// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DUALDESC_REVERSE_SEARCH_H_
#define SRC_POLY_POLY_DUALDESC_REVERSE_SEARCH_H_

// clang-format off
#include "Boost_bitset.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "POLY_Fundamental.h"
#include <limits>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_REVERSE_SEARCH
#endif

#ifdef DISABLE_DEBUG_REVERSE_SEARCH
#undef DEBUG_REVERSE_SEARCH
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_REVERSE_SEARCH
#endif

// A rewrite of the lrs dual description backend (POLY_DualDesc_lrslib.h),
// which is a line-by-line translation of David Avis' lrslib. The algorithm
// is the one of
//  * D. Avis, K. Fukuda, "A pivoting algorithm for convex hulls and vertex
//    enumeration of arrangements and polyhedra", Discrete Comput. Geom. 8
//    (1992) 295-313.
// with the standard lrslib refinements over the plain paper:
//  * lexicographic ratio test instead of Bland's rule, so the reverse
//    search runs on the lex-positive dictionaries of the (implicitly)
//    perturbed problem: the perturbed problem is nondegenerate, the root is
//    unique, and no special treatment of degenerate optimal dictionaries
//    (section 3.2 of the paper) is needed;
//  * integer pivoting (fraction-free Bareiss with the tracked determinant),
//    so the whole enumeration runs division-free over a ring;
//  * the lex-min test (Proposition 3.4) so each vertex / ray is reported
//    exactly once even on degenerate inputs.
//
// The formulation is the one of lrslib in H-mode: the rows of the input
// matrix EXT are read as inequalities EXT(i,.) . (1, x) >= 0 over
// x in R^{nbCol-1}, and the enumeration reports the vertices and extreme
// rays of that polyhedron, which are the facets of the cone spanned by the
// rows. The public entry points first prepend a zero column when column 0
// is not already zero (FirstColumnZero), so in production the system is
// homogeneous: the right-hand side column stays identically zero through
// every pivot, the only vertex is the apex (skipped), and every facet
// comes out as a ray.
//
// What is different from the old translation:
//  * flat contiguous tableau instead of a row-pointer array, std::vector
//    ownership throughout, no linked-list dictionary cache;
//  * backtracking through a stack of saved dictionaries -- one slot per
//    tree level, restored by O(1) buffer swaps -- so no re-pivoting on
//    backtrack (the old code cached only the last 10 dictionaries and
//    re-pivoted below that); beyond a depth cap the parent is recomputed
//    by a forward pivot as in the paper, which needs no stored data;
//  * the volume / triangulation entry points are not duplicated here, they
//    remain in the lrs namespace.
//
// The machine-integer fast path is the same as in the other backends: on a
// field input the kernel runs over the underlying ring, and first over
// TryInt64 with the deferred-overflow discipline; the overflow flag is
// checked before every emission, and on overflow the exact rerun skips the
// already-emitted prefix (the enumeration is deterministic).

namespace rev_search {

#ifdef POLY_LRS_TRY_SIMD
using RsTryInt = TrySimdInt64;
#else
using RsTryInt = TryCarryInt64;
#endif

// Above this depth the dictionary stack stops storing copies and the
// parent dictionary is recomputed by one forward pivot on backtrack.
constexpr int STACK_CAP = 4096;

template <typename T> inline int sign_int(T const &a) {
  if (a < 0)
    return -1;
  if (a > 0)
    return 1;
  return 0;
}

// Compare Na * Nb with Nc * Nd: 1 / 0 / -1 for > / = / <.
template <typename T>
inline int comprod(T const &Na, T const &Nb, T const &Nc, T const &Nd) {
  T prod1 = Na * Nb;
  T prod2 = Nc * Nd;
  if (prod1 > prod2)
    return 1;
  if (prod1 < prod2)
    return -1;
  return 0;
}

// The dictionary: tableau A of size (m+1) x (d+1) stored flat row-major,
// row 0 the cost row, column 0 the right-hand side. B/C are the basic and
// cobasic variable indices kept sorted ascending, Row/Col their tableau
// locations: the basic variable B[i] owns tableau row Row[i], the cobasic
// variable C[j] owns tableau column Col[j].
template <typename T> struct Dictionary {
  int m;   // number of input rows
  int d;   // current number of decision columns
  int depth;
  bool lexflag;
  T det;   // positive determinant of the current basis
  std::vector<T> Adata;      // (m+1) * (d+1)
  std::vector<int> B, Row;   // size m+1
  std::vector<int> C, Col;   // size d+1

  Dictionary(int _m, int _d) : m(_m), d(_d), depth(0), lexflag(true), det(1) {
    Adata.assign(static_cast<size_t>(m + 1) * (d + 1), T(0));
    B.resize(m + 1);
    Row.resize(m + 1);
    C.resize(d + 1);
    Col.resize(d + 1);
    for (int i = 0; i <= m; i++) {
      B[i] = (i == 0) ? 0 : d + i;
      Row[i] = i;
    }
    for (int j = 0; j < d; j++) {
      C[j] = j + 1;
      Col[j] = j + 1;
    }
    C[d] = m + d + 1;
    Col[d] = 0;
  }
  T *row(int i) { return Adata.data() + static_cast<size_t>(i) * (d + 1); }
  T const *row(int i) const {
    return Adata.data() + static_cast<size_t>(i) * (d + 1);
  }
  T &A(int i, int j) { return Adata[static_cast<size_t>(i) * (d + 1) + j]; }
  T const &A(int i, int j) const {
    return Adata[static_cast<size_t>(i) * (d + 1) + j];
  }
};

// The problem data that does not change during the search.
struct Problem {
  int m;           // input rows
  int n;           // input columns
  bool homogeneous;  // true if input column 0 is identically zero
  int lastdv;      // index of the last decision variable after preprocessing
  int nredundcol;  // number of dependent input columns
  std::vector<int> inequality;  // variable index - lastdv -> input row (1-based)
  std::vector<int> redundcol;   // input columns removed as dependent
  std::vector<int> minratio;    // scratch for the ratio test
};

// One saved dictionary for O(1) backtracking; buffers are swapped in and
// out so a slot's capacity is reused across the whole run.
template <typename T> struct SavedState {
  std::vector<T> Adata;
  std::vector<int> B, Row, C, Col;
  T det;
  int i, j;
  bool used = false;
};

// Fraction-free Bareiss pivot on basis position bas, cobasis position cob.
// The invariant hoisting mirrors the old code: the copies Ars / Ais / the
// row pointer break the aliasing so the compiler keeps them in registers.
template <typename T> void pivot(Dictionary<T> &dict, int bas, int cob) {
  int const r = dict.Row[bas];
  int const s = dict.Col[cob];
  int const d = dict.d;
  int const m = dict.m;
  T const Ars = dict.A(r, s);
  // The Bareiss divisor is the previous determinant carrying the sign of
  // the pivot element; the new determinant is |Ars|.
  T const Ndet = (Ars > 0) ? dict.det : -dict.det;
  T *Ar = dict.row(r);
  if (dict.det == 1) {
    // The division by the previous determinant is exact; when it is one
    // (frequent on small-coordinate input) it is skipped entirely, with
    // the sign folded into the update for a negative pivot element.
    if (Ars > 0) {
      for (int i = 0; i <= m; i++) {
        if (i == r)
          continue;
        T *Ai = dict.row(i);
        T const Ais = Ai[s];
        for (int j = 0; j < s; j++)
          Ai[j] = Ai[j] * Ars - Ais * Ar[j];
        for (int j = s + 1; j <= d; j++)
          Ai[j] = Ai[j] * Ars - Ais * Ar[j];
      }
    } else {
      for (int i = 0; i <= m; i++) {
        if (i == r)
          continue;
        T *Ai = dict.row(i);
        T const Ais = Ai[s];
        for (int j = 0; j < s; j++)
          Ai[j] = Ais * Ar[j] - Ai[j] * Ars;
        for (int j = s + 1; j <= d; j++)
          Ai[j] = Ais * Ar[j] - Ai[j] * Ars;
      }
    }
  } else {
    bool magic_done = false;
    if constexpr (is_try_int<T>::value) {
      // The Bareiss division is exact, so over the machine-integer type it
      // can be done as a shift and a wrapping multiplication with the
      // inverse of the odd part of the divisor modulo 2^64 (the classical
      // exact-division trick), which is several times cheaper than a
      // hardware division. The numerator products stay overflow-checked:
      // if they wrapped, the flag is already set and the garbage quotient
      // is discarded by the rerun, the same discipline as the checked
      // division it replaces -- with the advantage that this path cannot
      // trap on a garbage divisor.
      static_assert(sizeof(T) == sizeof(int64_t));
      int64_t const Draw = Ndet.get_const_val();
      uint64_t D = (Draw < 0) ? -static_cast<uint64_t>(Draw)
                              : static_cast<uint64_t>(Draw);
      if (D != 0) {
        int const shift = __builtin_ctzll(D);
        uint64_t const Dodd = D >> shift;
        uint64_t inv = Dodd;
        for (int it = 0; it < 5; it++)
          inv *= 2 - Dodd * inv;
        if (Draw < 0)
          inv = 0 - inv;
        auto exact_div = [&](T const &num) -> int64_t {
          int64_t const x = num.get_const_val();
          return static_cast<int64_t>(
              static_cast<uint64_t>(x >> shift) * inv);
        };
        for (int i = 0; i <= m; i++) {
          if (i == r)
            continue;
          T *Ai = dict.row(i);
          T const Ais = Ai[s];
          for (int j = 0; j < s; j++)
            Ai[j].get_val() = exact_div(Ai[j] * Ars - Ais * Ar[j]);
          for (int j = s + 1; j <= d; j++)
            Ai[j].get_val() = exact_div(Ai[j] * Ars - Ais * Ar[j]);
        }
        magic_done = true;
      }
    }
    if (!magic_done) {
      for (int i = 0; i <= m; i++) {
        if (i == r)
          continue;
        T *Ai = dict.row(i);
        T const Ais = Ai[s];
        for (int j = 0; j < s; j++)
          Ai[j] = (Ai[j] * Ars - Ais * Ar[j]) / Ndet;
        for (int j = s + 1; j <= d; j++)
          Ai[j] = (Ai[j] * Ars - Ais * Ar[j]) / Ndet;
      }
    }
  }
  if (Ars > 0) {
    for (int j = 0; j <= d; j++)
      if (Ar[j] != 0)
        Ar[j] = -Ar[j];
  } else {
    for (int i = 0; i <= m; i++) {
      T &Ais = dict.A(i, s);
      if (Ais != 0)
        Ais = -Ais;
    }
  }
  dict.A(r, s) = Ndet;
  dict.det = (Ars > 0) ? Ars : T(-Ars);
  // Deferred overflow check: no-op except for the TryInt64 types.
  terminate_in_arithmetic_error<T>();
}

// Restore the sorted order of an index array after one element changed,
// keeping the companion location array aligned.
inline void reorder1(std::vector<int> &a, std::vector<int> &b, int newone,
                     int range) {
  while (newone > 0 && a[newone] < a[newone - 1]) {
    std::swap(a[newone], a[newone - 1]);
    std::swap(b[newone], b[newone - 1]);
    newone--;
  }
  while (newone < range - 1 && a[newone] > a[newone + 1]) {
    std::swap(a[newone], a[newone + 1]);
    std::swap(b[newone], b[newone + 1]);
    newone++;
  }
}

// Update B/C after the pivot of B[i] with C[j]; on return i and j are the
// new positions of the entered and left variables.
template <typename T> void update(Dictionary<T> &dict, int *i, int *j) {
  int const leave = dict.B[*i];
  int const enter = dict.C[*j];
  dict.B[*i] = enter;
  reorder1(dict.B, dict.Row, *i, dict.m + 1);
  dict.C[*j] = leave;
  reorder1(dict.C, dict.Col, *j, dict.d);
  for (*i = 1; dict.B[*i] != enter; (*i)++) {
  }
  for (*j = 0; dict.C[*j] != leave; (*j)++) {
  }
}

// Lexicographic minimum ratio test on tableau column col. Returns the
// basis position of the unique leaving variable, or 0 when the column is
// nonnegative on the rows (a ray). The tie-breaking walks the columns of
// the basis inverse in basis order, which realizes the lexicographic
// perturbation.
template <typename T> int lex_ratio(Dictionary<T> &dict, Problem &prob, int col) {
  int const m = dict.m;
  int const d = dict.d;
  int const lastdv = prob.lastdv;
  std::vector<int> &minratio = prob.minratio;
  int degencount = 0;
  for (int j = lastdv + 1; j <= m; j++) {
    if (dict.A(dict.Row[j], col) < 0)
      minratio[degencount++] = j;
  }
  if (degencount == 0)
    return 0;
  T Nmin(0), Dmin(0);
  int ratiocol = 0;    // tableau column being compared, first the rhs
  int start = 0;       // window into minratio
  int bindex = d + 1;  // next basic variable index to consider
  int cindex = 0;      // next cobasic position to consider
  // The variable index whose basis-inverse column is tested; the value d
  // stands for the right-hand side pass. On homogeneous input the rhs
  // column is identically zero, every candidate ties, and the pass is
  // skipped outright.
  int basicindex = prob.homogeneous ? d + 1 : d;
  int nstart = 0, ndegencount = 0;
  while (degencount > 1) {
    if (dict.B[bindex] == basicindex) {
      // Identity column of the basis inverse: the candidate equal to
      // bindex leaves the tie immediately, everyone else stays.
      if (minratio[start] == bindex) {
        start++;
        degencount--;
      }
      bindex++;
    } else {
      bool firstime = true;
      if (basicindex != d)
        ratiocol = dict.Col[cindex++];
      for (int j = start; j < start + degencount; j++) {
        int const i = dict.Row[minratio[j]];
        int comp = 1;
        if (firstime) {
          firstime = false;
        } else {
          if (Nmin > 0 || dict.A(i, ratiocol) < 0) {
            if (Nmin < 0 || dict.A(i, ratiocol) > 0)
              comp = comprod(Nmin, dict.A(i, col), dict.A(i, ratiocol), Dmin);
            else
              comp = -1;
          } else if (Nmin == 0 && dict.A(i, ratiocol) == 0) {
            comp = 0;
          }
          if (ratiocol == 0)
            comp = -comp;  // signs reversed for the rhs column
        }
        if (comp == 1) {
          nstart = j;
          Nmin = dict.A(i, ratiocol);
          Dmin = dict.A(i, col);
          ndegencount = 1;
        } else if (comp == 0) {
          minratio[nstart + ndegencount++] = minratio[j];
        }
      }
      degencount = ndegencount;
      start = nstart;
    }
    basicindex++;
  }
  return minratio[start];
}

// Forward pivot selection: smallest cobasic variable with positive cost,
// leaving variable by the lexicographic ratio test.
template <typename T>
bool selectpivot(Dictionary<T> &dict, Problem &prob, int *r, int *s) {
  int const d = dict.d;
  *r = 0;
  *s = d;
  int j = 0;
  while (j < d && dict.A(0, dict.Col[j]) <= 0)
    j++;
  if (j < d) {
    *s = j;
    *r = lex_ratio(dict, prob, dict.Col[j]);
    if (*r != 0)
      return true;
  }
  return false;
}

// Reverse pivot test: true if pivoting B[*r] with C[s] gives a dictionary
// whose forward pivot is exactly that pair, i.e. the current dictionary is
// the child of that dictionary in the reverse search tree. The cost row of
// the would-be parent is examined by sign combinations, no pivot is done.
template <typename T>
bool reverse(Dictionary<T> &dict, Problem &prob, int *r, int s) {
  int const d = dict.d;
  int const col = dict.Col[s];
  if (dict.A(0, col) >= 0)
    return false;
  *r = lex_ratio(dict, prob, col);
  if (*r == 0)
    return false;
  int const row = dict.Row[*r];
  for (int i = 0; i < d && dict.C[i] < dict.B[*r]; i++) {
    if (i != s) {
      int const j = dict.Col[i];
      if (dict.A(0, j) > 0 || dict.A(row, j) < 0) {
        if ((dict.A(0, j) >= 0 && dict.A(row, j) <= 0) ||
            comprod(dict.A(0, j), dict.A(row, col), dict.A(0, col),
                    dict.A(row, j)) == -1)
          return false;
      }
    }
  }
  return true;
}

// True if A[r][s] realizes a minimum ratio for column s.
template <typename T> bool ismin(Dictionary<T> &dict, int r, int s) {
  int const m = dict.m;
  for (int i = 1; i <= m; i++) {
    if (i != r && dict.A(i, s) < 0 &&
        comprod(dict.A(i, 0), dict.A(r, s), dict.A(i, s), dict.A(r, 0)) != 0)
      return false;
  }
  return true;
}

// Lex-min test (Proposition 3.4): the basis is not the lexicographically
// minimal one for the current basic solution iff a basic variable with a
// zero value could be exchanged against a smaller cobasic one.
template <typename T> bool lexmin(Dictionary<T> &dict, Problem &prob, int col) {
  int const m = dict.m;
  int const d = dict.d;
  for (int i = prob.lastdv + 1; i <= m; i++) {
    int const r = dict.Row[i];
    if (dict.A(r, col) == 0) {
      for (int j = 0; j < d; j++) {
        int const s = dict.Col[j];
        if (dict.B[i] > dict.C[j]) {
          if (dict.A(r, 0) == 0) {
            if (dict.A(r, s) != 0)
              return false;
          } else if (dict.A(r, s) < 0 && ismin(dict, r, s)) {
            return false;
          }
        }
      }
    }
  }
  return true;
}

// Pivot the decision variables into the basis, trying to pivot out the
// slacks in the order row m, m-1, ..., 1. Decision variables still cobasic
// afterwards are dependent input columns; they are recorded and removed.
template <typename T>
void get_initial_basis(Dictionary<T> &dict, Problem &prob) {
  int const m = dict.m;
  int const d = dict.d;
  for (int jj = 0; jj < m; jj++) {
    int const target = d + (m - jj);  // slack of input row m - jj
    int i = 0;
    while (i <= m && dict.B[i] != target)
      i++;
    if (i <= m) {
      int k = 0;
      while (dict.C[k] <= d && dict.A(dict.Row[i], dict.Col[k]) == 0)
        k++;
      if (dict.C[k] <= d) {
        pivot(dict, i, k);
        update(dict, &i, &k);
      }
    }
  }
  int nredundcol = 0;
  for (int k = 0; k < d && dict.C[k] <= d; k++)
    prob.redundcol[nredundcol++] = dict.C[k];
  prob.nredundcol = nredundcol;
  prob.lastdv = d - nredundcol;
}

// Remove the cobasic variable at position k (a dependent decision column)
// from the problem, renumbering all larger variable indices down by one.
template <typename T> void remove_cobasic(Dictionary<T> &dict, int k) {
  int const m = dict.m;
  int const d = dict.d;
  int const cindex = dict.C[k];
  int const deloc = dict.Col[k];
  for (int i = 1; i <= m; i++)
    if (dict.B[i] > cindex)
      dict.B[i]--;
  for (int j = k; j < d; j++) {
    dict.C[j] = dict.C[j + 1] - 1;
    dict.Col[j] = dict.Col[j + 1];
  }
  if (deloc != d) {
    for (int i = 0; i <= m; i++)
      dict.A(i, deloc) = dict.A(i, d);
    int j = 0;
    while (dict.Col[j] != d)
      j++;
    dict.Col[j] = deloc;
  }
  // Shrink the flat tableau to the new stride d (columns 0..d-1 kept).
  int const dnew = d - 1;
  std::vector<T> Anew(static_cast<size_t>(m + 1) * (dnew + 1));
  for (int i = 0; i <= m; i++)
    for (int j = 0; j <= dnew; j++)
      Anew[static_cast<size_t>(i) * (dnew + 1) + j] = dict.A(i, j);
  dict.Adata = std::move(Anew);
  dict.d = dnew;
  dict.C.resize(dnew + 1);
  dict.Col.resize(dnew + 1);
}

// Dual-Bland pivots to primal feasibility. The cost row is zero at this
// point so no ratio test is needed. For the homogeneous inputs of the
// public entry points the loop body never executes.
template <typename T> bool primal_feasible(Dictionary<T> &dict, Problem &prob) {
  int const m = dict.m;
  while (true) {
    int i = prob.lastdv + 1;
    while (i <= m && dict.A(dict.Row[i], 0) >= 0)
      i++;
    if (i > m)
      return true;
    int j = 0;
    while (j < dict.d && dict.A(dict.Row[i], dict.Col[j]) <= 0)
      j++;
    if (j >= dict.d)
      return false;
    pivot(dict, i, j);
    update(dict, &i, &j);
  }
}

// Build the root dictionary: initial basis, removal of dependent columns,
// primal feasibility, the all-negative cost row that makes the root
// optimal, and the relabeling that makes the basis exactly {0, 1, ..., m}
// with inequality[] recording which input row each variable index carries.
template <typename T>
bool get_first_basis(Dictionary<T> &dict, Problem &prob) {
  int const m = dict.m;
  for (int j = 0; j <= dict.d; j++)
    dict.A(0, j) = T(0);
  get_initial_basis(dict, prob);
  for (int i = 1; i <= m; i++)
    prob.inequality[i] = i;
  for (int i = 0; i < prob.nredundcol; i++)
    remove_cobasic(dict, 0);
  if (!primal_feasible(dict, prob))
    return false;
  for (int j = 1; j <= dict.d; j++)
    dict.A(0, j) = -dict.det;
  dict.A(0, 0) = T(0);
  int const lastdv = prob.lastdv;
  while (dict.C[0] <= m) {
    int const i = dict.C[0];
    std::swap(prob.inequality[dict.B[i] - lastdv],
              prob.inequality[dict.C[0] - lastdv]);
    dict.C[0] = dict.B[i];
    dict.B[i] = i;
    reorder1(dict.C, dict.Col, 0, dict.d);
  }
  return true;
}

// The incidence of the solution in column col: the input rows whose slack
// is cobasic (except the ray column itself) or basic with value zero and a
// zero coefficient in the ray column.
template <typename T>
void set_face(Dictionary<T> const &dict, Problem const &prob, int col,
              Face &f) {
  int const nbRow = prob.m;
  for (int i = 0; i < nbRow; i++)
    f[i] = 0;
  for (int i = 0; i < dict.d; i++) {
    if (dict.Col[i] != col) {
      int const idx = prob.inequality[dict.C[i] - prob.lastdv] - 1;
      f[idx] = 1;
    }
  }
  for (int i = prob.lastdv + 1; i <= dict.m; i++) {
    int const iRow = dict.Row[i];
    if (dict.A(iRow, 0) == 0) {
      if (col == 0 || dict.A(iRow, col) == 0)
        f[iRow - 1] = 1;
    }
  }
}

// Extract the basic solution (col == 0) into output, skipping bases that
// are not lex-min. Dependent input columns get a zero coordinate.
template <typename T>
bool get_vertex(Dictionary<T> const &dict, Problem const &prob, T *output) {
  if (!dict.lexflag)
    return false;
  output[0] = dict.det;
  int i = 1;
  int ired = 0;
  for (int ind = 1; ind < prob.n; ind++) {
    if (ired < prob.nredundcol && prob.redundcol[ired] == ind) {
      output[ind] = T(0);
      ired++;
    } else {
      output[ind] = dict.A(dict.Row[i], 0);
      i++;
    }
  }
  return true;
}

// Extract the ray of cobasic column col into output.
template <typename T>
bool get_ray(Dictionary<T> const &dict, Problem const &prob, int col,
             T *output) {
  int i = 1;
  int ired = 0;
  for (int ind = 0; ind < prob.n; ind++) {
    if (ind == 0) {
      output[0] = T(0);
    } else if (ired < prob.nredundcol && prob.redundcol[ired] == ind) {
      output[ind] = T(0);
      ired++;
    } else {
      output[ind] = dict.A(dict.Row[i], col);
      i++;
    }
  }
  return true;
}

// Check whether column col of the current dictionary carries an output:
// the basic solution for col == 0, an extreme ray otherwise (negative cost
// coefficient, nonnegative column on the rows, lex-min basis).
template <typename T>
bool get_solution(Dictionary<T> &dict, Problem &prob, T *output, int col) {
  if (col == 0)
    return get_vertex(dict, prob, output);
  if (dict.A(0, col) >= 0)
    return false;
  int j = prob.lastdv + 1;
  while (j <= dict.m && dict.A(dict.Row[j], col) >= 0)
    j++;
  if (j <= dict.m)
    return false;
  if (lexmin(dict, prob, col))
    return get_ray(dict, prob, col, output);
  return false;
}

// The dictionary stack. push saves the current dictionary (by copy, the
// buffers of a previously popped slot are reused); pop restores by O(1)
// swaps. Beyond STACK_CAP no copy is stored and the pop recomputes the
// parent with one forward pivot, which is self-contained.
template <typename T> struct DictStack {
  std::vector<SavedState<T>> slots;
  int sp = 0;
  void push(Dictionary<T> const &dict, int i, int j) {
    if (sp >= STACK_CAP) {
      sp++;
      return;
    }
    if (sp >= static_cast<int>(slots.size()))
      slots.resize(sp + 1);
    SavedState<T> &sl = slots[sp];
    sl.Adata = dict.Adata;
    sl.B = dict.B;
    sl.Row = dict.Row;
    sl.C = dict.C;
    sl.Col = dict.Col;
    sl.det = dict.det;
    sl.i = i;
    sl.j = j;
    sl.used = true;
    sp++;
  }
  // Restore the parent; returns false if it must be recomputed by pivot.
  bool pop(Dictionary<T> &dict, int *i, int *j) {
    sp--;
    if (sp >= STACK_CAP)
      return false;
    SavedState<T> &sl = slots[sp];
    dict.Adata.swap(sl.Adata);
    dict.B.swap(sl.B);
    dict.Row.swap(sl.Row);
    dict.C.swap(sl.C);
    dict.Col.swap(sl.Col);
    dict.det = sl.det;
    *i = sl.i;
    *j = sl.j;
    return true;
  }
};

// One step of the depth-first reverse search: move to the next dictionary
// of the tree, false when the search is exhausted. On return the new
// dictionary's lexflag is set for the vertex test.
template <typename T>
bool next_basis(Dictionary<T> &dict, Problem &prob, DictStack<T> &stack,
                bool backtrack) {
  int i = 0, j = 0;
  int const d = dict.d;
  int const m = dict.m;
  if (backtrack && dict.depth == 0)
    return false;
  while (j < d || dict.B[m] != m) {
    if (backtrack) {
      backtrack = false;
      dict.depth--;
      if (!stack.pop(dict, &i, &j)) {
        selectpivot(dict, prob, &i, &j);
        pivot(dict, i, j);
        update(dict, &i, &j);
        // selectpivot returned the forward pivot to the parent; i, j are
        // now the positions of the variables that moved, which is exactly
        // the reverse pivot position we descended through.
      }
      j++;
    }
    while (j < d && !reverse(dict, prob, &i, j))
      j++;
    if (j == d) {
      backtrack = true;
      if (dict.depth == 0)
        return false;
    } else {
      stack.push(dict, i, j);
      dict.depth++;
      pivot(dict, i, j);
      update(dict, &i, &j);
      dict.lexflag = lexmin(dict, prob, 0);
      return true;
    }
  }
  return false;
}

// Run the full enumeration on EXT, calling f(dict, prob, col, output) for
// every output (vertex or ray) except the very first one, which is the
// basic solution of the root dictionary (the apex for homogeneous input).
template <typename T, typename F>
void Kernel_DualDescription(MyMatrix<T> const &EXT, F const &f) {
  int const nbRow = EXT.rows();
  int const nbCol = EXT.cols();
  int const m = nbRow;
  int const d = nbCol - 1;
  Dictionary<T> dict(m, d);
  Problem prob;
  prob.m = m;
  prob.n = nbCol;
  prob.homogeneous = [&]() -> bool {
    for (int iRow = 0; iRow < nbRow; iRow++)
      if (EXT(iRow, 0) != 0)
        return false;
    return true;
  }();
  prob.lastdv = d;
  prob.nredundcol = 0;
  prob.inequality.assign(m + 1, 0);
  prob.redundcol.assign(d + 1, 0);
  prob.minratio.assign(m + 1, 0);
  for (int iRow = 0; iRow < nbRow; iRow++)
    for (int j = 0; j < nbCol; j++)
      dict.A(iRow + 1, j) = EXT(iRow, j);
  if (!get_first_basis(dict, prob)) {
    std::cerr << "RS: failure in get_first_basis\n";
    throw TerminalException{1};
  }
  std::vector<T> output(prob.n + 1);
  DictStack<T> stack;
  bool is_first = true;
  bool backtrack = false;
  while (true) {
    for (int col = 0; col <= dict.d; col++) {
      if (get_solution(dict, prob, output.data(), col)) {
        if (!is_first)
          f(dict, prob, col, output.data());
        is_first = false;
      }
    }
    if (!next_basis(dict, prob, stack, backtrack))
      break;
    backtrack = false;
  }
}

template <typename T> MyMatrix<T> FirstColumnZero(MyMatrix<T> const &M) {
  int const nbRow = M.rows();
  for (int iRow = 0; iRow < nbRow; iRow++)
    if (M(iRow, 0) != 0)
      return AddFirstZeroColumn(M);
  return M;
}

template <typename T>
std::pair<MyMatrix<T>, int> FirstColumnZeroCond(MyMatrix<T> const &M) {
  int const nbRow = M.rows();
  for (int iRow = 0; iRow < nbRow; iRow++)
    if (M(iRow, 0) != 0)
      return {AddFirstZeroColumn(M), 1};
  return {M, 0};
}

template <typename T, typename Tw> inline T rs_output_convert(Tw const &val) {
  if constexpr (std::is_same_v<Tw, T>) {
    return val;
  } else if constexpr (is_try_int<Tw>::value) {
    return ConvertFromTryInt64<T>(val);
  } else {
    return UniversalScalarConversion<T, Tw>(val);
  }
}

// Ring kernel with the TryInt64 first attempt: the overflow flag is
// checked before every emission, and on overflow the exact rerun skips the
// n_emitted facets already handed out (the enumeration is deterministic
// for identical values).
template <typename Tring, typename Ffacet>
void Kernel_DualDescription_ring(MyMatrix<Tring> const &EXTring,
                                 Ffacet f_facet) {
  if constexpr (use_try_int64<Tring>::value) {
    size_t n_emitted = 0;
    try {
      MyMatrix<RsTryInt> EXTtry = ConvertMatrixToTryInt64<RsTryInt>(EXTring);
      auto f_try = [&](Dictionary<RsTryInt> &dict, Problem &prob,
                       int const &col, RsTryInt *out) -> void {
        terminate_in_arithmetic_error<RsTryInt>();
        f_facet(dict, prob, col, out);
        n_emitted++;
      };
      Kernel_DualDescription(EXTtry, f_try);
      terminate_in_arithmetic_error<RsTryInt>();
      return;
    } catch (TryIntException const &) {
    }
    size_t idx = 0;
    auto f_skip = [&](Dictionary<Tring> &dict, Problem &prob, int const &col,
                      Tring *out) -> void {
      if (idx >= n_emitted)
        f_facet(dict, prob, col, out);
      idx++;
    };
    Kernel_DualDescription(EXTring, f_skip);
  } else {
    Kernel_DualDescription(EXTring, f_facet);
  }
}

// Field-vs-ring dispatch: a field input is scaled row by row to the
// underlying ring, a ring input runs directly.
template <typename T, typename Ffacet>
void Kernel_DualDescription_process(MyMatrix<T> const &EXTwork,
                                    Ffacet f_facet) {
  if constexpr (is_ring_field<T>::value) {
    using Tring = typename underlying_ring<T>::ring_type;
    int const nbRow = EXTwork.rows();
    int const nbCol = EXTwork.cols();
    MyMatrix<Tring> EXTring(nbRow, nbCol);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      MyVector<T> eRow =
          NonUniqueScaleToIntegerVector(GetMatrixRow(EXTwork, iRow));
      AssignMatrixRow(EXTring, iRow, UniversalVectorConversion<Tring, T>(eRow));
    }
    Kernel_DualDescription_ring(EXTring, f_facet);
  } else {
    Kernel_DualDescription_ring(EXTwork, f_facet);
  }
}

template <typename T> vectface DualDescription_incd(MyMatrix<T> const &EXT) {
  MyMatrix<T> EXTwork = FirstColumnZero(EXT);
  size_t const nbRow = EXTwork.rows();
  vectface ListIncd(nbRow);
  Face face(nbRow);
  auto f_facet = [&](auto &dict, Problem &prob, int const &col,
                     [[maybe_unused]] auto *out) -> void {
    set_face(dict, prob, col, face);
    ListIncd.push_back(face);
  };
  Kernel_DualDescription_process(EXTwork, f_facet);
  return ListIncd;
}

template <typename T> MyMatrix<T> DualDescription(MyMatrix<T> const &EXT) {
  std::pair<MyMatrix<T>, int> pair = FirstColumnZeroCond(EXT);
  MyMatrix<T> const &EXTwork = pair.first;
  int const shift = pair.second;
  int const nbColRed = EXTwork.cols() - shift;
  std::vector<MyVector<T>> ListVect;
  MyVector<T> V(nbColRed);
  auto f_facet = [&]([[maybe_unused]] auto &dict, [[maybe_unused]] Problem &prob,
                     [[maybe_unused]] int const &col, auto *out) -> void {
    for (int i = 0; i < nbColRed; i++)
      V(i) = rs_output_convert<T>(out[i + shift]);
    ListVect.push_back(V);
  };
  Kernel_DualDescription_process(EXTwork, f_facet);
  return MatrixFromVectorFamily(ListVect);
}

// The pair handed to f_process is reused across facets.
template <typename T, typename Fprocess>
void DualDescriptionFaceIneq(MyMatrix<T> const &EXT, Fprocess f_process) {
  std::pair<MyMatrix<T>, int> ePair = FirstColumnZeroCond(EXT);
  MyMatrix<T> const &EXTwork = ePair.first;
  int const shift = ePair.second;
  int const nbRow = EXTwork.rows();
  int const nbColRed = EXTwork.cols() - shift;
  std::pair<Face, MyVector<T>> pair{Face(nbRow), MyVector<T>(nbColRed)};
  auto f_facet = [&](auto &dict, Problem &prob, int const &col,
                     auto *out) -> void {
    for (int i = 0; i < nbColRed; i++)
      pair.second(i) = rs_output_convert<T>(out[i + shift]);
    set_face(dict, prob, col, pair.first);
    f_process(pair);
  };
  Kernel_DualDescription_process(EXTwork, f_facet);
}

// clang-format off
}  // namespace rev_search
#endif  // SRC_POLY_POLY_DUALDESC_REVERSE_SEARCH_H_
// clang-format on

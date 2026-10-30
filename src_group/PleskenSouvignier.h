// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_PLESKENSOUVIGNIER_H_
#define SRC_GROUP_PLESKENSOUVIGNIER_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include <algorithm>
#include <limits>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_PLESKEN_SOUVIGNIER
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_PLESKEN_SOUVIGNIER
#endif

#ifdef TIMINGS
#define TIMINGS_PLESKEN_SOUVIGNIER
#endif

/*
  The Plesken-Souvignier backtrack algorithm for the automorphism group of
  a lattice and for isometry testing:

  W. Plesken, B. Souvignier, Computing isometries of lattices, Journal of
  Symbolic Computation (1997) 24, 327-334.

  The structure follows Bernd Souvignier's programs AUTO and ISOM through
  the port in Hecke.jl (src/QuadForm/Morphism.jl), with one deviation
  noted in ps_stab.

  The problem solved here is the one of PolytopeEquiStabInt.h --
  stabilizer and equivalence of a configuration of vectors under
  GL_n(Z) -- by a different method. There the full graph of pairwise
  scalar products is built and handed to a graph automorphism program,
  which is quadratic in the number of vectors before the graph code even
  starts. Here nothing quadratic in the family is ever built: images of
  one basis are searched in the family, and the pairwise products are
  computed only along the search. On families of thousands of vectors
  with moderate symmetry this is orders of magnitude cheaper; the graph
  method remains preferable when a canonical form is wanted, which the
  backtrack does not provide.

  Conventions, matching the rest of the repository:
  * Vectors are rows. An automorphism g satisfies g * M * g^T = M for
    every matrix M of the configuration and acts on a row vector v as
    v -> v * g.
  * The isometry returned by PleskenSouvignierIsometry satisfies
    P * M1 * P^T = M2, as TestIntEquivalence_ListMat_Vdiag does.

  The vector family S is given as SHVhalf, one representative per
  antipodal pair; the pair set is what matters and the sign of each
  representative is normalized internally. Requirements on the input:

  (1) ListMat[0] is symmetric positive definite and every other matrix is
      symmetric. Positive definiteness is what makes the pointwise
      stabilizer of a basis trivial, which the order computation uses.
  (2) The set {±v} is invariant under the group being computed. Any set
      of the form {v : v M v^T <= bound for the matrices M} is.
  (3) S has full rank. It does not have to contain the standard basis
      and does not even have to span Z^n.

  The searched basis is drawn from the FAMILY: the standard basis when S
  contains it (the vector set of the paper, {v : v ListMat[0] v^T <=
  max_i ListMat[0](i,i)}, always does, and an LLL-reduced Gram matrix
  keeps it small -- see src_latt/LatticePleskenSouvignier.h for the
  wrappers), and otherwise n independent rows of S, shortest first. The
  images of the basis under any automorphism are in S because S itself
  is, which is what the search needs.

  A complete assignment of images v_1, ..., v_n in S with
  (v_i, v_j)_M = (b_i, b_j)_M for every M defines the matrix
  g = B^{-1} V; the products against the full-rank basis force
  g M g^T = M globally, and |det g| = 1 follows from ListMat[0] being
  definite. When the basis spans Z^n (index 1, always the case for the
  standard basis) g is integral by construction and no test is needed;
  when it only spans a finite index sublattice, the leaf checks that
  B^{-1} V is integral and discards the rational solutions. That is the
  price of a general family: a family spanning a proper sublattice L'
  finds exactly Aut(Z^n, ListMat) all the same -- every integral
  automorphism preserves the invariant S, hence L', hence appears; the
  rational isometries of L' that do not preserve Z^n are the discarded
  leaves, and when they vastly outnumber the integral ones the subspace
  machinery of MatrixGroup.h (compute the rational group once, cut it to
  the integral subgroup) is the better tool.

  The search prunes with the fingerprint of the paper (Section 4): the
  basis is reordered so that the number of candidate images per level is
  minimal, and a partial assignment is abandoned unless the number of
  candidates at the next level is EXACTLY the fingerprint value, an
  automorphism mapping candidate sets bijectively. The automorphism
  group is assembled level by level along the stabilizer chain of the
  basis (Section 8), with the orbits of the already-found group pruning
  the candidates, and Schreier elements harvested from the orbit
  computations (ps_stab) supplementing the generators found by
  backtracking.

  The vector sums of Section 5 are implemented with a depth parameter,
  defaulting to n/10 rounded as in Hecke; without them the fingerprint
  alone loses rank-14 lattices of moderate symmetry to combinatorial
  explosion (the paper reports 50 h without result at depth 0 on a
  dimension-16 lattice whose depth-5 run takes 14 s). At the level where
  the images v_1, ..., v_l are chosen, every vector w of the family
  carries the tuple of its scalar products with the last `depth` chosen
  images, taken across all forms and normalized in sign. A partial
  automorphism maps the tuple classes of the identity prefix bijectively
  onto those of the chosen prefix, so the branch dies unless the tuples
  hit the precomputed identity classes with the same multiplicities; and
  it maps the per-class vector sums onto each other, so the Gram
  matrices of the sums have to agree as well. Not implemented: the
  Bacher polynomials (Section 6), an additional invariant for the
  remaining hard cases, pluggable into ps_cand the same way.
 */

// One vector per antipodal pair, sign-normalized so that the first
// non-zero coordinate is positive. A point of the search is a signed
// 1-based index: +k is row k-1 of V, -k its negative.
template <typename Tint> struct PleskenSouvignierVectorSet {
  int n;
  int m;
  MyMatrix<Tint> V;
  std::unordered_map<MyVector<Tint>, int> Map;
  // The signed index of w, or 0 when w is not in the family. A 0 from a
  // vector that an automorphism produced means the family is not
  // invariant, which the callers treat as an input error.
  int find_point(MyVector<Tint> const &w) const {
    int len = w.size();
    int ipos = 0;
    while (ipos < len && w(ipos) == 0) {
      ipos++;
    }
    if (ipos == len) {
      return 0;
    }
    if (w(ipos) > 0) {
      auto iter = Map.find(w);
      if (iter == Map.end()) {
        return 0;
      }
      return iter->second + 1;
    }
    MyVector<Tint> wNeg = -w;
    auto iter = Map.find(wNeg);
    if (iter == Map.end()) {
      return 0;
    }
    return -(iter->second + 1);
  }
};

template <typename Tint>
PleskenSouvignierVectorSet<Tint>
PleskenSouvignierBuildVectorSet(MyMatrix<Tint> const &SHVhalf) {
  int m = SHVhalf.rows();
  int n = SHVhalf.cols();
  PleskenSouvignierVectorSet<Tint> VS;
  VS.n = n;
  VS.m = m;
  VS.V = MyMatrix<Tint>(m, n);
  for (int j = 0; j < m; j++) {
    int ipos = 0;
    while (ipos < n && SHVhalf(j, ipos) == 0) {
      ipos++;
    }
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
    if (ipos == n) {
      std::cerr << "PS: the vector family contains the zero vector\n";
      throw TerminalException{1};
    }
#endif
    Tint sign = (SHVhalf(j, ipos) > 0) ? Tint(1) : Tint(-1);
    for (int i = 0; i < n; i++) {
      VS.V(j, i) = sign * SHVhalf(j, i);
    }
    VS.Map[GetMatrixRow(VS.V, j)] = j;
  }
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  if (static_cast<int>(VS.Map.size()) != m) {
    std::cerr << "PS: the vector family has duplicated antipodal pairs\n";
    throw TerminalException{1};
  }
#endif
  return VS;
}

/*
  The vector-sum data of one level (Section 5 of the paper), computed
  for the identity prefix: the distinct sign-normalized tuples of scalar
  products with the trailing dep basis vectors, their multiplicities,
  and the Gram matrices of the per-tuple vector sums. The search
  compares the same data computed for its chosen prefix.
 */
template <typename Tint> struct PleskenSouvignierVectorSumLevel {
  // Number of trailing positions used at this level, min(depth, level).
  int dep;
  std::unordered_map<MyVector<Tint>, int> TupleMap;
  std::vector<int> counts;
  // SumNorm[iMat][a] is the norm of the sum of the vectors of class a.
  std::vector<std::vector<Tint>> SumNorm;
  int n_tuple() const { return counts.size(); }
};

/*
  The search context of one lattice: the family, the product tables, the
  fingerprint, and the group data filled by the automorphism search.

  W[iMat] = V * ListMat[iMat], so that the scalar product of rows j and k
  is the dot product of row j of V with row k of W[iMat], and the product
  of row j with a standard basis vector e_p is the single entry
  W[iMat](j, p). Lengths(j, iMat) is the norm of row j. Nothing quadratic
  in the family size is ever built.
 */
template <typename Tint> struct PleskenSouvignierContext {
  std::vector<MyMatrix<Tint>> ListMat;
  PleskenSouvignierVectorSet<Tint> VS;
  std::vector<MyMatrix<Tint>> W;
  MyMatrix<Tint> Lengths;
  /*
    The searched basis, as rows of the family: bas holds the family row
    indices before the fingerprint reordering, Fmat[iMat] their pairwise
    products, so that the chain-position-i basis vector is row
    bas[per[i]] and its products are Fmat[iMat](per[i], per[k]). For the
    standard basis Fmat[iMat] is ListMat[iMat] itself.

    Badj and Bden describe the inverse of the position-ordered basis
    matrix B (row i = V.row(bas[per[i]])): B^{-1} = Badj / Bden with
    Badj integral and Bden = |det B|, so a leaf matrix is
    (Badj * A) / Bden with an exact-division check; Bden == 1 is the
    unimodular case where every leaf is integral.
   */
  std::vector<int> bas;
  std::vector<MyMatrix<Tint>> Fmat;
  MyMatrix<Tint> Badj;
  Tint Bden;
  /*
    The adjugate and determinant of ListMat[0], with adjM0 / detM0 its
    inverse and detM0 > 0 by positive definiteness. The inverse of an
    automorphism X is (M0 * X^T * adjM0) / detM0: every intermediate is
    bounded by the certificate of ps_certify_int64_kernel, while the
    adjugate route of Inverse(X) goes through minors of X that are
    exponentially larger than its entries and overflow a machine word.
   */
  MyMatrix<Tint> adjM0;
  Tint detM0;
  // The basis order per and the fingerprint diagonal (candidate count
  // per level under the identity prefix); std_basis[i] is the signed
  // index of the chain-position-i basis vector in the family.
  std::vector<int> per;
  std::vector<int> fp_diagonal;
  std::vector<int> std_basis;
  // g[i]: generators fixing e_{per[0]}, ..., e_{per[i-1]}; the first
  // nsg[i] of them are Schreier elements from ps_stab, the others come
  // from the backtrack. orders[i] is the orbit length of e_{per[i]}
  // under <g[i], ..., g[n-1]>, so the group order is the product.
  std::vector<std::vector<MyMatrix<Tint>>> g;
  std::vector<int> nsg;
  std::vector<int> orders;
  // The vector-sum invariant: depth 0 disables it, and vs_data[I] holds
  // the identity-prefix data of the level with I images assigned, for
  // I = 1, ..., n - 1.
  int depth;
  std::vector<PleskenSouvignierVectorSumLevel<Tint>> vs_data;
  // Adaptive-depth instrumentation. node_count is the number of backtrack
  // nodes in the CURRENT top-level subtree (the exploration rooted at one
  // candidate image of a basis vector); it is reset per subtree, and when
  // node_budget is nonzero and it is exceeded the search aborts
  // (PSBudgetExceeded) so the caller can raise the depth. Budgeting per
  // subtree rather than in total is what separates a depth-starved
  // explosion (one giant subtree) from a legitimately large group (many
  // small subtrees, e.g. Leech). total_nodes accumulates across subtrees
  // for diagnostics only.
  mutable size_t node_count = 0;
  mutable size_t total_nodes = 0;
  size_t node_budget = 0;
  /*
    Where a search aborted by the budget resumes. The steps before
    resume_step were completed, and what they established does not depend
    on the depth: their generators are automorphisms, and a candidate
    whose subtree was searched to the end without a leaf has no
    extension. So a higher depth resumes at resume_step with the group
    found so far, and with resume_bad, the candidates of that step already
    shown not to extend, instead of starting from -Id at step 0.
   */
  int resume_step = 0;
  std::vector<int> resume_bad;
  int n() const { return VS.n; }
};

// Thrown by the backtrack when the node budget is exceeded, signalling the
// caller to retry the search at a higher vector-sum depth.
struct PSBudgetExceeded {};

// The image of the point of signed index s under the matrix A, as a
// signed index. The family is required to be invariant, so a miss is an
// input error, not a data case.
template <typename Tint>
int ps_operate(PleskenSouvignierContext<Tint> const &ctx, int s,
               MyMatrix<Tint> const &A) {
  int j0 = (s > 0 ? s : -s) - 1;
  int n = ctx.n();
  MyVector<Tint> w(n);
  for (int i = 0; i < n; i++) {
    Tint sum(0);
    for (int k = 0; k < n; k++) {
      sum += ctx.VS.V(j0, k) * A(k, i);
    }
    w(i) = (s > 0) ? sum : -sum;
  }
  int im = ctx.VS.find_point(w);
  if (im == 0) {
    std::cerr << "PS: an automorphism maps a family vector outside the "
              << "family, so the family is not invariant\n";
    throw TerminalException{1};
  }
  return im;
}

// The scalar product of family rows j0 and k0 for the iMat-th form.
template <typename Tint>
Tint ps_scal_rows(PleskenSouvignierContext<Tint> const &ctx, int j0, int k0,
                  int iMat) {
  int n = ctx.n();
  Tint sum(0);
  for (int i = 0; i < n; i++) {
    sum += ctx.VS.V(j0, i) * ctx.W[iMat](k0, i);
  }
  return sum;
}

/*
  The count of the paper's fingerprint: the number of signed vectors
  whose norms match those of the basis vector of list index J and whose
  products with the first I chain-position basis vectors match those of
  that basis vector, both read off Fmat.
 */
template <typename Tint>
int ps_possible(PleskenSouvignierContext<Tint> const &ctx, int I, int J) {
  int m = ctx.VS.m;
  int nbMat = ctx.ListMat.size();
  int count = 0;
  for (int j = 0; j < m; j++) {
    bool good_length = true;
    for (int iMat = 0; iMat < nbMat; iMat++) {
      if (ctx.Lengths(j, iMat) != ctx.Fmat[iMat](J, J)) {
        good_length = false;
        break;
      }
    }
    if (!good_length) {
      continue;
    }
    bool okp = true;
    bool okm = true;
    for (int iMat = 0; iMat < nbMat && (okp || okm); iMat++) {
      for (int k = 0; k < I; k++) {
        Tint sc = ps_scal_rows(ctx, j, ctx.bas[ctx.per[k]], iMat);
        Tint const &target = ctx.Fmat[iMat](J, ctx.per[k]);
        if (sc != target) {
          okp = false;
        }
        if (sc != -target) {
          okm = false;
        }
        if (!okp && !okm) {
          break;
        }
      }
    }
    if (okp) {
      count++;
    }
    if (okm) {
      count++;
    }
  }
  return count;
}

/*
  The fingerprint of Section 4: the basis is reordered so that the
  candidate count of each level is minimal given the previous levels, and
  the diagonal counts are kept for the exact-count pruning of ps_cand.
 */
template <typename Tint>
void ps_fingerprint(PleskenSouvignierContext<Tint> &ctx) {
  int n = ctx.n();
  ctx.per.resize(n);
  for (int i = 0; i < n; i++) {
    ctx.per[i] = i;
  }
  std::vector<std::vector<int>> fp(n, std::vector<int>(n, 0));
  for (int i = 0; i < n; i++) {
    fp[0][i] = ps_possible(ctx, 0, i);
  }
  for (int i = 0; i < n - 1; i++) {
    int mini = i;
    for (int j = i + 1; j < n; j++) {
      if (fp[i][ctx.per[j]] < fp[i][ctx.per[mini]]) {
        mini = j;
      }
    }
    std::swap(ctx.per[mini], ctx.per[i]);
    for (int j = i + 1; j < n; j++) {
      fp[i + 1][ctx.per[j]] = ps_possible(ctx, i + 1, ctx.per[j]);
    }
  }
  ctx.fp_diagonal.resize(n);
  for (int i = 0; i < n; i++) {
    ctx.fp_diagonal[i] = fp[i][ctx.per[i]];
  }
}

/*
  The searched basis, as family row indices: the standard basis when the
  family contains it entirely -- it spans Z^n, so every leaf is integral
  by construction -- and otherwise n independent rows taken shortest
  first, short basis vectors having few candidate images.
 */
template <typename Tint>
std::vector<int> ps_select_basis(PleskenSouvignierContext<Tint> const &ctx) {
  using Tfield = typename overlying_field<Tint>::field_type;
  int n = ctx.n();
  int m = ctx.VS.m;
  int nbMat = ctx.ListMat.size();
  {
    std::vector<int> bas(n);
    bool has_std = true;
    for (int i = 0; i < n && has_std; i++) {
      MyVector<Tint> e = ZeroVector<Tint>(n);
      e(i) = 1;
      auto iter = ctx.VS.Map.find(e);
      if (iter == ctx.VS.Map.end()) {
        has_std = false;
      } else {
        bas[i] = iter->second;
      }
    }
    if (has_std) {
      return bas;
    }
  }
  if constexpr (std::is_same_v<Tint, int64_t>) {
    // The certificate admits int64 only when the family contains the
    // standard basis: the rational rank tool and the adjugate descent
    // below have intermediates no machine word holds.
    std::cerr << "PS: an int64 run without the standard basis in the "
              << "family\n";
    throw TerminalException{1};
  }
  std::vector<int> order(m);
  for (int j = 0; j < m; j++) {
    order[j] = j;
  }
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) -> bool {
    for (int iMat = 0; iMat < nbMat; iMat++) {
      if (ctx.Lengths(a, iMat) != ctx.Lengths(b, iMat)) {
        return ctx.Lengths(a, iMat) < ctx.Lengths(b, iMat);
      }
    }
    return a < b;
  });
  RankTool<Tfield> tool(n);
  std::vector<int> bas;
  MyVector<Tfield> Vf(n);
  for (int idx = 0; idx < m && static_cast<int>(bas.size()) < n; idx++) {
    int j = order[idx];
    for (int i = 0; i < n; i++) {
      Vf(i) = UniversalScalarConversion<Tfield, Tint>(ctx.VS.V(j, i));
    }
    int rank_prev = tool.get_rank();
    tool.insert_if_indep(Vf);
    if (tool.get_rank() != rank_prev) {
      bas.push_back(j);
    }
  }
  if (static_cast<int>(bas.size()) != n) {
    std::cerr << "PS: the vector family does not have full rank\n";
    throw TerminalException{1};
  }
  /*
    The greedy basis can span a proper sublattice even when the family
    contains unimodular subsets, and every unit of index multiplies the
    rational leaves the search has to discard. A swap descent on |det|
    fixes that: by Cramer, replacing row i by the vector w changes the
    determinant to sum_k w(k) * adj(k, i), a dot product against the
    adjugate, so a full scan over (row, position) is cheap and |det|
    strictly decreases at each accepted swap.
   */
  auto get_matrix = [&]() -> MyMatrix<Tint> {
    MyMatrix<Tint> B(n, n);
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        B(i, k) = ctx.VS.V(bas[i], k);
      }
    }
    return B;
  };
  // The adjugate and determinant together, fraction free (adj = det B^{-1}),
  // instead of a field inverse scaled by a separately computed determinant.
  MyMatrix<Tint> B = get_matrix();
  std::pair<MyMatrix<Tint>, Tint> adjdet = AdjugateDeterminant(B);
  Tint abs_det = T_abs(adjdet.second);
  while (abs_det > 1) {
    MyMatrix<Tint> const &adj = adjdet.first;
    bool improved = false;
    for (int idx = 0; idx < m && !improved; idx++) {
      int j = order[idx];
      for (int i = 0; i < n && !improved; i++) {
        Tint det_new(0);
        for (int k = 0; k < n; k++) {
          det_new += ctx.VS.V(j, k) * adj(k, i);
        }
        Tint abs_new = T_abs(det_new);
        if (abs_new != 0 && abs_new < abs_det) {
          bas[i] = j;
          improved = true;
        }
      }
    }
    if (!improved) {
      break;
    }
    B = get_matrix();
    adjdet = AdjugateDeterminant(B);
    abs_det = T_abs(adjdet.second);
  }
  return bas;
}

// Sign-normalizes the vector in place so that its first non-zero entry
// is positive; returns the sign applied, 0 for the zero vector.
template <typename Tint> int ps_normalize_sign(MyVector<Tint> &v) {
  int len = v.size();
  int ipos = 0;
  while (ipos < len && v(ipos) == 0) {
    ipos++;
  }
  if (ipos == len) {
    return 0;
  }
  if (v(ipos) > 0) {
    return 1;
  }
  for (int i = ipos; i < len; i++) {
    v(i) = -v(i);
  }
  return -1;
}

// The identity-prefix vector-sum data of every level (Section 5): the
// tuples come straight out of the W tables, the products with a basis
// vector being single entries.
template <typename Tint>
void ps_init_vector_sums(PleskenSouvignierContext<Tint> &ctx, int depth,
                         [[maybe_unused]] std::ostream &os) {
  ctx.depth = depth;
  if (depth == 0) {
    return;
  }
  int n = ctx.n();
  int m = ctx.VS.m;
  int nbMat = ctx.ListMat.size();
  // Rebuild from scratch: on a depth escalation this is called a second
  // time, and the per-level tuple maps and sums must not retain the
  // previous depth's data.
  ctx.vs_data.clear();
  ctx.vs_data.resize(n);
  for (int I = 1; I < n; I++) {
    PleskenSouvignierVectorSumLevel<Tint> &lev = ctx.vs_data[I];
    int dep = std::min(depth, I);
    lev.dep = dep;
    int len = nbMat * dep;
    std::vector<MyVector<Tint>> sums;
    for (int j = 0; j < m; j++) {
      MyVector<Tint> tup(len);
      int pos = 0;
      for (int iMat = 0; iMat < nbMat; iMat++) {
        for (int k = I - dep; k < I; k++) {
          tup(pos) = ps_scal_rows(ctx, j, ctx.bas[ctx.per[k]], iMat);
          pos++;
        }
      }
      int sign = ps_normalize_sign(tup);
      int idx;
      auto iter = lev.TupleMap.find(tup);
      if (iter == lev.TupleMap.end()) {
        idx = lev.counts.size();
        lev.TupleMap[tup] = idx;
        lev.counts.push_back(0);
        sums.push_back(ZeroVector<Tint>(n));
      } else {
        idx = iter->second;
      }
      lev.counts[idx]++;
      if (sign != 0) {
        for (int i = 0; i < n; i++) {
          sums[idx](i) += Tint(sign) * ctx.VS.V(j, i);
        }
      }
    }
    int n_tuple = lev.counts.size();
    lev.SumNorm.resize(nbMat);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      lev.SumNorm[iMat].resize(n_tuple);
      for (int idx = 0; idx < n_tuple; idx++) {
        lev.SumNorm[iMat][idx] =
            EvaluationQuadForm(ctx.ListMat[iMat], sums[idx]);
      }
    }
#ifdef DEBUG_PLESKEN_SOUVIGNIER
    {
      // The rank of the lattice spanned by the vector sums measures the
      // discriminating power of the depth-dep invariant at this level;
      // AUTO's guidance is to raise depth until this reaches full rank n.
      int rk = 0;
      if (n_tuple > 0) {
        MyMatrix<Tint> SM(n_tuple, n);
        for (int idx = 0; idx < n_tuple; idx++) {
          for (int i = 0; i < n; i++) {
            SM(idx, i) = sums[idx](i);
          }
        }
        rk = RankMat(SM);
      }
      os << "PS vsrank: level " << I << " dep " << dep << " n_tuple "
         << n_tuple << " rank " << rk << " / " << n << "\n";
    }
#endif
  }
}

/*
  The runtime side of the vector-sum invariant, at the level with I
  images assigned: the tuple of every family vector must hit the
  identity classes with the same multiplicities, and the per-class sums
  must have the identity norms. Sound because an automorphism extending
  the prefix maps the identity classes bijectively onto the runtime
  ones, tuple value to tuple value, and the class sums accordingly.
 */
template <typename Tint>
bool ps_vector_sum_check(PleskenSouvignierContext<Tint> const &Ci,
                         PleskenSouvignierContext<Tint> const &Co, int I,
                         std::vector<int> const &x) {
  PleskenSouvignierVectorSumLevel<Tint> const &lev = Ci.vs_data[I];
  int n = Ci.n();
  int m = Co.VS.m;
  int nbMat = Ci.ListMat.size();
  int dep = lev.dep;
  int len = nbMat * dep;
  int n_tuple = lev.counts.size();
  std::vector<int> counts_run(n_tuple, 0);
  std::vector<MyVector<Tint>> sums_run(n_tuple, ZeroVector<Tint>(n));
  MyVector<Tint> tup(len);
  for (int j = 0; j < m; j++) {
    int pos = 0;
    for (int iMat = 0; iMat < nbMat; iMat++) {
      for (int k = I - dep; k < I; k++) {
        int xk = x[k];
        int row = (xk > 0 ? xk : -xk) - 1;
        Tint sc(0);
        for (int i = 0; i < n; i++) {
          sc += Co.VS.V(j, i) * Co.W[iMat](row, i);
        }
        if (xk < 0) {
          sc = -sc;
        }
        tup(pos) = sc;
        pos++;
      }
    }
    int sign = ps_normalize_sign(tup);
    auto iter = lev.TupleMap.find(tup);
    if (iter == lev.TupleMap.end()) {
      return false;
    }
    int idx = iter->second;
    counts_run[idx]++;
    if (counts_run[idx] > lev.counts[idx]) {
      return false;
    }
    if (sign != 0) {
      for (int i = 0; i < n; i++) {
        sums_run[idx](i) += Tint(sign) * Co.VS.V(j, i);
      }
    }
  }
  // Every count is bounded by its identity value and the totals agree,
  // so the multisets agree.
  for (int iMat = 0; iMat < nbMat; iMat++) {
    for (int idx = 0; idx < n_tuple; idx++) {
      Tint nrm = EvaluationQuadForm(Co.ListMat[iMat], sums_run[idx]);
      if (nrm != lev.SumNorm[iMat][idx]) {
        return false;
      }
    }
  }
  return true;
}

/*
  The context. with_fingerprint is set for the lattice whose basis is
  searched (the only one for automorphisms, the first one for an
  isometry); the fingerprint, the location of the standard basis and the
  vector-sum data are only needed there. depth = -1 selects the default
  of n/10 rounded, as Hecke does; depth = 0 disables the vector sums.
 */
/*
  The searched basis of a context and the inverse data of its
  position-ordered matrix (see the context). Computed once over the big
  integers by a context built there, it lets an int64 context skip
  ps_select_basis and the adjugate of the basis, whose intermediates do
  not fit a machine word, while every quantity of the search itself is
  covered by the certificate.
 */
template <typename Tint> struct PSBasisData {
  std::vector<int> bas;
  MyMatrix<Tint> Badj;
  Tint Bden;
};

template <typename Tint>
PSBasisData<Tint> ps_basis_data(PleskenSouvignierContext<Tint> const &ctx) {
  return {ctx.bas, ctx.Badj, ctx.Bden};
}

template <typename Tint>
PleskenSouvignierContext<Tint> PleskenSouvignierBuildContext(
    std::vector<MyMatrix<Tint>> const &ListMat, MyMatrix<Tint> const &SHVhalf,
    bool with_fingerprint, int depth, [[maybe_unused]] std::ostream &os,
    std::optional<std::pair<MyMatrix<Tint>, Tint>> const &adjdet0 = {},
    std::optional<PSBasisData<Tint>> const &basis_data = {}) {
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  MicrosecondTime time;
#endif
  PleskenSouvignierContext<Tint> ctx;
  ctx.ListMat = ListMat;
  if (adjdet0) {
    ctx.adjM0 = adjdet0->first;
    ctx.detM0 = adjdet0->second;
  } else {
    if constexpr (std::is_same_v<Tint, int64_t>) {
      // The Bareiss minors of the adjugate overflow a machine word; an
      // int64 context has to receive the adjugate from the caller, who
      // computed it over the big integers.
      std::cerr << "PS: an int64 context needs the precomputed adjugate of "
                << "ListMat[0]\n";
      throw TerminalException{1};
    } else {
      std::pair<MyMatrix<Tint>, Tint> pr = AdjugateDeterminant(ListMat[0]);
      ctx.adjM0 = std::move(pr.first);
      ctx.detM0 = pr.second;
    }
  }
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  if (ctx.detM0 <= 0) {
    std::cerr << "PS: ListMat[0] has to be positive definite, so its "
              << "determinant has to be positive\n";
    throw TerminalException{1};
  }
#endif
  ctx.VS = PleskenSouvignierBuildVectorSet(SHVhalf);
  int n = ctx.VS.n;
  int m = ctx.VS.m;
  int nbMat = ListMat.size();
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  for (auto &eMat : ListMat) {
    if (!IsSymmetricMatrix(eMat)) {
      std::cerr << "PS: the matrices of the configuration have to be "
                << "symmetric\n";
      throw TerminalException{1};
    }
  }
#endif
  ctx.W.reserve(nbMat);
  for (int iMat = 0; iMat < nbMat; iMat++) {
    ctx.W.push_back(ctx.VS.V * ListMat[iMat]);
  }
  ctx.Lengths = MyMatrix<Tint>(m, nbMat);
  for (int j = 0; j < m; j++) {
    for (int iMat = 0; iMat < nbMat; iMat++) {
      Tint sum(0);
      for (int i = 0; i < n; i++) {
        sum += ctx.VS.V(j, i) * ctx.W[iMat](j, i);
      }
      ctx.Lengths(j, iMat) = sum;
    }
  }
  ctx.g.assign(n, {});
  ctx.nsg.assign(n, 0);
  ctx.orders.assign(n, 1);
  ctx.depth = 0;
  // -Id preserves every symmetric form and every antipodal family.
  ctx.g[0].push_back(-IdentityMat<Tint>(n));
  if (with_fingerprint) {
    ctx.bas = basis_data ? basis_data->bas : ps_select_basis(ctx);
    ctx.Fmat.clear();
    for (int iMat = 0; iMat < nbMat; iMat++) {
      MyMatrix<Tint> F(n, n);
      for (int a = 0; a < n; a++) {
        for (int b = 0; b < n; b++) {
          F(a, b) = ps_scal_rows(ctx, ctx.bas[a], ctx.bas[b], iMat);
        }
      }
      ctx.Fmat.push_back(F);
    }
    ps_fingerprint(ctx);
    ctx.std_basis.resize(n);
    for (int i = 0; i < n; i++) {
      // The basis vectors are family rows, which are stored as the
      // sign-normalized representatives: their signed index is positive.
      ctx.std_basis[i] = ctx.bas[ctx.per[i]] + 1;
    }
    // The inverse data of the position-ordered basis matrix, through
    // the adjugate so that the leaves stay in integer arithmetic.
    MyMatrix<Tint> B(n, n);
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        B(i, k) = ctx.VS.V(ctx.bas[ctx.per[i]], k);
      }
    }
    // The adjugate and determinant of the basis, fraction free
    // (adj = det B^{-1}), so the leaves stay in integer arithmetic with
    // B^{-1} = Badj / Bden and Bden > 0. The fingerprint order per is the
    // same as in the context that computed basis_data, the data being the
    // same, so its inverse data applies as is.
    // In int64 without basis_data the basis is the standard one (the
    // certificate requires it), a signed permutation matrix, whose
    // adjugate stays small.
    if (basis_data) {
      ctx.Badj = basis_data->Badj;
      ctx.Bden = basis_data->Bden;
    } else {
      std::pair<MyMatrix<Tint>, Tint> adjdet = AdjugateDeterminant(B);
      ctx.Badj = std::move(adjdet.first);
      Tint det = adjdet.second;
      if (det < 0) {
        det = -det;
        ctx.Badj = -ctx.Badj;
      }
      ctx.Bden = det;
    }
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
    if (B * ctx.Badj != ctx.Bden * IdentityMat<Tint>(n)) {
      std::cerr << "PS: the inverse data does not invert the basis\n";
      throw TerminalException{1};
    }
#endif
#ifdef DEBUG_PLESKEN_SOUVIGNIER
    if (ctx.Bden != 1) {
      os << "PS: the searched basis spans a sublattice of index " << ctx.Bden
         << ", rational leaves will be discarded\n";
    }
#endif
    if (depth == -1) {
      // Hecke's default round(n / 10): 1 from dimension 5 up, 2 from 15.
      depth = (n + 5) / 10;
    }
    ps_init_vector_sums(ctx, depth, os);
  }
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  os << "|PS: BuildContext|=" << time << "\n";
#endif
  return ctx;
}

/*
  The candidates for the image of the level-I basis vector, given the
  images x[0], ..., x[I-1] already chosen: the signed vectors whose norms
  match e_{per[I]} and whose products with the chosen images match the
  Gram entries. Returns false unless the count is EXACTLY the fingerprint
  value: an automorphism maps the candidate set of the identity prefix
  bijectively to this one, so any other count kills the branch (Section 4
  of the paper).

  Two contexts so that the same function serves the isometry test: the
  basis, fingerprint and target products come from Ci, the vectors from
  Co. For automorphisms both are the same context.
 */
template <typename Tint>
bool ps_cand(PleskenSouvignierContext<Tint> const &Ci,
             PleskenSouvignierContext<Tint> const &Co, int I,
             std::vector<int> const &x, std::vector<int> &candidates) {
  int m = Co.VS.m;
  int nbMat = Ci.ListMat.size();
  int target = Ci.fp_diagonal[I];
  candidates.clear();
  // The vector-sum invariant kills most of the prefixes that the
  // count of candidates alone cannot.
  if (I >= 1 && Ci.depth > 0) {
    if (!ps_vector_sum_check(Ci, Co, I, x)) {
      return false;
    }
  }
  for (int j = 0; j < m; j++) {
    bool okp = true;
    bool okm = true;
    for (int iMat = 0; iMat < nbMat; iMat++) {
      if (Co.Lengths(j, iMat) != Ci.Fmat[iMat](Ci.per[I], Ci.per[I])) {
        okp = false;
        okm = false;
        break;
      }
      for (int k = 0; k < I; k++) {
        int xk = x[k];
        int row = (xk > 0 ? xk : -xk) - 1;
        int n = Ci.n();
        Tint sc(0);
        for (int i = 0; i < n; i++) {
          sc += Co.VS.V(j, i) * Co.W[iMat](row, i);
        }
        if (xk < 0) {
          sc = -sc;
        }
        Tint const &tgt = Ci.Fmat[iMat](Ci.per[I], Ci.per[k]);
        if (sc != tgt) {
          okp = false;
        }
        if (sc != -tgt) {
          okm = false;
        }
        if (!okp && !okm) {
          break;
        }
      }
      if (!okp && !okm) {
        break;
      }
    }
    if (okp) {
      if (static_cast<int>(candidates.size()) >= target) {
        return false;
      }
      candidates.push_back(j + 1);
    }
    if (okm) {
      if (static_cast<int>(candidates.size()) >= target) {
        return false;
      }
      candidates.push_back(-(j + 1));
    }
  }
  return static_cast<int>(candidates.size()) == target;
}

/*
  The matrix of a complete assignment: the solution X of B * X = A, with
  B the position-ordered basis and A the stacked images, computed as
  (Badj * A) / Bden in integer arithmetic. Nothing when X is not
  integral, which only happens when the basis spans a proper sublattice
  and the assignment is a rational solution to be discarded.
 */
template <typename Tint>
std::optional<MyMatrix<Tint>>
ps_matgen_opt(PleskenSouvignierContext<Tint> const &Ci,
              PleskenSouvignierContext<Tint> const &Co,
              std::vector<int> const &x) {
  int n = Ci.n();
  MyMatrix<Tint> A(n, n);
  for (int i = 0; i < n; i++) {
    int xi = x[i];
    int row = (xi > 0 ? xi : -xi) - 1;
    for (int k = 0; k < n; k++) {
      A(i, k) = (xi > 0) ? Co.VS.V(row, k) : -Co.VS.V(row, k);
    }
  }
  MyMatrix<Tint> num = Ci.Badj * A;
  if (Ci.Bden == 1) {
    return num;
  }
  MyMatrix<Tint> X(n, n);
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < n; k++) {
      if (ResInt(num(i, k), Ci.Bden) != 0) {
        return {};
      }
      X(i, k) = QuoInt(num(i, k), Ci.Bden);
    }
  }
  return X;
}

// The same for the callers whose assignment comes from established
// automorphisms, where a rational solution is a programming error.
template <typename Tint>
MyMatrix<Tint> ps_matgen(PleskenSouvignierContext<Tint> const &Ci,
                         PleskenSouvignierContext<Tint> const &Co,
                         std::vector<int> const &x) {
  std::optional<MyMatrix<Tint>> opt = ps_matgen_opt(Ci, Co, x);
  return unfold_opt(opt, "the assignment should give an integral matrix");
}

/*
  The inverse of an isometry Q with Q * M2 * Q^T = M1: from
  M2 * Q^T = Q^{-1} * M1 one gets Q^{-1} = (M2 * Q^T * adjM1) / detM1,
  and the division is exact. For an automorphism M1 = M2. Unlike
  Inverse(Q), whose adjugate goes through minors of Q exponentially
  larger than its entries, the intermediates here are bounded by
  n^2 * |M2| * |adjM1| * |Q|, which the int64 certificate covers.
 */
template <typename Tint>
MyMatrix<Tint> ps_isometry_inverse(MyMatrix<Tint> const &M2,
                                   MyMatrix<Tint> const &adjM1,
                                   Tint const &detM1, MyMatrix<Tint> const &Q) {
  int n = Q.rows();
  MyMatrix<Tint> Z = M2 * Q.transpose() * adjM1;
  MyMatrix<Tint> Qinv(n, n);
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < n; k++) {
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
      if (ResInt(Z(i, k), detM1) != 0) {
        std::cerr << "PS: the isometry inverse is not integral\n";
        throw TerminalException{1};
      }
#endif
      Qinv(i, k) = QuoInt(Z(i, k), detM1);
    }
  }
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  if (Q * Qinv != IdentityMat<Tint>(n)) {
    std::cerr << "PS: ps_isometry_inverse did not invert the matrix\n";
    throw TerminalException{1};
  }
#endif
  return Qinv;
}

// The inverse of an automorphism of the context's configuration.
template <typename Tint>
MyMatrix<Tint> ps_aut_inverse(PleskenSouvignierContext<Tint> const &ctx,
                              MyMatrix<Tint> const &X) {
  return ps_isometry_inverse(ctx.ListMat[0], ctx.adjM0, ctx.detM0, X);
}

// The orbit of a set of signed points under a list of matrices.
template <typename Tint>
std::vector<int> ps_orbit(PleskenSouvignierContext<Tint> const &ctx,
                          std::vector<int> const &pts,
                          std::vector<MyMatrix<Tint>> const &ListGen) {
  int m = ctx.VS.m;
  std::vector<uint8_t> flag(2 * m + 1, 0);
  std::vector<int> orb;
  for (auto &pt : pts) {
    if (flag[pt + m] == 0) {
      flag[pt + m] = 1;
      orb.push_back(pt);
    }
  }
  size_t cnd = 0;
  while (cnd < orb.size()) {
    for (auto &eGen : ListGen) {
      int im = ps_operate(ctx, orb[cnd], eGen);
      if (flag[im + m] == 0) {
        flag[im + m] = 1;
        orb.push_back(im);
      }
    }
    cnd++;
  }
  return orb;
}

// The orbit length of one point, computed no further than cap: the
// caller only compares it with cap.
template <typename Tint>
int ps_orbitlen(PleskenSouvignierContext<Tint> const &ctx, int pt, int cap,
                std::vector<MyMatrix<Tint>> const &ListGen) {
  int m = ctx.VS.m;
  std::vector<uint8_t> flag(2 * m + 1, 0);
  std::vector<int> orb{pt};
  flag[pt + m] = 1;
  size_t cnd = 0;
  while (cnd < orb.size() && static_cast<int>(orb.size()) < cap) {
    for (auto &eGen : ListGen) {
      int im = ps_operate(ctx, orb[cnd], eGen);
      if (flag[im + m] == 0) {
        flag[im + m] = 1;
        orb.push_back(im);
        if (static_cast<int>(orb.size()) >= cap) {
          break;
        }
      }
    }
    cnd++;
  }
  return orb.size();
}

// Removes the members of orb from the ordered candidate list.
inline void ps_setdiff(std::vector<int> &candidates,
                       std::vector<int> const &orb, int m) {
  std::vector<uint8_t> flag(2 * m + 1, 0);
  for (auto &pt : orb) {
    flag[pt + m] = 1;
  }
  std::vector<int> reduced;
  reduced.reserve(candidates.size());
  for (auto &c : candidates) {
    if (flag[c + m] == 0) {
      reduced.push_back(c);
    }
  }
  candidates = std::move(reduced);
}

/*
  The element X1 * G * X2^{-1}: x1 is the assignment of an element X1
  mapping the basis prefix somewhere, x2 of an element X2 with the same
  image of e_{per[I]}, and G a generator mapping the image of X1 to the
  image of X2, so the product stabilizes every basis vector on which the
  two words agree. X2 is in GL_n(Z), so the inverse is integral.
 */
template <typename Tint>
MyMatrix<Tint> ps_stabil(PleskenSouvignierContext<Tint> const &ctx,
                         std::vector<int> const &x1,
                         std::vector<int> const &x2,
                         MyMatrix<Tint> const &G) {
  int n = ctx.n();
  std::vector<int> x(n);
  for (int i = 0; i < n; i++) {
    x[i] = ps_operate(ctx, x1[i], G);
  }
  MyMatrix<Tint> XG = ps_matgen(ctx, ctx, x);
  MyMatrix<Tint> X2 = ps_matgen(ctx, ctx, x2);
  MyMatrix<Tint> S = XG * ps_aut_inverse(ctx, X2);
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  if (S * X2 != XG) {
    std::cerr << "PS: the stabilizer element does not solve S * X2 = XG\n";
    throw TerminalException{1};
  }
#endif
  return S;
}

// A deterministic xorshift, standing in for the random prodding of the
// original: the repository is deterministic throughout, and canonical
// runs must not depend on a seed.
struct PSDeterministicRng {
  uint64_t state = UINT64_C(0x9E3779B97F4A7C15);
  uint64_t next() {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
  }
  int next_index(int siz) {
    return static_cast<int>(next() % static_cast<uint64_t>(siz));
  }
};

/*
  Schreier harvesting along the orbit of e_{per[I]} (Section 8 and the
  stab routine of AUTO): while computing the orbit under the generators
  fixing the earlier basis vectors, every coincidence of two words
  yields an element fixing e_{per[0]}, ..., e_{per[I]}, which is kept
  when it enlarges the orbit of a deeper basis vector. These elements
  are much cheaper than ones found by backtracking, and the heuristic
  budget (Maxfail, Rest) bounds the redundancy; the budget only affects
  how many generators are harvested, never correctness.

  Deviation from the Hecke port: the level j at which a harvested
  element is tested is the first index where the two words differ, per
  Souvignier's original code; the port breaks on the first index where
  they agree, which is always I itself.
 */
template <typename Tint>
void ps_stab(PleskenSouvignierContext<Tint> &ctx, int I,
             [[maybe_unused]] std::ostream &os) {
  int n = ctx.n();
  int m = ctx.VS.m;
  int Rest = 0;
  for (int i = I; i < n; i++) {
    if (ctx.fp_diagonal[i] > 1 && ctx.orders[i] < ctx.fp_diagonal[i]) {
      Rest++;
    }
  }
  int Maxfail = Rest;
  for (int i = 0; i < n; i++) {
    if (ctx.fp_diagonal[i] > 1) {
      Maxfail++;
    }
  }
  std::vector<MyMatrix<Tint>> H;
  for (int i = I; i < n; i++) {
    for (auto &eGen : ctx.g[i]) {
      H.push_back(eGen);
    }
  }
  if (H.empty()) {
    return;
  }
  std::vector<std::vector<int>> w(2 * m + 1);
  std::vector<uint8_t> flag(2 * m + 1, 0);
  std::vector<int> orb;
  orb.push_back(ctx.std_basis[I]);
  flag[orb[0] + m] = 1;
  w[orb[0] + m] = ctx.std_basis;
  size_t cnd = 0;
  int fail = 0;
  PSDeterministicRng rng;
  while (cnd < orb.size() && fail < Maxfail + Rest) {
    size_t nH = H.size();
    for (size_t i = 0; i < nH; i++) {
      if (fail >= Maxfail + Rest) {
        break;
      }
      size_t use_cnd = cnd;
      size_t use_i = i;
      if (fail >= Maxfail) {
        // The regular harvest failed Maxfail times: Rest more elements
        // are taken by prodding arbitrary orbit points, which often
        // helps the deeper levels.
        use_cnd = rng.next_index(orb.size());
        use_i = rng.next_index(nH);
      }
      int im = ps_operate(ctx, orb[use_cnd], H[use_i]);
      if (flag[im + m] == 0) {
        flag[im + m] = 1;
        std::vector<int> wim(n);
        for (int j = 0; j < n; j++) {
          wim[j] = ps_operate(ctx, w[orb[use_cnd] + m][j], H[use_i]);
        }
        w[im + m] = std::move(wim);
        orb.push_back(im);
      } else {
        // A coincidence: the first level where the two words differ is
        // where the Schreier element can act.
        int j = I;
        while (j < n) {
          if (ps_operate(ctx, w[orb[use_cnd] + m][j], H[use_i]) !=
              w[im + m][j]) {
            break;
          }
          j++;
        }
        if (j < n && (ctx.orders[j] < ctx.fp_diagonal[j] || fail >= Maxfail)) {
          std::vector<int> x1(n);
          for (int k = 0; k < n; k++) {
            x1[k] = w[orb[use_cnd] + m][k];
          }
          MyMatrix<Tint> S = ps_stabil(ctx, x1, w[im + m], H[use_i]);
          std::vector<MyMatrix<Tint>> Hj{S};
          for (int k = j; k < n; k++) {
            for (auto &eGen : ctx.g[k]) {
              Hj.push_back(eGen);
            }
          }
          int tmplen =
              ps_orbitlen(ctx, ctx.std_basis[j], ctx.fp_diagonal[j], Hj);
          if (tmplen > ctx.orders[j] || fail >= Maxfail) {
            ctx.orders[j] = tmplen;
            ctx.nsg[j]++;
            ctx.g[j].insert(ctx.g[j].begin() + (ctx.nsg[j] - 1), S);
            H.push_back(S);
            if (fail < Maxfail) {
              fail = 0;
            } else {
              fail++;
            }
          } else {
            fail++;
          }
        } else {
          fail++;
        }
      }
    }
    if (fail < Maxfail) {
      cnd++;
    }
  }
}

// The backtrack extension: whether the partial assignment x[0..step-1]
// with the candidate list for level step extends to a complete
// automorphism; on success x holds it.
template <typename Tint>
bool ps_aut_extend(PleskenSouvignierContext<Tint> const &ctx, int step,
                   std::vector<int> &x, std::vector<int> const &cand_step) {
  int n = ctx.n();
  for (auto &c : cand_step) {
    ctx.node_count++;
    ctx.total_nodes++;
    if (ctx.node_budget != 0 && ctx.node_count > ctx.node_budget) {
      throw PSBudgetExceeded{};
    }
    x[step] = c;
    if (step == n - 1) {
      // A complete product-matching assignment; with a unimodular basis
      // it is integral by construction, otherwise the rational
      // solutions are discarded here.
      if (ctx.Bden == 1 || ps_matgen_opt(ctx, ctx, x).has_value()) {
        return true;
      }
      continue;
    }
    std::vector<int> cand_next;
    if (ps_cand(ctx, ctx, step + 1, x, cand_next)) {
      if (ps_aut_extend(ctx, step + 1, x, cand_next)) {
        return true;
      }
    }
  }
  return false;
}

/*
  The automorphism search along the stabilizer chain of the basis
  (Section 8): at level step the group G_step fixing the earlier basis
  vectors is assembled by trying, for each candidate image of e_{per[step]}
  outside the orbit of the group found so far, to extend it to a full
  automorphism. Success gives a new generator and the orbit swallows a
  chunk of candidates; failure kills the whole orbit of the candidate,
  a group element mapping one dead end to another dead end.
 */
template <typename Tint>
void ps_auto(PleskenSouvignierContext<Tint> &ctx, std::ostream &os) {
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  MicrosecondTime time;
#endif
  int n = ctx.n();
  int m = ctx.VS.m;
  for (int step = ctx.resume_step; step < n; step++) {
    std::vector<MyMatrix<Tint>> H;
    for (int i = step; i < n; i++) {
      for (auto &eGen : ctx.g[i]) {
        H.push_back(eGen);
      }
    }
    std::vector<int> x(n, 0);
    for (int i = 0; i < step; i++) {
      x[i] = ctx.std_basis[i];
    }
    std::vector<int> cand_step;
    if (ctx.fp_diagonal[step] > 1) {
      [[maybe_unused]] bool test = ps_cand(ctx, ctx, step, x, cand_step);
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
      if (!test) {
        std::cerr << "PS: the identity prefix has a wrong candidate count, "
                  << "the fingerprint is inconsistent\n";
        throw TerminalException{1};
      }
#endif
    } else {
      cand_step.push_back(ctx.std_basis[step]);
    }
    std::vector<int> orb = ps_orbit(ctx, {ctx.std_basis[step]}, H);
    ctx.orders[step] = orb.size();
    ps_setdiff(cand_step, orb, m);
    std::vector<int> bad;
    if (step == ctx.resume_step && !ctx.resume_bad.empty()) {
      bad = ctx.resume_bad;
      std::vector<int> oc = ps_orbit(ctx, bad, H);
      ps_setdiff(cand_step, oc, m);
    }
    ctx.resume_step = step;
    ctx.resume_bad = bad;
#ifdef DEBUG_PLESKEN_SOUVIGNIER
    os << "PS: step " << step << ", fingerprint " << ctx.fp_diagonal[step]
       << ", candidates outside the orbit " << cand_step.size() << "\n";
#endif
    while (!cand_step.empty()) {
      int im = cand_step[0];
      x[step] = im;
      for (int i = step + 1; i < n; i++) {
        x[i] = 0;
      }
      // The budget is per top-level subtree: reset before exploring this
      // candidate so that a legitimately large group (many small subtrees)
      // is not mistaken for a depth-starved explosion (one huge subtree).
      ctx.node_count = 0;
      bool found = false;
      if (step < n - 1) {
        std::vector<int> cand_next;
        if (ps_cand(ctx, ctx, step + 1, x, cand_next)) {
          found = ps_aut_extend(ctx, step + 1, x, cand_next);
        }
      } else {
        found = (ctx.Bden == 1) || ps_matgen_opt(ctx, ctx, x).has_value();
      }
      if (!found) {
        std::vector<int> oc = ps_orbit(ctx, {im}, H);
        ps_setdiff(cand_step, oc, m);
        bad.push_back(im);
        ctx.resume_bad = bad;
      } else {
        ctx.g[step].push_back(ps_matgen(ctx, ctx, x));
        H.push_back(ctx.g[step].back());
        orb = ps_orbit(ctx, {ctx.std_basis[step]}, H);
        ctx.orders[step] = orb.size();
        ps_setdiff(cand_step, orb, m);
        if (!bad.empty()) {
          std::vector<int> oc = ps_orbit(ctx, bad, H);
          ps_setdiff(cand_step, oc, m);
        }
#ifdef DEBUG_PLESKEN_SOUVIGNIER
        os << "PS: step " << step << ", new generator, orbit "
           << ctx.orders[step] << ", candidates left " << cand_step.size()
           << "\n";
#endif
      }
    }
    if (step < n - 1 && ctx.orders[step] > 1) {
      ps_stab(ctx, step, os);
    }
    ctx.resume_step = step + 1;
    ctx.resume_bad.clear();
  }
#ifdef DEBUG_PLESKEN_SOUVIGNIER
  os << "PS: auto nodes=" << ctx.total_nodes << " depth=" << ctx.depth << "\n";
#endif
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  os << "|PS: auto|=" << time << "\n";
#endif
}

/*
  The result of the automorphism computation. The group order is the
  product of the orbit lengths along the stabilizer chain.
 */
template <typename Tint> struct PleskenSouvignierAutomResult {
  std::vector<MyMatrix<Tint>> ListGen;
  std::vector<int> ListOrbitSize;
  mpz_class order;
  // The vector-sum depth the adaptive search settled on and the number of
  // backtrack nodes it visited (diagnostics for tuning the depth rule).
  int depth_used = 0;
  size_t nodes_used = 0;
};

inline mpz_class
PleskenSouvignierGroupOrder(std::vector<int> const &ListOrbitSize) {
  mpz_class order(1);
  for (auto &siz : ListOrbitSize) {
    order *= siz;
  }
  return order;
}

// The conversions of the int64 fast path. The downward one is reached
// only behind the certificate, whose bounds dominate every entry.
template <typename Tint>
MyMatrix<int64_t> ps_matrix_to_int64(MyMatrix<Tint> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  MyMatrix<int64_t> Ms(nbRow, nbCol);
  for (int i = 0; i < nbRow; i++) {
    for (int j = 0; j < nbCol; j++) {
      Ms(i, j) = UniversalScalarConversion<int64_t, Tint>(M(i, j));
    }
  }
  return Ms;
}

template <typename Tint>
MyMatrix<Tint> ps_matrix_from_int64(MyMatrix<int64_t> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  MyMatrix<Tint> Mt(nbRow, nbCol);
  for (int i = 0; i < nbRow; i++) {
    for (int j = 0; j < nbCol; j++) {
      Mt(i, j) = UniversalScalarConversion<Tint, int64_t>(M(i, j));
    }
  }
  return Mt;
}

template <typename Tint>
PSBasisData<int64_t> ps_basis_data_to_int64(PSBasisData<Tint> const &bd) {
  return {bd.bas, ps_matrix_to_int64(bd.Badj),
          UniversalScalarConversion<int64_t, Tint>(bd.Bden)};
}

/*
  The measured maxima from which the int64 certificate is computed, all
  in the big type of the caller. has_std_basis records whether the family
  contains the whole standard basis: then the searched basis is a
  permutation matrix, its adjugate is free of blow-up and every leaf is
  integral, which is what the int64 engine requires.
 */
template <typename Tint> struct PSInt64Bounds {
  Tint n;
  Tint m;
  Tint Vmax;
  Tint Mmax;
  Tint maxMii;
  Tint adjMmax;
  Tint detM0;
  bool has_std_basis;
  // Set when the searched basis comes as precomputed PSBasisData, BadjMax
  // being the largest entry of its adjugate: the leaves then compute
  // Badj * A, bounded by E9 below.
  bool has_basis_data = false;
  Tint BadjMax = Tint(0);
};

template <typename Tint>
PSInt64Bounds<Tint>
ps_int64_bounds(std::vector<MyMatrix<Tint>> const &ListMat,
                MyMatrix<Tint> const &SHVhalf,
                std::pair<MyMatrix<Tint>, Tint> const &adjdet0) {
  int n = SHVhalf.cols();
  int m = SHVhalf.rows();
  PSInt64Bounds<Tint> b;
  b.n = Tint(n);
  b.m = Tint(m);
  b.Vmax = Tint(1);
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      Tint val = T_abs(SHVhalf(j, i));
      if (val > b.Vmax) {
        b.Vmax = val;
      }
    }
  }
  b.Mmax = Tint(1);
  for (auto &eMat : ListMat) {
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        Tint val = T_abs(eMat(i, k));
        if (val > b.Mmax) {
          b.Mmax = val;
        }
      }
    }
  }
  b.maxMii = Tint(1);
  for (int i = 0; i < n; i++) {
    if (ListMat[0](i, i) > b.maxMii) {
      b.maxMii = ListMat[0](i, i);
    }
  }
  b.adjMmax = Tint(1);
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < n; k++) {
      Tint val = T_abs(adjdet0.first(i, k));
      if (val > b.adjMmax) {
        b.adjMmax = val;
      }
    }
  }
  b.detM0 = adjdet0.second;
  b.has_std_basis = true;
  std::vector<uint8_t> found(n, 0);
  for (int j = 0; j < m; j++) {
    int nb_nonzero = 0;
    int pos = -1;
    for (int i = 0; i < n; i++) {
      if (SHVhalf(j, i) != 0) {
        nb_nonzero++;
        pos = i;
      }
    }
    if (nb_nonzero == 1 && T_abs(SHVhalf(j, pos)) == 1) {
      found[pos] = 1;
    }
  }
  for (int i = 0; i < n; i++) {
    if (found[i] == 0) {
      b.has_std_basis = false;
    }
  }
  return b;
}

/*
  When the family lacks the standard basis, the basis of the search is
  selected once over the big integers, by a context built without vector
  sums, and its data added to the bounds; the int64 search then reuses it.
 */
template <typename Tint>
std::optional<PSBasisData<Tint>>
ps_basis_data_for_int64(std::vector<MyMatrix<Tint>> const &ListMat,
                        MyMatrix<Tint> const &SHVhalf,
                        std::pair<MyMatrix<Tint>, Tint> const &adjdet0,
                        PSInt64Bounds<Tint> &bounds, std::ostream &os) {
  if (bounds.has_std_basis) {
    return {};
  }
  PleskenSouvignierContext<Tint> ctx = PleskenSouvignierBuildContext(
      ListMat, SHVhalf, true, 0, os,
      std::optional<std::pair<MyMatrix<Tint>, Tint>>(adjdet0));
  PSBasisData<Tint> bd = ps_basis_data(ctx);
  bounds.has_basis_data = true;
  bounds.BadjMax = Tint(1);
  for (int i = 0; i < bd.Badj.rows(); i++) {
    for (int k = 0; k < bd.Badj.cols(); k++) {
      Tint val = T_abs(bd.Badj(i, k));
      if (val > bounds.BadjMax) {
        bounds.BadjMax = val;
      }
    }
  }
  if (bd.Bden > bounds.BadjMax) {
    bounds.BadjMax = bd.Bden;
  }
  return bd;
}

// The merge for an isometry test, where the one search touches the
// vectors and automorphisms of both configurations. The determinants
// have to agree for an isometry to exist, which the caller has checked.
template <typename Tint>
PSInt64Bounds<Tint> ps_int64_bounds_merge(PSInt64Bounds<Tint> const &b1,
                                          PSInt64Bounds<Tint> const &b2) {
  PSInt64Bounds<Tint> b = b1;
  auto f_max = [](Tint &tgt, Tint const &other) -> void {
    if (other > tgt) {
      tgt = other;
    }
  };
  f_max(b.m, b2.m);
  f_max(b.Vmax, b2.Vmax);
  f_max(b.Mmax, b2.Mmax);
  f_max(b.maxMii, b2.maxMii);
  f_max(b.adjMmax, b2.adjMmax);
  // The basis is searched on the first configuration only.
  b.has_std_basis = b1.has_std_basis;
  return b;
}

/*
  Whether the whole search provably fits int64_t: every runtime
  expression class is bounded against 2^61, with room to spare below the
  word, so the engine then runs on raw machine integers with no per
  operation check (the point of the analysis: SafeInt64 pays for every
  single operation, int64_t for none).

  Gb bounds the entries of every automorphism or isometry X appearing:
  row i of X has M-norm M_ii, and a vector v with v M v^T = c has
  |v_j| <= sqrt(c * (M^-1)_jj) <= sqrt(maxMii * adjMmax / detM0).
  The expressions:
  * E1 = n * Vmax * (n * Vmax * Mmax): the scalar product scans against
    the W tables, the W entries themselves, norms, Gram targets and the
    vector-sum tuples.
  * E2 = n * Vmax * Gb: the point images of ps_operate.
  * E3 = n * Gb^2: products of two automorphisms.
  * E5 = n * Mmax * Gb and E6 = n * adjMmax * E5: the two steps of
    ps_isometry_inverse; E6 also dominates detM0 and the adjugate
    entries.
  * E7 = m * Vmax and E8 = n^2 * Mmax * E7^2: the vector sums of one
    tuple class and their norms.
  * E9 = n * BadjMax * Vmax: the leaf products Badj * A, when the basis
    is not the standard one and comes as precomputed PSBasisData.
 */
template <typename Tint>
bool ps_certify_int64(PSInt64Bounds<Tint> const &b,
                      [[maybe_unused]] std::ostream &os) {
  if (!b.has_std_basis && !b.has_basis_data) {
#ifdef DEBUG_PLESKEN_SOUVIGNIER
    os << "PS: no int64, the family does not contain the standard basis\n";
#endif
    return false;
  }
  Tint CAP = UniversalScalarConversion<Tint, int64_t>(int64_t(1) << 61);
  Tint num = b.maxMii * b.adjMmax;
  Tint Gb(1);
  while (Gb * Gb * b.detM0 < num) {
    Gb *= 2;
    if (Gb > CAP) {
      return false;
    }
  }
  Gb += 1;
  Tint Wmax = b.n * b.Vmax * b.Mmax;
  Tint E1 = b.n * b.Vmax * Wmax;
  Tint E2 = b.n * b.Vmax * Gb;
  Tint E3 = b.n * Gb * Gb;
  Tint E5 = b.n * b.Mmax * Gb;
  Tint E6 = b.n * b.adjMmax * E5;
  Tint E7 = b.m * b.Vmax;
  Tint E8 = b.n * b.n * b.Mmax * E7 * E7;
  bool fits = E1 <= CAP && E2 <= CAP && E3 <= CAP && E5 <= CAP &&
              E6 <= CAP && E7 <= CAP && E8 <= CAP;
  if (!b.has_std_basis) {
    // The leaves of a basis drawn from the family: Badj * A with A the
    // stacked images, and the quotient by Bden.
    Tint E9 = b.n * b.BadjMax * b.Vmax;
    fits = fits && E9 <= CAP && b.BadjMax <= CAP;
  }
#ifdef DEBUG_PLESKEN_SOUVIGNIER
  os << "PS: int64 certificate fits=" << fits << "\n";
#endif
  return fits;
}

// The backtrack node budget above which the adaptive depth is raised.
// A depth-appropriate search visits a few hundred nodes at most (the
// largest healthy count observed on the d351 genus was ~160); a
// depth-starved one visits 10^4 to 10^7. A tight cap just above the
// healthy regime both separates the two and minimises the work wasted on
// the doomed low-depth attempt before escalating -- which matters on a
// large vector family, where every node scans the whole family. Only
// consulted when the depth is chosen automatically (depth == -1).
static const size_t PS_NODE_BUDGET = 500;

// The highest vector-sum depth the adaptive escalation will reach. Beyond
// it the tuple maps grow costly for little pruning gain, and a search that
// still overruns the budget is not depth-starved but genuinely large, so
// the last attempt runs to completion unbudgeted rather than escalating
// into ever-larger vector-sum data.
static const int PS_MAX_ADAPTIVE_DEPTH = 8;

template <typename Tint>
PleskenSouvignierAutomResult<Tint> PleskenSouvignierAutomorphism(
    std::vector<MyMatrix<Tint>> const &ListMat, MyMatrix<Tint> const &SHVhalf,
    std::ostream &os, int depth = -1,
    std::optional<std::pair<MyMatrix<Tint>, Tint>> const &adjdet0_opt = {},
    std::optional<PSBasisData<Tint>> const &basis_data = {}) {
  std::pair<MyMatrix<Tint>, Tint> adjdet0 = [&]() {
    if (adjdet0_opt) {
      return *adjdet0_opt;
    }
    if constexpr (std::is_same_v<Tint, int64_t>) {
      std::cerr << "PS: an int64 automorphism run needs the precomputed "
                << "adjugate of ListMat[0]\n";
      throw TerminalException{1};
    } else {
      return AdjugateDeterminant(ListMat[0]);
    }
  }();
  if constexpr (!std::is_same_v<Tint, int64_t>) {
    // The certified int64 fast path: the preanalysis proves that no
    // operation of the search can overflow, so the engine runs on raw
    // machine integers with no per-operation check.
    PSInt64Bounds<Tint> bounds = ps_int64_bounds(ListMat, SHVhalf, adjdet0);
    std::optional<PSBasisData<Tint>> bd =
        ps_basis_data_for_int64(ListMat, SHVhalf, adjdet0, bounds, os);
    if (ps_certify_int64(bounds, os)) {
      std::vector<MyMatrix<int64_t>> ListMat64;
      ListMat64.reserve(ListMat.size());
      for (auto &eMat : ListMat) {
        ListMat64.push_back(ps_matrix_to_int64(eMat));
      }
      MyMatrix<int64_t> SHVhalf64 = ps_matrix_to_int64(SHVhalf);
      std::optional<std::pair<MyMatrix<int64_t>, int64_t>> adjdet64 =
          std::make_pair(ps_matrix_to_int64(adjdet0.first),
                         UniversalScalarConversion<int64_t, Tint>(
                             adjdet0.second));
      std::optional<PSBasisData<int64_t>> bd64;
      if (bd) {
        bd64 = ps_basis_data_to_int64(*bd);
      }
      PleskenSouvignierAutomResult<int64_t> res64 =
          PleskenSouvignierAutomorphism<int64_t>(ListMat64, SHVhalf64, os,
                                                 depth, adjdet64, bd64);
      PleskenSouvignierAutomResult<Tint> res;
      res.ListGen.reserve(res64.ListGen.size());
      for (auto &eGen : res64.ListGen) {
        res.ListGen.push_back(ps_matrix_from_int64<Tint>(eGen));
      }
      res.ListOrbitSize = std::move(res64.ListOrbitSize);
      res.order = res64.order;
      res.depth_used = res64.depth_used;
      res.nodes_used = res64.nodes_used;
      return res;
    }
  }
  bool adaptive = (depth == -1);
  PleskenSouvignierContext<Tint> ctx = PleskenSouvignierBuildContext(
      ListMat, SHVhalf, true, depth, os,
      std::optional<std::pair<MyMatrix<Tint>, Tint>>(adjdet0), basis_data);
  int d = ctx.depth;
  while (true) {
    // Once the cap is reached the search runs unbudgeted, so it always
    // terminates rather than escalating without bound.
    bool budgeted = adaptive && d < PS_MAX_ADAPTIVE_DEPTH;
    ctx.node_count = 0;
    ctx.total_nodes = 0;
    ctx.node_budget = budgeted ? PS_NODE_BUDGET : 0;
    try {
      ps_auto(ctx, os);
      break;
    } catch (PSBudgetExceeded const &) {
      // The group found so far and the completed steps are kept; the
      // search resumes at ctx.resume_step (see the context).
      d++;
      ps_init_vector_sums(ctx, d, os);
#ifdef DEBUG_PLESKEN_SOUVIGNIER
      os << "PS: auto node budget exceeded, raising depth to " << d
         << ", resuming at step " << ctx.resume_step << "\n";
#endif
    }
  }
  PleskenSouvignierAutomResult<Tint> result;
  for (int i = 0; i < ctx.n(); i++) {
    for (int j = ctx.nsg[i]; j < static_cast<int>(ctx.g[i].size()); j++) {
      result.ListGen.push_back(ctx.g[i][j]);
    }
  }
  result.ListOrbitSize = ctx.orders;
  result.order = PleskenSouvignierGroupOrder(result.ListOrbitSize);
  result.depth_used = ctx.depth;
  result.nodes_used = ctx.total_nodes;
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  for (auto &eGen : result.ListGen) {
    for (auto &eMat : ListMat) {
      if (eGen * eMat * eGen.transpose() != eMat) {
        std::cerr << "PS: a generator does not preserve the configuration\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  return result;
}

/*
  A bounded harvest of elements of <ListGen> stabilizing the point pt,
  from the coincidences of the orbit transversal. Used only for pruning
  the isometry search, where any subgroup of the true stabilizer gives
  valid orbit deletions, so the bounds trade pruning power for time and
  nothing else.
 */
template <typename Tint>
std::vector<MyMatrix<Tint>>
ps_isostab(PleskenSouvignierContext<Tint> const &ctx, int pt,
           std::vector<MyMatrix<Tint>> const &ListGen) {
  int m = ctx.VS.m;
  int n = ctx.n();
  const size_t max_gens = 16;
  const int max_fail = 32;
  std::vector<MyMatrix<Tint>> H;
  if (ListGen.empty()) {
    return H;
  }
  std::vector<int> orb{pt};
  std::vector<uint8_t> flag(2 * m + 1, 0);
  std::vector<MyMatrix<Tint>> w(1);
  std::vector<int> wpos(2 * m + 1, -1);
  flag[pt + m] = 1;
  wpos[pt + m] = 0;
  w[0] = IdentityMat<Tint>(n);
  size_t cnd = 0;
  int fail = 0;
  while (cnd < orb.size() && fail < max_fail && H.size() < max_gens) {
    for (auto &eGen : ListGen) {
      int im = ps_operate(ctx, orb[cnd], eGen);
      if (flag[im + m] == 0) {
        flag[im + m] = 1;
        wpos[im + m] = w.size();
        w.push_back(w[wpos[orb[cnd] + m]] * eGen);
        orb.push_back(im);
      } else {
        MyMatrix<Tint> B = w[wpos[orb[cnd] + m]] * eGen;
        MyMatrix<Tint> const &Bknown = w[wpos[im + m]];
        if (B != Bknown) {
          MyMatrix<Tint> S = B * ps_aut_inverse(ctx, Bknown);
          bool is_known = false;
          for (auto &eH : H) {
            if (eH == S) {
              is_known = true;
              break;
            }
          }
          if (is_known) {
            fail++;
          } else {
            H.push_back(S);
            if (H.size() >= max_gens) {
              break;
            }
          }
        } else {
          fail++;
        }
      }
    }
    cnd++;
  }
  return H;
}

/*
  The isometry search: candidates for the images of the basis of the
  first lattice among the vectors of the second, with the orbit of a
  failed candidate under the stabilizer of the earlier choices deleted
  wholesale. ListGenAut2 is any set of automorphisms of the second
  configuration; the group they generate powers the orbit deletions, so
  a caller that knows Aut of the second lattice should pass it, and an
  empty list still leaves -Id.
 */
template <typename Tint>
bool ps_iso_extend(PleskenSouvignierContext<Tint> const &Ci,
                   PleskenSouvignierContext<Tint> const &Co, int step,
                   std::vector<int> &x, std::vector<int> &cand_step,
                   std::vector<MyMatrix<Tint>> const &H) {
  int n = Ci.n();
  int m = Co.VS.m;
  while (!cand_step.empty()) {
    int im = cand_step[0];
    // The budget is per top-level subtree (see ps_auto): reset at the root
    // of each one so a large target group is not read as an explosion.
    if (step == 0) {
      Ci.node_count = 0;
    }
    Ci.node_count++;
    Ci.total_nodes++;
    if (Ci.node_budget != 0 && Ci.node_count > Ci.node_budget) {
      throw PSBudgetExceeded{};
    }
    x[step] = im;
    if (step == n - 1) {
      if (Ci.Bden == 1 || ps_matgen_opt(Ci, Co, x).has_value()) {
        return true;
      }
    } else {
      std::vector<int> cand_next;
      if (ps_cand(Ci, Co, step + 1, x, cand_next)) {
        std::vector<MyMatrix<Tint>> Hnext = ps_isostab(Co, im, H);
        if (ps_iso_extend(Ci, Co, step + 1, x, cand_next, Hnext)) {
          return true;
        }
      }
    }
    std::vector<int> oc = ps_orbit(Co, {im}, H);
    if (oc.size() == 1) {
      // ps_setdiff would only remove im itself; keep the common case
      // cheap.
      cand_step.erase(cand_step.begin());
    } else {
      ps_setdiff(cand_step, oc, m);
    }
  }
  return false;
}

// The isometry, as P with P * ListMat1[i] * P^T = ListMat2[i] for all i,
// or nothing. The two families have to be built by the same invariant
// rule: vectors of norm at most B for the SAME bound B on both sides
// (see the src_latt wrappers).
template <typename Tint>
std::optional<MyMatrix<Tint>> PleskenSouvignierIsometry(
    std::vector<MyMatrix<Tint>> const &ListMat1,
    MyMatrix<Tint> const &SHVhalf1,
    std::vector<MyMatrix<Tint>> const &ListMat2,
    MyMatrix<Tint> const &SHVhalf2,
    std::vector<MyMatrix<Tint>> const &ListGenAut2, std::ostream &os,
    int depth = -1,
    std::optional<std::pair<MyMatrix<Tint>, Tint>> const &adjdet1_opt = {},
    std::optional<std::pair<MyMatrix<Tint>, Tint>> const &adjdet2_opt = {},
    std::optional<PSBasisData<Tint>> const &basis_data1 = {}) {
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  MicrosecondTime time;
#endif
  if (ListMat1.size() != ListMat2.size()) {
    return {};
  }
  if (SHVhalf1.rows() != SHVhalf2.rows()) {
    return {};
  }
  // An isometry P has det P = +-1, so P M1 P^T = M2 forces equal
  // determinants. Without this, families of equal size on lattices of
  // different determinants could match through an integral P of |det| > 1,
  // whose inverse is not integral. The int64 recursion receives the
  // adjugates from the caller, who has already run the check; the
  // Bareiss minors of DeterminantMat do not fit a machine word anyway.
  if (!adjdet1_opt) {
    if constexpr (std::is_same_v<Tint, int64_t>) {
      std::cerr << "PS: an int64 isometry run needs the precomputed "
                << "adjugates\n";
      throw TerminalException{1};
    } else {
      for (size_t iMat = 0; iMat < ListMat1.size(); iMat++) {
        if (DeterminantMat(ListMat1[iMat]) != DeterminantMat(ListMat2[iMat])) {
          return {};
        }
      }
    }
  }
  std::pair<MyMatrix<Tint>, Tint> adjdet1 = [&]() {
    if (adjdet1_opt) {
      return *adjdet1_opt;
    }
    if constexpr (std::is_same_v<Tint, int64_t>) {
      throw TerminalException{1};
    } else {
      return AdjugateDeterminant(ListMat1[0]);
    }
  }();
  std::pair<MyMatrix<Tint>, Tint> adjdet2 = [&]() {
    if (adjdet2_opt) {
      return *adjdet2_opt;
    }
    if constexpr (std::is_same_v<Tint, int64_t>) {
      throw TerminalException{1};
    } else {
      return AdjugateDeterminant(ListMat2[0]);
    }
  }();
  if constexpr (!std::is_same_v<Tint, int64_t>) {
    // The certified int64 fast path over the maxima of both sides.
    PSInt64Bounds<Tint> bounds = ps_int64_bounds_merge(
        ps_int64_bounds(ListMat1, SHVhalf1, adjdet1),
        ps_int64_bounds(ListMat2, SHVhalf2, adjdet2));
    // The basis is searched on the first configuration only.
    std::optional<PSBasisData<Tint>> bd =
        ps_basis_data_for_int64(ListMat1, SHVhalf1, adjdet1, bounds, os);
    if (ps_certify_int64(bounds, os)) {
      std::vector<MyMatrix<int64_t>> ListMat1_64;
      std::vector<MyMatrix<int64_t>> ListMat2_64;
      std::vector<MyMatrix<int64_t>> ListGenAut2_64;
      for (auto &eMat : ListMat1) {
        ListMat1_64.push_back(ps_matrix_to_int64(eMat));
      }
      for (auto &eMat : ListMat2) {
        ListMat2_64.push_back(ps_matrix_to_int64(eMat));
      }
      for (auto &eGen : ListGenAut2) {
        ListGenAut2_64.push_back(ps_matrix_to_int64(eGen));
      }
      MyMatrix<int64_t> SHVhalf1_64 = ps_matrix_to_int64(SHVhalf1);
      MyMatrix<int64_t> SHVhalf2_64 = ps_matrix_to_int64(SHVhalf2);
      auto f_adj = [](std::pair<MyMatrix<Tint>, Tint> const &pr)
          -> std::optional<std::pair<MyMatrix<int64_t>, int64_t>> {
        return std::make_pair(
            ps_matrix_to_int64(pr.first),
            UniversalScalarConversion<int64_t, Tint>(pr.second));
      };
      std::optional<MyMatrix<int64_t>> opt64 =
          PleskenSouvignierIsometry<int64_t>(
              ListMat1_64, SHVhalf1_64, ListMat2_64, SHVhalf2_64,
              ListGenAut2_64, os, depth, f_adj(adjdet1), f_adj(adjdet2),
              bd ? std::optional<PSBasisData<int64_t>>(
                       ps_basis_data_to_int64(*bd))
                 : std::optional<PSBasisData<int64_t>>());
      if (!opt64) {
        return {};
      }
      return ps_matrix_from_int64<Tint>(*opt64);
    }
  }
  bool adaptive = (depth == -1);
  PleskenSouvignierContext<Tint> Ci = PleskenSouvignierBuildContext(
      ListMat1, SHVhalf1, true, depth, os,
      std::optional<std::pair<MyMatrix<Tint>, Tint>>(adjdet1), basis_data1);
  PleskenSouvignierContext<Tint> Co = PleskenSouvignierBuildContext(
      ListMat2, SHVhalf2, false, 0, os,
      std::optional<std::pair<MyMatrix<Tint>, Tint>>(adjdet2));
  int n = Ci.n();
  std::vector<int> x(n, 0);
  std::vector<int> cand0_base;
  // The level-0 candidates use only the fingerprint, not the vector sums,
  // so they are independent of the depth and computed once.
  if (!ps_cand(Ci, Co, 0, x, cand0_base)) {
    return {};
  }
  std::vector<MyMatrix<Tint>> H = ListGenAut2;
  H.push_back(-IdentityMat<Tint>(n));
  int d = Ci.depth;
  bool found = false;
  // Kept across depth escalations: ps_iso_extend removes a level-0
  // candidate (with its orbit) only once its subtree has been searched to
  // the end, and that refutation does not depend on the depth.
  std::vector<int> cand0 = cand0_base;
  while (true) {
    bool budgeted = adaptive && d < PS_MAX_ADAPTIVE_DEPTH;
    Ci.node_count = 0;
    Ci.total_nodes = 0;
    Ci.node_budget = budgeted ? PS_NODE_BUDGET : 0;
    std::fill(x.begin(), x.end(), 0);
    try {
      found = ps_iso_extend(Ci, Co, 0, x, cand0, H);
      break;
    } catch (PSBudgetExceeded const &) {
      d++;
      ps_init_vector_sums(Ci, d, os);
#ifdef DEBUG_PLESKEN_SOUVIGNIER
      os << "PS: isometry node budget exceeded, raising depth to " << d << "\n";
#endif
    }
  }
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  os << "|PS: isometry found=" << found << "|=" << time << "\n";
#endif
  if (!found) {
    return {};
  }
  // The assignment gives Q with Q * M2 * Q^T = M1; the repository
  // convention is the inverse direction, computed through the forms so
  // that the intermediates stay within the certified bounds.
  MyMatrix<Tint> Q = ps_matgen(Ci, Co, x);
  MyMatrix<Tint> P =
      ps_isometry_inverse(ListMat2[0], Ci.adjM0, Ci.detM0, Q);
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  for (size_t iMat = 0; iMat < ListMat1.size(); iMat++) {
    if (P * ListMat1[iMat] * P.transpose() != ListMat2[iMat]) {
      std::cerr << "PS: the isometry does not map ListMat1 to ListMat2\n";
      throw TerminalException{1};
    }
  }
#endif
  return P;
}

// clang-format off
#endif  // SRC_GROUP_PLESKENSOUVIGNIER_H_
// clang-format on

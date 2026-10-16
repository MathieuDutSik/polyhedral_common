// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLYNORM_POLYNORM_PACKINGCONTACTS_H_
#define SRC_POLYNORM_POLYNORM_PACKINGCONTACTS_H_

// clang-format off
#include "POLY_DirectDualDesc.h"
#include "POLY_LinearProgramming.h"
#include "PolytopeEquiStab.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_POLYNORM_PACKING_CONTACTS
#endif

/*
  The densest lattice packing of a polytope P of dimension n by the
  enumeration of the contact sets of a critical lattice of D = P - P.

  Let Lambda be a critical lattice of D (admissible, of minimal determinant),
  S the set of its points on the boundary of D, and b_1, ..., b_n in S
  independent with |det(b_1, ..., b_n)| maximal. Then:
  (1) conv(+-b_i) is a lattice cross-polytope (it lies in D, whose interior
      has no lattice point but 0), so in the coordinates of the b_i the
      lattice is one of the finitely many classes L >= Z^n of lattice
      cross-polytopes;
  (2) every s in S has coordinates in [-1, 1]^n (Cramer: replacing b_i by s
      does not increase |det|), and every n elements of S have |det| <= 1
      in these coordinates;
  (3) S has at least n(n+1)/2 pairs (Swinnerton-Dyer);
  (4) for every choice of a facet F(s) of D containing each s, the basis W
      (rows b_i) is a local minimum of |det| on {W : s W in H_F(s)}: moving
      in that space keeps every s on a supporting hyperplane, hence outside
      the interior of D, and the other lattice points stay outside by
      continuity, so the nearby lattices are admissible (the argument (2.6)
      of Betke and Henk). On a space of dimension d the second order
      condition needs d <= n(n-1)/2, that is a rank >= n(n+1)/2.
  So the critical determinant is the smallest |det W| / [L : Z^n] over the
  admissible local minima of |det| on the affine spaces of the
  configurations (L, S, F) satisfying (1)-(4). The configurations are
  enumerated by backtracking over the candidate contact vectors of each
  class, each excluded or put on a facet, with the exact pruning of the
  linear programs (the points in their facets), of the additive relations
  (c = a +- b forces (F_a +- F_b) cap F_c nonempty) and of the
  determinants; the affine spaces are collected without repetition.
 */

// The classes of lattice cross-polytopes of dimension n, by the enumeration
// along flags (see Paper_poly_norm/LatticeCrossPolytopes): the group L/Z^n
// as numerators modulo the common denominator E.
struct PCClass {
  int64_t E;
  std::vector<std::vector<int64_t>> elts;
};

inline std::vector<PCClass> PC_LatticeCrossPolytopes(int n) {
  using Vec = std::vector<int64_t>;
  struct Lat {
    int k;
    int64_t E;
    std::vector<Vec> elts;
  };
  auto fact = [](int k) -> int64_t {
    int64_t f = 1;
    for (int i = 2; i <= k; i++)
      f *= i;
    return f;
  };
  auto bound = [&](int k) -> int64_t {
    return k <= 2 ? fact(k) : fact(k) - 1;
  };
  auto signed_perms = [](int k) {
    std::vector<std::pair<std::vector<int>, std::vector<int>>> L;
    std::vector<int> p(k);
    for (int i = 0; i < k; i++)
      p[i] = i;
    do {
      for (int mask = 0; mask < (1 << k); mask++) {
        std::vector<int> s(k);
        for (int i = 0; i < k; i++)
          s[i] = (mask >> i) & 1 ? -1 : 1;
        L.push_back({p, s});
      }
    } while (std::next_permutation(p.begin(), p.end()));
    return L;
  };
  auto reduce = [](Lat &L) {
    int64_t g = L.E;
    for (auto &x : L.elts)
      for (int i = 0; i < L.k; i++)
        g = std::gcd(g, x[i]);
    if (g > 1) {
      for (auto &x : L.elts)
        for (int i = 0; i < L.k; i++)
          x[i] /= g;
      L.E /= g;
    }
  };
  auto canonical = [](Lat const &L, auto const &SP) {
    std::vector<int64_t> best;
    for (auto const &s : SP) {
      std::vector<Vec> img;
      for (auto const &x : L.elts) {
        Vec v(L.k);
        for (int i = 0; i < L.k; i++) {
          int64_t c = s.second[i] * x[s.first[i]] % L.E;
          if (c < 0)
            c += L.E;
          v[i] = c;
        }
        img.push_back(v);
      }
      std::sort(img.begin(), img.end());
      std::vector<int64_t> flat{L.E};
      for (auto &v : img)
        flat.insert(flat.end(), v.begin(), v.end());
      if (best.empty() || flat < best)
        best = flat;
    }
    return best;
  };
  auto basis_EL = [](Lat const &L) {
    int k = L.k;
    std::vector<Vec> rows;
    for (int i = 0; i < k; i++) {
      Vec v(k, 0);
      v[i] = L.E;
      rows.push_back(v);
    }
    for (auto &x : L.elts)
      rows.push_back(x);
    std::vector<Vec> H;
    for (int c = 0; c < k; c++) {
      while (true) {
        int piv = -1;
        for (size_t r = 0; r < rows.size(); r++)
          if (rows[r][c] != 0 &&
              (piv < 0 || std::llabs(rows[r][c]) < std::llabs(rows[piv][c])))
            piv = r;
        if (piv < 0)
          break;
        bool done = true;
        for (size_t r = 0; r < rows.size(); r++) {
          if (static_cast<int>(r) == piv || rows[r][c] == 0)
            continue;
          int64_t q = rows[r][c] / rows[piv][c];
          for (int i = 0; i < k; i++)
            rows[r][i] -= q * rows[piv][i];
          if (rows[r][c] != 0)
            done = false;
        }
        if (done) {
          H.push_back(rows[piv]);
          rows.erase(rows.begin() + piv);
          break;
        }
      }
    }
    return H;
  };
  std::vector<Lat> level{Lat{1, 1, {Vec(1, 0)}}};
  for (int k = 1; k < n; k++) {
    int64_t B = bound(k + 1);
    auto SP = signed_perms(k + 1);
    std::map<std::vector<int64_t>, Lat> found;
    for (auto const &L : level) {
      int64_t ord = L.elts.size();
      std::vector<Vec> H = basis_EL(L);
      bool is_Zk = ord == 1;
      for (int64_t q = 1; q * ord <= B; q++) {
        std::vector<int64_t> c(k, 0);
        int64_t cmax = is_Zk ? q / 2 : q - 1;
        while (true) {
          bool ok_sorted = true;
          if (is_Zk)
            for (int i = 0; i + 1 < k; i++)
              if (c[i] < c[i + 1])
                ok_sorted = false;
          if (ok_sorted) {
            int64_t E = L.E, D = q * E;
            Vec v(k + 1, 0);
            for (int j = 0; j < k; j++)
              for (int i = 0; i < k; i++)
                v[i] += c[j] * H[j][i];
            v[k] = E;
            for (int i = 0; i <= k; i++) {
              v[i] %= D;
              if (v[i] < 0)
                v[i] += D;
            }
            bool adm = true;
            std::vector<Vec> elts;
            for (int64_t j = 0; j < q && adm; j++) {
              for (auto const &g : L.elts) {
                Vec x(k + 1);
                int64_t nrm = 0;
                bool zero = true;
                for (int i = 0; i <= k; i++) {
                  int64_t a = ((i < k ? q * g[i] : 0) + j * v[i]) % D;
                  x[i] = a;
                  if (a != 0)
                    zero = false;
                  nrm += std::min(a, D - a);
                }
                if (!zero && nrm < D) {
                  adm = false;
                  break;
                }
                elts.push_back(x);
              }
            }
            if (adm) {
              Lat L2{k + 1, D, std::move(elts)};
              reduce(L2);
              auto key = canonical(L2, SP);
              if (!found.count(key))
                found.emplace(std::move(key), std::move(L2));
            }
          }
          int i = 0;
          while (i < k && c[i] == cmax)
            c[i++] = 0;
          if (i == k || q == 1)
            break;
          c[i]++;
        }
      }
    }
    level.clear();
    for (auto &kv : found)
      level.push_back(kv.second);
  }
  std::vector<PCClass> ListClass;
  for (auto &L : level)
    ListClass.push_back({L.E, L.elts});
  return ListClass;
}

// The difference body D = P - P.
template <typename T> struct PCBody {
  int n;
  int m;
  MyMatrix<T> Lmat;
  MyMatrix<double> Lmat_d;
  std::vector<int> opp;
  std::vector<MyVector<T>> ListPoint;
  std::vector<std::vector<int>> neighbors;
  std::vector<std::vector<double>> box_lo, box_hi;
  std::vector<int> orbit;
  std::vector<uint8_t> is_rep;
  size_t n_orbit;
  T vol_P;
  double r_max;
};

template <typename T, typename Tgroup>
PCBody<T> PC_BuildBody(MyMatrix<T> const &EXT_P, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  int n = EXT_P.cols() - 1;
  PCBody<T> body;
  body.n = n;
  body.vol_P = rev_search::Kernel_VolumePolytope(EXT_P);
  int n_vert = EXT_P.rows();
  std::set<std::vector<T>> set_diff;
  for (int i = 0; i < n_vert; i++)
    for (int j = 0; j < n_vert; j++)
      if (i != j) {
        std::vector<T> v(n);
        for (int k = 0; k < n; k++)
          v[k] = EXT_P(i, k + 1) - EXT_P(j, k + 1);
        set_diff.insert(v);
      }
  int n_diff = set_diff.size();
  MyMatrix<T> EXTdiff(n_diff, n + 1);
  int pos = 0;
  for (auto &v : set_diff) {
    EXTdiff(pos, 0) = 1;
    for (int k = 0; k < n; k++)
      EXTdiff(pos, k + 1) = v[k];
    MyVector<T> V(n);
    for (int k = 0; k < n; k++)
      V(k) = v[k];
    body.ListPoint.push_back(V);
    pos++;
  }
  vectface vf = DirectFacetComputationIncidence(EXTdiff, "cdd", os);
  int m = vf.size();
  body.m = m;
  body.Lmat = MyMatrix<T>(m, n);
  for (int j = 0; j < m; j++) {
    MyVector<T> w = FindFacetInequality(EXTdiff, vf[j]);
    for (int k = 0; k < n; k++)
      body.Lmat(j, k) = -w(k + 1) / w(0);
  }
  body.Lmat_d = UniversalMatrixConversion<double, T>(body.Lmat);
  body.opp.assign(m, -1);
  for (int j = 0; j < m; j++)
    for (int k = 0; k < m; k++) {
      bool same = true;
      for (int u = 0; u < n; u++)
        if (body.Lmat(j, u) != -body.Lmat(k, u))
          same = false;
      if (same)
        body.opp[j] = k;
    }
  // Points on each facet, bounding boxes, ridge neighbors.
  std::vector<std::vector<int>> fpts(m);
  body.box_lo.assign(m, std::vector<double>(n, 1e300));
  body.box_hi.assign(m, std::vector<double>(n, -1e300));
  body.r_max = 0;
  for (int i = 0; i < n_diff; i++) {
    double nrm = 0;
    for (int k = 0; k < n; k++) {
      double x = UniversalScalarConversion<double, T>(body.ListPoint[i](k));
      nrm += x * x;
    }
    body.r_max = std::max(body.r_max, std::sqrt(nrm));
    for (int j = 0; j < m; j++) {
      T s(0);
      for (int k = 0; k < n; k++)
        s += body.Lmat(j, k) * body.ListPoint[i](k);
      if (s == 1) {
        fpts[j].push_back(i);
        for (int k = 0; k < n; k++) {
          double x =
              UniversalScalarConversion<double, T>(body.ListPoint[i](k));
          body.box_lo[j][k] = std::min(body.box_lo[j][k], x);
          body.box_hi[j][k] = std::max(body.box_hi[j][k], x);
        }
      }
    }
  }
  body.neighbors.resize(m);
  for (int j = 0; j < m; j++)
    for (int k = 0; k < m; k++) {
      if (k == j)
        continue;
      std::vector<int> common;
      for (auto &i : fpts[j])
        if (std::find(fpts[k].begin(), fpts[k].end(), i) != fpts[k].end())
          common.push_back(i);
      if (static_cast<int>(common.size()) >= n - 1) {
        MyMatrix<T> M(common.size(), n + 1);
        for (size_t u = 0; u < common.size(); u++) {
          M(u, 0) = 1;
          for (int c = 0; c < n; c++)
            M(u, c + 1) = body.ListPoint[common[u]](c);
        }
        if (RankMat(M) == n - 1)
          body.neighbors[j].push_back(k);
      }
    }
  // Orbits of facets under the linear symmetries of D.
  std::vector<Telt> LGen = LinPolytope_Automorphism_LGen<T, Tgroup>(body.Lmat, os);
  body.orbit.assign(m, -1);
  body.is_rep.assign(m, 0);
  int n_orbit = 0;
  for (int j = 0; j < m; j++) {
    if (body.orbit[j] != -1)
      continue;
    body.is_rep[j] = 1;
    std::vector<int> stack{j};
    body.orbit[j] = n_orbit;
    while (!stack.empty()) {
      int u = stack.back();
      stack.pop_back();
      for (auto &g : LGen) {
        int v = g.at(u);
        if (body.orbit[v] == -1) {
          body.orbit[v] = n_orbit;
          stack.push_back(v);
        }
      }
    }
    n_orbit++;
  }
  body.n_orbit = n_orbit;
  return body;
}

// A configuration: an affine space {W : sum_k A_k W = b_k} in the n^2
// entries of W, for a class.
struct PCLeaf {
  int i_class;
  int n_contact;
  int rank;
  // The exact reduced row echelon form, as a key.
  std::vector<std::string> key;
};

template <typename T> struct PCSearch {
  PCBody<T> const &body;
  std::vector<PCClass> const &ListClass;
  std::ostream &os;
  int n;
  // The leaves, without repetition.
  std::map<std::pair<int, std::vector<std::string>>, size_t> map_leaf;
  std::vector<PCLeaf> ListLeaf;
  std::vector<std::vector<std::pair<MyVector<T>, T>>> ListLeafEq;
  std::mutex mtx;
  // Statistics.
  std::atomic<size_t> n_node{0}, n_lp{0}, n_leaf_raw{0}, n_triple_lp{0};
  // Memo of the triple tests (F_a + F_b) cap F_c.
  std::unordered_map<uint64_t, bool> memo_triple;
  std::mutex mtx_triple;

  PCSearch(PCBody<T> const &_body, std::vector<PCClass> const &_ListClass,
           std::ostream &_os)
      : body(_body), ListClass(_ListClass), os(_os), n(_body.n) {}

  // c W d_f = level as a row over the n^2 entries (c has denominator E).
  MyVector<T> RowOf(std::vector<int64_t> const &c, int64_t E,
                    MyVector<T> const &d) const {
    MyVector<T> R(n * n);
    for (int a = 0; a < n; a++)
      for (int b = 0; b < n; b++)
        R(a * n + b) = T(static_cast<long>(c[a])) / T(static_cast<long>(E)) * d(b);
    return R;
  }

  bool Feasible(std::vector<std::vector<int64_t>> const &pts,
                std::vector<int> const &facets, int64_t E) {
    n_lp++;
    std::vector<MyVector<T>> rows;
    auto push = [&](MyVector<T> const &R, T const &c0) {
      MyVector<T> X(n * n + 1);
      X(0) = c0;
      for (int u = 0; u < n * n; u++)
        X(u + 1) = R(u);
      rows.push_back(X);
    };
    for (size_t i = 0; i < pts.size(); i++) {
      MyVector<T> d = GetMatrixRow(body.Lmat, facets[i]);
      MyVector<T> R = RowOf(pts[i], E, d);
      push(R, T(-1));
      push(-R, T(1));
      for (auto &g : body.neighbors[facets[i]]) {
        MyVector<T> dg = GetMatrixRow(body.Lmat, g);
        push(-RowOf(pts[i], E, dg), T(1));
      }
    }
    MyMatrix<T> M(rows.size(), n * n + 1);
    for (size_t i = 0; i < rows.size(); i++)
      for (int u = 0; u <= n * n; u++)
        M(i, u) = rows[i](u);
    MyVector<T> obj = ZeroVector<T>(n * n + 1);
    LpSolution<T> sol = SIMPLEX_LinearProgramming(M, obj, os);
    return sol.DirectSolution.has_value() && sol.DualSolution.has_value();
  }

  // (F_a + F_b) cap F_c nonempty, by a box test and an exact program.
  bool Triple(int a, int b, int c) {
    uint64_t key = (uint64_t(a) * 1024 + b) * 1024 + c;
    {
      std::lock_guard<std::mutex> lock(mtx_triple);
      auto it = memo_triple.find(key);
      if (it != memo_triple.end())
        return it->second;
    }
    bool res = true;
    for (int k = 0; k < n; k++) {
      double lo = body.box_lo[a][k] + body.box_lo[b][k];
      double hi = body.box_hi[a][k] + body.box_hi[b][k];
      if (hi < body.box_lo[c][k] - 1e-9 || body.box_hi[c][k] < lo - 1e-9)
        res = false;
    }
    if (res) {
      n_triple_lp++;
      // Variables x, y in R^n: x in F_a, y in F_b, x + y in F_c.
      std::vector<std::vector<T>> rows;
      auto push_facet = [&](int f, int which) {
        // which: 0 x, 1 y, 2 x + y
        auto mk = [&](int g, int sign, T c0) {
          std::vector<T> R(2 * n + 1, T(0));
          R[0] = c0;
          for (int k = 0; k < n; k++) {
            if (which != 1)
              R[1 + k] = T(sign) * body.Lmat(g, k);
            if (which != 0)
              R[1 + n + k] = T(sign) * body.Lmat(g, k);
          }
          rows.push_back(R);
        };
        mk(f, 1, T(-1));
        mk(f, -1, T(1));
        for (auto &g : body.neighbors[f])
          mk(g, -1, T(1));
      };
      push_facet(a, 0);
      push_facet(b, 1);
      push_facet(c, 2);
      MyMatrix<T> M(rows.size(), 2 * n + 1);
      for (size_t i = 0; i < rows.size(); i++)
        for (int u = 0; u <= 2 * n; u++)
          M(i, u) = rows[i][u];
      MyVector<T> obj = ZeroVector<T>(2 * n + 1);
      LpSolution<T> sol = SIMPLEX_LinearProgramming(M, obj, os);
      res = sol.DirectSolution.has_value() && sol.DualSolution.has_value();
    }
    std::lock_guard<std::mutex> lock(mtx_triple);
    memo_triple[key] = res;
    return res;
  }

  static int64_t Det(std::vector<std::vector<int64_t>> const &M) {
    // Bareiss, exact for small entries.
    int k = M.size();
    std::vector<std::vector<__int128>> A(k, std::vector<__int128>(k));
    for (int i = 0; i < k; i++)
      for (int j = 0; j < k; j++)
        A[i][j] = M[i][j];
    __int128 prev = 1;
    int sign = 1;
    for (int c = 0; c < k - 1; c++) {
      if (A[c][c] == 0) {
        int sw = -1;
        for (int i = c + 1; i < k; i++)
          if (A[i][c] != 0) {
            sw = i;
            break;
          }
        if (sw < 0)
          return 0;
        std::swap(A[c], A[sw]);
        sign = -sign;
      }
      for (int i = c + 1; i < k; i++)
        for (int j = c + 1; j < k; j++)
          A[i][j] = (A[i][j] * A[c][c] - A[i][c] * A[c][j]) / prev;
      prev = A[c][c];
    }
    return static_cast<int64_t>(sign * A[k - 1][k - 1]);
  }

  // The search for one class, restricted to the first facet f0 (one
  // representative per orbit under the symmetries of D).
  void SearchClass(int i_class, int f0) {
    PCClass const &C = ListClass[i_class];
    int64_t E = C.E;
    // Candidates: L cap [-1,1]^n up to sign; the unit vectors first.
    std::set<std::vector<int64_t>> set_cand;
    for (auto const &g : C.elts) {
      std::vector<int> z(n, -2);
      while (true) {
        std::vector<int64_t> v(n);
        bool ok = true, zero = true;
        for (int i = 0; i < n; i++) {
          v[i] = g[i] + E * z[i];
          if (std::llabs(v[i]) > E)
            ok = false;
          if (v[i] != 0)
            zero = false;
        }
        if (ok && !zero) {
          int64_t f = 0;
          for (int i = 0; i < n; i++)
            if (v[i] != 0) {
              f = v[i];
              break;
            }
          if (f < 0)
            for (auto &x : v)
              x = -x;
          set_cand.insert(v);
        }
        int i = 0;
        while (i < n && z[i] == 1)
          z[i++] = -2;
        if (i == n)
          break;
        z[i]++;
      }
    }
    std::vector<std::vector<int64_t>> cand;
    for (int i = 0; i < n; i++) {
      std::vector<int64_t> e(n, 0);
      e[i] = E;
      cand.push_back(e);
      set_cand.erase(e);
    }
    for (auto &v : set_cand)
      cand.push_back(v);
    int nc = cand.size();
    int64_t det_bound = 1;
    for (int i = 0; i < n; i++)
      det_bound *= E;
    // The additive relations c_k = s_i c_i + s_j c_j with i, j < k.
    std::vector<std::vector<std::array<int, 4>>> rel(nc);
    for (int k = 0; k < nc; k++)
      for (int i = 0; i < k; i++)
        for (int j = i + 1; j < k; j++)
          for (int si : {1, -1})
            for (int sj : {1, -1}) {
              bool match = true;
              for (int u = 0; u < n; u++)
                if (cand[k][u] != si * cand[i][u] + sj * cand[j][u])
                  match = false;
              if (match)
                rel[k].push_back({i, j, si, sj});
            }
    int need = n * (n + 1) / 2;
    std::vector<int> in_facet(nc, -1);  // -1 excluded, else the facet
    std::vector<std::vector<int64_t>> pts;
    std::vector<int> facets;
    std::function<void(int)> rec = [&](int k) {
      n_node++;
      if (static_cast<int>(pts.size()) + (nc - k) < need)
        return;
      if (k == nc) {
        Leaf(i_class, pts, facets, E);
        return;
      }
      // Exclude (not for the unit vectors, which are contacts).
      if (k >= n) {
        in_facet[k] = -1;
        rec(k + 1);
      }
      // Include on a facet.
      for (int f = 0; f < body.m; f++) {
        if (k == 0 && f != f0)
          continue;
        // Determinant condition with every n-1 earlier points.
        bool ok = true;
        int np = pts.size();
        if (np >= n - 1) {
          std::vector<int> sub(n - 1);
          for (int u = 0; u < n - 1; u++)
            sub[u] = u;
          while (ok) {
            std::vector<std::vector<int64_t>> M;
            for (int u = 0; u < n - 1; u++)
              M.push_back(pts[sub[u]]);
            M.push_back(cand[k]);
            if (std::llabs(Det(M)) > det_bound)
              ok = false;
            int a = n - 2;
            while (a >= 0 && sub[a] == np - (n - 1) + a)
              a--;
            if (a < 0)
              break;
            sub[a]++;
            for (int b = a + 1; b < n - 1; b++)
              sub[b] = sub[b - 1] + 1;
          }
          if (!ok)
            break;  // independent of f
        }
        // Additive relations with included points.
        for (auto &r : rel[k]) {
          if (in_facet[r[0]] < 0 || in_facet[r[1]] < 0)
            continue;
          int fa = r[2] == 1 ? in_facet[r[0]] : body.opp[in_facet[r[0]]];
          int fb = r[3] == 1 ? in_facet[r[1]] : body.opp[in_facet[r[1]]];
          if (!Triple(fa, fb, f)) {
            ok = false;
            break;
          }
        }
        if (!ok)
          continue;
        pts.push_back(cand[k]);
        facets.push_back(f);
        in_facet[k] = f;
        if (Feasible(pts, facets, E))
          rec(k + 1);
        in_facet[k] = -1;
        pts.pop_back();
        facets.pop_back();
      }
    };
    rec(0);
  }

  // Stage 2: the facet assignments of a given contact set S (of the class
  // i_class, with the common denominator E), the facet of the first point
  // being f0 (one representative per orbit of facets of D).
  void SearchSet(int i_class, std::vector<std::vector<int64_t>> const &S,
                 int f0) {
    int64_t E = ListClass[i_class].E;
    int k = S.size();
    std::vector<std::vector<std::array<int, 4>>> rel(k);
    for (int c = 0; c < k; c++)
      for (int a = 0; a < c; a++)
        for (int b = a + 1; b < c; b++)
          for (int sa : {1, -1})
            for (int sb : {1, -1}) {
              bool match = true;
              for (int u = 0; u < n; u++)
                if (S[c][u] != sa * S[a][u] + sb * S[b][u])
                  match = false;
              if (match)
                rel[c].push_back({a, b, sa, sb});
            }
    std::vector<int> facets;
    std::vector<std::vector<int64_t>> pts;
    std::function<void(int)> rec = [&](int i) {
      n_node++;
      if (i == k) {
        Leaf(i_class, pts, facets, E);
        return;
      }
      for (int f = 0; f < body.m; f++) {
        if (i == 0 && f != f0)
          continue;
        bool ok = true;
        for (auto &r : rel[i]) {
          int fa = r[2] == 1 ? facets[r[0]] : body.opp[facets[r[0]]];
          int fb = r[3] == 1 ? facets[r[1]] : body.opp[facets[r[1]]];
          if (!Triple(fa, fb, f)) {
            ok = false;
            break;
          }
        }
        if (!ok)
          continue;
        pts.push_back(S[i]);
        facets.push_back(f);
        if (i < n || Feasible(pts, facets, E))
          rec(i + 1);
        pts.pop_back();
        facets.pop_back();
      }
    };
    rec(0);
  }

  void Leaf(int i_class, std::vector<std::vector<int64_t>> const &pts,
            std::vector<int> const &facets, int64_t E) {
    n_leaf_raw++;
    int k = pts.size();
    MyMatrix<T> A(k, n * n + 1);
    for (int i = 0; i < k; i++) {
      MyVector<T> R = RowOf(pts[i], E, GetMatrixRow(body.Lmat, facets[i]));
      for (int u = 0; u < n * n; u++)
        A(i, u) = R(u);
      A(i, n * n) = 1;
    }
    // Reduced row echelon form of the augmented matrix.
    MyMatrix<T> R = A;
    int rank = 0;
    for (int c = 0; c < n * n && rank < k; c++) {
      int piv = -1;
      for (int i = rank; i < k; i++)
        if (R(i, c) != 0) {
          piv = i;
          break;
        }
      if (piv < 0)
        continue;
      R.row(rank).swap(R.row(piv));
      T inv = 1 / R(rank, c);
      for (int u = 0; u <= n * n; u++)
        R(rank, u) *= inv;
      for (int i = 0; i < k; i++)
        if (i != rank && R(i, c) != 0) {
          T f = R(i, c);
          for (int u = 0; u <= n * n; u++)
            R(i, u) -= f * R(rank, u);
        }
      rank++;
    }
    // Inconsistent system: a row 0 = 1.
    for (int i = rank; i < k; i++)
      if (R(i, n * n) != 0)
        return;
    if (rank < n * (n + 1) / 2)
      return;
    std::vector<std::string> key;
    std::vector<std::pair<MyVector<T>, T>> eqs;
    for (int i = 0; i < rank; i++) {
      MyVector<T> row(n * n);
      for (int u = 0; u <= n * n; u++) {
        std::ostringstream s;
        s << R(i, u);
        key.push_back(s.str());
      }
      for (int u = 0; u < n * n; u++)
        row(u) = R(i, u);
      eqs.push_back({row, R(i, n * n)});
    }
    std::lock_guard<std::mutex> lock(mtx);
    auto kk = std::make_pair(i_class, key);
    if (map_leaf.count(kk))
      return;
    map_leaf[kk] = ListLeaf.size();
    ListLeaf.push_back({i_class, k, rank, key});
    ListLeafEq.push_back(eqs);
  }
};

/*
  Stage 3: the critical points of det on an affine space of dimension d
  (d <= n(n-1)/2), numerically. det is a polynomial of degree n in the
  parameters, its gradient d polynomials of degree n-1, solved by a total
  degree homotopy in projective coordinates (start system
  X_i^{n-1} - X_0^{n-1}, (n-1)^d paths) with the gamma trick and retries;
  the real parts of the non-real endpoints are kept when they are roots
  (positive dimensional components are real affine subspaces).
 */

// A polynomial in d variables as a list of (exponents, coefficient).
struct PCPoly {
  std::vector<std::vector<int>> expo;
  std::vector<double> coef;
};

inline PCPoly PC_PolyMult(PCPoly const &P, PCPoly const &Q, int d) {
  std::map<std::vector<int>, double> M;
  for (size_t i = 0; i < P.expo.size(); i++)
    for (size_t j = 0; j < Q.expo.size(); j++) {
      std::vector<int> e(d);
      for (int u = 0; u < d; u++)
        e[u] = P.expo[i][u] + Q.expo[j][u];
      M[e] += P.coef[i] * Q.coef[j];
    }
  PCPoly R;
  for (auto &kv : M) {
    R.expo.push_back(kv.first);
    R.coef.push_back(kv.second);
  }
  return R;
}

struct PCRoot {
  std::vector<double> x;
  bool singular;
};

// Solve F_i = 0 (i < d), each F_i of degree <= deg, by the projective
// homotopy. Returns the real roots and counts the paths stopped in the
// middle.
inline std::vector<PCRoot> PC_SolveSystem(std::vector<PCPoly> const &F,
                                          int d, int deg,
                                          size_t &n_stalled_mid) {
  using C = std::complex<double>;
  using CV = Eigen::Matrix<C, Eigen::Dynamic, 1>;
  using CM = Eigen::Matrix<C, Eigen::Dynamic, Eigen::Dynamic>;
  int N = d + 1;
  // Homogenized evaluation: monomial x^e -> X_0^{deg - |e|} X^e.
  auto evalFh = [&](CV const &X, CV &val, CM &jac) {
    val = CV::Zero(d);
    jac = CM::Zero(d, N);
    for (int i = 0; i < d; i++) {
      for (size_t t = 0; t < F[i].expo.size(); t++) {
        std::vector<int> const &e = F[i].expo[t];
        int tot = 0;
        for (int u = 0; u < d; u++)
          tot += e[u];
        int e0 = deg - tot;
        std::vector<int> full(N);
        full[0] = e0;
        for (int u = 0; u < d; u++)
          full[u + 1] = e[u];
        C mono = F[i].coef[t];
        for (int u = 0; u < N; u++)
          for (int r = 0; r < full[u]; r++)
            mono *= X(u);
        val(i) += mono;
        for (int u = 0; u < N; u++) {
          if (full[u] == 0)
            continue;
          C der = F[i].coef[t] * double(full[u]);
          for (int w = 0; w < N; w++) {
            int p = full[w] - (w == u ? 1 : 0);
            for (int r = 0; r < p; r++)
              der *= X(w);
          }
          jac(i, u) += der;
        }
      }
    }
  };
  auto evalAffine = [&](CV const &x, CV &val, CM &jac) {
    CV X(N);
    X(0) = 1;
    X.tail(d) = x;
    CM J;
    evalFh(X, val, J);
    jac = J.rightCols(d);
  };
  CV chart(N);
  for (int i = 0; i < N; i++)
    chart(i) = C(std::cos(0.7 * i + 0.3), std::sin(1.3 * i + 0.1));
  std::vector<C> ListGamma{{0.6129831, 0.7900987},
                           {-0.4161468, 0.9092974},
                           {0.9800666, -0.1986693},
                           {-0.8775826, -0.4794255}};
  auto attempt = [&](C gamma, size_t &mid) {
    std::vector<PCRoot> roots;
    // Start points: X_0 = 1, X_i cube (deg-1) roots of unity.
    int k = deg;
    std::vector<C> unity(k);
    for (int r = 0; r < k; r++)
      unity[r] = std::polar(1.0, 2 * M_PI * r / k);
    size_t n_start = 1;
    for (int i = 0; i < d; i++)
      n_start *= k;
    auto evalH = [&](CV const &X, double t, CV &H, CM &J) {
      CV Fv;
      CM JF;
      evalFh(X, Fv, JF);
      CV G(d);
      CM JG = CM::Zero(d, N);
      for (int i = 0; i < d; i++) {
        G(i) = std::pow(X(i + 1), k) - std::pow(X(0), k);
        JG(i, i + 1) = double(k) * std::pow(X(i + 1), k - 1);
        JG(i, 0) = -double(k) * std::pow(X(0), k - 1);
      }
      H = CV(N);
      H.head(d) = (1 - t) * gamma * G + t * Fv;
      H(d) = chart.dot(X) - 1.0;
      J = CM(N, N);
      J.topRows(d) = (1 - t) * gamma * JG + t * JF;
      for (int i = 0; i < N; i++)
        J(d, i) = std::conj(chart(i));
    };
    auto evalHt = [&](CV const &X) {
      CV Fv;
      CM JF;
      evalFh(X, Fv, JF);
      CV Ht = CV::Zero(N);
      for (int i = 0; i < d; i++)
        Ht(i) = Fv(i) - gamma * (std::pow(X(i + 1), k) - std::pow(X(0), k));
      return Ht;
    };
    for (size_t sidx = 0; sidx < n_start; sidx++) {
      CV X(N);
      X(0) = 1;
      size_t rem = sidx;
      for (int i = 0; i < d; i++) {
        X(i + 1) = unity[rem % k];
        rem /= k;
      }
      X /= chart.dot(X);
      double t = 0, dt = 0.02;
      int n_step = 0;
      while (t < 1 && n_step < 200000) {
        n_step++;
        double h = std::min(dt, 1 - t);
        auto deriv = [&](CV const &Y, double tt) -> CV {
          CV H;
          CM J;
          evalH(Y, tt, H, J);
          return -J.fullPivLu().solve(evalHt(Y));
        };
        CV k1 = deriv(X, t);
        CV k2 = deriv(X + 0.5 * h * k1, t + 0.5 * h);
        CV k3 = deriv(X + 0.5 * h * k2, t + 0.5 * h);
        CV k4 = deriv(X + h * k3, t + h);
        CV Y = X + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        double tn = t + h;
        double tol = tn > 0.5 ? 1e-8 : 1e-11;
        bool ok = false;
        for (int it = 0; it < 8; it++) {
          CV H;
          CM J;
          evalH(Y, tn, H, J);
          CV delta = J.fullPivLu().solve(-H);
          Y += delta;
          if (delta.norm() < tol * (1 + Y.norm())) {
            ok = true;
            break;
          }
        }
        if (!ok || !(Y.norm() < 1e12)) {
          dt /= 2;
          if (dt < 1e-13)
            break;
          continue;
        }
        X = Y;
        t = tn;
        dt = std::min(2 * dt, 0.05);
      }
      if (t <= 0.9)
        mid++;
      if (std::abs(X(0)) < 1e-8 * X.norm())
        continue;
      CV x = X.tail(d) / X(0);
      for (int it = 0; it < 30; it++) {
        CV val;
        CM jac;
        evalAffine(x, val, jac);
        CV delta = jac.completeOrthogonalDecomposition().solve(-val);
        x += delta;
        if (delta.norm() < 1e-14 * (1 + x.norm()))
          break;
      }
      if (x.norm() > 1e8)
        continue;
      double mi = 0;
      for (int i = 0; i < d; i++)
        mi = std::max(mi, std::fabs(x(i).imag()));
      CV val;
      CM jac;
      if (mi > 1e-6 * (1 + x.norm())) {
        CV xr(d);
        for (int i = 0; i < d; i++)
          xr(i) = x(i).real();
        evalAffine(xr, val, jac);
        if (val.norm() > 1e-7 * (1 + xr.squaredNorm()))
          continue;
        x = xr;
      }
      evalAffine(x, val, jac);
      if (val.norm() > 1e-7 * (1 + x.squaredNorm()))
        continue;
      PCRoot r;
      for (int i = 0; i < d; i++)
        r.x.push_back(x(i).real());
      Eigen::JacobiSVD<CM> svd(jac);
      r.singular = svd.singularValues()(d - 1) <
                   1e-7 * (1 + svd.singularValues()(0));
      roots.push_back(r);
    }
    return roots;
  };
  std::vector<PCRoot> best;
  size_t best_mid = std::numeric_limits<size_t>::max();
  for (auto &g : ListGamma) {
    size_t mid = 0;
    auto R = attempt(g, mid);
    if (mid < best_mid) {
      best = R;
      best_mid = mid;
    }
    if (best_mid == 0)
      break;
  }
  n_stalled_mid += best_mid;
  return best;
}

template <typename T> struct PCCandidate {
  int i_leaf;
  double value;  // |det W| / index
  double min_gauge;
  bool singular;
  MyMatrix<double> W;
};

// The critical points of det on the affine space of a leaf, with the test
// points (the contact set) in D, and their admissibility.
template <typename T>
std::vector<PCCandidate<T>>
PC_Stage3(PCBody<T> const &body, PCClass const &C,
          std::vector<std::pair<MyVector<T>, T>> const &eqs, int i_leaf,
          double value_bound, size_t &n_stalled_mid) {
  int n = body.n;
  int nn = n * n;
  int r = eqs.size();
  MyMatrix<T> A(r, nn);
  MyVector<T> b(r);
  for (int i = 0; i < r; i++) {
    for (int u = 0; u < nn; u++)
      A(i, u) = eqs[i].first(u);
    b(i) = eqs[i].second;
  }
  MyMatrix<T> AT = TransposedMat(A);
  std::optional<MyVector<T>> opt = SolutionMat(AT, b);
  std::vector<PCCandidate<T>> out;
  if (!opt)
    return out;
  MyMatrix<T> NSP = NullspaceTrMat(A);
  int d = NSP.rows();
  MyVector<double> w0 = UniversalVectorConversion<double, T>(*opt);
  MyMatrix<double> M = UniversalMatrixConversion<double, T>(NSP);
  // Scale the directions.
  double nrm0 = std::max(1.0, w0.cwiseAbs().maxCoeff());
  for (int i = 0; i < d; i++) {
    double nr = M.row(i).cwiseAbs().maxCoeff();
    if (nr > 0)
      M.row(i) *= nrm0 / nr;
  }
  auto get_W = [&](std::vector<double> const &x) {
    MyMatrix<double> W(n, n);
    for (int a = 0; a < n; a++)
      for (int c = 0; c < n; c++) {
        double v = w0(a * n + c);
        for (int i = 0; i < d; i++)
          v += x[i] * M(i, a * n + c);
        W(a, c) = v;
      }
    return W;
  };
  std::vector<std::vector<double>> roots;
  std::vector<bool> sing;
  if (d == 0) {
    roots.push_back({});
    sing.push_back(false);
  } else {
    // det as a polynomial in the d parameters.
    std::vector<std::vector<PCPoly>> Wp(n, std::vector<PCPoly>(n));
    for (int a = 0; a < n; a++)
      for (int c = 0; c < n; c++) {
        PCPoly P;
        P.expo.push_back(std::vector<int>(d, 0));
        P.coef.push_back(w0(a * n + c));
        for (int i = 0; i < d; i++) {
          std::vector<int> e(d, 0);
          e[i] = 1;
          P.expo.push_back(e);
          P.coef.push_back(M(i, a * n + c));
        }
        Wp[a][c] = P;
      }
    std::vector<int> perm(n);
    for (int i = 0; i < n; i++)
      perm[i] = i;
    std::map<std::vector<int>, double> detm;
    do {
      int sg = 1;
      for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
          if (perm[i] > perm[j])
            sg = -sg;
      PCPoly P = Wp[0][perm[0]];
      for (int a = 1; a < n; a++)
        P = PC_PolyMult(P, Wp[a][perm[a]], d);
      for (size_t t = 0; t < P.expo.size(); t++)
        detm[P.expo[t]] += sg * P.coef[t];
    } while (std::next_permutation(perm.begin(), perm.end()));
    std::vector<PCPoly> grad(d);
    for (auto &kv : detm)
      for (int i = 0; i < d; i++)
        if (kv.first[i] > 0) {
          std::vector<int> e = kv.first;
          e[i]--;
          grad[i].expo.push_back(e);
          grad[i].coef.push_back(kv.second * kv.first[i]);
        }
    for (int i = 0; i < d; i++) {
      double mx = 0;
      for (auto &c : grad[i].coef)
        mx = std::max(mx, std::fabs(c));
      if (mx > 0)
        for (auto &c : grad[i].coef)
          c /= mx;
    }
    auto R = PC_SolveSystem(grad, d, n - 1, n_stalled_mid);
    for (auto &rt : R) {
      roots.push_back(rt.x);
      sing.push_back(rt.singular);
    }
  }
  int64_t index = C.elts.size();
  for (size_t ir = 0; ir < roots.size(); ir++) {
    MyMatrix<double> W = get_W(roots[ir]);
    double value = std::fabs(W.determinant()) / double(index);
    if (value < 1e-9 || value > value_bound)
      continue;
    // Admissibility of Lambda = L W: the minimum of the gauge over the
    // nonzero points (g / E + z) W, z in a box.
    MyMatrix<double> Winv = W.inverse();
    double R = body.r_max * 1.0001;
    std::vector<int> bnd(n);
    bool too_big = false;
    for (int i = 0; i < n; i++) {
      bnd[i] = static_cast<int>(std::floor(R * Winv.col(i).norm())) + 2;
      if (bnd[i] > 60)
        too_big = true;
    }
    double mg = std::numeric_limits<double>::max();
    if (too_big) {
      mg = -1;
    } else {
      for (auto const &g : C.elts) {
        std::vector<int> z(n);
        for (int i = 0; i < n; i++)
          z[i] = -bnd[i];
        while (true) {
          std::vector<double> c(n);
          bool zero = true;
          for (int i = 0; i < n; i++) {
            c[i] = double(g[i]) / double(C.E) + z[i];
            if (std::fabs(c[i]) > 1e-12)
              zero = false;
          }
          if (!zero) {
            std::vector<double> x(n, 0);
            for (int a = 0; a < n; a++)
              for (int b2 = 0; b2 < n; b2++)
                x[b2] += c[a] * W(a, b2);
            double gmax = -1e300;
            for (int j = 0; j < body.m; j++) {
              double sc = 0;
              for (int b2 = 0; b2 < n; b2++)
                sc += body.Lmat_d(j, b2) * x[b2];
              gmax = std::max(gmax, sc);
            }
            mg = std::min(mg, gmax);
          }
          int i = 0;
          while (i < n && z[i] == bnd[i])
            z[i++] = -bnd[i - 1 + 1 - 1 + 1 - 1];
          if (i == n)
            break;
          z[i]++;
        }
      }
    }
    out.push_back({i_leaf, value, mg, static_cast<bool>(sing[ir]), W});
  }
  return out;
}

// clang-format off
#endif  // SRC_POLYNORM_POLYNORM_PACKINGCONTACTS_H_
// clang-format on

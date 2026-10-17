// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLYNORM_POLYNORM_KISSING_H_
#define SRC_POLYNORM_POLYNORM_KISSING_H_

// clang-format off
#include "PolyNorm_Basic.h"
#include "PolyNorm_PackingContacts.h"
#include "POLY_PolytopeInt.h"
#include "ClassicLLL.h"
#include <algorithm>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_POLYNORM_KISSING
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_POLYNORM_KISSING
#endif

/*
  The maximal kissing number of a polytope P of dimension n: the maximum,
  over the lattices Lambda with P + Lambda a packing, of the number of
  translates P + v, v in Lambda nonzero, touching P. The translates P and
  P + v touch exactly when v is on the boundary of the difference body
  D = P - P, so the problem is the maximum number of points of an
  admissible lattice (no nonzero point in the interior of D) on the
  boundary of D. The value is scale invariant, so no packing scalar enters:
  the lattice is scaled instead of the polytope.

  The enumeration is the one of PolyNorm_PackingContacts.h without the
  criticality. The maximum is attained by a lattice whose contact set S
  spans R^n: if S spans a proper subspace V, contracting the lattice
  transversally to V fixes S and moves the other points inward, and at the
  largest admissible contraction a new contact appears. So with b_1, ...,
  b_n in S independent of maximal |det|, the points (1) and (2) of the
  packing enumeration hold for any admissible lattice: the lattice is one
  of the classes of lattice cross-polytopes in the coordinates of the b_i,
  the contacts have coordinates in [-1,1]^n and every n of them have
  |det| <= 1. The points (3) and (4) use the criticality and do not apply.
  A configuration (L, S, F) of a class, a candidate contact set and a
  choice of a facet of D containing each contact is searched by the same
  backtracking, maximizing |S| instead of collecting affine spaces, with
  the same linear program, determinant and additive-relation pruning.

  A configuration larger than the best witness is then decided exactly:
  --- the feasibility is a linear program over the n^2 entries of the basis
      W: each contact s W on its facet and in D (all the facet
      inequalities, not only the neighbors). The feasible set Q is a
      bounded polyhedron since the rows of W are the contacts e_i W, which
      lie in D.
  --- the admissibility (some W in Q whose lattice L W has no nonzero
      point in the interior of D) is a branch and bound: at a feasible W,
      a nonzero c in L with ||c W||_D < 1 is searched, by the LLL
      reduction of the Gram matrix, by the rational kernel when
      det(W) = 0, and by the exact enumeration of the lattice points of D
      otherwise. If there is none, W is a witness and its boundary points
      are counted (the enumeration returns them). Otherwise the subtree
      branches on the facet of D that keeps c outside: the constraint
      ell_j(c W) >= 1 for each facet j, and a branch contains the witness
      lattice of the configuration if it has one. The node budget bounds
      the branch and bound; a configuration exhausting it is reported as
      undecided, so that the certified value and the remaining upper bound
      are both returned.
 */

template <typename T> struct KissingWitness {
  int i_class;
  // The number of nonzero lattice vectors on the boundary of D: the
  // kissing number of the witness lattice (contacts come in opposite
  // pairs, the count is even).
  int n_contact;
  // A basis of the witness lattice, as rows.
  MyMatrix<T> LattBasis;
  // The contact vectors: the differences between touching translates.
  std::vector<MyVector<T>> ListContact;
  // |det| of the lattice (the packing density of P + Lambda is
  // vol(P) / det).
  T det;
};

template <typename T> struct KissingResult {
  // The certified value: realized by the witness, and maximal when
  // n_undecided = 0.
  int n_contact;
  // The upper bound including the undecided configurations; equal to
  // n_contact when every configuration was decided.
  int n_contact_possible;
  size_t n_undecided;
  std::optional<KissingWitness<T>> witness;
};

// A configuration whose admissibility branch and bound ran out of budget.
struct KissUndecided {
  int i_class;
  std::vector<std::vector<int64_t>> pts;
  std::vector<int> facets;
};

template <typename T> struct KissingSearch {
  PCBody<T> const &body;
  std::vector<PCClass> const &ListClass;
  std::ostream &os;
  int n;
  // The packing-contacts search, for its relaxed linear program (Feasible)
  // and the facet-sum pruning (Triple).
  PCSearch<T> aux;
  // The branch and bound budget: a small slice during the enumeration, the
  // full budget to re-decide the surviving configurations afterwards.
  int64_t adm_budget_per_leaf;
  int best;
  std::optional<KissingWitness<T>> witness;
  std::vector<KissUndecided> ListUndecided;
  // Statistics.
  size_t n_node, n_leaf, n_adm_lp;

  KissingSearch(PCBody<T> const &_body, std::vector<PCClass> const &_ListClass,
                int64_t _adm_budget_per_leaf, std::ostream &_os)
      : body(_body), ListClass(_ListClass), os(_os), n(_body.n),
        aux(_body, _ListClass, _os),
        adm_budget_per_leaf(_adm_budget_per_leaf), best(0), n_node(0),
        n_leaf(0), n_adm_lp(0) {}

  // The row of ell_f(p W) over the n^2 entries of W.
  MyVector<T> RowPoint(MyVector<T> const &p, int f) const {
    MyVector<T> R(n * n);
    for (int a = 0; a < n; a++)
      for (int b = 0; b < n; b++)
        R(a * n + b) = p(a) * body.Lmat(f, b);
    return R;
  }

  // A basis of the lattice L of the class, as rows.
  MyMatrix<T> ClassBasis(int i_class) const {
    PCClass const &C = ListClass[i_class];
    int n_elt = C.elts.size();
    MyMatrix<T> Gens(n + n_elt, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        Gens(i, j) = T(i == j ? 1 : 0);
    for (int i = 0; i < n_elt; i++)
      for (int j = 0; j < n; j++)
        Gens(n + i, j) =
            T(static_cast<long>(C.elts[i][j])) / T(static_cast<long>(C.E));
    return GetZbasis(Gens);
  }

  // A feasible point of {W : c0 + R . W >= 0 for the rows}, towards the
  // relative interior: the midpoint of the minimizer and the maximizer of
  // the sum of the entries (the polyhedron is bounded, the rows of W are
  // in D).
  std::optional<MyMatrix<T>> SolveQ(std::vector<MyVector<T>> const &rows) {
    int nn = n * n;
    MyMatrix<T> M(rows.size(), nn + 1);
    for (size_t i = 0; i < rows.size(); i++)
      for (int u = 0; u <= nn; u++)
        M(i, u) = rows[i](u);
    MyVector<T> obj = ZeroVector<T>(nn + 1);
    for (int u = 0; u < nn; u++)
      obj(u + 1) = T(1);
    n_adm_lp++;
    LpSolution<T> sol1 = SIMPLEX_LinearProgramming(M, obj, os);
    if (!sol1.DirectSolution)
      return {};
    MyVector<T> w = *sol1.DirectSolution;
    for (int u = 0; u < nn; u++)
      obj(u + 1) = T(-1);
    n_adm_lp++;
    LpSolution<T> sol2 = SIMPLEX_LinearProgramming(M, obj, os);
    if (sol2.DirectSolution) {
      MyVector<T> w2 = *sol2.DirectSolution;
      for (int u = 0; u < nn; u++)
        w(u) = (w(u) + w2(u)) / T(2);
    }
    MyMatrix<T> W(n, n);
    for (int a = 0; a < n; a++)
      for (int b = 0; b < n; b++)
        W(a, b) = w(a * n + b);
    return W;
  }

  // Either a nonzero c in L with c W in the interior of D, or (when there
  // is none) the contact vectors c W, c in L nonzero, on the boundary.
  std::optional<MyVector<T>> FindViolator(MyMatrix<T> const &B,
                                          MyMatrix<T> const &W,
                                          std::vector<MyVector<T>> &contacts) {
    using Tint = typename underlying_z_ring<T>::ring_type;
    contacts.clear();
    MyMatrix<T> M = B * W;
    T det = DeterminantMat(M);
    auto comb = [&](MyVector<T> const &u, MyMatrix<T> const &Mat) {
      MyVector<T> c = ZeroVector<T>(n);
      for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++)
          c(b) += u(a) * Mat(a, b);
      return c;
    };
    if (det == 0) {
      // A rational kernel vector scaled to Z^n, which is contained in L:
      // its lattice point is 0, in the interior of D.
      MyMatrix<T> NSP = NullspaceTrMat(TransposedMat(M));
      MyVector<T> u = RemoveFractionVector(GetMatrixRow(NSP, 0));
      return comb(u, B);
    }
    // The LLL reduction of the Gram matrix: a short basis vector inside D
    // is a violator, and when none is, the determinant is bounded below
    // and the enumeration below is bounded.
    MyMatrix<T> Gram = M * TransposedMat(M);
    LLLreduction<T, Tint> lll = LLLreducedBasis<T, Tint>(Gram, os);
    MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(lll.Pmat);
    MyMatrix<T> Mred = P_T * M;
    for (int i = 0; i < n; i++) {
      MyVector<T> x = GetMatrixRow(Mred, i);
      if (PolyNorm_Gauge(body.Lmat, x) < 1)
        return comb(GetMatrixRow(P_T, i), B);
    }
    // The lattice points of D for the lattice of basis M: the points u M,
    // u integral, with every facet inequality at most 1.
    int m = body.m;
    MyMatrix<T> G(m, n);
    for (int j = 0; j < m; j++)
      for (int i = 0; i < n; i++) {
        T s(0);
        for (int k = 0; k < n; k++)
          s += body.Lmat(j, k) * M(i, k);
        G(j, i) = s;
      }
    MyMatrix<T> FAC(m, n + 1);
    for (int j = 0; j < m; j++) {
      FAC(j, 0) = T(1);
      for (int i = 0; i < n; i++)
        FAC(j, i + 1) = -G(j, i);
    }
    std::optional<MyVector<T>> viol;
    auto f_insert = [&](MyVector<Tint> const &u) -> bool {
      bool zero = true;
      for (int i = 0; i < n; i++)
        if (u(i) != 0)
          zero = false;
      if (zero)
        return true;
      MyVector<T> u_T = UniversalVectorConversion<T, Tint>(u);
      T gauge(0);
      for (int j = 0; j < m; j++) {
        T s(0);
        for (int i = 0; i < n; i++)
          s += G(j, i) * u_T(i);
        if (j == 0 || s > gauge)
          gauge = s;
      }
      if (gauge < 1) {
        viol = comb(u_T, B);
        return false;
      }
      contacts.push_back(comb(u_T, M));
      return true;
    };
    Kernel_GetListIntegralPoint_ITER<T, Tint, decltype(f_insert)>(FAC, f_insert,
                                                                  os);
    if (viol) {
      contacts.clear();
      return viol;
    }
    return {};
  }

  // The admissibility branch and bound on the polyhedron of the rows.
  // Returns 0 when no admissible W exists in it, 1 when a witness was
  // found (Wout and contacts are then set), 2 when the budget ran out.
  int AdmSearch(std::vector<MyVector<T>> &rows, MyMatrix<T> const &B,
                int64_t &budget, MyMatrix<T> &Wout,
                std::vector<MyVector<T>> &contacts) {
    if (budget <= 0)
      return 2;
    budget--;
    std::optional<MyMatrix<T>> optW = SolveQ(rows);
    if (!optW)
      return 0;
    MyMatrix<T> const &W = *optW;
    std::vector<MyVector<T>> cont;
    std::optional<MyVector<T>> viol = FindViolator(B, W, cont);
    if (!viol) {
      Wout = W;
      contacts = std::move(cont);
      return 1;
    }
    MyVector<T> const &c = *viol;
    // c W is in the interior of D: branch on the facet pushed to at least
    // 1, the nearest facets first.
    MyVector<T> x = ZeroVector<T>(n);
    for (int a = 0; a < n; a++)
      for (int b = 0; b < n; b++)
        x(b) += c(a) * W(a, b);
    std::vector<std::pair<T, int>> order;
    for (int j = 0; j < body.m; j++) {
      T s(0);
      for (int b = 0; b < n; b++)
        s += body.Lmat(j, b) * x(b);
      order.push_back({s, j});
    }
    std::sort(order.begin(), order.end(),
              [](auto const &a, auto const &b) { return a.first > b.first; });
    int nn = n * n;
    for (auto const &of : order) {
      MyVector<T> R = RowPoint(c, of.second);
      MyVector<T> X(nn + 1);
      X(0) = T(-1);
      for (int u = 0; u < nn; u++)
        X(u + 1) = R(u);
      rows.push_back(X);
      int st = AdmSearch(rows, B, budget, Wout, contacts);
      rows.pop_back();
      if (st != 0)
        return st;
    }
    return 0;
  }

  // The exact decision of a configuration with the given budget, and the
  // update of the best witness. Returns the status of AdmSearch.
  int DecideLeaf(int i_class, std::vector<std::vector<int64_t>> const &pts,
                 std::vector<int> const &facets, int64_t budget_leaf) {
    n_leaf++;
    int64_t E = ListClass[i_class].E;
    int nn = n * n;
    std::vector<MyVector<T>> rows;
    auto push = [&](T const &c0, MyVector<T> const &R) {
      MyVector<T> X(nn + 1);
      X(0) = c0;
      for (int u = 0; u < nn; u++)
        X(u + 1) = R(u);
      rows.push_back(X);
    };
    for (size_t i = 0; i < pts.size(); i++) {
      MyVector<T> p(n);
      for (int a = 0; a < n; a++)
        p(a) = T(static_cast<long>(pts[i][a])) / T(static_cast<long>(E));
      for (int j = 0; j < body.m; j++) {
        MyVector<T> R = RowPoint(p, j);
        push(T(1), -R);
        if (j == facets[i])
          push(T(-1), R);
      }
    }
    MyMatrix<T> B = ClassBasis(i_class);
    int64_t budget = budget_leaf;
    MyMatrix<T> W;
    std::vector<MyVector<T>> contacts;
    int st = AdmSearch(rows, B, budget, W, contacts);
    if (st == 1) {
      int cnt = contacts.size();
#ifdef SANITY_CHECK_POLYNORM_KISSING
      if (cnt % 2 != 0 || cnt < 2 * static_cast<int>(pts.size())) {
        std::cerr << "KISSING: the witness has " << cnt << " contacts, the "
                  << "configuration has " << pts.size() << " pairs\n";
        throw TerminalException{1};
      }
      for (size_t i = 0; i < pts.size(); i++) {
        MyVector<T> x = ZeroVector<T>(n);
        for (int a = 0; a < n; a++)
          for (int b = 0; b < n; b++)
            x(b) += T(static_cast<long>(pts[i][a])) /
                    T(static_cast<long>(E)) * W(a, b);
        if (PolyNorm_Gauge(body.Lmat, x) != 1) {
          std::cerr << "KISSING: a configured contact is not on the "
                    << "boundary of D\n";
          throw TerminalException{1};
        }
      }
#endif
      if (cnt > best) {
        best = cnt;
        MyMatrix<T> LattBasis = B * W;
        T det = T_abs(DeterminantMat(LattBasis));
        witness = KissingWitness<T>{i_class, cnt, std::move(LattBasis),
                                    std::move(contacts), det};
#ifdef DEBUG_POLYNORM_KISSING
        os << "KISSING: new best " << best << " contacts, class " << i_class
           << ", " << pts.size() << " configured pairs\n";
#endif
      }
    }
    return st;
  }

  // The enumeration-time decision: a slice of the budget, the exhausted
  // configurations are kept for the second phase.
  void Leaf(int i_class, std::vector<std::vector<int64_t>> const &pts,
            std::vector<int> const &facets) {
    int64_t slice = std::max<int64_t>(adm_budget_per_leaf / 32, 250);
    int st = DecideLeaf(i_class, pts, facets, slice);
    if (st == 2) {
      ListUndecided.push_back({i_class, pts, facets});
#ifdef DEBUG_POLYNORM_KISSING
      os << "KISSING: undecided configuration of " << 2 * pts.size()
         << " contacts in class " << i_class << "\n";
#endif
    }
  }

  // The second phase: the surviving configurations larger than the best
  // witness, decided with the full budget, the largest first.
  void ResolveUndecided() {
    std::sort(ListUndecided.begin(), ListUndecided.end(),
              [](KissUndecided const &a, KissUndecided const &b) {
                return a.pts.size() > b.pts.size();
              });
    std::vector<KissUndecided> remaining;
    for (auto &u : ListUndecided) {
      if (2 * static_cast<int>(u.pts.size()) <= best)
        continue;
      int st = DecideLeaf(u.i_class, u.pts, u.facets, adm_budget_per_leaf);
      if (st == 2)
        remaining.push_back(std::move(u));
    }
    ListUndecided = std::move(remaining);
  }

  // The lattice Z^n scaled to its first minimum: an admissible lattice
  // whose contacts seed the search.
  void SeedZn() {
    using Tint = typename underlying_z_ring<T>::ring_type;
    int ic0 = -1;
    for (size_t ic = 0; ic < ListClass.size(); ic++)
      if (ListClass[ic].elts.size() == 1)
        ic0 = static_cast<int>(ic);
    if (ic0 < 0)
      return;
    // The first minimum: the smallest gauge of a nonzero z in Z^n, among
    // the lattice points of t0 D for the gauge t0 of a unit vector.
    MyVector<T> e = ZeroVector<T>(n);
    e(0) = T(1);
    T t0 = PolyNorm_Gauge(body.Lmat, e);
    std::vector<MyVector<Tint>> ListZ =
        PolyNorm_LatticePoints<T, Tint>(body.Lmat, t0, os);
    T m0 = t0;
    for (auto &z : ListZ) {
      if (IsZeroVector(z))
        continue;
      MyVector<T> z_T = UniversalVectorConversion<T, Tint>(z);
      T g = PolyNorm_Gauge(body.Lmat, z_T);
      if (g < m0)
        m0 = g;
    }
    MyMatrix<T> W = ZeroMatrix<T>(n, n);
    for (int i = 0; i < n; i++)
      W(i, i) = T(1) / m0;
    MyMatrix<T> B = IdentityMat<T>(n);
    std::vector<MyVector<T>> contacts;
    std::optional<MyVector<T>> viol = FindViolator(B, W, contacts);
#ifdef SANITY_CHECK_POLYNORM_KISSING
    if (viol) {
      std::cerr << "KISSING: the seed lattice Z^n / m0 is not admissible\n";
      throw TerminalException{1};
    }
#else
    if (viol)
      return;
#endif
    int cnt = contacts.size();
    if (cnt > best) {
      best = cnt;
      T det = T_abs(DeterminantMat(W));
      witness =
          KissingWitness<T>{ic0, cnt, std::move(W), std::move(contacts), det};
#ifdef DEBUG_POLYNORM_KISSING
      os << "KISSING: seed Z^n gives " << best << " contacts\n";
#endif
    }
  }

  // The search for one class, the facet of the first contact being f0 (one
  // representative per orbit under the symmetries of D).
  void SearchClass(int i_class, int f0) {
    PCClass const &C = ListClass[i_class];
    int64_t E = C.E;
    std::vector<std::vector<int64_t>> cand = PC_Candidates(C, n);
    int nc = cand.size();
    int64_t det_bound = 1;
    for (int i = 0; i < n; i++)
      det_bound *= E;
    // Membership in the scaled lattice E L, by the residues mod E.
    std::set<std::vector<int64_t>> residues;
    for (auto const &g : C.elts) {
      std::vector<int64_t> r(n);
      for (int u = 0; u < n; u++) {
        int64_t x = g[u] % E;
        if (x < 0)
          x += E;
        r[u] = x;
      }
      residues.insert(r);
    }
    auto in_lattice = [&](std::vector<int64_t> const &v) -> bool {
      std::vector<int64_t> r(n);
      for (int u = 0; u < n; u++) {
        int64_t x = v[u] % E;
        if (x < 0)
          x += E;
        r[u] = x;
      }
      return residues.count(r) > 0;
    };
    std::map<std::vector<int64_t>, int> cand_index;
    for (int i = 0; i < nc; i++)
      cand_index[cand[i]] = i;
    // A candidate v with v / t in the lattice (t >= 2) is over an interior
    // lattice point (gauge 1 / t) and is never a contact. The unit vectors
    // are primitive: a fraction of one would be interior to the
    // cross-polytope of the class.
    std::vector<uint8_t> dead(nc, 0);
    for (int i = n; i < nc; i++) {
      int64_t g = 0;
      for (int u = 0; u < n; u++)
        g = std::gcd(g, std::llabs(cand[i][u]));
      for (int64_t t = 2; t <= g && !dead[i]; t++) {
        if (g % t != 0)
          continue;
        std::vector<int64_t> q(n);
        for (int u = 0; u < n; u++)
          q[u] = cand[i][u] / t;
        if (in_lattice(q))
          dead[i] = 1;
      }
    }
    // The closure of the contact set: for contacts a, b, a lattice point
    // (a +- b) / 2 is on the boundary (in D by convexity, not interior by
    // admissibility), so it is a contact too: pair_mid lists the forced
    // candidates of a pair, rev_mid the pairs forcing a candidate. A
    // lattice point (a +- b) / t with t >= 3 has gauge at most 2 / t < 1,
    // which is interior: the pair is forbidden (pair_forb).
    std::vector<std::vector<int>> pair_mid(nc * nc);
    std::vector<uint8_t> pair_forb(nc * nc, 0);
    std::vector<std::vector<std::pair<int, int>>> rev_mid(nc);
    for (int i = 0; i < nc; i++)
      for (int j = i + 1; j < nc; j++) {
        for (int s : {1, -1}) {
          std::vector<int64_t> w(n);
          int64_t g = 0;
          for (int u = 0; u < n; u++) {
            w[u] = cand[i][u] + s * cand[j][u];
            g = std::gcd(g, std::llabs(w[u]));
          }
          for (int64_t t = 2; t <= g; t++) {
            if (g % t != 0)
              continue;
            std::vector<int64_t> q(n);
            for (int u = 0; u < n; u++)
              q[u] = w[u] / t;
            if (!in_lattice(q))
              continue;
            if (t == 2) {
              int64_t f = 0;
              for (int u = 0; u < n; u++)
                if (q[u] != 0) {
                  f = q[u];
                  break;
                }
              if (f < 0)
                for (auto &x : q)
                  x = -x;
              int idx = cand_index.at(q);
              pair_mid[i * nc + j].push_back(idx);
              rev_mid[idx].push_back({i, j});
            } else {
              pair_forb[i * nc + j] = 1;
            }
          }
        }
      }
    // The alive candidates from each position on, for the bound.
    std::vector<int> rem_alive(nc + 1, 0);
    for (int k = nc - 1; k >= 0; k--)
      rem_alive[k] = rem_alive[k + 1] + (dead[k] ? 0 : 1);
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
    std::vector<int> in_facet(nc, -1);
    std::vector<std::vector<int64_t>> pts;
    std::vector<int> facets;
    std::function<void(int)> rec = [&](int k) {
      n_node++;
      if (2 * (static_cast<int>(pts.size()) + rem_alive[k]) <= best)
        return;
      if (k == nc) {
        if (2 * static_cast<int>(pts.size()) > best)
          Leaf(i_class, pts, facets);
        return;
      }
      // The closure of the contact set with the earlier candidates,
      // independent of the facet: a forbidden pair, or a forced midpoint
      // already excluded, rules the inclusion of k out.
      bool closure_ok = !dead[k];
      for (int i = 0; i < k && closure_ok; i++) {
        if (in_facet[i] < 0)
          continue;
        int a = std::min(i, k), b = std::max(i, k);
        if (pair_forb[a * nc + b]) {
          closure_ok = false;
          break;
        }
        for (int idx : pair_mid[a * nc + b])
          if (idx < k && in_facet[idx] < 0) {
            closure_ok = false;
            break;
          }
      }
      // Determinant condition with every n-1 earlier points, independent
      // of the facet.
      bool det_ok = closure_ok;
      int np = pts.size();
      if (det_ok && np >= n - 1) {
        std::vector<int> sub(n - 1);
        for (int u = 0; u < n - 1; u++)
          sub[u] = u;
        while (det_ok) {
          std::vector<std::vector<int64_t>> M;
          for (int u = 0; u < n - 1; u++)
            M.push_back(pts[sub[u]]);
          M.push_back(cand[k]);
          if (std::llabs(PCSearch<T>::Det(M)) > det_bound)
            det_ok = false;
          int a = n - 2;
          while (a >= 0 && sub[a] == np - (n - 1) + a)
            a--;
          if (a < 0)
            break;
          sub[a]++;
          for (int b = a + 1; b < n - 1; b++)
            sub[b] = sub[b - 1] + 1;
        }
      }
      // Include on a facet first: the maximum wants the contacts in.
      if (det_ok) {
        for (int f = 0; f < body.m; f++) {
          if (k == 0 && f != f0)
            continue;
          bool ok = true;
          for (auto &r : rel[k]) {
            if (in_facet[r[0]] < 0 || in_facet[r[1]] < 0)
              continue;
            int fa = r[2] == 1 ? in_facet[r[0]] : body.opp[in_facet[r[0]]];
            int fb = r[3] == 1 ? in_facet[r[1]] : body.opp[in_facet[r[1]]];
            if (!aux.Triple(fa, fb, f)) {
              ok = false;
              break;
            }
          }
          if (!ok)
            continue;
          pts.push_back(cand[k]);
          facets.push_back(f);
          in_facet[k] = f;
          if (aux.Feasible(pts, facets, E))
            rec(k + 1);
          in_facet[k] = -1;
          pts.pop_back();
          facets.pop_back();
        }
      }
      // Exclude (not for the unit vectors, which are contacts, nor the
      // forced midpoint of two included contacts).
      if (k >= n) {
        bool excl_ok = true;
        for (auto &pr : rev_mid[k]) {
          if (pr.first < k && pr.second < k && in_facet[pr.first] >= 0 &&
              in_facet[pr.second] >= 0) {
            excl_ok = false;
            break;
          }
        }
        if (excl_ok) {
          in_facet[k] = -1;
          rec(k + 1);
        }
      }
    };
    rec(0);
  }
};

// The maximal kissing number of the polytope of vertices EXT (homogeneous
// coordinates). The budget bounds the admissibility branch and bound of
// each configuration; a configuration exhausting it leaves n_undecided
// nonzero and widens n_contact_possible.
template <typename T, typename Tgroup>
KissingResult<T> ComputeMaxKissing(MyMatrix<T> const &EXT,
                                   int64_t adm_budget_per_leaf,
                                   std::ostream &os) {
  PCBody<T> body = PC_BuildBody<T, Tgroup>(EXT, os);
  int n = body.n;
  std::vector<PCClass> ListClass = PC_LatticeCrossPolytopes(n);
#ifdef DEBUG_POLYNORM_KISSING
  os << "KISSING: n=" << n << " facets of D: " << body.m << " ("
     << body.n_orbit << " orbits), " << ListClass.size() << " classes\n";
#endif
  KissingSearch<T> search(body, ListClass, adm_budget_per_leaf, os);
  search.SeedZn();
  // The bound of Hadwiger: no convex body has more than 3^n - 1 touching
  // translates, so the search stops once that many contacts are realized.
  int hadwiger = 1;
  for (int i = 0; i < n; i++)
    hadwiger *= 3;
  hadwiger -= 1;
  for (size_t ic = 0; ic < ListClass.size() && search.best < hadwiger; ic++)
    for (int f = 0; f < body.m && search.best < hadwiger; f++)
      if (body.is_rep[f])
        search.SearchClass(static_cast<int>(ic), f);
  search.ResolveUndecided();
#ifdef DEBUG_POLYNORM_KISSING
  os << "KISSING: nodes " << search.n_node << ", leaves " << search.n_leaf
     << ", admissibility linear programs " << search.n_adm_lp
     << ", undecided " << search.ListUndecided.size() << "\n";
#endif
  KissingResult<T> res;
  res.n_contact = search.best;
  res.n_contact_possible = search.best;
  res.n_undecided = search.ListUndecided.size();
  for (auto &u : search.ListUndecided)
    res.n_contact_possible = std::max(res.n_contact_possible,
                                      2 * static_cast<int>(u.pts.size()));
  res.witness = std::move(search.witness);
  return res;
}

// clang-format off
#endif  // SRC_POLYNORM_POLYNORM_KISSING_H_
// clang-format on

// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ISOMETRYGLUE_H_
#define SRC_GENERA_ISOMETRYGLUE_H_

// clang-format off
#include "IsometryMass.h"
#include "GenusEnumerationGroup.h"
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ISOMETRY_GLUE
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ISOMETRY_GLUE
#endif

/*
  The seeds of the enumeration of pairs by gluing. Given two pairs (L_1,
  sigma_1) and (L_2, sigma_2), even positive definite, the sigma-stable
  primitive extensions L of L_1 (+) L_2, that is the even lattices

      L_1 (+) L_2 subset L subset L_1^# (+) L_2^#

  stable under sigma = sigma_1 (+) sigma_2 with L cap (L_i (x) Q) = L_i,
  correspond to the subgroups Gamma = L / (L_1 (+) L_2) of D(L_1) (+) D(L_2)
  that are
  - sigma-stable,
  - isotropic for the discriminant quadratic form (L is then even), and
  - meet neither D(L_1) nor D(L_2) (L_i is primitive in L).
  Such a Gamma is the graph of a sigma-equivariant anti-isometry between a
  subgroup of D(L_1) and one of D(L_2), and it splits into its primary
  parts, which are enumerated prime by prime.

  The genus of a pair (L, sigma) obtained in this way depends only on the
  genera of the pairs (L_i, sigma_i) and on the local glue: the
  discriminant forms with the action of sigma are local invariants. So
  gluing fixed representatives of the genera of the (L_i, sigma_i) gives a
  seed in every genus of pairs (L, sigma) whose sublattices L cap V_i lie
  in those genera, and the enumeration from these seeds gives all the
  classes. Gluing more than two pieces is done one piece at a time, since
  L cap (V_1 (+) V_2) is a primitive extension of L_1 (+) L_2.

  Conventions as in GenusEnumerationGroup.h: the isometries g act on the
  columns, g^T A g = A, and internally P = g^T acts on the rows.
 */

struct IsometryGlueBounds {
  // the elements of one primary part of D(L_1) (+) D(L_2)
  int64_t max_module = 1 << 20;
  // the sigma-stable isotropic subgroups of one primary part
  size_t max_subgroup = 100000;
};

// The direct sum of two q-primary modules, the first r1 coordinates being
// those of M1.
inline PrimaryQuadModule DirectSumModule(PrimaryQuadModule const &M1,
                                         PrimaryQuadModule const &M2) {
  PrimaryQuadModule M;
  M.q = M1.q;
  M.N = std::max(M1.N, M2.N);
  M.d = 0;
  size_t r1 = M1.rank(), r2 = M2.rank(), r = r1 + r2;
  M.ord = M1.ord;
  M.ord.insert(M.ord.end(), M2.ord.begin(), M2.ord.end());
  M.bnum.assign(r, std::vector<int64_t>(r, 0));
  M.qnum.assign(r, 0);
  M.sigma_img.assign(r, ModuleElt(r, 0));
  auto insert = [&](PrimaryQuadModule const &Mi, size_t shift) -> void {
    int64_t scale = M.N / Mi.N;
    for (size_t i = 0; i < Mi.rank(); i++) {
      for (size_t j = 0; j < Mi.rank(); j++) {
        M.bnum[shift + i][shift + j] = ModPos64(Mi.bnum[i][j] * scale, M.N);
        M.sigma_img[shift + i][shift + j] = Mi.sigma_img[i][j];
      }
      M.qnum[shift + i] = ModPos64(Mi.qnum[i] * scale, 2 * M.N);
    }
  };
  insert(M1, 0);
  insert(M2, r1);
  return M;
}

// The element of index idx in the mixed radix of the orders, and back.
inline ModuleElt ModuleEltOfIndex(PrimaryQuadModule const &M, size_t idx) {
  ModuleElt x(M.rank());
  for (size_t i = 0; i < M.rank(); i++) {
    size_t o = static_cast<size_t>(M.ord[i]);
    x[i] = static_cast<int64_t>(idx % o);
    idx /= o;
  }
  return x;
}

inline size_t ModuleIndexOfElt(PrimaryQuadModule const &M, ModuleElt const &x) {
  size_t idx = 0;
  for (size_t i = M.rank(); i-- > 0;) {
    idx = idx * static_cast<size_t>(M.ord[i]) + static_cast<size_t>(x[i]);
  }
  return idx;
}

inline ModuleElt ModuleSum(PrimaryQuadModule const &M, ModuleElt const &x,
                           ModuleElt const &y) {
  ModuleElt z(M.rank());
  for (size_t i = 0; i < M.rank(); i++) {
    z[i] = ModPos64(x[i] + y[i], M.ord[i]);
  }
  return z;
}

/*
  A sigma-stable subgroup: its elements as sorted indices, and generators,
  from which the vectors of the glue are lifted.
 */
struct GlueSubgroup {
  std::vector<size_t> elements;
  std::vector<ModuleElt> gens;
};

/*
  The sigma-stable isotropic subgroups of M = M1 (+) M2 meeting neither
  summand, of order dividing max_order (0: any order), the trivial one
  included. They are grown one generator at a time from the candidates,
  the isotropic elements with both projections nonzero: a subgroup is the
  span of the sigma-orbits of its generators, and is kept when all its
  elements are isotropic (then the bilinear form vanishes on it too) and
  every nonzero element has both projections nonzero. Every admissible
  subgroup is reached, since its elements are candidates and it is the
  union of a chain of spans of them. Nothing is returned beyond the bounds.

  With group generators (automorphisms of M commuting with sigma), one
  subgroup per orbit is returned: the orbit of each new subgroup is
  marked as known, and only the representatives are extended, which is
  enough since phi(span(Gamma, x)) = span(phi(Gamma), phi(x)).
 */
inline std::optional<std::vector<GlueSubgroup>>
GlueSubgroups(PrimaryQuadModule const &M, size_t r1, int64_t max_order,
              std::vector<ModuleAut> const &group_gens,
              IsometryGlueBounds const &bounds) {
  if (M.size() > bounds.max_module) {
    return {};
  }
  size_t n_elt = static_cast<size_t>(M.size());
  size_t r = M.rank();
  std::vector<ModuleElt> elts(n_elt);
  std::vector<int64_t> qval(n_elt);
  std::vector<uint8_t> admissible(n_elt);
  for (size_t u = 0; u < n_elt; u++) {
    elts[u] = ModuleEltOfIndex(M, u);
    qval[u] = ModuleQ(M, elts[u]);
    bool nz1 = false, nz2 = false;
    for (size_t i = 0; i < r; i++) {
      if (elts[u][i] != 0) {
        (i < r1 ? nz1 : nz2) = true;
      }
    }
    admissible[u] = (u == 0) || (nz1 && nz2 && qval[u] == 0);
  }
  std::vector<size_t> candidates;
  for (size_t u = 1; u < n_elt; u++) {
    if (admissible[u]) {
      candidates.push_back(u);
    }
  }
  auto apply_sigma = [&](ModuleElt const &x) -> ModuleElt {
    return ModuleCombination(M, x, M.sigma_img);
  };
  // the span of the sigma-orbits of gens, or nothing if not admissible
  auto span = [&](std::vector<ModuleElt> const &gens)
      -> std::optional<std::vector<size_t>> {
    std::vector<ModuleElt> trans;
    for (auto const &g : gens) {
      ModuleElt h = g;
      for (int k = 0; k < 100000; k++) {
        trans.push_back(h);
        h = apply_sigma(h);
        if (h == g) {
          break;
        }
      }
    }
    std::set<size_t> seen{0};
    std::vector<size_t> queue{0};
    while (!queue.empty()) {
      size_t u = queue.back();
      queue.pop_back();
      for (auto const &h : trans) {
        size_t v = ModuleIndexOfElt(M, ModuleSum(M, elts[u], h));
        if (seen.insert(v).second) {
          if (!admissible[v]) {
            return {};
          }
          if (max_order != 0 &&
              static_cast<int64_t>(seen.size()) > max_order) {
            return {};
          }
          queue.push_back(v);
        }
      }
    }
    return std::vector<size_t>(seen.begin(), seen.end());
  };
  auto image = [&](std::vector<size_t> const &elements,
                   ModuleAut const &phi) -> std::vector<size_t> {
    std::vector<size_t> img;
    for (size_t u : elements) {
      img.push_back(ModuleIndexOfElt(M, ModuleCombination(M, elts[u], phi)));
    }
    std::sort(img.begin(), img.end());
    return img;
  };
  std::vector<GlueSubgroup> res{{{0}, {}}};
  std::set<std::vector<size_t>> known{{0}};
  // marks the orbit of a new subgroup, returns false if it was known
  auto insert_orbit = [&](std::vector<size_t> const &elements) -> bool {
    if (!known.insert(elements).second) {
      return false;
    }
    std::vector<std::vector<size_t>> queue{elements};
    while (!queue.empty()) {
      std::vector<size_t> x = queue.back();
      queue.pop_back();
      for (auto const &phi : group_gens) {
        std::vector<size_t> y = image(x, phi);
        if (known.insert(y).second) {
          queue.push_back(y);
        }
      }
    }
    return true;
  };
  for (size_t pos = 0; pos < res.size(); pos++) {
    std::vector<size_t> const elements = res[pos].elements;
    std::vector<ModuleElt> const gens = res[pos].gens;
    for (size_t u : candidates) {
      if (std::binary_search(elements.begin(), elements.end(), u)) {
        continue;
      }
      std::vector<ModuleElt> gens_new = gens;
      gens_new.push_back(elts[u]);
      std::optional<std::vector<size_t>> opt = span(gens_new);
      if (!opt) {
        continue;
      }
      int64_t order = static_cast<int64_t>(opt->size());
      if (max_order != 0 && max_order % order != 0) {
        continue;
      }
      if (insert_orbit(*opt)) {
        res.push_back({*opt, gens_new});
        if (res.size() > bounds.max_subgroup) {
          return {};
        }
      }
    }
  }
  return res;
}

// The floor of the square root of a nonnegative integer.
template <typename Tint> Tint IntegerSquareRoot(Tint const &a) {
  Tint lo(0), hi(1);
  while (hi * hi <= a) {
    hi *= 2;
  }
  // lo^2 <= a < hi^2
  while (hi - lo > 1) {
    Tint mid = (lo + hi) / 2;
    if (mid * mid <= a) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return lo;
}

template <typename T> struct GluedPair {
  MyMatrix<T> gram;
  // acting on the columns, g^T A g = A
  MyMatrix<T> isom;
  T det;
  T glue_order;
};

template <typename T> MyMatrix<T> BlockDiagonal(MyMatrix<T> const &X,
                                                MyMatrix<T> const &Y) {
  int a = X.rows(), b = Y.rows();
  MyMatrix<T> Z = ZeroMatrix<T>(a + b, a + b);
  for (int i = 0; i < a; i++) {
    for (int j = 0; j < a; j++) {
      Z(i, j) = X(i, j);
    }
  }
  for (int i = 0; i < b; i++) {
    for (int j = 0; j < b; j++) {
      Z(a + i, a + j) = Y(i, j);
    }
  }
  return Z;
}

/*
  The generators of the q-primary part of D(L) as vectors of L^# (rows, in
  the coordinates of L), in the order of BuildPrimaryQuadModule.
 */
template <typename T>
std::vector<MyVector<T>> PrimaryModuleGenerators(IsometryComponent<T> const &comp,
                                                 PrimaryQuadModule const &M) {
  std::vector<MyVector<T>> gens;
  int k = comp.gram.rows();
  for (size_t s = 0; s < M.rank(); s++) {
    int i = M.index[s];
    MyVector<T> g(k);
    for (int j = 0; j < k; j++) {
      g(j) = comp.row(i, j) / ScalarFromInt64<T>(M.ord[s]);
    }
    gens.push_back(g);
  }
  return gens;
}

/*
  The automorphism of the q-primary module M of D(L) induced by an isometry
  h of L acting on the rows (h A h^T = A), as the images of the generators.
 */
template <typename T, typename Tint>
ModuleAut ModuleActionOfIsometry(IsometryComponent<T> const &comp,
                                 PrimaryQuadModule const &M,
                                 std::vector<MyVector<T>> const &gens,
                                 MyMatrix<Tint> const &h) {
  MyMatrix<T> h_T = UniversalMatrixConversion<T, Tint>(h);
  ModuleAut phi;
  for (auto const &g : gens) {
    MyVector<T> img = h_T.transpose() * g;
    phi.push_back(ModuleCoordinates(M, comp, img));
  }
  return phi;
}

// phi on the first summand of M1 (+) M2 (shift 0) or the second.
inline ModuleAut ExtendModuleAut(ModuleAut const &phi, size_t r1, size_t r2,
                                 bool first) {
  size_t r = r1 + r2;
  ModuleAut res(r, ModuleElt(r, 0));
  size_t shift = first ? 0 : r1;
  size_t ri = first ? r1 : r2;
  for (size_t i = 0; i < r; i++) {
    if (i >= shift && i < shift + ri) {
      for (size_t j = 0; j < ri; j++) {
        res[i][shift + j] = phi[i - shift][j];
      }
    } else {
      res[i][i] = 1;
    }
  }
  return res;
}

/*
  The distinct classes of the sigma-stable primitive extensions of (L_1,
  sigma_1) (+) (L_2, sigma_2), of determinant det when det is given, as
  canonical pairs.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<GluedPair<T>>
GlueIsometryPairs(MyMatrix<T> const &A1, MyMatrix<T> const &g1,
                  MyMatrix<T> const &A2, MyMatrix<T> const &g2,
                  std::optional<T> const &det, IsometryGlueBounds const &bounds,
                  std::ostream &os) {
  CheckIsometryInput(A1, g1);
  CheckIsometryInput(A2, g2);
  int n1 = A1.rows(), n2 = A2.rows(), n = n1 + n2;
  MyMatrix<T> P1 = g1.transpose(), P2 = g2.transpose();
  IsometryComponent<T> comp1 = BuildIsometryComponent<T, Tint>(
      A1, IdentityMat<T>(n1), P1, IsometryOrder(P1));
  IsometryComponent<T> comp2 = BuildIsometryComponent<T, Tint>(
      A2, IdentityMat<T>(n2), P2, IsometryOrder(P2));
  T det1 = DeterminantMat(A1), det2 = DeterminantMat(A2);
  // the order of the glue: |Gamma|^2 = det1 det2 / det
  std::optional<Tint> target;
  if (det) {
    T sq = det1 * det2 / *det;
    if (!IsInteger(sq)) {
      return {};
    }
    Tint sq_i = UniversalScalarConversion<Tint, T>(sq);
    Tint root = IntegerSquareRoot(sq_i);
    if (root * root != sq_i) {
      return {};
    }
    target = root;
  }
  // the primes dividing both determinants
  Tint g_i = UniversalScalarConversion<Tint, T>(GcdPair(det1, det2));
  std::vector<Tint> primes;
  if (g_i > 1) {
    primes = FactorsInt(g_i);
    std::sort(primes.begin(), primes.end());
    primes.erase(std::unique(primes.begin(), primes.end()), primes.end());
  }
  if (target) {
    // the glue order is supported on these primes
    Tint rest = *target;
    for (auto const &q : primes) {
      while (ResInt(rest, q) == 0) {
        rest /= q;
      }
    }
    if (rest != 1) {
      return {};
    }
  }
  /*
    The glue up to the centralisers C_i of sigma_i in O(L_i): glue groups in
    the same orbit of C_1 x C_2 give isometric pairs. The orbits are taken
    at one prime only, the one with the largest module: the group acts on
    all the primes at once, and moving the part at that prime to its
    representative leaves the other parts arbitrary, which are all kept.
   */
  std::vector<MyMatrix<Tint>> ListGen1 =
      GetPairAutInfo<T, Tint, Tgroup>(A1, std::vector<MyMatrix<T>>{P1}, os)
          .ListGenMat;
  std::vector<MyMatrix<Tint>> ListGen2 =
      GetPairAutInfo<T, Tint, Tgroup>(A2, std::vector<MyMatrix<T>>{P2}, os)
          .ListGenMat;
  size_t i_reduced = 0;
  {
    int64_t best = -1;
    for (size_t iq = 0; iq < primes.size(); iq++) {
      int64_t q = UniversalScalarConversion<int64_t, Tint>(primes[iq]);
      int64_t size = BuildPrimaryQuadModule(comp1, q).size() *
                     BuildPrimaryQuadModule(comp2, q).size();
      if (size > best) {
        best = size;
        i_reduced = iq;
      }
    }
  }
  // per prime, the subgroups with the lifts of their generators
  std::vector<std::vector<std::vector<MyVector<T>>>> ListLifts;
  std::vector<std::vector<int64_t>> ListOrders;
  for (size_t iq = 0; iq < primes.size(); iq++) {
    Tint const &q_i = primes[iq];
    int64_t q = UniversalScalarConversion<int64_t, Tint>(q_i);
    PrimaryQuadModule M1 = BuildPrimaryQuadModule(comp1, q);
    PrimaryQuadModule M2 = BuildPrimaryQuadModule(comp2, q);
    PrimaryQuadModule M = DirectSumModule(M1, M2);
    std::vector<MyVector<T>> gens1 = PrimaryModuleGenerators(comp1, M1);
    std::vector<MyVector<T>> gens2 = PrimaryModuleGenerators(comp2, M2);
    std::vector<ModuleAut> group_gens;
    if (iq == i_reduced) {
      for (auto const &h : ListGen1) {
        group_gens.push_back(ExtendModuleAut(
            ModuleActionOfIsometry(comp1, M1, gens1, h), M1.rank(),
            M2.rank(), true));
      }
      for (auto const &h : ListGen2) {
        group_gens.push_back(ExtendModuleAut(
            ModuleActionOfIsometry(comp2, M2, gens2, h), M1.rank(),
            M2.rank(), false));
      }
    }
    int64_t max_order = 0;
    if (target) {
      Tint qe(1);
      Tint rest = *target;
      while (ResInt(rest, q_i) == 0) {
        rest /= q_i;
        qe *= q_i;
      }
      max_order = UniversalScalarConversion<int64_t, Tint>(qe);
    }
    std::optional<std::vector<GlueSubgroup>> opt =
        GlueSubgroups(M, M1.rank(), max_order, group_gens, bounds);
    if (!opt) {
      std::cerr << "ISOM_GLUE: the glue at q=" << q
                << " exceeds the bounds\n";
      throw TerminalException{1};
    }
    std::vector<std::vector<MyVector<T>>> lifts;
    std::vector<int64_t> orders;
    for (auto const &sg : *opt) {
      int64_t order = static_cast<int64_t>(sg.elements.size());
      if (max_order != 0 && order != max_order) {
        continue;
      }
      std::vector<MyVector<T>> vects;
      for (auto const &x : sg.gens) {
        MyVector<T> v = ZeroVector<T>(n);
        for (size_t i = 0; i < M1.rank(); i++) {
          for (int j = 0; j < n1; j++) {
            v(j) += ScalarFromInt64<T>(x[i]) * gens1[i](j);
          }
        }
        for (size_t i = 0; i < M2.rank(); i++) {
          for (int j = 0; j < n2; j++) {
            v(n1 + j) += ScalarFromInt64<T>(x[M1.rank() + i]) * gens2[i](j);
          }
        }
        vects.push_back(v);
      }
      lifts.push_back(vects);
      orders.push_back(order);
    }
#ifdef DEBUG_ISOMETRY_GLUE
    os << "ISOM_GLUE: q=" << q << " |D1|=" << M1.size() << " |D2|="
       << M2.size() << " subgroups=" << lifts.size() << "\n";
#endif
    ListLifts.push_back(lifts);
    ListOrders.push_back(orders);
  }
  MyMatrix<T> A0 = BlockDiagonal(A1, A2);
  MyMatrix<T> P0 = BlockDiagonal(P1, P2);
  std::vector<GluedPair<T>> res;
  std::set<std::pair<int, std::vector<MyMatrix<T>>>> known;
  // the product over the primes
  std::vector<size_t> choice(primes.size(), 0);
  while (true) {
    std::vector<MyVector<T>> rows;
    T glue_order(1);
    for (size_t iq = 0; iq < primes.size(); iq++) {
      for (auto const &v : ListLifts[iq][choice[iq]]) {
        rows.push_back(v);
      }
      glue_order *= ScalarFromInt64<T>(ListOrders[iq][choice[iq]]);
    }
    bool ok = true;
    for (size_t iq = 0; iq < primes.size(); iq++) {
      ok = ok && !ListLifts[iq].empty();
    }
    if (ok) {
      MyMatrix<T> Gen(n + static_cast<int>(rows.size()), n);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          Gen(i, j) = (i == j) ? T(1) : T(0);
        }
      }
      for (size_t k = 0; k < rows.size(); k++) {
        for (int j = 0; j < n; j++) {
          Gen(n + static_cast<int>(k), j) = rows[k](j);
        }
      }
      // a Z-basis of the overlattice, through integer coordinates
      std::vector<T> ListDen;
      for (int i = 0; i < Gen.rows(); i++) {
        for (int j = 0; j < n; j++) {
          ListDen.push_back(GetDenominator(Gen(i, j)));
        }
      }
      T den = LCMlist(ListDen);
      MyMatrix<Tint> Gen_i =
          UniversalMatrixConversion<Tint, T>(MyMatrix<T>(Gen * den));
      MyMatrix<Tint> B_i = GetZbasis(Gen_i);
      MyMatrix<T> B = UniversalMatrixConversion<T, Tint>(B_i) / den;
      MyMatrix<T> A = B * A0 * B.transpose();
      MyMatrix<T> P = B * P0 * Inverse(B);
      if (!IsIntegralMatrix(A) || !IsIntegralMatrix(P)) {
        std::cerr << "ISOM_GLUE: the glued lattice is not integral or not "
                     "stable\n";
        throw TerminalException{1};
      }
      for (int i = 0; i < n; i++) {
        if (!IsInteger(A(i, i) / T(2))) {
          std::cerr << "ISOM_GLUE: the glued lattice is not even\n";
          throw TerminalException{1};
        }
      }
      CanonicalPairData<T> can =
          CanonicalPair<T, Tint, Tgroup>(A, std::vector<MyMatrix<T>>{P}, os);
      std::pair<int, std::vector<MyMatrix<T>>> key{
          can.method, PairForms(can.Acan, can.ListPcan)};
      if (known.insert(key).second) {
        T det_glued = DeterminantMat(can.Acan);
        res.push_back({can.Acan, can.ListPcan[0].transpose(), det_glued,
                       glue_order});
      }
    }
    size_t iq = 0;
    while (iq < primes.size()) {
      choice[iq]++;
      if (choice[iq] < ListLifts[iq].size()) {
        break;
      }
      choice[iq] = 0;
      iq++;
    }
    if (iq == primes.size()) {
      break;
    }
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_ISOMETRYGLUE_H_
// clang-format on

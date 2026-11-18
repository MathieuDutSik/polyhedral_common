// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ISOMETRYGLUE_H_
#define SRC_GENERA_ISOMETRYGLUE_H_

// clang-format off
#include "ResidueArithmeticLegacy.h"
#include "GeneraArithmetic.h"
#include "IsometryMass.h"
#include "GenusEnumerationGroup.h"
#include <functional>
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
  size_t max_module = size_t(1) << 20;
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
  // the further generators, block diagonal as sigma
  size_t n_extra = M1.extra_sigma_img.size();
  if (M2.extra_sigma_img.size() != n_extra) {
    std::cerr << "ISOM_GLUE: the two modules have different numbers of "
                 "generators\n";
    throw TerminalException{1};
  }
  M.extra_sigma_img.assign(n_extra, std::vector<ModuleElt>(r, ModuleElt(r, 0)));
  auto insert = [&](PrimaryQuadModule const &Mi, size_t shift) -> void {
    int64_t scale = M.N / Mi.N;
    for (size_t i = 0; i < Mi.rank(); i++) {
      for (size_t j = 0; j < Mi.rank(); j++) {
        M.bnum[shift + i][shift + j] = ResMulI64(Mi.bnum[i][j], scale, M.N);
        M.sigma_img[shift + i][shift + j] = Mi.sigma_img[i][j];
        for (size_t e = 0; e < n_extra; e++) {
          M.extra_sigma_img[e][shift + i][shift + j] =
              Mi.extra_sigma_img[e][i][j];
        }
      }
      M.qnum[shift + i] = ResMulI64(Mi.qnum[i], scale, 2 * M.N);
    }
  };
  insert(M1, 0);
  insert(M2, r1);
  return M;
}

inline ModuleElt ModuleSum(PrimaryQuadModule const &M, ModuleElt const &x,
                           ModuleElt const &y) {
  ModuleElt z(M.rank());
  for (size_t i = 0; i < M.rank(); i++) {
    z[i] = ResAddI64(x[i], y[i], M.ord[i]);
  }
  return z;
}

/*
  A sigma-stable subgroup: its elements as sorted indices, generators (the
  subgroup is the span of their sigma-orbits) and the elements of these
  orbits, which generate it as a group and from which the vectors of the
  glue are lifted: lifting the generators alone gives a lattice that is not
  sigma-stable when sigma moves them, which happens when sigma acts on both
  summands by more than +-1.
 */
struct GlueSubgroup {
  std::vector<size_t> elements;
  std::vector<ModuleElt> gens;
  std::vector<ModuleElt> orbit_gens;
};

/*
  The sigma-stable isotropic subgroups of M = M1 (+) M2 meeting neither
  summand, of order dividing max_order (0: any order), the trivial one
  included, inside the elements killed by exponent (a power of q, 0: no
  bound; the glue is killed by the resultant of the minimal polynomials
  of the two sides, see GlueIsometryPairs). That submodule is closed under
  the sums, sigma and the automorphisms, so it is listed alone: on a
  generator of order q^k its coordinate is c q^(k - min(k, v)), c taken
  modulo q^min(k, v), for exponent = q^v. They are grown one generator at a time from the candidates,
  the isotropic elements with both projections nonzero: a subgroup is the
  span of the sigma-orbits of its generators, and is kept when all its
  elements are isotropic (then the bilinear form vanishes on it too) and
  every nonzero element has both projections nonzero. Every admissible
  subgroup is reached, since its elements are candidates and it is the
  union of a chain of spans of them. Nothing is returned beyond the bounds.

  With action generators (automorphisms of M commuting with sigma), all
  the subgroups are returned with the permutations of the list by these
  generators, for the orbits on several primes at once in
  GlueIsometryPairs.

  With group generators (automorphisms of M commuting with sigma), one
  subgroup per orbit is returned: the orbit of each new subgroup is
  marked as known, and only the representatives are extended, which is
  enough since phi(span(Gamma, x)) = span(phi(Gamma), phi(x)).
 */
struct GlueSubgroupList {
  std::vector<GlueSubgroup> subgroups;
  // with action generators: perms[k][i] is the index of the image of
  // subgroups[i] by the k-th generator
  std::vector<std::vector<size_t>> perms;
};

inline std::optional<GlueSubgroupList>
GlueSubgroups(PrimaryQuadModule const &M, size_t r1, int64_t max_order,
              int64_t exponent, std::vector<ModuleAut> const &group_gens,
              std::vector<ModuleAut> const &action_gens,
              IsometryGlueBounds const &bounds) {
  size_t r = M.rank();
  // the orders of the coordinates c in the torsion submodule, and the
  // multipliers from c to the coordinates in M
  std::vector<int64_t> sub_ord(r), mult(r);
  size_t n_elt = 1;
  for (size_t i = 0; i < r; i++) {
    sub_ord[i] = M.ord[i];
    if (exponent != 0 && exponent < M.ord[i]) {
      sub_ord[i] = exponent;
    }
    mult[i] = M.ord[i] / sub_ord[i];
    size_t o = static_cast<size_t>(sub_ord[i]);
    if (n_elt > bounds.max_module / o) {
      return {};
    }
    n_elt *= o;
  }
  auto elt_of_index = [&](size_t idx) -> ModuleElt {
    ModuleElt x(r);
    for (size_t i = 0; i < r; i++) {
      size_t o = static_cast<size_t>(sub_ord[i]);
      x[i] = static_cast<int64_t>(idx % o) * mult[i];
      idx /= o;
    }
    return x;
  };
  auto index_of_elt = [&](ModuleElt const &x) -> size_t {
    size_t idx = 0;
    for (size_t i = r; i-- > 0;) {
      idx = idx * static_cast<size_t>(sub_ord[i]) +
            static_cast<size_t>(x[i] / mult[i]);
    }
    return idx;
  };
  std::vector<ModuleElt> elts(n_elt);
  std::vector<int64_t> qval(n_elt);
  std::vector<uint8_t> admissible(n_elt);
  for (size_t u = 0; u < n_elt; u++) {
    elts[u] = elt_of_index(u);
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
  // the elements of the sigma-orbits of gens
  // the orbits under the group generated by sigma and the further
  // generators (a single cycle of sigma for one isometry)
  auto orbits = [&](std::vector<ModuleElt> const &gens)
      -> std::vector<ModuleElt> {
    std::vector<ModuleElt> trans;
    if (M.extra_sigma_img.empty()) {
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
      return trans;
    }
    std::set<ModuleElt> seen;
    for (auto const &g : gens) {
      if (!seen.insert(g).second) {
        continue;
      }
      std::vector<ModuleElt> queue{g};
      while (!queue.empty()) {
        ModuleElt h = queue.back();
        queue.pop_back();
        trans.push_back(h);
        std::vector<ModuleElt> imgs{apply_sigma(h)};
        for (auto const &img : M.extra_sigma_img) {
          imgs.push_back(ModuleCombination(M, h, img));
        }
        for (auto const &x : imgs) {
          if (seen.insert(x).second) {
            queue.push_back(x);
          }
        }
      }
    }
    return trans;
  };
  // the span of the sigma-orbits of gens, or nothing if not admissible
  auto span = [&](std::vector<ModuleElt> const &gens)
      -> std::optional<std::vector<size_t>> {
    std::vector<ModuleElt> trans = orbits(gens);
    std::set<size_t> seen{0};
    std::vector<size_t> queue{0};
    while (!queue.empty()) {
      size_t u = queue.back();
      queue.pop_back();
      for (auto const &h : trans) {
        size_t v = index_of_elt(ModuleSum(M, elts[u], h));
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
      img.push_back(index_of_elt(ModuleCombination(M, elts[u], phi)));
    }
    std::sort(img.begin(), img.end());
    return img;
  };
  std::vector<GlueSubgroup> res{{{0}, {}, {}}};
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
        res.push_back({*opt, gens_new, orbits(gens_new)});
        if (res.size() > bounds.max_subgroup) {
          return {};
        }
      }
    }
  }
  GlueSubgroupList out;
  if (!action_gens.empty()) {
    std::map<std::vector<size_t>, size_t> index;
    for (size_t i = 0; i < res.size(); i++) {
      index[res[i].elements] = i;
    }
    for (auto const &phi : action_gens) {
      std::vector<size_t> perm;
      for (auto const &sg : res) {
        auto iter = index.find(image(sg.elements, phi));
        if (iter == index.end()) {
          std::cerr << "ISOM_GLUE: the subgroups are not stable under the "
                       "centralizers\n";
          throw TerminalException{1};
        }
        perm.push_back(iter->second);
      }
      out.perms.push_back(perm);
    }
  }
  out.subgroups = std::move(res);
  return out;
}

/*
  Polynomials as the vectors of their coefficients, of increasing degree.
 */
template <typename T> std::vector<T> PolyProductCoeffs(std::vector<T> const &a,
                                                       std::vector<T> const &b) {
  std::vector<T> c(a.size() + b.size() - 1, T(0));
  for (size_t i = 0; i < a.size(); i++) {
    for (size_t j = 0; j < b.size(); j++) {
      c[i + j] += a[i] * b[j];
    }
  }
  return c;
}

// The monic minimal polynomial of the square matrix g.
template <typename T>
std::vector<T> MinimalPolynomialCoeffs(MyMatrix<T> const &g) {
  int n = g.rows();
  auto flat = [&](MyMatrix<T> const &M) -> MyVector<T> {
    MyVector<T> v(n * n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        v(i * n + j) = M(i, j);
      }
    }
    return v;
  };
  std::vector<MyVector<T>> pw{flat(IdentityMat<T>(n))};
  MyMatrix<T> M = IdentityMat<T>(n);
  while (true) {
    M = M * g;
    MyVector<T> w = flat(M);
    std::optional<MyVector<T>> opt =
        SolutionMat(MatrixFromVectorFamily(pw), w);
    if (opt) {
      // g^k = sum_i c_i g^i
      std::vector<T> mu;
      for (int i = 0; i < opt->size(); i++) {
        mu.push_back(-(*opt)(i));
      }
      mu.push_back(T(1));
      return mu;
    }
    pw.push_back(w);
  }
}

// The resultant of a and b, the determinant of their Sylvester matrix.
template <typename T>
T PolyResultant(std::vector<T> const &a, std::vector<T> const &b) {
  int da = static_cast<int>(a.size()) - 1;
  int db = static_cast<int>(b.size()) - 1;
  if (da == 0 && db == 0) {
    return T(1);
  }
  if (da == 0) {
    return MyPow<T>(a[0], db);
  }
  if (db == 0) {
    return MyPow<T>(b[0], da);
  }
  int m = da + db;
  MyMatrix<T> S = ZeroMatrix<T>(m, m);
  for (int i = 0; i < db; i++) {
    for (int k = 0; k <= da; k++) {
      S(i, i + k) = a[da - k];
    }
  }
  for (int i = 0; i < da; i++) {
    for (int k = 0; k <= db; k++) {
      S(db + i, i + k) = b[db - k];
    }
  }
  return DeterminantMat(S);
}

template <typename T> struct GluedPair {
  MyMatrix<T> gram;
  // acting on the columns, g^T A g = A
  MyMatrix<T> isom;
  T det;
  T glue_order;
  // the further generators of a group, acting on the columns (empty for a
  // single isometry)
  std::vector<MyMatrix<T>> extra_isom;
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
  sigma_1) (+) (L_2, sigma_2), of determinant det when det is given.

  An isometry between two of them maps the sublattices ker mu_i(sigma) onto
  each other, so restricts to an element of C_1 x C_2, the centralizers of
  the sigma_i in O(L_i): the classes are the C_1 x C_2-orbits of the glue
  groups. At one prime, GlueSubgroups gives one subgroup per orbit. At
  several primes, the group acts on the lists of all the subgroups of each
  prime: one subgroup per orbit at the first prime, then per orbit of its
  stabilizer at the second, and so on. Either way the glued pairs are
  distinct classes, returned in an LLL-reduced basis.

  For a group G, gens1 and gens2 are the actions of the same generators of
  G on L_1 and L_2: the glue is G-stable (spans of G-orbits), killed by
  the gcd over the generators of the resultants, and C_i is the
  centralizer of G in O(L_i).
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<GluedPair<T>>
GlueGroupPairs(MyMatrix<T> const &A1, std::vector<MyMatrix<T>> const &gens1,
               MyMatrix<T> const &A2, std::vector<MyMatrix<T>> const &gens2,
               std::optional<T> const &det,
               std::optional<Tint> const &glue_bound,
               IsometryGlueBounds const &bounds, std::ostream &os) {
  if (gens1.empty() || gens1.size() != gens2.size()) {
    std::cerr << "ISOM_GLUE: the two pairs need the same number of "
                 "generators\n";
    throw TerminalException{1};
  }
  size_t n_gen = gens1.size();
  std::vector<MyMatrix<T>> ListP1, ListP2;
  for (size_t k = 0; k < n_gen; k++) {
    CheckIsometryInput(A1, gens1[k]);
    CheckIsometryInput(A2, gens2[k]);
    ListP1.push_back(gens1[k].transpose());
    ListP2.push_back(gens2[k].transpose());
  }
  int n1 = A1.rows(), n2 = A2.rows(), n = n1 + n2;
  MyMatrix<T> const &P1 = ListP1[0];
  MyMatrix<T> const &P2 = ListP2[0];
  IsometryComponent<T> comp1 = BuildIsometryComponent<T, Tint>(
      A1, IdentityMat<T>(n1), P1, IsometryOrder(P1));
  IsometryComponent<T> comp2 = BuildIsometryComponent<T, Tint>(
      A2, IdentityMat<T>(n2), P2, IsometryOrder(P2));
  // the further generators (the basis is that of L_i: sigma_k = P_k)
  for (size_t k = 1; k < n_gen; k++) {
    comp1.extra_sigma.push_back(ListP1[k]);
    comp2.extra_sigma.push_back(ListP2[k]);
  }
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
  // the glue is killed by g_1 and g_2 acting on it, so by the resultant
  // of their minimal polynomials (u mu_1 + v mu_2 = Res, u, v in Z[x]),
  // when they are coprime: only its primes, and the elements killed by it.
  // For a group, by the gcd of these over the generators (those with a
  // non zero resultant).
  T res_mu(0);
  for (size_t k = 0; k < n_gen; k++) {
    T r = T_abs(PolyResultant(MinimalPolynomialCoeffs(ListP1[k]),
                              MinimalPolynomialCoeffs(ListP2[k])));
    if (r != 0) {
      res_mu = (res_mu == 0) ? r : GcdPair(res_mu, r);
    }
  }
  // a multiple of the order of the glue known to the caller (e.g. by the
  // characters of an abelian group, IsometryGlueFactors.h) bounds its
  // exponent as well
  if (glue_bound && *glue_bound != 0) {
    T b = UniversalScalarConversion<T, Tint>(*glue_bound);
    res_mu = (res_mu == 0) ? b : GcdPair(res_mu, b);
  }
  // the primes dividing both determinants (and the resultant)
  T g_T = GcdPair(det1, det2);
  if (res_mu != 0) {
    g_T = GcdPair(g_T, res_mu);
  }
  Tint g_i = UniversalScalarConversion<Tint, T>(g_T);
  std::vector<Tint> primes;
  if (g_i > 1) {
    primes = FactorsInt(g_i);
    std::sort(primes.begin(), primes.end());
    primes.erase(std::unique(primes.begin(), primes.end()), primes.end());
  }
  // the q-part of the resultant, 0 for no bound
  auto exponent_at = [&](Tint const &q_i) -> int64_t {
    if (res_mu == 0) {
      return 0;
    }
    Tint rest = UniversalScalarConversion<Tint, T>(res_mu);
    Tint qe(1);
    while (ResInt(rest, q_i) == 0) {
      rest /= q_i;
      qe *= q_i;
    }
    return UniversalScalarConversion<int64_t, Tint>(qe);
  };
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
    The glue up to the centralizers C_i of sigma_i in O(L_i): glue groups in
    the same orbit of C_1 x C_2 give isometric pairs. The orbits are taken
    at one prime only, the one with the largest module: the group acts on
    all the primes at once, and moving the part at that prime to its
    representative leaves the other parts arbitrary, which are all kept.
   */
  PairMethods methods = DefaultPairMethods();
  std::vector<MyMatrix<Tint>> ListGen1 =
      PairCentralizer<T, Tint, Tgroup>(A1, ListP1, methods.automorphism, os)
          .ListGen;
  std::vector<MyMatrix<Tint>> ListGen2 =
      PairCentralizer<T, Tint, Tgroup>(A2, ListP2, methods.automorphism, os)
          .ListGen;
  size_t i_reduced = 0;
  {
    size_t best = 0;
    for (size_t iq = 0; iq < primes.size(); iq++) {
      int64_t q = UniversalScalarConversion<int64_t, Tint>(primes[iq]);
      size_t size = SaturatedSizeProduct(BuildPrimaryQuadModule(comp1, q).size(),
                                         BuildPrimaryQuadModule(comp2, q).size());
      if (iq == 0 || size > best) {
        best = size;
        i_reduced = iq;
      }
    }
  }
  // the orbits are taken on all the primes at once when there are several
  bool exact_orbits = primes.size() > 1;
  // per prime, the subgroups with the lifts of their generators, and the
  // permutations of them by the generators of C_1 x C_2 (exact_orbits)
  std::vector<std::vector<std::vector<MyVector<T>>>> ListLifts;
  std::vector<std::vector<int64_t>> ListOrders;
  std::vector<std::vector<std::vector<size_t>>> ListPerms;
  for (size_t iq = 0; iq < primes.size(); iq++) {
    Tint const &q_i = primes[iq];
    int64_t q = UniversalScalarConversion<int64_t, Tint>(q_i);
    PrimaryQuadModule M1 = BuildPrimaryQuadModule(comp1, q);
    PrimaryQuadModule M2 = BuildPrimaryQuadModule(comp2, q);
    PrimaryQuadModule M = DirectSumModule(M1, M2);
    std::vector<MyVector<T>> gens1 = PrimaryModuleGenerators(comp1, M1);
    std::vector<MyVector<T>> gens2 = PrimaryModuleGenerators(comp2, M2);
    std::vector<ModuleAut> group_gens;
    if (iq == i_reduced || exact_orbits) {
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
    std::optional<GlueSubgroupList> opt = GlueSubgroups(
        M, M1.rank(), max_order, exponent_at(q_i),
        exact_orbits ? std::vector<ModuleAut>{} : group_gens,
        exact_orbits ? group_gens : std::vector<ModuleAut>{}, bounds);
    if (!opt) {
      std::cerr << "ISOM_GLUE: the glue at q=" << q
                << " exceeds the bounds\n";
      throw TerminalException{1};
    }
    std::vector<std::vector<MyVector<T>>> lifts;
    std::vector<int64_t> orders;
    // the kept subgroups (of the order asked for) are a union of orbits
    std::vector<size_t> new_index(opt->subgroups.size(), 0);
    for (size_t isg = 0; isg < opt->subgroups.size(); isg++) {
      GlueSubgroup const &sg = opt->subgroups[isg];
      int64_t order = static_cast<int64_t>(sg.elements.size());
      if (max_order != 0 && order != max_order) {
        continue;
      }
      new_index[isg] = lifts.size();
      std::vector<MyVector<T>> vects;
      for (auto const &x : sg.orbit_gens) {
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
    std::vector<std::vector<size_t>> perms;
    for (auto const &perm : opt->perms) {
      std::vector<size_t> p(lifts.size());
      for (size_t isg = 0; isg < opt->subgroups.size(); isg++) {
        int64_t order = static_cast<int64_t>(opt->subgroups[isg].elements.size());
        if (max_order == 0 || order == max_order) {
          p[new_index[isg]] = new_index[perm[isg]];
        }
      }
      perms.push_back(p);
    }
    ListLifts.push_back(lifts);
    ListOrders.push_back(orders);
    ListPerms.push_back(perms);
  }
  MyMatrix<T> A0 = BlockDiagonal(A1, A2);
  std::vector<MyMatrix<T>> ListP0;
  for (size_t k = 0; k < n_gen; k++) {
    ListP0.push_back(BlockDiagonal(ListP1[k], ListP2[k]));
  }
  std::vector<GluedPair<T>> res;
  // no subgroup at some prime (none of the order the determinant asks
  // for): no glued lattice, before the product reads ListLifts[iq][0]
  for (auto const &lifts : ListLifts) {
    if (lifts.empty()) {
      return res;
    }
  }
  auto build = [&](std::vector<size_t> const &choice) -> void {
    std::vector<MyVector<T>> rows;
    T glue_order(1);
    for (size_t iq = 0; iq < primes.size(); iq++) {
      for (auto const &v : ListLifts[iq][choice[iq]]) {
        rows.push_back(v);
      }
      glue_order *= ScalarFromInt64<T>(ListOrders[iq][choice[iq]]);
    }
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
    MyMatrix<T> Binv = Inverse(B);
    std::vector<MyMatrix<T>> ListP;
    for (auto const &P0 : ListP0) {
      ListP.push_back(B * P0 * Binv);
    }
    bool stable = IsIntegralMatrix(A);
    for (auto const &P : ListP) {
      stable = stable && IsIntegralMatrix(P);
    }
    if (!stable) {
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
    MyMatrix<T> B_T =
        UniversalMatrixConversion<T, Tint>(LLLreducedBasis<T, Tint>(A, os).Pmat);
    MyMatrix<T> Ared = B_T * A * B_T.transpose();
    MyMatrix<T> B_Tinv = Inverse(B_T);
    std::vector<MyMatrix<T>> extra;
    for (size_t k = 1; k < n_gen; k++) {
      extra.push_back(MyMatrix<T>(B_T * ListP[k] * B_Tinv).transpose());
    }
    MyMatrix<T> Pred = B_T * ListP[0] * B_Tinv;
    res.push_back(
        {Ared, Pred.transpose(), DeterminantMat(Ared), glue_order, extra});
  };
  if (!exact_orbits) {
    // at most one prime, its subgroups up to C_1 x C_2 already
    std::vector<size_t> choice(primes.size(), 0);
    if (primes.empty()) {
      build(choice);
    } else {
      for (size_t i = 0; i < ListLifts[0].size(); i++) {
        choice[0] = i;
        build(choice);
      }
    }
    return res;
  }
  // C_1 x C_2 acting on the disjoint union of the lists of the primes
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<size_t> offset;
  size_t n_point = 0;
  for (auto const &lifts : ListLifts) {
    offset.push_back(n_point);
    n_point += lifts.size();
  }
  std::vector<Telt> LGen;
  size_t n_perm_gen = ListPerms[0].size();
  for (size_t k = 0; k < n_perm_gen; k++) {
    std::vector<Tidx> v(n_point);
    for (size_t iq = 0; iq < primes.size(); iq++) {
      for (size_t i = 0; i < ListLifts[iq].size(); i++) {
        v[offset[iq] + i] = static_cast<Tidx>(offset[iq] + ListPerms[iq][k][i]);
      }
    }
    LGen.push_back(Telt(v));
  }
  Tgroup G(LGen, Telt(static_cast<Tidx>(n_point)));
  std::vector<size_t> choice(primes.size(), 0);
  std::function<void(size_t, Tgroup const &)> rec =
      [&](size_t iq, Tgroup const &H) -> void {
    if (iq == primes.size()) {
      build(choice);
      return;
    }
    std::vector<Telt> hgens = H.GeneratorsOfGroup();
    size_t sz = ListLifts[iq].size();
    std::vector<uint8_t> seen(sz, 0);
    for (size_t i = 0; i < sz; i++) {
      if (seen[i]) {
        continue;
      }
      // the orbit of i under H
      seen[i] = 1;
      std::vector<size_t> queue{i};
      while (!queue.empty()) {
        size_t u = queue.back();
        queue.pop_back();
        for (auto const &h : hgens) {
          size_t w = static_cast<size_t>(h.at(static_cast<Tidx>(offset[iq] + u))) -
                     offset[iq];
          if (!seen[w]) {
            seen[w] = 1;
            queue.push_back(w);
          }
        }
      }
      choice[iq] = i;
      if (iq + 1 < primes.size()) {
        rec(iq + 1, H.Stabilizer_OnPoints(static_cast<Tidx>(offset[iq] + i)));
      } else {
        rec(iq + 1, H);
      }
    }
  };
  rec(0, G);
  return res;
}

// The case of a single isometry on each side.
template <typename T, typename Tint, typename Tgroup>
std::vector<GluedPair<T>>
GlueIsometryPairs(MyMatrix<T> const &A1, MyMatrix<T> const &g1,
                  MyMatrix<T> const &A2, MyMatrix<T> const &g2,
                  std::optional<T> const &det, IsometryGlueBounds const &bounds,
                  std::ostream &os) {
  return GlueGroupPairs<T, Tint, Tgroup>(A1, std::vector<MyMatrix<T>>{g1}, A2,
                                         std::vector<MyMatrix<T>>{g2}, det,
                                         std::optional<Tint>{}, bounds, os);
}

// clang-format off
#endif  // SRC_GENERA_ISOMETRYGLUE_H_
// clang-format on

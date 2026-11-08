// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ISOMETRYMASS_H_
#define SRC_GENERA_ISOMETRYMASS_H_

// clang-format off
#include "LocalDensity.h"
#include "GenusMass.h"
#include "CyclotomicHermitianMass.h"
#include "MultiplicityLocalFactor.h"
#include "WeightMatrix.h"
#include "GRAPH_GraphicalBasic.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "factorizations.h"
#include <algorithm>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ISOMETRY_MASS
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ISOMETRY_MASS
#endif

/*
  The mass of the genus of a lattice with an isometry of finite order,

      mass(L, sigma) = sum_{(L', sigma)} 1 / |C_{O(L')}(sigma)|,

  the sum over the isometry classes of pairs in the genus of the pair (the
  sigma-stable lattices of V = L (x) Q locally isometric to L as pairs).

  Write C for the centraliser of sigma in O(V), K the stabiliser of L in
  C(A_f), L_d = L cap ker Phi_d(sigma) for the d dividing the order m, and
  L_0 = (+)_d L_d, of finite index in L. An element of C preserves every
  ker Phi_d(sigma), so K is contained in the stabiliser K_0 = prod_d K_d of
  L_0, and since C = prod_d C_d (orthogonal groups for d = 1, 2, unitary
  groups over Q(zeta_d) for d >= 3):

      mass(L, sigma) = prod_d mass(L_d, sigma_d) * [K_0 : K].

  The index is local: [K_0 : K] = prod_q |G_q . H_q|, the orbit of the glue
  H_q = (L / L_0) (x) Z_q in D(L_0)_q = (+)_d D(L_d)_q under the image G_q
  of K_{0,q}. H is a subgroup of the discriminant group and K_q is its
  stabiliser. Only the primes q dividing [L : L_0] contribute, and they
  divide m.

  The image G_q is the product over d of the images of the local groups
  of the L_d in the isometries of the discriminant form D(L_d)_q that
  commute with sigma. For d = 1, 2 (sigma = +-1 on L_d) O(L_{d,q}) maps
  onto O(D(L_d)_q) (Nikulin, Integral symmetric bilinear forms and some of
  their applications, Theorem 1.9.5), so the full group is the image. For
  d >= 3 the code uses all the isometries commuting with sigma; that they
  are the image of the local unitary group is not proved here, and the
  result says when it is relied on.

  The masses of the orthogonal L_d come from GenusMass. The genus of the
  pair (L_d, sigma_d), d >= 3, is the genus of L_d as a Hermitian lattice
  over Z[zeta_d], for the form h with b = Tr h. Its mass is computed
  exactly (CyclotomicHermitianMass.h); a mass supplied by the caller takes
  precedence.

  Conventions: the Gram matrix A of L, the isometry g acting on the column
  vectors, g^T A g = A. Internally the row vectors are acted on by P =
  g^T, as in the rest of the code.
 */

template <typename T>
MyMatrix<T> EvaluatePolynomialMatrix(std::vector<int64_t> const &coeff,
                                     MyMatrix<T> const &P) {
  int n = P.rows();
  MyMatrix<T> res = ZeroMatrix<T>(n, n);
  for (size_t k = coeff.size(); k-- > 0;) {
    res = res * P;
    for (int i = 0; i < n; i++) {
      res(i, i) += ScalarFromInt64<T>(coeff[k]);
    }
  }
  return res;
}

inline int64_t ModPos64(int64_t a, int64_t m) {
  a %= m;
  return (a < 0) ? a + m : a;
}

inline int64_t InverseMod64(int64_t u, int64_t m) {
  int64_t a = ModPos64(u, m), b = m, x0 = 1, x1 = 0;
  while (b != 0) {
    int64_t q = a / b;
    int64_t t = a - q * b;
    a = b;
    b = t;
    t = x0 - q * x1;
    x0 = x1;
    x1 = t;
  }
  if (a != 1) {
    std::cerr << "ISOM_MASS: " << u << " is not invertible modulo " << m << "\n";
    throw TerminalException{1};
  }
  return ModPos64(x0, m);
}

/*
  One component L_d: a basis (rows, coordinates in L), its Gram matrix, the
  action of sigma (sigma(basis_i) = sum_j sigma(i, j) basis_j), and the
  Smith form of the Gram matrix, which coordinatises D(L_d): with ROW G COL
  = diag(s), x in L_d^# (coordinates in the basis) is z ROW with z_i in
  (1/s_i) Z, and its residues are z_i s_i mod s_i.
 */
template <typename T> struct IsometryComponent {
  // d = 0: a component of a non-abelian group on which it does not act
  // through a cyclic quotient
  int d;
  MyMatrix<T> basis;
  MyMatrix<T> gram;
  MyMatrix<T> sigma;
  // the actions of the further generators of the group (non-abelian case)
  std::vector<MyMatrix<T>> extra_sigma;
  // End_G = Q: the lattice X = L cap (l (x) U) of the multiplicity space,
  // whose genus gives the mass up to the local factors at p | |G|
  std::optional<MyMatrix<T>> multiplicity_gram;
  // the product over p | |G| of [O(X_p) : Stab(L_p)], when computed
  std::optional<T> multiplicity_factor;
  MyMatrix<T> row;
  MyMatrix<T> row_inv;
  std::vector<int64_t> elem;
};

template <typename T>
std::vector<int64_t> DiscriminantResidues(IsometryComponent<T> const &comp,
                                          MyVector<T> const &x) {
  MyVector<T> z = comp.row_inv.transpose() * x;
  std::vector<int64_t> res;
  for (size_t i = 0; i < comp.elem.size(); i++) {
    T val = z(i) * ScalarFromInt64<T>(comp.elem[i]);
    if (!IsInteger(val)) {
      std::cerr << "ISOM_MASS: a vector is not in the dual lattice\n";
      throw TerminalException{1};
    }
    int64_t r = UniversalScalarConversion<int64_t, T>(val);
    res.push_back(ModPos64(r, comp.elem[i]));
  }
  return res;
}

/*
  The q-primary part of D(L_d): generators g_i (one per elementary divisor
  divisible by q) of orders ord_i, N the exponent, b(g_i, g_j) = bnum_ij / N
  mod 1, q(g_i) = qnum_i / N mod 2, and the images of the generators by
  sigma in coordinates.
 */
struct PrimaryQuadModule {
  int64_t q = 0;
  int64_t N = 1;
  int d = 0;
  std::vector<int> index;
  std::vector<int64_t> ord;
  std::vector<int64_t> tinv;
  std::vector<std::vector<int64_t>> bnum;
  std::vector<int64_t> qnum;
  std::vector<std::vector<int64_t>> sigma_img;
  // the images by the further generators
  std::vector<std::vector<std::vector<int64_t>>> extra_sigma_img;
  int64_t size() const {
    int64_t s = 1;
    for (int64_t o : ord) {
      s *= o;
    }
    return s;
  }
  size_t rank() const { return ord.size(); }
};

using ModuleElt = std::vector<int64_t>;

inline int64_t ModuleB(PrimaryQuadModule const &M, ModuleElt const &x,
                       ModuleElt const &y) {
  int64_t val = 0;
  for (size_t i = 0; i < M.rank(); i++) {
    if (x[i] == 0) {
      continue;
    }
    for (size_t j = 0; j < M.rank(); j++) {
      int64_t t = ModPos64(x[i] * y[j], M.N);
      val = ModPos64(val + ModPos64(t * M.bnum[i][j], M.N), M.N);
    }
  }
  return val;
}

inline int64_t ModuleQ(PrimaryQuadModule const &M, ModuleElt const &x) {
  int64_t N2 = 2 * M.N;
  int64_t val = 0;
  for (size_t i = 0; i < M.rank(); i++) {
    int64_t t = ModPos64(x[i] * x[i], N2);
    val = ModPos64(val + ModPos64(t * M.qnum[i], N2), N2);
    for (size_t j = i + 1; j < M.rank(); j++) {
      int64_t u = ModPos64(2 * ModPos64(x[i] * x[j], N2), N2);
      val = ModPos64(val + ModPos64(u * M.bnum[i][j], N2), N2);
    }
  }
  return val;
}

// sum_i a_i y_i, the y_i given in coordinates.
inline ModuleElt ModuleCombination(PrimaryQuadModule const &M,
                                   ModuleElt const &a,
                                   std::vector<ModuleElt> const &y) {
  ModuleElt res(M.rank(), 0);
  for (size_t i = 0; i < M.rank(); i++) {
    if (a[i] == 0) {
      continue;
    }
    for (size_t j = 0; j < M.rank(); j++) {
      res[j] = ModPos64(res[j] + ModPos64(a[i] * y[i][j], M.ord[j]), M.ord[j]);
    }
  }
  return res;
}

template <typename T>
ModuleElt ModuleCoordinates(PrimaryQuadModule const &M,
                            IsometryComponent<T> const &comp,
                            MyVector<T> const &x) {
  std::vector<int64_t> res = DiscriminantResidues(comp, x);
  ModuleElt a(M.rank());
  for (size_t i = 0; i < M.rank(); i++) {
    int64_t r = ModPos64(res[M.index[i]], M.ord[i]);
    a[i] = ModPos64(r * M.tinv[i], M.ord[i]);
  }
  return a;
}

template <typename T>
PrimaryQuadModule BuildPrimaryQuadModule(IsometryComponent<T> const &comp,
                                         int64_t q) {
  PrimaryQuadModule M;
  M.q = q;
  M.d = comp.d;
  std::vector<MyVector<T>> gens;
  int k = comp.gram.rows();
  for (size_t i = 0; i < comp.elem.size(); i++) {
    int64_t s = comp.elem[i];
    int64_t qe = 1;
    while (s % (qe * q) == 0) {
      qe *= q;
    }
    if (qe == 1) {
      continue;
    }
    M.index.push_back(static_cast<int>(i));
    M.ord.push_back(qe);
    M.tinv.push_back(InverseMod64(s / qe, qe));
    M.N = std::max(M.N, qe);
    MyVector<T> g(k);
    for (int j = 0; j < k; j++) {
      g(j) = comp.row(i, j) / ScalarFromInt64<T>(qe);
    }
    gens.push_back(g);
  }
  size_t r = gens.size();
  M.bnum.assign(r, std::vector<int64_t>(r, 0));
  M.qnum.assign(r, 0);
  T N_T = ScalarFromInt64<T>(M.N);
  for (size_t i = 0; i < r; i++) {
    for (size_t j = 0; j < r; j++) {
      T val = N_T * gens[i].dot(comp.gram * gens[j]);
      if (!IsInteger(val)) {
        std::cerr << "ISOM_MASS: the discriminant form is not of exponent N\n";
        throw TerminalException{1};
      }
      int64_t v = UniversalScalarConversion<int64_t, T>(val);
      M.bnum[i][j] = ModPos64(v, M.N);
      if (i == j) {
        M.qnum[i] = ModPos64(v, 2 * M.N);
      }
    }
  }
  for (size_t i = 0; i < r; i++) {
    MyVector<T> img = comp.sigma.transpose() * gens[i];
    M.sigma_img.push_back(ModuleCoordinates(M, comp, img));
  }
  for (auto const &sig : comp.extra_sigma) {
    std::vector<std::vector<int64_t>> imgs;
    for (size_t i = 0; i < r; i++) {
      MyVector<T> img = sig.transpose() * gens[i];
      imgs.push_back(ModuleCoordinates(M, comp, img));
    }
    M.extra_sigma_img.push_back(imgs);
  }
  return M;
}

// The group acts by +-1 (d = 1, 2), where Nikulin's surjectivity holds and
// every isometry of the module commutes with it.
inline bool IsScalarAction(PrimaryQuadModule const &M) {
  return M.d >= 1 && M.d <= 2;
}

/*
  All the isometries of the module commuting with sigma, as the images of
  the generators, by backtracking: the image of g_i has an order dividing
  ord_i, the same q value, and the same b with the images already chosen.
  Such a map is a homomorphism preserving the nondegenerate form, hence an
  automorphism. Nothing is returned when the module or the group exceeds
  the bounds.
 */
using ModuleAut = std::vector<ModuleElt>;

inline std::vector<ModuleElt> ListModuleElements(PrimaryQuadModule const &M) {
  std::vector<ModuleElt> res;
  ModuleElt x(M.rank(), 0);
  while (true) {
    res.push_back(x);
    size_t i = 0;
    while (i < M.rank()) {
      x[i]++;
      if (x[i] < M.ord[i]) {
        break;
      }
      x[i] = 0;
      i++;
    }
    if (i == M.rank()) {
      break;
    }
  }
  return res;
}

inline ModuleElt ApplyModuleAut(PrimaryQuadModule const &M, ModuleAut const &phi,
                                ModuleElt const &x) {
  return ModuleCombination(M, x, phi);
}

inline std::optional<std::vector<ModuleAut>>
EquivariantModuleIsometries(PrimaryQuadModule const &M, int64_t max_module,
                            size_t max_group) {
  if (M.size() > max_module) {
    return {};
  }
  std::vector<ModuleElt> elts = ListModuleElements(M);
  std::vector<int64_t> qval(elts.size());
  for (size_t u = 0; u < elts.size(); u++) {
    qval[u] = ModuleQ(M, elts[u]);
  }
  auto order_divides = [&](ModuleElt const &x, int64_t o) -> bool {
    for (size_t j = 0; j < M.rank(); j++) {
      if (ModPos64(o * x[j], M.ord[j]) != 0) {
        return false;
      }
    }
    return true;
  };
  size_t r = M.rank();
  std::vector<std::vector<size_t>> cand(r);
  for (size_t i = 0; i < r; i++) {
    for (size_t u = 0; u < elts.size(); u++) {
      if (qval[u] == M.qnum[i] && order_divides(elts[u], M.ord[i])) {
        cand[i].push_back(u);
      }
    }
  }
  std::vector<ModuleAut> res;
  ModuleAut phi(r);
  bool overflow = false;
  std::function<void(size_t)> rec = [&](size_t i) -> void {
    if (overflow) {
      return;
    }
    if (i == r) {
      if (!IsScalarAction(M)) {
        for (size_t j = 0; j < r; j++) {
          ModuleElt lhs = ApplyModuleAut(M, phi, M.sigma_img[j]);
          ModuleElt rhs = ModuleCombination(M, phi[j], M.sigma_img);
          if (lhs != rhs) {
            return;
          }
          for (auto const &img : M.extra_sigma_img) {
            if (ApplyModuleAut(M, phi, img[j]) !=
                ModuleCombination(M, phi[j], img)) {
              return;
            }
          }
        }
      }
      res.push_back(phi);
      if (res.size() > max_group) {
        overflow = true;
      }
      return;
    }
    for (size_t u : cand[i]) {
      ModuleElt const &y = elts[u];
      bool ok = true;
      for (size_t j = 0; j < i; j++) {
        if (ModuleB(M, y, phi[j]) != M.bnum[i][j]) {
          ok = false;
          break;
        }
      }
      if (!ok) {
        continue;
      }
      phi[i] = y;
      rec(i + 1);
    }
  };
  rec(0);
  if (overflow) {
    return {};
  }
  return res;
}

// A small generating set of a group given by its list of elements.
inline std::vector<ModuleAut>
GeneratingSubset(PrimaryQuadModule const &M, std::vector<ModuleAut> const &grp) {
  auto compose = [&](ModuleAut const &a, ModuleAut const &b) -> ModuleAut {
    // (a o b)(g_i) = a(b(g_i))
    ModuleAut c(M.rank());
    for (size_t i = 0; i < M.rank(); i++) {
      c[i] = ApplyModuleAut(M, a, b[i]);
    }
    return c;
  };
  std::vector<ModuleAut> gens;
  std::set<ModuleAut> closure;
  ModuleAut id(M.rank(), ModuleElt(M.rank(), 0));
  for (size_t i = 0; i < M.rank(); i++) {
    id[i][i] = 1;
  }
  closure.insert(id);
  for (auto const &g : grp) {
    if (closure.count(g) > 0) {
      continue;
    }
    gens.push_back(g);
    std::vector<ModuleAut> work(closure.begin(), closure.end());
    while (!work.empty()) {
      ModuleAut x = work.back();
      work.pop_back();
      for (auto const &h : gens) {
        ModuleAut y = compose(h, x);
        if (closure.insert(y).second) {
          work.push_back(y);
        }
      }
    }
    if (closure.size() == grp.size()) {
      break;
    }
  }
  return gens;
}

/*
  The orbit of the glue at one prime: the modules of the components on
  which the glue has a nonzero projection, the glue generators in each of
  them, and the groups. The subgroups are stored as the sorted list of the
  codes of their elements (mixed radix over all the coordinates).
 */
struct GlueOrbitResult {
  int64_t q;
  int64_t glue_order;
  std::vector<int> list_d;
  std::vector<int64_t> list_module_order;
  std::vector<mpz_class> list_group_order;
  mpz_class orbit_size;
  // "listing" or "graph"
  std::string method;
};

// The elements of the glue, as vectors of all the coordinates.
inline std::vector<ModuleElt>
GlueElements(std::vector<PrimaryQuadModule> const &mods,
             std::vector<std::vector<ModuleElt>> const &glue_gens) {
  size_t t = mods.size();
  std::vector<int64_t> radix;
  std::vector<size_t> offset;
  for (auto const &M : mods) {
    offset.push_back(radix.size());
    for (int64_t o : M.ord) {
      radix.push_back(o);
    }
  }
  size_t R = radix.size();
  std::set<ModuleElt> H;
  H.insert(ModuleElt(R, 0));
  std::vector<ModuleElt> work{ModuleElt(R, 0)};
  while (!work.empty()) {
    ModuleElt x = work.back();
    work.pop_back();
    for (auto const &g : glue_gens) {
      ModuleElt y(R);
      for (size_t s = 0; s < t; s++) {
        for (size_t i = 0; i < mods[s].rank(); i++) {
          size_t pos = offset[s] + i;
          y[pos] = ModPos64(x[pos] + g[s][i], radix[pos]);
        }
      }
      if (H.insert(y).second) {
        work.push_back(y);
      }
    }
  }
  return std::vector<ModuleElt>(H.begin(), H.end());
}

/*
  The orbit length of the glue without listing the groups. For a module M
  with a nondegenerate form, a bijection of M preserving b(x, y) for all x,
  y is additive (phi(x + y) - phi(x) - phi(y) is orthogonal to everything),
  so the isometries of M commuting with sigma are the automorphisms of the
  complete graph on the elements of M with the vertices coloured by q(x)
  and the edges by b(x, y) and by whether y = sigma(x). The glue H is added
  as one vertex per element, joined to its components: the automorphisms of
  that graph are the stabilizer of H in G = prod_s G_s (the module of each
  vertex is in its colour, so the modules are not exchanged, and an element
  of H is determined by its components). The orbit length is |G| / |Stab_G
  H|, the orders being those of the graph automorphism groups.
 */
inline std::optional<GlueOrbitResult>
GlueOrbitGraph(int64_t q, std::vector<PrimaryQuadModule> const &mods,
               std::vector<std::vector<ModuleElt>> const &glue_gens,
               size_t max_vertices, std::ostream &os) {
  using Tgr = GraphListAdj;
  using Tidx = uint32_t;
  using Tidx_value = int32_t;
  GlueOrbitResult res;
  res.q = q;
  res.method = "graph";
  size_t t = mods.size();
  std::vector<ModuleElt> H = GlueElements(mods, glue_gens);
  res.glue_order = static_cast<int64_t>(H.size());
  std::vector<std::vector<ModuleElt>> elts(t);
  size_t n_vert = H.size();
  int64_t BIG = 1;
  for (size_t s = 0; s < t; s++) {
    if (mods[s].size() > static_cast<int64_t>(max_vertices)) {
      return {};
    }
    elts[s] = ListModuleElements(mods[s]);
    n_vert += elts[s].size();
    BIG = std::max(BIG, 4 * mods[s].N + 4);
  }
  if (n_vert > max_vertices) {
    return {};
  }
  // the images by the generators of the group, one edge colour bit each
  std::vector<std::vector<std::vector<ModuleElt>>> sigma_img(t);
  size_t n_bit = 1;
  for (size_t s = 0; s < t; s++) {
    std::vector<std::vector<std::vector<int64_t>>> imgs{mods[s].sigma_img};
    imgs.insert(imgs.end(), mods[s].extra_sigma_img.begin(),
                mods[s].extra_sigma_img.end());
    n_bit = std::max(n_bit, imgs.size());
    for (auto const &img : imgs) {
      std::vector<ModuleElt> col;
      for (auto const &x : elts[s]) {
        col.push_back(ModuleCombination(mods[s], x, img));
      }
      sigma_img[s].push_back(col);
    }
  }
  // the vertices: (module s, index) for the modules, then the glue
  std::vector<std::pair<int, size_t>> vert;
  for (size_t s = 0; s < t; s++) {
    for (size_t i = 0; i < elts[s].size(); i++) {
      vert.push_back({static_cast<int>(s), i});
    }
  }
  size_t n_mod = vert.size();
  // the position of the coordinates of module s in an element of H
  std::vector<size_t> offset(t, 0);
  for (size_t s = 1; s < t; s++) {
    offset[s] = offset[s - 1] + mods[s - 1].rank();
  }
  auto weight = [&](size_t i, size_t j, bool with_glue) -> int64_t {
    (void)with_glue;
    bool gi = i >= n_mod, gj = j >= n_mod;
    if (!gi && !gj) {
      int s = vert[i].first, s2 = vert[j].first;
      ModuleElt const &x = elts[s][vert[i].second];
      if (s != s2) {
        return -2;
      }
      if (i == j) {
        return 4 * (s * BIG + ModuleQ(mods[s], x));
      }
      ModuleElt const &y = elts[s][vert[j].second];
      int64_t flag = 0;
      if (!IsScalarAction(mods[s])) {
        for (size_t b = 0; b < sigma_img[s].size(); b++) {
          if (sigma_img[s][b][vert[i].second] == y) {
            flag |= int64_t(1) << b;
          }
        }
      }
      return 4 * ((s * BIG + ModuleB(mods[s], x, y)) *
                      (int64_t(1) << n_bit) +
                  flag) +
             1;
    }
    if (gi && gj) {
      return (i == j) ? -1 : -5;
    }
    size_t ih = gi ? i - n_mod : j - n_mod;
    size_t im = gi ? j : i;
    int s = vert[im].first;
    ModuleElt const &x = elts[s][vert[im].second];
    ModuleElt const &h = H[ih];
    for (size_t k = 0; k < mods[s].rank(); k++) {
      if (h[offset[s] + k] != x[k]) {
        return -4;
      }
    }
    return -3;
  };
  auto group_order = [&](size_t n, auto f) -> mpz_class {
    WeightMatrix<false, int64_t, Tidx_value> WMat(n, f, os);
    return GetStabilizerWeightMatrix_Kernel<int64_t, Tgr, Tidx, mpz_class,
                                            Tidx_value, false>(WMat, os)
        .order;
  };
  mpz_class order_G(1);
  size_t base = 0;
  for (size_t s = 0; s < t; s++) {
    size_t ns = elts[s].size();
    auto f = [&](size_t i, size_t j) -> int64_t {
      return weight(base + i, base + j, false);
    };
    mpz_class o = group_order(ns, f);
    res.list_d.push_back(mods[s].d);
    res.list_module_order.push_back(mods[s].size());
    res.list_group_order.push_back(o);
    order_G *= o;
    base += ns;
  }
  auto f_all = [&](size_t i, size_t j) -> int64_t { return weight(i, j, true); };
  mpz_class order_stab = group_order(n_vert, f_all);
  if (order_G % order_stab != 0) {
    std::cerr << "ISOM_MASS: the stabilizer of the glue does not divide the "
                 "group\n";
    throw TerminalException{1};
  }
  res.orbit_size = order_G / order_stab;
  return res;
}

inline std::optional<GlueOrbitResult>
GlueOrbitListing(int64_t q, std::vector<PrimaryQuadModule> const &mods,
                 std::vector<std::vector<ModuleElt>> const &glue_gens,
                 int64_t max_module, size_t max_group, size_t max_orbit) {
  GlueOrbitResult res;
  res.q = q;
  res.method = "listing";
  size_t t = mods.size();
  std::vector<int64_t> radix;
  std::vector<size_t> offset;
  for (auto const &M : mods) {
    offset.push_back(radix.size());
    for (int64_t o : M.ord) {
      radix.push_back(o);
    }
  }
  size_t R = radix.size();
  {
    // the codes must fit in int64
    long double tot = 1;
    for (int64_t o : radix) {
      tot *= static_cast<long double>(o);
    }
    if (tot > 4.0e18L) {
      return {};
    }
  }
  auto encode = [&](ModuleElt const &x) -> int64_t {
    int64_t c = 0;
    for (size_t i = R; i-- > 0;) {
      c = c * radix[i] + x[i];
    }
    return c;
  };
  // the elements of the glue
  std::set<ModuleElt> H;
  H.insert(ModuleElt(R, 0));
  {
    std::vector<ModuleElt> work{ModuleElt(R, 0)};
    while (!work.empty()) {
      ModuleElt x = work.back();
      work.pop_back();
      for (auto const &g : glue_gens) {
        ModuleElt y(R);
        for (size_t s = 0; s < t; s++) {
          for (size_t i = 0; i < mods[s].rank(); i++) {
            size_t pos = offset[s] + i;
            y[pos] = ModPos64(x[pos] + g[s][i], radix[pos]);
          }
        }
        if (H.insert(y).second) {
          work.push_back(y);
        }
      }
    }
  }
  res.glue_order = static_cast<int64_t>(H.size());
  std::vector<std::vector<ModuleAut>> gens(t);
  for (size_t s = 0; s < t; s++) {
    std::optional<std::vector<ModuleAut>> opt =
        EquivariantModuleIsometries(mods[s], max_module, max_group);
    if (!opt) {
      return {};
    }
    res.list_d.push_back(mods[s].d);
    res.list_module_order.push_back(mods[s].size());
    res.list_group_order.push_back(mpz_class(static_cast<long>(opt->size())));
    gens[s] = GeneratingSubset(mods[s], *opt);
  }
  std::vector<ModuleElt> Hlist(H.begin(), H.end());
  auto key_of = [&](std::vector<ModuleElt> const &l) -> std::vector<int64_t> {
    std::vector<int64_t> k;
    for (auto const &x : l) {
      k.push_back(encode(x));
    }
    std::sort(k.begin(), k.end());
    return k;
  };
  std::set<std::vector<int64_t>> orbit;
  std::vector<std::vector<ModuleElt>> work{Hlist};
  orbit.insert(key_of(Hlist));
  while (!work.empty()) {
    std::vector<ModuleElt> cur = work.back();
    work.pop_back();
    for (size_t s = 0; s < t; s++) {
      for (auto const &phi : gens[s]) {
        std::vector<ModuleElt> img = cur;
        for (auto &x : img) {
          ModuleElt part(mods[s].rank());
          for (size_t i = 0; i < mods[s].rank(); i++) {
            part[i] = x[offset[s] + i];
          }
          ModuleElt part_img = ApplyModuleAut(mods[s], phi, part);
          for (size_t i = 0; i < mods[s].rank(); i++) {
            x[offset[s] + i] = part_img[i];
          }
        }
        if (orbit.insert(key_of(img)).second) {
          if (orbit.size() > max_orbit) {
            return {};
          }
          work.push_back(img);
        }
      }
    }
  }
  res.orbit_size = mpz_class(static_cast<long>(orbit.size()));
  return res;
}

/*
  The orbit of the glue: by listing the groups when they are small, and
  otherwise by the graph automorphism groups (GlueOrbitGraph).
 */
inline std::optional<GlueOrbitResult>
GlueOrbit(int64_t q, std::vector<PrimaryQuadModule> const &mods,
          std::vector<std::vector<ModuleElt>> const &glue_gens,
          int64_t max_module, size_t max_group, size_t max_orbit,
          size_t max_vertices, std::ostream &os) {
  std::optional<GlueOrbitResult> opt =
      GlueOrbitListing(q, mods, glue_gens, max_module, max_group, max_orbit);
  if (opt) {
#ifdef SANITY_CHECK_ISOMETRY_MASS
    std::optional<GlueOrbitResult> opt2 =
        GlueOrbitGraph(q, mods, glue_gens, max_vertices, os);
    if (opt2 && opt2->orbit_size != opt->orbit_size) {
      std::cerr << "ISOM_MASS: the glue orbit by listing is " << opt->orbit_size
                << " and by the graph " << opt2->orbit_size << "\n";
      throw TerminalException{1};
    }
#endif
    return opt;
  }
  return GlueOrbitGraph(q, mods, glue_gens, max_vertices, os);
}

template <typename T> struct IsometryComponentMass {
  int d;
  int rank;
  T det;
  std::optional<T> mass;
};

template <typename T> struct IsometryMassResult {
  int order;
  std::vector<IsometryComponent<T>> components;
  std::vector<IsometryComponentMass<T>> component_masses;
  T index;
  std::vector<GlueOrbitResult> glue;
  T glue_factor;
  // whether a component with d >= 3 enters the glue, in which case the
  // image of the local unitary group is taken to be all the equivariant
  // isometries of its discriminant form
  bool image_assumed;
  // a component with End_G = Q of dimension > 1 over the multiplicity
  // space: its mass is that of the multiplicity lattice, the local factors
  // at the primes dividing |G| not being computed
  bool local_correction_open = false;
  // a component of a type whose mass is not implemented (End_G neither Q
  // nor reached through a cyclic quotient)
  bool unsupported_component = false;
  std::optional<T> mass;
};

struct IsometryMassBounds {
  int64_t max_module = 1 << 12;
  size_t max_group = 100000;
  size_t max_orbit = 100000;
  // the graph method: the number of vertices (module elements and glue)
  size_t max_vertices = 6000;
};

template <typename T> int IsometryOrder(MyMatrix<T> const &P) {
  int n = P.rows();
  MyMatrix<T> Id = IdentityMat<T>(n);
  MyMatrix<T> pw = P;
  for (int k = 1; k <= 100000; k++) {
    if (pw == Id) {
      return k;
    }
    pw = pw * P;
  }
  std::cerr << "ISOM_MASS: the isometry is not of finite order <= 100000\n";
  throw TerminalException{1};
}

/*
  The component of the saturated sublattice with basis E (rows), on which
  P acts with the d-th cyclotomic polynomial as minimal polynomial.
 */
template <typename T, typename Tint>
IsometryComponent<T> BuildIsometryComponent(MyMatrix<T> const &A,
                                            MyMatrix<T> const &E,
                                            MyMatrix<T> const &P, int d) {
  IsometryComponent<T> comp;
  comp.d = d;
  comp.basis = E;
  comp.gram = E * A * E.transpose();
  comp.sigma = E * P * A * E.transpose() * Inverse(comp.gram);
  if (!IsIntegralMatrix(comp.sigma) || comp.sigma * E != E * P) {
    std::cerr << "ISOM_MASS: the component of d=" << d << " is not stable\n";
    throw TerminalException{1};
  }
  MyMatrix<Tint> G_i = UniversalMatrixConversion<Tint, T>(comp.gram);
  ResultSmithNormalForm<Tint> snf = SmithNormalForm(G_i);
  MyMatrix<Tint> D_i = snf.ROW * G_i * snf.COL;
  comp.row = UniversalMatrixConversion<T, Tint>(snf.ROW);
  comp.row_inv = Inverse(comp.row);
  for (int i = 0; i < D_i.rows(); i++) {
    Tint s = T_abs(D_i(i, i));
    comp.elem.push_back(UniversalScalarConversion<int64_t, Tint>(s));
  }
  return comp;
}

/*
  The rows x of the saturated lattice spanned by the rows of E with
  x Phi_d(P) = 0, in the coordinates of L.
 */
template <typename T, typename Tint>
MyMatrix<T> CyclotomicKernel(MyMatrix<T> const &E, MyMatrix<T> const &P,
                             int d) {
  MyMatrix<T> Phi = EvaluatePolynomialMatrix(CyclotomicPolynomial(d), P);
  MyMatrix<T> EPhi = E * Phi;
  MyMatrix<Tint> EPhi_i = UniversalMatrixConversion<Tint, T>(EPhi);
  // NullspaceIntTrMat gives the x with M x^T = 0
  MyMatrix<Tint> K_i = NullspaceIntTrMat(MyMatrix<Tint>(EPhi_i.transpose()));
  if (K_i.rows() == 0) {
    return MyMatrix<T>(0, E.cols());
  }
  MyMatrix<T> K = UniversalMatrixConversion<T, Tint>(K_i) * E;
  if (!IsZeroMatrix(MyMatrix<T>(K * Phi))) {
    std::cerr << "ISOM_MASS: the kernel of Phi_" << d << " is wrong\n";
    throw TerminalException{1};
  }
  return K;
}

template <typename T, typename Tint>
std::vector<IsometryComponent<T>> IsometryComponents(MyMatrix<T> const &A,
                                                     MyMatrix<T> const &P,
                                                     int m) {
  std::vector<IsometryComponent<T>> res;
  int n = A.rows();
  MyMatrix<T> Id = IdentityMat<T>(n);
  int tot = 0;
  for (int d = 1; d <= m; d++) {
    if (m % d != 0) {
      continue;
    }
    MyMatrix<T> E = CyclotomicKernel<T, Tint>(Id, P, d);
    if (E.rows() == 0) {
      continue;
    }
    tot += E.rows();
    res.push_back(BuildIsometryComponent<T, Tint>(A, E, P, d));
  }
  if (tot != n) {
    std::cerr << "ISOM_MASS: the components have total rank " << tot
              << " instead of " << n << "\n";
    throw TerminalException{1};
  }
  return res;
}

/*
  The elements of the group generated by the ListP (acting on the rows),
  by closure; nothing beyond max_group elements.
 */
template <typename T>
std::vector<MyMatrix<T>> GroupElementsClosure(std::vector<MyMatrix<T>> const &ListP,
                                              size_t max_group) {
  int n = ListP[0].rows();
  auto flatten = [&](MyMatrix<T> const &M) -> std::vector<T> {
    std::vector<T> v;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        v.push_back(M(i, j));
      }
    }
    return v;
  };
  std::vector<MyMatrix<T>> elts{IdentityMat<T>(n)};
  std::set<std::vector<T>> seen{flatten(elts[0])};
  for (size_t pos = 0; pos < elts.size(); pos++) {
    for (auto const &h : ListP) {
      MyMatrix<T> x = elts[pos] * h;
      if (seen.insert(flatten(x)).second) {
        elts.push_back(x);
        if (elts.size() > max_group) {
          std::cerr << "ISOM_MASS: the group has more than " << max_group
                    << " elements\n";
          throw TerminalException{1};
        }
      }
    }
  }
  return elts;
}

/*
  The components of a finite abelian group G. A rational irreducible
  representation of G is the Galois orbit of a character chi, which factors
  through the cyclic quotient G / ker chi of order e, faithfully; the
  orbit is determined by the orders of the chi(h), and V_chi is the
  intersection over h in G of the kernels of Phi_{ord chi(h)}(h). On V_chi
  the group acts through the powers of one element g_chi with chi(g_chi) of
  order e, so the centraliser of G in O(V_chi) is that of g_chi, and the
  component is that of the cyclic case for (L_chi, g_chi) and d = e.
 */
template <typename T, typename Tint>
std::vector<IsometryComponent<T>>
AbelianGroupComponents(MyMatrix<T> const &A,
                       std::vector<MyMatrix<T>> const &elts) {
  int n = A.rows();
  struct Block {
    MyMatrix<T> E;
    std::vector<int> sig;
  };
  std::vector<Block> blocks{{IdentityMat<T>(n), {}}};
  for (auto const &h : elts) {
    int o = IsometryOrder(h);
    std::vector<Block> next;
    for (auto const &blk : blocks) {
      for (int d = 1; d <= o; d++) {
        if (o % d != 0) {
          continue;
        }
        MyMatrix<T> K = CyclotomicKernel<T, Tint>(blk.E, h, d);
        if (K.rows() == 0) {
          continue;
        }
        std::vector<int> sig = blk.sig;
        sig.push_back(d);
        next.push_back({K, sig});
      }
    }
    blocks = std::move(next);
  }
  std::vector<IsometryComponent<T>> res;
  int tot = 0;
  for (auto const &blk : blocks) {
    size_t i_max = 0;
    for (size_t i = 0; i < blk.sig.size(); i++) {
      if (blk.sig[i] > blk.sig[i_max]) {
        i_max = i;
      }
    }
    IsometryComponent<T> comp = BuildIsometryComponent<T, Tint>(
        A, blk.E, elts[i_max], blk.sig[i_max]);
#ifdef SANITY_CHECK_ISOMETRY_MASS
    // every element acts on the component as a power of sigma
    for (auto const &h : elts) {
      MyMatrix<T> h_c = blk.E * h * A * blk.E.transpose() * Inverse(comp.gram);
      MyMatrix<T> pw = IdentityMat<T>(blk.E.rows());
      bool found = false;
      for (int k = 0; k < blk.sig[i_max] && !found; k++) {
        found = (pw == h_c);
        pw = pw * comp.sigma;
      }
      if (!found) {
        std::cerr << "ISOM_MASS: the group is not cyclic on a component\n";
        throw TerminalException{1};
      }
    }
#endif
    tot += blk.E.rows();
    res.push_back(comp);
  }
  if (tot != n) {
    std::cerr << "ISOM_MASS: the components have total rank " << tot
              << " instead of " << n << "\n";
    throw TerminalException{1};
  }
  return res;
}

/*
  The mass of the genus of (L, sigma). herm_mass gives the masses of the
  Hermitian components by d; the total mass is returned only when every
  component has one.
 */
template <typename T>
void CheckIsometryInput(MyMatrix<T> const &A, MyMatrix<T> const &g) {
  int n = A.rows();
  if (!IsSymmetricMatrix(A) || !IsIntegralMatrix(A) || !IsIntegralMatrix(g)) {
    std::cerr << "ISOM_MASS: A must be symmetric integral and g integral\n";
    throw TerminalException{1};
  }
  for (int i = 0; i < n; i++) {
    if (!IsInteger(A(i, i) / T(2))) {
      std::cerr << "ISOM_MASS: the lattice must be even\n";
      throw TerminalException{1};
    }
  }
  if (g.transpose() * A * g != A) {
    std::cerr << "ISOM_MASS: g^T A g is not A\n";
    throw TerminalException{1};
  }
  if (T_abs(DeterminantMat(g)) != T(1)) {
    std::cerr << "ISOM_MASS: g is not in GL_n(Z)\n";
    throw TerminalException{1};
  }
}

template <typename T, typename Tint>
IsometryMassResult<T>
IsometryGenusMassComponents(MyMatrix<T> const &A, int order,
                            std::vector<IsometryComponent<T>> const &components,
                            std::map<int, T> const &herm_mass,
                            IsometryMassBounds const &bounds, std::ostream &os);

template <typename T, typename Tint>
IsometryMassResult<T> IsometryGenusMass(MyMatrix<T> const &A,
                                        MyMatrix<T> const &g,
                                        std::map<int, T> const &herm_mass,
                                        IsometryMassBounds const &bounds,
                                        std::ostream &os) {
  CheckIsometryInput(A, g);
  MyMatrix<T> P = g.transpose();
  int order = IsometryOrder(P);
  return IsometryGenusMassComponents<T, Tint>(
      A, order, IsometryComponents<T, Tint>(A, P, order), herm_mass, bounds,
      os);
}

/*
  The mass from the components, for an isometry or an abelian group: the
  product of the masses of the components and of the glue orbits. herm_mass
  is by d, which identifies the component only for a single isometry.
 */
template <typename T, typename Tint>
IsometryMassResult<T>
IsometryGenusMassComponents(MyMatrix<T> const &A, int order,
                            std::vector<IsometryComponent<T>> const &components,
                            std::map<int, T> const &herm_mass,
                            IsometryMassBounds const &bounds,
                            [[maybe_unused]] std::ostream &os) {
  int n = A.rows();
  IsometryMassResult<T> res;
  res.order = order;
  res.components = components;
  res.image_assumed = false;
  // the component masses
  bool all_mass = true;
  T mass0(1);
  for (auto const &comp : res.components) {
    IsometryComponentMass<T> cm;
    cm.d = comp.d;
    cm.rank = comp.gram.rows();
    cm.det = DeterminantMat(comp.gram);
    if (comp.d == 0) {
      if (comp.multiplicity_gram) {
        cm.mass = GenusMass<T, Tint>(*comp.multiplicity_gram);
        // multiplicity one: the centraliser on the component is {+-1},
        // which preserves every lattice, so no local factor is missing;
        // otherwise the local factors at p | |G|, unless beyond the bounds
        if (comp.multiplicity_factor) {
          *cm.mass *= *comp.multiplicity_factor;
        } else if (comp.multiplicity_gram->rows() > 1) {
          res.local_correction_open = true;
        }
      } else {
        res.unsupported_component = true;
      }
    } else if (comp.d <= 2) {
      cm.mass = GenusMass<T, Tint>(comp.gram);
    } else {
      auto iter = herm_mass.find(comp.d);
      if (iter != herm_mass.end()) {
        cm.mass = iter->second;
      } else {
        // computed exactly except at a ramified prime of E/K in rank >= 2
        std::optional<mpq_class> opt = HermitianMassCyclotomicLattice(
            UniversalMatrixConversion<mpq_class, T>(comp.gram),
            UniversalMatrixConversion<mpq_class, T>(comp.sigma), comp.d);
        if (opt) {
          cm.mass = UniversalScalarConversion<T, mpq_class>(*opt);
        }
      }
    }
    if (cm.mass) {
      mass0 *= *cm.mass;
    } else {
      all_mass = false;
    }
    res.component_masses.push_back(cm);
  }
  // L_0 and the glue
  MyMatrix<T> Concat(n, n);
  {
    int pos = 0;
    for (auto const &comp : res.components) {
      for (int i = 0; i < comp.basis.rows(); i++) {
        for (int j = 0; j < n; j++) {
          Concat(pos, j) = comp.basis(i, j);
        }
        pos++;
      }
    }
  }
  res.index = T_abs(DeterminantMat(Concat));
  MyMatrix<T> ConcatInv = Inverse(Concat);
  // the projections of e_i to the components, in the coordinates of L_d
  std::vector<std::vector<MyVector<T>>> proj(n);
  for (int i = 0; i < n; i++) {
    int pos = 0;
    for (auto const &comp : res.components) {
      int k = comp.basis.rows();
      MyVector<T> v(k);
      for (int j = 0; j < k; j++) {
        v(j) = ConcatInv(i, pos + j);
      }
      proj[i].push_back(v);
      pos += k;
    }
  }
  res.glue_factor = T(1);
  Tint index_i = UniversalScalarConversion<Tint, T>(res.index);
  std::vector<Tint> primes;
  if (index_i > 1) {
    primes = FactorsInt(index_i);
    std::sort(primes.begin(), primes.end());
    primes.erase(std::unique(primes.begin(), primes.end()), primes.end());
  }
  for (Tint const &q_i : primes) {
    int64_t q = UniversalScalarConversion<int64_t, Tint>(q_i);
    std::vector<PrimaryQuadModule> mods;
    std::vector<size_t> comp_pos;
    for (size_t c = 0; c < res.components.size(); c++) {
      PrimaryQuadModule M = BuildPrimaryQuadModule(res.components[c], q);
      if (M.rank() == 0) {
        continue;
      }
      // keep the component only if the glue projects nontrivially to it
      bool nonzero = false;
      for (int i = 0; i < n; i++) {
        ModuleElt a = ModuleCoordinates(M, res.components[c], proj[i][c]);
        for (int64_t x : a) {
          if (x != 0) {
            nonzero = true;
          }
        }
      }
      if (nonzero) {
        mods.push_back(M);
        comp_pos.push_back(c);
      }
    }
    std::vector<std::vector<ModuleElt>> glue_gens;
    for (int i = 0; i < n; i++) {
      std::vector<ModuleElt> g_i;
      for (size_t s = 0; s < mods.size(); s++) {
        size_t c = comp_pos[s];
        g_i.push_back(ModuleCoordinates(mods[s], res.components[c], proj[i][c]));
      }
      glue_gens.push_back(g_i);
    }
    std::optional<GlueOrbitResult> opt =
        GlueOrbit(q, mods, glue_gens, bounds.max_module, bounds.max_group,
                  bounds.max_orbit, bounds.max_vertices, os);
    if (!opt) {
      std::cerr << "ISOM_MASS: the glue at q=" << q
                << " exceeds the bounds of both the listing and the graph\n";
      throw TerminalException{1};
    }
    for (auto const &M : mods) {
      if (!IsScalarAction(M)) {
        res.image_assumed = true;
      }
    }
#ifdef DEBUG_ISOMETRY_MASS
    os << "ISOM_MASS: q=" << q << " |H_q|=" << opt->glue_order
       << " orbit=" << opt->orbit_size << "\n";
#endif
    res.glue_factor *= UniversalScalarConversion<T, mpz_class>(opt->orbit_size);
    res.glue.push_back(*opt);
  }
#ifdef SANITY_CHECK_ISOMETRY_MASS
  {
    T prod(1);
    for (auto const &gl : res.glue) {
      prod *= ScalarFromInt64<T>(gl.glue_order);
    }
    if (prod != res.index) {
      std::cerr << "ISOM_MASS: the glue has order " << prod
                << " but [L : L_0] = " << res.index << "\n";
      throw TerminalException{1};
    }
  }
#endif
  if (all_mass) {
    res.mass = mass0 * res.glue_factor;
  }
  return res;
}

/*
  The mass of the genus of (L, G), G a finite abelian group of isometries
  given by generators g_i (on the columns, g_i^T A g_i = A):

      mass(L, G) = sum_{(L', G)} 1 / |C_{O(L')}(G)|.

  The decomposition of the cyclic case holds with the components L_chi of
  AbelianGroupComponents, and the glue isometries commute with G, that is
  with each g_chi. The order reported is |G|.
 */
template <typename T, typename Tint>
IsometryMassResult<T> AbelianGroupGenusMass(MyMatrix<T> const &A,
                                            std::vector<MyMatrix<T>> const &list_g,
                                            IsometryMassBounds const &bounds,
                                            std::ostream &os) {
  if (list_g.empty()) {
    std::cerr << "ISOM_MASS: the group needs at least one generator\n";
    throw TerminalException{1};
  }
  std::vector<MyMatrix<T>> ListP;
  for (auto const &g : list_g) {
    CheckIsometryInput(A, g);
    ListP.push_back(g.transpose());
  }
  for (size_t i = 0; i < ListP.size(); i++) {
    for (size_t j = i + 1; j < ListP.size(); j++) {
      if (ListP[i] * ListP[j] != ListP[j] * ListP[i]) {
        std::cerr << "ISOM_MASS: the generators " << i << " and " << j
                  << " do not commute\n";
        throw TerminalException{1};
      }
    }
  }
  std::vector<MyMatrix<T>> elts = GroupElementsClosure(ListP, bounds.max_group);
  std::vector<IsometryComponent<T>> components =
      AbelianGroupComponents<T, Tint>(A, elts);
  std::map<int, T> herm_mass;
  return IsometryGenusMassComponents<T, Tint>(
      A, static_cast<int>(elts.size()), components, herm_mass, bounds, os);
}

/*
  The rational isotypic components of a finite group G (elements elts, on
  the rows). The rational class sums R_K = sum of the classes of the g^k,
  k prime to the order of g, act on the isotypic part of a rational
  irreducible rho by an integer, and these integers tell the rho apart, so
  the components are the common eigenspaces of the R_K; their eigenvalues
  are found numerically and checked exactly.

  On a component where the generators commute, G acts through a cyclic
  quotient and the component is that of the cyclic case for the element
  acting with the largest order. Otherwise, End_G(V_rho) = Q is recognised
  as a centre of dimension 1, End_G of dimension m^2, and a common
  eigenspace F of some elements h for eigenvalues +-1 of dimension m: F is
  then the image l (x) U of an idempotent of Q[H] of rank one on W, and X
  = L cap F, with the restriction of the form, is the multiplicity lattice
  (Section on non-abelian groups of the paper). Other components get no
  mass.
 */
template <typename T, typename Tint>
std::vector<IsometryComponent<T>>
FiniteGroupComponents(MyMatrix<T> const &A, std::vector<MyMatrix<T>> const &ListP,
                      std::vector<MyMatrix<T>> const &elts,
                      [[maybe_unused]] std::ostream &os) {
  int n = A.rows();
  auto flatten = [&](MyMatrix<T> const &M) -> std::vector<T> {
    std::vector<T> v;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        v.push_back(M(i, j));
      }
    }
    return v;
  };
  std::map<std::vector<T>, size_t> index;
  for (size_t i = 0; i < elts.size(); i++) {
    index[flatten(elts[i])] = i;
  }
  std::vector<MyMatrix<T>> inv;
  for (auto const &x : elts) {
    inv.push_back(Inverse(x));
  }
  // the rational classes
  size_t n_elt = elts.size();
  std::vector<int> cls(n_elt, -1);
  int n_cls = 0;
  for (size_t i = 0; i < n_elt; i++) {
    if (cls[i] != -1) {
      continue;
    }
    std::vector<size_t> queue{i};
    cls[i] = n_cls;
    while (!queue.empty()) {
      size_t u = queue.back();
      queue.pop_back();
      std::vector<size_t> nb;
      for (size_t j = 0; j < n_elt; j++) {
        nb.push_back(index[flatten(MyMatrix<T>(inv[j] * elts[u] * elts[j]))]);
      }
      // the powers prime to the order
      int o = IsometryOrder(elts[u]);
      MyMatrix<T> pw = elts[u];
      for (int k = 2; k < o; k++) {
        pw = pw * elts[u];
        if (std::gcd(k, o) == 1) {
          nb.push_back(index[flatten(pw)]);
        }
      }
      for (size_t v : nb) {
        if (cls[v] == -1) {
          cls[v] = n_cls;
          queue.push_back(v);
        }
      }
    }
    n_cls++;
  }
  std::vector<MyMatrix<T>> ListR(n_cls, ZeroMatrix<T>(n, n));
  for (size_t i = 0; i < n_elt; i++) {
    ListR[cls[i]] += elts[i];
  }
  // the common eigenspaces, refined class sum by class sum
  struct Block {
    MyMatrix<T> E;
    std::vector<MyMatrix<T>> conditions;
  };
  std::vector<Block> blocks{{IdentityMat<T>(n), {}}};
  for (auto const &Rk : ListR) {
    std::vector<Block> next;
    for (auto const &blk : blocks) {
      int k = blk.E.rows();
      // the action of R_K on the block, numerically
      MyMatrix<T> act =
          blk.E * Rk * blk.E.transpose() *
          Inverse(MyMatrix<T>(blk.E * blk.E.transpose()));
      Eigen::MatrixXd act_d(k, k);
      for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
          act_d(i, j) = UniversalScalarConversion<double, T>(act(i, j));
        }
      }
      Eigen::EigenSolver<Eigen::MatrixXd> es(act_d, false);
      std::set<int64_t> ListLambda;
      for (int i = 0; i < k; i++) {
        ListLambda.insert(static_cast<int64_t>(std::llround(es.eigenvalues()(i).real())));
      }
      int tot = 0;
      for (int64_t lambda : ListLambda) {
        MyMatrix<T> M = Rk - ScalarFromInt64<T>(lambda) * IdentityMat<T>(n);
        MyMatrix<T> EM = blk.E * M;
        MyMatrix<Tint> EM_i = UniversalMatrixConversion<Tint, T>(EM);
        MyMatrix<Tint> K_i = NullspaceIntTrMat(MyMatrix<Tint>(EM_i.transpose()));
        if (K_i.rows() == 0) {
          continue;
        }
        MyMatrix<T> K = UniversalMatrixConversion<T, Tint>(K_i) * blk.E;
        std::vector<MyMatrix<T>> cond = blk.conditions;
        cond.push_back(M);
        next.push_back({K, cond});
        tot += K.rows();
      }
      if (tot != k) {
        std::cerr << "ISOM_MASS: the class sums are not diagonalizable with "
                     "integer eigenvalues on a block\n";
        throw TerminalException{1};
      }
    }
    blocks = std::move(next);
  }
  std::vector<IsometryComponent<T>> res;
  for (auto const &blk : blocks) {
    int k = blk.E.rows();
    MyMatrix<T> Gram = blk.E * A * blk.E.transpose();
    MyMatrix<T> GramInv = Inverse(Gram);
    auto action = [&](MyMatrix<T> const &P) -> MyMatrix<T> {
      return blk.E * P * A * blk.E.transpose() * GramInv;
    };
    std::vector<MyMatrix<T>> ListAct;
    for (auto const &P : ListP) {
      ListAct.push_back(action(P));
    }
    bool commuting = true;
    for (size_t i = 0; i < ListAct.size(); i++) {
      for (size_t j = i + 1; j < ListAct.size(); j++) {
        if (ListAct[i] * ListAct[j] != ListAct[j] * ListAct[i]) {
          commuting = false;
        }
      }
    }
    if (commuting) {
      // a cyclic quotient: the element of largest order on the block
      size_t i_best = 0;
      int o_best = 0;
      for (size_t i = 0; i < n_elt; i++) {
        int o = IsometryOrder(action(elts[i]));
        if (o > o_best) {
          o_best = o;
          i_best = i;
        }
      }
      res.push_back(
          BuildIsometryComponent<T, Tint>(A, blk.E, elts[i_best], o_best));
      continue;
    }
    IsometryComponent<T> comp =
        BuildIsometryComponent<T, Tint>(A, blk.E, ListP[0], 1);
    comp.d = 0;
    for (size_t i = 1; i < ListP.size(); i++) {
      comp.extra_sigma.push_back(action(ListP[i]));
    }
    // End_G of the block and its centre
    auto commutant = [&](std::vector<MyMatrix<T>> const &ListX)
        -> std::vector<MyMatrix<T>> {
      // the Y with X Y = Y X for all X, by the nullspace
      int kk = k * k;
      MyMatrix<T> Eq(kk, kk * static_cast<int>(ListX.size()));
      for (size_t ix = 0; ix < ListX.size(); ix++) {
        MyMatrix<T> const &X = ListX[ix];
        for (int u = 0; u < k; u++) {
          for (int v = 0; v < k; v++) {
            int col = static_cast<int>(ix) * kk + u * k + v;
            for (int w = 0; w < kk; w++) {
              Eq(w, col) = 0;
            }
            // (X Y - Y X)(u, v) = sum_w X(u,w) Y(w,v) - Y(u,w) X(w,v)
            for (int w = 0; w < k; w++) {
              Eq(w * k + v, col) += X(u, w);
              Eq(u * k + w, col) -= X(w, v);
            }
          }
        }
      }
      MyMatrix<T> ker = NullspaceMat(Eq);
      std::vector<MyMatrix<T>> res_c;
      for (int r = 0; r < ker.rows(); r++) {
        MyMatrix<T> Y(k, k);
        for (int u = 0; u < k; u++) {
          for (int v = 0; v < k; v++) {
            Y(u, v) = ker(r, u * k + v);
          }
        }
        res_c.push_back(Y);
      }
      return res_c;
    };
    std::vector<MyMatrix<T>> End = commutant(ListAct);
    // the centre of End_G: the elements of End_G commuting with End_G
    std::vector<MyMatrix<T>> ListBoth = ListAct;
    ListBoth.insert(ListBoth.end(), End.begin(), End.end());
    int c = static_cast<int>(commutant(ListBoth).size());
    int e = static_cast<int>(End.size());
    int m = static_cast<int>(std::llround(std::sqrt(static_cast<double>(e))));
    bool central_q = (c == 1);
    if (central_q && m * m == e) {
      // a common eigenspace of dimension m, for one or two elements
      std::optional<MyMatrix<T>> found;
      MyMatrix<T> found_basis;
      std::vector<std::pair<size_t, int>> cand;
      for (size_t i = 0; i < n_elt; i++) {
        cand.push_back({i, 1});
        cand.push_back({i, -1});
      }
      auto try_conditions =
          [&](std::vector<std::pair<size_t, int>> const &cs) -> bool {
        // the integral kernel of the conditions of the block and of the
        // eigenvalue conditions
        std::vector<MyMatrix<T>> all = blk.conditions;
        for (auto const &c : cs) {
          all.push_back(elts[c.first] -
                        ScalarFromInt64<T>(c.second) * IdentityMat<T>(n));
        }
        MyMatrix<T> Big(n, n * static_cast<int>(all.size()));
        for (size_t a = 0; a < all.size(); a++) {
          for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
              Big(i, static_cast<int>(a) * n + j) = all[a](i, j);
            }
          }
        }
        MyMatrix<Tint> Big_i = UniversalMatrixConversion<Tint, T>(Big);
        MyMatrix<Tint> K_i =
            NullspaceIntTrMat(MyMatrix<Tint>(Big_i.transpose()));
        if (K_i.rows() != m) {
          return false;
        }
        MyMatrix<T> X = UniversalMatrixConversion<T, Tint>(K_i);
        found = X * A * X.transpose();
        found_basis = X;
        return true;
      };
      for (size_t a = 0; a < cand.size() && !found; a++) {
        (void)try_conditions({cand[a]});
      }
      for (size_t a = 0; a < cand.size() && !found; a++) {
        for (size_t b = a + 1; b < cand.size() && !found; b++) {
          if (cand[a].first != cand[b].first) {
            (void)try_conditions({cand[a], cand[b]});
          }
        }
      }
      comp.multiplicity_gram = found;
      if (found) {
        MultiplicityFactorBounds mbounds;
        comp.multiplicity_factor = MultiplicityLocalFactor<T, Tint>(
            A, found_basis, blk.E, elts, mbounds, os);
      }
    }
#ifdef DEBUG_ISOMETRY_MASS
    os << "ISOM_MASS: non-cyclic component of rank " << k << ", dim End_G="
       << e << ", multiplicity lattice "
       << (comp.multiplicity_gram ? "found" : "not found") << "\n";
#endif
    res.push_back(comp);
  }
  return res;
}

/*
  The mass of the genus of (L, G) for a finite group G given by
  generators: AbelianGroupGenusMass when they commute, otherwise the
  components of FiniteGroupComponents with the glue commuting with all the
  generators. The result says when local factors are left out
  (local_correction_open) or a component has no mass.
 */
template <typename T, typename Tint>
IsometryMassResult<T> FiniteGroupGenusMass(MyMatrix<T> const &A,
                                           std::vector<MyMatrix<T>> const &list_g,
                                           IsometryMassBounds const &bounds,
                                           std::ostream &os) {
  std::vector<MyMatrix<T>> ListP;
  for (auto const &g : list_g) {
    CheckIsometryInput(A, g);
    ListP.push_back(g.transpose());
  }
  bool commuting = true;
  for (size_t i = 0; i < ListP.size(); i++) {
    for (size_t j = i + 1; j < ListP.size(); j++) {
      if (ListP[i] * ListP[j] != ListP[j] * ListP[i]) {
        commuting = false;
      }
    }
  }
  if (commuting) {
    return AbelianGroupGenusMass<T, Tint>(A, list_g, bounds, os);
  }
  std::vector<MyMatrix<T>> elts = GroupElementsClosure(ListP, bounds.max_group);
  std::vector<IsometryComponent<T>> components =
      FiniteGroupComponents<T, Tint>(A, ListP, elts, os);
  std::map<int, T> herm_mass;
  return IsometryGenusMassComponents<T, Tint>(
      A, static_cast<int>(elts.size()), components, herm_mass, bounds, os);
}

template <typename T>
void WriteIsometryMassResult(std::ostream &os, IsometryMassResult<T> const &res) {
  os << "return rec(order:=" << res.order << ",\n  components:=[";
  for (size_t c = 0; c < res.component_masses.size(); c++) {
    auto const &cm = res.component_masses[c];
    if (c > 0) {
      os << ",";
    }
    os << "\n    rec(d:=" << cm.d << ", rank:=" << cm.rank << ", det:=" << cm.det;
    if (cm.mass) {
      os << ", mass:=" << *cm.mass;
    }
    os << ", gram:=";
    WriteMatrixGAP(os, res.components[c].gram);
    os << ", sigma:=";
    WriteMatrixGAP(os, res.components[c].sigma);
    os << ")";
  }
  os << "],\n  index:=" << res.index << ",\n  glue:=[";
  for (size_t i = 0; i < res.glue.size(); i++) {
    auto const &gl = res.glue[i];
    if (i > 0) {
      os << ",";
    }
    os << "\n    rec(q:=" << gl.q << ", glue_order:=" << gl.glue_order
       << ", d:=" << StringStdVectorGAP(gl.list_d)
       << ", module_orders:=" << StringStdVectorGAP(gl.list_module_order)
       << ", group_orders:=" << StringStdVectorGAP(gl.list_group_order)
       << ", orbit:=" << gl.orbit_size << ", method:=\"" << gl.method << "\")";
  }
  os << "],\n  glue_factor:=" << res.glue_factor
     << ",\n  image_assumed:=" << GAP_logical(res.image_assumed)
     << ",\n  local_correction_open:="
     << GAP_logical(res.local_correction_open)
     << ",\n  unsupported_component:="
     << GAP_logical(res.unsupported_component);
  if (res.mass) {
    os << ",\n  mass:=" << *res.mass;
  }
  os << ");\n";
}

// clang-format off
#endif  // SRC_GENERA_ISOMETRYMASS_H_
// clang-format on

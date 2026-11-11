// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ISOMETRYMASS_H_
#define SRC_GENERA_ISOMETRYMASS_H_

// clang-format off
#include "LocalDensity.h"
#include "GenusMass.h"
#include "CyclotomicHermitianMass.h"
#include "MultiplicityLocalFactor.h"
#include "AbelianHermitianMass.h"
#include "AbelianQuadraticMass.h"
#include "QuaternionicMass.h"
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
  // End_G with a centre K other than Q: the multiplicity lattice X, with
  // the actions of a Z-basis of O_K, is Hermitian over K when K is CM
  // (herm_field) and quadratic over K when K is totally real (quad_field)
  std::optional<AbelianCMField> herm_field;
  std::optional<AbelianCMField> quad_field;
  MyMatrix<T> centre_gram;
  std::vector<MyMatrix<T>> centre_omega;
  // End_G = M_r(D), D a totally definite quaternion algebra over a real
  // centre: the mass of the component, computed directly
  // (QuaternionicMass.h) when it could be
  bool quaternionic = false;
  std::optional<T> quaternionic_mass;
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
  // the number of elements, for the listings: saturated at the largest
  // size_t, so that a module too large to be listed is never mistaken for a
  // small one
  size_t size() const {
    size_t s = 1;
    for (int64_t o : ord) {
      size_t uo = static_cast<size_t>(o);
      if (uo != 0 && s > std::numeric_limits<size_t>::max() / uo) {
        return std::numeric_limits<size_t>::max();
      }
      s *= uo;
    }
    return s;
  }
  // the order, an unbounded integer
  template <typename Tint> Tint order() const {
    Tint s(1);
    for (int64_t o : ord) {
      s *= UniversalScalarConversion<Tint, int64_t>(o);
    }
    return s;
  }
  size_t rank() const { return ord.size(); }
};

using ModuleElt = std::vector<int64_t>;
// an automorphism, as the images of the generators
using ModuleAut = std::vector<ModuleElt>;

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

/*
  All the isometries of the module commuting with sigma, as the images of
  the generators, by backtracking: the image of g_i has an order dividing
  ord_i, the same q value, and the same b with the images already chosen.
  Such a map is a homomorphism preserving the nondegenerate form, hence an
  automorphism. Nothing is returned when the module or the group exceeds
  the bounds.
 */

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
EquivariantModuleIsometries(PrimaryQuadModule const &M, size_t max_module,
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
template <typename Tint> struct GlueOrbitResult {
  int64_t q;
  // the number of elements of the glue, which are listed
  size_t glue_order;
  std::vector<int> list_d;
  std::vector<Tint> list_module_order;
  std::vector<Tint> list_group_order;
  // per module: "scalar" (Nikulin), "full" (the image of the local group
  // is all the equivariant isometries), "proper" (a proper subgroup, used
  // as such) or "assumed" (not computed, all the isometries used)
  std::vector<std::string> list_image;
  Tint orbit_size;
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
template <typename Tint>
inline std::optional<GlueOrbitResult<Tint>>
GlueOrbitGraph(int64_t q, std::vector<PrimaryQuadModule> const &mods,
               std::vector<std::vector<ModuleElt>> const &glue_gens,
               size_t max_vertices, std::ostream &os) {
  using Tgr = GraphListAdj;
  using Tidx = uint32_t;
  using Tidx_value = int32_t;
  GlueOrbitResult<Tint> res;
  res.q = q;
  res.method = "graph";
  size_t t = mods.size();
  std::vector<ModuleElt> H = GlueElements(mods, glue_gens);
  res.glue_order = H.size();
  std::vector<std::vector<ModuleElt>> elts(t);
  size_t n_vert = H.size();
  int64_t BIG = 1;
  for (size_t s = 0; s < t; s++) {
    if (mods[s].size() > max_vertices) {
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
  auto group_order = [&](size_t n, auto f) -> Tint {
    WeightMatrix<false, int64_t, Tidx_value> WMat(n, f, os);
    return GetStabilizerWeightMatrix_Kernel<int64_t, Tgr, Tidx, Tint,
                                            Tidx_value, false>(WMat, os)
        .order;
  };
  Tint order_G(1);
  size_t base = 0;
  for (size_t s = 0; s < t; s++) {
    size_t ns = elts[s].size();
    auto f = [&](size_t i, size_t j) -> int64_t {
      return weight(base + i, base + j, false);
    };
    Tint o = group_order(ns, f);
    res.list_d.push_back(mods[s].d);
    res.list_module_order.push_back(mods[s].template order<Tint>());
    res.list_group_order.push_back(o);
    order_G *= o;
    base += ns;
  }
  auto f_all = [&](size_t i, size_t j) -> int64_t { return weight(i, j, true); };
  Tint order_stab = group_order(n_vert, f_all);
  if (ResInt(order_G, order_stab) != 0) {
    std::cerr << "ISOM_MASS: the stabilizer of the glue does not divide the "
                 "group\n";
    throw TerminalException{1};
  }
  res.orbit_size = QuoInt(order_G, order_stab);
  return res;
}

template <typename Tint>
inline std::optional<GlueOrbitResult<Tint>>
GlueOrbitListing(int64_t q, std::vector<PrimaryQuadModule> const &mods,
                 std::vector<std::vector<ModuleElt>> const &glue_gens,
                 std::vector<std::optional<std::vector<ModuleAut>>> const
                     &image_groups,
                 size_t max_module, size_t max_group, size_t max_orbit) {
  GlueOrbitResult<Tint> res;
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
  res.glue_order = H.size();
  std::vector<std::vector<ModuleAut>> gens(t);
  for (size_t s = 0; s < t; s++) {
    // the image of the local group when it was computed, otherwise all
    // the equivariant isometries
    std::optional<std::vector<ModuleAut>> opt = image_groups[s];
    if (!opt) {
      opt = EquivariantModuleIsometries(mods[s], max_module, max_group);
    }
    if (!opt) {
      return {};
    }
    res.list_d.push_back(mods[s].d);
    res.list_module_order.push_back(mods[s].template order<Tint>());
    res.list_group_order.push_back(
        UniversalScalarConversion<Tint, int64_t>(static_cast<int64_t>(opt->size())));
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
  res.orbit_size =
      UniversalScalarConversion<Tint, int64_t>(static_cast<int64_t>(orbit.size()));
  return res;
}

/*
  The orbit of the glue: by listing the groups when they are small, and
  otherwise by the graph automorphism groups (GlueOrbitGraph).
 */
template <typename Tint>
inline std::optional<GlueOrbitResult<Tint>>
GlueOrbit(int64_t q, std::vector<PrimaryQuadModule> const &mods,
          std::vector<std::vector<ModuleElt>> const &glue_gens,
          std::vector<std::optional<std::vector<ModuleAut>>> const
              &image_groups,
          bool image_proper, size_t max_module, size_t max_group,
          size_t max_orbit, size_t max_vertices, std::ostream &os) {
  std::optional<GlueOrbitResult<Tint>> opt = GlueOrbitListing<Tint>(
      q, mods, glue_gens, image_groups, max_module, max_group, max_orbit);
  if (opt || image_proper) {
    // the graph method uses all the equivariant isometries, which is not
    // the group when an image is a proper subgroup
#ifdef SANITY_CHECK_ISOMETRY_MASS
    std::optional<GlueOrbitResult<Tint>> opt2 =
        GlueOrbitGraph<Tint>(q, mods, glue_gens, max_vertices, os);
    if (opt && !image_proper && opt2 && opt2->orbit_size != opt->orbit_size) {
      std::cerr << "ISOM_MASS: the glue orbit by listing is " << opt->orbit_size
                << " and by the graph " << opt2->orbit_size << "\n";
      throw TerminalException{1};
    }
#endif
    return opt;
  }
  return GlueOrbitGraph<Tint>(q, mods, glue_gens, max_vertices, os);
}

// The commutant of comp.sigma and comp.extra_sigma in M_k(Z), saturated.
template <typename T, typename Tint>
std::vector<MyMatrix<Tint>> IntegralCommutant(IsometryComponent<T> const &comp) {
  int k = comp.gram.rows();
  std::vector<MyMatrix<T>> ListX{comp.sigma};
  ListX.insert(ListX.end(), comp.extra_sigma.begin(), comp.extra_sigma.end());
  // the integral commutant: Y X = X Y, saturated
  int kk = k * k;
  MyMatrix<Tint> Eq(kk, kk * static_cast<int>(ListX.size()));
  for (int i = 0; i < Eq.rows(); i++) {
    for (int j = 0; j < Eq.cols(); j++) {
      Eq(i, j) = 0;
    }
  }
  for (size_t ix = 0; ix < ListX.size(); ix++) {
    MyMatrix<Tint> X = UniversalMatrixConversion<Tint, T>(ListX[ix]);
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        int col = static_cast<int>(ix) * kk + u * k + v;
        for (int w = 0; w < k; w++) {
          Eq(u * k + w, col) += X(w, v);
          Eq(w * k + v, col) -= X(u, w);
        }
      }
    }
  }
  MyMatrix<Tint> ker = NullspaceIntTrMat(MyMatrix<Tint>(Eq.transpose()));
  std::vector<MyMatrix<Tint>> res;
  for (int r = 0; r < ker.rows(); r++) {
    MyMatrix<Tint> Y(k, k);
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        Y(u, v) = ker(r, u * k + v);
      }
    }
    res.push_back(Y);
  }
  return res;
}

/*
  Elements of the local centraliser by the Cayley transform: for Y in the
  commutant with Y A + A Y^T = 0 and I + Y invertible modulo q, g = (I - Y)
  (I + Y)^{-1} is an isometry (Y A = -A Y^T gives (I + Y) A = A (I - Y^T))
  commuting with the group, with denominators prime to q. Returned as
  integer matrices congruent to g modulo q^c, which is all the action on
  D(L)_q sees. count random Y with small coefficients, from a fixed seed.
 */
template <typename T, typename Tint>
std::vector<MyMatrix<Tint>> CayleyLocalElements(IsometryComponent<T> const &comp,
                                                int64_t q, int c, size_t count) {
  int k = comp.gram.rows();
  std::vector<MyMatrix<Tint>> ListComm = IntegralCommutant<T, Tint>(comp);
  size_t dim = ListComm.size();
  MyMatrix<Tint> A = UniversalMatrixConversion<Tint, T>(comp.gram);
  // the skew part: x with sum x_i (Y_i A + A Y_i^T) = 0
  int n_eq = k * k;
  MyMatrix<Tint> Eq(static_cast<int>(dim), n_eq);
  for (size_t i = 0; i < dim; i++) {
    MyMatrix<Tint> Z = ListComm[i] * A + A * ListComm[i].transpose();
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        Eq(static_cast<int>(i), u * k + v) = Z(u, v);
      }
    }
  }
  MyMatrix<Tint> ker = NullspaceIntTrMat(MyMatrix<Tint>(Eq.transpose()));
  std::vector<MyMatrix<Tint>> ListSkew;
  for (int r = 0; r < ker.rows(); r++) {
    MyMatrix<Tint> Y = ZeroMatrix<Tint>(k, k);
    for (size_t i = 0; i < dim; i++) {
      Tint coef = ker(r, static_cast<int>(i));
      if (coef != 0) {
        MyMatrix<Tint> term = ListComm[i] * coef;
        Y += term;
      }
    }
    ListSkew.push_back(Y);
  }
  std::vector<MyMatrix<Tint>> res;
  if (ListSkew.empty()) {
    return res;
  }
  Tint q_T = UniversalScalarConversion<Tint, int64_t>(q);
  Tint qc(1);
  for (int i = 0; i < c; i++) {
    qc *= q_T;
  }
  std::mt19937_64 rng(static_cast<uint64_t>(q) * 2654435761ULL + k);
  std::uniform_int_distribution<int> unif(-2, 2);
  MyMatrix<Tint> Id = IdentityMat<Tint>(k);
  for (size_t iter = 0; iter < 20 * count && res.size() < count; iter++) {
    MyMatrix<Tint> Y = ZeroMatrix<Tint>(k, k);
    for (auto const &S : ListSkew) {
      int coef = unif(rng);
      if (coef != 0) {
        MyMatrix<Tint> term = S * Tint(coef);
        Y += term;
      }
    }
    MyMatrix<Tint> Ip = Id + Y;
    Tint det = DeterminantMat(Ip);
    Tint r = ((det % q_T) + q_T) % q_T;
    if (r == 0) {
      continue;
    }
    // (I - Y) adj(I + Y) det^{-1} modulo q^c
    MyMatrix<T> IpInv = Inverse(UniversalMatrixConversion<T, Tint>(Ip));
    MyMatrix<Tint> adj = UniversalMatrixConversion<Tint, T>(
        MyMatrix<T>(IpInv * UniversalScalarConversion<T, Tint>(det)));
    Tint dm = ((det % qc) + qc) % qc;
    int64_t dinv = InvModQ(UniversalScalarConversion<int64_t, Tint>(dm),
                           UniversalScalarConversion<int64_t, Tint>(qc));
    MyMatrix<Tint> ImY = Id - Y;
    MyMatrix<Tint> g = ImY * adj;
    Tint dinv_T = UniversalScalarConversion<Tint, int64_t>(dinv);
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        Tint x = g(u, v) * dinv_T;
        g(u, v) = ((x % qc) + qc) % qc;
      }
    }
    res.push_back(g);
  }
  return res;
}

/*
  The image in the automorphisms of the q-primary module M of D(L_d) of
  the stabiliser of L_d in the centraliser C(Q_q) of the group acting on
  the component: the elements of the commutant Lambda of comp.sigma and
  comp.extra_sigma that preserve the form, found modulo q^c (q^c the
  exponent of M) by LocalStabilizerImage, then moved to M. Nothing beyond
  the bounds.
 */
template <typename T, typename Tint>
std::optional<std::vector<ModuleAut>>
LocalImageOnModule(IsometryComponent<T> const &comp, PrimaryQuadModule const &M,
                   [[maybe_unused]] std::ostream &os) {
  int k = comp.gram.rows();
  int kk = k * k;
  std::vector<MyMatrix<Tint>> ListComm = IntegralCommutant<T, Tint>(comp);
  MyMatrix<Tint> ker(static_cast<int>(ListComm.size()), kk);
  for (size_t r = 0; r < ListComm.size(); r++) {
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        ker(static_cast<int>(r), u * k + v) = ListComm[r](u, v);
      }
    }
  }
  int64_t q = M.q;
  int c = 0;
  for (int64_t x = M.N; x > 1; x /= q) {
    c++;
  }
  int e = 0;
  {
    MyMatrix<T> Ainv = Inverse(comp.gram);
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < k; j++) {
        if (Ainv(i, j) != T(0)) {
          e = std::max(e, ValuationOfRational(GetDenominator(Ainv(i, j)), q));
        }
      }
    }
  }
  // exact integer matrices: the solutions must stay in Lambda
  auto exact = [&](MyMatrix<Tint> const &X) -> ModMat {
    ModMat R(k, std::vector<int64_t>(k));
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < k; j++) {
        R[i][j] = UniversalScalarConversion<int64_t, Tint>(X(i, j));
      }
    }
    return R;
  };
  std::vector<ModMat> basis;
  for (int r = 0; r < ker.rows(); r++) {
    MyMatrix<Tint> Y(k, k);
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        Y(u, v) = ker(r, u * k + v);
      }
    }
    basis.push_back(exact(Y));
  }
  ModMat A_mod = exact(UniversalMatrixConversion<Tint, T>(comp.gram));
  MultiplicityFactorBounds mbounds;
  std::optional<std::vector<ModMat>> opt =
      LocalStabilizerImage(A_mod, basis, q, c, e, mbounds);
  if (!opt) {
    return {};
  }
  std::vector<MyVector<T>> gens = PrimaryModuleGenerators(comp, M);
  std::set<ModuleAut> images;
  for (auto const &g : *opt) {
    MyMatrix<Tint> h(k, k);
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < k; j++) {
        h(i, j) = UniversalScalarConversion<Tint, int64_t>(g[i][j]);
      }
    }
    images.insert(ModuleActionOfIsometry(comp, M, gens, h));
  }
#ifdef DEBUG_ISOMETRY_MASS
  os << "ISOM_MASS: q=" << q << " d=" << comp.d << " dim Lambda="
     << basis.size() << " solutions mod q^c=" << opt->size()
     << " image on D=" << images.size() << "\n";
#endif
  return std::vector<ModuleAut>(images.begin(), images.end());
}

// The order of the group generated by the automorphisms gens of M, by
// closure, stopping beyond max_size.
inline size_t ModuleGroupClosureSize(PrimaryQuadModule const &M,
                                     std::vector<ModuleAut> const &gens,
                                     size_t max_size) {
  size_t r = M.rank();
  ModuleAut id(r, ModuleElt(r, 0));
  for (size_t i = 0; i < r; i++) {
    id[i][i] = 1;
  }
  std::set<ModuleAut> seen{id};
  std::vector<ModuleAut> queue{id};
  while (!queue.empty()) {
    ModuleAut x = queue.back();
    queue.pop_back();
    for (auto const &g : gens) {
      // x then g: the images of the generators by x, moved by g
      ModuleAut y;
      for (auto const &xi : x) {
        y.push_back(ApplyModuleAut(M, g, xi));
      }
      if (seen.insert(y).second) {
        if (seen.size() > max_size) {
          return seen.size();
        }
        queue.push_back(y);
      }
    }
  }
  return seen.size();
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
  std::vector<GlueOrbitResult<typename underlying_z_ring<T>::ring_type>> glue;
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
  size_t max_module = size_t(1) << 12;
  size_t max_group = 100000;
  size_t max_orbit = 100000;
  // the graph method: the number of vertices (module elements and glue)
  size_t max_vertices = 6000;
  // the L-values of the quadratic components with E not abelian
  QuadLValueMethod lvalue = QuadLValueMethod::Mpfr;
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
    if (comp.d == 0 && comp.quaternionic) {
      if (comp.quaternionic_mass) {
        cm.mass = *comp.quaternionic_mass;
      } else {
        res.unsupported_component = true;
      }
    } else if (comp.d == 0 && (comp.herm_field || comp.quad_field)) {
      std::optional<T> hm =
          comp.herm_field
              ? HermitianMassAbelianLattice(comp.centre_gram, comp.centre_omega,
                                            *comp.herm_field)
              : QuadraticMassAbelianLattice(comp.centre_gram, comp.centre_omega,
                                            *comp.quad_field, bounds.lvalue);
      if (hm && comp.multiplicity_factor) {
        cm.mass = *hm * *comp.multiplicity_factor;
      } else if (hm) {
        cm.mass = *hm;
        res.local_correction_open = true;
      } else {
        res.unsupported_component = true;
      }
    } else if (comp.d == 0) {
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
        std::optional<T> opt =
            HermitianMassCyclotomicLattice(comp.gram, comp.sigma, comp.d);
        if (opt) {
          cm.mass = *opt;
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
    // the images of the local groups on the non-scalar modules
    std::vector<std::optional<std::vector<ModuleAut>>> image_groups;
    std::vector<std::string> list_image;
    bool image_proper = false;
    for (size_t s = 0; s < mods.size(); s++) {
      if (IsScalarAction(mods[s])) {
        image_groups.push_back({});
        list_image.push_back("scalar");
        continue;
      }
      std::optional<std::vector<ModuleAut>> full = EquivariantModuleIsometries(
          mods[s], bounds.max_module, bounds.max_group);
      if (full) {
        // a lower bound first: -1 and, for a cyclic component, sigma are
        // in the local group; if they generate all the equivariant
        // isometries the image is full without a search
        IsometryComponent<T> const &comp = res.components[comp_pos[s]];
        std::vector<MyVector<T>> gens = PrimaryModuleGenerators(comp, mods[s]);
        int kc = comp.gram.rows();
        std::vector<ModuleAut> lower{ModuleActionOfIsometry(
            comp, mods[s], gens,
            MyMatrix<Tint>(-IdentityMat<Tint>(kc)))};
        if (comp.extra_sigma.empty()) {
          lower.push_back(ModuleActionOfIsometry(
              comp, mods[s], gens,
              UniversalMatrixConversion<Tint, T>(comp.sigma)));
        }
        bool is_full =
            ModuleGroupClosureSize(mods[s], lower, full->size()) ==
            full->size();
        if (!is_full) {
          // and elements of the local group by the Cayley transform
          int c = 0;
          for (int64_t x = mods[s].N; x > 1; x /= mods[s].q) {
            c++;
          }
          for (auto const &g :
               CayleyLocalElements<T, Tint>(comp, mods[s].q, c, 12)) {
            lower.push_back(ModuleActionOfIsometry(comp, mods[s], gens, g));
          }
          is_full = ModuleGroupClosureSize(mods[s], lower, full->size()) ==
                    full->size();
        }
        if (is_full) {
          image_groups.push_back(full);
          list_image.push_back("full");
          continue;
        }
      }
      std::optional<std::vector<ModuleAut>> img =
          LocalImageOnModule<T, Tint>(res.components[comp_pos[s]], mods[s], os);
      if (!img) {
        image_groups.push_back({});
        list_image.push_back("assumed");
        res.image_assumed = true;
        continue;
      }
      if (full && full->size() == img->size()) {
        // the whole group: the graph method may be used as well
        image_groups.push_back(full);
        list_image.push_back("full");
      } else {
        image_groups.push_back(img);
        list_image.push_back("proper");
        image_proper = true;
      }
    }
    std::optional<GlueOrbitResult<typename underlying_z_ring<T>::ring_type>> opt =
        GlueOrbit<typename underlying_z_ring<T>::ring_type>(
        q, mods, glue_gens, image_groups, image_proper, bounds.max_module,
        bounds.max_group, bounds.max_orbit, bounds.max_vertices, os);
    if (!opt) {
      std::cerr << "ISOM_MASS: the glue at q=" << q
                << " exceeds the bounds of both the listing and the graph\n";
      throw TerminalException{1};
    }
    opt->list_image = list_image;
#ifdef DEBUG_ISOMETRY_MASS
    os << "ISOM_MASS: q=" << q << " |H_q|=" << opt->glue_order
       << " orbit=" << opt->orbit_size << "\n";
#endif
    res.glue_factor *= UniversalScalarConversion<T, typename underlying_z_ring<T>::ring_type>(
        opt->orbit_size);
    res.glue.push_back(*opt);
  }
#ifdef SANITY_CHECK_ISOMETRY_MASS
  {
    T prod(1);
    for (auto const &gl : res.glue) {
      prod *= UniversalScalarConversion<T, int64_t>(static_cast<int64_t>(gl.glue_order));
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
  Reflections of X_p, X with Gram matrix A: for v in X with 2 (x, v) /
  (v, v) p-integral for all x in X, s_v : x -> x - 2 (x, v) / (v, v) v is
  in O(X_p); it supplies the elements of determinant -1 that the Cayley
  transforms miss. The v are taken with coordinates in {-1, 0, 1}, up to
  count of them; returned modulo p^c, rows convention.
 */
template <typename T, typename Tint>
std::vector<MyMatrix<Tint>> LocalReflections(MyMatrix<T> const &A, int64_t p,
                                             int c, size_t count) {
  int m = A.rows();
  T p_T = UniversalScalarConversion<T, int64_t>(p);
  Tint pc(1);
  for (int i = 0; i < c; i++) {
    pc *= UniversalScalarConversion<Tint, int64_t>(p);
  }
  auto p_integral = [&](T const &x) -> bool {
    Tint den = GetDenominator_z(x);
    return ResInt(den, UniversalScalarConversion<Tint, int64_t>(p)) != 0;
  };
  auto reduce = [&](T const &x) -> Tint {
    Tint nm = ResInt(GetNumerator_z(x), pc);
    Tint dn = ResInt(GetDenominator_z(x), pc);
    Tint prod = nm * mod_inv(dn, pc);
    return ResInt(prod, pc);
  };
  (void)p_T;
  std::vector<MyMatrix<Tint>> res;
  size_t n_cand = 1;
  for (int i = 0; i < m; i++) {
    n_cand *= 3;
  }
  for (size_t u = 1; u < n_cand && res.size() < count; u++) {
    MyVector<T> v(m);
    size_t t = u;
    for (int i = 0; i < m; i++) {
      v(i) = T(static_cast<int>(t % 3) - 1);
      t /= 3;
    }
    MyVector<T> Av = A * v;
    T norm = v.dot(Av);
    if (norm == T(0)) {
      continue;
    }
    bool ok = true;
    for (int i = 0; i < m && ok; i++) {
      ok = p_integral(T(2) * Av(i) / norm);
    }
    if (!ok) {
      continue;
    }
    // rows: x -> x - (2 x A v^T / (v A v^T)) v, matrix I - (2/norm) (A v) v
    MyMatrix<Tint> h(m, m);
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++) {
        T val = T(i == j ? 1 : 0) - T(2) * Av(i) * v(j) / norm;
        h(i, j) = reduce(val);
      }
    }
    res.push_back(h);
  }
  return res;
}

#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
#define MULT_FACTOR_FAIL(tag) \
  std::cerr << "MULT_FACTOR: gives up at check " << (tag) << "\n"
#else
#define MULT_FACTOR_FAIL(tag)
#endif

// The size of the orbit of the subgroup spanned by gens under the
// automorphisms group_gens of M; nothing beyond max_size subgroups.
inline std::optional<size_t>
ModuleSubgroupOrbitSize(PrimaryQuadModule const &M,
                        std::vector<ModuleElt> const &gens,
                        std::vector<ModuleAut> const &group_gens,
                        size_t max_size) {
  size_t r = M.rank();
  auto encode = [&](ModuleElt const &x) -> size_t {
    size_t idx = 0;
    for (size_t i = r; i-- > 0;) {
      idx = idx * static_cast<size_t>(M.ord[i]) + static_cast<size_t>(x[i]);
    }
    return idx;
  };
  auto span = [&](std::vector<ModuleElt> const &g) -> std::vector<size_t> {
    std::set<size_t> seen{0};
    std::vector<ModuleElt> queue{ModuleElt(r, 0)};
    std::vector<ModuleElt> elts;
    while (!queue.empty()) {
      ModuleElt x = queue.back();
      queue.pop_back();
      for (auto const &y : g) {
        ModuleElt z(r);
        for (size_t i = 0; i < r; i++) {
          z[i] = ModPos64(x[i] + y[i], M.ord[i]);
        }
        if (seen.insert(encode(z)).second) {
          queue.push_back(z);
        }
      }
    }
    return std::vector<size_t>(seen.begin(), seen.end());
  };
  // the orbit acts on generating sets: the image of a subgroup is spanned
  // by the images of its generators
  std::set<std::vector<size_t>> orbit;
  std::vector<std::vector<ModuleElt>> queue{gens};
  orbit.insert(span(gens));
  while (!queue.empty()) {
    std::vector<ModuleElt> g = queue.back();
    queue.pop_back();
    for (auto const &phi : group_gens) {
      std::vector<ModuleElt> img;
      for (auto const &x : g) {
        img.push_back(ApplyModuleAut(M, phi, x));
      }
      if (orbit.insert(span(img)).second) {
        if (orbit.size() > max_size) {
          return {};
        }
        queue.push_back(img);
      }
    }
  }
  return orbit.size();
}

/*
  The product over the primes of [O(X_p) : Stab(L_p)].
  A: the Gram matrix of L; Xb: a basis of X (rows, coordinates of L); E: a
  basis of L_rho = L cap V_rho; elts: the elements of G (on the rows).
  Nothing when the bounds are exceeded.
 */
template <typename T, typename Tint>
std::optional<T> MultiplicityLocalFactor(MyMatrix<T> const &A,
                                         MyMatrix<T> const &Xb,
                                         MyMatrix<T> const &E,
                                         std::vector<MyMatrix<T>> const &elts,
                                         std::vector<MyMatrix<T>> const &ListP,
                                         std::vector<MyMatrix<T>> const &OmegaX,
                                         std::vector<MyMatrix<T>> const &ListCentral,
                                         MultiplicityFactorBounds const &bounds,
                                         [[maybe_unused]] std::ostream &os) {
  int m = Xb.rows();
  int k = E.rows();
  int n = A.cols();
  if (k % m != 0) {
    return {};
  }
  int w = k / m;
  // the gamma_a with the gamma_a X in direct sum (V_rho = W (x)_K U and
  // X in l (x) U, l a line of W over the centre K)
  std::vector<size_t> ListGamma;
  {
    MyMatrix<T> rows(0, n);
    for (size_t ig = 0; ig < elts.size() && static_cast<int>(ListGamma.size()) < w; ig++) {
      MyMatrix<T> Y = Xb * elts[ig];
      MyMatrix<T> test(rows.rows() + m, n);
      for (int i = 0; i < rows.rows(); i++) {
        for (int j = 0; j < n; j++) {
          test(i, j) = rows(i, j);
        }
      }
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          test(rows.rows() + i, j) = Y(i, j);
        }
      }
      if (RankMat(test) == test.rows()) {
        rows = test;
        ListGamma.push_back(ig);
      }
    }
  }
  if (static_cast<int>(ListGamma.size()) != w) {
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    os << "MULT_FACTOR: only " << ListGamma.size() << " of the " << w
       << " translates of X found in direct sum\n";
#endif
    return {};
  }
  // the basis B_V = {gamma_a x_i}, index a * m + i
  MyMatrix<T> BV(k, n);
  for (int a = 0; a < w; a++) {
    MyMatrix<T> Y = Xb * elts[ListGamma[a]];
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        BV(a * m + i, j) = Y(i, j);
      }
    }
  }
  if (RankMat(BV) != k) {
    MULT_FACTOR_FAIL(1872);
    return {};
  }
  // coordinates in B_V of vectors of V_rho
  auto coords = [&](MyMatrix<T> const &M) -> MyMatrix<T> {
    // M = C BV, C = M BV^T (BV BV^T)^{-1}
    MyMatrix<T> BB = BV * BV.transpose();
    return M * BV.transpose() * Inverse(BB);
  };
  // N = Z[G] X
  MyMatrix<T> Ngen(static_cast<int>(elts.size()) * m, n);
  for (size_t ig = 0; ig < elts.size(); ig++) {
    MyMatrix<T> Y = Xb * elts[ig];
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        Ngen(static_cast<int>(ig) * m + i, j) = Y(i, j);
      }
    }
  }
  MyMatrix<Tint> Ngen_i = UniversalMatrixConversion<Tint, T>(Ngen);
  MyMatrix<T> Nb = UniversalMatrixConversion<T, Tint>(GetZbasis(Ngen_i));
  if (Nb.rows() != k) {
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    os << "MULT_FACTOR: Z[G] X has rank " << Nb.rows() << " instead of " << k
       << "\n";
#endif
    return {};
  }
  MyMatrix<T> GramN = Nb * A * Nb.transpose();
  MyMatrix<T> Nsharp = Inverse(GramN) * Nb;
  MyMatrix<T> N_c = coords(Nb);
  MyMatrix<T> Ns_c = coords(Nsharp);
  MyMatrix<T> L_c = coords(E);
  // N in the coordinates of N^#: R = Nb Nsharp^{-1}, integral
  MyMatrix<T> Rm = N_c * Inverse(Ns_c);
  MyMatrix<Tint> R_i = UniversalMatrixConversion<Tint, T>(Rm);
  ResultSmithNormalForm<Tint> snf = SmithNormalForm(R_i);
  MyMatrix<Tint> Dm = snf.ROW * R_i * snf.COL;
  MyMatrix<T> COL_T = UniversalMatrixConversion<T, Tint>(snf.COL);
  // f = COL^{-1} Nsharp: the coordinates of x in N^# are x_F COL
  MyMatrix<T> F_c = Inverse(COL_T) * Ns_c;
  MyMatrix<T> F_inv = Inverse(F_c);
  // the index [L_rho : N]
  T index = T_abs(DeterminantMat(MyMatrix<T>(N_c * Inverse(L_c))));
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
  os << "MULT_FACTOR: [L_rho : N] = " << index << ", multiplicity " << m
     << "\n";
#endif
  T factor(1);
  if (index == T(1)) {
    return factor;
  }
  std::vector<Tint> primes =
      FactorsInt(UniversalScalarConversion<Tint, T>(index));
  std::sort(primes.begin(), primes.end());
  primes.erase(std::unique(primes.begin(), primes.end()), primes.end());
  MyMatrix<T> GramX = Xb * A * Xb.transpose();
  for (auto const &p_i : primes) {
    int64_t p = UniversalScalarConversion<int64_t, Tint>(p_i);
    // the p-part of N^# / N: coordinates i with p | d_i, orders p^{v_i}
    std::vector<int> ListIdx;
    std::vector<int64_t> ListOrd;
    int c = 0;
    for (int i = 0; i < k; i++) {
      T d = T_abs(UniversalScalarConversion<T, Tint>(Dm(i, i)));
      if (d == T(0)) {
        MULT_FACTOR_FAIL(1933);
        return {};
      }
      int v = ValuationOfRational(d, p);
      if (v > 0) {
        ListIdx.push_back(i);
        int64_t o = 1;
        for (int u = 0; u < v; u++) {
          o *= p;
        }
        ListOrd.push_back(o);
        c = std::max(c, v);
      }
    }
    size_t r = ListIdx.size();
    // a vector of V_rho (B_V coordinates, p-integral in N^#) to its class
    auto residue = [&](MyVector<T> const &x) -> std::vector<int64_t> {
      MyVector<T> y = F_inv.transpose() * x;
      std::vector<int64_t> res(r);
      for (size_t s = 0; s < r; s++) {
        T val = y(ListIdx[s]);
        // p-integral: numerator times the inverse of the denominator
        Tint o_T = UniversalScalarConversion<Tint, int64_t>(ListOrd[s]);
        Tint nm = ResInt(GetNumerator_z(val), o_T);
        Tint dn = ResInt(GetDenominator_z(val), o_T);
        int64_t nm64 = UniversalScalarConversion<int64_t, Tint>(nm);
        int64_t dn64 = UniversalScalarConversion<int64_t, Tint>(dn);
        res[s] = MulModQ(nm64, InvModQ(dn64, ListOrd[s]), ListOrd[s]);
      }
      return res;
    };
    // |N^#/N|_p, saturated: only the element lists of the sanity check
    // need it to fit (the orbit itself is on Howell forms)
    size_t q_size = 1;
    bool q_fits = true;
    for (auto o : ListOrd) {
      size_t uo = static_cast<size_t>(o);
      if (uo != 0 && q_size > bounds.max_module / uo) {
        q_fits = false;
        q_size = bounds.max_module;
      } else {
        q_size *= uo;
      }
    }
    (void)q_fits;
    auto encode = [&](std::vector<int64_t> const &x) -> size_t {
      size_t idx = 0;
      for (size_t s = r; s-- > 0;) {
        idx = idx * static_cast<size_t>(ListOrd[s]) + static_cast<size_t>(x[s]);
      }
      return idx;
    };
    // the subgroup generated by vectors, as the sorted set of its elements
    auto span = [&](std::vector<std::vector<int64_t>> const &gens)
        -> std::vector<size_t> {
      std::set<size_t> seen;
      std::vector<std::vector<int64_t>> queue{std::vector<int64_t>(r, 0)};
      seen.insert(0);
      while (!queue.empty()) {
        std::vector<int64_t> x = queue.back();
        queue.pop_back();
        for (auto const &g : gens) {
          std::vector<int64_t> y(r);
          for (size_t s = 0; s < r; s++) {
            y[s] = (x[s] + g[s]) % ListOrd[s];
          }
          if (seen.insert(encode(y)).second) {
            queue.push_back(y);
          }
        }
      }
      return std::vector<size_t>(seen.begin(), seen.end());
    };
    // L / N
    std::vector<std::vector<int64_t>> Lgens;
    for (int i = 0; i < k; i++) {
      Lgens.push_back(residue(GetMatrixRow(L_c, i)));
    }
    /*
      The bounds first: the orbit of L/N under elements of the image (-1
      and Cayley transforms of O(X_p), acting by 1 (x) g) is at most the
      factor, and the orbit under all the G-equivariant isometries of the
      p-part of D(N), which contain the image, at least. When they agree
      that is the factor, without the search of the solutions.
     */
    auto bounds_factor = [&](size_t max_dn,
                             size_t max_group) -> std::optional<size_t> {
      IsometryComponent<T> compN =
          BuildIsometryComponent<T, Tint>(A, Nb, ListP[0], 1);
      compN.d = 0;
      MyMatrix<T> GramNinv = Inverse(GramN);
      for (size_t i = 1; i < ListP.size(); i++) {
        compN.extra_sigma.push_back(Nb * ListP[i] * A * Nb.transpose() *
                                    GramNinv);
      }
      for (auto const &Zc : ListCentral) {
        compN.extra_sigma.push_back(Nb * Zc * A * Nb.transpose() * GramNinv);
      }
      PrimaryQuadModule MN = BuildPrimaryQuadModule(compN, p);
      std::optional<std::vector<ModuleAut>> upper =
          EquivariantModuleIsometries(MN, max_dn, max_group);
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
      if (!upper) {
        os << "MULT_FACTOR: p=" << p << " the equivariant isometries of the "
           << "p-part of D(N), of exponent " << MN.N << " and order "
           << MN.size() << ", exceed the bound\n";
      }
#endif
      if (upper) {
        std::vector<MyVector<T>> gensMN = PrimaryModuleGenerators(compN, MN);
        MyMatrix<T> NbInv_c = Inverse(N_c);
        std::vector<ModuleElt> LgensMN;
        for (int i = 0; i < k; i++) {
          MyVector<T> x = NbInv_c.transpose() * GetMatrixRow(L_c, i);
          LgensMN.push_back(ModuleCoordinates(MN, compN, x));
        }
        // 1 (x) h on N, h on X (rows convention), as a matrix in N coords
        auto on_N = [&](MyMatrix<Tint> const &h) -> MyMatrix<Tint> {
          MyMatrix<T> h_T = UniversalMatrixConversion<T, Tint>(h);
          MyMatrix<T> img(k, k);
          for (int r = 0; r < k; r++) {
            for (int a = 0; a < w; a++) {
              for (int j = 0; j < m; j++) {
                T sum(0);
                for (int i = 0; i < m; i++) {
                  sum += N_c(r, a * m + i) * h_T(i, j);
                }
                img(r, a * m + j) = sum;
              }
            }
          }
          return UniversalMatrixConversion<Tint, T>(MyMatrix<T>(img * NbInv_c));
        };
        IsometryComponent<T> compX = BuildIsometryComponent<T, Tint>(
            A, Xb, IdentityMat<T>(A.rows()), 1);
        // the Cayley transforms commute with O_K when it acts on X
        compX.extra_sigma = OmegaX;
        std::vector<ModuleAut> lower{ModuleActionOfIsometry(
            compN, MN, gensMN, MyMatrix<Tint>(-IdentityMat<Tint>(k)))};
        int cN = 0;
        for (int64_t x = MN.N; x > 1; x /= p) {
          cN++;
        }
        for (auto const &h : CayleyLocalElements<T, Tint>(compX, p, cN, 16)) {
          lower.push_back(ModuleActionOfIsometry(compN, MN, gensMN, on_N(h)));
        }
        if (OmegaX.empty()) {
          for (auto const &h : LocalReflections<T, Tint>(GramX, p, cN, 24)) {
            lower.push_back(ModuleActionOfIsometry(compN, MN, gensMN, on_N(h)));
          }
        }
        std::optional<size_t> orb_low = ModuleSubgroupOrbitSize(
            MN, LgensMN, lower, bounds.max_module);
        std::optional<size_t> orb_up = ModuleSubgroupOrbitSize(
            MN, LgensMN, *upper, bounds.max_module);
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
        os << "MULT_FACTOR: p=" << p << " bounds: orbit under the lower "
           << "group " << (orb_low ? static_cast<long>(*orb_low) : -1)
           << ", under the upper " << (orb_up ? static_cast<long>(*orb_up) : -1)
           << "\n";
#endif
        if (orb_low && orb_up && *orb_low == *orb_up) {
          return *orb_low;
        }
      }
      return {};
    };
    // with a centre K other than Q the search over the primes of K over
    // p, pi-adically, first: it is exact and the fastest
    std::optional<std::vector<ModMat>> opt;
    bool rows_convention = !OmegaX.empty();
    if (!OmegaX.empty()) {
      opt = LocalImagePiAdic(GramX, OmegaX, p, c, bounds.max_solution);
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
      os << "MULT_FACTOR: p=" << p << " pi-adic image "
         << (opt ? static_cast<long>(opt->size()) : -1L) << " elements\n";
#endif
    }
    // the image of O(X_p) modulo p^c
    int e = 0;
    {
      MyMatrix<T> Ainv = Inverse(GramX);
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
          if (Ainv(i, j) != T(0)) {
            T den = GetDenominator(Ainv(i, j));
            e = std::max(e, ValuationOfRational(den, p));
          }
        }
      }
    }
    ModMat A_mod(m, std::vector<int64_t>(m));
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++) {
        A_mod[i][j] = UniversalScalarConversion<int64_t, T>(GramX(i, j));
      }
    }
    // which construction gave the image (for the sanity checks)
    bool gens_over_Q = false, gens_over_K = false;
    // End_G = Q: generators of the image, complete by construction, with
    // p^e A^{-1} modulo p^W for the Newton steps
    if (!opt && OmegaX.empty()) {
      int Wn = MultiplicityNewtonPrecision(p);
      Tint pW = UniversalScalarConversion<Tint, int64_t>(1);
      for (int i = 0; i < Wn; i++) {
        pW *= UniversalScalarConversion<Tint, int64_t>(p);
      }
      MyMatrix<T> Binv = Inverse(GramX);
      T pe(1);
      for (int i = 0; i < e; i++) {
        pe *= UniversalScalarConversion<T, int64_t>(p);
      }
      ModMat Be(m, std::vector<int64_t>(m));
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
          T x = Binv(i, j) * pe;
          Tint den = GetDenominator_z(x);
          Tint inv = mod_inv(ResInt(den, pW), pW);
          Tint v = ResInt(Tint(GetNumerator_z(x) * inv), pW);
          if (v < 0) {
            v += pW;
          }
          Be[i][j] = UniversalScalarConversion<int64_t, Tint>(v);
        }
      }
      LocalImageOrder order;
      opt = LocalOrthogonalImageGenerators(A_mod, Be, {}, {}, {}, p, c, e, bounds,
                                           &order);
      gens_over_Q = opt.has_value();
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
      os << "MULT_FACTOR: p=" << p << " generators of the image "
         << (opt ? static_cast<long>(opt->size()) : -1L) << ", |image|="
         << LocalImageOrderValue(order, p) << "\n";
#endif
    }
    // a centre K beyond the pi-adic search: the same generators within the
    // integral commutant of O_K (K-linear maps), on the columns
    if (!opt && !OmegaX.empty()) {
      int Wn = MultiplicityNewtonPrecision(p);
      Tint pW = UniversalScalarConversion<Tint, int64_t>(1);
      for (int i = 0; i < Wn; i++) {
        pW *= UniversalScalarConversion<Tint, int64_t>(p);
      }
      MyMatrix<T> Binv = Inverse(GramX);
      T pe(1);
      for (int i = 0; i < e; i++) {
        pe *= UniversalScalarConversion<T, int64_t>(p);
      }
      ModMat Be(m, std::vector<int64_t>(m));
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
          T x = Binv(i, j) * pe;
          Tint den = GetDenominator_z(x);
          Tint inv = mod_inv(ResInt(den, pW), pW);
          Tint v = ResInt(Tint(GetNumerator_z(x) * inv), pW);
          if (v < 0) {
            v += pW;
          }
          Be[i][j] = UniversalScalarConversion<int64_t, Tint>(v);
        }
      }
      auto to_cols = [&](MyMatrix<T> const &Y) -> ModMat {
        ModMat Ym(m, std::vector<int64_t>(m));
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++) {
            Ym[i][j] = UniversalScalarConversion<int64_t, T>(Y(j, i));
          }
        }
        return Ym;
      };
      std::vector<ModMat> Om, Lam;
      for (auto const &W2 : OmegaX) {
        Om.push_back(to_cols(W2));
      }
      IsometryComponent<T> compX = BuildIsometryComponent<T, Tint>(
          A, Xb, IdentityMat<T>(A.rows()), 1);
      compX.extra_sigma = OmegaX;
      for (auto const &Y : IntegralCommutant<T, Tint>(compX)) {
        Lam.push_back(to_cols(UniversalMatrixConversion<T, Tint>(Y)));
      }
      LocalImageOrder order;
      opt = LocalOrthogonalImageGenerators(A_mod, Be, Om, Lam, {}, p, c, e,
                                           bounds, &order);
      if (opt) {
        rows_convention = false;
        gens_over_K = true;
      }
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
      os << "MULT_FACTOR: p=" << p << " c=" << c << " e=" << e
         << " K-linear generators of the image "
         << (opt ? static_cast<long>(opt->size()) : -1L) << ", |image|="
         << LocalImageOrderValue(order, p) << "\n";
#endif
    }
    if (!opt) {
      if (std::optional<size_t> fb = bounds_factor(size_t(1) << 12, 100000)) {
        factor *= UniversalScalarConversion<T, int64_t>(static_cast<int64_t>(*fb));
        continue;
      }
    }
    if (opt) {
      // found above
    } else if (OmegaX.empty()) {
      opt = LocalOrthogonalImage(A_mod, p, c, e, bounds);
    } else {
      // the K-linear endomorphisms of X, exact
      IsometryComponent<T> compX = BuildIsometryComponent<T, Tint>(
          A, Xb, IdentityMat<T>(A.rows()), 1);
      compX.extra_sigma = OmegaX;
      std::vector<ModMat> basis;
      for (auto const &Y : IntegralCommutant<T, Tint>(compX)) {
        ModMat Ym(m, std::vector<int64_t>(m));
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++) {
            Ym[i][j] = UniversalScalarConversion<int64_t, Tint>(Y(i, j));
          }
        }
        basis.push_back(Ym);
      }
      opt = LocalStabilizerImage(A_mod, basis, p, c, e, bounds);
    }
    if (!opt) {
      // the search is beyond its bounds: the bounds again, with the p-part
      // of D(N) allowed to be larger
      if (std::optional<size_t> fb = bounds_factor(size_t(1) << 17, 2000000)) {
        factor *= UniversalScalarConversion<T, int64_t>(static_cast<int64_t>(*fb));
        continue;
      }
      MULT_FACTOR_FAIL(2128);
      return {};
    }
#ifdef SANITY_CHECK_MULT
    // the K-linear generators over Z_p against the pi-adic image (rows,
    // a list or generators): the generated groups must agree
    if (opt && !gens_over_K && rows_convention && !OmegaX.empty()) {
      int Wn = MultiplicityNewtonPrecision(p);
      Tint pW = UniversalScalarConversion<Tint, int64_t>(1);
      for (int i = 0; i < Wn; i++) {
        pW *= UniversalScalarConversion<Tint, int64_t>(p);
      }
      MyMatrix<T> Binv = Inverse(GramX);
      T pe(1);
      for (int i = 0; i < e; i++) {
        pe *= UniversalScalarConversion<T, int64_t>(p);
      }
      ModMat Be(m, std::vector<int64_t>(m));
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
          T x = Binv(i, j) * pe;
          Tint den = GetDenominator_z(x);
          Tint inv = mod_inv(ResInt(den, pW), pW);
          Tint v = ResInt(Tint(GetNumerator_z(x) * inv), pW);
          if (v < 0) {
            v += pW;
          }
          Be[i][j] = UniversalScalarConversion<int64_t, Tint>(v);
        }
      }
      auto to_cols = [&](MyMatrix<T> const &Y) -> ModMat {
        ModMat Ym(m, std::vector<int64_t>(m));
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++) {
            Ym[i][j] = UniversalScalarConversion<int64_t, T>(Y(j, i));
          }
        }
        return Ym;
      };
      std::vector<ModMat> Om, Lam;
      for (auto const &W2 : OmegaX) {
        Om.push_back(to_cols(W2));
      }
      IsometryComponent<T> compX3 = BuildIsometryComponent<T, Tint>(
          A, Xb, IdentityMat<T>(A.rows()), 1);
      compX3.extra_sigma = OmegaX;
      for (auto const &Y : IntegralCommutant<T, Tint>(compX3)) {
        Lam.push_back(to_cols(UniversalMatrixConversion<T, Tint>(Y)));
      }
      LocalImageOrder order_k;
      std::optional<std::vector<ModMat>> kg = LocalOrthogonalImageGenerators(
          A_mod, Be, Om, Lam, {}, p, c, e, bounds, &order_k);
      if (kg) {
        int64_t pcl = 1;
        for (int i = 0; i < c; i++) {
          pcl *= p;
        }
        auto closure_of = [&](std::vector<ModMat> const &gens, bool transpose,
                              size_t cap) -> std::set<ModMat> {
          ModMat Idm(m, std::vector<int64_t>(m, 0));
          for (int i = 0; i < m; i++) {
            Idm[i][i] = 1;
          }
          std::vector<ModMat> gs;
          for (auto const &g : gens) {
            ModMat h(m, std::vector<int64_t>(m));
            for (int i = 0; i < m; i++) {
              for (int j = 0; j < m; j++) {
                h[i][j] = ModPos64(transpose ? g[j][i] : g[i][j], pcl);
              }
            }
            gs.push_back(h);
          }
          std::set<ModMat> cl{Idm};
          std::vector<ModMat> queue{Idm};
          while (!queue.empty() && cl.size() <= cap) {
            ModMat x = queue.back();
            queue.pop_back();
            for (auto const &g : gs) {
              ModMat y(m, std::vector<int64_t>(m, 0));
              for (int a = 0; a < m; a++) {
                for (int k2 = 0; k2 < m; k2++) {
                  for (int b = 0; b < m; b++) {
                    y[a][b] = ModPos64(y[a][b] + x[a][k2] * g[k2][b], pcl);
                  }
                }
              }
              if (cl.insert(y).second) {
                queue.push_back(y);
              }
            }
          }
          return cl;
        };
        // when the pi-adic result is the list of the image: the K-linear
        // generators are in it and their group, of the order given by the
        // chain, has its order; else the closures are compared
        std::set<ModMat> listed_set;
        for (auto const &g : *opt) {
          ModMat h(m, std::vector<int64_t>(m));
          for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
              h[i][j] = ModPos64(g[j][i], pcl);
            }
          }
          listed_set.insert(h);
        }
        bool all_in = true;
        for (auto const &g : *kg) {
          ModMat h = g;
          for (auto &row : h) {
            for (auto &x : row) {
              x = ModPos64(x, pcl);
            }
          }
          if (listed_set.count(h) == 0) {
            all_in = false;
          }
        }
        if (all_in &&
            static_cast<size_t>(std::llround(LocalImageOrderValue(order_k, p))) == listed_set.size()) {
          std::cerr << "MULT_FACTOR: SANITY: the K-linear generators lie in the "
                    << "pi-adic image, of the same order " << listed_set.size()
                    << " (p=" << p << ")\n";
        } else {
        std::set<ModMat> c2 = closure_of(*opt, true, 2000000);
        std::set<ModMat> c1 = closure_of(*kg, false, c2.size());
        if (c1 != c2) {
          std::cerr << "MULT_FACTOR: SANITY: the K-linear generators give "
                    << c1.size() << " elements, the pi-adic image "
                    << c2.size() << " (p=" << p << ")\n";
          throw TerminalException{1};
        }
        std::cerr << "MULT_FACTOR: SANITY: the K-linear generators give the "
                  << c1.size() << " elements of the pi-adic image (p=" << p
                  << ")\n";
        }
      }
    }
    // the group generated by the generators is the image listed by the old
    // searches, when these fit in their bounds (on the columns)
    if (gens_over_Q || gens_over_K) {
      std::optional<std::vector<ModMat>> listed;
      if (gens_over_Q) {
        listed = LocalOrthogonalImage(A_mod, p, c, e, bounds);
      } else {
        IsometryComponent<T> compX2 = BuildIsometryComponent<T, Tint>(
            A, Xb, IdentityMat<T>(A.rows()), 1);
        compX2.extra_sigma = OmegaX;
        std::vector<ModMat> basis2;
        for (auto const &Y : IntegralCommutant<T, Tint>(compX2)) {
          ModMat Ym(m, std::vector<int64_t>(m));
          for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
              Ym[i][j] = UniversalScalarConversion<int64_t, Tint>(Y(i, j));
            }
          }
          basis2.push_back(Ym);
        }
        listed = LocalStabilizerImage(A_mod, basis2, p, c, e, bounds);
        if (listed) {
          // rows to columns
          for (auto &g : *listed) {
            ModMat h(m, std::vector<int64_t>(m));
            for (int i = 0; i < m; i++) {
              for (int j = 0; j < m; j++) {
                h[i][j] = g[j][i];
              }
            }
            g = h;
          }
        }
      }
      if (listed) {
        int64_t pcl = 1;
        for (int i = 0; i < c; i++) {
          pcl *= p;
        }
        std::set<ModMat> listed_set;
        for (auto const &g : *listed) {
          ModMat h = g;
          for (auto &row : h) {
            for (auto &x : row) {
              x = ModPos64(x, pcl);
            }
          }
          listed_set.insert(h);
        }
        ModMat Idm(m, std::vector<int64_t>(m, 0));
        for (int i = 0; i < m; i++) {
          Idm[i][i] = 1;
        }
        std::set<ModMat> closure{Idm};
        std::vector<ModMat> queue{Idm};
        while (!queue.empty() && closure.size() <= listed_set.size()) {
          ModMat x = queue.back();
          queue.pop_back();
          for (auto const &g : *opt) {
            ModMat y(m, std::vector<int64_t>(m, 0));
            for (int a = 0; a < m; a++) {
              for (int k2 = 0; k2 < m; k2++) {
                for (int b = 0; b < m; b++) {
                  y[a][b] = ModPos64(y[a][b] + x[a][k2] * ModPos64(g[k2][b], pcl), pcl);
                }
              }
            }
            if (closure.insert(y).second) {
              queue.push_back(y);
            }
          }
        }
        if (closure != listed_set) {
          std::cerr << "MULT_FACTOR: SANITY: the generators give "
                    << closure.size() << " elements, the list has "
                    << listed_set.size() << " (p=" << p << ")\n";
          throw TerminalException{1};
        }
        std::cerr << "MULT_FACTOR: SANITY: the generators give the "
                  << listed_set.size() << " elements of the list (p=" << p
                  << ")\n";
      }
    }
#endif
    // the action of an element on the generators f_s of the p-part: the
    // rows f_s (B_V coordinates) moved by 1 (x) g
    auto act_gens = [&](ModMat const &g) -> std::vector<std::vector<int64_t>> {
      std::vector<std::vector<int64_t>> res;
      for (size_t s = 0; s < r; s++) {
        MyVector<T> f = GetMatrixRow(F_c, ListIdx[s]);
        MyVector<T> img = ZeroVector<T>(k);
        for (int a = 0; a < w; a++) {
          for (int i = 0; i < m; i++) {
            T x = f(a * m + i);
            if (x == T(0)) {
              continue;
            }
            for (int j = 0; j < m; j++) {
              // rows convention: x_i -> sum_j g(j, i) x_j, g acting on the
              // columns of coordinates of X
              int64_t gij = rows_convention ? g[i][j] : g[j][i];
              img(a * m + j) += x * UniversalScalarConversion<T, int64_t>(gij);
            }
          }
        }
        res.push_back(residue(img));
      }
      return res;
    };
    std::vector<std::vector<std::vector<int64_t>>> ListAct;
    for (auto const &g : *opt) {
      ListAct.push_back(act_gens(g));
    }
    auto apply = [&](std::vector<std::vector<int64_t>> const &act,
                     std::vector<int64_t> const &x) -> std::vector<int64_t> {
      std::vector<int64_t> y(r, 0);
      for (size_t s = 0; s < r; s++) {
        if (x[s] == 0) {
          continue;
        }
        for (size_t u = 0; u < r; u++) {
          y[u] = (y[u] + MulModQ(x[s], act[s][u], ListOrd[u])) % ListOrd[u];
        }
      }
      return y;
    };
    auto decode = [&](size_t idx) -> std::vector<int64_t> {
      std::vector<int64_t> x(r);
      for (size_t s = 0; s < r; s++) {
        x[s] = static_cast<int64_t>(idx % static_cast<size_t>(ListOrd[s]));
        idx /= static_cast<size_t>(ListOrd[s]);
      }
      return x;
    };
    // the orbit of L/N, each subgroup by its Howell normal form in (Z /
    // p^C)^r, coordinate s scaled by p^{C - a_s} (ListOrd[s] = p^{a_s})
    int Cexp = 0;
    std::vector<int> aexp(r, 0);
    for (size_t s2 = 0; s2 < r; s2++) {
      for (int64_t o = ListOrd[s2]; o > 1; o /= p) {
        aexp[s2]++;
      }
      Cexp = std::max(Cexp, aexp[s2]);
    }
    auto powp = [&](int k) -> int64_t {
      int64_t x = 1;
      for (int i = 0; i < k; i++) {
        x *= p;
      }
      return x;
    };
    auto canon = [&](std::vector<std::vector<int64_t>> const &gens_orig)
        -> ModMat {
      ModMat sc;
      for (auto const &x : gens_orig) {
        std::vector<int64_t> y(r);
        for (size_t s2 = 0; s2 < r; s2++) {
          y[s2] = x[s2] * powp(Cexp - aexp[s2]);
        }
        sc.push_back(y);
      }
      return HowellFormPrimePower(sc, p, Cexp);
    };
    auto unscale = [&](std::vector<int64_t> const &y) -> std::vector<int64_t> {
      std::vector<int64_t> x(r);
      for (size_t s2 = 0; s2 < r; s2++) {
        x[s2] = (y[s2] / powp(Cexp - aexp[s2])) % ListOrd[s2];
      }
      return x;
    };
    // the element lists only serve the sanity check
    (void)span;
    (void)decode;
    ModMat L0c = canon(Lgens);
    std::set<ModMat> orbit{L0c};
    std::vector<ModMat> queue{L0c};
    while (!queue.empty()) {
      ModMat H = queue.back();
      queue.pop_back();
      for (auto const &act : ListAct) {
        std::vector<std::vector<int64_t>> img;
        for (auto const &row : H) {
          img.push_back(apply(act, unscale(row)));
        }
        ModMat Hc = canon(img);
        if (orbit.insert(Hc).second) {
          queue.push_back(Hc);
        }
      }
    }
    // |L/N| from the form: the product of the p^{C - v} of the pivots
    double logL = 0;
    for (auto const &row : L0c) {
      for (size_t s2 = 0; s2 < r; s2++) {
        if (row[s2] != 0) {
          int v = 0;
          for (int64_t x = row[s2]; x % p == 0; x /= p) {
            v++;
          }
          logL += (Cexp - v) * std::log(static_cast<double>(p));
          break;
        }
      }
    }
#ifdef SANITY_CHECK_MULT
    // the orbit on the sorted element lists, when L/N is small
    if (q_fits && std::exp(logL) <= 200000.5) {
      std::vector<size_t> L0 = span(Lgens);
      std::set<std::vector<size_t>> orbit_old{L0};
      std::vector<std::vector<size_t>> queue_old{L0};
      while (!queue_old.empty()) {
        std::vector<size_t> H = queue_old.back();
        queue_old.pop_back();
        for (auto const &act : ListAct) {
          std::vector<size_t> img;
          for (size_t idx : H) {
            img.push_back(encode(apply(act, decode(idx))));
          }
          std::sort(img.begin(), img.end());
          if (orbit_old.insert(img).second) {
            queue_old.push_back(img);
          }
        }
      }
      if (orbit_old.size() != orbit.size() ||
          static_cast<double>(L0.size()) != std::round(std::exp(logL))) {
        std::cerr << "MULT_FACTOR: SANITY: the Howell orbit " << orbit.size()
                  << " (|L/N| " << std::exp(logL) << ") against "
                  << orbit_old.size() << " (|L/N| " << L0.size() << ")\n";
        throw TerminalException{1};
      }
      os << "MULT_FACTOR: SANITY: the Howell orbit agrees, " << orbit.size()
         << "\n";
    }
#endif
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    os << "MULT_FACTOR: p=" << p << " c=" << c << " e=" << e
       << " |image of O(X_p)|=" << opt->size() << " |N^#/N|_p"
       << (q_fits ? "=" : ">=") << q_size
       << " |L/N|_p=" << std::exp(logL) << " orbit=" << orbit.size() << "\n";
#endif
    factor *= UniversalScalarConversion<T, int64_t>(
        static_cast<int64_t>(orbit.size()));
  }
  return factor;
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
    // the class sums, the exponent N of G and the a with sigma_a fixing
    // the centre K of End_G (the subfield of Q(zeta_N) it is)
    auto centre_field = [&]()
        -> std::tuple<std::vector<MyMatrix<T>>, int64_t, std::vector<int64_t>> {
      // the conjugacy classes and their sums on the block
      std::vector<int> ccls(n_elt, -1);
      int n_ccls = 0;
      for (size_t i = 0; i < n_elt; i++) {
        if (ccls[i] != -1) {
          continue;
        }
        for (size_t j = 0; j < n_elt; j++) {
          size_t u = index[flatten(MyMatrix<T>(inv[j] * elts[i] * elts[j]))];
          ccls[u] = n_ccls;
        }
        n_ccls++;
      }
      std::vector<MyMatrix<T>> ClassSum(n_ccls, ZeroMatrix<T>(n, n));
      for (size_t i = 0; i < n_elt; i++) {
        ClassSum[ccls[i]] += elts[i];
      }
      int64_t Nexp = 1;
      for (auto const &x : elts) {
        Nexp = std::lcm(Nexp, static_cast<int64_t>(IsometryOrder(x)));
      }
      std::vector<int64_t> Hgens;
      for (int64_t a = 1; a < Nexp; a++) {
        if (std::gcd(a, Nexp) != 1) {
          continue;
        }
        bool fixes = true;
        for (size_t i = 0; i < n_elt && fixes; i++) {
          MyMatrix<T> pw = IdentityMat<T>(n);
          for (int64_t t = 0; t < a; t++) {
            pw = pw * elts[i];
          }
          int ci = ccls[i];
          int ca = ccls[index[flatten(pw)]];
          if (ci != ca && action(ClassSum[ci]) != action(ClassSum[ca])) {
            fixes = false;
          }
        }
        if (fixes) {
          Hgens.push_back(a);
        }
      }
      return {ClassSum, Nexp, Hgens};
    };
    std::vector<MyMatrix<T>> End = commutant(ListAct);
    // the centre of End_G: the elements of End_G commuting with End_G
    std::vector<MyMatrix<T>> ListBoth = ListAct;
    ListBoth.insert(ListBoth.end(), End.begin(), End.end());
    int c = static_cast<int>(commutant(ListBoth).size());
    int e = static_cast<int>(End.size());
    int m = static_cast<int>(std::llround(std::sqrt(static_cast<double>(e))));
    bool central_q = (c == 1);
    {
      // End_G = M_r(D), D a quaternion algebra (Albert type III), in the
      // columns convention
      std::vector<MyMatrix<T>> EndC, CentreC, GensC, EltsC;
      for (auto const &Y : End) {
        EndC.push_back(Y.transpose());
      }
      for (auto const &Z : commutant(ListBoth)) {
        CentreC.push_back(Z.transpose());
      }
      int rq = QuaternionicRank(Gram, EndC, CentreC);
      if (rq > 0) {
        for (auto const &M : ListAct) {
          GensC.push_back(M.transpose());
        }
        for (auto const &x : elts) {
          EltsC.push_back(action(x).transpose());
        }
        for (auto const &C : GensC) {
          if (MyMatrix<T>(C.transpose() * Gram * C) != Gram) {
            std::cerr << "ISOM_MASS: a generator is not an isometry of the "
                         "block\n";
            throw TerminalException{1};
          }
        }
        std::optional<AbelianCMField> FKq;
        T dKq(1);
        if (c > 1) {
          auto [ClassSumq, Nexpq, Hgensq] = centre_field();
          (void)ClassSumq;
          FKq = MakeAbelianRealField(Nexpq, Hgensq);
          if (AbelianFieldDegree(*FKq) != c) {
            std::cerr << "ISOM_MASS: the centre has degree " << c
                      << " but the field found has degree "
                      << AbelianFieldDegree(*FKq) << "\n";
            throw TerminalException{1};
          }
          DirichletGroup Gd = ComputeDirichletGroup(FKq->n);
          for (auto const &t : FieldCharacters(Gd, FKq->H, false)) {
            dKq *= ScalarFromInt64<T>(PrimitiveCharacter(Gd, t).first);
          }
        }
        comp.quaternionic = true;
        QuaternionicMassBounds qbounds;
        comp.quaternionic_mass = QuaternionicComponentMass<T, Tint>(
            Gram, EltsC, GensC, CentreC, FKq, dKq, rq, qbounds, os);
#ifdef DEBUG_ISOMETRY_MASS
        os << "ISOM_MASS: quaternionic component of rank " << k << ", r=" << rq
           << ", centre of degree " << c << ", mass "
           << (comp.quaternionic_mass ? "found" : "not found") << "\n";
#endif
        res.push_back(comp);
        continue;
      }
    }
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
            A, found_basis, blk.E, elts, ListP, {}, {}, mbounds, os);
      }
    }
    // with a centre of degree c, End_G = M_m(K) has dimension m^2 c
    if (c > 1 && e % c == 0) {
      m = static_cast<int>(std::llround(std::sqrt(static_cast<double>(e / c))));
    }
    if (c > 1 && m * m * c == e) {
      /*
        End_G = M_m(K) with K = Q(chi) the centre, of degree c. The
        multiplicity space is X = L cap F, F the image of an idempotent of
        Q[H] of K-rank one on W, of dimension m c, with the action of the
        order generated by the class sums: Hermitian over K when K is CM
        (the adjoint acts nontrivially on the centre), quadratic over K
        when K is totally real. K is the subfield of Q(zeta_N), N the
        exponent of G, fixed by the a with C_g = C_{g^a} on V_rho for every
        class (sigma_a sends the central character of C_g to that of
        C_{g^a}). F is stable under End_G, so it has no dimension m c when
        End_G is not split (M_{m/2}(D), D a quaternion algebra over K, has
        no module of that dimension), and the component is then left
        without a mass.
       */
      std::vector<MyMatrix<T>> Centre = commutant(ListBoth);
      MyMatrix<T> Gblk = blk.E * A * blk.E.transpose();
      MyMatrix<T> GblkInv = Inverse(Gblk);
      bool is_cm = false;
      for (auto const &Z : Centre) {
        MyMatrix<T> Zs = (GblkInv * Z * Gblk).transpose();
        if (Zs != Z) {
          is_cm = true;
        }
      }
      auto [ClassSum, Nexp, Hgens] = centre_field();
      std::optional<MyMatrix<T>> found_basis;
      if (Nexp > 2) {
        std::vector<std::pair<size_t, int>> cand;
        for (size_t i = 0; i < n_elt; i++) {
          cand.push_back({i, 1});
          cand.push_back({i, -1});
        }
        auto try_conditions =
            [&](std::vector<std::pair<size_t, int>> const &cs) -> bool {
          std::vector<MyMatrix<T>> all = blk.conditions;
          for (auto const &cc : cs) {
            all.push_back(elts[cc.first] -
                          ScalarFromInt64<T>(cc.second) * IdentityMat<T>(n));
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
          if (K_i.rows() != m * c) {
            return false;
          }
          found_basis = UniversalMatrixConversion<T, Tint>(K_i);
          return true;
        };
        for (size_t a = 0; a < cand.size() && !found_basis; a++) {
          (void)try_conditions({cand[a]});
        }
        for (size_t a = 0; a < cand.size() && !found_basis; a++) {
          for (size_t b = a + 1; b < cand.size() && !found_basis; b++) {
            if (cand[a].first != cand[b].first) {
              (void)try_conditions({cand[a], cand[b]});
            }
          }
        }
      }
      if (found_basis) {
        MyMatrix<T> const &Xb = *found_basis;
        int mX = Xb.rows();
        MyMatrix<T> XX = Xb * Xb.transpose();
        MyMatrix<T> XXinv = Inverse(XX);
        // the class sums on X, then the ring they generate
        auto on_X = [&](MyMatrix<T> const &Mc) -> MyMatrix<T> {
          return Xb * Mc * Xb.transpose() * XXinv;
        };
        std::vector<MyMatrix<T>> gens{IdentityMat<T>(mX)};
        for (auto const &Cs : ClassSum) {
          gens.push_back(on_X(Cs));
        }
        auto zbasis = [&](std::vector<MyMatrix<T>> const &L)
            -> std::vector<MyMatrix<T>> {
          MyMatrix<Tint> F(static_cast<int>(L.size()), mX * mX);
          for (size_t i = 0; i < L.size(); i++) {
            for (int u = 0; u < mX; u++) {
              for (int v = 0; v < mX; v++) {
                F(static_cast<int>(i), u * mX + v) =
                    UniversalScalarConversion<Tint, T>(L[i](u, v));
              }
            }
          }
          MyMatrix<Tint> Bz = GetZbasis(F);
          std::vector<MyMatrix<T>> res_b;
          for (int i = 0; i < Bz.rows(); i++) {
            MyMatrix<T> Y(mX, mX);
            for (int u = 0; u < mX; u++) {
              for (int v = 0; v < mX; v++) {
                Y(u, v) = UniversalScalarConversion<T, Tint>(Bz(i, u * mX + v));
              }
            }
            res_b.push_back(Y);
          }
          return res_b;
        };
        std::vector<MyMatrix<T>> Omega = zbasis(gens);
        for (int iter = 0; iter < 8; iter++) {
          std::vector<MyMatrix<T>> prods = Omega;
          for (auto const &X1 : Omega) {
            for (auto const &X2 : Omega) {
              prods.push_back(X1 * X2);
            }
          }
          std::vector<MyMatrix<T>> next = zbasis(prods);
          if (next.size() == Omega.size()) {
            Omega = next;
            break;
          }
          Omega = next;
        }
        AbelianCMField Fk = is_cm ? MakeAbelianCMField(Nexp, Hgens)
                                  : MakeAbelianRealField(Nexp, Hgens);
        // the order is O_K when its discriminant is that of K (the product
        // of the conductors of the characters)
        bool maximal = (static_cast<int>(Omega.size()) == c);
        if (maximal) {
          MyMatrix<T> Tr(c, c);
          for (int i = 0; i < c; i++) {
            for (int j = 0; j < c; j++) {
              MyMatrix<T> P = Omega[i] * Omega[j];
              T t(0);
              for (int u = 0; u < mX; u++) {
                t += P(u, u);
              }
              Tr(i, j) = t / ScalarFromInt64<T>(m);
            }
          }
          T dO = T_abs(DeterminantMat(Tr));
          DirichletGroup Gd = ComputeDirichletGroup(Fk.n);
          T dK(1);
          for (auto const &t : FieldCharacters(Gd, Fk.H, false)) {
            dK *= ScalarFromInt64<T>(PrimitiveCharacter(Gd, t).first);
          }
          maximal = (dO == dK);
        }
        if (maximal) {
          if (is_cm) {
            comp.herm_field = Fk;
          } else {
            comp.quad_field = Fk;
          }
          comp.centre_gram = Xb * A * Xb.transpose();
          comp.centre_omega = Omega;
          MultiplicityFactorBounds mbounds;
          comp.multiplicity_factor = MultiplicityLocalFactor<T, Tint>(
              A, Xb, blk.E, elts, ListP, Omega, ClassSum, mbounds, os);
        }
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
       << ", orbit:=" << gl.orbit_size << ", method:=\"" << gl.method
       << "\", image:=[";
    for (size_t s = 0; s < gl.list_image.size(); s++) {
      os << (s > 0 ? ", " : "") << "\"" << gl.list_image[s] << "\"";
    }
    os << "])";
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

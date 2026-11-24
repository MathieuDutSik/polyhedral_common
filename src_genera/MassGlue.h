// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_MASSGLUE_H_
#define SRC_GENERA_MASSGLUE_H_

// clang-format off
#include "MassComponents.h"
#include "MassLocalImage.h"
// clang-format on

/*
  The glue of the components (MassIsometry.h): the q-primary parts of the
  discriminant forms of the components as modules (PrimaryQuadModule),
  the actions of the isometries on them, the images of the local
  centralizers (integral commutant, Cayley transforms, local
  reflections), and the orbits of the glue subgroups under them.
 */

/*
  The q-primary part of D(L_d): generators g_i (one per elementary divisor
  divisible by q) of orders ord_i, N the exponent, b(g_i, g_j) = bnum_ij / N
  mod 1, q(g_i) = qnum_i / N mod 2, and the images of the generators by
  sigma in coordinates. The coordinates are residues modulo ord_i, bnum
  modulo N and qnum modulo 2 N, all handled by ResidueArithmetic.h:
  BuildPrimaryQuadModule checks that 2 N is a valid modulus, after which
  no operation on the module can overflow.
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
    if (x[i] != 0) {
      for (size_t j = 0; j < M.rank(); j++) {
        int64_t t = ResMulI64(x[i], y[j], M.N);
        val = ResAddI64(val, ResMulI64(t, M.bnum[i][j], M.N), M.N);
      }
    }
  }
  return val;
}

inline int64_t ModuleQ(PrimaryQuadModule const &M, ModuleElt const &x) {
  int64_t N2 = 2 * M.N;
  int64_t val = 0;
  for (size_t i = 0; i < M.rank(); i++) {
    int64_t t = ResMulI64(x[i], x[i], N2);
    val = ResAddI64(val, ResMulI64(t, M.qnum[i], N2), N2);
    for (size_t j = i + 1; j < M.rank(); j++) {
      int64_t u = ResMulI64(2, ResMulI64(x[i], x[j], N2), N2);
      val = ResAddI64(val, ResMulI64(u, M.bnum[i][j], N2), N2);
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
    if (a[i] != 0) {
      for (size_t j = 0; j < M.rank(); j++) {
        res[j] = ResAddI64(res[j], ResMulI64(a[i], y[i][j], M.ord[j]), M.ord[j]);
      }
    }
  }
  return res;
}

/*
  The coordinates in M of the class of x in L_d^#: with z = x ROW^{-1},
  z_i s_i is an integer, and the coordinate on the generator of order
  ord_i = q^v (s_i = q^v t) is (z_i s_i) t^{-1} mod ord_i.
 */
template <typename T>
ModuleElt ModuleCoordinates(PrimaryQuadModule const &M,
                            IsometryComponent<T> const &comp,
                            MyVector<T> const &x) {
  MyVector<T> z = comp.row_inv.transpose() * x;
  ModuleElt a(M.rank());
  for (size_t i = 0; i < M.rank(); i++) {
    int ii = M.index[i];
    T val = z(ii) * UniversalScalarConversion<
                        T, typename underlying_z_ring<T>::ring_type>(
                        comp.elem[ii]);
    if (!IsInteger(val)) {
      std::cerr << "ISOM_MASS: a vector is not in the dual lattice\n";
      throw TerminalException{1};
    }
    int64_t r = ResModOfI64(GetNumerator_z(val), M.ord[i]);
    a[i] = ResMulI64(r, M.tinv[i], M.ord[i]);
  }
  return a;
}

template <typename T>
PrimaryQuadModule BuildPrimaryQuadModule(IsometryComponent<T> const &comp,
                                         int64_t q) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  PrimaryQuadModule M;
  M.q = q;
  M.d = comp.d;
  std::vector<MyVector<T>> gens;
  int k = comp.gram.rows();
  Tint q_T = UniversalScalarConversion<Tint, int64_t>(q);
  for (size_t i = 0; i < comp.elem.size(); i++) {
    // s = q^v t, t prime to q
    Tint t = comp.elem[i];
    int v = 0;
    while (ResInt(t, q_T) == Tint(0)) {
      t = QuoInt(t, q_T);
      v++;
    }
    if (v == 0) {
      continue;
    }
    int64_t qe = ResiduePrimePower(q, v);
    // the quadratic form is read modulo 2 N
    if (qe > ResToSigned(kMaxResidueModulus / 2)) {
      ResidueModulusTooLarge("2 * " + std::to_string(qe));
    }
    M.index.push_back(static_cast<int>(i));
    M.ord.push_back(qe);
    M.tinv.push_back(ResInvI64(ResModOfI64(t, qe), qe));
    M.N = std::max(M.N, qe);
    MyVector<T> g(k);
    for (int j = 0; j < k; j++) {
      g(j) = comp.row(i, j) / UniversalScalarConversion<T, int64_t>(qe);
    }
    gens.push_back(g);
  }
  size_t r = gens.size();
  M.bnum.assign(r, std::vector<int64_t>(r, 0));
  M.qnum.assign(r, 0);
  T N_T = UniversalScalarConversion<T, int64_t>(M.N);
  for (size_t i = 0; i < r; i++) {
    for (size_t j = 0; j < r; j++) {
      T val = N_T * gens[i].dot(comp.gram * gens[j]);
      if (!IsInteger(val)) {
        std::cerr << "ISOM_MASS: the discriminant form is not of exponent N\n";
        throw TerminalException{1};
      }
      M.bnum[i][j] = ResidueOfInteger(val, M.N);
      if (i == j) {
        M.qnum[i] = ResidueOfInteger(val, 2 * M.N);
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
      g(j) = comp.row(i, j) / UniversalScalarConversion<T, int64_t>(M.ord[s]);
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
      if (ResMulI64(ResModI64(o, M.ord[j]), x[j], M.ord[j]) != 0) {
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
  // the equivariance is only checked on the complete maps: a budget of
  // nodes, beyond which nothing is returned (the image of the local
  // centralizer then comes from its generators, LocalImageOnModule)
  size_t n_node = 0;
  size_t const max_node = 500000;
  std::function<void(size_t)> rec = [&](size_t i) -> void {
    if (overflow) {
      return;
    }
    if (++n_node > max_node) {
      overflow = true;
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
          y[pos] = ResAddI64(x[pos], g[s][i], radix[pos]);
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
  complete graph on the elements of M with the vertices colored by q(x)
  and the edges by b(x, y) and by whether y = sigma(x). The glue H is added
  as one vertex per element, joined to its components: the automorphisms of
  that graph are the stabilizer of H in G = prod_s G_s (the module of each
  vertex is in its color, so the modules are not exchanged, and an element
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
  // the images by the generators of the group, one edge color bit each
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
  // the edge colors are coded as 4 ((s BIG + b) 2^n_bit + flag) + 1 with
  // s < t, b < N <= BIG / 4, flag < 2^n_bit: the largest code must fit
  {
    __int128 code = static_cast<__int128>(t) * BIG + BIG;
    for (size_t b = 0; b < n_bit + 1 && code < (__int128(1) << 62); b++) {
      code *= 2;
    }
    code = 4 * code + 1;
    if (code >= (__int128(1) << 62)) {
      return {};
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
            y[pos] = ResAddI64(x[pos], g[s][i], radix[pos]);
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
  MyMatrix<Tint> Eq = ZeroMatrix<Tint>(kk, kk * static_cast<int>(ListX.size()));
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
  MyMatrix<Tint> ker = IntegralKernelOfConditions(Eq);
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
  Elements of the local centralizer by the Cayley transform: for Y in the
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
  Tint qc = MyPow<Tint>(q_T, c);
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
    int64_t qc64 = ResidueModulusOf(qc);
    int64_t dinv = ResInvI64(ResModOfI64(dm, qc64), qc64);
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
  the stabilizer of L_d in the centralizer C(Q_q) of the group acting on
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
          e = std::max(e, RationalValuation(GetDenominator(Ainv(i, j)), q));
        }
      }
    }
  }
  // the listing keeps its solutions by their coordinates in the basis of
  // Lambda (an integral kernel), so they stay in Lambda; it only sees the
  // basis and the Gram matrix modulo q^{K+1}, given as residues, whatever
  // the size of the exact entries. Beyond 2^62 the listing is left out.
  std::optional<int64_t> listing_mod = LocalImageModulus(q, c, e);
  bool listable = listing_mod.has_value();
  auto exact = [&](MyMatrix<Tint> const &X) -> ModRowsI64 {
    ModRowsI64 R(k, std::vector<int64_t>(k, 0));
    if (listable) {
      for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
          R[i][j] = CenteredResidueOfInteger(X(i, j), *listing_mod);
        }
      }
    }
    return R;
  };
  // the residues modulo q^W, for the certified generators
  int Wq = MultiplicityNewtonPrecision(M.q);
  Tint qW = MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(M.q), Wq);
  auto residues = [&](MyMatrix<Tint> const &X) -> ModRowsI64 {
    ModRowsI64 R(k, std::vector<int64_t>(k));
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < k; j++) {
        Tint v = ResInt(X(i, j), qW);
        if (v < 0) {
          v += qW;
        }
        R[i][j] = UniversalScalarConversion<int64_t, Tint>(v);
      }
    }
    return R;
  };
  std::vector<ModRowsI64> basis;
  for (int r = 0; r < ker.rows(); r++) {
    MyMatrix<Tint> Y(k, k);
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        Y(u, v) = ker(r, u * k + v);
      }
    }
    basis.push_back(exact(Y));
  }
  ModRowsI64 A_mod = exact(UniversalMatrixConversion<Tint, T>(comp.gram));
  MultiplicityFactorBounds mbounds;
  /*
    First the certified generators of the image modulo q^c
    (LocalOrthogonalImageGenerators, on the columns, commuting with the
    group and within its integral commutant), and the closure of their
    actions on M: the image is not listed. The listing of the image
    (LocalStabilizerImage) when they are beyond the bounds.
   */
  {
    int Wn = MultiplicityNewtonPrecision(q);
    Tint pW = MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(q), Wn);
    MyMatrix<T> Binv = Inverse(comp.gram);
    T pe = MyPow<T>(UniversalScalarConversion<T, int64_t>(q), e);
    ModRowsI64 Be(k, std::vector<int64_t>(k));
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < k; j++) {
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
    // the group on the columns: the transposes when the matrices act on
    // the rows (sigma A sigma^T = A)
    std::vector<MyMatrix<T>> ListX{comp.sigma};
    ListX.insert(ListX.end(), comp.extra_sigma.begin(), comp.extra_sigma.end());
    bool on_rows = true;
    for (auto const &X : ListX) {
      if (MyMatrix<T>(X * comp.gram * X.transpose()) != comp.gram) {
        on_rows = false;
      }
    }
    std::vector<ModRowsI64> Om, Lam;
    for (auto const &X : ListX) {
      MyMatrix<Tint> Xi = UniversalMatrixConversion<Tint, T>(
          on_rows ? MyMatrix<T>(X.transpose()) : X);
      Om.push_back(residues(Xi));
    }
    for (auto const &Y : ListComm) {
      Lam.push_back(residues(on_rows ? MyMatrix<Tint>(Y.transpose()) : Y));
    }
    std::vector<MyMatrix<Tint>> CommCols;
    for (auto const &Y : ListComm) {
      CommCols.push_back(on_rows ? MyMatrix<Tint>(Y.transpose()) : Y);
    }
    std::optional<std::vector<ModRowsI64>> Lie = ResidueMatrices(
        IntegralLieLatticeColumns<T, Tint>(comp.gram, CommCols),
        ResidueModulusOf(qW));
    LocalImageOrder ord;
    std::optional<std::vector<ModRowsI64>> ggens = LocalOrthogonalImageGenerators(
        residues(UniversalMatrixConversion<Tint, T>(comp.gram)), Be, Om, Lam,
        {}, Lie, q, c, e, mbounds, &ord);
    if (ggens) {
      std::vector<MyVector<T>> mgens = PrimaryModuleGenerators(comp, M);
      std::vector<ModuleAut> gimg;
      for (auto const &g : *ggens) {
        // back to the convention of ModuleActionOfIsometry (the rows)
        MyMatrix<Tint> h(k, k);
        for (int i = 0; i < k; i++) {
          for (int j = 0; j < k; j++) {
            h(i, j) = UniversalScalarConversion<Tint, int64_t>(
                on_rows ? g[j][i] : g[i][j]);
          }
        }
        gimg.push_back(ModuleActionOfIsometry(comp, M, mgens, h));
      }
      size_t r = M.rank();
      ModuleAut id(r, ModuleElt(r, 0));
      for (size_t i = 0; i < r; i++) {
        id[i][i] = 1;
      }
      std::set<ModuleAut> seen{id};
      std::vector<ModuleAut> queue{id};
      bool too_big = false;
      while (!queue.empty() && !too_big) {
        ModuleAut x = queue.back();
        queue.pop_back();
        for (auto const &gg : gimg) {
          ModuleAut y;
          for (auto const &xi : x) {
            y.push_back(ApplyModuleAut(M, gg, xi));
          }
          if (seen.insert(y).second) {
            queue.push_back(y);
            if (seen.size() > 1000000) {
              too_big = true;
              break;
            }
          }
        }
      }
      if (!too_big) {
#ifdef DEBUG_ISOMETRY_MASS
        os << "ISOM_MASS: q=" << q << " d=" << comp.d << " image on D="
           << seen.size() << " from " << ggens->size() << " generators\n";
#endif
#ifdef SANITY_CHECK_ISOMETRY_MASS
        // the listed image, when within the bounds, is the same set
        std::optional<std::vector<ModRowsI64>> lst;
        if (listable) {
          lst = LocalStabilizerImage(A_mod, basis, q, c, e, mbounds);
        }
        if (lst) {
          std::set<ModuleAut> listed;
          for (auto const &g : *lst) {
            MyMatrix<Tint> h(k, k);
            for (int i = 0; i < k; i++) {
              for (int j = 0; j < k; j++) {
                h(i, j) = UniversalScalarConversion<Tint, int64_t>(g[i][j]);
              }
            }
            listed.insert(ModuleActionOfIsometry(comp, M, mgens, h));
          }
          if (listed != seen) {
            std::cerr << "ISOM_MASS: the image on D from the generators ("
                      << seen.size() << ") differs from the listed one ("
                      << listed.size() << ")\n";
            throw TerminalException{1};
          }
        }
#endif
        return std::vector<ModuleAut>(seen.begin(), seen.end());
      }
    }
  }
  if (!listable) {
    return {};
  }
  std::optional<std::vector<ModRowsI64>> opt =
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
  Tint pc = MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(p), c);
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
  size_t n_cand = MyPow<size_t>(3, m);
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
          z[i] = ResAddI64(x[i], y[i], M.ord[i]);
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

// clang-format off
#endif  // SRC_GENERA_MASSGLUE_H_
// clang-format on

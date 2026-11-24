// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_MASSISOMETRY_H_
#define SRC_GENERA_MASSISOMETRY_H_

// clang-format off
#include "LocalDensity.h"
#include "MassGenus.h"
#include "MassGroupComponents.h"
// clang-format on

#ifdef SANITY_CHECK
# define SANITY_CHECK_ISOMETRY_MASS
#endif

/*
  The mass of the genus of a lattice with an isometry of finite order,

      mass(L, sigma) = sum_{(L', sigma)} 1 / |C_{O(L')}(sigma)|,

  the sum over the isometry classes of pairs in the genus of the pair (the
  sigma-stable lattices of V = L (x) Q locally isometric to L as pairs).

  Write C for the centralizer of sigma in O(V), K the stabilizer of L in
  C(A_f), L_d = L cap ker Phi_d(sigma) for the d dividing the order m, and
  L_0 = (+)_d L_d, of finite index in L. An element of C preserves every
  ker Phi_d(sigma), so K is contained in the stabilizer K_0 = prod_d K_d of
  L_0, and since C = prod_d C_d (orthogonal groups for d = 1, 2, unitary
  groups over Q(zeta_d) for d >= 3):

      mass(L, sigma) = prod_d mass(L_d, sigma_d) * [K_0 : K].

  The index is local: [K_0 : K] = prod_q |G_q . H_q|, the orbit of the glue
  H_q = (L / L_0) (x) Z_q in D(L_0)_q = (+)_d D(L_d)_q under the image G_q
  of K_{0,q}. H is a subgroup of the discriminant group and K_q is its
  stabilizer. Only the primes q dividing [L : L_0] contribute, and they
  divide m.

  The image G_q is the product over d of the images of the local groups
  of the L_d in the isometries of the discriminant form D(L_d)_q that
  commute with sigma. For d = 1, 2 (sigma = +-1 on L_d) O(L_{d,q}) maps
  onto O(D(L_d)_q) (Nikulin, Integral symmetric bilinear forms and some of
  their applications, Theorem 1.9.5), so the full group is the image. For
  d >= 3 the image of the local unitary group is computed (or proved to be
  all the isometries commuting with sigma by elements of it); that it is
  always all of them is not proved (Question q:surj of the paper), and it
  is assumed only when IsometryMassBounds::allow_image_assumption is set,
  the result then saying image_assumed.

  The masses of the orthogonal L_d come from GenusMass. The genus of the
  pair (L_d, sigma_d), d >= 3, is the genus of L_d as a Hermitian lattice
  over Z[zeta_d], for the form h with b = Tr h. Its mass is computed
  exactly (MassCyclotomicHermitian.h); a mass supplied by the caller takes
  precedence.

  Conventions: the Gram matrix A of L, the isometry g acting on the column
  vectors, g^T A g = A. Internally the row vectors are acted on by P =
  g^T, as in the rest of the code.
 */


template <typename T> struct IsometryComponentMass {
  int d;
  int rank;
  T det;
  std::optional<T> mass;
  // the multiplicity factors at p | |G| are beyond the bounds: mass is
  // only that of the multiplicity lattice X, a lower bound
  bool local_factor_missing = false;
};

template <typename T> struct IsometryMassResult {
  int order;
  std::vector<IsometryComponent<T>> components;
  std::vector<IsometryComponentMass<T>> component_masses;
  T index;
  std::vector<GlueOrbitResult<typename underlying_z_ring<T>::ring_type>> glue;
  T glue_factor;
  // the image of a local centralizer on a glue module was beyond the
  // bounds and taken to be all the equivariant isometries (Question q:surj
  // of the paper); only with IsometryMassBounds::allow_image_assumption
  bool image_assumed;
  // the multiplicity factors of a component are beyond the bounds
  // (IsometryComponentMass::local_factor_missing)
  bool local_correction_open = false;
  // a component of a type whose mass is not implemented
  bool unsupported_component = false;
  /*
    The mass of the genus of the pair. It is set only when it is proved:
    every component has its exact mass, with its local factors. When local
    factors are missing the product is mass_without_local_factors instead,
    a lower bound of the mass (each factor is an index, at least 1), which
    must never be used as the target of an enumeration: the enumeration
    could stop before it is complete.
  */
  std::optional<T> mass;
  std::optional<T> mass_without_local_factors;
};

struct IsometryMassBounds {
  size_t max_module = size_t(1) << 12;
  size_t max_group = 100000;
  size_t max_orbit = 100000;
  // the graph method: the number of vertices (module elements and glue)
  size_t max_vertices = 6000;
  // the L-values of the quadratic components with E not abelian
  QuadLValueMethod lvalue = QuadLValueMethod::Mpfr;
  // when the image of the local centralizer on a glue module is beyond
  // the bounds: stop (false), or take it to be all the equivariant
  // isometries (Question q:surj of the paper), the result then saying
  // image_assumed. A wrong assumption gives a wrong mass, which can stop
  // an enumeration before it is complete, so it is not done by default.
  bool allow_image_assumption = false;
};

/*
  The mass of a component according to the kind of its centralizer
  (ComponentKind, MassComponents.h), times the local factor at the primes
  dividing |G| when it was computed (multiplicity_factor). When it was not,
  the mass is returned without it and local_factor_missing is set, except
  for a lattice of rank one or a cyclic component, where there is none.
 */
template <typename T> struct ComponentMassValue {
  std::optional<T> mass;
  bool local_factor_missing = false;
};

template <typename T, typename Tint>
ComponentMassValue<T> ComponentMass(IsometryComponent<T> const &comp,
                                    IsometryMassBounds const &bounds) {
  std::optional<T> base;
  bool needs_factor = true;
  switch (comp.kind) {
  case ComponentKind::Orthogonal:
    return {GenusMass<T, Tint>(comp.gram), false};
  case ComponentKind::CyclotomicHermitian:
    return {HermitianMassCyclotomicLattice(comp.gram, comp.sigma, comp.d),
            false};
  case ComponentKind::Multiplicity:
    if (!comp.multiplicity_gram) {
      return {};
    }
    base = GenusMass<T, Tint>(*comp.multiplicity_gram);
    // multiplicity one: the centralizer on the component is {+-1}, which
    // preserves every lattice, so no local factor is missing
    needs_factor = comp.multiplicity_gram->rows() > 1;
    break;
  case ComponentKind::AbelianHermitian:
    base = HermitianMassAbelianLattice(comp.center_gram, comp.center_omega,
                                       *comp.center_field);
    break;
  case ComponentKind::AbelianQuadratic:
    base = QuadraticMassAbelianLattice(comp.center_gram, comp.center_omega,
                                       *comp.center_field, bounds.lvalue);
    break;
  case ComponentKind::Precomputed:
    base = comp.precomputed_mass;
    break;
  case ComponentKind::Unsupported:
    return {};
  }
  if (!base) {
    return {};
  }
  if (comp.multiplicity_factor) {
    return {*base * *comp.multiplicity_factor, false};
  }
  return {*base, needs_factor};
}

/*
  The mass from the components, for an isometry or an abelian group: the
  product of the masses of the components and of the glue orbits.
 */
template <typename T, typename Tint>
IsometryMassResult<T>
IsometryGenusMassComponents(MyMatrix<T> const &A, int order,
                            std::vector<IsometryComponent<T>> const &components,
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
    ComponentMassValue<T> v = ComponentMass<T, Tint>(comp, bounds);
    cm.mass = v.mass;
    if (v.local_factor_missing) {
      cm.local_factor_missing = true;
      res.local_correction_open = true;
    }
    if (!v.mass && comp.d == 0) {
      res.unsupported_component = true;
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
      {
        /*
          A cyclic component (d >= 3, one element acting) is a lattice over
          Z[sigma] = Z[zeta_d], the maximal order of Q(zeta_d), with an even
          trace form: the image of its local unitary group is the whole
          centralizer of sigma in O(D(L_d)_q) (Brandhorst-Hofmann, Finite
          subgroups of automorphisms of K3 surfaces, Theorem 6.12 with
          Proposition 6.9). An empty entry means all the equivariant
          isometries, for the listing and for the graph.
         */
        IsometryComponent<T> const &comp_c = res.components[comp_pos[s]];
        if (comp_c.d >= 3 && comp_c.extra_sigma.empty()) {
#ifdef SANITY_CHECK_ISOMETRY_MASS
          std::optional<std::vector<ModuleAut>> all_eq =
              EquivariantModuleIsometries(mods[s], bounds.max_module,
                                          bounds.max_group);
          std::optional<std::vector<ModuleAut>> img_c =
              LocalImageOnModule<T, Tint>(comp_c, mods[s], os);
          if (all_eq && img_c && all_eq->size() != img_c->size()) {
            std::cerr << "ISOM_MASS: SANITY: the image on a cyclic component "
                         "has " << img_c->size() << " elements, not all the "
                      << all_eq->size() << " equivariant isometries\n";
            throw TerminalException{1};
          }
#endif
          image_groups.push_back({});
          list_image.push_back("full");
          continue;
        }
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
        if (!bounds.allow_image_assumption) {
          std::cerr << "ISOM_MASS: the image of the local centralizer on the "
                       "glue module at q=" << q << " is beyond the bounds and "
                       "the image assumption (Question q:surj) is not "
                       "allowed, please insert the code\n";
          throw TerminalException{1};
        }
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
#ifdef SANITY_CHECK_ISOMETRY_MASS
    if (!opt) {
      std::cerr << "ISOM_MASS: the glue at q=" << q
                << " exceeds the bounds of both the listing and the graph\n";
      throw TerminalException{1};
    }
#endif
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
    if (res.local_correction_open) {
      res.mass_without_local_factors = mass0 * res.glue_factor;
    } else {
      res.mass = mass0 * res.glue_factor;
    }
  }
  return res;
}

/*
  The mass of the genus of (L, G), G the finite group generated by the
  isometries list_g of L (g^T A g = A, on the column vectors):

      mass(L, G) = sum_{(L', G)} 1 / |C_{O(L')}(G)|.

  One generator: the components are the L_d = L cap ker Phi_d(g)
  (IsometryComponents). Commuting generators: the components L_chi of
  AbelianGroupComponents, the glue isometries commuting with each g_chi.
  Otherwise: the components of FiniteGroupComponents, the glue commuting
  with all the generators. The order reported is |G|. The result says
  when local factors are left out (local_correction_open) or a component
  has no mass.
 */
template <typename T, typename Tint>
IsometryMassResult<T> GroupGenusMass(MyMatrix<T> const &A,
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
  if (ListP.size() == 1) {
    int order = IsometryOrder(ListP[0]);
    return IsometryGenusMassComponents<T, Tint>(
        A, order, IsometryComponents<T, Tint>(A, ListP[0], order), bounds, os);
  }
  bool commuting = true;
  for (size_t i = 0; i < ListP.size(); i++) {
    for (size_t j = i + 1; j < ListP.size(); j++) {
      if (ListP[i] * ListP[j] != ListP[j] * ListP[i]) {
        commuting = false;
      }
    }
  }
  std::vector<MyMatrix<T>> elts = GroupElementsClosure(ListP, bounds.max_group);
  std::vector<IsometryComponent<T>> components =
      commuting ? AbelianGroupComponents<T, Tint>(A, elts)
                : FiniteGroupComponents<T, Tint>(A, ListP, elts, os);
  return IsometryGenusMassComponents<T, Tint>(
      A, static_cast<int>(elts.size()), components, bounds, os);
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
      os << (cm.local_factor_missing ? ", mass_without_local_factors:="
                                     : ", mass:=")
         << *cm.mass;
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
  if (res.mass_without_local_factors) {
    os << ",\n  mass_without_local_factors:="
       << *res.mass_without_local_factors;
  }
  os << ");\n";
}

// clang-format off
#endif  // SRC_GENERA_MASSISOMETRY_H_
// clang-format on
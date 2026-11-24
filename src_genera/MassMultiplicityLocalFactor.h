// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_MASSMULTIPLICITYLOCALFACTOR_H_
#define SRC_GENERA_MASSMULTIPLICITYLOCALFACTOR_H_

// clang-format off
#include "MassGlue.h"
// clang-format on

/*
  The local factor of the mass of a component at the primes dividing |G|,
  prod_p [O(X_p) : Stab(L_p)] for the multiplicity lattice X
  (MultiplicityLocalFactor), from the images of the local orthogonal
  groups (MassLocalImage.h) acting on the glue modules (MassGlue.h).
 */
/*
  The product over the primes of [O(X_p) : Stab(L_p)].
  A: the Gram matrix of L; Xb: a basis of X (rows, coordinates of L); E: a
  basis of L_rho = L cap V_rho; elts: the elements of G (on the rows);
  OmegaX: a Z-basis of the order acting on X (O_K, or an order of a
  quaternion algebra D when quaternionic: the pi-adic search, which needs
  a commutative O_K, is then skipped, and the forms B(w x, y), w in the
  trace dual, are imposed in the search). Nothing when the bounds are
  exceeded.
 */
template <typename T, typename Tint>
std::optional<T> MultiplicityLocalFactor(MyMatrix<T> const &A,
                                         MyMatrix<T> const &Xb,
                                         MyMatrix<T> const &E,
                                         std::vector<MyMatrix<T>> const &elts,
                                         std::vector<MyMatrix<T>> const &ListP,
                                         std::vector<MyMatrix<T>> const &OmegaX,
                                         std::vector<MyMatrix<T>> const &ListCentral,
                                         bool quaternionic,
                                         MultiplicityFactorBounds const &bounds,
                                         [[maybe_unused]] std::ostream &os) {
  int m = Xb.rows();
  int k = E.rows();
  int n = A.cols();
  /*
    The failures below are not bounds: they mean that the multiplicity
    lattice does not have the expected structure, a case the code does not
    handle (or a bug of the decomposition). They stop the program rather
    than leave the factor out (IsometryComponentMass::local_factor_missing
    is for the bounds only).
  */
  auto structural_failure = [](std::string const &msg) -> void {
    std::cerr << "MULT_FACTOR: " << msg << ", please insert the code\n";
    throw TerminalException{1};
  };
  if (k % m != 0) {
    structural_failure("the rank of X does not divide that of the component");
  }
  int w = k / m;
  // the gamma_a with the gamma_a X in direct sum (V_rho = W (x)_K U and
  // X in l (x) U, l a line of W over the center K)
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
    structural_failure("only " + std::to_string(ListGamma.size()) + " of the " +
                       std::to_string(w) +
                       " translates of X found in direct sum by the greedy "
                       "choice");
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
    structural_failure("the translates of X do not span the component");
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
    structural_failure("Z[G] X has rank " + std::to_string(Nb.rows()) +
                       " instead of " + std::to_string(k));
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
        structural_failure("N^# / N is infinite");
      }
      int v = RationalValuation(d, p);
      if (v > 0) {
        ListIdx.push_back(i);
        ListOrd.push_back(ResiduePrimePower(p, v));
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
        int64_t o = ListOrd[s];
        res[s] = ResMulI64(ResModOfI64(GetNumerator_z(val), o),
                        ResInvI64(ResModOfI64(GetDenominator_z(val), o), o), o);
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
    // with a center K other than Q the search over the primes of K over
    // p, pi-adically, first: it is exact and the fastest
    std::optional<std::vector<ModRowsI64>> opt;
    bool rows_convention = !OmegaX.empty() && !quaternionic;
    if (!OmegaX.empty() && !quaternionic) {
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
            e = std::max(e, RationalValuation(den, p));
          }
        }
      }
    }
    ModRowsI64 A_mod(m, std::vector<int64_t>(m));
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++) {
        // exact (the generators read it modulo p^W, the listings modulo
        // p^{K+1}); the conversion throws beyond 64 bits
        A_mod[i][j] = UniversalScalarConversion<int64_t, T>(GramX(i, j));
      }
    }
    // which construction gave the image (for the sanity checks)
    bool gens_over_Q = false, gens_over_K = false;
    (void)gens_over_Q;
    (void)gens_over_K;
    // End_G = Q: generators of the image, complete by construction, with
    // p^e A^{-1} modulo p^W for the Newton steps
    if (!opt && OmegaX.empty()) {
      int Wn = MultiplicityNewtonPrecision(p);
      Tint pW = UniversalScalarConversion<Tint, int64_t>(1);
      pW *= MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(p), Wn);
      MyMatrix<T> Binv = Inverse(GramX);
      T pe = MyPow<T>(UniversalScalarConversion<T, int64_t>(p), e);
      ModRowsI64 Be(m, std::vector<int64_t>(m));
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
      std::optional<std::vector<ModRowsI64>> Lie = ResidueMatrices(
          IntegralLieLatticeColumns<T, Tint>(GramX, {}), ResidueModulusOf(pW));
      opt = LocalOrthogonalImageGenerators(A_mod, Be, {}, {}, {}, Lie, p, c, e,
                                           bounds, &order);
      gens_over_Q = opt.has_value();
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
      os << "MULT_FACTOR: p=" << p << " generators of the image "
         << (opt ? static_cast<long>(opt->size()) : -1L) << ", |image|="
         << LocalImageOrderValue(order, p) << "\n";
#endif
    }
    // a center K beyond the pi-adic search: the same generators within the
    // integral commutant of O_K (K-linear maps), on the columns
    if (!opt && !OmegaX.empty()) {
      int Wn = MultiplicityNewtonPrecision(p);
      Tint pW = UniversalScalarConversion<Tint, int64_t>(1);
      pW *= MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(p), Wn);
      MyMatrix<T> Binv = Inverse(GramX);
      T pe = MyPow<T>(UniversalScalarConversion<T, int64_t>(p), e);
      ModRowsI64 Be(m, std::vector<int64_t>(m));
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
      auto to_cols = [&](MyMatrix<T> const &Y) -> ModRowsI64 {
        ModRowsI64 Ym(m, std::vector<int64_t>(m));
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++) {
            // residues modulo p^W, which is all the generators read
            Tint v = ResInt(UniversalScalarConversion<Tint, T>(Y(j, i)), pW);
            if (v < 0) {
              v += pW;
            }
            Ym[i][j] = UniversalScalarConversion<int64_t, Tint>(v);
          }
        }
        return Ym;
      };
      std::vector<ModRowsI64> Om, Lam;
      for (auto const &W2 : OmegaX) {
        Om.push_back(to_cols(W2));
      }
      IsometryComponent<T> compX = BuildIsometryComponent<T, Tint>(
          A, Xb, IdentityMat<T>(A.rows()), 1);
      compX.extra_sigma = OmegaX;
      std::vector<MyMatrix<Tint>> CommColsX;
      for (auto const &Y : IntegralCommutant<T, Tint>(compX)) {
        Lam.push_back(to_cols(UniversalMatrixConversion<T, Tint>(Y)));
        CommColsX.push_back(MyMatrix<Tint>(Y.transpose()));
      }
      std::optional<std::vector<ModRowsI64>> LieX = ResidueMatrices(
          IntegralLieLatticeColumns<T, Tint>(GramX, CommColsX),
          ResidueModulusOf(pW));
      // for a quaternion order, the forms B(w x, y) = x W GramX y^T, w in
      // the trace dual (the integral points of their span), as residues
      // modulo p^W
      std::vector<ModRowsI64> Forms;
      if (quaternionic) {
        std::vector<MyMatrix<T>> Ls;
        for (auto const &W2 : OmegaX) {
          Ls.push_back(MyMatrix<T>(W2 * GramX));
        }
        for (auto const &N : QuatIntegralPoints<T, Tint>(Ls, m)) {
          ModRowsI64 Nm(m, std::vector<int64_t>(m));
          for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
              Tint v = ResInt(UniversalScalarConversion<Tint, T>(N(i, j)), pW);
              if (v < 0) {
                v += pW;
              }
              Nm[i][j] = UniversalScalarConversion<int64_t, Tint>(v);
            }
          }
          Forms.push_back(Nm);
        }
      }
      LocalImageOrder order;
      opt = LocalOrthogonalImageGenerators(A_mod, Be, Om, Lam, Forms, LieX, p,
                                           c, e, bounds, &order);
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
      // the basis as residues modulo the modulus of the listing
      std::optional<int64_t> listing_mod = LocalImageModulus(p, c, e);
      if (listing_mod) {
        std::vector<ModRowsI64> basis;
        for (auto const &Y : IntegralCommutant<T, Tint>(compX)) {
          ModRowsI64 Ym(m, std::vector<int64_t>(m));
          for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
              Ym[i][j] = CenteredResidueOfInteger(Y(i, j), *listing_mod);
            }
          }
          basis.push_back(Ym);
        }
        opt = LocalStabilizerImage(A_mod, basis, p, c, e, bounds);
      }
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
      pW *= MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(p), Wn);
      MyMatrix<T> Binv = Inverse(GramX);
      T pe = MyPow<T>(UniversalScalarConversion<T, int64_t>(p), e);
      ModRowsI64 Be(m, std::vector<int64_t>(m));
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
      auto to_cols = [&](MyMatrix<T> const &Y) -> ModRowsI64 {
        ModRowsI64 Ym(m, std::vector<int64_t>(m));
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++) {
            // residues modulo p^W, which is all the generators read
            Tint v = ResInt(UniversalScalarConversion<Tint, T>(Y(j, i)), pW);
            if (v < 0) {
              v += pW;
            }
            Ym[i][j] = UniversalScalarConversion<int64_t, Tint>(v);
          }
        }
        return Ym;
      };
      std::vector<ModRowsI64> Om, Lam;
      for (auto const &W2 : OmegaX) {
        Om.push_back(to_cols(W2));
      }
      IsometryComponent<T> compX3 = BuildIsometryComponent<T, Tint>(
          A, Xb, IdentityMat<T>(A.rows()), 1);
      compX3.extra_sigma = OmegaX;
      std::vector<MyMatrix<Tint>> CommCols3;
      for (auto const &Y : IntegralCommutant<T, Tint>(compX3)) {
        Lam.push_back(to_cols(UniversalMatrixConversion<T, Tint>(Y)));
        CommCols3.push_back(MyMatrix<Tint>(Y.transpose()));
      }
      std::optional<std::vector<ModRowsI64>> Lie3 = ResidueMatrices(
          IntegralLieLatticeColumns<T, Tint>(GramX, CommCols3),
          ResidueModulusOf(pW));
      LocalImageOrder order_k;
      std::optional<std::vector<ModRowsI64>> kg = LocalOrthogonalImageGenerators(
          A_mod, Be, Om, Lam, {}, Lie3, p, c, e, bounds, &order_k);
      if (kg) {
        int64_t pcl = ResiduePrimePower(p, c);
        auto closure_of = [&](std::vector<ModRowsI64> const &gens, bool transpose,
                              size_t cap) -> std::set<ModRowsI64> {
          ModRowsI64 Idm(m, std::vector<int64_t>(m, 0));
          for (int i = 0; i < m; i++) {
            Idm[i][i] = 1;
          }
          std::vector<ModRowsI64> gs;
          for (auto const &g : gens) {
            ModRowsI64 h(m, std::vector<int64_t>(m));
            for (int i = 0; i < m; i++) {
              for (int j = 0; j < m; j++) {
                h[i][j] = ResModI64(transpose ? g[j][i] : g[i][j], pcl);
              }
            }
            gs.push_back(h);
          }
          std::set<ModRowsI64> cl{Idm};
          std::vector<ModRowsI64> queue{Idm};
          while (!queue.empty() && cl.size() <= cap) {
            ModRowsI64 x = queue.back();
            queue.pop_back();
            for (auto const &g : gs) {
              ModRowsI64 y(m, std::vector<int64_t>(m, 0));
              for (int a = 0; a < m; a++) {
                for (int k2 = 0; k2 < m; k2++) {
                  for (int b = 0; b < m; b++) {
                    y[a][b] = ResAddI64(y[a][b], ResMulI64(x[a][k2], g[k2][b], pcl), pcl);
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
        std::set<ModRowsI64> listed_set;
        for (auto const &g : *opt) {
          ModRowsI64 h(m, std::vector<int64_t>(m));
          for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
              h[i][j] = ResModI64(g[j][i], pcl);
            }
          }
          listed_set.insert(h);
        }
        bool all_in = true;
        for (auto const &g : *kg) {
          ModRowsI64 h = g;
          for (auto &row : h) {
            for (auto &x : row) {
              x = ResModI64(x, pcl);
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
        std::set<ModRowsI64> c2 = closure_of(*opt, true, 2000000);
        std::set<ModRowsI64> c1 = closure_of(*kg, false, c2.size());
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
      std::optional<std::vector<ModRowsI64>> listed;
      if (gens_over_Q) {
        listed = LocalOrthogonalImage(A_mod, p, c, e, bounds);
      } else {
        IsometryComponent<T> compX2 = BuildIsometryComponent<T, Tint>(
            A, Xb, IdentityMat<T>(A.rows()), 1);
        compX2.extra_sigma = OmegaX;
        std::vector<ModRowsI64> basis2;
        std::optional<int64_t> listing_mod = LocalImageModulus(p, c, e);
        for (auto const &Y : listing_mod ? IntegralCommutant<T, Tint>(compX2)
                                         : std::vector<MyMatrix<Tint>>{}) {
          ModRowsI64 Ym(m, std::vector<int64_t>(m));
          for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
              Ym[i][j] = CenteredResidueOfInteger(Y(i, j), *listing_mod);
            }
          }
          basis2.push_back(Ym);
        }
        if (listing_mod) {
          listed = LocalStabilizerImage(A_mod, basis2, p, c, e, bounds);
        }
        if (listed) {
          // rows to columns
          for (auto &g : *listed) {
            ModRowsI64 h(m, std::vector<int64_t>(m));
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
        int64_t pcl = ResiduePrimePower(p, c);
        std::set<ModRowsI64> listed_set;
        for (auto const &g : *listed) {
          ModRowsI64 h = g;
          for (auto &row : h) {
            for (auto &x : row) {
              x = ResModI64(x, pcl);
            }
          }
          listed_set.insert(h);
        }
        ModRowsI64 Idm(m, std::vector<int64_t>(m, 0));
        for (int i = 0; i < m; i++) {
          Idm[i][i] = 1;
        }
        std::set<ModRowsI64> closure{Idm};
        std::vector<ModRowsI64> queue{Idm};
        while (!queue.empty() && closure.size() <= listed_set.size()) {
          ModRowsI64 x = queue.back();
          queue.pop_back();
          for (auto const &g : *opt) {
            ModRowsI64 y(m, std::vector<int64_t>(m, 0));
            for (int a = 0; a < m; a++) {
              for (int k2 = 0; k2 < m; k2++) {
                for (int b = 0; b < m; b++) {
                  y[a][b] = ResAddI64(y[a][b], ResMulI64(x[a][k2], ResModI64(g[k2][b], pcl), pcl), pcl);
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
    auto act_gens = [&](ModRowsI64 const &g) -> std::vector<std::vector<int64_t>> {
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
        if (x[s] != 0) {
          for (size_t u = 0; u < r; u++) {
            y[u] = ResAddI64(y[u], ResMulI64(x[s], act[s][u], ListOrd[u]), ListOrd[u]);
          }
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
    // the subgroups live in (Z / p^C)^r: p^C is the modulus, and x p^{C -
    // a_s} < p^C for x < p^{a_s}
    int64_t pC = ResiduePrimePower(p, Cexp);
    auto powp = [&](int k) -> int64_t { return ResiduePrimePower(p, k); };
    auto canon = [&](std::vector<std::vector<int64_t>> const &gens_orig)
        -> ModRowsI64 {
      ModRowsI64 sc;
      for (auto const &x : gens_orig) {
        std::vector<int64_t> y(r);
        for (size_t s2 = 0; s2 < r; s2++) {
          y[s2] = ResMulI64(x[s2], powp(Cexp - aexp[s2]), pC);
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
    ModRowsI64 L0c = canon(Lgens);
    std::set<ModRowsI64> orbit{L0c};
    std::vector<ModRowsI64> queue{L0c};
    while (!queue.empty()) {
      ModRowsI64 H = queue.back();
      queue.pop_back();
      for (auto const &act : ListAct) {
        std::vector<std::vector<int64_t>> img;
        for (auto const &row : H) {
          img.push_back(apply(act, unscale(row)));
        }
        ModRowsI64 Hc = canon(img);
        if (orbit.insert(Hc).second) {
          queue.push_back(Hc);
        }
      }
    }
    // |L/N| from the form: the product of the p^{C - v} of the pivots
    double logL = 0;
    (void)logL;
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

// clang-format off
#endif  // SRC_GENERA_MASSMULTIPLICITYLOCALFACTOR_H_
// clang-format on

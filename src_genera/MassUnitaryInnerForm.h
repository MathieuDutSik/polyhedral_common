// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_MASSUNITARYINNERFORM_H_
#define SRC_GENERA_MASSUNITARYINNERFORM_H_

// clang-format off
#include "MassAlbertTypeTwo.h"
#include <optional>
#include <set>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
# define SANITY_CHECK_UNITARY_INNER_FORM_MASS
#endif

#ifdef DEBUG
# define DEBUG_UNITARY_INNER_FORM_MASS
#endif

/*
  The mass of a component whose centralizer is a unitary group over a CM
  center K, U(A, tau) with A = End_G(V_rho) central simple of degree n over
  K and tau the adjoint involution, of the second kind: a form over K^+ of
  U_n, inner when A is not split (A = D a quaternion division algebra over
  K, a complex character of Schur index 2, n = 2). U is connected with a
  simply connected derived group, so its Tamagawa number is 2 (that of
  U_1); with Ono's measure, whose convergence factors are the L_v(1, chi)
  of chi = chi_{K/K^+},

      mass(X) = 2 prod_{i=1}^n L(i, chi^i) / (vol(U(R)) prod_p rho_p),

  L(i, chi^i) = zeta_{K^+}(i) for i even. At a prime v of K^+ where U is
  smooth, vol(U(O_v)) = |U(F_v)| q_v^{-n^2} = prod_{i=1}^n (1 - chi(v)^i
  q_v^{-i}), and the product of these over all v is 1 / prod L(i, chi^i);
  rho_p is the ratio of the computed volume at p (QuatLocalImageOrders, the
  measure of the trace form, as for type III) to that product over the v
  above p. At the real places U(R) = U(n)^{[K^+:Q]}, U(n) acting on X_R =
  (C^n)^mu, so that the metric -tr_X is s = 2 mu times -tr on C^n, and
  vol(U(n)) = s^{n^2/2} (2 pi)^{n(n+1)/2} / prod_{k<n} k!.

  The L-values from the functional equations, with A = sqrt(d_K / d_{K^+}):
  L(1, chi) = pi^c L(0, chi) / A, zeta_{K^+}(2) = 2^c pi^{2c} |zeta_{K^+}(-1)|
  / d_{K^+}^{3/2}, c = [K^+:Q]; the square of the mass is rational. n <= 2.

  X is the whole block (G spans the commutant of A): gens are matrices on
  X, columns convention, whose commutant is A.
 */
template <typename T, typename Tint>
std::optional<T> UnitaryInnerFormMass(MyMatrix<T> const &B,
                                      std::vector<MyMatrix<T>> const &gens,
                                      AbelianField const &FK, int n,
                                      QuaternionicMassBounds const &bounds,
                                      [[maybe_unused]] std::ostream &os) {
  if (n < 1 || n > 2) {
    return {};
  }
  int k = B.rows();
  AbelianField FKp = MakeAbelianRealField(FK.n, FK.H);
  int c = static_cast<int>(AbelianFieldDegree(FKp));
  int dim_u = n * n * c;
  std::vector<MyMatrix<T>> O0 = QuatIntegralPoints<T, Tint>(gens, k);
  std::vector<MyMatrix<T>> u = QuatIntegralCommutant<T, Tint>(B, O0, true);
#ifdef SANITY_CHECK_UNITARY_INNER_FORM_MASS
  if (static_cast<int>(u.size()) != dim_u) {
    std::cerr << "UNIT_MASS: the Lie lattice has rank " << u.size()
              << ", not " << dim_u << "\n";
    throw TerminalException{1};
  }
#endif
  T disc_u = QuatTraceDiscriminant(u);
  std::vector<MyMatrix<T>> Lam = QuatIntegralCommutant<T, Tint>(B, O0, false);
  std::vector<MyMatrix<T>> Forms;
  {
    std::vector<MyMatrix<T>> L;
    for (auto const &w : O0) {
      L.push_back(MyMatrix<T>(w.transpose() * B));
    }
    Forms = QuatIntegralPoints<T, Tint>(L, k);
  }
  // the discriminants of K and K^+, products of the conductors
  auto disc_of = [&](AbelianField const &F) -> T {
    DirichletGroup Gd = ComputeDirichletGroup(F.n);
    T d(1);
    for (auto const &t : FieldCharacters(Gd, F.H, false)) {
      d *= UniversalScalarConversion<T, int64_t>(PrimitiveCharacter(Gd, t).first);
    }
    return d;
  };
  T dK = disc_of(FK);
  T dKp = disc_of(FKp);
  // the bad primes
  auto prime_set = [&](T const &x) -> std::set<int64_t> {
    std::set<int64_t> Sp;
    for (Tint xi : {Tint(T_abs(GetNumerator_z(x))), Tint(GetDenominator_z(x))}) {
      if (xi > 1) {
        for (Tint const &q : FactorsInt(xi)) {
          Sp.insert(UniversalScalarConversion<int64_t, Tint>(q));
        }
      }
    }
    return Sp;
  };
  std::set<int64_t> primes{2};
  for (T const &x : {DeterminantMat(B), disc_u, dK}) {
    for (int64_t q : prime_set(x)) {
      primes.insert(q);
    }
  }
  {
    // the discriminant of the order of A on X, for its non-maximal primes
    int nl = static_cast<int>(Lam.size());
    MyMatrix<T> Gt(nl, nl);
    for (int i = 0; i < nl; i++) {
      for (int j = 0; j < nl; j++) {
        Gt(i, j) = QuatTrace(MyMatrix<T>(Lam[i] * Lam[j]));
      }
    }
    for (int64_t q : prime_set(DeterminantMat(Gt))) {
      primes.insert(q);
    }
  }
  // the generic volume prod_{v | p} prod_i (1 - chi(v)^i q_v^{-i})
  auto g_factor = [&](int64_t p) -> T {
    auto [eK, fK, gK] = AbelianPrimeDecomposition(FK.n, FK.H, p);
    auto [ep, fp, gp] = AbelianPrimeDecomposition(FKp.n, FKp.H, p);
    int chi = 0;
    if (gK == 2 * gp) {
      chi = 1;
    } else if (fK == 2 * fp) {
      chi = -1;
    } else if (eK == 2 * ep) {
      chi = 0;
    }
#ifdef SANITY_CHECK_UNITARY_INNER_FORM_MASS
    else {
      std::cerr << "UNIT_MASS: the decomposition of " << p << " in K and K^+ "
                << "does not fit\n";
      throw TerminalException{1};
    }
#endif
    T q = MyPow<T>(UniversalScalarConversion<T, int64_t>(p), fp);
    T res(1);
    for (int v = 0; v < gp; v++) {
      // chi^i, the trivial character for i even (also at a ramified v)
      T qi(1);
      for (int i = 1; i <= n; i++) {
        qi *= q;
        int ci = (i % 2 == 0) ? 1 : chi;
        res *= T(1) - UniversalScalarConversion<T, int64_t>(ci) / qi;
      }
    }
    return res;
  };
  std::vector<TamagawaLocalData<T>> loc;
  for (int64_t p : primes) {
    int cc = (p == 2) ? 2 : 1;
    std::optional<std::pair<Tint, Tint>> o = QuatLocalImageOrders<T, Tint>(
        B, O0, Lam, Forms, p, cc, bounds.image, os);
    if (!o) {
      return {};
    }
    TamagawaLocalData<T> ld;
    ld.p = p;
    ld.c = cc;
    ld.n_over_idx = UniversalScalarConversion<T, Tint>(o->first);
    ld.v_disc_u = RationalValuation(disc_u, p);
    ld.g_p = g_factor(p);
#ifdef DEBUG_UNITARY_INNER_FORM_MASS
    os << "UNIT_MASS: n=" << n << " p=" << p << " |U(X_p) mod p^c|=" << o->first
       << " v(disc u)=" << ld.v_disc_u << " g_p=" << ld.g_p << "\n";
#endif
    loc.push_back(ld);
  }
  // the L-values over pi^.., squared: L(1, chi)^2 / pi^{2c} =
  // L(0, chi)^2 d_{K^+} / d_K, zeta_{K^+}(2)^2 / pi^{4c} = 4^c zeta(-1)^2 /
  // d_{K^+}^3
  T L0 = StandardHermitianMassAbelian<T>(FK, 1);
  for (int i = 0; i < c - 1; i++) {
    L0 *= T(2);
  }
  T Lsq = L0 * L0 * dKp / dK;
  if (n >= 2) {
    T z = T_abs(DedekindZetaAbelian<T>(FKp, 2));
    T z2 = z * z / (dKp * dKp * dKp);
    for (int i = 0; i < c; i++) {
      z2 *= T(4);
    }
    Lsq *= z2;
  }
  // vol(U(n))^2 / pi^.. at each real place: s^{n^2} (2^{n(n+1)/2} / prod
  // k!)^2, s = 2 mu, mu = k / (2 n c)
  T s = UniversalScalarConversion<T, int64_t>(k) / UniversalScalarConversion<T, int64_t>(n * c);
  T V2(1);
  for (int pl = 0; pl < c; pl++) {
    for (int i = 0; i < n * n; i++) {
      V2 *= s;
    }
    T r(1);
    for (int i = 0; i < n * (n + 1) / 2; i++) {
      r *= T(2);
    }
    for (int i = 1; i < n; i++) {
      r /= UniversalScalarConversion<T, int64_t>(i);
    }
    V2 *= r * r;
  }
  // TamagawaMassSquared squares tau * Lprod / Vinf: pass 1 and multiply
  T m2 = TamagawaMassSquared<T>(T(2), T(1), T(1), dim_u, loc) * Lsq / V2;
  std::optional<T> mass = RationalSquareRoot<T, Tint>(m2);
#ifdef DEBUG_UNITARY_INNER_FORM_MASS
  os << "UNIT_MASS: n=" << n << " c=" << c << " d_K=" << dK << " d_K+=" << dKp
     << " mass^2=" << m2 << "\n";
#endif
#ifdef SANITY_CHECK_UNITARY_INNER_FORM_MASS
  if (!mass) {
    std::cerr << "UNIT_MASS: the squared mass " << m2 << " is not a square\n";
    throw TerminalException{1};
  }
#endif
  return mass;
}

// clang-format off
#endif  // SRC_GENERA_MASSUNITARYINNERFORM_H_
// clang-format on

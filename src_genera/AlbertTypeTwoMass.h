// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ALBERTTYPETWOMASS_H_
#define SRC_GENERA_ALBERTTYPETWOMASS_H_

// clang-format off
#include "QuaternionicMass.h"
#include "CyclotomicHermitianMass.h"
#include "GeneraArithmetic.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include <optional>
#include <set>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ALBERT_TYPE_TWO_MASS
#endif

/*
  The mass of a component V_rho of a finite group G whose End_G(V_rho) is
  a quaternion division algebra D over Q with an adjoint involution of
  orthogonal type (Albert type II, r = 1: a rational character with
  Frobenius-Schur indicator +1 and Schur index 2).

  X is a lattice in V_rho = D-module of rank one, D acting by the span of
  G (or by the e g e on X = L cap e V_rho), and the centralizer is the
  group N of the x in End_D(X) = D' (a quaternion algebra) with x^* x = 1,
  x^* = mu^{-1} x-bar mu for the skew-adjoint mu of D', mu^2 = -delta,
  delta > 0. N is an inner form of O_2: its identity component is the
  norm-one torus T of E = Q(mu), imaginary quadratic, and the other
  component is E j, j mu = -mu j. An element u = z j of that component has
  u^* u = -nrd(u) = N(z) j^2, so it exists over Q_p iff j^2 is a norm of
  E_p, iff D'_p is split; D' being a division algebra, N(Q) = T(Q).

  The classes of the genus of X under N are N(Q) \ N(A_f) / K_N, K_N the
  stabilizer of X. Since N(Q) = T(Q) lies in T(A_f), they split by the
  image in N(A_f) / T(A_f) K_N, which is the product of a Z / 2 over the
  primes p not in Ram(D') where K_{N, p} = K_{T, p}; over each coset the
  classes are those of the T-genus of a lattice with the same stabilizer
  in T (an element of the other component normalizes T and inverts it,
  so it preserves K_T), with the same finite groups T(Q) cap K. So

      mass_N(X) = 2^{s(X)} mass_T(X),
      mass_T(X) = mass_T(S) prod_p [T(O_{E, p}) : K_{T, p}(X)],

  s(X) the number of those primes, and S any lattice whose stabilizer in
  T is the maximal compact T(O_E-hat): S = O_E X. mass_T(S) is the mass
  of a Hermitian lattice of rank one over E (its unitary group is T, its
  stabilizer the maximal compact), RankOneHermitianMassAbelian.

  The local quantities come from the image orders of type III
  (QuatLocalImageOrders): |K mod p^c| for the stabilizer of the lattice in
  N (commuting with D) and in T (commuting with D and mu); the congruence
  subgroups lie in T, so [K_N : K_T] is the ratio of the two orders, and
  [T(O_{E, p}) : K_{T, p}(X)] is the ratio of the volumes of the T-
  stabilizers of S and X, with the discriminants of the Lie lattices. The
  primes are those dividing 2, det X, det S, disc D', disc E and the
  index of the Lie lattices; elsewhere X_p is unimodular and O_{E,p}-
  stable with an element of the other component in its stabilizer, which
  is checked at a few other primes under SANITY_CHECK_ALBERT_TYPE_TWO_MASS.
 */

struct AlbertTypeTwoMassBounds {
  QuaternionicMassBounds quat;
};

// The squarefree part of a positive integer.
template <typename Tint> Tint SquarefreePart(Tint const &x) {
  Tint res(1);
  if (x > Tint(1)) {
    std::map<Tint, size_t> fac = FactorsIntMap(x);
    for (auto const &[q, mult] : fac) {
      if (mult % 2 == 1) {
        res *= q;
      }
    }
  }
  return res;
}

/*
  The imaginary quadratic field Q(sqrt(-d)), d > 0 squarefree, as the
  abelian CM field of conductor |disc| fixed by the kernel of the
  Kronecker character of disc.
 */
inline AbelianCMField ImaginaryQuadraticField(int64_t d) {
  int64_t disc = (d % 4 == 3) ? -d : -4 * d;
  int64_t n = -disc;
  std::vector<int64_t> H;
  for (int64_t a = 1; a < n; a++) {
    if (std::gcd(a, n) == 1 && kronecker_symbol<int64_t>(disc, a) == 1) {
      H.push_back(a);
    }
  }
  return MakeAbelianCMField(n, H);
}

/*
  Whether End_G of the block (Gram matrix B, End_G and its center in the
  columns convention) is a quaternion division algebra over Q with an
  adjoint involution of orthogonal type: dimension 4 over the center Q,
  a skew-adjoint part of dimension 1 (3 for the symplectic type III), and
  a maximal order of discriminant > 1 (M_2(Q), type I with multiplicity
  two, has the same dimensions).
 */
template <typename T, typename Tint>
bool IsAlbertTypeTwoRankOne(MyMatrix<T> const &B,
                            std::vector<MyMatrix<T>> const &End,
                            std::vector<MyMatrix<T>> const &Center,
                            QuaternionicMassBounds const &bounds,
                            std::ostream &os) {
  int k = B.rows();
  if (Center.size() != 1 || End.size() != 4 || k % 4 != 0) {
    return false;
  }
  MyMatrix<T> Binv = Inverse(B);
  MyMatrix<T> M(4, k * k);
  for (int i = 0; i < 4; i++) {
    MyMatrix<T> Sy = End[i] + Binv * End[i].transpose() * B;
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        M(i, u * k + v) = Sy(u, v);
      }
    }
  }
  if (4 - RankMat(M) != 1) {
    return false;
  }
  std::vector<MyMatrix<T>> O = QuatIntegralPoints<T, Tint>(End, k);
  std::optional<std::vector<MyMatrix<T>>> opt_O =
      QuatMaximalOrder<T, Tint>(O, k / 4, T(1), bounds, os);
  return opt_O && QuatOrderDisc(*opt_O, k / 4) > T(1);
}

/*
  The mass of the component, columns convention: B the Gram matrix of X
  (rank 4), Dspan matrices spanning D on X over Q. Nothing when D is split
  (the component is then of type I, M_2(Q), and handled as such), or
  beyond the bounds of the image searches.
 */
template <typename T, typename Tint>
std::optional<T> AlbertTypeTwoRankOneMass(MyMatrix<T> const &B,
                                          std::vector<MyMatrix<T>> const &Dspan,
                                          AlbertTypeTwoMassBounds const &bounds,
                                          std::ostream &os) {
  int k = B.rows();
  if (k != 4) {
    std::cerr << "ALB2_MASS: X has rank " << k << ", not 4\n";
    throw TerminalException{1};
  }
  MyMatrix<T> Binv = Inverse(B);
  auto adj = [&](MyMatrix<T> const &Y) -> MyMatrix<T> {
    return Binv * Y.transpose() * B;
  };
  // the multiplier ring of X in D, and the maximal order for Ram(D)
  std::vector<MyMatrix<T>> O0 = QuatIntegralPoints<T, Tint>(Dspan, k);
  if (O0.size() != 4) {
    std::cerr << "ALB2_MASS: D has dimension " << O0.size() << "\n";
    throw TerminalException{1};
  }
  std::optional<std::vector<MyMatrix<T>>> opt_O =
      QuatMaximalOrder<T, Tint>(O0, 1, T(1), bounds.quat, os);
  if (!opt_O) {
    return {};
  }
  T discD = QuatOrderDisc(*opt_O, 1);
  Tint discD_i = UniversalScalarConversion<Tint, T>(discD);
  std::set<int64_t> ram;
  if (discD_i > 1) {
    for (Tint const &q : FactorsInt(discD_i)) {
      ram.insert(UniversalScalarConversion<int64_t, Tint>(q));
    }
  }
  if (ram.empty()) {
    // D = M_2(Q)
    return {};
  }
  // the D-linear endomorphisms and the skew-adjoint mu among them
  std::vector<MyMatrix<T>> uX = QuatIntegralCommutant<T, Tint>(B, O0, true);
  if (uX.size() != 1) {
    std::cerr << "ALB2_MASS: the skew-adjoint D-linear maps have rank "
              << uX.size() << ", not 1\n";
    throw TerminalException{1};
  }
  MyMatrix<T> const &mu = uX[0];
  MyMatrix<T> mu2 = mu * mu;
  T delta0 = -mu2(0, 0);
  if (mu2 != MyMatrix<T>(-delta0 * IdentityMat<T>(k)) || delta0 <= 0) {
    std::cerr << "ALB2_MASS: mu^2 is not a negative scalar\n";
    throw TerminalException{1};
  }
  if (adj(mu) != MyMatrix<T>(-mu)) {
    std::cerr << "ALB2_MASS: mu is not skew-adjoint\n";
    throw TerminalException{1};
  }
  // E = Q(sqrt(-d)), d the squarefree part of delta0 (an integer: mu is
  // integral on X)
  if (!IsInteger(delta0)) {
    std::cerr << "ALB2_MASS: mu^2 is not integral\n";
    throw TerminalException{1};
  }
  Tint delta_i = UniversalScalarConversion<Tint, T>(delta0);
  Tint d_i = SquarefreePart(delta_i);
  int64_t d = UniversalScalarConversion<int64_t, Tint>(d_i);
  Tint s2 = delta_i / d_i;
  Tint s_i = UniversalScalarConversion<Tint, double>(
      std::round(std::sqrt(UniversalScalarConversion<double, Tint>(s2))));
  if (s_i * s_i != s2) {
    std::cerr << "ALB2_MASS: delta / d is not a square\n";
    throw TerminalException{1};
  }
  AbelianCMField FE = ImaginaryQuadraticField(d);
  T massT = RankOneHermitianMassAbelian<T>(FE);
  // sqrt(-d) = mu / s and the generator omega of O_E
  MyMatrix<T> sqrtmd = mu / UniversalScalarConversion<T, Tint>(s_i);
  MyMatrix<T> omega = (d % 4 == 3)
                          ? MyMatrix<T>((IdentityMat<T>(k) + sqrtmd) / T(2))
                          : sqrtmd;
  // S = X + omega X, P its basis (columns) in the coordinates of X
  MyMatrix<T> P;
  {
    MyMatrix<T> Gen(2 * k, k);
    MyMatrix<T> omT = omega.transpose();
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < k; j++) {
        Gen(i, j) = (i == j) ? T(1) : T(0);
        Gen(k + i, j) = omT(i, j);
      }
    }
    P = GetZbasis(Gen).transpose();
  }
  MyMatrix<T> Pinv = Inverse(P);
  MyMatrix<T> BS = P.transpose() * B * P;
  auto to_S = [&](MyMatrix<T> const &Y) -> MyMatrix<T> {
    MyMatrix<T> Z = Pinv * Y * P;
    if (!IsIntegralMatrix(Z)) {
      std::cerr << "ALB2_MASS: an element does not preserve S = O_E X\n";
      throw TerminalException{1};
    }
    return Z;
  };
  std::vector<MyMatrix<T>> O0S;
  for (auto const &w : O0) {
    O0S.push_back(to_S(w));
  }
  MyMatrix<T> muS = to_S(mu);
  // the data of a lattice for the image searches, for N and for T
  struct LatticeData {
    MyMatrix<T> B;
    std::vector<MyMatrix<T>> gensN, gensT;
    std::vector<MyMatrix<T>> LamN, LamT, FormsN, FormsT;
    T disc_u;
  };
  auto forms = [&](MyMatrix<T> const &Bz, std::vector<MyMatrix<T>> const &W)
      -> std::vector<MyMatrix<T>> {
    std::vector<MyMatrix<T>> L;
    for (auto const &w : W) {
      L.push_back(MyMatrix<T>(w.transpose() * Bz));
    }
    return QuatIntegralPoints<T, Tint>(L, Bz.rows());
  };
  auto make_data = [&](MyMatrix<T> const &Bz, std::vector<MyMatrix<T>> const &gens,
                       MyMatrix<T> const &muz) -> LatticeData {
    LatticeData ld;
    ld.B = Bz;
    ld.gensN = gens;
    ld.gensT = gens;
    ld.gensT.push_back(muz);
    ld.LamN = QuatIntegralCommutant<T, Tint>(Bz, ld.gensN, false);
    ld.LamT = QuatIntegralCommutant<T, Tint>(Bz, ld.gensT, false);
    ld.FormsN = forms(Bz, gens);
    std::vector<MyMatrix<T>> WT = gens;
    for (auto const &w : gens) {
      WT.push_back(MyMatrix<T>(w * muz));
    }
    ld.FormsT = forms(Bz, WT);
    std::vector<MyMatrix<T>> u = QuatIntegralCommutant<T, Tint>(Bz, gens, true);
    if (u.size() != 1) {
      std::cerr << "ALB2_MASS: a Lie lattice of rank " << u.size() << "\n";
      throw TerminalException{1};
    }
    ld.disc_u = QuatTraceDiscriminant(u);
    return ld;
  };
  LatticeData dX = make_data(B, O0, mu);
  LatticeData dS = make_data(BS, O0S, muS);
  if (dS.LamT.size() != 2 || dX.LamT.size() != 2) {
    std::cerr << "ALB2_MASS: the commutant of D and mu has rank "
              << dX.LamT.size() << ", " << dS.LamT.size() << ", not 2\n";
    throw TerminalException{1};
  }
  // the primes of the numerator and of the denominator of x
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
  for (T const &x : {DeterminantMat(B), DeterminantMat(BS), discD,
                     QuatOrderDisc(O0, 1), delta0, dX.disc_u, dS.disc_u}) {
    for (int64_t q : prime_set(x)) {
      primes.insert(q);
    }
  }
#ifdef DEBUG_ALBERT_TYPE_TWO_MASS
  os << "ALB2_MASS: E = Q(sqrt(-" << d << ")), mass_T(max) = " << massT
     << ", Ram(D) =";
  for (int64_t q : ram) {
    os << " " << q;
  }
  os << ", det X = " << DeterminantMat(B) << ", det S = " << DeterminantMat(BS)
     << "\n";
#endif
  // the local factor at p, or nothing beyond the bounds
  auto local_factor = [&](int64_t p) -> std::optional<T> {
    int cc = (p == 2) ? 2 : 1;
    auto order = [&](LatticeData const &ld, bool torus) -> std::optional<Tint> {
      std::optional<std::pair<Tint, Tint>> o = QuatLocalImageOrders<T, Tint>(
          ld.B, torus ? ld.gensT : ld.gensN, torus ? ld.LamT : ld.LamN,
          torus ? ld.FormsT : ld.FormsN, p, cc, bounds.quat.image, os);
      if (!o) {
        return {};
      }
      return o->first;
    };
    std::optional<Tint> aX = order(dX, false);
    std::optional<Tint> bX = order(dX, true);
    std::optional<Tint> bS = order(dS, true);
    if (!aX || !bX || !bS) {
      return {};
    }
    // [K_N(X_p) : K_T(X_p)]
    if (*aX != *bX && *aX != Tint(2) * *bX) {
      std::cerr << "ALB2_MASS: at p=" << p << " the stabilizers in N and T "
                << "have images " << *aX << " and " << *bX << "\n";
      throw TerminalException{1};
    }
    bool in_ram = ram.count(p) == 1;
    if (in_ram && *aX != *bX) {
      std::cerr << "ALB2_MASS: at p=" << p << " in Ram(D) the stabilizer "
                << "meets the other component\n";
      throw TerminalException{1};
    }
    // [T(O_{E,p}) : K_T(X_p)] = vol_T(S_p) / vol_T(X_p)
    int vX = RationalValuation(dX.disc_u, p);
    int vS = RationalValuation(dS.disc_u, p);
    if ((vX - vS) % 2 != 0) {
      std::cerr << "ALB2_MASS: the Lie discriminants differ by an odd power "
                   "of " << p << "\n";
      throw TerminalException{1};
    }
    T index_T = UniversalScalarConversion<T, Tint>(*bS) /
                UniversalScalarConversion<T, Tint>(*bX);
    int ex = (vX - vS) / 2;
    T p_T = UniversalScalarConversion<T, int64_t>(p);
    for (int i = 0; i < (ex >= 0 ? ex : -ex); i++) {
      index_T = (ex >= 0) ? T(index_T * p_T) : T(index_T / p_T);
    }
    if (!IsInteger(index_T) || index_T < 1) {
      std::cerr << "ALB2_MASS: [T(O_E) : K_T(X)] = " << index_T << " at p="
                << p << "\n";
      throw TerminalException{1};
    }
    T factor = index_T;
    if (!in_ram && *aX == *bX) {
      factor *= 2;
    }
#ifdef DEBUG_ALBERT_TYPE_TWO_MASS
    os << "ALB2_MASS: p=" << p << (in_ram ? " (ramified in D)" : "")
       << " |K_N(X)|=" << *aX << " |K_T(X)|=" << *bX << " |K_T(S)|=" << *bS
       << " v(disc u)=" << vX << "," << vS << " [T(O_E):K_T(X)]=" << index_T
       << " factor=" << factor << "\n";
#endif
    return factor;
  };
  T mass = massT;
  for (int64_t p : primes) {
    std::optional<T> f = local_factor(p);
    if (!f) {
      return {};
    }
    mass *= *f;
  }
#ifdef SANITY_CHECK_ALBERT_TYPE_TWO_MASS
  {
    int64_t p = 2;
    int n_checked = 0;
    while (n_checked < 3) {
      p++;
      if (primes.count(p) == 1 || !IsPrime(p)) {
        continue;
      }
      std::optional<T> f = local_factor(p);
      if (f && *f != T(1)) {
        std::cerr << "ALB2_MASS: the factor at the good prime " << p << " is "
                  << *f << "\n";
        throw TerminalException{1};
      }
      n_checked++;
    }
  }
#endif
  return mass;
}

// clang-format off
#endif  // SRC_GENERA_ALBERTTYPETWOMASS_H_
// clang-format on

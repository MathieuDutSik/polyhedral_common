// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ALBERTTYPETWOMASS_H_
#define SRC_GENERA_ALBERTTYPETWOMASS_H_

// clang-format off
#include "QuaternionicMass.h"
#include "CyclotomicHermitianMass.h"
#include "GeneraArithmetic.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "MAT_MatrixModInt64.h"
#include <optional>
#include <set>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ALBERT_TYPE_TWO_MASS
#endif

/*
  The mass of a component V_rho of a finite group G whose End_G(V_rho) is
  M_r(D), D a quaternion division algebra with an adjoint involution of
  orthogonal type (Albert type II: a real character with Frobenius-Schur
  indicator +1 and Schur index 2). Rank r = 1 over a totally real abelian
  center K (below over Q, the changes over K with AlbertTypeTwoRankOneMass)
  and rank r >= 2 over Q (AlbertTypeTwoOrthogonalMass, by the Tamagawa
  number).

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
  // L(0, chi_{E/K}) for a non-abelian E
  QuadLValueMethod lvalue = QuadLValueMethod::Mpfr;
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
inline AbelianField ImaginaryQuadraticField(int64_t d) {
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
  Whether End_G = D (rank one), given in the columns convention with its
  center of degree c and the discriminant dK of the center, is a quaternion
  division algebra: a maximal order of its integral points has |disc| >
  dK^4 (M_2(K), type I with multiplicity two, has the same dimensions).
 */
template <typename T, typename Tint>
bool EndIsQuaternionDivision(std::vector<MyMatrix<T>> const &End, int k, int c,
                             T const &dK, QuaternionicMassBounds const &bounds,
                             std::ostream &os) {
  std::vector<MyMatrix<T>> O = QuatIntegralPoints<T, Tint>(End, k);
  if (static_cast<int>(O.size()) != 4 * c || k % (4 * c) != 0) {
    return false;
  }
  std::optional<std::vector<MyMatrix<T>>> opt_O =
      QuatMaximalOrder<T, Tint>(O, k / (4 * c), dK, bounds, os);
  return opt_O && QuatOrderDisc(*opt_O, k / (4 * c)) > dK * dK * dK * dK;
}

/*
  The rank r of a candidate of type II: End_G of dimension 4 c r^2 over a
  self-adjoint (totally real) center of degree c, with a skew-adjoint
  part of dimension c r (2 r - 1). The split M_{2r}(K) has the same
  dimensions: it is told apart by EndIsQuaternionDivision (r = 1) or on X
  (QuaternionSpanIsDivision). 0 otherwise.
 */
template <typename T>
int AlbertTypeTwoRankCandidate(MyMatrix<T> const &B,
                               std::vector<MyMatrix<T>> const &End,
                               std::vector<MyMatrix<T>> const &Center) {
  int k = B.rows();
  int e = static_cast<int>(End.size());
  int c = static_cast<int>(Center.size());
  if (c == 0 || e % (4 * c) != 0) {
    return 0;
  }
  int r = static_cast<int>(
      std::llround(std::sqrt(static_cast<double>(e / (4 * c)))));
  if (r < 1 || 4 * c * r * r != e) {
    return 0;
  }
  MyMatrix<T> Binv = Inverse(B);
  for (auto const &Z : Center) {
    if (MyMatrix<T>(Binv * Z.transpose() * B) != Z) {
      return 0;
    }
  }
  MyMatrix<T> M(e, k * k);
  for (int i = 0; i < e; i++) {
    MyMatrix<T> Sy = End[i] + Binv * End[i].transpose() * B;
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        M(i, u * k + v) = Sy(u, v);
      }
    }
  }
  if (e - RankMat(M) != c * r * (2 * r - 1)) {
    return 0;
  }
  return r;
}

// Whether the span of Dspan (on F of dimension 4 r) is a quaternion
// division algebra over Q, by the discriminant of a maximal order.
template <typename T, typename Tint>
bool QuaternionSpanIsDivision(std::vector<MyMatrix<T>> const &Dspan, int k,
                              int r, QuaternionicMassBounds const &bounds,
                              std::ostream &os) {
  std::vector<MyMatrix<T>> O = QuatIntegralPoints<T, Tint>(Dspan, k);
  if (O.size() != 4) {
    return false;
  }
  std::optional<std::vector<MyMatrix<T>>> opt_O =
      QuatMaximalOrder<T, Tint>(O, r, T(1), bounds, os);
  return opt_O && QuatOrderDisc(*opt_O, r) > T(1);
}


/*
  The maximal order of a commutative Q-algebra E (a field) given by a
  Z-basis O of an order (k x k integral matrices), by Round 2 at each p
  with p^2 dividing the discriminant: I the p-radical, the kernel of the
  Frobenius x -> x^q (q = p^j >= rank) on O / pO, which is additive, and
  O' = {x in E : x I subset I}, contained in p^{-1} O; O is p-maximal when
  O' = O. The elements are coordinate vectors on the given basis O, with
  its integral structure constants.
 */
template <typename T, typename Tint>
std::vector<MyMatrix<T>> CommutativeMaximalOrder(std::vector<MyMatrix<T>> const &O0,
                                                 int k) {
  int n = static_cast<int>(O0.size());
  MyMatrix<T> F(n, k * k);
  for (int i = 0; i < n; i++) {
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        F(i, u * k + v) = O0[i](u, v);
      }
    }
  }
  auto coords0 = [&](MyMatrix<T> const &Y) -> MyVector<T> {
    MyVector<T> y(k * k);
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        y(u * k + v) = Y(u, v);
      }
    }
    std::optional<MyVector<T>> s = SolutionMat(F, y);
    if (!s) {
      std::cerr << "ALB2_MASS: a product is not in the algebra\n";
      throw TerminalException{1};
    }
    return *s;
  };
  // the structure constants on O0: e_i e_j = sum_l C0[i][j](l) e_l
  std::vector<std::vector<MyVector<T>>> C0(n, std::vector<MyVector<T>>(n));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      C0[i][j] = coords0(MyMatrix<T>(O0[i] * O0[j]));
    }
  }
  auto mult = [&](MyVector<T> const &a, MyVector<T> const &b) -> MyVector<T> {
    MyVector<T> c = ZeroVector<T>(n);
    for (int i = 0; i < n; i++) {
      if (a(i) != 0) {
        for (int j = 0; j < n; j++) {
          if (b(j) != 0) {
            c += a(i) * b(j) * C0[i][j];
          }
        }
      }
    }
    return c;
  };
  // the current order: rows of Bo, coordinates on O0
  MyMatrix<T> Bo = IdentityMat<T>(n);
  auto disc = [&]() -> T {
    MyMatrix<T> G(n, n);
    std::vector<MyMatrix<T>> mats;
    for (int i = 0; i < n; i++) {
      MyMatrix<T> Y = ZeroMatrix<T>(k, k);
      for (int l = 0; l < n; l++) {
        Y += Bo(i, l) * O0[l];
      }
      mats.push_back(Y);
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        G(i, j) = QuatTrace(MyMatrix<T>(mats[i] * mats[j]));
      }
    }
    return T_abs(DeterminantMat(G));
  };
  T D = disc();
  std::vector<int64_t> primes;
  {
    Tint num = GetNumerator_z(D);
    if (num > 1) {
      for (auto const &[q, mult_q] : FactorsIntMap(num)) {
        if (mult_q >= 2) {
          primes.push_back(UniversalScalarConversion<int64_t, Tint>(q));
        }
      }
    }
  }
  for (int64_t p : primes) {
    T_small_prime pu = static_cast<T_small_prime>(p);
    T pT = ScalarFromInt64<T>(p);
    while (true) {
      MyMatrix<T> Boinv = Inverse(Bo);
      auto row = [&](int i) -> MyVector<T> {
        MyVector<T> v(n);
        for (int l = 0; l < n; l++) {
          v(l) = Bo(i, l);
        }
        return v;
      };
      auto in_basis = [&](MyVector<T> const &a, MyMatrix<T> const &Binv) -> MyVector<T> {
        return Binv.transpose() * a;
      };
      // the Frobenius on O / p O, by repeated multiplication of the basis
      int64_t q = p;
      while (q < n) {
        q *= p;
      }
      ModMatrix Frob(n, n);
      for (int i = 0; i < n; i++) {
        MyVector<T> x = row(i);
        MyVector<T> y = x;
        for (int64_t t = 1; t < q; t++) {
          y = mult(y, x);
          // reduce modulo p O to keep the entries small
          MyVector<T> cy = in_basis(y, Boinv);
          for (int l = 0; l < n; l++) {
            cy(l) = UniversalScalarConversion<T, T_residue>(
                ResModOf(GetNumerator_z(cy(l)), pu));
          }
          y = Bo.transpose() * cy;
        }
        MyVector<T> cy = in_basis(y, Boinv);
        for (int l = 0; l < n; l++) {
          if (!IsInteger(cy(l))) {
            std::cerr << "ALB2_MASS: the order is not closed\n";
            throw TerminalException{1};
          }
          Frob(i, l) = ResModOf(GetNumerator_z(cy(l)), pu);
        }
      }
      // the radical: the a with sum_i a_i x_i^q = 0 mod p
      ModMatrix rad = LeftNullspaceModP(Frob, pu);
      std::vector<MyVector<T>> gensI;
      for (int t = 0; t < rad.rows(); t++) {
        MyVector<T> a = ZeroVector<T>(n);
        for (int i = 0; i < n; i++) {
          a += UniversalScalarConversion<T, T_residue>(rad(t, i)) * row(i);
        }
        gensI.push_back(a);
      }
      for (int i = 0; i < n; i++) {
        gensI.push_back(MyVector<T>(pT * row(i)));
      }
      MyMatrix<T> Gi(static_cast<int>(gensI.size()), n);
      for (size_t t = 0; t < gensI.size(); t++) {
        for (int l = 0; l < n; l++) {
          Gi(static_cast<int>(t), l) = gensI[t](l);
        }
      }
      MyMatrix<T> BI = GetZbasis(Gi);
      MyMatrix<T> BIinv = Inverse(BI);
      // the a modulo p with (sum a_i x_i) y in p I for the basis y of I
      ModMatrix Cond(n, n * n);
      for (int i = 0; i < n; i++) {
        for (int t = 0; t < n; t++) {
          MyVector<T> yt(n);
          for (int l = 0; l < n; l++) {
            yt(l) = BI(t, l);
          }
          MyVector<T> c = in_basis(mult(row(i), yt), BIinv);
          for (int l = 0; l < n; l++) {
            if (!IsInteger(c(l))) {
              std::cerr << "ALB2_MASS: O I is not in I\n";
              throw TerminalException{1};
            }
            Cond(i, t * n + l) = ResModOf(GetNumerator_z(c(l)), pu);
          }
        }
      }
      ModMatrix ker = LeftNullspaceModP(Cond, pu);
      if (ker.rows() == 0) {
        break;
      }
      MyMatrix<T> Gn(n + ker.rows(), n);
      for (int i = 0; i < n; i++) {
        for (int l = 0; l < n; l++) {
          Gn(i, l) = Bo(i, l);
        }
      }
      for (int t = 0; t < ker.rows(); t++) {
        MyVector<T> a = ZeroVector<T>(n);
        for (int i = 0; i < n; i++) {
          a += UniversalScalarConversion<T, T_residue>(ker(t, i)) * row(i);
        }
        for (int l = 0; l < n; l++) {
          Gn(n + t, l) = a(l) / pT;
        }
      }
      MyMatrix<T> Bnew = GetZbasis(Gn);
      if (T_abs(DeterminantMat(Bnew)) == T_abs(DeterminantMat(Bo))) {
        break;
      }
      Bo = Bnew;
    }
  }
  std::vector<MyMatrix<T>> res;
  for (int i = 0; i < n; i++) {
    MyMatrix<T> Y = ZeroMatrix<T>(k, k);
    for (int l = 0; l < n; l++) {
      Y += Bo(i, l) * O0[l];
    }
    res.push_back(Y);
  }
  return res;
}

/*
  The mass of a Hermitian lattice of rank one over the CM field E, from
  the quadratic lattice M over K of rank two (Gram matrix BM, O_K = OmegaM
  on the rows) whose discriminant field is E: 2^{1 - [K:Q]} |L(0,
  chi_{E/K})| 2^{-t}, t the number of finite primes of K ramified in E.
  E abelian: RankOneHermitianMassAbelian of E as a subfield of a
  cyclotomic field (K(sqrt D) for a rational D, or as found by
  AbelianQuadraticExtensionField); otherwise L(0, chi_{E/K}) numerically,
  rounded as in QuadraticMassNonAbelianE.
 */
template <typename T>
std::optional<T> RankOneHermitianMassOverK(MyMatrix<T> const &BM,
                                           std::vector<MyMatrix<T>> const &OmegaM,
                                           AbelianField const &FK,
                                           std::set<int64_t> const &ram_E,
                                           QuadLValueMethod method,
                                           [[maybe_unused]] std::ostream &os) {
  std::optional<std::pair<NumberField<T>, ElemE<T>>> Kd =
      QuadDiscriminantK(BM, OmegaM);
  if (!Kd) {
    return {};
  }
  NumberField<T> const &K = Kd->first;
  ElemE<T> const &d = Kd->second;
  if (SqrtInK(d, K)) {
    std::cerr << "ALB2_MASS: K(mu) is not a field\n";
    throw TerminalException{1};
  }
  std::optional<AbelianField> FE;
  {
    // the primes ramified in E: a rational D with K(sqrt d) = K(sqrt D)
    // is a product of them (and of those of 2 n)
    using Tint = typename underlying_z_ring<T>::ring_type;
    std::set<int64_t> ps = ram_E;
    for (Tint const &p :
         FactorsInt(UniversalScalarConversion<Tint, int64_t>(2 * FK.n))) {
      ps.insert(UniversalScalarConversion<int64_t, Tint>(p));
    }
    std::vector<int64_t> odd(ps.begin(), ps.end());
    odd.erase(std::remove(odd.begin(), odd.end(), 2), odd.end());
    for (uint64_t sub = 0; sub < (uint64_t(1) << odd.size()) && !FE; sub++) {
      int64_t Dodd = 1;
      for (size_t i = 0; i < odd.size(); i++) {
        if ((sub >> i) & 1) {
          Dodd *= (odd[i] % 4 == 1) ? odd[i] : -odd[i];
        }
      }
      for (int64_t two : {1, -4, 8, -8}) {
        int64_t Dq = Dodd * two;
        if (Dq >= 0) {
          continue;
        }
        if (!SqrtInK(EScal(d, UniversalScalarConversion<T, int64_t>(Dq)), K)) {
          continue;
        }
        int64_t N = std::lcm(FK.n, std::abs(Dq));
        AbelianField FK2 = LiftAbelianField(FK, N);
        std::vector<int64_t> HE;
        for (int64_t a : FK2.H) {
          if (kronecker_symbol<int64_t>(Dq, a) == 1) {
            HE.push_back(a);
          }
        }
        FE = AbelianField{N, HE};
        break;
      }
    }
  }
  if (!FE) {
    FE = AbelianQuadraticExtensionFieldPrimes(BM, OmegaM, FK, d, K, 2, ram_E);
  }
  if (FE) {
#ifdef DEBUG_ALBERT_TYPE_TWO_MASS
    os << "ALB2_MASS: E abelian, of conductor " << FE->n << "\n";
#endif
    return RankOneHermitianMassAbelian<T>(*FE);
  }
  // E not abelian: t from the Hecke character, L(0) numerically
  auto f_none = [](int64_t, std::vector<QuadHeckePrime> const &) -> void {};
  std::optional<QuadHeckeCharacter<T>> C =
      ComputeQuadHeckeCharacter(BM, OmegaM, FK, d, K, 2, 2, false, f_none);
  if (!C) {
    return {};
  }
  int t = 0;
  for (auto const &[p, lp] : C->primes) {
    for (auto const &P : lp) {
      if (P.chi == 0) {
        t++;
      }
    }
  }
  T local(2);
  for (int i = 0; i < t; i++) {
    local /= T(2);
  }
#ifdef DEBUG_ALBERT_TYPE_TWO_MASS
  os << "ALB2_MASS: E not abelian, " << t << " ramified primes of K\n";
#endif
  return QuadraticMassNonAbelianE(BM, OmegaM, FK, d, K, 2, local, method);
}

/*
  The mass of the component, columns convention: B the Gram matrix of X
  (rank 4 c), Dspan matrices spanning D on X over Q, FK the totally real
  center K when c = [K:Q] > 1 and dK its discriminant. Over K the primes P
  of K over p outside Ram(D) each contribute a 2 when the stabilizer at P
  does not meet the other component: 2^{g'_p - a_p} at p, g'_p their
  number and 2^{a_p} = [K_{N,p} : K_{T,p}]. Nothing when D is split (the
  component is then of type I and handled as such), or beyond the bounds.
 */
template <typename T, typename Tint>
std::optional<T> AlbertTypeTwoRankOneMass(MyMatrix<T> const &B,
                                          std::vector<MyMatrix<T>> const &Dspan,
                                          std::optional<AbelianField> const &FK,
                                          T const &dK,
                                          AlbertTypeTwoMassBounds const &bounds,
                                          std::ostream &os) {
  int k = B.rows();
  int c = FK ? static_cast<int>(AbelianFieldDegree(*FK)) : 1;
  if (k != 4 * c) {
    std::cerr << "ALB2_MASS: X has rank " << k << ", not 4 c = " << 4 * c
              << "\n";
    throw TerminalException{1};
  }
  MyMatrix<T> Binv = Inverse(B);
  auto adj = [&](MyMatrix<T> const &Y) -> MyMatrix<T> {
    return Binv * Y.transpose() * B;
  };
  // the multiplier ring of X in D, and the maximal order for Ram(D)
  std::vector<MyMatrix<T>> O0 = QuatIntegralPoints<T, Tint>(Dspan, k);
  if (static_cast<int>(O0.size()) != 4 * c) {
    std::cerr << "ALB2_MASS: D has dimension " << O0.size() << "\n";
    throw TerminalException{1};
  }
  std::optional<std::vector<MyMatrix<T>>> opt_O =
      QuatMaximalOrder<T, Tint>(O0, 1, dK, bounds.quat, os);
  if (!opt_O) {
    return {};
  }
  // N(d(D))^2, and the number of primes of K over p in Ram(D)
  T dK4 = dK * dK * dK * dK;
  T discD = QuatOrderDisc(*opt_O, 1) / dK4;
  Tint discD_i = UniversalScalarConversion<Tint, T>(discD);
  std::set<int64_t> ram;
  std::map<int64_t, int> g_ram;
  auto decomposition = [&](int64_t p) -> std::tuple<int, int, int> {
    if (!FK) {
      return {1, 1, 1};
    }
    return AbelianPrimeDecomposition(FK->n, FK->H, p);
  };
  if (discD_i > 1) {
    for (Tint const &q : FactorsInt(discD_i)) {
      int64_t p = UniversalScalarConversion<int64_t, Tint>(q);
      auto [ee, ff, gg] = decomposition(p);
      (void)ee;
      int v = RationalValuation(discD, p);
      if (v % (2 * ff) != 0 || v / (2 * ff) > gg) {
        std::cerr << "ALB2_MASS: v_p(N(d(D))^2)=" << v << " at p=" << p
                  << " does not fit f=" << ff << " g=" << gg << "\n";
        throw TerminalException{1};
      }
      ram.insert(p);
      g_ram[p] = v / (2 * ff);
    }
  }
  if (ram.empty()) {
    // D = M_2(K)
    return {};
  }
  // the D-linear endomorphisms and the skew-adjoint mu among them
  std::vector<MyMatrix<T>> uX = QuatIntegralCommutant<T, Tint>(B, O0, true);
  if (static_cast<int>(uX.size()) != c) {
    std::cerr << "ALB2_MASS: the skew-adjoint D-linear maps have rank "
              << uX.size() << ", not " << c << "\n";
    throw TerminalException{1};
  }
  MyMatrix<T> const &mu = uX[0];
  MyMatrix<T> mu2 = mu * mu;
  T delta0 = -mu2(0, 0);
  if (adj(mu) != MyMatrix<T>(-mu)) {
    std::cerr << "ALB2_MASS: mu is not skew-adjoint\n";
    throw TerminalException{1};
  }
  if (c == 1 && (mu2 != MyMatrix<T>(-delta0 * IdentityMat<T>(k)) ||
                 delta0 <= 0)) {
    std::cerr << "ALB2_MASS: mu^2 is not a negative scalar\n";
    throw TerminalException{1};
  }
  // O_E (generators on X, 1 first) and the mass of the maximal compact
  std::vector<MyMatrix<T>> OE;
  T massT;
  std::string E_name;
  if (c == 1) {
    // E = Q(sqrt(-d)), d the squarefree part of delta0 (an integer: mu
    // is integral on X)
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
    massT = RankOneHermitianMassAbelian<T>(ImaginaryQuadraticField(d));
    // sqrt(-d) = mu / s and the generator omega of O_E
    MyMatrix<T> sqrtmd = mu / UniversalScalarConversion<T, Tint>(s_i);
    MyMatrix<T> omega = (d % 4 == 3)
                            ? MyMatrix<T>((IdentityMat<T>(k) + sqrtmd) / T(2))
                            : sqrtmd;
    OE = {IdentityMat<T>(k), omega};
    E_name = "Q(sqrt(-" + std::to_string(d) + "))";
  } else {
    // O_E by Round 2 from the integral points of E = K(mu) on X
    std::vector<MyMatrix<T>> gE = O0;
    gE.push_back(mu);
    std::vector<MyMatrix<T>> LamE = QuatIntegralCommutant<T, Tint>(B, gE, false);
    if (static_cast<int>(LamE.size()) != 2 * c) {
      std::cerr << "ALB2_MASS: E has rank " << LamE.size() << ", not "
                << 2 * c << "\n";
      throw TerminalException{1};
    }
    OE = CommutativeMaximalOrder<T, Tint>(LamE, k);
    // M = O_E x, a Hermitian lattice of rank one over E, as the quadratic
    // lattice of rank two over K whose discriminant field is E. Its
    // discriminant carries N(h(x, x))^2, whose primes would all enter the
    // conductor bound of the search for E: x is taken among e_i and e_i +-
    // e_j so that det(BM) has the smallest largest prime factor.
    int nE = 2 * c;
    auto make_Pm = [&](MyVector<T> const &x) -> MyMatrix<T> {
      MyMatrix<T> Pm(k, nE);
      for (int i = 0; i < nE; i++) {
        MyVector<T> y = OE[i] * x;
        for (int a = 0; a < k; a++) {
          Pm(a, i) = y(a);
        }
      }
      return Pm;
    };
    auto largest_prime = [&](T const &x) -> Tint {
      Tint res(1);
      for (Tint xi : {Tint(T_abs(GetNumerator_z(x))), Tint(GetDenominator_z(x))}) {
        if (xi > 1) {
          for (Tint const &q : FactorsInt(xi)) {
            if (q > res) {
              res = q;
            }
          }
        }
      }
      return res;
    };
    std::optional<MyMatrix<T>> best_Pm;
    Tint best_lp(0);
    T best_det(0);
    for (int i = 0; i < k; i++) {
      for (int j = i; j < k; j++) {
        for (int sgn : {1, -1}) {
          if (j == i && sgn == -1) {
            continue;
          }
          MyVector<T> x = ZeroVector<T>(k);
          x(i) = T(1);
          if (j != i) {
            x(j) = ScalarFromInt64<T>(sgn);
          }
          MyMatrix<T> Pm = make_Pm(x);
          T det = T_abs(DeterminantMat(MyMatrix<T>(Pm.transpose() * B * Pm)));
          if (det == 0) {
            continue;
          }
          Tint lp = largest_prime(det);
          if (!best_Pm || lp < best_lp || (lp == best_lp && det < best_det)) {
            best_Pm = Pm;
            best_lp = lp;
            best_det = det;
          }
        }
      }
    }
    MyMatrix<T> const &Pm = *best_Pm;
    MyMatrix<T> BM = Pm.transpose() * B * Pm;
#ifdef DEBUG_ALBERT_TYPE_TWO_MASS
    os << "ALB2_MASS: M = O_E x with det " << best_det << ", largest prime "
       << best_lp << "\n";
#endif
    std::vector<MyMatrix<T>> CenterK;
    {
      std::vector<MyMatrix<T>> Zc;
      int d0 = static_cast<int>(O0.size());
      int kk = k * k;
      MyMatrix<T> Eq(d0, kk * d0);
      for (int i = 0; i < d0; i++) {
        for (int j = 0; j < d0; j++) {
          MyMatrix<T> Cm = O0[i] * O0[j] - O0[j] * O0[i];
          for (int u = 0; u < k; u++) {
            for (int v = 0; v < k; v++) {
              Eq(i, j * kk + u * k + v) = Cm(u, v);
            }
          }
        }
      }
      MyMatrix<T> ker = NullspaceMat(Eq);
      for (int t = 0; t < ker.rows(); t++) {
        MyMatrix<T> Z = ZeroMatrix<T>(k, k);
        for (int i = 0; i < d0; i++) {
          Z += ker(t, i) * O0[i];
        }
        Zc.push_back(Z);
      }
      CenterK = CommutativeMaximalOrder<T, Tint>(
          QuatIntegralPoints<T, Tint>(Zc, k), k);
    }
    // O_K on M, on the rows
    std::vector<MyMatrix<T>> OmegaM;
    for (auto const &z : CenterK) {
      MyMatrix<T> Zm(nE, nE);
      for (int i = 0; i < nE; i++) {
        MyVector<T> col(k);
        for (int a = 0; a < k; a++) {
          col(a) = T(0);
          for (int b2 = 0; b2 < k; b2++) {
            col(a) += z(a, b2) * Pm(b2, i);
          }
        }
        std::optional<MyVector<T>> sol =
            SolutionMat(MyMatrix<T>(Pm.transpose()), col);
        if (!sol) {
          std::cerr << "ALB2_MASS: O_K does not preserve O_E x\n";
          throw TerminalException{1};
        }
        for (int j = 0; j < nE; j++) {
          Zm(j, i) = (*sol)(j);
        }
      }
      OmegaM.push_back(MyMatrix<T>(Zm.transpose()));
    }
    // the primes ramified in E, those of disc(O_E)
    std::set<int64_t> ram_E;
    {
      MyMatrix<T> Gt(nE, nE);
      for (int i = 0; i < nE; i++) {
        for (int j = 0; j < nE; j++) {
          Gt(i, j) = QuatTrace(MyMatrix<T>(OE[i] * OE[j]));
        }
      }
      T dE = T_abs(DeterminantMat(Gt));
      Tint dEi = GetNumerator_z(dE);
      // tr_X = (k / 2c) Tr_{E/Q}: only the primes of k / 2c come in besides
      if (dEi > 1) {
        for (Tint const &q : FactorsInt(dEi)) {
          ram_E.insert(UniversalScalarConversion<int64_t, Tint>(q));
        }
      }
    }
    std::optional<T> mE = RankOneHermitianMassOverK<T>(BM, OmegaM, *FK, ram_E,
                                                       bounds.lvalue, os);
    if (!mE) {
      return {};
    }
    massT = *mE;
    E_name = "K(mu) of degree " + std::to_string(2 * c);
  }
  // S = O_E X, P its basis (columns) in the coordinates of X
  MyMatrix<T> P;
  {
    int nO = static_cast<int>(OE.size());
    MyMatrix<T> Gen(nO * k, k);
    for (int t = 0; t < nO; t++) {
      MyMatrix<T> wT = OE[t].transpose();
      for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
          Gen(t * k + i, j) = wT(i, j);
        }
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
    if (static_cast<int>(u.size()) != c) {
      std::cerr << "ALB2_MASS: a Lie lattice of rank " << u.size() << "\n";
      throw TerminalException{1};
    }
    ld.disc_u = QuatTraceDiscriminant(u);
    return ld;
  };
  LatticeData dX = make_data(B, O0, mu);
  LatticeData dS = make_data(BS, O0S, muS);
  if (static_cast<int>(dS.LamT.size()) != 2 * c ||
      static_cast<int>(dX.LamT.size()) != 2 * c) {
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
  for (T const &x : {DeterminantMat(B), DeterminantMat(BS), discD, dK,
                     QuatOrderDisc(O0, 1), DeterminantMat(mu), dX.disc_u,
                     dS.disc_u}) {
    for (int64_t q : prime_set(x)) {
      primes.insert(q);
    }
  }
#ifdef DEBUG_ALBERT_TYPE_TWO_MASS
  os << "ALB2_MASS: E = " << E_name << ", mass_T(max) = " << massT
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
    // [K_N(X_p) : K_T(X_p)] = 2^{a_p}, at most the number g'_p of the
    // primes over p outside Ram(D), where the other component has points
    bool in_ram = ram.count(p) == 1;
    int g_p = std::get<2>(decomposition(p));
    int g_good = g_p - (in_ram ? g_ram[p] : 0);
    int a_p = 0;
    {
      Tint q = *aX;
      if (q % *bX != 0) {
        std::cerr << "ALB2_MASS: at p=" << p << " the image in T does not "
                  << "divide that in N\n";
        throw TerminalException{1};
      }
      q /= *bX;
      while (q > 1 && q % 2 == 0) {
        q /= 2;
        a_p++;
      }
      if (q != 1 || a_p > g_good) {
        std::cerr << "ALB2_MASS: at p=" << p << " the stabilizers in N and T "
                  << "have images " << *aX << " and " << *bX << " for "
                  << g_good << " primes outside Ram(D)\n";
        throw TerminalException{1};
      }
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
    for (int i = a_p; i < g_good; i++) {
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

/*
  Rank r >= 2 over D, by the Tamagawa number. N is the group of the x in
  A = End_D(X) = M_r(D') with x^* x = 1, an inner form of O_{2r}; its
  identity component N^0, an inner form of SO_{2r}, is semisimple with
  Tamagawa number 2 (Ono: |Pic| = 2, Sha trivial), so with the Haar
  measures of the trace form tr_X(Y Z) on the Lie algebra u (the measure
  of QuatLocalImageOrders at the finite primes, whose product with the
  real one is the Tamagawa measure by the product formula),

      mass_{N^0}(X) = 2 / (vol(N^0(R)) prod_p vol(K_{N^0, p})).

  N^0(R) = SO(2r) acting on X_R = (R^{2r})^2, so the metric -tr_X is 4
  times -tr/2 and vol(N^0(R)) = 2^{dim u} prod_{j=1}^{2r-1} vol(S^j). At a
  prime where N^0 is smooth over Z_p, vol(K_p) = |N^0(F_p)| p^{-dim u} =
  prod_{i<r} (1 - p^{-2i}) (1 - chi(p) p^{-r}), chi the character of the
  discriminant field F = Q(sqrt(disc sigma)), disc sigma = (-1)^r Nrd(a)
  for a skew-adjoint a; the product of these factors over all p is
  1 / (prod_{i<r} zeta(2i) L(r, chi)), and the bad primes are corrected by
  the ratios of the computed volumes to them. The square of the mass is
  rational: sqrt(f) in L(r, chi) and the half powers |disc u|_p^{1/2}.

  From N^0 to N as in rank one: an improper isometry exists over Q_p only
  where A_p is split, and not over Q (A is not split), so the mass of N
  is 2^{|B|} mass_{N^0}, B the primes outside Ram(D) where the stabilizer
  is proper. Properness is read off the reduced norm, det of g on e X / p X
  for a primitive idempotent e of D modulo p, on the generators of the
  image of the stabilizer (p odd; at p = 2 outside Ram(D) nothing is
  returned).
 */

// vol(S^j) / pi^{ceil(j / 2)}
template <typename T> T SphereVolumeRational(int j) {
  T res(2);
  if (j % 2 == 1) {
    int m = (j + 1) / 2;
    for (int i = 2; i < m; i++) {
      res /= ScalarFromInt64<T>(i);
    }
  } else {
    for (int i = 0; i < j / 2; i++) {
      res *= T(2);
    }
    for (int i = j - 1; i >= 1; i -= 2) {
      res /= ScalarFromInt64<T>(i);
    }
  }
  return res;
}

// zeta(2 i) / pi^{2 i}
template <typename T> T ZetaEvenOverPi(int i) {
  std::vector<T> Bern = BernoulliNumbers<T>(2 * i);
  T res = T_abs(Bern[2 * i]);
  for (int t = 0; t < 2 * i; t++) {
    res *= T(2);
  }
  res /= T(2);
  for (int t = 2; t <= 2 * i; t++) {
    res /= ScalarFromInt64<T>(t);
  }
  return res;
}

// (L(r, chi_Df) / pi^r)^2, chi_Df of parity r (Df = 1: zeta(r), r even)
template <typename T> T LValueOverPiSquared(int r, int64_t Df) {
  int64_t f = (Df < 0) ? -Df : Df;
  T Brc = GeneralizedBernoulli<T>(r, Df);
  T fT = ScalarFromInt64<T>(f);
  T res = fT / T(4);
  for (int t = 0; t < r; t++) {
    res *= T(4) / (fT * fT);
  }
  T fact(1);
  for (int t = 2; t <= r; t++) {
    fact *= ScalarFromInt64<T>(t);
  }
  return res * Brc * Brc / (fact * fact);
}

// The fundamental discriminant of Q(sqrt(x)), x a nonzero integer.
template <typename Tint> int64_t FundamentalDiscriminantOf(Tint const &x) {
  Tint ax = (x < 0) ? Tint(-x) : x;
  Tint sq = SquarefreePart(ax);
  int64_t d = UniversalScalarConversion<int64_t, Tint>(sq);
  if (x < 0) {
    d = -d;
  }
  if (d == 1) {
    return 1;
  }
  int64_t m4 = ((d % 4) + 4) % 4;
  return (m4 == 1) ? d : 4 * d;
}

// The exact square root of a nonnegative rational, if it is a square.
template <typename T, typename Tint>
std::optional<T> RationalSquareRoot(T const &x) {
  if (x < 0) {
    return {};
  }
  Tint num = GetNumerator_z(x), den = GetDenominator_z(x);
  auto isqrt = [](Tint const &n) -> std::optional<Tint> {
    Tint s = UniversalScalarConversion<Tint, double>(
        std::floor(std::sqrt(UniversalScalarConversion<double, Tint>(n))));
    for (Tint t = s - 2; t <= s + 2; t += 1) {
      if (t >= 0 && t * t == n) {
        return t;
      }
    }
    return {};
  };
  std::optional<Tint> a = isqrt(num), b = isqrt(den);
  if (!a || !b) {
    return {};
  }
  return UniversalScalarConversion<T, Tint>(*a) /
         UniversalScalarConversion<T, Tint>(*b);
}

/*
  The reduced norms modulo the odd prime p, not ramified in D, of the
  image generators (columns, residues modulo p^W) of a stabilizer
  commuting with the order O0 of D: the determinants on e X / p X, e a
  primitive idempotent of O0 / p = M_2(F_p).
 */
template <typename T>
std::optional<std::vector<T_residue>>
ReducedNormsModP(std::vector<MyMatrix<T>> const &O0,
                 std::vector<ModRowsI64> const &image_gens, int k,
                 T_small_prime p) {
  // a primitive idempotent from an element with two distinct roots
  std::optional<ModMatrix> e;
  ModMatrix Id = ModMatrixIdentity(k);
  int n_o = static_cast<int>(O0.size());
  std::vector<ModMatrix> O0p;
  for (auto const &w : O0) {
    O0p.push_back(ModMatrixOfMatrix(w, p));
  }
  std::vector<int> cf(n_o, 0);
  for (int iter = 0; iter < 4096 && !e; iter++) {
    int t = 0;
    while (t < n_o && cf[t] == 4) {
      cf[t] = 0;
      t++;
    }
    if (t == n_o) {
      break;
    }
    cf[t]++;
    ModMatrix x = ZeroMatrix<T_residue>(k, k);
    for (int i = 0; i < n_o; i++) {
      for (int a = 0; a < k; a++) {
        for (int b = 0; b < k; b++) {
          x(a, b) = ResAdd(x(a, b), ResMul(static_cast<T_residue>(cf[i]),
                                           O0p[i](a, b), p), p);
        }
      }
    }
    // x^2 = tr x - nr, from the flattened matrices
    ModMatrix basis(2, k * k);
    basis.row(0) = Eigen::Map<Eigen::Matrix<T_residue, 1, Eigen::Dynamic>>(
        FlattenModMatrix(x).data(), k * k);
    basis.row(1) = Eigen::Map<Eigen::Matrix<T_residue, 1, Eigen::Dynamic>>(
        FlattenModMatrix(Id).data(), k * k);
    if (RankModPrime(basis, p) != 2) {
      continue;
    }
    std::optional<ModVector> sol =
        SolveRowsModP(basis, FlattenModMatrix(ModMatrixProduct(x, x, p)), p);
    if (!sol) {
      continue;
    }
    T_residue tr = (*sol)[0];
    T_residue nr = ResNeg((*sol)[1], p);
    T_residue disc = ResSub(ResMul(tr, tr, p), ResMul(T_residue(4), nr, p), p);
    if (disc == 0 || KroneckerResidue(disc, p) != 1) {
      continue;
    }
    T_residue s = SqrtModPrime(disc, p);
    T_residue inv2 = ResInv(T_residue(2), p);
    T_residue l1 = ResMul(ResAdd(tr, s, p), inv2, p);
    T_residue l2 = ResMul(ResSub(tr, s, p), inv2, p);
    T_residue inv12 = ResInv(ResSub(l1, l2, p), p);
    ModMatrix ee(k, k);
    for (int a = 0; a < k; a++) {
      for (int b = 0; b < k; b++) {
        T_residue v = ResSub(x(a, b), (a == b) ? l2 : T_residue(0), p);
        ee(a, b) = ResMul(v, inv12, p);
      }
    }
    if (ModMatrixProduct(ee, ee, p) == ee) {
      e = ee;
    }
  }
  if (!e) {
    return {};
  }
  // W = e X / p X, by the columns of e
  RowSpaceModP W = ComputeRowSpaceModP(ModMatrix(e->transpose()), p);
  int w = W.basis.rows();
  if (2 * w != k) {
    std::cerr << "ALB2_MASS: the idempotent has rank " << w << "\n";
    throw TerminalException{1};
  }
  std::vector<T_residue> res;
  for (auto const &g64 : image_gens) {
    ModMatrix g = ModMatrixOfModRowsI64(g64, k);
    for (int a = 0; a < k; a++) {
      for (int b = 0; b < k; b++) {
        g(a, b) %= p;
      }
    }
    std::vector<int64_t> M(w * w);
    for (int i = 0; i < w; i++) {
      ModVector col(k);
      for (int a = 0; a < k; a++) {
        T_residue sum = 0;
        for (int b = 0; b < k; b++) {
          sum = ResAdd(sum, ResMul(g(a, b), W.basis(i, b), p), p);
        }
        col[a] = sum;
      }
      ModVector cw = CoordinatesInRowSpaceModP(W, col);
      for (int j = 0; j < w; j++) {
        M[j * w + i] = ResToSigned(cw[j]);
      }
    }
    int64_t det = DeterminantMatModInt64(M, w, ResToSigned(p), true);
    res.push_back(ResModSigned(det, p));
  }
  return res;
}

/*
  The common part of the masses by the Tamagawa number, squared: tau^2
  times the squared L-product over the squared real volume, divided by the
  squared ratios rho_p = vol(K_p) / g_p at the bad primes, from the image
  orders n_p of the connected stabilizer.
 */
template <typename T> struct TamagawaLocalData {
  int64_t p;
  int c;
  // |K_{N^0, p} mod p^c|: the image order divided by [K_N : K_{N^0}]
  T n_over_idx;
  int v_disc_u;
  T g_p;
};

template <typename T>
T TamagawaMassSquared(T const &tau, T const &Lprod_over_pi, T const &Vinf_over_pi,
                      int dim_u, std::vector<TamagawaLocalData<T>> const &loc) {
  T num = tau * Lprod_over_pi;
  T res = num * num / (Vinf_over_pi * Vinf_over_pi);
  for (auto const &ld : loc) {
    T pT = ScalarFromInt64<T>(ld.p);
    // vol^2 = n^2 p^{-2 c dim} p^{-v}
    T vol2 = ld.n_over_idx;
    vol2 *= vol2;
    int ex = 2 * ld.c * dim_u + ld.v_disc_u;
    for (int i = 0; i < ex; i++) {
      vol2 /= pT;
    }
    T g = ld.g_p;
    res *= g * g / vol2;
  }
  return res;
}

template <typename T, typename Tint>
std::optional<T> AlbertTypeTwoOrthogonalMass(MyMatrix<T> const &B,
                                             std::vector<MyMatrix<T>> const &Dspan,
                                             int r,
                                             AlbertTypeTwoMassBounds const &bounds,
                                             std::ostream &os) {
  int k = B.rows();
  if (k != 4 * r || r < 2) {
    std::cerr << "ALB2_MASS: X of rank " << k << " for r=" << r << "\n";
    throw TerminalException{1};
  }
  int dim_u = r * (2 * r - 1);
  std::vector<MyMatrix<T>> O0 = QuatIntegralPoints<T, Tint>(Dspan, k);
  if (O0.size() != 4) {
    std::cerr << "ALB2_MASS: D has dimension " << O0.size() << "\n";
    throw TerminalException{1};
  }
  std::optional<std::vector<MyMatrix<T>>> opt_O =
      QuatMaximalOrder<T, Tint>(O0, r, T(1), bounds.quat, os);
  if (!opt_O) {
    return {};
  }
  T discD = QuatOrderDisc(*opt_O, r);
  std::set<int64_t> ram;
  {
    Tint dd = UniversalScalarConversion<Tint, T>(discD);
    if (dd > 1) {
      for (Tint const &q : FactorsInt(dd)) {
        ram.insert(UniversalScalarConversion<int64_t, Tint>(q));
      }
    }
  }
  if (ram.empty()) {
    return {};
  }
  std::vector<MyMatrix<T>> u = QuatIntegralCommutant<T, Tint>(B, O0, true);
  if (static_cast<int>(u.size()) != dim_u) {
    std::cerr << "ALB2_MASS: the Lie lattice has rank " << u.size()
              << ", not " << dim_u << "\n";
    throw TerminalException{1};
  }
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
  // the discriminant field, from an invertible skew-adjoint element
  int64_t Df = 0;
  for (int t = 0; t < static_cast<int>(u.size()) && Df == 0; t++) {
    for (int s = -1; s <= 1 && Df == 0; s += 2) {
      MyMatrix<T> Y = u[t];
      for (int t2 = 0; t2 < static_cast<int>(u.size()); t2++) {
        if (t2 != t) {
          Y += ScalarFromInt64<T>(s * (t2 + 1)) * u[t2];
        }
      }
      T det = DeterminantMat(Y);
      if (det == 0) {
        continue;
      }
      std::optional<T> nrd = RationalSquareRoot<T, Tint>(det);
      if (!nrd) {
        std::cerr << "ALB2_MASS: det_X of a skew element is not a square\n";
        throw TerminalException{1};
      }
      T discs = (r % 2 == 0) ? *nrd : T(-*nrd);
      Tint num = GetNumerator_z(discs) * GetDenominator_z(discs);
      Df = FundamentalDiscriminantOf<Tint>(num);
    }
  }
  if (Df == 0) {
    std::cerr << "ALB2_MASS: no invertible skew-adjoint element found\n";
    throw TerminalException{1};
  }
  if ((Df > 0) != (r % 2 == 0)) {
    std::cerr << "ALB2_MASS: the discriminant " << Df << " has the wrong sign "
              << "for r=" << r << "\n";
    throw TerminalException{1};
  }
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
  for (T const &x : {DeterminantMat(B), disc_u, QuatOrderDisc(O0, r), discD,
                     ScalarFromInt64<T>(Df)}) {
    for (int64_t q : prime_set(x)) {
      primes.insert(q);
    }
  }
  auto g_factor = [&](int64_t p) -> T {
    T pinv = T(1) / ScalarFromInt64<T>(p);
    T res(1), pw(1);
    for (int i = 1; i < r; i++) {
      pw *= pinv * pinv;
      res *= T(1) - pw;
    }
    T pr(1);
    for (int i = 0; i < r; i++) {
      pr *= pinv;
    }
    int ch = kronecker_symbol<int64_t>(Df, p);
    return res * (T(1) - ScalarFromInt64<T>(ch) * pr);
  };
  std::vector<TamagawaLocalData<T>> loc;
  int nB = 0;
  auto local = [&](int64_t p) -> std::optional<std::pair<TamagawaLocalData<T>, int>> {
    int cc = (p == 2) ? 2 : 1;
    std::vector<ModRowsI64> img;
    std::optional<std::pair<Tint, Tint>> o = QuatLocalImageOrdersGenerators<T, Tint>(
        B, O0, Lam, Forms, p, cc, bounds.quat.image, &img, os);
    if (!o) {
      return {};
    }
    int idx = 1;
    bool in_ram = ram.count(p) == 1;
    if (!in_ram) {
      if (p == 2) {
        return {};
      }
      std::optional<std::vector<T_residue>> nrds =
          ReducedNormsModP(O0, img, k, static_cast<T_small_prime>(p));
      if (!nrds) {
        return {};
      }
      for (T_residue x : *nrds) {
        if (x == static_cast<T_residue>(p - 1)) {
          idx = 2;
        } else if (x != 1) {
          std::cerr << "ALB2_MASS: a reduced norm " << x << " modulo " << p
                    << "\n";
          throw TerminalException{1};
        }
      }
    }
    TamagawaLocalData<T> ld;
    ld.p = p;
    ld.c = cc;
    ld.n_over_idx = UniversalScalarConversion<T, Tint>(o->first) /
                    ScalarFromInt64<T>(idx);
    ld.v_disc_u = RationalValuation(disc_u, p);
    ld.g_p = g_factor(p);
#ifdef DEBUG_ALBERT_TYPE_TWO_MASS
    os << "ALB2_MASS: r=" << r << " p=" << p << (in_ram ? " (ramified in D)" : "")
       << " |K_N mod p^c|=" << o->first << " [K_N:K_N0]=" << idx
       << " v(disc u)=" << ld.v_disc_u << " g_p=" << ld.g_p << "\n";
#endif
    return std::make_pair(ld, (!in_ram && idx == 1) ? 1 : 0);
  };
  for (int64_t p : primes) {
    auto res = local(p);
    if (!res) {
      return {};
    }
    loc.push_back(res->first);
    nB += res->second;
  }
  T Lprod = ScalarFromInt64<T>(1);
  for (int i = 1; i < r; i++) {
    Lprod *= ZetaEvenOverPi<T>(i);
  }
  T L2 = LValueOverPiSquared<T>(r, Df);
  T Vinf(1);
  for (int i = 0; i < dim_u; i++) {
    Vinf *= T(2);
  }
  for (int j = 1; j <= 2 * r - 1; j++) {
    Vinf *= SphereVolumeRational<T>(j);
  }
  // the square of mass_{N^0}: L(r, chi)^2 enters squared already
  T m2 = TamagawaMassSquared<T>(T(2), Lprod, Vinf, dim_u, loc) * L2;
  std::optional<T> m0 = RationalSquareRoot<T, Tint>(m2);
  if (!m0) {
    std::cerr << "ALB2_MASS: the squared mass " << m2 << " is not a square\n";
    throw TerminalException{1};
  }
  T mass = *m0;
  for (int i = 0; i < nB; i++) {
    mass *= T(2);
  }
#ifdef DEBUG_ALBERT_TYPE_TWO_MASS
  os << "ALB2_MASS: r=" << r << " Df=" << Df << " mass(N0)=" << *m0
     << " |B|=" << nB << " mass=" << mass << "\n";
#endif
#ifdef SANITY_CHECK_ALBERT_TYPE_TWO_MASS
  {
    int64_t p = 2;
    int n_checked = 0;
    while (n_checked < 2) {
      p++;
      if (primes.count(p) == 1 || !IsPrime(p)) {
        continue;
      }
      auto res = local(p);
      if (res) {
        TamagawaLocalData<T> const &ld = res->first;
        T vol = ld.n_over_idx;
        for (int i = 0; i < ld.c * dim_u + ld.v_disc_u / 2; i++) {
          vol /= ScalarFromInt64<T>(p);
        }
        if (ld.v_disc_u != 0 || vol != ld.g_p || res->second != 0) {
          std::cerr << "ALB2_MASS: the good prime " << p << " has vol " << vol
                    << " for g_p " << ld.g_p << "\n";
          throw TerminalException{1};
        }
      }
      n_checked++;
    }
  }
#endif
  return mass;
}

/*
  The same machinery on Albert type III of rank one over Q (U = D^1 =
  SL_1(D), tau = 1, U(R) = SU(2) on X_R = H^{k/4}, metric -tr_X = k |Y|^2
  on the pure quaternions, so vol = 2 pi^2 k^{3/2}; g_p = 1 - p^{-2}), as
  a check of the measures against QuaternionicComponentMass.
 */
template <typename T, typename Tint>
std::optional<T> TypeThreeRankOneTamagawaMass(MyMatrix<T> const &B,
                                              std::vector<MyMatrix<T>> const &Dspan,
                                              AlbertTypeTwoMassBounds const &bounds,
                                              std::ostream &os) {
  int k = B.rows();
  std::vector<MyMatrix<T>> O0 = QuatIntegralPoints<T, Tint>(Dspan, k);
  std::vector<MyMatrix<T>> u = QuatIntegralCommutant<T, Tint>(B, O0, true);
  if (u.size() != 3) {
    return {};
  }
  T disc_u = QuatTraceDiscriminant(u);
  std::vector<MyMatrix<T>> Lam = QuatIntegralCommutant<T, Tint>(B, O0, false);
  std::vector<MyMatrix<T>> L;
  for (auto const &w : O0) {
    L.push_back(MyMatrix<T>(w.transpose() * B));
  }
  std::vector<MyMatrix<T>> Forms = QuatIntegralPoints<T, Tint>(L, k);
  std::set<int64_t> primes{2};
  for (T const &x : {DeterminantMat(B), disc_u, QuatOrderDisc(O0, k / 4)}) {
    for (Tint xi : {Tint(T_abs(GetNumerator_z(x))), Tint(GetDenominator_z(x))}) {
      if (xi > 1) {
        for (Tint const &q : FactorsInt(xi)) {
          primes.insert(UniversalScalarConversion<int64_t, Tint>(q));
        }
      }
    }
  }
  std::vector<TamagawaLocalData<T>> loc;
  for (int64_t p : primes) {
    int cc = (p == 2) ? 2 : 1;
    std::optional<std::pair<Tint, Tint>> o = QuatLocalImageOrders<T, Tint>(
        B, O0, Lam, Forms, p, cc, bounds.quat.image, os);
    if (!o) {
      return {};
    }
    TamagawaLocalData<T> ld;
    ld.p = p;
    ld.c = cc;
    ld.n_over_idx = UniversalScalarConversion<T, Tint>(o->first);
    ld.v_disc_u = RationalValuation(disc_u, p);
    ld.g_p = T(1) - T(1) / ScalarFromInt64<T>(p * p);
    loc.push_back(ld);
  }
  // vol^2 = 4 pi^4 k^3
  T Vinf2 = T(4) * ScalarFromInt64<T>(k) * ScalarFromInt64<T>(k) *
            ScalarFromInt64<T>(k);
  T m2 = TamagawaMassSquared<T>(T(1), ZetaEvenOverPi<T>(1), T(1), 3, loc) / Vinf2;
  return RationalSquareRoot<T, Tint>(m2);
}

// clang-format off
#endif  // SRC_GENERA_ALBERTTYPETWOMASS_H_
// clang-format on

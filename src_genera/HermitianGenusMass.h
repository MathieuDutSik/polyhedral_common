// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_HERMITIANGENUSMASS_H_
#define SRC_GENERA_HERMITIANGENUSMASS_H_

// clang-format off
#include "HermitianGenusSeed.h"
#include "ClassicalGroupOrders.h"
#include "MassQuadraticHeckeL.h"
#include <map>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_HERMITIAN_GENUS_MASS
#endif

/*
  The local factors of the mass of a genus of positive definite Hermitian
  lattices over the order O_E of any CM field (CMOrder.h), as local_factor
  of Hecke.jl (QuadForm/Herm/Mass.jl, after Gan-Hanke-Yu and Cho), from the
  local symbols of HermitianGenusSymbol.h. The mass is the standard mass
  times the product of the local factors at the primes of the symbol.

  At a prime p of K_0 that is not dyadic and ramified, with q = N(p) and the
  Jordan blocks (s_i, m_i, d_i) of L_p, m = sum m_i, r_i = sum_{j > i} m_j:

    q^N f,  f = G(m) / prod_i G_i(m_i),

  G the group order of Hecke.jl (ClassicalGroupOrders.h): GL at a split
  prime, U at an inert one, Sp(2 floor(m / 2)) at a ramified one, with G_i
  = GL, U, and at a ramified prime Sp(m_i) for an odd s_i, O+ or O- for an
  even s_i as the discriminant (-1)^{m_i (m_i - 1) / 2} det of the block is
  a norm or not; and N = sum_i (- c_i - d s_i m_i r_i + d s_i m_i m / 2),
  d = 2 and c_i = s_i m_i^2 at an unramified prime, d = 1 and c_i =
  floor(s_i / 2) m_i^2 (plus (m_i + 1) floor(m_i / 2) for an odd s_i) at a
  ramified one, plus m / 2 at a ramified prime when m is even. The factor
  is 1 at an unramified prime with a single block.

  At a dyadic ramified prime, the chain method of Hecke.jl: with the
  lattices whose norm has the p-valuation at least that of L (a-maximal for
  a = n(L)_p, which is the maximal integral case of Hecke.jl after its
  rescaling), the chain L = L_0 in L_1 in ... in L_k of overlattices M +
  P^{-1} x up to an a-maximal L_k gives the factor
  f_max prod_i #{M in L_{i+1} of index N(P), M ~ L_i} / #{N over L_i of
  index N(P), N ~ L_{i+1}}, the lattices compared by their local symbols
  (complete invariants), f_max the factor of a maximal lattice: 1/2 for an
  odd m, and when the discriminant of the space is not a norm (q^m - 1) /
  (2 q + 2) (ramified), 1 otherwise.
 */

// The norm of the prime of K_0 under P: N(P) when p splits or ramifies,
// its square root when p is inert.
template <typename Tint> int64_t HermK0Norm(HermIdealPrime<Tint> const &P) {
  int f = (P.type == HermIdealType::Inert) ? P.f / 2 : P.f;
  int64_t q = 1;
  for (int i = 0; i < f; i++) {
    q *= P.ell;
  }
  return q;
}

template <typename T, typename Tint, typename Tgroup>
T HermLocalFactor(CMOrder<T, Tint> const &ord,
                  RealifiedRing<T, Tint> const &ring, MyMatrix<T> const &A,
                  std::vector<HermIdealPrime<Tint>> const &l_prime, size_t i_P) {
  HermIdealPrime<Tint> const &P = l_prime[i_P];
  HermLatticeData<T, Tint> L = HermLatticeOfBasis(ord, ring, A);
  int m = L.n() / ord.deg;
  T q = UniversalScalarConversion<T, int64_t>(HermK0Norm(P));
  bool ram = P.type == HermIdealType::Ramified;
  bool split = P.type == HermIdealType::Split;
  auto disc_norm_class = [&]() -> int {
    // the discriminant (-1)^{m (m - 1) / 2} det of the space
    MyMatrix<T> H(m, m);
    std::vector<MyVector<Tint>> l_v = HermLocalBasis(L, P.ell);
    (void)H;
    // det h on the local basis, in K_0
    int d = ord.deg;
    MyMatrix<T> Bt = ord.basis_theta;
    MyMatrix<T> Bt_inv = Inverse(Bt);
    NumberField<T> const &K = ord.field;
    SmallGram<T> G(m, std::vector<ElemE<T>>(m));
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++) {
        MyVector<T> th = Bt.transpose() * L.herm(l_v[i], l_v[j]);
        ElemE<T> x(d, T(0));
        for (int k = 0; k < d; k++) {
          x[k] = th(k);
        }
        G[i][j] = x;
      }
    }
    // the determinant by elimination over E
    ElemE<T> det = EOne(K);
    for (int c = 0; c < m; c++) {
      int piv = -1;
      for (int r = c; r < m; r++) {
        if (!G[r][c].is_zero()) {
          piv = r;
          break;
        }
      }
      if (piv != c) {
        std::swap(G[piv], G[c]);
        det = -det;
      }
      det = EMul(det, G[c][c], K);
      ElemE<T> inv = EInv(G[c][c], K);
      for (int r = c + 1; r < m; r++) {
        ElemE<T> f = EMul(G[r][c], inv, K);
        for (int k = c; k < m; k++) {
          G[r][k] = G[r][k] - EMul(f, G[c][k], K);
        }
      }
    }
    if ((m * (m - 1) / 2) % 2 == 1) {
      det = -det;
    }
    MyVector<T> th(d);
    for (int k = 0; k < d; k++) {
      th(k) = det[k];
    }
    return HermLocalNormClass(ord, l_prime, i_P, MyVector<T>(Bt_inv.transpose() * th));
  };
  if (ram && P.ell == 2) {
    // the chain method
    auto f_max = [&]() -> T {
      if (m % 2 == 1) {
        return T(1) / T(2);
      }
      if (disc_norm_class() == -1) {
        T qm(1);
        for (int i = 0; i < m; i++) {
          qm *= q;
        }
        return (qm - 1) / (2 * q + 2);
      }
      return T(1);
    };
    int a = HermNormValAt(L, P);
    int diff_val = HermDifferentVal(ord, l_prime, i_P);
    auto symbol = [&](MyMatrix<T> const &B) -> HermLocalSymbol {
      return HermLocalSymbolOfLattice(HermLatticeOfBasis(ord, ring, B),
                                      l_prime, i_P);
    };
    auto same = [&](HermLocalSymbol const &x, HermLocalSymbol const &y) -> bool {
      return HermLocalSymbolEquivalent(x, y, true, true, diff_val);
    };
    MyMatrix<T> P_inv = HermIdealInverse(ord, P);
    auto admissible = [&](MyMatrix<T> const &B) -> bool {
      MyMatrix<T> G = B * UniversalMatrixConversion<T, Tint>(ring.B) *
                      B.transpose();
      if (!IsIntegralMatrix(G)) {
        return false;
      }
      return HermNormValAt(HermLatticeOfBasis(ord, ring, B), P) >= a;
    };
    // the overlattices M + P^{-1} x of the lattice of basis B
    auto overlattices = [&](MyMatrix<T> const &B) -> std::vector<MyMatrix<T>> {
      HermLatticeData<T, Tint> LB = HermLatticeOfBasis(ord, ring, B);
      HermLineSpace ls = HermGetLineSpace(LB, P);
      std::vector<MyMatrix<T>> res;
      for (size_t idx = 0; idx < ls.n_line(); idx++) {
        std::vector<int64_t> z = ls.line_vector(idx);
        MyVector<Tint> x(LB.n());
        for (int k = 0; k < LB.n(); k++) {
          x(k) = UniversalScalarConversion<Tint, int64_t>(z[k]);
        }
        res.push_back(MyMatrix<T>(HermAddIdealTimes(LB, P_inv, x) * B));
      }
      return res;
    };
    auto sublattices = [&](MyMatrix<T> const &B) -> std::vector<MyMatrix<T>> {
      HermLatticeData<T, Tint> LB = HermLatticeOfBasis(ord, ring, B);
      HermLineSpace ls = HermGetLineSpace(LB, P);
      std::vector<MyMatrix<T>> res;
      for (size_t idx = 0; idx < ls.n_line(); idx++) {
        std::vector<int> c =
            ProjectiveLineOfIndex(ls.m, static_cast<int>(ls.F.q), idx);
        res.push_back(MyMatrix<T>(HermHyperplaneSublattice(LB, ls, c) * B));
      }
      return res;
    };
    std::vector<MyMatrix<T>> chain{A};
    while (true) {
      std::optional<MyMatrix<T>> next;
      for (auto const &Bn : overlattices(chain.back())) {
        if (admissible(Bn)) {
          next = Bn;
          break;
        }
      }
      if (!next) {
        break;
      }
      chain.push_back(*next);
    }
    T f = f_max();
    for (size_t i = 0; i + 1 < chain.size(); i++) {
      HermLocalSymbol si = symbol(chain[i]), sn = symbol(chain[i + 1]);
      int64_t n_sub = 0, n_over = 0;
      for (auto const &B : sublattices(chain[i + 1])) {
        if (same(symbol(B), si)) {
          n_sub++;
        }
      }
      for (auto const &B : overlattices(chain[i])) {
        if (admissible(B) && same(symbol(B), sn)) {
          n_over++;
        }
      }
      f *= UniversalScalarConversion<T, int64_t>(n_sub);
      f /= UniversalScalarConversion<T, int64_t>(n_over);
    }
    return f;
  }
  HermLocalSymbol sym = HermLocalSymbolOfLattice(L, l_prime, i_P);
  if (sym.blocks.size() == 1 && !ram) {
    return T(1);
  }
  int type = split ? 0 : 1;
  T f = ram ? HeckeGroupOrderSp<T>(2 * (m / 2), q)
            : HermitianGroupOrderNormalized<T>(type, m, q);
  int d = ram ? 1 : 2;
  int64_t N = 0;
  int neg_one = ram ? HermLocalNormClass(ord, l_prime, i_P,
                                         MyVector<T>([&]() {
                                           MyVector<T> v = ZeroVector<T>(ord.deg);
                                           v(0) = -1;
                                           return v;
                                         }()))
                    : 1;
  size_t t = sym.blocks.size();
  for (size_t i = 0; i < t; i++) {
    int64_t s = sym.blocks[i][0];
    int64_t mi = sym.blocks[i][1];
    int64_t ri = 0;
    for (size_t j = i + 1; j < t; j++) {
      ri += sym.blocks[j][1];
    }
    if (ram) {
      int64_t fl = (s >= 0) ? s / 2 : -((-s + 1) / 2);
      N -= fl * mi * mi;
      if (s % 2 != 0) {
        N -= (mi + 1) * (mi / 2);
        f /= HeckeGroupOrderSp<T>(static_cast<int>(mi), q);
      } else {
        int cls = sym.blocks[i][2];
        if ((mi * (mi - 1) / 2) % 2 == 1) {
          cls *= neg_one;
        }
        f /= HeckeGroupOrderO<T>(static_cast<int>(mi), cls, q);
      }
    } else {
      N -= s * mi * mi;
      f /= HermitianGroupOrderNormalized<T>(type, static_cast<int>(mi), q);
    }
    N -= d * s * mi * ri;
    N += (d * s * mi * m) / 2;
  }
  if (ram && m % 2 == 0) {
    N += m / 2;
  }
  T qN(1);
  for (int64_t i = 0; i < std::abs(N); i++) {
    qN *= q;
  }
  if (N < 0) {
    qN = T(1) / qN;
  }
  return qN * f;
}

/*
  The product of the local factors of the lattice of basis A at the primes
  of its genus symbol (the primes of HermGenusSymbolOfLattice).
 */
template <typename T, typename Tint, typename Tgroup>
T HermLocalMass(CMOrder<T, Tint> const &ord, RealifiedRing<T, Tint> const &ring,
                MyMatrix<T> const &A, std::ostream &os) {
  HermLatticeData<T, Tint> L = HermLatticeOfBasis(ord, ring, A);
  Tint det = T_abs(DeterminantMat(L.G));
  std::set<int64_t> l_ell{2};
  Tint rest = det;
  for (int64_t p = 2;; p++) {
    Tint p_T = UniversalScalarConversion<Tint, int64_t>(p);
    if (p_T * p_T > rest) {
      break;
    }
    while (rest % p_T == 0) {
      l_ell.insert(p);
      rest /= p_T;
    }
  }
  if (rest > 1) {
    l_ell.insert(UniversalScalarConversion<int64_t, Tint>(rest));
  }
  T res(1);
  for (int64_t ell : l_ell) {
    std::vector<HermIdealPrime<Tint>> l_prime = HermIdealPrimesAbove(ord, ell);
    for (size_t j = 0; j < l_prime.size(); j++) {
      if (HermK0Representative(l_prime, j) != j) {
        continue;
      }
      T lf = HermLocalFactor<T, Tint, Tgroup>(ord, ring, A, l_prime, j);
#ifdef DEBUG_HERMITIAN_GENUS_MASS
      os << "HERMMASS: local factor " << lf << " at a prime above " << ell
         << "\n";
#endif
      res *= lf;
    }
  }
  (void)os;
  return res;
}

/*
  The value zeta_{K_0}(1 - 2k) of the Dedekind zeta function of the real
  subfield K_0 of E (degree n >= 2), exactly, by Siegel's formula
  (Berechnung von Zetafunktionen an ganzzahligen Stellen, 1969). The
  restriction to the diagonal of the Hilbert Eisenstein series of parallel
  weight 2k,

    F = zeta_{K_0}(1 - 2k) / 2^n + sum_{j >= 1} c_j q^j,
    c_j = sum_{nu in d^{-1}, nu >> 0, Tr nu = j} sigma_{2k-1}(nu d),

  is in M_h(SL_2(Z)), h = 2 k n. For T in M^!_{2-h}, F T has weight 2 and no
  constant term; with T = E_w Delta^{-N}, w = 12 N + 2 - h in {0, 4, 6, 8,
  10, 14}, scaled to the constant term 1, this gives c_0 = - sum_{j=1}^{N}
  t_{-j} c_j. The nu are the short vectors of the inverse different d^{-1}
  (the dual of O_{K_0} for the trace) with Tr(nu^2) <= j^2, kept when
  totally positive (Tr(nu y^2) > 0) of trace j; sigma_{2k-1}(a) = prod_p
  sum_{i <= v_p(a)} N(p)^{i (2k - 1)}, the valuations of a = nu d at the
  primes of K_0 being those at the primes of E above them.
 */
template <typename T, typename Tint>
T DedekindZetaTotallyRealExact(CMOrder<T, Tint> const &ord, int k,
                               std::ostream &os) {
  int d = ord.deg;
  int n = d / 2;
  if (n == 1) {
    // K_0 = Q: zeta(1 - 2k) = -B_{2k} / (2k), the Bernoulli numbers by
    // sum_{j <= m} binom(m + 1, j) B_j = 0
    std::vector<T> Bn{T(1)};
    for (int mm = 1; mm <= 2 * k; mm++) {
      T sum(0);
      T binom(1);
      for (int j = 0; j < mm; j++) {
        sum += binom * Bn[j];
        binom = binom * T(mm + 1 - j) / T(j + 1);
      }
      Bn.push_back(-sum / T(mm + 1));
    }
    return -Bn[2 * k] / T(2 * k);
  }
  int h = 2 * k * n;
  // T = E_w Delta^{-N}
  int N = std::max(0, (h - 2 + 11) / 12);
  auto w_of = [&](int NN) -> int { return 12 * NN + 2 - h; };
  std::set<int> good_w{0, 4, 6, 8, 10, 14};
  while (!good_w.count(w_of(N))) {
    N++;
  }
  int w = w_of(N);
  int prec = N + 1;
  // the q-expansions up to q^N, with integer coefficients
  auto mul = [&](std::vector<Tint> const &x, std::vector<Tint> const &y) {
    std::vector<Tint> z(prec, Tint(0));
    for (int i = 0; i < prec; i++) {
      for (int j = 0; i + j < prec; j++) {
        z[i + j] += x[i] * y[j];
      }
    }
    return z;
  };
  auto eis = [&](int64_t c, int e) -> std::vector<Tint> {
    // 1 + c sum sigma_e(m) q^m
    std::vector<Tint> z(prec, Tint(0));
    z[0] = 1;
    for (int m = 1; m < prec; m++) {
      Tint sg(0);
      for (int dv = 1; dv <= m; dv++) {
        if (m % dv == 0) {
          Tint p(1);
          for (int i = 0; i < e; i++) {
            p *= Tint(dv);
          }
          sg += p;
        }
      }
      z[m] = UniversalScalarConversion<Tint, int64_t>(c) * sg;
    }
    return z;
  };
  std::vector<Tint> E4 = eis(240, 3), E6 = eis(-504, 5);
  std::vector<Tint> Ew(prec, Tint(0));
  Ew[0] = 1;
  if (w == 4) {
    Ew = E4;
  } else if (w == 6) {
    Ew = E6;
  } else if (w == 8) {
    Ew = mul(E4, E4);
  } else if (w == 10) {
    Ew = mul(E4, E6);
  } else if (w == 14) {
    Ew = mul(mul(E4, E4), E6);
  }
  // prod_{m >= 1} (1 - q^m)^{-24 N}: the inverse of prod (1 - q^m)^{24 N}
  std::vector<Tint> P(prec, Tint(0));
  P[0] = 1;
  for (int m = 1; m < prec; m++) {
    for (int t = 0; t < 24 * N; t++) {
      // multiply by (1 - q^m)
      for (int i = prec - 1; i >= m; i--) {
        P[i] -= P[i - m];
      }
    }
  }
  std::vector<Tint> Pinv(prec, Tint(0));
  Pinv[0] = 1;
  for (int i = 1; i < prec; i++) {
    Tint sum(0);
    for (int j = 1; j <= i; j++) {
      sum += P[j] * Pinv[i - j];
    }
    Pinv[i] = -sum;
  }
  std::vector<Tint> G = mul(Ew, Pinv);
  // T = q^{-N} G: t_{-j} = G[N - j], t_0 = G[N]
  if (G[N] == 0) {
    std::cerr << "HERMMASS: the form of weight 2 - h has no constant term\n";
    throw TerminalException{1};
  }
  // O_{K_0}: the fixed points of the conjugation on O_E
  MyMatrix<Tint> B0 = NullspaceIntMat(MyMatrix<Tint>(ord.conj - IdentityMat<Tint>(d)));
  MyMatrix<Tint> Q = ord.trace_norm_form();
  // Tr_{K_0}(x y) = Tr_E(x y) / 2 on K_0
  MyMatrix<T> TK = UniversalMatrixConversion<T, Tint>(MyMatrix<Tint>(B0 * Q * B0.transpose())) / T(2);
  MyMatrix<T> TKinv = Inverse(TK);
  MyMatrix<T> B0_T = UniversalMatrixConversion<T, Tint>(B0);
  // the dual basis, on e
  MyMatrix<T> Dual = TKinv * B0_T;
  // the primes of K_0 of the different: v_p(D) = - min v_p(dual basis)
  auto elem_of = [&](MyVector<T> const &y) -> MyVector<T> {
    return Dual.transpose() * y;
  };
  T c0(0);
  // the c_j, j <= N
  std::vector<T> cj(N + 1, T(0));
  if (N >= 1) {
    CVPSolver<T, Tint> solver(TKinv, os);
    T bound = T(N) * T(N);
    std::vector<MyVector<Tint>> l_v = solver.at_most_norm_vectors(bound);
    std::map<int64_t, std::vector<HermIdealPrime<Tint>>> cache;
    auto primes = [&](int64_t ell) -> std::vector<HermIdealPrime<Tint>> const & {
      auto it = cache.find(ell);
      if (it == cache.end()) {
        it = cache.emplace(ell, HermIdealPrimesAbove(ord, ell)).first;
      }
      return it->second;
    };
    for (auto const &v0 : l_v) {
      for (int sgn : {1, -1}) {
        MyVector<T> y = UniversalVectorConversion<T, Tint>(v0) * T(sgn);
        MyVector<T> nu = elem_of(y);
        T tr(0);
        for (int a = 0; a < d; a++) {
          tr += nu(a) * UniversalScalarConversion<T, Tint>(ord.trace[a]);
        }
        tr /= T(2);
        if (tr <= 0 || tr > T(N) || !IsInteger(tr)) {
          continue;
        }
        // totally positive: Tr(nu y^2) > 0 on O_{K_0}
        MyMatrix<T> Mnu = ZeroMatrix<T>(d, d);
        for (int a = 0; a < d; a++) {
          if (nu(a) != 0) {
            Mnu += nu(a) * UniversalMatrixConversion<T, Tint>(ord.mult[a]);
          }
        }
        MyMatrix<T> Tnu = B0_T * Mnu * UniversalMatrixConversion<T, Tint>(Q) * B0_T.transpose();
        if (!IsPositiveDefinite(Tnu, os)) {
          continue;
        }
        int j = UniversalScalarConversion<int, T>(tr);
        // N_{K_0}(nu)^2 = N_E(nu); a = nu d has the norm N(nu) |d_{K_0}|
        T NE = DeterminantMat(Mnu);
        T dK = DeterminantMat(TK);
        T Na2 = NE * dK * dK;
        // the primes of a: the divisors of the numerator of N(a)^2
        Tint num = GetNumerator_z(Na2);
        std::vector<int64_t> l_ell;
        for (int64_t p = 2;; p++) {
          Tint p_T = UniversalScalarConversion<Tint, int64_t>(p);
          if (p_T * p_T > num) {
            break;
          }
          if (num % p_T == 0) {
            l_ell.push_back(p);
            while (num % p_T == 0) {
              num /= p_T;
            }
          }
        }
        if (num > 1) {
          l_ell.push_back(UniversalScalarConversion<int64_t, Tint>(num));
        }
        T sigma(1);
        for (int64_t ell : l_ell) {
          std::vector<HermIdealPrime<Tint>> const &lp = primes(ell);
          for (size_t i = 0; i < lp.size(); i++) {
            if (HermK0Representative(lp, i) != i) {
              continue;
            }
            int e_rel = (lp[i].type == HermIdealType::Ramified) ? 2 : 1;
            int v_nu = HermIdealVal(ord, lp[i], nu) / e_rel;
            int v_D = 1 << 28;
            for (int r = 0; r < Dual.rows(); r++) {
              v_D = std::min(v_D, HermIdealVal(ord, lp[i], MyVector<T>(GetMatrixRow(Dual, r))) / e_rel);
            }
            int v = v_nu - v_D;
            T Np = UniversalScalarConversion<T, int64_t>(HermK0Norm(lp[i]));
            T term(0), pw(1), step(1);
            for (int t = 0; t < 2 * k - 1; t++) {
              step *= Np;
            }
            for (int t = 0; t <= v; t++) {
              term += pw;
              pw *= step;
            }
            sigma *= term;
          }
        }
        cj[j] += sigma;
      }
    }
  }
  for (int j = 1; j <= N; j++) {
    c0 -= UniversalScalarConversion<T, Tint>(G[N - j]) * cj[j];
  }
  c0 /= UniversalScalarConversion<T, Tint>(G[N]);
  T res = c0;
  for (int i = 0; i < n; i++) {
    res *= T(2);
  }
  return res;
}

// The minimal polynomial (integer, by increasing degree) of the element x
// of O_E (coordinates on e), from its powers on the basis e.
template <typename T, typename Tint>
std::vector<int64_t> HermMinimalPolynomial(CMOrder<T, Tint> const &ord,
                                           MyVector<Tint> const &x) {
  int d = ord.deg;
  std::vector<MyVector<Tint>> pows;
  MyVector<Tint> pw = ZeroVector<Tint>(d);
  pw(0) = 1;
  for (int k = 0; k <= d; k++) {
    pows.push_back(pw);
    MyMatrix<T> M(k + 1, d);
    for (int i = 0; i <= k; i++) {
      for (int j = 0; j < d; j++) {
        M(i, j) = UniversalScalarConversion<T, Tint>(pows[i](j));
      }
    }
    MyMatrix<T> NSP = NullspaceMat(M);
    if (NSP.rows() > 0) {
      std::vector<int64_t> res;
      T lead = NSP(0, k);
      for (int i = 0; i <= k; i++) {
        T c = NSP(0, i) / lead;
        if (!IsInteger(c)) {
          std::cerr << "HERMMASS: the minimal polynomial is not integral\n";
          throw TerminalException{1};
        }
        res.push_back(UniversalScalarConversion<int64_t, T>(c));
      }
      return res;
    }
    pw = ord.product(pw, x);
  }
  std::cerr << "HERMMASS: no minimal polynomial found\n";
  throw TerminalException{1};
}

// The number of the roots modulo the prime q of the integer polynomial f.
inline int HermRootCountModQ(std::vector<int64_t> const &f, int64_t q) {
  T_small_prime qu = static_cast<T_small_prime>(q);
  PolyQ fq = PolyReduceSigned(f, qu);
  if (fq.degree() <= 0) {
    return 0;
  }
  PolyQ x = PolyQ::X(qu);
  PolyQ xq = PolyPowModQ(x, static_cast<uint64_t>(q), fq);
  PolyQ g = PolyGcdQ(xq - x, fq);
  return std::max(0, g.degree());
}

/*
  The mass of the genus of the lattice of basis A over any CM field: the
  standard mass 2^{1 - n m} prod_{i even <= m} zeta_{K_0}(1 - i) prod_{i odd
  <= m} L(1 - i, chi_{E/K_0}) times the local factors (HermLocalMass). The
  zeta values are exact (DedekindZetaTotallyRealExact); the L-values are
  computed numerically (QuadHeckeLNegative of MassQuadraticHeckeL.h, by
  Dokchitser's method, gamma factor Gamma_R(s + 1)^n, conductor A^2 =
  |d_E| / |d_{K_0}|) from the primes of HermIdealPrimesAbove, the
  coefficient c_1(q) at a good prime being #roots(Phi_E mod q) -
  #roots(Phi_{K_0} mod q) (a degree one prime of K_0 splits into two of
  degree one, or stays inert with one of degree two). The mass is made
  exact as in Hecke.jl: every U(L) is a finite subgroup of GL_{d m}(Z), so
  its order divides Minkowski's bound M(d m) and b = M(d m) mass is an
  integer; the product is evaluated in Arb balls (with flint) until one
  integer only is in the ball, or else rounded at two precisions that must
  agree.
 */
template <typename T, typename Tint, typename Tgroup>
T HermMassAnyField(CMOrder<T, Tint> const &ord, RealifiedRing<T, Tint> const &ring,
                   MyMatrix<T> const &A, std::ostream &os) {
  int d = ord.deg;
  int n = d / 2;
  int m = ring.m;
  T local = HermLocalMass<T, Tint, Tgroup>(ord, ring, A, os);
  // the exact part
  T exact(1);
  for (int i = 0; i < n * m - 1; i++) {
    exact /= T(2);
  }
  for (int i = 2; i <= m; i += 2) {
    exact *= DedekindZetaTotallyRealExact<T, Tint>(ord, i / 2, os);
  }
  // Minkowski's bound for the finite subgroups of GL_N(Z)
  int Nmink = d * m;
  Tint M(1);
  for (int64_t p = 2; p <= Nmink + 1; p++) {
    bool is_prime = true;
    for (int64_t q = 2; q * q <= p; q++) {
      if (p % q == 0) {
        is_prime = false;
      }
    }
    if (!is_prime) {
      continue;
    }
    int64_t e = 0;
    for (int64_t pj = 1; Nmink / (pj * (p - 1)) > 0; pj *= p) {
      e += Nmink / (pj * (p - 1));
    }
    for (int64_t i = 0; i < e; i++) {
      M *= UniversalScalarConversion<Tint, int64_t>(p);
    }
  }
  T factor = exact * local * UniversalScalarConversion<T, Tint>(M);
  // the character
  MyVector<Tint> theta = ord.theta;
  std::vector<int64_t> PhiE = ord.field.Phi;
  std::optional<std::vector<int64_t>> PhiK0;
  for (int t = 0; t < 4 && !PhiK0; t++) {
    MyVector<Tint> x;
    if (t == 0) {
      x = theta + ord.conjugate(theta);
    } else if (t == 1) {
      x = ord.product(theta, ord.conjugate(theta));
    } else {
      MyVector<Tint> y = theta;
      for (int i = 1; i < t; i++) {
        y = ord.product(y, theta);
      }
      x = y + ord.conjugate(y) + ord.product(theta, ord.conjugate(theta));
    }
    std::vector<int64_t> f = HermMinimalPolynomial(ord, x);
    if (static_cast<int>(f.size()) - 1 == n) {
      PhiK0 = f;
    }
  }
  if (!PhiK0) {
    std::cerr << "HERMMASS: no generator of K_0 found\n";
    throw TerminalException{1};
  }
  auto disc_of = [](std::vector<int64_t> const &f) -> Tint {
    MyMatrix<T> PT = PowerTraceMatrix<T>(f);
    return UniversalScalarConversion<Tint, T>(DeterminantMat(PT));
  };
  std::set<int64_t> bad{2};
  for (Tint D : {disc_of(PhiE), disc_of(*PhiK0)}) {
    D = T_abs(D);
    for (int64_t p = 2;; p++) {
      Tint p_T = UniversalScalarConversion<Tint, int64_t>(p);
      if (p_T * p_T > D) {
        break;
      }
      while (D % p_T == 0) {
        bad.insert(p);
        D /= p_T;
      }
    }
    if (D > 1) {
      bad.insert(UniversalScalarConversion<int64_t, Tint>(D));
    }
  }
  auto primes_fn = [&](int64_t p) -> std::vector<QuadHeckePrime> {
    std::vector<HermIdealPrime<Tint>> lp = HermIdealPrimesAbove(ord, p);
    std::vector<QuadHeckePrime> res;
    for (size_t i = 0; i < lp.size(); i++) {
      if (HermK0Representative(lp, i) != i) {
        continue;
      }
      if (lp[i].type == HermIdealType::Split) {
        res.push_back({lp[i].f, 1, 0});
      } else if (lp[i].type == HermIdealType::Inert) {
        res.push_back({lp[i].f / 2, -1, 0});
      } else {
        res.push_back({lp[i].f, 0, 0});
      }
    }
    return res;
  };
  std::vector<int64_t> PhiK0v = *PhiK0;
  auto c1_fn = [PhiE, PhiK0v](int64_t q) -> int {
    return HermRootCountModQ(PhiE, q) - HermRootCountModQ(PhiK0v, q);
  };
  MyMatrix<Tint> Q = ord.trace_norm_form();
  MyMatrix<Tint> B0 = NullspaceIntMat(MyMatrix<Tint>(ord.conj - IdentityMat<Tint>(d)));
  T dE = T_abs(UniversalScalarConversion<T, Tint>(DeterminantMat(Q)));
  T dK = T_abs(DeterminantMat(MyMatrix<T>(
      UniversalMatrixConversion<T, Tint>(MyMatrix<Tint>(B0 * Q * B0.transpose())) / T(2))));
  T A2_T = dE / dK;
  if (!IsInteger(A2_T)) {
    std::cerr << "HERMMASS: |d_E| / |d_K0| is not an integer\n";
    throw TerminalException{1};
  }
  Tint A2 = UniversalScalarConversion<Tint, T>(A2_T);
  std::vector<QuadHeckeLData<T>> l_data;
  for (int i = 1; i <= m; i += 2) {
    l_data.push_back(MakeQuadHeckeLDataGeneric<T>(n, 1, A2, bad, primes_fn, c1_fn, i));
  }
  using mpfr::mpreal;
  // the size of b
  mpreal::set_default_prec(256);
  mpreal prod0 = mpreal_from_rational(factor, 256);
  for (auto &Ld : l_data) {
    std::optional<mpreal> L0 = QuadHeckeLNegative<mpreal>(Ld, 24, false);
    if (!L0) {
      std::cerr << "HERMMASS: an L-value cannot be computed\n";
      throw TerminalException{1};
    }
    prod0 *= *L0;
  }
  double log2_b = mpfr::log2(mpfr::abs(prod0) + 1).toDouble();
  int bits = std::max(48, static_cast<int>(std::ceil(log2_b)) + 32);
  for (int iter = 0; iter < 6; iter++, bits *= 2) {
#ifdef ENABLE_FLINT_SUPPORT
    ArbPrecisionGuard g(bits + 128);
    ArbReal prod = arb_from_rational(factor);
    bool ok = true;
    for (auto &Ld : l_data) {
      std::optional<ArbReal> La = QuadHeckeLNegative<ArbReal>(Ld, bits, false);
      if (!La) {
        ok = false;
        break;
      }
      prod = prod * *La;
    }
    if (!ok) {
      continue;
    }
    std::optional<Tint> b = arb_unique_integer<Tint>(prod);
    if (!b) {
      continue;
    }
#else
    auto rounded = [&](int prec) -> std::optional<Tint> {
      mpfr_prec_t wp = prec + 64;
      mpreal bf = mpreal_from_rational(factor, wp);
      for (auto &Ld : l_data) {
        std::optional<mpreal> L = QuadHeckeLNegative<mpreal>(Ld, prec, false);
        if (!L) {
          return {};
        }
        bf *= *L;
      }
      Tint bi = mpreal_nearest_integer<Tint>(bf);
      if (mpfr::abs(bf - mpreal_from_integer(bi, wp)) > mpfr::pow(mpreal(2, wp), -8)) {
        return {};
      }
      return bi;
    };
    std::optional<Tint> b = rounded(bits);
    std::optional<Tint> b2 = rounded(bits + 16);
    if (!b || !b2 || *b != *b2) {
      continue;
    }
#endif
    T mass = T_abs(UniversalScalarConversion<T, Tint>(*b)) /
             UniversalScalarConversion<T, Tint>(M);
#ifdef DEBUG_HERMITIAN_GENUS_MASS
    os << "HERMMASS: mass " << mass << " (local " << local << ", b = " << *b
       << ")\n";
#endif
    return mass;
  }
  std::cerr << "HERMMASS: the mass could not be made exact\n";
  throw TerminalException{1};
}

// clang-format off
#endif  // SRC_GENERA_HERMITIANGENUSMASS_H_
// clang-format on

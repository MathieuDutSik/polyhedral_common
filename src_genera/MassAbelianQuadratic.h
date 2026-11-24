// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_MASSABELIANQUADRATIC_H_
#define SRC_GENERA_MASSABELIANQUADRATIC_H_

// clang-format off
#include "MassAbelianQuadraticLocal.h"
#include "MassQuadraticHeckeL.h"
// clang-format on

#ifdef SANITY_CHECK
# define SANITY_CHECK_ABELIAN_QUADRATIC_MASS
#endif

/*
  The mass of a quadratic lattice over a totally real abelian field K
  (QuadraticMassAbelianLattice): the field and the discriminant of the
  lattice, the standard mass when E = K(sqrt disc) is abelian
  (MassAbelianQuadraticLocal.h) or from L(1 - r, chi_{E/K}) otherwise
  (MassQuadraticHeckeL.h), times the local factors.
 */

/*
  The field K = Q(theta) of a generator theta of the center, with a minimal
  polynomial f in Z[x], and the discriminant d = (-1)^{m(m-1)/2} det of the
  K-bilinear form on a K-basis of the lattice, on the powers of theta: E =
  K(sqrt d) is determined by its square class. Nothing when no small prime
  has a local generator.
 */
template <typename T>
inline std::optional<std::pair<NumberField<T>, ElemE<T>>>
QuadDiscriminantK(MyMatrix<T> const &B, std::vector<MyMatrix<T>> const &Omega) {
  int deg = static_cast<int>(Omega.size());
  int m = B.rows() / deg;
  auto is_prime = [](int64_t x) -> bool {
    for (int64_t p = 2; p * p <= x; p++) {
      if (x % p == 0) {
        return false;
      }
    }
    return x >= 2;
  };
  std::optional<LocalGenerator<T>> lg;
  int64_t q0 = 2;
  while (!lg) {
    if (is_prime(q0)) {
      lg = FindLocalGenerator(B, Omega, q0);
    }
    if (!lg) {
      q0++;
    }
    if (q0 > 1000) {
      return {};
    }
  }
  NumberField<T> K;
  K.d = 0;
  K.deg = deg;
  K.Phi = lg->f;
  K.conj = lg->conj;
  ElemE<T> d = DetK(LocalOKGram(B, *lg, q0), K);
  if ((m * (m - 1) / 2) % 2 == 1) {
    d = -d;
  }
  return std::make_pair(K, d);
}

template <typename T>
std::optional<AbelianField>
AbelianQuadraticExtensionField(MyMatrix<T> const &B,
                               std::vector<MyMatrix<T>> const &Omega,
                               AbelianField const &FK, ElemE<T> const &d,
                               NumberField<T> const &K, int m) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  T Nd = DeterminantMat(EMulMatrix(d, K));
  std::set<int64_t> ps;
  for (Tint x : {Tint(T_abs(GetNumerator_z(Nd))), Tint(GetDenominator_z(Nd))}) {
    if (x > 1) {
      for (Tint const &p : FactorsInt(x)) {
        ps.insert(UniversalScalarConversion<int64_t, Tint>(p));
      }
    }
  }
  return AbelianQuadraticExtensionFieldPrimes(B, Omega, FK, d, K, m, ps);
}

/*
  The mass of the genus of the quadratic lattice over K (the standard mass
  times the local factors), or nothing when a local factor is not computed
  or, for m even with a discriminant d not a square in K, when K(sqrt d) is
  not K(sqrt D) for a fundamental discriminant D (built from the primes
  dividing 2 n N(d)); E = K(sqrt D) is then the subfield of Q(zeta_N), N =
  LCMpair(n, |D|), fixed by the elements of H with Kronecker symbol (D / a) =
  1.
 */
template <typename T>
inline std::optional<T>
QuadraticMassAbelianLattice(MyMatrix<T> const &B,
                            std::vector<MyMatrix<T>> const &Omega,
                            AbelianField const &FK, QuadLValueMethod method) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int deg = static_cast<int>(Omega.size());
  int m = B.rows() / deg;
  if (m == 1) {
    return T(1) / T(2);
  }
  T local(1);
  for (auto const &kv : QuadraticLocalFactors(B, Omega, FK)) {
    if (!kv.second) {
      return {};
    }
    local *= *kv.second;
  }
  std::optional<T> sm;
  if (m % 2 == 1) {
    sm = StandardQuadraticMassAbelian<T>(FK, m, false, {});
  } else {
    std::optional<std::pair<NumberField<T>, ElemE<T>>> Kd =
        QuadDiscriminantK(B, Omega);
    if (!Kd) {
      return {};
    }
    NumberField<T> const &K = Kd->first;
    ElemE<T> const &d = Kd->second;
#ifdef DEBUG_ABQUAD
    std::cerr << "ABQUAD: disc =";
    for (auto const &x : d) {
      std::cerr << " " << x;
    }
    std::cerr << " (powers of theta, f =";
    for (auto const &x : K.Phi) {
      std::cerr << " " << x;
    }
    std::cerr << ")\n";
#endif
    if (SqrtInK(d, K)) {
      sm = StandardQuadraticMassAbelian<T>(FK, m, true, {});
    } else {
      std::set<int64_t> ps;
      T Nd = DeterminantMat(EMulMatrix(d, K));
      for (Tint x : {Tint(abs(GetNumerator_z(Nd))), Tint(GetDenominator_z(Nd)),
                          UniversalScalarConversion<Tint, int64_t>(2 * FK.n)}) {
        if (x > 1) {
          for (Tint const &p : FactorsInt(x)) {
            ps.insert(UniversalScalarConversion<int64_t, Tint>(p));
          }
        }
      }
      std::vector<int64_t> odd(ps.begin(), ps.end());
      odd.erase(std::remove(odd.begin(), odd.end(), 2), odd.end());
      for (uint64_t sub = 0; sub < (uint64_t(1) << odd.size()) && !sm; sub++) {
        int64_t Dodd = 1;
        for (size_t i = 0; i < odd.size(); i++) {
          if ((sub >> i) & 1) {
            Dodd *= (odd[i] % 4 == 1) ? odd[i] : -odd[i];
          }
        }
        for (int64_t two : {1, -4, 8, -8}) {
          int64_t Dq = Dodd * two;
          if (Dq == 1) {
            continue;
          }
          if (!SqrtInK(d * UniversalScalarConversion<T, int64_t>(Dq), K)) {
            continue;
          }
#ifdef DEBUG_ABQUAD
          std::cerr << "ABQUAD: E = K(sqrt " << Dq << ")\n";
#endif
          int64_t N = std::lcm(FK.n, std::abs(Dq));
          AbelianField FK2 = LiftAbelianField(FK, N);
          std::vector<int64_t> HE;
          for (int64_t a : FK2.H) {
            if (kronecker_symbol<int64_t>(Dq, a) == 1) {
              HE.push_back(a);
            }
          }
          sm = StandardQuadraticMassAbelian<T>(FK2, m, false,
                                            AbelianField{N, HE});
          break;
        }
      }
      if (!sm) {
        // E abelian but not K(sqrt D)
        std::optional<AbelianField> FE =
            AbelianQuadraticExtensionField(B, Omega, FK, d, K, m);
        if (!FE) {
          // E not abelian: L(1 - r, chi_{E/K}) numerically, rounded
          return QuadraticMassNonAbelianE(B, Omega, FK, d, K, m, local,
                                          method);
        }
        sm = StandardQuadraticMassAbelian<T>(LiftAbelianField(FK, FE->n), m,
                                             false, *FE);
      }
    }
  }
  if (!sm) {
    return {};
  }
  return *sm * local;
}

/*
  The image modulo p^c of the K-linear isometries of the trace form on X_p
  (the centralizer of a group component with center K, see
  MultiplicityLocalFactor), computed pi-adically. With X_p = (+)_P X_P over
  the primes P of K (or of the CM field) over p, all stable under the
  involution, the image is the product of the images of the U(X_P), each
  found in M_m(O_P / P^j) one pi-digit at a time:

  - h on m vectors x_i spanning X_p (LocalOKGramBasis), s_0 the least
    valuation of its entries; level j: g modulo pi^j with v(g H g^* - H) >=
    j + s_0;
  - level 1: the rows y_i of g modulo pi, the conditions y_i Hn y_k^T =
    Hn_{ik} and y_k Hn y_i^T = Hn_{ki} (k < i, Hn the residue of
    H / pi^{s_0}) being linear in y_i, the one for k = i quadratic;
  - g -> g + pi^j Delta, Delta modulo pi solving res(E / pi^{j+s_0}) +
    Delta A_1 + u_j A_2 Delta^T = 0, A_1 = res(H g^* / pi^{s_0}), A_2 =
    res(g H / pi^{s_0}), u_j the residue of conj(pi)^j / pi^j (the
    involution is trivial on the residue field);
  - the solutions modulo pi^{e c} are kept when they lift to pi^K, K =
    max(e c + e_H + t, e_H + 2 t + 1), e_H the exponent of the denominator
    of (H / pi^{s_0})^{-1} and t = v_P(2), by the Newton step g -> g -
    (1/2) E g^{-*} H^{-1}.

  Each is returned as the matrix on X modulo p^c (rows convention) of the
  map acting by g on X_P and as the identity on the other X_P'. Nothing
  when a prime is split over the real subfield or inert in a CM field,
  when there is no generator of O_K prime to the index, or beyond the
  bound on the solutions.
 */
template <typename T>
inline std::optional<std::vector<std::vector<std::vector<int64_t>>>>
LocalImagePiAdic(MyMatrix<T> const &B,
                 std::vector<MyMatrix<T>> const &Omega, int64_t p,
                 int c, size_t max_solution) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using IMat = std::vector<std::vector<int64_t>>;
  std::optional<LocalGenerator<T>> lg = FindLocalGenerator(B, Omega, p);
  if (!lg) {
    return {};
  }
  int deg = static_cast<int>(lg->f.size()) - 1;
  int r = B.rows();
  std::vector<std::pair<PolyQ, int>> fac =
      FactorModQWithMultiplicity(PolyReduceSigned(lg->f, p));
  {
    std::vector<AbelianPrime> aps = AbelianPrimesOverQ(*lg, p, 1);
    for (size_t i = 0; i < aps.size(); i++) {
      if (aps[i].conj_index != i) {
        return {};
      }
    }
  }
  LocalOKBasisGram<T> bg = LocalOKGramBasis(B, *lg, p);
  std::vector<std::vector<ElemE<T>>> const &H = bg.H;
  int m = static_cast<int>(H.size());
  // W: the rows x_i theta^t, index i deg + t
  std::vector<MyMatrix<T>> Spow{IdentityMat<T>(r)};
  for (int t = 1; t < deg; t++) {
    Spow.push_back(Spow.back() * lg->S);
  }
  MyMatrix<T> W(r, r);
  for (int i = 0; i < m; i++) {
    for (int t = 0; t < deg; t++) {
      MyVector<T> w = Spow[t].transpose() * bg.basis[i];
      for (int a = 0; a < r; a++) {
        W(i * deg + t, a) = w(a);
      }
    }
  }
  MyMatrix<T> Winv = Inverse(W);
  Tint pc = MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(p), c);
  auto mod_pc = [&](T const &x) -> int64_t {
    Tint den = GetDenominator_z(x);
    Tint inv = mod_inv(ResInt(den, pc), pc);
#ifdef SANITY_CHECK_ABELIAN_QUADRATIC_MASS
    if (inv == 0) {
      std::cerr << "ABQUAD: a denominator divisible by p\n";
      throw TerminalException{1};
    }
#endif
    Tint prod = GetNumerator_z(x) * inv;
    Tint v = ResInt(prod, pc);
    return UniversalScalarConversion<int64_t, Tint>(v);
  };

  std::vector<IMat> result;
  for (auto const &[g, e] : fac) {
    LocalPrimeK<T> P = MakeLocalPrimeK(*lg, g, e, p);
    NumberField<T> const &K = P.K;
    FiniteField F = MakeFiniteField(static_cast<T_small_prime>(p), g);
    int fr = P.f_res;
    // the conjugate of pi; the involution must be trivial on F_P
    ElemE<T> pibar = EConj(P.pi, K);
    {
      ElemE<T> th = EZetaPower(1, K);
      ElemE<T> d = EConj(th, K) - th;
      if (!d.is_zero() && LocalVal(d, P) < 1) {
        return {};
      }
    }
    int s0 = 1 << 28;
    for (auto const &row : H) {
      for (auto const &x : row) {
        if (!x.is_zero()) {
          s0 = std::min(s0, LocalVal(x, P));
        }
      }
    }
    auto res_shift = [&](ElemE<T> const &x, int k) -> FqElt {
      if (x.is_zero()) {
        return {};
      }
      return LocalResidue(LocalShift(x, k, P), P);
    };
    // the residue of H / pi^{s0}, and e_H
    std::vector<std::vector<FqElt>> Hn(m, std::vector<FqElt>(m));
    for (int a = 0; a < m; a++) {
      for (int b = 0; b < m; b++) {
        Hn[a][b] = res_shift(H[a][b], s0);
      }
    }
    int eH = 0;
    {
      // (H / pi^{s0})^{-1} = pi^{s0} H^{-1}
      std::vector<std::vector<ElemE<T>>> G = H;
      std::vector<std::vector<ElemE<T>>> Inv(
          m, std::vector<ElemE<T>>(m, ElemE<T>(deg, T(0))));
      for (int a = 0; a < m; a++) {
        Inv[a][a] = EZetaPower(0, K);
      }
      for (int col = 0; col < m; col++) {
        int piv = -1;
        for (int a = col; a < m; a++) {
          if (!G[a][col].is_zero()) {
            piv = a;
            break;
          }
        }
#ifdef SANITY_CHECK_ABELIAN_QUADRATIC_MASS
        if (piv < 0) {
          std::cerr << "ABQUAD: a degenerate local form\n";
          throw TerminalException{1};
        }
#endif
        std::swap(G[piv], G[col]);
        std::swap(Inv[piv], Inv[col]);
        ElemE<T> iv = EInv(G[col][col], K);
        for (int b = 0; b < m; b++) {
          G[col][b] = EMul(G[col][b], iv, K);
          Inv[col][b] = EMul(Inv[col][b], iv, K);
        }
        for (int a = 0; a < m; a++) {
          if (a == col || G[a][col].is_zero()) {
            continue;
          }
          ElemE<T> fct = G[a][col];
          for (int b = 0; b < m; b++) {
            G[a][b] -= EMul(fct, G[col][b], K);
            Inv[a][b] -= EMul(fct, Inv[col][b], K);
          }
        }
      }
      for (auto const &row : Inv) {
        for (auto const &x : row) {
          if (!x.is_zero()) {
            eH = std::max(eH, -(LocalVal(x, P) + s0));
          }
        }
      }
    }
    int tt = (p == 2) ? e : 0;
    int kT = e * c;
    int Kpi = std::max(kT + eH + tt, eH + 2 * tt + 1);
    // F_q-linear algebra on F_P = F_q^fr
    auto coeffs = [&](FqElt const &a) -> std::vector<int64_t> {
      std::vector<int64_t> v(fr, 0);
      for (size_t i = 0; i < a.size() && static_cast<int>(i) < fr; i++) {
        v[i] = a[i];
      }
      return v;
    };
    auto mulmat = [&](FqElt const &a) -> IMat {
      // column t: a x^t
      IMat M(fr, std::vector<int64_t>(fr, 0));
      for (int t = 0; t < fr; t++) {
        std::vector<int64_t> v = coeffs(FqMul(F, a, FqXPower(F, t)));
        for (int u = 0; u < fr; u++) {
          M[u][t] = v[u];
        }
      }
      return M;
    };
    // affine solutions of rows [coefficients | rhs] over F_q
    auto affine = [&](IMat const &aug, size_t nvar)
        -> std::optional<std::pair<std::vector<int64_t>, IMat>> {
      T_small_prime p_u = static_cast<T_small_prime>(p);
      int n_row = static_cast<int>(aug.size());
      int n_var = static_cast<int>(nvar);
      ModMatrix A(n_row, n_var + 1);
      for (int i = 0; i < n_row; i++) {
        for (int j = 0; j <= n_var; j++) {
          A(i, j) = ResModSigned(aug[i][j], p_u);
        }
      }
      // a pivot in the column of the right-hand side: no solution
      RowSpaceModP rs = ComputeRowSpaceModP(std::move(A), p_u);
      std::vector<int64_t> part(nvar, 0);
      for (size_t i = 0; i < rs.pivots.size(); i++) {
        if (rs.pivots[i] == nvar) {
          return {};
        }
        part[rs.pivots[i]] = ResToSigned(rs.basis(static_cast<int>(i), n_var));
      }
      ModMatrix K = FreeKernelBasisModP(rs, n_var, p_u);
      IMat ker(K.rows(), std::vector<int64_t>(nvar));
      for (int i = 0; i < K.rows(); i++) {
        for (int j = 0; j < n_var; j++) {
          ker[i][j] = ResToSigned(K(i, j));
        }
      }
      return std::make_pair(part, ker);
    };
    auto from_coeffs = [&](std::vector<int64_t> const &v, size_t off) -> FqElt {
      return FqFromCoeffs(F, std::vector<int64_t>(
                                 v.begin() + static_cast<long>(off),
                                 v.begin() + static_cast<long>(off + fr)));
    };
    using GMat = std::vector<std::vector<ElemE<T>>>;
    // level 1: the rows over F_P, all of them when there are not too many
    // (the list, when the chain below cannot be used, and in the sanity
    // checks)
    std::vector<GMat> level;
    bool level_over = false;
    auto list_level1 = [&]() -> void {
      std::vector<std::vector<FqElt>> ys;
      size_t n_node = 0;
      bool over = false;
      auto bil = [&](std::vector<FqElt> const &x, std::vector<FqElt> const &y) {
        FqElt v;
        for (int a = 0; a < m; a++) {
          for (int b = 0; b < m; b++) {
            v = FqAdd(F, v, FqMul(F, FqMul(F, x[a], Hn[a][b]), y[b]));
          }
        }
        return v;
      };
      std::function<void(int)> rec = [&](int i) {
        if (over) {
          return;
        }
        if (i == m) {
          GMat gm(m, std::vector<ElemE<T>>(m));
          for (int a = 0; a < m; a++) {
            for (int b = 0; b < m; b++) {
              gm[a][b] = LiftResidue(ys[a][b], K);
            }
          }
          level.push_back(gm);
          if (level.size() > max_solution) {
            over = true;
          }
          return;
        }
        size_t nvar = static_cast<size_t>(m * fr);
        IMat aug;
        for (int k = 0; k < i; k++) {
          // y_i Hn y_k^T = Hn_ik: coefficients w_b = (Hn y_k^T)_b
          for (int side = 0; side < 2; side++) {
            std::vector<FqElt> w(m);
            for (int b = 0; b < m; b++) {
              FqElt v;
              for (int d = 0; d < m; d++) {
                v = (side == 0) ? FqAdd(F, v, FqMul(F, Hn[b][d], ys[k][d]))
                                : FqAdd(F, v, FqMul(F, ys[k][d], Hn[d][b]));
              }
              w[b] = v;
            }
            FqElt rhs = (side == 0) ? Hn[i][k] : Hn[k][i];
            std::vector<IMat> blocks;
            for (int b = 0; b < m; b++) {
              blocks.push_back(mulmat(w[b]));
            }
            std::vector<int64_t> rc = coeffs(rhs);
            for (int u = 0; u < fr; u++) {
              std::vector<int64_t> row(nvar + 1, 0);
              for (int b = 0; b < m; b++) {
                for (int t = 0; t < fr; t++) {
                  row[b * fr + t] = blocks[b][u][t];
                }
              }
              row[nvar] = rc[u];
              aug.push_back(row);
            }
          }
        }
        std::optional<std::pair<std::vector<int64_t>, IMat>> sol;
        if (aug.empty()) {
          IMat ker;
          for (size_t v = 0; v < nvar; v++) {
            std::vector<int64_t> e1(nvar, 0);
            e1[v] = 1;
            ker.push_back(e1);
          }
          sol = std::make_pair(std::vector<int64_t>(nvar, 0), ker);
        } else {
          sol = affine(aug, nvar);
        }
        if (!sol) {
          return;
        }
        size_t n_y = 1;
        for (size_t u = 0; u < sol->second.size(); u++) {
          n_y *= static_cast<size_t>(p);
          if (n_y > 50000000) {
            over = true;
            return;
          }
        }
        for (size_t u = 0; u < n_y && !over; u++) {
          if (++n_node > 50000000) {
            over = true;
            return;
          }
          std::vector<int64_t> v = sol->first;
          size_t t2 = u;
          for (auto const &kv : sol->second) {
            int64_t cc = static_cast<int64_t>(t2 % static_cast<size_t>(p));
            t2 /= static_cast<size_t>(p);
            if (cc != 0) {
              for (size_t a = 0; a < nvar; a++) {
                v[a] = (v[a] + cc * kv[a]) % p;
              }
            }
          }
          std::vector<FqElt> y(m);
          for (int b = 0; b < m; b++) {
            y[b] = from_coeffs(v, static_cast<size_t>(b * fr));
          }
          if (bil(y, y) != Hn[i][i]) {
            continue;
          }
          // independence over F_P
          FqMat Mrows;
          for (auto const &yy : ys) {
            Mrows.push_back(yy);
          }
          Mrows.push_back(y);
          if (!FqLeftKernel(F, Mrows, Mrows.size(), static_cast<size_t>(m)).empty()) {
            continue;
          }
          ys.push_back(y);
          rec(i + 1);
          ys.pop_back();
        }
      };
      rec(0);
      if (over) {
        level_over = true;
        level.clear();
      }
    };
    // the lifts of g from pi^j to pi^{j+1}
    auto matmul = [&](GMat const &X, GMat const &Y) -> GMat {
      GMat Z(m, std::vector<ElemE<T>>(m, ElemE<T>(deg, T(0))));
      for (int a = 0; a < m; a++) {
        for (int b = 0; b < m; b++) {
          for (int d = 0; d < m; d++) {
            if (!X[a][d].is_zero() && !Y[d][b].is_zero()) {
              Z[a][b] += EMul(X[a][d], Y[d][b], K);
            }
          }
        }
      }
      return Z;
    };
    auto star = [&](GMat const &X) -> GMat {
      GMat Z(m, std::vector<ElemE<T>>(m));
      for (int a = 0; a < m; a++) {
        for (int b = 0; b < m; b++) {
          Z[a][b] = EConj(X[b][a], K);
        }
      }
      return Z;
    };
    // the residues of conj(pi)^j / pi^j
    std::vector<FqElt> ures;
    {
      ElemE<T> u1 = EMul(pibar, P.pi_inv, K);
      ElemE<T> uj = EZetaPower(0, K);
      for (int j = 0; j <= Kpi + 1; j++) {
        ures.push_back(LocalResidue(uj, P));
        uj = EMul(uj, u1, K);
      }
    }
    // the affine space of the Delta lifting g from pi^j to pi^{j+1}
    auto lift_space = [&](GMat const &gm, int j)
        -> std::optional<std::pair<std::vector<int64_t>, IMat>> {
      GMat gs = star(gm);
      GMat E = matmul(matmul(gm, H), gs);
      for (int a = 0; a < m; a++) {
        for (int b = 0; b < m; b++) {
          E[a][b] -= H[a][b];
#ifdef SANITY_CHECK_ABQUAD
          if (!E[a][b].is_zero() && LocalVal(E[a][b], P) < j + s0) {
            std::cerr << "ABQUAD: not a solution at level " << j << "\n";
            throw TerminalException{1};
          }
#endif
        }
      }
      GMat A1 = matmul(H, gs), A2 = matmul(gm, H);
      FqElt uj_res = ures[j];
      size_t nvar = static_cast<size_t>(m * m * fr);
      IMat aug;
      for (int a = 0; a < m; a++) {
        for (int b = 0; b < m; b++) {
          std::vector<IMat> coef(static_cast<size_t>(m * m));
          std::vector<bool> has(static_cast<size_t>(m * m), false);
          auto add = [&](int cidx, FqElt const &cf) {
            if (cf.is_zero()) {
              return;
            }
            IMat Mm = mulmat(cf);
            if (!has[cidx]) {
              coef[cidx] = Mm;
              has[cidx] = true;
            } else {
              for (int u = 0; u < fr; u++) {
                for (int t = 0; t < fr; t++) {
                  coef[cidx][u][t] = (coef[cidx][u][t] + Mm[u][t]) % p;
                }
              }
            }
          };
          for (int d = 0; d < m; d++) {
            // Delta_{a d} A1_{d b}
            add(a * m + d, res_shift(A1[d][b], s0));
            // u_j A2_{a d} Delta_{b d}
            add(b * m + d, FqMul(F, uj_res, res_shift(A2[a][d], s0)));
          }
          std::vector<int64_t> rc = coeffs(res_shift(E[a][b], j + s0));
          for (int u = 0; u < fr; u++) {
            std::vector<int64_t> row(nvar + 1, 0);
            for (int cidx = 0; cidx < m * m; cidx++) {
              if (!has[cidx]) {
                continue;
              }
              for (int t = 0; t < fr; t++) {
                row[static_cast<size_t>(cidx * fr + t)] = coef[cidx][u][t];
              }
            }
            row[nvar] = (p - rc[u]) % p;
            aug.push_back(row);
          }
        }
      }
      return affine(aug, nvar);
    };
    auto n_lifts = [&](std::pair<std::vector<int64_t>, IMat> const &sp) -> size_t {
      size_t n_y = MyPow<size_t>(static_cast<size_t>(p), sp.second.size());
      return n_y;
    };
    auto make_lift = [&](GMat const &gm, int j,
                         std::pair<std::vector<int64_t>, IMat> const &sp,
                         size_t u) -> GMat {
      size_t nvar = sp.first.size();
      std::vector<int64_t> v = sp.first;
      size_t t2 = u;
      for (auto const &kv : sp.second) {
        int64_t cc = static_cast<int64_t>(t2 % static_cast<size_t>(p));
        t2 /= static_cast<size_t>(p);
        if (cc != 0) {
          for (size_t a = 0; a < nvar; a++) {
            v[a] = (v[a] + cc * kv[a]) % p;
          }
        }
      }
      ElemE<T> pij = PiPower(j, P);
      GMat g2 = gm;
      for (int a = 0; a < m; a++) {
        for (int b = 0; b < m; b++) {
          FqElt dlt = from_coeffs(v, static_cast<size_t>((a * m + b) * fr));
          if (!dlt.is_zero()) {
            g2[a][b] += EMul(LiftResidue(dlt, K), pij, K);
          }
        }
      }
      return g2;
    };
    auto list_levels = [&]() -> void {
    list_level1();
    for (int j = 1; j < kT && !level_over; j++) {
      std::vector<GMat> next;
      for (auto const &gm : level) {
        std::optional<std::pair<std::vector<int64_t>, IMat>> sp = lift_space(gm, j);
        if (!sp) {
          continue;
        }
        if (sp->second.size() > 40) {
          level_over = true;
          break;
        }
        size_t n_y = n_lifts(*sp);
        if (next.size() + n_y > max_solution) {
          level_over = true;
          break;
        }
        for (size_t u = 0; u < n_y; u++) {
          next.push_back(make_lift(gm, j, *sp, u));
        }
      }
      level = std::move(next);
    }
    if (level_over) {
      level.clear();
    }
    };
    // one lift at a time, the first that lifts further ending the search
    std::function<bool(GMat const &, int)> liftable = [&](GMat const &gm,
                                                         int j) -> bool {
      if (j >= Kpi) {
        return true;
      }
      std::optional<std::pair<std::vector<int64_t>, IMat>> sp = lift_space(gm, j);
      if (!sp) {
        return false;
      }
      size_t n_y = n_lifts(*sp);
      for (size_t u = 0; u < n_y; u++) {
        if (liftable(make_lift(gm, j, *sp, u), j + 1)) {
          return true;
        }
      }
      return false;
    };
    // the matrices on X modulo p^c
    int Nprec = c + 1;
    Tint Qp = MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(p), Nprec);
    PolyZ<Tint> eps = LocalEps(P, Nprec);
    auto to_X = [&](GMat const &gm) -> IMat {
      MyMatrix<T> Mf = ZeroMatrix<T>(r, r);
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
          auto [yz, dz] = ClearDenominators(gm[i][j]);
          Tint mod = Qp * dz;
          PolyZ<Tint> t = PolyZReduce(PolyZMul(eps, yz), lg->f, mod);
          ElemE<T> tot(deg, T(0));
          for (int u = 0; u < deg; u++) {
            tot[u] = T(t[u]) / T(dz);
          }
          if (i == j) {
            // + (1 - eps)
            for (int u = 0; u < deg; u++) {
              tot[u] -= T(eps[u]);
            }
            tot[0] += 1;
          }
          MyMatrix<T> Mm = EMulMatrix(tot, K);
          for (int t1 = 0; t1 < deg; t1++) {
            for (int t2 = 0; t2 < deg; t2++) {
              Mf(i * deg + t1, j * deg + t2) = Mm(t1, t2);
            }
          }
        }
      }
      MyMatrix<T> Me = Winv * Mf * W;
      IMat out(r, std::vector<int64_t>(r));
      for (int a = 0; a < r; a++) {
        for (int b = 0; b < r; b++) {
          out[a][b] = mod_pc(Me(a, b));
        }
      }
      return out;
    };
    std::vector<IMat> listed;
    auto list_image = [&]() -> void {
      list_levels();
      for (auto const &gm : level) {
        if (liftable(gm, kT)) {
          listed.push_back(to_X(gm));
        }
      }
    };
    /*
      When the list overflows: generators of the image of U(X_P), complete
      by construction, as LocalOrthogonalImageGenerators does over Z_p but
      on the rows over F_P (the residue field of P) and the pi-adic levels:
      the stabilizer chain modulo pi of the basis rows e_0, ..., e_{m-1}
      (one element searched for each new image of e_i, whose rows are
      completed one at a time on their affine spaces and lifted to
      pi^{Kpi} by the lifts of lift_space, which by the Lemma comes from an
      isometry), then the layers of the elements congruent to 1 modulo pi^j,
      j < kT, F_q-spaces given by a basis of the solutions of lift_space at
      1 when all of its elements lift, else by enumerating them (bounded).
     */
    auto run_chain = [&](std::vector<IMat> &out_gens) -> bool {
      size_t nd = static_cast<size_t>(m * fr);
      size_t n_pts = 1;
      for (size_t u = 0; u < nd; u++) {
        n_pts *= static_cast<size_t>(p);
        if (n_pts > (size_t(1) << 24)) {
          return false;
        }
      }
      using Row = std::vector<FqElt>;
      auto encode_row = [&](Row const &y) -> size_t {
        size_t idx = 0;
        for (int b = m - 1; b >= 0; b--) {
          std::vector<int64_t> cf = coeffs(y[b]);
          for (int u = fr - 1; u >= 0; u--) {
            idx = idx * static_cast<size_t>(p) + static_cast<size_t>(((cf[u] % p) + p) % p);
          }
        }
        return idx;
      };
      auto decode_row = [&](size_t idx) -> Row {
        Row y(m);
        for (int b = 0; b < m; b++) {
          std::vector<int64_t> cf(fr);
          for (int u = 0; u < fr; u++) {
            cf[u] = static_cast<int64_t>(idx % static_cast<size_t>(p));
            idx /= static_cast<size_t>(p);
          }
          y[b] = FqFromCoeffs(F, cf);
        }
        return y;
      };
      auto bil2 = [&](Row const &x, Row const &y) -> FqElt {
        FqElt v;
        for (int a2 = 0; a2 < m; a2++) {
          for (int b2 = 0; b2 < m; b2++) {
            v = FqAdd(F, v, FqMul(F, FqMul(F, x[a2], Hn[a2][b2]), y[b2]));
          }
        }
        return v;
      };
      // the affine space of the rows y with y Hn x_k^T = Hn_ik and x_k Hn
      // y^T = Hn_ki for the rows x_k of the prefix
      auto row_aff = [&](std::vector<Row> const &prefix, int i)
          -> std::optional<std::pair<std::vector<int64_t>, IMat>> {
        size_t nvar = nd;
        IMat aug;
        for (size_t k = 0; k < prefix.size(); k++) {
          for (int side = 0; side < 2; side++) {
            Row w(m);
            for (int b2 = 0; b2 < m; b2++) {
              FqElt v;
              for (int d2 = 0; d2 < m; d2++) {
                v = (side == 0) ? FqAdd(F, v, FqMul(F, Hn[b2][d2], prefix[k][d2]))
                                : FqAdd(F, v, FqMul(F, prefix[k][d2], Hn[d2][b2]));
              }
              w[b2] = v;
            }
            FqElt rhs = (side == 0) ? Hn[i][k] : Hn[k][i];
            std::vector<IMat> blocks;
            for (int b2 = 0; b2 < m; b2++) {
              blocks.push_back(mulmat(w[b2]));
            }
            std::vector<int64_t> rc = coeffs(rhs);
            for (int u = 0; u < fr; u++) {
              std::vector<int64_t> row(nvar + 1, 0);
              for (int b2 = 0; b2 < m; b2++) {
                for (int t3 = 0; t3 < fr; t3++) {
                  row[b2 * fr + t3] = blocks[b2][u][t3];
                }
              }
              row[nvar] = rc[u];
              aug.push_back(row);
            }
          }
        }
        if (aug.empty()) {
          IMat ker;
          for (size_t v = 0; v < nvar; v++) {
            std::vector<int64_t> e1(nvar, 0);
            e1[v] = 1;
            ker.push_back(e1);
          }
          return std::make_pair(std::vector<int64_t>(nvar, 0), ker);
        }
        return affine(aug, nvar);
      };
      size_t n_node2 = 0;
      bool over2 = false;
      const size_t max_node2 = 5000000;
      // f(y) for the rows y of the affine space with y Hn y^T = Hn_ii and
      // independent of the prefix, until f is true
      auto for_each_row = [&](std::vector<Row> const &prefix, int i,
                              std::function<bool(Row const &)> const &f) -> void {
        std::optional<std::pair<std::vector<int64_t>, IMat>> sol =
            row_aff(prefix, i);
        if (!sol) {
          return;
        }
        size_t dim = sol->second.size();
        std::vector<int64_t> cc(dim, 0);
        while (!over2) {
          if (++n_node2 > max_node2) {
            over2 = true;
            return;
          }
          std::vector<int64_t> v = sol->first;
          for (size_t u = 0; u < dim; u++) {
            if (cc[u] != 0) {
              for (size_t a2 = 0; a2 < nd; a2++) {
                v[a2] = (v[a2] + cc[u] * sol->second[u][a2]) % p;
              }
            }
          }
          Row y(m);
          for (int b2 = 0; b2 < m; b2++) {
            y[b2] = from_coeffs(v, static_cast<size_t>(b2 * fr));
          }
          if (bil2(y, y) == Hn[i][i]) {
            FqMat Mrows(prefix.begin(), prefix.end());
            Mrows.push_back(y);
            if (FqLeftKernel(F, Mrows, Mrows.size(), static_cast<size_t>(m))
                    .empty() &&
                f(y)) {
              return;
            }
          }
          size_t u = 0;
          while (u < dim && cc[u] == p - 1) {
            cc[u] = 0;
            u++;
          }
          if (u == dim) {
            return;
          }
          cc[u]++;
        }
      };
      // g + pi^j (sum of the coefficient vector times the kernel), lazily
      auto make_lift_vec = [&](GMat const &gm, int j,
                               std::pair<std::vector<int64_t>, IMat> const &sp,
                               std::vector<int64_t> const &cvec) -> GMat {
        size_t nvar = sp.first.size();
        std::vector<int64_t> v = sp.first;
        for (size_t u = 0; u < sp.second.size(); u++) {
          if (cvec[u] != 0) {
            for (size_t a2 = 0; a2 < nvar; a2++) {
              v[a2] = (v[a2] + cvec[u] * sp.second[u][a2]) % p;
            }
          }
        }
        ElemE<T> pij = PiPower(j, P);
        GMat g2 = gm;
        for (int a2 = 0; a2 < m; a2++) {
          for (int b2 = 0; b2 < m; b2++) {
            FqElt dlt = from_coeffs(v, static_cast<size_t>((a2 * m + b2) * fr));
            if (!dlt.is_zero()) {
              g2[a2][b2] += EMul(LiftResidue(dlt, K), pij, K);
            }
          }
        }
        return g2;
      };
      // a solution modulo pi^{Kpi} extending the solution gm modulo pi^j,
      // the first found; by the Lemma it is congruent modulo pi^{kT} to an
      // element of the image
      std::function<std::optional<GMat>(GMat const &, int)> path =
          [&](GMat const &gm, int j) -> std::optional<GMat> {
        if (j >= Kpi) {
          return gm;
        }
        if (over2 || ++n_node2 > max_node2) {
          over2 = true;
          return {};
        }
        std::optional<std::pair<std::vector<int64_t>, IMat>> sp =
            lift_space(gm, j);
        if (!sp) {
          return {};
        }
        size_t dim = sp->second.size();
        std::vector<int64_t> cv(dim, 0);
        while (!over2) {
          std::optional<GMat> res = path(make_lift_vec(gm, j, *sp, cv), j + 1);
          if (res) {
            return res;
          }
          size_t u = 0;
          while (u < dim && cv[u] == p - 1) {
            cv[u] = 0;
            u++;
          }
          if (u == dim) {
            return {};
          }
          cv[u]++;
        }
        return {};
      };
      auto act_row = [&](Row const &y, std::vector<Row> const &g) -> Row {
        Row z(m);
        for (int b2 = 0; b2 < m; b2++) {
          FqElt v;
          for (int a2 = 0; a2 < m; a2++) {
            v = FqAdd(F, v, FqMul(F, y[a2], g[a2][b2]));
          }
          z[b2] = v;
        }
        return z;
      };
      std::vector<std::vector<Row>> gens1;
      std::vector<GMat> elems;
      [[maybe_unused]] double log_order = 0;
      std::vector<Row> rows(m);
      std::function<std::optional<GMat>(int)> complete =
          [&](int j) -> std::optional<GMat> {
        if (j == m) {
          GMat gm(m, std::vector<ElemE<T>>(m));
          for (int a2 = 0; a2 < m; a2++) {
            for (int b2 = 0; b2 < m; b2++) {
              gm[a2][b2] = LiftResidue(rows[a2][b2], K);
            }
          }
          return path(gm, 1);
        }
        std::vector<Row> prefix(rows.begin(), rows.begin() + j);
        std::optional<GMat> res;
        for_each_row(prefix, j, [&](Row const &y) -> bool {
          rows[j] = y;
          res = complete(j + 1);
          return res.has_value();
        });
        return res;
      };
      for (int i = m - 1; i >= 0; i--) {
        std::vector<Row> prefix(i);
        for (int k = 0; k < i; k++) {
          Row ek(m);
          ek[k] = FqOne(F);
          prefix[k] = ek;
        }
        Row ei(m);
        ei[i] = FqOne(F);
        std::vector<char> in_orbit(n_pts, 0);
        size_t orbit_len = 0;
        auto close = [&]() -> void {
          std::vector<size_t> queue;
          if (orbit_len == 0) {
            in_orbit[encode_row(ei)] = 1;
            orbit_len = 1;
          }
          for (size_t idx = 0; idx < n_pts; idx++) {
            if (in_orbit[idx]) {
              queue.push_back(idx);
            }
          }
          while (!queue.empty()) {
            size_t idx = queue.back();
            queue.pop_back();
            Row y = decode_row(idx);
            for (auto const &g : gens1) {
              size_t iz = encode_row(act_row(y, g));
              if (!in_orbit[iz]) {
                in_orbit[iz] = 1;
                orbit_len++;
                queue.push_back(iz);
              }
            }
          }
        };
        close();
        [[maybe_unused]] size_t n_cand = 0, n_fail = 0;
        for_each_row(prefix, i, [&](Row const &y) -> bool {
          n_cand++;
          if (in_orbit[encode_row(y)]) {
            return false;
          }
          for (int k = 0; k < i; k++) {
            rows[k] = prefix[k];
          }
          rows[i] = y;
          std::optional<GMat> g = complete(i + 1);
          if (over2) {
            return true;
          }
          if (g) {
            std::vector<Row> g1(m, Row(m));
            for (int a2 = 0; a2 < m; a2++) {
              for (int b2 = 0; b2 < m; b2++) {
                g1[a2][b2] = LocalResidue((*g)[a2][b2], P);
              }
            }
            gens1.push_back(g1);
            elems.push_back(*g);
            close();
          } else {
            n_fail++;
          }
          return false;
        });
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
        std::cerr << "ABQUAD: pi-adic chain at row " << i << ": candidates "
                  << n_cand << " failed " << n_fail << " orbit " << orbit_len
                  << " nodes " << n_node2 << "\n";
#endif
        if (over2) {
          return false;
        }
        log_order += std::log(static_cast<double>(orbit_len));
      }
      // the layers
      GMat Id(m, std::vector<ElemE<T>>(m, ElemE<T>(deg, T(0))));
      for (int a2 = 0; a2 < m; a2++) {
        Id[a2][a2] = EZetaPower(0, K);
      }
      for (int j = 1; j < kT; j++) {
        std::optional<std::pair<std::vector<int64_t>, IMat>> sp =
            lift_space(Id, j);
        if (!sp) {
          return false;
        }
        size_t dim = sp->second.size();
        std::vector<GMat> lifts;
        bool all = true;
        for (size_t u = 0; u < dim; u++) {
          std::vector<int64_t> cv(dim, 0);
          cv[u] = 1;
          std::optional<GMat> g = path(make_lift_vec(Id, j, *sp, cv), j + 1);
          if (over2) {
            return false;
          }
          if (!g) {
            all = false;
            break;
          }
          lifts.push_back(*g);
        }
        if (all) {
          elems.insert(elems.end(), lifts.begin(), lifts.end());
          log_order += static_cast<double>(dim) * std::log(static_cast<double>(p));
        } else {
          // the liftable subspace, by enumerating the space
          size_t n_el = 1;
          for (size_t u = 0; u < dim; u++) {
            n_el *= static_cast<size_t>(p);
            if (n_el > 100000) {
              return false;
            }
          }
          IMat basis_cv;
          size_t n_ok = 0;
          for (size_t s2 = 0; s2 < n_el; s2++) {
            std::vector<int64_t> cv(dim);
            size_t tt = s2;
            for (size_t u = 0; u < dim; u++) {
              cv[u] = static_cast<int64_t>(tt % static_cast<size_t>(p));
              tt /= static_cast<size_t>(p);
            }
            std::optional<GMat> g = path(make_lift_vec(Id, j, *sp, cv), j + 1);
            if (over2) {
              return false;
            }
            if (g) {
              IMat test = basis_cv;
              test.push_back(cv);
              size_t rk = 0;
              {
                IMat aug2 = test;
                for (auto &row : aug2) {
                  row.push_back(0);
                }
                std::optional<std::pair<std::vector<int64_t>, IMat>> kk =
                    affine(aug2, dim);
                rk = (kk) ? dim - kk->second.size() : 0;
              }
              if (rk == test.size()) {
                basis_cv = test;
                elems.push_back(*g);
                n_ok++;
              }
            }
          }
          log_order += static_cast<double>(n_ok) * std::log(static_cast<double>(p));
        }
      }
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
      std::cerr << "ABQUAD: pi-adic chain: " << elems.size()
                << " generators, |image of U(X_P)| = " << std::exp(log_order)
                << "\n";
#endif
      for (auto const &gm : elems) {
        out_gens.push_back(to_X(gm));
      }
      return true;
    };
    // the chain first; the list when the chain cannot be used
    std::vector<IMat> chain_gens;
    bool chain_done = run_chain(chain_gens);
#ifndef SANITY_CHECK_MULT
    if (chain_done) {
      result.insert(result.end(), chain_gens.begin(), chain_gens.end());
      continue;
    }
#endif
    list_image();
    if (level_over) {
      if (!chain_done) {
        return {};
      }
      result.insert(result.end(), chain_gens.begin(), chain_gens.end());
      continue;
    }
#ifdef SANITY_CHECK_MULT
    {
      // the group generated by the chain is the listed image
      if (chain_done) {
        std::set<IMat> listed_set(listed.begin(), listed.end());
        std::set<IMat> closure;
        IMat Idm(r, std::vector<int64_t>(r, 0));
        for (int a2 = 0; a2 < r; a2++) {
          Idm[a2][a2] = 1;
        }
        int64_t pcl = UniversalScalarConversion<int64_t, Tint>(pc);
        std::vector<IMat> queue{Idm};
        closure.insert(Idm);
        while (!queue.empty() && closure.size() <= listed_set.size()) {
          IMat x = queue.back();
          queue.pop_back();
          for (auto const &g : chain_gens) {
            IMat y(r, std::vector<int64_t>(r, 0));
            for (int a2 = 0; a2 < r; a2++) {
              for (int k2 = 0; k2 < r; k2++) {
                if (x[a2][k2] != 0) {
                  for (int b2 = 0; b2 < r; b2++) {
                    y[a2][b2] = (y[a2][b2] + x[a2][k2] * g[k2][b2]) % pcl;
                  }
                }
              }
            }
            if (closure.insert(y).second) {
              queue.push_back(y);
            }
          }
        }
        if (closure != listed_set) {
          std::cerr << "ABQUAD: SANITY: the pi-adic chain generates "
                    << closure.size() << " elements, the list has "
                    << listed_set.size() << "\n";
          throw TerminalException{1};
        }
        std::cerr << "ABQUAD: SANITY: the pi-adic chain generates the "
                  << listed_set.size() << " elements of the list\n";
      }
    }
#endif
    (void)chain_done;
    result.insert(result.end(), listed.begin(), listed.end());
  }
  return result;
}

// clang-format off
#endif  // SRC_GENERA_MASSABELIANQUADRATIC_H_
// clang-format on

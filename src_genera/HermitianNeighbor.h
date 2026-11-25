// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_HERMITIANNEIGHBOR_H_
#define SRC_GENERA_HERMITIANNEIGHBOR_H_

// clang-format off
#include "CMOrder.h"
#include "ResidueModInteger.h"
#include "ScalarFiniteField.h"
#include "ScalarPolynomialModQ.h"
#include "FiniteFieldProjective.h"
#include "Basic_random.h"
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_HERMITIAN_NEIGHBOR
#endif

/*
  The Kneser neighbors of a Hermitian lattice over the order O_E of a CM
  field (CMOrder.h) at any prime P of O_E, as neighbours of Hecke.jl
  (QuadForm/Herm/GenusRep.jl).

  The primes. Over a rational prime ell, A = O_E / ell O_E is a finite
  commutative F_ell-algebra. Its radical R is the kernel of x -> x^(ell^j),
  ell^j >= [E : Q], and A / R is the product of the residue fields of the
  primes above ell. An element z (from a generator of fixed seed, the same
  on every run) whose minimal polynomial has
  distinct irreducible factors g_i of total degree dim A / R generates every
  residue field: the primary idempotents of F_ell[z] are then the primitive
  idempotents of A, the prime P_i is the kernel of O_E -> F_ell[x] / g_i, z
  -> x, its residue degree is deg g_i and its ramification index is dim
  eps_i A / deg g_i. This needs neither a generator of O_E at ell nor the
  index of Z[theta]. The conjugate of P is the prime containing the
  conjugates of its elements, and P is, over the real subfield K_0 with
  p = P cap K_0: split when conj(P) != P, ramified when conj(P) = P and the
  conjugation is trivial on O_E / P, inert otherwise.

  The neighbors, for L unimodular at the primes above ell (h(L, L) in O_E
  and L equal to its dual there: v_ell(det Tr h) = m v_ell(disc O_E)), at a
  prime P over p, C = conj(P) when p splits and C = P otherwise, for a
  vector x of L:

    M = { y in L : h(x, y) in C },   L' = M + C^{-1} x,

  which is an O_E-lattice, integral and unimodular at p, of the covolume
  of L:
  - p split: x = pi w, pi in P but not in P^2 nor in C, for any line w of
    L / C L; every line gives one neighbor;
  - p inert: w isotropic, h(w, w) in p, and x = w + pi s e_k, pi a
    uniformizer of p in K_0 and e_k a basis vector with h(e_k, w) not in
    P, for the q_p elements s of O_E / P with h(x, x) in p^2 (an affine
    condition on s: Tr(s h(e_k, w)) = -h(w, w) / pi modulo P);
  - p ramified: w isotropic, x = w + pi s e_k, pi a uniformizer of P, for
    the q elements s of O_E / P.
  The odd scales at a ramified prime (special in Hecke.jl) and the dyadic
  ramified primes are not treated: such a prime is not used.
 */

// The algebra O_E / ell O_E on the basis e of the order.
struct HermResidueAlgebra {
  int64_t ell;
  int d;
  // mult[a][b * d + c]: the multiplication by e_a, row b column c, mod ell
  std::vector<std::vector<int64_t>> mult;
  std::vector<int64_t> product(std::vector<int64_t> const &x,
                               std::vector<int64_t> const &y) const {
    std::vector<int64_t> res(d, 0);
    for (int a = 0; a < d; a++) {
      if (x[a] == 0) {
        continue;
      }
      for (int b = 0; b < d; b++) {
        if (y[b] == 0) {
          continue;
        }
        int64_t coef = x[a] * y[b] % ell;
        for (int c = 0; c < d; c++) {
          res[c] = (res[c] + coef * mult[a][b * d + c]) % ell;
        }
      }
    }
    return res;
  }
  std::vector<int64_t> one() const {
    std::vector<int64_t> v(d, 0);
    v[0] = 1;
    return v;
  }
  std::vector<int64_t> power(std::vector<int64_t> x, uint64_t e) const {
    std::vector<int64_t> res = one();
    while (e > 0) {
      if (e & 1) {
        res = product(res, x);
      }
      x = product(x, x);
      e >>= 1;
    }
    return res;
  }
  std::vector<int64_t> poly_at(PolyQ const &g,
                               std::vector<int64_t> const &z) const {
    std::vector<int64_t> res(d, 0);
    for (int k = g.degree(); k >= 0; k--) {
      res = product(res, z);
      res[0] = (res[0] + static_cast<int64_t>(g.coeff(k))) % ell;
    }
    return res;
  }
};

template <typename Tint>
MyMatrix<Tint> HermRowsToMatrix(std::vector<std::vector<int64_t>> const &rows,
                                int n_col) {
  MyMatrix<Tint> M(rows.size(), n_col);
  for (size_t i = 0; i < rows.size(); i++) {
    for (int j = 0; j < n_col; j++) {
      M(i, j) = UniversalScalarConversion<Tint, int64_t>(rows[i][j]);
    }
  }
  return M;
}

/*
  A prime P of O_E above ell: the residue field F = F_ell[x] / g with z ->
  x, the residues of the basis e, the Z-basis of P and of P^2 (rows, on e),
  e and f over Q, the index of conj(P) and the type over K_0.
 */
enum class HermIdealType { Split, Inert, Ramified };

template <typename Tint> struct HermIdealPrime {
  int64_t ell;
  int e;
  int f;
  FiniteField F;
  // the residue of e_a, and the lift of the generator x of F in O_E
  std::vector<FqElt> res_basis;
  MyVector<Tint> z;
  MyMatrix<Tint> ideal;
  MyMatrix<Tint> ideal_sq;
  size_t conj_index;
  HermIdealType type;
  std::string name;
  // P^1, P^2, ... as Z-bases, computed on demand by HermIdealVal
  mutable std::vector<MyMatrix<Tint>> powers;
  uint64_t q() const { return F.q; }
};

template <typename Tint>
FqElt HermResidue(HermIdealPrime<Tint> const &P, MyVector<Tint> const &a) {
  FqElt res = FqFromCoeffs(P.F, std::vector<int64_t>{});
  for (int k = 0; k < a.size(); k++) {
    int64_t c = order_residue(a(k), P.ell);
    if (c != 0) {
      FqElt term = FqMul(P.F, FqFromCoeffs(P.F, std::vector<int64_t>{c}),
                         P.res_basis[k]);
      res = FqAdd(P.F, res, term);
    }
  }
  return res;
}

// The residue of an element of E with denominators prime to ell.
template <typename T, typename Tint>
FqElt HermResidueRat(HermIdealPrime<Tint> const &P, MyVector<T> const &a) {
  FqElt res = FqFromCoeffs(P.F, std::vector<int64_t>{});
  for (int k = 0; k < a.size(); k++) {
    int64_t c = order_residue_rat(a(k), P.ell);
    if (c != 0) {
      FqElt term = FqMul(P.F, FqFromCoeffs(P.F, std::vector<int64_t>{c}),
                         P.res_basis[k]);
      res = FqAdd(P.F, res, term);
    }
  }
  return res;
}

// Whether the integral vector a (on e) is in the Z-lattice of the rows of B.
template <typename T, typename Tint>
bool HermInLattice(MyMatrix<Tint> const &B, MyVector<Tint> const &a) {
  std::optional<MyVector<T>> opt =
      SolutionMat(UniversalMatrixConversion<T, Tint>(B),
                  UniversalVectorConversion<T, Tint>(a));
  return opt && IsIntegralVector(*opt);
}

/*
  The element of P with the valuation 1 at P (and, with not_in, outside the
  other primes listed), among the elements of the Z-lattice of the rows of
  B (the ideal P, or P cap K_0), by random small combinations.
 */
template <typename T, typename Tint>
MyVector<Tint>
HermFindUniformizer(MyMatrix<Tint> const &B, HermIdealPrime<Tint> const &P,
                    std::vector<HermIdealPrime<Tint>> const &l_other) {
  int k = B.rows();
  int d = B.cols();
  for (int iter = 0; iter < 10000; iter++) {
    int range = 1 + iter / 1000;
    MyVector<Tint> a = ZeroVector<Tint>(d);
    for (int i = 0; i < k; i++) {
      int c = random_int(-range, range);
      if (c != 0) {
        a += Tint(c) * GetMatrixRow(B, i);
      }
    }
    if (HermInLattice<T, Tint>(P.ideal_sq, a)) {
      continue;
    }
    bool ok = true;
    for (auto const &Q : l_other) {
      // a must not be in the other primes above ell
      if (HermResidue(Q, a).is_zero()) {
        ok = false;
        break;
      }
    }
    if (ok) {
      return a;
    }
  }
  std::cerr << "HERMNB: no uniformizer found at " << P.name << "\n";
  throw TerminalException{1};
}

/*
  The valuation v_P of an element of E (rational coordinates on e): the
  powers of ell are removed (v_P(ell) = e), then the membership in the
  successive powers of P, kept in P.powers.
 */
template <typename T, typename Tint>
int HermIdealVal(CMOrder<T, Tint> const &ord, HermIdealPrime<Tint> const &P,
                 MyVector<T> const &a) {
  int d = ord.deg;
  if (IsZeroVector(a)) {
    return 1 << 28;
  }
  Tint ell_T = UniversalScalarConversion<Tint, int64_t>(P.ell);
  Tint den(1);
  for (int k = 0; k < d; k++) {
    den = LCMpair(den, GetDenominator_z(a(k)));
  }
  // a = num / den: v_P(a) = v_P(num) - e v_ell(den)
  MyVector<Tint> num(d);
  for (int k = 0; k < d; k++) {
    num(k) = UniversalScalarConversion<Tint, T>(a(k) * UniversalScalarConversion<T, Tint>(den));
  }
  int val = 0;
  while (den % ell_T == 0) {
    den /= ell_T;
    val -= P.e;
  }
  while (true) {
    bool all = true;
    for (int k = 0; k < d && all; k++) {
      all = (num(k) % ell_T == 0);
    }
    if (!all) {
      break;
    }
    num /= ell_T;
    val += P.e;
  }
  if (P.powers.empty()) {
    P.powers.push_back(P.ideal);
  }
  int k = 0;
  while (true) {
    if (static_cast<int>(P.powers.size()) <= k) {
      MyMatrix<Tint> const &A = P.powers.back();
      MyMatrix<Tint> Gen(A.rows() * P.ideal.rows(), d);
      for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < P.ideal.rows(); j++) {
          MyVector<Tint> pr = ord.product(GetMatrixRow(A, i), GetMatrixRow(P.ideal, j));
          for (int c = 0; c < d; c++) {
            Gen(i * P.ideal.rows() + j, c) = pr(c);
          }
        }
      }
      P.powers.push_back(GetZbasis(Gen));
    }
    if (!HermInLattice<T, Tint>(P.powers[k], num)) {
      return val + k;
    }
    k++;
  }
}

/*
  The primes of O_E above ell.
 */
template <typename T, typename Tint>
std::vector<HermIdealPrime<Tint>> HermIdealPrimesAbove(CMOrder<T, Tint> const &ord,
                                             int64_t ell) {
  int d = ord.deg;
  HermResidueAlgebra A;
  A.ell = ell;
  A.d = d;
  for (int a = 0; a < d; a++) {
    std::vector<int64_t> M(d * d);
    for (int b = 0; b < d; b++) {
      for (int c = 0; c < d; c++) {
        M[b * d + c] = order_residue(ord.mult[a](b, c), ell);
      }
    }
    A.mult.push_back(M);
  }
  Tint ell_T = UniversalScalarConversion<Tint, int64_t>(ell);
  // the radical: the kernel of the Frobenius power
  uint64_t qj = static_cast<uint64_t>(ell);
  while (qj < static_cast<uint64_t>(d)) {
    qj *= static_cast<uint64_t>(ell);
  }
  std::vector<std::vector<int64_t>> frob;
  for (int a = 0; a < d; a++) {
    std::vector<int64_t> ea(d, 0);
    ea[a] = 1;
    frob.push_back(A.power(ea, qj));
  }
  MyMatrix<Tint> Rad = NullspaceMatMod(HermRowsToMatrix<Tint>(frob, d), ell_T);
  int dim_ss = d - static_cast<int>(Rad.rows());
  // a generator z of A / R: distinct irreducible factors of total degree
  // dim A / R in its minimal polynomial
  T_small_prime ell_u = static_cast<T_small_prime>(ell);
  // The candidates come from a linear congruential generator of fixed seed,
  // so that every process of an MPI computation gets the same primes, in
  // the same order, and the same line spaces.
  uint64_t state = 0x9E3779B97F4A7C15ULL;
  for (int iter = 0; iter < 1000; iter++) {
    std::vector<int64_t> z(d);
    for (int a = 0; a < d; a++) {
      state = state * 6364136223846793005ULL + 1442695040888963407ULL;
      z[a] = static_cast<int64_t>((state >> 33) % static_cast<uint64_t>(ell));
    }
    std::vector<std::vector<int64_t>> pows{A.one()};
    std::optional<PolyQ> opt_mu;
    while (!opt_mu) {
      pows.push_back(A.product(pows.back(), z));
      MyMatrix<Tint> NSP =
          NullspaceMatMod(HermRowsToMatrix<Tint>(pows, d), ell_T);
      if (NSP.rows() > 0) {
        int k = static_cast<int>(pows.size()) - 1;
        int64_t lead = order_residue(NSP(0, k), ell);
        int64_t inv = InverseModPrime(lead, ell);
        std::vector<T_residue> c(k + 1);
        for (int i = 0; i <= k; i++) {
          c[i] = static_cast<T_residue>(order_residue(NSP(0, i), ell) * inv %
                                        ell);
        }
        opt_mu = PolyQ(ell_u, c);
      }
    }
    std::vector<std::pair<PolyQ, int>> fac =
        FactorModQWithMultiplicity(*opt_mu);
    int tot = 0;
    for (auto &pr : fac) {
      tot += pr.first.degree();
    }
    if (tot != dim_ss) {
      continue;
    }
    std::vector<HermIdealPrime<Tint>> l_prime;
    MyVector<Tint> z_lift(d);
    for (int a = 0; a < d; a++) {
      z_lift(a) = UniversalScalarConversion<Tint, int64_t>(z[a]);
    }
    for (auto &[g, mult] : fac) {
      // the idempotent of the g-primary part of F_ell[z]
      PolyQ gk(ell_u, {1});
      for (int i = 0; i < mult; i++) {
        gk = gk * g;
      }
      PolyQ rest = (*opt_mu) / gk;
      PolyQ u = rest * PolyInverseModQ(rest % gk, gk);
      std::vector<int64_t> eps = A.poly_at(u, z);
      // eps A, of dimension e f
      std::vector<std::vector<int64_t>> rows;
      for (int a = 0; a < d; a++) {
        std::vector<int64_t> ea(d, 0);
        ea[a] = 1;
        rows.push_back(A.product(eps, ea));
      }
      int dim_loc = d - static_cast<int>(
                            NullspaceMatMod(HermRowsToMatrix<Tint>(rows, d),
                                            ell_T)
                                .rows());
      HermIdealPrime<Tint> P;
      P.ell = ell;
      P.f = g.degree();
      P.e = dim_loc / P.f;
      P.F = MakeFiniteField(ell_u, g);
      P.z = z_lift;
      // eps R, then the residues: a eps = sum_k c_k z^k eps + (eps R)
      std::vector<std::vector<int64_t>> eps_rad;
      for (int i = 0; i < Rad.rows(); i++) {
        std::vector<int64_t> r(d);
        for (int a = 0; a < d; a++) {
          r[a] = order_residue(Rad(i, a), ell);
        }
        eps_rad.push_back(A.product(eps, r));
      }
      std::vector<std::vector<int64_t>> zk_eps;
      std::vector<int64_t> zk = A.one();
      for (int k = 0; k < P.f; k++) {
        zk_eps.push_back(A.product(zk, eps));
        zk = A.product(zk, z);
      }
      for (int a = 0; a < d; a++) {
        std::vector<int64_t> ea(d, 0);
        ea[a] = 1;
        std::vector<std::vector<int64_t>> sys = zk_eps;
        for (auto &r : eps_rad) {
          sys.push_back(r);
        }
        sys.push_back(A.product(ea, eps));
        MyMatrix<Tint> NSP =
            NullspaceMatMod(HermRowsToMatrix<Tint>(sys, d), ell_T);
        int last = static_cast<int>(sys.size()) - 1;
        std::optional<std::vector<int64_t>> coef;
        for (int i = 0; i < NSP.rows() && !coef; i++) {
          int64_t lv = order_residue(NSP(i, last), ell);
          if (lv != 0) {
            int64_t inv = InverseModPrime(ell - lv, ell);
            std::vector<int64_t> c(P.f);
            for (int k = 0; k < P.f; k++) {
              c[k] = order_residue(NSP(i, k), ell) * inv % ell;
            }
            coef = c;
          }
        }
        if (!coef) {
          std::cerr << "HERMNB: no residue for a basis element\n";
          throw TerminalException{1};
        }
        P.res_basis.push_back(FqFromCoeffs(P.F, *coef));
      }
      // P: the kernel of the residue map, and ell O_E
      std::vector<std::vector<int64_t>> res_rows;
      for (int a = 0; a < d; a++) {
        std::vector<int64_t> r(P.f, 0);
        for (int k = 0; k < P.f; k++) {
          r[k] = P.res_basis[a].coeff(k);
        }
        res_rows.push_back(r);
      }
      MyMatrix<Tint> Ker =
          NullspaceMatMod(HermRowsToMatrix<Tint>(res_rows, P.f), ell_T);
      MyMatrix<Tint> Gen(Ker.rows() + d, d);
      for (int i = 0; i < Ker.rows(); i++) {
        for (int a = 0; a < d; a++) {
          Gen(i, a) = Ker(i, a);
        }
      }
      for (int a = 0; a < d; a++) {
        for (int b = 0; b < d; b++) {
          Gen(Ker.rows() + a, b) = (a == b) ? ell_T : Tint(0);
        }
      }
      P.ideal = GetZbasis(Gen);
      MyMatrix<Tint> Gen2(d * d, d);
      for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
          MyVector<Tint> pr = ord.product(GetMatrixRow(P.ideal, i),
                                          GetMatrixRow(P.ideal, j));
          for (int a = 0; a < d; a++) {
            Gen2(i * d + j, a) = pr(a);
          }
        }
      }
      P.ideal_sq = GetZbasis(Gen2);
      l_prime.push_back(P);
    }
    // the conjugate primes and the types
    for (size_t i = 0; i < l_prime.size(); i++) {
      HermIdealPrime<Tint> &P = l_prime[i];
      std::optional<size_t> opt_j;
      for (size_t j = 0; j < l_prime.size() && !opt_j; j++) {
        bool all = true;
        for (int r = 0; r < P.ideal.rows() && all; r++) {
          MyVector<Tint> c = ord.conjugate(GetMatrixRow(P.ideal, r));
          all = HermResidue(l_prime[j], c).is_zero();
        }
        if (all) {
          opt_j = j;
        }
      }
      if (!opt_j) {
        std::cerr << "HERMNB: the conjugate of a prime is not found\n";
        throw TerminalException{1};
      }
      P.conj_index = *opt_j;
      if (*opt_j != i) {
        P.type = HermIdealType::Split;
      } else {
        FqElt rz = HermResidue(P, P.z);
        FqElt rcz = HermResidue(P, ord.conjugate(P.z));
        P.type = (rz == rcz) ? HermIdealType::Ramified : HermIdealType::Inert;
      }
      P.name = "(" + std::to_string(ell) + "," + std::to_string(i) + ")";
    }
    return l_prime;
  }
  std::cerr << "HERMNB: no generator of O_E / rad modulo " << ell << "\n";
  throw TerminalException{1};
}

/*
  A Hermitian lattice L in the coordinates of a Z-basis: the trace form G
  and the multiplication matrices Omega_a of the basis e of O_E (y -> y
  Omega_a). h(y, x) on e is the solution z of Tr(e_a z) = y Omega_a G x^T.
 */
template <typename T, typename Tint> struct HermLatticeData {
  CMOrder<T, Tint> const *ord;
  MyMatrix<Tint> G;
  std::vector<MyMatrix<Tint>> Omega;
  MyMatrix<T> tr_inv;
  int n() const { return G.rows(); }
  MyMatrix<Tint> mult_of(MyVector<Tint> const &a) const {
    MyMatrix<Tint> M = ZeroMatrix<Tint>(n(), n());
    for (int k = 0; k < a.size(); k++) {
      if (a(k) != 0) {
        M += a(k) * Omega[k];
      }
    }
    return M;
  }
  MyVector<T> herm(MyVector<Tint> const &y, MyVector<Tint> const &x) const {
    int d = ord->deg;
    MyVector<Tint> Gx = G * x;
    MyVector<T> t(d);
    for (int a = 0; a < d; a++) {
      MyVector<Tint> yO = Omega[a].transpose() * y;
      t(a) = UniversalScalarConversion<T, Tint>(yO.dot(Gx));
    }
    return tr_inv * t;
  }
};

template <typename T, typename Tint>
HermLatticeData<T, Tint>
HermLatticeOfBasis(CMOrder<T, Tint> const &ord, RealifiedRing<T, Tint> const &ring,
                   MyMatrix<T> const &A) {
  HermLatticeData<T, Tint> L;
  L.ord = &ord;
  MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(ring.B);
  MyMatrix<T> G_T = A * B_T * A.transpose();
  L.G = UniversalMatrixConversion<Tint, T>(G_T);
  MyMatrix<T> Ainv = Inverse(A);
  for (auto &R : ring.l_ring_basis) {
    MyMatrix<T> O = A * UniversalMatrixConversion<T, Tint>(R) * Ainv;
    if (!IsIntegralMatrix(O)) {
      std::cerr << "HERMNB: the lattice is not an O_E-module\n";
      throw TerminalException{1};
    }
    L.Omega.push_back(UniversalMatrixConversion<Tint, T>(O));
  }
  int d = ord.deg;
  MyMatrix<T> Tr(d, d);
  for (int a = 0; a < d; a++) {
    for (int b = 0; b < d; b++) {
      MyVector<Tint> ea = ZeroVector<Tint>(d), eb = ZeroVector<Tint>(d);
      ea(a) = 1;
      eb(b) = 1;
      Tr(a, b) = UniversalScalarConversion<T, Tint>(
          ord.trace_of(ord.product(ea, eb)));
    }
  }
  L.tr_inv = Inverse(Tr);
  return L;
}

// Whether L is unimodular at every prime above ell.
template <typename T, typename Tint>
bool HermIsUnimodularAt(HermLatticeData<T, Tint> const &L, int64_t ell) {
  int d = L.ord->deg;
  int m = L.n() / d;
  Tint det = DeterminantMat(L.G);
  Tint disc = DeterminantMat(L.ord->trace_norm_form());
  Tint ell_T = UniversalScalarConversion<Tint, int64_t>(ell);
  int v_det = 0, v_disc = 0;
  while (det != 0 && det % ell_T == 0) {
    det /= ell_T;
    v_det++;
  }
  while (disc % ell_T == 0) {
    disc /= ell_T;
    v_disc++;
  }
  return v_det == m * v_disc;
}

/*
  The lines of L / C L over F_q = O_E / C: the rows of C L / ell L, then an
  F_q-basis v_j given as the z^t v_j over F_ell, inverted modulo ell. The
  interface is that of a projective space over F_q: codes of F_q, the
  lines of P^{m-1}(F_q) by ProjectiveLineOfIndex.
 */
struct HermLineSpace {
  int n = 0;
  int m = 0;
  int n_sub = 0;
  FiniteField F;
  std::vector<int64_t> basis;
  std::vector<int64_t> basis_inv;
  int64_t ell() const { return static_cast<int64_t>(F.p); }
  size_t n_line() const {
    return NumberProjectiveLines(m, static_cast<int>(F.q));
  }
  std::vector<uint64_t> coordinates(std::vector<int64_t> const &z) const {
    int64_t p = ell();
    std::vector<int64_t> y(n, 0);
    for (int j = 0; j < n; j++) {
      int64_t sum = 0;
      for (int k = 0; k < n; k++) {
        sum = (sum + z[k] * basis_inv[k * n + j]) % p;
      }
      y[j] = sum;
    }
    std::vector<uint64_t> c(m, 0);
    for (int j = 0; j < m; j++) {
      uint64_t code = 0, pw = 1;
      for (int t = 0; t < F.f; t++) {
        code += static_cast<uint64_t>(y[n_sub + j * F.f + t]) * pw;
        pw *= static_cast<uint64_t>(p);
      }
      c[j] = code;
    }
    return c;
  }
  std::optional<size_t> line_index(std::vector<int64_t> const &z) const {
    std::vector<uint64_t> c = coordinates(z);
    int ipos = 0;
    while (ipos < m && c[ipos] == 0) {
      ipos++;
    }
    if (ipos == m) {
      return {};
    }
    uint64_t scal = FqCodeInv(F, c[ipos]);
    std::vector<int> v(m);
    for (int j = 0; j < m; j++) {
      v[j] = static_cast<int>(FqCodeMul(F, c[j], scal));
    }
    return IndexOfProjectiveLine(v, static_cast<int>(F.q));
  }
  std::vector<int64_t> line_vector(size_t const &idx) const {
    int64_t p = ell();
    std::vector<int> c = ProjectiveLineOfIndex(m, static_cast<int>(F.q), idx);
    std::vector<int64_t> z(n, 0);
    for (int j = 0; j < m; j++) {
      uint64_t code = static_cast<uint64_t>(c[j]);
      for (int t = 0; t < F.f; t++) {
        int64_t a = static_cast<int64_t>(code % static_cast<uint64_t>(p));
        code /= static_cast<uint64_t>(p);
        if (a == 0) {
          continue;
        }
        int row = n_sub + j * F.f + t;
        for (int k = 0; k < n; k++) {
          z[k] = (z[k] + a * basis[row * n + k]) % p;
        }
      }
    }
    return z;
  }
};

// The Z-basis (rows, in the coordinates of L) of I L, I given by a Z-basis.
template <typename T, typename Tint>
MyMatrix<Tint> HermIdealTimesLattice(HermLatticeData<T, Tint> const &L,
                                     MyMatrix<Tint> const &I) {
  int n = L.n();
  MyMatrix<Tint> Gen(I.rows() * n, n);
  for (int r = 0; r < I.rows(); r++) {
    MyMatrix<Tint> M = L.mult_of(GetMatrixRow(I, r));
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        Gen(r * n + i, j) = M(i, j);
      }
    }
  }
  return GetZbasis(Gen);
}

template <typename T, typename Tint>
HermLineSpace HermGetLineSpace(HermLatticeData<T, Tint> const &L,
                               HermIdealPrime<Tint> const &C) {
  int n = L.n();
  int64_t ell = C.ell;
  Tint ell_T = UniversalScalarConversion<Tint, int64_t>(ell);
  HermLineSpace ls;
  ls.n = n;
  ls.F = C.F;
  MyMatrix<Tint> CL = HermIdealTimesLattice(L, C.ideal);
  // the rows of C L independent modulo ell
  std::vector<std::vector<int64_t>> rows;
  auto rank_of = [&](std::vector<std::vector<int64_t>> const &R) -> int {
    if (R.empty()) {
      return 0;
    }
    return static_cast<int>(R.size()) -
           static_cast<int>(
               NullspaceMatMod(HermRowsToMatrix<Tint>(R, n), ell_T).rows());
  };
  for (int i = 0; i < CL.rows(); i++) {
    std::vector<int64_t> r(n);
    for (int k = 0; k < n; k++) {
      r[k] = order_residue(CL(i, k), ell);
    }
    rows.push_back(r);
    if (rank_of(rows) < static_cast<int>(rows.size())) {
      rows.pop_back();
    }
  }
  ls.n_sub = static_cast<int>(rows.size());
  int f = C.f;
  if ((n - ls.n_sub) % f != 0) {
    std::cerr << "HERMNB: L / C L has no F_q-structure\n";
    throw TerminalException{1};
  }
  ls.m = (n - ls.n_sub) / f;
  MyMatrix<Tint> Mz = L.mult_of(C.z);
  for (int k = 0; k < n && static_cast<int>(rows.size()) < n; k++) {
    std::vector<std::vector<int64_t>> cand = rows;
    MyVector<Tint> v = ZeroVector<Tint>(n);
    v(k) = 1;
    for (int t = 0; t < f; t++) {
      std::vector<int64_t> r(n);
      for (int j = 0; j < n; j++) {
        r[j] = order_residue(v(j), ell);
      }
      cand.push_back(r);
      v = Mz.transpose() * v;
    }
    if (rank_of(cand) == static_cast<int>(cand.size())) {
      rows = cand;
    }
  }
  if (static_cast<int>(rows.size()) != n) {
    std::cerr << "HERMNB: no F_q-basis of L / C L\n";
    throw TerminalException{1};
  }
  ls.basis.resize(n * n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      ls.basis[i * n + j] = rows[i][j];
    }
  }
  std::optional<std::vector<int64_t>> opt =
      InverseMatModInt64(ls.basis, n, ell, true);
  if (!opt) {
    std::cerr << "HERMNB: the adapted basis is not invertible\n";
    throw TerminalException{1};
  }
  ls.basis_inv = *opt;
  return ls;
}

/*
  The fractional ideal C^{-1} = ell^{-1} { b in O_E : b C in ell O_E }, as
  the rows of a Z-basis on e.
 */
template <typename T, typename Tint>
MyMatrix<T> HermIdealInverse(CMOrder<T, Tint> const &ord,
                             HermIdealPrime<Tint> const &C) {
  int d = ord.deg;
  int64_t ell = C.ell;
  Tint ell_T = UniversalScalarConversion<Tint, int64_t>(ell);
  MyMatrix<Tint> Cond(d, d * C.ideal.rows());
  for (int a = 0; a < d; a++) {
    MyVector<Tint> ea = ZeroVector<Tint>(d);
    ea(a) = 1;
    for (int r = 0; r < C.ideal.rows(); r++) {
      MyVector<Tint> pr = ord.product(ea, GetMatrixRow(C.ideal, r));
      for (int c = 0; c < d; c++) {
        Cond(a, r * d + c) = pr(c);
      }
    }
  }
  MyMatrix<Tint> Ker = NullspaceMatMod(Cond, ell_T);
  MyMatrix<Tint> Gen(Ker.rows() + d, d);
  for (int r = 0; r < Ker.rows(); r++) {
    for (int a = 0; a < d; a++) {
      Gen(r, a) = Ker(r, a);
    }
  }
  for (int a = 0; a < d; a++) {
    for (int b = 0; b < d; b++) {
      Gen(Ker.rows() + a, b) = (a == b) ? ell_T : Tint(0);
    }
  }
  MyMatrix<Tint> D = GetZbasis(Gen);
  return UniversalMatrixConversion<T, Tint>(D) /
         UniversalScalarConversion<T, int64_t>(ell);
}

/*
  The data of the neighbors at P: the prime C of the line space, the
  uniformizer, the Z-basis of C^{-1} (rows on e, rational).
 */
template <typename T, typename Tint> struct HermNeighborPrime {
  HermIdealPrime<Tint> P;
  HermIdealPrime<Tint> C;
  MyVector<Tint> pi;
  MyMatrix<T> C_inv;
};

template <typename T, typename Tint>
HermNeighborPrime<T, Tint>
HermGetNeighborPrime(CMOrder<T, Tint> const &ord,
                     std::vector<HermIdealPrime<Tint>> const &l_prime, size_t i) {
  int d = ord.deg;
  HermNeighborPrime<T, Tint> np;
  np.P = l_prime[i];
  np.C = l_prime[np.P.conj_index];
  int64_t ell = np.P.ell;
  Tint ell_T = UniversalScalarConversion<Tint, int64_t>(ell);
  std::vector<HermIdealPrime<Tint>> l_other;
  for (size_t j = 0; j < l_prime.size(); j++) {
    if (j != i) {
      l_other.push_back(l_prime[j]);
    }
  }
  if (np.P.type == HermIdealType::Inert) {
    // a uniformizer of p in K_0: P cap K_0, the elements fixed by conj
    MyMatrix<Tint> Pc = np.P.ideal;
    MyMatrix<Tint> Diff(Pc.rows(), d);
    for (int r = 0; r < Pc.rows(); r++) {
      MyVector<Tint> v = GetMatrixRow(Pc, r);
      MyVector<Tint> w = ord.conjugate(v) - v;
      for (int a = 0; a < d; a++) {
        Diff(r, a) = w(a);
      }
    }
    MyMatrix<Tint> Ker = NullspaceIntMat(Diff);
    MyMatrix<Tint> Fixed = Ker * Pc;
    np.pi = HermFindUniformizer<T, Tint>(Fixed, np.P, l_other);
  } else {
    np.pi = HermFindUniformizer<T, Tint>(np.P.ideal, np.P, l_other);
  }
  np.C_inv = HermIdealInverse(ord, np.C);
  return np;
}

/*
  The neighbors of L at P for the line of the residue vector z of L / C L,
  as Z-bases in the coordinates of L.
 */
template <typename T, typename Tint>
std::vector<MyMatrix<T>>
HermNeighbors(HermLatticeData<T, Tint> const &L,
              HermNeighborPrime<T, Tint> const &np,
              std::vector<int64_t> const &z) {
  CMOrder<T, Tint> const &ord = *L.ord;
  int n = L.n();
  int d = ord.deg;
  int64_t ell = np.P.ell;
  Tint ell_T = UniversalScalarConversion<Tint, int64_t>(ell);
  HermIdealPrime<Tint> const &P = np.P;
  HermIdealPrime<Tint> const &C = np.C;
  MyVector<Tint> w(n);
  for (int k = 0; k < n; k++) {
    w(k) = UniversalScalarConversion<Tint, int64_t>(z[k]);
  }
  MyMatrix<Tint> Mpi = L.mult_of(np.pi);
  // the vectors x of the neighbors
  std::vector<MyVector<Tint>> l_x;
  if (P.type == HermIdealType::Split) {
    l_x.push_back(Mpi.transpose() * w);
  } else {
    MyVector<T> hww = L.herm(w, w);
    if (!HermResidueRat(P, hww).is_zero()) {
      return {};
    }
    // a basis vector e_k with h(e_k, w) not in P
    std::optional<int> opt_k;
    FqElt lambda;
    for (int k = 0; k < n && !opt_k; k++) {
      MyVector<Tint> ek = ZeroVector<Tint>(n);
      ek(k) = 1;
      FqElt r = HermResidueRat(P, L.herm(ek, w));
      if (!r.is_zero()) {
        opt_k = k;
        lambda = r;
      }
    }
    if (!opt_k) {
      std::cerr << "HERMNB: w is in the radical, L is not unimodular\n";
      throw TerminalException{1};
    }
    // the elements s of O_E / P, as polynomials in z
    uint64_t q = P.F.q;
    std::optional<FqElt> target;
    int half = P.f / 2;
    if (P.type == HermIdealType::Inert) {
      // Tr(s lambda) = -h(w, w) / pi modulo P
      MyMatrix<T> Mpi_e = UniversalMatrixConversion<T, Tint>(ord.mult_matrix(np.pi));
      MyVector<T> hww_pi = Inverse(Mpi_e).transpose() * hww;
      target = FqNeg(P.F, HermResidueRat(P, hww_pi));
    }
    MyMatrix<Tint> Mz = L.mult_of(P.z);
    for (uint64_t code = 0; code < q; code++) {
      FqElt s = FqDecode(P.F, code);
      if (target) {
        FqElt sl = FqMul(P.F, s, lambda);
        FqElt tr = FqAdd(P.F, sl, FqFrobenius(P.F, sl, half));
        if (!(tr == *target)) {
          continue;
        }
      }
      // pi s e_k, s = sum_t s_t z^t
      MyVector<Tint> v = ZeroVector<Tint>(n);
      v(*opt_k) = 1;
      MyVector<Tint> sv = ZeroVector<Tint>(n);
      for (int t = 0; t < P.f; t++) {
        int64_t c = s.coeff(t);
        if (c != 0) {
          sv += UniversalScalarConversion<Tint, int64_t>(c) * v;
        }
        v = Mz.transpose() * v;
      }
      l_x.push_back(MyVector<Tint>(w + Mpi.transpose() * sv));
    }
  }
  std::vector<MyMatrix<T>> l_basis;
  for (auto const &x : l_x) {
    // M: the kernel of y -> h(x, y) = conj(h(y, x)) modulo C
    std::vector<std::vector<int64_t>> func;
    for (int k = 0; k < n; k++) {
      MyVector<Tint> ek = ZeroVector<Tint>(n);
      ek(k) = 1;
      MyVector<T> hyx = L.herm(ek, x);
      MyVector<T> hxy =
          UniversalMatrixConversion<T, Tint>(ord.conj).transpose() * hyx;
      FqElt r = HermResidueRat(C, hxy);
      std::vector<int64_t> row(C.f, 0);
      for (int t = 0; t < C.f; t++) {
        row[t] = r.coeff(t);
      }
      func.push_back(row);
    }
    MyMatrix<Tint> Ker = NullspaceMatMod(HermRowsToMatrix<Tint>(func, C.f), ell_T);
    MyMatrix<T> Gens(Ker.rows() + n + np.C_inv.rows(), n);
    int pos = 0;
    for (int i = 0; i < Ker.rows(); i++) {
      for (int k = 0; k < n; k++) {
        Gens(pos, k) = UniversalScalarConversion<T, Tint>(Ker(i, k));
      }
      pos++;
    }
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        Gens(pos, k) = (i == k) ? UniversalScalarConversion<T, int64_t>(ell) : T(0);
      }
      pos++;
    }
    MyVector<T> x_T = UniversalVectorConversion<T, Tint>(x);
    for (int r = 0; r < np.C_inv.rows(); r++) {
      MyMatrix<T> Mc = ZeroMatrix<T>(n, n);
      for (int a = 0; a < d; a++) {
        if (np.C_inv(r, a) != 0) {
          Mc += np.C_inv(r, a) * UniversalMatrixConversion<T, Tint>(L.Omega[a]);
        }
      }
      MyVector<T> cx = Mc.transpose() * x_T;
      for (int k = 0; k < n; k++) {
        Gens(pos, k) = cx(k);
      }
      pos++;
    }
    MyMatrix<T> Basis = GetZbasis(Gens);
#ifdef SANITY_CHECK_HERMITIAN_NEIGHBOR
    if (Basis.rows() != n || T_abs(DeterminantMat(Basis)) != 1) {
      std::cerr << "HERMNB: the neighbor does not have the covolume of L\n";
      throw TerminalException{1};
    }
    MyMatrix<T> Binv = Inverse(Basis);
    for (auto &O : L.Omega) {
      if (!IsIntegralMatrix(MyMatrix<T>(
              Basis * UniversalMatrixConversion<T, Tint>(O) * Binv))) {
        std::cerr << "HERMNB: the neighbor is not an O_E-module\n";
        throw TerminalException{1};
      }
    }
    if (!IsIntegralMatrix(MyMatrix<T>(
            Basis * UniversalMatrixConversion<T, Tint>(L.G) *
            Basis.transpose()))) {
      std::cerr << "HERMNB: the neighbor is not integral\n";
      throw TerminalException{1};
    }
#endif
    l_basis.push_back(Basis);
  }
  return l_basis;
}

// clang-format off
#endif  // SRC_GENERA_HERMITIANNEIGHBOR_H_
// clang-format on

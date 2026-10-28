// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ORDERLATTICENEIGHBOR_H_
#define SRC_GENERA_ORDERLATTICENEIGHBOR_H_

// clang-format off
#include "RealifiedRing.h"
#include "GenusEnumeration.h"
#include "MAT_MatrixMod.h"
#include "MAT_MatrixModInt64.h"
#include "SiegelSpecialValues.h"
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ORDER_LATTICE
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ORDER_LATTICE
#endif

/*
  Lattices over R = O_K (K imaginary quadratic) with a Hermitian form, and
  their Kneser neighbours at a prime ideal P of O_K.

  A lattice L is held by a Z-basis, the rows of a matrix A in the realified
  coordinates of RealifiedRing.h. In the coordinates of that basis the form
  is G = A B A^T and the multiplication by omega is Omega = A W A^{-1}, an
  integral matrix since L is an O_K-module. The Hermitian form is recovered
  from the pair (G, G Omega^T): for z = h(y, x) = a + b omega,
    t1 = y G x^T       = Tr z       = 2 a + s b,
    t2 = y G Omega^T x^T = Tr(conj(omega) z) = s a + 2 n b,
  with omega^2 = s omega - n, so that a = (2 n t1 - s t2) / |D| and b = (2 t2 -
  s t1) / |D|, D = s^2 - 4 n. The same pair is the family whose canonical
  form and automorphism group are those of the Hermitian lattice (see
  LatticeStabEquiCan.h).

  The quotient L / P L is a vector space of dimension m over O_K / P = F_q,
  and its lines are indexed through an adapted basis:
  - P = (ell, omega - r) split (q = ell): L / ell L is the sum of the two
    eigenspaces of Omega mod ell, for r and conj(r) = s - r, and the
    projection onto the first one along the second is L / P L,
  - P = (ell) inert (q = ell^2): a basis e_1, ..., e_m of L / ell L over F_q,
    that is e_j, e_j Omega a basis over F_ell, the element a + b omega of F_q
    being coded as the integer a + ell b,
  - over R = Z, P = (ell): L / ell L itself.
  The lines are then those of F_q^m, indexed by ProjectiveLineOfIndex and
  IndexOfProjectiveLine of GenusEnumeration.h, which only need the codes of
  0 and 1 to be 0 and 1.

  The neighbour of L at P for the line of x:
  - split: x adjusted modulo P L so that ell divides h(x, x), M = {y in L :
    h(y, x) in conj(P)} and L' = M + P^{-1} x, P^{-1} = conj(P) / ell, which
    is Z x + Z (omega - conj(r)) x / ell. Every line gives a neighbour, and
    there are (q^m - 1)/(q - 1) of them.
  - inert: x isotropic modulo ell, adjusted modulo ell L so that ell^2
    divides h(x, x), M = {y in L : h(y, x) in ell O_K} and L' = M + O_K x /
    ell. The admissible lifts x + ell y form an affine line over F_ell modulo
    M, so each isotropic line gives ell neighbours, and T(ell) has degree
    ell (ell^{2m-3} + ...) (ell (ell^3 + 1) for m = 3).
  - over Z (L even, h(x, x) = Q(x) = B(x, x) / 2): x isotropic modulo ell,
    adjusted modulo ell L so that ell^2 divides Q(x), M = {y in L : B(y, x)
    = 0 mod ell} and L' = M + Z x / ell, the classical Kneser neighbour.
  They are unimodular at ell when L is, and in the genus of L. The ramified
  primes are not treated.
 */

// The number types are T (a field, mpq_class or fmpq_class) and Tint (its
// ring of integers); int64_t goes through ScalarFromInt64.
template <typename Tring> Tring amf_int(int64_t const &x) { return ScalarFromInt64<Tring>(x); }

/*
  A prime ideal P of R: the rational prime ell under it and the residue
  degree f; for f = 1 over O_K the residue r of omega (omega = r mod P).
 */
struct OrderPrime {
  std::string name;
  int64_t ell = 0;
  int f = 1;
  int64_t r = 0;
  int64_t q() const { return f == 1 ? ell : ell * ell; }
};

inline bool order_is_prime(int64_t const &x) {
  if (x < 2) {
    return false;
  }
  for (int64_t d = 2; d * d <= x; d++) {
    if (x % d == 0) {
      return false;
    }
  }
  return true;
}

/*
  The prime of R named by e: over Z an integer prime; over O_K an integer
  prime inert in O_K, or a generator a + b omega (b != 0) of a split prime.
 */
template <typename T, typename Tint>
OrderPrime get_order_prime(RealifiedRing<T, Tint> const &ring,
                                  RingElementName const &e) {
  OrderPrime P;
  P.name = e.name;
  if (ring.d == 1 || e.b == 0) {
    if (e.b != 0) {
      std::cerr << "AMF: the prime " << e.name << " should be an integer over Z\n";
      throw TerminalException{1};
    }
    P.ell = e.a < 0 ? -e.a : e.a;
    if (!order_is_prime(P.ell)) {
      std::cerr << "AMF: " << e.name << " is not a prime\n";
      throw TerminalException{1};
    }
    if (ring.d == 1) {
      P.f = 1;
      return P;
    }
    // inert iff x^2 - s x + n has no root mod ell (and ell does not ramify)
    int64_t D = ring.s * ring.s - 4 * ring.p;
    if (ResidueModPrime(D, P.ell) == 0) {
      std::cerr << "AMF: the prime " << e.name << " ramifies, not treated\n";
      throw TerminalException{1};
    }
    for (int64_t x = 0; x < P.ell; x++) {
      int64_t val = ResidueModPrime(x * x - ring.s * x + ring.p, P.ell);
      if (val == 0) {
        std::cerr << "AMF: " << e.name << " splits in O_K, give a generator "
                  << "a+bw of a prime over it\n";
        throw TerminalException{1};
      }
    }
    P.f = 2;
    return P;
  }
  int64_t norm = e.a * e.a + ring.s * e.a * e.b + ring.p * e.b * e.b;
  if (!order_is_prime(norm)) {
    std::cerr << "AMF: the norm " << norm << " of " << e.name
              << " is not a prime (a split prime of degree one is needed)\n";
    throw TerminalException{1};
  }
  int64_t D = ring.s * ring.s - 4 * ring.p;
  if (ResidueModPrime(D, norm) == 0) {
    std::cerr << "AMF: the prime " << e.name << " ramifies, not treated\n";
    throw TerminalException{1};
  }
  P.ell = norm;
  P.f = 1;
  // a + b omega divides omega - r iff a + b r = 0 mod ell
  P.r = ResidueModPrime(-ResidueModPrime(e.a, norm) *
                            InverseModPrime(ResidueModPrime(e.b, norm), norm),
                        norm);
#ifdef SANITY_CHECK_ORDER_LATTICE
  if (ResidueModPrime(P.r * P.r - ring.s * P.r + ring.p, norm) != 0) {
    std::cerr << "AMF: the residue of omega is not a root\n";
    throw TerminalException{1};
  }
#endif
  return P;
}

/*
  The primes given by name, a comma separated list: an integer, a generator
  a+bw (see get_order_prime), or ell@r for the prime (ell, omega - r) of
  degree one, which need not be principal.
 */
template <typename T, typename Tint>
std::vector<OrderPrime> parse_order_primes(RealifiedRing<T, Tint> const &ring, std::string const &str) {
  std::vector<OrderPrime> l_prime;
  for (auto &tok : STRING_Split(str, ",")) {
    size_t pos = tok.find('@');
    if (pos == std::string::npos) {
      l_prime.push_back(get_order_prime(ring, parse_ring_elements(tok)[0]));
      continue;
    }
    OrderPrime P;
    P.name = tok;
    P.ell = ParseScalar<int64_t>(tok.substr(0, pos));
    P.r = ResidueModPrime(ParseScalar<int64_t>(tok.substr(pos + 1)), P.ell);
    P.f = 1;
    int64_t D = ring.s * ring.s - 4 * ring.p;
    if (ring.d != 2 || !order_is_prime(P.ell) || ResidueModPrime(D, P.ell) == 0 ||
        ResidueModPrime(P.r * P.r - ring.s * P.r + ring.p, P.ell) != 0) {
      std::cerr << "AMF: " << tok << " is not a split prime (ell, omega - r) of O_K\n";
      throw TerminalException{1};
    }
    l_prime.push_back(P);
  }
  return l_prime;
}

template <typename Tint> int64_t order_residue(Tint const &x, int64_t const &ell) {
  Tint ell_T = ScalarFromInt64<Tint>(ell);
  Tint r = x % ell_T;
  if (r < 0) {
    r += ell_T;
  }
  return UniversalScalarConversion<int64_t, Tint>(r);
}

// The residue of a rational whose denominator is prime to ell.
template <typename T> int64_t order_residue_rat(T const &x, int64_t const &ell) {
  int64_t num = order_residue(GetNumerator_z(x), ell);
  int64_t den = order_residue(GetDenominator_z(x), ell);
  if (den == 0) {
    std::cerr << "AMF: a denominator is divisible by " << ell << "\n";
    throw TerminalException{1};
  }
  return (num * InverseModPrime(den, ell)) % ell;
}

/*
  The arithmetic of F_q = O_K / P on the codes: for f = 1 the residues mod
  ell, for f = 2 the code a + ell b of a + b omega, omega^2 = s omega - n.
 */
struct ResidueField {
  int64_t ell;
  int f;
  int64_t s;
  int64_t n;
  int64_t q() const { return f == 1 ? ell : ell * ell; }
  int64_t mul(int64_t const &x, int64_t const &y) const {
    if (f == 1) {
      return (x * y) % ell;
    }
    int64_t a1 = x % ell, b1 = x / ell, a2 = y % ell, b2 = y / ell;
    int64_t bb = (b1 * b2) % ell;
    int64_t a = ResidueModPrime(a1 * a2 - n * bb, ell);
    int64_t b = ResidueModPrime(a1 * b2 + a2 * b1 + s * bb, ell);
    return a + ell * b;
  }
  int64_t inv(int64_t const &x) const {
    if (f == 1) {
      return InverseModPrime(x, ell);
    }
    // conj(a + b omega) = (a + s b) - b omega, divided by the norm
    int64_t a = x % ell, b = x / ell;
    int64_t norm = ResidueModPrime(a * a + s * a % ell * b + n * b % ell * b, ell);
    int64_t norm_inv = InverseModPrime(norm, ell);
    int64_t ca = ResidueModPrime((a + s * b) * norm_inv, ell);
    int64_t cb = ResidueModPrime(-b * norm_inv, ell);
    return ca + ell * cb;
  }
};

// The rank mod ell of the rows of M.
template <typename Tint>
int order_rank_mod(MyMatrix<Tint> const &M, int64_t const &ell) {
  MyMatrix<Tint> NSP = NullspaceMatMod(M, amf_int<Tint>(ell));
  return M.rows() - NSP.rows();
}

/*
  The lines of L / P L for a lattice of matrix Omega (empty over Z): the
  adapted basis and its inverse modulo ell, and the conversions between
  the vectors of L / ell L (residues) and the normalized coordinates in
  F_q^m.
 */
struct LineSpace {
  int n = 0;
  int m = 0;
  ResidueField F;
  // row major n x n, mod ell: the rows e_1, ..., e_m (split, Z), or e_1,
  // e_1 Omega, ..., e_m, e_m Omega (inert), then the complement
  std::vector<int64_t> basis;
  std::vector<int64_t> basis_inv;
  size_t n_line() const { return NumberProjectiveLines(m, static_cast<int>(F.q())); }
  // The coordinates in F_q^m of the residue vector z (not normalized).
  std::vector<int64_t> coordinates(std::vector<int64_t> const &z) const {
    std::vector<int64_t> y(n, 0);
    for (int j = 0; j < n; j++) {
      int64_t sum = 0;
      for (int k = 0; k < n; k++) {
        sum = (sum + z[k] * basis_inv[k * n + j]) % F.ell;
      }
      y[j] = sum;
    }
    std::vector<int64_t> c(m);
    for (int j = 0; j < m; j++) {
      c[j] = F.f == 1 ? y[j] : y[2 * j] + F.ell * y[2 * j + 1];
    }
    return c;
  }
  // The index of the line of the residue vector z, or nothing for z in P L.
  std::optional<size_t> line_index(std::vector<int64_t> const &z) const {
    std::vector<int64_t> c = coordinates(z);
    int ipos = 0;
    while (ipos < m && c[ipos] == 0) {
      ipos++;
    }
    if (ipos == m) {
      return {};
    }
    int64_t scal = F.inv(c[ipos]);
    std::vector<int> v(m);
    for (int j = 0; j < m; j++) {
      v[j] = static_cast<int>(F.mul(c[j], scal));
    }
    return IndexOfProjectiveLine(v, static_cast<int>(F.q()));
  }
  // A residue vector of L / ell L on the line of index idx.
  std::vector<int64_t> line_vector(size_t const &idx) const {
    std::vector<int> c = ProjectiveLineOfIndex(m, static_cast<int>(F.q()), idx);
    std::vector<int64_t> z(n, 0);
    for (int j = 0; j < m; j++) {
      int64_t a = F.f == 1 ? c[j] : c[j] % F.ell;
      int64_t b = F.f == 1 ? 0 : c[j] / F.ell;
      int row_a = F.f == 1 ? j : 2 * j;
      for (int k = 0; k < n; k++) {
        z[k] = (z[k] + a * basis[row_a * n + k]) % F.ell;
        if (b != 0) {
          z[k] = (z[k] + b * basis[(row_a + 1) * n + k]) % F.ell;
        }
      }
    }
    return z;
  }
};

template <typename T>
std::vector<int64_t> order_residue_vector(MyVector<T> const &v, int64_t const &ell) {
  std::vector<int64_t> z(v.size());
  for (int k = 0; k < v.size(); k++) {
    z[k] = order_residue_rat(v(k), ell);
  }
  return z;
}

template <typename T, typename Tint>
LineSpace get_line_space(RealifiedRing<T, Tint> const &ring,
                                MyMatrix<Tint> const &Omega,
                                OrderPrime const &P) {
  LineSpace ls;
  int n = ring.n();
  ls.n = n;
  ls.m = ring.m;
  ls.F = ResidueField{P.ell, P.f, ring.s, ring.p};
  int64_t ell = P.ell;
  MyMatrix<Tint> Bas(n, n);
  if (ring.d == 1) {
    Bas = IdentityMat<Tint>(n);
  } else if (P.f == 1) {
    // the eigenspaces of Omega mod ell for r and conj(r) = s - r
    int64_t r_conj = ResidueModPrime(ring.s - P.r, ell);
    MyMatrix<Tint> Id = IdentityMat<Tint>(n);
    MyMatrix<Tint> V1 = NullspaceMatMod(MyMatrix<Tint>(Omega - amf_int<Tint>(P.r) * Id), amf_int<Tint>(ell));
    MyMatrix<Tint> V2 = NullspaceMatMod(MyMatrix<Tint>(Omega - amf_int<Tint>(r_conj) * Id), amf_int<Tint>(ell));
    if (V1.rows() != ring.m || V2.rows() != ring.m) {
      std::cerr << "AMF: the eigenspaces of omega mod " << ell << " have dimensions " << V1.rows() << ", "
                << V2.rows() << " instead of " << ring.m << "\n";
      throw TerminalException{1};
    }
    for (int j = 0; j < ring.m; j++) {
      for (int k = 0; k < n; k++) {
        Bas(j, k) = V1(j, k);
        Bas(ring.m + j, k) = V2(j, k);
      }
    }
  } else {
    int n_found = 0;
    for (int k = 0; k < n && n_found < ring.m; k++) {
      MyMatrix<Tint> Test(2 * n_found + 1, n);
      for (int i = 0; i < 2 * n_found; i++) {
        for (int j = 0; j < n; j++) {
          Test(i, j) = Bas(i, j);
        }
      }
      for (int j = 0; j < n; j++) {
        Test(2 * n_found, j) = (j == k) ? 1 : 0;
      }
      if (order_rank_mod(Test, ell) == 2 * n_found + 1) {
        for (int j = 0; j < n; j++) {
          Bas(2 * n_found, j) = (j == k) ? 1 : 0;
          Bas(2 * n_found + 1, j) = Omega(k, j);
        }
        n_found++;
      }
    }
    if (n_found != ring.m) {
      std::cerr << "AMF: no basis of L / ell L over F_q found\n";
      throw TerminalException{1};
    }
  }
  ls.basis.resize(n * n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      ls.basis[i * n + j] = order_residue(Bas(i, j), ell);
    }
  }
  std::optional<std::vector<int64_t>> opt = InverseMatModInt64(ls.basis, n, ell, true);
  if (!opt) {
    std::cerr << "AMF: the adapted basis is not invertible mod " << ell << "\n";
    throw TerminalException{1};
  }
  ls.basis_inv = *opt;
  return ls;
}

/*
  The Hermitian product h(e_k, x) for the basis vectors e_k of L, as the
  pairs (a_k, b_k), h = a + b omega; over Z simply B(e_k, x) (b_k = 0).
 */
template <typename T, typename Tint>
std::vector<std::pair<Tint, Tint>>
order_hermitian_products(RealifiedRing<T, Tint> const &ring, MyMatrix<Tint> const &G,
                         MyMatrix<Tint> const &Omega, MyVector<Tint> const &x) {
  int n = G.rows();
  MyVector<Tint> t1 = G * x;
  std::vector<std::pair<Tint, Tint>> l_h(n);
  if (ring.d == 1) {
    for (int k = 0; k < n; k++) {
      l_h[k] = {t1(k), Tint(0)};
    }
    return l_h;
  }
  MyVector<Tint> xO = Omega.transpose() * x;
  MyVector<Tint> t2 = G * xO;
  Tint absD = 4 * ring.p - ring.s * ring.s;
  for (int k = 0; k < n; k++) {
    Tint a_num = 2 * ring.p * t1(k) - ring.s * t2(k);
    Tint b_num = 2 * t2(k) - ring.s * t1(k);
#ifdef SANITY_CHECK_ORDER_LATTICE
    if (a_num % absD != 0 || b_num % absD != 0) {
      std::cerr << "AMF: the Hermitian product is not in O_K\n";
      throw TerminalException{1};
    }
#endif
    l_h[k] = {a_num / absD, b_num / absD};
  }
  return l_h;
}

// h(x, x) = G(x, x) / 2 for R = O_K.
template <typename Tint>
Tint order_hermitian_norm(MyMatrix<Tint> const &G, MyVector<Tint> const &x) {
  Tint val = x.dot(G * x);
#ifdef SANITY_CHECK_ORDER_LATTICE
  if (val % 2 != 0) {
    std::cerr << "AMF: Tr h(x, x) is odd\n";
    throw TerminalException{1};
  }
#endif
  return val / 2;
}

/*
  The neighbours of L at P for the residue vector z of a line of L / P L,
  each as the rows of a Z-basis in the coordinates of L: one for a split
  prime or over Z, ell for an inert prime (see below), none when the line
  is not isotropic.
 */
template <typename T, typename Tint>
std::vector<MyMatrix<T>>
order_neighbor(RealifiedRing<T, Tint> const &ring, MyMatrix<Tint> const &G,
               MyMatrix<Tint> const &Omega, OrderPrime const &P, std::vector<int64_t> const &z) {
  // Over Z (an even lattice, h(x, x) = Q(x) = B(x, x) / 2) the construction
  // is the one of an inert prime with the single functional B(y, x): L' = M
  // + Z x / ell, M = {y : B(y, x) = 0 mod ell}.
  bool split = ring.d == 2 && P.f == 1;
  int n = G.rows();
  int64_t ell = P.ell;
  Tint ell_T = amf_int<Tint>(ell);
  MyVector<Tint> x(n);
  for (int k = 0; k < n; k++) {
    x(k) = amf_int<Tint>(z[k]);
  }
  MyMatrix<Tint> Id = IdentityMat<Tint>(n);
  if (split) {
    // x modulo P L, P L generated by ell L and the rows of Omega - r.
    Tint hx = order_hermitian_norm(G, x);
    int64_t hx_res = order_residue(hx, ell);
    if (hx_res != 0) {
      MyMatrix<Tint> OmR = Omega - amf_int<Tint>(P.r) * Id;
      bool done = false;
      for (int k = 0; k < n && !done; k++) {
        MyVector<Tint> y = GetMatrixRow(OmR, k);
        int64_t bxy = order_residue(Tint(x.dot(G * y)), ell);
        if (bxy != 0) {
          int64_t c = ResidueModPrime(-hx_res * InverseModPrime(bxy, ell), ell);
          x += amf_int<Tint>(c) * y;
          done = true;
        }
      }
      if (!done) {
        std::cerr << "AMF: the line cannot be lifted to h(x, x) = 0 mod ell\n";
        throw TerminalException{1};
      }
    }
  } else {
    Tint hx = order_hermitian_norm(G, x);
    if (order_residue(hx, ell) != 0) {
      return {};
    }
    // x + ell c e_k: h = h(x, x) + ell c Tr h(x, e_k) mod ell^2
    int64_t hx_red = order_residue(Tint(hx / ell_T), ell);
    if (hx_red != 0) {
      MyVector<Tint> Gx = G * x;
      bool done = false;
      for (int k = 0; k < n && !done; k++) {
        int64_t bxe = order_residue(Gx(k), ell);
        if (bxe != 0) {
          int64_t c = ResidueModPrime(-hx_red * InverseModPrime(bxe, ell), ell);
          x(k) += ell_T * amf_int<Tint>(c);
          done = true;
        }
      }
      if (!done) {
        std::cerr << "AMF: L is not unimodular at " << ell << "\n";
        throw TerminalException{1};
      }
    }
  }
#ifdef SANITY_CHECK_ORDER_LATTICE
  {
    Tint hx = order_hermitian_norm(G, x);
    Tint mod = split ? ell_T : ell_T * ell_T;
    if (hx % mod != 0) {
      std::cerr << "AMF: the adjusted vector is not isotropic enough\n";
      throw TerminalException{1};
    }
  }
#endif
  // M: the kernel of y -> h(y, x) mod conj(P) (split) or mod ell (inert, Z).
  std::vector<std::pair<Tint, Tint>> l_h = order_hermitian_products(ring, G, Omega, x);
  int n_func = ring.d == 1 || split ? 1 : 2;
  MyMatrix<Tint> Func(n, n_func);
  int64_t r_conj = ResidueModPrime(ring.s - P.r, ell);
  for (int k = 0; k < n; k++) {
    int64_t a = order_residue(l_h[k].first, ell);
    int64_t b = order_residue(l_h[k].second, ell);
    if (split) {
      // omega = conj(r) mod conj(P)
      Func(k, 0) = amf_int<Tint>((a + b * r_conj) % ell);
    } else {
      Func(k, 0) = amf_int<Tint>(a);
      if (ring.d == 2) {
        Func(k, 1) = amf_int<Tint>(b);
      }
    }
  }
  MyMatrix<Tint> NSP = NullspaceMatMod(Func, ell_T);
  int n_gen = NSP.rows() + n + 2;
  MyMatrix<T> Gens(n_gen, n);
  int pos = 0;
  for (int i = 0; i < NSP.rows(); i++) {
    for (int k = 0; k < n; k++) {
      Gens(pos, k) = UniversalScalarConversion<T, Tint>(NSP(i, k));
    }
    pos++;
  }
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < n; k++) {
      Gens(pos, k) = (i == k) ? amf_int<T>(ell) : T(0);
    }
    pos++;
  }
  /*
    The lifts of the line. Split prime, Z: x alone. Inert prime: L / M is
    O / ell = F_{ell^2} through y -> h(y, x), and the lifts x + ell y with
    ell^2 | h(x + ell y, x + ell y) form an affine line, y being determined
    modulo M up to the elements with Tr h(x, y) = 0 mod ell. So there are ell
    of them, x + ell t y1 for t in F_ell, with h(y1, x) = 2 omega - s mod ell
    (of trace 0), and as many neighbours of the line.
   */
  std::vector<MyVector<Tint>> l_lift{x};
  if (!split && ring.d == 2) {
    // y1 = c1 e_k1 + c2 e_k2 for two basis vectors whose h(e_k, x) are
    // independent mod ell
    int64_t target_a = ResidueModPrime(-ring.s, ell), target_b = 2 % ell;
    std::optional<MyVector<Tint>> opt_y1;
    for (int k1 = 0; k1 < n && !opt_y1; k1++) {
      for (int k2 = k1 + 1; k2 < n && !opt_y1; k2++) {
        int64_t a1 = order_residue(l_h[k1].first, ell), b1 = order_residue(l_h[k1].second, ell);
        int64_t a2 = order_residue(l_h[k2].first, ell), b2 = order_residue(l_h[k2].second, ell);
        int64_t det = ResidueModPrime(a1 * b2 - a2 * b1, ell);
        if (det == 0) {
          continue;
        }
        int64_t det_inv = InverseModPrime(det, ell);
        int64_t c1 = ResidueModPrime((target_a * b2 - a2 * target_b) % ell * det_inv, ell);
        int64_t c2 = ResidueModPrime((a1 * target_b - target_a * b1) % ell * det_inv, ell);
        MyVector<Tint> y1 = ZeroVector<Tint>(n);
        y1(k1) = amf_int<Tint>(c1);
        y1(k2) = amf_int<Tint>(c2);
        opt_y1 = y1;
      }
    }
    if (!opt_y1) {
      std::cerr << "AMF: y -> h(y, x) is not onto O / ell\n";
      throw TerminalException{1};
    }
    for (int64_t t = 1; t < ell; t++) {
      Tint coef = ell_T * amf_int<Tint>(t);
      l_lift.push_back(MyVector<Tint>(x + coef * (*opt_y1)));
    }
  }
  std::vector<MyMatrix<T>> l_basis;
  for (auto &xl : l_lift) {
#ifdef SANITY_CHECK_ORDER_LATTICE
    {
      Tint hxl = order_hermitian_norm(G, xl);
      Tint mod = split ? ell_T : ell_T * ell_T;
      if (hxl % mod != 0) {
        std::cerr << "AMF: a lift is not isotropic enough\n";
        throw TerminalException{1};
      }
    }
#endif
    // P^{-1} x: x and (omega - conj(r)) x / ell (split); x / ell and omega x
    // / ell (inert); x / ell (Z).
    MyVector<Tint> xO = ring.d == 2 ? MyVector<Tint>(Omega.transpose() * xl) : xl;
    for (int k = 0; k < n; k++) {
      if (split) {
        Gens(pos, k) = UniversalScalarConversion<T, Tint>(xl(k));
        Gens(pos + 1, k) = UniversalScalarConversion<T, Tint>(xO(k) - amf_int<Tint>(r_conj) * xl(k)) / amf_int<T>(ell);
      } else {
        Gens(pos, k) = UniversalScalarConversion<T, Tint>(xl(k)) / amf_int<T>(ell);
        Gens(pos + 1, k) = UniversalScalarConversion<T, Tint>(xO(k)) / amf_int<T>(ell);
      }
    }
    MyMatrix<T> Basis = GetZbasis(Gens);
#ifdef SANITY_CHECK_ORDER_LATTICE
    if (Basis.rows() != n || T_abs(DeterminantMat(Basis)) != 1) {
      std::cerr << "AMF: the neighbour does not have the covolume of L\n";
      throw TerminalException{1};
    }
#endif
    l_basis.push_back(Basis);
  }
  return l_basis;
}

// clang-format off
#endif  // SRC_GENERA_ORDERLATTICENEIGHBOR_H_
// clang-format on

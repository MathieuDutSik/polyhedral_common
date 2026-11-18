// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_HERMITIANLATTICE_H_
#define SRC_GENERA_HERMITIANLATTICE_H_

// clang-format off
#include "HermitianGenus.h"
#include "GenusSymbol.h"
#include "MAT_MatrixInt.h"
#include "MAT_MatrixPolynomial.h"
#include "ClassicLLL.h"
#include <cmath>
#include <deque>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
// clang-format on

/*
  Hermitian lattices over E = Q(zeta_d) seen as pairs (L, sigma): L a
  Z-lattice (Gram matrix G on Z^n) and sigma an isometry of minimal
  polynomial Phi_d acting on the rows (x -> x P), zeta acting as sigma.

  The local structure. At a prime l, the l-part of D(L) = L^# / L is a
  module over O_E (x) Z_l = prod_{P | l} O_P, a sum over the P of cyclic
  modules O_P / P^k. Its multiplicities r_{P,k} are read from the orders of
  kernels on D(L) of polynomials in sigma (lattice intersections):
  - a uniformizer at all the P above l: l when l does not divide d, 1 -
    zeta^m for d = l^a m (zeta^m is a primitive l^a-th root of unity, and
    E / Q(zeta_{l^a}) is not ramified above l);
  - the P-part: the kernel of a high power of g(zeta'), for zeta' = zeta^(l^a)
    a primitive m-th root of unity and g a factor of Phi_m mod l (a unit at
    the other primes above l, of positive valuation at P), when Phi_m has
    several factors mod l;
  then |D_P[P^i]| = q^(sum_k min(i, k) r_k), q the size of the residue
  field of P.

  The invariant of a pair: the genus symbol of L (Conway-Sloane) and, per
  prime l of det L, the multiset over the P of the multiplicities r_{P,k}.
  Two pairs of the same hermitian genus have the same invariant; whether
  the invariant separates the genera is tested on the pieces of Hecke.jl.
 */

// d with d = 2 mod 4 replaced by d / 2, and the matrix of a primitive
// root of unity of that order: sigma, or -sigma for d = 2 mod 4.
template <typename T>
MyMatrix<T> HermPrimitiveRoot(MyMatrix<T> const &P, int64_t d) {
  if (d % 4 == 2) {
    return -P;
  }
  return P;
}



/*
  The order of the kernel on D(L) = L^# / L of the endomorphisms x -> x M_t
  (M_t integral, invertible over Q, preserving L): with x = y G^-1, the y of
  Z^n with y G^-1 M_t integral, a lattice containing Z^n G, of index the
  order sought.
 */
// A polynomial mod l (coefficients in [0, l)) as one over Z.
inline std::vector<int64_t> HermLiftPoly(PolyQ const &g) {
  std::vector<int64_t> res;
  for (auto c : g) {
    res.push_back(static_cast<int64_t>(c));
  }
  return res;
}

/*
  A basis (rows) of the lattice of the y in the lattice of basis S (rows,
  rational) with y A_t integral for all t (A_t invertible over Q).
 */
template <typename T, typename Tint>
MyMatrix<T> HermIntegralPreimage(MyMatrix<T> const &S,
                                 std::vector<MyMatrix<T>> const &As) {
  int n = S.cols();
  std::vector<T> ListDenS;
  for (int i = 0; i < S.rows(); i++) {
    for (int j = 0; j < n; j++) {
      ListDenS.push_back(GetDenominator(S(i, j)));
    }
  }
  T denS = LCMlist(ListDenS);
  // K: the current lattice scaled by denS, integral
  MyMatrix<Tint> K = UniversalMatrixConversion<Tint, T>(MyMatrix<T>(S * denS));
  for (auto const &A : As) {
    // {y : y A in Z^n} = Z^n A^-1, scaled by denS: denS Z^n A^-1
    MyMatrix<T> Ainv = Inverse(A);
    std::vector<T> ListDen;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        ListDen.push_back(GetDenominator(Ainv(i, j) * denS));
      }
    }
    T den = LCMlist(ListDen);
    MyMatrix<Tint> B =
        UniversalMatrixConversion<Tint, T>(MyMatrix<T>(Ainv * denS * den));
    Tint den_i = UniversalScalarConversion<Tint, T>(den);
    MyMatrix<Tint> I = IntersectionLattice(MyMatrix<Tint>(K * den_i), B);
    for (int i = 0; i < I.rows(); i++) {
      for (int j = 0; j < n; j++) {
        I(i, j) /= den_i;
      }
    }
    K = I;
  }
  return UniversalMatrixConversion<T, Tint>(K) / denS;
}

template <typename T, typename Tint>
Tint HermKernelOrder(MyMatrix<T> const &G, std::vector<MyMatrix<T>> const &Ms) {
  // the x of L^# (basis G^-1) with x M integral, over L = Z^n
  MyMatrix<T> K = HermIntegralPreimage<T, Tint>(Inverse(G), Ms);
  return UniversalScalarConversion<Tint, T>(T(1) / T_abs(DeterminantMat(K)));
}

// The valuation of the positive integer x at l.
template <typename Tint> int HermValuation(Tint x, int64_t l) {
  int v = 0;
  Tint l_i(static_cast<long>(l));
  while (ResInt(x, l_i) == 0) {
    x /= l_i;
    v++;
  }
  return v;
}

/*
  The multiplicities r_k (k = 1, 2, ...) of the O_P-module D_P from the
  log_q of the orders of the D_P[P^i], i = 1, 2, ..., until they stop
  growing: log_q |D_P[P^i]| = sum_k min(i, k) r_k.
 */
inline std::vector<int> HermMultiplicities(std::vector<int> const &logs) {
  // first differences: number of summands with k >= i
  std::vector<int> ge;
  int prev = 0;
  for (int x : logs) {
    ge.push_back(x - prev);
    prev = x;
  }
  std::vector<int> r;
  for (size_t i = 0; i < ge.size(); i++) {
    int next = (i + 1 < ge.size()) ? ge[i + 1] : 0;
    r.push_back(ge[i] - next);
  }
  while (!r.empty() && r.back() == 0) {
    r.pop_back();
  }
  return r;
}

/*
  Per prime l of det G, the multiplicities r_{P,k} of D(L) at the primes P
  above l, as a sorted list over the P.
 */
template <typename T, typename Tint>
std::map<int64_t, std::vector<std::vector<int>>>
HermitianLocalStructure(MyMatrix<T> const &G, MyMatrix<T> const &P,
                        int64_t d) {
  int n = G.rows();
  int64_t dn = NormalizeCyclotomicIndex(d);
  MyMatrix<T> Z = HermPrimitiveRoot(P, d);
  Tint detG = UniversalScalarConversion<Tint, T>(T_abs(DeterminantMat(G)));
  std::map<int64_t, std::vector<std::vector<int>>> res;
  if (detG == 1) {
    return res;
  }
  for (auto const &[l_i, mult] : FactorsIntMap(detG)) {
    int64_t l = UniversalScalarConversion<int64_t, Tint>(l_i);
    int64_t la = 1, m = dn;
    while (m % l == 0) {
      m /= l;
      la *= l;
    }
    // the uniformizer at the primes above l
    MyMatrix<T> U;
    if (la == 1) {
      U = T(static_cast<long>(l)) * IdentityMat<T>(n);
    } else {
      U = IdentityMat<T>(n) - MatrixPower(Z, m);
    }
    // the cuts to the primes P above l
    std::vector<MyMatrix<T>> cuts;
    int64_t f = MultiplicativeOrder(l, m);
    if (m > 2 && EulerPhi(m) / f > 1) {
      MyMatrix<T> Zp = MatrixPower(Z, la);
      int vmax = static_cast<int>(EulerPhi(la) * static_cast<int64_t>(mult));
      for (auto const &g : HermCyclotomicFactorsModL(m, l)) {
        cuts.push_back(MatrixPower(EvaluatePolynomialMatrix(HermLiftPoly(g), Zp), vmax + 1));
      }
    } else {
      cuts.push_back(IdentityMat<T>(n) * T(0));
    }
    int64_t q = 1;
    for (int64_t k = 0; k < f; k++) {
      q *= l;
    }
    std::vector<std::vector<int>> structs;
    for (auto const &cut : cuts) {
      std::vector<int> logs;
      MyMatrix<T> Ui = IdentityMat<T>(n);
      int last = -1;
      // the P-adic lengths are at most e v_l(det), e = phi(l^a) the
      // ramification index of P over l
      int imax = static_cast<int>(EulerPhi(la) * static_cast<int64_t>(mult)) + 1;
      for (int i = 1; i <= imax; i++) {
        Ui = Ui * U;
        std::vector<MyMatrix<T>> Ms{Ui};
        bool has_cut = !IsZeroMatrix(cut);
        if (has_cut) {
          Ms.push_back(cut);
        }
        Tint ord = HermKernelOrder<T, Tint>(G, Ms);
        int v = HermValuation(ord, l);
        if (v % static_cast<int>(f) != 0) {
          std::cerr << "HERM_LATTICE: the kernel order is not a power of q\n";
          throw TerminalException{1};
        }
        int lg = v / static_cast<int>(f);
        logs.push_back(lg);
        if (lg == last) {
          break;
        }
        last = lg;
      }
      structs.push_back(HermMultiplicities(logs));
    }
    std::sort(structs.begin(), structs.end());
    res[l] = structs;
  }
  return res;
}

/*
  The uniformizer at the primes above l and the cuts to each of them
  (empty: a single prime), the cuts of power cut_power.
 */
template <typename T> struct HermPrimeCuts {
  MyMatrix<T> U;
  std::vector<MyMatrix<T>> cuts;
};

template <typename T>
HermPrimeCuts<T> HermGetPrimeCuts(MyMatrix<T> const &Z, int64_t dn, int64_t l,
                                  int cut_power) {
  int n = Z.rows();
  int64_t la = 1, m = dn;
  while (m % l == 0) {
    m /= l;
    la *= l;
  }
  HermPrimeCuts<T> res;
  if (la == 1) {
    res.U = T(static_cast<long>(l)) * IdentityMat<T>(n);
  } else {
    res.U = IdentityMat<T>(n) - MatrixPower(Z, m);
  }
  int64_t f = MultiplicativeOrder(l, m);
  if (m > 2 && EulerPhi(m) / f > 1) {
    MyMatrix<T> Zp = MatrixPower(Z, la);
    for (auto const &g : HermCyclotomicFactorsModL(m, l)) {
      res.cuts.push_back(MatrixPower(EvaluatePolynomialMatrix(HermLiftPoly(g), Zp), cut_power));
    }
  }
  return res;
}

/*
  Representatives of the elements of B / S, B and S lattices (bases as
  rows, rational) with S in B.
 */
template <typename T, typename Tint>
std::vector<MyVector<T>> HermQuotientReps(MyMatrix<T> const &B,
                                          MyMatrix<T> const &S) {
  int n = B.rows();
  MyMatrix<T> C = S * Inverse(B);
  MyMatrix<Tint> C_i = UniversalMatrixConversion<Tint, T>(C);
  ResultSmithNormalForm<Tint> snf = SmithNormalForm(C_i);
  MyMatrix<Tint> D = snf.ROW * C_i * snf.COL;
  MyMatrix<T> Bp = Inverse(UniversalMatrixConversion<T, Tint>(snf.COL)) * B;
  // the orders of the cyclic factors, enumerated so bounded in size
  std::vector<size_t> ord;
  Tint total(1);
  for (int i = 0; i < n; i++) {
    Tint o = T_abs(D(i, i));
    total *= o;
    if (total > Tint(static_cast<long>(1) << 24)) {
      std::cerr << "HERM_LATTICE: the quotient is too large to enumerate\n";
      throw TerminalException{1};
    }
    ord.push_back(UniversalScalarConversion<size_t, Tint>(o));
  }
  std::vector<MyVector<T>> res;
  std::vector<size_t> c(n, 0);
  while (true) {
    MyVector<T> v = ZeroVector<T>(B.cols());
    for (int i = 0; i < n; i++) {
      if (c[i] != 0) {
        for (int j = 0; j < B.cols(); j++) {
          v(j) += UniversalScalarConversion<T, size_t>(c[i]) * Bp(i, j);
        }
      }
    }
    res.push_back(v);
    int i = 0;
    while (i < n) {
      c[i]++;
      if (c[i] < ord[i]) {
        break;
      }
      c[i] = 0;
      i++;
    }
    if (i == n) {
      break;
    }
  }
  return res;
}


// The pair on the lattice of basis W (rows, in the coordinates of L),
// LLL reduced.
template <typename T, typename Tint>
std::pair<MyMatrix<T>, MyMatrix<T>>
HermRestrictPair(MyMatrix<T> const &G, MyMatrix<T> const &P,
                 MyMatrix<T> const &W) {
  MyMatrix<T> G1 = W * G * W.transpose();
  MyMatrix<T> P1 = W * P * Inverse(W);
  LLLreduction<T, Tint> red = LLLreducedBasis<T, Tint>(G1, std::cerr);
  MyMatrix<T> R = UniversalMatrixConversion<T, Tint>(red.Pmat);
  return {R * G1 * R.transpose(), R * P1 * Inverse(R)};
}

template <typename T, typename Tint>
std::string HermitianPairInvariant(MyMatrix<T> const &G, MyMatrix<T> const &P,
                                   int64_t d) {
  std::ostringstream os;
  os << GenusSymbolString(GenusSymbolOfGram<T, Tint>(G));
  for (auto const &[l, structs] : HermitianLocalStructure<T, Tint>(G, P, d)) {
    os << " | " << l << ":";
    for (auto const &s : structs) {
      os << " [";
      for (size_t k = 0; k < s.size(); k++) {
        os << (k > 0 ? "," : "") << s[k];
      }
      os << "]";
    }
  }
  return os.str();
}

/*
  The neighbors of (L, sigma) at the primes above l: the overlattices L +
  O_E v, v in D(L) killed by a prime P above l (integral and even), and
  the sublattices L' with L / L' = O_E / P, the duals of the overlattices
  L^# + O_E w of L^#, w killed by P modulo L^#.
 */
template <typename T, typename Tint>
std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>>
HermNeighbors(MyMatrix<T> const &G, MyMatrix<T> const &P, int64_t d,
              int64_t l, bool over, bool under, bool first_only) {
  int n = G.rows();
  int64_t dn = NormalizeCyclotomicIndex(d);
  MyMatrix<T> Z = HermPrimitiveRoot(P, d);
  int phi = static_cast<int>(EulerPhi(dn));
  HermPrimeCuts<T> pc = HermGetPrimeCuts(Z, dn, l, 1);
  std::vector<MyMatrix<T>> parts = pc.cuts;
  if (parts.empty()) {
    parts.push_back(MyMatrix<T>());
  }
  MyMatrix<T> Ginv = Inverse(G);
  MyMatrix<T> Id = IdentityMat<T>(n);
  std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>> res;
  auto orbit_rows = [&](MyMatrix<T> const &base,
                        MyVector<T> const &v) -> MyMatrix<T> {
    MyMatrix<T> M(base.rows() + phi, n);
    for (int i = 0; i < base.rows(); i++) {
      for (int j = 0; j < n; j++) {
        M(i, j) = base(i, j);
      }
    }
    MyVector<T> w = v;
    for (int k = 0; k < phi; k++) {
      for (int j = 0; j < n; j++) {
        M(base.rows() + k, j) = w(j);
      }
      w = P.transpose() * w;
    }
    return M;
  };
  // keeps the lattice of basis W when even, true when the search may stop
  auto keep_if_even = [&](MyMatrix<T> const &W) -> bool {
    MyMatrix<T> G1 = W * G * W.transpose();
    bool ok = IsIntegralMatrix(G1);
    for (int i = 0; ok && i < n; i++) {
      ok = IsInteger(G1(i, i) / T(2));
    }
    if (ok) {
      res.push_back(HermRestrictPair<T, Tint>(G, P, W));
    }
    return ok && first_only;
  };
  for (auto const &cut : parts) {
    if (over) {
      // B = {x in L^# : x U in L, x cut in L}, in the coordinates of L
      std::vector<MyMatrix<T>> As{pc.U};
      if (cut.rows() > 0) {
        As.push_back(cut);
      }
      MyMatrix<T> B = HermIntegralPreimage<T, Tint>(Ginv, As);
      for (auto const &v : HermQuotientReps<T, Tint>(B, Id)) {
        if (v.isZero()) {
          continue;
        }
        if (keep_if_even(GetZbasis(orbit_rows(Id, v)))) {
          return res;
        }
      }
    }
    if (under) {
      // B = {w in L^# U^-1 : w cut in L^#} (x in L^# iff x G is
      // integral), over L^#
      std::vector<MyMatrix<T>> As;
      if (cut.rows() > 0) {
        As.push_back(MyMatrix<T>(cut * G));
      }
      MyMatrix<T> B = HermIntegralPreimage<T, Tint>(
          MyMatrix<T>(Ginv * Inverse(pc.U)), As);
      for (auto const &w : HermQuotientReps<T, Tint>(B, Ginv)) {
        if (w.isZero()) {
          continue;
        }
        MyMatrix<T> Wd = GetZbasis(orbit_rows(Ginv, w));
        // the dual of the lattice of basis Wd: (Wd G Wd^T)^-1 Wd
        if (keep_if_even(GetZbasis(MyMatrix<T>(
                Inverse(MyMatrix<T>(Wd * G * Wd.transpose())) * Wd)))) {
          return res;
        }
      }
    }
  }
  return res;
}


/*
  The pair of O_E^b with the hermitian form diag(1, ..., 1, a), a = sum_k
  c_k t^k in K_0 (t = zeta + zeta^-1), trace form Tr_{E/Q}: on the basis
  zeta^k of each copy of O_E, Tr(a_i zeta^(k - l)); sigma the
  multiplication by zeta (by -zeta for d = 2 mod 4).
 */
template <typename T>
std::pair<MyMatrix<T>, MyMatrix<T>>
HermStandardPair(int64_t d, int b, std::vector<int64_t> const &a) {
  int64_t dn = NormalizeCyclotomicIndex(d);
  int phi = static_cast<int>(EulerPhi(dn));
  int n = b * phi;
  // a as a Laurent polynomial in zeta: sum_k c_k (zeta + zeta^-1)^k
  std::map<int64_t, int64_t> lau;
  for (size_t k = 0; k < a.size(); k++) {
    // binomial expansion
    int64_t binom = 1;
    for (int64_t s = 0; s <= static_cast<int64_t>(k); s++) {
      lau[static_cast<int64_t>(k) - 2 * s] += a[k] * binom;
      binom = binom * (static_cast<int64_t>(k) - s) / (s + 1);
    }
  }
  auto tr = [&](bool with_a, int64_t j) -> int64_t {
    if (!with_a) {
      return RamanujanSum(dn, j);
    }
    int64_t s = 0;
    for (auto const &[e, c] : lau) {
      s += c * RamanujanSum(dn, j + e);
    }
    return s;
  };
  MyMatrix<T> G = ZeroMatrix<T>(n, n);
  for (int blk = 0; blk < b; blk++) {
    bool with_a = (blk == b - 1);
    for (int k = 0; k < phi; k++) {
      for (int l = 0; l < phi; l++) {
        G(blk * phi + k, blk * phi + l) = T(static_cast<long>(tr(with_a, k - l)));
      }
    }
  }
  std::vector<int64_t> cyc = CyclotomicPolynomial(static_cast<int>(dn));
  MyMatrix<T> Z = ZeroMatrix<T>(n, n);
  for (int blk = 0; blk < b; blk++) {
    int o = blk * phi;
    for (int k = 0; k + 1 < phi; k++) {
      Z(o + k, o + k + 1) = 1;
    }
    for (int k = 0; k < phi; k++) {
      Z(o + phi - 1, o + k) = T(static_cast<long>(-cyc[k]));
    }
  }
  if (d % 4 == 2) {
    Z = -Z;
  }
  return {G, Z};
}

/*
  Starts for the search: the a = c_0 + c_1 t + c_2 t^2 (t = zeta + zeta^-1,
  |c_i| <= bound, distinct as elements of K_0 of degree phi(d) / 2) that are
  totally positive (positive at t = 2 cos(2 pi k / d), k prime to d) with a
  norm N_{K_0/Q}(a) whose primes divide 2 d N: a determinant class -1 at
  another prime, inert, is no piece's (a piece is unimodular there), so the
  other a are dropped by the descent anyway. Sorted by norm, at most
  max_starts of them. They give the hermitian spaces of the determinant
  classes at the primes of 2 d N.
 */
inline std::vector<std::vector<int64_t>>
HermDefaultStarts(int64_t d, int64_t N_small, int64_t bound,
                  size_t max_starts) {
  int64_t dn = NormalizeCyclotomicIndex(d);
  int64_t deg = EulerPhi(dn) / 2;
  int max_deg = static_cast<int>(std::min<int64_t>(2, deg - 1));
  std::vector<double> conj;
  for (int64_t k = 1; 2 * k < dn || (dn <= 2 && k == 1); k++) {
    if (std::gcd(k, dn) == 1) {
      conj.push_back(2 * std::cos(2 * M_PI * static_cast<double>(k) /
                                  static_cast<double>(dn)));
    }
  }
  int64_t rel = 2 * dn * N_small;
  std::vector<std::pair<int64_t, std::vector<int64_t>>> cand;
  std::set<std::vector<int64_t>> seen;
  for (int64_t c0 = 1; c0 <= bound; c0++) {
    for (int64_t c1 = (max_deg >= 1 ? -bound : 0);
         c1 <= (max_deg >= 1 ? bound : 0); c1++) {
      for (int64_t c2 = (max_deg >= 2 ? -bound : 0);
           c2 <= (max_deg >= 2 ? bound : 0); c2++) {
        double norm = 1;
        bool pos = true;
        for (double t : conj) {
          double v = static_cast<double>(c0) + static_cast<double>(c1) * t +
                     static_cast<double>(c2) * t * t;
          pos = pos && v > 1e-9;
          norm *= v;
        }
        if (!pos) {
          continue;
        }
        int64_t nm = static_cast<int64_t>(std::llround(norm));
        // the primes of the norm must divide 2 d N
        bool relevant = true;
        if (nm > 1) {
          for (auto const &[p, mult] : FactorsIntMap(nm)) {
            relevant = relevant && (rel % p == 0);
          }
        }
        if (!relevant) {
          continue;
        }
        std::vector<int64_t> a{c0};
        if (max_deg >= 1) {
          a.push_back(c1);
        }
        if (max_deg >= 2) {
          a.push_back(c2);
        }
        while (a.size() > 1 && a.back() == 0) {
          a.pop_back();
        }
        if (seen.insert(a).second) {
          cand.push_back({nm, a});
        }
      }
    }
  }
  std::sort(cand.begin(), cand.end());
  std::vector<std::vector<int64_t>> res;
  for (auto const &[nm, a] : cand) {
    if (res.size() < max_starts) {
      res.push_back(a);
    }
  }
  return res;
}


/*
  The pieces: one pair (L, sigma) per invariant (so per hermitian genus,
  see above) of rank b over O_E, even, with det L dividing N, in the
  hermitian spaces of the standard pairs of the starts a (forms diag(1,
  ..., 1, a)). One pair per invariant is enough: the neighbors of the pairs
  of a genus have the same invariants, locally at the prime of the step.
  - At a prime l not dividing N the pieces are maximal even lattices
    locally (det prime to l), and these form one local genus in a given
    space: any chain of overlattice neighbors at l ends there. Done
    greedily for each start, a start that cannot lose l being dropped.
  - At the primes of N: a piece X lies in a maximal even lattice M, and the
    lattices between X and M have determinants dividing det X. So from any
    lattice of the space, overlattice neighbors lead to the genus of M and
    sublattice neighbors from there to that of X, all of determinant
    dividing N once below it: a search by neighbors both ways, among the
    pairs of determinant dividing N, after a greedy descent to it.
 */
template <typename T, typename Tint>
std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>>
HermitianPiecesSearch(int64_t d, int b, Tint const &N,
                      std::vector<std::vector<int64_t>> const &starts,
                      std::ostream &os) {
  using Tpair = std::pair<MyMatrix<T>, MyMatrix<T>>;
  auto det_of = [](MyMatrix<T> const &G) -> Tint {
    return UniversalScalarConversion<Tint, T>(T_abs(DeterminantMat(G)));
  };
  std::vector<int64_t> primes;
  if (N > 1) {
    for (auto const &[l_i, mult] : FactorsIntMap(N)) {
      primes.push_back(UniversalScalarConversion<int64_t, Tint>(l_i));
    }
  }
  // the greedy descent: one overlattice neighbor at a time at a prime l of
  // the determinant, while l does not divide N or the power of l in det
  // exceeds that in N
  auto descend = [&](Tpair pr) -> std::optional<Tpair> {
    while (true) {
      Tint det = det_of(pr.first);
      if (ResInt(N, det) == 0) {
        return pr;
      }
      bool moved = false;
      for (auto const &[l_i, mult] : FactorsIntMap(det)) {
        int64_t l = UniversalScalarConversion<int64_t, Tint>(l_i);
        int vN = HermValuation(N, l);
        if (static_cast<int>(mult) <= vN) {
          continue;
        }
        std::vector<Tpair> nbs = HermNeighbors<T, Tint>(
            pr.first, pr.second, d, l, true, false, true);
        if (!nbs.empty()) {
          pr = nbs[0];
          moved = true;
          break;
        }
      }
      if (!moved) {
        return {};
      }
    }
  };
  std::map<std::string, Tpair> visited;
  std::deque<Tpair> queue;
  size_t n_dropped = 0;
  // the invariants are taken after the descent, where the determinants
  // are small (those of the starts are large, and so are the powers of the
  // cuts)
  for (auto const &a : starts) {
    std::optional<Tpair> opt = descend(HermStandardPair<T>(d, b, a));
    if (!opt) {
      n_dropped++;
      continue;
    }
    std::string key = HermitianPairInvariant<T, Tint>(opt->first, opt->second, d);
    if (visited.emplace(key, *opt).second) {
      queue.push_back(*opt);
    }
  }
  os << "HERM_PIECES: d=" << d << " b=" << b << " N=" << N << ": "
     << starts.size() << " starts, " << n_dropped << " dropped, "
     << visited.size() << " invariants after the descent\n";
  while (!queue.empty()) {
    Tpair pr = queue.front();
    queue.pop_front();
    for (int64_t l : primes) {
      for (auto const &nb :
           HermNeighbors<T, Tint>(pr.first, pr.second, d, l, true, true,
                                  false)) {
        if (ResInt(N, det_of(nb.first)) != 0) {
          continue;
        }
        std::string key =
            HermitianPairInvariant<T, Tint>(nb.first, nb.second, d);
        if (visited.emplace(key, nb).second) {
          queue.push_back(nb);
        }
      }
    }
  }
  std::vector<Tpair> res;
  for (auto const &[key, pr] : visited) {
    res.push_back(pr);
  }
  os << "HERM_PIECES: " << res.size() << " pieces\n";
  return res;
}

/*
  The pieces of rank b over Z[zeta_d] with even trace lattice of
  determinant dividing N, from the default starts (coefficients up to 4,
  at most 40 of them, enough in all the tests).
 */
template <typename T, typename Tint>
std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>>
HermitianPieces(int64_t d, int b, Tint const &N, std::ostream &os) {
  int64_t N_small = UniversalScalarConversion<int64_t, Tint>(N);
  return HermitianPiecesSearch<T, Tint>(
      d, b, N, HermDefaultStarts(d, N_small, 4, 40), os);
}

// clang-format off
#endif  // SRC_GENERA_HERMITIANLATTICE_H_
// clang-format on

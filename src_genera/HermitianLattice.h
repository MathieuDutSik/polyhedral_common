// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_HERMITIANLATTICE_H_
#define SRC_GENERA_HERMITIANLATTICE_H_

// clang-format off
#include "HermitianGenus.h"
#include "IsotypicPieces.h"
#include "GenusSymbol.h"
#include "MAT_MatrixInt.h"
#include "MAT_MatrixPolynomial.h"
#include "ClassicLLL.h"
#include <cmath>
#include <deque>
#include <functional>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
# define SANITY_CHECK_HERMITIAN_LATTICE
#endif

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
  The uniformizer at the primes above l and the cuts to each of them
  (empty: a single prime), the cuts of power cut_power.
 */
template <typename T> struct HermPrimeCuts {
  MyMatrix<T> U;
  std::vector<MyMatrix<T>> cuts;
  // the residue degree f and the ramification index e of the primes
  // above l, d = l^a m: f = ord(l mod m), e = phi(l^a)
  int64_t f;
  int64_t e;
};

template <typename T>
HermPrimeCuts<T> HermGetPrimeCuts(MyMatrix<T> const &Z, int64_t dn, int64_t l,
                                  int cut_power) {
  int n = Z.rows();
  int64_t la = HermPrimePowerPart(dn, l);
  int64_t m = dn / la;
  HermPrimeCuts<T> res;
  if (la == 1) {
    res.U = T(static_cast<long>(l)) * IdentityMat<T>(n);
  } else {
    res.U = IdentityMat<T>(n) - MatrixPower(Z, m);
  }
  res.f = MultiplicativeOrder(l, m);
  res.e = EulerPhi(la);
  if (m > 2 && EulerPhi(m) / res.f > 1) {
    MyMatrix<T> Zp = MatrixPower(Z, la);
    for (auto const &g : HermCyclotomicFactorsModL(m, l)) {
      res.cuts.push_back(MatrixPower(EvaluatePolynomialMatrix(HermLiftPoly(g), Zp), cut_power));
    }
  }
  return res;
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
    // the P-adic lengths are at most e v_l(det): a cut of that power kills
    // the parts at the other primes
    int64_t e = EulerPhi(HermPrimePowerPart(dn, l));
    int imax = static_cast<int>(e * static_cast<int64_t>(mult)) + 1;
    HermPrimeCuts<T> pc = HermGetPrimeCuts(Z, dn, l, imax);
    std::vector<std::optional<MyMatrix<T>>> cuts;
    for (auto const &c : pc.cuts) {
      cuts.push_back(c);
    }
    if (cuts.empty()) {
      cuts.push_back(std::nullopt);
    }
    int f = static_cast<int>(pc.f);
    std::vector<std::vector<int>> structs;
    for (auto const &cut : cuts) {
      std::vector<int> logs;
      MyMatrix<T> Ui = IdentityMat<T>(n);
      int last = -1;
      for (int i = 1; i <= imax; i++) {
        Ui = Ui * pc.U;
        std::vector<MyMatrix<T>> Ms{Ui};
        if (cut) {
          Ms.push_back(*cut);
        }
        Tint ord = HermKernelOrder<T, Tint>(G, Ms);
        int v = HermValuation(ord, l);
#ifdef SANITY_CHECK_HERMITIAN_LATTICE
        if (v % f != 0) {
          std::cerr << "HERM_LATTICE: the kernel order is not a power of q\n";
          throw TerminalException{1};
        }
#endif
        int lg = v / f;
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
  Representatives of the elements of B / S, B and S lattices (bases as
  rows, rational) with S in B, given in turn to f until it returns true.
 */
template <typename T, typename Tint, typename F>
void HermQuotientRepsForEach(MyMatrix<T> const &B, MyMatrix<T> const &S,
                             F const &f) {
  int n = B.rows();
  MyMatrix<T> C = S * Inverse(B);
  MyMatrix<Tint> C_i = UniversalMatrixConversion<Tint, T>(C);
  ResultSmithNormalForm<Tint> snf = SmithNormalForm(C_i);
  MyMatrix<Tint> D = snf.ROW * C_i * snf.COL;
  MyMatrix<T> Bp = Inverse(UniversalMatrixConversion<T, Tint>(snf.COL)) * B;
  // the orders of the cyclic factors (the enumeration may stop early, so
  // only bounded to fit the counters)
  std::vector<size_t> ord;
  Tint total(1);
  for (int i = 0; i < n; i++) {
    Tint o = T_abs(D(i, i));
    total *= o;
    if (total > Tint(static_cast<long>(1) << 62)) {
      std::cerr << "HERM_LATTICE: the quotient is too large to enumerate\n";
      throw TerminalException{1};
    }
    ord.push_back(UniversalScalarConversion<size_t, Tint>(o));
  }
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
    if (f(v)) {
      return;
    }
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
  // L + O_E v (v in L^#) is even exactly when q(v) is even and the
  // b(v, sigma^k v), k < phi, are integers: the b(v, sigma^k v) for the
  // other k are integral combinations of these, and q(lambda v) for lambda
  // in O_E is then even
  auto over_even = [&](MyVector<T> const &v) -> bool {
    MyVector<T> Gv = G * v;
    if (!IsInteger(v.dot(Gv) / T(2))) {
      return false;
    }
    MyVector<T> w = v;
    for (int k = 1; k < phi; k++) {
      w = P.transpose() * w;
      if (!IsInteger(w.dot(Gv))) {
        return false;
      }
    }
    return true;
  };
  for (auto const &cut : parts) {
    if (over) {
      // B = {x in L^# : x U in L, x cut in L}, in the coordinates of L
      std::vector<MyMatrix<T>> As{pc.U};
      if (cut.rows() > 0) {
        As.push_back(cut);
      }
      MyMatrix<T> B = HermIntegralPreimage<T, Tint>(Ginv, As);
      // B / L is a vector space over O_E / P; when L + O_E v = B for one v
      // not in L it has dimension one, and every such v gives B
      bool cyclic = false;
      for (int i = 0; i < B.rows() && !cyclic; i++) {
        MyVector<T> v = GetMatrixRow(B, i);
        if (IsIntegralVector(v)) {
          continue;
        }
        MyMatrix<T> W = GetZbasis(orbit_rows(Id, v));
        if (T_abs(DeterminantMat(W)) != T_abs(DeterminantMat(B))) {
          break;
        }
        cyclic = true;
        if (over_even(v) && keep_if_even(W)) {
          return res;
        }
      }
      if (!cyclic) {
        bool stop = false;
        HermQuotientRepsForEach<T, Tint>(
            B, Id, [&](MyVector<T> const &v) -> bool {
              if (v.isZero() || !over_even(v)) {
                return false;
              }
              stop = keep_if_even(GetZbasis(orbit_rows(Id, v)));
              return stop;
            });
        if (stop) {
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
      bool stop = false;
      HermQuotientRepsForEach<T, Tint>(
          B, Ginv, [&](MyVector<T> const &w) -> bool {
            if (w.isZero()) {
              return false;
            }
            MyMatrix<T> Wd = GetZbasis(orbit_rows(Ginv, w));
            // the dual of the lattice of basis Wd: (Wd G Wd^T)^-1 Wd
            stop = keep_if_even(GetZbasis(MyMatrix<T>(
                Inverse(MyMatrix<T>(Wd * G * Wd.transpose())) * Wd)));
            return stop;
          });
      if (stop) {
        return res;
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
  The pieces: the pairs (L, sigma) of rank b over O_E, even, with det L
  dividing N, in the hermitian spaces of the standard pairs of the starts a
  (forms diag(1, ..., 1, a)), one per class (the canonical form of the
  pair): the search walks classes, as the neighbors of isometric pairs are
  isometric, so that no genus is missed, whether or not the invariant
  (HermitianPairInvariant) separates the genera; it is only counted for
  the report (several classes may share a genus, and so an invariant).
  The starts must meet every hermitian space that contains a piece
  (HermHermiteStarts).
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
template <typename T, typename Tint, typename Tgroup>
std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>>
HermitianPiecesSearch(
    int64_t d, int b, Tint const &N,
    std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>> const &starts,
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
  // the classes, by the canonical form of the pair (the Gram matrix and
  // the form of sigma)
  LatticeMethodChoice choice = DefaultPairMethods().canonical;
  auto class_key = [&](Tpair const &pr) -> std::vector<MyMatrix<T>> {
    return PairCanonicalForms<T, Tint, Tgroup>(
        pr.first, std::vector<MyMatrix<T>>{pr.second}, choice, os);
  };
  std::map<std::vector<MyMatrix<T>>, Tpair> visited;
  std::deque<Tpair> queue;
  size_t n_dropped = 0;
  // the invariants are taken after the descent, where the determinants
  // are small (those of the starts are large, and so are the powers of the
  // cuts)
  for (auto const &start : starts) {
    std::optional<Tpair> opt = descend(start);
    if (!opt) {
      n_dropped++;
      continue;
    }
    if (visited.emplace(class_key(*opt), *opt).second) {
      queue.push_back(*opt);
    }
  }
  os << "HERM_PIECES: d=" << d << " b=" << b << " N=" << N << ": "
     << starts.size() << " starts, " << n_dropped << " dropped, "
     << visited.size() << " classes after the descent\n";
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
        if (visited.emplace(class_key(nb), nb).second) {
          queue.push_back(nb);
        }
      }
    }
  }
  std::vector<Tpair> res;
  std::set<std::string> invariants;
  for (auto const &[key, pr] : visited) {
    res.push_back(pr);
    invariants.insert(HermitianPairInvariant<T, Tint>(pr.first, pr.second, d));
  }
  os << "HERM_PIECES: " << res.size() << " classes, " << invariants.size()
     << " invariants\n";
  return res;
}

/*
  Starts meeting every hermitian space of a piece. Let X be a piece: rank
  b over O_E, trace form q even of determinant D <= N, rank R = b phi
  over Z, with successive minima lambda_1 <= ... <= lambda_R (values of q).
  Among the first (i - 1) phi + 1 of the vectors realising them one lies
  outside the E-span of i - 1 others, so X contains E-independent y_1,
  ..., y_b with q(y_i) <= lambda_{(i - 1) phi + 1}, and by Minkowski's
  second theorem prod_i q(y_i)^phi <= prod_j lambda_j <= gamma_R^R D. Let
  H = (h(y_i, y_j)) and delta = det H, totally positive in K_0. At each
  embedding sigma of E, Hadamard gives sigma(delta) <= prod_i sigma(h(y_i,
  y_i)), so Tr_{E/Q}(delta) <= prod_i Tr_{E/Q}(h(y_i, y_i)) = prod_i q(y_i)
  <= (gamma_R^R D)^{1/phi}. As q is integral, h(X, X) is in the inverse
  different D_E^{-1}, so delta is in D_E^{-b} cap K_0. By Landherr's
  theorem (a definite hermitian space over E / K_0 is determined by its
  rank and its determinant modulo the norms) the space of X is that of
  diag(1, ..., 1, delta). So the starts: the delta of D_E^{-b} cap K_0,
  totally positive, of trace at most (gamma_R^R N)^{1/phi}; the twist of
  delta s^2 (s a positive integer, so the same class, making it integral
  and even) raised to a maximal even lattice of its space; one per class
  of these, summed with b - 1 copies of O_E (h = 1). From any lattice of
  a space the search reaches every genus of the space with determinant
  dividing N. No unit or class group of K_0 enters.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>>
HermHermiteStarts(int64_t d, int b, Tint const &N, std::ostream &os) {
  using Tpair = std::pair<MyMatrix<T>, MyMatrix<T>>;
  int64_t dn = NormalizeCyclotomicIndex(d);
  int phi = static_cast<int>(EulerPhi(dn));
  Tpair unit = HermStandardPair<T>(d, 1, std::vector<int64_t>{1});
  MyMatrix<T> const &A = unit.first;
  MyMatrix<T> const &Z = unit.second;
  // E in the basis zeta^k (row coordinates): the multiplication by zeta,
  // the traces Tr_{E/Q}(zeta^j)
  std::vector<int64_t> cyc = CyclotomicPolynomial(static_cast<int>(dn));
  MyMatrix<T> Mz = ZeroMatrix<T>(phi, phi);
  for (int k = 0; k + 1 < phi; k++) {
    Mz(k, k + 1) = 1;
  }
  for (int k = 0; k < phi; k++) {
    Mz(phi - 1, k) = T(static_cast<long>(-cyc[k]));
  }
  std::vector<T> trace_table(dn);
  for (int64_t j = 0; j < dn; j++) {
    trace_table[j] = T(static_cast<long>(RamanujanSum(dn, j)));
  }
  auto trace_zeta = [&](int64_t j) -> T {
    return trace_table[((j % dn) + dn) % dn];
  };
  // the real embeddings of K_0: zeta -> exp(2 i pi k / dn), k prime to dn
  std::vector<std::vector<double>> cosines;
  for (int64_t k = 1; 2 * k <= dn; k++) {
    if (std::gcd(k, dn) == 1) {
      std::vector<double> row(phi);
      for (int j = 0; j < phi; j++) {
        row[j] = std::cos(2 * M_PI * static_cast<double>(k * j) /
                          static_cast<double>(dn));
      }
      cosines.push_back(row);
    }
  }
  // the multiplication by x = sum_k c_k zeta^k, on the row coordinates
  auto mult_matrix = [&](MyVector<T> const &c) -> MyMatrix<T> {
    MyMatrix<T> M = ZeroMatrix<T>(phi, phi);
    MyMatrix<T> P = IdentityMat<T>(phi);
    for (int k = 0; k < phi; k++) {
      M += c(k) * P;
      P = P * Mz;
    }
    return M;
  };
  // D_E^{-1} = {x : Tr(x zeta^l) integral}, then D_E^{-b} as products
  MyMatrix<T> Tm(phi, phi);
  for (int k = 0; k < phi; k++) {
    for (int l = 0; l < phi; l++) {
      Tm(k, l) = trace_zeta(k + l);
    }
  }
  MyMatrix<T> Dinv = Inverse(Tm);
  MyMatrix<T> I = Dinv;
  for (int e = 1; e < b; e++) {
    MyMatrix<T> Gen(I.rows() * phi, phi);
    int r = 0;
    for (int u = 0; u < I.rows(); u++) {
      for (int v = 0; v < phi; v++) {
        MyVector<T> y = GetMatrixRow(Dinv, v);
        MyVector<T> x = GetMatrixRow(I, u);
        MyVector<T> xy = mult_matrix(y).transpose() * x;
        for (int w = 0; w < phi; w++) {
          Gen(r, w) = xy(w);
        }
        r++;
      }
    }
    I = GetZbasis(Gen);
  }
  // the intersection with K_0: x fixed by zeta -> zeta^-1
  MyMatrix<T> Cj(phi, phi);
  {
    MyMatrix<T> Mzinv = Inverse(Mz);
    MyMatrix<T> P = IdentityMat<T>(phi);
    for (int k = 0; k < phi; k++) {
      for (int w = 0; w < phi; w++) {
        Cj(k, w) = P(0, w);
      }
      P = P * Mzinv;
    }
  }
  MyMatrix<T> Kb;
  {
    MyMatrix<T> Eq = I * (Cj - IdentityMat<T>(phi));
    MyMatrix<Tint> Eq_i =
        UniversalMatrixConversion<Tint, T>(RemoveFractionMatrix(Eq));
    MyMatrix<T> U = UniversalMatrixConversion<T, Tint>(NullspaceIntMat(Eq_i));
    Kb = U * I;
  }
  int m = Kb.rows();
  // the bound
  int R = b * phi;
  double bound_tr = std::exp(
      (std::log(HermiteConstantPower(R)) +
       std::log(UniversalScalarConversion<double, Tint>(N))) /
      phi);
  T Btr = UniversalScalarConversion<T, int64_t>(
      static_cast<int64_t>(std::ceil(bound_tr * (1 + 1e-9))) + 1);
  T bound2 = Btr * Btr;
  // Tr_{E/Q} on the coordinates, and the form Tr(x y) on Kb
  MyVector<T> tvec(phi);
  for (int k = 0; k < phi; k++) {
    tvec(k) = trace_zeta(k);
  }
  auto trace_of = [&](MyVector<T> const &c) -> T { return c.dot(tvec); };
  auto form_of = [&](MyMatrix<T> const &basis) -> MyMatrix<T> {
    int k = basis.rows();
    MyMatrix<T> F(k, k);
    for (int u = 0; u < k; u++) {
      MyMatrix<T> Mu = mult_matrix(GetMatrixRow(basis, u));
      for (int v = 0; v < k; v++) {
        MyVector<T> prod = Mu.transpose() * GetMatrixRow(basis, v);
        F(u, v) = trace_of(prod);
      }
    }
    return F;
  };
  {
    LLLreduction<T, Tint> red = LLLreducedBasis<T, Tint>(form_of(Kb), os);
    Kb = UniversalMatrixConversion<T, Tint>(red.Pmat) * Kb;
  }
  MyMatrix<T> Q = form_of(Kb);
  // Fincke-Pohst for Tr(delta^2) <= Btr^2 (the conjugates of delta are
  // positive): Q(v) = sum_i qd_i (v_i + sum_{j > i} mu_ij v_j)^2, in doubles
  // with a margin for the pruning only, the conditions checked exactly on
  // the leaves
  std::vector<double> qd(m);
  std::vector<std::vector<double>> mu(m, std::vector<double>(m, 0));
  {
    std::vector<std::vector<double>> a(m, std::vector<double>(m));
    for (int u = 0; u < m; u++) {
      for (int v = 0; v < m; v++) {
        a[u][v] = UniversalScalarConversion<double, T>(Q(u, v));
      }
    }
    for (int u = 0; u < m; u++) {
      qd[u] = a[u][u];
      for (int v = u + 1; v < m; v++) {
        mu[u][v] = a[u][v] / qd[u];
      }
      for (int v = u + 1; v < m; v++) {
        for (int w = u + 1; w < m; w++) {
          a[v][w] -= mu[u][v] * qd[u] * mu[u][w];
        }
      }
    }
  }
  double bound2_d = UniversalScalarConversion<double, T>(bound2);
  double margin = 1e-6 * (1 + bound2_d);
  // the twist of delta: Tr(delta zeta^(k - l)), scaled by s^2 to be
  // integral and even
  auto twist_of = [&](MyVector<T> const &c) -> MyMatrix<T> {
    MyMatrix<T> G(phi, phi);
    for (int k = 0; k < phi; k++) {
      for (int l = 0; l < phi; l++) {
        T x(0);
        for (int j = 0; j < phi; j++) {
          x += c(j) * trace_zeta(j + k - l);
        }
        G(k, l) = x;
      }
    }
    for (int64_t sv = 1;; sv++) {
      T s2 = UniversalScalarConversion<T, int64_t>(sv * sv);
      MyMatrix<T> Gs = s2 * G;
      bool ok = IsIntegralMatrix(Gs);
      for (int k = 0; ok && k < phi; k++) {
        ok = IsInteger(Gs(k, k) / T(2));
      }
      if (ok) {
        return Gs;
      }
    }
  };
  auto det_of = [](MyMatrix<T> const &G) -> Tint {
    return UniversalScalarConversion<Tint, T>(T_abs(DeterminantMat(G)));
  };
  // up to a maximal even lattice, by steps up
  auto ascend = [&](Tpair pr) -> Tpair {
    while (true) {
      Tint det = det_of(pr.first);
      bool moved = false;
      if (det > 1) {
        for (auto const &[l_i, mult] : FactorsIntMap(det)) {
          int64_t l = UniversalScalarConversion<int64_t, Tint>(l_i);
          std::vector<Tpair> nbs = HermNeighbors<T, Tint>(
              pr.first, pr.second, d, l, true, false, true);
          if (!nbs.empty()) {
            pr = nbs[0];
            moved = true;
            break;
          }
        }
      }
      if (!moved) {
        return pr;
      }
    }
  };
  LatticeMethodChoice choice = DefaultPairMethods().canonical;
  auto class_key = [&](Tpair const &pr) -> std::vector<MyMatrix<T>> {
    return PairCanonicalForms<T, Tint, Tgroup>(
        pr.first, std::vector<MyMatrix<T>>{pr.second}, choice, os);
  };
  std::map<std::vector<MyMatrix<T>>, Tpair> rank_one;
  std::set<std::vector<MyMatrix<T>>> seen_twists;
  size_t n_delta = 0;
  std::vector<int64_t> v(m);
  std::function<void(int, double)> rec = [&](int i, double rest) -> void {
    if (i < 0) {
      MyVector<T> c = ZeroVector<T>(phi);
      for (int u = 0; u < m; u++) {
        if (v[u] != 0) {
          c += UniversalScalarConversion<T, int64_t>(v[u]) * GetMatrixRow(Kb, u);
        }
      }
      T tr = trace_of(c);
      if (tr <= 0 || tr > Btr) {
        return;
      }
      // totally positive: screened in doubles, then exactly
      std::vector<double> cd(phi);
      for (int j = 0; j < phi; j++) {
        cd[j] = UniversalScalarConversion<double, T>(c(j));
      }
      for (auto const &row : cosines) {
        double x = 0;
        for (int j = 0; j < phi; j++) {
          x += cd[j] * row[j];
        }
        if (x < -1e-6) {
          return;
        }
      }
      MyMatrix<T> G = twist_of(c);
      if (!IsPositiveDefinite(G, std::cerr)) {
        return;
      }
      n_delta++;
      // the twists of delta and delta u u-bar, u a unit, are isometric
      if (seen_twists.insert(class_key(Tpair{G, Z})).second) {
        Tpair top = ascend(Tpair{G, Z});
        rank_one.emplace(class_key(top), top);
      }
      return;
    }
    double cc = 0;
    for (int w = i + 1; w < m; w++) {
      cc += mu[i][w] * static_cast<double>(v[w]);
    }
    double r = std::sqrt(std::max(0.0, (rest + margin) / qd[i]));
    int64_t lo = static_cast<int64_t>(std::ceil(-cc - r - 1e-9));
    int64_t hi = static_cast<int64_t>(std::floor(-cc + r + 1e-9));
    for (int64_t x = lo; x <= hi; x++) {
      double y = static_cast<double>(x) + cc;
      double used = qd[i] * y * y;
      if (used > rest + margin) {
        continue;
      }
      v[i] = x;
      rec(i - 1, rest - used);
    }
  };
  rec(m - 1, bound2_d);
  // the starts: b - 1 copies of O_E, and the rank one lattice
  std::vector<Tpair> starts;
  int n = b * phi;
  for (auto const &[key, top] : rank_one) {
    MyMatrix<T> G = ZeroMatrix<T>(n, n);
    MyMatrix<T> P = ZeroMatrix<T>(n, n);
    for (int u = 0; u < b; u++) {
      MyMatrix<T> const &Gu = (u + 1 < b) ? A : top.first;
      MyMatrix<T> const &Pu = (u + 1 < b) ? Z : top.second;
      for (int r = 0; r < phi; r++) {
        for (int c = 0; c < phi; c++) {
          G(u * phi + r, u * phi + c) = Gu(r, c);
          P(u * phi + r, u * phi + c) = Pu(r, c);
        }
      }
    }
    starts.push_back({G, P});
  }
  os << "HERM_PIECES: Hermite starts, trace bound " << Btr << ", " << n_delta
     << " determinants, " << seen_twists.size() << " classes of twists, "
     << starts.size() << " rank one maximal classes\n";
  return starts;
}

/*
  The pieces of rank b over Z[zeta_d] with even trace lattice of
  determinant dividing N: the search from the Hermite starts, complete.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>>
HermitianPieces(int64_t d, int b, Tint const &N, std::ostream &os) {
  return HermitianPiecesSearch<T, Tint, Tgroup>(
      d, b, N, HermHermiteStarts<T, Tint, Tgroup>(d, b, N, os), os);
}

// clang-format off
#endif  // SRC_GENERA_HERMITIANLATTICE_H_
// clang-format on

// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_HERMITIANGENUSSYMBOL_H_
#define SRC_GENERA_HERMITIANGENUSSYMBOL_H_

// clang-format off
#include "HermitianNeighbor.h"
#include "HermitianGenus.h"
#include "LocalJordan.h"
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_HERMITIAN_GENUS_SYMBOL
#endif

/*
  The genus symbols of Hermitian lattices over the order O_E of a CM field
  E (CMOrder.h), relative to its real subfield K_0, as genus(L) of Hecke.jl
  (QuadForm/Herm/Genus.jl, after Kirschmer).

  The local symbol at a prime p of K_0 (HermLocalSymbol of HermitianGenus.h)
  is the list of the Jordan blocks (s, r, d) of L_p: s the valuation of the
  scale at the prime P of E above p, r the rank, d = 1 or -1 as the
  determinant of the block is a local norm from E_P or not. The global
  symbol lists the local symbols at the primes of K_0 where L is not
  unimodular, the primes ramified in E / K_0 and the dyadic primes (the
  bad primes of Hecke.jl, with discriminant and dyadic), the signatures
  being 0 for the positive definite lattices treated here.

  The Jordan blocks come from LocalJordanSplit (LocalJordan.h) on the
  Hermitian Gram matrix of a local basis at ell: m vectors of L whose
  O_E-span has an index prime to ell in L. The norm test of a in K_0 at p,
  at a prime that is not dyadic and ramified:
  - p split in E: every element is a norm;
  - p inert: the norms are the elements of even valuation;
  - p ramified: a = N(Pi)^k u, Pi a uniformizer of P, k = v_p(a), and a is
    a norm iff the unit u is a square in the residue field (the norms of
    the units are the squares of the residues, the residue extension being
    trivial).
  At a dyadic prime ramified in E / K_0 a block also has the p-valuation n
  of its norm (the ideal of the h(x, x)), written (s,r,+,n), and the norm
  test is HermDyadicNormClass.

  The notation of a prime of K_0: ell when it is the only one above ell,
  otherwise (ell,[c_0,...,c_{d-1}]) with an element sum_k c_k theta^k of
  O_{K_0} that is in it and in no other prime of K_0 above ell. A global
  symbol is written

    m=<rank>; <prime>: (s,r,+)(s,r,-)...; <prime>: ...

  the blocks as Hecke.jl prints them.
 */

// The primes of K_0 above ell: the orbits of the conjugation on the primes
// of E, each given by its first prime.
template <typename Tint> struct HermK0Prime {
  size_t i_P;
  HermIdealType type;
  std::string name;
};

// An element of O_{K_0} in the prime P of E and in no prime of E of l_other.
template <typename T, typename Tint>
MyVector<Tint> HermK0ElementIn(CMOrder<T, Tint> const &ord,
                               HermIdealPrime<Tint> const &P,
                               std::vector<HermIdealPrime<Tint>> const &l_other) {
  int d = ord.deg;
  MyMatrix<Tint> Pc = P.ideal;
  MyMatrix<Tint> Diff(Pc.rows(), d);
  for (int r = 0; r < Pc.rows(); r++) {
    MyVector<Tint> v = GetMatrixRow(Pc, r);
    MyVector<Tint> w = ord.conjugate(v) - v;
    for (int a = 0; a < d; a++) {
      Diff(r, a) = w(a);
    }
  }
  MyMatrix<Tint> Fixed = NullspaceIntMat(Diff) * Pc;
  for (int iter = 0; iter < 10000; iter++) {
    int range = 1 + iter / 1000;
    MyVector<Tint> a = ZeroVector<Tint>(d);
    for (int i = 0; i < Fixed.rows(); i++) {
      int c = random_int(-range, range);
      if (c != 0) {
        a += Tint(c) * GetMatrixRow(Fixed, i);
      }
    }
    bool ok = !IsZeroVector(a);
    for (auto const &Q : l_other) {
      if (ok && HermResidue(Q, a).is_zero()) {
        ok = false;
      }
    }
    if (ok) {
      return a;
    }
  }
  std::cerr << "HERMSYM: no element of O_K0 separating the primes\n";
  throw TerminalException{1};
}

template <typename T, typename Tint>
std::string HermThetaString(CMOrder<T, Tint> const &ord, MyVector<Tint> const &a) {
  MyVector<T> a_T = UniversalVectorConversion<T, Tint>(a);
  MyVector<T> th = ord.basis_theta.transpose() * a_T;
  std::ostringstream os;
  os << "[";
  for (int k = 0; k < th.size(); k++) {
    os << (k > 0 ? "," : "") << th(k);
  }
  os << "]";
  return os.str();
}

template <typename T, typename Tint>
std::vector<HermK0Prime<Tint>>
HermK0PrimesAbove(CMOrder<T, Tint> const &ord,
                  std::vector<HermIdealPrime<Tint>> const &l_prime) {
  std::vector<HermK0Prime<Tint>> res;
  for (size_t i = 0; i < l_prime.size(); i++) {
    if (l_prime[i].conj_index < i) {
      continue;
    }
    res.push_back({i, l_prime[i].type, ""});
  }
  int64_t ell = l_prime[0].ell;
  for (auto &p : res) {
    if (res.size() == 1) {
      p.name = std::to_string(ell);
      continue;
    }
    std::vector<HermIdealPrime<Tint>> l_other;
    for (size_t j = 0; j < l_prime.size(); j++) {
      if (j != p.i_P && j != l_prime[p.i_P].conj_index) {
        l_other.push_back(l_prime[j]);
      }
    }
    MyVector<Tint> a = HermK0ElementIn(ord, l_prime[p.i_P], l_other);
    p.name = "(" + std::to_string(ell) + "," + HermThetaString(ord, a) + ")";
  }
  return res;
}

// The valuation v_P of the different of E / K_0 at a ramified prime P:
// O_{E,P} = O_{K_0,p}[Pi] for a uniformizer Pi (Eisenstein), so the
// different is generated by Pi - conj(Pi).
template <typename T, typename Tint>
int HermDifferentVal(CMOrder<T, Tint> const &ord,
                     std::vector<HermIdealPrime<Tint>> const &l_prime,
                     size_t i_P) {
  HermIdealPrime<Tint> const &P = l_prime[i_P];
  std::vector<HermIdealPrime<Tint>> l_other;
  for (size_t j = 0; j < l_prime.size(); j++) {
    if (j != i_P) {
      l_other.push_back(l_prime[j]);
    }
  }
  MyVector<Tint> Pi = HermFindUniformizer<T, Tint>(P.ideal, P, l_other);
  MyVector<Tint> diff = Pi - ord.conjugate(Pi);
  return HermIdealVal(ord, P, UniversalVectorConversion<T, Tint>(diff));
}

/*
  The norm test at a dyadic prime ramified in E / K_0, e_0 = v_p(2): 1 +
  4 pi O_{K_0,p} consists of squares, so a is a norm iff N(y) / a is in it
  for some y = Pi^v w, v = v_P(a) / 2 and w a unit; and since Tr(P^M) is in
  p^{floor(M / 2)}, w matters modulo P^M, M = 2 (2 e_0 + 1), only. The
  units w are taken in the box of the coordinates in [0, ell^k), ell^k O_E
  in P^M.
 */
template <typename T, typename Tint>
int HermDyadicNormClass(CMOrder<T, Tint> const &ord,
                        std::vector<HermIdealPrime<Tint>> const &l_prime,
                        size_t i_P, MyVector<T> const &a) {
  HermIdealPrime<Tint> const &P = l_prime[i_P];
  int d = ord.deg;
  int e0 = P.e / 2;
  int M = 2 * (2 * e0 + 1);
  int k = (M + P.e - 1) / P.e;
  int64_t box = 1;
  for (int i = 0; i < k; i++) {
    box *= P.ell;
  }
  std::vector<HermIdealPrime<Tint>> l_other;
  for (size_t j = 0; j < l_prime.size(); j++) {
    if (j != i_P) {
      l_other.push_back(l_prime[j]);
    }
  }
  // the multiplication matrix of an element of E (rational coordinates)
  auto mult = [&](MyVector<T> const &x) -> MyMatrix<T> {
    MyMatrix<T> R = ZeroMatrix<T>(d, d);
    for (int c = 0; c < d; c++) {
      if (x(c) != 0) {
        R += x(c) * UniversalMatrixConversion<T, Tint>(ord.mult[c]);
      }
    }
    return R;
  };
  MyMatrix<T> Cj = UniversalMatrixConversion<T, Tint>(ord.conj).transpose();
  MyVector<T> Pi = UniversalVectorConversion<T, Tint>(
      HermFindUniformizer<T, Tint>(P.ideal, P, l_other));
  int v = HermIdealVal(ord, P, a) / 2;
  // Pi^v
  MyVector<T> Pv = ZeroVector<T>(d);
  Pv(0) = 1;
  MyMatrix<T> M_Pi = mult(Pi).transpose();
  MyMatrix<T> M_Pi_inv = Inverse(M_Pi);
  for (int i = 0; i < std::abs(v); i++) {
    Pv = (v > 0) ? MyVector<T>(M_Pi * Pv) : MyVector<T>(M_Pi_inv * Pv);
  }
  MyMatrix<T> M_a_inv = Inverse(MyMatrix<T>(mult(a).transpose()));
  MyMatrix<T> M_Pv = mult(Pv).transpose();
  MyVector<T> one = ZeroVector<T>(d);
  one(0) = 1;
  int64_t total = 1;
  for (int c = 0; c < d; c++) {
    total *= box;
  }
  for (int64_t idx = 0; idx < total; idx++) {
    MyVector<Tint> w(d);
    int64_t t = idx;
    for (int c = 0; c < d; c++) {
      w(c) = UniversalScalarConversion<Tint, int64_t>(t % box);
      t /= box;
    }
    if (HermResidue(P, w).is_zero()) {
      continue;
    }
    MyVector<T> y = M_Pv * UniversalVectorConversion<T, Tint>(w);
    MyVector<T> Ny = mult(y).transpose() * (Cj * y);
    MyVector<T> x = M_a_inv * Ny;
    if (HermIdealVal(ord, P, MyVector<T>(x - one)) >= M) {
      return 1;
    }
  }
  return -1;
}

/*
  The residue at P = l_prime[i_P] of an element u of E with v_P(u) >= 0,
  whatever its valuations at the other primes above ell: u is multiplied by
  y, a product of powers of uniformizers of those primes taken outside P,
  so that y u is integral at every prime above ell (its coordinates then
  have no ell in their denominators), and res(u) = res(y u) / res(y).
 */
template <typename T, typename Tint>
FqElt HermResidueAtP(CMOrder<T, Tint> const &ord,
                     std::vector<HermIdealPrime<Tint>> const &l_prime,
                     size_t i_P, MyVector<T> const &u) {
  HermIdealPrime<Tint> const &P = l_prime[i_P];
  int d = ord.deg;
  MyVector<T> yu = u;
  MyVector<Tint> y = ZeroVector<Tint>(d);
  y(0) = 1;
  for (size_t j = 0; j < l_prime.size(); j++) {
    if (j == i_P) {
      continue;
    }
    int v = HermIdealVal(ord, l_prime[j], yu);
    if (v >= 0) {
      continue;
    }
    std::vector<HermIdealPrime<Tint>> l_other;
    for (size_t k = 0; k < l_prime.size(); k++) {
      if (k != j) {
        l_other.push_back(l_prime[k]);
      }
    }
    MyVector<Tint> pi =
        HermFindUniformizer<T, Tint>(l_prime[j].ideal, l_prime[j], l_other);
    MyMatrix<T> Mpi =
        UniversalMatrixConversion<T, Tint>(ord.mult_matrix(pi)).transpose();
    for (int t = 0; t < -v; t++) {
      yu = Mpi * yu;
      y = ord.product(y, pi);
    }
  }
  FqElt num = HermResidueRat(P, yu);
  FqElt den = HermResidue(P, y);
  return FqMul(P.F, num, FqInv(P.F, den));
}

/*
  Whether a, an element of K_0 (coordinates on e), is a local norm from E
  at the prime of K_0 under P: 1 or -1.
 */
template <typename T, typename Tint>
int HermLocalNormClass(CMOrder<T, Tint> const &ord,
                       std::vector<HermIdealPrime<Tint>> const &l_prime,
                       size_t i_P, MyVector<T> const &a) {
  HermIdealPrime<Tint> const &P = l_prime[i_P];
  if (P.type == HermIdealType::Split) {
    return 1;
  }
  int v = HermIdealVal(ord, P, a);
  if (P.type == HermIdealType::Inert) {
    return (v % 2 == 0) ? 1 : -1;
  }
  if (P.ell == 2) {
    return HermDyadicNormClass(ord, l_prime, i_P, a);
  }
  std::vector<HermIdealPrime<Tint>> l_other;
  for (size_t j = 0; j < l_prime.size(); j++) {
    if (j != i_P) {
      l_other.push_back(l_prime[j]);
    }
  }
  MyVector<Tint> Pi = HermFindUniformizer<T, Tint>(P.ideal, P, l_other);
  MyVector<Tint> NPi = ord.product(Pi, ord.conjugate(Pi));
  MyMatrix<T> M = UniversalMatrixConversion<T, Tint>(ord.mult_matrix(NPi));
  MyMatrix<T> Minv = Inverse(M).transpose();
  MyMatrix<T> Mt = M.transpose();
  // v is even (a is in K_0); a negative valuation is that of a lattice
  // of negative scale
  MyVector<T> u = a;
  for (int k = 0; k < v / 2; k++) {
    u = Minv * u;
  }
  for (int k = 0; k < -v / 2; k++) {
    u = Mt * u;
  }
  FqElt r = HermResidueAtP(ord, l_prime, i_P, u);
  if (r.is_zero()) {
    std::cerr << "HERMSYM: the unit part has a zero residue at a prime above "
              << P.ell << " (e=" << P.e << ", f=" << P.f << "), v_P=" << v
              << "\n";
    throw TerminalException{1};
  }
  return FqQuadraticCharacter(P.F, r);
}

/*
  A local basis of L at ell: m vectors whose O_E-span has index prime to
  ell, by random small combinations of the basis of L.
 */
template <typename T, typename Tint>
std::vector<MyVector<Tint>> HermLocalBasis(HermLatticeData<T, Tint> const &L,
                                           int64_t ell) {
  int n = L.n();
  int d = L.ord->deg;
  int m = n / d;
  Tint ell_T = UniversalScalarConversion<Tint, int64_t>(ell);
  for (int iter = 0; iter < 10000; iter++) {
    int range = 1 + iter / 1000;
    std::vector<MyVector<Tint>> l_v;
    MyMatrix<Tint> Gen(m * d, n);
    for (int i = 0; i < m; i++) {
      MyVector<Tint> v = ZeroVector<Tint>(n);
      for (int k = 0; k < n; k++) {
        v(k) = Tint(random_int(-range, range));
      }
      l_v.push_back(v);
      for (int a = 0; a < d; a++) {
        MyVector<Tint> w = L.Omega[a].transpose() * v;
        for (int k = 0; k < n; k++) {
          Gen(i * d + a, k) = w(k);
        }
      }
    }
    Tint det = DeterminantMat(Gen);
    if (det == 0 || det % ell_T == 0) {
      continue;
    }
    return l_v;
  }
  std::cerr << "HERMSYM: no local basis found at " << ell << "\n";
  throw TerminalException{1};
}

/*
  The local symbol of L at the prime of K_0 under l_prime[i_P].
 */
template <typename T, typename Tint>
HermLocalSymbol HermLocalSymbolOfLattice(
    HermLatticeData<T, Tint> const &L,
    std::vector<HermIdealPrime<Tint>> const &l_prime, size_t i_P) {
  CMOrder<T, Tint> const &ord = *L.ord;
  NumberField<T> const &K = ord.field;
  HermIdealPrime<Tint> const &P = l_prime[i_P];
  bool dyadic_ram = P.type == HermIdealType::Ramified && P.ell == 2;
  std::vector<MyVector<Tint>> l_v = HermLocalBasis(L, P.ell);
  int m = static_cast<int>(l_v.size());
  MyMatrix<T> Bt = ord.basis_theta;
  MyMatrix<T> Bt_inv = Inverse(Bt);
  auto to_theta = [&](MyVector<T> const &a_e) -> ElemE<T> {
    MyVector<T> th = Bt.transpose() * a_e;
    ElemE<T> x(K.deg, T(0));
    for (int k = 0; k < K.deg; k++) {
      x[k] = th(k);
    }
    return x;
  };
  auto to_e = [&](ElemE<T> const &x) -> MyVector<T> {
    MyVector<T> th(K.deg);
    for (int k = 0; k < K.deg; k++) {
      th(k) = x[k];
    }
    return Bt_inv.transpose() * th;
  };
  SmallGram<T> H(m, std::vector<ElemE<T>>(m));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      H[i][j] = to_theta(L.herm(l_v[i], l_v[j]));
    }
  }
  auto val = [&](ElemE<T> const &x) -> int {
    return x.is_zero() ? (1 << 28) : HermIdealVal(ord, P, to_e(x));
  };
  if (P.type == HermIdealType::Split) {
    /*
      At a split prime O_E (x) O_{K_0,p} = O_P x O_conj(P) is not a local
      ring: L_p is a pair of lattices paired by h, and the symbol is given
      by the elementary divisors of that pairing at P, the blocks (s, r, +)
      (every element is a norm). They are those of H over the localization
      of O_E at P, a discrete valuation ring: a Smith elimination with
      independent row and column operations, the pivot of least v_P.
     */
    SmallGram<T> M = H;
    std::map<int, int> count;
    int sz = m;
    while (sz > 0) {
      int bi = 0, bj = 0, bv = 1 << 28;
      for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz; j++) {
          int v = val(M[i][j]);
          if (v < bv) {
            bv = v;
            bi = i;
            bj = j;
          }
        }
      }
      if (bv == (1 << 28)) {
        std::cerr << "HERMSYM: the form is degenerate at a split prime\n";
        throw TerminalException{1};
      }
      count[bv]++;
      ElemE<T> inv = EInv(M[bi][bj], K);
      // clear the column bj and the row bi, then remove them
      for (int i = 0; i < sz; i++) {
        if (i == bi || M[i][bj].is_zero()) {
          continue;
        }
        ElemE<T> c = EMul(M[i][bj], inv, K);
        for (int j = 0; j < sz; j++) {
          M[i][j] = M[i][j] - EMul(c, M[bi][j], K);
        }
      }
      SmallGram<T> M2;
      for (int i = 0; i < sz; i++) {
        if (i == bi) {
          continue;
        }
        std::vector<ElemE<T>> row;
        for (int j = 0; j < sz; j++) {
          if (j != bj) {
            row.push_back(M[i][j]);
          }
        }
        M2.push_back(row);
      }
      M = M2;
      sz--;
    }
    HermLocalSymbol sym;
    for (auto const &[sv, r] : count) {
      sym.blocks.push_back({sv, r, 1});
    }
    return sym;
  }
  std::vector<LocalJordanComponent<T>> comps = LocalJordanSplitVal<T>(
      H, K, val, P.ell == 2, JordanForm::Hermitian, nullptr);
  HermLocalSymbol sym;
  for (auto const &C : comps) {
    MyVector<T> det_e = to_e(LocalComponentDet(C, K));
#ifdef SANITY_CHECK_HERMITIAN_GENUS_SYMBOL
    if (UniversalMatrixConversion<T, Tint>(ord.conj).transpose() * det_e !=
        det_e) {
      std::cerr << "HERMSYM: the determinant of a block is not in K_0\n";
      throw TerminalException{1};
    }
#endif
    int cls = HermLocalNormClass(ord, l_prime, i_P, det_e);
    sym.blocks.push_back({C.s, C.rank(), cls});
    if (dyadic_ram) {
      // the norm: the diagonal entries and the traces Tr_{E/K_0}(x O_E) of
      // the others (_get_normsuation_from_gram_matrix of Hecke.jl)
      int nv = 1 << 28;
      MyMatrix<T> Cj = UniversalMatrixConversion<T, Tint>(ord.conj).transpose();
      for (auto const &G : C.pieces) {
        int sz = static_cast<int>(G.size());
        for (int i = 0; i < sz; i++) {
          nv = std::min(nv, HermIdealVal(ord, P, to_e(G[i][i])));
          for (int j = i + 1; j < sz; j++) {
            MyVector<T> x = to_e(G[i][j]);
            for (int a = 0; a < K.deg; a++) {
              MyVector<Tint> ea = ZeroVector<Tint>(K.deg);
              ea(a) = 1;
              MyVector<T> xe = UniversalMatrixConversion<T, Tint>(
                                   ord.mult_matrix(ea)).transpose() * x;
              nv = std::min(nv, HermIdealVal(ord, P, MyVector<T>(xe + Cj * xe)));
            }
          }
        }
      }
      if (nv % 2 != 0) {
        std::cerr << "HERMSYM: the norm valuation is odd at a ramified prime\n";
        throw TerminalException{1};
      }
      sym.norms.push_back(nv / 2);
    }
  }
  return sym;
}

// The global symbol: the prime of K_0, its rational prime and its local
// symbol.
struct HermGlobalSymbolEntry {
  int64_t ell;
  std::string name;
  HermLocalSymbol sym;
};

struct HermGlobalSymbol {
  int m;
  std::vector<HermGlobalSymbolEntry> entries;
};

inline std::string HermLocalSymbolString(HermLocalSymbol const &sym) {
  std::ostringstream os;
  for (size_t i = 0; i < sym.blocks.size(); i++) {
    auto const &b = sym.blocks[i];
    os << "(" << b[0] << "," << b[1] << "," << (b[2] == 1 ? "+" : "-");
    if (!sym.norms.empty()) {
      os << "," << sym.norms[i];
    }
    os << ")";
  }
  return os.str();
}

inline std::string HermGlobalSymbolString(HermGlobalSymbol const &G) {
  std::ostringstream os;
  os << "m=" << G.m;
  for (auto const &e : G.entries) {
    os << "; " << e.name << ": " << HermLocalSymbolString(e.sym);
  }
  return os.str();
}

/*
  The global symbol of L: the primes of K_0 above the divisors of 2 det of
  the trace form (which contains those of the discriminant of O_E), kept
  when dyadic, ramified in E / K_0 or with a symbol that is not unimodular.
 */
template <typename T, typename Tint>
HermGlobalSymbol HermGenusSymbolOfLattice(HermLatticeData<T, Tint> const &L) {
  CMOrder<T, Tint> const &ord = *L.ord;
  int m = L.n() / ord.deg;
  Tint det = T_abs(DeterminantMat(L.G));
  std::vector<int64_t> l_ell{2};
  Tint rest = det;
  while (rest % 2 == 0) {
    rest /= 2;
  }
  for (int64_t p = 3;; p += 2) {
    Tint p_T = UniversalScalarConversion<Tint, int64_t>(p);
    if (p_T * p_T > rest) {
      break;
    }
    if (rest % p_T == 0) {
      l_ell.push_back(p);
      while (rest % p_T == 0) {
        rest /= p_T;
      }
    }
  }
  if (rest > 1) {
    l_ell.push_back(UniversalScalarConversion<int64_t, Tint>(rest));
  }
  HermGlobalSymbol G;
  G.m = m;
  HermLocalSymbol unimod;
  unimod.blocks.push_back({0, m, 1});
  for (int64_t ell : l_ell) {
    std::vector<HermIdealPrime<Tint>> l_prime = HermIdealPrimesAbove(ord, ell);
    for (auto const &p : HermK0PrimesAbove(ord, l_prime)) {
      HermLocalSymbol sym = HermLocalSymbolOfLattice(L, l_prime, p.i_P);
      if (ell == 2 || p.type == HermIdealType::Ramified || !(sym == unimod)) {
        G.entries.push_back({ell, p.name, sym});
      }
    }
  }
  return G;
}

// clang-format off
#endif  // SRC_GENERA_HERMITIANGENUSSYMBOL_H_
// clang-format on

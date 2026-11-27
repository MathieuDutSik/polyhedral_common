// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_HERMITIANENUMERATION_H_
#define SRC_GENERA_HERMITIANENUMERATION_H_

// clang-format off
#include "HermitianGenusSeed.h"
#include "MassAbelianHermitian.h"
#include "CyclotomicArithmetic.h"
#include "ScalarAbelianField.h"
#include "LatticeStabEquiCan.h"
#include "ClassicLLL.h"
#include "Namelist.h"
#include <map>
#include <queue>
#include <tuple>
#include <optional>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_HERMITIAN_ENUMERATION
#endif

#ifdef SANITY_CHECK
# define SANITY_CHECK_HERMITIAN_ENUMERATION
#endif

#ifdef TIMINGS
#define TIMINGS_HERMITIAN_ENUMERATION
#endif

/*
  The enumeration of a genus of positive definite Hermitian lattices over
  the order O_E of a CM field (CMOrder.h) by Kneser neighbors
  (HermitianNeighbor.h), from a seed (HermitianGenusSeed.h).

  A class is kept by a canonical basis: the canonical form of the pair
  (G, G Omega^T) of the trace form and the multiplication by theta (an
  E-linear map preserving Tr h preserves h), on the LLL-reduced lattice. Its
  automorphism group is that of the same pair, the unitary group U(L), whose
  order enters the mass sum 1 / |U(L)|.

  A class is expanded at a prime P of O_E where it is unimodular (the
  neighbors of HermitianNeighbor.h: split, inert and ramified, the dyadic
  ramified ones included): one neighbor computation per orbit of U(L) on the lines of
  L / C L, the automorphisms being O_E-linear so that they preserve C L.
  The neighbor primes are the primes of E above the rational primes where
  the seed is unimodular, by increasing norm N(P); all the classes are
  expanded at the first one, and while the mass is not reached the next
  one is brought in and all the classes are expanded at it. The neighbor
  graphs of several primes connect the genus when their classes generate
  the relevant quotient of the class group, which the mass certifies: the
  enumeration is complete exactly when the sum of the 1 / |U(L)| reaches
  the mass of the genus.
 */

/*
  The abelian CM field E fixed by H in Q(zeta_n) as Q(theta): theta =
  Tr_{Q(zeta_n)/E}(beta) for beta = zeta^k (then zeta + zeta^k, ...) taken
  when its conjugates, the sigma_c(theta) for the cosets c of H, are
  distinct; its minimal polynomial is prod_c (x - sigma_c(theta)), computed
  in Q(zeta_n), and its conjugate sigma_{-1}(theta) is written on the powers
  of theta. With this field the mass of the genus is exact
  (HermitianMassAbelianLattice of MassAbelianHermitian.h).
 */
template <typename T>
NumberField<T> NumberFieldOfAbelianCM(AbelianField const &F) {
  int64_t n = F.n;
  NumberField<T> Z;
  Z.d = n;
  Z.Phi = CyclotomicPolynomial(static_cast<int>(n));
  Z.deg = static_cast<int>(Z.Phi.size()) - 1;
  Z.conj = ElemE<T>(Z.deg, T(0));
  std::set<int64_t> Hs(F.H.begin(), F.H.end());
  // the coset representatives of H in (Z/n)^*
  std::vector<int64_t> reps;
  std::set<int64_t> covered;
  for (int64_t c = 1; c < std::max<int64_t>(n, 2); c++) {
    if (std::gcd(c, n) != 1 || covered.count(c % n)) {
      continue;
    }
    reps.push_back(c);
    for (int64_t h : F.H) {
      covered.insert((c * h) % n);
    }
  }
  int deg = static_cast<int>(reps.size());
  // sigma_c(Tr(beta)) for beta = sum_k b_k zeta^k
  auto conjugate = [&](std::vector<int64_t> const &beta, int64_t c) -> ElemE<T> {
    ElemE<T> x(Z.deg, T(0));
    for (int64_t h : F.H) {
      for (size_t k = 0; k < beta.size(); k++) {
        if (beta[k] != 0) {
          x = x + UniversalScalarConversion<T, int64_t>(beta[k]) * EZetaPower(static_cast<int64_t>(k) * c * h, Z);
        }
      }
    }
    return x;
  };
  for (int64_t trial = 1; trial < 4 * n + 4; trial++) {
    std::vector<int64_t> beta(n, 0);
    if (trial < n) {
      beta[trial] = 1;
    } else {
      beta[1] = 1;
      beta[(trial % (n - 1)) + 1] += 1 + trial / n;
    }
    std::vector<ElemE<T>> l_conj;
    for (int64_t c : reps) {
      l_conj.push_back(conjugate(beta, c));
    }
    bool distinct = true;
    for (int i = 0; i < deg && distinct; i++) {
      for (int j = i + 1; j < deg && distinct; j++) {
        distinct = !(l_conj[i] == l_conj[j]);
      }
    }
    if (!distinct) {
      continue;
    }
    // prod_c (x - theta_c), coefficients in Q(zeta_n)
    std::vector<ElemE<T>> poly{EOne(Z)};
    for (auto const &tc : l_conj) {
      std::vector<ElemE<T>> next(poly.size() + 1, ElemE<T>(Z.deg, T(0)));
      for (size_t k = 0; k < poly.size(); k++) {
        next[k + 1] = next[k + 1] + poly[k];
        next[k] = next[k] - EMul(poly[k], tc, Z);
      }
      poly = next;
    }
    NumberField<T> K;
    K.d = 0;
    K.deg = deg;
    for (auto const &c : poly) {
      for (int k = 1; k < Z.deg; k++) {
        if (c[k] != 0) {
          std::cerr << "HERMENUM: the minimal polynomial is not rational\n";
          throw TerminalException{1};
        }
      }
      if (!IsInteger(c[0])) {
        std::cerr << "HERMENUM: the minimal polynomial is not integral\n";
        throw TerminalException{1};
      }
      K.Phi.push_back(UniversalScalarConversion<int64_t, T>(c[0]));
    }
    // sigma_{-1}(theta) on the powers of theta
    ElemE<T> theta = l_conj[0];
    MyMatrix<T> Pow(deg, Z.deg);
    ElemE<T> pw = EOne(Z);
    for (int k = 0; k < deg; k++) {
      for (int j = 0; j < Z.deg; j++) {
        Pow(k, j) = pw[j];
      }
      pw = EMul(pw, theta, Z);
    }
    ElemE<T> cth = conjugate(beta, n - 1);
    MyVector<T> v(Z.deg);
    for (int j = 0; j < Z.deg; j++) {
      v(j) = cth[j];
    }
    std::optional<MyVector<T>> sol = SolutionMat(Pow, v);
    if (!sol) {
      std::cerr << "HERMENUM: the conjugate is not in Q(theta)\n";
      throw TerminalException{1};
    }
    K.conj = ElemE<T>(deg, T(0));
    for (int k = 0; k < deg; k++) {
      K.conj[k] = (*sol)(k);
    }
    return K;
  }
  std::cerr << "HERMENUM: no generator of the fixed field found\n";
  throw TerminalException{1};
}

/*
  The pair (G, G Omega^T) of the lattice of basis A (rows in the ambient
  space): G its trace form, Omega the multiplication by theta on the basis.
  The integral matrices preserving both are the unitary isometries, so the
  automorphism group, the canonical form and the isometries of the pair are
  those of the Hermitian lattice.
 */
template <typename T, typename Tint>
std::vector<MyMatrix<T>> HermPairOfBasis(RealifiedRing<T, Tint> const &ring,
                                         MyMatrix<T> const &A) {
  MyMatrix<T> G = A * UniversalMatrixConversion<T, Tint>(ring.B) * A.transpose();
  MyMatrix<T> Om = A * UniversalMatrixConversion<T, Tint>(ring.W) * Inverse(A);
  return {G, MyMatrix<T>(G * Om.transpose())};
}

/*
  The field of a namelist block (GENUS_HermitianEnumerate,
  GENUS_HermitianGenera): AbelianConductor and AbelianGenerators, or
  FieldPolynomial and FieldConjugate. The error message, for the program to
  report, when the input is not valid.
 */
template <typename T> struct HermFieldInput {
  NumberField<T> K;
  std::optional<AbelianField> opt_F;
};

template <typename T>
std::optional<HermFieldInput<T>> HermFieldFromBlock(SingleBlock const &B,
                                                    std::string &error) {
  HermFieldInput<T> res;
  int conductor = B.get_int("AbelianConductor");
  if (conductor > 0) {
    std::vector<int64_t> gens;
    for (int g : B.get_list_int("AbelianGenerators")) {
      gens.push_back(g);
    }
    res.opt_F = MakeAbelianCMField(conductor, gens);
    res.K = NumberFieldOfAbelianCM<T>(*res.opt_F);
    return res;
  }
  res.K.d = 0;
  for (int c : B.get_list_int("FieldPolynomial")) {
    res.K.Phi.push_back(c);
  }
  res.K.deg = static_cast<int>(res.K.Phi.size()) - 1;
  std::vector<std::string> l_conj = B.get_list_string("FieldConjugate");
  if (res.K.deg < 1 || static_cast<int>(l_conj.size()) != res.K.deg) {
    error = "the field is given by AbelianConductor or by FieldPolynomial "
            "and FieldConjugate (" + std::to_string(res.K.deg) +
            " coefficients)";
    return {};
  }
  res.K.conj = ElemE<T>(res.K.deg, T(0));
  for (int k = 0; k < res.K.deg; k++) {
    res.K.conj[k] = ParseScalar<T>(l_conj[k]);
  }
  return res;
}

template <typename T, typename Tint> struct HermEnumClass {
  // the canonical basis, rows in the coordinates of the ambient space
  MyMatrix<T> basis;
  mpz_class aut_order;
  // the generators of U(L) on the coordinates of basis
  std::vector<MyMatrix<Tint>> l_gen;
};

template <typename T, typename Tint> struct HermEnumResult {
  std::vector<HermEnumClass<T, Tint>> l_class;
  T mass;
  std::optional<T> target_mass;
  bool complete = false;
  std::vector<int64_t> primes_used;
};

/*
  The orbits of the group generated by l_gen (on the coordinates of L) on
  the lines of the line space: a representative and the size of each.
 */
template <typename Tint>
std::vector<std::pair<size_t, size_t>>
HermLineOrbits(HermLineSpace const &ls, std::vector<MyMatrix<Tint>> const &l_gen) {
  int n = ls.n;
  int64_t ell = ls.ell();
  size_t n_line = ls.n_line();
  std::vector<std::vector<int64_t>> gen_red;
  for (auto const &g : l_gen) {
    std::vector<int64_t> gr(n * n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        gr[i * n + j] = order_residue(g(i, j), ell);
      }
    }
    gen_red.push_back(gr);
  }
  std::vector<char> seen(n_line, 0);
  std::vector<std::pair<size_t, size_t>> res;
  std::vector<size_t> stack;
  std::vector<int64_t> w(n);
  for (size_t idx = 0; idx < n_line; idx++) {
    if (seen[idx]) {
      continue;
    }
    seen[idx] = 1;
    stack.push_back(idx);
    size_t size = 0;
    while (!stack.empty()) {
      size_t cur = stack.back();
      stack.pop_back();
      size++;
      std::vector<int64_t> z = ls.line_vector(cur);
      for (auto const &gr : gen_red) {
        for (int j = 0; j < n; j++) {
          int64_t s = 0;
          for (int i = 0; i < n; i++) {
            s = (s + z[i] * gr[i * n + j]) % ell;
          }
          w[j] = s;
        }
        std::optional<size_t> opt = ls.line_index(w);
        if (opt && !seen[*opt]) {
          seen[*opt] = 1;
          stack.push_back(*opt);
        }
      }
    }
    res.push_back({idx, size});
  }
  return res;
}

template <typename T, typename Tint, typename Tgroup>
HermEnumResult<T, Tint>
HermitianEnumeration(CMOrder<T, Tint> const &ord, HermSeed<T, Tint> const &seed,
                     std::optional<T> const &target_mass, int max_prime,
                     std::ostream &os) {
  RealifiedRing<T, Tint> const &ring = seed.ring;
  MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(ring.B);
  MyMatrix<T> W_T = UniversalMatrixConversion<T, Tint>(ring.W);
  HermEnumResult<T, Tint> res;
  res.mass = T(0);
  res.target_mass = target_mass;
  std::map<std::pair<MyMatrix<T>, MyMatrix<T>>, size_t> known;
  auto closed = [&]() -> bool {
    return target_mass && res.mass == *target_mass;
  };
  auto insert = [&](MyMatrix<T> const &A_in) -> void {
    LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(
        MyMatrix<T>(A_in * B_T * A_in.transpose()), os);
    MyMatrix<T> A = UniversalMatrixConversion<T, Tint>(rec.Pmat) * A_in;
    std::vector<MyMatrix<T>> l_mat = HermPairOfBasis(ring, A);
    MyMatrix<Tint> Pc =
        ComputeCanonicalFormMultiple<T, Tint, Tgroup>(l_mat, os);
    MyMatrix<T> Pc_T = UniversalMatrixConversion<T, Tint>(Pc);
    MyMatrix<T> Ac = Pc_T * A;
    std::pair<MyMatrix<T>, MyMatrix<T>> key{Ac * B_T * Ac.transpose(),
                                            Ac * W_T * Inverse(Ac)};
    if (known.count(key)) {
      return;
    }
    IntAutomorphismGroup<Tint, Tgroup> aut =
        ArithmeticAutomorphismGroupMultiple<T, Tint, Tgroup>(l_mat, os);
    HermEnumClass<T, Tint> cl;
    cl.basis = Ac;
    cl.aut_order = aut.order;
    MyMatrix<T> Pc_inv = Inverse(Pc_T);
    for (auto const &g : aut.ListGen) {
      MyMatrix<T> gc = Pc_T * UniversalMatrixConversion<T, Tint>(g) * Pc_inv;
      cl.l_gen.push_back(UniversalMatrixConversion<Tint, T>(gc));
    }
    known[key] = res.l_class.size();
    res.l_class.push_back(cl);
    res.mass += T(1) / UniversalScalarConversion<T, mpz_class>(aut.order);
#ifdef DEBUG_HERMITIAN_ENUMERATION
    os << "HERMENUM: class " << res.l_class.size() << " |U(L)|=" << aut.order
       << " mass " << res.mass << "\n";
#endif
  };
  insert(seed.basis);
  HermLatticeData<T, Tint> L0 = HermLatticeOfBasis(ord, ring, seed.basis);
  // The primes P of E above ell <= max_prime, ell where the genus is
  // unimodular, by increasing norm N(P): the cost of the neighbors at P
  // grows with N(P) (the lines of L / P L), and an inert prime of norm ell^2
  // never changes the special genus, so a split prime of small norm comes
  // before it. As N(P) >= ell, the primes above every ell <= N(P) are known
  // before P is taken.
  std::map<int64_t, std::vector<HermIdealPrime<Tint>>> primes_above;
  // (N(P), ell, index of P above ell), smallest first
  std::priority_queue<std::tuple<uint64_t, int64_t, size_t>,
                      std::vector<std::tuple<uint64_t, int64_t, size_t>>,
                      std::greater<std::tuple<uint64_t, int64_t, size_t>>>
      pending;
  int64_t next_ell = 2;
  auto decompose_up_to = [&](uint64_t bound) -> void {
    while (next_ell <= max_prime && static_cast<uint64_t>(next_ell) <= bound) {
      int64_t ell = next_ell++;
      bool is_prime = true;
      for (int64_t q = 2; q * q <= ell; q++) {
        if (ell % q == 0) {
          is_prime = false;
        }
      }
      if (!is_prime || !HermIsUnimodularAt(L0, ell)) {
        continue;
      }
      std::vector<HermIdealPrime<Tint>> l_p = HermIdealPrimesAbove(ord, ell);
#ifdef DEBUG_HERMITIAN_ENUMERATION
      for (auto const &P : l_p) {
        os << "HERMENUM: above " << ell << " the prime " << P.name
           << " of norm " << P.q() << "\n";
      }
#endif
      for (size_t ip = 0; ip < l_p.size(); ip++) {
        // one prime of each conjugate pair: the neighbors at conj(P) are
        // the reversed edges of those at P (Hecke.jl uses a single prime
        // too)
        if (l_p[ip].conj_index >= ip) {
          pending.push({l_p[ip].q(), ell, ip});
        }
      }
      primes_above[ell] = l_p;
    }
  };
  while (!closed()) {
    if (pending.empty()) {
      if (next_ell > max_prime) {
        break;
      }
      decompose_up_to(static_cast<uint64_t>(next_ell));
      continue;
    }
    // every ell <= N(P) must be decomposed before the smallest P is taken
    uint64_t top_norm = std::get<0>(pending.top());
    if (next_ell <= max_prime && static_cast<uint64_t>(next_ell) <= top_norm) {
      decompose_up_to(top_norm);
      continue;
    }
    auto [norm, ell, ip] = pending.top();
    pending.pop();
    std::vector<HermIdealPrime<Tint>> const &l_prime = primes_above[ell];
    {
      HermNeighborPrime<T, Tint> np = HermGetNeighborPrime(ord, l_prime, ip);
      res.primes_used.push_back(ell);
#ifdef DEBUG_HERMITIAN_ENUMERATION
      os << "HERMENUM: the neighbors at a prime above " << ell << ", of norm "
         << norm << "\n";
#endif
      for (size_t i = 0; i < res.l_class.size() && !closed(); i++) {
        MyMatrix<T> A = res.l_class[i].basis;
        std::vector<MyMatrix<Tint>> l_gen = res.l_class[i].l_gen;
        HermLatticeData<T, Tint> L = HermLatticeOfBasis(ord, ring, A);
        HermLineSpace ls = HermGetLineSpace(L, np.C);
        for (auto const &[idx, size] : HermLineOrbits(ls, l_gen)) {
          (void)size;
          for (auto const &Bn : HermNeighbors(L, np, ls.line_vector(idx))) {
#ifdef SANITY_CHECK_HERMITIAN_ENUMERATION
            // the neighbor has the local symbol of L at the primes above ell
            HermLatticeData<T, Tint> Ln =
                HermLatticeOfBasis(ord, ring, MyMatrix<T>(Bn * A));
            for (size_t jp = 0; jp < l_prime.size(); jp++) {
              if (HermK0Representative(l_prime, jp) != jp) {
                continue;
              }
              bool ram = l_prime[jp].type == HermIdealType::Ramified;
              bool dy = ram && ell == 2;
              int dv = dy ? HermDifferentVal(ord, l_prime, jp) : 0;
              if (!HermLocalSymbolEquivalent(
                      HermLocalSymbolOfLattice(Ln, l_prime, jp),
                      HermLocalSymbolOfLattice(L, l_prime, jp), ram, dy, dv)) {
                std::cerr << "HERMENUM: a neighbor at " << l_prime[ip].name
                          << " is not in the genus\n";
                throw TerminalException{1};
              }
            }
#endif
            insert(MyMatrix<T>(Bn * A));
            if (closed()) {
              break;
            }
          }
          if (closed()) {
            break;
          }
        }
      }
    }
  }
  res.complete = closed();
  return res;
}

// clang-format off
#endif  // SRC_GENERA_HERMITIANENUMERATION_H_
// clang-format on

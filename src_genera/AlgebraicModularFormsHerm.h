// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ALGEBRAICMODULARFORMSHERM_H_
#define SRC_GENERA_ALGEBRAICMODULARFORMSHERM_H_

// clang-format off
#include "AlgebraicModularFormsZ.h"
#include "HermitianEnumeration.h"
#include "HermitianGenusMass.h"
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
# define DEBUG_ALGEBRAIC_MODULAR_FORMS_HERM
#endif

#ifdef TIMINGS
# define TIMINGS_ALGEBRAIC_MODULAR_FORMS_HERM
#endif

#ifdef SANITY_CHECK
# define SANITY_CHECK_ALGEBRAIC_MODULAR_FORMS_HERM
#endif

/*
  Algebraic modular forms of trivial weight for the unitary group of a
  genus of positive definite Hermitian lattices over a CM field E, from the
  list of its classes (GENUS_HermitianEnumerate, OutFormat CPP, read by
  GENUS_MPI_AlgebraicModularForms), with a level structure at an ideal n
  of O_E or none: the analogue over E of AlgebraicModularFormsZ.h, whose
  result type and checks it shares.

  The states are the pairs (L_j, O), O an orbit of U(L_j) on the level
  lines of L_j / n L_j, a single state per class without a level. A level
  line is a free summand O x + n L of rank one of L / n L (the submodule
  O x + n L / n L is isomorphic to O / n) that is the reduction of an
  isotropic vector, the analogue of the cusp of Gamma_0(n): h(x, x) is in
  the trace ideal Tr_{E/K_0}(n) at the primes under n, the ideal by which
  h(x, x) moves when x moves by n L. At a split prime the condition is
  empty (all the lines), at an inert or odd ramified prime it is h(x, x)
  in P; at a dyadic ramified prime the trace ideal can be smaller than P
  cap K_0. At a prime level the lines are those of the projective space of
  L / P L (HermLineSpace); otherwise they are enumerated in L / n L and
  keyed by the Hermite normal form of O x + n L.

  The Hecke operator T_P at a prime P of E where the genus is unimodular:
  the entry (L_i, O) -> (L_j, O') counts the pairs (x, l), x a line of
  L_i / C L_i and M one of its neighbors (HermitianNeighbor.h: one for a
  split prime, several for an inert or ramified one) with M ~ L_j, l a line
  of O carried into O' by the isometry M -> L_j (M and L_i agree away from
  P, so M / N M = L_i / N L_i). Summed over the orbits X of U(L_i) on the
  lines, it is |X| / |O| times the number of l in O sent to O'. The
  neighbors are recognized by the decision tree of LatticeInvariantFamily.h
  on the pairs (G, G Omega^T) of HermPairOfBasis, which separates the
  Hermitian classes. The adjoint of T_P for the inner product |O| / |U(L)|
  of the states is T_conj(P): T_P is self-adjoint at an inert or ramified
  prime, and T_P, T_conj(P) are both computed at a split one.
 */

// The pair of HermPairOfBasis in Tint (both are integral).
template <typename T, typename Tint>
std::pair<MyMatrix<Tint>, MyMatrix<Tint>>
HermPairOfBasisInt(RealifiedRing<T, Tint> const &ring, MyMatrix<T> const &A) {
  std::vector<MyMatrix<T>> l_mat = HermPairOfBasis(ring, A);
  return {UniversalMatrixConversion<Tint, T>(l_mat[0]),
          UniversalMatrixConversion<Tint, T>(l_mat[1])};
}

/*
  A level n = prod P_i^{e_i}: the Z-basis of n (in the coordinates of the
  basis of O_E), its norm N(n), which kills O_E / n, and for each P_i the
  valuation of Tr_{E/K_0}(n) O_E at P_i.
 */
template <typename Tint> struct HermLevel {
  std::vector<HermIdealPrime<Tint>> l_prime;
  std::vector<int> l_exp;
  MyMatrix<Tint> ideal;
  int64_t norm = 1;
  std::vector<int> l_tau;
  std::string name;
  bool is_prime() const { return l_prime.size() == 1 && l_exp[0] == 1; }
};

template <typename T, typename Tint>
HermLevel<Tint> HermGetLevel(CMOrder<T, Tint> const &ord,
                             std::vector<HermIdealPrime<Tint>> const &l_prime,
                             std::vector<int> const &l_exp) {
  int d = ord.deg;
  HermLevel<Tint> lev;
  lev.l_prime = l_prime;
  lev.l_exp = l_exp;
  lev.ideal = IdentityMat<Tint>(d);
  for (size_t i = 0; i < l_prime.size(); i++) {
    lev.name += (i > 0 ? " " : "") + l_prime[i].name;
    if (l_exp[i] > 1) {
      lev.name += "^" + std::to_string(l_exp[i]);
    }
    for (int k = 0; k < l_exp[i]; k++) {
      MyMatrix<Tint> const &P = l_prime[i].ideal;
      MyMatrix<Tint> Gen(lev.ideal.rows() * P.rows(), d);
      for (int a = 0; a < lev.ideal.rows(); a++) {
        for (int b = 0; b < P.rows(); b++) {
          MyVector<Tint> pr =
              ord.product(GetMatrixRow(lev.ideal, a), GetMatrixRow(P, b));
          for (int c = 0; c < d; c++) {
            Gen(a * P.rows() + b, c) = pr(c);
          }
        }
      }
      lev.ideal = GetZbasis(Gen);
      lev.norm *= static_cast<int64_t>(l_prime[i].q());
    }
  }
  // Tr(n) is the Z-span of the traces of a Z-basis (it is an O_{K_0}
  // module), so its valuation is the least one of these traces.
  for (auto const &P : l_prime) {
    int tau = 1 << 28;
    for (int a = 0; a < lev.ideal.rows(); a++) {
      MyVector<Tint> b = GetMatrixRow(lev.ideal, a);
      MyVector<Tint> tr = b + ord.conjugate(b);
      tau = std::min(tau, HermIdealVal<T, Tint>(ord, P,
                                                UniversalVectorConversion<T, Tint>(tr)));
    }
    lev.l_tau.push_back(tau);
  }
  return lev;
}

// Whether h(x, x) is in Tr(n) at the primes of n: x is isotropic mod n.
template <typename T, typename Tint>
bool HermLevelIsotropic(HermLatticeData<T, Tint> const &L,
                        HermLevel<Tint> const &lev, std::vector<int64_t> const &z) {
  int n = L.n();
  MyVector<Tint> x(n);
  for (int k = 0; k < n; k++) {
    x(k) = UniversalScalarConversion<Tint, int64_t>(z[k]);
  }
  MyVector<T> hxx = L.herm(x, x);
  for (size_t i = 0; i < lev.l_prime.size(); i++) {
    if (HermIdealVal<T, Tint>(*L.ord, lev.l_prime[i], hxx) < lev.l_tau[i]) {
      return false;
    }
  }
  return true;
}

/*
  The lines of L / n L: at a prime level the projective space of L / P L
  (every line, the level ones found by HermLevelIsotropic), otherwise the
  level lines only, enumerated over the N(n)^m vectors of L / n L and keyed
  by the Hermite normal form of O x + n L (HnfKeyModM, modulo N(n)).
 */
struct HermLevelSpace {
  bool general = false;
  HermLineSpace ls;
  int n = 0;
  int64_t M = 1;
  // n L in the coordinates of L, the multiplications by the basis of O_E,
  // modulo M
  std::vector<std::vector<int64_t>> rows_nL;
  std::vector<std::vector<int64_t>> l_mult;
  std::vector<std::vector<int64_t>> l_line;
  std::map<std::vector<int64_t>, size_t> map_line;
  int64_t modulus() const { return general ? M : ls.ell(); }
  size_t n_line() const { return general ? l_line.size() : ls.n_line(); }
  std::vector<int64_t> line_vector(size_t const &idx) const {
    return general ? l_line[idx] : ls.line_vector(idx);
  }
  std::vector<int64_t> line_key(std::vector<int64_t> const &z) const {
    std::vector<std::vector<int64_t>> rows = rows_nL;
    for (auto const &Om : l_mult) {
      std::vector<int64_t> zO(n, 0);
      for (int j = 0; j < n; j++) {
        int64_t sum = 0;
        for (int k = 0; k < n; k++) {
          sum = (sum + ResidueModPrime(z[k], M) * Om[k * n + j]) % M;
        }
        zO[j] = sum;
      }
      rows.push_back(zO);
    }
    return HnfKeyModM(rows, n, M);
  }
  std::optional<size_t> line_index(std::vector<int64_t> const &z) const {
    if (!general) {
      return ls.line_index(z);
    }
    auto iter = map_line.find(line_key(z));
    if (iter == map_line.end()) {
      return {};
    }
    return iter->second;
  }
};

template <typename T, typename Tint>
HermLevelSpace HermGetLevelSpace(HermLatticeData<T, Tint> const &L,
                                 HermLevel<Tint> const &lev) {
  HermLevelSpace res;
  int n = L.n();
  res.n = n;
  if (lev.is_prime()) {
    res.ls = HermGetLineSpace(L, lev.l_prime[0]);
    return res;
  }
  res.general = true;
  int64_t M = lev.norm;
  res.M = M;
  MyMatrix<Tint> nL = HermIdealTimesLattice(L, lev.ideal);
  for (int i = 0; i < nL.rows(); i++) {
    std::vector<int64_t> r(n);
    for (int j = 0; j < n; j++) {
      r[j] = order_residue(nL(i, j), M);
    }
    res.rows_nL.push_back(r);
  }
  for (auto const &Om : L.Omega) {
    std::vector<int64_t> om(n * n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        om[i * n + j] = order_residue(Om(i, j), M);
      }
    }
    res.l_mult.push_back(om);
  }
  // a transversal of L / n L from the diagonal of the key of n L, and the
  // free lines: the index of O x + n L is N(n)^{m - 1}
  std::vector<int64_t> H = HnfKeyModM(res.rows_nL, n, M);
  int d = L.ord->deg;
  int m = n / d;
  int64_t target = 1;
  for (int k = 0; k < m - 1; k++) {
    target *= M;
  }
  std::vector<int64_t> diag(n);
  int64_t total = 1;
  for (int i = 0; i < n; i++) {
    diag[i] = H[i * n + i];
    total *= diag[i];
  }
  std::vector<int64_t> z(n, 0);
  for (int64_t idx = 0; idx < total; idx++) {
    int64_t w = idx;
    for (int i = 0; i < n; i++) {
      z[i] = w % diag[i];
      w /= diag[i];
    }
    std::vector<int64_t> key = res.line_key(z);
    int64_t index = 1;
    for (int i = 0; i < n; i++) {
      index *= key[i * n + i];
    }
    if (index != target || res.map_line.count(key) == 1) {
      continue;
    }
    if (!HermLevelIsotropic(L, lev, z)) {
      continue;
    }
    res.map_line[key] = res.l_line.size();
    res.l_line.push_back(z);
  }
  return res;
}

/*
  The level orbits of a lattice: for each line of ls (the lines of L / n
  L), its orbit under the generators, -1 for a line that is not a level
  line, and the size of every orbit; no orbit when L / N L has no level
  line (the caller reports it).
 */
template <typename T, typename Tint>
AmfZLevelData HermAmfLevelOrbits(HermLatticeData<T, Tint> const &L,
                                 HermLevelSpace const &ls,
                                 HermLevel<Tint> const &lev,
                                 std::vector<MyMatrix<Tint>> const &l_gen) {
  int n = L.n();
  int64_t ell = ls.modulus();
  size_t n_line = ls.n_line();
  // -2: a level line whose orbit is not known yet (the lines of a general
  // level space are all level lines)
  AmfZLevelData data;
  data.orbit_of_line.assign(n_line, -1);
  for (size_t idx = 0; idx < n_line; idx++) {
    if (ls.general || HermLevelIsotropic(L, lev, ls.line_vector(idx))) {
      data.orbit_of_line[idx] = -2;
    }
  }
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
  std::vector<size_t> stack;
  std::vector<int64_t> w(n);
  for (size_t idx = 0; idx < n_line; idx++) {
    if (data.orbit_of_line[idx] != -2) {
      continue;
    }
    int64_t o = static_cast<int64_t>(data.orbit_size.size());
    data.orbit_of_line[idx] = o;
    stack.push_back(idx);
    size_t size = 0;
    while (!stack.empty()) {
      size_t cur = stack.back();
      stack.pop_back();
      size++;
      std::vector<int64_t> z = ls.line_vector(cur);
      for (auto const &gr : gen_red) {
        for (int j = 0; j < n; j++) {
          int64_t sum = 0;
          for (int i = 0; i < n; i++) {
            sum = (sum + z[i] * gr[i * n + j]) % ell;
          }
          w[j] = sum;
        }
        std::optional<size_t> opt = ls.line_index(w);
        if (opt && data.orbit_of_line[*opt] == -2) {
          data.orbit_of_line[*opt] = o;
          stack.push_back(*opt);
        }
      }
    }
    data.orbit_size.push_back(size);
  }
  return data;
}

/*
  The row of the class i at the Hecke prime np. l_basis are the bases of
  the classes, family their decision tree (on the pairs), l_gen_i the
  generators of U(L_i). With a level, N is its ideal, level_ls the line
  spaces of L_j / N L_j and level the level orbits of every class.
 */
template <typename T, typename Tint, typename Tgroup>
AmfZRow HermAmfComputeRow(CMOrder<T, Tint> const &ord,
                          RealifiedRing<T, Tint> const &ring,
                          std::vector<MyMatrix<T>> const &l_basis, size_t i,
                          HermNeighborPrime<T, Tint> const &np,
                          InvariantFamily<Tint> const &family,
                          std::vector<MyMatrix<Tint>> const &l_gen_i,
                          std::optional<HermLevel<Tint>> const &N,
                          std::vector<HermLevelSpace> const &level_ls,
                          std::vector<AmfZLevelData> const &level,
                          std::ostream &os) {
  MyMatrix<T> const &A_i = l_basis[i];
  HermLatticeData<T, Tint> L = HermLatticeOfBasis(ord, ring, A_i);
  int n = L.n();
#ifdef TIMINGS_ALGEBRAIC_MODULAR_FORMS_HERM
  MicrosecondTime time;
#endif
  HermLineSpace ls = HermGetLineSpace(L, np.C);
  std::vector<std::pair<size_t, size_t>> l_orbit = HermLineOrbits(ls, l_gen_i);
#ifdef TIMINGS_ALGEBRAIC_MODULAR_FORMS_HERM
  os << "|AMFH: class " << i << " " << np.P.name << " line orbits|=" << time
     << "\n";
#endif
  // the sums of |X| times the number of lines of O sent to O2, divided by
  // |O| at the end
  std::map<std::tuple<size_t, size_t, size_t>, int64_t> sums;
  AmfZRow row;
  for (auto const &[idx, size] : l_orbit) {
    int64_t size_X = static_cast<int64_t>(size);
    for (auto const &Bn : HermNeighbors(L, np, ls.line_vector(idx))) {
      row.n_neighbor += size;
      MyMatrix<T> A_M = Bn * A_i;
      std::pair<MyMatrix<Tint>, MyMatrix<Tint>> GS = HermPairOfBasisInt(ring, A_M);
      size_t j =
          family_invariant_index<Tint, Tgroup>(family, GS.first, os, GS.second);
      if (!N) {
        sums[{0, j, 0}] += size_X;
        continue;
      }
      // An isometry Q: M -> L_j, Q {G_j, S_j} Q^T = {G_M, S_M}, then the
      // level lines of L_i in the coordinates of L_j: z -> z Bn^{-1} Q, the
      // denominators of Bn^{-1} being prime to N.
      std::vector<MyMatrix<T>> l_j = HermPairOfBasis(ring, l_basis[j]);
      std::vector<MyMatrix<T>> l_M = HermPairOfBasis(ring, A_M);
      std::optional<MyMatrix<Tint>> opt_iso =
          ArithmeticEquivalenceMultiple<T, Tint, Tgroup>(l_j, l_M, os);
#ifdef SANITY_CHECK_ALGEBRAIC_MODULAR_FORMS_HERM
      if (!opt_iso) {
        std::cerr << "AMFH: the neighbor recognized as the class " << j
                  << " is not isometric to it\n";
        throw TerminalException{1};
      }
#endif
      MyMatrix<T> Q = UniversalMatrixConversion<T, Tint>(*opt_iso);
      auto maps = [&](MyMatrix<T> const &U) -> bool {
        return U * l_j[0] * U.transpose() == l_M[0] &&
               U * l_j[1] * U.transpose() == l_M[1];
      };
      if (!maps(Q)) {
        // the isometry the other way
        Q = Inverse(Q);
#ifdef SANITY_CHECK_ALGEBRAIC_MODULAR_FORMS_HERM
        if (!maps(Q)) {
          std::cerr << "AMFH: the isometry does not map the class to the "
                       "neighbor\n";
          throw TerminalException{1};
        }
#endif
      }
      MyMatrix<T> R = Inverse(Bn) * Q;
      int64_t ellN = level_ls[i].modulus();
      std::vector<int64_t> R_red(n * n);
      for (int a = 0; a < n; a++) {
        for (int b = 0; b < n; b++) {
          R_red[a * n + b] = order_residue_rat(R(a, b), ellN);
        }
      }
      AmfZLevelData const &level_i = level[i];
      AmfZLevelData const &level_j = level[j];
      std::map<std::pair<size_t, size_t>, int64_t> count;
      std::vector<int64_t> w(n);
      for (size_t li = 0; li < level_i.orbit_of_line.size(); li++) {
        int64_t o = level_i.orbit_of_line[li];
        if (o < 0) {
          continue;
        }
        std::vector<int64_t> z = level_ls[i].line_vector(li);
        for (int b = 0; b < n; b++) {
          int64_t sum = 0;
          for (int a = 0; a < n; a++) {
            sum = (sum + z[a] * R_red[a * n + b]) % ellN;
          }
          w[b] = sum;
        }
        std::optional<size_t> opt = level_ls[j].line_index(w);
#ifdef SANITY_CHECK_ALGEBRAIC_MODULAR_FORMS_HERM
        if (!opt || level_j.orbit_of_line[*opt] < 0) {
          std::cerr << "AMFH: a level line is not sent to a level line\n";
          throw TerminalException{1};
        }
#endif
        count[{static_cast<size_t>(o),
               static_cast<size_t>(level_j.orbit_of_line[*opt])}]++;
      }
      for (auto const &[key, c] : count) {
        sums[{key.first, j, key.second}] += size_X * c;
      }
    }
  }
  for (auto const &[key, s] : sums) {
    int64_t size_O =
        static_cast<int64_t>(level[i].orbit_size[std::get<0>(key)]);
#ifdef SANITY_CHECK_ALGEBRAIC_MODULAR_FORMS_HERM
    if (s % size_O != 0) {
      std::cerr << "AMFH: a Hecke matrix entry is not an integer\n";
      throw TerminalException{1};
    }
#endif
    row.entries[key] = s / size_O;
  }
#ifdef TIMINGS_ALGEBRAIC_MODULAR_FORMS_HERM
  os << "|AMFH: class " << i << " " << np.P.name << ", " << l_orbit.size()
     << " orbits|=" << time << "\n";
#endif
#ifdef DEBUG_ALGEBRAIC_MODULAR_FORMS_HERM
  os << "AMFH: class " << i << " " << np.P.name << " done, "
     << l_orbit.size() << " orbits of lines, " << row.n_neighbor
     << " neighbors\n";
#endif
  return row;
}

// clang-format off
#endif  // SRC_GENERA_ALGEBRAICMODULARFORMSHERM_H_
// clang-format on

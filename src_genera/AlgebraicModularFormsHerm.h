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
  GENUS_MPI_AlgebraicModularForms), with a
  level structure at a prime N of E or none: the analogue over E of
  AlgebraicModularFormsZ.h, whose result type and checks it shares.

  The states are the pairs (L_j, O), O an orbit of U(L_j) on the level
  lines of L_j / N L_j: all the lines when N is split over K_0 (U(L_j) acts
  there through GL_m), the isotropic ones (h(x, x) in N) when N = conj(N);
  a single state per class without a level.

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
  The level orbits of a lattice: for each line of ls (the lines of L / N
  L), its orbit under the generators, -1 for a line that is not a level
  line, and the size of every orbit; no orbit when L / N L has no level
  line (the caller reports it).
 */
template <typename T, typename Tint>
AmfZLevelData HermAmfLevelOrbits(HermLatticeData<T, Tint> const &L,
                                 HermLineSpace const &ls,
                                 HermIdealPrime<Tint> const &N,
                                 std::vector<MyMatrix<Tint>> const &l_gen) {
  int n = L.n();
  int64_t ell = N.ell;
  bool all_lines = N.type == HermIdealType::Split;
  size_t n_line = ls.n_line();
  // -2: a level line whose orbit is not known yet
  AmfZLevelData data;
  data.orbit_of_line.assign(n_line, -1);
  MyVector<Tint> zv(n);
  for (size_t idx = 0; idx < n_line; idx++) {
    bool level = all_lines;
    if (!level) {
      std::vector<int64_t> z = ls.line_vector(idx);
      for (int k = 0; k < n; k++) {
        zv(k) = UniversalScalarConversion<Tint, int64_t>(z[k]);
      }
      level = HermResidueRat(N, L.herm(zv, zv)).is_zero();
    }
    if (level) {
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
  generators of U(L_i). With a level, N is its prime, level_ls the line
  spaces of L_j / N L_j and level the level orbits of every class.
 */
template <typename T, typename Tint, typename Tgroup>
AmfZRow HermAmfComputeRow(CMOrder<T, Tint> const &ord,
                          RealifiedRing<T, Tint> const &ring,
                          std::vector<MyMatrix<T>> const &l_basis, size_t i,
                          HermNeighborPrime<T, Tint> const &np,
                          InvariantFamily<Tint> const &family,
                          std::vector<MyMatrix<Tint>> const &l_gen_i,
                          std::optional<HermIdealPrime<Tint>> const &N,
                          std::vector<HermLineSpace> const &level_ls,
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
      int64_t ellN = N->ell;
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

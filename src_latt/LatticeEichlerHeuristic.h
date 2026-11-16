// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEEICHLERHEURISTIC_H_
#define SRC_LATT_LATTICEEICHLERHEURISTIC_H_

// clang-format off
#include "LatticeEichlerDecomposition.h"
#include "LatticeMethods.h"
#include <algorithm>
#include <map>
#include <optional>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_LATTICE_EICHLER_HEURISTIC
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_LATTICE_EICHLER_HEURISTIC
#endif

/*
  The canonical form, the automorphism group and the isometry test of a
  positive definite lattice through its Eichler decomposition
  (LatticeEichlerDecomposition.h), the method heuristic of LatticeMethods.h
  being applied to each indecomposable summand. An option, not used by
  default.

  The decomposition is unique, so it is an isometry invariant, and the
  results of isometric lattices remain comparable:
  * the canonical form: the summands are put in canonical form by
    heuristic, ordered by their canonical Gram matrices (rank first), and
    stacked; the Gram matrix of L is the block diagonal one of the sorted
    canonical summands. For an indecomposable L it is the canonical form of
    heuristic itself;
  * the automorphism group: an isometry of L permutes the summands, so it is
    the product over the isometry classes of summands K, occurring m times,
    of the wreath products Aut(K) wr S_m: order prod |Aut K|^m m!, the
    generators being those of Aut(K) on the first copy and the permutations
    of the copies, which the canonical bases make identical;
  * the isometry: P = E_2^{-1} E_1 for the canonical bases E_i, when the
    canonical forms agree.

  The gain is on the lattices where heuristic chooses a method for the
  whole lattice that fails on one summand: Leech + A_1 has roots, so
  heuristic takes the Weyl-chamber method, whose family on the root-free
  part explodes, while the summands Leech and A_1 are each easy.
 */

template <typename T, typename Tint> struct EichlerCanonicalSummand {
  // rows: a canonical basis of the summand, in the coordinates of L
  MyMatrix<Tint> CanBasis;
  MyMatrix<T> CanGram;
};

template <typename T>
bool EichlerGramLess(MyMatrix<T> const &A, MyMatrix<T> const &B) {
  if (A.rows() != B.rows()) {
    return A.rows() < B.rows();
  }
  int n = A.rows();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (A(i, j) != B(i, j)) {
        return A(i, j) < B(i, j);
      }
    }
  }
  return false;
}

/*
  The summands with an LLL-reduced basis. The basis given by the
  decomposition is a Hermite normal form, and the methods are much slower
  from it than from a reduced one, although they reduce their input
  themselves (measured on the classes of rank 25 and determinant 2).
 */
template <typename T, typename Tint>
std::vector<EichlerComponent<T, Tint>>
EichlerReducedComponents(MyMatrix<T> const &G, std::ostream &os) {
  std::vector<EichlerComponent<T, Tint>> ListComp =
      EichlerDecomposition<T, Tint>(G, os);
  if (ListComp.size() == 1) {
    // indecomposable: the callers then work on G itself
    return ListComp;
  }
  for (auto &comp : ListComp) {
    LLLreduction<T, Tint> red = LLLreducedBasis<T, Tint>(comp.Gram, os);
    comp.Basis = red.Pmat * comp.Basis;
    comp.Gram = red.GramMatRed;
  }
#ifdef DEBUG_LATTICE_EICHLER_HEURISTIC
  os << "EICHLER_HEURISTIC: " << ListComp.size() << " summands of ranks";
  for (auto const &comp : ListComp) {
    os << " " << comp.Gram.rows();
  }
  os << "\n";
#endif
  return ListComp;
}

template <typename T, typename Tint, typename Tgroup>
EichlerCanonicalSummand<T, Tint>
EichlerCanonicalOfComponent(EichlerComponent<T, Tint> const &comp,
                            std::ostream &os) {
  MyMatrix<Tint> C = LatticeCanonicalBasis<T, Tint, Tgroup>(
      comp.Gram, LatticeStabCanIsomMethod::Heuristic, os);
  EichlerCanonicalSummand<T, Tint> s;
  s.CanBasis = C * comp.Basis;
  MyMatrix<T> C_T = UniversalMatrixConversion<T, Tint>(C);
  s.CanGram = C_T * comp.Gram * C_T.transpose();
  return s;
}

// The stacked bases: rows, a basis of L.
template <typename Tint>
MyMatrix<Tint> EichlerStackedBasis(std::vector<MyMatrix<Tint>> const &ListB,
                                   int n) {
  MyMatrix<Tint> E(n, n);
  int pos = 0;
  for (auto const &B : ListB) {
    for (int i = 0; i < B.rows(); i++) {
      for (int j = 0; j < n; j++) {
        E(pos + i, j) = B(i, j);
      }
    }
    pos += B.rows();
  }
  return E;
}

/*
  A canonical basis: rows B with B G B^T canonical. An indecomposable
  lattice gets the canonical basis of heuristic, computed from G itself;
  otherwise the summands are put in canonical form, sorted by their
  canonical Gram matrices and stacked.
 */
template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint> EichlerHeuristicCanonicalBasis(MyMatrix<T> const &G,
                                              std::ostream &os) {
  std::vector<EichlerComponent<T, Tint>> ListComp =
      EichlerReducedComponents<T, Tint>(G, os);
  if (ListComp.size() == 1) {
    return LatticeCanonicalBasis<T, Tint, Tgroup>(
        G, LatticeStabCanIsomMethod::Heuristic, os);
  }
  std::vector<EichlerCanonicalSummand<T, Tint>> ListCan;
  for (auto const &comp : ListComp) {
    ListCan.push_back(EichlerCanonicalOfComponent<T, Tint, Tgroup>(comp, os));
  }
  std::stable_sort(ListCan.begin(), ListCan.end(),
                   [](EichlerCanonicalSummand<T, Tint> const &a,
                      EichlerCanonicalSummand<T, Tint> const &b) -> bool {
                     return EichlerGramLess(a.CanGram, b.CanGram);
                   });
  std::vector<MyMatrix<Tint>> ListB;
  for (auto const &s : ListCan) {
    ListB.push_back(s.CanBasis);
  }
  return EichlerStackedBasis(ListB, G.rows());
}

/*
  The automorphism group. An indecomposable lattice gets that of heuristic.
  Otherwise the summands are grouped in isometry classes: those alone of
  their rank and determinant form a class by themselves and keep their
  reduced basis; the others are put in canonical form, and the equal
  canonical Gram matrices make the classes. Each class K, of multiplicity
  m, contributes Aut(K) on its first copy and the permutations of the
  copies, which have identical Gram matrices in their bases.
 */
template <typename T, typename Tint, typename Tgroup>
IntAutomorphismGroup<Tint, Tgroup>
EichlerHeuristicAutomorphismGroup(MyMatrix<T> const &G, std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
  int n = G.rows();
  std::vector<EichlerComponent<T, Tint>> ListComp =
      EichlerReducedComponents<T, Tint>(G, os);
  if (ListComp.size() == 1) {
    std::vector<MyMatrix<T>> ListMat{G};
    return LatticeAutomorphismGroup<T, Tint, Tgroup>(
        ListMat, LatticeStabCanIsomMethod::Heuristic, os);
  }
  // the isometry classes: a list of copies (basis, Gram), all copies of one
  // class with the same Gram matrix
  std::map<std::pair<int, T>, std::vector<size_t>> by_rank_det;
  for (size_t i = 0; i < ListComp.size(); i++) {
    by_rank_det[{static_cast<int>(ListComp[i].Gram.rows()),
                 DeterminantMat(ListComp[i].Gram)}]
        .push_back(i);
  }
  std::vector<std::pair<MyMatrix<T>, std::vector<MyMatrix<Tint>>>> classes;
  for (auto const &[key, members] : by_rank_det) {
    if (members.size() == 1) {
      EichlerComponent<T, Tint> const &comp = ListComp[members[0]];
      classes.push_back({comp.Gram, {comp.Basis}});
      continue;
    }
    std::vector<EichlerCanonicalSummand<T, Tint>> ListCan;
    for (size_t i : members) {
      ListCan.push_back(
          EichlerCanonicalOfComponent<T, Tint, Tgroup>(ListComp[i], os));
    }
    std::stable_sort(ListCan.begin(), ListCan.end(),
                     [](EichlerCanonicalSummand<T, Tint> const &a,
                        EichlerCanonicalSummand<T, Tint> const &b) -> bool {
                       return EichlerGramLess(a.CanGram, b.CanGram);
                     });
    size_t i = 0;
    while (i < ListCan.size()) {
      std::vector<MyMatrix<Tint>> copies{ListCan[i].CanBasis};
      size_t j = i + 1;
      while (j < ListCan.size() && ListCan[j].CanGram == ListCan[i].CanGram) {
        copies.push_back(ListCan[j].CanBasis);
        j++;
      }
      classes.push_back({ListCan[i].CanGram, copies});
      i = j;
    }
  }
  std::vector<MyMatrix<Tint>> ListB;
  for (auto const &cl : classes) {
    for (auto const &B : cl.second) {
      ListB.push_back(B);
    }
  }
  MyMatrix<Tint> E = EichlerStackedBasis(ListB, n);
  MyMatrix<Tint> Einv = Inverse(E);
  IntAutomorphismGroup<Tint, Tgroup> res;
  res.order = TintGroup(1);
  // g_E, preserving the block diagonal Gram matrix E G E^T, gives
  // E^{-1} g_E E preserving G
  auto push = [&](MyMatrix<Tint> const &gE) -> void {
    res.ListGen.push_back(Einv * gE * E);
  };
  int pos = 0;
  for (auto const &[GramK, copies] : classes) {
    int r = GramK.rows();
    size_t m = copies.size();
    std::vector<MyMatrix<T>> ListMat{GramK};
    IntAutomorphismGroup<Tint, Tgroup> aut =
        LatticeAutomorphismGroup<T, Tint, Tgroup>(
            ListMat, LatticeStabCanIsomMethod::Heuristic, os);
    for (size_t k = 0; k < m; k++) {
      res.order *= aut.order;
    }
    for (size_t k = 2; k <= m; k++) {
      res.order *= TintGroup(static_cast<long>(k));
    }
    // Aut(K) on the first copy
    for (auto const &h : aut.ListGen) {
      MyMatrix<Tint> gE = IdentityMat<Tint>(n);
      for (int a = 0; a < r; a++) {
        for (int b = 0; b < r; b++) {
          gE(pos + a, pos + b) = h(a, b);
        }
      }
      push(gE);
    }
    // the permutations of the copies: a transposition and an m-cycle
    auto block_permutation = [&](std::vector<size_t> const &sigma) -> void {
      MyMatrix<Tint> gE = IdentityMat<Tint>(n);
      for (size_t k = 0; k < m; k++) {
        int src = pos + static_cast<int>(k) * r;
        int dst = pos + static_cast<int>(sigma[k]) * r;
        for (int a = 0; a < r; a++) {
          gE(src + a, src + a) = 0;
        }
        for (int a = 0; a < r; a++) {
          gE(src + a, dst + a) = 1;
        }
      }
      push(gE);
    };
    if (m >= 2) {
      std::vector<size_t> transposition(m), cycle(m);
      for (size_t k = 0; k < m; k++) {
        transposition[k] = k;
        cycle[k] = (k + 1) % m;
      }
      std::swap(transposition[0], transposition[1]);
      block_permutation(transposition);
      if (m >= 3) {
        block_permutation(cycle);
      }
    }
    pos += static_cast<int>(m) * r;
  }
#ifdef SANITY_CHECK_LATTICE_EICHLER_HEURISTIC
  for (auto const &g : res.ListGen) {
    MyMatrix<T> g_T = UniversalMatrixConversion<T, Tint>(g);
    if (g_T * G * g_T.transpose() != G) {
      std::cerr << "EICHLER_HEURISTIC: a generator does not preserve G\n";
      throw TerminalException{1};
    }
  }
#endif
  return res;
}

/*
  An isometry P with P G1 P^T = G2, if there is one. Two indecomposable
  lattices go to the isometry test of heuristic; otherwise the canonical
  bases E_i give P = E_2^{-1} E_1 when the canonical forms agree.
 */
template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
EichlerHeuristicIsometry(MyMatrix<T> const &G1, MyMatrix<T> const &G2,
                         std::ostream &os) {
  if (G1.rows() != G2.rows()) {
    return {};
  }
  std::vector<EichlerComponent<T, Tint>> ListComp1 =
      EichlerReducedComponents<T, Tint>(G1, os);
  std::vector<EichlerComponent<T, Tint>> ListComp2 =
      EichlerReducedComponents<T, Tint>(G2, os);
  if (ListComp1.size() != ListComp2.size()) {
    return {};
  }
  if (ListComp1.size() == 1) {
    std::vector<MyMatrix<T>> LM1{G1}, LM2{G2};
    return LatticeIsometry<T, Tint, Tgroup>(
        LM1, LM2, {}, LatticeStabCanIsomMethod::Heuristic, os);
  }
  MyMatrix<Tint> E1 = EichlerHeuristicCanonicalBasis<T, Tint, Tgroup>(G1, os);
  MyMatrix<Tint> E2 = EichlerHeuristicCanonicalBasis<T, Tint, Tgroup>(G2, os);
  MyMatrix<T> E1_T = UniversalMatrixConversion<T, Tint>(E1);
  MyMatrix<T> E2_T = UniversalMatrixConversion<T, Tint>(E2);
  if (E1_T * G1 * E1_T.transpose() != E2_T * G2 * E2_T.transpose()) {
    return {};
  }
  MyMatrix<Tint> P = Inverse(E2) * E1;
#ifdef SANITY_CHECK_LATTICE_EICHLER_HEURISTIC
  MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(P);
  if (P_T * G1 * P_T.transpose() != G2) {
    std::cerr << "EICHLER_HEURISTIC: the isometry is wrong\n";
    throw TerminalException{1};
  }
#endif
  return P;
}

// clang-format off
#endif  // SRC_LATT_LATTICEEICHLERHEURISTIC_H_
// clang-format on

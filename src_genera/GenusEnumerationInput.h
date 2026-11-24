// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GENUSENUMERATIONINPUT_H_
#define SRC_GENERA_GENUSENUMERATIONINPUT_H_

// clang-format off
#include "GenusEnumeration.h"
#include "LocalDensity.h"
#include "MassGenus.h"
#include "GenusEnumerationOdd.h"
#include "GenusSymbol.h"
#include "GenusSpinor.h"
#include <optional>
#include <string>
#include <vector>
// clang-format on

/*
  The input of an enumeration, shared by GENUS_Enumerate and
  GENUS_MPI_Enumerate: the genus from its symbol, a seed built from it,
  the mass, the spinor genera of the enumerated genus and, from them, the
  prime of out_prime and the jump primes of two_plus_jumps. Every step is
  deterministic, so all the processes of the MPI program compute the same
  input.
 */
template <typename T, typename Tgroup> struct GenusEnumerationInput {
  using Tmass = typename GenusEnumerationResult<T, Tgroup>::Tmass;
  GenusSpec spec;
  std::vector<MyMatrix<T>> ListSeed;
  Tmass TotalMass;
  int prime;
  bool even;
  std::string OutFormat;
  std::string OutFile;
  int max_runtime_second;
};

template <typename T, typename Tint, typename Tgroup>
GenusEnumerationInput<T, Tgroup>
ReadGenusEnumerationInput(FullNamelist const &eFull, std::ostream &os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  std::string GenusSymbolStr = BlockDATA.get_string("GenusSymbol");
  std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  int max_runtime_second = BlockSYSTEM.get_int("max_runtime_second");
  GenusSpec spec = ReadGenusSpec(eFull);
  // The input is the genus, as a symbol (GenusSymbol.h); a lattice is
  // turned into its symbol by GENUS_GetGenusFromLattice. The seed is built
  // from the symbol, in mpq_class / mpz_class whatever T is.
  if (GenusSymbolStr == "unset") {
    std::cerr << "GENUS_Enumerate: GenusSymbol is needed (GENUS_"
                 "GetGenusFromLattice gives the symbol of a lattice)\n";
    throw TerminalException{1};
  }
  GenusSymbol symbol = ParseGenusSymbol<mpz_class>(GenusSymbolStr);
  MyMatrix<mpq_class> Seed =
      GenusSymbolRepresentative<mpq_class, mpz_class>(symbol, os);
  os << "GENUS_Enumerate: the seed built from the genus symbol "
     << GenusSymbolString(symbol) << " is\n";
  WriteMatrix(os, Seed);
  std::vector<MyMatrix<T>> ListSeed{UniversalMatrixConversion<T, mpq_class>(Seed)};
  // A sum of 1 / |Aut(L)|, so a fraction over the type the groups count
  // with, which has nothing to do with the type of the Gram matrices.
  using Tmass = typename GenusEnumerationResult<T, Tgroup>::Tmass;
  int rank = symbol.rank;
  T det = DeterminantMat(ListSeed[0]);
  // Even genera are enumerated directly; odd ones through the genus of the
  // even sublattices (GenusEnumerationOdd.h).
  bool even = symbol.even;
  if (!even) {
    os << "GENUS_Enumerate: the genus is odd, it is enumerated through the "
          "genus of the even sublattices\n";
  }
  // The mass, computed from the first seed by the mass formula of Conway
  // and Sloane (MassGenus.h); the enumeration stops when sum 1/|Aut|
  // reaches it, and that equality certifies completeness. Tmass is
  // mpq_class, the type GenusMass computes in. The genus of an odd L has the
  // classes, automorphism groups and mass of the genus of L(2), which is
  // even.
  MyMatrix<mpq_class> GramMass =
      UniversalMatrixConversion<mpq_class, T>(ListSeed[0]);
  if (!even) {
    GramMass *= mpq_class(2);
  }
  Tmass TotalMass = GenusMass<mpq_class, mpz_class>(GramMass);
  os << "GENUS_Enumerate: the mass of the genus is " << TotalMass << "\n";
  // The spinor genera of the enumerated genus (GenusSpinor.h). For odd
  // seeds the enumerated genus is that of the even sublattices, of
  // determinant 4 det. The rank must be at least 3 for Kneser's theorem.
  T det_enum = even ? det : T(4) * det;
  auto is_odd_prime = [](int l) -> bool {
    for (int f = 3; f * f <= l; f += 2) {
      if (l % f == 0) {
        return false;
      }
    }
    return l >= 3 && l % 2 == 1;
  };
  std::optional<SpinorGenusData> spinor;
  if (rank >= 3) {
    MyMatrix<mpq_class> G_enum = UniversalMatrixConversion<mpq_class, T>(
        even ? ListSeed[0] : EvenSublatticeGram<T, Tint>(ListSeed[0]));
    spinor = GenusSpinorData<mpq_class, mpz_class>(
        GenusSymbolOfGram<mpq_class, mpz_class>(G_enum), G_enum);
    os << "GENUS_Enumerate: " << spinor->n_spinor_genera()
       << " spinor genera in the enumerated genus\n";
  }
  int prime = spec.prime;
  if (spec.strategy == NeighborStrategy::OutPrime && prime == 0) {
    prime = ChooseNeighborPrime<T>(det_enum);
    // When there are several spinor genera, the l-neighbour graph of the
    // genus is connected only for some primes l: the smallest odd one is
    // taken, if there is one, so that a single seed suffices.
    if (spinor && spinor->n_spinor_genera() > 1) {
      std::optional<int> connected;
      for (int l = 3; l < 1000 && !connected; l += 2) {
        if (is_odd_prime(l) && ResInt(det_enum, T(l)) != 0 &&
            spinor->neighbor_components<mpz_class>(l) == 1) {
          connected = l;
        }
      }
      if (connected) {
        prime = *connected;
        os << "GENUS_Enumerate: the prime " << prime
           << " connects them, it is used\n";
      } else {
        os << "GENUS_Enumerate: no single prime connects them, use "
              "two_plus_jumps, iterate or random_walk\n";
      }
    }
  }
  if (spec.strategy == NeighborStrategy::TwoPlusJumps) {
    // The jump primes: the smallest odd primes not dividing det whose
    // classes Delta(l) generate the group A of the spinor genera modulo the
    // improper kernel, chosen greedily (the independent sets form a
    // matroid, so this minimizes the largest one). Their neighbour graphs
    // together connect the genus whatever the 2-neighbour graph does. Below
    // rank 3, the smallest odd prime not dividing det.
    std::vector<int> jump_primes;
    if (spinor) {
      SpinorSubspace span = spinor->kernel();
      for (int l = 3; span.rank() < spinor->A.dim() && l < 100000; l += 2) {
        if (is_odd_prime(l) && ResInt(det_enum, T(l)) != 0 &&
            span.insert(SpinorDelta(spinor->A, mpz_class(l)))) {
          jump_primes.push_back(l);
        }
      }
    }
    // With one spinor genus (or below rank 3) any prime not dividing 2 det
    // connects the genus: the smallest one, needed when 2 divides det since
    // nothing is known of the 2-neighbour graph then.
    if (jump_primes.empty()) {
      int l = 3;
      while (!is_odd_prime(l) || ResInt(det_enum, T(l)) == 0) {
        l += 2;
      }
      jump_primes.push_back(l);
    }
    spec.rw.jump_primes = jump_primes;
    os << "GENUS_Enumerate: two_plus_jumps, the jump primes are";
    for (int l : jump_primes) {
      os << " " << l;
    }
    os << "\n";
  }
  GenusEnumerationInput<T, Tgroup> input;
  input.spec = spec;
  input.ListSeed = ListSeed;
  input.TotalMass = TotalMass;
  input.prime = prime;
  input.even = even;
  input.OutFormat = OutFormat;
  input.OutFile = OutFile;
  input.max_runtime_second = max_runtime_second;
  return input;
}

// clang-format off
#endif  // SRC_GENERA_GENUSENUMERATIONINPUT_H_
// clang-format on

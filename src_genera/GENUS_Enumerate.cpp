// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "GenusEnumeration.h"
#include "LocalDensity.h"
#include "GenusMass.h"
#include "GenusEnumerationOdd.h"
#include "GenusSymbol.h"
#include "GenusSpinor.h"
#include "SignatureSymmetric.h"
#include "Arithmetic_dispatch_types.h"
#include <optional>
#include <string>
#include <vector>
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void ProcessGenus(FullNamelist const &eFull, std::ostream &os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  std::string FileLattice = BlockDATA.get_string("FileLattice");
  std::string GenusSymbolStr = BlockDATA.get_string("GenusSymbol");
  std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  int max_runtime_second = BlockSYSTEM.get_int("max_runtime_second");
  GenusSpec spec = ReadGenusSpec(eFull);
  // The seeds: read from FileLattice, or, when it is unset, one seed built
  // from the genus symbol (GenusSymbol.h). The symbol is handled in
  // mpq_class / mpz_class whatever T is.
  std::optional<GenusSymbol> symbol;
  if (GenusSymbolStr != "unset") {
    symbol = ParseGenusSymbol<mpz_class>(GenusSymbolStr);
  }
  std::vector<MyMatrix<T>> ListSeed;
  if (FileLattice == "unset") {
    if (!symbol) {
      std::cerr << "GENUS_Enumerate: FileLattice and GenusSymbol are both "
                   "unset, one of them is needed\n";
      throw TerminalException{1};
    }
    MyMatrix<mpq_class> Seed =
        GenusSymbolRepresentative<mpq_class, mpz_class>(*symbol, os);
    os << "GENUS_Enumerate: the seed built from the genus symbol "
       << GenusSymbolString(*symbol) << " is\n";
    WriteMatrix(os, Seed);
    ListSeed.push_back(UniversalMatrixConversion<T, mpq_class>(Seed));
  } else {
    ListSeed = ReadListMatrixFile<T>(FileLattice);
  }
  if (ListSeed.empty()) {
    std::cerr << "GENUS_Enumerate: no seed lattice in " << FileLattice << "\n";
    throw TerminalException{1};
  }
  // A sum of 1 / |Aut(L)|, so a fraction over the type the groups count
  // with, which has nothing to do with the type of the Gram matrices.
  using Tmass = typename GenusEnumerationResult<T, Tgroup>::Tmass;
  // The seeds of one genus have the same rank and the same determinant.
  int rank = ListSeed[0].rows();
  T det = DeterminantMat(ListSeed[0]);
  for (auto &eSeed : ListSeed) {
    if (!IsSymmetricMatrix(eSeed) || !IsPositiveDefinite(eSeed, std::cerr)) {
      std::cerr << "GENUS_Enumerate: a seed Gram matrix in " << FileLattice
                << " is not symmetric positive definite\n";
      throw TerminalException{1};
    }
    // The seeds are integral. Even seeds are enumerated directly; odd ones
    // through their even sublattices (GenusEnumerationOdd.h).
    if (!IsIntegralMatrix(eSeed)) {
      std::cerr << "GENUS_Enumerate: a seed Gram matrix in " << FileLattice
                << " is not integral\n";
      throw TerminalException{1};
    }
    if (eSeed.rows() != rank) {
      std::cerr << "GENUS_Enumerate: the seeds in " << FileLattice
                << " have the ranks " << rank << " and " << eSeed.rows()
                << "\n";
      throw TerminalException{1};
    }
    if (DeterminantMat(eSeed) != det) {
      std::cerr << "GENUS_Enumerate: the seeds in " << FileLattice
                << " have the determinants " << det << " and "
                << DeterminantMat(eSeed) << "\n";
      throw TerminalException{1};
    }
  }
  if (symbol && FileLattice != "unset") {
    for (auto &eSeed : ListSeed) {
      GenusSymbol s_seed = GenusSymbolOfGram<mpq_class, mpz_class>(
          UniversalMatrixConversion<mpq_class, T>(eSeed));
      if (s_seed != *symbol) {
        std::cerr << "GENUS_Enumerate: a seed in " << FileLattice
                  << " has the genus " << GenusSymbolString(s_seed)
                  << ", not GenusSymbol = " << GenusSymbolString(*symbol)
                  << "\n";
        throw TerminalException{1};
      }
    }
  }
  // Even or odd: all the seeds of one genus have the same parity.
  auto is_even = [](MyMatrix<T> const &M) -> bool {
    for (int i = 0; i < M.rows(); i++) {
      if (!IsInteger(M(i, i) / T(2))) {
        return false;
      }
    }
    return true;
  };
  bool even = is_even(ListSeed[0]);
  for (auto &eSeed : ListSeed) {
    if (is_even(eSeed) != even) {
      std::cerr << "GENUS_Enumerate: the seeds in " << FileLattice
                << " are not all even or all odd\n";
      throw TerminalException{1};
    }
  }
  if (!even) {
    os << "GENUS_Enumerate: the seeds are odd, the genus is enumerated "
          "through the genus of their even sublattices\n";
  }
  // The mass, computed from the first seed by the mass formula of Conway
  // and Sloane (GenusMass.h); the enumeration stops when sum 1/|Aut|
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
  int prime = spec.prime;
  if (spec.strategy == NeighborStrategy::OutPrime && prime == 0) {
    // For odd seeds the enumerated genus is that of the even sublattices,
    // of determinant 4 det.
    T det_enum = even ? det : T(4) * det;
    prime = ChooseNeighborPrime<T>(det_enum);
    // The spinor genera of the enumerated genus (GenusSpinor.h). When there
    // are several, the l-neighbour graph of the genus is connected only for
    // some primes l: the smallest odd one is taken, if there is one, so that
    // a single seed suffices.
    if (rank >= 3) {
      MyMatrix<mpq_class> G_enum = UniversalMatrixConversion<mpq_class, T>(
          even ? ListSeed[0] : EvenSublatticeGram<T, Tint>(ListSeed[0]));
      SpinorGenusData spinor = GenusSpinorData<mpq_class, mpz_class>(
          GenusSymbolOfGram<mpq_class, mpz_class>(G_enum), G_enum);
      int n_spinor = spinor.n_spinor_genera();
      os << "GENUS_Enumerate: " << n_spinor
         << " spinor genera in the enumerated genus\n";
      if (n_spinor > 1) {
        std::optional<int> connected;
        for (int l = 3; l < 1000 && !connected; l += 2) {
          bool is_prime = true;
          for (int f = 3; f * f <= l; f += 2) {
            if (l % f == 0) {
              is_prime = false;
            }
          }
          if (is_prime && ResInt(det_enum, T(l)) != 0 &&
              spinor.neighbor_components<mpz_class>(l) == 1) {
            connected = l;
          }
        }
        if (connected) {
          prime = *connected;
          os << "GENUS_Enumerate: the prime " << prime
             << " connects them, it is used\n";
        } else if (static_cast<int>(ListSeed.size()) < n_spinor) {
          os << "GENUS_Enumerate: no single prime connects them; give one "
                "seed per spinor genus, or use iterate or random_walk\n";
        }
      }
    }
  }
  GenusEnumerationResult<T, Tgroup> result =
      even ? GenusEnumeration<T, Tint, Tgroup>(
                 ListSeed, TotalMass, prime, spec.methods, spec.strategy,
                 spec.rw, spec.hecke, spec.order, max_runtime_second, os)
           : GenusEnumerationOdd<T, Tint, Tgroup>(
                 ListSeed, TotalMass, prime, spec.methods, spec.strategy,
                 spec.rw, spec.hecke, spec.order, max_runtime_second, os);
  // With the iterate and random_walk strategies there is no single prime; 0
  // marks that in the output record.
  int prime_out = (spec.strategy == NeighborStrategy::OutPrime) ? prime : 0;
  auto f = [&](std::ostream &os_out) -> void {
    WriteGenusEnumerationResult(os_out, result, OutFormat, prime_out);
  };
  FILE_PrintStderrStdoutFile(OutFile, f);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    // The group orders are counted in TintGroup and the mass lives in
    // overlying_field<TintGroup>. Both are unrelated to the arithmetic of
    // the lattice coordinates, so Tgroup is fixed once here while the
    // coordinate types T / Tint are chosen per arithmetic below.
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    //
    FullNamelist eFull = NAMELIST_GetStandard_GENUS_ENUMERATE();
    ArithmeticDispatch dispatch;
    add_rational_integer_arithmetics_by_library(
        dispatch, [&]<typename T, typename Tint>() -> void {
          ProcessGenus<T, Tint, Tgroup>(eFull, std::cerr);
        });
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GENUS_Enumerate [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      std::cerr << "arithmetic values:\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    std::string arith = GetNamelistStringEntry(eFull, "DATA", "arithmetic");
    dispatch.run(arith);
    std::cerr << "Normal termination of GENUS_Enumerate\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_Enumerate\n";
    exit(e.eVal);
  }
  runtime(time);
}

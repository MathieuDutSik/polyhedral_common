// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "GenusEnumeration.h"
#include "SignatureSymmetric.h"
#include <string>
#include <vector>
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void ProcessGenus(std::string const &FileGenus, std::string const &FileLattice,
                  std::string const &FileMass, std::string const &OutFormat,
                  GenusScheme scheme, std::ostream &os) {
  GenusSpec<T> spec = ReadGenusSpecFile<T>(FileGenus);
  std::vector<MyMatrix<T>> ListSeed = ReadListMatrixFile<T>(FileLattice);
  if (ListSeed.empty()) {
    std::cerr << "GENUS_Enumerate: no seed lattice in " << FileLattice << "\n";
    throw TerminalException{1};
  }
  // A sum of 1 / |Aut(L)|, so a fraction over the type the groups count
  // with, which has nothing to do with the type of the Gram matrices.
  using Tmass = typename GenusEnumerationResult<T, Tgroup>::Tmass;
  Tmass TotalMass = ReadMassFile<Tmass>(FileMass);
  for (auto &eSeed : ListSeed) {
    if (!IsSymmetricMatrix(eSeed) || !IsPositiveDefinite(eSeed, std::cerr)) {
      std::cerr << "GENUS_Enumerate: a seed Gram matrix in " << FileLattice
                << " is not symmetric positive definite\n";
      throw TerminalException{1};
    }
    if (eSeed.rows() != spec.rank) {
      std::cerr << "GENUS_Enumerate: a seed lattice has rank " << eSeed.rows()
                << " but the genus has rank " << spec.rank << "\n";
      throw TerminalException{1};
    }
    if (DeterminantMat(eSeed) != spec.det) {
      std::cerr << "GENUS_Enumerate: a seed lattice has determinant "
                << DeterminantMat(eSeed) << " but the genus has determinant "
                << spec.det << "\n";
      throw TerminalException{1};
    }
  }
  int prime = spec.prime;
  if (prime == 0) {
    prime = ChooseNeighborPrime<T>(spec.det);
  }
  GenusEnumerationResult<T, Tgroup> result =
      GenusEnumeration<T, Tint, Tgroup>(ListSeed, TotalMass, prime, scheme,
                                        std::cerr);
  WriteGenusEnumerationResult(os, result, OutFormat, prime);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 5 && argc != 7 && argc != 8) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GENUS_Enumerate [arith] [Genus] [Lattice] [Mass]\n";
      std::cerr << "    or\n";
      std::cerr << "GENUS_Enumerate [arith] [Genus] [Lattice] [Mass] "
                << "[OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "GENUS_Enumerate [arith] [Genus] [Lattice] [Mass] "
                << "[OutFormat] [OutFile] [Strategy]\n";
      std::cerr << "\n";
      std::cerr << "Strategy (optional): canonical (default), schemeA / ps "
                << "(Plesken-Souvignier isomorphism), or schemeB / matrixgroup "
                << "(MatrixGroup isomorphism)\n";
      std::cerr << "\n";
      std::cerr << "Genus   (input) : the description of the genus, as the\n";
      std::cerr << "                  key/value lines\n";
      std::cerr << "                     rank <n>\n";
      std::cerr << "                     det <determinant>\n";
      std::cerr << "                     prime <p>     (optional, 0 = auto)\n";
      std::cerr << "Lattice (input) : one or more Gram matrices, in the\n";
      std::cerr << "                  ListMatrix format. One seed per SPINOR\n";
      std::cerr << "                  genus is needed: the p-neighbour graph\n";
      std::cerr << "                  is connected on a spinor genus, not on\n";
      std::cerr << "                  a genus, so a single seed may reach\n";
      std::cerr << "                  only part of the genus\n";
      std::cerr << "Mass    (input) : the mass of the genus, as \"num den\"\n";
      std::cerr << "                  or a single integer. The enumeration\n";
      std::cerr << "                  stops when sum 1/|Aut| reaches it, and\n";
      std::cerr << "                  that equality certifies completeness\n";
      std::cerr << "OutFormat values:\n";
      std::cerr << "  GAP     : a record with the Gram matrices, the group\n";
      std::cerr << "            orders and the mass (default)\n";
      std::cerr << "  CPP     : the list of Gram matrices\n";
      std::cerr << "  Summary : the class number and the mass check only\n";
      return -1;
    }
    // The group orders are counted in TintGroup and the mass lives in
    // overlying_field<TintGroup>. Both are unrelated to the arithmetic of
    // the lattice coordinates, so Tgroup is fixed once here while the
    // coordinate types T / Tint are chosen per arithmetic below.
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    //
    std::string arith = argv[1];
    std::string FileGenus = argv[2];
    std::string FileLattice = argv[3];
    std::string FileMass = argv[4];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc >= 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    GenusScheme scheme = GenusScheme::CanonicalForm;
    if (argc == 8) {
      std::string s = argv[7];
      if (s == "schemeA" || s == "ps" || s == "pleskensouvignier" ||
          s == "isomorphism" || s == "iso") {
        scheme = GenusScheme::PleskenSouvignier;
      } else if (s == "schemeB" || s == "matrixgroup" || s == "matgrp") {
        scheme = GenusScheme::MatrixGroup;
      } else if (s != "canonical") {
        std::cerr << "Strategy must be 'canonical', 'schemeA'/'ps', or "
                  << "'schemeB'/'matrixgroup', got " << s << "\n";
        throw TerminalException{1};
      }
    }
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
        using T = boost::multiprecision::mpz_int;
        using Tint = boost::multiprecision::mpz_int;
#else
        using T = mpz_class;
        using Tint = mpz_class;
#endif
        return ProcessGenus<T, Tint, Tgroup>(FileGenus, FileLattice, FileMass,
                                             OutFormat, scheme, os);
      }
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpz_class;
        using Tint = fmpz_class;
        return ProcessGenus<T, Tint, Tgroup>(FileGenus, FileLattice, FileMass,
                                             OutFormat, scheme, os);
      }
#endif
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Allowed values: gmp, flint\n";
#else
      std::cerr << "Allowed values: gmp (build with ENABLE_FLINT_SUPPORT=1 "
                << "for flint)\n";
#endif
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of GENUS_Enumerate\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_Enumerate\n";
    exit(e.eVal);
  }
  runtime(time);
}

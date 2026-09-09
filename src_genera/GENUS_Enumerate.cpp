// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
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
                  std::ostream &os) {
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
      GenusEnumeration<T, Tint, Tgroup>(ListSeed, TotalMass, prime, std::cerr);
  WriteGenusEnumerationResult(os, result, OutFormat, prime);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 5 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GENUS_Enumerate [arith] [Genus] [Lattice] [Mass]\n";
      std::cerr << "    or\n";
      std::cerr << "GENUS_Enumerate [arith] [Genus] [Lattice] [Mass] "
                << "[OutFormat] [OutFile]\n";
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
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using T = boost::multiprecision::mpz_int;
    using Tint = boost::multiprecision::mpz_int;
#else
    using T = mpz_class;
    using Tint = mpz_class;
#endif
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
    if (argc == 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        return ProcessGenus<T, Tint, Tgroup>(FileGenus, FileLattice, FileMass,
                                             OutFormat, os);
      }
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
      std::cerr << "Allowed values: gmp\n";
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

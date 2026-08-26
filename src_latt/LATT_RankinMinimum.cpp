// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "Enumeration_k_space.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename Tgroup, typename T, typename Tint>
void compute_rankin_minimum(int const &k, std::string const &mode,
                            std::string const &FileGram,
                            std::string const &OutFormat, std::ostream &os) {
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileGram);
  if (mode == "all") {
    T tol(0);
    ResultKRankinMin<T, Tint> result =
        Rankin_k_minimum<T, Tint>(GramMat, k, tol, std::cerr);
    if (OutFormat == "GAP") {
      os << "return rec(min:=" << result.min << ", ListSubspace:=";
      WriteListMatrixGAP(os, result.l_spaces);
      os << ");\n";
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
              << "\n";
    throw TerminalException{1};
  }
  if (mode == "orbit") {
    ResultKRankinMinOrbits<T, Tint> result =
        Rankin_k_minimum_orbits<T, Tint, Tgroup>(GramMat, k, std::cerr);
    std::vector<MyMatrix<Tint>> l_repr;
    for (auto &eOrbit : result.l_orbits) {
      l_repr.push_back(eOrbit.representative);
    }
    if (OutFormat == "GAP") {
      os << "return rec(min:=" << result.min << ", ListSubspace:=";
      WriteListMatrixGAP(os, l_repr);
      os << ", ListOrbitSize:=[";
      bool IsFirst = true;
      for (auto &eOrbit : result.l_orbits) {
        if (!IsFirst)
          os << ",";
        IsFirst = false;
        os << eOrbit.orbit_size;
      }
      os << "]);\n";
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
              << "\n";
    throw TerminalException{1};
  }
  std::cerr << "Failed to find a matching entry for mode=" << mode << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 5 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_RankinMinimum [arith] [k] [mode] [FileGram]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_RankinMinimum [arith] [k] [mode] [FileGram] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "arith: gmp, gmp_boost, multi_boost\n";
      std::cerr << "k: The dimension of the sublattices\n";
      std::cerr << "mode: all (all the sublattices attaining the minimum)\n";
      std::cerr << "      or orbit (their orbits under the automorphism "
                   "group)\n";
      std::cerr << "FileGram: The Gram matrix on input\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string strK = argv[2];
    std::string mode = argv[3];
    std::string FileGram = argv[4];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    int k = ParseScalar<int>(strK);
    //
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint_grp = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint_grp>;
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return compute_rankin_minimum<Tgroup, T, Tint>(k, mode, FileGram,
                                                       OutFormat, os);
      }
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return compute_rankin_minimum<Tgroup, T, Tint>(k, mode, FileGram,
                                                       OutFormat, os);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return compute_rankin_minimum<Tgroup, T, Tint>(k, mode, FileGram,
                                                       OutFormat, os);
      }
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of LATT_RankinMinimum\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_RankinMinimum\n";
    exit(e.eVal);
  }
  runtime(time);
}

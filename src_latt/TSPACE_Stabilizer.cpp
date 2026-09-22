// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Tspace_General.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T>
void write_group(std::vector<MyMatrix<T>> const &LGen,
                 std::string const &OutFormat, std::ostream &os) {
  if (OutFormat == "count") {
    os << "number of generators=" << LGen.size() << "\n";
    return;
  }
  if (OutFormat == "GAP") {
    os << "return Group([";
    bool IsFirst = true;
    for (auto &eGen : LGen) {
      if (!IsFirst) {
        os << ",\n";
      }
      os << StringMatrixGAP_line(eGen);
      IsFirst = false;
    }
    os << "]);\n";
    return;
  }
  std::cerr << "Failed to find a matching format\n";
  throw TerminalException{1};
}

template <typename T, typename Tint>
void ComputeStabilizer(std::string const &FileTspace,
                       std::string const &FileGram,
                       std::string const &OutFormat, std::ostream &os) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  LinSpaceMatrix<T> LinSpa = ReadLinSpaceFile<T>(FileTspace, std::cerr);
  MyMatrix<T> eMat = ReadMatrixFile<T>(FileGram);
  std::optional<MyMatrix<T>> CommonGramMat;
  std::vector<MyMatrix<T>> ListGen =
      LINSPA_ComputeStabilizer<T, Tint, Tgroup>(LinSpa, eMat, CommonGramMat,
                                                std::cerr);
  write_group(ListGen, OutFormat, os);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TSPACE_Stabilizer [arith] [FileTspace] [FileGram]\n";
      std::cerr << "or\n";
      std::cerr << "TSPACE_Stabilizer [arith] [FileTspace] [FileGram] "
                   "[OutFormat] [FileOut]\n";
      std::cerr << "\n";
      std::cerr << "arith values:\n";
      std::cerr << "  gmp         : mpq_class / mpz_class (default choice)\n";
      std::cerr << "  gmp_boost   : the boost bindings to the gmp types\n";
      std::cerr << "  multi_boost : the boost multiprecision types\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "  flint       : fmpq_class / fmpz_class\n";
#endif
      std::cerr << "OutFormat values:\n";
      std::cerr << "  count : the number of generators only (default)\n";
      std::cerr << "  GAP   : the group generators, GAP readable\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileTspace = argv[2];
    std::string FileGram = argv[3];
    std::string OutFormat = "count";
    std::string FileOut = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      FileOut = argv[5];
    }
    //
    auto f = [&](std::ostream &os_out) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return ComputeStabilizer<T, Tint>(FileTspace, FileGram, OutFormat,
                                          os_out);
      }
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return ComputeStabilizer<T, Tint>(FileTspace, FileGram, OutFormat,
                                          os_out);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return ComputeStabilizer<T, Tint>(FileTspace, FileGram, OutFormat,
                                          os_out);
      }
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        using Tint = fmpz_class;
        return ComputeStabilizer<T, Tint>(FileTspace, FileGram, OutFormat,
                                          os_out);
      }
#endif
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Available possibilities: gmp, gmp_boost, "
                << "multi_boost, flint\n";
#else
      std::cerr << "Available possibilities: gmp, gmp_boost, "
                << "multi_boost (build with ENABLE_FLINT_SUPPORT=1 "
                << "for flint)\n";
#endif
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(FileOut, f);
    std::cerr << "Normal termination of TSPACE_Stabilizer\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TSPACE_Stabilizer\n";
    exit(e.eVal);
  }
  runtime(time);
}

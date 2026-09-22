// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "Tspace_General.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T>
void write_result(std::optional<MyMatrix<T>> const &opt,
                  std::string const &OutFormat, std::ostream &os) {
  if (OutFormat == "simple") {
    if (opt) {
      os << "It is equivalence\n";
    } else {
      os << "It is not equivalence\n";
    }
    return;
  }
  if (OutFormat == "GAP") {
    if (opt) {
      os << "return " << StringMatrixGAP_line(*opt) << ";\n";
    } else {
      os << "return fail;\n";
    }
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

template <typename T, typename Tint>
void ComputeEquivalence(std::string const &FileTspace,
                        std::string const &FileGram1,
                        std::string const &FileGram2,
                        std::string const &OutFormat, std::ostream &os) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  LinSpaceMatrix<T> LinSpa = ReadLinSpaceFile<T>(FileTspace, std::cerr);
  MyMatrix<T> eMat1 = ReadMatrixFile<T>(FileGram1);
  MyMatrix<T> eMat2 = ReadMatrixFile<T>(FileGram2);
  std::optional<MyMatrix<Tint>> opt =
      LINSPA_TestEquivalenceGramMatrix<T, Tint, Tgroup>(LinSpa, eMat1, eMat2,
                                                        std::cerr);
  write_result(opt, OutFormat, os);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 5 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TSPACE_Equivalence [arith] [FileLinSpa] [FileMat1] "
                   "[FileMat2]\n";
      std::cerr << "or\n";
      std::cerr << "TSPACE_Equivalence [arith] [FileLinSpa] [FileMat1] "
                   "[FileMat2] [OutFormat] [FileOut]\n";
      std::cerr << "\n";
      std::cerr << "arith values:\n";
      std::cerr << "  gmp         : mpq_class / mpz_class (default choice)\n";
      std::cerr << "  gmp_boost   : the boost bindings to the gmp types\n";
      std::cerr << "  multi_boost : the boost multiprecision types\n";
      std::cerr << "OutFormat values:\n";
      std::cerr << "  simple : a human readable one line answer (default)\n";
      std::cerr << "  GAP    : the equivalence matrix or fail, GAP readable\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileTspace = argv[2];
    std::string FileGram1 = argv[3];
    std::string FileGram2 = argv[4];
    std::string OutFormat = "simple";
    std::string FileOut = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      FileOut = argv[6];
    }
    //
    auto f = [&](std::ostream &os_out) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return ComputeEquivalence<T, Tint>(FileTspace, FileGram1, FileGram2,
                                           OutFormat, os_out);
      }
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return ComputeEquivalence<T, Tint>(FileTspace, FileGram1, FileGram2,
                                           OutFormat, os_out);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return ComputeEquivalence<T, Tint>(FileTspace, FileGram1, FileGram2,
                                           OutFormat, os_out);
      }
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: gmp, gmp_boost, multi_boost\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(FileOut, f);
    std::cerr << "Normal termination of TSPACE_Equivalence\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TSPACE_Equivalence\n";
    exit(e.eVal);
  }
  runtime(time);
}

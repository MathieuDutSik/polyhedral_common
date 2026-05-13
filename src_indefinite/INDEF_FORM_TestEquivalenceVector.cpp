// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "CombinedAlgorithms.h"
#include "Group.h"
#include "Permutation.h"
// clang-format on

template <typename T, typename Tint>
void process(std::string const &FileM, std::string const &FileV1,
             std::string const &FileV2, std::string const &OutFormat,
             std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = Tint;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> Q = ReadMatrixFile<T>(FileM);
  MyMatrix<Tint> v1 = ReadVectorFile<Tint>(FileV1);
  MyMatrix<Tint> v2 = ReadVectorFile<Tint>(FileV2);
  IndefiniteCombinedAlgo<T, Tint, Tgroup> comb(std::cerr);
  std::optional<MyMatrix<Tint>> opt =
      comb.INDEF_FORM_EquivalenceVector(Q, Q, v1, v2);
  if (OutFormat == "PYTHON") {
    if (opt) {
      WriteMatrixPYTHON(os_out, *opt);
    } else {
      os_out << "None";
    }
    os_out << "\n";
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return ";
    if (opt) {
      WriteMatrixGAP(os_out, *opt);
    } else {
      os_out << "fail";
    }
    os_out << ";\n";
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 5 && argc != 7) {
      std::cerr << "INDEF_FORM_TestEquivalenceVector [arith] [FileM] [FileV1] "
                   "[FileV2]\n";
      std::cerr << "or\n";
      std::cerr << "INDEF_FORM_TestEquivalenceVector [arith] [FileM] [FileV1] "
                   "[FileV2] [OutFormat] [OutFile]\n";
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string FileM = argv[2];
    std::string FileV1 = argv[3];
    std::string FileV2 = argv[4];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint>(FileM, FileV1, FileV2, OutFormat, os);
      }
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return process<T, Tint>(FileM, FileV1, FileV2, OutFormat, os);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return process<T, Tint>(FileM, FileV1, FileV2, OutFormat, os);
      }
      std::cerr << "Failed to find matching type for arith\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(OutFile, f);
    //
    std::cerr << "Normal termination of INDEF_FORM_TestEquivalence\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in INDEF_FORM_TestEquivalence, runtime=" << time
              << "\n";
    exit(e.eVal);
  }
  runtime(time);
}

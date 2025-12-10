// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "NumberTheorySafeInt.h"
#include "DeterminantMinimization.h"
// clang-format on

template <typename T>
void process(std::string const &FileI, std::string const &OutFormat,
             std::ostream &os) {
  MyMatrix<T> Q1 = ReadMatrixFile<T>(FileI);
  //
  auto f_print = [&](MyMatrix<T> const &Qinp) -> void {
    T det1 = DeterminantMat(Qinp);
    ResultDetMin<T> res = DeterminantMinimization(Qinp, true, std::cerr);
    T det2 = DeterminantMat(res.Mred);
    std::cerr << "det1=" << det1 << " det2=" << det2 << "\n";
    if (OutFormat == "GAP") {
      os << "return rec(P:=";
      WriteMatrixGAP(os, res.P);
      os << ", Mred:=";
      WriteMatrixGAP(os, res.Mred);
      os << ");\n";
      return;
    }
    if (OutFormat == "CPP") {
      WriteMatrix(os, res.Mred);
      return;
    }
    std::cerr << "Failed to find a matching OutFormat\n";
    throw TerminalException{1};
  };
  MyMatrix<T> Qinv = Inverse(Q1);
  MyMatrix<T> Q2 = RemoveFractionMatrix(Qinv);
  f_print(Q1);
  f_print(Q2);
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "LATT_Automorphyreduction arithmetic [FileI] [OutFormat] "
                   "[FileO]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_AutomorphyReduction arithmetic [FileI]\n";
      std::cerr << "or\n";
      std::cerr << "Determining if two matrices are equivalent and computing "
                   "the automorphism group\n";
      std::cerr << "are tough problems. Two tools may help:\n";
      std::cerr << "* The equivalence for the inverse may be easier\n";
      std::cerr << "* The determinant may be reduceable and the problem easier "
                   "to solve for the reduced form\n";
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileO = argv[4];
    }

    //
    auto f = [&](std::ostream &os) -> void {
      if (arith == "rational") {
        using T = mpq_class;
        return process<T>(FileI, OutFormat, os);
      }
      /*
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(FileI, OutFormat, os);
      }
      */

      /*
      if (arith == "cpp_rational") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(FileI, OutFormat, os);
      }
      */

      /*
      if (arith == "mpq_rational") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(FileI, OutFormat, os);
      }
      */

      std::cerr << "Failed to find matching type for arith. Possibilities: "
                   "rational\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(FileO, f);
    std::cerr << "Normal termination of LATT_AutomorphyReduction\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_AutomorphyReduction\n";
    exit(e.eVal);
  }
  runtime(time);
}

// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "NumberTheorySafeInt.h"
#include "Isotropic.h"
// clang-format on

template <typename T>
void process(std::string const &FileI, std::string const &OutFormat,
             std::ostream &os_out) {
  MyMatrix<T> Q = ReadMatrixFile<T>(FileI);
  //
  std::optional<MyVector<T>> opt = FindIsotropic(Q, std::cerr);
  if (OutFormat == "GAP") {
    if (opt) {
      MyVector<T> const &eV = *opt;
      T val = EvaluationQuadForm(Q, eV);
      if (val != 0) {
        std::cerr << "LATT_FindIsotropic: eV is not an isotropic vector\n";
        throw TerminalException{1};
      }
      os_out << "return rec(has_isotropic:=true, ";
      os_out << "V:=" << StringVectorGAP(eV) << ");\n";
    } else {
      os_out << "return rec(has_isotropic:=false);\n";
    }
    return;
  }
  if (OutFormat == "PYTHON") {
    if (opt) {
      MyVector<T> const &eV = *opt;
      WriteVectorPYTHON(os_out, eV);
    } else {
      os_out << "None";
    }
    os_out << "\n";
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "LATT_FindIsotropic [arith] [FileI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_FindIsotropic [arith] [FileI]\n";
      std::cerr << "\n";
      std::cerr << "Possibilities for [arith]: rational\n";
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
      std::cerr << "Failed to find matching type for arith\n";
      throw TerminalException{1};
    };
    if (FileO == "stderr") {
      f(std::cerr);
    } else {
      if (FileO == "stdout") {
        f(std::cout);
      } else {
        std::ofstream os(FileO);
        f(os);
      }
    }
    std::cerr << "Normal termination of LATT_FindIsotropic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_FindIsotropic\n";
    exit(e.eVal);
  }
  runtime(time);
}

// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "NumberTheorySafeInt.h"
#include "Indefinite_LLL.h"
// clang-format on

template <typename T, typename Tint>
void process(std::string const &FileI, std::string const &OutFormat,
             std::ostream &os) {
  MyMatrix<T> M = ReadMatrixFile<T>(FileI);
  std::cerr << "We have M\n";

  ResultIndefiniteLLL<T, Tint> res = Indefinite_LLL<T, Tint>(M);
  std::cerr << "B_T=\n";
  WriteMatrix(std::cerr, res.B);
  std::cerr << "Mred=\n";
  WriteMatrix(std::cerr, res.Mred);
  if (OutFormat == "GAP") {
    os << "return rec(B:=" << StringMatrixGAP(res.B);
    os << ", Mred:=" << StringMatrixGAP(res.Mred);
    if (res.Xisotrop) {
      MyVector<T> const& Xisotrop = *res.Xisotrop;
      os << ", Xisotrop:=" << StringVectorGAP(Xisotrop);
    }
    os << ");\n";
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 3 && argc != 2) {
      std::cerr
          << "LATT_IndefiniteLLL arithmetic [FileI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_IndefiniteLLL arithmetic [FileI]\n";
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
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint>(FileI, OutFormat, os);
      }
      if (arith == "safe") {
        using T = Rational<SafeInt64>;
        using Tint = SafeInt64;
        return process<T, Tint>(FileI, OutFormat, os);
      }
      if (arith == "boost_cpp") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return process<T, Tint>(FileI, OutFormat, os);
      }
      /*
      if (arith == "boost_gmp") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return process<T,Tint>(FileI, OutFormat, os);
      }
      */
      std::cerr << "Failed to find a matching type\n";
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
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in IndefiniteLLL\n";
    exit(e.eVal);
  }
  runtime(time1);
}

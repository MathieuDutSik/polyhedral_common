// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Indefinite_LLL.h"
// clang-format on

template <typename T>
void process(std::string const &FileI, std::string const &OutFormat,
             std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  MyMatrix<T> M = ReadMatrixFile<T>(FileI);

  ResultIndefiniteLLL<T, Tint> res = Indefinite_LLL<T, Tint>(M, std::cerr);
  if (OutFormat == "GAP") {
    os << "return rec(B:=" << StringMatrixGAP(res.B);
    os << ", Mred:=" << StringMatrixGAP(res.Mred);
    if (res.Xisotrop) {
      MyVector<T> const &Xisotrop = *res.Xisotrop;
      os << ", Xisotrop:=" << StringVectorGAP(Xisotrop);
    }
    os << ");\n";
    return;
  }
  if (OutFormat == "CPP") {
    os << "B_T=\n";
    WriteMatrix(os, res.B);
    os << "Mred=\n";
    WriteMatrix(os, res.Mred);
    if (res.Xisotrop) {
      MyVector<T> const &Xisotrop = *res.Xisotrop;
      WriteVector(os, Xisotrop);
    } else {
      os << "no isotrop vector found\n";
    }
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
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
        return process<T>(FileI, OutFormat, os);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(FileI, OutFormat, os);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(FileI, OutFormat, os);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        return process<T>(FileI, OutFormat, os);
      }
#endif
      std::cerr << "Failed to find a matching type for arith\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Possibilities: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << ", flint\n";
#else
      std::cerr << "Possibilities: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << " (build with ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(FileO, f);
    std::cerr << "Normal termination of LATT_IndefiniteLLL\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_IndefiniteLLL\n";
    exit(e.eVal);
  }
  runtime(time);
}

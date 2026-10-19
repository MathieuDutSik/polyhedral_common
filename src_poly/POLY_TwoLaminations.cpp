// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Laminations.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
MyMatrix<T> ReordListPoint(const std::vector<MyVector<T>> &ListPoint) {
  std::set<MyVector<T>> e_set;
  for (auto &ePt : ListPoint)
    e_set.insert(ePt);
  std::vector<MyVector<T>> ListPt;
  for (auto &ePt : e_set)
    ListPt.push_back(ePt);
  return MatrixFromVectorFamily(ListPt);
}

template <typename T>
void process(std::string const &opt, std::string const &FileM,
             std::string const &OutFormat, std::ostream &os) {
  MyMatrix<T> M = ReadMatrixFile<T>(FileM);
  if (opt == "all") {
    vectface vf = compute_all_two_laminations(M);
    if (OutFormat == "GAP") {
      os << "return ";
      WriteListFaceGAP(os, vf);
      os << ";\n";
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
              << "\n";
    throw TerminalException{1};
  }
  if (opt == "one") {
    std::optional<Face> opt = compute_one_two_laminations(M);
    if (OutFormat == "GAP") {
      os << "return ";
      if (opt) {
        WriteFaceGAP(os, *opt);
      } else {
        os << "fail";
      }
      os << ";\n";
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
              << "\n";
    throw TerminalException{1};
  }
  std::cerr << "Failed to find a matching entry for opt=" << opt << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string opt, FileM;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(opt, FileM, OutFormat, *os_ptr);
    };
    add_rational_arithmetics(dispatch, f_arith);
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "POLY_TwoLaminations [arith] [opt] [FileM]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_TwoLaminations [arith] [opt] [FileM] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "\n";
      std::cerr << "arith: " << dispatch.choices() << "\n";
      std::cerr << "opt: one or all\n";
      std::cerr << "FileM: File of the matrix\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    opt = argv[2];
    FileM = argv[3];
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    //
    auto f = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    //
    std::cerr << "Normal termination of POLY_TwoLaminations\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_TwoLaminations\n";
    exit(e.eVal);
  }
  runtime(time1);
}

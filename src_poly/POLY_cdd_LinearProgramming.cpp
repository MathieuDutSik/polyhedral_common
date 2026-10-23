// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "NumberTheoryRealField.h"
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "POLY_LinearProgramming.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const &eFileFAC, std::string const &eFileIneq,
             std::string const &OutFormat, std::ostream &os) {
  MyMatrix<T> TheEXT = ReadMatrixFile<T>(eFileFAC);
  MyVector<T> eVect = ReadVectorFile<T>(eFileIneq);
  if (TheEXT.cols() != eVect.size()) {
    std::cerr << "|TheEXT|=" << TheEXT.rows() << " / " << TheEXT.cols() << "\n";
    std::cerr << "|eVect|=" << eVect.size() << "\n";
    std::cerr << "This is inconsistentt\n";
    throw TerminalException{1};
  }
  LpSolution<T> eSol = SIMPLEX_LinearProgramming(TheEXT, eVect, std::cerr);
  if (OutFormat == "GAP") {
    os << "return rec(";
    os << "OptimalValue:=" << eSol.OptimalValue;
    if (eSol.DirectSolution) {
      os << ",\n primal_solution:=" << StringVectorGAP(*eSol.DirectSolution);
    }
    if (eSol.DualSolution) {
      os << ",\n dual_solution:=" << StringVectorGAP(*eSol.DualSolution);
    }
    if (eSol.DirectSolution && eSol.DualSolution) {
      os << ", face:=";
      Face eFace = ComputeFaceLpSolution(TheEXT, eSol);
      WriteFaceGAP(os, eFace);
    }
    os << ");\n";
    return;
  }
  if (OutFormat == "Oscar") {
    os << eSol.OptimalValue << "\n";
    if (eSol.DirectSolution) {
      WriteVector(os, *eSol.DirectSolution);
    } else {
      os << "0\n";
    }
    if (eSol.DualSolution) {
      WriteVector(os, *eSol.DualSolution);
    } else {
      os << "0\n";
    }
    return;
  }
  std::cerr << "Failed to find a matching for OutFormat=" << OutFormat << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using Trat = boost::multiprecision::mpq_rational;
#else
    using Trat = mpq_class;
#endif
    std::string eFileFAC;
    std::string eFileIneq;
    std::string OutFormat = "GAP";
    std::string eFileO = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(eFileFAC, eFileIneq, OutFormat, *os_ptr);
    };
    add_rational_arithmetics<Trat>(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "POLY_cdd_LinearProgramming arith [FileFAC] [FileIneq] "
             "[OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_cdd_LinearProgramming arith [FileFAC] [FileIneq]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "arith     : the chosen arithmetic\n";
      std::cerr << "FileFAC   : The list of inequalities\n";
      std::cerr << "FileIneq  : The inequality to be minimized\n";
      std::cerr << "OutFormat : The formatting of the output, GAP or Oscar\n";
      std::cerr << "DATAOUT   : The file of output (if present, otherwise "
                   "std::cerr)\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    //
    std::string arith = argv[1];
    eFileFAC = argv[2];
    eFileIneq = argv[3];
    if (argc == 6) {
      OutFormat = argv[4];
      eFileO = argv[5];
    }
    auto call_lp = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(eFileO, call_lp);
    std::cerr << "Normal termination of POLY_cdd_LinearProgramming\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_cdd_LinearProgramming\n";
    exit(e.eVal);
  }
  runtime(time1);
}

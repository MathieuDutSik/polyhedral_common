// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "POLY_RedundancyElimination.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process_A(std::string const &eFileI, std::string const &eFileO,
               std::string const &method, std::string const &OutFormat,
               std::ostream &os) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileI);
  auto get_list_irred = [&]() -> std::vector<int> {
    if (method == "Clarkson") {
      return SIMPLEX_RedundancyReductionClarksonExt(EXT, os);
    }
    if (method == "HitAndRun") {
      return EliminationByRedundance_HitAndRun(EXT, os);
    }
    std::cerr << "Failed to find a matching method\n";
    std::cerr << "Allowed ones: Clarkson and HitAndRun\n";
    throw TerminalException{1};
  };
  std::vector<int> ListIrred = get_list_irred();
  int nbIrred = ListIrred.size();
  std::cerr << "nbIrred=" << nbIrred << "\n";
  auto print_result = [&](std::ostream &os_out) -> void {
    if (OutFormat == "GAP") {
      os_out << "return [";
      for (int i = 0; i < nbIrred; i++) {
        if (i > 0)
          os_out << ",";
        int eVal = ListIrred[i] + 1;
        os_out << eVal;
      }
      os_out << "];\n";
      return;
    }
    if (OutFormat == "Python") {
      for (int i = 0; i < nbIrred; i++) {
        if (i > 0)
          os_out << " ";
        int eVal = ListIrred[i];
        os_out << eVal;
      }
      os_out << "\n";
      return;
    }
    std::cerr << "Failed to find a matching entry\n";
    throw TerminalException{1};
  };
  FILE_PrintStderrStdoutFile(eFileO, print_result);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string method, eFileI;
    std::string OutFormat = "GAP";
    std::string eFileO = "stderr";
    ArithmeticDispatch dispatch;
    auto f = [&]<typename T>() -> void {
      process_A<T>(eFileI, eFileO, method, OutFormat, std::cerr);
    };
    add_rational_arithmetics(dispatch, f);
    add_algebraic_field_arithmetics(dispatch, f);
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "POLY_redundancy method arith [DATAIN] [OutFormat] [DATAOUT]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_redundancy method arith [DATAIN]\n";
      std::cerr << "\n";
      std::cerr << "choice : the choice for presenting the output\n";
      std::cerr << "arith : the chosen arithmetic for the output\n";
      std::cerr << "DATAIN : The polyhedral cone inequalities\n";
      std::cerr << "DATAOUT : The list of irredundant facets\n";
      std::cerr << "\n";
      std::cerr << "     ---- method ----\n";
      std::cerr << "\n";
      std::cerr << "Clarkson  : For the Clarkson method\n";
      std::cerr << "HitAndRun : For the hit and run method\n";
      std::cerr << "\n";
      std::cerr << "     ---- OutFormat ----\n";
      std::cerr << "\n";
      std::cerr << "GAP : For having a gap readable file\n";
      std::cerr << "Python : For having a python readable file\n";
      std::cerr << "\n";
      std::cerr << "     ---- arith ----\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    method = argv[1];
    std::string arith = argv[2];
    eFileI = argv[3];
    if (argc == 6) {
      OutFormat = argv[4];
      eFileO = argv[5];
    }
    dispatch.run(arith);
    std::cerr << "Normal termination of POLY_redundancy\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_redundancy\n";
    exit(e.eVal);
  }
  runtime(time1);
}

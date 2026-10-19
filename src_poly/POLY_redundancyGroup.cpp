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
#include "GRP_GroupFct.h"
#include "GRP_GroupFile.h"
#include "Group.h"
#include "POLY_RedundancyElimination.h"
#include "Permutation.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T, typename Tgroup>
void process(std::string const &FileEXT, std::string const &FileGRP,
             std::string const &method, std::string const &OutFormat,
             std::ostream &os_out, std::ostream &os) {
  std::cerr << "process: start\n";
  MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
  std::cerr << "process: We have EXT\n";
  size_t nbRow = EXT.rows();
  Tgroup GRP = ReadGroupFile<Tgroup>(FileGRP);
  std::cerr << "|GRP|=" << GRP.size() << " nbRow=" << nbRow << "\n";
  if (GRP.n_act() != nbRow) {
    std::cerr << "The group acts on n_act=" << GRP.n_act()
              << " points but the matrix has nbRow=" << nbRow << " rows\n";
    std::cerr << "The group file must be for the action on the rows\n";
    throw TerminalException{1};
  }
  auto get_list_irred = [&]() -> std::vector<int> {
    if (method == "ClarksonBlock") {
      vectface vfo = DecomposeOrbitPoint_Full(GRP);
      size_t n_orbit = vfo.size();
      std::vector<int> BlockBelong(nbRow);
      for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
        Face f = vfo[i_orbit];
        for (size_t i = 0; i < nbRow; i++) {
          if (f[i] == 1) {
            BlockBelong[i] = i_orbit;
          }
        }
      }
      return SIMPLEX_RedundancyReductionClarksonBlocksExt(EXT, BlockBelong, os);
    }
    if (method == "Equivariant") {
      std::cerr << "process: before GetNonRedundant_Equivariant\n";
      return GetNonRedundant_Equivariant(EXT, GRP, os);
    }
    std::cerr << "Failed to find a relevant method\n";
    std::cerr << "Allowed ones: ClarksonBlock and Equivariant\n";
    throw TerminalException{1};
  };
  std::vector<int> ListIrred = get_list_irred();
  int nbIrred = ListIrred.size();
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
  if (OutFormat == "MyVector") {
    os_out << nbIrred << "\n";
    for (int i = 0; i < nbIrred; i++) {
      os_out << " ";
      int eVal = ListIrred[i];
      os_out << eVal;
    }
    return;
  }
  std::cerr << "No meaningful choice for OutFormat\n";
  std::cerr << "Allowed options are GAP and MyVector\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string method;
    std::string FileEXT;
    std::string FileGRP;
    std::string OutFormat = "GAP";
    std::string FileOut = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T, Tgroup>(FileEXT, FileGRP, method, OutFormat, *os_ptr,
                                  std::cerr);
    };
    add_rational_arithmetics(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
    if (argc != 5 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_redundancyGroup method arith [FileEXT] [FileGRP] "
                   "[OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_redundancyGroup method arith [FileEXT] [FileGRP]\n";
      std::cerr << "\n";
      std::cerr << "FileEXT   : The polyhedral cone inequalities (or "
                   "generators of vertices/extreme rays)\n";
      std::cerr << "FileGRP   : The group being used\n";
      std::cerr
          << "OutFormat : Format for output, GAP or MyVector are allowed\n";
      std::cerr << "DATAOUT   : The list of irredundant facets (if absent then "
                   "std::cerr)\n";
      std::cerr << "\n";
      std::cerr << "        --- method ---\n";
      std::cerr << "\n";
      std::cerr << "ClarksonBlock : The Clarkson method with orbit-wise "
                   "propagation of the conclusions\n";
      std::cerr << "Equivariant   : For using the equivariant method\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    //
    method = argv[1];
    std::string arith = argv[2];
    FileEXT = argv[3];
    FileGRP = argv[4];
    if (argc == 7) {
      OutFormat = argv[5];
      FileOut = argv[6];
    }
    auto compute_redundancy = [&](std::ostream &os_out) -> void {
      os_ptr = &os_out;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileOut, compute_redundancy);
    std::cerr << "Normal termination of POLY_redundancyGroup\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_redundancyGroup\n";
    exit(e.eVal);
  }
  runtime(time1);
}

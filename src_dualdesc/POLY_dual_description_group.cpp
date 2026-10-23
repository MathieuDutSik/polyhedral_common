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
#include "POLY_Fundamental.h"
#include "POLY_DirectDualDesc.h"
#include "GRP_DoubleCoset.h"
#include "GRP_GroupFct.h"
#include "GRP_GroupFile.h"
#include "Group.h"
#include "Permutation.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T, typename Tgroup>
void process(std::string const &eFileI, std::string const &eFileG,
             std::string const &ansProg, std::string const &input,
             std::string const &OutFormat, std::ostream &os) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileI);
  auto get_vf = [&]() -> vectface {
    if (input == "general") {
      return DirectFacetComputationIncidence_General(EXT, ansProg, os);
    }
    if (input == "vertices") {
      return DirectFacetComputationIncidence(EXT, ansProg, os);
    }
    std::cerr << "POLY_dual_description_group : input=" << input
              << " but allowed possibilities are general and vertices\n";
    throw TerminalException{1};
  };
  vectface vf = get_vf();
  Tgroup GRP = ReadGroupFile<Tgroup>(eFileG);
  std::unordered_set<Face> SetFace(vf.begin(), vf.end());
  vectface vfo = OrbitSplittingSet_T(SetFace, GRP);
  if (OutFormat == "GAP") {
    os << "return ";
    VectVectInt_Gap_Print(os, vfo);
    os << ";\n";
    return;
  }
  if (OutFormat == "Oscar") {
    MyMatrix<int> M = VectfaceAsMatrix(vfo);
    WriteMatrix(os, M);
    return;
  }
  std::cerr << "POLY_dual_description_group : Failed to find a matching entry "
               "for OutFormat="
            << OutFormat << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using Tint = boost::multiprecision::mpz_int;
    using Trat = boost::multiprecision::mpq_rational;
#else
    using Tint = mpz_class;
    using Trat = mpq_class;
#endif
    //
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string command;
    std::string input;
    std::string eFileI;
    std::string eFileG;
    std::string OutFormat = "GAP";
    std::string eFileO = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T, Tgroup>(eFileI, eFileG, command, input, OutFormat, *os_ptr);
    };
    add_rational_arithmetics<Trat>(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
    if (argc != 6 && argc != 8) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_dual_description_group arith command input [FileEXT] "
                   "[FileGRP] [OutFormat] [FileOUT]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_dual_description_group arith command input [FileEXT] "
                   "[FileGRP]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "arith     : the chosen arithmetic\n";
      std::cerr << "command   : the program used for computing the dual "
                   "description\n";
      std::cerr << "input     : vertices if the rows of FileEXT are distinct "
                   "extreme rays,\n";
      std::cerr << "            general otherwise (the incidences are then "
                   "over all the rows)\n";
      std::cerr << "FileEXT   : The polyhedral cone inequalities\n";
      std::cerr << "FileGRP   : The permutation group used for reduction\n";
      std::cerr << "OutFormat : The file format in output, GAP or Oscar\n";
      std::cerr << "FileOUT   : The file of output (if present, otherwise "
                   "std::cerr)\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      std::cerr << "\n";
      std::cerr << "        --- command ---\n";
      std::cerr << "\n";
      std::cerr << "cdd      : the cdd program\n";
      std::cerr << "lrs      : the lrs program\n";
      std::cerr
          << "pd_lrs   : the primal dual technique with lrs used for test\n";
      std::cerr << "normaliz : the external program normaliz\n";
      return -1;
    }
    std::string arith = argv[1];
    command = argv[2];
    input = argv[3];
    eFileI = argv[4];
    eFileG = argv[5];
    if (argc == 8) {
      OutFormat = argv[6];
      eFileO = argv[7];
    }
    auto call_dualdesc = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(eFileO, call_dualdesc);
    std::cerr << "Normal termination of POLY_dual_description_group\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_dual_description_group\n";
    exit(e.eVal);
  }
  runtime(time1);
}

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
#include "POLY_DirectDualDesc.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const &eFileI, std::string const &ansProg,
             std::string const &choice, std::ostream &os_out) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileI);
  if (choice == "control") {
    MyMatrix<T> FAC =
      DirectFacetComputationInequalities(EXT, ansProg, std::cerr);
    vectface vf = DirectFacetComputationIncidence(EXT, ansProg, std::cerr);
    os_out << "Obtained results:\n";
    os_out << "|vf|=" << vf.n << " / " << vf.n_face << "\n";
    size_t pos = 0;
    for (Face f : vf) {
      os_out << "pos=" << pos << " f=" << StringFace(f) << " |f|=" << f.count()
         << "\n";
    }
    return WriteMatrix(os_out, FAC);
  }
  if (choice == "GAP") {
    MyMatrix<T> FAC =
      DirectFacetComputationInequalities(EXT, ansProg, std::cerr);
    os_out << "return ";
    WriteMatrixGAP(os_out, FAC);
    os_out << ";\n";
    return;
  }
  if (choice == "GAPincidence") {
    vectface vf = DirectFacetComputationIncidence(EXT, ansProg, std::cerr);
    os_out << "return ";
    VectVectInt_Gap_Print(os_out, vf);
    os_out << ";\n";
    return;
  }
  if (choice == "CPP") {
    MyMatrix<T> FAC =
      DirectFacetComputationInequalities(EXT, ansProg, std::cerr);
    return WriteMatrix(os_out, FAC);
  }
  if (choice == "Number") {
    MyMatrix<T> FAC =
      DirectFacetComputationInequalities(EXT, ansProg, std::cerr);
    os_out << "|FAC|=" << FAC.rows() << "\n";
    return;
  }
  std::cerr << "We have choice=" << choice << "\n";
  std::cerr << "But allowed possibilities are control, CPP, GAP and Number\n";
  std::cerr << "Failed to find a matching entry\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string command;
    std::string choice;
    std::string eFileI;
    std::string eFileO = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(eFileI, command, choice, *os_ptr);
    };
    dispatch.add("mpz_class", "integer arithmetic over GMP mpz_class",
                 [f_arith]() -> void { f_arith.template operator()<mpz_class>(); });
    add_rational_arithmetics(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
    if (argc != 5 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "POLY_dual_description arith command choice [DATAIN] [DATAOUT]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_dual_description arith command choice [DATAIN]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "arith   : the chosen arithmetic\n";
      std::cerr
          << "command : the program used for computing the dual description\n";
      std::cerr << "DATAIN  : The polyhedral cone inequalities\n";
      std::cerr
          << "DATAOUT : The file of output (if present, otherwise std::cerr)\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      std::cerr << "\n";
      std::cerr << "        --- command ---\n";
      std::cerr << "\n";
      std::cerr << "cdd      : the cdd program\n";
      std::cerr << "lrs      : the lrs program\n";
      std::cerr << "pd_lrs   : the pd with lrs used for checks\n";
      std::cerr << "bb       : the beneath-and-beyond method\n";
      std::cerr << "normaliz : the native port of the Normaliz primal "
                   "algorithm\n";
      std::cerr << "\n";
      std::cerr << "        --- choice ---\n";
      std::cerr << "\n";
      std::cerr << "control : the full data set for control and debugging\n";
      std::cerr << "CPP     : the matrix for output (also used in Oscar)\n";
      std::cerr << "Number  : the number of rays found\n";
      std::cerr << "GAP     : returns a GAP readable file (except for "
                   "algebraic field)\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    command = argv[2];
    choice = argv[3];
    eFileI = argv[4];
    if (argc == 6)
      eFileO = argv[5];
    auto dual_desc = [&](std::ostream &os_out) -> void {
      os_ptr = &os_out;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(eFileO, dual_desc);
    std::cerr << "Normal termination of POLY_dual_description\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_dual_description\n";
    exit(e.eVal);
  }
  runtime(time1);
}

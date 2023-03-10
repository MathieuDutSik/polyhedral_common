// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "QuadField.h"
#include "POLY_DirectDualDesc.h"
// clang-format on

template <typename T>
void process(std::string const &eFileI, std::string const& ansProg, std::ostream &os) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileI);
  vectface vf = DirectFacetComputationIncidence(EXT, ansProg, os);
  MyMatrix<T> FAC = DirectFacetComputationInequalities(EXT, ansProg, os);
  os << "Obtained results:\n";
  os << "|vf|=" << vf.n << " / " << vf.n_face << "\n";
  size_t pos = 0;
  for (Face f : vf) {
    os << "pos=" << pos << " f=" << StringFace(f) << "\n";
  }
  WriteMatrix(os, FAC);
}

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 4 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_dual_description arith command [DATAIN] [DATAOUT]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_dual_description arith command [DATAIN]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "choice  : the chosen processing option\n";
      std::cerr << "arith   : the chosen arithmetic\n";
      std::cerr << "DATAIN  : The polyhedral cone inequalities\n";
      std::cerr << "DATAOUT : The file of output (if present, otherwise std::cerr)\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "rational : rational arithmetic on input\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr << "RealAlgebraic=FileDesc  : For the real algebraic case of a "
                   "field whose description is in FileDesc\n";
      std::cerr << "\n";
      std::cerr << "        --- command ---\n";
      std::cerr << "\n";
      std::cerr << "cdd      : the cdd program\n";
      std::cerr << "lrs      : the lrs program\n";
      std::cerr << "lrs_ring : the lrs program but reduced to ring computation (remove denominators)\n";
      std::cerr << "glrs     : the external program glrs\n";
      std::cerr << "ppl_ext  : the external program ppl_ext\n";
      std::cerr << "cdd_ext  : the external program cdd_ext\n";
      std::cerr << "normaliz : the external program normaliz\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string command = argv[2];
    std::string eFileI = argv[3];
    std::string eFileO = "stderr";
    if (argc == 5)
      eFileO = argv[4];
    auto call_lrs = [&](std::ostream &os) -> void {
      if (arith == "rational") {
        return process<mpq_class>(eFileI, command, os);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return process<T>(eFileI, command, os);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return process<T>(eFileI, command, os);
      }
      std::optional<std::string> opt_realalgebraic =
          get_postfix(arith, "RealAlgebraic=");
      if (opt_realalgebraic) {
        using T_rat = mpq_class;
        std::string FileAlgebraicField = *opt_realalgebraic;
        if (!IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                    << " is missing\n";
          throw TerminalException{1};
        }
        HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        using T = RealField<idx_real_algebraic_field>;
        return process<T>(eFileI, command, os);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
                   "RealAlgebraic\n";
      throw TerminalException{1};
    };
    if (eFileO == "stderr") {
      call_lrs(std::cerr);
    } else {
      if (eFileO == "stdout") {
        call_lrs(std::cout);
      } else {
        std::ofstream os(eFileO);
        call_lrs(os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_dual_description\n";
    exit(e.eVal);
  }
  runtime(time1);
}

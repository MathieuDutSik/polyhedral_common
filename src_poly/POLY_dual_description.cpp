// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "POLY_DirectDualDesc.h"
// clang-format on

template <typename T>
void process(std::string const &eFileI, std::string const &ansProg,
             std::string const &choice, std::ostream &os) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileI);
  if (choice == "control") {
    vectface vf = DirectFacetComputationIncidence(EXT, ansProg, std::cerr);
    MyMatrix<T> FAC =
        DirectFacetComputationInequalities(EXT, ansProg, std::cerr);
    os << "Obtained results:\n";
    os << "|vf|=" << vf.n << " / " << vf.n_face << "\n";
    size_t pos = 0;
    for (Face f : vf) {
      os << "pos=" << pos << " f=" << StringFace(f) << " |f|=" << f.count()
         << "\n";
    }
    return WriteMatrix(os, FAC);
  }
  if (choice == "GAP") {
    MyMatrix<T> FAC =
        DirectFacetComputationInequalities(EXT, ansProg, std::cerr);
    os << "return ";
    WriteMatrixGAP(os, FAC);
    os << ";\n";
    return;
  }
  if (choice == "CPP") {
    MyMatrix<T> FAC =
        DirectFacetComputationInequalities(EXT, ansProg, std::cerr);
    return WriteMatrix(os, FAC);
  }
  std::cerr << "choice=" << choice
            << " but allowed possibilities are control and CPP\n";
  std::cerr << "Failed to find a matching entry\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
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
      std::cerr << "safe_rational          : rational arithmetic based on "
                   "int64_t that fails\n";
      std::cerr << "    gracefully on overflowing\n";
      std::cerr << "cpp_rational           : rational arithmetic based on "
                   "boost header library\n";
      std::cerr << "mpq_rational           : rational arithmetic based on "
                   "boost mpq data type\n";
      std::cerr << "rational               : rational arithmetic on input\n";
      std::cerr
          << "Qsqrt2                 : arithmetic over the field Q(sqrt(2))\n";
      std::cerr
          << "Qsqrt5                 : arithmetic over the field Q(sqrt(5))\n";
      std::cerr
          << "RealAlgebraic=FileDesc : For the real algebraic case of a\n";
      std::cerr << "    field whose description is in FileDesc\n";
      std::cerr << "\n";
      std::cerr << "        --- command ---\n";
      std::cerr << "\n";
      std::cerr << "cdd      : the cdd program\n";
      std::cerr << "lrs      : the lrs program\n";
      std::cerr << "pd_lrs   : the pd with lrs used for checks\n";
      std::cerr << "lrs_ring : the lrs program but reduced to ring computation "
                   "(remove denominators)\n";
      std::cerr << "glrs     : the external program glrs\n";
      std::cerr << "ppl_ext  : the external program ppl_ext\n";
      std::cerr << "cdd_ext  : the external program cdd_ext\n";
      std::cerr << "normaliz : the external program normaliz\n";
      std::cerr << "\n";
      std::cerr << "        --- choice ---\n";
      std::cerr << "\n";
      std::cerr << "control : the full data set for control and debugging\n";
      std::cerr << "CPP     : the matrix for output (also used in Oscar)\n";
      std::cerr << "GAP     : returns a GAP readable file (except for "
                   "algebraic field)\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string command = argv[2];
    std::string choice = argv[3];
    std::string eFileI = argv[4];
    std::string eFileO = "stderr";
    if (argc == 6)
      eFileO = argv[5];
    auto dual_desc = [&](std::ostream &os) -> void {
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(eFileI, command, choice, os);
      }
      if (arith == "cpp_rational") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(eFileI, command, choice, os);
      }
      if (arith == "mpq_rational") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(eFileI, command, choice, os);
      }
      if (arith == "rational") {
        using T = mpq_class;
        return process<T>(eFileI, command, choice, os);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return process<T>(eFileI, command, choice, os);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return process<T>(eFileI, command, choice, os);
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
        return process<T>(eFileI, command, choice, os);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
                   "RealAlgebraic\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(eFileO, dual_desc);
    std::cerr << "Normal termination of POLY_dual_description\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_dual_description\n";
    exit(e.eVal);
  }
  runtime(time1);
}

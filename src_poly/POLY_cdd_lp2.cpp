// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "QuadField.h"
#include "POLY_LinearProgramming.h"
// clang-format on

template <typename T>
void process(std::string const &eFileFAC, std::string const& eFileIneq, std::ostream &os) {
  MyMatrix<T> TheEXT = ReadMatrixFile<T>(eFileFAC);
  MyVector<T> eVect = ReadVectorFile<T>(eFileIneq);
  if (TheEXT.cols() != eVect.size()) {
    std::cerr << "|TheEXT|=" << TheEXT.rows() << " / " << TheEXT.cols() << "\n";
    std::cerr << "|eVect|=" << eVect.size() << "\n";
    std::cerr << "This is inconsistentt\n";
    throw TerminalException{1};
  }
  LpSolution<T> eSol = CDD_LinearProgramming(TheEXT, eVect);
  os << "return rec(";
  os << "answer:=\"" << eSol.Answer << "\",\n";
  os << "OptimalValue:=" << eSol.OptimalValue;
  if (eSol.PrimalDefined) {
    os << ",\n primal_solution:=" << StringVectorGAP(eSol.DirectSolution);
  }
  if (eSol.DualDefined) {
    os << ",\n dual_solution:=" << StringVectorGAP(eSol.DualSolution);
  }
  os << ");\n";
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 4 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_cdd_lp2 arith [FileFAC] [FileIneq] [DATAOUT]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_cdd_lp2 arith [FileFAC] [FileIneq]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "arith    : the chosen arithmetic\n";
      std::cerr << "FileFAC  : The list of inequalities\n";
      std::cerr << "FileIneq : The inequality to be minimized\n";
      std::cerr << "DATAOUT  : The file of output (if present, otherwise std::cerr)\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "rational : rational arithmetic on input\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr << "RealAlgebraic=FileDesc  : For the real algebraic case of a";
      std::cerr << "  field whose description is in FileDesc\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string eFileFAC = argv[2];
    std::string eFileIneq = argv[3];
    std::string eFileO = "stderr";
    if (argc == 5)
      eFileO = argv[4];
    auto call_lp = [&](std::ostream &os) -> void {
      if (arith == "rational") {
        using T = mpq_class;
        return process<T>(eFileFAC, eFileIneq, os);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return process<T>(eFileFAC, eFileIneq, os);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return process<T>(eFileFAC, eFileIneq, os);
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
        return process<T>(eFileFAC, eFileIneq, os);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
                   "RealAlgebraic\n";
      throw TerminalException{1};
    };
    if (eFileO == "stderr") {
      call_lp(std::cerr);
    } else {
      if (eFileO == "stdout") {
        call_lp(std::cout);
      } else {
        std::ofstream os(eFileO);
        call_lp(os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_dual_description\n";
    exit(e.eVal);
  }
  runtime(time1);
}
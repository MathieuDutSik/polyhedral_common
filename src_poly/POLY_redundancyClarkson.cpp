// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "QuadField.h"
#include "POLY_cddlib.h"
// clang-format on

template<typename T>
void process_A(std::string const& eFileI, std::string const& eFileO, std::string const& choice) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileI);
  std::vector<int> ListIrred = cdd::RedundancyReductionClarkson(EXT);
  auto print_result=[&](std::ostream& os) -> void {
    if (choice == "GAP") {
      os << "return [";
      int nbIrred = ListIrred.size();
      for (int i = 0; i < nbIrred; i++) {
        if (i > 0)
          os << ",";
        int eVal = ListIrred[i] + 1;
        os << eVal;
      }
      os << "];\n";
    }
    if (choice == "Python") {
      int nbIrred = ListIrred.size();
      for (int i = 0; i < nbIrred; i++) {
        if (i > 0)
          os << " ";
        int eVal = ListIrred[i];
        os << eVal;
      }
    }
    std::cerr << "Failed to find a matching entry\n";
    throw TerminalException{1};
  };
  if (eFileO == "stderr")
    return print_result(std::cerr);
  if (eFileO == "stdout")
    return print_result(std::cout);
  std::ofstream os(eFileO);
  return print_result(os);
}



void process_B(std::string const& eFileI, std::string const& eFileO, std::string const& choice, std::string const& arith) {
  if (arith == "rational") {
    return process_A<mpq_class>(eFileI, eFileO, choice);
  }
  if (arith == "Qsqrt5") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return process_A<T>(eFileI, eFileO, choice);
  }
  if (arith == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return process_A<T>(eFileI, eFileO, choice);
  }
  std::optional<std::string> opt_realalgebraic = get_postfix(arith, "RealAlgebraic=");
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
    return process_A<T>(eFileI, eFileO, choice);
  }
  std::cerr << "Failed to find a matching field for arith=" << arith << "\n";
  std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
    "RealAlgebraic\n";
  throw TerminalException{1};
}



int main(int argc, char *argv[]) {
  try {
    if (argc != 4 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_redundancyClarskson choice arith [DATAIN] [DATAOUT]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_redundancyClarskson choice arith [DATAIN]\n";
      std::cerr << "\n";
      std::cerr << "choice : the choice for presenting the output\n";
      std::cerr << "arith : the chosen arithmetic for the output\n";
      std::cerr << "DATAIN : The polyhedral cone inequalities\n";
      std::cerr << "DATAOUT : The list of irredundant facets\n";
      std::cerr << "\n";
      std::cerr << "     ---- choice ----\n";
      std::cerr << "\n";
      std::cerr << "GAP : For having a gap readiable file (via ReadAsFunction)\n";
      std::cerr << "Python : For having a python readable file\n";
      std::cerr << "\n";
      std::cerr << "     ---- arith ----\n";
      std::cerr << "\n";
      std::cerr << "rational : rational arithmetic on input\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr << "RealAlgebraic=FileDesc  : For the real algebraic case of a field whose description is in FileDesc\n";
      return -1;
    }
    if (argc == 4)
      process_B(argv[3], "stderr", argv[1], argv[2]);
    if (argc == 5)
      process_B(argv[3], argv[4], argv[1], argv[2]);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}

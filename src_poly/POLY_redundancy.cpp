// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "POLY_cddlib.h"
#include "POLY_lrslib.h"
#include "POLY_RedundancyElimination.h"
// clang-format on

template <typename T>
void process_A(std::string const &eFileI, std::string const &eFileO,
               std::string const &method, std::string const &OutFormat,
               std::ostream& os) {
  MyMatrix<T> preEXT = ReadMatrixFile<T>(eFileI);
  MyMatrix<T> EXT = lrs::FirstColumnZeroCond(preEXT).first;
  auto get_list_irred = [&]() -> std::vector<int> {
    if (method == "Clarkson") {
      return cdd::RedundancyReductionClarkson(EXT, os);
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
  if (eFileO == "stderr")
    return print_result(std::cerr);
  if (eFileO == "stdout")
    return print_result(std::cout);
  std::ofstream osF(eFileO);
  return print_result(osF);
}

void process_B(std::string const &eFileI, std::string const &eFileO,
               std::string const &method, std::string const &OutFormat,
               std::string const &arith, std::ostream& os) {
  if (arith == "safe_rational") {
    using T = Rational<SafeInt64>;
    return process_A<T>(eFileI, eFileO, method, OutFormat, os);
  }
  if (arith == "rational") {
    using T = mpq_class;
    return process_A<T>(eFileI, eFileO, method, OutFormat, os);
  }
  if (arith == "Qsqrt5") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return process_A<T>(eFileI, eFileO, method, OutFormat, os);
  }
  if (arith == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return process_A<T>(eFileI, eFileO, method, OutFormat, os);
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
    return process_A<T>(eFileI, eFileO, method, OutFormat, os);
  }
  std::cerr << "Failed to find a matching field for arith=" << arith << "\n";
  std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
               "RealAlgebraic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
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
      std::cerr << "safe_rational : rational arithmetic based on int64_t that "
                   "fails\n";
      std::cerr << "    gracefully on overflow\n";
      std::cerr << "rational : rational arithmetic on input\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr
          << "RealAlgebraic=FileDesc  : For the real algebraic case of a\n";
      std::cerr << "    field whose description is in FileDesc\n";
      return -1;
    }
    std::string method = argv[1];
    std::string arith = argv[2];
    std::string eFileI = argv[3];
    std::string OutFormat = "GAP";
    std::string eFileO = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      eFileO = argv[5];
    }
    process_B(eFileI, eFileO, method, OutFormat, arith, std::cerr);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happenned somewhere\n";
    exit(e.eVal);
  }
  runtime(time1);
}

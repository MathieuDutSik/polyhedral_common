// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "GRP_GroupFct.h"
#include "Group.h"
#include "POLY_cddlib.h"
#include "POLY_lrslib.h"
#include "POLY_RedundancyElimination.h"
#include "Permutation.h"
// clang-format on

template <typename T, typename Tgroup>
void process(std::string const &FileEXT, std::string const &FileGRP,
             std::string const &method, std::string const &OutFormat,
             std::ostream &os_out, std::ostream& os) {
  MyMatrix<T> preEXT = ReadMatrixFile<T>(FileEXT);
  MyMatrix<T> EXT = lrs::FirstColumnZeroCond(preEXT).first;
  size_t nbRow = EXT.rows();
  Tgroup GRP = ReadGroupFile<Tgroup>(FileGRP);
  std::cerr << "|GRP|=" << GRP.size() << " nbRow=" << nbRow << "\n";
  auto get_list_irred = [&]() -> std::vector<int> {
    if (method == "ClarksonBlock") {
      vectface vfo = DecomposeOrbitPoint_Full(GRP);
      size_t n_orbit = vfo.size();
      std::vector<int> BlockBelong(nbRow);
      for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
        Face f = vfo[i_orbit];
        //    std::cerr << "i_orbit=" << i_orbit << "/" << n_orbit << " |f|=" <<
        //    f.size() << "/" << f.count() << "\n";
        for (size_t i = 0; i < nbRow; i++) {
          if (f[i] == 1) {
            BlockBelong[i] = i_orbit;
          }
        }
      }
      return cdd::RedundancyReductionClarksonBlocks(EXT, BlockBelong, os);
    }
    if (method == "Equivariant") {
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
  HumanTime time1;
  try {
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
      std::cerr << "ClarksonBlock : For using Clarkson method wwith blocks\n";
      std::cerr << "Equivariant   : For using the equivariant method\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "safe_rational  : rational arithmetic based on int64_t that "
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
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    std::string method = argv[1];
    std::string arith = argv[2];
    std::string FileEXT = argv[3];
    std::string FileGRP = argv[4];
    std::string OutFormat = "GAP";
    std::string FileOut = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      FileOut = argv[6];
    }
    auto compute_redundancy = [&](std::ostream &os_out) -> void {
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T, Tgroup>(FileEXT, FileGRP, method, OutFormat, os_out, std::cerr);
      }
      if (arith == "rational") {
        using T = mpq_class;
        return process<T, Tgroup>(FileEXT, FileGRP, method, OutFormat, os_out, std::cerr);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return process<T, Tgroup>(FileEXT, FileGRP, method, OutFormat, os_out, std::cerr);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return process<T, Tgroup>(FileEXT, FileGRP, method, OutFormat, os_out, std::cerr);
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
        return process<T, Tgroup>(FileEXT, FileGRP, method, OutFormat, os_out, std::cerr);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
                   "RealAlgebraic\n";
      throw TerminalException{1};
    };
    if (FileOut == "stderr") {
      compute_redundancy(std::cerr);
    } else {
      if (FileOut == "stdout") {
        compute_redundancy(std::cout);
      } else {
        std::ofstream osF(FileOut);
        compute_redundancy(osF);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_redundancyClarksonBlocks\n";
    exit(e.eVal);
  }
  runtime(time1);
}

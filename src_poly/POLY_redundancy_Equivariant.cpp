// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheorySafeInt.h"
#include "Group.h"
#include "Permutation.h"
#include "POLY_RedundancyElimination.h"
// clang-format on

template<typename T, typename Tgroup>
void process(std::string const& FileExtI, std::string const& FileGrpI, std::string const& OutFormat, std::ostream & os) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(FileExtI);
  Tgroup GRP = ReadGroupFile<Tgroup>(FileGrpI);
  Face status = GetNonRedundant_Equivariant(EXT, GRP);
  auto print = [&]() -> void {
    size_t nbIrred = status.count();
    if (OutFormat == "GAP") {
      os << "return [";
      std::cerr << "nbIrred=" << nbIrred << "\n";
      boost::dynamic_bitset<>::size_type pos = status.find_first();
      for (size_t i = 0; i < nbIrred; i++) {
        if (i > 0)
          os << ",";
        int eVal = pos + 1;
        os << eVal;
        pos = status.find_next(pos);
      }
      os << "];\n";
      return;
    }
    if (OutFormat == "Raw") {
      os << nbIrred << "\n";
      boost::dynamic_bitset<>::size_type pos = status.find_first();
      for (size_t i = 0; i < nbIrred; i++) {
        if (i > 0)
          os << ",";
        int eVal = pos;
        os << " " << eVal;
        pos = status.find_next(pos);
      }
      os << "\n";
      return;
    }
    if (OutFormat == "NB") {
      os << "Total number of non-redundant facets=" << nbIrred << "\n";
      return;
    }
    std::cerr << "Failed to find a matching OutFormat. OutFormat=" << OutFormat << "\n";
    std::cerr << "Allowed values: GAP, Raw, NB\n";
    throw TerminalException{1};
  };
  print();
}


int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_redundancy_Equivariant arith [FileExtI] [FileGrpI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_redundancy_Equivariant arith [FileExtI] [FileGrpI]\n";
      std::cerr << "\n";
      std::cerr << "FileExtI  : The polyhedral cone inequalities (or vertices)\n";
      std::cerr << "FileGrpI  : The permutation group on the inequalities (or vertices)\n";
      std::cerr << "OutFormat : The format of the Output. Possibilities: GAP, Raw, NB\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "safe_rational          : rational arithmetic based on int64_t that fails\n";
      std::cerr << "    gracefully on overflow\n";
      std::cerr << "rational               : rational arithmetic based on gmp\n";
      std::cerr << "Qsqrt2                 : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5                 : arithmetic over the field Q(sqrt(5))\n";
      std::cerr << "RealAlgebraic=FileDesc : For the real algebraic case of a\n";
      std::cerr << "    field whose description is in FileDesc\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string FileExtI = argv[2];
    std::string FileGrpI = argv[3];
    std::string OutFormat = "NB";
    std::string FileO = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      FileO = argv[5];
    }
    auto treat=[&](std::ostream & os) -> void {
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T,Tgroup>(FileExtI, FileGrpI, OutFormat, os);
      }
      if (arith == "rational") {
        using T = mpq_class;
        return process<T,Tgroup>(FileExtI, FileGrpI, OutFormat, os);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return process<T,Tgroup>(FileExtI, FileGrpI, OutFormat, os);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return process<T,Tgroup>(FileExtI, FileGrpI, OutFormat, os);
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
        return process<T,Tgroup>(FileExtI, FileGrpI, OutFormat, os);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
                   "RealAlgebraic\n";
      throw TerminalException{1};
    };
    if (FileO == "stderr") {
      treat(std::cerr);
    } else {
      if (FileO == "stdout") {
        treat(std::cout);
      } else {
        std::ofstream os(FileO);
        treat(os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_redundancy_Equivariant\n";
    exit(e.eVal);
  }
  runtime(time1);
}

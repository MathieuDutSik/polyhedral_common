// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "Group.h"
#include "POLY_Kskeletton.h"
#include "Permutation.h"
// clang-format on

template <typename T, typename Tgroup>
void MainFunctionFaceLattice(std::string const &FACfile,
                             std::string const &GRPfile, int const &LevSearch,
                             std::string const &OutFormat,
                             std::ostream &os_out) {
  MyMatrix<T> FAC = ReadMatrixFile<T>(FACfile);
  Tgroup GRP = ReadGroupFile<Tgroup>(GRPfile);
  std::string method_spann = "LinearProgramming";
  std::string method_final = "all";
  MyMatrix<T> EXT; // Unused since we are using linear programming
  bool ComputeTotalNumberFaces = false;
  std::vector<vectface> TheOutput =
      EnumerationFaces(GRP, FAC, EXT, LevSearch, method_spann, method_final,
                       ComputeTotalNumberFaces, std::cerr);
  OutputFaces_stream(TheOutput, os_out, OutFormat);
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    //
    if (argc != 5 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "POLY_DirectFaceLattice [arith] [FACfile] [GRPfile] [LevSearch]\n";
      std::cerr << "     os\n";
      std::cerr << "POLY_DirectFaceLattice [arith] [FACfile] [GRPfile] "
                   "[LevSearch] [OutFormat] [OutFile]\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FACfile = argv[2];
    std::string GRPfile = argv[3];
    std::string LevSearch_str = argv[4];
    int LevSearch = ParseScalar<int>(LevSearch_str);
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    auto f = [&](std::ostream &os_out) -> void {
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return MainFunctionFaceLattice<T, Tgroup>(FACfile, GRPfile, LevSearch,
                                                  OutFormat, os_out);
      }
      if (arith == "rational") {
        using T = mpq_class;
        return MainFunctionFaceLattice<T, Tgroup>(FACfile, GRPfile, LevSearch,
                                                  OutFormat, os_out);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return MainFunctionFaceLattice<T, Tgroup>(FACfile, GRPfile, LevSearch,
                                                  OutFormat, os_out);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return MainFunctionFaceLattice<T, Tgroup>(FACfile, GRPfile, LevSearch,
                                                  OutFormat, os_out);
      }
      std::optional<std::string> opt_realalgebraic =
          get_postfix(arith, "RealAlgebraic=");
      if (opt_realalgebraic) {
        std::string const &FileAlgebraicField = *opt_realalgebraic;
        if (!IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                    << " is missing\n";
          throw TerminalException{1};
        }
        using T_rat = mpq_class;
        HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        using T = RealField<idx_real_algebraic_field>;
        return MainFunctionFaceLattice<T, Tgroup>(FACfile, GRPfile, LevSearch,
                                                  OutFormat, os_out);
      }
      std::cerr << "Failed to find a matching arithmetic arith=" << arith
                << "\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(OutFile, f);
    //
    std::cerr << "Normal termination of POLY_DirectFaceLattice\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_DirectFaceLattice\n";
    exit(e.eVal);
  }
  runtime(time);
}

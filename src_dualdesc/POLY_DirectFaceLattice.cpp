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
#include "Group.h"
#include "POLY_Kskeletton.h"
#include "GRP_GroupFile.h"
#include "Permutation.h"
#include "Arithmetic_dispatch_types.h"
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
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    //
    std::string FACfile;
    std::string GRPfile;
    int LevSearch = 0;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      MainFunctionFaceLattice<T, Tgroup>(FACfile, GRPfile, LevSearch,
                                         OutFormat, *os_ptr);
    };
    add_rational_arithmetics(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
    if (argc != 5 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "POLY_DirectFaceLattice [arith] [FACfile] [GRPfile] [LevSearch]\n";
      std::cerr << "     os\n";
      std::cerr << "POLY_DirectFaceLattice [arith] [FACfile] [GRPfile] "
                   "[LevSearch] [OutFormat] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    std::string arith = argv[1];
    FACfile = argv[2];
    GRPfile = argv[3];
    std::string LevSearch_str = argv[4];
    LevSearch = ParseScalar<int>(LevSearch_str);
    if (argc == 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    auto f = [&](std::ostream &os_out) -> void {
      os_ptr = &os_out;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    //
    std::cerr << "Normal termination of POLY_DirectFaceLattice\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_DirectFaceLattice\n";
    exit(e.eVal);
  }
  runtime(time);
}

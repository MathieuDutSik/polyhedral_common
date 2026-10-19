// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "GRP_GroupFct.h"
#include "PolytopeEquiStab.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const &FileExt1, std::string const &FileExt2,
             std::string const &OutFormat, std::string const &FileO) {
  using Tidx = uint32_t;
  MyMatrix<T> EXT1 = ReadMatrixFile<T>(FileExt1);
  MyMatrix<T> EXT2 = ReadMatrixFile<T>(FileExt2);
  size_t nbRow = EXT1.rows();
  //
  std::optional<std::vector<Tidx>> equiv =
      LinPolytope_Isomorphism<T, Tidx>(EXT1, EXT2, std::cerr);
  //
  auto print_info = [&](std::ostream &os) -> void {
    if (OutFormat == "Oscar") {
      if (equiv) {
        const std::vector<Tidx> &V = *equiv;
        os << nbRow << "\n";
        for (size_t iRow = 0; iRow < nbRow; iRow++) {
          if (iRow > 0)
            os << " ";
          int eVal = V[iRow] + 1;
          os << eVal;
        }
      } else {
        os << "-1\n";
      }
      return;
    }
    if (OutFormat == "GAP") {
      if (equiv) {
        const std::vector<Tidx> &V = *equiv;
        os << "return [";
        for (size_t iRow = 0; iRow < nbRow; iRow++) {
          if (iRow > 0)
            os << ",";
          int eVal = V[iRow] + 1;
          os << eVal;
        }
        os << "];\n";
      } else {
        os << "return fail;\n";
      }
      return;
    }
    std::cerr << "Failed to find a matching entry\n";
    throw TerminalException{1};
  };
  FILE_PrintStderrStdoutFile(FileO, print_info);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string FileExt1;
    std::string FileExt2;
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(FileExt1, FileExt2, OutFormat, FileO);
    };
    add_rational_arithmetics(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
#ifdef QUADRATIC_FIELDS
    dispatch.add("Qsqrt3", "arithmetic over the field Q(sqrt(3))",
                 [&]() -> void {
                   f_arith.template operator()<QuadField<Trat_quad_field, 3>>();
                 });
#endif
    if (argc != 6 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "GRP_LinPolytope_Isomorphism [arith] [EXT1] [EXT2] [OutFormat] "
             "[OutEquiv]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytope_Isomorphism [arith] [EXT1] [EXT2]\n";
      std::cerr << "\n";
      std::cerr << "OutEquiv : The equivalence information file (otherwise "
                   "printed to screen)\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    //
    std::string arith = argv[1];
    FileExt1 = argv[2];
    FileExt2 = argv[3];
    if (argc == 6) {
      OutFormat = argv[4];
      FileO = argv[5];
    }
    dispatch.run(arith);
    std::cerr << "Normal termination of GRP_LinPolytope_Isomorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Isomorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}

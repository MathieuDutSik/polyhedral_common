// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GRP_GroupFct.h"
#include "Group.h"
#include "Permutation.h"
#include "PolytopeEquiStab.h"
// clang-format on

template <typename Tint>
void process_A(std::string const &FileExt1, std::string const &FileExt2,
               std::string const &OutFormat, std::ostream &os) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tgroup = permutalib::Group<Telt, Tint>;
  using Tidx_value = uint32_t;
  //    using Tgr = GraphBitset;
  using Tgr = GraphListAdj;
  MyMatrix<Tint> EXT1 = ReadMatrixFile<Tint>(FileExt1);
  MyMatrix<Tint> EXT2 = ReadMatrixFile<Tint>(FileExt2);
  size_t nbCol = EXT1.cols();
  size_t nbRow = EXT1.rows();
  std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
  //
  std::optional<MyMatrix<Tint>> equiv =
      LinPolytopeIntegral_Isomorphism<Tint, Tidx, Tgroup, Tidx_value, Tgr>(EXT1, EXT2, std::cerr);
  if (OutFormat == "GAP") {
    if (equiv) {
      os << "return ";
      WriteMatrixGAP(os, *equiv);
      os << ";\n";
    } else {
      os << "return fail;\n";
    }
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 6 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinPolytopeIntegral_Isomorphism arith [EXT1] [EXT2] "
                   "[OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytopeIntegral_Isomorphism arith [EXT1] [EXT2]\n";
      std::cerr << "\n";
      std::cerr << "         ------ arith -------\n";
      std::cerr << "\n";
      std::cerr << "mpz_class   : The GMP class with the GMP binder\n";
      std::cerr << "mpq_integer : The GMP class embedded into boost\n";
      std::cerr << "\n";
      std::cerr << "EXT1 and EXT2 are files of the input files\n";
      std::cerr << "\n";
      std::cerr << "         ----- OutFormat ------\n";
      std::cerr << "\n";
      std::cerr << "GAP: The output in the GAP readable file\n";
      std::cerr << "\n";
      std::cerr << "         ----- FileOut -----\n";
      std::cerr << "\n";
      std::cerr << "If stderr, or stdout, then output to standard error or "
                   "standrd output\n";
      std::cerr << "Other output to the designated file name\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string FileExt1 = argv[2];
    std::string FileExt2 = argv[3];
    std::string OutFormat = "GAP";
    std::string FileOut = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      FileOut = argv[5];
    }
    //
    auto process_B = [&](std::ostream &os) -> void {
      if (arith == "mpz_class") {
        using Tint = mpz_class;
        return process_A<Tint>(FileExt1, FileExt2, OutFormat, os);
      }
      if (arith == "mpq_integer") {
        using Tint = mpz_class;
        return process_A<Tint>(FileExt1, FileExt2, OutFormat, os);
      }
      std::cerr << "Failed to find a matching type for arith\n";
      throw TerminalException{1};
    };
    if (FileOut == "stderr") {
      process_B(std::cerr);
    } else {
      if (FileOut == "stdout") {
        process_B(std::cout);
      } else {
        std::ofstream os(FileOut);
        process_B(os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytopeIntegral_Isomorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}

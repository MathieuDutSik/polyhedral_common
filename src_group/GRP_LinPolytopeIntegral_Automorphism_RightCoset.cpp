// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GRP_GroupFct.h"
#include "Group.h"
#include "Permutation.h"
#include "Temp_PolytopeEquiStab.h"
// clang-format on

template<typename Tint>
void process_A(std::string const& FileExt, std::string const& OutFormat, std::ostream & os) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tgroup = permutalib::Group<Telt, Tint>;
  using Tidx_value = uint32_t;
  //    using Tgr = GraphBitset;
  using Tgr = GraphListAdj;
  MyMatrix<Tint> EXT = ReadMatrixFile<Tint>(FileExt);
  size_t nbCol = EXT.cols();
  size_t nbRow = EXT.rows();
  std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
  //
  //    const bool use_scheme = true;
  const bool use_scheme = false;
  std::pair<Tgroup,std::vector<Telt>> pair = LinPolytopeIntegral_Automorphism_RightCoset<Tint, Tidx, Tgroup, Tidx_value, Tgr, use_scheme>(EXT, std::cerr);
  Tgroup GRP = pair.first;
  if (OutFormat == "GAP") {
    os << "return " << GRP.GapString() << ";\n";
    return;
  }
  if (OutFormat == "RecGAP") {
    std::string strGAPmatr = "[";
    bool IsFirst = true;
    for (auto & eGen : GRP.GeneratorsOfGroup()) {
      MyMatrix<Tint> M = RepresentVertexPermutation(EXT, EXT, eGen);
      if (!IsFirst)
        strGAPmatr += ",";
      strGAPmatr += StringMatrixGAP(M);
    }
    strGAPmatr += "]";
    os << "return rec(GAPperm:=" << GRP.GapString() << ", GAPmatr:=" << strGAPmatr << ");";
    return;
  }
  if (OutFormat == "Oscar") {
    WriteGroup(os, GRP);
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}


int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinPolytopeIntegral_Automorphism arith [EXT] [OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytopeIntegral_Automorphism arith [EXT]\n";
      std::cerr << "\n";
      std::cerr << "         ------ arith -------\n";
      std::cerr << "\n";
      std::cerr << "mpz_class   : The GMP class with the GMP binder\n";
      std::cerr << "mpq_integer : The GMP class embedded into boost\n";
      std::cerr << "\n";
      std::cerr << "EXT is file of the input polytope\n";
      std::cerr << "\n";
      std::cerr << "         ----- OutFormat ------\n";
      std::cerr << "\n";
      std::cerr << "GAP    : The permutation group as a GAP readable file\n";
      std::cerr << "RecGAP : The permutation group and the matrix group as a GAP readable file\n";
      std::cerr << "Oscar  : The output to the Oscar readable file\n";
      std::cerr << "\n";
      std::cerr << "         ----- FileOut -----\n";
      std::cerr << "\n";
      std::cerr << "If stderr, or stdout, then output to standard error or standrd output\n";
      std::cerr << "Other output to the designated file name\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string FileExt = argv[2];
    std::string OutFormat = "GAP";
    std::string FileOut = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileOut = argv[4];
    }
    //
    auto process_B=[&](std::ostream & os) -> void {
      if (arith == "rational") {
        using Tint = mpz_class;
        return process_A<Tint>(FileExt, OutFormat, os);
      }
      if (arith == "mpq_integer") {
        using Tint = mpz_class;
        return process_A<Tint>(FileExt, OutFormat, os);
      }
      std::cerr << "Failed to find a matching type for arith\n";
      std::cerr << "arith=" << arith << " allowed = rational and mpq_integer\n";
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

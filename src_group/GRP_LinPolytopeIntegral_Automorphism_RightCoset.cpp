// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GRP_GroupFct.h"
#include "Group.h"
#include "Permutation.h"
#include "PolytopeEquiStab.h"
// clang-format on

template <typename Tint>
void process_A(std::string const &FileExt, std::string const &OutFormat,
               std::ostream &os) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tgroup = permutalib::Group<Telt, Tint>;
  using Tidx_value = uint32_t;
  using Tfield = typename overlying_field<Tint>::field_type;
  MyMatrix<Tint> EXT = ReadMatrixFile<Tint>(FileExt);
  size_t nbCol = EXT.cols();
  size_t nbRow = EXT.rows();
  std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
  //
  std::pair<Tgroup, std::vector<Telt>> pair =
      LinPolytopeIntegral_Automorphism_RightCoset<Tint, Tgroup>(EXT, std::cerr);
  MyMatrix<Tfield> EXT_T = UniversalMatrixConversion<Tfield, Tint>(EXT);
  Tgroup GRPisom =
      LinPolytope_Automorphism<Tfield, Tgroup>(EXT_T, os);
  std::cerr << "|GRPisom|=" << GRPisom.size()
            << " |pair.first|=" << pair.first.size()
            << " |pair.second|=" << pair.second.size() << "\n";
  Tgroup GRP = pair.first;
  std::vector<Telt> l_elt = GRP.get_all_element();
  std::set<Telt> s_elt;
  for (auto &e_elt : l_elt) {
    for (auto &f_elt : pair.second) {
      Telt prod = e_elt * f_elt;
      if (!GRPisom.isin(prod)) {
        std::cerr << "The element is not in the product\n";
        throw TerminalException{1};
      }
      s_elt.insert(prod);
    }
  }
  Tint s_elt_size = s_elt.size();
  if (s_elt_size != GRPisom.size()) {
    std::cerr << "|s_elt|=" << s_elt.size() << " |GRPisom|=" << GRPisom.size()
              << "\n";
    std::cerr << "s_elt has the wrong size\n";
    throw TerminalException{1};
  }
  auto get_as_string = [&](std::vector<Telt> const &l_elt) -> std::string {
    std::string strGAPmatr = "[";
    bool IsFirst = true;
    for (auto &eElt : l_elt) {
      MyMatrix<Tint> M = RepresentVertexPermutation(EXT, EXT, eElt);
      if (!IsFirst)
        strGAPmatr += ",";
      strGAPmatr += StringMatrixGAP(M);
    }
    strGAPmatr += "]";
    return strGAPmatr;
  };
  if (OutFormat == "GAP") {
    os << "return " << GRP.GapString() << ";\n";
    return;
  }
  if (OutFormat == "RecGAP") {
    std::string strGAPgroup =
        "Group(" + get_as_string(GRP.GeneratorsOfGroup()) + ")";
    std::string strCoset = get_as_string(pair.second);
    os << "return rec(GAPperm:=" << GRP.GapString()
       << ", GAPmatr:=" << strGAPgroup << ", ListCoset:=" << strCoset << ");";
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
      std::cerr << "GRP_LinPolytopeIntegral_Automorphism arith [EXT] "
                   "[OutFormat] [FileOut]\n";
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
      std::cerr << "RecGAP : The permutation group and the matrix group as a "
                   "GAP readable file\n";
      std::cerr << "Oscar  : The output to the Oscar readable file\n";
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
    std::string FileExt = argv[2];
    std::string OutFormat = "GAP";
    std::string FileOut = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileOut = argv[4];
    }
    //
    auto process_B = [&](std::ostream &os) -> void {
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

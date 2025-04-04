// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GRP_GroupFct.h"
#include "Group.h"
#include "Permutation.h"
#include "PolytopeEquiStabInt.h"
// clang-format on

template <typename Tint>
void process_A(std::string const &FileExt, std::string const &OutFormat,
               std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tgroup = permutalib::Group<Telt, Tint>;
  using Tfield = typename overlying_field<Tint>::field_type;
  MyMatrix<Tint> EXT = ReadMatrixFile<Tint>(FileExt);
  //
  std::pair<Tgroup, std::vector<Telt>> pair =
      LinPolytopeIntegral_Automorphism_RightCoset<Tint, Tgroup>(EXT, std::cerr);
  MyMatrix<Tfield> EXT_T = UniversalMatrixConversion<Tfield, Tint>(EXT);
  Tgroup GRPisom = LinPolytope_Automorphism<Tfield, Tgroup>(EXT_T, std::cerr);
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
  if (OutFormat == "GAP") {
    os_out << "return " << GRP.GapString() << ";\n";
    return;
  }
  if (OutFormat == "RecGAP") {
    std::string strGAPgroup =
      "Group(" + get_matrs_as_string(EXT,  GRP.GeneratorsOfGroup()) + ")";
    std::string strCoset = get_matrs_as_string(EXT, pair.second);
    os_out << "return rec(GAPperm:=" << GRP.GapString()
           << ", GAPmatr:=" << strGAPgroup << ", ListCoset:=" << strCoset << ");";
    return;
  }
  if (OutFormat == "Oscar") {
    WriteGroup(os_out, GRP);
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
      std::cerr << "GRP_LinPolytopeIntegral_Automorphism_RightCoset arith [EXT] "
                   "[OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytopeIntegral_Automorphism_RightCoset arith [EXT]\n";
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
    print_stderr_stdout_file(FileOut, process_B);
    std::cerr << "Normal termination of GRP_LinPolytopeIntegral_Automorphism_RightCoset\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytopeIntegral_Automorphism_RightCoset\n";
    exit(e.eVal);
  }
  runtime(time1);
}

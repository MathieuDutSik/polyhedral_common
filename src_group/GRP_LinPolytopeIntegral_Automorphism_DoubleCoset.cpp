// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GRP_GroupFct.h"
#include "Group.h"
#include "Permutation.h"
#include "PolytopeEquiStabInt.h"
// clang-format on

template <typename Tint>
void process_A(std::string const &FileExt, std::string const& FileGrpV,
               std::string const &OutFormat, std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tgroup = permutalib::Group<Telt, Tint>;
  using Tfield = typename overlying_field<Tint>::field_type;
  MyMatrix<Tint> EXT = ReadMatrixFile<Tint>(FileExt);
  MyMatrix<Tfield> EXT_field = UniversalMatrixConversion<Tfield,Tint>(EXT);
  Tgroup GrpV = ReadGroupFile<Tgroup>(FileGrpV);
  //
  std::pair<Tgroup, std::vector<Telt>> pair =
    LinPolytopeIntegral_Automorphism_DoubleCoset<Tint, Tgroup>(EXT, GrpV, std::cerr);
  MyMatrix<Tfield> EXT_T = UniversalMatrixConversion<Tfield, Tint>(EXT);
  Tgroup const& GRP = pair.first;
  auto get_perms_as_string = [&](std::vector<Telt> const &l_elt) -> std::string {
    std::string strGAPperm = "[";
    bool IsFirst = true;
    for (auto &eElt : l_elt) {
      if (!IsFirst)
        strGAPperm += ",";
      IsFirst = false;
      strGAPperm += std::to_string(eElt);
    }
    strGAPperm += "]";
    return strGAPperm;
  };
  if (OutFormat == "GAP") {
    os_out << "return " << GRP.GapString() << ";\n";
    return;
  }
  if (OutFormat == "RecGAP") {
    std::string strGAPgroupMatr =
      "Group(" + get_matrs_as_string(EXT, GRP.GeneratorsOfGroup()) + ")";
    std::string strCosetMatr = get_matrs_as_string(EXT_field, pair.second);
    std::string strCosetPerm = get_perms_as_string(pair.second);
    os_out << "return rec(GAPperm:=" << GRP.GapString()
           << ", GAPmatr:=" << strGAPgroupMatr
           << ", DoubleCosetsMatr:=" << strCosetMatr
           << ", DoubleCosetsPerm:=" << strCosetPerm << ");";
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
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinPolytopeIntegral_Automorphism_DoubleCoset arith [EXT] [GRP_V] "
                   "[OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytopeIntegral_Automorphism_DoubleCoset arith [EXT] [GRP_V]\n";
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
    std::string FileGrpV = argv[3];
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
        return process_A<Tint>(FileExt, FileGrpV, OutFormat, os);
      }
      std::cerr << "Failed to find a matching type for arith\n";
      std::cerr << "arith=" << arith << " allowed = mpz_class\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(FileOut, process_B);
    std::cerr << "Normal termination of GRP_LinPolytopeIntegral_Automorphism_DoubleCoset\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytopeIntegral_Automorphism_DoubleCoset\n";
    exit(e.eVal);
  }
  runtime(time1);
}

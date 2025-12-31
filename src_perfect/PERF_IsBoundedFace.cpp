// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "perfect_complex.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

FullNamelist NAMELIST_get_is_bounded_face() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, std::string> ListStringValues1;
  ListStringValues1["arithmetic"] = "gmp";
  ListStringValues1["FileSHV"] = "unset";
  SingleBlock BlockDATA;
  BlockDATA.setListStringValues(ListStringValues1);
  ListBlock["DATA"] = BlockDATA;
  // SYSTEM
  ListBlock["SYSTEM"] = SINGLEBLOCK_Get_System();
  // TSPACE
  ListBlock["TSPACE"] = SINGLEBLOCK_Get_Tspace_Description();
  // Merging all data
  return FullNamelist(ListBlock);
}


template <typename T, typename Tint, typename Tgroup>
void process(FullNamelist const &eFull) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa =
      ReadTspace<T, Tint, Tgroup>(BlockTSPACE, std::cerr);
  //
  // SHV and is_bounded
  //
  std::string FileSHV = BlockDATA.get_string("FileSHV");
  MyMatrix<Tint> SHV = ReadMatrixFile<Tint>(FileSHV);
  bool result = is_bounded_face_iterative<T,Tint>(LinSpa, SHV, std::cerr);
  //
  // Output the data
  //
  std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  auto f_print = [&](std::ostream &os_out) -> void {
    if (OutFormat == "GAP") {
      os_out << "return " << GAP_logical(result) << ";\n";
      return;
    }
    if (OutFormat == "CPP") {
      os_out << result << "\n";
      return;
    }
    std::cerr << "PERF: No matching output for OutFormat=" << OutFormat << "\n";
    throw TerminalException{1};
  };
  print_stderr_stdout_file(OutFile, f_print);
}

int main(int argc, char *argv[]) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_get_is_bounded_face();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "PERF_IsBoundedFace [file.nml]\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    SingleBlock const &BlockDATA = eFull.get_block("DATA");
    std::string arithmetic = BlockDATA.get_string("arithmetic");
    auto f=[&]() -> void {
      if (arithmetic == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T,Tint,Tgroup>(eFull);
      }
      std::cerr << "Failed to find a matching entry for arithmetic=" << arithmetic << "\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of PERF_IsBoundedFace\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PERF_IsBoundedFace\n";
    exit(e.eVal);
  }
  runtime(time);
}

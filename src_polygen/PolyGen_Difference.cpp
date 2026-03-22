// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "generalized_polytopes.h"
// clang-format on

template <typename T>
GeneralizedPolytope<T> list_ext_to_generalizedpolytope(std::vector<MyMatrix<T>> const& list_ext) {
  std::vector<SinglePolytope<T>> polytopes;
  for (size_t i=0; i<list_ext.size(); i++) {
    MyMatrix<T> FAC = DirectDualDescription(EXT, std::cerr);
    SinglePolytope<T> sp = get_single_polytope(FAC, EXT);
    polytopes.push_back(sp);
  }
  GeneralizedPolytope<T> gp{polytopes};
  return gp;
}



template <typename T, typename Tint>
void process(std::string const &GenPolyFile1,
             std::string const &GenPolyFile2,
             std::string const &OutFormat,
             std::string const &OutFile) {
  std::vector<MyMatrix<T>> list_ext1 = ReadListMatrixFile<T>(GenPolyFile1);
  std::vector<MyMatrix<T>> list_ext2 = ReadListMatrixFile<T>(GenPolyFile2);
  GeneralizedPolytope<T> gp1 = list_ext_to_generalizedpolytope(list_ext1);
  GeneralizedPolytope<T> gp2 = list_ext_to_generalizedpolytope(list_ext2);
  GeneralizedPolytope<T> gp_d12 = difference_gp_gp(gp1, gp2, std::cerr);
  GeneralizedPolytope<T> gp_d21 = difference_gp_gp(gp1, gp2, std::cerr);
  GeneralizedPolytope<T> gp_int = intersection_gp_gp(gp1, gp2, std::cerr);


  T vol1 = volume_gp(gp1, std::cerr);
  T vol2 = volume_gp(gp2, std::cerr);
  T vol_int = volume_gp(gp_int, std::cerr);
  T vol_d12 = volume_gp(gp_d12, std::cerr);
  T vol_d21 = volume_gp(gp_d21, std::cerr);
  T vol1_sum = vol_int + vol_d12;
  T vol2_sum = vol_int + vol_d21;
  if (vol1 != vol1_sum) {
    std::cerr << "Inconsistent sums. vol1=" << vol1 << " vol1_sum=" << vol1_sum << "\n";
    throw TerminalException{1};
  }
  if (vol2 != vol2_sum) {
    std::cerr << "Inconsistent sums. vol2=" << vol2 << " vol2_sum=" << vol2_sum << "\n";
    throw TerminalException{1};
  }
  auto f_print=[&](std::ostream& os_out) -> void {
    if (OutFormat == "GAP") {
      os_out << "return rec(vol1:=" << vol1
             << ", vol2:=" << vol2
             << ", vol_int:=" << vol_int
             << ", vol_d12:=" << vol_d12
             << ", vol_d21:=" << vol_d21 << ");\n";
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
    throw TerminalException{1};
  };
  print_stderr_stdout_file(OutFile, f_print);
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 5 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "PolyGen_difference arith [ListEXT1] [ListEXT2] [OutFormat] [OutFile]\n";
      std::cerr << "       or\n";
      std::cerr << "PolyGen_difference arith [ListEXT1] [ListEXT2]\n";
      std::cerr << "\n";
      std::cerr << "allowed choices:\n";
      std::cerr << "arithmetic: gmp\n";
      std::cerr << "OutFormat: GAP\n";
      std::cerr << "OutFile: stderr, stdout, my_file\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string PolyFile1 = argv[2];
    std::string PolyFile2 = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    auto f=[&]() -> void {
      if (arith == "mpq_class") {
        using T = mpq_class;
        return process<T>(PolyFile1, PolyFile2, OutFormat, OutFile);
      }
      std::cerr << "Error for the template parameter arith=" << arith << "\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination PolyGen_difference\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PolyGen_difference\n";
    exit(e.eVal);
  }
  runtime(time);
}

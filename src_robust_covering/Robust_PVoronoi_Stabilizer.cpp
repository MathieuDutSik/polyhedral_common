// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "Permutation.h"
#include "Group.h"
#include "NumberTheory.h"
#include "enum_robust_covering.h"
// clang-format on

template <typename T, typename Tint>
void process(std::string const &MatFile,
             std::string const &PVoronoiFile,
             std::string const &OutFormat,
             std::string const &OutFile) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> GramMat = ReadMatrixFile<T>(MatFile);
  PVoronoi<T, Tint> pv = ReadEntryCPP_PVoronoi_File<T,Tint>(PVoronoiFile);

  int dimEXT = GramMat.rows() + 1;
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, std::cerr);
  DataLattice<T, Tint, Tgroup> eData =
      GetDataLattice<T, Tint, Tgroup>(GramMat, AllArr, std::cerr);
  std::vector<MyMatrix<Tint>> l_gens = get_p_voronoi_stabilizer(eData, pv);
  auto f_print = [&](std::ostream &os_out) -> void {
    if (OutFormat == "GAP") {
      os_out << "return ";
      WriteListMatrixGAP(os_out, l_gens);
      os_out << ";\n";
      return;
    }
    if (OutFormat == "CPP") {
      WriteListMatrix(os_out, l_gens);
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat\n";
    std::cerr << "Allowed choices: GAP, CPP\n";
    throw TerminalException{1};
  };
  print_stderr_stdout_file(OutFile, f_print);
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 6 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Robust_PVoronoi_Stabilizer [arith] [MatFile] [PVoronoiFile] [OutFormat] [OutFile]\n";
      std::cerr << "       or\n";
      std::cerr << "Robust_PVoronoi_Stabilizer [arith] [MatFile] [PVoronoiFile]\n";
      std::cerr << "allowed choices:\n";
      std::cerr << "arithmetic: gmp\n";
      std::cerr << "OutFormat: GAP\n";
      std::cerr << "OutFile: stderr, stdout, my_file\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    std::string MatFile = argv[2];
    std::string PVoronoiFile = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    auto f=[&]() -> void {
      if (arithmetic == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint>(MatFile, PVoronoiFile, OutFormat, OutFile);
      }
      std::cerr << "failure: No matching entry for arithmetic\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of Robust_PVoronoi_Stabilizer\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Robust_PVoronoi_Stabilizer\n";
    exit(e.eVal);
  }
  runtime(time);
}

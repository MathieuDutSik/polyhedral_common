// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "Permutation.h"
#include "Group.h"
#include "NumberTheory.h"
#include "enum_robust_covering.h"
// clang-format on

template <typename T, typename Tint>
void process_B(std::string const &MatFile,
               std::string const &VFile,
               std::string const &OutFormat,
               std::string const &OutFile) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> GramMat = ReadMatrixFile<T>(MatFile);
  auto get_p_voronoi=[&]() -> PVoronoi<T, Tint> {
    CVPSolver<T, Tint> solver(GramMat, std::cerr);
    int dimEXT = GramMat.rows() + 1;
    PolyHeuristicSerial<TintGroup> AllArr =
    AllStandardHeuristicSerial<T, TintGroup>(dimEXT, std::cerr);
    DataLattice<T, Tint, Tgroup> eData =
    GetDataLattice<T, Tint, Tgroup>(GramMat, AllArr, std::cerr);
    if (!IsExistingFile(VFile)) {
      return initial_p_polytope<T,Tint>(eData);
    }
    MyVector<T> eV = ReadVectorFile<T>(VFile);
    std::optional<PVoronoi<T, Tint>> opt = find_p_voronoi(eData, eV);
    if (!opt) {
      std::cerr << "ROBUST: The find_p_voronoi failed\n";
      throw TerminalException{1};
    }
    return *opt;
  };
  PVoronoi<T, Tint> p_voronoi = get_p_voronoi();
  auto f_print = [&](std::ostream &os_out) -> void {
    if (OutFormat == "GAP") {
      os_out << "return ";
      WriteEntryGAP(os_out, p_voronoi);
      os_out << ";\n";
      return;
    }
    if (OutFormat == "CPP") {
      WriteEntryCPP(os_out, p_voronoi);
      return;
    }
    if (OutFormat == "boost") {
      boost::archive::text_oarchive oa(os_out);
      oa << p_voronoi;
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat\n";
    std::cerr << "Allowed choices: GAP, CPP, boost\n";
    throw TerminalException{1};
  };
  print_stderr_stdout_file(OutFile, f_print);
}

void process_A(std::string const &arithmetic, std::string const &MatFile, std::string const &VFile,
               std::string const &OutFormat, std::string const &OutFile) {
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_B<T, Tint>(MatFile, VFile, OutFormat, OutFile);
  }
  std::cerr << "process_A failure: No matching entry for arithmetic_mat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 6 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Robust_GetPVoronoi [arith] [MatFile] [VFile] [OutFormat] [OutFile]\n";
      std::cerr << "       or\n";
      std::cerr << "Robust_GetPVoronoi [arith] [MatFile] [VFile]\n";
      std::cerr << "allowed choices:\n";
      std::cerr << "arithmetic: gmp\n";
      std::cerr << "OutFormat: GAP, CPP, boost\n";
      std::cerr << "OutFile: stderr, stdout, my_file\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    std::string MatFile = argv[2];
    std::string VFile = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    process_A(arithmetic, MatFile, VFile, OutFormat, OutFile);
    std::cerr << "Normal termination of Robust_GetPVoronoi\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Robust_GetPVoronoi\n";
    exit(e.eVal);
  }
  runtime(time);
}

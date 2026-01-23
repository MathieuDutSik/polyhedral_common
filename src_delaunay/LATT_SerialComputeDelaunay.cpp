// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeDelaunay.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void process(std::string const &GramFile, std::string const &OutFormat,
             std::ostream &os_out) {
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> GramMat = ReadMatrixFile<T>(GramFile);

  int dimEXT = GramMat.rows() + 1;
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, std::cerr);
  DataLattice<T, Tint, Tgroup> data_lattice =
      GetDataLattice<T, Tint, Tgroup>(GramMat, AllArr, std::cerr);
  auto f_incorrect =
      [&]([[maybe_unused]] Delaunay_Obj<T, Tgroup> const &x) -> bool {
    return false;
  };
  int max_runtime_second = 0;
  std::optional<DelaunayTesselation<T, Tgroup>> opt =
      EnumerationDelaunayPolytopes<T, Tint, Tgroup, decltype(f_incorrect)>(
          data_lattice, f_incorrect, max_runtime_second);
  DelaunayTesselation<T, Tgroup> DT =
      unfold_opt(opt, "The Delaunay tesselation");
  WriteDelaunayTesselation(OutFormat, os_out, GramMat, DT);
}

int main(int argc, char *argv[]) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_SerialComputeDelaunay [arith] [GramFile]\n";
      std::cerr << "     or\n";
      std::cerr << "LATT_SerialComputeDelaunay [arith] [GramFile] [OutFormat] "
                   "[OutFile]\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string GramFile = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint, Tgroup>(GramFile, OutFormat, os);
      }
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
      std::cerr << "Allowed arith: gmp\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(OutFile, f);
    std::cerr << "Normal termination of LATT_SerialComputeDelaunay\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_SerialComputeDelaunay\n";
    exit(e.eVal);
  }
  runtime(time);
}

// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "IsoDelaunayDomains.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint>
void process(std::string const& FileListMat, std::string const& OutFormat, std::ostream& os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::vector<MyMatrix<T>> ListMat = ReadListMatrixFile<T>(FileListMat);

  LinSpaceMatrix<T> LinSpa = BuildLinSpaceMatrix<T,Tint>(ListMat, std::cerr);
  //
  int dimEXT = LinSpa.n + 1;
  PolyHeuristicSerial<TintGroup> AllArr =
    AllStandardHeuristicSerial<T, TintGroup>(dimEXT, std::cerr);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, std::cerr);
  //
  std::optional<MyMatrix<T>> CommonGramMat;
  DataIsoDelaunayDomains<T, Tint, Tgroup> data{LinSpa, std::move(rddo),
                                               CommonGramMat};
  //
  using Tdata = DataIsoDelaunayDomainsFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_Serial<Tobj, TadjO>;
  int max_runtime_second = 0;
  auto f_incorrect = [&]([[maybe_unused]] Tobj const &x) -> bool {
    return false;
  };
  std::optional<std::vector<Tout>> opt = EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
      data_func, f_incorrect, max_runtime_second);
  if (!opt) {
    std::cerr << "Failed to terminate the enumeration, which is abnormal\n";
    throw TerminalException{1};
  }
  std::vector<Tout> const& l_tot = *opt;
  bool result = WriteFamilyObjects(data, OutFormat, os_out, l_tot, std::cerr);
  if (result) {
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
    throw TerminalException{1};
  }
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_SerialLattice_IsoDelaunayDomain [arith] [FileListMat]\n";
      std::cerr << "      or\n";
      std::cerr << "LATT_SerialLattice_IsoDelaunayDomain [arith] [FileListMat] [OutFormat] [OutFile]\n";
      std::cerr << "with arith the arithmetic type\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileListMat = argv[2];
    std::string OutFormat = "NumberGAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    auto f=[&](std::ostream& os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T,Tint>(FileListMat, OutFormat, os);
      }
      std::cerr << "Failed to find a matching entry for arith=" << arith << "\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(OutFile, f);
    std::cerr << "Normal termination of LATT_SerialLattice_IsoDelaunayDomain\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_SerialLattice_IsoDelaunayDomain\n";
    exit(e.eVal);
  }
  runtime(time);
}

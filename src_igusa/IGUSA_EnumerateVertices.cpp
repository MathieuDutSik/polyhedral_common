// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "igusa.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void process(FullNamelist const &eFull, std::ostream &os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint, Tgroup>(BlockTSPACE, os);
  IgusaParameters<T> params{BlockDATA.get_int("NormBound"),
                            BlockDATA.get_int("NbRandomFacet"),
                            BlockDATA.get_int("NbCleanRound"),
                            BlockDATA.get_string("IlpMethod")};
  IgusaSpace<T, Tint> space = build_igusa_space<T, Tint>(LinSpa, params, os);
  using Tdata = DataIgusaFunc<T, Tint, Tgroup>;
  Tdata data{std::move(space), os};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_Serial<Tobj, TadjO>;
  auto f_incorrect = [&]([[maybe_unused]] Tobj const &x) -> bool {
    return false;
  };
  int max_runtime_second = BlockDATA.get_int("max_runtime_second");
  std::optional<std::vector<Tout>> opt_l_tot =
      EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
          data, f_incorrect, max_runtime_second);
  std::vector<Tout> l_tot =
      unfold_opt(opt_l_tot, "EnumerateAndStore_Serial (Igusa vertices)");
  os << "IGUSA: number of orbits of vertices=" << l_tot.size() << "\n";
  //
  std::string OutFormat = BlockDATA.get_string("OutFormat");
  std::string OutFile = BlockDATA.get_string("OutFile");
  auto f_print = [&](std::ostream &os_out) -> void {
    auto write_list = [&](auto f_write) -> void {
      os_out << "[";
      bool IsFirst = true;
      for (auto &ent : l_tot) {
        if (!IsFirst) {
          os_out << ",\n";
        }
        IsFirst = false;
        f_write(ent);
      }
      os_out << "]";
    };
    if (OutFormat == "GAP") {
      os_out << "return ";
      write_list([&](Tout const &ent) -> void { WriteEntryGAP(os_out, ent); });
      os_out << ";\n";
      return;
    }
    if (OutFormat == "PYTHON") {
      write_list(
          [&](Tout const &ent) -> void { WriteEntryPYTHON(os_out, ent); });
      os_out << "\n";
      return;
    }
    std::cerr << "IGUSA: Unknown OutFormat=" << OutFormat << "\n";
    std::cerr << "IGUSA: Allowed formats: GAP, PYTHON\n";
    throw TerminalException{1};
  };
  if (OutFile == "stderr") {
    return f_print(std::cerr);
  }
  if (OutFile == "stdout") {
    return f_print(std::cout);
  }
  std::ofstream os_out(OutFile);
  f_print(os_out);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_ENUMERATE_IGUSA_TSPACE();
    if (argc != 2) {
      std::cerr << "IGUSA_EnumerateVertices [file.nml]\n";
      std::cerr << "with file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    SingleBlock const &BlockDATA = eFull.get_block("DATA");
    std::string arithmetic = BlockDATA.get_string("arithmetic");
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    if (arithmetic == "gmp") {
      using T = mpq_class;
      using Tint = mpz_class;
      using Tgroup = permutalib::Group<Telt, Tint>;
      process<T, Tint, Tgroup>(eFull, std::cerr);
    } else {
      std::cerr << "IGUSA: Unknown arithmetic=" << arithmetic << "\n";
      std::cerr << "IGUSA: Allowed arithmetic: gmp\n";
      throw TerminalException{1};
    }
    std::cerr << "Normal termination of IGUSA_EnumerateVertices\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in IGUSA_EnumerateVertices\n";
    exit(e.eVal);
  }
  runtime(time);
}

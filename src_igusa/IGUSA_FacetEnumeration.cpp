// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "igusa_facet_enum.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void process(FullNamelist const &eFull, std::ostream &os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint, Tgroup>(BlockTSPACE, os);
  int n = LinSpa.n;
  if (static_cast<int>(LinSpa.ListMat.size()) != n * (n + 1) / 2) {
    std::cerr << "IGUSA_FACET_ENUM: only the full space of symmetric "
                 "matrices (Classic) is supported\n";
    throw TerminalException{1};
  }
  IgusaParameters<T> params{};
  params.IlpMethod = BlockDATA.get_string("IlpMethod");
  IgusaSpace<T, Tint> space = build_igusa_space<T, Tint>(LinSpa, params, os);
  igusa_insert_initial_cuts(space);
  std::vector<MyMatrix<T>> ListF =
      ReadListMatrixFile<T>(BlockDATA.get_string("FileInitialFacets"));
  std::ifstream is_rhs(BlockDATA.get_string("FileInitialRhs"));
  std::vector<std::pair<MyMatrix<T>, T>> ListInit;
  for (auto &F : ListF) {
    std::string str;
    is_rhs >> str;
    ListInit.push_back({F, ParseScalar<T>(str)});
  }
  int max_orbit = BlockDATA.get_int("MaxOrbit");
  IgusaFacetEnumResult<T> l_orbit =
      igusa_facet_enumeration<T, Tint, Tgroup>(space, ListInit, max_orbit, os);
  os << "IGUSA_FACET_ENUM: number of orbits of full rank facets="
     << l_orbit.l_facet.size() << " of bounded ridges in facets of lower rank="
     << l_orbit.l_ridge.size() << "\n";
  std::string OutFile = BlockDATA.get_string("OutFile");
  if (OutFile == "stderr") {
    return WriteFacetOrbitsGAP(std::cerr, l_orbit);
  }
  if (OutFile == "stdout") {
    return WriteFacetOrbitsGAP(std::cout, l_orbit);
  }
  std::ofstream os_out(OutFile);
  WriteFacetOrbitsGAP(os_out, l_orbit);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_IGUSA_FACET_ENUMERATION();
    if (argc != 2) {
      std::cerr << "IGUSA_FacetEnumeration [file.nml]\n";
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
    auto f_process = [&]() -> void {
      if (arithmetic == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        using Tgroup = permutalib::Group<Telt, mpz_class>;
        return process<T, Tint, Tgroup>(eFull, std::cerr);
      }
#ifdef ENABLE_FLINT_SUPPORT
      if (arithmetic == "flint") {
        using T = fmpq_class;
        using Tint = fmpz_class;
        using Tgroup = permutalib::Group<Telt, mpz_class>;
        return process<T, Tint, Tgroup>(eFull, std::cerr);
      }
#endif
      std::cerr << "IGUSA_FACET_ENUM: Unknown arithmetic=" << arithmetic
                << "\n";
      throw TerminalException{1};
    };
    f_process();
    std::cerr << "Normal termination of IGUSA_FacetEnumeration\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in IGUSA_FacetEnumeration\n";
    exit(e.eVal);
  }
  runtime(time);
}

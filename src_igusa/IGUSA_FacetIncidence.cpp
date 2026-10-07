// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "igusa_facet.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void process(FullNamelist const &eFull, std::ostream &os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint, Tgroup>(BlockTSPACE, os);
  IgusaParameters<T> params{};
  params.IlpMethod = BlockDATA.get_string("IlpMethod");
  IgusaSpace<T, Tint> space = build_igusa_space<T, Tint>(LinSpa, params, os);
  igusa_insert_initial_cuts(space);
  std::string FileFacet = BlockDATA.get_string("FileFacet");
  MyMatrix<T> F = ReadMatrixFile<T>(FileFacet);
  T rhs = ParseScalar<T>(BlockDATA.get_string("FacetRhs"));
  IgusaFacetIncidence<T> inc = igusa_facet_incidence(space, F, rhs, os);
  os << "IGUSA_FACET: rhs=" << inc.rhs << " min over I=" << inc.igusa_min
     << " min over the Ryshkov polyhedron=" << inc.ryshkov_min
     << " slack=" << inc.slack
     << " |Inc|=" << inc.l_incident.size() << " rank=" << inc.rank
     << " (dim=" << space.dim << ")\n";
  std::string OutFile = BlockDATA.get_string("OutFile");
  auto f_print = [&](std::ostream &os_out) -> void {
    os_out << "return rec(F:=";
    WriteMatrixGAP(os_out, inc.F);
    os_out << ",\n rhs:=" << inc.rhs << ",\n RyshkovMin:=" << inc.ryshkov_min
           << ",\n slack:=" << inc.slack
           << ",\n RyshkovMinimizer:=";
    WriteMatrixGAP(os_out, inc.ryshkov_minimizer);
    os_out << ",\n rank:=" << inc.rank << ",\n ListIncident:=[";
    bool is_first = true;
    for (auto &X : inc.l_incident) {
      if (!is_first) {
        os_out << ",\n";
      }
      is_first = false;
      WriteMatrixGAP(os_out, X);
    }
    os_out << "]);\n";
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
    FullNamelist eFull = NAMELIST_GetStandard_IGUSA_FACET_INCIDENCE();
    if (argc != 2) {
      std::cerr << "IGUSA_FacetIncidence [file.nml]\n";
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
      std::cerr << "IGUSA_FACET: Unknown arithmetic=" << arithmetic << "\n";
      throw TerminalException{1};
    };
    f_process();
    std::cerr << "Normal termination of IGUSA_FacetIncidence\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in IGUSA_FacetIncidence\n";
    exit(e.eVal);
  }
  runtime(time);
}

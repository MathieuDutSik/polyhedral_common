// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "CtypeMPI_enumeration.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

int main(int argc, char *argv[]) {
  maybe_install_tls_gmp_pool();
  boost::mpi::environment env(boost::mpi::threading::serialized);
  if (env.thread_level() < boost::mpi::threading::serialized) {
    env.abort(-1);
  }
  boost::mpi::communicator world;
  HumanTime time;
  using T = mpq_class;
  // int is not enough: the canonicalization of the principal domain in
  // dimension 6 overflows Rational<int> arithmetic.
  using Tint = int64_t;
  using Tidx = uint16_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_COMPUTE_LATTICE_IsoEdgeDomains();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CTYP_MPI_AdjScheme [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    ComputeLatticeIsoEdgeDomains<T, Tint, Tgroup>(world, eFull);
    std::cerr << "Normal termination in CTYP_MPI_AdjScheme\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CTYP_MPI_AdjScheme\n";
    exit(e.eVal);
  }
  runtime(time);
}

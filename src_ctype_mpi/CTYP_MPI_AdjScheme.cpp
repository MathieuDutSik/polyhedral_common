// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "CtypeMPI_enumeration.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

int main(int argc, char *argv[]) {
  boost::mpi::environment env(boost::mpi::threading::serialized);
  if (env.thread_level() < boost::mpi::threading::serialized) {
    env.abort(-1);
  }
  boost::mpi::communicator world;
  HumanTime time;
  using T = mpq_class;
  using Tint = int;
  using Tidx = uint16_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  try {
    Eigen::initParallel();
    FullNamelist eFull = NAMELIST_GetStandard_COMPUTE_LATTICE_IsoEdgeDomains();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CTYP_MPI_AdjScheme [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    ComputeLatticeIsoEdgeDomains<T,Tint,Tgroup>(world, eFull);
  } catch (TerminalException const &e) {
    std::cerr << "Error in CTYP_MPI_AdjScheme\n";
    exit(e.eVal);
  }
  runtime(time);
}

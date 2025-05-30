// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeDelaunay_mpi.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint>
void process_C(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  return ComputeDelaunayPolytope_MPI<T, Tint, Tgroup>(comm, eFull);
}

template <typename T>
void process_B(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  std::string arithmetic_Tint =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic_Tint");
  if (arithmetic_Tint == "gmp_integer") {
    using Tint = mpz_class;
    return process_C<T, Tint>(comm, eFull);
  }
  std::cerr << "LATT_MPI_ComputeDelaunay B: Failed to find a matching type for "
               "arithmetic_Tint="
            << arithmetic_Tint << "\n";
  std::cerr << "Available types: gmp_integer\n";
  throw TerminalException{1};
}

void process_A(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  std::string arithmetic_T =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic_T");
  if (arithmetic_T == "gmp_rational") {
    using T = mpq_class;
    return process_B<T>(comm, eFull);
  }
  std::cerr << "LATT_MPI_ComputeDelaunay A: Failed to find a matching type for "
               "arithmetic_T="
            << arithmetic_T << "\n";
  std::cerr << "Available types: gmp_rational\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  boost::mpi::environment env(boost::mpi::threading::serialized);
  if (env.thread_level() < boost::mpi::threading::serialized) {
    env.abort(-1);
  }
  boost::mpi::communicator world;
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_COMPUTE_DELAUNAY();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_MPI_ComputeDelaunay [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    //
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_A(world, eFull);
    std::cerr << "Normal termination of LATT_MPI_ComputeDelaunay\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_MPI_ComputeDelaunay\n";
    exit(e.eVal);
  }
  runtime(time1);
}

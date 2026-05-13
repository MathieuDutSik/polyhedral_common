// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "perfect_tspace_mpi.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint>
void process_C(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  return ComputePerfectTspace_mpi<T, Tint, Tgroup>(comm, eFull);
}

void process_A(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  std::string arithmetic =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic");
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_C<T, Tint>(comm, eFull);
  }
  if (arithmetic == "gmp_boost") {
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
    return process_C<T, Tint>(comm, eFull);
  }
  if (arithmetic == "multi_boost") {
    using T = boost::multiprecision::cpp_rational;
    using Tint = boost::multiprecision::cpp_int;
    return process_C<T, Tint>(comm, eFull);
  }
  std::cerr
      << "PERF_MPI_EnumeratePerfectCones: Failed to find matching type for "
      << "arithmetic=" << arithmetic << "\n";
  std::cerr << "Available types: gmp, gmp_boost, multi_boost\n";
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
    FullNamelist eFull = NAMELIST_GetStandard_ENUMERATE_PERFECT_TSPACE();

    if (argc != 2) {
      std::cerr << "Number of arguments = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "PERF_MPI_EnumeratePerfectCones [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }

    std::string eFileName = argv[1];
    std::cerr << "Reading namelist file: " << eFileName << "\n";
    NAMELIST_ReadNamelistFile(eFileName, eFull);

    process_A(world, eFull);
    std::cerr << "Normal termination of PERF_MPI_EnumeratePerfectCones\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PERF_MPI_EnumeratePerfectCones\n";
    exit(e.eVal);
  }
  runtime(time1);
}

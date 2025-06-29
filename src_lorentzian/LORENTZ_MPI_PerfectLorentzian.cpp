// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "lorentzian_perfect_mpi.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint>
void process_C(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  return ComputePerfectLorentzian<T, Tint, Tgroup>(comm, eFull);
}

template <typename T>
void process_B(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  std::string arithmetic_Tint =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic_Tint");
  if (arithmetic_Tint == "gmp_integer") {
    using Tint = mpz_class;
    return process_C<T, Tint>(comm, eFull);
  }
  std::cerr << "LORENTZ_MPI_PerfectLorentzian B: Failed to find a matching "
               "type for arithmetic_Tint="
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
  std::cerr << "LORENTZ_MPI_PerfectLorentzian A: Failed to find a matching "
               "type for arithmetic_T="
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
    std::cerr << "LORENTZ_MPI_PerfectLorentzian, step 1\n";
    FullNamelist eFull = NAMELIST_GetStandard_COMPUTE_PERFECT_LORENTZIAN();
    SingleBlock const& BlockDATA = eFull.get_block("DATA");
    int max_runtime_second = BlockDATA.get_int("max_runtime_second");
    std::cerr << "LORENTZ_MPI_PerfectLorentzian, step 1, max_runtime_second=" << max_runtime_second << "\n";
    std::cerr << "LORENTZ_MPI_PerfectLorentzian, step 2\n";
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LORENTZ_MPI_PerfectLorentzian [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    //
    std::string eFileName = argv[1];
    std::cerr << "LORENTZ_MPI_PerfectLorentzian, step 3\n";
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    std::cerr << "LORENTZ_MPI_PerfectLorentzian, step 4\n";
    process_A(world, eFull);
    std::cerr << "Normal termination of LORENTZ_MPI_PerfectLorentzian\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LORENTZ_MPI_PerfectLorentzian\n";
    exit(e.eVal);
  }
  runtime(time1);
}

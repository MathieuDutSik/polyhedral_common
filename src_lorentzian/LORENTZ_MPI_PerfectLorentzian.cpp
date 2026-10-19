// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "lorentzian_perfect_mpi.h"
#include "Permutation.h"
#include "Group.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T, typename Tint>
void process_C(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  return ComputePerfectLorentzian_mpi<T, Tint, Tgroup>(comm, eFull);
}

void process_A(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  std::string arithmetic =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic");
  ArithmeticDispatch dispatch;
  add_rational_integer_arithmetics_by_library(
      dispatch, [&]<typename T, typename Tint>() -> void {
        process_C<T, Tint>(comm, eFull);
      });
  dispatch.run(arithmetic);
}

int main(int argc, char *argv[]) {
  maybe_install_tls_gmp_pool();
  boost::mpi::environment env(boost::mpi::threading::serialized);
  if (env.thread_level() < boost::mpi::threading::serialized) {
    env.abort(-1);
  }
  boost::mpi::communicator world;
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_COMPUTE_PERFECT_LORENTZIAN();
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
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_A(world, eFull);
    std::cerr << "Normal termination of LORENTZ_MPI_PerfectLorentzian\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LORENTZ_MPI_PerfectLorentzian\n";
    exit(e.eVal);
  }
  runtime(time1);
}

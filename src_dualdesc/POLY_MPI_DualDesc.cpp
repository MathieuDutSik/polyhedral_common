// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "Permutation.h"
#include "Group.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "POLY_RecursiveDualDesc_MPI.h"
// clang-format on

template <typename T, typename Tidx>
void Process_eFull(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint>;
  //    using Tidx_value = int16_t;
  using Tidx_value = int32_t;
  MPI_MainFunctionDualDesc<T, Tgroup, Tidx_value>(comm, eFull);
}

int main(int argc, char *argv[]) {
  // The construction is relatively subtle.
  // ---We need to have the env and comm on top.
  // ---This is because throwing an exception that goes
  // outside of the lifetime of world creates a MPI_ABORT.
  // ---We should also avoid using exit when terminating
  // because it makes an unscheduled destruction of the communicator.
  boost::mpi::environment env;
  boost::mpi::communicator world;
  SingletonTime start;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_RecursiveDualDescription();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_MPI_DualDesc [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    using T = mpq_class;
    //    using T = boost::multiprecision::cpp_rational;
    //    using T = boost::multiprecision::mpq_rational;
    MyMatrix<T> EXT = GetEXT_from_efull<T>(eFull);
    //
    auto process = [&]() -> void {
      if (size_t(EXT.rows()) < std::numeric_limits<uint8_t>::max())
        return Process_eFull<T, uint8_t>(world, eFull);
      if (size_t(EXT.rows()) < std::numeric_limits<uint16_t>::max())
        return Process_eFull<T, uint16_t>(world, eFull);
      if (size_t(EXT.rows()) < std::numeric_limits<uint32_t>::max())
        return Process_eFull<T, uint32_t>(world, eFull);
#if !defined __APPLE__
      if (size_t(EXT.rows()) < std::numeric_limits<uint64_t>::max())
        return Process_eFull<T, uint64_t>(world, eFull);
#endif
      std::cerr << "Failed to find a numeric type that matches\n";
      throw TerminalException{1};
    };
    process();
    //
    std::cerr << "Normal termination of the program runtime=" << si(start)
              << "\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program runtime=" << si(start)
              << "\n";
    exit(e.eVal);
  } catch (RuntimeException const &e) {
    std::cerr << "Runtime termination of the program runtime=" << si(start)
              << "\n";
  }
}

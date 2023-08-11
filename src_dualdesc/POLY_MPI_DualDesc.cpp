// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "Permutation.h"
#include "Group.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheorySafeInt.h"
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

template<typename T>
void Process_eFull_select_type(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  MyMatrix<T> EXT = GetEXT_from_efull<T>(eFull);
  //
  if (size_t(EXT.rows()) < std::numeric_limits<uint8_t>::max())
    return Process_eFull<T, uint8_t>(comm, eFull);
  if (size_t(EXT.rows()) < std::numeric_limits<uint16_t>::max())
    return Process_eFull<T, uint16_t>(comm, eFull);
  if (size_t(EXT.rows()) < std::numeric_limits<uint32_t>::max())
    return Process_eFull<T, uint32_t>(comm, eFull);
#if !defined __APPLE__
  if (size_t(EXT.rows()) < std::numeric_limits<uint64_t>::max())
    return Process_eFull<T, uint64_t>(comm, eFull);
#endif
  std::cerr << "Failed to find a numeric type that matches\n";
  throw TerminalException{1};
}



int main(int argc, char *argv[]) {
  // The construction is relatively subtle.
  // ---We need to have the env and comm on top.
  // ---This is because throwing an exception that goes
  // outside of the lifetime of world creates a MPI_ABORT.
  // ---We should also avoid using exit when terminating
  // because it makes an unscheduled destruction of the communicator.
  // ---Mpi is not thread-safe by default, we pass the threading::serialized
  // opton to allow mpi calls from all threads (but not at the same time).
  // this is used by the communication threads.
  boost::mpi::environment env(boost::mpi::threading::serialized);
  if (env.thread_level() < boost::mpi::threading::serialized) {
    env.abort(-1);
  }
  boost::mpi::communicator world;
  HumanTime start;
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
    std::string NumericalType = GetNumericalType(eFull);
    auto process=[&]() -> void {
      if (NumericalType == "rational") {
        using T = mpq_class;
        return Process_eFull_select_type<T>(world, eFull);
      }
      if (NumericalType == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return Process_eFull_select_type<T>(world, eFull);
      }
      if (NumericalType == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return Process_eFull_select_type<T>(world, eFull);
      }
      if (NumericalType == "RealAlgebraic") {
        using T_rat = mpq_class;
        SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
        std::string FileAlgebraicField =
          BlockDATA.ListStringValues.at("FileAlgebraicField");
        if (!IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                    << " is missing\n";
          throw TerminalException{1};
        }
        HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        using T = RealField<idx_real_algebraic_field>;
        return Process_eFull_select_type<T>(world, eFull);
      }
      std::cerr << "Failed to find a matching type entry\n";
      throw TerminalException{1};
    };
    process();
    //
    //    using T = boost::multiprecision::cpp_rational;
    //    using T = boost::multiprecision::mpq_rational;
    //
    std::cerr << "Normal termination of the program runtime=" << start
              << "\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program runtime=" << start
              << "\n";
    exit(e.eVal);
  } catch (RuntimeException const &e) {
    std::cerr << "Runtime termination of the program runtime=" << start
              << "\n";
  }
}

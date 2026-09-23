// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "lorentzian_perfect_serial.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint>
void process_B(FullNamelist const &eFull) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  return ComputePerfectLorentzian_serial<T, Tint, Tgroup>(eFull, std::cerr);
}

void process_A(FullNamelist const &eFull) {
  std::string arithmetic =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic");
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_B<T, Tint>(eFull);
  }
#ifdef ENABLE_BOOST_TYPES
  if (arithmetic == "gmp_boost") {
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
    return process_B<T, Tint>(eFull);
  }
  if (arithmetic == "multi_boost") {
    using T = boost::multiprecision::cpp_rational;
    using Tint = boost::multiprecision::cpp_int;
    return process_B<T, Tint>(eFull);
  }
#endif
#ifdef ENABLE_FLINT_SUPPORT
  if (arithmetic == "flint") {
    using T = fmpq_class;
    using Tint = fmpz_class;
    return process_B<T, Tint>(eFull);
  }
#endif
  std::cerr << "LORENTZ_Serial_PerfectLorentzian: Failed to find a matching "
               "type for arithmetic="
            << arithmetic << "\n";
#ifdef ENABLE_FLINT_SUPPORT
  std::cerr << "Available types: gmp, gmp_boost, multi_boost, flint\n";
#else
  std::cerr << "Available types: gmp, gmp_boost, multi_boost (build "
            << "with ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_COMPUTE_PERFECT_LORENTZIAN();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LORENTZ_Serial_PerfectLorentzian [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    //
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_A(eFull);
    std::cerr << "Normal termination of LORENTZ_Serial_PerfectLorentzian\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LORENTZ_Serial_PerfectLorentzian\n";
    exit(e.eVal);
  }
  runtime(time1);
}

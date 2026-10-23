// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Shvec_exact.h"
#include "LatticeStabEquiCan.h"
#include "OrbitsVectorBasis.h"
#include "Permutation.h"
#include "Group.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename Tgroup, typename T, typename Tint>
void compute_orbit_basis(std::string const &FileM, std::string const &OutFormat,
                         std::ostream &os) {
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileM);
  MyMatrix<Tint> SHV = T_ShortestVector<T, Tint>(GramMat, std::cerr).SHV;
  std::vector<MyMatrix<Tint>> ListGen =
      ArithmeticAutomorphismGroup<T, Tint, Tgroup>(GramMat, std::cerr);
  vectface vf = EnumerateOrbitBasis<Tgroup, Tint>(SHV, ListGen, std::cerr);
  if (OutFormat == "GAP") {
    os << "return rec(SHV:=";
    WriteMatrixGAP(os, SHV);
    os << ", vf:=";
    WriteListFaceGAP(os, vf);
    os << ");\n";
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_ComputeShortestOrbits [arith] [FileM]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_ComputeShortestOrbits [arith] [FileM] "
                   "[OutFormat] [OutFille]\n";
      std::cerr << "\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "arith: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << ", flint\n";
#else
      std::cerr << "arith: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << " (build with ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
      std::cerr << "FileM: The Gram matrix on input\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string FileM = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint_grp = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint_grp>;
    auto f = [&](std::ostream &os) -> void {
      ArithmeticDispatch dispatch;
      add_rational_integer_arithmetics_by_library(
          dispatch, [&]<typename T, typename Tint>() -> void {
            compute_orbit_basis<Tgroup, T, Tint>(FileM, OutFormat, os);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of LATT_ComputeShortestOrbits\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_ComputeShortestOrbits\n";
    exit(e.eVal);
  }
  runtime(time);
}

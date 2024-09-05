// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Shvec_exact.h"
#include "LatticeStabEquiCan.h"
#include "OrbitsVectorBasis.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template<typename Tgroup, typename T, typename Tint>
void compute_orbit_basis(std::string const& FileM, std::string const& OutFormat, std::ostream& os) {
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileM);
  T_shvec_info<T, Tint> info = computeMinimum_GramMat<T, Tint>(GramMat);
  MyMatrix<Tint> SHV = MatrixFromVectorFamily(info.short_vectors);
  std::vector<MyMatrix<Tint>> ListGen = ArithmeticAutomorphismGroup<T,Tint>(GramMat, std::cerr);
  vectface vf = EnumerateOrbitBasis<Tgroup,Tint>(SHV, ListGen, std::cerr);
  if (OutFormat == "GAP") {
    os << "return rec(SHV:=";
    WriteMatrixGAP(os, SHV);
    os << ", vf:=";
    WriteListFaceGAP(os, vf);
    os << ");\n";
    return ;
  }
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_ComputeShortestOrbitsBasis [arith] [FileM]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_ComputeShortestOrbitsBasis [arith] [FileM] [OutFormat] [OutFille]\n";
      std::cerr << "\n";
      std::cerr << "arith: gmp\n";
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
    auto f=[&](std::ostream& os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return compute_orbit_basis<Tgroup,T,Tint>(FileM, OutFormat, os);
      }
      std::cerr << "Failed to find a matching entry for arith=" << arith << "\n";
      throw TerminalException{1};
    };
    if (OutFile == "stderr") {
      f(std::cerr);
    } else {
      if (OutFile == "stdout") {
        f(std::cout);
      } else {
        std::ofstream os(OutFile);
        f(os);
      }
    }
    std::cerr << "Normal termination of LATT_ComputeShortestOrbitsBasis\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_ComputeShortestOrbitsBasis\n";
    exit(e.eVal);
  }
  runtime(time);
}

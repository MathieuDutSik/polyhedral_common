// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "GenusEnumeration.h"
#include "SignatureSymmetric.h"
#include "Arithmetic_dispatch_types.h"
#include <string>
#include <vector>
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void ProcessGenus(FullNamelist const &eFull, std::ostream &os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  std::string FileLattice = BlockDATA.get_string("FileLattice");
  std::string FileMass = BlockDATA.get_string("FileMass");
  std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  GenusSpec spec = ReadGenusSpec(eFull);
  std::vector<MyMatrix<T>> ListSeed = ReadListMatrixFile<T>(FileLattice);
  if (ListSeed.empty()) {
    std::cerr << "GENUS_Enumerate: no seed lattice in " << FileLattice << "\n";
    throw TerminalException{1};
  }
  // A sum of 1 / |Aut(L)|, so a fraction over the type the groups count
  // with, which has nothing to do with the type of the Gram matrices.
  using Tmass = typename GenusEnumerationResult<T, Tgroup>::Tmass;
  Tmass TotalMass = ReadMassFile<Tmass>(FileMass);
  // The seeds of one genus have the same rank and the same determinant.
  int rank = ListSeed[0].rows();
  T det = DeterminantMat(ListSeed[0]);
  for (auto &eSeed : ListSeed) {
    if (!IsSymmetricMatrix(eSeed) || !IsPositiveDefinite(eSeed, std::cerr)) {
      std::cerr << "GENUS_Enumerate: a seed Gram matrix in " << FileLattice
                << " is not symmetric positive definite\n";
      throw TerminalException{1};
    }
    if (eSeed.rows() != rank) {
      std::cerr << "GENUS_Enumerate: the seeds in " << FileLattice
                << " have the ranks " << rank << " and " << eSeed.rows()
                << "\n";
      throw TerminalException{1};
    }
    if (DeterminantMat(eSeed) != det) {
      std::cerr << "GENUS_Enumerate: the seeds in " << FileLattice
                << " have the determinants " << det << " and "
                << DeterminantMat(eSeed) << "\n";
      throw TerminalException{1};
    }
  }
  int prime = spec.prime;
  if (spec.strategy == NeighborStrategy::OutPrime && prime == 0) {
    prime = ChooseNeighborPrime<T>(det);
  }
  GenusEnumerationResult<T, Tgroup> result = GenusEnumeration<T, Tint, Tgroup>(
      ListSeed, TotalMass, prime, spec.methods, spec.strategy, spec.rw,
      spec.hecke, spec.order, os);
  // With the iterate and random_walk strategies there is no single prime; 0
  // marks that in the output record.
  int prime_out = (spec.strategy == NeighborStrategy::OutPrime) ? prime : 0;
  auto f = [&](std::ostream &os_out) -> void {
    WriteGenusEnumerationResult(os_out, result, OutFormat, prime_out);
  };
  FILE_PrintStderrStdoutFile(OutFile, f);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    // The group orders are counted in TintGroup and the mass lives in
    // overlying_field<TintGroup>. Both are unrelated to the arithmetic of
    // the lattice coordinates, so Tgroup is fixed once here while the
    // coordinate types T / Tint are chosen per arithmetic below.
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    //
    FullNamelist eFull = NAMELIST_GetStandard_GENUS_ENUMERATE();
    ArithmeticDispatch dispatch;
    add_rational_integer_arithmetics_by_library(
        dispatch, [&]<typename T, typename Tint>() -> void {
          ProcessGenus<T, Tint, Tgroup>(eFull, std::cerr);
        });
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GENUS_Enumerate [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      std::cerr << "arithmetic values:\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    std::string arith = GetNamelistStringEntry(eFull, "DATA", "arithmetic");
    dispatch.run(arith);
    std::cerr << "Normal termination of GENUS_Enumerate\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_Enumerate\n";
    exit(e.eVal);
  }
  runtime(time);
}

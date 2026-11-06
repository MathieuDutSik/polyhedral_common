// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "FiniteGroupNormalizer.h"
#include "Permutation.h"
#include "Group.h"
#include <string>
#include <vector>
// clang-format on

/*
  Whether two finite groups G1, G2 of GL_n(Z) are conjugate in GL_n(Z)
  (Z-equivalent), and a matrix X with X G1 X^{-1} = G2 if they are, from
  the perfect forms of the spaces of invariant forms (see
  FiniteGroupNormalizer.h).
 */

template <typename T, typename Tint, typename Tgroup>
void process(std::string const &FileGens1, std::string const &FileGens2,
             std::ostream &os) {
  std::vector<MyMatrix<Tint>> gens1 = UniversalStdVectorMatrixConversion<Tint, T>(
      ReadListMatrixFile<T>(FileGens1));
  std::vector<MyMatrix<Tint>> gens2 = UniversalStdVectorMatrixConversion<Tint, T>(
      ReadListMatrixFile<T>(FileGens2));
  FiniteGroupNormalizerBounds bounds;
  FiniteGroupConjugacyResult<Tint> res =
      FiniteGroupConjugacy<T, Tint, Tgroup>(gens1, gens2, bounds, std::cerr);
  if (res.conjugator) {
    os << "return rec(conjugate:=true, conjugator:=";
    WriteMatrixGAP(os, *res.conjugator);
    os << ");\n";
  } else {
    os << "return rec(conjugate:=false, reason:=\"" << res.reason << "\");\n";
  }
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    if (argc != 3 && argc != 4) {
      std::cerr << "PERF_FiniteGroupConjugacy [FileGens1] [FileGens2] "
                   "[OutFile]\n";
      std::cerr << "\n";
      std::cerr << "FileGens1, FileGens2 : the generators of two finite\n";
      std::cerr << "           groups G1, G2 of GL_n(Z), in the ListMatrix\n";
      std::cerr << "           format\n";
      std::cerr << "OutFile  : stderr (default), stdout or a file; a GAP\n";
      std::cerr << "           record: conjugate, and the conjugator X with\n";
      std::cerr << "           X G1 X^{-1} = G2 when there is one\n";
      return -1;
    }
    std::string FileGens1 = argv[1];
    std::string FileGens2 = argv[2];
    std::string OutFile = (argc == 4) ? argv[3] : "stderr";
    auto f = [&](std::ostream &os) -> void {
      process<T, Tint, Tgroup>(FileGens1, FileGens2, os);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of PERF_FiniteGroupConjugacy\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PERF_FiniteGroupConjugacy\n";
    exit(e.eVal);
  }
  runtime(time);
}

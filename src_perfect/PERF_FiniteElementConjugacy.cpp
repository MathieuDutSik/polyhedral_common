// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "FiniteGroupNormalizer.h"
#include "Permutation.h"
#include "Group.h"
#include <string>
// clang-format on

/*
  Whether two matrices g, h of finite order of GL_n(Z) are conjugate in
  GL_n(Z), and a matrix X with X g X^{-1} = h if they are: the cyclic
  groups by the perfect forms of their spaces of invariant forms, then the
  orbit of the element under the normalizer (see FiniteGroupNormalizer.h).
 */

template <typename T, typename Tint, typename Tgroup>
void process(std::string const &FileG, std::string const &FileH,
             std::ostream &os) {
  MyMatrix<Tint> g = UniversalMatrixConversion<Tint, T>(ReadMatrixFile<T>(FileG));
  MyMatrix<Tint> h = UniversalMatrixConversion<Tint, T>(ReadMatrixFile<T>(FileH));
  FiniteGroupNormalizerBounds bounds;
  FiniteGroupConjugacyResult<Tint> res =
      FiniteElementConjugacy<T, Tint, Tgroup>(g, h, bounds, std::cerr);
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
      std::cerr << "PERF_FiniteElementConjugacy [FileG] [FileH] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "FileG, FileH : two matrices g, h of finite order of\n";
      std::cerr << "           GL_n(Z), in the Matrix format\n";
      std::cerr << "OutFile  : stderr (default), stdout or a file; a GAP\n";
      std::cerr << "           record: conjugate, and the conjugator X with\n";
      std::cerr << "           X g X^{-1} = h when there is one\n";
      return -1;
    }
    std::string FileG = argv[1];
    std::string FileH = argv[2];
    std::string OutFile = (argc == 4) ? argv[3] : "stderr";
    auto f = [&](std::ostream &os) -> void {
      process<T, Tint, Tgroup>(FileG, FileH, os);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of PERF_FiniteElementConjugacy\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PERF_FiniteElementConjugacy\n";
    exit(e.eVal);
  }
  runtime(time);
}

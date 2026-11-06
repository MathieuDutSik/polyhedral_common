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
  The normalizer in GL_n(Z) of a finite group G of GL_n(Z), given by
  generators acting on the forms by M -> g M g^T, from the G-perfect forms
  of the space of invariant forms (see FiniteGroupNormalizer.h).
 */

template <typename T, typename Tint, typename Tgroup>
void process(std::string const &FileGens, std::ostream &os) {
  std::vector<MyMatrix<T>> gens_T = ReadListMatrixFile<T>(FileGens);
  std::vector<MyMatrix<Tint>> gens =
      UniversalStdVectorMatrixConversion<Tint, T>(gens_T);
  FiniteGroupNormalizerBounds bounds;
  FiniteGroupNormalizerResult<Tint> res =
      FiniteGroupNormalizer<T, Tint, Tgroup>(gens, bounds, std::cerr);
  os << "return rec(n:=" << res.n << ", order_G:=" << res.order_G
     << ", order_Bravais:=" << res.order_B << ", dim_forms:=" << res.dim_space
     << ", n_perfect:=" << res.n_perfect
     << ", index_in_normalizer_Bravais:=" << res.index_in_NB
     << ",\n  generators:=";
  WriteListMatrixGAP(os, res.gens);
  os << ");\n";
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
    if (argc != 2 && argc != 3) {
      std::cerr << "PERF_FiniteGroupNormalizer [FileGens] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "FileGens : the generators of a finite group G of GL_n(Z),\n";
      std::cerr << "           in the ListMatrix format\n";
      std::cerr << "OutFile  : stderr (default), stdout or a file; a GAP\n";
      std::cerr << "           record with the generators of the normalizer\n";
      std::cerr << "           of G in GL_n(Z)\n";
      return -1;
    }
    std::string FileGens = argv[1];
    std::string OutFile = (argc == 3) ? argv[2] : "stderr";
    auto f = [&](std::ostream &os) -> void {
      process<T, Tint, Tgroup>(FileGens, os);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of PERF_FiniteGroupNormalizer\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PERF_FiniteGroupNormalizer\n";
    exit(e.eVal);
  }
  runtime(time);
}

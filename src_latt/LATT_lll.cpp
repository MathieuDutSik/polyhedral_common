// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "LatticeReduction.h"
// clang-format on

static void PrintUsage() {
  std::cerr << R"(
LATT_lll [arith] [method] [FileI] [OutFormat] [FileO]
or
LATT_lll [arith] [method] [FileI]

Reduces a positive definite Gram matrix, returning the reduced form and the
unimodular transformation that produced it.

 ------- arith --------

  gmp          mpq_class over mpz_class
)";
#ifdef ENABLE_BOOST_TYPES
  std::cerr << R"(  gmp_boost    boost mpq_rational over mpz_int
  multi_boost  boost cpp_rational over cpp_int
)";
#endif
#ifdef ENABLE_FLINT_SUPPORT
  std::cerr << R"(  flint        fmpq_class over fmpz_class
)";
#endif
  std::cerr << R"(
 ------- method --------

  direct       LLL
  dual         LLL through the adjugate of the form. On every measure applied
               in this package this is the weaker of the two LLL variants,
               sometimes by orders of magnitude
  seysen       Seysen, minimising sum_i |b_i|^2 |dual b_i|^2. Treats every
               index symmetrically and reduces the form and its dual in one
               descent; the one to use when what matters is that every
               coefficient be small rather than that one vector be short
  seysen_best  Seysen taking the best move of each sweep. Slower, and in
               measurement no better; its trajectory does not depend on the
               order pairs are visited in, which makes it a diagnostic
  seysen_lll   Seysen and LLL alternated while the measure improves
  deep         Schnorr-Euchner deep insertion: LLL's move set widened from
               the adjacent swap to an insertion at any earlier position
  deep5        deep insertion restricted to depth 5
  deep10       deep insertion restricted to depth 10
  bkz4         BKZ at block size 4
  bkz8         BKZ at block size 8
  bkz12        BKZ at block size 12
  slide4       Gama-Nguyen slide reduction, block size at most 4
  slide8       slide reduction, block size at most 8
  best         run all of the above and keep whichever minimises the squared
               orthogonality defect prod_i G_ii / det G, ties broken by the
               sum of absolute entries. The input is among the candidates, so
               the result is never worse than it

The block methods bkz and slide ask at each index for a shortest vector of a
projected block, so they cost more than the rest, superexponentially in the
block size, and the quality they buy improves slowly; above block size eight
the gains measured in this package were negligible. Slide reduction differs
from BKZ in using two families of conditions on non-overlapping blocks, which
buys a polynomial bound on the number of steps that BKZ has no analogue of.

Note that "best" is a choice and not a canonical notion. It minimises the
orthogonality defect, which is scale free and equals one exactly for an
orthogonal basis, and breaks ties on the size of the integers. A caller
wanting a short first vector, a balanced Gram-Schmidt profile or a small
Seysen measure should name the method instead, because the best by one of
those need not be the best by another.

 ------- OutFormat --------

  GAP          the reduced Gram matrix and the transformation, GAP readable
  CPP_G        just the reduced Gram matrix
  CPP_P        just the transformation matrix
               Default: GAP

 ------- FileO --------

  stdout       write to std::cout
  stderr       write to std::cerr
  anything     else is taken as a file name to write to
               Default: stderr
)";
}

/*
  The measures of the form, printed before and after so that the pair can be
  read at a glance. The orthogonality defect is scale free and equals one
  exactly for an orthogonal basis; the L1 norm is the size of the integers
  everything downstream will compute with.
 */
template <typename T>
static void PrintMeasures(std::string const &label, MyMatrix<T> const &G) {
  LatticeReductionQuality<T> q = ComputeLatticeReductionQuality(G);
  std::cerr << label << "\n";
  std::cerr << "  defect^2           : " << q.orth_defect_sq << "\n";
  std::cerr << "  L1 norm            : " << q.l1_norm << "\n";
  std::cerr << "  Linf norm          : " << Linfinity_norm_mat(G) << "\n";
}

template <typename T, typename Tint>
void process(std::string const &FileI, std::string const &method,
             std::string const &OutFormat, std::ostream &os) {
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileI);
  std::cerr << "input: dimension " << GramMat.rows() << ", method=" << method
            << "\n";
  PrintMeasures("before:", GramMat);
  LatticeReductionResult<T, Tint> res =
      LatticeReducedByName<T, Tint>(GramMat, method, std::cerr);
  PrintMeasures("after:", res.GramMatRed);
  std::cerr << "method used: " << res.method << "\n";
  if (OutFormat == "GAP") {
    os << "return rec(GramMat:=";
    WriteMatrixGAP(os, res.GramMatRed);
    os << ", Pmat:=";
    WriteMatrixGAP(os, res.Pmat);
    os << ", method:=\"" << res.method << "\");\n";
    return;
  }
  if (OutFormat == "CPP_G") {
    WriteMatrix(os, res.GramMatRed);
    return;
  }
  if (OutFormat == "CPP_P") {
    WriteMatrix(os, res.Pmat);
    return;
  }
  std::cerr << "Failed to find a matching OutFormat=" << OutFormat
            << ". Allowed: GAP, CPP_G, CPP_P\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 4 && argc != 6) {
      PrintUsage();
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string method = argv[2];
    std::string FileI = argv[3];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      FileO = argv[5];
    }
    //
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint>(FileI, method, OutFormat, os);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return process<T, Tint>(FileI, method, OutFormat, os);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return process<T, Tint>(FileI, method, OutFormat, os);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        using Tint = fmpz_class;
        return process<T, Tint>(FileI, method, OutFormat, os);
      }
#endif
      std::cerr << "Failed to find a matching type for arith=" << arith
                << "\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Allowed: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << ", flint\n";
#else
      std::cerr << "Allowed: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << " (build with ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(FileO, f);
    std::cerr << "Normal termination of LATT_lll\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_lll\n";
    exit(e.eVal);
  }
  runtime(time);
}

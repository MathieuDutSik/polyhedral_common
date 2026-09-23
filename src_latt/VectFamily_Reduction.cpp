// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "VectFamilyReduction.h"
#include "norms.h"
// clang-format on

static void PrintUsage() {
  std::cerr << R"(
VectFamily_Reduction [FileI] [method] [OutFormat] [FileO]
or
VectFamily_Reduction [FileI] [method]

Reduces a family of vectors: a change of coordinates in the ambient space
that makes the coefficients of the family small. This is what is applied to
an EXT matrix before a dual description.

 ------- method --------

What is reduced is the Gram matrix of the COLUMNS. Note that this is not the
same problem as reducing a lattice basis: what a consumer of the output pays
for is the size of EVERY coefficient, and for the dual description the
quantity that governs it is the Hadamard estimate on the facet coefficients,
which depends on all the vectors symmetrically. A reduction aiming at one
short vector, which is what LLL does, is therefore not obviously the right
choice here, and in measurement the winner varies by instance.

  direct       LLL
  dual         LLL through the adjugate of the form
  seysen       Seysen, minimising sum_i |b_i|^2 |dual b_i|^2. Treats every
               coordinate symmetrically and reduces the form and its dual in
               one descent; the natural first thing to try here
  seysen_best  Seysen taking the best move of each sweep. Slower, and in
               measurement no better; useful as a diagnostic, its trajectory
               not depending on the order pairs are visited in
  seysen_lll   Seysen and LLL alternated while the measure improves
  deep         Schnorr-Euchner deep insertion, unrestricted depth
  deep5        deep insertion restricted to depth 5
  deep10       deep insertion restricted to depth 10
  bkz4         BKZ at block size 4
  bkz8         BKZ at block size 8
  bkz12        BKZ at block size 12
  slide4       Gama-Nguyen slide reduction, block size at most 4
  slide8       slide reduction, block size at most 8
  best         run all of the above and keep whichever actually minimises the
               facet coefficient estimate. The unreduced input is among the
               candidates, so the result is never worse than what was handed
               in. Recommended unless the cost of the reduction itself
               matters: it is a small multiple of one reduction, and that is
               nothing against the dual description that follows

 ------- OutFormat --------

  GAP          write a GAP readable file
  CPP          write a CPP polyhedral readable file

 ------- FileO --------

  stderr       write to std::cerr
  stdout       write to std::cout
  anything     else is taken as a file name to write to
)";
}

/*
  The measures, aligned in one column so that a before and after pair can be
  read at a glance. The facet estimate is a SQUARED quantity, so its square
  root is what should be compared against the range of the integer type the
  dual description will run in; that comparison is the reason the estimate is
  printed at all.
 */
template <typename T>
static void PrintMeasures(std::string const &label, MyMatrix<T> const &M) {
  T facet_sqr = sqr_estimate_facet_coefficients(M);
  double facet_d = UniversalScalarConversion<double, T>(facet_sqr);
  double int64_max =
      UniversalScalarConversion<double, int64_t>(
          std::numeric_limits<int64_t>::max());
  std::cerr << label << "\n";
  std::cerr << "  facet estimate (squared) : " << facet_sqr << "\n";
  std::cerr << "  facet estimate           : " << sqrt(facet_d)
            << "   (int64 max " << int64_max << ")\n";
  std::cerr << "  L1 norm                  : " << L1_norm_mat(M) << "\n";
  std::cerr << "  Linf norm                : " << Linfinity_norm_mat(M) << "\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    using T = mpq_class;
    if (argc != 3 && argc != 5) {
      PrintUsage();
      throw TerminalException{1};
    }
    std::string FileInput = argv[1];
    std::string method = argv[2];
    std::string OutFormat = "CPP";
    std::string FileO = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileO = argv[4];
    }
    MyMatrix<T> M = ReadMatrixFile<T>(FileInput);
    std::cerr << "input: " << M.rows() << " vectors in dimension " << M.cols()
              << ", method=" << method << "\n";
    PrintMeasures("before:", M);
    VectFamilyReductionResult<T> res =
        ReduceVectorFamilyGeneral(M, method, std::cerr);
    PrintMeasures("after:", res.Mred);
    std::cerr << "method used: " << res.method << "\n";
    auto print_mat = [&](std::ostream &os_out) -> void {
      if (OutFormat == "GAP") {
        os_out << "return ";
        WriteMatrixGAP(os_out, res.Mred);
        os_out << ";\n";
        return;
      }
      if (OutFormat == "CPP") {
        return WriteMatrix(os_out, res.Mred);
      }
      std::cerr << "No matching format in print_mat. Allowed: GAP, CPP\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(FileO, print_mat);
    std::cerr << "Normal termination of VectFamily_Reduction\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in VectFamily_Reduction\n";
    exit(e.eVal);
  }
  runtime(time1);
}

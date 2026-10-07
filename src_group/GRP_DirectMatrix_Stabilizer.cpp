// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "GRP_GroupFct.h"
#include "Group.h"
#include "NumberTheory.h"
#include "Permutation.h"
#include "WeightMatrix.h"
// clang-format on

/*
  The stabilizer of a square matrix M given directly: the permutations s of
  the indices with M(s(i), s(j)) = M(i, j) for all i, j, the diagonal
  included. M need not be symmetric; a non-symmetric M goes through the
  non-symmetric weight matrix and its graph on 2 n + 1 vertices, which is
  what this program exists to exercise.
 */
template <typename T, typename Tidx, typename Tint, typename Tidx_value>
StabGeneratorsOrder<Tidx, Tint>
DirectMatrix_Stabilizer_Tidx_value(MyMatrix<T> const &M, bool is_symm,
                                   std::ostream &os) {
  using Tgr = GraphListAdj;
  if (is_symm) {
    WeightMatrix<true, T, Tidx_value> WMat =
        WeightedMatrixFromMyMatrix<true, T, Tidx_value>(M, os);
    return GetStabilizerWeightMatrix_Kernel<T, Tgr, Tidx, Tint, Tidx_value,
                                            true>(WMat, os);
  }
  WeightMatrix<false, T, Tidx_value> WMat =
      WeightedMatrixFromMyMatrix<false, T, Tidx_value>(M, os);
  return GetStabilizerWeightMatrix_Kernel<T, Tgr, Tidx, Tint, Tidx_value,
                                          false>(WMat, os);
}

template <typename T, typename Tidx, typename Tint>
StabGeneratorsOrder<Tidx, Tint> DirectMatrix_Stabilizer(MyMatrix<T> const &M,
                                                        std::ostream &os) {
  size_t n = M.rows();
  bool is_symm = IsSymmetricMatrix(M);
  // A non-symmetric matrix can have one distinct value per entry.
  size_t max_poss_val = weightmatrix_get_nb(is_symm, n);
  auto f_dispatch = [&]<typename Tidx_value>() {
    return DirectMatrix_Stabilizer_Tidx_value<T, Tidx, Tint, Tidx_value>(
        M, is_symm, os);
  };
  return call_with_smallest_unsigned(
      max_poss_val, "DirectMatrix_Stabilizer", f_dispatch);
}

template <typename T, typename Tgroup>
void process(std::string const &FileMat, std::string const &OutFormat,
             std::ostream &os_out) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  MyMatrix<T> M = ReadMatrixFile<T>(FileMat);
  if (M.rows() != M.cols()) {
    std::cerr << "GRP_DirectMatrix_Stabilizer: the matrix should be square, "
                 "we have |M|="
              << M.rows() << " / " << M.cols() << "\n";
    throw TerminalException{1};
  }
  size_t n = M.rows();
  StabGeneratorsOrder<Tidx, typename Tgroup::Tint> gens_order =
      DirectMatrix_Stabilizer<T, Tidx, typename Tgroup::Tint>(M, std::cerr);
  Tgroup GRP = GroupFromStabGeneratorsOrder<Tgroup>(gens_order, n);
  std::cerr << "n=" << n << " is_symmetric=" << IsSymmetricMatrix(M)
            << " |GRP|=" << GRP.size() << "\n";
  if (OutFormat == "GAP") {
    os_out << "return " << GRP.GapString() << ";\n";
    return;
  }
  if (OutFormat == "Oscar") {
    WriteGroup(os_out, GRP);
    return;
  }
  std::cerr << "GRP_DirectMatrix_Stabilizer: No matching entry for "
               "OutFormat="
            << OutFormat << ". Allowed values: GAP, Oscar\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_DirectMatrix_Stabilizer [arith] [FileMatrix] "
                   "[OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_DirectMatrix_Stabilizer [arith] [FileMatrix]\n";
      std::cerr << "\n";
      std::cerr << "It computes the permutations s with M(s(i),s(j)) = "
                   "M(i,j) for the square\n";
      std::cerr << "matrix M of FileMatrix, which does not have to be "
                   "symmetric.\n";
      std::cerr << "\n";
      std::cerr << "arith     : rational\n";
      std::cerr << "OutFormat : GAP or Oscar (default GAP)\n";
      std::cerr << "FileOut   : The output file (default stderr)\n";
      return -1;
    }
    using Tint = mpz_class;
    using Telt = permutalib::SingleSidedPerm<uint32_t>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string arith = argv[1];
    std::string FileMat = argv[2];
    std::string OutFormat = "GAP";
    std::string FileOut = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileOut = argv[4];
    }
    auto f = [&](std::ostream &os_out) -> void {
      if (arith == "rational") {
        using T = mpq_class;
        return process<T, Tgroup>(FileMat, OutFormat, os_out);
      }
      std::cerr << "GRP_DirectMatrix_Stabilizer: No matching entry for arith="
                << arith << ". Allowed values: rational\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(FileOut, f);
    std::cerr << "Normal termination of GRP_DirectMatrix_Stabilizer\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_DirectMatrix_Stabilizer\n";
    exit(e.eVal);
  }
  runtime(time1);
}

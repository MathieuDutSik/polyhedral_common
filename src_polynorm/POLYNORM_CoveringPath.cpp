// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "PolyNorm_CoveringPath.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
MyMatrix<T> ReadPolytopeFile(std::string const &FileEXT, int n_expected) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
  int n = EXT.cols() - 1;
  if (n_expected >= 0 && n != n_expected) {
    std::cerr << "The polytopes should be of the same dimension, " << FileEXT
              << " has dimension " << n << " instead of " << n_expected << "\n";
    throw TerminalException{1};
  }
  for (int i = 0; i < EXT.rows(); i++) {
    if (EXT(i, 0) != 1) {
      std::cerr << "The vertices of " << FileEXT << " should be in "
                << "homogeneous coordinates with first entry 1\n";
      throw TerminalException{1};
    }
  }
  if (EXT.rows() == 0 || RankMat(EXT) != n + 1) {
    std::cerr << "The polytope of " << FileEXT
              << " should be full dimensional\n";
    throw TerminalException{1};
  }
  return EXT;
}

template <typename T>
void ComputeCoveringPath(std::string const &FileEXT1,
                         std::string const &FileEXT2,
                         std::string const &OutFormat, std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> EXT1 = ReadPolytopeFile<T>(FileEXT1, -1);
  MyMatrix<T> EXT2 = ReadPolytopeFile<T>(FileEXT2, EXT1.cols() - 1);
  PolyNormCoveringPath<T> res =
      PolyNorm_ComputeCoveringPath<T, Tint, Tgroup>(EXT1, EXT2, std::cerr);
  if (OutFormat == "text") {
    os_out << "The function g(lambda) = 1 / mu((1 - lambda) P1 + lambda P2) "
           << "on [0,1], " << res.ListPiece.size() << " affine pieces, "
           << "all proved\n";
    for (auto &P : res.ListPiece) {
      os_out << "  lambda in [" << P.alpha << ", " << P.beta
             << "]: g = " << P.s0 << " + " << P.s1 << " lambda,"
             << " mu from " << T(1) / PolyNorm_PathPieceValue(P, P.alpha)
             << " to " << T(1) / PolyNorm_PathPieceValue(P, P.beta) << "\n";
    }
    os_out << "Search: " << res.n_evaluation << " exact evaluations, "
           << res.n_certificate << " certificates, " << res.n_node << " nodes, "
           << res.n_lp << " linear programs\n";
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return [";
    bool IsFirst = true;
    for (auto &P : res.ListPiece) {
      if (!IsFirst) {
        os_out << ",";
      }
      IsFirst = false;
      os_out << "rec(alpha:=" << P.alpha << ", beta:=" << P.beta
             << ", s0:=" << P.s0 << ", s1:=" << P.s1 << ")";
    }
    os_out << "];\n";
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
    std::string FileEXT1, FileEXT2;
    std::string OutFormat = "text";
    std::string OutFile = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics_by_library(dispatch, [&]<typename T>() -> void {
      ComputeCoveringPath<T>(FileEXT1, FileEXT2, OutFormat, *os_ptr);
    });
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLYNORM_CoveringPath [arith] [FileP1] [FileP2] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "POLYNORM_CoveringPath [arith] [FileP1] [FileP2]\n";
      std::cerr << "\n";
      std::cerr << "arith values:\n";
      dispatch.print_descriptions(std::cerr);
      std::cerr
          << "FileP1, FileP2 : the vertices of the polytopes P1 and P2,\n";
      std::cerr << "                 one per row in homogeneous coordinates\n";
      std::cerr << "OutFormat      : text (default) or GAP\n";
      std::cerr << "OutFile        : the output file, or stderr / stdout\n";
      std::cerr << "\n";
      std::cerr << "Computes exactly the function lambda -> 1 / mu((1 -\n";
      std::cerr << "lambda) P1 + lambda P2) on [0,1], mu the covering radius\n";
      std::cerr << "for Z^n: its affine pieces, each one proved\n";
      return -1;
    }
    std::string arith = argv[1];
    FileEXT1 = argv[2];
    FileEXT2 = argv[3];
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    auto prt = [&](std::ostream &os_out) -> void {
      os_ptr = &os_out;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, prt);
    std::cerr << "Normal termination of POLYNORM_CoveringPath\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLYNORM_CoveringPath\n";
    exit(e.eVal);
  }
  runtime(time);
}

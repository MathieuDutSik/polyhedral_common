// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "LatticeMethods.h"
#include "SignatureSymmetric.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

/*
  The conjugacy in O(L) of two isometries g and h of the lattice L of Gram
  matrix A: a U with U A U^T = A and U g U^{-1} = h. The matrices act on
  the row vectors, v -> v g, with g A g^T = A.

  U (A g^T) U^T = U A U^T (U g U^{-1})^T = A (U g U^{-1})^T for U in O(L),
  so U conjugates g to h if and only if it maps the tuple of forms
  (A, A g^T) to (A, A h^T): the conjugacy is the isometry test of the two
  tuples, by any method that applies to several forms. For the conjugacy
  in GL_n(Z), without the lattice, see PERF_FiniteElementConjugacy.
 */

template <typename T>
void ComputeConjugation(std::string const &FileGram, std::string const &FileG,
                        std::string const &FileH, std::string const &OutFormat,
                        std::string const &method, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> A = ReadMatrixFile<T>(FileGram);
  MyMatrix<T> g = ReadMatrixFile<T>(FileG);
  MyMatrix<T> h = ReadMatrixFile<T>(FileH);
  int n = A.rows();
  if (!IsSymmetricMatrix(A) || !IsPositiveDefinite(A, std::cerr)) {
    std::cerr << "LATT_TestConjugation: the Gram matrix in " << FileGram
              << " is not symmetric positive definite\n";
    throw TerminalException{1};
  }
  for (auto const &[name, M] :
       std::vector<std::pair<std::string, MyMatrix<T>>>{{"g", g}, {"h", h}}) {
    if (M.rows() != n || M.cols() != n || !IsIntegralMatrix(M) ||
        M * A * M.transpose() != A) {
      std::cerr << "LATT_TestConjugation: " << name
                << " is not an integral isometry of A (" << name << " A "
                << name << "^T = A expected)\n";
      throw TerminalException{1};
    }
  }
  LatticeMethodChoice choice = ParseLatticeMethodChoice(method);
  if (!LatticeMethodAppliesToMultipleForms(choice.method)) {
    std::cerr << "LATT_TestConjugation: the method " << method
              << " reads the first form only and does not apply\n";
    throw TerminalException{1};
  }
  std::vector<MyMatrix<T>> ListMat1{A, MyMatrix<T>(A * g.transpose())};
  std::vector<MyMatrix<T>> ListMat2{A, MyMatrix<T>(A * h.transpose())};
  std::vector<MyMatrix<Tint>> ListGenAut2;
  std::optional<MyMatrix<Tint>> opt = LatticeIsometry<T, Tint, Tgroup>(
      ListMat1, ListMat2, ListGenAut2, choice, std::cerr);
  if (opt) {
    MyMatrix<T> U = UniversalMatrixConversion<T, Tint>(*opt);
    if (U * A * U.transpose() != A || U * g != h * U) {
      std::cerr << "LATT_TestConjugation: the conjugating matrix is wrong\n";
      throw TerminalException{1};
    }
  }
  if (OutFormat == "GAP") {
    if (opt) {
      os << "return ";
      WriteMatrixGAP(os, *opt);
      os << ";\n";
    } else {
      os << "return false;\n";
    }
    return;
  }
  if (OutFormat == "Oscar") {
    if (opt) {
      WriteMatrix(os, *opt);
    } else {
      os << "0 0\n";
    }
    return;
  }
  std::cerr << "LATT_TestConjugation: unknown OutFormat=" << OutFormat
            << " (GAP or Oscar)\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 5 && argc != 7 && argc != 8) {
      std::cerr << "LATT_TestConjugation [arith] [Gram] [g] [h]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_TestConjugation [arith] [Gram] [g] [h] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_TestConjugation [arith] [Gram] [g] [h] [OutFormat] "
                   "[OutFile] [method]\n";
      std::cerr << "\n";
      std::cerr << "Whether the isometries g and h of the lattice of Gram\n";
      std::cerr << "matrix A are conjugate in O(L): a U with U A U^T = A and\n";
      std::cerr << "U g U^{-1} = h. The matrices act on the row vectors,\n";
      std::cerr << "g A g^T = A. Gram, g, h are matrix files.\n";
      std::cerr << "\n";
      std::cerr << "arith values: gmp (mpq_class / mpz_class)";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << ", flint";
#endif
      std::cerr << "\n";
      std::cerr << "OutFormat values:\n";
      std::cerr << "  GAP   : U or false, GAP readable (default)\n";
      std::cerr << "  Oscar : U in the Oscar format, 0 0 for none\n";
      std::cerr << "OutFile : stderr (default), stdout or a file\n";
      std::cerr << "method values (default heuristic), those that apply to "
                   "several forms:\n";
      for (auto &entry : ListLatticeStabCanIsomMethod()) {
        if (entry.several_forms) {
          std::cerr << "  " << entry.name << " : " << entry.description
                    << "\n";
        }
      }
      return -1;
    }
    std::string arith = argv[1];
    std::string FileGram = argv[2];
    std::string FileG = argv[3];
    std::string FileH = argv[4];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    std::string method = "heuristic";
    if (argc >= 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    if (argc == 8) {
      method = argv[7];
    }
    auto prt = [&](std::ostream &os) -> void {
      ArithmeticDispatch dispatch;
      auto f_arith = [&]<typename T>() -> void {
        ComputeConjugation<T>(FileGram, FileG, FileH, OutFormat, method, os);
      };
      add_rational_arithmetics_by_library(dispatch, f_arith);
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, prt);
    std::cerr << "Normal termination of LATT_TestConjugation\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_TestConjugation\n";
    exit(e.eVal);
  }
  runtime(time);
}

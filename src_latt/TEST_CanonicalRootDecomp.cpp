// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "LatticeRootDecomposition.h"
#include "Positivity.h"
#include <set>
// clang-format on

// A fixed family of GL_n(Z) elements to conjugate by: several shears and
// a sign flip, deterministic so the test is reproducible.
template <typename T>
std::vector<MyMatrix<T>> test_unimodulars(int n) {
  std::vector<MyMatrix<T>> res;
  for (int variant = 0; variant < 3; variant++) {
    MyMatrix<T> U = IdentityMat<T>(n);
    for (int i = 0; i + 1 < n; i++) {
      U(i, i + 1) = T(1 + ((i + variant) % 3));
    }
    MyMatrix<T> L = IdentityMat<T>(n);
    for (int i = 1; i < n; i++) {
      L(i, 0) = T(1 - 2 * ((i + variant) % 2));
    }
    MyMatrix<T> D = IdentityMat<T>(n);
    D(variant % n, variant % n) = T(-1);
    res.push_back(U * L * D);
  }
  return res;
}

template <typename T>
void process(std::string const &FileListGram, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::vector<MyMatrix<T>> ListGram = ReadListMatrixFile<T>(FileListGram);
  int n_error = 0;
  std::set<MyMatrix<T>> distinct;
  for (size_t i = 0; i < ListGram.size(); i++) {
    MyMatrix<T> const &eMat = ListGram[i];
    int n = eMat.rows();
    MyMatrix<Tint> B =
        ComputeCanonicalFormRootDecomposed<T, Tint, Tgroup>(eMat, os);
    MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
    MyMatrix<T> can = B_T * eMat * B_T.transpose();
    distinct.insert(can);
    // (a) isometry invariance: conjugating by a unimodular U must not
    // change the canonical form.
    bool inv_ok = true;
    for (auto &U : test_unimodulars<T>(n)) {
      MyMatrix<T> eMatU = U * eMat * U.transpose();
      MyMatrix<Tint> BU =
          ComputeCanonicalFormRootDecomposed<T, Tint, Tgroup>(eMatU, os);
      MyMatrix<T> BU_T = UniversalMatrixConversion<T, Tint>(BU);
      MyMatrix<T> canU = BU_T * eMatU * BU_T.transpose();
      if (canU != can) {
        inv_ok = false;
        break;
      }
    }
    os << "class " << (i + 1) << ": invariance=" << (inv_ok ? "yes" : "NO")
       << "\n";
    if (!inv_ok) {
      std::cerr << "CANROOTDEC TEST: canonical form not invariant on class "
                << (i + 1) << "\n";
      n_error++;
    }
  }
  os << "CANROOTDEC TEST: " << ListGram.size() << " lattices, "
     << distinct.size() << " distinct canonical forms\n";
  // On a genus (pairwise non-isometric classes) the count of distinct
  // canonical forms must equal the number of lattices.
  if (distinct.size() != ListGram.size()) {
    std::cerr << "CANROOTDEC TEST: " << ListGram.size() << " non-isometric "
              << "lattices collapsed to " << distinct.size()
              << " canonical forms\n";
    n_error++;
  }
  if (n_error > 0) {
    std::cerr << "CANROOTDEC TEST: " << n_error << " errors\n";
    throw TerminalException{1};
  }
  os << "CANROOTDEC TEST: all checks passed\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 2) {
      std::cerr << "TEST_CanonicalRootDecomp [FileListGram]\n";
      return -1;
    }
    std::string FileListGram = argv[1];
    using T = mpq_class;
    process<T>(FileListGram, std::cerr);
    std::cerr << "Normal termination of TEST_CanonicalRootDecomp\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_CanonicalRootDecomp\n";
    exit(e.eVal);
  }
  runtime(time);
}

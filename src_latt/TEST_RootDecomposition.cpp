// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "LatticeRootDecomposition.h"
#include "Positivity.h"
// clang-format on

template <typename T>
void process(std::string const &FileListGram, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::vector<MyMatrix<T>> ListGram = ReadListMatrixFile<T>(FileListGram);
  for (size_t i = 0; i < ListGram.size(); i++) {
    MyMatrix<T> const &eMat = ListGram[i];
    RootDecomposition<T, Tint> dec =
        ComputeRootDecomposition<T, Tint, Tgroup>(eMat, os);
    // Validation: rank(R) + rank(R^perp) = n, and the components carry
    // the right Gram. Round-trip check: R (+) R^perp is a sublattice of
    // index glue_index in L, and its Gram has determinant det(L) *
    // glue_index^2.
    int n = eMat.rows();
    bool ok = (dec.root_rank + dec.PerpBasis.rows() == n);
    std::pair<MyMatrix<T>, MyMatrix<T>> can =
        CanonicalComponentGrams<T, Tint, Tgroup>(dec, os);
    os << "class " << (i + 1) << ": n_roots=" << dec.n_roots
       << " root_rank=" << dec.root_rank
       << " perp_rank=" << dec.PerpBasis.rows()
       << " glue_index=" << dec.glue_index
       << " decomp_ok=" << (ok ? "yes" : "NO") << "\n";
    if (!ok) {
      std::cerr << "ROOTDEC TEST: decomposition rank check failed on class "
                << (i + 1) << "\n";
      throw TerminalException{1};
    }
    // det check: det(RootGram) * det(PerpGram) = det(L) * glue_index^2
    if (dec.root_rank > 0 && dec.PerpBasis.rows() > 0) {
      T dR = DeterminantMat(dec.RootGram);
      T dP = DeterminantMat(dec.PerpGram);
      T dL = DeterminantMat(eMat);
      T gi = UniversalScalarConversion<T, Tint>(dec.glue_index);
      if (dR * dP != dL * gi * gi) {
        std::cerr << "ROOTDEC TEST: det(R)*det(Rperp)=" << (dR * dP)
                  << " != det(L)*glue^2=" << (dL * gi * gi) << " on class "
                  << (i + 1) << "\n";
        throw TerminalException{1};
      }
    }
  }
  os << "ROOTDEC TEST: all " << ListGram.size() << " decompositions valid\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 2) {
      std::cerr << "TEST_RootDecomposition [FileListGram]\n";
      return -1;
    }
    std::string FileListGram = argv[1];
    using T = mpq_class;
    process<T>(FileListGram, std::cerr);
    std::cerr << "Normal termination of TEST_RootDecomposition\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_RootDecomposition\n";
    exit(e.eVal);
  }
  runtime(time);
}

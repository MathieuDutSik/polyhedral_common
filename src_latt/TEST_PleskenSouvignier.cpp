// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "LatticePleskenSouvignier.h"
#include "LatticeStabEquiCan.h"
#include "Positivity.h"
// clang-format on

/*
  Validation of the Plesken-Souvignier engine against the graph-based
  machinery of LatticeStabEquiCan.h, on a list of Gram matrices:

  (1) For every matrix M the automorphism group is computed by both
      methods and the orders are compared, each order being read off the
      permutation action of the generators on the same invariant vector
      family. The product of the orbit lengths of the stabilizer chain
      is checked against the permutation-group order as well.
  (2) For every matrix, the isometry test must find an equivalence from
      M to U * M * U^T for a fixed unimodular U, and the returned P is
      checked to realize it.
  (3) For every consecutive pair of DISTINCT matrices of the list the
      isometry test must fail; the input is expected to be a list of
      pairwise non-isometric lattices, such as the classes of a genus.

  Exits non-zero on the first mismatch.
 */

template <typename T, typename Tint, typename Tgroup>
typename Tgroup::Tint OrderFromGenerators(
    MyMatrix<T> const &GramMat, std::vector<MyMatrix<Tint>> const &ListGen,
    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  // The permutation action on the full antipodal family of the vectors
  // of norm at most the diagonal bound of the LLL-reduced matrix, an
  // invariant of the lattice: faithful, a form-preserving map fixing a
  // full rank family being the identity.
  LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(GramMat, os);
  T bound = PleskenSouvignierBound(rec.GramMatRed);
  MyMatrix<Tint> SHVhalf_red =
      PleskenSouvignierVectorFamily<T, Tint>(rec.GramMatRed, bound, os);
  int n_pair = SHVhalf_red.rows();
  int n = SHVhalf_red.cols();
  MyMatrix<Tint> SHV(2 * n_pair, n);
  for (int i = 0; i < n_pair; i++) {
    MyVector<Tint> V = GetMatrixRow(SHVhalf_red, i);
    // A row y in the reduced basis has original coordinates y * Pmat.
    MyVector<Tint> Vorig = rec.Pmat.transpose() * V;
    for (int k = 0; k < n; k++) {
      SHV(i, k) = Vorig(k);
      SHV(n_pair + i, k) = -Vorig(k);
    }
  }
  int n_row = SHV.rows();
  std::unordered_map<MyVector<Tint>, Tidx> MapRow;
  for (int i = 0; i < n_row; i++) {
    MapRow[GetMatrixRow(SHV, i)] = static_cast<Tidx>(i);
  }
  std::vector<Telt> ListPermGens;
  for (auto &eGen : ListGen) {
    MyMatrix<Tint> Mtr = eGen.transpose();
    std::vector<Tidx> ePerm(n_row);
    for (int i = 0; i < n_row; i++) {
      MyVector<Tint> w = Mtr * GetMatrixRow(SHV, i);
      auto iter = MapRow.find(w);
      if (iter == MapRow.end()) {
        std::cerr << "TEST_PS: a generator does not permute the family\n";
        throw TerminalException{1};
      }
      ePerm[i] = iter->second;
    }
    ListPermGens.push_back(Telt(ePerm));
  }
  Tgroup grp(ListPermGens, n_row);
  return grp.size();
}

template <typename T>
MyMatrix<T> get_test_unimodular(int n) {
  // A fixed deterministic element of GL_n(Z): upper shear, lower shear,
  // and a sign flip, mildly mixing every coordinate.
  MyMatrix<T> U = IdentityMat<T>(n);
  for (int i = 0; i + 1 < n; i++) {
    U(i, i + 1) = T(1 + (i % 3));
  }
  MyMatrix<T> L = IdentityMat<T>(n);
  for (int i = 1; i < n; i++) {
    L(i, 0) = T(1 - 2 * (i % 2));
  }
  MyMatrix<T> D = IdentityMat<T>(n);
  D(0, 0) = T(-1);
  return U * L * D;
}

template <typename T, typename Tint>
void process(std::string const &FileListGram, std::ostream &os) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::vector<MyMatrix<T>> ListGram = ReadListMatrixFile<T>(FileListGram);
  int n_error = 0;
  for (size_t i = 0; i < ListGram.size(); i++) {
    MyMatrix<T> const &eMat = ListGram[i];
    if (!IsSymmetricMatrix(eMat) || !IsPositiveDefinite(eMat, std::cerr)) {
      std::cerr << "TEST_PS: matrix " << i << " is not symmetric positive "
                << "definite\n";
      throw TerminalException{1};
    }
    std::vector<MyMatrix<T>> ListMat{eMat};
    // (1) automorphisms, both methods
    PleskenSouvignierAutomResult<Tint> res_ps =
        PleskenSouvignierLatticeAutomorphism<T, Tint>(ListMat, os);
    mpz_class order_chain =
        PleskenSouvignierGroupOrder<mpz_class>(res_ps.ListOrbitSize);
    mpz_class order_ps = OrderFromGenerators<T, Tint, Tgroup>(
        eMat, res_ps.ListGen, os);
    std::vector<MyMatrix<Tint>> ListGenRef =
        ArithmeticAutomorphismGroup<T, Tint, Tgroup>(eMat, os);
    mpz_class order_ref = OrderFromGenerators<T, Tint, Tgroup>(
        eMat, ListGenRef, os);
    std::cerr << "TEST_PS: matrix " << i << " |Aut| chain=" << order_chain
              << " perm=" << order_ps << " reference=" << order_ref << "\n";
    if (order_chain != order_ps || order_ps != order_ref) {
      std::cerr << "TEST_PS: ERROR, the orders disagree on matrix " << i
                << "\n";
      n_error++;
    }
    // (1b) The same with the searched basis forced to come from the
    // family rows, exercising the general-basis path (shortest-row
    // selection, and the leaf integrality checks when that basis spans
    // a proper sublattice).
    PleskenSouvignierAutomResult<Tint> res_forced =
        PleskenSouvignierLatticeAutomorphism<T, Tint>(ListMat, os, -1, true);
    mpz_class order_forced =
        PleskenSouvignierGroupOrder<mpz_class>(res_forced.ListOrbitSize);
    if (order_forced != order_ref) {
      std::cerr << "TEST_PS: ERROR, the forced-family-basis order "
                << order_forced << " disagrees on matrix " << i << "\n";
      n_error++;
    }
    // (2) isometry against a transformed copy. The automorphisms of the
    // second configuration are the U-conjugates of those of the first,
    // which also exercises the orbit pruning of the search.
    MyMatrix<T> U = get_test_unimodular<T>(eMat.rows());
    MyMatrix<T> eMatB = U * eMat * U.transpose();
    std::vector<MyMatrix<T>> ListMatB{eMatB};
    MyMatrix<Tint> U_int = UniversalMatrixConversion<Tint, T>(U);
    MyMatrix<Tint> Uinv_int = Inverse(U_int);
    std::vector<MyMatrix<Tint>> ListGenAut2;
    for (auto &eGen : res_ps.ListGen) {
      ListGenAut2.push_back(U_int * eGen * Uinv_int);
    }
    std::optional<MyMatrix<Tint>> opt =
        PleskenSouvignierLatticeIsometry<T, Tint>(ListMat, ListMatB,
                                                  ListGenAut2, os);
    if (!opt) {
      std::cerr << "TEST_PS: ERROR, no isometry found from matrix " << i
                << " to its transformed copy\n";
      n_error++;
    } else {
      MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(*opt);
      if (P_T * eMat * P_T.transpose() != eMatB) {
        std::cerr << "TEST_PS: ERROR, the returned isometry is wrong on "
                  << "matrix " << i << "\n";
        n_error++;
      }
    }
    // (3) non-isometry of distinct classes
    if (i + 1 < ListGram.size()) {
      std::vector<MyMatrix<T>> ListMatNext{ListGram[i + 1]};
      std::vector<MyMatrix<Tint>> no_gens;
      std::optional<MyMatrix<Tint>> optNeg =
          PleskenSouvignierLatticeIsometry<T, Tint>(ListMat, ListMatNext,
                                                    no_gens, os);
      if (optNeg) {
        std::cerr << "TEST_PS: ERROR, matrices " << i << " and " << (i + 1)
                  << " found isometric, they should be distinct classes\n";
        n_error++;
      }
    }
  }
  if (n_error > 0) {
    std::cerr << "TEST_PS: " << n_error << " errors\n";
    throw TerminalException{1};
  }
  std::cerr << "TEST_PS: all tests passed on " << ListGram.size()
            << " matrices\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 2) {
      std::cerr << "This program is used as\n";
      std::cerr << "TEST_PleskenSouvignier [FileListGram]\n";
      std::cerr << "\n";
      std::cerr << "FileListGram: a list of pairwise non-isometric symmetric "
                << "positive definite Gram matrices, e.g. the classes of a "
                << "genus\n";
      return -1;
    }
    std::string FileListGram = argv[1];
    using T = mpq_class;
    using Tint = mpz_class;
    process<T, Tint>(FileListGram, std::cerr);
    std::cerr << "Normal termination of TEST_PleskenSouvignier\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_PleskenSouvignier\n";
    exit(e.eVal);
  }
  runtime(time);
}

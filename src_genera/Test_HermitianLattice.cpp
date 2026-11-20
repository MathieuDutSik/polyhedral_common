// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "HermitianLattice.h"
#include <set>
#include <string>
#include <vector>
// clang-format on

/*
  invariants FilePieces d : the invariant (HermitianLattice.h) of each pair
  of a file of pieces (Gram matrix, then isometry acting on the columns, in
  turn), one per hermitian genus over Q(zeta_d) by construction, and
  whether they are distinct.
 */

using T = mpq_class;
using Tint = mpz_class;
using Tidx = uint32_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using Tgroup = permutalib::Group<Telt, Tint>;

int main(int argc, char *argv[]) {
  try {
    if (argc >= 2 && std::string(argv[1]) == "compare") {
      // compare FilePieces d b N : the pieces of the search against those
      // of the file (Hecke): every invariant of the file must be found
      std::vector<MyMatrix<T>> lst = ReadListMatrixFile<T>(argv[2]);
      int64_t d = std::stoll(argv[3]);
      int b = std::stoi(argv[4]);
      Tint N(argv[5]);
      std::set<std::string> ours;
      for (auto const &pr : HermitianPieces<T, Tint, Tgroup>(d, b, N, std::cerr)) {
        ours.insert(HermitianPairInvariant<T, Tint>(pr.first, pr.second, d));
      }
      size_t n_missing = 0, n_file = 0;
      for (size_t i = 0; i + 1 < lst.size(); i += 2) {
        n_file++;
        std::string inv = HermitianPairInvariant<T, Tint>(
            lst[i], MyMatrix<T>(lst[i + 1].transpose()), d);
        if (ours.count(inv) == 0) {
          n_missing++;
          std::cout << "missing: " << inv << "\n";
        }
      }
      std::cout << "file=" << n_file << " ours=" << ours.size()
                << " missing=" << n_missing << "\n";
      return 0;
    }
    if (argc != 4 || std::string(argv[1]) != "invariants") {
      std::cerr << "Test_HermitianLattice invariants FilePieces d\n";
      std::cerr << "Test_HermitianLattice compare FilePieces d b N\n";
      return -1;
    }
    std::vector<MyMatrix<T>> lst = ReadListMatrixFile<T>(argv[2]);
    int64_t d = std::stoll(argv[3]);
    std::set<std::string> seen;
    size_t n_dup = 0;
    for (size_t i = 0; i + 1 < lst.size(); i += 2) {
      MyMatrix<T> const &G = lst[i];
      MyMatrix<T> P = lst[i + 1].transpose();
      std::string inv = HermitianPairInvariant<T, Tint>(G, P, d);
      std::cout << inv << "\n";
      if (!seen.insert(inv).second) {
        n_dup++;
      }
    }
    std::cout << "pieces=" << lst.size() / 2 << " distinct=" << seen.size()
              << " duplicates=" << n_dup << "\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_HermitianLattice\n";
    exit(e.eVal);
  }
}

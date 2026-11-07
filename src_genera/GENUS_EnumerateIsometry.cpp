// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "GenusEnumerationGroup.h"
#include <string>
#include <utility>
#include <vector>
// clang-format on

/*
  The classes of pairs (L, g) in the genus of a pair, g an isometry of
  finite order, by the sigma-stable Kneser neighbours, until the mass of
  the genus of the pair (GENUS_IsometryMass) is reached. See
  GenusEnumerationGroup.h.
 */

template <typename T, typename Tint, typename Tgroup>
void process(std::string const &FileGram, std::string const &FileIsom,
             std::string const &FileMass, std::string const &OutFormat,
             int prime, size_t max_prime, std::ostream &os) {
  std::vector<MyMatrix<T>> ListGram = ReadListMatrixFile<T>(FileGram);
  std::vector<MyMatrix<T>> ListIsom = ReadListMatrixFile<T>(FileIsom);
  if (ListGram.empty() || ListGram.size() != ListIsom.size()) {
    std::cerr << "GENUS_EnumerateIsometry: the files must hold the same "
                 "positive number of Gram matrices and isometries\n";
    throw TerminalException{1};
  }
  std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>> ListSeed;
  for (size_t i = 0; i < ListGram.size(); i++) {
    MyMatrix<T> const &A = ListGram[i];
    MyMatrix<T> const &g = ListIsom[i];
    if (!IsSymmetricMatrix(A) || !IsPositiveDefinite(A, std::cerr)) {
      std::cerr << "GENUS_EnumerateIsometry: a Gram matrix is not symmetric "
                   "positive definite\n";
      throw TerminalException{1};
    }
    if (g.transpose() * A * g != A) {
      std::cerr << "GENUS_EnumerateIsometry: g^T A g is not A for a seed\n";
      throw TerminalException{1};
    }
    ListSeed.push_back({A, g});
  }
  using Tmass = typename GenusEnumerationGroupResult<T, Tgroup>::Tmass;
  Tmass TotalMass = ReadMassFile<Tmass>(FileMass);
  GenusEnumerationGroupResult<T, Tgroup> result =
      GenusEnumerationGroup<T, Tint, Tgroup>(ListSeed, TotalMass, prime,
                                             max_prime, std::cerr);
  WriteGenusEnumerationGroupResult(os, result, OutFormat);
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
    if (argc < 4 || argc > 8) {
      std::cerr << "GENUS_EnumerateIsometry [Gram] [Isometry] [Mass] "
                   "[OutFormat] [OutFile] [prime] [max_prime]\n";
      std::cerr << "\n";
      std::cerr << "Gram     (input) : the Gram matrices of the seeds, in the\n";
      std::cerr << "                   ListMatrix format\n";
      std::cerr << "Isometry (input) : their isometries g (g^T A g = A), in\n";
      std::cerr << "                   the ListMatrix format, one per seed\n";
      std::cerr << "Mass     (input) : the mass of the genus of the pair, as\n";
      std::cerr << "                   \"num den\" (GENUS_IsometryMass)\n";
      std::cerr << "OutFormat        : Summary (default) or GAP\n";
      std::cerr << "OutFile          : stderr (default), stdout or a file\n";
      std::cerr << "prime            : 0 (default): the good primes in\n";
      std::cerr << "                   increasing order; otherwise only that\n";
      std::cerr << "                   prime, which must not divide m det\n";
      std::cerr << "max_prime        : the number of primes brought in at\n";
      std::cerr << "                   most when prime = 0 (default 6)\n";
      return -1;
    }
    std::string FileGram = argv[1];
    std::string FileIsom = argv[2];
    std::string FileMass = argv[3];
    std::string OutFormat = (argc >= 5) ? argv[4] : "Summary";
    std::string OutFile = (argc >= 6) ? argv[5] : "stderr";
    int prime = (argc >= 7) ? std::stoi(argv[6]) : 0;
    size_t max_prime = (argc >= 8) ? std::stoul(argv[7]) : 6;
    auto f = [&](std::ostream &os) -> void {
      process<T, Tint, Tgroup>(FileGram, FileIsom, FileMass, OutFormat, prime,
                               max_prime, os);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of GENUS_EnumerateIsometry\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_EnumerateIsometry\n";
    exit(e.eVal);
  }
  runtime(time);
}

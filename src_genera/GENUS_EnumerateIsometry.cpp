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
  GenusEnumerationGroup.h. With k isometries per seed, the classes of
  (L, G) for the finite abelian group G = <g_1, ..., g_k>, the mass then
  given by GENUS_IsometryMass group; generators that do not commute give
  the classes of (L, G) for a non-abelian G, by the G-stable neighbours of
  the irreducible isotropic submodules.

  Exit status: 0 when the mass closed, 2 when the classes found weigh less
  than the target (incomplete), 3 when they weigh more (the target is
  wrong); see GenusEnumerationGroupExitStatus.
 */

template <typename T, typename Tint, typename Tgroup>
int process(std::string const &FileGram, std::string const &FileIsom,
             std::string const &FileMass, std::string const &OutFormat,
             int prime, size_t max_prime, std::ostream &os) {
  std::vector<MyMatrix<T>> ListGram = ReadListMatrixFile<T>(FileGram);
  std::vector<MyMatrix<T>> ListIsom = ReadListMatrixFile<T>(FileIsom);
  if (ListGram.empty() || ListIsom.empty() ||
      ListIsom.size() % ListGram.size() != 0) {
    std::cerr << "GENUS_EnumerateIsometry: the isometry file must hold k "
                 "isometries per Gram matrix, k >= 1\n";
    throw TerminalException{1};
  }
  size_t k = ListIsom.size() / ListGram.size();
  std::vector<std::pair<MyMatrix<T>, std::vector<MyMatrix<T>>>> ListSeed;
  for (size_t i = 0; i < ListGram.size(); i++) {
    MyMatrix<T> const &A = ListGram[i];
    if (!IsSymmetricMatrix(A) || !IsPositiveDefinite(A, std::cerr)) {
      std::cerr << "GENUS_EnumerateIsometry: a Gram matrix is not symmetric "
                   "positive definite\n";
      throw TerminalException{1};
    }
    std::vector<MyMatrix<T>> ListG;
    for (size_t j = 0; j < k; j++) {
      MyMatrix<T> const &g = ListIsom[i * k + j];
      if (g.transpose() * A * g != A) {
        std::cerr << "GENUS_EnumerateIsometry: g^T A g is not A for a seed\n";
        throw TerminalException{1};
      }
      ListG.push_back(g);
    }
    ListSeed.push_back({A, ListG});
  }
  using Tmass = typename GenusEnumerationGroupResult<T, Tgroup>::Tmass;
  Tmass TotalMass = ReadMassFile<Tmass>(FileMass);
  GenusEnumerationGroupResult<T, Tgroup> result =
      GenusEnumerationAbelianGroup<T, Tint, Tgroup>(ListSeed, TotalMass, prime,
                                                    max_prime, std::cerr);
  WriteGenusEnumerationGroupResult(os, result, OutFormat);
  return GenusEnumerationGroupExitStatus(result);
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
      std::cerr << "                   the ListMatrix format, k per seed\n";
      std::cerr << "                   (consecutive), generating a finite\n";
      std::cerr << "                   group\n";
      std::cerr << "Mass     (input) : the mass of the genus of the pair, as\n";
      std::cerr << "                   \"num den\" (GENUS_IsometryMass)\n";
      std::cerr << "OutFormat        : Summary (default) or GAP\n";
      std::cerr << "OutFile          : stderr (default), stdout or a file\n";
      std::cerr << "prime            : 0 (default): the good primes in\n";
      std::cerr << "                   increasing order; otherwise only that\n";
      std::cerr << "                   prime, which must not divide m det\n";
      std::cerr << "                   (m the exponent of the group)\n";
      std::cerr << "max_prime        : the number of primes brought in at\n";
      std::cerr << "                   most when prime = 0 (default 6)\n";
      std::cerr << "\n";
      std::cerr << "Exit status: 0 complete, 2 incomplete (the classes found\n";
      std::cerr << "weigh less than the mass), 3 the mass is exceeded (the\n";
      std::cerr << "target is wrong), 1 an error.\n";
      return -1;
    }
    std::string FileGram = argv[1];
    std::string FileIsom = argv[2];
    std::string FileMass = argv[3];
    std::string OutFormat = (argc >= 5) ? argv[4] : "Summary";
    std::string OutFile = (argc >= 6) ? argv[5] : "stderr";
    int prime = (argc >= 7) ? std::stoi(argv[6]) : 0;
    size_t max_prime = (argc >= 8) ? std::stoul(argv[7]) : 6;
    int status = 0;
    auto f = [&](std::ostream &os) -> void {
      status = process<T, Tint, Tgroup>(FileGram, FileIsom, FileMass,
                                        OutFormat, prime, max_prime, os);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    if (status == GENUS_GROUP_EXIT_MASS_EXCEEDED) {
      std::cerr << "GENUS_EnumerateIsometry: ERROR: the classes found weigh "
                   "more than the target mass, which is therefore wrong\n";
      runtime(time);
      exit(status);
    }
    if (status == GENUS_GROUP_EXIT_INCOMPLETE) {
      std::cerr << "GENUS_EnumerateIsometry: INCOMPLETE: the mass did not "
                   "close; another prime or another seed is needed\n";
      runtime(time);
      exit(status);
    }
    std::cerr << "Normal termination of GENUS_EnumerateIsometry\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_EnumerateIsometry\n";
    exit(e.eVal);
  }
  runtime(time);
}

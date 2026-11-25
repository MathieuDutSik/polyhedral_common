// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "HermitianLattice.h"
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
// clang-format on

/*
  The pieces (L_2, zeta) of GENUS_GlueIsometry factors: one pair per
  hermitian genus over Z[zeta_d] of rank b, positive definite, whose trace
  lattice is even with det dividing N (HermitianLattice.h), written as a
  list of Gram matrices and isometries acting on the columns (g^T G g = G),
  in turn, as HermitianPieces.jl does with Hecke.jl.

    GENUS_HermitianPieces d b N OutFile
    GENUS_HermitianPieces --jobs FileJobs   (lines "d b N E OutFile", E
                                              unused, as for
                                              HermitianPieces.jl)

  and, for the comparisons with hermitian_genera of Hecke.jl:

    GENUS_HermitianPieces genera d b N
        a line "d b N number": the hermitian genera of rank b over
        Q(zeta_d) / Q(zeta_d)^+ with scales between O_E and N O_E
        (HermitianGenus.h)
    GENUS_HermitianPieces compare FilePieces d b N
        every invariant (HermitianPairInvariant) of the pieces of the file
        (Gram matrices and isometries acting on the columns, in turn) is
        one of the pieces of the search; the missing ones, then a line
        "file= ours= missing="
    GENUS_HermitianPieces invariants FilePieces d
        the invariant of each pair of the file, and whether they are
        distinct
    GENUS_HermitianPieces primes dmax lmax
        lines "d l kind primes e0": for 3 <= d <= dmax, d not 2 mod 4, and
        the primes l < lmax, the type of the primes of Q(zeta_d)^+ over l
        in Q(zeta_d) (0 split, 1 inert, 2 ramified), their number and
        their ramification index (HermPrimesAbove)
 */

using T = mpq_class;
using Tint = mpz_class;
using Tidx = uint32_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using Tgroup = permutalib::Group<Telt, Tint>;

void WritePieces(int64_t d, int b, Tint const &N, std::string const &out) {
  std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>> pieces =
      HermitianPieces<T, Tint, Tgroup>(d, b, N, std::cerr);
  std::ofstream os(out);
  os << 2 * pieces.size() << "\n";
  for (auto const &[G, P] : pieces) {
    WriteMatrix(os, G);
    WriteMatrix(os, MyMatrix<T>(P.transpose()));
  }
  std::cout << out << ": " << pieces.size() << " pieces" << std::endl;
}

void PrintPrimes(int64_t dmax, int64_t lmax) {
  for (int64_t d = 3; d <= dmax; d++) {
    if (d % 4 == 2) {
      continue;
    }
    for (int64_t l = 2; l < lmax; l++) {
      if (!IsPrime(l)) {
        continue;
      }
      HermPrimeData pd = HermPrimesAbove(d, l);
      std::cout << d << " " << l << " " << pd.kind << " " << pd.n_prime << " "
                << pd.e0 << "\n";
    }
  }
}

void PrintGenera(int64_t d, int b, int64_t N) {
  std::vector<HermPrimeSymbols> ps = HermitianPrimeSymbols(d, b, N);
  std::cout << d << " " << b << " " << N << " "
            << HermitianGlobalGenera(ps).size() << "\n";
}

void ComparePieces(std::string const &file, int64_t d, int b,
                   Tint const &N) {
  std::vector<MyMatrix<T>> lst = ReadListMatrixFile<T>(file);
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
}

void PrintInvariants(std::string const &file, int64_t d) {
  std::vector<MyMatrix<T>> lst = ReadListMatrixFile<T>(file);
  std::set<std::string> seen;
  size_t n_dup = 0;
  for (size_t i = 0; i + 1 < lst.size(); i += 2) {
    std::string inv = HermitianPairInvariant<T, Tint>(
        lst[i], MyMatrix<T>(lst[i + 1].transpose()), d);
    std::cout << inv << "\n";
    if (!seen.insert(inv).second) {
      n_dup++;
    }
  }
  std::cout << "pieces=" << lst.size() / 2 << " distinct=" << seen.size()
            << " duplicates=" << n_dup << "\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    std::string mode = (argc > 1) ? argv[1] : "";
    if (mode == "genera" && argc == 5) {
      PrintGenera(std::stoll(argv[2]), std::stoi(argv[3]), std::stoll(argv[4]));
      return 0;
    }
    if (mode == "compare" && argc == 6) {
      ComparePieces(argv[2], std::stoll(argv[3]), std::stoi(argv[4]),
                    Tint(argv[5]));
      return 0;
    }
    if (mode == "invariants" && argc == 4) {
      PrintInvariants(argv[2], std::stoll(argv[3]));
      return 0;
    }
    if (mode == "primes" && argc == 4) {
      PrintPrimes(std::stoll(argv[2]), std::stoll(argv[3]));
      return 0;
    }
    if (argc == 3 && std::string(argv[1]) == "--jobs") {
      std::ifstream is(argv[2]);
      std::string line;
      while (std::getline(is, line)) {
        std::istringstream ls(line);
        int64_t d;
        int b;
        std::string N_s, E_s, out;
        if (line.empty() || line[0] == '#' ||
            !(ls >> d >> b >> N_s >> E_s >> out)) {
          continue;
        }
        WritePieces(d, b, Tint(N_s), out);
      }
    } else if (argc == 5) {
      WritePieces(std::stoll(argv[1]), std::stoi(argv[2]), Tint(argv[3]),
                  argv[4]);
    } else {
      std::cerr << "GENUS_HermitianPieces d b N OutFile\n";
      std::cerr << "GENUS_HermitianPieces --jobs FileJobs\n";
      std::cerr << "GENUS_HermitianPieces genera d b N\n";
      std::cerr << "GENUS_HermitianPieces compare FilePieces d b N\n";
      std::cerr << "GENUS_HermitianPieces invariants FilePieces d\n";
      std::cerr << "GENUS_HermitianPieces primes dmax lmax\n";
      return -1;
    }
    std::cerr << "Normal termination of GENUS_HermitianPieces\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_HermitianPieces\n";
    exit(e.eVal);
  }
  runtime(time);
}

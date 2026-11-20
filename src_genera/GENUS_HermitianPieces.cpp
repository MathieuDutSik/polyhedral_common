// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "HermitianLattice.h"
#include <fstream>
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

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
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
      return -1;
    }
    std::cerr << "Normal termination of GENUS_HermitianPieces\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_HermitianPieces\n";
    exit(e.eVal);
  }
  runtime(time);
}

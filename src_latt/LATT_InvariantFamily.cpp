// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "LatticeInvariantFamily.h"
#include <chrono>
// clang-format on

/*
  LATT_InvariantFamily: the recognition of lattices in a known family of
  pairwise non-isometric lattices (LatticeInvariantFamily.h).

  index: the decision tree of the lattices of FileFamily is built, then the
    index in the family of every lattice of FileQuery is printed, one per
    line; a lattice in no class of the family is an error. With a query
    file of the classes under changes of basis, the expected output is
    their positions.
  cost: for every lattice of FileQuery, the line
      n N_bv t_bv t_can
    with n the rank, N_bv the vectors (up to sign) of the BV invariant, t_bv
    its time and t_can the time of the canonical form (heuristic), in
    seconds, after the reduction of the family: the data of the cost model
    InvariantFamilyCostBV / InvariantFamilyCostCanonical.

  The files are in the ListMatrix format.
 */

using Tidx = uint32_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using Tgroup = permutalib::Group<Telt, mpz_class>;

template <typename Tint>
void process_index(std::string const &FileFamily, std::string const &FileQuery,
                   std::ostream &os_out) {
  std::vector<MyMatrix<Tint>> ListGram = ReadListMatrixFile<Tint>(FileFamily);
  std::vector<MyMatrix<Tint>> ListQuery = ReadListMatrixFile<Tint>(FileQuery);
  InvariantFamily<Tint> family =
      get_invariant_family<Tint, Tgroup>(ListGram, std::cerr);
  std::map<std::string, size_t> count_step;
  for (auto const &node : family.nodes) {
    count_step[NameInvariantFamilyStep(node.step)]++;
  }
  std::cerr << "LATT_InvariantFamily: " << ListGram.size() << " lattices, "
            << family.nodes.size() << " nodes:";
  for (auto const &[name, c] : count_step) {
    std::cerr << " " << name << "=" << c;
  }
  std::cerr << "\n";
  for (auto const &G : ListQuery) {
    os_out << family_invariant_index<Tint, Tgroup>(family, G, std::cerr)
           << "\n";
  }
}

template <typename Tint>
void process_cost(std::string const &FileQuery, std::ostream &os_out) {
  using T = typename overlying_field<Tint>::field_type;
  std::vector<MyMatrix<Tint>> ListGram = ReadListMatrixFile<Tint>(FileQuery);
  auto now = []() { return std::chrono::steady_clock::now(); };
  auto sec = [](auto a, auto b) {
    return std::chrono::duration<double>(b - a).count();
  };
  std::ostream &os = std::cerr;
  for (auto const &G_in : ListGram) {
    MyMatrix<Tint> G_int = InvariantFamilyReduced(G_in, os);
    MyMatrix<T> G = UniversalMatrixConversion<T, Tint>(G_int);
    LatticeBVOptions opt;
    auto t0 = now();
    auto [M, SHV] = BVShortVectors<T, Tint>(G, opt, os);
    (void)M;
    (void)LatticeBV<T, Tint>(G, opt, os);
    auto t1 = now();
    MyMatrix<mpq_class> G_q = InvariantFamilyGmpMatrix(G_int);
    auto t2 = now();
    (void)LatticeCanonicalBasis<mpq_class, mpz_class, Tgroup>(
        G_q, LatticeStabCanIsomMethod::Heuristic, os);
    auto t3 = now();
    os_out << G.rows() << " " << SHV.size() << " " << sec(t0, t1) << " "
           << sec(t2, t3) << std::endl;
  }
}

template <typename Tint>
void process(std::string const &mode, std::vector<std::string> const &files) {
  if (mode == "index" && files.size() == 2) {
    process_index<Tint>(files[0], files[1], std::cout);
    return;
  }
  if (mode == "cost" && files.size() == 1) {
    process_cost<Tint>(files[0], std::cout);
    return;
  }
  std::cerr << "LATT_InvariantFamily: unknown mode " << mode << " with "
            << files.size() << " files\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  if (argc != 4 && argc != 5) {
    std::cerr << "LATT_InvariantFamily [arith] index [FileFamily] [FileQuery]\n";
    std::cerr << "LATT_InvariantFamily [arith] cost [FileQuery]\n";
    std::cerr << "with arith gmp or flint (with ENABLE_FLINT_SUPPORT)\n";
    return -1;
  }
  std::string arith = argv[1];
  std::string mode = argv[2];
  std::vector<std::string> files;
  for (int i = 3; i < argc; i++) {
    files.push_back(argv[i]);
  }
  try {
    if (arith == "gmp") {
      process<mpz_class>(mode, files);
    }
#ifdef ENABLE_FLINT_SUPPORT
    else if (arith == "flint") {
      process<fmpz_class>(mode, files);
    }
#endif
    else {
      std::cerr << "LATT_InvariantFamily: unknown arithmetic " << arith << "\n";
      return -1;
    }
    std::cerr << "Normal termination of LATT_InvariantFamily\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_InvariantFamily\n";
    return e.eVal;
  }
  runtime(time);
}

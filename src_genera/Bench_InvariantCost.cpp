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
  The costs of the two expensive invariants of LatticeInvariantFamily.h, the
  BV invariant and the canonical form, on the lattices of a ListMatrix file
  (every stride-th one, at most max_count), to build a model of their costs.
  One line per lattice:
    n N_bv t_bv N_can t_family t_can
  with N_bv the vectors (up to sign) of the BV invariant, N_can those of
  the canonical vector family, t_family the time of that family and t_can
  the time of the whole canonical form, in seconds. The lattices are first
  LLL-reduced, as in the invariant family.
 */

using Tidx = uint32_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using Tgroup = permutalib::Group<Telt, mpz_class>;

template <typename Tint>
void run_bench(std::string const &FileI, size_t stride, size_t max_count) {
  using T = typename overlying_field<Tint>::field_type;
  std::vector<MyMatrix<Tint>> ListGram = ReadListMatrixFile<Tint>(FileI);
  auto now = []() { return std::chrono::steady_clock::now(); };
  auto sec = [](auto a, auto b) {
    return std::chrono::duration<double>(b - a).count();
  };
  std::ostream &os = std::cerr;
  size_t count = 0;
  for (size_t i = 0; i < ListGram.size() && count < max_count; i += stride) {
    MyMatrix<Tint> G_int = InvariantFamilyReduced(ListGram[i], os);
    MyMatrix<T> G = UniversalMatrixConversion<T, Tint>(G_int);
    LatticeBVOptions opt;
    auto t0 = now();
    auto [M, SHV] = BVShortVectors<T, Tint>(G, opt, os);
    (void)M;
    (void)LatticeBV<T, Tint>(G, opt, os);
    auto t1 = now();
    MyMatrix<mpq_class> G_q = InvariantFamilyGmpMatrix(G_int);
    auto t2 = now();
    CanonicVectorFamily<mpz_class> fam =
        GetCanonicVectorFamily<mpq_class, mpz_class>(G_q, os);
    auto t3 = now();
    (void)LatticeCanonicalBasis<mpq_class, mpz_class, Tgroup>(
        G_q, LatticeStabCanIsomMethod::Heuristic, os);
    auto t4 = now();
    std::cout << G.rows() << " " << SHV.size() << " " << sec(t0, t1) << " "
              << fam.SHVhalf.rows() << " " << sec(t2, t3) << " "
              << sec(t3, t4) << std::endl;
    count++;
  }
}

int main(int argc, char *argv[]) {
  if (argc != 5) {
    std::cerr << "Bench_InvariantCost [gmp|flint] [FileListMatrix] [stride] "
                 "[max_count]\n";
    return -1;
  }
  std::string arith = argv[1];
  size_t stride = std::stoul(argv[3]);
  size_t max_count = std::stoul(argv[4]);
  try {
    if (arith == "gmp") {
      run_bench<mpz_class>(argv[2], stride, max_count);
      return 0;
    }
#ifdef ENABLE_FLINT_SUPPORT
    if (arith == "flint") {
      run_bench<fmpz_class>(argv[2], stride, max_count);
      return 0;
    }
#endif
    std::cerr << "Bench_InvariantCost: the arithmetic is gmp or flint\n";
    return -1;
  } catch (TerminalException const &e) {
    std::cerr << "Error in Bench_InvariantCost\n";
    return e.eVal;
  }
}

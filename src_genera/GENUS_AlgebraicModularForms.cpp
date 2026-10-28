// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "AlgebraicModularForms.h"
#include "Arithmetic_dispatch_types.h"
#include <fstream>
#include <string>
#include <vector>
// clang-format on

namespace {

template <typename T, typename Tint>
void write_amf_result(std::ostream &os, AmfResult<T, Tint> const &res, std::string const &level) {
  os << "return rec(level:=\"" << level << "\", n_class:=" << res.l_class.size() << ", n_state:=" << res.n_state
     << ",\nListAutOrder:=[";
  for (size_t j = 0; j < res.l_class.size(); j++) {
    os << (j > 0 ? ", " : "") << res.l_class[j].aut_order;
  }
  os << "],\nListNumberOrbit:=[";
  for (size_t j = 0; j < res.l_class.size(); j++) {
    os << (j > 0 ? ", " : "") << res.l_class[j].l_orbit_repr.size();
  }
  os << "],\nmass:=" << res.mass << ",\nmass_formula:=";
  if (res.mass_formula) {
    os << *res.mass_formula << ", complete:=" << (*res.mass_formula == res.mass ? "true" : "false");
  } else {
    os << "fail, complete:=fail";
  }
  os << ",\nListHecke:=[";
  for (size_t ip = 0; ip < res.l_hecke.size(); ip++) {
    AmfHeckeMatrix const &hm = res.l_hecke[ip];
    os << (ip > 0 ? ",\n" : "") << "rec(prime:=\"" << hm.P.name << "\", n_neighbor:=" << hm.n_neighbor
       << ", matrix:=[";
    for (size_t r = 0; r < res.n_state; r++) {
      os << (r > 0 ? ",\n" : "") << "[";
      for (size_t c = 0; c < res.n_state; c++) {
        auto iter = hm.entries.find({r, c});
        os << (c > 0 ? "," : "") << (iter == hm.entries.end() ? 0 : iter->second);
      }
      os << "]";
    }
    os << "])";
  }
  os << "]);\n";
}

template <typename T, typename Tint>
void process_amf(std::string const &file_desc, std::string const &str_level, std::string const &str_primes,
                 std::string const &file_out) {
  DescriptionTokenReader<Tint> tr{read_description_tokens(file_desc)};
  RealifiedRing<T, Tint> ring = read_realified_ring<T, Tint>(tr);
  if (ring.d == 1 && str_level != "1") {
    std::cerr << "GENUS_AlgebraicModularForms: over Z only the trivial level for now\n";
    throw TerminalException{1};
  }
  std::optional<OrderPrime> N;
  if (str_level != "1") {
    N = parse_order_primes(ring, str_level)[0];
  }
  std::vector<OrderPrime> l_prime;
  for (auto &P : parse_order_primes(ring, str_primes)) {
    if (N && N->ell == P.ell) {
      std::cerr << "GENUS_AlgebraicModularForms: the prime " << P.name << " lies over the same rational prime as "
                << "the level\n";
      throw TerminalException{1};
    }
    l_prime.push_back(P);
  }
  AmfResult<T, Tint> res = compute_algebraic_modular_forms(ring, N, l_prime, std::cerr);
  std::cerr << "GENUS_AlgebraicModularForms: " << res.l_class.size() << " classes, mass " << res.mass << ", "
            << res.n_state << " states\n";
  if (res.mass_formula) {
    std::cerr << "GENUS_AlgebraicModularForms: mass formula " << *res.mass_formula << ", the genus is "
              << (*res.mass_formula == res.mass ? "complete" : "NOT complete") << "\n";
  } else {
    std::cerr << "GENUS_AlgebraicModularForms: no mass formula for this genus, completeness not certified\n";
  }
  std::ofstream os(file_out);
  write_amf_result(os, res, str_level);
}

}  // namespace

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    std::string file_desc, str_level, str_primes, file_out;
    ArithmeticDispatch dispatch;
    add_rational_integer_arithmetics_by_library(dispatch, [&]<typename T, typename Tint>() -> void {
      process_amf<T, Tint>(file_desc, str_level, str_primes, file_out);
    });
    if (argc != 6) {
      std::cerr << "GENUS_AlgebraicModularForms [arith] [description] [level] [primes] [output]\n";
      std::cerr << "\n";
      std::cerr << "Algebraic modular forms of trivial weight on the unitary group of a\n";
      std::cerr << "definite Hermitian lattice over O_K, or the orthogonal group of a\n";
      std::cerr << "positive definite even lattice over Z, and their Hecke operators.\n";
      std::cerr << "  arith values:\n";
      dispatch.print_descriptions(std::cerr);
      std::cerr << "  description: the ring, the dimension and the form of the seed lattice\n";
      std::cerr << "               R^m (see RealifiedRing.h)\n";
      std::cerr << "  level: 1, or a prime of O_K: an inert integer prime, a generator a+bw\n";
      std::cerr << "         of a split prime (w = omega), or ell@r for the split prime\n";
      std::cerr << "         (ell, omega - r), which need not be principal\n";
      std::cerr << "  primes: comma separated list of the Hecke primes, given in the same way\n";
      std::cerr << "  output: the GAP file of the classes and the Hecke matrices\n";
      return -1;
    }
    std::string arith = argv[1];
    file_desc = argv[2];
    str_level = argv[3];
    str_primes = argv[4];
    file_out = argv[5];
    dispatch.run(arith);
    std::cerr << "Normal termination of GENUS_AlgebraicModularForms\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_AlgebraicModularForms\n";
    exit(e.eVal);
  }
  runtime(time);
}

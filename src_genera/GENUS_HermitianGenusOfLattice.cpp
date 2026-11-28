// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "HermitianGenusSymbol.h"
#include <string>
// clang-format on

/*
  GENUS_HermitianGenusOfLattice: the genus symbol (HermitianGenusSymbol.h)
  of a Hermitian lattice over a CM field, given by a description file
  (read_cm_realified_ring of CMOrder.h):

    ring NumberField deg c_0 ... c_deg conj t_0 ... t_{deg-1}
    dimension m
    form [h_0] ... [h_{deg-1}]
    generators r [g_0] ... [g_{deg-1}]   (optional)

  Without generators the lattice is O_E^m. With them it is the O_E-span of
  the r rows of g = sum_k theta^k g_k (r x m matrices, rational entries),
  which need not be free: a lattice of rank m that is not free is given by
  m + 1 generators. The symbol is the input of GENUS_HermitianEnumerate.
 */

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    using T = mpq_class;
    using Tint = mpz_class;
    if (argc != 2 && argc != 3) {
      std::cerr << "GENUS_HermitianGenusOfLattice [description] [output]\n";
      return -1;
    }
    DescriptionTokenReader<Tint> tr{read_description_tokens(argv[1])};
    auto [ring, opt_ord] = read_cm_realified_ring<T, Tint>(tr, std::cerr);
    if (!opt_ord) {
      std::cerr << "GENUS_HermitianGenusOfLattice: the ring must be given as "
                   "a NumberField\n";
      throw TerminalException{1};
    }
    CMOrder<T, Tint> const &ord = *opt_ord;
    MyMatrix<T> A = read_herm_lattice_basis<T, Tint>(tr, ord, ring);
    HermLatticeData<T, Tint> L = HermLatticeOfBasis(ord, ring, A);
    std::string sym = HermGlobalSymbolString(HermGenusSymbolOfLattice(L));
    std::string out = (argc == 3) ? argv[2] : "stdout";
    auto f = [&](std::ostream &os) -> void { os << sym << "\n"; };
    FILE_PrintStderrStdoutFile(out, f);
    std::cerr << "Normal termination of GENUS_HermitianGenusOfLattice\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_HermitianGenusOfLattice\n";
    exit(e.eVal);
  }
  runtime(time);
}

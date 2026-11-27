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
    int d = ord.deg;
    int m = ring.m;
    int n = ring.n();
    MyMatrix<T> A = IdentityMat<T>(n);
    if (!tr.at_end()) {
      tr.expect("generators");
      int r = tr.next_int();
      // the matrices g_k, read with rational entries
      std::vector<MyMatrix<T>> l_g;
      for (int k = 0; k < d; k++) {
        int n_row = tr.next_int();
        int n_col = tr.next_int();
        if (n_row != r || n_col != m) {
          std::cerr << "GENUS_HermitianGenusOfLattice: the generator matrices "
                       "must be " << r << " x " << m << "\n";
          throw TerminalException{1};
        }
        MyMatrix<T> g(r, m);
        for (int i = 0; i < r; i++) {
          for (int j = 0; j < m; j++) {
            g(i, j) = ParseScalar<T>(tr.next());
          }
        }
        l_g.push_back(g);
      }
      // a generator in the realified coordinates (e_a times the i-th unit
      // vector at a m + i), then its multiples by the basis e of O_E
      MyMatrix<T> Binv = Inverse(ord.basis_theta);
      std::vector<MyMatrix<T>> l_mult;
      for (auto const &M : ring.l_ring_basis) {
        l_mult.push_back(UniversalMatrixConversion<T, Tint>(M));
      }
      MyMatrix<T> Gens(r * d, n);
      for (int i = 0; i < r; i++) {
        MyVector<T> v = ZeroVector<T>(n);
        for (int j = 0; j < m; j++) {
          MyVector<T> c(d);
          for (int k = 0; k < d; k++) {
            c(k) = l_g[k](i, j);
          }
          MyVector<T> x = Binv.transpose() * c;
          for (int a = 0; a < d; a++) {
            v(a * m + j) = x(a);
          }
        }
        for (int b = 0; b < d; b++) {
          MyVector<T> w = l_mult[b].transpose() * v;
          for (int k = 0; k < n; k++) {
            Gens(i * d + b, k) = w(k);
          }
        }
      }
      A = GetZbasis(Gens);
      if (A.rows() != n) {
        std::cerr << "GENUS_HermitianGenusOfLattice: the generators span a "
                     "lattice of rank " << A.rows() / d << ", not " << m
                  << "\n";
        throw TerminalException{1};
      }
    }
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

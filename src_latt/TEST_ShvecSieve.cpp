// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeReductionBench.h"
#include "Shvec_sieve.h"
// clang-format on

/*
  The Gauss sieve of Shvec_sieve.h against the enumeration of Shvec_exact.h.

  Each instance is a lattice of the benchmark families Zn, An, Dn, E8 and
  random, presented through a random unimodular change of basis. On each of
  them are run T_ShortestVector, the enumeration, and T_ShortestVectorSieve,
  the sieve followed by the certifying enumeration. The two must return the
  same minimum and the same set of minimal vectors, and any difference is an
  error. What is also reported is how often the sieve alone, uncertified,
  found the minimum, and the times of the three.

  Usage:
    TEST_ShvecSieve [dim] [n_iter] [seed]
 */

template <typename Tint>
std::vector<MyVector<Tint>> SortedRows(MyMatrix<Tint> const &M) {
  std::vector<MyVector<Tint>> l_row;
  for (int i = 0; i < M.rows(); i++) {
    l_row.push_back(GetMatrixRow(M, i));
  }
  std::sort(l_row.begin(), l_row.end(),
            [](MyVector<Tint> const &a, MyVector<Tint> const &b) {
              return std::lexicographical_compare(
                  a.data(), a.data() + a.size(), b.data(), b.data() + b.size());
            });
  return l_row;
}

template <typename T>
void process(int dim, int n_iter, unsigned long seed, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  std::mt19937_64 rng(seed);
  std::vector<std::string> l_name{"Zn", "An", "Dn", "E8", "random"};
  int n_ops = 40;
  int n_case = 0, n_sieve_min = 0;
  double ms_enum = 0, ms_sieve = 0, ms_certified = 0;
  for (auto &name : l_name) {
    int n_use = (name == "E8") ? 8 : dim;
    if (name == "Dn" && n_use < 4) {
      continue;
    }
    for (int i_iter = 0; i_iter < n_iter; i_iter++) {
      HiddenBasisInstance<T, Tint> inst =
          MakeHiddenBasisInstance<T, Tint>(name, n_use, n_ops, rng);
      MyMatrix<T> const &G = inst.GramBad;
      SieveOptions<T> opts;
      opts.seed = seed + static_cast<unsigned long>(n_case);
      //
      MicrosecondTime time_enum;
      Tshortest<T, Tint> rec_enum = T_ShortestVector<T, Tint>(G, os);
      double t_enum = static_cast<double>(time_enum.const_eval_int64()) / 1000;
      //
      MicrosecondTime time_sieve;
      SieveResult<T, Tint> res = GaussSieve<T, Tint>(G, opts, os);
      double t_sieve =
          static_cast<double>(time_sieve.const_eval_int64()) / 1000;
      //
      MicrosecondTime time_cert;
      Tshortest<T, Tint> rec_sieve =
          T_ShortestVectorSieve<T, Tint>(G, opts, os);
      double t_cert = static_cast<double>(time_cert.const_eval_int64()) / 1000;
      //
      if (rec_enum.min != rec_sieve.min ||
          SortedRows(rec_enum.SHV) != SortedRows(rec_sieve.SHV)) {
        std::cerr << "TEST_ShvecSieve: mismatch on " << name << " n=" << inst.n
                  << " iter=" << i_iter << " enumeration min=" << rec_enum.min
                  << " |SHV|=" << rec_enum.SHV.rows()
                  << " sieve min=" << rec_sieve.min
                  << " |SHV|=" << rec_sieve.SHV.rows() << "\n";
        std::cerr << "G=\n";
        WriteMatrix(std::cerr, G);
        throw TerminalException{1};
      }
      bool sieve_found = (res.norms[0] == rec_enum.min);
      if (sieve_found) {
        n_sieve_min++;
      }
      os << "case " << name << " n=" << inst.n << " iter=" << i_iter
         << " min=" << rec_enum.min << " |SHV|=" << rec_enum.SHV.rows()
         << " sieve_first=" << res.norms[0] << " |L|=" << res.list.size()
         << " samples=" << res.n_samples << " enum=" << t_enum
         << "ms sieve=" << t_sieve << "ms certified=" << t_cert << "ms\n";
      ms_enum += t_enum;
      ms_sieve += t_sieve;
      ms_certified += t_cert;
      n_case++;
    }
  }
  os << "\nTEST_ShvecSieve: " << n_case
     << " cases, certified sieve = enumeration on all of them\n";
  os << "  the sieve alone found the minimum on " << n_sieve_min << "/"
     << n_case << "\n";
  os << "  total times: enumeration=" << ms_enum << "ms sieve=" << ms_sieve
     << "ms certified=" << ms_certified << "ms\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TEST_ShvecSieve [dim] [n_iter] [seed]\n";
      return -1;
    }
    int dim = std::stoi(argv[1]);
    int n_iter = std::stoi(argv[2]);
    unsigned long seed = std::stoul(argv[3]);
    process<mpq_class>(dim, n_iter, seed, std::cout);
    std::cerr << "Normal termination of TEST_ShvecSieve\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_ShvecSieve\n";
    exit(e.eVal);
  }
  runtime(time);
}

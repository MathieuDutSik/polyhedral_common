// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeReductionBench.h"
#include "MinkowskiReduction.h"
// clang-format on

/*
  Validates the Minkowski reduction.

  The decisive check is IsMinkowskiReduced, which redoes the enumeration from
  the output matrix and asks, at every index, whether any admissible vector is
  strictly shorter than the one in place. It shares no state with the descent,
  so it tests the condition rather than restating the algorithm.

  Two further properties are asserted. A Minkowski reduced basis has
  |b_1| = lambda_1(L), the first index being admissible for every primitive
  vector, so the leading diagonal entry must equal the true minimum, which is
  computed here independently by T_ShortestVector. And |b_i| <= |b_i + b_j| for
  every j < i, since b_i + b_j is admissible at index i whenever b_i is: that
  gives 2|<b_i, b_j>| <= |b_j|^2, which is size reduction with respect to the
  basis vectors themselves.

  The dimensions are kept small on purpose. The cost of the method is
  exponential in the dimension, and inherently so.

  Usage:
    TEST_MinkowskiReduction [dim] [n_iter] [seed]
 */

template <typename T, typename Tint>
void process(int dim, int n_iter, unsigned long seed, std::ostream &os) {
  std::mt19937_64 rng(seed);
  std::vector<std::string> l_name{"Zn", "An", "Dn", "random"};
  std::vector<int> l_nops{10, 40};
  int n_case = 0, n_error = 0;
  for (auto &name : l_name) {
    int n_use = dim;
    if (name == "Dn" && n_use < 4) {
      continue;
    }
    for (auto &n_ops : l_nops) {
      for (int i_iter = 0; i_iter < n_iter; i_iter++) {
        HiddenBasisInstance<T, Tint> inst =
            MakeHiddenBasisInstance<T, Tint>(name, n_use, n_ops, rng);
        int n = inst.n;
        LLLreduction<T, Tint> res =
            MinkowskiReducedBasis<T, Tint>(inst.GramBad, os);
        std::string tag = name + "/" + std::to_string(n_ops) + "/" +
                          std::to_string(i_iter);
        //
        // The transformation is unimodular and produces the stated form.
        //
        MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(res.Pmat);
        if (P_T * inst.GramBad * P_T.transpose() != res.GramMatRed) {
          std::cerr << "MINKOWSKI TEST: " << tag
                    << ": the transformation does not produce the reduced "
                       "form\n";
          n_error++;
        }
        T det_P = DeterminantMat(P_T);
        if (det_P != 1 && det_P != -1) {
          std::cerr << "MINKOWSKI TEST: " << tag
                    << ": the transformation is not unimodular, det=" << det_P
                    << "\n";
          n_error++;
        }
        //
        // The condition itself.
        //
        if (!IsMinkowskiReduced<T, Tint>(res.GramMatRed, os)) {
          std::cerr << "MINKOWSKI TEST: " << tag
                    << ": the output is not Minkowski reduced\n";
          n_error++;
        }
        //
        // |b_1|^2 is the true minimum of the lattice.
        //
        Tshortest<T, Tint> shv =
            T_ShortestVector<T, Tint>(res.GramMatRed, os);
        if (res.GramMatRed(0, 0) != shv.min) {
          std::cerr << "MINKOWSKI TEST: " << tag << ": |b_1|^2 is "
                    << res.GramMatRed(0, 0) << " but the minimum is "
                    << shv.min << "\n";
          n_error++;
        }
        //
        // 2|<b_i,b_j>| <= |b_j|^2 for j < i, a consequence of the condition
        // at index i applied to b_i +- b_j.
        //
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < i; j++) {
            T lhs = 2 * T_abs(res.GramMatRed(i, j));
            if (lhs > res.GramMatRed(j, j)) {
              std::cerr << "MINKOWSKI TEST: " << tag << ": 2|G(" << i << ","
                        << j << ")|=" << lhs << " exceeds G(" << j << "," << j
                        << ")=" << res.GramMatRed(j, j) << "\n";
              n_error++;
            }
          }
        }
        os << "MINKOWSKI TEST: " << tag << " min=" << shv.min
           << " diagonal";
        for (int i = 0; i < n; i++) {
          os << " " << res.GramMatRed(i, i);
        }
        os << "\n";
        n_case++;
      }
    }
  }
  if (n_error > 0) {
    std::cerr << "MINKOWSKI TEST: " << n_error << " errors over " << n_case
              << " cases\n";
    throw TerminalException{1};
  }
  os << "MINKOWSKI TEST: all " << n_case << " cases pass\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 1 && argc != 4) {
      std::cerr << "TEST_MinkowskiReduction [dim] [n_iter] [seed]\n";
      std::cerr << "or\n";
      std::cerr << "TEST_MinkowskiReduction\n";
      std::cerr << "for the default dim=5, n_iter=2, seed=1\n";
      return -1;
    }
    int dim = 5;
    int n_iter = 2;
    unsigned long seed = 1;
    if (argc == 4) {
      dim = std::stoi(argv[1]);
      n_iter = std::stoi(argv[2]);
      seed = std::stoul(argv[3]);
    }
    using T = mpq_class;
    using Tint = mpz_class;
    process<T, Tint>(dim, n_iter, seed, std::cerr);
    std::cerr << "Normal termination of TEST_MinkowskiReduction\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_MinkowskiReduction\n";
    exit(e.eVal);
  }
  runtime(time);
}

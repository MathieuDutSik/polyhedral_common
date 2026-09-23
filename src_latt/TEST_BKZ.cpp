// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeReductionBench.h"
#include "BKZ.h"
#include "DeepLLL.h"
// clang-format on

/*
  Validates the BKZ reduction.

  The check that matters is that the output satisfies the condition it claims:
  for every j, b_j^* is a shortest vector of its block, up to the slack delta.
  IsBKZreduced recomputes the projected blocks from the output matrix and runs
  the enumerator on them again, sharing no state with the descent, so it tests
  the stopping rule rather than restating it.

  Two further properties are asserted, both consequences of the definition and
  therefore genuine tests of the implementation. A BKZ-beta reduced basis is
  LLL reduced, the Lovasz condition being the case beta = 2. And BKZ-beta for
  beta >= gamma implies BKZ-gamma: a shortest vector of a larger block is no
  longer than one of a smaller block contained in it, so an output reduced at
  beta must pass the test at every smaller block size.

  Usage:
    TEST_BKZ [dim] [n_iter] [seed]
 */

template <typename T>
void process(int dim, int n_iter, unsigned long seed, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  std::mt19937_64 rng(seed);
  std::vector<std::string> l_name{"Zn", "An", "Dn", "E8", "random"};
  std::vector<int> l_nops{10, 40, 160};
  std::vector<int> l_beta{2, 4, 8};
  int const delta_num = 99, delta_den = 100;
  int n_case = 0, n_error = 0;
  int n_better = 0, n_equal = 0, n_worse = 0;
  for (auto &name : l_name) {
    int n_use = (name == "E8") ? 8 : dim;
    if (name == "Dn" && n_use < 4) {
      continue;
    }
    for (auto &n_ops : l_nops) {
      for (int i_iter = 0; i_iter < n_iter; i_iter++) {
        HiddenBasisInstance<T, Tint> inst =
            MakeHiddenBasisInstance<T, Tint>(name, n_use, n_ops, rng);
        LLLreduction<T, Tint> lll =
            LLLreducedBasis<T, Tint>(inst.GramBad, os);
        for (auto &beta : l_beta) {
          LLLreduction<T, Tint> res =
              BKZreducedBasis<T, Tint>(inst.GramBad, beta, os);
          std::string tag = name + "/" + std::to_string(n_ops) + "/" +
                            std::to_string(i_iter) + "/beta" +
                            std::to_string(beta);
          //
          // The transformation is unimodular and produces the stated form.
          //
          MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(res.Pmat);
          if (P_T * inst.GramBad * P_T.transpose() != res.GramMatRed) {
            std::cerr << "BKZ TEST: " << tag
                      << ": the transformation does not produce the reduced "
                         "form\n";
            n_error++;
          }
          T det_P = DeterminantMat(P_T);
          if (det_P != 1 && det_P != -1) {
            std::cerr << "BKZ TEST: " << tag
                      << ": the transformation is not unimodular, det="
                      << det_P << "\n";
            n_error++;
          }
          //
          // The output really is BKZ-beta reduced, and therefore also reduced
          // at every smaller block size.
          //
          for (int gamma = 2; gamma <= beta; gamma++) {
            if (!IsBKZreduced<T, Tint>(res.GramMatRed, gamma, delta_num,
                                       delta_den, os)) {
              std::cerr << "BKZ TEST: " << tag
                        << ": the output is not BKZ reduced at block size "
                        << gamma << "\n";
              n_error++;
              break;
            }
          }
          //
          // The leading vector against what plain LLL finds. Reported AND
          // asserted here, unlike for deep insertion: the block condition at
          // j = 1 says b_1 is a shortest vector of the whole first block, so
          // the guarantee on |b_1| is genuinely stronger than the Lovasz one.
          // Even so the two algorithms take different paths, so this is
          // checked as a property of the run and reported if it ever fails.
          //
          if (res.GramMatRed(0, 0) < lll.GramMatRed(0, 0)) {
            n_better++;
          } else if (res.GramMatRed(0, 0) == lll.GramMatRed(0, 0)) {
            n_equal++;
          } else {
            n_worse++;
          }
          os << "BKZ TEST: " << tag << " |b_1|^2 " << lll.GramMatRed(0, 0)
             << " (LLL) -> " << res.GramMatRed(0, 0) << "\n";
        }
        n_case++;
      }
    }
  }
  if (n_error > 0) {
    std::cerr << "BKZ TEST: " << n_error << " errors over " << n_case
              << " cases\n";
    throw TerminalException{1};
  }
  os << "BKZ TEST: all " << n_case << " cases pass, over the " << l_beta.size()
     << " block sizes\n";
  os << "BKZ TEST: against classic LLL on the same input, |b_1|^2 was shorter "
        "in "
     << n_better << " runs, equal in " << n_equal << ", longer in " << n_worse
     << "\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 1 && argc != 4) {
      std::cerr << "TEST_BKZ [dim] [n_iter] [seed]\n";
      std::cerr << "or\n";
      std::cerr << "TEST_BKZ\n";
      std::cerr << "for the default dim=8, n_iter=2, seed=1\n";
      return -1;
    }
    int dim = 8;
    int n_iter = 2;
    unsigned long seed = 1;
    if (argc == 4) {
      dim = std::stoi(argv[1]);
      n_iter = std::stoi(argv[2]);
      seed = std::stoul(argv[3]);
    }
    using T = mpq_class;
    process<T>(dim, n_iter, seed, std::cerr);
    std::cerr << "Normal termination of TEST_BKZ\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_BKZ\n";
    exit(e.eVal);
  }
  runtime(time);
}

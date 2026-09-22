// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeReductionBench.h"
#include "SlideReduction.h"
// clang-format on

/*
  Validates the slide reduction.

  The decisive check is IsSlideReduced, which rebuilds every block from the
  output and re-tests both families: each disjoint primal block HKZ-reduced,
  and each shifted dual block having its last Gram-Schmidt norm maximal. It
  shares no state with the descent and re-runs the enumerator itself, so it
  tests the stopping rule rather than restating it.

  The dual half is the part worth testing hardest, being the part that is not
  simply BKZ on a block: it goes through the reversed dual J adj(G) J, finds a
  shortest vector there, and maps the transformation back by J U^{-T} J. An
  error in any of those three would leave the last Gram-Schmidt norm
  unimproved, which the verifier would catch, or would break the lattice, which
  the unimodularity check would catch.

  Usage:
    TEST_SlideReduction [dim] [n_iter] [seed]
 */

template <typename T, typename Tint>
void process(int dim, int n_iter, unsigned long seed, std::ostream &os) {
  std::mt19937_64 rng(seed);
  std::vector<std::string> l_name{"Zn", "An", "Dn", "E8", "random"};
  std::vector<int> l_nops{10, 40, 160};
  std::vector<int> l_kmax{2, 4, 6};
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
        for (auto &k_max : l_kmax) {
          int k = SlideBlockSize(inst.n, k_max);
          LLLreduction<T, Tint> res =
              SlideReducedBasis<T, Tint>(inst.GramBad, k, os);
          std::string tag = name + "/" + std::to_string(n_ops) + "/" +
                            std::to_string(i_iter) + "/k" + std::to_string(k);
          //
          // The transformation is unimodular and produces the stated form.
          //
          MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(res.Pmat);
          if (P_T * inst.GramBad * P_T.transpose() != res.GramMatRed) {
            std::cerr << "SLIDE TEST: " << tag
                      << ": the transformation does not produce the reduced "
                         "form\n";
            n_error++;
          }
          T det_P = DeterminantMat(P_T);
          if (det_P != 1 && det_P != -1) {
            std::cerr << "SLIDE TEST: " << tag
                      << ": the transformation is not unimodular, det="
                      << det_P << "\n";
            n_error++;
          }
          //
          // Both families of conditions hold on the output.
          //
          if (!IsSlideReduced<T, Tint>(res.GramMatRed, k, delta_num, delta_den,
                                       os)) {
            std::cerr << "SLIDE TEST: " << tag
                      << ": the output is not slide reduced\n";
            n_error++;
          }
          //
          // A slide reduced basis is in particular LLL reduced: the primal
          // block condition contains the Lovasz condition inside each block,
          // and the dual condition is what controls the boundaries.
          //
          if (!IsLLLreduced(res.GramMatRed, delta_num, delta_den, os)) {
            std::cerr << "SLIDE TEST: " << tag
                      << ": the output is not LLL reduced\n";
            n_error++;
          }
          if (res.GramMatRed(0, 0) < lll.GramMatRed(0, 0)) {
            n_better++;
          } else if (res.GramMatRed(0, 0) == lll.GramMatRed(0, 0)) {
            n_equal++;
          } else {
            n_worse++;
          }
          os << "SLIDE TEST: " << tag << " |b_1|^2 " << lll.GramMatRed(0, 0)
             << " (LLL) -> " << res.GramMatRed(0, 0) << "\n";
        }
        n_case++;
      }
    }
  }
  if (n_error > 0) {
    std::cerr << "SLIDE TEST: " << n_error << " errors over " << n_case
              << " cases\n";
    throw TerminalException{1};
  }
  os << "SLIDE TEST: all " << n_case << " cases pass, over the "
     << l_kmax.size() << " block sizes\n";
  os << "SLIDE TEST: against classic LLL on the same input, |b_1|^2 was "
        "shorter in "
     << n_better << " runs, equal in " << n_equal << ", longer in " << n_worse
     << "\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 1 && argc != 4) {
      std::cerr << "TEST_SlideReduction [dim] [n_iter] [seed]\n";
      std::cerr << "or\n";
      std::cerr << "TEST_SlideReduction\n";
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
    using Tint = mpz_class;
    process<T, Tint>(dim, n_iter, seed, std::cerr);
    std::cerr << "Normal termination of TEST_SlideReduction\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_SlideReduction\n";
    exit(e.eVal);
  }
  runtime(time);
}

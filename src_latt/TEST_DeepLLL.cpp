// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeReductionBench.h"
#include "DeepLLL.h"
// clang-format on

/*
  Validates the Schnorr-Euchner deep insertion.

  The check that matters is that the output satisfies the deep condition it
  claims: size reduced, and delta |b_i^*|^2 <= |pi_i(b_k)|^2 at every
  admissible pair i < k. That is verified by IsDeepLLLreduced, which
  recomputes the Gram-Schmidt data of the output from scratch and shares no
  state with the descent, so it tests the stopping rule rather than restating
  it.

  One further property is asserted: the deep condition contains the Lovasz
  condition as its case i = k-1, and also the condition at i = 1, so a deep
  reduced basis is in particular LLL reduced and the check at depth 1 must
  pass.

  The comparison of |b_1|^2 against what classic LLL returns on the same input
  is reported and NOT asserted. It is tempting to assert it, the deep
  condition being a strengthening of the Lovasz one, but that is a statement
  about the two outputs and not about the two conditions: the algorithms take
  different trajectories through GL_n(Z) and neither endpoint dominates the
  other by any argument known to the author. It happens to hold on every case
  tried here, which is worth recording as an observation and would be a bug to
  record as an invariant.

  Usage:
    TEST_DeepLLL [dim] [n_iter] [seed]
 */

template <typename T,
          typename Tint = typename underlying_z_ring<T>::ring_type>
void process(int dim, int n_iter, unsigned long seed, std::ostream &os) {
  std::mt19937_64 rng(seed);
  std::vector<std::string> l_name{"Zn", "An", "Dn", "E8", "random"};
  std::vector<int> l_nops{10, 40, 160};
  std::vector<std::pair<std::string, int>> l_method{
      {"full", 0}, {"depth5", 5}, {"depth10", 10}};
  int const delta_num = 99, delta_den = 100;
  int n_case = 0, n_error = 0;
  int n_shorter = 0, n_longer = 0;
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
        for (auto &method : l_method) {
          LLLreduction<T, Tint> res = DeepLLLreducedGeneral<T, Tint>(
              inst.GramBad, method.first, os);
          std::string tag = name + "/" + std::to_string(n_ops) + "/" +
                            std::to_string(i_iter) + "/" + method.first;
          //
          // The transformation is unimodular and produces the stated form.
          //
          MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(res.Pmat);
          if (P_T * inst.GramBad * P_T.transpose() != res.GramMatRed) {
            std::cerr << "DEEPLLL TEST: " << tag
                      << ": the transformation does not produce the reduced "
                         "form\n";
            n_error++;
          }
          T det_P = DeterminantMat(P_T);
          if (det_P != 1 && det_P != -1) {
            std::cerr << "DEEPLLL TEST: " << tag
                      << ": the transformation is not unimodular, det="
                      << det_P << "\n";
            n_error++;
          }
          //
          // The output really is deep reduced at the depth it was asked for.
          //
          if (!IsDeepLLLreduced(res.GramMatRed, method.second, delta_num,
                                delta_den, os)) {
            std::cerr << "DEEPLLL TEST: " << tag
                      << ": the output is not deep reduced\n";
            n_error++;
          }
          //
          // A deep reduced basis is LLL reduced, the Lovasz condition being
          // the case i = k-1 of the deep one. Checked at depth 0, which is
          // what "LLL reduced" means here.
          //
          if (!IsDeepLLLreduced(res.GramMatRed, 1, delta_num, delta_den, os)) {
            std::cerr << "DEEPLLL TEST: " << tag
                      << ": the output fails the Lovasz condition, which the "
                         "deep condition contains\n";
            n_error++;
          }
          //
          // Reported, not asserted; see the head of the file.
          //
          if (res.GramMatRed(0, 0) < lll.GramMatRed(0, 0)) {
            n_shorter++;
          }
          if (res.GramMatRed(0, 0) > lll.GramMatRed(0, 0)) {
            n_longer++;
          }
          os << "DEEPLLL TEST: " << tag << " |b_1|^2 " << lll.GramMatRed(0, 0)
             << " (LLL) -> " << res.GramMatRed(0, 0) << "\n";
        }
        n_case++;
      }
    }
  }
  if (n_error > 0) {
    std::cerr << "DEEPLLL TEST: " << n_error << " errors over " << n_case
              << " cases\n";
    throw TerminalException{1};
  }
  os << "DEEPLLL TEST: all " << n_case << " cases pass, over the "
     << l_method.size() << " depths\n";
  os << "DEEPLLL TEST: against classic LLL on the same input, |b_1|^2 was "
        "shorter in "
     << n_shorter << " runs and longer in " << n_longer << " out of "
     << (n_case * static_cast<int>(l_method.size()))
     << ". This is an observation and not an invariant.\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 1 && argc != 4) {
      std::cerr << "TEST_DeepLLL [dim] [n_iter] [seed]\n";
      std::cerr << "or\n";
      std::cerr << "TEST_DeepLLL\n";
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
    std::cerr << "Normal termination of TEST_DeepLLL\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_DeepLLL\n";
    exit(e.eVal);
  }
  runtime(time);
}

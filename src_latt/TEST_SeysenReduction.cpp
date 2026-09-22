// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeReductionBench.h"
#include "SeysenReduction.h"
// clang-format on

/*
  Validates the Seysen reduction.

  The decisive test is the last one: the output is checked to be a local
  minimum of Seysen's measure by brute force, every transvection (i, j, lambda)
  with |lambda| <= K being tried and required not to improve. This tests the
  closed form for the optimal coefficient rather than trusting it. If
  lambda_opt or the gain formula were wrong the descent would stop away from a
  local minimum and the scan would find the improvement that was missed. It is
  a genuinely independent check: the scan recomputes the measure from the
  inverse of the transformed Gram matrix and shares no code with the descent.

  Usage:
    TEST_SeysenReduction [dim] [n_iter] [seed]
 */

template <typename T, typename Tint>
void process(int dim, int n_iter, unsigned long seed, std::ostream &os) {
  std::mt19937_64 rng(seed);
  std::vector<std::string> l_name{"Zn", "An", "Dn", "E8", "random"};
  std::vector<int> l_nops{10, 40, 160};
  std::vector<std::string> l_method{"first", "best", "seysen_lll"};
  // The radius of the brute-force scan. Seysen steps are rounded quotients
  // and are small once the form is anywhere near reduced, so a radius of 4
  // covers every step the descent could have taken and a margin besides.
  int const K = 4;
  int n_case = 0, n_error = 0;
  for (auto &name : l_name) {
    int n_use = (name == "E8") ? 8 : dim;
    if (name == "Dn" && n_use < 4) {
      continue;
    }
    for (auto &n_ops : l_nops) {
      for (int i_iter = 0; i_iter < n_iter; i_iter++) {
        HiddenBasisInstance<T, Tint> inst =
            MakeHiddenBasisInstance<T, Tint>(name, n_use, n_ops, rng);
        int n = inst.n;
        T meas_in = SeysenMeasure(inst.GramBad);
        for (auto &method : l_method) {
          LLLreduction<T, Tint> res =
              SeysenReducedGeneral<T, Tint>(inst.GramBad, method, os);
          std::string tag = name + "/" + std::to_string(n_ops) + "/" +
                            std::to_string(i_iter) + "/" + method;
          //
          // The transformation is unimodular and is the one that produced the
          // stated reduced form.
          //
          MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(res.Pmat);
          if (P_T * inst.GramBad * P_T.transpose() != res.GramMatRed) {
            std::cerr << "SEYSEN TEST: " << tag
                      << ": the transformation does not produce the reduced "
                         "form\n";
            n_error++;
          }
          T det_P = DeterminantMat(P_T);
          if (det_P != 1 && det_P != -1) {
            std::cerr << "SEYSEN TEST: " << tag
                      << ": the transformation is not unimodular, det="
                      << det_P << "\n";
            n_error++;
          }
          //
          // The measure did not increase, and is never below the dimension.
          //
          T meas_out = SeysenMeasure(res.GramMatRed);
          if (meas_out > meas_in) {
            std::cerr << "SEYSEN TEST: " << tag
                      << ": the measure increased, from " << meas_in << " to "
                      << meas_out << "\n";
            n_error++;
          }
          if (meas_out < T(n)) {
            std::cerr << "SEYSEN TEST: " << tag << ": the measure " << meas_out
                      << " is below the dimension " << n
                      << ", which Cauchy-Schwarz forbids\n";
            n_error++;
          }
          //
          // The integral potential and the measure agree, the first being the
          // second times the determinant.
          //
          std::pair<MyMatrix<T>, T> pair = AdjugateDeterminant(res.GramMatRed);
          T pot = SeysenIntegralPotential(res.GramMatRed, pair.first);
          if (pot != meas_out * pair.second) {
            std::cerr << "SEYSEN TEST: " << tag
                      << ": the integral potential and the measure disagree\n";
            n_error++;
          }
          //
          // The output is a local minimum. This is the test that matters; see
          // the comment at the head of the file. The "seysen_lll" method is
          // excluded from it, its last pass being an LLL one that is free to
          // leave the Seysen minimum.
          //
          if (method != "seysen_lll") {
            for (int i = 0; i < n && n_error == 0; i++) {
              for (int j = 0; j < n && n_error == 0; j++) {
                if (i == j) {
                  continue;
                }
                for (int lambda = -K; lambda <= K; lambda++) {
                  if (lambda == 0) {
                    continue;
                  }
                  MyMatrix<T> U = IdentityMat<T>(n);
                  U(i, j) = T(lambda);
                  MyMatrix<T> Gnew = U * res.GramMatRed * U.transpose();
                  T meas_new = SeysenMeasure(Gnew);
                  if (meas_new < meas_out) {
                    std::cerr << "SEYSEN TEST: " << tag
                              << ": the output is not a local minimum, the "
                                 "move (i="
                              << i << ", j=" << j << ", lambda=" << lambda
                              << ") lowers the measure from " << meas_out
                              << " to " << meas_new << "\n";
                    n_error++;
                    break;
                  }
                }
              }
            }
          }
          os << "SEYSEN TEST: " << tag << " measure " << meas_in << " -> "
             << meas_out << "\n";
        }
        n_case++;
      }
    }
  }
  if (n_error > 0) {
    std::cerr << "SEYSEN TEST: " << n_error << " errors over " << n_case
              << " cases\n";
    throw TerminalException{1};
  }
  os << "SEYSEN TEST: all " << n_case << " cases pass, over the " << l_method.size()
     << " methods\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 1 && argc != 4) {
      std::cerr << "TEST_SeysenReduction [dim] [n_iter] [seed]\n";
      std::cerr << "or\n";
      std::cerr << "TEST_SeysenReduction\n";
      std::cerr << "for the default dim=6, n_iter=2, seed=1\n";
      return -1;
    }
    int dim = 6;
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
    std::cerr << "Normal termination of TEST_SeysenReduction\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_SeysenReduction\n";
    exit(e.eVal);
  }
  runtime(time);
}

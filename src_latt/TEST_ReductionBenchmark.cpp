// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeReductionBench.h"
#include "BKZ.h"
#include "SlideReduction.h"
#include "DeepLLL.h"
// clang-format on

/*
  The hidden-good-basis benchmark of the KAN reduction study.

  For each lattice in the list and each destruction strength n_ops, a random
  unimodular matrix is applied to a known-good Gram matrix and the reducers are
  run on the result. What is printed per case is, for each algorithm, whether
  the hidden presentation was recovered up to signed permutation, whether the
  reduced form is at least as good as the hidden one, and the quality measures
  themselves.

  The point of the harness is that it separates four choices that are easy to
  confuse when an experiment disappoints: the representation, the objective,
  the move set and the search strategy. A reducer that loses here has lost on
  one of them, and the measures are reported separately so that one can tell
  which.

  Usage:
    TEST_ReductionBenchmark [dim] [n_iter] [seed]

  with dim the ambient dimension used for the families Zn, An and Dn (E8 is
  fixed at 8), n_iter the number of random instances per case, and seed the
  generator seed, so that a reported failure can be replayed exactly.
 */

template <typename T, typename Tint>
void process(int dim, int n_iter, unsigned long seed, std::ostream &os) {
  std::mt19937_64 rng(seed);
  std::vector<std::string> l_name{"Zn", "An", "Dn", "E8", "random"};
  std::vector<int> l_nops{10, 40, 160};
  //
  // The reducers under test. LLLreducedBasis is the baseline that everything
  // else has to beat; the dual variant is included because
  // SublatticeBasisReductionKernel already alternates the two, which is an
  // ad hoc form of a primal-dual reduction and worth measuring on its own.
  //
  using Freduce = std::function<LLLreduction<T, Tint>(MyMatrix<T> const &,
                                                      std::ostream &)>;
  std::vector<std::pair<std::string, Freduce>> l_algo;
  l_algo.push_back({"none", [](MyMatrix<T> const &G,
                               [[maybe_unused]] std::ostream &os_i) {
                      return LLLnoreduction<T, Tint>(G);
                    }});
  l_algo.push_back({"LLL", [](MyMatrix<T> const &G, std::ostream &os_i) {
                      return LLLreducedBasis<T, Tint>(G, os_i);
                    }});
  l_algo.push_back({"LLLdual", [](MyMatrix<T> const &G, std::ostream &os_i) {
                      return LLLreducedBasisDual<T, Tint>(G, os_i);
                    }});
  l_algo.push_back({"Deep", [](MyMatrix<T> const &G, std::ostream &os_i) {
                      return DeepLLLreducedBasis<T, Tint>(G, os_i);
                    }});
  l_algo.push_back({"Deep5", [](MyMatrix<T> const &G, std::ostream &os_i) {
                      return DeepLLLreducedBasisDepth<T, Tint>(G, 5, os_i);
                    }});
  l_algo.push_back({"BKZ4", [](MyMatrix<T> const &G, std::ostream &os_i) {
                      return BKZreducedBasis<T, Tint>(G, 4, os_i);
                    }});
  l_algo.push_back({"BKZ8", [](MyMatrix<T> const &G, std::ostream &os_i) {
                      return BKZreducedBasis<T, Tint>(G, 8, os_i);
                    }});
  l_algo.push_back({"BKZ12", [](MyMatrix<T> const &G, std::ostream &os_i) {
                      return BKZreducedBasis<T, Tint>(G, 12, os_i);
                    }});
  l_algo.push_back({"Slide4", [](MyMatrix<T> const &G, std::ostream &os_i) {
                      return SlideReducedBasisAuto<T, Tint>(G, 4, os_i);
                    }});
  l_algo.push_back({"Slide8", [](MyMatrix<T> const &G, std::ostream &os_i) {
                      return SlideReducedBasisAuto<T, Tint>(G, 8, os_i);
                    }});
  l_algo.push_back({"Seysen", [](MyMatrix<T> const &G, std::ostream &os_i) {
                      return SeysenReducedBasis<T, Tint>(G, os_i);
                    }});
  l_algo.push_back({"SeysenB", [](MyMatrix<T> const &G, std::ostream &os_i) {
                      return SeysenReducedBasisBest<T, Tint>(G, os_i);
                    }});
  l_algo.push_back({"SeyLLL", [](MyMatrix<T> const &G, std::ostream &os_i) {
                      return SeysenLLLreducedBasis<T, Tint>(G, os_i);
                    }});
  //
  int n_case = 0;
  std::map<std::string, int> n_recovered, n_good_quality;
  std::map<std::string, double> total_ms;
  for (auto &name : l_name) {
    int n_use = (name == "E8") ? 8 : dim;
    if (name == "Dn" && n_use < 4) {
      continue;
    }
    for (auto &n_ops : l_nops) {
      for (int i_iter = 0; i_iter < n_iter; i_iter++) {
        HiddenBasisInstance<T, Tint> inst =
            MakeHiddenBasisInstance<T, Tint>(name, n_use, n_ops, rng);
        ReductionQuality<T> q_good = ComputeReductionQuality(inst.GramGood);
        os << "case " << name << " n=" << inst.n << " n_ops=" << n_ops
           << " iter=" << i_iter << "\n";
        PrintReductionQuality(os, "hidden", q_good);
        for (auto &algo : l_algo) {
          ReductionOutcome<T> out = RunOneReduction<T, Tint>(
              inst, algo.first, algo.second, os);
          std::string label = algo.first;
          label += std::string(8 - std::min<size_t>(8, label.size()), ' ');
          label += (out.recovered ? " rec=yes" : " rec=no ");
          label += (out.matches_good_quality ? " qual=ok " : " qual=BAD");
          PrintReductionQuality(os, label, out.quality);
          if (out.recovered) {
            n_recovered[algo.first]++;
          }
          if (out.matches_good_quality) {
            n_good_quality[algo.first]++;
          }
          total_ms[algo.first] += out.runtime_ms;
        }
        n_case++;
      }
    }
  }
  os << "\nREDUCTION_BENCH: summary over " << n_case << " cases\n";
  for (auto &algo : l_algo) {
    os << "  " << algo.first << ": recovered=" << n_recovered[algo.first]
       << "/" << n_case << " quality_ok=" << n_good_quality[algo.first] << "/"
       << n_case << " total_time=" << total_ms[algo.first] << " ms\n";
  }
  //
  // The one assertion of the harness. Everything above is measurement, to be
  // read rather than asserted, but a reducer that does not reduce at all is a
  // bug and not a disappointing experiment.
  //
  if (n_good_quality["LLL"] == 0 && n_case > 0) {
    std::cerr << "REDUCTION_BENCH: LLL never reached the quality of the "
                 "hidden basis, this is not a plausible outcome\n";
    throw TerminalException{1};
  }
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 1 && argc != 4) {
      std::cerr << "TEST_ReductionBenchmark [dim] [n_iter] [seed]\n";
      std::cerr << "or\n";
      std::cerr << "TEST_ReductionBenchmark\n";
      std::cerr << "for the default dim=8, n_iter=3, seed=1\n";
      return -1;
    }
    int dim = 8;
    int n_iter = 3;
    unsigned long seed = 1;
    if (argc == 4) {
      dim = std::stoi(argv[1]);
      n_iter = std::stoi(argv[2]);
      seed = std::stoul(argv[3]);
    }
    using T = mpq_class;
    using Tint = mpz_class;
    process<T, Tint>(dim, n_iter, seed, std::cerr);
    std::cerr << "Normal termination of TEST_ReductionBenchmark\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_ReductionBenchmark\n";
    exit(e.eVal);
  }
  runtime(time);
}

// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "LocalDelaunayConfig.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

/*
  Whether a configuration of unimodular Delaunay simplices can be extended,
  facet by facet, to depth_limit more simplices for a form of covering
  density (or packing-covering constant) at most the record. See
  ExtendConfiguration in LocalDelaunayConfig.h.
 */

// The neighbours of sigma_0 across the facet opposite to the origin that
// survive the pair bound, all of them (not one per orbit).
template <typename T, typename Tint>
std::vector<MyVector<T>> ComputeP0(int n, bool is_pc, double threshold,
                                   double dstar) {
  using namespace local_delaunay;
  MyMatrix<T> Sigma0 = ZeroMatrix<T>(n + 1, n);
  for (int i = 0; i < n; i++) {
    Sigma0(i + 1, i) = T(1);
  }
  SimplexData<T> sd = GetSimplexData(Sigma0);
  std::vector<MyVector<T>> P0;
  std::optional<double> opt_lambda = InertiaLowerBound(n, 1, dstar);
  if (!opt_lambda) {
    return P0;
  }
  double EdgeBound = 4.0 / *opt_lambda;
  std::vector<NeighborCandidate<T>> ListCand =
      EnumerateNeighborCandidates<T, Tint>(sd, {0}, EdgeBound, 1, std::cerr);
  std::vector<AffineMap<T>> ListMap = SimplexStabilizer(sd);
  std::vector<size_t> ListRep = OrbitRepresentatives(ListCand, ListMap);
  std::vector<std::optional<PairBound>> ListBound(ListCand.size());
  for (size_t i = 0; i < ListCand.size(); i++) {
    size_t r = ListRep[i];
    if (!ListBound[r]) {
      ListBound[r] =
          is_pc ? BoundPairPackCov<T, Tint>(Sigma0, ListCand[r], false, 3,
                                            std::cerr)
                : BoundPair<T, Tint>(Sigma0, ListCand[r], false, 0, std::cerr);
    }
    if (!ListBound[r]->success || ListBound[r]->bound >= threshold) {
      P0.push_back(ListCand[i].w);
    }
  }
  return P0;
}

template <typename T, typename Tint>
void process(std::string const &FileConfig, std::string const &mode, double p1,
             double p2, int Vmax, int depth_limit, size_t node_limit,
             std::ostream &os_out) {
  using namespace local_delaunay;
  std::ifstream is(FileConfig);
  int n_simp;
  is >> n_simp;
  std::vector<MyMatrix<T>> ListStart;
  for (int i = 0; i < n_simp; i++) {
    ListStart.push_back(ReadMatrix<T>(is));
  }
  if (ListStart.empty()) {
    std::cerr << "The configuration file FileConfig=" << FileConfig
              << " has no simplex\n";
    throw TerminalException{1};
  }
  int n = ListStart[0].cols();
  bool is_pc = mode == "pc";
  if (mode != "cov" && mode != "pc") {
    std::cerr << "mode should be cov or pc\n";
    throw TerminalException{1};
  }
  double ThetaRec = is_pc ? p2 * std::pow(p1, n) : p1;
  double lambda_min = is_pc ? 4.0 / (p1 * p1) : 0;
  double dstar = DeterminantThreshold(n, ThetaRec);
  double threshold = is_pc ? lambda_min : dstar;
  // Vmax = 0: all the simplices are unimodular and the neighbours are
  // given by P0, the neighbours of sigma_0 across the facet opposite to the
  // origin that survive the pair bound. Vmax > 0: the neighbours of each
  // simplex, of volume at most Vmax, are enumerated (general mode).
  std::vector<MyVector<T>> P0;
  if (Vmax == 0) {
    P0 = ComputeP0<T, Tint>(n, is_pc, threshold, dstar);
  }
  os_out << "n=" << n << " mode=" << mode << " threshold=" << threshold
         << " Vmax=" << Vmax << " |P0|=" << P0.size()
         << " |start|=" << ListStart.size() << std::endl;
  HumanTime time;
  ExtendResult res =
      ExtendConfiguration<T, Tint>(ListStart, P0, is_pc, threshold, depth_limit,
                                   node_limit, Vmax, dstar, os_out, std::cerr);
  os_out << "n_node=" << res.n_node << " n_bound=" << res.n_bound
         << " max_depth=" << res.max_depth << " complete=" << res.complete
         << " time=" << time << "\n";
  if (res.refuted && res.complete) {
    os_out << "REFUTED\n";
  } else if (!res.refuted) {
    os_out << "NOT REFUTED: a branch reached depth " << depth_limit
           << " or closed\n";
  } else {
    os_out << "INCOMPLETE\n";
  }
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    if (argc != 9 && argc != 10) {
      std::cerr
          << "LATT_LocalDelaunayExtend [arith] [FileConfig] [cov|pc] [p1] [p2] "
             "[Vmax] [depth_limit] [node_limit] [OutFile]\n";
      std::cerr << "Vmax = 0: unimodular simplices only. Vmax > 0: "
                   "neighbours of volume at most Vmax\n";
      std::cerr << "FileConfig: the number of simplices, then their vertex "
                   "matrices\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileConfig = argv[2];
    std::string mode = argv[3];
    double p1 = ParseScalar<double>(argv[4]);
    double p2 = ParseScalar<double>(argv[5]);
    int Vmax = ParseScalar<int>(argv[6]);
    int depth_limit = ParseScalar<int>(argv[7]);
    size_t node_limit = ParseScalar<size_t>(argv[8]);
    std::string OutFile = "stdout";
    if (argc == 10) {
      OutFile = argv[9];
    }
    auto f = [&](std::ostream &os) -> void {
      ArithmeticDispatch dispatch;
      add_rational_integer_arithmetics_by_library<false>(
          dispatch, [&]<typename T, typename Tint>() -> void {
            process<T, Tint>(FileConfig, mode, p1, p2, Vmax, depth_limit,
                             node_limit, os);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of LATT_LocalDelaunayExtend\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_LocalDelaunayExtend\n";
    exit(e.eVal);
  }
  runtime(time1);
}

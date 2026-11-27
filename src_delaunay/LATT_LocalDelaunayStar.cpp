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
  Whether a lattice simplex S can occur in the Delaunay triangulation of a
  form of covering density at most ThetaRec: the neighbours of S across each
  facet are listed and filtered by the pair bound, then the stars of S (one
  neighbour per facet) are searched. See LocalDelaunayConfig.h.
 */

template <typename T, typename Tint>
void process(std::string const &FileS, std::string const &mode, double p1,
             double p2, int Vmax, size_t node_limit, size_t max_print,
             std::ostream &os_out) {
  using namespace local_delaunay;
  MyMatrix<T> Vert = ReadMatrixFile<T>(FileS);
  int n = Vert.cols();
  if (Vert.rows() != n + 1) {
    std::cerr << "The file FileS=" << FileS
              << " should have n+1 rows, the vertices of the simplex\n";
    throw TerminalException{1};
  }
  SimplexData<T> sd = GetSimplexData(Vert);
  // cov: p1 = Theta_rec. pc: p1 = gamma_rec, p2 = delta_n, and the
  // threshold is on lambda for R = 1 (see LATT_LocalDelaunayNeighbors).
  bool is_pc = mode == "pc";
  if (mode != "cov" && mode != "pc") {
    std::cerr << "mode should be cov or pc\n";
    throw TerminalException{1};
  }
  double ThetaRec = is_pc ? p2 * std::pow(p1, n) : p1;
  double lambda_min = is_pc ? 4.0 / (p1 * p1) : 0;
  double dstar = DeterminantThreshold(n, ThetaRec);
  double threshold = is_pc ? lambda_min : dstar;
  std::optional<double> opt_lambda = InertiaLowerBound(n, sd.V, dstar);
  os_out << "n=" << n << " V=" << sd.V << " dstar=" << dstar << "\n";
  if (!opt_lambda) {
    os_out << "EXCLUDED: no simplex of this volume below the threshold\n";
    return;
  }
  double EdgeBound = 4.0 / *opt_lambda;
  HumanTime time;
  std::vector<int> ListFacet;
  for (int k = 0; k <= n; k++) {
    ListFacet.push_back(k);
  }
  std::vector<NeighborCandidate<T>> ListCand =
      EnumerateNeighborCandidates<T, Tint>(sd, ListFacet, EdgeBound, Vmax,
                                           std::cerr);
  std::vector<AffineMap<T>> ListMap = SimplexStabilizer(sd);
  std::vector<size_t> ListRep = OrbitRepresentatives(ListCand, ListMap);
  std::vector<std::optional<PairBound>> ListBound(ListCand.size());
  std::vector<std::vector<StarCandidate<T>>> ListCandFacet(n + 1);
  for (size_t i_cand = 0; i_cand < ListCand.size(); i_cand++) {
    NeighborCandidate<T> const &cand = ListCand[i_cand];
    size_t i_rep = ListRep[i_cand];
    if (!ListBound[i_rep]) {
      ListBound[i_rep] =
          is_pc
              ? BoundPairPackCov<T, Tint>(Vert, ListCand[i_rep], false, 3,
                                          std::cerr)
              : BoundPair<T, Tint>(Vert, ListCand[i_rep], false, 0, std::cerr);
    }
    PairBound pb = *ListBound[i_rep];
    if (pb.success && pb.bound < threshold) {
      continue;
    }
    MyMatrix<T> Vert2 = Vert;
    for (int j = 0; j < n; j++) {
      Vert2(cand.k, j) = cand.w(j);
    }
    ListCandFacet[cand.k].push_back(
        {cand.k, cand.w, cand.beta, Vert2, pb.bound, {}, {}});
  }
  os_out << "Pair filter done, time=" << time << std::endl;
  int root_facet = 0;
  for (int k = 0; k <= n; k++) {
    os_out << "facet " << k << ": " << ListCandFacet[k].size()
           << " neighbours after the pair bound\n";
    if (ListCandFacet[k].size() < ListCandFacet[root_facet].size()) {
      root_facet = k;
    }
  }
  if (ListCandFacet[root_facet].empty()) {
    os_out << "EXCLUDED: facet " << root_facet << " has no neighbour\n";
    return;
  }
  // The root candidates up to the stabilizer of S fixing the root facet.
  os_out << "|Stab(S)|=" << ListMap.size() << "\n";
  std::vector<StarCandidate<T>> RootList;
  std::set<std::vector<T>> SetSeen;
  for (auto &c : ListCandFacet[root_facet]) {
    std::vector<T> canon;
    for (auto &am : ListMap) {
      if (am.perm[root_facet] != root_facet) {
        continue;
      }
      MyVector<T> img = ApplyAffine(am, c.w);
      std::vector<T> key(img.data(), img.data() + n);
      if (canon.empty() || key < canon) {
        canon = key;
      }
    }
    if (SetSeen.insert(canon).second) {
      RootList.push_back(c);
    }
  }
  os_out << "root facet " << root_facet << ": " << RootList.size()
         << " orbits of neighbours" << std::endl;
  StarSearchResult res =
      StarSearch<T, Tint>(sd, ListCandFacet, RootList, root_facet, threshold,
                          is_pc, node_limit, max_print, os_out, std::cerr);
  os_out << "n_node=" << res.n_node << " n_maxdet=" << res.n_maxdet
         << " n_full=" << res.n_full << " complete=" << res.complete
         << " min_full_bound=" << res.min_full_bound
         << " max_full_bound=" << res.max_full_bound
         << " threshold=" << threshold << " time=" << time << "\n";
  if (res.complete && res.n_full == 0) {
    os_out << "EXCLUDED: no star of S survives\n";
  } else if (res.complete) {
    os_out << "NOT EXCLUDED: " << res.n_full << " stars survive\n";
  } else {
    os_out << "INCOMPLETE: node limit reached\n";
  }
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    if (argc != 9 && argc != 10) {
      std::cerr
          << "LATT_LocalDelaunayStar [arith] [FileS] [cov|pc] [p1] [p2] [Vmax] "
             "[node_limit] [max_print] [OutFile]\n";
      std::cerr << "cov: p1 = Theta_rec, p2 unused. pc: p1 = gamma_rec, p2 = "
                   "delta_n\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileS = argv[2];
    std::string mode = argv[3];
    double p1 = ParseScalar<double>(argv[4]);
    double p2 = ParseScalar<double>(argv[5]);
    int Vmax = ParseScalar<int>(argv[6]);
    size_t node_limit = ParseScalar<size_t>(argv[7]);
    size_t max_print = ParseScalar<size_t>(argv[8]);
    std::string OutFile = "stdout";
    if (argc == 10) {
      OutFile = argv[9];
    }
    auto f = [&](std::ostream &os) -> void {
      ArithmeticDispatch dispatch;
      add_rational_integer_arithmetics_by_library<false>(
          dispatch, [&]<typename T, typename Tint>() -> void {
            process<T, Tint>(FileS, mode, p1, p2, Vmax, node_limit, max_print,
                             os);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of LATT_LocalDelaunayStar\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_LocalDelaunayStar\n";
    exit(e.eVal);
  }
  runtime(time1);
}

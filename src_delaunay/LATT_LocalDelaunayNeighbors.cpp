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
  For a lattice simplex S, the possible neighbours of S across its facets in
  the Delaunay triangulation of a form with covering density at most
  ThetaRec, each with the determinant maximization bound of the pair.
  See LocalDelaunayConfig.h.
 */

template <typename T, typename Tint>
void process(std::string const &FileS, std::string const &mode, double p1,
             double p2, int Vmax, std::string const &FacetChoice,
             bool emptiness, int max_cut, std::ostream &os_out) {
  using namespace local_delaunay;
  MyMatrix<T> Vert = ReadMatrixFile<T>(FileS);
  int n = Vert.cols();
  if (Vert.rows() != n + 1) {
    std::cerr << "The file FileS=" << FileS
              << " should have n+1 rows, the vertices of the simplex\n";
    throw TerminalException{1};
  }
  SimplexData<T> sd = GetSimplexData(Vert);
  if (sd.V == 0) {
    std::cerr << "The simplex is degenerate\n";
    throw TerminalException{1};
  }
  // For the covering p1 is Theta_rec. For the packing-covering p1 is
  // gamma_rec and p2 the optimal lattice packing density delta_n, so that
  // Theta(Q) = delta(Q) gamma(Q)^n <= delta_n gamma_rec^n and, with R = 1,
  // lambda(Q) >= 4 / gamma_rec^2.
  bool is_pc = mode == "pc";
  if (mode != "cov" && mode != "pc") {
    std::cerr << "mode should be cov or pc\n";
    throw TerminalException{1};
  }
  double ThetaRec = is_pc ? p2 * std::pow(p1, n) : p1;
  double lambda_min = is_pc ? 4.0 / (p1 * p1) : 0;
  double dstar = DeterminantThreshold(n, ThetaRec);
  std::optional<double> opt_lambda = InertiaLowerBound(n, sd.V, dstar);
  os_out << "n=" << n << " V=" << sd.V << " mode=" << mode
         << " ThetaRec=" << ThetaRec << " dstar=" << dstar;
  if (is_pc) {
    os_out << " lambda_min=" << lambda_min;
  }
  os_out << "\n";
  if (!opt_lambda) {
    os_out << "No simplex of volume " << sd.V
           << " is possible below the threshold\n";
    return;
  }
  double lambda = *opt_lambda;
  double EdgeBound = 4.0 / lambda;
  os_out << "lambda_S > " << lambda << " EdgeBound=" << EdgeBound << "\n";
  std::vector<int> ListFacet;
  if (FacetChoice == "all") {
    for (int k = 0; k <= n; k++) {
      ListFacet.push_back(k);
    }
  } else {
    ListFacet.push_back(ParseScalar<int>(FacetChoice));
  }
  HumanTime time;
  std::vector<NeighborCandidate<T>> ListCand =
      EnumerateNeighborCandidates<T, Tint>(sd, ListFacet, EdgeBound, Vmax,
                                           std::cerr);
  os_out << "|ListCand|=" << ListCand.size() << " time=" << time << "\n";
  // One bound per orbit of the stabilizer of S, the counts being over
  // all the candidates.
  std::vector<AffineMap<T>> ListMap = SimplexStabilizer(sd);
  std::vector<size_t> ListRep = OrbitRepresentatives(ListCand, ListMap);
  os_out << "|Stab(S)|=" << ListMap.size() << " time=" << time << std::endl;
  std::vector<std::optional<PairBound>> ListBound(ListCand.size());
  std::map<int, std::pair<size_t, size_t>> map_facet;
  std::map<int, double> map_facet_max;
  size_t n_solved = 0;
  for (size_t i_cand = 0; i_cand < ListCand.size(); i_cand++) {
    NeighborCandidate<T> const &cand = ListCand[i_cand];
    size_t i_rep = ListRep[i_cand];
    bool is_rep = i_rep == i_cand;
    if (!ListBound[i_rep]) {
      ListBound[i_rep] =
          is_pc ? BoundPairPackCov<T, Tint>(Vert, ListCand[i_rep], emptiness,
                                            max_cut, std::cerr)
                : BoundPair<T, Tint>(Vert, ListCand[i_rep], emptiness, max_cut,
                                     std::cerr);
      n_solved++;
    }
    PairBound pb = *ListBound[i_rep];
    // For the packing-covering the bound is on lambda, to be compared with
    // 4 / gamma_rec^2.
    double threshold = is_pc ? lambda_min : dstar;
    bool accepted = !pb.success || pb.bound >= threshold;
    map_facet[cand.k].first++;
    if (accepted) {
      map_facet[cand.k].second++;
      if (map_facet_max.count(cand.k) == 0 ||
          pb.bound > map_facet_max[cand.k]) {
        map_facet_max[cand.k] = pb.bound;
      }
    }
    if (!is_rep) {
      continue;
    }
    os_out << "k=" << cand.k << " w=";
    for (int i = 0; i < n; i++) {
      os_out << " " << cand.w(i);
    }
    os_out << " V'=" << cand.V << " det=" << pb.det << " bound=" << pb.bound
           << " n_cut=" << pb.n_cut << " success=" << pb.success << " "
           << (accepted ? "ACCEPTED" : "rejected") << "\n";
  }
  os_out << "Summary: n_solved=" << n_solved << " time=" << time << "\n";
  for (auto &kv : map_facet) {
    os_out << "facet " << kv.first << ": candidates=" << kv.second.first
           << " accepted=" << kv.second.second;
    if (map_facet_max.count(kv.first) > 0) {
      os_out << " max_bound=" << map_facet_max[kv.first];
    }
    os_out << "\n";
  }
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    if (argc != 10 && argc != 11) {
      std::cerr
          << "LATT_LocalDelaunayNeighbors [arith] [FileS] [cov|pc] [p1] [p2] "
             "[Vmax] [facet|all] [emptiness 0|1] [max_cut] [OutFile]\n";
      std::cerr << "cov: p1 = Theta_rec, p2 unused. pc: p1 = gamma_rec, p2 = "
                   "delta_n\n";
      std::cerr << "FileS: the n+1 vertices of a lattice simplex, one per "
                   "row\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileS = argv[2];
    std::string mode = argv[3];
    double p1 = ParseScalar<double>(argv[4]);
    double p2 = ParseScalar<double>(argv[5]);
    int Vmax = ParseScalar<int>(argv[6]);
    std::string FacetChoice = argv[7];
    bool emptiness = ParseScalar<int>(argv[8]) == 1;
    int max_cut = ParseScalar<int>(argv[9]);
    std::string OutFile = "stdout";
    if (argc == 11) {
      OutFile = argv[10];
    }
    auto f = [&](std::ostream &os) -> void {
      ArithmeticDispatch dispatch;
      add_rational_integer_arithmetics_by_library<false>(
          dispatch, [&]<typename T, typename Tint>() -> void {
            process<T, Tint>(FileS, mode, p1, p2, Vmax, FacetChoice, emptiness,
                             max_cut, os);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of LATT_LocalDelaunayNeighbors\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_LocalDelaunayNeighbors\n";
    exit(e.eVal);
  }
  runtime(time1);
}

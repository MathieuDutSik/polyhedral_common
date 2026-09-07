// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_COVERINGRECORDSEARCH_H_
#define SRC_DELAUNAY_COVERINGRECORDSEARCH_H_

// clang-format off
#include "CoveringMaxdet.h"
#include "PeriodicDelaunay.h"
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>
// clang-format on

/*
  Random-walk search for an iso-Delaunay domain whose covering optimum beats
  a record.

  The covering density of a point set is the minimum, over the iso-Delaunay
  domains of its parameter space, of the determinant maximization optimum of
  CoveringMaxdet.h over each domain. Enumerating every domain answers the
  question outright, but the number of domains explodes with the dimension
  and with the coset structure, so for the cases of interest it is out of
  reach.

  What is in reach is a walk. The strategy is the one of
  LookForFullRankRayDomain (itself modelled on
  src_ctype/CTYP_LookForNoFreeVector.cpp), which is how the record covering
  in dimension 6 was found:

  * evaluate the objective on the current domain and on every adjacent one;
  * if some adjacent domain is strictly better, move to a uniformly random
    one among those tying for the best;
  * otherwise the walk sits in a local minimum, so take n_walk_steps
    unbiased random jumps in the adjacency graph to leave the basin;
  * stop on a hit, on the runtime budget, or when the graph runs out.

  The objective here is the covering density optimum over the domain, which
  is what has to be pushed below the record.

  --- Failure is the expected outcome ---

  In dimensions 3 to 5 the best known covering is a lattice one, and it is
  conjectured to be the best covering altogether, so a periodic point set is
  not expected to beat it. The search therefore treats "no record" as a
  normal outcome: it returns the best density it saw, together with the
  domain attaining it, and the caller reports that rather than failing. A
  walk that runs out of budget, that reaches a domain with no flippable wall
  or whose optimization does not converge is likewise reported, not thrown.

  --- What a "record" means here ---

  The optimization is numerical (see CoveringMaxdet.h), so a density that
  sits within a few units of the last digit of the record proves nothing.
  GetKnownLatticeCoveringRecord returns the exact closed form for the A_n^*
  of dimensions up to 5 and the published decimal beyond, which is rounded;
  a hit that close to the record has to be confirmed by hand.
 */

#ifdef DEBUG
#define DEBUG_COVERING_RECORD
#endif

#ifdef TIMINGS
#define TIMINGS_COVERING_RECORD
#endif

namespace covering_record {

// ---------------------------------------------------------------------------
// The records to beat
// ---------------------------------------------------------------------------

/*
  The volume kappa_n of the unit ball, taken from the covering density
  computation itself so that the record and the densities compared with it
  are built on the very same constant.
 */
inline double GetVolumeUnitBall(int n) {
  ResultCov<double> rc =
      ComputeCoveringDensityFromDimDetCov<double>(n, 1.0, 1.0);
  return rc.VolumeBall;
}

/*
  The covering density of the dual root lattice A_n^*, in closed form:

      Theta(A_n^*) = kappa_n sqrt(n+1) ( n(n+2) / (12(n+1)) )^{n/2}.

  A_n^* is the best known lattice covering for n <= 5, and the best possible
  one for n <= 5 by the classification of the iso-Delaunay domains.
 */
inline double GetCoveringDensityAnStar(int n) {
  double dn = static_cast<double>(n);
  double base = dn * (dn + 2) / (12 * (dn + 1));
  return GetVolumeUnitBall(n) * std::sqrt(dn + 1) * std::pow(base, dn / 2);
}

/*
  The least dense known lattice covering of dimension n, which is what a
  periodic covering has to beat to be a record.

  Up to dimension 5 this is A_n^* and the closed form above is used, so the
  value is exact to the last bit a double carries. From dimension 6 on it is
  the L^c_n of Table 2 of Dutour Sikiric, Schuermann and Vallentin, Duke
  Math. J. 142 (2008), whose entries are rounded decimals -- in particular
  the dimension 6 record 2.464801 is the lattice L^c_6 of Vallentin. Beyond
  dimension 8 nothing is returned and the caller has to state the record.
 */
inline std::optional<double> GetKnownLatticeCoveringRecord(int n) {
  if (n >= 1 && n <= 5) {
    return GetCoveringDensityAnStar(n);
  }
  static const std::map<int, double> map_record{
      {6, 2.464801}, {7, 2.900024}, {8, 3.142202}};
  std::map<int, double>::const_iterator iter = map_record.find(n);
  if (iter == map_record.end()) {
    return {};
  }
  return iter->second;
}

// ---------------------------------------------------------------------------
// The objective
// ---------------------------------------------------------------------------

/*
  The covering optimum over one iso-Delaunay domain: the whole
  CoveringMaxdet.h pipeline, from the defining inequalities of the domain to
  the barrier solve.

  A solve that does not reach the tolerance is kept: it still returns a
  strictly feasible point, whose density is a covering density the point set
  attains and therefore an upper bound on the optimum of the domain. That is
  what the walk compares and what a record claim rests on, so throwing those
  away would both waste the work and risk discarding the very domain that
  holds a record. Only a domain that yields no feasible point at all is
  dropped, and the caller treats that as "this domain is not usable" rather
  than as an error.

  How often that happens depends sharply on the point set: it is zero on
  most configurations and 95% on some, when the forms interior to the
  domains are anisotropic enough that a double cannot hold the circumradius
  blocks of the different orbits at once (see GetStartingPoint). The count
  is reported as n_domain_failed, and a result carrying a large one is the
  best over a small part of what the walk visited rather than over all of
  it.
 */
template <typename T, typename Tint, typename Tgroup>
std::optional<covering_maxdet::MaxdetResult<double>>
GetDomainCoveringOptimum(IsoDelaunayDomain<T, Tint, Tgroup> const &x,
                         LinSpaceMatrix<T> const &LinSpa,
                         T const &point_density, std::ostream &os) {
#ifdef TIMINGS_COVERING_RECORD
  MicrosecondTime time;
#endif
  std::vector<FullAdjInfo<Tint>> ListIneq =
      ComputeListIneqFromTesselationIneq(x.DT);
  MyMatrix<T> FAC = UniversalMatrixConversion<T, Tint>(GetFACineq(ListIneq));
  std::vector<int> ListIrred = get_non_redundant_indices(FAC, os);
  MyMatrix<T> FACred = SelectRow(FAC, ListIrred);
  covering_maxdet::CoveringData<T> cd = covering_maxdet::BuildCoveringData(
      x, LinSpa, FACred, point_density, os);
  covering_maxdet::CoveringData<double> cd_f =
      covering_maxdet::ConvertCoveringData<double, T>(cd);
  std::optional<MyVector<double>> opt_start =
      covering_maxdet::GetStartingPoint<T, double>(cd, cd_f, os);
  if (!opt_start) {
    os << "COVERING_RECORD: the domain that could not be optimized has a Gram "
          "matrix of largest entry "
       << x.GramMat.cwiseAbs().maxCoeff() << "\n";
    return {};
  }

  covering_maxdet::MaxdetOptions<double> opts =
      covering_maxdet::GetDefaultMaxdetOptions<double>();
  covering_maxdet::MaxdetResult<double> res =
      covering_maxdet::SolveCoveringMaxdet(cd_f, *opt_start, opts, os);
#ifdef TIMINGS_COVERING_RECORD
  os << "|COVERING_RECORD: GetDomainCoveringOptimum|=" << time << "\n";
#endif
  if (!res.has_point) {
    return {};
  }
  return res;
}

/*
  How skewed a representative is: the largest entry of its Gram matrix. The
  domains of a point set have bounded coordinates, so a value far above what
  a fresh domain shows is the drift of the walk rather than a property of
  the domain.
 */
template <typename T, typename Tint, typename Tgroup>
double GetDomainSkew(IsoDelaunayDomain<T, Tint, Tgroup> const &x) {
  Tint val = x.GramMat.cwiseAbs().maxCoeff();
  return UniversalScalarConversion<double, Tint>(val);
}

// ---------------------------------------------------------------------------
// The search
// ---------------------------------------------------------------------------

struct RecordSearchOptions {
  // The covering density to beat: a domain whose optimum is strictly below
  // this is a hit and ends the search.
  double record;
  // The number of unbiased random jumps taken to leave a local minimum.
  int n_walk_steps;
  // The wall-clock budget in seconds; 0 means no limit.
  int max_runtime_second;
  // The walk is restarted from a fresh domain once the Gram matrix of the
  // current one has an entry above this; 0 disables the guard. See the
  // comment on the drift in LookForRecordCovering.
  double max_gram_entry;
  // The directory where the Gram matrix of every improvement is written.
  // Empty writes nothing.
  std::string Prefix;
};

template <typename Tint> struct RecordSearchResult {
  // Whether a domain beating the record was found. False is the expected
  // outcome in the dimensions where the lattice covering is optimal.
  bool found_record = false;
  // Whether any domain at all could be optimized.
  bool has_best = false;
  // The best density seen, the Gram matrix of the domain attaining it (an
  // interior point of that domain, over the ring) and the optimizer itself.
  double best_density = 0;
  MyMatrix<Tint> best_domain_gram;
  MyMatrix<double> best_gram;
  double best_covering_radius_sq = 0;
  // How the search ended, and what it cost.
  std::string message;
  int n_iter = 0;
  int n_walk = 0;
  int n_restart = 0;
  int n_domain_evaluated = 0;
  int n_domain_failed = 0;
  int runtime_second = 0;
};

/*
  Take n_iter unbiased random steps in the adjacency graph, the escape move
  of the walk. Same as RandomWalkIsoDelaunay of IsoDelaunayDomains.h, except
  that the generators of the stabilizer are supplied by the caller: the
  periodic point set needs its own, and the wall and flip kernels below them
  are shared.
 */
template <typename T, typename Tint, typename Tgroup, typename Fstab>
IsoDelaunayDomain<T, Tint, Tgroup>
RandomWalkStab(IsoDelaunayDomain<T, Tint, Tgroup> const &x,
               DataIsoDelaunayDomains<T, Tint, Tgroup> &data,
               Fstab f_stab_gens, int n_iter,
               [[maybe_unused]] std::ostream &os) {
  IsoDelaunayDomain<T, Tint, Tgroup> Work = x;
  for (int iter = 0; iter < n_iter; iter++) {
    // Only the combinatorial part is needed to know which walls are
    // flippable; the flip itself is done for the single chosen neighbour.
    PreResultDelaunayAdj<T, Tint, Tgroup> pre =
        get_pre_result_delaunay_adj_kernel<T, Tint, Tgroup>(Work, data,
                                                            f_stab_gens(Work));
    std::vector<size_t> l_flippable;
    for (size_t k = 0; k < pre.l_test.size(); k++) {
      if (pre.l_test[k]) {
        l_flippable.push_back(k);
      }
    }
    size_t n_adj = l_flippable.size();
    if (n_adj == 0) {
#ifdef DEBUG_COVERING_RECORD
      os << "COVERING_RECORD: RandomWalkStab, no adjacent domain at iter="
         << iter << ", stopping early\n";
#endif
      break;
    }
    size_t pos = random() % n_adj;
    Work = get_adjacent(Work, data, pre.ListIneqRed[l_flippable[pos]]).DT_gram;
  }
  return Work;
}

/*
  The walk itself. See the header comment for the strategy; f_stab_gens
  supplies the generators of the stabilizer of the object being walked on
  (the periodic subgroup for a periodic point set), and f_restart a fresh
  starting domain.

  --- The drift, and why f_restart is needed ---

  A flip re-expresses the tessellation in the same lattice basis, and
  nothing brings it back. Composing thousands of them makes the coordinates
  grow without bound: on Z^3 + {0, (1/2,1/2,1/2), (3/4,3/4,3/4)}, whose 46
  domains have Gram matrices of largest entry at most 108, a twenty minute
  walk was reaching representatives of largest entry 1e15, with a median of
  1e12. Those are the same domains seen in ever more skewed coordinates: the
  covering optimum is unchanged, being an invariant, but the forms interior
  to them become so anisotropic that the optimization cannot be carried out
  in double, and 95% of the evaluations of that run were lost that way.

  The enumeration does not suffer from this because it maps every domain it
  reaches to a canonical representative; a walk has no such step. Nor can
  the representative simply be reduced here: a unimodular reduction of the
  form would move the cosets with it, giving a different point set unless it
  is taken in the subgroup preserving the one at hand.

  So the walk restarts instead. Once the current domain is more skewed than
  max_gram_entry, or once its optimization fails, it is dropped for a fresh
  domain from f_restart, whose coordinates are small by construction. The
  best found so far is kept across restarts, so nothing is lost but the
  position.
 */
template <typename T, typename Tint, typename Tgroup, typename Fstab,
          typename Frestart>
RecordSearchResult<Tint>
LookForRecordCovering(DataIsoDelaunayDomains<T, Tint, Tgroup> &data,
                      IsoDelaunayDomain<T, Tint, Tgroup> const &start,
                      LinSpaceMatrix<T> const &LinSpa, T const &point_density,
                      Fstab f_stab_gens, Frestart f_restart,
                      RecordSearchOptions const &opts, std::ostream &os) {
  SingletonTime time_start;
  RecordSearchResult<Tint> res;
  // Below this relative gain an improvement is the noise of the numerical
  // solve rather than a different domain, and is not worth reporting.
  double relative_noise = 1e-9;
  auto out_of_time = [&]() -> bool {
    return opts.max_runtime_second > 0 &&
           si(time_start) > opts.max_runtime_second;
  };
  auto finish = [&](std::string const &message) -> RecordSearchResult<Tint> {
    res.message = message;
    res.runtime_second = si(time_start);
    return res;
  };
  // Evaluating a domain and keeping it if it improves on everything seen so
  // far. The Gram matrix of an improvement is written out immediately: a
  // long walk that is later killed still leaves its best find on disk.
  auto register_domain =
      [&](IsoDelaunayDomain<T, Tint, Tgroup> const &y,
          covering_maxdet::MaxdetResult<double> const &opt) -> void {
    if (res.has_best && opt.cov_density >= res.best_density) {
      return;
    }
    // The walk returns to the same domain over and over, and the optimizer
    // lands on a slightly different point every time, so a bare "strictly
    // better" would report and write a file on every revisit. Only an
    // improvement larger than the noise of the solve is worth saying, the
    // best is tracked either way, and a record is always worth saying.
    bool is_record = opt.cov_density < opts.record;
    bool is_real_gain =
        !res.has_best ||
        opt.cov_density < res.best_density * (1 - relative_noise);
    res.has_best = true;
    res.best_density = opt.cov_density;
    res.best_domain_gram = y.GramMat;
    res.best_gram = opt.Q;
    res.best_covering_radius_sq = opt.cov_radius_sq;
    if (is_record) {
      res.found_record = true;
    }
    if (!is_real_gain && !is_record) {
      return;
    }
    os << "COVERING_RECORD: new best density=" << opt.cov_density
       << " record=" << opts.record << "\n";
    if (opts.Prefix.size() > 0) {
      std::string FileOut = FILE_FindAvailableFileFromPrefix(opts.Prefix);
      WriteMatrixFile(FileOut, y.GramMat);
      os << "COVERING_RECORD: wrote the domain Gram matrix to " << FileOut
         << "\n";
    }
  };
  auto evaluate = [&](IsoDelaunayDomain<T, Tint, Tgroup> const &y)
      -> std::optional<covering_maxdet::MaxdetResult<double>> {
    std::optional<covering_maxdet::MaxdetResult<double>> opt =
        GetDomainCoveringOptimum(y, LinSpa, point_density, os);
    if (opt) {
      res.n_domain_evaluated++;
    } else {
      res.n_domain_failed++;
    }
    return opt;
  };
  //
  IsoDelaunayDomain<T, Tint, Tgroup> Work = start;
  double curr = 0;
  // Move to a fresh domain and take its value, the escape from a drifted or
  // unusable position. Returns false when even a fresh domain cannot be
  // optimized, which stops the walk.
  auto restart = [&]() -> bool {
    res.n_restart++;
    for (int i_try = 0; i_try < 10; i_try++) {
      Work = f_restart();
      std::optional<covering_maxdet::MaxdetResult<double>> opt =
          evaluate(Work);
      if (opt) {
        register_domain(Work, *opt);
        curr = opt->cov_density;
        return true;
      }
    }
    return false;
  };
  auto is_drifted = [&](IsoDelaunayDomain<T, Tint, Tgroup> const &y) -> bool {
    return opts.max_gram_entry > 0 && GetDomainSkew(y) > opts.max_gram_entry;
  };
  std::optional<covering_maxdet::MaxdetResult<double>> opt_curr =
      evaluate(Work);
  if (!opt_curr) {
    if (!restart()) {
      return finish("no domain could be optimized, even after restarting");
    }
  } else {
    register_domain(Work, *opt_curr);
    curr = opt_curr->cov_density;
  }
  if (res.found_record) {
    return finish("the starting domain already beats the record");
  }
  while (true) {
    if (out_of_time()) {
      return finish("the runtime budget ran out");
    }
    res.n_iter++;
    ResultDelaunayAdj<T, Tint, Tgroup> result =
        get_result_delaunay_adj_kernel<T, Tint, Tgroup>(Work, data,
                                                        f_stab_gens(Work));
    size_t n_adj = result.l_adj.size();
    if (n_adj == 0) {
      // The parameter space is a single domain, so the walk has nowhere to
      // go and the answer is already the global one.
      return finish("the domain has no flippable wall, the walk is over");
    }
    std::vector<size_t> ListIdx;
    double the_min = 0;
    bool has_min = false;
    for (size_t i = 0; i < n_adj; i++) {
      if (out_of_time()) {
        break;
      }
      std::optional<covering_maxdet::MaxdetResult<double>> opt =
          evaluate(result.l_adj[i].DT_gram);
      if (!opt) {
        continue;
      }
      register_domain(result.l_adj[i].DT_gram, *opt);
      if (res.found_record) {
        return finish("a domain beating the record was found");
      }
      if (!has_min || opt->cov_density < the_min) {
        the_min = opt->cov_density;
        has_min = true;
        ListIdx.clear();
      }
      if (opt->cov_density == the_min) {
        ListIdx.push_back(i);
      }
    }
#ifdef DEBUG_COVERING_RECORD
    os << "COVERING_RECORD: n_iter=" << res.n_iter << " n_adj=" << n_adj
       << " curr=" << curr << " the_min=" << (has_min ? the_min : -1) << "\n";
#endif
    if (out_of_time()) {
      // The neighbourhood was only partly evaluated, so neither the descent
      // nor the escape move would be based on the whole of it. Stopping
      // here also keeps the budget honest, an escape move being
      // n_walk_steps flips long.
      return finish("the runtime budget ran out");
    }
    if (!has_min || the_min >= curr) {
      // A local minimum, or a neighbourhood none of whose domains could be
      // optimized: jump out of the basin.
      Work = RandomWalkStab<T, Tint, Tgroup, Fstab>(Work, data, f_stab_gens,
                                                    opts.n_walk_steps, os);
      res.n_walk++;
      if (is_drifted(Work)) {
        if (!restart()) {
          return finish("no domain could be optimized, even after restarting");
        }
      } else {
        std::optional<covering_maxdet::MaxdetResult<double>> opt =
            evaluate(Work);
        if (!opt) {
          // Not usable, and staying would only walk deeper into the same
          // region, so start again from a fresh domain.
          if (!restart()) {
            return finish(
                "no domain could be optimized, even after restarting");
          }
        } else {
          register_domain(Work, *opt);
          curr = opt->cov_density;
        }
      }
      if (res.found_record) {
        return finish("a domain beating the record was found");
      }
    } else {
      size_t pos = random() % ListIdx.size();
      Work = result.l_adj[ListIdx[pos]].DT_gram;
      curr = the_min;
      if (is_drifted(Work) && !restart()) {
        return finish("no domain could be optimized, even after restarting");
      }
      if (res.found_record) {
        return finish("a domain beating the record was found");
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------------

/*
  The outcome as a GAP record. Written whether or not a record was found:
  not finding one is the expected result, and the best density seen is the
  useful part of it.
 */
template <typename Tint>
void WriteRecordSearchGAP(std::string const &FileName,
                          RecordSearchResult<Tint> const &res,
                          double const &record) {
  std::ofstream os_out(FileName);
  os_out << std::setprecision(17) << std::showpoint;
  os_out << "return rec(found_record:="
         << (res.found_record ? "true" : "false");
  os_out << ", has_best:=" << (res.has_best ? "true" : "false");
  os_out << ", record:=" << record;
  os_out << ", best_density:=" << res.best_density;
  os_out << ", best_covering_radius_sq:=" << res.best_covering_radius_sq;
  os_out << ", message:=\"" << res.message << "\"";
  os_out << ", n_iter:=" << res.n_iter;
  os_out << ", n_walk:=" << res.n_walk;
  os_out << ", n_restart:=" << res.n_restart;
  os_out << ", n_domain_evaluated:=" << res.n_domain_evaluated;
  os_out << ", n_domain_failed:=" << res.n_domain_failed;
  os_out << ", runtime_second:=" << res.runtime_second;
  if (res.has_best) {
    os_out << ", best_domain_gram:=" << StringMatrixGAP(res.best_domain_gram);
    os_out << ", best_gram:=[";
    for (int i = 0; i < res.best_gram.rows(); i++) {
      if (i > 0) {
        os_out << ", ";
      }
      os_out << "[";
      for (int j = 0; j < res.best_gram.cols(); j++) {
        if (j > 0) {
          os_out << ", ";
        }
        os_out << res.best_gram(i, j);
      }
      os_out << "]";
    }
    os_out << "]";
  }
  os_out << ");\n";
}

// clang-format off
}  // namespace covering_record
// clang-format on

// clang-format off
#endif  // SRC_DELAUNAY_COVERINGRECORDSEARCH_H_
// clang-format on

// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Thompson-sampled variants of DirectDualDescription_{mat,vf,f}.
//
// The plain DirectDualDescription_* functions in POLY_DirectDualDesc.h
// pick their method through get_dual_desc_method, which currently just
// returns "lrs" unconditionally. That leaves on the table the choice
// between "lrs", "cdd" and "bb" (beneath-and-beyond), each of which can
// be genuinely faster than the others depending on the shape of EXT. The
// _ts variants here let a caller supply a ThompsonSamplingHeuristic that
// samples between "lrs", "cdd" and "bb", records the wall time of each
// attempt, and progressively converges on the empirically-fastest choice
// for the caller's workload.
//
// Contract on the sampler:
//   * The Thompson state consumed by the _ts functions must list the
//     methods "lrs", "cdd" and "bb" as options in its ListDescription
//     (e.g. "lrs:distri1 cdd:distri1 bb:distri1"), so that all three are
//     Thompson-sampled at every level. Any answer other than
//     "lrs"/"cdd"/"bb" is defensively mapped to "lrs" (which works for
//     both fields and rings). Note that "cdd" and "bb" require a field.
//   * The caller owns the sampler's lifetime. Its posterior distribution
//     is stateful; reuse the same instance across calls to accumulate
//     learning.
//
// The feature map passed to get_eval() uses the standard PolytopeField
// key names, matching what the recursive-dual-desc layer already does
// via PolytopeInputInfo::named_map(). At this free-standing layer we
// have no group or recursion context, so we set groupsize=1 and
// delta/level/time=0. The sampler's key-compression configuration
// decides which features actually condition the choice; unused fields
// are silently ignored.

#ifndef SRC_DUALDESC_POLY_DIRECTDUALDESC_TS_H_
#define SRC_DUALDESC_POLY_DIRECTDUALDESC_TS_H_

// clang-format off
#include "POLY_DirectDualDesc.h"
#include "POLY_Heuristics.h"
#include <map>
#include <string>
// clang-format on

namespace poly_direct_dual_desc_ts_detail {

// Build the feature map the Thompson sampler consumes. We reuse
// PolytopeInputInfo::named_map() so the set of keys stays in sync
// automatically if PolytopeField ever grows a member.
//
// The workload-discriminating feature is `delta = incidence - rank`
// (nbRow - nbCol), the same quantity ComputeInitialInfo uses in the
// recursive-dual-desc layer. It measures how many vertices/rays the
// polytope has beyond the ambient dimension — small delta favors one
// method, large delta the other. Fields we can't populate at this
// layer (groupsize / level / time) get sentinel defaults; the
// sampler's key-compression config decides which fields actually
// condition the choice.
template <typename Tint>
inline std::map<std::string, Tint> ts_features(int nbRow, int nbCol) {
  PolytopeInputInfo<Tint> info;
  info.groupsize = Tint(1);
  info.incidence = nbRow;
  info.rank = nbCol;
  info.delta = nbRow - nbCol;
  info.level = 0;
  info.time = 0;
  return info.named_map();
}

// Sample a program and map the string to a DualDescProgram. Only "lrs",
// "cdd" and "bb" are exposed; anything else is a namelist-configuration bug
// and we fall back to "lrs" (which is a safe default because it works for
// both fields and rings).
template <typename T, typename Tint>
inline DualDescProgram
ts_pick_program(MyMatrix<T> const &EXT,
                ThompsonSamplingHeuristic<Tint> &eTS) {
  std::string choice = eTS.get_eval(ts_features<Tint>(EXT.rows(), EXT.cols()));
  if (choice == "cdd")
    return DualDescProgram::cdd;
  if (choice == "lrs")
    return DualDescProgram::lrs;
  if (choice == "bb")
    return DualDescProgram::beneath_beyond;
  return DualDescProgram::lrs;
}

}  // namespace poly_direct_dual_desc_ts_detail

// Zero-configuration factory. Returns a ThompsonSamplingHeuristic<Tint>
// preconfigured for DirectDualDescription_*_ts:
//
//   * A single Thompson state ("state_opts") offers "lrs", "cdd" and "bb"
//     as options, so all three are Thompson-sampled at every level.
//   * Only one feature is used to separate the sampler's blocks: `delta`
//     (= nbRow - nbCol, the polytope's excess-vertices count).
//   * `delta` is binned in intervals of width 5: [0,4], [5,9], …,
//     [95,99], then a tail bin [100, infinity) so very wide polytopes
//     don't create an unbounded number of empty bins.
//   * The empirical prior is dirac at 0.0 with a weak seed (Nstart = 2)
//     and modest memory (Nmax = 25). All answers start indistinguishable;
//     the sampler explores each and lets timings differentiate them.
//   * No classical-heuristic override — the sampler drives every choice.
//
// The returned sampler is stateful: **reuse the same instance across calls**
// to accumulate learning. Constructing a fresh sampler per call throws away
// everything the previous ones learned.
template <typename Tint>
inline ThompsonSamplingHeuristic<Tint>
MakeDualDescThompsonSampler(std::ostream &os) {
  // Build the KEY_COMPRESSION.ListDescription string: 20 finite width-5
  // bins covering [0, 99], plus one open-ended tail bin [100, infinity).
  std::string bin_desc;
  int const step = 5;
  int const n_finite_bins = 20;
  for (int i = 0; i < n_finite_bins; i++) {
    int lo = i * step;
    int hi = lo + step - 1;
    if (!bin_desc.empty())
      bin_desc += ",";
    bin_desc += std::to_string(lo) + "-" + std::to_string(hi);
  }
  bin_desc += "," + std::to_string(n_finite_bins * step) + "-infinity";

  std::vector<std::string> lstr_proba = {
      "&PROBABILITY_DISTRIBUTIONS",
      " ListName = \"distri1\"",
      " ListNmax = 25",
      " ListNstart = 2",
      " ListNature = \"dirac\"",
      " ListDescription = \"0.0\"",
      "/",
  };
  // A single Thompson state offering *all three* methods. Because "lrs",
  // "cdd" and "bb" appear in the same ListDescription, get_lowest_sampling
  // draws from each option's empirical distribution and returns the
  // fastest — i.e. the three options are genuinely sampled against each
  // other at every level. (ListAnswer is a per-state label; it only
  // matters for noprior states, which we do not use here, but its length
  // must match ListName / ListDescription.)
  std::vector<std::string> lstr_thompson_prior = {
      "&THOMPSON_PRIOR",
      " ListAnswer = \"lrs_cdd_bb\"",
      " ListName = \"state_opts\"",
      " ListDescription = \"lrs:distri1 cdd:distri1 bb:distri1\"",
      "/",
  };
  std::vector<std::string> lstr_key = {
      "&KEY_COMPRESSION",
      " ListKey = \"delta\"",
      " ListDescription = \"" + bin_desc + "\"",
      "/",
  };
  // Every delta bucket resolves to the single "state_opts" state, which
  // offers "lrs", "cdd" and "bb". delta = nbRow - nbCol is >= 0 for any
  // well-posed dual-desc input, so the "delta >= 0" condition is always
  // true; together with the identical DefaultPrior it guarantees that all
  // three options are Thompson-sampled at every level, whatever delta is.
  std::vector<std::string> lstr_heuristic_prior = {
      "&HEURISTIC_PRIOR",
      " DefaultPrior = \"state_opts\"",
      " ListFullCond = \"delta >= 0\"",
      " ListConclusion = \"state_opts\"",
      "/",
  };
  std::vector<std::string> lstr_io = {
      "&IO",
      " name = \"dualdesc_ts\"",
      " WriteLog = F",
      " ProcessExistingDataIfExist = F",
      " LogFileToProcess = \"input_logfile\"",
      "/",
  };

  std::vector<std::string> lstr;
  for (auto const &e_lstr : {lstr_proba, lstr_thompson_prior, lstr_key,
                             lstr_heuristic_prior, lstr_io}) {
    lstr.push_back("");
    for (auto &estr : e_lstr)
      lstr.push_back(estr);
  }
  FullNamelist eFull = NAMELIST_ThompsonSamplingRuntime();
  NAMELIST_ReadListString(eFull, lstr);
  return ThompsonSamplingHeuristic<Tint>(eFull, os);
}

// Returns the facet inequalities as an inequality matrix.
template <typename T, typename Tint>
MyMatrix<T>
DirectDualDescription_mat_ts(MyMatrix<T> const &EXT,
                             ThompsonSamplingHeuristic<Tint> &eTS,
                             std::ostream &os) {
  DualDescProgram prog =
      poly_direct_dual_desc_ts_detail::ts_pick_program(EXT, eTS);
  MyMatrix<T> ret = DirectFacetComputationInequalities(EXT, prog, os);
  eTS.pop(os);
  return ret;
}

// Returns the facet incidences (per-facet Face bitmask) as a vectface.
template <typename T, typename Tint>
vectface
DirectDualDescription_vf_ts(MyMatrix<T> const &EXT,
                            ThompsonSamplingHeuristic<Tint> &eTS,
                            std::ostream &os) {
  DualDescProgram prog =
      poly_direct_dual_desc_ts_detail::ts_pick_program(EXT, eTS);
  vectface ret = DirectFacetComputationIncidence(EXT, prog, os);
  eTS.pop(os);
  return ret;
}

// Streams each facet as a (Face, MyVector<T>) pair to `f_process`.
template <typename T, typename Tint, typename Fprocess>
void
DirectDualDescription_f_ts(MyMatrix<T> const &EXT,
                           ThompsonSamplingHeuristic<Tint> &eTS,
                           Fprocess f_process, std::ostream &os) {
  DualDescProgram prog =
      poly_direct_dual_desc_ts_detail::ts_pick_program(EXT, eTS);
  DirectFacetComputationFaceIneq(EXT, prog, f_process, os);
  eTS.pop(os);
}

// clang-format off
#endif  // SRC_DUALDESC_POLY_DIRECTDUALDESC_TS_H_
// clang-format on

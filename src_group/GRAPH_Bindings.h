// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_GRAPH_BINDINGS_H_
#define SRC_GROUP_GRAPH_BINDINGS_H_

// This is bindings to the bliss or traces code.
//
// The default is TRACES. If not selected then
// the BLISS code is used.

#ifndef USE_BLISS
#define USE_TRACES
#endif

#ifdef USE_BLISS
#include "GRAPH_bliss.h"
#endif

#ifdef USE_TRACES
#include "GRAPH_traces.h"
#endif

template<typename Tgr, typename TidxC>
std::vector<TidxC> GRAPH_GetCanonicalOrdering(Tgr const& eGR, [[maybe_unused]] std::ostream& os) {
#ifdef USE_BLISS
  return BLISS_GetCanonicalOrdering<Tgr, TidxC>(eGR);
#endif
#ifdef USE_TRACES
  return TRACES_GetCanonicalOrdering<Tgr, TidxC>(eGR, os);
#endif
}

template<typename Tgr, typename TidxC, typename TidxG>
std::pair<std::vector<TidxC>, std::vector<std::vector<TidxG>>>
GRAPH_GetCanonicalOrdering_ListGenerators(Tgr const& eGR, size_t const& nbRow, [[maybe_unused]] std::ostream& os) {
#ifdef USE_BLISS
  return BLISS_GetCanonicalOrdering_ListGenerators<Tgr, TidxC, TidxG>(eGR, nbRow);
#endif
#ifdef USE_TRACES
  return TRACES_GetCanonicalOrdering_ListGenerators<Tgr, TidxC, TidxG>(eGR, nbRow, os);
#endif
}

template<typename Tgr, typename TidxG>
std::vector<std::vector<TidxG>> GRAPH_GetListGenerators(Tgr const& eGR, size_t const& nbRow, [[maybe_unused]] std::ostream& os) {
#ifdef USE_BLISS
  return BLISS_GetListGenerators<Tgr, TidxG>(eGR, nbRow);
#endif
#ifdef USE_TRACES
  return TRACES_GetListGenerators<Tgr, TidxG>(eGR, nbRow, os);
#endif
}



// clang-format off
#endif  // SRC_GROUP_GRAPH_BINDINGS_H_
// clang-format on


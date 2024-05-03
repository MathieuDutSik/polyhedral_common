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

struct SimplifiedVectexColoredGraph {
  size_t nbVert;
  std::vector<int> d;
  std::vector<size_t> v;
  std::vector<int> e;
  std::vector<int> ListBlockSize;
};


#ifdef USE_BLISS
GraphListAdj GetGraphListAdj_from_simplified(SimplifiedVectexColoredGraph const& x) {
  size_t nbVert = d.nbVert;
  GraphListAdj eGR(nbVert);
  size_t pos = 0;
  size_t i_color = 0;
  for (auto & blk_size : x.ListBlockSize) {
    size_t blk_size_z = blk_size;
    for (size_t u=0; u<blk_size_z; u++) {
      eGR.SetColor(pos, i_color);
      pos++;
    }
    i_color++;
  }
  int pos = 0;
  for (size_t iVert=0; iVert<nbVert; iVert++) {
    int nbAdj = x.d[iVert];
    for (int u=0; u<nbAdj; u++) {
      size_t jVert = x.e[pos];
      eGR.AddAdjacent(iVert, jVert);
      pos++;
    }
  }
  return eGR;
}
#endif

#ifdef USE_TRACES
DataTraces GetDataTraces(SimplifiedVectexColoredGraph const& x) {
  size_t nbVert = x.nbVert;
  int nbAdjacent = 0;
  for (size_t i=0; i<nbVert; i++) {
    nbAdjacent += x.d[i];
  }
  DataTraces DT(nbVert, nbAdjacent);
  int pos = 0;
  for (size_t i=0; i<nbVert; i++) {
    int nbAdj = x.d[i];
    DT.sg1.d[i] = nbAdj;
    DT.sg1.v[i] = pos;
    pos += nbAdj;
  }
  for (size_t i = 0; i < nbVert; i++) {
    DT.lab1[i] = i;
    DT.ptn[i] = NAUTY_INFINITY;
  }
  int tot_size = 0;
  for (auto & blk_size : x.ListBlockSize) {
    tot_size += blk_size;
    DT.ptn[tot_size-1] = 0;
  }
  for (int i_ent=0; i_ent<nbAdjacent; i_ent++) {
    DT.sg1.e[i_ent] = x.e[i_ent];
  }
  return DT;
}
#endif


template<typename TidxC, typename TidxG>
std::pair<std::vector<TidxC>, std::vector<std::vector<TidxG>>>
GRAPH_GetCanonicalOrdering_ListGenerators_Simp(SimplifiedVectexColoredGraph const& s, size_t const& nbRow, std::ostream& os) {
#ifdef USE_BLISS
  GraphListAdj eGR = GetGraphListAdj_from_simplified(s);
  return BLISS_GetCanonicalOrdering_ListGenerators<GraphListAdj, TidxC, TidxG>(eGR, nbRow);
#endif
#ifdef USE_TRACES
  DataTraces DT = GetDataTraces(s);
  return TRACES_GetCanonicalOrdering_ListGenerators_Arr<TidxC, TidxG>(DT, nbRow, os);
#endif
}




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


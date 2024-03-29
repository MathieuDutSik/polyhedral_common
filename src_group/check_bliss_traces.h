// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_CHECK_BLISS_TRACES_H_
#define SRC_GROUP_CHECK_BLISS_TRACES_H_

// clang-format off
#include "GRAPH_bliss.h"
#include "GRAPH_traces.h"
#include <vector>
// clang-format on

template <typename Tgr, typename Tgroup>
void PrintStabilizerGroupSizes(std::ostream &os, Tgr const &eGR) {
  bliss::Graph g = GetBlissGraphFromGraph(eGR);
  size_t nbVert = eGR.GetNbVert();
  std::vector<std::vector<unsigned int>> ListGen1 =
      BLISS_GetListGenerators(eGR);
  std::vector<std::vector<unsigned int>> ListGen2 =
      TRACES_GetListGenerators(eGR);
  auto siz1 = GetGroupListGen<Tgroup>(ListGen1, nbVert).size();
  auto siz2 = GetGroupListGen<Tgroup>(ListGen2, nbVert).size();
  bool test1 = CheckListGenerators(ListGen1, eGR);
  bool test2 = CheckListGenerators(ListGen2, eGR);
  os << "|GRP bliss|=" << siz1 << " |GRP traces|=" << siz2 << " test1=" << test1
     << " test2=" << test2 << "\n";
}

// clang-format off
#endif  // SRC_GROUP_CHECK_BLISS_TRACES_H_
// clang-format on

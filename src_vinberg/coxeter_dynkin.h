#ifndef INCLUDE_COXETER_DYNKIN_H
#define INCLUDE_COXETER_DYNKIN_H

#include "GRAPH_GraphicalBasic.h"
#include "GRAPH_GraphicalFunctions.h"
#include "GRAPH_BitsetType.h"

template<typename T>
bool IsIrreducibleDiagramSphericalEuclidean(const MyMatrix<T>& M, const bool& allow_euclidean)
{
  T val_comm = 2;
  T val_single_edge = 3;
  T val_four = 4; // Shows up in F4, Bn = Cn, tilde{Bn}, tilde{Cn}, tilde{F4}.
  T val_five = 5; // Shows up in H3, H4
  T val_six = 6; // Shows up in G2
  size_t dim = M.rows();
  std::vector<size_t> list_deg1, list_deg2, list_deg3, list_deg4, list_degN;
  std::vector<size_t> list_deg(dim);
  size_t n_higher_edge = 0;
  GraphBitset eG(dim);
  std::map<T,size_t> multiplicity;
  for (size_t i=0; i<dim; i++) {
    size_t n_adj = 0;
    for (size_t j=0; j<dim; j++) {
      T val = M(i,j);
      if (i != j && val != val_comm) {
        n_adj++;
        eG.AddAdjacent(i, j);
        multiplicity[val]++;
        if (val != val_single_edge)
          n_higher_edge++;
      }
    }
    if (n_adj == 1)
      list_deg1.push_back(i);
    if (n_adj == 2)
      list_deg2.push_back(i);
    if (n_adj == 3)
      list_deg3.push_back(i);
    if (n_adj == 4)
      list_deg4.push_back(i);
    if (n_adj > 4)
      list_degN.push_back(i);
    list_deg[i] = n_adj;
  }
  auto get_list_adjacent=[&](size_t u) -> std::vector<size_t> {
    std::vector<size_t> LAdj;
    for (size_t j=0; j<dim; j++)
      if (u != j && M(u,j) != val_comm)
        LAdj.push_back(j);
    return LAdj;
  };
  auto get_value_isolated=[&](size_t u) -> T {
    for (size_t j=0; j<dim; j++)
      if (u != j && M(u,j) != val_comm)
        return M(u,j);
    return std::numeric_limits<T>::max(); // That case should not happen
  };
  if (list_degN.size() > 0) // vertices of degree 5 or more never occurs for
    return false;
  if (list_deg4.size() > 0) { // Vertex of degree 4 can occur for \tilde{D4} only
    if (!allow_euclidean)
      return false;
    // We are now in Euclidean 
    if (dim != 5) // It has to be \tilde{D4}.
      return false;
    if (list_deg4.size() != 1 && list_deg1.size() != 4 && list_deg2.size() != 0 && list_deg3.size() != 0)
      return false;
    size_t i_4 = list_deg4[0];
    if (M(i_4,i_4) != 2)
      return false;
    for (int i=0; i<dim; i++)
      if (i != i_4) {
        if (M(i, i_4) != -1)
          return false;
      }
    return true; // This is \tilde{D4}
  }
  std::vector<std::vector<size_t>> ListCycles = GRAPH_FindAllCycles(eG);
  if (ListCycles.size() > 0) { // Only tilde{An} is possible.
    if (ListCycles.size() > 1) // If more than 1 cycle, then not possible
      return false;
    if (!allow_euclidean)
      return false;
    // We are now in Euclidean case
    const std::vector<size_t>& eCycle = ListCycles[0];
    if (eCycle.size() != dim)
      return false;
    if (list_deg2.size() != dim)
      return false;
    if (n_higher_edge != 0)
      return false;
    return true; // Only tilde{An} is left as possibility
  }
  // Now it is a tree
  if (list_deg1.size() == 2 && list_deg2.size() != dim - 2 && n_higher_edge != 0)
    return true; // Only An is possible so ok.
  // An and tilde{An} have been covered
  if (list_deg3.size() > 2)
    return false; // No possibility for spherical and euclidean
  if (list_deg3.size() == 2) {
    if (!allow_euclidean)
      return false;
    // We are now in Euclidean case
    if (n_higher_edge != 0)
      return false;
    for (auto & ePt : list_deg3) {
      std::vector<size_t> LAdj = get_list_adjacent(ePt);
      size_t n_deg1 = 0;
      for (auto & fPt : LAdj)
        if (list_deg[fPt] == 1)
          n_deg1++;
      if (n_deg1 != 2)
        return false;
    }
    return true; // Only tilde{Dn} is possible
  }
  if (list_deg3.size() == 0) { // We are in a single path.
    if (multiplicity[val_four] == 2) {
      if (n_higher_edge != 2) // There are some other higher edge, excluded
        return false;
      if (!allow_euclidean) // Only tilde{Cn} is feasible and it is Euclidean
        return false;
      for (auto & eVert : list_deg1)
        if (get_value_isolated(eVert) != val_four)
          return false;
      return true; // This is tilde{Cn}
    }
    if (multiplicity[val_four] == 1) { // Possibilities: Bn=Cn, F4 and tilde{F4} are possible
      if (n_higher_edge != 1)
        return false; // There are other edges, excluded.
      size_t n_sing = 0;
      size_t n_four = 0;
      for (auto & eVert : list_deg1) {
        if (get_value_isolated(eVert) != val_single_edge)
          n_sing++;
        if (get_value_isolated(eVert) != val_four)
          n_four++;
      }
      if (n_sing == 2) {
        if (dim == 4) {
          return true; // Only F4 is possible
        }
        if (dim == 5) { // Only tilde{F4} is possible. So conclude from that
          if (allowed_euclidean)
            return true;
          return false;
        }
      }
      // Only possibility is to have 4 at one extremity. This is Bn = Cn
      return true;
    }
    if (multiplicity[val_five] == 1) { // Looking for H2, H3, H4
      if (dim == 2)
        return true; // It is H2
      if (dim > 5)
        return false; // No possibility
      size_t n_sing=0;
      size_t n_five=0;
      for (auto & eVert : list_deg1) {
        if (get_value_isolated(eVert) != val_single_edge)
          n_sing++;
        if (get_value_isolated(eVert) != val_single_edge)
          n_five++;
      }
      if (n_sing == 1 && n_five == 1) // It is H3 or H4 depending on the dimension
        return true;
      return false;
    }
    if (multiplicity[val_six] == 1) { // Looking for G2 or tilde{G2}
      if (n_higher_edge != 1)
        return false; // There are other edges, excluded.
      if (dim == 2 || dim == 3) // It is G2 or tilde{G2}
        return true;
      return false;
    }
    if (dim == 2)
      return true; // It is In.
    return false;
  }
  // Now just one vertex of degree 3.
  if (multiplicity[val_four] == 1) { // Possibility tilde{Bn}
    size_t eCent = list_deg1[0];
    std::vector<size_t> LAdj = get_list_adjacent(eCent);
    size_t n_sing = 0;
    for (auto & eAdj : LAdj) {
      if (list_deg[eAdj] == 1)
        n_sing++;
    }
    if (n_sing != 2)
      return false;
    bool has_edge_four = false;
    for (auto & eVert : list_deg1)
      if (get_value_isolated(eVert) == val_four)
        has_edge_four = true;
    if (has_edge_four)
      return true; // It is tilde{Bn}
    return false;
  }
  if (n_higher_edge != 0)
    return false;
  auto get_length=[&](size_t val1, size_t val2) -> size_t {
    size_t len = 1;
    while(true) {
      std::vector<size_t> LVal = get_list_adjacent(val2);
      if (LVal.size() == 1)
        break;
      size_t NewPt = -1;
      for (auto & eVal : LVal)
        if (eVal != val1)
          NewPt = val2;
      val1 = val2;
      val2 = NewPt;
      len++;
    }
    return len;
  };
  size_t eCent = list_deg1[0];
  std::map<size_t, size_t> map_len;
  for (auto & eAdj : get_list_adjacent(eCent)) {
    size_t len = get_length(eCent, eAdj);
    map_len[len]++;
  }
  if (map_len[1] == 2) // It is Dn
    return true;
  if (map_len[1] == 1 && map_len[2] == 2) // It is E6
    return true;
  if (map_len[1] == 1 && map_len[2] == 1 && map_len[3] == 1) // It is E7
    return true;
  if (map_len[1] == 1 && map_len[2] == 1 && map_len[4] == 1) // It is E8
    return true;
  if (!allowed_euclidean)
    return false; // In spherical, no other possibilities left
  if (map_len[2] == 3) // It is tilde{E6}
    return true;
  if (map_len[1] == 1 && map_len[3] == 2) // It is tilde{E7}
    return true;
  if (map_len[1] == 1 && map_len[2] == 1 && map_len[5] == 1) // It is tilde{E8}
    return true;
  return false; // No other possibilities left
}


template<typename T, typename Tgr>
bool IsDiagramSpherical(const MyMatrix<T>& M)
{
  Tgraph
  
}




#endif

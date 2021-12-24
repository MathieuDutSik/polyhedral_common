#ifndef INCLUDE_COXETER_DYNKIN_H
#define INCLUDE_COXETER_DYNKIN_H

#include "GRAPH_GraphicalBasic.h"
#include "GRAPH_GraphicalFunctions.h"
#include "GRAPH_BitsetType.h"


/*
  The Coxeter Dynkin diagram are build in the following way:
  --- The values for i != j are the exponent m such that (g_ig_j)^m = Id
  --- The values M(i,i) are not used.
 */


template<typename T>
bool IsIrreducibleDiagramSphericalEuclidean(const MyMatrix<T>& M, const bool& allow_euclidean)
{
  T val_comm = 2;
  T val_single_edge = 3;
  T val_four = 4; // Shows up in F4, Bn = Cn, tilde{Bn}, tilde{Cn}, tilde{F4}.
  T val_five = 5; // Shows up in H3, H4
  T val_six = 6; // Shows up in G2
  size_t dim = M.rows();
  if (dim == 1) // Case of A1
    return true;
  std::vector<size_t> list_deg1, list_deg2, list_deg3, list_deg4, list_degN;
  std::vector<size_t> list_deg(dim, 0);
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
  if (list_degN.size() > 0) // vertices of degree 5 or more never occurs.
    return false;
  if (list_deg4.size() > 0) { // Vertex of degree 4 can occur for \tilde{D4} only
    if (!allow_euclidean) // Only possibilities is not allowed, exit.
      return false;
    // We are now in Euclidean 
    if (dim != 5) // It has to be \tilde{D4}.
      return false;
    if (list_deg4.size() != 1 || list_deg1.size() != 4 || list_deg2.size() != 0 || list_deg3.size() != 0)
      return false;
    size_t i_4 = list_deg4[0];
    for (size_t i=0; i<dim; i++)
      if (i != i_4)
        if (M(i, i_4) != 3)
          return false;
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
  if (list_deg1.size() == 2 && list_deg2.size() == dim - 2 && n_higher_edge == 0)
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
          if (allow_euclidean)
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
  if (!allow_euclidean)
    return false; // In spherical, no other possibilities left
  if (map_len[2] == 3) // It is tilde{E6}
    return true;
  if (map_len[1] == 1 && map_len[3] == 2) // It is tilde{E7}
    return true;
  if (map_len[1] == 1 && map_len[2] == 1 && map_len[5] == 1) // It is tilde{E8}
    return true;
  return false; // No other possibilities left
}


template<typename T>
bool IsDiagramSpherical(const MyMatrix<T>& M, const bool& allow_euclidean)
{
  T val_comm = 2;
  size_t dim = M.rows();
  std::cerr << "IsDiagramSpherical dim=" << dim << "\n";
  GraphBitset eG(dim);
  for (size_t i=0; i<dim; i++) {
    for (size_t j=i+1; j<dim; j++) {
      if (M(i,j) != val_comm) {
        eG.AddAdjacent(i,j);
        eG.AddAdjacent(j,i);
      }
    }
  }
  std::cerr << "eG is built\n";
  std::vector<std::vector<size_t>> LConn = ConnectedComponents_set(eG);
  std::cerr << "LConn is built\n";
  for (auto & eConn : LConn) {
    size_t dim_res=eConn.size();
    std::cerr << "dim_res=" << dim_res << "\n";
    MyMatrix<T> Mres(dim_res, dim_res);
    for (size_t i=0; i<dim_res; i++)
      for (size_t j=0; j<dim_res; j++)
        Mres(i,j) = M(eConn[i], eConn[j]);
    std::cerr << "Before IsIrreducibleDiagramSphericalEuclidean\n";
    bool test = IsIrreducibleDiagramSphericalEuclidean(M, allow_euclidean);
    std::cerr << "After IsIrreducibleDiagramSphericalEuclidean\n";
    if (!test)
      return false;
  }
  return true;
}


template<typename T>
std::vector<MyVector<T>> FindDiagramExtensions(const MyMatrix<T>& M, const bool& allow_euclidean)
{
  std::cerr << "FindDiagramExtensions, step 1\n";
  std::set<MyVector<T>> SetExtensions;
  T val_comm = 2;
  T val_single_edge = 3;
  T val_six = 6;
  size_t dim = M.rows();
  std::vector<size_t> list_deg(dim);
  std::vector<size_t> list_isolated;
  for (size_t i=0; i<dim; i++) {
    size_t n_adj = 0;
    for (size_t j=0; j<dim; j++) {
      T val = M(i,j);
      if (i != j && val != val_comm)
        n_adj++;
    }
    list_deg[i] = n_adj;
    if (n_adj == 0)
      list_isolated.push_back(i);
  }
  // Consider the case of adding unconnected vector
  std::cerr << "FindDiagramExtensions, step 2\n";
  MyVector<T> V_basic(dim);
  for (size_t i=0; i<dim; i++)
    V_basic(i) = val_comm;
  std::cerr << "FindDiagramExtensions, step 3\n";
  auto test_vector_and_insert=[&](const MyVector<T>& V) -> void {
    MyMatrix<T> Mtest(dim+1,dim+1);
    for (size_t i=0; i<dim; i++)
      for (size_t j=0; j<dim; j++)
        Mtest(i,j) = M(i,j);
    for (size_t i=0; i<dim; i++) {
      Mtest(i,dim) = V(i);
      Mtest(dim,i) = V(i);
    }
    std::cerr << "Mtest built\n";
    if (IsDiagramSpherical(Mtest, allow_euclidean))
      SetExtensions.insert(V);
  };
  test_vector_and_insert(V_basic);
  std::cerr << "FindDiagramExtensions, step 4\n";
  // Considering the case of just one edge
  for (size_t i=0; i<dim; i++) {
    // Here we have an arbitrary value
    for (T val=val_single_edge; val<128; val++) {
      std::cerr << "i=" << i << " val=" << val << "\n";
      MyVector<T> V = V_basic;
      V(i) = val;
      test_vector_and_insert(V);
    }
  }
  std::cerr << "FindDiagramExtensions, step 5\n";
  // Considering the case of 2 edges
  for (size_t i=0; i<dim; i++) {
    for (size_t j=0; j<dim; j++) {
      if (i != j) {
        for (T val=val_single_edge; val<=val_six; val++) {
          MyVector<T> V = V_basic;
          V(i) = val_single_edge;
          V(j) = val;
          test_vector_and_insert(V);
        }
      }
    }
  }
  std::cerr << "FindDiagramExtensions, step 6\n";
  // Considering the case of 3 edges. All have to be single edges
  SetCppIterator SCI_A(dim,3);
  for (auto & eV : SCI_A) {
    MyVector<T> V = V_basic;
    for (auto & eVal : eV)
      V(eVal) = val_single_edge;
    test_vector_and_insert(V);
  }
  // Considering the case of 4 edges. Only tilde{D4} is possible
  std::cerr << "FindDiagramExtensions, step 7\n";
  size_t n_isolated = list_isolated.size();
  SetCppIterator SCI_B(n_isolated,4);
  for (auto & eV : SCI_B) {
    MyVector<T> V = V_basic;
    for (auto & eVal : eV)
      V(list_isolated[eVal]) = val_single_edge;
    test_vector_and_insert(V);
  }
  std::cerr << "FindDiagramExtensions, step 8\n";
  std::vector<MyVector<T>> ListExtensions;
  for (auto &eEnt : SetExtensions)
    ListExtensions.push_back(eEnt);
  std::cerr << "FindDiagramExtensions, step 9\n";
  return ListExtensions;
}



template<typename T, typename Tint>
std::pair<MyMatrix<T>,MyMatrix<T>> ComputeCoxeterMatrix(MyMatrix<T> const& G, std::vector<MyVector<Tint>> const& l_root)
{
  int n = G.rows();
  auto get_scal=[&](int i, int j) -> T {
    T sum=0;
    for (int u=0; u<n; u++)
      for (int v=0; v<n; v++)
        sum += G(u,v) * l_root[i](u) * l_root[j](v);
    return sum;
  };
  T val3 = T(1) / T(4);
  T val4 = T(1) / T(2);
  T val6 = T(3) / T(4);
  auto get_val=[&](int i, int j) -> std::pair<T,T> {
    if (i == j) {
      T scal12 = get_scal(i, i);
      return {scal12, scal12};
    } else {
      T scal12 = get_scal(i, j);
      if (scal12 == 0)
        return {2,scal12};
      T scal11 = get_scal(i, i);
      T scal22 = get_scal(j, j);
      T quot = (scal12 * scal12) / (scal11 * scal22);
      if (quot == val3)
        return {3,scal12};
      if (quot == val4)
        return {4,scal12};
      if (quot == val6)
        return {6,scal12};
      std::cerr << "Failed to find matching entry\n";
      throw TerminalException{1};
    }
  };
  size_t n_root=l_root.size();
  MyMatrix<T> CoxMat(n_root,n_root);
  MyMatrix<T> ScalMat(n_root,n_root);
  for (size_t i=0; i<n_root; i++)
    for (size_t j=0; j<n_root; j++) {
      std::pair<T,T> ep = get_val(i,j);
      CoxMat(i,j) = ep.first;
      ScalMat(i,j) = ep.second;
    }
  return {CoxMat, ScalMat};
}



template<typename T>
struct Possible_Extension {
  MyVector<T> u_component;
  T res_norm;
  T e_norm;
};

template<typename T, typename Tint>
std::vector<Possible_Extension<T>> ComputePossibleExtensions(MyMatrix<T> const& G, std::vector<MyVector<Tint>> const& l_root, std::vector<T> const& l_norm, bool allow_euclidean)
{
  std::cerr << "ComputePossibleExtensions, step 1\n";
  std::pair<MyMatrix<T>,MyMatrix<T>> ep = ComputeCoxeterMatrix(G, l_root);
  std::cerr << "ComputePossibleExtensions, step 2\n";
  const MyMatrix<T> & CoxMat = ep.first;
  const MyMatrix<T> & ScalMat = ep.second;
  int dim = G.rows();
  int dim_cox = l_root.size();
  std::cerr << "ComputePossibleExtensions, step 3\n";
  std::vector<MyVector<T>> l_vect = FindDiagramExtensions(CoxMat, allow_euclidean);
  std::cerr << "ComputePossibleExtensions, step 4\n";
  T val3 = T(1) / T(4);
  T val4 = T(1) / T(2);
  T val6 = T(3) / T(4);
  auto get_cos_square=[&](T val) -> T {
    if (val == 3)
      return val3;
    if (val == 4)
      return val4;
    if (val == 6)
      return val6;
    std::cerr << "Failed to find a matching entry\n";
    throw TerminalException{1};
  };
  std::vector<Possible_Extension<T>> l_extensions;
  auto get_entry=[&](MyVector<T> const& e_vect, T const& e_norm) -> void {
    MyVector<T> l_scal(dim_cox);
    for (int i=0; i<dim_cox; i++) {
      T val = e_vect(i);
      T cos_square = get_cos_square(val);
      T scal_square = cos_square * G(i,i) * e_norm;
      std::optional<T> opt = UniversalSquareRoot(scal_square);
      if (!opt)
        return;
      T scal = - *opt;
      l_scal(i) = scal;
    }
    /* So, we have computed alpha.dot.ui = u.dot.ui
       If u = sum wi u_i then w = G^{-1} l_scal
       eNorm = w.dot.w  is the Euclidean norm of u.
     */
    MyVector<T> w = Inverse(ScalMat) * l_scal;
    T eNorm = l_scal.dot(w);
    T res_norm = e_norm - eNorm;
    if (res_norm <= 0)
      return;
    MyVector<T> u_component = ZeroVector<T>(dim);
    for (int i=0; i<dim_cox; i++)
      u_component += w(i) * UniversalVectorConversion<T,Tint>(l_root[i]);
    l_extensions.push_back({u_component, res_norm, e_norm});
  };
  std::cerr << "ComputePossibleExtensions, step 5\n";
  for (auto & e_norm : l_norm)
    for (auto & e_vect : l_vect)
      get_entry(e_vect, e_norm);
  std::cerr << "ComputePossibleExtensions, step 6\n";
  return l_extensions;
}

#endif

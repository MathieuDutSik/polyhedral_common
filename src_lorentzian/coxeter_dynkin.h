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


struct IrrCoxDyn {
  std::string type;
  size_t dim;
  int param; // For In only
};


bool IsDiagramSpherical(IrrCoxDyn const& cd)
{
  std::string type = cd.type;
  if (type == "tildeA" || type == "tildeB" || type == "tildeC" || type == "tildeD" || type == "tildeE" || type == "tildeF" || type == "tildeG")
    return false;
  return true;
}


bool IsDiagramADE(IrrCoxDyn const& cd)
{
  std::string type = cd.type;
  if (type == "A" || type == "D" || type == "E")
    return true;
  return false;
}


int GetNrVertices(IrrCoxDyn const& cd)
{
  std::string type = cd.type;
  int dim = cd.dim;
  if (type.size() > 5) {
    std::string s_res = type.substr(0,5);
    if (s_res == "tilde")
      return dim+1;
  }
  return dim;
}


std::string IrrCoxDyn_to_string(IrrCoxDyn const& cd)
{
  std::string type = cd.type;
  if (type == "A" || type == "B" || type == "C" || type == "D" || type == "E" || type == "F" || type == "G" || type == "H")
    return cd.type + "_{" + std::to_string(cd.dim) + "}";
  if (type == "I") {
    return std::string("I_2(") + std::to_string(cd.param) + ")";
  }
  if (type == "tildeA" || type == "tildeB" || type == "tildeC" || type == "tildeD" || type == "tildeE" || type == "tildeF" || type == "tildeG") {
    std::string type_red = type.substr(5,1);
    return std::string("\\tilde{") + type_red + "_{" + std::to_string(cd.dim) + "}}";
  }
  std::cerr << "cd  type=" << cd.type << " dim=" << cd.dim << " param=" << cd.param << "\n";
  std::cerr << "Failed to matching entry. Maybe bug or non-conforming input\n";
  throw TerminalException{1};
}


IrrCoxDyn string_to_IrrCoxDyn(std::string const& s)
{
  std::string s_work = s;
  auto remove_char=[&](char ec) -> void {
    s_work.erase(remove(s_work.begin(), s_work.end(), ec), s_work.end());
  };
  remove_char('_');
  remove_char('\\');
  remove_char('{');
  remove_char('}');
  remove_char('(');
  remove_char(')');
  auto recognize=[&]() -> IrrCoxDyn {
    std::vector<std::string> LS{"A", "B", "C", "D", "E", "F", "G", "H", "tildeA", "tildeB", "tildeC", "tildeD", "tildeE", "I2"};
    for (auto & eLS : LS) {
      size_t len1 = s_work.size();
      size_t len2 = eLS.size();
      if (len1 > len2) {
        std::string s_red = s_work.substr(0,len2);
        if (s_red == eLS) {
          std::string s_rem = s_work.substr(len2,len1-len2);
          int val = atoi(s_rem.c_str());
          if (eLS == "I2") {
            return IrrCoxDyn{"I", 2, val};
          } else {
            return IrrCoxDyn{eLS,size_t(val),0};
          }
        }
      }
    }
    std::cerr << "s_work=" << s_work << "\n";
    throw TerminalException{1};
  };
  IrrCoxDyn cd = recognize();
  std::string s_map = IrrCoxDyn_to_string(cd);
  if (s_map != s) {
    std::cerr << "Initial string in input is s=" << s << "\n";
    std::cerr << "Found matching to be type=" << cd.type << " dim=" << cd.dim << " param=" << cd.param << "\n";
    std::cerr << "Mapped string is s_map=" << s_map << "\n";
    throw TerminalException{1};
  }
  return cd;
}



template<typename T>
MyMatrix<T> Kernel_IrrCoxDyn_to_matrix(IrrCoxDyn const& cd)
{
  T val_comm = 2;
  T val_single_edge = 3;
  T val_four = 4; // Shows up in F4, Bn = Cn, tilde{Bn}, tilde{Cn}, tilde{F4}.
  T val_five = 5; // Shows up in H3, H4
  T val_six = 6; // Shows up in G2
  std::string type = cd.type;
  int dim = cd.dim;
  int n_vert = GetNrVertices(cd);
  MyMatrix<T> M(n_vert,n_vert);
  for (int i=0; i<n_vert; i++)
    for (int j=0; j<n_vert; j++)
      M(i,j) = val_comm;
  size_t n_assign = 0;
  auto set_v=[&](int i, int j, T val) -> void {
    M(i,j) = val;
    M(j,i) = val;
    n_assign++;
  };
  if (type == "A") {
    for (int i=1; i<dim; i++)
      set_v(i-1,i,val_single_edge);
  }
  if (type == "B" || type == "C") {
    for (int i=1; i<dim-1; i++)
      set_v(i-1,i,val_single_edge);
    set_v(dim-2,dim-1,val_four);
  }
  if (type == "D") {
    for (int i=1; i<dim-1; i++)
      set_v(i-1,i,val_single_edge);
    set_v(1,dim-1,val_single_edge);
  }
  auto set_from_triple=[&](int l1, int l2, int l3) -> void {
    auto set_line=[&](int l, int shift) -> void {
      for (int i=0; i<l; i++) {
        if (i == 0)
          set_v(0, shift + i + 1, val_single_edge);
        else
          set_v(shift + i, shift + i + 1, val_single_edge);
      }
    };
    set_line(l1, 0);
    set_line(l2, l1);
    set_line(l3, l1+l2);
  };
  if (type == "E") {
    if (dim == 6)
      set_from_triple(1,2,2);
    if (dim == 7)
      set_from_triple(1,2,3);
    if (dim == 8)
      set_from_triple(1,2,4);
  }
  if (type == "F" && dim == 4) {
    set_v(0, 1, val_single_edge);
    set_v(1, 2, val_four);
    set_v(2, 3, val_single_edge);
  }
  if (type == "G" && dim == 2) {
    set_v(0, 1, val_six);
  }
  if (type == "H") {
    for (int i=2; i<dim; i++)
      set_v(i-1,i,val_single_edge);
    set_v(0,1,val_five);
  }
  if (type == "I") {
    set_v(0,1, cd.param);
  }
  // Now the euclidean ones
  if (type == "tildeA") {
    for (int i=1; i<n_vert; i++)
      set_v(i-1,i,val_single_edge);
    set_v(0,n_vert-1,val_single_edge);
  }
  if (type == "tildeB") {
    for (int i=1; i<n_vert-2; i++)
      set_v(i-1,i,val_single_edge);
    set_v(0,n_vert-2,val_four);
    set_v(n_vert-4,n_vert-1,val_single_edge);
  }
  if (type == "tildeC") {
    for (int i=1; i<n_vert-2; i++)
      set_v(i-1,i,val_single_edge);
    set_v(0,n_vert-2,val_four);
    set_v(n_vert-3,n_vert-1,val_four);
  }
  if (type == "tildeD") {
    for (int i=1; i<n_vert-2; i++)
      set_v(i-1,i,val_single_edge);
    set_v(1,n_vert-2,val_single_edge);
    set_v(n_vert-4,n_vert-1,val_single_edge);
  }
  if (type == "tildeE") {
    if (dim == 6)
      set_from_triple(2,2,2);
    if (dim == 7)
      set_from_triple(1,3,3);
    if (dim == 8)
      set_from_triple(1,2,5);
  }
  if (type == "tildeF" && dim == 4) {
    set_v(0, 1, val_single_edge);
    set_v(1, 2, val_four);
    set_v(2, 3, val_single_edge);
    set_v(3, 4, val_single_edge);
  }
  if (type == "tildeG" && dim == 4) {
    set_v(0, 1, val_single_edge);
    set_v(1, 2, val_six);
  }
  if (type == "tildeI" && dim == 2) {
    std::cerr << "For tildeI, we should have value infinity which we cannot represent\n";
    throw TerminalException{1};
  }
  if (n_assign == 0) {
    std::cerr << "We assign 0 edges. Likely that case was not covered\n";
    std::cerr << "type=" << type << " dim=" << dim << " param=" << cd.param << "\n";
    throw TerminalException{1};
  }
  return M;
}





template<typename T>
std::optional<IrrCoxDyn> IsIrreducibleDiagramSphericalEuclidean(const MyMatrix<T>& M, bool allow_euclidean)
{
  std::cerr << "IsIrreducibleDiagramSphericalEuclidean allow_euclidean=" << allow_euclidean << " M=\n";
  WriteMatrix(std::cerr, M);
  T val_comm = 2;
  T val_single_edge = 3;
  T val_four = 4; // Shows up in F4, Bn = Cn, tilde{Bn}, tilde{Cn}, tilde{F4}.
  T val_five = 5; // Shows up in H3, H4
  T val_six = 6; // Shows up in G2
  size_t n_vert = M.rows();
  if (n_vert == 1) // Case of A1
    return IrrCoxDyn{"A",1,0};
  std::cerr << "IsIrreducibleDiagramSphericalEuclidean, step 1\n";
  std::vector<size_t> list_deg1, list_deg2, list_deg3, list_deg4, list_degN;
  std::vector<size_t> list_deg(n_vert, 0);
  size_t n_higher_edge = 0;
  GraphBitset eG(n_vert);
  std::map<T,size_t> multiplicity;
  for (size_t i=0; i<n_vert; i++) {
    size_t n_adj = 0;
    for (size_t j=0; j<n_vert; j++) {
      T val = M(i,j);
      if (i != j && val != val_comm) {
        n_adj++;
        eG.AddAdjacent(i, j);
        if (j > i) {
          multiplicity[val]++;
          if (val != val_single_edge)
            n_higher_edge++;
        }
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
  std::cerr << "IsIrreducibleDiagramSphericalEuclidean, step 2\n";
  auto get_list_adjacent=[&](size_t u) -> std::vector<size_t> {
    std::vector<size_t> LAdj;
    for (size_t j=0; j<n_vert; j++)
      if (u != j && M(u,j) != val_comm)
        LAdj.push_back(j);
    return LAdj;
  };
  auto get_value_isolated=[&](size_t u) -> T {
    for (size_t j=0; j<n_vert; j++)
      if (u != j && M(u,j) != val_comm)
        return M(u,j);
    return std::numeric_limits<T>::max(); // That case should not happen
  };
  if (list_degN.size() > 0) // vertices of degree 5 or more never occurs.
    return {};
  std::cerr << "IsIrreducibleDiagramSphericalEuclidean, step 3\n";
  if (list_deg4.size() > 0) { // Vertex of degree 4 can occur for \tilde{D4} only
    if (!allow_euclidean) // Only possibilities is not allowed, exit.
      return {};
    // We are now in Euclidean
    if (n_vert != 5) // It has to be \tilde{D4}.
      return {};
    if (list_deg4.size() != 1 || list_deg1.size() != 4 || list_deg2.size() != 0 || list_deg3.size() != 0)
      return {};
    size_t i_4 = list_deg4[0];
    for (size_t i=0; i<n_vert; i++)
      if (i != i_4)
        if (M(i, i_4) != 3)
          return {};
    return IrrCoxDyn{"tildeD", 4, 0}; // This is \tilde{D4}
  }
  std::cerr << "IsIrreducibleDiagramSphericalEuclidean, step 4\n";
  std::vector<std::vector<size_t>> ListCycles = GRAPH_FindAllCycles(eG);
  std::cerr << "|ListCycles|=" << ListCycles.size() << "\n";
  if (ListCycles.size() > 0) { // Only tilde{An} is possible.
    if (ListCycles.size() > 1) // If more than 1 cycle, then not possible
      return {};
    if (!allow_euclidean)
      return {};
    // We are now in Euclidean case
    const std::vector<size_t>& eCycle = ListCycles[0];
    if (eCycle.size() != n_vert)
      return {};
    if (list_deg2.size() != n_vert)
      return {};
    if (n_higher_edge != 0)
      return {};
    return IrrCoxDyn{"tildeA", n_vert-1, 0}; // Only tilde{An} is left as possibility
  }
  std::cerr << "IsIrreducibleDiagramSphericalEuclidean, step 5\n";
  // Now it is a tree
  if (list_deg1.size() == 2 && list_deg2.size() == n_vert - 2 && n_higher_edge == 0)
    return IrrCoxDyn{"A",n_vert,0}; // Only An is possible so ok.
  // An and tilde{An} have been covered
  if (list_deg3.size() > 2)
    return {}; // No possibility for spherical and euclidean
  if (list_deg3.size() == 2) {
    if (!allow_euclidean)
      return {};
    // We are now in Euclidean case
    if (n_higher_edge != 0)
      return {};
    for (auto & ePt : list_deg3) {
      std::vector<size_t> LAdj = get_list_adjacent(ePt);
      size_t n_deg1 = 0;
      for (auto & fPt : LAdj)
        if (list_deg[fPt] == 1)
          n_deg1++;
      if (n_deg1 != 2)
        return {};
    }
    return IrrCoxDyn{"tildeD",n_vert-1, 0}; // Only tilde{Dn} is possible
  }
  std::cerr << "IsIrreducibleDiagramSphericalEuclidean, step 6\n";
  if (list_deg3.size() == 0) { // We are in a single path.
    if (multiplicity[val_four] == 2) {
      if (n_higher_edge != 2) // There are some other higher edge, excluded
        return {};
      if (!allow_euclidean) // Only tilde{Cn} is feasible and it is Euclidean
        return {};
      for (auto & eVert : list_deg1)
        if (get_value_isolated(eVert) != val_four)
          return {};
      return IrrCoxDyn{"tildeC",n_vert-1,0}; // This is tilde{Cn}
    }
    if (multiplicity[val_four] == 1) { // Possibilities: Bn=Cn, F4 and tilde{F4} are possible
      if (n_higher_edge != 1)
        return {}; // There are other edges, excluded.
      size_t n_sing = 0;
      size_t n_four = 0;
      for (auto & eVert : list_deg1) {
        if (get_value_isolated(eVert) == val_single_edge)
          n_sing++;
        if (get_value_isolated(eVert) == val_four)
          n_four++;
      }
      if (n_sing == 2) {
        if (n_vert == 4) {
          return IrrCoxDyn{"F", 4,0}; // Only F4 is possible
        }
        if (n_vert == 5) { // Only tilde{F4} is possible. So conclude from that
          if (allow_euclidean)
            return IrrCoxDyn{"tildeF", 4, 0};
          return {};
        }
      }
      if (n_four == 1 && n_sing == 1) {
        // Only possibility is to have 4 at one extremity. This is Bn = Cn
        return IrrCoxDyn{"B",n_vert,0};
      }
      // No other possibilities
      return {};
    }
    std::cerr << "multiplicity[val_five]=" << multiplicity[val_five] << "\n";
    if (multiplicity[val_five] == 1) { // Looking for H2, H3, H4
      std::cerr << "n_vert=" << n_vert << "\n";
      if (n_vert == 2)
        return IrrCoxDyn{"H", 2,0}; // It is H2
      if (n_vert > 5)
        return {}; // No possibility
      size_t n_sing=0;
      size_t n_five=0;
      std::cerr << "list_deg1 =";
      for (auto & eVert : list_deg1) {
        std::cerr << " " << eVert;
        if (get_value_isolated(eVert) == val_single_edge)
          n_sing++;
        if (get_value_isolated(eVert) == val_five)
          n_five++;
      }
      std::cerr << "\n";
      std::cerr << "n_sing=" << n_sing << " n_five=" << n_five << "\n";
      if (n_sing == 1 && n_five == 1) { // It is H3 or H4 depending on the dimension
        if (n_vert == 3)
          return IrrCoxDyn{"H", 3, 0};
        if (n_vert == 4)
          return IrrCoxDyn{"H", 4, 0};
      }
      return {};
    }
    if (multiplicity[val_six] == 1) { // Looking for G2 or tilde{G2}
      if (n_higher_edge != 1)
        return {}; // There are other edges, excluded.
      if (n_vert == 2 || n_vert == 3) { // It is G2 or tilde{G2}
        if (n_vert == 2)
          return IrrCoxDyn{"G", 2, 0};
        if (!allow_euclidean)
          return {};
        if (n_vert == 3)
          return IrrCoxDyn{"tildeG", 2, 0};
      }
      return {};
    }
    if (n_vert == 2) {
      int param = UniversalScalarConversion<int,T>(M(0,1));
      return IrrCoxDyn{"I", 2, param}; // It is I2(n)
    }
    return {};
  }
  std::cerr << "IsIrreducibleDiagramSphericalEuclidean, step 7\n";
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
      return {};
    bool has_edge_four = false;
    for (auto & eVert : list_deg1)
      if (get_value_isolated(eVert) == val_four)
        has_edge_four = true;
    if (has_edge_four)
      return IrrCoxDyn{"tildeB",n_vert-1,0}; // It is tilde{Bn}
    return {};
  }
  std::cerr << "IsIrreducibleDiagramSphericalEuclidean, step 8\n";
  if (n_higher_edge != 0)
    return {};
  auto get_length=[&](size_t val1, size_t val2) -> size_t {
    size_t len = 1;
    size_t iter=0;
    while(true) {
      std::cerr << "get_length, passing iter=" << iter << " val1=" << val1 << " val2=" << val2 << "\n";
      iter++;
      std::vector<size_t> LVal = get_list_adjacent(val2);
      std::cerr << "val2=" << val2 << "    LVal =";
      for (auto & u : LVal)
        std::cerr << " " << u;
      std::cerr << "\n";
      if (LVal.size() == 1)
        break;
      size_t NewPt = -1;
      for (auto & eVal : LVal)
        if (eVal != val1)
          NewPt = eVal;
      val1 = val2;
      val2 = NewPt;
      len++;
    }
    return len;
  };
  std::cerr << "|list_deg3|=" << list_deg3.size() << "\n";
  size_t eCent = list_deg3[0];
  std::cerr << "eCent=" << eCent << "\n";
  std::map<size_t, size_t> map_len;
  for (auto & eAdj : get_list_adjacent(eCent)) {
    size_t len = get_length(eCent, eAdj);
    std::cerr << "eAdj=" << eAdj << " len=" << len << "\n";
    map_len[len]++;
  }
  for (auto & kv : map_len) {
    std::cerr << "kv=" << kv.first << " / " << kv.second << "\n";
  }
  std::cerr << "IsIrreducibleDiagramSphericalEuclidean, step 9\n";
  if (map_len[1] == 3) // It is D4
    return IrrCoxDyn{"D",4,0};
  if (map_len[1] == 2) // It is Dn
    return IrrCoxDyn{"D",n_vert,0};
  if (map_len[1] == 1 && map_len[2] == 2) // It is E6
    return IrrCoxDyn{"E",6,0};
  if (map_len[1] == 1 && map_len[2] == 1 && map_len[3] == 1) // It is E7
    return IrrCoxDyn{"E",7,0};
  if (map_len[1] == 1 && map_len[2] == 1 && map_len[4] == 1) // It is E8
    return IrrCoxDyn{"E",8,0};
  if (!allow_euclidean)
    return {}; // In spherical, no other possibilities left
  if (map_len[2] == 3) // It is tilde{E6}
    return IrrCoxDyn{"tildeE",6,0};
  if (map_len[1] == 1 && map_len[3] == 2) // It is tilde{E7}
    return IrrCoxDyn{"tildeE",7,0};
  if (map_len[1] == 1 && map_len[2] == 1 && map_len[5] == 1) // It is tilde{E8}
    return IrrCoxDyn{"tildeE",8,0};
  return {}; // No other possibilities left
}


template<typename T>
MyMatrix<T> IrrCoxDyn_to_matrix(IrrCoxDyn const& cd)
{
  MyMatrix<T> M = Kernel_IrrCoxDyn_to_matrix<T>(cd);
  bool allow_euclidean = true;
  std::optional<IrrCoxDyn> opt = IsIrreducibleDiagramSphericalEuclidean(M, allow_euclidean);
  if (opt) {
    IrrCoxDyn cd2 = *opt;
    if (cd.type != cd2.type || cd.dim != cd2.dim || cd.param != cd2.param) {
      std::cerr << "M=\n";
      WriteMatrix(std::cerr, M);
      std::cerr << "cd   type=" << cd.type << " dim" << cd.dim << " param=" << cd.param << "\n";
      std::cerr << "cd2  type=" << cd2.type << " dim" << cd2.dim << " param=" << cd2.param << "\n";
      std::cerr << "The recognition of the matrix did not yield the original Coxeter-Dynkin diagram\n";
      throw TerminalException{1};
    }
    return M;
  }
  std::cerr << "The created matrix was not recognized. Some bug somewhere\n";
  throw TerminalException{1};
}




template<typename T>
std::optional<std::vector<IrrCoxDyn>> IsDiagramSphericalEuclidean(const MyMatrix<T>& M, const bool& allow_euclidean)
{
  T val_comm = 2;
  size_t dim = M.rows();
  std::cerr << "IsDiagramSphericalEuclidean dim=" << dim << "\n";
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
  std::vector<IrrCoxDyn> l_cd;
  for (auto & eConn : LConn) {
    size_t dim_res=eConn.size();
    std::cerr << "dim_res=" << dim_res << "\n";
    MyMatrix<T> Mres(dim_res, dim_res);
    for (size_t i=0; i<dim_res; i++)
      for (size_t j=0; j<dim_res; j++)
        Mres(i,j) = M(eConn[i], eConn[j]);
    std::cerr << "Before IsIrreducibleDiagramSphericalEuclidean\n";
    std::optional<IrrCoxDyn> opt = IsIrreducibleDiagramSphericalEuclidean(Mres, allow_euclidean);
    std::cerr << "After IsIrreducibleDiagramSphericalEuclidean\n";
    if (opt) {
      IrrCoxDyn cd = *opt;
      l_cd.push_back(cd);
      std::cerr << "symb=" << IrrCoxDyn_to_string(cd) << "\n";
    } else {
      std::cerr << "Answer is false\n";
      return {};
    }
  }
  return l_cd;
}



template<typename T>
MyMatrix<T> string_to_coxdyn_matrix(std::string const& str)
{
  std::vector<std::string> LStr = STRING_Split(str, "+");
  std::vector<MyMatrix<T>> LMat;
  size_t n_vert = 0;
  for (auto & s1 : LStr) {
    std::string s2 = s1;
    s2.erase(remove(s2.begin(), s2.end(), ' '), s2.end());
    IrrCoxDyn cd = string_to_IrrCoxDyn(s2);
    MyMatrix<T> M = IrrCoxDyn_to_matrix<T>(cd);
    LMat.push_back(M);
    n_vert += M.rows();
  }
  MyMatrix<T> Mret(n_vert,n_vert);
  for (int i=0; i<n_vert; i++)
    for (int j=0; j<n_vert; j++)
      Mret(i,j) = 2;
  size_t shift = 0;
  for (auto & M : LMat) {
    size_t len = M.rows();
    for (size_t i=0; i<len; i++)
      for (size_t j=0; j<len; j++)
        Mret(shift + i, shift + j) = M(i,j);
    shift += len;
  }
  return Mret;
}

template<typename T>
std::string coxdyn_matrix_to_string(MyMatrix<T> const& M)
{
  bool allow_euclidean = true;
  std::optional<std::vector<IrrCoxDyn>> opt = IsDiagramSphericalEuclidean(M, allow_euclidean);
  if (opt) {
    const std::vector<IrrCoxDyn> & l_irr = *opt;
    std::string str = IrrCoxDyn_to_string(l_irr[0]);
    for (int i=1; i<l_irr.size(); i++)
      str += "+" + IrrCoxDyn_to_string(l_irr[i]);
    return str;
  }
  std::cerr << "The recognition failed so coxdyn_matrix_to_string cannot work\n";
  throw TerminalException{1};
}





template<typename T>
MyMatrix<T> ExtendMatrix(MyMatrix<T> const& M, MyVector<T> const& V)
{
  size_t len = M.rows();
  MyMatrix<T> Mret(len+1, len+1);
  for (int i=0; i<len; i++)
    for (int j=0; j<len; j++)
      Mret(i,j) = M(i,j);
  for (int i=0; i<len; i++) {
    Mret(len,i) = V(i);
    Mret(i,len) = V(i);
  }
  return Mret;
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
    if (IsDiagramSphericalEuclidean(Mtest, allow_euclidean))
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

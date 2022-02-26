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

//#define DEBUG_COXETER_DYNKIN_COMBINATORICS


template<typename T>
struct IrrCoxDyn {
  std::string type;
  size_t dim;
  T param; // For In only
};


template<typename T>
bool IsDiagramSpherical(IrrCoxDyn<T> const& cd)
{
  std::string type = cd.type;
  if (type == "tildeA" || type == "tildeB" || type == "tildeC" || type == "tildeD" || type == "tildeE" || type == "tildeF" || type == "tildeG")
    return false;
  return true;
}


template<typename T>
bool IsDiagramSimplyLaced(IrrCoxDyn<T> const& cd)
{
  std::string type = cd.type;
  if (type == "A" || type == "D" || type == "E")
    return true;
  return false;
}


template<typename T>
bool IsDiagramIntegerLorentzianFeasible(IrrCoxDyn<T> const& cd)
{
  std::string type = cd.type;
  if (type == "I") {
    if (cd.dim == 1 && cd.param == practical_infinity<T>())
      return true;
    return false;
  }
  if (type == "H")
    return false;
  return true;
}


template<typename T>
int GetNrVertices(IrrCoxDyn<T> const& cd)
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


template<typename T>
std::string IrrCoxDyn_to_string(IrrCoxDyn<T> const& cd)
{
  std::string type = cd.type;
  if (type == "A" || type == "B" || type == "C" || type == "D" || type == "E" || type == "F" || type == "G" || type == "H")
    return cd.type + "_{" + std::to_string(cd.dim) + "}";
  if (type == "I") {
    if (cd.param == practical_infinity<T>())
      return "I_2(infinity)";
    return std::string("I_2(") + std::to_string(cd.param) + ")";
  }
  if (type == "tildeA" || type == "tildeB" || type == "tildeC" || type == "tildeD" || type == "tildeE" || type == "tildeF" || type == "tildeG") {
    std::string type_red = type.substr(5,1);
    return std::string("\\tilde{") + type_red + "_{" + std::to_string(cd.dim) + "}}";
  }
  if (type == "tildeI")
    return std::string("\\tilde{I_{1}(infinity)}");
  std::cerr << "cd  type=" << cd.type << " dim=" << cd.dim << " param=" << cd.param << "\n";
  std::cerr << "Failed to matching entry. Maybe bug or non-conforming input\n";
  throw TerminalException{1};
}


template<typename T>
IrrCoxDyn<T> string_to_IrrCoxDyn(std::string const& s)
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
  auto recognize=[&]() -> IrrCoxDyn<T> {
    std::vector<std::string> LS{"A", "B", "C", "D", "E", "F", "G", "H", "tildeA", "tildeB", "tildeC", "tildeD", "tildeE", "I2", "tildeI1"};
    for (auto & eLS : LS) {
      size_t len1 = s_work.size();
      size_t len2 = eLS.size();
      if (len1 > len2) {
        std::string s_red = s_work.substr(0,len2);
        if (s_red == eLS) {
          std::string s_rem = s_work.substr(len2,len1-len2);
          if (eLS == "tildeI1") {
            if (s_rem == "infinity")
              return IrrCoxDyn<T>{"tildeI", 1, practical_infinity<T>()};
          } else {
            if (eLS == "I2") {
              T val = ParseScalar<T>(s_rem);
              return IrrCoxDyn<T>{"I", 2, val};
            } else {
              size_t val = ParseScalar<size_t>(s_rem);
              return IrrCoxDyn<T>{eLS,val,0};
            }
          }
        }
      }
    }
    std::cerr << "s_work=" << s_work << "\n";
    throw TerminalException{1};
  };
  IrrCoxDyn<T> cd = recognize();
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
MyMatrix<T> Kernel_IrrCoxDyn_to_matrix(IrrCoxDyn<T> const& cd)
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
    if (dim == 1) // That one has zero edges. So, we need to increase so as to pass.
      n_assign++;
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
  if (type == "H" && dim >= 2 && dim <= 4) {
    for (int i=2; i<dim; i++)
      set_v(i-1,i,val_single_edge);
    set_v(0,1,val_five);
  }
  if (type == "I" && dim == 2) {
    if (cd.param < 7) {
      std::cerr << "For the I case the values 2, 3, 4, 5, 6 are covered by A1+A1 , A2 , B2 , H2 and G2\n";
      throw TerminalException{1};
    }
    if (cd.param == practical_infinity<T>()) {
      std::cerr << "The I2(infinity) is not allowed\n";
      throw TerminalException{1};
    }
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
  if (type == "tildeG" && dim == 2) {
    set_v(0, 1, val_single_edge);
    set_v(1, 2, val_six);
  }
  if (type == "tildeI" && dim == 1) {
    set_v(0, 1, practical_infinity<T>());
  }
  if (n_assign == 0) {
    std::cerr << "We assign 0 edges. Likely that case was not covered\n";
    std::cerr << "type=" << type << " dim=" << dim << " param=" << cd.param << "\n";
    throw TerminalException{1};
  }
  return M;
}












template<typename T>
std::pair<std::vector<size_t>,IrrCoxDyn<T>> FindExtensionVerticesIrreducibleSphericalEuclideanDiagram(const MyMatrix<T>& M)
{
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram M=\n";
  WriteMatrix(std::cerr, M);
#endif
  T val_comm = 2;
  T val_single_edge = 3;
  T val_four = 4; // Shows up in F4, Bn = Cn, tilde{Bn}, tilde{Cn}, tilde{F4}.
  T val_five = 5; // Shows up in H3, H4
  T val_six = 6; // Shows up in G2
  size_t n_vert = M.rows();
  if (n_vert == 1) // Case of A1
    return { {0}, IrrCoxDyn<T>{"A",1,0} };
  auto error=[&]() -> void {
    std::cerr << "The diagram is not spherical or euclidean. Contrary to hypothesis\n";
    throw TerminalException{1};
  };
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 1\n";
#endif
  std::vector<size_t> list_deg1, list_deg2, list_deg3, list_deg4, list_degN;
  std::vector<size_t> list_deg(n_vert, 0);
  size_t n_higher_edge = 0;
  GraphBitset eG(n_vert);
  std::map<T,size_t> multiplicity;
  std::vector<std::vector<size_t>> LLAdj;
  std::vector<T> list_isolated_adjacent(n_vert);
  for (size_t i=0; i<n_vert; i++) {
    size_t n_adj = 0;
    std::vector<size_t> LAdj;
    T e_isolated_adjacent = std::numeric_limits<T>::max();
    for (size_t j=0; j<n_vert; j++) {
      T val = M(i,j);
      if (i != j && val != val_comm) {
        n_adj++;
        eG.AddAdjacent(i, j);
        LAdj.push_back(j);
        e_isolated_adjacent = val;
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
    LLAdj.push_back(LAdj);
    list_isolated_adjacent[i] = e_isolated_adjacent;
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 2\n";
#endif
  if (list_degN.size() > 0) // vertices of degree 5 or more never occurs.
    error();
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 3\n";
#endif
  if (list_deg4.size() > 0) { // Vertex of degree 4 can occur for \tilde{D4} only
    // We are now in Euclidean
    if (n_vert != 5) // It has to be \tilde{D4}.
      error();
    if (list_deg4.size() != 1 || list_deg1.size() != 4 || list_deg2.size() != 0 || list_deg3.size() != 0)
      error();
    size_t i_4 = list_deg4[0];
    for (size_t i=0; i<n_vert; i++)
      if (i != i_4)
        if (M(i, i_4) != 3)
          error();;
    return { {}, IrrCoxDyn<T>{"tildeD", 4, 0} }; // This is \tilde{D4}. No extension vertices.
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 4\n";
#endif
  std::vector<std::vector<size_t>> ListCycles = GRAPH_FindAllCycles(eG);
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "|ListCycles|=" << ListCycles.size() << "\n";
#endif
  if (ListCycles.size() > 0) { // Only tilde{An} is possible.
    if (ListCycles.size() > 1) // If more than 1 cycle, then not possible
      error();
    // We are now in Euclidean case
    const std::vector<size_t>& eCycle = ListCycles[0];
    if (eCycle.size() != n_vert)
      error();
    if (list_deg2.size() != n_vert)
      error();
    if (n_higher_edge != 0)
      error();
    return { {}, IrrCoxDyn<T>{"tildeA", n_vert-1, 0} }; // This is tilde{An}. No extension possible
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 5\n";
#endif
  // Now it is a tree
  if (list_deg1.size() == 2 && list_deg2.size() == n_vert - 2 && n_higher_edge == 0) {
    // This is A_n. But finding the possible vertices for extension is not going to be easy:
    std::vector<size_t> AllExtens;
    if (n_vert <= 7) { // In that case all the vertices can potentially be used
      for (size_t i=0; i<n_vert; i++)
        AllExtens.push_back(i);
      return { AllExtens, IrrCoxDyn<T>{"A",n_vert,0} };
    }
    auto append=[&](size_t e_vert, size_t e_lev) -> void {
      size_t evert1 = e_vert;
      size_t evert2 = LLAdj[evert1][0];
      for (size_t i_lev=0; i_lev<e_lev; i_lev++) {
        AllExtens.push_back(evert1);
        size_t evert2sav = evert2;
        auto set_evert2=[&]() -> void { // It all collapses for A2, but it does not matter here
          for (auto & eAdj : LLAdj[evert2]) {
            if (eAdj != evert1) {
              evert2 = eAdj;
              return;
            }
          }
        };
        set_evert2();
        evert1 = evert2sav;
      }
    };
    if (n_vert == 8) { // We have the extension to D9, but also to tilde{E8}
      append(list_deg1[0], 3);
      append(list_deg1[1], 3);
    } else { // Only A(n+1) and D(n+1) are possible
      append(list_deg1[0], 2);
      append(list_deg1[1], 2);
    }
    return { AllExtens, IrrCoxDyn<T>{"A",n_vert,0} };
  }
  // An and tilde{An} have been covered
  if (list_deg3.size() > 2)
    error(); // No possibility for spherical and euclidean
  if (list_deg3.size() == 2) {
    // We are now in Euclidean case
    if (n_higher_edge != 0)
      error();
    for (auto & ePt : list_deg3) {
      std::vector<size_t> const& LAdj = LLAdj[ePt];
      size_t n_deg1 = 0;
      for (auto & fPt : LAdj)
        if (list_deg[fPt] == 1)
          n_deg1++;
      if (n_deg1 != 2)
        error();
    }
    return { {}, IrrCoxDyn<T>{"tildeD",n_vert-1, 0} }; // For tilde{Dn} no extension is possible
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 6 |list_deg3|=" << list_deg3.size() << "\n";
#endif
  if (list_deg3.size() == 0) { // We are in a single path.
    if (multiplicity[val_four] == 2) {
      if (n_higher_edge != 2) // There are some other higher edge, excluded
        error();
      // Only tilde{Cn} is feasible and it is Euclidean
      for (auto & eVert : list_deg1)
        if (list_isolated_adjacent[eVert] != val_four)
          error();
      return { {}, IrrCoxDyn<T>{"tildeC",n_vert-1,0} }; // No extension possible for tilde{Cn}
    }
    if (multiplicity[val_four] == 1) { // Possibilities: Bn=Cn, F4, tilde{F4}, and I2(4) are possible
      if (n_vert == 2) {
        return { {0,1}, IrrCoxDyn<T>{"B", 2, 0} }; // For I2(4), There is extension on both vertices
      }
      if (n_higher_edge != 1)
        error(); // There are other edges, excluded.
      size_t n_sing = 0;
      size_t n_four = 0;
      size_t e_vert_sing;
      size_t e_vert_four;
      for (auto & eVert : list_deg1) {
        if (list_isolated_adjacent[eVert] == val_single_edge) {
          n_sing++;
          e_vert_sing = eVert;
        }
        if (list_isolated_adjacent[eVert] == val_four) {
          n_four++;
          e_vert_four = eVert;
        }
      }
      if (n_sing == 2) {
        if (n_vert == 4) {
          return { list_deg1, IrrCoxDyn<T>{"F", 4,0} }; // For F4 there are extension on both extremities.
        }
        if (n_vert == 5) { // Only tilde{F4} is possible. So conclude from that
          return { {}, IrrCoxDyn<T>{"tildeF", 4, 0} }; // No extension is possible here
        }
      }
      if (n_four == 1 && n_sing == 1) {
        // Only possibility is to have 4 at one extremity. This is Bn = Cn
        // But this is complicated to cover all cases.
        if (n_vert <= 3) { // All vertices are possible then
          return { {0,1,2}, IrrCoxDyn<T>{"B",n_vert,0} };
        }
        std::vector<size_t> AllExtens;
        if (n_vert == 4) { // For B4, we can extend to tilde{F4}
          AllExtens.push_back(e_vert_four);
        }
        AllExtens.push_back(e_vert_sing); // Extension to B(n+1)
        size_t e_adj = LLAdj[e_vert_sing][0];
        AllExtens.push_back(e_adj); // Extension to tilde(Bn)
        return { AllExtens, IrrCoxDyn<T>{"B",n_vert,0} };
      }
      // No other possibilities
      error();
    }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
    std::cerr << "multiplicity[val_five]=" << multiplicity[val_five] << "\n";
#endif
    if (multiplicity[val_five] == 1) { // Looking for H2, H3, H4
      if (n_vert == 2)
        return { {0,1}, IrrCoxDyn<T>{"H", 2,0} }; // Extensions to H3 possible on both ends.
      if (n_vert > 5)
        error(); // No possibility
      size_t n_sing=0;
      size_t n_five=0;
      size_t e_vert_sing;
      for (auto & eVert : list_deg1) {
        if (list_isolated_adjacent[eVert] == val_single_edge) {
          n_sing++;
          e_vert_sing = eVert;
        }
        if (list_isolated_adjacent[eVert] == val_five)
          n_five++;
      }
      if (n_sing == 1 && n_five == 1) { // It is H3 or H4 depending on the dimension
        if (n_vert == 3)
          return { {e_vert_sing}, IrrCoxDyn<T>{"H", 3, 0} }; // Extension to H4 is the only possibility
        if (n_vert == 4)
          return { {}, IrrCoxDyn<T>{"H", 4, 0} }; // H4 cannot be extended.
      }
      error();
    }
    if (multiplicity[val_six] == 1) { // Looking for G2 or tilde{G2}
      if (n_higher_edge != 1)
        error(); // There are other edges, excluded.
      if (n_vert == 2)
        return { {0,1}, IrrCoxDyn<T>{"G", 2, 0} }; // Extension to tilde{G2} possible on both vertices
      if (n_vert == 3)
        return { {}, IrrCoxDyn<T>{"tildeG", 2, 0} }; // tilde{G2} cannot be extended
      error();
    }
    if (n_vert == 2) {
      T param = M(0,1);
      if (param == practical_infinity<T>()) {
        return { {}, IrrCoxDyn<T>{"tildeI", 1, param} }; // I1(infinity) cannot be extended
      }
      return { {}, IrrCoxDyn<T>{"I", 2, param} }; // We have now param > 6 and no extension is possible for I2(param) in that case
    }
    error();
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 7\n";
#endif
  if (list_deg3.size() != 1) {
    std::cerr << "We should hqve just one vertex of degree 3\n";
    throw TerminalException{1};
  }
  size_t eCent = list_deg3[0];
  // Now just one vertex of degree 3.
  if (multiplicity[val_four] == 1) { // Possibility tilde{Bn}
    std::vector<size_t> const& LAdj = LLAdj[eCent];
    size_t n_sing = 0;
    for (auto & eAdj : LAdj)
      if (list_deg[eAdj] == 1)
        n_sing++;
    if (n_sing == 3) { // It is actualla tilde{B3}. No extension possible
      return { {}, IrrCoxDyn<T>{"tildeB",3,0} };
    }
    if (n_sing != 2)
      error();
    bool has_edge_four = false;
    for (auto & eVert : list_deg1)
      if (list_isolated_adjacent[eVert] == val_four)
        has_edge_four = true;
    if (has_edge_four) {
      return { {}, IrrCoxDyn<T>{"tildeB",n_vert-1,0} }; // tilde{Bn} cannot be extended
    }
    error();
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 8\n";
#endif
  if (n_higher_edge != 0)
    error();
  auto get_length=[&](size_t val1, size_t val2) -> std::pair<size_t,size_t> {
    size_t len = 1;
    size_t iter=0;
    while(true) {
      iter++;
      std::vector<size_t> const& LVal = LLAdj[val2];
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
    return {len,val2};
  };
  // map from the length to the list of vertices of that length
  std::map<size_t, std::vector<size_t>> map_len;
  for (auto & eAdj : LLAdj[eCent]) {
    std::pair<size_t,size_t> ep = get_length(eCent, eAdj);
    map_len[ep.first].push_back(ep.second);
  }
  if (map_len[1].size() == 3) { // It is D4
    return { {0,1,2,3}, IrrCoxDyn<T>{"D",4,0} }; // Possible extensions are to D5 from the 3 vertices of degree 1 and to tilde{D4} from the vertex of degree 3
  }
  std::vector<size_t> AllExtens;
  if (map_len[1].size() == 2) { // It is Dn with n >= 5
    size_t oth_len = n_vert - 3;
    size_t oth_vert = map_len[oth_len][0];
    AllExtens.push_back(oth_vert); // Extension to D(n+1) always possible
    size_t oth_vert_b = LLAdj[oth_vert][0];
    AllExtens.push_back(oth_vert_b); // Extension to tilde(Dn) always possible
    if (n_vert <= 8) { // Extension to E6, E7, E8, tilde{E8} are possible.
      AllExtens.push_back(map_len[1][0]);
      AllExtens.push_back(map_len[1][1]);
    }
    return { AllExtens, IrrCoxDyn<T>{"D",n_vert,0} };
  }
  if (map_len[1].size() == 1 && map_len[2].size() == 2) { // It is E6
    return { {map_len[1][0], map_len[2][0], map_len[2][1]}, IrrCoxDyn<T>{"E",6,0} }; // Extensions to E7 and to tilde{E6}
  }
  if (map_len[1].size() == 1 && map_len[2].size() == 1 && map_len[3].size() == 1) { // It is E7
    return { {map_len[3][0], map_len[2][0]}, IrrCoxDyn<T>{"E",7,0} }; // Extensions to E8 and to tilde{E7}
  }
  if (map_len[1].size() == 1 && map_len[2].size() == 1 && map_len[4].size() == 1) { // It is E8
    return { {map_len[4][0]}, IrrCoxDyn<T>{"E",8,0} }; // Extension to tilde{E8} is the only possibility
  }
  // In spherical, no other possibilities left
  if (map_len[2].size() == 3) { // It is tilde{E6}
    return { {}, IrrCoxDyn<T>{"tildeE",6,0} }; // No extension is possible
  }
  if (map_len[1].size() == 1 && map_len[3].size() == 2) { // It is tilde{E7}
    return { {}, IrrCoxDyn<T>{"tildeE",7,0} }; // No extension is possible
  }
  if (map_len[1].size() == 1 && map_len[2].size() == 1 && map_len[5].size() == 1) { // It is tilde{E8}
    return { {}, IrrCoxDyn<T>{"tildeE",8,0} }; // No extension is possible
  }
  error(); // No other possibilities left
}






template<typename T>
std::optional<IrrCoxDyn<T>> RecognizeIrreducibleSphericalEuclideanDiagram(const MyMatrix<T>& M)
{
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram M=\n";
  WriteMatrix(std::cerr, M);
#endif
  T val_comm = 2;
  T val_single_edge = 3;
  T val_four = 4; // Shows up in F4, Bn = Cn, tilde{Bn}, tilde{Cn}, tilde{F4}.
  T val_five = 5; // Shows up in H3, H4
  T val_six = 6; // Shows up in G2
  size_t n_vert = M.rows();
  if (n_vert == 1) // Case of A1
    return IrrCoxDyn<T>{"A",1,0};
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 1\n";
#endif
  std::vector<size_t> list_deg1, list_deg2, list_deg3, list_deg4, list_degN;
  std::vector<size_t> list_deg(n_vert, 0);
  size_t n_higher_edge = 0;
  GraphBitset eG(n_vert);
  std::map<T,size_t> multiplicity;
  std::vector<std::vector<size_t>> LLAdj;
  std::vector<T> list_isolated_adjacent(n_vert);
  for (size_t i=0; i<n_vert; i++) {
    size_t n_adj = 0;
    std::vector<size_t> LAdj;
    T e_isolated_adjacent = std::numeric_limits<T>::max();
    for (size_t j=0; j<n_vert; j++) {
      T val = M(i,j);
      if (i != j && val != val_comm) {
        n_adj++;
        eG.AddAdjacent(i, j);
        LAdj.push_back(j);
        e_isolated_adjacent = val;
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
    LLAdj.push_back(LAdj);
    list_isolated_adjacent[i] = e_isolated_adjacent;
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 2\n";
#endif
  if (list_degN.size() > 0) // vertices of degree 5 or more never occurs.
    return {};
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 3\n";
#endif
  if (list_deg4.size() > 0) { // Vertex of degree 4 can occur for \tilde{D4} only
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
    return IrrCoxDyn<T>{"tildeD", 4, 0}; // This is \tilde{D4}
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 4\n";
#endif
  std::vector<std::vector<size_t>> ListCycles = GRAPH_FindAllCycles(eG);
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "|ListCycles|=" << ListCycles.size() << "\n";
#endif
  if (ListCycles.size() > 0) { // Only tilde{An} is possible.
    if (ListCycles.size() > 1) // If more than 1 cycle, then not possible
      return {};
    // We are now in Euclidean case
    const std::vector<size_t>& eCycle = ListCycles[0];
    if (eCycle.size() != n_vert)
      return {};
    if (list_deg2.size() != n_vert)
      return {};
    if (n_higher_edge != 0)
      return {};
    return IrrCoxDyn<T>{"tildeA", n_vert-1, 0}; // Only tilde{An} is left as possibility
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 5\n";
#endif
  // Now it is a tree
  if (list_deg1.size() == 2 && list_deg2.size() == n_vert - 2 && n_higher_edge == 0)
    return IrrCoxDyn<T>{"A",n_vert,0}; // Only An is possible so ok.
  // An and tilde{An} have been covered
  if (list_deg3.size() > 2)
    return {}; // No possibility for spherical and euclidean
  if (list_deg3.size() == 2) {
    // We are now in Euclidean case
    if (n_higher_edge != 0)
      return {};
    for (auto & ePt : list_deg3) {
      std::vector<size_t> const& LAdj = LLAdj[ePt];
      size_t n_deg1 = 0;
      for (auto & fPt : LAdj)
        if (list_deg[fPt] == 1)
          n_deg1++;
      if (n_deg1 != 2)
        return {};
    }
    return IrrCoxDyn<T>{"tildeD",n_vert-1, 0}; // Only tilde{Dn} is possible
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 6 |list_deg3|=" << list_deg3.size() << "\n";
#endif
  if (list_deg3.size() == 0) { // We are in a single path.
    if (multiplicity[val_four] == 2) {
      if (n_higher_edge != 2) // There are some other higher edge, excluded
        return {};
      // Only tilde{Cn} is feasible and it is Euclidean
      for (auto & eVert : list_deg1)
        if (list_isolated_adjacent[eVert] != val_four)
          return {};
      return IrrCoxDyn<T>{"tildeC",n_vert-1,0}; // This is tilde{Cn}
    }
    if (multiplicity[val_four] == 1) { // Possibilities: Bn=Cn, F4, tilde{F4}, and B2 are possible
      if (n_vert == 2) {
        T param = 4;
        return IrrCoxDyn<T>{"B", 2, 0}; // It is I2(4)
      }
      if (n_higher_edge != 1)
        return {}; // There are other edges, excluded.
      size_t n_sing = 0;
      size_t n_four = 0;
      for (auto & eVert : list_deg1) {
        if (list_isolated_adjacent[eVert] == val_single_edge)
          n_sing++;
        if (list_isolated_adjacent[eVert] == val_four)
          n_four++;
      }
      if (n_sing == 2) {
        if (n_vert == 4) {
          return IrrCoxDyn<T>{"F", 4,0}; // Only F4 is possible
        }
        if (n_vert == 5) { // Only tilde{F4} is possible. So conclude from that
          return IrrCoxDyn<T>{"tildeF", 4, 0};
        }
      }
      if (n_four == 1 && n_sing == 1) {
        // Only possibility is to have 4 at one extremity. This is Bn = Cn
        return IrrCoxDyn<T>{"B",n_vert,0};
      }
      // No other possibilities
      return {};
    }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
    std::cerr << "multiplicity[val_five]=" << multiplicity[val_five] << "\n";
#endif
    if (multiplicity[val_five] == 1) { // Looking for H2, H3, H4
      if (n_vert == 2)
        return IrrCoxDyn<T>{"H", 2,0}; // It is H2
      if (n_vert > 5)
        return {}; // No possibility
      size_t n_sing=0;
      size_t n_five=0;
      for (auto & eVert : list_deg1) {
        if (list_isolated_adjacent[eVert] == val_single_edge)
          n_sing++;
        if (list_isolated_adjacent[eVert] == val_five)
          n_five++;
      }
      if (n_sing == 1 && n_five == 1) { // It is H3 or H4 depending on the dimension
        if (n_vert == 3)
          return IrrCoxDyn<T>{"H", 3, 0};
        if (n_vert == 4)
          return IrrCoxDyn<T>{"H", 4, 0};
      }
      return {};
    }
    if (multiplicity[val_six] == 1) { // Looking for G2 or tilde{G2}
      if (n_higher_edge != 1)
        return {}; // There are other edges, excluded.
      if (n_vert == 2)
        return IrrCoxDyn<T>{"G", 2, 0};
      if (n_vert == 3)
        return IrrCoxDyn<T>{"tildeG", 2, 0};
      return {};
    }
    if (n_vert == 2) {
      T param = M(0,1);
      if (param == practical_infinity<T>()) {
        return IrrCoxDyn<T>{"tildeI", 1, param}; // It is I1(infinity)
      }
      return IrrCoxDyn<T>{"I", 2, param}; // It is I2(n)
    }
    return {};
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 7\n";
#endif
  if (list_deg3.size() != 1) {
    std::cerr << "We should hqve just one vertex of degree 3\n";
    throw TerminalException{1};
  }
  size_t eCent = list_deg3[0];
  // Now just one vertex of degree 3.
  if (multiplicity[val_four] == 1) { // Possibility tilde{Bn}
    std::vector<size_t> const& LAdj = LLAdj[eCent];
    size_t n_sing = 0;
    for (auto & eAdj : LAdj)
      if (list_deg[eAdj] == 1)
        n_sing++;
    if (n_sing == 3) { // All neighbors are single. Only one possibility
      return IrrCoxDyn<T>{"tildeB",3,0}; // It is tilde{B3}
    }
    if (n_sing != 2)
      return {};
    bool has_edge_four = false;
    for (auto & eVert : list_deg1)
      if (list_isolated_adjacent[eVert] == val_four)
        has_edge_four = true;
    if (has_edge_four)
      return IrrCoxDyn<T>{"tildeB",n_vert-1,0}; // It is tilde{Bn}
    return {};
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 8\n";
#endif
  if (n_higher_edge != 0)
    return {};
  auto get_length=[&](size_t val1, size_t val2) -> size_t {
    size_t len = 1;
    size_t iter=0;
    while(true) {
      iter++;
      std::vector<size_t> const& LVal = LLAdj[val2];
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
  std::map<size_t, size_t> map_len;
  for (auto & eAdj : LLAdj[eCent]) {
    size_t len = get_length(eCent, eAdj);
    map_len[len]++;
  }
  if (map_len[1] == 3) // It is D4
    return IrrCoxDyn<T>{"D",4,0};
  if (map_len[1] == 2) // It is Dn
    return IrrCoxDyn<T>{"D",n_vert,0};
  if (map_len[1] == 1 && map_len[2] == 2) // It is E6
    return IrrCoxDyn<T>{"E",6,0};
  if (map_len[1] == 1 && map_len[2] == 1 && map_len[3] == 1) // It is E7
    return IrrCoxDyn<T>{"E",7,0};
  if (map_len[1] == 1 && map_len[2] == 1 && map_len[4] == 1) // It is E8
    return IrrCoxDyn<T>{"E",8,0};
  // In spherical, no other possibilities left
  if (map_len[2] == 3) // It is tilde{E6}
    return IrrCoxDyn<T>{"tildeE",6,0};
  if (map_len[1] == 1 && map_len[3] == 2) // It is tilde{E7}
    return IrrCoxDyn<T>{"tildeE",7,0};
  if (map_len[1] == 1 && map_len[2] == 1 && map_len[5] == 1) // It is tilde{E8}
    return IrrCoxDyn<T>{"tildeE",8,0};
  return {}; // No other possibilities left
}


template<typename T>
MyMatrix<T> IrrCoxDyn_to_matrix(IrrCoxDyn<T> const& cd)
{
  MyMatrix<T> M = Kernel_IrrCoxDyn_to_matrix<T>(cd);
  std::optional<IrrCoxDyn<T>> opt = RecognizeIrreducibleSphericalEuclideanDiagram(M);
  if (opt) {
    IrrCoxDyn<T> cd2 = *opt;
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
std::vector<std::vector<size_t>> GetIrreducibleComponents(const MyMatrix<T>& M)
{
  T val_comm = 2;
  size_t dim = M.rows();
  GraphBitset eG(dim);
  //  std::cerr << "LEdge =";
  for (size_t i=0; i<dim; i++) {
    for (size_t j=i+1; j<dim; j++) {
      if (M(i,j) != val_comm) {
        eG.AddAdjacent(i,j);
        eG.AddAdjacent(j,i);
        //        std::cerr << " [" << i << "/" << j << "]";
      }
    }
  }
  //  std::cerr << "\n";
  return ConnectedComponents_set(eG);
}

template<typename T>
std::optional<std::vector<IrrCoxDyn<T>>> RecognizeSphericalEuclideanDiagram(const MyMatrix<T>& M)
{
  std::vector<std::vector<size_t>> LConn = GetIrreducibleComponents(M);
  std::vector<IrrCoxDyn<T>> l_cd;
  for (auto & eConn : LConn) {
    size_t dim_res=eConn.size();
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
    std::cerr << "eConn =";
    for (auto & val : eConn)
      std::cerr << " " << val;
    std::cerr << "\n";
#endif
    MyMatrix<T> Mres(dim_res, dim_res);
    for (size_t i=0; i<dim_res; i++)
      for (size_t j=0; j<dim_res; j++)
        Mres(i,j) = M(eConn[i], eConn[j]);
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
    std::cerr << "Mres=\n";
    WriteMatrix(std::cerr, Mres);
#endif
    std::optional<IrrCoxDyn<T>> opt = RecognizeIrreducibleSphericalEuclideanDiagram(Mres);
    if (opt) {
      IrrCoxDyn<T> cd = *opt;
      l_cd.push_back(cd);
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
      std::cerr << "symb=" << IrrCoxDyn_to_string(cd) << "\n";
#endif
    } else {
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
      std::cerr << "Answer is false\n";
#endif
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
    IrrCoxDyn<T> cd = string_to_IrrCoxDyn<T>(s2);
    MyMatrix<T> M = IrrCoxDyn_to_matrix<T>(cd);
    LMat.push_back(M);
    n_vert += M.rows();
  }
  MyMatrix<T> Mret(n_vert,n_vert);
  for (size_t i=0; i<n_vert; i++)
    for (size_t j=0; j<n_vert; j++)
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
  std::optional<std::vector<IrrCoxDyn<T>>> opt = RecognizeSphericalEuclideanDiagram(M);
  if (opt) {
    const std::vector<IrrCoxDyn<T>> & l_irr = *opt;
    std::vector<std::string> l_str;
    for (auto & icd : l_irr)
      l_str.push_back(IrrCoxDyn_to_string(icd));
    std::sort(l_str.begin(), l_str.end());
    //
    std::string str = l_str[0];
    for (size_t i=1; i<l_str.size(); i++)
      str += "+" + l_str[i];
    return str;
  }
  std::cerr << "M=\n";
  WriteMatrix(std::cerr, M);
  std::cerr << "The recognition failed so coxdyn_matrix_to_string cannot work\n";
  throw TerminalException{1};
}


struct DiagramSelector {
  bool OnlySimplyLaced;
  bool OnlyLorentzianAdmissible;
  bool OnlySpherical;
};


template<typename T>
bool CheckIrreducibleDiagram(const MyMatrix<T>& M, DiagramSelector const& DS)
{
  int n = M.rows();
  T valInfinity = practical_infinity<T>();
  for (int i=0; i<n; i++)
    for (int j=i+1; j<n; j++) {
      T val = M(i,j);
      if (DS.OnlySimplyLaced)
        if (val != 2 && val != 3)
          return false;
      if (DS.OnlyLorentzianAdmissible)
        if (val != 2 && val != 3 && val != 4 && val != 6 && val != valInfinity)
          return false;
    }
  std::optional<IrrCoxDyn<T>> opt = RecognizeIrreducibleSphericalEuclideanDiagram(M);
  if (opt) {
    if (DS.OnlySpherical) {
      const IrrCoxDyn<T>& cd = *opt;
      return IsDiagramSpherical(cd);
    }
    return true;
  }
  return false;
}


template<typename T>
bool CheckDiagram(const MyMatrix<T>& M, DiagramSelector const& DS)
{
  std::vector<std::vector<size_t>> LConn = GetIrreducibleComponents(M);
  //  std::cerr << "CheckDiagram M=\n";
  //  WriteMatrix(std::cerr, M);
  std::vector<IrrCoxDyn<T>> l_cd;
  for (auto & eConn : LConn) {
    size_t dim_res=eConn.size();
    /*
    std::cerr << "eConn =";
    for (auto & val : eConn)
      std::cerr << " " << val;
      std::cerr << "\n";*/
    MyMatrix<T> Mres(dim_res, dim_res);
    for (size_t i=0; i<dim_res; i++)
      for (size_t j=0; j<dim_res; j++)
        Mres(i,j) = M(eConn[i], eConn[j]);
    bool test = CheckIrreducibleDiagram(Mres, DS);
    /*
    std::cerr << "  test=" << test << " eConn=";
    for (auto & v : eConn)
      std::cerr << v << " ";
    std::cerr << "Mres=\n";
    WriteMatrix(std::cerr, Mres);*/
    if (!test)
      return false;
  }
  return true;
}


template<typename T>
MyMatrix<T> ExtendMatrix(MyMatrix<T> const& M, MyVector<T> const& V)
{
  size_t len = M.rows();
  MyMatrix<T> Mret(len+1, len+1);
  for (size_t i=0; i<len; i++)
    for (size_t j=0; j<len; j++)
      Mret(i,j) = M(i,j);
  for (size_t i=0; i<len; i++) {
    Mret(len,i) = V(i);
    Mret(i,len) = V(i);
  }
  return Mret;
}


template<typename T>
MyMatrix<T> ExtendMatrixNorm(MyMatrix<T> const& M, MyVector<T> const& V, T const& norm)
{
  size_t len = M.rows();
  MyMatrix<T> Mret(len+1, len+1);
  for (size_t i=0; i<len; i++)
    for (size_t j=0; j<len; j++)
      Mret(i,j) = M(i,j);
  for (size_t i=0; i<len; i++) {
    Mret(len,i) = V(i);
    Mret(i,len) = V(i);
  }
  Mret(len,len) = norm;
  return Mret;
}


template<typename T>
std::vector<MyVector<T>> FindDiagramExtensions(const MyMatrix<T>& M, const DiagramSelector& DS)
{
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "FindDiagramExtensions, step 1\n";
#endif
  std::vector<std::vector<size_t>> LConn = GetIrreducibleComponents(M);
  std::set<MyVector<T>> SetExtensions;
  T val_comm = 2;
  T val_single_edge = 3;
  T val_four = 4;
  T val_six = 6;
  T val_inf = practical_infinity<T>(); // For supporting I1(infinity)
  std::vector<T> allowed_vals;
  if (DS.OnlyLorentzianAdmissible) {
    allowed_vals.push_back(val_single_edge);
    allowed_vals.push_back(val_four);
    allowed_vals.push_back(val_six);
    if (!DS.OnlySpherical)
      allowed_vals.push_back(val_inf);
  } else {
    for (T val=val_single_edge; val<128; val++)
      allowed_vals.push_back(val);
    if (!DS.OnlySpherical)
      allowed_vals.push_back(val_inf);
  }
  size_t dim = M.rows();
  std::vector<size_t> list_deg(dim);
  std::vector<size_t> list_n_higher(dim);
  std::vector<size_t> list_isolated;
  for (size_t i=0; i<dim; i++) {
    size_t n_adj = 0;
    size_t n_higher = 0;
    for (size_t j=0; j<dim; j++) {
      T val = M(i,j);
      if (i != j && val != val_comm)
        n_adj++;
      if (i != j && val != val_comm && val != val_single_edge)
        n_higher++;
    }
    list_deg[i] = n_adj;
    list_n_higher[i] = n_higher;
    if (n_adj == 0)
      list_isolated.push_back(i);
  }
  size_t n_isolated = list_isolated.size();
  std::vector<size_t> list_deg01_for_triple_vertex; // Part of En, Dn, tilde or not
  for (auto &eConn : LConn) {
    size_t dim_res=eConn.size();
    MyMatrix<T> Mres(dim_res, dim_res);
    for (size_t i=0; i<dim_res; i++)
      for (size_t j=0; j<dim_res; j++)
        Mres(i,j) = M(eConn[i], eConn[j]);
    std::optional<IrrCoxDyn<T>> opt = RecognizeIrreducibleSphericalEuclideanDiagram(Mres);
    if (!opt) {
      std::cerr << "The diagram should have been recognized\n";
      throw TerminalException{1};
    }
    const IrrCoxDyn<T>& cd = *opt;
    // The diagrams that can be part of triple points are:
    //   An, Dn, Bn
    // The diagrams that cannot be part of triple points are:
    //   En, tildeEn, tildeDn, tildeBn, tildeCn, F4, tildeF4, G2, tildeG2
    if (cd.type == "A" || cd.type == "D" || cd.type == "B") {
      for (auto & eVert : eConn) {
        if (list_deg[eVert] <= 1) // Case 0 corresponds to A1
          list_deg01_for_triple_vertex.push_back(eVert);
      }
    }
  }
  // Consider the case of adding unconnected vector
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "FindDiagramExtensions, step 2\n";
#endif
  MyVector<T> V_basic(dim);
  for (size_t i=0; i<dim; i++)
    V_basic(i) = val_comm;
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "FindDiagramExtensions, step 3\n";
#endif
  size_t n_diagram_considered = 0;
  size_t n_diagram_match = 0;
  auto test_vector_and_insert=[&](const MyVector<T>& V) -> void {
    MyMatrix<T> Mtest(dim+1,dim+1);
    for (size_t i=0; i<dim; i++)
      for (size_t j=0; j<dim; j++)
        Mtest(i,j) = M(i,j);
    for (size_t i=0; i<dim; i++) {
      Mtest(i,dim) = V(i);
      Mtest(dim,i) = V(i);
    }
    //    std::cerr << "V=" << StringVectorGAP(V) << "\n";
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
    std::cerr << "Mtest built\n";
#endif
    n_diagram_considered++;
    bool test = CheckDiagram(Mtest, DS);
    //    std::cerr << "test=" << test << "\n";
    if (test) {
      SetExtensions.insert(V);
      n_diagram_match++;
    }
  };
  test_vector_and_insert(V_basic); // Adding just an A1, always works.
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "FindDiagramExtensions, step 4\n";
#endif
  // Considering the case of just one edge
  if (DS.OnlyLorentzianAdmissible) {
    for (size_t i=0; i<dim; i++) {
      MyVector<T> V = V_basic;
      V(i) = val_single_edge;
      test_vector_and_insert(V);
      V(i) = val_four;
      test_vector_and_insert(V);
      V(i) = val_six; // This covers G2 and tilde{G2}
      test_vector_and_insert(V);
    }
    for (auto & eIsol : list_isolated) {
      MyVector<T> V = V_basic;
      if (!DS.OnlySpherical) {
        V(eIsol) = val_inf;
        test_vector_and_insert(V); // I1(infinity), always works.
      }
    }
  } else {
    for (size_t i=0; i<dim; i++) {
      // Here we have an arbitrary value
      for (auto & val : allowed_vals) {
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
        std::cerr << "i=" << i << " val=" << val << "\n";
#endif
        MyVector<T> V = V_basic;
        V(i) = val;
        test_vector_and_insert(V);
      }
    }
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "FindDiagramExtensions, step 5\n";
#endif
  // Considering the case of 2 edges
  if (DS.OnlyLorentzianAdmissible) {
    for (size_t i=0; i<dim; i++) {
      for (size_t j=i+1; j<dim; j++) {
        MyVector<T> V = V_basic;
        V(i) = val_single_edge;
        V(j) = val_single_edge;
        test_vector_and_insert(V);
      }
    }
    for (size_t i=0; i<dim; i++) {
      for (size_t j=0; j<dim; j++) {
        if (i != j) {
          MyVector<T> V = V_basic;
          V(i) = val_single_edge;
          V(j) = val_four;
          test_vector_and_insert(V);
        }
      }
    }
    if (!DS.OnlySpherical) { // Only tildeG2 and tilde{C2} are possible here
      for (size_t i=0; i<n_isolated; i++) {
        for (size_t j=0; j<n_isolated; j++) {
          if (i != j) {
            MyVector<T> V = V_basic;
            V(list_isolated[i]) = val_single_edge;
            V(list_isolated[j]) = val_six;
            test_vector_and_insert(V); // For tilde{G2}
          }
        }
      }
      for (size_t i=0; i<n_isolated; i++) {
        for (size_t j=i+1; j<n_isolated; j++) {
          MyVector<T> V = V_basic;
          V(list_isolated[i]) = val_four;
          V(list_isolated[j]) = val_four;
          test_vector_and_insert(V); // For tilde{C2}
        }
      }
    }
  } else {
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
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "FindDiagramExtensions, step 6\n";
#endif
  // Considering the case of 3 edges. Considering first part with single edges
  size_t n_cand_triple = list_deg01_for_triple_vertex.size();
  SetCppIterator SCI_A(n_cand_triple,3);
  for (auto & eV : SCI_A) {
    MyVector<T> V = V_basic;
    for (auto & eVal : eV)
      V(list_deg01_for_triple_vertex[eVal]) = val_single_edge;
    test_vector_and_insert(V);
  }
  if (!DS.OnlySpherical) { // Considering now the tilde{B3} cases
    SetCppIterator SCI_B(n_isolated,3);
    T val;
    for (auto & eV : SCI_B) {
      MyVector<T> V = V_basic;
      for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
          val = val_single_edge;
          if (i == j)
            val = val_four;
          V(list_isolated[eV[j]]) = val;
        }
        test_vector_and_insert(V);
      }
    }
  }
  // Considering the case of 4 edges. Only tilde{D4} is possible
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "FindDiagramExtensions, step 7\n";
#endif
  if (!DS.OnlySpherical) { // Only tildeD4 is feasible, and it is not euclidean
    SetCppIterator SCI_B(n_isolated,4);
    for (auto & eV : SCI_B) {
      MyVector<T> V = V_basic;
      for (auto & eVal : eV)
        V(list_isolated[eVal]) = val_single_edge;
      test_vector_and_insert(V);
    }
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "FindDiagramExtensions, step 8\n";
#endif
  std::vector<MyVector<T>> ListExtensions;
  for (auto &eEnt : SetExtensions)
    ListExtensions.push_back(eEnt);
  std::cerr << "Stats : |ListExtensions|=" << ListExtensions.size() << " n_diagram_considered=" << n_diagram_considered << " n_diagram_match=" << n_diagram_match << "\n";
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
  T cossqr_val3 = T(1) / T(4);
  T cossqr_val4 = T(1) / T(2);
  T cossqr_val6 = T(3) / T(4);
  T cossqr_valInf = 1;
  T pr_inf = practical_infinity<T>();
  auto get_cossqr_scal=[&](int i, int j) -> std::pair<T,T> {
    T scal12 = get_scal(i, j);
    if (i == j) {
      return {scal12, scal12};
    } else {
      if (scal12 == 0)
        return {2,scal12};
      T scal11 = get_scal(i, i);
      T scal22 = get_scal(j, j);
      T quot = (scal12 * scal12) / (scal11 * scal22);
      if (quot == cossqr_val3)
        return {3,scal12};
      if (quot == cossqr_val4)
        return {4,scal12};
      if (quot == cossqr_val6)
        return {6,scal12};
      if (quot == cossqr_valInf)
        return {pr_inf,scal12};
      std::cerr << "i=" << i << " j=" << j << "\n";
      std::cerr << "scal12=" << scal12 << " scal11=" << scal11 << " scal22=" << scal22 << "\n";
      std::cerr << "l_root=\n";
      WriteMatrix(std::cerr, MatrixFromVectorFamily(l_root));
      std::cerr << "Failed to find matching entry quot=" << quot << "\n";
      throw TerminalException{1};
    }
  };
  size_t n_root=l_root.size();
  MyMatrix<T> CoxMat(n_root,n_root);
  MyMatrix<T> ScalMat(n_root,n_root);
  for (size_t i=0; i<n_root; i++)
    for (size_t j=0; j<n_root; j++) {
      std::pair<T,T> ep = get_cossqr_scal(i,j);
      CoxMat(i,j) = ep.first;
      ScalMat(i,j) = ep.second;
    }
  return {std::move(CoxMat), std::move(ScalMat)};
}









template<typename T>
struct Possible_Extension {
  MyVector<T> u_component;
  T res_norm;
  T e_norm;
};

template<typename T, typename Tint>
std::vector<Possible_Extension<T>> ComputePossibleExtensions(MyMatrix<T> const& G, std::vector<MyVector<Tint>> const& l_root, std::vector<T> const& l_norm, bool only_spherical)
{
  std::cerr << "------------------------------------ ComputePossibleExtension ---------------------------------\n";
  DiagramSelector DS;
  DS.OnlySimplyLaced = false;
  DS.OnlyLorentzianAdmissible = true;
  DS.OnlySpherical = only_spherical;
  //
  std::cerr << "ComputePossibleExtensions, step 1\n";
  //  std::cerr << "G=\n";
  //  WriteMatrixGAP(std::cerr, G);
  //  std::cerr << "\n";
  //  std::cerr << "l_root=\n";
  //  for (auto & e_root : l_root)
  //    std::cerr << "e_root=" << StringVectorGAP(e_root) << "\n";

  std::pair<MyMatrix<T>,MyMatrix<T>> ep = ComputeCoxeterMatrix(G, l_root);
  const MyMatrix<T> & CoxMat = ep.first;
  const MyMatrix<T> & ScalMat = ep.second;
  MyMatrix<T> ScalMatInv = Inverse(ScalMat);
  std::cerr << "ScalMat=\n"; WriteMatrix(std::cerr, ScalMat);
  std::cerr << "CoxMat=\n"; WriteMatrix(std::cerr, CoxMat);
  std::cerr << "Symbol of M=" << coxdyn_matrix_to_string(CoxMat) << "\n";
  int dim = G.rows();
  int n_root = l_root.size();
  std::cerr << "ComputePossibleExtensions, step 3\n";
  std::vector<MyVector<T>> l_vect = FindDiagramExtensions(CoxMat, DS);
  std::cerr << "|l_vect|=" << l_vect.size() << "\n";
  T val2 = 0;
  T val3 = T(1) / T(4);
  T val4 = T(1) / T(2);
  T val6 = T(3) / T(4);
  T valInfinity = 1;
  auto get_cos_square=[&](T val) -> T {
    if (val == 2)
      return val2;
    if (val == 3)
      return val3;
    if (val == 4)
      return val4;
    if (val == 6)
      return val6;
    if (val == practical_infinity<T>())
      return valInfinity;
    std::cerr << "Failed to find a matching entry val=" << val << "\n";
    throw TerminalException{1};
  };
  std::vector<Possible_Extension<T>> l_extensions;
  auto get_entry=[&](MyVector<T> const& e_vect, T const& e_norm) -> void {
    //    std::cerr << "---------------- e_norm=" << e_norm << " e_vect=" <<  StringVectorGAP( e_vect) << "\n";
    MyVector<T> l_scal(n_root);
    for (int i=0; i<n_root; i++) {
      T val = e_vect(i);
      T cos_square = get_cos_square(val);
      T scal_square = cos_square * CoxMat(i,i) * e_norm;
      std::optional<T> opt = UniversalSquareRoot(scal_square);
      //      std::cerr << "i=" << i << " cos_square=" << cos_square << " CoxMat(i,i)=" << CoxMat(i,i) << " e_norm=" << e_norm << " scal_square=" << scal_square << "\n";
      //      std::cerr << "i=" << i << " scal_square=" << scal_square << "\n";
      if (!opt) {
        //        std::cerr << "   Failed to match\n";
        return;
      }
      T scal = - *opt;
      //      std::cerr << "     scal=" << scal << "\n";
      l_scal(i) = scal;
    }
    std::cerr << "---------------- e_norm=" << e_norm << " e_vect=" << StringVectorGAP(e_vect) << "\n";
    //    std::cerr << "Scalar products found : l_scal = " << StringVectorGAP(l_scal) << "\n";
    /* So, we have computed l_scal(i) = alpha.dot.ui = u.dot.ui
       If u = sum wi u_i then w = G^{-1} l_scal
       eNorm = w.dot.w  is the Euclidean norm of u.
     */
    MyVector<T> w = ScalMatInv * l_scal;
    //    std::cerr << "w = " << StringVectorGAP(w) << "\n";
    T eNorm = l_scal.dot(w);
    T res_norm = e_norm - eNorm;
    std::cerr << "  eNorm=" << eNorm << " res_norm=" << res_norm << "\n";
    MyVector<T> u_component = ZeroVector<T>(dim);
    for (int i=0; i<n_root; i++)
      u_component += w(i) * UniversalVectorConversion<T,Tint>(l_root[i]);
#ifdef SANITY_CHECK
    if (res_norm == 0) {
      std::vector<MyVector<T>> l_root_tot;
      for (auto & ev_tint : l_root)
        l_root_tot.push_back(UniversalVectorConversion<T,Tint>(ev_tint));
      l_root_tot.push_back(u_component);
      MyMatrix<T> CoxMatExt = ComputeCoxeterMatrix(G, l_root_tot).first;
      MyMatrix<T> CoxMatBld = ExtendMatrixNorm(CoxMat, e_vect, e_norm);
      if (CoxMatExt != CoxMatBld) {
        std::cerr << "CoxMatExt=\n";
        WriteMatrix(std::cerr, CoxMatExt);
        std::cerr << "CoxMatBld=\n";
        WriteMatrix(std::cerr, CoxMatBld);
        std::cerr << "The matrices should be equal\n";
        throw TerminalException{1};
      }
      //      std::cerr << "Symbol of CoxMatExt=" << coxdyn_matrix_to_string(CoxMatExt) << "\n";
    }
#endif
    std::cerr << "  u_component="; WriteVectorGAP(std::cerr, u_component); std::cerr << "\n";
    l_extensions.push_back({u_component, res_norm, e_norm});
  };
  std::cerr << "ComputePossibleExtensions, step 5\n";
  size_t len = l_norm.size();
  Face status_norm(len);
  for (auto & e_vect : l_vect) {
    for (size_t i_norm=0; i_norm<len; i_norm++)
      status_norm[i_norm] = 1;
    for (int i=0; i<n_root; i++) {
      T norm = CoxMat(i,i);
      for (size_t i_norm=0; i_norm<len; i_norm++) {
        if (e_vect(i) == 3) {
          if (l_norm[i_norm] != norm)
            status_norm[i_norm] = 0;
        }
        if (e_vect(i) == 4) {
          if (2 * l_norm[i_norm] != norm && l_norm[i_norm] != 2 * norm)
            status_norm[i_norm] = 0;
        }
        if (e_vect(i) == 6) {
          if (3 * l_norm[i_norm] != norm && l_norm[i_norm] != 3 * norm)
            status_norm[i_norm] = 0;
        }
      }
    }
    for (size_t i_norm=0; i_norm<len; i_norm++)
      if (status_norm[i_norm] == 1)
        get_entry(e_vect, l_norm[i_norm]);
  }
  std::cerr << "ComputePossibleExtensions, step 6 |l_extensions|=" << l_extensions.size() << "\n";
  return l_extensions;
}

#endif

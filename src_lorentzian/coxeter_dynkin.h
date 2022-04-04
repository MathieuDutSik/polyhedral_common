#ifndef SRC_LORENTZIAN_COXETER_DYNKIN_H_
#define SRC_LORENTZIAN_COXETER_DYNKIN_H_

#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalBasic.h"
#include "GRAPH_GraphicalFunctions.h"
#include <limits>
#include <map>
#include <algorithm>
#include <utility>
#include <set>
#include <string>
#include <vector>

/*
  The Coxeter Dynkin diagram are build in the following way:
  --- The values for i != j are the exponent m such that (g_ig_j)^m = Id
  --- The values M(i,i) are not used.
 */

//#define DEBUG_COXETER_DYNKIN_COMBINATORICS

//#define CHECK_EFFICIENT_ENUMERATION

struct DiagramSelector {
  bool OnlySimplyLaced;
  bool OnlyLorentzianAdmissible;
  bool OnlySpherical;
};

template <typename T> struct IrrCoxDyn {
  std::string type;
  size_t dim;
  T param; // For In only
};

template <typename T> bool IsDiagramSpherical(IrrCoxDyn<T> const &cd) {
  std::string type = cd.type;
  if (type == "tildeA" || type == "tildeB" || type == "tildeC" ||
      type == "tildeD" || type == "tildeE" || type == "tildeF" ||
      type == "tildeG")
    return false;
  return true;
}

template <typename T> bool IsDiagramSimplyLaced(IrrCoxDyn<T> const &cd) {
  std::string type = cd.type;
  if (type == "A" || type == "D" || type == "E")
    return true;
  return false;
}

template <typename T>
bool IsDiagramIntegerLorentzianFeasible(IrrCoxDyn<T> const &cd) {
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

template <typename T> int GetNrVertices(IrrCoxDyn<T> const &cd) {
  std::string type = cd.type;
  int dim = cd.dim;
  if (type.size() > 5) {
    std::string s_res = type.substr(0, 5);
    if (s_res == "tilde")
      return dim + 1;
  }
  return dim;
}

template <typename T> std::string IrrCoxDyn_to_string(IrrCoxDyn<T> const &cd) {
  std::string type = cd.type;
  if (type == "A" || type == "B" || type == "C" || type == "D" || type == "E" ||
      type == "F" || type == "G" || type == "H")
    return cd.type + "_{" + std::to_string(cd.dim) + "}";
  if (type == "I") {
    if (cd.param == practical_infinity<T>())
      return "I_2(infinity)";
    return std::string("I_2(") + std::to_string(cd.param) + ")";
  }
  if (type == "tildeA" || type == "tildeB" || type == "tildeC" ||
      type == "tildeD" || type == "tildeE" || type == "tildeF" ||
      type == "tildeG") {
    std::string type_red = type.substr(5, 1);
    return std::string("\\tilde{") + type_red + "_{" + std::to_string(cd.dim) +
           "}}";
  }
  if (type == "tildeI")
    return std::string("\\tilde{I_{1}(infinity)}");
  std::cerr << "cd  type=" << cd.type << " dim=" << cd.dim
            << " param=" << cd.param << "\n";
  std::cerr << "Failed to matching entry. Maybe bug or non-conforming input\n";
  throw TerminalException{1};
}

template <typename T> IrrCoxDyn<T> string_to_IrrCoxDyn(std::string const &s) {
  std::string s_work = s;
  auto remove_char = [&](char ec) -> void {
    s_work.erase(remove(s_work.begin(), s_work.end(), ec), s_work.end());
  };
  remove_char('_');
  remove_char('\\');
  remove_char('{');
  remove_char('}');
  remove_char('(');
  remove_char(')');
  auto recognize = [&]() -> IrrCoxDyn<T> {
    std::vector<std::string> LS{
        "A",      "B",      "C",      "D",      "E",      "F",  "G",      "H",
        "tildeA", "tildeB", "tildeC", "tildeD", "tildeE", "I2", "tildeI1"};
    for (auto &eLS : LS) {
      size_t len1 = s_work.size();
      size_t len2 = eLS.size();
      if (len1 > len2) {
        std::string s_red = s_work.substr(0, len2);
        if (s_red == eLS) {
          std::string s_rem = s_work.substr(len2, len1 - len2);
          if (eLS == "tildeI1") {
            if (s_rem == "infinity")
              return IrrCoxDyn<T>{"tildeI", 1, practical_infinity<T>()};
          } else {
            if (eLS == "I2") {
              T val = ParseScalar<T>(s_rem);
              return IrrCoxDyn<T>{"I", 2, val};
            } else {
              size_t val = ParseScalar<size_t>(s_rem);
              return IrrCoxDyn<T>{eLS, val, 0};
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
    std::cerr << "Found matching to be type=" << cd.type << " dim=" << cd.dim
              << " param=" << cd.param << "\n";
    std::cerr << "Mapped string is s_map=" << s_map << "\n";
    throw TerminalException{1};
  }
  return cd;
}

template <typename T>
MyMatrix<T> Kernel_IrrCoxDyn_to_matrix(IrrCoxDyn<T> const &cd) {
  T val_comm = 2;
  T val_single_edge = 3;
  T val_four = 4; // Shows up in F4, Bn = Cn, tilde{Bn}, tilde{Cn}, tilde{F4}.
  T val_five = 5; // Shows up in H3, H4
  T val_six = 6;  // Shows up in G2
  std::string type = cd.type;
  int dim = cd.dim;
  int n_vert = GetNrVertices(cd);
  MyMatrix<T> M(n_vert, n_vert);
  for (int i = 0; i < n_vert; i++)
    for (int j = 0; j < n_vert; j++)
      M(i, j) = val_comm;
  size_t n_assign = 0;
  auto set_v = [&](int i, int j, T val) -> void {
    M(i, j) = val;
    M(j, i) = val;
    n_assign++;
  };
  if (type == "A") {
    for (int i = 1; i < dim; i++)
      set_v(i - 1, i, val_single_edge);
    if (dim ==
        1) // That one has zero edges. So, we need to increase so as to pass.
      n_assign++;
  }
  if (type == "B" || type == "C") {
    for (int i = 1; i < dim - 1; i++)
      set_v(i - 1, i, val_single_edge);
    set_v(dim - 2, dim - 1, val_four);
  }
  if (type == "D") {
    for (int i = 1; i < dim - 1; i++)
      set_v(i - 1, i, val_single_edge);
    set_v(1, dim - 1, val_single_edge);
  }
  auto set_from_triple = [&](int l1, int l2, int l3) -> void {
    auto set_line = [&](int l, int shift) -> void {
      for (int i = 0; i < l; i++) {
        if (i == 0)
          set_v(0, shift + i + 1, val_single_edge);
        else
          set_v(shift + i, shift + i + 1, val_single_edge);
      }
    };
    set_line(l1, 0);
    set_line(l2, l1);
    set_line(l3, l1 + l2);
  };
  if (type == "E") {
    if (dim == 6)
      set_from_triple(1, 2, 2);
    if (dim == 7)
      set_from_triple(1, 2, 3);
    if (dim == 8)
      set_from_triple(1, 2, 4);
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
    for (int i = 2; i < dim; i++)
      set_v(i - 1, i, val_single_edge);
    set_v(0, 1, val_five);
  }
  if (type == "I" && dim == 2) {
    if (cd.param < 7) {
      std::cerr << "For the I case the values 2, 3, 4, 5, 6 are covered by "
                   "A1+A1 , A2 , B2 , H2 and G2\n";
      throw TerminalException{1};
    }
    if (cd.param == practical_infinity<T>()) {
      std::cerr << "The I2(infinity) is not allowed\n";
      throw TerminalException{1};
    }
    set_v(0, 1, cd.param);
  }
  // Now the euclidean ones
  if (type == "tildeA") {
    for (int i = 1; i < n_vert; i++)
      set_v(i - 1, i, val_single_edge);
    set_v(0, n_vert - 1, val_single_edge);
  }
  if (type == "tildeB") {
    for (int i = 1; i < n_vert - 2; i++)
      set_v(i - 1, i, val_single_edge);
    set_v(0, n_vert - 2, val_four);
    set_v(n_vert - 4, n_vert - 1, val_single_edge);
  }
  if (type == "tildeC") {
    for (int i = 1; i < n_vert - 2; i++)
      set_v(i - 1, i, val_single_edge);
    set_v(0, n_vert - 2, val_four);
    set_v(n_vert - 3, n_vert - 1, val_four);
  }
  if (type == "tildeD") {
    for (int i = 1; i < n_vert - 2; i++)
      set_v(i - 1, i, val_single_edge);
    set_v(1, n_vert - 2, val_single_edge);
    set_v(n_vert - 4, n_vert - 1, val_single_edge);
  }
  if (type == "tildeE") {
    if (dim == 6)
      set_from_triple(2, 2, 2);
    if (dim == 7)
      set_from_triple(1, 3, 3);
    if (dim == 8)
      set_from_triple(1, 2, 5);
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
    std::cerr << "type=" << type << " dim=" << dim << " param=" << cd.param
              << "\n";
    throw TerminalException{1};
  }
  return M;
}

template <typename T>
std::optional<IrrCoxDyn<T>>
RecognizeIrreducibleSphericalEuclideanDiagram(const MyMatrix<T> &M) {
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram M=\n";
  WriteMatrix(std::cerr, M);
#endif
  T val_comm = 2;
  T val_single_edge = 3;
  T val_four = 4; // Shows up in F4, Bn = Cn, tilde{Bn}, tilde{Cn}, tilde{F4}.
  T val_five = 5; // Shows up in H3, H4
  T val_six = 6;  // Shows up in G2
  size_t n_vert = M.rows();
  if (n_vert == 1) // Case of A1
    return IrrCoxDyn<T>{"A", 1, 0};
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 1\n";
#endif
  std::vector<size_t> list_deg1, list_deg2, list_deg3, list_deg4, list_degN;
  std::vector<size_t> list_deg(n_vert, 0);
  size_t n_higher_edge = 0;
  GraphBitset eG(n_vert);
  std::map<T, size_t> multiplicity;
  std::vector<std::vector<size_t>> LLAdj;
  std::vector<T> list_isolated_adjacent(n_vert);
  for (size_t i = 0; i < n_vert; i++) {
    size_t n_adj = 0;
    std::vector<size_t> LAdj;
    T e_isolated_adjacent = std::numeric_limits<T>::max();
    for (size_t j = 0; j < n_vert; j++) {
      T val = M(i, j);
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
  if (list_deg4.size() >
      0) { // Vertex of degree 4 can occur for \tilde{D4} only
    // We are now in Euclidean
    if (n_vert != 5) // It has to be \tilde{D4}.
      return {};
    if (list_deg4.size() != 1 || list_deg1.size() != 4 ||
        list_deg2.size() != 0 || list_deg3.size() != 0)
      return {};
    size_t i_4 = list_deg4[0];
    for (size_t i = 0; i < n_vert; i++)
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
    const std::vector<size_t> &eCycle = ListCycles[0];
    if (eCycle.size() != n_vert)
      return {};
    if (list_deg2.size() != n_vert)
      return {};
    if (n_higher_edge != 0)
      return {};
    return IrrCoxDyn<T>{"tildeA", n_vert - 1,
                        0}; // Only tilde{An} is left as possibility
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 5\n";
#endif
  // Now it is a tree
  if (list_deg1.size() == 2 && list_deg2.size() == n_vert - 2 &&
      n_higher_edge == 0)
    return IrrCoxDyn<T>{"A", n_vert, 0}; // Only An is possible so ok.
  // An and tilde{An} have been covered
  if (list_deg3.size() > 2)
    return {}; // No possibility for spherical and euclidean
  if (list_deg3.size() == 2) {
    // We are now in Euclidean case
    if (n_higher_edge != 0)
      return {};
    for (auto &ePt : list_deg3) {
      std::vector<size_t> const &LAdj = LLAdj[ePt];
      size_t n_deg1 = 0;
      for (auto &fPt : LAdj)
        if (list_deg[fPt] == 1)
          n_deg1++;
      if (n_deg1 != 2)
        return {};
    }
    return IrrCoxDyn<T>{"tildeD", n_vert - 1, 0}; // Only tilde{Dn} is possible
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr
      << "RecognizeIrreducibleSphericalEuclideanDiagram, step 6 |list_deg3|="
      << list_deg3.size() << "\n";
#endif
  if (list_deg3.size() == 0) { // We are in a single path.
    if (multiplicity[val_four] == 2) {
      if (n_higher_edge != 2) // There are some other higher edge, excluded
        return {};
      // Only tilde{Cn} is feasible and it is Euclidean
      for (auto &eVert : list_deg1)
        if (list_isolated_adjacent[eVert] != val_four)
          return {};
      return IrrCoxDyn<T>{"tildeC", n_vert - 1, 0}; // This is tilde{Cn}
    }
    if (multiplicity[val_four] ==
        1) { // Possibilities: Bn=Cn, F4, tilde{F4}, and B2 are possible
      if (n_vert == 2) {
        T param = 4;
        return IrrCoxDyn<T>{"B", 2, 0}; // It is I2(4)
      }
      if (n_higher_edge != 1)
        return {}; // There are other edges, excluded.
      size_t n_sing = 0;
      size_t n_four = 0;
      for (auto &eVert : list_deg1) {
        if (list_isolated_adjacent[eVert] == val_single_edge)
          n_sing++;
        if (list_isolated_adjacent[eVert] == val_four)
          n_four++;
      }
      if (n_sing == 2) {
        if (n_vert == 4) {
          return IrrCoxDyn<T>{"F", 4, 0}; // Only F4 is possible
        }
        if (n_vert == 5) { // Only tilde{F4} is possible. So conclude from that
          return IrrCoxDyn<T>{"tildeF", 4, 0};
        }
      }
      if (n_four == 1 && n_sing == 1) {
        // Only possibility is to have 4 at one extremity. This is Bn = Cn
        return IrrCoxDyn<T>{"B", n_vert, 0};
      }
      // No other possibilities
      return {};
    }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
    std::cerr << "multiplicity[val_five]=" << multiplicity[val_five] << "\n";
#endif
    if (multiplicity[val_five] == 1) { // Looking for H2, H3, H4
      if (n_vert == 2)
        return IrrCoxDyn<T>{"H", 2, 0}; // It is H2
      if (n_vert > 5)
        return {}; // No possibility
      size_t n_sing = 0;
      size_t n_five = 0;
      for (auto &eVert : list_deg1) {
        if (list_isolated_adjacent[eVert] == val_single_edge)
          n_sing++;
        if (list_isolated_adjacent[eVert] == val_five)
          n_five++;
      }
      if (n_sing == 1 &&
          n_five == 1) { // It is H3 or H4 depending on the dimension
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
      T param = M(0, 1);
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
    std::cerr << "We should have just one vertex of degree 3\n";
    throw TerminalException{1};
  }
  size_t eCent = list_deg3[0];
  // Now just one vertex of degree 3.
  if (multiplicity[val_four] == 1) { // Possibility tilde{Bn}
    std::vector<size_t> const &LAdj = LLAdj[eCent];
    size_t n_sing = 0;
    size_t n_sing_simple = 0;
    for (auto &eAdj : LAdj)
      if (list_deg[eAdj] == 1) {
        n_sing++;
        if (list_isolated_adjacent[eAdj] == val_single_edge)
          n_sing_simple++;
      }
    if (n_sing == 3) { // All neighbors are single. Only one possibility
      if (n_sing_simple != 2)
        return {};
      return IrrCoxDyn<T>{"tildeB", 3, 0}; // It is tilde{B3}
    }
    if (n_sing != 2 || n_sing_simple != 2)
      return {};
    bool has_edge_four = false;
    for (auto &eVert : list_deg1)
      if (list_isolated_adjacent[eVert] == val_four)
        has_edge_four = true;
    if (has_edge_four) {
      return IrrCoxDyn<T>{"tildeB", n_vert - 1, 0}; // It is tilde{Bn}
    }
    return {};
  }
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "RecognizeIrreducibleSphericalEuclideanDiagram, step 8\n";
#endif
  if (n_higher_edge != 0)
    return {};
  auto get_length = [&](size_t val1, size_t val2) -> size_t {
    size_t len = 1;
    size_t iter = 0;
    while (true) {
      iter++;
      std::vector<size_t> const &LVal = LLAdj[val2];
      if (LVal.size() == 1)
        break;
      size_t NewPt = -1;
      for (auto &eVal : LVal)
        if (eVal != val1)
          NewPt = eVal;
      val1 = val2;
      val2 = NewPt;
      len++;
    }
    return len;
  };
  std::map<size_t, size_t> map_len;
  for (auto &eAdj : LLAdj[eCent]) {
    size_t len = get_length(eCent, eAdj);
    map_len[len]++;
  }
  if (map_len[1] == 3) // It is D4
    return IrrCoxDyn<T>{"D", 4, 0};
  if (map_len[1] == 2) // It is Dn
    return IrrCoxDyn<T>{"D", n_vert, 0};
  if (map_len[1] == 1 && map_len[2] == 2) // It is E6
    return IrrCoxDyn<T>{"E", 6, 0};
  if (map_len[1] == 1 && map_len[2] == 1 && map_len[3] == 1) // It is E7
    return IrrCoxDyn<T>{"E", 7, 0};
  if (map_len[1] == 1 && map_len[2] == 1 && map_len[4] == 1) // It is E8
    return IrrCoxDyn<T>{"E", 8, 0};
  // In spherical, no other possibilities left
  if (map_len[2] == 3) // It is tilde{E6}
    return IrrCoxDyn<T>{"tildeE", 6, 0};
  if (map_len[1] == 1 && map_len[3] == 2) // It is tilde{E7}
    return IrrCoxDyn<T>{"tildeE", 7, 0};
  if (map_len[1] == 1 && map_len[2] == 1 && map_len[5] == 1) // It is tilde{E8}
    return IrrCoxDyn<T>{"tildeE", 8, 0};
  return {}; // No other possibilities left
}

template <typename T>
bool CheckIrreducibleDiagram(const MyMatrix<T> &M, DiagramSelector const &DS) {
  int n = M.rows();
  T valInfinity = practical_infinity<T>();
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++) {
      T val = M(i, j);
      if (DS.OnlySimplyLaced)
        if (val != 2 && val != 3)
          return false;
      if (DS.OnlyLorentzianAdmissible)
        if (val != 2 && val != 3 && val != 4 && val != 6 && val != valInfinity)
          return false;
    }
  std::optional<IrrCoxDyn<T>> opt =
      RecognizeIrreducibleSphericalEuclideanDiagram(M);
  if (opt) {
    if (DS.OnlySpherical) {
      const IrrCoxDyn<T> &cd = *opt;
      return IsDiagramSpherical(cd);
    }
    return true;
  }
  return false;
}

template <typename T>
std::vector<std::vector<size_t>>
GetIrreducibleComponents(const MyMatrix<T> &M) {
  T val_comm = 2;
  size_t dim = M.rows();
  GraphBitset eG(dim);
  //  std::cerr << "LEdge =";
  for (size_t i = 0; i < dim; i++) {
    for (size_t j = i + 1; j < dim; j++) {
      if (M(i, j) != val_comm) {
        eG.AddAdjacent(i, j);
        eG.AddAdjacent(j, i);
        //        std::cerr << " [" << i << "/" << j << "]";
      }
    }
  }
  //  std::cerr << "\n";
  return ConnectedComponents_set(eG);
}

template <typename T>
bool CheckDiagram(const MyMatrix<T> &M, DiagramSelector const &DS) {
  std::vector<std::vector<size_t>> LConn = GetIrreducibleComponents(M);
  //  std::cerr << "CheckDiagram M=\n";
  //  WriteMatrix(std::cerr, M);
  std::vector<IrrCoxDyn<T>> l_cd;
  for (auto &eConn : LConn) {
    size_t dim_res = eConn.size();
    /*
    std::cerr << "eConn =";
    for (auto & val : eConn)
      std::cerr << " " << val;
      std::cerr << "\n";*/
    MyMatrix<T> Mres(dim_res, dim_res);
    for (size_t i = 0; i < dim_res; i++)
      for (size_t j = 0; j < dim_res; j++)
        Mres(i, j) = M(eConn[i], eConn[j]);
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

template <typename T>
std::vector<MyVector<T>>
FindDiagramExtensions_Efficient(const MyMatrix<T> &M,
                                const DiagramSelector &DS) {
  if (!DS.OnlyLorentzianAdmissible) {
    std::cerr << "We work only with Lorentzian admissible lattices\n";
    throw TerminalException{1};
  }
  std::vector<std::vector<size_t>> LConn = GetIrreducibleComponents(M);
  std::vector<MyVector<T>> ListExtensions;
#ifdef CHECK_EFFICIENT_ENUMERATION
  std::set<MyVector<T>> SetExtensions;
#endif
  T val_comm = 2;
  T val_single_edge = 3;
  T val_four = 4;
  T val_six = 6;
  T val_inf = practical_infinity<T>(); // For supporting I1(infinity)
  size_t n_vert = M.rows();
  std::vector<size_t> list_deg(n_vert);
  std::vector<size_t> list_n_higher(n_vert);
  std::vector<size_t> list_isolated;
  std::vector<T> list_isolated_adjacent_value(n_vert);
  std::vector<size_t> list_isolated_adjacent_index(n_vert);
  std::vector<std::vector<size_t>> LLAdj;
  for (size_t i = 0; i < n_vert; i++) {
    size_t n_adj = 0;
    size_t n_higher = 0;
    std::vector<size_t> LAdj;
    for (size_t j = 0; j < n_vert; j++) {
      T val = M(i, j);
      if (i != j && val != val_comm) {
        n_adj++;
        LAdj.push_back(j);
        list_isolated_adjacent_value[i] = val;
        list_isolated_adjacent_index[i] = j;
      }
      if (i != j && val != val_comm && val != val_single_edge)
        n_higher++;
    }
    list_deg[i] = n_adj;
    list_n_higher[i] = n_higher;
    LLAdj.push_back(LAdj);
    if (n_adj == 0)
      list_isolated.push_back(i);
  }
  size_t n_isolated = list_isolated.size();
  std::vector<std::vector<size_t>> list_extremal_A2;
  std::vector<std::vector<size_t>> list_extremal_A3;
  std::vector<std::vector<size_t>> list_extremal_AN;
  std::vector<size_t> list_vert_A2;
  std::vector<size_t> list_vert_G2;
  std::vector<size_t> list_middle_A3;
  std::vector<size_t> list_ends_A3;
  std::vector<size_t> list_ends_A4;
  std::vector<size_t> list_ends_A5;
  std::vector<size_t> list_ends_F4;
  std::vector<size_t> list_cent_D4;
  std::vector<size_t> list_expand_Bn; // For B2 this is the two vertices, for Bn
                                      // (n > 2) this is the expanding vertex
  std::vector<size_t>
      list_expand_m1_Bn; // For B2 this is the two vertices, for Bn (n > 2) this
                         // is the expanding vertex
  std::vector<size_t>
      list_non_expand_Bn; // Only for n > 2. This is the vertex adjacent with
                          // weight 4 which cannot be extended to B(n+1)
  std::vector<size_t> list_expand_Dn; // For D4 this is the 3 vertices, For Dn
                                      // (n > 4) this is the expanding vertex
  std::vector<size_t>
      list_expand_m1_Dn; // For D4 this is the 3 vertices, For Dn (n > 4) this
                         // is the expanding vertex
  std::vector<size_t>
      list_non_expand_Dn; // Only for n > 4, this is the two vertices which
                          // cannot be expanded to D(n+1)
  std::vector<size_t> VertToConn(
      n_vert); // Mapping from vertices to connected component. Useful for sums
               // like Ak + A3
  std::vector<size_t> VertToLocDim(
      n_vert); // Mapping from vertices to connected component. Useful for sums
               // like Ak + A3
  std::vector<size_t>
      list_extm1_AN; // For An the vertices from which we can expand to D(n+1)
  std::vector<size_t> list_extm2_AN; // Inspired by above, useful o get E6, E7,
                                     // E8, tilde{{E8} from A5, A6, A7 and A8.
  std::vector<size_t>
      list_extm3_AN; // Same as above, only useful to get tilde{E7} from A7
  std::vector<size_t>
      list_dist2_extrem_E6; // The two vertices of degree 1 of E6 diagram
  std::vector<size_t>
      list_dist1_extrem_E6; // The single vertices of degree 1 at distance 1
                            // from the vertex of E6 diagram
  std::vector<size_t> list_dist3_extrem_E7; // The single vertex of degree 1 at
                                            // distance 3 of E7 diagram
  std::vector<size_t> list_dist2_extrem_E7; // The single vertex of degree 1 at
                                            // distance 2 of E7 diagram
  std::vector<size_t> list_dist4_extrem_E8; // The single vertex of degree 1 at
                                            // distance 4 of E8 diagram
  size_t iConn = 0;
  for (auto &eConn : LConn) {
    size_t dim_res = eConn.size();
    for (auto &eVert : eConn) {
      VertToConn[eVert] = iConn;
      VertToLocDim[eVert] = dim_res;
    }
    iConn++;
    MyMatrix<T> Mres(dim_res, dim_res);
    for (size_t i = 0; i < dim_res; i++)
      for (size_t j = 0; j < dim_res; j++)
        Mres(i, j) = M(eConn[i], eConn[j]);
    std::optional<IrrCoxDyn<T>> opt =
        RecognizeIrreducibleSphericalEuclideanDiagram(Mres);
    if (!opt) {
      std::cerr
          << "The recognition did not get something Euclidean or Spherical\n";
      throw TerminalException{1};
    }
    const IrrCoxDyn<T> &cd = *opt;
    if (cd.type == "A" && cd.dim == 2) {
      list_extremal_A2.push_back(eConn);
      list_extremal_AN.push_back(eConn);
      for (auto &eVert : eConn)
        list_vert_A2.push_back(eVert);
    }
    if (cd.type == "G" && cd.dim == 2) {
      for (auto &eVert : eConn)
        list_vert_G2.push_back(eVert);
    }
    if (cd.type == "F" && cd.dim == 4) {
      for (auto &eVert : eConn)
        if (list_deg[eVert] == 1)
          list_ends_F4.push_back(eVert);
    }
    if (cd.type == "A" && 2 < cd.dim) {
      std::vector<size_t> Lext;
      for (auto &eVert : eConn)
        if (list_deg[eVert] == 1)
          Lext.push_back(eVert);
      if (cd.dim == 3)
        list_extremal_A3.push_back(Lext);
      if (cd.dim == 4) {
        for (auto &eVert : Lext)
          list_ends_A4.push_back(eVert);
      }
      if (cd.dim == 5) {
        for (auto &eVert : Lext)
          list_ends_A5.push_back(eVert);
      }
      list_extremal_AN.push_back(Lext);
      std::vector<size_t> Lextm1;
      bool IsFirst = true;
      for (auto &eVert : Lext) {
        if (IsFirst || cd.dim > 3) {
          size_t nVert = list_isolated_adjacent_index[eVert];
          Lextm1.push_back(nVert);
          list_extm1_AN.push_back(nVert);
        }
        IsFirst = false;
      }
      if (cd.dim > 4) {
        std::vector<size_t> Lextm2;
        bool IsFirst = true;
        for (auto &eVert : Lextm1) {
          for (auto &fVert : LLAdj[eVert]) {
            if (PositionVect(Lext, fVert) == -1) {
              if (IsFirst || cd.dim > 5) {
                Lextm2.push_back(fVert);
                list_extm2_AN.push_back(fVert);
              }
              IsFirst = false;
            }
          }
        }
        if (cd.dim == 7) { // only that case makes sense for A7 and its
                           // extension to tilde{E7}
          bool IsFirst = true;
          for (auto &eVert : Lextm2) {
            for (auto &fVert : LLAdj[eVert]) {
              if (PositionVect(Lextm1, fVert) == -1) {
                if (IsFirst) {
                  list_extm3_AN.push_back(fVert);
                }
                IsFirst = false;
              }
            }
          }
        }
      }
    }
    if (cd.type == "A" && cd.dim == 3) {
      for (auto &eVert : eConn) {
        if (list_deg[eVert] == 2)
          list_middle_A3.push_back(eVert);
        if (list_deg[eVert] == 1)
          list_ends_A3.push_back(eVert);
      }
    }
    if (cd.type == "B") {
      if (cd.dim == 2) {
        for (auto &eVert : eConn)
          list_expand_Bn.push_back(eVert);
      } else {
        for (auto &eVert : eConn) {
          if (list_deg[eVert] == 1) {
            if (list_isolated_adjacent_value[eVert] == val_single_edge) {
              list_expand_Bn.push_back(eVert);
              list_expand_m1_Bn.push_back(list_isolated_adjacent_index[eVert]);
            }
            if (list_isolated_adjacent_value[eVert] == val_four)
              list_non_expand_Bn.push_back(eVert);
          }
        }
      }
    }
    if (cd.type == "D") {
      if (cd.dim == 4) {
        for (auto &eVert : eConn) {
          if (list_deg[eVert] == 1)
            list_expand_Dn.push_back(eVert);
          if (list_deg[eVert] == 3)
            list_cent_D4.push_back(eVert);
        }
      } else {
        auto f = [&]() -> void {
          for (auto &eVert : eConn) {
            if (list_deg[eVert] == 3) {
              for (auto &fVert : eConn) {
                if (list_deg[fVert] == 1) {
                  if (list_isolated_adjacent_index[fVert] != eVert) {
                    list_expand_Dn.push_back(fVert);
                    list_expand_m1_Dn.push_back(
                        list_isolated_adjacent_index[fVert]);
                  }
                  if (list_isolated_adjacent_index[fVert] == eVert)
                    list_non_expand_Dn.push_back(fVert);
                }
              }
              return;
            }
          }
          std::cerr << "We should never reach that stage\n";
          throw TerminalException{1};
        };
        f();
      }
    }
    if (cd.type == "E") {
      auto get_length = [&](size_t val1,
                            size_t val2) -> std::pair<size_t, size_t> {
        size_t len = 1;
        size_t iter = 0;
        while (true) {
          iter++;
          std::vector<size_t> const &LVal = LLAdj[val2];
          if (LVal.size() == 1)
            break;
          size_t NewPt = -1;
          for (auto &eVal : LVal)
            if (eVal != val1)
              NewPt = eVal;
          val1 = val2;
          val2 = NewPt;
          len++;
        }
        return {len, val2};
      };
      for (auto &eVert : eConn) {
        if (list_deg[eVert] == 3) {
          for (auto &fVert : LLAdj[eVert]) {
            std::pair<size_t, size_t> ep = get_length(eVert, fVert);
            if (cd.dim == 6) {
              if (ep.first == 2)
                list_dist2_extrem_E6.push_back(ep.second);
              if (ep.first == 1)
                list_dist1_extrem_E6.push_back(ep.second);
            }
            if (cd.dim == 7) {
              if (ep.first == 3)
                list_dist3_extrem_E7.push_back(ep.second);
              if (ep.first == 2)
                list_dist2_extrem_E7.push_back(ep.second);
            }
            if (cd.dim == 8) {
              if (ep.first == 4)
                list_dist4_extrem_E8.push_back(ep.second);
            }
          }
        }
      }
    }
  }
  size_t n_A2 = list_extremal_A2.size();
  size_t n_A3 = list_extremal_A3.size();
  size_t n_AN = list_extremal_AN.size();
  // Consider the case of adding unconnected vector
  MyVector<T> V_basic(n_vert);
  for (size_t i = 0; i < n_vert; i++)
    V_basic(i) = val_comm;
  auto test_vector_and_insert = [&](const MyVector<T> &V) -> void {
#ifdef CHECK_EFFICIENT_ENUMERATION
    MyMatrix<T> Mtest(n_vert + 1, n_vert + 1);
    for (size_t i = 0; i < n_vert; i++)
      for (size_t j = 0; j < n_vert; j++)
        Mtest(i, j) = M(i, j);
    for (size_t i = 0; i < n_vert; i++) {
      Mtest(i, n_vert) = V(i);
      Mtest(n_vert, i) = V(i);
    }
    bool test = CheckDiagram(Mtest, DS);
    if (!test) {
      std::cerr << "V=" << StringVectorGAP(V) << "\n";
      std::cerr << "The proposed diagram extension is not adequate. We want to "
                   "avoid that\n";
      throw TerminalException{1};
    }
    if (SetExtensions.count(V) != 0) {
      std::cerr << "V=" << StringVectorGAP(V) << "\n";
      std::cerr << "The diagram is already present. We want to avoid that\n";
      throw TerminalException{1};
    }
    SetExtensions.insert(V);
#endif
    ListExtensions.push_back(V);
  };
  //
  // Case of 0 edge. Always valid
  //
  test_vector_and_insert(V_basic); // Adding just an A1, always works.
  //
  // Considering the case of just one edge
  //
  auto f_single = [&](size_t v) -> void {
    MyVector<T> V = V_basic;
    V(v) = val_single_edge;
    test_vector_and_insert(V);
  };
  // An obtained from A(n-1)
  for (auto &Lext : list_extremal_AN) {
    for (auto &v : Lext)
      f_single(v);
  }
  // A2 obtained from A1
  for (auto &v : list_isolated)
    f_single(v);
  // Bn obtained from A(n-1)
  for (auto &Lext : list_extremal_AN) {
    for (auto &v : Lext) {
      MyVector<T> V = V_basic;
      V(v) = val_four;
      test_vector_and_insert(V);
    }
  }
  // B2 obtained from A1
  for (auto &v : list_isolated) {
    MyVector<T> V = V_basic;
    V(v) = val_four;
    test_vector_and_insert(V);
  }
  // Bn obtained from B(n-1)
  for (auto &v : list_expand_Bn)
    f_single(v);
  // Dn from D(n-1)
  for (auto &v : list_expand_Dn)
    f_single(v);
  // Dn from A(n-1)
  for (auto &v : list_extm1_AN)
    f_single(v);
  // E6 obtained from A5
  for (auto &v : list_extm2_AN) {
    if (VertToLocDim[v] == 5)
      f_single(v);
  }
  // E6 obtained from D5
  for (auto &v : list_non_expand_Dn) {
    if (VertToLocDim[v] == 5)
      f_single(v);
  }
  // E7 obtained from A6
  for (auto &v : list_extm2_AN) {
    if (VertToLocDim[v] == 6)
      f_single(v);
  }
  // E7 obtained from D6
  for (auto &v : list_non_expand_Dn) {
    if (VertToLocDim[v] == 6)
      f_single(v);
  }
  // E7 obtained from E6
  for (auto &v : list_dist2_extrem_E6)
    f_single(v);
  // E8 from A7
  for (auto &v : list_extm2_AN) {
    if (VertToLocDim[v] == 7)
      f_single(v);
  }
  // E8 from D7
  for (auto &v : list_non_expand_Dn) {
    if (VertToLocDim[v] == 7)
      f_single(v);
  }
  // E8 from E7
  for (auto &v : list_dist3_extrem_E7)
    f_single(v);
  // G2 from A1
  for (auto &v : list_isolated) {
    MyVector<T> V = V_basic;
    V(v) = val_six;
    test_vector_and_insert(V);
  }
  // F4 from B3
  for (auto &v : list_non_expand_Bn) {
    if (VertToLocDim[v] == 3)
      f_single(v);
  }
  if (!DS.OnlySpherical) {
    // tilde{Bn} from Dn
    for (auto &v : list_expand_Dn) {
      MyVector<T> V = V_basic;
      V(v) = val_four;
      test_vector_and_insert(V);
    }
    // tilde{Bn} from Bn
    for (auto &v : list_expand_m1_Bn)
      f_single(v);
    // tilde{B3} from A3
    for (auto &v : list_middle_A3) {
      MyVector<T> V = V_basic;
      V(v) = val_four;
      test_vector_and_insert(V);
    }
    // tilde{Cn} from Bn
    for (auto &v : list_expand_Bn) {
      MyVector<T> V = V_basic;
      V(v) = val_four;
      test_vector_and_insert(V);
    }
    // tilde{Dn} from Dn
    for (auto &v : list_expand_m1_Dn)
      f_single(v);
    // tilde{D4} from D4
    for (auto &v : list_cent_D4)
      f_single(v);
    // tilde{E6} from E6
    for (auto &v : list_dist1_extrem_E6)
      f_single(v);
    // tilde{E7} from A7
    for (auto &v : list_extm3_AN)
      f_single(v);
    // tilde{E7} from E7
    for (auto &v : list_dist2_extrem_E7)
      f_single(v);
    // tilde{E8} from D8
    for (auto &v : list_non_expand_Dn)
      if (VertToLocDim[v] == 8)
        f_single(v);
    // tilde{E8} from A8
    for (auto &v : list_extm2_AN)
      if (VertToLocDim[v] == 8)
        f_single(v);
    // tilde{E8} from E8
    for (auto &v : list_dist4_extrem_E8)
      f_single(v);
    // tilde{G2} for G2
    for (auto &v : list_vert_G2)
      f_single(v);
    // tilde{G2} for A2
    for (auto &v : list_vert_A2) {
      MyVector<T> V = V_basic;
      V(v) = val_six;
      test_vector_and_insert(V);
    }
    // tilde{F4} from F4
    for (auto &v : list_ends_F4)
      f_single(v);
    // tilde{F4} from B4
    for (auto &v : list_non_expand_Bn)
      if (VertToLocDim[v] == 4)
        f_single(v);
    // I1(infinity) from A1
    for (auto &v : list_isolated) {
      MyVector<T> V = V_basic;
      V(v) = val_inf;
      test_vector_and_insert(V);
    }
  }
  //
  // Considering the case of 2 edges
  //
  auto f_pair = [&](size_t const &v1, size_t const &v2) -> void {
#ifdef CHECK_EFFICIENT_ENUMERATION
    if (v1 == v2) {
      std::cerr << "We should have v1 != v2\n";
      throw TerminalException{1};
    }
#endif
    MyVector<T> V = V_basic;
    V(v1) = val_single_edge;
    V(v2) = val_single_edge;
    test_vector_and_insert(V);
  };
  // An formed from Ak + Al with k+l = n-1 , k >= 2 , l >= 2
  SetCppIterator SCI_Ak_Al(n_AN, 2);
  for (auto &eV : SCI_Ak_Al) {
    for (auto &v1 : list_extremal_AN[eV[0]]) {
      for (auto &v2 : list_extremal_AN[eV[1]])
        f_pair(v1, v2);
    }
  }
  // An formed from A(n-2) + A1
  for (auto &v1 : list_isolated) {
    for (auto &LTerm : list_extremal_AN) {
      for (auto &v2 : LTerm)
        f_pair(v1, v2);
    }
  }
  // A3 formed from A1+A1
  SetCppIterator SCI_A3(n_isolated, 2);
  for (auto &eV : SCI_A3) {
    size_t v1 = list_isolated[eV[0]];
    size_t v2 = list_isolated[eV[1]];
    f_pair(v1, v2);
  }
  // B3 obtained from A1 + A1
  for (auto &v1 : list_isolated) {
    for (auto &v2 : list_isolated) {
      if (v1 != v2) {
        MyVector<T> V = V_basic;
        V(v1) = val_four;
        V(v2) = val_single_edge;
        test_vector_and_insert(V);
      }
    }
  }
  // Bn formed from A(n-2) + A1
  for (auto &v1 : list_isolated) {
    for (auto &LTerm : list_extremal_AN) {
      for (auto &v2 : LTerm) {
        MyVector<T> V = V_basic;
        V(v1) = val_four;
        V(v2) = val_single_edge;
        test_vector_and_insert(V);
      }
    }
  }
  // Bn formed from Bk + Al with k+l = n-1 , k>= 2 , l >= 2
  for (auto &v1 : list_expand_Bn) {
    for (auto &LTerm : list_extremal_AN) {
      for (auto &v2 : LTerm)
        f_pair(v1, v2);
    }
  }
  // Bn formed from B(n-2) + A1
  for (auto &v1 : list_expand_Bn) {
    for (auto &v2 : list_isolated)
      f_pair(v1, v2);
  }
  // Dn formed from A3 + Ak
  for (auto &v1 : list_middle_A3) {
    for (auto &LTerm : list_extremal_AN) {
      for (auto &v2 : LTerm) {
        if (VertToConn[v1] != VertToConn[v2])
          f_pair(v1, v2);
      }
    }
  }
  // D5 formed from A3 + A1
  for (auto &v1 : list_middle_A3) {
    for (auto &v2 : list_isolated) {
      f_pair(v1, v2);
    }
  }
  // Dn formed from Dk + Al with k+l = n-1 , k >= 4 , l >= 2
  for (auto &v1 : list_expand_Dn) {
    for (auto &LTerm : list_extremal_AN) {
      for (auto &v2 : LTerm)
        f_pair(v1, v2);
    }
  }
  // Dn formed from D(n-2) + A1
  for (auto &v1 : list_expand_Dn) {
    for (auto &v2 : list_isolated)
      f_pair(v1, v2);
  }
  // E6 formed as A4 + A1
  for (auto &v1 : list_isolated) {
    for (auto &v2 : list_extm1_AN) {
      if (VertToLocDim[v2] == 4)
        f_pair(v1, v2);
    }
  }
  // E7 formed as A5 + A1
  for (auto &v1 : list_isolated) {
    for (auto &v2 : list_extm1_AN) {
      if (VertToLocDim[v2] == 5)
        f_pair(v1, v2);
    }
  }
  // E7 formed as A4 + A2
  for (auto &v1 : list_vert_A2) {
    for (auto &v2 : list_extm1_AN) {
      if (VertToLocDim[v2] == 4)
        f_pair(v1, v2);
    }
  }
  // E7 formed as D5 + A1
  for (auto &v1 : list_isolated) {
    for (auto &v2 : list_non_expand_Dn) {
      if (VertToLocDim[v2] == 5)
        f_pair(v1, v2);
    }
  }
  // E8 from E6 + A1
  for (auto &v1 : list_dist2_extrem_E6) {
    for (auto &v2 : list_isolated)
      f_pair(v1, v2);
  }
  // E8 formed as A6 + A1
  for (auto &v1 : list_isolated) {
    for (auto &v2 : list_extm1_AN) {
      if (VertToLocDim[v2] == 6)
        f_pair(v1, v2);
    }
  }
  // E8 formed as A4 + A3
  for (auto &v1 : list_extm1_AN) {
    if (VertToLocDim[v1] == 4) {
      for (auto &v2 : list_ends_A3)
        f_pair(v1, v2);
    }
  }
  // E8 formed as D5 + A2
  for (auto &v1 : list_non_expand_Dn) {
    if (VertToLocDim[v1] == 5) {
      for (auto &v2 : list_vert_A2)
        f_pair(v1, v2);
    }
  }
  // F4 formed as A2 + A1
  for (auto &v1 : list_vert_A2) {
    for (auto &v2 : list_isolated) {
      MyVector<T> V = V_basic;
      V(v1) = val_four;
      V(v2) = val_single_edge;
      test_vector_and_insert(V);
    }
  }
  if (!DS.OnlySpherical) {
    // tilde{An} formed from An
    for (auto &Lext : list_extremal_AN)
      f_pair(Lext[0], Lext[1]);
    // tilde{Bn} obtained from A3 + B(n-3)
    for (auto &v1 : list_expand_Bn) {
      for (auto &v2 : list_middle_A3)
        f_pair(v1, v2);
    }
    // tilde{B4} from A3 + A1
    for (auto &v1 : list_isolated) {
      for (auto &v2 : list_middle_A3) {
        MyVector<T> V = V_basic;
        V(v1) = val_four;
        V(v2) = val_single_edge;
        test_vector_and_insert(V);
      }
    }
    // tilde{Bn} obtained from Dk + Bl with k+l = n , k >= 4 , l >= 2
    for (auto &v1 : list_expand_Bn) {
      for (auto &v2 : list_expand_Dn)
        f_pair(v1, v2);
    }
    // tilde{Bn} obtained from D(n-1) + A1
    for (auto &v1 : list_isolated) {
      for (auto &v2 : list_expand_Dn) {
        MyVector<T> V = V_basic;
        V(v1) = val_four;
        V(v2) = val_single_edge;
        test_vector_and_insert(V);
      }
    }
    // tilde{Cn} obtained from Bk + Bl with k+l = n , k >= 2 , l >= 2
    size_t n_expand_Bn = list_expand_Bn.size();
    SetCppIterator SCI_Bk_Bl(n_expand_Bn, 2);
    for (auto &eV : SCI_Bk_Bl) {
      size_t v1 = list_expand_Bn[eV[0]];
      size_t v2 = list_expand_Bn[eV[1]];
      if (VertToConn[v1] != VertToConn[v2])
        f_pair(v1, v2);
    }
    // tilde{Cn} obtained from C(n-1) + A1
    for (auto &v1 : list_isolated) {
      for (auto &v2 : list_expand_Bn) {
        MyVector<T> V = V_basic;
        V(v1) = val_four;
        V(v2) = val_single_edge;
        test_vector_and_insert(V);
      }
    }
    // tilde{C2} obtained from A1 + A1
    SetCppIterator SCI_A1_A1(n_isolated, 2);
    for (auto &eV : SCI_A1_A1) {
      size_t v1 = list_isolated[eV[0]];
      size_t v2 = list_isolated[eV[1]];
      MyVector<T> V = V_basic;
      V(v1) = val_four;
      V(v2) = val_four;
      test_vector_and_insert(V);
    }
    // tilde{G2} formed from A1+A1
    for (auto &v1 : list_isolated) {
      for (auto &v2 : list_isolated) {
        if (VertToConn[v1] != VertToConn[v2]) {
          MyVector<T> V = V_basic;
          V(v1) = val_single_edge;
          V(v2) = val_six;
          test_vector_and_insert(V);
        }
      }
    }
    // tilde{D6} obtained from A3 + A3
    size_t n_middle_A3 = list_middle_A3.size();
    SetCppIterator SCI_A3_A3(n_middle_A3, 2);
    for (auto &eV : SCI_A3_A3) {
      f_pair(list_middle_A3[eV[0]], list_middle_A3[eV[1]]);
    }
    // tilde{Dn} obtained from A3 + D(n-3)
    for (auto &v1 : list_middle_A3) {
      for (auto &v2 : list_expand_Dn)
        f_pair(v1, v2);
    }
    // tilde{Dn} obtained from Dk + Dl with k+l = n , k >= 4, l >= 4
    size_t n_expand_Dn = list_expand_Dn.size();
    SetCppIterator SCI_exp_Dk_Dl(n_expand_Dn, 2);
    for (auto &eV : SCI_exp_Dk_Dl) {
      size_t v1 = list_expand_Dn[eV[0]];
      size_t v2 = list_expand_Dn[eV[1]];
      if (VertToConn[v1] != VertToConn[v2])
        f_pair(v1, v2);
    }
    // tilde{E6} from A5+A1
    for (auto &v1 : list_extm2_AN) {
      if (VertToLocDim[v1] == 5) {
        for (auto &v2 : list_isolated)
          f_pair(v1, v2);
      }
    }
    // tilde{E7} from A5 + A2
    for (auto &v1 : list_extm1_AN) {
      if (VertToLocDim[v1] == 5) {
        for (auto &v2 : list_vert_A2)
          f_pair(v1, v2);
      }
    }
    // tilde{E7} from D6 + A1
    for (auto &v1 : list_non_expand_Dn) {
      if (VertToLocDim[v1] == 6) {
        for (auto &v2 : list_isolated)
          f_pair(v1, v2);
      }
    }
    // tilde{E8} from A7 + A1
    for (auto &v1 : list_extm1_AN) {
      if (VertToLocDim[v1] == 7) {
        for (auto &v2 : list_isolated)
          f_pair(v1, v2);
      }
    }
    // tilde{E8} from A4 + A4
    for (auto &v1 : list_extm1_AN) {
      if (VertToLocDim[v1] == 4) {
        for (auto &v2 : list_ends_A4) {
          if (VertToConn[v1] != VertToConn[v2])
            f_pair(v1, v2);
        }
      }
    }
    // tilde{E8} from D5 + A3
    for (auto &v1 : list_non_expand_Dn) {
      if (VertToLocDim[v1] == 5) {
        for (auto &v2 : list_ends_A3)
          f_pair(v1, v2);
      }
    }
    // tilde{E8} from E6 + A2
    for (auto &v1 : list_dist2_extrem_E6) {
      for (auto &v2 : list_vert_A2)
        f_pair(v1, v2);
    }
    // tilde{E8} from E6 + A2
    for (auto &v1 : list_dist3_extrem_E7) {
      for (auto &v2 : list_isolated)
        f_pair(v1, v2);
    }
    // tilde{F4} from A3 + A1
    for (auto &v1 : list_ends_A3) {
      for (auto &v2 : list_isolated) {
        MyVector<T> V = V_basic;
        V(v1) = val_four;
        V(v2) = val_single_edge;
        test_vector_and_insert(V);
      }
    }
    // tilde{F4} from A2 + A2
    for (auto &v1 : list_vert_A2) {
      for (auto &v2 : list_vert_A2) {
        if (VertToConn[v1] != VertToConn[v2]) {
          MyVector<T> V = V_basic;
          V(v1) = val_four;
          V(v2) = val_single_edge;
          test_vector_and_insert(V);
        }
      }
    }
    // tilde{F4} from tilde{B2} + A1
    for (auto &v1 : list_non_expand_Bn) {
      if (VertToLocDim[v1] == 3) {
        for (auto &v2 : list_isolated)
          f_pair(v1, v2);
      }
    }
  }
  //
  // Considering the case of 3 edges.
  //
  auto f_triple = [&](size_t w1, size_t w2, size_t w3) -> void {
#ifdef CHECK_EFFICIENT_ENUMERATION
    if (w1 == w2 || w1 == w3 || w2 == w3) {
      std::cerr << "We should have w1, w2, w3 all different\n";
      throw TerminalException{1};
    }
#endif
    MyVector<T> V = V_basic;
    V(w1) = val_single_edge;
    V(w2) = val_single_edge;
    V(w3) = val_single_edge;
    test_vector_and_insert(V);
  };
  // E6 formed from A2 + A2 + A1
  SetCppIterator SCI_A2_A2(n_A2, 2);
  for (auto &eV : SCI_A2_A2) {
    std::vector<size_t> const &vA2_1 = list_extremal_A2[eV[0]];
    std::vector<size_t> const &vA2_2 = list_extremal_A2[eV[1]];
    for (auto &v1 : vA2_1) {
      for (auto &v2 : vA2_2) {
        for (auto &v3 : list_isolated)
          f_triple(v1, v2, v3);
      }
    }
  }
  // E7 formed from A3 + A2 + A1
  for (auto &v1 : list_vert_A2) {
    for (auto &v2 : list_ends_A3) {
      for (auto &v3 : list_isolated)
        f_triple(v1, v2, v3);
    }
  }
  // E8 formed from A4 + A2 + A1
  for (auto &v1 : list_vert_A2) {
    for (auto &v2 : list_ends_A4) {
      for (auto &v3 : list_isolated)
        f_triple(v1, v2, v3);
    }
  }
  // D4 formed from A1 + A1 + A1
  SetCppIterator SCI_D4(n_isolated, 3);
  for (auto &eV : SCI_D4) {
    size_t v1 = list_isolated[eV[0]];
    size_t v2 = list_isolated[eV[1]];
    size_t v3 = list_isolated[eV[2]];
    f_triple(v1, v2, v3);
  }
  // Dn (for n > 4) formed from A(n-3) + A1 + A1
  SetCppIterator SCI_A(n_isolated, 2);
  for (auto &eV : SCI_A) {
    size_t v1 = list_isolated[eV[0]];
    size_t v2 = list_isolated[eV[1]];
    for (auto &LTerm : list_extremal_AN) {
      for (auto &v3 : LTerm)
        f_triple(v1, v2, v3);
    }
  }
  if (!DS.OnlySpherical) { // Doing the tilde{E6}; tilde{E7} and tilde{E8}
    // tilde{E6} from A2 + A2 + A2
    SetCppIterator SCI_A(n_A2, 3);
    for (auto &eV : SCI_A) {
      for (auto &v1 : list_extremal_A2[eV[0]]) {
        for (auto &v2 : list_extremal_A2[eV[1]]) {
          for (auto &v3 : list_extremal_A2[eV[2]])
            f_triple(v1, v2, v3);
        }
      }
    }
    // tilde{E7} from A3 + A3 + A1
    SetCppIterator SCI_tildeE7(n_A3, 2);
    for (auto &eV : SCI_tildeE7) {
      std::vector<size_t> const &vA3_1 = list_extremal_A3[eV[0]];
      std::vector<size_t> const &vA3_2 = list_extremal_A3[eV[1]];
      for (auto &v1 : vA3_1) {
        for (auto &v2 : vA3_2) {
          for (auto &v3 : list_isolated)
            f_triple(v1, v2, v3);
        }
      }
    }
    // tilde{E8} from A2 + A5 + A1
    for (auto &v1 : list_vert_A2) {
      for (auto &v2 : list_ends_A5) {
        for (auto &v3 : list_isolated)
          f_triple(v1, v2, v3);
      }
    }
    // tilde{B3} from A1 + A1 + A1
    SetCppIterator SCI_tildeB3(n_isolated, 3);
    T val;
    for (auto &eV : SCI_tildeB3) {
      MyVector<T> V = V_basic;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          val = val_single_edge;
          if (i == j)
            val = val_four;
          V(list_isolated[eV[j]]) = val;
        }
        test_vector_and_insert(V);
      }
    }
    // The tilde{Bn} for n >= 4 cases from B(n-2) + A1 + A1
    SetCppIterator SCI_tildeBn(n_isolated, 2);
    for (auto &eV : SCI_tildeBn) {
      size_t v1 = list_isolated[eV[0]];
      size_t v2 = list_isolated[eV[1]];
      for (auto &v3 : list_expand_Bn)
        f_triple(v1, v2, v3);
    }
    // tilde{Dn} formed from D(n-2) + A1 + A1
    SetCppIterator SCI_tildeDn(n_isolated, 2);
    for (auto &eV : SCI_tildeDn) {
      size_t v1 = list_isolated[eV[0]];
      size_t v2 = list_isolated[eV[1]];
      for (auto &v3 : list_expand_Dn)
        f_triple(v1, v2, v3);
    }
    // tilde{D5} formed from A3 + A1 + A1
    SetCppIterator SCI_tildeD5(n_isolated, 2);
    for (auto &eV : SCI_tildeD5) {
      size_t v1 = list_isolated[eV[0]];
      size_t v2 = list_isolated[eV[1]];
      for (auto &v3 : list_middle_A3)
        f_triple(v1, v2, v3);
    }
  }
  // Considering the case of 4 edges. Only tilde{D4} is possible
  if (!DS.OnlySpherical) { // Only tildeD4 is feasible, and it is not euclidean
    SetCppIterator SCI_B(n_isolated, 4);
    for (auto &eV : SCI_B) {
      MyVector<T> V = V_basic;
      for (auto &eVal : eV)
        V(list_isolated[eVal]) = val_single_edge;
      test_vector_and_insert(V);
    }
  }
  return ListExtensions;
}

template <typename T> MyMatrix<T> IrrCoxDyn_to_matrix(IrrCoxDyn<T> const &cd) {
  MyMatrix<T> M = Kernel_IrrCoxDyn_to_matrix<T>(cd);
  std::optional<IrrCoxDyn<T>> opt =
      RecognizeIrreducibleSphericalEuclideanDiagram(M);
  if (opt) {
    IrrCoxDyn<T> cd2 = *opt;
    if (cd.type != cd2.type || cd.dim != cd2.dim || cd.param != cd2.param) {
      std::cerr << "M=\n";
      WriteMatrix(std::cerr, M);
      std::cerr << "cd   type=" << cd.type << " dim" << cd.dim
                << " param=" << cd.param << "\n";
      std::cerr << "cd2  type=" << cd2.type << " dim" << cd2.dim
                << " param=" << cd2.param << "\n";
      std::cerr << "The recognition of the matrix did not yield the original "
                   "Coxeter-Dynkin diagram\n";
      throw TerminalException{1};
    }
    return M;
  }
  std::cerr << "The created matrix was not recognized. Some bug somewhere\n";
  throw TerminalException{1};
}

template <typename T>
std::optional<std::vector<IrrCoxDyn<T>>>
RecognizeSphericalEuclideanDiagram(const MyMatrix<T> &M) {
  std::vector<std::vector<size_t>> LConn = GetIrreducibleComponents(M);
  std::vector<IrrCoxDyn<T>> l_cd;
  for (auto &eConn : LConn) {
    size_t dim_res = eConn.size();
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
    std::cerr << "eConn =";
    for (auto &val : eConn)
      std::cerr << " " << val;
    std::cerr << "\n";
#endif
    MyMatrix<T> Mres(dim_res, dim_res);
    for (size_t i = 0; i < dim_res; i++)
      for (size_t j = 0; j < dim_res; j++)
        Mres(i, j) = M(eConn[i], eConn[j]);
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
    std::cerr << "Mres=\n";
    WriteMatrix(std::cerr, Mres);
#endif
    std::optional<IrrCoxDyn<T>> opt =
        RecognizeIrreducibleSphericalEuclideanDiagram(Mres);
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

template <typename T>
MyMatrix<T> string_to_coxdyn_matrix(std::string const &str) {
  std::vector<std::string> LStr = STRING_Split(str, "+");
  std::vector<MyMatrix<T>> LMat;
  size_t n_vert = 0;
  for (auto &s1 : LStr) {
    std::string s2 = s1;
    s2.erase(remove(s2.begin(), s2.end(), ' '), s2.end());
    IrrCoxDyn<T> cd = string_to_IrrCoxDyn<T>(s2);
    MyMatrix<T> M = IrrCoxDyn_to_matrix<T>(cd);
    LMat.push_back(M);
    n_vert += M.rows();
  }
  MyMatrix<T> Mret(n_vert, n_vert);
  for (size_t i = 0; i < n_vert; i++)
    for (size_t j = 0; j < n_vert; j++)
      Mret(i, j) = 2;
  size_t shift = 0;
  for (auto &M : LMat) {
    size_t len = M.rows();
    for (size_t i = 0; i < len; i++)
      for (size_t j = 0; j < len; j++)
        Mret(shift + i, shift + j) = M(i, j);
    shift += len;
  }
  return Mret;
}

template <typename T>
std::string coxdyn_matrix_to_string(MyMatrix<T> const &M) {
  std::optional<std::vector<IrrCoxDyn<T>>> opt =
      RecognizeSphericalEuclideanDiagram(M);
  if (opt) {
    const std::vector<IrrCoxDyn<T>> &l_irr = *opt;
    std::vector<std::string> l_str;
    for (auto &icd : l_irr)
      l_str.push_back(IrrCoxDyn_to_string(icd));
    std::sort(l_str.begin(), l_str.end());
    //
    std::string str = l_str[0];
    for (size_t i = 1; i < l_str.size(); i++)
      str += "+" + l_str[i];
    return str;
  }
  std::cerr << "M=\n";
  WriteMatrix(std::cerr, M);
  std::cerr
      << "The recognition failed so coxdyn_matrix_to_string cannot work\n";
  throw TerminalException{1};
}

template <typename T>
MyMatrix<T> ExtendMatrix(MyMatrix<T> const &M, MyVector<T> const &V) {
  size_t len = M.rows();
  MyMatrix<T> Mret(len + 1, len + 1);
  for (size_t i = 0; i < len; i++)
    for (size_t j = 0; j < len; j++)
      Mret(i, j) = M(i, j);
  for (size_t i = 0; i < len; i++) {
    Mret(len, i) = V(i);
    Mret(i, len) = V(i);
  }
  return Mret;
}

template <typename T>
MyMatrix<T> ExtendMatrixNorm(MyMatrix<T> const &M, MyVector<T> const &V,
                             T const &norm) {
  size_t len = M.rows();
  MyMatrix<T> Mret(len + 1, len + 1);
  for (size_t i = 0; i < len; i++)
    for (size_t j = 0; j < len; j++)
      Mret(i, j) = M(i, j);
  for (size_t i = 0; i < len; i++) {
    Mret(len, i) = V(i);
    Mret(i, len) = V(i);
  }
  Mret(len, len) = norm;
  return Mret;
}

template <typename T>
std::vector<MyVector<T>> FindDiagramExtensions(const MyMatrix<T> &M,
                                               const DiagramSelector &DS) {
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
    for (T val = val_single_edge; val < 128; val++)
      allowed_vals.push_back(val);
    if (!DS.OnlySpherical)
      allowed_vals.push_back(val_inf);
  }
  size_t dim = M.rows();
  std::vector<size_t> list_deg(dim);
  std::vector<size_t> list_n_higher(dim);
  std::vector<size_t> list_isolated;
  for (size_t i = 0; i < dim; i++) {
    size_t n_adj = 0;
    size_t n_higher = 0;
    for (size_t j = 0; j < dim; j++) {
      T val = M(i, j);
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
  std::vector<size_t>
      list_cand_for_triple_vertex; // Part of En, Dn, tilde or not
  for (auto &eConn : LConn) {
    size_t dim_res = eConn.size();
    MyMatrix<T> Mres(dim_res, dim_res);
    for (size_t i = 0; i < dim_res; i++)
      for (size_t j = 0; j < dim_res; j++)
        Mres(i, j) = M(eConn[i], eConn[j]);
    std::optional<IrrCoxDyn<T>> opt =
        RecognizeIrreducibleSphericalEuclideanDiagram(Mres);
    if (!opt) {
      std::cerr << "The diagram should have been recognized\n";
      throw TerminalException{1};
    }
    const IrrCoxDyn<T> &cd = *opt;
    // The diagrams that can be part of triple points are:
    //   An, Dn, Bn
    // The diagrams that cannot be part of triple points are:
    //   En, tildeEn, tildeDn, tildeBn, tildeCn, F4, tildeF4, G2, tildeG2
    if (cd.type == "A" || cd.type == "D" || cd.type == "B") {
      for (auto &eVert : eConn) {
        if (list_deg[eVert] <= 1) // Case 0 corresponds to A1
          list_cand_for_triple_vertex.push_back(eVert);
      }
    }
    if (cd.type == "A" && cd.dim == 3) {
      for (auto &eVert : eConn) {
        if (list_deg[eVert] ==
            2) // The middle vertex can match and get us tilde{D5}
          list_cand_for_triple_vertex.push_back(eVert);
      }
    }
  }
  // Consider the case of adding unconnected vector
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "FindDiagramExtensions, step 2\n";
#endif
  MyVector<T> V_basic(dim);
  for (size_t i = 0; i < dim; i++)
    V_basic(i) = val_comm;
#ifdef DEBUG_COXETER_DYNKIN_COMBINATORICS
  std::cerr << "FindDiagramExtensions, step 3\n";
#endif
  size_t n_diagram_considered = 0;
  size_t n_diagram_match = 0;
  auto test_vector_and_insert = [&](const MyVector<T> &V) -> void {
    MyMatrix<T> Mtest(dim + 1, dim + 1);
    for (size_t i = 0; i < dim; i++)
      for (size_t j = 0; j < dim; j++)
        Mtest(i, j) = M(i, j);
    for (size_t i = 0; i < dim; i++) {
      Mtest(i, dim) = V(i);
      Mtest(dim, i) = V(i);
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
    for (size_t i = 0; i < dim; i++) {
      MyVector<T> V = V_basic;
      V(i) = val_single_edge;
      test_vector_and_insert(V);
      V(i) = val_four;
      test_vector_and_insert(V);
      V(i) = val_six; // This covers G2 and tilde{G2}
      test_vector_and_insert(V);
    }
    for (auto &eIsol : list_isolated) {
      MyVector<T> V = V_basic;
      if (!DS.OnlySpherical) {
        V(eIsol) = val_inf;
        test_vector_and_insert(V); // I1(infinity), always works.
      }
    }
  } else {
    for (size_t i = 0; i < dim; i++) {
      // Here we have an arbitrary value
      for (auto &val : allowed_vals) {
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
    for (size_t i = 0; i < dim; i++) {
      for (size_t j = i + 1; j < dim; j++) {
        MyVector<T> V = V_basic;
        V(i) = val_single_edge;
        V(j) = val_single_edge;
        test_vector_and_insert(V);
      }
    }
    for (size_t i = 0; i < dim; i++) {
      for (size_t j = 0; j < dim; j++) {
        if (i != j) {
          MyVector<T> V = V_basic;
          V(i) = val_single_edge;
          V(j) = val_four;
          test_vector_and_insert(V);
        }
      }
    }
    if (!DS.OnlySpherical) { // Only tildeG2 and tilde{C2} are possible here
      for (size_t i = 0; i < n_isolated; i++) {
        for (size_t j = 0; j < n_isolated; j++) {
          if (i != j) {
            MyVector<T> V = V_basic;
            V(list_isolated[i]) = val_single_edge;
            V(list_isolated[j]) = val_six;
            test_vector_and_insert(V); // For tilde{G2}
          }
        }
      }
      for (size_t i = 0; i < n_isolated; i++) {
        for (size_t j = i + 1; j < n_isolated; j++) {
          MyVector<T> V = V_basic;
          V(list_isolated[i]) = val_four;
          V(list_isolated[j]) = val_four;
          test_vector_and_insert(V); // For tilde{C2}
        }
      }
    }
  } else {
    for (size_t i = 0; i < dim; i++) {
      for (size_t j = 0; j < dim; j++) {
        if (i != j) {
          for (T val = val_single_edge; val <= val_six; val++) {
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
  size_t n_cand_triple = list_cand_for_triple_vertex.size();
  SetCppIterator SCI_A(n_cand_triple, 3);
  for (auto &eV : SCI_A) {
    MyVector<T> V = V_basic;
    for (auto &eVal : eV)
      V(list_cand_for_triple_vertex[eVal]) = val_single_edge;
    test_vector_and_insert(V);
  }
  if (!DS.OnlySpherical) { // Considering now the tilde{B3} cases
    SetCppIterator SCI_B(n_isolated, 3);
    T val;
    for (auto &eV : SCI_B) {
      MyVector<T> V = V_basic;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
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
    SetCppIterator SCI_B(n_isolated, 4);
    for (auto &eV : SCI_B) {
      MyVector<T> V = V_basic;
      for (auto &eVal : eV)
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
  std::cerr << "Stats : |ListExtensions|=" << ListExtensions.size()
            << " n_diagram_considered=" << n_diagram_considered
            << " n_diagram_match=" << n_diagram_match << "\n";
  return ListExtensions;
}

template <typename T, typename Tint>
std::pair<MyMatrix<T>, MyMatrix<T>>
ComputeCoxeterMatrix(MyMatrix<T> const &G,
                     std::vector<MyVector<Tint>> const &l_root) {
  int n = G.rows();
  auto get_scal = [&](int i, int j) -> T {
    T sum = 0;
    for (int u = 0; u < n; u++)
      for (int v = 0; v < n; v++)
        sum += G(u, v) * l_root[i](u) * l_root[j](v);
    return sum;
  };
  T cossqr_val3 = T(1) / T(4);
  T cossqr_val4 = T(1) / T(2);
  T cossqr_val6 = T(3) / T(4);
  T cossqr_valInf = 1;
  T pr_inf = practical_infinity<T>();
  auto get_cossqr_scal = [&](int i, int j) -> std::pair<T, T> {
    T scal12 = get_scal(i, j);
    if (i == j) {
      return {scal12, scal12};
    } else {
      if (scal12 == 0)
        return {2, scal12};
      T scal11 = get_scal(i, i);
      T scal22 = get_scal(j, j);
      T quot = (scal12 * scal12) / (scal11 * scal22);
      if (quot == cossqr_val3)
        return {3, scal12};
      if (quot == cossqr_val4)
        return {4, scal12};
      if (quot == cossqr_val6)
        return {6, scal12};
      if (quot == cossqr_valInf)
        return {pr_inf, scal12};
      std::cerr << "i=" << i << " j=" << j << "\n";
      std::cerr << "scal12=" << scal12 << " scal11=" << scal11
                << " scal22=" << scal22 << "\n";
      std::cerr << "l_root=\n";
      WriteMatrix(std::cerr, MatrixFromVectorFamily(l_root));
      std::cerr << "Failed to find matching entry quot=" << quot << "\n";
      throw TerminalException{1};
    }
  };
  size_t n_root = l_root.size();
  MyMatrix<T> CoxMat(n_root, n_root);
  MyMatrix<T> ScalMat(n_root, n_root);
  for (size_t i = 0; i < n_root; i++)
    for (size_t j = 0; j < n_root; j++) {
      std::pair<T, T> ep = get_cossqr_scal(i, j);
      CoxMat(i, j) = ep.first;
      ScalMat(i, j) = ep.second;
    }
  return {std::move(CoxMat), std::move(ScalMat)};
}

template <typename T> struct Possible_Extension {
  MyVector<T> u_component;
  T res_norm;
  T e_norm;
};

template <typename T, typename Tint>
std::vector<Possible_Extension<T>>
ComputePossibleExtensions(MyMatrix<T> const &G,
                          std::vector<MyVector<Tint>> const &l_root,
                          std::vector<T> const &l_norm, bool only_spherical) {
#ifdef DEBUG_COMPUTE_POSSIBLE_EXTENSIONS
  std::cerr << "------------------------------------ ComputePossibleExtension "
               "---------------------------------\n";
  std::cerr << "only_spherical=" << only_spherical << "\n";
#endif
  DiagramSelector DS;
  DS.OnlySimplyLaced = false;
  DS.OnlyLorentzianAdmissible = true;
  DS.OnlySpherical = only_spherical;
  //
  //  std::cerr << "G=\n";
  //  WriteMatrixGAP(std::cerr, G);
  //  std::cerr << "\n";
  //  std::cerr << "l_root=\n";
  //  for (auto & e_root : l_root)
  //    std::cerr << "e_root=" << StringVectorGAP(e_root) << "\n";

  std::pair<MyMatrix<T>, MyMatrix<T>> ep = ComputeCoxeterMatrix(G, l_root);
  const MyMatrix<T> &CoxMat = ep.first;
  const MyMatrix<T> &ScalMat = ep.second;
  MyMatrix<T> ScalMatInv = Inverse(ScalMat);
#ifdef DEBUG_COMPUTE_POSSIBLE_EXTENSIONS
  std::cerr << "ScalMat=\n";
  WriteMatrix(std::cerr, ScalMat);
  std::cerr << "CoxMat=\n";
  WriteMatrix(std::cerr, CoxMat);
  std::cerr << "Symbol of M=" << coxdyn_matrix_to_string(CoxMat) << "\n";
#endif
  int dim = G.rows();
  int n_root = l_root.size();
  std::vector<MyVector<T>> l_vect = FindDiagramExtensions_Efficient(CoxMat, DS);
#ifdef CHECK_EFFICIENT_ENUMERATION
  std::vector<MyVector<T>> l_vect_B = FindDiagramExtensions(CoxMat, DS);
  if (l_vect.size() != l_vect_B.size()) {
    std::cerr << "The two enumeration codes return different results\n";
    std::set<MyVector<T>> s_vect;
    for (auto &eV : l_vect)
      s_vect.insert(eV);
    std::set<MyVector<T>> s_vect_B;
    for (auto &eV : l_vect_B)
      s_vect_B.insert(eV);
    std::cerr << "In s_vect but not in s_vect_B=\n";
    auto prt = [&](MyVector<T> const &eV) -> void {
      MyMatrix<T> CoxMatExp(n_root + 1, n_root + 1);
      for (int i = 0; i < n_root; i++)
        for (int j = 0; j < n_root; j++)
          CoxMatExp(i, j) = CoxMat(i, j);
      for (int i = 0; i < n_root; i++) {
        CoxMatExp(i, n_root) = eV(i);
        CoxMatExp(n_root, i) = eV(i);
      }
      CoxMatExp(n_root, n_root) = 2;
      std::string symb = coxdyn_matrix_to_string(CoxMatExp);
      std::cerr << "V=" << StringVectorGAP(eV) << " symb=" << symb << "\n";
    };
    for (auto &eV : s_vect)
      if (s_vect_B.count(eV) == 0)
        prt(eV);
    std::cerr << "In s_vect_B but not in s_vect=\n";
    for (auto &eV : s_vect_B)
      if (s_vect.count(eV) == 0)
        prt(eV);
    std::cerr << "only_spherical = " << only_spherical << "\n";
    std::cerr << "|l_vect|=" << l_vect.size()
              << " |l_vect_B|=" << l_vect_B.size() << "\n";
    std::cerr << "|s_vect|=" << s_vect.size()
              << " |s_vect_B|=" << s_vect_B.size() << "\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_COMPUTE_POSSIBLE_EXTENSIONS
  std::cerr << "|l_vect|=" << l_vect.size() << "\n";
#endif
  T val2 = 0;
  T val3 = T(1) / T(4);
  T val4 = T(1) / T(2);
  T val6 = T(3) / T(4);
  T valInfinity = 1;
  auto get_cos_square = [&](T val) -> T {
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
  auto get_entry = [&](MyVector<T> const &e_vect, T const &e_norm) -> void {
    //    std::cerr << "---------------- e_norm=" << e_norm << " e_vect=" <<
    //    StringVectorGAP( e_vect) << "\n";
    MyVector<T> l_scal(n_root);
    for (int i = 0; i < n_root; i++) {
      T val = e_vect(i);
      T cos_square = get_cos_square(val);
      T scal_square = cos_square * CoxMat(i, i) * e_norm;
      std::optional<T> opt = UniversalSquareRoot(scal_square);
      //      std::cerr << "i=" << i << " cos_square=" << cos_square << "
      //      CoxMat(i,i)=" << CoxMat(i,i) << " e_norm=" << e_norm << "
      //      scal_square=" << scal_square << "\n"; std::cerr << "i=" << i << "
      //      scal_square=" << scal_square << "\n";
      if (!opt) {
        //        std::cerr << "   Failed to match\n";
        return;
      }
      T scal = -*opt;
      //      std::cerr << "     scal=" << scal << "\n";
      l_scal(i) = scal;
    }
#ifdef DEBUG_COMPUTE_POSSIBLE_EXTENSIONS
    std::cerr << "---------------- e_norm=" << e_norm
              << " e_vect=" << StringVectorGAP(e_vect) << "\n";
#endif
    //    std::cerr << "Scalar products found : l_scal = " <<
    //    StringVectorGAP(l_scal) << "\n";
    /* So, we have computed l_scal(i) = alpha.dot.ui = u.dot.ui
       If u = sum wi u_i then w = G^{-1} l_scal
       eNorm = w.dot.w  is the Euclidean norm of u.
     */
    MyVector<T> w = ScalMatInv * l_scal;
    //    std::cerr << "w = " << StringVectorGAP(w) << "\n";
    T eNorm = l_scal.dot(w);
    T res_norm = e_norm - eNorm;
#ifdef DEBUG_COMPUTE_POSSIBLE_EXTENSIONS
    std::cerr << "  eNorm=" << eNorm << " res_norm=" << res_norm << "\n";
#endif
    MyVector<T> u_component = ZeroVector<T>(dim);
    for (int i = 0; i < n_root; i++)
      u_component += w(i) * UniversalVectorConversion<T, Tint>(l_root[i]);
#ifdef SANITY_CHECK
    if (res_norm == 0) {
      std::vector<MyVector<T>> l_root_tot;
      for (auto &ev_tint : l_root)
        l_root_tot.push_back(UniversalVectorConversion<T, Tint>(ev_tint));
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
      //      std::cerr << "Symbol of CoxMatExt=" <<
      //      coxdyn_matrix_to_string(CoxMatExt) << "\n";
    }
#endif
#ifdef DEBUG_COMPUTE_POSSIBLE_EXTENSIONS
    std::cerr << "  u_component=";
    WriteVectorGAP(std::cerr, u_component);
    std::cerr << "\n";
#endif
    l_extensions.push_back({u_component, res_norm, e_norm});
  };
#ifdef DEBUG_COMPUTE_POSSIBLE_EXTENSIONS
  std::cerr << "ComputePossibleExtensions, step 5\n";
#endif
  size_t len = l_norm.size();
  Face status_norm(len);
  for (auto &e_vect : l_vect) {
    for (size_t i_norm = 0; i_norm < len; i_norm++)
      status_norm[i_norm] = 1;
    for (int i = 0; i < n_root; i++) {
      T norm = CoxMat(i, i);
      for (size_t i_norm = 0; i_norm < len; i_norm++) {
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
    for (size_t i_norm = 0; i_norm < len; i_norm++)
      if (status_norm[i_norm] == 1)
        get_entry(e_vect, l_norm[i_norm]);
  }
#ifdef DEBUG_COMPUTE_POSSIBLE_EXTENSIONS
  std::cerr << "ComputePossibleExtensions, step 6 |l_extensions|="
            << l_extensions.size() << "\n";
#endif
  return l_extensions;
}

#endif  // SRC_LORENTZIAN_COXETER_DYNKIN_H_

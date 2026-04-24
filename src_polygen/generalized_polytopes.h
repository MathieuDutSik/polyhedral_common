// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_GENERALIZED_POLYTOPES_H_
#define SRC_ROBUST_COVERING_GENERALIZED_POLYTOPES_H_

// clang-format off
#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalBasic.h"
#include "POLY_LinearProgramming.h"
#include "POLY_DirectDualDesc.h"
#include "POLY_Fundamental.h"
#include "POLY_RedundancyElimination.h"
#include "boost_serialization.h"
#include <boost/dynamic_bitset/serialization.hpp>
// clang-format on

/*
  This is a framework for dealing with objects defined as
  unions of full dimensional polytopes. We call those
  generalized polytopes.
  ---
  For classical polytopes, the description by vertices
  and facets is both unique and useful. But that is not
  going to be the case for geometric union.
  This indicates two kinds of description:
  * A description as a finite union
  * A description as unique set of faces.
  ---
  The description as a finite union is going to be
  fundamental for our work.
  We store the facets and vertices of each of the
  components. The umion should be disjoint.
  What can we compute as a description as finite union
  * Do checks:
    + Check pairwise intersections.
    + Check for possible unions. Just heuristics. We
      cannot have a canonical schema.
  * Intersection [not useful right now, but could be]
    + Compute the pairwise intersections.
    + Keep the non-empty ones.
  * Difference A - B [useful for computing]
    + The difference can be done cell by cell.
    + So, we need to compute P - Q for P and Q polytopes.
    + Select a violated facet (f(x) >= 0) of Q, write
      P = P1 union P2 with f(x) doing the splitting.
      P1 with f(x) <= 0 is one component.
      Then iterate with P2 - Q until finished.
  * A decomposition in connected components.
    + Compute the facets of each component.
    + We can look for the pairwise intersections are
      correct.
    + Adjacency is variable. It can be sharing a facet,
      or just a vertex.
  * Get the facet as a polytope union. Very useful.
    Relatively easy to compute.
  * Get an interior point to a polytope union.
    It just suffices to have one for one component.
    But maybe we could select a better one.
  ---
  The functionalities that we want from the unique
  description:
  * Compute the vertices.
    x The computation of each components gets us
    an initial list of candidate vertices.
    x So, we have to test them. Which means that we have
    to test whether the corresponding facets are indeed
    real facets.
  * Compute the face lattice. But each face can be a
    generalized polytope.
  * Be able to test for equality of generalized polytopes.
  * Be able to compute the stabilizer of a generalized
    polytope and compute equivalences.
  ---
 */

#ifdef DEBUG
#define DEBUG_GENERALIZED_POLYTOPE
#endif


template <typename T> struct SinglePolytope {
  MyMatrix<T> EXT;
  MyMatrix<T> FAC;
  vectface facets;
  SinglePolytope(SinglePolytope<T> const& x) : EXT(x.EXT), FAC(x.FAC), facets(x.facets.clone()) {
  }
  SinglePolytope(MyMatrix<T> const& _EXT, MyMatrix<T> const& _FAC, vectface const& _facets) : EXT(_EXT), FAC(_FAC), facets(_facets.clone()) {
  }
  SinglePolytope() : facets() {
  }
  SinglePolytope<T>& operator=(SinglePolytope<T> const& x) {
    EXT = x.EXT;
    FAC = x.FAC;
    facets = x.facets.clone();
    return *this;
  }
  std::string ext_string_gap() const {
    return StringMatrixGAP(EXT);
  }
};

template <typename T>
void WriteEntryGAP(std::ostream &os_out, SinglePolytope<T> const &sp) {
  os_out << "rec(EXT:=";
  WriteMatrixGAP(os_out, sp.EXT);
  os_out << ", FAC:=";
  WriteMatrixGAP(os_out, sp.FAC);
  os_out << ", facets:=" << StringVectfaceGAP(sp.facets) << ")";
}

template <typename T>
SinglePolytope<T> get_single_polytope(MyMatrix<T> const &FAC, MyMatrix<T> const &EXT) {
  int n_fac = FAC.rows();
  int n_ext = EXT.rows();
  int dim = EXT.cols();
  vectface facets(n_ext);
  for (int i_fac = 0; i_fac < n_fac; i_fac++) {
    Face f(n_ext);
    for (int i_ext = 0; i_ext < n_ext; i_ext++) {
      T scal(0);
      for (int i = 0; i < dim; i++) {
        scal += FAC(i_fac, i) * EXT(i_ext, i);
      }
      if (scal == 0) {
        f[i_ext] = 1;
      }
    }
    facets.push_back(f);
  }
  MyMatrix<T> EXTret(n_ext, dim);
  for (int i_ext=0; i_ext<n_ext; i_ext++) {
    T val0(1);
    if (EXT(i_ext, 0) > 0) {
      val0 = EXT(i_ext, 0);
    }
    for (int i=0; i<dim; i++) {
      EXTret(i_ext, i) = EXT(i_ext, i) / val0;
    }
  }
  MyMatrix<T> FACret(n_fac, dim);
  for (int i_fac=0; i_fac<n_fac; i_fac++) {
    MyVector<T> eFAC = GetMatrixRow(FAC, i_fac);
    MyVector<T> fFAC = ScalarCanonicalizationVector(eFAC);
    AssignMatrixRow(FACret, i_fac, fFAC);
  }
  return SinglePolytope<T>(EXTret, FACret, std::move(facets));
}

template <typename T>
SinglePolytope<T> mat_product(SinglePolytope<T> const& sp, MyMatrix<T> const& P) {
  MyMatrix<T> EXTnew = sp.EXT * P;
  MyMatrix<T> Pcgr = CongrMap(P);
  MyMatrix<T> FACnew = sp.FAC * Pcgr;
  vectface facets = sp.facets.clone();
  return SinglePolytope<T>(EXTnew, FACnew, std::move(facets));
}


template<typename T>
std::optional<SinglePolytope<T>> singlepolytope_halfspace_int(SinglePolytope<T> const& sp, MyVector<T> const& eFAC, std::ostream &os) {
  MyMatrix<T> FACtot = ConcatenateMatVec(sp.FAC, eFAC);
  if (!IsFullDimensional(FACtot, os)) {
    return {};
  }
  std::vector<int> ListIrred = get_non_redundant_indices(FACtot, os);
  MyMatrix<T> FAC = SelectRow(FACtot, ListIrred);
  MyMatrix<T> EXT = DirectDualDescription_mat(FAC, os);
  return get_single_polytope(FAC, EXT);
}

template<typename T>
std::optional<SinglePolytope<T>> singlepolytope_halfspaces_int(SinglePolytope<T> const& sp, MyMatrix<T> const& FAC_inp, std::ostream &os) {
  MyMatrix<T> FACtot = Concatenate(sp.FAC, FAC_inp);
  if (!IsFullDimensional(FACtot, os)) {
    return {};
  }
  std::vector<int> ListIrred = get_non_redundant_indices(FACtot, os);
  MyMatrix<T> FAC = SelectRow(FACtot, ListIrred);
  MyMatrix<T> EXT = DirectDualDescription_mat(FAC, os);
  return get_single_polytope(FAC, EXT);
}






template<typename T>
struct ConvexBoundary {
  MyVector<T> V; // The orthogonal
  MyMatrix<T> NSP;
  SinglePolytope<T> sp; // The polytope in the face.
  MyVector<T> get_isobarycenter([[maybe_unused]] std::ostream &os) const {
#ifdef DEBUG_GENERALIZED_POLYTOPE
    os << "GP: get_isobarycenter, Before Isobarycenter\n";
#endif
    MyVector<T> eIso1 = Isobarycenter(sp.EXT);
#ifdef DEBUG_GENERALIZED_POLYTOPE
    os << "GP: get_isobarycenter, eIso1=" << StringVectorGAP(eIso1) << "\n";
#endif
    MyVector<T> eIso2 = NSP.transpose() * eIso1;
#ifdef DEBUG_GENERALIZED_POLYTOPE
    os << "GP: get_isobarycenter, eIso2=" << StringVectorGAP(eIso2) << "\n";
#endif
    T val = eIso2(0);
    for (int u=0; u<eIso2.size(); u++) {
      eIso2(u) = eIso2(u) / val;
    }
#ifdef DEBUG_GENERALIZED_POLYTOPE
    os << "GP: get_isobarycenter, eIso2_B=" << StringVectorGAP(eIso2) << "\n";
#endif
    return eIso2;
  }
  MyVector<T> random_interior_point(int const&N, std::ostream& os) const {
#ifdef DEBUG_GENERALIZED_POLYTOPE
    os << "GP: random_interior_point, Before Isobarycenter\n";
#endif
    MyVector<T> eIso1 = random_interior_pt(sp.EXT, N, os);
#ifdef DEBUG_GENERALIZED_POLYTOPE
    os << "GP: random_interior_point, eIso1=" << StringVectorGAP(eIso1) << "\n";
#endif
    MyVector<T> eIso2 = NSP.transpose() * eIso1;
#ifdef DEBUG_GENERALIZED_POLYTOPE
    os << "GP: random_interior_point, eIso2_A=" << StringVectorGAP(eIso2) << "\n";
#endif
    T val = eIso2(0);
    for (int u=0; u<eIso2.size(); u++) {
      eIso2(u) = eIso2(u) / val;
    }
#ifdef DEBUG_GENERALIZED_POLYTOPE
    os << "GP: random_interior_point, eIso2_B=" << StringVectorGAP(eIso2) << "\n";
#endif
    return eIso2;
  }
  std::vector<MyVector<T>> get_list_vertices() const {
    std::vector<MyVector<T>> l_vertices;
    MyMatrix<T> EXT2 = sp.EXT * NSP;
    int dim = EXT2.cols();
    for (int i_row=0; i_row<EXT2.rows(); i_row++) {
      MyVector<T> V = GetMatrixRow(EXT2, i_row);
      T val = V(0);
      for (int u=0; u<dim; u++) {
        V(u) = V(u) / val;
      }
      l_vertices.emplace_back(std::move(V));
    }
    return l_vertices;
  }
};

template <typename T>
void WriteEntryGAP(std::ostream &os_out, ConvexBoundary<T> const &cb) {
  os_out << "rec(V:=" << StringVectorGAP(cb.V);
  os_out << ", NSP:=";
  WriteMatrixGAP(os_out, cb.NSP);
  os_out << ", sp:=";
  WriteEntryGAP(os_out, cb.sp);
  os_out << ")";
}

template<typename T>
std::vector<int> get_adjacent_facet_indices(SinglePolytope<T> const& sp, int const& i_fac) {
  int dim = sp.FAC.cols();
  int n_fac = sp.FAC.rows();
  int n_ext = sp.EXT.rows();
  std::vector<int> l_idx_facet;
  Face f1 = sp.facets[i_fac];
  for (int j_fac=0; j_fac<n_fac; j_fac++) {
    if (i_fac != j_fac) {
      Face f2 = sp.facets[j_fac];
      int n_incd = 0;
      for (int i_ext=0; i_ext<n_ext; i_ext++) {
        if (f1[i_ext] == 1 && f2[i_ext] == 1) {
          n_incd += 1;
        }
      }
      bool is_facet = false;
      if (n_incd >= dim - 2) {
        MyMatrix<T> EXTincd(n_incd, dim);
        int i_incd = 0;
        for (int i_ext=0; i_ext<n_ext; i_ext++) {
          if (f1[i_ext] == 1 && f2[i_ext] == 1) {
            for (int i=0; i<dim; i++) {
              EXTincd(i_incd, i) = sp.EXT(i_ext, i);
            }
            i_incd += 1;
          }
        }
        int rnk = RankMat(EXTincd);
        if (rnk == dim - 2) {
          is_facet = true;
        }
      }
      if (is_facet) {
        l_idx_facet.push_back(j_fac);
      }
    }
  }
  return l_idx_facet;
}



template<typename T>
ConvexBoundary<T> get_convex_boundary(SinglePolytope<T> const& sp, int const& i_fac, [[maybe_unused]] std::ostream &os) {
  int dim = sp.FAC.cols();
  MyVector<T> V = GetMatrixRow(sp.FAC, i_fac);
  MyMatrix<T> NSP = NullspaceMatSingleVectorExt(V);
#ifdef DEBUG_GENERALIZED_POLYTOPE_DISABLE
  os << "GP: get_convex_boundary, NSP=\n";
  WriteMatrix(os, NSP);
#endif
  int n_ext = sp.EXT.rows();
  std::vector<int> l_idx_facet = get_adjacent_facet_indices(sp, i_fac);
  Face f1 = sp.facets[i_fac];
  std::vector<int> f1_map(n_ext, -1);
  int pos = 0;
  for (int i_ext=0; i_ext<n_ext; i_ext++) {
    if (f1[i_ext] == 1) {
      f1_map[i_ext] = pos;
      pos += 1;
    }
  }
  MyMatrix<T> EXT2 = SelectRow(sp.EXT, f1);
  SolutionMatRepetitive<T> smr(NSP);
  int n_ext_ret = EXT2.rows();
  MyMatrix<T> EXT_ret(n_ext_ret, dim-1);
  for (int i_ext=0; i_ext<n_ext_ret; i_ext++) {
    MyVector<T> v1 = GetMatrixRow(EXT2, i_ext);
    std::optional<MyVector<T>> opt = smr.GetSolution(v1);
    if (!opt) {
      std::cerr << "Failed to find a solution to the system\n";
      throw TerminalException{1};
    }
    AssignMatrixRow(EXT_ret, i_ext, *opt);
  }
  int n_fac_ret = l_idx_facet.size();
  MyMatrix<T> FAC_ret(n_fac_ret, dim-1);
  for (int j_fac_ret=0; j_fac_ret<n_fac_ret; j_fac_ret++) {
    int j_fac = l_idx_facet[j_fac_ret];
    Face f2 = sp.facets[j_fac];
    Face f(n_ext_ret);
    for (int i_ext=0; i_ext<n_ext; i_ext++) {
      if (f1[i_ext] == 1 && f2[i_ext] == 1) {
        int pos = f1_map[i_ext];
        f[pos] = 1;
      }
    }
    MyVector<T> eFAC = FindFacetInequality(EXT_ret, f);
    AssignMatrixRow(FAC_ret, j_fac_ret, eFAC);
  }
  SinglePolytope<T> sp_ret = get_single_polytope(FAC_ret, EXT_ret);
  return {V, NSP, std::move(sp_ret)};
}

template<typename T>
std::optional<ConvexBoundary<T>> convexboundary_halfspace_int(ConvexBoundary<T> const& cb, MyVector<T> const& eFAC, std::ostream &os) {
  MyVector<T> eFAC_call = cb.NSP * eFAC;
  std::optional<SinglePolytope<T>> opt = singlepolytope_halfspace_int(cb.sp, eFAC_call, os);
  if (!opt) {
    return {};
  }
  SinglePolytope<T> sp = *opt;
  ConvexBoundary<T> cb_ret{cb.V, cb.NSP, sp};
  return cb_ret;
}

template<typename T>
std::optional<ConvexBoundary<T>> convexboundary_halfspaces_int(ConvexBoundary<T> const& cb, MyMatrix<T> const& FAC, std::ostream &os) {
  MyMatrix<T> FAC_call = FAC * cb.NSP.transpose();
  std::optional<SinglePolytope<T>> opt = singlepolytope_halfspaces_int(cb.sp, FAC_call, os);
  if (!opt) {
    return {};
  }
  SinglePolytope<T> sp = *opt;
  ConvexBoundary<T> cb_ret{cb.V, cb.NSP, sp};
  return cb_ret;
}

template<typename T>
int get_matching_face_position(MyMatrix<T> const& FAC, MyVector<T> const& eFAC) {
  int n_fac = FAC.rows();
  int dim = FAC.cols();
  auto f_is_match=[&](int const& j_fac, int sign) -> bool {
    for (int i=0; i<dim; i++) {
      T val = sign * eFAC(i);
      if (FAC(j_fac, i) != val) {
        return false;
      }
    }
    return true;
  };
  for (int i_fac=0; i_fac<n_fac; i_fac++) {
    if (f_is_match(i_fac, 1)) {
      return i_fac;
    }
    if (f_is_match(i_fac, -1)) {
      return i_fac;
    }
  }
  return -1;
}

template <typename T>
MyVector<T> get_interior_facet_pt(SinglePolytope<T> const& sp, int i_facet) {
  int dim = sp.EXT.cols();
  int n_ext = sp.EXT.rows();
  MyVector<T> V = ZeroVector<T>(dim);
  Face f = sp.facets[i_facet];
  int incd = 0;
  for (int i_ext=0; i_ext<n_ext; i_ext++) {
    if (f[i_ext] == 1) {
      for (int i=0; i<dim; i++) {
        V(i) += sp.EXT(i_ext, i);
      }
      incd += 1;
    }
  }
  for (int i=0; i<dim; i++) {
    V(i) = V(i) / incd;
  }
  return V;
}



template <typename T>
SinglePolytope<T> generate_single_polytope(MyMatrix<T> const &FACinput,
                                           std::ostream &os) {
#ifdef DEBUG_GENERALIZED_POLYTOPE_DISABLE
  os << "GP: generate_single_polytope FACinput=\n";
  WriteMatrix(os, FACinput);
#endif
  std::vector<int> ListIrred = get_non_redundant_indices(FACinput, os);
#ifdef DEBUG_GENERALIZED_POLYTOPE_DISABLE
  os << "GP: generate_single_polytope |ListIrred|=" << ListIrred.size() << "\n";
#endif
  MyMatrix<T> FAC = SelectRow(FACinput, ListIrred);
#ifdef DEBUG_GENERALIZED_POLYTOPE_DISABLE
  os << "GP: generate_single_polytope FAC=\n";
  WriteMatrix(os, FAC);
#endif
  MyMatrix<T> EXT = DirectDualDescription_mat(FAC, os);
  return get_single_polytope(FAC, EXT);
}

template <typename T> struct GeneralizedPolytope {
  int dim;
  std::vector<SinglePolytope<T>> polytopes;
  size_t size() const {
    return polytopes.size();
  }
  bool empty() const {
    return polytopes.empty();
  }
  std::string ext_string_gap() const {
    std::ostringstream os;
    bool is_first=true;
    os << "[";
    for (auto & sp: polytopes) {
      if (!is_first) {
        os << ",";
      }
      is_first=false;
      os << sp.ext_string_gap();
    }
    os << "]";
    return os.str();
  }
};

template <typename T>
void WriteEntryGAP(std::ostream &os_out, GeneralizedPolytope<T> const &gp) {
  os_out << "rec(dim:=" << gp.dim << ", polytopes:=[";
  bool IsFirst = true;
  for (auto &sp : gp.polytopes) {
    if (!IsFirst) {
      os_out << ",";
    }
    IsFirst = false;
    WriteEntryGAP(os_out, sp);
  }
  os_out << "])";
}

template <typename T>
GeneralizedPolytope<T> list_ext_to_generalizedpolytope(int const& dim, std::vector<MyMatrix<T>> const& list_ext) {
  std::vector<SinglePolytope<T>> polytopes;
  for (size_t i=0; i<list_ext.size(); i++) {
    MyMatrix<T> const& EXT = list_ext[i];
    MyMatrix<T> FAC = DirectDualDescription_mat(EXT, std::cerr);
    SinglePolytope<T> sp = get_single_polytope(FAC, EXT);
    polytopes.push_back(sp);
  }
  GeneralizedPolytope<T> gp{dim, polytopes};
  return gp;
}

template <typename T>
GeneralizedPolytope<T> mat_product(GeneralizedPolytope<T> const& gp, MyMatrix<T> const& P) {
  std::vector<SinglePolytope<T>> polytopes;
  for (auto &sp : gp.polytopes) {
    SinglePolytope<T> sp_new = mat_product(sp, P);
    polytopes.push_back(sp_new);
  }
  int dim = gp.dim;
  GeneralizedPolytope<T> gp_new{dim, polytopes};
  return gp_new;
}


template <typename T>
bool check_pairwise_intersection(GeneralizedPolytope<T> const &gp,
                                 std::ostream &os) {
  size_t n_comp = gp.size();
  for (size_t i = 0; i < n_comp; i++) {
    for (size_t j = i + 1; j < n_comp; j++) {
      MyMatrix<T> FAC = Concatenate(gp.polytopes[i].FAC, gp.polytopes[j].FAC);
      bool test = IsFullDimensional(FAC, os);
      if (test) {
        return false;
      }
    }
  }
  return true;
}

template <typename T>
std::unordered_map<MyVector<T>, size_t>
get_map_total_vertices(GeneralizedPolytope<T> const &gp,
                       [[maybe_unused]] std::ostream &os) {
  std::unordered_map<MyVector<T>, size_t> map_total_vertices;
  size_t pos = 1;
  for (auto &polytope : gp.polytopes) {
    int n_ext = polytope.EXT.rows();
    for (int i_ext = 0; i_ext < n_ext; i_ext++) {
      MyVector<T> eEXT = GetMatrixRow(polytope.EXT, i_ext);
      size_t &find_pos = map_total_vertices[eEXT];
      if (find_pos == 0) {
        find_pos = pos;
        pos += 1;
      }
    }
  }
  for (auto &[vertex, pos] : map_total_vertices) {
    pos -= 1;
  }
  return map_total_vertices;
}

template <typename T>
size_t simplify_generalized_polytopes(GeneralizedPolytope<T> &gp,
                                      std::ostream &os) {
  // Determine all the vertices
  std::unordered_map<MyVector<T>, size_t> map_total_vertices =
      get_map_total_vertices(gp, os);
  size_t n_ext_tot = map_total_vertices.size();
  // Finding the common facet
  struct EntryContain {
    size_t i_polytope;
    int i_fac;
    MyVector<T> eFAC;
  };
  std::unordered_map<Face, std::vector<EntryContain>> map_faces;
  size_t n_polytope_ins = 0;
  auto f_insert_polytope = [&](SinglePolytope<T> const &polytope) -> void {
    size_t i_polytope = n_polytope_ins;
    int n_ext = polytope.EXT.rows();
    int n_fac = polytope.FAC.rows();
    for (int i_fac = 0; i_fac < n_fac; i_fac++) {
      Face fac = polytope.facets[i_fac];
      MyVector<T> eFAC = GetMatrixRow(polytope.FAC, i_fac);
      MyVector<T> fFAC = ScalarCanonicalizationVector(eFAC);
      Face f_big(n_ext_tot);
      for (int i_ext = 0; i_ext < n_ext; i_ext++) {
        if (fac[i_ext] == 1) {
          MyVector<T> eEXT = GetMatrixRow(polytope.EXT, i_ext);
          size_t pos = map_total_vertices.at(eEXT);
          f_big[pos] = 1;
        }
      }
      EntryContain ec{i_polytope, i_fac, fFAC};
      map_faces[f_big].push_back(ec);
    }
    n_polytope_ins += 1;
  };
  auto drop_from_list = [&](size_t i_polytope_delete) -> void {
    std::vector<Face> l_delete;
    for (auto &[face, entries] : map_faces) {
      std::vector<EntryContain> l_ec;
      for (auto &ec : entries) {
        if (ec.i_polytope < i_polytope_delete) {
          l_ec.push_back(ec);
        }
        if (ec.i_polytope > i_polytope_delete) {
          ec.i_polytope -= 1;
          l_ec.push_back(ec);
        }
      }
      entries = l_ec;
      if (entries.empty()) {
        l_delete.push_back(face);
      }
    }
    for (auto &f_del : l_delete) {
      map_faces.erase(f_del);
    }
    gp.polytopes.erase(gp.polytopes.begin() + i_polytope_delete);
    n_polytope_ins -= 1;
  };
  auto test_mergeness = [&](MyMatrix<T> const &EXT1, MyMatrix<T> const &FAC2,
                            int const &i_fac2) -> bool {
    int n_fac2 = FAC2.rows();
    int n_ext1 = EXT1.rows();
    int dim = EXT1.cols();
    for (int j_fac2 = 0; j_fac2 < n_fac2; j_fac2++) {
      if (j_fac2 != i_fac2) {
        for (int i_ext1 = 0; i_ext1 < n_ext1; i_ext1++) {
          T scal(0);
          for (int i = 0; i < dim; i++) {
            scal += EXT1(i_ext1, i) * FAC2(j_fac2, i);
          }
          if (scal < 0) {
            return false;
          }
        }
      }
    }
    return true;
  };
  auto get_merged = [&](MyMatrix<T> const &FAC1, int const &i_fac1,
                        MyMatrix<T> const &FAC2,
                        int const &i_fac2) -> SinglePolytope<T> {
    std::unordered_set<MyVector<T>> set_facet;
    auto insert_block = [&](MyMatrix<T> const &FAC, int const &i_fac) -> void {
      int n_fac = FAC.rows();
      for (int j_fac = 0; j_fac < n_fac; j_fac++) {
        if (i_fac != j_fac) {
          MyVector<T> eFAC = GetMatrixRow(FAC, j_fac);
          MyVector<T> fFAC = ScalarCanonicalizationVector(eFAC);
          set_facet.insert(fFAC);
        }
      }
    };
    insert_block(FAC1, i_fac1);
    insert_block(FAC2, i_fac2);
    std::vector<MyVector<T>> list_facet;
    for (auto &facet : set_facet) {
      list_facet.push_back(facet);
    }
    MyMatrix<T> FACcand = MatrixFromVectorFamily(list_facet);
    return generate_single_polytope(FACcand, os);
  };
  auto look_for_merge = [&]() -> bool {
    for (auto &kv : map_faces) {
      if (kv.second.size() > 2) {
        std::cerr << "GP: That scenario should never happen\n";
        throw TerminalException{1};
      }
      if (kv.second.size() == 2) {
        EntryContain &ec0 = kv.second[0];
        EntryContain &ec1 = kv.second[1];
        MyVector<T> sumFAC = ec0.eFAC + ec1.eFAC;
        if (!IsZeroVector(sumFAC)) {
          std::cerr << "GP: The facet vectors should be opposite\n";
          throw TerminalException{1};
        }
        size_t i_polytope0 = ec0.i_polytope;
        size_t i_polytope1 = ec1.i_polytope;
        int i_fac0 = ec0.i_fac;
        int i_fac1 = ec1.i_fac;
        bool test = test_mergeness(gp.polytopes[i_polytope0].EXT,
                                   gp.polytopes[i_polytope1].FAC, i_fac1);
        if (!test) {
          test = test_mergeness(gp.polytopes[i_polytope1].EXT,
                                gp.polytopes[i_polytope0].FAC, i_fac0);
        }
        if (test) {
          SinglePolytope<T> poly =
              get_merged(gp.polytopes[i_polytope0].FAC, i_fac0,
                         gp.polytopes[i_polytope1].FAC, i_fac1);
          if (i_polytope1 > i_polytope0) {
            drop_from_list(i_polytope1);
            drop_from_list(i_polytope0);
          } else {
            drop_from_list(i_polytope0);
            drop_from_list(i_polytope1);
          }
          f_insert_polytope(poly);
          return true;
        }
      }
    }
    return false;
  };
  auto iterative_merge = [&]() -> void {
    while (true) {
      bool test = look_for_merge();
      if (!test) {
        break;
      }
    }
  };
  // Insert the entries.
  for (auto &poly : gp.polytopes) {
    f_insert_polytope(poly);
  }
  // Simplify them
  iterative_merge();
  return gp.polytopes.size();
}

template <typename T>
std::pair<MyVector<T>, int> get_face_can(MyVector<T> const &eFAC) {
  MyVector<T> fFAC = ScalarCanonicalizationVector(eFAC);
  int dim = fFAC.size();
  for (int i = 0; i < dim; i++) {
    T const &val = fFAC(i);
    if (val != 0) {
      if (val > 0) {
        return {fFAC, 1};
      } else {
        MyVector<T> gFAC = -fFAC;
        return {gFAC, -1};
      }
    }
  }
  std::cerr << "GP: Failed to find a non-zero index\n";
  throw TerminalException{1};
};

template <typename T>
MyMatrix<T> get_fac_subspace(SinglePolytope<T> const& sp, int const &i_fac,
                             MyMatrix<T> const &NSP) {
  MyMatrix<T> const& FAC = sp.FAC;
  int dim = FAC.cols();
  std::vector<int> l_idx_facet = get_adjacent_facet_indices(sp, i_fac);
  int n_fac = l_idx_facet.size();
  MyMatrix<T> FACsub(n_fac, dim - 1);
  int pos = 0;
  for (auto& j_fac: l_idx_facet) {
    MyVector<T> v = GetMatrixRow(FAC, j_fac);
    MyVector<T> eLine = NSP * v;
    for (int i = 0; i < dim - 1; i++) {
      FACsub(pos, i) = eLine(i);
    }
    pos += 1;
  }
  return FACsub;
}

// Two generalized polytopes are adjacent if they share a facet
template <typename T>
std::vector<GeneralizedPolytope<T>>
connected_components_decomposition(GeneralizedPolytope<T> const &gp,
                                   std::ostream &os) {
  int dim = gp.dim;
  struct TrackInfo {
    size_t i_poly;
    int sign;
    MyMatrix<T> FACrel;
  };
  struct AllTrackInfo {
    MyMatrix<T> NSP;
    std::vector<TrackInfo> l_ti;
  };
  std::unordered_map<MyVector<T>, AllTrackInfo> full_track;
  for (size_t i_poly = 0; i_poly < gp.size(); i_poly++) {
    int n_fac = gp.polytopes[i_poly].FAC.rows();
    for (int i_fac = 0; i_fac < n_fac; i_fac++) {
      MyVector<T> eFAC = GetMatrixRow(gp.polytopes[i_poly].FAC, i_fac);
      std::pair<MyVector<T>, int> pair = get_face_can(eFAC);
      AllTrackInfo &rec = full_track[pair.first];
      if (rec.NSP.rows() == 0) {
        rec.NSP = NullspaceMatSingleVectorExt(eFAC);
#ifdef DEBUG_GENERALIZED_POLYTOPE_DISABLE
        os << "GP: connected_components_decomposition, rec.NSP=\n";
        WriteMatrix(os, rec.NSP);
#endif
      }
      MyMatrix<T> FAC = get_fac_subspace(gp.polytopes[i_poly], i_fac, rec.NSP);
      TrackInfo ti{i_poly, pair.second, FAC};
      rec.l_ti.push_back(ti);
    }
  }
  size_t n_polytopes = gp.size();
  GraphBitset GR(n_polytopes);
  for (auto &kv : full_track) {
    std::vector<TrackInfo> &l_ti = kv.second.l_ti;
    size_t n_match = l_ti.size();
    for (size_t i_match = 0; i_match < n_match; i_match++) {
      for (size_t j_match = i_match + 1; j_match < n_match; j_match++) {
        MyMatrix<T> FACtot = Concatenate(l_ti[i_match].FACrel, l_ti[j_match].FACrel);
        bool test = IsFullDimensional(FACtot, os);
        if (test) {
          int sign_prod = l_ti[i_match].sign * l_ti[j_match].sign;
          if (sign_prod == 1) {
            std::cerr << "GP: If that occurs, then it means that we have non-empty "
                         "intersection\n";
            throw TerminalException{1};
          }
          size_t i_poly = l_ti[i_match].i_poly;
          size_t j_poly = l_ti[j_match].i_poly;
          GR.AddAdjacent(i_poly, j_poly);
          GR.AddAdjacent(j_poly, i_poly);
        }
      }
    }
  }
  std::vector<std::vector<size_t>> LConn = ConnectedComponents_set(GR);
  std::vector<GeneralizedPolytope<T>> l_gp;
  for (auto &eConn : LConn) {
    std::vector<SinglePolytope<T>> polytopes;
    for (auto &i_poly : eConn) {
      polytopes.push_back(gp.polytopes[i_poly]);
    }
    GeneralizedPolytope<T> gp_ins{dim, polytopes};
    l_gp.push_back(std::move(gp_ins));
  }
  return l_gp;
}

template <typename T>
GeneralizedPolytope<T>
intersection_gp_gp(GeneralizedPolytope<T> const &gp1,
                   GeneralizedPolytope<T> const &gp2,
                   std::ostream &os) {
  int dim = gp1.dim;
  size_t n_polytope1 = gp1.size();
  size_t n_polytope2 = gp2.size();
  std::vector<SinglePolytope<T>> polytopes;
  for (size_t i1 = 0; i1 < n_polytope1; i1++) {
    for (size_t i2 = 0; i2 < n_polytope2; i2++) {
      MyMatrix<T> const &FAC1 = gp1.polytopes[i1].FAC;
      MyMatrix<T> const &FAC2 = gp2.polytopes[i2].FAC;
      MyMatrix<T> FAC = Concatenate(FAC1, FAC2);
      bool test = IsFullDimensional(FAC, os);
      if (test) {
        polytopes.push_back(generate_single_polytope(FAC, os));
      }
    }
  }
  GeneralizedPolytope<T> gp_ret{dim, polytopes};
  return gp_ret;
}

template <typename T>
bool is_contained_p_vert(SinglePolytope<T> const &p, MyVector<T> const &eEXT,
                         [[maybe_unused]] std::ostream &os) {
  int dim = p.FAC.cols();
  int n_fac = p.FAC.rows();
  for (int i_fac = 0; i_fac < n_fac; i_fac++) {
    T scal(0);
    for (int i = 0; i < dim; i++) {
      scal += p.FAC(i_fac, i) * eEXT(i);
    }
    if (scal < 0) {
      return false;
    }
  }
  return true;
}

template <typename T>
bool is_interior_gp_vert(GeneralizedPolytope<T> const &gp,
                         MyVector<T> const &eEXT, std::ostream &os) {
  for (size_t i_polytope = 0; i_polytope < gp.size(); i_polytope++) {
    bool test = is_contained_p_vert(gp.polytopes[i_polytope], eEXT, os);
    if (test) {
      return true;
    }
  }
  return false;
}

template <typename T>
std::optional<MyVector<T>>
get_interior_point_gp(GeneralizedPolytope<T> const &gp,
                      [[maybe_unused]] std::ostream &os) {
  for (size_t i_polytope = 0; i_polytope < gp.size(); i_polytope++) {
    MyVector<T> eEXT = Isobarycenter(gp.polytopes[i_polytope].EXT);
    return eEXT;
  }
  return {};
}

template <typename T>
std::optional<MyVector<T>>
get_random_interior_point_gp(GeneralizedPolytope<T> const &gp,
                             int const& N,
                             std::ostream &os) {
  for (size_t i_polytope = 0; i_polytope < gp.size(); i_polytope++) {
    MyVector<T> eEXT = random_interior_pt(gp.polytopes[i_polytope].EXT, N, os);
    return eEXT;
  }
  return {};
}

// Tests whether p_sma is contained in p_big.
template <typename T>
bool is_contained_p_p(SinglePolytope<T> const &p_big,
                      SinglePolytope<T> const &p_sma, std::ostream &os) {
  int n_ext = p_sma.EXT.rows();
  for (int i_ext = 0; i_ext < n_ext; i_ext++) {
    MyVector<T> eEXT = GetMatrixRow(p_sma.EXT, i_ext);
    bool test = is_contained_p_vert(p_big, eEXT, os);
    if (!test) {
      return false;
    }
  }
  return true;
}

// Computes p1 - p2.
template <typename T>
GeneralizedPolytope<T> difference_p_p(SinglePolytope<T> const &p1,
                                      SinglePolytope<T> const &p2,
                                      std::ostream &os) {
  std::vector<SinglePolytope<T>> new_polytopes;
  int dim = p1.EXT.cols();
  MyMatrix<T> FACconcat = Concatenate(p1.FAC, p2.FAC);
  if (!IsFullDimensional(FACconcat, os)) {
    // p2 is not contained in p1, so returning just p1.
    new_polytopes.push_back(p1);
    GeneralizedPolytope<T> gp_new{dim, new_polytopes};
    return gp_new;
  }
  if (is_contained_p_p(p2, p1, os)) {
    // p1 totally contained in p2. Returning empty
    GeneralizedPolytope<T> gp_new{dim, new_polytopes};
    return gp_new;
  }
  int n_fac2 = p2.FAC.rows();
  std::vector<MyVector<T>> l_ineq;
  for (int i_fac1 = 0; i_fac1 < p1.FAC.rows(); i_fac1++) {
    MyVector<T> eFAC = GetMatrixRow(p1.FAC, i_fac1);
    l_ineq.push_back(eFAC);
  }
  for (int i_fac2 = 0; i_fac2 < n_fac2; i_fac2++) {
    MyVector<T> eFAC2 = GetMatrixRow(p2.FAC, i_fac2);
    std::vector<MyVector<T>> cand_ineqs = l_ineq;
    cand_ineqs.push_back(-eFAC2);
    MyMatrix<T> FACinput = MatrixFromVectorFamily(cand_ineqs);
    if (IsFullDimensional(FACinput, os)) {
      SinglePolytope<T> sp = generate_single_polytope(FACinput, os);
      new_polytopes.push_back(sp);
    }
    l_ineq.push_back(eFAC2);
  }
  GeneralizedPolytope<T> gp_ret{dim, new_polytopes};
  return gp_ret;
}




template<typename T>
std::vector<ConvexBoundary<T>> convec_boundary_minus_sp(ConvexBoundary<T> const& cb, SinglePolytope<T> const& sp,  std::ostream &os) {
  MyVector<T> eFAC = ScalarCanonicalizationVector(cb.V);
  int i_fac = get_matching_face_position(sp.FAC, eFAC);
  if (i_fac == -1) {
    std::cerr << "GP: The eFAC is not a facet\n";
    throw TerminalException{1};
  }
  std::vector<int> l_idx_facet = get_adjacent_facet_indices(sp, i_fac);
  MyMatrix<T> FAC1 = SelectRow(sp.FAC, l_idx_facet);
  MyMatrix<T> FAC2 = FAC1 * cb.NSP.transpose();
  MyMatrix<T> EXT2 = DirectDualDescription_mat(FAC2, os);
  SinglePolytope<T> sp2 = get_single_polytope(FAC2, EXT2);
  GeneralizedPolytope<T> gp = difference_p_p(cb.sp, sp2, os);
  std::vector<ConvexBoundary<T>> l_cb;
  for (auto & sp_ent: gp.polytopes) {
    ConvexBoundary<T> cb_new{cb.V, cb.NSP, sp_ent};
    l_cb.emplace_back(std::move(cb_new));
  }
  return l_cb;
}




template <typename T>
GeneralizedPolytope<T> difference_gp_p(GeneralizedPolytope<T> const &gp,
                                       SinglePolytope<T> const &p,
                                       std::ostream &os) {
  int dim = gp.dim;
  std::vector<SinglePolytope<T>> new_polytopes;
  for (size_t i = 0; i < gp.size(); i++) {
    GeneralizedPolytope<T> gp_out = difference_p_p(gp.polytopes[i], p, os);
    for (auto &poly : gp_out.polytopes) {
      new_polytopes.push_back(poly);
    }
  }
  GeneralizedPolytope<T> gp_ret{dim, new_polytopes};
  return gp_ret;
}

// Compute the difference gp1 - gp2.
template <typename T>
GeneralizedPolytope<T> difference_gp_gp(GeneralizedPolytope<T> const &gp1,
                                        GeneralizedPolytope<T> const &gp2,
                                        std::ostream &os) {
  GeneralizedPolytope<T> ret_gp = gp1;
  for (size_t i2 = 0; i2 < gp2.size(); i2++) {
    ret_gp = difference_gp_p(ret_gp, gp2.polytopes[i2], os);
  }
  return ret_gp;
}

template <typename T>
bool is_equal(GeneralizedPolytope<T> const &gp1,
              GeneralizedPolytope<T> const &gp2,
              std::ostream &os) {
  GeneralizedPolytope<T> gp_21 = difference_gp_gp(gp2, gp1, os);
  if (!gp_21.empty()) {
    return false;
  }
  GeneralizedPolytope<T> gp_12 = difference_gp_gp(gp1, gp2, os);
  if (!gp_12.empty()) {
    return false;
  }
  return true;
}



template <typename T> struct DataFacetPlusMinus {
  MyMatrix<T> NSP;
  GeneralizedPolytope<T> gp_plus;
  GeneralizedPolytope<T> gp_minus;
};

template <typename T> struct BoundaryGeneralizedPolytope {
  int n;
  std::unordered_map<MyVector<T>, DataFacetPlusMinus<T>> full_data_facets;
  size_t size() const {
    return full_data_facets.size();
  }
  bool is_empty() const {
    return full_data_facets.size() == 0;
  }
};


template <typename T>
void write_generalized_polytope(GeneralizedPolytope<T> const&gp,
                                 std::ostream& os_out) {
  size_t n_p = gp.polytopes.size();
  os_out << "GP: n_p=" << n_p << "\n";
  for (size_t i_p=0; i_p<n_p; i_p++) {
    os_out << "i_p=" << i_p << " EXT=\n";
    WriteMatrix(os_out, gp.polytopes[i_p].EXT);
    os_out << "      FAC=\n";
    WriteMatrix(os_out, gp.polytopes[i_p].FAC);
  }
}





template <typename T>
void write_boundary_generalized_polytope(BoundaryGeneralizedPolytope<T> const& bnd,
                                         std::ostream& os_out) {
  size_t siz = bnd.full_data_facets.size();
  os_out << "n=" << bnd.n << " siz=" << siz << "\n";
  int iter = 0;
  auto f_prt=[&](MyMatrix<T> const& NSP, GeneralizedPolytope<T> const& gp) -> void {
    for (size_t i_p=0; i_p<gp.polytopes.size(); i_p++) {
      MyMatrix<T> EXTw = gp.polytopes[i_p].EXT * NSP;
      os_out << "    EXT" << i_p << "=\n";
      WriteMatrix(os_out, EXTw);
    }
  };
  for (auto &[k, v]: bnd.full_data_facets) {
    os_out << "iter=" << iter << " V=" << StringVectorGAP(k) << "\n";
    os_out << "  gp_minus\n";
    f_prt(v.NSP, v.gp_minus);
    os_out << "  gp_plus\n";
    f_prt(v.NSP, v.gp_plus);
    iter += 1;
  }
}






template <typename T>
BoundaryGeneralizedPolytope<T>
find_generalized_polytope_boundary(GeneralizedPolytope<T> const &gp,
                                   std::ostream &os) {
  int n = gp.polytopes[0].FAC.cols();
  std::unordered_map<MyVector<T>, DataFacetPlusMinus<T>> full_data_facets;
#ifdef DEBUG_GENERALIZED_POLYTOPE_DISABLE
  os << "GP:  find_generalized_polytope_boundary(fgpb) start\n";
#endif
  for (size_t i = 0; i < gp.size(); i++) {
#ifdef DEBUG_GENERALIZED_POLYTOPE_DISABLE
    os << "GP: fgpb, polytope " << i << " EXT=\n";
    WriteMatrix(os, gp.polytopes[i].EXT);
    os << "GP:     FAC=\n";
    WriteMatrix(os, gp.polytopes[i].FAC);
#endif
    int n_fac = gp.polytopes[i].FAC.rows();
    for (int i_fac = 0; i_fac < n_fac; i_fac++) {
      MyVector<T> eFAC = GetMatrixRow(gp.polytopes[i].FAC, i_fac);
      std::pair<MyVector<T>, int> pair = get_face_can(eFAC);
#ifdef DEBUG_GENERALIZED_POLYTOPE
      os << "GP: find_generalized_polytope_boundary, i_fac=" << i_fac << " eFAC=" << StringVectorGAP(eFAC) << "\n";
#endif
      DataFacetPlusMinus<T> &rec = full_data_facets[pair.first];
      if (rec.NSP.rows() == 0) {
        rec.NSP = NullspaceMatSingleVectorExt(eFAC);
#ifdef DEBUG_GENERALIZED_POLYTOPE_DISABLE
        os << "GP: find_generalized_polytope_boundary, rec.NSP=\n";
        WriteMatrix(os, rec.NSP);
#endif
      }
      MyMatrix<T> FAC = get_fac_subspace(gp.polytopes[i], i_fac, rec.NSP);
#ifdef DEBUG_GENERALIZED_POLYTOPE
      os << "GP: fgpb, polytope " << i << " i_fac=" << i_fac << " FAC=\n";
      WriteMatrix(os, FAC);
#endif
      SinglePolytope<T> sp = generate_single_polytope(FAC, os);
#ifdef DEBUG_GENERALIZED_POLYTOPE
      os << "GP: fgpb, polytope " << i << " i_fac=" << i_fac << " EXT=\n";
      WriteMatrix(os, sp.EXT);
#endif
      if (pair.second == 1) {
        rec.gp_plus.polytopes.push_back(sp);
      } else {
        rec.gp_minus.polytopes.push_back(sp);
      }
    }
  }
  std::vector<MyVector<T>> to_remove;
  for (auto &kv : full_data_facets) {
    GeneralizedPolytope<T> diff_p_m =
        difference_gp_gp(kv.second.gp_plus, kv.second.gp_minus, os);
    GeneralizedPolytope<T> diff_m_p =
        difference_gp_gp(kv.second.gp_minus, kv.second.gp_plus, os);
    if (diff_p_m.empty() && diff_m_p.empty()) {
      to_remove.push_back(kv.first);
    } else {
      kv.second.gp_plus = diff_p_m;
      kv.second.gp_minus = diff_m_p;
    }
  }
  for (auto &vect : to_remove) {
    full_data_facets.erase(vect);
  }
  return {n, full_data_facets};
}


template<typename T>
struct InteriorPtDir {
  MyVector<T> pt;
  MyVector<T> FacIneq;
  MyVector<T> get_point(T const& shift) const {
    int dim = pt.size();
    MyVector<T> V(dim);
    V(0) = pt(0);
    for (int u=1; u<dim; u++) {
      V(u) = pt(u) + shift * FacIneq(u);
    }
    return V;
  }
  std::string to_string() const {
    std::string ret = " pt=" + StringVectorGAP(pt) + " FacIneq=" + StringVectorGAP(FacIneq);
    return ret;
  }
};




template<typename T>
InteriorPtDir<T> ipd_opposite(InteriorPtDir<T> const& ipd) {
  MyVector<T> FacIneq = - ipd.FacIneq;
  return {ipd.pt, FacIneq};
}


template <typename T>
std::optional<InteriorPtDir<T>>
get_random_interior_point_bnd(BoundaryGeneralizedPolytope<T> const &bnd,
                              int const& N,
                              std::ostream &os) {
  for (auto &kv : bnd.full_data_facets) {
    auto get_opt = [&]() -> std::optional<InteriorPtDir<T>> {
      std::optional<MyVector<T>> opt1 = get_random_interior_point_gp(kv.second.gp_plus, N, os);
      if (opt1) {
        MyVector<T> const& FacIneq = kv.first;
        InteriorPtDir<T> ipd{*opt1, FacIneq};
        return ipd;
      }
      std::optional<MyVector<T>> opt2 = get_random_interior_point_gp(kv.second.gp_minus, N, os);
      if (opt2) {
        MyVector<T> FacIneq = - kv.first;
        InteriorPtDir<T> ipd{*opt2, std::move(FacIneq)};
        return ipd;
      }
      std::cerr << "GP: the size should be non-zero. Since otherwise, it "
                   "should be removed from the list\n";
      throw TerminalException{1};
    };
    std::optional<InteriorPtDir<T>> opt = get_opt();
    if (opt) {
      InteriorPtDir<T> const &sol_A = *opt;
      MyVector<T> eIso_B = kv.second.NSP.transpose() * sol_A.pt;
      T val = eIso_B(0);
      for (int u=0; u<eIso_B.size(); u++) {
        eIso_B(u) = eIso_B(u) / val;
      }
      InteriorPtDir<T> sol_B{eIso_B, sol_A.FacIneq};
      return sol_B;
    }
  }
  return {};
}



template <typename T>
std::optional<InteriorPtDir<T>>
get_interior_point_bnd(BoundaryGeneralizedPolytope<T> const &bnd,
                       std::ostream &os) {
  for (auto &kv : bnd.full_data_facets) {
    auto get_opt = [&]() -> std::optional<InteriorPtDir<T>> {
      std::optional<MyVector<T>> opt1 =
          get_interior_point_gp(kv.second.gp_plus, os);
      if (opt1) {
        MyVector<T> const& FacIneq = kv.first;
        InteriorPtDir<T> ipd{*opt1, FacIneq};
        return ipd;
      }
      std::optional<MyVector<T>> opt2 =
          get_interior_point_gp(kv.second.gp_minus, os);
      if (opt2) {
        MyVector<T> FacIneq = - kv.first;
        InteriorPtDir<T> ipd{*opt2, std::move(FacIneq)};
        return ipd;
      }
      std::cerr << "GP: the size should be non-zero. Since otherwise, it "
                   "should be removed from the list\n";
      throw TerminalException{1};
    };
    std::optional<InteriorPtDir<T>> opt = get_opt();
    if (opt) {
      InteriorPtDir<T> const &sol_A = *opt;
      MyVector<T> eIso_B = kv.second.NSP.transpose() * sol_A.pt;
      T val = eIso_B(0);
      for (int u=0; u<eIso_B.size(); u++) {
        eIso_B(u) = eIso_B(u) / val;
      }
      InteriorPtDir<T> sol_B{eIso_B, sol_A.FacIneq};
      return sol_B;
    }
  }
  return {};
}

template <typename T>
bool is_boundary_point(InteriorPtDir<T> const& ipd, BoundaryGeneralizedPolytope<T> const &bnd, std::ostream &os) {
  MyVector<T> V = ScalarCanonicalizationVector(ipd.FacIneq);
  std::pair<MyVector<T>, int> pair = get_face_can(V);
  if (!bnd.full_data_facets.contains(pair.first)) {
#ifdef DEBUG_GENERALIZED_POLYTOPE
    os << "GP: is_boundary_point, false case 1\n";
#endif
    return false;
  }
  DataFacetPlusMinus<T> const& dfpm = bnd.full_data_facets.at(pair.first);
  std::optional<MyVector<T>> opt = SolutionMat(dfpm.NSP, ipd.pt);
  if (!opt) {
#ifdef DEBUG_GENERALIZED_POLYTOPE
    os << "GP: is_boundary_point, false case 2\n";
#endif
    return false;
  }
  MyVector<T> const& pt = *opt;
#ifdef DEBUG_GENERALIZED_POLYTOPE
  os << "GP: is_boundary_point, determining via the is_interior_gp_vert\n";
#endif
  if (pair.second == 1) {
    return is_interior_gp_vert(dfpm.gp_plus, pt, os);
  } else {
    return is_interior_gp_vert(dfpm.gp_minus, pt, os);
  }
}




template <typename T>
void reduce_boundary_generalized_polytope(BoundaryGeneralizedPolytope<T> &bnd,
                                          GeneralizedPolytope<T> const &gp,
                                          std::ostream &os) {
  for (size_t i = 0; i < gp.size(); i++) {
    int n_fac = gp.polytopes[i].FAC.rows();
    for (int i_fac = 0; i_fac < n_fac; i_fac++) {
      MyVector<T> eFAC = GetMatrixRow(gp.polytopes[i].FAC, i_fac);
      std::pair<MyVector<T>, int> pair = get_face_can(eFAC);
      if (bnd.full_data_facets.contains(pair.first)) {
        DataFacetPlusMinus<T> &df = bnd.full_data_facets[pair.first];
        MyMatrix<T> FAC = get_fac_subspace(gp.polytopes[i], i_fac, df.NSP);
        SinglePolytope<T> sp = generate_single_polytope(FAC, os);
        if (pair.second == 1) {
          df.gp_minus = difference_gp_p(df.gp_minus, sp, os);
        } else {
          df.gp_plus = difference_gp_p(df.gp_plus, sp, os);
        }
      }
    }
  }
  std::vector<MyVector<T>> to_remove;
  for (auto &kv : bnd.full_data_facets) {
    if (kv.second.gp_plus.empty() && kv.second.gp_minus.empty()) {
      to_remove.push_back(kv.first);
    }
  }
  for (auto &vect : to_remove) {
    bnd.full_data_facets.erase(vect);
  }
}

template <typename T>
MyMatrix<T> get_vertices_gp_bnd(GeneralizedPolytope<T> const &gp,
                                BoundaryGeneralizedPolytope<T> const &bnd,
                                std::ostream &os) {
  if (gp.empty()) {
    return {};
  }
  int dim = gp.dim;
  std::unordered_set<MyVector<T>> set_vertices;
  for (size_t i = 0; i < gp.size(); i++) {
    int n_ext = gp.polytopes[i].EXT.rows();
    for (int i_ext = 0; i_ext < n_ext; i_ext++) {
      MyVector<T> eEXT = GetMatrixRow(gp.polytopes[i].EXT, i_ext);
      set_vertices.insert(eEXT);
    }
  }
#ifdef DEBUG_GENERALIZED_POLYTOPE
  os << "GP: get_vertices_gp_bnd, |gp|=" << gp.size() << " |set_vertices|=" << set_vertices.size() << "\n";
  write_generalized_polytope(gp, os);
  write_boundary_generalized_polytope(bnd, os);
#endif
  std::vector<MyVector<T>> l_vertices;
  for (auto &eEXT : set_vertices) {
    std::vector<MyVector<T>> l_fac;
    for (auto &kv : bnd.full_data_facets) {
      MyVector<T> const &eFAC = kv.first;
      T scal = eFAC.dot(eEXT);
      if (scal == 0) {
        std::optional<MyVector<T>> opt = SolutionMat(kv.second.NSP, eEXT);
        if (!opt) {
          std::cerr << "GP: It should be in the subspace\n";
          throw TerminalException{1};
        }
        MyVector<T> const &eEXTred = *opt;
        bool test1 = is_interior_gp_vert(kv.second.gp_minus, eEXTred, os);
        bool test2 = is_interior_gp_vert(kv.second.gp_plus, eEXTred, os);
        if (test1 || test2) {
          l_fac.push_back(eFAC);
        }
      }
    }
    int rnk = RankMat(MatrixFromVectorFamilyDim(dim, l_fac));
    if (rnk == dim - 1) {
      l_vertices.push_back(eEXT);
    }
  }
  return MatrixFromVectorFamilyDim(dim, l_vertices);
}

template <typename T>
MyMatrix<T> get_vertices_gp(GeneralizedPolytope<T> const &gp, std::ostream &os) {
  BoundaryGeneralizedPolytope<T> bnd = find_generalized_polytope_boundary(gp, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: convert_p_voronoi_part, step 3\n";
#endif
  return get_vertices_gp_bnd(gp, bnd, os);
}

template <typename T>
MyMatrix<T> get_vertices_bnd(BoundaryGeneralizedPolytope<T> const& bnd, std::ostream &os) {
  int dim = bnd.n;
  std::unordered_set<MyVector<T>> set_vert;
  for (auto & kv: bnd.full_data_facets) {
    auto f_insert=[&](GeneralizedPolytope<T> const& gp) -> void {
      MyMatrix<T> EXT1 = get_vertices_gp(gp, os);
      MyMatrix<T> EXT2 = EXT1 * kv.second.NSP;
      int n_ext = EXT2.rows();
      for (int i_ext=0; i_ext<n_ext; i_ext++) {
        MyVector<T> V = GetMatrixRow(EXT2, i_ext);
        set_vert.insert(V);
      }
    };
    f_insert(kv.second.gp_plus);
    f_insert(kv.second.gp_minus);
  }
  std::vector<MyVector<T>> l_vert(set_vert.begin(), set_vert.end());
  return MatrixFromVectorFamilyDim(dim, l_vert);
}







template <typename T>
size_t ComputeHash_gp(GeneralizedPolytope<T> const &gp, size_t seed_in, std::ostream& os) {
  MyMatrix<T> M1 = get_vertices_gp(gp, os);
  MyMatrix<T> M2 = reorder_matrix(M1);
  return Matrix_Hash(M2, seed_in);
}


template<typename T>
std::vector<GeneralizedPolytope<T>> get_p_voronoi_orbit(std::vector<MyMatrix<T>> const& LGen,
                                                        GeneralizedPolytope<T> const& gp,
                                                        std::ostream &os) {
  auto f_prod=[&](GeneralizedPolytope<T> const& gp,
                  MyMatrix<T> const& P) -> GeneralizedPolytope<T> {
    return mat_product(gp, P);
  };
  auto f_hash=[&](GeneralizedPolytope<T> const& gp) -> size_t {
    size_t seed_in = 34;
    return ComputeHash_gp(gp, seed_in, os);
  };
  auto f_equal=[&](GeneralizedPolytope<T> const& gp1, GeneralizedPolytope<T> const& gp2) -> bool {
    return is_equal(gp1, gp2, os);
  };
  return OrbitComputationGen(LGen, gp, f_prod, f_hash, f_equal, os);
}








template<typename T>
T volume_gp(GeneralizedPolytope<T> const& gp,  [[maybe_unused]] std::ostream& os) {
  T volume(0);
  for (auto & sp: gp.polytopes) {
    T vol = lrs::Kernel_VolumePolytope(sp.EXT);
    volume += vol;
  }
  return volume;
}

template <typename T>
void WriteEntryCPP(std::ostream &os, SinglePolytope<T> const &sp) {
  WriteMatrix(os, sp.EXT);
  WriteMatrix(os, sp.FAC);
  WriteListFace(os, sp.facets);
}

template <typename T>
SinglePolytope<T> ReadEntryCPP_SinglePolytope(std::istream &is) {
  MyMatrix<T> EXT = ReadMatrix<T>(is);
  MyMatrix<T> FAC = ReadMatrix<T>(is);
  vectface facets = ReadListFace(is);
  return SinglePolytope<T>(EXT, FAC, facets);
}

template <typename T>
void WriteEntryCPP(std::ostream &os, ConvexBoundary<T> const &cb) {
  WriteVector(os, cb.V);
  WriteMatrix(os, cb.NSP);
  WriteEntryCPP(os, cb.sp);
}

template <typename T>
ConvexBoundary<T> ReadEntryCPP_ConvexBoundary(std::istream &is) {
  MyVector<T> V = ReadVector<T>(is);
  MyMatrix<T> NSP = ReadMatrix<T>(is);
  SinglePolytope<T> sp = ReadEntryCPP_SinglePolytope<T>(is);
  return {V, NSP, sp};
}

template <typename T>
void WriteEntryCPP(std::ostream &os, GeneralizedPolytope<T> const &gp) {
  os << gp.dim << "\n";
  size_t n = gp.polytopes.size();
  os << n << "\n";
  for (size_t i = 0; i < n; i++) {
    WriteEntryCPP(os, gp.polytopes[i]);
  }
}

template <typename T>
GeneralizedPolytope<T> ReadEntryCPP_GeneralizedPolytope(std::istream &is) {
  int dim;
  is >> dim;
  size_t n;
  is >> n;
  std::vector<SinglePolytope<T>> polytopes;
  for (size_t i = 0; i < n; i++) {
    polytopes.push_back(ReadEntryCPP_SinglePolytope<T>(is));
  }
  return {dim, polytopes};
}

namespace boost::serialization {

template <class Archive, typename T>
inline void serialize(Archive &ar, SinglePolytope<T> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("EXT", val.EXT);
  ar &make_nvp("FAC", val.FAC);
  ar &make_nvp("facets", val.facets);
}

template <class Archive, typename T>
inline void serialize(Archive &ar, ConvexBoundary<T> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("V", val.V);
  ar &make_nvp("NSP", val.NSP);
  ar &make_nvp("sp", val.sp);
}

template <class Archive, typename T>
inline void serialize(Archive &ar, GeneralizedPolytope<T> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("polytopes", val.polytopes);
}

} // namespace boost::serialization

// clang-format off
#endif  // SRC_ROBUST_COVERING_GENERALIZED_POLYTOPES_H_
// clang-format on

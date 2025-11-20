// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_GENERALIZED_POLYTOPES_H_
#define SRC_ROBUST_COVERING_GENERALIZED_POLYTOPES_H_

// clang-format off
#include "POLY_LinearProgramming.h"
#include "POLY_DirectDualDesc.h"
#include "POLY_Fundamental.h"
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


template<typename T>
struct SinglePolytope {
  MyMatrix<T> EXT;
  MyMatrix<T> FAC;
  vectface facets;
}

template<typename T>
SinglePolytope<T> generate_single_polytope(MyMatrix<T> const& FACinput, std::ostream& os) {
  std::vector<int> ListIrred = cdd::RedundancyReductionClarkson(FACinput, os);
  MyMatrix<T> FAC = SelectRow(FACcand, ListIrred);
  MyMatrix<T> EXT = DirectDualDescription(FAC, os);
  int n_fac = FAC.rows();
  int n_ext = EXT.rows();
  int dim = EXT.cols();
  vectface facets(n_ext);
  for (int i_fac=0; i_fac<n_fac; i_fac++) {
    Face f(n_ext);
    for (int i_ext=0; i_ext<n_ext; i_ext++) {
      T scal(0);
      for (int i=0; i<dim; i++) {
        scal += FAC(i_fac, i) * EXT(i_ext, i);
      }
      if (scal == 0) {
        f[i_ext] = 1;
      }
    }
    facets.push_back(f);
  }
  return {std::move(EXT), std::move(FAC), std::move(facets)};
}

template<typename T>
struct GeneralizedPolytope {
  SinglePolytope<T> polytopes;
}

template<typename T>
bool check_pairwise_intersection(GeneralizedPolytope<T> const& gp, std::ostream& os) {
  size_t n_comp = gp.polytopes.size();
  for (size_t i=0; i<n_comp; i++) {
    for (size_t j=i+1; j<n_comp; j++) {
      MyMatrix<T> FAC = Concatenate(gp.polytopes[i].FAC, gp.polytopes[j].FAC);
      bool test = IsFullDimensional(FAC, os);
      if (test) {
        return false;
      }
    }
  }
  return true;
}

template<typename T>
std::unordered_map<MyVector<T>, size_t> get_map_total_vertices(GeneralizedPolytope<T> const& gp, std::ostream& os) {
  std::unordered_map<MyVector<T>, size_t> map_total_vertices;
  size_t pos = 1;
  for (auto & polytope: gp.polytopes) {
    int n_ext = polytope.EXT.rows();
    for (int i_ext=0; i_ext<n_ext; i_ext++) {
      MyVector<T> eEXT = GetMatrixRow(polytope.EXT, i_ext);
      size_t& find_pos = map_total_vertices[eEXT];
      if (find_pos == 0) {
        find_pos = pos;
        pos += 1;
      }
    }
  }
  for (auto &kv: map_total_vertices) {
    kv.second -= 1;
  }
  return map_total_vertices;
}



template<typename T>
size_t simplify_generalized_polytopes(GeneralizedPolytope<T> const& gp, std::ostream& os) {
  size_t n_polytope = gp.polytopes.size();
  // Determine all the vertices
  std::unordered_map<MyVector<T>, size_t> map_total_vertices = get_map_total_vertices(gp, os);
  size_t n_ext_tot = map_total_vertices.size();
  // Finding the common facet
  struct EntryContain {
    size_t i_polytope;
    int i_fac;
    MyVector<T> eFAC;
  };
  std::unordered_map<Face, std::vector<EntryContain>> map_faces;
  size_t n_polytope_ins = 0;
  auto f_insert_polytope=[&](SinglePolytope<T> const& polytope) -> void {
    size_t i_polytope = n_polytope_ins;
    int n_ext = polytope.EXT.rows();
    int n_fac = polytope.FAC.rows();
    for (int i_fac=0; i_ext<n_ext; i_ext++) {
      Face fac = polytope.facets[i_fac];
      MyVector<T> eFAC = GetMatrixRow(polytope.FAC, i_fac);
      MyVector<T> fFAC = ScalarCanonicalizationVector(eFAC);
      Face f_big(n_ext_tot);
      for (int i_ext=0; i_ext<n_ext; i_ext++) {
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
  auto drop_from_list=[&](size_t i_polytope_delete) -> void {
    std::vector<Face> l_delete;
    for (auto &kv: map_faces) {
      std::vector<EntryContain> l_ec;
      for (auto &ec: kv.second) {
        if (ec.i_polytope < i_polytope_delete) {
          l_ec.push_back(ec);
        }
        if (ec.i_polytope > i_polytope_delete) {
          ec.i_polytope -= 1;
          l_ec.push_back(ec);
        }
      }
      kv.second = l_ec;
      if (kv.second.size() == 0) {
        l_delete.push_back(kv.first);
      }
    }
    for (auto & f_del: l_delete) {
      map_faces.delete(f_del);
    }
    gp.polytopes.remove(i_polytope_delete);
    n_polytope_ins -= 1;
  };
  auto test_mergeness=[&](MyMatrix<T> const& EXT1, MyMatrix<T> const& FAC2, int const& i_fac2) -> bool {
    int n_fac2 = FAC2.rows();
    int n_ext1 = EXT1.rows();
    int dim = EXT1.cols();
    for (int j_fac2=0; j_fac2<n_fac2; j_fac2++) {
      if (j_fac2 != i_fac2) {
        for (int i_ext1=0; i_ext1<n_ext1; i_ext1++) {
          T scal(0);
          for (int i=0; i<dim; i++) {
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
  auto get_merged=[&](MyMatrix<T> const& FAC1, int const& i_fac1, MyMatrix<T> const& FAC2, int const& i_fac2) -> SinglePolytope<T> {
    std::unordered_set<MyVector<T>> set_facet;
    auto insert_block=[&](MyMatrix<T> const& FAC, int const& i_fac) -> void {
      int n_fac = FAC.rows();
      for (int j_fac=0; j_fac<n_fac; j_fac++) {
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
    for (auto & facet: set_facet) {
      list_facet.push_back(facet);
    }
    MyMatrix<T> FACcand = MatrixFromVectorFamily(list_facet);
    return generate_single_polytope(FACcand, os);
  };
  auto look_for_merge=[&]() -> bool {
    size_t n_polytope = n_polytope_ins;
    for (auto & kv: map_faces) {
      if (kv.second.size() > 2) {
        std::cerr << "GP: That scenario should never happen\n";
        throw TerminalException{1};
      }
      if (kv.second.size() == 2) {
        EntryContain& ec0 = kv.second[0];
        EntryContain& ec1 = kv.second[1];
        MyVector<T> sumFAC = ec0.eFAC + ec1.eFAC;
        if (!IsZeroVector(sumFAC)) {
          std::cerr << "GP: The facet vectors should be opposite\n";
          throw TerminalException{1};
        }
        size_t i_polytope0 = ec0.i_polytope;
        size_t i_polytope1 = ec1.i_polytope;
        int i_fac0 = ec0.i_fac;
        int i_fac1 = ec1.i_fac;
        bool test = test_mergeness(gp.polytopes[i_polytope0].EXT, gp.polytopes[i_polytope1].FAC, i_fac1);
        if (!test) {
          test = test_mergeness(gp.polytopes[i_polytope1].EXT, gp.polytopes[i_polytope0].FAC, i_fac0);
        }
        if (test) {
          SinglePolytope<T> poly = get_merged(gp.polytopes[i_polytope0].FAC, i_fac0, gp.polytopes[i_polytope1].FAC, i_fac1);
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
  auto iterative_merge=[&]() -> void {
    while(true) {
      bool test = look_for_merge();
      if (!test) {
        break;
      }
    }
  };
  // Insert the entries.
  for (auto & poly: gp.polytopes) {
    f_insert_polytope(poly);
  }
  // Simplify them
  iterative_merge();
}


template<typename T>
std::pair<MyVector<T>, int> get_face_can(MyVector<T> const& eFAC) {
  MyVector<T> fFAC = ScalarCanonicalizationVector(eFAC);
  int dim = fFAC.size();
  for (int i=0; i<dim; i++) {
    T const& val = fFAC(i);
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

template<typename T>
MyMatrix<T> get_fac_subspace(MyMatrix<T> const& FAC, int const& i_fac, MyMatrix<T> const& NSP) {
  int n_fac = FAC.rows();
  int dim = FAC.cols();
  MyMatrix<T> FACsub(n_fac - 1, dim-1);
  int pos = 0;
  for (int j_fac=0; j_fac<n_fac; j_fac++) {
    if (i_fac != j_fac) {
      MyVector<T> v = GetMatrixRow(FAC, j_fac);
      MyVector<T> eLine = rec.NSP * v;
      for (int i=0; i<dim-1; i++) {
        FACsub(pos, i) = eLine(i);
      }
      pos += 1;
    }
  }
  return FACsub;
}



// Two generalized polytopes are adjacent if they share a facet
template<typename T>
std::vector<GeneralizedPolytope<T>> connected_components_decomposition(GeneralizedPolytope<T> const& gp, std::ostream& os) {
  int dim = gp.polytopes[0].FAC.cols();
  //  std::unordered_map<MyVector<T>, size_t> map_total_vertices = get_map_total_vertices(gp, os);

  struct TrackInfo {
    size_t i_polytope;
    int sign;
    MyMatrix<T> FACrel;
  };
  struct AllTrackInfo {
    MyMatrix<T> NSP;
    std::vector<TrackInfo> l_ti;
  };
  std::unordered_map<MyVector<T>, AllTrackInfo> full_track;
  for (size_t i_polytope=0; i_polytope<gp.polytopes.size(); i_polytope++) {
    int n_fac = gp.polytopes[i_polytope].FAC.rows();
    for (int i_fac=0; i_fac<n_fac; i_fac++) {
      MyVector<T> eFAC = GetMatrixRow(gp.polytopes[i_polytope].FAC, i_fac);
      std::pair<MyVector<T>, int> pair = get_face_can(eFAC);
      AllTrackInfo& rec = full_track[pair.first];
      if (rec.NSP.rows() == 0) {
        rec.NSP = NullspaceVector(eFAC);
      }
      MyMatrix<T> FAC = get_fac_subspace(gp.polytopes[i_polytope].FAC, i_fac, rec.NSP);
      TrackInfo ti{i_polytope, pair.second, FAC};
      rec.l_ti.push_back(ti);
    }
  }
  GraphBitset GR(n_polytopes);
  for (auto & kv: full_track) {
    std::vector<TrackInfo> & l_ti = kv.second.l_ti;
    size_t n_match = l_ti.size();
    for (size_t i_match=0; i_match<n_match; i_match++) {
      for (size_t j_match=i_match+1; j_match<n_match; j_match++) {
        MyMatrix<T> FACtot = Concatenate(l_ti[i_match].FAC, l_ti[j_match].FAC);
        bool test = IsFullDimensional(FACtot, os);
        if (test) {
          int sign_prod = l_ti[i_match].sign * l_ti[j_match].sign;
          if (sign_prod == 1) {
            std::cerr << "GP: If that occurs, then it means that we have empty intersection\n";
            throw TerminalException{1};
          }
          size_t i_poly = l_ti[i_match].i_polytope;
          size_t j_poly = l_ti[j_match].i_polytope;
          GR.AddAdjacent(i_poly, j_poly);
          GR.AddAdjacent(j_poly, i_poly);
        }
      }
    }
  }
  std::vector<std::vector<size_t>> LConn = ConnectedComponents_set(GR);
  std::vector<GeneralizedPolytope<T>> l_gp;
  for (auto & eConn: LConn) {
    GeneralizedPolytope<T> gp_ins;
    for (auto & i_polytope: eConn) {
      gp_ins.polytopes.push_back(gp.polytopes[i_polytope]);
    }
  }
  return l_gp;
}

template<typename T>
GeneralizedPolytope<T> intersection_generalized_polytope(GeneralizedPolytope<T> const& gp1, GeneralizedPolytope<T> const& gp2, std::ostream& os) {
  size_t n_polytope1 = gp1.polytopes.size();
  size_t n_polytope2 = gp1.polytopes.size();
  std::vector<SinglePolytope<T>> polytopes;
  for (size_t i1=0; i1<n_polytope1; i1++) {
    for (size_t i2=0; i2<n_polytope2; i2++) {
      MyMatrix<T> const& FAC1 = gp1.polytopes[i1].FAC;
      MyMatrix<T> const& FAC2 = gp1.polytopes[i2].FAC;
      MyMatrix<T> FAC = Concatenate(FAC1, FAC2);
      bool test = IsFullDimensional(FAC, os);
      if (test) {
        polytopes.push_back(generate_single_polytope(FACcand, os));
      }
    }
  }
  return {polytopes};
}

// Computes p1 - p2.
template<typename T>
GeneralizedPolytope<T> difference_p_p(SinglePolytope<T> const& p1, SinglePolytope<T> const& p2, std::ostream& os) {
  std::vector<SinglePolytope<T>> new_polytopes;
  int n_fac2 = p2.FAC.rows();
  std::vector<MyVector<T>> l_ineq;
  for (int i_fac1=0; i_fac1<p1.FAC.rows(); i_fac1++) {
    MyVector<T> eFAC = GetMatrixRow(p1.FAC, i_fac1);
    l_ineq.push_back(eFAC);
  }
  for (int i_fac2=0; i_fac2<n_fac2; i_fac2++) {
    MyVector<T> eFAC2 = GetMatrixRow(p2.FAC, i_fac2);
    std::vector<MyVector<T>> cand_ineqs = l_ineq;
    cand_ineqs.push_back(-eFAC2);
    MyMatrix<T> FACinput = MatrixFromVectorFamily(cand_ineqs);
    if (IsFullDimensional(FACinput)) {
      SinglePolytope<T> sp = generate_single_polytope(FACinput, os);
      new_polytopes.push_back(sp);
    }
    l_ineq.push_back(eFAC2);
  }
  return {new_polytopes};
}

template<typename T>
GeneralizedPolytope<T> difference_gp_p(GeneralizedPolytope<T> const& gp, SinglePolytope<T> const& p, std::ostream& os) {
  std::vector<SinglePolytope<T>> new_polytopes;
  for (size_t i=0; i<gp.polytopes.size(); i++) {
    GeneralizedPolytope<T> gp = difference_p_p(gp.polytopes[i], p, os);
    for (auto & poly: gp.polytopes) {
      new_polytopes.push_back(poly);
    }
  }
  return {new_polytopes};
}


// Compute the difference gp1 - gp2.
template<typename T>
GeneralizedPolytope<T> difference_gp_gp(GeneralizedPolytope<T> const& gp1, GeneralizedPolytope<T> const& gp2, std::ostream& os) {
  GeneralizedPolytope<T> ret_gp = gp1;
  for (size_t i2=0; i2<gp2.polytope.size(); i2++) {
    ret_gp = difference_gp_p(ret_gp, gp2.polytopes[i2], os);
  }
  return ret_gp;
}

template<typename T>
struct DataFacet {
  MyMatrix<T> NSP;
  GeneralizedPolytope<T> gp_plus;
  GeneralizedPolytope<T> gp_minus;
};

template<typename T>
struct BoundaryGeneralizedPolytope {
  std::unordered_map<MyVector<T>, DataFacet<T>> full_data_facets;
};

template<typename T>
ListFacetGeneralizedPolytope<T> find_generalized_polytope_boundary(GeneralizedPolytope<T> const& gp, std::ostream& os) {
  std::unordered_map<MyVector<T>, DataFacet> full_data_facets;
  for (size_t i=0; i<gp.polytopes.size(); i++) {
    int n_fac = gp.polytopes[i].FAC.rows();
    for (int i_fac=0; i_fac<n_fac; i_fac++) {
      MyVector<T> eFAC = GetMatrixRow(gp.polytopes[i].FAC, i_fac);
      std::pair<MyVector<T>, int> pair = get_face_can(eFAC);
      DataFacet& rec = full_data_facets[pair.first];
      if (rec.NSP.rows() == 0) {
        rec.NSP = NullspaceVec(eFAC);
      }
      MyMatrix<T> FAC = get_fac_subspace(gp.polytopes[i].FAC, i_fac, rec.NSP);
      SinglePolytope<T> sp = generate_single_polytope(FAC, os);
      if (pair.second == 1) {
        gp_plus.polytopes.push_back(sp);
      } else {
        gp_minus.polytopes.push_back(sp);
      }
    }
  }
  ListFacetGeneralizedPolytope<T> return_val;
  for (auto & kv: full_data_facets) {
    GeneralizedPolytope<T> diff_p_m = difference_gp_gp(kv.second.gp_plus, kv.second.gp_minus, os);
    GeneralizedPolytope<T> diff_m_p = difference_gp_gp(kv.second.gp_minus, kv.second.gp_plus, os);
    kv.second.gp_plus = diff_p_m;
    kv.second.gp_minus = diff_m_p;
  }
  return {full_data_facets};
}

template<typename T>
void reduce_boundary_generalized_polytope(BoundaryGeneralizedPolytope<T> & bnd, GeneralizedPolytope<T> const& gp, std::ostream& os) {
  for (size_t i=0; i<gp.polytopes.size(); i++) {
    int n_fac = gp.polytopes[i].FAC.rows();
    for (int i_fac=0; i_fac<n_fac; i_fac++) {
      MyVector<T> eFAC = GetMatrixRow(gp.polytopes[i].FAC, i_fac);
      std::pair<MyVector<T>, int> pair = get_face_can(eFAC);
      if (bnd.full_data_facets.count(pair.first) == 1) {
        DataFacet<T>& df = bnd.full_data_facets[pair.first];
        MyMatrix<T> FAC = get_fac_subspace(gp.polytopes[i].FAC, i_fac, df.NSP);
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
  for (auto & kv: bnd.full_data_facets) {
    if (kv.second.gp_plus.size() == 0 && kv.second.gp_minus.size() == 0) {
      to_remove.push_back(kv.first);
    }
  }
  for (auto & vect: to_remove) {
    bnd.full_data_facets.remove(vect);
  }
}






// clang-format off
#endif  // SRC_ROBUST_COVERING_GENERALIZED_POLYTOPES_H_
// clang-format on

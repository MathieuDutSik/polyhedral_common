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
      or just a vertex
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
size_t simplify_generalized_polytopes(GeneralizedPolytope<T> const& gp, std::ostream& os) {
  size_t n_polytope = gp.polytopes.size();
  // Determine all the vertices
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
          size_t pos = map_total_vertices.at(eEXT) - 1;
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

// Compute the difference gp1 - gp2.
template<typename T>
GeneralizedPolytope<T> difference_generalized_polytope(GeneralizedPolytope<T> const& gp1, GeneralizedPolytope<T> const& gp2, std::ostream& os) {
  std::vector<SinglePolytope<T>> polytopes = gp1.polytopes;
  size_t n_polytope2 = gp2.polytopes.size();
  for (size_t i2=0; i2<n_polytope2; i2++) {
    size_t n_polytope1 = gp1.polytopes.size();
  }
  for (size_t i1=0; i1<n_polytope1; i1++) {
  }


}







// clang-format off
#endif  // SRC_ROBUST_COVERING_GENERALIZED_POLYTOPES_H_
// clang-format on

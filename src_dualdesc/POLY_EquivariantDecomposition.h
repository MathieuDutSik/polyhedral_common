// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_EQUIVARIANTTRIANGULATION_H_
#define SRC_DUALDESC_POLY_EQUIVARIANTTRIANGULATION_H_

// clang-format off
#include "POLY_RecursiveDualDesc.h"
#include <unordered_map>
#include <vector>
// clang-format on

// The computation of several things require using triangulations.
//
// We have essentially two methods for computing them:
// * The LRS can produce a triangulation.
// * For a symmetric polytope, compute the barycenter.
//   Then from the list of orbit of facets, compute the
//   decomposition from the isobarycenter.
//
// We cannot compute an equivariant triangulation in a
// simple way but fortunately, we do not really neeed that.
// We do not store the triangulations, that can be done
// by the functions but it is not a requirement.
//
// Possible usage for this:
// * Computing the vector inside a cone of fixed norm in a
//   Polyhedral cone using Copositive programming.
// * Computing an integral over a polytope.
//
// The templatized functions must look like:
// * The f_lrs that takes a triangle in lrs (A mapping of the function
//     f used for Kernel_Simplices_cond.
// * A type T_merge that is updated.
// * A function f_compilation that compiles everything that has
//   been computed.
//   --For the integral of a polynomial, it is just an average over
//     a group.
//   --For the fixed norm vectors, it is just the compilation of the
//     database, so nothing particular.
// * The function is then

template<typename TintGroup>
struct PolyHeuristicDecomposition {
  TheHeuristic<TintGroup> Splitting;
  TheHeuristic<TintGroup> InitialFacetSet;
};



// The equivalence:
// * The matrix getting the equivalence.
// * The f corresponding to the facet.
// * The corresponding orbit, if i_orbit = miss_val then
//   the orbit is not already known.
template<typename T>
struct DecompositionEquiv {
  Face f;
  std::option<std::pair<i_orbit, MyMatrix<T>>> index_mat;
};

// A component of the decomposition
// * The list of vertices.
// * vf: The triangulation of the object
// * ListEquiv: The list of transformation mapping to the adjacent domains.
template<typename T, typename Tgroup>
struct ComponentDecomposition {
  MyMatrix<T> EXT;
  Tgroup GRPperm;
  vectface vf_trig;
  std::vector<DecompositionEquiv<T,Tint>> ListEquiv;
};

template<typename T>
MyMatrix<T> GetFacetIso(MyMatrix<T> const& EXT, Face const& f) {
  int nbCol = EXT.cols();
  MyVector<T> eVect = Isobarycenter(EXT);
  MyMatrix<T> EXTret(f.count() + 1, nbCol);
  int iRow = 0;
  for (auto & eRow : FaceToVector<int>(f)) {
    for (int iCol=0; iCol<nbCol; iCol++) {
      EXTret(iRow, iCol) = EXT(eRow, iCol);
    }
    iRow += 1;
  }
  for (int iCol=0; iCol<nbCol; iCol++) {
    EXTret(iRow, iCol) = eVect(iCol);
  }
  return EXTret;
}

vectface vf_add_one_included_vertex(vectface const& vf) {
  size_t n = vf.get_n();
  vectface vf_ret(n+1);
  Face fnew(n+1);
  fnew[n] = 1;
  for (auto & eFace : vf) {
    for (size_t i=0; i<n; i++) {
      fnew[i] = eFace[i];
    }
    vf_ret.push_back(fnew);
  }
  return vf_ret;
}

Face face_add_one_included_vertex(Face const& eFace) {
  size_t n = eFace.size();
  Face fnew(n + 1);
  fnew[n] = 1;
  for (size_t i=0; i<n; i++) {
    fnew[i] = eFace[i];
  }
  return fnew;
}

template<typename Tgroup>
Tgroup group_add_one_vertex(Tgroup const& GRP) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n_act = GRP.n_act();
  std::vector<Telt> ListGen;
  for (auto & eGen : GRP.GeneratorsOfGroup()) {
    std::vector<Tidx> eList(n_act + 1);
    for (Tidx i=0; i<n_act; i++) {
      eList[i] = eGen.at(i);
    }
    eList[n_act] = n_act;
    Telt ePerm(std::move(eList));
    ListGen.emplace_back(std::move(ePerm));
  }
  return Tgroup(ListGen, n_act + 1);
}




template<typename T, typename Tgroup>
std::vector<ComponentDecomposition<T,Tgroup>> get_full_decomposition(MyMatrix<T> const& EXT, Tgroup const& GRP,
                                                                     PolyHeuristicDecomposition<typename Tgroup::Tint> const& AllArr, std::ostream& os) {
  using TintGroup = typename Tgroup::Tint;
  std::map<std::string, TintGroup> TheMap = ComputeInitialMap(EXT, GRP, AllArr.dimEXT);
  std::vector<ComponentDecomposition<T>> l_comp_decomp;
  std::string ansSplit = HeuristicEvaluation(TheMap, AllArr.Splitting);
  if (ansSplit != "split") {
    std::pair<vectface, vectface> pair = GetTriangulationFacet(EXT);
    vectface ListOrbFacet = OrbitSplittingSet(pair.second, GRP);
    std::vector<DecompositionEquiv<T,Tint>> ListEquiv;
    for (auto & f : ListOrbFacet) {
      DecompositionEquiv<T> de{f, {}};
      ListEquiv.push_back(de);
    }
    ComponentDecomposition<T,Tgroup> cd{EXT, GRP, pair.first, ListEquiv};
    std::vector<ComponentDecomposition<T,Tgroup>> vect_cd{cd};
    return vect_cd;
  }
  //
  // Decomposition the space
  // We have to deal with:
  // * The complexity of the different 
  //
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  MyMatrix<T> EXT = ColumnReduction(EXT);
  MyVector<T> eVectRed = Isobarycenter(EXT);
  std::string ansSamp = HeuristicEvaluation(TheMap, AllArr.InitialFacetSet);
  vectface vf_samp = DirectComputationInitialFacetSet(EXT, ansSamp, os);
  Face f_first = vf_samp[0];
  // The ways of the managing of the indices.
  struct pairIdx {
    size_t i_facet;
    size_t i_orbit;
  };
  std::map<pairIdx, size_t> map_pair;
  std::vector<pairIdx> vect_pair;

  std::vector<ComponentDecomposition<T,Tgroup>> vect_cd;
  //
  std::vector<Tgroup> list_stab1;
  std::vector<Tgroup> list_stab2;
  std::vector<Tgroup> list_EXT2;
  vectface vf_facet(nbRow);
  std::vector<size_t> list_siz_n_comp;
  std::vector<size_t> list_shift_n_comp;
  auto f_insert_direct=[&](Face const& f) -> void {
    Tgroup eStab1 = GRP.Stabilizer_OnSets(f);
    Tgroup eStab2 = ReducedGroupActionAddPtFace(eStab1, f);
    MyMatrix<T> EXT1 = SelectRow(EXT, f);
    MyMatrix<T> EXT2 = ColumnReduction(EXT2);
    list_stab1.push_back(eStab1);
    list_stab2.push_back(eStab2);
    list_EXT2.push_back(EXT2);
    vf_facet.push_back(f);
  };
  size_t n_done = 0;
  while(true) {
    size_t len = list_stab2.len();
    for (size_t pos=n_done; pos<len; pos++) {
      Face facet = vf_facet[pos];
      Tgroup eStab2 = list_stab2[pos];
      MyMatrix<T> EXT2 = list_EXT2[pos];
      std::vector<ComponentDecomposition<T,Tgroup>> vec_cd = get_full_decomposition(EXT2, eStab2, AllArr, os);
      for (auto & ecd : vec_cd) {
        vectface vf_trig = vf_add_one_included_vertex(ecd.vf_trig);
        MyMatrix<T> EXTcomp = GetFacetIso(EXT, facet);
        Tgroup GRPperm = group_add_one_vertex(ecd.GRPperm);
        std::vector<DecompositionEquiv<T>> ListEquiv;
        for (auto & equiv : ecd.ListEquiv) {
          Face f = vf_add_one_included_vertex(equiv.f);
          auto iife_equiv=[&]() -> std::optional<std::pair<size_t, MyMatrix<T>>> {
            if (equiv) {
              std::pair<size_t, MyMatrix<T>> const& pair = *equiv;

            } else {
              return {};
            }
          };
          std::optional<size_t, MyMatrix<T>> new_equiv = iife_equiv();
          ListEquiv.push_back(new_equiv);
        }
        
      }
    }

    
  }



  
}



// clang-format off
#endif  // SRC_DUALDESC_POLY_EQUIVARIANTTRIANGULATION_H_
// clang-format on

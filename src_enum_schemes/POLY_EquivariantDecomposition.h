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
MyMatrix<T> get_transformation(MyMatrix<T> const& EXTfull, MyMatrix<T> const& EXT, MyVector<T> const& eVect, ComponentDecomposition<T, Tgroup> const& cd_adj, ComponentDecomposition<T, Tgroup> const& cd, Face const& f, MyMatrix<T> const& P) {
  //
  // The polytopes and their containers to find the indiced
  //
  int nbColFull = EXTfull.cols();
  ContainerMatrix<T> ContEXT(EXT);
  MyMatrix<T> cd_adjEXTimg = cd_adj.EXT * P;
  ContainerMatrix<T> Cont_cd_adjEXTimg(cd_adjEXTimg);
  MyMatrix<T> cd_EXT = SelectRow(cd.EXT, f);
  //
  // Building the corresponding matrices
  //
  int siz = f.count();
  MyMatrix<T> EXTfull1(siz + 1, nbColFull);
  MyMatrix<T> EXTfull2(siz + 1, nbColFull);
  for (int i=0; i<siz; i++) {
    MyVector<T> V1 = GetMatrixRow(cd_EXT, i);
    std::optional<size_t> opt1 = ContEXT.GetIdx_v(V1);
    size_t pos1 = unfold_opt(opt1, "cd_EXT should be contained in EXT");
    MyVector<T> V2 = GetMatrixRow(EXTfull, pos1);
    AssignMatrixRow(EXTfull1, i, V2);
    //
    std::optional<size_t> opt2 = Cont_cd_adjEXTimg.GetIdx_v(V1);
    size_t pos2 = unfold_opt(opt2, "cd_EXT should be contained in cd_adjEXTimg");
    MyVector<T> V3 = GetMatrixRow(cd_adj.EXT, pos2);
    std::optional<size_t> opt3 = ContEXT.GetIdx_v(V3);
    size_t pos3 = unfold_opt(opt3, "cd_adj.eXT should be contained in EXT");
    MyVector<T> V4 = GetMatrixRow(EXTfull, pos3);
    AssignMatrixRow(EXTfull2, i, V4);
  }
  AssignMatrixRow(EXTfull1, siz, eVect);
  AssignMatrixRow(EXTfull2, siz, eVect);
  //
  // Finding the relevant transformation
  //
  return FindTransformationId(EXTfull2, EXTfull1);
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
  // * The complexity of the different indices from the subcomponents (hence the pairIdx)
  // * We decompose the faces into smaller ones. The problem is that we need to have the
  //   faces matching. And if we have decomposed one side but the opposite is not decomposed
  //   then we have a problem. Therefore, the decomposition has to use only intrinsic
  //   information
  //   - The number of rays.
  //   - The dimension.
  //   - The size of the automorphism group.
  // * The above problem being addressed, we need anyway to compute the mapping from a
  //   component of a facet to a 
  //
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  MyMatrix<T> EXT = ColumnReduction(EXT);
  MyVector<T> eVect = Isobarycenter(EXT);
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
  auto f_insert_pair_idx=[&](size_t const& i_facet, size_t const& i_orbit) -> void {
    pairIdx pi{i_facet, i_orbit};
    size_t len = vect_pair.size();
    map_pair[pi] = len + 1;
    vect_pair.push_back(pi);
  };
  auto get_idx_from_pair=[&](size_t const& i_facet, size_t const& i_orbit) -> size_t {
    pairIdx pi{i_facet, i_orbit};
    size_t pos = map_pair[pi];
    if (pos == 0) {
      std::cerr << "The entry is actually is missing\n";
      throw TerminalException{1};
    }
    return pos - 1;
  };
  //
  // The returned data
  //
  std::vector<ComponentDecomposition<T,Tgroup>> vec_cd;
  //
  // The array of data corresponding to the facet
  //
  std::vector<Tgroup> list_stab1;
  std::vector<Tgroup> list_stab2;
  std::vector<MyMatrix<T>> list_EXT1;
  std::vector<MyMatrix<T>> list_EXT2;
  vectface vf_facet(nbRow);
  std::vector<size_t> list_siz_n_orbit;
  std::vector<size_t> list_shift_n_orbit;
  auto f_insert_direct=[&](Face const& f) -> void {
    Tgroup eStab1 = GRP.Stabilizer_OnSets(f);
    Tgroup eStab2 = ReducedGroupActionAddPtFace(eStab1, f);
    MyMatrix<T> EXT1 = SelectRow(EXT, f);
    MyMatrix<T> EXT2 = ColumnReduction(EXT2);
    list_stab1.push_back(eStab1);
    list_stab2.push_back(eStab2);
    list_EXT1.push_back(EXT1);
    list_EXT2.push_back(EXT2);
    vf_facet.push_back(f);
  };
  f_insert_direct(f_first);
  //
  // The main iteration.
  //
  struct FacetEntry {
    Face facet;
    std::vector<Face> l_face;
  };
  std::vector<std::vector<FacetEntry>> ll_facetentry;

  size_t n_done = 0;
  while(true) {
    size_t len = list_stab2.len();
    for (size_t pos=n_done; pos<len; pos++) {
      Face facet = vf_facet[pos];
      Tgroup eStab2 = list_stab2[pos];
      MyMatrix<T> EXT1 = list_EXT1[pos];
      MyMatrix<T> EXT2 = list_EXT2[pos];
      ContainerMatrix<T> ContEXT2(EXT2);
      std::vector<ComponentDecomposition<T,Tgroup>> local_vec_cd = get_full_decomposition(EXT2, eStab2, AllArr, os);
      size_t n_orbit = local_vec_cd.size();
      int nbVect = EXT1.rows();
      list_siz_n_orbit.push_back(n_comp);
      list_shift_n_orbit.push_back(vec_cd.size());
      for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++) {
        f_insert_pair_idx(pos, i_orbit);
      }
      vectface vf_facet(nbRow);
      for (auto & ecd : local_vec_cd) {
        vectface vf_trig = vf_add_one_included_vertex(ecd.vf_trig);
        MyMatrix<T> EXTcomp = GetFacetIso(EXT, facet);
        Tgroup GRPperm = group_add_one_vertex(ecd.GRPperm);
        std::vector<DecompositionEquiv<T>> ListEquiv;
        for (auto & equiv : ecd.ListEquiv) {
          Face f = vf_add_one_included_vertex(equiv.f);
          auto iife_equiv=[&]() -> std::optional<std::pair<size_t, MyMatrix<T>>> {
            if (equiv) {
              std::pair<size_t, MyMatrix<T>> const& pair = *equiv;
              size_t const& j_orbit = pair.first;
              MyMatrix<T> const& P = pair.second;
              ComponentDecomposition<T,Tgroup> const& cd_adj = local_vec_cd[j_orbit];
              MyMatrix<T> P_final = get_transformation(EXT1, EXT2, eVect, cd_adj, ecd, equiv.f, P);
              size_t j_orbit_final = get_idx_from_pair(pos, j_orbit);
              std::pair<size_t, MyMatrix<T>> pair_final{j_orbit_final, P_final};
              return pair_final;
            } else {
              Face new_facet(nbRow);
              for (auto & val : FaceToVector<size_t>(equiv.f)) {
                MyVector<T> V1 = GetMatrixRow(ecd.EXT, val);
                std::optional<size_t> opt1 = ContEXT2.GetIdx_v(V1);
                size_t pos1 = unfold_opt(opt1, "ecd.EXT should be contained in EXT2");
                new_facet[pos1] = 1;
              }
              vf_facet.push_back(new_facet);
              return {};
            }
          };
          std::optional<size_t, MyMatrix<T>> new_equiv = iife_equiv();
          DecompositionEquiv<T> de{f, new_equiv};
          ListEquiv.push_back(de);
        }
        Face f(facet.count() + 1);
        for (size_t u=0; u<facet.count(); u++) {
          f[u] = 1;
        }
        DecompositionEquiv<T> de{f, {}};
        ListEquiv.push_back(de);
        //
        ComponentDecomposition<T,Tgroup> fcd{EXTcomp, GRPperm, vf_trig, ListEquiv};
        vec_cd.push_back(fcd);
      }
      //
      for (auto & eFace
    }
    n_done = len;
    if (list_stab2.size() == len) {
      break;
    }
  }



  
}



// clang-format off
#endif  // SRC_DUALDESC_POLY_EQUIVARIANTTRIANGULATION_H_
// clang-format on

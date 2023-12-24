// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ISODELAUNAYDOMAINS_H_
#define SRC_LATT_ISODELAUNAYDOMAINS_H_

// clang-format off
#include "POLY_LinearProgramming.h"
#include "ShortestUniversal.h"
#include "Temp_Positivity.h"
#include "LatticeDelaunay.h"
#include <string>
#include <vector>
// clang-format on

/*
  Code for the L-type domains.

  Two main use case:
  ---Lattice case: Then the Tvert is actually a Tint and can be mpz_class, int32_t, etc.
  ---Periodic structure case: Then the coordinates are no longer integral.
    Also the equivalence are no longer integral. Sure the matrix transformation is
    integral, but the translation vector is not necessarily so.

  We use the definitions from LatticeDelaunay.h
  They are for the lattice case but can be generalized for the periodic case.
 */


template<typename T>
struct VoronoiInequalityPreComput {
  //  MyMatrix<Tvert> VertBasis;
  MyMatrix<T> VertBasis_T;
  MyMatrix<T> VertBasisInv_T;
  std::vector<MyVector<T>> VertBasisRed_T;
};

template<typename T>
VoronoiInequalityPreComput<T> BuildVoronoiIneqPreCompute(MyMatrix<Tvert> const& EXT) {
  int n = EXT.cols() - 1;
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T,Tvert>(EXT);
  MyMatrix<T> VertBasis_T = RowReduction(EXT_T);
  MyMatrix<T> VertBasisInv_T = TransposedMat(Inverse(VertBasis_T));
  std::vector<MyVector<T>> VertBasisRed_T;
  for (int i=0; i<=n; i++) {
    MyVector<T> V(n);
    for (int u=0; u<n; u++) {
      V(u) = VertBasis_T(i,u+1);
    }
    VertBasisRed_T.push_back(V);
  }
  return {std::move(VertBasis_T), std::move(VertBasisInv_T), std::move(VertBasisRed_T)};
}


template<typename T, typename Tvert>
MyVector<T> VoronoiLinearInequality(VoronoiInequalityPreComput<T> const& vipc, MyVector<Tvert> const& TheVert, std::vector<std::vector<T>> const& ListGram) {
  int n = ListGram[0].rows();
  int dimSpace = ListGram.size();
  MyVector<T> TheVert_T = UniversalVectorConversion<T,Tvert>(TheVert);
  MyVector<T> B = vipc.VertBasisInv_T * TheVert_T.transpose();
  MyVector<T> Ineq(dimSpace);
  MyVector<T> TheVertRed(n);
  for (int u=0; u<n; u++) {
    TheVertRed(u) = TheVert_T(u+1);
  }
  int iGram = 0;
  for (auto & eLineMat : ListGram) {
    T val = EvaluateLineVector(eLineMat, TheVertRed);
    for (int k=0; k<=n; k++) {
      val += B(k) * EvaluateLineVector(eLineMat, vipc.VertBasisRed_T[k]);
    }
    Ineq(iGram) = val;
    iGram++;
  }
  return Ineq;
}

struct AdjInfo {
  int iOrb;
  int i_adj;
};


/*
  Compute the defining inequalities of an iso-Delaunay domain
 */
template<typename T, typename Tvert, typename Tgroup>
std::unordered_map<MyVector<T>,std::vector<AdjInfo>> ComputeDefiningIneqIsoDelaunayDomain(DelaunayTesselation<Tvert, Tgroup> const& DT, std::vector<std::vector<T>> const& ListGram) {
  std::unordered_map<MyVector<T>,std::vector<AdjInfo>> map;
  int n_del = DT.l_dels.size();
  for (int i_del=0; i_del<n_del; i_del++) {
    int n_adj = DT.l_dels.l_adj.size();
    VoronoiInequalityPreComput<T> vipc = BuildVoronoiIneqPreCompute(DT.l_dels[i_del].obj);
    ContainerMatrix<Tvert> cont(DT.l_dels[i_del].obj);
    auto get_ineq=[&](int const& i_adj) -> MyVector<T> {
      Delaunay_AdjO<Tvert> adj = DT.l_dels[i_del].l_adj[i_adj];
      int j_del = adj.iOrb;
      MyMatrix<Tvert> EXTadj = DT.l_dels[j_del].obj * adj.P;
      int len = EXTadj.rows();
      for (int u=0; u<len; u++) {
        MyVector<Tvert> TheVert = GetMatrixRow(EXTadj, u);
        std::optional<size_t> opt = cont.GetIdx_v(TheVert);
        if (!opt) {
          return VoronoiLinearInequality(vipc, TheVert, ListGram);
        }
      }
      std::cerr << "Failed to find a matching entry\n";
      throw TerminalException{1};
    };
    for (int i_adj=0; i_adj<n_adj; i_adj++) {
      MyVector<T> V = get_ineq(i_adj);
      MyVector<T> V_red = ScalarCanonicalizationVector(V);
      AdjInfo eAdj{i_del, i_adj};
      map[V_red].push_back(eAdj);
    }
  }
  return map;
}







// clang-format off
#endif  // SRC_LATT_ISODELAUNAYDOMAINS_H_
// clang-format on

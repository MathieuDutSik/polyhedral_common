// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_RIGIDITY_H_
#define SRC_DELAUNAY_RIGIDITY_H_

// clang-format off
#include "IsoDelaunayDomains.h"
#include "LatticeDelaunay.h"
#include "LatticeStabEquiCan.h"
#include "Tspace_Canonical.h"
#include "MAT_MatrixNullspace.h"
#include <vector>
// clang-format on

/*
  C++ port of MyPolyhedral/lib/LatticeDelaunays.g::Lspace and its caller
  GetRigidityDegree.

  GAP reference (annotated):
    Lspace := function(DataLattice, DelaunayDatabase)
      n := DataLattice.n;
      TheCan := CanonicalSymmetricBasis(n);     # full n(n+1)/2-dim Sym^n
      TheBasisMatrix := TheCan.TheBasisMatrix;
      TheVectorSpace := TheCan.TheVectorSpace;
      NewMatrixGRP := GetSym2Representation(n, DataLattice.PointGRP);
      for iOrb in [1..DelaunayDatabase.FuncDelaunayGetNumber()] do
        EXT := DelaunayDatabase.FuncDelaunayGetEXT(iOrb);
        ListEqua := ListEqualitiesDelaunay(EXT, TheBasisMatrix);
        NSPvector := NullspaceMat(TransposedMat(ListEqua));      # right null
        TheVectorSpaceProv := ...(NSPvector expressed in ambient);
        TheVectorSpace := IntersectionOrbitSubspace(TheVectorSpaceProv,
                                                    NewMatrixGRP);
        TheBasisMatrix := List(TheVectorSpace.ListVect,
                               x -> VectorToSymmetricMatrix(x, n));
      od;
      return TheBasisMatrix;
    end;
    GetRigidityDegree() := Length(Lspace(...))

  Key subtlety: IntersectionOrbitSubspace(V, G) returns the LARGEST G-stable
  subspace contained in V (not the G-fixed subspace). For Z^n the diagonal
  forms are permuted but stay diagonal under the hyperoctahedral group, so
  the diagonal subspace is G-stable; the G-fixed subspace would only be the
  ray R*I_n, which gives the wrong answer.
 */

template <typename T> MyVector<T> SymMatrixToVecRig(MyMatrix<T> const &M) {
  int n = M.rows();
  int dim = (n * (n + 1)) / 2;
  MyVector<T> v(dim);
  int pos = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      v(pos) = M(i, j);
      pos++;
    }
  }
  return v;
}

/*
  Intersection of two subspaces of Sym^n given by bases A and B (each entry
  an n x n symmetric matrix). Returns a linearly independent basis of A ∩ B.

  Each n x n symmetric matrix is encoded as an n(n+1)/2 vector of its
  upper-triangular entries. With M of shape (|A| + |B|, n(n+1)/2) holding A
  on top and B on bottom, the left null space of M parametrises pairs (a, b)
  satisfying a · A + b · B = 0 (i.e. -a · A = b · B lies in A ∩ B).
 */
template <typename T>
std::vector<MyMatrix<T>>
IntersectionSubspaceMat(std::vector<MyMatrix<T>> const &A,
                        std::vector<MyMatrix<T>> const &B, int n) {
  if (A.empty() || B.empty()) {
    return {};
  }
  int dim_amb = (n * (n + 1)) / 2;
  int kA = A.size();
  int kB = B.size();
  MyMatrix<T> M(kA + kB, dim_amb);
  for (int i = 0; i < kA; i++) {
    MyVector<T> v = SymMatrixToVecRig(A[i]);
    AssignMatrixRow(M, i, v);
  }
  for (int i = 0; i < kB; i++) {
    MyVector<T> v = SymMatrixToVecRig(B[i]);
    AssignMatrixRow(M, kA + i, v);
  }
  // Left null space: rows w with w · M = 0. Split each w into (a, b);
  // then u = -sum_j a_j A[j] (= sum_j b_j B[j]) lies in A ∩ B.
  MyMatrix<T> NSP = NullspaceMat(M);
  int dimInt = NSP.rows();
  std::vector<MyMatrix<T>> result;
  result.reserve(dimInt);
  for (int i = 0; i < dimInt; i++) {
    MyMatrix<T> uMat = ZeroMatrix<T>(n, n);
    for (int j = 0; j < kA; j++) {
      uMat -= NSP(i, j) * A[j];
    }
    result.push_back(uMat);
  }
  return result;
}

/*
  Image of a subspace V of Sym^n under a single generator g acting as
  M -> g * M * g^T. Returns the same number of basis matrices.
 */
template <typename T>
std::vector<MyMatrix<T>> ApplyGenSym2(std::vector<MyMatrix<T>> const &V,
                                      MyMatrix<T> const &g) {
  std::vector<MyMatrix<T>> result;
  result.reserve(V.size());
  MyMatrix<T> gT = g.transpose();
  for (auto &M : V) {
    result.push_back(g * M * gT);
  }
  return result;
}

/*
  Port of MyPolyhedral/lib/BasicGeomNbr.g::IntersectionOrbitSubspace:
  iteratively replace V by ∩_{g in gens} (V ∩ g·V) until stable. The fixed
  point is the largest G-stable subspace inside the initial V.
 */
template <typename T>
std::vector<MyMatrix<T>>
LargestGStableSubspace(std::vector<MyMatrix<T>> const &V_init,
                       std::vector<MyMatrix<T>> const &ListGen, int n) {
  std::vector<MyMatrix<T>> V = V_init;
  while (true) {
    std::vector<MyMatrix<T>> NewV = V;
    for (auto &g : ListGen) {
      std::vector<MyMatrix<T>> gV = ApplyGenSym2(V, g);
      NewV = IntersectionSubspaceMat(NewV, gV, n);
      if (NewV.empty()) {
        return NewV;
      }
    }
    if (NewV.size() == V.size()) {
      return V;
    }
    V = NewV;
  }
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<T>>
ComputeLspace(MyMatrix<T> const &GramMat,
              DelaunayTesselation<T, Tgroup> const &DT, std::ostream &os) {
  int n = GramMat.rows();
  // Initial basis: full canonical Sym^n (the n*(n+1)/2 elementary
  // symmetric matrices), matching TheCan.TheBasisMatrix in the GAP code.
  std::vector<MyMatrix<T>> Lbasis =
      TSPACE_canonical_get_list_matrices<T>(n);
  // Lattice automorphism generators acting on R^n; their Sym^2
  // representation is the M -> g M g^T action used by ApplyGenSym2.
  std::vector<MyMatrix<Tint>> ListGen_int =
      ArithmeticAutomorphismGroup<T, Tint, Tgroup>(GramMat, os);
  std::vector<MyMatrix<T>> ListGen_T;
  for (auto &eGen : ListGen_int) {
    ListGen_T.push_back(UniversalMatrixConversion<T, Tint>(eGen));
  }
  for (auto &eDel : DT.l_dels) {
    if (Lbasis.empty()) {
      break;
    }
    // Equalities of this Delaunay orbit expressed in current Lbasis coords.
    std::vector<std::vector<T>> ListLineMat;
    for (auto &eMat : Lbasis) {
      ListLineMat.push_back(GetLineVector(eMat));
    }
    VoronoiInequalityPreComput<T> vipc =
        BuildVoronoiIneqPreCompute<T>(eDel.EXT, os);
    std::vector<MyVector<T>> ListEqua;
    int n_vert = eDel.EXT.rows();
    for (int u = 0; u < n_vert; u++) {
      MyVector<T> TheVert = GetMatrixRow(eDel.EXT, u);
      MyVector<T> eEq =
          VoronoiLinearInequality(vipc, TheVert, ListLineMat, os);
      ListEqua.push_back(eEq);
    }
    MyMatrix<T> MatEqua = MatrixFromVectorFamily(ListEqua);
    // Right null space (GAP's NullspaceMat(TransposedMat(...))).
    MyMatrix<T> NSP = NullspaceTrMat(MatEqua);
    int new_dim = NSP.rows();
    std::vector<MyMatrix<T>> ProvLbasis;
    ProvLbasis.reserve(new_dim);
    for (int i = 0; i < new_dim; i++) {
      MyMatrix<T> NewMat = ZeroMatrix<T>(n, n);
      int cur_dim = Lbasis.size();
      for (int j = 0; j < cur_dim; j++) {
        NewMat += NSP(i, j) * Lbasis[j];
      }
      ProvLbasis.push_back(NewMat);
    }
    // Apply G-stability AFTER each Delaunay (mirrors the GAP loop body).
    Lbasis = LargestGStableSubspace(ProvLbasis, ListGen_T, n);
  }
  return Lbasis;
}

template <typename T, typename Tint, typename Tgroup>
int ComputeRigidityDegreeLattice(MyMatrix<T> const &GramMat,
                                 DelaunayTesselation<T, Tgroup> const &DT,
                                 std::ostream &os) {
  std::vector<MyMatrix<T>> Lbasis =
      ComputeLspace<T, Tint, Tgroup>(GramMat, DT, os);
  return static_cast<int>(Lbasis.size());
}

#endif  // SRC_DELAUNAY_RIGIDITY_H_

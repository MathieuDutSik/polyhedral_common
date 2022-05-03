#ifndef SRC_POLY_POLY_CDD_GRAPH_H_
#define SRC_POLY_POLY_CDD_GRAPH_H_

#include "GRAPH_BitsetType.h"
#include "POLY_cddlib.h"

namespace cdd {
template <typename T>
GraphBitset GetRidgeGraph(dd_polyhedradata<T> const *poly, int const &nbFac) {
  dd_raydata<T> *RayPtr1;
  dd_raydata<T> *RayPtr2;
  long pos1, pos2;
  bool adj;
  GraphBitset RidgeGraph(nbFac);
  poly->child->LastRay->Next = nullptr;
  for (RayPtr1 = poly->child->FirstRay, pos1 = 1; RayPtr1 != nullptr;
       RayPtr1 = RayPtr1->Next, pos1++)
    for (RayPtr2 = poly->child->FirstRay, pos2 = 1; RayPtr2 != nullptr;
         RayPtr2 = RayPtr2->Next, pos2++)
      if (RayPtr1 != RayPtr2) {
        dd_CheckAdjacency(poly->child, &RayPtr1, &RayPtr2, &adj);
        if (adj)
          RidgeGraph.AddAdjacent(pos1 - 1, pos2 - 1);
      }
  return RidgeGraph;
}

template <typename T> GraphBitset GetSkeletonGraph(dd_polyhedradata<T> *poly) {
  dd_rowrange i, j;
  if (poly->AincGenerated == false)
    dd_ComputeAinc(poly);
  int nbVert = poly->m1;
  GraphBitset SkelGraph(nbVert);
  set_type common;
  long lastn = 0;
  for (i = 1; i <= nbVert; i++)
    for (j = 1; j <= nbVert; j++)
      if (i != j && dd_InputAdjacentQ<T>(common, lastn, poly, i, j))
        SkelGraph.AddAdjacent(i - 1, j - 1);
  return SkelGraph;
}

template <typename T> struct DDA {
  MyMatrix<T> EXT;
  GraphBitset SkelGraph;
  GraphBitset RidgeGraph;
};

template <typename T>
DDA<T> DualDescriptionAdjacencies(MyMatrix<T> const &TheEXT) {
  dd_polyhedradata<T> *poly;
  dd_matrixdata<T> *M;
  dd_ErrorType err;
  int nbCol = TheEXT.cols();
  M = MyMatrix_PolyFile2Matrix(TheEXT);
  poly = dd_DDMatrix2Poly(M, &err);

  MyMatrix<T> TheFAC = FAC_from_poly(poly, nbCol);
  int nbFac = TheFAC.rows();
  GraphBitset eRidgeGraph = GetRidgeGraph(poly, nbFac);
  GraphBitset eSkelGraph = GetSkeletonGraph(poly);
  dd_FreePolyhedra(poly);
  dd_FreeMatrix(M);
  return DDA<T>{TheFAC, eSkelGraph, eRidgeGraph};
}
} // namespace cdd

// clang-format off
#endif  // SRC_POLY_POLY_CDD_GRAPH_H_
// clang-format on

// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#define GMPRATIONAL
#ifdef USE_CDDLIB

#include "Boost_bitset_kernel.h"
#include "cdd.h"
#include "gmpxx.h"
#include "setoper.h"
#include <boost/multiprecision/gmp.hpp>
#include <vector>

#include <boost/dynamic_bitset.hpp>
typedef boost::dynamic_bitset<> Face;

#include "MatrixTypes.h"

namespace cbased_cdd {

template <typename T, typename Fset>
vectface DualDescription_incd_T(MyMatrix<T> const &TheEXT, Fset fset) {
  dd_ErrorType err;
  size_t nbCol = TheEXT.cols();
  size_t nbRow = TheEXT.rows();
  vectface ListIncd(nbRow);
  dd_set_global_constants();
  dd_MatrixPtr M = nullptr;
  dd_rowrange m_input, i;
  dd_colrange d_input, j;
  dd_RepresentationType rep;
  m_input = TheEXT.rows();
  d_input = TheEXT.cols();
  rep = dd_Generator; /* using dd_Inequality led to horrible bugs */
  M = dd_CreateMatrix(m_input, d_input);
  M->representation = rep;

  for (i = 0; i < m_input; i++) {
    for (j = 0; j < d_input; j++) {
      T val_T = TheEXT(i, j);
      fset(M->matrix[i][j], val_T);
    }
  }
  M->representation = dd_Generator;
  //
  dd_polyhedradata *poly = dd_DDMatrix2Poly(M, &err);
  //
  dd_raydata *RayPtr = poly->child->FirstRay;
  std::vector<T> V(nbCol);
  T eScal;
  auto isincd = [&](size_t iRow) -> bool {
    eScal = 0;
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      eScal += V[iCol] * TheEXT(iRow, iCol);
    return eScal == 0;
  };
  while (RayPtr != nullptr) {
    if (RayPtr->feasible) {
      for (size_t iCol = 0; iCol < nbCol; iCol++) {
        T eVal(RayPtr->Ray[iCol]);
        V[iCol] = eVal;
      }
      ListIncd.InsertFaceRef(isincd);
    }
    RayPtr = RayPtr->Next;
  }
  dd_FreePolyhedra(poly);
  dd_FreeMatrix(M);
  dd_free_global_constants();
  return ListIncd;
}

vectface DualDescription_incd_mpq_class(MyMatrix<mpq_class> const &TheEXT) {
  auto fset = [](mpq_t &ptr, mpq_class const &val) -> void {
    mpq_set(ptr, val.get_mpq_t());
  };
  return DualDescription_incd_T(TheEXT, fset);
}

vectface DualDescription_incd_boost_mpq_rational(
    MyMatrix<boost::multiprecision::mpq_rational> const &TheEXT) {
  auto fset = [](mpq_t &ptr,
                 boost::multiprecision::mpq_rational const &val) -> void {
    mpq_set(ptr, val.backend().data());
  };
  return DualDescription_incd_T(TheEXT, fset);
}

// clang-format off
}  // namespace cbased_cdd
// clang-format on
#endif

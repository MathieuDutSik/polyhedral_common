#ifndef SRC_POLY_POLY_C_CDDLIB_H_
#define SRC_POLY_POLY_C_CDDLIB_H_

#ifdef USE_CDDLIB

// clang-format off
#include "setoper.h"
#include "POLY_c_cddlib_mpq.h"
#include "cdd.h"
#include <type_traits>
#include <vector>
// clang-format on

namespace cbased_cdd {

//
// Removal of equation using double computation
//

template <typename T>
dd_MatrixPtr MyMatrix_PolyFile2Matrix_double(MyMatrix<T> const &TheEXT) {
  dd_MatrixPtr M = nullptr;
  dd_rowrange m_input, i;
  dd_colrange d_input, j;
  dd_RepresentationType rep;
  m_input = TheEXT.rows();
  d_input = TheEXT.cols();
  rep = dd_Generator; /* using dd_Inequality led to horrible bugs */
  M = dd_CreateMatrix(m_input, d_input);
  M->representation = rep;

  for (i = 0; i < m_input; i++)
    for (j = 0; j < d_input; j++) {
      T val_T = TheEXT(i, j);
      double val_d = UniversalScalarConversion<double, T>(val_T);
      dd_set(M->matrix[i][j], &val_d);
    }
  return M;
}

template <typename T>
std::vector<int> RedundancyReductionClarkson(MyMatrix<T> const &TheEXT) {
  dd_ErrorType err;
  int nbRow = TheEXT.rows();
  dd_set_global_constants();
  dd_MatrixPtr M = MyMatrix_PolyFile2Matrix_double(TheEXT);
  M->representation = dd_Inequality; // The choice between dd_Inequality and
                                     // dd_Generator is not clear.
  dd_rowset redset = dd_RedundantRowsViaShooting(M, &err);
  std::vector<int> ListIdx;
  for (int i_row = 0; i_row < nbRow; i_row++) {
    bool isin = set_member(i_row + 1, redset);
    if (!isin)
      ListIdx.push_back(i_row);
  }
  dd_FreeMatrix(M);
  set_free(redset);
  dd_free_global_constants();
  return ListIdx;
}

//
// Removal of equation using mpq_class
//

template <typename T>
dd_MatrixPtr MyMatrix_PolyFile2Matrix_T(MyMatrix<T> const &TheEXT) {
  dd_MatrixPtr M = nullptr;
  dd_rowrange m_input, i;
  dd_colrange d_input, j;
  dd_RepresentationType rep;
  m_input = TheEXT.rows();
  d_input = TheEXT.cols();
  rep = dd_Generator; /* using dd_Inequality led to horrible bugs */
  M = dd_CreateMatrix(m_input, d_input);
  M->representation = rep;

  for (i = 0; i < m_input; i++)
    for (j = 0; j < d_input; j++) {
      T val_T = TheEXT(i, j);
      mpq_set(M->matrix[i][j], val_T.get_mpq_t());
    }
  return M;
}

template <typename T>
inline typename std::enable_if<is_mpq_class<T>::value, vectface>::type
DualDescription_incd(MyMatrix<T> const &TheEXT) {
  return DualDescription_incd_mpq_class(TheEXT);
}

template <typename T>
inline typename std::enable_if<is_boost_mpq_rational<T>::value, vectface>::type
DualDescription_incd(MyMatrix<T> const &TheEXT) {
  return DualDescription_incd_boost_mpq_rational(TheEXT);
}

template <typename T>
inline typename std::enable_if<(not is_mpq_class<T>::value) &&
                                   (not is_boost_mpq_rational<T>::value),
                               vectface>::type
DualDescription_incd(MyMatrix<T> const &TheEXT) {
  MyMatrix<mpq_class> EXT_mpq = UniversalMatrixConversion<mpq_class, T>(TheEXT);
  return DualDescription_incd_mpq_class(EXT_mpq);
}

// clang-format off
}  // namespace cbased_cdd
#endif  // USE_CDDLIB
#endif  // SRC_POLY_POLY_C_CDDLIB_H_
// clang-format on

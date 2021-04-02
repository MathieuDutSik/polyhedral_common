#define GMPRATIONAL

#include "gmpxx.h"
#include <type_traits>
#include "setoper.h"
#include "cdd.h"
#include "Boost_bitset.h"
#include "MAT_Matrix.h"
#include <vector>



namespace cbased_cdd {

//
// Removal of equation using mpq_class
//


dd_MatrixPtr MyMatrix_PolyFile2Matrix_mpq(MyMatrix<mpq_class> const& TheEXT)
{
  using T = mpq_class;
  dd_MatrixPtr M=nullptr;
  dd_rowrange m_input, i;
  dd_colrange d_input, j;
  dd_RepresentationType rep;
  m_input=TheEXT.rows();
  d_input=TheEXT.cols();
  rep=dd_Generator; /* using dd_Inequality led to horrible bugs */
  M=dd_CreateMatrix(m_input, d_input);
  M->representation=rep;

  for (i = 0; i < m_input; i++)
    for (j = 0; j < d_input; j++) {
      T val_T = TheEXT(i, j);
      mpq_set(M->matrix[i][j], val_T.get_mpq_t());
    }
  return M;
}


std::vector<Face> DualDescription_incd_mpq(MyMatrix<mpq_class> const& TheEXT)
{
  using T = mpq_class;
  dd_ErrorType err;
  size_t nbCol=TheEXT.cols();
  size_t nbRow=TheEXT.rows();
  dd_set_global_constants();
  dd_MatrixPtr M=MyMatrix_PolyFile2Matrix_mpq(TheEXT);
  M->representation = dd_Inequality; // The choice between dd_Inequality and dd_Generator is not clear.
  //
  dd_polyhedradata *poly;
  poly=dd_DDMatrix2Poly(M, &err);
  //
  std::vector<Face> ListIncd;
  dd_raydata* RayPtr = poly->child->FirstRay;
  while (RayPtr != nullptr) {
    if (RayPtr->feasible) {
      Face V(nbRow);
      for (size_t iRow=0; iRow<nbRow; iRow++) {
        T eScal=0;
        for (size_t iCol=0; iCol<nbCol; iCol++) {
          mpq_class eVal(RayPtr->Ray[iCol]);
          eScal += eVal * TheEXT(iRow,iCol);
        }
        if (eScal == 0)
          V[iRow]=1;
      }
      ListIncd.push_back(V);
    }
    RayPtr = RayPtr->Next;
  }
  dd_FreePolyhedra(poly);
  dd_FreeMatrix(M);
  dd_free_global_constants();
  return ListIncd;
}

}

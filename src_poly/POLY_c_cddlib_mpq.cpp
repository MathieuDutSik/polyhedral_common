#define GMPRATIONAL

#include "gmpxx.h"
#include "setoper.h"
#include "cdd.h"
#include <vector>
#include "Boost_bitset_kernel.h"

#include <boost/dynamic_bitset.hpp>
typedef boost::dynamic_bitset<> Face;



#include <Eigen/Dense>
template <typename T>
using MyMatrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;




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



vectface DualDescription_incd_mpq(MyMatrix<mpq_class> const& TheEXT)
{
  dd_ErrorType err;
  size_t nbCol=TheEXT.cols();
  size_t nbRow=TheEXT.rows();
  vectface ListIncd(nbRow);
  dd_set_global_constants();
  dd_MatrixPtr M=MyMatrix_PolyFile2Matrix_mpq(TheEXT);
  M->representation = dd_Generator;
  //
  dd_polyhedradata* poly = dd_DDMatrix2Poly(M, &err);
  //
  dd_raydata* RayPtr = poly->child->FirstRay;
  std::vector<mpq_class> V(nbCol);
  mpq_class eScal;
  auto process=[&](std::vector<mpq_class> const& V) {
    auto isincd=[&](size_t iRow) -> bool {
      eScal=0;
      for (size_t iCol=0; iCol<nbCol; iCol++)
        eScal += V[iCol] * TheEXT(iRow,iCol);
      return eScal == 0;
    };
    ListIncd.InsertFace(isincd);
  };
  while (RayPtr != nullptr) {
    if (RayPtr->feasible) {
      for (size_t iCol=0; iCol<nbCol; iCol++) {
        mpq_class eVal(RayPtr->Ray[iCol]);
        V[iCol] = eVal;
      }
      process(V);
    }
    RayPtr = RayPtr->Next;
  }
  dd_FreePolyhedra(poly);
  dd_FreeMatrix(M);
  dd_free_global_constants();
  return ListIncd;
}

}

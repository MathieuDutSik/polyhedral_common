#define GMPRATIONAL

#include "gmpxx.h"
#include "setoper.h"
#include "cdd.h"
#include <vector>
#include "Boost_bitset_kernel.h"

#include <boost/dynamic_bitset.hpp>
typedef boost::dynamic_bitset<> Face;


#include "MatrixTypes.h"


namespace cbased_cdd {

//
// Removal of equation using mpq_class
//



template<typename T, typename F>
vectface DualDescription_incd_T(MyMatrix<T> const& TheEXT, F f)
{
  dd_ErrorType err;
  size_t nbCol=TheEXT.cols();
  size_t nbRow=TheEXT.rows();
  vectface ListIncd(nbRow);
  dd_set_global_constants();
  dd_MatrixPtr M=nullptr;
  dd_rowrange m_input, i;
  dd_colrange d_input, j;
  dd_RepresentationType rep;
  m_input=TheEXT.rows();
  d_input=TheEXT.cols();
  rep=dd_Generator; /* using dd_Inequality led to horrible bugs */
  M=dd_CreateMatrix(m_input, d_input);
  M->representation=rep;

  for (i = 0; i < m_input; i++) {
    for (j = 0; j < d_input; j++) {
      T val_T = TheEXT(i, j);
      mpq_set(M->matrix[i][j], f(val_T));
    }
  }
  M->representation = dd_Generator;
  //
  dd_polyhedradata* poly = dd_DDMatrix2Poly(M, &err);
  //
  dd_raydata* RayPtr = poly->child->FirstRay;
  std::vector<T> V(nbCol);
  T eScal;
  auto isincd=[&](size_t iRow) -> bool {
    eScal=0;
    for (size_t iCol=0; iCol<nbCol; iCol++)
      eScal += V[iCol] * TheEXT(iRow,iCol);
    return eScal == 0;
  };
  while (RayPtr != nullptr) {
    if (RayPtr->feasible) {
      for (size_t iCol=0; iCol<nbCol; iCol++) {
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

vectface DualDescription_incd_mpq_class(MyMatrix<mpq_class> const& TheEXT, F f)
{
  auto f=[](mpq_class const& val) -> mpq_t {
    return val.get_mpq_t();
  };
  return DualDescription_incd_T(TheEXT, f);
}

#ifdef INCLUDE_NUMBER_THEORY_BOOST_GMP_INT

vectface DualDescription_incd_mpq_class(MyMatrix<mpq_class> const& TheEXT)
{
  auto f=[](mpq_class const& val) -> mpq_t {
    return val.backend().data();
  };
  return DualDescription_incd_T(TheEXT, f);
}


#endif

}

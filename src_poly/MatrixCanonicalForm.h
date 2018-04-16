#ifndef MATRIC_CANONICAL_FORM_H
#define MATRIC_CANONICAL_FORM_H


#include "Shvec_double.h"
#include "Temp_PolytopeEquiStab.h"



template<typename T,typename Tint>
MyMatrix<T> ComputeCanonicalForm(MyMatrix<T> const& inpMat)
{
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis(inpMat);
  int nbRow=SHV.rows();
  int n=SHV.cols();
  MyMatrix<Tint> SHVred(nbRow,n);
  std::vector<MyVector<T>> OldBasis;
  std::vector<MyVector<T>> NewBasis;
  int TheDim=0;
  for (int iRow=0; iRow<nbRow; iRow++) {
    MyVector<Tint> eRow=GetMatrixRow(SHV, iRow);
    MyVector<T> eRow_T = ConvertVectorUniversal<T,Tint>(eRow);
    
  }


  
}






#endif

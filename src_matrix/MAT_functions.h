#ifndef MATRIX_FUNCTIONS
#define MATRIX_FUNCTIONS


#include "MAT_Matrix.h"



//
// This routine is adequate for simplifying systems used
// for linear programming.
//
template<typename T>
MyMatrix<T> SortUnicizeMatrix(MyMatrix<T> const& M)
{
  int nbRow=M.rows();
  int nbCol=M.cols();
  auto comp=[](MyVector<T> const& V1, MyVector<T> const& V2) -> bool {
    return IsLower(V1, V2);
  };
  std::set<MyVector<T>, std::function<bool(MyVector<T>,MyVector<T>)> > eListVect(comp);
  //  eListVect=std::set<MyVector<T>, decltype(comp)> (comp);
  for (int iRow=0; iRow<nbRow; iRow++) {
    MyVector<T> V=GetMatrixRow(M, iRow);
    /*    std::cerr << "V=";
    WriteVector(std::cerr, V);
    std::cerr << "Before CanonicalizeVector\n";*/
    MyVector<T> Vcan=CanonicalizeVector(V);
    //    std::cerr << " After CanonicalizeVector\n";
    eListVect.insert(Vcan);
  }
  int nbVect=eListVect.size();
  MyMatrix<T> Mret(nbVect, nbCol);
  int idx=0;
  for (auto & eV : eListVect) {
    for (int iCol=0; iCol<nbCol; iCol++)
      Mret(idx,iCol)=eV(iCol);
    idx++;
  }
  return Mret;
}



#endif

#ifndef DEFINE_INVARIANT_VECTOR_FAMILY_H
#define DEFINE_INVARIANT_VECTOR_FAMILY_H

#include "ShortestUniversal.h"


template<typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamily(MyMatrix<T> const& eMat, std::function<bool(MyMatrix<Tint> const&)> const& fCorrect)
{
  T MaxNorm=MaximumDiagonal(eMat);
  MyMatrix<Tint> SHVall=T_ShortVector<T,Tint>(eMat, MaxNorm);
  //  std::cerr << "MaxNorm = " << MaxNorm << " eMat =\n";
  //  WriteMatrix(std::cerr, eMat);
  //  std::cerr << "|SHVall|=" << SHVall.rows() << "\n";
  //  WriteMatrix(std::cerr, SHVall);
  std::set<T> SetNorm;
  int nbSHV=SHVall.rows();
  std::vector<T> ListNorm(nbSHV);
  for (int iSHV=0; iSHV<nbSHV; iSHV++) {
    MyVector<Tint> eRow=GetMatrixRow(SHVall, iSHV);
    T eNorm=EvaluationQuadForm(eMat, eRow);
    ListNorm[iSHV]=eNorm;
    SetNorm.insert(eNorm);
  }
  std::vector<MyVector<Tint>> ListVect;
  for (auto const& eNorm : SetNorm) {
    //    std::cerr << "eNorm=" << eNorm << "\n";
    for (int iSHV=0; iSHV<nbSHV; iSHV++)
      if (ListNorm[iSHV] == eNorm) {
	MyVector<Tint> eRow=GetMatrixRow(SHVall, iSHV);
	ListVect.push_back(eRow);
      }
    MyMatrix<Tint> SHVret=MatrixFromVectorFamily(ListVect);
    if (fCorrect(SHVret))
      return SHVret;
  }
  std::cerr << "We should never reach that stage\n";
  throw TerminalException{1};
}

template<typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyFullRank(MyMatrix<T> const& eMat)
{
  int n=eMat.rows();
  std::function<bool(MyMatrix<Tint> const&)> fCorrect=[&](MyMatrix<Tint> const& M) -> bool {
    if (RankMat(M) == n)
      return true;
    return false;
  };
  return ExtractInvariantVectorFamily(eMat, fCorrect);
}





template<typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyZbasis_Kernel(MyMatrix<T> const& eMat)
{
  int n=eMat.rows();
  std::function<bool(MyMatrix<Tint> const&)> fCorrect=[&](MyMatrix<Tint> const& M) -> bool {
    if (RankMat(M) < n)
      return false;
    //    std::cerr << "Before Int_IndexLattice computation\n";
    Tint indx=Int_IndexLattice(M);
    //    std::cerr << "indx=" << indx << "\n";
    if (T_Norm(indx) == 1)
      return true;
    return false;
  };
  return ExtractInvariantVectorFamily(eMat, fCorrect);
}




template<typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyZbasis(MyMatrix<T> const& eMat)
{
  LLLreduction<T,Tint> recLLL = LLLreducedBasis<T,Tint>(eMat);
  //  std::cerr << "recLLL.GramMatRed=\n";
  //  WriteMatrix(std::cerr, recLLL.GramMatRed);
  MyMatrix<T> Pmat_T = ConvertMatrixUniversal<T,Tint>(recLLL.Pmat);
  //
  //  std::cerr << "eMat=\n";
  //  WriteMatrix(std::cerr, eMat);
  MyMatrix<T> eMatRed=Pmat_T * eMat * TransposedMat(Pmat_T);
  //  std::cerr << "eMatRed=\n";
  //  WriteMatrix(std::cerr, eMatRed);
  //
  MyMatrix<Tint> SHVred = ExtractInvariantVectorFamilyZbasis_Kernel<T,Tint>(eMatRed);
  MyMatrix<Tint> SHVret = SHVred * recLLL.Pmat;
  return SHVret;
}



#endif

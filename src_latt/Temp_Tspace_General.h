#ifndef TEMP_TSPACE_GENERAL_INCLUDE
#define TEMP_TSPACE_GENERAL_INCLUDE

#include "POLY_PolytopeFct.h"
#include "Temp_PolytopeEquiStab.h"
#include "Temp_Positivity.h"
#include "ShortestUniversal.h"


// Here we use SuperMat as an array because if we declare
// it directly then we need to give a size and we do not
// want that
template<typename T>
struct LinSpaceMatrix {
  int n;
  MyMatrix<T> SuperMat;
  std::vector<MyMatrix<T>> ListMat;
  std::vector<MyMatrix<T>> ListComm;
};


template<typename T>
std::istream& operator>>(std::istream& is, LinSpaceMatrix<T>& obj)
{
  MyMatrix<T> SuperMat=ReadMatrix<T>(is);
  int n=SuperMat.rows();
  //
  int nbMat;
  is >> nbMat;
  std::vector<MyMatrix<T>> ListMat(nbMat);
  for (int iMat=0; iMat<nbMat; iMat++) {
    MyMatrix<T> eMat=ReadMatrix<T>(is);
    ListMat[iMat]=eMat;
  }
  //
  int nbComm;
  is >> nbComm;
  std::vector<MyMatrix<T>> ListComm(nbComm);
  for (int iComm=0; iComm<nbComm; iComm++) {
    MyMatrix<T> eComm=ReadMatrix<T>(is);
    ListComm[iComm]=eComm;
  }
  obj={n, SuperMat, ListMat, ListComm};
  return is;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, LinSpaceMatrix<T> const& obj)
{
  WriteMatrix(os, obj.SuperMat);
  //
  int nbMat=obj.ListMat.size();
  os << " " << nbMat << "\n";
  for (int iMat=0; iMat<nbMat; iMat++)
    WriteMatrix(os, obj.ListMat[iMat]);
  //
  int nbComm=obj.ListComm.size();
  os << " " << nbComm << "\n";
  for (int iComm=0; iComm<nbComm; iComm++)
    WriteMatrix(os, obj.ListComm[iComm]);
  return os;
}




template<typename T>
LinSpaceMatrix<T> ComputeCanonicalSpace(int const& n)
{
  std::vector<MyMatrix<T>> ListMat;
  int i, j;
  T eAss;
  eAss=1;
  for (i=0; i<n; i++)
    for (j=i; j<n; j++) {
      MyMatrix<T> eMat(n, n);
      ZeroAssignation(eMat);
      eMat(i, j)=eAss;
      eMat(j, i)=eAss;
      ListMat.push_back(eMat);
    }
  MyMatrix<T> SuperMat(n, n);
  ZeroAssignation(SuperMat);
  for (i=0; i<n; i++)
    SuperMat(i, i)=eAss;
  return {n, SuperMat, ListMat, {}};
}


template<typename T>
MyMatrix<T> __RealQuadMatSpace(MyMatrix<T> const& eMatB, MyMatrix<T> const& eMatC, int n, T const& eSum, T const& eProd)
{
  int i2, j2;
  T eVal1, eVal1b, eVal2;
  T bVal, cVal;
  MyMatrix<T> eMatN(2*n, 2*n);
  ZeroAssignation(eMatN);
  for (i2=0; i2<n; i2++)
    for (j2=0; j2<n; j2++) {
      bVal=eMatB(i2, j2);
      cVal=eMatC(i2, j2);
      eVal1=2*bVal + eSum*cVal;
      eVal1b=(eSum*eSum - 2*eProd)*bVal + eSum*eProd*cVal;
      eMatN(i2  , j2  )=eVal1;
      eMatN(i2+n, j2+n)=eVal1b;
      eVal2=eSum*bVal + 2*eProd*cVal;
      eMatN(i2  , j2+n)=eVal2;
      eMatN(i2+n, j2  )=eVal2;
    }
  return eMatN;
}


template<typename T>
LinSpaceMatrix<T> ComputeRealQuadraticSpace(int n, T const& eSum, T const& eProd)
{
  std::vector<MyMatrix<T>> ListMat;
  int nbMat, i, j;
  T eOne;
  MyMatrix<T> eMatB(n,n);
  MyMatrix<T> eMatC(n,n);
  nbMat=n*(n+1);
  eOne=1;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      ZeroAssignation(eMatB);
      ZeroAssignation(eMatC);
      eMatB(i, j)=eOne;
      eMatB(j, i)=eOne;
      //
      MyMatrix<T> eMat=__RealQuadMatSpace(eMatB, eMatC, n, eSum, eProd);
      ListMat.push_back(eMat);
      //
      ZeroAssignation(eMatB);
      ZeroAssignation(eMatC);
      eMatC(i, j)=eOne;
      eMatC(j, i)=eOne;
      //
      MyMatrix<T> fMat=__RealQuadMatSpace(eMatB, eMatC, n, eSum, eProd);
      ListMat.push_back(fMat);
    }
  eMatB=IdentityMat<int>(n);
  ZeroAssignation(eMatC);
  MyMatrix<T> SuperMat=__RealQuadMatSpace(eMatB, eMatC, n, eSum, eProd);
  MyMatrix<T> eComm(2*n, 2*n);
  for (i=0; i<n; i++) {
    eComm(i, n+i)=eOne;
    eComm(n+i, i)=-eProd;
    eComm(n+i, n+i)=eSum;
  }
  return {2*n, SuperMat, ListMat, {eComm}};
}





template<typename T>
LinSpaceMatrix<T> ComputeImagQuadraticSpace(int n, T const& eSum, T const& eProd)
{
  std::vector<MyMatrix<T>> ListMat;
  int i, j;
  T eOne, eVal;
  T Discriminant;
  Discriminant=eSum*eSum - 4*eProd;
  if (Discriminant >= 0) {
    std::cerr << "Discriminant is positive\n";
    std::cerr << "Impossible for an imaginary space\n";
    throw TerminalException{1};
  }
  eOne=1;
  for (i=0; i<n; i++)
    for (j=i; j<n; j++) {
      MyMatrix<T> eMat(2*n, 2*n);
      ZeroAssignation(eMat);
      eMat(i, j)=eOne;
      eMat(j, i)=eOne;
      eVal=eSum/2;
      eMat(n+i, j)=eVal;
      eMat(n+j, i)=eVal;
      eMat(i, n+j)=eVal;
      eMat(j, n+i)=eVal;
      //
      eMat(n+i, n+j)=eProd;
      eMat(n+j, n+i)=eProd;
      ListMat.push_back(eMat);
    }
  for (i=0; i<n-1; i++)
    for (j=i+1; j<n; j++) {
      MyMatrix<T> eMat(2*n, 2*n);
      ZeroAssignation(eMat);
      eMat(n+i, j)=eOne;
      eMat(j, n+i)=eOne;
      eVal=-1;
      eMat(n+j, i)=eVal;
      eMat(i, n+j)=eVal;
      ListMat.push_back(eMat);
    }
  MyMatrix<T> SuperMat(2*n, 2*n);
  ZeroAssignation(SuperMat);
  for (i=0; i<n; i++) {
    SuperMat(i, i)=eOne;
    eVal=eSum/2;
    SuperMat(n+i, i)=eVal;
    SuperMat(i, n+i)=eVal;
    //
    SuperMat(n+i, n+i)=eProd;
  }
  MyMatrix<T> eComm(2*n, 2*n);
  ZeroAssignation(eComm);
  for (i=0; i<n; i++) {
    eComm(i, n+i)=eOne;
    eVal=-eProd;
    eComm(n+i, i)=eVal;
    eComm(n+i, n+i)=eSum;
  }
  return {2*n, SuperMat, ListMat, {eComm}};
}


template<typename T>
MyMatrix<T> GetMatrixFromBasis(std::vector<MyMatrix<T>> const& ListMat, MyVector<T> const& eVect)
{
  int n=ListMat[0].rows();
  MyMatrix<T> RetMat=ZeroMatrix<T>(n, n);
  int nbMat=ListMat.size();
  int siz=eVect.size();
  if (siz != nbMat) {
    std::cerr << "Error in GetMatrixFromBasis\n";
    std::cerr << "  siz=" << siz << "\n";
    std::cerr << "nbMat=" << nbMat << "\n";
    std::cerr << "But they should be both equal\n";
    throw TerminalException{1};
  }
  for (int iMat=0; iMat<nbMat; iMat++)
    RetMat += eVect(iMat) * ListMat[iMat];
  return RetMat;
}



template<typename T>
MyMatrix<T> LINSPA_GetMatrixInTspace(LinSpaceMatrix<T> const& LinSpa, MyVector<T> const& eVect)
{
  return GetMatrixFromBasis(LinSpa.ListMat, eVect);
}


template<typename T>
MyVector<T> LINSPA_GetVectorOfMatrixExpression(LinSpaceMatrix<T> const& LinSpa, MyMatrix<T> const& eMat)
{
  MyVector<T> eMatVect=SymmetricMatrixToVector(eMat);
  int dimSymm=eMatVect.size();
  int dimLinSpa=LinSpa.ListMat.size();
  MyMatrix<T> TotalMatrix(dimLinSpa, dimSymm);
  for (int iLinSpa=0; iLinSpa<dimLinSpa; iLinSpa++) {
    MyVector<T> V=SymmetricMatrixToVector(LinSpa.ListMat[iLinSpa]);
    AssignMatrixRow(TotalMatrix, iLinSpa, V);
  }
  SolMatResult<T> RecSol=SolutionMat(TotalMatrix, eMatVect);
  if (!RecSol.result) {
    std::cerr << "The matrix does not belong to the linear space of matrix. Exclude it\n";
    throw TerminalException{1};
  }
  return RecSol.eSol;
}





template<typename T>
T T_GRAM_GetUpperBound(MyMatrix<T> const&TheMat)
{
  T TheLowEst=0;
  int n=TheMat.rows();
  for (int i=0; i<n; i++) {
    T eVal=TheMat(i, i);
    if (eVal > TheLowEst)
      TheLowEst=eVal;
  }
  T MaxNorm=TheLowEst;
  return MaxNorm;
}


template<typename T>
void T_UpdateListValue(MyMatrix<T> const&eMat, T const & TheTol, std::vector<T> &ListVal)
{
  int nbRow=eMat.nbRow;
  int nbCol=eMat.nbCol;
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      bool WeFound=false;
      int nb=ListVal.size();
      T fVal=eMat(iRow, iCol);
      for (int iVal=0; iVal<nb; iVal++) {
	T hVal=ListVal[iVal];
	if (T_abs(hVal - fVal) < TheTol)
	  WeFound=true;
      }
      if (!WeFound)
	ListVal.push_back(fVal);
    }
}


template<typename T>
MyMatrix<T> T_GRAM_GetScalProdMat(MyMatrix<T> const& eMat, MyMatrix<int> const& ListShort)
{
  int nbShort=ListShort.rows();
  int n=ListShort.cols();
  MyMatrix<T> ScalProdMat(nbShort, nbShort);
  for (int iShort=0; iShort<nbShort; iShort++)
    for (int jShort=0; jShort<nbShort; jShort++) {
      T eScal=0;
      for (int i=0; i<n; i++)
	for (int j=0; j<n; j++) {
	  int eVal12=ListShort(jShort, j) * ListShort(iShort, i);
	  T eVal=eMat(i, j);
	  eScal += eVal12*eVal;
	}
      ScalProdMat(iShort, jShort)=eScal;
    }
  return ScalProdMat;
}


template<typename T, typename Tint, typename Tgroup>
void T_GetGramMatrixAutomorphismGroup(MyMatrix<T> const& eMat, T const& TheTol, Tgroup & GRPperm, std::vector<MyMatrix<Tint>> &ListMatrGens)
{
  using Tidx_value = int16_t;
  T MaxDet=T_GRAM_GetUpperBound(eMat);
  MyMatrix<Tint> ListShort=T_ShortVector(eMat, MaxDet);
  MyMatrix<T> ListShort_T=UniversalMatrixConversion<T,Tint>(ListShort);
  WeightMatrix<true, T, Tidx_value> WMat=T_TranslateToMatrix_QM_SHV<T,Tint,Tidx_value>(eMat, ListShort);
  GRPperm=GetStabilizerWeightMatrix(WMat);
  ListMatrGens.clear();
  for (auto & eGen : GRPperm.group->S) {
    MyMatrix<T> M3_T=RepresentVertexPermutation(ListShort_T, ListShort_T, *eGen);
    MyMatrix<Tint> M3_I=UniversalMatrixConversion<Tint,T>(M3_T);
    ListMatrGens.push_back(M3_I);
  }
}

template<typename T, typename Tint, typename Telt>
bool T_TestGramMatrixEquivalence(MyMatrix<T> const& eMat1, MyMatrix<T> const& eMat2, T const &TheTol)
{
  using Tidx_value = int16_t;
  T MaxDet1=T_GRAM_GetUpperBound(eMat1);
  T MaxDet2=T_GRAM_GetUpperBound(eMat2);
  T MaxDet=MaxDet1;
  if (MaxDet1 > MaxDet2)
    MaxDet=MaxDet2;
  MyMatrix<Tint> ListShort1=T_ShortVector(eMat1, MaxDet);
  MyMatrix<Tint> ListShort2=T_ShortVector(eMat2, MaxDet);
  WeightMatrix<true, T, Tidx_value> WMat1 = T_TranslateToMatrix_QM_SHV<T,Tint,Tidx_value>(eMat1, ListShort1);
  WeightMatrix<true, T, Tidx_value> WMat2 = T_TranslateToMatrix_QM_SHV<T,Tint,Tidx_value>(eMat2, ListShort2);
  EquivTest<Telt> eResEquiv=TestEquivalenceWeightMatrix<T,Telt>(WMat1, WMat2);
  return eResEquiv.TheReply;
}


//
// We search for the set of matrices satisfying g M g^T = M for all g in ListGen
//
template<typename T>
std::vector<MyMatrix<T>> BasisInvariantForm(int const& n, std::vector<MyMatrix<T>> const& ListGen)
{
  std::vector<std::vector<int>> ListCoeff;
  for (int iLin=0; iLin<n; iLin++)
    for (int iCol=0; iCol<=iLin; iCol++)
      ListCoeff.push_back({iLin,iCol});
  int nbCoeff=ListCoeff.size();
  auto FuncPos=[&](int const&i, int const&j) -> int {
    if (i < j)
      return PositionVect(ListCoeff, {j,i});
    return PositionVect(ListCoeff, {i,j});
  };
  std::vector<MyVector<T>> ListEquations;
  for (auto & eGen : ListGen)
    for (int i=0; i<n; i++)
      for (int j=0; j<=i; j++) {
	MyVector<T> TheEquation=ZeroVector<T>(nbCoeff);
	TheEquation(FuncPos(i,j)) += 1;
	for (int k=0; k<n; k++)
	  for (int l=0; l<n; l++) {
	    int pos=FuncPos(k,l);
	    TheEquation(pos) += - eGen(i,k)*eGen(j,l);
	  }
	ListEquations.push_back(TheEquation);
      }
  MyMatrix<T> MatEquations=MatrixFromVectorFamily(ListEquations);
  std::cerr << "Before call to NullspaceTrMat nbGen=" << ListGen.size() << " MatEquations(rows/cols)=" << MatEquations.rows() << " / " << MatEquations.cols() << "\n";
  MyMatrix<T> NSP=NullspaceTrMat(MatEquations);
  std::cerr << "After call to NullspaceTrMat NSP.rows=" << NSP.rows() << " cols=" << NSP.cols() << "\n";
  int dimSpa=NSP.rows();
  std::vector<MyMatrix<T>> TheBasis(dimSpa);
  for (int iDim=0; iDim<dimSpa; iDim++) {
    MyVector<T> eRow=GetMatrixRow(NSP,iDim);
    MyMatrix<T> eMat(n,n);
    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
	eMat(i,j)=eRow(FuncPos(i,j));
    TheBasis[iDim]=eMat;
  }
  return TheBasis;
}


template<typename T>
struct PairV {
  T val;
  int nb;
};

template<typename T>
bool operator==(PairV<T> const& x, PairV<T> const& y)
{
  if (x.val != y.val)
    return false;
  if (x.nb != y.nb)
    return false;
  return true;
}

template<typename T>
bool operator<(PairV<T> const& x, PairV<T> const& y)
{
  if (x.val < y.val)
    return true;
  if (x.val > y.val)
    return false;
  return x.nb < y.nb;
}

template<typename T>
bool operator>(PairV<T> const& x, PairV<T> const& y)
{
  return y<x;
}



template<typename T>
struct GramInvVect {
  T eDet;
  std::vector<PairV<T>> ListValOff;
  std::vector<PairV<T>> ListValDiag;
};


template<typename T>
bool operator<(GramInvVect<T> const& x, GramInvVect<T> const& y)
{
  if (x.eDet < y.eDet)
    return true;
  if (x.eDet > y.eDet)
    return false;
  int xSizOff=x.ListValOff.size();
  int ySizOff=y.ListValOff.size();
  if (xSizOff < ySizOff)
    return true;
  if (xSizOff > ySizOff)
    return false;
  for (int i=0; i<xSizOff; i++) {
    if (x.ListValOff[i] < y.ListValOff[i])
      return true;
    if (x.ListValOff[i] > y.ListValOff[i])
      return false;
  }
  int xSizDiag=x.ListValDiag.size();
  int ySizDiag=y.ListValDiag.size();
  if (xSizDiag < ySizDiag)
    return true;
  if (xSizDiag > ySizDiag)
    return false;
  for (int i=0; i<xSizDiag; i++) {
    if (x.ListValDiag[i] < y.ListValDiag[i])
      return true;
    if (x.ListValDiag[i] > y.ListValDiag[i])
      return false;
  }
  return false;
}

template<typename T>
std::istream& operator>>(std::istream& is, GramInvVect<T>& obj)
{
  T eDet;
  is >> eDet;
  int nbOff;
  is >> nbOff;
  std::vector<PairV<T>> ListValOff(nbOff);
  for (int i=0; i<nbOff; i++) {
    T eVal;
    int nb;
    is >> eVal;
    is >> nb;
    ListValOff[i]={eVal,nb};
  }
  int nbDiag;
  is >> nbDiag;
  std::vector<PairV<T>> ListValDiag(nbDiag);
  for (int i=0; i<nbOff; i++) {
    T eVal;
    int nb;
    is >> eVal;
    is >> nb;
    ListValDiag[i]={eVal,nb};
  }
  obj={eDet, ListValOff, ListValDiag};
  return is;
}


template<typename T>
std::ostream& operator<<(std::ostream& os, GramInvVect<T> const& obj)
{
  os << obj.eDet << "\n";
  int nbOff=obj.ListValOff.size();
  os << nbOff << "\n";
  for (int i=0; i<nbOff; i++)
    os << obj.ListValOff[i].val << " " << obj.ListValOff[i].nb << "\n";
  int nbDiag=obj.ListValDiag.size();
  os << nbDiag << "\n";
  for (int i=0; i<nbDiag; i++)
    os << obj.ListValDiag[i].val << " " << obj.ListValDiag[i].nb << "\n";
  return os;
}

template<typename T>
bool operator==(GramInvVect<T> const& x, GramInvVect<T> const& y)
{
  if (x.eDet != y.eDet)
    return false;
  if (x.ListValOff != y.ListValOff)
    return false;
  if (x.ListValDiag != y.ListValDiag)
    return false;
  return true;
}










template<typename T, typename Tint>
GramInvVect<T> GetInvariantGram(MyMatrix<T> const& eGram, MyMatrix<Tint> const& SHV)
{
  Tshortest<T,Tint> eRec=T_ShortestVector<T,Tint>(eGram);
  MyMatrix<T> eGramRef=eGram / eRec.eMin;
  T eDet=DeterminantMat(eGramRef);
  int nbVect=SHV.rows();
  //  std::cerr << "eDet=" << eDet << " nbVect=" << nbVect << "\n";
  std::vector<T> ListValOff;
  std::vector<int> ListNbOff;
  std::vector<T> ListValDiag;
  std::vector<int> ListNbDiag;
  for (int iVect=0; iVect<nbVect; iVect++) {
    MyVector<Tint> eRow=GetMatrixRow(SHV,iVect);
    T eNorm=EvaluationQuadForm(eGramRef, eRow);
    int pos=PositionVect(ListValDiag, eNorm);
    if (pos == -1) {
      ListValDiag.push_back(eNorm);
      ListNbDiag.push_back(1);
    }
    else {
      ListNbDiag[pos]++;
    }
    for (int jVect=iVect+1; jVect<nbVect; jVect++) {
      MyVector<Tint> fRow=GetMatrixRow(SHV, jVect);
      T eScal=ScalarProductQuadForm(eGramRef, eRow, fRow);
      int posB=PositionVect(ListValOff, eScal);
      //      std::cerr << "iVect=" << iVect << " jVect=" << jVect << " eScal=" << eScal << " posB=" << posB << "\n";
      if (posB == -1) {
	ListValOff.push_back(eScal);
	ListNbOff.push_back(1);
      }
      else {
	ListNbOff[posB]++;
      }
    }
  }
  /*
  std::cerr << "sum(ListNbOff)=" << VectorSum(ListNbOff) << "\n";
  std::cerr << "LComb =";
  for (int u=0; u<int(ListNbOff.size()); u++) {
    std::cerr << " [" << ListValOff[u] << "," << ListNbOff[u] << "]";
  }
  std::cerr << "\n";*/
  std::function<bool(PairV<T> const&,PairV<T> const&)> Comp=[](PairV<T> const& x, PairV<T> const& y) -> bool {
    if (x.val < y.val)
      return true;
    if (x.val > y.val)
      return false;
    return x.nb < y.nb;
  };
  std::set<PairV<T>,std::function<bool(PairV<T> const&,PairV<T> const&)>> ListValOff_B(Comp);
  std::set<PairV<T>,std::function<bool(PairV<T> const&,PairV<T> const&)>> ListValDiag_B(Comp);
  int nbOff=ListValOff.size();
  int nbDiag=ListValDiag.size();
  for (int i=0; i<nbOff; i++)
    ListValOff_B.insert({ListValOff[i], ListNbOff[i]});
  for (int i=0; i<nbDiag; i++)
    ListValDiag_B.insert({ListValDiag[i], ListNbDiag[i]});
  std::vector<PairV<T>> ListValOff_C;
  for (auto & x : ListValOff_B)
    ListValOff_C.push_back(x);
  std::vector<PairV<T>> ListValDiag_C;
  for (auto & x : ListValDiag_B)
    ListValDiag_C.push_back(x);
  return {eDet, ListValOff_C, ListValDiag_C};
}



template<typename T>
std::vector<MyMatrix<T>> StandardSymmetricBasis(int n)
{
  std::vector<MyMatrix<T>> ListMat;
  T eOne=1;
  for (int i=0; i<n; i++)
    for (int j=0; j<=i; j++) {
      MyMatrix<T> eMat=ZeroMatrix<T>(n,n);
      eMat(i,j)=eOne;
      eMat(j,i)=eOne;
      ListMat.push_back(eMat);
    }
  return ListMat;
}


#endif

#ifndef TEMP_ENUMERATION_CONFIG
#define TEMP_ENUMERATION_CONFIG

#include "MAT_MatrixInt.h"
#include "Shvec_double.h"
#include "Temp_Positivity.h"

template<typename T>
T GetHermitePower(int const& n)
{
  if (n == 1)
    return 1;
  if (n == 2)
    return 4/3;
  if (n == 3)
    return 2;
  if (n == 4)
    return 4;
  if (n == 5)
    return 8;
  if (n == 6)
    return 64/3;
  if (n == 7)
    return 64;
  if (n == 8)
    return 256;
  T nDoubl=n;
  T h=nDoubl*(nDoubl-1)/2;
  T eVal=((double)64)/((double)3);
  return pow(eVal, h);
}



template<typename T>
MyMatrix<T> GetMatrixAv(MyMatrix<T> const& TheMat, MyVector<int> const& TheV, MyMatrix<int> const& TheCompl)
{
  int n=TheMat.rows();
  T rNorm=EvaluationQuadForm<T,int>(TheMat, TheV);
  MyMatrix<T> TheProj=MyMatrix<T>(n-1,n);
  for (int j=0; j<n-1; j++) {
    MyVector<int> eCompl=TheCompl.row(j);
    T eScal=ScalarProductQuadForm<T,int>(TheMat, TheV, eCompl);
    for (int i=0; i<n; i++)
      {
	T eVal1=TheCompl(j,i);
	T eVal2=TheV(i);
	T hVal=eVal1 - (eScal/rNorm)*eVal2;
	TheProj(j,i)=hVal;
      }
  }
  MyMatrix<T> RetMat(n-1,n-1);
  for (int i=0; i<n-1; i++) {
    for (int j=0; j<n-1; j++) {
      MyVector<T> V1=TheProj.row(i);
      MyVector<T> V2=TheProj.row(i);
      T eScal=ScalarProductQuadForm(TheMat, V1, V2);
      RetMat(i,j)=eScal;
    }
  }
  return RetMat;
}


template<typename T>
void EnumerationKspaces(MyMatrix<T> const& TheMat, int const& k, T const& MaxDet, std::vector<MyMatrix<int> > &ListSubspace, std::vector<T> &ListDet)
{
  ListSubspace.clear();
  ListDet.clear();
  int n=TheMat.rows();
  //  fprintf(stdout, "EnumerationKspace k=%d  MaxDet=%lg\n", k, MaxDet);
  //  PrintMat_Doubl(stdout, TheMat);
  if (k == 1) {
    MyMatrix<int> ListShort=T_ShortVector(TheMat, MaxDet);
    int nbShort=ListShort.rows();
    for (int iShort=0; iShort<nbShort; iShort++) {
      MyVector<int> eShort=ListShort.row(iShort);
      if (IsVectorPrimitive(eShort) == 1) {
	MyMatrix<int> eSubspace=MyMatrix<int>(1,n);
	AssignMatrixRow(eSubspace, 0, eShort);
	ListSubspace.push_back(eSubspace);
      }
    }
    return;
  }
  T kDoubl=k;
  T eExpo=GetHermitePower<T>(k)*MaxDet;
  T MaxNorm=pow(eExpo, 1/kDoubl);
  MyMatrix<int> ListShort=T_ShortVector(TheMat, MaxNorm);
  int nbShort=ListShort.rows();
  for (int j=0; j<nbShort; j++) {
    MyVector<int> eShort=ListShort.row(j);
    if (IsVectorPrimitive(eShort) == 1) {
      T eNorm=EvaluationQuadForm<T,int>(TheMat, eShort);
      MyMatrix<int> TheCompl=ComplementToBasis<int>(eShort);
      MyMatrix<T> RetMat=GetMatrixAv(TheMat, eShort, TheCompl);
      T NewMaxDet=MaxDet/eNorm;
      std::vector<MyMatrix<int> > RedListSubspace;
      std::vector<T> RedListDet;
      EnumerationKspaces(RetMat, k-1, NewMaxDet, RedListSubspace, RedListDet);
      int nbSubspace=RedListSubspace.size();
      for (int iSub=0; iSub<nbSubspace; iSub++) {
	MyMatrix<int> eSubspace=RedListSubspace[iSub];
	MyMatrix<int> NewSubspace=MyMatrix<int>(k,n);
	AssignMatrixRow(NewSubspace, 0, eShort);
	for (int i=0; i<k-1; i++)
	  for (int r=0; r<n; r++) {
	    int eSum=0;
	    for (int l=0; l<n-1; l++) {
	      int eVal1=eSubspace(i,l);
	      int eVal2=TheCompl(l,r);
	      eSum += eVal1*eVal2;
	    }
	    NewSubspace(i+1,r)=eSum;
	  }
	T NewDet=eNorm*RedListDet[iSub];
	FuncInsertSubspace(ListSubspace, ListDet, NewSubspace, NewDet);
      }
    }
  }
}


template<typename T>
void UpperEstimateKMinimum(MyMatrix<T> const& TheMat, int const& k, MyMatrix<int> & OneSubspace, T &TheEstimate)
{
  int n=TheMat.rows();
  T MinDiag=TheMat(0,0);
  for (int i=1; i<n; i++) {
    T eVal=TheMat(i,i);
    if (eVal < MinDiag)
      MinDiag=eVal;
  }
  MyMatrix<int> ListShort=T_ShortVector(TheMat, MinDiag);
  int nbShort=ListShort.rows();
  int IsFirst=1;
  T MinNorm=0;
  MyVector<int> FoundMinVect=MyVector<int>(n);
  for (int iShort=0; iShort<nbShort; iShort++) {
    MyVector<int> eShort=ListShort.row(iShort);
    T rNorm=EvaluationQuadForm(TheMat, eShort);
    if (IsFirst == 1) {
      IsFirst=0;
      MinNorm=rNorm;
      FoundMinVect=eShort;
    }
    else {
      if (rNorm < MinNorm) {
	MinNorm=rNorm;
	FoundMinVect=eShort;
      }
    }
  }
  if (IsVectorPrimitive(FoundMinVect) == 0) {
    std::cerr << "FoundMinVect=";
    WriteVector(std::cerr, FoundMinVect);
    std::cerr << "TheMat=\n";
    WriteMatrix(std::cerr, TheMat);
    std::cerr << "Deep incoherency, PANIC\n";
    throw TerminalException{1};
  }
  if (k == 1) {
    TheEstimate=MinNorm;
    return;
  }
  MyMatrix<int> TheCompl=ComplementToBasis<int>(FoundMinVect);
  MyMatrix<T> RetMat=GetMatrixAv(TheMat, FoundMinVect, TheCompl);
  MyMatrix<int> RedSubspace=MyMatrix<int>(k-1,n);
  T RedEstimate;
  UpperEstimateKMinimum(RetMat, k-1, RedSubspace, RedEstimate);
  AssignMatrixRow(OneSubspace, 0, FoundMinVect);
  for (int i=0; i<k-1; i++)
    for (int r=0; r<n; r++) {
      int eSum=0;
      for (int l=0; l<n-1; l++) {
	int eVal1=RedSubspace(i,l);
	int eVal2=TheCompl(l,r);
	eSum += eVal1*eVal2;
      }
      OneSubspace(i+1,r)=eSum;
    }
  TheEstimate=RedEstimate*MinNorm;
}


#endif

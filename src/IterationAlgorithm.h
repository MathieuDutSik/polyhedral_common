#ifndef ITERATION_ALGORITHM_VORONOI_RANKINE
#define ITERATION_ALGORITHM_VORONOI_RANKINE

#include "EnumerationConfig.h"
#include "MAT_MatrixFLT.h"
#include "Temp_Tspace_General.h"

template<typename T>
struct AllTol {
  int nbIteration;
  T tolDirMat;
  T deltaInc;
  T IsomTol;
  T MaxDetTol;
  T MaxDetTargetTol;
  T MaxMoveError;
};

                                        


template<typename T>
MyMatrix<T> EvaluatorMatrix(MyMatrix<T> const& RetMat,
			    MyMatrix<int> const& eSubspace, 
			    T & TheDet)
{
  MyMatrix<int> eSubspaceTr=eSubspace.transpose();
  MyMatrix<T> eProd1=MatrixProduct_T_Oth(RetMat, eSubspaceTr);
  MyMatrix<T> eProd2=MatrixProduct_Oth_T(eSubspace, eProd1);
  MyMatrix<T> eInv=Inverse(eProd2);
  MyMatrix<T> eProd3=MatrixProduct_Oth_T(eSubspaceTr, eInv);
  MyMatrix<T> eEvalMat=MatrixProduct_T_Oth(eProd3, eSubspace);
  TheDet=DeterminantMat(eProd2);
  return eEvalMat;
}

template<typename T>
T ComputeError(MyMatrix<T> const& InpMat, std::vector<MyMatrix<int> > const& ListSubspace)
{
  T TotalError=0;
  int nbSpace=ListSubspace.size();
  for (int iSub=0; iSub<nbSpace; iSub++) {
    T TheDet;
    MyMatrix<T> eEvalMat=EvaluatorMatrix(InpMat, ListSubspace[iSub], TheDet);
    TotalError += std::abs(TheDet-1);
  }
  return TotalError;
}

template<typename T>
T ComputeErrorDebug(MyMatrix<T> const& InpMat, std::vector<MyMatrix<int> > const& ListSubspace)
{
  T TotalError=0;
  int nbSpace=ListSubspace.size();
  for (int iSub=0; iSub<nbSpace; iSub++) {
    T TheDet;
    MyMatrix<T> eValMat=EvaluatorMatrix(InpMat, ListSubspace[iSub], TheDet);
    std::cerr << "iSub=" << iSub << " TheDet=" << TheDet << "\n";
    TotalError += std::abs(TheDet-1);
  }
  return TotalError;
}

template<typename T>
MyMatrix<T> ReductionMatrix(MyMatrix<T> const& InpMat, std::vector<MyMatrix<T> > const& ListDirMat, std::vector<MyMatrix<int> > const& ListSubspace, int nbIteration)
{
  int n=InpMat.rows();
  int nbSpace=ListSubspace.size();
  int nbDir=ListDirMat.size();
  MyVector<T> ListDiff(nbSpace);
  MyMatrix<T> MatScal=ZeroMatrix<T>(nbDir, nbSpace);
  MyMatrix<T> MatGram=ZeroMatrix<T>(nbDir, nbDir);

  MyVector<T> bSide=MyVector<T>(nbDir);
  MyMatrix<T> RetMat=MyMatrix<T>(n,n);
  RetMat=InpMat;
  for (int iIter=0; iIter<nbIteration; iIter++) {
    T TotalError=0;
    for (int iSub=0; iSub<nbSpace; iSub++) {
      T TheDet;
      MyMatrix<T> eEvalMat=EvaluatorMatrix(RetMat, ListSubspace[iSub], TheDet);
      T eVal=1 - TheDet;
      ListDiff(iSub)=eVal;
      TotalError += std::abs(TheDet-1);
      for (int iDir=0; iDir<nbDir; iDir++)
	MatScal(iDir, iSub)=TheDet*MatrixScalarProduct(eEvalMat, ListDirMat[iDir]);
    }
    for (int iDir=0; iDir<nbDir; iDir++)
      for (int jDir=0; jDir<nbDir; jDir++) {
	MyVector<T> V1=MatScal.row(iDir);
	MyVector<T> V2=MatScal.row(jDir);
	T eScal=ScalarProduct(V1, V2);
	MatGram(iDir, jDir)=eScal;
      }
    for (int iDir=0; iDir<nbDir; iDir++) {
      MyVector<T> V=MatScal.row(iDir);
      T eScal=ScalarProduct(V, ListDiff);
      bSide(iDir)=eScal;
    }
    MyMatrix<T> InvMatGram=Inverse(MatGram);
    MyVector<T> xSol=VectorMatrix(bSide, InvMatGram);
    T ErrSol=0;
    for (int i=0; i<n; i++) {
      T eSum=0;
      for (int j=0; j<n; j++) {
	T eVal=MatGram(i,j);
	eSum += xSol(j)*eVal;
      }
      T eDiff=bSide(i) - eSum;
      ErrSol += std::abs(eDiff);
    }
    ZeroAssignation(RetMat);
    for (int iDir=0; iDir<nbDir; iDir++) {
      T eCoeff=xSol(iDir);
      RetMat += eCoeff*ListDirMat[iDir];
    }
  }
  return RetMat;
}


template<typename T>
void KernelMatrix(MyMatrix<T> const& TheMat,
		  std::vector<MyMatrix<int> > const& ListSubspace,
		  MyMatrix<T> const& DirMat, 
		  std::vector<MyMatrix<T> > & ListOrthBasis,
		  MyMatrix<T> & NewDirMat,
		  MyVector<T> & ListEigVal)
{
  int nbSubspace=ListSubspace.size();
  int n=TheMat.rows();
  int nSym=n*(n+1)/2;
  MyMatrix<T> SpaceMatrix=ZeroMatrix<T>(nSym, nSym);
  for (int iSub=0; iSub<nbSubspace; iSub++) {
    T TheDet;
    MyMatrix<T> eEval=EvaluatorMatrix(TheMat, ListSubspace[iSub], TheDet);
    MyVector<T> TheV=SymmetricMatrixToVector(eEval);
    MyMatrix<T> rnkOne=RankOneMatrix(TheV);
    SpaceMatrix += rnkOne;
  }
  MyVector<T> TheVdir=SymmetricMatrixToVectorB(DirMat);
  //  T eNorm1=ScalarProductQuadForm(SpaceMatrix, TheVdir, TheVdir);
  //  T eNorm2=ScalarProduct(TheVdir, TheVdir);
  MyMatrix<T> ListEigVect=MyMatrix<T>(nSym, nSym);
  jacobi_double(SpaceMatrix, ListEigVal, ListEigVect);
  T MinEigVal=ListEigVal(0);
  int idxSelect=0;
  for (int i=1; i<nSym; i++) {
    T eEig=ListEigVal(i);
    if (eEig < MinEigVal) {
      MinEigVal=eEig;
      idxSelect=i;
    }
  }
  ListOrthBasis.clear();
  for (int i=0; i<nSym; i++) {
    MyVector<T> eEigVect=ListEigVect.row(i);
    MyMatrix<T> eMatB=VectorToSymmetricMatrixB(eEigVect, n);
    if (i == idxSelect)
      NewDirMat=eMatB;
    else
      ListOrthBasis.push_back(eMatB);
  }
  T ScalNN=MatrixScalarProduct(NewDirMat, NewDirMat);
  T ScalNO=MatrixScalarProduct(DirMat, NewDirMat);
  T ScalOO=MatrixScalarProduct(DirMat, DirMat);
  T eQuot=ScalNO/sqrt(ScalNN*ScalOO);
  T eVal;
  if (eQuot < 0)
    eVal=-1/sqrt(ScalNN);
  else
    eVal=1/sqrt(ScalNN);
  NewDirMat=eVal*NewDirMat;
}


template<typename T>
MyMatrix<T> KFlippingOperation(MyMatrix<T> const& OrigMat, 
			       int const& k, 
			       MyMatrix<T> const& DirMat, 
			       AllTol<T> const& RecTol)
{
  std::vector<MyMatrix<int> > NewListSubspace, ListSubspace, ListRelSubspace;
  std::vector<MyMatrix<int> > TestListSubspace, ListSubspaceMin;
  std::vector<T> ListDet, NewListDet, TestListDet;
  int n=OrigMat.rows();
  T kDoubl=k;
  std::cerr << "OrigMat=\n";
  WriteMatrix(std::cerr, OrigMat);
  MyMatrix<int> OneSubspace=MyMatrix<int>(k,n);
  T TheEstimate;
  UpperEstimateKMinimum(OrigMat, k, OneSubspace, TheEstimate);
  std::cerr << "TheEstimate=" << TheEstimate << "\n";
  T MaxDet=2*TheEstimate;
  std::cerr << "MaxDet=" << MaxDet << "\n";
  EnumerationKspaces(OrigMat, k, MaxDet, ListSubspace, ListDet);
  
  T MinDet=ListDet[0];
  int nbSubspace=ListSubspace.size();
  for (int iSub=1; iSub<nbSubspace; iSub++)
    if (ListDet[iSub] < MinDet)
      MinDet=ListDet[iSub];
  std::cerr << "MinDet=" << MinDet << "\n";
  T eScal=pow(MinDet, -1/kDoubl);
  MyMatrix<T> StartMat=eScal*OrigMat;
  T MaxTolDet=MinDet*RecTol.MaxDetTol;
  for (int iSub=0; iSub<nbSubspace; iSub++)
    if (ListDet[iSub] < MaxTolDet)
      {
	T TheDet;
	ListSubspaceMin.push_back(ListSubspace[iSub]);
	MyMatrix<T> eEval=EvaluatorMatrix(OrigMat, ListSubspace[iSub], TheDet);
	T eScalB=MatrixScalarProduct(eEval, DirMat);
	if (eScalB < RecTol.tolDirMat)
	  ListRelSubspace.push_back(ListSubspace[iSub]);
      }
  int nbMinSubspace=ListSubspaceMin.size();
  std::cerr << "nbMinSubspace=" << nbMinSubspace << "\n";
  int nbRelSubspace=ListRelSubspace.size();
  std::cerr << "nbRelSubspace=" << nbRelSubspace << "\n";
  MyMatrix<T> MovingMat=DirMat;
  std::cerr << "MovingMat=\n";
  WriteMatrix(std::cerr, MovingMat);
  MyMatrix<T> RetMat=StartMat;
  int WeFoundNewSubspace=0;
  T deltaIncNew=3; /* irrelevant value, just to avoid warning */
  while(true) {
    T deltaInc=RecTol.deltaInc;
    int nbIteration=RecTol.nbIteration;
    T deltaMove=deltaInc;
    int nbIterReducing=0;
    MyMatrix<T> PreRetMat=MyMatrix<T>(n,n);
    while(true) {
      nbIterReducing++;
      PreRetMat = RetMat + deltaMove * MovingMat;
      T TotalError=ComputeError(PreRetMat, ListRelSubspace);
      std::cerr << "deltaMove=" << deltaMove << " TotalError=" << TotalError << " WeFNewSub=" << WeFoundNewSubspace << " nbIter=" << nbIterReducing << "\n";
      if (TotalError < RecTol.MaxMoveError)
	break;
      deltaMove=deltaMove/2;
      if (nbIterReducing > 1500) {
	std::cerr << "We have a clear bug to solve\n";
	throw TerminalException{1};
      }
    }
    if (WeFoundNewSubspace == 1) {
      std::cerr << "deltaMove=" << deltaMove << " deltaIncNew=" << deltaIncNew << "\n";
      if (deltaIncNew < deltaMove) {
	deltaMove=deltaIncNew;
	MyMatrix<T> PreRetMatB = RetMat + deltaMove*MovingMat;
	std::cerr << "Now doing the canonicalization\n";
	std::vector<MyMatrix<T> > ListDirMat=StandardSymmetricBasis<T>(n);
	MyMatrix<T> FinalMat=ReductionMatrix(PreRetMatB, ListDirMat, NewListSubspace, nbIteration);
	std::cerr << "StartMat=\n";
	WriteMatrix(std::cerr, StartMat);
	std::cerr << "FinalMat=\n";
	WriteMatrix(std::cerr, FinalMat);
	return FinalMat;
      }
    }
    deltaMove=deltaInc;
    std::cerr << "deltaMove=" << deltaMove << "\n";
    int nSym=n*(n+1)/2;
    std::vector<MyMatrix<T> > ListOrthBasis;
    MyMatrix<T> NewDirMat=MyMatrix<T>(n,n);
    MyVector<T> ListEigVal(nSym);
    KernelMatrix(PreRetMat, ListRelSubspace, MovingMat, 
		 ListOrthBasis, NewDirMat, ListEigVal);
    RetMat=ReductionMatrix(PreRetMat, ListOrthBasis, ListRelSubspace, nbIteration);
    MaxDet=RecTol.MaxDetTol;
    EnumerationKspaces(RetMat, k, MaxDet, NewListSubspace, NewListDet);
    int NewNbSubspace=NewListSubspace.size();
    std::cerr << "NewNbSubspace=" << NewNbSubspace << " nbRelSubspace=" << nbRelSubspace << "\n";
    int testEqua=TestInclusionFamilySubspace(ListSubspaceMin, NewListSubspace);
    if (NewNbSubspace > nbRelSubspace && testEqua == 0) {
      std::cerr << "MaxDet=" << MaxDet << "\n";
      WeFoundNewSubspace=1;
      MyMatrix<int> eMatNew=GetNoncontainedSubspace(NewListSubspace, ListSubspaceMin);
      T TheDetNew;
      MyMatrix<T> eEvalNew=EvaluatorMatrix(RetMat, eMatNew, TheDetNew);
      T eScalNew=MatrixScalarProduct(eEvalNew, NewDirMat);
      deltaIncNew=(1-TheDetNew)/eScalNew;
    }
    MovingMat=NewDirMat;
    MaxDet=2;
    EnumerationKspaces(RetMat, k, MaxDet, TestListSubspace, TestListDet);
    ReorderSubspaceDet(TestListSubspace, TestListDet);
  }
}

#endif

#ifndef VORONOI_RANKINE
#define VORONOI_RANKINE

#include "IterationAlgorithm.h"
#include "Temp_PolytopeEquiStab.h"
#include "Temp_Tspace_General.h"
#include "POLY_PolytopeFct.h"
#include "POLY_cddlib.h"

template<typename T>
MyMatrix<T> GetLatticeAdjA4(int n, int k)
{
  T eFact=sqrt(8);
  MyMatrix<T> RetMat=ZeroMatrix<T>(n, n);
  RetMat(0,0)=4/eFact;
  RetMat(0,1)=-2/eFact;

  RetMat(1,0)=-2/eFact;
  RetMat(1,1)=3/eFact;
  RetMat(1,2)=-1/eFact;
  RetMat(1,3)=-1/eFact;

  RetMat(2,1)=-1/eFact;
  RetMat(2,2)=3/eFact;
  RetMat(2,3)=-1/eFact;

  RetMat(3,1)=-1/eFact;
  RetMat(3,2)=-1/eFact;
  RetMat(3,3)=3/eFact;
  return RetMat;
}

template<typename T>
MyMatrix<T> GetLatticePbDim5(int n, int k)
{
  if (n != 5 || k != 2) {
    std::cerr << "Error in function call\n";
    throw TerminalException{1};
  }
  MyMatrix<T> RetMat=ZeroMatrix<T>(n,n);
  T eFact=sqrt(3);

  RetMat(0,0)=2;
  RetMat(0,1)=0;
  RetMat(0,2)=-1;
  RetMat(0,3)=0;
  RetMat(0,4)=0;

  RetMat(1,0)=0;
  RetMat(1,1)=3/2;
  RetMat(1,2)=-3/2;
  RetMat(1,3)=0;
  RetMat(1,4)=0;

  RetMat(2,0)=-1;
  RetMat(2,1)=-3/2;
  RetMat(2,2)=7/2;
  RetMat(2,3)=-1;
  RetMat(2,4)=0;

  RetMat(3,0)=0;
  RetMat(3,1)=0;
  RetMat(3,2)=-1;
  RetMat(3,3)=2;
  RetMat(3,4)=-1;

  RetMat(4,0)=0;
  RetMat(4,1)=0;
  RetMat(4,2)=0;
  RetMat(4,3)=-1;
  RetMat(4,4)=2;
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) {
      T eVal=RetMat(i,j);
      RetMat(i,j)=eVal/eFact;
    }
  return RetMat;
}


template<typename T>
MyMatrix<T> GetLatticeAn(size_t n, size_t k)
{
  MyMatrix<T> RetMat=ZeroMatrix<T>(n,n);
  int eVal;
  T eFact=exp(log(double(k+1))/double(k));
  for (size_t i=0; i<n; i++)
    for (size_t j=0; j<n; j++) {
      if (i == j)
	eVal=2;
      else {
	if (i == j-1 || i == j+1)
	  eVal=-1;
	else
	  eVal=0;
      }
      T fVal=eVal/eFact;
      RetMat(i,j)=eVal/eFact;
    }
  return RetMat;
}


template<typename T>
void FuncInsertMatrixInList(MyMatrix<T> const& eMat, 
			    std::vector<MyMatrix<T> > & ListMat, 
			    std::vector<int> &ListStatus, 
			    AllTol<T> const& RecTol)
{
  int nbMat=ListStatus.size();
  for (int iMat=0; iMat<nbMat; iMat++) {
    std::cerr << "FuncInsertMatrix iMat=" << iMat << "/" << nbMat << "\n";
    bool test=T_TestGramMatrixEquivalence(eMat, ListMat[iMat], RecTol.IsomTol);
    if (test)
      return;
  }
  ListMat.push_back(eMat);
  ListStatus.push_back(0);
}

template<typename T>
std::vector<int> GetPermutationOnSubspace(std::vector<MyMatrix<T> > const& ListSubspace, MyMatrix<T> const& eMatEquiv)
{
  int nbSub=ListSubspace.size();
  std::vector<int> eList(nbSub);
  for (int iSub=0; iSub<nbSub; iSub++) {
    MyMatrix<T> NewSub=ListSubspace[iSub]*eMatEquiv;
    int pos=PositionSubspace(ListSubspace, NewSub);
    eList[iSub]=pos;
  }
  return eList;
}

template<typename T>
std::vector<MyMatrix<T> > OrbitEnumerationDirMat(int n, int k, 
						 MyMatrix<T> const& eMat, 
						 AllTol<T> const& RecTol)
{
  std::vector<std::vector<int> > ListVectInt;
  TheGroupFormat GRPperm;
  std::vector<MyMatrix<int> > ListMatrGens;
  T_GetGramMatrixAutomorphismGroup(eMat, RecTol.IsomTol, GRPperm, ListMatrGens);
  T MaxDet=RecTol.MaxDetTol;
  std::vector<T> ListDet;
  std::vector<MyMatrix<int> > ListSubspace;
  EnumerationKspaces(eMat, k, MaxDet, ListSubspace, ListDet);
  int nbSub=ListSubspace.size();
  int nbGen=ListMatrGens.size();
  std::vector<permlib::dom_int> v(nbSub);
  std::list<permlib::Permutation::ptr> generatorList;
  for (int iGen=0; iGen<nbGen; iGen++) {
    MyMatrix<int> eMatEquiv=ListMatrGens[iGen];
    std::vector<int> eList=GetPermutationOnSubspace(ListSubspace, eMatEquiv);
    for (int iSub=0; iSub<nbSub; iSub++)
      v[iSub]=eList[iSub];
    generatorList.push_back(permlib::Permutation::ptr(new permlib::Permutation(v)));
  }
  TheGroupFormat GRPsub;
  GRPsub.n=nbSub;
  GRPsub.group=construct(nbSub, generatorList.begin(), generatorList.end());
  std::cerr << "We have GRPsub\n";
  std::vector<MyMatrix<T> > ListEval;
  
  int nSym=n*(n+1)/2;
  MyMatrix<T> ThePerfectCone=MyMatrix<T>(nSym, nbSub);
  for (int iSub=0; iSub<nbSub; iSub++) {
    T TheDet;
    MyMatrix<T> eEvalMat=EvaluatorMatrix(eMat, ListSubspace[iSub], TheDet);
    ListEval.push_back(eEvalMat);
    MyVector<T> TheV=SymmetricMatrixToVector(eEvalMat);
    AssignMatrixRow(ThePerfectCone, iSub, TheV);
  }
  T smallVal=RecTol.tolDirMat;
  MyMatrix<T> TheFAC=cdd::DualDescription(ThePerfectCone, smallVal);
  int nbFac=TheFAC.rows();
  std::cerr << "nbFac=" << nbFac << "\n";
  std::vector<Face> ListTotal(nbFac);
  std::vector<MyMatrix<T> > RawListDirMat;
  for (int iFac=0; iFac<nbFac; iFac++) {
    MyVector<T> eVect=TheFAC.row(iFac);
    MyMatrix<T> eDirMat=VectorToSymmetricMatrixB(eVect, n);
    T eNorm=sqrt(MatrixScalarProduct(eDirMat, eDirMat));
    MyMatrix<T> fDirMat = (1/eNorm) * eDirMat;
    Face eFace(nbSub);
    for (int iSub=0; iSub<nbSub; iSub++) {
      T eScal=MatrixScalarProduct(ListEval[iSub], fDirMat);
      if (eScal < RecTol.tolDirMat)
	eFace[iSub]=1;
    }
    ListTotal[iFac]=eFace;
    RawListDirMat[iFac]=fDirMat;
  }
  int sizTotal=ListTotal.size();
  std::cerr << "sizTotal=" << sizTotal << "\n";
  std::vector<Face> TheReturn=OrbitSplittingSet(ListTotal, GRPsub);
  int nbOrb=TheReturn.size();
  std::cerr << "nbOrb=" << nbOrb << "\n";
  std::vector<MyMatrix<T> > ListDirMat;
  for (int iOrb=0; iOrb<nbOrb; iOrb++) {
    Face eSet=TheReturn[iOrb];
    for (int iFac=0; iFac<nbFac; iFac++)
      if (ListTotal[iFac] == eSet)
	ListDirMat.push_back(RawListDirMat[iFac]);
  }
  std::cerr << "Leaving OrbitEnumerationDirMat\n";
  return ListDirMat;
}


template<typename T>
void PrintKPerfectMatrix(std::ostream & os, int n, int k, 
			 MyMatrix<T> const& ePerfect,
			 AllTol<T> const& RecTol)
{
  std::vector<T> ListDet;
  std::vector<MyMatrix<int> > ListSubspace;
  TheGroupFormat GRPperm;
  std::vector<MyMatrix<int> > ListMatrGens;
  os << "n=" << n << " k=" << k << "\n";
  WriteMatrix(os, ePerfect);
  T MaxDet=RecTol.MaxDetTol;
  EnumerationKspaces(ePerfect, k, MaxDet, ListSubspace, ListDet);
  int nbSub=ListDet.size();
  for (int iSub=0; iSub<nbSub; iSub++) {
    os << "iSub=" << iSub << "/" << nbSub << " det=" << ListDet[iSub] << "\n";
    WriteMatrix(os, ListSubspace[iSub]);
  }
  T_GetGramMatrixAutomorphismGroup(ePerfect, RecTol.IsomTol, GRPperm, ListMatrGens);
  int nbGen=ListMatrGens.size();
  for (int iGen=0; iGen<nbGen; iGen++) {
    os << "iGen=" << iGen << "/" << nbGen << "\n";
    WriteMatrix(os, ListMatrGens[iGen]);
  }
}

template<typename T>
void PrintFamilyKPerfectMatrices(std::ostream & os,
				 int n, int k, 
				 std::vector<MyMatrix<T> > const& ListPerfect, 
				 AllTol<T> const& RecTol)
{
  int nbPerfect=ListPerfect.size();
  os << "nbPerfect=" << nbPerfect << "\n";
  for (int iPerfect=0; iPerfect<nbPerfect; iPerfect++) {
    os << "iPerfect=" << iPerfect << "/" << nbPerfect << "\n";
    MyMatrix<T> ePerfect=ListPerfect[iPerfect];
    PrintKPerfectMatrix(os, n, k, ePerfect, RecTol);
  }
}



template<typename T>
void PrintKPerfectMatrixGAP(std::ostream & os,
			    int n, int k, 
			    MyMatrix<T> const&  ePerfect,
			    AllTol<T> const& RecTol)
{
  os << "rec(n:=" << n << ", k:=" << k << ",\n";
  T MaxDet=RecTol.MaxDetTol;
  std::vector<MyMatrix<int> > ListSubspace;
  std::vector<T> ListDet;
  EnumerationKspaces(ePerfect, k, MaxDet, ListSubspace, ListDet);
  int nbSub=ListDet.size();
  os << "ListSubspace:=[\n";
  for (int iSub=0; iSub<nbSub; iSub++) {
    WriteMatrixGAP(os, ListSubspace[iSub]);
    if (iSub < nbSub-1)
      os << ",\n";
  }
  os << "], \n";
  TheGroupFormat GRPperm;
  std::vector<MyMatrix<int> > ListMatrGens;
  T_GetGramMatrixAutomorphismGroup(ePerfect, RecTol.IsomTol, GRPperm, ListMatrGens);
  int nbGen=ListMatrGens.size();
  os << "ListMatrGen:=[\n";
  for (int iGen=0; iGen<nbGen; iGen++) {
    WriteMatrixGAP(os, ListMatrGens[iGen]);
    if (iGen < nbGen-1)
      os << ",\n";
  }
  os << "])\n";
}

template<typename T>
void PrintFamilyKPerfectMatricesGAP(std::ostream & os,
				    int n, int k, 
				    std::vector<MyMatrix<T> > const& ListPerfect, 
				    AllTol<T> const& RecTol)
{
  os << "return [";
  int nbPerfect=ListPerfect.size();
  for (int iPerfect=0; iPerfect<nbPerfect; iPerfect++) {
    if (iPerfect > 0)
      os << ",\n";
    MyMatrix<T> ePerfect=ListPerfect[iPerfect];
    PrintKPerfectMatrixGAP(os, n, k, ePerfect, RecTol);
  }
  os << "];";
}



template<typename T>
std::vector<MyMatrix<T> > VoronoiEnumerationRankin(int const& n, int const& k, 
	  AllTol<T> const& RecTol)
{
  std::vector<int> ListStatus;
  std::vector<MyMatrix<T> > ListPerfect;
  //  eMat=GetLatticeAn(n, k);
  //  eMat=GetLatticeAdjA4(n, k);
  MyMatrix<T> eMat=GetLatticePbDim5<T>(n,k);
  std::cerr << "RecTol.IsomTol=" << RecTol.IsomTol << "\n";
  FuncInsertMatrixInList(eMat, ListPerfect, ListStatus, RecTol);
  while(true) {
    int nbMat=ListPerfect.size();
    int IsFinished=1;
    std::cerr << "nbMat=" << nbMat << "\n";
    for (int iMat=0; iMat<nbMat; iMat++)
      if (ListStatus[iMat] == 0) {
	std::cerr << "iMat=" << iMat << "/" << nbMat << "\n";
	ListStatus[iMat]=1;
	IsFinished=0;
	MyMatrix<T> fMat=ListPerfect[iMat];
	std::cerr << "eMat=\n";
	WriteMatrix(std::cerr, fMat);
	std::cerr << "Before OrbitEnumerationDirMat\n";
	std::vector<MyMatrix<T> > ListDirMat=OrbitEnumerationDirMat(n, k,
								    fMat, RecTol);
	int nbDir=ListDirMat.size();
	std::cerr << "nbDir=" << nbDir << "\n";
	for (int iDir=0; iDir<nbDir; iDir++) {
	  std::cerr << "iDir=" << iDir << "/" << nbDir << "\n";
	  MyMatrix<T> DirMat=ListDirMat[iDir];
	  std::cerr << "DirMat=\n";
	  WriteMatrix(std::cerr, DirMat);
	  MyMatrix<T> NewMat=KFlippingOperation(fMat, k, DirMat, RecTol);
	  std::cerr << "NewMat=\n";
	  WriteMatrix(std::cerr, NewMat);
	  std::cerr << "Eigenvalues of PreRetMat\n";
	  PrintEigenvalues(std::cerr, NewMat);
	  std::cerr << "Before inserting into ListPerfect\n";
	  FuncInsertMatrixInList(NewMat, ListPerfect, ListStatus, RecTol);
	  std::cerr << "After inserting into ListPerfect\n";
	}
      }
    if (IsFinished == 1)
      break;
  }
  std::cerr << "Enumeration is finished\n";
  return ListPerfect;
}

template<typename T>
AllTol<T> ReadToleranceFile(std::istream & is)
{
  AllTol<T> RecTol;
  int nbIteration;
  is >> nbIteration;
  RecTol.nbIteration=nbIteration;
  std::cerr << "nbIteration=" << nbIteration << "\n";
  //
  T tolDirMat;
  is >> tolDirMat;
  RecTol.tolDirMat=tolDirMat;
  std::cerr << "  tolDirMat=" << tolDirMat << "\n";
  //
  T deltaInc;
  is >> deltaInc;
  RecTol.deltaInc=deltaInc;
  std::cerr << "  deltaInc=" << deltaInc << "\n";
  //
  T IsomTol;
  is >> IsomTol;
  RecTol.IsomTol=IsomTol;
  std::cerr << "  IsomTol=" << IsomTol << "\n";
  //
  T MaxDetTol;
  is >> MaxDetTol;
  RecTol.MaxDetTol=MaxDetTol;
  std::cerr << "  MaxDetTol=" << MaxDetTol << "\n";
  //
  T MaxDetTargetTol;
  is >> MaxDetTargetTol;
  RecTol.MaxDetTargetTol=MaxDetTargetTol;
  std::cerr << "  MaxDetTargetTol=" << MaxDetTargetTol << "\n";
  //
  T MaxMoveError;
  is >> MaxMoveError;
  RecTol.MaxMoveError=MaxMoveError;
  std::cerr << "  MaxMoveError=" << MaxMoveError << "\n";
  //
  return RecTol;
}
#endif

#ifndef INCLUDE_CVP_NIEMEIER_ALGORITHM
#define INCLUDE_CVP_NIEMEIER_ALGORITHM

#include "LatticeDefinitions.h"
#include "NiemeierN23_24A1.h"
#include "NumberTheory.h"
#include "COMB_Combinatorics.h"


//
// We have eBasis*GramMat*TransposedMat(eBasis) = q0 
// with q0 the diagonal matrix having 2 on the diagonal.
// Therefore, we have GramMat = Inv(eBasis)*Q0*TransposedMat(Inv(eBasis))
// So, we have the equation GramMat[eV - x] to solve with
// x= eCos + y*eBasis
// Therefore, the equation becomes
// q0[ (eV - eCos)*Inv(eBasis) - y] to solve.
// Writing fV=eV*Inv(eBasis) and eCosInv=eCos*Inv(eBasis)
// we are led to
// q0[ fV - eCosInv - y]
//
template<typename T>
resultCVP<T> CVP_N23_24A1_Version1(MyVector<T> const& eV)
{
  
  MyMatrix<T> eBasis=Get_N23_eBasis<T>();
  MyMatrix<T> InvBasis=Get_N23_InvBasis<T>();
  MyMatrix<T> TotBasisQuotInv=Get_N23_TotBasisQuotInv<T>();
  BlockIteration BlCoset(12, 2);
  MyVector<T> fV=ProductVectorMatrix(eV, InvBasis);
  std::vector<MyVector<T>> ListNearest;
  T MinNorm=0;
  bool IsFirst=true;
  T eOne=1;
  T eTwo=2;
  T eHalf=eOne/eTwo;
  MyVector<T> eTransCoset(24);
  MyVector<T> ListLow(24);
  while(true) {
    //    std::cerr << "First loop\n";
    std::vector<int> LVal=BlCoset.GetVect();
    /*    for (int iDim=0; iDim<12; iDim++)
      std::cerr << " " << LVal[iDim];
      std::cerr << "\n";*/
    for (int i=0; i<24; i++)
      eTransCoset(i)=0;
    for (int iDim=0; iDim<12; iDim++)
      if (LVal[iDim] == 1) {
	for (int i=0; i<24; i++)
	  eTransCoset(i) += TotBasisQuotInv(iDim,i);
      }
    MyVector<T> eSumVect=fV - eTransCoset;
    std::vector<int> ListDimPair;
    T eNorm=0;
    for (int i=0; i<24; i++) {
      T eVal=eSumVect(i);
      T eNear=UniversalNearestInteger<T,T>(eVal);
      ListLow(i)=eNear;
      T eDiff=eVal - eNear;
      //      std::cerr << "i=" << i << " eDiff=" << eDiff << "\n";
      eNorm += eDiff*eDiff;
      if (IsInteger(eDiff - eHalf)) 
	ListDimPair.push_back(i);
    }
    bool WeAppend=false;
    if (IsFirst) {
      WeAppend=true;
      IsFirst=false;
      MinNorm=eNorm;
    }
    else {
      if (eNorm == MinNorm) {
	WeAppend=true;
      }
      else {
	if (eNorm < MinNorm) {
	  WeAppend=true;
	  MinNorm=eNorm;
	  ListNearest.clear();
	}
      }
    }
    //    std::cerr << "eNorm=" << eNorm << " IsFirst=" << IsFirst << " MinNorm=" << MinNorm << "\n";
    if (WeAppend) {
      int nbPair=ListDimPair.size();
      //      std::cerr << "nbPair=" << nbPair << " MinNorm=" << eNorm << "\n";
      BlockIteration BlCVP(nbPair, 2);
      while(true) {
	std::vector<int> Wvect=BlCVP.GetVect();
	MyVector<T> W=ListLow;
	for (int iPair=0; iPair<nbPair; iPair++)
	  if (Wvect[iPair] == 1) {
	    int eIdx=ListDimPair[iPair];
	    W[eIdx]++;
	  }
	MyVector<T> fVnear= W + eTransCoset;
	ListNearest.push_back(fVnear);
	int val=BlCVP.IncrementShow();
	if (val == -1)
	  break;
      }
    }
    int val=BlCoset.IncrementShow();
    if (val == -1)
      break;
  }
  int nbFinal=ListNearest.size();
  MyMatrix<T> MatNearest(nbFinal,24);
  for (int iFin=0; iFin<nbFinal; iFin++) {
    MyVector<T> fVnear=ListNearest[iFin];
    MyVector<T> eVnear=ProductVectorMatrix(fVnear, eBasis);
    for (int i=0; i<24; i++)
      MatNearest(iFin,i)=eVnear(i);
  }
  T NearestNorm=2*MinNorm;
  //  std::cerr << "End of CVP_N23_24A1\n";
  return {NearestNorm, MatNearest};
}


// This version has a more ordered iteration
// over cosets. This turned to give no gain. Now 9m32s vs 9m30 for Version1.
template<typename T>
resultCVP<T> CVP_N23_24A1_Version2(MyVector<T> const& eV)
{
  
  MyMatrix<T> eBasis=Get_N23_eBasis<T>();
  MyMatrix<T> InvBasis=Get_N23_InvBasis<T>();
  MyMatrix<T> TotBasisQuotInv=Get_N23_TotBasisQuotInv<T>();
  BlockIteration BlCoset(12, 2);
  MyVector<T> fV=ProductVectorMatrix(eV, InvBasis);
  T MinNorm=0;
  bool IsFirst=true;
  T eOne=1;
  T eTwo=2;
  T eHalf=eOne/eTwo;
  MyVector<T> eTransCoset(24);
  int nbFinal=0;
  std::vector<MyVector<T>> ListTransCoset;
  while(true) {
    std::vector<int> LVal=BlCoset.GetVect();
    for (int i=0; i<24; i++)
      eTransCoset(i)=0;
    //    std::cerr << "Before Assignation of eTransCoset\n";
    for (int iDim=0; iDim<12; iDim++)
      if (LVal[iDim] == 1) {
	for (int i=0; i<24; i++)
	  eTransCoset(i) += TotBasisQuotInv(iDim,i);
      }
    //    std::cerr << "After Assignation of eTransCoset\n";
    MyVector<T> eSumVect=fV - eTransCoset;
    int nbPair=0;
    T eNorm=0;
    for (int i=0; i<24; i++) {
      T eVal=eSumVect(i);
      T eNear=UniversalNearestInteger<T,T>(eVal);
      T eDiff=eVal - eNear;
      //      std::cerr << "i=" << i << " eDiff=" << eDiff << "\n";
      eNorm += eDiff*eDiff;
      if (IsInteger(eDiff - eHalf))
	nbPair++;
    }
    bool WeAppend=false;
    if (IsFirst) {
      WeAppend=true;
      IsFirst=false;
      MinNorm=eNorm;
    }
    else {
      if (eNorm == MinNorm) {
	WeAppend=true;
      }
      else {
	if (eNorm < MinNorm) {
	  WeAppend=true;
	  MinNorm=eNorm;
	  nbFinal=0;
	  ListTransCoset.clear();
	}
      }
    }
    //    std::cerr << "eNorm=" << eNorm << " IsFirst=" << IsFirst << " MinNorm=" << MinNorm << "\n";
    if (WeAppend) {
      int eSiz=1;
      for (int iPair=0; iPair<nbPair; iPair++)
	eSiz *=2;
      nbFinal += eSiz;
      ListTransCoset.push_back(eTransCoset);
    }
    int val=BlCoset.IncrementShow();
    if (val == -1)
      break;
  }
  //  std::cerr << "nbFinal=" << nbFinal << "\n";
  int iFin=0;
  MyMatrix<T> MatNearest(nbFinal,24);
  MyVector<T> ListLow(24);
  for (auto & fTransCoset : ListTransCoset) {
    //    std::cerr << "Iterate, step 1\n";
    MyVector<T> eSumVect=fV - fTransCoset;
    std::vector<int> ListDimPair;
    //    T eNorm=0;
    //    std::cerr << "Iterate, step 2\n";
    for (int i=0; i<24; i++) {
      T eVal=eSumVect(i);
      T eNear=UniversalNearestInteger<T,T>(eVal);
      ListLow(i)=eNear;
      T eDiff=eVal - eNear;
      //      std::cerr << "i=" << i << " eDiff=" << eDiff << "\n";
      //      eNorm += eDiff*eDiff;
      if (IsInteger(eDiff - eHalf)) 
	ListDimPair.push_back(i);
    }
    //    std::cerr << "eNorm=" << eNorm << "\n";
    int nbPair=ListDimPair.size();
    //    std::cerr << "nbPair=" << nbPair << "\n";
    BlockIteration BlCVP(nbPair, 2);
    //    std::cerr << "Iterate, step 3\n";
    while(true) {
      std::vector<int> Wvect=BlCVP.GetVect();
      //      std::cerr << "Loop, step 1\n";
      MyVector<T> W=ListLow;
      for (int iPair=0; iPair<nbPair; iPair++)
	if (Wvect[iPair] == 1) {
	  int eIdx=ListDimPair[iPair];
	  W[eIdx]++;
	}
      MyVector<T> fVnear= W + fTransCoset;
      MyVector<T> eVnear=ProductVectorMatrix(fVnear, eBasis);
      for (int i=0; i<24; i++)
	MatNearest(iFin,i)=eVnear(i);
      iFin++;
      int val=BlCVP.IncrementShow();
      if (val == -1)
	break;
      //      std::cerr << "Loop, step 7\n";
    }
    //    std::cerr << "Iterate, step 4\n";
  }
  T NearestNorm=2*MinNorm;
  //  std::cerr << "End of CVP_N23_24A1\n";
  return {NearestNorm, MatNearest};
}



template<typename T>
MyMatrix<T> GetShiftedMatrix(MyMatrix<T> const& eMat, int const& eSize)
{
  int nbRow=eMat.rows();
  int nbCol=eMat.cols();
  MyMatrix<T> RetMat(nbRow, nbCol);
  for (int iRow=0; iRow<nbRow; iRow++) {
    RetMat.row(iRow)=eMat.row(iRow);
    T alpha=-(eSize-1);
    for (int jRow=0; jRow<iRow; jRow++) {
      RetMat.row(iRow)=RetMat.row(iRow) + alpha*RetMat.row(jRow);
    }
  }
  return RetMat;
}


// This version shows some improvement with respect
// to the iteration over he cosets. Definitely less
// operations are needed.
// Time is 9m18s vs 9m32s before. Slight improvement.
template<typename T>
resultCVP<T> CVP_N23_24A1_Version3(MyVector<T> const& eV)
{
  MyMatrix<T> eBasis=Get_N23_eBasis<T>();
  MyMatrix<T> InvBasis=Get_N23_InvBasis<T>();
  MyMatrix<T> TotBasisQuotInv=Get_N23_TotBasisQuotInv<T>();
  MyMatrix<T> TotBasisQuotInvShift=GetShiftedMatrix(TotBasisQuotInv, 2);
  BlockIteration BlCoset(12, 2);
  MyVector<T> fV=ProductVectorMatrix(eV, InvBasis);
  T MinNorm=0;
  bool IsFirst=true;
  T eOne=1;
  T eTwo=2;
  T eHalf=eOne/eTwo;
  MyVector<T> eTransCoset(24);
  for (int i=0; i<24; i++)
    eTransCoset(i)=0;
  int nbFinal=0;
  std::vector<MyVector<T>> ListTransCoset;
  while(true) {
    std::vector<int> LVal=BlCoset.GetVect();
    MyVector<T> eSumVect=fV - eTransCoset;
    int nbPair=0;
    T eNorm=0;
    for (int i=0; i<24; i++) {
      T eVal=eSumVect(i);
      T eNear=UniversalNearestInteger<T,T>(eVal);
      T eDiff=eVal - eNear;
      //      std::cerr << "i=" << i << " eDiff=" << eDiff << "\n";
      eNorm += eDiff*eDiff;
      if (IsInteger(eDiff - eHalf))
	nbPair++;
    }
    bool WeAppend=false;
    if (IsFirst) {
      WeAppend=true;
      IsFirst=false;
      MinNorm=eNorm;
    }
    else {
      if (eNorm == MinNorm) {
	WeAppend=true;
      }
      else {
	if (eNorm < MinNorm) {
	  WeAppend=true;
	  MinNorm=eNorm;
	  nbFinal=0;
	  ListTransCoset.clear();
	}
      }
    }
    //    std::cerr << "eNorm=" << eNorm << " IsFirst=" << IsFirst << " MinNorm=" << MinNorm << "\n";
    if (WeAppend) {
      int eSiz=1;
      for (int iPair=0; iPair<nbPair; iPair++)
	eSiz *=2;
      nbFinal += eSiz;
      ListTransCoset.push_back(eTransCoset);
    }
    int val=BlCoset.IncrementShow();
    if (val == -1)
      break;
    for (int i=0; i<24; i++)
      eTransCoset(i) += TotBasisQuotInvShift(val,i);
  }
  //  std::cerr << "nbFinal=" << nbFinal << "\n";
  int iFin=0;
  MyMatrix<T> MatNearest(nbFinal,24);
  MyVector<T> ListLow(24);
  for (auto & fTransCoset : ListTransCoset) {
    //    std::cerr << "Iterate, step 1\n";
    MyVector<T> eSumVect=fV - fTransCoset;
    std::vector<int> ListDimPair;
    //    T eNorm=0;
    //    std::cerr << "Iterate, step 2\n";
    for (int i=0; i<24; i++) {
      T eVal=eSumVect(i);
      T eNear=UniversalNearestInteger<T,T>(eVal);
      ListLow(i)=eNear;
      T eDiff=eVal - eNear;
      //      std::cerr << "i=" << i << " eDiff=" << eDiff << "\n";
      //      eNorm += eDiff*eDiff;
      if (IsInteger(eDiff - eHalf)) 
	ListDimPair.push_back(i);
    }
    //    std::cerr << "eNorm=" << eNorm << "\n";
    int nbPair=ListDimPair.size();
    //    std::cerr << "nbPair=" << nbPair << "\n";
    BlockIteration BlCVP(nbPair, 2);
    //    std::cerr << "Iterate, step 3\n";
    while(true) {
      std::vector<int> Wvect=BlCVP.GetVect();
      //      std::cerr << "Loop, step 1\n";
      MyVector<T> W=ListLow;
      for (int iPair=0; iPair<nbPair; iPair++)
	if (Wvect[iPair] == 1) {
	  int eIdx=ListDimPair[iPair];
	  W[eIdx]++;
	}
      MyVector<T> fVnear= W + fTransCoset;
      MyVector<T> eVnear=ProductVectorMatrix(fVnear, eBasis);
      for (int i=0; i<24; i++)
	MatNearest(iFin,i)=eVnear(i);
      iFin++;
      int val=BlCVP.IncrementShow();
      if (val == -1)
	break;
      //      std::cerr << "Loop, step 7\n";
    }
    //    std::cerr << "Iterate, step 4\n";
  }
  T NearestNorm=2*MinNorm;
  //  std::cerr << "End of CVP_N23_24A1\n";
  return {NearestNorm, MatNearest};
}


// We change the rounding function.
// hopefully faster.
template<typename T>
resultCVP<T> CVP_N23_24A1(MyVector<T> const& eV)
{
  MyMatrix<T> eBasis=Get_N23_eBasis<T>();
  MyMatrix<T> InvBasis=Get_N23_InvBasis<T>();
  MyMatrix<T> TotBasisQuotInv=Get_N23_TotBasisQuotInv<T>();
  MyMatrix<T> TotBasisQuotInvShift=GetShiftedMatrix(TotBasisQuotInv, 2);
  BlockIteration BlCoset(12, 2);
  MyVector<T> fV=ProductVectorMatrix(eV, InvBasis);
  T MinNorm=0;
  bool IsFirst=true;
  T eOne=1;
  T eTwo=2;
  T eHalf=eOne/eTwo;
  MyVector<T> eTransCoset(24);
  for (int i=0; i<24; i++)
    eTransCoset(i)=0;
  int nbFinal=0;
  std::vector<MyVector<T>> ListTransCoset;
  while(true) {
    std::vector<int> LVal=BlCoset.GetVect();
    MyVector<T> eSumVect=fV - eTransCoset;
    int nbPair=0;
    T eNorm=0;
    for (int i=0; i<24; i++) {
      T eVal=eSumVect(i);
      T eNear=NearestInteger_rpi(eVal);
      T eDiff=eVal - eNear;
      //      std::cerr << "i=" << i << " eDiff=" << eDiff << "\n";
      eNorm += eDiff*eDiff;
      if (IsInteger(eDiff - eHalf))
	nbPair++;
    }
    bool WeAppend=false;
    if (IsFirst) {
      WeAppend=true;
      IsFirst=false;
      MinNorm=eNorm;
    }
    else {
      if (eNorm == MinNorm) {
	WeAppend=true;
      }
      else {
	if (eNorm < MinNorm) {
	  WeAppend=true;
	  MinNorm=eNorm;
	  nbFinal=0;
	  ListTransCoset.clear();
	}
      }
    }
    //    std::cerr << "eNorm=" << eNorm << " IsFirst=" << IsFirst << " MinNorm=" << MinNorm << "\n";
    if (WeAppend) {
      int eSiz=1;
      for (int iPair=0; iPair<nbPair; iPair++)
	eSiz *=2;
      nbFinal += eSiz;
      ListTransCoset.push_back(eTransCoset);
    }
    int val=BlCoset.IncrementShow();
    if (val == -1)
      break;
    for (int i=0; i<24; i++)
      eTransCoset(i) += TotBasisQuotInvShift(val,i);
  }
  //  std::cerr << "nbFinal=" << nbFinal << "\n";
  int iFin=0;
  MyMatrix<T> MatNearest(nbFinal,24);
  MyVector<T> ListUpp(24);
  for (auto & fTransCoset : ListTransCoset) {
    //    std::cerr << "Iterate, step 1\n";
    MyVector<T> eSumVect=fV - fTransCoset;
    std::vector<int> ListDimPair;
    //    T eNorm=0;
    //    std::cerr << "Iterate, step 2\n";
    for (int i=0; i<24; i++) {
      T eVal=eSumVect(i);
      T eNear=NearestInteger_rpi(eVal);
      ListUpp(i)=eNear;
      T eDiff=eVal - eNear;
      //      std::cerr << "i=" << i << " eDiff=" << eDiff << "\n";
      //      eNorm += eDiff*eDiff;
      if (IsInteger(eDiff - eHalf)) 
	ListDimPair.push_back(i);
    }
    //    std::cerr << "eNorm=" << eNorm << "\n";
    int nbPair=ListDimPair.size();
    //    std::cerr << "nbPair=" << nbPair << "\n";
    BlockIteration BlCVP(nbPair, 2);
    //    std::cerr << "Iterate, step 3\n";
    while(true) {
      std::vector<int> Wvect=BlCVP.GetVect();
      //      std::cerr << "Loop, step 1\n";
      MyVector<T> W=ListUpp;
      for (int iPair=0; iPair<nbPair; iPair++)
	if (Wvect[iPair] == 1) {
	  int eIdx=ListDimPair[iPair];
	  W[eIdx]--;
	}
      MyVector<T> fVnear= W + fTransCoset;
      MyVector<T> eVnear=ProductVectorMatrix(fVnear, eBasis);
      for (int i=0; i<24; i++)
	MatNearest(iFin,i)=eVnear(i);
      iFin++;
      int val=BlCVP.IncrementShow();
      if (val == -1)
	break;
      //      std::cerr << "Loop, step 7\n";
    }
    //    std::cerr << "Iterate, step 4\n";
  }
  T NearestNorm=2*MinNorm;
  //  std::cerr << "End of CVP_N23_24A1\n";
  return {NearestNorm, MatNearest};
}



#endif

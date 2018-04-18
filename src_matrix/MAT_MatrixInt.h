#ifndef TEMP_MATRIX_INTEGRAL_RING
#define TEMP_MATRIX_INTEGRAL_RING

#include "MAT_Matrix.h"
#include "NumberTheory.h"

template<typename T>
MyMatrix<int> TMat_ConvertToInt(MyMatrix<T> const&eMatI)
{
  double eVal_d;
  int nbRow=eMatI.rows();
  int nbCol=eMatI.cols();
  MyMatrix<int> eMatO(nbRow, nbCol);
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      T eVal=eMatI(iRow, iCol);
      GET_DOUBLE(eVal, eVal_d);
      int eVal_i=int(round(eVal_d));
      eMatO(iRow, iCol)=eVal_i;
    }
  return eMatO;
}





// Now declarations of generic code.
// The code below generally requires the field T to be the ring (or fraction ring) of
// a factorial ring.
// Operations may work for fields and rings as well.
template<typename T>
bool IsIntegralVector(MyVector<T> const& V)
{
  int nbCol=V.size();
  for (int i=0; i<nbCol; i++)
    if (!IsInteger(V(i)))
      return false;
  return true;
}


template<typename T>
bool IsIntegralMatrix(MyMatrix<T> const& M)
{
  int nbRow=M.rows();
  int nbCol=M.cols();
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++)
      if (!IsInteger(M(iRow, iCol)))
	return false;
  return true;
}


// T must be an integral domain with a norm function
// for computing GCD.
// eMat is a m x n matrix which defines a submodule L of T^n.
// The function computes the index i of L in T^n.
template<typename T>
T Int_IndexLattice(MyMatrix<T> const& eMat)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain");
  int iRowF=-1, iColF=-1;
  MyMatrix<T> eMatW=eMat;
  int nbCol=eMat.cols();
  int nbRow=eMat.rows();
  std::vector<int> colStat(nbCol,1);
  std::vector<int> rowStat(nbRow,1);
  int nbDone=0;
  T TheIndex=1;
  while(true) {
    bool IsFirst=true;
    int MinPivot=0;
    for (int iCol=0; iCol<nbCol; iCol++)
      if (colStat[iCol] == 1)
	for (int iRow=0; iRow<nbRow; iRow++)
	  if (rowStat[iRow] == 1) {
	    T eVal=eMatW(iRow, iCol);
	    if (eVal != 0) {
	      int eValA=T_Norm(eVal);
	      if (IsFirst) {
		iRowF=iRow;
		iColF=iCol;
		MinPivot=eValA;
	      }
	      else {
		if (eValA < MinPivot) {
		  iRowF=iRow;
		  iColF=iCol;
		  MinPivot=eValA;
		}
	      }
	      IsFirst=false;
	    }
	  }
    if (IsFirst)
      return 0;
    std::cerr << "MinPivot=" << MinPivot << " iRowF=" << iRowF << " iColF=" << iColF << "\n";
    if (MinPivot == 0) {
      std::cerr << "Clear error in the code of IndexLattice\n";
      throw TerminalException{1};
    }
    std::cerr << "Before row operations\n";
    T ThePivot=eMatW(iRowF, iColF);
    bool IsFinished=true;
    for (int iRow=0; iRow<nbRow; iRow++) {
      std::cerr << "iRow=" << iRow << " nbRow=" << nbRow << "\n";
      if (rowStat[iRow] == 1 && iRow != iRowF) {
	T eVal=eMatW(iRow, iColF);
	std::cerr << "eVal=" << eVal << "\n";
	if (eVal != 0) {
	  IsFinished=false;
	  T TheQ=QuoInt(eVal, ThePivot);
	  std::cerr << "eVal=" << eVal << " ThePivot=" << ThePivot << " TheQ=" << TheQ << "\n";
	  eMatW.row(iRow) -= TheQ*eMatW.row(iRowF);
	}
      }
    }
    std::cerr << "After row operations IsFinished=" << IsFinished << "\n";
    if (IsFinished) {
      colStat[iColF]=0;
      rowStat[iRowF]=0;
      nbDone++;
      TheIndex=TheIndex*ThePivot;
    }
    std::cerr << "Now updated index\n";
    if (nbDone == nbCol)
      return TheIndex;
    std::cerr << "Continuing loop\n";
  }
}

template<typename T>
struct GCD_int {
  T gcd;
  std::vector<T> ListA;
  std::vector<std::vector<T>> ListListK;
  MyMatrix<T> Pmat;
};


template<typename T>
void WriteGCD_int(std::ostream & os, GCD_int<T> const& eGCD)
{
  int siz=eGCD.ListListK.size();
  auto ThePrint=[&](std::vector<T> const& eV) -> void {
    os << "[";
    int len=eV.size();
    for (int i=0; i<len; i++) {
      if (i>0)
	os << ",";
      os << eV[i];
    }
    os << "]";
  };
  os << "gcd=" << eGCD.gcd << "\n";
  os << "A=";
  ThePrint(eGCD.ListA);
  os << "\n";
  for (int iVect=0; iVect<siz; iVect++) {
    os << "iVect=" << iVect << " V=";
    ThePrint(eGCD.ListListK[iVect]);
    os << "\n";
  }
}


template<typename T>
GCD_int<T> ComputePairGcd(T const& m, T const& n)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain");
  T f, g, h, fm, gm, hm, q;
  if (n == 0 && m == 0) {
    f=0;
    std::vector<T> eVect{1, 0};
    std::vector<T> fVect{0, 1};
    std::vector<std::vector<T> > ListListK = {fVect};
    MyMatrix<T> Pmat=IdentityMat<T>(2);
    GCD_int<T> Case2{f, eVect, ListListK, Pmat};
    return Case2;
  }
  if (0 <= m) {
    f=m; fm=1;
  }
  else {
    f=-m; fm=-1;
  }
  if (0 <= n) {
    g=n; gm=0;
  }
  else {
    g=-n; gm=0;
  }
  while (g != 0) {
    q = QuoInt( f, g );
    h = g;          hm = gm;
    g = f - q * g;  gm = fm - q * gm;
    f = h;          fm = hm;
  }
  T eCoeff1, eCoeff2;
  if (n == 0) {
    eCoeff1=fm;
    eCoeff2=0;
  }
  else {
    eCoeff1=fm;
    eCoeff2=(f - fm * m) / n;
  }
  //  std::cerr << "n=" << n << " m=" << m << "\n";
  std::vector<T> eVect{eCoeff1, eCoeff2};
  std::vector<T> fVect{-n/f, m/f};
  std::vector<std::vector<T> > ListListK = {fVect};
  MyMatrix<T> Pmat(2,2);
  Pmat(0,0)=eVect[0];
  Pmat(1,0)=eVect[1];
  Pmat(0,1)=fVect[0];
  Pmat(1,1)=fVect[1];
  GCD_int<T> Case2{f, eVect, ListListK, Pmat};
  return Case2;
}

template<typename T>
T NakedGcdPair(T const& a, T const& b)
{
  GCD_int<T> eGCD=ComputePairGcd(a, b);
  return eGCD.gcd;
}

template<typename T>
inline typename std::enable_if<is_totally_ordered<T>::value,T>::type GcdPair(T const& a, T const& b)
{
  T eGCD=NakedGcdPair(a,b);
  if (eGCD < 0) {
    return -eGCD;
  }
  return eGCD;
}

template<typename T>
inline typename std::enable_if<(not is_totally_ordered<T>::value),T>::type GcdPair(T const& a, T const& b)
{
  T eGCD=NakedGcdPair(a,b);
  return eGCD;
}


template<typename T>
T LCMpair(T const& a, T const& b)
{
  if (a == 0)
    return b;
  if (b == 0)
    return a;
  return a*b/GcdPair(a,b);
}



template<typename T>
T ComputeLCM(std::vector<T> const& eVect)
{
  int siz=eVect.size();
  T eLCM=1;
  for (int i=0; i<siz; i++) {
    //    std::cerr << "i=" << i << " eLCM=" << eLCM << " eVect[i]=" << eVect[i] << "\n";
    eLCM=LCMpair(eLCM, eVect[i]);
    //    std::cerr << "After LCMpair operation eLCM=" << eLCM << "\n";
  }
  return eLCM;
}

template<typename T>
struct FractionMatrix {
  T TheMult;
  MyMatrix<T> TheMat;
};


template<typename T>
FractionMatrix<T> RemoveFractionMatrixPlusCoeff(MyMatrix<T> const& M)
{
  int nbRow=M.rows();
  int nbCol=M.cols();
  int idx=0;
  std::vector<T> eVect(nbRow*nbCol);
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      T eDen=GetDenominator(M(iRow,iCol));
      eVect[idx]=eDen;
      idx++;
    }
  T eLCM=ComputeLCM(eVect);
  MyMatrix<T> Mret(nbRow, nbCol);
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      T eVal=M(iRow,iCol)*eLCM;
      Mret(iRow, iCol)=eVal;
    }
  return {eLCM, Mret};
}

template<typename T>
MyMatrix<T> RemoveFractionMatrix(MyMatrix<T> const& M)
{
  return RemoveFractionMatrixPlusCoeff(M).TheMat;
}

template<typename T>
struct FractionVector {
  T TheMult;
  MyVector<T> TheVect;
};


template<typename T>
FractionVector<T> RemoveFractionVectorPlusCoeff(MyVector<T> const& V)
{
  //  std::cerr << "RemoveFractionVectorPlusCoeff, step 1\n";
  int n=V.size();
  std::vector<T> eVect(n);
  //  std::cerr << "RemoveFractionVectorPlusCoeff, step 2\n";
  for (int i=0; i<n; i++) {
    T eDen=GetDenominator(V(i));
    eVect[i]=eDen;
  }
  //  std::cerr << "RemoveFractionVectorPlusCoeff, step 3\n";
  T eLCM=ComputeLCM(eVect);
  //  std::cerr << "RemoveFractionVectorPlusCoeff, step 4\n";
  MyVector<T> Vret(n);
  for (int i=0; i<n; i++) {
    T eVal=V(i)*eLCM;
    Vret(i)=eVal;
  }
  //  std::cerr << "RemoveFractionVectorPlusCoeff, step 5\n";
  return {eLCM, Vret};
}

template<typename T>
inline typename std::enable_if<is_implementation_of_Z<T>::value,MyVector<T>>::type CanonicalizeVector(MyVector<T> const& V)
{
  //  std::cerr << "MAT_MatrixInt : CanonicalizeVector\n";
  return RemoveFractionVectorPlusCoeff(V).TheVect;
}


template<typename T>
inline typename std::enable_if<(not is_float_arithmetic<T>::value),MyVector<T>>::type RemoveFractionVector(MyVector<T> const& V)
{
  return RemoveFractionVectorPlusCoeff(V).TheVect;
}












int IsVectorPrimitive(MyVector<int> const& TheV)
{
  int n=TheV.size();
  int TheGCD=TheV(0);
  for (int i=1; i<n; i++) {
    int eValI=TheV(i);
    GCD_int<int> eRec=ComputePairGcd(TheGCD, eValI);
    TheGCD=eRec.gcd;
  }
  if (abs(TheGCD) == 1)
    return 1;
  return 0;
}








template<typename T>
GCD_int<T> ComputeGCD_information(std::vector<T> const& ListX)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain");
  int siz=ListX.size();
  if (siz == 1) {
    T gcd=ListX[0];
    std::vector<T> ListA={1};
    MyMatrix<T> Pmat=IdentityMat<T>(1);
    GCD_int<T> eGCD_int{gcd, ListA, {}, Pmat};
    return eGCD_int;
  }
  if (siz == 2)
    return ComputePairGcd(ListX[0], ListX[1]);
  std::vector<T> ListXred(ListX.begin(), ListX.end()-1);
  GCD_int<T> eGCD_int=ComputeGCD_information(ListXred);
  GCD_int<T> eGCD2=ComputePairGcd(eGCD_int.gcd, ListX[siz-1]);
  std::vector<T> NewListA;
  for (int i=0; i<siz-1; i++) {
    T eVal=eGCD2.ListA[0]*eGCD_int.ListA[i];
    NewListA.push_back(eVal);
  }
  NewListA.push_back(eGCD2.ListA[1]);
  //
  std::vector<std::vector<T> > NewListListK(siz-1);
  for (int i=0; i<siz-2; i++) {
    std::vector<T> eVect=eGCD_int.ListListK[i];
    eVect.push_back(0);
    NewListListK[i]=eVect;
  }
  std::vector<T> fVect(siz);
  for (int i=0; i<siz-1; i++) {
    T eVal=eGCD2.ListListK[0][0]*eGCD_int.ListA[i];
    fVect[i]=eVal;
  }
  fVect[siz-1]=eGCD2.ListListK[0][1];
  NewListListK[siz-2]=fVect;
  MyMatrix<T> Pmat=MyMatrix<T>(siz,siz);
  for (int i=0; i<siz; i++)
    Pmat(i,0)=NewListA[i];
  for (int iCol=1; iCol<siz; iCol++)
    for (int iRow=0; iRow<siz; iRow++)
      Pmat(iRow, iCol)=NewListListK[iCol-1][iRow];
  //
  GCD_int<T> retGCD{eGCD2.gcd, NewListA, NewListListK, Pmat};
  return retGCD;
}


template<typename T>
void SwitchRow(MyMatrix<T> & eMat, int const& iRow, int const& jRow)
{
  int nbCol=eMat.cols();
  if (iRow == jRow)
    return;
  for (int iCol=0; iCol<nbCol; iCol++) {
    T eVal1=eMat(iRow, iCol);
    T eVal2=eMat(jRow, iCol);
    eMat(iRow, iCol)=eVal2;
    eMat(jRow, iCol)=eVal1;
  }
}

template<typename T>
void INT_ClearColumn(MyMatrix<T> & eMat, int const& iCol, int const& MinAllowedRow, int & iRowFound)
{
  using Treal=typename underlying_totally_ordered_ring<T>::real_type;
  int nbRow=eMat.rows();
  while(true) {
    Treal MinVal=-1;
    int nbFound=0;
    for (int iRow=MinAllowedRow; iRow<nbRow; iRow++) {
      T eVal=eMat(iRow, iCol);
      if (eVal != 0) {
	Treal AbsEVal=T_NormGen(eVal);
	if (nbFound == 0) {
	  MinVal=AbsEVal;
	  iRowFound=iRow;
	}
	else {
	  if (AbsEVal < MinVal) {
	    MinVal=AbsEVal;
	    iRowFound=iRow;
	  }
	}
	nbFound++;
      }
    }
    if (nbFound == 0) {
      std::cerr << "The column is zero. No work possible\n";
      throw TerminalException{1};
    }
    T ThePivot=eMat(iRowFound, iCol);
    for (int iRow=0; iRow<nbRow; iRow++)
      if (iRow != iRowFound) {
	T eVal=eMat(iRow, iCol);
	T TheQ=QuoInt(eVal, ThePivot);
	eMat.row(iRow) -= TheQ*eMat.row(iRowFound);
      }
    if (nbFound == 1)
      return;
  }
}


template<typename T>
bool IsColumnNonEmpty(MyMatrix<T> const& eMat, int const& minAllowed, int const& iCol)
{
  int nbRow=eMat.rows();
  for (int iRow=minAllowed; iRow<nbRow; iRow++) {
    T eVal=eMat(iRow, iCol);
    if (eVal != 0)
      return true;
  }
  return false;
}

// T should be a ring of integers with factorization,
// e.g. Gaussian Integers, Eisenstein Integers, etc.
// but we need to provide QuoInt and T_Norm functions.
//
// The matrix in output is the matrix NullspaceIntMat(TransposedMat(M))
template<typename T>
MyMatrix<T> NullspaceIntTrMat(MyMatrix<T> const& eMat)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain");
  MyMatrix<T> eMatW=eMat;
  int nbCol=eMat.cols();
  int nbRow=eMat.rows();
  std::vector<int> ListIndex;
  std::vector<int> ListNonIndex;
  int eRank=0;
  for (int iCol=0; iCol<nbCol; iCol++)
    if (IsColumnNonEmpty(eMatW, eRank, iCol)) {
      ListIndex.push_back(iCol);
      int iRowFound=444;
      INT_ClearColumn(eMatW, iCol, eRank, iRowFound);
      SwitchRow(eMatW, eRank, iRowFound);
      eRank++;
    }
    else {
      ListNonIndex.push_back(iCol);
    }
  int dimSpace=ListNonIndex.size();
  std::vector<std::vector<T> > TheBasis;
  for (int i=0; i<dimSpace; i++) {
    std::vector<T> eVect;
    for (int j=0; j<dimSpace; j++) {
      if (i == j) {
	eVect.push_back(1);
      }
      else {
	eVect.push_back(0);
      }
    }
    TheBasis.push_back(eVect);
  }
  for (int iRank=0; iRank<eRank; iRank++) {
    int iRow=eRank-1-iRank;
    int iCol=ListIndex[iRow];
    std::vector<T> ListX;
    T eVal=eMatW(iRow, iCol);
    ListX.push_back(eVal);
    std::vector<int> ListRelIndex;
    for (int jRow=iRow+1; jRow<eRank; jRow++) {
      int jCol=ListIndex[jRow];
      ListRelIndex.push_back(jCol);
    }
    for (auto & jCol : ListNonIndex) {
      ListRelIndex.push_back(jCol);
    }
    int sizRelIndex=ListRelIndex.size();
    for (int iVect=0; iVect<dimSpace; iVect++) {
      std::vector<T> eVect=TheBasis[iVect];
      T eSum=0;
      for (int iRel=0; iRel<sizRelIndex; iRel++) {
	int jCol=ListRelIndex[iRel];
	T fVal=eMatW(iRow, jCol);
	eSum += eVect[iRel]*fVal;
      }
      ListX.push_back(eSum);
    }
    GCD_int<T> eGCD=ComputeGCD_information(ListX);
    std::vector<std::vector<T> > NewBasis;
    for (int iVect=0; iVect<dimSpace; iVect++) {
      std::vector<T> kerVect=eGCD.ListListK[iVect];
      std::vector<T> eVectNew(sizRelIndex+1,0);
      eVectNew[0]=kerVect[0];
      for (int i=1; i<=dimSpace; i++) {
	T fVal=kerVect[i];
	std::vector<T> basVect=TheBasis[i-1];
	for (int j=0; j<sizRelIndex; j++)
	  eVectNew[j+1] += fVal*basVect[j];
      }
      NewBasis.push_back(eVectNew);
    }
    TheBasis=NewBasis;
  }
  MyMatrix<T> retNSP(dimSpace,nbCol);
  for (int iVect=0; iVect<dimSpace; iVect++) {
    std::vector<T> eVect=TheBasis[iVect];
    int idx=0;
    for (int iRank=0; iRank<eRank; iRank++) {
      int iCol=ListIndex[iRank];
      retNSP(iVect, iCol)=eVect[idx];
      idx++;
    }
    for (int iDim=0; iDim<dimSpace; iDim++) {
      int iCol=ListNonIndex[iDim];
      retNSP(iVect, iCol)=eVect[idx];
      idx++;
    }
  }
  for (int iVect=0; iVect<dimSpace; iVect++)
    for (int iRow=0; iRow<nbRow; iRow++) {
      T eSum=0;
      for (int iCol=0; iCol<nbCol; iCol++)
	eSum += eMat(iRow, iCol) * retNSP(iVect, iCol);
      if (eSum != 0) {
	std::cerr << "There are remaining errors in NullspaceIntTrMat\n";
	throw TerminalException{1};
      }
    }
  return retNSP;
}

template<typename T>
MyMatrix<T> NullspaceIntMat(MyMatrix<T> const& eMat)
{
  return NullspaceIntTrMat(TransposedMat(eMat));
}



// Given a vector v of T^n which is primitive
// the function returns a family of vector v1, ...., v(n-1)
// such that (v1, ...., v(n-1), v) is a T-basis of T^n
template<typename T>
MyMatrix<T> ComplementToBasis(MyVector<T> const& TheV)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain");
  using Treal=typename underlying_totally_ordered_ring<T>::real_type;
  std::vector<int> OperCol1;
  std::vector<int> OperCol2;
  std::vector<T> OperCoef;
  Treal AbsVal=-400;
  MyVector<T> TheVcopy=TheV;
  int n=TheV.size();
  int idxSelect;
  while(true) {
    bool IsFirst=true;
    int nbDiffZero=0;
    idxSelect=-1;
    for (int i=0; i<n; i++)
      if (TheVcopy(i) != 0) {
	nbDiffZero++;
	Treal hVal=T_NormGen(TheVcopy(i));
	if (IsFirst) {
	  IsFirst=false;
	  AbsVal=hVal;
	  idxSelect=i;
	}
	else {
	  if (hVal < AbsVal)
	    {
	      idxSelect=i;
	      AbsVal=hVal;
	    }
	}
      }
    if (idxSelect == -1) {
      std::cerr << "Inconsistency in computation of value\n";
      throw TerminalException{1};
    }
    if (nbDiffZero == 1) {
      if (AbsVal != 1) {
	std::cerr << "Wrong value for AbsVal\n";
	throw TerminalException{1};
      }
      break;
    }
    for (int j=0; j<n; j++)
      if (TheVcopy(j) != 0 && idxSelect !=j)
	{
	  T eVal1=TheVcopy(idxSelect);
	  T eVal2=TheVcopy(j);
	  T TheQ=QuoInt(eVal2, eVal1);
	  T res=eVal2 - TheQ*eVal1;
	  TheVcopy(j)=res;
	  OperCol1.push_back(idxSelect);
	  OperCol2.push_back(j);
	  OperCoef.push_back(TheQ);
	}
  }
  int nbOper=OperCol1.size();
  MyMatrix<T> TheReturn=MyMatrix<T>(n-1, n);
  int idx=0;
  for (int i=0; i<n; i++)
    if (i != idxSelect) {
      for (int j=0; j<n; j++) {
	T eVal=0;
	if (j == i)
	  eVal=1;
	TheReturn(idx, j)=eVal;
      }
      idx++;
    }
  for (int iOper=0; iOper<nbOper; iOper++) {
    int jOper=nbOper-1-iOper;
    int i1=OperCol1[jOper];
    int i2=OperCol2[jOper];
    T eCoeff=OperCoef[jOper];
    for (int iRow=0; iRow <n-1; iRow++) {
      T eVal=TheReturn(iRow, i2) + eCoeff*TheReturn(iRow, i1);
      TheReturn(iRow, i2)=eVal;
    }
    T eVal=TheVcopy(i2) + eCoeff*TheVcopy(i1);
    TheVcopy(i2)=eVal;
  }
  if (!TestEquality(TheVcopy, TheV)) {
    std::cerr << "TheVcopy =";
    WriteVector(std::cerr, TheVcopy);
    std::cerr << "TheV =";
    WriteVector(std::cerr, TheV);
    std::cerr << "Bookkeeping error\n";
    throw TerminalException{1};
  }
  return TheReturn;
}


// We have two matrices M1 and M2 and we check if they define
// the same subspace of T^n 
template<typename T>
bool TestEqualitySpaces(MyMatrix<T> const& M1, MyMatrix<T> const& M2)
{
  using Treal=typename underlying_totally_ordered_ring<T>::real_type;
  int idxSelect=-1;
  int k=M1.rows();
  int n=M1.cols();
  MyMatrix<T> M1copy=M1;
  MyMatrix<T> M2copy=M2;
  std::vector<int> StatusRow(k, 0);
  int idxSearch=0;
  for (int iK=0; iK<k; iK++) {
    while(true) {
      int nbDiff=0;
      for (int j=0; j<k; j++)
	if (StatusRow[j] == 0)
	  if (M1copy(j, idxSearch) != 0)
	    nbDiff++;
      if (nbDiff > 0)
	break;
      idxSearch++;
    }
    while(true) {
      int nbDiff=0;
      bool IsFirst=true;
      Treal AbsVal=0;
      for (int j=0; j<k; j++) {
	T eVal=M1copy(j,idxSearch);
	if (eVal != 0 && StatusRow[j] == 0) {
	  nbDiff++;
	  Treal hVal=T_NormGen(eVal);
	  if (IsFirst) {
	    IsFirst=false;
	    AbsVal=hVal;
	    idxSelect=j;
	  }
	  else {
	    if (hVal < AbsVal) {
	      idxSelect=j;
	      AbsVal=hVal;
	    }
	  }
	}
      }
      if (nbDiff == 1)
	break;
      for (int j=0; j<k; j++)
	if (j != idxSelect) {
	  T eVal1=M1copy(idxSelect, idxSearch);
	  T eVal2=M1copy(j, idxSearch);
	  T TheQ=QuoInt(eVal2, eVal1);
	  M1copy.row(j)=M1copy.row(j) - TheQ*M1copy.row(idxSelect);
	}
    }
    StatusRow[idxSelect]=1;
    T eVal1=M1copy(idxSelect, idxSearch);
    for (int j=0; j<k; j++) {
      T eVal2=M2copy(j, idxSearch);
      T TheQ=QuoInt(eVal2, eVal1);
      T res=eVal2 - TheQ*eVal1;
      if (res != 0)
	return false;
      M2copy.row(j)=M2copy.row(j) - TheQ*M1copy.row(idxSelect);
    }
    idxSearch++;
  }
  for (int j=0; j<k; j++)
    for (int i=0; i<n; i++)
      if (M2copy(j, i) != 0)
	return false;
  return true;
}


template<typename T>
int PositionSubspace(std::vector<MyMatrix<T> > const& ListSubspace, MyMatrix<T> const& OneSubspace)
{
  int nbSpace=ListSubspace.size();
  for (int iSub=0; iSub<nbSpace; iSub++) {
    bool test=TestEqualitySpaces(ListSubspace[iSub], OneSubspace);
    if (test)
      return iSub;
  }
  return -1;
}

template<typename Ti, typename Td>
void FuncInsertSubspace(std::vector<MyMatrix<Ti> > & ListSubspace, std::vector<Td> & ListDet, MyMatrix<Ti> const& OneSubspace, Td const& OneDet)
{
  if (PositionSubspace(ListSubspace, OneSubspace) == -1) {
    ListSubspace.push_back(OneSubspace);
    ListDet.push_back(OneDet);
  }
}


// T1 is integer type and T2 is real kind type
template<typename Ti, typename Td>
void ReorderSubspaceDet(std::vector<MyMatrix<Ti> > &ListSubspace, std::vector<Td> &ListDet)
{
  int nbSub=ListSubspace.size();
  for (int i=0; i<nbSub-1; i++)
    for (int j=i+1; j<nbSub; j++)
      if (ListDet[i] > ListDet[j])
	{
	  MyMatrix<Ti> eSub=ListSubspace[i];
	  Td eDet=ListDet[i];
	  ListSubspace[i]=ListSubspace[j];
	  ListDet[i]=ListDet[j];
	  ListSubspace[j]=eSub;
	  ListDet[j]=eDet;
	}
}

/* test if ListSub2 is a subset of ListSub1 */
template<typename T>
bool TestInclusionFamilySubspace(std::vector<MyMatrix<T> > const& ListSub1, std::vector<MyMatrix<T> > const& ListSub2)
{
  int nbSub2=ListSub2.size();
  for (int i=0; i<nbSub2; i++)
    if (PositionSubspace(ListSub1, ListSub2[i]) == -1)
      return false;
  return true;
}


template<typename T>
bool TestEqualityFamilySubspace(std::vector<MyMatrix<T> > const& ListSub1, std::vector<MyMatrix<T> > const& ListSub2)
{
  int nbSub1=ListSub1.size();
  int nbSub2=ListSub2.size();
  if (nbSub1 != nbSub2)
    return false;
  return TestInclusionFamilySubspace(ListSub1, ListSub2);
}



template<typename T>
MyMatrix<T> GetNoncontainedSubspace(std::vector<MyMatrix<T> > const& ListSubBig, std::vector<MyMatrix<T> > const& ListSubSma)
{
  int nbSubBig=ListSubBig.size();
  int nbSubSma=ListSubSma.size();
  for (int i=0; i<nbSubBig; i++)
    if (PositionSubspace(ListSubSma, ListSubBig[i]) == -1)
      return ListSubBig[i];
  std::cerr << "We did not find any subspace\n";
  std::cerr << "nbSubBig=" << nbSubBig << " nbSumSma=" << nbSubSma << "\n";
  throw TerminalException{1};
}


template<typename T>
struct ResultSolutionIntMat {
  bool TheRes;
  MyVector<T> eSol;
};



// Find an integral solution of the equation Y = X A
// if it exists.
template<typename T>
ResultSolutionIntMat<T> SolutionIntMat(MyMatrix<T> const& TheMat, MyVector<T> const& TheVect)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain");
  using Treal=typename underlying_totally_ordered_ring<T>::real_type;
  int nbDiff;
  int nbVect=TheMat.rows();
  int nbCol=TheMat.cols();
  if (nbVect == 0) {
    MyVector<T> eSol;
    if (IsZeroVector(TheVect)) {
      return {true, eSol};
    }
    else {
      return {false, {}};
    }
  }
  //  std::cerr << "nbVect=" << nbVect << " nbCol=" << nbCol << "\n";
  MyVector<T> eSol=ZeroVector<T>(nbVect);
  MyMatrix<T> eEquivMat=IdentityMat<T>(nbVect);
  MyMatrix<T> TheMatWork=TheMat;
  MyVector<T> TheVectWork=TheVect;
  std::vector<int> VectStatus(nbVect,1);
  for (int i=0; i<nbCol; i++) {
    //    std::cerr << "i=" << i << " nbCol=" << nbCol << "\n";
    int iVectFound=-1;
    while(true) {
      bool IsFirst=true;
      Treal MinValue=0;
      nbDiff=0;
      for (int iVect=0; iVect<nbVect; iVect++)
        if (VectStatus[iVect] == 1) {
          T prov1=TheMatWork(iVect, i);
          Treal eNorm=T_NormGen(prov1);
	  //	  std::cerr << "iVect=" << iVect << " eNorm=" << eNorm << "\n";
	  if (prov1 != 0) {
	    nbDiff++;
	    if (IsFirst) {
	      IsFirst=false;
	      MinValue=eNorm;
	      iVectFound=iVect;
	    }
	    else {
	      if (eNorm < MinValue) {
		MinValue=eNorm;
		iVectFound=iVect;
	      }
	    }
	  }
	}
      if (nbDiff == 1 || nbDiff == 0)
	break;
      if (MinValue == 0) {
	std::cerr << "MinValue should not be zero\n";
	throw TerminalException{1};
      }
      //      std::cerr << "i=" << i << " iVectFound=" << iVectFound << " nbDiff=" << nbDiff << " MinValue=" << MinValue << "\n";
      for (int iVect=0; iVect<nbVect; iVect++)
	if (VectStatus[iVect] == 1 && iVect != iVectFound) {
	  T prov1b=TheMatWork(iVectFound, i);
	  T prov2=TheMatWork(iVect, i);
	  T TheQ=QuoInt(prov2, prov1b);
	  //	  std::cerr << "iVect=" << iVect << " prov1b=" << prov1b << " prov2=" << prov2 << " q=" << TheQ << "\n";
	  TheMatWork.row(iVect) -= TheQ*TheMatWork.row(iVectFound);
	  eEquivMat.row(iVect) -= TheQ*eEquivMat.row(iVectFound);
	}
      /*      std::cerr << "Now TheMatWork=\n";
	      WriteMatrix(std::cerr, TheMatWork);*/
      /*      MyMatrix<T> eProdMat=eEquivMat*TheMat;
      if (TheMatWork != eProdMat) {
	std::cerr << "Matrix error in the construction of eEquivMat and TheMatWork\n";
	throw TerminalException{1};
	}*/
    }
    if (nbDiff == 1) {
      if (iVectFound == -1) {
        std::cerr << "Clear error in the program\n";
	throw TerminalException{1};
      }
      VectStatus[iVectFound]=0;
      /*      std::cerr << "TheMatWork=\n";
	      WriteMatrix(std::cerr, TheMatWork);*/
      T prov1=TheVectWork(i);
      T prov2=TheMatWork(iVectFound, i);
      T TheQ=QuoInt(prov1, prov2);
      for (int j=0; j<nbCol; j++)
	TheVectWork(j) -= TheQ*TheMatWork(iVectFound, j);
      for (int iVect=0; iVect<nbVect; iVect++)
	eSol(iVect) += TheQ*eEquivMat(iVectFound,iVect);
    }
    if (TheVectWork(i) != 0)
      return {false, {}};
  }
  /*
  MyVector<T> eProd=ProductVectorMatrix(eSol,TheMat);
  for (int i=0; i<nbCol; i++)
    if (eProd(i) != TheVect(i)) {
      std::cerr << "eProd=\n";
      WriteVector(std::cerr, eProd);
      std::cerr << "TheVect=\n";
      WriteVector(std::cerr, TheVect);
      std::cerr << "SolutionIntMat: Error in the solutioning algorithm\n";
      throw TerminalException{1};
      }*/
  return {true, eSol};
}


template<typename T>
struct CanSolIntMat {
  std::vector<int> ListRow;
  MyMatrix<T> TheMatWork;
  MyMatrix<T> eEquivMat;
};




template<typename T>
CanSolIntMat<T> ComputeCanonicalFormFastReduction(MyMatrix<T> const& TheMat)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain");
  using Treal=typename underlying_totally_ordered_ring<T>::real_type;
  int nbDiff;
  int nbVect=TheMat.rows();
  int nbCol=TheMat.cols();
  if (nbVect == 0) {
    std::cerr << "Need to write the code here\n";
    throw TerminalException{1};
  }
  MyMatrix<T> eEquivMat=IdentityMat<T>(nbVect);
  MyMatrix<T> TheMatWork=TheMat;
  std::vector<int> VectStatus(nbVect,1);
  std::vector<int> ListRow(nbCol);
  for (int i=0; i<nbCol; i++) {
    //    std::cerr << "i=" << i << " nbCol=" << nbCol << "\n";
    int iVectFound=-1;
    while(true) {
      bool IsFirst=true;
      Treal MinValue=0;
      nbDiff=0;
      for (int iVect=0; iVect<nbVect; iVect++)
        if (VectStatus[iVect] == 1) {
          T prov1=TheMatWork(iVect, i);
          Treal eNorm=T_NormGen(prov1);
	  //	  std::cerr << "iVect=" << iVect << " eNorm=" << eNorm << "\n";
	  if (prov1 != 0) {
	    nbDiff++;
	    if (IsFirst) {
	      IsFirst=false;
	      MinValue=eNorm;
	      iVectFound=iVect;
	    }
	    else {
	      if (eNorm < MinValue) {
		MinValue=eNorm;
		iVectFound=iVect;
	      }
	    }
	  }
	}
      if (nbDiff == 1 || nbDiff == 0)
	break;
      if (MinValue == 0) {
	std::cerr << "MinValue should not be zero\n";
	throw TerminalException{1};
      }
      for (int iVect=0; iVect<nbVect; iVect++)
	if (VectStatus[iVect] == 1 && iVect != iVectFound) {
	  T prov1b=TheMatWork(iVectFound, i);
	  T prov2=TheMatWork(iVect, i);
	  T TheQ=QuoInt(prov2, prov1b);
	  TheMatWork.row(iVect) -= TheQ*TheMatWork.row(iVectFound);
	  eEquivMat.row(iVect) -= TheQ*eEquivMat.row(iVectFound);
	}
      /*      MyMatrix<T> eProdMat=eEquivMat*TheMat;
      if (TheMatWork != eProdMat) {
	std::cerr << "Matrix error in the construction of eEquivMat and TheMatWork\n";
	throw TerminalException{1};
	}*/
    }
    int eVal;
    if (nbDiff == 1) {
      eVal=iVectFound;
      //      ListRow.push_back(iVectFound);
      //      ListCol.push_back(i);
      if (iVectFound == -1) {
        std::cerr << "Clear error in the program\n";
	throw TerminalException{1};
      }
      VectStatus[iVectFound]=0;
    }
    else {
      eVal=-1;
    }
    ListRow[i]=eVal;
  }
  return {ListRow, TheMatWork, eEquivMat};
}

template<typename T>
bool CanTestSolutionIntMat(CanSolIntMat<T> const& eCan, MyVector<T> const& TheVect)
{
  int nbCol=eCan.TheMatWork.cols();
  MyVector<T> TheVectWork=TheVect;
  for (int i=0; i<nbCol; i++) {
    int iRow=eCan.ListRow[i];
    if (iRow >= 0) {
      T prov1=TheVectWork(i);
      T prov2=eCan.TheMatWork(iRow, i);
      T TheQ=QuoInt(prov1, prov2);
      for (int j=0; j<nbCol; j++)
	TheVectWork(j) -= TheQ*eCan.TheMatWork(iRow, j);
    }
    if (TheVectWork(i) != 0)
      return false;
  }
  return true;
}

template<typename T>
ResultSolutionIntMat<T> CanSolutionIntMat(CanSolIntMat<T> const& eCan, MyVector<T> const& TheVect)
{
  int nbVect=eCan.TheMatWork.rows();
  int nbCol=eCan.TheMatWork.cols();
  MyVector<T> TheVectWork=TheVect;
  MyVector<T> eSol=ZeroVector<T>(nbVect);
  for (int i=0; i<nbCol; i++) {
    int iRow=eCan.ListRow[i];
    if (iRow >= 0) {
      T prov1=TheVectWork(i);
      T prov2=eCan.TheMatWork(iRow, i);
      T TheQ=QuoInt(prov1, prov2);
      for (int j=0; j<nbCol; j++)
	TheVectWork(j) -= TheQ*eCan.TheMatWork(iRow, j);
      for (int iVect=0; iVect<nbVect; iVect++)
	eSol(iVect) += TheQ*eCan.eEquivMat(iRow,iVect);
    }
    if (TheVectWork(i) != 0)
      return {false,{}};
  }
  return {true,eSol};
}



template<typename T>
struct BasisReduction {
  MyMatrix<T> TheBasisReduced;
  MyMatrix<T> Pmat;
  std::vector<int> IdxVector;
  MyMatrix<T> TheBasisReord;
};


// We need TheBasis to be of full rank.
// This code is for the Copositivity.
//
// Following operation is done on the matrix
// The matrix TheBasis is transformed by operations
// on the columns in order to diminish the size
// of the coefficients.
//
// Specifically, the matrix TheBasisReord is changed so that
// it is of the form
// (a11 0 ....... 0)
// ( x  a22 0 ... 0)
//     .
//     .
// ( x  .... x  ann)
// We have aII > 0
// This is done in two ways
// Applying an integral matrix to the column
// and then a permutation on the rows.
template<typename T>
BasisReduction<T> ComputeBasisReduction(MyMatrix<T> const& TheBasis)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain");
  using Treal=typename underlying_totally_ordered_ring<T>::real_type;
  int nbCol=TheBasis.cols();
  int nbRow=TheBasis.rows();
  std::vector<int> colStat(nbCol, 1);
  std::vector<int> rowStat(nbRow, 1);
  MyMatrix<T> TheBasisReduced = TheBasis;
  MyMatrix<T> Pmat=IdentityMat<T>(nbCol);
  MyMatrix<T> TheBasisReord(nbRow, nbCol);
  //  std::cerr << "Before writing TheBasis\n";
  //  WriteMatrix(std::cerr, TheBasis);
  //  std::cerr << " After writing TheBasis\n";
  std::vector<int> IdxVector;
  auto FindMinGCDrow=[&](int const& iRank) -> int {
    int iRowSearch=-1;
    bool IsFirst=true;
    Treal AbsVal=0;
    for (int iRow=0; iRow<nbRow; iRow++)
      if (rowStat[iRow] == 1) {
	std::vector<T> eRowRed(nbCol - iRank);
	//	std::cerr << "eRowRed=";
	for (int iCol=0; iCol<nbCol - iRank; iCol++) {
	  T eVal=TheBasisReduced(iRow, iCol+iRank);
	  //	  std::cerr << " " << eVal;
	  eRowRed[iCol]=eVal;
	}
	//	std::cerr << "\n";
	//	std::cerr << "Before GCD computation\n";
	GCD_int<T> eGCD=ComputeGCD_information(eRowRed);
	//	std::cerr << "After GCD computation\n";
	//	WriteGCD_int(std::cerr, eGCD);
	Treal eNorm=T_NormGen(eGCD.gcd);
	if (IsFirst) {
	  IsFirst=false;
	  AbsVal=eNorm;
	  iRowSearch=iRow;
	}
	else {
	  if (eNorm < AbsVal) {
	    AbsVal=eNorm;
	    iRowSearch=iRow;
	  }
	}
      }
    return iRowSearch;
  };
  auto SingleMultiplicationUpdate=[&](MyMatrix<T> const& PartMat) -> void {
    //    std::cerr << "SingleMultiplicationUpdate\n";
    //    std::cerr << "PartMat=\n";
    //    WriteMatrix(std::cerr, PartMat);
    Pmat = Pmat*PartMat;
    //    std::cerr << "Now Pmat=\n";
    //    WriteMatrix(std::cerr, Pmat);
    TheBasisReduced = TheBasisReduced*PartMat;
    //    std::cerr << "Now TheBasisReduced=\n";
    //    WriteMatrix(std::cerr, TheBasisReduced);
  };
  auto UpdateMatrices=[&](int const& iRank, int const& iRowSearch) -> void {
    rowStat[iRowSearch]=0;
    IdxVector.push_back(iRowSearch);
    std::vector<T> eRowRed(nbCol - iRank);
    for (int iCol=0; iCol<nbCol - iRank; iCol++)
      eRowRed[iCol]=TheBasisReduced(iRowSearch, iCol+iRank);
    /*    std::cerr << "eRowRed=";
    for (auto & eVal : eRowRed)
      std::cerr << " " << eVal;
      std::cerr << "\n";*/
    GCD_int<T> eGCD=ComputeGCD_information(eRowRed);
    //    std::cerr << "eGCD.Pmat=\n";
    //    WriteMatrix(std::cerr, eGCD.Pmat);
    MyMatrix<T> PartMat=IdentityMat<T>(nbCol);
    //    std::cerr << "Before operations: iRank=" << iRank << "\n";
    for (int iCol=iRank; iCol<nbCol; iCol++)
      for (int iRow=iRank; iRow<nbCol; iRow++) {
	//	std::cerr << "iRow=" << iRow << " iCol=" << iCol << "\n";
	T eVal=eGCD.Pmat(iRow-iRank, iCol-iRank);
	//	std::cerr << "We have eVal\n";
	PartMat(iRow, iCol)=eVal;
	//	std::cerr << "Assign PartMat\n";
      }
    //    std::cerr << "Expression of PartMat:\n";
    //    WriteMatrix(std::cerr, PartMat);

    SingleMultiplicationUpdate(PartMat);
    for (int iCol=0; iCol<iRank; iCol++) {
      T a=TheBasisReduced(iRowSearch, iCol);
      T b=TheBasisReduced(iRowSearch, iRank);
      T q=QuoInt(a, b);
      //      std::cerr << "q=" << q << "\n";
      MyMatrix<T> PartMatB=IdentityMat<T>(nbCol);
      PartMatB(iRank, iCol)=-q;
      SingleMultiplicationUpdate(PartMatB);
    }
    if (TheBasisReduced(iRowSearch, iRank) < 0) {
      MyMatrix<T> PartMatC=IdentityMat<T>(nbCol);
      PartMatC(iRank, iRank)=-1;
      SingleMultiplicationUpdate(PartMatC);
    }
    //    std::cerr << "iRank=" << iRank << " DiagVal=" << TheBasisReduced(iRowSearch, iRank) << "\n";
    TheBasisReord.row(iRank)=TheBasisReduced.row(iRowSearch);
  };
  int TheRank=std::min(nbCol, nbRow);
  for (int iRank=0; iRank<TheRank; iRank++) {
    //    std::cerr << "iRank=" << iRank << "\n";
    int iRowSearch=FindMinGCDrow(iRank);
    //    std::cerr << "iRowSearch=" << iRowSearch << "\n";
    UpdateMatrices(iRank, iRowSearch);
  }
  //  std::cerr << "Pmat=\n";
  //  WriteMatrix(std::cerr, Pmat);
  BasisReduction<T> eRed{TheBasisReduced, Pmat, IdxVector, TheBasisReord};
  return eRed;
}









struct AffineBasisResult {
  bool result;
  std::vector<int> ListIdx;
};


// Given a family of points EXT, find a subset (v1, ...., vN) of EXT
// such that for every point of v of EXT there exist lambda1, ..., lambdaN in T
// with v=lambda1 v1 + .... + lambdaN vN
// This may not exist 
template<typename T>
AffineBasisResult Kernel_ComputeAffineBasis(MyMatrix<T> const& EXT)
{
  static_assert(is_ring_field<T>::value, "Requires T to have inverses");
  int nbRow=EXT.rows();
  int nbCol=EXT.cols();
  int n=nbCol;
  MyMatrix<T> ListExp=ZeroMatrix<T>(nbRow, nbCol);
  MyMatrix<T> EXTwork=EXT;
  std::vector<int> RowStatus(nbRow, 0);
  std::vector<int> ColumnStatus(nbCol,1);
  std::cerr << "Starting Kernel_ComputeAffineBasis\n";
  MyVector<T> V1(nbCol);
  MyVector<T> V2(nbCol);
  MyVector<T> eExpr(nbCol);
  auto fInsertValue=[&](int const& idx, int const& iVect) -> bool {
    int eCol=-1;
    for (int iCol=0; iCol<nbCol; iCol++)
      if (eCol == -1 && EXTwork(iVect, iCol) != 0 && ColumnStatus[iCol] == 1)
	eCol=iCol;
    std::cerr << "eCol=" << eCol << "\n";
    if (eCol == -1) {
      std::cerr << "This should not be selected\n";
      std::cerr << "nbCol=" << nbCol << "\n";
      for (int iCol=0; iCol<nbCol; iCol++) {
	std::cerr << " iCol=" << iCol << " stat=" << ColumnStatus[iCol] << " val=" << EXTwork(iVect,iCol) << "\n";
      }
      throw TerminalException{1};
    }
    V1=EXTwork.row(iVect);
    std::cerr << "V1 rows=" << V1.rows() << " cols=" << V1.cols() << "\n";
    ListExp(iVect, idx)=1;
    //    std::cerr << "fInsertValue, step 1\n";
    for (int iRow=0; iRow<nbRow; iRow++)
      if (RowStatus[iRow] == 0) {
	V2=EXTwork.row(iRow);
	bool test=IsVectorMultiple(V1, V2);
	if (test) {
	  T eQuot=V2(eCol)/V1(eCol);
	  eExpr=ListExp.row(iRow) - eQuot*ListExp.row(iVect);
	  for (int iCol=0; iCol<nbCol; iCol++)
	    if (!IsInteger(eExpr(iCol)))
	      return false;
	}
      }
    ColumnStatus[eCol]=0;
    //    std::cerr << "fInsertValue, step 2\n";
    RowStatus[iVect]=1;
    for (int iRow=0; iRow<nbRow; iRow++)
      if (RowStatus[iRow] == 0) {
	//	std::cerr << "iRow=" << iRow << "\n";
	V2=EXTwork.row(iRow);
	//	std::cerr << "V2 extracted\n";
	T eQuot=V2(eCol)/V1(eCol);
	//	std::cerr << "eQuot extracted\n";
	ListExp.row(iRow)=ListExp.row(iRow) - eQuot*ListExp.row(iVect);
	//	std::cerr << "After ListExp assignation\n";
	EXTwork.row(iRow)=EXTwork.row(iRow) - eQuot*EXTwork.row(iVect);
	//	std::cerr << "After EXTwork assignation\n";
	V2=EXTwork.row(iRow);
	//	std::cerr << "After eVect extraction\n";
	bool IsZero=IsZeroVector(V2);
	int iRowPrint=-2249;
	if (iRow == iRowPrint) {
	  std::cerr << iRowPrint << ": eVect\n";
	  for (int jCol=0; jCol<nbCol; jCol++)
	    std::cerr << " jCol=" << jCol << " stat=" << ColumnStatus[jCol] << " val=" << EXTwork(iRowPrint,jCol) << "\n";
	  std::cerr << iRowPrint << ": IsZero=" << IsZero << "\n";
	}
	//	std::cerr << "After IsZero check\n";
	if (IsZero)
	  RowStatus[iRow]=1;
      }
    int nbFinished=0;
    for (int iRow=0; iRow<nbRow; iRow++)
      if (RowStatus[iRow] == 1)
	nbFinished++;
    std::cerr << "nbFinished=" << nbFinished << "\n";
    /*
    for (int iRow=0; iRow<nbRow; iRow++)
      for (int iCol=0; iCol<nbCol; iCol++)
	if (ColumnStatus[iCol] == 0 && RowStatus[iRow] == 0)
	  if (EXTwork(iRow, iCol) != 0) {
	    std::cerr << "Inconsistency error\n";
	    std::cerr << "eCol=" << eCol << "\n";
	    std::cerr << "idx=" << idx << "\n";
	    std::cerr << "iVect=" << iVect << "\n";
	    for (int jCol=0; jCol<nbCol; jCol++)
	      std::cerr << " jCol=" << jCol << " stat=" << ColumnStatus[jCol] << " val=" << EXTwork(iRow,jCol) << "\n";
	    std::cerr << "nbRow=" << nbRow << "\n";
	    std::cerr << "nbCol=" << nbCol << "\n";
	    throw TerminalException{1};
	    }
    */
    //    std::cerr << "fInsertValue, step 3\n";
    return true;
  };
  int nbIter=1000;
  std::vector<int> ListIdx(n);
  std::vector<int> UsedNumber(nbRow, 0);
  auto GetRandomNumber=[&]() -> int {
    for (int iter=0; iter<nbIter; iter++) {
      int eVal=rand() % nbRow;
      if (UsedNumber[eVal] == 0 && RowStatus[eVal] == 0)
	return eVal;
    }
    return -1;
  };
  auto SetLocallyCorrectIndex=[&](int const& idx) -> int {
    for (int iter=0; iter<nbIter; iter++) {
      int eVal=GetRandomNumber();
      //      std::cerr << "SetLocallyCorrectIndex eVal=" << eVal << "\n";
      if (eVal == -1)
	return -1;
      //      std::cerr << "Before fInsertValue eVal=" << eVal << "\n";
      bool res=fInsertValue(idx, eVal);
      //      std::cerr << "After fInsertValue\n";
      UsedNumber[eVal]=1;
      if (res) {
	ListIdx[idx]=eVal;
	return 0;
      }
    }
    return -1;
  };
  for (int i=0; i<n; i++) {
    std::cerr << "i=" << i << "\n";
    //    std::cerr << "Before SetLocallyCorectIndex\n";
    //    WriteMatrix(std::cerr, EXTwork);
    int eVal=SetLocallyCorrectIndex(i);
    //    std::cerr << "i=" << i << "\n";
    //    std::cerr << "Reduced Matrix form\n";
    //    WriteMatrix(std::cerr, EXTwork);
    if (eVal == -1)
      return {false, {}};
  }
  return {true, ListIdx};
}

template<typename T> 
AffineBasisResult ComputeAffineBasis(MyMatrix<T> const& EXT)
{
  int nbIter=1000;
  for (int iter=0; iter<nbIter; iter++) {
    AffineBasisResult eAffRes=Kernel_ComputeAffineBasis(EXT);
    if (eAffRes.result)
      return eAffRes;
  }
  return {false, {}};
}


// Compute the translation classes.
// Two classes eV and fV are equivalent if there exists a vector w integer such that
// eV - fV = w M
// that is we need to compute (eV - fV) M^(-1)
template<typename T>
std::vector<MyVector<int>> ComputeTranslationClasses(MyMatrix<T> const& M)
{
  int n=M.rows();
  MyMatrix<T> eInv=Inverse(M);
  std::vector<MyVector<int>> ListClasses;
  std::vector<int> ListStatus;
  auto IsEquivalent=[&](MyVector<int> const& eV, MyVector<int> const& fV) -> bool {
    MyVector<int> diff=eV - fV;
    for (int i=0; i<n; i++) {
      T eVal=0;
      for (int j=0; j<n; j++)
	eVal += diff(j) * eInv(j,i);
      if (!IsInteger(eVal))
	return false;
    }
    return true;
  };
  auto FuncInsert=[&](MyVector<int> const& eV) -> void {
    for (auto & fV : ListClasses) {
      if (IsEquivalent(eV, fV))
	return;
    }
    ListClasses.push_back(eV);
    ListStatus.push_back(1);
  };
  MyVector<int> zerV(n);
  for (int i=0; i<n; i++)
    zerV(i)=0;
  FuncInsert(zerV);
  while(true) {
    bool IsFinished=true;
    int nbClass=ListClasses.size();
    for (int iClass=0; iClass<nbClass; iClass++) {
      if (ListStatus[iClass] == 1) {
	ListStatus[iClass]=0;
	IsFinished=false;
	MyVector<int> eClass=ListClasses[iClass];
	for (int i=0; i<n; i++) {
	  MyVector<int> fClass=eClass;
	  fClass[i]++;
	  FuncInsert(fClass);
	}
      }
    }
    if (IsFinished)
      break;
  }
  return ListClasses;

  
}




// Given a family of vector v1, ....., vN
// Find an integral family of vector w1, ...., wR
// such that all vi are expressed integrally in term of wi
// and uniquely.
// This always exist as opposed to affine basis.
template<typename T>
MyMatrix<T> GetZbasis(MyMatrix<T> const& ListElement)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain");
  using Treal=typename underlying_totally_ordered_ring<T>::real_type;
  int TheDim=ListElement.cols();
  MyMatrix<T> ListEqua;
  MyMatrix<T> InvMatrix;
  MyMatrix<T> InvMatrixTr;
  MyMatrix<T> TheBasis;
  std::vector<int> eSet;
  auto fGetOneBasis=[&](MyVector<T> const& eSol) -> MyMatrix<T> {
    int DimLoc=TheBasis.rows();
    //    std::cerr << "DimLoc=" << DimLoc << " |eSol|=" << eSol.size() << "\n";
    MyMatrix<T> TheRedMat=ZeroMatrix<T>(DimLoc+1, DimLoc);
    for (int i=0; i<DimLoc; i++)
      TheRedMat(i,i)=1;
    for (int i=0; i<DimLoc; i++)
      TheRedMat(DimLoc,i)=eSol(i);
    //    std::cerr << "After TheRedMat construction\n";
    MyMatrix<T> NSP=NullspaceIntMat(TheRedMat);
    //    std::cerr << "We have NSP\n";
    if (NSP.rows() != 1) {
      std::cerr << "|NSP|=" << NSP.rows() << " when it should be 1\n";
      std::cerr << "TheRedMat:\n";
      WriteMatrix(std::cerr, TheRedMat);
      std::cerr << "TheRedMat:\n";
      WriteMatrix(std::cerr, NSP);
      std::cerr << "Inconsistency that needs to be corrected\n";
      throw TerminalException{1};
    }
    MyVector<T> eVect=GetMatrixRow(NSP,0);
    /*    std::cerr << "eVect=";
	  WriteVector(std::cerr, eVect);*/
    int n2=DimLoc+1;
    while(true) {
      std::vector<int> ListIdxZ;
      std::vector<int> ListIdxNZ;
      for (int i=0; i<n2; i++) {
	if (eVect(i) == 0) {
	  ListIdxZ.push_back(i);
	}
	else {
	  ListIdxNZ.push_back(i);
	}
      }
      if (ListIdxNZ.size() == 1)
	return SelectRow(TheRedMat, ListIdxZ);
      std::vector<int> AbsList;
      bool IsFirst=true;
      int ThePivot=-1;
      Treal TheMin=-1;
      for (auto & eVal : ListIdxNZ) {
	T eVal_T=eVect(eVal);
	Treal eAbs=T_NormGen(eVal_T);
	if (IsFirst) {
	  TheMin=eAbs;
	  ThePivot=eVal;
	}
	else {
	  if (eAbs < TheMin) {
	    TheMin=eAbs;
	    ThePivot=eVal;
	  }
	}
	IsFirst=false;
      }
      for (int iCol=0; iCol<n2; iCol++)
	if (iCol != ThePivot) {
	  T q=QuoInt(eVect(iCol), eVect(ThePivot));
	  TheRedMat.row(ThePivot)=TheRedMat.row(ThePivot) + q*TheRedMat.row(iCol);
	  eVect(iCol)=eVect(iCol) - q*eVect(ThePivot);
	}
    }
  };
  auto fComputeSpeed=[&]() -> void {
    int dimSpace=TheBasis.rows();
    if (dimSpace == 0) {
      ListEqua=IdentityMat<T>(TheDim);
      eSet={};
    }
    else {
      /*      std::cerr << "TheBasis=\n";
	      WriteMatrix(std::cerr, TheBasis);*/
      ListEqua=NullspaceTrMat(TheBasis);
      eSet=ColumnReductionSet(TheBasis);
      //      std::cerr << "|eSet|=" << eSet.size() << "\n";
      InvMatrix=Inverse(SelectColumn(TheBasis, eSet));
      InvMatrixTr=InvMatrix.transpose();
    }
  };
  auto IsInSpace=[&](MyVector<T> const& eElt) -> bool {
    int nbEqua=ListEqua.rows();
    //    std::cerr << "nbEqua=" << nbEqua << "\n";
    for (int iEqua=0; iEqua<nbEqua; iEqua++) {
      T eSum=0;
      for (int i=0; i<TheDim; i++)
	eSum += ListEqua(iEqua,i)*eElt(i);
      if (eSum != 0)
	return false;
    }
    return true;
  };
  auto fInsert=[&](MyVector<T> const& eElt) -> void {
    bool test=IsInSpace(eElt);
    //    std::cerr << "test=" << test << "\n";
    if (!test) {
      TheBasis=ConcatenateMatVec(TheBasis, eElt);
      fComputeSpeed();
    }
    else {
      if (TheBasis.rows() == 0)
	return;
      MyVector<T> eEltRed=SelectColumnVector(eElt, eSet);
      MyVector<T> eSol=InvMatrixTr*eEltRed;
      if (IsIntegralVector(eSol))
	return;
      MyMatrix<T> NewBasis=fGetOneBasis(eSol);
      TheBasis=NewBasis*TheBasis;
      fComputeSpeed();
    }
  };
  fComputeSpeed();
  int nbElt=ListElement.rows();
  //  std::cerr << "nbElt=" << nbElt << "\n";
  for (int iElt=0; iElt<nbElt; iElt++) {
    //    std::cerr << "iElt=" << iElt << "\n";
    MyVector<T> eElt=GetMatrixRow(ListElement, iElt);
    fInsert(eElt);
    //    std::cerr << "After fInsert\n";
  }
  bool DoCheck=true;
  if (DoCheck) {
    int DimSpace=TheBasis.rows();
    for (int iBas=0; iBas<DimSpace; iBas++) {
      MyVector<T> eLine=GetMatrixRow(TheBasis, iBas);
      //      std::cerr << "Before SolutionIntMat, iBas=" << iBas << "\n";
      ResultSolutionIntMat<T> eResIntMat=SolutionIntMat(ListElement, eLine);
      /*      std::cerr << "ListElement=\n";
      WriteMatrixGAP(std::cerr, ListElement);
      std::cerr << "eLine=\n";
      WriteVectorGAP(std::cerr, eLine);
      std::cerr << "After SolutionIntMat 1\n";*/
      if (!eResIntMat.TheRes) {
	std::cerr << "Error in GetZbasis 1\n";
	throw TerminalException{1};
      }
    }
    for (int iElt=0; iElt<nbElt; iElt++) {
      MyVector<T> eElt=GetMatrixRow(ListElement, iElt);
      //      std::cerr << "Before SolutionIntMat, iElt=" << iElt << "\n";
      ResultSolutionIntMat<T> eResIntMat=SolutionIntMat(TheBasis, eElt);
      /*      std::cerr << "TheBasis=\n";
      WriteMatrixGAP(std::cerr, TheBasis);
      std::cerr << "eElt=\n";
      WriteVectorGAP(std::cerr, eElt);
      std::cerr << "After SolutionIntMat 2 eResIntMat.TheRes=" << eResIntMat.TheRes << "\n";*/
      if (!eResIntMat.TheRes) {
	std::cerr << "Error in GetZbasis 2\n";
	throw TerminalException{1};
      }
    }
  }
  return TheBasis;
}








#endif

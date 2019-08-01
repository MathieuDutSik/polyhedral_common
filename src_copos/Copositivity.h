#ifndef COPOSITIVITY_INCLUDE
#define COPOSITIVITY_INCLUDE

#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "Temp_Positivity.h"
#include "POLY_LinearProgramming.h"
#include "LatticeDefinitions.h"



template<typename T>
bool TestCopositivityByPositivityCoeff(MyMatrix<T> const& eSymmMatB)
{
  int n=eSymmMatB.rows();
  for (int i=0; i<n; i++) {
    T eNorm=eSymmMatB(i,i);
    if (eNorm <= 0)
      return false;
  }
  for (int i=0; i<n-1; i++)
    for (int j=i+1; j<n; j++) {
      T eScal=eSymmMatB(i,j);
      if (eScal < 0)
	return false;
    }
  return true;
}





template<typename T>
struct CopositivityInfoReduction {
  MyMatrix<T> TheBasisReduced;
  MyMatrix<T> eSymmMatReduced;
  MyMatrix<T> Pmat;
  std::vector<int> IdxVector;
  std::vector<int> IdxColumn;
};



//#define DEBUG_UNDER_POS_COND

// Given C >0 a >=0 and b>0 find the maximum
// value of x>=0 such that a x + b x^2 <= C
// reduced form is a2 x + x^2 <= C2
// polynomial is   x^2 + a2 x - C2 = 0
//               a x^2 +  b x + c  = 0
template<typename T>
int FindLargest(T const& a, T const& b, T const& C)
{
  T C2=C/b;
  T a2=a/b;
  double a2_doubl = UniversalTypeConversion<double,T>(a2);
  double C2_doubl = UniversalTypeConversion<double,T>(C2);
#ifdef DEBUG_UNDER_POS_COND
  if (b <= 0 || C <= 0) {
    std::cerr << "b should be strictly positive. b=" << b << "\n";
    std::cerr << "C should be strictly positive. C=" << C << "\n";
    throw TerminalException{1};
  }
#endif
  double delta=a2_doubl*a2_doubl + 4*C2_doubl;
  double x1=0.5*(-a2_doubl + sqrt(delta));
  double eD1=floor(x1);
  long int eD2=lround(eD1);
  int eReturn=eD2;
  auto f=[&](int const& x) -> bool {
    T eDiff=C2 - a2*x - x*x;
    if (eDiff >= 0)
      return true;
    return false;
  };
  while(true) {
    bool test1=f(eReturn);
    bool test2=f(eReturn+1);
    if (test1 && !test2)
      break;
    if (!test1)
      eReturn--;
    if (test2)
      eReturn++;
  }
  return eReturn;
}




// This function is the kernel of our methodology
// We have a symmetric matrix eSymmMat and 
// a family of vectors in TheBasis.
//
// Assumptions
// ---TheBasis is a set of vectors (v1, ..., vn)
// ---We have vi in the positive orthant R_{\geq 0}^n
// ---(v1, ...., vn) spans a simplicial cone. We have no constraint on det(v1, ..., vn). Just be non-zero
// ---For all 1<= i,j <= n  we have vi eSymmMat * vj >= 0
template<typename T, typename Tint>
std::vector<MyVector<Tint>> EnumerateShortVectorInCone_UnderPositivityCond(MyMatrix<T> const& eSymmMat, MyMatrix<Tint> const& TheBasis, T const& MaxNorm)
{
#ifdef DEBUG_UNDER_POS_COND
  MyMatrix<T> tstTheBasis_T=ConvertMatrixUniversal<T,Tint>(TheBasis);
  MyMatrix<T> tstSymmMatB=TheBasis_T*eSymmMat*TheBasis_T.transpose();
  bool test1=TestCopositivityByPositivityCoeff(tstSymmMatB);
  if (!test1) {
    std::cerr << "Inconsistency in the positivity of coefficients\n";
    throw TerminalException{1};
  }
#endif
  BasisReduction<Tint> eRedBas=ComputeBasisReduction(TheBasis);
  MyMatrix<Tint> TheBasisReord=eRedBas.TheBasisReord;
  MyMatrix<T> Pmat_T=ConvertMatrixUniversal<T,Tint>(eRedBas.Pmat);
  MyMatrix<T> Pinv=Pmat_T.inverse().eval();
  MyMatrix<Tint> Pinv_int=ConvertMatrixUniversal<Tint,T>(Pinv);
  MyMatrix<Tint> Pinv_cgr=Pinv_int.transpose();
  MyMatrix<T> eSymmMatB=Pinv*eSymmMat*Pinv.transpose();
#ifdef DEBUG_UNDER_POS_COND
  std::cerr << "eRedBas.Pmat=\n";
  WriteMatrix(std::cerr, eRedBas.Pmat);
  std::cerr << "eSymmMat=\n";
  WriteMatrix(std::cerr, eSymmMat);
  std::cerr << "TheBasis=\n";
  WriteMatrix(std::cerr, TheBasis);
  //
  std::cerr << "eSymmMatB=\n";
  WriteMatrix(std::cerr, eSymmMatB);
  std::cerr << "TheBasisReord=\n";
  WriteMatrix(std::cerr, TheBasisReord);
#endif
  MyMatrix<T> TheBasisReord_T=ConvertMatrixUniversal<T,Tint>(TheBasisReord);
  MyMatrix<T> eSymmMatC=TheBasisReord_T*eSymmMatB*TheBasisReord_T.transpose();
#ifdef DEBUG_UNDER_POS_COND
  bool test2=TestCopositivityByPositivityCoeff(eSymmMatC);
  if (!test2) {
    std::cerr << "Inconsistency in the positivity of coefficients\n";
    std::cerr << "and a matrix error as well (We should have test1=test2)\n";
    throw TerminalException{1};
  }
#endif
  int n=eSymmMat.cols();
  std::vector<T> AII(n);
  for (int i=0; i<n; i++)
    AII[i]=TheBasisReord(i,i);
#ifdef DEBUG_UNDER_POS_COND
  std::cerr << "AII=";
  for (int i=0; i<n; i++)
    std::cerr << " " << AII[i];
  std::cerr << "\n";
#endif
  //
  // Now the variables whose value will change over the iteration
  //
  // x = sum_i lambda_i v_i
  // with lambda_i = h_i + z_i/a_ii  with 0 <= z_i <= z_max(i)
  //
  // We define the partial quadratic forms
  // q_i = q[sum_{j=i}^{n-1} lambda_j v_j].
  // this expression depends only on the X[j] for j>=i.
  // q_0 is the full quadratic form.
  // q_i(Z_i, Z_{i+1}, Z_{n-1}) = q_ii Z_i^2 + q_{i+1}(Z_{i+1}, Z_{i+2}, ..., Z_{n-1})
  //        + W_{i} Z_i
  //

  // Partial norms of vectors
  std::vector<T> ListNorm(n);
  // The sought after vector
  MyVector<Tint> X(n);
  // Xmax values
  std::vector<Tint> Zmax(n);
  // X values over which we iterate
  std::vector<Tint> Z(n, 0);
  // List of H: smallest values of lambda
  std::vector<T> ListH(n);
  // values of the Lambda_i:
  std::vector<T> Lambda(n);
  // PartialSum is defined as
  // PartialSum[i][j]
  MyMatrix<T> PartialSum=ZeroMatrix<T>(n+1,n+1);
  auto ComputeHXZmaxNorm=[&](int const& nupdt)-> void {
    // First we can compute the ListH just from the X
#ifdef DEBUG_UNDER_POS_COND
    std::cerr << "ComputeHXZmaxNorm : Before computing ListH, Lambda and X\n";
#endif
    for (int i=nupdt-1; i>=0; i--) {
      T eNumber=0;
      for (int j=i+1; j<n; j++)
	eNumber += TheBasisReord(j,i)*Lambda[j];
      T eFr=FractionalPart(-eNumber) / AII[i];
      ListH[i]=eFr;
      T eZ=Z[i];
      Lambda[i]=eZ/AII[i] + ListH[i];
      T eX_q=eNumber + Lambda[i]*AII[i];
#ifdef DEBUG_UNDER_POS_COND
      mpz_class eDen=eX_q.get_den();
      if (eDen > 1) {
	std::cerr << "eDen should be equal to 1\n";
	std::cerr << "eZ=" << eZ << "\n";
	std::cerr << "eX_q=" << eX_q << "\n";
	throw TerminalException{1};
      }
#endif
      X(i)=UniversalTypeConversion<Tint,T>(eX_q);
    }
#ifdef DEBUG_UNDER_POS_COND
    std::cerr << "X=";
    for (int i=0; i<n; i++)
      std::cerr << " " << X(i);
    std::cerr << "\n";
    std::cerr << "Z=";
    for (int i=0; i<n; i++)
      std::cerr << " " << Z[i];
    std::cerr << "\n";
    std::cerr << "Lambda=";
    for (int i=0; i<n; i++)
      std::cerr << " " << Lambda[i];
    std::cerr << "\n";
    std::cerr << "ListH=";
    for (int i=0; i<n; i++)
      std::cerr << " " << ListH[i];
    std::cerr << "\n";
    std::cerr << "Computing norms and Zmax\n";
    std::cerr << "n=" << n << "\n";
#endif
    for (int i=nupdt-1; i>=0; i--) {
#ifdef DEBUG_UNDER_POS_COND
      std::cerr << "i=" << i << " / " << nupdt << "\n";
#endif
      T Qii=eSymmMatC(i,i);
      T eNorm=Qii*Lambda[i]*Lambda[i];
      T eW=0;
#ifdef DEBUG_UNDER_POS_COND
      std::cerr << "Norm compute step 1\n";
#endif
      for (int j=i+1; j<n; j++)
	eW += 2*eSymmMatC(i,j)*Lambda[j];
#ifdef DEBUG_UNDER_POS_COND
      std::cerr << "Norm compute step 2\n";
#endif
      T NextNorm;
      if (i<n-1)
	NextNorm=ListNorm[i+1];
      else
	NextNorm=0;
      eNorm += eW*Lambda[i] + NextNorm;
      ListNorm[i]=eNorm;
      T eH=ListH[i];
      T eVal0=NextNorm + eW*eH + Qii*eH*eH;
      if (eVal0 >= MaxNorm) {
	Zmax[i] = 0;
      }
      else {
	T eVal1=(eW + 2*eH*Qii)/AII[i];
	T eVal2=Qii/(AII[i]*AII[i]);
	T C=MaxNorm - eVal0;
	Zmax[i] = FindLargest(eVal1, eVal2, C);
      }
    }
#ifdef DEBUG_UNDER_POS_COND
    std::cerr << "Zmax=";
    for (int i=0; i<n; i++)
      std::cerr << " " << Zmax[i];
    std::cerr << "\n";
    std::cerr << "ListNorm=";
    for (int i=0; i<n; i++)
      std::cerr << " " << ListNorm[i];
    std::cerr << "\n";
#endif
  };
  //
  // The initialization
  //
  ComputeHXZmaxNorm(n);
#ifdef DEBUG_UNDER_POS_COND
  std::cerr << "After ComputeHXZmaxNorm\n";
#endif
  //
  // The tree search function
  //
  auto NextInTree=[&]() -> bool {
#ifdef DEBUG_UNDER_POS_COND
    for (int i=0; i<n; i++)
      std::cerr << "i=" << i << " Zmax=" << Zmax[i] << "\n";
#endif
    for (int i=0; i<n; i++) {
      if (Z[i] < Zmax[i]) {
#ifdef DEBUG_UNDER_POS_COND
	std::cerr << "Clear update here Z[i]=" << Z[i] << "\n";
#endif
	Z[i]++;
#ifdef DEBUG_UNDER_POS_COND
	std::cerr << "i=" << i << " Z[i]=" << Z[i] << "\n";
#endif
	for (int j=0; j<i; j++)
	  Z[j]=0;
	ComputeHXZmaxNorm(i+1);
	return true;
      }
    }
    return false;
  };
  std::vector<MyVector<Tint>> ListVectRet;
  while(true) {
    if (ListNorm[0] <= MaxNorm) {
#ifdef DEBUG_UNDER_POS_COND
      std::cerr << "Before computing NewX\n";
      std::cerr << "Pinv_int=\n";
      WriteMatrix(std::cerr, Pinv_int);
#endif
      MyVector<Tint> NewX=Pinv_cgr*X;
#ifdef DEBUG_UNDER_POS_COND
      T eVal=EvaluationQuadForm(eSymmMat, NewX);
      if (eVal != ListNorm[0]) {
	std::cerr << "Norm inconsistency in the code\n";
	std::cerr << "eVal=" << eVal << " ListNorm[0]=" << ListNorm[0] << "\n";
	std::cerr << "NewX=";
	WriteVector(std::cerr, NewX);
	std::cerr << "X=";
	WriteVector(std::cerr, X);
	std::cerr << "\n";
	for (int i=0; i<n; i++) {
	  std::cerr << "i=" << i << " norm=" << ListNorm[i] << "\n";
	}
	throw TerminalException{1};
      }
      std::cerr << "After computing NewX. NewX=\n";
      WriteVector(std::cerr, NewX);
      if (NewX.minCoeff() < 0) {
	std::cerr << "X=";
	WriteVector(std::cerr, X);
	std::cerr << "Error. We found a negative vector\n";
	std::cerr << "eSymmMat=\n";
	WriteMatrix(std::cerr, eSymmMat);
	std::cerr << "TheBasis=\n";
	WriteMatrix(std::cerr, TheBasis);
	std::cerr << "MaxNorm=" << MaxNorm << "\n";
	throw TerminalException{1};
      }
#endif
      ListVectRet.push_back(NewX);
    }
#ifdef DEBUG_UNDER_POS_COND
    std::cerr << "X=[";
    for (int i=0; i<n; i++)
      std::cerr << " " << X(i);
    std::cerr << "]\n";
    std::cerr << "Z=[";
    for (int i=0; i<n; i++)
      std::cerr << " " << Z[i];
    std::cerr << "]\n";
    std::cerr << "Lambda=[";
    for (int i=0; i<n; i++)
      std::cerr << " " << Lambda[i];
    std::cerr << "]\n";
#endif
    bool test=NextInTree();
    if (!test)
      break;
  }
  return ListVectRet;
}

template<typename Tint>
struct SingleTestResult {
  bool test;
  std::string strNature;
  MyVector<Tint> eVectResult1;
};


template<typename Tint>
struct CopositivityEnumResult {
  bool test;
  int nbCone;
  std::vector<MyMatrix<Tint>> ListBasis;
  std::vector<MyVector<Tint>> TotalListVect;
  SingleTestResult<Tint> eResult;
};


template<typename T, typename Tint>
std::vector<MyMatrix<Tint>> PairDecomposition(MyMatrix<T> const& eSymmMat, MyMatrix<Tint> const& TheBasis)
{
  int n=eSymmMat.rows();
  MyMatrix<T> TheBasis_T=ConvertMatrixUniversal<T,Tint>(TheBasis);
  MyMatrix<T> eSymmMatB=TheBasis_T*eSymmMat*TheBasis_T.transpose();
  std::pair<Tint,Tint> PairXY_found{-1,-1};
  std::pair<int,int> PairIJ_found{-1,-1};
  // We form a matrix 
  // | a b |
  // | b c |
  // and we look for an integer vector (x,y) with x >= 0 and y >= 0 such that
  // a x + b y >= 0 and a b + c y >= 0
  // In practice it means that when we have a negative value
  auto GetSplit=[](T const& a, T const& b, T const& c) -> std::pair<Tint,Tint> {
    int eSum=2;
    while(true) {
      for (int x=1; x<eSum; x++) {
	int y=eSum-x;
	T val1 = a * x + b * y;
	T val2 = b * x + c * y;
	if (val1 >= 0 && val2 >= 0)
	  return {Tint(x),Tint(y)};
      }
      eSum++;
    }
  };
  int OptMinimum=1;
  bool IsFirst=true;
  Tint MinSumXY=0;
  T MaxValQuot=0;
  for (int i=0; i<n-1; i++)
    for (int j=i+1; j<n; j++) {
      T a=eSymmMatB(i,i);
      T b=eSymmMatB(i,j);
      T c=eSymmMatB(j,j);
      if (b < 0) {
	std::pair<Tint,Tint> CandXY;
	Tint sumCoef=0;
	T eQuot=0;
	if (OptMinimum == 1) {
	  CandXY = {1, 1};
	  eQuot = (b*b) / (a*c);
	}
	if (OptMinimum == 2) {
	  CandXY = GetSplit(a, b, c);
	  sumCoef = CandXY.first + CandXY.second;
	}
	auto InsertValue=[&]() -> bool {
	  if (IsFirst) {
	    IsFirst=false;
	    MaxValQuot=eQuot;
	    MinSumXY=sumCoef;
	    return true;
	  }
	  if (OptMinimum == 1 && eQuot > MaxValQuot) {
	    MaxValQuot = eQuot;
	    return true;
	  }
	  if (OptMinimum == 2 && sumCoef < MinSumXY) {
	    MinSumXY=sumCoef;
	    return true;
	  }
	  return false;
	};
	bool test=InsertValue();
	if (test) {
	  PairXY_found = CandXY;
	  PairIJ_found = {i, j};
	}
      }
    }
  int iFound = PairIJ_found.first;
  int jFound = PairIJ_found.second;
  if (iFound == -1) {
    std::cerr << "iFound=" << iFound << " jFound=" << jFound << " x=" << PairXY_found.first << " y=" << PairXY_found.second << "\n";
    std::cerr << "eSymmMatB=\n";
    WriteMatrix(std::cerr, eSymmMatB);
    std::cerr << "TheBasis=\n";
    WriteMatrix(std::cerr, TheBasis);
    std::cerr << "Looks like the matrix is all non-negative. No need for splitting. Most likely a bug\n";
    throw TerminalException{1};
  }
#ifdef PRINT_SPLIT_CONE
  std::cerr << "iFound=" << iFound << " jFound=" << jFound << " PairXY_found=" << PairXY_found.first << " / " << PairXY_found.second << "\n";
#endif
  
  MyVector<Tint> eVectMid(n);
  for (int i=0; i<n; i++)
    eVectMid(i)=PairXY_found.first * TheBasis(iFound, i) + PairXY_found.second * TheBasis(jFound, i);
  MyMatrix<Tint> TheBasis1(n,n);
  MyMatrix<Tint> TheBasis2(n,n);
  for (int i=0; i<n; i++) {
    MyVector<Tint> eVect1;
    MyVector<Tint> eVect2;
    if (i == iFound) {
      eVect1=eVectMid;
    }
    else {
      eVect1=TheBasis.row(i);
    }
    //
    if (i == jFound) {
      eVect2=eVectMid;
    }
    else {
      eVect2=TheBasis.row(i);
    }
    TheBasis1.row(i) = eVect1;
    TheBasis2.row(i) = eVect2;
  }
  MyMatrix<T> TheBasis1_T=ConvertMatrixUniversal<T,Tint>(TheBasis1);
  MyMatrix<T> TheBasis2_T=ConvertMatrixUniversal<T,Tint>(TheBasis2);
  MyMatrix<T> eSymmMatB1=TheBasis1_T*eSymmMat*TheBasis1_T.transpose();
  MyMatrix<T> eSymmMatB2=TheBasis2_T*eSymmMat*TheBasis2_T.transpose();
  if (OptMinimum == 2) {
    if (eSymmMatB1(iFound, jFound) < 0 || eSymmMatB2(iFound, jFound) < 0) {
      std::cerr << "We clearly have a bug in our computation\n";
      throw TerminalException{1};
    }
  }
  return {TheBasis1, TheBasis2};
}




template<typename T, typename Tint>
SingleTestResult<Tint> SingleTestStrictCopositivity(MyMatrix<T> const& eSymmMat, MyMatrix<Tint> const& TheBasis, MyMatrix<T> const& eSymmMatB)
{
  int n=eSymmMat.rows();
  MyVector<Tint> eVectZero;
  for (int i=0; i<n; i++)
    if (eSymmMatB(i,i) <= 0)
      return {false, "Not Positive diag", TheBasis.row(i)};
  for (int i=0; i<n; i++)
    for (int j=i+1; j<n; j++) {
      T a=eSymmMatB(i,i);
      T b=eSymmMatB(i,j);
      T c=eSymmMatB(j,j);
      if (b < 0) {
	if (b*b == a*c) {
	  MyVector<T> V(2);
	  V(0) = -b;
	  V(1) = a;
	  MyVector<T> Vred = RemoveFractionVector(V);
	  MyVector<Tint> Vint = ConvertVectorUniversal<Tint,T>(Vred);
	  MyVector<Tint> eVect = Vint(0) * TheBasis.row(i) + Vint(1) * TheBasis.row(j);
	  return {false, "zero vector detected", eVect};
	}
	if (b*b > a*c) {
	  int eSum=2;
	  while(true) {
	    for (int x=1; x<eSum; x++) {
	      int y=eSum - x;
	      T qVal = a * x*x + 2*b * x*y + c * y*y;
	      if (qVal < 0) {
		MyVector<Tint> eVect = x * TheBasis.row(i) + y * TheBasis.row(j);
		return {false, "Off diagonal violation", eVect};
	      }
	    }
	    eSum++;
	  }
	}
      }
    }
  return {true, "on the surface ok", eVectZero};
}


template<typename T, typename Tint>
SingleTestResult<Tint> SingleTestCopositivity(MyMatrix<T> const& eSymmMat, MyMatrix<Tint> const& TheBasis, MyMatrix<T> const& eSymmMatB)
{
  int n=eSymmMat.rows();
  MyVector<Tint> eVectZero;  
  for (int i=0; i<n; i++)
    if (eSymmMatB(i,i) < 0)
      return {false, "Not Positive diag", TheBasis.row(i)};
  for (int i=0; i<n; i++)
    for (int j=i+1; j<n; j++) {
      T a=eSymmMatB(i,i);
      T b=eSymmMatB(i,j);
      T c=eSymmMatB(j,j);
      if (b < 0 && b*b > a*c) {
	int eSum=2;
	while(true) {
	  for (int x=1; x<eSum; x++) {
	    int y=eSum - x;
	    T qVal = a * x*x + 2*b * x*y + c * y*y;
	    if (qVal < 0) {
	      MyVector<Tint> eVect = x * TheBasis.row(i) + y * TheBasis.row(j);
	      return {false, "Off diagonal violation", eVect};
	    }
	  }
	  eSum++;
	}
      }
    }
  return {true, "on the surface ok", eVectZero};
}





template<typename T, typename Tint>
SingleTestResult<Tint> SearchByZeroInKernel(MyMatrix<T> const& eSymmMat)
{
  MyVector<Tint> eVectZero;  
  MyMatrix<T> NSP = NullspaceMat(eSymmMat);
  int nbCol=NSP.cols();
  int nbNSP=NSP.rows();
  if (nbNSP == 0)
    return {true, "on the surface ok", eVectZero};
  MyMatrix<T> FAC(nbCol+1, nbNSP+1);
  for (int iCol=0; iCol<nbCol; iCol++) {
    FAC(iCol, 0) = 0;
    for (int iNSP=0; iNSP<nbNSP; iNSP++)
      FAC(iCol, iNSP+1) = NSP(iNSP, iCol);
  }
  MyVector<T> ToBeMinimized(nbNSP+1);
  ToBeMinimized(0) = 0;
  FAC(nbCol, 0) = -1;
  for (int iNSP=0; iNSP<nbNSP; iNSP++) {
    T eSum=0;
    for (int iCol=0; iCol<nbCol; iCol++)
      eSum += NSP(iNSP, iCol);
    ToBeMinimized(iNSP+1) = eSum;
    FAC(nbCol, iNSP+1) = eSum;
  }
  LpSolution<T> eSol=CDD_LinearProgramming(FAC, ToBeMinimized);
  if (!eSol.PrimalDefined)
    return {true, "on the surface ok", eVectZero};
  MyVector<T> eVect1(nbCol);
  for (int iCol=0; iCol<nbCol; iCol++) {
    T eSum=0;
    for (int iNSP=0; iNSP<nbNSP; iNSP++)
      eSum += eSol.DirectSolution(iNSP) * NSP(iNSP,iCol);
    eVect1(iCol) = eSum;
  }
  MyVector<T> eVect2=RemoveFractionVector(eVect1);
  MyVector<Tint> eVect3=ConvertVectorUniversal<Tint,T>(eVect2);
  return {false, "Zero vector from kernel", eVect3};
}



template<typename T>
struct RequestCopositivity {
  T MaxNorm;
  bool DoListCone;
};


template<typename T, typename Tint>
CopositivityEnumResult<Tint> KernelEnumerateShortVectorInCone(MyMatrix<T> const& eSymmMat, MyMatrix<Tint> const& TheBasis, RequestCopositivity<T> const& CopoReq)
{
  std::vector<MyVector<Tint>> TotalList;
  //
  // First clear up with all signs being positive case
  //
  MyMatrix<T> TheBasis_T=ConvertMatrixUniversal<T,Tint>(TheBasis);
  MyMatrix<T> eSymmMatB=TheBasis_T*eSymmMat*TheBasis_T.transpose();
  bool test=TestCopositivityByPositivityCoeff(eSymmMatB);
  if (test) {
    TotalList = EnumerateShortVectorInCone_UnderPositivityCond<T,Tint>(eSymmMat, TheBasis, CopoReq.MaxNorm);
    if (!CopoReq.DoListCone)
      return {true, 1, {}, TotalList, {}};
    else
      return {true, 1, {TheBasis}, TotalList, {}};
  }
  //
  // Second see if the vectors allow to conclude directly.
  //
  SingleTestResult<Tint> eResult=SingleTestStrictCopositivity(eSymmMat, TheBasis, eSymmMatB);
  if (!eResult.test)
    return {false, 0, {}, {}, eResult};
  //
  // Third, split the domain.
  //
  std::vector<MyMatrix<Tint>> ListBasis=PairDecomposition(eSymmMat, TheBasis);
  int nbCone=0;
  std::vector<MyMatrix<Tint>> ListBasisInt;
  for (auto & eBasis : ListBasis) {
    CopositivityEnumResult<Tint> fEnumResult=KernelEnumerateShortVectorInCone(eSymmMat, eBasis, CopoReq);
    nbCone += fEnumResult.nbCone;
    if (CopoReq.DoListCone)
      ListBasisInt.insert(ListBasisInt.end(), fEnumResult.ListBasis.begin(), fEnumResult.ListBasis.end());
    if (!fEnumResult.test)
      return {false, nbCone, {}, {}, fEnumResult.eResult};
    TotalList.insert(TotalList.end(), fEnumResult.TotalListVect.begin(), fEnumResult.TotalListVect.end());
  }
  return {true, nbCone, ListBasisInt, TotalList, {}};
}







// This is supposed to be a tree version and so faster than KernelEnumerateShortVectorInCone
// The other program is recursive and so has a problem with too large level.
template<typename T, typename Tint>
SingleTestResult<Tint> EnumerateCopositiveShortVector_Kernel(MyMatrix<T> const& eSymmMat, std::function<bool(MyMatrix<Tint> const&,MyMatrix<T> const&)> const& fInsert, std::function<SingleTestResult<Tint>(MyMatrix<Tint> const&, MyMatrix<T> const&)> const& fSingleTest)
{
  std::cerr << "Begin of EnumerateCopositiveShortVector_Kernel\n";
  int n=eSymmMat.rows();
  struct DataPair {
    int idx;
    MyMatrix<Tint> TheBasis0;
    MyMatrix<Tint> TheBasis1;
  };
  std::vector<DataPair> ListPair;
  auto SetValue=[&](int const& pos, DataPair const& ePair) -> void {
    int siz=ListPair.size();
    if (pos == siz) {
      ListPair.push_back(ePair);
    }
    else {
      ListPair[pos] = ePair;
    }
  };
  //
  SingleTestResult<Tint> eResult{true, "all ok", {}};
  //
  int position=-1;
  auto GetBasis=[&]() -> MyMatrix<Tint> {
    if (position == -1)
      return IdentityMat<Tint>(n);
    int idx=ListPair[position].idx;
    if (idx == 0)
      return ListPair[position].TheBasis0;
    return ListPair[position].TheBasis1;
  };
  auto GoUpNextInTree=[&]() -> bool {
    while(true) {
      if (position == -1)
	return false;
      int idx=ListPair[position].idx;
      if (idx == 0) {
	ListPair[position].idx = 1;
	return true;
      }
      position--;
    }
  };
  auto GetDataPair=[&](std::vector<MyMatrix<Tint>> const& ListBasis) -> DataPair {
    return {0, ListBasis[0], ListBasis[1]};
  };
  auto NextInTree=[&]() -> bool {
    MyMatrix<Tint> TheBasis = GetBasis();
    //
    // If we have vi A vj >= 0 then direct approach works.
    //
    MyMatrix<T> TheBasis_T=ConvertMatrixUniversal<T,Tint>(TheBasis);
    MyMatrix<T> eSymmMatB=TheBasis_T*eSymmMat*TheBasis_T.transpose();
    bool testB=fInsert(TheBasis, eSymmMatB);
    if (testB) {
      return GoUpNextInTree();
    }
    //
    // If we have some negative values then we exit by failure
    //
    SingleTestResult<Tint> fResult=fSingleTest(TheBasis, eSymmMatB);
    if (!fResult.test) {
      eResult = fResult;
      return false;
    }
    //
    // Now doing the split
    //
    std::vector<MyMatrix<Tint>> ListBasis=PairDecomposition(eSymmMat, TheBasis);
    position++;
    SetValue(position, GetDataPair(ListBasis));
    return true;
  };
  while(true) {
    bool res = NextInTree();
    if (!res)
      break;
  }
  return eResult;
}


template<typename T, typename Tint>
Tshortest<T,Tint> T_CopositiveShortestVector(MyMatrix<T> const& eSymmMat)
{
  SingleTestResult<Tint> kerResult = SearchByZeroInKernel<T,Tint>(eSymmMat);
  if (!kerResult.test) {
    std::cerr << "Inconsistency in the run\n";
    throw TerminalException{1};
  }
  int n=eSymmMat.rows();
  T MinNorm=MinimumDiagonal(eSymmMat);
  int nbCone=0;
  std::set<MyVector<Tint>> TotalListVect_set;
  std::function<bool(MyMatrix<Tint> const&, MyMatrix<T> const&)> fInsert=[&](MyMatrix<Tint> const& TheBasis, MyMatrix<T> const& eSymmMatB) -> bool {
    bool test=TestCopositivityByPositivityCoeff(eSymmMatB);
    if (!test)
      return false;
    std::vector<MyVector<Tint>> TotalList=EnumerateShortVectorInCone_UnderPositivityCond(eSymmMat, TheBasis, MinNorm);
    for (auto & eVect : TotalList) {
      T eNorm=0;
      for (int i=0; i<n; i++) {
	for (int j=0; j<n; j++) {
	  T val1=eVect(i);
	  T val2=eVect(j);
	  eNorm += eSymmMat(i,j) * val1 * val2;
	}
      }
      if (eNorm > 0) {
	if (eNorm < MinNorm) {
	  MinNorm=eNorm;
	  TotalListVect_set.clear();
	  TotalListVect_set.insert(eVect);
	}
	else {
	  if (eNorm == MinNorm)
	    TotalListVect_set.insert(eVect);
	}
      }
    }
    nbCone++;
    return true;
  };
  std::function<SingleTestResult<Tint>(MyMatrix<Tint> const&, MyMatrix<T> const&)> fSingleTest=[&](MyMatrix<Tint> const& TheBasis, MyMatrix<T> const& eSymmMatB) -> SingleTestResult<Tint> {
    return SingleTestStrictCopositivity(eSymmMat, TheBasis, eSymmMatB);
  };
  SingleTestResult<Tint> eResult = EnumerateCopositiveShortVector_Kernel(eSymmMat, fInsert, fSingleTest);
#ifdef PRINT_NBCONE
  std::cerr << "nbCone=" << nbCone << "\n";
#endif
  if (!eResult.test) {
    std::cerr << "We should not have non-copositivity in T_CopositiveShortestVector\n";
    throw TerminalException{1};
  }
  int nbVect=TotalListVect_set.size();
  int iVect=0;
  MyMatrix<Tint> SHV(nbVect, n);
  for (auto & eVect : TotalListVect_set) {
    AssignMatrixRow(SHV, iVect, eVect);
    iVect++;
  }
  return {MinNorm, SHV};
}




template<typename T, typename Tint>
CopositivityEnumResult<Tint> EnumerateCopositiveShortVector_V2(MyMatrix<T> const& eSymmMat, RequestCopositivity<T> const& CopoReq)
{
  SingleTestResult<Tint> kerResult = SearchByZeroInKernel<T,Tint>(eSymmMat);
  if (!kerResult.test) {
    return {false, 0, {}, {}, kerResult};
  }
  int nbCone=0;
  std::vector<MyMatrix<Tint>> ListBasis;
  std::set<MyVector<Tint>> TotalListVect_set;
  std::function<bool(MyMatrix<Tint> const&, MyMatrix<T> const&)> fInsert=[&](MyMatrix<Tint> const& TheBasis, MyMatrix<T> const& eSymmMatB) -> bool {
    bool test=TestCopositivityByPositivityCoeff(eSymmMatB);
    if (!test)
      return false;
    std::vector<MyVector<Tint>> TotalList=EnumerateShortVectorInCone_UnderPositivityCond(eSymmMat, TheBasis, CopoReq.MaxNorm);
    if (CopoReq.DoListCone)
      ListBasis.push_back(TheBasis);
    for (auto & eVect : TotalList)
      TotalListVect_set.insert(eVect);
    nbCone++;
    return true;
  };
  std::function<SingleTestResult<Tint>(MyMatrix<Tint> const&, MyMatrix<T> const&)> fSingleTest=[&](MyMatrix<Tint> const& TheBasis, MyMatrix<T> const& eSymmMatB) -> SingleTestResult<Tint> {
    return SingleTestStrictCopositivity(eSymmMat, TheBasis, eSymmMatB);
  };
  SingleTestResult<Tint> eResult = EnumerateCopositiveShortVector_Kernel(eSymmMat, fInsert, fSingleTest);
  if (!eResult.test) {
    return {false, 0, {}, {}, eResult};
  }
  std::vector<MyVector<Tint>> TotalListVect;
  for (auto & eVect : TotalListVect_set)
    TotalListVect.push_back(eVect);
  return {true, nbCone, ListBasis, TotalListVect, eResult};
}



template<typename T, typename Tint>
SingleTestResult<Tint> TestCopositivity(MyMatrix<T> const& eSymmMat)
{
  int n=eSymmMat.rows();
  int nbCone=0;
  std::function<bool(MyMatrix<int> const&, MyMatrix<T> const&)> fInsert=[&](MyMatrix<int> const& TheBasis, MyMatrix<T> const& eSymmMatB) -> bool {
    for (int i=0; i<n; i++)
      for (int j=i+1; j<n; j++)
	if (eSymmMatB(i,j) < 0)
	  return false;
    nbCone++;
    return true;
  };
  std::function<SingleTestResult<Tint>(MyMatrix<int> const&, MyMatrix<T> const&)> fSingleTest=[&](MyMatrix<int> const& TheBasis, MyMatrix<T> const& eSymmMatB) -> SingleTestResult<Tint> {
    return SingleTestCopositivity(eSymmMat, TheBasis, eSymmMatB);
  };
  SingleTestResult<Tint> eResult = EnumerateCopositiveShortVector_Kernel(eSymmMat, fInsert, fSingleTest);
  std::cerr << "nbCone=" << nbCone << "\n";
  return eResult;
}







template<typename T, typename Tint>
CopositivityEnumResult<Tint> EnumerateCopositiveShortVector_V1(MyMatrix<T> const& eSymmMat, RequestCopositivity<T> const& CopoReq)
{
  int n=eSymmMat.rows();
  MyMatrix<int> TheBasis=IdentityMat<int>(n);
  CopositivityEnumResult<Tint> CopoRes=KernelEnumerateShortVectorInCone(eSymmMat, TheBasis, CopoReq);
  std::set<std::vector<int>> TheSet;
  // We have to do something non-clever to remove duplicates
  for (auto & eVect : CopoRes.TotalListVect) {
    std::vector<int> eVectB(n);
    for (int i=0; i<n; i++)
      eVectB[i]=eVect(i);
    TheSet.insert(eVectB);
  }
  std::vector<MyVector<int>> TotalListRed;
  for (auto & eVect : TheSet) {
    MyVector<int> eVectB(n);
    for (int i=0; i<n; i++)
      eVectB(i)=eVect[i];
    TotalListRed.push_back(eVectB);
  }
  return {CopoRes.test, CopoRes.nbCone, {}, TotalListRed, CopoRes.eResult};
}


template<typename T, typename Tint>
CopositivityEnumResult<Tint> EnumerateCopositiveShortVector(MyMatrix<T> const& eSymmMat, RequestCopositivity<T> const& CopoReq)
{
//#define FULL_DEBUG
#ifdef FULL_DEBUG
  CopositivityEnumResult<Tint> res1 = EnumerateCopositiveShortVector_V1(eSymmMat, CopoReq);
  CopositivityEnumResult<Tint> res2 = EnumerateCopositiveShortVector_V2(eSymmMat, CopoReq);
  if (res1.test != res2.test) {
    std::cerr << "Inconsistency in computing res1.test and res2.test\n";
    throw TerminalException{1};
  }
  if (!res1.test)
    return res1;
  if (res1.TotalListVect != res2.TotalListVect) {
    std::cerr << "res1.TotalListVect, |res1.TotalListVect|=" << res1.TotalListVect.size() << "\n";
    for (auto & eVect : res1.TotalListVect)
      WriteVector(std::cerr, eVect);
    //
    std::cerr << "res2.TotalListVect, |res2.TotalListVect|=" << res2.TotalListVect.size() << "\n";
    for (auto & eVect : res2.TotalListVect)
      WriteVector(std::cerr, eVect);
    //
    throw TerminalException{1};
  }
  return res1;
#else
  return EnumerateCopositiveShortVector_V2<T,Tint>(eSymmMat, CopoReq);
#endif
}



template<typename T, typename Tint>
Tshortest<T,Tint> T_CopositiveShortestVector_V1(MyMatrix<T> const& eSymmMat)
{
  int n=eSymmMat.rows();
  T MaxNorm=MinimumDiagonal(eSymmMat);
  RequestCopositivity<T> CopoReq{MaxNorm, false};
  CopositivityEnumResult<Tint> CopoRes=EnumerateCopositiveShortVector(eSymmMat, CopoReq);
  if (!CopoRes.test) {
    std::cerr << "Inconsistency in result\n";
    throw TerminalException{1};
  }
  int nbVect=CopoRes.TotalListVect.size();
  MyMatrix<int> SHV(nbVect, n);
  for (int iVect=0; iVect<nbVect; iVect++) {
    MyVector<Tint> V = CopoRes.TotalListVect[iVect];
    AssignMatrixRow(SHV, iVect, V);
  }
  return SelectShortestVector<T,Tint>(eSymmMat, SHV);
}





#endif

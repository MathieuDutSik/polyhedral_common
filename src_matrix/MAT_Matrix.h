#ifndef TEMP_MATRIX_INCLUDE
#define TEMP_MATRIX_INCLUDE

// For reference over Eigen, see for example
// http://eigen.tuxfamily.org/dox/AsciiQuickReference.txt
//

#include "Temp_common.h"
#include "Basic_string.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>

template <typename T>
using MyVector = Eigen::Matrix<T,Eigen::Dynamic,1>;

template <typename T>
using MyMatrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;

template <typename T>
using MySparseMatrix = Eigen::SparseMatrix<T,Eigen::ColMajor>;


//
// Matrix type conversion
//

template<typename T1, typename T2>
MyVector<T2> ConvertVector(MyVector<T1> const& V, std::function<T2(T1 const&)> const& f)
{
  int n=V.size();
  MyVector<T2> eRet(n);
  for (int i=0; i<n; i++) {
    T1 eVal1=V(i);
    T2 eVal2=f(eVal1);
    eRet(i)=eVal2;
  }
  return eRet;
}

template<typename T2, typename T1>
MyVector<T2> ConvertVectorUniversal(MyVector<T1> const& V)
{
  int n=V.size();
  MyVector<T2> eRet(n);
  for (int i=0; i<n; i++) {
    T1 eVal1=V(i);
    T2 eVal2=UniversalTypeConversion<T2,T1>(eVal1);
    eRet(i)=eVal2;
  }
  return eRet;
}





template<typename T1, typename T2>
MyMatrix<T2> ConvertMatrix(MyMatrix<T1> const& M, std::function<T2(const T1&)> const& f)
{
  int eta_rho=M.rows();
  int xi_rho=M.cols();
  MyMatrix<T2> eRet(eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      T1 eVal1=M(i,j);
      T2 eVal2=f(eVal1);
      eRet(i,j)=eVal2;
  }
  return eRet;
}

template<typename T2, typename T1>
MyMatrix<T2> ConvertMatrixUniversal(MyMatrix<T1> const& M)
{
  int eta_rho=M.rows();
  int xi_rho=M.cols();
  MyMatrix<T2> eRet(eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      T1 eVal1=M(i,j);
      T2 eVal2=UniversalTypeConversion<T2,T1>(eVal1);
      eRet(i,j)=eVal2;
    }
  return eRet;
}


template<typename T>
MyVector<T> VectorWithIdenticalEntries(int const& len, T const& val)
{
  MyVector<T> V(len);
  for (int i=0; i<len; i++)
    V(i) = val;
  return V;
}



template<typename T1, typename T2>
MySparseMatrix<T2> ConvertSparseMatrix(MySparseMatrix<T1> const& M, std::function<T2(T1)> const& f)
{
  int nbRow=M.rows();
  int nbCol=M.cols();
  int nnz=M.nonZeros();
  using Ttrip = Eigen::Triplet<T2>;
  std::vector<Ttrip> tripletList(nnz);
  int nb=0;
  for (int k=0; k<M.outerSize(); ++k)
    for (typename MySparseMatrix<T1>::InnerIterator it(M,k); it; ++it) {
      T1 eVal1=it.value();
      T2 eVal2=f(eVal1);
      int iRow=it.row();   // row index
      int iCol=it.col();   // col index (here it is equal to k)
      tripletList[nb]=Ttrip(iRow,iCol,eVal2);
      nb++;
    }
  MySparseMatrix<T2> SpMat=MySparseMatrix<T2>(nbRow, nbCol);
  SpMat.setFromTriplets(tripletList.begin(), tripletList.end());
  return SpMat;
}

//
// ReadWriteMatrices and vectors
//



template<typename T>
MyMatrix<T> ReadMatrix(std::istream &is)
{
  if (!is.good()) {
    std::cerr << "ReadMatrix operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  T eVal;
  int nbRow, nbCol;
  int iRow, iCol;
  is >> nbRow >> nbCol;
  std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
  if (nbRow < 0 || nbCol < 0) {
    std::cerr << "We should have nbRow > 0 and nbCol > 0\n";
    std::cerr << "But we have nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    throw TerminalException{1};
  }
  //std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
  MyMatrix<T> TheMat(nbRow, nbCol);
  for (iRow=0; iRow<nbRow; iRow++)
    for (iCol=0; iCol<nbCol; iCol++) {
      is >> eVal;
      TheMat(iRow, iCol)=eVal;
    }
  return TheMat;
}



template<typename T>
MyMatrix<T> ReadMatrixLrsCdd(std::istream &is)
{
  if (!is.good()) {
    std::cerr << "ReadMatrixLrs operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  std::string header;
  is >> header;
  if (header != "V-representation" && header != "H-representation") {
    std::cerr << "Error while reading header\n";
    std::cerr << "header = " << header << "\n";
    throw TerminalException{1};
  }
  //
  std::string head2;
  is >> head2;
  if (head2 != "begin") {
    std::cerr << "Error while reading\n";
    std::cerr << "head2 = " << head2 << "\n";
    throw TerminalException{1};
  }
  //
  int nbRow, nbCol;
  std::string str;
  is >> nbRow >> nbCol >> str;
  if (str != "integer" && str != "rational" && str != "double") {
    std::cerr << "Error while reading\n";
    std::cerr << "str = " << str << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
  MyMatrix<T> TheMat(nbRow, nbCol);
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      T eVal;
      is >> eVal;
      TheMat(iRow, iCol)=eVal;
    }
  //
  std::string footer;
  is >> footer;
  if (footer != "end") {
    std::cerr << "Error while reading\n";
    std::cerr << "footer = " << footer << "\n";
    throw TerminalException{1};
  }
  return TheMat;
}


template<typename T>
bool IsEqualSizeMatrices(MyMatrix<T> const& X, MyMatrix<T> const& Y)
{
  if (X.rows() != Y.rows())
    return false;
  if (X.cols() != Y.cols())
    return false;
  return true;
}


template<typename T>
std::string StringSizeMatrix(MyMatrix<T> const& X)
{
  int nbRow=X.rows();
  int nbCol=X.cols();
  return IntToString(nbRow) + " / " + IntToString(nbCol);
}


template<typename T>
std::string MinMaxMatrix(MyMatrix<T> const& X)
{
  T minV=X.minCoeff();
  T maxV=X.maxCoeff();
  std::stringstream s;
  s << "min/max=" << minV << " / " << maxV;
  std::string converted(s.str());
  return converted;
}




template<typename T>
MyVector<T> ReadVector(std::istream &is)
{
  if (!is.good()) {
    std::cerr << "ReadVector operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  T eVal;
  size_t nbRow;
  is >> nbRow;
  //  std::cerr << "nbRow=" << nbRow << "\n";
  MyVector<T> eVect=MyVector<T>(nbRow);
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    is >> eVal;
    eVect(iRow)=eVal;
  }
  return eVect;
}

template<typename T>
std::vector<T> ReadStdVector(std::istream &is)
{
  if (!is.good()) {
    std::cerr << "ReadStdVector operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  T eVal;
  size_t nbRow;
  is >> nbRow;
  //  std::cerr << "nbRow=" << nbRow << "\n";
  std::vector<T> eVect(nbRow);
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    is >> eVal;
    eVect[iRow]=eVal;
  }
  return eVect;
}



template<typename T>
void WriteMatrix(std::ostream &os, MyMatrix<T> const&TheMat)
{
  int nbRow=TheMat.rows();
  int nbCol=TheMat.cols();
  //  TerminalEnding();
  os << nbRow << " " << nbCol << "\n";
  for (int iRow=0; iRow<nbRow; iRow++) {
    for (int iCol=0; iCol<nbCol; iCol++) {
      T eVal=TheMat(iRow, iCol);
      os << " " << eVal;
    }
    os << "\n";
  }
}

template<typename T>
void WriteMatrixMatlab(std::ostream &os, MyMatrix<T> const&TheMat)
{
  int nbRow=TheMat.rows();
  int nbCol=TheMat.cols();
  for (int iRow=0; iRow<nbRow; iRow++) {
    for (int iCol=0; iCol<nbCol; iCol++) {
      T eVal=TheMat(iRow, iCol);
      os << " " << eVal;
    }
    os << "\n";
  }
}



template<typename T>
void WriteSparseMatrixGAP(std::ostream &os, MySparseMatrix<T> const&TheMat)
{
  long nbRow=TheMat.rows();
  long nbCol=TheMat.cols();
  struct PairCV {
    int iCol;
    T eVal;
  };
  std::vector<std::vector<PairCV> > LLPair(nbRow);
  for (int k=0; k<TheMat.outerSize(); ++k)
    for (typename MySparseMatrix<T>::InnerIterator it(TheMat,k); it; ++it) {
      T eVal=it.value();
      int iRow=it.row();
      int iCol=it.col();
      PairCV ePair{iCol, eVal};
      LLPair[iRow].push_back(ePair);
    }
  os << "rec(nbRow:=" << nbRow << ", nbCol:=" << nbCol << ", ListEntries:=[";
  for (int iRow=0; iRow<nbRow; iRow++) {
    if (iRow > 0)
      os << ",\n";
    int len=LLPair[iRow].size();
    os << "rec(ListCol:=[";
    for (int i=0; i<len; i++) {
      if (i>0)
	os << ", ";
      int eCol=LLPair[iRow][i].iCol + 1;
      os << eCol;
    }
    os << "], ListVal:=[";
    for (int i=0; i<len; i++) {
      if (i>0)
	os << ", ";
      os << LLPair[iRow][i].eVal;
    }
    os << "])";
  }
  os << "])";
}





template<typename T>
void WriteMatrixGAP(std::ostream &os, MyMatrix<T> const&TheMat)
{
  long nbRow=TheMat.rows();
  long nbCol=TheMat.cols();
  os << "[ ";
  for (long iRow=0; iRow<nbRow; iRow++) {
    if (iRow > 0)
      os << ",\n";
    os << "[";
    for (long iCol=0; iCol<nbCol; iCol++) {
      T eVal=TheMat(iRow, iCol);
      if (iCol > 0)
	os << ",";
      os << " " << eVal;
    }
    os << "]";
  }
  os << "]";
}


template<typename T>
void WriteMatrixGAPfile(std::string const& eFile, MyMatrix<T> const&TheMat)
{
  std::ofstream os(eFile);
  os << "return ";
  WriteMatrixGAP(os, TheMat);
  os << ";\n";
}



template<typename T>
void WriteVector(std::ostream &os, MyVector<T> const & TheVec)
{
  int n=TheVec.size();
  for (int i=0; i<n; i++)
    os << " " << TheVec(i);
  os << "\n";
}

template<typename T>
void WriteVectorGAP(std::ostream &os, MyVector<T> const & TheVec)
{
  int n=TheVec.size();
  os << "[";
  for (int i=0; i<n; i++) {
    if (i>0)
      os << ", ";
    os << TheVec(i);
  }
  os << "]";
}




template<typename T>
T VectorSum(MyVector<T> const& eVect)
{
  T eSum=0;
  int siz=eVect.size();
  for (int i=0; i<siz; i++)
    eSum += eVect(i);
  return eSum;
}


template<typename T>
T ScalarProduct(MyVector<T> const& V1, MyVector<T> const & V2)
{
  if (V1.size() != V2.size()) {
    std::cerr << "Vectors of wrong sizes\n";
    throw TerminalException{1};
  }
  size_t siz=V1.size();
  T eSum=0;
  for (size_t i=0; i<siz; i++)
    eSum += V1(i)*V2(i);
  return eSum;
}

// compute the scalar product of V with the lines of eMat
template<typename T>
MyVector<T> ListScalarProduct(MyVector<T> const& V, MyMatrix<T> const& eMat)
{
  int dim=V.size();
  int nbRow=eMat.rows();
  int nbCol=eMat.cols();
  if (nbCol != dim) {
    std::cerr << "dim should equal nbCol\n";
    throw TerminalException{1};
  }
  MyVector<T> retVect(nbRow);
  for (int iRow=0; iRow<nbRow; iRow++) {
    T eSum=0;
    for (int iCol=0; iCol<dim; iCol++) {
      T eVal1=V(iCol);
      T eVal2=eMat(iRow, iCol);
      eSum += eVal1*eVal2;
    }
    retVect(iRow)=eSum;
  }
  return retVect;
}



template<typename T>
MyMatrix<T> ZeroMatrix(int const& nbRow, int const& nbCol)
{
  MyMatrix<T> retMat(nbRow, nbCol);
  T eZero;
  eZero=0;
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++)
      retMat(iRow, iCol)=eZero;
  return retMat;
}


template<typename T>
MyMatrix<T> ConstantMatrix(int const& nbRow, int const& nbCol, T const& eVal)
{
  MyMatrix<T> retMat(nbRow, nbCol);
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++)
      retMat(iRow, iCol)=eVal;
  return retMat;
}




template<typename T>
MyVector<T> ZeroVector(int const& nbRow)
{
  MyVector<T> retVect(nbRow);
  T eZero;
  eZero=0;
  for (int iRow=0; iRow<nbRow; iRow++)
    retVect(iRow)=eZero;
  return retVect;
}




template<typename T>
void TVec_ZeroAssignation(MyVector<T> &TheVect)
{
  int i;
  T eZero;
  eZero=0;
  for (i=0; i<TheVect.n; i++)
    TheVect(i)=eZero;
}




/* eMatO should already be allocated with same ssizes as eMatI */
template<typename T>
MyMatrix<T> TMat_ConvertFromInt(MyMatrix<int> const&eMatI)
{
  int nbRow=eMatI.rows();
  int nbCol=eMatI.cols();
  MyMatrix<T> eMatO(nbRow, nbCol);
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      int eVal=eMatI(iRow, iCol);
      T eVal_T=eVal;
      eMatO(iRow, iCol)=eVal_T;
    }
  return eMatO;
}






template<typename T>
void ZeroAssignation(MyMatrix<T> &TheMat)
{
  int nbRow=TheMat.rows();
  int nbCol=TheMat.cols();
  T eVal;
  eVal=0;
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++)
      TheMat(iRow, iCol)=eVal;
}

template<typename T>
MyMatrix<T> TransposedMat(MyMatrix<T> const&TheMat)
{
  int nbCol=TheMat.cols();
  int nbRow=TheMat.rows();
  MyMatrix<T> TheTrans(nbCol, nbRow);
  for (int iCol=0; iCol<nbCol; iCol++)
    for (int iRow=0; iRow<nbRow; iRow++)
      TheTrans(iCol, iRow)=TheMat(iRow, iCol);
  return TheTrans;
}


// Compute the product XM
template<typename T>
MyVector<T> ProductVectorMatrix(MyVector<T> const& X, MyMatrix<T> const& M)
{
  int nbCol=M.cols();
  int nbRow=M.rows();
  if (X.size() != nbRow) {
    std::cerr << "Error in the product X A\n";
    throw TerminalException{1};
  }
  MyVector<T> Vret(nbCol);
  for (int iCol=0; iCol<nbCol; iCol++) {
    T sum=0;
    for (int iRow=0; iRow<nbRow; iRow++)
      sum += M(iRow,iCol)*X(iRow);
    Vret(iCol)=sum;
  }
  return Vret;
}


template<typename T>
MyMatrix<T> RankOneMatrix(MyVector<T> const& V)
{
  int n=V.size();
  MyMatrix<T> retMat(n,n);
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      retMat(i,j)=V(i) * V(j);
  return retMat;
}


template<typename T1, typename T2, typename T3>
MyMatrix<T3> MatrixProduct_T1_T2_T3(MyMatrix<T1> const&M1, MyMatrix<T2> const&M2)
{
  int nbCol1=M1.cols();
  int nbRow1=M1.rows();
  int nbCol2=M2.cols();
  int nbRow2=M2.rows();
  if (nbCol1 != nbRow2) {
    std::cerr << "Error in matrix sizes\n";
    throw TerminalException{1};
  }
  MyMatrix<T3> TheProd(nbRow1, nbCol2);
  for (int iCol=0; iCol<nbCol2; iCol++)
    for (int iRow=0; iRow<nbRow1; iRow++) {
      T3 eSum=0;
      for (int i=0; i<nbCol1; i++) {
	T1 prov1=M1(iRow, i);
	T2 prov2=M2(i, iCol);
	T3 prov12=prov1*prov2;
	eSum += prov12;
      }
      TheProd(iRow, iCol)=eSum;
    }
  return TheProd;
}




template<typename T1, typename T2>
MyMatrix<T1> MatrixProduct_T_Oth(MyMatrix<T1> const&M1, MyMatrix<T2> const&M2)
{
  int nbCol1=M1.cols();
  int nbRow1=M1.rows();
  int nbCol2=M2.cols();
  int nbRow2=M2.rows();
  if (nbCol1 != nbRow2) {
    std::cerr << "Error in matrix sizes\n";
    throw TerminalException{1};
  }
  MyMatrix<T1> TheProd(nbRow1, nbCol2);
  for (int iCol=0; iCol<nbCol2; iCol++)
    for (int iRow=0; iRow<nbRow1; iRow++) {
      T1 eSum=0;
      for (int i=0; i<nbCol1; i++) {
	T1 prov1=M1(iRow, i);
	T2 prov2=M2(i, iCol);
	T1 prov12=prov1*prov2;
	eSum += prov12;
      }
      TheProd(iRow, iCol)=eSum;
    }
  return TheProd;
}







template<typename T1, typename T2>
MyMatrix<T2> MatrixProduct_Oth_T(MyMatrix<T1> const&M1, MyMatrix<T2> const&M2)
{
  int nbCol1=M1.cols();
  int nbRow1=M1.rows();
  int nbCol2=M2.cols();
  int nbRow2=M2.rows();
  if (nbCol1 != nbRow2) {
    std::cerr << "Error in matrix sizes\n";
    throw TerminalException{1};
  }
  MyMatrix<T2> TheProd(nbRow1, nbCol2);
  for (int iCol=0; iCol<nbCol2; iCol++)
    for (int iRow=0; iRow<nbRow1; iRow++) {
      T2 eSum=0;
      for (int i=0; i<nbCol1; i++)
	eSum += M1(iRow, i) * M2(i, iCol);
      TheProd(iRow, iCol)=eSum;
    }
  return TheProd;
}


template<typename T>
std::vector<T> RootPolynomial(std::vector<T> const& eVect)
{
  int nbEnt=eVect.size();
  MyMatrix<T> CompMat=ZeroMatrix<T>(nbEnt,nbEnt);
  for (int i=0; i<nbEnt-1; i++)
    CompMat(i+1,i)=1;
  for (int i=0; i<nbEnt; i++)
    CompMat(i,nbEnt-1)=-eVect[i];
  Eigen::EigenSolver<MyMatrix<T>> eig(CompMat);
  MyVector<T> ListEig=eig.eigenvalues();
  std::vector<T> ListRoot(nbEnt);
  for (int i=0; i<nbEnt; i++)
    ListRoot[i]=ListEig(i);
  return ListRoot;
}




// form the product x*A with x a line vector and A a matrix
template<typename T>
MyVector<T> VectorMatrix(MyVector<T> const& eVect, MyMatrix<T> const& eMat)
{
  int nbCol=eMat.cols();
  int nbRow=eMat.rows();
  int n=eVect.size();
  if (n != nbRow) {
    std::cerr << "n should be equal to nbRow\n";
    throw TerminalException{1};
  }
  MyVector<T> rVect(nbCol);
  for (int iCol=0; iCol<nbCol; iCol++) {
    T eSum=0;
    for (int iRow=0; iRow<nbRow; iRow++)
      eSum += eMat(iRow, iCol) * eVect(iRow);
    rVect(iCol)=eSum;
  }
  return rVect;
}

template<typename T>
void SwapValues(T& val1, T& val2)
{
  T prov;
  prov=val1;
  val1=val2;
  val2=prov;
}

template<typename T>
void AssignMatrixRow(MyMatrix<T> &eMat, int const& iRow, MyVector<T> const& eVect)
{
  int nbCol=eMat.cols();
  for (int iCol=0; iCol<nbCol; iCol++)
    eMat(iRow, iCol)=eVect(iCol);
}



template<typename T>
bool IsVectorPositiveMultiple(MyVector<T> const&eVec1, MyVector<T> const&eVec2)
{
  int n=eVec1.size();
  bool IsAssign=false;
  T val2save, val1save;
  for (int i=0; i<n; i++) {
    T eVal1=eVec1(i);
    T eVal2=eVec2(i);
    if (!IsAssign) {
      if (eVal1 != 0) {
	if (eVal2 == 0)
	  return false;
	val2save=eVal2;
	val1save=eVal1;
	IsAssign=true;
      }
      else {
	if (eVal2 != 0)
	  return false;
      }
    }
    else {
      T eDiff=val2save * eVal1 - val1save * eVal2;
      if (eDiff != 0)
	return false;
    }
  }
  if (!IsAssign) {
    std::cerr << "Vectors eVec1 is 0!\n";
    throw TerminalException{1};
  }
  if (val2save > 0 && val1save > 0)
    return true;
  if (val2save < 0 && val1save < 0)
    return true;
  return false;
}


template<typename T>
bool IsVectorMultiple(MyVector<T> const&eVec1, MyVector<T> const&eVec2)
{
  int n=eVec1.size();
  bool IsAssign=false;
  T val2save, val1save;
  for (int i=0; i<n; i++) {
    T eVal1=eVec1(i);
    T eVal2=eVec2(i);
    if (!IsAssign) {
      if (eVal1 != 0) {
	if (eVal2 == 0)
	  return false;
	val2save=eVal2;
	val1save=eVal1;
	IsAssign=true;
      }
      else {
	if (eVal2 != 0)
	  return false;
      }
    }
    else {
      T eDiff=val2save*eVal1 - val1save*eVal2;
      if (eDiff != 0)
	return false;
    }
  }
  if (!IsAssign) {
    std::cerr << "Vectors eVec1 is 0!\n";
    throw TerminalException{1};
  }
  return true;
}


template<typename T>
T DivideVector(MyVector<T> const& V1, MyVector<T> const& V2)
{
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int n=V1.size();
  for (int i=0; i<n; i++)
    if (V1(i) != 0)
      return V1(i)/V2(i);
  std::cerr << "Error in DivideVector\n";
  throw TerminalException{1};
}



template<typename T>
struct Inverse_exception {
  std::string errmsg;
  T pivot;
};

template<typename T>
void TMat_Inverse_destroy(MyMatrix<T> &Input, MyMatrix<T> &Output)
{
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int iCol, iRow;
  int iColB;
  int nbRow=Input.rows();
  int nbCol=Input.cols();
  T prov1, prov2, eVal;
  if (nbRow != nbCol) {
    std::cerr << "Error on nbRow, nbCol in TMat_Inverse_destroy";
    throw TerminalException{1};
  }
  for (iRow=0; iRow<nbRow; iRow++)
    for (iCol=0; iCol<nbRow; iCol++) {
      if (iRow == iCol)
	prov1=1;
      else
	prov1=0;
      Output(iRow,iCol)=prov1;
    }
  int iRowFound=-400;
  for (iCol=0; iCol<nbCol; iCol++) {
    prov1=0;
    for (iRow=iCol; iRow<nbRow; iRow++)
      if (prov1 == 0) {
	eVal=Input(iRow,iCol);
	if (eVal != 0) {
	  iRowFound=iRow;
	  prov1=1/eVal;
	}
      }
    if (prov1 == 0) {
      Inverse_exception<T> eExcept;
      eExcept.errmsg="Error in matrix inversion";
      eExcept.pivot=0;
      throw eExcept;
    }
    for (iColB=0; iColB<nbCol; iColB++) {
      Input(iRowFound,iColB) *= prov1;
      Output(iRowFound,iColB) *= prov1;
    }
    for (iRow=0; iRow<nbRow; iRow++)
      if (iRow != iRowFound) {
	prov2=Input(iRow, iCol);
	if (prov2 != 0) {
	  for (iColB=0; iColB<nbCol; iColB++) {
	    Input(iRow, iColB) -= prov2*Input(iRowFound,iColB);
	    Output(iRow,iColB) -= prov2*Output(iRowFound,iColB);
	  }
	}
      }
    if (iRowFound != iCol)
      for (iColB=0; iColB<nbCol; iColB++) {
	SwapValues(Input(iRowFound, iColB), Input(iCol, iColB));
	SwapValues(Output(iRowFound, iColB), Output(iCol, iColB));
      }
  }
}

template<typename T>
MyMatrix<T> InverseKernel(MyMatrix<T> const&Input)
{
  int nbRow=Input.rows();
  MyMatrix<T> provMat=Input;
  MyMatrix<T> Output(nbRow, nbRow);
  TMat_Inverse_destroy(provMat, Output);
  return Output;
}


template<typename T>
inline typename std::enable_if<is_ring_field<T>::value,MyMatrix<T>>::type Inverse(MyMatrix<T> const& Input)
{
  return InverseKernel(Input);
}

template<typename T>
inline typename std::enable_if<(not is_ring_field<T>::value),MyMatrix<T>>::type Inverse(MyMatrix<T> const& Input)
{
  using Tfield=typename overlying_field<T>::field_type;
  MyMatrix<Tfield> InputF = ConvertMatrixUniversal<Tfield,T>(Input);
  MyMatrix<Tfield> OutputF=InverseKernel(InputF);
  return ConvertMatrixUniversal<T,Tfield>(OutputF);
}









/* This function is for rank calculation.
   Of course, it can be used for many other purpose:
   1> Selecting specific sets of rows and columns for reduction 
   *   2> Computing a set of generators of the zero set.
   Initial matrix is of the form
   Input (nbRow, nbCol)
   The zero matrix is of the form
   NSP (dimKer, nbCol)
*/
template<typename T>
struct SelectionRowCol {
  int TheRank;
  MyMatrix<T> NSP;
  std::vector<int> ListColSelect;
  std::vector<int> ListRowSelect;
};




// The NSP array is assigned to NullspaceMat(TransposedMat(Input)) 
template<typename T>
SelectionRowCol<T> TMat_SelectRowCol(MyMatrix<T> const&Input)
{
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int nbRow=Input.rows();
  int nbCol=Input.cols();
  int maxRank=nbRow;
  if (nbCol < maxRank)
    maxRank=nbCol;
  int sizMat=maxRank+1;
  MyMatrix<T> provMat(sizMat, nbCol);
  std::vector<int> ListColSelect;
  std::vector<int> ListRowSelect;
  std::vector<int> ListColSelect01(nbCol,0);
  int eRank=0;
  for (int iRow=0; iRow<nbRow; iRow++) {
    for (int iCol=0; iCol<nbCol; iCol++)
      provMat(eRank, iCol)=Input(iRow, iCol);
    for (int iRank=0; iRank<eRank; iRank++) {
      int eCol=ListColSelect[iRank];
      T eVal1=provMat(eRank, eCol);
      if (eVal1 != 0) {
	for (int iCol=0; iCol<nbCol; iCol++)
	  provMat(eRank, iCol) -= eVal1*provMat(iRank,iCol);
      }
    }
    int FirstNonZeroCol=-1;
    for (int iCol=0; iCol<nbCol; iCol++)
      if (FirstNonZeroCol == -1) {
	T eVal=provMat(eRank, iCol);
	if (eVal != 0)
	  FirstNonZeroCol=iCol;
      }
    if (FirstNonZeroCol != -1) {
      ListColSelect.push_back(FirstNonZeroCol);
      ListRowSelect.push_back(iRow);
      ListColSelect01[FirstNonZeroCol]=1;
      T eVal=provMat(eRank, FirstNonZeroCol);
      T eVal2=1/eVal;
      for (int iCol=0; iCol<nbCol; iCol++)
	provMat(eRank, iCol) *= eVal2;
      for (int iRank=0; iRank<eRank; iRank++) {
	T eVal1=provMat(iRank, FirstNonZeroCol);
	if (eVal1 != 0) {
	  for (int iCol=0; iCol<nbCol; iCol++)
	    provMat(iRank, iCol) -= eVal1*provMat(eRank, iCol);
	}
      }
      eRank++;
    }
  }
  int nbVectZero=nbCol-eRank;
  MyMatrix<T> NSP=ZeroMatrix<T>(nbVectZero, nbCol);
  int nbVect=0;
  for (int iCol=0; iCol<nbCol; iCol++)
    if (ListColSelect01[iCol] == 0) {
      NSP(nbVect, iCol)=1;
      for (int iRank=0; iRank<eRank; iRank++) {
	int eCol=ListColSelect[iRank];
	NSP(nbVect, eCol)=-provMat(iRank, iCol);
      }
      nbVect++;
    }
  SelectionRowCol<T> retSelect{eRank, NSP, ListColSelect, ListRowSelect};
  return retSelect;
}

template<typename T>
MyMatrix<T> NullspaceTrMat(MyMatrix<T> const& M)
{
  return TMat_SelectRowCol(M).NSP;
}

template<typename T>
MyMatrix<T> NullspaceMat(MyMatrix<T> const& M)
{
  return TMat_SelectRowCol(TransposedMat(M)).NSP;
}





template<typename T>
int RankMatKernel(MyMatrix<T> const& Input)
{
  SelectionRowCol<T> eSelect=TMat_SelectRowCol(Input);
  return eSelect.TheRank;
}


template<typename T>
inline typename std::enable_if<is_ring_field<T>::value,int>::type RankMat(MyMatrix<T> const& Input)
{
  return RankMatKernel(Input);
}

template<typename T>
inline typename std::enable_if<(not is_ring_field<T>::value),int>::type RankMat(MyMatrix<T> const& Input)
{
  using Tfield=typename overlying_field<T>::field_type;
  MyMatrix<Tfield> InputF = ConvertMatrixUniversal<Tfield,T>(Input);
  return RankMatKernel(InputF);
}






template<typename T>
MyMatrix<T> IdentityMat(int n)
{
  T t;
  MyMatrix<T> TheMat(n, n);
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) {
      if (i == j)
	t=1;
      else
	t=0;
      TheMat(i, j)=t;
    }
  return TheMat;
}



template<typename T>
void TMat_ImageIntVector(MyVector<T> &eVect, MyMatrix<T> &TheMat, MyVector<T> &eVectImg)
{
  T t, prov1, prov2;
  int iCol, nbCol, iRow, nbRow, n;
  n=eVect->n;
  nbRow=TheMat.rows();
  nbCol=TheMat.cols();
  if (n != nbRow) {
    std::cerr << "Error in ImageIntVector\n";
    std::cerr << "n=" << n << " nbRow=" << nbRow << "\n";
    throw TerminalException{1};
  }
  for (iCol=0; iCol<nbCol; iCol++) {
    t=0;
    for (iRow=0; iRow<n; iRow++) {
      prov1=TheMat(iRow, iCol);
      prov2=prov1*eVect(iRow);
      t=t + prov2;
    }
    eVectImg(iCol)=t;
  }
}

template<typename T>
MyMatrix<T> TMat_CongrMap(MyMatrix<T> const&eMat)
{
  MyMatrix<T> TheInv=Inverse(eMat);
  return TransposedMat(TheInv);
}

template<typename T>
T DeterminantMatKernel(MyMatrix<T> const&TheMat)
{
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int n, i, j, k, jSel;
  T TheDet, hVal, alpha;
  T eVal1, eVal2, nVal;
  n=TheMat.rows();
  MyMatrix<T> WorkMat = TheMat;
  std::vector<int> eVectPos(n,-1);
  TheDet=1;
  for (i=0; i<n; i++) {
    jSel=-1;
    for (j=0; j<n; j++)
      if (eVectPos[j] == -1) {
	hVal=WorkMat(i, j);
	if (hVal != 0)
	  jSel=j;
      }
    if (jSel == -1) {
      TheDet=0;
      return TheDet;
    }
    eVectPos[jSel]=i;
    for (j=0; j<n; j++)
      if (j != jSel) {
	eVal1=WorkMat(i, j);
	eVal2=WorkMat(i, jSel);
	alpha=eVal1/eVal2;
	for (k=0; k<n; k++) {
	  nVal=WorkMat(k, j) - alpha*WorkMat(k, jSel);
	  WorkMat(k, j)=nVal;
	}
      }
    TheDet=TheDet*WorkMat(i, jSel);
  }
  for (i=0; i<n-1; i++)
    for (j=i+1; j<n; j++)
      if (eVectPos[i] > eVectPos[j])
	TheDet=-TheDet;
  return TheDet;
}

template<typename T>
inline typename std::enable_if<is_ring_field<T>::value,T>::type DeterminantMat(MyMatrix<T> const& Input)
{
  return DeterminantMatKernel(Input);
}

template<typename T>
inline typename std::enable_if<(not is_ring_field<T>::value),T>::type DeterminantMat(MyMatrix<T> const& Input)
{
  using Tfield=typename overlying_field<T>::field_type;
  MyMatrix<Tfield> InputF=ConvertMatrixUniversal<Tfield,T>(Input);
  Tfield eDet_field=DeterminantMatKernel(InputF);
  T eDet=UniversalTypeConversion<T,Tfield>(eDet_field);
  return eDet;
}



template<typename Tint, typename Tfloat>
std::pair<Tfloat,MyVector<Tint>> FindBestIntegerApproximation(MyVector<Tfloat> const& V, int const& N)
{
  int n=V.size();
  Tfloat SumAbsVal=0;
  for (int i=0; i<n; i++)
    SumAbsVal += T_abs(V(i));
  Tfloat MinErr=SumAbsVal;
  MyVector<Tint> MinimalApprox(n);
  //
  for (int iter=1; iter<N; iter++) {
    MyVector<Tint> CandApprox(n);
    Tfloat SumErr=0;
    Tfloat TheMult = iter;
    for (int i=0; i<n; i++) {
      Tfloat val= TheMult * V(i);
      Tint val_i = UniversalNearestInteger<Tint,Tfloat>(val);
      CandApprox(i) = val_i;
      Tfloat valApprox_d = UniversalTypeConversion<Tfloat,Tint>(val_i) / TheMult;
      SumErr += T_abs(valApprox_d - V(i));
    }
    //
    if (SumErr < MinErr) {
      MinErr=SumErr;
      MinimalApprox = CandApprox;
    }
  }
  return {MinErr,MinimalApprox};
}






template<typename T>
MyVector<T> OrthogonalHyperplane(MyMatrix<T> const& M)
{
  int nbVert=M.rows();
  MyVector<T> eVect(nbVert);
  MyMatrix<T> Mred(nbVert-1,nbVert-1);
  int eCoeff=1;
  for (int iVert=0; iVert<nbVert; iVert++) {
    int iRow=0; 
    for (int iLine=0; iLine<nbVert; iLine++)
      if (iLine != iVert) {
	for (int iCol=0; iCol<nbVert-1; iCol++)
	  Mred(iRow,iCol) = M(iLine,iCol);
	iRow++;
      }
    eVect(iVert) = eCoeff * DeterminantMat(Mred);
    eCoeff *= -1;
  }
  return eVect;
}



template<typename T>
MyMatrix<T> SelectRow(MyMatrix<T> const&TheMat, std::vector<int> const& eList)
{
  int nbRowRed=eList.size();
  int nbCol=TheMat.cols();
  MyMatrix<T> TheProv(nbRowRed, nbCol);
  for (int iRow=0; iRow<nbRowRed; iRow++) {
    int jRow=eList[iRow];
    for (int iCol=0; iCol<nbCol; iCol++)
      TheProv(iRow, iCol)=TheMat(jRow, iCol);
  }
  return TheProv;
}



template<typename T>
MyMatrix<T> SelectColumn(MyMatrix<T> const& TheMat, std::vector<int> const& eList)
{
  int nbRow=TheMat.rows();
  int nbColRed=eList.size();
  MyMatrix<T> TheProv(nbRow, nbColRed);
  for (int iCol=0; iCol<nbColRed; iCol++) {
    int jCol=eList[iCol];
    for (int iRow=0; iRow<nbRow; iRow++)
      TheProv(iRow, iCol)=TheMat(iRow, jCol);
  }
  return TheProv;
}


template<typename T>
MyVector<T> SelectColumnVector(MyVector<T> const& TheV, std::vector<int> const& eList)
{
  int nbColRed=eList.size();
  MyVector<T> TheProv(nbColRed);
  for (int iCol=0; iCol<nbColRed; iCol++) {
    int jCol=eList[iCol];
    TheProv(iCol)=TheV(jCol);
  }
  return TheProv;
}




template<typename T>
bool TestEquality(MyVector<T> const& V1, MyVector<T> const& V2)
{
  int n1=V1.size();
  int n2=V2.size();
  if (n1 != n2)
    return false;
  for (int i=0; i<n1; i++) {
    T eVal1=V1(i);
    T eVal2=V2(i);
    if (eVal1 != eVal2)
      return false;
  }
  return true;
}


template<typename T>
bool TestEqualityMatrix(MyMatrix<T> const& M1, MyMatrix<T> const& M2)
{
  int n1=M1.size();
  int n2=M2.size();
  if (n1 != n2)
    return false;
  for (int i=0; i<n1; i++)
    if (M1(i) != M2(i))
      return false;
  return true;
}







// In input a matrix
// in output the lines that spann its rank in sequential order
template<typename T>
MyMatrix<T> RowReduction(MyMatrix<T> const&eMatIn)
{
  SelectionRowCol<T> eSelect=TMat_SelectRowCol(eMatIn);
  std::vector<int> ListRowSelect=eSelect.ListRowSelect;
  return SelectRow(eMatIn, ListRowSelect);
}


template<typename T>
struct SolMatResult {
  bool result;
  MyVector<T> eSol;
};


// Given the equation Y = XA, we find one solution X if it exists.
// 
template<typename T>
SolMatResult<T> SolutionMat(MyMatrix<T> const& eMat, MyVector<T> const& eVect)
{
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  if (eMat.rows() == 0) {
    if (!IsZeroVector(eVect))
      return {false,{}};
    MyVector<T> eSol(0);
    return {true, eSol};
  }
  int nbRow=eMat.rows();
  int nbCol=eMat.cols();
  SelectionRowCol<T> eSelect=TMat_SelectRowCol(eMat);
  int eRank=eSelect.TheRank;
  std::vector<int> ListRowSelect=eSelect.ListRowSelect;
  std::vector<int> ListColSelect=eSelect.ListColSelect;
  MyMatrix<T> eMat2=SelectRow(eMat, ListRowSelect);
  MyMatrix<T> eMat3=SelectColumn(eMat2, ListColSelect);
  MyVector<T> eVectB=SelectColumnVector(eVect, ListColSelect);
  MyMatrix<T> eMatInv=Inverse(eMat3);
  MyVector<T> eSol=ProductVectorMatrix(eVectB, eMatInv);
  MyVector<T> eProd=ProductVectorMatrix(eSol, eMat2);
  for (int iCol=0; iCol<nbCol; iCol++)
    if (eProd(iCol) != eVect(iCol))
      return {false, {}};
  MyVector<T> eRetSol=ZeroVector<T>(nbRow);
  for (int iRank=0; iRank<eRank; iRank++) {
    int iRow=ListRowSelect[iRank];
    eRetSol(iRow)=eSol(iRank);
  }
  return {true, eRetSol};
}


template<typename T>
MyMatrix<T> SelectNonZeroRows(MyMatrix<T> const& EXT)
{
  int nbRow=EXT.rows();
  int nbCol=EXT.cols();
  std::vector<int> ListIdx;
  for (int iRow=0; iRow<nbRow; iRow++) {
    bool IsZero=true;
    for (int iCol=0; iCol<nbCol; iCol++)
      if (EXT(iRow,iCol) != 0)
	IsZero=false;
    if (!IsZero)
      ListIdx.push_back(iRow);
  }
  return SelectRow(EXT, ListIdx);
}




template<typename T>
std::vector<int> ColumnReductionSet(MyMatrix<T> const& eMatIn)
{
  SelectionRowCol<T> eSelect=TMat_SelectRowCol(eMatIn);
  return eSelect.ListColSelect;
}

template<typename T>
MyMatrix<T> ColRowSymmetricMatrix(MyMatrix<T> const& M, std::vector<int> const& LSel)
{
  int siz=LSel.size();
  MyMatrix<T> Mred(siz, siz);
  for (int i=0; i<siz; i++)
    for (int j=0; j<siz; j++) {
      int iFull=LSel[i];
      int jFull=LSel[j];
      Mred(i,j)=M(iFull,jFull);
    }
  return Mred;
}





template<typename T>
MyMatrix<T> ColumnReduction(MyMatrix<T> const& eMatIn)
{
  SelectionRowCol<T> eSelect=TMat_SelectRowCol(eMatIn);
  return SelectColumn(eMatIn, eSelect.ListColSelect);
}





template<typename T>
MyMatrix<T> Concatenate(MyMatrix<T> const& eMat1, MyMatrix<T> const& eMat2)
{
  int nbRow1=eMat1.rows();
  int nbRow2=eMat2.rows();
  int nbCol1=eMat1.cols();
  int nbCol2=eMat2.cols();
  if (nbCol1 != nbCol2) {
    std::cerr << "nbCol1=" << nbCol1 << " nbCol2=" << nbCol2 << "\n";
    std::cerr << "Error in the number of columns\n";
    throw TerminalException{1};
  }
  int nbCol=nbCol1;
  MyMatrix<T> eMatRet(nbRow1+nbRow2, nbCol);
  for (int iRow=0; iRow<nbRow1; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      T eVal=eMat1(iRow, iCol);
      eMatRet(iRow, iCol)=eVal;
    }
  for (int iRow=0; iRow<nbRow2; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      T eVal=eMat2(iRow, iCol);
      eMatRet(iRow+nbRow1, iCol)=eVal;
    }
  return eMatRet;
}

template<typename T>
MyMatrix<T> ConcatenateMatVec(MyMatrix<T> const& M, MyVector<T> const& V)
{
  int nbRow=M.rows();
  int nbColM=M.cols();
  int nbCol=V.size();
  if (nbRow != 0) {
    if (nbColM != nbCol) {
      std::cerr << "Error in ConcatenateMatVec\n";
      std::cerr << "We have nbCol=" << nbCol << " nbColM=" << nbColM << "\n";
      throw TerminalException{1};
    }
  }
  //  std::cerr << "M(rows/cols)=" << M.rows() << "/" << M.cols() << "\n";
  //  std::cerr << "V(size)=" << V.size() << "\n";
  MyMatrix<T> Mret(nbRow+1,nbCol);
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++)
      Mret(iRow,iCol)=M(iRow,iCol);
  for (int iCol=0; iCol<nbCol; iCol++)
    Mret(nbRow,iCol)=V(iCol);
  return Mret;
}





template<typename T>
MySparseMatrix<T> ReadSparseMatrix(std::istream &is)
{
  if (!is.good()) {
    std::cerr << "ReadSparseMatrix operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  T eVal;
  int nbRow, nbCol, nnz;
  int iRow, iCol;
  is >> nbRow >> nbCol >> nnz;
  //  std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << " nnz=" << nnz << "\n";
  using T2=Eigen::Triplet<T>;
  std::vector<T2> tripletList(nnz);
  for (int iNNZ=0; iNNZ<nnz; iNNZ++) {
    is >> iRow >> iCol >> eVal;
    tripletList[iNNZ]=T2(iRow,iCol,eVal);
    if (iRow >= nbRow) {
      std::cerr << "iRow is too large\n";
      std::cerr << " iRow=" << iRow << "\n";
      std::cerr << "nbRow=" << nbRow << "\n";
      throw TerminalException{1};
    }
    if (iCol >= nbCol) {
      std::cerr << "iCol is too large\n";
      std::cerr << " iCol=" << iCol << "\n";
      std::cerr << "nbCol=" << nbCol << "\n";
      throw TerminalException{1};
    }
  }
  //  std::cerr << "After the tripletList assignment\n";
  MySparseMatrix<T> SpMat(nbRow, nbCol);
  //  std::cerr << "Creation of the sparse matrix\n";
  SpMat.setFromTriplets(tripletList.begin(), tripletList.end());
  //  std::cerr << "After the setFromTriplets\n";
  //  std::cerr << "(rows/cols)SpMat = " << SpMat.rows() << " / " << SpMat.cols() << "\n";
  return SpMat;
}


template<typename T>
void WriteSparseMatrix(std::ostream &os, MySparseMatrix<T> const& eMat)
{
  int nbRow=eMat.rows();
  int nbCol=eMat.cols();
  int nnz=eMat.nonZeros();
  os << nbRow << " " << nbCol << " " << nnz;
  int nb=0;
  for (int k=0; k<eMat.outerSize(); ++k)
    for (typename MySparseMatrix<T>::InnerIterator it(eMat,k); it; ++it) {
      T eVal=it.value();
      int iRow=it.row();   // row index
      int iCol=it.col();   // col index (here it is equal to k)
      os << iRow << " " << iCol << " " << eVal << "\n";
      nb++;
    }
  if (nb != nnz) {
    std::cerr << "Logical error in our understanding of \n";
    throw TerminalException{1};
  }
}



template<typename T>
MyMatrix<T> MyMatrixFromSparseMatrix(MySparseMatrix<T> const& eMat)
{
  int nbRow=eMat.rows();
  int nbCol=eMat.cols();
  int nnz=eMat.nonZeros();
  MyMatrix<T> M=ZeroMatrix<T>(nbRow,nbCol);
  for (int k=0; k<eMat.outerSize(); ++k)
    for (typename MySparseMatrix<T>::InnerIterator it(eMat,k); it; ++it) {
      T eVal=it.value();
      int iRow=it.row();   // row index
      int iCol=it.col();   // col index (here it is equal to k)
      M(iRow,iCol) = eVal;
    }
  return M;
}











template<typename T>
T MatrixScalarProduct(MyMatrix<T> const& M1, MyMatrix<T> const& M2)
{
  int n1=M1.rows();
  int p1=M1.cols();
  int n2=M2.rows();
  int p2=M2.cols();
  if (n1 != n2 || p1 != p2) {
    std::cerr << "Incoherency\n";
    std::cerr << "n1=" << n1 << " n2=" << n2 << "\n";
    std::cerr << "p1=" << p1 << " p1=" << p2 << "\n";
    throw TerminalException{1};
  }
  T eSum=0;
  for (int i=0; i<n1; i++)
    for (int j=0; j<p1; j++)
      eSum += M1(i,j)*M2(i,j);
  return eSum;
}

template<typename T>
MyMatrix<T> VectorToSymmetricMatrix(MyVector<T> const& V, int const& n)
{
  MyMatrix<T> RetMat=ZeroMatrix<T>(n, n);
  int idx=0;
  for (int i=0; i<n; i++)
    for (int j=0; j<=i; j++) {
      RetMat(i,j)=V(idx);
      RetMat(j,i)=V(idx);
      idx++;
    }
  return RetMat;
}

template<typename T>
MyMatrix<T> VectorToSymmetricMatrixB(MyVector<T> const& V, int const& n)
{
  MyMatrix<T> RetMat(n,n);
  int idx=0;
  for (int i=0; i<n; i++)
    for (int j=0; j<=i; j++) {
      T eVal=V(idx);
      if (i != j)
	eVal=eVal/2;
      RetMat(i,j)=eVal;
      RetMat(j,i)=eVal;
      idx++;
    }
  return RetMat;
}


template<typename T>
MyVector<T> SymmetricMatrixToVector(MyMatrix<T> const& M)
{
  int n=M.rows();
  int dim=(n*(n+1))/2;
  MyVector<T> eVect(dim);
  int idx=0;
  for (int i=0; i<n; i++)
    for (int j=0; j<=i; j++) {
      eVect(idx)=M(i,j);
      idx++;
    }
  return eVect;
}

template<typename T>
MyVector<T> GetSymmetricMatrixWeightVector(int const& n)
{
  int dim=(n*(n+1))/2;
  MyVector<T> eVect(dim);
  int idx=0;
  for (int i=0; i<n; i++)
    for (int j=0; j<=i; j++) {
      T eVal;
      if (i == j)
	eVal = 1;
      else
	eVal = 2;
      eVect(idx)=eVal;
      idx++;
    }
  return eVect;
}





template<typename T>
MyVector<T> SymmetricMatrixToVectorB(MyMatrix<T> const& M)
{
  int n=M.rows();
  int dim=(n*(n+1))/2;
  MyVector<T> eVect(dim);
  int idx=0;
  T eVal;
  for (int i=0; i<n; i++)
    for (int j=0; j<=i; j++) {
      eVal=M(i,j);
      if (i != j)
	eVal=2*eVal;
      eVect(idx)=eVal;
    }
  return eVect;
}



template<typename T>
void PrintEigenvalueDefect(MyMatrix<T> const& Sinp, MyVector<T> const& ListEigVal, MyMatrix<T> const& ListEigVect, std::ostream & os)
{
  int n=Sinp.rows();
  for (int i=0; i<n; i++) {
    T eDelta=0;
    for (int j=0; j<n; j++) {
      T eSum=ListEigVal(i) * ListEigVect(i,j);
      for (int k=0; k<n; k++)
	eSum -= ListEigVect(i,k) * Sinp(j,k);
      eDelta += T_abs(eSum);
    }
    os << "i=" << i << " err=" << eDelta << "\n";
  }
}

template<typename T>
MyVector<T> SolveConjGrad(MyMatrix<T> const& A, MyVector<T> const& b)
{
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int n=b.size();
  MyVector<T> r(n);
  MyVector<T> p(n);
  MyVector<T> x(n);
  MyVector<T> Ap(n);
  double rsold, rsnew, alpha;
  r=b;
  p=r;
  rsold=ScalarProduct(r, r);
  T eZer=0;
  for (int i=0; i<n; i++)
    x(i)=eZer;
  int nbOper=4*n;
  for (int i=0; i<nbOper; i++) {
    Ap=VectorMatrix(p, A);
    alpha=rsold/ScalarProductQuadForm(A, p, p);
    x += alpha*p;
    r -= alpha*Ap;
    rsnew=ScalarProduct_Doubl(r, r);
    p = r + (rsnew/rsold)*p;
    rsold=rsnew;
  }
  return x;
}

template<typename T>
MyVector<T> Concatenation(MyVector<T> const& V1, MyVector<T> const& V2)
{
  int dim1=V1.size();
  int dim2=V2.size();
  MyVector<T> V(dim1 + dim2);
  for (int i=0; i<dim1; i++)
    V(i)=V1(i);
  for (int i=0; i<dim2; i++)
    V(i+dim1)=V2(i);
  return V;
}

template<typename T>
MyVector<T> Isobarycenter(MyMatrix<T> const& eMat)
{
  int nbRow=eMat.rows();
  int nbCol=eMat.cols();
  MyVector<T> eVect(nbCol);
  T nbRow_T=nbRow;
  for (int iCol=0; iCol<nbCol; iCol++) {
    T eSum=0;
    for (int iRow=0; iRow<nbRow; iRow++)
      eSum += eMat(iRow, iCol);
    T eVal=eSum/nbRow_T;
    eVect(iCol)=eVal;
  }
  return eVect;
}

template<typename T>
bool IsZeroMatrix(MyMatrix<T> const& M)
{
  int nbRow=M.rows();
  int nbCol=M.cols();
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++)
      if (M(iRow, iCol) != 0)
	return false;
  return true;
}

template<typename T>
bool IsZeroVector(MyVector<T> const& V)
{
  int n=V.size();
  for (int i=0; i<n; i++)
    if (V(i) != 0)
      return false;
  return true;
}


template<typename T>
inline typename std::enable_if<is_ring_field<T>::value,MyVector<T>>::type CanonicalizeVector(MyVector<T> const& V)
{
  int n=V.size();
  T TheMin=0;
  int iSelect=-1;
  for (int i=0; i<n; i++) {
    T eVal=V(i);
    if (eVal != 0) {
      T eAbs=T_abs(eVal);
      if (iSelect == -1) {
	TheMin=eAbs;
	iSelect=i;
      }
      else {
	if (eAbs < TheMin) {
	  TheMin=eAbs;
	  iSelect=i;
	}
      }
    }
  }
  if (iSelect == -1)
    return V;
  MyVector<T> Vret(n);
  T eQuot=1/TheMin;
  for (int i=0; i<n; i++)
    Vret(i)=V(i)*eQuot;
  return Vret;
}






/* return true if V1 < V2 according to lexicographic order */
template<typename T>
bool IsLower(MyVector<T> const& V1, MyVector<T> const& V2)
{
  int n=V1.size();
  for (int i=0; i<n; i++) {
    if (V1(i) != V2(i)) {
      if (V1(i) < V2(i)) {
	return true;
      }
      else {
	return false;
      }
    }
  }
  return false;
}

template<typename T>
bool operator==(MyVector<T> const& V1, MyVector<T> const& V2)
{
  int n=V1.size();
  if (V2.size() != n) {
    std::cerr << "We should not have different sizes\n";
    throw TerminalException{1};
  }
  for (int i=0; i<n; i++) {
    if (V1(i) != V2(i))
      return false;
  }
  return true;
}

template<typename T>
bool operator<(MyMatrix<T> const& M1, MyMatrix<T> const& M2)
{
  int nbRow=M1.rows();
  int nbCol=M1.cols();
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      if (M1(iRow, iCol) < M2(iRow,iCol))
	return true;
      if (M1(iRow, iCol) > M2(iRow,iCol))
	return false;
    }
  return false;
}

template<typename T>
bool operator<(MyVector<T> const& V1, MyVector<T> const& V2)
{
  int siz=V1.size();
  for (int i=0; i<siz; i++) {
    if (V1(i) < V2(i))
      return true;
    if (V1(i) > V2(i))
      return false;
  }
  return false;
}



namespace std {
  template<typename T>
  struct less<MyVector<T>> {
    bool operator()(MyVector<T> const& V1, MyVector<T> const& V2) const
    {
      int siz=V1.size();
      for (int i=0; i<siz; i++) {
	if (V1(i) < V2(i))
	  return true;
	if (V2(i) < V1(i))
	  return false;
      }
      return false;
    }
  };
}




template<typename T>
T L1_norm(MyVector<T> const& V)
{
  int siz=V.size();
  T norm=0;
  for (int i=0; i<siz; i++)
    norm += T_abs(V(i));
  return norm;
}




// Just sorting (no unicization)
//
template<typename T>
MyMatrix<T> SortMatrix(MyMatrix<T> const& M)
{
  int nbRow=M.rows();
  int nbCol=M.cols();
  auto comp=[](MyVector<T> const& V1, MyVector<T> const& V2) -> bool {return IsLower(V1, V2);};
  std::set<MyVector<T>, std::function<bool(MyVector<T>,MyVector<T>)> > eListVect(comp);
  //  eListVect=std::set<MyVector<T>, decltype(comp)> (comp);
  for (int iRow=0; iRow<nbRow; iRow++) {
    MyVector<T> V(nbCol);
    for (int iCol=0; iCol<nbCol; iCol++)
      V(iCol)=M(iRow,iCol);
    eListVect.insert(V);
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





template<typename T>
MyVector<T> GetMatrixRow(MyMatrix<T> const& M, int const& iRow)
{
  int nbCol=M.cols();
  MyVector<T> V(nbCol);
  for (int iCol=0; iCol<nbCol; iCol++)
    V(iCol)=M(iRow,iCol);
  return V;
}

template<typename T>
MyVector<T> GetMatrixColumn(MyMatrix<T> const& M, int const& iCol)
{
  int nbRow=M.cols();
  MyVector<T> V(nbRow);
  for (int iRow=0; iRow<nbRow; iRow++)
    V(iRow)=M(iRow,iCol);
  return V;
}



template<typename T>
std::vector<T> ConcatenateVect(std::vector<T> const& ListV1, std::vector<T> const& ListV2)
{
  std::vector<T> ListV=ListV1;
  for (auto & eV : ListV2)
    ListV.push_back(eV);
  return ListV;
}




template<typename T>
MyMatrix<T> MatrixFromVectorFamily(std::vector<MyVector<T>> const& ListVect)
{
  int nbVect=ListVect.size();
  if (nbVect == 0) {
    std::cerr << "Error in MatrixFromVectorFamily. ListVect is empty\n";
    std::cerr << "We cannot create the matrix since we cannot know the dimension\n";
    throw TerminalException{1};
  }
  int dim=ListVect[0].size();
  MyMatrix<T> M(nbVect,dim);
  for (int iVect=0; iVect<nbVect; iVect++)
    for (int i=0; i<dim; i++)
      M(iVect,i)=ListVect[iVect](i);
  return M;
}

template<typename T>
MyVector<T> SumMatrix(MyMatrix<T> const& M)
{
  int nbRow=M.rows();
  int nbCol=M.cols();
  MyVector<T> V=ZeroVector<T>(nbCol);
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++)
      V(iCol) += M(iRow,iCol);
  return V;
}



template<typename T>
std::vector<T> StdVectorFromVector(MyVector<T> const& eV)
{
  int siz=eV.size();
  std::vector<T> eVect(siz);
  for (int i=0; i<siz; i++)
    eVect[i]=eV(i);
  return eVect;
}

template<typename T>
MyVector<T> VectorFromStdVector(std::vector<T> const& eList)
{
  int siz=eList.size();
  MyVector<T> eVect(siz);
  for (int i=0; i<siz; i++)
    eVect(i)=eList[i];
  return eVect;
}

template<typename T>
inline typename std::enable_if<is_float_arithmetic<T>::value,MyVector<T>>::type RemoveFractionVector(MyVector<T> const& V)
{
  return V;
}


template<typename T, typename Tint>
T EvaluationQuadForm(MyMatrix<T> const& eMat, MyVector<Tint> const& eVect)
{
  int n=eVect.size();
  T eSum=0;
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) {
      Tint eVal12=eVect(i)*eVect(j);
      T eVal=eMat(i,j);
      eSum += eVal12*eVal;
    }
  return eSum;
}

//
// The matrix is supposed to describe a vector space. We want
// to canonicalize it in order to get a better representation
// of the space.
// 
template<typename T>
MyMatrix<T> CanonicalizeBasisVectorSpace(MyMatrix<T> const& inputMat)
{
  MyMatrix<T> WorkMat = inputMat;
  int nbRow=WorkMat.rows();
  int nbCol=WorkMat.cols();
  std::vector<int> ColStatus(nbCol,1);
  std::vector<int> RowStatus(nbRow,1);
  for (int iRow=0; iRow<nbRow; iRow++) {
    int FoundRow=-1;
    int FoundCol=-1;
    T MaxValue=0;
    for (int eRow=0; eRow<nbRow; eRow++) {
      for (int eCol=0; eCol<nbCol; eCol++) {
	if (ColStatus[eCol] == 1 && RowStatus[eRow] == 1) {
	  T aVal = T_abs(WorkMat(eRow, eCol));
	  if (FoundRow == -1) {
	    MaxValue=aVal;
	    FoundRow=eRow;
	    FoundCol=eCol;
	  }
	  else {
	    if (aVal > MaxValue) {
	      MaxValue=aVal;
	      FoundRow=eRow;
	      FoundCol=eCol;
	    }
	  }
	}
      }
    }
    //
    ColStatus[FoundCol]=0;
    RowStatus[FoundRow]=0;
    for (int eRow=0; eRow<nbRow; eRow++) {
      if (eRow != FoundRow) {
	T alpha=WorkMat(eRow, FoundCol) / WorkMat(FoundRow, FoundCol);
	for (int iCol=0; iCol<nbCol; iCol++)
	  WorkMat(eRow, iCol) -= alpha * WorkMat(FoundRow, iCol);
      }
    }
  }
  return WorkMat;
}





template<typename T, typename Tint>
T ScalarProductQuadForm(MyMatrix<T> const& eMat, MyVector<Tint> const& V1, MyVector<Tint> const& V2)
{
  int n=V1.size();
  T eSum=0;
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) {
      Tint eVal12=V1(i) * V2(j);
      T eVal=eMat(i,j);
      eSum += eVal12*eVal;
    }
  return eSum;
}

template<typename T>
MyVector<T> SolutionMat_LeastSquare(MyMatrix<T> const& M, MyVector<double> const& V)
{
  int nbRow=M.rows();
  int nbCol=M.cols();
  MyMatrix<double> Msqr=ZeroMatrix<double>(nbCol, nbCol);
  for (int j1=0; j1<nbCol; j1++)
    for (int j2=0; j2<nbCol; j2++)
      for (int i=0; i<nbRow; i++)
	Msqr(j1,j2) += M(i,j1) * M(i,j2);
  MyVector<double> B=ZeroVector<double>(nbCol);
  for (int j=0; j<nbCol; j++)
    for (int i=0; i<nbRow; i++)
      B(i,j) += V(i) * M(i,j);
  //
  // We have now a linear system to solve
  //
  Eigen::FullPivLU<MyMatrix<double>> solver;
  solver.compute(Msqr);
  MyVector<double> eSol=solver.solve(B);
  //
  return eSol;
}

template<typename T>
bool IsSymmetricMatrix(MyMatrix<T> const& M)
{
  if (M.rows() != M.cols())
    return false;
  int n=M.rows();
  for (int i=0; i<n; i++)
    for (int j=i; j<n; j++)
      if (M(i,j) != M(j,i))
	return false;
  return true;
}



#endif

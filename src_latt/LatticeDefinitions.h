#ifndef INCLUDE_LATTICE_DEFINITIONS
#define INCLUDE_LATTICE_DEFINITIONS

#include "MAT_Matrix.h"
#include "Boost_bitset.h"
#include "COMB_Combinatorics.h"

template<typename T, typename Tint>
struct Tshortest {
  T eMin;
  MyMatrix<Tint> SHV;
};

#undef DEBUG

template<typename T, typename Tint>
Tshortest<T, Tint> SelectShortestVector(MyMatrix<T> const& eMat, MyMatrix<Tint> const& SHV)
{
  int n=eMat.rows();
  int nbRow=SHV.rows();
  std::vector<int> ListStatus(nbRow, 0);
  int nbShort=0;
  T MinNorm = -1;
  for (int iRow=0; iRow<nbRow; iRow++) {
    MyVector<Tint> eVect=SHV.row(iRow);
    T eNorm=EvaluationQuadForm(eMat, eVect);
    if (eNorm > 0) {
      if (nbShort == 0) {
	MinNorm=eNorm;
	nbShort=1;
	ListStatus[iRow]=1;
      }
      else {
	if (eNorm < MinNorm) {
	  MinNorm=eNorm;
	  nbShort=1;
	  for (int jRow=0; jRow<iRow; jRow++)
	    ListStatus[jRow]=0;
	  ListStatus[iRow]=1;
	}
	else {
	  if (eNorm == MinNorm) {
	    nbShort++;
	    ListStatus[iRow]=1;
	  }
	}
      }
    }
  }
  MyMatrix<Tint> TheSHV(nbShort, n);
  int iShort=0;
  for (int iRow=0; iRow<nbRow; iRow++)
    if (ListStatus[iRow] == 1) {
      for (int i=0; i<n; i++)
	TheSHV(iShort, i)=SHV(iRow, i);
      iShort++;
    }
  //  std::cerr << "MinNorm=" << MinNorm << "\n";
  return {MinNorm, TheSHV};
}



template<typename T, typename Tint>
struct resultCVP {
  T TheNorm;
  MyMatrix<Tint> ListVect;
};

template<typename T,typename Tint>
bool operator==(resultCVP<T,Tint> const& x, resultCVP<T,Tint> const& y)
{
  if (x.TheNorm != y.TheNorm)
    return false;
  MyMatrix<Tint> ListVectSort1=SortMatrix(x.ListVect);
  MyMatrix<Tint> ListVectSort2=SortMatrix(y.ListVect);
  return ListVectSort1 == ListVectSort2;
}

template<typename T,typename Tint>
bool operator!=(resultCVP<T,Tint> const& x, resultCVP<T,Tint> const& y)
{
  if (x == y)
    return false;
  return true;
}



template<typename Tint>
void EnumerateOrbitPrimitiveVector(std::vector<MyMatrix<int>> const& ListMat, int const& pPrime, int const& dim, std::function<void(std::vector<int> const&, int const&)> const& FCT)
{
  Tint expo=1;
  for (int i=0; i<dim; i++)
    expo *= pPrime;
  Tint TotalNb=(expo - 1) / (pPrime -1);
  Face StatusVect(TotalNb);
  for (Tint i=0; i<TotalNb; i++)
    StatusVect[i]=1;
  //
  // The modulo p Operations
  //
  auto Canonicalization=[&](int const& x) -> int {
    int res=x % pPrime;
    return res;
  };
  auto ProdModP=[&](int const& a, int const& b) -> int {
    return Canonicalization(a*b);
  };
  std::vector<int> ListInverse(pPrime);
  for (int i=0; i<pPrime; i++) {
    int eInv=-400;
    for (int j=0; j<pPrime; j++) {
      int res = ProdModP(i,j);
      if (res == 1)
	eInv=j;
    }
    ListInverse[i]=eInv;
    //    std::cerr << "i=" << i << " ListInverse=" << eInv << "\n";
  }
  auto PrimitivizeVector=[&](std::vector<int> & V) -> void {
    /*    std::cerr << "PrimitiveVector Vi=";
    for (int i=0; i<dim; i++)
      std::cerr << " " << V[i];
      std::cerr << "\n";*/
    int ifound=-1;
    for (int i=0; i<dim; i++)
      if (V[i] != 0)
	ifound=i;
    int eVal=V[ifound];
    int eInv=ListInverse[eVal];
    //    std::cerr << "eVal=" << eVal << " eInv=" << eInv << "\n";
    for (int i=0; i<=ifound; i++)
      V[i]=ProdModP(eInv, V[i]);
    /*    std::cerr << "PrimitiveVector Vo=";
    for (int i=0; i<dim; i++)
      std::cerr << " " << V[i];
      std::cerr << "\n";*/
  };
  auto ImagePrimitiveVector=[&](MyMatrix<int> const& M, std::vector<int> const& V) -> std::vector<int> {
    std::vector<int> Vret(dim);
    for (int i=0; i<dim; i++) {
      int sum=0;
      for (int j=0; j<dim; j++)
	sum += M(j,i) * V[j];
      int res = Canonicalization(sum);
      //      std::cerr << "sum=" << sum << " pPrime=" << pPrime << " res=" << res << "\n";
      Vret[i]=res;
    }
    PrimitivizeVector(Vret);
    return Vret;
  };
  auto ImageNumber=[&](MyMatrix<int> const& M, Tint const& valIn) -> Tint {
    std::vector<int> Vin=ConvertNumberToPrimitiveVector(valIn, pPrime, dim);
    int NewVal=ConvertPrimitiveVectorToNumber<Tint>(Vin, pPrime);
    if (NewVal != valIn) {
      std::cerr << "Inconsistency to be solved\n";
      throw TerminalException{1};
    }
    /*    std::cerr << "Vin=";
    for (int i=0; i<dim; i++)
      std::cerr << " " << Vin[i];
    std::cerr << "\n";
    std::cerr << "M=\n";
    WriteMatrix(std::cerr, M);*/
    std::vector<int> Vout=ImagePrimitiveVector(M, Vin);
    /*    std::cerr << "Vout=";
    for (int i=0; i<dim; i++)
      std::cerr << " " << Vout[i];
      std::cerr << "\n";*/
    return ConvertPrimitiveVectorToNumber<Tint>(Vout, pPrime);
  };
  //
  // Combinatorial enumeration code
  //
  std::vector<std::vector<int>> ListRepresentative;
  Tint cnt = TotalNb;
  StackStorage<Tint> ActivePoints;
  int ThePos=0;
  auto GetPosition=[&]() -> int {
    while(true) {
      if (StatusVect[ThePos] == 1)
	return ThePos;
      ThePos++;
    }
  };
  int nbOrbit=0;
  while(true) {
    if (cnt == 0)
      break;
    int OrbSize=0;
    int ePos=GetPosition();
    std::vector<int> eRepr=ConvertNumberToPrimitiveVector(ePos, pPrime, dim);
    auto FuncInsert=[&](Tint const& val) -> void {
      if (StatusVect[val] == 1) {
	StatusVect[val]=0;
	ActivePoints.push_back(val);
	cnt--;
	OrbSize++;
      }
    };
    FuncInsert(ePos);
    while(true) {
      if (ActivePoints.size() == 0)
	break;
      Tint ePt=ActivePoints.pop();
      for (auto & M : ListMat) {
	Tint fPt=ImageNumber(M, ePt);
	FuncInsert(fPt);
      }
    }
    FCT(eRepr, OrbSize);
    std::cerr << "Find an orbit of size " << OrbSize << "\n";
    nbOrbit++;
  }
  std::cerr << "nbOrbit=" << nbOrbit << "\n";
}


template<typename T, typename Tint>
struct LLLreduction {
  MyMatrix<T> GramMatRed;
  MyMatrix<Tint> Pmat;
};



//
// Adapted from LLLReducedBasis   in zlattice.gi GAP code
//
template<typename Tmat, typename Tint>
LLLreduction<Tmat,Tint> LLLreducedBasis(MyMatrix<Tmat> const & GramMat)
{
  using Tfield=typename overlying_field<Tmat>::field_type;
  MyMatrix<Tmat> gram = GramMat;
  int nbRow=gram.rows();
  int nbCol=gram.cols();
  //  std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
  if (nbRow != nbCol) {
    throw TerminalException{1};
  }
  int n = nbRow;
  int k = 1;
  int kmax = 0;
  MyMatrix<Tfield> mue = ZeroMatrix<Tfield>(n,n);
  int r = -1;
  MyVector<Tfield> ak(n);
  MyMatrix<Tint> H = IdentityMat<Tint>(n);
  auto RED=[&](int const& l) -> void {
    //    std::cerr << "k=" << k << " l=" << l << " mue(k,l)=" << mue(k,l) << "\n";
    if (1 < mue(k,l) * 2 || mue(k,l) * 2 < -1) {
      Tint q=UniversalNearestInteger<Tint,Tfield>(mue(k,l));
      //      std::cerr << "RED, before oper q=" << q << "\n";
      //      WriteMatrix(std::cerr, gram);
      gram(k,k) -= q * gram(k,l);
      //      std::cerr << "RED step 1\n";
      for (int i=r+1; i<=l; i++)
        gram(k,i) -= q * gram(l,i);
      //      std::cerr << "RED step 2\n";
      for (int i=l+1; i<=k; i++)
        gram(k,i) -= q * gram(i,l);
      //      std::cerr << "RED step 3 n=" << n << "\n";
      for (int i=k+1; i<n; i++)
        gram(i,k) -= q * gram(i,l);
      //      std::cerr << "After gram Oper\n";
      //      WriteMatrix(std::cerr, gram);
      //      std::cerr << "RED step 4\n";
      mue(k,l) = mue(k,l) - q;
      //      std::cerr << "RED step 5\n";
      for (int i=r+1; i<=l-1; i++)
        mue(k,i) -= q * mue(l,i);
      //      std::cerr << "RED step 6\n";
      H.row(k) -= q * H.row(l);
      //      std::cerr << "RED step 7\n";
    }
  };
  Tfield y=Tfield(99)/Tfield(100);
  //  std::cerr << " y=" << y << "\n";
  int i=0;
  while(true) {
    Tmat eVal = gram(i,i);
    //    std::cerr << "i=" << i << " eVal=" << eVal << "\n";
    if (eVal > 0)
      break;
    i++;
    if (i == n)
      break;
  }
  //  std::cerr << "After the while loop i=" << i << "\n";
  if (i >= n) {
    r = n-1;
    k = n;
  }
  else {
    if (i > 0) {
      for (int j=i+1; j<n; j++) {
        gram(j,0) = gram(j,i);
        gram(j,i) = 0;
      }
      gram(0,0) = gram(i,i);
      gram(i,i)=0;
      H.row(i).swap(H.row(i));
    }
  }
  //  std::cerr << "After the if test r=" << r << " k=" << k << " kmax=" << kmax << "\n";
  MyVector<Tfield> B(n);
  B(0) = UniversalTypeConversion<Tfield,Tmat>(gram(0,0));
  //  std::cerr << "Before while loop\n";
  while (k < n) {
    //    std::cerr << "While loop, step 1 k=" << k << " kmax=" << kmax << "\n";
    if (k > kmax) {
      kmax = k;
      B(k) = UniversalTypeConversion<Tfield,Tmat>(gram(k,k));
      for (int u=0; u<n; u++)
	mue(k,u) = 0;
      for (int j=r+1; j<=k-1; j++) {
        ak(j) = UniversalTypeConversion<Tfield,Tmat>(gram(k,j));
        for (int i=r+1; i<=j-1; i++)
          ak(j) -= mue(j,i) * ak(i);
        mue(k,j) = ak(j) / B(j);
	//	std::cerr << "j=" << j << " ak(j)=" << ak(j) << " B(j) = " << B(j) << " mue(k,l)=" << mue(k,j) << "\n";
        B(k) -= mue(k,j) * ak(j);
      }
    }
    //    std::cerr << "While loop, step 2 k=" << k << "\n";
    RED(k-1);
    //    bool test=B(k) < ( y - mue(k,k-1) * mue(k,k-1) ) * B(k-1);
    //    std::cerr << "While loop, step 3 y=" << y << " mue(k,k-1)=" << mue(k,k-1) << " B(k)=" << B(k) << " B(k-1)=" << B(k-1) << " test=" << test << "\n";
    while (B(k) < ( y - mue(k,k-1) * mue(k,k-1) ) * B(k-1)) {
      //      std::cerr << "Sec While loop, step 1\n";
      H.row(k).swap(H.row(k-1));
      //      std::cerr << "Sec While loop, step 2\n";
      for (int j=r+1; j<=k-2; j++)
        std::swap(gram(k,j), gram(k-1,j));
      //      std::cerr << "Sec While loop, step 3\n";
      for (int j=k+1; j<n; j++)
        std::swap(gram(j,k), gram(j,k-1));
      //      std::cerr << "Sec While loop, step 4\n";
      std::swap(gram(k-1,k-1), gram(k,k));
      //      std::cerr << "Sec While loop, step 5\n";
      for (int j=r+1; j<=k-2; j++)
        std::swap(mue(k,j), mue(k-1,j));
      //      std::cerr << "Sec While loop, step 6\n";
      Tfield mmue = mue(k,k-1);
      Tfield BB= B(k) + mmue * mmue * B(k-1);
      //      std::cerr << "Sec While loop, step 7 BB=" << BB << "\n";
      if (BB == 0) {
	B(k)   = B(k-1);
	B(k-1) = 0;
	for (int i=k+1; k<=kmax; k++)
	  std::swap(mue(i,k), mue(i,k-1));
      }
      else {
	//	std::cerr << "Sec While loop, step 7.2 B(k) = " << B(k) << " mmue=" << mmue << " k=" << k << "\n";
	if (B(k) == 0 && mmue != 0) {
	  B(k-1)= BB;
	  mue(k,k-1) = 1 / mmue;
	  for (int i=k+1; k<=kmax; k++)
	    mue(i,k-1) /= mmue;
	}
	else {
	  //	  std::cerr << "Sec While loop, step 7.2.1\n";
	  Tfield q= B(k-1) / BB;
	  //	  std::cerr << "Sec While loop, step 7.2.2\n";
	  mue(k,k-1) = mmue * q;
	  //	  std::cerr << "Sec While loop, step 7.2.3\n";
	  B(k) *= q;
	  //	  std::cerr << "Sec While loop, step 7.2.4\n";
	  B(k-1)= BB;
	  //	  std::cerr << "Sec While loop, step 7.2.5 kmax=" << kmax << " n=" << n << "\n";
	  for (int i=k+1; i<=kmax; i++) {
	    //	    std::cerr << "Sec While loop, step 7.2.5.1\n";
	    Tfield q= mue(i,k);
	    //	    std::cerr << "Sec While loop, step 7.2.5.2\n";
	    mue(i,k)= mue(i,k-1) - mmue * q;
	    //	    std::cerr << "Sec While loop, step 7.2.5.3\n";
	    mue(i,k-1)= q + mue(k,k-1) * mue(i,k);
	    //	    std::cerr << "Sec While loop, step 7.2.5.4\n";
	  }
	  //	  std::cerr << "Sec While loop, step 7.2.6\n";
	}
      }
      //      std::cerr << "Sec While loop, step 8\n";
      if (k>1) k--;
      RED( k-1 );
      //      std::cerr << "Sec While loop, step 9\n";
    }
    //    std::cerr << "While loop, step 4\n";
    if (B(r+1) == 0) r++;
    //    std::cerr << "While loop, step 5 k=" << k << " r=" << r << "\n";
    for (int l=k-2; l>=r+1; l--) {
      //      std::cerr << "  While loop, step 5.1 l=" << l << "\n";
      RED(l);
    }
    //    std::cerr << "While loop, step 6\n";
    k++;
  }
  for (int i=1; i<n; i++)
    for (int j=0; j<i; j++)
      gram(j,i) = gram(i,j);
#ifdef DEBUG
  MyMatrix<Tmat> H_T = ConvertMatrixUniversal<Tmat,Tint>(H);
  MyMatrix<Tmat> eProd = H_T * GramMat * H_T.transpose();
  if (!TestEqualityMatrix(eProd, gram)) {
    std::cerr << "The needed equality is not satisfied in LLL. DEBUG!!!!\n";
    throw TerminalException{1};
  }
#endif
  return {gram, H};
}



#endif

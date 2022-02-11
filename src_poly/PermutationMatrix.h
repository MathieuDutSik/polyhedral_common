#ifndef DEFINE_PERMUTALIB_PERMUTATION_MATRIX_H
#define DEFINE_PERMUTALIB_PERMUTATION_MATRIX_H

#include <vector>
#include <stdlib.h>
#include <string>
#include <iostream>
#include "exception.h"
#include "MAT_Matrix.h"
#include "Face_basic.h"
#include "hash_fct.h"

namespace permutalib {



template<typename Tidx_inp, typename T_impl>
struct PermutationMatrix {
public:
  using Tidx = Tidx_inp;
  using T = T_impl;
  //
  // The constructors
  //
  PermutationMatrix ()
  {
    siz=0;
    ListVal = {};
    dim = 0;
  }
  PermutationMatrix (Tidx const& _siz, int const& _dim) : siz(_siz), dim(_dim)
  {
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (n >= std::numeric_limits<Tidx>::max() - 1) {
      std::cerr << "Tidx is too small for representing the vector\n";
      throw PermutalibException{1};
    }
#endif
    ListVal = std::vector<Tidx>(siz);
    for (Tidx i=0; i<siz; i++)
      ListVal[i]=i;
    M = IdentityMat<T>(dim);
  }
  PermutationMatrix(std::vector<Tidx> && v, MyMatrix<T> && _M)
  {
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (v.size() >= std::numeric_limits<Tidx>::max() - 1) {
      std::cerr << "Tidx is too small for representing the vector\n";
      throw PermutalibException{1};
    }
#endif
    ListVal = v;
    siz = Tidx(v.size());
    M = std::move(_M);
    dim = M.rows();
  }
  PermutationMatrix(std::vector<Tidx> const& v, MyMatrix<T> const& _M)
  {
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (v.size() >= std::numeric_limits<Tidx>::max() - 1) {
      std::cerr << "Tidx is too small for representing the vector\n";
      throw PermutalibException{1};
    }
#endif
    ListVal=v;
    siz=Tidx(v.size());
    M = _M;
    dim = M.rows();
  }
  PermutationMatrix(PermutationMatrix const& ePermMatr)
  {
    siz     = ePermMatr.siz;
    ListVal = ePermMatr.ListVal;
    dim     = ePermMatr.dim;
    M       = ePermMatr.M;
  }
  PermutationMatrix(PermutationMatrix&& ePermMatr)
  {
    siz       = ePermMatr.siz;
    ListVal   = std::move(ePermMatr.ListVal);
    dim       = ePermMatr.dim;
    M         = std::move(ePermMatr.M);
    ePermMatr.siz = 0;
    ePermMatr.dim = 0;
  }
  //
  // Copy operator
  //
  PermutationMatrix<Tidx,T> operator=(PermutationMatrix<Tidx,T> const& ePermMatr)
  {
    siz     = ePermMatr.siz;
    ListVal = ePermMatr.ListVal;
    dim     = ePermMatr.dim;
    M       = ePermMatr.M;
    return *this;
  }
  PermutationMatrix<Tidx,T> operator=(PermutationMatrix<Tidx,T>&& ePermMatr)
  {
    siz     = ePermMatr.siz;
    ListVal = std::move(ePermMatr.ListVal);
    dim     = ePermMatr.dim;
    M       = std::move(ePermMatr.M);
    ePermMatr.siz = 0;
    ePermMatr.dim = 0;
    return *this;
  }
  //
  // The destructor
  //
  ~PermutationMatrix()
  {
  }
  //
  // The destructor
  //
  bool isIdentity() const
  {
    for (Tidx i=0; i<siz; i++)
      if (ListVal[i] != i)
	return false;
    for (int i=0; i<dim; i++)
      for (int j=0; j<dim; j++) {
        if (i == j && M(i,j) != 1)
          return false;
        if (i != j && M(i,j) != 0)
          return false;
      }
    return true;
  }
  Tidx at(Tidx const& i) const
  {
    return ListVal[i];
  }
  Tidx atRev(Tidx const& i) const
  {
    Tidx i1 = i, i2 = ListVal[i];
    while(true) {
      if (i2 == i)
        return i1;
      i1 = i2;
      i2 = ListVal[i2];
    }
    /*
    for (Tidx j=0; j<siz; j++)
      if (ListVal[j] == i)
        return j;
    return std::numeric_limits<Tidx>::max();
    */
  }
  const Tidx* getPtr() const
  {
    return ListVal.data();
  }
  const std::vector<Tidx>& getListVal() const
  {
    return ListVal;
  }
  const MyMatrix<T>& getMatrix() const
  {
    return M;
  }
  Tidx operator[](Tidx const& i) const
  {
    return ListVal[i];
  }
  Tidx size() const
  {
    return siz;
  }
  //
public: // Should be private in a more classic construction
  Tidx siz;
  std::vector<Tidx> ListVal;
  int dim;
  MyMatrix<T> M;
};


template<typename Tidx, typename T>
bool operator==(PermutationMatrix<Tidx,T> const& v1, PermutationMatrix<Tidx,T> const& v2)
{
  Tidx siz=v1.size();
  if (siz != v2.size() )
    return false;
  for (Tidx i=0; i<siz; i++)
    if (v1.at(i) != v2.at(i))
      return false;
  return v1.M == v2.M;
}


template<typename Tidx, typename T>
bool operator!=(PermutationMatrix<Tidx,T> const& v1, PermutationMatrix<Tidx,T> const& v2)
{
  Tidx siz=v1.size();
  if (siz != v2.size() )
    return true;
  for (Tidx i=0; i<siz; i++)
    if (v1.at(i) != v2.at(i))
      return true;
  return v1.M != v2.M;
}


template<typename Tidx, typename T>
bool operator<(PermutationMatrix<Tidx,T> const& v1, PermutationMatrix<Tidx,T> const& v2)
{
  Tidx siz1=v1.size();
  Tidx siz2=v2.size();
  if (siz1 != siz2)
    return siz1 < siz2;
  Tidx siz=siz1;
  for (Tidx i=0; i<siz; i++) {
    if (v1.at(i) != v2.at(i))
      return v1.at(i) < v2.at(i);
  }
  return v1.M < v2.M;
}

template<typename Tidx, typename T>
PermutationMatrix<Tidx,T> operator~(PermutationMatrix<Tidx,T> const& ePermMatr)
{
  Tidx siz = ePermMatr.size();
  const std::vector<Tidx>& LVal = ePermMatr.getListVal();
  std::vector<Tidx> v(siz);
  for (Tidx i=0; i<siz; i++)
    v[LVal[i]] = i;
  MyMatrix<T> Minv = Inverse(ePermMatr.M);
  return PermutationMatrix<Tidx,T>(std::move(v), std::move(Minv));
}







// Form the product v1 * v2
template<typename Tidx, typename T>
PermutationMatrix<Tidx,T> operator*(PermutationMatrix<Tidx,T> const& v1, PermutationMatrix<Tidx,T> const& v2)
{
  Tidx siz=v1.size();
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (siz != v2.size() ) {
    std::cerr << "Error in the DoubleSidedPerm product\n";
    throw PermutalibException{1};
  }
#endif
  const std::vector<Tidx>& LVal1 = v1.getListVal();
  const std::vector<Tidx>& LVal2 = v2.getListVal();
  std::vector<Tidx> vVal(siz);
  for (Tidx i=0; i<siz; i++) {
    Tidx j=LVal1[i];
    Tidx k=LVal2[j];
    vVal[i]=k;
  }
  //  return SingleSidedPerm<Tidx>(vVal);
  MyMatrix<T> Mprod = v1.M * v2.M;
  return SingleSidedPerm<Tidx>(std::move(vVal), std::move(Mprod));
}


template<typename Tidx, typename T>
void operator*=(PermutationMatrix<Tidx,T> & v1, PermutationMatrix<Tidx,T> const& v2)
{
  Tidx siz_idx = Tidx(v1.size());
  for (Tidx i=0; i<siz_idx; i++) {
    Tidx j=v1.at(i);
    Tidx k=v2.at(j);
    v1.ListVal[i]=k;
  }
  v1.M *= v2.M;
}






template<typename Tidx, typename T>
PermutationMatrix<Tidx,T> Conjugation(PermutationMatrix<Tidx,T> const& v1, PermutationMatrix<Tidx,T> const& v2)
{
  Tidx siz=v1.size();
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (siz != v2.size() ) {
    std::cerr << "Error in the DoubleSidedPerm conjugation\n";
    throw PermutalibException{1};
  }
#endif
  std::vector<Tidx> v(siz);
  for (Tidx i=0; i<siz; i++) {
    Tidx j=v1[i];
    Tidx i2=v2[i];
    Tidx j2=v2[j];
    v[i2]=j2;
  }
  MyMatrix<T> Mret = Inverse(v2.M) * v1.M * v2.M;
  return SingleSidedPerm<Tidx>(std::move(v), std::move(Mret));
}



template<typename Tidx, typename T>
Tidx PowAct(Tidx const& i, PermutationMatrix<Tidx,T> const& g)
{
  return g.at(i);
}



template<typename Tidx, typename T>
Tidx SlashAct(Tidx const& i, PermutationMatrix<Tidx,T> const& g)
{
  return g.atRev(i);
}


// LeftQuotient(a,b) = a^{-1} * b in the list.gi file
template<typename Tidx, typename T>
PermutationMatrix<Tidx,T> LeftQuotient(PermutationMatrix<Tidx,T> const& a, PermutationMatrix<Tidx,T> const& b)
{
  Tidx siz=Tidx(a.size());
  const std::vector<Tidx>& Val_A = a.getListVal();
  const std::vector<Tidx>& Val_B = b.getListVal();
  std::vector<Tidx> ListVal(siz);
  for (Tidx i=0; i<siz; i++)
    ListVal[Val_A[i]] = Val_B[i];
  MyMatrix<T> Mret = Inverse(a.M) * b.M;
  return SingleSidedPerm<Tidx>(std::move(ListVal), std::move(Mret));
}





template<typename Tidx,typename T>
PermutationMatrix<Tidx,T> Inverse(PermutationMatrix<Tidx,T> const& ePermMatr)
{
  return ~ePermMatr;
}




template<typename Tidx, typename T>
std::string GapStyleStringShift(PermutationMatrix<Tidx,T> const& ePermMatr, int const& eShift)
{
  Tidx n=ePermMatr.size();
  Face ListStat(n);
  std::string eRet;

  for (Tidx i=0; i<n; i++) {
    if (ListStat[i] == 0) {
      Tidx eFirst=i;
      Tidx eCurr=i;
      std::string ePart = "(";
      bool IsFirst=true;
      Tidx len=0;
      while(true) {
	if (!IsFirst)
	  ePart += ",";
	IsFirst=false;
	ePart += std::to_string(eCurr + eShift);
	ListStat[eCurr] = 1;
	Tidx eNext = ePermMatr.at(eCurr);
	len++;
	if (eNext == eFirst)
	  break;
	eCurr = eNext;
      }
      ePart += ")";
      if (len > 1)
	eRet += ePart;
    }
  }
  if (eRet.size() > 0)
    return eRet;
  return "()";
}

template<typename Tidx, typename T>
std::string GapStyleString(PermutationMatrix<Tidx,T> const& ePermMatr)
{
  return GapStyleStringShift(ePermMatr, 1);
}


template<typename Tidx, typename T>
std::ostream& operator<<(std::ostream& os, PermutationMatrix<Tidx,T> const& ePermMatr)
{
  os << GapStyleStringShift(ePermMatr, 1);
  return os;
}


}







namespace std {
  template<typename Tidx, typename T>
  struct hash<permutalib::PermutationMatrix<Tidx,T>>
  {
    std::size_t operator()(const permutalib::PermutationMatrix<Tidx,T> & e_val) const
    {
      uint32_t seed1 = 0x1b873540;
      const Tidx* ptr_tidx = e_val.getPtr();
      const uint8_t* ptr_i = (const uint8_t*)ptr_tidx;
      size_t len = sizeof(Tidx) * e_val.size();
      size_t seed2 = permutalib::robin_hood_hash_bytes(ptr_i, len, seed1);
      return Matrix_Hash(e_val.M, seed2);
    }
  };

}





#endif

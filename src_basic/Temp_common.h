#ifndef TEMP_COMMON_TO_ALL
#define TEMP_COMMON_TO_ALL

// Standard includes

// All the code here should have common usage with
// ---Oceanography
// ---Chemistry
// ---Mathematics
// If it applies to only one of those fields, then it should
// be out.

// Only standard includes are allowed. Forbidden are:
// ---any boost include
// ---any TBB include
// ---anything that requires a non-trivial installation

// C-style includes

#include <ctype.h>
#include <malloc.h>
#include <unistd.h>
#include <getopt.h>
#include <chrono>
#include <ctime>


#include <math.h>

#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

// STL facilities

#include <exception>
#include <vector>
#include <list>
#include <set>
#include <map>

// Functional code

#include <functional>
#include <algorithm>

// IO streams

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

// Boost serialization

#include <boost/archive/tmpdir.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>

// synonyms

typedef unsigned long ulong;
typedef unsigned int uint;

// types for exception

struct TerminalException {
  int eVal;
};

// This is guaranteed to trigger an end.
// Also it gives something that can be used for having the stacktrace via gdb.
void TerminalEnding()
{
  int val1=4;
  assert(val1 == 5);
}



typedef std::vector<std::vector<int> > VectVectInt;


//All the definitions of special fields are in other include.
//Nothing of this should depend on GMP or MPREAL or FLINT or whatever.
//
// For example, all GMP are in NumberTheory.h



// Trait definition for subset of integers

template <typename T>
struct is_implementation_of_Z {
};

template<>
struct is_implementation_of_Z<double> {
  static const bool value = false;
};

template<>
struct is_implementation_of_Z<float> {
  static const bool value = false;
};

template<>
struct is_implementation_of_Z<int> {
  static const bool value = true;
};

template<>
struct is_implementation_of_Z<long> {
  static const bool value = true;
};

// is mpreal

template <typename T>
struct is_mpreal {
  static const bool value = false;
};


// Trait definition for exactness

template <typename T>
struct is_exact_arithmetic {
};

template<>
struct is_exact_arithmetic<double> {
  static const bool value = false;
};

template<>
struct is_exact_arithmetic<float> {
  static const bool value = false;
};

// type mappings for the rings and fields:
template<typename T>
struct overlying_field {
};

template<typename T>
struct underlying_ring {
};






// Trait definition for fields.

template <typename T>
struct is_ring_field {
};

template <>
struct is_ring_field<short> {
  static const bool value = false;
};

template <>
struct is_ring_field<long> {
  static const bool value = false;
};

template <>
struct is_ring_field<int> {
  static const bool value = false;
};

template <>
struct is_ring_field<long long> {
  static const bool value = false;
};

template <>
struct is_ring_field<double> {
  static const bool value = true;
};

template <>
struct is_ring_field<float> {
  static const bool value = true;
};

// Trait of totally ordered set

template <typename T>
struct is_totally_ordered {
  static const bool value = false;
};

template <>
struct is_totally_ordered<short> {
  static const bool value = true;
};

template <>
struct is_totally_ordered<long> {
  static const bool value = true;
};

template <>
struct is_totally_ordered<int> {
  static const bool value = true;
};

template <>
struct is_totally_ordered<long long> {
  static const bool value = true;
};

template <>
struct is_totally_ordered<double> {
  static const bool value = true;
};

template <>
struct is_totally_ordered<float> {
  static const bool value = true;
};

// is double type

template <typename T>
struct is_double_type {
  static const bool value = false;
};

template <>
struct is_double_type<double> {
  static const bool value = true;
};

// is double type

template <typename T>
struct is_int_type {
  static const bool value = false;
};

template <>
struct is_int_type<int> {
  static const bool value = true;
};

// is floating arithmetic

template <typename T>
struct is_float_arithmetic {
  static const bool value = false;
};

template<>
struct is_float_arithmetic<float> {
  static const bool value = true;
};

template<>
struct is_float_arithmetic<double> {
  static const bool value = true;
};

// Trait definition for is_mpq

template <typename T>
struct is_mpq_class {
  static const bool value = false;
};

// Trait definition for is_mpz

template <typename T>
struct is_mpz_class {
  static const bool value = false;
};

// basic numeric code

int IntFloor(double const& x)
{
  return int(floor(x));
}

void GET_DOUBLE(double const& eQ, double & eD)
{
  eD=eQ;
}

template<typename T>
struct is_graphsparseimmutable_class {
  static const bool value = false;
};

template<typename T>
struct is_graphbitset_class {
  static const bool value = false;
};

template<typename T1, typename T2>
T1 UniversalTypeConversion(T2 const& a)
{
  T1 ret;
  TYPE_CONVERSION(a, ret);
  return ret;
}


template<typename T>
T T_abs(T const& eVal)
{
  if (eVal > 0)
    return eVal;
  T fVal= - eVal;
  return fVal;
}



void NearestInteger_double_int(double const& xI, int & xO)
{
  //  std::cerr << "Temp_common : NearestInteger\n";
  double xRnd_d=round(xI);
  int xRnd_z=int(xRnd_d);
  //  std::cerr << "xI=" << xI << "\n";
  auto GetErr=[&](int const& u) -> double {
    double diff = double(u) - xI;
    return T_abs(diff);
  };
  double err=GetErr(xRnd_z);
  //  std::cerr << "err=" << err << "\n";
  while(true) {
    bool IsOK=true;
    for (int i=0; i<2; i++) {
      int shift=2*i -1;
      int xTest = xRnd_z + shift;
      double TheErr=GetErr(xTest);
      //      std::cerr << "i=" << i << " shift=" << shift << " xTest=" << xTest << " TheErr=" << TheErr << "\n";
      if (TheErr < err) {
	IsOK=false;
	xRnd_z=xTest;
      }
    }
    if (IsOK)
      break;
  }
  xO=xRnd_z;
}





template<typename To>
void NearestInteger_double_To(double const& xI, To & xO)
{
  //  std::cerr << "Temp_common : NearestInteger\n";
  double xRnd_d=round(xI);
  int xRnd_i=int(xRnd_d);
  To xRnd_To=xRnd_i;
  //  std::cerr << "xI=" << xI << "\n";
  auto GetErr=[&](To const& u) -> double {
    double u_doubl;
    GET_DOUBLE(u, u_doubl);
    double diff = u_doubl - xI;
    return T_abs(diff);
  };
  double err=GetErr(xRnd_To);
  //  std::cerr << "err=" << err << "\n";
  while(true) {
    bool IsOK=true;
    for (int i=0; i<2; i++) {
      int shift=2*i -1;
      To xTest = xRnd_To + shift;
      double TheErr=GetErr(xTest);
      //      std::cerr << "i=" << i << " shift=" << shift << " xTest=" << xTest << " TheErr=" << TheErr << "\n";
      if (TheErr < err) {
	IsOK=false;
	xRnd_To=xTest;
      }
    }
    if (IsOK)
      break;
  }
  xO=xRnd_To;
}





template<typename To, typename Ti>
inline typename std::enable_if<(not is_mpreal<Ti>::value) && (not is_double_type<Ti>::value),To>::type UniversalNearestInteger(Ti const& a)
{
  To ret;
  NearestInteger(a, ret);
  //  std::cerr << "Temp_common a=" << a << " ret=" << ret << "\n";
  return ret;
}


template<typename To, typename Ti>
inline typename std::enable_if<is_double_type<Ti>::value,To>::type UniversalNearestInteger(Ti const& a)
{
  To ret;
  NearestInteger_double_To<To>(a, ret);
  //  std::cerr << "Temp_common(I) a=" << a << " ret=" << ret << "\n";
  return ret;
}






template<typename T>
T VectorSum(std::vector<T> const& eVect)
{
  T eSum=0;
  for (T eVal : eVect)
    eSum += eVal;
  return eSum;
}

template<typename T>
T VectorMin(std::vector<T> const& eVect)
{
  T eMin=eVect[0];
  for (T eVal : eVect)
    if (eVal < eMin)
      eMin=eVal;
  return eMin;
}

template<typename T>
T VectorMax(std::vector<T> const& eVect)
{
  T eMax=eVect[0];
  for (T eVal : eVect)
    if (eVal > eMax)
      eMax=eVal;
  return eMax;
}


/*
use instead std::max
template<typename T>
T T_max(T const& eVal1, T const& eVal2)  Use std::max instead
{
  if (eVal1 > eVal2)
    return eVal1;
  return eVal2;
  }*/


/*
use std::min instead
template<typename T>
T T_min(T const& eVal1, T const& eVal2)  Use std::min instead
{
  if (eVal1 > eVal2)
    return eVal2;
  return eVal1;
  }*/


template<typename T>
int T_sign(T const& eVal)
{
  if (eVal > 0)
    return 1;
  if (eVal < 0)
    return -1;
  return 0;
}

std::string random_string( size_t length )
{
  srand ( time(NULL) );
  auto randchar = []() -> char {
    const char charset[] =
    "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";
    const size_t max_index = (sizeof(charset) - 1);
    return charset[ rand() % max_index ];
  };
  std::string str(length,0);
  std::generate_n( str.begin(), length, randchar );
  return str;
}



std::string random_string_restricted( size_t length )
{
  srand ( time(NULL) );
  auto randchar = []() -> char {
    const char charset[] = "abcdefghijklmnopqrstuvwxyz";
    const size_t max_index = (sizeof(charset) - 1);
    return charset[ rand() % max_index ];
  };
  std::string str(length,0);
  std::generate_n( str.begin(), length, randchar );
  return str;
}


std::string GAP_logical(bool const& x)
{
  if (x)
    return "true";
  return "false";
}



template<typename T>
void WriteStdVectorStdVectorGAP(std::ostream & os, std::vector<std::vector<T> > const& ListVect)
{
  os << "[";
  int IsFirstVect=true;
  for (std::vector<int> const& eVect : ListVect) {
    if (!IsFirstVect)
      os << ",\n";
    IsFirstVect=false;
    os << "[";
    int siz=eVect.size();
    for (int i=0; i<siz; i++) {
      if (i>0)
        os << ",";
      os << eVect[i];
    }
    os << "]";
  }
  os << "]";
}

template<typename T>
void WriteStdVector(std::ostream& os, std::vector<T> const& V)
{
  for (auto & eVal : V)
    os << " " << eVal;
  os << "\n";
}


template<typename T>
void WriteStdVectorGAP(std::ostream& os, std::vector<T> const& V)
{
  os << "[";
  bool IsFirst=true;
  for (auto & eVal : V) {
    if (!IsFirst)
      os << ",";
    IsFirst=false;
    os << eVal;
  }
  os << "]";
}

template<typename T>
std::istream& operator>>(std::istream& is, std::vector<T>& obj)
{
  int n;
  is >> n;
  obj.resize(n);
  for (int i=0; i<n; i++) {
    T eVal;
    is >> eVal;
    obj[i]=eVal;
  }
  return is;
}
template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> const& ListVal)
{
  int n=ListVal.size();
  os << " " << n;
  for (int i=0; i<n; i++)
    os << " " << ListVal[i];
  return os;
}






template<typename T>
struct CollectedResult {
  std::vector<T> LVal;
  std::vector<int> LMult;
};

template<typename T>
CollectedResult<T> Collected(std::vector<T> const& eVect)
{
  std::set<T> SetVal;
  for (auto & eVal : eVect)
    SetVal.insert(eVal);
  std::vector<T> LVal;
  for (auto & eVal : SetVal)
    LVal.push_back(eVal);
  int eSize=LVal.size();
  std::vector<int> LMult(eSize,0);
  auto UpPosition=[&](T const& eVal) -> void {
    for (int i=0; i<eSize; i++)
      if (LVal[i] == eVal) {
	LMult[i] += 1;
	return;
      }
    std::cerr << "Should never reach that stage\n";
    throw TerminalException{1};
  };
  for (auto & eVal : eVect)
    UpPosition(eVal);
  return {LVal, LMult};
}

int NextIdx(int const& len,int const& i)
{
  if (i == len-1)
    return 0;
  return i+1;
}

int PrevIdx(int const& len,int const& i)
{
  if (i == 0)
    return len-1;
  return i-1;
}


template<typename T>
T MyPow(T const& eVal, int const& n)
{
  T eRet=1;
  for (int i=0; i<n; i++)
    eRet *= eVal;
  return eRet;
}

template<typename Tequiv>
struct EquivTest {
  bool TheReply;
  Tequiv TheEquiv;
};


std::vector<int> StdVectorFirstNentries(int const& N)
{
  std::vector<int> eList(N);
  for (int i=0; i<N; i++)
    eList[i]=i;
  return eList;
}

template<typename T>
int PositionVect(std::vector<T> const& V, T const& eVal)
{
  int len=V.size();
  for (int i=0; i<len; i++)
    if (V[i] == eVal)
      return i;
  return -1;
}

template<typename T>
bool IsVectorConstant(std::vector<T> const& V)
{
  if (V.size() == 0)
    return true;
  T eVal=V[0];
  for (auto & fVal : V) {
    if (eVal != fVal)
      return false;
  }
  return true;
}



void WriteVectorInt_GAP(std::ostream &os, std::vector<int> const& OneInc)
{
  int siz=OneInc.size();
  os << "[";
  for (int i=0; i<siz; i++) {
    if (i>0)
      os << ",";
    int eVal=OneInc[i]+1;
    os << eVal;
  }
  os << "]";
}


std::vector<int> DivideListPosition(int const& len, int const& nbBlock)
{
  std::vector<int> ListVal;
  for (int i=0; i<=nbBlock; i++) {
    double pos_d = (double(i) / double(nbBlock)) * double(len);
    int pos_i;
    NearestInteger_double_int(pos_d, pos_i);
    pos_i = std::max(0, pos_i);
    pos_i = std::min(len, pos_i);
    ListVal.push_back(pos_i);
  }
  return ListVal;
}


template<typename T>
std::vector<T> ConcatenationTwo(std::vector<T> const& L1, std::vector<T> const& L2)
{
  std::vector<T> ret = L1;
  for (auto & eVal : L2)
    ret.push_back(eVal);
  return ret;
}



/*
template<typename U>
std::vector<U> operator+(std::vector<U> const& L1, std::vector<U> const& L2)
{
  std::vector<U> ret = L1;
  for (auto & eVal : L2)
    ret.push_back(eVal);
  return ret;
}
*/


/*
  Solution using fold-expression.
  It works, but we are forced to having to use operator +. No way to use ConcatenationTwo
  as it is forbidden by the standard that only allows one of the 32 operators.
  No way to put the operator+ inside of Concatenation. Since this is explicit forbidden to
  write routine inside a routine.

template<typename... T>
auto Concatenation(T... s)
{
  return (... + s);
}
*/

template<typename T>
std::vector<T> Concatenation(std::vector<T> const& L)
{
  return L;
}


template<typename T, typename... Args>
std::vector<T> Concatenation(std::vector<T> const& first, Args... args)
{
  std::vector<T> ret = first;
  for (auto & eVal : Concatenation(args...))
    ret.push_back(eVal);
  return ret;
}


template<typename T>
int PositionProperty(std::vector<T> const& V, std::function<bool(T const&)> const& f)
{
  int len=V.size();
  for (int i=0; i<len; i++)
    if (f(V[i]))
      return i;
  return -1;
}

template<typename T>
bool ForAll(std::vector<T> const& V, std::function<bool(T const&)> const& f)
{
  for (auto & eVal : V)
    if (!f(eVal))
      return false;
  return true;
}

template<typename T>
std::vector<T> Filtered(std::vector<T> const& V, std::function<bool(T const&)> const& f)
{
  std::vector<T> LRet;
  for (auto & eVal : V)
    if (f(eVal))
      LRet.push_back(eVal);
  return LRet;
}

// Returned fields {v1,v2}
// v1 is the field returning the sorted index to the original index
// v2 is the field returning the original index to the sorted index
template<typename T>
std::pair<std::vector<int>,std::vector<int>> SortingList(std::vector<T> const & ListV)
{
  struct PairData {
    std::size_t i;
    T x;
  };
  std::size_t len=ListV.size();
  std::vector<PairData> ListPair(len);
  for (std::size_t i=0; i<len; i++) {
    PairData ePair{i, ListV[i]};
    ListPair[i]=ePair;
  }
  sort(ListPair.begin(), ListPair.end(),
       [](PairData const & x1, PairData const& x2) -> bool {
         if (x1.x < x2.x)
           return true;
         if (x2.x < x1.x)
           return false;
         return x1.i< x2.i;
       });
  std::vector<int> v1(len);
  std::vector<int> v2(len);
  for (std::size_t i=0; i<len; i++) {
    int eIdx=ListPair[i].i;
    v1[i]=eIdx;
    v2[eIdx]=i;
  }
  return {v1,v2};
}



#endif

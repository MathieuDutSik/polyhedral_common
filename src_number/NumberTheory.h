#ifndef INCLUDE_NUMBER_THEORY
#define INCLUDE_NUMBER_THEORY

#include "Temp_common.h"
#include "gmpxx.h"


// is an implementation of Z 

template<>
struct is_implementation_of_Z<mpz_class> {
  static const bool value = true;
};

template<>
struct is_implementation_of_Z<mpq_class> {
  static const bool value = false;
};



// is_euclidean_domain property

template <typename T>
struct is_euclidean_domain {
  static const bool value = false;
};

template <>
struct is_euclidean_domain<short> {
  static const bool value = true;
};

template <>
struct is_euclidean_domain<long> {
  static const bool value = true;
};

template <>
struct is_euclidean_domain<int> {
  static const bool value = true;
};

template <>
struct is_euclidean_domain<long long> {
  static const bool value = true;
};

template <>
struct is_euclidean_domain<mpz_class> {
  static const bool value = true;
};

template <>
struct is_euclidean_domain<mpq_class> {
  static const bool value = true;
};

// is_ring_field (i.e. non-zero elements are invertible)

template <>
struct is_ring_field<mpz_class> {
  static const bool value = false;
};

template <>
struct is_ring_field<mpq_class> {
  static const bool value = true;
};

// is_totally_ordered (i.e. not a complex field or ring)

template <>
struct is_totally_ordered<mpz_class> {
  static const bool value = true;
};

template <>
struct is_totally_ordered<mpq_class> {
  static const bool value = true;
};

// is exact arithmetic
// --- This excludes floating point types such as float/double
// --- This excludes limited types such as int/long

template<>
struct is_exact_arithmetic<mpz_class> {
  static const bool value = true;
};

template<>
struct is_exact_arithmetic<mpq_class> {
  static const bool value = true;
};

// is_mpq_class

template <>
struct is_mpq_class<mpq_class> {
  static const bool value = true;
};

// is_mpz_class

template <>
struct is_mpz_class<mpz_class> {
  static const bool value = true;
};

//
// Underlying ring
// For some operations, we do not need divisions
// but we need some ways to convert from one setting to another
//

template<>
struct underlying_ring<mpq_class> {
  typedef mpz_class ring_type;
};


mpz_class GetRingElement(mpq_class const& eVal)
{
  return eVal.get_num();
}


//
// Overlying field
// For some operations, we do need divisions and we need
// a canonical way to do the conversion
//

template<>
struct overlying_field<mpz_class> {
  typedef mpq_class field_type;
};

template<>
struct overlying_field<int> {
  typedef mpq_class field_type;
};

template<>
struct overlying_field<long> {
  typedef mpq_class field_type;
};

template<typename T>
struct underlying_totally_ordered_ring {
};

template<>
struct underlying_totally_ordered_ring<mpq_class> {
  typedef mpq_class real_type;
};

template<>
struct underlying_totally_ordered_ring<mpz_class> {
  typedef mpz_class real_type;
};

template<>
struct underlying_totally_ordered_ring<int> {
  typedef int real_type;
};



namespace boost { namespace serialization {

    // mpq_class
    
    template<class Archive>
      inline void load(Archive & ar,
		       mpq_class & val,
		       const unsigned int version)
    {
      std::cerr << "load(mpq_class), step 1\n";
      std::string str;
      ar & make_nvp("mpq", str);
      std::istringstream is(str);
      is >> val;
      std::cerr << "load(mpq_class), step 2\n";
    }


    
    template<class Archive>
      inline void save(Archive & ar,
		       mpq_class const& val,
		       const unsigned int version)
    {
      std::cerr << "save(mpq_class), step 1\n";
      std::ostringstream os;
      os << val;
      std::string str=os.str();
      ar & make_nvp("mpq", str);
      std::cerr << "save(mpq_class), step 2\n";
    }

    template<class Archive>
      inline void serialize(Archive & ar,
			    mpq_class & val,
			    const unsigned int version)
    {
      std::cerr << "split_free(mpq_class), step 1\n";
      split_free(ar, val, version);
      std::cerr << "split_free(mpq_class), step 2\n";
    }

    // mpz_class

    template<class Archive>
      inline void load(Archive & ar,
		       mpz_class & val,
		       const unsigned int version)
    {
      std::cerr << "load(mpz_class), step 1\n";
      std::string str;
      ar & make_nvp("mpz", str);
      std::istringstream is(str);
      is >> val;
      std::cerr << "load(mpz_class), step 2\n";
    }


    
    template<class Archive>
      inline void save(Archive & ar,
		       mpz_class const& val,
		       const unsigned int version)
    {
      std::cerr << "save(mpz_class), step 1\n";
      std::ostringstream os;
      os << val;
      std::string str=os.str();
      ar & make_nvp("mpz", str);
      std::cerr << "save(mpz_class), step 2\n";
    }

    template<class Archive>
      inline void serialize(Archive & ar,
			    mpz_class & val,
			    const unsigned int version)
    {
      std::cerr << "split_free(mpz_class), step 1\n";
      split_free(ar, val, version);
      std::cerr << "split_free(mpz_class), step 2\n";
    }
}}

template<typename T>
struct Tplusinfinity {
  Tplusinfinity(bool const& val1, T const& val2)
  {
    IsInfinity = val1;
    value = val2;
  }
  Tplusinfinity(T const& val)
  {
    IsInfinity=false;
    value = val;
  }
  void SetToInfinity()
  {
    IsInfinity=true;
  }
  Tplusinfinity<T> operator=(Tplusinfinity<T> const& x)
  {
    IsInfinity = x.IsInfinity;
    value = x.value;
  }
  Tplusinfinity<T> operator=(T const& val)
  {
    IsInfinity = false;
    value = val;
  }
  bool GetInfinity() const
  {
    return IsInfinity;
  }
  T GetValue() const
  {
    return value;
  }
private:
  bool IsInfinity;
  T value;
};


template<typename T>
bool operator==(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.GetInfinity() && y.GetInfinity())
    return true;
  if (x.GetInfinity() != y.GetInfinity())
    return false;
  return x.GetValue() == y.GetValue();
}


template<typename T>
bool operator<(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.GetInfinity() && y.GetInfinity())
    return false;
  if (!x.GetInfinity() && y.GetInfinity())
    return true;
  if (x.GetInfinity() && !y.GetInfinity())
    return false;
  return x.GetValue() < y.GetValue();
}



template<typename T>
bool operator<=(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.GetInfinity() && y.GetInfinity())
    return true;
  if (!x.GetInfinity() && y.GetInfinity())
    return true;
  if (x.GetInfinity() && !y.GetInfinity())
    return false;
  return x.GetValue() <= y.GetValue();
}


template<typename T>
bool operator>(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.GetInfinity() && y.GetInfinity())
    return false;
  if (!x.GetInfinity() && y.GetInfinity())
    return false;
  if (x.GetInfinity() && !y.GetInfinity())
    return true;
  return x.GetValue() > y.GetValue();
}


template<typename T>
Tplusinfinity<T> operator+(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.GetInfinity() || y.GetInfinity())
    return Tplusinfinity<T>(true, 0);
  return Tplusinfinity<T>(false, x.GetValue() + y.GetValue());
}

template<typename T>
Tplusinfinity<T> operator+(Tplusinfinity<T> const& x, T const& y)
{
  return Tplusinfinity<T>(x.GetInfinity(), x.GetValue() + y);
}

template<typename T>
Tplusinfinity<T> operator+(T const& x, Tplusinfinity<T> const& y)
{
  return Tplusinfinity<T>(y.GetInfinity(), x + y.GetValue());
}

template<typename T>
Tplusinfinity<T> operator*(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.GetInfinity() || y.GetInfinity())
    return Tplusinfinity<T>(true, 0);
  return Tplusinfinity<T>(false, x.GetValue() * y.GetValue());
}

template<typename T>
Tplusinfinity<T> operator*(Tplusinfinity<T> const& x, T const& y)
{
  return Tplusinfinity<T>(x.GetInfinity(), x.GetValue() * y);
}

template<typename T>
Tplusinfinity<T> operator*(T const& x, Tplusinfinity<T> const& y)
{
  return Tplusinfinity<T>(y.GetInfinity(), x * y.GetValue());
}






// Specific functions for number type.
// We basically want to do the PID (Principal Ideal Domain) case.
// What are needed for that are two functions:
// ---QuotInt for the quotient
// ---T_Norm for the norm of elements.
// Result of QuoInt(a,b) is an integer q such that q is the quotient.
// We then have a = q b + r
int QuoInt(int const& a, int const& b)
{
  int res=a % b;
  int TheQ=(a - res)/b;
  return TheQ;
}


mpz_class QuoInt(mpz_class const& a, mpz_class const& b)
{
  mpz_class res_z=a % b;
  mpz_class b_abs;
  if (b > 0)
    b_abs = b;
  else
    b_abs = -b;
  while(true) {
    if (res_z >= 0 && res_z < b_abs)
      break;
    if (res_z < 0)
      res_z += b_abs;
    if (res_z >= b_abs)
      res_z -= b_abs;
  }
  mpz_class quot_z=(a - res_z)/b;
  return quot_z;
}


mpq_class QuoInt(mpq_class const& a, mpq_class const& b)
{
  mpz_class a_den=a.get_den();
  mpz_class b_den=b.get_den();
  mpz_class eGcd;
  mpz_gcd(eGcd.get_mpz_t(), a_den.get_mpz_t(), b_den.get_mpz_t());
  mpz_class eLCM=a_den*b_den/eGcd;
  mpq_class aProd=a*eLCM;
  mpq_class bProd=b*eLCM;
  mpz_class a_num=aProd.get_num();
  mpz_class b_num=bProd.get_num();
  mpz_class b_num_pos;
  if (b_num < 0) {
    b_num_pos=-b_num;
  }
  else {
    b_num_pos=b_num;
  }
  mpz_class res_z=a_num % b_num_pos;
  while(true) {
    if (res_z >= 0 && res_z < b_num_pos)
      break;
    if (res_z < 0)
      res_z += b_num_pos;
    if (res_z >= b_num_pos)
      res_z -= b_num_pos;
    //    std::cerr << "res_z=" << res_z << " b_num_pos=" << b_num_pos << "\n";
  }
  //std::cerr << "a_num=" << a_num << " b_num=" << b_num << " res_z=" << res_z << "\n";
  mpq_class res_q=res_z;
  mpq_class quot_q=(a -res_q)/b;
  return quot_q;
}






// T_Norm should always return an integer, whatever the input type
int T_Norm(int const& eVal)
{
  return abs(eVal);
}


int T_Norm(mpq_class const& x)
{
  mpz_class eDen=x.get_den();
  if (eDen != 1) {
    std::cerr << "Denominator is not 1 as wished\n";
    std::cerr << "x=" << x << " eDen=" << eDen << "\n";
    std::cerr << "Error in T_Norm computation\n";
    throw TerminalException{1};
  }
  double x_d=x.get_d();
  int eValI=int(round(x_d));
  if (eValI > 0)
    return eValI;
  return -eValI;
}

mpq_class T_NormGen(mpq_class const& x)
{
  return T_abs(x);
}

mpz_class T_NormGen(mpz_class const& x)
{
  return T_abs(x);
}

int T_NormGen(int const& x)
{
  return abs(x);
}







bool IsInteger(mpq_class const& x)
{
  mpz_class eDen=x.get_den();
  if (eDen == 1)
    return true;
  return false;
}



mpq_class GetDenominator(mpq_class const& x)
{
  mpz_class eDen=x.get_den();
  mpq_class eDen_q=eDen;
  return eDen_q;
}

int GetDenominator(int const& x)
{
  return 1;
}

long GetDenominator(long const& x)
{
  return 1;
}

mpz_class GetDenominator(mpz_class const& x)
{
  return 1;
}





std::vector<mpz_class> FactorsInt(mpz_class const& x)
{
  if (!IsInteger(x)) {
    std::cerr << "FactorsInt requires the input to be integral\n";
    throw TerminalException{1};
  }
  std::vector<mpz_class> ListDiv;
  mpz_class TheWork=x;
  std::function<void(void)> GetOneDivisor=[&]() -> void {
    mpz_class eVal=2;
    while(true) {
      mpz_class TheRes=TheWork % eVal;
      if (TheRes == 0) {
	ListDiv.push_back(eVal);
	TheWork = TheWork / eVal;
	return;
      }
      eVal += 1;
    }
  };
  while(true) {
    if (TheWork == 1)
      break;
    GetOneDivisor();
  }
  return ListDiv;
}

std::vector<mpq_class> FactorsInt(mpq_class const& x)
{
  mpz_class x_z=x.get_num();
  std::vector<mpz_class> LFact=FactorsInt(x_z);
  std::vector<mpq_class> LFact_q;
  for (auto & eVal : LFact) {
    mpq_class eVal_q=eVal;
    LFact_q.push_back(eVal_q);
  }
  return LFact_q;
}




mpq_class GetFieldElement(mpz_class const& eVal)
{
  return eVal;
}


mpq_class GetFieldElement(long const& eVal)
{
  return eVal;
}




void GET_DOUBLE(mpq_class const& eQ, double & eD)
{
  eD=eQ.get_d();
}

void TYPE_CONVERSION(mpq_class const& a1, double & a2)
{
  a2=a1.get_d();
}


void TYPE_CONVERSION(mpq_class const& a1, mpz_class & a2)
{
  a2=a1.get_num();
  mpz_class a1_den=a1.get_den();
  if (a1_den != 1) {
    std::cerr << "a1 is not an integer\n";
    throw TerminalException{1};
  }
}




void TYPE_CONVERSION(int const& a1, mpq_class & a2)
{
  a2=a1;
}

void TYPE_CONVERSION(int const& a1, mpz_class & a2)
{
  a2=a1;
}

void TYPE_CONVERSION(mpq_class const& a1, int & a2)
{
  if (!IsInteger(a1)) {
    std::cerr << "a1=" << a1 << "\n";
    std::cerr << "Error in TYPE_CONVERSION the value is not integer\n";
    std::cerr << "no conversion possible\n";
    TerminalEnding();
    //    throw TerminalException{1};
  }
  mpz_class a1_z=a1.get_num();
  long a1_long=a1_z.get_si();
  a2=a1_long;
}

void TYPE_CONVERSION(mpz_class const& a1, mpq_class & a2)
{
  a2=a1;
}

void TYPE_CONVERSION(mpq_class const& a1, long & a2)
{
  mpz_class a1_z=a1.get_num();
  a2=a1_z.get_si();
}

void TYPE_CONVERSION(mpz_class const& a1, int & a2)
{
  long eVal_long=a1.get_si();
  a2=eVal_long;
}


//
// Nearest integer and similar stuff.
// 
mpq_class FractionalPart(mpq_class const& x)
{
  mpz_class eNum=x.get_num();
  mpz_class eDen=x.get_den();
  //  std::cerr << "FRAC eNum=" << eNum << " eDen=" << eDen << "\n";
  mpz_class res;
  mpz_mod(res.get_mpz_t(), eNum.get_mpz_t(), eDen.get_mpz_t());
  //  std::cerr << " res=" << res << "\n";
  mpq_class eRet=mpq_class(res, eDen);
  //  std::cerr << "x=" << x << " eRet=" << eRet << "\n";
  return eRet;
}

mpq_class Floor(mpq_class const& x)
{
  mpq_class eFrac=FractionalPart(x);
  return x-eFrac;
}

mpq_class Ceil(mpq_class const& x)
{
  mpq_class eFrac=FractionalPart(x);
  if (eFrac == 0)
    return x;
  return 1 + x - eFrac;
}



// return the nearest integer to x.
// If x is of the form y + 1/2 then it returns y.
mpq_class NearestInteger_rni(mpq_class const& x)
{
  mpq_class eFrac=FractionalPart(x);
  mpq_class eDiff1=eFrac;
  mpq_class eDiff2=1-eFrac;
  mpq_class RetVal=x-eFrac;
  if (eDiff1 <= eDiff2) {
    return RetVal;
  }
  else {
    return RetVal+1;
  }
}




void NearestInteger(mpq_class const& xI, mpq_class & xO)
{
  //  std::cerr << "NearestInteger mpq -> mpq\n";
  xO=NearestInteger_rni(xI);
}


void NearestInteger(mpq_class const& xI, mpz_class & xO)
{
  //  std::cerr << "NearestInteger mpq -> mpz\n";
  mpq_class xO_q=NearestInteger_rni(xI);
  xO=xO_q.get_num();
}



void NearestInteger(int const& xI, mpq_class & xO)
{
  xO=xI;
}



void NearestInteger(mpq_class const& xI, int & xO)
{
  mpz_class a1_z=xI.get_num();
  long a1_long=a1_z.get_si();
  xO=a1_long;
  auto GetPenalty=[&](int const& eInt) -> mpq_class {
    mpq_class delta=mpq_class(eInt) - xI;
    if (delta < 0)
      return -delta;
    return delta;
  };
  while(true) {
    mpq_class ePenalty=GetPenalty(xO);
    mpq_class ePenaltyP=GetPenalty(xO+1);
    mpq_class ePenaltyM=GetPenalty(xO-1);
    bool DoSomething=false;
    if (ePenaltyP < ePenalty) {
      DoSomething=true;
      xO++;
    }
    if (ePenaltyM < ePenalty) {
      DoSomething=true;
      xO--;
    }
    if (!DoSomething)
      return;
  }
}








// return the nearest integer to x.
// If x is of the form y + 1/2 then it returns y+1 and not y.
// rpi: "Rounding towards Positive Integers"
// See https://en.wikipedia.org/wiki/Floor_and_ceiling_functions#Rounding
mpq_class NearestInteger_rpi(mpq_class const& x)
{
  //  std::cerr << "--------------------------------------\n";
  mpq_class eFrac=FractionalPart(x);
  mpq_class eOne=1;
  mpq_class eTwo=2;
  //  std::cerr << "We have eOne, eTwo\n";
  mpq_class eHalf=eOne/eTwo;
  //  std::cerr << "We have eHalf\n";
  mpq_class x2=x + eHalf;
  //  std::cerr << "We have x=" << x << " eHalf=" << eHalf << " x2=" << x2 << "\n";
  mpq_class x3=Floor(x2);
  //  std::cerr << "We have x2=" << x2 << " x3=" << x3 << "\n";
  return x3;
  /*
  mpq_class residual=x3 - eHalf;
  std::cerr << "We have x3=" << x3 << " x4=" << residual << "\n";
  std::cerr << "x=" << x << " x4=" << residual << "\n";
  if (residual < -eHalf || residual >= eHalf) {
    std::cerr << "inconsistency error in residual computation\n";
    throw TerminalException{1};
  }
  mpq_class x5=x-residual;
  std::cerr << "x=" << x << " nearestInt=" << x5 << "\n";
  std::cerr << "--------------------------------------\n";
  return x5;*/
}







#endif

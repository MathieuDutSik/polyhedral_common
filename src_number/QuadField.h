#ifndef TEMP_QUAD_FIELD
#define TEMP_QUAD_FIELD

#include "Temp_common.h"
#include "NumberTheory.h"

template<typename T, int d>
class QuadField {
 private:
  T a;
  T b;
 public:
  QuadField(QuadField<T,d> const& x)  // copy constructor
    {
      a=x.a;
      b=x.b;
    }
  //  QuadField<T,d>& operator=(QuadField<T,d> const&); // assignment operator
  //  QuadField<T,d>& operator=(T const&); // assignment operator from T
  //  QuadField<T,d>& operator=(int const&); // assignment operator from T
  QuadField() // constructor from nothing
    {
      a=0;
      b=0;
    }
  QuadField(int const& u) // constructor from integer
    {
      a=u;
      b=0;
    }
  QuadField(T const& u) // constructor from T
    {
      a=u;
      b=0;
    }
  QuadField<T,d> operator=(int const& u) // assignment operator from int
    {
      QuadField<T,d> x;
      x.a=u;
      x.b=0;
      return x;
    }
  QuadField<T,d> operator=(QuadField<T,d> const& x) // assignment operator
    {
      QuadField<T,d> y;
      y.a=x.a;
      y.b=x.b;
      return y;
    }
  //
  // Arithmetic operators below:
  void operator+=(QuadField<T,d> const&x)
    {
      a += x.a;
      b += x.b;
    }
  void operator-=(QuadField<T,d> const&x)
    {
      a -= x.a;
      b -= x.b;
    }
  friend QuadField<T,d> operator+(QuadField<T,d> const&x, QuadField<T,d> const&y)
    {
      QuadField<T,d> z;
      z.a=x.a + y.a;
      z.b=x.b + y.b;
      return z;
    }
  friend QuadField<T,d> operator-(QuadField<T,d> const&x, QuadField<T,d> const&y)
    {
      QuadField<T,d> z;
      z.a=x.a - y.a;
      z.b=x.b - y.b;
      return z;
    }
  friend QuadField<T,d> operator-(QuadField<T,d> const&x, T const&y)
    {
      QuadField<T,d> z;
      z.a=x.a - y;
      z.b=x.b;
      return z;
    }
  friend QuadField<T,d> operator-(QuadField<T,d> const&x)
    {
      QuadField<T,d> z;
      z.a=- x.a;
      z.b=- x.b;
      return z;
    }
  friend QuadField<T,d> operator/(int const&x, QuadField<T,d> const&y)
    {
      T disc;
      QuadField<T,d> z;
      disc=y.a*y.a - d*y.b*y.b;
      z.a=x*y.a/disc;
      z.b=-x*y.b/disc;
      return z;
    }
  friend QuadField<T,d> operator/(QuadField<T,d> const&x, QuadField<T,d> const&y)
    {
      T disc;
      QuadField<T,d> z;
      disc=y.a*y.a - d*y.b*y.b;
      z.a=(x.a*y.a - d*x.b*y.b)/disc;
      z.b=(x.b*y.a - x.a*y.b)/disc;
      return z;
    }
  double get_d() const
  {
    T hA, hB;
    double hA_d, hB_d, doubl_d;
    hA_d=a.get_d();
    hB_d=b.get_d();
    doubl_d=(double)d;
    return hA_d + sqrt(doubl_d)*hB_d;
  }
  void operator*=(QuadField<T,d> const& x)
    {
      T hA;
      hA=a*x.a + d*  b*x.b;
      b =a*x.b + b*x.a;
      a=hA;
    }
  friend QuadField<T,d> operator*(QuadField<T,d> const&x, QuadField<T,d> const&y)
    {
      QuadField<T,d> z;
      z.a=x.a*y.a + d*  x.b*y.b;
      z.b=x.a*y.b + x.b*y.a;
      return z;
    }
  friend QuadField<T,d> operator*(int const&x, QuadField<T,d> const&y)
    {
      QuadField<T,d> z;
      z.a=x*y.a;
      z.b=x*y.b;
      return z;
    }
  friend std::ostream& operator<<(std::ostream& os, QuadField<T,d> const &v)
  {
    return os << v.a << " " << v.b;
  }
  friend std::istream& operator>>(std::istream &is, QuadField<T,d> &v)
  {
    T tmpA, tmpB;
    is >> tmpA;
    is >> tmpB;
    v.a=tmpA;
    v.b=tmpB;
    return is;
  }
  friend bool operator == (QuadField<T,d> const&x, QuadField<T,d> const&y)
  {
    if (x.a != y.a)
      return false;
    if (x.b != y.b)
      return false;
    return true;
  }
  friend bool operator != (QuadField<T,d> const&x, QuadField<T,d> const&y)
  {
    if (x.a != y.a)
      return true;
    if (x.b != y.b)
      return true;
    return false;
  }
  friend bool operator != (QuadField<T,d> const&x, int const&y)
  {
    if (x.a != y)
      return true;
    if (x.b != 0)
      return true;
    return false;
  }
  friend bool IsNonNegative(QuadField<T,d> const& x)
    {
      T disc;
      if (x.a == 0 && x.b == 0)
	return true;
      if (x.a >= 0 && x.b >= 0)
	return true;
      if (x.a<=0 && x.b <= 0)
	return false;
      disc=x.a*x.a - d*x.b*x.b;
      if (disc > 0)
	{
	  if (x.a >= 0 && x.b <= 0)
	    return true;
	  if (x.a <=0 && x.b >= 0)
	    return false;
	}
      else
	{
	  if (x.a >= 0 && x.b <= 0)
	    return false;
	  if (x.a <=0 && x.b >= 0)
	    return true;
	}
      std::cerr << "Major errors in the code\n";
      return false;
    }
  friend bool operator>=(QuadField<T,d> const& x, QuadField<T,d> const& y)
  {
    QuadField<T,d> z;
    z=x-y;
    return IsNonNegative(z);
  }
  friend bool operator <= (QuadField<T,d> const& x, QuadField<T,d> const& y)
  {
    QuadField<T,d> z;
    z=y-x;
    return IsNonNegative(z);
  }
  friend bool operator <= (QuadField<T,d> const& x, int const& y)
  {
    QuadField<T,d> z;
    z=y-x;
    return IsNonNegative(z);
  }
  friend bool operator > (QuadField<T,d> const& x, QuadField<T,d> const& y)
  {
    QuadField<T,d> z;
    z=x-y;
    if (z.a == 0 && z.b == 0)
      return false;
    return IsNonNegative(z);
  }
  friend bool operator < (QuadField<T,d> const& x, QuadField<T,d> const& y)
  {
    QuadField<T,d> z;
    z=y-x;
    if (z.a == 0 && z.b == 0)
      return false;
    return IsNonNegative(z);
  }
  friend bool operator < (QuadField<T,d> const& x, int const& y)
  {
    QuadField<T,d> z;
    z=y-x;
    if (z.a == 0 && z.b == 0)
      return false;
    return IsNonNegative(z);
  }
  
};

template<typename T, int d>
void GET_DOUBLE(QuadField<T, d> const& eQ, double & eD)
{
  eD=eQ.get_d();
}

template<typename T, int d>
struct is_totally_ordered<QuadField<T,d> > {
  static const bool value = true;
};

template<typename T, int d>
struct is_ring_field<QuadField<T,d> > {
  static const bool value = is_ring_field<T>::value;
};

template<typename T, int d>
struct is_exact_arithmetic<QuadField<T,d> > {
  static const bool value = true;
};




#endif

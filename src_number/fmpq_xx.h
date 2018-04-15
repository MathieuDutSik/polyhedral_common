#ifndef FLINT_FMPQ_CPP_CONTAINER
#define FLINT_FMPQ_CPP_CONTAINER

#include "fmpz.h"
#include "fmpq.h"

double fmpq_get_d(fmpq_t const& a)
{
  mpq_t eQ;
  mpq_init(eQ);
  fmpq_get_mpq(eQ, a);
  double x=mpq_get_d(eQ);
  mpq_clear(eQ);
  return x;
}


class fmpq_class {
 private:
  fmpq_t a;
 public:
  ~fmpq_class() {
    fmpq_clear(a);
  }
  fmpq_class(fmpq_class const& x)  // copy constructor
    {
      fmpq_init(a);
      fmpq_set(a, x.a);
    }
  //  fmpq_class<T,d>& operator=(fmpq_class<T,d> const&); // assignment operator
  //  fmpq_class<T,d>& operator=(T const&); // assignment operator from T
  //  fmpq_class<T,d>& operator=(int const&); // assignment operator from T
  fmpq_class() // constructor from nothing
    {
      fmpq_init(a);
    }
  fmpq_class(int const& u) // constructor from int
    {
      fmpq_init(a);
      fmpq_set_si(a, u, 1);
    }
  fmpq_class operator=(int const& u) // assignment operator from int
    {
      fmpq_class x;
      fmpq_set_si(x.a, u, 1);
      return x;
    }
  fmpq_class operator=(fmpq_class const& x) // assignment operator
    {
      fmpq_class y;
      fmpq_set(y.a, x.a);
      return y;
    }
  //
  // Arithmetic operators below:
  friend void operator+=(fmpq_class & x, fmpq_class const& y)
    {
      fmpq_add(x.a, x.a, y.a);
    }
  friend fmpq_class operator+(fmpq_class const&x, fmpq_class const&y)
    {
      fmpq_class z;
      fmpq_add(z.a, x.a, y.a);
      return z;
    }
  friend fmpq_class operator-(fmpq_class const&x, fmpq_class const&y)
    {
      fmpq_class z;
      fmpq_sub(z.a, x.a, y.a);
      return z;
    }
  /*
  friend fmpq_class<T,d> operator-(fmpq_class<T,d> const&x)
    {
      fmpq_class<T,d> z;
      z.a = -x.a;
      z.b = -x.b;
      return z;
      }*/
  friend const fmpq_class operator-(int const&x, fmpq_class const&y)
    {
      fmpq_class z;
      fmpq_t b;
      fmpq_init(b);
      fmpq_set_si(b, x, 1);
      fmpq_sub(z.a, b, y.a);
      fmpq_clear(b);
      return z;
    }
  friend fmpq_class operator-(fmpq_class const&x, int const&y)
    {
      fmpq_class z, b;
      b=y;
      fmpq_sub(z.a, x.a, b.a);
      return z;
    }
  friend fmpq_class operator-(fmpq_class const&x)
    {
      fmpq_class z;
      fmpq_neg(z.a, x.a);
      return z;
    }
  friend fmpq_class operator/(int const&x, fmpq_class const& y)
    {
      fmpq_class z, b;
      b=x;
      fmpq_div(z.a, b.a, y.a);
      return z;
    }
  friend fmpq_class operator/(fmpq_class const&x, fmpq_class const&y)
    {
      fmpq_class z;
      fmpq_div(z.a, x.a, y.a);
      return z;
    }
  double get_d() const
  {
    double x=fmpq_get_d(a);
    return x;
  }
  friend fmpq_class operator/(fmpq_class const&x, int const&y)
    {
      fmpq_class z, b;
      b=y;
      fmpq_div(z.a, x.a, b.a);
      return z;
    }
  friend void operator*=(fmpq_class &x, fmpq_class const& y)
    {
      fmpq_mul(x.a, x.a, y.a);
    }
  friend fmpq_class operator*(fmpq_class const&x, fmpq_class const&y)
    {
      fmpq_class z;
      fmpq_mul(z.a, x.a, y.a);
      return z;
    }
  friend fmpq_class operator*(int const&x, fmpq_class const&y)
    {
      fmpq_class z, b;
      b=x;
      fmpq_mul(z.a, b.a, y.a);
      return z;
    }
  friend std::ostream& operator<<(std::ostream& os, fmpq_class const &v)
    {
      char *str;
      str=fmpq_get_str(NULL, 10, v.a);
      std::string eStr=str;
      free(str);
      return os<<eStr;
    }
  friend std::istream& operator>>(std::istream &is, fmpq_class &v)
    {
      mpq_class eQ;
      is >> eQ;
      fmpq_set_mpq(v.a, eQ.get_mpq_t());
      return is;
    }
  friend bool operator == (fmpq_class const&x, fmpq_class const&y)
  {
    int eVal=fmpq_equal(x.a, y.a);
    if (eVal == 0)
      return true;
    return false;
  }
  friend bool operator != (fmpq_class const&x, fmpq_class const&y)
  {
    int eVal=fmpq_equal(x.a, y.a);
    if (eVal == 0)
      return false;
    return true;
  }
  friend bool operator != (fmpq_class const&x, int const&y)
  {
    fmpq_class y_f;
    y_f=y;
    int eVal=fmpq_cmp(x.a, y_f.a);
    if (eVal == 0)
      return false;
    return true;
  }
  friend bool operator>=(fmpq_class const& x, fmpq_class const& y)
  {
    int eVal=fmpq_cmp(x.a, y.a);
    // true case is x >= y
    // false case is x < y
    if (eVal < 0)
      return false;
    return true;
  }
  friend bool operator <= (fmpq_class const& x, fmpq_class const& y)
  {
    int eVal=fmpq_cmp(x.a, y.a);
    // true case is x <= y
    // false case is x > y
    if (eVal > 0)
      return false;
    return true;
  }
  friend bool operator <= (fmpq_class const& x, int const& y)
  {
    fmpq_class y_f;
    y_f=y;
    int eVal=fmpq_cmp(x.a, y_f.a);
    // true case is x <= y
    // false case is x > y
    if (eVal > 0)
      return false;
    return true;
  }
  friend bool operator > (fmpq_class const& x, fmpq_class const& y)
  {
    int eVal=fmpq_cmp(x.a, y.a);
    if (eVal > 0)
      return true;
    return false;
  }
  friend bool operator < (fmpq_class const& x, fmpq_class const& y)
  {
    int eVal=fmpq_cmp(x.a, y.a);
    if (eVal <0)
      return true;
    return false;
  }
  friend bool operator < (fmpq_class const& x, int const& y)
  {
    fmpq_class y_f;
    y_f=y;
    int eVal=fmpq_cmp(x.a, y_f.a);
    // true case is y > x
    if (eVal <0)
      return true;
    return false;
  }
};

void GET_DOUBLE(fmpq_class const& eQ, double & eD)
{
  eD=eQ.get_d();
}



#endif

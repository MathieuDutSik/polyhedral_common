class double_poly : public double {
 public:
  double h;
  void FilePrint(FILE *f)
  {
    fprintf(f, "%lg", h);
  }
  double toDouble()
  {
    return h;
  }
}


class QuadInteger_poly {
 public:
  mpq_class h1, h2;
  void FilePrint(FILE *f)
  {
    mpq_out_str(f, 10, h1.get_mpq_t());
    fprintf(f, " ");
    mpq_out_str(f, 10, h2.get_mpq_t());
  }
  void FileRead(FILE *f)
  {
    mpq_t x, y;
    mpq_init(x);
    mpq_init(y);
    mpq_inp_str(x, f, 10);
    mpq_inp_str(y, f, 10);
    h1=mpq_class(x);
    h2=mpq_class(y);
    mpq_clear(x);
    mpq_clear(y);
  }
  double toDouble()
  {
    double aX, aY, b, rVal;
    b=sqrt(5);
    aX=mpq_get_d(h1.get_mpq_t());
    aY=mpq_get_d(h2.get_mpq_t());
    rVal=aX+b*aY;
    return rVal;
  }
  QuadInteger_poly& operator=(const mpq_class u)
    {
      h1=u;
      h2=0;
      return *this;
    }
  QuadInteger_poly& operator=(const int u)
    {
      h1=u;
      h2=0;
      return *this;
    }
  QuadInteger_poly& operator+=(const QuadInteger_poly& u)
    {
      h1=h1 + u.h1;
      h2=h2 + u.h2;
      return *this;
    }
  QuadInteger_poly& operator-=(const QuadInteger_poly& u)
    {
      h1=h1 - u.h1;
      h2=h2 - u.h2;
      return *this;
    }
  QuadInteger_poly& operator*=(const QuadInteger_poly& u)
    {
      mpq_class x, y;
      x=h1*u.h1 + 5*(h2*u.h2);
      y=h1*u.h2 + h2*u.h1;
      h1=x;
      h2=y;
      return *this;
    }
  QuadInteger_poly& inv(const QuadInteger_poly& u)
    {
      mpq_class x, y, delta;
      delta=u.h1*u.h1 - 5*(u.h2*u.h2);
      x=u.h1/delta;
      y=-y/delta;
      h1=x;
      h2=y;
      return *this;
    }
  QuadInteger_poly& operator/=(const QuadInteger_poly& u)
    {
      QuadInteger_poly uInv;
      uInv=inv(u);
      *this
      mpq_class x, y;
      x=h1*u.h1 + 5*(h2*u.h2);
      y=h1*u.h2 + h2*u.h1;
      h1=x;
      h2=y;
      return *this;
    }
  
  




}

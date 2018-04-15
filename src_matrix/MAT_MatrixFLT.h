#ifndef TEMP_MATRIX_FLT
#define TEMP_MATRIX_FLT
#include <vector>
#include <stdio.h>
#include <stdlib.h>


template<typename T>
int jacobi_kernel(T **a, int n, T *d, T **v)
{
  auto ROTATE=[](T **a, int const& i, int const&j, int const& k, int const& l, T const& tau, T const& s) -> void {
    T g=a[i][j];
    T h=a[k][l];
    a[i][j]=g-s*(h+g*tau);
    a[k][l]=h+s*(g-h*tau);
  };
  int j,iq,ip,i;
  T tresh,theta,tau,t,sm,s,h,g,c;
  std::vector<T> b(n+1);
  std::vector<T> z(n+1,0);
  for (ip=1;ip<=n;ip++) {
    for (iq=1;iq<=n;iq++)
      v[ip][iq]=T(0);
    v[ip][ip]=T(1);
  }
  for (ip=1;ip<=n;ip++) {
    b[ip]=a[ip][ip];
    d[ip]=a[ip][ip];
  }
  int nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++)
	sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0)
      return nrot;
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++) {
	g=100.0*fabs(a[ip][iq]);
	if (i > 4 && T_abs(d[ip])+g == T_abs(d[ip]) && T_abs(d[iq])+g == T_abs(d[iq]))
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if (fabs(h)+g == fabs(h))
	    t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;
	  for (j=1;j<=ip-1;j++) {
	    ROTATE(a,j,ip,j,iq, tau, s);
	  }
	  for (j=ip+1;j<=iq-1;j++) {
	    ROTATE(a,ip,j,j,iq, tau, s);
	  }
	  for (j=iq+1;j<=n;j++) {
	    ROTATE(a,ip,j,iq,j, tau, s);
	  }
	  for (j=1;j<=n;j++) {
	    ROTATE(v,j,ip,j,iq, tau, s);
	  }
	  nrot++;
	}
      }
    }
    for (ip=1;ip<=n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  std::cerr << "Too many iterations in routine jacobi\n";
  throw TerminalException{1};
}

template<typename T>
void jacobi_double(MyMatrix<T> const& OrigMat, MyVector<T> & ListEigVal, MyMatrix<T> & ListEigVect)
{
  int n;
  int i;
  T **a, **v, *d;
  std::vector<T> eVect;
  int j;
  n=OrigMat.rows();
  a=new T*[n+1];
  for (i=0; i<n+1; i++)
    a[i]=new T[n+1];
  v=new T*[n+1];
  for (i=0; i<n+1; i++)
    v[i]=new T[n+1];
  d=new T[n+1];
  for (i=0; i<n+1; i++)
    d[i]=0;
  for (i=0; i<n+1; i++)
    for (j=0; j<n+1; j++)
      a[i][j]=0;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      a[i+1][j+1]=OrigMat(i,j);
  int nrot=jacobi_kernel(a, n, d, v);
  if (nrot < 0) {
    std::cerr << "Wild inconsistency\n";
    throw TerminalException{1};
  }
  for (i=0; i<n; i++)
    ListEigVal(i)=d[i+1];
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      ListEigVect(i,j)=v[j+1][i+1];
  for (i=0; i<n+1; i++)
    delete [] a[i];
  delete [] a;
  for (i=0; i<n+1; i++)
    delete [] v[i];
  delete [] v;
  delete [] d;
}

template<typename T>
std::vector<std::vector<T> > InverseSquareMatrix(std::vector<std::vector<T> > TheMat)
{
  int n, i, j, k, jSel;
  T MaxVal, hVal, eVal1, eVal2, eVal;
  std::vector<std::vector<T> > WorkMat, InvMat;
  std::vector<T> eVect;
  WorkMat=TheMat;
  n=TheMat.size();
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      if (i == j)
	eVal=1;
      else
	eVal=0;
      eVect.push_back(eVal);
    }
    InvMat.push_back(eVect);
    eVect.clear();
  }
  for (i=0; i<n; i++) {
    MaxVal=abs(WorkMat[i][i]);
    jSel=i;
    for (j=i+1; j<n; j++) {
      hVal=abs(WorkMat[i][j]);
      if (hVal > MaxVal) {
	jSel=j;
	MaxVal=hVal;
      }
    }
    if (i != jSel)
      for (k=0; k<n; k++) {
	eVal1=WorkMat[k][i];
	eVal2=WorkMat[k][jSel];
	WorkMat[k][i]=eVal1;
	WorkMat[k][jSel]=eVal2;
	eVal1=InvMat[k][i];
	eVal2=InvMat[k][jSel];
	InvMat[k][i]=eVal1;
	InvMat[k][jSel]=eVal2;
      }
    T alpha=1/WorkMat[i][i];
    for (k=0; k<n; k++) {
      WorkMat[k][i] = alpha*WorkMat[k][i];
      InvMat[k][i] = alpha*InvMat[k][i];
    }
    for (j=0; j<n; j++)
      if (i != j) {
	T alpha2=WorkMat[i][j];
	for (k=0; k<n; k++) {
	  WorkMat[k][j] = WorkMat[k][j] - WorkMat[k][i]*alpha2;
	  InvMat[k][j]  = InvMat[k][j]  - InvMat[k][i]*alpha2;
	}
      }
  }
  return InvMat;
}

template<typename T>
void PrintEigenvalues(std::ostream & os, MyMatrix<T> const& eMat)
{
  int n=eMat.rows();
  MyVector<T> ListEigVal(n);
  MyMatrix<T> ListEigVect(n,n);
  jacobi_double(eMat, ListEigVal, ListEigVect);
  for (int i=0; i<n; i++) {
    T eEig=ListEigVal(i);
    os << "i=" << i << "/" << n << " eig=" << eEig << "\n";
  }
}


#endif

#include "IterationAlgorithm.h"


#include <math.h>
/*#include "nrutil.h"*/
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

void jacobi(double **a, int n, double d[], double **v, int *nrot)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
	
	fprintf(stderr, "jacobi, step 1\n");
	if ((b = (double*)malloc((n+1)*sizeof(double))) == 0)
	  exit(EXIT_FAILURE);
	if ((z = (double*)malloc((n+1)*sizeof(double))) == 0)
	  exit(EXIT_FAILURE);
	fprintf(stderr, "jacobi, step 2\n");
	for (ip=1;ip<=n;ip++)
	  {
	    for (iq=1;iq<=n;iq++)
	      v[ip][iq]=0.0;
	    v[ip][ip]=1.0;
	  }
	fprintf(stderr, "jacobi, step 2.1\n");
	for (ip=1;ip<=n;ip++)
	  {
	    b[ip]=a[ip][ip];
	    d[ip]=a[ip][ip];
	    z[ip]=0.0;
	  }
	fprintf(stderr, "jacobi, step 3\n");
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
		  free(z);
		  free(b);
		  return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
					&& (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
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
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	fprintf(stderr, "Too many iterations in routine jacobi");
}
#undef ROTATE

int main(int argc, char *argv[])
{
  FILE *DATAINP=NULL;
  MatDOUBL OrigMat, NewMat;
  int eRet;
  int n;
  int i;
  double **a, **v, *d;
  VectDOUBL eVectReal;
  VectDOUBL ListEigVal;
  MatDOUBL ListEigVect;
  int j, nrot;
  if (argc != 2)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "StandaloneEigenvalue [DATAINP]\n");
      fprintf(stderr, "DATAINP: The input data of the form\n");
      fprintf(stderr, "n\n");
      fprintf(stderr, "a11 .... an1\n");
      fprintf(stderr, ".          .\n");
      fprintf(stderr, ".          .\n");
      fprintf(stderr, "a1n .... ann\n");
      return -1;
    }
  fprintf(stderr, "Reading input\n");
  DATAINP=fopen(argv[1], "r");
  if (DATAINP == NULL)
    {
      fprintf(stderr,"input: The file %s was not found\n",argv[1]);
      return -1;
    }
  fprintf(stderr, "  input file:%s\n", argv[1]);

  eRet=fscanf(DATAINP, "%d", &n);
  if ((a = (double**)malloc((n+1)*sizeof(double*))) == 0)
    exit(EXIT_FAILURE);
  for (i=0; i<n+1; i++)
    if ((a[i] = (double*)malloc((n+1)*sizeof(double))) == 0)
      exit(EXIT_FAILURE);
  if ((v = (double**)malloc((n+1)*sizeof(double*))) == 0)
    exit(EXIT_FAILURE);
  for (i=0; i<n+1; i++)
    if ((v[i] = (double*)malloc((n+1)*sizeof(double))) == 0)
      exit(EXIT_FAILURE);
  if ((d = (double*)malloc((n+1)*sizeof(double))) == 0)
    exit(EXIT_FAILURE);
  OrigMat=ReadMat_Doubl(DATAINP, n, n);
  for (i=0; i<n+1; i++)
    for (j=0; j<n+1; j++)
      a[i][j]=0;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      a[i+1][j+1]=OrigMat[i][j];
  nrot=0;
  fprintf(stderr, "Before jacobi work\n");
  jacobi(a, n, d, v, &nrot);
  fprintf(stderr, "After jacobi work\n");
  fprintf(stderr, "n=%d\n", n);
  for (i=0; i<n; i++)
    {
      fprintf(stderr, "i=%d eEig=%lg\n", i, d[i+1]);
    }
  fprintf(stderr, "Completion of the program\n");
  for (i=0; i<n+1; i++)
    free(a[i]);
  free(a);
  for (i=0; i<n+1; i++)
    free(v[i]);
  free(v);
  free(d);
}

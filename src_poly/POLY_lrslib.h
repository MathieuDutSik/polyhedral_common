#ifndef LRS_TEMP_H_
#define LRS_TEMP_H_

/* A templatized version of David Avis' lrs.
 * Code is mainly taken from lrslib.c
 *
 * Copyright: David Avis 2003,2006 avis@cs.mcgill.ca
 */

/* This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


#include "MAT_Matrix.h"
#include "Boost_bitset.h"
#include "MAT_MatrixInt.h"

namespace lrs {
  // some #defines and global variables from the original lrs code
  namespace globals {
    const long POS = 1L;
    const long NEG = -1L;
    const long FALSE = 0L;
    const long TRUE = 1L;
    const long MAXIMIZE = 1L;  /* maximize the lp  */
    const long MINIMIZE = 0L;  /* maximize the lp  */
    const long GE = 1L;        /* constraint is >= */
    const long EQ = 0L;        /* constraint is linearity */
    const unsigned long dict_limit = 10;
  } // end NS lrs::globals
  

template<typename T>
struct lrs_dic {
  T** A;
  long m;                     /* A has m+1 rows, row 0 is cost row            */
  long m_A;                   /* =m or m-d if nonnegative flag set            */
  long d;                     /* A has d+1 columns, col 0 is b-vector         */
  long d_orig;                /* value of d as A was allocated  (E.G.)        */
  long lexflag;               /* true if lexmin basis for this vertex         */
  long depth;                 /* depth of basis/vertex in reverse search tree */
  long i, j;                  /* last pivot row and column pivot indices      */
  T det;                 /* current determinant of basis                 */
  long *B, *Row;              /* basis, row location indices                  */
  long *C, *Col;              /* cobasis, column location indices             */
  lrs_dic *prev, *next;
};

template<typename T>
struct lrs_dat {
  long unbounded;             /* lp unbounded */
  
  long *inequality;           /* indices of inequalities corr. to cobasic ind */
    /* initially holds order used to find starting  */
    /* basis, default: m,m-1,...,2,1                */
  long *facet;                /* cobasic indices for restart in needed        */
  long *redundcol;            /* holds columns which are redundant            */
  long *linearity;            /* holds cobasic indices of input linearities   */
  long *minratio;             /* used for lexicographic ratio test            */
  long *temparray;            /* for sorting indices, dimensioned to d        */
  long inputd;                /* input dimension: n-1 for H-rep, n for V-rep  */
  
  long m;                     /* number of rows in input file                 */
  long n;                     /* number of columns in input file              */
  long lastdv;                /* index of last dec. variable after preproc    */
  /* given by inputd-nredundcol                   */
  long count[10];             /* count[0]=rays [1]=verts. [2]=base [3]=pivots */
  /* count[4]=integer vertices                    */
  long deepest;               /* max depth ever reached in search             */
  long nredundcol;            /* number of redundant columns                  */
  long nlinearity;            /* number of input linearities                  */
  long totalnodes;            /* count total number of tree nodes evaluated   */
  long runs;                  /* probes for estimate function                 */
  long seed;                  /* seed for random number generator             */
  double cest[10];            /* ests: 0=rays,1=vert,2=bases,3=vol,4=int vert */
  /**** flags  **********                         */
  long allbases;              /* globals::TRUE if all bases should be printed          */
  long bound;                 /* globals::TRUE if upper/lower bound on objective given */
  long dualdeg;               /* globals::TRUE if start dictionary is dual degenerate  */
  long frequency;             /* frequency to print cobasis indices           */
  long geometric;             /* globals::TRUE if incident vertex prints after each ray */
  long getvolume;             /* do volume calculation                        */
  long givenstart;            /* globals::TRUE if a starting cobasis is given          */
  long homogeneous;           /* globals::TRUE if all entries in column one are zero   */
  long hull;                  /* do convex hull computation if globals::TRUE           */
  long incidence;             /* print all tight inequalities (vertices/rays) */
  long lponly;                /* true if only lp solution wanted              */
  long maxdepth;              /* max depth to search to in treee              */
  long maximize;              /* flag for LP maximization                     */
  long maxoutput;             /* if positive, maximum number of output lines  */
  long minimize;              /* flag for LP minimization                     */
  long mindepth;              /* do not backtrack above mindepth              */
  long nonnegative;           /* globals::TRUE if last d constraints are nonnegativity */
  long polytope;              /* globals::TRUE for facet computation of a polytope     */
  long printcobasis;          /* globals::TRUE if all cobasis should be printed        */
  long printslack;            /* globals::TRUE if indices of slack inequal. printed    */
  long truncate;              /* globals::TRUE: truncate tree when moving from opt vert*/
  long verbose;               /* globals::FALSE for minimalist output                  */
  long restart;               /* globals::TRUE if restarting from some cobasis         */
  
  /* Variables for saving/restoring cobasis,  db */
  
  long id;                    /* numbered sequentially */
  
  long saved_count[3];        /* How often to print out current cobasis */
  long *saved_C;
  long saved_depth;
  long saved_d;
  
  long saved_flag;            /* There is something in the saved cobasis */
  
  /* Variables for cacheing dictionaries, db */
  lrs_dic<T> *Qhead, *Qtail;
};


template<typename T>
void storesign(T &a, long sa)
{
  T eProd=a*sa;
  if (eProd < 0)
    a=-a;
}

template<typename T>
long sign(T const& a)
{
  if (a < 0)
    return globals::NEG;
  if (a > 0)
    return globals::POS;
  return 0;
}

template<typename T>
long zero(T const & a)
{
  if (a == 0)
    return 1;
  return 0;
}

template<typename T>
long one(T const& a)
{
  if (a == 1)
    return 1;
  return 0;
}




template<typename T>
long lrs_greater(T const& a, T const& b)
{
  if (a > b)
    return 1;
  return 0;
}

template<typename T>
long positive(T const& a)
{
  if (a > 0)
    return 1;
  return 0;
}

template<typename T>
long negative(T const& a)
{
  if (a < 0)
    return 1;
  return 0;
}



template<typename T>
long lrs_getsolution (lrs_dic<T> * P, lrs_dat<T> * Q, T* &output, long col)
   /* check if column indexed by col in this dictionary */
   /* contains output                                   */
   /* col=0 for vertex 1....d for ray/facet             */
{
  long j;			/* cobasic index     */
  T **A = P->A;
  long *Row = P->Row;
  //  std::cerr << "col=" << col << " A[0][col]=" << A[0][col] << "\n";
  if (col == 0) {
    //    std::cerr << "col=" << col << " : lrs_getsolution, exit case 1\n";
    return lrs_getvertex (P, Q, output);
  }

  if (Q->lponly)
    {
      if (!positive(A[0][col])) {
	//	std::cerr << "col=" << col << " : lrs_getsolution, exit case 2\n";
	return globals::FALSE;
      }
    }
  else if (!negative(A[0][col])) {
    //    std::cerr << "col=" << col << " : lrs_getsolution, exit case 3\n";
    return globals::FALSE;
  }

  j = Q->lastdv + 1;
  while (j <= P->m && !negative(A[Row[j]][col]))
    j++;

  if (j <= P->m) {
    //    std::cerr << "col=" << col << " : lrs_getsolution, exit case 4\n";
    return globals::FALSE;
  }

  if (Q->geometric || Q->allbases || lexmin (P, Q, col) || Q->lponly)
    return lrs_getray (P, Q, col, Q->n, output);
  //  std::cerr << "col=" << col << " : lrs_getsolution, exit case 5\n";
  return globals::FALSE;			/* no more output in this dictionary */
}


/***********************************/
/* allocate and initialize lrs_dat */
/***********************************/
template<typename T>
lrs_dat<T> * lrs_alloc_dat ()
{
  lrs_dat<T> *Q;
  long i;
  Q = new lrs_dat<T>;

/* initialize variables */
  Q->m = 0L;
  Q->n = 0L;
  Q->inputd = 0L;
  Q->deepest = 0L;
  Q->nlinearity = 0L;
  Q->nredundcol = 0L;
  Q->runs = 0L;
  Q->seed = 1234L;
  Q->totalnodes = 0L;
  for (i = 0; i < 10; i++)
    {
      Q->count[i] = 0L;
      Q->cest[i] = 0.0;
    }
  Q->count[2] = 1L;		/* basis counter */
/* initialize flags */
  Q->allbases = globals::FALSE;
  Q->bound = globals::FALSE;            /* upper/lower bound on objective function given */
  Q->frequency = 0L;
  Q->geometric = globals::FALSE;
  Q->homogeneous = globals::TRUE;
  Q->hull = globals::FALSE;
  Q->incidence = globals::FALSE;
  Q->lponly = globals::FALSE;
  Q->maxdepth = std::numeric_limits<long>::max();
  Q->mindepth = std::numeric_limits<long>::min();
  Q->maxoutput = 0L;
  Q->nonnegative = globals::FALSE;
  Q->printcobasis = globals::FALSE;
  Q->printslack = globals::FALSE;
  Q->truncate = globals::FALSE;          /* truncate tree when moving from opt vertex        */
  Q->verbose=globals::FALSE;
  Q->maximize = globals::FALSE;		/*flag for LP maximization                          */
  Q->minimize = globals::FALSE;		/*flag for LP minimization                          */
  Q->restart = globals::FALSE;		/* globals::TRUE if restarting from some cobasis             */
  Q->givenstart = globals::FALSE;	/* globals::TRUE if a starting cobasis is given              */
  Q->saved_flag = 0;		/* no cobasis saved initially, db */
  return Q;
}


void reorder1 (long a[], long b[], long newone, long range)
/*reorder array a in increasing order with one misplaced element at index newone */
/*elements of array b are updated to stay aligned with a */
{
  long temp;
  while (newone > 0 && a[newone] < a[newone - 1])
    {
      temp = a[newone];
      a[newone] = a[newone - 1];
      a[newone - 1] = temp;
      temp = b[newone];
      b[newone] = b[newone - 1];
      b[--newone] = temp;
    }
  while (newone < range - 1 && a[newone] > a[newone + 1])
    {
      temp = a[newone];
      a[newone] = a[newone + 1];
      a[newone + 1] = temp;
      temp = b[newone];
      b[newone] = b[newone + 1];
      b[++newone] = temp;
    }
}


template<typename T>
long lrs_getfirstbasis (lrs_dic<T> ** D_p, lrs_dat<T> * Q, T** &Lin, long no_output)
/* gets first basis, globals::FALSE if none              */
/* P may get changed if lin. space Lin found    */
/* no_output is globals::TRUE supresses output headers   */
{
  //  lrs_dic<T> * D = *D_p;
  T scale, Temp;

  long i, j, k;

/* assign local variables to structures */

  T **A;
  long *B, *C, *Col;
  long *inequality;
  long *linearity;
  long hull = Q->hull;
  long m, d, lastdv, nlinearity, nredundcol;
  if (Q->lponly)
    no_output = globals::TRUE;
  m = (*D_p)->m;
  d = (*D_p)->d;
  lastdv = Q->lastdv;
  //  std::cerr << "lrs_getfirstbasis lastdv = " << lastdv << "\n";

  nredundcol = 0L;		/* will be set after getabasis        */
  nlinearity = Q->nlinearity;	/* may be reset if new linearity read */
  linearity = Q->linearity;

  A = (*D_p)->A;
  B = (*D_p)->B;
  C = (*D_p)->C;
  Col = (*D_p)->Col;
  inequality = Q->inequality;

/* default is to look for starting cobasis using linearies first, then     */
/* filling in from last rows of input as necessary                         */
/* linearity array is assumed sorted here                                  */
/* note if restart/given start inequality indices already in place         */
/* from nlinearity..d-1                                                    */
  for (i = 0; i < nlinearity; i++)	/* put linearities first in the order */
    inequality[i] = linearity[i];

  k = 0;			/* index for linearity array   */
  if (Q->givenstart)
    k = d;
  else
    k = nlinearity;
  for (i = m; i >= 1; i--)
    {
      j = 0;
      while (j < k && inequality[j] != i)
	j++;			/* see if i is in inequality  */
      if (j == k)
	inequality[k++] = i;
    }
  //  std::cerr << "Q->maximize/minimize = " << Q->maximize << " / " << Q->minimize << "\n";

/* for voronoi convert to h-description using the transform                  */
/* a_0 .. a_d-1 -> (a_0^2 + ... a_d-1 ^2)-2a_0x_0-...-2a_d-1x_d-1 + x_d >= 0 */
/* note constant term is stored in column d, and column d-1 is all ones      */
/* the other coefficients are multiplied by -2 and shifted one to the right  */
  if (!Q->maximize && !Q->minimize)
    for (j = 0; j <= d; j++)
      A[0][j]=0;

/* Now we pivot to standard form, and then find a primal feasible basis       */
/* Note these steps MUST be done, even if restarting, in order to get         */
/* the same index/inequality correspondance we had for the original prob.     */
/* The inequality array is used to give the insertion order                   */
/* and is defaulted to the last d rows when givenstart=globals::FALSE                  */

  if(Q->nonnegative) {
/* no need for initial pivots here, labelling already done */
    Q->lastdv = d;
    Q->nredundcol = 0;
  }
  else {
    if (!getabasis (*D_p, Q, inequality)) {
      //      std::cerr << "Exit case 1\n";
      return globals::FALSE;
    }
  }
  nredundcol = Q->nredundcol;
  lastdv = Q->lastdv;
  d = (*D_p)->d;
/* Reset up the inequality array to remember which index is which input inequality */
/* inequality[B[i]-lastdv] is row number of the inequality with index B[i]              */
/* inequality[C[i]-lastdv] is row number of the inequality with index C[i]              */

  for (i = 1; i <= m; i++)
    inequality[i] = i;
  if (nlinearity > 0)		/* some cobasic indices will be removed */
    {
      for (i = 0; i < nlinearity; i++)	/* remove input linearity indices */
	inequality[linearity[i]] = 0;
      k = 1;			/* counter for linearities         */
      for (i = 1; i <= m - nlinearity; i++)
	{
	  while (k <= m && inequality[k] == 0)
	    k++;		/* skip zeroes in corr. to linearity */
	  inequality[i] = inequality[k++];
	}
    }				/* end if linearity */
  if (nredundcol > 0)
    {
      Lin = new T*[nredundcol+1];
      for (i = 0; i < nredundcol; i++)
	{
	  Lin[i] = new T[Q->n+1];
	  if (!(Q->homogeneous && Q->hull && i == 0))	/* skip redund col 1 for homog. hull */
	    lrs_getray (*D_p, Q, Col[0], (*D_p)->C[0] + i - hull, Lin[i]);		/* adjust index for deletions */

	  if (!removecobasicindex (*D_p, Q, 0L)) {
	    //	    std::cerr << "Exit case 2\n";
	    return globals::FALSE;
	  }
	}
    }
  if (!primalfeasible (*D_p, Q)) {
    //    std::cerr << "Exit case 3\n";
    return globals::FALSE;
  }

/* Now solve LP if objective function was given */
  if (Q->maximize || Q->minimize)
    {
      Q->unbounded = !lrs_solvelp (*D_p, Q, Q->maximize);
      if (Q->lponly) {
	//	std::cerr << "Exit case 4\n";
	return globals::TRUE;
      }
      else                         /* check to see if objective is dual degenerate */
       {
	  j = 1;
	  while (j <= d && A[0][j] != 0)
	    j++;
	}
    }
  else
    {
      for (j = 1; j <= d; j++)
	{
	  A[0][j] = (*D_p)->det;
	  storesign (A[0][j], globals::NEG);
	}
      A[0][0]=0;
    }

  while (C[0] <= m)
    {
      i = C[0];
      j = inequality[B[i] - lastdv];
      inequality[B[i] - lastdv] = inequality[C[0] - lastdv];
      inequality[C[0] - lastdv] = j;
      C[0] = B[i];
      B[i] = i;
      reorder1 (C, Col, 0, d);
    }
/* Check to see if necessary to resize */
  if (Q->inputd > (*D_p)->d)
    *D_p = resize (*D_p, Q);

  //  std::cerr << "Exit case 6\n";
  return globals::TRUE;
}


/*****************************************/
/* getnextbasis in reverse search order  */
/*****************************************/


template<typename T>
long lrs_getnextbasis (lrs_dic<T> ** D_p, lrs_dat<T> * Q, long backtrack, unsigned long &dict_count)
	 /* gets next reverse search tree basis, globals::FALSE if none  */
	 /* switches to estimator if maxdepth set               */
	 /* backtrack globals::TRUE means backtrack from here            */
{
  /* assign local variables to structures */
  long i = 0L, j = 0L;
  long m = (*D_p)->m;
  long d = (*D_p)->d;


  if (backtrack && (*D_p)->depth == 0) {
    //    std::cerr << "lrs_getnextbasis, exit case 1\n";
    return globals::FALSE;                       /* cannot backtrack from root      */
  }

  if (Q->maxoutput > 0 && Q->count[0]+Q->count[1]-Q->hull >= Q->maxoutput) {
    //    std::cerr << "lrs_getnextbasis, exit case 2\n";
    return globals::FALSE;
  }
  //  std::cerr << "d=" << d << " m=" << m << "\n";
  //  PrintP(*D_p, "Before while loop");
  while ((j < d) || ((*D_p)->B[m] != m))
    {
      //      std::cerr << "j=" << j << " D-B[m]=" << (*D_p)->B[m] << "\n";
      if ((*D_p)->depth >= Q->maxdepth)
	{
	  backtrack = globals::TRUE;
	  if (Q->maxdepth == 0)	{ /* estimate only */
	    //	    std::cerr << "lrs_getnextbasis, exit case 3\n";
	    return globals::FALSE;
	  }
	}

      if ( Q->truncate && negative((*D_p)->A[0][0]))   /* truncate when moving from opt. vertex */
          backtrack = globals::TRUE;

      //      PrintP(*D_p, "Before backtrack test");
      if (backtrack)		/* go back to prev. dictionary, restore i,j */
	{
	  backtrack = globals::FALSE;

	  if (!check_cache (D_p, Q, &i, &j))
	    {
	      (*D_p)->depth--;
	      selectpivot (*D_p, Q, &i, &j);
	      pivot (*D_p, Q, i, j);
	      update (*D_p, Q, &i, &j);	/*Update B,C,i,j */
	    }

	  j++;			/* go to next column */
	}			/* end of if backtrack  */
      //      PrintP(*D_p, "After backtrack test");

      if ((*D_p)->depth < Q->mindepth)
	break;
      while ((j < d) && !reverse (*D_p, Q, &i, j))
	j++;

      //      PrintP(*D_p, "Before exiting test");
      if (j == d)
	backtrack = globals::TRUE;
      else
	{
	  cache_dict (D_p, Q, i, j, dict_count);
	  //	  PrintP(*D_p, "After cache_dict");
	  /* Note that the next two lines must come _after_ the 
	     call to cache_dict */

	  (*D_p)->depth++;
	  if ((*D_p)->depth > Q->deepest)
	    Q->deepest++;

	  pivot (*D_p, Q, i, j);
	  //	  PrintP(*D_p, "After pivot");
	  update (*D_p, Q, &i, &j);	/*Update B,C,i,j */
	  //	  PrintP(*D_p, "After update");

	  (*D_p)->lexflag = lexmin (*D_p, Q, 0);	/* see if lexmin basis */
	  Q->count[2]++;
	  Q->totalnodes++;

	  save_basis (*D_p, Q);
	  //	  PrintP(*D_p, "After save_basis");
	  //	  std::cerr << "lrs_getnextbasis, exit case 4\n";
	  return globals::TRUE;
	}
    }				/* end of main while loop for getnextbasis */
  //  std::cerr << "lrs_getnextbasis, exit case 5\n";
  return globals::FALSE;
}

/*************************************/
/* print out one line of output file */
/*************************************/
template<typename T>
long lrs_getvertex (lrs_dic<T> * P, lrs_dat<T> * Q, T* &output)
/*Print out current vertex if it is lexmin and return it in output */
/* return globals::FALSE if no output generated  */
{

  long i;
  long ind;			/* output index                                  */
  long ired;			/* counts number of redundant columns            */
/* assign local variables to structures */
  long *redundcol = Q->redundcol;

  long hull;
  long lexflag;

  hull = Q->hull;
  lexflag = P->lexflag;
  if (lexflag || Q->allbases)
    ++(Q->count[1]);


  if (hull)
    return globals::FALSE;		/* skip printing the origin */

  if (!lexflag && !Q->allbases && !Q->lponly)	/* not lexmin, and not printing forced */
    return globals::FALSE;


  /* copy column 0 to output */

  i = 1;
  ired = 0;
  output[0] = P->det;

  for (ind = 1; ind < Q->n; ind++)	/* extract solution */

    if ((ired < Q->nredundcol) && (redundcol[ired] == ind))
      /* column was deleted as redundant */
      {
	output[ind]=0;
	ired++;
      }
    else
      /* column not deleted as redundant */
      {
	getnextoutput (P, Q, i, 0, output[ind]);
	i++;
      }
  if (lexflag && one(output[0]))
      ++Q->count[4];               /* integer vertex */
  return globals::TRUE;
}

template<typename T>
long lrs_getray (lrs_dic<T> * P, lrs_dat<T> * Q, long col, long redcol, T* &output)
/*Print out solution in col and return it in output   */
/*redcol =n for ray/facet 0..n-1 for linearity column */
/*hull=1 implies facets will be recovered             */
/* return globals::FALSE if no output generated in column col  */
{
  long i;
  long ind;			/* output index                                  */
  long ired;			/* counts number of redundant columns            */
/* assign local variables to structures */
  long *redundcol = Q->redundcol;
  long *count = Q->count;
  long hull = Q->hull;
  long n = Q->n;
  T dettest;
  if (redcol == n)
    ++count[0];
  i = 1;
  ired = 0;

  for (ind = 0; ind < n; ind++)	/* print solution */
    {
      if (ind == 0 && !hull)	/* must have a ray, set first column to zero */
	output[0]=0;

      else if ((ired < Q->nredundcol) && (redundcol[ired] == ind))
	/* column was deleted as redundant */
	{
	  if (redcol == ind)	/* true for linearity on this cobasic index */
	    /* we print reduced determinant instead of zero */
	    output[ind] = P->det;
	  else
	    output[ind]=0;
	  ired++;
	}
      else
	/* column not deleted as redundant */
	{
	  getnextoutput (P, Q, i, col, output[ind]);
	  i++;
	}

    }
  return globals::TRUE;
}

template<typename T>
void getnextoutput(lrs_dic<T> * P, lrs_dat<T> * Q, long i, long col, T &out)
{
  long row;
  long m = P->m;
  long d = P->d;
  long lastdv = Q->lastdv;
  T **A = P->A;
  long *B = P->B;
  long *Row = P->Row;
  long j;
  row = Row[i];
  if (Q->nonnegative)
    {
      for (j = lastdv+ 1; j <= m; j++)
	{
	  if ( Q->inequality[B[j]-lastdv] == m-d+i )
	    {
	      out = A[Row[j]][col];
	      return;
	    }
	}
/* did not find inequality m-d+i in basis */
      if ( i == col )
	out=P->det;
      else
	out=0;
    }
  else
    out=A[row][col];
}

template<typename T>
int comprod(T Na, T Nb, T Nc, T Nd)
{
  if (Na*Nb > Nc*Nd)
    return 1;
  if (Na*Nb < Nc*Nd)
    return -1;
  return 0;
}


template<typename T>
long lrs_ratio (lrs_dic<T> *P, lrs_dat<T> *Q, long col)	/*find lex min. ratio */
		  /* find min index ratio -aig/ais, ais<0 */
		  /* if multiple, checks successive basis columns */
		  /* recoded Dec 1997                     */
{
  long i, j, comp, ratiocol, basicindex, start, nstart, cindex, bindex;
  long firstime;		/*For ratio test, true on first pass,else false */
  T Nmin, Dmin;
  long degencount, ndegencount;
/* assign local variables to structures */
  T** A = P->A;
  long *B = P->B;
  long *Row = P->Row;
  long *Col = P->Col;
  long *minratio = Q->minratio;
  long m, d, lastdv;

  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;


  nstart=0;
  ndegencount=0;
  degencount = 0;
  for (j = lastdv + 1; j <= m; j++)
    {
      /* search rows with negative coefficient in dictionary */
      /*  minratio contains indices of min ratio cols        */
      if (negative(A[Row[j]][col]))
	minratio[degencount++] = j;
    }				/* end of for loop */
  if (degencount == 0)
    return (degencount);	/* non-negative pivot column */

  ratiocol = 0;			/* column being checked, initially rhs */
  start = 0;			/* starting location in minratio array */
  bindex = d + 1;		/* index of next basic variable to consider */
  cindex = 0;			/* index of next cobasic variable to consider */
  basicindex = d;		/* index of basis inverse for current ratio test, except d=rhs test */
  while (degencount > 1)	/*keep going until unique min ratio found */
    {
      if (B[bindex] == basicindex)	/* identity col in basis inverse */
	{
	  if (minratio[start] == bindex)
	    /* remove this index, all others stay */
	    {
	      start++;
	      degencount--;
	    }
	  bindex++;
	}
      else
	/* perform ratio test on rhs or column of basis inverse */
	{
	  firstime = globals::TRUE;
	  /*get next ratio column and increment cindex */
	  if (basicindex != d)
	    ratiocol = Col[cindex++];
	  for (j = start; j < start + degencount; j++)
	    {
	      i = Row[minratio[j]];	/* i is the row location of the next basic variable */
	      comp = 1;		/* 1:  lhs>rhs;  0:lhs=rhs; -1: lhs<rhs */
	      if (firstime)
		firstime = globals::FALSE;	/*force new min ratio on first time */
	      else
		{
		  if (positive(Nmin) || negative(A[i][ratiocol]))
		    {
		      if (negative(Nmin) || positive(A[i][ratiocol]))
			comp = comprod (Nmin, A[i][col], A[i][ratiocol], Dmin);
		      else
			comp = -1;
		    }

		  else if (zero(Nmin) && zero(A[i][ratiocol]))
		    comp = 0;

		  if (ratiocol == 0)
		    comp = -comp;	/* all signs reversed for rhs */
		}
	      if (comp == 1)
		{		/*new minimum ratio */
		  nstart = j;
		  Nmin = A[i][ratiocol];
		  Dmin = A[i][col];
		  ndegencount = 1;
		}
	      else if (comp == 0)	/* repeated minimum */
		minratio[nstart + ndegencount++] = minratio[j];

	    }			/* end of  for (j=start.... */
	  degencount = ndegencount;
	  start = nstart;
	}			/* end of else perform ratio test statement */
      basicindex++;		/* increment column of basis inverse to check next */
    }				/*end of while loop */
  return (minratio[start]);
}				/* end of ratio */





template<typename T>
long reverse(lrs_dic<T> *P, lrs_dat<T> *Q, long *r, long s)
/*  find reverse indices  */
/* globals::TRUE if B[*r] C[s] is a reverse lexicographic pivot */
{
  long i, j, row, col;

/* assign local variables to structures */
  T** A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long d = P->d;

  col = Col[s];
  if (!negative(A[0][col]))
    return globals::FALSE;

  *r = lrs_ratio<T>(P, Q, col);
  if (*r == 0)			/* we have a ray */
    return (globals::FALSE);

  row = Row[*r];

/* check cost row after "pivot" for smaller leaving index    */
/* ie. j s.t.  A[0][j]*A[row][col] < A[0][col]*A[row][j]     */
/* note both A[row][col] and A[0][col] are negative          */

  for (i = 0; i < d && C[i] < B[*r]; i++)
    if (i != s)
      {
	j = Col[i];
	if (positive(A[0][j]) || negative(A[row][j]))		/*or else sign test fails trivially */
	  if ((!negative(A[0][j]) && !positive(A[row][j])) ||
	      comprod (A[0][j], A[row][col], A[0][col], A[row][j]) == -1)
	    return (globals::FALSE);
      }
  return (globals::TRUE);
}

template<typename T>
long selectpivot (lrs_dic<T> *P, lrs_dat<T> *Q, long *r, long *s)
/* select pivot indices using lexicographic rule   */
/* returns globals::TRUE if pivot found else globals::FALSE          */
/* pivot variables are B[*r] C[*s] in locations Row[*r] Col[*s] */
{
  long j, col;
/* assign local variables to structures */
  T** A = P->A;
  long *Col = P->Col;
  long d = P->d;

  *r = 0;
  *s = d;
  j = 0;
  while ((j < d) && (!positive(A[0][Col[j]])))
    j++;

  if (j < d) {
    *s = j;
    col = Col[j];
    *r = lrs_ratio<T>(P, Q, col);
    if (*r != 0)
      return (globals::TRUE);
  }
  return (globals::FALSE);
}

template<typename T>
void  pivot (lrs_dic<T> * P, lrs_dat<T> * Q, long bas, long cob)
		     /* Qpivot routine for array A              */
		     /* indices bas, cob are for Basis B and CoBasis C    */
		     /* corresponding to row Row[bas] and column       */
		     /* Col[cob]   respectively                       */
{
  long r, s;
  long i, j;
  T Ns, Nt, Ars;
/* assign local variables to structures */

  T** A = P->A;
  long *Row = P->Row;
  long *Col = P->Col;
  long d, m_A;

  d = P->d;
  m_A = P->m_A;
  Q->count[3]++;    /* count the pivot */

  r = Row[bas];
  s = Col[cob];

/* Ars=A[r][s]    */
  Ars = A[r][s];
  storesign (P->det, sign (Ars));	/*adjust determinant to new sign */


  for (i = 0; i <= m_A; i++)
    if (i != r)
      for (j = 0; j <= d; j++)
	if (j != s)
	  {
/*          A[i][j]=(A[i][j]*Ars-A[i][s]*A[r][j])/P->det; */
	    A[i][j]=(A[i][j]*Ars-A[i][s]*A[r][j])/P->det;
	  }

  if (sign (Ars) == globals::POS)
    {
      for (j = 0; j <= d; j++)	/* no need to change sign if Ars neg */
	/*   A[r][j]=-A[r][j];              */
	if (!zero(A[r][j]))
	  A[r][j]=-A[r][j];
    }
  else
    for (i = 0; i <= m_A; i++)
      if (!zero(A[i][s]))
	A[i][s]=-A[i][s];

  A[r][s] = P->det;
  P->det = Ars;
  storesign (P->det, globals::POS);
}



template<typename T>
long primalfeasible (lrs_dic<T> * P, lrs_dat<T> * Q)
/* Do dual pivots to get primal feasibility */
/* Note that cost row is all zero, so no ratio test needed for Dual Bland's rule */
{
  long primalinfeasible = globals::TRUE;
  long i, j;
/* assign local variables to structures */
  T** A = P->A;
  long *Row = P->Row;
  long *Col = P->Col;
  long m, d, lastdv;
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;
  while (primalinfeasible)
    {
      i=lastdv+1;
      //      std::cerr << "i=" << i << " Row=" << Row[i] << "\n";
      //      std::cerr << "val=" << A[Row[i]][0] << "\n";
      while (i <= m && !negative(A[Row[i]][0]))
        i++;
      if (i <= m )
	{
	  j = 0;		/*find a positive entry for in row */
	  while (j < d && !positive(A[Row[i]][Col[j]])) {
	    //	    std::cerr << "j=" << j << " A[R][C]=" << A[Row[i]][Col[j]] << "\n";
	    j++;
	  }
	  if (j >= d)
	    return (globals::FALSE);	/* no positive entry */
	  pivot (P, Q, i, j);
	  update (P, Q, &i, &j);
        }
      else
         primalinfeasible = globals::FALSE;
    }				/* end of while primalinfeasibile */
  return (globals::TRUE);
}				/* end of primalfeasible */


template<typename T>
long lrs_solvelp (lrs_dic<T> * P, lrs_dat<T> * Q, long maximize)
/* Solve primal feasible lp by Dantzig`s rule and lexicographic ratio test */
/* return globals::TRUE if bounded, globals::FALSE if unbounded                              */
{
  long i, j;
/* assign local variables to structures */
  long d = P->d;

  while (dan_selectpivot (P, Q, &i, &j))
    {
      Q->count[3]++;
      pivot (P, Q, i, j);
      update (P, Q, &i, &j);	/*Update B,C,i,j */
    }

  if (j < d && i == 0)		/* selectpivot gives information on unbounded solution */
    {
      return globals::FALSE;
    }
  return globals::TRUE;
}				/* end of lrs_solvelp  */

template<typename T>
long getabasis (lrs_dic<T> * P, lrs_dat<T> * Q, long order[])
/* Pivot Ax<=b to standard form */
/*Try to find a starting basis by pivoting in the variables x[1]..x[d]        */
/*If there are any input linearities, these appear first in order[]           */
/* Steps: (a) Try to pivot out basic variables using order                    */
/*            Stop if some linearity cannot be made to leave basis            */
/*        (b) Permanently remove the cobasic indices of linearities           */
/*        (c) If some decision variable cobasic, it is a linearity,           */
/*            and will be removed.                                            */
{
  long i, j, k;
/* assign local variables to structures */
  T** A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long *linearity = Q->linearity;
  long *redundcol = Q->redundcol;
  long m, d, nlinearity;
  long nredundcol = 0L;		/* will be calculated here */
  nlinearity = Q->nlinearity;
  m = P->m;
  d = P->d;

  for (j = 0; j < m; j++)
    {
      i = 0;
      while (i <= m && B[i] != d + order[j])
	i++; /* find leaving basis index i */
      if (j < nlinearity && i > m)	/* cannot pivot linearity to cobasis */
	return globals::FALSE;
      if (i <= m)
	{			/* try to do a pivot */
	  k = 0;
	  while (C[k] <= d && zero(A[Row[i]][Col[k]]))
	    k++;
	  if (C[k] <= d)
	    {
	      pivot (P, Q, i, k);
	      update (P, Q, &i, &k);
	    }
	  else if (j < nlinearity)
	    {			/* cannot pivot linearity to cobasis */
	      if (zero(A[Row[i]][0]))
		linearity[j] = 0;
	      else
		return globals::FALSE;
	    }			/* end if j < nlinearity */
	}			/* end of if i <= m .... */
    }				/* end of for   */
/* update linearity array to get rid of redundancies */
  i = 0;
  k = 0;			/* counters for linearities         */
  while (k < nlinearity)
    {
      while (k < nlinearity && linearity[k] == 0)
	k++;
      if (k < nlinearity)
	linearity[i++] = linearity[k++];
    }
  nlinearity = i;
/* column dependencies now can be recorded  */
/* redundcol contains input column number 0..n-1 where redundancy is */
  k = 0;
  while (k < d && C[k] <= d) {
    if (C[k] <= d)		/* decision variable still in cobasis */
      redundcol[nredundcol++] = C[k] - Q->hull;	/* adjust for hull indices */
    k++;
  }
/* now we know how many decision variables remain in problem */
  Q->nredundcol = nredundcol;
  Q->lastdv = d - nredundcol;

/* Remove linearities from cobasis for rest of computation */
/* This is done in order so indexing is not screwed up */

  for (i = 0; i < nlinearity; i++)
    {				/* find cobasic index */
      k = 0;
      while (k < d && C[k] != linearity[i] + d)
	k++;
      if (k >= d)
	{
	  //	  std::cerr << "Error removing linearity\n";
	  return globals::FALSE;
	}
      if (!removecobasicindex (P, Q, k))
	return globals::FALSE;
      d = P->d;
    }

/* Check feasability */
  if (Q->givenstart)
    {
      i = Q->lastdv + 1;
      while (i <= m && !negative(A[Row[i]][0]))
	i++;
      if (i <= m)
	fprintf (stderr, "\n*Infeasible startingcobasis - will be modified");
    }
  return globals::TRUE;
}				/*  end of getabasis */

template<typename T>
long removecobasicindex (lrs_dic<T> * P, lrs_dat<T> * Q, long k)
/* remove the variable C[k] from the problem */
/* used after detecting column dependency    */
{
  long i, j, cindex, deloc;
/* assign local variables to structures */
  T** A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Col = P->Col;
  long m, d;
  m = P->m;
  d = P->d;

  cindex = C[k];		/* cobasic index to remove              */
  deloc = Col[k];		/* matrix column location to remove     */

  for (i = 1; i <= m; i++)	/* reduce basic indices by 1 after index */
    if (B[i] > cindex)
      B[i]--;
  
  for (j = k; j < d; j++)	/* move down other cobasic variables    */
    {
      C[j] = C[j + 1] - 1;	/* cobasic index reduced by 1           */
      Col[j] = Col[j + 1];
    }

  //  std::cerr << "deloc=" << deloc << "\n";
  if (deloc != d)               
    {
  /* copy col d to deloc */
      for (i = 0; i <= m; i++)
        A[i][deloc] = A[i][d];

  /* reassign location for moved column */
      j = 0;
      while (Col[j] != d)
        j++;
      
    Col[j] = deloc;
  }

  P->d--;
  return globals::TRUE;
}				/* end of removecobasicindex */

template<typename T>
lrs_dic<T> *new_lrs_dic (long m, long d, long m_A)
{
  lrs_dic<T> *p;
  int i;

  p = new lrs_dic<T>;

  p->B   = new long[m+1];
  p->Row = new long[m+1];
  p->C   = new long[d+1];
  p->Col = new long[d+1];

  p->d_orig=d;

  p->A=new T*[m_A+1];
  for (i=0; i<=m_A; i++)
    p->A[i]=new T[d+1];
  return p;
}



template<typename T>
lrs_dic<T> *resize (lrs_dic<T> * P, lrs_dat<T> * Q)
	/* resize the dictionary after some columns are deleted, ie. inputd>d */
	/* a new lrs_dic record is created with reduced size, and items copied over */
{
  lrs_dic<T> *P1;			/* to hold new dictionary in case of resizing */

  long i, j;
  long m, d, m_A;


  m = P->m;
  d = P->d;
  m_A = P->m_A;

/* get new dictionary record */

  P1 = new_lrs_dic<T> (m, d, m_A);

/* copy data from P to P1    */
  P1->i = P->i;
  P1->j = P->j;
  P1->depth = P->depth;
  P1->m = P->m;
  P1->d = P1->d_orig = d;
  P1->lexflag = P->lexflag;
  P1->m_A = P->m_A;
  P1->det = P->det;

  for (i = 0; i <= m; i++)
    {
      P1->B[i] = P->B[i];
      P1->Row[i] = P->Row[i];
    }
  for (i = 0; i <= m_A; i++)
    {
      for (j = 0; j <= d; j++)
	P1->A[i][j] = P->A[i][j];
    }


  for (j = 0; j <= d; j++)
    {
      P1->Col[j] = P->Col[j];
      P1->C[j] = P->C[j];
    }

  lrs_free_dic (P,Q);

/* Reassign cache pointers */

  Q->Qhead = P1;
  Q->Qtail = P1;
  P1->next = P1;
  P1->prev = P1;

  return P1;

}


template<typename T>
long restartpivots (lrs_dic<T> * P, lrs_dat<T> * Q)
/* facet contains a list of the inequalities in the cobasis for the restart */
/* inequality contains the relabelled inequalities after initialization     */
{
  long i, j, k;
  long *Cobasic;		/* when restarting, Cobasic[j]=1 if j is in cobasis */
/* assign local variables to structures */
  T** A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long *inequality = Q->inequality;
  long *facet = Q->facet;
  long nlinearity = Q->nlinearity;
  long m, d, lastdv;
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

  Cobasic = new long[m + d + 2];

  /* set Cobasic flags */
  for (i = 0; i < m + d + 1; i++)
    Cobasic[i] = 0;
  for (i = 0; i < d; i++)	/* find index corresponding to facet[i] */
    {
      j = 1;
      while (facet[i + nlinearity] != inequality[j])
	j++;
      Cobasic[j + lastdv] = 1;
    }

  /* Note that the order of doing the pivots is important, as */
  /* the B and C vectors are reordered after each pivot       */

/* Suggested new code from db starts */
  i=m;
  while (i>d){
    while(Cobasic[B[i]]){
      k = d - 1;
      while (k >= 0 && (zero(A[Row[i]][Col[k]]) || Cobasic[C[k]])){
       k--;
      }
      if (k >= 0)  {
           /*db asks: should i really be modified here? (see old code) */
           /*da replies: modifying i only makes is larger, and so      */
           /*the second while loop will put it back where it was       */
           /*faster (and safer) as done below                          */
       long  ii=i;
       pivot (P, Q, ii, k);
       update (P, Q, &ii, &k);
      } else {
	delete [] Cobasic;
       return globals::FALSE;
      }
    }
    i--;
  }
/* Suggested new code from db ends */

  if (lexmin (P, Q, 0))
    --Q->count[1];		/* decrement vertex count if lexmin */
/* check restarting from a primal feasible dictionary               */
  for (i = lastdv + 1; i <= m; i++)
    if (negative(A[Row[i]][0]))
      {
        delete [] Cobasic;
	return globals::FALSE;
      }
  delete [] Cobasic;
  return globals::TRUE;
}



template<typename T>
long lexmin (lrs_dic<T> * P, lrs_dat<T> * Q, long col)
  /*test if basis is lex-min for vertex or ray, if so globals::TRUE */
  /* globals::FALSE if a_r,g=0, a_rs !=0, r > s          */
{
/*do lexmin test for vertex if col=0, otherwise for ray */
  long r, s, i, j;
/* assign local variables to structures */
  T** A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long m = P->m;
  long d = P->d;

  for (i = Q->lastdv + 1; i <= m; i++)
    {
      r = Row[i];
      if (zero(A[r][col]))	/* necessary for lexmin to fail */
	for (j = 0; j < d; j++)
	  {
	    s = Col[j];
	    if (B[i] > C[j])	/* possible pivot to reduce basis */
	      {
		if (zero (A[r][0]))	/* no need for ratio test, any pivot feasible */
		  {
		    if (!zero (A[r][s]))
		      return (globals::FALSE);
		  }
		else if (negative(A[r][s]) && ismin (P, Q, r, s))
		  {
		    return (globals::FALSE);
		  }
	      }			/* end of if B[i] ... */
	  }
    }
  return (globals::TRUE);
}

template<typename T>
long ismin (lrs_dic<T> * P, lrs_dat<T> * Q, long r, long s)
/*test if A[r][s] is a min ratio for col s */
{
  long i;
/* assign local variables to structures */
  T** A = P->A;
  long m_A = P->m_A;

  for (i = 1; i <= m_A; i++)
    if ((i != r) &&
	negative(A[i][s]) && comprod (A[i][0], A[r][s], A[i][s], A[r][0]))
      {
	return (globals::FALSE);
      }

  return (globals::TRUE);
}

template<typename T>
void update (lrs_dic<T> * P, lrs_dat<T> * Q, long *i, long *j)
 /*update the B,C arrays after a pivot */
 /*   involving B[bas] and C[cob]           */
{

  long leave, enter;
/* assign local variables to structures */
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long m = P->m;
  long d = P->d;

  leave = B[*i];
  enter = C[*j];
  B[*i] = enter;
  reorder1 (B, Row, *i, m + 1);
  C[*j] = leave;
  reorder1 (C, Col, *j, d);
/* restore i and j to new positions in basis */
  for (*i = 1; B[*i] != enter; (*i)++);		/*Find basis index */
  for (*j = 0; C[*j] != leave; (*j)++);		/*Find co-basis index */
}				/* end of update */

template<typename T>
long lrs_degenerate (lrs_dic<T> * P, lrs_dat<T> * Q)
/* globals::TRUE if the current dictionary is primal degenerate */
/* not thoroughly tested   2000/02/15                  */
{
  long i;
  long *B, *Row;

  T** A = P->A;
  long d = P->d;
  long m = P->m;

  B = P->B;
  Row = P->Row;

  for (i = d + 1; i <= m; i++)
    if (zero (A[Row[i]][0]))
      return globals::TRUE;

  return globals::FALSE;
}


/*********************************************************/
/*                 Miscellaneous                         */
/******************************************************* */

void 
reorder (long a[], long range)
/*reorder array in increasing order with one misplaced element */
{
  long i, temp;
  for (i = 0; i < range - 1; i++)
    if (a[i] > a[i + 1])
      {
	temp = a[i];
	a[i] = a[i + 1];
	a[i + 1] = temp;
      }
  for (i = range - 2; i >= 0; i--)
    if (a[i] > a[i + 1])
      {
	temp = a[i];
	a[i] = a[i + 1];
	a[i + 1] = temp;
      }

}				/* end of reorder */


template<typename T>
long checkredund (lrs_dic<T> * P, lrs_dat<T> * Q)
/* Solve primal feasible lp by least subscript and lex min basis method */
/* to check redundancy of a row in objective function                   */
/* returns globals::TRUE if redundant, else globals::FALSE                                */
{
  T Ns, Nt;
  long i, j;
  long r, s;

/* assign local variables to structures */
  T** A = P->A;
  long *B, *C, *Row, *Col;
  long d = P->d;

  B = P->B;
  C = P->C;
  Row = P->Row;
  Col = P->Col;

  while (selectpivot (P, Q, &i, &j))
    {
      Q->count[2]++;

/* sign of new value of A[0][0]            */
/* is      A[0][s]*A[r][0]-A[0][0]*A[r][s] */

      r = Row[i];
      s = Col[j];
      Ns=A[0][s] * A[r][0];
      Nt=A[0][0] * A[r][s];
      if (lrs_greater(Ns,Nt))
	return globals::FALSE;		/* non-redundant */

      pivot (P, Q, i, j);
      update (P, Q, &i, &j);	/*Update B,C,i,j */
    }
  return !(j < d && i == 0);	/* unbounded is also non-redundant */
}				/* end of checkredund  */

template<typename T>
long checkcobasic (lrs_dic<T> * P, lrs_dat<T> * Q, long index)
/* globals::TRUE if index is cobasic and nonredundant                         */
/* globals::FALSE if basic, or degen. cobasic, where it will get pivoted out  */

{

/* assign local variables to structures */

  T** A = P->A;
  long *B, *C, *Row, *Col;
  long d = P->d;
  long m = P->m;
  long i = 0;
  long j = 0;
  long s;

  B = P->B;
  C = P->C;
  Row = P->Row;
  Col = P->Col;


  while ((j < d) && C[j] != index)
    j++;

  if (j == d)
    return globals::FALSE;		/* not cobasic index */


/* index is cobasic */

  s = Col[j];
  i = Q->lastdv + 1;

  while ((i <= m) &&
	 (zero (A[Row[i]][s]) || !zero (A[Row[i]][0])))
    i++;

  if (i > m)
    {
      return globals::TRUE;
    }

  pivot (P, Q, i, j);
  update (P, Q, &i, &j);	/*Update B,C,i,j */

  return globals::FALSE;			/*index is no longer cobasic */

}				/* end of checkcobasic */

template<typename T>
long checkindex (lrs_dic<T> * P, lrs_dat<T> * Q, long index)
/* 0 if index is non-redundant inequality */
/* 1 if index is redundant     inequality */
/* 2 if index is input linearity          */
/*NOTE: row is returned all zero if redundant!! */
{
  long i, j;

  T** A = P->A;
  long *Row = P->Row;
  long *B = P->B;
  long d = P->d;
  long m = P->m;


/* each slack index must be checked for redundancy */
/* if in cobasis, it is pivoted out if degenerate */
/* else it is non-redundant                       */

  if (checkcobasic (P, Q, index))
    return 0;

/* index is basic   */
  j = 1;
  while ((j <= m) && (B[j] != index))
    j++;

  i = Row[j];

  /* copy row i to cost row, and set it to zero */

  for (j = 0; j <= d; j++)
    {
      A[0][j] = A[i][j];
      A[0][j] = - A[0][j];
      A[i][j]=0;
    }


  if (checkredund (P, Q))
    return 1L;

/* non-redundant, copy back and change sign */

  for (j = 0; j <= d; j++)
    A[i][j]=-A[0][j];
  return 0;
}				/* end of checkindex */

/***************************************************************/
/*                                                             */
/*     Routines for caching, allocating etc.                   */
/*                                                             */
/***************************************************************/

/* From here mostly Bremner's handiwork */

template<typename T>
void cache_dict (lrs_dic<T> ** D_p, lrs_dat<T> * global, long i, long j, unsigned long &dict_count)
{
  if (globals::dict_limit > 1)
    {
      (*D_p)->i = i;
      (*D_p)->j = j;
      pushQ (global, (*D_p)->m, (*D_p)->d, (*D_p)->m_A, dict_count);
      copy_dict (global, global->Qtail, *D_p);
    }
  *D_p = global->Qtail;
}

template<typename T>
void copy_dict (lrs_dat<T> * global, lrs_dic<T> * dest, lrs_dic<T> * src)
{
  long m = src->m;
  long m_A = src->m_A;        /* number of rows in A */
  long d = src->d;
  long r,s;

  for ( r=0;r<=m_A;r++)
    for( s=0;s<=d;s++)
       dest->A[r][s]=src->A[r][s];

  dest->i = src->i;
  dest->j = src->j;
  dest->m = m;
  dest->d = d;
  dest->m_A  = src->m_A;

  dest->depth = src->depth;
  dest->lexflag = src->lexflag;

  dest->det = src->det;

  for (int u=0; u<m+1; u++) {
    dest->B[u] = src->B[u];
    dest->Row[u] = src->Row[u];
  }
  for (int u=0; u<d+1; u++) {
    dest->C[u] = src->C[u];
    dest->Col[u] = src->Col[u];
  }
}

/* 
 * pushQ(lrs_dat *globals,m,d):
 * this routine ensures that Qtail points to a record that 
 * may be copied into.
 *
 * It may create a new record, or it may just move the head pointer
 * forward so that know that the old record has been overwritten.
 */
template<typename T>
void pushQ (lrs_dat<T> * global, long m, long d ,long m_A, unsigned long & dict_count)
{
  if ((global->Qtail->next) == global->Qhead)
    {
      if (dict_count < globals::dict_limit)
	{
	  lrs_dic<T> *p;
	  p = new_lrs_dic<T> (m, d, m_A);
	  if (p)
	    {
	      p->next = global->Qtail->next;
	      (global->Qtail->next)->prev = p;
	      (global->Qtail->next) = p;
	      p->prev = global->Qtail;
	      dict_count++;
	      global->Qtail = p;
	    }
	  else
	    {
	      /* virtual memory exhausted. bummer */
	      global->Qhead = global->Qhead->next;
	      global->Qtail = global->Qtail->next;
	    }
	}
      else
	{
	  /*
	   * user defined limit reached. start overwriting the
	   * beginning of Q 
	   */
	  global->Qhead = global->Qhead->next;
	  global->Qtail = global->Qtail->next;
	}
    }

  else
    {
      global->Qtail = global->Qtail->next;
    }
}

template<typename T>
lrs_dic<T> *lrs_getdic(lrs_dat<T> *Q)
{
  lrs_dic<T> *p;
  long m;
  m = Q->m;
  if(Q->nonnegative)
    m = m+Q->inputd;
  p = new_lrs_dic<T> (m, Q->inputd, Q->m);
  if (!p)
    return nullptr;

  p->next = p;
  p->prev = p;
  Q->Qhead = p;
  Q->Qtail = p;

  return p;
}



template<typename T>
void lrs_free_dic (lrs_dic<T> * P, lrs_dat<T> *Q)
{
/* do the same steps as for allocation, but backwards */
/* gmp variables cannot be cleared using free: use lrs_clear_mp* */
  lrs_dic<T> *P1;
  int i;

/* repeat until cache is empty */
  do
  {
    /* I moved these here because I'm not certain the cached dictionaries
       need to be the same size. Well, it doesn't cost anything to be safe. db */

  long m_A = P->m_A;

  for (i=0; i<=m_A; i++)
    delete [] P->A[i];
  delete [] P->A;

  delete [] P->Row;
  delete [] P->Col;
  delete [] P->C;
  delete [] P->B;

/* go to next record in cache if any */
  P1 =P->next;
  delete P;
  P=P1;

  }  while (Q->Qhead != P );


}

template<typename T>
void lrs_free_dat ( lrs_dat<T> *Q )
{
/* most of these items were allocated in lrs_alloc_dic */
  delete [] Q->inequality;
  delete [] Q->linearity;
  delete [] Q->facet;
  delete [] Q->redundcol;
  delete [] Q->minratio;
  delete [] Q->temparray;
  delete [] Q->saved_C;
  delete Q;
}


template<typename T>
long check_cache (lrs_dic<T> ** D_p, lrs_dat<T> * global, long *i_p, long *j_p)
{
  if (global->Qtail == global->Qhead)
    return 0;
  else
    {
      global->Qtail = global->Qtail->prev;

      *D_p = global->Qtail;

      *i_p = global->Qtail->i;
      *j_p = global->Qtail->j;

      return 1;
    }
}


template<typename T>
lrs_dic<T> * lrs_alloc_dic (lrs_dat<T> * Q)
{

  lrs_dic<T> *p;
  long i, j;
  long m, d, m_A;
  //  std::cerr << "Q->hull=" << Q->hull << "\n";
  if (Q->hull)                       /* d=col dimension of A */
    Q->inputd = Q->n;                /* extra column for hull */
  else
    Q->inputd = Q->n - 1;

  m = Q->m;
  d = Q->inputd;
  //  std::cerr << "m=" << m << " d=" << d << "\n";
  m_A = m;

/* nonnegative flag set means that problem is d rows "bigger"     */
/* since nonnegative constraints are not kept explicitly          */

  if(Q->nonnegative)
    m = m+d;

  p = new_lrs_dic<T> (m, d, m_A);
  if (!p)
    return nullptr;

  p->next = p;
  p->prev = p;
  Q->Qhead = p;
  Q->Qtail = p;

/* Initializations */

  p->d = p->d_orig = d;
  p->m = m;
  p->m_A  = m_A;
  p->depth = 0L;
  p->lexflag = globals::TRUE;
  p->det=1;

/*m+d+1 is the number of variables, labelled 0,1,2,...,m+d  */
/*  initialize array to zero   */
  for (i = 0; i <= m_A; i++)
    for (j = 0; j <= d; j++)
      p->A[i][j]=0;

  Q->inequality = new long[m+1];
  for (i=0; i<=m; i++)
    Q->inequality[i]=0;
  if (Q->nlinearity == 0)   /* linearity may already be allocated */
    Q->linearity = new long[m+1];
  for (i=0; i<=m; i++)
    Q->linearity[i]=0;

  Q->facet = new long[d+1];
  Q->redundcol = new long[d+1];
  Q->minratio = new long[m+1];
  Q->temparray = new long[d+1];

  Q->inequality[0] = 2L;
  Q->saved_C = new long[d+1];

  Q->lastdv = d;      /* last decision variable may be decreased */
                      /* if there are redundant columns          */

/*initialize basis and co-basis indices, and row col locations */
/*if nonnegative, we label differently to avoid initial pivots */
/* set basic indices and rows */
 if(Q->nonnegative)
  for (i = 0; i <= m; i++)
    {
      p->B[i] = i;
      if (i <= d )
	p->Row[i]=0; /* no row for decision variables */
      else 
	p->Row[i]=i-d;
    }
 else
   for (i = 0; i <= m; i++)
    {
      if (i == 0 )
          p->B[0]=0;
      else
          p->B[i] = d + i;
      p->Row[i] = i;
    }

  for (j = 0; j < d; j++)
    {
      if(Q->nonnegative)
          p->C[j] = m+j+1;
      else
          p->C[j] = j + 1;
      p->Col[j] = j + 1;
    }
  p->C[d] = m + d + 1;
  p->Col[d] = 0;
  return p;
}

/* 
   this routine makes a copy of the information needed to restart, 
   so that we can guarantee that if a signal is received, we 
   can guarantee that nobody is messing with it.
   This as opposed to adding all kinds of critical regions in 
   the main line code.

   It is also used to make sure that in case of overflow, we
   have a valid cobasis to restart from.
 */
template<typename T>
void save_basis (lrs_dic<T> * P, lrs_dat<T> * Q)
{
  int i;
/* assign local variables to structures */
  long *C = P->C;
  long d;
  d = P->d;
  Q->saved_flag = 1;
  for (i = 0; i < 3; i++)
    Q->saved_count[i] = Q->count[i];
  for (i = 0; i < d + 1; i++)
    Q->saved_C[i] = C[i];
  Q->saved_d = P->d;
  Q->saved_depth = P->depth;
}



template<typename T>
void lrs_set_row_mp(lrs_dic<T> *P, lrs_dat<T> *Q, long row, T* num, long ineq)
/* set row of dictionary using num and den arrays for rational input */
/* ineq = 1 (globals::GE)   - ordinary row  */
/*      = 0 (globals::EQ)   - linearity     */
{
  T Temp, mpone;

  long i, j;

/* assign local variables to structures */

  T** A;
  long hull;
  long d;
  hull = Q->hull;
  A = P->A;
  d = P->d;
  i=row;
  for (j = hull; j <= d; j++)       /* hull data copied to cols 1..d */
    A[i][j] = num[j-hull];
  if (hull)
    A[i][0]=0;
  if (!zero (A[i][hull]))   /* for H-rep, are zero in column 0     */
    Q->homogeneous = globals::FALSE; /* for V-rep, all zero in column 1     */
  if ( ineq == globals::EQ )        /* input is linearity */
    {
      Q->linearity[Q->nlinearity]=row;
      Q->nlinearity++;
    }
}

template<typename T>
void lrs_set_obj_mp(lrs_dic<T> *P, lrs_dat<T> *Q, T* num, long max)
{
  long i;

  if (max == globals::MAXIMIZE)
    Q->maximize=globals::TRUE;
  else
    {
      Q->minimize=globals::TRUE;
      for(i=0;i<=P->d;i++)
	num[i]=-num[i];
    }
  lrs_set_row_mp(P,Q,0L,num,globals::GE);
}

template<typename T>
long lrs_solve_lp(lrs_dic<T> *P, lrs_dat<T> *Q)
/* user callable function to solve lp only */
{
  T** Lin;		/* holds input linearities if any are found             */

  Q->lponly = globals::TRUE;

  if (!lrs_getfirstbasis (&P, Q, Lin, globals::FALSE))
    return globals::FALSE;
  return globals::TRUE;
} /* end of lrs_solve_lp */

template<typename T>
long dan_selectpivot (lrs_dic<T> * P, lrs_dat<T> * Q, long *r, long *s)
/* select pivot indices using dantzig simplex method             */
/* largest coefficient with lexicographic rule to avoid cycling  */
/* Bohdan Kaluzny's handiwork                                    */
/* returns globals::TRUE if pivot found else globals::FALSE                        */
/* pivot variables are B[*r] C[*s] in locations Row[*r] Col[*s]  */
{
  long j,k,col;
  T coeff;
/* assign local variables to structures */
  T **A = P->A;
  long *Col = P->Col;
  long d = P->d;

  *r = 0;
  *s = d;
  j = 0;
  k = 0;
 
  coeff=0;
/*find positive cost coef */
  while (k < d) 
     {
       if (lrs_greater(A[0][Col[k]], coeff))
        {
          j = k;
          coeff = A[0][Col[j]];
	}
      k++;
     }

  if (positive(coeff))
    {
      *s = j;
      col = Col[j];

      /*find min index ratio */
      *r = lrs_ratio<T> (P, Q, col);
      if (*r != 0)
        {
	return (globals::TRUE);		/* unbounded */
        }
    }
  return (globals::FALSE);
}

template<typename T>
long phaseone (lrs_dic<T> * P, lrs_dat<T> * Q)
/* Do a dual pivot to get primal feasibility (pivot in X_0)*/
/* Bohdan Kaluzny's handiwork                                    */
{
  long i, j, k;
/* assign local variables to structures */
  T *A = P->A;
  long *Row = P->Row;
  long *Col = P->Col;
  long m, d;
  T b_vector;
  m = P->m;
  d = P->d;
  i = 0;
  k = d+1;

  b_vector=0;
  while (k <= m)
     {
       if(lrs_greater(b_vector, A[Row[k]][0]))
        {
          i = k;
          b_vector=A[Row[i]][0];
        }
      k++;
     }

  if (negative(b_vector))                       /* pivot row found! */
    {
      j = 0;            /*find a positive entry for in row */
      while (j < d && !positive(A[Row[i]][Col[j]]))
        j++;
      if (j >= d)
        {
          return (globals::FALSE);       /* no positive entry */
        }
      pivot (P, Q, i, j);
      update (P, Q, &i, &j);
    }
  return (globals::TRUE);
}

/*
long 
lrs_read_dic (lrs_dic * P, lrs_dat * Q)
{
  lrs_mp Temp,Tempn,Tempd, mpone, mpten;
  lrs_mp_vector oD;

  long i, j;
  char name[100];
  int c;

  lrs_mp_matrix A;
  lrs_mp_vector Gcd, Lcm;
  long hull = Q->hull;
  long m, d;
  long dualperturb=FALSE;

  lrs_alloc_mp(Temp); lrs_alloc_mp(mpone);
  lrs_alloc_mp(Tempn); lrs_alloc_mp(Tempd); lrs_alloc_mp(mpten);
  A = P->A;
  m = Q->m;
  d = Q->inputd;
  Gcd = Q->Gcd;
  Lcm = Q->Lcm;

  oD = lrs_alloc_mp_vector (d);


  itomp (ONE, mpone);
  itomp (ONE, A[0][0]);
  itomp (ONE, Lcm[0]);
  itomp (ONE, Gcd[0]);

  for (i = 1; i <= m; i++)
    {
      itomp (ONE, Lcm[i]);
      itomp (ZERO, Gcd[i]);
      for (j = hull; j <= d; j++)
	{
	  if (readrat (A[i][j], A[0][j]))
	    lcm (Lcm[i], A[0][j]);
	  copy (Temp, A[i][j]);
	  gcd (Gcd[i], Temp);
	}
      if (!zero (A[i][hull]))
	Q->homogeneous = FALSE;
      storesign (Gcd[i], POS);
      storesign (Lcm[i], POS);
      if (greater (Gcd[i], mpone) || greater (Lcm[i], mpone))
	for (j = 0; j <= d; j++)
	  {
	    exactdivint (A[i][j], Gcd[i], Temp);
	    mulint (Lcm[i], Temp, Temp);
	    exactdivint (Temp, A[0][j], A[i][j]);
	  }

    }
  return TRUE;
}
*/

/* readrat is the function that reads off the values from the input file
   The entries in the first line are supposed to be the lcm or the Linear programming 
   The function lrs_read_dic remains potentially very dangerous 
   and we have to be very wary of its potential to create trouble 
   The readrat(Na, Da) returns the value Na/Da 
   Thus the read readrat (A[i][j], A[0][j])) puts the numerator in A[i][j] */
  
 

 
template<typename T>
void fillModelLRS(MyMatrix<T> const& EXT, lrs_dic<T> *P, lrs_dat<T> *Q)
{
  int j;
  int iRow, nbRow, nbCol;
  long n;
  long ineq;
  T* num;
  nbRow=EXT.rows();
  nbCol=EXT.cols();
  n=nbCol;
  num = new T[n+1];
  ineq=1;
  for (iRow=0; iRow<nbRow; iRow++) {
    for (j=0; j<nbCol; ++j)
      num[j]=EXT(iRow, j);
    lrs_set_row_mp(P, Q, iRow+1, num, ineq);
  }
  for (j=0; j<nbCol; j++)
    P->A[0][j]=1;
  delete [] num;
}


template<typename T>
void PrintP(lrs_dic<T>* & P, std::string const& message)
{
  std::cerr << "message = " << message << "\n";
  std::cerr << "P->A=\n";
  for (int i=0; i<=P->m_A; i++) {
    for (int j=0; j<=P->d; j++)
      std::cerr << " " << P->A[i][j];
    std::cerr << "\n";
  }

  std::cerr << "B=";
  for (int i=0; i<=P->m; i++)
    std::cerr << P->B[i] << " ";
  std::cerr << "\n";
  std::cerr << "Row=";
  for (int i=0; i<=P->m; i++)
    std::cerr << P->Row[i] << " ";
  std::cerr << "\n";
  //
  std::cerr << "C=";
  for (int i=0; i<=P->d; i++)
    std::cerr << P->C[i] << " ";
  std::cerr << "\n";
  std::cerr << "Col=";
  for (int i=0; i<=P->d; i++)
    std::cerr << P->Col[i] << " ";
  std::cerr << "\n";
}

 
template<typename T>
void initLRS(MyMatrix<T> const& EXT, lrs_dic<T>* & P, lrs_dat<T>* & Q)
{
  T** Lin;
  Q = lrs_alloc_dat<T> ();
  //  std::cerr << "After lrs_alloc_dat\n";
  if (Q == nullptr) {
    throw TerminalException{1};
  }
  int nbrow=EXT.rows();
  int nbcol=EXT.cols();
  Q->n = nbcol;
  Q->m = nbrow;
  //  std::cerr << "no need to call a lrs_read_dat a priori\n";
  P = lrs_alloc_dic (Q);
  //  PrintP(P, "after lrs_alloc_dic");
  if (P == nullptr) {
    std::cerr << "We failed allocation, let's die\n";
    throw TerminalException{1};
  }
  fillModelLRS(EXT, P, Q);
  //  PrintP(P, "Before lrs_getfirstbasis");

  if (!lrs_getfirstbasis (&P, Q, Lin, globals::FALSE)) {
    std::cerr << "Error in call to lrs_getfirstbasis\n";
    throw TerminalException{1};
  }
  //  std::cerr << "After lrs_getfirstbasis\n";
}

template<typename T>
void freeLRS(lrs_dic<T>* & P, lrs_dat<T>* & Q)
{
  lrs_free_dic(P,Q);
  lrs_free_dat(Q);
}

template<typename T>
void Kernel_DualDescription(MyMatrix<T> const& EXT, std::function<void(T*)> const& f)
{
  lrs_dic<T> *P;
  lrs_dat<T> *Q;
  T* output;
  int col;
  //  std::cerr << "Before call to initLRS in DualDescription_temp_incd\n";
  initLRS(EXT, P, Q);
  output=new T[Q->n+1];
  unsigned long dict_count = 1;
  /*
  int nbCol=EXT.cols();
  int nbRow=EXT.rows();
  int nbIter=0;
  std::cerr << "nbRow=" << nbRow << "\n";
  std::cerr << "P->d=" << P->d << "\n";
  std::cerr << "Q->hull=" << Q->hull << "  Q->n=" << Q->n << "\n";*/
  do {
    /*
      nbIter++;
      std::cerr << "nbIter=" << nbIter << " nbCol=" << nbCol << "\n";
      PrintP(P, "Before lrs_getsolution loop");*/
    for (col = 0; col <= P->d; col++) {
      if (lrs_getsolution (P, Q, output, col)) {
	f(output);
      }
    }
  }
  while (lrs_getnextbasis (&P, Q, globals::FALSE, dict_count));
  delete [] output;
  lrs_free_dic (P,Q);
  lrs_free_dat (Q);
}



template<typename T>
void Kernel_DualDescription_limited(MyMatrix<T> const& EXT, std::function<void(T*)> const& f, int const& UpperLimit)
{
  lrs_dic<T> *P;
  lrs_dat<T> *Q;
  T* output;
  int col;
  initLRS(EXT, P, Q);
  output=new T[Q->n+1];
  unsigned long dict_count = 1;
  int nbDone=0;
  do {
    for (col = 0; col <= P->d; col++)
      if (lrs_getsolution (P, Q, output, col)) {
	f(output);
	nbDone++;
      }
    if (nbDone > UpperLimit)
      break;
  }
  while (lrs_getnextbasis (&P, Q, globals::FALSE, dict_count));
  delete [] output;
  lrs_free_dic (P,Q);
  lrs_free_dat (Q);
}

 

template<typename T>
MyMatrix<T> FirstColumnZero(MyMatrix<T> const& M)
{
  int nbRow=M.rows();
  int nbCol=M.cols();
  for (int iRow=0; iRow<nbRow; iRow++) {
    T eVal=M(iRow,0);
    if (eVal != 0) {
      MyMatrix<T> Mret(nbRow,nbCol+1);
      for (int jRow=0; jRow<nbRow; jRow++) {
	Mret(jRow,0)=0;
	for (int iCol=0; iCol<nbCol; iCol++)
	  Mret(jRow,iCol+1)=M(jRow,iCol);
      }
      return Mret;
    }
  }
  return M;
}

 
template<typename T>
std::vector<Face> DualDescription_temp_incd(MyMatrix<T> const& EXT)
{
  MyMatrix<T> EXTwork=FirstColumnZero(EXT);
  int nbCol=EXTwork.cols();
  int nbRow=EXTwork.rows();
  std::vector<Face> ListIncd;
  bool IsFirst=true;
  std::function<void(T*)> f=[&](T* out) -> void {
    if (!IsFirst) {
      Face V(nbRow);
      for (int iRow=0; iRow<nbRow; iRow++) {
	T eScal=0;
	for (int iCol=0; iCol<nbCol; iCol++)
	  eScal += out[iCol]*EXTwork(iRow,iCol);
	if (eScal == 0)
	  V[iRow]=1;
      }
      ListIncd.push_back(V);
    }
    IsFirst=false;
  };
  Kernel_DualDescription(EXTwork, f);
  return ListIncd;
}



template<typename T>
std::vector<Face> DualDescription_temp_incd_limited(MyMatrix<T> const& EXT, int const& UpperLimit)
{
  MyMatrix<T> EXTwork=FirstColumnZero(EXT);
  int nbCol=EXTwork.cols();
  int nbRow=EXTwork.rows();
  std::vector<Face> ListIncd;
  bool IsFirst=true;
  std::function<void(T*)> f=[&](T* out) -> void {
    if (!IsFirst) {
      Face V(nbRow);
      for (int iRow=0; iRow<nbRow; iRow++) {
	T eScal=0;
	for (int iCol=0; iCol<nbCol; iCol++)
	  eScal += out[iCol]*EXTwork(iRow,iCol);
	if (eScal == 0)
	  V[iRow]=1;
      }
      ListIncd.push_back(V);
    }
    IsFirst=false;
  };
  Kernel_DualDescription_limited(EXTwork, f, UpperLimit);
  return ListIncd;
}




template<typename T>
std::vector<Face> DualDescription_temp_incd_reduction(MyMatrix<T> const& EXT)
{
  MyMatrix<T> EXTwork=FirstColumnZero(EXT);
  using Tring = typename underlying_ring<T>::ring_type;
  //  typedef typename underlying_ring<T>::ring_type Tring;
  int nbCol=EXTwork.cols();
  int nbRow=EXTwork.rows();
  MyMatrix<Tring> EXTring(nbRow,nbCol);
  for (int iRow=0; iRow<nbRow; iRow++) {
    MyVector<T> eRow1=GetMatrixRow(EXTwork, iRow);
    MyVector<T> eRow2=RemoveFractionVector(eRow1);
    MyVector<Tring> eRow3=ConvertVectorUniversal<Tring,T>(eRow2);
    AssignMatrixRow(EXTring, iRow, eRow3);
  }
  std::vector<Face> ListIncd;
  bool IsFirst=true;
  std::function<void(Tring*)> f=[&](Tring* out) -> void {
    if (!IsFirst) {
      Face V(nbRow);
      for (int iRow=0; iRow<nbRow; iRow++) {
	T eScal=0;
	for (int iCol=0; iCol<nbCol; iCol++)
	  eScal += out[iCol]*EXTring(iRow,iCol);
	if (eScal == 0)
	  V[iRow]=1;
      }
      ListIncd.push_back(V);
    }
    IsFirst=false;
  };
  Kernel_DualDescription(EXTring, f);
  return ListIncd;
}



 
}
#endif

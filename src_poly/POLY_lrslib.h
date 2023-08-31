#ifndef SRC_POLY_POLY_LRSLIB_H_
#define SRC_POLY_POLY_LRSLIB_H_
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

// clang-format off
#include "Boost_bitset.h"
#include "COMB_Stor.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include <limits>
#include <string>
// clang-format on

#ifdef PRINT
# define PRINT_LRS_ANALYSIS
#endif

namespace lrs {
// some #defines and global variables from the original lrs code
namespace globals {
const int64_t POS = 1L;
const int64_t NEG = -1L;
const int64_t L_FALSE = 0L;
const int64_t L_TRUE = 1L;
const int64_t MAXIMIZE = 1L; /* maximize the lp  */
const int64_t MINIMIZE = 0L; /* maximize the lp  */
const int64_t GE = 1L;       /* constraint is >= */
const int64_t EQ = 0L;       /* constraint is linearity */
const uint64_t dict_limit = 10;
} // namespace globals

template <typename T> struct lrs_dic {
  T **A;
  int64_t m;        /* A has m+1 rows, row 0 is cost row            */
  int64_t m_A;      /* =m or m-d if nonnegative flag set            */
  int64_t d;        /* A has d+1 columns, col 0 is b-vector         */
  int64_t d_orig;   /* value of d as A was allocated  (E.G.)        */
  int64_t lexflag;  /* true if lexmin basis for this vertex         */
  int64_t depth;    /* depth of basis/vertex in reverse search tree */
  int64_t i, j;     /* last pivot row and column pivot indices      */
  T det;            /* current determinant of basis                 */
  int64_t *B, *Row; /* basis, row location indices                  */
  int64_t *C, *Col; /* cobasis, column location indices             */
  lrs_dic *prev, *next;
};

template <typename T> struct lrs_dat {
  int64_t unbounded; /* lp unbounded */

  int64_t *inequality; /* indices of inequalities corr. to cobasic ind */
  /* initially holds order used to find starting  */
  /* basis, default: m,m-1,...,2,1                */
  int64_t *facet;     /* cobasic indices for restart in needed        */
  int64_t *redundcol; /* holds columns which are redundant            */
  int64_t *linearity; /* holds cobasic indices of input linearities   */
  int64_t *minratio;  /* used for lexicographic ratio test            */
  int64_t inputd;     /* input dimension: n-1 for H-rep, n for V-rep  */

  int64_t m;      /* number of rows in input file                 */
  int64_t n;      /* number of columns in input file              */
  int64_t lastdv; /* index of last dec. variable after preproc    */
  /* given by inputd-nredundcol                   */
  int64_t count[10]; /* count[0]=rays [1]=verts. [2]=base [3]=pivots */
  /* count[4]=integer vertices                    */
  int64_t nredundcol; /* number of redundant columns                  */
  int64_t nlinearity; /* number of input linearities                  */
  int64_t runs;       /* probes for estimate function                 */
  int64_t seed;       /* seed for random number generator             */
  /**** flags  **********                         */
  int64_t bound;     /* globals::TRUE if upper/lower bound on objective given */
  int64_t dualdeg;   /* globals::TRUE if start dictionary is dual degenerate  */
  int64_t getvolume; /* do volume calculation                        */
  int64_t givenstart;  /* globals::TRUE if a starting cobasis is given  */
  int64_t homogeneous; /* globals::TRUE if all entries in column one are zero */
  int64_t hull;      /* do convex hull computation if globals::TRUE           */
  int64_t lponly;    /* true if only lp solution wanted              */
  int64_t maxdepth;  /* max depth to search to in treee              */
  int64_t maximize;  /* flag for LP maximization                     */
  int64_t minimize;  /* flag for LP minimization                     */
  int64_t mindepth;  /* do not backtrack above mindepth              */
  int64_t
      nonnegative;  /* globals::TRUE if last d constraints are nonnegativity */
  int64_t polytope; /* globals::TRUE for facet computation of a polytope     */
  int64_t truncate; /* globals::TRUE: truncate tree when moving from opt vert*/
  int64_t restart;  /* globals::TRUE if restarting from some cobasis         */

  /* Variables for saving/restoring cobasis,  db */

  int64_t id; /* numbered sequentially */

  /* Variables for cacheing dictionaries, db */
  lrs_dic<T> *Qhead, *Qtail;
};

template <typename T> inline void storesign(T &a, int64_t const &sa) {
  if (a > 0) {
    if (sa < 0)
      a = -a;
  } else {
    if (a < 0) {
      if (sa > 0)
        a = -a;
    }
  }
}

template <typename T> inline int64_t sign(T const &a) {
  if (a < 0)
    return globals::NEG;
  if (a > 0)
    return globals::POS;
  return 0;
}

template <typename T> inline int comprod(T const& Na, T const& Nb, T const& Nc, T const& Nd) {
  T prod1 = Na * Nb;
  T prod2 = Nc * Nd;
  if (prod1 > prod2)
    return 1;
  if (prod1 < prod2)
    return -1;
  return 0;
}

template <typename T>
int64_t lrs_getsolution(lrs_dic<T> *P, lrs_dat<T> *Q, T *&output, int64_t col)
/* check if column indexed by col in this dictionary */
/* contains output                                   */
/* col=0 for vertex 1....d for ray/facet             */
{
  int64_t j; /* cobasic index     */
  T **A = P->A;
  int64_t *Row = P->Row;
  if (col == 0) {
    return lrs_getvertex(P, Q, output);
  }

  if (Q->lponly) {
    if (A[0][col] <= 0) {
      return globals::L_FALSE;
    }
  } else if (A[0][col] >= 0) {
    return globals::L_FALSE;
  }

  j = Q->lastdv + 1;
  while (j <= P->m && A[Row[j]][col] >= 0)
    j++;

  if (j <= P->m) {
    return globals::L_FALSE;
  }

  if (lexmin(P, Q, col) || Q->lponly)
    return lrs_getray(P, Q, col, Q->n, output);
  return globals::L_FALSE; /* no more output in this dictionary */
}

/***********************************/
/* allocate and initialize lrs_dat */
/***********************************/
template <typename T> lrs_dat<T> *lrs_alloc_dat() {
  lrs_dat<T> *Q;
  Q = new lrs_dat<T>;

  /* initialize variables */
  Q->m = 0L;
  Q->n = 0L;
  Q->inputd = 0L;
  Q->nlinearity = 0L;
  Q->nredundcol = 0L;
  Q->runs = 0L;
  Q->seed = 1234L;
  Q->bound =
      globals::L_FALSE; /* upper/lower bound on objective function given */
  Q->homogeneous = globals::L_TRUE;
  Q->hull = globals::L_FALSE;
  Q->lponly = globals::L_FALSE;
  Q->maxdepth = std::numeric_limits<int64_t>::max();
  Q->mindepth = std::numeric_limits<int64_t>::min();
  Q->nonnegative = globals::L_FALSE;
  Q->truncate =
      globals::L_FALSE; /* truncate tree when moving from opt vertex        */
  Q->maximize = globals::L_FALSE; /*flag for LP maximization */
  Q->minimize = globals::L_FALSE; /*flag for LP minimization */
  Q->restart =
      globals::L_FALSE; /* globals::TRUE if restarting from some cobasis */
  Q->givenstart =
      globals::L_FALSE; /* globals::TRUE if a starting cobasis is given */
  return Q;
}

void reorder1(int64_t a[], int64_t b[], int64_t newone, int64_t range)
/*reorder array a in increasing order with one misplaced element at index newone
 */
/*elements of array b are updated to stay aligned with a */
{
  int64_t temp;
  while (newone > 0 && a[newone] < a[newone - 1]) {
    temp = a[newone];
    a[newone] = a[newone - 1];
    a[newone - 1] = temp;
    temp = b[newone];
    b[newone] = b[newone - 1];
    b[--newone] = temp;
  }
  while (newone < range - 1 && a[newone] > a[newone + 1]) {
    temp = a[newone];
    a[newone] = a[newone + 1];
    a[newone + 1] = temp;
    temp = b[newone];
    b[newone] = b[newone + 1];
    b[++newone] = temp;
  }
}

template <typename T>
int64_t lrs_getfirstbasis(lrs_dic<T> **D_p, lrs_dat<T> *Q, T **&Lin)
/* gets first basis, globals::FALSE if none              */
/* P may get changed if lin. space Lin found    */
/* no_output is globals::TRUE supresses output headers   */
{
  int64_t i, j, k;

  /* assign local variables to structures */

  T **A;
  int64_t *B, *C, *Col;
  int64_t *inequality;
  int64_t *linearity;
  int64_t hull = Q->hull;
  int64_t m, d, lastdv, nlinearity, nredundcol;
  m = (*D_p)->m;
  d = (*D_p)->d;
  lastdv = Q->lastdv;
  //  std::cerr << "lrs_getfirstbasis lastdv = " << lastdv << "\n";

  nredundcol = 0L;            /* will be set after getabasis        */
  nlinearity = Q->nlinearity; /* may be reset if new linearity read */
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
  for (i = 0; i < nlinearity; i++) /* put linearities first in the order */
    inequality[i] = linearity[i];

  k = 0; /* index for linearity array   */
  if (Q->givenstart)
    k = d;
  else
    k = nlinearity;
  for (i = m; i >= 1; i--) {
    j = 0;
    while (j < k && inequality[j] != i)
      j++; /* see if i is in inequality  */
    if (j == k)
      inequality[k++] = i;
  }
  //  std::cerr << "Q->maximize/minimize = " << Q->maximize << " / " <<
  //  Q->minimize << "\n";

  /* for voronoi convert to h-description using the transform */
  /* a_0 .. a_d-1 -> (a_0^2 + ... a_d-1 ^2)-2a_0x_0-...-2a_d-1x_d-1 + x_d >= 0
   */
  /* note constant term is stored in column d, and column d-1 is all ones */
  /* the other coefficients are multiplied by -2 and shifted one to the right */
  if (!Q->maximize && !Q->minimize)
    for (j = 0; j <= d; j++)
      A[0][j] = 0;

  /* Now we pivot to standard form, and then find a primal feasible basis */
  /* Note these steps MUST be done, even if restarting, in order to get */
  /* the same index/inequality correspondance we had for the original prob. */
  /* The inequality array is used to give the insertion order */
  /* and is defaulted to the last d rows when givenstart=globals::FALSE */

  if (Q->nonnegative) {
    /* no need for initial pivots here, labelling already done */
    Q->lastdv = d;
    Q->nredundcol = 0;
  } else {
    if (!getabasis(*D_p, Q, inequality)) {
      //      std::cerr << "Exit case 1\n";
      return globals::L_FALSE;
    }
  }
  nredundcol = Q->nredundcol;
  lastdv = Q->lastdv;
  d = (*D_p)->d;
  /* Reset up the inequality array to remember which index is which input
   * inequality */
  /* inequality[B[i]-lastdv] is row number of the inequality with index B[i] */
  /* inequality[C[i]-lastdv] is row number of the inequality with index C[i] */

  for (i = 1; i <= m; i++)
    inequality[i] = i;
  if (nlinearity > 0) {              /* some cobasic indices will be removed */
    for (i = 0; i < nlinearity; i++) /* remove input linearity indices */
      inequality[linearity[i]] = 0;
    k = 1; /* counter for linearities         */
    for (i = 1; i <= m - nlinearity; i++) {
      while (k <= m && inequality[k] == 0)
        k++; /* skip zeroes in corr. to linearity */
      inequality[i] = inequality[k++];
    }
  } /* end if linearity */
  if (nredundcol > 0) {
    Lin = new T *[nredundcol + 1];
    for (i = 0; i < nredundcol; i++) {
      Lin[i] = new T[Q->n + 1];
      if (!(Q->homogeneous && Q->hull &&
            i == 0)) /* skip redund col 1 for homog. hull */
        lrs_getray(*D_p, Q, Col[0], (*D_p)->C[0] + i - hull,
                   Lin[i]); /* adjust index for deletions */
      if (!removecobasicindex(*D_p, 0)) {
        //	    std::cerr << "Exit case 2\n";
        return globals::L_FALSE;
      }
    }
  }
  if (!primalfeasible(*D_p, Q)) {
    //    std::cerr << "Exit case 3\n";
    return globals::L_FALSE;
  }

  /* Now solve LP if objective function was given */
  if (Q->maximize || Q->minimize) {
    Q->unbounded = !lrs_solvelp(*D_p, Q);
    if (Q->lponly) {
      //	std::cerr << "Exit case 4\n";
      return globals::L_TRUE;
    } else { /* check to see if objective is dual degenerate */
      j = 1;
      while (j <= d && A[0][j] != 0)
        j++;
    }
  } else {
    for (j = 1; j <= d; j++) {
      A[0][j] = (*D_p)->det;
      storesign(A[0][j], globals::NEG);
    }
    A[0][0] = 0;
  }
  while (C[0] <= m) {
    i = C[0];
    j = inequality[B[i] - lastdv];
    inequality[B[i] - lastdv] = inequality[C[0] - lastdv];
    inequality[C[0] - lastdv] = j;
    C[0] = B[i];
    B[i] = i;
    reorder1(C, Col, 0, d);
  }
  /* Check to see if necessary to resize */
  if (Q->inputd > (*D_p)->d)
    *D_p = resize(*D_p, Q);

  //  std::cerr << "Exit case 6\n";
  return globals::L_TRUE;
}

/*****************************************/
/* getnextbasis in reverse search order  */
/*****************************************/

template <typename T>
int64_t lrs_getnextbasis(lrs_dic<T> **D_p, lrs_dat<T> *Q, int64_t backtrack,
                         uint64_t &dict_count)
/* gets next reverse search tree basis, globals::FALSE if none  */
/* switches to estimator if maxdepth set               */
/* backtrack globals::TRUE means backtrack from here            */
{
  /* assign local variables to structures */
  int64_t i = 0L, j = 0L;
  int64_t m = (*D_p)->m;
  int64_t d = (*D_p)->d;

  if (backtrack && (*D_p)->depth == 0) {
    return globals::L_FALSE; /* cannot backtrack from root      */
  }

  //  PrintP(*D_p, "Before while loop");
  while ((j < d) || ((*D_p)->B[m] != m)) {
    if ((*D_p)->depth >= Q->maxdepth) {
      backtrack = globals::L_TRUE;
      if (Q->maxdepth == 0) { /* estimate only */
        return globals::L_FALSE;
      }
    }
    if (Q->truncate &&
        (*D_p)->A[0][0] < 0) /* truncate when moving from opt. vertex */
      backtrack = globals::L_TRUE;

    //      PrintP(*D_p, "Before backtrack test");
    if (backtrack) { /* go back to prev. dictionary, restore i,j */
      backtrack = globals::L_FALSE;

      if (!check_cache(D_p, Q, &i, &j)) {
        (*D_p)->depth--;
        selectpivot(*D_p, Q, &i, &j);
        pivot(*D_p, i, j);
        update(*D_p, &i, &j);
      }

      j++; /* go to next column */
    }
    //      PrintP(*D_p, "After backtrack test");

    if ((*D_p)->depth < Q->mindepth)
      break;
    while ((j < d) && !reverse(*D_p, Q, &i, j))
      j++;

    //      PrintP(*D_p, "Before exiting test");
    if (j == d)
      backtrack = globals::L_TRUE;
    else {
      cache_dict(D_p, Q, i, j, dict_count);
      //	  PrintP(*D_p, "After cache_dict");
      /* Note that the next two lines must come _after_ the
         call to cache_dict */

      (*D_p)->depth++;

      pivot(*D_p, i, j);
      //	  PrintP(*D_p, "After pivot");
      update(*D_p, &i, &j);
      //	  PrintP(*D_p, "After update");

      (*D_p)->lexflag = lexmin(*D_p, Q, 0); /* see if lexmin basis */
      return globals::L_TRUE;
    }
  } /* end of main while loop for getnextbasis */
  return globals::L_FALSE;
}

/*************************************/
/* print out one line of output file */
/*************************************/
template <typename T>
int64_t lrs_getvertex(lrs_dic<T> *P, lrs_dat<T> *Q, T *&output)
/*Print out current vertex if it is lexmin and return it in output */
/* return globals::FALSE if no output generated  */
{

  int64_t i;
  int64_t ind;  /* output index                                  */
  int64_t ired; /* counts number of redundant columns            */
                /* assign local variables to structures */
  int64_t *redundcol = Q->redundcol;

  int64_t hull;
  int64_t lexflag;

  hull = Q->hull;
  lexflag = P->lexflag;
  if (hull)
    return globals::L_FALSE; /* skip printing the origin */

  if (!lexflag && !Q->lponly) /* not lexmin, and not printing forced */
    return globals::L_FALSE;

  /* copy column 0 to output */

  i = 1;
  ired = 0;
  output[0] = P->det;

  for (ind = 1; ind < Q->n; ind++) { /* extract solution */
    if ((ired < Q->nredundcol) &&
        (redundcol[ired] == ind)) { /* column was deleted as redundant */
      output[ind] = 0;
      ired++;
    } else { /* column not deleted as redundant */
      getnextoutput(P, Q, i, 0, output[ind]);
      i++;
    }
  }
  return globals::L_TRUE;
}

template <typename T>
int64_t lrs_getray(lrs_dic<T> *P, lrs_dat<T> *Q, int64_t col, int64_t redcol,
                   T *&output)
/*Print out solution in col and return it in output   */
/*redcol =n for ray/facet 0..n-1 for linearity column */
/*hull=1 implies facets will be recovered             */
/* return globals::FALSE if no output generated in column col  */
{
  int64_t i;
  int64_t ind;  /* output index                                  */
  int64_t ired; /* counts number of redundant columns            */
                /* assign local variables to structures */
  int64_t *redundcol = Q->redundcol;
  int64_t hull = Q->hull;
  int64_t n = Q->n;
  i = 1;
  ired = 0;

  for (ind = 0; ind < n; ind++) { /* print solution */
    if (ind == 0 && !hull) /* must have a ray, set first column to zero */
      output[0] = 0;
    else if ((ired < Q->nredundcol) &&
             (redundcol[ired] == ind)) { /* column was deleted as redundant */
      if (redcol == ind) /* true for linearity on this cobasic index */
        /* we print reduced determinant instead of zero */
        output[ind] = P->det;
      else
        output[ind] = 0;
      ired++;
    } else { /* column not deleted as redundant */
      getnextoutput(P, Q, i, col, output[ind]);
      i++;
    }
  }
  return globals::L_TRUE;
}

template <typename T>
void getnextoutput(lrs_dic<T> *P, lrs_dat<T> *Q, int64_t i, int64_t col,
                   T &out) {
  int64_t row;
  int64_t m = P->m;
  int64_t d = P->d;
  int64_t lastdv = Q->lastdv;
  T **A = P->A;
  int64_t *B = P->B;
  int64_t *Row = P->Row;
  int64_t j;
  row = Row[i];
  if (Q->nonnegative) {
    for (j = lastdv + 1; j <= m; j++) {
      if (Q->inequality[B[j] - lastdv] == m - d + i) {
        out = A[Row[j]][col];
        return;
      }
    }
    /* did not find inequality m-d+i in basis */
    if (i == col)
      out = P->det;
    else
      out = 0;
  } else {
    out = A[row][col];
  }
}

template <typename T>
int64_t lrs_ratio(lrs_dic<T> *P, lrs_dat<T> *Q,
                  int64_t col) /*find lex min. ratio */
/* find min index ratio -aig/ais, ais<0 */
/* if multiple, checks successive basis columns */
/* recoded Dec 1997                     */
{
  int64_t i, j, comp, ratiocol, basicindex, start, nstart, cindex, bindex;
  int64_t firstime; /*For ratio test, true on first pass,else false */
  T Nmin, Dmin;
  int64_t degencount, ndegencount;
  /* assign local variables to structures */
  T **A = P->A;
  int64_t *B = P->B;
  int64_t *Row = P->Row;
  int64_t *Col = P->Col;
  int64_t *minratio = Q->minratio;
  int64_t m, d, lastdv;

  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

  nstart = 0;
  ndegencount = 0;
  degencount = 0;
  for (j = lastdv + 1; j <= m; j++) {
    /* search rows with negative coefficient in dictionary */
    /*  minratio contains indices of min ratio cols        */
    if (A[Row[j]][col] < 0)
      minratio[degencount++] = j;
  } /* end of for loop */
  if (degencount == 0)
    return (degencount); /* non-negative pivot column */

  ratiocol = 0;   /* column being checked, initially rhs */
  start = 0;      /* starting location in minratio array */
  bindex = d + 1; /* index of next basic variable to consider */
  cindex = 0;     /* index of next cobasic variable to consider */
  basicindex =
      d; /* index of basis inverse for current ratio test, except d=rhs test */
  while (degencount > 1) {         /*keep going until unique min ratio found */
    if (B[bindex] == basicindex) { /* identity col in basis inverse */
      if (minratio[start] == bindex) { /* remove this index, all others stay */
        start++;
        degencount--;
      }
      bindex++;
    } else { /* perform ratio test on rhs or column of basis inverse */
      firstime = globals::L_TRUE;
      /*get next ratio column and increment cindex */
      if (basicindex != d)
        ratiocol = Col[cindex++];
      for (j = start; j < start + degencount; j++) {
        i = Row[minratio[j]]; /* i is the row location of the next basic
                                 variable */
        comp = 1;             /* 1:  lhs>rhs;  0:lhs=rhs; -1: lhs<rhs */
        if (firstime)
          firstime = globals::L_FALSE; /*force new min ratio on first time */
        else {
          if (Nmin > 0 || A[i][ratiocol] < 0) {
            if (Nmin < 0 || A[i][ratiocol] > 0)
              comp = comprod(Nmin, A[i][col], A[i][ratiocol], Dmin);
            else
              comp = -1;
          } else if (Nmin == 0 && A[i][ratiocol] == 0)
            comp = 0;
          if (ratiocol == 0)
            comp = -comp; /* all signs reversed for rhs */
        }
        if (comp == 1) { /*new minimum ratio */
          nstart = j;
          Nmin = A[i][ratiocol];
          Dmin = A[i][col];
          ndegencount = 1;
        } else if (comp == 0) /* repeated minimum */
          minratio[nstart + ndegencount++] = minratio[j];
      } /* end of  for (j=start.... */
      degencount = ndegencount;
      start = nstart;
    }             /* end of else perform ratio test statement */
    basicindex++; /* increment column of basis inverse to check next */
  }               /*end of while loop */
  return minratio[start];
} /* end of ratio */

template <typename T>
int64_t reverse(lrs_dic<T> *P, lrs_dat<T> *Q, int64_t *r, int64_t s)
/*  find reverse indices  */
/* globals::TRUE if B[*r] C[s] is a reverse lexicographic pivot */
{
  int64_t i, j, row, col;

  /* assign local variables to structures */
  T **A = P->A;
  int64_t *B = P->B;
  int64_t *C = P->C;
  int64_t *Row = P->Row;
  int64_t *Col = P->Col;
  int64_t d = P->d;

  col = Col[s];
  if (A[0][col] >= 0)
    return globals::L_FALSE;

  *r = lrs_ratio<T>(P, Q, col);
  if (*r == 0) /* we have a ray */
    return globals::L_FALSE;

  row = Row[*r];

  /* check cost row after "pivot" for smaller leaving index    */
  /* ie. j s.t.  A[0][j]*A[row][col] < A[0][col]*A[row][j]     */
  /* note both A[row][col] and A[0][col] are negative          */

  for (i = 0; i < d && C[i] < B[*r]; i++) {
    if (i != s) {
      j = Col[i];
      if (A[0][j] > 0 || A[row][j] < 0) /*or else sign test fails trivially */
        if ((A[0][j] >= 0 && A[row][j] <= 0) ||
            comprod(A[0][j], A[row][col], A[0][col], A[row][j]) == -1)
          return globals::L_FALSE;
    }
  }
  return globals::L_TRUE;
}

template <typename T>
int64_t selectpivot(lrs_dic<T> *P, lrs_dat<T> *Q, int64_t *r, int64_t *s)
/* select pivot indices using lexicographic rule   */
/* returns globals::TRUE if pivot found else globals::FALSE          */
/* pivot variables are B[*r] C[*s] in locations Row[*r] Col[*s] */
{
  int64_t j, col;
  /* assign local variables to structures */
  T **A = P->A;
  int64_t *Col = P->Col;
  int64_t d = P->d;

  *r = 0;
  *s = d;
  j = 0;
  while (j < d && A[0][Col[j]] <= 0)
    j++;

  if (j < d) {
    *s = j;
    col = Col[j];
    *r = lrs_ratio<T>(P, Q, col);
    if (*r != 0)
      return globals::L_TRUE;
  }
  return globals::L_FALSE;
}

template <typename T>
void pivot(lrs_dic<T> *P, int64_t bas, int64_t cob)
/* Qpivot routine for array A              */
/* indices bas, cob are for Basis B and CoBasis C    */
/* corresponding to row Row[bas] and column       */
/* Col[cob]   respectively                       */
{
  int64_t r, s;
  int64_t i, j;
  T Ars;
  /* assign local variables to structures */

  T **A = P->A;
  int64_t *Row = P->Row;
  int64_t *Col = P->Col;
  int64_t d, m_A;

  d = P->d;
  m_A = P->m_A;

  r = Row[bas];
  s = Col[cob];

  /* Ars=A[r][s]    */
  Ars = A[r][s];
  storesign(P->det, sign(Ars)); /*adjust determinant to new sign */

  for (i = 0; i <= m_A; i++)
    if (i != r)
      for (j = 0; j <= d; j++)
        if (j != s) {
          /*        A[i][j]=(A[i][j]*Ars-A[i][s]*A[r][j])/P->det; */
          A[i][j] = (A[i][j] * Ars - A[i][s] * A[r][j]) / P->det;
        }

  if (Ars > 0) {
    for (j = 0; j <= d; j++) /* no need to change sign if Ars neg */
      /*   A[r][j]=-A[r][j];              */
      if (A[r][j] != 0)
        A[r][j] = -A[r][j];
  } else {
    for (i = 0; i <= m_A; i++)
      if (A[i][s] != 0)
        A[i][s] = -A[i][s];
  }

  A[r][s] = P->det;
  P->det = Ars;
  storesign(P->det, globals::POS);
}

template <typename T>
int64_t primalfeasible(lrs_dic<T> *P, lrs_dat<T> *Q)
/* Do dual pivots to get primal feasibility */
/* Note that cost row is all zero, so no ratio test needed for Dual Bland's rule
 */
{
  int64_t primalinfeasible = globals::L_TRUE;
  int64_t i, j;
  /* assign local variables to structures */
  T **A = P->A;
  int64_t *Row = P->Row;
  int64_t *Col = P->Col;
  int64_t m, d, lastdv;
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;
  while (primalinfeasible) {
    i = lastdv + 1;
    while (i <= m && A[Row[i]][0] >= 0)
      i++;
    if (i <= m) {
      j = 0; /*find a positive entry for in row */
      while (j < d && A[Row[i]][Col[j]] <= 0) {
        j++;
      }
      if (j >= d)
        return globals::L_FALSE; /* no positive entry */
      pivot(P, i, j);
      update(P, &i, &j);
    } else
      primalinfeasible = globals::L_FALSE;
  } /* end of while primalinfeasibile */
  return globals::L_TRUE;
} /* end of primalfeasible */

template <typename T>
int64_t lrs_solvelp(lrs_dic<T> *P, lrs_dat<T> *Q)
/* Solve primal feasible lp by Dantzig`s rule and lexicographic ratio test */
/* return globals::TRUE if bounded, globals::FALSE if unbounded */
{
  int64_t i, j;
  /* assign local variables to structures */
  int64_t d = P->d;

  while (dan_selectpivot(P, Q, &i, &j)) {
    pivot(P, i, j);
    update(P, &i, &j);
  }

  if (j < d &&
      i == 0) { /* selectpivot gives information on unbounded solution */
    return globals::L_FALSE;
  }
  return globals::L_TRUE;
}

template <typename T>
int64_t getabasis(lrs_dic<T> *P, lrs_dat<T> *Q, int64_t order[])
/* Pivot Ax<=b to standard form */
/*Try to find a starting basis by pivoting in the variables x[1]..x[d]        */
/*If there are any input linearities, these appear first in order[]           */
/* Steps: (a) Try to pivot out basic variables using order                    */
/*            Stop if some linearity cannot be made to leave basis            */
/*        (b) Permanently remove the cobasic indices of linearities           */
/*        (c) If some decision variable cobasic, it is a linearity,           */
/*            and will be removed.                                            */
{
  int64_t i, j, k;
  /* assign local variables to structures */
  T **A = P->A;
  int64_t *B = P->B;
  int64_t *C = P->C;
  int64_t *Row = P->Row;
  int64_t *Col = P->Col;
  int64_t *linearity = Q->linearity;
  int64_t *redundcol = Q->redundcol;
  int64_t m, d, nlinearity;
  int64_t nredundcol = 0L; /* will be calculated here */
  nlinearity = Q->nlinearity;
  m = P->m;
  d = P->d;

  for (j = 0; j < m; j++) {
    i = 0;
    while (i <= m && B[i] != d + order[j])
      i++;                       /* find leaving basis index i */
    if (j < nlinearity && i > m) /* cannot pivot linearity to cobasis */
      return globals::L_FALSE;
    if (i <= m) { /* try to do a pivot */
      k = 0;
      while (C[k] <= d && A[Row[i]][Col[k]] == 0)
        k++;
      if (C[k] <= d) {
        pivot(P, i, k);
        update(P, &i, &k);
      } else if (j < nlinearity) { /* cannot pivot linearity to cobasis */
        if (A[Row[i]][0] == 0)
          linearity[j] = 0;
        else
          return globals::L_FALSE;
      } /* end if j < nlinearity */
    }   /* end of if i <= m .... */
  }     /* end of for   */
        /* update linearity array to get rid of redundancies */
  i = 0;
  k = 0; /* counters for linearities         */
  while (k < nlinearity) {
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
    if (C[k] <= d) /* decision variable still in cobasis */
      redundcol[nredundcol++] = C[k] - Q->hull; /* adjust for hull indices */
    k++;
  }
  /* now we know how many decision variables remain in problem */
  Q->nredundcol = nredundcol;
  Q->lastdv = d - nredundcol;

  /* Remove linearities from cobasis for rest of computation */
  /* This is done in order so indexing is not screwed up */

  for (i = 0; i < nlinearity; i++) { /* find cobasic index */
    k = 0;
    while (k < d && C[k] != linearity[i] + d)
      k++;
    if (k >= d) {
      return globals::L_FALSE;
    }
    if (!removecobasicindex(P, k))
      return globals::L_FALSE;
    d = P->d;
  }

  /* Check feasability */
  if (Q->givenstart) {
    i = Q->lastdv + 1;
    while (i <= m && A[Row[i]][0] >= 0)
      i++;
    if (i <= m)
      std::cerr << "Infeasible startingcobasis - will be modified";
  }
  return globals::L_TRUE;
} /*  end of getabasis */

template <typename T>
int64_t removecobasicindex(lrs_dic<T> *P, int64_t k)
/* remove the variable C[k] from the problem */
/* used after detecting column dependency    */
{
  int64_t i, j, cindex, deloc;
  /* assign local variables to structures */
  T **A = P->A;
  int64_t *B = P->B;
  int64_t *C = P->C;
  int64_t *Col = P->Col;
  int64_t m, d;
  m = P->m;
  d = P->d;

  cindex = C[k];  /* cobasic index to remove              */
  deloc = Col[k]; /* matrix column location to remove     */

  for (i = 1; i <= m; i++) /* reduce basic indices by 1 after index */
    if (B[i] > cindex)
      B[i]--;

  for (j = k; j < d; j++) { /* move down other cobasic variables    */
    C[j] = C[j + 1] - 1;    /* cobasic index reduced by 1           */
    Col[j] = Col[j + 1];
  }

  if (deloc != d) {
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
  return globals::L_TRUE;
}

template <typename T>
lrs_dic<T> *new_lrs_dic(int64_t m, int64_t d, int64_t m_A) {
  lrs_dic<T> *p = new lrs_dic<T>;

  p->B = new int64_t[m + 1];
  p->Row = new int64_t[m + 1];
  p->C = new int64_t[d + 1];
  p->Col = new int64_t[d + 1];

  p->d_orig = d;

  p->A = new T *[m_A + 1];
  for (int i = 0; i <= m_A; i++)
    p->A[i] = new T[d + 1];
  return p;
}

template <typename T>
lrs_dic<T> *resize(lrs_dic<T> *P, lrs_dat<T> *Q)
/* resize the dictionary after some columns are deleted, ie. inputd>d */
/* a new lrs_dic record is created with reduced size, and items copied over */
{
  lrs_dic<T> *P1; /* to hold new dictionary in case of resizing */

  int64_t i, j;
  int64_t m, d, m_A;

  m = P->m;
  d = P->d;
  m_A = P->m_A;

  /* get new dictionary record */

  P1 = new_lrs_dic<T>(m, d, m_A);

  /* copy data from P to P1    */
  P1->i = P->i;
  P1->j = P->j;
  P1->depth = P->depth;
  P1->m = P->m;
  P1->d = P1->d_orig = d;
  P1->lexflag = P->lexflag;
  P1->m_A = P->m_A;
  P1->det = P->det;

  for (i = 0; i <= m; i++) {
    P1->B[i] = P->B[i];
    P1->Row[i] = P->Row[i];
  }
  for (i = 0; i <= m_A; i++) {
    for (j = 0; j <= d; j++)
      P1->A[i][j] = P->A[i][j];
  }

  for (j = 0; j <= d; j++) {
    P1->Col[j] = P->Col[j];
    P1->C[j] = P->C[j];
  }

  lrs_free_dic(P, Q);

  /* Reassign cache pointers */

  Q->Qhead = P1;
  Q->Qtail = P1;
  P1->next = P1;
  P1->prev = P1;
  return P1;
}

template <typename T>
int64_t restartpivots(lrs_dic<T> *P, lrs_dat<T> *Q)
/* facet contains a list of the inequalities in the cobasis for the restart */
/* inequality contains the relabelled inequalities after initialization     */
{
  int64_t i, j, k;
  int64_t *Cobasic; /* when restarting, Cobasic[j]=1 if j is in cobasis */
                    /* assign local variables to structures */
  T **A = P->A;
  int64_t *B = P->B;
  int64_t *C = P->C;
  int64_t *Row = P->Row;
  int64_t *Col = P->Col;
  int64_t *inequality = Q->inequality;
  int64_t *facet = Q->facet;
  int64_t nlinearity = Q->nlinearity;
  int64_t m, d, lastdv;
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

  Cobasic = new int64_t[m + d + 2];

  /* set Cobasic flags */
  for (i = 0; i < m + d + 1; i++)
    Cobasic[i] = 0;
  for (i = 0; i < d; i++) { /* find index corresponding to facet[i] */
    j = 1;
    while (facet[i + nlinearity] != inequality[j])
      j++;
    Cobasic[j + lastdv] = 1;
  }

  /* Note that the order of doing the pivots is important, as */
  /* the B and C vectors are reordered after each pivot       */

  /* Suggested new code from db starts */
  i = m;
  while (i > d) {
    while (Cobasic[B[i]]) {
      k = d - 1;
      while (k >= 0 && (A[Row[i]][Col[k]] == 0 || Cobasic[C[k]])) {
        k--;
      }
      if (k >= 0) {
        /*db asks: should i really be modified here? (see old code) */
        /*da replies: modifying i only makes is larger, and so      */
        /*the second while loop will put it back where it was       */
        /*faster (and safer) as done below                          */
        int64_t ii = i;
        pivot(P, ii, k);
        update(P, &ii, &k);
      } else {
        delete[] Cobasic;
        return globals::L_FALSE;
      }
    }
    i--;
  }
  /* Suggested new code from db ends */

  /* check restarting from a primal feasible dictionary               */
  for (i = lastdv + 1; i <= m; i++)
    if (A[Row[i]][0] < 0) {
      delete[] Cobasic;
      return globals::L_FALSE;
    }
  delete[] Cobasic;
  return globals::L_TRUE;
}

template <typename T>
int64_t lexmin(lrs_dic<T> *P, lrs_dat<T> *Q, int64_t col)
/*test if basis is lex-min for vertex or ray, if so globals::TRUE */
/* globals::FALSE if a_r,g=0, a_rs !=0, r > s          */
{
  /*do lexmin test for vertex if col=0, otherwise for ray */
  int64_t r, s, i, j;
  /* assign local variables to structures */
  T **A = P->A;
  int64_t *B = P->B;
  int64_t *C = P->C;
  int64_t *Row = P->Row;
  int64_t *Col = P->Col;
  int64_t m = P->m;
  int64_t d = P->d;

  for (i = Q->lastdv + 1; i <= m; i++) {
    r = Row[i];
    if (A[r][col] == 0) /* necessary for lexmin to fail */
      for (j = 0; j < d; j++) {
        s = Col[j];
        if (B[i] > C[j]) {    /* possible pivot to reduce basis */
          if (A[r][0] == 0) { /* no need for ratio test, any pivot feasible */
            if (A[r][s] != 0)
              return globals::L_FALSE;
          } else if (A[r][s] < 0 && ismin(P, r, s)) {
            return globals::L_FALSE;
          }
        } /* end of if B[i] ... */
      }
  }
  return globals::L_TRUE;
}

template <typename T>
int64_t ismin(lrs_dic<T> *P, int64_t r, int64_t s)
/*test if A[r][s] is a min ratio for col s */
{
  int64_t i;
  /* assign local variables to structures */
  T **A = P->A;
  int64_t m_A = P->m_A;

  for (i = 1; i <= m_A; i++)
    if (i != r && A[i][s] < 0 && comprod(A[i][0], A[r][s], A[i][s], A[r][0])) {
      return globals::L_FALSE;
    }

  return globals::L_TRUE;
}

template <typename T>
void update(lrs_dic<T> *P, int64_t *i, int64_t *j)
/*update the B,C arrays after a pivot */
/*   involving B[bas] and C[cob]           */
{

  int64_t leave, enter;
  /* assign local variables to structures */
  int64_t *B = P->B;
  int64_t *C = P->C;
  int64_t *Row = P->Row;
  int64_t *Col = P->Col;
  int64_t m = P->m;
  int64_t d = P->d;

  leave = B[*i];
  enter = C[*j];
  B[*i] = enter;
  reorder1(B, Row, *i, m + 1);
  C[*j] = leave;
  reorder1(C, Col, *j, d);
  /* restore i and j to new positions in basis */
  for (*i = 1; B[*i] != enter; (*i)++)
    ; /*Find basis index */
  for (*j = 0; C[*j] != leave; (*j)++)
    ; /*Find co-basis index */
}

/*********************************************************/
/*                 Miscellaneous                         */
/******************************************************* */
/*reorder array in increasing order with one misplaced element */
void reorder(int64_t a[], int64_t range) {
  int64_t i, temp;
  for (i = 0; i < range - 1; i++)
    if (a[i] > a[i + 1]) {
      temp = a[i];
      a[i] = a[i + 1];
      a[i + 1] = temp;
    }
  for (i = range - 2; i >= 0; i--)
    if (a[i] > a[i + 1]) {
      temp = a[i];
      a[i] = a[i + 1];
      a[i + 1] = temp;
    }
} /* end of reorder */

template <typename T>
int64_t checkredund(lrs_dic<T> *P, lrs_dat<T> *Q)
/* Solve primal feasible lp by least subscript and lex min basis method */
/* to check redundancy of a row in objective function                   */
/* returns globals::TRUE if redundant, else globals::FALSE */
{
  T Ns, Nt;
  int64_t i, j;
  int64_t r, s;

  /* assign local variables to structures */
  T **A = P->A;
  int64_t *Row, *Col;
  int64_t d = P->d;

  Row = P->Row;
  Col = P->Col;

  while (selectpivot(P, Q, &i, &j)) {

    /* sign of new value of A[0][0]            */
    /* is      A[0][s]*A[r][0]-A[0][0]*A[r][s] */

    r = Row[i];
    s = Col[j];
    Ns = A[0][s] * A[r][0];
    Nt = A[0][0] * A[r][s];
    if (Ns > Nt)
      return globals::L_FALSE; /* non-redundant */

    pivot(P, i, j);
    update(P, &i, &j);
  }
  return !(j < d && i == 0); /* unbounded is also non-redundant */
} /* end of checkredund  */

template <typename T>
int64_t checkcobasic(lrs_dic<T> *P, lrs_dat<T> *Q, int64_t index)
/* globals::TRUE if index is cobasic and nonredundant                         */
/* globals::FALSE if basic, or degen. cobasic, where it will get pivoted out  */

{

  /* assign local variables to structures */

  T **A = P->A;
  int64_t *C, *Row, *Col;
  int64_t d = P->d;
  int64_t m = P->m;
  int64_t i = 0;
  int64_t j = 0;
  int64_t s;

  C = P->C;
  Row = P->Row;
  Col = P->Col;

  while ((j < d) && C[j] != index)
    j++;

  if (j == d)
    return globals::L_FALSE; /* not cobasic index */

  /* index is cobasic */

  s = Col[j];
  i = Q->lastdv + 1;

  while (i <= m && (A[Row[i]][s] == 0 || A[Row[i]][0] != 0))
    i++;

  if (i > m)
    return globals::L_TRUE;

  pivot(P, i, j);
  update(P, &i, &j);

  return globals::L_FALSE; /*index is no longer cobasic */
}

/***************************************************************/
/*                                                             */
/*     Routines for caching, allocating etc.                   */
/*                                                             */
/***************************************************************/

/* From here mostly Bremner's handiwork */

template <typename T>
void cache_dict(lrs_dic<T> **D_p, lrs_dat<T> *global, int64_t i, int64_t j,
                uint64_t &dict_count) {
  if (globals::dict_limit > 1) {
    (*D_p)->i = i;
    (*D_p)->j = j;
    pushQ(global, (*D_p)->m, (*D_p)->d, (*D_p)->m_A, dict_count);
    copy_dict(global->Qtail, *D_p);
  }
  *D_p = global->Qtail;
}

template <typename T> void copy_dict(lrs_dic<T> *dest, lrs_dic<T> *src) {
  int64_t m = src->m;
  int64_t m_A = src->m_A; /* number of rows in A */
  int64_t d = src->d;
  int64_t r, s;

  for (r = 0; r <= m_A; r++)
    for (s = 0; s <= d; s++)
      dest->A[r][s] = src->A[r][s];

  dest->i = src->i;
  dest->j = src->j;
  dest->m = m;
  dest->d = d;
  dest->m_A = src->m_A;

  dest->depth = src->depth;
  dest->lexflag = src->lexflag;

  dest->det = src->det;

  for (int u = 0; u < m + 1; u++) {
    dest->B[u] = src->B[u];
    dest->Row[u] = src->Row[u];
  }
  for (int u = 0; u < d + 1; u++) {
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
template <typename T>
void pushQ(lrs_dat<T> *global, int64_t m, int64_t d, int64_t m_A,
           uint64_t &dict_count) {
  if ((global->Qtail->next) == global->Qhead) {
    if (dict_count < globals::dict_limit) {
      lrs_dic<T> *p;
      p = new_lrs_dic<T>(m, d, m_A);
      p->next = global->Qtail->next;
      (global->Qtail->next)->prev = p;
      (global->Qtail->next) = p;
      p->prev = global->Qtail;
      dict_count++;
      global->Qtail = p;
    } else {
      /*
       * user defined limit reached. start overwriting the
       * beginning of Q
       */
      global->Qhead = global->Qhead->next;
      global->Qtail = global->Qtail->next;
    }
  } else {
    global->Qtail = global->Qtail->next;
  }
}

template <typename T> lrs_dic<T> *lrs_getdic(lrs_dat<T> *Q) {
  lrs_dic<T> *p;
  int64_t m;
  m = Q->m;
  if (Q->nonnegative)
    m = m + Q->inputd;
  p = new_lrs_dic<T>(m, Q->inputd, Q->m);

  p->next = p;
  p->prev = p;
  Q->Qhead = p;
  Q->Qtail = p;

  return p;
}

template <typename T> void lrs_free_dic(lrs_dic<T> *P, lrs_dat<T> *Q) {
  /* do the same steps as for allocation, but backwards */
  /* gmp variables cannot be cleared using free: use lrs_clear_mp* */
  lrs_dic<T> *P1;
  int i;

  /* repeat until cache is empty */
  do {
    /* I moved these here because I'm not certain the cached dictionaries
       need to be the same size. Well, it doesn't cost anything to be safe. db
     */

    int64_t m_A = P->m_A;

    for (i = 0; i <= m_A; i++)
      delete[] P->A[i];
    delete[] P->A;

    delete[] P->Row;
    delete[] P->Col;
    delete[] P->C;
    delete[] P->B;

    /* go to next record in cache if any */
    P1 = P->next;
    delete P;
    P = P1;

  } while (Q->Qhead != P);
}

template <typename T> void lrs_free_dat(lrs_dat<T> *Q) {
  /* most of these items were allocated in lrs_alloc_dic */
  delete[] Q->inequality;
  delete[] Q->linearity;
  delete[] Q->facet;
  delete[] Q->redundcol;
  delete[] Q->minratio;
  delete Q;
}

template <typename T>
int64_t check_cache(lrs_dic<T> **D_p, lrs_dat<T> *global, int64_t *i_p,
                    int64_t *j_p) {
  if (global->Qtail == global->Qhead)
    return 0;
  else {
    global->Qtail = global->Qtail->prev;

    *D_p = global->Qtail;

    *i_p = global->Qtail->i;
    *j_p = global->Qtail->j;

    return 1;
  }
}

template <typename T> lrs_dic<T> *lrs_alloc_dic(lrs_dat<T> *Q) {

  lrs_dic<T> *p;
  int64_t i, j;
  int64_t m, d, m_A;
  //  std::cerr << "Q->hull=" << Q->hull << "\n";
  if (Q->hull)        /* d=col dimension of A */
    Q->inputd = Q->n; /* extra column for hull */
  else
    Q->inputd = Q->n - 1;

  m = Q->m;
  d = Q->inputd;
  //  std::cerr << "m=" << m << " d=" << d << "\n";
  m_A = m;

  /* nonnegative flag set means that problem is d rows "bigger"     */
  /* since nonnegative constraints are not kept explicitly          */

  if (Q->nonnegative)
    m = m + d;

  p = new_lrs_dic<T>(m, d, m_A);

  p->next = p;
  p->prev = p;
  Q->Qhead = p;
  Q->Qtail = p;

  /* Initializations */

  p->d = p->d_orig = d;
  p->m = m;
  p->m_A = m_A;
  p->depth = 0L;
  p->lexflag = globals::L_TRUE;
  p->det = 1;

  /*m+d+1 is the number of variables, labelled 0,1,2,...,m+d  */
  /*  initialize array to zero   */
  for (i = 0; i <= m_A; i++)
    for (j = 0; j <= d; j++)
      p->A[i][j] = 0;

  Q->inequality = new int64_t[m + 1];
  for (i = 0; i <= m; i++)
    Q->inequality[i] = 0;
  if (Q->nlinearity == 0) /* linearity may already be allocated */
    Q->linearity = new int64_t[m + 1];
  for (i = 0; i <= m; i++)
    Q->linearity[i] = 0;

  Q->facet = new int64_t[d + 1];
  Q->redundcol = new int64_t[d + 1];
  Q->minratio = new int64_t[m + 1];

  Q->inequality[0] = 2L;

  Q->lastdv = d; /* last decision variable may be decreased */
                 /* if there are redundant columns          */

  /*initialize basis and co-basis indices, and row col locations */
  /*if nonnegative, we label differently to avoid initial pivots */
  /* set basic indices and rows */
  if (Q->nonnegative) {
    for (i = 0; i <= m; i++) {
      p->B[i] = i;
      if (i <= d)
        p->Row[i] = 0; /* no row for decision variables */
      else
        p->Row[i] = i - d;
    }
  } else {
    for (i = 0; i <= m; i++) {
      if (i == 0)
        p->B[0] = 0;
      else
        p->B[i] = d + i;
      p->Row[i] = i;
    }
  }
  for (j = 0; j < d; j++) {
    if (Q->nonnegative)
      p->C[j] = m + j + 1;
    else
      p->C[j] = j + 1;
    p->Col[j] = j + 1;
  }
  p->C[d] = m + d + 1;
  p->Col[d] = 0;
  return p;
}

template <typename T>
void lrs_set_row_mp(lrs_dic<T> *P, lrs_dat<T> *Q, int64_t row, T *num,
                    int64_t ineq)
/* set row of dictionary using num and den arrays for rational input */
/* ineq = 1 (globals::GE)   - ordinary row  */
/*      = 0 (globals::EQ)   - linearity     */
{
  int64_t i, j;

  /* assign local variables to structures */

  T **A;
  int64_t hull;
  int64_t d;
  hull = Q->hull;
  A = P->A;
  d = P->d;
  i = row;
  for (j = hull; j <= d; j++) /* hull data copied to cols 1..d */
    A[i][j] = num[j - hull];
  if (hull)
    A[i][0] = 0;
  if (A[i][hull] != 0)                 /* for H-rep, are zero in column 0     */
    Q->homogeneous = globals::L_FALSE; /* for V-rep, all zero in column 1     */
  if (ineq == globals::EQ)             /* input is linearity */
  {
    Q->linearity[Q->nlinearity] = row;
    Q->nlinearity++;
  }
}

template <typename T>
void lrs_set_obj_mp(lrs_dic<T> *P, lrs_dat<T> *Q, T *num, int64_t max) {
  int64_t i;

  if (max == globals::MAXIMIZE)
    Q->maximize = globals::L_TRUE;
  else {
    Q->minimize = globals::L_TRUE;
    for (i = 0; i <= P->d; i++)
      num[i] = -num[i];
  }
  lrs_set_row_mp(P, Q, 0L, num, globals::GE);
}

template <typename T>
int64_t dan_selectpivot(lrs_dic<T> *P, lrs_dat<T> *Q, int64_t *r, int64_t *s)
/* select pivot indices using dantzig simplex method             */
/* largest coefficient with lexicographic rule to avoid cycling  */
/* Bohdan Kaluzny's handiwork                                    */
/* returns globals::TRUE if pivot found else globals::L_FALSE    */
/* pivot variables are B[*r] C[*s] in locations Row[*r] Col[*s]  */
{
  int64_t j, k, col;
  T coeff;
  /* assign local variables to structures */
  T **A = P->A;
  int64_t *Col = P->Col;
  int64_t d = P->d;

  *r = 0;
  *s = d;
  j = 0;
  k = 0;

  coeff = 0;
  /*find positive cost coef */
  while (k < d) {
    if (A[0][Col[k]] > coeff) {
      j = k;
      coeff = A[0][Col[j]];
    }
    k++;
  }

  if (coeff > 0) {
    *s = j;
    col = Col[j];
    /*find min index ratio */
    *r = lrs_ratio<T>(P, Q, col);
    if (*r != 0)
      return globals::L_TRUE; /* unbounded */
  }
  return globals::L_FALSE;
}

template <typename T>
void fillModelLRS(MyMatrix<T> const &EXT, lrs_dic<T> *P, lrs_dat<T> *Q) {
  int j;
  int iRow, nbRow, nbCol;
  int64_t n;
  int64_t ineq;
  T *num;
  nbRow = EXT.rows();
  nbCol = EXT.cols();
  n = nbCol;
  num = new T[n + 1];
  ineq = 1;
  for (iRow = 0; iRow < nbRow; iRow++) {
    for (j = 0; j < nbCol; ++j)
      num[j] = EXT(iRow, j);
    lrs_set_row_mp(P, Q, iRow + 1, num, ineq);
  }
  for (j = 0; j < nbCol; j++)
    P->A[0][j] = 1;
  delete[] num;
}

template <typename T> void PrintP(lrs_dic<T> *&P, std::string const &message) {
  std::cerr << "message = " << message << "\n";
  std::cerr << "P->A=\n";
  for (int i = 0; i <= P->m_A; i++) {
    for (int j = 0; j <= P->d; j++)
      std::cerr << " " << P->A[i][j];
    std::cerr << "\n";
  }

  std::cerr << "B=";
  for (int i = 0; i <= P->m; i++)
    std::cerr << P->B[i] << " ";
  std::cerr << "\n";
  std::cerr << "Row=";
  for (int i = 0; i <= P->m; i++)
    std::cerr << P->Row[i] << " ";
  std::cerr << "\n";
  //
  std::cerr << "C=";
  for (int i = 0; i <= P->d; i++)
    std::cerr << P->C[i] << " ";
  std::cerr << "\n";
  std::cerr << "Col=";
  for (int i = 0; i <= P->d; i++)
    std::cerr << P->Col[i] << " ";
  std::cerr << "\n";
}

template <typename T>
void initLRS(MyMatrix<T> const &EXT, lrs_dic<T> *&P, lrs_dat<T> *&Q) {
  T **Lin;
  Q = lrs_alloc_dat<T>();
  //  std::cerr << "After lrs_alloc_dat\n";
  if (Q == nullptr) {
    throw TerminalException{1};
  }
  int nbrow = EXT.rows();
  int nbcol = EXT.cols();
  Q->n = nbcol;
  Q->m = nbrow;
  //  std::cerr << "no need to call a lrs_read_dat a priori\n";
  P = lrs_alloc_dic(Q);
  //  PrintP(P, "after lrs_alloc_dic");
  if (P == nullptr) {
    std::cerr << "We failed allocation, let's die\n";
    throw TerminalException{1};
  }
  fillModelLRS(EXT, P, Q);
  //  PrintP(P, "Before lrs_getfirstbasis");

  if (!lrs_getfirstbasis(&P, Q, Lin)) {
    std::cerr << "Error in call to lrs_getfirstbasis\n";
    throw TerminalException{1};
  }
  //  std::cerr << "After lrs_getfirstbasis\n";
}

template <typename T> void freeLRS(lrs_dic<T> *&P, lrs_dat<T> *&Q) {
  lrs_free_dic(P, Q);
  lrs_free_dat(Q);
}


template <typename T>
void set_face(lrs_dic<T> *P, lrs_dat<T> *Q, int const& col, Face & f) {
  int nbRow = Q->m;
  for (int i=0; i<nbRow; i++)
    f[i] = 0;
  for (int i=0; i<P->d; i++) {
    int the_col = P->Col[i];
    if (the_col != col) {
      int idx1 = P->C[i];
      int idx2 = Q->lastdv;
      int idx = Q->inequality[idx1 - idx2] - 1;
      f[idx] = 1;
    }
  }
  for (int i=Q->lastdv+1; i<=P->m; i++) {
    int iRow = P->Row[i];
    if (P->A[iRow][0] == 0) {
      if (col == 0 || P->A[iRow][col] == 0) {
        int idx = iRow - 1;
        f[idx] = 1;
      }
    }
  }
}

template <typename T, typename F>
void Kernel_DualDescription(MyMatrix<T> const &EXT, F const &f) {
  lrs_dic<T> *P;
  lrs_dat<T> *Q;
  int col;
  initLRS(EXT, P, Q);
  T *output = new T[Q->n + 1];
  uint64_t dict_count = 1;
  bool is_first = true;
#ifdef PRINT_LRS_ANALYSIS
  size_t n_entry = 0;
  size_t max_n_error = 0;
#endif
  do {
    for (col = 0; col <= P->d; col++) {
      if (lrs_getsolution(P, Q, output, col)) {
        if (!is_first) {
#ifdef PRINT_LRS_ANALYSIS
          size_t nbCol = EXT.cols();
          size_t nbRow = EXT.rows();
          size_t real_incidence = 0;
          std::cerr << "------------ Entry " << n_entry << " col=" << col << " ------------\n";
          std::cerr << "    ScalProd =";
          Face real_incd(nbRow);
          for (size_t iRow = 0; iRow < nbRow; iRow++) {
            T eScal(0);
            for (size_t iCol = 0; iCol < nbCol; iCol++)
              eScal += output[iCol] * EXT(iRow, iCol);
            if (eScal == 0) {
              std::cerr << " " << iRow;
              real_incidence += 1;
              real_incd[iRow] = 1;
            }
          }
          std::cerr << "\n";
          size_t lrs_incidence = P->d - 1;
          size_t n_error = 0;
          size_t n_correct = 0;
          std::cerr << "    Lrs_Dict =";
          Face lrs_incd(nbRow);
          int max_iRow = std::numeric_limits<int>::min();
          int min_iRow = std::numeric_limits<int>::max();
          for (int i=0; i<P->d; i++) {
            int idx1 = P->C[i];
            int idx2 = Q->lastdv;
            //          std::cerr << "idx1=" << idx1 << " idx2=" << idx2 << "\n";
            int idx = Q->inequality[idx1 - idx2] - 1;
            int the_col = P->Col[i];
            if (the_col != col) {
              std::cerr << " " << idx;
              lrs_incd[idx] = 1;
              if (real_incd[idx] == 1) {
                n_correct++;
              } else {
                n_error++;
              }
            }
          }
          for (int i=Q->lastdv+1; i<=P->m; i++) {
            int iRow = P->Row[i];
            if (iRow < min_iRow)
              min_iRow = iRow;
            if (iRow > max_iRow)
              max_iRow = iRow;
            if (P->A[iRow][0] == 0) {
              if (col == 0 || P->A[iRow][col] == 0) {
                int idx = iRow - 1;
                std::cerr << " " << idx;
                lrs_incidence += 1;
                lrs_incd[idx] = 1;
                if (real_incd[idx] == 1) {
                  n_correct++;
                } else {
                  n_error++;
                }
              }
            }
          }
          if (n_error > max_n_error) {
            max_n_error = n_error;
          }
          std::cerr << "\n";
          std::cerr << "    real_incidence=" << real_incidence << " n_correct=" << n_correct << " n_error=" << n_error << "\n";
          std::cerr << "    min_iRow=" << min_iRow << " max_iRow=" << max_iRow << "\n";
          if (real_incidence != lrs_incidence) {
            std::cerr << "The incidence are different\n";
            std::cerr << "real_incidence=" << real_incidence << " lrs_incidence=" << lrs_incidence << "\n";
            throw TerminalException{1};
          }
          if (real_incd != lrs_incd) {
            std::cerr << "The real_incd is not equal to lrs_incd\n";
            throw TerminalException{1};
          }
          n_entry += 1;
#endif
          f(P, Q, col, output);
        }
        is_first = false;
      }
    }
  } while (lrs_getnextbasis(&P, Q, globals::L_FALSE, dict_count));
#ifdef PRINT_LRS_ANALYSIS
  std::cerr << "max_n_error=" << max_n_error << "\n";
#endif
  delete[] output;
  lrs_free_dic(P, Q);
  lrs_free_dat(Q);
}

template <typename T, typename F>
void Kernel_DualDescription_cond(MyMatrix<T> const &EXT, F const &f) {
  lrs_dic<T> *P;
  lrs_dat<T> *Q;
  int col;
  initLRS(EXT, P, Q);
  T *output = new T[Q->n + 1];
  uint64_t dict_count = 1;
  bool is_first = true;
  do {
    bool is_finished = false;
    for (col = 0; col <= P->d; col++)
      if (lrs_getsolution(P, Q, output, col)) {
        if (!is_first) {
          bool test = f(P, Q, col, output);
          if (!test)
            is_finished = true;
        }
        is_first = false;
      }
    if (is_finished)
      break;
  } while (lrs_getnextbasis(&P, Q, globals::L_FALSE, dict_count));
  delete[] output;
  lrs_free_dic(P, Q);
  lrs_free_dat(Q);
}

template <typename T> MyMatrix<T> FirstColumnZero(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  for (int iRow = 0; iRow < nbRow; iRow++) {
    T eVal = M(iRow, 0);
    if (eVal != 0) {
      MyMatrix<T> Mret(nbRow, nbCol + 1);
      for (int jRow = 0; jRow < nbRow; jRow++) {
        Mret(jRow, 0) = 0;
        for (int iCol = 0; iCol < nbCol; iCol++)
          Mret(jRow, iCol + 1) = M(jRow, iCol);
      }
      return Mret;
    }
  }
  return M;
}

template <typename T>
std::pair<MyMatrix<T>, int> FirstColumnZeroCond(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  for (int iRow = 0; iRow < nbRow; iRow++) {
    T eVal = M(iRow, 0);
    if (eVal != 0) {
      MyMatrix<T> Mret(nbRow, nbCol + 1);
      for (int jRow = 0; jRow < nbRow; jRow++) {
        Mret(jRow, 0) = 0;
        for (int iCol = 0; iCol < nbCol; iCol++)
          Mret(jRow, iCol + 1) = M(jRow, iCol);
      }
      return {Mret, 1};
    }
  }
  return {M, 0};
}

template <typename T> vectface DualDescription_incd(MyMatrix<T> const &EXT) {
  MyMatrix<T> EXTwork = FirstColumnZero(EXT);
  size_t nbRow = EXTwork.rows();
  vectface ListIncd(nbRow);
  T eScal;
#if !defined USE_ISINCD
  Face face(nbRow);
#endif
  auto f = [&](lrs_dic<T> *P, lrs_dat<T> *Q, int const& col, [[maybe_unused]] T *out) -> void {
    set_face(P, Q, col, face);
    ListIncd.push_back(face);
  };
  Kernel_DualDescription(EXTwork, f);
  return ListIncd;
}

template <typename T> MyMatrix<T> DualDescription(MyMatrix<T> const &EXT) {
  std::pair<MyMatrix<T>, int> pair = FirstColumnZeroCond(EXT);
  MyMatrix<T> const &EXTwork = pair.first;
  int shift = pair.second;
  int nbCol = EXTwork.cols();
  int nbColRed = nbCol - shift;
  std::vector<MyVector<T>> ListVect;
  MyVector<T> V(nbColRed);
  auto f = [&]([[maybe_unused]] lrs_dic<T> *P, [[maybe_unused]] lrs_dat<T> *Q, [[maybe_unused]] int const& col, T *out) -> void {
    for (int i = 0; i < nbColRed; i++)
      V(i) = out[i + shift];
    ListVect.push_back(V);
  };
  Kernel_DualDescription(EXTwork, f);
  return MatrixFromVectorFamily(ListVect);
}

template <typename T, typename Fprocess>
void DualDescriptionFaceIneq(MyMatrix<T> const &EXT, Fprocess f_process) {
  std::pair<MyMatrix<T>, int> ePair = FirstColumnZeroCond(EXT);
  MyMatrix<T> const &EXTwork = ePair.first;
  int shift = ePair.second;
  int nbCol = EXTwork.cols();
  int nbRow = EXTwork.rows();
  int nbColRed = nbCol - shift;
  std::pair<Face, MyVector<T>> pair{Face(nbRow), MyVector<T>(nbColRed)};
  T eScal;
  auto f = [&](lrs_dic<T> *P, lrs_dat<T> *Q, int const& col, T *out) -> void {
    for (int i = 0; i < nbColRed; i++)
      pair.second(i) = out[i + shift];
    set_face(P, Q, col, pair.first);
    f_process(pair);
  };
  Kernel_DualDescription(EXTwork, f);
}

template <typename T>
vectface DualDescription_incd_limited(MyMatrix<T> const &EXT,
                                      int const &UpperLimit) {
  MyMatrix<T> EXTwork = FirstColumnZero(EXT);
  size_t nbRow = EXTwork.rows();
  vectface ListIncd(nbRow);
  T eScal;
  int nbFound = 0;
  Face face(nbRow);
  auto f = [&](lrs_dic<T> *P, lrs_dat<T> *Q, int const& col, [[maybe_unused]] T *out) -> bool {
    set_face(P, Q, col, face);
    ListIncd.push_back(face);
    nbFound++;
    return nbFound != UpperLimit;
  };
  Kernel_DualDescription_cond(EXTwork, f);
  return ListIncd;
}

template <typename T>
vectface DualDescription_incd_reduction(MyMatrix<T> const &EXT) {
  MyMatrix<T> EXTwork = FirstColumnZero(EXT);
  using Tring = typename underlying_ring<T>::ring_type;
  size_t nbCol = EXTwork.cols();
  size_t nbRow = EXTwork.rows();
  MyMatrix<Tring> EXTring(nbRow, nbCol);
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    MyVector<T> eRow1 = GetMatrixRow(EXTwork, iRow);
    MyVector<T> eRow2 = NonUniqueScaleToIntegerVector(eRow1);
    MyVector<Tring> eRow3 = UniversalVectorConversion<Tring, T>(eRow2);
    AssignMatrixRow(EXTring, iRow, eRow3);
  }
  vectface ListIncd(nbRow);
  Tring eScal;
  Face face(nbRow);
  auto f = [&](lrs_dic<Tring> *P, lrs_dat<Tring> *Q, int const& col, [[maybe_unused]] Tring *out) -> void {
    set_face(P, Q, col, face);
    ListIncd.push_back(face);
  };
  Kernel_DualDescription(EXTring, f);
  return ListIncd;
}

template <typename T>
MyMatrix<T> DualDescription_reduction(MyMatrix<T> const &EXT) {
  std::pair<MyMatrix<T>, int> pair = FirstColumnZeroCond(EXT);
  MyMatrix<T> const &EXTwork = pair.first;
  int shift = pair.second;
  using Tring = typename underlying_ring<T>::ring_type;
  int nbCol = EXTwork.cols();
  int nbRow = EXTwork.rows();
  int nbColRed = nbCol - shift;
  MyMatrix<Tring> EXTring(nbRow, nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<T> eRow1 = GetMatrixRow(EXTwork, iRow);
    MyVector<T> eRow2 = NonUniqueScaleToIntegerVector(eRow1);
    MyVector<Tring> eRow3 = UniversalVectorConversion<Tring, T>(eRow2);
    AssignMatrixRow(EXTring, iRow, eRow3);
  }
  std::vector<MyVector<T>> ListVect;
  MyVector<T> V(nbColRed);
  auto f = [&]([[maybe_unused]] lrs_dic<Tring> *P, [[maybe_unused]] lrs_dat<Tring> *Q, [[maybe_unused]] int const& col, Tring *out) -> void {
    for (int i = 0; i < nbColRed; i++)
      V(i) = out[i + shift];
    ListVect.push_back(V);
  };
  Kernel_DualDescription(EXTring, f);
  return MatrixFromVectorFamily(ListVect);
}

template <typename T, typename Fprocess>
void DualDescriptionFaceIneq_reduction(MyMatrix<T> const &EXT, Fprocess f_process) {
  std::pair<MyMatrix<T>, int> ePair = FirstColumnZeroCond(EXT);
  MyMatrix<T> const &EXTwork = ePair.first;
  int shift = ePair.second;
  using Tring = typename underlying_ring<T>::ring_type;
  int nbCol = EXTwork.cols();
  int nbRow = EXTwork.rows();
  int nbColRed = nbCol - shift;
  MyMatrix<Tring> EXTring(nbRow, nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<T> eRow1 = GetMatrixRow(EXTwork, iRow);
    MyVector<T> eRow2 = NonUniqueScaleToIntegerVector(eRow1);
    MyVector<Tring> eRow3 = UniversalVectorConversion<Tring, T>(eRow2);
    AssignMatrixRow(EXTring, iRow, eRow3);
  }
  std::pair<Face, MyVector<T>> pair{Face(nbRow), MyVector<T>(nbColRed)};
  Tring eScal;
  auto f = [&](lrs_dic<Tring> *P, lrs_dat<Tring> *Q, int const& col, Tring *out) -> void {
    for (int i = 0; i < nbColRed; i++)
      pair.second(i) = out[i + shift];
    set_face(P, Q, col, pair.first);
    f_process(pair);
  };
  Kernel_DualDescription(EXTring, f);
}

// clang-format off
}  // namespace lrs
#endif  // SRC_POLY_POLY_LRSLIB_H_
// clang-format on

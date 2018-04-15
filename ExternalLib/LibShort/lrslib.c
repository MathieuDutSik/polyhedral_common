/* lrslib.c     library code for lrs                     */

/* last modified: June 19, 2001                          */
/* Copyright: David Avis 2001, avis@cs.mcgill.ca         */

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

#include <stdio.h>
#include <string.h>
#include "lrslib.h"

/* Variables and functions global to this file only */
static long lrs_checkpoint_seconds = 0;

static long lrs_global_count = 0;	/* Track how many lrs_dat records are 
					   allocated */

static lrs_dat_p *lrs_global_list[MAX_LRS_GLOBALS + 1];


static void lrs_dump_state ();


/*******************************/
/* functions  for external use */
/*******************************/

/*******************************************************/
/* lrs_main is driver for lrs.c does H/V enumeration   */
/* showing function calls intended for public use      */
/*******************************************************/
long
lrs_main (int argc, char *argv[])

{

  lrs_dic *P;			/* structure for holding current dictionary and indices */
  lrs_dat *Q;			/* structure for holding static problem data            */

  lrs_mp_vector output;		/* holds one line of output; ray,vertex,facet,linearity */
  lrs_mp_matrix Lin;		/* holds input linearities if any are found             */
  long col;			/* output column index for dictionary                   */
  long startcol = 0;
  long prune = FALSE;		/* if TRUE, getnextbasis will prune tree and backtrack  */

/* global variables lrs_ifp and lrs_ofp are file pointers for input and output   */
/* they default to stdin and stdout, but may be overidden by command line parms. */


/***************************************************
 Step 0: 
  Do some global initialization that should only be done once,
  no matter how many lrs_dat records are allocated. db

***************************************************/

  if ( !lrs_init ("\n*lrs:"))
    return 1;
  printf(AUTHOR);

/*********************************************************************************/
/* Step 1: Allocate lrs_dat, lrs_dic and set up the problem                      */
/*********************************************************************************/

  Q = lrs_alloc_dat ("LRS globals");	/* allocate and init structure for static problem data */

  if (Q == NULL)
    return 1;

  if (!lrs_read_dat (Q, argc, argv))	/* read first part of problem data to get dimensions   */
    return 1;                   	/* and problem type: H- or V- input representation     */

  P = lrs_alloc_dic (Q);	/* allocate and initialize lrs_dic                     */
  if (P == NULL)
    return 1;

  if (!lrs_read_dic (P, Q))	/* read remainder of input to setup P and Q            */
    return 1;

  output = lrs_alloc_mp_vector (Q->n);	/* output holds one line of output from dictionary     */


/*********************************************************************************/
/* Step 2: Find a starting cobasis from default of specified order               */
/*         P is created to hold  active dictionary data and may be cached        */
/*         Lin is created if necessary to hold linearity space                   */
/*         Print linearity space if any, and retrieve output from first dict.    */
/*********************************************************************************/

  if (!lrs_getfirstbasis (&P, Q, &Lin, FALSE))
    return 1;

  /* Pivot to a starting dictionary                      */
  /* There may have been column redundancy               */
  /* If so the linearity space is obtained and redundant */
  /* columns are removed. User can access linearity space */
  /* from lrs_mp_matrix Lin dimensions nredundcol x d+1  */



  if (Q->homogeneous && Q->hull)
    startcol++;			/* col zero not treated as redundant   */

  for (col = startcol; col < Q->nredundcol; col++)	/* print linearity space               */
    lrs_printoutput (Q, Lin[col]);	/* Array Lin[][] holds the coeffs.     */



/*********************************************************************************/
/* Step 3: Terminate if lponly option set, otherwise initiate a reverse          */
/*         search from the starting dictionary. Get output for each new dict.    */
/*********************************************************************************/



  /* We initiate reverse search from this dictionary       */
  /* getting new dictionaries until the search is complete */
  /* User can access each output line from output which is */
  /* vertex/ray/facet from the lrs_mp_vector output         */
  /* prune is TRUE if tree should be pruned at current node */
  do
    {
      for (col = 0; col <= P->d; col++)
	if (lrs_getsolution (P, Q, output, col))
	  lrs_printoutput (Q, output);
    }
  while (!Q->lponly && lrs_getnextbasis (&P, Q, prune));

  if (Q->lponly)
    fprintf (lrs_ofp, "\n*LP solution only requested\n");
  else
    lrs_printtotals (P, Q);	/* print final totals, including estimates       */

  lrs_clear_mp_vector(output, Q->n);
  lrs_free_dic (P,Q);           /* deallocate lrs_dic */

  lrs_free_dat (Q);             /* deallocate lrs_dat */

  lrs_close ("lrs:");

  return 0;
}
/*********************************************/
/* end of model test program for lrs library */
/*********************************************/
/*******************************************************/
/* redund_main is driver for redund.c, removes all     */
/* redundant rows from an H or V-representation        */
/* showing function calls intended for public use      */
/*******************************************************/
long
redund_main (int argc, char *argv[])

{
  lrs_mp_matrix Ain;		/* holds a copy of the input matrix to output at the end */

  long *redineq;		/* redineq[i]=0 if ineq i non-red,1 if red,2 linearity  */
  long ineq;			/* input inequality number of current index             */

  lrs_dic *P;			/* structure for holding current dictionary and indices */
  lrs_dat *Q;			/* structure for holding static problem data            */

  lrs_mp_matrix Lin;		/* holds input linearities if any are found             */

  long i, j, d, m;
  long nlinearity;		/* number of linearities in input file                  */
  long nredund;			/* number of redundant rows in input file               */
  long lastdv;
  long debug;
  long index;			/* basic index for redundancy test */

/* global variables lrs_ifp and lrs_ofp are file pointers for input and output   */
/* they default to stdin and stdout, but may be overidden by command line parms. */
/* Lin is global 2-d array for linearity space if it is found (redund columns)   */

  lrs_ifp = stdin;
  lrs_ofp = stdout;
/***************************************************
 Step 0: 
  Do some global initialization that should only be done once,
  no matter how many lrs_dat records are allocated. db

***************************************************/

  if ( !lrs_init ("\n*redund:"))
    return 1;

  printf (AUTHOR);

/*********************************************************************************/
/* Step 1: Allocate lrs_dat, lrs_dic and set up the problem                      */
/*********************************************************************************/

  Q = lrs_alloc_dat ("LRS globals");	/* allocate and init structure for static problem data */

  if (Q == NULL)
    return 1;

  if (!lrs_read_dat (Q, argc, argv))	/* read first part of problem data to get dimensions   */
    return 1;                         	/* and problem type: H- or V- input representation     */

  P = lrs_alloc_dic (Q);	/* allocate and initialize lrs_dic                     */
  if (P == NULL)
    return 1;

  if (!lrs_read_dic (P, Q))	/* read remainder of input to setup P and Q            */
    return 1;

/* if non-negative flag is set, non-negative constraints are not input */
/* explicitly, and are not checked for redundancy                      */

  m = P->m_A;              /* number of rows of A matrix */   
  d = P->d;
  debug = Q->debug;

  redineq = calloc ((m + 1), sizeof (long));
  Ain = lrs_alloc_mp_matrix (m, d);	/* make a copy of A matrix for output later            */

  for (i = 1; i <= m; i++)
    {
      for (j = 0; j <= d; j++)
	copy (Ain[i][j], P->A[i][j]);

      if (debug)
	lrs_printrow ("*", Q, Ain[i], d);
    }

/*********************************************************************************/
/* Step 2: Find a starting cobasis from default of specified order               */
/*         Lin is created if necessary to hold linearity space                   */
/*********************************************************************************/

  if (!lrs_getfirstbasis (&P, Q, &Lin, TRUE))
    return 1;

  /* Pivot to a starting dictionary                      */
  /* There may have been column redundancy               */
  /* If so the linearity space is obtained and redundant */
  /* columns are removed. User can access linearity space */
  /* from lrs_mp_matrix Lin dimensions nredundcol x d+1  */


/*********************************************************************************/
/* Step 3: Test each row of the dictionary to see if it is redundant             */
/*********************************************************************************/

/* note some of these may have been changed in getting initial dictionary        */
  m = P->m_A;
  d = P->d;
  nlinearity = Q->nlinearity;
  lastdv = Q->lastdv;

/* linearities are not considered for redundancy */

  for (i = 0; i < nlinearity; i++)
    redineq[Q->linearity[i]] = 2L;

/* rows 0..lastdv are cost, decsion variables, or linearities  */
/* other rows need to be tested                                */

  for (index = lastdv + 1; index <= m + d; index++)
    {
      ineq = Q->inequality[index - lastdv];	/* the input inequality number corr. to this index */

      redineq[ineq] = checkindex (P, Q, index);
      if (debug)
	fprintf (lrs_ofp, "\ncheck index=%ld, inequality=%ld, redineq=%ld", index, ineq, redineq[ineq]);
      if (redineq[ineq] == ONE)
	fprintf (lrs_ofp, "\n*row %ld was redundant and removed", ineq);

    }				/* end for index ..... */

  if (debug)
    {
      fprintf (lrs_ofp, "\n*redineq:");
      for (i = 1; i <= m; i++)
	fprintf (lrs_ofp, " %ld", redineq[i]);
    }

  /*
  if (!Q->hull)
    fprintf (lrs_ofp, "\nH-representation");
  else
    fprintf (lrs_ofp, "\nV-representation");
  */

/* linearities will be printed first in output */

  if (nlinearity > 0)
    {
      fprintf (lrs_ofp, "\nlinearity %ld", nlinearity);
      for (i = 1; i <= nlinearity; i++)
	fprintf (lrs_ofp, " %ld", i);

    }
  nredund = nlinearity;		/* count number of non-redundant inequalities */
  for (i = 1; i <= m; i++)
    if (redineq[i] == 0)
      nredund++;
  fprintf (lrs_ofp, "\nbegin");
  fprintf (lrs_ofp, "\n%ld %ld rational", nredund, Q->n);

/* print the linearities first */

  for (i = 0; i < nlinearity; i++)
    lrs_printrow ("", Q, Ain[Q->linearity[i]], Q->inputd);

  for (i = 1; i <= m; i++)
    if (redineq[i] == 0)
      lrs_printrow ("", Q, Ain[i], Q->inputd);
  fprintf (lrs_ofp, "\nend");
  fprintf (lrs_ofp, "\n*Input had %ld rows and %ld columns", m, Q->n);
  fprintf (lrs_ofp, ": %ld row(s) redundant", m - nredund);

  lrs_free_dic (P,Q);           /* deallocate lrs_dic */
  lrs_free_dat (Q);             /* deallocate lrs_dat */

  lrs_close ("redund:");

  return 0;
}
/*********************************************/
/* end of redund.c                           */
/*********************************************/
/*******************/
/* lrs_printoutput */
/*******************/
void 
lrs_printoutput (lrs_dat * Q, lrs_mp_vector output)
{
  long i;

  fprintf (lrs_ofp, "\n");
  if (Q->hull || zero (output[0]))	/*non vertex */
    {
      for (i = 0; i < Q->n; i++)
	pmp ("", output[i]);

    }
  else
    {				/* vertex   */
      fprintf (lrs_ofp, " 1 ");
      for (i = 1; i < Q->n; i++)
	prat ("", output[i], output[0]);
    }
  fflush(lrs_ofp);

}
/**************************/
/* end of lrs_printoutput */
/**************************/

void 
lrs_printrow (char name[], lrs_dat * Q, lrs_mp_vector output, long rowd)
/* print a row of A matrix in output in "original" form  */
/* rowd+1 is the dimension of output vector                */
/* if input is H-rep. output[0] contains the RHS      */
/* if input is V-rep. vertices are scaled by 1/output[1] */
{
  long i;

  fprintf (lrs_ofp, "\n%s", name);
  if (!Q->hull)			/* input was inequalities, print directly */

    {

      for (i = 0; i <= rowd; i++)
	pmp ("", output[i]);
      return;
    }

/* input was vertex/ray */

  if (zero (output[1]))		/*non-vertex */
    {
      for (i = 1; i <= rowd; i++)
	pmp ("", output[i]);

    }
  else
    {				/* vertex */
      fprintf (lrs_ofp, " 1 ");
      for (i = 2; i <= rowd; i++)
	prat ("", output[i], output[1]);
    }

  return;

}				/* end of lrs_printrow */

long 
lrs_getsolution (lrs_dic * P, lrs_dat * Q, lrs_mp_vector output, long col)
   /* check if column indexed by col in this dictionary */
   /* contains output                                   */
   /* col=0 for vertex 1....d for ray/facet             */
{
  long j;			/* cobasic index     */

  lrs_mp_matrix A = P->A;
  long *Row = P->Row;

  if (col == ZERO)		/* check for lexmin vertex */
    return lrs_getvertex (P, Q, output);

/*  check for rays: negative in row 0 , positive if lponly */

  if (Q->lponly)
    {
      if (!positive (A[0][col]))
	return FALSE;
    }

  else if (!negative (A[0][col]))
    return FALSE;


/*  and non-negative for all basic non decision variables */

  j = Q->lastdv + 1;
  while (j <= P->m && !negative (A[Row[j]][col]))
    j++;

  if (j <= P->m)
    return FALSE;

  if (Q->geometric || Q->allbases || lexmin (P, Q, col) || Q->lponly)

    return lrs_getray (P, Q, col, Q->n, output);

  return FALSE;			/* no more output in this dictionary */

}				/* end of lrs_getsolution */


long
lrs_init (char *name)       /* returns TRUE if successful, else FALSE */
{
  /*
  printf ("%s", name);
  printf (TITLE);
  printf (VERSION);
  printf ("(");
  printf (BIT);
  printf (",");
  printf (ARITH);
  */
  if (!lrs_mp_init (ZERO, stdin, stdout))  /* initialize arithmetic */
    return FALSE;
  /*
  printf (")");
  */


  lrs_global_count = 0;
  lrs_checkpoint_seconds = 0;
#ifdef SIGNALS
  setup_signals ();
#endif
  return TRUE;
}

void 
lrs_close (char *name)
{

  /*
  fprintf (lrs_ofp, "\n*%s", name);
  fprintf (lrs_ofp, TITLE);
  fprintf (lrs_ofp, VERSION);
  fprintf (lrs_ofp, "(");
  fprintf (lrs_ofp, BIT);
  fprintf (lrs_ofp, ",");
  fprintf (lrs_ofp, ARITH);
  fprintf (lrs_ofp, ")");

#ifdef MP   
  fprintf (lrs_ofp, " max digits=%ld/%ld", DIG2DEC (lrs_record_digits), DIG2DEC (lrs_digits));
#endif

#ifdef TIMES
  ptimes ();
#endif
  */
  fclose (lrs_ifp);
  if (lrs_ofp != stdout)
    fclose (lrs_ofp);
}

/***********************************/
/* allocate and initialize lrs_dat */
/***********************************/
lrs_dat *
lrs_alloc_dat (char *name)
{
  lrs_dat *Q;
  long i;


  if (lrs_global_count >= MAX_LRS_GLOBALS)
    {
      fprintf (stderr,
	       "Fatal: Attempt to allocate more than %ld global data blocks\n", MAX_LRS_GLOBALS);
      exit (1);

    }

  Q = (lrs_dat *) malloc (sizeof (lrs_dat));
  if (Q == NULL)
    return Q;			/* failure to allocate */

  lrs_global_list[lrs_global_count] = Q;
  Q->id = lrs_global_count;
  lrs_global_count++;
  Q->name=(char *) CALLOC ((unsigned) strlen(name)+1, sizeof (char));
  strcpy(Q->name,name); 

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
  for (i = 0; i < 4; i++)
    {
      Q->count[i] = 0L;
      Q->cest[i] = 0.0;
    }
  Q->count[2] = 1L;		/* basis counter */
/* initialize flags */
  Q->allbases = FALSE;
  Q->debug = FALSE;
  Q->frequency = 0L;
  Q->dualdeg = TRUE;
  Q->geometric = FALSE;
  Q->getvolume = FALSE;
  Q->homogeneous = TRUE;
  Q->polytope = FALSE;
  Q->hull = FALSE;
  Q->incidence = FALSE;
  Q->lponly = FALSE;
  Q->maxdepth = MAXD;
  Q->mindepth = -MAXD;
  Q->nonnegative = FALSE;
  Q->printcobasis = FALSE;
  Q->truncate = FALSE;          /* truncate tree when moving from opt vertex        */
  Q->verbose=FALSE;
  Q->voronoi = FALSE;
  Q->maximize = FALSE;		/*flag for LP maximization                          */
  Q->minimize = FALSE;		/*flag for LP minimization                          */
  Q->restart = FALSE;		/* TRUE if restarting from some cobasis             */
  Q->givenstart = FALSE;	/* TRUE if a starting cobasis is given              */
  Q->strace = -1L;		/* turn on  debug at basis # strace */
  Q->etrace = -1L;		/* turn off debug at basis # etrace */

  Q->saved_flag = 0;		/* no cobasis saved initially, db */
  lrs_alloc_mp (Q->Nvolume);
  lrs_alloc_mp (Q->Dvolume);
  lrs_alloc_mp (Q->sumdet);
  lrs_alloc_mp (Q->objnum);
  lrs_alloc_mp (Q->objden);
  lrs_alloc_mp (Q->saved_det);
  itomp (ZERO, Q->Nvolume);
  itomp (ONE, Q->Dvolume);
  itomp (ZERO, Q->objnum);
  itomp (ONE, Q->objden);
  itomp (ZERO, Q->sumdet);

  return Q;
}				/* end of allocate and initialize lrs_dat */

/*******************************/
/*  lrs_read_dat               */
/*******************************/
long 
lrs_read_dat (lrs_dat * Q, int argc, char *argv[])
{
  char name[100];
  long dec_digits = 0;
  long firstline = TRUE;	/*flag for picking off name at line 1               */

  char c;			/* for fgetc */


  if (argc > 1)			/* command line argument overides stdin        */
    {
      if ((lrs_ifp = fopen (argv[1], "r")) == NULL)
	{
	  printf ("\nBad input file name");
	  return (FALSE);
	}
      else
	printf ("\n*Input taken from file %s", argv[1]);
    }

  if (argc > 2)			/* command line argument overides stdin        */
    {
      if ((lrs_ofp = fopen (argv[2], "a")) == NULL)
	{
	  printf ("\nBad output file name");
	  return (FALSE);
	}
      else
	printf ("\n*Output sent to file %s\n", argv[2]);
    }


/* process input file */
  fscanf (lrs_ifp, "%s", name);

  while (strcmp (name, "begin") != 0)	/*skip until "begin" found processing options */
    {
      if (strncmp (name, "*", 1) == 0)	/* skip any line beginning with * */
	{
	  c = name[0];
	  while (c != EOF && c != '\n')
	    c = fgetc (lrs_ifp);
	}

      else if (strcmp (name, "H-representation") == 0)
	Q->hull = FALSE;
      else if ((strcmp (name, "hull") == 0) || (strcmp (name, "V-representation") == 0))
	Q->hull = TRUE;
      else if (strcmp (name, "digits") == 0)
	{
	  if (fscanf (lrs_ifp, "%ld", &dec_digits) == EOF)
	    {
	      fprintf (lrs_ofp, "\nNo begin line");
	      return (FALSE);
	    }
          if  (!lrs_set_digits(dec_digits))
             return (FALSE);
	}
      else if (strcmp (name, "linearity") == 0)
	{
	  if (!readlinearity (Q))
	    return FALSE;
	}
      else if (strcmp (name, "nonnegative") == 0)
	{
	  Q->nonnegative = TRUE;
	}
      else if (firstline)
	{
	  stringcpy (Q->fname, name);
	  fprintf (lrs_ofp, "\n%s", Q->fname);
	  firstline = FALSE;
	}

      if (fscanf (lrs_ifp, "%s", name) == EOF)
	{
	  fprintf (lrs_ofp, "\nNo begin line");
	  return (FALSE);
	}

    }				/* end of while */


  if (fscanf (lrs_ifp, "%ld %ld %s", &Q->m, &Q->n, name) == EOF)
    {
      fprintf (lrs_ofp, "\nNo data in file");
      return (FALSE);
    }
  if (strcmp (name, "integer") != 0 && strcmp (name, "rational") != 0)
    {
      fprintf (lrs_ofp, "\nData type must be integer of rational");
      return (FALSE);
    }


  if (Q->m == 0)
    {
      fprintf (lrs_ofp, "\nNo input given");	/* program dies ungracefully */
      return (FALSE);
    }



  /* inputd may be reduced in preprocessing of linearities and redund cols */

  return TRUE;
}				/* end of lrs_read_dat */
/****************************/
/* set up lrs_dic structure */
/****************************/
long 
lrs_read_dic (lrs_dic * P, lrs_dat * Q)
/* read constraint matrix and set up problem and dictionary  */

{
  lrs_mp Temp, mpone;
  lrs_mp_vector oD;		/* Denom for objective function */

  long i, j;
  char name[100];
  char c;

/* assign local variables to structures */

  lrs_mp_matrix A;
  lrs_mp_vector Gcd, Lcm;
  long hull = Q->hull;
  long m, d;
  lrs_alloc_mp(Temp); lrs_alloc_mp(mpone);
  A = P->A;
  m = Q->m;
  d = Q->inputd;
  Gcd = Q->Gcd;
  Lcm = Q->Lcm;
  Gcd = Q->Gcd;
  Lcm = Q->Lcm;

  oD = lrs_alloc_mp_vector (d);

  Q->polytope = hull;		/* will be updated as input read                   */
  itomp (ONE, mpone);
  itomp (ONE, A[0][0]);
  itomp (ONE, Lcm[0]);
  itomp (ONE, Gcd[0]);

  for (i = 1; i <= m; i++)	/* read in input matrix row by row                 */
    {
      itomp (ONE, Lcm[i]);	/* Lcm of denominators */
      itomp (ZERO, Gcd[i]);	/* Gcd of numerators */
      for (j = hull; j <= d; j++)	/* hull data copied to cols 1..d */
	{
	  if (readrat (A[i][j], A[0][j]))
	    lcm (Lcm[i], A[0][j]);	/* update lcm of denominators */
	  copy (Temp, A[i][j]);
	  gcd (Gcd[i], Temp);	/* update gcd of numerators   */
	}

      if (hull)
	{
	  itomp (ZERO, A[i][0]);	/*for hull, we have to append an extra column of zeroes */
	  if (!one (A[i][1]) || !one (A[0][1]))		/* all rows must have a one in column one */
	    Q->polytope = FALSE;
	}
      if (!zero (A[i][hull]))	/* for H-rep, are zero in column 0     */
	Q->homogeneous = FALSE;	/* for V-rep, all zero in column 1     */

      storesign (Gcd[i], POS);
      storesign (Lcm[i], POS);
      if (greater (Gcd[i], mpone) || greater (Lcm[i], mpone))
	for (j = 0; j <= d; j++)
	  {
	    exactdivint (A[i][j], Gcd[i], Temp);	/*reduce numerators by Gcd  */
	    mulint (Lcm[i], Temp, Temp);	/*remove denominators */
	    exactdivint (Temp, A[0][j], A[i][j]);	/*reduce by former denominator */
	  }

    }				/* end of for i=       */


  if (Q->homogeneous && Q->verbose)
    {
      fprintf (lrs_ofp, "\n*Input is homogeneous, column 1 not treated as redundant");
    }


/* read in flags */
  while (fscanf (lrs_ifp, "%s", name) != EOF)
    {
      if (strncmp (name, "*", 1) == 0)	/* skip any line beginning with * */
	{
	  c = name[0];
	  while (c != EOF && c != '\n')
	    c = fgetc (lrs_ifp);
	}


      if (strcmp (name, "checkpoint") == 0)
	{
	  long seconds;

	  fscanf (lrs_ifp, "%ld", &seconds);
#ifdef SIGNALS
	  if (seconds > 0)
	    {
	      lrs_checkpoint_seconds = seconds;
	      errcheck ("signal", signal (SIGALRM, timecheck));
	      alarm (lrs_checkpoint_seconds);
	    }
#endif
	}

      if (strcmp (name, "debug") == 0)
	{
	  fscanf (lrs_ifp, "%ld %ld", &Q->strace, &Q->etrace);
	  fprintf (lrs_ofp, "\n*%s from B#%ld to B#%ld", name, Q->strace, Q->etrace);
          Q->verbose=TRUE;
	  if (Q->strace <= 1)
	    Q->debug = TRUE;
	}
      if (strcmp (name, "startingcobasis") == 0)
	{
          if(Q->nonnegative)
	      fprintf (lrs_ofp, "\n*startingcobasis incompatible with nonnegative option:skipped");
          else
            {    
	      fprintf (lrs_ofp, "\n*startingcobasis");
	      Q->givenstart = TRUE;
	      if (!readfacets (Q, Q->inequality))
	          return FALSE;
	    }
        }

      if (strcmp (name, "restart") == 0)
	{
	  Q->restart = TRUE;
          if(Q->voronoi)
           {
             fscanf (lrs_ifp, "%ld %ld %ld %ld", &Q->count[1], &Q->count[0], &Q->count[2], &P->depth);
             fprintf (lrs_ofp, "\n*%s V#%ld R#%ld B#%ld h=%ld data points", name, Q->count[1], Q->count[0], Q->count[2], P->depth);
            }
          else if(hull)
            {
	     fscanf (lrs_ifp, "%ld %ld %ld", &Q->count[0], &Q->count[2], &P->depth);
	     fprintf (lrs_ofp, "\n*%s F#%ld B#%ld h=%ld vertices/rays", name, Q->count[0], Q->count[2], P->depth);
            }
          else
            {
	     fscanf (lrs_ifp, "%ld %ld %ld %ld", &Q->count[1], &Q->count[0], &Q->count[2], &P->depth);
	     fprintf (lrs_ofp, "\n*%s V#%ld R#%ld B#%ld h=%ld facets", name, Q->count[1], Q->count[0], Q->count[2], P->depth);
            }
	  if (!readfacets (Q, Q->facet))
	    return FALSE;
	}			/* end of restart */

/* The next flag request a LP solution only */
      if (strcmp (name, "lponly") == 0)
	{
	  if (Q->hull)
	    fprintf (lrs_ofp, "\n*lponly  option not valid for V-representation-skipped");
	  else
	    Q->lponly = TRUE;
	}


/* The LP will be solved after initialization to get starting vertex   */
/* Used also with lponly flag                                          */
      if (strcmp (name, "maximize") == 0 || strcmp (name, "minimize") == 0)
	{
	  if (Q->hull)
	    fprintf (lrs_ofp, "\n*%s option not valid for V-representation-skipped", name);
	  else
	    {
	      {
		if (strcmp (name, "maximize") == 0)
		  Q->maximize = TRUE;
		else
		  Q->minimize = TRUE;
	      }
	      printf ("\n*%s", name);

	      for (j = 0; j <= d; j++)
		{
		  if (readrat (A[0][j], oD[j]))
		    {
		      reduce (A[0][j], oD[j]);
		      lcm (Q->Lcm[0], oD[j]);	/* update lcm of denominators */
		    }
		  prat ("", A[0][j], oD[j]);
		  if (!Q->maximize)
		    changesign (A[0][j]);
		}
	      storesign (Q->Lcm[0], POS);
	      if (greater (Q->Lcm[0], mpone))
		for (j = 0; j <= d; j++)
		  {
		    mulint (Q->Lcm[0], A[0][j], A[0][j]);	/*remove denominators */
		    copy (Temp, A[0][j]);
		    exactdivint (Temp, oD[j], A[0][j]);
		  }
	      if (Q->debug)
		printA (P, Q);
	    }
	}			/* end of LP setup */
      if (strcmp (name, "volume") == 0)
	{
	  fprintf (lrs_ofp, "\n*%s", name);
	  Q->getvolume = TRUE;
	}
      if (strcmp (name, "geometric") == 0)
	{
	  fprintf (lrs_ofp, "\n*%s", name);
	  if (hull & !Q->voronoi)
	    fprintf (lrs_ofp, " - option for H-representation or voronoi only, skipped");
	  else
	    Q->geometric = TRUE;
	}
      if (strcmp (name, "allbases") == 0)
	{
	  fprintf (lrs_ofp, "\n*%s", name);
	  Q->allbases = TRUE;
	}

      if (strcmp (name, "incidence") == 0)
	{
	  fprintf (lrs_ofp, "\n*%s", name);
	  Q->incidence = TRUE;
	}

      if (strcmp (name, "#incidence") == 0) /* number of incident inequalities only */
	{
	  Q->printcobasis = TRUE;
	}

      if (strcmp (name, "printcobasis") == 0)
	{
	  fscanf (lrs_ifp, "%ld", &Q->frequency);
	  fprintf (lrs_ofp, "\n*%s", name);
          if (Q->frequency > 0)
            fprintf(lrs_ofp," %ld", Q->frequency);
	  Q->printcobasis = TRUE;
	}

      if (strcmp (name, "cache") == 0)
	{
	  fscanf (lrs_ifp, "%ld", &dict_limit);
	  fprintf (lrs_ofp, "\n*cache %ld", dict_limit);
	  if (dict_limit < 1)
	    dict_limit = 1;
	}
      if (strcmp (name, "linearity") == 0)
	{
	  if (!readlinearity (Q))
	    return FALSE;
	}

      if (strcmp (name, "maxdepth") == 0)
	{
	  fscanf (lrs_ifp, "%ld", &Q->maxdepth);
	  fprintf (lrs_ofp, "\n*%s  %ld", name, Q->maxdepth);
	}
      if (strcmp (name, "mindepth") == 0)
	{
	  fscanf (lrs_ifp, "%ld", &Q->mindepth);
	  fprintf (lrs_ofp, "\n*%s  %ld", name, Q->mindepth);
	}
      if (strcmp (name, "truncate") == 0)
        {
          fprintf (lrs_ofp, "\n*%s", name);
          if (!hull)
            Q->truncate = TRUE;
          else
            fprintf (lrs_ofp, " - option for H-representation only, skipped");
        }


      if (strcmp (name, "verbose") == 0)
          Q->verbose = TRUE;

      if (strcmp (name, "nonnegative") == 0)
	{
	  fprintf (lrs_ofp, "\n*%s", name);
	  fprintf (lrs_ofp, " - option must come before begin line - skipped");
	}

      if (strcmp (name, "seed") == 0)
	{
	  fscanf (lrs_ifp, "%ld", &Q->seed);
	  fprintf (lrs_ofp, "\n*seed= %ld ", Q->seed);
	}

      if (strcmp (name, "estimates") == 0)
	{
	  fscanf (lrs_ifp, "%ld", &Q->runs);
	  fprintf (lrs_ofp, "\n*%ld %s", Q->runs, name);
	}

      if ((strcmp (name, "voronoi") == 0) || (strcmp (name, "Voronoi") == 0))
	{
	  if (!hull)
	    fprintf (lrs_ofp, "\n*voronoi requires V-representation - option skipped");
	  else
	    {
	      Q->voronoi = TRUE;
	      Q->polytope = FALSE;
	    }
	}

    }				/* end of while for reading flags */

  if (Q->polytope)
    Q->getvolume = TRUE;	/* might as well get volume, it doesn't cost much */


/* Certain options are incompatible, this is fixed here */

  if (Q->restart)
    Q->getvolume = FALSE;       /* otherwise incorrect volume reported            */

  if (Q->incidence)
    {
      Q->printcobasis = TRUE;
      Q->frequency    = ZERO;
    }

  if (Q->debug)
    {
      printA (P, Q);
      fprintf (lrs_ofp, "\nexiting lrs_read_dic");
    }
  lrs_clear_mp(Temp); lrs_clear_mp(mpone);
  lrs_clear_mp_vector (oD,d);
  return TRUE;

}
/********* end of lrs_read_dic ***************/

/* In lrs_getfirstbasis and lrs_getnextbasis we use D instead of P */
/* since the dictionary P may change, ie. &P in calling routine    */

#define D (*D_p)

long 
lrs_getfirstbasis (lrs_dic ** D_p, lrs_dat * Q, lrs_mp_matrix * Lin, long no_output)
/* gets first basis, FALSE if none              */
/* P may get changed if lin. space Lin found    */
/* no_output is TRUE supresses output headers   */
{
  lrs_mp scale, Temp;

  long i, j, k;

/* assign local variables to structures */

  lrs_mp_matrix A;
  long *B, *C, *Row, *Col;
  long *inequality;
  long *linearity;
  long hull = Q->hull;
  long m, d, lastdv, nlinearity, nredundcol;

  lrs_alloc_mp(Temp); lrs_alloc_mp(scale);

  if (Q->lponly)
    no_output = TRUE;
  m = D->m;
  d = D->d;
  lastdv = Q->lastdv;

  nredundcol = 0L;		/* will be set after getabasis        */
  nlinearity = Q->nlinearity;	/* may be reset if new linearity read */
  linearity = Q->linearity;

  A = D->A;
  B = D->B;
  C = D->C;
  Row = D->Row;
  Col = D->Col;
  inequality = Q->inequality;


  if (Q->nlinearity > 0 && Q->nonnegative)
   {
    fprintf (lrs_ofp, "\n*linearity and nonnegative options incompatible");
    fprintf (lrs_ofp, " - all linearities are skipped");
    fprintf (lrs_ofp, "\n*add nonnegative constraints explicitly and ");
    fprintf (lrs_ofp, " remove nonnegative option");
   }

  if (Q->nlinearity && Q->voronoi)
    fprintf (lrs_ofp, "\n*linearity and Voronoi options set - results unpredictable");
  if (Q->lponly && !Q->maximize && !Q->minimize)
    fprintf (lrs_ofp, "\n*LP has no objective function given - assuming all zero");


  if (Q->runs > 0)		/* arrays for estimator */
    {
      Q->isave = (long *) CALLOC ((unsigned) (m * d), sizeof (long));
      Q->jsave = (long *) CALLOC ((unsigned) (m * d), sizeof (long));
    }
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
  if (Q->debug)
    {
      fprintf (lrs_ofp, "\n*Starting cobasis uses input row order");
      for (i = 0; i < m; i++)
	fprintf (lrs_ofp, " %ld", inequality[i]);
    }
/* for voronoi convert to h-description using the transform                  */
/* a_0 .. a_d-1 -> (a_0^2 + ... a_d-1 ^2)-2a_0x_0-...-2a_d-1x_d-1 + x_d >= 0 */
/* note constant term is stored in column d, and column d-1 is all ones      */
/* the other coefficients are multiplied by -2 and shifted one to the right  */
  if (Q->debug)
    printA (D, Q);
  if (Q->voronoi)
    {
      Q->hull = FALSE;
      hull = FALSE;
      for (i = 1; i <= m; i++)
	{
	  if (zero (A[i][1]))
	    {
	      fprintf (lrs_ofp, "\nWith voronoi option column one must be all one");
	      return (FALSE);
	    }
	  copy (scale, A[i][1]);	/*adjust for scaling to integers of rationals */
	  itomp (ZERO, A[i][0]);
	  for (j = 2; j <= d; j++)	/* transform each input row */
	    {
	      copy (Temp, A[i][j]);
	      mulint (A[i][j], Temp, Temp);
	      linint (A[i][0], ONE, Temp, ONE);
	      linint (A[i][j - 1], ZERO, A[i][j], -TWO);
	      mulint (scale, A[i][j - 1], A[i][j - 1]);
	    }			/* end of for (j=1;..) */
	  copy (A[i][d], scale);
	  mulint (scale, A[i][d], A[i][d]);
	}			/* end of for (i=1;..) */
      if (Q->debug)
	printA (D, Q);
    }				/* end of if(voronoi)     */
  if (!Q->maximize && !Q->minimize)
    for (j = 0; j <= d; j++)
      itomp (ZERO, A[0][j]);

/* Now we pivot to standard form, and then find a primal feasible basis       */
/* Note these steps MUST be done, even if restarting, in order to get         */
/* the same index/inequality correspondance we had for the original prob.     */
/* The inequality array is used to give the insertion order                   */
/* and is defaulted to the last d rows when givenstart=FALSE                  */

  if(Q->nonnegative) 
   {
/* no need for initial pivots here, labelling already done */
     Q->lastdv = d;
     Q->nredundcol = 0;
   }
  else
  {
     if (!getabasis (D, Q, inequality))
          return FALSE;
  }
  if(Q->debug)
  {
    printf("\nafter getabasis");
    printA(D, Q);
  }
  nredundcol = Q->nredundcol;
  lastdv = Q->lastdv;
  d = D->d;

/********************************************************************/
/* now we start printing the output file  unless no output requested */
/********************************************************************/
  if (!no_output || Q->debug)
    {
      /*
      if (Q->voronoi)
	fprintf (lrs_ofp, "\n*Voronoi Diagram: Voronoi vertices and rays are output");
      if (hull)
	fprintf (lrs_ofp, "\nH-representation");
      else
	fprintf (lrs_ofp, "\nV-representation");
      */

/* Print linearity space                 */
/* Don't print linearity if first column zero in hull computation */

      if (hull && Q->homogeneous)
	k = 1;			/* 0 normally, 1 for homogeneous case     */
      else
	k = 0;

      if (nredundcol > k)
	{
	  fprintf (lrs_ofp, "\nlinearity %ld ", nredundcol - k);	/*adjust nredundcol for homog. */
	  for (i = 1; i <= nredundcol - k; i++)
	    fprintf (lrs_ofp, " %ld", i);
	}			/* end print of linearity space */
      /*
      fprintf (lrs_ofp, "\nbegin");
      fprintf (lrs_ofp, "\n***** %ld rational", Q->n);
      */

    }				/* end of if !no_output .......   */

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
  if (Q->debug)
    {
      fprintf (lrs_ofp, "\ninequality array initialization:");
      for (i = 1; i <= m - nlinearity; i++)
	fprintf (lrs_ofp, " %ld", inequality[i]);
    }
  if (nredundcol > 0)
    {
      *Lin = lrs_alloc_mp_matrix (nredundcol, Q->n);

      for (i = 0; i < nredundcol; i++)
	{
	  if (!(Q->homogeneous && Q->hull && i == 0))	/* skip redund col 1 for homog. hull */
	    {
	      lrs_getray (D, Q, Col[0], D->C[0] + i - hull, (*Lin)[i]);		/* adjust index for deletions */
	    }

	  if (!removecobasicindex (D, Q, 0L))
	    return FALSE;
	}
    }				/* end if nredundcol > 0 */

  if (Q->lponly)
      if (Q->verbose)
      {
      fprintf (lrs_ofp, "\nNumber of pivots for starting dictionary: %ld",Q->count[3]);
                     printA (D, Q);
       }

/* Do dual pivots to get primal feasibility */
  if (!primalfeasible (D, Q))
    {
      fprintf (lrs_ofp, "\nNo feasible solution");
      return FALSE;
    }

  if (Q->lponly)
      if (Q->verbose)
     {
      fprintf (lrs_ofp, "\nNumber of pivots for feasible solution: %ld",Q->count[3]);
      printA (D, Q);
     }


/* Now solve LP if objective function was given */
  if (Q->maximize || Q->minimize)
    {
      lrs_solvelp (D, Q, Q->maximize);
      if (Q->lponly)		
        {
	  /*
         fprintf (lrs_ofp, "\nNumber of pivots for optimum solution: %ld",Q->count[3]);
	  */
         if (Q->verbose)
         {
           printA (D, Q);
          }
          lrs_clear_mp(Temp); lrs_clear_mp(scale);
          return TRUE;
        }

      else                         /* check to see if objective is dual degenerate */
       {
	  j = 1;
	  while (j <= d && !zero (A[0][j]))
	    j++;
	  if (j <= d)
	    Q->dualdeg = FALSE;
	}
    }
  else
/* re-initialize cost row to -det */
    {
      for (j = 1; j <= d; j++)
	{
	  copy (A[0][j], D->det);
	  storesign (A[0][j], NEG);
	}

      itomp (ZERO, A[0][0]);	/* zero optimum objective value */
    }


/* reindex basis to 0..m if necessary */
/* we use the fact that cobases are sorted by index value */
  if (Q->debug)
    printA (D, Q);
  while (C[0] <= m)
    {
      i = C[0];
      j = inequality[B[i] - lastdv];
      inequality[B[i] - lastdv] = inequality[C[0] - lastdv];
      inequality[C[0] - lastdv] = j;
      C[0] = B[i];
      B[i] = i;
      reorder1 (C, Col, ZERO, d);
    }

  if (Q->debug)
    {
      fprintf (lrs_ofp, "\n*Inequality numbers for indices %ld .. %ld : ", lastdv + 1, m + d);
      for (i = 1; i <= m - nlinearity; i++)
	fprintf (lrs_ofp, " %ld ", inequality[i]);
      printA (D, Q);
    }



  if (Q->restart)
    {
      if (Q->debug)
	fprintf (lrs_ofp, "\nPivoting to restart co-basis");
      if (!restartpivots (D, Q))
	return FALSE;
      D->lexflag = lexmin (D, Q, ZERO);		/* see if lexmin basis */
      if (Q->debug)
	printA (D, Q);
    }
/* Check to see if necessary to resize */
  if (Q->inputd > D->d)
    *D_p = resize (D, Q);

  lrs_clear_mp(Temp); lrs_clear_mp(scale);
  return TRUE;
}
/********* end of lrs_getfirstbasis  ***************/


/*****************************************/
/* getnextbasis in reverse search order  */
/*****************************************/


long 
lrs_getnextbasis (lrs_dic ** D_p, lrs_dat * Q, long backtrack)
	 /* gets next reverse search tree basis, FALSE if none  */
	 /* switches to estimator if maxdepth set               */
	 /* backtrack TRUE means backtrack from here            */

{
  /* assign local variables to structures */
  long i = 0L, j = 0L;
  long m = D->m;
  long d = D->d;

  while ((j < d) || (D->B[m] != m))	/*main while loop for getnextbasis */
    {
      if (D->depth >= Q->maxdepth)
	{
	  backtrack = TRUE;
	  if (Q->runs > 0)	/*get an estimate of remaining tree */
	    lrs_estimate (D, Q);
	  if (Q->maxdepth == 0)	/* estimate only */
	    return FALSE;	/* no nextbasis  */
	}

      if ( Q->truncate && negative(D->A[0][0]))   /* truncate when moving from opt. vertex */
          backtrack = TRUE;

      if (backtrack)		/* go back to prev. dictionary, restore i,j */
	{
	  backtrack = FALSE;

	  if (check_cache (D_p, Q, &i, &j))
	    {
	      if (Q->debug)
		fprintf (lrs_ofp, "\n Cached Dict. restored to depth %ld\n", D->depth);
	    }
	  else
	    {
	      D->depth--;
	      selectpivot (D, Q, &i, &j);
	      pivot (D, Q, i, j);
	      update (D, Q, &i, &j);	/*Update B,C,i,j */
	    }

	  if (Q->debug)
	    {
	      fprintf (lrs_ofp, "\n Backtrack Pivot: indices i=%ld j=%ld depth=%ld", i, j, D->depth);
	      printA (D, Q);
	    };

	  j++;			/* go to next column */
	}			/* end of if backtrack  */

      if (D->depth < Q->mindepth)
	break;

      /* try to go down tree */

      while ((j < d) && !reverse (D, Q, &i, j))
	j++;
      if (j == d)
	backtrack = TRUE;

      else
	/*reverse pivot found */
	{
	  cache_dict (D_p, Q, i, j);
	  /* Note that the next two lines must come _after_ the 
	     call to cache_dict */

	  D->depth++;
	  if (D->depth > Q->deepest)
	    Q->deepest++;

	  pivot (D, Q, i, j);
	  update (D, Q, &i, &j);	/*Update B,C,i,j */

	  D->lexflag = lexmin (D, Q, ZERO);	/* see if lexmin basis */
	  Q->count[2]++;
	  Q->totalnodes++;

	  save_basis (*D_p, Q);
	  if (Q->strace == Q->count[2])
	    Q->debug = TRUE;
	  if (Q->etrace == Q->count[2])
	    Q->debug = FALSE;
	  return TRUE;		/*return new dictionary */
	}


    }				/* end of main while loop for getnextbasis */
  return FALSE;			/* done, no more bases */
}				/*end of lrs_getnextbasis */

/*************************************/
/* print out one line of output file */
/*************************************/
long 
lrs_getvertex (lrs_dic * P, lrs_dat * Q, lrs_mp_vector output)
/*Print out current vertex if it is lexmin and return it in output */
/* return FALSE if no output generated  */
{
  long i;
  long ind;			/* output index                                  */
  long ired;			/* counts number of redundant columns            */
/* assign local variables to structures */
  long *redundcol = Q->redundcol;
  long *count = Q->count;
  long hull;
  long m, d;
  long lexflag;

  hull = Q->hull;
  m = P->m;
  d = P->d;
  lexflag = P->lexflag;

  if (lexflag || Q->allbases)
    ++(Q->count[1]);
  if (Q->debug)
    printA (P, Q);

  linint (Q->sumdet, 1, P->det, 1);
  if (Q->getvolume)
    updatevolume (P, Q);


  /*print cobasis if printcobasis=TRUE and count[2] a multiple of frequency */
  /* or for lexmin basis, except origin for hull computation - ugly!        */

  if (Q->printcobasis)
    if ((lexflag && !hull)  || ((Q->frequency > 0) && (count[2] == (count[2] / Q->frequency) * Q->frequency)))
      lrs_printcobasis (P, Q, ZERO);

  if (hull)
    return FALSE;		/* skip printing the origin */

  if (!lexflag && !Q->allbases && !Q->lponly)	/* not lexmin, and not printing forced */
    return FALSE;


  /* copy column 0 to output */

  i = 1;
  ired = 0;
  copy (output[0], P->det);

  for (ind = 1; ind < Q->n; ind++)	/* extract solution */

    if ((ired < Q->nredundcol) && (redundcol[ired] == ind))
      /* column was deleted as redundant */
      {
	itomp (ZERO, output[ind]);
	ired++;
      }
    else
      /* column not deleted as redundant */
      {
	getnextoutput (P, Q, i, ZERO, output[ind]);
	i++;
      }

  reducearray (output, Q->n);
  return TRUE;
}				/* end of lrs_getvertex */

long 
lrs_getray (lrs_dic * P, lrs_dat * Q, long col, long redcol, lrs_mp_vector output)
/*Print out solution in col and return it in output   */
/*redcol =n for ray/facet 0..n-1 for linearity column */
/*hull=1 implies facets will be recovered             */
/* return FALSE if no output generated in column col  */
{
  long i;
  long ind;			/* output index                                  */
  long ired;			/* counts number of redundant columns            */
/* assign local variables to structures */
  long *redundcol = Q->redundcol;
  long *count = Q->count;
  long hull = Q->hull;
  long n = Q->n;

  if (Q->debug)
    {
      printA (P, Q);
      for (i = 0; i < Q->nredundcol; i++)
	fprintf (lrs_ofp, " %ld", redundcol[i]);
      fflush(lrs_ofp);
    }

  if (redcol == n)
    {
      ++count[0];
      if (Q->printcobasis)
	lrs_printcobasis (P, Q, col);
    }

  i = 1;
  ired = 0;

  for (ind = 0; ind < n; ind++)	/* print solution */
    {
      if (ind == 0 && !hull)	/* must have a ray, set first column to zero */
	itomp (ZERO, output[0]);

      else if ((ired < Q->nredundcol) && (redundcol[ired] == ind))
	/* column was deleted as redundant */
	{
	  if (redcol == ind)	/* true for linearity on this cobasic index */
	    /* we print reduced determinant instead of zero */
	    copy (output[ind], P->det);
	  else
	    itomp (ZERO, output[ind]);
	  ired++;
	}
      else
	/* column not deleted as redundant */
	{
	  getnextoutput (P, Q, i, col, output[ind]);
	  i++;
	}

    }
  reducearray (output, n);
  return TRUE;
}				/* end of lrs_getray */


void 
getnextoutput (lrs_dic * P, lrs_dat * Q, long i, long col, lrs_mp out)
      /* get A[B[i][col] and copy to out */
{
  long row;
  long m = P->m;
  long d = P->d;
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *Row = P->Row;
  static long j;

  if (i == d && Q->voronoi)
    return;			/* skip last column if voronoi set */

  row = Row[i];

  if (Q->nonnegative)		/* if m+i basic get correct value from dictionary */
    {
      if (i == 1)		/* set up j */
	for (j = m - d + 1; j <= m && B[j] <= m; j++);

      if (j <= m && B[j] == m + i)
	copy (out, A[Row[j++]][col]);
      else
	{
         if ( i == col )
            copy(out,P->det);
         else
            itomp(ZERO,out);
	}
    }
  else
    copy (out, A[row][col]);

}				/* end of getnextoutput */

void 
lrs_printcobasis (lrs_dic * P, lrs_dat * Q, long col)
/* col is output column being printed */
{
  long i;
  long rflag;			/* used to find inequality number for ray column */
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Col = P->Col;
  long *Row = P->Row;
  long *inequality = Q->inequality;
  long *temparray = Q->temparray;
  long *count = Q->count;
  long hull = Q->hull;
  long d = P->d;
  long lastdv = Q->lastdv;
  long m=P->m;
  long firstime=TRUE;
  long nincidence;       /* count number of tight inequalities */

  if (hull)
    fprintf (lrs_ofp, "\nF#%ld B#%ld h=%ld vertices/rays ", count[0], count[2], P->depth);
  else if (Q->voronoi)
    fprintf (lrs_ofp, "\nV#%ld R#%ld B#%ld h=%ld data points ", count[1], count[0], count[2], P->depth);
  else
    fprintf (lrs_ofp, "\nV#%ld R#%ld B#%ld h=%ld facets ", count[1], count[0], count[2], P->depth);

  rflag = (-1);
  for (i = 0; i < d; i++)
    {
      temparray[i] = inequality[C[i] - lastdv];
      if (Col[i] == col)
	rflag = temparray[i];	/* look for ray index */
    }
  for (i = 0; i < d; i++)
    reorder (temparray, d);
  for (i = 0; i < d; i++)
    {
      fprintf (lrs_ofp, " %ld", temparray[i]);

      if (!(col == ZERO) && (rflag == temparray[i])) /* missing cobasis element for ray */
           fprintf (lrs_ofp, "*");

    }

/* get and print incidence information */
  if ( col == 0 )
     nincidence = d;
  else
     nincidence = d-1;

  for(i=lastdv+1;i<=m;i++)
    if ( zero (A[Row[i]][0] ))
       if( ( col == ZERO ) || zero (A[Row[i]] [col]) )
          { 
            nincidence++;
            if( Q->incidence )
              {
                if (firstime)
                  {
                    fprintf (lrs_ofp," :");
                    firstime = FALSE;
                   }
                fprintf(lrs_ofp," %ld",inequality[B[i] - lastdv ] );
              }
           }
         
  fprintf(lrs_ofp," I#%ld",nincidence);

  pmp (" det=", P->det);

}				/* end of lrs_printcobasis */

/*********************/
/* print final totals */
/*********************/
void 
lrs_printtotals (lrs_dic * P, lrs_dat * Q)
{
  long i;
  double x;
/* local assignments */
  double *cest = Q->cest;
  long *count = Q->count;
  long *inequality = Q->inequality;
  long *linearity = Q->linearity;
  long *temparray = Q->temparray;

  long *C = P->C;

  long hull = Q->hull;
  long homogeneous = Q->homogeneous;
  long nlinearity = Q->nlinearity;
  long nredundcol = Q->nredundcol;
  long m, d, lastdv;
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

  fprintf (lrs_ofp, "\nend");
  if (!Q->dualdeg)
    {
      fprintf (lrs_ofp, "\n*Warning: Starting dictionary is dual degenerate");
      fprintf (lrs_ofp, "\n*Complete enumeration may not have been produced\n");
    }
  if (Q->truncate)
    fprintf(lrs_ofp,"\n*Tree truncated at each new vertex");
  if (Q->maxdepth < MAXD)
    fprintf (lrs_ofp, "\n*Tree truncated at depth %ld", Q->maxdepth);


#ifdef LONG
  fprintf (lrs_ofp, "\n*Caution: no overflow checking with long integer arithemtic");
#else
  if( Q->verbose)
    {
      fprintf (lrs_ofp, "\n*Sum of det(B)=");
      pmp ("", Q->sumdet);
    }
#endif

/* next block with volume rescaling must come before estimates are printed */

  if (Q->getvolume)
    {
      rescalevolume (P, Q, Q->Nvolume, Q->Dvolume);

      if (Q->polytope)
	prat ("\n*Volume=", Q->Nvolume, Q->Dvolume);
      else
	prat ("\n*Pseudovolume=", Q->Nvolume, Q->Dvolume);
    }

  if (hull)     /* output things that are specific to hull computation */
    {
      fprintf (lrs_ofp, "\n*Totals: facets=%ld bases=%ld", count[0], count[2]);

      if (nredundcol > homogeneous)	/* don't count column 1 as redundant if homogeneous */
       {
        fprintf (lrs_ofp, " linearities=%ld", nredundcol - homogeneous);
        fprintf (lrs_ofp, " facets+linearities=%ld",nredundcol-homogeneous+count[0]);
       }


      if ((cest[2] > 0) || (cest[0] > 0))
      {
	fprintf (lrs_ofp, "\n*Estimates: facets=%g bases=%g", count[0] + cest[0], count[2] + cest[2]);
        if (Q->getvolume)
	  {
	    rattodouble (Q->Nvolume, Q->Dvolume, &x);
	    for (i = 2; i < d; i++)
	      cest[3] = cest[3] / i;	/*adjust for dimension */
	    fprintf (lrs_ofp, " volume=%g", cest[3] + x);
	  }
      fprintf (lrs_ofp, "\n*Total number of tree nodes evaluated: %ld", Q->totalnodes);

      }

      if ( Q-> restart || Q->allbases || (count[0] > 1 && !Q->homogeneous && !Q->polytope))
	    fprintf (lrs_ofp, "\n*Note! Duplicate facets may be present");

    }
  else         /* output things specific to vertex/ray computation */
    {
      fprintf (lrs_ofp, "\n*Totals: vertices=%ld rays=%ld bases=%ld", count[1], count[0], count[2]);

      if (nredundcol > 0)
        fprintf (lrs_ofp, " linearities=%ld", nredundcol);
      if ( count[0] + nredundcol > 0 )
         {
           fprintf (lrs_ofp, " vertices+rays");
           if ( nredundcol > 0 )
              fprintf (lrs_ofp, "+linearities");
           fprintf (lrs_ofp, "=%ld",nredundcol+count[0]+count[1]);
         }

      if ((cest[2] > 0) || (cest[0] > 0))
        {
	fprintf (lrs_ofp, "\n*Estimates: vertices=%g rays=%g bases=%g", count[1]+cest[1], count[0]+cest[0], count[2]+cest[2]);

         if (Q->getvolume)
	   {
	     rattodouble (Q->Nvolume, Q->Dvolume, &x);
	     for (i = 2; i <= d-homogeneous; i++)
	       cest[3] = cest[3] / i;	/*adjust for dimension */
	     fprintf (lrs_ofp, " pseudovolume=%g", cest[3] + x);
	   }
         fprintf (lrs_ofp, "\n*Total number of tree nodes evaluated: %ld", Q->totalnodes);
        }

      if (Q->restart || Q->allbases)        /* print warning  */
          fprintf (lrs_ofp, "\n*Note! Duplicate vertices/rays may be present");

      else if ( (count[0] > 1 && !Q->homogeneous))
          fprintf (lrs_ofp, "\n*Note! Duplicate rays may be present");

    }				/* end of output for vertices/rays */


  if(!Q->verbose)
     return;

  fprintf (lrs_ofp, "\n*Input size m=%ld rows n=%ld columns", P->m, Q->n);
  if (hull)
    fprintf (lrs_ofp, " working dimension=%ld", d - 1 + homogeneous);
  else
    fprintf (lrs_ofp, " working dimension=%ld", d);

  fprintf (lrs_ofp, "\n*Starting cobasis defined by input rows");
  for (i = 0; i < nlinearity; i++)
    temparray[i] = linearity[i];
  for (i = nlinearity; i < lastdv; i++)
    temparray[i] = inequality[C[i - nlinearity] - lastdv];
  for (i = 0; i < lastdv; i++)
    reorder (temparray, lastdv);
  for (i = 0; i < lastdv; i++)
    fprintf (lrs_ofp, " %ld", temparray[i]);

  fprintf (lrs_ofp, "\n*Dictionary Cache: max size= %ld misses= %ld/%ld   Tree Depth= %ld", dict_count, cache_misses, cache_tries, Q->deepest);


}				/* end of lrs_printtotals */
/************************/
/*  Estimation function */
/************************/
void 
lrs_estimate (lrs_dic * P, lrs_dat * Q)
		   /*get estimate of tree size from current node    */
		   /*current node is not counted.                   */
		   /*cest[0]rays [1]vertices [2]bases [3]volume     */
{

  lrs_mp Nvol, Dvol;		/* hold volume of current basis */
  long estdepth = 0;		/* depth of basis/vertex in subtree for estimate */
  long i = 0, j = 0, k, nchild, runcount, col;
  double prod = 0.0;
  double cave[] =
  {0.0, 0.0, 0.0, 0.0};
  double nvertices, nbases, nrays, nvol;
  long rays = 0;
  double newvol = 0.0;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *isave = Q->isave;
  long *jsave = Q->jsave;
  double *cest = Q->cest;
  long d = P->d;
  lrs_alloc_mp(Nvol); lrs_alloc_mp(Dvol);
/* Main Loop of Estimator */

  for (runcount = 1; runcount <= Q->runs; runcount++)
    {				/* runcount counts number of random probes */
      j = 0;
      nchild = 1;
      prod = 1;
      nvertices = 0.0;
      nbases = 0.0;
      nrays = 0.0;
      nvol = 0.0;
      while (nchild != 0)	/* while not finished yet */
	{

	  nchild = 0;
	  while (j < d)
	    {
	      if (reverse (P, Q, &i, j))
		{
		  isave[nchild] = i;
		  jsave[nchild] = j;
		  nchild++;
		}
	      j++;
	    }

	  if (estdepth == 0 && nchild == 0)
	    {
	      cest[0] = cest[0] + rays;		/* may be some rays here */
              lrs_clear_mp(Nvol); lrs_clear_mp(Dvol);
	      return;		/*subtree is a leaf */
	    }

	  prod = prod * nchild;
	  nbases = nbases + prod;
	  if (Q->debug)
	    {
	      fprintf (lrs_ofp, "   degree= %ld ", nchild);
	      fprintf (lrs_ofp, "\nPossible reverse pivots: i,j=");
	      for (k = 0; k < nchild; k++)
		fprintf (lrs_ofp, "%ld,%ld ", isave[k], jsave[k]);
	    }

	  if (nchild > 0)	/*reverse pivot found choose random child */
	    {
	      k = myrandom (Q->seed, nchild);
	      Q->seed = myrandom (Q->seed, 977L);
	      i = isave[k];
	      j = jsave[k];
	      if (Q->debug)
		fprintf (lrs_ofp, "  selected pivot k=%ld seed=%ld ", k, Q->seed);
	      estdepth++;
	      Q->totalnodes++;	/* calculate total number of nodes evaluated */
	      pivot (P, Q, i, j);
	      update (P, Q, &i, &j);	/*Update B,C,i,j */
	      if (lexmin (P, Q, ZERO))	/* see if lexmin basis for vertex */
		nvertices = nvertices + prod;

	      rays = 0;
	      for (col = 1; col <= d; col++)
		if (negative (A[0][col]) && (ratio (P, Q, col) == 0) && lexmin (P, Q, col))
		  rays++;
	      nrays = nrays + prod * rays;	/* update ray info */

	      if (Q->getvolume)
		{
		  rescaledet (P, Q, Nvol, Dvol);	/* scales determinant in case input rational */
		  rattodouble (Nvol, Dvol, &newvol);
		  nvol = nvol + newvol * prod;	/* adjusts volume for degree              */
		}
	      j = 0;
	    }
	}
      cave[0] = cave[0] + nrays;
      cave[1] = cave[1] + nvertices;
      cave[2] = cave[2] + nbases;
      cave[3] = cave[3] + nvol;

/*  backtrack to root and do it again */

      while (estdepth > 0)
	{
	  estdepth = estdepth - 1;
	  selectpivot (P, Q, &i, &j);
	  pivot (P, Q, i, j);
	  update (P, Q, &i, &j);	/*Update B,C,i,j */
	  /*fprintf(lrs_ofp,"\n0  +++"); */
	  if (Q->debug)
	    {
	      fprintf (lrs_ofp, "\n Backtrack Pivot: indices i,j %ld %ld ", i, j);
	      printA (P, Q);
	    }
	  j++;
	}

    }				/* end of for loop on runcount */
  for (i = 0; i < 4; i++)
    cest[i] = cave[i] / Q->runs + cest[i];
  lrs_clear_mp(Nvol); lrs_clear_mp(Dvol);
}				/* end of lrs_estimate  */


/*********************************/
/* Internal functions            */
/*********************************/
/* Basic Dictionary functions    */
/******************************* */

long 
reverse (lrs_dic * P, lrs_dat * Q, long *r, long s)
/*  find reverse indices  */
/* TRUE if B[*r] C[s] is a reverse lexicographic pivot */
{
  long i, j, enter, row, col;

/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long d = P->d;

  enter = C[s];
  col = Col[s];
  if (Q->debug)
    {
      fprintf (lrs_ofp, "\n+reverse: col index %ld C %ld Col %ld ", s, enter, col);
      fflush (stdout);
    }
  if (!negative (A[0][col]))
    {
      if (Q->debug)
	fprintf (lrs_ofp, " Pos/Zero Cost Coeff");
      return (FALSE);
    }

  *r = ratio (P, Q, col);
  if (*r == 0)			/* we have a ray */
    {
      if (Q->debug)
	fprintf (lrs_ofp, " Pivot col non-negative:  ray found");
      return (FALSE);
    }

  row = Row[*r];

/* check cost row after "pivot" for smaller leaving index    */
/* ie. j s.t.  A[0][j]*A[row][col] < A[0][col]*A[row][j]     */
/* note both A[row][col] and A[0][col] are negative          */

  for (i = 0; i < d && C[i] < B[*r]; i++)
    if (i != s)
      {
	j = Col[i];
	if (positive (A[0][j]) || negative (A[row][j]))		/*or else sign test fails trivially */
	  if ((!negative (A[0][j]) && !positive (A[row][j])) ||
	      comprod (A[0][j], A[row][col], A[0][col], A[row][j]) == -1)
	    {			/*+ve cost found */
	      if (Q->debug)
               {
		fprintf (lrs_ofp, "\nPositive cost found: index %ld C %ld Col %ld", i, C[i], j);
                fflush(lrs_ofp);
               }
	      return (FALSE);
	    }
      }
  if (Q->debug)
    {
      fprintf (lrs_ofp, "\n+end of reverse : indices r %ld s %ld \n", *r, s);
      fflush (stdout);
    }
  return (TRUE);
}				/* end of reverse */

long 
selectpivot (lrs_dic * P, lrs_dat * Q, long *r, long *s)
/* select pivot indices using lexicographic rule   */
/* returns TRUE if pivot found else FALSE          */
/* pivot variables are B[*r] C[*s] in locations Row[*r] Col[*s] */
{
  long j, col;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Col = P->Col;
  long d = P->d;

  *r = 0;
  *s = d;
  j = 0;

/*find positive cost coef */
  while ((j < d) && (!positive (A[0][Col[j]])))
    j++;

  if (j < d)			/* pivot column found! */
    {
      *s = j;
      col = Col[j];

      /*find min index ratio */
      *r = ratio (P, Q, col);
      if (*r != 0)
	return (TRUE);		/* unbounded */
    }
  return (FALSE);
}				/* end of selectpivot        */
/******************************************************* */

void 
pivot (lrs_dic * P, lrs_dat * Q, long bas, long cob)
		     /* Qpivot routine for array A              */
		     /* indices bas, cob are for Basis B and CoBasis C    */
		     /* corresponding to row Row[bas] and column       */
		     /* Col[cob]   respectively                       */
{
  long r, s;
  long i, j;
  lrs_mp Ns, Nt, Ars;
/* assign local variables to structures */

  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long m, d, m_A;

  lrs_alloc_mp(Ns); lrs_alloc_mp(Nt); lrs_alloc_mp(Ars);

  m = P->m;
  d = P->d;
  m_A = P->m_A;
  Q->count[3]++;    /* count the pivot */

  r = Row[bas];
  s = Col[cob];

/* Ars=A[r][s]    */
  if (Q->debug)
    {
      fprintf (lrs_ofp, "\n pivot  B[%ld]=%ld  C[%ld]=%ld ", bas, B[bas], cob, C[cob]);
      printA(P,Q);
      fflush (stdout);
    }
  copy (Ars, A[r][s]);
  storesign (P->det, sign (Ars));	/*adjust determinant to new sign */


  for (i = 0; i <= m_A; i++)
    if (i != r)
      for (j = 0; j <= d; j++)
	if (j != s)
	  {
/*          A[i][j]=(A[i][j]*Ars-A[i][s]*A[r][j])/P->det; */

	    mulint (A[i][j], Ars, Nt);
	    mulint (A[i][s], A[r][j], Ns);
	    decint (Nt, Ns);
	    exactdivint (Nt, P->det, A[i][j]);
	  }			/* end if j ....  */

  if (sign (Ars) == POS)
    {
      for (j = 0; j <= d; j++)	/* no need to change sign if Ars neg */
	/*   A[r][j]=-A[r][j];              */
	if (!zero (A[r][j]))
	  changesign (A[r][j]);
    }				/* watch out for above "if" when removing this "}" ! */
  else
    for (i = 0; i <= m_A; i++)
      if (!zero (A[i][s]))
	changesign (A[i][s]);

  /*  A[r][s]=P->det;                  */

  copy (A[r][s], P->det);		/* restore old determinant */
  copy (P->det, Ars);
  storesign (P->det, POS);		/* always keep positive determinant */


  if (Q->debug)
    {
      fprintf (lrs_ofp, " depth=%ld ", P->depth);
      pmp ("det=", P->det);
      fflush(stdout);
    }
/* set the new rescaled objective function value */

  mulint (P->det, Q->Lcm[0], Q->objden);
  mulint (Q->Gcd[0], A[0][0], Q->objnum);

  if (!Q->maximize)
        changesign (Q->objnum);
  if (zero (Q->objnum))
        storesign (Q->objnum, POS);

   reduce (Q->objnum,Q->objden);


  lrs_clear_mp(Ns); lrs_clear_mp(Nt); lrs_clear_mp(Ars);
}				/* end of pivot */

long 
primalfeasible (lrs_dic * P, lrs_dat * Q)
/* Do dual pivots to get primal feasibility */
/* Note that cost row is all zero, so no ratio test needed for Dual Bland's rule */
{
  long primalinfeasible = TRUE;
  long i, j;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Row = P->Row;
  long *Col = P->Col;
  long m, d, lastdv;
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

/*temporary: try to get new start after linearity */

  while (primalinfeasible)
    {
      i=lastdv+1;
      while (i <= m && !negative (A[Row[i]][0]) )
        i++;
      if (i <= m )
	{
	      j = 0;		/*find a positive entry for in row */
	      while (j < d && !positive (A[Row[i]][Col[j]]))
		j++;
	      if (j >= d)
		return (FALSE);	/* no positive entry */
	      pivot (P, Q, i, j);
	      update (P, Q, &i, &j);
        }
      else
         primalinfeasible = FALSE;
    }				/* end of while primalinfeasibile */
  return (TRUE);
}				/* end of primalfeasible */


long 
lrs_solvelp (lrs_dic * P, lrs_dat * Q, long maximize)
/* Solve primal feasible lp by least subscript and lex min basis method */
/* return TRUE if bounded, FALSE if unbounded                           */
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
  if (Q->debug)
    printA (P, Q);

  if (j < d && i == 0)		/* selectpivot gives information on unbounded solution */
    {
      if (Q->lponly)
	fprintf (lrs_ofp, "\n*Unbounded solution");
      return FALSE;
    }
   if (Q->verbose)
        prat ("\n*Objective function has value ", Q->objnum, Q->objden);
  return TRUE;
}				/* end of lrs_solvelp  */

long 
getabasis (lrs_dic * P, lrs_dat * Q, long order[])

/* Pivot Ax<=b to standard form */
/*Try to find a starting basis by pivoting in the variables x[1]..x[d]        */
/*If there are any input linearities, these appear first in order[]           */
/* Steps: (a) Try to pivot out basic variables using order                    */
/*            Stop if some linearity cannot be made to leave basis            */
/*        (b) Permanently remove the cobasic indices of linearities           */
/*        (c) If some decision variable cobasic, it is a linearity,           */
/*            and will be removed.                                            */
/*        (d) Dictionary is pivoted to get a feasible basis                   */

{
  long i, j, k;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long *linearity = Q->linearity;
  long *redundcol = Q->redundcol;
  long m, d, nlinearity;
  long nredundcol = 0L;		/* will be calculated here */
  m = P->m;
  d = P->d;
  nlinearity = Q->nlinearity;

  if (Q->debug)
    {
      fprintf (lrs_ofp, "\ngetabasis from inequalities given in order");
      for (i = 0; i < m; i++)
	fprintf (lrs_ofp, " %ld", order[i]);
    }
  for (j = 0; j < m; j++)
    {
      i = 0;
      while (i <= m && B[i] != d + order[j])
	i++;			/* find leaving basis index i */
      if (j < nlinearity && i > m)	/* cannot pivot linearity to cobasis */
	{
	  if (Q->debug)
	    printA (P, Q);
	  fprintf (lrs_ofp, "\nCannot find linearity in the basis");
	  return FALSE;
	}
      if (i <= m)
	{			/* try to do a pivot */
	  k = 0;
	  while (C[k] <= d && zero (A[Row[i]][Col[k]]))
	    k++;

	  if (C[k] <= d)
	    {
	      pivot (P, Q, i, k);
	      update (P, Q, &i, &k);
	    }
	  else if (j < nlinearity)
	    {			/* cannot pivot linearity to cobasis */
	      if (zero (A[Row[i]][0]))
		{
		  fprintf (lrs_ofp, "\n*Input linearity in row %ld is redundant--skipped", order[j]);
		  linearity[j] = 0;
		}
	      else
		{
		  if (Q->debug)
		    printA (P, Q);
		  fprintf (lrs_ofp, "\nInconsistent linearities");
		  return FALSE;
		}
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
  while (k < d && C[k] <= d)
    {
      if (C[k] <= d)		/* decision variable still in cobasis */
	redundcol[nredundcol++] = C[k] - Q->hull;	/* adjust for hull indices */
      k++;
    }

/* now we know how many decision variables remain in problem */
  Q->nredundcol = nredundcol;
  Q->lastdv = d - nredundcol;
  if (Q->debug)
    {
      fprintf (lrs_ofp, "\nend of first phase of getabasis: ");
      fprintf (lrs_ofp, "lastdv=%ld nredundcol=%ld", Q->lastdv, Q->nredundcol);
      fprintf (lrs_ofp, "\nredundant cobases:");
      for (i = 0; i < nredundcol; i++)
	fprintf (lrs_ofp, " %ld", redundcol[i]);
      printA (P, Q);
    }
/* Remove linearities from cobasis for rest of computation */
/* This is done in order so indexing is not screwed up */
  for (i = 0; i < nlinearity; i++)
    {				/* find cobasic index */
      k = 0;
      while (k < d && C[k] != linearity[i] + d)
	k++;
      if (k >= d)
	{
	  fprintf (lrs_ofp, "\nError removing linearity");
	  return FALSE;
	}
      if (!removecobasicindex (P, Q, k))
	return FALSE;
      d = P->d;
    }
  if (Q->debug && nlinearity > 0)
    printA (P, Q);
/* set index value for first slack variable */

/* Check feasability */
  if (Q->givenstart)
    {
      i = Q->lastdv + 1;
      while (i <= m && !negative (A[Row[i]][0]))
	i++;
      if (i <= m)
	fprintf (lrs_ofp, "\n*Infeasible startingcobasis - will be modified");
    }
  return TRUE;
}				/*  end of getabasis */

long 
removecobasicindex (lrs_dic * P, lrs_dat * Q, long k)
/* remove the variable C[k] from the problem */
/* used after detecting column dependency    */
{
  long i, j, cindex, deloc;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Col = P->Col;
  long m, d;
  m = P->m;
  d = P->d;

  if (Q->debug)
    fprintf (lrs_ofp, "\nremoving cobasic index k=%ld C[k]=%ld", k, C[k]);
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

  if (deloc != d)               
    {
  /* copy col d to deloc */
      for (i = 0; i <= m; i++)
        copy (A[i][deloc], A[i][d]);

  /* reassign location for moved column */
      j = 0;
      while (Col[j] != d)
        j++;
      
    Col[j] = deloc;
  }

  P->d--;
  if (Q->debug)
    printA (P, Q);
  return TRUE;
}				/* end of removecobasicindex */

lrs_dic *
resize (lrs_dic * P, lrs_dat * Q)
	/* resize the dictionary after some columns are deleted, ie. inputd>d */
	/* a new lrs_dic record is created with reduced size, and items copied over */
{
  lrs_dic *P1;			/* to hold new dictionary in case of resizing */

  long i, j;
  long m, d, m_A;


  m = P->m;
  d = P->d;
  m_A = P->m_A;

/* get new dictionary record */

  P1 = new_lrs_dic (m, d, m_A);

/* copy data from P to P1    */
  P1->i = P->i;
  P1->j = P->j;
  P1->depth = P->depth;
  P1->m = P->m;
  P1->d = P->d;
  P1->lexflag = P->lexflag;
  P1->m_A = P->m_A;
  copy (P1->det, P->det);

  for (i = 0; i <= m; i++)
    {
      P1->B[i] = P->B[i];
      P1->Row[i] = P->Row[i];
    }
  for (i = 0; i <= m_A; i++)
    {
      for (j = 0; j <= d; j++)
	copy (P1->A[i][j], P->A[i][j]);
    }


  for (j = 0; j <= d; j++)
    {
      P1->Col[j] = P->Col[j];
      P1->C[j] = P->C[j];
    }

  if (Q->debug)
    {
      fprintf (lrs_ofp, "\nDictionary resized from d=%ld to d=%ld", Q->inputd, P->d);
      printA (P1, Q);
    }

  lrs_free_dic (P,Q);

/* Reassign cache pointers */

  Q->Qhead = P1;
  Q->Qtail = P1;
  P1->next = P1;
  P1->prev = P1;

  return P1;

}
/********* resize                    ***************/


long 
restartpivots (lrs_dic * P, lrs_dat * Q)
/* facet contains a list of the inequalities in the cobasis for the restart */
/* inequality contains the relabelled inequalities after initialization     */
{
  long i, j, k;
  long *Cobasic;		/* when restarting, Cobasic[j]=1 if j is in cobasis */
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
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

  Cobasic = (long *) CALLOC ((unsigned) m + d + 2, sizeof (long));

  if (Q->debug)
    fprintf(lrs_ofp,"\nCobasic flags in restartpivots");
  /* set Cobasic flags */
  for (i = 0; i < m + d + 1; i++)
    Cobasic[i] = 0;
  for (i = 0; i < d; i++)	/* find index corresponding to facet[i] */
    {
      j = 1;
      while (facet[i + nlinearity] != inequality[j])
	j++;
      Cobasic[j + lastdv] = 1;
      if (Q->debug)
        fprintf(lrs_ofp," %ld %ld;",facet[i+nlinearity],j+lastdv);
    }

  /* Note that the order of doing the pivots is important, as */
  /* the B and C vectors are reordered after each pivot       */

  for (i = m; i >= d + 1; i--)	/*see if a basic variable should leave */
    if (Cobasic[B[i]])		/* basic variable must leave  */
      {
	/* find a Cobasic variable that must enter basis */
	k = d - 1;
	while ((k >= 0) &&
	       (zero (A[Row[i]][Col[k]]) || Cobasic[C[k]]))
	  k--;
	if (k >= 0)
	  {
	    pivot (P, Q, i, k);
	    update (P, Q, &i, &k);
	  }
	else
	  {
	    fprintf (lrs_ofp, "\nInvalid Co-basis - does not have correct rank");
            free(Cobasic);
	    return FALSE;
	  }
      }				/* end of if(Cobasic..  */
  if (lexmin (P, Q, ZERO))
    --Q->count[1];		/* decrement vertex count if lexmin */
/* check restarting from a primal feasible dictionary               */
  for (i = lastdv + 1; i <= m; i++)
    if (negative (A[Row[i]][0]))
      {
	fprintf (lrs_ofp, "\nTrying to restart from infeasible dictionary");
        free(Cobasic);
	return FALSE;
      }
  free(Cobasic);
  return TRUE;

}				/* end of restartpivots */


long 
ratio (lrs_dic * P, lrs_dat * Q, long col)	/*find lex min. ratio */
		  /* find min index ratio -aig/ais, ais<0 */
		  /* if multiple, checks successive basis columns */
		  /* recoded Dec 1997                     */
{
  long i, j, comp, ratiocol, basicindex, start, nstart, cindex, bindex;
  long firstime;		/*For ratio test, true on first pass,else false */
  lrs_mp Nmin, Dmin;
  long degencount, ndegencount;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
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
      if (negative (A[Row[j]][col]))
	minratio[degencount++] = j;
    }				/* end of for loop */
  if (Q->debug)
    {
      fprintf (lrs_ofp, "  Min ratios: ");
      for (i = 0; i < degencount; i++)
	fprintf (lrs_ofp, " %ld ", B[minratio[i]]);
    }
  if (degencount == 0)
    return (degencount);	/* non-negative pivot column */

  lrs_alloc_mp(Nmin); lrs_alloc_mp(Dmin);
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
	  firstime = TRUE;
	  /*get next ratio column and increment cindex */
	  if (basicindex != d)
	    ratiocol = Col[cindex++];
	  for (j = start; j < start + degencount; j++)
	    {
	      i = Row[minratio[j]];	/* i is the row location of the next basic variable */
	      comp = 1;		/* 1:  lhs>rhs;  0:lhs=rhs; -1: lhs<rhs */
	      if (firstime)
		firstime = FALSE;	/*force new min ratio on first time */
	      else
		{
		  if (positive (Nmin) || negative (A[i][ratiocol]))
		    {
		      if (negative (Nmin) || positive (A[i][ratiocol]))
			comp = comprod (Nmin, A[i][col], A[i][ratiocol], Dmin);
		      else
			comp = -1;
		    }

		  else if (zero (Nmin) && zero (A[i][ratiocol]))
		    comp = 0;

		  if (ratiocol == ZERO)
		    comp = -comp;	/* all signs reversed for rhs */
		}
	      if (comp == 1)
		{		/*new minimum ratio */
		  nstart = j;
		  copy (Nmin, A[i][ratiocol]);
		  copy (Dmin, A[i][col]);
		  ndegencount = 1;
		}
	      else if (comp == 0)	/* repeated minimum */
		minratio[nstart + ndegencount++] = minratio[j];

	    }			/* end of  for (j=start.... */
	  degencount = ndegencount;
	  start = nstart;
	}			/* end of else perform ratio test statement */
      basicindex++;		/* increment column of basis inverse to check next */
      if (Q->debug)
	{
	  fprintf (lrs_ofp, " ratiocol=%ld degencount=%ld ", ratiocol, degencount);
	  fprintf (lrs_ofp, "  Min ratios: ");
	  for (i = start; i < start + degencount; i++)
	    fprintf (lrs_ofp, " %ld ", B[minratio[i]]);
	}
    }				/*end of while loop */
  lrs_clear_mp(Nmin); lrs_clear_mp(Dmin);
  return (minratio[start]);
}				/* end of ratio */



long 
lexmin (lrs_dic * P, lrs_dat * Q, long col)
  /*test if basis is lex-min for vertex or ray, if so TRUE */
  /* FALSE if a_r,g=0, a_rs !=0, r > s          */
{
/*do lexmin test for vertex if col=0, otherwise for ray */
  long r, s, i, j;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long m = P->m;
  long d = P->d;

  for (i = Q->lastdv + 1; i <= m; i++)
    {
      r = Row[i];
      if (zero (A[r][col]))	/* necessary for lexmin to fail */
	for (j = 0; j < d; j++)
	  {
	    s = Col[j];
	    if (B[i] > C[j])	/* possible pivot to reduce basis */
	      {
		if (zero (A[r][0]))	/* no need for ratio test, any pivot feasible */
		  {
		    if (!zero (A[r][s]))
		      return (FALSE);
		  }
		else if (negative (A[r][s]) && ismin (P, Q, r, s))
		  {
		    return (FALSE);
		  }
	      }			/* end of if B[i] ... */
	  }
    }
  if ((col != ZERO) && Q->debug)
    {
      fprintf (lrs_ofp, "\n lexmin ray in col=%ld ", col);
      printA (P, Q);
    }
  return (TRUE);
}				/* end of lexmin */

long 
ismin (lrs_dic * P, lrs_dat * Q, long r, long s)
/*test if A[r][s] is a min ratio for col s */
{
  long i;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long m_A = P->m_A;

  for (i = 1; i <= m_A; i++)
    if ((i != r) &&
	negative (A[i][s]) && comprod (A[i][0], A[r][s], A[i][s], A[r][0]))
      {
	return (FALSE);
      }

  return (TRUE);
}

void 
update (lrs_dic * P, lrs_dat * Q, long *i, long *j)
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

long 
lrs_degenerate (lrs_dic * P, lrs_dat * Q)
/* TRUE if the current dictionary is primal degenerate */
/* not thoroughly tested   2000/02/15                  */
{
  long i;
  long *B, *Row;

  lrs_mp_matrix A = P->A;
  long d = P->d;
  long m = P->m;

  B = P->B;
  Row = P->Row;

  for (i = d + 1; i <= m; i++)
    if (zero (A[Row[i]][0]))
      return TRUE;

  return FALSE;
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

void 
reorder1 (long a[], long b[], long newone, long range)
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
}				/* end of reorder1 */


void 
rescaledet (lrs_dic * P, lrs_dat * Q, lrs_mp Vnum, lrs_mp Vden)
  /* rescale determinant to get its volume */
  /* Vnum/Vden is volume of current basis  */
{
  lrs_mp gcdprod;		/* to hold scale factors */
  long i;
/* assign local variables to structures */
  long *B = P->B;
  long *C = P->C;
  long m, d, lastdv;

  lrs_alloc_mp(gcdprod);
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

  itomp (ONE, gcdprod);
  itomp (ONE, Vden);
  for (i = 0; i < d; i++)
    if (B[i] <= m)
      {
	mulint (Q->Gcd[Q->inequality[C[i] - lastdv]], gcdprod, gcdprod);
	mulint (Q->Lcm[Q->inequality[C[i] - lastdv]], Vden, Vden);
      }
  mulint (P->det, gcdprod, Vnum);
  reduce (Vnum, Vden);
  lrs_clear_mp(gcdprod);
}				/* end rescaledet */

void 
rescalevolume (lrs_dic * P, lrs_dat * Q, lrs_mp Vnum, lrs_mp Vden)
/* adjust volume for dimension */
{
  lrs_mp temp, dfactorial;
/* assign local variables to structures */
  long lastdv = Q->lastdv;

  lrs_alloc_mp(temp); lrs_alloc_mp(dfactorial);

  /*reduce Vnum by d factorial  */
  getfactorial (dfactorial, lastdv);
  mulint (dfactorial, Vden, Vden);
  if (Q->hull && !Q->homogeneous)
    {				/* For hull option multiply by d to correct for lifting */
      itomp (lastdv, temp);
      mulint (temp, Vnum, Vnum);
    }

  reduce (Vnum, Vden);
  lrs_clear_mp(temp); lrs_clear_mp(dfactorial);
}


void 
updatevolume (lrs_dic * P, lrs_dat * Q)		/* rescale determinant and update the volume */
{
  lrs_mp tN, tD, Vnum, Vden;
  lrs_alloc_mp(tN); lrs_alloc_mp(tD); lrs_alloc_mp(Vnum); lrs_alloc_mp(Vden);
  rescaledet (P, Q, Vnum, Vden);
  copy (tN, Q->Nvolume);
  copy (tD, Q->Dvolume);
  linrat (tN, tD, ONE, Vnum, Vden, ONE, Q->Nvolume, Q->Dvolume);
  if (Q->debug)
    {
      prat ("\n*Volume=", Q->Nvolume, Q->Dvolume);
      pmp (" Vnum=", Vnum);
      pmp (" Vden=", Vden);
    }
  lrs_clear_mp(tN); lrs_clear_mp(tD); lrs_clear_mp(Vnum); lrs_clear_mp(Vden);

}				/* end of updatevolume */



/***************************************************/
/* Routines for redundancy checking                */
/***************************************************/

long 
checkredund (lrs_dic * P, lrs_dat * Q)
/* Solve primal feasible lp by least subscript and lex min basis method */
/* to check redundancy of a row in objective function                   */
/* returns TRUE if redundant, else FALSE                                */
{
  lrs_mp Ns, Nt;
  long i, j;
  long r, s;

/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B, *C, *Row, *Col;
  long d = P->d;

  lrs_alloc_mp(Ns); lrs_alloc_mp(Nt);
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

      mulint (A[0][s], A[r][0], Ns);
      mulint (A[0][0], A[r][s], Nt);

      if (greater (Ns, Nt))
        {
          lrs_clear_mp(Ns); lrs_clear_mp(Nt);
	  return FALSE;		/* non-redundant */
        }

      pivot (P, Q, i, j);
      update (P, Q, &i, &j);	/*Update B,C,i,j */

    }

  lrs_clear_mp(Ns); lrs_clear_mp(Nt);

  return !(j < d && i == 0);	/* unbounded is also non-redundant */


}				/* end of checkredund  */

long 
checkcobasic (lrs_dic * P, lrs_dat * Q, long index)
/* TRUE if index is cobasic and nonredundant                         */
/* FALSE if basic, or degen. cobasic, where it will get pivoted out  */

{

/* assign local variables to structures */

  lrs_mp_matrix A = P->A;
  long *B, *C, *Row, *Col;
  long d = P->d;
  long m = P->m;
  long debug = Q->debug;
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
    return FALSE;		/* not cobasic index */


/* index is cobasic */

  if (debug)
    fprintf (lrs_ofp, "\nindex=%ld cobasic", index);
/* not debugged for new LOC 
   s=LOC[index];
 */
  s = Col[j];
  i = Q->lastdv + 1;

  while ((i <= m) &&
	 (zero (A[Row[i]][s]) || !zero (A[Row[i]][0])))
    i++;

  if (i > m)
    {
      if (debug)
	fprintf (lrs_ofp, " is non-redundant");
      return TRUE;
    }
  if (debug)
    fprintf (lrs_ofp, " is degenerate B[i]=%ld", B[i]);

  pivot (P, Q, i, j);
  update (P, Q, &i, &j);	/*Update B,C,i,j */

  return FALSE;			/*index is no longer cobasic */

}				/* end of checkcobasic */

long 
checkindex (lrs_dic * P, lrs_dat * Q, long index)
/* 0 if index is non-redundant inequality */
/* 1 if index is redundant     inequality */
/* 2 if index is input linearity          */
/*NOTE: row is returned all zero if redundant!! */
{
  long i, j;

  lrs_mp_matrix A = P->A;
  long *Row = P->Row;
  long *B = P->B;
  long d = P->d;
  long m = P->m;

  if (Q->debug)
    printA (P, Q);

/* each slack index must be checked for redundancy */
/* if in cobasis, it is pivoted out if degenerate */
/* else it is non-redundant                       */

  if (checkcobasic (P, Q, index))
    return ZERO;

/* index is basic   */
/* not debugged for new LOC 
   i=LOC[index];
 */
  j = 1;
  while ((j <= m) && (B[j] != index))
    j++;

  i = Row[j];

  /* copy row i to cost row, and set it to zero */

  for (j = 0; j <= d; j++)
    {
      copy (A[0][j], A[i][j]);
      changesign (A[0][j]);
      itomp (ZERO, A[i][j]);
    }


  if (checkredund (P, Q))
    return ONE;

/* non-redundant, copy back and change sign */

  for (j = 0; j <= d; j++)
    {
      copy (A[i][j], A[0][j]);
      changesign (A[i][j]);
    }

  return ZERO;

}				/* end of checkindex */

/***************************************************************/
/*                                                             */
/*            Package of I/O routines                          */
/*                                                             */
/***************************************************************/

void
lprat (const char *name, long Nt, long Dt)
/*print the long precision rational Nt/Dt without reducing  */
{
  if ( Nt > 0 ) 
    fprintf (lrs_ofp, " ");
  fprintf (lrs_ofp, "%s%ld", name, Nt);
  if (Dt != 1)   
    fprintf (lrs_ofp, "/%ld", Dt);
  fprintf (lrs_ofp, " ");
}                               /* lprat */

long
lreadrat (long *Num, long *Den)
 /* read a rational string and convert to long    */
 /* returns true if denominator is not one        */
{
  char in[MAXINPUT], num[MAXINPUT], den[MAXINPUT];
  fscanf (lrs_ifp, "%s", in);
  atoaa (in, num, den);         /*convert rational to num/dem strings */
  *Num = atol (num);
  if (den[0] == '\0')
    {
      *Den = 1L;
      return (FALSE);
    }
  *Den = atol (den);
  return (TRUE);
}

void
lrs_getinput(lrs_dic *P,lrs_dat *Q,long *num,long *den, long m, long d)
/* code for reading data matrix in lrs/cdd format */
{
  long j,row;

  printf("\nEnter each row: b_i  a_ij j=1..%ld",d);
  for (row=1;row<=m;row++)
    {
      printf("\nEnter row %ld: ",row );
      for(j=0;j<=d;j++)
       {
         lreadrat(&num[j],&den[j]);
         lprat(" ",num[j],den[j]);
       }

       lrs_set_row(P,Q,row,num,den,GE);
     }

 printf("\nEnter objective row c_j j=1..%ld: ",d);
 num[0]=0; den[0]=1;
 for(j=1;j<=d;j++)
       {
         lreadrat(&num[j],&den[j]);
         lprat(" ",num[j],den[j]);
       }

 lrs_set_obj(P,Q,num,den,MAXIMIZE);
}


long 
readlinearity (lrs_dat * Q)	/* read in and check linearity list */
{
  long i, j;
  long nlinearity;
  fscanf (lrs_ifp, "%ld", &nlinearity);
  if (nlinearity < 1)
    {
      fprintf (lrs_ofp, "\nLinearity option invalid, indices must be positive");
      return (FALSE);
    } 

  Q->linearity  = CALLOC ((nlinearity + 1), sizeof (long));

  for (i = 0; i < nlinearity; i++)
    {
      fscanf (lrs_ifp, "%ld", &j);
      Q->linearity[i] = j;	

    }
  for (i = 1; i < nlinearity; i++)	/*sort in order */
    reorder (Q->linearity, nlinearity);

  Q->nlinearity = nlinearity;
  Q->polytope = FALSE;
  return TRUE;
}				/* end readlinearity */


long 
readfacets (lrs_dat * Q, long facet[])
/* read and check facet list for obvious errors during start/restart */
/* this must be done after linearity option is processed!! */
{
  long i, j;
/* assign local variables to structures */
  long m, d;
  long *linearity = Q->linearity;
  m = Q->m;
  d = Q->inputd;

  for (j = Q->nlinearity; j < d; j++)	/* note we place these after the linearity indices */
    {
      fscanf (lrs_ifp, "%ld", &facet[j]);
      fprintf (lrs_ofp, " %ld", facet[j]);
      if (facet[j] < 1 || facet[j] > m)
	{
	  fprintf (lrs_ofp, "\n Start/Restart cobasic indices must be in range 1 .. %ld ", m);
	  return FALSE;
	}
      for (i = 0; i < Q->nlinearity; i++)
	if (linearity[i] == facet[j])
	  {
	    fprintf (lrs_ofp, "\n Start/Restart cobasic indices should not include linearities");
	    return FALSE;
	  }
      for (i = 0; i < j; i++)
	if (facet[i] == facet[j])
	  {
	    fprintf (lrs_ofp, "\n  Start/Restart cobasic indices must be distinct");
	    return FALSE;
	  }
    }
  return TRUE;
}				/* end of readfacets */

void 
printA (lrs_dic * P, lrs_dat * Q)	/* print the integer m by n array A 
					   with B,C,Row,Col vectors         */
{
  long i, j;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long m, d, lastdv;
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

  fprintf (lrs_ofp, "\n Basis    ");
  for (i = 0; i <= m; i++)
    fprintf (lrs_ofp, "%ld ", B[i]);
  fprintf (lrs_ofp, " Row ");
  for (i = 0; i <= m; i++)
    fprintf (lrs_ofp, "%ld ", Row[i]);
  fprintf (lrs_ofp, "\n Co-Basis ");
  for (i = 0; i <= d; i++)
    fprintf (lrs_ofp, "%ld ", C[i]);
  fprintf (lrs_ofp, " Column ");
  for (i = 0; i <= d; i++)
    fprintf (lrs_ofp, "%ld ", Col[i]);
  pmp (" det=", P->det);
  fprintf (lrs_ofp, "\n");
  i=0;
  while ( i<= m )
    {
      for (j = 0; j <= d; j++)
	pimat (P, i, j, A[Row[i]][Col[j]], "A");
      fprintf (lrs_ofp, "\n");
      if (i==0 && Q->nonnegative)  /* skip basic rows - don't exist! */
          i=d;
      i++;
  fflush (stdout);
    }
  fflush (stdout);
}


void 
pimat (lrs_dic * P, long r, long s, lrs_mp Nt, char name[])
 /*print the long precision integer in row r col s of matrix A */
{
  long *B = P->B;
  long *C = P->C;
  if (s == 0)
    fprintf (lrs_ofp, "%s[%ld][%ld]=", name, B[r], C[s]);
  else
    fprintf (lrs_ofp, "[%ld]=", C[s]);
  pmp ("", Nt);

}

/***************************************************************/
/*                                                             */
/*     Routines for caching, allocating etc.                   */
/*                                                             */
/***************************************************************/

/* From here mostly Bremner's handiwork */

void 
cache_dict (lrs_dic ** D_p, lrs_dat * global, long i, long j)
{

  if (dict_limit > 1)
    {
      /* save row, column indicies */
      (*D_p)->i = i;
      (*D_p)->j = j;

/* Make a new, blank spot at the end of the queue to copy into */ 

      pushQ (global, (*D_p)->m, (*D_p)->d, (*D_p)->m_A);

      copy_dict (global, global->Qtail, *D_p);	/* Copy current dictionary */
    }
  *D_p = global->Qtail;

}

void 
copy_dict (lrs_dat * global, lrs_dic * dest, lrs_dic * src)
{
  long m = src->m;
  long m_A = src->m_A;        /* number of rows in A */
  long d = src->d;

#ifdef GMP 
  long r,s;
  for ( r=0;r<=m_A;r++)
    for( s=0;s<=d;s++)
       copy(dest->A[r][s],src->A[r][s]);

#else            /* fast copy for MP and LONG arithmetic */
  /* Note that the "A" pointer trees need not be copied, since they
     always point to the same places within the corresponding space
*/
  memcpy (dest->A[0][0], (global->Qtail->prev)->A[0][0],
          (d + 1) * (lrs_digits + 1) * (m_A + 1) * sizeof (long));

#endif

  dest->i = src->i;
  dest->j = src->j;
  dest->m = m;
  dest->d = d;
  dest->m_A  = src->m_A;

  dest->depth = src->depth;
  dest->lexflag = src->lexflag;

  copy (dest->det, src->det);

  if (global->debug)
    fprintf (lrs_ofp, "\nSaving dict at depth %ld\n", src->depth);

  memcpy (dest->B, src->B, (m + 1) * sizeof (long));
  memcpy (dest->C, src->C, (d + 1) * sizeof (long));
  memcpy (dest->Row, src->Row, (m + 1) * sizeof (long));
  memcpy (dest->Col, src->Col, (d + 1) * sizeof (long));
}

/* 
 * pushQ(lrs_dat *globals,m,d):
 * this routine ensures that Qtail points to a record that 
 * may be copied into.
 *
 * It may create a new record, or it may just move the head pointer
 * forward so that know that the old record has been overwritten.
 */
#if 0
#define TRACE(s) fprintf(stderr,"\n%s %p %p\n",s,global->Qhead,global->Qtail);
#else
#define TRACE(s)
#endif

void 
pushQ (lrs_dat * global, long m, long d ,long m_A)
{

  if ((global->Qtail->next) == global->Qhead)
    {
      /* the Queue is full */
      if (dict_count < dict_limit)
	{
	  /* but we are allowed to create more */
	  lrs_dic *p;

	  p = new_lrs_dic (m, d, m_A);

	  if (p)
	    {

	      /* we successfully created another record */

	      p->next = global->Qtail->next;
	      (global->Qtail->next)->prev = p;
	      (global->Qtail->next) = p;
	      p->prev = global->Qtail;

	      dict_count++;
	      global->Qtail = p;

	      TRACE ("Added new record to Q");

	    }
	  else
	    {
	      /* virtual memory exhausted. bummer */
	      global->Qhead = global->Qhead->next;
	      global->Qtail = global->Qtail->next;

	      TRACE ("VM exhausted");
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
	  TRACE ("User  limit");

	}
    }

  else
    {
      global->Qtail = global->Qtail->next;
      TRACE ("Reusing");
    }
}

#define NULLRETURN(e) if (!(e)) return NULL;

lrs_dic *
new_lrs_dic (long m, long d, long m_A)
{
  lrs_dic *p;

  NULLRETURN (p = (lrs_dic *) malloc (sizeof (lrs_dic)));


  NULLRETURN (p->B = calloc ((m + 1), sizeof (long)));
  NULLRETURN (p->Row = calloc ((m + 1), sizeof (long)));

  NULLRETURN (p->C = calloc ((d + 1), sizeof (long)));
  NULLRETURN (p->Col = calloc ((d + 1), sizeof (long)));

#ifdef GMP
  lrs_alloc_mp(p->det);
#endif

  p->A=lrs_alloc_mp_matrix(m_A,d);

  return p;
}

void 
lrs_free_dic (lrs_dic * P, lrs_dat *Q)
{
/* do the same steps as for allocation, but backwards */
/* gmp variables cannot be cleared using free: use lrs_clear_mp* */
  lrs_dic *P1;
  long d = P->d;
  long m_A = P->m_A;

/* repeat until cache is empty */

  do
  {
  lrs_clear_mp_matrix (P->A,m_A,d);

/* "it is a ghastly error to free something not assigned my malloc" KR167 */
/* so don't try: free (P->det);                                           */

  lrs_clear_mp (P->det);      

  free (P->Row);
  free (P->Col);
  free (P->C);
  free (P->B);

/* go to next record in cache if any */
  P1 =P->next;
  free (P);
  P=P1;

  }  while (Q->Qhead != P );


}

void
lrs_free_dat ( lrs_dat *Q )
{
   long m=Q->m;

/* most of these items were allocated in lrs_alloc_dic */

  lrs_clear_mp_vector (Q->Gcd,m);
  lrs_clear_mp_vector (Q->Lcm,m);

  lrs_clear_mp (Q->sumdet);
  lrs_clear_mp (Q->Nvolume);
  lrs_clear_mp (Q->Dvolume);
  lrs_clear_mp (Q->objnum);
  lrs_clear_mp (Q->objden);
  lrs_clear_mp (Q->saved_det);

  free (Q->inequality);
  free (Q->facet);
  free (Q->redundcol);
  free (Q->linearity);
  free (Q->minratio);
  free (Q->temparray);

  free (Q->name);  
  free (Q->saved_C);

  lrs_global_count--;

  free(Q);
}


long 
check_cache (lrs_dic ** D_p, lrs_dat * global, long *i_p, long *j_p)
{
/* assign local variables to structures */



  cache_tries++;

  if (global->Qtail == global->Qhead)
    {
      TRACE ("cache miss");
      /* Q has only one element */
      cache_misses++;
      return 0;

    }
  else
    {
      global->Qtail = global->Qtail->prev;

      *D_p = global->Qtail;

      *i_p = global->Qtail->i;
      *j_p = global->Qtail->j;

      TRACE ("restoring dict");
      return 1;
    }
}


lrs_dic *
lrs_alloc_dic (lrs_dat * Q)
/* allocate and initialize lrs_dic */
{

  lrs_dic *p;
  long i, j;
  long m, d, m_A;

  if (Q->hull)                       /* d=col dimension of A */
    Q->inputd = Q->n;                /* extra column for hull */
  else
    Q->inputd = Q->n - 1;

  m = Q->m;
  d = Q->inputd;
  m_A = m;   /* number of rows in A */

/* nonnegative flag set means that problem is d rows "bigger"     */
/* since nonnegative constraints are not kept explicitly          */

  if(Q->nonnegative)
    m = m+d;

  p = new_lrs_dic (m, d, m_A);
  if (!p)
    return NULL;

  p->next = p;
  p->prev = p;
  Q->Qhead = p;
  Q->Qtail = p;


  dict_count = 1;
  dict_limit = 10;
  cache_tries = 0;
  cache_misses = 0;

/* Initializations */

  p->d = d;
  p->m = m;
  p->m_A  = m_A;
  p->depth = 0L;
  p->lexflag = TRUE;
  itomp (ONE, p->det);

/*m+d+1 is the number of variables, labelled 0,1,2,...,m+d  */
/*  initialize array to zero   */
  for (i = 0; i <= m_A; i++)
    for (j = 0; j <= d; j++)
      itomp (ZERO, p->A[i][j]);

  Q->inequality = CALLOC ((m + 1), sizeof (long));
  if (Q->nlinearity == ZERO)   /* linearity may already be allocated */
      Q->linearity  = CALLOC ((m + 1), sizeof (long));

  Q->facet = (long *) CALLOC ((unsigned) d + 1, sizeof (long));
  Q->redundcol = CALLOC ((d + 1), sizeof (long));
  Q->minratio = CALLOC ((m + 1), sizeof (long));
  Q->temparray = CALLOC ((unsigned) d + 1, sizeof (long));

  Q->inequality[0] = 2L;
  Q->Gcd = lrs_alloc_mp_vector(m);
  Q->Lcm = lrs_alloc_mp_vector(m);
  Q->saved_C = CALLOC (d + 1, sizeof (long));

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
}				/* end of lrs_alloc_dic */


/* 
   this routine makes a copy of the information needed to restart, 
   so that we can guarantee that if a signal is received, we 
   can guarantee that nobody is messing with it.
   This as opposed to adding all kinds of critical regions in 
   the main line code.

   It is also used to make sure that in case of overflow, we
   have a valid cobasis to restart from.
 */
void 
save_basis (lrs_dic * P, lrs_dat * Q)
{
  int i;
/* assign local variables to structures */
  long *C = P->C;
  long d;

#ifdef SIGNALS
  sigset_t oset, blockset;
  sigemptyset (&blockset);
  sigaddset (&blockset, SIGTERM);
  sigaddset (&blockset, SIGHUP);
  sigaddset (&blockset, SIGUSR1);

  errcheck ("sigprocmask", sigprocmask (SIG_BLOCK, &blockset, &oset));
#endif
  d = P->d;

  Q->saved_flag = 1;

  for (i = 0; i < 3; i++)
    Q->saved_count[i] = Q->count[i];

  for (i = 0; i < d + 1; i++)
    Q->saved_C[i] = C[i];

  copy (Q->saved_det, P->det);

  Q->saved_d = P->d;
  Q->saved_depth = P->depth;

#ifdef SIGNALS
  errcheck ("sigprocmask", sigprocmask (SIG_SETMASK, &oset, 0));
#endif
}

/* digits overflow is a call from lrs_mp package */

void 
digits_overflow ()
{
  fprintf (lrs_ofp, "\nOverflow at digits=%ld", DIG2DEC (lrs_digits));
  fprintf (lrs_ofp, "\nRerun with option: digits n, where n > %ld\n", DIG2DEC (lrs_digits));
  lrs_dump_state ();

  notimpl("");
}

static void 
lrs_dump_state ()
{
  long i;

  fprintf (stderr, "\n\nlrs_lib: checkpointing:\n");

  fprintf (stderr, "lrs_lib: Current digits at %ld out of %ld\n",
	   DIG2DEC (lrs_record_digits),
	   DIG2DEC (lrs_digits));

  for (i = 0; i < lrs_global_count; i++)
    {
      print_basis (stderr, lrs_global_list[i]);
    }
  fprintf (stderr, "lrs_lib: checkpoint finished\n");
}


/* print out the saved copy of the basis */
void 
print_basis (FILE * fp, lrs_dat * global)
{
  int i;
/* assign local variables to structures */
  fprintf (fp, "lrs_lib: State #%ld: (%s)\t", global->id, global->name);

  if (global->saved_flag)
    {

      fprintf (fp, "V#%ld R#%ld B#%ld h=%ld facets ",
	       global->saved_count[1],
	       global->saved_count[0],
	       global->saved_count[2],
	       global->saved_depth);
      for (i = 0; i < global->saved_d; i++)
	fprintf (fp, "%ld ",
		 global->inequality[global->saved_C[i] - global->lastdv]);
      pmp (" det=", global->saved_det);
      fprintf (fp, "\n");

    }
  else
    {
      fprintf (fp, "lrs_lib: Computing initial basis\n");
    }


  fflush (fp);
}

#ifdef SIGNALS

/*
   If given a signal
   USR1            print current cobasis and continue
   TERM            print current cobasis and terminate
   INT (ctrl-C) ditto
   HUP                     ditto
 */
void 
setup_signals ()
{
  errcheck ("signal", signal (SIGTERM, die_gracefully));
  errcheck ("signal", signal (SIGALRM, timecheck));
  errcheck ("signal", signal (SIGHUP, die_gracefully));
  errcheck ("signal", signal (SIGINT, die_gracefully));
  errcheck ("signal", signal (SIGUSR1, checkpoint));
}

void 
timecheck ()
{
  lrs_dump_state ();
  errcheck ("signal", signal (SIGALRM, timecheck));
  alarm (lrs_checkpoint_seconds);
}

void 
checkpoint ()
{
  lrs_dump_state ();
  errcheck ("signal", signal (SIGUSR1, checkpoint));
}

void 
die_gracefully ()
{
  lrs_dump_state ();

  exit (1);
}

#endif

#ifdef TIMES
/* 
 * Not sure about the portability of this yet, 
 *              - db
 */
#include <sys/resource.h>
#define double_time(t) ((double)(t.tv_sec)+(double)(t.tv_usec)/1000000)

void 
ptimes ()
{
  struct rusage rusage;
  getrusage (RUSAGE_SELF, &rusage);
  fprintf (lrs_ofp, "\n*%0.3fu %0.3fs %ldKb %ld flts %ld swaps %ld blks-in %ld blks-out \n",
	   double_time (rusage.ru_utime),
	   double_time (rusage.ru_stime),
	   rusage.ru_maxrss, rusage.ru_majflt, rusage.ru_nswap,
	   rusage.ru_inblock, rusage.ru_oublock);

}

#endif

/* Routines based on lp_solve */

void
lrs_set_row(lrs_dic *P, lrs_dat *Q, long row, long num[], long den[], long ineq)
/* convert to lrs_mp then call lrs_set_row */
{
 lrs_mp_vector Num, Den;
 long m,d;
 long j;
 
  m = P->m;
  d = P->d;

  Num=lrs_alloc_mp_vector(d+1);
  Den=lrs_alloc_mp_vector(d+1);

  for (j=0;j<=d;j++)
  {
  itomp(num[j],Num[j]);
  itomp(den[j],Den[j]);
  }
 
  lrs_set_row_mp(P,Q,row,Num,Den,ineq);

  lrs_clear_mp_vector(Num,d+1);
  lrs_clear_mp_vector(Den,d+1);

}

void
lrs_set_row_mp(lrs_dic *P, lrs_dat *Q, long row, lrs_mp_vector num, lrs_mp_vector den, long ineq)
/* set row of dictionary using num and den arrays for rational input */
/* ineq = 1 (GE)   - ordinary row  */
/*      = 0 (EQ)   - linearity     */
{
  lrs_mp Temp, mpone;
  lrs_mp_vector oD;             /* denominator for row  */

  long i, j;

/* assign local variables to structures */

  lrs_mp_matrix A;
  lrs_mp_vector Gcd, Lcm;
  long hull;
  long m, d;
  lrs_alloc_mp(Temp); lrs_alloc_mp(mpone);
  hull = Q->hull;
  A = P->A;
  m = P->m;
  d = P->d;
  Gcd = Q->Gcd;
  Lcm = Q->Lcm;

  oD = lrs_alloc_mp_vector (d);
  itomp (ONE, mpone);
  itomp (ONE, oD[0]);

  i=row;
  itomp (ONE, Lcm[i]);      /* Lcm of denominators */
  itomp (ZERO, Gcd[i]);     /* Gcd of numerators */
  for (j = hull; j <= d; j++)       /* hull data copied to cols 1..d */
        {
          copy( A[i][j],num[j-hull]);
          copy(oD[j],den[j-hull]);
          if (!one(oD[j]))
            lcm (Lcm[i], oD[j]);      /* update lcm of denominators */
          copy (Temp, A[i][j]);
          gcd (Gcd[i], Temp);   /* update gcd of numerators   */
        }

  if (hull)
        {
          itomp (ZERO, A[i][0]);        /*for hull, we have to append an extra column of zeroes */
          if (!one (A[i][1]) || !one (oD[1]))         /* all rows must have a one in column one */
            Q->polytope = FALSE;
        }
  if (!zero (A[i][hull]))   /* for H-rep, are zero in column 0     */
        Q->homogeneous = FALSE; /* for V-rep, all zero in column 1     */

  storesign (Gcd[i], POS);
  storesign (Lcm[i], POS);
  if (greater (Gcd[i], mpone) || greater (Lcm[i], mpone))
        for (j = 0; j <= d; j++)
          {
            exactdivint (A[i][j], Gcd[i], Temp);        /*reduce numerators by Gcd  */
            mulint (Lcm[i], Temp, Temp);        /*remove denominators */
            exactdivint (Temp, oD[j], A[i][j]);       /*reduce by former denominator */
          }

  if ( ineq == EQ )        /* input is linearity */
     {
      Q->linearity[Q->nlinearity]=row;
      Q->nlinearity++;
     }


  lrs_clear_mp_vector (oD,d);
  lrs_clear_mp(Temp); lrs_clear_mp(mpone);
}          /* end of lrs_set_row_mp */

void 
lrs_set_obj(lrs_dic *P, lrs_dat *Q, long num[], long den[], long max)
{
  long i;

  if (max == MAXIMIZE)
       Q->maximize=TRUE;
  else
       {
       Q->minimize=TRUE;
       for(i=0;i<=P->d;i++)
         num[i]=-num[i];
       }

  lrs_set_row(P,Q,0L,num,den,GE);
}

void
lrs_set_obj_mp(lrs_dic *P, lrs_dat *Q, lrs_mp_vector num, lrs_mp_vector den, long max)
{
  long i;

  if (max == MAXIMIZE)
       Q->maximize=TRUE;
  else
       {
       Q->minimize=TRUE;
       for(i=0;i<=P->d;i++)
         changesign(num[i]);
       }

  lrs_set_row_mp(P,Q,0L,num,den,GE);
}


long
lrs_solve_lp(lrs_dic *P, lrs_dat *Q)
/* user callable function to solve lp only */
{
  lrs_mp_matrix Lin;		/* holds input linearities if any are found             */
  long col;

  Q->lponly = TRUE;

  if (!lrs_getfirstbasis (&P, Q, &Lin, FALSE))
    return FALSE;

/* There may have been column redundancy                */
/* If so the linearity space is obtained and redundant  */
/* columns are removed. User can access linearity space */
/* from lrs_mp_matrix Lin dimensions nredundcol x d+1   */

  for (col = 0; col < Q->nredundcol; col++)	/* print linearity space               */
    lrs_printoutput (Q, Lin[col]);   	        /* Array Lin[][] holds the coeffs.     */

  return TRUE;
} /* end of lrs_solve_lp */


long 
dan_selectpivot (lrs_dic * P, lrs_dat * Q, long *r, long *s)
/* select pivot indices using dantzig simplex method             */
/* largest coefficient with lexicographic rule to avoid cycling  */
/* Bohdan Kaluzny's handiwork                                    */
/* returns TRUE if pivot found else FALSE                        */
/* pivot variables are B[*r] C[*s] in locations Row[*r] Col[*s]  */
{
  long j,k,col;
  lrs_mp coeff;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Col = P->Col;
  long d = P->d;

  lrs_alloc_mp (coeff);
  *r = 0;
  *s = d;
  j = 0;
  k = 0;
 
  itomp(0,coeff);
/*find positive cost coef */
  while (k < d) 
     {
       if(greater(A[0][Col[k]],coeff))
        {
          j = k;
          copy(coeff,A[0][Col[j]]);
	}
      k++;
     }

  if (positive(coeff))			/* pivot column found! */
    {
      *s = j;
      col = Col[j];

      /*find min index ratio */
      *r = ratio (P, Q, col);
      if (*r != 0)
        {
        lrs_clear_mp(coeff);
	return (TRUE);		/* unbounded */
        }
    }
  lrs_clear_mp(coeff);
  return (FALSE);
}				/* end of dan_selectpivot        */

long
lrs_set_digits(long dec_digits)
{
/* convert user specified decimal digits to mp digits */

  fprintf (lrs_ofp, "\n*digits %ld", dec_digits);
  if (dec_digits > 0)
    lrs_digits = DEC2DIG (dec_digits);
  if (lrs_digits > MAX_DIGITS)
    {
      fprintf (lrs_ofp, "\nDigits must be at most %ld\nChange MAX_DIGITS and recompile",
	       DIG2DEC (MAX_DIGITS));
      fflush(stdout);
      return (FALSE);
    }
  return (TRUE);
}

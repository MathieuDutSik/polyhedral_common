// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_CDDLIB_H_
#define SRC_POLY_POLY_CDDLIB_H_
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

/* The original author of the CDD code is Komei Fukuda.
   This C++ version was adapted by Mathieu Dutour Sikiric
   from the original C version.
*/

// clang-format off
#include "Boost_bitset.h"
#include "POLY_LinearProgrammingFund.h"
#include "POLY_Fundamental.h"
#include "MAT_Matrix.h"
#include <cassert>
#include <vector>
#include <utility>
#include <unordered_map>
#include <string>
// clang-format on

#ifdef DEBUG
#define DEBUG_CDD
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_CDD
#endif

#ifdef TIMINGS
#define TIMINGS_CDD
#endif

// templatized code of the CDD library

namespace cdd {
typedef enum { dd_CrissCross, dd_DualSimplex } dd_LPSolverType;
const dd_LPSolverType dd_choiceLPSolverDefault =
    dd_DualSimplex; /* Default LP solver Algorithm */
const dd_LPSolverType dd_choiceRedcheckAlgorithm =
    dd_DualSimplex; /* Redundancy Checking Algorithm */
const bool dd_choiceLexicoPivotQ =
    true; /* whether to use the lexicographic pivot */
typedef unsigned char set_card_lut_t;

template <typename T> void dd_WriteT(std::ostream &os, T *a, int d) {
  os << "CDD: dd_WriteT a=";
  for (int i = 0; i < d; i++)
    os << "  " << a[i];
  os << "\n";
}

template <typename T> inline bool dd_LargerFrac(T val1, T q1, T val2, T q2) {
#ifdef DEBUG_CDD
  assert(q1 > 0);
  assert(q2 > 0);
#endif
  return val1 * q2 > val2 * q1;
}

template <typename T>
inline bool dd_SmallerFrac(T val1, T q1, T val2, T q2)
// We want to have val1 / q1 < val2 / q2 which is equivalent to val1 * q2 < val2
// * q1
{
#ifdef DEBUG_CDD
  assert(q1 > 0);
  assert(q2 > 0);
#endif
  return val1 * q2 < val2 * q1;
}

template <typename T> inline bool dd_EqualFrac(T val1, T q1, T val2, T q2) {
#ifdef DEBUG_CDD
  assert(q1 > 0);
  assert(q2 > 0);
#endif
  return val1 * q2 < val2 * q1;
}

template <typename T> inline void dd_abs(T &absval, T val) {
  if (val < 0)
    absval = -val;
  else
    absval = val;
}

typedef long dd_rowrange;
typedef long dd_colrange;
typedef long dd_bigrange;

typedef long *dd_rowindex;
typedef int *dd_rowflag;
typedef long *dd_colindex;

// API for the set systems.
typedef unsigned long *set_type; /* set type definition */

#define SETBITS (sizeof(long) * 8)

#define LUTBLOCKS(set)                                                         \
  (((set[0] - 1) / SETBITS + 1) * (sizeof(long) / sizeof(set_card_lut_t)))

/* It is a vector of length 256. For each entry we give
   the number of entries that are matching. This is very
   useful for the computation of cardinality.
 */
static unsigned char set_card_lut[] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};
/* End of Definitions for optimized set_card */

inline void set_free(set_type set) { delete[] set; }

inline unsigned long set_blocks(long len) {
  long blocks = 1L;

  if (len > 0) {
    blocks = (len - 1) / SETBITS + 2;
  }
  return blocks;
}

void set_initialize(set_type *setp, long length)
/* Make a set with a given bit lengths  */
{
  long i, forlim1, len;

  if (length <= 0)
    len = 1;
  else
    len = length;
  /* if negative length is requested, it generates the shortest length */

  forlim1 = set_blocks(len);
  using ulong = unsigned long;
  *setp = new ulong[forlim1];
  (*setp)[0] = ulong(len); /* size of the ground set */
  for (i = 1; i < forlim1; i++)
    (*setp)[i] = 0;
}

void reset_initialize(set_type *setp, long length)
/* Make a set with a given bit lengths  */
{
  long i, forlim1, len;
  if (length <= 0)
    len = 1;
  else
    len = length;

  forlim1 = set_blocks(len);
  using ulong = unsigned long;
  (*setp)[0] = ulong(len); /* size of the ground set */
  for (i = 1; i < forlim1; i++)
    (*setp)[i] = 0;
}

void set_emptyset(set_type set)
/* Set set to be the emptyset  */
{
  long forlim = set_blocks(set[0]) - 1;
  for (long i = 1; i <= forlim; i++)
    set[i] = 0;
}

void set_copy(set_type setcopy, set_type set)
/* Copy the set set[] to setcopy[] with setcopy[] length */
{
  long forlim = set_blocks(setcopy[0]) - 1;
  for (long i = 1; i <= forlim; i++)
    setcopy[i] = set[i];
}

inline void set_addelem(set_type set, long elem)
/* add elem only if it is within the set[] range */
{
  long i, j;
  unsigned long change;
  unsigned long one = 1U;
  if (elem <= (long)set[0]) {
    i = (elem - 1) / SETBITS + 1;
    j = (elem - 1) % SETBITS;
    change = one << j; /* put 1 in jth position */
    set[i] = set[i] | change;
  }
}

inline void set_delelem(set_type set, long elem)
/* delete elem only if it is within the set[] range */
{
  long i, j;
  unsigned long change;
  unsigned long one = 1U;
  if (elem <= (long)set[0]) {
    i = (elem - 1) / SETBITS + 1;
    j = (elem - 1) % SETBITS;
    change = one << j; /* put 1 in jth position */
    set[i] = (set[i] | change) ^ change;
  }
}

inline void set_int(set_type set, set_type set1, set_type set2)
/* Set intersection, assuming set1 and set2 have the same length as set */
{
  long forlim = set_blocks(set[0]) - 1;
  for (long i = 1; i <= forlim; i++)
    set[i] = (set1[i] & set2[i]);
}

inline void set_uni(set_type set, set_type set1, set_type set2)
/* Set union,assuming set1 and set2 have the same length as set */
{
  long forlim = set_blocks(set[0]) - 1;
  for (long i = 1; i <= forlim; i++)
    set[i] = set1[i] | set2[i];
}

inline void set_diff(set_type set, set_type set1, set_type set2)
/* Set difference se1/set2, assuming set1 and set2 have the same length as set
 */
{
  long forlim = set_blocks(set[0]) - 1;
  for (long i = 1; i <= forlim; i++)
    set[i] = set1[i] & (~set2[i]);
}

inline void set_compl(set_type set, set_type set1)
/* set[] will be set to the complement of set1[] */
{
  long i, j, l, forlim;
  unsigned long change;
  unsigned long one = 1U;

  forlim = set_blocks(set[0]) - 1;
  for (i = 1; i <= forlim; i++)
    set[i] = ~set1[i];
  l = (set[0] - 1) %
      SETBITS; /* the position of the last elem in the last block */
  for (j = l + 1; j <= (long)SETBITS - 1; j++) {
    change = one << j;
    set[forlim] = (set[forlim] | change) ^ change;
  }
}

inline bool set_subset(set_type set1, set_type set2)
/* Set containment check, set1 <= set2 */
{
  bool reply = true;
  long i, forlim;
  forlim = set_blocks(set2[0]) - 1;
  for (i = 1; i <= forlim && reply; i++)
    if ((set1[i] | set2[i]) != set2[i])
      reply = false;
  return reply;
}

inline bool set_member(long elem, set_type set)
/* Set membership check, elem in set */
{
  bool reply = false;
  long i, j;
  unsigned long testset;
  unsigned long one = 1U;
  if (elem <= (long)set[0]) {
    i = (elem - 1) / SETBITS + 1;
    j = (elem - 1) % SETBITS;
    testset = set[i] | (one << j); /* add elem to set[i] */
    if (testset == set[i])
      reply = true;
  }
  return reply;
}

/*set cardinality, modified by David Bremner bremner@cs.mcgill.ca
   to optimize for speed.*/
inline long set_card(set_type set) {
  unsigned long block;
  long car = 0;
  set_card_lut_t *p;
  p = (set_card_lut_t *)&set[1];
  for (block = 0; block < LUTBLOCKS(set); block++)
    car += set_card_lut[p[block]];
  return car;
}

typedef set_type dd_rowset;
typedef set_type dd_colset;
typedef set_type *dd_SetVector;
typedef set_type *dd_Aincidence;
typedef set_type rowset; /* set_type defined in setoper.h */
typedef set_type colset;

// More functionalities

template <typename T>
inline void dd_InnerProduct(T &prod, dd_colrange d, T *v1, T *v2) {
  dd_colrange j;
  prod = 0;
  for (j = 0; j < d; j++)
    prod += v1[j] * v2[j];
}

template <typename T> struct dd_raydata {
  T *Ray;
  dd_rowset ZeroSet;
  dd_rowrange FirstInfeasIndex; /* the first inequality the ray violates */
  bool feasible;                /* flag to store the feasibility */
  T ARay;                       /* temporary area to store some row of A*Ray */
  dd_raydata<T> *Next;
};

template <typename T> struct dd_adjacencydata {
  dd_raydata<T> *Ray1;
  dd_raydata<T> *Ray2;
  dd_adjacencydata<T> *Next;
};

typedef enum { dd_Combinatorial, dd_Algebraic } dd_AdjacencyTestType;

typedef enum {
  dd_MaxIndex,
  dd_MinIndex,
  dd_MinCutoff,
  dd_MaxCutoff,
  dd_MixCutoff,
  dd_LexMin,
  dd_LexMax,
  dd_RandomRow
} dd_RowOrderType;

typedef enum {
  dd_Unspecified = 0,
  dd_Inequality,
  dd_Generator
} dd_RepresentationType;

typedef enum {
  dd_IneToGen,
  dd_GenToIne,
  dd_LPMax,
  dd_LPMin,
  dd_InteriorFind
} dd_ConversionType;

typedef enum {
  dd_IncOff = 0,
  dd_IncCardinality,
  dd_IncSet
} dd_IncidenceOutputType;

typedef enum {
  dd_AdjOff = 0,
  dd_AdjacencyList,
  dd_AdjacencyDegree
} dd_AdjacencyOutputType;

typedef enum { dd_Auto, dd_SemiAuto, dd_Manual } dd_FileInputModeType;
/* Auto if a input filename is specified by command arguments */

typedef enum {
  dd_DimensionTooLarge,
  dd_ImproperInputFormat,
  dd_NegativeMatrixSize,
  dd_EmptyVrepresentation,
  dd_EmptyHrepresentation,
  dd_EmptyRepresentation,
  dd_IFileNotFound,
  dd_OFileNotOpen,
  dd_NoLPObjective,
  dd_NoRealNumberSupport,
  dd_NotAvailForH,
  dd_NotAvailForV,
  dd_CannotHandleLinearity,
  dd_RowIndexOutOfRange,
  dd_ColIndexOutOfRange,
  dd_LPCycling,
  dd_NumericallyInconsistent,
  dd_NonZeroLinearity,
  dd_NoError
} dd_ErrorType;

void dd_WriteErrorMessages(std::ostream &os, dd_ErrorType Error) {
  switch (Error) {

  case dd_NonZeroLinearity:
    os << "*Input Error: Input entry has a non-zero linearity:\n";
    os << "*The function does not support this.\n";
    break;

  case dd_DimensionTooLarge:
    os << "*Input Error: Input matrix is too large:\n";
    os << "*Please increase MMAX and/or NMAX in the source code and "
          "recompile.\n";
    break;

  case dd_IFileNotFound:
    os << "*Input Error: Specified input file does not exist.\n";
    break;

  case dd_OFileNotOpen:
    os << "*Output Error: Specified output file cannot be opened.\n";
    break;

  case dd_NegativeMatrixSize:
    os << "*Input Error: Input matrix has a negative size:\n";
    os << "*Please check rowsize or colsize.\n";
    break;

  case dd_ImproperInputFormat:
    os << "*Input Error: Input format is not correct.\n";
    os << "*Format:\n";
    os << " begin\n";
    os << "   m   n  NumberType(real, rational or integer)\n";
    os << "   b  -A\n";
    os << " end\n";
    break;

  case dd_EmptyVrepresentation:
    os << "*Input Error: V-representation is empty:\n";
    os << "*cddlib does not accept this trivial case for which output can be "
          "any inconsistent system.\n";
    break;

  case dd_EmptyHrepresentation:
    os << "*Input Error: H-representation is empty.\n";
    break;

  case dd_EmptyRepresentation:
    os << "*Input Error: Representation is empty.\n";
    break;

  case dd_NoLPObjective:
    os << "*LP Error: No LP objective (max or min) is set.\n";
    break;

  case dd_NoRealNumberSupport:
    os << "*LP Error: The binary (with GMP Rational) does not support Real "
          "number input.\n";
    os << "         : Use a binary compiled without -DGMPRATIONAL option.\n";
    break;

  case dd_NotAvailForH:
    os << "*Error: A function is called with H-rep which does not support an "
          "H-representation.\n";
    break;

  case dd_NotAvailForV:
    os << "*Error: A function is called with V-rep which does not support an "
          "V-representation.\n";
    break;

  case dd_CannotHandleLinearity:
    os << "*Error: The function called cannot handle linearity.\n";
    break;

  case dd_RowIndexOutOfRange:
    os << "*Error: Specified row index is out of range\n";
    break;

  case dd_ColIndexOutOfRange:
    os << "*Error: Specified column index is out of range\n";
    break;

  case dd_LPCycling:
    os << "*Error: Possibly an LP cycling occurs.  Use the Criss-Cross "
          "method.\n";
    break;

  case dd_NumericallyInconsistent:
    os << "*Error: Numerical inconsistency is found.  Use the GMP exact "
          "arithmetic.\n";
    break;

  case dd_NoError:
    os << "*No Error found.\n";
    break;
  }
}

typedef enum { dd_InProgress, dd_AllFound, dd_RegionEmpty } dd_CompStatusType;

/* --- LP types ---- */

typedef enum { dd_LPnone = 0, dd_LPmax, dd_LPmin } dd_LPObjectiveType;

typedef enum {
  dd_LPSundecided,
  dd_Optimal,
  dd_Inconsistent,
  dd_DualInconsistent,
  dd_StrucInconsistent,
  dd_StrucDualInconsistent,
  dd_Unbounded,
  dd_DualUnbounded,
  dd_TooManyIterations
} dd_LPStatusType;

template <typename T> struct dd_lpdata {
  //  dd_DataFileType filename;
  dd_LPObjectiveType objective;
  dd_LPSolverType solver;
  dd_rowrange m;
  dd_colrange d;
  T **A;
  T **B;
  dd_rowrange objrow;
  dd_colrange rhscol;
  dd_rowset equalityset;

  bool redcheck_extensive; /* Apply the extensive redundancy check. */
  dd_rowset redset_extra; /* a set of rows that are newly recognized redundan by
                             the extensive search. */
  dd_rowset redset_accum; /* the accumulated set of rows that are recognized
                             redundant */
  dd_rowset posset_extra; /* a set of rows that are recognized non-linearity */

  bool lexicopivot; /* flag to use the lexicogrphic pivot rule (symbolic
                       perturbation). */

  dd_LPStatusType LPS; /* the current solution status */
  dd_rowrange m_alloc; /* the allocated row size of matrix A */
  dd_colrange d_alloc; /* the allocated col size of matrix A */
  T optvalue;          /* optimal value */
  T *sol;              /* primal solution */
  T *dsol;             /* dual solution */
  dd_colindex nbindex; /* current basis represented by nonbasic indices */
  dd_rowrange re; /* row index as a certificate in the case of inconsistency */
  dd_colrange
      se; /* col index as a certificate in the case of dual inconsistency */
  long pivots[5];
  /* pivots[0]=setup (to find a basis), pivots[1]=PhaseI or Criss-Cross,
     pivots[2]=Phase II, pivots[3]=Anticycling, pivots[4]=GMP postopt. */
  long total_pivots;
  int use_given_basis;       /* switch to indicate the use of the given basis */
  dd_colindex given_nbindex; /* given basis represented by nonbasic indices */
};

template <typename T>
void dd_WriteAmatrix(std::ostream &os, T **A, long rowmax, long colmax) {
  long i, j;
  if (A == NULL) {
    os << "WriteAmatrix: The requested matrix is empty\n";
    return;
  }
  os << "begin\n";
  os << " " << rowmax << " " << colmax << " rational\n";
  for (i = 1; i <= rowmax; i++) {
    for (j = 1; j <= colmax; j++) {
      os << " " << A[i - 1][j - 1];
    }
    os << "\n";
  }
  os << "end\n";
}

template <typename T> void dd_WriteArow(std::ostream &os, T *a, dd_colrange d) {
  dd_colrange j;
  for (j = 0; j < d; j++)
    os << " " << a[j];
  os << "\n";
}

template <typename T> void dd_WriteLP(std::ostream &os, dd_lpdata<T> *lp) {
  if (lp == nullptr) {
    os << "WriteLP: The requested lp is empty\n";
    return;
  }
  os << "H-representation\n";

  dd_WriteAmatrix(os, lp->A, (lp->m) - 1, lp->d);
  if (lp->objective != dd_LPnone) {
    if (lp->objective == dd_LPmax)
      os << "maximize\n";
    else
      os << "minimize\n";
    dd_WriteArow(os, lp->A[lp->objrow - 1], lp->d);
  }
}

template <typename T>
void dd_WriteLPResult(std::ostream &os, dd_lpdata<T> *lp, dd_ErrorType err) {
  long j;
  os << "* cdd LP solver result\n";

  if (err != dd_NoError) {
    dd_WriteErrorMessages(os, err);
    return;
  }

  os << "* #constraints = " << (lp->m - 1) << "\n";
  os << "* #variables   = " << (lp->d - 1) << "\n";

  switch (lp->solver) {
  case dd_DualSimplex:
    os << "* Algorithm: dual simplex algorithm\n";
    break;
  case dd_CrissCross:
    os << "* Algorithm: criss-cross method\n";
    break;
  }

  switch (lp->objective) {
  case dd_LPmax:
    os << "* maximization is chosen\n";
    break;
  case dd_LPmin:
    os << "* minimization is chosen\n";
    break;
  case dd_LPnone:
    os << "* no objective type (max or min) is chosen\n";
    break;
  }

  if (lp->objective == dd_LPmax || lp->objective == dd_LPmin) {
    os << "* Objective function is\n";
    for (j = 0; j < lp->d; j++) {
      if (j > 0 && lp->A[lp->objrow - 1][j] >= 0)
        os << " +";
      if (j > 0 && (j % 5) == 0)
        os << "\n";
      os << lp->A[lp->objrow - 1][j];
      if (j > 0)
        os << " X[" << j << "]";
    }
    os << "\n";
  }

  switch (lp->LPS) {
  case dd_Optimal:
    os << "* LP status: a dual pair (x,y) of optimal solutions found.\n";
    os << "begin\n";
    os << "  primal_solution\n";
    for (j = 1; j < lp->d; j++) {
      os << "  " << j << " : " << lp->sol[j] << "\n";
    }
    os << "  dual_solution\n";
    for (j = 1; j < lp->d; j++)
      if (lp->nbindex[j + 1] > 0)
        os << "  " << lp->nbindex[j + 1] << " : " << lp->dsol[j] << "\n";
    os << "  optimal_value : " << lp->optvalue << "\n";
    os << "end\n";
    break;

  case dd_Inconsistent:
    os << "* LP status: LP is inconsistent.\n";
    os << "* The positive combination of original inequalities with\n";
    os << "* the following coefficients will prove the inconsistency.\n";
    os << "begin\n";
    os << "  dual_direction\n";
    os << "  " << lp->re << " : 1\n";
    for (j = 1; j < lp->d; j++)
      if (lp->nbindex[j + 1] > 0)
        os << "  " << lp->nbindex[j + 1] << " : " << lp->dsol[j] << "\n";
    os << "end\n";
    break;
  case dd_DualInconsistent:
  case dd_StrucDualInconsistent:
    os << "* LP status: LP is dual inconsistent.\n";
    os << "* The linear combination of columns with\n";
    os << "* the following coefficients will prove the dual inconsistency.\n";
    os << "* (It is also an unbounded direction for the primal LP.)\n";
    os << "begin\n";
    os << "  primal_direction\n";
    for (j = 1; j < lp->d; j++)
      os << "  " << j << " : " << lp->sol[j] << "\n";
    os << "end\n";
    break;

  default:
    break;
  }
  os << "* number of pivot operations = " << lp->total_pivots
     << " (ph0 = " << lp->pivots[0] << ", ph1 = " << lp->pivots[1]
     << ", ph2 = " << lp->pivots[2] << ", ph3 = " << lp->pivots[3]
     << ", ph4 = " << lp->pivots[4] << ")\n";
}

/*----  end of LP Types ----- */

template <typename T> struct data_temp_simplex {
  T *rcost;
  dd_colset tieset, stieset;
  T *Rtemp;
  dd_colindex nbindex_ref;
  dd_rowindex bflag;
  dd_rowindex OrderVector; /* the permutation vector to store a preordered row
                              indeces */
  dd_colindex nbindex_ref_ds; /* to be used to store the initial feasible basis
                                 for lexico rule */
};

template <typename T>
data_temp_simplex<T> *allocate_data_simplex(dd_rowrange m_size,
                                            dd_colrange d_size) {
  //  std::cout << "allocate_data_simplex with d_size=" << d_size << "\n";
  data_temp_simplex<T> *data = new data_temp_simplex<T>;
  data->rcost = new T[d_size];
  set_initialize(&data->tieset, d_size);
  set_initialize(&data->stieset, d_size);
  data->Rtemp = new T[d_size];
  data->nbindex_ref = new long[d_size + 1];
  //
  data->OrderVector = new long[m_size + 1];
  data->bflag = new long[m_size + 2];
  data->nbindex_ref_ds = new long[d_size + 1];
  return data;
}

template <typename T> void free_data_simplex(data_temp_simplex<T> *data) {
  delete[] data->rcost;
  set_free(data->tieset);
  set_free(data->stieset);
  delete[] data->Rtemp;
  delete[] data->nbindex_ref;
  delete[] data->OrderVector;
  delete[] data->bflag;
  delete[] data->nbindex_ref_ds;
  delete data;
}

/* Containing data sets */

template <typename T> void dd_FreeAmatrix(dd_rowrange m, T **A) {
  dd_rowrange i;

  if (A != nullptr) {
    for (i = 0; i < m; i++)
      delete[] A[i];
    delete[] A;
  }
}

template <typename T> void dd_FreeBmatrix(dd_colrange d, T **B) {
  dd_colrange j;

  if (B != nullptr) {
    for (j = 0; j < d; j++)
      delete[] B[j];
    delete[] B;
  }
}

template <typename T> struct dd_matrixdata {
  dd_rowrange rowsize;
  dd_rowset linset;
  /*  a subset of rows of linearity (ie, generators of
      linearity space for V-representation, and equations
      for H-representation. */
  dd_colrange colsize;
  dd_RepresentationType representation;
  T **matrix;
  dd_LPObjectiveType objective;
  T *rowvec;
};

template <typename T> void dd_FreeArow(T *a) { delete[] a; }

template <typename T> void dd_FreeMatrix(dd_matrixdata<T> *M) {
  dd_rowrange m1;

  if (M != nullptr) {
    if (M->rowsize <= 0)
      m1 = 1;
    else
      m1 = M->rowsize;
    dd_FreeAmatrix(m1, M->matrix);
    dd_FreeArow(M->rowvec);
    set_free(M->linset);
    delete M;
  }
}

template <typename T> void dd_FreeShallowMatrix(dd_matrixdata<T> *M) {
  if (M != nullptr) {
    delete[] M->matrix;
    dd_FreeArow(M->rowvec);
    set_free(M->linset);
    delete M;
  }
}

struct dd_setfamily {
  dd_bigrange famsize;
  dd_bigrange setsize;
  dd_SetVector set;
};

void dd_FreeSetFamily(dd_setfamily *F) {
  dd_bigrange i, f1;

  if (F != nullptr) {
    if (F->famsize <= 0)
      f1 = 1;
    else
      f1 = F->famsize;
    /* the smallest created size is one */
    for (i = 0; i < f1; i++)
      set_free(F->set[i]);
    delete[] F->set;
    delete[] F;
  }
}

typedef struct dd_nodedata *dd_NodePtr;
typedef struct dd_nodedata {
  dd_bigrange key;
  dd_NodePtr next;
} dd_NodeType;

typedef struct dd_graphdata *dd_GraphPtr;
typedef struct dd_graphdata {
  dd_bigrange vsize;
  dd_NodePtr *adjlist; /* should be initialized to have vsize components */
} dd_GraphType;

// typedef struct dd_polyhedradata *dd_PolyhedraPtr;
// typedef struct dd_conedata *dd_ConePtr;
template <typename T> struct dd_polyhedradata;

template <typename T> struct dd_conedata {
  dd_RepresentationType representation;
  dd_rowrange m;
  dd_colrange d;
  T **A;
  dd_polyhedradata<T> *parent; /* pointing to the original polyhedra data */
  dd_rowrange m_alloc;         /* allocated row size of matrix A */
  dd_colrange d_alloc;         /* allocated col size of matrix A */

  /* CONTROL: variables to control computation */
  dd_rowrange Iteration;
  dd_RowOrderType HalfspaceOrder;
  dd_raydata<T> *FirstRay;
  dd_raydata<T> *LastRay;
  dd_raydata<T> *ArtificialRay;
  dd_raydata<T> *PosHead;
  dd_raydata<T> *ZeroHead;
  dd_raydata<T> *NegHead;
  dd_raydata<T> *PosLast;
  dd_raydata<T> *ZeroLast;
  dd_raydata<T> *NegLast;
  dd_adjacencydata<T> **Edges; /* adjacency relation storage for iteration k */
  unsigned int rseed;          /* random seed for random row permutation */

  bool ColReduced; /* flag to indicate that a column basis is computed and
                      reduced */
  dd_bigrange LinearityDim;
  /*  the dimension of the linearity space (when input is H), and
      the size of a minimal system of equations to determine the space (when V).
   */
  dd_colrange d_orig; /* the size d of the original matrix A */
  dd_colindex newcol; /* the size d of the original matrix A */

  dd_colindex InitialRayIndex; /* InitialRayIndex[s] (s>=1) stores the corr. row
                                  index */
  dd_rowindex OrderVector;
  bool RecomputeRowOrder;
  bool PreOrderedRun;
  dd_rowset GroundSet, EqualitySet, NonequalitySet, AddedHalfspaces,
      WeaklyAddedHalfspaces, InitialHalfspaces;
  long RayCount, FeasibleRayCount, WeaklyFeasibleRayCount, TotalRayCount,
      ZeroRayCount;
  long EdgeCount, TotalEdgeCount;
  long count_int, count_int_good,
      count_int_bad; /* no. of intersection operations */

  T **B;
  T **Bsave; /* a copy of the dual basis inverse used to reduce the matrix A */

  dd_ErrorType Error;
  dd_CompStatusType CompStatus; /* Computation Status */
};

template <typename T> struct dd_polyhedradata {
  dd_RepresentationType representation; /* given representation */
  bool homogeneous;
  dd_colrange d;
  dd_rowrange m;
  T **A;                 /* Inequality System:  m times d matrix */
  dd_conedata<T> *child; /* pointing to the homogenized cone data */
  dd_rowrange m_alloc;   /* allocated row size of matrix A */
  dd_colrange d_alloc;   /* allocated col size of matrix A */
  T *c;                  /* cost vector */

  dd_rowflag EqualityIndex;
  /* ith component is 1 if it is equality, -1 if it is strict inequality, 0
   * otherwise. */

  bool IsEmpty; /* This is to tell whether the set is empty or not */

  bool NondegAssumed;
  bool InitBasisAtBottom;
  bool RestrictedEnumeration;
  bool RelaxedEnumeration;

  dd_rowrange m1;
  /* = m or m+1 (when representation=Inequality && !homogeneous)
     This data is written after dd_ConeDataLoad is called.  This
     determines the size of Ainc. */
  bool AincGenerated;
  /* Indicates whether Ainc, Ared, Adom are all computed.
     All the variables below are valid only when this is TRUE */
  dd_colrange ldim; /* linearity dimension */
  dd_bigrange n;
  /* the size of output = total number of rays
     in the computed cone + linearity dimension */
  dd_Aincidence Ainc;

  /* incidence of input and output */
  dd_rowset Ared;
  /* redundant set of rows whose removal results in a minimal system */
  dd_rowset Adom;
  /* dominant set of rows (those containing all rays). */
};

dd_setfamily *dd_CreateSetFamily(dd_bigrange fsize, dd_bigrange ssize) {
  dd_bigrange i, f0, f1, s0, s1;

  if (fsize <= 0) {
    f0 = 0;
    f1 = 1;
    /* if fsize<=0, the fsize is set to zero and the created size is one */
  } else {
    f0 = fsize;
    f1 = fsize;
  }
  if (ssize <= 0) {
    s0 = 0;
    s1 = 1;
    /* if ssize<=0, the ssize is set to zero and the created size is one */
  } else {
    s0 = ssize;
    s1 = ssize;
  }

  dd_setfamily *F = new dd_setfamily;
  F->set = new set_type[f1];
  for (i = 0; i < f1; i++)
    set_initialize(&(F->set[i]), s1);
  F->famsize = f0;
  F->setsize = s0;
  return F;
}

template <typename T>
void dd_AppendMatrix2Poly(dd_polyhedradata<T> **poly, dd_matrixdata<T> *M) {
  dd_matrixdata<T> *Mpoly, Mnew = nullptr;
  dd_ErrorType err;

  if ((*poly) != nullptr && (*poly)->m >= 0 && (*poly)->d >= 0 &&
      (*poly)->d == M->colsize && M->rowsize > 0) {
    Mpoly = dd_CopyInput(*poly);
    Mnew = dd_AppendMatrix(Mpoly, M);
    dd_FreePolyhedra(*poly);
    *poly = dd_DDMatrix2Poly(Mnew, &err);
    dd_FreeMatrix(Mpoly);
    dd_FreeMatrix(Mnew);
  }
}

// Basic Initialization procedures

template <typename T>
void dd_AllocateAmatrix(dd_rowrange m, dd_colrange d, T ***A) {
  dd_rowrange i;

  if (m > 0)
    (*A) = new T *[m];
  for (i = 0; i < m; i++)
    (*A)[i] = new T[d];
}

template <typename T> void dd_AllocateArow(dd_colrange d, T **a) {
  if (d > 0)
    *a = new T[d];
}

template <typename T> void dd_AllocateBmatrix(dd_colrange d, T ***B) {
  dd_colrange j;

  (*B) = new T *[d];
  for (j = 0; j < d; j++)
    (*B)[j] = new T[d];
}

// More sophisticated creation functionalities.

template <typename T>
dd_matrixdata<T> *dd_CreateMatrix(dd_rowrange m_size, dd_colrange d_size) {
  dd_matrixdata<T> *M;
  dd_rowrange m0, m1;
  dd_colrange d0, d1;

  if (m_size <= 0) {
    m0 = 0;
    m1 = 1;
    /* if m_size <=0, the number of rows is set to zero, the actual size is 1 */
  } else {
    m0 = m_size;
    m1 = m_size;
  }
  if (d_size <= 0) {
    d0 = 0;
    d1 = 1;
    /* if d_size <=0, the number of cols is set to zero, the actual size is 1 */
  } else {
    d0 = d_size;
    d1 = d_size;
  }
  M = new dd_matrixdata<T>;
  dd_AllocateAmatrix(m1, d1, &(M->matrix));
  dd_AllocateArow(d1, &(M->rowvec));
  M->rowsize = m0;
  set_initialize(&(M->linset), m1);
  M->colsize = d0;
  M->objective = dd_LPnone;
  M->representation = dd_Unspecified;
  return M;
}

/* This is a kind of shallow copy. The matrix entries are not assigned but
   supposed to be obtained via pointer assignment from another matrix. However,
   we may lose some speed due to alignment stuff, so it is unsure if this makes
   sense at all.
*/
template <typename T>
dd_matrixdata<T> *dd_CreateShallowMatrix(dd_rowrange m_size,
                                         dd_colrange d_size) {
  dd_matrixdata<T> *M;
  dd_rowrange m0, m1;
  dd_colrange d0, d1;

  if (m_size <= 0) {
    m0 = 0;
    m1 = 1;
    /* if m_size <=0, the number of rows is set to zero, the actual size is 1 */
  } else {
    m0 = m_size;
    m1 = m_size;
  }
  if (d_size <= 0) {
    d0 = 0;
    d1 = 1;
    /* if d_size <=0, the number of cols is set to zero, the actual size is 1 */
  } else {
    d0 = d_size;
    d1 = d_size;
  }
  M = new dd_matrixdata<T>;
  M->matrix = new T *[m1];
  dd_AllocateArow(d1, &(M->rowvec));
  M->rowsize = m0;
  set_initialize(&(M->linset), m1);
  M->colsize = d0;
  M->objective = dd_LPnone;
  M->representation = dd_Unspecified;
  return M;
}

template <typename T>
dd_polyhedradata<T> *dd_CreatePolyhedraData(dd_rowrange m, dd_colrange d) {
  dd_rowrange i;
  dd_polyhedradata<T> *poly = nullptr;

  poly = new dd_polyhedradata<T>;
  poly->child = nullptr; /* this links the homogenized cone data */
  poly->c = nullptr;     /* Added to avoid some ambiguity */
  poly->m = m;
  poly->d = d;
  poly->n = -1;          /* the size of output is not known */
  poly->m_alloc = m + 2; /* the allocated row size of matrix A */
  poly->d_alloc = d;     /* the allocated col size of matrix A */
  poly->ldim = 0;        /* initialize the linearity dimension */
  dd_AllocateAmatrix(poly->m_alloc, poly->d_alloc, &(poly->A));
  dd_AllocateArow(d, &(poly->c)); /* cost vector */
  poly->representation = dd_Inequality;
  poly->homogeneous = false;

  poly->EqualityIndex = new int[m + 2];
  /* size increased to m+2 in 092b because it is used by the child cone,
     This is a bug fix suggested by Thao Dang. */
  /* ith component is 1 if it is equality, -1 if it is strict inequality, 0
   * otherwise. */
  for (i = 0; i <= m + 1; i++)
    poly->EqualityIndex[i] = 0;

  poly->IsEmpty =
      -1; /* initially set to -1, neither TRUE nor FALSE, meaning unknown */

  poly->NondegAssumed = false;
  poly->InitBasisAtBottom = false;
  poly->RestrictedEnumeration = false;
  poly->RelaxedEnumeration = false;

  poly->AincGenerated =
      false; /* Ainc is a set array to store the input incidence. */

  return poly;
}

template <typename T>
void dd_InitializeConeData(dd_rowrange m, dd_colrange d,
                           dd_conedata<T> **cone) {
  dd_colrange j;

  (*cone) = new dd_conedata<T>;

  /* INPUT: A given representation of a cone: inequality */
  (*cone)->m = m;
  (*cone)->d = d;
  (*cone)->m_alloc = m + 2; /* allocated row size of matrix A */
  (*cone)->d_alloc = d;     /* allocated col size of matrix A, B and Bsave */
  (*cone)->parent = nullptr;

  /* CONTROL: variables to control computation */
  (*cone)->Iteration = 0;

  (*cone)->HalfspaceOrder = dd_LexMin;

  (*cone)->ArtificialRay = nullptr;
  (*cone)->FirstRay = nullptr;
  (*cone)->LastRay = nullptr; /* The second description: Generator */
  (*cone)->PosHead = nullptr;
  (*cone)->ZeroHead = nullptr;
  (*cone)->NegHead = nullptr;
  (*cone)->PosLast = nullptr;
  (*cone)->ZeroLast = nullptr;
  (*cone)->NegLast = nullptr;
  (*cone)->RecomputeRowOrder = true;
  (*cone)->PreOrderedRun = false;
  set_initialize(&((*cone)->GroundSet), (*cone)->m_alloc);
  set_initialize(&((*cone)->EqualitySet), (*cone)->m_alloc);
  set_initialize(&((*cone)->NonequalitySet), (*cone)->m_alloc);
  set_initialize(&((*cone)->AddedHalfspaces), (*cone)->m_alloc);
  set_initialize(&((*cone)->WeaklyAddedHalfspaces), (*cone)->m_alloc);
  set_initialize(&((*cone)->InitialHalfspaces), (*cone)->m_alloc);
  (*cone)->RayCount = 0;
  (*cone)->FeasibleRayCount = 0;
  (*cone)->WeaklyFeasibleRayCount = 0;
  (*cone)->TotalRayCount = 0;
  (*cone)->ZeroRayCount = 0;
  (*cone)->EdgeCount = 0;
  (*cone)->TotalEdgeCount = 0;
  (*cone)->count_int = 0;
  (*cone)->count_int_good = 0;
  (*cone)->count_int_bad = 0;
  (*cone)->rseed = 1; /* random seed for random row permutation */

  dd_AllocateBmatrix((*cone)->d_alloc, &((*cone)->B));
  dd_AllocateBmatrix((*cone)->d_alloc, &((*cone)->Bsave));
  dd_AllocateAmatrix((*cone)->m_alloc, (*cone)->d_alloc, &((*cone)->A));

  (*cone)->Edges = new dd_adjacencydata<T> *[(*cone)->m_alloc];
  for (j = 0; j < (*cone)->m_alloc; j++)
    (*cone)->Edges[j] = nullptr;
  (*cone)->InitialRayIndex = new long[d + 1];
  for (int i = 0; i <= d; i++)
    (*cone)->InitialRayIndex[i] = 0;
  (*cone)->OrderVector = new long[(*cone)->m_alloc + 1];
  for (int i = 0; i <= (*cone)->m_alloc; i++)
    (*cone)->OrderVector[i] = 0;

  (*cone)->newcol = new long[((*cone)->d) + 1];
  for (int i = 0; i <= (*cone)->d; i++)
    (*cone)->newcol[i] = 0;
  for (j = 0; j <= (*cone)->d; j++)
    (*cone)->newcol[j] = j;   /* identity map, initially */
  (*cone)->LinearityDim = -2; /* -2 if it is not computed */
  (*cone)->ColReduced = false;
  (*cone)->d_orig = d;
}

template <typename T>
dd_conedata<T> *dd_ConeDataLoad(dd_polyhedradata<T> *poly) {
  dd_conedata<T> *cone = nullptr;
  dd_colrange d, j;
  dd_rowrange m, i;

  m = poly->m;
  d = poly->d;
  if (!(poly->homogeneous) && poly->representation == dd_Inequality)
    m = poly->m + 1;
  poly->m1 = m;

  dd_InitializeConeData(m, d, &cone);
  cone->representation = poly->representation;

  /* Points to the original polyhedra data, and reversely */
  cone->parent = poly;
  poly->child = cone;

  for (i = 1; i <= poly->m; i++)
    for (j = 1; j <= cone->d; j++)
      cone->A[i - 1][j - 1] = poly->A[i - 1][j - 1];

  if (poly->representation == dd_Inequality && !poly->homogeneous) {
    cone->A[m - 1][0] = 1;
    for (j = 1; j < d; j++)
      cone->A[m - 1][j] = 0;
  }

  return cone;
}

template <typename T>
dd_lpdata<T> *dd_CreateLPData(dd_rowrange m, dd_colrange d) {
  dd_lpdata<T> *lp;
  lp = new dd_lpdata<T>;
  lp->solver = dd_choiceLPSolverDefault; /* set the default lp solver */
  lp->d = d;
  lp->m = m;
  lp->objrow = m;
  lp->rhscol = 1L;
  lp->objective = dd_LPnone;
  lp->LPS = dd_LPSundecided;

  lp->nbindex = new long[d + 1];
  lp->given_nbindex = new long[d + 1];
  for (int i = 0; i <= d; i++)
    lp->nbindex[i] = 0;
  for (int i = 0; i <= d; i++)
    lp->given_nbindex[i] = 0;
  set_initialize(&(lp->equalityset), m);
  /* i must be in the set iff i-th row is equality . */

  lp->redcheck_extensive = false; /* this is on only for RedundantExtensive */
  set_initialize(&(lp->redset_extra), m);
  /* i is in the set if i-th row is newly recognized redundant (during the
   * checking the row ired). */
  set_initialize(&(lp->redset_accum), m);
  /* i is in the set if i-th row is recognized redundant (during the checking
   * the row ired). */
  set_initialize(&(lp->posset_extra), m);
  /* i is in the set if i-th row is recognized non-linearity (during the course
   * of computation). */
  lp->lexicopivot = dd_choiceLexicoPivotQ; /* dd_choice... is set in
                                              dd_set_global_constants() */

  lp->m_alloc = lp->m + 2;
  lp->d_alloc = lp->d + 2;
  dd_AllocateBmatrix(lp->d_alloc, &(lp->B));
  dd_AllocateAmatrix(lp->m_alloc, lp->d_alloc, &(lp->A));
  dd_AllocateArow(lp->d_alloc, &(lp->sol));
  dd_AllocateArow(lp->d_alloc, &(lp->dsol));
  return lp;
}

template <typename T> void dd_LPData_reset_m(dd_rowrange m, dd_lpdata<T> *lp) {
  dd_colrange d = lp->d;
  lp->m = m;
  lp->objrow = m;
  for (int i = 0; i <= d; i++)
    lp->nbindex[i] = 0;
  for (int i = 0; i <= d; i++)
    lp->given_nbindex[i] = 0;
  reset_initialize(&(lp->equalityset), m);
  reset_initialize(&(lp->redset_extra), m);
  reset_initialize(&(lp->redset_accum), m);
  reset_initialize(&(lp->posset_extra), m);
}

template <typename T> inline dd_rowrange get_m_size(dd_matrixdata<T> *M) {
  dd_rowrange linc = set_card(M->linset);
  return M->rowsize + 1 + linc;
}

template <typename T> inline dd_rowrange get_d_size(dd_matrixdata<T> *M) {
  dd_colrange d;
  if (M->representation == dd_Generator) {
    std::cout << " Generator case\n";
    d = M->colsize + 1;
  } else {
    d = M->colsize;
  }
  return d;
}

template <typename T>
dd_lpdata<T> *dd_CreateLPData_from_M(dd_matrixdata<T> *M) {
  bool localdebug = false;
  dd_colrange d = get_d_size(M);
  if (localdebug)
    std::cout << "dd_CreateLPData_from_M d=" << d << "\n";
  dd_rowrange m = get_m_size(M);
  return dd_CreateLPData<T>(m, d);
}

template <typename T> void dd_CopyArow(T *acopy, T *a, dd_colrange d) {
  dd_colrange j;
  for (j = 0; j < d; j++) {
    acopy[j] = a[j];
  }
}

template <typename T>
void dd_CopyAmatrix(T **Acopy, T **A, dd_rowrange m, dd_colrange d) {
  dd_rowrange i;

  for (i = 0; i < m; i++) {
    dd_CopyArow(Acopy[i], A[i], d);
  }
}

template <typename T> dd_matrixdata<T> *dd_MatrixCopy(dd_matrixdata<T> *M) {
  dd_matrixdata<T> *Mcopy = nullptr;
  dd_rowrange m;
  dd_colrange d;

  m = M->rowsize;
  d = M->colsize;
  if (m >= 0 && d >= 0) {
    Mcopy = dd_CreateMatrix<T>(m, d);
    dd_CopyAmatrix(Mcopy->matrix, M->matrix, m, d);
    dd_CopyArow(Mcopy->rowvec, M->rowvec, d);
    set_copy(Mcopy->linset, M->linset);
    Mcopy->representation = M->representation;
    Mcopy->objective = M->objective;
  }
  return Mcopy;
}

template <typename T> dd_matrixdata<T> *dd_CopyMatrix(dd_matrixdata<T> *M) {
  return dd_MatrixCopy(M);
}

template <typename T>
dd_matrixdata<T> *dd_MatrixNormalizedCopy(dd_matrixdata<T> *M) {
  dd_matrixdata<T> *Mcopy = nullptr;
  dd_rowrange m;
  dd_colrange d;

  m = M->rowsize;
  d = M->colsize;
  if (m >= 0 && d >= 0) {
    Mcopy = dd_CreateMatrix<T>(m, d);
    dd_CopyNormalizedAmatrix(Mcopy->matrix, M->matrix, m, d);
    dd_CopyArow(Mcopy->rowvec, M->rowvec, d);
    set_copy(Mcopy->linset, M->linset);
    Mcopy->representation = M->representation;
    Mcopy->objective = M->objective;
  }
  return Mcopy;
}

template <typename T>
dd_matrixdata<T> *dd_MatrixAppend(dd_matrixdata<T> *M1, dd_matrixdata<T> *M2) {
  dd_matrixdata<T> *M = nullptr;
  dd_rowrange i, m, m1, m2;
  dd_colrange j, d, d1, d2;

  m1 = M1->rowsize;
  d1 = M1->colsize;
  m2 = M2->rowsize;
  d2 = M2->colsize;

  m = m1 + m2;
  d = d1;

  if (d1 >= 0 && d1 == d2 && m1 >= 0 && m2 >= 0) {
    M = dd_CreateMatrix<T>(m, d);
    dd_CopyAmatrix(M->matrix, M1->matrix, m1, d);
    dd_CopyArow(M->rowvec, M1->rowvec, d);
    for (i = 0; i < m1; i++) {
      if (set_member(i + 1, M1->linset))
        set_addelem(M->linset, i + 1);
    }
    for (i = 0; i < m2; i++) {
      for (j = 0; j < d; j++)
        M->matrix[m1 + i][j] = M2->matrix[i][j];
      /* append the second matrix */
      if (set_member(i + 1, M2->linset))
        set_addelem(M->linset, m1 + i + 1);
    }
  }
  return M;
}

void dd_RandomPermutation(dd_rowindex OV, long t, unsigned int seed) {
  long k, j, ovj;
  srand(seed);
  for (j = t; j > 1; j--) {
    k = 1 + random() % t;
    //    std::cerr << "j=" << j << " r=" << r << " u=" << u << " xk=" << xk <<
    //    " k=" << k << "\n"; std::cerr << "j=" << j << " k=" << k << "\n";
    ovj = OV[j];
    OV[j] = OV[k];
    OV[k] = ovj;
  }
}

template <typename T>
bool dd_LexSmaller(T *v1, T *v2,
                   long dmax) { /* dmax is the size of vectors v1,v2 */
  bool determined, smaller;
  dd_colrange j;

  smaller = false;
  determined = false;
  j = 1;
  do {
    if (v1[j - 1] != v2[j - 1]) {
      if (v1[j - 1] < v2[j - 1])
        smaller = true;
      determined = true;
    } else {
      j++;
    }
  } while (!(determined) && (j <= dmax));
  return smaller;
}

template <typename T>
bool dd_LexSmallerFrac(T *v1, T q1, T *v2, T q2,
                       long dmax) { /* dmax is the size of vectors v1,v2 */
  bool determined, smaller;
  dd_colrange j;

  smaller = false;
  determined = false;
  j = 1;
  do {
    if (!dd_EqualFrac(v1[j - 1], q1, v2[j - 1], q2)) {    /* 086 */
      if (dd_SmallerFrac(v1[j - 1], q1, v2[j - 1], q2)) { /*086 */
        smaller = true;
      }
      determined = true;
    } else {
      j++;
    }
  } while (!(determined) && (j <= dmax));
  return smaller;
}

template <typename T> bool dd_LexLarger(T *v1, T *v2, long dmax) {
  return dd_LexSmaller(v2, v1, dmax);
}

template <typename T>
long dd_Partition(dd_rowindex OV, long p, long r, T **A, long dmax) {
  T *x;
  long i, j, ovi;

  x = A[OV[p] - 1];

  i = p - 1;
  j = r + 1;
  while (true) {
    do {
      j--;
    } while (dd_LexLarger(A[OV[j] - 1], x, dmax));
    do {
      i++;
    } while (dd_LexSmaller(A[OV[i] - 1], x, dmax));
    if (i < j) {
      ovi = OV[i];
      OV[i] = OV[j];
      OV[j] = ovi;
    } else {
      return j;
    }
  }
  return -417;
}

template <typename T>
void dd_QuickSort(dd_rowindex OV, long p, long r, T **A, long dmax) {
  long q;

  if (p < r) {
    q = dd_Partition(OV, p, r, A, dmax);
    dd_QuickSort(OV, p, q, A, dmax);
    dd_QuickSort(OV, q + 1, r, A, dmax);
  }
}

template <typename T>
dd_matrixdata<T> *dd_MatrixNormalizedSortedCopy(dd_matrixdata<T> *M,
                                                dd_rowindex *newpos) /* 094 */
{
  /* Sort the rows of Amatrix lexicographically, and return a link to this
  sorted copy. The vector newpos is allocated, where newpos[i] returns the new
  row index of the original row i (i=1,...,M->rowsize). */
  dd_matrixdata<T> *Mcopy = nullptr, Mnorm = nullptr;
  dd_rowrange m, i;
  dd_colrange d;

  /* if (newpos!=nullptr) free(newpos); */
  m = M->rowsize;
  d = M->colsize;
  std::vector<long> roworder(m + 1, 0);
  *newpos = new long[m + 1];
  for (i = 0; i <= m; i++)
    (*newpos)[i] = 0;
  if (m >= 0 && d >= 0) {
    Mnorm = dd_MatrixNormalizedCopy(M);
    Mcopy = dd_CreateMatrix<T>(m, d);
    for (i = 1; i <= m; i++)
      roworder[i] = i;

    dd_RandomPermutation(roworder.data(), m, 123);
    dd_QuickSort(roworder, 1, m, Mnorm->matrix, d);

    dd_PermuteCopyAmatrix(Mcopy->matrix, Mnorm->matrix, m, d, roworder);
    dd_CopyArow(Mcopy->rowvec, M->rowvec, d);
    for (i = 1; i <= m; i++) {
      if (set_member(roworder[i], M->linset))
        set_addelem(Mcopy->linset, i);
      (*newpos)[roworder[i]] = i;
    }
    Mcopy->representation = M->representation;
    Mcopy->objective = M->objective;
    dd_FreeMatrix(Mnorm);
  }
  return Mcopy;
}

template <typename T>
dd_matrixdata<T> *dd_MatrixUniqueCopy(dd_matrixdata<T> *M,
                                      dd_rowindex *newpos) {
  /* Remove row duplicates, and return a link to this sorted copy.
     Linearity rows have priority over the other rows.
     It is better to call this after sorting with dd_MatrixNormalizedSortedCopy.
     The vector newpos is allocated, where *newpos[i] returns the new row index
     of the original row i (i=1,...,M->rowsize).  *newpos[i] is negative if the
     original row is dominated by -*newpos[i] and eliminated in the new copy.
  */
  dd_matrixdata<T> *Mcopy = nullptr;
  dd_rowrange m, i, uniqrows;
  dd_rowset preferredrows;
  dd_colrange d;
  dd_rowindex roworder;

  m = M->rowsize;
  d = M->colsize;
  preferredrows = M->linset;
  roworder = new long[m + 1];
  for (i = 0; i <= m; i++)
    roworder[i] = 0;
  if (m >= 0 && d >= 0) {
    for (i = 1; i <= m; i++)
      roworder[i] = i;
    dd_UniqueRows(roworder, 1, m, M->matrix, d, preferredrows, &uniqrows);

    Mcopy = dd_CreateMatrix<T>(uniqrows, d);
    dd_PermutePartialCopyAmatrix(Mcopy->matrix, M->matrix, m, d, roworder);
    dd_CopyArow(Mcopy->rowvec, M->rowvec, d);
    for (i = 1; i <= m; i++) {
      if (roworder[i] > 0 && set_member(i, M->linset))
        set_addelem(Mcopy->linset, roworder[i]);
    }
    Mcopy->representation = M->representation;
    Mcopy->objective = M->objective;
  }
  *newpos = roworder;
  return Mcopy;
}

template <typename T>
dd_matrixdata<T> *
dd_MatrixNormalizedSortedUniqueCopy(dd_matrixdata<T> *M,
                                    dd_rowindex *newpos) /* 094 */
{
  /* Sort and remove row duplicates, and return a link to this sorted copy.
     Linearity rows have priority over the other rows.
     It is better to call this after sorting with dd_MatrixNormalizedSortedCopy.
     The vector newpos is allocated, where *newpos[i] returns the new row index
     of the original row i (i=1,...,M->rowsize).  *newpos[i] is negative if the
     original row is dominated by -*newpos[i] and eliminated in the new copy.
  */
  dd_matrixdata<T> *M1 = nullptr, M2 = nullptr;
  dd_rowrange m, i;
  dd_colrange d;
  dd_rowindex newpos1 = nullptr, newpos2 = nullptr;

  m = M->rowsize;
  d = M->colsize;
  *newpos = new long[m + 1];
  for (i = 0; i <= m; i++)
    (*newpos)[i] = 0;
  std::vector<long> newpos1r(m + 1, 0);
  if (m >= 0 && d >= 0) {
    M1 = dd_MatrixNormalizedSortedCopy(M, &newpos1);
    for (i = 1; i <= m; i++)
      newpos1r[newpos1[i]] = i; /* reverse of newpos1 */
    M2 = dd_MatrixUniqueCopy(M1, &newpos2);
    set_emptyset(M2->linset);
    for (i = 1; i <= m; i++) {
      if (newpos2[newpos1[i]] > 0) {
        printf("newpos1[%ld]=%ld, newpos2[newpos1[%ld]]=%ld\n", i, newpos1[i],
               i, newpos2[newpos1[i]]);
        if (set_member(i, M->linset))
          set_addelem(M2->linset, newpos2[newpos1[i]]);
        (*newpos)[i] = newpos2[newpos1[i]];
      } else {
        (*newpos)[i] = -newpos1r[-newpos2[newpos1[i]]];
      }
    }
    dd_FreeMatrix(M1);
    delete[] newpos1;
    delete[] newpos2;
  }
  return M2;
}

template <typename T>
dd_matrixdata<T> *dd_MatrixSortedUniqueCopy(dd_matrixdata<T> *M,
                                            dd_rowindex *newpos) /* 094 */
{
  /* Same as dd_MatrixNormalizedSortedUniqueCopy except that it returns a
     unnormalized origial data with original ordering.
  */
  dd_matrixdata<T> *M1 = nullptr, M2 = nullptr;
  dd_rowrange m, i, k, ii;
  dd_colrange d;
  dd_rowindex newpos1 = nullptr, newpos2 = nullptr;

  m = M->rowsize;
  d = M->colsize;
  *newpos = new long[m + 1];
  for (i = 0; i <= m; i++)
    (*newpos)[i] = 0;
  std::vector<long> newpos1r(m + 1, 0);
  if (m >= 0 && d >= 0) {
    M1 = dd_MatrixNormalizedSortedCopy(M, &newpos1);
    for (i = 1; i <= m; i++)
      newpos1r[newpos1[i]] = i; /* reverse of newpos1 */
    M2 = dd_MatrixUniqueCopy(M1, &newpos2);
    dd_FreeMatrix(M1);
    set_emptyset(M2->linset);
    for (i = 1; i <= m; i++) {
      if (newpos2[newpos1[i]] > 0) {
        if (set_member(i, M->linset))
          set_addelem(M2->linset, newpos2[newpos1[i]]);
        (*newpos)[i] = newpos2[newpos1[i]];
      } else {
        (*newpos)[i] = -newpos1r[-newpos2[newpos1[i]]];
      }
    }

    ii = 0;
    set_emptyset(M2->linset);
    for (i = 1; i <= m; i++) {
      k = (*newpos)[i];
      if (k > 0) {
        ii += 1;
        (*newpos)[i] = ii;
        dd_CopyArow(M2->matrix[ii - 1], M->matrix[i - 1], d);
        if (set_member(i, M->linset))
          set_addelem(M2->linset, ii);
      }
    }

    delete[] newpos1;
    delete[] newpos2;
  }
  return M2;
}

template <typename T>
dd_matrixdata<T> *dd_AppendMatrix(dd_matrixdata<T> *M1, dd_matrixdata<T> *M2) {
  return dd_MatrixAppend(M1, M2);
}

template <typename T>
bool dd_MatrixAppendTo(dd_matrixdata<T> **M1, dd_matrixdata<T> *M2) {
  dd_matrixdata<T> *M = nullptr;
  dd_rowrange i, m, m1, m2;
  dd_colrange j, d, d1, d2;
  bool success = false;

  m1 = (*M1)->rowsize;
  d1 = (*M1)->colsize;
  m2 = M2->rowsize;
  d2 = M2->colsize;

  m = m1 + m2;
  d = d1;

  if (d1 >= 0 && d1 == d2 && m1 >= 0 && m2 >= 0) {
    M = dd_CreateMatrix<T>(m, d);
    dd_CopyAmatrix(M->matrix, (*M1)->matrix, m1, d);
    dd_CopyArow(M->rowvec, (*M1)->rowvec, d);
    for (i = 0; i < m1; i++) {
      if (set_member(i + 1, (*M1)->linset))
        set_addelem(M->linset, i + 1);
    }
    for (i = 0; i < m2; i++) {
      for (j = 0; j < d; j++)
        M->matrix[m1 + i][j] = M2->matrix[i][j];
      /* append the second matrix */
      if (set_member(i + 1, M2->linset))
        set_addelem(M->linset, m1 + i + 1);
    }
    dd_FreeMatrix(*M1);
    *M1 = M;
    success = true;
  }
  return success;
}

template <typename T>
void dd_MatrixRowRemove(dd_matrixdata<T> **M, dd_rowrange r) /* 092 */
{
  dd_rowrange i, m;
  m = (*M)->rowsize;

  if (r >= 1 && r <= m) {
    (*M)->rowsize = m - 1;
    dd_FreeArow((*M)->matrix[r - 1]);
    set_delelem((*M)->linset, r);
    /* slide the row headers */
    for (i = r; i < m; i++) {
      (*M)->matrix[i - 1] = (*M)->matrix[i];
      if (set_member(i + 1, (*M)->linset)) {
        set_delelem((*M)->linset, i + 1);
        set_addelem((*M)->linset, i);
      }
    }
  }
}

template <typename T>
void dd_MatrixRowRemove2(dd_matrixdata<T> **M, dd_rowrange r) /* 094 */
{
  dd_rowrange i, m;

  m = (*M)->rowsize;

  if (r >= 1 && r <= m) {
    (*M)->rowsize = m - 1;
    dd_FreeArow((*M)->matrix[r - 1]);
    set_delelem((*M)->linset, r);
    /* slide the row headers */
    for (i = r; i < m; i++) {
      (*M)->matrix[i - 1] = (*M)->matrix[i];
      if (set_member(i + 1, (*M)->linset)) {
        set_delelem((*M)->linset, i + 1);
        set_addelem((*M)->linset, i);
      }
    }
  }
}

template <typename T>
dd_matrixdata<T> *dd_MatrixSubmatrix(dd_matrixdata<T> *M,
                                     dd_rowset delset) /* 092 */
{
  dd_matrixdata<T> *Msub = nullptr;
  dd_rowrange i, isub = 1, m, msub;
  dd_colrange d;

  m = M->rowsize;
  d = M->colsize;
  msub = m;
  if (m >= 0 && d >= 0) {
    for (i = 1; i <= m; i++) {
      if (set_member(i, delset))
        msub -= 1;
    }
    Msub = dd_CreateMatrix<T>(msub, d);
    for (i = 1; i <= m; i++) {
      if (!set_member(i, delset)) {
        dd_CopyArow(Msub->matrix[isub - 1], M->matrix[i - 1], d);
        if (set_member(i, M->linset)) {
          set_addelem(Msub->linset, isub);
        }
        isub++;
      }
    }
    dd_CopyArow(Msub->rowvec, M->rowvec, d);
    Msub->representation = M->representation;
    Msub->objective = M->objective;
  }
  return Msub;
}

template <typename T>
dd_matrixdata<T> *dd_MatrixSubmatrix2(dd_matrixdata<T> *M, dd_rowset delset,
                                      dd_rowindex *newpos) /* 092 */
{ /* returns a pointer to a new matrix which is a submatrix of M with rows in
  delset removed.  *newpos[i] returns the position of the original row i in the
  new matrix. It is -1 if and only if it is deleted.
  */

  dd_matrixdata<T> *Msub = nullptr;
  dd_rowrange i, isub = 1, m, msub;
  dd_colrange d;
  dd_rowindex roworder;

  m = M->rowsize;
  d = M->colsize;
  msub = m;
  if (m >= 0 && d >= 0) {
    roworder = new long[m + 1];
    for (i = 0; i <= m; i++)
      roworder[i] = 0;
    for (i = 1; i <= m; i++) {
      if (set_member(i, delset))
        msub -= 1;
    }
    Msub = dd_CreateMatrix<T>(msub, d);
    for (i = 1; i <= m; i++) {
      if (set_member(i, delset)) {
        roworder[i] = 0; /* zero means the row i is removed */
      } else {
        dd_CopyArow(Msub->matrix[isub - 1], M->matrix[i - 1], d);
        if (set_member(i, M->linset)) {
          set_addelem(Msub->linset, isub);
        }
        roworder[i] = isub;
        isub++;
      }
    }
    *newpos = roworder;
    dd_CopyArow(Msub->rowvec, M->rowvec, d);
    Msub->representation = M->representation;
    Msub->objective = M->objective;
  }
  return Msub;
}

template <typename T>
dd_matrixdata<T> *dd_MatrixSubmatrix2L(dd_matrixdata<T> *M, dd_rowset delset,
                                       dd_rowindex *newpos) /* 094 */
{ /* This is same as dd_MatrixSubmatrix2 except that the linearity rows will be
     shifted up so that they are at the top of the matrix.
  */
  dd_matrixdata<T> *Msub = nullptr;
  dd_rowrange i, iL, iI, m, msub;
  dd_colrange d;
  dd_rowindex roworder;

  m = M->rowsize;
  d = M->colsize;
  msub = m;
  if (m >= 0 && d >= 0) {
    roworder = new long[m + 1];
    for (i = 0; i <= m; i++)
      roworder[i] = 0;
    for (i = 1; i <= m; i++) {
      if (set_member(i, delset))
        msub -= 1;
    }
    Msub = dd_CreateMatrix<T>(msub, d);
    iL = 1;
    iI = set_card(M->linset) + 1; /* starting positions */
    for (i = 1; i <= m; i++) {
      if (set_member(i, delset)) {
        roworder[i] = 0; /* zero means the row i is removed */
      } else {
        if (set_member(i, M->linset)) {
          dd_CopyArow(Msub->matrix[iL - 1], M->matrix[i - 1], d);
          set_delelem(Msub->linset, i);
          set_addelem(Msub->linset, iL);
          roworder[i] = iL;
          iL += 1;
        } else {
          dd_CopyArow(Msub->matrix[iI - 1], M->matrix[i - 1], d);
          roworder[i] = iI;
          iI += 1;
        }
      }
    }
    *newpos = roworder;
    dd_CopyArow(Msub->rowvec, M->rowvec, d);
    Msub->representation = M->representation;
    Msub->objective = M->objective;
  }
  return Msub;
}

template <typename T>
void dd_MatrixRowsRemove(dd_matrixdata<T> **M, dd_rowset delset) /* 094 */
{
  dd_matrixdata<T> *Msub = nullptr;

  Msub = dd_MatrixSubmatrix(*M, delset);
  dd_FreeMatrix(*M);
  *M = Msub;
}

template <typename T>
void dd_MatrixRowsRemove2(dd_matrixdata<T> **M, dd_rowset delset,
                          dd_rowindex *newpos) /* 094 */
{
  dd_matrixdata<T> *Msub = nullptr;

  Msub = dd_MatrixSubmatrix2(*M, delset, newpos);
  dd_FreeMatrix(*M);
  *M = Msub;
}

template <typename T>
void dd_MatrixShiftupLinearity(dd_matrixdata<T> **M,
                               dd_rowindex *newpos) /* 094 */
{
  dd_matrixdata<T> *Msub = nullptr;
  dd_rowset delset;

  set_initialize(&delset, (*M)->rowsize); /* emptyset */
  Msub = dd_MatrixSubmatrix2L(*M, delset, newpos);
  dd_FreeMatrix(*M);
  *M = Msub;

  delete[] delset;
}

template <typename T> void dd_SetLinearity(dd_matrixdata<T> *M, char *line) {
  int i = 0;
  dd_rowrange eqsize, var;
  char *next;
  const char ct[] = ", "; /* allows separators "," and " ". */

  next = strtok(line, ct);
  eqsize = atol(next);
  while (i < eqsize && (next = strtok(nullptr, ct)) != nullptr) {
    var = atol(next);
    set_addelem(M->linset, var);
    i++;
  }
  if (i != eqsize) {
    std::cout << "* Warning: there are inconsistencies in linearity setting.\n";
  }
  return;
}

template <typename T>
dd_polyhedradata<T> *dd_DDMatrix2Poly(dd_matrixdata<T> *M, dd_ErrorType *err, size_t const& maxiter, std::ostream& os) {
  dd_rowrange i;
  dd_colrange j;

  *err = dd_NoError;
  if (M->rowsize < 0 || M->colsize < 0) {
    *err = dd_NegativeMatrixSize;
    return nullptr;
  }
  dd_polyhedradata<T> *poly = dd_CreatePolyhedraData<T>(M->rowsize, M->colsize);
  poly->representation = M->representation;
  poly->homogeneous = true;

  for (i = 1; i <= M->rowsize; i++) {
    if (set_member(i, M->linset))
      poly->EqualityIndex[i] = 1;
    for (j = 1; j <= M->colsize; j++) {
      poly->A[i - 1][j - 1] = M->matrix[i - 1][j - 1];
      if (j == 1 && M->matrix[i - 1][j - 1] != 0)
        poly->homogeneous = false;
    }
  }
  dd_DoubleDescription(poly, err, maxiter, os);
  return poly;
}

template <typename T>
dd_polyhedradata<T> *dd_DDMatrix2Poly2(dd_matrixdata<T> *M,
                                       dd_RowOrderType horder,
                                       dd_ErrorType *err, size_t const& maxiter, std::ostream& os) {
  dd_rowrange i;
  dd_colrange j;

  *err = dd_NoError;
  if (M->rowsize < 0 || M->colsize < 0) {
    *err = dd_NegativeMatrixSize;
    return nullptr;
  }
  dd_polyhedradata<T> *poly = dd_CreatePolyhedraData<T>(M->rowsize, M->colsize);
  poly->representation = M->representation;
  poly->homogeneous = true;

  for (i = 1; i <= M->rowsize; i++) {
    if (set_member(i, M->linset))
      poly->EqualityIndex[i] = 1;
    for (j = 1; j <= M->colsize; j++) {
      poly->A[i - 1][j - 1] = M->matrix[i - 1][j - 1];
      if (j == 1 && M->matrix[i - 1][j - 1] != 0)
        poly->homogeneous = false;
    }
  }
  dd_DoubleDescription2(poly, horder, err, maxiter, os);
  return poly;
}

template <typename T>
void dd_CopyRay(T *a, dd_colrange d_origsize, dd_raydata<T> *RR,
                dd_colindex reducedcol) {
  long j, j1;
  for (j = 1; j <= d_origsize; j++) {
    j1 = reducedcol[j];
    if (j1 > 0) {
      a[j - 1] = RR->Ray[j1 - 1];
      /* the original column j is mapped to j1, and thus
         copy the corresponding component */
    } else {
      a[j - 1] = 0;
      /* original column is redundant and removed for computation */
    }
  }
}

template <typename T> void dd_ComputeAinc(dd_polyhedradata<T> *poly) {
  /* This generates the input incidence array poly->Ainc, and
     two sets: poly->Ared, poly->Adom.
  */
  dd_bigrange k;
  dd_rowrange i, m1;
  dd_colrange j;
  bool redundant;
  dd_matrixdata<T> *M = nullptr;
  T sum, temp;

  if (poly->AincGenerated == true)
    return;

  M = dd_CopyOutput(poly);
  poly->n = M->rowsize;
  m1 = poly->m1;
  /* this number is same as poly->m, except when
     poly is given by nonhomogeneous inequalty:
     !(poly->homogeneous) && poly->representation==Inequality,
     it is poly->m+1.   See dd_ConeDataLoad.
  */
  poly->Ainc = new set_type[m1];
  for (i = 1; i <= m1; i++)
    set_initialize(&(poly->Ainc[i - 1]), poly->n);
  set_initialize(&(poly->Ared), m1);
  set_initialize(&(poly->Adom), m1);

  for (k = 1; k <= poly->n; k++) {
    for (i = 1; i <= poly->m; i++) {
      sum = 0;
      for (j = 1; j <= poly->d; j++)
        sum += poly->A[i - 1][j - 1] * M->matrix[k - 1][j - 1];
      if (sum == 0) {
        set_addelem(poly->Ainc[i - 1], k);
      }
    }
    if (!(poly->homogeneous) && poly->representation == dd_Inequality) {
      if (M->matrix[k - 1][0] == 0) {
        set_addelem(poly->Ainc[m1 - 1],
                    k); /* added infinity inequality (1,0,0,...,0) */
      }
    }
  }

  for (i = 1; i <= m1; i++) {
    if (set_card(poly->Ainc[i - 1]) == M->rowsize) {
      set_addelem(poly->Adom, i);
    }
  }
  for (i = m1; i >= 1; i--) {
    if (set_card(poly->Ainc[i - 1]) == 0) {
      redundant = true;
      set_addelem(poly->Ared, i);
    } else {
      redundant = false;
      for (k = 1; k <= m1; k++) {
        if (k != i && !set_member(k, poly->Ared) &&
            !set_member(k, poly->Adom) &&
            set_subset(poly->Ainc[i - 1], poly->Ainc[k - 1])) {
          if (!redundant) {
            redundant = true;
          }
          set_addelem(poly->Ared, i);
        }
      }
    }
  }
  dd_FreeMatrix(M);
  poly->AincGenerated = true;
}

template <typename T>
bool dd_InputAdjacentQ(set_type &common, long &lastn, dd_polyhedradata<T> *poly,
                       dd_rowrange i1, dd_rowrange i2)
/* Before calling this function, RedundantSet must be
   a set of row indices whose removal results in a minimal
   nonredundant system to represent the input polyhedron,
   DominantSet must be the set of row indices which are
   active at every extreme points/rays.
*/
{
  dd_rowrange i;

  if (poly->AincGenerated == false)
    dd_ComputeAinc<T>(poly);
  if (lastn != poly->n) {
    if (lastn > 0)
      set_free(common);
    set_initialize(&common, poly->n);
    lastn = poly->n;
  }
  if (set_member(i1, poly->Ared) || set_member(i2, poly->Ared))
    return false;
  // dominant inequality is considered adjacencent to all others.
  if (set_member(i1, poly->Adom) || set_member(i2, poly->Adom))
    return true;
  set_int(common, poly->Ainc[i1 - 1], poly->Ainc[i2 - 1]);
  i = 0;
  while (i < poly->m1) {
    i++;
    if (i != i1 && i != i2 && !set_member(i, poly->Ared) &&
        !set_member(i, poly->Adom) && set_subset(common, poly->Ainc[i - 1]))
      return false;
  }
  return true;
}

template <typename T>
dd_setfamily *dd_CopyIncidence(dd_polyhedradata<T> *poly) {
  dd_bigrange k;
  dd_rowrange i;

  if (poly->child == nullptr || poly->child->CompStatus != dd_AllFound)
    return nullptr;
  if (poly->AincGenerated == false)
    dd_ComputeAinc(poly);
  dd_setfamily *F = dd_CreateSetFamily(poly->n, poly->m1);
  for (i = 1; i <= poly->m1; i++)
    for (k = 1; k <= poly->n; k++)
      if (set_member(k, poly->Ainc[i - 1]))
        set_addelem(F->set[k - 1], i);
  return F;
}

template <typename T>
dd_setfamily *dd_CopyInputIncidence(dd_polyhedradata<T> *poly) {
  dd_rowrange i;

  if (poly->child == nullptr || poly->child->CompStatus != dd_AllFound)
    return nullptr;
  if (poly->AincGenerated == false)
    dd_ComputeAinc(poly);
  dd_setfamily *F = dd_CreateSetFamily(poly->m1, poly->n);
  for (i = 0; i < poly->m1; i++)
    set_copy(F->set[i], poly->Ainc[i]);
  return F;
}

template <typename T>
dd_setfamily *dd_CopyAdjacency(dd_polyhedradata<T> *poly) {
  dd_raydata<T> *RayPtr1;
  dd_raydata<T> *RayPtr2;
  dd_setfamily *F = nullptr;
  long pos1, pos2;
  dd_bigrange lstart, k, n;
  set_type linset, allset;
  bool adj;

  if (poly->n == 0 && poly->homogeneous &&
      poly->representation == dd_Inequality) {
    n = 1; /* the origin (the unique vertex) should be output. */
  } else {
    n = poly->n;
  }
  set_initialize(&linset, n);
  set_initialize(&allset, n);
  if (poly->child == nullptr || poly->child->CompStatus != dd_AllFound)
    goto _L99;
  F = dd_CreateSetFamily(n, n);
  if (n <= 0)
    goto _L99;
  poly->child->LastRay->Next = nullptr;
  for (RayPtr1 = poly->child->FirstRay, pos1 = 1; RayPtr1 != nullptr;
       RayPtr1 = RayPtr1->Next, pos1++) {
    for (RayPtr2 = poly->child->FirstRay, pos2 = 1; RayPtr2 != nullptr;
         RayPtr2 = RayPtr2->Next, pos2++) {
      if (RayPtr1 != RayPtr2) {
        dd_CheckAdjacency(poly->child, &RayPtr1, &RayPtr2, &adj);
        if (adj) {
          set_addelem(F->set[pos1 - 1], pos2);
        }
      }
    }
  }
  lstart = poly->n - poly->ldim + 1;
  set_compl(allset, allset); /* allset is set to the ground set. */
  for (k = lstart; k <= poly->n; k++) {
    set_addelem(linset, k); /* linearity set */
    set_copy(F->set[k - 1],
             allset); /* linearity generator is adjacent to all */
  }
  for (k = 1; k < lstart; k++) {
    set_uni(F->set[k - 1], F->set[k - 1], linset);
    /* every generator is adjacent to all linearity generators */
  }
_L99:;
  set_free(allset);
  set_free(linset);
  return F;
}

template <typename T>
dd_setfamily *dd_CopyInputAdjacency(dd_polyhedradata<T> *poly) {
  dd_rowrange i, j;
  set_type common;
  long lastn = 0;
  if (poly->child == nullptr || poly->child->CompStatus != dd_AllFound)
    return nullptr;
  if (poly->AincGenerated == false)
    dd_ComputeAinc(poly);
  dd_setfamily *F = dd_CreateSetFamily(poly->m1, poly->m1);
  for (i = 1; i <= poly->m1; i++)
    for (j = 1; j <= poly->m1; j++)
      if (i != j && dd_InputAdjacentQ<T>(common, lastn, poly, i, j))
        set_addelem(F->set[i - 1], j);
  return F;
}

template <typename T>
dd_matrixdata<T> *dd_CopyOutput(dd_polyhedradata<T> *poly) {
  dd_raydata<T> *RayPtr;
  dd_rowrange i = 0;
  dd_rowrange total;
  dd_colrange j, j1;
  bool outputorigin = false;

  total = poly->child->LinearityDim + poly->child->FeasibleRayCount;

  if (poly->child->d <= 0 || poly->child->newcol[1] == 0)
    total = total - 1;
  if (total == 0 && poly->homogeneous &&
      poly->representation == dd_Inequality) {
    total = 1;
    outputorigin = true;
    // the origin (the unique vertex) should be output.
  }
  if (poly->child == nullptr || poly->child->CompStatus != dd_AllFound)
    return nullptr;

  dd_matrixdata<T> *M = dd_CreateMatrix<T>(total, poly->d);
  RayPtr = poly->child->FirstRay;
  while (RayPtr != nullptr) {
    if (RayPtr->feasible) {
      dd_CopyRay(M->matrix[i], poly->d, RayPtr, poly->child->newcol);
      i++;
    }
    RayPtr = RayPtr->Next;
  }
  for (j = 2; j <= poly->d; j++) {
    if (poly->child->newcol[j] == 0) {
      for (j1 = 0; j1 < poly->d; j1++)
        M->matrix[i][j1] = poly->child->Bsave[j1][j - 1];
      set_addelem(M->linset, i + 1);
      i++;
    }
  }
  if (outputorigin) {
    // output the origin for homogeneous H-polyhedron with no rays.
    M->matrix[0][0] = 1;
    for (j = 1; j < poly->d; j++)
      M->matrix[0][j] = 0;
  }
  //  dd_MatrixIntegerFilter(M);
  if (poly->representation == dd_Inequality)
    M->representation = dd_Generator;
  else
    M->representation = dd_Inequality;
  return M;
}

template <typename T>
dd_matrixdata<T> *dd_CopyInput(dd_polyhedradata<T> *poly) {
  dd_matrixdata<T> *M = nullptr;
  dd_rowrange i;

  M = dd_CreateMatrix<T>(poly->m, poly->d);
  dd_CopyAmatrix(M->matrix, poly->A, poly->m, poly->d);
  for (i = 1; i <= poly->m; i++)
    if (poly->EqualityIndex[i] == 1)
      set_addelem(M->linset, i);
  // dd_MatrixIntegerFilter(M);
  if (poly->representation == dd_Generator)
    M->representation = dd_Generator;
  else
    M->representation = dd_Inequality;
  return M;
}

template <typename T>
dd_matrixdata<T> *dd_CopyGenerators(dd_polyhedradata<T> *poly) {
  dd_matrixdata<T> *M = nullptr;

  if (poly->representation == dd_Generator)
    M = dd_CopyInput(poly);
  else
    M = dd_CopyOutput(poly);
  return M;
}

template <typename T>
dd_matrixdata<T> *dd_CopyInequalities(dd_polyhedradata<T> *poly) {
  dd_matrixdata<T> *M = nullptr;

  if (poly->representation == dd_Inequality)
    M = dd_CopyInput(poly);
  else
    M = dd_CopyOutput(poly);
  return M;
}

template <typename T>
dd_matrixdata<T> *dd_BlockElimination(dd_matrixdata<T> *M, dd_colset delset,
                                      dd_ErrorType *error)
/* Eliminate the variables (columns) delset by
   the Block Elimination with dd_DoubleDescription algorithm.

   Given (where y is to be eliminated):
   c1 + A1 x + B1 y >= 0
   c2 + A2 x + B2 y =  0

   1. First construct the dual system:  z1^T B1 + z2^T B2 = 0, z1 >= 0.
   2. Compute the generators of the dual.
   3. Then take the linear combination of the original system with each
   generator.
   4. Remove redundant inequalies.

*/
{
  dd_matrixdata<T> *Mdual = nullptr, Mproj = nullptr, Gdual = nullptr;
  dd_rowrange i, h, m, mproj, mdual, linsize;
  dd_colrange j, k, d, dproj, ddual, delsize;
  T temp, prod;
  dd_polyhedradata<T> *dualpoly;
  dd_ErrorType err = dd_NoError;
  bool localdebug = false;

  *error = dd_NoError;
  m = M->rowsize;
  d = M->colsize;
  std::vector<long> delindex(d);
  k = 0;
  delsize = 0;
  for (j = 1; j <= d; j++) {
    if (set_member(j, delset)) {
      delsize++;
      delindex[k] = j; /* stores the kth deletion column index */
      k++;
    }
  }

  linsize = set_card(M->linset);
  ddual = m + 1;
  mdual = delsize + m - linsize; /* #equalitions + dimension of z1 */

  /* setup the dual matrix */
  Mdual = dd_CreateMatrix<T>(mdual, ddual);
  Mdual->representation = dd_Inequality;
  for (i = 1; i <= delsize; i++) {
    set_addelem(Mdual->linset, i); /* equality */
    for (j = 1; j <= m; j++)
      Mdual->matrix[i - 1][j] = M->matrix[j - 1][delindex[i - 1] - 1];
  }

  k = 0;
  for (i = 1; i <= m; i++) {
    if (!set_member(i, M->linset)) {
      /* set nonnegativity for the dual variable associated with
         each non-linearity inequality. */
      k++;
      Mdual->matrix[delsize + k - 1][i] = 1;
    }
  }

  /* 2. Compute the generators of the dual system. */
  dualpoly = dd_DDMatrix2Poly(Mdual, &err);
  Gdual = dd_CopyGenerators(dualpoly);

  /* 3. Take the linear combination of the original system with each generator.
   */
  dproj = d - delsize;
  mproj = Gdual->rowsize;
  Mproj = dd_CreateMatrix<T>(mproj, dproj);
  Mproj->representation = dd_Inequality;
  set_copy(Mproj->linset, Gdual->linset);

  for (i = 0; i < mproj; i++) {
    k = 0;
    for (j = 1; j <= d; j++) {
      if (!set_member(j, delset)) {
        prod = 0;
        for (h = 1; h <= m; h++)
          prod += M->matrix[h - 1][j - 1] * Gdual->matrix[i][h];
        Mproj->matrix[i][k] = prod;
        k++; /* new index of the variable x_j  */
      }
    }
  }
  if (localdebug)
    printf("Size of the projection system: %ld x %ld\n", mproj, dproj);

  dd_FreePolyhedra(dualpoly);
  dd_FreeMatrix(Mdual);
  dd_FreeMatrix(Gdual);
  return Mproj;
}

template <typename T>
dd_matrixdata<T> *dd_FourierElimination(dd_matrixdata<T> *M,
                                        dd_ErrorType *error)
/* Eliminate the last variable (column) from the given H-matrix using
   the standard Fourier Elimination.
 */
{
  dd_rowrange i, inew, ip, in, iz, m, mpos = 0, mneg = 0, mzero = 0, mnew;
  dd_colrange j, d, dnew;
  T temp1, temp2;
  bool localdebug = false;

  *error = dd_NoError;
  m = M->rowsize;
  d = M->colsize;
  if (d <= 1) {
    *error = dd_ColIndexOutOfRange;
    if (localdebug)
      printf(
          "The number of column is too small: %ld for Fourier's Elimination.\n",
          d);
    return nullptr;
  }

  if (M->representation == dd_Generator) {
    *error = dd_NotAvailForV;
    if (localdebug)
      printf("Fourier's Elimination cannot be applied to a V-polyhedron.\n");
    return nullptr;
  }

  if (set_card(M->linset) > 0) {
    *error = dd_CannotHandleLinearity;
    if (localdebug)
      printf("The Fourier Elimination function does not handle equality in "
             "this version.\n");
    return nullptr;
  }

  /* Create temporary spaces to be removed at the end of this function */
  std::vector<long> posrowindex(m, 0);
  std::vector<long> negrowindex(m, 0);
  std::vector<long> zerorowindex(m, 0);

  for (i = 0; i < m; i++) {
    if (M->matrix[i][d - 1] > 0) {
      posrowindex[mpos] = i;
      mpos++;
    } else {
      if (M->matrix[i][d - 1] < 0) {
        negrowindex[mneg] = i;
        mneg++;
      } else {
        zerorowindex[mzero] = i;
        mzero++;
      }
    }
  }

  /* The present code generates so many redundant inequalities and thus
     is quite useless, except for very small examples
  */
  mnew = mzero + mpos * mneg; /* the total number of rows after elimination */
  dnew = d - 1;

  dd_matrixdata<T> *Mnew = dd_CreateMatrix<T>(mnew, dnew);
  dd_CopyArow(Mnew->rowvec, M->rowvec, dnew);
  Mnew->representation = M->representation;
  Mnew->objective = M->objective;

  /* Copy the inequalities independent of x_d to the top of the new matrix. */
  for (iz = 0; iz < mzero; iz++)
    for (j = 0; j < dnew; j++)
      Mnew->matrix[iz][j] = M->matrix[zerorowindex[iz]][j];

  /* Create the new inequalities by combining x_d positive and negative ones. */
  inew = mzero; /* the index of the last x_d zero inequality */
  for (ip = 0; ip < mpos; ip++) {
    for (in = 0; in < mneg; in++) {
      temp1 = -M->matrix[negrowindex[in]][d - 1];
      temp2 = M->matrix[posrowindex[ip]][d - 1];
      for (j = 0; j < dnew; j++)
        Mnew->matrix[inew][j] = M->matrix[posrowindex[ip]][j] * temp1 +
                                M->matrix[negrowindex[in]][j] * temp2;
      inew++;
    }
  }

  return Mnew;
}

template <typename T> dd_lpdata<T> *dd_Matrix2LP(dd_matrixdata<T> *M) {
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_lpdata<T> *lp;

  linc = set_card(M->linset);
  m = M->rowsize + 1 + linc;
  /* We represent each equation by two inequalities.
     This is not the best way but makes the code simple. */
  d = M->colsize;

  lp = dd_CreateLPData<T>(m, d);
  lp->objective = M->objective;

  irev = M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    if (set_member(i, M->linset)) {
      irev++;
      set_addelem(lp->equalityset, i); /* it is equality. */
      /* the reversed row irev is not in the equality set. */
      for (j = 0; j < M->colsize; j++)
        lp->A[irev - 1][j] = -M->matrix[i - 1][j];
    }
    for (j = 0; j < M->colsize; j++)
      lp->A[i - 1][j] = M->matrix[i - 1][j];
  }
  for (j = 0; j < M->colsize; j++)
    lp->A[m - 1][j] = M->rowvec[j]; /* objective row */
  return lp;
}

template <typename T>
dd_lpdata<T> *dd_Matrix2Feasibility(dd_matrixdata<T> *M, dd_ErrorType *err)
/* Load a matrix to create an LP object for feasibility.   It is
   essentially the dd_Matrix2LP except that the objject function
   is set to identically ZERO (maximization).

*/
/*  094 */
{
  dd_rowrange m, linc;
  dd_colrange j;
  dd_lpdata<T> *lp;

  *err = dd_NoError;
  linc = set_card(M->linset);
  m = M->rowsize + 1 + linc;
  /* We represent each equation by two inequalities.
     This is not the best way but makes the code simple. */

  lp = dd_Matrix2LP(M);
  lp->objective =
      dd_LPmax; /* since the objective is zero, this is not important */
  for (j = 1; j <= M->colsize; j++) {
    lp->A[m - 1][j - 1] = 0; /* set the objective to zero. */
  } /*of j*/

  return lp;
}

template <typename T>
dd_lpdata<T> *dd_Matrix2Feasibility2(dd_matrixdata<T> *M, dd_rowset R,
                                     dd_rowset S, dd_ErrorType *err)
/* Load a matrix to create an LP object for feasibility with additional equality
   and strict inequality constraints given by R and S.  There are three types of
   inequalities:

   b_r + A_r x =  0     Linearity (Equations) specified by M
   b_s + A_s x >  0     Strict Inequalities specified by row index set S
   b_t + A_t x >= 0     The rest inequalities in M

   Where the linearity is considered here as the union of linearity specified by
   M and the additional set R.  When S contains any linearity rows, those
   rows are considered linearity (equation).  Thus S does not overlide
   linearity. To find a feasible solution, we set an LP

   maximize  z
   subject to
   b_r + A_r x     =  0      all r in Linearity
   b_s + A_s x - z >= 0      for all s in S
   b_t + A_t x     >= 0      for all the rest rows t
   1           - z >= 0      to make the LP bounded.

   Clearly, the feasibility problem has a solution iff the LP has a positive
   optimal value. The variable z will be the last variable x_{d+1}.

*/
/*  094 */
{
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_lpdata<T> *lp;
  dd_rowset L;

  *err = dd_NoError;
  set_initialize(&L, M->rowsize);
  set_uni(L, M->linset, R);
  linc = set_card(L);
  m = M->rowsize + 1 + linc + 1;
  /* We represent each equation by two inequalities.
     This is not the best way but makes the code simple. */
  d = M->colsize + 1;

  lp = dd_CreateLPData<T>(m, d);
  lp->objective = dd_LPmax;

  irev = M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    if (set_member(i, L)) {
      irev++;
      set_addelem(lp->equalityset, i); /* it is equality. */
      /* the reversed row irev is not in the equality set. */
      for (j = 1; j <= M->colsize; j++)
        lp->A[irev - 1][j - 1] = -M->matrix[i - 1][j - 1];
    } else {
      if (set_member(i, S)) {
        lp->A[i - 1][M->colsize] = -1;
      }
    }
    for (j = 1; j <= M->colsize; j++) {
      lp->A[i - 1][j - 1] = M->matrix[i - 1][j - 1];
    } /*of j*/
  } /*of i*/
  for (j = 1; j <= d; j++) {
    lp->A[m - 2][j - 1] = 0; /* initialize */
  } /*of j*/
  lp->A[m - 2][0] = 1;           /* the bounding constraint. */
  lp->A[m - 2][M->colsize] = -1; /* the bounding constraint. */
  for (j = 1; j <= d; j++) {
    lp->A[m - 1][j - 1] = 0; /* initialize */
  } /*of j*/
  lp->A[m - 1][M->colsize] = 1;

  set_free(L);
  return lp;
}

template <typename T> void dd_FreeLPData(dd_lpdata<T> *lp) {
  if ((lp) != nullptr) {
    dd_FreeArow(lp->dsol);
    dd_FreeArow(lp->sol);
    dd_FreeBmatrix(lp->d_alloc, lp->B);
    dd_FreeAmatrix(lp->m_alloc, lp->A);
    set_free(lp->equalityset);
    set_free(lp->redset_extra);
    set_free(lp->redset_accum);
    set_free(lp->posset_extra);
    delete[] lp->nbindex;
    delete[] lp->given_nbindex;
    delete lp;
  }
}

template <typename T> bool dd_LPReverseRow(dd_lpdata<T> *lp, dd_rowrange i) {
  dd_colrange j;
  bool success = false;

  if (i >= 1 && i <= lp->m) {
    lp->LPS = dd_LPSundecided;
    for (j = 1; j <= lp->d; j++) {
      lp->A[i - 1][j - 1] = -lp->A[i - 1][j - 1];
      /* negating the i-th constraint of A */
    }
    success = true;
  }
  return success;
}

template <typename T>
bool dd_LPReplaceRow(dd_lpdata<T> *lp, dd_rowrange i, T *a) {
  dd_colrange j;
  bool success = false;

  if (i >= 1 && i <= lp->m) {
    lp->LPS = dd_LPSundecided;
    for (j = 1; j <= lp->d; j++) {
      lp->A[i - 1][j - 1] = a[j - 1];
      /* replacing the i-th constraint by a */
    }
    success = true;
  }
  return success;
}

template <typename T>
void dd_TableauEntry(T &x, dd_colrange d_size, T **X, T **Ts, dd_rowrange r,
                     dd_colrange s)
/* Compute the (r,s) entry of X.T   */
{
  dd_colrange j;
  x = 0;
  for (j = 0; j < d_size; j++)
    x += X[r - 1][j] * Ts[j][s - 1];
}

template <typename T>
void dd_GetRedundancyInformation(dd_rowrange m_size, dd_colrange d_size, T **A,
                                 T **Ts, dd_rowindex bflag, dd_rowset redset)
/* Some basic variables that are forced to be nonnegative will be output.  These
   are variables whose dictionary row components are all nonnegative.   */
{
  dd_colrange j;
  dd_rowrange i;
  T x;
  bool red = false;

  for (i = 1; i <= m_size; i++) {
    red = true;
    for (j = 1; j <= d_size; j++) {
      dd_TableauEntry(x, d_size, A, Ts, i, j);
      if (red && x < 0)
        red = false;
    }
    if (bflag[i] < 0 && red) {
      set_addelem(redset, i);
    }
  }
}

template <typename T>
void dd_SelectDualSimplexPivot(dd_rowrange m_size, dd_colrange d_size,
                               int Phase1, T **A, T **Ts,
                               dd_colindex nbindex_ref, dd_rowindex bflag,
                               dd_rowrange objrow, dd_colrange rhscol,
                               bool lexicopivot, dd_rowrange *r, dd_colrange *s,
                               bool *selected, dd_LPStatusType *lps,
                               data_temp_simplex<T> *data) {
  /* selects a dual simplex pivot (*r,*s) if the current
     basis is dual feasible and not optimal. If not dual feasible,
     the procedure returns *selected=false and *lps=LPSundecided.
     If Phase1=true, the RHS column will be considered as the negative
     of the column of the largest variable (==m_size).  For this case, it is
     assumed that the caller used the auxiliary row (with variable m_size) to
     make the current dictionary dual feasible before calling this routine so
     that the nonbasic column for m_size corresponds to the auxiliary variable.
  */
  bool colselected = false, rowselected = false, dualfeasible = true;
  dd_rowrange i, iref;
  dd_colrange j, k;
  T val, valn, minval(0), rat, minrat(0), minrat_q(1);
  //  T* rcost;
  //  dd_colset tieset;
  //  dd_colset stieset;  /* store the column indices with tie */
  //  rcost=new T[d_size];
  //  set_initialize(&tieset,d_size);
  //  set_initialize(&stieset,d_size);

  *r = 0;
  *s = 0;
  *selected = false;
  *lps = dd_LPSundecided;
  for (j = 1; j <= d_size; j++) {
    if (j != rhscol) {
      dd_TableauEntry(data->rcost[j - 1], d_size, A, Ts, objrow, j);
      if (data->rcost[j - 1] > 0) {
        dualfeasible = false;
      }
    }
  }
  if (dualfeasible) {
    while ((*lps == dd_LPSundecided) && (!rowselected) && (!colselected)) {
      for (i = 1; i <= m_size; i++) {
        if (i != objrow && bflag[i] == -1) { /* i is a basic variable */
          if (Phase1) {
            dd_TableauEntry(val, d_size, A, Ts, i, bflag[m_size]);
            val = -val;
            /* for dual Phase I.  The RHS (dual objective) is the negative of
             * the auxiliary variable column. */
          } else {
            dd_TableauEntry(val, d_size, A, Ts, i, rhscol);
          }
          if (val < minval) {
            *r = i;
            minval = val;
          }
        }
      }
      if (minval >= 0) {
        *lps = dd_Optimal;
      } else {
        rowselected = true;
        set_emptyset(data->tieset);
        for (j = 1; j <= d_size; j++) {
          dd_TableauEntry(val, d_size, A, Ts, *r, j);
          if (j != rhscol && val > 0) {
            bool is_field = true;
            if (is_field) {
              rat = -data->rcost[j - 1] / val;
              if (*s == 0 || rat < minrat) {
                minrat = rat;
                *s = j;
                set_emptyset(data->tieset);
                set_addelem(data->tieset, j);
              } else {
                if (rat == minrat) {
                  set_addelem(data->tieset, j);
                }
              }
            } else {
              // ring case
              rat = -data->rcost[j - 1];
              if (*s == 0 || dd_SmallerFrac(rat, val, minrat, minrat_q)) {
                minrat = rat;
                minrat_q = val;
                *s = j;
                set_emptyset(data->tieset);
                set_addelem(data->tieset, j);
              } else {
                if (dd_EqualFrac(rat, val, minrat, minrat_q)) {
                  set_addelem(data->tieset, j);
                }
              }
            }
          }
        }
        if (*s > 0) {
          if (!lexicopivot || set_card(data->tieset) == 1) {
            colselected = true;
            *selected = true;
          } else { /* lexicographic rule with respect to the given reference
                      cobasis.  */
            *s = 0;
            k = 2; /* k runs through the column indices except RHS. */
            do {
              iref = nbindex_ref[k]; /* iref runs though the reference basic
                                        indices */
              if (iref > 0) {
                j = bflag[iref];
                if (j > 0) {
                  if (set_member(j, data->tieset) &&
                      set_card(data->tieset) == 1) {
                    *s = j;
                    colselected = true;
                  } else {
                    set_delelem(data->tieset, j);
                    /* iref is cobasic, and the corresponding col is not the
                     * pivot column except it is the last one. */
                  }
                } else {
                  *s = 0;
                  for (j = 1; j <= d_size; j++) {
                    if (set_member(j, data->tieset)) {
                      dd_TableauEntry(val, d_size, A, Ts, *r, j);
                      dd_TableauEntry(valn, d_size, A, Ts, iref, j);
                      if (j != rhscol && val > 0) {
                        rat = valn / val;
                        if (*s == 0 || rat < minrat) {
                          minrat = rat;
                          *s = j;
                          set_emptyset(data->stieset);
                          set_addelem(data->stieset, j);
                        } else {
                          if (rat == minrat) {
                            set_addelem(data->stieset, j);
                          }
                        }
                      }
                    }
                  }
                  set_copy(data->tieset, data->stieset);
                  if (set_card(data->tieset) == 1)
                    colselected = true;
                }
              }
              k += 1;
            } while (!colselected && k <= d_size);
            *selected = true;
          }
        } else {
          *lps = dd_Inconsistent;
        }
      }
    } /* end of while */
  }
}

void dd_SelectPreorderedNext2(dd_rowrange m_size, rowset excluded,
                              dd_rowindex OV, dd_rowrange *hnext) {
  dd_rowrange i, k;

  *hnext = 0;
  for (i = 1; i <= m_size && *hnext == 0; i++) {
    k = OV[i];
    if (!set_member(k, excluded))
      *hnext = k;
  }
}

template <typename T>
void dd_SelectPivot2(dd_rowrange m_size, dd_colrange d_size, T **A, T **Ts,
                     dd_rowindex ordervec, rowset equalityset,
                     dd_rowrange rowmax, rowset NopivotRow, colset NopivotCol,
                     dd_rowrange *r, dd_colrange *s, bool *selected)
/* Select a position (*r,*s) in the matrix A.T such that (A.T)[*r][*s] is
   nonzero The choice is feasible, i.e., not on NopivotRow and NopivotCol, and
   best with respect to the specified roworder
 */
{
  int stop;
  dd_rowrange i, rtemp;
  rowset rowexcluded;
  T Xtemp;

  stop = false;
  set_initialize(&rowexcluded, m_size);
  set_copy(rowexcluded, NopivotRow);
  for (i = rowmax + 1; i <= m_size; i++)
    set_addelem(rowexcluded, i); /* cannot pivot on any row > rmax */
  *selected = false;
  do {
    rtemp = 0;
    i = 1;
    while (i <= m_size &&
           rtemp == 0) { /* equalityset vars have highest priorities */
      if (set_member(i, equalityset) && !set_member(i, rowexcluded)) {
        rtemp = i;
      }
      i++;
    }
    if (rtemp == 0)
      dd_SelectPreorderedNext2(m_size, rowexcluded, ordervec, &rtemp);
    ;
    if (rtemp >= 1) {
      *r = rtemp;
      *s = 1;
      while (*s <= d_size && !*selected) {
        dd_TableauEntry(Xtemp, d_size, A, Ts, *r, *s);
        if (!set_member(*s, NopivotCol) && Xtemp != 0) {
          *selected = true;
          stop = true;
        } else {
          (*s)++;
        }
      }
      if (!*selected) {
        set_addelem(rowexcluded, rtemp);
      }
    } else {
      *r = 0;
      *s = 0;
      stop = true;
    }
  } while (!stop);
  set_free(rowexcluded);
}

template <typename T>
void dd_GaussianColumnPivot(dd_colrange d_size, T **X, T **Ts, dd_rowrange r,
                            dd_colrange s, T *Rtemp)
/* Update the Transformation matrix T with the pivot operation on (r,s)
   This procedure performs a implicit pivot operation on the matrix X by
   updating the dual basis inverse  T.
 */
{
  dd_colrange j, j1;
  T Xtemp0, Xtemp;

  for (j = 1; j <= d_size; j++) {
    dd_TableauEntry(Rtemp[j - 1], d_size, X, Ts, r, j);
  }
  bool is_field = true;
  if (is_field) {
    Xtemp0 = Rtemp[s - 1];
    for (j = 1; j <= d_size; j++) {
      if (j != s) {
        Xtemp = Rtemp[j - 1] / Xtemp0;
        for (j1 = 1; j1 <= d_size; j1++)
          Ts[j1 - 1][j - 1] -= Xtemp * Ts[j1 - 1][s - 1];
      }
    }
    for (j = 1; j <= d_size; j++)
      Ts[j - 1][s - 1] /= Xtemp0;
  } else {
    // ring case now
    for (j = 1; j <= d_size; j++) {
      if (j != s) {
        for (j1 = 1; j1 <= d_size; j1++) {
          Ts[j1 - 1][j - 1] = Rtemp[s - 1] * Ts[j1 - 1][j - 1] -
                              Rtemp[j - 1] * Ts[j1 - 1][s - 1];
        }
      }
    }
    Xtemp0 = Rtemp[s - 1];
    for (j = 1; j <= d_size; j++)
      Ts[j - 1][s - 1] /= Xtemp0;
    for (j = 1; j <= d_size; j++) {
      if (j != s) {
        T alpha;
        dd_TableauEntry(alpha, d_size, X, Ts, r, j);
        if (alpha != 0)
          std::cout << "j=" << j << " alpha=" << alpha << "\n";
      }
    }
  }
}

template <typename T>
void dd_GaussianColumnPivot2(dd_colrange d_size, T **A, T **Ts,
                             dd_colindex nbindex, dd_rowindex bflag,
                             dd_rowrange r, dd_colrange s, T *Rtemp)
/* Update the Transformation matrix T with the pivot operation on (r,s)
   This procedure performs a implicit pivot operation on the matrix A by
   updating the dual basis inverse  T.
 */
{
  long entering;

  dd_GaussianColumnPivot(d_size, A, Ts, r, s, Rtemp);
  entering = nbindex[s];
  bflag[r] = s;   /* the nonbasic variable r corresponds to column s */
  nbindex[s] = r; /* the nonbasic variable on s column is r */

  if (entering > 0)
    bflag[entering] = -1;
  /* original variables have negative index and should not affect the row index
   */
}

template <typename T> void dd_SetToIdentity(dd_colrange d_size, T **Ts) {
  dd_colrange j1, j2;

  for (j1 = 1; j1 <= d_size; j1++) {
    for (j2 = 1; j2 <= d_size; j2++) {
      if (j1 == j2)
        Ts[j1 - 1][j2 - 1] = 1;
      else
        Ts[j1 - 1][j2 - 1] = 0;
    }
  }
}

template <typename T>
void dd_ResetTableau(dd_rowrange m_size, dd_colrange d_size, T **Ts,
                     dd_colindex nbindex, dd_rowindex bflag, dd_rowrange objrow,
                     dd_colrange rhscol) {
  dd_rowrange i;
  dd_colrange j;

  /* Initialize T and nbindex */
  for (j = 1; j <= d_size; j++)
    nbindex[j] = -j;
  nbindex[rhscol] = 0;
  /* RHS is already in nonbasis and is considered to be associated
     with the zero-th row of input. */
  dd_SetToIdentity(d_size, Ts);

  /* Set the bflag according to nbindex */
  for (i = 1; i <= m_size; i++)
    bflag[i] = -1;
  /* all basic variables have index -1 */
  bflag[objrow] = 0;
  /* bflag of the objective variable is 0,
     different from other basic variables which have -1 */
  for (j = 1; j <= d_size; j++)
    if (nbindex[j] > 0)
      bflag[nbindex[j]] = j;
  /* bflag of a nonbasic variable is its column number */
}

template <typename T>
void dd_SelectCrissCrossPivot(dd_rowrange m_size, dd_colrange d_size, T **A,
                              T **Ts, dd_rowindex bflag, dd_rowrange objrow,
                              dd_colrange rhscol, dd_rowrange *r,
                              dd_colrange *s, bool *selected,
                              dd_LPStatusType *lps) {
  int colselected = false, rowselected = false;
  dd_rowrange i;
  T val;

  *selected = false;
  *lps = dd_LPSundecided;
  while ((*lps == dd_LPSundecided) && (!rowselected) && (!colselected)) {
    for (i = 1; i <= m_size; i++) {
      if (i != objrow && bflag[i] == -1) { /* i is a basic variable */
        dd_TableauEntry(val, d_size, A, Ts, i, rhscol);
        if (val < 0) {
          rowselected = true;
          *r = i;
          break;
        }
      } else {
        if (bflag[i] > 0) { /* i is nonbasic variable */
          dd_TableauEntry(val, d_size, A, Ts, objrow, bflag[i]);
          if (val > 0) {
            colselected = true;
            *s = bflag[i];
            break;
          }
        }
      }
    }
    if ((!rowselected) && (!colselected)) {
      *lps = dd_Optimal;
      return;
    } else if (rowselected) {
      for (i = 1; i <= m_size; i++) {
        if (bflag[i] > 0) { /* i is nonbasic variable */
          dd_TableauEntry(val, d_size, A, Ts, *r, bflag[i]);
          if (val > 0) {
            colselected = true;
            *s = bflag[i];
            *selected = true;
            break;
          }
        }
      }
    } else {
      if (colselected) {
        for (i = 1; i <= m_size; i++) {
          if (i != objrow && bflag[i] == -1) { /* i is a basic variable */
            dd_TableauEntry(val, d_size, A, Ts, i, *s);
            if (val < 0) {
              rowselected = true;
              *r = i;
              *selected = true;
              break;
            }
          }
        }
      }
    }
    if (!rowselected) {
      *lps = dd_DualInconsistent;
    } else {
      if (!colselected) {
        *lps = dd_Inconsistent;
      }
    }
  }
}

template <typename T>
void dd_CrissCrossSolve(dd_lpdata<T> *lp, dd_ErrorType *err,
                        data_temp_simplex<T> *data, size_t const& maxiter, std::ostream& os) {
  switch (lp->objective) {
  case dd_LPmax:
    dd_CrissCrossMaximize(lp, err, data, maxiter, os);
    break;

  case dd_LPmin:
    dd_CrissCrossMinimize(lp, err, data, maxiter, os);
    break;

  case dd_LPnone:
    *err = dd_NoLPObjective;
    break;
  }
}

template <typename T>
void dd_DualSimplexSolve(dd_lpdata<T> *lp, dd_ErrorType *err,
                         data_temp_simplex<T> *data,
                         size_t const& maxiter, std::ostream& os) {
  bool localdebug = false;
  if (localdebug)
    std::cout << "Running dd_DualSimplexSolve\n";
  switch (lp->objective) {
  case dd_LPmax:
#ifdef DEBUG_CDD
    os << "CDD: Before dd_DualSimplexMaximize in dd_DualSimplexSolve\n";
#endif
    dd_DualSimplexMaximize(lp, err, data, maxiter, os);
#ifdef DEBUG_CDD
    os << "CDD: After dd_DualSimplexMaximize in dd_DualSimplexSolve\n";
#endif
    break;

  case dd_LPmin:
#ifdef DEBUG_CDD
    os << "CDD: Before dd_DualSimplexMinimize in dd_DualSimplexSolve\n";
#endif
    dd_DualSimplexMinimize(lp, err, data, maxiter, os);
#ifdef DEBUG_CDD
    os << "CDD: After dd_DualSimplexMinimize in dd_DualSimplexSolve\n";
#endif
    break;

  case dd_LPnone:
    *err = dd_NoLPObjective;
    break;
  }
}

template <typename T>
void dd_FindLPBasis(dd_rowrange m_size, dd_colrange d_size, T **A, T **Ts,
                    dd_rowindex OV, dd_rowset equalityset, dd_colindex nbindex,
                    dd_rowindex bflag, dd_rowrange objrow, dd_colrange rhscol,
                    dd_colrange *cs, bool *found, dd_LPStatusType *lps,
                    long *pivot_no, data_temp_simplex<T> *data) {
  bool chosen, stop;
  long pivots_p0 = 0, rank;
  colset ColSelected;
  rowset RowSelected;
  T val;

  dd_rowrange r;
  dd_colrange j, s;

  *found = false;
  *cs = 0;
  rank = 0;
  stop = false;
  *lps = dd_LPSundecided;

  set_initialize(&RowSelected, m_size);
  set_initialize(&ColSelected, d_size);
  set_addelem(RowSelected, objrow);
  set_addelem(ColSelected, rhscol);

  stop = false;
  do { /* Find a LP basis */
    dd_SelectPivot2(m_size, d_size, A, Ts, OV, equalityset, m_size, RowSelected,
                    ColSelected, &r, &s, &chosen);
    if (chosen) {
      set_addelem(RowSelected, r);
      set_addelem(ColSelected, s);
      //      std::cout << "dd_GaussianColumnPivot2 call 1\n";
      dd_GaussianColumnPivot2(d_size, A, Ts, nbindex, bflag, r, s, data->Rtemp);
      pivots_p0++;
      rank++;
    } else {
      for (j = 1; j <= d_size && *lps == dd_LPSundecided; j++) {
        if (j != rhscol && nbindex[j] < 0) {
          dd_TableauEntry(val, d_size, A, Ts, objrow, j);
          if (val != 0) { /* dual inconsistent */
            *lps = dd_StrucDualInconsistent;
            *cs = j;
            /* dual inconsistent because the nonzero reduced cost */
          }
        }
      }
      if (*lps == dd_LPSundecided)
        *found = true;
      /* dependent columns but not dual inconsistent. */
      stop = true;
    }
    /* printf("d_size=%ld, rank=%ld\n",d_size,rank); */
    if (rank == d_size - 1) {
      stop = true;
      *found = true;
    }
  } while (!stop);

  *pivot_no = pivots_p0;
  set_free(RowSelected);
  set_free(ColSelected);
}

template <typename T>
void dd_FindLPBasis2(dd_rowrange m_size, dd_colrange d_size, T **A, T **Ts,
                     dd_rowindex OV, dd_rowset equalityset, dd_colindex nbindex,
                     dd_rowindex bflag, dd_rowrange objrow, dd_colrange rhscol,
                     dd_colrange *cs, bool *found, long *pivot_no,
                     data_temp_simplex<T> *data) {
  /* Similar to dd_FindLPBasis but it is much simpler.  This tries to recompute
  T for the specified basis given by nbindex.  It will return *found=false if
  the specified basis is not a basis.
  */
  bool chosen, stop;
  long pivots_p0 = 0, rank;
  dd_colset ColSelected, DependentCols;
  dd_rowset RowSelected, NopivotRow;
  T val;
  bool localdebug = false;

  dd_rowrange r, negcount = 0;
  dd_colrange j, s;

  *found = false;
  *cs = 0;
  rank = 0;

  set_initialize(&RowSelected, m_size);
  set_initialize(&DependentCols, d_size);
  set_initialize(&ColSelected, d_size);
  set_initialize(&NopivotRow, m_size);
  set_addelem(RowSelected, objrow);
  set_addelem(ColSelected, rhscol);
  set_compl(NopivotRow, NopivotRow); /* set NopivotRow to be the groundset */

  for (j = 2; j <= d_size; j++)
    if (nbindex[j] > 0) {
      set_delelem(NopivotRow, nbindex[j]);
    } else {
      if (nbindex[j] < 0) {
        negcount++;
        set_addelem(DependentCols, -nbindex[j]);
        set_addelem(ColSelected, -nbindex[j]);
      }
    }

  set_uni(
      RowSelected, RowSelected,
      NopivotRow); /* RowSelected is the set of rows not allowed to poviot on */

  stop = false;
  do { /* Find a LP basis */
    dd_SelectPivot2(m_size, d_size, A, Ts, OV, equalityset, m_size, RowSelected,
                    ColSelected, &r, &s, &chosen);
    if (chosen) {
      set_addelem(RowSelected, r);
      set_addelem(ColSelected, s);
      //      std::cout << "dd_GaussianColumnPivot2 call 2\n";
      dd_GaussianColumnPivot2(d_size, A, Ts, nbindex, bflag, r, s, data->Rtemp);
      pivots_p0++;
      rank++;
    } else {
      *found = false; /* cannot pivot on any of the spacified positions. */
      stop = true;
    }
    if (rank == d_size - 1 - negcount) {
      if (negcount) {
        /* Now it tries to pivot on rows that are supposed to be dependent. */
        set_diff(ColSelected, ColSelected, DependentCols);
        dd_SelectPivot2(m_size, d_size, A, Ts, OV, equalityset, m_size,
                        RowSelected, ColSelected, &r, &s, &chosen);
        if (chosen)
          *found = false; /* not supposed to be independent */
        else
          *found = true;
        if (localdebug) {
          printf("Try to check the dependent cols:");
          if (chosen)
            printf("They are not dependent.  Can still pivot on (%ld, %ld)\n",
                   r, s);
          else
            printf("They are indeed dependent.\n");
        }
      } else {
        *found = true;
      }
      stop = true;
    }
  } while (!stop);

  for (j = 1; j <= d_size; j++)
    if (nbindex[j] > 0)
      bflag[nbindex[j]] = j;
  *pivot_no = pivots_p0;
  set_free(RowSelected);
  set_free(ColSelected);
  set_free(NopivotRow);
  set_free(DependentCols);
}

template <typename T>
void dd_FindDualFeasibleBasis(dd_rowrange m_size, dd_colrange d_size, T **A,
                              T **Ts, dd_colindex nbindex, dd_rowindex bflag,
                              dd_rowrange objrow, dd_colrange rhscol,
                              bool lexicopivot, dd_colrange *s,
                              dd_ErrorType *err, dd_LPStatusType *lps,
                              long *pivot_no, long maxpivots,
                              data_temp_simplex<T> *data) {
  /* Find a dual feasible basis using Phase I of Dual Simplex method.
     If the problem is dual feasible,
     the procedure returns *err=NoError, *lps=LPSundecided and a dual feasible
     basis.   If the problem is dual infeasible, this returns
     *err=NoError, *lps=DualInconsistent and the evidence column *s.
     Caution: matrix A must have at least one extra row:  the row space
     A[m_size] must have been allocated.
  */
  bool phase1, dualfeasible = true;
  bool localdebug = false, chosen, stop;
  dd_LPStatusType LPSphase1;
  long pivots_p1 = 0;
  dd_rowrange i, r_val;
  dd_colrange j, l, ms = 0, s_val, local_m_size;
  T x, val, maxcost(0), axvalue, maxratio(0), maxratio_q(1);

  T scaling; /* random scaling mytype value */
  T minval(0);

  *err = dd_NoError;
  *lps = dd_LPSundecided;
  *s = 0;
  local_m_size = m_size + 1; /* increase m_size by 1 */

  ms =
      0; /* ms will be the index of column which has the largest reduced cost */
  for (j = 1; j <= d_size; j++) {
    if (j != rhscol) {
      dd_TableauEntry(data->rcost[j - 1], d_size, A, Ts, objrow, j);
      if (data->rcost[j - 1] > maxcost) {
        maxcost = data->rcost[j - 1];
        ms = j;
      }
    }
  }
  if (maxcost > 0)
    dualfeasible = false;

  if (!dualfeasible) {
    for (j = 1; j <= d_size; j++) {
      A[local_m_size - 1][j - 1] = 0;
      for (l = 1; l <= d_size; l++) {
        if (nbindex[l] > 0) {
          scaling = l + 10;
          A[local_m_size - 1][j - 1] -= A[nbindex[l] - 1][j - 1] * scaling;
          /* To make the auxiliary row (0,-11,-12,...,-d-10).
             It is likely to be better than  (0, -1, -1, ..., -1)
             to avoid a degenerate LP.  Version 093c. */
        }
      }
    }

    ms = 0;
    /* Ratio Test: ms will be now the index of column which has the largest
       reduced cost over the auxiliary row entry */
    for (j = 1; j <= d_size; j++) {
      if ((j != rhscol) && data->rcost[j - 1] > 0) {
        dd_TableauEntry(axvalue, d_size, A, Ts, local_m_size, j);
        if (axvalue >= 0) {
          *err = dd_NumericallyInconsistent;
          /* This should not happen as they are set negative above.  Quit the
           * phase I.*/
          goto _L99;
        }
        // So now axvalue < 0
        axvalue = -axvalue;
        // So now axvalue > 0
        bool is_field = true;
        if (is_field) {
          axvalue =
              data->rcost[j - 1] / axvalue; /* axvalue is the negative of ratio
                                               that is to be maximized. */
          if (axvalue > maxratio) {
            maxratio = axvalue;
            ms = j;
          }
        } else {
          if (dd_LargerFrac(data->rcost[j - 1], axvalue, maxratio,
                            maxratio_q)) {
            maxratio = data->rcost[j - 1];
            maxratio_q = axvalue;
            ms = j;
          }
        }
      }
    }

    if (ms == 0) {
      *err = dd_NumericallyInconsistent; /* This should not happen. Quit the
                                            phase I.*/
      goto _L99;
    }

    /* Pivot on (local_m_size,ms) so that the dual basic solution becomes
     * feasible */
    //    std::cout << "dd_GaussianColumnPivot2 call 3\n";
    dd_GaussianColumnPivot2(d_size, A, Ts, nbindex, bflag, local_m_size, ms,
                            data->Rtemp);
    pivots_p1++;
    if (localdebug) {
      printf("\ndd_FindDualFeasibleBasis: Pivot on %ld %ld.\n", local_m_size,
             ms);
    }

    for (j = 1; j <= d_size; j++)
      data->nbindex_ref[j] = nbindex[j];
    /* set the reference basis to be the current feasible basis. */

    phase1 = true;
    stop = false;
    do { /* Dual Simplex Phase I */
      chosen = false;
      LPSphase1 = dd_LPSundecided;
      if (pivots_p1 > maxpivots) {
        *err = dd_LPCycling;
        goto _L99; /* failure due to max no. of pivots performed */
      }
      dd_SelectDualSimplexPivot(
          local_m_size, d_size, phase1, A, Ts, data->nbindex_ref, bflag, objrow,
          rhscol, lexicopivot, &r_val, &s_val, &chosen, &LPSphase1, data);
      if (!chosen) {
        /* The current dictionary is terminal.  There are two cases:
           dd_TableauEntry(local_m_size,d_size,A,T,objrow,ms) is negative or
           zero. The first case implies dual infeasible, and the latter implies
           dual feasible but local_m_size is still in nonbasis. We must pivot in
           the auxiliary variable local_m_size.
        */
        dd_TableauEntry(x, d_size, A, Ts, objrow, ms);
        if (x < 0) {
          *err = dd_NoError;
          *lps = dd_DualInconsistent;
          *s = ms;
        }

        r_val = 0;
        for (i = 1; i <= local_m_size; i++) {
          if (bflag[i] < 0) {
            /* i is basic and not the objective variable */
            dd_TableauEntry(val, d_size, A, Ts, i, ms); /* auxiliary column*/
            if (val < minval) {
              r_val = i;
              minval = val;
            }
          }
        }

        if (r_val == 0) {
          *err = dd_NumericallyInconsistent; /* This should not happen. Quit the
                                                phase I.*/
          goto _L99;
        }

        //        std::cout << "dd_GaussianColumnPivot2 call 4\n";
        dd_GaussianColumnPivot2(d_size, A, Ts, nbindex, bflag, r_val, ms,
                                data->Rtemp);
        stop = true;
      } else {
        //        std::cout << "dd_GaussianColumnPivot2 call 5\n";
        dd_GaussianColumnPivot2(d_size, A, Ts, nbindex, bflag, r_val, s_val,
                                data->Rtemp);
        if (bflag[local_m_size] < 0)
          stop = true;
      }
      pivots_p1++;
    } while (!stop);
  }
_L99:
  *pivot_no = pivots_p1;
}

template <typename T>
void dd_DualSimplexMinimize(dd_lpdata<T> *lp, dd_ErrorType *err,
                            data_temp_simplex<T> *data,
                            size_t const& maxiter, [[maybe_unused]] std::ostream& os) {
  dd_colrange j;
  *err = dd_NoError;
  for (j = 1; j <= lp->d; j++)
    lp->A[lp->objrow - 1][j - 1] = -lp->A[lp->objrow - 1][j - 1];
  dd_DualSimplexMaximize(lp, err, data, maxiter, os);
  lp->optvalue = -lp->optvalue;
  for (j = 1; j <= lp->d; j++) {
    if (lp->LPS != dd_Inconsistent)
      lp->dsol[j - 1] = -lp->dsol[j - 1];
    lp->A[lp->objrow - 1][j - 1] = -lp->A[lp->objrow - 1][j - 1];
  }
}

template <typename T>
void dd_DualSimplexMaximize(dd_lpdata<T> *lp, dd_ErrorType *err,
                            data_temp_simplex<T> *data,
                            size_t const& maxiter, std::ostream& os)
/*
When LP is inconsistent then lp->re returns the evidence row.
When LP is dual-inconsistent then lp->se returns the evidence column.
*/
{
  bool stop, chosen, phase1, found;
  long pivots_ds = 0, pivots_p0 = 0, pivots_p1 = 0, pivots_pc = 0, maxpivots,
       maxpivfactor = 20;
  bool localdebug1 = false;
  dd_rowrange i, r;
  dd_colrange j, s;
  size_t n_iter;

  unsigned int rseed = 1;
#ifdef DEBUG_CDD
  os << "CDD: Beginning of dd_DualSimplexMaximize with maxiter=" << maxiter << "\n";
#endif
  // r value assigned designed to create segfault in case it is not set later
  r = -40000;

  /* *err=dd_NoError; */
  set_emptyset(lp->redset_extra);
  for (i = 0; i <= 4; i++)
    lp->pivots[i] = 0;
  maxpivots =
      maxpivfactor *
      lp->d; /* maximum pivots to be performed before cc pivot is applied. */
  /* Initializing control variables. */
#ifdef DEBUG_CDD
  os << "CDD: dd_DualSimplexMaximize, step 1\n";
#endif
  dd_ComputeRowOrderVector2(lp->m, lp->d, lp->A, data->OrderVector, dd_MinIndex,
                            rseed);
#ifdef DEBUG_CDD
  os << "CDD: dd_DualSimplexMaximize, step 2\n";
#endif

  lp->re = 0;
  lp->se = 0;

  dd_ResetTableau(lp->m, lp->d, lp->B, lp->nbindex, data->bflag, lp->objrow,
                  lp->rhscol);
#ifdef DEBUG_CDD
  os << "CDD: dd_DualSimplexMaximize, step 3\n";
#endif

  dd_FindLPBasis(lp->m, lp->d, lp->A, lp->B, data->OrderVector, lp->equalityset,
                 lp->nbindex, data->bflag, lp->objrow, lp->rhscol, &s, &found,
                 &(lp->LPS), &pivots_p0, data);
#ifdef DEBUG_CDD
  os << "CDD: dd_DualSimplexMaximize, step 4\n";
#endif
  lp->pivots[0] = pivots_p0;

  if (!found) {
    lp->se = s;
    goto _L99;
    /* No LP basis is found, and thus Inconsistent.
    Output the evidence column. */
  }

  dd_FindDualFeasibleBasis(lp->m, lp->d, lp->A, lp->B, lp->nbindex, data->bflag,
                           lp->objrow, lp->rhscol, lp->lexicopivot, &s, err,
                           &(lp->LPS), &pivots_p1, maxpivots, data);
#ifdef DEBUG_CDD
  os << "CDD: dd_DualSimplexMaximize, step 5\n";
#endif
  lp->pivots[1] = pivots_p1;

  /* set the reference basis to be the current feasible basis. */
  for (j = 1; j <= lp->d; j++)
    data->nbindex_ref_ds[j] = lp->nbindex[j];

#ifdef DEBUG_CDD
  os << "CDD: dd_DualSimplexMaximize, step 6\n";
#endif
  if (*err == dd_LPCycling || *err == dd_NumericallyInconsistent) {
#ifdef DEBUG_CDD
    os << "CDD: dd_DualSimplexMaximize, calling dd_CrissCrossMaximize\n";
#endif
    dd_CrissCrossMaximize(lp, err, data, maxiter, os);
    return;
  }
#ifdef DEBUG_CDD
  os << "CDD: dd_DualSimplexMaximize, step 7\n";
#endif

  if (lp->LPS == dd_DualInconsistent) {
    lp->se = s;
    goto _L99;
    /* No dual feasible basis is found, and thus DualInconsistent.
    Output the evidence column. */
  }
#ifdef DEBUG_CDD
  os << "CDD: dd_DualSimplexMaximize, step 8\n";
#endif

  /* Dual Simplex Method */
  stop = false;
  n_iter = 0;
  do {
    n_iter += 1;
#ifdef DEBUG_CDD_DISABLE
    os << "CDD: dd_DualSimplexMaximize n_iter=" << n_iter << "\n";
#endif
    if (maxiter != 0) {
      if (n_iter == maxiter) {
#ifdef DEBUG_CDD
        os << "CDD: Exiting from dd_DualSimplexMaximize at maxiter=" << maxiter << " (too many iterations)\n";
#endif
        lp->LPS = dd_TooManyIterations;
        break;
      }
    }
    chosen = false;
    lp->LPS = dd_LPSundecided;
    phase1 = false;
    if (pivots_ds < maxpivots) {
      dd_SelectDualSimplexPivot(lp->m, lp->d, phase1, lp->A, lp->B,
                                data->nbindex_ref_ds, data->bflag, lp->objrow,
                                lp->rhscol, lp->lexicopivot, &r, &s, &chosen,
                                &(lp->LPS), data);
    }
    if (chosen) {
      pivots_ds++;
      if (lp->redcheck_extensive) {
        dd_GetRedundancyInformation(lp->m, lp->d, lp->A, lp->B, data->bflag,
                                    lp->redset_extra);
        set_uni(lp->redset_accum, lp->redset_accum, lp->redset_extra);
      }
    }
    if (!chosen && lp->LPS == dd_LPSundecided) {
      if (localdebug1) {
        std::cout
            << "Warning: an emergency CC pivot in Phase II is performed\n";
        /* In principle this should not be executed because we already have dual
           feasibility attained and dual simplex pivot should have been chosen.
           This might occur under floating point computation, or the case of
           cycling.
        */
      }

      dd_SelectCrissCrossPivot(lp->m, lp->d, lp->A, lp->B, data->bflag,
                               lp->objrow, lp->rhscol, &r, &s, &chosen,
                               &(lp->LPS));
      if (chosen)
        pivots_pc++;
    }
    if (chosen) {
      //      std::cout << "dd_GaussianColumnPivot2 call 6\n";
      dd_GaussianColumnPivot2(lp->d, lp->A, lp->B, lp->nbindex, data->bflag, r,
                              s, data->Rtemp);
    } else {
      switch (lp->LPS) {
      case dd_Inconsistent:
        lp->re = r;
        break;
      case dd_DualInconsistent:
        lp->se = s;
        break;
      case dd_LPSundecided:
        break;
      case dd_Unbounded:
        break;
      case dd_DualUnbounded:
        break;
      case dd_StrucInconsistent:
        break;
      case dd_StrucDualInconsistent:
        break;
      case dd_Optimal:
        break;
      case dd_TooManyIterations:
#ifdef DEBUG_CDD
        os << "CDD: We should not pass here\n";
        throw TerminalException{1};
#endif
        break;
      }
      stop = true;
    }
  } while (!stop);
#ifdef DEBUG_CDD
  os << "CDD: Terminates at maxiter=" << maxiter << "\n";
#endif
_L99:
  lp->pivots[2] = pivots_ds;
  lp->pivots[3] = pivots_pc;
#ifdef DEBUG_CDD
  os << "CDD: Before dd_SetSolutions in dd_DualSimplexMaximize\n";
#endif
  dd_SetSolutions(lp->m, lp->d, lp->A, lp->B, lp->objrow, lp->rhscol, lp->LPS,
                  lp->optvalue, lp->sol, lp->dsol, lp->posset_extra, lp->re,
                  lp->se, data->bflag, os);
}

template <typename T>
void dd_CrissCrossMinimize(dd_lpdata<T> *lp, dd_ErrorType *err,
                           data_temp_simplex<T> *data, size_t const& maxiter, std::ostream& os) {
  dd_colrange j;

  *err = dd_NoError;
  for (j = 1; j <= lp->d; j++)
    lp->A[lp->objrow - 1][j - 1] = -lp->A[lp->objrow - 1][j - 1];
  dd_CrissCrossMaximize(lp, err, data, maxiter, os);
  lp->optvalue = -lp->optvalue;
  for (j = 1; j <= lp->d; j++) {
    if (lp->LPS != dd_Inconsistent) {
      /* Inconsistent certificate stays valid for minimization, 0.94e */
      lp->dsol[j - 1] = -lp->dsol[j - 1];
    }
    lp->A[lp->objrow - 1][j - 1] = -lp->A[lp->objrow - 1][j - 1];
  }
}

template <typename T>
void dd_CrissCrossMaximize(dd_lpdata<T> *lp, dd_ErrorType *err,
                           data_temp_simplex<T> *data,
                           size_t const& maxiter, std::ostream& os)
/*
When LP is inconsistent then lp->re returns the evidence row.
When LP is dual-inconsistent then lp->se returns the evidence column.
*/
{
  bool stop, chosen, found;
  long pivots0, pivots1;

  dd_rowrange i, r;
  dd_colrange s;
  unsigned int rseed = 1;
  size_t n_iter;

  *err = dd_NoError;
  std::vector<long> nbtemp(lp->d + 1, 0);
  for (i = 0; i <= 4; i++)
    lp->pivots[i] = 0;
  /* Initializing control variables. */
  dd_ComputeRowOrderVector2(lp->m, lp->d, lp->A, data->OrderVector, dd_MinIndex,
                            rseed);

  lp->re = 0;
  lp->se = 0;
  pivots1 = 0;

  dd_ResetTableau(lp->m, lp->d, lp->B, lp->nbindex, data->bflag, lp->objrow,
                  lp->rhscol);

  dd_FindLPBasis(lp->m, lp->d, lp->A, lp->B, data->OrderVector, lp->equalityset,
                 lp->nbindex, data->bflag, lp->objrow, lp->rhscol, &s, &found,
                 &(lp->LPS), &pivots0, data);
  lp->pivots[0] += pivots0;

  if (!found) {
    lp->se = s;
    goto _L99;
    /* No LP basis is found, and thus Inconsistent.
    Output the evidence column. */
  }

  stop = false;
  n_iter = 0;
  do { /* Criss-Cross Method */
    n_iter += 1;
#ifdef DEBUG_CDD_DISABLE
    os << "CDD: dd_CrissCrossMaximize n_iter=" << n_iter << "\n";
#endif
    if (maxiter != 0) {
      if (n_iter == maxiter) {
#ifdef DEBUG_CDD
        os << "CDD: Exiting from dd_CrissCrossMaximize at maxiter=" << maxiter << " (too many iterations)\n";
#endif
        lp->LPS = dd_TooManyIterations;
        break;
      }
    }

    dd_SelectCrissCrossPivot(lp->m, lp->d, lp->A, lp->B, data->bflag,
                             lp->objrow, lp->rhscol, &r, &s, &chosen,
                             &(lp->LPS));
    if (chosen) {
      //      std::cout << "dd_GaussianColumnPivot2 call 7\n";
      dd_GaussianColumnPivot2(lp->d, lp->A, lp->B, lp->nbindex, data->bflag, r,
                              s, data->Rtemp);
      pivots1++;
    } else {
      switch (lp->LPS) {
      case dd_Inconsistent:
        lp->re = r;
        break;
      case dd_DualInconsistent:
        lp->se = s;
        break;
      case dd_Optimal:
        break;
      case dd_LPSundecided:
        break;
      case dd_Unbounded:
        break;
      case dd_DualUnbounded:
        break;
      case dd_StrucInconsistent:
        break;
      case dd_StrucDualInconsistent:
        break;
      case dd_TooManyIterations:
#ifdef DEBUG_CDD
        std::cerr << "CDD: That case should not match\n";
        throw TerminalException{1};
#endif
        break;
      }
      stop = true;
    }
  } while (!stop);

_L99:
  lp->pivots[1] += pivots1;
#ifdef DEBUG_CDD
  os << "CDD: Before dd_SetSolutions in dd_CrissCrossMaximize\n";
#endif
  dd_SetSolutions(lp->m, lp->d, lp->A, lp->B, lp->objrow, lp->rhscol, lp->LPS,
                  lp->optvalue, lp->sol, lp->dsol, lp->posset_extra, lp->re,
                  lp->se, data->bflag, os);
}

template <typename T>
void dd_SetSolutions(dd_rowrange m_size, dd_colrange d_size, T **A, T **Ts,
                     dd_rowrange objrow, dd_colrange rhscol,
                     dd_LPStatusType LPS, T &optvalue, T *sol, T *dsol,
                     dd_rowset posset, dd_rowrange re, dd_colrange se,
                     dd_rowindex bflag, [[maybe_unused]] std::ostream& os)
/*
Assign the solution vectors to sol,dsol,*optvalue after solving
the LP.
*/
{
  dd_rowrange i;
  dd_colrange j;
  T x, sw;
  int localdebug = false;

  if (LPS == dd_TooManyIterations) {
    return;
  }

  if (localdebug)
    std::cout << "SetSolutions:\n";
  switch (LPS) {
  case dd_Optimal:
    for (j = 1; j <= d_size; j++) {
      sol[j - 1] = Ts[j - 1][rhscol - 1];
      dd_TableauEntry(x, d_size, A, Ts, objrow, j);
      dsol[j - 1] = -x;
      dd_TableauEntry(optvalue, d_size, A, Ts, objrow, rhscol);
    }
    for (i = 1; i <= m_size; i++) {
      if (bflag[i] == -1) { /* i is a basic variable */
        dd_TableauEntry(x, d_size, A, Ts, i, rhscol);
        if (x > 0)
          set_addelem(posset, i);
      }
    }

    break;
  case dd_Inconsistent:
    if (localdebug)
      std::cout << "SetSolutions: LP is inconsistent.\n";
    for (j = 1; j <= d_size; j++) {
      sol[j - 1] = Ts[j - 1][rhscol - 1];
      dd_TableauEntry(x, d_size, A, Ts, re, j);
      dsol[j - 1] = -x;
    }
    break;

  case dd_LPSundecided:
    std::cout
        << "Case dd_LPSundecided has not been programmed in dd_SetSolutions\n";
    throw TerminalException{1};

  case dd_StrucInconsistent:
    std::cout << "Case dd_StrucInconsistent has not been programmed in "
                 "dd_SetSolutions\n";
    throw TerminalException{1};

  case dd_Unbounded:
    std::cout
        << "Case dd_Unbounded has not been programmed in dd_SetSolutions\n";
    throw TerminalException{1};

  case dd_DualUnbounded:
    std::cout
        << "Case dd_DualUnbounded has not been programmed in dd_SetSolutions\n";
    throw TerminalException{1};

  case dd_DualInconsistent:
    if (localdebug)
      printf("SetSolutions: LP is dual inconsistent.\n");
    for (j = 1; j <= d_size; j++) {
      sol[j - 1] = Ts[j - 1][se - 1];
      dd_TableauEntry(x, d_size, A, Ts, objrow, j);
      dsol[j - 1] = -x;
    }
    break;

  case dd_TooManyIterations:
#ifdef DEBUG_CDD
    os << "CDD: That case should not be met here\n";
    throw TerminalException{1};
#endif
    break;

  case dd_StrucDualInconsistent:
    dd_TableauEntry(x, d_size, A, Ts, objrow, se);
    if (x > 0)
      sw = 1;
    else
      sw = -1;
    for (j = 1; j <= d_size; j++) {
      sol[j - 1] = sw * Ts[j - 1][se - 1];
      dd_TableauEntry(x, d_size, A, Ts, objrow, j);
      dsol[j - 1] = -x;
    }
    if (localdebug)
      std::cout << "SetSolutions: LP is dual inconsistent.\n";
    break;
  }

}

template <typename T>
void dd_ComputeRowOrderVector2(dd_rowrange m_size, dd_colrange d_size, T **A,
                               dd_rowindex OV, dd_RowOrderType ho,
                               unsigned int rseed) {
  long i, itemp;

  OV[0] = 0;
  switch (ho) {
  case dd_MaxIndex:
    for (i = 1; i <= m_size; i++)
      OV[i] = m_size - i + 1;
    break;

  case dd_LexMin:
    for (i = 1; i <= m_size; i++)
      OV[i] = i;
    dd_QuickSort(OV, 1, m_size, A, d_size);
    break;

  case dd_LexMax:
    for (i = 1; i <= m_size; i++)
      OV[i] = i;
    dd_QuickSort(OV, 1, m_size, A, d_size);
    for (i = 1; i <= m_size / 2; i++) { /* just reverse the order */
      itemp = OV[i];
      OV[i] = OV[m_size - i + 1];
      OV[m_size - i + 1] = itemp;
    }
    break;

  case dd_RandomRow:
    for (i = 1; i <= m_size; i++)
      OV[i] = i;
    if (rseed <= 0)
      rseed = 1;
    dd_RandomPermutation(OV, m_size, rseed);
    break;

  case dd_MinIndex:
    for (i = 1; i <= m_size; i++)
      OV[i] = i;
    break;

  case dd_MinCutoff:
  case dd_MaxCutoff:
  case dd_MixCutoff:
    for (i = 1; i <= m_size; i++)
      OV[i] = i;
    break;
  }
}

template <typename T> dd_lpdata<double> *dd_LPgmp2LPf(dd_lpdata<T> *lp) {
  dd_rowrange i;
  dd_colrange j;
  dd_lpdata<double> *lpf;
  bool localdebug = false;

  if (localdebug)
    std::cout << "Converting a GMP-LP to a float-LP.\n";

  lpf = dd_CreateLPData<double>(lp->m, lp->d);
  lpf->objective = lp->objective;

  for (i = 1; i <= lp->m; i++) {
    if (set_member(i, lp->equalityset))
      set_addelem(lpf->equalityset, i);
    /* it is equality. Its reversed row will not be in this set */
    for (j = 1; j <= lp->d; j++)
      lpf->A[i - 1][j - 1] =
          UniversalScalarConversion<double, T>(lp->A[i - 1][j - 1]);
  }
  return lpf;
}

template <typename T>
inline bool dd_LPSolve_data(dd_lpdata<T> *lp, dd_LPSolverType solver,
                            dd_ErrorType *err, data_temp_simplex<T> *data,
                            size_t const& maxiter, std::ostream& os)
/*
The current version of dd_LPSolve that solves an LP with floating-arithmetics
first and then with the specified arithimetics if it is GMP.

When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  int i;
  bool found = false;

  *err = dd_NoError;
  lp->solver = solver;

  // There is a bug when using USE_DOUBLE_FIRST.
  switch (lp->solver) {
  case dd_CrissCross:
    dd_CrissCrossSolve(lp, err, data, maxiter, os);
    break;
  case dd_DualSimplex:
#ifdef DEBUG_CDD
    os << "CDD: Before dd_DualSimplexSolve in dd_LPSolve_data\n";
#endif
    dd_DualSimplexSolve(lp, err, data, maxiter, os);
#ifdef DEBUG_CDD
    os << "CDD: After dd_DualSimplexSolve in dd_LPSolve_data\n";
#endif
    break;
  }

  lp->total_pivots = 0;
  for (i = 0; i <= 4; i++)
    lp->total_pivots += lp->pivots[i];
  if (*err == dd_NoError)
    found = true;
  return found;
}

template <typename T>
bool dd_LPSolve(dd_lpdata<T> *lp, dd_LPSolverType solver, dd_ErrorType *err, size_t const& maxiter, std::ostream& os) {
#ifdef DEBUG_CDD
  os << "CDD: dd_LPSolve, maxiter=" << maxiter << "\n";
#endif
  data_temp_simplex<T> *data = allocate_data_simplex<T>(lp->m, lp->d);
  bool test = dd_LPSolve_data(lp, solver, err, data, maxiter, os);
  free_data_simplex(data);
  return test;
}

template <typename T>
bool dd_LPSolve0(dd_lpdata<T> *lp, dd_LPSolverType solver, dd_ErrorType *err, size_t const& maxiter, std::ostream& os)
/*
The original version of dd_LPSolve that solves an LP with specified
arithimetics.

When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  int i;
  bool found = false;

  *err = dd_NoError;
  lp->solver = solver;

  data_temp_simplex<T> *data = allocate_data_simplex<T>(lp->m, lp->d);
  switch (lp->solver) {
  case dd_CrissCross:
    dd_CrissCrossSolve(lp, err, data, maxiter, os);
    break;
  case dd_DualSimplex:
#ifdef DEBUG_CDD
    os << "CDD: Before dd_DualSimplexSolve in dd_LPSolve0\n";
#endif
    dd_DualSimplexSolve(lp, err, data, maxiter, os);
#ifdef DEBUG_CDD
    os << "CDD: After dd_DualSimplexSolve in dd_LPSolve0\n";
#endif
    break;
  }
  free_data_simplex(data);

  lp->total_pivots = 0;
  for (i = 0; i <= 4; i++)
    lp->total_pivots += lp->pivots[i];
  if (*err == dd_NoError)
    found = true;
  return found;
}

template <typename T>
dd_lpdata<T> *dd_MakeLPforInteriorFinding(dd_lpdata<T> *lp)
/* Delete the objective row,
   add an extra column with -1's to the matrix A,
   add an extra row with (bceil, 0,...,0,-1),
   add an objective row with (0,...,0,1), and
   rows & columns, and change m_size and d_size accordingly, to output new_A.
  This sets up the LP:
  maximize      x_{d+1}
  s.t.    A x + x_{d+1}  <=  b
                x_{d+1}  <=  bm * bmax,
  where bm is set to 2 by default, and bmax=max{1, b[1],...,b[m_size]}.
  Note that the equalitions (linearity) in the input lp will be ignored.
*/
{
  dd_rowrange m;
  dd_colrange d;
  dd_lpdata<T> *lpnew;
  dd_rowrange i;
  dd_colrange j;
  T bm, bmax, bceil;

  bm = 2;
  bmax = 1;
  m = lp->m + 1;
  d = lp->d + 1;

  lpnew = dd_CreateLPData<T>(m, d);
  lpnew->objective = dd_LPmax;

  for (i = 1; i <= lp->m; i++)
    if (lp->A[i - 1][lp->rhscol - 1] > bmax)
      bmax = lp->A[i - 1][lp->rhscol - 1];
  bceil = bm * bmax;

  for (i = 1; i <= lp->m; i++)
    for (j = 1; j <= lp->d; j++)
      lpnew->A[i - 1][j - 1] = lp->A[i - 1][j - 1];

  for (i = 1; i <= lp->m; i++)
    lpnew->A[i - 1][lp->d] = -1;

  for (j = 1; j <= lp->d; j++)
    lpnew->A[m - 2][j - 1] = 0; /* new row (bceil, 0,...,0,-1) */
  lpnew->A[m - 2][0] = bceil;   /* new row (bceil, 0,...,0,-1) */

  for (j = 1; j <= d - 1; j++)
    lpnew->A[m - 1][j - 1] = 0; /* new obj row with (0,...,0,1) */
  lpnew->A[m - 1][d - 1] = 1;
  return lpnew;
}

template <typename T>
dd_lpdata<T> *dd_CreateLP_H_ImplicitLinearity(dd_matrixdata<T> *M) {
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_lpdata<T> *lp;
  bool localdebug = false;

  linc = set_card(M->linset);
  m = M->rowsize + 1 + linc + 1;
  /* We represent each equation by two inequalities.
     This is not the best way but makes the code simple. */
  d = M->colsize + 1;

  lp = dd_CreateLPData<T>(m, d);
  lp->objective = dd_LPmax;
  lp->redcheck_extensive = false; /* this is default */

  irev = M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    if (set_member(i, M->linset)) {
      irev++;
      set_addelem(lp->equalityset, i); /* it is equality. */
      /* the reversed row irev is not in the equality set. */
      for (j = 1; j <= M->colsize; j++) {
        lp->A[irev - 1][j - 1] = -M->matrix[i - 1][j - 1];
      }
    } else {
      lp->A[i - 1][d - 1] = -1; /* b_I + A_I x - 1 z >= 0  (z=x_d) */
    }
    for (j = 1; j <= M->colsize; j++) {
      lp->A[i - 1][j - 1] = M->matrix[i - 1][j - 1];
    } /*of j*/
  } /*of i*/
  lp->A[m - 2][0] = 1;
  lp->A[m - 2][d - 1] = -1;
  /* make the LP bounded.  */

  lp->A[m - 1][d - 1] = 1;
  /* objective is to maximize z.  */

  if (localdebug) {
    std::cout << "dd_CreateLP_H_ImplicitLinearity: an new lp is\n";
    dd_WriteLP(std::cout, lp);
  }

  return lp;
}

template <typename T>
dd_lpdata<T> *dd_CreateLP_V_ImplicitLinearity(dd_matrixdata<T> *M) {
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_lpdata<T> *lp;
  bool localdebug = false;

  linc = set_card(M->linset);
  m = M->rowsize + 1 + linc + 1;
  /* We represent each equation by two inequalities.
     This is not the best way but makes the code simple. */
  d = (M->colsize) + 2;
  /* Two more columns.  This is different from the H-reprentation case */

  /* The below must be modified for V-representation!!!  */

  lp = dd_CreateLPData<T>(m, d);
  lp->objective = dd_LPmax;
  lp->redcheck_extensive = false; /* this is default */

  irev = M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    lp->A[i - 1][0] = 0; /* It is almost completely degerate LP */
    if (set_member(i, M->linset)) {
      irev++;
      set_addelem(lp->equalityset, i); /* it is equality. */
      /* the reversed row irev is not in the equality set. */
      for (j = 2; j <= d - 1; j++) {
        lp->A[irev - 1][j - 1] = -M->matrix[i - 1][j - 2];
      } /*of j*/
      if (localdebug)
        fprintf(stdout, "equality row %ld generates the reverse row %ld.\n", i,
                irev);
    } else {
      lp->A[i - 1][d - 1] = -1;
    }
    for (j = 2; j <= d - 1; j++) {
      lp->A[i - 1][j - 1] = M->matrix[i - 1][j - 2];
    } /*of j*/
  } /*of i*/
  lp->A[m - 2][0] = 1;
  lp->A[m - 2][d - 1] = -1;
  /* make the LP bounded.  */
  lp->A[m - 1][d - 1] = 1;
  /* objective is to maximize z.  */

  if (localdebug) {
    fprintf(stdout, "dd_CreateLP_V_ImplicitLinearity: an new lp is\n");
    dd_WriteLP(std::cout, lp);
  }

  return lp;
}

template <typename T>
void dd_CreateLP_H_Redundancy(dd_matrixdata<T> *M, dd_rowrange itest,
                              dd_lpdata<T> *lp) {
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;

  linc = set_card(M->linset);
  m = M->rowsize + 1 + linc;
  dd_LPData_reset_m(m, lp);
  /* We represent each equation by two inequalities.
     This is not the best way but makes the code simple. */
  d = M->colsize;

  lp->objective = dd_LPmin;
  lp->redcheck_extensive = false; /* this is default */

  irev = M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    if (set_member(i, M->linset)) {
      irev++;
      set_addelem(lp->equalityset, i); /* it is equality. */
      for (j = 1; j <= d; j++)
        lp->A[irev - 1][j - 1] = -M->matrix[i - 1][j - 1];
    }
    for (j = 1; j <= d; j++)
      lp->A[i - 1][j - 1] = M->matrix[i - 1][j - 1];
  }
  for (j = 1; j <= d; j++)
    lp->A[m - 1][j - 1] = M->matrix[itest - 1][j - 1];
  lp->A[itest - 1][0] += T(1); /* relax the original inequality by one */
}

template <typename T>
void dd_CreateLP_V_Redundancy(dd_matrixdata<T> *M, dd_rowrange itest,
                              dd_lpdata<T> *lp) {
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;

  linc = set_card(M->linset);
  m = M->rowsize + 1 + linc;
  dd_LPData_reset_m(m, lp);
  /* We represent each equation by two inequalities.
     This is not the best way but makes the code simple. */
  d = M->colsize + 1;
  /* One more column.  This is different from the H-reprentation case */

  /* The below must be modified for V-representation!!!  */

  lp->objective = dd_LPmin;
  lp->redcheck_extensive = false; /* this is default */

  irev = M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    if (i == itest)
      lp->A[i - 1][0] =
          1; /* this is to make the LP bounded, ie. the min >= -1 */
    else
      lp->A[i - 1][0] = 0; /* It is almost completely degerate LP */
    if (set_member(i, M->linset)) {
      irev++;
      set_addelem(lp->equalityset, i); /* it is equality. */
      /* the reversed row irev is not in the equality set. */
      for (j = 2; j <= d; j++)
        lp->A[irev - 1][j - 1] = -M->matrix[i - 1][j - 2];
    }
    for (j = 2; j <= d; j++)
      lp->A[i - 1][j - 1] = M->matrix[i - 1][j - 2];
  }
  /* objective is to violate the inequality in question.  */
  for (j = 2; j <= d; j++)
    lp->A[m - 1][j - 1] = M->matrix[itest - 1][j - 2];
  lp->A[m - 1][0] = 0; /* the constant term for the objective is zero */
}

template <typename T>
dd_lpdata<T> *dd_CreateLP_V_SRedundancy(dd_matrixdata<T> *M,
                                        dd_rowrange itest) {
  /*
       V-representation (=boundary problem)
         g* = maximize
           1^T b_{I-itest} x_0 + 1^T A_{I-itest}    (the sum of slacks)
         subject to
           b_itest x_0     + A_itest x      =  0 (the point has to lie on the
     boundary) b_{I-itest} x_0 + A_{I-itest} x >=  0 (all nonlinearity
     generators in one side) 1^T b_{I-itest} x_0 + 1^T A_{I-itest} x <=  1 (to
     make an LP bounded) b_L x_0         + A_L x = 0.  (linearity generators)

      The redundant row is strongly redundant if and only if g* is zero.
  */

  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_lpdata<T> *lp;

  linc = set_card(M->linset);
  m = M->rowsize + 1 + linc + 2;
  /* We represent each equation by two inequalities.
     This is not the best way but makes the code simple.
     Two extra constraints are for the first equation and the bouding
     inequality.
     */
  d = M->colsize + 1;
  /* One more column.  This is different from the H-reprentation case */

  /* The below must be modified for V-representation!!!  */

  lp = dd_CreateLPData<T>(m, d);
  lp->objective = dd_LPmax;

  irev = M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    if (i == itest)
      lp->A[i - 1][0] = 0; /* this is a half of the boundary constraint. */
    else
      lp->A[i - 1][0] = 0; /* It is almost completely degerate LP */
    if (set_member(i, M->linset) || i == itest) {
      irev++;
      set_addelem(lp->equalityset, i); /* it is equality. */
      /* the reversed row irev is not in the equality set. */
      for (j = 2; j <= d; j++)
        lp->A[irev - 1][j - 1] = -M->matrix[i - 1][j - 2];
    }
    for (j = 2; j <= d; j++) {
      lp->A[i - 1][j - 1] = M->matrix[i - 1][j - 2];
      lp->A[m - 1][j - 1] +=
          lp->A[i - 1][j - 1]; /* the objective is the sum of all ineqalities */
    }
  }
  for (j = 2; j <= d; j++) {
    lp->A[m - 2][j - 1] = -lp->A[m - 1][j - 1];
    /* to make an LP bounded.  */
  }
  lp->A[m - 2][0] = 1; /* the constant term for the bounding constraint is 1 */

  return lp;
}

template <typename T>
bool dd_Redundant(dd_matrixdata<T> *M, dd_rowrange itest, T *certificate,
                  dd_ErrorType *error, data_temp_simplex<T> *data, size_t const& maxiter, std::ostream& os) {
  /* Checks whether the row itest is redundant for the representation.
     All linearity rows are not checked and considered NONredundant.
     This code works for both H- and V-representations.  A certificate is
     given in the case of non-redundancy, showing a solution x violating only
    the itest inequality for H-representation, a hyperplane RHS and normal (x_0,
    x) that separates the itest from the rest.  More explicitly, the LP to be
    setup is

     H-representation
       f* = minimize
         b_itest     + A_itest x
       subject to
         b_itest + 1 + A_itest x     >= 0 (relaxed inequality to make an LP
    bounded) b_{I-itest} + A_{I-itest} x >= 0 (all inequalities except for
    itest) b_L         + A_L x = 0.  (linearity)

     V-representation (=separation problem)
       f* = minimize
         b_itest x_0     + A_itest x
       subject to
         b_itest x_0     + A_itest x     >= -1 (to make an LP bounded)
         b_{I-itest} x_0 + A_{I-itest} x >=  0 (all nonlinearity generators
    except for itest in one side) b_L x_0         + A_L x = 0.  (linearity
    generators)

    Here, the input matrix is considered as (b, A), i.e. b corresponds to the
    first column of input and the row indices of input is partitioned into I and
    L where L is the set of linearity. In both cases, the itest data is
    nonredundant if and only if the optimal value f* is negative. The
    certificate has dimension one more for V-representation case.
  */

  dd_colrange j;
  dd_lpdata<T> *lp;
  dd_ErrorType err = dd_NoError;
  bool answer = false, localdebug = false;

  *error = dd_NoError;
  if (set_member(itest, M->linset)) {
    if (localdebug)
      printf(
          "The %ld th row is linearity and redundancy checking is skipped.\n",
          itest);
    return answer;
  }

  /* Create an LP data for redundancy checking */
  lp = dd_CreateLPData_from_M<T>(M);
  if (M->representation == dd_Generator) {
    dd_CreateLP_V_Redundancy(M, itest, lp);
  } else {
    dd_CreateLP_H_Redundancy(M, itest, lp);
  }

  dd_LPSolve_data(lp, dd_choiceRedcheckAlgorithm, &err, data, maxiter, os);
  if (err != dd_NoError) {
    *error = err;
  } else {
    for (j = 0; j < lp->d; j++)
      certificate[j] = lp->sol[j];

    if (lp->optvalue < 0) {
      answer = false;
      if (localdebug)
        fprintf(stdout, "==> %ld th row is nonredundant.\n", itest);
    } else {
      answer = true;
      if (localdebug)
        fprintf(stdout, "==> %ld th row is redundant.\n", itest);
    }
  }
  dd_FreeLPData(lp);
  return answer;
}

template <typename T>
bool dd_RedundantExtensive(dd_matrixdata<T> *M, dd_rowrange itest,
                           T *certificate, dd_rowset *redset,
                           dd_ErrorType *error, size_t const& maxiter, std::ostream& os) {
  /* This uses the same LP construction as dd_Redundant.  But, while it is
     checking the redundancy of itest, it also tries to find some other variable
     that are redundant (i.e. forced to be nonnegative).  This is expensive as
     it used the complete tableau information at each DualSimplex pivot.  The
     redset must be initialized before this function is called.
  */

  dd_colrange j;
  dd_lpdata<T> *lp;
  dd_ErrorType err = dd_NoError;
  bool answer = false, localdebug = false;

  *error = dd_NoError;
  if (set_member(itest, M->linset)) {
    if (localdebug)
      printf(
          "The %ld th row is linearity and redundancy checking is skipped.\n",
          itest);
    return answer;
  }

  /* Create an LP data for redundancy checking */
  lp = dd_CreateLPData_from_M<T>(M);
  if (M->representation == dd_Generator) {
    dd_CreateLP_V_Redundancy(M, itest, lp);
  } else {
    dd_CreateLP_H_Redundancy(M, itest, lp);
  }

  lp->redcheck_extensive = true;

  dd_LPSolve0(lp, dd_DualSimplex, &err, maxiter, os);
  if (err != dd_NoError) {
    *error = err;
  } else {
    set_copy(*redset, lp->redset_extra);
    set_delelem(*redset, itest);
    /* itest row might be redundant in the lp but this has nothing to do with
    its redundancy in the original system M.   Thus we must delete it.  */
    for (j = 0; j < lp->d; j++) {
      certificate[j] = lp->sol[j];
    }

    if (lp->optvalue < 0) {
      answer = false;
      if (localdebug)
        fprintf(stdout, "==> %ld th row is nonredundant.\n", itest);
    } else {
      answer = true;
      if (localdebug)
        fprintf(stdout, "==> %ld th row is redundant.\n", itest);
    }
  }
  dd_FreeLPData(lp);
  return answer;
}

template <typename T>
dd_rowset dd_RedundantRows(dd_matrixdata<T> *M, dd_ErrorType *error, size_t const& maxiter, std::ostream& os) {
  dd_rowrange i, m;
  dd_rowset redset;
  dd_matrixdata<T> *Mcopy;
  T *cvec; /* certificate */
  bool localdebug = false;

  m = M->rowsize;
  dd_colrange d = get_d_size(M);

  Mcopy = dd_MatrixCopy(M);
  dd_AllocateArow(d, &cvec);
  set_initialize(&redset, m);
  data_temp_simplex<T> *data = allocate_data_simplex<T>(get_m_size(M), d);
  for (i = m; i >= 1; i--) {
    if (dd_Redundant(Mcopy, i, cvec, error, data, maxiter, os)) {
      if (localdebug)
        printf("dd_RedundantRows: the row %ld is redundant.\n", i);
      set_addelem(redset, i);
      dd_MatrixRowRemove(&Mcopy, i);
    } else {
      if (localdebug)
        printf("dd_RedundantRows: the row %ld is essential.\n", i);
    }
    if (*error != dd_NoError)
      goto _L99;
  }
_L99:
  free_data_simplex(data);
  dd_FreeMatrix(Mcopy);
  dd_FreeArow(cvec);
  return redset;
}

template <typename T>
bool dd_MatrixRedundancyRemove(dd_matrixdata<T> **M, dd_rowset *redset,
                               dd_rowindex *newpos, dd_ErrorType *error, size_t const& maxiter, std::ostream& os) {
  /* It returns the set of all redundant rows.  This should be called after all
     implicit linearity are recognized with dd_MatrixCanonicalizeLinearity.
  */

  dd_rowrange i, k, m, m1;
  dd_colrange d;
  dd_rowset redset1;
  dd_rowindex newpos1;
  dd_matrixdata<T> *M1 = nullptr;
  T *cvec; /* certificate */
  bool success = false;
  bool localdebug = false;

  m = (*M)->rowsize;
  set_initialize(redset, m);
  M1 = dd_MatrixSortedUniqueCopy(*M, newpos);
  for (i = 1; i <= m; i++) {
    if ((*newpos)[i] <= 0)
      set_addelem(*redset, i);
    if (localdebug)
      printf(" %ld:%ld", i, (*newpos)[i]);
  }
  if (localdebug)
    printf("\n");

  if ((*M)->representation == dd_Generator) {
    d = (*M)->colsize + 1;
  } else {
    d = (*M)->colsize;
  }
  m1 = M1->rowsize;
  if (localdebug) {
    fprintf(stdout,
            "dd_MatrixRedundancyRemove: By sorting, %ld rows have been "
            "removed.  The remaining has %ld rows.\n",
            m - m1, m1);
    /* dd_WriteMatrix(stdout,M1);  */
  }
  dd_AllocateArow(d, &cvec);
  set_initialize(&redset1, M1->rowsize);
  k = 1;
  do {
    if (dd_RedundantExtensive(M1, k, cvec, &redset1, error, maxiter, os)) {
      set_addelem(redset1, k);
      dd_MatrixRowsRemove2(&M1, redset1, &newpos1);
      for (i = 1; i <= m; i++) {
        if ((*newpos)[i] > 0) {
          if (set_member((*newpos)[i], redset1)) {
            set_addelem(*redset, i);
            (*newpos)[i] = 0; /* now the original row i is recognized redundant
                                 and removed from M1 */
          } else {
            (*newpos)[i] =
                newpos1[(*newpos)[i]]; /* update the new pos vector */
          }
        }
      }
      set_free(redset1);
      set_initialize(&redset1, M1->rowsize);
      if (localdebug) {
        printf("dd_MatrixRedundancyRemove: the row %ld is redundant. The new "
               "matrix has %ld rows.\n",
               k, M1->rowsize);
        /* dd_WriteMatrix(stdout, M1);  */
      }
      delete[] newpos1;
    } else {
      if (set_card(redset1) > 0) {
        dd_MatrixRowsRemove2(&M1, redset1, &newpos1);
        for (i = 1; i <= m; i++) {
          if ((*newpos)[i] > 0) {
            if (set_member((*newpos)[i], redset1)) {
              set_addelem(*redset, i);
              (*newpos)[i] = 0; /* now the original row i is recognized
                                   redundant and removed from M1 */
            } else {
              (*newpos)[i] =
                  newpos1[(*newpos)[i]]; /* update the new pos vector */
            }
          }
        }
        set_free(redset1);
        set_initialize(&redset1, M1->rowsize);
        delete[] newpos1;
      }
      if (localdebug) {
        printf("dd_MatrixRedundancyRemove: the row %ld is essential. The new "
               "matrix has %ld rows.\n",
               k, M1->rowsize);
        /* dd_WriteMatrix(stdout, M1);  */
      }
      k++;
    }
    if (*error != dd_NoError)
      goto _L99;
  } while (k <= M1->rowsize);
  if (localdebug)
    dd_WriteMatrix(stdout, M1);
  success = true;

_L99:
  dd_FreeMatrix(*M);
  *M = M1;
  dd_FreeArow(cvec);
  set_free(redset1);
  return success;
}

template <typename T>
bool dd_SRedundant(dd_matrixdata<T> *M, dd_rowrange itest, T *certificate,
                   dd_ErrorType *error, size_t const& maxiter, std::ostream& os) {
  /* Checks whether the row itest is strongly redundant for the representation.
     A row is strongly redundant in H-representation if every point in
     the polyhedron satisfies it with strict inequality.
     A row is strongly redundant in V-representation if this point is in
     the interior of the polyhedron.

     All linearity rows are not checked and considered NOT strongly redundant.
     This code works for both H- and V-representations.  A certificate is
     given in the case of non-redundancy, showing a solution x violating only
    the itest inequality for H-representation, a hyperplane RHS and normal (x_0,
    x) that separates the itest from the rest.  More explicitly, the LP to be
    setup is

     H-representation
       f* = minimize
         b_itest     + A_itest x
       subject to
         b_itest + 1 + A_itest x     >= 0 (relaxed inequality to make an LP
    bounded) b_{I-itest} + A_{I-itest} x >= 0 (all inequalities except for
    itest) b_L         + A_L x = 0.  (linearity)

     V-representation (=separation problem)
       f* = minimize
         b_itest x_0     + A_itest x
       subject to
         b_itest x_0     + A_itest x     >= -1 (to make an LP bounded)
         b_{I-itest} x_0 + A_{I-itest} x >=  0 (all nonlinearity generators
    except for itest in one side) b_L x_0         + A_L x = 0.  (linearity
    generators)

    Here, the input matrix is considered as (b, A), i.e. b corresponds to the
    first column of input and the row indices of input is partitioned into I and
    L where L is the set of linearity. In H-representation, the itest data is
    strongly redundant if and only if the optimal value f* is positive. In
    V-representation, the itest data is redundant if and only if the optimal
    value f* is zero (as the LP is homogeneous and the optimal value is always
    non-positive).  To recognize strong redundancy, one can set up a second LP

     V-representation (=boundary problem)
       g* = maximize
         1^T b_{I-itest} x_0 + 1^T A_{I-itest}    (the sum of slacks)
       subject to
         b_itest x_0     + A_itest x      =  0 (the point has to lie on the
    boundary) b_{I-itest} x_0 + A_{I-itest} x >=  0 (all nonlinearity generators
    in one side) 1^T b_{I-itest} x_0 + 1^T A_{I-itest} x <=  1 (to make an LP
    bounded) b_L x_0         + A_L x = 0.  (linearity generators)

    The redundant row is strongly redundant if and only if g* is zero.

    The certificate has dimension one more for V-representation case.
  */

  dd_colrange j;
  dd_lpdata<T> *lp;
  dd_ErrorType err = dd_NoError;
  bool answer = false;
  bool localdebug = false;

  *error = dd_NoError;
  if (set_member(itest, M->linset)) {
    if (localdebug)
      printf("The %ld th row is linearity and strong redundancy checking is "
             "skipped.\n",
             itest);
    goto _L99;
  }

  /* Create an LP data for redundancy checking */
  lp = dd_CreateLPData_from_M<T>(M);
  if (M->representation == dd_Generator) {
    dd_CreateLP_V_Redundancy(M, itest, lp);
  } else {
    dd_CreateLP_H_Redundancy(M, itest, lp);
  }

  dd_LPSolve(lp, dd_choiceRedcheckAlgorithm, &err, maxiter, os);
  if (err != dd_NoError) {
    *error = err;
    goto _L999;
  } else {

    for (j = 0; j < lp->d; j++) {
      certificate[j] = lp->sol[j];
    }

    if (M->representation == dd_Inequality) {
      if (lp->optvalue > 0) {
        answer = true;
        if (localdebug)
          fprintf(stdout, "==> %ld th inequality is strongly redundant.\n",
                  itest);
      } else {
        answer = false;
        if (localdebug)
          fprintf(stdout, "==> %ld th inequality is not strongly redundant.\n",
                  itest);
      }
    } else {
      if (lp->optvalue < 0) {
        answer = false;
        if (localdebug)
          fprintf(stdout, "==> %ld th point is not strongly redundant.\n",
                  itest);
      } else {
        /* for V-representation, we have to solve another LP */
        dd_FreeLPData(lp);
        lp = dd_CreateLP_V_SRedundancy(M, itest);
        dd_LPSolve(lp, dd_DualSimplex, &err, maxiter, os);
        if (localdebug)
          dd_WriteLPResult(std::cout, lp, err);
        if (lp->optvalue > 0) {
          answer = false;
          if (localdebug)
            fprintf(stdout, "==> %ld th point is not strongly redundant.\n",
                    itest);
        } else {
          answer = true;
          if (localdebug)
            fprintf(stdout, "==> %ld th point is strongly redundant.\n", itest);
        }
      }
    }
  }
_L999:
  dd_FreeLPData(lp);
_L99:
  return answer;
}

template <typename T>
dd_rowset dd_SRedundantRows(dd_matrixdata<T> *M, dd_ErrorType *error, size_t const& maxiter, std::ostream& os) {
  dd_rowrange i, m;
  dd_colrange d;
  dd_rowset redset;
  dd_matrixdata<T> *Mcopy;
  T *cvec; /* certificate */
  bool localdebug = false;

  m = M->rowsize;
  if (M->representation == dd_Generator) {
    d = M->colsize + 1;
  } else {
    d = M->colsize;
  }
  Mcopy = dd_MatrixCopy(M);
  dd_AllocateArow(d, &cvec);
  set_initialize(&redset, m);
  for (i = m; i >= 1; i--) {
    if (dd_SRedundant(Mcopy, i, cvec, error, maxiter, os)) {
      if (localdebug)
        printf("dd_SRedundantRows: the row %ld is strongly redundant.\n", i);
      set_addelem(redset, i);
      dd_MatrixRowRemove(&Mcopy, i);
    } else {
      if (localdebug)
        printf("dd_SRedundantRows: the row %ld is not strongly redundant.\n",
               i);
    }
    if (*error != dd_NoError)
      goto _L99;
  }
_L99:
  dd_FreeMatrix(Mcopy);
  dd_FreeArow(cvec);
  return redset;
}

template <typename T>
dd_rowset dd_RedundantRowsViaShooting(dd_matrixdata<T> *M, dd_ErrorType *error,
                                      size_t const& maxiter, std::ostream &os) {
#ifdef DEBUG_CDD
  os << "CDD: dd_RedundantRowsViaShooting, step 1 *error=" << *error << "\n";
#endif
  /*
     For H-representation only and not quite reliable,
     especially when floating-point arithmetic is used.
     Use the ordinary (slower) method dd_RedundantRows.
  */

  dd_rowrange i, m, ired;
  dd_colrange j, k, d;
  T *shootdir;
  dd_LPSolverType solver = dd_DualSimplex;
  bool localdebug = false;

  m = M->rowsize;
  d = M->colsize;
  dd_rowset redset;
  set_initialize(&redset, m);
  dd_AllocateArow(d, &shootdir);
  if (localdebug) {
    std::cout << "ViaShooting : M->colsize=" << M->colsize << "\n";
  }
#ifdef DEBUG_CDD
  os << "CDD: dd_RedundantRowsViaShooting, step 2 *error=" << *error << "\n";
#endif

  if (set_card(M->linset) != 0) {
    std::cerr << "This code works only in the absence of linearity relations\n";
    *error = dd_NonZeroLinearity;
#ifdef DEBUG_CDD
    os << "CDD: Error dd_NonZeroLinearity\n";
#endif
    return redset;
  }
#ifdef DEBUG_CDD
  os << "CDD: dd_RedundantRowsViaShooting, step 3 *error=" << *error << "\n";
#endif

  dd_lpdata<T> *lpw = dd_CreateLPData_from_M<T>(M);
  lpw->objective = dd_LPmin;
  lpw->redcheck_extensive = false; /* this is default */
  dd_LPData_reset_m(1, lpw);

  dd_ErrorType err = dd_NoError;
  data_temp_simplex<T> *data =
      allocate_data_simplex<T>(get_m_size(M), get_d_size(M));
#ifdef DEBUG_CDD
  os << "CDD: dd_RedundantRowsViaShooting, step 4 *error=" << *error << "\n";
#endif
  auto dd_Redundant_loc = [&]() -> bool {
    dd_colrange j;
    dd_rowrange mi = lpw->m;
    for (j = 0; j < d; j++)
      lpw->A[mi - 1][j] = lpw->A[mi - 2][j];
    lpw->A[mi - 2][0] += T(1);
    if (localdebug) {
      std::cout << "Hyperplane case\n";
      std::cout << "dd_Redundant_loc: lpw->m=" << lpw->m << " lpw=\n";
      dd_WriteLP(std::cout, lpw);
    }
    dd_LPSolve_data(lpw, dd_choiceRedcheckAlgorithm, &err, data, maxiter, os);
    lpw->A[mi - 2][0] -= T(1);
    if (lpw->optvalue < 0)
      return false;
    else
      return true;
  };
  auto set_entry_in_lpw = [&](dd_rowrange irow) -> void {
    dd_rowrange mi = lpw->m;
    for (k = 0; k < d; k++)
      lpw->A[mi - 2][k] = M->matrix[irow - 1][k];
  };
  auto insert_entry_in_lpw = [&](dd_rowrange irow) -> void {
    dd_rowrange mi = lpw->m;
    for (k = 0; k < d; k++)
      lpw->A[mi - 1][k] = M->matrix[irow - 1][k];
    mi++;
    dd_LPData_reset_m(mi, lpw);
  };
  auto decrement_entry_in_lpw = [&]() -> void {
    dd_rowrange mi = lpw->m;
    mi--;
    dd_LPData_reset_m(mi, lpw);
  };

  /* Whether we have reached a conclusion in any way on the code */
  dd_rowset is_decided;
  set_initialize(&is_decided, m);
#ifdef DEBUG_CDD
  os << "CDD: dd_RedundantRowsViaShooting, step 5 *error=" << *error << "\n";
#endif

  /* First find some (likely) nonredundant inequalities by Interior Point Find.
   */
  dd_lpdata<T> *lp0 = dd_Matrix2LP(M);
  dd_lpdata<T> *lp = dd_MakeLPforInteriorFinding(lp0);
  dd_FreeLPData(lp0);
  dd_LPSolve(lp, solver, &err, maxiter, os);
  if (localdebug) {
    std::cout << "lp->sol=";
    dd_WriteT(std::cout, lp->sol, d);
  }
#ifdef DEBUG_CDD
  os << "CDD: dd_RedundantRowsViaShooting, step 6 *error=" << *error << "\n";
#endif

  if (lp->optvalue > 0) {
    if (localdebug)
      std::cout << "dd_Positive=T case\n";
    /* An interior point is found.  Use rayshooting to find some nonredundant
       inequalities. */
    for (k = 0; k < d; k++)
      shootdir[k] = 0;
    for (j = 1; j < d; j++) {
      shootdir[j] = 1; /* j-th unit vector */
      if (localdebug)
        dd_WriteT(std::cout, shootdir, d);
      ired = dd_RayShooting(M, lp->sol, shootdir);
      if (localdebug)
        printf("nonredundant row %3ld found by shooting.\n", ired);
      if (ired > 0 && !set_member(ired, is_decided)) {
        set_addelem(is_decided, ired);
        insert_entry_in_lpw(ired);
      }
      shootdir[j] = -1; /* negative of the j-th unit vector */
      ired = dd_RayShooting(M, lp->sol, shootdir);
      if (localdebug)
        printf("nonredundant row %3ld found by shooting.\n", ired);
      if (ired > 0 && !set_member(ired, is_decided)) {
        set_addelem(is_decided, ired);
        insert_entry_in_lpw(ired);
      }
      shootdir[j] = 0; /* restore to 0 */
    }

    if (localdebug) {
      printf("The initial nonredundant set is:");
      for (i = 1; i <= m; i++)
        if (set_member(i, is_decided))
          printf(" %ld", i);
      printf("\n");
    }

    i = 1;
    while (i <= m) {
      if (localdebug)
        std::cout << "i=" << i << " is_decided=" << set_member(i, is_decided)
                  << "\n";
      if (!set_member(i,
                      is_decided)) { /* the ith inequality is not yet checked */
        if (localdebug)
          std::cout << "Checking redundancy of " << i << " th inequality\n";
        insert_entry_in_lpw(i);
        if (localdebug) {
          std::cout << "M->matrix[i-1]=";
          dd_WriteT(std::cout, M->matrix[i - 1], d);
        }
        if (!dd_Redundant_loc()) {
          for (k = 0; k < d; k++)
            shootdir[k] = lpw->sol[k] - lp->sol[k];
          if (localdebug) {
            std::cout << "shootdir=";
            dd_WriteT(std::cout, shootdir, d);
          }
          ired = dd_RayShooting(M, lp->sol, shootdir);
          set_addelem(is_decided, ired);
          set_entry_in_lpw(ired);
          if (localdebug) {
            fprintf(stdout,
                    "The %ld th inequality is nonredundant for the subsystem\n",
                    i);
            fprintf(stdout,
                    "The nonredundancy of %ld th inequality is found by "
                    "shooting.\n",
                    ired);
            dd_WriteT(std::cout, M->matrix[ired - 1], d);
          }
        } else {
          if (localdebug)
            fprintf(stdout,
                    "The %ld th inequality is redundant for the subsystem and "
                    "thus for the whole.\n",
                    i);
          set_addelem(is_decided, i);
          decrement_entry_in_lpw();
          set_addelem(redset, i);
          i++;
        }
      } else {
        if (localdebug)
          std::cout << "Case already decided\n";
        i++;
      }
    } /* endwhile */
  } else {
    if (localdebug)
      std::cout << "dd_Positive=F case\n";
    /* No interior point is found.  Apply the standard LP technique.  */
    set_free(redset);
    redset = dd_RedundantRows(M, error, maxiter, os);
#ifdef DEBUG_CDD
    if (*error != dd_NoError) {
      os << "CDD: Error in dd_RedundantRows\n";
    }
#endif
  }
#ifdef DEBUG_CDD
  os << "CDD: dd_RedundantRowsViaShooting, step 7 *error=" << *error << "\n";
#endif
  free_data_simplex(data);

  dd_FreeLPData(lp);
  dd_FreeLPData(lpw);

  dd_FreeArow(shootdir);
  set_free(is_decided);
  return redset;
}

template <typename T>
dd_rowset
dd_RedundantRowsViaShootingBlocks(dd_matrixdata<T> *M, dd_ErrorType *error,
                                  std::vector<int> const &BlockBelong, size_t const& maxiter, std::ostream& os) {
  dd_rowrange i, m, ired;
  dd_colrange j, k, d;
  T *shootdir;
  dd_LPSolverType solver = dd_DualSimplex;
  bool localdebug = true;
  //  bool localdebug = false;

  m = M->rowsize;
  d = M->colsize;
  std::cerr << "m=" << m << " d=" << d << "\n";
  dd_rowset redset;
  set_initialize(&redset, m);
  dd_AllocateArow(d, &shootdir);
  if (localdebug) {
    std::cout << "ViaShooting : M->colsize=" << M->colsize << "\n";
  }

  if (set_card(M->linset) != 0) {
    std::cerr << "This code works only in the absence of linearity relations\n";
    *error = dd_NonZeroLinearity;
    return redset;
  }

  dd_lpdata<T> *lpw = dd_CreateLPData_from_M<T>(M);
  lpw->objective = dd_LPmin;
  lpw->redcheck_extensive = false; /* this is default */
  dd_LPData_reset_m(1, lpw);

  dd_ErrorType err = dd_NoError;
  data_temp_simplex<T> *data =
      allocate_data_simplex<T>(get_m_size(M), get_d_size(M));
  auto dd_Redundant_loc = [&]() -> bool {
    dd_colrange j;
    dd_rowrange mi = lpw->m;
    for (j = 0; j < d; j++)
      lpw->A[mi - 1][j] = lpw->A[mi - 2][j];
    lpw->A[mi - 2][0] += T(1);
    if (localdebug) {
      std::cout << "Hyperplane case\n";
      std::cout << "dd_Redundant_loc: lpw->m=" << lpw->m << " lpw=\n";
      dd_WriteLP(std::cout, lpw);
    }
    dd_LPSolve_data(lpw, dd_choiceRedcheckAlgorithm, &err, data, maxiter, os);
    lpw->A[mi - 2][0] -= T(1);
    return lpw->optvalue >= 0;
  };
  auto set_entry_in_lpw = [&](dd_rowrange irow) -> void {
    dd_rowrange mi = lpw->m;
    for (k = 0; k < d; k++)
      lpw->A[mi - 2][k] = M->matrix[irow - 1][k];
  };
  auto insert_entry_in_lpw = [&](dd_rowrange irow) -> void {
    dd_rowrange mi = lpw->m;
    for (k = 0; k < d; k++)
      lpw->A[mi - 1][k] = M->matrix[irow - 1][k];
    mi++;
    dd_LPData_reset_m(mi, lpw);
  };
  auto decrement_entry_in_lpw = [&]() -> void {
    dd_rowrange mi = lpw->m;
    mi--;
    dd_LPData_reset_m(mi, lpw);
  };

  /* Whether we have reached a conclusion in any way on the code */
  dd_rowset is_decided;
  set_initialize(&is_decided, m);
  int e_max = 0;
  for (auto &e_val : BlockBelong) {
    if (e_val > e_max)
      e_max = e_val;
  }
  size_t n_block = e_max + 1;
  std::vector<std::vector<dd_rowrange>> list_blocks(n_block);
  for (size_t i = 0; i < BlockBelong.size(); i++) {
    int iBlock = BlockBelong[i];
    dd_rowrange iredw = i + 1;
    list_blocks[iBlock].push_back(iredw);
  }
  auto get_block = [&](dd_rowrange const &pos) -> std::vector<dd_rowrange> {
    int iBlock = BlockBelong[pos - 1];
    return list_blocks[iBlock];
  };

  /* First find some (likely) nonredundant inequalities by Interior Point Find.
   */
  dd_lpdata<T> *lp0 = dd_Matrix2LP(M);
  dd_lpdata<T> *lp = dd_MakeLPforInteriorFinding(lp0);
  dd_FreeLPData(lp0);
  dd_LPSolve(lp, solver, &err, maxiter, os);
  if (localdebug) {
    std::cout << "lp->sol=";
    dd_WriteT(std::cout, lp->sol, d);
  }

  if (lp->optvalue > 0) {
    if (localdebug)
      std::cout << "dd_Positive=T case\n";
    /* An interior point is found.  Use rayshooting to find some nonredundant
       inequalities. */
    for (k = 0; k < d; k++)
      shootdir[k] = 0;
    for (j = 1; j < d; j++) {
      shootdir[j] = 1; /* j-th unit vector */
      if (localdebug)
        dd_WriteT(std::cout, shootdir, d);
      ired = dd_RayShooting(M, lp->sol, shootdir);
      if (localdebug)
        std::cout << "nonredundant row " << ired << " found by shooting\n";
      if (ired > 0 && !set_member(ired, is_decided)) {
        std::vector<dd_rowrange> eBlock = get_block(ired);
        std::cout << "ired=" << ired << " |eBlock|=" << eBlock.size() << "\n";
        for (auto &jred : eBlock) {
          set_addelem(is_decided, jred);
          std::cout << "1: Deciding " << jred << "\n";
          insert_entry_in_lpw(jred);
        }
      }
      shootdir[j] = -1; /* negative of the j-th unit vector */
      ired = dd_RayShooting(M, lp->sol, shootdir);
      if (localdebug)
        printf("nonredundant row %3ld found by shooting.\n", ired);
      if (ired > 0 && !set_member(ired, is_decided)) {
        std::vector<dd_rowrange> eBlock = get_block(ired);
        std::cout << "ired=" << ired << " |eBlock|=" << eBlock.size() << "\n";
        for (auto &jred : eBlock) {
          set_addelem(is_decided, jred);
          std::cout << "2: Deciding " << jred << "\n";
          insert_entry_in_lpw(jred);
        }
      }
      shootdir[j] = 0; /* restore to 0 */
    }

    if (localdebug) {
      printf("The initial nonredundant set is:");
      for (i = 1; i <= m; i++)
        if (set_member(i, is_decided))
          printf(" %ld", i);
      printf("\n");
    }

    i = 1;
    while (i <= m) {
      if (localdebug)
        std::cout << "i=" << i << " is_decided=" << set_member(i, is_decided)
                  << "\n";
      if (!set_member(i,
                      is_decided)) { /* the ith inequality is not yet checked */
        if (localdebug)
          std::cout << "Checking redundancy of " << i << " th inequality\n";
        insert_entry_in_lpw(i);
        if (localdebug) {
          std::cout << "M->matrix[i-1]=";
          dd_WriteT(std::cout, M->matrix[i - 1], d);
        }
        if (!dd_Redundant_loc()) {
          for (k = 0; k < d; k++)
            shootdir[k] = lpw->sol[k] - lp->sol[k];
          if (localdebug) {
            std::cout << "shootdir=";
            dd_WriteT(std::cout, shootdir, d);
          }
          ired = dd_RayShooting(M, lp->sol, shootdir);
          if (!set_member(ired, is_decided)) {
            std::vector<dd_rowrange> eBlock = get_block(ired);
            std::cout << "ired=" << ired << " |eBlock|=" << eBlock.size()
                      << "\n";
            for (auto &jred : eBlock) {
              set_addelem(is_decided, jred);
              std::cout << "3: Deciding " << jred << "\n";
              set_entry_in_lpw(jred);
            }
          }
          if (localdebug) {
            std::cout << "The " << i
                      << " inequality is nonredundant for the subsystem\n";
            std::cout << "The nonredundancy of " << ired
                      << "% inequality is found by shooting\n";
            dd_WriteT(std::cout, M->matrix[ired - 1], d);
          }
        } else {
          if (localdebug)
            std::cout << "The " << i
                      << " inequality is redundant for the subsystem and thus "
                         "for the whole\n";
          decrement_entry_in_lpw();
          if (!set_member(i, is_decided)) {
            std::vector<dd_rowrange> eBlock = get_block(i);
            std::cout << "i=" << i << " |eBlock|=" << eBlock.size() << "\n";
            for (auto &jred : eBlock) {
              set_addelem(is_decided, jred);
              std::cout << "4: Deciding " << jred << "\n";
              set_addelem(redset, jred);
            }
          }
          i++;
        }
      } else {
        if (localdebug)
          std::cout << "Case already decided\n";
        i++;
      }
    } /* endwhile */
  } else {
    if (localdebug)
      std::cout << "dd_Positive=F case\n";
    /* No interior point is found.  Apply the standard LP technique.  */
    set_free(redset);
    redset = dd_RedundantRows(M, error, maxiter, os);
  }
  free_data_simplex(data);

  dd_FreeLPData(lp);
  dd_FreeLPData(lpw);

  dd_FreeArow(shootdir);
  return redset;
}

template <typename T>
dd_setfamily *dd_Matrix2Adjacency(dd_matrixdata<T> *M, dd_ErrorType *error, size_t const& maxiter, std::ostream& os) {
  /* This is to generate the (facet) graph of a polyheron (H) V-represented by M
     using LPs. Since it does not use the representation conversion, it should
     work for a large scale problem.
  */
  dd_rowrange i, m;
  dd_colrange d;
  dd_rowset redset;
  dd_matrixdata<T> *Mcopy;
  dd_setfamily *F = nullptr;

  m = M->rowsize;
  d = M->colsize;
  if (m <= 0 || d <= 0) {
    *error = dd_EmptyRepresentation;
    goto _L999;
  }
  Mcopy = dd_MatrixCopy(M);
  F = dd_CreateSetFamily(m, m);
  for (i = 1; i <= m; i++) {
    if (!set_member(i, M->linset)) {
      set_addelem(Mcopy->linset, i);
      redset = dd_RedundantRows(Mcopy, error, maxiter, os);
      /* redset should contain all nonadjacent ones */
      set_uni(redset, redset,
              Mcopy->linset); /* all linearity elements should be nonadjacent */
      set_compl(F->set[i - 1], redset); /* set the adjacency list of vertex i */
      set_delelem(Mcopy->linset, i);
      set_free(redset);
      if (*error != dd_NoError)
        goto _L99;
    }
  }
_L99:
  dd_FreeMatrix(Mcopy);
_L999:
  return F;
}

template <typename T>
dd_setfamily *dd_Matrix2WeakAdjacency(dd_matrixdata<T> *M,
                                      dd_ErrorType *error, size_t const& maxiter, std::ostream& os) {
  /* This is to generate the weak-adjacency (facet) graph of a polyheron (H)
     V-represented by M using LPs. Since it does not use the representation
     conversion, it should work for a large scale problem.
  */
  dd_rowrange i, m;
  dd_colrange d;
  dd_rowset redset;
  dd_matrixdata<T> *Mcopy;
  dd_setfamily *F = nullptr;

  m = M->rowsize;
  d = M->colsize;
  if (m <= 0 || d <= 0) {
    *error = dd_EmptyRepresentation;
    goto _L999;
  }
  Mcopy = dd_MatrixCopy(M);
  F = dd_CreateSetFamily(m, m);
  for (i = 1; i <= m; i++) {
    if (!set_member(i, M->linset)) {
      set_addelem(Mcopy->linset, i);
      redset = dd_SRedundantRows(Mcopy, error, maxiter, os);
      /* redset should contain all weakly nonadjacent ones */
      set_uni(redset, redset,
              Mcopy->linset); /* all linearity elements should be nonadjacent */
      set_compl(F->set[i - 1], redset); /* set the adjacency list of vertex i */
      set_delelem(Mcopy->linset, i);
      set_free(redset);
      if (*error != dd_NoError)
        goto _L99;
    }
  }
_L99:
  dd_FreeMatrix(Mcopy);
_L999:
  return F;
}

template <typename T>
bool dd_ImplicitLinearity(dd_matrixdata<T> *M, dd_rowrange itest,
                          T *certificate, dd_ErrorType *error,
                          size_t const& maxiter, std::ostream& os)
/* 092 */
{
  /* Checks whether the row itest is implicit linearity for the representation.
     All linearity rows are not checked and considered non implicit linearity
    (false). This code works for both H- and V-representations.  A certificate
    is given in the case of false, showing a feasible solution x satisfying the
    itest strict inequality for H-representation, a hyperplane RHS and normal
    (x_0, x) that separates the itest from the rest.  More explicitly, the LP to
    be setup is the same thing as redundancy case but with maximization:

     H-representation
       f* = maximize
         b_itest     + A_itest x
       subject to
         b_itest + 1 + A_itest x     >= 0 (relaxed inequality. This is not
    necessary but kept for simplicity of the code) b_{I-itest} + A_{I-itest} x
    >= 0 (all inequalities except for itest) b_L         + A_L x = 0.
    (linearity)

     V-representation (=separation problem)
       f* = maximize
         b_itest x_0     + A_itest x
       subject to
         b_itest x_0     + A_itest x     >= -1 (again, this is not necessary but
    kept for simplicity.) b_{I-itest} x_0 + A_{I-itest} x >=  0 (all
    nonlinearity generators except for itest in one side) b_L x_0         + A_L
    x = 0.  (linearity generators)

    Here, the input matrix is considered as (b, A), i.e. b corresponds to the
    first column of input and the row indices of input is partitioned into I and
    L where L is the set of linearity. In both cases, the itest data is implicit
    linearity if and only if the optimal value f* is nonpositive. The
    certificate has dimension one more for V-representation case.
  */

  dd_colrange j;
  dd_lpdata<T> *lp;
  dd_ErrorType err = dd_NoError;
  bool answer = false, localdebug = false;

  *error = dd_NoError;
  if (set_member(itest, M->linset)) {
    if (localdebug)
      printf(
          "The %ld th row is linearity and redundancy checking is skipped.\n",
          itest);
    return answer;
  }

  /* Create an LP data for redundancy checking */
  lp = dd_CreateLPData_from_M<T>(M);
  if (M->representation == dd_Generator) {
    dd_CreateLP_V_Redundancy(M, itest, lp);
  } else {
    dd_CreateLP_H_Redundancy(M, itest, lp);
  }

  lp->objective = dd_LPmax; /* the lp->objective is set by CreateLP* to LPmin */
  dd_LPSolve(lp, dd_choiceRedcheckAlgorithm, &err, maxiter, os);
  if (err != dd_NoError) {
    *error = err;
  } else {
    for (j = 0; j < lp->d; j++) {
      certificate[j] = lp->sol[j];
    }

    if (lp->LPS == dd_Optimal && lp->optvalue == 0) {
      answer = true;
      if (localdebug)
        fprintf(stdout, "==> %ld th data is an implicit linearity.\n", itest);
    } else {
      answer = false;
      if (localdebug)
        fprintf(stdout, "==> %ld th data is not an implicit linearity.\n",
                itest);
    }
  }
  dd_FreeLPData(lp);
  return answer;
}

template <typename T>
void dd_FreeOfImplicitLinearity(dd_matrixdata<T> *M, T *certificate,
                                dd_rowset *imp_linrows, dd_ErrorType *error,
                                size_t const& maxiter, std::ostream& os)
/* 092 */
{
  /* Checks whether the matrix M constains any implicit linearity at all.
  It returns 1 if it is free of any implicit linearity.  This means that
  the present linearity rows define the linearity correctly.  It returns
  nonpositive values otherwise.


     H-representation
       f* = maximize    z
       subject to
         b_I  + A_I x - 1 z >= 0
         b_L  + A_L x = 0  (linearity)
         z <= 1.

     V-representation (=separation problem)
       f* = maximize    z
       subject to
         b_I x_0 + A_I x - 1 z >= 0 (all nonlinearity generators in one side)
         b_L x_0 + A_L x  = 0  (linearity generators)
         z <= 1.

    Here, the input matrix is considered as (b, A), i.e. b corresponds to the
  first column of input and the row indices of input is partitioned into I and L
  where L is the set of linearity. In both cases, any implicit linearity exists
  if and only if the optimal value f* is nonpositive. The certificate has
  dimension one more for V-representation case.
  */

  dd_lpdata<T> *lp;
  dd_rowrange i, m;
  dd_colrange j, d1;
  dd_ErrorType err = dd_NoError;
  T *cvec; /* certificate for implicit linearity */

  bool localdebug = false;
  int answer = 0; // Used to return the value in the old CDD code.

  *error = dd_NoError;
  /* Create an LP data for redundancy checking */
  if (M->representation == dd_Generator) {
    lp = dd_CreateLP_V_ImplicitLinearity(M);
  } else {
    lp = dd_CreateLP_H_ImplicitLinearity(M);
  }

  dd_LPSolve(lp, dd_choiceRedcheckAlgorithm, &err, maxiter, os);
  if (err != dd_NoError) {
    *error = err;
    goto _L999;
  } else {

    for (j = 0; j < lp->d; j++) {
      certificate[j] = lp->sol[j];
    }

    if (localdebug)
      dd_WriteLPResult(std::cout, lp, err);

    /* *posset contains a set of row indices that are recognized as
     * nonlinearity.  */

    if (M->representation == dd_Generator) {
      d1 = M->colsize + 1;
    } else {
      d1 = M->colsize;
    }
    m = M->rowsize;
    dd_AllocateArow(d1, &cvec);
    set_initialize(imp_linrows, m);

    if (lp->LPS == dd_Optimal) {
      if (lp->optvalue > 0) {
        answer = 1;
        if (localdebug)
          fprintf(stdout, "==> The matrix has no implicit linearity.\n");
      } else {
        if (lp->optvalue < 0) {
          answer = -1;
          if (localdebug)
            fprintf(stdout, "==> The matrix defines the trivial system.\n");
        } else {
          answer = 0;
          if (localdebug)
            fprintf(stdout, "==> The matrix has some implicit linearity.\n");
        }
      }
    } else {
      answer = -2;
      if (localdebug)
        fprintf(stdout, "==> The LP fails.\n");
    }
    if (answer == 0) {
      /* List the implicit linearity rows */
      for (i = m; i >= 1; i--) {
        if (!set_member(i, lp->posset_extra)) {
          if (dd_ImplicitLinearity(M, i, cvec, error, maxiter, os)) {
            set_addelem(*imp_linrows, i);
            if (localdebug) {
              fprintf(stdout, " row %ld is implicit linearity\n", i);
              fprintf(stdout, "\n");
            }
          }
          if (*error != dd_NoError)
            goto _L999;
        }
      }
    } /* end of if (answer==0) */
    if (answer == -1) {
      for (i = m; i >= 1; i--)
        set_addelem(*imp_linrows, i);
    } /* all rows are considered implicit linearity */

    dd_FreeArow(cvec);
  }
_L999:
  dd_FreeLPData(lp);
}

template <typename T>
dd_rowset dd_ImplicitLinearityRows(dd_matrixdata<T> *M,
                                   dd_ErrorType *error,
                                   size_t const& maxiter, std::ostream& os)
{
  dd_colrange d;
  dd_rowset imp_linset;
  T *cvec; /* certificate */

  if (M->representation == dd_Generator) {
    d = M->colsize + 2;
  } else {
    d = M->colsize + 1;
  }

  dd_AllocateArow(d, &cvec);
  dd_FreeOfImplicitLinearity(M, cvec, &imp_linset, error, maxiter, os);

  dd_FreeArow(cvec);
  return imp_linset;
}

template <typename T>
bool dd_MatrixCanonicalizeLinearity(dd_matrixdata<T> **M,
                                    dd_rowset *impl_linset, dd_rowindex *newpos,
                                    dd_ErrorType *error,
                                    size_t const& maxiter, std::ostream& os)
{
  /* This is to recongnize all implicit linearities, and put all linearities at
     the top of the matrix.    All implicit linearities will be returned by
     *impl_linset.
  */
  dd_rowset linrows, ignoredrows, basisrows;
  dd_colset ignoredcols, basiscols;
  dd_rowrange i, k, m;
  dd_rowindex newpos1;

  linrows = dd_ImplicitLinearityRows(*M, error, maxiter, os);
  if (*error != dd_NoError)
    return false;

  m = (*M)->rowsize;

  set_uni((*M)->linset, (*M)->linset, linrows);
  /* add the implicit linrows to the explicit linearity rows */

  /* To remove redundancy of the linearity part,
     we need to compute the rank and a basis of the linearity part. */
  set_initialize(&ignoredrows, (*M)->rowsize);
  set_initialize(&ignoredcols, (*M)->colsize);
  set_compl(ignoredrows, (*M)->linset);
  (void)dd_MatrixRank(*M, ignoredrows, ignoredcols, &basisrows, &basiscols);
  set_diff(ignoredrows, (*M)->linset, basisrows);
  dd_MatrixRowsRemove2(M, ignoredrows, newpos);

  dd_MatrixShiftupLinearity(M, &newpos1);

  for (i = 1; i <= m; i++) {
    k = (*newpos)[i];
    if (k > 0) {
      (*newpos)[i] = newpos1[k];
    }
  }

  *impl_linset = linrows;
  delete[] newpos1;
  set_free(basisrows);
  set_free(basiscols);
  set_free(ignoredrows);
  set_free(ignoredcols);
  return true;
}

template <typename T>
bool dd_MatrixCanonicalize(dd_matrixdata<T> **M, dd_rowset *impl_linset,
                           dd_rowset *redset, dd_rowindex *newpos,
                           dd_ErrorType *error, size_t const& maxiter, std::ostream& os) {
  /* This is to find a canonical representation of a matrix *M by
     recognizing all implicit linearities and all redundancies.
     All implicit linearities will be returned by *impl_linset and
     redundancies will be returned by *redset.
  */
  dd_rowrange i, k, m;
  dd_rowindex newpos1, revpos;
  dd_rowset redset1;
  bool success = true;

  m = (*M)->rowsize;
  set_initialize(redset, m);
  revpos = new long[m + 1];
  for (i = 0; i <= m; i++)
    revpos[i] = 0;

  success = dd_MatrixCanonicalizeLinearity(M, impl_linset, newpos, error);

  if (!success)
    goto _L99;

  for (i = 1; i <= m; i++) {
    k = (*newpos)[i];
    if (k > 0)
      revpos[k] = i; /* inverse of *newpos[] */
  }

  success = dd_MatrixRedundancyRemove(M, &redset1, &newpos1, error, maxiter, os);

  if (!success)
    goto _L99;

  for (i = 1; i <= m; i++) {
    k = (*newpos)[i];
    if (k > 0) {
      (*newpos)[i] = newpos1[k];
      if (newpos1[k] < 0)
        (*newpos)[i] = -revpos[-newpos1[k]]; /* update the certificate of its
                                                duplicate removal. */
      if (set_member(k, redset1))
        set_addelem(*redset, i);
    }
  }

_L99:
  set_free(redset1);
  delete[] newpos1;
  delete[] revpos;
  return success;
}

template <typename T>
bool dd_ExistsRestrictedFace(dd_matrixdata<T> *M, dd_rowset R, dd_rowset S,
                             dd_ErrorType *err, size_t const& maxiter, std::ostream& os)
/* 0.94 */
{
  /* This function checkes if there is a point that satifies all the constraints
  of the matrix M (interpreted as an H-representation) with additional equality
  contraints specified by R and additional strict inequality constraints
  specified by S. The set S is supposed to be disjoint from both R and
  M->linset.   When it is not, the set S will be considered as S\(R U
  M->linset).
  */
  bool answer = false;
  dd_lpdata<T> *lp = nullptr;

  lp = dd_Matrix2Feasibility2(M, R, S, err);

  if (*err != dd_NoError)
    return answer;

  /* Solve the LP by cdd LP solver. */
  dd_LPSolve(lp, dd_DualSimplex, err, maxiter, os);
  if (*err != dd_NoError)
    return answer;
  if (lp->LPS == dd_Optimal && lp->optvalue > 0) {
    answer = true;
  }

  dd_FreeLPData(lp);
  return answer;
}

template <typename T>
dd_rowrange dd_RayShooting(dd_matrixdata<T> *M, T *p, T *r) {
  dd_rowrange imin = -1, i, m;
  dd_colrange j, d;
  T *vecmin, *vec;
  T min, t1, t2, alpha, t1min;
  bool started = false;
  bool localdebug = false;
  T dd_one;
  dd_one = 1;
  m = M->rowsize;
  d = M->colsize;
  if (dd_one != p[0]) {
    fprintf(stdout, "Warning: RayShooting is called with a point with first "
                    "coordinate not 1.\n");
    p[0] = dd_one;
  }
  if (r[0] != 0) {
    fprintf(stdout, "Warning: RayShooting is called with a direction with "
                    "first coordinate not 0.\n");
    r[0] = 0;
  }

  dd_AllocateArow(d, &vecmin);
  dd_AllocateArow(d, &vec);

  for (i = 1; i <= m; i++) {
    dd_InnerProduct(t1, d, M->matrix[i - 1], p);
    if (localdebug)
      std::cout << "dd_RayShooting: i=" << i << " t1=" << t1 << "\n";
    if (t1 > 0) {
      dd_InnerProduct(t2, d, M->matrix[i - 1], r);
      bool is_field = true;
      if (is_field) {
        alpha = t2 / t1;
        if (!started) {
          imin = i;
          min = alpha;
          t1min = t1; /* store the denominator. */
          started = true;
          if (localdebug)
            std::cout << "dd_RayShooting: Level 1: imin = " << imin
                      << " and min = " << min << "\n";
        } else {
          if (alpha < min) {
            imin = i;
            min = alpha;
            t1min = t1; /* store the denominator. */
            if (localdebug)
              std::cout << "dd_RayShootni: Level 2: imin = " << imin
                        << " and min = " << min << "\n";
          } else {
            if (alpha == min) { /* tie break */
              for (j = 1; j <= d; j++) {
                vecmin[j - 1] = M->matrix[imin - 1][j - 1] / t1min;
                vec[j - 1] = M->matrix[i - 1][j - 1] / t1;
              }
              if (dd_LexSmaller(vec, vecmin, d)) {
                imin = i;
                min = alpha;
                t1min = t1; /* store the denominator. */
                if (localdebug)
                  std::cout << "dd_RayShooting: Level 3: imin = " << imin
                            << " and min = " << min << "\n";
              }
            }
          }
        }
      } else {
        if (!started) {
          imin = i;
          min = t2;
          t1min = t1; /* store the denominator. */
          started = true;
          if (localdebug)
            std::cout << "dd_RayShooting: Level 1: imin = " << imin
                      << " and min = " << min << "\n";
        } else {
          if (dd_SmallerFrac(t2, t1, min, t1min)) {
            imin = i;
            min = t2;
            t1min = t1; /* store the denominator. */
            if (localdebug)
              std::cout << "dd_RayShootni: Level 2: imin = " << imin
                        << " and min = " << min << "\n";
          } else {
            if (dd_EqualFrac(t2, t1, min, t1min)) { /* tie break */
              if (dd_LexSmallerFrac(M->matrix[i - 1], t1, M->matrix[imin - 1],
                                    t1min, d)) {
                imin = i;
                min = t2;
                t1min = t1; /* store the denominator. */
                if (localdebug)
                  std::cout << "dd_RayShooting: Level 3: imin = " << imin
                            << " and min = " << min << "\n";
              }
            }
          }
        }
      }
    }
  }

  dd_FreeArow(vecmin);
  dd_FreeArow(vec);
  return imin;
}

template <typename T>
void dd_BasisStatusMaximize(dd_rowrange m_size, dd_colrange d_size, T **A,
                            T **Ts, dd_rowset equalityset, dd_rowrange objrow,
                            dd_colrange rhscol, dd_LPStatusType LPS,
                            T &optvalue, T *sol, T *dsol, dd_rowset posset,
                            dd_colindex nbindex, dd_rowrange re, dd_colrange se,
                            dd_colrange *nse, long *pivots, bool *found,
                            bool *LPScorrect, data_temp_simplex<T> *data,
                            std::ostream & os)
/*  This is just to check whether the status LPS of the basis given by
nbindex with extra certificates se or re is correct.  It is done
by recomputing the basis inverse matrix T.  It does not solve the LP
when the status *LPS is undecided.  Thus the input is
m_size, d_size, A, equalityset, LPS, nbindex, re and se.
Other values will be recomputed from scratch.

The main purpose of the function is to verify the correctness
of the result of floating point computation with the GMP rational
arithmetics.
*/
{
  long pivots0, fbasisrank;
  dd_rowrange i, is;
  dd_colrange s, senew, j;
  unsigned int rseed = 1;
  T val;

  std::vector<long> nbtemp(d_size + 1, 0);
  for (i = 0; i <= 4; i++)
    pivots[i] = 0;

  /* Initializing control variables. */
  dd_ComputeRowOrderVector2(m_size, d_size, A, data->OrderVector, dd_MinIndex,
                            rseed);

  dd_ResetTableau(m_size, d_size, Ts, nbtemp.data(), data->bflag, objrow,
                  rhscol);

  is = nbindex[se];

  fbasisrank = d_size - 1;
  for (j = 1; j <= d_size; j++) {
    if (nbindex[j] < 0)
      fbasisrank = fbasisrank - 1;
    /* fbasisrank=the basis rank computed by floating-point */
  }

  if (fbasisrank < d_size - 1) {
    *found = false;
    return;
    /* Suspicious case.  Rerun the LP solver with GMP. */
  }

  dd_FindLPBasis2(m_size, d_size, A, Ts, data->OrderVector, equalityset,
                  nbindex, data->bflag, objrow, rhscol, &s, found, &pivots0,
                  data);

  /* set up the new se column and corresponding variable */
  senew = data->bflag[is];
  is = nbindex[senew];

  pivots[4] = pivots0; /*GMP postopt pivots */

  if (!(*found)) {
    return;
  }

  /* Check whether a recomputed basis is of the type specified by LPS */
  *LPScorrect = true;
  switch (LPS) {
  case dd_Optimal:
    for (i = 1; i <= m_size; i++) {
      if (i != objrow && data->bflag[i] == -1) { /* i is a basic variable */
        dd_TableauEntry(val, d_size, A, Ts, i, rhscol);
        if (val < 0) {
          *LPScorrect = false;
          break;
        }
      } else {
        if (data->bflag[i] > 0) { /* i is nonbasic variable */
          dd_TableauEntry(val, d_size, A, Ts, objrow, data->bflag[i]);
          if (val > 0) {
            *LPScorrect = false;
            break;
          }
        }
      }
    }
    break;
  case dd_Inconsistent:
    for (j = 1; j <= d_size; j++) {
      dd_TableauEntry(val, d_size, A, Ts, re, j);
      if (j == rhscol) {
        if (val >= 0) {
          *LPScorrect = false;
          break;
        }
      } else {
        if (val > 0) {
          *LPScorrect = false;
          break;
        }
      }
    }
    break;
  case dd_LPSundecided:
    break;
  case dd_StrucInconsistent:
    break;
  case dd_StrucDualInconsistent:
    break;
  case dd_Unbounded:
    break;
  case dd_DualUnbounded:
    break;
  case dd_TooManyIterations:
#ifdef DEBUG_CDD
    os << "CDD: It should not pass by this entry\n";
    throw TerminalException{1};
#endif
    break;

  case dd_DualInconsistent:
    for (i = 1; i <= m_size; i++) {
      dd_TableauEntry(val, d_size, A, Ts, i, data->bflag[is]);
      if (i == objrow) {
        if (val <= 0) {
          *LPScorrect = false;
          break;
        }
      } else {
        if (val < 0) {
          *LPScorrect = false;
          break;
        }
      }
    };
    break;
  }
#ifdef DEBUG_CDD
  os << "CDD: Before dd_SetSolutions in dd_BasisStatusMaximize\n";
#endif
  dd_SetSolutions(m_size, d_size, A, Ts, objrow, rhscol, LPS, optvalue, sol,
                  dsol, posset, re, senew, data->bflag, os);
  *nse = senew;
}

template <typename T>
void dd_BasisStatusMinimize(dd_rowrange m_size, dd_colrange d_size, T **A,
                            T **Ts, dd_rowset equalityset, dd_rowrange objrow,
                            dd_colrange rhscol, dd_LPStatusType LPS,
                            T &optvalue, T *sol, T *dsol, dd_rowset posset,
                            dd_colindex nbindex, dd_rowrange re, dd_colrange se,
                            dd_colrange *nse, long *pivots, bool *found,
                            bool *LPScorrect, data_temp_simplex<T> *data,
                            std::ostream & os) {
  dd_colrange j;

  for (j = 1; j <= d_size; j++)
    A[objrow - 1][j - 1] = -A[objrow - 1][j - 1];
  dd_BasisStatusMaximize(m_size, d_size, A, Ts, equalityset, objrow, rhscol,
                         LPS, optvalue, sol, dsol, posset, nbindex, re, se, nse,
                         pivots, found, LPScorrect, data, os);
  optvalue = -optvalue;
  for (j = 1; j <= d_size; j++) {
    if (LPS != dd_Inconsistent) {
      /* Inconsistent certificate stays valid for minimization, 0.94e */
      dsol[j - 1] = -dsol[j - 1];
    }
    A[objrow - 1][j - 1] = -A[objrow - 1][j - 1];
  }
}

template <typename T>
void dd_BasisStatus(dd_lpdata<double> *lpf, dd_lpdata<T> *lp, bool *LPScorrect,
                    data_temp_simplex<T> *data, std::ostream& os) {
  int i;
  dd_colrange se, j;
  bool basisfound;

  switch (lp->objective) {
  case dd_LPmax:
    dd_BasisStatusMaximize<T>(lp->m, lp->d, lp->A, lp->B, lp->equalityset,
                              lp->objrow, lp->rhscol, lpf->LPS, lp->optvalue,
                              lp->sol, lp->dsol, lp->posset_extra, lpf->nbindex,
                              lpf->re, lpf->se, &se, lp->pivots, &basisfound,
                              LPScorrect, data, os);
    if (*LPScorrect) {
      /* printf("BasisStatus Check: the current basis is verified with GMP\n");
       */
      lp->LPS = lpf->LPS;
      lp->re = lpf->re;
      lp->se = se;
      for (j = 1; j <= lp->d; j++)
        lp->nbindex[j] = lpf->nbindex[j];
    }
    for (i = 1; i <= 5; i++)
      lp->pivots[i - 1] += lpf->pivots[i - 1];
    break;
  case dd_LPmin:
    dd_BasisStatusMinimize<T>(lp->m, lp->d, lp->A, lp->B, lp->equalityset,
                              lp->objrow, lp->rhscol, lpf->LPS, lp->optvalue,
                              lp->sol, lp->dsol, lp->posset_extra, lpf->nbindex,
                              lpf->re, lpf->se, &se, lp->pivots, &basisfound,
                              LPScorrect, data, os);
    if (*LPScorrect) {
      /* printf("BasisStatus Check: the current basis is verified with GMP\n");
       */
      lp->LPS = lpf->LPS;
      lp->re = lpf->re;
      lp->se = se;
      for (j = 1; j <= lp->d; j++) {
        lp->nbindex[j] = lpf->nbindex[j];
      }
    }
    for (i = 1; i <= 5; i++)
      lp->pivots[i - 1] += lpf->pivots[i - 1];
    break;
  case dd_LPnone:
    break;
  }
}

/* end of cddlp.c */

/* cddcore.c:  Core Procedures for cddlib
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.94g, March 23, 2012
*/

/* cddlib : C-library of the double description method for
   computing all vertices and extreme rays of the polyhedron
   P= {x :  b - A x >= 0}.
   Please read COPYING (GNU General Public Licence) and
   the manual cddlibman.tex for detail.
*/

template <typename T>
void dd_CheckAdjacency(dd_conedata<T> *cone, dd_raydata<T> **RP1,
                       dd_raydata<T> **RP2, bool *adjacent) {
  dd_raydata<T> *TempRay;
  bool localdebug = false;
  dd_rowset Face, Face1;
  set_initialize(&Face, cone->m);
  set_initialize(&Face1, cone->m);

  *adjacent = true;
  set_int(Face1, (*RP1)->ZeroSet, (*RP2)->ZeroSet);
  set_int(Face, Face1, cone->AddedHalfspaces);
  if (set_card(Face) < cone->d - 2) {
    *adjacent = false;
    if (localdebug) {
      fprintf(stdout, "non adjacent: set_card(face) %ld < %ld = cone->d.\n",
              set_card(Face), cone->d);
    }
    set_free(Face);
    set_free(Face1);
    return;
  } else {
    if (cone->parent->NondegAssumed) {
      *adjacent = true;
      set_free(Face);
      set_free(Face1);
      return;
    }
  }
  TempRay = cone->FirstRay;
  while (TempRay != nullptr && *adjacent) {
    if (TempRay != *RP1 && TempRay != *RP2) {
      set_int(Face1, TempRay->ZeroSet, cone->AddedHalfspaces);
      if (set_subset(Face, Face1))
        *adjacent = false;
    }
    TempRay = TempRay->Next;
  }
  set_free(Face);
  set_free(Face1);
}

template <typename T>
void dd_Eliminate(dd_conedata<T> *cone, dd_raydata<T> **Ptr) {
  /*eliminate the record pointed by Ptr->Next*/
  dd_raydata<T> *TempPtr;

  TempPtr = (*Ptr)->Next;
  (*Ptr)->Next = (*Ptr)->Next->Next;
  if (TempPtr == cone->FirstRay) /*Update the first pointer*/
    cone->FirstRay = (*Ptr)->Next;
  if (TempPtr == cone->LastRay) /*Update the last pointer*/
    cone->LastRay = *Ptr;

  delete[] TempPtr->Ray;      /* free the ray vector memory */
  set_free(TempPtr->ZeroSet); /* free the ZeroSet memory */
  delete TempPtr;             /* free the dd_Ray structure memory */
  cone->RayCount--;
}

template <typename T> void dd_SetInequalitySets(dd_conedata<T> *cone) {
  dd_rowrange i;

  set_emptyset(cone->GroundSet);
  set_emptyset(cone->EqualitySet);
  set_emptyset(cone->NonequalitySet);
  for (i = 1; i <= (cone->parent->m); i++) {
    set_addelem(cone->GroundSet, i);
    if (cone->parent->EqualityIndex[i] == 1)
      set_addelem(cone->EqualitySet, i);
    if (cone->parent->EqualityIndex[i] == -1)
      set_addelem(cone->NonequalitySet, i);
  }
}

template <typename T>
void dd_AValue(T *val, dd_colrange d_size, T **A, T *p, dd_rowrange i) {
  /*return the ith component of the vector  A x p */
  dd_colrange j;

  *val = 0;
  /* Changed by Marc Pfetsch 010219 */

  for (j = 0; j < d_size; j++)
    *val += A[i - 1][j] * p[j];
}

template <typename T>
void dd_StoreRay1(dd_conedata<T> *cone, T *p,
                  bool *feasible) { /* Original ray storing routine when
                                       RelaxedEnumeration is false */
  dd_rowrange i, k, fii = cone->m + 1;
  dd_colrange j;
  T temp;
  dd_raydata<T> *RR;
  bool localdebug = false;

  RR = cone->LastRay;
  *feasible = true;
  set_initialize(&(RR->ZeroSet), cone->m);
  for (j = 0; j < cone->d; j++) {
    RR->Ray[j] = p[j];
  }
  for (i = 1; i <= cone->m; i++) {
    k = cone->OrderVector[i];
    dd_AValue(&temp, cone->d, cone->A, p, k);
    if (temp == 0) {
      set_addelem(RR->ZeroSet, k);
      if (localdebug) {
        fprintf(stdout, "recognized zero!\n");
      }
    }
    if (temp < 0) {
      if (localdebug) {
        fprintf(stdout, "recognized negative!\n");
      }
      *feasible = false;
      if (fii > cone->m)
        fii = i; /* the first violating inequality index */
    }
  }
  RR->FirstInfeasIndex = fii;
  RR->feasible = *feasible;
}

template <typename T>
void dd_StoreRay2(dd_conedata<T> *cone, T *p, bool *feasible,
                  bool *weaklyfeasible)
/* Ray storing routine when RelaxedEnumeration is true.
    weaklyfeasible is true iff it is feasible with
    the strict_inequality conditions deleted. */
{
  dd_raydata<T> *RR;
  dd_rowrange i, k, fii = cone->m + 1;
  dd_colrange j;
  T temp;

  RR = cone->LastRay;
  *feasible = true;
  *weaklyfeasible = true;
  set_initialize(&(RR->ZeroSet), cone->m);
  for (j = 0; j < cone->d; j++) {
    RR->Ray[j] = p[j];
  }
  for (i = 1; i <= cone->m; i++) {
    k = cone->OrderVector[i];
    dd_AValue(&temp, cone->d, cone->A, p, k);
    if (temp == 0) {
      set_addelem(RR->ZeroSet, k);
      if (cone->parent->EqualityIndex[k] == -1)
        *feasible = false; /* strict inequality required */
    }
    if (temp < 0) {
      *feasible = false;
      if (fii > cone->m && cone->parent->EqualityIndex[k] >= 0) {
        fii = i; /* the first violating inequality index */
        *weaklyfeasible = false;
      }
    }
  }
  RR->FirstInfeasIndex = fii;
  RR->feasible = *feasible;
}

template <typename T> void dd_AddRay(dd_conedata<T> *cone, T *p) {
  bool feasible, weaklyfeasible;
  bool localdebug = false;

  if (cone->FirstRay == nullptr) {
    cone->FirstRay = new dd_raydata<T>;
    cone->FirstRay->Ray = new T[cone->d];
    if (localdebug)
      fprintf(stdout, "Create the first ray pointer\n");
    cone->LastRay = cone->FirstRay;
    cone->ArtificialRay->Next = cone->FirstRay;
  } else {
    cone->LastRay->Next = new dd_raydata<T>;
    cone->LastRay->Next->Ray = new T[cone->d];
    if (localdebug)
      fprintf(stdout, "Create a new ray pointer\n");
    cone->LastRay = cone->LastRay->Next;
  }
  cone->LastRay->Next = nullptr;
  cone->RayCount++;
  cone->TotalRayCount++;
  if (localdebug) {
    if (cone->TotalRayCount % 100 == 0) {
      fprintf(stdout,
              "*Rays (Total, Currently Active, Feasible) =%8ld%8ld%8ld\n",
              cone->TotalRayCount, cone->RayCount, cone->FeasibleRayCount);
    }
  }
  if (cone->parent->RelaxedEnumeration) {
    dd_StoreRay2(cone, p, &feasible, &weaklyfeasible);
    if (weaklyfeasible)
      cone->WeaklyFeasibleRayCount++;
  } else {
    dd_StoreRay1(cone, p, &feasible);
    if (feasible)
      cone->WeaklyFeasibleRayCount++;
    /* weaklyfeasible is equiv. to feasible in this case. */
  }
  if (feasible)
    cone->FeasibleRayCount++;
}

template <typename T> void dd_AddArtificialRay(dd_conedata<T> *cone) {
  bool localdebug = false;

  if (cone->ArtificialRay != nullptr) {
    fprintf(stdout, "Warning !!!  FirstRay in not nil.  Illegal Call\n");
    return;
  }

  T *zerovector;
  dd_colrange d1;
  if (cone->d <= 0)
    d1 = 1;
  else
    d1 = cone->d;
  dd_AllocateArow(d1, &zerovector);
  cone->ArtificialRay = new dd_raydata<T>;
  cone->ArtificialRay->Ray = new T[d1];

  if (localdebug)
    fprintf(stdout, "Create the artificial ray pointer\n");

  cone->LastRay = cone->ArtificialRay;
  bool feasible;
  dd_StoreRay1(cone, zerovector, &feasible);
  cone->ArtificialRay->Next = nullptr;
  delete[] zerovector; /* 086 */
}

template <typename T>
void dd_ConditionalAddEdge(dd_conedata<T> *cone, dd_raydata<T> *Ray1,
                           dd_raydata<T> *Ray2, dd_raydata<T> *ValidFirstRay) {
  long it, it_row, fii1, fii2, fmin, fmax;
  bool adjacent, lastchance;
  dd_raydata<T> *TempRay;
  dd_raydata<T> *Rmin;
  dd_raydata<T> *Rmax;
  dd_adjacencydata<T> *NewEdge;
  bool localdebug = false;
  dd_rowset ZSmin, ZSmax;
  dd_rowset Face, Face1;

  set_initialize(&Face, cone->m);
  set_initialize(&Face1, cone->m);

  fii1 = Ray1->FirstInfeasIndex;
  fii2 = Ray2->FirstInfeasIndex;
  if (fii1 < fii2) {
    fmin = fii1;
    fmax = fii2;
    Rmin = Ray1;
    Rmax = Ray2;
  } else {
    fmin = fii2;
    fmax = fii1;
    Rmin = Ray2;
    Rmax = Ray1;
  }
  ZSmin = Rmin->ZeroSet;
  ZSmax = Rmax->ZeroSet;
  if (localdebug) {
    fprintf(stdout, "dd_ConditionalAddEdge: FMIN = %ld (row%ld)   FMAX=%ld\n",
            fmin, cone->OrderVector[fmin], fmax);
  }
  if (fmin == fmax) {
    if (localdebug)
      fprintf(stdout,
              "dd_ConditionalAddEdge: equal FII value-> No edge added\n");
  } else if (set_member(cone->OrderVector[fmin], ZSmax)) {
    if (localdebug)
      fprintf(stdout,
              "dd_ConditionalAddEdge: No strong separation -> No edge added\n");
  } else { /* the pair will be separated at the iteration fmin */
    lastchance = true;
    /* flag to check it will be the last chance to store the edge candidate */
    set_int(Face1, ZSmax, ZSmin);
    cone->count_int++;
    if (localdebug) {
      fprintf(stdout, "Face: ");
      for (it = 1; it <= cone->m; it++) {
        it_row = cone->OrderVector[it];
        if (set_member(it_row, Face1))
          fprintf(stdout, "%ld ", it_row);
      }
      fprintf(stdout, "\n");
    }
    for (it = cone->Iteration + 1; it < fmin && lastchance; it++) {
      it_row = cone->OrderVector[it];
      if (cone->parent->EqualityIndex[it_row] >= 0 &&
          set_member(it_row, Face1)) {
        lastchance = false;
        (cone->count_int_bad)++;
      }
    }
    if (lastchance) {
      adjacent = true;
      (cone->count_int_good)++;
      /* adjacent checking */
      set_int(Face, Face1, cone->AddedHalfspaces);
      if (set_card(Face) < cone->d - 2) {
        adjacent = false;
      } else if (cone->parent->NondegAssumed) {
        adjacent = true;
      } else {
        TempRay = ValidFirstRay; /* the first ray for adjacency checking */
        while (TempRay != nullptr && adjacent) {
          if (TempRay != Ray1 && TempRay != Ray2) {
            set_int(Face1, TempRay->ZeroSet, cone->AddedHalfspaces);
            if (set_subset(Face, Face1)) {
              adjacent = false;
            }
          }
          TempRay = TempRay->Next;
        }
      }
      if (adjacent) {
        NewEdge = new dd_adjacencydata<T>;
        NewEdge->Ray1 =
            Rmax; /* save the one remains in iteration fmin in the first */
        NewEdge->Ray2 =
            Rmin; /* save the one deleted in iteration fmin in the second */
        NewEdge->Next = nullptr;
        (cone->EdgeCount)++;
        (cone->TotalEdgeCount)++;
        if (cone->Edges[fmin] == nullptr) {
          cone->Edges[fmin] = NewEdge;
        } else {
          NewEdge->Next = cone->Edges[fmin];
          cone->Edges[fmin] = NewEdge;
        }
      }
    }
  }
  set_free(Face);
  set_free(Face1);
}

template <typename T> void dd_CreateInitialEdges(dd_conedata<T> *cone) {
  dd_raydata<T> *Ptr1;
  dd_raydata<T> *Ptr2;
  dd_rowrange fii1, fii2;
  bool adj;

  cone->Iteration = cone->d; /* CHECK */
  if (cone->FirstRay == nullptr || cone->LastRay == nullptr) {
    /* fprintf(stdout,"Warning: dd_ CreateInitialEdges called with nullptr
     * pointer(s)\n"); */
    return;
  }
  Ptr1 = cone->FirstRay;
  while (Ptr1 != cone->LastRay && Ptr1 != nullptr) {
    fii1 = Ptr1->FirstInfeasIndex;
    Ptr2 = Ptr1->Next;
    while (Ptr2 != nullptr) {
      fii2 = Ptr2->FirstInfeasIndex;
      dd_CheckAdjacency(cone, &Ptr1, &Ptr2, &adj);
      if (fii1 != fii2 && adj)
        dd_ConditionalAddEdge(cone, Ptr1, Ptr2, cone->FirstRay);
      Ptr2 = Ptr2->Next;
    }
    Ptr1 = Ptr1->Next;
  }
}

template <typename T>
void dd_UpdateEdges(dd_conedata<T> *cone, dd_raydata<T> *RRbegin,
                    dd_raydata<T> *RRend)
/* This procedure must be called after the ray list is sorted
   by dd_EvaluateARay2 so that FirstInfeasIndex's are monotonically
   increasing.
*/
{
  dd_raydata<T> *Ptr1;
  dd_raydata<T> *Ptr2;
  dd_raydata<T> *Ptr2begin;
  dd_rowrange fii1;
  bool ptr2found, quit;
  long pos1, pos2;
  float workleft, prevworkleft = 110.0, totalpairs;
  bool localdebug = false;

  totalpairs = (cone->ZeroRayCount - 1.0) * (cone->ZeroRayCount - 2.0) + 1.0;
  Ptr2begin = nullptr;
  if (RRbegin == nullptr || RRend == nullptr) {
    if (1)
      fprintf(stdout,
              "Warning: dd_UpdateEdges called with nullptr pointer(s)\n");
    return;
  }
  Ptr1 = RRbegin;
  pos1 = 1;
  do {
    ptr2found = false;
    quit = false;
    fii1 = Ptr1->FirstInfeasIndex;
    pos2 = 2;
    for (Ptr2 = Ptr1->Next; !ptr2found && !quit; Ptr2 = Ptr2->Next, pos2++) {
      if (Ptr2->FirstInfeasIndex > fii1) {
        Ptr2begin = Ptr2;
        ptr2found = true;
      } else {
        if (Ptr2 == RRend)
          quit = true;
      }
    }
    if (ptr2found) {
      quit = false;
      for (Ptr2 = Ptr2begin; !quit; Ptr2 = Ptr2->Next) {
        dd_ConditionalAddEdge(cone, Ptr1, Ptr2, RRbegin);
        if (Ptr2 == RRend || Ptr2->Next == nullptr)
          quit = true;
      }
    }
    Ptr1 = Ptr1->Next;
    pos1++;
    workleft = 100.0 * (cone->ZeroRayCount - pos1) *
               (cone->ZeroRayCount - pos1 - 1.0) / totalpairs;
    if (localdebug) {
      if (cone->ZeroRayCount >= 500 && pos1 % 10 == 0 &&
          prevworkleft - workleft >= 10) {
        fprintf(stdout,
                "*Work of iteration %5ld(/%ld): %4ld/%4ld => %4.1f%% left\n",
                cone->Iteration, cone->m, pos1, cone->ZeroRayCount, workleft);
        prevworkleft = workleft;
      }
    }
  } while (Ptr1 != RRend && Ptr1 != nullptr);
}

template <typename T> void dd_FreeDDMemory0(dd_conedata<T> *cone) {
  dd_raydata<T> *Ptr;
  dd_raydata<T> *PrevPtr;

  /* THIS SHOULD BE REWRITTEN carefully */
  PrevPtr = cone->ArtificialRay;
  if (PrevPtr != nullptr) {
    for (Ptr = cone->ArtificialRay->Next; Ptr != nullptr; Ptr = Ptr->Next) {
      delete[] PrevPtr->Ray;
      delete[] PrevPtr->ZeroSet;
      delete PrevPtr;
      PrevPtr = Ptr;
    }
    cone->FirstRay = nullptr;

    delete[] cone->LastRay->Ray;
    cone->LastRay->Ray = nullptr;
    set_free(cone->LastRay->ZeroSet);
    cone->LastRay->ZeroSet = nullptr;
    delete cone->LastRay;
    cone->LastRay = nullptr;
    cone->ArtificialRay = nullptr;
  }
  /* must add (by Sato) */
  delete[] cone->Edges;

  set_free(cone->GroundSet);
  set_free(cone->EqualitySet);
  set_free(cone->NonequalitySet);
  set_free(cone->AddedHalfspaces);
  set_free(cone->WeaklyAddedHalfspaces);
  set_free(cone->InitialHalfspaces);
  delete[] cone->InitialRayIndex;
  delete[] cone->OrderVector;
  delete[] cone->newcol;

  /* Fixed by Shawn Rusaw.  Originally it was cone->d instead of cone->d_alloc
   */
  dd_FreeBmatrix(cone->d_alloc, cone->B);
  dd_FreeBmatrix(cone->d_alloc, cone->Bsave);

  /* Fixed by Marc Pfetsch 010219*/
  dd_FreeAmatrix(cone->m_alloc, cone->A);
  cone->A = nullptr;

  delete cone;
}

template <typename T> void dd_FreeDDMemory(dd_polyhedradata<T> *poly) {
  dd_FreeDDMemory0(poly->child);
  poly->child = nullptr;
}

template <typename T> void dd_FreePolyhedra(dd_polyhedradata<T> *poly) {
  dd_bigrange i;

  if ((poly)->child != nullptr)
    dd_FreeDDMemory(poly);
  dd_FreeAmatrix(poly->m_alloc, poly->A);
  if (poly->c != nullptr)
    dd_FreeArow(poly->c);
  delete[] poly->EqualityIndex;
  if (poly->AincGenerated) {
    for (i = 1; i <= poly->m1; i++) {
      set_free(poly->Ainc[i - 1]);
    }
    delete[] poly->Ainc;
    set_free(poly->Ared);
    set_free(poly->Adom);
    poly->Ainc = nullptr;
  }

  delete poly;
}

template <typename T>
void dd_ZeroIndexSet(dd_rowrange m_size, dd_colrange d_size, T **A, T *x,
                     dd_rowset ZS) {
  dd_rowrange i;
  T temp;

  /* Changed by Marc Pfetsch 010219 */
  set_emptyset(ZS);
  for (i = 1; i <= m_size; i++) {
    dd_AValue(&temp, d_size, A, x, i);
    if (temp == 0)
      set_addelem(ZS, i);
  }
}

template <typename T>
void dd_CopyBmatrix(dd_colrange d_size, T **Ts, T **TCOPY) {
  for (dd_rowrange i = 0; i < d_size; i++)
    for (dd_colrange j = 0; j < d_size; j++)
      TCOPY[i][j] = Ts[i][j];
}

template <typename T>
void dd_CopyNormalizedArow(T *acopy, T *a, dd_colrange d) {
  dd_CopyArow(acopy, a, d);
}

template <typename T>
void dd_CopyNormalizedAmatrix(T **Acopy, T **A, dd_rowrange m, dd_colrange d) {
  dd_rowrange i;

  for (i = 0; i < m; i++) {
    dd_CopyNormalizedArow(Acopy[i], A[i], d);
  }
}

template <typename T>
void dd_PermuteCopyAmatrix(T **Acopy, T **A, dd_rowrange m, dd_colrange d,
                           dd_rowindex roworder) {
  dd_rowrange i;

  for (i = 1; i <= m; i++) {
    dd_CopyArow(Acopy[i - 1], A[roworder[i] - 1], d);
  }
}

template <typename T>
void dd_PermutePartialCopyAmatrix(T **Acopy, T **A, dd_rowrange m,
                                  dd_colrange d, dd_rowindex roworder) {
  /* copy the rows of A whose roworder is positive.  roworder[i] is the row
   * index of the copied row. */
  for (dd_rowrange i = 1; i <= m; i++)
    if (roworder[i] > 0)
      dd_CopyArow(Acopy[roworder[i] - 1], A[i - 1], d);
}

template <typename T> void dd_ColumnReduce(dd_conedata<T> *cone) {
  dd_colrange j, j1 = 0;
  dd_rowrange i;

  for (j = 1; j <= cone->d; j++) {
    if (cone->InitialRayIndex[j] > 0) {
      j1++;
      if (j1 < j) {
        for (i = 1; i <= cone->m; i++)
          cone->A[i - 1][j1 - 1] = cone->A[i - 1][j - 1];
        cone->newcol[j] = j1;
      }
    } else {
      cone->newcol[j] = 0;
    }
  }
  cone->d = j1; /* update the dimension. cone->d_orig remembers the old. */
  dd_CopyBmatrix(cone->d_orig, cone->B, cone->Bsave);
  /* save the dual basis inverse as Bsave.  This matrix contains the linearity
   * space generators. */
  cone->ColReduced = true;
}

template <typename T>
long dd_MatrixRank(dd_matrixdata<T> *M, dd_rowset ignoredrows,
                   dd_colset ignoredcols, dd_rowset *rowbasis,
                   dd_colset *colbasis) {
  bool stop, chosen;
  dd_rowset NopivotRow, PriorityRow;
  dd_colset ColSelected;
  T **B;
  dd_rowrange r;
  dd_colrange s;
  long rank;
  bool localdebug = false;
  std::vector<T> Rtemp(M->colsize);

  rank = 0;
  stop = false;
  set_initialize(&ColSelected, M->colsize);
  set_initialize(&NopivotRow, M->rowsize);
  set_initialize(rowbasis, M->rowsize);
  set_initialize(colbasis, M->colsize);
  set_initialize(&PriorityRow, M->rowsize);
  set_copy(NopivotRow, ignoredrows);
  set_copy(ColSelected, ignoredcols);
  dd_AllocateBmatrix(M->colsize, &B);
  dd_SetToIdentity(M->colsize, B);
  std::vector<long> roworder(M->rowsize + 1);
  for (r = 0; r < M->rowsize; r++)
    roworder[r + 1] = r + 1;
  roworder[M->rowsize] = 0;

  do { /* Find a set of rows for a basis */
    dd_SelectPivot2(M->rowsize, M->colsize, M->matrix, B, roworder.data(),
                    PriorityRow, M->rowsize, NopivotRow, ColSelected, &r, &s,
                    &chosen);
    if (localdebug && chosen)
      fprintf(stdout, "Procedure dd_MatrixRank: pivot on (r,s) =(%ld, %ld).\n",
              r, s);
    if (chosen) {
      set_addelem(NopivotRow, r);
      set_addelem(*rowbasis, r);
      set_addelem(ColSelected, s);
      set_addelem(*colbasis, s);
      rank++;
      //        std::cout << "dd_GaussianColumnPivot call 8\n";
      dd_GaussianColumnPivot(M->colsize, M->matrix, B, r, s, Rtemp.data());
    } else {
      stop = true;
    }
    if (rank == M->colsize)
      stop = true;
  } while (!stop);
  dd_FreeBmatrix(M->colsize, B);
  set_free(ColSelected);
  set_free(NopivotRow);
  set_free(PriorityRow);
  return rank;
}

template <typename T> void dd_FindBasis(dd_conedata<T> *cone, long *rank) {
  bool stop, chosen;
  dd_rowset NopivotRow;
  dd_colset ColSelected;
  dd_rowrange r;
  dd_colrange j, s;
  bool localdebug = false;
  std::vector<T> Rtemp(cone->d);

  *rank = 0;
  stop = false;
  for (j = 0; j <= cone->d; j++)
    cone->InitialRayIndex[j] = 0;
  set_emptyset(cone->InitialHalfspaces);
  set_initialize(&ColSelected, cone->d);
  set_initialize(&NopivotRow, cone->m);
  set_copy(NopivotRow, cone->NonequalitySet);
  dd_SetToIdentity(cone->d, cone->B);
  do { /* Find a set of rows for a basis */
    dd_SelectPivot2(cone->m, cone->d, cone->A, cone->B, cone->OrderVector,
                    cone->EqualitySet, cone->m, NopivotRow, ColSelected, &r, &s,
                    &chosen);
    if (localdebug && chosen)
      fprintf(stdout, "Procedure dd_FindBasis: pivot on (r,s) =(%ld, %ld).\n",
              r, s);
    if (chosen) {
      set_addelem(cone->InitialHalfspaces, r);
      set_addelem(NopivotRow, r);
      set_addelem(ColSelected, s);
      cone->InitialRayIndex[s] =
          r; /* cone->InitialRayIndex[s] stores the corr. row index */
      (*rank)++;
      //        std::cout << "dd_GaussianColumnPivot call 9\n";
      dd_GaussianColumnPivot(cone->d, cone->A, cone->B, r, s, Rtemp.data());
    } else {
      stop = true;
    }
    if (*rank == cone->d)
      stop = true;
  } while (!stop);
  set_free(ColSelected);
  set_free(NopivotRow);
}

template <typename T>
void dd_FindInitialRays(dd_conedata<T> *cone, bool *found) {
  dd_rowset CandidateRows;
  dd_rowrange i;
  long rank;
  dd_RowOrderType roworder_save = dd_LexMin;
  bool localdebug = false;

  *found = false;
  set_initialize(&CandidateRows, cone->m);
  if (cone->parent->InitBasisAtBottom == true) {
    roworder_save = cone->HalfspaceOrder;
    cone->HalfspaceOrder = dd_MaxIndex;
    cone->PreOrderedRun = false;
  } else {
    cone->PreOrderedRun = true;
  }
  for (i = 1; i <= cone->m; i++)
    if (!set_member(i, cone->NonequalitySet))
      set_addelem(CandidateRows, i);
  /*all rows not in NonequalitySet are candidates for initial cone*/
  dd_FindBasis(cone, &rank);
  cone->LinearityDim = cone->d - rank;
  if (localdebug)
    fprintf(stdout, "Linearity Dimension = %ld\n", cone->LinearityDim);
  if (cone->LinearityDim > 0) {
    dd_ColumnReduce(cone);
    dd_FindBasis(cone, &rank);
  }
  if (!set_subset(cone->EqualitySet, cone->InitialHalfspaces))
    if (localdebug)
      fprintf(
          stdout,
          "Equality set is dependent. Equality Set and an initial basis:\n");
  *found = true;
  set_free(CandidateRows);
  if (cone->parent->InitBasisAtBottom == true)
    cone->HalfspaceOrder = roworder_save;
  if (cone->HalfspaceOrder == dd_MaxCutoff ||
      cone->HalfspaceOrder == dd_MinCutoff ||
      cone->HalfspaceOrder == dd_MixCutoff) {
    cone->PreOrderedRun = false;
  } else {
    cone->PreOrderedRun = true;
  }
}

template <typename T>
void dd_CheckEquality(dd_colrange d_size, dd_raydata<T> **RP1,
                      dd_raydata<T> **RP2, bool *equal) {
  long j;
  bool localdebug = false;

  if (localdebug)
    fprintf(stdout, "Check equality of two rays\n");
  *equal = true;
  j = 1;
  while (j <= d_size && *equal) {
    if ((*RP1)->Ray[j - 1] != (*RP2)->Ray[j - 1])
      *equal = false;
    j++;
  }
  if (*equal)
    fprintf(stdout, "Equal records found !!!!\n");
}

template <typename T>
void dd_CreateNewRay(dd_conedata<T> *cone, dd_raydata<T> *Ptr1,
                     dd_raydata<T> *Ptr2, dd_rowrange ii) {
  /*Create a new ray by taking a linear combination of two rays*/
  dd_colrange j;
  T a1, a2, v1, v2;
  T *NewRay;
  NewRay = new T[cone->d];

  dd_AValue(&a1, cone->d, cone->A, Ptr1->Ray, ii);
  dd_AValue(&a2, cone->d, cone->A, Ptr2->Ray, ii);
  v1 = T_abs(a1);
  v2 = T_abs(a2);
  for (j = 0; j < cone->d; j++)
    NewRay[j] = Ptr1->Ray[j] * v2 + Ptr2->Ray[j] * v1;
  dd_AddRay(cone, NewRay);
  delete[] NewRay;
}

template <typename T>
void dd_EvaluateARay1(dd_rowrange i, dd_conedata<T> *cone)
/* Evaluate the ith component of the vector  A x RD.Ray
    and rearrange the linked list so that
    the infeasible rays with respect to  i  will be
    placed consecutively from First
 */
{
  dd_colrange j;
  T temp;
  dd_raydata<T> *Ptr;
  dd_raydata<T> *PrevPtr;
  dd_raydata<T> *TempPtr;

  Ptr = cone->FirstRay;
  PrevPtr = cone->ArtificialRay;
  if (PrevPtr->Next != Ptr) {
    std::cerr << "Error.  Artificial Ray does not point to FirstRay!!!\n";
    throw TerminalException{1};
  }
  while (Ptr != nullptr) {
    temp = 0;
    for (j = 0; j < cone->d; j++)
      temp += cone->A[i - 1][j] * Ptr->Ray[j];
    Ptr->ARay = temp;
    if (temp < 0 && Ptr != cone->FirstRay) {
      /* fprintf(stdout,"Moving an infeasible record w.r.t. %ld to
       * FirstRay\n",i); */
      if (Ptr == cone->LastRay)
        cone->LastRay = PrevPtr;
      TempPtr = Ptr;
      Ptr = Ptr->Next;
      PrevPtr->Next = Ptr;
      cone->ArtificialRay->Next = TempPtr;
      TempPtr->Next = cone->FirstRay;
      cone->FirstRay = TempPtr;
    } else {
      PrevPtr = Ptr;
      Ptr = Ptr->Next;
    }
  }
}

template <typename T>
void dd_EvaluateARay2(dd_rowrange i, dd_conedata<T> *cone)
/* Evaluate the ith component of the vector  A x RD.Ray
   and rearrange the linked list so that
   the infeasible rays with respect to  i  will be
   placed consecutively from First. Also for all feasible rays,
   "positive" rays and "zero" rays will be placed consecutively.
 */
{
  dd_colrange j;
  T temp;
  dd_raydata<T> *Ptr;
  dd_raydata<T> *NextPtr;
  bool zerofound = false, negfound = false, posfound = false;

  if (cone == nullptr || cone->TotalRayCount <= 0)
    return;
  cone->PosHead = nullptr;
  cone->ZeroHead = nullptr;
  cone->NegHead = nullptr;
  cone->PosLast = nullptr;
  cone->ZeroLast = nullptr;
  cone->NegLast = nullptr;
  Ptr = cone->FirstRay;
  while (Ptr != nullptr) {
    NextPtr = Ptr->Next; /* remember the Next record */
    Ptr->Next = nullptr; /* then clear the Next pointer */
    temp = 0;
    for (j = 0; j < cone->d; j++)
      temp += cone->A[i - 1][j] * Ptr->Ray[j];
    Ptr->ARay = temp;
    if (temp < 0) {
      if (!negfound) {
        negfound = true;
        cone->NegHead = Ptr;
        cone->NegLast = Ptr;
      } else {
        Ptr->Next = cone->NegHead;
        cone->NegHead = Ptr;
      }
    } else {
      if (temp > 0) {
        if (!posfound) {
          posfound = true;
          cone->PosHead = Ptr;
          cone->PosLast = Ptr;
        } else {
          Ptr->Next = cone->PosHead;
          cone->PosHead = Ptr;
        }
      } else {
        if (!zerofound) {
          zerofound = true;
          cone->ZeroHead = Ptr;
          cone->ZeroLast = Ptr;
        } else {
          Ptr->Next = cone->ZeroHead;
          cone->ZeroHead = Ptr;
        }
      }
    }
    Ptr = NextPtr;
  }
  /* joining three neg, pos and zero lists */
  if (negfound) { /* -list nonempty */
    cone->FirstRay = cone->NegHead;
    if (posfound) { /* -list & +list nonempty */
      cone->NegLast->Next = cone->PosHead;
      if (zerofound) { /* -list, +list, 0list all nonempty */
        cone->PosLast->Next = cone->ZeroHead;
        cone->LastRay = cone->ZeroLast;
      } else { /* -list, +list nonempty but  0list empty */
        cone->LastRay = cone->PosLast;
      }
    } else {           /* -list nonempty & +list empty */
      if (zerofound) { /* -list,0list nonempty & +list empty */
        cone->NegLast->Next = cone->ZeroHead;
        cone->LastRay = cone->ZeroLast;
      } else { /* -list nonempty & +list,0list empty */
        cone->LastRay = cone->NegLast;
      }
    }
  } else if (posfound) { /* -list empty & +list nonempty */
    cone->FirstRay = cone->PosHead;
    if (zerofound) { /* -list empty & +list,0list nonempty */
      cone->PosLast->Next = cone->ZeroHead;
      cone->LastRay = cone->ZeroLast;
    } else { /* -list,0list empty & +list nonempty */
      cone->LastRay = cone->PosLast;
    }
  } else { /* -list,+list empty & 0list nonempty */
    cone->FirstRay = cone->ZeroHead;
    cone->LastRay = cone->ZeroLast;
  }
  cone->ArtificialRay->Next = cone->FirstRay;
  cone->LastRay->Next = nullptr;
}

template <typename T>
void dd_DeleteNegativeRays(dd_conedata<T> *cone)
/* Eliminate the infeasible rays with respect to  i  which
   are supposed to be consecutive from the head of the dd_Ray list,
   and sort the zero list assumed to be consecutive at the
   end of the list.
 */
{
  dd_rowrange fii, fiitest;
  T temp;
  dd_raydata<T> *Ptr;
  dd_raydata<T> *PrevPtr;
  dd_raydata<T> *NextPtr;
  dd_raydata<T> *ZeroPtr1;
  dd_raydata<T> *ZeroPtr0;
  bool found, completed, zerofound = false, negfound = false, posfound = false;

  cone->PosHead = nullptr;
  cone->ZeroHead = nullptr;
  cone->NegHead = nullptr;
  cone->PosLast = nullptr;
  cone->ZeroLast = nullptr;
  cone->NegLast = nullptr;

  /* Delete the infeasible rays  */
  PrevPtr = cone->ArtificialRay;
  Ptr = cone->FirstRay;
  if (PrevPtr->Next != Ptr)
    fprintf(stdout, "Error at dd_DeleteNegativeRays: ArtificialRay does not "
                    "point the FirstRay.\n");
  completed = false;
  while (Ptr != nullptr && !completed) {
    if (Ptr->ARay < 0) {
      dd_Eliminate(cone, &PrevPtr);
      Ptr = PrevPtr->Next;
    } else {
      completed = true;
    }
  }

  /* Sort the zero rays */
  Ptr = cone->FirstRay;
  cone->ZeroRayCount = 0;
  while (Ptr != nullptr) {
    NextPtr = Ptr->Next; /* remember the Next record */
    temp = Ptr->ARay;
    if (temp < 0) {
      if (!negfound) {
        fprintf(stdout, "Error: An infeasible ray found after their removal\n");
        negfound = true;
      }
    } else {
      if (temp > 0) {
        if (!posfound) {
          posfound = true;
          cone->PosHead = Ptr;
          cone->PosLast = Ptr;
        } else {
          cone->PosLast = Ptr;
        }
      } else {
        (cone->ZeroRayCount)++;
        if (!zerofound) {
          zerofound = true;
          cone->ZeroHead = Ptr;
          cone->ZeroLast = Ptr;
          cone->ZeroLast->Next = nullptr;
        } else { /* Find a right position to store the record sorted w.r.t.
                    FirstInfeasIndex */
          fii = Ptr->FirstInfeasIndex;
          found = false;
          ZeroPtr1 = nullptr;
          for (ZeroPtr0 = cone->ZeroHead; !found && ZeroPtr0 != nullptr;
               ZeroPtr0 = ZeroPtr0->Next) {
            fiitest = ZeroPtr0->FirstInfeasIndex;
            if (fiitest >= fii) {
              found = true;
            } else {
              ZeroPtr1 = ZeroPtr0;
            }
          }
          /* fprintf(stdout,"insert position found \n %d  index %ld\n",found,
           * fiitest); */
          if (!found) { /* the new record must be stored at the end of list */
            cone->ZeroLast->Next = Ptr;
            cone->ZeroLast = Ptr;
            cone->ZeroLast->Next = nullptr;
          } else {
            if (ZeroPtr1 == nullptr) { /* store the new one at the head, and
                                          update the head ptr */
              /* fprintf(stdout,"Insert at the head\n"); */
              Ptr->Next = cone->ZeroHead;
              cone->ZeroHead = Ptr;
            } else { /* store the new one inbetween ZeroPtr1 and 0 */
              /* fprintf(stdout,"Insert inbetween\n");  */
              Ptr->Next = ZeroPtr1->Next;
              ZeroPtr1->Next = Ptr;
            }
          }
          /*
            Ptr->Next=cone->ZeroHead;
            cone->ZeroHead=Ptr;
          */
        }
      }
    }
    Ptr = NextPtr;
  }
  /* joining the pos and zero lists */
  if (posfound) { /* -list empty & +list nonempty */
    cone->FirstRay = cone->PosHead;
    if (zerofound) { /* +list,0list nonempty */
      cone->PosLast->Next = cone->ZeroHead;
      cone->LastRay = cone->ZeroLast;
    } else { /* 0list empty & +list nonempty */
      cone->LastRay = cone->PosLast;
    }
  } else { /* +list empty & 0list nonempty */
    cone->FirstRay = cone->ZeroHead;
    cone->LastRay = cone->ZeroLast;
  }
  cone->ArtificialRay->Next = cone->FirstRay;
  cone->LastRay->Next = nullptr;
}

template <typename T>
void dd_FeasibilityIndices(long *fnum, long *infnum, dd_rowrange i,
                           dd_conedata<T> *cone) {
  /*Evaluate the number of feasible rays and infeasible rays*/
  /*  w.r.t the hyperplane  i*/
  dd_colrange j;
  T temp;
  dd_raydata<T> *Ptr;

  *fnum = 0;
  *infnum = 0;
  Ptr = cone->FirstRay;
  while (Ptr != nullptr) {
    temp = 0;
    for (j = 0; j < cone->d; j++)
      temp += cone->A[i - 1][j] * Ptr->Ray[j];
    if (temp >= 0)
      (*fnum)++;
    else
      (*infnum)++;
    Ptr = Ptr->Next;
  }
}

template <typename T>
bool dd_LexEqual(T *v1, T *v2,
                 long dmax) { /* dmax is the size of vectors v1,v2 */
  bool determined, equal;
  dd_colrange j;

  equal = true;
  determined = false;
  j = 1;
  do {
    if (v1[j - 1] != v2[j - 1]) { /* 093c */
      equal = false;
      determined = true;
    } else {
      j++;
    }
  } while (!(determined) && (j <= dmax));
  return equal;
}

template <typename T>
void dd_AddNewHalfspace1(dd_conedata<T> *cone, dd_rowrange hnew)
/* This procedure 1 must be used with PreorderedRun=false
   This procedure is the most elementary implementation of
   DD and can be used with any type of ordering, including
   dynamic ordering of rows, e.g. MaxCutoff, MinCutoff.
   The memory requirement is minimum because it does not
   store any adjacency among the rays.
*/
{
  dd_raydata<T> *RayPtr0;
  dd_raydata<T> *RayPtr1;
  dd_raydata<T> *RayPtr2;
  dd_raydata<T> *RayPtr2s;
  dd_raydata<T> *RayPtr3;
  long pos1, pos2;
  double prevprogress, progress;
  T value1, value2;
  bool adj, equal, completed;
  bool localdebug = false;

  dd_EvaluateARay1(hnew, cone);
  /*Check feasibility of rays w.r.t. hnew
    and put all infeasible ones consecutively */

  RayPtr0 = cone->ArtificialRay; /*Pointer pointing RayPrt1*/
  RayPtr1 = cone->FirstRay; /*1st hnew-infeasible ray to scan and compare with
                               feasible rays*/
  value1 = cone->FirstRay->ARay;
  if (value1 >= 0) {
    if (cone->RayCount == cone->WeaklyFeasibleRayCount)
      cone->CompStatus = dd_AllFound;
    return; /* Sicne there is no hnew-infeasible ray and nothing to do */
  } else {
    RayPtr2s = RayPtr1->Next; /* RayPtr2s must point the first feasible ray */
    pos2 = 1;
    while (RayPtr2s != nullptr && RayPtr2s->ARay < 0) {
      RayPtr2s = RayPtr2s->Next;
      pos2++;
    }
  }
  if (RayPtr2s == nullptr) {
    cone->FirstRay = nullptr;
    cone->ArtificialRay->Next = cone->FirstRay;
    cone->RayCount = 0;
    cone->CompStatus = dd_AllFound;
    return; /* All rays are infeasible, and the computation must stop */
  }
  RayPtr2 = RayPtr2s;      /*2nd feasible ray to scan and compare with 1st*/
  RayPtr3 = cone->LastRay; /*Last feasible for scanning*/
  prevprogress = -10.0;
  pos1 = 1;
  completed = false;
  while ((RayPtr1 != RayPtr2s) && !completed) {
    value1 = RayPtr1->ARay;
    value2 = RayPtr2->ARay;
    dd_CheckEquality(cone->d, &RayPtr1, &RayPtr2, &equal);
    if ((value1 > 0 && value2 < 0) || (value1 < 0 && value2 > 0)) {
      dd_CheckAdjacency(cone, &RayPtr1, &RayPtr2, &adj);
      if (adj)
        dd_CreateNewRay(cone, RayPtr1, RayPtr2, hnew);
    }
    if (RayPtr2 != RayPtr3) {
      RayPtr2 = RayPtr2->Next;
      continue;
    }
    if (value1 < 0 || equal) {
      dd_Eliminate(cone, &RayPtr0);
      RayPtr1 = RayPtr0->Next;
      RayPtr2 = RayPtr2s;
    } else {
      completed = true;
    }
    pos1++;
    progress =
        100.0 * (static_cast<double>(pos1) / pos2) * (2.0 * pos2 - pos1) / pos2;
    if (progress - prevprogress >= 10 && pos1 % 10 == 0 && localdebug) {
      fprintf(
          stdout,
          "*Progress of iteration %5ld(/%ld):   %4ld/%4ld => %4.1f%% done\n",
          cone->Iteration, cone->m, pos1, pos2, progress);
      prevprogress = progress;
    }
  }
  if (cone->RayCount == cone->WeaklyFeasibleRayCount)
    cone->CompStatus = dd_AllFound;
}

template <typename T>
void dd_AddNewHalfspace2(dd_conedata<T> *cone, dd_rowrange hnew)
/* This procedure must be used under PreOrderedRun mode */
{
  bool localdebug = false;
  dd_raydata<T> *RayPtr1;
  dd_raydata<T> *RayPtr2;

  dd_adjacencydata<T> *EdgePtr, *EdgePtr0;
  dd_rowrange fii1, fii2;

  dd_EvaluateARay2(hnew, cone);
  /* Check feasibility of rays w.r.t. hnew
     and sort them. ( -rays, +rays, 0rays)*/

  if (cone->PosHead == nullptr && cone->ZeroHead == nullptr) {
    cone->FirstRay = nullptr;
    cone->ArtificialRay->Next = cone->FirstRay;
    cone->RayCount = 0;
    cone->CompStatus = dd_AllFound;
    return; /* All rays are infeasible, and the computation must stop */
  }
  if (cone->ZeroHead == nullptr)
    cone->ZeroHead = cone->LastRay;

  EdgePtr = cone->Edges[cone->Iteration];
  while (EdgePtr != nullptr) {
    RayPtr1 = EdgePtr->Ray1;
    RayPtr2 = EdgePtr->Ray2;
    fii1 = RayPtr1->FirstInfeasIndex;
    dd_CreateNewRay(cone, RayPtr1, RayPtr2, hnew);
    fii2 = cone->LastRay->FirstInfeasIndex;
    if (fii1 != fii2)
      dd_ConditionalAddEdge(cone, RayPtr1, cone->LastRay, cone->PosHead);
    EdgePtr0 = EdgePtr;
    EdgePtr = EdgePtr->Next;
    delete EdgePtr0;
    (cone->EdgeCount)--;
  }
  cone->Edges[cone->Iteration] = nullptr;

  dd_DeleteNegativeRays(cone);

  set_addelem(cone->AddedHalfspaces, hnew);

  if (cone->Iteration < cone->m) {
    if (cone->ZeroHead != nullptr && cone->ZeroHead != cone->LastRay) {
      if (cone->ZeroRayCount > 200 && localdebug)
        fprintf(stdout, "*New edges being scanned...\n");
      dd_UpdateEdges(cone, cone->ZeroHead, cone->LastRay);
    }
  }

  if (cone->RayCount == cone->WeaklyFeasibleRayCount)
    cone->CompStatus = dd_AllFound;
}

template <typename T>
dd_rowrange dd_SelectNextHalfspace0(dd_conedata<T> *cone, dd_rowset excluded) {
  /*A natural way to choose the next hyperplane.  Simply the largest index*/
  long i;
  bool determined;

  i = cone->m;
  determined = false;
  do {
    if (set_member(i, excluded))
      i--;
    else
      determined = true;
  } while (!determined && i >= 1);
  if (determined)
    return i;
  else
    return 0;
}

template <typename T>
dd_rowrange dd_SelectNextHalfspace1(dd_conedata<T> *cone, dd_rowset excluded) {
  /*Natural way to choose the next hyperplane.  Simply the least index*/
  long i;
  bool determined;

  i = 1;
  determined = false;
  do {
    if (set_member(i, excluded))
      i++;
    else
      determined = true;
  } while (!determined && i <= cone->m);
  if (determined)
    return i;
  else
    return 0;
}

template <typename T>
dd_rowrange dd_SelectNextHalfspace2(dd_conedata<T> *cone, dd_rowset excluded) {
  long i, fea, inf, infmin; /*feasibility and infeasibility numbers*/

  infmin = cone->RayCount + 1;
  dd_rowrange retidx = 0;
  for (i = 1; i <= cone->m; i++) {
    if (!set_member(i, excluded)) {
      dd_FeasibilityIndices(&fea, &inf, i, cone);
      if (inf < infmin) {
        infmin = inf;
        retidx = i;
      }
    }
  }
  return retidx;
}

template <typename T>
dd_rowrange dd_SelectNextHalfspace3(dd_conedata<T> *cone, dd_rowset excluded) {
  /*Choose the next hyperplane with maximum infeasibility*/
  long i, fea, inf, infmax; /*feasibility and infeasibility numbers*/

  infmax = -1;
  dd_rowrange retidx = 0;
  for (i = 1; i <= cone->m; i++)
    if (!set_member(i, excluded)) {
      dd_FeasibilityIndices(&fea, &inf, i, cone);
      if (inf > infmax) {
        infmax = inf;
        retidx = i;
      }
    }
  return retidx;
}

template <typename T>
dd_rowrange dd_SelectNextHalfspace4(dd_conedata<T> *cone, dd_rowset excluded) {
  long i, fea, inf, max, tmax, fi = 0, infi = 0;
  bool localdebug = false;

  max = -1;
  dd_rowrange retidx = 0;
  for (i = 1; i <= cone->m; i++) {
    if (!set_member(i, excluded)) {
      dd_FeasibilityIndices(&fea, &inf, i, cone);
      if (fea <= inf)
        tmax = inf;
      else
        tmax = fea;
      if (tmax > max) {
        max = tmax;
        fi = fea;
        infi = inf;
        retidx = i;
      }
    }
  }
  if (localdebug) {
    if (max == fi) {
      fprintf(stdout, "*infeasible rays (min) =%5ld, #feas rays =%5ld\n", infi,
              fi);
    } else {
      fprintf(stdout, "*infeasible rays (max) =%5ld, #feas rays =%5ld\n", infi,
              fi);
    }
  }
  return retidx;
}

template <typename T>
dd_rowrange dd_SelectNextHalfspace5(dd_conedata<T> *cone, dd_rowset excluded) {
  /*Choose the next hyperplane which is lexico-min*/
  long i, minindex;
  T *v1, *v2;

  minindex = 0;
  v1 = nullptr;
  for (i = 1; i <= cone->m; i++) {
    if (!set_member(i, excluded)) {
      v2 = cone->A[i - 1];
      if (minindex == 0) {
        minindex = i;
        v1 = v2;
      } else {
        if (dd_LexSmaller(v2, v1, cone->d)) {
          minindex = i;
          v1 = v2;
        }
      }
    }
  }
  return minindex;
}

template <typename T>
dd_rowrange dd_SelectNextHalfspace6(dd_conedata<T> *cone, dd_rowset excluded) {
  /*Choose the next hyperplane which is lexico-max*/
  long i, maxindex;
  T *v1, *v2;

  maxindex = 0;
  v1 = nullptr;
  for (i = 1; i <= cone->m; i++) {
    if (!set_member(i, excluded)) {
      v2 = cone->A[i - 1];
      if (maxindex == 0) {
        maxindex = i;
        v1 = v2;
      } else {
        if (dd_LexLarger(v2, v1, cone->d)) {
          maxindex = i;
          v1 = v2;
        }
      }
    }
  }
  return maxindex;
}

template <typename T>
void dd_UniqueRows(dd_rowindex OV, long p, long r, T **A, long dmax,
                   dd_rowset preferred, long *uniqrows) {
  /* Select a subset of rows of A (with range [p, q] up to dimension dmax) by
     removing duplicates. When a row subset preferred is nonempty, those row
     indices in the set have priority.  If two priority rows define the same row
     vector, one is chosen. For a selected unique row i, OV[i] returns a new
     position of the unique row i. For other nonuniqu row i, OV[i] returns a
     negative of the original row j dominating i. Thus the original contents of
     OV[p..r] will be rewritten.  Other components remain the same. *uniqrows
     returns the number of unique rows.
 */
  long i, iuniq, j;
  T *x;

  if (p <= 0 || p > r)
    return;
  iuniq = p;
  j = 1; /* the first unique row candidate */
  x = A[p - 1];
  OV[p] = j; /* tentative row index of the candidate */
  for (i = p + 1; i <= r; i++) {
    if (!dd_LexEqual(x, A[i - 1], dmax)) {
      /* a new row vector found. */
      iuniq = i;
      j = j + 1;
      OV[i] = j; /* Tentatively register the row i.  */
      x = A[i - 1];
    } else {
      /* rows are equal */
      if (set_member(i, preferred) && !set_member(iuniq, preferred)) {
        OV[iuniq] = -i; /* the row iuniq is dominated by the row i */
        iuniq = i;      /* the row i is preferred.  Change the candidate. */
        OV[i] = j;      /* the row i is tentatively registered. */
        x = A[i - 1];
      } else {
        OV[i] = -iuniq; /* the row iuniq is dominated by the row i */
      }
    }
  }
  *uniqrows = j;
}

template <typename T> void dd_ComputeRowOrderVector(dd_conedata<T> *cone) {
  long i, itemp;

  cone->OrderVector[0] = 0;
  switch (cone->HalfspaceOrder) {
  case dd_MaxIndex:
    for (i = 1; i <= cone->m; i++)
      cone->OrderVector[i] = cone->m - i + 1;
    break;

  case dd_MinIndex:
    for (i = 1; i <= cone->m; i++)
      cone->OrderVector[i] = i;
    break;

  case dd_LexMin:
  case dd_MinCutoff:
  case dd_MixCutoff:
  case dd_MaxCutoff:
    for (i = 1; i <= cone->m; i++)
      cone->OrderVector[i] = i;
    dd_RandomPermutation(cone->OrderVector, cone->m, cone->rseed);
    dd_QuickSort(cone->OrderVector, 1, cone->m, cone->A, cone->d);
    break;

  case dd_LexMax:
    for (i = 1; i <= cone->m; i++)
      cone->OrderVector[i] = i;
    dd_RandomPermutation(cone->OrderVector, cone->m, cone->rseed);
    dd_QuickSort(cone->OrderVector, 1, cone->m, cone->A, cone->d);
    for (i = 1; i <= cone->m / 2; i++) { /* just reverse the order */
      itemp = cone->OrderVector[i];
      cone->OrderVector[i] = cone->OrderVector[cone->m - i + 1];
      cone->OrderVector[cone->m - i + 1] = itemp;
    }
    break;

  case dd_RandomRow:
    for (i = 1; i <= cone->m; i++)
      cone->OrderVector[i] = i;
    dd_RandomPermutation(cone->OrderVector, cone->m, cone->rseed);
    break;
  }
}

template <typename T>
void dd_UpdateRowOrderVector(dd_conedata<T> *cone, dd_rowset PriorityRows)
/* Update the RowOrder vector to shift selected rows
in highest order.
*/
{
  dd_rowrange i, j, k, j1 = 0, oj = 0;
  long rr;
  bool found;

  found = true;
  rr = set_card(PriorityRows);
  for (i = 1; i <= rr; i++) {
    found = false;
    for (j = i; j <= cone->m && !found; j++) {
      oj = cone->OrderVector[j];
      if (set_member(oj, PriorityRows)) {
        found = true;
        j1 = j;
      }
    }
    if (found) {
      if (j1 > i) {
        /* shift everything lower: ov[i]->cone->ov[i+1]..ov[j1-1]->cone->ov[j1]
         */
        for (k = j1; k >= i; k--)
          cone->OrderVector[k] = cone->OrderVector[k - 1];
        cone->OrderVector[i] = oj;
      }
    } else {
      fprintf(stdout, "UpdateRowOrder: Error.\n");
      return;
    }
  }
}

template <typename T>
dd_rowrange dd_SelectPreorderedNext(dd_conedata<T> *cone, dd_rowset excluded) {
  dd_rowrange i, k;

  dd_rowrange retidx = 0;
  for (i = 1; i <= cone->m && retidx == 0; i++) {
    k = cone->OrderVector[i];
    if (!set_member(k, excluded))
      retidx = k;
  }
  return retidx;
}

template <typename T>
dd_rowrange dd_SelectNextHalfspace(dd_conedata<T> *cone, dd_rowset excluded) {
  if (cone->PreOrderedRun) {
    return dd_SelectPreorderedNext(cone, excluded);
  } else {
    switch (cone->HalfspaceOrder) {

    case dd_MaxIndex:
      return dd_SelectNextHalfspace0(cone, excluded);

    case dd_MinIndex:
      return dd_SelectNextHalfspace1(cone, excluded);

    case dd_MinCutoff:
      return dd_SelectNextHalfspace2(cone, excluded);

    case dd_MaxCutoff:
      return dd_SelectNextHalfspace3(cone, excluded);

    case dd_MixCutoff:
      return dd_SelectNextHalfspace4(cone, excluded);

    case dd_LexMin:
    case dd_LexMax:
    case dd_RandomRow:
      return dd_SelectNextHalfspace0(cone, excluded);
    }
  }
  // That case should never happen
  return -1;
}

template <typename T> void dd_DDInit(dd_conedata<T> *cone) {
  cone->Error = dd_NoError;
  cone->CompStatus = dd_InProgress;
  cone->RayCount = 0;
  cone->TotalRayCount = 0;
  cone->FeasibleRayCount = 0;
  cone->WeaklyFeasibleRayCount = 0;
  cone->EdgeCount = 0;      /* active edge count */
  cone->TotalEdgeCount = 0; /* active edge count */
  dd_SetInequalitySets(cone);
  dd_ComputeRowOrderVector(cone);
  cone->RecomputeRowOrder = false;
}

template <typename T> void dd_DDMain(dd_conedata<T> *cone) {
  dd_rowrange itemp, otemp;
  bool localdebug = false;

  auto clean = [&]() -> void {
    if (cone->d <= 0 ||
        cone->newcol[1] == 0) { /* fixing the number of output */
      cone->parent->n = cone->LinearityDim + cone->FeasibleRayCount - 1;
      cone->parent->ldim = cone->LinearityDim - 1;
    } else {
      cone->parent->n = cone->LinearityDim + cone->FeasibleRayCount;
      cone->parent->ldim = cone->LinearityDim;
    }
  };
  if (cone->d <= 0) {
    cone->Iteration = cone->m;
    cone->FeasibleRayCount = 0;
    cone->CompStatus = dd_AllFound;
    return clean();
  }
  while (cone->Iteration <= cone->m) {
    dd_rowrange hh = dd_SelectNextHalfspace(cone, cone->WeaklyAddedHalfspaces);
    if (set_member(hh, cone->NonequalitySet)) { /* Skip the row hh */
      if (localdebug) {
        fprintf(stdout,
                "*The row # %3ld should be inactive and thus skipped.\n", hh);
      }
      set_addelem(cone->WeaklyAddedHalfspaces, hh);
    } else {
      if (cone->PreOrderedRun) {
        dd_AddNewHalfspace2(cone, hh);
      } else {
        dd_AddNewHalfspace1(cone, hh);
      }
      set_addelem(cone->AddedHalfspaces, hh);
      set_addelem(cone->WeaklyAddedHalfspaces, hh);
    }
    if (!cone->PreOrderedRun) {
      otemp = -400;
      for (itemp = 1; cone->OrderVector[itemp] != hh; itemp++)
        otemp = cone->OrderVector[cone->Iteration];
      cone->OrderVector[cone->Iteration] = hh;
      /* store the dynamic ordering in ordervec */
      cone->OrderVector[itemp] = otemp;
      /* store the dynamic ordering in ordervec */
    }
    if (localdebug) {
      fprintf(stdout,
              "(Iter, Row, #Total, #Curr, #Feas)= %5ld %5ld %9ld %6ld %6ld\n",
              cone->Iteration, hh, cone->TotalRayCount, cone->RayCount,
              cone->FeasibleRayCount);
    }
    if (cone->CompStatus == dd_AllFound || cone->CompStatus == dd_RegionEmpty) {
      set_addelem(cone->AddedHalfspaces, hh);
      return clean();
    }
    (cone->Iteration)++;
  }
}

template <typename T> void dd_InitialDataSetup(dd_conedata<T> *cone) {
  long j, r;
  dd_rowset ZSet;

  std::vector<T> Vector1(cone->d);
  std::vector<T> Vector2(cone->d);

  cone->RecomputeRowOrder = false;
  cone->ArtificialRay = nullptr;
  cone->FirstRay = nullptr;
  cone->LastRay = nullptr;
  set_initialize(&ZSet, cone->m);
  dd_AddArtificialRay(cone);
  set_copy(cone->AddedHalfspaces, cone->InitialHalfspaces);
  set_copy(cone->WeaklyAddedHalfspaces, cone->InitialHalfspaces);
  dd_UpdateRowOrderVector(cone, cone->InitialHalfspaces);
  for (r = 1; r <= cone->d; r++) {
    for (j = 0; j < cone->d; j++) {
      Vector1[j] = cone->B[j][r - 1];
      Vector2[j] = -cone->B[j][r - 1];
    }
    dd_ZeroIndexSet(cone->m, cone->d, cone->A, Vector1.data(), ZSet);
    if (set_subset(cone->EqualitySet, ZSet)) {
      dd_AddRay(cone, Vector1.data());
      if (cone->InitialRayIndex[r] == 0)
        dd_AddRay(cone, Vector2.data());
    }
  }
  dd_CreateInitialEdges(cone);
  cone->Iteration = cone->d + 1;
  if (cone->Iteration > cone->m)
    cone->CompStatus = dd_AllFound; /* 0.94b  */
  set_free(ZSet);
}

template <typename T>
bool dd_CheckEmptiness(dd_polyhedradata<T> *poly, dd_ErrorType *err, size_t const& maxiter, std::ostream& os) {
  bool answer = false;

  *err = dd_NoError;

  if (poly->representation == dd_Inequality) {
    dd_matrixdata<T> *M = dd_CopyInequalities(poly);
    dd_rowset R, S;
    set_initialize(&R, M->rowsize);
    set_initialize(&S, M->rowsize);
    if (!dd_ExistsRestrictedFace(M, R, S, err, maxiter, os)) {
      poly->child->CompStatus = dd_AllFound;
      poly->IsEmpty = true;
      poly->n = 0;
      answer = true;
    }
    set_free(R);
    set_free(S);
    dd_FreeMatrix(M);
  } else {
    if (poly->representation == dd_Generator && poly->m <= 0) {
      *err = dd_EmptyVrepresentation;
      poly->IsEmpty = true;
      poly->child->CompStatus = dd_AllFound;
      answer = true;
      poly->child->Error = *err;
    }
  }

  return answer;
}

template <typename T>
bool dd_DoubleDescription(dd_polyhedradata<T> *poly, dd_ErrorType *err, size_t const& maxiter, std::ostream& os) {
  dd_conedata<T> *cone = nullptr;
  bool found = false;

  *err = dd_NoError;

  if (poly != nullptr &&
      (poly->child == nullptr || poly->child->CompStatus != dd_AllFound)) {
    cone = dd_ConeDataLoad(poly);
    /* create a cone associated with poly by homogenization */
    dd_DDInit(cone);
    if (poly->representation == dd_Generator && poly->m <= 0) {
      *err = dd_EmptyVrepresentation;
      cone->Error = *err;
      return found;
    }
    /* Check emptiness of the polyhedron */
    dd_CheckEmptiness(poly, err, maxiter, os);

    if (cone->CompStatus != dd_AllFound) {
      dd_FindInitialRays(cone, &found);
      if (found) {
        dd_InitialDataSetup(cone);
        if (cone->CompStatus == dd_AllFound)
          return found;
        dd_DDMain(cone);
        if (cone->FeasibleRayCount != cone->RayCount)
          *err = dd_NumericallyInconsistent; /* cddlib-093d */
      }
    }
  }
  return found;
}

template <typename T>
bool dd_DoubleDescription2(dd_polyhedradata<T> *poly, dd_RowOrderType horder,
                           dd_ErrorType *err, size_t const& maxiter, std::ostream& os) {
  dd_conedata<T> *cone = nullptr;
  bool found = false;

  *err = dd_NoError;

  if (poly != nullptr &&
      (poly->child == nullptr || poly->child->CompStatus != dd_AllFound)) {
    cone = dd_ConeDataLoad(poly);
    // create a cone associated with poly by homogenization
    cone->HalfspaceOrder = horder;
    dd_DDInit(cone);
    if (poly->representation == dd_Generator && poly->m <= 0) {
      *err = dd_EmptyVrepresentation;
      cone->Error = *err;
      return found;
    }
    // Check emptiness of the polyhedron
    dd_CheckEmptiness(poly, err, maxiter, os);

    if (cone->CompStatus != dd_AllFound) {
      dd_FindInitialRays(cone, &found);
      if (found) {
        dd_InitialDataSetup(cone);
        if (cone->CompStatus == dd_AllFound)
          return found;
        dd_DDMain(cone);
        if (cone->FeasibleRayCount != cone->RayCount)
          *err = dd_NumericallyInconsistent;
      }
    }
  }
  return found;
}

template <typename T>
bool dd_DDInputAppend(dd_polyhedradata<T> **poly, dd_matrixdata<T> *M,
                      dd_ErrorType *err, size_t const& maxiter, std::ostream& os) {
  /* This is imcomplete.  It simply solves the problem from scratch.  */

  if ((*poly)->child != nullptr)
    dd_FreeDDMemory(*poly);
  dd_AppendMatrix2Poly(poly, M);
  (*poly)->representation = dd_Inequality;
  return dd_DoubleDescription(*poly, err, maxiter, os);
}

template <typename T>
dd_matrixdata<T> *MyMatrix_PolyFile2Matrix(MyMatrix<T> const &TheEXT) {
  dd_matrixdata<T> *M = nullptr;
  dd_rowrange m_input, i;
  dd_colrange d_input, j;
  dd_RepresentationType rep;
  bool localdebug = false;

  m_input = TheEXT.rows();
  d_input = TheEXT.cols();

  rep = dd_Generator; /* using dd_Inequality led to horrible bugs */
  M = dd_CreateMatrix<T>(m_input, d_input);
  M->representation = rep;

  for (i = 0; i < m_input; i++)
    for (j = 0; j < d_input; j++) {
      M->matrix[i][j] = TheEXT(i, j);
      if (localdebug)
        std::cout << "i=" << i << " j=" << j << " value=" << TheEXT(i, j)
                  << "\n";
    }
  return M;
}

template <typename T>
dd_matrixdata<T> *MyMatrix_PolyFile2MatrixExt(MyMatrix<T> const &TheEXT) {
  dd_matrixdata<T> *M = nullptr;
  dd_rowrange m_input, i;
  dd_colrange d_input, j;
  dd_RepresentationType rep;
  bool localdebug = false;

  m_input = TheEXT.rows();
  d_input = TheEXT.cols();

  rep = dd_Generator; /* using dd_Inequality led to horrible bugs */
  M = dd_CreateMatrix<T>(m_input, d_input + 1);
  M->representation = rep;

  for (i = 0; i < m_input; i++) {
    M->matrix[i][0] = 0;
    for (j = 0; j < d_input; j++) {
      M->matrix[i][j + 1] = TheEXT(i, j);
      if (localdebug)
        std::cout << "i=" << i << " j=" << j << " value=" << TheEXT(i, j)
                  << "\n";
    }
  }
  return M;
}

template <typename T>
MyMatrix<T> FAC_from_poly(dd_polyhedradata<T> const *poly, int const &nbCol) {
  dd_raydata<T> *RayPtr = poly->child->FirstRay;
  int nbRay = 0;
  while (RayPtr != nullptr) {
    if (RayPtr->feasible)
      nbRay++;
    RayPtr = RayPtr->Next;
  }
  MyMatrix<T> TheFAC(nbRay, nbCol);
  int iRay = 0;
  RayPtr = poly->child->FirstRay;
  while (RayPtr != nullptr) {
    if (RayPtr->feasible) {
      for (int iCol = 0; iCol < nbCol; iCol++)
        TheFAC(iRay, iCol) = RayPtr->Ray[iCol];
      iRay++;
    }
    RayPtr = RayPtr->Next;
  }
  return TheFAC;
}

template <typename T>
vectface ListIncd_from_poly(dd_polyhedradata<T> const *poly,
                            MyMatrix<T> const &EXT) {
  size_t nbRow = EXT.rows();
  vectface ListIncd(nbRow);
  dd_raydata<T> *RayPtr = poly->child->FirstRay;
  T eScal;
  Face f(nbRow);
  while (RayPtr != nullptr) {
    if (RayPtr->feasible) {
      for (size_t iRow = 0; iRow < nbRow; iRow++) {
        long elem = iRow + 1;
        f[iRow] = set_member(elem, RayPtr->ZeroSet);
      }
      ListIncd.push_back(f);
    }
    RayPtr = RayPtr->Next;
  }
  return ListIncd;
}

template <typename T, typename Fprocess>
void ListFaceIneq_from_poly(dd_polyhedradata<T> const *poly,
                            MyMatrix<T> const &EXT, Fprocess f_process) {
  size_t nbCol = EXT.cols();
  size_t nbRow = EXT.rows();
  dd_raydata<T> *RayPtr = poly->child->FirstRay;
  T eScal;
  std::pair<Face, MyVector<T>> pair{Face(nbRow), MyVector<T>(nbCol)};
  while (RayPtr != nullptr) {
    if (RayPtr->feasible) {
      for (size_t iRow = 0; iRow < nbRow; iRow++) {
        long elem = iRow + 1;
        pair.first[iRow] = set_member(elem, RayPtr->ZeroSet);
      }
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        pair.second(iCol) = RayPtr->Ray[iCol];
      f_process(pair);
    }
    RayPtr = RayPtr->Next;
  }
}

template <typename T>
std::vector<int> RedundancyReductionClarkson(MyMatrix<T> const &TheEXT,
                                             std::ostream &os) {
  dd_ErrorType err = dd_NoError;
  int nbRow = TheEXT.rows();
  dd_matrixdata<T> *M = MyMatrix_PolyFile2Matrix(TheEXT);
  M->representation = dd_Inequality;
  //  M->representation = dd_Generator;
  size_t maxiter = 0;
  dd_rowset redset = dd_RedundantRowsViaShooting(M, &err, maxiter, os);
#ifdef DEBUG_CDD
  os << "CDD: RedundancyReductionClarkson err=" << err << "\n";
#endif
  if (err != dd_NoError) {
    std::cerr << "TheEXT=\n";
    WriteMatrix(std::cerr, TheEXT);
    std::cerr << "RedundancyReductionClarkson internal CDD error\n";
    throw TerminalException{1};
  }
  std::vector<int> ListIdx;
  for (int i_row = 0; i_row < nbRow; i_row++) {
    bool isin = set_member(i_row + 1, redset);
    if (!isin)
      ListIdx.push_back(i_row);
  }
  dd_FreeMatrix(M);
  set_free(redset);
  return ListIdx;
}

template <typename T>
std::vector<int>
RedundancyReductionClarksonBlocks(MyMatrix<T> const &TheEXT,
                                  std::vector<int> const &BlockBelong,
                                  [[maybe_unused]] std::ostream &os) {
  dd_ErrorType err = dd_NoError;
  int nbRow = TheEXT.rows();
  dd_matrixdata<T> *M = MyMatrix_PolyFile2Matrix(TheEXT);
  M->representation = dd_Inequality;
  //  M->representation = dd_Generator;
  size_t maxiter = 0;
  dd_rowset redset = dd_RedundantRowsViaShootingBlocks(M, &err, BlockBelong, maxiter, os);
  if (err != dd_NoError) {
    std::cerr << "RedundancyReductionClarksonBlocks internal CDD error\n";
    throw TerminalException{1};
  }
  std::vector<int> ListIdx;
  for (int i_row = 0; i_row < nbRow; i_row++) {
    bool isin = set_member(i_row + 1, redset);
    if (!isin)
      ListIdx.push_back(i_row);
  }
  dd_FreeMatrix(M);
  set_free(redset);
  return ListIdx;
}

template <typename T>
std::pair<MyMatrix<T>, Face>
KernelLinearDeterminedByInequalitiesAndIndices_DirectLP(
    MyMatrix<T> const &FAC, size_t const& maxiter, std::ostream &os) {
  dd_ErrorType err = dd_NoError;
  int nbRow = FAC.rows();
  int nbCol = FAC.cols();
  dd_matrixdata<T> *M = MyMatrix_PolyFile2MatrixExt(FAC);
  M->representation = dd_Inequality;
  dd_rowset linset = dd_ImplicitLinearityRows(M, &err, maxiter, os);
  if (err != dd_NoError) {
    std::cerr << "DualDescription_incd internal CDD error\n";
    throw TerminalException{1};
  }
  Face f(nbRow);
  std::vector<MyVector<T>> ListV;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    long elem = iRow + 1;
    bool test = set_member(elem, linset);
    f[iRow] = test;
    if (test) {
      MyVector<T> eFAC = GetMatrixRow(FAC, iRow);
      ListV.push_back(eFAC);
    }
  }
#ifdef DEBUG_CDD
  os << "CDD: |ListV|=" << ListV.size() << "\n";
#endif
  auto get_nsp = [&]() -> MyMatrix<T> {
    if (ListV.size() > 0) {
      MyMatrix<T> MatEqua = MatrixFromVectorFamily(ListV);
      return NullspaceTrMat(MatEqua);
    } else {
      return IdentityMat<T>(nbCol);
    }
  };
  MyMatrix<T> NSP = get_nsp();
  dd_FreeMatrix(M);
  set_free(linset);
  return {std::move(NSP), f};
}

/*
  We apply the preceding algorithm. But when we find an equality, we use it
  to reduce the dimensionality.
  ---
  The reasoning of the method is that any gain in dimensionality has to
  exploited to the full. Note that the gain is not just about the dimensionality
  but also the number of inequalities since many inequalities become multiple
  of each other when going to a lower dimensional subspace.
 */
template <typename T>
std::pair<MyMatrix<T>, Face>
KernelLinearDeterminedByInequalitiesAndIndices_LPandNullspace(
    MyMatrix<T> const &FAC, std::ostream &os) {
  int nbRow = FAC.rows();
  int nbCol = FAC.cols();
  dd_matrixdata<T> *M = MyMatrix_PolyFile2MatrixExt(FAC);
  M->representation = dd_Inequality;
  size_t maxiter = 0;
  int d1 = M->colsize; // We are in the inequality case.
  std::vector<T> cvec(d1);
#ifdef DEBUG_CDD
  os << "CDD: LPandNullspace, nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
#endif
  auto get_linear_entry = [&]() -> std::optional<int> {
    dd_ErrorType err = dd_NoError;
    for (int iRow = 0; iRow < nbRow; iRow++) {
      long i = iRow + 1;
      bool test = dd_ImplicitLinearity(M, i, cvec.data(), &err, maxiter, os);
      if (err != dd_NoError) {
        std::cerr << "LinearDeterminedByInequalitiesAndIndices_LPandNullspace "
                     "internal CDD error\n";
        throw TerminalException{1};
      }
      if (test) {
        return iRow;
      }
    }
    return {};
  };
  std::optional<int> opt = get_linear_entry();
  if (opt) {
    int idx = *opt;
#ifdef DEBUG_CDD
    os << "CDD: LPandNullspace, idx=" << idx << "\n";
#endif
    MyVector<T> Vlin = GetMatrixRow(FAC, idx);
    MyMatrix<T> NSP = NullspaceMatSingleVector(Vlin);
#ifdef DEBUG_CDD
    os << "CDD: Vlin=";
    WriteVectorNoDim(os, Vlin);
    os << "CDD: FAC=\n";
    WriteMatrix(os, FAC);
    os << "CDD: NSP=\n";
    WriteMatrix(os, NSP);
#endif
    std::vector<int> l_zeros;
    std::unordered_map<MyVector<T>, std::vector<int>> map;
    for (int iRow = 0; iRow < nbRow; iRow++) {
      MyVector<T> V = GetMatrixRow(FAC, iRow);
      MyVector<T> ScalV = NSP * V;
#ifdef DEBUG_CDD
      os << "CDD:    ScalV=";
      WriteVectorNoDim(os, ScalV);
#endif
      if (IsZeroVector(ScalV)) {
        l_zeros.push_back(iRow);
      } else {
        MyVector<T> ScalVcan = CanonicalizeVector(ScalV);
#ifdef DEBUG_CDD
        os << "CDD: ScalVcan=";
        WriteVectorNoDim(os, ScalVcan);
#endif
        map[ScalVcan].push_back(iRow);
      }
    }
    int siz = map.size();
#ifdef DEBUG_CDD
    os << "CDD: siz=" << siz << " nbRow=" << nbRow << "\n";
#endif
    MyMatrix<T> FACred(siz, nbCol - 1);
    std::vector<std::vector<int>> ll_idx(siz);
    int pos = 0;
    for (auto &kv : map) {
      for (int u = 0; u < nbCol - 1; u++) {
        FACred(pos, u) = kv.first(u);
      }
      ll_idx[pos] = kv.second;
      pos++;
    }
#ifdef DEBUG_CDD
    os << "CDD: FACred=\n";
    WriteMatrix(os, FACred);
    os << "CDD: |FACred|=" << FACred.rows() << " / " << FACred.cols() << "\n";
    if (FACred.rows() > 0) {
      size_t maxiter = 0;
      std::pair<MyMatrix<T>, Face> pair_loc =
        KernelLinearDeterminedByInequalitiesAndIndices_DirectLP(FACred, maxiter, os);
      os << "CDD: pair_loc.first=\n";
      WriteMatrix(os, pair_loc.first);
    }
#endif
    std::pair<MyMatrix<T>, Face> pairRed =
        KernelLinearDeterminedByInequalitiesAndIndices_LPandNullspace(FACred,
                                                                      os);
#ifdef DEBUG_CDD
    os << "CDD: We have pairRed\n";
#endif
    MyMatrix<T> NSPnew = pairRed.first * NSP;
    Face f(nbRow);
    for (auto &idx : l_zeros) {
      f[idx] = 1;
    }
    for (int u = 0; u < siz; u++) {
      if (pairRed.second[u]) {
        for (auto &pos : ll_idx[u]) {
          f[pos] = 1;
        }
      }
    }
#ifdef DEBUG_CDD
    os << "CDD: We have NSPnew / f\n";
#endif
    return {std::move(NSPnew), std::move(f)};
  } else {
#ifdef DEBUG_CDD
    os << "CDD: None case\n";
#endif
    MyMatrix<T> Spa = IdentityMat<T>(nbCol);
    Face f(nbRow);
    return {std::move(Spa), std::move(f)};
  }
}

template <typename T>
MyMatrix<T> DualDescription(MyMatrix<T> const &TheEXT, std::ostream& os) {
  dd_ErrorType err = dd_NoError;
  int nbCol = TheEXT.cols();
  dd_matrixdata<T> *M = MyMatrix_PolyFile2Matrix(TheEXT);
  size_t maxiter = 0;
  dd_polyhedradata<T> *poly = dd_DDMatrix2Poly(M, &err, maxiter, os);
  if (err != dd_NoError) {
    std::cerr << "DualDescription internal CDD error\n";
    throw TerminalException{1};
  }
  MyMatrix<T> TheFAC = FAC_from_poly(poly, nbCol);
  dd_FreePolyhedra(poly);
  dd_FreeMatrix(M);
  return TheFAC;
}

template <typename T>
vectface DualDescription_incd(MyMatrix<T> const &TheEXT, std::ostream& os) {
  dd_ErrorType err = dd_NoError;
  dd_matrixdata<T> *M = MyMatrix_PolyFile2Matrix(TheEXT);
  size_t maxiter = 0;
  dd_polyhedradata<T> *poly = dd_DDMatrix2Poly(M, &err, maxiter, os);
  if (err != dd_NoError) {
    std::cerr << "DualDescription_incd internal CDD error\n";
    throw TerminalException{1};
  }
  vectface ListIncd = ListIncd_from_poly(poly, TheEXT);
  dd_FreePolyhedra(poly);
  dd_FreeMatrix(M);
  return ListIncd;
}

template <typename T, typename Fprocess>
void DualDescriptionFaceIneq(MyMatrix<T> const &TheEXT, Fprocess f_process, std::ostream& os) {
  dd_ErrorType err = dd_NoError;
  dd_matrixdata<T> *M = MyMatrix_PolyFile2Matrix(TheEXT);
  size_t maxiter = 0;
  dd_polyhedradata<T> *poly = dd_DDMatrix2Poly(M, &err, maxiter, os);
  if (err != dd_NoError) {
    std::cerr << "DualDescriptionFaceIneq internal CDD error\n";
    throw TerminalException{1};
  }
  ListFaceIneq_from_poly(poly, TheEXT, f_process);
  dd_FreePolyhedra(poly);
  dd_FreeMatrix(M);
}

// clang-format off
}  // namespace cdd
// clang-format on

template <typename T>
void WriteInputFileCdd(std::string const &FileName, MyMatrix<T> const &ListIneq,
                       MyVector<T> const &ToBeMinimized) {
  int nbRow = ListIneq.rows();
  int nbCol = ListIneq.cols();
  std::ofstream os(FileName);
  os << "H-representation\n";
  os << "begin\n";
  os << " " << nbRow << " " << nbCol << " integer\n";
  for (int iRow = 0; iRow < nbRow; iRow++) {
    for (int iCol = 0; iCol < nbCol; iCol++)
      os << " " << ListIneq(iRow, iCol);
    os << "\n";
  }
  os << "end\n";
  os << "minimize\n";
  WriteVectorNoDim(os, ToBeMinimized);
}

template <typename T>
LpSolution<T>
CDD_LinearProgramming_External(MyMatrix<T> const &InequalitySet,
                               MyVector<T> const &ToBeMinimized,
                               [[maybe_unused]] std::ostream &os) {
  std::cerr << "Begin CDD_LinearProgramming_External\n";
  std::string eStr = random_string(20);
  std::string FileIne = "/tmp/LP_" + eStr + ".ine";
  std::string FileLps = "/tmp/LP_" + eStr + ".lps";
  std::string FileErr = "/tmp/LP_" + eStr + ".error";
  std::string FileDdl = "/tmp/LP_" + eStr + ".ddl";
  std::string FileLog = "/tmp/LP_" + eStr + ".log";
  std::string FileCpp = "/tmp/LP_" + eStr + ".cpp";
  auto CleanFile = [&]() -> void {
    RemoveFileIfExist(FileIne);
    RemoveFileIfExist(FileLps);
    RemoveFileIfExist(FileErr);
    RemoveFileIfExist(FileDdl);
    RemoveFileIfExist(FileLog);
    RemoveFileIfExist(FileCpp);
  };
  CleanFile();
  //
  WriteInputFileCdd(FileIne, InequalitySet, ToBeMinimized);
  //
  std::string FileTestlp2 = "testlp2_gmp";
  std::string eComm1 =
      FileTestlp2 + " " + FileIne + " 2> " + FileErr + " > " + FileLog;
  int iret1 = system(eComm1.c_str());
  if (iret1 != 0) {
    std::cerr << "iret1=" << iret1 << "\n";
    std::cerr << "Call to " << FileTestlp2 << " failed\n";
    std::cerr << "eComm1=" << eComm1 << "\n";
    throw TerminalException{1};
  }
  std::string FilelpcddcleanerCpp = "lpcddcleanerCpp";
  std::string eComm2 = FilelpcddcleanerCpp + " < " + FileLog + " > " + FileCpp;
  int iret2 = system(eComm2.c_str());
  if (iret2 != 0) {
    std::cerr << "iret2=" << iret2 << "\n";
    std::cerr << "Call to " << FilelpcddcleanerCpp << "\n";
    std::cerr << "eComm2=" << eComm2 << "\n";
    throw TerminalException{1};
  }
  //
  std::ifstream is(FileCpp);
  LpSolution<T> eSol;
  bool PrimalDefined;
  is >> PrimalDefined;
  eSol.PrimalDefined = PrimalDefined;
  //
  bool DualDefined;
  is >> DualDefined;
  eSol.DualDefined = DualDefined;
  //
  MyVector<T> DualSolution = ReadVector<T>(is);
  eSol.DualSolution = DualSolution;
  //
  T OptimalValue;
  is >> OptimalValue;
  eSol.OptimalValue = OptimalValue;
  //
  MyVector<T> DirectSolution = ReadVector<T>(is);
  eSol.DirectSolution = DirectSolution;
  //
  CleanFile();
  return eSol;
}

template <typename T>
std::optional<LpSolution<T>> GetLpSolutionFromLpData(
    MyMatrix<T> const &EXT, [[maybe_unused]] MyVector<T> const &eVect,
    cdd::dd_lpdata<T> *lp, [[maybe_unused]] std::ostream &os) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  cdd::dd_colrange j;
  int idx;
  bool PrimalDefined = false;
  bool DualDefined = false;
  switch (lp->LPS) {
  case cdd::dd_Optimal:
    // Correspondence in cddlp.c: We have dual_solution, primal_solution and
    // optimal_value
    PrimalDefined = true;
    DualDefined = true;
    break;
  case cdd::dd_Inconsistent:
    // Corrrespondence in cddlp.c: We have dual_direction
    PrimalDefined = false;
    DualDefined = true;
    break;
  case cdd::dd_DualInconsistent:
    // Correspondence in cddlp.c: Nothing actually.
    PrimalDefined = true;
    DualDefined = false;
    break;
  case cdd::dd_StrucDualInconsistent:
    // Correspondence in cddlp.c: We have primal_direction
    PrimalDefined = true;
    DualDefined = false;
    break;
  case cdd::dd_LPSundecided:
    // Programming of this case is missing. Please work from cdd code
#ifdef DEBUG_CDD
    os << "Programming of the case dd_LPSundecided is missing\n";
#endif
    return {};
  case cdd::dd_StrucInconsistent:
#ifdef DEBUG_CDD
    os << "Programming of the case dd_StructInconsistent is missing\n";
#endif
    return {};
  case cdd::dd_Unbounded:
#ifdef DEBUG_CDD
    os << "Programming of the case dd_Unbounded is missing\n";
#endif
    return {};
  case cdd::dd_DualUnbounded:
#ifdef DEBUG_CDD
    os << "Programming of the case dd_Unbounded is missing\n";
#endif
    return {};
  case cdd::dd_TooManyIterations:
    std::cerr << "CDD: Too many iterations done\n";
    throw TerminalException{1};
  }
  LpSolution<T> eSol;
  if (PrimalDefined) {
    MyVector<T> DirectSolution(nbCol - 1);
    for (j = 1; j < lp->d; j++) {
      DirectSolution(j - 1) = lp->sol[j];
    }
    eSol.PrimalDefined = true;
    eSol.DirectSolution = DirectSolution;
  }
  MyVector<T> DualSolution = ZeroVector<T>(nbRow);
  if (DualDefined) {
    for (j = 1; j < lp->d; j++) {
      idx = lp->nbindex[j + 1];
      if (idx > 0) {
        DualSolution(idx - 1) = lp->dsol[j];
      }
    }
    eSol.DualDefined = true;
    eSol.DualSolution = DualSolution;
  }
  if (PrimalDefined && DualDefined) {
    eSol.OptimalValue = lp->optvalue;
  }
#ifdef SANITY_CHECK_CDD
  if (PrimalDefined && DualDefined) {
    bool test = CheckDualSolutionGetOptimal(EXT, eVect, eSol, os);
    if (!test) {
      std::cerr << "This is not a valid dual solution\n";
      throw TerminalException{1};
    }
  }
#endif
  return eSol;
}

template<typename T>
bool is_lifting_possible(cdd::dd_lpdata<T> *lp) {
  if (lp->LPS != cdd::dd_Optimal) {
    // If we do not have an optimal solution, then we cannot have a configuration
    // of inequalities that get us the optimal vertex.
    return false;
  }
  cdd::dd_colrange j, d = lp->d;
  for (j = 1; j < d; j++) {
    long idx = lp->nbindex[j + 1];
    if (idx <= 0) {
      // A negative index means that we do not have a full configuration of vectors
      // and so the solution of the linear system strategy will not work.
      return false;
    }
  }
  return true;
}


template <typename T, typename Tfloat>
std::optional<LpSolution<T>>
LiftFloatingPointSolution(MyMatrix<T> const &EXT, MyVector<T> const &eVect,
                          cdd::dd_lpdata<Tfloat> *lp,
                          [[maybe_unused]] std::ostream &os) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  if (lp->LPS != cdd::dd_Optimal) {
#ifdef DEBUG_CDD
    os << "CDD: LIFT ERROR, We did not get an optinal solution. Therefore we "
          "cannot lift solution\n";
#endif
    return {};
  }
  // Now using the dual solution to find the exact vertex.
  cdd::dd_colrange i, j;
  cdd::dd_colrange d = lp->d;
  MyVector<T> V(d - 1);
  MyMatrix<T> M(d - 1, d - 1);
#ifdef DEBUG_CDD
  os << "CDD: Definition of V and M (step 1)\n";
#endif
  for (j = 1; j < d; j++) {
    long idx = lp->nbindex[j + 1];
#ifdef DEBUG_CDD_DISABLE
    os << "CDD: j=" << static_cast<int>(j) << "/" << static_cast<int>(d) << " idx=" << static_cast<int>(idx) << " |EXT|=" << EXT.rows() << " / " << EXT.cols() << "\n";
#endif
    V(j - 1) = EXT(idx - 1, 0);
    for (i = 0; i < d - 1; i++) {
      M(i, j - 1) = -EXT(idx - 1, i + 1);
    }
#ifdef DEBUG_CDD_DISABLR
    os << "CDD: After V and M sets\n";
#endif
  }
#ifdef DEBUG_CDD
  os << "CDD: Definition of V and M (step 2)\n";
#endif
  std::optional<MyVector<T>> optA = SolutionMat(M, V);
#ifdef DEBUG_CDD
  os << "CDD: We have optA\n";
#endif
  if (!optA) {
#ifdef DEBUG_CDD
    os << "CDD: LIFT ERROR, Could not find a solution to SolutionMat(M, V)\n";
#endif
    return {};
  }
  //
  // Getting the direct solution and testing it.
  //
  MyVector<T> const &DirectSolution = *optA;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    T eSum = EXT(iRow, 0);
    for (int iCol = 0; iCol < nbCol - 1; iCol++) {
      eSum += DirectSolution(iCol) * EXT(iRow, iCol + 1);
    }
    if (eSum < 0) {
#ifdef DEBUG_CDD
      os << "CDD: LIFT ERROR, Not an interior point at iRow=" << iRow
         << " eSum=" << eSum << "\n";
#endif
      return {};
    }
  }
  T objDirect = eVect(0);
  for (int iCol = 0; iCol < nbCol - 1; iCol++) {
    objDirect += DirectSolution(iCol) * eVect(iCol + 1);
  }
#ifdef DEBUG_CDD
  os << "CDD: We have objDirect\n";
#endif
  //
  // Getting the dual solution
  //
  MyVector<T> eVectRed(nbCol - 1);
  for (int iCol = 0; iCol < nbCol - 1; iCol++) {
    eVectRed(iCol) = eVect(iCol + 1);
  }
#ifdef DEBUG_CDD
  os << "CDD: We have eVectRed\n";
#endif
  MyMatrix<T> M2(nbCol - 1, nbCol - 1);
  for (j = 1; j < d; j++) {
    long idx = lp->nbindex[j + 1];
    for (int iCol = 0; iCol < nbCol - 1; iCol++) {
      M2(j - 1, iCol) = EXT(idx - 1, iCol + 1);
    }
  }
#ifdef DEBUG_CDD
  os << "CDD: We have M2\n";
#endif
  std::optional<MyVector<T>> optB = SolutionMat(M2, eVectRed);
#ifdef DEBUG_CDD
  os << "CDD: We have optB\n";
#endif
  if (!optB) {
#ifdef DEBUG_CDD
    os << "CDD: LIFT ERROR: No solution found for SolutionMat(M2, eVectRed)\n";
#endif
    return {};
  }
  MyVector<T> const &partDualSolution = *optB;
  MyVector<T> DualSolution = ZeroVector<T>(nbRow);
  T objDual = eVect(0);
  for (j = 1; j < d; j++) {
    long idx = lp->nbindex[j + 1];
    T scal = -partDualSolution(j - 1);
    DualSolution(idx - 1) = scal;
    if (scal > 0) {
#ifdef DEBUG_CDD
      os << "CDD: LIFT ERROR, scal=" << scal << "\n";
#endif
      return {};
    }
    objDual += scal * EXT(idx - 1, 0);
  }
#ifdef DEBUG_CDD
  os << "CDD: We have objDual\n";
#endif
  if (objDual != objDirect) {
#ifdef DEBUG_CDD
    os << "CDD: LIFT ERROR, objDual=" << objDual << " objDirect=" << objDirect
       << "\n";
#endif
    return {};
  }
#ifdef DEBUG_CDD
  os << "CDD: An apparently valid solution has been found\n";
#endif
  LpSolution<T> eSol;
  eSol.PrimalDefined = true;
  eSol.DualDefined = true;
  eSol.DualSolution = DualSolution;
  eSol.OptimalValue = objDirect;
  eSol.DirectSolution = DirectSolution;
  return eSol;
}

template <typename T>
LpSolution<T>
CDD_LinearProgramming_exact_V1(MyMatrix<T> const &EXT, MyVector<T> const &eVect,
                               std::ostream &os) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  cdd::dd_ErrorType error = cdd::dd_NoError;
  cdd::dd_matrixdata<T> *M;
  cdd::dd_LPSolverType solver =
      cdd::dd_DualSimplex; /* either DualSimplex or CrissCross */
  cdd::dd_lpdata<T>
      *lp; /* pointer to LP data structure that is not visible by user. */
  M = cdd::MyMatrix_PolyFile2Matrix(EXT);
  M->representation = cdd::dd_Inequality;
  cdd::dd_colrange j;
  cdd::dd_colrange d_input = EXT.cols();
  for (j = 0; j < d_input; j++) {
    M->rowvec[j] = eVect(j);
  }
  lp = cdd::dd_Matrix2LP(M);
  lp->objective = cdd::dd_LPmin;
  size_t maxiter = 0;
  dd_LPSolve(lp, solver, &error, maxiter, os);
  std::optional<LpSolution<T>> optA =
      GetLpSolutionFromLpData(EXT, eVect, lp, os);
  if (optA) {
    LpSolution<T> const &eSolA = *optA;
#ifdef DEBUG_CDD
    if (is_lifting_possible(lp)) {
      std::optional<LpSolution<T>> optB =
        LiftFloatingPointSolution(EXT, eVect, lp, os);
      if (optB) {
        LpSolution<T> const &eSolB = *optB;
        if (eSolB.OptimalValue != eSolA.OptimalValue) {
          std::cerr << "CDD: We should have the same optimal value\n";
          throw TerminalException{1};
        }
        if (eSolA.DualSolution != eSolB.DualSolution) {
          std::cerr << "CDD: DualSolution(A)=" << StringVector(eSolA.DualSolution)
                    << "\n";
          std::cerr << "CDD: DualSolution(B)=" << StringVector(eSolB.DualSolution)
                    << "\n";
          throw TerminalException{1};
        }
        if (eSolA.DirectSolution != eSolB.DirectSolution) {
          std::cerr << "CDD: DirectSolution(A)=" << StringVector(eSolA.DirectSolution)
                    << "\n";
          std::cerr << "CDD_ DirectSolution(B)=" << StringVector(eSolB.DirectSolution)
                    << "\n";
          throw TerminalException{1};
        }
        os << "CDD: DualSolution(A)=" << StringVector(eSolA.DualSolution) << "\n";
        os << "CDD: DualSolution(B)=" << StringVector(eSolB.DualSolution) << "\n";
        os << "CDD: DirectSolution(A)=" << StringVector(eSolA.DirectSolution) << "\n";
        os << "CDD: DirectSolution(B)=" << StringVector(eSolB.DirectSolution) << "\n";
      } else {
        std::cerr << "CDD: We should have been able to lift the solution\n";
        throw TerminalException{1};
      }
    }
#endif
    dd_FreeMatrix(M);
    dd_FreeLPData(lp);
    return eSolA;
  } else {
    throw TerminalException{1};
  }
}

template <typename T, typename Tfloat>
LpSolution<T>
CDD_LinearProgramming_exact_V2(MyMatrix<T> const &EXT, MyVector<T> const &eVect,
                               [[maybe_unused]] std::ostream &os) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  MyMatrix<Tfloat> EXT_float = UniversalMatrixConversion<Tfloat, T>(EXT);
  cdd::dd_ErrorType error = cdd::dd_NoError;
  cdd::dd_matrixdata<Tfloat> *M;
  cdd::dd_LPSolverType solver =
      cdd::dd_DualSimplex; /* either DualSimplex or CrissCross */
  cdd::dd_lpdata<Tfloat>
      *lp; /* pointer to LP data structure that is not visible by user. */
  M = cdd::MyMatrix_PolyFile2Matrix(EXT_float);
  M->representation = cdd::dd_Inequality;
  cdd::dd_colrange j;
  cdd::dd_colrange d_input = EXT.cols();
  for (j = 0; j < d_input; j++) {
    M->rowvec[j] = UniversalScalarConversion<Tfloat, T>(eVect(j));
  }
  lp = cdd::dd_Matrix2LP(M);
  lp->objective = cdd::dd_LPmin;
#ifdef DEBUG_CDD
  os << "CDD: Before dd_LPSolve\n";
#endif
  size_t maxiter = 10 * (EXT.rows() + EXT.cols() + 3);
  dd_LPSolve(lp, solver, &error, maxiter, os);
#ifdef DEBUG_CDD
  os << "CDD: After dd_LPSolve\n";
#endif
  if (lp->LPS == cdd::dd_TooManyIterations) {
#ifdef DEBUG_CDD
    os << "CDD: Error is TooManyIterations, calling CDD_LinearProgramming_exact_V1\n";
#endif
    return CDD_LinearProgramming_exact_V1(EXT, eVect, os);
  }
#ifdef DEBUG_CDD
  os << "CDD: After dd_LPSolve passed \n";
#endif
  if (is_lifting_possible(lp)) {
#ifdef DEBUG_CDD
    os << "CDD: is_lifting_possible = true\n";
#endif
    std::optional<LpSolution<T>> optB =
      LiftFloatingPointSolution<T, Tfloat>(EXT, eVect, lp, os);
#ifdef DEBUG_CDD
    os << "CDD: We have optB\n";
#endif
    if (optB) {
#ifdef DEBUG_CDD
      os << "CDD: The lifing of floating point solution went nicely\n";
#endif
      return *optB;
    }
  }
#ifdef DEBUG_CDD
  os << "CDD: The lifting scheme failed, now using the direct approach\n";
#endif
  return CDD_LinearProgramming_exact_V1(EXT, eVect, os);
}

template <typename T>
LpSolution<T> CDD_LinearProgramming(MyMatrix<T> const &EXT, MyVector<T> const &eVect, std::ostream &os) {
  if (EXT.cols() < 4) {
    return CDD_LinearProgramming_exact_V1(EXT, eVect, os);
  }
  return CDD_LinearProgramming_exact_V2<T, double>(EXT, eVect, os);
}

template <typename T>
LpSolution<T> CDD_LinearProgramming_BugSearch(MyMatrix<T> const &TheEXT,
                                              MyVector<T> const &eVect,
                                              std::ostream &os) {
  LpSolution<T> eSol1 = CDD_LinearProgramming(TheEXT, eVect, os);
  LpSolution<T> eSol2 = CDD_LinearProgramming_External(TheEXT, eVect, os);
  if (eSol1.PrimalDefined != eSol2.PrimalDefined ||
      eSol1.DualDefined != eSol2.DualDefined) {
    WriteInputFileCdd("bugSearch.ine", TheEXT, eVect);
    std::cerr << "We find the bug we were after\n";
    throw TerminalException{1};
  }
  return eSol1;
}

// clang-format off
#endif  // SRC_POLY_POLY_CDDLIB_H_
// clang-format on

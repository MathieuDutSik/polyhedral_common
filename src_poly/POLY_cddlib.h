#ifndef CDD_TEMP_H_
#define CDD_TEMP_H_
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
   This C++ version is adapted by Mathieu Dutour Sikiric
*/

#include "MAT_Matrix.h"
#include "Boost_bitset.h"
#include<cassert>

#define DEBUG_CDD

namespace cdd {
  typedef enum {
    dd_CrissCross, dd_DualSimplex
  } dd_LPSolverType;
  const dd_LPSolverType dd_choiceLPSolverDefault=dd_DualSimplex;  /* Default LP solver Algorithm */
  const dd_LPSolverType dd_choiceRedcheckAlgorithm=dd_DualSimplex;  /* Redundancy Checking Algorithm */
  const bool dd_choiceLexicoPivotQ=true;    /* whether to use the lexicographic pivot */
  typedef unsigned char set_card_lut_t;

  template<typename T>
  inline void dd_set(T& a, T b)
  {
    a=b;
  }

  template<typename T>
  inline void dd_set_si(T& a, int b)
  {
    a=b;
  }

template<typename T>
inline void dd_set_si2(T& a, int b)
{
  a=b;
}

template<typename T>
inline void dd_add(T &a, T b, T c)
{
  a=b + c;
}

template<typename T>
inline void dd_sub(T &a, T b, T c)
{
  a=b - c;
}

template<typename T>
inline void dd_mul(T &a, T b, T c)
{
  a=b * c;
}

template<typename T>
inline void dd_div(T &a, T b, T c)
{
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  a=b / c;
}

template<typename T>
inline void dd_neg(T &a, T b)
{
  a=-b;
}

template<typename T>
inline void dd_inv(T &a, T b)
{
  a=1/b;
}

template<typename T>
inline int dd_pos(T a, T b)
{
  if (a>b)
    return 1;
  if (a<b)
    return -1;
  return 0;
}


template<typename T>
void dd_WriteT(std::ostream & os, T* a, int d)
{
  os << "dd_Write T a=";
  for (int i=0; i<d; i++)
    os << " " << a[i];
  os << "\n";
}

template<typename T>
inline bool dd_Nonnegative(T val)
{
  return val >= 0;
}

template<typename T>
inline bool dd_Nonpositive(T val)
{
  return val <= 0;
}

template<typename T>
inline bool dd_Positive(T val)
{
  return val > 0;
}

template<typename T>
inline bool dd_Negative(T val)
{
  return val < 0;
}

template<typename T>
inline bool dd_EqualToZero(T val)
{
  return val == 0;
}

template<typename T>
inline bool dd_Nonzero(T val)
{
  return val != 0;
}

template<typename T>
inline bool dd_Larger(T val1,T val2)
{
  return val1 > val2;
}

template<typename T>
inline bool dd_LargerFrac(T val1, T q1, T val2, T q2)
{
#ifdef DEBUG_CDD
  assert(q1>0);
  assert(q2>0);
#endif
  return val1 * q2 > val2 * q1;
}

template<typename T>
inline bool dd_Smaller(T val1,T val2)
{
  return val1 < val2;
}

template<typename T>
inline bool dd_SmallerFrac(T val1, T q1, T val2, T q2)
// We want to have val1 / q1 < val2 / q2 which is equivalent to val1 * q2 < val2 * q1
{
#ifdef DEBUG_CDD
  assert(q1>0);
  assert(q2>0);
#endif
  return val1 * q2 < val2 * q1;
}

template<typename T>
inline bool dd_Equal(T val1,T val2)
{
  return val1 == val2;
}

template<typename T>
inline bool dd_EqualFrac(T val1,T q1, T val2, T q2)
{
#ifdef DEBUG_CDD
  assert(q1>0);
  assert(q2>0);
#endif
  return val1 * q2 < val2 * q1;
}

template<typename T>
inline void dd_abs(T &absval, T val)
{
  if (val < 0) dd_neg(absval,val);
  else dd_set(absval,val);
}



typedef long dd_rowrange;
typedef long dd_colrange;
typedef long dd_bigrange;

typedef long *dd_rowindex;
typedef int *dd_rowflag;
typedef long *dd_colindex;

// API for the set systems.
typedef unsigned long *set_type;   /* set type definition */

inline void set_free(set_type set)
/* Free the space created with the set pointer set*/
{
    delete [] set;
}


#define SETBITS (sizeof(long) * 8)

#define LUTBLOCKS(set) (((set[0]-1)/SETBITS+1)*(sizeof(long)/sizeof(set_card_lut_t)))

static unsigned char set_card_lut[]={
0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};
/* End of Definitions for optimized set_card */

unsigned long set_blocks(long len)
{
	long blocks=1L;

	if (len>0) blocks=((long)len-1)/SETBITS+2;
	return blocks;
}

void set_initialize(set_type *setp, long length)
/* Make a set with a given bit lengths  */
{
    long i,forlim1,len;

    if (length<=0)
      len=1;
    else
      len=length;
     /* if negative length is requested, it generates the shortest length */

    forlim1=set_blocks(len);
    using ulong=unsigned long;
    *setp=new ulong[forlim1];
    (*setp)[0]=ulong(len);  /* size of the ground set */
    for (i=1; i<forlim1; i++)
      (*setp)[i]=0U;
}

void set_emptyset(set_type set)
/* Set set to be the emptyset  */
{
	long i,forlim;

	forlim=set_blocks(set[0])-1;
	for (i=1; i<=forlim; i++)
		set[i]=0U;
}

void set_copy(set_type setcopy,set_type set)
/* Copy the set set[] to setcopy[] with setcopy[] length */
{
	long i,forlim;

	forlim=set_blocks(setcopy[0])-1;
	for (i=1; i<=forlim; i++)
		setcopy[i]=set[i];
}

void set_addelem(set_type set, long elem)
/* add elem only if it is within the set[] range */
{
  long i,j;
  unsigned long change;
  unsigned long one=1U;

  if (elem<=(long)set[0])
    {
      i=(elem-1)/SETBITS+1;
      j=(elem-1)%SETBITS;
      change= one << j;  /* put 1 in jth position */
      set[i]=set[i] | change;
    }
}

void set_delelem(set_type set, long elem)
/* delete elem only if it is within the set[] range */
{
	long  i,j;
	unsigned long change;
	unsigned long one=1U;

	if (elem<=(long)set[0])
	{
		i=(elem-1)/SETBITS+1;
		j=(elem-1)%SETBITS;
		change=one << j;  /* put 1 in jth position */
		set[i]=(set[i] | change) ^ change;
	}
}

void set_int(set_type set,set_type set1,set_type set2)
/* Set intersection, assuming set1 and set2 have the same length as set */
{
  long  i,forlim;

  forlim=set_blocks(set[0])-1;
  for (i=1;i<=forlim;i++)
    set[i]=(set1[i] & set2[i]);
}

void set_uni(set_type set,set_type set1,set_type set2)
/* Set union,assuming set1 and set2 have the same length as set */
{
	long  i,forlim;

	forlim=set_blocks(set[0])-1;
	for (i=1;i<=forlim;i++)
		set[i]=set1[i] | set2[i];
}

void set_diff(set_type set,set_type set1,set_type set2)
/* Set difference se1/set2, assuming set1 and set2 have the same length as set */
{
	long  i,forlim;

	forlim=set_blocks(set[0])-1;
	for (i=1;i<=forlim;i++)
		set[i]=set1[i] & (~set2[i]);
}

void set_compl(set_type set,set_type set1)
/* set[] will be set to the complement of set1[] */
{
	long  i,j,l,forlim;
	unsigned long change;
	unsigned long one=1U;

	forlim=set_blocks(set[0])-1;
	for (i=1;i<=forlim;i++)
		set[i]= ~set1[i];
	l=(set[0]-1)%SETBITS; /* the position of the last elem in the last block */
    	for (j=l+1; j<=(long)SETBITS-1; j++){
		change=one << j;
		set[forlim]=(set[forlim] | change) ^ change;
    	}
}

int set_subset(set_type set1,set_type set2)
/* Set containment check, set1 <= set2 */
{
	int  yes=1;
	long i,forlim;

	forlim=set_blocks(set2[0])-1;
	for (i=1;i<=forlim && yes;i++)
		if ((set1[i] | set2[i])!=set2[i])
			yes=0;
	return yes;
}

int set_member(long elem, set_type set)
/* Set membership check, elem in set */
{
	int  yes=0;
	long  i,j;
	unsigned long testset;
	unsigned long one=1U;

	if (elem<=(long)set[0])
	{
		i=(elem-1)/SETBITS+1;
		j=(elem-1)%SETBITS;
		testset=set[i] | (one<<j);   /* add elem to set[i] */
		if (testset==set[i])
			yes=1;
	}
	return yes;
}

/*set cardinality, modified by David Bremner bremner@cs.mcgill.ca
   to optimize for speed.*/
long set_card(set_type set)
{
  unsigned long block;
  long car=0;
  set_card_lut_t *p;

  p=(set_card_lut_t *)&set[1];
  for (block=0; block< LUTBLOCKS(set);block++) {
    car+=set_card_lut[p[block]];
  }
  return car;
}


typedef set_type dd_rowset;
typedef set_type dd_colset;
typedef set_type *dd_SetVector;
typedef set_type *dd_Aincidence;
typedef set_type rowset;  /* set_type defined in setoper.h */
typedef set_type colset;


// More functionalities



template<typename T>
inline void dd_InnerProduct(T& prod, dd_colrange d, T* v1, T* v2)
{
  T temp;
  dd_colrange j;
  dd_set_si(prod, 0);
  for (j = 0; j < d; j++){
    dd_mul(temp,v1[j],v2[j]);
    dd_add(prod,prod,temp);
  }
}


template<typename T>
struct dd_raydata {
  T *Ray;
  dd_rowset ZeroSet;
  dd_rowrange FirstInfeasIndex;  /* the first inequality the ray violates */
  bool feasible;  /* flag to store the feasibility */
  T ARay;   /* temporary area to store some row of A*Ray */
  dd_raydata<T> *Next;
};


template<typename T>
struct dd_adjacencydata {
  dd_raydata<T> *Ray1;
  dd_raydata<T> *Ray2;
  dd_adjacencydata<T> *Next;
};

typedef enum {
  dd_Combinatorial, dd_Algebraic
} dd_AdjacencyTestType;

typedef enum {
  dd_MaxIndex, dd_MinIndex, dd_MinCutoff, dd_MaxCutoff, dd_MixCutoff,
   dd_LexMin, dd_LexMax, dd_RandomRow
} dd_RowOrderType;

typedef enum {
  dd_Unspecified=0, dd_Inequality, dd_Generator
} dd_RepresentationType;

typedef enum {
  dd_IneToGen, dd_GenToIne, dd_LPMax, dd_LPMin, dd_InteriorFind
} dd_ConversionType;

typedef enum {
  dd_IncOff=0, dd_IncCardinality, dd_IncSet
} dd_IncidenceOutputType;

typedef enum {
  dd_AdjOff=0, dd_AdjacencyList,  dd_AdjacencyDegree
} dd_AdjacencyOutputType;

typedef enum {
  dd_Auto, dd_SemiAuto, dd_Manual
} dd_FileInputModeType;
   /* Auto if a input filename is specified by command arguments */

typedef enum {
  dd_DimensionTooLarge, dd_ImproperInputFormat,
  dd_NegativeMatrixSize, dd_EmptyVrepresentation, dd_EmptyHrepresentation, dd_EmptyRepresentation,
  dd_IFileNotFound, dd_OFileNotOpen, dd_NoLPObjective, dd_NoRealNumberSupport,
  dd_NotAvailForH, dd_NotAvailForV, dd_CannotHandleLinearity,
  dd_RowIndexOutOfRange, dd_ColIndexOutOfRange,
  dd_LPCycling, dd_NumericallyInconsistent,
  dd_NoError
} dd_ErrorType;

typedef enum {
  dd_InProgress, dd_AllFound, dd_RegionEmpty
} dd_CompStatusType;

/* --- LP types ---- */

typedef enum {
  dd_LPnone=0, dd_LPmax, dd_LPmin
} dd_LPObjectiveType;


typedef enum {
  dd_LPSundecided, dd_Optimal, dd_Inconsistent, dd_DualInconsistent,
  dd_StrucInconsistent, dd_StrucDualInconsistent,
  dd_Unbounded, dd_DualUnbounded
} dd_LPStatusType;


template<typename T>
struct dd_lpsolution {
  //  dd_DataFileType filename;
  dd_LPObjectiveType objective;
  dd_LPSolverType solver;
  dd_rowrange m;
  dd_colrange d;

  dd_LPStatusType LPS;  /* the current solution status */
  T optvalue;  /* optimal value */
  T* sol;   /* primal solution */
  T* dsol;  /* dual solution */
  dd_colindex nbindex;  /* current basis represented by nonbasic indices */
  dd_rowrange re;  /* row index as a certificate in the case of inconsistency */
  dd_colrange se;  /* col index as a certificate in the case of dual inconsistency */
  long pivots[5];
   /* pivots[0]=setup (to find a basis), pivots[1]=PhaseI or Criss-Cross,
      pivots[2]=Phase II, pivots[3]=Anticycling, pivots[4]=GMP postopt. */
  long total_pivots;
};


template<typename T>
struct dd_lpdata {
  //  dd_DataFileType filename;
  dd_LPObjectiveType objective;
  dd_LPSolverType solver;
  bool Homogeneous;
     /* The first column except for the obj row is all zeros. */
  dd_rowrange m;
  dd_colrange d;
  T** A;
  T** B;
  dd_rowrange objrow;
  dd_colrange rhscol;
  dd_rowrange eqnumber;  /* the number of equalities */
  dd_rowset equalityset;

  bool redcheck_extensive;  /* Apply the extensive redundancy check. */
  dd_rowrange ired; /* the row index for the redundancy checking */
  dd_rowset redset_extra;  /* a set of rows that are newly recognized redundan by the extensive search. */
  dd_rowset redset_accum;  /* the accumulated set of rows that are recognized redundant */
  dd_rowset posset_extra;  /* a set of rows that are recognized non-linearity */

  bool lexicopivot;  /* flag to use the lexicogrphic pivot rule (symbolic perturbation). */

  dd_LPStatusType LPS;  /* the current solution status */
  dd_rowrange m_alloc; /* the allocated row size of matrix A */
  dd_colrange d_alloc; /* the allocated col size of matrix A */
  T optvalue;  /* optimal value */
  T* sol;   /* primal solution */
  T* dsol;  /* dual solution */
  dd_colindex nbindex;  /* current basis represented by nonbasic indices */
  dd_rowrange re;  /* row index as a certificate in the case of inconsistency */
  dd_colrange se;  /* col index as a certificate in the case of dual inconsistency */
  long pivots[5];
   /* pivots[0]=setup (to find a basis), pivots[1]=PhaseI or Criss-Cross,
      pivots[2]=Phase II, pivots[3]=Anticycling, pivots[4]=GMP postopt. */
  long total_pivots;
  int use_given_basis;  /* switch to indicate the use of the given basis */
  dd_colindex given_nbindex;  /* given basis represented by nonbasic indices */
};


/*----  end of LP Types ----- */

template<typename T>
struct data_temp_simplex {
  T* rcost;
  dd_colset tieset, stieset;
  T* Rtemp;
  dd_colindex nbindex_ref;
};

template<typename T>
data_temp_simplex<T>* allocate_data_simplex(dd_colrange d_size)
{
  //  std::cout << "allocate_data_simplex with d_size=" << d_size << "\n";
  data_temp_simplex<T>* data = new data_temp_simplex<T>;
  data->rcost = new T[d_size];
  set_initialize(&data->tieset,d_size);
  set_initialize(&data->stieset,d_size);
  data->Rtemp = new T[d_size];
  data->nbindex_ref = new long[d_size+1];
  return data;
}

template<typename T>
void free_data_simplex(data_temp_simplex<T>* data)
{
  delete [] data->rcost;
  set_free(data->tieset);
  set_free(data->stieset);
  delete [] data->Rtemp;
  delete [] data->nbindex_ref;
}






/* Containing data sets */


template<typename T>
void dd_FreeAmatrix(dd_rowrange m, T** A)
{
  dd_rowrange i;

  if (A!=nullptr) {
    for (i = 0; i < m; i++) {
      delete [] A[i];
    }
    delete [] A;
  }
}

template<typename T>
void dd_FreeBmatrix(dd_colrange d,T** B)
{
  dd_colrange j;

  if (B!=nullptr) {
    for (j = 0; j < d; j++) {
      delete [] B[j];
    }
    delete [] B;
  }
}



template<typename T>
struct dd_matrixdata {
  dd_rowrange rowsize;
  dd_rowset linset;
    /*  a subset of rows of linearity (ie, generators of
        linearity space for V-representation, and equations
        for H-representation. */
  dd_colrange colsize;
  dd_RepresentationType representation;
  T** matrix;
  dd_LPObjectiveType objective;
  T* rowvec;
};

template<typename T>
void dd_FreeArow(T* a)
{
  delete [] a;
}

template<typename T>
void dd_FreeMatrix(dd_matrixdata<T> *M)
{
  dd_rowrange m1;

  if (M!=nullptr) {
    if (M->rowsize<=0) m1=1; else m1=M->rowsize;
    if (M!=nullptr) {
      dd_FreeAmatrix(m1, M->matrix);
      dd_FreeArow(M->rowvec);
      set_free(M->linset);
      delete M;
    }
  }
}



template<typename T>
void dd_FreeShallowMatrix(dd_matrixdata<T> *M)
{
  dd_rowrange m1;

  if (M!=nullptr) {
    if (M!=nullptr) {
      delete [] M->matrix;
      dd_FreeArow(M->rowvec);
      set_free(M->linset);
      delete M;
    }
  }
}

typedef struct dd_setfamily *dd_SetFamilyPtr;
typedef struct dd_setfamily {
  dd_bigrange famsize;
  dd_bigrange setsize;
  dd_SetVector set;
} dd_SetFamilyType;



void dd_FreeSetFamily(dd_SetFamilyPtr F)
{
  dd_bigrange i,f1;

  if (F!=nullptr){
    if (F->famsize<=0) f1=1; else f1=F->famsize;
      /* the smallest created size is one */
    for (i=0; i<f1; i++) {
      set_free(F->set[i]);
    }
    delete [] F->set;
    delete [] F;
  }
}




typedef struct dd_nodedata *dd_NodePtr;
typedef struct dd_nodedata {dd_bigrange key; dd_NodePtr next;} dd_NodeType;

typedef struct dd_graphdata *dd_GraphPtr;
typedef struct dd_graphdata {
  dd_bigrange vsize;
  dd_NodePtr *adjlist;  /* should be initialized to have vsize components */
} dd_GraphType;


//typedef struct dd_polyhedradata *dd_PolyhedraPtr;
//typedef struct dd_conedata *dd_ConePtr;
template<typename T>
struct dd_polyhedradata;


template<typename T>
struct dd_conedata {
  dd_RepresentationType representation;
  dd_rowrange m;
  dd_colrange d;
  T** A;
  dd_polyhedradata<T> *parent;  /* pointing to the original polyhedra data */
  dd_rowrange m_alloc; /* allocated row size of matrix A */
  dd_colrange d_alloc; /* allocated col size of matrix A */

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
  dd_adjacencydata<T> **Edges;  /* adjacency relation storage for iteration k */
  unsigned int rseed;  /* random seed for random row permutation */

  bool ColReduced;  /* flag to indicate that a column basis is computed and reduced */
  dd_bigrange LinearityDim;
           /*  the dimension of the linearity space (when input is H), and
               the size of a minimal system of equations to determine the space (when V). */
  dd_colrange d_orig;  /* the size d of the original matrix A */
  dd_colindex newcol;  /* the size d of the original matrix A */

  dd_colindex InitialRayIndex;   /* InitialRayIndex[s] (s>=1) stores the corr. row index */
  dd_rowindex OrderVector;
  bool RecomputeRowOrder;
  bool PreOrderedRun;
  dd_rowset GroundSet, EqualitySet, NonequalitySet,
       AddedHalfspaces, WeaklyAddedHalfspaces, InitialHalfspaces;
  long RayCount, FeasibleRayCount, WeaklyFeasibleRayCount,
       TotalRayCount, ZeroRayCount;
  long EdgeCount, TotalEdgeCount;
  long count_int,count_int_good,count_int_bad; /* no. of intersection operations */

  T** B;
  T** Bsave;  /* a copy of the dual basis inverse used to reduce the matrix A */

  dd_ErrorType Error;
  dd_CompStatusType CompStatus;  /* Computation Status */
};


template<typename T>
struct dd_polyhedradata {
  dd_RepresentationType representation;  /* given representation */
  bool homogeneous;
  dd_colrange d;
  dd_rowrange m;
  T** A;   /* Inequality System:  m times d matrix */
  dd_conedata<T> *child;  /* pointing to the homogenized cone data */
  dd_rowrange m_alloc; /* allocated row size of matrix A */
  dd_colrange d_alloc; /* allocated col size of matrix A */
  T* c;           /* cost vector */

  dd_rowflag EqualityIndex;
    /* ith component is 1 if it is equality, -1 if it is strict inequality, 0 otherwise. */

  bool IsEmpty;  /* This is to tell whether the set is empty or not */

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
  dd_colrange ldim;   /* linearity dimension */
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


dd_SetFamilyPtr dd_CreateSetFamily(dd_bigrange fsize, dd_bigrange ssize)
{
  dd_SetFamilyPtr F;
  dd_bigrange i,f0,f1,s0,s1;

  if (fsize<=0) {
    f0=0; f1=1;
    /* if fsize<=0, the fsize is set to zero and the created size is one */
  } else {
    f0=fsize; f1=fsize;
  }
  if (ssize<=0) {
    s0=0; s1=1;
    /* if ssize<=0, the ssize is set to zero and the created size is one */
  } else {
    s0=ssize; s1=ssize;
  }

  F = new dd_SetFamilyType;
  F->set = new set_type[f1];
  for (i=0; i<f1; i++)
    set_initialize(&(F->set[i]), s1);
  F->famsize=f0;
  F->setsize=s0;
  return F;
}

template<typename T>
bool dd_AppendMatrix2Poly(dd_polyhedradata<T> **poly, dd_matrixdata<T> *M)
{
  bool success=false;
  dd_matrixdata<T> *Mpoly,Mnew=nullptr;
  dd_ErrorType err;

  if ((*poly)!=nullptr && (*poly)->m >=0 && (*poly)->d>=0 &&
      (*poly)->d==M->colsize && M->rowsize>0){
    Mpoly=dd_CopyInput(*poly);
    Mnew=dd_AppendMatrix(Mpoly, M);
    dd_FreePolyhedra(*poly);
    *poly=dd_DDMatrix2Poly(Mnew,&err);
    dd_FreeMatrix(Mpoly);
    dd_FreeMatrix(Mnew);
    if (err==dd_NoError) success=true;
  }
  return success;
}

// Basic Initialization procedures

template<typename T>
void dd_InitializeAmatrix(dd_rowrange m,dd_colrange d,T** *A)
{
  dd_rowrange i;

  if (m>0) (*A)=new T*[m];
  for (i = 0; i < m; i++)
    (*A)[i] =new T[d];
}

template<typename T>
void dd_InitializeArow(dd_colrange d, T* *a)
{
  if (d>0) *a=new T[d];
}

template<typename T>
void dd_InitializeBmatrix(dd_colrange d,T** *B)
{
  dd_colrange j;

  (*B)=new T*[d];
  for (j = 0; j < d; j++) {
    (*B)[j]=new T[d];
  }
}

// More sophisticated creation functionalities.

template<typename T>
dd_matrixdata<T> *dd_CreateMatrix(dd_rowrange m_size,dd_colrange d_size)
{
  dd_matrixdata<T> *M;
  dd_rowrange m0,m1;
  dd_colrange d0,d1;

  if (m_size<=0){
    m0=0; m1=1;
    /* if m_size <=0, the number of rows is set to zero, the actual size is 1 */
  } else {
    m0=m_size; m1=m_size;
  }
  if (d_size<=0){
    d0=0; d1=1;
    /* if d_size <=0, the number of cols is set to zero, the actual size is 1 */
  } else {
    d0=d_size; d1=d_size;
  }
  M=new dd_matrixdata<T>;
  dd_InitializeAmatrix(m1,d1,&(M->matrix));
  dd_InitializeArow(d1,&(M->rowvec));
  M->rowsize=m0;
  set_initialize(&(M->linset), m1);
  M->colsize=d0;
  M->objective=dd_LPnone;
  M->representation=dd_Unspecified;
  return M;
}


/* This is a kind of shallow copy. The matrix entries are not assigned but supposed
   to be obtained via pointer assignment from another matrix.
   However, we may lose some speed due to alignment stuff, so it is unsure
   if this makes sense at all.
*/
template<typename T>
dd_matrixdata<T> *dd_CreateShallowMatrix(dd_rowrange m_size,dd_colrange d_size)
{
  dd_matrixdata<T> *M;
  dd_rowrange m0,m1;
  dd_colrange d0,d1;

  if (m_size<=0){
    m0=0; m1=1;
    /* if m_size <=0, the number of rows is set to zero, the actual size is 1 */
  } else {
    m0=m_size; m1=m_size;
  }
  if (d_size<=0){
    d0=0; d1=1;
    /* if d_size <=0, the number of cols is set to zero, the actual size is 1 */
  } else {
    d0=d_size; d1=d_size;
  }
  M=new dd_matrixdata<T>;
  M->matrix = new T*[m1];
  dd_InitializeArow(d1,&(M->rowvec));
  M->rowsize=m0;
  set_initialize(&(M->linset), m1);
  M->colsize=d0;
  M->objective=dd_LPnone;
  M->representation=dd_Unspecified;
  return M;
}




template<typename T>
dd_polyhedradata<T> *dd_CreatePolyhedraData(dd_rowrange m, dd_colrange d)
{
  dd_rowrange i;
  dd_polyhedradata<T> *poly=nullptr;

  poly=new dd_polyhedradata<T>;
  poly->child       =nullptr; /* this links the homogenized cone data */
  poly->m           =m;
  poly->d           =d;
  poly->n           =-1;  /* the size of output is not known */
  poly->m_alloc     =m+2; /* the allocated row size of matrix A */
  poly->d_alloc     =d;   /* the allocated col size of matrix A */
  poly->ldim		=0;   /* initialize the linearity dimension */
  dd_InitializeAmatrix(poly->m_alloc,poly->d_alloc,&(poly->A));
  dd_InitializeArow(d,&(poly->c));           /* cost vector */
  poly->representation       =dd_Inequality;
  poly->homogeneous =false;

  poly->EqualityIndex = new int[m+2];
    /* size increased to m+2 in 092b because it is used by the child cone,
       This is a bug fix suggested by Thao Dang. */
    /* ith component is 1 if it is equality, -1 if it is strict inequality, 0 otherwise. */
  for (i = 0; i <= m+1; i++) poly->EqualityIndex[i]=0;

  poly->IsEmpty                 = -1;  /* initially set to -1, neither TRUE nor FALSE, meaning unknown */

  poly->NondegAssumed           = false;
  poly->InitBasisAtBottom       = false;
  poly->RestrictedEnumeration   = false;
  poly->RelaxedEnumeration      = false;

  poly->AincGenerated=false;  /* Ainc is a set array to store the input incidence. */

  return poly;
}

template<typename T>
bool dd_InitializeConeData(dd_rowrange m, dd_colrange d, dd_conedata<T> **cone)
{
  bool success=true;
  dd_colrange j;

  (*cone)=new dd_conedata<T>;

/* INPUT: A given representation of a cone: inequality */
  (*cone)->m=m;
  (*cone)->d=d;
  (*cone)->m_alloc=m+2; /* allocated row size of matrix A */
  (*cone)->d_alloc=d;   /* allocated col size of matrix A, B and Bsave */
  (*cone)->parent=nullptr;

/* CONTROL: variables to control computation */
  (*cone)->Iteration=0;

  (*cone)->HalfspaceOrder=dd_LexMin;

  (*cone)->ArtificialRay=nullptr;
  (*cone)->FirstRay=nullptr;
  (*cone)->LastRay=nullptr; /* The second description: Generator */
  (*cone)->PosHead=nullptr;
  (*cone)->ZeroHead=nullptr;
  (*cone)->NegHead=nullptr;
  (*cone)->PosLast=nullptr;
  (*cone)->ZeroLast=nullptr;
  (*cone)->NegLast=nullptr;
  (*cone)->RecomputeRowOrder  = true;
  (*cone)->PreOrderedRun      = false;
  set_initialize(&((*cone)->GroundSet),(*cone)->m_alloc);
  set_initialize(&((*cone)->EqualitySet),(*cone)->m_alloc);
  set_initialize(&((*cone)->NonequalitySet),(*cone)->m_alloc);
  set_initialize(&((*cone)->AddedHalfspaces),(*cone)->m_alloc);
  set_initialize(&((*cone)->WeaklyAddedHalfspaces),(*cone)->m_alloc);
  set_initialize(&((*cone)->InitialHalfspaces),(*cone)->m_alloc);
  (*cone)->RayCount=0;
  (*cone)->FeasibleRayCount=0;
  (*cone)->WeaklyFeasibleRayCount=0;
  (*cone)->TotalRayCount=0;
  (*cone)->ZeroRayCount=0;
  (*cone)->EdgeCount=0;
  (*cone)->TotalEdgeCount=0;
  (*cone)->count_int=0;
  (*cone)->count_int_good=0;
  (*cone)->count_int_bad=0;
  (*cone)->rseed=1;  /* random seed for random row permutation */

  dd_InitializeBmatrix((*cone)->d_alloc, &((*cone)->B));
  dd_InitializeBmatrix((*cone)->d_alloc, &((*cone)->Bsave));
  dd_InitializeAmatrix((*cone)->m_alloc,(*cone)->d_alloc,&((*cone)->A));

  (*cone)->Edges=new dd_adjacencydata<T>*[(*cone)->m_alloc];
  for (j=0; j<(*cone)->m_alloc; j++)
    (*cone)->Edges[j]=nullptr;
  (*cone)->InitialRayIndex = new long[d+1];
  for (int i=0; i<=d; i++)
    (*cone)->InitialRayIndex[i] = 0;
  (*cone)->OrderVector = new long[(*cone)->m_alloc+1];
  for (int i=0; i<=(*cone)->m_alloc; i++)
    (*cone)->OrderVector[i] = 0;

  (*cone)->newcol = new long[((*cone)->d)+1];
  for (int i=0; i<=(*cone)->d; i++)
    (*cone)->newcol[i] = 0;
  for (j=0; j<=(*cone)->d; j++) (*cone)->newcol[j]=j;  /* identity map, initially */
  (*cone)->LinearityDim = -2; /* -2 if it is not computed */
  (*cone)->ColReduced   = false;
  (*cone)->d_orig = d;

/* STATES: variables to represent current state. */
/*(*cone)->Error;
  (*cone)->CompStatus;
  (*cone)->starttime;
  (*cone)->endtime;
*/

  return success;
}

template<typename T>
dd_conedata<T> *dd_ConeDataLoad(dd_polyhedradata<T> *poly)
{
  dd_conedata<T> *cone=nullptr;
  dd_colrange d,j;
  dd_rowrange m,i;

  m=poly->m;
  d=poly->d;
  if (!(poly->homogeneous) && poly->representation==dd_Inequality){
    m=poly->m+1;
  }
  poly->m1=m;

  dd_InitializeConeData(m, d, &cone);
  cone->representation=poly->representation;

/* Points to the original polyhedra data, and reversely */
  cone->parent=poly;
  poly->child=cone;

  for (i=1; i<=poly->m; i++)
    for (j=1; j<=cone->d; j++)
      dd_set(cone->A[i-1][j-1],poly->A[i-1][j-1]);

  if (poly->representation==dd_Inequality && !(poly->homogeneous)){
    cone->A[m-1][0]=1;
    for (j=2; j<=d; j++) cone->A[m-1][j-1]=0;
  }

  return cone;
}

template<typename T>
dd_lpsolution<T> *dd_CopyLPSolution(dd_lpdata<T> *lp)
{
  dd_lpsolution<T> *lps;
  dd_colrange j;

  lps=new dd_lpsolution<T>;
  //  for (i=1; i<=globals::dd_filenamelen; i++) lps->filename[i-1]=lp->filename[i-1];
  lps->objective=lp->objective;
  lps->solver=lp->solver;
  lps->m=lp->m;
  lps->d=lp->d;

  lps->LPS=lp->LPS;  /* the current solution status */

  dd_set(lps->optvalue,lp->optvalue);  /* optimal value */
  dd_InitializeArow(lp->d+1,&(lps->sol));
  dd_InitializeArow(lp->d+1,&(lps->dsol));
  lps->nbindex = new long[(lp->d)+1];
  for (j=0; j<=lp->d; j++){
    dd_set(lps->sol[j],lp->sol[j]);
    dd_set(lps->dsol[j],lp->dsol[j]);
    lps->nbindex[j]=lp->nbindex[j];
  }
  lps->pivots[0]=lp->pivots[0];
  lps->pivots[1]=lp->pivots[1];
  lps->pivots[2]=lp->pivots[2];
  lps->pivots[3]=lp->pivots[3];
  lps->pivots[4]=lp->pivots[4];
  lps->total_pivots=lp->total_pivots;

  return lps;
}



template<typename T>
dd_lpdata<T> *dd_CreateLPData(dd_rowrange m,dd_colrange d)
{
  dd_lpdata<T> *lp;
  lp=new dd_lpdata<T>;
  lp->solver=dd_choiceLPSolverDefault;  /* set the default lp solver */
  lp->d=d;
  lp->m=m;
  lp->objrow=m;
  lp->rhscol=1L;
  lp->objective=dd_LPnone;
  lp->LPS=dd_LPSundecided;
  lp->eqnumber=0;  /* the number of equalities */

  lp->nbindex = new long[d+1];
  for (int i=0; i<=d; i++)
    lp->nbindex[i] = 0;
  lp->given_nbindex = new long[d+1];
  for (int i=0; i<=d; i++)
    lp->given_nbindex[i] = 0;
  set_initialize(&(lp->equalityset),m);
    /* i must be in the set iff i-th row is equality . */

  lp->redcheck_extensive=false; /* this is on only for RedundantExtensive */
  lp->ired=0;
  set_initialize(&(lp->redset_extra),m);
    /* i is in the set if i-th row is newly recognized redundant (during the checking the row ired). */
  set_initialize(&(lp->redset_accum),m);
    /* i is in the set if i-th row is recognized redundant (during the checking the row ired). */
  set_initialize(&(lp->posset_extra),m);
    /* i is in the set if i-th row is recognized non-linearity (during the course of computation). */
  lp->lexicopivot=dd_choiceLexicoPivotQ;  /* dd_choice... is set in dd_set_global_constants() */

  lp->m_alloc=lp->m+2;
  lp->d_alloc=lp->d+2;
  dd_InitializeBmatrix(lp->d_alloc,&(lp->B));
  dd_InitializeAmatrix(lp->m_alloc,lp->d_alloc,&(lp->A));
  dd_InitializeArow(lp->d_alloc,&(lp->sol));
  dd_InitializeArow(lp->d_alloc,&(lp->dsol));
  return lp;
}



template<typename T>
dd_lpdata<T> *dd_CreateLPData_from_M(dd_matrixdata<T> *M)
{
  dd_rowrange m, linc;
  dd_colrange d;
  linc=set_card(M->linset);
  m=M->rowsize+1+linc;
  if (M->representation==dd_Generator){
    d=(M->colsize)+1;
  } else {
    d=M->colsize;
  }
  return dd_CreateLPData<T>(m, d);
}


















template<typename T>
void dd_CopyArow(T *acopy, T *a, dd_colrange d)
{
  dd_colrange j;

  for (j = 0; j < d; j++) {
    dd_set(acopy[j],a[j]);
  }
}


template<typename T>
void dd_CopyAmatrix(T **Acopy, T **A, dd_rowrange m, dd_colrange d)
{
  dd_rowrange i;

  for (i = 0; i< m; i++) {
    dd_CopyArow(Acopy[i],A[i],d);
  }
}

template<typename T>
dd_matrixdata<T> *dd_MatrixCopy(dd_matrixdata<T> *M)
{
  dd_matrixdata<T> *Mcopy=nullptr;
  dd_rowrange m;
  dd_colrange d;

  m = M->rowsize;
  d = M->colsize;
  if (m >=0 && d >=0){
    Mcopy=dd_CreateMatrix<T>(m, d);
    dd_CopyAmatrix(Mcopy->matrix, M->matrix, m, d);
    dd_CopyArow(Mcopy->rowvec, M->rowvec, d);
    set_copy(Mcopy->linset,M->linset);
    Mcopy->representation=M->representation;
    Mcopy->objective=M->objective;
  }
  return Mcopy;
}

template<typename T>
dd_matrixdata<T> *dd_CopyMatrix(dd_matrixdata<T> *M)
{
  return dd_MatrixCopy(M);
}

template<typename T>
dd_matrixdata<T> *dd_MatrixNormalizedCopy(dd_matrixdata<T> *M)
{
  dd_matrixdata<T> *Mcopy=nullptr;
  dd_rowrange m;
  dd_colrange d;

  m= M->rowsize;
  d= M->colsize;
  if (m >=0 && d >=0){
    Mcopy=dd_CreateMatrix<T>(m, d);
    dd_CopyNormalizedAmatrix(Mcopy->matrix, M->matrix, m, d);
    dd_CopyArow(Mcopy->rowvec, M->rowvec, d);
    set_copy(Mcopy->linset,M->linset);
    Mcopy->representation=M->representation;
    Mcopy->objective=M->objective;
  }
  return Mcopy;
}


template<typename T>
dd_matrixdata<T> *dd_MatrixAppend(dd_matrixdata<T> *M1, dd_matrixdata<T> *M2)
{
  dd_matrixdata<T> *M=nullptr;
  dd_rowrange i, m,m1,m2;
  dd_colrange j, d,d1,d2;

  m1=M1->rowsize;
  d1=M1->colsize;
  m2=M2->rowsize;
  d2=M2->colsize;

  m=m1+m2;
  d=d1;

  if (d1>=0 && d1==d2 && m1>=0 && m2>=0){
    M=dd_CreateMatrix<T>(m, d);
    dd_CopyAmatrix(M->matrix, M1->matrix, m1, d);
    dd_CopyArow(M->rowvec, M1->rowvec, d);
    for (i=0; i<m1; i++){
      if (set_member(i+1,M1->linset)) set_addelem(M->linset,i+1);
    }
    for (i=0; i<m2; i++){
       for (j=0; j<d; j++)
         dd_set(M->matrix[m1+i][j],M2->matrix[i][j]);
         /* append the second matrix */
       if (set_member(i+1,M2->linset)) set_addelem(M->linset,m1+i+1);
    }
  }
  return M;
}

void dd_RandomPermutation(dd_rowindex OV, long t, unsigned int seed)
{
  long k,j,ovj;
  double u,xk,r,rand_max=(double) RAND_MAX;

  srand(seed);
  for (j=t; j>1 ; j--) {
    r=rand();
    u=r/rand_max;
    xk=(double)(j*u +1);
    k=(long)xk;
    ovj=OV[j];
    OV[j]=OV[k];
    OV[k]=ovj;
  }
}

void dd_RandomPermutation2(dd_rowindex OV,long t,unsigned int seed)
{
  long k,j,ovj;
  double u,xk,r,rand_max=(double) RAND_MAX;

  srand(seed);
  for (j=t; j>1 ; j--) {
    r=rand();
    u=r/rand_max;
    xk=(double)(j*u +1);
    k=(long)xk;
    ovj=OV[j];
    OV[j]=OV[k];
    OV[k]=ovj;
  }
}

template<typename T>
bool dd_LexSmaller(T *v1, T *v2, long dmax)
{ /* dmax is the size of vectors v1,v2 */
  bool determined, smaller;
  dd_colrange j;

  smaller = false;
  determined = false;
  j = 1;
  do {
    if (!dd_Equal(v1[j - 1],v2[j - 1])) {  /* 086 */
      if (dd_Smaller(v1[j - 1],v2[j - 1])) {  /*086 */
	    smaller = true;
	  }
      determined = true;
    } else
      j++;
  } while (!(determined) && (j <= dmax));
  return smaller;
}

template<typename T>
bool dd_LexSmallerFrac(T *v1, T q1, T *v2, T q2, long dmax)
{ /* dmax is the size of vectors v1,v2 */
  bool determined, smaller;
  dd_colrange j;

  smaller = false;
  determined = false;
  j = 1;
  do {
    if (!dd_EqualFrac(v1[j - 1], q1, v2[j - 1], q2)) {  /* 086 */
      if (dd_SmallerFrac(v1[j - 1], q1, v2[j - 1], q2)) {  /*086 */
	    smaller = true;
	  }
      determined = true;
    } else
      j++;
  } while (!(determined) && (j <= dmax));
  return smaller;
}

template<typename T>
bool dd_LexLarger(T *v1, T *v2, long dmax)
{
  return dd_LexSmaller(v2, v1, dmax);
}

template<typename T>
long dd_Partition(dd_rowindex OV, long p, long r, T** A, long dmax)
{
  T *x;
  long i,j,ovi;

  x=A[OV[p]-1];

  i=p-1;
  j=r+1;
  while (true)
    {
      do
	{
	  j--;
	} while (dd_LexLarger(A[OV[j]-1],x,dmax));
      do
	{
	  i++;
	} while (dd_LexSmaller(A[OV[i]-1],x,dmax));
      if (i<j)
	{
	  ovi=OV[i];
	  OV[i]=OV[j];
	  OV[j]=ovi;
	}
      else
	{
	  return j;
	}
    }
  return -417;
}



template<typename T>
void dd_QuickSort(dd_rowindex OV, long p, long r, T** A, long dmax)
{
  long q;

  if (p < r){
    q = dd_Partition(OV, p, r, A, dmax);
    dd_QuickSort(OV, p, q, A, dmax);
    dd_QuickSort(OV, q+1, r, A, dmax);
  }
}

template<typename T>
dd_matrixdata<T> *dd_MatrixNormalizedSortedCopy(dd_matrixdata<T> *M,dd_rowindex *newpos)  /* 094 */
{
  /* Sort the rows of Amatrix lexicographically, and return a link to this sorted copy.
  The vector newpos is allocated, where newpos[i] returns the new row index
  of the original row i (i=1,...,M->rowsize). */
  dd_matrixdata<T> *Mcopy=nullptr,Mnorm=nullptr;
  dd_rowrange m,i;
  dd_colrange d;
  dd_rowindex roworder;

  /* if (newpos!=nullptr) free(newpos); */
  m= M->rowsize;
  d= M->colsize;
  roworder = new long[m+1];
  for (int i=0; i<=m; i++)
    roworder[i] = 0;
  *newpos = new long[m+1];
  for (int i=0; i<=m; i++)
    (*newpos)[i] = 0;
  if (m >=0 && d >=0){
    Mnorm=dd_MatrixNormalizedCopy(M);
    Mcopy=dd_CreateMatrix<T>(m, d);
    for(i=1; i<=m; i++) roworder[i]=i;

    dd_RandomPermutation(roworder, m, 123);
    dd_QuickSort(roworder,1,m,Mnorm->matrix,d);

    dd_PermuteCopyAmatrix(Mcopy->matrix, Mnorm->matrix, m, d, roworder);
    dd_CopyArow(Mcopy->rowvec, M->rowvec, d);
    for(i=1; i<=m; i++) {
      if (set_member(roworder[i],M->linset)) set_addelem(Mcopy->linset, i);
      (*newpos)[roworder[i]]=i;
    }
    Mcopy->representation=M->representation;
    Mcopy->objective=M->objective;
    dd_FreeMatrix(Mnorm);
  }
  delete [] roworder;
  return Mcopy;
}

template<typename T>
dd_matrixdata<T> *dd_MatrixUniqueCopy(dd_matrixdata<T> *M,dd_rowindex *newpos)
{
  /* Remove row duplicates, and return a link to this sorted copy.
     Linearity rows have priority over the other rows.
     It is better to call this after sorting with dd_MatrixNormalizedSortedCopy.
     The vector newpos is allocated, where *newpos[i] returns the new row index
     of the original row i (i=1,...,M->rowsize).  *newpos[i] is negative if the original
     row is dominated by -*newpos[i] and eliminated in the new copy.
  */
  dd_matrixdata<T> *Mcopy=nullptr;
  dd_rowrange m,i,uniqrows;
  dd_rowset preferredrows;
  dd_colrange d;
  dd_rowindex roworder;

  m= M->rowsize;
  d= M->colsize;
  preferredrows=M->linset;
  roworder = new long[m+1];
  for (int i=0; i<=m; i++)
    roworder[i] = 0;
  if (m >=0 && d >=0){
    for(i=1; i<=m; i++) roworder[i]=i;
    dd_UniqueRows(roworder, 1, m, M->matrix, d,preferredrows, &uniqrows);

    Mcopy=dd_CreateMatrix<T>(uniqrows, d);
    dd_PermutePartialCopyAmatrix(Mcopy->matrix, M->matrix, m, d, roworder);
    dd_CopyArow(Mcopy->rowvec, M->rowvec, d);
    for(i=1; i<=m; i++) {
      if (roworder[i]>0 && set_member(i,M->linset)) set_addelem(Mcopy->linset, roworder[i]);
    }
    Mcopy->representation=M->representation;
    Mcopy->objective=M->objective;
  }
  *newpos=roworder;
  return Mcopy;
}


template<typename T>
dd_matrixdata<T> *dd_MatrixNormalizedSortedUniqueCopy(dd_matrixdata<T> *M,dd_rowindex *newpos) /* 094 */
{
  /* Sort and remove row duplicates, and return a link to this sorted copy.
     Linearity rows have priority over the other rows.
     It is better to call this after sorting with dd_MatrixNormalizedSortedCopy.
     The vector newpos is allocated, where *newpos[i] returns the new row index
     of the original row i (i=1,...,M->rowsize).  *newpos[i] is negative if the original
     row is dominated by -*newpos[i] and eliminated in the new copy.
  */
  dd_matrixdata<T> *M1=nullptr,M2=nullptr;
  dd_rowrange m,i;
  dd_colrange d;
  dd_rowindex newpos1=nullptr,newpos1r=nullptr,newpos2=nullptr;

  m= M->rowsize;
  d= M->colsize;
  *newpos = new long[m+1];
  for (int i=0; i<=m; i++)
    (*newpos)[i] = 0;
  newpos1r = new long[m+1];
  for (int i=0; i<=m; i++)
    newpos1r[i] = 0;
  if (m>=0 && d>=0){
    M1=dd_MatrixNormalizedSortedCopy(M,&newpos1);
    for (i=1; i<=m;i++) newpos1r[newpos1[i]]=i;  /* reverse of newpos1 */
    M2=dd_MatrixUniqueCopy(M1,&newpos2);
    set_emptyset(M2->linset);
    for(i=1; i<=m; i++) {
      if (newpos2[newpos1[i]]>0){
         printf("newpos1[%ld]=%ld, newpos2[newpos1[%ld]]=%ld\n",i,newpos1[i], i,newpos2[newpos1[i]]);
         if (set_member(i,M->linset)) set_addelem(M2->linset, newpos2[newpos1[i]]);
         (*newpos)[i]=newpos2[newpos1[i]];
      } else {
         (*newpos)[i]=-newpos1r[-newpos2[newpos1[i]]];
      }
    }
    dd_FreeMatrix(M1);
    delete [] newpos1;
    delete [] newpos2;
  }
  delete [] newpos1r;
  return M2;
}

template<typename T>
dd_matrixdata<T> *dd_MatrixSortedUniqueCopy(dd_matrixdata<T> *M,dd_rowindex *newpos)  /* 094 */
{
  /* Same as dd_MatrixNormalizedSortedUniqueCopy except that it returns a unnormalized origial data
     with original ordering.
  */
  dd_matrixdata<T> *M1=nullptr,M2=nullptr;
  dd_rowrange m,i,k,ii;
  dd_colrange d;
  dd_rowindex newpos1=nullptr,newpos1r=nullptr,newpos2=nullptr;

  m= M->rowsize;
  d= M->colsize;
  *newpos = new long[m+1];
  for (int i=0; i<=m; i++)
    (*newpos)[i] = 0;
  newpos1r = new long[m+1];
  for (int i=0; i<=m; i++)
    newpos1r[i] = 0;
  if (m>=0 && d>=0){
    M1=dd_MatrixNormalizedSortedCopy(M,&newpos1);
    for (i=1; i<=m;i++) newpos1r[newpos1[i]]=i;  /* reverse of newpos1 */
    M2=dd_MatrixUniqueCopy(M1,&newpos2);
    dd_FreeMatrix(M1);
    set_emptyset(M2->linset);
    for(i=1; i<=m; i++) {
      if (newpos2[newpos1[i]]>0){
         if (set_member(i,M->linset)) set_addelem(M2->linset, newpos2[newpos1[i]]);
         (*newpos)[i]=newpos2[newpos1[i]];
      } else {
         (*newpos)[i]=-newpos1r[-newpos2[newpos1[i]]];
      }
    }

    ii=0;
    set_emptyset(M2->linset);
    for (i = 1; i<=m ; i++) {
      k=(*newpos)[i];
      if (k>0) {
        ii+=1;
        (*newpos)[i]=ii;
        dd_CopyArow(M2->matrix[ii-1],M->matrix[i-1],d);
        if (set_member(i,M->linset)) set_addelem(M2->linset, ii);
      }
    }

    delete [] newpos1;
    delete [] newpos2;
  }
  delete [] newpos1r;

  return M2;
}

template<typename T>
dd_matrixdata<T> *dd_AppendMatrix(dd_matrixdata<T> *M1, dd_matrixdata<T> *M2)
{
  return dd_MatrixAppend(M1,M2);
}

template<typename T>
int dd_MatrixAppendTo(dd_matrixdata<T> **M1, dd_matrixdata<T> *M2)
{
  dd_matrixdata<T> *M=nullptr;
  dd_rowrange i, m,m1,m2;
  dd_colrange j, d,d1,d2;
  bool success=0;

  m1=(*M1)->rowsize;
  d1=(*M1)->colsize;
  m2=M2->rowsize;
  d2=M2->colsize;

  m=m1+m2;
  d=d1;

  if (d1>=0 && d1==d2 && m1>=0 && m2>=0){
    M=dd_CreateMatrix<T>(m, d);
    dd_CopyAmatrix(M->matrix, (*M1)->matrix, m1, d);
    dd_CopyArow(M->rowvec, (*M1)->rowvec, d);
    for (i=0; i<m1; i++){
      if (set_member(i+1,(*M1)->linset)) set_addelem(M->linset,i+1);
    }
    for (i=0; i<m2; i++){
       for (j=0; j<d; j++)
         dd_set(M->matrix[m1+i][j],M2->matrix[i][j]);
         /* append the second matrix */
       if (set_member(i+1,M2->linset)) set_addelem(M->linset,m1+i+1);
    }
    dd_FreeMatrix(*M1);
    *M1=M;
    success=1;
  }
  return success;
}

template<typename T>
int dd_MatrixRowRemove(dd_matrixdata<T> **M, dd_rowrange r) /* 092 */
{
  dd_rowrange i,m;
  bool success=0;
  m=(*M)->rowsize;

  if (r >= 1 && r <=m){
    (*M)->rowsize=m-1;
    dd_FreeArow((*M)->matrix[r-1]);
    set_delelem((*M)->linset,r);
    /* slide the row headers */
    for (i=r; i<m; i++){
      (*M)->matrix[i-1]=(*M)->matrix[i];
      if (set_member(i+1, (*M)->linset)){
        set_delelem((*M)->linset,i+1);
        set_addelem((*M)->linset,i);
      }
    }
    success=1;
  }
  return success;
}

template<typename T>
int dd_MatrixRowRemove2(dd_matrixdata<T> **M, dd_rowrange r) /* 094 */
{
  dd_rowrange i,m;
  dd_colrange d;
  bool success=0;

  m=(*M)->rowsize;
  d=(*M)->colsize;

  if (r >= 1 && r <=m){
    (*M)->rowsize=m-1;
    dd_FreeArow((*M)->matrix[r-1]);
    set_delelem((*M)->linset,r);
    /* slide the row headers */
    for (i=r; i<m; i++){
      (*M)->matrix[i-1]=(*M)->matrix[i];
      if (set_member(i+1, (*M)->linset)){
        set_delelem((*M)->linset,i+1);
        set_addelem((*M)->linset,i);
      }
    }
    success=1;
  }
  return success;
}

template<typename T>
dd_matrixdata<T> *dd_MatrixSubmatrix(dd_matrixdata<T> *M, dd_rowset delset) /* 092 */
{
  dd_matrixdata<T> *Msub=nullptr;
  dd_rowrange i,isub=1, m,msub;
  dd_colrange d;

  m= M->rowsize;
  d= M->colsize;
  msub=m;
  if (m >=0 && d >=0){
    for (i=1; i<=m; i++) {
       if (set_member(i,delset)) msub-=1;
    }
    Msub=dd_CreateMatrix<T>(msub, d);
    for (i=1; i<=m; i++){
      if (!set_member(i,delset)){
        dd_CopyArow(Msub->matrix[isub-1], M->matrix[i-1], d);
        if (set_member(i, M->linset)){
          set_addelem(Msub->linset,isub);
        }
        isub++;
      }
    }
    dd_CopyArow(Msub->rowvec, M->rowvec, d);
    Msub->representation=M->representation;
    Msub->objective=M->objective;
  }
  return Msub;
}

template<typename T>
dd_matrixdata<T> *dd_MatrixSubmatrix2(dd_matrixdata<T> *M, dd_rowset delset,dd_rowindex *newpos) /* 092 */
{ /* returns a pointer to a new matrix which is a submatrix of M with rows in delset
  removed.  *newpos[i] returns the position of the original row i in the new matrix.
  It is -1 if and only if it is deleted.
  */

  dd_matrixdata<T> *Msub=nullptr;
  dd_rowrange i,isub=1, m,msub;
  dd_colrange d;
  dd_rowindex roworder;

  m= M->rowsize;
  d= M->colsize;
  msub=m;
  if (m >=0 && d >=0){
    roworder = new long[m+1];
    for (int i=0; i<=m; i++)
      roworder[i] = 0;
    for (i=1; i<=m; i++) {
       if (set_member(i,delset)) msub-=1;
    }
    Msub=dd_CreateMatrix<T>(msub, d);
    for (i=1; i<=m; i++){
      if (set_member(i,delset)){
        roworder[i]=0; /* zero means the row i is removed */
      } else {
        dd_CopyArow(Msub->matrix[isub-1], M->matrix[i-1], d);
        if (set_member(i, M->linset)){
          set_addelem(Msub->linset,isub);
        }
        roworder[i]=isub;
        isub++;
      }
    }
    *newpos=roworder;
    dd_CopyArow(Msub->rowvec, M->rowvec, d);
    Msub->representation=M->representation;
    Msub->objective=M->objective;
  }
  return Msub;
}


template<typename T>
dd_matrixdata<T> *dd_MatrixSubmatrix2L(dd_matrixdata<T> *M, dd_rowset delset,dd_rowindex *newpos) /* 094 */
{ /* This is same as dd_MatrixSubmatrix2 except that the linearity rows will be shifted up
     so that they are at the top of the matrix.
  */
  dd_matrixdata<T> *Msub=nullptr;
  dd_rowrange i,iL, iI, m,msub;
  dd_colrange d;
  dd_rowindex roworder;

  m= M->rowsize;
  d= M->colsize;
  msub=m;
  if (m >=0 && d >=0){
    roworder = new long[m+1];
    for (int i=0; i<=m; i++)
      roworder[i] = 0;
    for (i=1; i<=m; i++) {
       if (set_member(i,delset)) msub-=1;
    }
    Msub=dd_CreateMatrix<T>(msub, d);
    iL=1; iI=set_card(M->linset)+1;  /* starting positions */
    for (i=1; i<=m; i++){
      if (set_member(i,delset)){
        roworder[i]=0; /* zero means the row i is removed */
      } else {
        if (set_member(i,M->linset)){
          dd_CopyArow(Msub->matrix[iL-1], M->matrix[i-1], d);
          set_delelem(Msub->linset,i);
          set_addelem(Msub->linset,iL);
          roworder[i]=iL;
          iL+=1;
        } else {
          dd_CopyArow(Msub->matrix[iI-1], M->matrix[i-1], d);
          roworder[i]=iI;
          iI+=1;
        }
       }
    }
    *newpos=roworder;
    dd_CopyArow(Msub->rowvec, M->rowvec, d);
    Msub->representation=M->representation;
    Msub->objective=M->objective;
  }
  return Msub;
}

template<typename T>
int dd_MatrixRowsRemove(dd_matrixdata<T> **M, dd_rowset delset) /* 094 */
{
  dd_matrixdata<T> *Msub=nullptr;
  int success;

  Msub=dd_MatrixSubmatrix(*M, delset);
  dd_FreeMatrix(*M);
  *M=Msub;
  success=1;
  return success;
}

template<typename T>
int dd_MatrixRowsRemove2(dd_matrixdata<T> **M, dd_rowset delset,dd_rowindex *newpos) /* 094 */
{
  dd_matrixdata<T> *Msub=nullptr;
  int success;

  Msub=dd_MatrixSubmatrix2(*M, delset,newpos);
  dd_FreeMatrix(*M);
  *M=Msub;
  success=1;
  return success;
}

template<typename T>
int dd_MatrixShiftupLinearity(dd_matrixdata<T> **M,dd_rowindex *newpos) /* 094 */
{
  dd_matrixdata<T> *Msub=nullptr;
  int success;
  dd_rowset delset;

  set_initialize(&delset,(*M)->rowsize);  /* emptyset */
  Msub=dd_MatrixSubmatrix2L(*M, delset,newpos);
  dd_FreeMatrix(*M);
  *M=Msub;

  delete [] delset;
  success=1;
  return success;
}

template<typename T>
void dd_SetLinearity(dd_matrixdata<T> *M, char *line)
{
  int i=0;
  dd_rowrange eqsize,var;
  char *next;
  const char ct[]=", ";  /* allows separators "," and " ". */

  next=strtok(line,ct);
  eqsize=atol(next);
  while (i < eqsize && (next=strtok(nullptr,ct))!=nullptr) {
     var=atol(next);
     set_addelem(M->linset,var); i++;
  }
  if (i!=eqsize) {
    fprintf(stderr,"* Warning: there are inconsistencies in linearity setting.\n");
  }
  return;
}

template<typename T>
dd_polyhedradata<T> *dd_DDMatrix2Poly(dd_matrixdata<T> *M, dd_ErrorType *err)
{
  dd_rowrange i;
  dd_colrange j;
  dd_polyhedradata<T> *poly=nullptr;

  *err=dd_NoError;
  if (M->rowsize<0 || M->colsize<0){
    *err=dd_NegativeMatrixSize;
    goto _L99;
  }
  poly=dd_CreatePolyhedraData<T>(M->rowsize, M->colsize);
  poly->representation=M->representation;
  poly->homogeneous=true;

  for (i = 1; i <= M->rowsize; i++)
    {
      if (set_member(i, M->linset))
	poly->EqualityIndex[i]=1;
      for (j = 1; j <= M->colsize; j++) {
	dd_set(poly->A[i-1][j-1], M->matrix[i-1][j-1]);
	if (j==1 && dd_Nonzero(M->matrix[i-1][j-1])) poly->homogeneous = false;
      }
    }
  dd_DoubleDescription(poly,err);
_L99:
  return poly;
}

template<typename T>
dd_polyhedradata<T> *dd_DDMatrix2Poly2(dd_matrixdata<T> *M, dd_RowOrderType horder, dd_ErrorType *err)
{
  dd_rowrange i;
  dd_colrange j;
  dd_polyhedradata<T> *poly=nullptr;

  *err=dd_NoError;
  if (M->rowsize<0 || M->colsize<0){
    *err=dd_NegativeMatrixSize;
    goto _L99;
  }
  poly=dd_CreatePolyhedraData<T>(M->rowsize, M->colsize);
  poly->representation=M->representation;
  poly->homogeneous=true;

  for (i = 1; i <= M->rowsize; i++) {
    if (set_member(i, M->linset)) {
      poly->EqualityIndex[i]=1;
    }
    for (j = 1; j <= M->colsize; j++) {
      dd_set(poly->A[i-1][j-1], M->matrix[i-1][j-1]);
      if (j==1 && dd_Nonzero(M->matrix[i-1][j-1])) poly->homogeneous = false;
    }
  }
  dd_DoubleDescription2(poly, horder, err);
_L99:
  return poly;
}

template<typename T>
void dd_CopyRay(T *a, dd_colrange d_origsize, dd_raydata<T> *RR,
		dd_RepresentationType rep, dd_colindex reducedcol)
{
  long j,j1;
  T b;
  for (j = 1; j <= d_origsize; j++){
    j1=reducedcol[j];
    if (j1>0){
      dd_set(a[j-1],RR->Ray[j1-1]);
        /* the original column j is mapped to j1, and thus
           copy the corresponding component */
    } else {
      a[j-1]=0;
        /* original column is redundant and removed for computation */
    }
  }
}



template<typename T>
void dd_ComputeAinc(dd_polyhedradata<T> *poly)
{
/* This generates the input incidence array poly->Ainc, and
   two sets: poly->Ared, poly->Adom.
*/
  dd_bigrange k;
  dd_rowrange i,m1;
  dd_colrange j;
  bool redundant;
  dd_matrixdata<T> *M=nullptr;
  T sum,temp;

  if (poly->AincGenerated==true) return;

  M=dd_CopyOutput(poly);
  poly->n=M->rowsize;
  m1=poly->m1;
   /* this number is same as poly->m, except when
      poly is given by nonhomogeneous inequalty:
      !(poly->homogeneous) && poly->representation==Inequality,
      it is poly->m+1.   See dd_ConeDataLoad.
   */
  poly->Ainc = new set_type[m1];
  for(i=1; i<=m1; i++) set_initialize(&(poly->Ainc[i-1]),poly->n);
  set_initialize(&(poly->Ared), m1);
  set_initialize(&(poly->Adom), m1);

  for (k=1; k<=poly->n; k++){
    for (i=1; i<=poly->m; i++){
      sum=0;
      for (j=1; j<=poly->d; j++){
        dd_mul(temp,poly->A[i-1][j-1],M->matrix[k-1][j-1]);
        dd_add(sum,sum,temp);
      }
      if (dd_EqualToZero(sum)) {
        set_addelem(poly->Ainc[i-1], k);
      }
    }
    if (!(poly->homogeneous) && poly->representation==dd_Inequality){
      if (dd_EqualToZero(M->matrix[k-1][0])) {
        set_addelem(poly->Ainc[m1-1], k);  /* added infinity inequality (1,0,0,...,0) */
      }
    }
  }

  for (i=1; i<=m1; i++){
    if (set_card(poly->Ainc[i-1])==M->rowsize){
      set_addelem(poly->Adom, i);
    }
  }
  for (i=m1; i>=1; i--){
    if (set_card(poly->Ainc[i-1])==0){
      redundant=true;
      set_addelem(poly->Ared, i);
    }else {
      redundant=false;
      for (k=1; k<=m1; k++) {
        if (k!=i && !set_member(k, poly->Ared)  && !set_member(k, poly->Adom) &&
            set_subset(poly->Ainc[i-1], poly->Ainc[k-1])){
          if (!redundant){
            redundant=true;
          }
          set_addelem(poly->Ared, i);
        }
      }
    }
  }
  dd_FreeMatrix(M);
  poly->AincGenerated=true;
}

template<typename T>
bool dd_InputAdjacentQ(set_type &common,
			     long & lastn,
			     dd_polyhedradata<T> *poly,
			     dd_rowrange i1, dd_rowrange i2)
/* Before calling this function, RedundantSet must be
   a set of row indices whose removal results in a minimal
   nonredundant system to represent the input polyhedron,
   DominantSet must be the set of row indices which are
   active at every extreme points/rays.
*/
{
  bool adj=true;
  dd_rowrange i;

  if (poly->AincGenerated==false) dd_ComputeAinc<T>(poly);
  if (lastn!=poly->n){
    if (lastn >0)
      set_free(common);
    set_initialize(&common, poly->n);
    lastn=poly->n;
  }
  if (set_member(i1, poly->Ared) || set_member(i2, poly->Ared)){
    adj=false;
    goto _L99;
  }
  if (set_member(i1, poly->Adom) || set_member(i2, poly->Adom)){
  // dominant inequality is considered adjacencent to all others.
    adj=true;
    goto _L99;
  }
  set_int(common, poly->Ainc[i1-1], poly->Ainc[i2-1]);
  i=0;
  while (i<poly->m1 && adj==true){
    i++;
    if (i!=i1 && i!=i2 && !set_member(i, poly->Ared) &&
        !set_member(i, poly->Adom) && set_subset(common,poly->Ainc[i-1])){
      adj=false;
    }
  }
_L99:;
  return adj;
}


template<typename T>
dd_SetFamilyPtr dd_CopyIncidence(dd_polyhedradata<T> *poly)
{
  dd_SetFamilyPtr F=nullptr;
  dd_bigrange k;
  dd_rowrange i;

  if (poly->child==nullptr || poly->child->CompStatus!=dd_AllFound) goto _L99;
  if (poly->AincGenerated==false) dd_ComputeAinc(poly);
  F=dd_CreateSetFamily(poly->n, poly->m1);
  for (i=1; i<=poly->m1; i++)
    for (k=1; k<=poly->n; k++)
      if (set_member(k,poly->Ainc[i-1])) set_addelem(F->set[k-1],i);
_L99:;
  return F;
}

template<typename T>
dd_SetFamilyPtr dd_CopyInputIncidence(dd_polyhedradata<T> *poly)
{
  dd_rowrange i;
  dd_SetFamilyPtr F=nullptr;

  if (poly->child==nullptr || poly->child->CompStatus!=dd_AllFound) goto _L99;
  if (poly->AincGenerated==false) dd_ComputeAinc(poly);
  F=dd_CreateSetFamily(poly->m1, poly->n);
  for(i=0; i< poly->m1; i++){
    set_copy(F->set[i], poly->Ainc[i]);
  }
_L99:;
  return F;
}

template<typename T>
dd_SetFamilyPtr dd_CopyAdjacency(dd_polyhedradata<T> *poly)
{
  dd_raydata<T> *RayPtr1;
  dd_raydata<T> *RayPtr2;
  dd_SetFamilyPtr F=nullptr;
  long pos1, pos2;
  dd_bigrange lstart,k,n;
  set_type linset,allset;
  bool adj;

  if (poly->n==0 && poly->homogeneous && poly->representation==dd_Inequality){
    n=1; /* the origin (the unique vertex) should be output. */
  } else n=poly->n;
  set_initialize(&linset, n);
  set_initialize(&allset, n);
  if (poly->child==nullptr || poly->child->CompStatus!=dd_AllFound) goto _L99;
  F=dd_CreateSetFamily(n, n);
  if (n<=0) goto _L99;
  poly->child->LastRay->Next=nullptr;
  for (RayPtr1=poly->child->FirstRay, pos1=1;
       RayPtr1 != nullptr;
       RayPtr1 = RayPtr1->Next, pos1++) {
    for (RayPtr2=poly->child->FirstRay, pos2=1;
         RayPtr2 != nullptr;
         RayPtr2 = RayPtr2->Next, pos2++) {
      if (RayPtr1 != RayPtr2) {
        dd_CheckAdjacency(poly->child, &RayPtr1, &RayPtr2, &adj);
        if (adj){
          set_addelem(F->set[pos1-1], pos2);
        }
      }
    }
  }
  lstart=poly->n - poly->ldim + 1;
  set_compl(allset,allset);  /* allset is set to the ground set. */
  for (k=lstart; k<=poly->n; k++){
    set_addelem(linset,k);     /* linearity set */
    set_copy(F->set[k-1],allset);  /* linearity generator is adjacent to all */
  }
  for (k=1; k<lstart; k++){
    set_uni(F->set[k-1],F->set[k-1],linset);
     /* every generator is adjacent to all linearity generators */
  }
_L99:;
  set_free(allset); set_free(linset);
  return F;
}

template<typename T>
dd_SetFamilyPtr dd_CopyInputAdjacency(dd_polyhedradata<T> *poly)
{
  dd_rowrange i,j;
  dd_SetFamilyPtr F=nullptr;
  set_type common;
  long lastn=0;
  if (poly->child==nullptr || poly->child->CompStatus!=dd_AllFound) goto _L99;
  if (poly->AincGenerated==false) dd_ComputeAinc(poly);
  F=dd_CreateSetFamily(poly->m1, poly->m1);
  for (i=1; i<=poly->m1; i++)
    for (j=1; j<=poly->m1; j++)
      if (i!=j && dd_InputAdjacentQ<T>(common, lastn, poly, i, j))
        set_addelem(F->set[i-1],j);
_L99:;
  return F;
}


template<typename T>
dd_matrixdata<T> *dd_CopyOutput(dd_polyhedradata<T> *poly)
{
  dd_raydata<T> *RayPtr;
  dd_matrixdata<T> *M=nullptr;
  dd_rowrange i=0;
  dd_rowrange total;
  dd_colrange j, j1;
  T b;
  dd_RepresentationType outputrep=dd_Inequality;
  bool outputorigin=false;

  total=poly->child->LinearityDim + poly->child->FeasibleRayCount;

  if (poly->child->d<=0 || poly->child->newcol[1]==0) total=total-1;
  if (poly->representation==dd_Inequality) outputrep=dd_Generator;
  if (total==0 && poly->homogeneous && poly->representation==dd_Inequality){
    total=1;
    outputorigin=true;
    // the origin (the unique vertex) should be output.
  }
  if (poly->child==nullptr || poly->child->CompStatus!=dd_AllFound) goto _L99;

  M=dd_CreateMatrix<T>(total, poly->d);
  RayPtr = poly->child->FirstRay;
  while (RayPtr != nullptr) {
    if (RayPtr->feasible) {
      dd_CopyRay(M->matrix[i], poly->d, RayPtr, outputrep, poly->child->newcol);
      i++;
    }
    RayPtr = RayPtr->Next;
  }
  for (j=2; j<=poly->d; j++){
    if (poly->child->newcol[j]==0){
      for (j1=0; j1<poly->d; j1++)
        dd_set(M->matrix[i][j1],poly->child->Bsave[j1][j-1]);
      set_addelem(M->linset, i+1);
      i++;
    }
  }
  if (outputorigin){
    // output the origin for homogeneous H-polyhedron with no rays.
    M->matrix[0][0]=1;
    for (j=1; j<poly->d; j++){
      M->matrix[0][j]=0;
    }
  }
  //  dd_MatrixIntegerFilter(M);
  if (poly->representation==dd_Inequality)
    M->representation=dd_Generator;
  else
    M->representation=dd_Inequality;
_L99:;
  return M;
}


template<typename T>
dd_matrixdata<T> *dd_CopyInput(dd_polyhedradata<T> *poly)
{
  dd_matrixdata<T> *M=nullptr;
  dd_rowrange i;

  M=dd_CreateMatrix<T>(poly->m, poly->d);
  dd_CopyAmatrix(M->matrix, poly->A, poly->m, poly->d);
  for (i=1; i<=poly->m; i++)
    if (poly->EqualityIndex[i]==1) set_addelem(M->linset,i);
  // dd_MatrixIntegerFilter(M);
  if (poly->representation==dd_Generator)
    M->representation=dd_Generator;
  else
    M->representation=dd_Inequality;
  return M;
}

template<typename T>
dd_matrixdata<T> *dd_CopyGenerators(dd_polyhedradata<T> *poly)
{
  dd_matrixdata<T> *M=nullptr;

  if (poly->representation==dd_Generator){
    M=dd_CopyInput(poly);
  } else {
    M=dd_CopyOutput(poly);
  }
  return M;
}

template<typename T>
dd_matrixdata<T> *dd_CopyInequalities(dd_polyhedradata<T> *poly)
{
  dd_matrixdata<T> *M=nullptr;

  if (poly->representation==dd_Inequality){
    M=dd_CopyInput(poly);
  } else {
    M=dd_CopyOutput(poly);
  }
  return M;
}


long set_groundsize(set_type set)
{
	return set[0];
}



template<typename T>
dd_matrixdata<T> *dd_BlockElimination(dd_matrixdata<T> *M, dd_colset delset, dd_ErrorType *error)
/* Eliminate the variables (columns) delset by
   the Block Elimination with dd_DoubleDescription algorithm.

   Given (where y is to be eliminated):
   c1 + A1 x + B1 y >= 0
   c2 + A2 x + B2 y =  0

   1. First construct the dual system:  z1^T B1 + z2^T B2 = 0, z1 >= 0.
   2. Compute the generators of the dual.
   3. Then take the linear combination of the original system with each generator.
   4. Remove redundant inequalies.

*/
{
  dd_matrixdata<T> *Mdual=nullptr, Mproj=nullptr, Gdual=nullptr;
  dd_rowrange i,h,m,mproj,mdual,linsize;
  dd_colrange j,k,d,dproj,ddual,delsize;
  dd_colindex delindex;
  T temp,prod;
  dd_polyhedradata<T> *dualpoly;
  dd_ErrorType err=dd_NoError;
  bool localdebug=false;

  *error=dd_NoError;
  m= M->rowsize;
  d= M->colsize;
  delindex = new long[d+1];
  for (int i=0; i<=d; i++)
    delindex[i] = 0;
  k=0; delsize=0;
  for (j=1; j<=d; j++){
    if (set_member(j, delset)){
      k++;  delsize++;
      delindex[k]=j;  /* stores the kth deletion column index */
    }
  }

  linsize=set_card(M->linset);
  ddual=m+1;
  mdual=delsize + m - linsize;  /* #equalitions + dimension of z1 */

  /* setup the dual matrix */
  Mdual=dd_CreateMatrix<T>(mdual, ddual);
  Mdual->representation=dd_Inequality;
  for (i = 1; i <= delsize; i++){
    set_addelem(Mdual->linset,i);  /* equality */
    for (j = 1; j <= m; j++) {
      dd_set(Mdual->matrix[i-1][j], M->matrix[j-1][delindex[i]-1]);
    }
  }

  k=0;
  for (i = 1; i <= m; i++){
    if (!set_member(i, M->linset)){
      /* set nonnegativity for the dual variable associated with
         each non-linearity inequality. */
      k++;
      Mdual->matrix[delsize+k-1][i]=1;
    }
  }

  /* 2. Compute the generators of the dual system. */
  dualpoly=dd_DDMatrix2Poly(Mdual, &err);
  Gdual=dd_CopyGenerators(dualpoly);

  /* 3. Take the linear combination of the original system with each generator.  */
  dproj=d-delsize;
  mproj=Gdual->rowsize;
  Mproj=dd_CreateMatrix<T>(mproj, dproj);
  Mproj->representation=dd_Inequality;
  set_copy(Mproj->linset, Gdual->linset);

  for (i=1; i<=mproj; i++){
    k=0;
    for (j=1; j<=d; j++){
      if (!set_member(j, delset)){
        k++;  /* new index of the variable x_j  */
        prod=0;
        for (h = 1; h <= m; h++){
          dd_mul(temp,M->matrix[h-1][j-1],Gdual->matrix[i-1][h]);
          dd_add(prod,prod,temp);
        }
        dd_set(Mproj->matrix[i-1][k-1],prod);
      }
    }
  }
  if (localdebug) printf("Size of the projection system: %ld x %ld\n", mproj, dproj);

  dd_FreePolyhedra(dualpoly);
  delete [] delindex;
  dd_FreeMatrix(Mdual);
  dd_FreeMatrix(Gdual);
  return Mproj;
}



template<typename T>
inline void dd_LinearComb(T &lc, T v1, T c1, T v2, T c2)
{
  T temp;
  dd_mul(lc,v1,c1);
  dd_mul(temp,v2,c2);
  dd_add(lc,lc,temp);
}

template<typename T>
dd_matrixdata<T> *dd_FourierElimination(dd_matrixdata<T> *M,dd_ErrorType *error)
/* Eliminate the last variable (column) from the given H-matrix using
   the standard Fourier Elimination.
 */
{
  dd_matrixdata<T> *Mnew=nullptr;
  dd_rowrange i,inew,ip,in,iz,m,mpos=0,mneg=0,mzero=0,mnew;
  dd_colrange j,d,dnew;
  dd_rowindex posrowindex, negrowindex,zerorowindex;
  T temp1,temp2;
  bool localdebug=false;

  *error=dd_NoError;
  m= M->rowsize;
  d= M->colsize;
  if (d<=1){
    *error=dd_ColIndexOutOfRange;
    if (localdebug) {
      printf("The number of column is too small: %ld for Fourier's Elimination.\n",d);
    }
    goto _L99;
  }

  if (M->representation==dd_Generator){
    *error=dd_NotAvailForV;
    if (localdebug) {
      printf("Fourier's Elimination cannot be applied to a V-polyhedron.\n");
    }
    goto _L99;
  }

  if (set_card(M->linset)>0){
    *error=dd_CannotHandleLinearity;
    if (localdebug) {
      printf("The Fourier Elimination function does not handle equality in this version.\n");
    }
    goto _L99;
  }

  /* Create temporary spaces to be removed at the end of this function */
  posrowindex = new long[m+1];
  for (int i=0; i<=m; i++)
    posrowindex[i] = 0;
  negrowindex = new long[m+1];
  for (int i=0; i<=m; i++)
    negrowindex[i] = 0;
  zerorowindex = new long[m+1];
  for (int i=0; i<=m; i++)
    zerorowindex[i] = 0;

  for (i = 1; i <= m; i++) {
    if (dd_Positive(M->matrix[i-1][d-1])){
      mpos++;
      posrowindex[mpos]=i;
    } else if (dd_Negative(M->matrix[i-1][d-1])) {
      mneg++;
      negrowindex[mneg]=i;
    } else {
      mzero++;
      zerorowindex[mzero]=i;
    }
  }  /*of i*/


  /* The present code generates so many redundant inequalities and thus
     is quite useless, except for very small examples
  */
  mnew=mzero+mpos*mneg;  /* the total number of rows after elimination */
  dnew=d-1;

  Mnew=dd_CreateMatrix<T>(mnew, dnew);
  dd_CopyArow(Mnew->rowvec, M->rowvec, dnew);
/*  set_copy(Mnew->linset,M->linset);  */
  Mnew->representation=M->representation;
  Mnew->objective=M->objective;


  /* Copy the inequalities independent of x_d to the top of the new matrix. */
  for (iz = 1; iz <= mzero; iz++){
    for (j = 1; j <= dnew; j++) {
      dd_set(Mnew->matrix[iz-1][j-1], M->matrix[zerorowindex[iz]-1][j-1]);
    }
  }

  /* Create the new inequalities by combining x_d positive and negative ones. */
  inew=mzero;  /* the index of the last x_d zero inequality */
  for (ip = 1; ip <= mpos; ip++){
    for (in = 1; in <= mneg; in++){
      inew++;
      dd_neg(temp1, M->matrix[negrowindex[in]-1][d-1]);
      for (j = 1; j <= dnew; j++) {
        dd_LinearComb(temp2,M->matrix[posrowindex[ip]-1][j-1],temp1,\
          M->matrix[negrowindex[in]-1][j-1],\
          M->matrix[posrowindex[ip]-1][d-1]);
        dd_set(Mnew->matrix[inew-1][j-1],temp2);
      }
    }
  }


  delete [] posrowindex;
  delete [] negrowindex;
  delete [] zerorowindex;

 _L99:
  return Mnew;
}


template<typename T>
dd_lpdata<T> *dd_Matrix2LP(dd_matrixdata<T> *M, dd_ErrorType *err)
{
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_lpdata<T> *lp;

  *err=dd_NoError;
  linc=set_card(M->linset);
  m=M->rowsize+1+linc;
     /* We represent each equation by two inequalities.
        This is not the best way but makes the code simple. */
  d=M->colsize;

  lp=dd_CreateLPData<T>(m, d);
  lp->Homogeneous = true;
  lp->objective = M->objective;
  lp->eqnumber=linc;  /* this records the number of equations */

  irev=M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    if (set_member(i, M->linset)) {
      irev=irev+1;
      set_addelem(lp->equalityset,i);    /* it is equality. */
                                         /* the reversed row irev is not in the equality set. */
      for (j = 1; j <= M->colsize; j++) {
        dd_neg(lp->A[irev-1][j-1],M->matrix[i-1][j-1]);
      }  /*of j*/
    }
    for (j = 1; j <= M->colsize; j++) {
      dd_set(lp->A[i-1][j-1],M->matrix[i-1][j-1]);
      if (j==1 && i<M->rowsize && dd_Nonzero(M->matrix[i-1][j-1])) lp->Homogeneous = false;
    }  /*of j*/
  }  /*of i*/
  for (j = 1; j <= M->colsize; j++) {
    dd_set(lp->A[m-1][j-1],M->rowvec[j-1]);  /* objective row */
  }  /*of j*/

  return lp;
}

template<typename T>
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

  *err=dd_NoError;
  linc=set_card(M->linset);
  m=M->rowsize+1+linc;
     /* We represent each equation by two inequalities.
        This is not the best way but makes the code simple. */

  lp=dd_Matrix2LP(M, err);
  lp->objective = dd_LPmax;   /* since the objective is zero, this is not important */
  for (j = 1; j <= M->colsize; j++) {
    lp->A[m-1][j-1]=0;  /* set the objective to zero. */
  }  /*of j*/

  return lp;
}

template<typename T>
dd_lpdata<T> *dd_Matrix2Feasibility2(dd_matrixdata<T> *M, dd_rowset R, dd_rowset S, dd_ErrorType *err)
/* Load a matrix to create an LP object for feasibility with additional equality and
   strict inequality constraints given by R and S.  There are three types of inequalities:

   b_r + A_r x =  0     Linearity (Equations) specified by M
   b_s + A_s x >  0     Strict Inequalities specified by row index set S
   b_t + A_t x >= 0     The rest inequalities in M

   Where the linearity is considered here as the union of linearity specified by
   M and the additional set R.  When S contains any linearity rows, those
   rows are considered linearity (equation).  Thus S does not overlide linearity.
   To find a feasible solution, we set an LP

   maximize  z
   subject to
   b_r + A_r x     =  0      all r in Linearity
   b_s + A_s x - z >= 0      for all s in S
   b_t + A_t x     >= 0      for all the rest rows t
   1           - z >= 0      to make the LP bounded.

   Clearly, the feasibility problem has a solution iff the LP has a positive optimal value.
   The variable z will be the last variable x_{d+1}.

*/
/*  094 */
{
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_lpdata<T> *lp;
  dd_rowset L;

  *err=dd_NoError;
  set_initialize(&L, M->rowsize);
  set_uni(L,M->linset,R);
  linc=set_card(L);
  m=M->rowsize+1+linc+1;
     /* We represent each equation by two inequalities.
        This is not the best way but makes the code simple. */
  d=M->colsize+1;

  lp=dd_CreateLPData<T>(m, d);
  lp->Homogeneous = true;
  lp->objective = dd_LPmax;
  lp->eqnumber=linc;  /* this records the number of equations */

  irev=M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    if (set_member(i, L)) {
      irev=irev+1;
      set_addelem(lp->equalityset,i);    /* it is equality. */
                                         /* the reversed row irev is not in the equality set. */
      for (j = 1; j <= M->colsize; j++) {
        dd_neg(lp->A[irev-1][j-1],M->matrix[i-1][j-1]);
      }  /*of j*/
    } else if (set_member(i, S)) {
	  lp->A[i-1][M->colsize]=-1;
    }
    for (j = 1; j <= M->colsize; j++) {
      dd_set(lp->A[i-1][j-1],M->matrix[i-1][j-1]);
      if (j==1 && i<M->rowsize && dd_Nonzero(M->matrix[i-1][j-1])) lp->Homogeneous = false;
    }  /*of j*/
  }  /*of i*/
  for (j = 1; j <= d; j++) {
    lp->A[m-2][j-1]=0;  /* initialize */
  }  /*of j*/
  lp->A[m-2][0]=1;  /* the bounding constraint. */
  lp->A[m-2][M->colsize]=-1;  /* the bounding constraint. */
  for (j = 1; j <= d; j++) {
    lp->A[m-1][j-1]=0;  /* initialize */
  }  /*of j*/
  lp->A[m-1][M->colsize]=1;

  set_free(L);
  return lp;
}

template<typename T>
void dd_FreeLPData(dd_lpdata<T> *lp)
{
  if ((lp)!=nullptr){
    dd_FreeArow(lp->dsol);
    dd_FreeArow(lp->sol);
    dd_FreeBmatrix(lp->d_alloc,lp->B);
    dd_FreeAmatrix(lp->m_alloc,lp->A);
    set_free(lp->equalityset);
    set_free(lp->redset_extra);
    set_free(lp->redset_accum);
    set_free(lp->posset_extra);
    delete [] lp->nbindex;
    delete [] lp->given_nbindex;
    delete lp;
  }
}

template<typename T>
void dd_FreeLPSolution(dd_lpsolution<T> *lps)
{
  if (lps!=nullptr){
    delete [] lps->nbindex;
    dd_FreeArow(lps->dsol);
    dd_FreeArow(lps->sol);

    delete lps;
  }
}

template<typename T>
int dd_LPReverseRow(dd_lpdata<T> *lp, dd_rowrange i)
{
  dd_colrange j;
  int success=0;

  if (i>=1 && i<=lp->m){
    lp->LPS=dd_LPSundecided;
    for (j=1; j<=lp->d; j++) {
      dd_neg(lp->A[i-1][j-1],lp->A[i-1][j-1]);
      /* negating the i-th constraint of A */
    }
    success=1;
  }
  return success;
}

template<typename T>
int dd_LPReplaceRow(dd_lpdata<T> *lp, dd_rowrange i, T* a)
{
  dd_colrange j;
  int success=0;

  if (i>=1 && i<=lp->m){
    lp->LPS=dd_LPSundecided;
    for (j=1; j<=lp->d; j++) {
      dd_set(lp->A[i-1][j-1],a[j-1]);
      /* replacing the i-th constraint by a */
    }
    success=1;
  }
  return success;
}


template<typename T>
void dd_TableauEntry(T & x, dd_colrange d_size, T** X, T** Ts,
                     dd_rowrange r, dd_colrange s)
/* Compute the (r,s) entry of X.T   */
{
  dd_colrange j;
  T temp;

  x=0;
  for (j=0; j< d_size; j++) {
    dd_mul(temp,X[r-1][j], Ts[j][s-1]);
    x += temp;
  }
}

template<typename T>
void dd_GetRedundancyInformation(dd_rowrange m_size,dd_colrange d_size,T** A,T** Ts,
				 dd_rowindex bflag, dd_rowset redset)
/* Some basic variables that are forced to be nonnegative will be output.  These are
   variables whose dictionary row components are all nonnegative.   */
{
  dd_colrange j;
  dd_rowrange i;
  T x;
  bool red=false;
  long numbred=0;

  for (i=1; i<= m_size; i++) {
    red=true;
    for (j=1; j<= d_size; j++) {
      dd_TableauEntry(x, d_size,A,Ts,i,j);
      if (red && dd_Negative(x)) red=false;
    }
    if (bflag[i]<0 && red) {
      numbred+=1;
      set_addelem(redset,i);
    }
  }
}


template<typename T>
void dd_SelectDualSimplexPivot(dd_rowrange m_size,dd_colrange d_size,
    int Phase1, T** A, T** Ts,
    dd_colindex nbindex_ref, dd_rowindex bflag,
    dd_rowrange objrow,dd_colrange rhscol, bool lexicopivot,
    dd_rowrange *r, dd_colrange *s, bool *selected, dd_LPStatusType *lps,
    data_temp_simplex<T>* data)
{
  /* selects a dual simplex pivot (*r,*s) if the current
     basis is dual feasible and not optimal. If not dual feasible,
     the procedure returns *selected=false and *lps=LPSundecided.
     If Phase1=true, the RHS column will be considered as the negative
     of the column of the largest variable (==m_size).  For this case, it is assumed
     that the caller used the auxiliary row (with variable m_size) to make the current
     dictionary dual feasible before calling this routine so that the nonbasic
     column for m_size corresponds to the auxiliary variable.
  */
  bool colselected=false,rowselected=false, dualfeasible=true;
  dd_rowrange i,iref;
  dd_colrange j,k;
  T val, valn, minval=0, rat, minrat=0, minrat_q = 1;
  //  T* rcost;
  //  dd_colset tieset;
  //  dd_colset stieset;  /* store the column indices with tie */
  //  rcost=new T[d_size];
  //  set_initialize(&tieset,d_size);
  //  set_initialize(&stieset,d_size);

  *r=0;
  *s=0;
  *selected=false;
  *lps=dd_LPSundecided;
  for (j=1; j<=d_size; j++){
    if (j!=rhscol){
      dd_TableauEntry(data->rcost[j-1], d_size,A,Ts,objrow,j);
      if (dd_Positive(data->rcost[j-1])) {
        dualfeasible=false;
      }
    }
  }
  if (dualfeasible){
    while ((*lps==dd_LPSundecided) && (!rowselected) && (!colselected)) {
      for (i=1; i<=m_size; i++) {
        if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          if (Phase1){
            dd_TableauEntry(val, d_size,A,Ts,i,bflag[m_size]);
            dd_neg(val,val);
            /* for dual Phase I.  The RHS (dual objective) is the negative of the auxiliary variable column. */
          }
          else {dd_TableauEntry(val, d_size,A,Ts,i,rhscol);}
          if (dd_Smaller(val,minval)) {
            *r=i;
            dd_set(minval,val);
          }
        }
      }
      if (dd_Nonnegative(minval)) {
        *lps=dd_Optimal;
      }
      else {
        rowselected=true;
        set_emptyset(data->tieset);
        for (j=1; j<=d_size; j++){
          dd_TableauEntry(val, d_size,A,Ts,*r,j);
          if (j!=rhscol && dd_Positive(val)) {
            bool is_field=true;
            if (is_field) {
                rat = - data->rcost[j-1] / val;
                if (*s==0 || dd_Smaller(rat, minrat)) {
                  dd_set(minrat,rat);
                  *s=j;
                  set_emptyset(data->tieset);
                  set_addelem(data->tieset, j);
                } else if (dd_Equal(rat, minrat)) {
                  set_addelem(data->tieset,j);
                }
              } else { // ring case
              rat = - data->rcost[j-1];
              if (*s==0 || dd_SmallerFrac(rat, val, minrat, minrat_q)){
                dd_set(minrat  ,rat);
                dd_set(minrat_q,val);
                *s=j;
                set_emptyset(data->tieset);
                set_addelem(data->tieset, j);
              } else if (dd_EqualFrac(rat, val, minrat, minrat_q)){
                set_addelem(data->tieset,j);
              }
            }
          }
        }
        if (*s>0) {
          if (!lexicopivot || set_card(data->tieset)==1){
            colselected=true; *selected=true;
          } else { /* lexicographic rule with respect to the given reference cobasis.  */
            *s=0;
            k=2; /* k runs through the column indices except RHS. */
            do {
              iref=nbindex_ref[k];  /* iref runs though the reference basic indices */
              if (iref>0) {
                j=bflag[iref];
                if (j>0) {
                  if (set_member(j,data->tieset) && set_card(data->tieset)==1) {
                    *s=j;
                     colselected=true;
                  } else {
                    set_delelem(data->tieset, j);
                    /* iref is cobasic, and the corresponding col is not the pivot column except it is the last one. */
                  }
                } else {
                  *s=0;
                  for (j=1; j<=d_size; j++){
                    if (set_member(j,data->tieset)) {
                      dd_TableauEntry(val, d_size,A,Ts,*r,j);
                      dd_TableauEntry(valn, d_size,A,Ts,iref,j);
                      if (j!=rhscol && dd_Positive(val)) {
                        rat = valn / val;
                        if (*s==0 || dd_Smaller(rat,minrat)) {
                          dd_set(minrat,rat);
                          *s=j;
                          set_emptyset(data->stieset);
                          set_addelem(data->stieset, j);
                        } else if (dd_Equal(rat,minrat)){
                          set_addelem(data->stieset,j);
                        }
                      }
                    }
                  }
                  set_copy(data->tieset,data->stieset);
                  if (set_card(data->tieset)==1) colselected=true;
                }
              }
              k+=1;
            } while (!colselected && k<=d_size);
            *selected=true;
          }
        } else *lps=dd_Inconsistent;
      }
    } /* end of while */
  }
}


void dd_SelectPreorderedNext2(dd_rowrange m_size,
    rowset excluded,dd_rowindex OV,dd_rowrange *hnext)
{
  dd_rowrange i,k;

  *hnext=0;
  for (i=1; i<=m_size && *hnext==0; i++){
    k=OV[i];
    if (!set_member(k,excluded)) *hnext=k ;
  }
}



template<typename T>
void dd_SelectPivot2(dd_rowrange m_size,dd_colrange d_size,T** A,T** Ts,
		     dd_rowindex ordervec, rowset equalityset,
		     dd_rowrange rowmax,rowset NopivotRow,
		     colset NopivotCol,dd_rowrange *r,dd_colrange *s,
		     bool *selected)
/* Select a position (*r,*s) in the matrix A.T such that (A.T)[*r][*s] is nonzero
   The choice is feasible, i.e., not on NopivotRow and NopivotCol, and
   best with respect to the specified roworder
 */
{
  int stop;
  dd_rowrange i,rtemp;
  rowset rowexcluded;
  T Xtemp;

  stop = false;
  set_initialize(&rowexcluded,m_size);
  set_copy(rowexcluded,NopivotRow);
  for (i=rowmax+1;i<=m_size;i++) {
    set_addelem(rowexcluded,i);   /* cannot pivot on any row > rmax */
  }
  *selected = false;
  do {
    rtemp=0; i=1;
    while (i<=m_size && rtemp==0) {  /* equalityset vars have highest priorities */
      if (set_member(i,equalityset) && !set_member(i,rowexcluded)){
        rtemp=i;
      }
      i++;
    }
    if (rtemp==0) dd_SelectPreorderedNext2(m_size, rowexcluded,ordervec,&rtemp);;
    if (rtemp>=1) {
      *r=rtemp;
      *s=1;
      while (*s <= d_size && !*selected) {
        dd_TableauEntry(Xtemp, d_size,A,Ts,*r,*s);
        if (!set_member(*s,NopivotCol) && dd_Nonzero(Xtemp)) {
          *selected = true;
          stop = true;
        } else {
          (*s)++;
        }
      }
      if (!*selected) {
        set_addelem(rowexcluded,rtemp);
      }
    }
    else {
      *r = 0;
      *s = 0;
      stop = true;
    }
  } while (!stop);
  set_free(rowexcluded);
}

template<typename T>
void dd_GaussianColumnPivot(dd_colrange d_size,
    T** X, T** Ts, dd_rowrange r, dd_colrange s, T* Rtemp)
/* Update the Transformation matrix T with the pivot operation on (r,s)
   This procedure performs a implicit pivot operation on the matrix X by
   updating the dual basis inverse  T.
 */
{
  dd_colrange j, j1;
  T Xtemp0, Xtemp1, Xtemp;

  for (j=1; j<=d_size; j++) {
    dd_TableauEntry(Rtemp[j-1], d_size, X, Ts, r, j);
  }
  bool is_field = true;
  if (is_field) {
      dd_set(Xtemp0, Rtemp[s-1]);
      for (j = 1; j <= d_size; j++) {
        if (j != s) {
          dd_div(Xtemp, Rtemp[j-1], Xtemp0);
          Xtemp1=0;
          for (j1 = 1; j1 <= d_size; j1++){
            dd_mul(Xtemp1, Xtemp, Ts[j1-1][s-1]);
            dd_sub(Ts[j1-1][j-1], Ts[j1-1][j-1], Xtemp1);
          }
        }
      }
      for (j = 1; j <= d_size; j++)
        dd_div(Ts[j-1][s - 1], Ts[j-1][s - 1], Xtemp0);
    } else { // ring case now
    for (j = 1; j <= d_size; j++) {
      if (j != s) {
        for (j1 = 1; j1 <= d_size; j1++) {
          Ts[j1-1][j-1] = Rtemp[s-1] * Ts[j1-1][j-1] - Rtemp[j-1] * Ts[j1-1][s-1];
        }
      }
    }
    dd_set(Xtemp0, Rtemp[s-1]);
    for (j = 1; j <= d_size; j++)
      dd_div(Ts[j-1][s - 1], Ts[j-1][s - 1], Xtemp0);
    for (j=1; j<=d_size; j++) {
      if (j != s) {
        T alpha;
        dd_TableauEntry(alpha, d_size, X, Ts, r, j);
        if (alpha != 0)
          std::cout << "j=" << j << " alpha=" << alpha << "\n";
      }
    }
  }
}

template<typename T>
void dd_GaussianColumnPivot2(dd_colrange d_size,
    T** A,T** Ts,dd_colindex nbindex,dd_rowindex bflag,dd_rowrange r,dd_colrange s, T* Rtemp)
/* Update the Transformation matrix T with the pivot operation on (r,s)
   This procedure performs a implicit pivot operation on the matrix A by
   updating the dual basis inverse  T.
 */
{
  long entering;

  dd_GaussianColumnPivot(d_size, A, Ts, r, s, Rtemp);
  entering=nbindex[s];
  bflag[r]=s;     /* the nonbasic variable r corresponds to column s */
  nbindex[s]=r;   /* the nonbasic variable on s column is r */

  if (entering>0) bflag[entering]=-1;
     /* original variables have negative index and should not affect the row index */

}



template<typename T>
void dd_SetToIdentity(dd_colrange d_size, T** Ts)
{
  dd_colrange j1, j2;

  for (j1 = 1; j1 <= d_size; j1++) {
    for (j2 = 1; j2 <= d_size; j2++) {
      if (j1 == j2)
        Ts[j1 - 1][j2 - 1]=1;
      else
        Ts[j1 - 1][j2 - 1]=0;
    }
  }
}


template<typename T>
void dd_ResetTableau(dd_rowrange m_size,dd_colrange d_size,T** Ts,
    dd_colindex nbindex,dd_rowindex bflag,dd_rowrange objrow,dd_colrange rhscol)
{
  dd_rowrange i;
  dd_colrange j;

  /* Initialize T and nbindex */
  for (j=1; j<=d_size; j++) nbindex[j]=-j;
  nbindex[rhscol]=0;
    /* RHS is already in nonbasis and is considered to be associated
       with the zero-th row of input. */
  dd_SetToIdentity(d_size,Ts);

  /* Set the bflag according to nbindex */
  for (i=1; i<=m_size; i++) bflag[i]=-1;
    /* all basic variables have index -1 */
  bflag[objrow]= 0;
    /* bflag of the objective variable is 0,
       different from other basic variables which have -1 */
  for (j=1; j<=d_size; j++) if (nbindex[j]>0) bflag[nbindex[j]]=j;
    /* bflag of a nonbasic variable is its column number */

}

template<typename T>
void dd_SelectCrissCrossPivot(dd_rowrange m_size,dd_colrange d_size,
			      T** A,T** Ts,
			      dd_rowindex bflag,dd_rowrange objrow,dd_colrange rhscol,
			      dd_rowrange *r,dd_colrange *s,
			      bool *selected,dd_LPStatusType *lps)
{
  int colselected=false,rowselected=false;
  dd_rowrange i;
  T val;

  *selected=false;
  *lps=dd_LPSundecided;
  while ((*lps==dd_LPSundecided) && (!rowselected) && (!colselected)) {
    for (i=1; i<=m_size; i++) {
      if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
        dd_TableauEntry(val, d_size, A, Ts, i, rhscol);
        if (dd_Negative(val)) {
          rowselected=true;
          *r=i;
          break;
        }
      }
      else if (bflag[i] >0) { /* i is nonbasic variable */
        dd_TableauEntry(val, d_size,A,Ts,objrow,bflag[i]);
        if (dd_Positive(val)) {
          colselected=true;
          *s=bflag[i];
          break;
        }
      }
    }
    if  ((!rowselected) && (!colselected)) {
      *lps=dd_Optimal;
      return;
    }
    else if (rowselected) {
     for (i=1; i<=m_size; i++) {
       if (bflag[i] >0) { /* i is nonbasic variable */
          dd_TableauEntry(val, d_size,A,Ts,*r,bflag[i]);
          if (dd_Positive(val)) {
            colselected=true;
            *s=bflag[i];
            *selected=true;
            break;
          }
        }
      }
    }
    else if (colselected) {
      for (i=1; i<=m_size; i++) {
        if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          dd_TableauEntry(val, d_size,A,Ts,i,*s);
          if (dd_Negative(val)) {
            rowselected=true;
            *r=i;
            *selected=true;
            break;
          }
        }
      }
    }
    if (!rowselected) {
      *lps=dd_DualInconsistent;
    }
    else if (!colselected) {
      *lps=dd_Inconsistent;
    }
  }
}

template<typename T>
void dd_CrissCrossSolve(dd_lpdata<T> *lp, dd_ErrorType *err)
{
  data_temp_simplex<T>* data = allocate_data_simplex<T>(lp->d);
  switch (lp->objective) {
    case dd_LPmax:
      dd_CrissCrossMaximize(lp, err, data);
      break;

    case dd_LPmin:
      dd_CrissCrossMinimize(lp, err, data);
      break;

    case dd_LPnone: *err=dd_NoLPObjective; break;
  }
  free_data_simplex(data);
}

template<typename T>
void dd_DualSimplexSolve(dd_lpdata<T> *lp, dd_ErrorType *err)
{
  data_temp_simplex<T>* data = allocate_data_simplex<T>(lp->d);
  switch (lp->objective) {
    case dd_LPmax:
      dd_DualSimplexMaximize(lp, err, data);
      break;

    case dd_LPmin:
      dd_DualSimplexMinimize(lp, err, data);
      break;

    case dd_LPnone: *err=dd_NoLPObjective; break;
  }
  free_data_simplex(data);
}


template<typename T>
void dd_FindLPBasis(dd_rowrange m_size,dd_colrange d_size,
		    T** A, T** Ts,dd_rowindex OV,dd_rowset equalityset, dd_colindex nbindex,
		    dd_rowindex bflag,dd_rowrange objrow,dd_colrange rhscol,
		    dd_colrange *cs, bool *found,dd_LPStatusType *lps,long *pivot_no,
                    data_temp_simplex<T>* data)
{
  bool chosen,stop;
  long pivots_p0=0,rank;
  colset ColSelected;
  rowset RowSelected;
  T val;

  dd_rowrange r;
  dd_colrange j,s;

  *found=false; *cs=0; rank=0;
  stop=false;
  *lps=dd_LPSundecided;

  set_initialize(&RowSelected,m_size);
  set_initialize(&ColSelected,d_size);
  set_addelem(RowSelected,objrow);
  set_addelem(ColSelected,rhscol);

  stop=false;
  do {   /* Find a LP basis */
    dd_SelectPivot2(m_size,d_size,A,Ts,OV,equalityset,
		    m_size,RowSelected,ColSelected,&r,&s,&chosen);
    if (chosen) {
      set_addelem(RowSelected,r);
      set_addelem(ColSelected,s);
      //      std::cout << "dd_GaussianColumnPivot2 call 1\n";
      dd_GaussianColumnPivot2(d_size,A,Ts,nbindex,bflag,r,s, data->Rtemp);
      pivots_p0++;
      rank++;
    } else {
      for (j=1;j<=d_size  && *lps==dd_LPSundecided; j++) {
        if (j!=rhscol && nbindex[j]<0){
          dd_TableauEntry(val, d_size,A,Ts,objrow,j);
          if (dd_Nonzero(val)){  /* dual inconsistent */
            *lps=dd_StrucDualInconsistent;
            *cs=j;
            /* dual inconsistent because the nonzero reduced cost */
          }
        }
      }
      if (*lps==dd_LPSundecided) *found=true;
         /* dependent columns but not dual inconsistent. */
      stop=true;
    }
    /* printf("d_size=%ld, rank=%ld\n",d_size,rank); */
    if (rank==d_size-1) {
      stop = true;
      *found=true;
    }
  } while (!stop);

  *pivot_no=pivots_p0;
  set_free(RowSelected);
  set_free(ColSelected);
}


template<typename T>
void dd_FindLPBasis2(dd_rowrange m_size,dd_colrange d_size,
                     T** A, T** Ts,dd_rowindex OV,dd_rowset equalityset, dd_colindex nbindex,
                     dd_rowindex bflag,dd_rowrange objrow,dd_colrange rhscol,
		     dd_colrange *cs, bool *found,long *pivot_no, data_temp_simplex<T>* data)
{
  /* Similar to dd_FindLPBasis but it is much simpler.  This tries to recompute T for
  the specified basis given by nbindex.  It will return *found=false if the specified
  basis is not a basis.
  */
  int chosen,stop;
  long pivots_p0=0,rank;
  dd_colset ColSelected,DependentCols;
  dd_rowset RowSelected, NopivotRow;
  T val;
  bool localdebug=false;

  dd_rowrange r,negcount=0;
  dd_colrange j,s;

  *found=false; *cs=0; rank=0;

  set_initialize(&RowSelected,m_size);
  set_initialize(&DependentCols,d_size);
  set_initialize(&ColSelected,d_size);
  set_initialize(&NopivotRow,m_size);
  set_addelem(RowSelected,objrow);
  set_addelem(ColSelected,rhscol);
  set_compl(NopivotRow, NopivotRow);  /* set NopivotRow to be the groundset */

  for (j=2; j<=d_size; j++)
    if (nbindex[j]>0)
       set_delelem(NopivotRow, nbindex[j]);
    else if (nbindex[j]<0){
       negcount++;
       set_addelem(DependentCols, -nbindex[j]);
       set_addelem(ColSelected, -nbindex[j]);
    }

  set_uni(RowSelected, RowSelected, NopivotRow);  /* RowSelected is the set of rows not allowed to poviot on */

  stop=false;
  do {   /* Find a LP basis */
    dd_SelectPivot2(m_size,d_size,A,Ts,OV,equalityset, m_size,RowSelected,ColSelected,&r,&s,&chosen);
    if (chosen) {
      set_addelem(RowSelected,r);
      set_addelem(ColSelected,s);
      //      std::cout << "dd_GaussianColumnPivot2 call 2\n";
      dd_GaussianColumnPivot2(d_size,A,Ts,nbindex,bflag,r,s, data->Rtemp);
      pivots_p0++;
      rank++;
    } else{
      *found=false;   /* cannot pivot on any of the spacified positions. */
      stop=true;
    }
    if (rank==d_size-1-negcount) {
      if (negcount){
        /* Now it tries to pivot on rows that are supposed to be dependent. */
        set_diff(ColSelected, ColSelected, DependentCols);
        dd_SelectPivot2(m_size,d_size,A,Ts,OV,equalityset, m_size,RowSelected,ColSelected,&r,&s,&chosen);
        if (chosen) *found=false;  /* not supposed to be independent */
        else *found=true;
        if (localdebug){
          printf("Try to check the dependent cols:");
          if (chosen) printf("They are not dependent.  Can still pivot on (%ld, %ld)\n",r, s);
          else printf("They are indeed dependent.\n");
        }
      } else {
        *found=true;
     }
     stop = true;
    }
  } while (!stop);

  for (j=1; j<=d_size; j++) if (nbindex[j]>0) bflag[nbindex[j]]=j;
  *pivot_no=pivots_p0;
  set_free(RowSelected);
  set_free(ColSelected);
  set_free(NopivotRow);
  set_free(DependentCols);
}

template<typename T>
void dd_FindDualFeasibleBasis(dd_rowrange m_size,dd_colrange d_size,
                              T** A,T** Ts,
                              dd_colindex nbindex,dd_rowindex bflag,dd_rowrange objrow,
                              dd_colrange rhscol, bool lexicopivot,
                              dd_colrange *s,dd_ErrorType *err,dd_LPStatusType *lps,long *pivot_no,
                              long maxpivots, data_temp_simplex<T>* data)
{
  /* Find a dual feasible basis using Phase I of Dual Simplex method.
     If the problem is dual feasible,
     the procedure returns *err=NoError, *lps=LPSundecided and a dual feasible
     basis.   If the problem is dual infeasible, this returns
     *err=NoError, *lps=DualInconsistent and the evidence column *s.
     Caution: matrix A must have at least one extra row:  the row space A[m_size] must
     have been allocated.
  */
  bool phase1,dualfeasible=true;
  bool localdebug=false,chosen,stop;
  dd_LPStatusType LPSphase1;
  long pivots_p1=0;
  dd_rowrange i,r_val;
  dd_colrange j,l,ms=0,s_val,local_m_size;
  T x, val, maxcost=0, axvalue, maxratio=0, maxratio_q=1;

  T scaling,svalue;  /* random scaling mytype value */
  T minval=0;

  *err=dd_NoError; *lps=dd_LPSundecided; *s=0;
  local_m_size=m_size+1;  /* increase m_size by 1 */

  ms=0;  /* ms will be the index of column which has the largest reduced cost */
  for (j=1; j<=d_size; j++){
    if (j!=rhscol){
      dd_TableauEntry(data->rcost[j-1], d_size,A,Ts,objrow,j);
      if (dd_Larger(data->rcost[j-1],maxcost)) {dd_set(maxcost,data->rcost[j-1]); ms = j;}
    }
  }
  if (dd_Positive(maxcost)) dualfeasible=false;

  if (!dualfeasible){
    for (j=1; j<=d_size; j++){
      A[local_m_size-1][j-1]=0;
      for (l=1; l<=d_size; l++){
        if (nbindex[l]>0) {
          dd_set_si(scaling,l+10);
          dd_mul(svalue,A[nbindex[l]-1][j-1],scaling);
          dd_sub(A[local_m_size-1][j-1],A[local_m_size-1][j-1],svalue);
          /* To make the auxiliary row (0,-11,-12,...,-d-10).
             It is likely to be better than  (0, -1, -1, ..., -1)
             to avoid a degenerate LP.  Version 093c. */
        }
      }
    }

    ms=0;
     /* Ratio Test: ms will be now the index of column which has the largest reduced cost
        over the auxiliary row entry */
    for (j=1; j<=d_size; j++) {
      if ((j!=rhscol) && dd_Positive(data->rcost[j-1])){
        dd_TableauEntry(axvalue, d_size,A,Ts,local_m_size,j);
        if (dd_Nonnegative(axvalue)) {
          *err=dd_NumericallyInconsistent;
           /* This should not happen as they are set negative above.  Quit the phase I.*/
          goto _L99;
        }
        // So now axvalue < 0
        dd_neg(axvalue,axvalue);
        // So now axvalue > 0
        bool is_field=true;
        if (is_field) {
            dd_div(axvalue,data->rcost[j-1],axvalue);  /* axvalue is the negative of ratio that is to be maximized. */
            if (dd_Larger(axvalue, maxratio)) {
              dd_set(maxratio,axvalue);
              ms = j;
            }
          } else {
          if (dd_LargerFrac(data->rcost[j-1], axvalue, maxratio, maxratio_q)) {
            dd_set(maxratio,data->rcost[j-1]);
            dd_set(maxratio_q,axvalue);
            ms = j;
          }
        }
      }
    }

    if (ms==0) {
      *err=dd_NumericallyInconsistent; /* This should not happen. Quit the phase I.*/
      goto _L99;
    }

    /* Pivot on (local_m_size,ms) so that the dual basic solution becomes feasible */
    //    std::cout << "dd_GaussianColumnPivot2 call 3\n";
    dd_GaussianColumnPivot2(d_size,A,Ts,nbindex,bflag,local_m_size,ms, data->Rtemp);
    pivots_p1=pivots_p1+1;
    if (localdebug) {
      printf("\ndd_FindDualFeasibleBasis: Pivot on %ld %ld.\n",local_m_size,ms);
    }

  for (j=1; j<=d_size; j++) data->nbindex_ref[j]=nbindex[j];
     /* set the reference basis to be the current feasible basis. */

    phase1=true; stop=false;
    do {   /* Dual Simplex Phase I */
      chosen=false; LPSphase1=dd_LPSundecided;
      if (pivots_p1>maxpivots) {
        *err=dd_LPCycling;
        goto _L99;  /* failure due to max no. of pivots performed */
      }
      dd_SelectDualSimplexPivot(local_m_size, d_size, phase1, A, Ts, data->nbindex_ref, bflag,
				objrow, rhscol, lexicopivot, &r_val, &s_val, &chosen, &LPSphase1,
                                data);
      if (!chosen) {
        /* The current dictionary is terminal.  There are two cases:
           dd_TableauEntry(local_m_size,d_size,A,T,objrow,ms) is negative or zero.
           The first case implies dual infeasible,
           and the latter implies dual feasible but local_m_size is still in nonbasis.
           We must pivot in the auxiliary variable local_m_size.
        */
        dd_TableauEntry(x, d_size,A,Ts,objrow,ms);
        if (dd_Negative(x)){
          *err=dd_NoError; *lps=dd_DualInconsistent;  *s=ms;
        }

        r_val=0;
        for (i=1; i<=local_m_size; i++){
          if (bflag[i]<0) {
             /* i is basic and not the objective variable */
            dd_TableauEntry(val, d_size,A,Ts,i,ms);  /* auxiliary column*/
            if (dd_Smaller(val, minval)) {
              r_val=i;
              dd_set(minval,val);
            }
          }
        }

        if (r_val==0) {
          *err=dd_NumericallyInconsistent; /* This should not happen. Quit the phase I.*/
          goto _L99;
        }

        //        std::cout << "dd_GaussianColumnPivot2 call 4\n";
        dd_GaussianColumnPivot2(d_size,A,Ts,nbindex,bflag,r_val,ms, data->Rtemp);
        pivots_p1++;
        stop=true;
      } else {
        //        std::cout << "dd_GaussianColumnPivot2 call 5\n";
        dd_GaussianColumnPivot2(d_size,A,Ts,nbindex,bflag,r_val,s_val, data->Rtemp);
        pivots_p1=pivots_p1+1;
        if (bflag[local_m_size]<0) {
          stop=true;
        }
      }
    } while(!stop);
  }
_L99:
  *pivot_no=pivots_p1;
}

template<typename T>
void dd_DualSimplexMinimize(dd_lpdata<T> *lp, dd_ErrorType *err, data_temp_simplex<T>* data)
{
   dd_colrange j;
   *err=dd_NoError;
   for (j=1; j<=lp->d; j++)
     dd_neg(lp->A[lp->objrow-1][j-1], lp->A[lp->objrow-1][j-1]);
   dd_DualSimplexMaximize(lp, err, data);
   dd_neg(lp->optvalue, lp->optvalue);
   for (j=1; j<=lp->d; j++) {
     if (lp->LPS!=dd_Inconsistent)
       dd_neg(lp->dsol[j-1], lp->dsol[j-1]);
     dd_neg(lp->A[lp->objrow-1][j-1], lp->A[lp->objrow-1][j-1]);
   }
}

template<typename T>
void dd_DualSimplexMaximize(dd_lpdata<T> *lp,dd_ErrorType *err, data_temp_simplex<T>* data)
/*
When LP is inconsistent then lp->re returns the evidence row.
When LP is dual-inconsistent then lp->se returns the evidence column.
*/
{
  bool stop,chosen,phase1,found;
  long pivots_ds=0,pivots_p0=0,pivots_p1=0,pivots_pc=0,maxpivots,maxpivfactor=20;
  bool localdebug1=false;
  dd_rowrange i,r;
  dd_colrange j,s;
  dd_rowindex bflag;
  dd_rowindex OrderVector;  /* the permutation vector to store a preordered row indeces */
  dd_colindex nbindex_ref; /* to be used to store the initial feasible basis for lexico rule */

  unsigned int rseed=1;
  r=-40000; // values designed to create segfault in case it is not set later

  /* *err=dd_NoError; */
  set_emptyset(lp->redset_extra);
  for (i=0; i<= 4; i++) lp->pivots[i]=0;
  maxpivots=maxpivfactor*lp->d;  /* maximum pivots to be performed before cc pivot is applied. */
  OrderVector = new long[lp->m+1];
  for (int i=0; i<=lp->m; i++)
    OrderVector[i] = 0;
  bflag = new long[lp->m+2];
  for (int i=0; i<=lp->m+1; i++)
    bflag[i] = 0;
  nbindex_ref = new long[lp->d+1];
  for (int i=0; i<=lp->d; i++)
    nbindex_ref[i] = 0;
  /* Initializing control variables. */
  dd_ComputeRowOrderVector2(lp->m,lp->d,lp->A,OrderVector,dd_MinIndex,rseed);

  lp->re=0; lp->se=0;

  dd_ResetTableau(lp->m,lp->d,lp->B,lp->nbindex,bflag,lp->objrow,lp->rhscol);

  dd_FindLPBasis(lp->m,lp->d,lp->A,lp->B,OrderVector,lp->equalityset,lp->nbindex,bflag,
		 lp->objrow,lp->rhscol,&s,&found,&(lp->LPS),&pivots_p0, data);
  lp->pivots[0]=pivots_p0;

  if (!found){
     lp->se=s;
     goto _L99;
     /* No LP basis is found, and thus Inconsistent.
     Output the evidence column. */
  }

  dd_FindDualFeasibleBasis(lp->m,lp->d,lp->A,lp->B,lp->nbindex,bflag,
			   lp->objrow,lp->rhscol,lp->lexicopivot,&s, err,&(lp->LPS),&pivots_p1, maxpivots,
                           data);
  lp->pivots[1]=pivots_p1;

  for (j=1; j<=lp->d; j++) nbindex_ref[j]=lp->nbindex[j];
     /* set the reference basis to be the current feasible basis. */

  if (*err==dd_LPCycling || *err==dd_NumericallyInconsistent) {
    dd_CrissCrossMaximize(lp, err, data);
    delete [] OrderVector;
    delete [] bflag;
    delete [] nbindex_ref;
    return;
  }

  if (lp->LPS==dd_DualInconsistent){
     lp->se=s;
     goto _L99;
     /* No dual feasible basis is found, and thus DualInconsistent.
     Output the evidence column. */
  }

  /* Dual Simplex Method */
  stop=false;
  do {
    chosen=false; lp->LPS=dd_LPSundecided; phase1=false;
    if (pivots_ds<maxpivots) {
      dd_SelectDualSimplexPivot(lp->m,lp->d,phase1,lp->A,lp->B,nbindex_ref,bflag,
				lp->objrow,lp->rhscol,lp->lexicopivot,&r,&s,&chosen,&(lp->LPS),
                                data);
    }
    if (chosen) {
      pivots_ds=pivots_ds+1;
      if (lp->redcheck_extensive) {
        dd_GetRedundancyInformation(lp->m,lp->d,lp->A,lp->B, bflag, lp->redset_extra);
        set_uni(lp->redset_accum, lp->redset_accum,lp->redset_extra);
      }
    }
    if (!chosen && lp->LPS==dd_LPSundecided) {
      if (localdebug1){
         fprintf(stderr,"Warning: an emergency CC pivot in Phase II is performed\n");
         /* In principle this should not be executed because we already have dual feasibility
            attained and dual simplex pivot should have been chosen.  This might occur
            under floating point computation, or the case of cycling.
         */
      }


      dd_SelectCrissCrossPivot(lp->m,lp->d,lp->A,lp->B,bflag,
			       lp->objrow,lp->rhscol,&r,&s,&chosen,&(lp->LPS));
      if (chosen) pivots_pc=pivots_pc+1;
    }
    if (chosen) {
      //      std::cout << "dd_GaussianColumnPivot2 call 6\n";
      dd_GaussianColumnPivot2(lp->d,lp->A,lp->B,lp->nbindex,bflag,r,s, data->Rtemp);
    } else {
      switch (lp->LPS){
      case dd_Inconsistent: lp->re=r; break;
      case dd_DualInconsistent: lp->se=s; break;
      case dd_LPSundecided: break;
      case dd_Unbounded: break;
      case dd_DualUnbounded: break;
      case dd_StrucInconsistent: break;
      case dd_StrucDualInconsistent: break;
      case dd_Optimal: break;
      }
      stop=true;
    }
  } while(!stop);
_L99:
  lp->pivots[2]=pivots_ds;
  lp->pivots[3]=pivots_pc;
  dd_SetSolutions(lp->m,lp->d,lp->A,lp->B,lp->objrow,lp->rhscol,lp->LPS,lp->optvalue,lp->sol,lp->dsol,lp->posset_extra, lp->re,lp->se,bflag);
  delete [] OrderVector;
  delete [] bflag;
  delete [] nbindex_ref;
}



template<typename T>
void dd_CrissCrossMinimize(dd_lpdata<T> *lp,dd_ErrorType *err, data_temp_simplex<T>* data)
{
   dd_colrange j;

   *err=dd_NoError;
   for (j=1; j<=lp->d; j++)
     dd_neg(lp->A[lp->objrow-1][j-1],lp->A[lp->objrow-1][j-1]);
   dd_CrissCrossMaximize(lp, err, data);
   dd_neg(lp->optvalue, lp->optvalue);
   for (j=1; j<=lp->d; j++){
     if (lp->LPS!=dd_Inconsistent) {
	   /* Inconsistent certificate stays valid for minimization, 0.94e */
	   dd_neg(lp->dsol[j-1],lp->dsol[j-1]);
	 }
     dd_neg(lp->A[lp->objrow-1][j-1],lp->A[lp->objrow-1][j-1]);
   }
}

template<typename T>
void dd_CrissCrossMaximize(dd_lpdata<T> *lp,dd_ErrorType *err, data_temp_simplex<T>* data)
/*
When LP is inconsistent then lp->re returns the evidence row.
When LP is dual-inconsistent then lp->se returns the evidence column.
*/
{
  bool stop,chosen,found;
  long pivots0,pivots1;

  dd_rowrange i,r;
  dd_colrange s;
  dd_rowindex bflag;
  dd_rowindex OrderVector;  /* the permutation vector to store a preordered row indeces */
  unsigned int rseed=1;
  dd_colindex nbtemp;

  *err=dd_NoError;
  nbtemp = new long[lp->d+1];
  for (int i=0; i<=lp->d; i++)
    nbtemp[i] = 0;
  for (i=0; i<= 4; i++)
    lp->pivots[i]=0;
  bflag = new long[lp->m+1];
  for (int i=0; i<=lp->m; i++)
    bflag[i] = 0;
  OrderVector = new long[lp->m+1];
  for (int i=0; i<=lp->m; i++)
    OrderVector[i] = 0;
  /* Initializing control variables. */
  dd_ComputeRowOrderVector2(lp->m,lp->d,lp->A,OrderVector,dd_MinIndex,rseed);

  lp->re=0; lp->se=0; pivots1=0;

  dd_ResetTableau(lp->m,lp->d,lp->B,lp->nbindex,bflag,lp->objrow,lp->rhscol);

  dd_FindLPBasis(lp->m,lp->d,lp->A,lp->B,OrderVector,lp->equalityset,
		 lp->nbindex,bflag,lp->objrow,lp->rhscol,&s,&found,&(lp->LPS),&pivots0, data);
  lp->pivots[0]+=pivots0;

  if (!found){
     lp->se=s;
     goto _L99;
     /* No LP basis is found, and thus Inconsistent.
     Output the evidence column. */
  }

  stop=false;
  do {   /* Criss-Cross Method */

    dd_SelectCrissCrossPivot(lp->m,lp->d,lp->A,lp->B,bflag,
			     lp->objrow,lp->rhscol,&r,&s,&chosen,&(lp->LPS));
    if (chosen) {
      //      std::cout << "dd_GaussianColumnPivot2 call 7\n";
      dd_GaussianColumnPivot2(lp->d,lp->A,lp->B,lp->nbindex,bflag,r,s, data->Rtemp);
      pivots1++;
    } else {
      switch (lp->LPS){
      case dd_Inconsistent: lp->re=r; break;
      case dd_DualInconsistent: lp->se=s; break;
      case dd_Optimal: break;
      case dd_LPSundecided: break;
      case dd_Unbounded: break;
      case dd_DualUnbounded: break;
      case dd_StrucInconsistent: break;
      case dd_StrucDualInconsistent: break;
      }
      stop=true;
    }
  } while(!stop);

_L99:
  lp->pivots[1]+=pivots1;
  dd_SetSolutions(lp->m,lp->d,lp->A,lp->B,
		  lp->objrow,lp->rhscol,lp->LPS,lp->optvalue,lp->sol,lp->dsol,lp->posset_extra, lp->re,lp->se,bflag);
  delete [] nbtemp;
  delete [] bflag;   /* called previously with different lp->m */
  delete [] OrderVector;
}

template<typename T>
void dd_SetSolutions(dd_rowrange m_size,dd_colrange d_size,
   T** A,T** Ts,
   dd_rowrange objrow,dd_colrange rhscol,dd_LPStatusType LPS,
   T & optvalue,T* sol,T* dsol,dd_rowset posset,
		     dd_rowrange re,dd_colrange se,dd_rowindex bflag)
/*
Assign the solution vectors to sol,dsol,*optvalue after solving
the LP.
*/
{
  dd_rowrange i;
  dd_colrange j;
  T x,sw;
  int localdebug=false;

  if (localdebug) fprintf(stderr,"SetSolutions:\n");
  switch (LPS){
  case dd_Optimal:
    for (j=1;j<=d_size; j++) {
      dd_set(sol[j-1],Ts[j-1][rhscol-1]);
      dd_TableauEntry(x, d_size,A,Ts,objrow,j);
      dd_neg(dsol[j-1],x);
      dd_TableauEntry(optvalue, d_size,A,Ts,objrow,rhscol);
    }
    for (i=1; i<=m_size; i++) {
      if (bflag[i]==-1) {  /* i is a basic variable */
        dd_TableauEntry(x, d_size,A,Ts,i,rhscol);
        if (dd_Positive(x)) set_addelem(posset, i);
      }
    }

    break;
  case dd_Inconsistent:
    if (localdebug) fprintf(stderr,"SetSolutions: LP is inconsistent.\n");
    for (j=1;j<=d_size; j++) {
      dd_set(sol[j-1],Ts[j-1][rhscol-1]);
      dd_TableauEntry(x, d_size,A,Ts,re,j);
      dd_neg(dsol[j-1],x);
    }
    break;

  case dd_LPSundecided:
    std::cerr << "Case dd_LPSundecided has not been programmed in dd_SetSolutions\n";
    throw TerminalException{1};

  case dd_StrucInconsistent:
    std::cerr << "Case dd_StrucInconsistent has not been programmed in dd_SetSolutions\n";
    throw TerminalException{1};

  case dd_Unbounded:
    std::cerr << "Case dd_Unbounded has not been programmed in dd_SetSolutions\n";
    throw TerminalException{1};

  case dd_DualUnbounded:
    std::cerr << "Case dd_DualUnbounded has not been programmed in dd_SetSolutions\n";
    throw TerminalException{1};


  case dd_DualInconsistent:
    if (localdebug) printf( "SetSolutions: LP is dual inconsistent.\n");
    for (j=1;j<=d_size; j++) {
      dd_set(sol[j-1],Ts[j-1][se-1]);
      dd_TableauEntry(x, d_size,A,Ts,objrow,j);
      dd_neg(dsol[j-1],x);
    }
	break;

  case dd_StrucDualInconsistent:
    dd_TableauEntry(x, d_size,A,Ts,objrow,se);
    if (dd_Positive(x)) sw=1;
    else sw=-1;
    for (j=1;j<=d_size; j++) {
      dd_mul(sol[j-1],sw,Ts[j-1][se-1]);
      dd_TableauEntry(x, d_size,A,Ts,objrow,j);
      dd_neg(dsol[j-1],x);
    }
    if (localdebug) fprintf(stderr,"SetSolutions: LP is dual inconsistent.\n");
    break;

  }
}



template<typename T>
void dd_ComputeRowOrderVector2(dd_rowrange m_size,dd_colrange d_size,T** A,
			       dd_rowindex OV,dd_RowOrderType ho,unsigned int rseed)
{
  long i,itemp;

  OV[0]=0;
  switch (ho){
  case dd_MaxIndex:
    for(i=1; i<=m_size; i++) OV[i]=m_size-i+1;
    break;

  case dd_LexMin:
    for(i=1; i<=m_size; i++) OV[i]=i;
    dd_QuickSort(OV,1,m_size,A,d_size);
   break;

  case dd_LexMax:
    for(i=1; i<=m_size; i++) OV[i]=i;
    dd_QuickSort(OV,1,m_size,A,d_size);
    for(i=1; i<=m_size/2;i++){   /* just reverse the order */
      itemp=OV[i];
      OV[i]=OV[m_size-i+1];
      OV[m_size-i+1]=itemp;
    }
    break;

  case dd_RandomRow:
    for(i=1; i<=m_size; i++) OV[i]=i;
    if (rseed<=0) rseed=1;
    dd_RandomPermutation2(OV,m_size,rseed);
    break;

  case dd_MinIndex:
    for(i=1; i<=m_size; i++) OV[i]=i;
    break;

  case dd_MinCutoff: case dd_MaxCutoff: case dd_MixCutoff:
    for(i=1; i<=m_size; i++) OV[i]=i;
    break;
 }
}

template<typename T>
bool dd_LPSolve(dd_lpdata<T> *lp,dd_LPSolverType solver,dd_ErrorType *err)
/*
The current version of dd_LPSolve that solves an LP with floating-arithmetics first
and then with the specified arithimetics if it is GMP.

When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  int i;
  bool found=false;

  *err=dd_NoError;
  lp->solver=solver;

#ifndef GMPRATIONAL
  switch (lp->solver) {
    case dd_CrissCross:
      dd_CrissCrossSolve(lp,err);
      break;
    case dd_DualSimplex:
      dd_DualSimplexSolve(lp,err);
      break;
  }
#else
  lpf=dd_LPgmp2LPf(lp);
  switch (lp->solver) {
    case dd_CrissCross:
      ddf_CrissCrossSolve(lpf,&errf);    /* First, run with double float. */
      if (errf==ddf_NoError){   /* 094a:  fix for a bug reported by Dima Pasechnik */
        dd_BasisStatus(lpf,lp, &LPScorrect);    /* Check the basis. */
      } else {LPScorrect=false;}
      if (!LPScorrect) {
        if (localdebug) printf("BasisStatus: the current basis is NOT verified with GMP. Rerun with GMP.\n");
        dd_CrissCrossSolve(lp,err);  /* Rerun with GMP if fails. */
      } else {
        if (localdebug) printf("BasisStatus: the current basis is verified with GMP. The LP Solved.\n");
      }
      break;
    case dd_DualSimplex:
      ddf_DualSimplexSolve(lpf,&errf);    /* First, run with double float. */
      if (errf==ddf_NoError){   /* 094a:  fix for a bug reported by Dima Pasechnik */
        dd_BasisStatus(lpf,lp, &LPScorrect);    /* Check the basis. */
      } else {LPScorrect=false;}
      if (!LPScorrect){
	if (localdebug) printf("BasisStatus: the current basis is NOT verified with GMP. Rerun with GMP.\n");
	dd_DualSimplexSolve(lp,err);  /* Rerun with GMP if fails. */
	if (localdebug){
	  printf("*total number pivots = %ld (ph0 = %ld, ph1 = %ld, ph2 = %ld, ph3 = %ld, ph4 = %ld)\n",
		 lp->total_pivots,lp->pivots[0],lp->pivots[1],lp->pivots[2],lp->pivots[3],lp->pivots[4]);
	  ddf_WriteLPResult(stdout, lpf, errf);
	  dd_WriteLP(stdout, lp);
	}
      } else {
         if (localdebug) printf("BasisStatus: the current basis is verified with GMP. The LP Solved.\n");
      }
      break;
  }
  ddf_FreeLPData(lpf);
#endif

  lp->total_pivots=0;
  for (i=0; i<=4; i++) lp->total_pivots+=lp->pivots[i];
  if (*err==dd_NoError) found=true;
  return found;
}


template<typename T>
bool dd_LPSolve0(dd_lpdata<T> *lp,dd_LPSolverType solver,dd_ErrorType *err)
/*
The original version of dd_LPSolve that solves an LP with specified arithimetics.

When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  int i;
  bool found=false;

  *err=dd_NoError;
  lp->solver=solver;

  switch (lp->solver) {
    case dd_CrissCross:
      dd_CrissCrossSolve(lp,err);
      break;
    case dd_DualSimplex:
      dd_DualSimplexSolve(lp,err);
      break;
  }

  lp->total_pivots=0;
  for (i=0; i<=4; i++) lp->total_pivots+=lp->pivots[i];
  if (*err==dd_NoError) found=true;
  return found;
}


template<typename T>
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
  T bm,bmax,bceil;

  bm=2;
  bmax=1;
  m=lp->m+1;
  d=lp->d+1;

  lpnew=dd_CreateLPData<T>(m, d);
  lpnew->objective = dd_LPmax;

  for (i=1; i<=lp->m; i++) {
    if (dd_Larger(lp->A[i-1][lp->rhscol-1],bmax))
      dd_set(bmax,lp->A[i-1][lp->rhscol-1]);
  }
  dd_mul(bceil,bm,bmax);

  for (i=1; i <= lp->m; i++) {
    for (j=1; j <= lp->d; j++) {
      dd_set(lpnew->A[i-1][j-1],lp->A[i-1][j-1]);
    }
  }

  for (i=1;i<=lp->m; i++){
    lpnew->A[i-1][lp->d]=-1;
  }

  for (j=1;j<=lp->d;j++){
    lpnew->A[m-2][j-1]=0;   /* new row (bceil, 0,...,0,-1) */
  }
  dd_set(lpnew->A[m-2][0],bceil);  /* new row (bceil, 0,...,0,-1) */

  for (j=1;j<= d-1;j++) {
    lpnew->A[m-1][j-1]=0;  /* new obj row with (0,...,0,1) */
  }
  lpnew->A[m-1][d-1]=1;

  return lpnew;
}


template<typename T>
dd_lpdata<T> *dd_CreateLP_H_ImplicitLinearity(dd_matrixdata<T> *M)
{
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_lpdata<T> *lp;
  bool localdebug=false;

  linc=set_card(M->linset);
  m=M->rowsize+1+linc+1;
     /* We represent each equation by two inequalities.
        This is not the best way but makes the code simple. */
  d=M->colsize+1;

  lp=dd_CreateLPData<T>(m, d);
  lp->Homogeneous = true;
  lp->objective = dd_LPmax;
  lp->eqnumber=linc;  /* this records the number of equations */
  lp->redcheck_extensive=false;  /* this is default */

  irev=M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    if (set_member(i, M->linset)) {
      irev=irev+1;
      set_addelem(lp->equalityset,i);    /* it is equality. */
            /* the reversed row irev is not in the equality set. */
      for (j = 1; j <= M->colsize; j++) {
        dd_neg(lp->A[irev-1][j-1],M->matrix[i-1][j-1]);
      }  /*of j*/
    } else {
      lp->A[i-1][d-1]=-1;  /* b_I + A_I x - 1 z >= 0  (z=x_d) */
    }
    for (j = 1; j <= M->colsize; j++) {
      dd_set(lp->A[i-1][j-1],M->matrix[i-1][j-1]);
      if (j==1 && i<M->rowsize && dd_Nonzero(M->matrix[i-1][j-1])) lp->Homogeneous = false;
    }  /*of j*/
  }  /*of i*/
  lp->A[m-2][0]=1;
  lp->A[m-2][d-1]=-1;
      /* make the LP bounded.  */

  lp->A[m-1][d-1]=1;
      /* objective is to maximize z.  */

  if (localdebug) {
    fprintf(stderr,"dd_CreateLP_H_ImplicitLinearity: an new lp is\n");
    dd_WriteLP(stderr,lp);
  }

  return lp;
}

template<typename T>
dd_lpdata<T> *dd_CreateLP_V_ImplicitLinearity(dd_matrixdata<T> *M)
{
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_lpdata<T> *lp;
  bool localdebug=false;

  linc=set_card(M->linset);
  m=M->rowsize+1+linc+1;
     /* We represent each equation by two inequalities.
        This is not the best way but makes the code simple. */
  d=(M->colsize)+2;
     /* Two more columns.  This is different from the H-reprentation case */

/* The below must be modified for V-representation!!!  */

  lp=dd_CreateLPData<T>(m, d);
  lp->Homogeneous = false;
  lp->objective = dd_LPmax;
  lp->eqnumber = linc;  /* this records the number of equations */
  lp->redcheck_extensive = false;  /* this is default */

  irev=M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    lp->A[i-1][0]=0;  /* It is almost completely degerate LP */
    if (set_member(i, M->linset)) {
      irev=irev+1;
      set_addelem(lp->equalityset,i);    /* it is equality. */
            /* the reversed row irev is not in the equality set. */
      for (j = 2; j <= d-1; j++) {
        dd_neg(lp->A[irev-1][j-1],M->matrix[i-1][j-2]);
      }  /*of j*/
      if (localdebug) fprintf(stderr,"equality row %ld generates the reverse row %ld.\n",i,irev);
    } else {
      lp->A[i-1][d-1]=-1;
    }
    for (j = 2; j <= d-1; j++) {
      dd_set(lp->A[i-1][j-1],M->matrix[i-1][j-2]);
    }  /*of j*/
  }  /*of i*/
  lp->A[m-2][0]=1;
  lp->A[m-2][d-1]=-1;
  /* make the LP bounded.  */
  lp->A[m-1][d-1]=1;
  /* objective is to maximize z.  */

  if (localdebug) {
    fprintf(stderr,"dd_CreateLP_V_ImplicitLinearity: an new lp is\n");
    dd_WriteLP(stderr,lp);
  }

  return lp;
}


template<typename T>
dd_lpdata<T> *dd_CreateLP_H_Redundancy(dd_matrixdata<T> *M, dd_rowrange itest)
{
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_lpdata<T> *lp;

  linc=set_card(M->linset);
  m=M->rowsize+1+linc;
     /* We represent each equation by two inequalities.
        This is not the best way but makes the code simple. */
  d=M->colsize;

  lp=dd_CreateLPData_from_M<T>(M);
  lp->Homogeneous = true;
  lp->objective = dd_LPmin;
  lp->eqnumber=linc;  /* this records the number of equations */
  lp->redcheck_extensive=false;  /* this is default */

  irev=M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    if (set_member(i, M->linset)) {
      irev=irev+1;
      set_addelem(lp->equalityset,i);    /* it is equality. */
            /* the reversed row irev is not in the equality set. */
      for (j = 1; j <= M->colsize; j++) {
        dd_neg(lp->A[irev-1][j-1],M->matrix[i-1][j-1]);
      }  /*of j*/
    }
    for (j = 1; j <= d; j++) {
      dd_set(lp->A[i-1][j-1],M->matrix[i-1][j-1]);
      if (j==1 && i<M->rowsize && dd_Nonzero(M->matrix[i-1][j-1])) lp->Homogeneous = false;
    }  /*of j*/
  }  /*of i*/
  for (j = 1; j <= d; j++) {
    dd_set(lp->A[m-1][j-1],M->matrix[itest-1][j-1]);
      /* objective is to violate the inequality in question.  */
  }  /*of j*/
  dd_add<T>(lp->A[itest-1][0],lp->A[itest-1][0],1); /* relax the original inequality by one */

  return lp;
}


template<typename T>
dd_lpdata<T> *dd_CreateLP_V_Redundancy(dd_matrixdata<T> *M, dd_rowrange itest)
{
  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_lpdata<T> *lp;

  linc=set_card(M->linset);
  m=M->rowsize+1+linc;
     /* We represent each equation by two inequalities.
        This is not the best way but makes the code simple. */
  d=(M->colsize)+1;
     /* One more column.  This is different from the H-reprentation case */

/* The below must be modified for V-representation!!!  */

  lp=dd_CreateLPData_from_M<T>(M);
  lp->Homogeneous = false;
  lp->objective = dd_LPmin;
  lp->eqnumber=linc;  /* this records the number of equations */
  lp->redcheck_extensive=false;  /* this is default */

  irev=M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    if (i==itest){
      lp->A[i-1][0]=1; /* this is to make the LP bounded, ie. the min >= -1 */
    } else {
      lp->A[i-1][0]=0;  /* It is almost completely degerate LP */
    }
    if (set_member(i, M->linset)) {
      irev=irev+1;
      set_addelem(lp->equalityset,i);    /* it is equality. */
            /* the reversed row irev is not in the equality set. */
      for (j = 2; j <= (M->colsize)+1; j++) {
        dd_neg(lp->A[irev-1][j-1],M->matrix[i-1][j-2]);
      }  /*of j*/
    }
    for (j = 2; j <= d; j++) {
      dd_set(lp->A[i-1][j-1],M->matrix[i-1][j-2]);
    }  /*of j*/
  }  /*of i*/
  for (j = 2; j <= d; j++) {
    dd_set(lp->A[m-1][j-1],M->matrix[itest-1][j-2]);
      /* objective is to violate the inequality in question.  */
  }  /*of j*/
  lp->A[m-1][0]=0;   /* the constant term for the objective is zero */

  return lp;
}


template<typename T>
dd_lpdata<T> *dd_CreateLP_V_SRedundancy(dd_matrixdata<T> *M, dd_rowrange itest)
{
/*
     V-representation (=boundary problem)
       g* = maximize
         1^T b_{I-itest} x_0 + 1^T A_{I-itest}    (the sum of slacks)
       subject to
         b_itest x_0     + A_itest x      =  0 (the point has to lie on the boundary)
         b_{I-itest} x_0 + A_{I-itest} x >=  0 (all nonlinearity generators in one side)
         1^T b_{I-itest} x_0 + 1^T A_{I-itest} x <=  1 (to make an LP bounded)
         b_L x_0         + A_L x = 0.  (linearity generators)

    The redundant row is strongly redundant if and only if g* is zero.
*/

  dd_rowrange m, i, irev, linc;
  dd_colrange d, j;
  dd_lpdata<T> *lp;

  linc=set_card(M->linset);
  m=M->rowsize+1+linc+2;
     /* We represent each equation by two inequalities.
        This is not the best way but makes the code simple.
        Two extra constraints are for the first equation and the bouding inequality.
        */
  d=(M->colsize)+1;
     /* One more column.  This is different from the H-reprentation case */

/* The below must be modified for V-representation!!!  */

  lp=dd_CreateLPData<T>(m, d);
  lp->Homogeneous = false;
  lp->objective = dd_LPmax;
  lp->eqnumber=linc;  /* this records the number of equations */

  irev=M->rowsize; /* the first row of the linc reversed inequalities. */
  for (i = 1; i <= M->rowsize; i++) {
    if (i==itest){
      lp->A[i-1][0]=0;  /* this is a half of the boundary constraint. */
    } else {
      lp->A[i-1][0]=0;  /* It is almost completely degerate LP */
    }
    if (set_member(i, M->linset) || i==itest) {
      irev=irev+1;
      set_addelem(lp->equalityset,i);    /* it is equality. */
            /* the reversed row irev is not in the equality set. */
      for (j = 2; j <= (M->colsize)+1; j++) {
        dd_neg(lp->A[irev-1][j-1],M->matrix[i-1][j-2]);
      }  /*of j*/
    }
    for (j = 2; j <= (M->colsize)+1; j++) {
      dd_set(lp->A[i-1][j-1],M->matrix[i-1][j-2]);
      dd_add(lp->A[m-1][j-1],lp->A[m-1][j-1],lp->A[i-1][j-1]);  /* the objective is the sum of all ineqalities */
    }  /*of j*/
  }  /*of i*/
  for (j = 2; j <= (M->colsize)+1; j++) {
    dd_neg(lp->A[m-2][j-1],lp->A[m-1][j-1]);
      /* to make an LP bounded.  */
  }  /*of j*/
  lp->A[m-2][0]=1;   /* the constant term for the bounding constraint is 1 */

  return lp;
}

template<typename T>
bool dd_Redundant(dd_matrixdata<T> *M, dd_rowrange itest, T* certificate, dd_ErrorType *error)
  /* 092 */
{
  /* Checks whether the row itest is redundant for the representation.
     All linearity rows are not checked and considered NONredundant.
     This code works for both H- and V-representations.  A certificate is
     given in the case of non-redundancy, showing a solution x violating only the itest
     inequality for H-representation, a hyperplane RHS and normal (x_0, x) that
     separates the itest from the rest.  More explicitly, the LP to be setup is

     H-representation
       f* = minimize
         b_itest     + A_itest x
       subject to
         b_itest + 1 + A_itest x     >= 0 (relaxed inequality to make an LP bounded)
         b_{I-itest} + A_{I-itest} x >= 0 (all inequalities except for itest)
         b_L         + A_L x = 0.  (linearity)

     V-representation (=separation problem)
       f* = minimize
         b_itest x_0     + A_itest x
       subject to
         b_itest x_0     + A_itest x     >= -1 (to make an LP bounded)
         b_{I-itest} x_0 + A_{I-itest} x >=  0 (all nonlinearity generators except for itest in one side)
         b_L x_0         + A_L x = 0.  (linearity generators)

    Here, the input matrix is considered as (b, A), i.e. b corresponds to the first column of input
    and the row indices of input is partitioned into I and L where L is the set of linearity.
    In both cases, the itest data is nonredundant if and only if the optimal value f* is negative.
    The certificate has dimension one more for V-representation case.
  */

  dd_colrange j;
  dd_lpdata<T> *lp;
  dd_ErrorType err=dd_NoError;
  bool answer=false,localdebug=false;

  *error=dd_NoError;
  if (set_member(itest, M->linset)){
    if (localdebug) printf("The %ld th row is linearity and redundancy checking is skipped.\n",itest);
    goto _L99;
  }

  /* Create an LP data for redundancy checking */
  if (M->representation==dd_Generator){
    lp=dd_CreateLP_V_Redundancy(M, itest);
  } else {
    lp=dd_CreateLP_H_Redundancy(M, itest);
  }

  dd_LPSolve(lp,dd_choiceRedcheckAlgorithm,&err);
  if (err!=dd_NoError){
    *error=err;
    goto _L999;
  } else {
    for (j=0; j<lp->d; j++) {
      dd_set(certificate[j], lp->sol[j]);
    }

    if (dd_Negative(lp->optvalue)){
      answer=false;
      if (localdebug) fprintf(stderr,"==> %ld th row is nonredundant.\n",itest);
    } else {
      answer=true;
      if (localdebug) fprintf(stderr,"==> %ld th row is redundant.\n",itest);
    }
  }
  _L999:
  dd_FreeLPData(lp);
_L99:
  return answer;
}

template<typename T>
bool dd_RedundantExtensive(dd_matrixdata<T> *M, dd_rowrange itest, T* certificate,
                           dd_rowset *redset,dd_ErrorType *error)
{
  /* This uses the same LP construction as dd_Redundant.  But, while it is checking
     the redundancy of itest, it also tries to find some other variable that are
     redundant (i.e. forced to be nonnegative).  This is expensive as it used
     the complete tableau information at each DualSimplex pivot.  The redset must
     be initialized before this function is called.
  */

  dd_colrange j;
  dd_lpdata<T> *lp;
  dd_ErrorType err=dd_NoError;
  bool answer=false,localdebug=false;

  *error=dd_NoError;
  if (set_member(itest, M->linset)){
    if (localdebug) printf("The %ld th row is linearity and redundancy checking is skipped.\n",itest);
    goto _L99;
  }

  /* Create an LP data for redundancy checking */
  if (M->representation==dd_Generator){
    lp=dd_CreateLP_V_Redundancy(M, itest);
  } else {
    lp=dd_CreateLP_H_Redundancy(M, itest);
  }

  lp->redcheck_extensive=true;

  dd_LPSolve0(lp,dd_DualSimplex,&err);
  if (err!=dd_NoError){
    *error=err;
    goto _L999;
  } else {
    set_copy(*redset,lp->redset_extra);
    set_delelem(*redset, itest);
    /* itest row might be redundant in the lp but this has nothing to do with its redundancy
    in the original system M.   Thus we must delete it.  */
    for (j=0; j<lp->d; j++) {
      dd_set(certificate[j], lp->sol[j]);
    }

    if (dd_Negative(lp->optvalue)){
      answer=false;
      if (localdebug) fprintf(stderr,"==> %ld th row is nonredundant.\n",itest);
    } else {
      answer=true;
      if (localdebug) fprintf(stderr,"==> %ld th row is redundant.\n",itest);
    }
  }
  _L999:
  dd_FreeLPData(lp);
_L99:
  return answer;
}

template<typename T>
dd_rowset dd_RedundantRows(dd_matrixdata<T> *M, dd_ErrorType *error)
{
  dd_rowrange i,m;
  dd_colrange d;
  dd_rowset redset;
  dd_matrixdata<T> *Mcopy;
  T* cvec; /* certificate */
  bool localdebug=false;

  m=M->rowsize;
  if (M->representation==dd_Generator){
    d=(M->colsize)+1;
  } else {
    d=M->colsize;
  }
  Mcopy=dd_MatrixCopy(M);
  dd_InitializeArow(d,&cvec);
  set_initialize(&redset, m);
  for (i=m; i>=1; i--) {
    if (dd_Redundant(Mcopy, i, cvec, error)) {
      if (localdebug) printf("dd_RedundantRows: the row %ld is redundant.\n", i);
      set_addelem(redset, i);
      dd_MatrixRowRemove(&Mcopy, i);
    } else {
      if (localdebug) printf("dd_RedundantRows: the row %ld is essential.\n", i);
    }
    if (*error!=dd_NoError) goto _L99;
  }
_L99:
  dd_FreeMatrix(Mcopy);
  dd_FreeArow(cvec);
  return redset;
}


template<typename T>
bool dd_MatrixRedundancyRemove(dd_matrixdata<T> **M, dd_rowset *redset,dd_rowindex *newpos, dd_ErrorType *error)
{
  /* It returns the set of all redundant rows.  This should be called after all
     implicit linearity are recognized with dd_MatrixCanonicalizeLinearity.
  */


  dd_rowrange i,k,m,m1;
  dd_colrange d;
  dd_rowset redset1;
  dd_rowindex newpos1;
  dd_matrixdata<T> *M1=nullptr;
  T* cvec; /* certificate */
  bool success=false, localdebug=false;

  m=(*M)->rowsize;
  set_initialize(redset, m);
  M1=dd_MatrixSortedUniqueCopy(*M,newpos);
  for (i=1; i<=m; i++){
    if ((*newpos)[i]<=0) set_addelem(*redset,i);
    if (localdebug) printf(" %ld:%ld",i,(*newpos)[i]);
  }
  if (localdebug) printf("\n");

  if ((*M)->representation==dd_Generator){
    d=((*M)->colsize)+1;
  } else {
    d=(*M)->colsize;
  }
  m1=M1->rowsize;
  if (localdebug){
    fprintf(stderr,"dd_MatrixRedundancyRemove: By sorting, %ld rows have been removed.  The remaining has %ld rows.\n",m-m1,m1);
    /* dd_WriteMatrix(stdout,M1);  */
  }
  dd_InitializeArow(d,&cvec);
  set_initialize(&redset1, M1->rowsize);
  k=1;
  do {
    if (dd_RedundantExtensive(M1, k, cvec, &redset1,error)) {
      set_addelem(redset1, k);
      dd_MatrixRowsRemove2(&M1,redset1,&newpos1);
      for (i=1; i<=m; i++){
        if ((*newpos)[i]>0){
          if  (set_member((*newpos)[i],redset1)){
            set_addelem(*redset,i);
            (*newpos)[i]=0;  /* now the original row i is recognized redundant and removed from M1 */
          } else {
            (*newpos)[i]=newpos1[(*newpos)[i]];  /* update the new pos vector */
          }
        }
      }
      set_free(redset1);
      set_initialize(&redset1, M1->rowsize);
      if (localdebug) {
        printf("dd_MatrixRedundancyRemove: the row %ld is redundant. The new matrix has %ld rows.\n", k, M1->rowsize);
        /* dd_WriteMatrix(stderr, M1);  */
      }
      delete [] newpos1;
    } else {
      if (set_card(redset1)>0) {
        dd_MatrixRowsRemove2(&M1,redset1,&newpos1);
        for (i=1; i<=m; i++){
          if ((*newpos)[i]>0){
            if  (set_member((*newpos)[i],redset1)){
              set_addelem(*redset,i);
              (*newpos)[i]=0;  /* now the original row i is recognized redundant and removed from M1 */
            } else {
              (*newpos)[i]=newpos1[(*newpos)[i]];  /* update the new pos vector */
            }
          }
        }
        set_free(redset1);
        set_initialize(&redset1, M1->rowsize);
        delete [] newpos1;
      }
      if (localdebug) {
        printf("dd_MatrixRedundancyRemove: the row %ld is essential. The new matrix has %ld rows.\n", k, M1->rowsize);
        /* dd_WriteMatrix(stderr, M1);  */
      }
      k=k+1;
    }
    if (*error!=dd_NoError) goto _L99;
  } while  (k<=M1->rowsize);
  if (localdebug) dd_WriteMatrix(stderr, M1);
  success=true;

_L99:
  dd_FreeMatrix(*M);
  *M=M1;
  dd_FreeArow(cvec);
  set_free(redset1);
  return success;
}


template<typename T>
bool dd_SRedundant(dd_matrixdata<T> *M, dd_rowrange itest, T* certificate, dd_ErrorType *error)
{
  /* Checks whether the row itest is strongly redundant for the representation.
     A row is strongly redundant in H-representation if every point in
     the polyhedron satisfies it with strict inequality.
     A row is strongly redundant in V-representation if this point is in
     the interior of the polyhedron.

     All linearity rows are not checked and considered NOT strongly redundant.
     This code works for both H- and V-representations.  A certificate is
     given in the case of non-redundancy, showing a solution x violating only the itest
     inequality for H-representation, a hyperplane RHS and normal (x_0, x) that
     separates the itest from the rest.  More explicitly, the LP to be setup is

     H-representation
       f* = minimize
         b_itest     + A_itest x
       subject to
         b_itest + 1 + A_itest x     >= 0 (relaxed inequality to make an LP bounded)
         b_{I-itest} + A_{I-itest} x >= 0 (all inequalities except for itest)
         b_L         + A_L x = 0.  (linearity)

     V-representation (=separation problem)
       f* = minimize
         b_itest x_0     + A_itest x
       subject to
         b_itest x_0     + A_itest x     >= -1 (to make an LP bounded)
         b_{I-itest} x_0 + A_{I-itest} x >=  0 (all nonlinearity generators except for itest in one side)
         b_L x_0         + A_L x = 0.  (linearity generators)

    Here, the input matrix is considered as (b, A), i.e. b corresponds to the first column of input
    and the row indices of input is partitioned into I and L where L is the set of linearity.
    In H-representation, the itest data is strongly redundant if and only if the optimal value f* is positive.
    In V-representation, the itest data is redundant if and only if the optimal value f* is zero (as the LP
    is homogeneous and the optimal value is always non-positive).  To recognize strong redundancy, one
    can set up a second LP

     V-representation (=boundary problem)
       g* = maximize
         1^T b_{I-itest} x_0 + 1^T A_{I-itest}    (the sum of slacks)
       subject to
         b_itest x_0     + A_itest x      =  0 (the point has to lie on the boundary)
         b_{I-itest} x_0 + A_{I-itest} x >=  0 (all nonlinearity generators in one side)
         1^T b_{I-itest} x_0 + 1^T A_{I-itest} x <=  1 (to make an LP bounded)
         b_L x_0         + A_L x = 0.  (linearity generators)

    The redundant row is strongly redundant if and only if g* is zero.

    The certificate has dimension one more for V-representation case.
  */

  dd_colrange j;
  dd_lpdata<T> *lp;
  dd_ErrorType err=dd_NoError;
  bool answer=false,localdebug=false;

  *error=dd_NoError;
  if (set_member(itest, M->linset)){
    if (localdebug) printf("The %ld th row is linearity and strong redundancy checking is skipped.\n",itest);
    goto _L99;
  }

  /* Create an LP data for redundancy checking */
  if (M->representation==dd_Generator){
    lp=dd_CreateLP_V_Redundancy(M, itest);
  } else {
    lp=dd_CreateLP_H_Redundancy(M, itest);
  }

  dd_LPSolve(lp,dd_choiceRedcheckAlgorithm,&err);
  if (err!=dd_NoError){
    *error=err;
    goto _L999;
  } else {

    for (j=0; j<lp->d; j++) {
      dd_set(certificate[j], lp->sol[j]);
    }

    if (M->representation==dd_Inequality){
      if (dd_Positive(lp->optvalue)){
          answer=true;
          if (localdebug) fprintf(stderr,"==> %ld th inequality is strongly redundant.\n",itest);
        } else {
          answer=false;
          if (localdebug) fprintf(stderr,"==> %ld th inequality is not strongly redundant.\n",itest);
        }
    } else {
      if (dd_Negative(lp->optvalue)){
          answer=false;
          if (localdebug) fprintf(stderr,"==> %ld th point is not strongly redundant.\n",itest);
        } else {
          /* for V-representation, we have to solve another LP */
          dd_FreeLPData(lp);
          lp=dd_CreateLP_V_SRedundancy(M, itest);
          dd_LPSolve(lp,dd_DualSimplex,&err);
          if (localdebug) dd_WriteLPResult(stdout,lp,err);
          if (dd_Positive(lp->optvalue)){
            answer=false;
            if (localdebug) fprintf(stderr,"==> %ld th point is not strongly redundant.\n",itest);
          } else {
            answer=true;
            if (localdebug) fprintf(stderr,"==> %ld th point is strongly redundant.\n",itest);
          }
       }
    }
  }
  _L999:
  dd_FreeLPData(lp);
_L99:
  return answer;
}

template<typename T>
dd_rowset dd_SRedundantRows(dd_matrixdata<T> *M, dd_ErrorType *error)
{
  dd_rowrange i,m;
  dd_colrange d;
  dd_rowset redset;
  dd_matrixdata<T> *Mcopy;
  T* cvec; /* certificate */
  bool localdebug=false;

  m=M->rowsize;
  if (M->representation==dd_Generator){
    d=(M->colsize)+1;
  } else {
    d=M->colsize;
  }
  Mcopy=dd_MatrixCopy(M);
  dd_InitializeArow(d,&cvec);
  set_initialize(&redset, m);
  for (i=m; i>=1; i--) {
    if (dd_SRedundant(Mcopy, i, cvec, error)) {
      if (localdebug) printf("dd_SRedundantRows: the row %ld is strongly redundant.\n", i);
      set_addelem(redset, i);
      dd_MatrixRowRemove(&Mcopy, i);
    } else {
      if (localdebug) printf("dd_SRedundantRows: the row %ld is not strongly redundant.\n", i);
    }
    if (*error!=dd_NoError) goto _L99;
  }
_L99:
  dd_FreeMatrix(Mcopy);
  dd_FreeArow(cvec);
  return redset;
}

template<typename T>
dd_rowset dd_RedundantRowsViaShooting(dd_matrixdata<T> *M, dd_ErrorType *error)
{
  /*
     For H-representation only and not quite reliable,
     especially when floating-point arithmetic is used.
     Use the ordinary (slower) method dd_RedundantRows.
  */

  dd_rowrange i,m, ired, irow=0;
  dd_colrange j,k,d;
  dd_rowset redset;
  dd_rowindex rowflag;
    /* ith comp is negative if the ith inequality (i-1 st row) is redundant.
                   zero     if it is not decided.
                   k > 0    if it is nonredundant and assigned to the (k-1)th row of M1.
    */
  dd_matrixdata<T> *M1;
  T* shootdir;
  T* cvec=nullptr;
  dd_lpdata<T>* lp0;
  dd_lpdata<T>* lp;
  dd_ErrorType err;
  dd_LPSolverType solver=dd_DualSimplex;
  bool localdebug=false;

  m=M->rowsize;
  d=M->colsize;
  M1=dd_CreateMatrix<T>(m,d);
  M1->rowsize=0;  /* cheat the rowsize so that smaller matrix can be stored */
  set_initialize(&redset, m);
  dd_InitializeArow(d, &shootdir);
  dd_InitializeArow(d, &cvec);

  rowflag = new long[m+1];
  for (int i=0; i<=m; i++)
    rowflag[i] = 0;

  /* First find some (likely) nonredundant inequalities by Interior Point Find. */
  lp0=dd_Matrix2LP(M, &err);
  lp=dd_MakeLPforInteriorFinding(lp0);
  dd_FreeLPData(lp0);
  dd_LPSolve(lp, solver, &err);  /* Solve the LP */

  if (dd_Positive(lp->optvalue)){
    if (localdebug) std::cout << "dd_Positive=T case\n";
    /* An interior point is found.  Use rayshooting to find some nonredundant
       inequalities. */
    for (j=1; j<d; j++){
      for (k=1; k<=d; k++) shootdir[k-1]=0;
      shootdir[j]=1;  /* j-th unit vector */
      if (localdebug) dd_WriteT(std::cerr, shootdir, d);
      ired=dd_RayShooting(M, lp->sol, shootdir);
      if (localdebug) printf("nonredundant row %3ld found by shooting.\n", ired);
      if (ired>0 && rowflag[ired]<=0) {
        irow++;
        rowflag[ired]=irow;
        for (k=1; k<=d; k++) dd_set(M1->matrix[irow-1][k-1], M->matrix[ired-1][k-1]);
      }

      shootdir[j]=-1;  /* negative of the j-th unit vector */
      ired=dd_RayShooting(M, lp->sol, shootdir);
      if (localdebug) printf("nonredundant row %3ld found by shooting.\n", ired);
      if (ired>0 && rowflag[ired]<=0) {
        irow++;
        rowflag[ired]=irow;
        for (k=1; k<=d; k++) dd_set(M1->matrix[irow-1][k-1], M->matrix[ired-1][k-1]);
      }
    }

    M1->rowsize=irow;
    if (localdebug) {
      printf("The initial nonredundant set is:");
      for (i=1; i<=m; i++) if (rowflag[i]>0) printf(" %ld", i);
      printf("\n");
    }

    i=1;
    while(i<=m){
      if (localdebug) std::cout << "i=" << i << " rowflag=" << rowflag[i] << "\n";
      if (rowflag[i]==0){ /* the ith inequality is not yet checked */
        if (localdebug) fprintf(stderr, "Checking redundancy of %ld th inequality\n", i);
        irow++;  M1->rowsize=irow;
        for (k=1; k<=d; k++) dd_set(M1->matrix[irow-1][k-1], M->matrix[i-1][k-1]);
        if (!dd_Redundant(M1, irow, cvec, &err)){
          for (k=1; k<=d; k++) dd_sub(shootdir[k-1], cvec[k-1], lp->sol[k-1]);
          ired=dd_RayShooting(M, lp->sol, shootdir);
          rowflag[ired]=irow;
          for (k=1; k<=d; k++) dd_set(M1->matrix[irow-1][k-1], M->matrix[ired-1][k-1]);
          if (localdebug) {
            fprintf(stderr, "The %ld th inequality is nonredundant for the subsystem\n", i);
            fprintf(stderr, "The nonredundancy of %ld th inequality is found by shooting.\n", ired);
          }
        } else {
          if (localdebug) fprintf(stderr, "The %ld th inequality is redundant for the subsystem and thus for the whole.\n", i);
          rowflag[i]=-1;
          set_addelem(redset, i);
          i++;
        }
      } else {
        if (localdebug) std::cout << "Case rowflag != 0\n";
        i++;
      }
    } /* endwhile */
  } else {
    if (localdebug) std::cout << "dd_Positive=F case\n";
    /* No interior point is found.  Apply the standard LP technique.  */
    set_free(redset);
    redset=dd_RedundantRows(M, error);
  }

  dd_FreeLPData(lp);

  M1->rowsize=m; M1->colsize=d;  /* recover the original sizes. Needed for the deallocation */
  dd_FreeMatrix(M1);
  dd_FreeArow(shootdir);
  dd_FreeArow(cvec);
  delete [] rowflag;
  return redset;
}

template<typename T>
dd_SetFamilyPtr dd_Matrix2Adjacency(dd_matrixdata<T> *M, dd_ErrorType *error)
{
  /* This is to generate the (facet) graph of a polyheron (H) V-represented by M using LPs.
     Since it does not use the representation conversion, it should work for a large
     scale problem.
  */
  dd_rowrange i,m;
  dd_colrange d;
  dd_rowset redset;
  dd_matrixdata<T> *Mcopy;
  dd_SetFamilyPtr F=nullptr;

  m=M->rowsize;
  d=M->colsize;
  if (m<=0 ||d<=0) {
    *error=dd_EmptyRepresentation;
    goto _L999;
  }
  Mcopy=dd_MatrixCopy(M);
  F=dd_CreateSetFamily(m, m);
  for (i=1; i<=m; i++) {
    if (!set_member(i, M->linset)){
      set_addelem(Mcopy->linset, i);
      redset=dd_RedundantRows(Mcopy, error);  /* redset should contain all nonadjacent ones */
      set_uni(redset, redset, Mcopy->linset); /* all linearity elements should be nonadjacent */
      set_compl(F->set[i-1], redset); /* set the adjacency list of vertex i */
      set_delelem(Mcopy->linset, i);
      set_free(redset);
      if (*error!=dd_NoError) goto _L99;
    }
  }
_L99:
  dd_FreeMatrix(Mcopy);
_L999:
  return F;
}

template<typename T>
dd_SetFamilyPtr dd_Matrix2WeakAdjacency(dd_matrixdata<T> *M, dd_ErrorType *error)
{
  /* This is to generate the weak-adjacency (facet) graph of a polyheron (H) V-represented by M using LPs.
     Since it does not use the representation conversion, it should work for a large
     scale problem.
  */
  dd_rowrange i,m;
  dd_colrange d;
  dd_rowset redset;
  dd_matrixdata<T> *Mcopy;
  dd_SetFamilyPtr F=nullptr;

  m=M->rowsize;
  d=M->colsize;
  if (m<=0 ||d<=0) {
    *error=dd_EmptyRepresentation;
    goto _L999;
  }
  Mcopy=dd_MatrixCopy(M);
  F=dd_CreateSetFamily(m, m);
  for (i=1; i<=m; i++) {
    if (!set_member(i, M->linset)){
      set_addelem(Mcopy->linset, i);
      redset=dd_SRedundantRows(Mcopy, error);  /* redset should contain all weakly nonadjacent ones */
      set_uni(redset, redset, Mcopy->linset); /* all linearity elements should be nonadjacent */
      set_compl(F->set[i-1], redset); /* set the adjacency list of vertex i */
      set_delelem(Mcopy->linset, i);
      set_free(redset);
      if (*error!=dd_NoError) goto _L99;
    }
  }
_L99:
  dd_FreeMatrix(Mcopy);
_L999:
  return F;
}


template<typename T>
bool dd_ImplicitLinearity(dd_matrixdata<T> *M, dd_rowrange itest, T* certificate, dd_ErrorType *error)
  /* 092 */
{
  /* Checks whether the row itest is implicit linearity for the representation.
     All linearity rows are not checked and considered non implicit linearity (false).
     This code works for both H- and V-representations.  A certificate is
     given in the case of false, showing a feasible solution x satisfying the itest
     strict inequality for H-representation, a hyperplane RHS and normal (x_0, x) that
     separates the itest from the rest.  More explicitly, the LP to be setup is
     the same thing as redundancy case but with maximization:

     H-representation
       f* = maximize
         b_itest     + A_itest x
       subject to
         b_itest + 1 + A_itest x     >= 0 (relaxed inequality. This is not necessary but kept for simplicity of the code)
         b_{I-itest} + A_{I-itest} x >= 0 (all inequalities except for itest)
         b_L         + A_L x = 0.  (linearity)

     V-representation (=separation problem)
       f* = maximize
         b_itest x_0     + A_itest x
       subject to
         b_itest x_0     + A_itest x     >= -1 (again, this is not necessary but kept for simplicity.)
         b_{I-itest} x_0 + A_{I-itest} x >=  0 (all nonlinearity generators except for itest in one side)
         b_L x_0         + A_L x = 0.  (linearity generators)

    Here, the input matrix is considered as (b, A), i.e. b corresponds to the first column of input
    and the row indices of input is partitioned into I and L where L is the set of linearity.
    In both cases, the itest data is implicit linearity if and only if the optimal value f* is nonpositive.
    The certificate has dimension one more for V-representation case.
  */

  dd_colrange j;
  dd_lpdata<T> *lp;
  dd_lpsolution<T> *lps;
  dd_ErrorType err=dd_NoError;
  bool answer=false,localdebug=false;

  *error=dd_NoError;
  if (set_member(itest, M->linset)){
    if (localdebug) printf("The %ld th row is linearity and redundancy checking is skipped.\n",itest);
    goto _L99;
  }

  /* Create an LP data for redundancy checking */
  if (M->representation==dd_Generator){
    lp=dd_CreateLP_V_Redundancy(M, itest);
  } else {
    lp=dd_CreateLP_H_Redundancy(M, itest);
  }

  lp->objective = dd_LPmax;  /* the lp->objective is set by CreateLP* to LPmin */
  dd_LPSolve(lp,dd_choiceRedcheckAlgorithm,&err);
  if (err!=dd_NoError){
    *error=err;
    goto _L999;
  } else {
    for (j=0; j<lp->d; j++) {
      dd_set(certificate[j], lp->sol[j]);
    }

    if (lp->LPS==dd_Optimal && dd_EqualToZero(lp->optvalue)){
      answer=true;
      if (localdebug) fprintf(stderr,"==> %ld th data is an implicit linearity.\n",itest);
    } else {
      answer=false;
      if (localdebug) fprintf(stderr,"==> %ld th data is not an implicit linearity.\n",itest);
    }
  }
  _L999:
  dd_FreeLPData(lp);
_L99:
  return answer;
}


template<typename T>
int dd_FreeOfImplicitLinearity(dd_matrixdata<T> *M, T* certificate, dd_rowset *imp_linrows, dd_ErrorType *error)
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

    Here, the input matrix is considered as (b, A), i.e. b corresponds to the first column of input
    and the row indices of input is partitioned into I and L where L is the set of linearity.
    In both cases, any implicit linearity exists if and only if the optimal value f* is nonpositive.
    The certificate has dimension one more for V-representation case.
  */

  dd_lpdata<T> *lp;
  dd_rowrange i,m;
  dd_colrange j,d1;
  dd_ErrorType err=dd_NoError;
  T* cvec; /* certificate for implicit linearity */

  int answer=0,localdebug=false;

  *error=dd_NoError;
  /* Create an LP data for redundancy checking */
  if (M->representation==dd_Generator){
    lp=dd_CreateLP_V_ImplicitLinearity(M);
  } else {
    lp=dd_CreateLP_H_ImplicitLinearity(M);
  }

  dd_LPSolve(lp,dd_choiceRedcheckAlgorithm,&err);
  if (err!=dd_NoError){
    *error=err;
    goto _L999;
  } else {

    for (j=0; j<lp->d; j++) {
      dd_set(certificate[j], lp->sol[j]);
    }

    if (localdebug) dd_WriteLPResult(stderr,lp,err);

    /* *posset contains a set of row indices that are recognized as nonlinearity.  */

    if (M->representation==dd_Generator){
      d1=(M->colsize)+1;
    } else {
      d1=M->colsize;
    }
    m=M->rowsize;
    dd_InitializeArow(d1,&cvec);
    set_initialize(imp_linrows,m);

    if (lp->LPS==dd_Optimal){
      if (dd_Positive(lp->optvalue)){
        answer=1;
        if (localdebug) fprintf(stderr,"==> The matrix has no implicit linearity.\n");
      } else if (dd_Negative(lp->optvalue)) {
          answer=-1;
          if (localdebug) fprintf(stderr,"==> The matrix defines the trivial system.\n");
        } else {
            answer=0;
            if (localdebug) fprintf(stderr,"==> The matrix has some implicit linearity.\n");
          }
    } else {
          answer=-2;
          if (localdebug) fprintf(stderr,"==> The LP fails.\n");
    }
    if (answer==0){
      /* List the implicit linearity rows */
      for (i=m; i>=1; i--) {
        if (!set_member(i,lp->posset_extra)) {
          if (dd_ImplicitLinearity(M, i, cvec, error)) {
            set_addelem(*imp_linrows, i);
            if (localdebug) {
              fprintf(stderr," row %ld is implicit linearity\n",i);
              fprintf(stderr,"\n");
            }
          }
          if (*error!=dd_NoError) goto _L999;
        }
      }
    }  /* end of if (answer==0) */
    if (answer==-1) {
      for (i=m; i>=1; i--) set_addelem(*imp_linrows, i);
    } /* all rows are considered implicit linearity */

    dd_FreeArow(cvec);
  }
  _L999:
  dd_FreeLPData(lp);

  return answer;
}


template<typename T>
dd_rowset dd_ImplicitLinearityRows(dd_matrixdata<T> *M, dd_ErrorType *error)  /* 092 */
{
  dd_colrange d;
  dd_rowset imp_linset;
  T* cvec; /* certificate */
  int foi;
  bool localdebug=false;

  if (M->representation==dd_Generator){
    d=(M->colsize)+2;
  } else {
    d=M->colsize+1;
  }

  dd_InitializeArow(d,&cvec);
  if (localdebug) fprintf(stdout, "\ndd_ImplicitLinearityRows: Check whether the system contains any implicit linearity.\n");
  foi=dd_FreeOfImplicitLinearity(M, cvec, &imp_linset, error);
  if (localdebug){
    switch (foi) {
      case 1:
        fprintf(stdout, "  It is free of implicit linearity.\n");
        break;

      case 0:
        fprintf(stdout, "  It is not free of implicit linearity.\n");
        break;

    case -1:
        fprintf(stdout, "  The input system is trivial (i.e. the empty H-polytope or the V-rep of the whole space.\n");
        break;

    default:
        fprintf(stdout, "  The LP was not solved correctly.\n");
        break;

    }
  }

  dd_FreeArow(cvec);
  return imp_linset;
}

template<typename T>
bool dd_MatrixCanonicalizeLinearity(dd_matrixdata<T> **M, dd_rowset *impl_linset,dd_rowindex *newpos,
					  dd_ErrorType *error) /* 094 */
{
/* This is to recongnize all implicit linearities, and put all linearities at the top of
   the matrix.    All implicit linearities will be returned by *impl_linset.
*/
  dd_rowrange rank;
  dd_rowset linrows,ignoredrows,basisrows;
  dd_colset ignoredcols,basiscols;
  dd_rowrange i,k,m;
  dd_rowindex newpos1;
  bool success=false;

  linrows=dd_ImplicitLinearityRows(*M, error);
  if (*error!=dd_NoError) goto _L99;

  m=(*M)->rowsize;

  set_uni((*M)->linset, (*M)->linset, linrows);
      /* add the implicit linrows to the explicit linearity rows */

  /* To remove redundancy of the linearity part,
     we need to compute the rank and a basis of the linearity part. */
  set_initialize(&ignoredrows,  (*M)->rowsize);
  set_initialize(&ignoredcols,  (*M)->colsize);
  set_compl(ignoredrows,  (*M)->linset);
  rank=dd_MatrixRank(*M,ignoredrows,ignoredcols,&basisrows,&basiscols);
  set_diff(ignoredrows,  (*M)->linset, basisrows);
  dd_MatrixRowsRemove2(M,ignoredrows,newpos);

  dd_MatrixShiftupLinearity(M,&newpos1);

  for (i=1; i<=m; i++){
    k=(*newpos)[i];
    if (k>0) {
      (*newpos)[i]=newpos1[k];
    }
  }

  *impl_linset=linrows;
  success=true;
  delete [] newpos1;
  set_free(basisrows);
  set_free(basiscols);
  set_free(ignoredrows);
  set_free(ignoredcols);
_L99:
  return success;
}

template<typename T>
bool dd_MatrixCanonicalize(dd_matrixdata<T> **M, dd_rowset *impl_linset, dd_rowset *redset,
				 dd_rowindex *newpos, dd_ErrorType *error) /* 094 */
{
/* This is to find a canonical representation of a matrix *M by
   recognizing all implicit linearities and all redundancies.
   All implicit linearities will be returned by *impl_linset and
   redundancies will be returned by *redset.
*/
  dd_rowrange i,k,m;
  dd_rowindex newpos1, revpos;
  dd_rowset redset1;
  bool success=true;

  m=(*M)->rowsize;
  set_initialize(redset, m);
  revpos = new long[m+1];
  for (int i=0; i<=m; i++)
    revpos[i] = 0;

  success=dd_MatrixCanonicalizeLinearity(M, impl_linset, newpos, error);

  if (!success) goto _L99;

  for (i=1; i<=m; i++){
    k=(*newpos)[i];
    if (k>0) revpos[k]=i;  /* inverse of *newpos[] */
  }

  success=dd_MatrixRedundancyRemove(M, &redset1, &newpos1, error);  /* 094 */

  if (!success) goto _L99;

  for (i=1; i<=m; i++){
    k=(*newpos)[i];
    if (k>0) {
      (*newpos)[i]=newpos1[k];
      if (newpos1[k]<0) (*newpos)[i]=-revpos[-newpos1[k]];  /* update the certificate of its duplicate removal. */
      if (set_member(k,redset1)) set_addelem(*redset, i);
    }
  }

_L99:
  set_free(redset1);
  delete [] newpos1;
  delete [] revpos;
  return success;
}


template<typename T>
bool dd_ExistsRestrictedFace(dd_matrixdata<T> *M, dd_rowset R, dd_rowset S, dd_ErrorType *err)
/* 0.94 */
{
/* This function checkes if there is a point that satifies all the constraints of
the matrix M (interpreted as an H-representation) with additional equality contraints
specified by R and additional strict inequality constraints specified by S.
The set S is supposed to be disjoint from both R and M->linset.   When it is not,
the set S will be considered as S\(R U M->linset).
*/
  bool answer=false;
  dd_lpdata<T> *lp=nullptr;

/*
  printf("\n--- ERF ---\n");
  printf("R = "); set_write(R);
  printf("S = "); set_write(S);
*/

  lp=dd_Matrix2Feasibility2(M, R, S, err);

  if (*err!=dd_NoError) goto _L99;

/* Solve the LP by cdd LP solver. */
  dd_LPSolve(lp, dd_DualSimplex, err);  /* Solve the LP */
  if (*err!=dd_NoError) goto _L99;
  if (lp->LPS==dd_Optimal && dd_Positive(lp->optvalue)) {
    answer=true;
  }

  dd_FreeLPData(lp);
_L99:
  return answer;
}

template<typename T>
bool dd_ExistsRestrictedFace2(dd_matrixdata<T> *M, dd_rowset R, dd_rowset S, dd_lpsolution<T> **lps, dd_ErrorType *err)
/* 0.94 */
{
/* This function checkes if there is a point that satifies all the constraints of
the matrix M (interpreted as an H-representation) with additional equality contraints
specified by R and additional strict inequality constraints specified by S.
The set S is supposed to be disjoint from both R and M->linset.   When it is not,
the set S will be considered as S\(R U M->linset).

This function returns a certificate of the answer in terms of the associated LP solutions.
*/
  bool answer=false;
  dd_lpdata<T> *lp=nullptr;

/*
  printf("\n--- ERF ---\n");
  printf("R = "); set_write(R);
  printf("S = "); set_write(S);
*/

  lp=dd_Matrix2Feasibility2(M, R, S, err);

  if (*err!=dd_NoError) goto _L99;

/* Solve the LP by cdd LP solver. */
  dd_LPSolve(lp, dd_DualSimplex, err);  /* Solve the LP */
  if (*err!=dd_NoError) goto _L99;
  if (lp->LPS==dd_Optimal && dd_Positive(lp->optvalue)) {
    answer=true;
  }


  (*lps)=dd_CopyLPSolution(lp);
  dd_FreeLPData(lp);
_L99:
  return answer;
}

template<typename T>
bool dd_FindRelativeInterior(dd_matrixdata<T> *M, dd_rowset *ImL, dd_rowset *Lbasis, dd_lpsolution<T> **lps, dd_ErrorType *err)
/* 0.94 */
{
/* This function computes a point in the relative interior of the H-polyhedron given by M.
Even the representation is V-representation, it simply interprete M as H-representation.
lps returns the result of solving an LP whose solution is a relative interior point.
ImL returns all row indices of M that are implicit linearities, i.e. their inqualities
are satisfied by equality by all points in the polyhedron.  Lbasis returns a row basis
of the submatrix of M consisting of all linearities and implicit linearities.  This means
that the dimension of the polyhedron is M->colsize - set_card(Lbasis) -1.
*/

  dd_rowset S;
  dd_colset Tc, Lbasiscols;
  bool success=false;
  dd_rowrange i;
  dd_colrange rank;


  *ImL=dd_ImplicitLinearityRows(M, err);

  if (*err!=dd_NoError) goto _L99;

  set_initialize(&S, M->rowsize);   /* the empty set */
  for (i=1; i <=M->rowsize; i++) {
	if (!set_member(i, M->linset) && !set_member(i, *ImL)){
	  set_addelem(S, i);  /* all nonlinearity rows go to S  */
	}
  }
  if (dd_ExistsRestrictedFace2(M, *ImL, S, lps, err)){
    /* printf("a relative interior point found\n"); */
    success=true;
  }

  set_initialize(&Tc,  M->colsize); /* empty set */
  rank=dd_MatrixRank(M,S,Tc,Lbasis,&Lbasiscols); /* the rank of the linearity submatrix of M.  */

  set_free(S);
  set_free(Tc);
  set_free(Lbasiscols);

_L99:
  return success;
}


template<typename T>
dd_rowrange dd_RayShooting(dd_matrixdata<T> *M, T* p, T* r)
{
  dd_rowrange imin=-1, i, m;
  dd_colrange j, d;
  T *vecmin, *vec;
  T min, t1, t2, alpha, t1min;
  bool started=false;
  bool localdebug=false;
  T dd_one;
  dd_one=1;
  m=M->rowsize;
  d=M->colsize;
  if (!dd_Equal(dd_one, p[0])){
    fprintf(stderr, "Warning: RayShooting is called with a point with first coordinate not 1.\n");
    dd_set(p[0],dd_one);
  }
  if (!dd_EqualToZero(r[0])){
    fprintf(stderr, "Warning: RayShooting is called with a direction with first coordinate not 0.\n");
    r[0]=0;
  }

  dd_InitializeArow(d,&vecmin);
  dd_InitializeArow(d,&vec);

  for (i=1; i<=m; i++){
    dd_InnerProduct(t1, d, M->matrix[i-1], p);
    if (localdebug) std::cerr << "dd_RayShooting: i=" << i << " t1=" << t1 << "\n";
    if (dd_Positive(t1)) {
      dd_InnerProduct(t2, d, M->matrix[i-1], r);
      bool is_field=true;
      if (is_field) {
          dd_div(alpha, t2, t1);
          if (!started){
            imin=i;  dd_set(min, alpha);
            dd_set(t1min, t1);  /* store the denominator. */
            started=true;
            if (localdebug) std::cerr << "dd_RayShooting: Level 1: imin = "<< imin << " and min = " << min << "\n";
          } else {
            if (dd_Smaller(alpha, min)){
              imin=i;  dd_set(min, alpha);
              dd_set(t1min, t1);  /* store the denominator. */
              if (localdebug) std::cerr << "dd_RayShootni: Level 2: imin = " << imin << " and min = " << min << "\n";
            } else {
              if (dd_Equal(alpha, min)) { /* tie break */
                for (j=1; j<= d; j++) {
                  dd_div(vecmin[j-1], M->matrix[imin-1][j-1], t1min);
                  dd_div(vec[j-1], M->matrix[i-1][j-1], t1);
                }
                if (dd_LexSmaller(vec,vecmin, d)){
                  imin=i;  dd_set(min, alpha);
                  dd_set(t1min, t1);  /* store the denominator. */
                  if (localdebug) std::cerr << "dd_RayShooting: Level 3: imin = " << imin << " and min = " << min << "\n";
                }
              }
            }
          }
        } else {
          if (!started){
            imin=i;
            dd_set(min, t2);
            dd_set(t1min, t1);  /* store the denominator. */
            started=true;
            if (localdebug) std::cerr << "dd_RayShooting: Level 1: imin = "<< imin << " and min = " << min << "\n";
          } else {
            if (dd_SmallerFrac(t2, t1, min, t1min)){
              imin=i;
              dd_set(min, t2);
              dd_set(t1min, t1);  /* store the denominator. */
              if (localdebug) std::cerr << "dd_RayShootni: Level 2: imin = " << imin << " and min = " << min << "\n";
            } else {
              if (dd_EqualFrac(t2, t1, min, t1min)) { /* tie break */
                if (dd_LexSmallerFrac(M->matrix[i-1], t1, M->matrix[imin-1], t1min, d)){
                  imin=i;
                  dd_set(min, t2);
                  dd_set(t1min, t1);  /* store the denominator. */
                  if (localdebug) std::cerr << "dd_RayShooting: Level 3: imin = " << imin << " and min = " << min << "\n";
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

template<typename T>
void dd_BasisStatusMaximize(dd_rowrange m_size,dd_colrange d_size,
			    T** A,T** Ts,dd_rowset equalityset,
			    dd_rowrange objrow,dd_colrange rhscol,dd_LPStatusType LPS,
			    T &optvalue,T* sol,T* dsol,dd_rowset posset, dd_colindex nbindex,
			    dd_rowrange re,dd_colrange se, dd_colrange *nse, long *pivots,
			    int *found, int *LPScorrect, data_temp_simplex<T>* data)
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
  long pivots0,pivots1,fbasisrank;
  dd_rowrange i,is;
  dd_colrange s,senew,j;
  dd_rowindex bflag;
  dd_rowindex OrderVector;  /* the permutation vector to store a preordered row indices */
  unsigned int rseed=1;
  T val;
  dd_colindex nbtemp;
  //  dd_LPStatusType ddlps;
  bool localdebug=false;

  nbtemp = new long[d_size+1];
  for (int i=0; i<=d_size; i++)
    nbtemp[i] = 0;
  for (i=0; i<= 4; i++)
    pivots[i]=0;
  bflag = new long[m_size+1];
  for (int i=0; i<=m_size; i++)
    bflag[i] = 0;
  OrderVector = new long[m_size+1];
  for (int i=0; i<=m_size; i++)
    OrderVector[i] = 0;

  /* Initializing control variables. */
  dd_ComputeRowOrderVector2(m_size,d_size,A,OrderVector,dd_MinIndex,rseed);

  pivots1=0;

  dd_ResetTableau(m_size,d_size,Ts,nbtemp,bflag,objrow,rhscol);


  is=nbindex[se];

  fbasisrank=d_size-1;
  for (j=1; j<=d_size; j++){
    if (nbindex[j]<0) fbasisrank=fbasisrank-1;
	/* fbasisrank=the basis rank computed by floating-point */
  }

  if (fbasisrank<d_size-1) {
    *found=false;
    goto _L99;
    /* Suspicious case.  Rerun the LP solver with GMP. */
  }



  dd_FindLPBasis2(m_size,d_size,A,Ts,OrderVector, equalityset,nbindex,bflag,
		  objrow,rhscol,&s,found,&pivots0, data);

/* set up the new se column and corresponding variable */
  senew=bflag[is];
  is=nbindex[senew];

  pivots[4]=pivots0;  /*GMP postopt pivots */

  if (!(*found)){
    goto _L99;
  }


  /* Check whether a recomputed basis is of the type specified by LPS */
  *LPScorrect=true;
  switch (LPS){
     case dd_Optimal:
       for (i=1; i<=m_size; i++) {
         if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
            dd_TableauEntry(val, d_size,A,Ts,i,rhscol);
            if (dd_Negative(val)) {
               *LPScorrect=false;
               break;
            }
          } else if (bflag[i] >0) { /* i is nonbasic variable */
            dd_TableauEntry(val, d_size,A,Ts,objrow,bflag[i]);
            if (dd_Positive(val)) {
               *LPScorrect=false;
               break;
            }
          }
       };
       break;
     case dd_Inconsistent:
       for (j=1; j<=d_size; j++){
          dd_TableauEntry(val, d_size,A,Ts,re,j);
          if (j==rhscol){
	    if (dd_Nonnegative(val)){
               *LPScorrect=false;
               break;
             }
	  } else if (dd_Positive(val)){
               *LPScorrect=false;
               break;
           }
       };
       break;
     case dd_LPSundecided: break;
     case dd_StrucInconsistent: break;
     case dd_StrucDualInconsistent: break;
     case dd_Unbounded: break;
     case dd_DualUnbounded: break;

     case dd_DualInconsistent:
        for (i=1; i<=m_size; i++){
          dd_TableauEntry(val, d_size,A,Ts,i,bflag[is]);
          if (i==objrow){
	    if (dd_Nonpositive(val)){
               *LPScorrect=false;
               break;
             }
	  } else if (dd_Negative(val)){
               *LPScorrect=false;
               break;
           }
       };
       break;
;
  }

  dd_SetSolutions(m_size,d_size,A,Ts,
   objrow,rhscol,LPS,optvalue,sol,dsol,posset,re,senew,bflag);
  *nse=senew;


_L99:
  delete [] nbtemp;
  delete [] bflag;
  delete [] OrderVector;
}

template<typename T>
void dd_BasisStatusMinimize(dd_rowrange m_size,dd_colrange d_size,
    T** A,T** Ts,dd_rowset equalityset,
    dd_rowrange objrow,dd_colrange rhscol,dd_LPStatusType LPS,
    T & optvalue,T* sol,T* dsol, dd_rowset posset, dd_colindex nbindex,
    dd_rowrange re,dd_colrange se,dd_colrange *nse,long *pivots, int *found, int *LPScorrect)
{
   dd_colrange j;

   for (j=1; j<=d_size; j++) dd_neg(A[objrow-1][j-1],A[objrow-1][j-1]);
   dd_BasisStatusMaximize(m_size,d_size,A,Ts,equalityset, objrow,rhscol,
     LPS,optvalue,sol,dsol,posset,nbindex,re,se,nse,pivots,found,LPScorrect);
   dd_neg(*optvalue,*optvalue);
   for (j=1; j<=d_size; j++){
	if (LPS!=dd_Inconsistent) {
	   /* Inconsistent certificate stays valid for minimization, 0.94e */
       dd_neg(dsol[j-1],dsol[j-1]);
	 }
     dd_neg(A[objrow-1][j-1],A[objrow-1][j-1]);
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


template<typename T>
void dd_CheckAdjacency(dd_conedata<T> *cone,
    dd_raydata<T> **RP1, dd_raydata<T> **RP2, bool *adjacent)
{
  dd_raydata<T> *TempRay;
  bool localdebug=false;
  dd_rowset Face, Face1;
  set_initialize(&Face, cone->m);
  set_initialize(&Face1, cone->m);

  *adjacent = true;
  set_int(Face1, (*RP1)->ZeroSet, (*RP2)->ZeroSet);
  set_int(Face, Face1, cone->AddedHalfspaces);
  if (set_card(Face)< cone->d - 2) {
    *adjacent = false;
    if (localdebug) {
      fprintf(stderr,"non adjacent: set_card(face) %ld < %ld = cone->d.\n",
        set_card(Face),cone->d);
    }
    set_free(Face);
    set_free(Face1);
    return;
  }
  else if (cone->parent->NondegAssumed) {
    *adjacent = true;
    set_free(Face);
    set_free(Face1);
    return;
  }
  TempRay = cone->FirstRay;
  while (TempRay != nullptr && *adjacent) {
    if (TempRay != *RP1 && TempRay != *RP2) {
    	set_int(Face1, TempRay->ZeroSet, cone->AddedHalfspaces);
      	if (set_subset(Face, Face1)) *adjacent = false;
    }
    TempRay = TempRay->Next;
  }
  set_free(Face);
  set_free(Face1);
}

template<typename T>
void dd_Eliminate(dd_conedata<T> *cone, dd_raydata<T> **Ptr)
{
  /*eliminate the record pointed by Ptr->Next*/
  dd_raydata<T> *TempPtr;

  TempPtr = (*Ptr)->Next;
  (*Ptr)->Next = (*Ptr)->Next->Next;
  if (TempPtr == cone->FirstRay)   /*Update the first pointer*/
    cone->FirstRay = (*Ptr)->Next;
  if (TempPtr == cone->LastRay)   /*Update the last pointer*/
    cone->LastRay = *Ptr;

  delete [] TempPtr->Ray;          /* free the ray vector memory */
  set_free(TempPtr->ZeroSet);  /* free the ZeroSet memory */
  delete TempPtr;   /* free the dd_Ray structure memory */
  cone->RayCount--;
}

template<typename T>
void dd_SetInequalitySets(dd_conedata<T> *cone)
{
  dd_rowrange i;

  set_emptyset(cone->GroundSet);
  set_emptyset(cone->EqualitySet);
  set_emptyset(cone->NonequalitySet);
  for (i = 1; i <= (cone->parent->m); i++){
    set_addelem(cone->GroundSet, i);
    if (cone->parent->EqualityIndex[i]==1) set_addelem(cone->EqualitySet,i);
    if (cone->parent->EqualityIndex[i]==-1) set_addelem(cone->NonequalitySet,i);
  }
}


template<typename T>
void dd_AValue(T *val, dd_colrange d_size, T** A, T *p, dd_rowrange i)
{
  /*return the ith component of the vector  A x p */
  dd_colrange j;
  T x;

  *val=0;
 /* Changed by Marc Pfetsch 010219 */

  for (j = 0; j < d_size; j++){
    T aVal=A[i - 1][j];
    T bVal=p[j];
    dd_mul(x, aVal, bVal);
    dd_add(*val, *val, x);
  }
}

template<typename T>
void dd_StoreRay1(dd_conedata<T> *cone, T *p, bool *feasible)
{  /* Original ray storing routine when RelaxedEnumeration is false */
  dd_rowrange i,k,fii=cone->m+1;
  dd_colrange j;
  T temp;
  dd_raydata<T> *RR;
  bool localdebug=false;

  RR=cone->LastRay;
  *feasible = true;
  set_initialize(&(RR->ZeroSet),cone->m);
  for (j = 0; j < cone->d; j++){
    dd_set(RR->Ray[j],p[j]);
  }
  for (i = 1; i <= cone->m; i++) {
    k=cone->OrderVector[i];
    dd_AValue(&temp, cone->d, cone->A, p, k);
    if (dd_EqualToZero(temp)) {
      set_addelem(RR->ZeroSet, k);
      if (localdebug) {
        fprintf(stderr,"recognized zero!\n");
      }
    }
    if (dd_Negative(temp)){
      if (localdebug) {
        fprintf(stderr,"recognized negative!\n");
      }
      *feasible = false;
      if (fii>cone->m) fii=i;  /* the first violating inequality index */
    }
  }
  RR->FirstInfeasIndex=fii;
  RR->feasible = *feasible;
}

template<typename T>
void dd_StoreRay2(dd_conedata<T> *cone, T *p,
		  bool *feasible, bool *weaklyfeasible)
   /* Ray storing routine when RelaxedEnumeration is true.
       weaklyfeasible is true iff it is feasible with
       the strict_inequality conditions deleted. */
{
  dd_raydata<T> *RR;
  dd_rowrange i,k,fii=cone->m+1;
  dd_colrange j;
  T temp;

  RR=cone->LastRay;
  *feasible = true;
  *weaklyfeasible = true;
  set_initialize(&(RR->ZeroSet),cone->m);
  for (j = 0; j < cone->d; j++){
    dd_set(RR->Ray[j],p[j]);
  }
  for (i = 1; i <= cone->m; i++) {
    k=cone->OrderVector[i];
    dd_AValue(&temp, cone->d, cone->A, p, k);
    if (dd_EqualToZero(temp)){
      set_addelem(RR->ZeroSet, k);
      if (cone->parent->EqualityIndex[k]==-1)
        *feasible=false;  /* strict inequality required */
    }
/*    if (temp < -zero){ */
    if (dd_Negative(temp)){
      *feasible = false;
      if (fii>cone->m && cone->parent->EqualityIndex[k]>=0) {
        fii=i;  /* the first violating inequality index */
        *weaklyfeasible=false;
      }
    }
  }
  RR->FirstInfeasIndex=fii;
  RR->feasible = *feasible;
}


template<typename T>
void dd_AddRay(dd_conedata<T> *cone, T *p)
{
  bool feasible, weaklyfeasible;
  bool localdebug=false;

  if (cone->FirstRay == nullptr) {
    cone->FirstRay = new dd_raydata<T>;
    cone->FirstRay->Ray = new T[cone->d];
    if (localdebug)
      fprintf(stderr,"Create the first ray pointer\n");
    cone->LastRay = cone->FirstRay;
    cone->ArtificialRay->Next = cone->FirstRay;
  } else {
    cone->LastRay->Next = new dd_raydata<T>;
    cone->LastRay->Next->Ray = new T[cone->d];
    if (localdebug) fprintf(stderr,"Create a new ray pointer\n");
    cone->LastRay = cone->LastRay->Next;
  }
  cone->LastRay->Next = nullptr;
  cone->RayCount++;
  cone->TotalRayCount++;
  if (localdebug) {
    if (cone->TotalRayCount % 100 == 0) {
      fprintf(stderr,"*Rays (Total, Currently Active, Feasible) =%8ld%8ld%8ld\n",
	 cone->TotalRayCount, cone->RayCount, cone->FeasibleRayCount);
    }
  }
  if (cone->parent->RelaxedEnumeration){
    dd_StoreRay2(cone, p, &feasible, &weaklyfeasible);
    if (weaklyfeasible) (cone->WeaklyFeasibleRayCount)++;
  } else {
    dd_StoreRay1(cone, p, &feasible);
    if (feasible) (cone->WeaklyFeasibleRayCount)++;
    /* weaklyfeasible is equiv. to feasible in this case. */
  }
  if (!feasible) return;
  else {
    (cone->FeasibleRayCount)++;
  }
}

template<typename T>
void dd_AddArtificialRay(dd_conedata<T> *cone)
{
  T* zerovector;
  dd_colrange d1;
  bool feasible;
  bool localdebug=false;

  if (cone->d<=0) d1=1; else d1=cone->d;
  dd_InitializeArow(d1, &zerovector);
  if (cone->ArtificialRay != nullptr) {
    fprintf(stderr,"Warning !!!  FirstRay in not nil.  Illegal Call\n");
    delete [] zerovector; /* 086 */
    return;
  }
  cone->ArtificialRay = new dd_raydata<T>;
  cone->ArtificialRay->Ray = new T[d1];

  if (localdebug) fprintf(stderr,"Create the artificial ray pointer\n");

  cone->LastRay=cone->ArtificialRay;
  dd_StoreRay1(cone, zerovector, &feasible);
  cone->ArtificialRay->Next = nullptr;
  delete [] zerovector; /* 086 */
}

template<typename T>
void dd_ConditionalAddEdge(dd_conedata<T> *cone,
    dd_raydata<T> *Ray1, dd_raydata<T> *Ray2, dd_raydata<T> *ValidFirstRay)
{
  long it,it_row,fii1,fii2,fmin,fmax;
  bool adjacent,lastchance;
  dd_raydata<T> *TempRay;
  dd_raydata<T> *Rmin;
  dd_raydata<T> *Rmax;
  dd_adjacencydata<T> *NewEdge;
  bool localdebug=false;
  dd_rowset ZSmin, ZSmax;
  dd_rowset Face, Face1;

  set_initialize(&Face, cone->m);
  set_initialize(&Face1, cone->m);

  fii1=Ray1->FirstInfeasIndex;
  fii2=Ray2->FirstInfeasIndex;
  if (fii1<fii2){
    fmin=fii1; fmax=fii2;
    Rmin=Ray1;
    Rmax=Ray2;
  }
  else{
    fmin=fii2; fmax=fii1;
    Rmin=Ray2;
    Rmax=Ray1;
  }
  ZSmin = Rmin->ZeroSet;
  ZSmax = Rmax->ZeroSet;
  if (localdebug) {
    fprintf(stderr,"dd_ConditionalAddEdge: FMIN = %ld (row%ld)   FMAX=%ld\n",
      fmin, cone->OrderVector[fmin], fmax);
  }
  if (fmin==fmax){
    if (localdebug) fprintf(stderr,"dd_ConditionalAddEdge: equal FII value-> No edge added\n");
  }
  else if (set_member(cone->OrderVector[fmin],ZSmax)){
    if (localdebug) fprintf(stderr,"dd_ConditionalAddEdge: No strong separation -> No edge added\n");
  }
  else {  /* the pair will be separated at the iteration fmin */
    lastchance=true;
    /* flag to check it will be the last chance to store the edge candidate */
    set_int(Face1, ZSmax, ZSmin);
    (cone->count_int)++;
    if (localdebug){
      fprintf(stderr,"Face: ");
      for (it=1; it<=cone->m; it++) {
        it_row=cone->OrderVector[it];
        if (set_member(it_row, Face1)) fprintf(stderr,"%ld ",it_row);
      }
      fprintf(stderr,"\n");
    }
    for (it=cone->Iteration+1; it < fmin && lastchance; it++){
      it_row=cone->OrderVector[it];
      if (cone->parent->EqualityIndex[it_row]>=0 && set_member(it_row, Face1)){
        lastchance=false;
        (cone->count_int_bad)++;
      }
    }
    if (lastchance){
      adjacent = true;
      (cone->count_int_good)++;
      /* adjacent checking */
      set_int(Face, Face1, cone->AddedHalfspaces);
      if (set_card(Face)< cone->d - 2) {
        adjacent = false;
      }
      else if (cone->parent->NondegAssumed) {
    	adjacent = true;
      }
      else{
        TempRay = ValidFirstRay;  /* the first ray for adjacency checking */
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
      if (adjacent){
        NewEdge=new dd_adjacencydata<T>;
        NewEdge->Ray1=Rmax;  /* save the one remains in iteration fmin in the first */
        NewEdge->Ray2=Rmin;  /* save the one deleted in iteration fmin in the second */
        NewEdge->Next=nullptr;
        (cone->EdgeCount)++;
        (cone->TotalEdgeCount)++;
        if (cone->Edges[fmin] == nullptr){
          cone->Edges[fmin]=NewEdge;
        }else{
          NewEdge->Next=cone->Edges[fmin];
          cone->Edges[fmin]=NewEdge;
        }
      }
    }
  }
  set_free(Face);
  set_free(Face1);
}

template<typename T>
void dd_CreateInitialEdges(dd_conedata<T> *cone)
{
  dd_raydata<T> *Ptr1;
  dd_raydata<T> *Ptr2;
  dd_rowrange fii1,fii2;
  long count=0;
  bool adj;

  cone->Iteration=cone->d;  /* CHECK */
  if (cone->FirstRay ==nullptr || cone->LastRay==nullptr){
    /* fprintf(stderr,"Warning: dd_ CreateInitialEdges called with nullptr pointer(s)\n"); */
    goto _L99;
  }
  Ptr1=cone->FirstRay;
  while(Ptr1!=cone->LastRay && Ptr1!=nullptr){
    fii1=Ptr1->FirstInfeasIndex;
    Ptr2=Ptr1->Next;
    while(Ptr2!=nullptr){
      fii2=Ptr2->FirstInfeasIndex;
      count++;
      dd_CheckAdjacency(cone, &Ptr1, &Ptr2, &adj);
      if (fii1 != fii2 && adj)
        dd_ConditionalAddEdge(cone, Ptr1, Ptr2, cone->FirstRay);
      Ptr2=Ptr2->Next;
    }
    Ptr1=Ptr1->Next;
  }
_L99:;
}


template<typename T>
void dd_UpdateEdges(dd_conedata<T> *cone, dd_raydata<T> *RRbegin, dd_raydata<T> *RRend)
/* This procedure must be called after the ray list is sorted
   by dd_EvaluateARay2 so that FirstInfeasIndex's are monotonically
   increasing.
*/
{
  dd_raydata<T> *Ptr1;
  dd_raydata<T> *Ptr2;
  dd_raydata<T> *Ptr2begin;
  dd_rowrange fii1;
  bool ptr2found,quit;
  long count=0,pos1, pos2;
  float workleft,prevworkleft=110.0,totalpairs;
  bool localdebug=false;

  totalpairs=(cone->ZeroRayCount-1.0)*(cone->ZeroRayCount-2.0)+1.0;
  Ptr2begin = nullptr;
  if (RRbegin ==nullptr || RRend==nullptr){
    if (1) fprintf(stderr,"Warning: dd_UpdateEdges called with nullptr pointer(s)\n");
    goto _L99;
  }
  Ptr1=RRbegin;
  pos1=1;
  do{
    ptr2found=false;
    quit=false;
    fii1=Ptr1->FirstInfeasIndex;
    pos2=2;
    for (Ptr2=Ptr1->Next; !ptr2found && !quit; Ptr2=Ptr2->Next,pos2++){
      if  (Ptr2->FirstInfeasIndex > fii1){
        Ptr2begin=Ptr2;
        ptr2found=true;
      }
      else if (Ptr2==RRend) quit=true;
    }
    if (ptr2found){
      quit=false;
      for (Ptr2=Ptr2begin; !quit ; Ptr2=Ptr2->Next){
        count++;
        dd_ConditionalAddEdge(cone, Ptr1,Ptr2,RRbegin);
        if (Ptr2==RRend || Ptr2->Next==nullptr) quit=true;
      }
    }
    Ptr1=Ptr1->Next;
    pos1++;
    workleft = 100.0 * (cone->ZeroRayCount-pos1) * (cone->ZeroRayCount - pos1-1.0) / totalpairs;
    if (localdebug) {
      if (cone->ZeroRayCount>=500 && pos1%10 ==0 && prevworkleft-workleft>=10 ) {
        fprintf(stderr,"*Work of iteration %5ld(/%ld): %4ld/%4ld => %4.1f%% left\n",
                cone->Iteration, cone->m, pos1, cone->ZeroRayCount, workleft);
        prevworkleft=workleft;
      }
    }
  } while (Ptr1!=RRend && Ptr1!=nullptr);
_L99:;
}

template<typename T>
void dd_FreeDDMemory0(dd_conedata<T> *cone)
{
  dd_raydata<T> *Ptr;
  dd_raydata<T> *PrevPtr;
  long count;

  /* THIS SHOULD BE REWRITTEN carefully */
  PrevPtr=cone->ArtificialRay;
  if (PrevPtr!=nullptr){
    count=0;
    for (Ptr=cone->ArtificialRay->Next; Ptr!=nullptr; Ptr=Ptr->Next){
      delete [] PrevPtr->Ray;
      delete [] PrevPtr->ZeroSet;
      delete PrevPtr;
      count++;
      PrevPtr=Ptr;
    };
    cone->FirstRay=nullptr;

    delete [] cone->LastRay->Ray;
    cone->LastRay->Ray = nullptr;
    set_free(cone->LastRay->ZeroSet);
    cone->LastRay->ZeroSet = nullptr;
    delete cone->LastRay;
    cone->LastRay = nullptr;
    cone->ArtificialRay=nullptr;
  }
/* must add (by Sato) */
  delete [] cone->Edges;

  set_free(cone->GroundSet);
  set_free(cone->EqualitySet);
  set_free(cone->NonequalitySet);
  set_free(cone->AddedHalfspaces);
  set_free(cone->WeaklyAddedHalfspaces);
  set_free(cone->InitialHalfspaces);
  delete [] cone->InitialRayIndex;
  delete [] cone->OrderVector;
  delete [] cone->newcol;

/* Fixed by Shawn Rusaw.  Originally it was cone->d instead of cone->d_alloc */
  dd_FreeBmatrix(cone->d_alloc,cone->B);
  dd_FreeBmatrix(cone->d_alloc,cone->Bsave);

/* Fixed by Marc Pfetsch 010219*/
  dd_FreeAmatrix(cone->m_alloc, cone->A);
  cone->A = nullptr;

  delete cone;
}

template<typename T>
void dd_FreeDDMemory(dd_polyhedradata<T> *poly)
{
  dd_FreeDDMemory0(poly->child);
  poly->child=nullptr;
}

template<typename T>
void dd_FreePolyhedra(dd_polyhedradata<T> *poly)
{
  dd_bigrange i;

  if ((poly)->child != nullptr) dd_FreeDDMemory(poly);
  dd_FreeAmatrix((poly)->m_alloc, poly->A);
  dd_FreeArow((poly)->c);
  delete [] (poly)->EqualityIndex;
  if (poly->AincGenerated){
    for (i=1; i<=poly->m1; i++){
      set_free(poly->Ainc[i-1]);
    }
    delete [] poly->Ainc;
    set_free(poly->Ared);
    set_free(poly->Adom);
    poly->Ainc=nullptr;
  }

  delete poly;
}

template<typename T>
void dd_ZeroIndexSet(dd_rowrange m_size, dd_colrange d_size, T** A, T *x, dd_rowset ZS)
{
  dd_rowrange i;
  T temp;

  /* Changed by Marc Pfetsch 010219 */
  set_emptyset(ZS);
  for (i = 1; i <= m_size; i++) {
    dd_AValue(&temp, d_size, A, x, i);
    if (dd_EqualToZero(temp)) set_addelem(ZS, i);
  }

}

template<typename T>
void dd_CopyBmatrix(dd_colrange d_size, T** Ts, T** TCOPY)
{
  dd_rowrange i;
  dd_colrange j;

  for (i=0; i < d_size; i++) {
    for (j=0; j < d_size; j++) {
      dd_set(TCOPY[i][j],Ts[i][j]);
    }
  }
}


template<typename T>
void dd_CopyNormalizedArow(T *acopy, T *a, dd_colrange d)
{
  dd_CopyArow(acopy, a, d);
}


template<typename T>
void dd_CopyNormalizedAmatrix(T **Acopy, T **A, dd_rowrange m, dd_colrange d)
{
  dd_rowrange i;

  for (i = 0; i< m; i++) {
    dd_CopyNormalizedArow(Acopy[i],A[i],d);
  }
}

template<typename T>
void dd_PermuteCopyAmatrix(T **Acopy, T **A, dd_rowrange m, dd_colrange d, dd_rowindex roworder)
{
  dd_rowrange i;

  for (i = 1; i<= m; i++) {
    dd_CopyArow(Acopy[i-1],A[roworder[i]-1],d);
  }
}

template<typename T>
void dd_PermutePartialCopyAmatrix(T **Acopy, T **A, dd_rowrange m, dd_colrange d, dd_rowindex roworder)
{
 /* copy the rows of A whose roworder is positive.  roworder[i] is the row index of the copied row. */
  dd_rowrange i,k;

  k=0;
  for (i = 1; i<= m; i++) {
    if (roworder[i]>0) dd_CopyArow(Acopy[roworder[i]-1],A[i-1],d);
  }
}







template<typename T>
void dd_ColumnReduce(dd_conedata<T> *cone)
{
  dd_colrange j,j1=0;
  dd_rowrange i;

  for (j=1;j<=cone->d;j++) {
    if (cone->InitialRayIndex[j]>0){
      j1=j1+1;
      if (j1<j) {
        for (i=1; i<=cone->m; i++) dd_set(cone->A[i-1][j1-1],cone->A[i-1][j-1]);
        cone->newcol[j]=j1;
      }
    } else {
      cone->newcol[j]=0;
    }
  }
  cone->d=j1;  /* update the dimension. cone->d_orig remembers the old. */
  dd_CopyBmatrix(cone->d_orig, cone->B, cone->Bsave);
    /* save the dual basis inverse as Bsave.  This matrix contains the linearity space generators. */
  cone->ColReduced=true;
}

template<typename T>
long dd_MatrixRank(dd_matrixdata<T> *M, dd_rowset ignoredrows, dd_colset ignoredcols, dd_rowset *rowbasis, dd_colset *colbasis)
{
  bool stop, chosen;
  dd_rowset NopivotRow,PriorityRow;
  dd_colset ColSelected;
  T** B;
  dd_rowindex roworder;
  dd_rowrange r;
  dd_colrange s;
  long rank;
  bool localdebug=false;
  T* Rtemp;
  Rtemp = new T[M->colsize];

  rank = 0;
  stop = false;
  set_initialize(&ColSelected, M->colsize);
  set_initialize(&NopivotRow, M->rowsize);
  set_initialize(rowbasis, M->rowsize);
  set_initialize(colbasis, M->colsize);
  set_initialize(&PriorityRow, M->rowsize);
  set_copy(NopivotRow,ignoredrows);
  set_copy(ColSelected,ignoredcols);
  dd_InitializeBmatrix(M->colsize, &B);
  dd_SetToIdentity(M->colsize, B);
  roworder = new long[M->rowsize+1];
  for (r=0; r<M->rowsize; r++)
    roworder[r+1]=r+1;
  roworder[M->rowsize] = 0;

  do {   /* Find a set of rows for a basis */
      dd_SelectPivot2(M->rowsize, M->colsize, M->matrix, B, roworder,
		      PriorityRow,M->rowsize, NopivotRow, ColSelected, &r, &s, &chosen);
      if (localdebug && chosen)
        fprintf(stderr,"Procedure dd_MatrixRank: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(NopivotRow, r);
        set_addelem(*rowbasis, r);
        set_addelem(ColSelected, s);
        set_addelem(*colbasis, s);
        rank++;
        //        std::cout << "dd_GaussianColumnPivot call 8\n";
        dd_GaussianColumnPivot(M->colsize, M->matrix, B, r, s, Rtemp);
      } else {
        stop=true;
      }
      if (rank==M->colsize) stop = true;
  } while (!stop);
  dd_FreeBmatrix(M->colsize,B);
  delete [] roworder;
  set_free(ColSelected);
  set_free(NopivotRow);
  set_free(PriorityRow);
  delete [] Rtemp;
  return rank;
}


template<typename T>
void dd_FindBasis(dd_conedata<T> *cone, long *rank)
{
  bool stop, chosen;
  dd_rowset NopivotRow;
  dd_colset ColSelected;
  dd_rowrange r;
  dd_colrange j,s;
  bool localdebug=false;
  T* Rtemp;
  Rtemp = new T[cone->d];

  *rank = 0;
  stop = false;
  for (j=0;j<=cone->d;j++) cone->InitialRayIndex[j]=0;
  set_emptyset(cone->InitialHalfspaces);
  set_initialize(&ColSelected, cone->d);
  set_initialize(&NopivotRow, cone->m);
  set_copy(NopivotRow,cone->NonequalitySet);
  dd_SetToIdentity(cone->d, cone->B);
  do {   /* Find a set of rows for a basis */
      dd_SelectPivot2(cone->m, cone->d, cone->A, cone->B, cone->OrderVector,
		      cone->EqualitySet,cone->m, NopivotRow, ColSelected, &r, &s, &chosen);
      if (localdebug && chosen)
        fprintf(stderr,"Procedure dd_FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(cone->InitialHalfspaces, r);
        set_addelem(NopivotRow, r);
        set_addelem(ColSelected, s);
        cone->InitialRayIndex[s]=r;    /* cone->InitialRayIndex[s] stores the corr. row index */
        (*rank)++;
        //        std::cout << "dd_GaussianColumnPivot call 9\n";
        dd_GaussianColumnPivot(cone->d, cone->A, cone->B, r, s, Rtemp);
      } else {
        stop=true;
      }
      if (*rank==cone->d) stop = true;
  } while (!stop);
  set_free(ColSelected);
  set_free(NopivotRow);
  delete [] Rtemp;
}


template<typename T>
void dd_FindInitialRays(dd_conedata<T> *cone, bool *found)
{
  dd_rowset CandidateRows;
  dd_rowrange i;
  long rank;
  dd_RowOrderType roworder_save=dd_LexMin;
  bool localdebug=false;

  *found = false;
  set_initialize(&CandidateRows, cone->m);
  if (cone->parent->InitBasisAtBottom==true) {
    roworder_save=cone->HalfspaceOrder;
    cone->HalfspaceOrder=dd_MaxIndex;
    cone->PreOrderedRun=false;
  }
  else cone->PreOrderedRun=true;
  for (i = 1; i <= cone->m; i++)
    if (!set_member(i,cone->NonequalitySet)) set_addelem(CandidateRows, i);
    /*all rows not in NonequalitySet are candidates for initial cone*/
  dd_FindBasis(cone, &rank);
  cone->LinearityDim=cone->d - rank;
  if (localdebug) fprintf(stderr,"Linearity Dimension = %ld\n", cone->LinearityDim);
  if (cone->LinearityDim > 0) {
     dd_ColumnReduce(cone);
     dd_FindBasis(cone, &rank);
  }
  if (!set_subset(cone->EqualitySet,cone->InitialHalfspaces)) {
    if (localdebug) {
      fprintf(stderr,"Equality set is dependent. Equality Set and an initial basis:\n");
    };
  }
  *found = true;
  set_free(CandidateRows);
  if (cone->parent->InitBasisAtBottom==true) {
    cone->HalfspaceOrder=roworder_save;
  }
  if (cone->HalfspaceOrder==dd_MaxCutoff||
      cone->HalfspaceOrder==dd_MinCutoff||
      cone->HalfspaceOrder==dd_MixCutoff){
    cone->PreOrderedRun=false;
  } else cone->PreOrderedRun=true;
}

template<typename T>
void dd_CheckEquality(dd_colrange d_size, dd_raydata<T> **RP1, dd_raydata<T> **RP2, bool *equal)
{
  long j;
  bool localdebug=false;

  if (localdebug)
    fprintf(stderr,"Check equality of two rays\n");
  *equal = true;
  j = 1;
  while (j <= d_size && *equal) {
    if (!dd_Equal((*RP1)->Ray[j - 1],(*RP2)->Ray[j - 1]))
      *equal = false;
    j++;
  }
  if (*equal)
    fprintf(stderr,"Equal records found !!!!\n");
}

template<typename T>
void dd_CreateNewRay(dd_conedata<T> *cone,
		     dd_raydata<T> *Ptr1, dd_raydata<T> *Ptr2, dd_rowrange ii)
{
  /*Create a new ray by taking a linear combination of two rays*/
  dd_colrange j;
  T a1, a2, v1, v2;
  T* NewRay;
  NewRay=new T[cone->d];

  dd_AValue(&a1, cone->d, cone->A, Ptr1->Ray, ii);
  dd_AValue(&a2, cone->d, cone->A, Ptr2->Ray, ii);
  dd_abs(v1,a1);
  dd_abs(v2,a2);
  for (j = 0; j < cone->d; j++){
    dd_LinearComb(NewRay[j], Ptr1->Ray[j],v2,Ptr2->Ray[j],v1);
  }
  dd_AddRay(cone, NewRay);
  delete [] NewRay;
}

template<typename T>
void dd_EvaluateARay1(dd_rowrange i, dd_conedata<T> *cone)
/* Evaluate the ith component of the vector  A x RD.Ray
    and rearrange the linked list so that
    the infeasible rays with respect to  i  will be
    placed consecutively from First
 */
{
  dd_colrange j;
  T temp,tnext;
  dd_raydata<T> *Ptr;
  dd_raydata<T> *PrevPtr;
  dd_raydata<T> *TempPtr;

  Ptr = cone->FirstRay;
  PrevPtr = cone->ArtificialRay;
  if (PrevPtr->Next != Ptr) {
    fprintf(stderr,"Error.  Artificial Ray does not point to FirstRay!!!\n");
  }
  while (Ptr != nullptr) {
    temp=0;
    for (j = 0; j < cone->d; j++){
      dd_mul(tnext,cone->A[i - 1][j],Ptr->Ray[j]);
      dd_add(temp,temp,tnext);
    }
    dd_set(Ptr->ARay,temp);
/*    if ( temp <= -zero && Ptr != cone->FirstRay) {*/
    if ( dd_Negative(temp) && Ptr != cone->FirstRay) {
      /* fprintf(stderr,"Moving an infeasible record w.r.t. %ld to FirstRay\n",i); */
      if (Ptr==cone->LastRay) cone->LastRay=PrevPtr;
      TempPtr=Ptr;
      Ptr = Ptr->Next;
      PrevPtr->Next = Ptr;
      cone->ArtificialRay->Next = TempPtr;
      TempPtr->Next = cone->FirstRay;
      cone->FirstRay = TempPtr;
    }
    else {
      PrevPtr = Ptr;
      Ptr = Ptr->Next;
    }
  }
}

template<typename T>
void dd_EvaluateARay2(dd_rowrange i, dd_conedata<T> *cone)
/* Evaluate the ith component of the vector  A x RD.Ray
   and rearrange the linked list so that
   the infeasible rays with respect to  i  will be
   placed consecutively from First. Also for all feasible rays,
   "positive" rays and "zero" rays will be placed consecutively.
 */
{
  dd_colrange j;
  T temp,tnext;
  dd_raydata<T> *Ptr;
  dd_raydata<T> *NextPtr;
  bool zerofound=false,negfound=false,posfound=false;

  if (cone==nullptr || cone->TotalRayCount<=0) goto _L99;
  cone->PosHead=nullptr;cone->ZeroHead=nullptr;cone->NegHead=nullptr;
  cone->PosLast=nullptr;cone->ZeroLast=nullptr;cone->NegLast=nullptr;
  Ptr = cone->FirstRay;
  while (Ptr != nullptr) {
    NextPtr=Ptr->Next;  /* remember the Next record */
    Ptr->Next=nullptr;     /* then clear the Next pointer */
    temp=0;
    for (j = 0; j < cone->d; j++){
      T aVal=cone->A[i - 1][j];
      T bVal=Ptr->Ray[j];
      dd_mul(tnext,aVal,bVal);
      dd_add(temp,temp,tnext);
    }
    dd_set(Ptr->ARay,temp);
/*    if ( temp < -zero) {*/
    if ( dd_Negative(temp)) {
      if (!negfound){
        negfound=true;
        cone->NegHead=Ptr;
        cone->NegLast=Ptr;
      }
      else{
        Ptr->Next=cone->NegHead;
        cone->NegHead=Ptr;
      }
    }
/*    else if (temp > zero){*/
    else if (dd_Positive(temp)){
      if (!posfound){
        posfound=true;
        cone->PosHead=Ptr;
        cone->PosLast=Ptr;
      }
      else{
        Ptr->Next=cone->PosHead;
        cone->PosHead=Ptr;
       }
    }
    else {
      if (!zerofound){
        zerofound=true;
        cone->ZeroHead=Ptr;
        cone->ZeroLast=Ptr;
      }
      else{
        Ptr->Next=cone->ZeroHead;
        cone->ZeroHead=Ptr;
      }
    }
    Ptr=NextPtr;
  }
  /* joining three neg, pos and zero lists */
  if (negfound){                 /* -list nonempty */
    cone->FirstRay=cone->NegHead;
    if (posfound){               /* -list & +list nonempty */
      cone->NegLast->Next=cone->PosHead;
      if (zerofound){            /* -list, +list, 0list all nonempty */
        cone->PosLast->Next=cone->ZeroHead;
        cone->LastRay=cone->ZeroLast;
      }
      else{                      /* -list, +list nonempty but  0list empty */
        cone->LastRay=cone->PosLast;
      }
    }
    else{                        /* -list nonempty & +list empty */
      if (zerofound){            /* -list,0list nonempty & +list empty */
        cone->NegLast->Next=cone->ZeroHead;
        cone->LastRay=cone->ZeroLast;
      }
      else {                      /* -list nonempty & +list,0list empty */
        cone->LastRay=cone->NegLast;
      }
    }
  }
  else if (posfound){            /* -list empty & +list nonempty */
    cone->FirstRay=cone->PosHead;
    if (zerofound){              /* -list empty & +list,0list nonempty */
      cone->PosLast->Next=cone->ZeroHead;
      cone->LastRay=cone->ZeroLast;
    }
    else{                        /* -list,0list empty & +list nonempty */
      cone->LastRay=cone->PosLast;
    }
  }
  else{                          /* -list,+list empty & 0list nonempty */
    cone->FirstRay=cone->ZeroHead;
    cone->LastRay=cone->ZeroLast;
  }
  cone->ArtificialRay->Next=cone->FirstRay;
  cone->LastRay->Next=nullptr;
  _L99:;
}

template<typename T>
void dd_DeleteNegativeRays(dd_conedata<T> *cone)
/* Eliminate the infeasible rays with respect to  i  which
   are supposed to be consecutive from the head of the dd_Ray list,
   and sort the zero list assumed to be consecutive at the
   end of the list.
 */
{
  dd_rowrange fii,fiitest;
  T temp;
  dd_raydata<T> *Ptr;
  dd_raydata<T> *PrevPtr;
  dd_raydata<T> *NextPtr;
  dd_raydata<T> *ZeroPtr1;
  dd_raydata<T> *ZeroPtr0;
  bool found, completed, zerofound=false,negfound=false,posfound=false;

  cone->PosHead=nullptr;cone->ZeroHead=nullptr;cone->NegHead=nullptr;
  cone->PosLast=nullptr;cone->ZeroLast=nullptr;cone->NegLast=nullptr;

  /* Delete the infeasible rays  */
  PrevPtr= cone->ArtificialRay;
  Ptr = cone->FirstRay;
  if (PrevPtr->Next != Ptr)
    fprintf(stderr,"Error at dd_DeleteNegativeRays: ArtificialRay does not point the FirstRay.\n");
  completed=false;
  while (Ptr != nullptr && !completed){
/*    if ( (Ptr->ARay) < -zero ){ */
    if ( dd_Negative(Ptr->ARay)){
      dd_Eliminate(cone, &PrevPtr);
      Ptr=PrevPtr->Next;
    }
    else{
      completed=true;
    }
  }

  /* Sort the zero rays */
  Ptr = cone->FirstRay;
  cone->ZeroRayCount=0;
  while (Ptr != nullptr) {
    NextPtr=Ptr->Next;  /* remember the Next record */
    dd_set(temp,Ptr->ARay);
    if ( dd_Negative(temp)) {
      if (!negfound){
        fprintf(stderr,"Error: An infeasible ray found after their removal\n");
        negfound=true;
      }
    }
/*    else if (temp > zero){*/
    else if (dd_Positive(temp)){
      if (!posfound){
        posfound=true;
        cone->PosHead=Ptr;
        cone->PosLast=Ptr;
      }
      else{
        cone->PosLast=Ptr;
       }
    }
    else {
      (cone->ZeroRayCount)++;
      if (!zerofound){
        zerofound=true;
        cone->ZeroHead=Ptr;
        cone->ZeroLast=Ptr;
        cone->ZeroLast->Next=nullptr;
      }
      else{/* Find a right position to store the record sorted w.r.t. FirstInfeasIndex */
        fii=Ptr->FirstInfeasIndex;
        found=false;
        ZeroPtr1=nullptr;
        for (ZeroPtr0=cone->ZeroHead; !found && ZeroPtr0!=nullptr ; ZeroPtr0=ZeroPtr0->Next){
          fiitest=ZeroPtr0->FirstInfeasIndex;
          if (fiitest >= fii){
            found=true;
          }
          else ZeroPtr1=ZeroPtr0;
        }
        /* fprintf(stderr,"insert position found \n %d  index %ld\n",found, fiitest); */
        if (!found){           /* the new record must be stored at the end of list */
          cone->ZeroLast->Next=Ptr;
          cone->ZeroLast=Ptr;
          cone->ZeroLast->Next=nullptr;
        }
        else{
          if (ZeroPtr1==nullptr){ /* store the new one at the head, and update the head ptr */
            /* fprintf(stderr,"Insert at the head\n"); */
            Ptr->Next=cone->ZeroHead;
            cone->ZeroHead=Ptr;
          }
          else{                /* store the new one inbetween ZeroPtr1 and 0 */
            /* fprintf(stderr,"Insert inbetween\n");  */
            Ptr->Next=ZeroPtr1->Next;
            ZeroPtr1->Next=Ptr;
          }
        }
        /*
        Ptr->Next=cone->ZeroHead;
        cone->ZeroHead=Ptr;
        */
      }
    }
    Ptr=NextPtr;
  }
  /* joining the pos and zero lists */
  if (posfound){            /* -list empty & +list nonempty */
    cone->FirstRay=cone->PosHead;
    if (zerofound){              /* +list,0list nonempty */
      cone->PosLast->Next=cone->ZeroHead;
      cone->LastRay=cone->ZeroLast;
    }
    else{                        /* 0list empty & +list nonempty */
      cone->LastRay=cone->PosLast;
    }
  }
  else{                          /* +list empty & 0list nonempty */
    cone->FirstRay=cone->ZeroHead;
    cone->LastRay=cone->ZeroLast;
  }
  cone->ArtificialRay->Next=cone->FirstRay;
  cone->LastRay->Next=nullptr;
}

template<typename T>
void dd_FeasibilityIndices(long *fnum, long *infnum, dd_rowrange i, dd_conedata<T> *cone)
{
  /*Evaluate the number of feasible rays and infeasible rays*/
  /*  w.r.t the hyperplane  i*/
  dd_colrange j;
  T temp, tnext;
  dd_raydata<T> *Ptr;

  *fnum = 0;
  *infnum = 0;
  Ptr = cone->FirstRay;
  while (Ptr != nullptr) {
    temp=0;
    for (j = 0; j < cone->d; j++){
      dd_mul(tnext, cone->A[i - 1][j],Ptr->Ray[j]);
      dd_add(temp, temp, tnext);
    }
    if (dd_Nonnegative(temp))
      (*fnum)++;
    else
      (*infnum)++;
    Ptr = Ptr->Next;
  }
}


template<typename T>
bool dd_LexEqual(T *v1, T *v2, long dmax)
{ /* dmax is the size of vectors v1,v2 */
  bool determined, equal;
  dd_colrange j;

  equal = true;
  determined = false;
  j = 1;
  do {
    if (!dd_Equal(v1[j - 1],v2[j - 1])) {  /* 093c */
	equal = false;
        determined = true;
    } else {
      j++;
    }
  } while (!(determined) && (j <= dmax));
  return equal;
}

template<typename T>
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
  bool localdebug=false;

  dd_EvaluateARay1(hnew, cone);
   /*Check feasibility of rays w.r.t. hnew
     and put all infeasible ones consecutively */

  RayPtr0 = cone->ArtificialRay;   /*Pointer pointing RayPrt1*/
  RayPtr1 = cone->FirstRay;        /*1st hnew-infeasible ray to scan and compare with feasible rays*/
  dd_set(value1,cone->FirstRay->ARay);
  if (dd_Nonnegative(value1)) {
    if (cone->RayCount==cone->WeaklyFeasibleRayCount) cone->CompStatus=dd_AllFound;
    goto _L99;        /* Sicne there is no hnew-infeasible ray and nothing to do */
  }
  else {
    RayPtr2s = RayPtr1->Next;/* RayPtr2s must point the first feasible ray */
    pos2=1;
    while (RayPtr2s!=nullptr && dd_Negative(RayPtr2s->ARay)) {
      RayPtr2s = RayPtr2s->Next;
      pos2++;
    }
  }
  if (RayPtr2s==nullptr) {
    cone->FirstRay=nullptr;
    cone->ArtificialRay->Next=cone->FirstRay;
    cone->RayCount=0;
    cone->CompStatus=dd_AllFound;
    goto _L99;   /* All rays are infeasible, and the computation must stop */
  }
  RayPtr2 = RayPtr2s;   /*2nd feasible ray to scan and compare with 1st*/
  RayPtr3 = cone->LastRay;    /*Last feasible for scanning*/
  prevprogress=-10.0;
  pos1 = 1;
  completed=false;
  while ((RayPtr1 != RayPtr2s) && !completed) {
    dd_set(value1,RayPtr1->ARay);
    dd_set(value2,RayPtr2->ARay);
    dd_CheckEquality(cone->d, &RayPtr1, &RayPtr2, &equal);
    if ((dd_Positive(value1) && dd_Negative(value2)) || (dd_Negative(value1) && dd_Positive(value2))){
      dd_CheckAdjacency(cone, &RayPtr1, &RayPtr2, &adj);
      if (adj) dd_CreateNewRay(cone, RayPtr1, RayPtr2, hnew);
    }
    if (RayPtr2 != RayPtr3) {
      RayPtr2 = RayPtr2->Next;
      continue;
    }
    if (dd_Negative(value1) || equal) {
      dd_Eliminate(cone, &RayPtr0);
      RayPtr1 = RayPtr0->Next;
      RayPtr2 = RayPtr2s;
    } else {
      completed=true;
    }
    pos1++;
    progress = 100.0 * ((double)pos1 / pos2) * (2.0 * pos2 - pos1) / pos2;
    if (progress-prevprogress>=10 && pos1%10==0 && localdebug) {
      fprintf(stderr,"*Progress of iteration %5ld(/%ld):   %4ld/%4ld => %4.1f%% done\n",
	     cone->Iteration, cone->m, pos1, pos2, progress);
      prevprogress=progress;
    }
  }
  if (cone->RayCount==cone->WeaklyFeasibleRayCount) cone->CompStatus=dd_AllFound;
  _L99:;
}

template<typename T>
void dd_AddNewHalfspace2(dd_conedata<T> *cone, dd_rowrange hnew)
/* This procedure must be used under PreOrderedRun mode */
{
  bool localdebug=false;
  //  dd_raydata<T> *RayPtr0;
  dd_raydata<T> *RayPtr1;
  dd_raydata<T> *RayPtr2;

  dd_adjacencydata<T> *EdgePtr, *EdgePtr0;
  dd_rowrange fii1, fii2;

  dd_EvaluateARay2(hnew, cone);
   /* Check feasibility of rays w.r.t. hnew
      and sort them. ( -rays, +rays, 0rays)*/

  if (cone->PosHead==nullptr && cone->ZeroHead==nullptr) {
    cone->FirstRay=nullptr;
    cone->ArtificialRay->Next=cone->FirstRay;
    cone->RayCount=0;
    cone->CompStatus=dd_AllFound;
    goto _L99;   /* All rays are infeasible, and the computation must stop */
  }
  if (cone->ZeroHead==nullptr) cone->ZeroHead=cone->LastRay;

  EdgePtr=cone->Edges[cone->Iteration];
  while (EdgePtr != nullptr){
    RayPtr1=EdgePtr->Ray1;
    RayPtr2=EdgePtr->Ray2;
    fii1=RayPtr1->FirstInfeasIndex;
    dd_CreateNewRay(cone, RayPtr1, RayPtr2, hnew);
    fii2=cone->LastRay->FirstInfeasIndex;
    if (fii1 != fii2)
      dd_ConditionalAddEdge(cone,RayPtr1,cone->LastRay,cone->PosHead);
    EdgePtr0=EdgePtr;
    EdgePtr=EdgePtr->Next;
    delete EdgePtr0;
    (cone->EdgeCount)--;
  }
  cone->Edges[cone->Iteration]=nullptr;

  dd_DeleteNegativeRays(cone);

  set_addelem(cone->AddedHalfspaces, hnew);

  if (cone->Iteration<cone->m){
    if (cone->ZeroHead!=nullptr && cone->ZeroHead!=cone->LastRay){
      if (cone->ZeroRayCount>200 && localdebug) fprintf(stderr,"*New edges being scanned...\n");
      dd_UpdateEdges(cone, cone->ZeroHead, cone->LastRay);
    }
  }

  if (cone->RayCount==cone->WeaklyFeasibleRayCount) cone->CompStatus=dd_AllFound;
_L99:;
}


template<typename T>
void dd_SelectNextHalfspace0(dd_conedata<T> *cone, dd_rowset excluded, dd_rowrange *hnext)
{
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
  } while (!determined && i>=1);
  if (determined)
    *hnext = i;
  else
    *hnext = 0;
}

template<typename T>
void dd_SelectNextHalfspace1(dd_conedata<T> *cone, dd_rowset excluded, dd_rowrange *hnext)
{
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
  } while (!determined && i<=cone->m);
  if (determined)
    *hnext = i;
  else
    *hnext=0;
}

template<typename T>
void dd_SelectNextHalfspace2(dd_conedata<T> *cone, dd_rowset excluded, dd_rowrange *hnext)
{
  long i, fea, inf, infmin;   /*feasibility and infeasibility numbers*/

  infmin = cone->RayCount + 1;
  for (i = 1; i <= cone->m; i++) {
    if (!set_member(i, excluded)) {
      dd_FeasibilityIndices(&fea, &inf, i, cone);
      if (inf < infmin) {
	infmin = inf;
	*hnext = i;
      }
    }
  }
}

template<typename T>
void dd_SelectNextHalfspace3(dd_conedata<T> *cone, dd_rowset excluded, dd_rowrange *hnext)
{
  /*Choose the next hyperplane with maximum infeasibility*/
  long i, fea, inf, infmax;   /*feasibility and infeasibility numbers*/

  infmax = -1;
  for (i = 1; i <= cone->m; i++)
    if (!set_member(i, excluded)) {
      dd_FeasibilityIndices(&fea, &inf, i, cone);
      if (inf > infmax) {
	infmax = inf;
	*hnext = i;
      }
    }
}

template<typename T>
void dd_SelectNextHalfspace4(dd_conedata<T> *cone, dd_rowset excluded, dd_rowrange *hnext)
{
  long i, fea, inf, max, tmax, fi=0, infi=0;
  bool localdebug=false;

  max = -1;
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
        *hnext = i;
      }
    }
  }
  if (localdebug) {
    if (max == fi) {
      fprintf(stderr,"*infeasible rays (min) =%5ld, #feas rays =%5ld\n", infi, fi);
    } else {
      fprintf(stderr,"*infeasible rays (max) =%5ld, #feas rays =%5ld\n", infi, fi);
    }
  }
}

template<typename T>
void dd_SelectNextHalfspace5(dd_conedata<T> *cone, dd_rowset excluded, dd_rowrange *hnext)
{
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
        v1=v2;
      } else if (dd_LexSmaller(v2,v1,cone->d)) {
        minindex = i;
        v1=v2;
      }
    }
  }
  *hnext = minindex;
}


template<typename T>
void dd_SelectNextHalfspace6(dd_conedata<T> *cone, dd_rowset excluded, dd_rowrange *hnext)
{
  /*Choose the next hyperplane which is lexico-max*/
  long i, maxindex;
  T *v1, *v2;

  maxindex = 0;
  v1 = nullptr;
  for (i = 1; i <= cone->m; i++) {
    if (!set_member(i, excluded)) {
      v2= cone->A[i - 1];
      if (maxindex == 0) {
        maxindex = i;
        v1=v2;
      } else if (dd_LexLarger(v2, v1, cone->d)) {
        maxindex = i;
        v1=v2;
     }
    }
  }
  *hnext = maxindex;
}

template<typename T>
void dd_UniqueRows(dd_rowindex OV, long p, long r, T** A, long dmax, dd_rowset preferred, long *uniqrows)
{
 /* Select a subset of rows of A (with range [p, q] up to dimension dmax) by removing duplicates.
    When a row subset preferred is nonempty, those row indices in the set have priority.  If
    two priority rows define the same row vector, one is chosen.
    For a selected unique row i, OV[i] returns a new position of the unique row i.
    For other nonuniqu row i, OV[i] returns a negative of the original row j dominating i.
    Thus the original contents of OV[p..r] will be rewritten.  Other components remain the same.
    *uniqrows returns the number of unique rows.
*/
  long i,iuniq,j;
  T *x;

  if (p<=0 || p > r) goto _L99;
  iuniq=p; j=1;  /* the first unique row candidate */
  x=A[p-1];
  OV[p]=j;  /* tentative row index of the candidate */
  for (i=p+1;i<=r; i++){
    if (!dd_LexEqual(x,A[i-1],dmax)) {
      /* a new row vector found. */
      iuniq=i;
      j=j+1;
      OV[i]=j;    /* Tentatively register the row i.  */
      x=A[i-1];
    } else {
      /* rows are equal */
      if (set_member(i, preferred) && !set_member(iuniq, preferred)){
        OV[iuniq]=-i;  /* the row iuniq is dominated by the row i */
        iuniq=i;  /* the row i is preferred.  Change the candidate. */
        OV[i]=j;  /* the row i is tentatively registered. */
        x=A[i-1];
      } else {
        OV[i]=-iuniq;  /* the row iuniq is dominated by the row i */
      }
    }
  }
  *uniqrows=j;
  _L99:;
}



#ifndef RAND_MAX
#define RAND_MAX 32767
#endif


template<typename T>
void dd_ComputeRowOrderVector(dd_conedata<T> *cone)
{
  long i,itemp;

  cone->OrderVector[0]=0;
  switch (cone->HalfspaceOrder){
  case dd_MaxIndex:
    for(i=1; i<=cone->m; i++) cone->OrderVector[i]=cone->m-i+1;
    break;

  case dd_MinIndex:
    for(i=1; i<=cone->m; i++) cone->OrderVector[i]=i;
    break;

  case dd_LexMin: case dd_MinCutoff: case dd_MixCutoff: case dd_MaxCutoff:
    for(i=1; i<=cone->m; i++) cone->OrderVector[i]=i;
    dd_RandomPermutation(cone->OrderVector, cone->m, cone->rseed);
    dd_QuickSort(cone->OrderVector, 1, cone->m, cone->A, cone->d);
    break;

  case dd_LexMax:
    for(i=1; i<=cone->m; i++) cone->OrderVector[i]=i;
    dd_RandomPermutation(cone->OrderVector, cone->m, cone->rseed);
    dd_QuickSort(cone->OrderVector, 1, cone->m, cone->A, cone->d);
    for(i=1; i<=cone->m/2;i++){   /* just reverse the order */
      itemp=cone->OrderVector[i];
      cone->OrderVector[i]=cone->OrderVector[cone->m-i+1];
      cone->OrderVector[cone->m-i+1]=itemp;
    }
    break;

  case dd_RandomRow:
    for(i=1; i<=cone->m; i++) cone->OrderVector[i]=i;
    dd_RandomPermutation(cone->OrderVector, cone->m, cone->rseed);
    break;

  }
}

template<typename T>
void dd_UpdateRowOrderVector(dd_conedata<T> *cone, dd_rowset PriorityRows)
/* Update the RowOrder vector to shift selected rows
in highest order.
*/
{
  dd_rowrange i,j,k,j1=0,oj=0;
  long rr;
  bool found;

  found=true;
  rr=set_card(PriorityRows);
  for (i=1; i<=rr; i++){
    found=false;
    for (j=i; j<=cone->m && !found; j++){
      oj=cone->OrderVector[j];
      if (set_member(oj, PriorityRows)){
        found=true;
        j1=j;
      }
    }
    if (found){
      if (j1>i) {
        /* shift everything lower: ov[i]->cone->ov[i+1]..ov[j1-1]->cone->ov[j1] */
        for (k=j1; k>=i; k--) cone->OrderVector[k]=cone->OrderVector[k-1];
        cone->OrderVector[i]=oj;
      }
    } else {
      fprintf(stderr,"UpdateRowOrder: Error.\n");
      goto _L99;
    }
  }
_L99:;
}

template<typename T>
void dd_SelectPreorderedNext(dd_conedata<T> *cone, dd_rowset excluded, dd_rowrange *hh)
{
  dd_rowrange i,k;

  *hh=0;
  for (i=1; i<=cone->m && *hh==0; i++){
    k=cone->OrderVector[i];
    if (!set_member(k, excluded)) *hh=k ;
  }
}

template<typename T>
void dd_SelectNextHalfspace(dd_conedata<T> *cone, dd_rowset excluded, dd_rowrange *hh)
{
  if (cone->PreOrderedRun){
    dd_SelectPreorderedNext(cone, excluded, hh);
  }
  else {

    switch (cone->HalfspaceOrder) {

    case dd_MaxIndex:
      dd_SelectNextHalfspace0(cone, excluded, hh);
      break;

    case dd_MinIndex:
      dd_SelectNextHalfspace1(cone, excluded, hh);
      break;

    case dd_MinCutoff:
      dd_SelectNextHalfspace2(cone, excluded, hh);
      break;

    case dd_MaxCutoff:
      dd_SelectNextHalfspace3(cone, excluded, hh);
      break;

    case dd_MixCutoff:
      dd_SelectNextHalfspace4(cone, excluded, hh);
      break;

    case dd_LexMin: case dd_LexMax: case dd_RandomRow:
      dd_SelectNextHalfspace0(cone, excluded, hh);
      break;
    }
  }
}

/* end of cddcore.c */


/* cddlib.c: cdd library  (library version of cdd)
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.94g, March 23, 2012
   Standard ftp site: ftp.ifor.math.ethz.ch, Directory: pub/fukuda/cdd
*/

/* cdd : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron
   P= {x :  b - A x >= 0}.
   Please read COPYING (GNU General Public Licence) and
   the manual cddlibman.tex for detail.
*/

/*  This program is free software; you can redistribute it and/or modify
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

/* The first version C0.21 was created on November 10,1993
   with Dave Gillespie's p2c translator
   from the Pascal program pdd.p written by Komei Fukuda.
*/


template<typename T>
void dd_DDInit(dd_conedata<T> *cone)
{
  cone->Error=dd_NoError;
  cone->CompStatus=dd_InProgress;
  cone->RayCount = 0;
  cone->TotalRayCount = 0;
  cone->FeasibleRayCount = 0;
  cone->WeaklyFeasibleRayCount = 0;
  cone->EdgeCount=0; /* active edge count */
  cone->TotalEdgeCount=0; /* active edge count */
  dd_SetInequalitySets(cone);
  dd_ComputeRowOrderVector(cone);
  cone->RecomputeRowOrder=false;
}

template<typename T>
void dd_DDMain(dd_conedata<T> *cone)
{
  dd_rowrange hh, itemp, otemp;
  bool localdebug=false;

  if (cone->d<=0){
    cone->Iteration=cone->m;
    cone->FeasibleRayCount=0;
    cone->CompStatus=dd_AllFound;
    goto _L99;
  }
  while (cone->Iteration <= cone->m) {
    dd_SelectNextHalfspace(cone, cone->WeaklyAddedHalfspaces, &hh);
    if (set_member(hh,cone->NonequalitySet)){  /* Skip the row hh */
      if (localdebug) {
        fprintf(stderr,"*The row # %3ld should be inactive and thus skipped.\n", hh);
      }
      set_addelem(cone->WeaklyAddedHalfspaces, hh);
    } else {
      if (cone->PreOrderedRun)
        dd_AddNewHalfspace2(cone, hh);
      else{
        dd_AddNewHalfspace1(cone, hh);
      }
      set_addelem(cone->AddedHalfspaces, hh);
      set_addelem(cone->WeaklyAddedHalfspaces, hh);
    }
    if (!cone->PreOrderedRun){
      otemp=-400;
      for (itemp=1; cone->OrderVector[itemp]!=hh; itemp++)
        otemp=cone->OrderVector[cone->Iteration];
      cone->OrderVector[cone->Iteration]=hh;
        /* store the dynamic ordering in ordervec */
      cone->OrderVector[itemp]=otemp;
        /* store the dynamic ordering in ordervec */
    }
    if (localdebug) {
      fprintf(stderr,"(Iter, Row, #Total, #Curr, #Feas)= %5ld %5ld %9ld %6ld %6ld\n",
        cone->Iteration, hh, cone->TotalRayCount, cone->RayCount,
        cone->FeasibleRayCount);
    }
    if (cone->CompStatus==dd_AllFound||cone->CompStatus==dd_RegionEmpty) {
      set_addelem(cone->AddedHalfspaces, hh);
      goto _L99;
    }
    (cone->Iteration)++;
  }
  _L99:;
  if (cone->d<=0 || cone->newcol[1]==0){ /* fixing the number of output */
     cone->parent->n=cone->LinearityDim + cone->FeasibleRayCount -1;
     cone->parent->ldim=cone->LinearityDim - 1;
  } else {
    cone->parent->n=cone->LinearityDim + cone->FeasibleRayCount;
    cone->parent->ldim=cone->LinearityDim;
  }
}


template<typename T>
void dd_InitialDataSetup(dd_conedata<T> *cone)
{
  long j, r;
  dd_rowset ZSet;
  T* Vector1;
  T* Vector2;

  Vector1=new T[cone->d];
  Vector2=new T[cone->d];

  cone->RecomputeRowOrder=false;
  cone->ArtificialRay = nullptr;
  cone->FirstRay = nullptr;
  cone->LastRay = nullptr;
  set_initialize(&ZSet,cone->m);
  dd_AddArtificialRay(cone);
  set_copy(cone->AddedHalfspaces, cone->InitialHalfspaces);
  set_copy(cone->WeaklyAddedHalfspaces, cone->InitialHalfspaces);
  dd_UpdateRowOrderVector(cone, cone->InitialHalfspaces);
  for (r = 1; r <= cone->d; r++) {
    for (j = 0; j < cone->d; j++) {
      dd_set(Vector1[j], cone->B[j][r-1]);
      dd_neg(Vector2[j], cone->B[j][r-1]);
    }
    dd_ZeroIndexSet(cone->m, cone->d, cone->A, Vector1, ZSet);
    if (set_subset(cone->EqualitySet, ZSet)) {
      dd_AddRay(cone, Vector1);
      if (cone->InitialRayIndex[r]==0)
        dd_AddRay(cone, Vector2);
    }
  }
  dd_CreateInitialEdges(cone);
  cone->Iteration = cone->d + 1;
  if (cone->Iteration > cone->m) cone->CompStatus=dd_AllFound; /* 0.94b  */
  set_free(ZSet);
  delete [] Vector1;
  delete [] Vector2;
}

template<typename T>
bool dd_CheckEmptiness(dd_polyhedradata<T> *poly, dd_ErrorType *err)
{
  dd_rowset R, S;
  dd_matrixdata<T> *M=nullptr;
  bool answer=false;

  *err=dd_NoError;

  if (poly->representation==dd_Inequality){
    M=dd_CopyInequalities(poly);
	set_initialize(&R, M->rowsize);
	set_initialize(&S, M->rowsize);
	if (!dd_ExistsRestrictedFace(M, R, S, err)){
	  poly->child->CompStatus=dd_AllFound;
	  poly->IsEmpty=true;
	  poly->n=0;
	  answer=true;
	}
	set_free(R);
	set_free(S);
	dd_FreeMatrix(M);
  } else if (poly->representation==dd_Generator && poly->m<=0){
	*err=dd_EmptyVrepresentation;
	poly->IsEmpty=true;
	poly->child->CompStatus=dd_AllFound;
	answer=true;
	poly->child->Error=*err;
  }

  return answer;
}


template<typename T>
bool dd_DoubleDescription(dd_polyhedradata<T> *poly, dd_ErrorType *err)
{
  dd_conedata<T> *cone=nullptr;
  bool found=false;

  *err=dd_NoError;

  if (poly!=nullptr && (poly->child==nullptr || poly->child->CompStatus!=dd_AllFound)){
    cone=dd_ConeDataLoad(poly);
    /* create a cone associated with poly by homogenization */
    dd_DDInit(cone);
    if (poly->representation==dd_Generator && poly->m<=0){
       *err=dd_EmptyVrepresentation;
       cone->Error=*err;
       goto _L99;
    }
    /* Check emptiness of the polyhedron */
    dd_CheckEmptiness(poly,err);

    if (cone->CompStatus!=dd_AllFound){
      dd_FindInitialRays(cone, &found);
      if (found) {
	dd_InitialDataSetup(cone);
	if (cone->CompStatus==dd_AllFound) goto _L99;
	dd_DDMain(cone);
	if (cone->FeasibleRayCount!=cone->RayCount) *err=dd_NumericallyInconsistent; /* cddlib-093d */
      }
    }
  }
_L99: ;
  return found;
}

template<typename T>
bool dd_DoubleDescription2(dd_polyhedradata<T> *poly, dd_RowOrderType horder, dd_ErrorType *err)
{
  dd_conedata<T> *cone=nullptr;
  bool found=false;

  *err=dd_NoError;

  if (poly!=nullptr && (poly->child==nullptr || poly->child->CompStatus!=dd_AllFound)){
    cone=dd_ConeDataLoad(poly);
    // create a cone associated with poly by homogenization
    cone->HalfspaceOrder=horder;  // set the row order
    dd_DDInit(cone);
    if (poly->representation==dd_Generator && poly->m<=0){
      *err=dd_EmptyVrepresentation;
      cone->Error=*err;
      goto _L99;
    }
    // Check emptiness of the polyhedron
    dd_CheckEmptiness(poly,err);

    if (cone->CompStatus!=dd_AllFound){
      dd_FindInitialRays(cone, &found);
      if (found) {
        dd_InitialDataSetup(cone);
        if (cone->CompStatus==dd_AllFound) goto _L99;
        dd_DDMain(cone);
        if (cone->FeasibleRayCount!=cone->RayCount) *err=dd_NumericallyInconsistent;
      }
    }
  }
_L99: ;
  return found;
}

template<typename T>
bool dd_DDInputAppend(dd_polyhedradata<T> **poly, dd_matrixdata<T> *M,
			    dd_ErrorType *err)
{
  /* This is imcomplete.  It simply solves the problem from scratch.  */
  bool found;
  dd_ErrorType error;

  if ((*poly)->child!=nullptr) dd_FreeDDMemory(*poly);
  dd_AppendMatrix2Poly(poly, M);
  (*poly)->representation=dd_Inequality;
  found=dd_DoubleDescription(*poly, &error);
  *err=error;
  return found;
}



template<typename T>
dd_matrixdata<T> *MyMatrix_PolyFile2Matrix(MyMatrix<T> const&TheEXT)
{
  dd_matrixdata<T> *M=nullptr;
  dd_rowrange m_input, i;
  dd_colrange d_input, j;
  dd_RepresentationType rep;
  bool localdebug=false;
  T value;

  m_input=TheEXT.rows();
  d_input=TheEXT.cols();

  rep=dd_Generator; /* using dd_Inequality led to horrible bugs */
  /*  NT=dd_dd_Integer;*/
  M=dd_CreateMatrix<T>(m_input, d_input);
  M->representation=rep;

  for (i = 1; i <= m_input; i++)
    for (j = 1; j <= d_input; j++)
      {
	value=TheEXT(i-1, j-1);
	M->matrix[i-1][j - 1]=value;
        if (localdebug) std::cerr << "i=" << i << " j=" << j << " value=" << value << "\n";
      }
  return M;
}


template<typename T>
MyMatrix<T> FAC_from_poly(dd_polyhedradata<T> const *poly, int const& nbCol)
{
  dd_raydata<T>* RayPtr = poly->child->FirstRay;
  int nbRay=0;
  while (RayPtr != nullptr) {
    if (RayPtr->feasible)
      nbRay++;
    RayPtr = RayPtr->Next;
  }
  MyMatrix<T> TheFAC=MyMatrix<T>(nbRay, nbCol);
  int iRay=0;
  RayPtr = poly->child->FirstRay;
  while (RayPtr != nullptr) {
    if (RayPtr->feasible) {
      for (int iCol=0; iCol<nbCol; iCol++)
        TheFAC(iRay, iCol)=RayPtr->Ray[iCol];
      iRay++;
    }
    RayPtr = RayPtr->Next;
  }
  return TheFAC;
}



template<typename T>
std::vector<Face> ListIncd_from_poly(dd_polyhedradata<T> const *poly, MyMatrix<T> const& EXT)
{
  std::vector<Face> ListIncd;
  size_t nbCol=EXT.cols();
  size_t nbRow=EXT.rows();
  MyVector<T> eFac(nbCol);
  dd_raydata<T>* RayPtr = poly->child->FirstRay;
  while (RayPtr != nullptr) {
    if (RayPtr->feasible) {
      for (size_t iCol=0; iCol<nbCol; iCol++)
	eFac(iCol)=RayPtr->Ray[iCol];
      MyVector<T> ListScal=ListScalarProduct(eFac, EXT);
      Face V(nbRow);
      for (size_t iRow=0; iRow<nbRow; iRow++) {
	T eVal=ListScal(iRow);
	if (eVal == 0)
	  V[iRow]=1;
      }
      ListIncd.push_back(V);
    }
    RayPtr = RayPtr->Next;
  }
  return ListIncd;
}



template<typename T>
std::vector<int> RedundancyReductionClarkson(MyMatrix<T> const&TheEXT)
{
  dd_ErrorType err;
  int nbRow=TheEXT.rows();
  std::cerr << "TheEXT=\n";
  WriteMatrix(std::cout, TheEXT);
  dd_matrixdata<T>* M=MyMatrix_PolyFile2Matrix(TheEXT);
  T smallVal = 0;
  dd_rowset rows = dd_RedundantRowsViaShooting(M, &err);
  std::vector<int> ListIdx;
  for (int i_row=0; i_row<nbRow; i_row++) {
    bool isin = set_member(i_row+1, rows);
    if (!isin) ListIdx.push_back(i_row);
  }
  dd_FreeMatrix(M);
  set_free(rows);
  return ListIdx;
}



template<typename T>
MyMatrix<T> DualDescription(MyMatrix<T> const&TheEXT)
{
  dd_polyhedradata<T> *poly;
  dd_matrixdata<T> *M;
  dd_ErrorType err;
  int nbCol=TheEXT.cols();
  M=MyMatrix_PolyFile2Matrix(TheEXT);
  poly=dd_DDMatrix2Poly(M, &err);
  MyMatrix<T> TheFAC=FAC_from_poly(poly, nbCol);
  dd_FreePolyhedra(poly);
  dd_FreeMatrix(M);
  return TheFAC;
}


template<typename T>
std::vector<Face> DualDescription_incd(MyMatrix<T> const&TheEXT)
{
  dd_polyhedradata<T> *poly;
  dd_matrixdata<T> *M;
  dd_ErrorType err;
  M=MyMatrix_PolyFile2Matrix(TheEXT);
  poly=dd_DDMatrix2Poly(M, &err);
  std::vector<Face> ListIncd=ListIncd_from_poly(poly, TheEXT);
  dd_FreePolyhedra(poly);
  dd_FreeMatrix(M);
  return ListIncd;
}

}

#endif

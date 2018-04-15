#if !defined DUALDESC_INCLUDE_SYMPOL
#include "MPQ_Matrix.h"
#include "GroupFct.h"
#include "Heuristic_fct.h"
#include "MPI_functions.h"
#include "MPQ_CDD_LRS.h"
#include "PolytopeFct.h"
#include "PolytopeEquiStab.h"
#include <vector>
#include "gmpxx.h"
using namespace std;

typedef vector<vector<int> > VectVectInt;

typedef struct {
  vector<vector<int> > ListRepresentative;
  vector<int> ListStatus;
  vector<int> ListSplittable;
  vector<mpz_class> OrbitSizes;
  int MaxAllowedUndone;
  int nbVert;
  mpz_t GRPsize;
} EnumerationOrbit;

void VectVectInt_Magma_Print(FILE *f, VectVectInt &ListOrbit);
void DUALDESC_MPI_SelectAndFlipComputation(MyMatrix *EXT, 
					   vector<int> eSelect, 
					   TheGroupFormat &GRP, 
					   TheHeuristic *TheHeuSplitting, 
					   int TheLevel, 
					   int NbProc, 
					   int *ListProc, 
					   VectVectInt &TheOutput, 
					   MPI_Comm comm);
void DUALDESC_MPI_AwaitingOrders(TheHeuristic *TheHeuSplitting, 
				 MPI_Comm comm);
void DUALDESC_MPI_SelectComputation(MyMatrix *EXT, 
				    TheGroupFormat &GRP, 
				    TheHeuristic *TheHeuSplitting, 
				    int TheLevel, 
				    int NbProc, 
				    int *ListProc, 
				    VectVectInt &TheOutput, 
				    MPI_Comm comm);
#endif
#define DUALDESC_INCLUDE_SYMPOL

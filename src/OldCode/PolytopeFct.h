#if !defined POLYTOPE_INCLUDE_SYMPOL
#include "MPQ_Matrix.h"
#include "GroupFct.h"
#include "gmpxx.h"
#include <vector>
#include <list>

using namespace std;

typedef struct {
  int nbDirect;
  int nbDual;
  mpq_t OptimalValue;
  int *DualSolutionPos;
  mpq_t *DualSolutionVal;
  int *DirectSolutionPos;
  mpq_t *DirectSolutionVal;
  int PrimalDefined;
  int DualDefined;
} LpSolution;


void ReadPolytope(FILE *f, MyMatrix *TheMat);
void WritePolytope(FILE *f, MyMatrix *TheMat);
void MatrixReduction(MyMatrix *TheMat, vector<int> SelectedRows, MyMatrix *TheMatRed);
set<vector<long> > GetCompressedFromPolytope(MyMatrix *TheEXT, MyMatrix *TheFAC);
void OrbitSplittingSet(set<vector<long> > ListTotal, 
		       TheGroupFormat &TheGRP, 
		       vector<vector<int> > &TheReturn);
void OrbitSplitting(MyMatrix *TheEXT, MyMatrix *TheFAC, 
		    TheGroupFormat &TheGRP, 
		    vector<vector<int> > &TheReturn);
vector<int> ComputeFlipping(MyMatrix *TheEXT, vector<int> OneInc, vector<int> sInc);
void PrintListOrbit(FILE *f, vector<vector<int> > &ListOrbit);
void TestFacetness(MyMatrix *TheEXT, vector<int> OneInc);
void CompressedToExpanded(vector<long>& TheCompressed, vector<int>& TheExpanded, long& nbExt);
void ExpandedToCompressed(vector<long>& TheCompressed, vector<int>& TheExpanded, long& nbExt);
#endif
#define POLYTOPE_INCLUDE_SYMPOL

#if !defined POLYTOPEEQUISTAB_INCLUDE_SYMPOL
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <vector>

#include "gmpxx.h"
#include "GroupFct.h"
#include "MPQ_Matrix.h"

#include "defs.hh"
#include "graph.hh"
#include "partition.hh"
#include "timer.hh"
#include "utils.hh"

using namespace std;
using namespace bliss;

typedef struct{
  int nbVertices;
  int hSize;
  int *AdjMat;
  int *ListColor;
} TheRecGraph;

typedef struct {
  int nbRow;
  int *TheMat;
  vector<mpq_class> ListWeight;
} WeightMatrix;

void GetGraphAutomorphismGroup(Graph *g, unsigned int nof_vertices, TheGroupFormat &GRP);
Graph* ReadGraphFromFile(FILE *f, unsigned int &nof_vertices);
void DefineWeightMatrix(WeightMatrix &WMat, int nbRow);
void ReadWeightedMatrix(FILE *f, WeightMatrix &WMat);
void PrintWeightedMatrix(FILE *f, WeightMatrix &WMat);
vector<int> GetLocalInvariantWeightMatrix(WeightMatrix &WMat, vector<int> eSet);
void GetWeightMatrix(WeightMatrix &WMat, MyMatrix *TheEXT);
void GetStabilizerWeightMatrix(WeightMatrix WMat, TheGroupFormat &GRP);
void TestEquivalenceWeightMatrix(WeightMatrix WMat1, WeightMatrix &WMat2, int *TheReply, vector<int> & TheEquiv);
#endif
#define POLYTOPEEQUISTAB_INCLUDE_SYMPOL

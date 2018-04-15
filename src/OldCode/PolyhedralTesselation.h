#if !defined POLYHEDRAL_TESSELATION_INCLUDE

#include "Matrix.h"
#include "PolytopeFct.h"
#include "gmpxx.h"
#include <vector>
#include <list>
using namespace std;

typedef struct {
  int iRank;
  int lenBoundary;
  int *ListSign;
  int *ListOrbit;
  MyMatrix *ListElt;
} TheBoundary;

typedef struct {
  MyMatrix *EXT;
  MyVector *IsoEXT;
  list<list<int> > ListIncd;
  /* The thing below is computed by the program */
  MyMatrix *TheSpann;
  MyMatrix *EXTexpSpann;
  MyMatrix *FACspann;
} OneCell;

typedef struct {
  int nbFac;
  MyMatrix *FAC;
  int *ListAdjacentDomain;
  MyMatrixPtr *ListMatrix;
} AdjacencyStructure;

typedef struct {
  int TheRank;
  int n;
  int *ListNbOrbit;
  OneCell **ListCells;
  TheBoundary **ListBoundary;
  AdjacencyStructure *ListAdjStruct;
  MyMatrix *QuadForm;
} TotalInformationTesselation;

typedef struct {
  int TheStatus;
  int iDomain;
  MyMatrix *TheMat;
  MyVector *TheVect;
} SingleConf;

typedef SingleConf* SingleConfPtr;

/* By complex face, we mean a face that is image under the
   Hecke operator */
typedef struct {
  int iRank;
  int iOrbit;
  MyMatrix *eActMat;
} ImageFace;

#endif
#define POLYHEDRAL_TESSELATION_INCLUDE

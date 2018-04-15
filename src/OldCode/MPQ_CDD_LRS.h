#if !defined MPQ_CDD_LRS_INCLUDE_SYMPOL
using namespace std;
#include <vector>
#include "MPQ_Matrix.h"
#include "PolytopeFct.h"

#include "setoper.h"
#include "cdd.h"
#include "lrslib.h"

void LPSOL_Free(LpSolution *eSol);
void CDD_LinearProgramming(MyMatrix *TheEXT, MyVector *eVect, LpSolution *eSol);
void CDD_DualDescription(MyMatrix *TheEXT, MyMatrix *TheFAC);
void LRS_DualDescription(MyMatrix *TheEXT, MyMatrix *TheFAC);
void FindOneInitialVertex(MyMatrix *EXT, vector<int> &eInc);
void Initialize_All();
void Free_All();
#endif
#define MPQ_CDD_LRS_INCLUDE_SYMPOL

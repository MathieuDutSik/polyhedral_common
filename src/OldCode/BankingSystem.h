#if !defined BANKING_INCLUDE_SYMPOL
#include "MPQ_Matrix.h"
#include "MPI_functions.h"
#include "PolytopeEquiStab.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "vector"

using namespace std;

typedef vector<vector<int> > VectVectInt;

typedef struct {
  vector<MyMatrixPtr> ListEXT;
  vector<VectVectInt> ListDualDesc;
} BankingDescription;

void BANK_OperatingPolyhedralBank(BankingDescription &TheBank, 
				  MPI_Comm comm);

#endif
#define BANKING_INCLUDE_SYMPOL

#if !defined LATTICE_ISOM_INCLUDE_SYMPOL
#include "MPQ_Matrix.h"
#include "MatrixOperation.h"
#include "PolytopeEquiStab.h"
#include "EnumerationConfig.h"
void DOUBL_GetGramMatrixAutomorphismGroup(MatDOUBL eMat, double TheTol, TheGroupFormat &GRPperm, vector<MatINT> &ListMatrGens);
void DOUBL_TestGramMatrixEquivalence(MatDOUBL eMat1, MatDOUBL eMat2, double TheTol, int &test);
#endif
#define LATTICE_ISOM_INCLUDE_SYMPOL

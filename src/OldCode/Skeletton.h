#if !defined SKELETTON_INCLUDE_SYMPOL
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <vector>

#include "MPQ_Matrix.h"
#include "MPQ_CDD_LRS.h"
#include "PolytopeEquiStab.h"


list<vector<int> > SPAN_face_LinearProgramming(vector<int> face,
					       TheGroupFormat StabFace,
					       MyMatrix* FAC,
					       TheGroupFormat FullGRP);
int TestPositiveRelationSimple(MyMatrix *ListVect);
void SearchPositiveRelationSimple(MyMatrix *ListVect, int *eTestExist, MyVector *TheRelat);
vector<vector<vector<int> > > EnumerationFaces(TheGroupFormat TheGRP, MyMatrix* FAC, int LevSearch);
vector<vector<int> > DoMemoryEfficientEnum(TheGroupFormat TheGRP, 
					   MyMatrix* FAC, int LevSearch);
void PrintListListOrb_IntGAP(FILE *f, vector<vector<vector<int> > > ListListOrb);
void PrintListOrb_GAP(FILE *f, vector<vector<int> > ListOrb);
vector<vector<vector<int> > > DoMemoryEfficientAllEnum(TheGroupFormat TheGRP, MyMatrix* FAC, 
						       int LevSearch);
#endif
#define SKELETTON_INCLUDE_SYMPOL

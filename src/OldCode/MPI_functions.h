#if !defined MPI_INCLUDE_SYMPOL
#include "GroupFct.h"
#include "MPQ_Matrix.h"
#include "PolytopeFct.h"
#include <stdlib.h>
#include <mpi.h>
#include <vector>
#include "gmpxx.h"

using namespace std;

typedef vector<mpq_class> VectGMP;

void MPI_SEND_mpq_vector(VectGMP &a, int dest, int tag, MPI_Comm comm);
void MPI_RECV_mpq_vector(VectGMP &a, int src, int tag, MPI_Comm comm);
void MPI_SEND_MyMatrix(MyMatrix *eMat, int dest, int tag, MPI_Comm comm);
void MPI_RECV_MyMatrix(MyMatrix *eMat, int src, int tag, MPI_Comm comm);
void MPI_SEND_Group(TheGroupFormat &GRP, int dest, int tag, MPI_Comm comm);
void MPI_RECV_Group(TheGroupFormat &GRP, int src, int tag, MPI_Comm comm);
void MPI_SEND_VectVectInt(vector<vector<int> > &eListListI, int dest, int tag, MPI_Comm comm);
void MPI_RECV_VectVectInt(vector<vector<int> > &eListListI, int src, int tag, MPI_Comm comm);
void MPI_SEND_VectInt(vector<int> &eListI, int dest, int tag, MPI_Comm comm);
void MPI_RECV_VectInt(vector<int> &eListI,  int src, int tag, MPI_Comm comm);
#endif
#define MPI_INCLUDE_SYMPOL

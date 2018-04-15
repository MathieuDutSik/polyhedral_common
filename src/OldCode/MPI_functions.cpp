#include "GroupFct.h"
#include "MPQ_Matrix.h"
#include "MPI_functions.h"



void MPI_SEND_mpq_vector(VectGMP &a, int dest, int tag, MPI_Comm comm)
{
  vector<int> eVectSend;
  mpz_t eNum, eDen, n_z, d_z, q_z, r_z;
  int *eVectSendC;
  int eRes, i, idx, eSign;
  int eSize, eVal, ret;
  eRes=147;
  mpz_init(n_z);
  mpz_init(q_z);
  mpz_init(r_z);
  mpz_init(eNum);
  mpz_init(eDen);
  mpz_init(d_z);
  mpz_set_ui(d_z, eRes);
  eSize=a.size();
  eVectSend.push_back(eSize);
  for (i=0; i<eSize; i++)
    {
      eSign=mpq_sgn(a[i].get_mpq_t());
      eVectSend.push_back(eSign);
      if (eSign != 0)
	{
	  mpq_get_num(eNum, a[i].get_mpq_t());
	  mpq_get_den(eDen, a[i].get_mpq_t());
	  if (eSign == -1)
	    mpz_neg(eNum, eNum);
	  mpz_set(n_z, eNum);
	  while(1)
	    {
	      mpz_fdiv_qr(q_z, r_z, n_z, d_z);
	      eVal=mpz_get_ui(r_z);
	      eVectSend.push_back(eVal);
	      eSign=mpz_sgn(q_z);
	      if (eSign == 0)
		{
		  eVectSend.push_back(-1);
		  break;
		}
	      mpz_set(n_z, q_z);
	    }
	  mpz_set(n_z, eDen);
	  while(1)
	    {
	      mpz_fdiv_qr(q_z, r_z, n_z, d_z);
	      eVal=mpz_get_ui(r_z);
	      eVectSend.push_back(eVal);
	      eSign=mpz_sgn(q_z);
	      if (eSign == 0)
		{
		  eVectSend.push_back(-1);
		  break;
		}
	      mpz_set(n_z, q_z);
	    }
	}
    }
  if ((eVectSendC = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  eSize=eVectSend.size();
  eVectSendC[0]=eSize;
  ret=MPI_Send(eVectSendC, 1, MPI_INT, dest, tag, comm);
  free(eVectSendC);
  //
  if ((eVectSendC = (int*)malloc(eVectSend.size()*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  for (idx=0; idx<eSize; idx++)
    eVectSendC[idx]=eVectSend[idx];
  
  ret=MPI_Send(eVectSendC, eSize, MPI_INT, dest, 1002, comm);
  free(eVectSendC);
  mpz_clear(n_z);
  mpz_clear(q_z);
  mpz_clear(r_z);
  mpz_clear(eNum);
  mpz_clear(eDen);
  mpz_clear(d_z);
}


void MPI_RECV_mpq_vector(VectGMP &a, int src, int tag, MPI_Comm comm)
{
  int *eVectRecvC;
  mpq_class tcla;
  mpq_t t_q1, t_q2, t_q3;
  mpz_t eNum, eDen, eExpo, d_z, t;
  int eVal, eSize, siz, nbEnt;
  MPI_Status status;
  int eRes, idx, eSign, idxVect;
  if ((eVectRecvC = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  MPI_Recv(eVectRecvC, 1, MPI_INT, src, tag, comm, &status);
  eSize=eVectRecvC[0];
  free(eVectRecvC);

  if ((eVectRecvC = (int*)malloc(eSize*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  MPI_Recv(eVectRecvC, eSize, MPI_INT, src, 1002, comm, &status);
  eRes=147;
  mpq_init(t_q1);
  mpq_init(t_q2);
  mpq_init(t_q3);
  mpz_init(eNum);
  mpz_init(eDen);
  mpz_init(eExpo);
  mpz_init(d_z);
  mpz_init(t);
  mpz_set_ui(d_z, eRes);
  idx=0;
  nbEnt=eVectRecvC[idx];
  idx++;
  idxVect=0;
  while(1)
    {
      eSign=eVectRecvC[idx];
      idx++;
      if (eSign == 0)
	{
	  mpq_set_ui(t_q3, 0, 1);
	  tcla=mpq_class(t_q3);
	  a.push_back(tcla);
	}
      else
	{
	  mpz_set_ui(eNum,0);
	  mpz_set_ui(eExpo,1);
	  while(1)
	    {
	      eVal=eVectRecvC[idx];
	      idx++;
	      if (eVal == -1)
		break;
	      mpz_set_si(t, eVal);
	      mpz_mul(t, t, eExpo);
	      mpz_mul(eExpo, eExpo, d_z);
	      mpz_add(eNum, eNum, t);
	    }
	  mpz_set_ui(eDen,0);
	  mpz_set_ui(eExpo,1);
	  while(1)
	    {
	      eVal=eVectRecvC[idx];
	      idx++;
	      if (eVal == -1)
		break;
	      mpz_set_ui(t, eVal);
	      mpz_mul(t, t, eExpo);
	      mpz_mul(eExpo, eExpo, d_z);
	      mpz_add(eDen, eDen, t);
	    }
	  mpq_set_z(t_q1,eNum);
	  if (eSign == -1)
	    mpq_neg(t_q1, t_q1);
	  mpq_set_z(t_q2,eDen);
	  mpq_div(t_q3, t_q1, t_q2);
	  tcla=mpq_class(t_q3);
	  a.push_back(tcla);
	}
      idxVect++;
      if (idx == eSize)
	break;
    }
  siz=a.size();
  if (siz != nbEnt)
    {
      fprintf(stderr, "Error in the number of entries\n");
      fprintf(stderr, "siz=%d nbEnt=%d\n", siz, nbEnt);
      exit(1);
    }
  free(eVectRecvC);
  mpq_clear(t_q1);
  mpq_clear(t_q2);
  mpq_clear(t_q3);
  mpz_clear(eNum);
  mpz_clear(eDen);
  mpz_clear(eExpo);
  mpz_clear(d_z);
  mpz_clear(t);
}


void MPI_SEND_MyMatrix(MyMatrix *eMat, int dest, int tag, MPI_Comm comm)
{
  vector<mpq_class> eVect;
  int *eVectSendC;
  int nbRow, nbCol, ret, nbEnt, i;
  mpq_class tcla;
  if ((eVectSendC = (int*)malloc(2*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  nbRow=eMat->nbRow;
  nbCol=eMat->nbCol;
  eVectSendC[0]=nbRow;
  eVectSendC[1]=nbCol;
  ret=MPI_Send(eVectSendC, 2, MPI_INT, dest, 1292, comm);
  free(eVectSendC);

  nbEnt=nbRow*nbCol;
  for (i=0; i<nbEnt; i++)
    {
      tcla=mpq_class(eMat->ListElt[i]);
      eVect.push_back(tcla);
    }
  MPI_SEND_mpq_vector(eVect, dest, tag, comm);
}



void MPI_RECV_MyMatrix(MyMatrix *eMat, int src, int tag, MPI_Comm comm)
{
  VectGMP eVect;
  int *eVectRecvC;
  MPI_Status status;
  int nbRow, nbCol, nbEnt, i;
  if ((eVectRecvC = (int*)malloc(2*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  MPI_Recv(eVectRecvC, 2, MPI_INT, src, 1292, comm, &status);
  nbRow=eVectRecvC[0];
  nbCol=eVectRecvC[1];
  free(eVectRecvC);

  MATRIX_Allocate(eMat, nbRow, nbCol);
  MPI_RECV_mpq_vector(eVect, src, tag, comm);
  nbEnt=nbRow*nbCol;
  for (i=0; i<nbEnt; i++)
    mpq_set(eMat->ListElt[i], eVect[i].get_mpq_t());
}



void MPI_SEND_Group(TheGroupFormat &GRP, int dest, int tag, MPI_Comm comm)
{
  list<Permutation::ptr> ListGen;
  list<Permutation::ptr>::iterator iter;
  int n, nbGen, h, idx, i, ret, eVal;
  int *eVectSendC;
  ListGen=GRP.group->S;
  if ((eVectSendC = (int*)malloc(2*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  nbGen=ListGen.size();
  n=GRP.n;
  eVectSendC[0]=nbGen;
  eVectSendC[1]=n;
  ret=MPI_Send(eVectSendC, 2, MPI_INT, dest, tag, comm);
  free(eVectSendC);
  
  h=n*nbGen;
  if ((eVectSendC = (int*)malloc(h*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  idx=0;
  iter=ListGen.begin();
  while(iter != ListGen.end())
    {
      for (i=0; i<n; i++)
	{
	  eVal=(*iter)->at(i);
	  eVectSendC[idx]=eVal;
	  idx++;
	}
      iter++;
    }
  ret=MPI_Send(eVectSendC, h, MPI_INT, dest, tag+4, comm);
  free(eVectSendC);
}


void MPI_RECV_Group(TheGroupFormat &GRP, int src, int tag, MPI_Comm comm)
{
  list<Permutation::ptr> generatorList;
  int i, n, h, iGen, nbGen, idx;
  int *eVectRecvC;
  int *g;
  MPI_Status status;
  if ((eVectRecvC = (int*)malloc(2*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  MPI_Recv(eVectRecvC, 2, MPI_INT, src, tag, comm, &status);
  nbGen=eVectRecvC[0];
  n=eVectRecvC[1];
  free(eVectRecvC);

  h=n*nbGen;
  if ((eVectRecvC = (int*)malloc(h*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  MPI_Recv(eVectRecvC, h, MPI_INT, src, tag+4, comm, &status);
  if ((g = (int*)malloc(n*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  idx=0;
  for (iGen=0; iGen<nbGen; iGen++)
    {
      for (i=0; i<n; i++)
	{
	  g[i]=eVectRecvC[idx];
	  idx++;
	}
      vector<dom_int> v(g, g+n);
      generatorList.push_back(Permutation::ptr(new Permutation(v)));
    }
  free(eVectRecvC);
  free(g);
  GRP.n=n;
  GRP.group=construct(n, generatorList.begin(), generatorList.end());
}

void MPI_SEND_VectVectInt(vector<vector<int> > &eListListI, int dest, int tag, MPI_Comm comm)
{
  int iCase, nbCase, nbNeed, i, siz, idx, ret;
  int *eVectSendC;
  nbCase=eListListI.size();
  nbNeed=1+nbCase;
  for (iCase=0; iCase<nbCase; iCase++)
    {
      siz=eListListI[iCase].size();
      nbNeed=nbNeed + siz;
    }
  if ((eVectSendC = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  eVectSendC[0]=nbNeed;
  ret=MPI_Send(eVectSendC, 1, MPI_INT, dest, tag, comm);
  free(eVectSendC);

  if ((eVectSendC = (int*)malloc(nbNeed*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  idx=0;
  eVectSendC[idx]=nbCase;
  idx++;
  for (iCase=0; iCase<nbCase; iCase++)
    {
      siz=eListListI[iCase].size();
      eVectSendC[idx]=siz;
      idx++;
      for (i=0; i<siz; i++)
	{
	  eVectSendC[idx]=eListListI[iCase][i];
	  idx++;
	}
    }
  ret=MPI_Send(eVectSendC, nbNeed, MPI_INT, dest, tag+5, comm);
  free(eVectSendC);
}

void MPI_RECV_VectVectInt(vector<vector<int> > &eListListI, int src, int tag, MPI_Comm comm)
{
  vector<int> eListI;
  int *eVectRecvC;
  int idx, iCase, nbCase, i, siz, eVal;
  MPI_Status status;
  int nbNeed;
  eListListI.clear();
  if ((eVectRecvC = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  MPI_Recv(eVectRecvC, 1, MPI_INT, src, tag, comm, &status);
  nbNeed=eVectRecvC[0];
  free(eVectRecvC);
  
  if ((eVectRecvC = (int*)malloc(nbNeed*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  MPI_Recv(eVectRecvC, nbNeed, MPI_INT, src, tag+5, comm, &status);
  idx=0;
  nbCase=eVectRecvC[idx];
  idx++;
  for (iCase=0; iCase<nbCase; iCase++)
    {
      siz=eVectRecvC[idx];
      idx++;
      eListI.clear();
      for (i=0; i<siz; i++)
	{
	  eVal=eVectRecvC[idx];
	  idx++;
	  eListI.push_back(eVal);
	}
      eListListI.push_back(eListI);
    }
  free(eVectRecvC);
}



void MPI_SEND_VectInt(vector<int> &eListI, int dest, int tag, MPI_Comm comm)
{
  int len, i, ret;
  int *eVectSendC;
  len=eListI.size();
  if ((eVectSendC = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  eVectSendC[0]=len;
  ret=MPI_Send(eVectSendC, 1, MPI_INT, dest, tag, comm);
  free(eVectSendC);

  if ((eVectSendC = (int*)malloc(len*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  for (i=0; i<len; i++)
    eVectSendC[i]=eListI[i];
  ret=MPI_Send(eVectSendC, len, MPI_INT, dest, tag+5, comm);
  free(eVectSendC);
}

void MPI_RECV_VectInt(vector<int> &eListI, int src, int tag, MPI_Comm comm)
{
  int *eVectRecvC;
  MPI_Status status;
  int len, i;
  if ((eVectRecvC = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  MPI_Recv(eVectRecvC, 1, MPI_INT, src, tag, comm, &status);
  len=eVectRecvC[0];
  free(eVectRecvC);
  
  if ((eVectRecvC = (int*)malloc(len*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  MPI_Recv(eVectRecvC, len, MPI_INT, src, tag+5, comm, &status);
  eListI.clear();
  for (i=0; i<len; i++)
    eListI.push_back(eVectRecvC[i]);
  free(eVectRecvC);
}

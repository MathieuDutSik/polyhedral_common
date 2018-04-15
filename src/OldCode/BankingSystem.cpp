#include "BankingSystem.h"


int GetIntFromVectChar(vector<char> TheInput)
{
  char *eStr;
  int siz, eVal, ret, i;
  siz=TheInput.size();
  if ((eStr = (char*)malloc(siz*sizeof(char))) == 0)
    exit (EXIT_FAILURE);
  for (i=0; i<siz; i++)
    eStr[i]=TheInput[i];
  ret=sscanf(eStr, "%d", &eVal);
  free(eStr);
  return eVal;
}


void GAP_ReadListInc(FILE *f, vector<vector<int> > &ListInc)
{
  char eChar;
  vector<char> ListChar;
  vector<int> eInc;
  int TheLevel, Act;
  TheLevel=0;
  while(1)
    {
      eChar=(char)fgetc(f);
      if (strcmp(&eChar,";") == 0)
	break;
      Act=0;
      if (strcmp(&eChar,"[") == 0)
	{
	  Act=1;
	  TheLevel++;
	}
      if (strcmp(&eChar,"]") == 0)
	{
	  Act=1;
	  if (TheLevel == 1)
	    {
	      ListInc.push_back(eInc);
	      eInc.clear();
	    }
	  if (TheLevel == 2)
	    {
	      eInc.push_back(GetIntFromVectChar(ListChar));
	      ListChar.clear();
	    }
	  TheLevel--;
	}
      if (TheLevel < 0)
	{
	  fprintf(stderr, "TheLevel=%d\n", TheLevel);
	  fprintf(stderr, "While it should be non-negative always\n");
	  exit(1);
	}
      if (strcmp(&eChar,"r") == 0 || strcmp(&eChar,"e") == 0 || strcmp(&eChar,"t") == 0 || strcmp(&eChar,"u") == 0 || strcmp(&eChar,"n") == 0)
	Act=1;
      fprintf(stderr, "TheLevel=%d\n", TheLevel);
      if (strcmp(&eChar, " ") != 0 && Act == 0)
	{
	  if (strcmp(&eChar,",") == 0)
	    {
	      if (TheLevel == 2)
		{
		  eInc.push_back(GetIntFromVectChar(ListChar));
		  ListChar.clear();
		}
	      if (TheLevel == 1)
		{
		  ListInc.push_back(eInc);
		  eInc.clear();
		}
	    }
	  else
	    ListChar.push_back(eChar);
	}
    }
  if (TheLevel != 0)
    {
      fprintf(stderr, "inconsistency in the code or in the file\n");
      exit(1);
    }
}


mpq_class GetMpqFromVectChar(vector<char> TheInput)
{
  char *eStr;
  mpq_t eVal_q;
  int siz, ret, i;
  mpq_class tcla;
  mpq_init(eVal_q);
  siz=TheInput.size();
  if ((eStr = (char*)malloc(siz*sizeof(char))) == 0)
    exit (EXIT_FAILURE);
  for (i=0; i<siz; i++)
    eStr[i]=TheInput[i];
  ret=mpq_set_str(eVal_q, eStr, 10);
  if (ret != 0)
    {
      fprintf(stderr, "Error in reading of mpq value\n");
      exit(1);
    }
  free(eStr);
  mpq_clear(eVal_q);
  return tcla;
}


void GAP_ReadEXT(FILE *f, vector<vector<mpq_class> > &EXT)
{
  char eChar;
  vector<char> ListChar;
  vector<mpq_class> eEXT;
  int TheLevel, Act;
  TheLevel=0;
  while(1)
    {
      eChar=(char)fgetc(f);
      if (strcmp(&eChar,";") == 0)
	break;
      Act=0;
      if (strcmp(&eChar,"[") == 0)
	{
	  Act=1;
	  TheLevel++;
	}
      if (strcmp(&eChar,"]") == 0)
	{
	  Act=1;
	  if (TheLevel == 1)
	    {
	      EXT.push_back(eEXT);
	      eEXT.clear();
	    }
	  if (TheLevel == 2)
	    {
	      eEXT.push_back(GetMpqFromVectChar(ListChar));
	      ListChar.clear();
	    }
	  TheLevel--;
	}
      if (TheLevel < 0)
	{
	  fprintf(stderr, "TheLevel=%d\n", TheLevel);
	  fprintf(stderr, "While it should be non-negative always\n");
	  exit(1);
	}
      if (strcmp(&eChar,"r") == 0 || strcmp(&eChar,"e") == 0 || strcmp(&eChar,"t") == 0 || strcmp(&eChar,"u") == 0 || strcmp(&eChar,"n") == 0)
	Act=1;
      fprintf(stderr, "TheLevel=%d\n", TheLevel);
      if (strcmp(&eChar, " ") != 0 && Act == 0)
	{
	  if (strcmp(&eChar,",") == 0)
	    {
	      if (TheLevel == 2)
		{
		  eEXT.push_back(GetMpqFromVectChar(ListChar));
		  ListChar.clear();
		}
	      if (TheLevel == 1)
		{
		  EXT.push_back(eEXT);
		  eEXT.clear();
		}
	    }
	  else
	    ListChar.push_back(eChar);
	}
    }
  if (TheLevel != 0)
    {
      fprintf(stderr, "inconsistency in the code or in the file\n");
      exit(1);
    }
}


void BANK_IsPresent(BankingDescription &TheBank,
		    MyMatrix *EXTask,
		    int *TheReply, vector<int> & TheEquiv)
{
  WeightMatrix WMatAsk;
  WeightMatrix WMatWork;
  int TheReplyInternal;
  int iEntry, nbEntry;
  GetWeightMatrix(WMatAsk, EXTask);
  nbEntry=TheBank.ListEXT.size();
  for (iEntry=0; iEntry<nbEntry; iEntry++)
    {
      GetWeightMatrix(WMatWork, TheBank.ListEXT[iEntry]);
      TestEquivalenceWeightMatrix(WMatAsk, WMatWork, &TheReplyInternal, TheEquiv);
      if (TheReplyInternal == 1)
	{
	  *TheReply=iEntry;
	  return;
	}
      free(WMatWork.TheMat);
    }
  free(WMatAsk.TheMat);
  *TheReply=-1;
}


void BANK_Clearing(BankingDescription &TheBank)
{
  int nbEntry, iEntry;
  nbEntry=TheBank.ListDualDesc.size();
  for (iEntry=0; iEntry<nbEntry; iEntry++)
    {
      MATRIX_Free(TheBank.ListEXT[iEntry]);
      free(TheBank.ListEXT[iEntry]);
    }
}


void BANK_ProcessRequest(BankingDescription &TheBank,
			 MyMatrix *EXTask,
			 int *TheReply, vector<vector<int> > &ListRepr)
{
  vector<int> TheEquiv;
  vector<int> eListI;
  int siz, i, iOrbit, nbOrbit;
  int iVert, jVert;
  BANK_IsPresent(TheBank, EXTask, TheReply, TheEquiv);
  if (*TheReply != -1)
    {
      ListRepr.clear();
      nbOrbit=TheBank.ListDualDesc[*TheReply].size();
      for (iOrbit=0; iOrbit<nbOrbit; iOrbit++)
	{
	  eListI.clear();
	  siz=TheBank.ListDualDesc[*TheReply][iOrbit].size();
	  for (i=0; i<siz; i++)
	    {
	      iVert=TheBank.ListDualDesc[*TheReply][iOrbit][i];
	      jVert=TheEquiv[iVert];
	      eListI.push_back(jVert);
	    }
	  ListRepr.push_back(eListI);
	}
    }
}


void BANK_InsertEntry(BankingDescription &TheBank, 
		      MyMatrix *EXTnew, 
		      vector<vector<int> > &ListReprNew)
{
  MyMatrix *eMat;
  int nbRow, nbCol;
  if ((eMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  nbRow=EXTnew->nbRow;
  nbCol=EXTnew->nbCol;
  MATRIX_Allocate(eMat, nbRow, nbCol);
  MATRIX_Copy(EXTnew, eMat);
  TheBank.ListEXT.push_back(eMat);
  TheBank.ListDualDesc.push_back(ListReprNew);
}


void BANK_OperatingPolyhedralBank(BankingDescription &TheBank, 
				  MPI_Comm comm)
{
  MPI_Status status;
  int *eVectOrig;
  int *eVectAction;
  int *eVectReply;
  MyMatrix *TheEXTask;
  vector<int> TheEquiv;
  vector<vector<int> > ListRepr;
  int TheSender, TheAct, TheReply, ret;
  if ((eVectOrig = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eVectAction = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eVectReply = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((TheEXTask = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  while(1)
    {
      MPI_Recv(eVectOrig, 1, MPI_INT, MPI_ANY_SOURCE, 1023, comm, &status);
      TheSender=eVectOrig[0];

      MPI_Recv(eVectAction, 1, MPI_INT, TheSender, 1492, comm, &status);
      TheAct=eVectAction[0];
      if (TheAct == 1)
	{
	  BANK_Clearing(TheBank);
	  break;
	}
      if (TheAct == 2)
	{
	  MPI_RECV_MyMatrix(TheEXTask, TheSender, 242, comm);
	  BANK_ProcessRequest(TheBank,
			      TheEXTask,
			      &TheReply, ListRepr);
	  eVectReply[0]=TheReply;
	  ret=MPI_Send(eVectReply, 1, MPI_INT, TheSender, 1897, comm);
	  if (TheReply == 1)
	    MPI_SEND_VectVectInt(ListRepr, TheSender, 1917, comm);
	  MATRIX_Free(TheEXTask);
	}
      if (TheAct == 3)
	{
	  MPI_RECV_MyMatrix(TheEXTask, TheSender, 702, comm);
	  BANK_IsPresent(TheBank,
			 TheEXTask,
			 &TheReply, TheEquiv);
	  eVectReply[0]=0;
	  ret=MPI_Send(eVectReply, 1, MPI_INT, TheSender, 2036, comm);
	  if (TheReply == 0)
	    {
	      MPI_RECV_VectVectInt(ListRepr, TheSender, 745, comm);
	      BANK_InsertEntry(TheBank, 
			       TheEXTask, 
			       ListRepr);
	    }
	  MATRIX_Free(TheEXTask);
	}
    }
  free(TheEXTask);
  free(eVectOrig);
  free(eVectAction);
  free(eVectReply);
}

#ifndef TEMP_BANKING_SYSTEM
#define TEMP_BANKING_SYSTEM


int GetIntFromVectChar(vector<char> const& TheInput)
{
  char *eStr;
  int siz, eVal, ret, i;
  siz=TheInput.size();
  if ((eStr = (char*)malloc(siz*sizeof(char))) == 0) {
    throw TerminalException{1};
  }
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
  while(1) {
    eChar=(char)fgetc(f);
    if (strcmp(&eChar,";") == 0)
      break;
    Act=0;
    if (strcmp(&eChar,"[") == 0) {
      Act=1;
      TheLevel++;
    }
    if (strcmp(&eChar,"]") == 0) {
      Act=1;
      if (TheLevel == 1) {
	ListInc.push_back(eInc);
	eInc.clear();
      }
      if (TheLevel == 2) {
	eInc.push_back(GetIntFromVectChar(ListChar));
	ListChar.clear();
      }
      TheLevel--;
    }
    if (TheLevel < 0) {
      std::cerr << "TheLevel=" << TheLevel << "\n";
      std::Cerr << "While it should be non-negative always\n";
      throw TerminalException{1};
    }
    if (strcmp(&eChar,"r") == 0 || strcmp(&eChar,"e") == 0 || strcmp(&eChar,"t") == 0 || strcmp(&eChar,"u") == 0 || strcmp(&eChar,"n") == 0)
      Act=1;
    std::cerr << "TheLevel=" << TheLevel << "\n";
    if (strcmp(&eChar, " ") != 0 && Act == 0) {
      if (strcmp(&eChar,",") == 0) {
	if (TheLevel == 2) {
	  eInc.push_back(GetIntFromVectChar(ListChar));
	  ListChar.clear();
	}
	if (TheLevel == 1) {
	  ListInc.push_back(eInc);
	  eInc.clear();
	}
      }
      else
	ListChar.push_back(eChar);
    }
  }
  if (TheLevel != 0) {
    std::cerr << "inconsistency in the code or in the file\n";
    throw TerminalException{1};
  }
}


template<typename T>
void BANK_IsPresent(BankingDescription &TheBank,
		    MyMatrix<T> const& EXTask,
		    int *TheReply, vector<int> & TheEquiv)
{
  WeightMatrix WMatAsk;
  WeightMatrix WMatWork;
  int TheReplyInternal;
  int iEntry, nbEntry;
  T_GetWeightMatrix(WMatAsk, EXTask);
  int nbEntry=TheBank.ListEXT.size();
  for (int iEntry=0; iEntry<nbEntry; iEntry++) {
    T_GetWeightMatrix(WMatWork, TheBank.ListEXT[iEntry]);
    ResultEquivWeightMatrix eResEquiv=TestEquivalenceWeightMatrix(WMatAsk, WMatWork);
    TheEquiv=eResEquiv.TheEquiv;
    if (eResEquiv.TheReply == true) {
      *TheReply=iEntry;
      return;
    }
    free(WMatWork.TheMat);
  }
  free(WMatAsk.TheMat);
  *TheReply=-1;
}

template<typename T>
void BANK_Clearing(BankingDescription &TheBank)
{
  int nbEntry, iEntry;
  nbEntry=TheBank.ListDualDesc.size();
  for (iEntry=0; iEntry<nbEntry; iEntry++) {
    TMat_Free(TheBank.ListEXT[iEntry]);
    free(TheBank.ListEXT[iEntry]);
  }
}


template<typename T>
void BANK_ProcessRequest(BankingDescription &TheBank,
			 MyMatrix<T> const& EXTask,
			 int *TheReply, vector<vector<int> > &ListRepr)
{
  vector<int> TheEquiv;
  vector<int> eListI;
  int siz, i, iOrbit, nbOrbit;
  int iVert, jVert;
  BANK_IsPresent(TheBank, EXTask, TheReply, TheEquiv);
  if (*TheReply != -1) {
    ListRepr.clear();
    nbOrbit=TheBank.ListDualDesc[*TheReply].size();
    for (iOrbit=0; iOrbit<nbOrbit; iOrbit++) {
      eListI.clear();
      siz=TheBank.ListDualDesc[*TheReply][iOrbit].size();
      for (i=0; i<siz; i++) {
	iVert=TheBank.ListDualDesc[*TheReply][iOrbit][i];
	jVert=TheEquiv[iVert];
	eListI.push_back(jVert);
      }
      ListRepr.push_back(eListI);
    }
  }
}


template<typename T>
void BANK_InsertEntry(BankingDescription &TheBank, 
		      MyMatrix<T> const& EXTnew, 
		      vector<vector<int> > &ListReprNew)
{
  MyMatrix<T> *eMat;
  int nbRow, nbCol;
  if ((eMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0) {
    throw TerminalException{1};
  }
  nbRow=EXTnew.GetNbRow();
  nbCol=EXTnew.GetNbCol();
  eMat=EXTnew;
  TMat_Copy(EXTnew, eMat);
  TheBank.ListEXT.push_back(eMat);
  TheBank.ListDualDesc.push_back(ListReprNew);
}


template<typename T>
void BANK_OperatingPolyhedralBank(BankingDescription &TheBank, 
				  MPI_Comm comm)
{
  MPI_Status status;
  int *eVectOrig;
  int *eVectAction;
  int *eVectReply;
  MyMatrix<T> *TheEXTask;
  vector<int> TheEquiv;
  vector<vector<int> > ListRepr;
  int TheSender, TheAct, TheReply, ret;
  if ((eVectOrig = (int*)malloc(1*sizeof(int))) == 0) {
    throw TerminalException{1};
  }
  if ((eVectAction = (int*)malloc(1*sizeof(int))) == 0) {
    throw TerminalException{1};
  }
  if ((eVectReply = (int*)malloc(1*sizeof(int))) == 0) {
    throw TerminalException{1};
  }
  if ((TheEXTask = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0) {
    throw TerminalException{1};
  }
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

#endif


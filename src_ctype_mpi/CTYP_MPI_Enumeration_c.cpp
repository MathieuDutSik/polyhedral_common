#include "CtypeMPI_types.h"
#include "NumberTheory.h"
#include "Namelist.h"
#include <unordered_map>

#include <boost/mpi.hpp>
#include "hash_functions.h"
namespace mpi = boost::mpi;


//#define TIMINGS_HASH
//#define ERR_LOG

/*
  Possible parallel schemes:
  ---We have a single file in input that is read at the beginning
  and data is dispatched to all the processors.
     PLUS: Simplicity. Only one exchange from one process to another.
     Can switch from one number of processors to another easily
     MINUS: Need for a merging function. At start, time spent on dispactching.
  ---We have a fixed number of entries.
     PLUS: No need for dispatching.
     MINUS: More exchanges during operation. More unwieldly database.
     When changing to another number of processor, great work needed.

 */
FullNamelist NAMELIST_GetStandard_ENUMERATE_CTYPE_MPI()
{
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListIntValues1["n"]=6;
  ListIntValues1["MaxNumberFlyingMessage"]=100;
  ListIntValues1["MaxStoredUnsentMatrices"]=1000;
  ListIntValues1["MaxRunTimeSecond"]=-1;
  ListIntValues1["TimeForDeclaringItOver"]=120;
  ListBoolValues1["StopWhenFinished"]=false;
  ListStringValues1["ListMatrixInput"] = "ListMatrix";
  //  ListStringValues1["PrefixDataSave"]="Output_";
  SingleBlock BlockDATA;
  BlockDATA.ListIntValues = ListIntValues1;
  BlockDATA.ListStringValues = ListStringValues1;
  BlockDATA.ListBoolValues = ListBoolValues1;
  ListBlock["DATA"]=BlockDATA;
  // Merging all data
  return {ListBlock, "undefined"};
}





static int tag_new_form = 37;
static int tag_termination = 157;


int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int irank_i, n_pes_i;
  MPI_Comm_size(MPI_COMM_WORLD, &n_pes_i);
  MPI_Comm_rank(MPI_COMM_WORLD,&irank_i);
  size_t irank=irank_i;
  size_t n_pes=n_pes_i;
#ifdef ERR_LOG
  std::cerr << "irank=" << irank << " n_pes=" << n_pes << "\n";
#endif
  //
  using Tint=int;
  //
  // The input file
  //
  FullNamelist eFull = NAMELIST_GetStandard_ENUMERATE_CTYPE_MPI();
  if (argc != 2) {
    std::cerr << "Number of argument is = " << argc << "\n";
    std::cerr << "This program is used as\n";
    std::cerr << "CTYP_MPI_Enumeration_c [file.nml]\n";
    std::cerr << "With file.nml a namelist file\n";
    NAMELIST_WriteNamelistFile(std::cerr, eFull);
    return -1;
  }
  std::string eFileName=argv[1];
  NAMELIST_ReadNamelistFile(eFileName, eFull);
  //
  // Parsing the input file
  //
  SingleBlock BlDATA = eFull.ListBlock["DATA"];
  //  int n=BlDATA.ListIntValues.at("n");
  int n = BlDATA.ListIntValues.at("n");
  int MaxNumberFlyingMessage = BlDATA.ListIntValues.at("MaxNumberFlyingMessage");
  int MaxStoredUnsentMatrices = BlDATA.ListIntValues.at("MaxStoredUnsentMatrices");
  int MaxRunTimeSecond = BlDATA.ListIntValues.at("MaxRunTimeSecond");
  int TimeForDeclaringItOver = BlDATA.ListIntValues.at("TimeForDeclaringItOver");
  bool StopWhenFinished = BlDATA.ListBoolValues.at("StopWhenFinished");
  std::string FileMatrix = BlDATA.ListStringValues.at("ListMatrixInput");
  //
  int n_vect = std::pow(2, n) - 1;
  int totalsiz_exch = 2*sizeof(int) + n_vect * n * sizeof(Tint) + sizeof(size_t) + 2*sizeof(int);
  std::string eFileO="LOG_" + std::to_string(irank);
  std::ofstream log(eFileO);
  log << "Initial log entry" << std::endl;
  //
  struct KeyData {
    int idxMatrix;
  };
  uint32_t seed= 0x1b873540;
  //
  // The list of requests.
  //
  std::vector<MPI_Request> ListRequest(MaxNumberFlyingMessage);
  std::vector<int> RequestStatus(MaxNumberFlyingMessage, 0);
  std::vector<std::vector<char>> ListMesg(MaxNumberFlyingMessage, std::vector<char>(totalsiz_exch));
  int nbRequest = 0;
  auto GetFreeIndex=[&]() -> int {
#ifdef ERR_LOG
    std::cerr << "Beginning of GetFreeIndex\n";
#endif
    for (int u=0; u<MaxNumberFlyingMessage; u++) {
#ifdef ERR_LOG
      std::cerr << "GetFreeIndex u=" << u << "\n";
#endif
      if (RequestStatus[u] == 0) {
#ifdef ERR_LOG
        std::cerr << "GetFreeIndex, returning u=" << u << "\n";
#endif
	return u;
      }
#ifdef ERR_LOG
      std::cerr << "Testing and getting a request\n";
#endif
      int flag;
      MPI_Status status1;
      int ierr1 = MPI_Test(&ListRequest[u], &flag, &status1);
      if (ierr1 != MPI_SUCCESS) {
        std::cerr << "Failing at MPI_Test\n";
        throw TerminalException{1};
      }
#ifdef ERR_LOG
      std::cerr << "flag=" << flag << "\n";
#endif
      if (flag) { // that request has ended. Let's read it.
        // As it turns out the MPI_Test does not set up the status1.MPI_ERROR
        // Thus the test should not be checked or it would led us to more strange error
        // that actually do not occur.
        // See https://www.mpich.org/static/docs/v3.2/www3/MPI_Test.html
        /*
	if (status1.MPI_ERROR != MPI_SUCCESS) {
	  std::cerr << "status1.tag = " << status1.MPI_TAG << "\n";
	  std::cerr << "status1.source = " << status1.MPI_SOURCE << " irank=" << irank << "\n";
	  std::cerr << "status1.error = " << status1.MPI_ERROR << "\n";
          char error_string[10000];
          int length_of_error_string;
          MPI_Error_string(status1.MPI_ERROR, error_string, &length_of_error_string);
          fprintf(stderr, "err: %s\n", error_string);
	  throw TerminalException{1};
	}
        */
	RequestStatus[u] = 0;
        nbRequest--;
#ifdef ERR_LOG
        std::cerr << "GetFreeIndex, clearing u=" << u << " returning it\n";
#endif
	return u;
      }
    }
    return -1;
  };
  //
  // The list of matrices being treated
  //
  std::unordered_map<TypeCtypeExch<Tint>,KeyData> ListCasesNotDone;
  std::unordered_map<TypeCtypeExch<Tint>,KeyData> ListCasesDone;
  int idxMatrixCurrent=0;
  auto fInsert=[&](PairExch<Tint> const& ePair) -> void {
    TypeCtypeExch<Tint> eCtype = ePair.eCtype;
#ifdef TIMINGS_HASH
    std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
    auto it1 = ListCasesDone.find(eCtype);
#ifdef TIMINGS_HASH
    std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
    std::cerr << "|HashMap1|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
    if (it1 != ListCasesDone.end()) {
      log << "Processed entry=" << ePair.eIndex << "END" << std::endl;
      return;
    }
#ifdef TIMINGS_HASH
    std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
#endif
    KeyData& eData = ListCasesNotDone[eCtype];
#ifdef TIMINGS_HASH
    std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
    std::cerr << "|HashMap2|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
#endif
    if (eData.idxMatrix != 0) {
      log << "Processed entry=" << ePair.eIndex << "END" << std::endl;
      return;
    }
    eData.idxMatrix = idxMatrixCurrent + 1;
    log << "Inserting New ctype" << ePair.eCtype << " idxMatrixCurrent=" << idxMatrixCurrent << " Obtained from " << ePair.eIndex << "END" << std::endl;
#ifdef ERR_LOG
    std::cerr << "Inserting new form, now we have |ListCasesNotDone|=" << ListCasesNotDone.size() << " |ListCasesDone|=" << ListCasesDone.size() << "\n";
#endif
    //    std::cerr << "idxMatrixCurrent=" << idxMatrixCurrent << " eCtype = " << ePair.eCtype << "\n";
    idxMatrixCurrent++;
  };
  auto GetUndoneEntry=[&]() -> boost::optional<std::pair<TypeCtypeExch<Tint>,int>> {
    auto it1 = ListCasesNotDone.begin();
    if (it1 != ListCasesNotDone.end()) {
      std::pair<TypeCtypeExch<Tint>,int> ePair = {it1->first, it1->second.idxMatrix-1};
      return boost::optional<std::pair<TypeCtypeExch<Tint>,int>>(ePair);
    }
    return {};
  };
  auto SetMatrixAsDone=[&](TypeCtypeExch<Tint> const& TheMat) -> void {
    KeyData eKey = ListCasesNotDone.at(TheMat);
    ListCasesNotDone.erase(TheMat);
    ListCasesDone[TheMat] = eKey;
  };
  //
  // The system for sending matrices
  //
  std::vector<std::pair<PairExch<Tint>, size_t>> ListMatrixUnsent;
  auto ClearUnsentAsPossible=[&]() -> void {
    int pos=ListMatrixUnsent.size() - 1;
    while(true) {
      if (pos == -1)
	break;
      int idx = GetFreeIndex();
      if (idx == -1)
	break;
#ifdef ERR_LOG
      std::cerr << "Assigning the request idx=" << idx << "\n";
      std::cerr << "world.isend to target =" << ListMatrixUnsent[pos].second << "\n";
#endif
      size_t iProc = ListMatrixUnsent[pos].second;
      PairExch_to_vectorchar(ListMatrixUnsent[pos].first, ListMesg[idx]);
      char* ptr = ListMesg[idx].data();
      MPI_Request* ereq_ptr = &ListRequest[idx];
      int ierr1 = MPI_Isend(ptr, totalsiz_exch, MPI_SIGNED_CHAR, iProc, tag_new_form, MPI_COMM_WORLD, ereq_ptr);
      if (ierr1 != MPI_SUCCESS) {
        std::cerr << "ierr1 wrongly set\n";
        throw TerminalException{1};
      }
      RequestStatus[idx] = 1;
      nbRequest++;
#ifdef ERR_LOG
      std::cerr << "Sending matrix Ctype=" << ListMatrixUnsent[pos].first.eCtype << " index=" << ListMatrixUnsent[pos].first.eIndex << "\n";
#endif
      ListMatrixUnsent.pop_back();
      pos--;
    }
  };
  auto fInsertUnsent=[&](PairExch<Tint> const& ePair) -> void {
    size_t e_hash = Matrix_Hash(ePair.eCtype.eMat, seed);
    size_t res = e_hash % n_pes;
    //    std::cerr << "fInsertUnsent e_hash=" << e_hash << " res=" << res << "\n";
    if (res == irank) {
      fInsert(ePair);
    }
    else {
      ListMatrixUnsent.push_back({ePair, res});
      ClearUnsentAsPossible();
    }
  };
  //
  // Reading the initial file
  //
  {
    std::ifstream is(FileMatrix);
#ifdef ERR_LOG
    std::cerr << "Beginning reading file=" << FileMatrix << "\n";
#endif
    int nbMatrixStart;
    is >> nbMatrixStart;
    for (int iMatStart=0; iMatStart<nbMatrixStart; iMatStart++) {
      int eStatus;
      is >> eStatus;
      MyMatrix<Tint> TheMat = ReadMatrix<Tint>(is);
      TypeCtypeExch<Tint> eRecMat{TheMat};
      size_t e_hash = Matrix_Hash(TheMat, seed);
      size_t res = e_hash % n_pes;
#ifdef ERR_LOG
      std::cerr << "iMatStart=" << iMatStart << " e_hash=" << e_hash << " res=" << res << "\n";
#endif
      if (res == irank) {
        KeyData eData{idxMatrixCurrent+1};
        if (eStatus == 0) {
          ListCasesNotDone[eRecMat] = eData;
        }
        else {
          ListCasesDone[eRecMat] = eData;
        }
        log << "Reading existing matrix=" << eRecMat << " idxMatrixCurrent=" << idxMatrixCurrent << "END" << std::endl;
        idxMatrixCurrent++;
      }
    }
  }
#ifdef ERR_LOG
  std::cerr << "Reading finished, we have |ListCasesDone|=" << ListCasesDone.size() << " |ListCasesNotDone|=" << ListCasesNotDone.size() << "\n";
#endif
  //
  // The main loop itself.
  //
  int iVal_synchronization = 72;
  bool TerminationNoticeSent = false;
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
  std::chrono::time_point<std::chrono::system_clock> last_timeoper = start;
  std::vector<int> StatusNeighbors(n_pes, 0);
  while(true) {
    // The reference time used for the comparisons
    std::chrono::time_point<std::chrono::system_clock> ref_time = std::chrono::system_clock::now();
    int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(ref_time - start).count();
    // Now the operations themselves
#ifdef ERR_LOG
    std::cerr << "Begin while, we have |ListCasesNotDone|=" << ListCasesNotDone.size() << " |ListCasesDone|=" << ListCasesDone.size() << " elapsed_time=" << elapsed_seconds << "\n";
#endif
    MPI_Status status1;
    int flag;
    int ierr1 = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag,&status1);
    if (ierr1 != MPI_SUCCESS) {
      std::cerr << "ierr1 wrongly set\n";
      throw TerminalException{1};
    }
    if (flag) {
      if (status1.MPI_ERROR != MPI_SUCCESS) {
        std::cerr << "Having an MPI_Error. Immediate death\n";
        throw TerminalException{1};
      }
      if (status1.MPI_TAG == tag_new_form) {
        StatusNeighbors[status1.MPI_SOURCE] = 0; // Getting a message pretty much means it is alive
        std::vector<char> eVect_c(totalsiz_exch);
        char* ptr = eVect_c.data();
        MPI_Status status2;
        int ierr2 = MPI_Recv(ptr, totalsiz_exch, MPI_SIGNED_CHAR, status1.MPI_SOURCE, tag_new_form,
                             MPI_COMM_WORLD, &status2);
        if (status2.MPI_ERROR != MPI_SUCCESS || ierr2 != MPI_SUCCESS) {
          std::cerr << "Failed status2 or ierr2\n";
          throw TerminalException{1};
        }
	PairExch<Tint> ePair = vectorchar_to_PairExch<Tint>(eVect_c);
#ifdef ERR_LOG
	std::cerr << "Receiving a matrix ePair=" << ePair.eCtype << " index=" << ePair.eIndex << "\n";
#endif
        fInsert(ePair);
        // Now the timings
        last_timeoper = std::chrono::system_clock::now();
      }
      if (status1.MPI_TAG == tag_termination) {
        StatusNeighbors[status1.MPI_SOURCE] = 1; // This is the termination message
        // Below is just customary. We are not really interested in the received value.
        int RecvInt;
        int* ptr_i = &RecvInt;
        MPI_Status status2;
        int ierr2 = MPI_Recv(ptr_i, 1, MPI_INT, status1.MPI_SOURCE, tag_termination,
                             MPI_COMM_WORLD, &status2);
        if (status2.MPI_ERROR != MPI_SUCCESS || ierr2 != MPI_SUCCESS) {
          std::cerr << "Failed status2 or ierr\n";
          throw TerminalException{1};
        }
        if (RecvInt != iVal_synchronization) {
          std::cerr << "The received integer is not what we expect\n";
          throw TerminalException{1};
        }
      }
    } else {
#ifdef ERR_LOG
      std::cerr << "irank=" << irank << " |ListMatrixUnsent|=" << ListMatrixUnsent.size() << " MaxStoredUnsentMatrices=" << MaxStoredUnsentMatrices << "\n";
#endif
      bool DoSomething = false;
      if (int(ListMatrixUnsent.size()) < MaxStoredUnsentMatrices) {
        if (MaxRunTimeSecond < 0 || elapsed_seconds < MaxRunTimeSecond) {
          // We pass the first test with respect to runtime
          if (!StopWhenFinished || ListCasesNotDone.size() > 0) {
            // We pass the test of being non-empty
            DoSomething = true;
          }
        }
      }
      // Now finding the adjacent if we indeed do something.
#ifdef ERR_LOG
      std::cerr << "DoSomething = " << DoSomething << "\n";
#endif
      if (DoSomething) {
	boost::optional<std::pair<TypeCtypeExch<Tint>,int>> eReq=GetUndoneEntry();
	if (eReq) {
          StatusNeighbors[irank] = 0;
          //          std::cerr << "eReq is non zero\n";
	  SetMatrixAsDone(eReq->first);
          //          std::cerr << "irank=" << irank << " eCtype=" << eReq->first << "\n";
          int idxMatrixF = eReq->second;
#ifdef ERR_LOG
          std::cerr << "Starting Adjacent Form Method\n";
#endif
          std::vector<TypeCtypeExch<Tint>> ListAdjacentObject = CTYP_GetAdjacentCanonicCtypes<Tint>(eReq->first);
          int nbAdjacent = ListAdjacentObject.size();
          log << "Number of Adjacent for idxMatrixF=" << idxMatrixF << " nbAdjacent=" << nbAdjacent << " END" << std::endl;
#ifdef ERR_LOG
          std::cerr << "Number of Adjacent for idxMatrixF=" << idxMatrixF << " nbAdjacent=" << nbAdjacent << " END\n";
#endif
          int iAdj=0;
	  for (auto & eObj1 : ListAdjacentObject) {
            TypeIndex eIndex{irank, idxMatrixF, iAdj};
            PairExch<Tint> ePair{eObj1, eIndex};
	    fInsertUnsent(ePair);
            iAdj++;
	  }
          // Now the timings
          last_timeoper = std::chrono::system_clock::now();
	}
      }
    }
    ClearUnsentAsPossible();
    //
    // Sending messages for the synchronization of ending the run
    //
    int TimeClear = std::chrono::duration_cast<std::chrono::seconds>(ref_time - last_timeoper).count();
#ifdef ERR_LOG
    std::cerr << "TimeClear=" << TimeClear << " TimeForDeclaringItOver=" << TimeForDeclaringItOver << "\n";
#endif
    if (TimeClear > TimeForDeclaringItOver && ListMatrixUnsent.size() == 0 && !TerminationNoticeSent) {
      TerminationNoticeSent = true;
      for (size_t i_pes=0; i_pes<n_pes; i_pes++) {
#ifdef ERR_LOG
        std::cerr << "Before world.isend for termination i_pes=" << i_pes << "\n";
#endif
        if (i_pes == irank) {
          StatusNeighbors[i_pes] = 1;
        } else {
          int idx = GetFreeIndex();
          if (idx == -1) {
            std::cerr << "We should be able to have an entry\n";
            throw TerminalException{1};
          }
          int* ptr_i = &iVal_synchronization;
          MPI_Request* ereq_ptr = &ListRequest[idx];
          int ierr1 = MPI_Isend(ptr_i, 1, MPI_INT, i_pes, tag_termination, MPI_COMM_WORLD, ereq_ptr);
          if (ierr1 != MPI_SUCCESS) {
            std::cerr << "ierr1 wrongly set\n";
            throw TerminalException{1};
          }
          RequestStatus[idx] = 1;
          nbRequest++;
        }
      }
    }
    //
    // We now clear with the messages
    //
    size_t nb_finished = 0;
    for (size_t i_pes=0; i_pes<n_pes; i_pes++) {
      nb_finished += StatusNeighbors[i_pes];
    }
    if (nb_finished == n_pes) {
#ifdef ERR_LOG
      std::cerr << "Before the all_reduce operation\n";
#endif
      int val_i=1;
      int val_o;
      int ierr1 = MPI_Allreduce(&val_i, &val_o, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if (ierr1 != MPI_SUCCESS) {
        std::cerr << "ierr1 wrongly set\n";
        throw TerminalException{1};
      }
      if (val_o == 1) {
        std::cerr << "Receive the termination message. All Exiting\n";
        std::cerr << "irank=" << irank << " |ListCasesDone|=" << ListCasesDone.size() << " |ListCasesNotDone|=" << ListCasesNotDone.size() << "\n";
        break;
      }
    }
  }
  std::cerr << "Normal termination of the program\n";
  MPI_Finalize();
}

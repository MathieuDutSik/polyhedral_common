#include "CtypeMPI_types.h"
#include "NumberTheory.h"
#include "Namelist.h"
#include <unordered_map>

#include <boost/mpi.hpp>
#include "hash_functions.h"
#include <netcdf>
namespace mpi = boost::mpi;


//#define TIMINGS_HASH
#define ERR_LOG

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
  It turns out that the runing of the system requires the second system.
  This is because the time for merging can be similar to the runtime which is
  unacceptable.

  ---

  The plus point of the system is that by exiting cleanly when the computation
  finifhsed, we do not need to keep track of just everything. The TypeIndex
  is actually not really needed anymore.

  ---

  The optimization of the hash table by using TSL. This could save some memory.
  The maximum of the coefficients has to be determined and move to int8_t

  ---

  With the estimate of N = 10^8 It gives us
  We have sizeof(size_t) = 8 in the platform used for this computation.
  or more precisely 4 because we use murmurhash32.
  So, this gets us 2^32 possibilities.
  Looking at the birthday problem https://en.wikipedia.org/wiki/Birthday_problem
  this gets us:
  n = 10^8 = 2^26 ...
  d = 2^32
  So the probability of no collision is going to be
  exp(-n^2 / 2d) = exp(-2^(2*26 - 1 - 32)) = exp(-2^19)
  So collision WILL definitely happen (if we use 32 bit hash).
  If we use a 64 bit hash then the probability is going to be zero.
  Therefore, we will need a 64 bit hash, whatever we decide to do.

  ---

  If we store the data on disk then we may save some memory.
  But the expense is going to be considerable. Thus this strategy
  is a last ressort one.

  ---

  The MPI exchanges need to be buffered:
  ---Exchanging the data with (nb, blk_1, ...., blk_nb, 0 ... 0)
  Some characters to exchange.
  ---We Need to have a variable like MpiBufferSize for the number
  of entries to exchange. Size 20 as a starting point.
  ---The working size to exchange has to vary in time. This is because
  we cannot exit the program without having all buffers being cleared.
  ---So, we have a variable CurrentBufferSize which satisfies
       1 <= CurrentBufferSize <= MpiBufferSize
  ---The list of matrices to exchange has to be partitionned by block.
  When the size of one becomes greater than CurrentBufferSize then
  we send data.
  ---If no operation has been done (exchange or adjacency) then we
  decrease CurrentBufferSize (if lower than 1).
  ---If some operation is done then it should increase because we
  aggressively want to have big blocks.

 */
FullNamelist NAMELIST_GetStandard_ENUMERATE_CTYPE_MPI()
{
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListIntValues1["n"]=6;
  ListStringValues1["WorkingPrefix"] = "LOG_";
  ListIntValues1["MaxNumberFlyingMessage"]=100;
  ListIntValues1["MaxStoredUnsentMatrices"]=1000;
  ListIntValues1["MaxRunTimeSecond"]=-1;
  ListIntValues1["TimeForDeclaringItOver"]=120;
  ListIntValues1["MpiBufferSize"]=20;
  ListBoolValues1["StopWhenFinished"]=false;
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


template<typename T, typename Tint>
void NC_ReadMatrix_T(netCDF::NcVar & varCtype, MyMatrix<int> & M, size_t const& n_vect, size_t const& n, int const& pos)
{
  std::vector<size_t> start2{size_t(pos), 0, 0};
  std::vector<size_t> count2{1, n_vect, n};
  std::vector<T> V(n_vect * n);
  varCtype.getVar(start2, count2, V.data());
  int idx=0;
  for (size_t i_vect=0; i_vect<n_vect; i_vect++)
    for (size_t i=0; i<n; i++) {
      M(i_vect, i) = V[idx];
      idx++;
    }
}


template<typename T, typename Tint>
void NC_WriteMatrix_T(netCDF::NcVar & varCtype, MyMatrix<int> const& M, size_t const& n_vect, size_t const& n, int const& pos)
{
  std::vector<size_t> start2{size_t(pos), 0, 0};
  std::vector<size_t> count2{1, n_vect, n};
  std::vector<T> V(n_vect * n);
  int idx=0;
  for (size_t i_vect=0; i_vect<n_vect; i_vect++)
    for (size_t i=0; i<n; i++) {
      V[idx] = M(i_vect, i);
      idx++;
    }
  varCtype.putVar(start2, count2, V.data());
}






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
  size_t MaxStoredUnsentMatrices = BlDATA.ListIntValues.at("MaxStoredUnsentMatrices");
  int MaxRunTimeSecond = BlDATA.ListIntValues.at("MaxRunTimeSecond");
  size_t MpiBufferSize = BlDATA.ListIntValues.at("MpiBufferSize");
  int TimeForDeclaringItOver = BlDATA.ListIntValues.at("TimeForDeclaringItOver");
  bool StopWhenFinished = BlDATA.ListBoolValues.at("StopWhenFinished");
  std::string WorkingPrefix = BlDATA.ListStringValues.at("WorkingPrefix");
  //
  // The basic sizes
  //
  int n_vect = std::pow(2, n) - 1;
  int siz_pairexch = n_vect * n * sizeof(Tint);
  int totalsiz_exch = sizeof(int) + MpiBufferSize * siz_pairexch;
  //
  // The netcdf interface
  //
  std::string WorkFile=WorkingPrefix + std::to_string(irank) + ".nc";
  netCDF::NcFile dataFile(WorkFile, netCDF::NcFile::write);
  netCDF::NcVar varCtype=dataFile.getVar("Ctype");
  int n_read = varCtype.getDim(2).getSize();
  if (n_read != n) {
    std::cerr << "n_read=" << n_read << " n=" << n << "\n";
    return 0;
  }
  netCDF::NcType eType=varCtype.getType();
  netCDF::NcVar varNbAdj=dataFile.getVar("nb_adjacent");
  int curr_nb_matrix = varNbAdj.getDim(0).getSize();
  auto NC_GetNbAdjacent=[&](int const& pos) -> int {
    std::vector<size_t> start{size_t(pos)};
    std::vector<size_t> count{1};
    int nbAdjacent;
    varNbAdj.getVar(start, count, &nbAdjacent);
    return nbAdjacent;
  };
  auto NC_WriteNbAdjacent=[&](int const& pos, int const& nbAdjacent) -> void {
    std::vector<size_t> start{size_t(pos)};
    std::vector<size_t> count{1};
    varNbAdj.putVar(start, count, &nbAdjacent);
  };
  auto NC_ReadMatrix=[&](int const& pos) -> TypeCtypeExch<Tint> {
    MyMatrix<Tint> M(n_vect, n);
    if (eType == netCDF::NcType::nc_BYTE)
      NC_ReadMatrix_T<int8_t,Tint>(varCtype, M, n_vect, n, pos);
    if (eType == netCDF::NcType::nc_SHORT)
      NC_ReadMatrix_T<int16_t,Tint>(varCtype, M, n_vect, n, pos);
    if (eType == netCDF::NcType::nc_INT)
      NC_ReadMatrix_T<int32_t,Tint>(varCtype, M, n_vect, n, pos);
    if (eType == netCDF::NcType::nc_INT64)
      NC_ReadMatrix_T<int64_t,Tint>(varCtype, M, n_vect, n, pos);
    return {M};
  };
  auto NC_AppendMatrix=[&](MyMatrix<Tint> const& M) -> void {
    if (eType == netCDF::NcType::nc_BYTE)
      NC_WriteMatrix_T<int8_t,Tint>(varCtype, M, n_vect, n, curr_nb_matrix);
    if (eType == netCDF::NcType::nc_SHORT)
      NC_WriteMatrix_T<int16_t,Tint>(varCtype, M, n_vect, n, curr_nb_matrix);
    if (eType == netCDF::NcType::nc_INT)
      NC_WriteMatrix_T<int32_t,Tint>(varCtype, M, n_vect, n, curr_nb_matrix);
    if (eType == netCDF::NcType::nc_INT64)
      NC_WriteMatrix_T<int64_t,Tint>(varCtype, M, n_vect, n, curr_nb_matrix);
    NC_WriteNbAdjacent(curr_nb_matrix, 0);
    curr_nb_matrix++;
  };
  //
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
  auto fctHash=[](size_t const& val) -> size_t {
    return val;
  };
  auto fctEqual=[](size_t const& val1, size_t const& val2) -> bool {
    return val1 == val2;
  };
  std::unordered_map<size_t,std::vector<int>,decltype(fctHash), decltype(fctEqual)> MapIndexByHash({}, fctHash, fctEqual);
  std::vector<int> ListUndoneIndex;
  int idxMatrixCurrent=0;
  auto fInsert=[&](TypeCtypeExch<Tint> const& eCtype) -> void {
#ifdef TIMINGS_HASH
    std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
    size_t e_hash = std::hash<TypeCtypeExch<Tint>>()(eCtype);
#ifdef ERR_LOG
    std::cerr << "e_hash=" << e_hash << "\n";
#endif
    std::vector<int> & eList = MapIndexByHash[e_hash];
#ifdef TIMINGS_HASH
    std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
    std::cerr << "|HashMap|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
    for (auto iIdx : eList) {
      TypeCtypeExch<Tint> fCtype = NC_ReadMatrix(iIdx);
      if (eCtype == fCtype)
        return;
    }
    eList.push_back(idxMatrixCurrent);
    ListUndoneIndex.push_back(idxMatrixCurrent);
    NC_AppendMatrix(eCtype.eMat);
    idxMatrixCurrent++;
#ifdef ERR_LOG
    int nb_collision=0;
    for (auto & kv : MapIndexByHash) {
      if (kv.second.size() > 1)
        nb_collision++;
    }
    std::cerr << "nb_collision=" << nb_collision << "\n";
#endif
  };
  auto GetUndoneEntry=[&]() -> boost::optional<std::pair<TypeCtypeExch<Tint>,int>> {
    size_t len = ListUndoneIndex.size();
    if (len > 0) {
      int idx = ListUndoneIndex[len - 1];
      ListUndoneIndex.pop_back();
      TypeCtypeExch<Tint> eCtype = NC_ReadMatrix(idx);
      std::pair<TypeCtypeExch<Tint>,int> ePair = {eCtype, idx};
      return boost::optional<std::pair<TypeCtypeExch<Tint>,int>>(ePair);
    }
    return {};
  };
  //
  // The system for sending matrices
  //
  std::vector<std::vector<TypeCtypeExch<Tint>>> ListListMatrixUnsent(n_pes);
  size_t CurrentBufferSize = 1;
  auto GetSendableIndex=[&]() -> int {
    for (size_t i_pes=0; i_pes<n_pes; i_pes++) {
      if (ListListMatrixUnsent[i_pes].size() >= CurrentBufferSize)
        return i_pes;
    }
    return -1;
  };
  auto AreBufferFullEnough=[&]() -> bool {
    size_t nb_unsend = 0;
#ifdef ERR_LOG
    std::cerr << "List|ListMatrixUnsent| =";
#endif
    for (size_t i_pes=0; i_pes<n_pes; i_pes++) {
      size_t the_siz = ListListMatrixUnsent[i_pes].size();
#ifdef ERR_LOG
      std::cerr << " " << the_siz;
#endif
      nb_unsend += the_siz;
    }
#ifdef ERR_LOG
    std::cerr << " nb_unsend=" << nb_unsend << "\n";
#endif
    return nb_unsend > MaxStoredUnsentMatrices;
  };
  auto ClearUnsentAsPossible=[&]() -> void {
    while(true) {
      int i_pes=GetSendableIndex();
      if (i_pes == -1)
	break;
      int idx = GetFreeIndex();
      if (idx == -1)
	break;
#ifdef ERR_LOG
      std::cerr << "Assigning the request idx=" << idx << " to processor " << i_pes << "\n";
#endif
      std::vector<TypeCtypeExch<Tint>>& eList = ListListMatrixUnsent[i_pes];
      char* ptr_o = ListMesg[idx].data();
      std::memcpy(ptr_o, (char*)(&CurrentBufferSize), sizeof(int));
      ptr_o += sizeof(int);
      int pos = eList.size();
      for (size_t i_mat=0; i_mat<CurrentBufferSize; i_mat++) {
        pos--;
        PairExch_to_vectorchar(eList[pos], n_vect, n, ptr_o);
        ptr_o += siz_pairexch;
#ifdef ERR_LOG
        std::cerr << "Appending matrix Ctype=" << eList[pos] << "\n";
#endif
        eList.pop_back();
      }
      char* ptr_send = ListMesg[idx].data();
      MPI_Request* ereq_ptr = &ListRequest[idx];
      int ierr1 = MPI_Isend(ptr_send, totalsiz_exch, MPI_SIGNED_CHAR, i_pes, tag_new_form, MPI_COMM_WORLD, ereq_ptr);
      if (ierr1 != MPI_SUCCESS) {
        std::cerr << "ierr1 wrongly set\n";
        throw TerminalException{1};
      }
      RequestStatus[idx] = 1;
      nbRequest++;
    }
  };
  auto fInsertUnsent=[&](TypeCtypeExch<Tint> const& eCtype) -> void {
    size_t e_hash = Matrix_Hash(eCtype.eMat, seed);
    size_t res = e_hash % n_pes;
    //    std::cerr << "fInsertUnsent e_hash=" << e_hash << " res=" << res << "\n";
    if (res == irank) {
      fInsert(eCtype);
    }
    else {
      ListListMatrixUnsent[res].push_back(eCtype);
      ClearUnsentAsPossible();
    }
  };
  //
  // Reading the initial file in memory
  //
#ifdef ERR_LOG
  std::cerr << "Beginning reading WorkFile=" << WorkFile << "\n";
#endif
  for (int iCurr=0; iCurr<curr_nb_matrix; iCurr++) {
    TypeCtypeExch<Tint> eCtype = NC_ReadMatrix(iCurr);
    size_t e_hash = std::hash<TypeCtypeExch<Tint>>()(eCtype);
    int nbAdjacent = NC_GetNbAdjacent(iCurr);
#ifdef ERR_LOG
    std::cerr << "iCurr=" << iCurr << "\n";
#endif
    std::vector<int>& eList= MapIndexByHash[e_hash];
    eList.push_back(idxMatrixCurrent);
    if (nbAdjacent == 0)
      ListUndoneIndex.push_back(idxMatrixCurrent);
    idxMatrixCurrent++;
  }
#ifdef ERR_LOG
  std::cerr << "Reading finished, we have |ListCases|=" << idxMatrixCurrent << " |ListCasesNotDone|=" << ListUndoneIndex.size() << "\n";
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
    bool DidSomething = false;
    std::chrono::time_point<std::chrono::system_clock> ref_time = std::chrono::system_clock::now();
    int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(ref_time - start).count();
    // Now the operations themselves
#ifdef ERR_LOG
    std::cerr << "Begin while, we have |ListCases|=" << idxMatrixCurrent << " |ListCasesNotDone|=" << ListUndoneIndex.size() << " elapsed_time=" << elapsed_seconds << "\n";
#endif
    MPI_Status status1;
    int flag;
    int ierr1 = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag,&status1);
    if (ierr1 != MPI_SUCCESS) {
      std::cerr << "ierr1 wrongly set\n";
      throw TerminalException{1};
    }
    if (flag) {
      DidSomething=true;
      if (status1.MPI_ERROR != MPI_SUCCESS) {
        std::cerr << "Having an MPI_Error. Immediate death\n";
        throw TerminalException{1};
      }
      if (status1.MPI_TAG == tag_new_form) {
        StatusNeighbors[status1.MPI_SOURCE] = 0; // Getting a message pretty much means it is alive
        std::vector<char> eVect_c(totalsiz_exch);
        char* ptr_recv = eVect_c.data();
        MPI_Status status2;
        int ierr2 = MPI_Recv(ptr_recv, totalsiz_exch, MPI_SIGNED_CHAR, status1.MPI_SOURCE, tag_new_form,
                             MPI_COMM_WORLD, &status2);
        if (status2.MPI_ERROR != MPI_SUCCESS || ierr2 != MPI_SUCCESS) {
          std::cerr << "Failed status2 or ierr2\n";
          throw TerminalException{1};
        }
        int nbRecv;
        std::memcpy((char*)(&nbRecv), ptr_recv, sizeof(int));
        ptr_recv += sizeof(int);
#ifdef ERR_LOG
	std::cerr << "Receiving nbRecv=" << nbRecv << " matrices\n";
#endif
        for (int iRecv=0; iRecv<nbRecv; iRecv++) {
          TypeCtypeExch<Tint> eCtype = vectorchar_to_PairExch<Tint>(ptr_recv, n_vect, n);
#ifdef ERR_LOG
          std::cerr << "iRecv=" << iRecv << " ctype=" << eCtype << "\n";
#endif
          ptr_recv += siz_pairexch;
          fInsert(eCtype);
        }
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
      std::cerr << "irank=" << irank <<  " MaxStoredUnsentMatrices=" << MaxStoredUnsentMatrices << "\n";
#endif
      bool DoSomething = false;
      if (!AreBufferFullEnough()) {
        if (MaxRunTimeSecond < 0 || elapsed_seconds < MaxRunTimeSecond) {
          // We pass the first test with respect to runtime
          if (!StopWhenFinished || ListUndoneIndex.size() > 0) {
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
          DidSomething = true;
          StatusNeighbors[irank] = 0;
          int idxMatrixF = eReq->second;
#ifdef ERR_LOG
          std::cerr << "Starting Adjacent Form Method\n";
          std::cerr << "eReq->first=" << eReq->first << "\n";
          //          WriteMatrix(std::cerr, eReq->first.eMat);
#endif
          std::vector<TypeCtypeExch<Tint>> ListAdjacentObject = CTYP_GetAdjacentCanonicCtypes<Tint>(eReq->first);
#ifdef ERR_LOG
          std::cerr << "We have ListAdjacentObject\n";
#endif
          int nbAdjacent = ListAdjacentObject.size();
          NC_WriteNbAdjacent(idxMatrixF, nbAdjacent);
#ifdef ERR_LOG
          std::cerr << "Number of Adjacent for idxMatrixF=" << idxMatrixF << " nbAdjacent=" << nbAdjacent << " END\n";
#endif
	  for (auto & eObj1 : ListAdjacentObject)
	    fInsertUnsent(eObj1);
          // Now the timings
          last_timeoper = std::chrono::system_clock::now();
	}
      }
    }
    ClearUnsentAsPossible();
    //
    if (DidSomething) {
      if (CurrentBufferSize < MpiBufferSize)
        CurrentBufferSize++;
    } else {
      if (CurrentBufferSize > 1)
        CurrentBufferSize--;
    }
    //
    // Sending messages for the synchronization of ending the run
    //
    int TimeClear = std::chrono::duration_cast<std::chrono::seconds>(ref_time - last_timeoper).count();
    size_t nb_unsent = 0;
    for (size_t i_pes=0; i_pes<n_pes; i_pes++)
      nb_unsent += ListListMatrixUnsent[i_pes].size();
#ifdef ERR_LOG
    std::cerr << "TimeClear=" << TimeClear << " TimeForDeclaringItOver=" << TimeForDeclaringItOver << "\n";
#endif
    if (TimeClear > TimeForDeclaringItOver && nb_unsent == 0 && CurrentBufferSize == 1 && !TerminationNoticeSent) {
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
      int64_t nbAll = idxMatrixCurrent;
      int64_t nbNotDone = ListUndoneIndex.size();
      int64_t nbAll_tot, nbNotDone_tot;
      int ierr1 = MPI_Allreduce(&nbAll, &nbAll_tot, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
      if (ierr1 != MPI_SUCCESS) {
        std::cerr << "ierr1 wrongly set\n";
        throw TerminalException{1};
      }
      int ierr2 = MPI_Allreduce(&nbNotDone, &nbNotDone_tot, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
      if (ierr2 != MPI_SUCCESS) {
        std::cerr << "ierr2 wrongly set\n";
        throw TerminalException{1};
      }
      std::cerr << "Receive the termination message. All Exiting\n";
      std::cerr << "FINAL irank=" << irank << " local |ListCases|=" << nbAll << " |ListCasesNotDone|=" << nbNotDone << "\n";
      std::cerr << "FINAL irank=" << irank << " total |ListCases|=" << nbAll_tot << " |ListCasesNotDone|=" << nbNotDone_tot << "\n";
      break;
    }
  }
  std::cerr << "Normal termination of the program\n";
  MPI_Finalize();
}

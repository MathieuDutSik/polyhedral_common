#include "CtypeMPI_types.h"
#include "NumberTheory.h"
#include "Namelist.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PerfectForm.h"
#include <unordered_map>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include "md5sum.h"
namespace mpi = boost::mpi;


FullNamelist NAMELIST_GetStandard_ENUMERATE_CTYPE_MPI()
{
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListIntValues1["n"]=9;
  ListIntValues1["MaxNumberFlyingMessage"]=100;
  ListIntValues1["MaxStoredUnsentMatrices"]=1000;
  ListIntValues1["MaxRunTimeSecond"]=-1;
  ListStringValues1["ListMatrixInput"] = "ListMatrix";
  //  ListStringValues1["PrefixDataSave"]="Output_";
  SingleBlock BlockDATA;
  BlockDATA.ListIntValues=ListIntValues1;
  BlockDATA.ListBoolValues=ListBoolValues1;
  BlockDATA.ListDoubleValues=ListDoubleValues1;
  BlockDATA.ListStringValues=ListStringValues1;
  BlockDATA.ListListStringValues=ListListStringValues1;
  ListBlock["DATA"]=BlockDATA;
  // Merging all data
  return {ListBlock, "undefined"};
}


template<typename T, typename Tint>
std::vector<TypeCtypeExch<Tint>> GetAdjacentObjects(TypeCtypeExch<Tint> const& eObjIn)
{
}








static int tag_new_form = 37;
static int tag_written_form = 38;


int main()
{
  using T=mpq_class;
  using Tint=long;
  //
  FullNamelist eFull = NAMELIST_GetStandard_ENUMERATE_CTYPE_MPI();
  std::string eFileName = "perfectenum.nml";
  NAMELIST_ReadNamelistFile(eFileName, eFull);
  SingleBlock BlDATA = eFull.ListBlock["DATA"];
  //  int n=BlDATA.ListIntValues.at("n");
  int MaxNumberFlyingMessage = BlDATA.ListIntValues.at("MaxNumberFlyingMessage");
  int MaxStoredUnsentMatrices = BlDATA.ListIntValues.at("MaxStoredUnsentMatrices");
  int MaxRunTimeSecond = BlDATA.ListIntValues.at("MaxRunTimeSecond");
  std::string FileMatrix = BlDATA.ListStringValues.at("ListMatrixInput");
  //
  boost::mpi::environment env;
  boost::mpi::communicator world;
  int irank=world.rank();
  int n_pes=world.size();
  std::string eFileO="LOG_" + IntToString(irank);
  std::ofstream log(eFileO);
  log << "Initial log entry" << std::endl;
  //
  std::ifstream is(FileMatrix);
  int nbMatrixStart;
  is >> nbMatrixStart;
  struct KeyData {
    int idxMatrix;
    int nbAdjacent;
    int nbProcessed;
  };
  // int StatusTreatedForm; // 0: untreated, 1: treated but status not written on disk, 2: done and treated
  //
  // The list of requests.
  //
  std::vector<boost::mpi::request> ListRequest(MaxNumberFlyingMessage);
  std::vector<int> RequestStatus(MaxNumberFlyingMessage, 0);
  auto GetFreeIndex=[&]() -> int {
    for (int u=0; u<MaxNumberFlyingMessage; u++) {
      if (RequestStatus[u] == 0)
	return u;
      boost::optional<boost::mpi::status> stat = ListRequest[u].test();
      if (stat) { // that request has ended. Let's read it.
	if (stat->error() != 0) {
	  std::cerr << "something went wrong in the MPI" << std::endl;
	  throw TerminalException{1};
	}
	RequestStatus[u] = 0;
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
    TypeCtypeExch<Tint> ePerfect = ePair.ePerfect;
    auto it1 = ListCasesDone.find(ePerfect);
    if (it1 != ListCasesDone.end()) {
      log << "Processed entry=" << ePair.eIndex << "END" << std::endl;
      return;
    }
    KeyData& eData = ListCasesNotDone[ePerfect];
    if (eData.idxMatrix != 0) {
      log << "Processed entry=" << ePair.eIndex << "END" << std::endl;
      return;
    }
    eData.idxMatrix = idxMatrixCurrent + 1;
    log << "Inserting New perfect form" << ePair.ePerfect << " idxMatrixCurrent=" << idxMatrixCurrent << " Obtained from " << ePair.eIndex << "END" << std::endl;
    std::cerr << "Inserting new form, now we have |ListCasesNotDone[pos]|=" << ListCasesNotDone.size() << " |ListCasesDone|=" << ListCasesDone.size() << "\n";
    std::cerr << "idxMatrixCurrent=" << idxMatrixCurrent << " ePerfect = " << ePair.ePerfect << "\n";
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
  auto fSendMatrix=[&](PairExch<Tint> const& ePair, int const& u) -> void {
    int res=IntegerDiscriminantInvariant(ePair.ePerfect.eMat, n_pes);
    ListRequest[u] = world.isend(res, tag_new_form, ePair);
    RequestStatus[u] = 1;
  };
  std::vector<PairExch<Tint>> ListMatrixUnsent;
  auto ClearUnsentAsPossible=[&]() -> void {
    int pos=ListMatrixUnsent.size() - 1;
    while(true) {
      if (pos == -1)
	break;
      int idx = GetFreeIndex();
      if (idx == -1)
	break;
      fSendMatrix(ListMatrixUnsent[pos], idx);
      ListMatrixUnsent.pop_back();
      pos--;
    }
  };
  auto fInsertUnsent=[&](PairExch<Tint> const& ePair) -> void {
    int res=IntegerDiscriminantInvariant(ePair.ePerfect.eMat, n_pes);
    if (res == irank) {
      fInsert(ePair);
    }
    else {
      ListMatrixUnsent.push_back(ePair);
      ClearUnsentAsPossible();
    }
  };
  int nbCaseNotDone=0;
  for (int iMatStart=0; iMatStart<nbMatrixStart; iMatStart++) {
    int eStatus;
    is >> eStatus;
    MyMatrix<Tint> TheMat = ReadMatrix<Tint>(is);
    TypeCtypeExch<Tint> eRecMat{TheMat};
    int res=IntegerDiscriminantInvariant(TheMat, n_pes);
    if (res == irank) {
      KeyData eData{idxMatrixCurrent+1};
      if (eStatus == 0) {
        ListCasesNotDone[eRecMat] = eData;
        nbCaseNotDone++;
      }
      else {
        ListCasesDone[eRecMat] = eData;
      }
      log << "Reading existing matrix=" << eRecMat << " idxMatrixCurrent=" << idxMatrixCurrent << "END" << std::endl;
      idxMatrixCurrent++;
    }
  }
  std::cerr << "Reading finished, we have |ListCasesDone|=" << ListCasesDone.size() << " nbCaseNotDone=" << nbCaseNotDone << "\n";
  std::cerr << " |ListCasesNotDone|=" << ListCasesNotDone.size() << "\n";
  //
  // The main loop itself.
  //
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
  while(true) {
    boost::optional<boost::mpi::status> prob = world.iprobe();
    if (prob) {
      std::cerr << "We are probing something\n";
      if (prob->tag() == tag_new_form) {
	PairExch<Tint> ePair;
	world.recv(prob->source(), prob->tag(), ePair);
        fInsert(ePair);
      }
    }
    else {
      std::cerr << "irank=" << irank << " |ListMatrixUnsent|=" << ListMatrixUnsent.size() << " MaxStoredUnsentMatrices=" << MaxStoredUnsentMatrices << "\n";
      if (int(ListMatrixUnsent.size()) < MaxStoredUnsentMatrices) {
	boost::optional<std::pair<TypeCtypeExch<Tint>,int>> eReq=GetUndoneEntry();
	if (eReq) {
          std::cerr << "irank=" << irank << " eReq is non zero\n";
	  SetMatrixAsDone(eReq->first);
          std::cerr << "irank=" << irank << " ePerfect=" << eReq->first << "\n";
          int idxMatrixF = eReq->second;
          std::cerr << "irank=" << irank << " Starting Adjacent Form Method\n";
          std::vector<TypeCtypeExch<Tint>> ListAdjacentObject = GetAdjacentObjects<T,Tint>(eReq->first);
          int nbAdjacent = ListAdjacentObject.size();
          log << "Number of Adjacent for idxMatrixF=" << idxMatrixF << " nbAdjacent=" << nbAdjacent << " END" << std::endl;
          std::cerr << "irank=" << irank << " Number of Adjacent for idxMatrixF=" << idxMatrixF << " nbAdjacent=" << nbAdjacent << " END\n";
          int iAdj=0;
	  for (auto & eObj1 : ListAdjacentObject) {
            TypeIndex eIndex{irank, idxMatrixF, iAdj};
            PairExch<Tint> ePair{eObj1, eIndex};
	    fInsertUnsent(ePair);
            iAdj++;
	  }
	}
      }
    }
    std::cerr << "irank=" << irank << " Before ClearUnsentAsPossible\n";
    ClearUnsentAsPossible();
    std::cerr << "irank=" << irank << " After ClearUnsentAsPossible\n";
    //
    // Checking for termination of the program
    //
    std::chrono::time_point<std::chrono::system_clock> curr = std::chrono::system_clock::now();
    int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(curr - start).count();
    if (MaxRunTimeSecond  > 0) {
      if (elapsed_seconds > MaxRunTimeSecond) {
        std::cerr << "Exiting because the runtime is higher than the one expected\n";
        break;
      }
    }
  }
  std::cerr << "Normal termination of the program\n";
}

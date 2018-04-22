#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
namespace mpi = boost::mpi;

#include "MAT_Matrix.h"

static tag_new_form = 37;


int main()
{
  using T=mpq_class;
  using Tint=mpz_class;
  struct TypePerfectExch {
    int incd; // the number of shortest vectors divided by 2
    MyMatrix<T> eMat;
  };
  //
  FullNamelist eFull = NAMELIST_GetStandard_ENUMERATE_PERFECT();
  std::string eFileName = "perfectenum.nml";
  NAMELIST_ReadNamelistFile(eFileName, eFull);
  SingleBlock BlDATA = eFull.ListBlock["DATA"];
  int n=BlDATA.ListIntValues.at("n");
  int MaxNumberFlyingMessage=BlDATA.ListIntValues.at("MaxNumerFlyingMessage");
  int MaxIncidenceTreating=BlDATA.ListIntValues.at("MaxIncidenceTreating");
  int MaxStoredUnsentMatrices=BlDATA.ListIntValues.at("MaxStoredUnsentMatrices");
  //
  boost::mpi::environment env;
  boost::mpi::communicator world;
  int irank=world.rank();
  int size=world.size();
  std::string eFileO="LOG_" + irank;
  std::ofstream log(eFileO);
  log << "Initial log entry\n";
  //
  std::string FileMatrix = "ListMatrix";
  std::ifstream is(FileMatrix);
  int nbMatrixStart;
  is >> nbMatrixStart;
  struct KeyData {
    int StatusTreatedForm; // 0: untreated, 1: treated but status not written on disk, 2: done and treated
  };
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
	  std::cerr << "something went wrong in the MPI\n";
	  throw TerminalExeption{1};
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
  std::map<TypePerfectExch,KeyData> ListCasesNotDone;
  std::map<TypePerfectExch,KeyData> ListCasesDone;
  auto fInsert=[&](TypePerfectExch const& NewMat) -> void {
    auto it1 = ListCasesDone.find(NewMat);
    if (it1 != ListCasesDone.end())
      return;
    auto it2 = ListCasesNotDone.find(NewMat);
    if (it2 != ListCasesNotDone.end())
      return;
    ListCasesNotDone[NewMat] = 0;
    os << "Inserting new form, now we have |ListCases|=" << ListCases.size() << "\n";
  };
  auto GetLowestIncidenceUndone=[&]() -> boost::optional<TypePerfectExch> {
    auto it1 = ListCasesNotDone.begin();
    if (it1 == ListCasesNotDone.end())
      return {};
    if (it1->incd > MaxIncidenceTreating)
      return {};
    return boost::optional(*it1);
  };
  auto SetMatrixAsDone=[&](TypePerfectExch const& TheMat) -> void {
    ListCasesNotDone.erase(TheMat);
    ListCasesDone[TheMat] = 1;
  };
  //
  // The system for sending matrices
  //
  auto fSendMatrix=[&](TypePerfectMatrix const& NewMat, int const& u) -> void {
    int KeyInv=IntegerDiscriminantInvariant(NewMat);
    int res=KeyInv % size;
    ListRequest[u] = world.isend(res, tag_new_form, NewMat);
    RequestStatus[u] = 1;
  };
  std::vector<TypePerfectMatrix> ListMatrixUnsent;
  auto ClearUnsentAsPossible=[&]() -> void {
    int pos=ListMatrixUnsent.size() - 1;
    while(true) {
      if (pos == -1)
	break;
      int idx = GetFreeIndex();
      if (idx == -1)
	break;
      fSendMatrix(ListMatrixUnsent[pos]);
      ListMatrixUnsent.pop_back();
      pos--;
    }
  };
  auto fInsertUnsent=[&](TypePerfectMatrix const& eMat) -> void {
    int KeyInv=IntegerDiscriminantInvariant(eMat);
    int res=KeyInv % size;
    if (res == irank) {
      fInsert(eMat);
    }
    else {
      ListMatrixUnsent.push_back(eMat);
      ClearUnsentAsPossible();
    }
  };
  for (int iMatStart=0; iMatStart<nbMatrixStart; iMatStart++) {
    int eStatus;
    is >> eStatus;
    MyMatrix<Tint> TheMat = ReadMatrix<Tint>(is);
    KeyData eData{eStatus,TheMat};
    int KeyInv=IntegerDiscriminantInvariant(TheMat);
    int res=KeyInv % size;
    if (res == irank)
      ListCases.push_back(eData);
  }
  os << "Reading finished, we have |ListCases|=" << ListCases.size() << "\n";
  //
  // The main loop itself.
  //
  while(true) {
    boost::optional<boost::mpi::status> prob = world.iprobe();
    if (prob) {
      if (prob->tag() == tag_new_form) {
	std::vector<MyMatrix<Tint>> data;
	world.recv(msg->source(), msg->tag(), data);
	for (auto & eMat : data)
	  fInsert(eMat);
      }
    }
    else {
      if (ListMatrixUnsent.size() < MaxStoredUnsentMatrices) {
	boost::optional<TypePerfectExch> eReq=GetLowestIncidenceUndone();
	if (eReq) {
	  SetMatrixAsDone(*eReq);
	  std::vector<MyMatrix<Tint>> ListAdjacent = GetAdjacentFormDirectMethod(eReq->eMat);
	  for (auto & eMat : ListAdjacent) {
	    MyMatrix<T> eMatCan = ComputeCanonicalForm<T,Tint>(eMat);
	    Tshortest<T,Tint> eRec = T_ShortestVector(eMatCan);
	    int incd = (eRec.SHV.rows()) / 2;
	    TypePerfectMatrix RecMat{incd, eMatCan};
	    fInsertUnsent(RecMat);
	  }
	}
      }
    }
    ClearUnsentAsPossible();
  }
  
}

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
    MyMatrix<Tint> eMat;
    int incd; // the number of shortest vectors divided by 2
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
  std::map<TypePerfectExch,KeyData> ListCases;
  auto fInsert=[&](TypePerfectExch const& NewMat) -> void {
    auto it = ListCases.find(NewMat);
    if (it != ListCases.end())
      return;
    ListCases[NewMat] = 0;
    os << "Inserting new form, now we have |ListCases|=" << ListCases.size() << "\n";
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
  auto fInsertUnsent=[&](TypePerfectMatrix const& eMat) -> void {
    ListMatrixUnsent.push_back(eMat);
  };
  auto ClearUnsentAsPossible=[&]() -> void {
    int pos=ListMatrixUnsent.size() - 1;
    while(true) {
      if (pos == -1)
	break;
      int idx = GetFreeIndex();
      if (idx == -1)
	break;
      fSendMatrix(ListMatrixUnsent[pos]);
      pos--;
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
      

      
      std::vector<MyMatrix<Tint>> ListAdjacent = 
    }
    ClearUnsentAsPossible();
  }
  
}

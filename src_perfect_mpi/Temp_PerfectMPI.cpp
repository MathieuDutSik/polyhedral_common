#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "Namelist.h"
#include "MatrixCanonicalForm.h"


#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
namespace mpi = boost::mpi;


FullNamelist NAMELIST_GetStandard_ENUMERATE_PERFECT_MPI()
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
  ListIntValues1["MaxIncidenceTreating"]=45 + 20;
  ListIntValues1["MaxStoredUnsentMatrices"]=1000;
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

template<typename T>
int IntegerDiscriminantInvariant(MyMatrix<T> const& NewMat)
{
  T TheDet=DeterminantMat(NewMat);
  int TheDet_i = UniversalTypeConversion<int,T>(TheDet);
  return TheDet_i;
}




static int tag_new_form = 37;


int main()
{
  using T=mpq_class;
  using Tint=mpz_class;
  struct TypePerfectExch {
    int incd; // the number of shortest vectors divided by 2
    MyMatrix<Tint> eMat;
  };
  //
  FullNamelist eFull = NAMELIST_GetStandard_ENUMERATE_PERFECT_MPI();
  std::string eFileName = "perfectenum.nml";
  NAMELIST_ReadNamelistFile(eFileName, eFull);
  SingleBlock BlDATA = eFull.ListBlock["DATA"];
  int n=BlDATA.ListIntValues.at("n");
  int MaxNumberFlyingMessage=BlDATA.ListIntValues.at("MaxNumberFlyingMessage");
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
  std::map<TypePerfectExch,KeyData> ListCasesNotDone;
  std::map<TypePerfectExch,KeyData> ListCasesDone;
  auto fInsert=[&](TypePerfectExch const& NewMat) -> void {
    auto it1 = ListCasesDone.find(NewMat);
    if (it1 != ListCasesDone.end())
      return;
    auto it2 = ListCasesNotDone.find(NewMat);
    if (it2 != ListCasesNotDone.end())
      return;
    ListCasesNotDone[NewMat] = {0};
    log << "Inserting new form, now we have |ListCasesNotDone|=" << ListCasesNotDone.size() << " |ListCasesDone|=" << ListCasesDone.size() << "\n";
  };
  auto GetLowestIncidenceUndone=[&]() -> boost::optional<TypePerfectExch> {
    auto it1 = ListCasesNotDone.begin();
    if (it1 == ListCasesNotDone.end())
      return {};
    if (it1->first.incd > MaxIncidenceTreating)
      return {};
    return boost::optional<TypePerfectExch>(it1->first);
  };
  auto SetMatrixAsDone=[&](TypePerfectExch const& TheMat) -> void {
    ListCasesNotDone.erase(TheMat);
    ListCasesDone[TheMat] = {1};
  };
  //
  // The system for sending matrices
  //
  auto fSendMatrix=[&](TypePerfectExch const& eRecMat, int const& u) -> void {
    int KeyInv=IntegerDiscriminantInvariant(eRecMat.eMat);
    int res=KeyInv % size;
    ListRequest[u] = world.isend(res, tag_new_form, eRecMat);
    RequestStatus[u] = 1;
  };
  std::vector<TypePerfectExch> ListMatrixUnsent;
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
  auto fInsertUnsent=[&](TypePerfectExch const& eRecMat) -> void {
    int KeyInv=IntegerDiscriminantInvariant(eRecMat.eMat);
    int res=KeyInv % size;
    if (res == irank) {
      fInsert(eRecMat);
    }
    else {
      ListMatrixUnsent.push_back(eRecMat);
      ClearUnsentAsPossible();
    }
  };
  for (int iMatStart=0; iMatStart<nbMatrixStart; iMatStart++) {
    int eStatus;
    is >> eStatus;
    int incd;
    is >> incd;
    MyMatrix<Tint> TheMat = ReadMatrix<Tint>(is);
    TypePerfectExch eRecMat{incd, TheMat};
    KeyData eData{eStatus};
    int KeyInv=IntegerDiscriminantInvariant(TheMat);
    int res=KeyInv % size;
    if (res == irank) {
      if (eStatus == 0) {
        ListCasesNotDone[eRecMat] = eData;
      }
      else {
        ListCasesDone[eRecMat] = eData;
      }
    }
  }
  log << "Reading finished, we have |ListCasesDone|=" << ListCasesDone.size() << " |ListCasesNotDone|=" << ListCasesNotDone.size() << "\n";
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
      if (int(ListMatrixUnsent.size()) < MaxStoredUnsentMatrices) {
	boost::optional<TypePerfectExch> eReq=GetLowestIncidenceUndone();
	if (eReq) {
	  SetMatrixAsDone(*eReq);
	  std::vector<MyMatrix<Tint>> ListAdjacent = GetAdjacentFormDirectMethod(eReq->eMat);
	  for (auto & eMat : ListAdjacent) {
	    MyMatrix<T> eMatCan = ComputeCanonicalForm<T,Tint>(eMat);
	    Tshortest<T,Tint> eRec = T_ShortestVector(eMatCan);
	    int incd = (eRec.SHV.rows()) / 2;
	    TypePerfectExch RecMat{incd, eMatCan};
	    fInsertUnsent(RecMat);
	  }
	}
      }
    }
    ClearUnsentAsPossible();
  }
}

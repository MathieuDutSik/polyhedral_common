#include "PerfectMPI_types.h"
#include "NumberTheory.h"
#include "Namelist.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PerfectForm.h"
#include <unordered_map>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include "md5sum.h"
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
  ListIntValues1["MinIncidenceRealized"]=0;
  ListIntValues1["MaxIncidenceRealized"]=1000;
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
std::vector<TypePerfectExch<Tint>> GetAdjacentObjects(TypePerfectExch<Tint> const& eObjIn)
{
  MyMatrix<T> eMat_T = ConvertMatrixUniversal<T,Tint>(eObjIn.eMat);
  Tshortest<T,Tint> eRec = T_ShortestVector<T,Tint>(eMat_T);
  int n=eRec.SHV.cols();
  int nbShort=eRec.SHV.rows() / 2;
  int dimSymm=n*(n+1)/2;
  MyMatrix<Tint> SHVred(nbShort, n);
  for (int iShort=0; iShort<nbShort; iShort++)
    for (int i=0; i<n; i++)
      SHVred(iShort,i) = eRec.SHV(2*iShort,i);
  MyMatrix<T> ConeClassical = GetNakedPerfectConeClassical<T,Tint>(SHVred);
  std::vector<Face> ListIncd = lrs::DualDescription_temp_incd(ConeClassical);
  MyVector<T> Wvect = GetSymmetricMatrixWeightVector<T>(n);
  std::vector<TypePerfectExch<Tint>> ListAdjMat;
  for (auto & eIncd : ListIncd) {
    MyVector<T> eFacet=FindFacetInequality(ConeClassical, eIncd);
    MyVector<T> Vexpand(dimSymm);
    for (int i=0; i<dimSymm; i++)
      Vexpand(i) = eFacet(i) / Wvect(i);
    MyMatrix<T> eMatDir = VectorToSymmetricMatrix(Vexpand, n);
    std::pair<MyMatrix<T>,Tshortest<T,Tint>> ePairAdj = Flipping_Perfect<T,Tint>(eMat_T, eMatDir);
    int incd = ePairAdj.second.SHV.rows() / 2;
    //
    MyMatrix<T> eMat2 = ComputeCanonicalForm<T,Tint>(ePairAdj.first).Mat;
    MyMatrix<T> eMat3 = RemoveFractionMatrix(eMat2);
    MyMatrix<Tint> eMat4 = ConvertMatrixUniversal<Tint,T>(eMat3);
    TypePerfectExch<Tint> RecMat{incd, eMat4};
    ListAdjMat.emplace_back(RecMat);
  }
  return ListAdjMat;
}





template<typename T>
int IntegerDiscriminantInvariant(MyMatrix<T> const& NewMat, int const& size)
{
  int TheChoice=3;
  if (TheChoice == 1) {
    T TheDet=DeterminantMat(NewMat);
    int TheDet_i = UniversalTypeConversion<int,T>(TheDet);
    return TheDet_i % size;
  }
  if (TheChoice == 2) {
    int nRow=NewMat.rows();
    int nCol=NewMat.cols();
    int TheInvariant=0;
    for (int iRow=0; iRow<nRow; iRow++)
      for (int iCol=0; iCol<nCol; iCol++) {
        int eCoeff1 = 3 + 5*iRow + 7*iCol;
        int eVal=NewMat(iRow,iCol);
        int eCoeff2;
        if (eVal < 0) {
          eCoeff2 = -2*eVal;
        }
        else {
          eCoeff2 = 2*eVal + 1;
        }
        TheInvariant += eCoeff1 * eCoeff2;
        std::cerr << "iRow=" << iRow << " iCol=" << iCol << " eCoeff1=" << eCoeff1 << " eCoeff2=" << eCoeff2 << " eVal=" << eVal << " TheInvariant=" << TheInvariant << "\n";
      }
    std::cerr << "TheInvariant=" << TheInvariant << "\n";
    return TheInvariant % size;
  }
  if (TheChoice == 3) {
    std::stringstream s;
    int nbRow=NewMat.rows();
    int nbCol=NewMat.cols();
    for (int iRow=0; iRow<nbRow; iRow++)
      for (int iCol=0; iCol<nbCol; iCol++)
        s << " " << NewMat(iRow, iCol);
    std::string converted(s.str());
    mpz_class e_hash_mpz = MD5_hash_mpz(converted);
    mpz_class e_modulo = size;
    mpz_class e_res = e_hash_mpz % e_modulo;
    int residue = UniversalTypeConversion<int,mpz_class>(e_res);
    return residue;
  }
  return 0;
}







static int tag_new_form = 37;


int main()
{
  using T=mpq_class;
  using Tint=long;
  //
  FullNamelist eFull = NAMELIST_GetStandard_ENUMERATE_PERFECT_MPI();
  std::string eFileName = "perfectenum.nml";
  NAMELIST_ReadNamelistFile(eFileName, eFull);
  SingleBlock BlDATA = eFull.ListBlock["DATA"];
  //  int n=BlDATA.ListIntValues.at("n");
  int MaxNumberFlyingMessage = BlDATA.ListIntValues.at("MaxNumberFlyingMessage");
  int MaxIncidenceTreating = BlDATA.ListIntValues.at("MaxIncidenceTreating");
  int MaxStoredUnsentMatrices = BlDATA.ListIntValues.at("MaxStoredUnsentMatrices");
  int MinIncidenceRealized = BlDATA.ListIntValues.at("MinIncidenceRealized");
  int MaxIncidenceRealized = BlDATA.ListIntValues.at("MaxIncidenceRealized");
  int MaxRunTimeSecond = BlDATA.ListIntValues.at("MaxRunTimeSecond");
  std::string FileMatrix = BlDATA.ListStringValues.at("ListMatrixInput");
  //
  boost::mpi::environment env;
  boost::mpi::communicator world;
  int irank=world.rank();
  int size=world.size();
  std::string eFileO="LOG_" + IntToString(irank);
  std::ofstream log(eFileO);
  log << "Initial log entry" << std::endl;
  //
  std::ifstream is(FileMatrix);
  int nbMatrixStart;
  is >> nbMatrixStart;
  struct KeyData {
    int idxMatrix;
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
  int nbCaseIncidence=MaxIncidenceRealized + 1 - MinIncidenceRealized;
  std::vector<std::unordered_map<TypePerfectExch<Tint>,KeyData>> ListCasesNotDone(nbCaseIncidence);
  std::unordered_map<TypePerfectExch<Tint>,KeyData> ListCasesDone;
  int idxMatrixCurrent=0;
  auto fInsert=[&](PairExch<Tint> const& ePair) -> void {
    TypePerfectExch<Tint> ePerfect = ePair.ePerfect;
    auto it1 = ListCasesDone.find(ePerfect);
    if (it1 != ListCasesDone.end()) {
      log << "Processed entry=" << ePair.eIndex << "END" << std::endl;
      return;
    }
    int pos = ePerfect.incd - MinIncidenceRealized;
    auto it2 = ListCasesNotDone[pos].find(ePerfect);
    if (it2 != ListCasesNotDone[pos].end()) {
      log << "Processed entry=" << ePair.eIndex << "END" << std::endl;
      return;
    }
    ListCasesNotDone[pos][ePerfect] = {idxMatrixCurrent};
    log << "Inserting New perfect form" << ePair.ePerfect << " idxMatrixCurrent=" << idxMatrixCurrent << " Obtained from " << ePair.eIndex << "END" << std::endl;
    std::cerr << "Inserting new form, now we have pos=" << pos << " |ListCasesNotDone[pos]|=" << ListCasesNotDone[pos].size() << " |ListCasesDone|=" << ListCasesDone.size() << "\n";
    std::cerr << "idxMatrixCurrent=" << idxMatrixCurrent << " ePerfect = " << ePair.ePerfect << "\n";
    idxMatrixCurrent++;
  };
  auto GetLowestIncidenceUndone=[&]() -> boost::optional<std::pair<TypePerfectExch<Tint>,int>> {
    for (int iCaseIncidence=0; iCaseIncidence<nbCaseIncidence; iCaseIncidence++) {
      int incd = iCaseIncidence + MinIncidenceRealized;
      if (incd <= MaxIncidenceTreating) {
        auto it1 = ListCasesNotDone[iCaseIncidence].begin();
        if (it1 != ListCasesNotDone[iCaseIncidence].end()) {
          std::pair<TypePerfectExch<Tint>,int> ePair = {it1->first, it1->second.idxMatrix};
          return boost::optional<std::pair<TypePerfectExch<Tint>,int>>(ePair);
        }
      }
    }
    return {};
  };
  auto SetMatrixAsDone=[&](TypePerfectExch<Tint> const& TheMat) -> void {
    int pos = TheMat.incd - MinIncidenceRealized;
    KeyData eKey = ListCasesNotDone[pos].at(TheMat);
    ListCasesNotDone[pos].erase(TheMat);
    ListCasesDone[TheMat] = eKey;
  };
  //
  // The system for sending matrices
  //
  auto fSendMatrix=[&](PairExch<Tint> const& ePair, int const& u) -> void {
    int res=IntegerDiscriminantInvariant(ePair.ePerfect.eMat, size);
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
    int res=IntegerDiscriminantInvariant(ePair.ePerfect.eMat, size);
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
    int incd;
    is >> incd;
    MyMatrix<Tint> TheMat = ReadMatrix<Tint>(is);
    TypePerfectExch<Tint> eRecMat{incd, TheMat};
    int res=IntegerDiscriminantInvariant(TheMat, size);
    if (res == irank) {
      KeyData eData{idxMatrixCurrent};
      if (eStatus == 0) {
        int pos = incd - MinIncidenceRealized;
        ListCasesNotDone[pos][eRecMat] = eData;
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
  for (int iCaseIncidence=0; iCaseIncidence<nbCaseIncidence; iCaseIncidence++) {
    int incd = iCaseIncidence + MinIncidenceRealized;
    std::cerr << "iCaseIncidence=" << iCaseIncidence << " incd=" << incd << " |ListCasesNotDone[pos]|=" << ListCasesNotDone[iCaseIncidence].size() << "\n";
  }
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
	boost::optional<std::pair<TypePerfectExch<Tint>,int>> eReq=GetLowestIncidenceUndone();
	if (eReq) {
          std::cerr << "irank=" << irank << " eReq is non zero\n";
	  SetMatrixAsDone(eReq->first);
          std::cerr << "irank=" << irank << " ePerfect=" << eReq->first << "\n";
          MyMatrix<T> eMat_T = ConvertMatrixUniversal<T,Tint>(eReq->first.eMat);
          int idxMatrixF = eReq->second;
          std::cerr << "irank=" << irank << " Starting Adjacent Form Method\n";
          std::vector<TypePerfectExch<Tint>> ListAdjacentObject = GetAdjacentObjects<T,Tint>(eReq->first);
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

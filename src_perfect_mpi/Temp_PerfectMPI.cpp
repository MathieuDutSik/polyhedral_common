#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "Namelist.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PerfectForm.h"


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
std::vector<MyMatrix<T>> GetAdjacentFormDirectMethod(MyMatrix<T> const& eMatIn)
{
  Tshortest<T,Tint> eRec = T_ShortestVector<T,Tint>(eMatIn);
  int n=eRec.SHV.cols();
  int nbShort=eRec.SHV.rows() / 2;
  int dimSymm=n*(n+1)/2;
  MyMatrix<Tint> SHVred(nbShort, n);
  for (int iShort=0; iShort<nbShort; iShort++) {
    for (int i=0; i<n; i++)
      SHVred(iShort,i) = eRec.SHV(2*iShort,i);
  }
  MyMatrix<T> ConeClassical = GetNakedPerfectConeClassical<T,Tint>(SHVred);
  std::vector<Face> ListIncd = lrs::DualDescription_temp_incd(ConeClassical);
  MyVector<T> Wvect = GetSymmetricMatrixWeightVector<T>(n);
  std::vector<MyMatrix<T>> ListAdjMat;
  for (auto & eIncd : ListIncd) {
    MyVector<T> eFacet=FindFacetInequality(ConeClassical, eIncd);
    MyVector<T> Vexpand(dimSymm);
    for (int i=0; i<dimSymm; i++)
      Vexpand(i) = eFacet(i) / Wvect(i);
    MyMatrix<T> eMatDir=VectorToSymmetricMatrix(Vexpand, n);
    MyMatrix<T> eMatAdj = Flipping_Perfect(eMatIn, eMatDir);
    ListAdjMat.push_back(eMatAdj);
  }
  return ListAdjMat;
}



struct TypeIndex {
  int iProc;
  int idxMatrix;
  int iAdj;
};


template<typename T>
struct TypePerfectExch {
  int incd; // the number of shortest vectors divided by 2
  MyMatrix<T> eMat;
};


template<typename T>
struct PairExch {
  TypePerfectExch<T> ePerfect;
  TypeIndex eIndex;
};


template<typename T>
std::ostream& operator<<(std::ostream& os, TypePerfectExch<T> const& obj)
{
  os << obj.incd;
  int nbRow=obj.eMat.rows();
  for (int iRow=0; iRow<nbRow; iRow++) {
    for (int iCol=0; iCol<nbRow; iCol++)
      os << " " << obj.eMat(iRow, iCol);
    os << "\n";
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, TypeIndex const& obj)
{
  os << obj.iProc << " , " << obj.idxMatrix << " , " << obj.iAdj;
  return os;
}






namespace std {
  template<typename T>
  struct less<TypePerfectExch<T>> {
    bool operator()(TypePerfectExch<T> const& eTPE1, TypePerfectExch<T> const& eTPE2) const
    {
      if (eTPE1.incd < eTPE2.incd)
        return true;
      if (eTPE1.incd > eTPE2.incd)
        return false;
      //
      int nbRow=eTPE1.eMat.rows();
      for (int iRow=0; iRow<nbRow; iRow++)
        for (int iCol=0; iCol<nbRow; iCol++) {
          if (eTPE1.eMat(iRow,iCol) < eTPE2.eMat(iRow,iCol))
            return true;
          if (eTPE1.eMat(iRow,iCol) > eTPE2.eMat(iRow,iCol))
            return false;
        }
      return false;
    }
  };
}


template<typename T>
int IntegerDiscriminantInvariant(MyMatrix<T> const& NewMat)
{
  T TheDet=DeterminantMat(NewMat);
  int TheDet_i = UniversalTypeConversion<int,T>(TheDet);
  return TheDet_i;
}

namespace boost { namespace serialization {
    // TypePerfectExch
    template<class Archive, typename T>
      inline void serialize(Archive & ar,
                            TypePerfectExch<T> & eRecMat,
                            const unsigned int version)
      {
        ar & make_nvp("incd", eRecMat.incd);
        int rows = eRecMat.eMat.rows();
        int cols = eRecMat.eMat.cols();
        ar & make_nvp("rows", rows);
        ar & make_nvp("cols", cols);
        eRecMat.eMat.resize(rows, cols);
        for(int r = 0; r < rows; ++r)
          for(int c = 0; c < cols; ++c)
            ar & make_nvp("val", eRecMat.eMat(r,c));
      }

    // TypePerfectExch
    template<class Archive>
      inline void serialize(Archive & ar,
                            TypeIndex & eTypIdx,
                            const unsigned int version)
      {
        ar & make_nvp("iProc", eTypIdx.iProc);
        ar & make_nvp("idxMatrix", eTypIdx.idxMatrix);
        ar & make_nvp("iAdj", eTypIdx.iAdj);
      }

    // PairExch
    template<class Archive, typename T>
      inline void serialize(Archive & ar,
                            PairExch<T> & ePair,
                            const unsigned int version)
      {
        ar & make_nvp("perfect", ePair.ePerfect);
        ar & make_nvp("index"  , ePair.eIndex);
      }


  }}






static int tag_new_form = 37;


int main()
{
  using T=mpq_class;
  using Tint=int;
  //
  FullNamelist eFull = NAMELIST_GetStandard_ENUMERATE_PERFECT_MPI();
  std::string eFileName = "perfectenum.nml";
  NAMELIST_ReadNamelistFile(eFileName, eFull);
  SingleBlock BlDATA = eFull.ListBlock["DATA"];
  //  int n=BlDATA.ListIntValues.at("n");
  int MaxNumberFlyingMessage = BlDATA.ListIntValues.at("MaxNumberFlyingMessage");
  int MaxIncidenceTreating = BlDATA.ListIntValues.at("MaxIncidenceTreating");
  int MaxStoredUnsentMatrices = BlDATA.ListIntValues.at("MaxStoredUnsentMatrices");
  std::string FileMatrix = BlDATA.ListStringValues.at("ListMatrixInput");
  //
  boost::mpi::environment env;
  boost::mpi::communicator world;
  int irank=world.rank();
  int size=world.size();
  std::string eFileO="LOG_" + irank;
  std::ofstream log(eFileO);
  log << "Initial log entry\n";
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
  std::map<TypePerfectExch<Tint>,KeyData> ListCasesNotDone;
  std::map<TypePerfectExch<Tint>,KeyData> ListCasesDone;
  int idxMatrixCurrent=0;
  auto fInsert=[&](PairExch<Tint> const& ePair) -> void {
    TypePerfectExch<Tint> ePerfect = ePair.ePerfect;
    auto it1 = ListCasesDone.find(ePerfect);
    if (it1 != ListCasesDone.end())
      return;
    auto it2 = ListCasesNotDone.find(ePerfect);
    if (it2 != ListCasesNotDone.end())
      return;
    ListCasesNotDone[ePerfect] = {idxMatrixCurrent};
    idxMatrixCurrent++;
    log << "Inserting New perfect form\n";
    log << ePair.ePerfect;
    log << "Obtained from " << ePair.eIndex << "\n";
    log << "COMMENT\n";
    log << "Inserting new form, now we have |ListCasesNotDone|=" << ListCasesNotDone.size() << " |ListCasesDone|=" << ListCasesDone.size() << "\n";
  };
  auto GetLowestIncidenceUndone=[&]() -> boost::optional<std::pair<TypePerfectExch<Tint>,int>> {
    auto it1 = ListCasesNotDone.begin();
    if (it1 == ListCasesNotDone.end())
      return {};
    if (it1->first.incd > MaxIncidenceTreating)
      return {};
    std::pair<TypePerfectExch<Tint>,int> ePair {it1->first, it1->second.idxMatrix};
    return boost::optional<std::pair<TypePerfectExch<Tint>,int>>(ePair);
  };
  auto SetMatrixAsDone=[&](TypePerfectExch<Tint> const& TheMat) -> void {
    KeyData eKey = ListCasesNotDone.at(TheMat);
    ListCasesNotDone.erase(TheMat);
    ListCasesDone[TheMat] = eKey;
  };
  //
  // The system for sending matrices
  //
  auto fSendMatrix=[&](PairExch<Tint> const& ePair, int const& u) -> void {
    int KeyInv=IntegerDiscriminantInvariant(ePair.ePerfect.eMat);
    int res=KeyInv % size;
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
    int KeyInv=IntegerDiscriminantInvariant(ePair.ePerfect.eMat);
    int res=KeyInv % size;
    if (res == irank) {
      fInsert(ePair);
    }
    else {
      ListMatrixUnsent.push_back(ePair);
      ClearUnsentAsPossible();
    }
  };
  for (int iMatStart=0; iMatStart<nbMatrixStart; iMatStart++) {
    int eStatus;
    is >> eStatus;
    int incd;
    is >> incd;
    MyMatrix<Tint> TheMat = ReadMatrix<Tint>(is);
    TypePerfectExch<Tint> eRecMat{incd, TheMat};
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
	PairExch<Tint> ePair;
	world.recv(prob->source(), prob->tag(), ePair);
        fInsert(ePair);
      }
    }
    else {
      if (int(ListMatrixUnsent.size()) < MaxStoredUnsentMatrices) {
	boost::optional<std::pair<TypePerfectExch<Tint>,int>> eReq=GetLowestIncidenceUndone();
	if (eReq) {
	  SetMatrixAsDone(eReq->first);
          MyMatrix<T> eMat_T = ConvertMatrixUniversal<T,Tint>(eReq->first.eMat);
          int idxMatrixF = eReq->second;
	  std::vector<MyMatrix<T>> ListAdjacent = GetAdjacentFormDirectMethod<T,Tint>(eMat_T);
          log << "Number of Adjacent for idxMatrixF=" << idxMatrixF << " nbAdjacent=" << ListAdjacent.size() << "\n";
          int iAdj=0;
	  for (auto & eMat1 : ListAdjacent) {
	    MyMatrix<T> eMat2 = ComputeCanonicalForm<T,Tint>(eMat1).second;
	    Tshortest<T,Tint> eRec = T_ShortestVector<T,Tint>(eMat2);
	    int incd = (eRec.SHV.rows()) / 2;
            MyMatrix<T> eMat3 = RemoveFractionMatrix(eMat2);
            MyMatrix<Tint> eMat4 = ConvertMatrixUniversal<Tint,T>(eMat3);
	    TypePerfectExch<Tint> RecMat{incd, eMat4};
            TypeIndex eIndex{irank, idxMatrixF, iAdj};
            PairExch<Tint> ePair{RecMat, eIndex};
	    fInsertUnsent(ePair);
            iAdj++;
	  }
	}
      }
    }
    ClearUnsentAsPossible();
  }
}

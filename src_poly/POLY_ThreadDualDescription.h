#ifndef TEMP_THREAD_DUAL_DESCRIPTION_INCLUDE
#define TEMP_THREAD_DUAL_DESCRIPTION_INCLUDE

#include "Basic_string.h"
#include "Basic_file.h"
#include "Namelist.h"

#include "POLY_Heuristics.h"

#include "ThreadManagement.h"
#include "Parallel_Classes.h"


#include "POLY_DirectDualDesc.h"
#include "POLY_SamplingFacet.h"


#include "Temp_PolytopeEquiStab.h"
#include "MatrixGroup.h"

//
// PolyhedralInv
//

template<typename T>
struct PolyhedralInv {
  int dim;
  T eVal;
};

template<typename T>
std::istream& operator>>(std::istream& is, PolyhedralInv<T>& obj)
{
  int dim;
  T eVal;
  is >> dim;
  is >> eVal;
  //
  obj={dim, eVal};
  return is;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, PolyhedralInv<T> const& obj)
{
  os << " " << obj.dim;
  os << " " << obj.eVal;
  return os;
}

template<typename T>
bool operator==(PolyhedralInv<T> const& x, PolyhedralInv<T> const& y)
{
  if (x.dim != y.dim)
    return false;
  if (x.eVal != y.eVal)
    return false;
  return true;
}


//
// PolyhedralEntry
//

template<typename T, typename Tgroup>
struct PolyhedralEntry {
  MyMatrix<T> EXT;
  Tgroup GRP;
  std::vector<Face> ListFace;
};

// We only want the full symmetry in the entry.
// so as to avoid any ambiguity in the process
template<typename T, typename Tgroup>
PolyhedralEntry<T,Tgroup> CanonicalizationPolyEntry(PolyhedralEntry<T,Tgroup> const& eEnt, std::ostream & os)
{
  MyMatrix<T> EXTred=ColumnReduction(eEnt.EXT);
  WeightMatrix<T,T> WMat=GetWeightMatrix(EXTred);
  os << "Canonicalization, we have WMat\n";
  Tgroup GRPlin=GetStabilizerWeightMatrix<T,T,Tgroup>(WMat);
  os << "Canonicalization, GRPlin.size=" << GRPlin.size() << " eEnt.GRP.size=" << eEnt.GRP.size() << "\n";
  if (GRPlin.size() == eEnt.GRP.size())
    return eEnt;
  os << "|eEnt.ListFace|=" << eEnt.ListFace.size() << "\n";
  {
    std::ofstream os1("CANONIC_GRPlin");
    std::ofstream os2("CANONIC_GRPlin.gap");
    std::ofstream os3("CANONIC_ListFace");
    std::ofstream os4("CANONIC_EXT");
    WriteGroup(os1, GRPlin);
    WriteGroupGAP(os2, GRPlin);
    WriteListFace(os3, eEnt.ListFace);
    WriteMatrix(os4, EXTred);
  }
  WeightMatrix<int,int> WMatInt=WeightMatrixFromPairOrbits<int,int>(GRPlin, os);
  os << "We have WMatInt\n";
  LocalInvInfo LocalInv=ComputeLocalInvariantStrategy(WMatInt, GRPlin, "pairinv", os);
  os << "We have LocalInv\n";
  struct Local {
    Face eFace;
    std::vector<int> eInv;
  };
  std::vector<Local> ListLocal;
  auto FuncInsert=[&](Face const& eFace) -> void {
    std::vector<int> eInv=GetLocalInvariantWeightMatrix_Enhanced<int>(LocalInv, eFace);
    for (auto & eRec : ListLocal)
      if (eRec.eInv == eInv) {
	{
	  std::ofstream os1("CANONIC_PairFace");
	  std::ofstream os2("CANONIC_PairFace.gap");
	  WriteListFace(os1, {eFace, eRec.eFace});
	  WriteListFaceGAP(os2, {eFace, eRec.eFace});
	}
	bool test=TestEquivalence(GRPlin, eFace, eRec.eFace);
	if (test)
	  return;
      }
    ListLocal.push_back({eFace,eInv});
  };
  int iter=0;
  for (auto & eFace : eEnt.ListFace) {
    os << "iter=" << iter << " Now |ListLocal|=" << ListLocal.size() << "\n";
    FuncInsert(eFace);
    os << "After FuncInsert\n";
    iter++;
  }
  os << "After FuncInsert\n";
  std::vector<Face> RetListRepr;
  for (auto & eRec : ListLocal)
    RetListRepr.push_back(eRec.eFace);
  PolyhedralEntry<T,Tgroup> fEnt{eEnt.EXT, GRPlin, RetListRepr};
  return fEnt;
}

template<typename T,typename Tgroup>
std::istream& operator>>(std::istream& is, PolyhedralEntry<T,Tgroup>& obj)
{
  MyMatrix<T> EXT=ReadMatrix<T>(is);
  Tgroup GRP=ReadGroup<Tgroup>(is);
  std::vector<Face> ListFace=ReadListFace(is);
  PolyhedralInv<T> eInv;
  //
  obj={EXT, GRP, ListFace};
  return is;
}

template<typename T,typename Tgroup>
std::ostream& operator<<(std::ostream& os, PolyhedralEntry<T,Tgroup> const& eEnt)
{
  WriteMatrix<T>(os, eEnt.EXT);
  WriteGroup(os, eEnt.GRP);
  WriteListFace(os, eEnt.ListFace);
  return os;
}





//
// SimpleOrbitFacet
//

template<typename T>
struct SimpleOrbitFacet {
  Face eRepr;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, SimpleOrbitFacet<T> const& eEnt)
{
  WriteFace(os, eEnt.eRepr);
  return os;
}


template<typename T>
std::istream& operator>>(std::istream& is, SimpleOrbitFacet<T> & eEnt)
{
  Face eSet=ReadFace(is);
  eEnt={eSet};
  return is;
}

//
// SimpleOrbitFacetInv related defines
//

using IntInv=short;

template<typename T>
struct SimpleOrbitFacetInv {
  int incd;
  mpz_class eOrbitSize;
  std::vector<IntInv> ListNb;
};

template<typename T>
bool operator==(SimpleOrbitFacetInv<T> const& x, SimpleOrbitFacetInv<T> const& y)
{
  if (x.incd != y.incd)
    return false;
  if (x.eOrbitSize != y.eOrbitSize)
    return false;
  return x.ListNb == y.ListNb;
}

template<typename T>
std::istream& operator>>(std::istream& is, SimpleOrbitFacetInv<T>& obj)
{
  int incd;
  mpz_class eOrbitSize;
  std::vector<IntInv> ListNb;
  is >> incd;
  is >> eOrbitSize;
  is >> ListNb;
  obj={incd, eOrbitSize, ListNb};
  return is;
}


template<typename T>
std::ostream& operator<<(std::ostream& os, SimpleOrbitFacetInv<T> const& obj)
{
  os << obj.incd << " ";
  os << obj.eOrbitSize;
  os << obj.ListNb;
  return os;
}

template<typename T>
bool operator<(SimpleOrbitFacetInv<T> const& eOrb, SimpleOrbitFacetInv<T> const& fOrb)
{
  // We prefer to treat small incidence first
  if (eOrb.incd < fOrb.incd)
    return true;
  if (eOrb.incd > fOrb.incd)
    return false;
  // At equal incidence, it is best to treat first the large orbits
  if (eOrb.eOrbitSize > fOrb.eOrbitSize)
    return true;
  if (eOrb.eOrbitSize < fOrb.eOrbitSize)
    return false;
  int eLen=eOrb.ListNb.size();
  int fLen=fOrb.ListNb.size();
  if (eLen < fLen)
    return true;
  if (eLen > fLen)
    return false;
  for (int i=0; i<eLen; i++) {
    if (eOrb.ListNb[i] < fOrb.ListNb[i])
      return true;
    if (eOrb.ListNb[i] > fOrb.ListNb[i])
      return false;
  }
  return false;
}

//
// PolyhedralBalinski
//

struct PolyhedralBalinski {
  bool final=false;
  bool IsFirst=true;
  bool IsComplete=true;
  int nbOrbitUnsolved=0;
  mpz_class nbUnsolved=0;
  std::vector<int> rList={};
};

std::ostream& operator<<(std::ostream& os, PolyhedralBalinski const& obj)
{
  os << "PolyhedralBalinski : final =" << obj.final << "\n";
  os << "  IsFirst=" << obj.IsFirst << " IsComplete=" << obj.IsComplete << "\n";
  os << "  nbUnsolved=" << obj.nbUnsolved << "\n";
  os << "  nbOrbitUnvoled=" << obj.nbOrbitUnsolved << "\n";
  int siz=obj.rList.size();
  int eSum=0;
  for (int i=0; i<siz; i++)
    eSum += obj.rList[i];
  os << "  siz=" << siz << " eSum=" << eSum << "\n";
  return os;
}



//
// template type mappings
//

template <typename T,typename Tgroup>
struct invariant_info<PolyhedralEntry<T,Tgroup> > {
  typedef PolyhedralInv<T> invariant_type;
};

template <typename T>
struct invariant_info<SimpleOrbitFacet<T>> {
  typedef SimpleOrbitFacetInv<T> invariant_type;
};

template <typename T,typename Tgroup>
struct equiv_info<PolyhedralEntry<T,Tgroup>> {
  typedef permlib::Permutation equiv_type;
};

template <typename T>
struct equiv_info<SimpleOrbitFacet<T>> {
  typedef permlib::Permutation equiv_type;
};

template <typename T>
struct balinski_info<SimpleOrbitFacet<T>> {
  typedef PolyhedralBalinski balinski_type;
};

//
// the proper functionality now
//

template<typename T,typename Tgroup>
FctsDataBank<PolyhedralEntry<T,Tgroup>> GetRec_FctsDataBank()
{
  using Telt=typename Tgroup::Telt;
  std::function<EquivTest<Telt>(PolyhedralEntry<T,Tgroup> const&,PolyhedralEntry<T,Tgroup> const&)> fEquiv=[](PolyhedralEntry<T,Tgroup> const& eRec1, PolyhedralEntry<T,Tgroup> const& eRec2) -> EquivTest<Telt> {
    MyMatrix<T> EXTred1=ColumnReduction(eRec1.EXT);
    MyMatrix<T> EXTred2=ColumnReduction(eRec2.EXT);
    WeightMatrix<T,T> WMat1=GetWeightMatrix(EXTred1);
    WeightMatrix<T,T> WMat2=GetWeightMatrix(EXTred2);
    EquivTest<Telt> eResEquiv=TestEquivalenceWeightMatrix<T,T,Telt>(WMat1, WMat2);
    if (!eResEquiv.TheReply)
      return {false, {}};
    return {true, eResEquiv.TheEquiv};
  };
  std::function<int(PolyhedralEntry<T,Tgroup> const&)> fSize=[](PolyhedralEntry<T,Tgroup> const& eRec) -> int {
    int siz=eRec.EXT.rows();
    return siz;
  };
  return {fEquiv, fSize};
}


template<typename T,typename Tgroup>
std::vector<Face> DUALDESC_THR_AdjacencyDecomposition(
         MainProcessor &MProc, int const& TheId,
	 DataBank<PolyhedralEntry<T,Tgroup>> &TheBank,
	 MyMatrix<T> const& EXT,
	 Tgroup const& GRP,
	 PolyHeuristic<mpz_class> const& AllArr,
	 std::string const& ePrefix,
	 int const& TheLevel)
{
  using Tint=typename Tgroup::Tint;
  using Telt=typename Tgroup::Telt;
  MyMatrix<T> EXTred=ColumnReduction(EXT);
  int nbRow=EXTred.rows();
  int eRank=EXTred.cols();
  WeightMatrix<T, T> WMat;
  bool HaveWMat=false;
  auto ComputeWMat=[&]() -> void {
    if (HaveWMat)
      return;
    WMat=GetWeightMatrix(EXTred);
    int nbWeight=WMat.GetWeightSize();
    MProc.GetO(TheId) << "nbWeight=" << nbWeight << "\n";
  };
  int TheMinSize=TheBank.GetMinSize();
  if (TheMinSize != -1 && nbRow >= TheMinSize) {
    ComputeWMat();
    T eValInv=GetInvariantWeightMatrix(WMat);
    PolyhedralInv<T> eInv{nbRow, eValInv};
    PolyhedralEntry<T,Tgroup> eEnt{EXT, GRP, {}};
    DataBank_ResultQuery<PolyhedralEntry<T,Tgroup>> eResBank=TheBank.ProcessRequest(eEnt, eInv, MProc.GetO(TheId));
    if (eResBank.test) {
      MProc.GetO(TheId) << "Begin the use of bank data\n";
      std::vector<Face> ListReprTrans;
      for (auto const& eOrbit : eResBank.eEnt.ListFace) {
	Face eListJ=OnFace(eOrbit, eResBank.TheEquiv);
	ListReprTrans.push_back(eListJ);
      }
      MProc.GetO(TheId) << "Before the orbit splitting |ListReprTrans|=" << ListReprTrans.size() << "\n";
      Tgroup GRPconj=ConjugateGroup(eResBank.eEnt.GRP, eResBank.TheEquiv);
      std::vector<Face> ListFaceRet=OrbitSplittingListOrbit(GRPconj, GRP, ListReprTrans, MProc.GetO(TheId));
      MProc.GetO(TheId) << "After the OrbitSplitting\n";
      for (auto & eFace : ListFaceRet) {
	TestFacetness(EXT, eFace);
      }
      MProc.GetO(TheId) << "After the checks\n";
      return ListFaceRet;
    }
  }
  Tgroup TheGRPrelevant;
  std::map<std::string, mpz_class> TheMap;
  int nbVert=EXT.rows();
  int delta=nbVert - eRank;
  TheMap["groupsize"]=GRP.size();
  TheMap["incidence"]=nbVert;
  TheMap["level"]=TheLevel;
  TheMap["rank"]=eRank;
  TheMap["delta"]=delta;
  std::string ansSplit=HeuristicEvaluation(TheMap, AllArr.Splitting);
  MProc.GetO(TheId) << "Not found in bank TheLevel=" << TheLevel << " ansSplit=" << ansSplit << "\n";
  //
  // The computations themselves
  //
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  std::vector<Face> ListOrbitFaces;
  std::string ansSymm;
  if (ansSplit != "split") {
    TheGRPrelevant=GRP;
    std::string ansProg=HeuristicEvaluation(TheMap, AllArr.DualDescriptionProgram);
    ListOrbitFaces=DirectFacetOrbitComputation(EXTred, GRP, ansProg, MProc.GetO(TheId));
    ansSymm="no";
  }
  else {
    ComputeWMat();
    ansSymm=HeuristicEvaluation(TheMap, AllArr.AdditionalSymmetry);
    MProc.GetO(TheId) << "ansSymm=" << ansSymm << "   |EXT|=" << EXT.rows() << "\n";
    if (ansSymm == "yes")
      TheGRPrelevant=GetStabilizerWeightMatrix(WMat);
    else
      TheGRPrelevant=GRP;
    if (TheGRPrelevant.size() == GRP.size())
      ansSymm="no";
    TheMap["groupsizerelevant"]=TheGRPrelevant.size();
    std::string ansGRP=HeuristicEvaluation(TheMap, AllArr.StabEquivFacet);
    std::string ansStratLocInv=HeuristicEvaluation(TheMap, AllArr.InvariantQuality);
    LocalInvInfo eLocalInv=ComputeLocalInvariantStrategy(WMat, TheGRPrelevant, ansStratLocInv, MProc.GetO(TheId));
    MProc.GetO(TheId) << "nbDiagCoeff=" << eLocalInv.nbDiagCoeff << " nbOffCoeff=" << eLocalInv.nbOffCoeff << "\n";
    Tint QuotSize=TheGRPrelevant.size() / GRP.size();
    MProc.GetO(TheId) << "ansSymm=" << ansSymm << " ansGRP=" << ansGRP << " |TheGRPrelevant|=" << TheGRPrelevant.size() << " |GRP|=" << GRP.size() << " QuotSize=" << QuotSize << "\n";
    mpz_class MaxAllowedUndone=eRank-2;
    std::function<bool(SimpleOrbitFacetInv<T> const&,SimpleOrbitFacetInv<T> const&)> CompFCT=[](SimpleOrbitFacetInv<T> const& x, SimpleOrbitFacetInv<T> const& y) -> bool {
      return x < y;
    };
    std::function<void(PolyhedralBalinski &,SimpleOrbitFacet<T> const&,SimpleOrbitFacetInv<T> const&,std::ostream&)> UpgradeBalinskiStat=[&](PolyhedralBalinski & eStat, SimpleOrbitFacet<T> const& fEnt, SimpleOrbitFacetInv<T> const& fInv, std::ostream& os) -> void {
      if (eStat.final)
	return;
      eStat.nbUnsolved += fInv.eOrbitSize;
      eStat.nbOrbitUnsolved++;
      //
      std::vector<int> gList=FaceTo01vector(fEnt.eRepr);
      std::vector<int> rList = OrbitIntersection(GRP, gList);
      if (eStat.IsFirst) {
	eStat.rList=rList;
	eStat.IsFirst=false;
      }
      else {
	for (int iVert=0; iVert<nbVert; iVert++)
	  eStat.rList[iVert] *= rList[iVert];
      }
      int eSum=0;
      for (int iVert=0; iVert<nbVert; iVert++)
	eSum += eStat.rList[iVert];
      if (eSum == 0 && eStat.nbUnsolved > MaxAllowedUndone) {
	eStat.final=true;
	eStat.IsComplete=false;
      }
    };
    int NewLevel=TheLevel+1;
    std::function<EquivTest<Telt>(SimpleOrbitFacet<T> const&,SimpleOrbitFacet<T> const&)> fEquiv;
    std::function<PairT_Tinv<SimpleOrbitFacet<T>>(Face const&, std::ostream&)> GetRecord;
    {
      std::ofstream os1("DEBUG_GRPrelevant");
      std::ofstream os2("DEBUG_GRPrelevant.gap");
      WriteGroup(os1, TheGRPrelevant);
      WriteGroupGAP(os2, TheGRPrelevant);
    }
    if (ansGRP == "classic") {
      fEquiv=[&](SimpleOrbitFacet<T> const& x, SimpleOrbitFacet<T> const& y) -> EquivTest<Telt> {
	std::chrono::time_point<std::chrono::system_clock> startLoc, endLoc;
	startLoc = std::chrono::system_clock::now();
	auto eReply=TheGRPrelevant.RepresentativeAction_OnSets(x.eRepr, y.eRepr);
	endLoc = std::chrono::system_clock::now();
	int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(endLoc-startLoc).count();
	MProc.GetO(TheId) << "CLASSIC: After the test time = " << elapsed_seconds << "\n";
	return eReply;
      };
      GetRecord=[&](Face const& eOrb, std::ostream &os) -> PairT_Tinv<SimpleOrbitFacet<T>> {
        Tgroup TheStab=TheGRPrelevant.Stabilizer_OnSets(eOrb);
	int siz=eOrb.count();
	Tint eOrbitSize=TheGRPrelevant.size() / TheStab.size();
	SimpleOrbitFacet<T> eOrbF{eOrb};
	std::vector<IntInv> ListSingleInv=GetLocalInvariantWeightMatrix_Enhanced<IntInv>(eLocalInv,eOrb);
	SimpleOrbitFacetInv<T> eInv{siz, eOrbitSize, ListSingleInv};
	return {eOrbF, eInv};
      };
    }
    if (ansGRP == "partition") {
      fEquiv=[&TheGRPrelevant,&MProc,&TheId,&WMat](SimpleOrbitFacet<T> const& x, SimpleOrbitFacet<T> const& y) -> EquivTest<Telt> {
	std::chrono::time_point<std::chrono::system_clock> startloc, endloc;
	startloc = std::chrono::system_clock::now();
	{
	  std::ofstream os1("DEBUG_PairFace");
	  std::ofstream os2("DEBUG_PairFace.gap");
	  WriteListFace(os1, {x.eRepr, y.eRepr});
	  WriteListFaceGAP(os2, {x.eRepr, y.eRepr});
	}
	//	MProc.GetO(TheId) << "TheGRPrelevant.n=" << TheGRPrelevant.n << "\n";
	//	MProc.GetO(TheId) << "|x.eRepr|=" << x.eRepr.size() << " |y.eRepr|=" << y.eRepr.size() << "\n";
	auto eReply=TestEquivalenceGeneralNaked(TheGRPrelevant, x.eRepr, y.eRepr, 1);
	endloc = std::chrono::system_clock::now();
	int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(endloc-startloc).count();
	MProc.GetO(TheId) << "PARTITION: After the test time = " << elapsed_seconds << "\n";
	if (elapsed_seconds > 60) {
	  //
	  std::chrono::time_point<std::chrono::system_clock> start_C, end_C;
	  start_C = std::chrono::system_clock::now();
	  auto eReplyB=TestEquivalenceSubset(WMat, x.eRepr, y.eRepr);
	  end_C = std::chrono::system_clock::now();
	  int elapsed_seconds_C = std::chrono::duration_cast<std::chrono::seconds>(end_C-start_C).count();
	  MProc.GetO(TheId) << "Second method (bliss) runtime = " << elapsed_seconds_C << "\n";
	  if (eReply.TheReply != eReplyB.TheReply) {
	    std::cerr << "We have a major bug to solve\n";
	    throw TerminalException{1};
	  }
	}
	return eReply;
      };
      GetRecord=[&](Face const& eOrb, std::ostream &os) -> PairT_Tinv<SimpleOrbitFacet<T>> {
        Tgroup TheStab=TheGRPrelevant.Stabilizer_OnSets(eOrb);
	int siz=eOrb.count();
	Tint eOrbitSize=TheGRPrelevant.size() / TheStab.size();
	SimpleOrbitFacet<T> eOrbF{eOrb};
	std::vector<IntInv> ListSingleInv=GetLocalInvariantWeightMatrix_Enhanced<IntInv>(eLocalInv,eOrb);
	SimpleOrbitFacetInv<T> eInv{siz, eOrbitSize, ListSingleInv};
	return {eOrbF, eInv};
      };
    }
    if (ansGRP == "exhaustive") {
      // we choose here to discard the element realizing the equivalence
      fEquiv=[&](SimpleOrbitFacet<T> const& x, SimpleOrbitFacet<T> const& y) -> EquivTest<Telt> {
	bool test=x.eRepr==y.eRepr;
	return {test, {}};
      };
      OrbitMinimumArr<Tint> ArrMin=GetInitialMinimumArray(TheGRPrelevant);
      GetRecord=[ArrMin](Face const& eOrb, std::ostream &os) -> PairT_Tinv<SimpleOrbitFacet<T>> {
        ResultMinimum<Tint> ResMin=GetMinimumRecord(ArrMin, eOrb);
	int siz=eOrb.count();
	SimpleOrbitFacet<T> eOrbF{ResMin.eMin};
	SimpleOrbitFacetInv<T> eInv{siz, ResMin.OrbitSize, {}};
	return {eOrbF, eInv};
      };
    }
    bool Saving=AllArr.Saving;
    bool eMemory=AllArr.eMemory;
    NewEnumerationWork<SimpleOrbitFacet<T>> ListOrbit(Saving, eMemory, ePrefix, CompFCT, UpgradeBalinskiStat, fEquiv, MProc.GetO(TheId));
    mpz_class TotalNumberFacet=0;
    auto FuncInsert=[&](Face const& eOrb, std::ostream &os) -> int {
      PairT_Tinv<SimpleOrbitFacet<T>> eRec=GetRecord(eOrb,os);
      int eVal=ListOrbit.InsertEntry(eRec, os);
      if (eVal == -1)
	TotalNumberFacet += eRec.xInv.eOrbitSize;
      return eVal;
    };
    int nbPresentOrbit=ListOrbit.GetNbEntry();
    if (nbPresentOrbit == 0) {
      std::string ansSamp=HeuristicEvaluation(TheMap, AllArr.InitialFacetSet);
      MProc.GetO(TheId) << "Before InitialFacetComputation ansSamp=" << ansSamp << "\n";
      std::vector<Face> ListFace=DirectComputationInitialFacetSet(EXTred, ansSamp);
      MProc.GetO(TheId) << " After InitialFacetComputation\n";
      for (auto & eInc : ListFace) {
	int RetVal=FuncInsert(eInc, MProc.GetO(TheId));
	MProc.GetO(TheId) << "Inserting |eInc|=" << eInc.count() <<" : RetVal=" << RetVal << "\n";
      }
    }
    std::atomic<int> nbSpannThread(0);
    std::atomic<int> nbStuckThread(0);
    std::condition_variable cv;
    std::mutex mtx_cv;
    auto WaitStuck=[&](int const& opt, int const& MyId) -> void {
      MProc.GetO(MyId) << "WaitStuck, nbSpannThread=" << nbSpannThread << " stuck=" << ListOrbit.IsStuck() << "\n";
      MProc.decNRT(MyId);
      nbStuckThread++;
      std::unique_lock<std::mutex> lk(mtx_cv);
      cv.wait(lk, [&]{if (opt == 0)
	    return ListOrbit.IsStuck() == false;
	  else
	    return nbSpannThread == 0;});
      MProc.incNRT(MyId);
      nbStuckThread--;
      MProc.GetO(MyId) << "Exiting WaitStuck\n";
    };
    std::function<void(void)> SpannNewThread;
    auto TreatDatabase=[&](int const& MyId) -> void {
      int nbWork=0;
      std::ostream & os = MProc.GetO(MyId);
      nbSpannThread++;
      while(true) {
	bool testStuck=ListOrbit.IsStuck();
	os << "nbEntry=" << ListOrbit.GetNbEntry() << " testStuck=" << testStuck << "\n";
	if (testStuck)
	  WaitStuck(0, MyId);
	bool IsComplete=ListOrbit.GetCompleteStatus();
	if (IsComplete)
	  break;
	int eEntry=ListOrbit.GetNonTreatedOrbit(os);
	if (eEntry == -1)
	  break;
	Face eListI=ListOrbit.GetRepresentative(eEntry).eRepr;
	MyMatrix<T> EXTredFace=SelectRow(EXT, eListI);
	Tgroup TheStab=TheGRPrelevant.Stabilizer_OnSets(eListI);
	Tint OrbSize=TheGRPrelevant.size() / TheStab.size();
	os << "eEntry=" << eEntry << " |TheGRPrelevant|=" << TheGRPrelevant.size() << "  |TheStab|=" << TheStab.size() << " |O|=" << OrbSize << "\n";
	Tgroup GRPred=ReducedGroupAction(TheStab, eListI);
	CondTempDirectory eDir(AllArr.Saving, ePrefix + "ADM" + IntToString(eEntry) + "/");
	std::vector<Face> TheOutput=DUALDESC_THR_AdjacencyDecomposition(MProc, MyId, TheBank, EXTredFace, GRPred, AllArr, eDir.str(), NewLevel);
	os << "TreatDatabase, NewLevel=" << NewLevel << "  |EXT|=" << EXTredFace.rows() << "  eRank=" << eRank << "  |TheOutput|=" << TheOutput.size() << " |GRPred|=" << GRPred.size() << "\n";
	int iter=0;
	for (auto& eOrbB : TheOutput) {
	  Face eFlip=ComputeFlipping(EXTred, eListI, eOrbB);
	  os << " iter=" << iter << " |eFlip|=" << eFlip.count() << "\n";
	  int eVal=FuncInsert(eFlip, MProc.GetO(MyId));
	  os << " After FuncInsert\n";
	  if (eVal == -1) {
	    if (nbStuckThread > 0) {
	      cv.notify_one();
	    }
	    else {
	      if (MProc.MPU_NumberFree() > 0)
		SpannNewThread();
	    }
	  }
	  iter++;
	}
	os << "\n";
	ListOrbit.SetEntryAsDone(eEntry, os);
	os << "TreatDatabase: After SetEntryAsDone\n";
	nbWork++;
      }
      os << "TreatDatabase: nbWork=" << nbWork << "\n";
      if (MyId != TheId)
	MProc.MPU_Terminate(MyId);
      nbSpannThread--;
      if (nbWork > 0)
	cv.notify_all();
    };
    std::vector<std::thread> ListThreads;
    SpannNewThread=[&]() -> void {
      int NewId=MProc.MPU_GetId();
      MProc.GetO(TheId) << "SpannNewThread NewId=" << NewId << "\n";
      if (NewId != -1) {
	ListThreads.push_back(std::thread(TreatDatabase, NewId));
	ListThreads[ListThreads.size()-1].detach();
      }
    };
    int NbThr=MProc.MPU_NumberFree();
    MProc.GetO(TheId) << "We have NbThr=" << NbThr << "\n";
    while(true) {
      bool IsCompleteSpann=ListOrbit.GetCompleteStatus();
      MProc.GetO(TheId) << "IsCompleteSpann=" << IsCompleteSpann << "\n";
      if (IsCompleteSpann)
	WaitStuck(1, TheId);
      else {
	for (int iThr=0; iThr<NbThr; iThr++)
	  SpannNewThread();
	MProc.GetO(TheId) << "Before TreatDatabase\n";
	TreatDatabase(TheId);
	MProc.GetO(TheId) << "After my own call to TreatDatabase\n";
      }
      if (nbSpannThread == 0)
	break;
    }
    bool IsComplete=ListOrbit.GetCompleteStatus();
    int nbOrbitFacet=ListOrbit.GetNbEntry();
    MProc.GetO(TheId) << "IsComplete=" << IsComplete << "\n";
    MProc.GetO(TheId) << "TotalNumberFacet=" << TotalNumberFacet << "  nbOrbitFacet=" << nbOrbitFacet << "\n";
    if (!IsComplete) {
      std::cerr << "Major error in the code. We should be complete now\n";
      throw TerminalException{1};
    }
    ListOrbitFaces.reserve(nbOrbitFacet);
    for (int iOF=0; iOF<nbOrbitFacet; iOF++) {
      SimpleOrbitFacet<T> x=ListOrbit.GetRepresentative(iOF);
      ListOrbitFaces.push_back(x.eRepr);
    }
    ListOrbit.FuncClear();
  }
  end = std::chrono::system_clock::now();
  int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
  TheMap["time"]=elapsed_seconds;
  std::string ansBank=HeuristicEvaluation(TheMap, AllArr.BankSave);
  MProc.GetO(TheId) << "elapsed_seconds = " << elapsed_seconds << " ansBank = " << ansBank << "\n";
  if (ansBank == "yes") {
    MProc.GetO(TheId) << "BANK work, step 1\n";
    PolyhedralEntry<T,Tgroup> eEntry{EXT, TheGRPrelevant, ListOrbitFaces};
    MProc.GetO(TheId) << "BANK work, step 2\n";
    PolyhedralEntry<T,Tgroup> eEntryCan=CanonicalizationPolyEntry(eEntry, MProc.GetO(TheId));
    MProc.GetO(TheId) << "BANK work, step 3\n";
    T eValInv=GetInvariantWeightMatrix(WMat);
    MProc.GetO(TheId) << "BANK work, step 4\n";
    PolyhedralInv<T> eInv{nbRow, eValInv};
    MProc.GetO(TheId) << "BANK work, step 5\n";
    TheBank.InsertEntry(eEntryCan, eInv);
    MProc.GetO(TheId) << "BANK work, step 6\n";
  }
  MProc.GetO(TheId) << "Bank entry processed\n";
  std::vector<Face> ListOrbitReturn;
  if (ansSymm == "yes") {
    MProc.GetO(TheId) << "|TheGRPrelevant|=" << TheGRPrelevant.size() << " |GRP|=" << GRP.size() << "\n";
    {
      std::ofstream os1("SPLIT_GRPrelevant");
      std::ofstream os2("SPLIT_GRP");
      std::ofstream os3("SPLIT_ListOrbitFaces");
      WriteGroup(os1, TheGRPrelevant);
      WriteGroup(os2, GRP);
      WriteListFace(os3, ListOrbitFaces);
    }
    ListOrbitReturn=OrbitSplittingListOrbit(TheGRPrelevant, GRP, ListOrbitFaces, MProc.GetO(TheId));
  }
  else
    ListOrbitReturn=ListOrbitFaces;
  MProc.GetO(TheId) << "Return list has been computed\n";
  return ListOrbitReturn;
}


FullNamelist NAMELIST_GetStandard_TEMP_THREADED_ADM()
{
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string> > ListListStringValues1;
  ListStringValues1["EXTfile"]="unset.ext";
  ListStringValues1["GRPfile"]="unset.grp";
  ListStringValues1["OUTfile"]="unset.out";
  SingleBlock BlockDATA;
  BlockDATA.ListIntValues=ListIntValues1;
  BlockDATA.ListBoolValues=ListBoolValues1;
  BlockDATA.ListDoubleValues=ListDoubleValues1;
  BlockDATA.ListStringValues=ListStringValues1;
  BlockDATA.ListListStringValues=ListListStringValues1;
  ListBlock["DATA"]=BlockDATA;
  // HEURISTIC
  std::map<std::string, std::string> ListStringValuesH;
  ListStringValuesH["SplittingHeuristicFile"]="unset.heu";
  ListStringValuesH["AdditionalSymmetryHeuristicFile"]="unset.heu";
  ListStringValuesH["DualDescriptionHeuristicFile"]="unset.heu";
  ListStringValuesH["StabEquivFacetHeuristicFile"]="unset.heu";
  ListStringValuesH["MethodInitialFacetSetFile"]="unset.heu";
  ListStringValuesH["BankSaveHeuristicFile"]="unset.heu";
  ListStringValuesH["MethodInvariantQualityFile"]="unset.heu";
  SingleBlock BlockHEURIS;
  BlockHEURIS.ListStringValues=ListStringValuesH;
  ListBlock["HEURISTIC"]=BlockHEURIS;
  // METHOD
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string> > ListListStringValues2;
  ListIntValues2["NPROC"]=1;
  ListBoolValues2["Saving"]=false;
  ListStringValues2["Prefix"]="/irrelevant/";
  ListBoolValues2["FullDataInMemory"]=true;
  SingleBlock BlockMETHOD;
  BlockMETHOD.ListIntValues=ListIntValues2;
  BlockMETHOD.ListBoolValues=ListBoolValues2;
  BlockMETHOD.ListDoubleValues=ListDoubleValues2;
  BlockMETHOD.ListStringValues=ListStringValues2;
  BlockMETHOD.ListListStringValues=ListListStringValues2;
  ListBlock["METHOD"]=BlockMETHOD;
  // BANK
  std::map<std::string, int> ListIntValues3;
  std::map<std::string, bool> ListBoolValues3;
  std::map<std::string, double> ListDoubleValues3;
  std::map<std::string, std::string> ListStringValues3;
  std::map<std::string, std::vector<std::string> > ListListStringValues3;
  ListStringValues3["Prefix"]="./unset/";
  ListBoolValues3["Saving"]=false;
  ListBoolValues3["FullDataInMemory"]=true;
  SingleBlock BlockBANK;
  BlockBANK.ListIntValues=ListIntValues3;
  BlockBANK.ListBoolValues=ListBoolValues3;
  BlockBANK.ListDoubleValues=ListDoubleValues3;
  BlockBANK.ListStringValues=ListStringValues3;
  BlockBANK.ListListStringValues=ListListStringValues3;
  ListBlock["BANK"]=BlockBANK;
  // Merging all data
  return {ListBlock, "undefined"};
}

template<typename T, typename Tgroup>
void MainFunctionComputeDualDesc(FullNamelist const& eFull)
{
  SingleBlock BlockBANK=eFull.ListBlock.at("BANK");
  bool BANK_IsSaving=BlockBANK.ListBoolValues.at("Saving");
  bool BANK_Memory=BlockBANK.ListBoolValues.at("FullDataInMemory");
  std::string BANK_Prefix=BlockBANK.ListStringValues.at("Prefix");
  FctsDataBank<PolyhedralEntry<T,Tgroup>> recFct=GetRec_FctsDataBank<T,Tgroup>();
  DataBank<PolyhedralEntry<T,Tgroup>> TheBank(BANK_IsSaving, BANK_Memory, BANK_Prefix, recFct);
  //
  std::cerr << "Reading DATA\n";
  SingleBlock BlockDATA=eFull.ListBlock.at("DATA");
  std::string EXTfile=BlockDATA.ListStringValues.at("EXTfile");
  IsExistingFileDie(EXTfile);
  std::cerr << "EXTfile=" << EXTfile << "\n";
  std::string GRPfile=BlockDATA.ListStringValues.at("GRPfile");
  IsExistingFileDie(GRPfile);
  std::cerr << "GRPfile=" << GRPfile << "\n";
  std::string OUTfile=BlockDATA.ListStringValues.at("OUTfile");
  std::cerr << "OUTfile=" << OUTfile << "\n";
  std::ifstream EXTfs(EXTfile);
  MyMatrix<T> EXT=ReadMatrix<T>(EXTfs);
  std::ifstream GRPfs(GRPfile);
  Tgroup GRP=ReadGroup<Tgroup>(GRPfs);
  //
  std::cerr << "Creating MPROC\n";
  SingleBlock BlockMETHOD=eFull.ListBlock.at("METHOD");
  int NbThr=BlockMETHOD.ListIntValues.at("NPROC");
  MainProcessor MProc(NbThr);
  int TheId=MProc.MPU_GetId();
  //
  PolyHeuristic<mpz_class> AllArr=AllStandardHeuristic<mpz_class>();
  //
  auto SetHeuristic=[&](std::string const& NamelistEnt, TheHeuristic<mpz_class> & eHeu) {
    SingleBlock BlockHEU=eFull.ListBlock.at("HEURISTIC");
    std::string NamelistEntFile=BlockHEU.ListStringValues.at(NamelistEnt);
    if (NamelistEntFile != "unset.heu") {
      IsExistingFileDie(NamelistEntFile);
      std::ifstream is(NamelistEntFile);
      eHeu=ReadHeuristic<mpz_class>(is);
    }
  };
  SetHeuristic("SplittingHeuristicFile", AllArr.Splitting);
  SetHeuristic("AdditionalSymmetryHeuristicFile", AllArr.AdditionalSymmetry);
  SetHeuristic("DualDescriptionHeuristicFile", AllArr.DualDescriptionProgram);
  SetHeuristic("StabEquivFacetHeuristicFile", AllArr.StabEquivFacet);
  SetHeuristic("MethodInitialFacetSetFile", AllArr.InitialFacetSet);
  SetHeuristic("MethodInvariantQualityFile", AllArr.InvariantQuality);
  SetHeuristic("BankSaveHeuristicFile", AllArr.BankSave);
  //
  bool DD_Saving=BlockMETHOD.ListBoolValues.at("Saving");
  bool DD_Memory=BlockMETHOD.ListBoolValues.at("FullDataInMemory");
  std::string DD_Prefix=BlockMETHOD.ListStringValues.at("Prefix");
  AllArr.Saving=DD_Saving;
  AllArr.eMemory=DD_Memory;
  //
  int TheLevel=0;
  std::vector<Face> TheOutput=DUALDESC_THR_AdjacencyDecomposition(MProc,
								  TheId,
								  TheBank, EXT,
								  GRP,
								  AllArr,
								  DD_Prefix, TheLevel);
  std::cerr << "We now have TheOutput\n";
  //
  std::ofstream OUTfs(OUTfile);
  VectVectInt_Magma_Print(OUTfs, TheOutput);
  //
}

#endif

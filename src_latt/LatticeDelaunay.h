#ifndef LATTICE_DELAUNAY_INCLUDE
#define LATTICE_DELAUNAY_INCLUDE

#include "FundamentalDelaunay.h"
#include "GRP_GroupFct.h"
#include "MatrixGroup.h"
#include "POLY_ThreadDualDescription.h"
#include "Namelist.h"

template<typename T, typename Tint>
struct DataLattice {
  int n;
  MyMatrix<T> GramMat;
  MyMatrix<Tint> SHV;
  bool SavingDelaunay;
  bool FullDataInMemory;
  std::string PrefixDelaunay;
  std::string PrefixPolyhedral;
  bool ReturnAll;
  mpz_class UpperLimitMethod4;
  std::string CVPmethod;
};



// The types for the Delaunay polytopes
template<typename T,typename Tint>
struct Delaunay {
  MyMatrix<Tint> EXT;
};

template<typename T,typename Tint>
std::istream& operator>>(std::istream& is, Delaunay<T,Tint>& obj)
{
  MyMatrix<Tint> EXT=ReadMatrix<Tint>(is);
  obj={EXT};
  return is;
}


template<typename T,typename Tint>
std::ostream& operator<<(std::ostream& os, Delaunay<T,Tint> const& obj)
{
  WriteMatrix(os, obj.EXT);
  return os;
}

template<typename T,typename Tint>
struct DelaunayInv {
  int nbVert;
  Tint eIndex;
  T ePolyInv;
};

template<typename T,typename Tint>
  bool operator==(DelaunayInv<T,Tint> const& x, DelaunayInv<T,Tint> const& y)
{
  if (x.nbVert != y.nbVert)
    return false;
  if (x.eIndex != y.eIndex)
    return false;
  if (x.ePolyInv != y.ePolyInv)
    return false;
  return true;
}



template<typename T,typename Tint>
std::istream& operator>>(std::istream& is, DelaunayInv<T,Tint>& obj)
{
  int nbVert;
  Tint eIndex;
  T ePolyInv;
  is >> nbVert;
  is >> eIndex;
  is >> ePolyInv;
  obj={nbVert, eIndex, ePolyInv};
  return is;
}


template<typename T,typename Tint>
std::ostream& operator<<(std::ostream& os, DelaunayInv<T,Tint> const& obj)
{
  os << " " << obj.nbVert;
  os << " " << obj.eIndex;
  os << " " << obj.ePolyInv;
  return os;
}


template<typename T,typename Tint>
bool operator<(DelaunayInv<T,Tint> const& x, DelaunayInv<T,Tint> const& y)
{
  if (x.nbVert < y.nbVert)
    return true;
  if (x.nbVert > y.nbVert)
    return false;
  if (x.eIndex < y.eIndex)
    return true;
  if (x.eIndex > y.eIndex)
    return false;
  return x.ePolyInv < y.ePolyInv;
}


template<typename T,typename Tint>
struct invariant_info<Delaunay<T,Tint>> {
  typedef DelaunayInv<T,Tint> invariant_type;
};

template<typename T,typename Tint>
  struct equiv_info<Delaunay<T,Tint>> {
  typedef MyMatrix<Tint> equiv_type;
};










template<typename T>
WeightMatrix<T, T> GetWeightMatrixFromGramEXT(MyMatrix<T> const& EXT, MyMatrix<T> const& GramMat, MyMatrix<T> const& SHV)
{
  int n=GramMat.rows();
  CP<T> eCP=CenterRadiusDelaunayPolytopeGeneral(GramMat, EXT);
  MyVector<T> TheCenter=eCP.eCent;
  std::vector<MyMatrix<T> > CharPair=CharacterizingPair(GramMat, TheCenter);
  MyMatrix<T> Qmat=CharPair[0];
  int nbVect=SHV.rows();
  int nbVert=EXT.rows();
  MyMatrix<T> EXText(nbVect + nbVert, n+1);
  for (int iVert=0; iVert<nbVert; iVert++)
    EXText.row(iVert)=EXT.row(iVert);
  for (int iVect=0; iVect<nbVect; iVect++) {
    EXText(nbVert+iVect,0)=0;
    for (int i=0; i<n; i++)
      EXText(nbVert+iVect,i+1)=SHV(iVect,i);
  }
  return GetSimpleWeightMatrix(EXText, Qmat);
}

template<typename T, typename Tgroup>
bool IsGroupCorrect(MyMatrix<T> const& EXT_T, Tgroup const& eGRP)
{
  using Telt=typename Tgroup::Telt;
  std::vector<Telt> LGen = eGRP.GeneratorsOfGroup();
  for (auto & eGen : LGen) {
    MyMatrix<T> eMat=FindTransformation(EXT_T, EXT_T, eGen);
    if (!IsIntegralMatrix(eMat))
      return false;
  }
  return true;
}




template<typename T, typename Tint, typename Tgroup>
Tgroup Delaunay_Stabilizer(DataLattice<T,Tint> const& eData, MyMatrix<Tint> const & EXT, Tint const& eIndex)
{
  MyMatrix<T> EXT_T=ConvertMatrixUniversal<T,Tint>(EXT);
  WeightMatrix<T,T> WMatRed=GetWeightMatrixFromGramEXT(EXT_T, eData.GramMat, {});
  Tgroup GRPisomRed=GetStabilizerWeightMatrix<T,T,Tgroup>(WMatRed);
  if (IsGroupCorrect(EXT_T, GRPisomRed))
    return GRPisomRed;
  if (eIndex == 1) {
    std::cerr << "Inconsistent. If the index is 1, then GRPisomRed\n";
    std::cerr << "should be the full group\n";
    throw TerminalException{1};
  }
  //
  // Now extending with the SHV vector set
  //
  MyMatrix<T> SHV_T=ConvertMatrixUniversal<T,Tint>(eData.SHV);
  WeightMatrix<T, T> WMat=GetWeightMatrixFromGramEXT(EXT_T, eData.GramMat, SHV_T);
  int nbVert=EXT_T.rows();
  int nbSHV=SHV_T.rows();
  Face eFace(nbVert + nbSHV);
  for (int iVert=0; iVert<nbVert; iVert++)
    eFace[iVert]=1;
  for (int iSHV=0; iSHV<nbSHV; iSHV++)
    eFace[nbVert+iSHV]=0;
  Tgroup PreGRPisom=GetStabilizerWeightMatrix<T,T,Tgroup>(WMat);
  Tgroup GRPisom=ReducedGroupAction(PreGRPisom, eFace);
  if (IsGroupCorrect(EXT_T, GRPisom)) {
    return GRPisom;
  }
  std::function<bool(MyMatrix<T>)> IsMatrixCorrect=[](MyMatrix<T> const& M) -> bool {
    return IsIntegralMatrix(M);
  };
  if (GRPisom.size() < eData.UpperLimitMethod4) {
    return LinPolytopeIntegral_Stabilizer_Method4(EXT_T, GRPisom, IsMatrixCorrect);
  }
  else {
    return LinPolytopeIntegral_Stabilizer_Method8(EXT_T, GRPisom);
  }
}



template<typename T, typename Tint, typename Tgroup>
EquivTest<MyMatrix<Tint>> Delaunay_TestEquivalence(DataLattice<T,Tint> const& eData, Delaunay<T,Tint> const& RecEXT1, Delaunay<T,Tint> const& RecEXT2, Tint const& eIndex)
{
  using Telt=typename Tgroup::Telt;
  std::cerr << "Begin Delaunay_TestEquivalence\n";
  auto ConvertEquiv=[](EquivTest<MyMatrix<T>> const& eEq) -> EquivTest<MyMatrix<Tint>> {
    if (!eEq.TheReply)
      return {false, {}};
    MyMatrix<Tint> eMat_I=ConvertMatrixUniversal<Tint,T>(eEq.TheEquiv);
    return {true, eMat_I};
  };
  MyMatrix<T> EXT1_T=ConvertMatrixUniversal<T,Tint>(RecEXT1.EXT);
  MyMatrix<T> EXT2_T=ConvertMatrixUniversal<T,Tint>(RecEXT2.EXT);
  WeightMatrix<T, T> WMatRed1=GetWeightMatrixFromGramEXT(EXT1_T, eData.GramMat, {});
  WeightMatrix<T, T> WMatRed2=GetWeightMatrixFromGramEXT(EXT2_T, eData.GramMat, {});
  EquivTest<Telt> eResRed=TestEquivalenceWeightMatrix<T,T,Telt>(WMatRed1, WMatRed2);
  if (!eResRed.TheReply) {
    std::cerr << "Leaving Delaunay_TestEquivalence with false\n";
    return {false, {}};
  }
  MyMatrix<T> MatEquivRed_T=FindTransformation(EXT1_T, EXT2_T, eResRed.TheEquiv);
  if (IsIntegralMatrix(MatEquivRed_T)) {
    MyMatrix<Tint> MatEquiv_I=ConvertMatrixUniversal<Tint,T>(MatEquivRed_T);
    std::cerr << "Leaving Delaunay_TestEquivalence with true\n";
    return {true, MatEquiv_I};
  };
  if (eIndex == 1) {
    std::cerr << "We should not reach that stage. If eIndex=1 then\n";
    std::cerr << "the basic algorithm should work\n";
    throw TerminalException{1};
  }
  //
  // Now extending by adding more vectors.
  //
  MyMatrix<T> SHV_T=ConvertMatrixUniversal<T,Tint>(eData.SHV);
  WeightMatrix<T, T> WMat1=GetWeightMatrixFromGramEXT(EXT1_T, eData.GramMat, SHV_T);
  WeightMatrix<T, T> WMat2=GetWeightMatrixFromGramEXT(EXT2_T, eData.GramMat, SHV_T);
  EquivTest<Telt> eRes=TestEquivalenceWeightMatrix<T,T,Telt>(WMat1, WMat2);
  if (!eRes.TheReply) {
    std::cerr << "Leaving Delaunay_TestEquivalence with false\n";
    return {false, {}};
  }
  MyMatrix<T> MatEquiv_T=FindTransformation(EXT1_T, EXT2_T, eRes.TheEquiv);
  if (IsIntegralMatrix(MatEquiv_T)) {
    MyMatrix<Tint> MatEquiv_I=ConvertMatrixUniversal<Tint,T>(MatEquiv_T);
    std::cerr << "Leaving Delaunay_TestEquivalence with true\n";
    return {true, MatEquiv_I};
  };
  std::cerr << "Trying other strategies\n";
  Tgroup GRP1=GetStabilizerWeightMatrix<T,T,Tgroup>(WMat1);
  if (GRP1.size() < eData.UpperLimitMethod4) {
    std::function<bool(MyMatrix<T>)> IsMatrixCorrect=[](MyMatrix<T> const& M) -> bool {
      return IsIntegralMatrix(M);
    };
    return ConvertEquiv(LinPolytopeIntegral_Isomorphism_Method4(EXT1_T, EXT2_T, GRP1, eRes.TheEquiv, IsMatrixCorrect));
  }
  return ConvertEquiv(LinPolytopeIntegral_Isomorphism_Method8(EXT1_T, EXT2_T, GRP1, eRes.TheEquiv));
}


template<typename T,typename Tint>
DelaunayInv<T,Tint> ComputeInvariantDelaunay(MyMatrix<T> const& GramMat, Delaunay<T,Tint> const& eDel)
{
  int nbVert=eDel.EXT.rows();
  int n=eDel.EXT.cols() - 1;
  Tint PreIndex=Int_IndexLattice(eDel.EXT);
  Tint eIndex=T_abs(PreIndex);
  WeightMatrix<T,T> WMat(nbVert,0);
  for (int i=0; i<nbVert; i++)
    for (int j=0; j<nbVert; j++) {
      MyVector<T> eDiff(n);
      for (int iCol=0; iCol<n; iCol++)
	eDiff(iCol)=eDel.EXT(i,iCol+1) - eDel.EXT(j,iCol+1);
      T eNorm=EvaluationQuadForm<T,T>(GramMat, eDiff);
      WMat.Update(i,j,eNorm);
    }
  T ePolyInv_T=GetInvariantWeightMatrix(WMat);
  DelaunayInv<T,Tint> eInv{nbVert, eIndex, ePolyInv_T};
  return eInv;
}










template<typename T,typename Tint, typename Tgroup>
std::vector<Delaunay<T,Tint>> EnumerationDelaunayPolytopes(
	    MainProcessor &MProc, int const& TheId,
	    DataBank<PolyhedralEntry<T,Tgroup>> &TheBank,
	    DataLattice<T,Tint> const& eData, 
	    PolyHeuristic<mpz_class> const& AllArr)
{
  std::function<bool(DelaunayInv<T,Tint> const&,DelaunayInv<T,Tint> const&)> CompFCT=[](DelaunayInv<T,Tint> const& x, DelaunayInv<T,Tint> const& y) -> bool {
    return x < y;
  };
  std::function<void(TrivialBalinski &,Delaunay<T,Tint> const&,DelaunayInv<T,Tint> const&,std::ostream&)> UpgradeBalinskiStat=[](TrivialBalinski const& eStat, Delaunay<T,Tint> const& eEnt, DelaunayInv<T,Tint> const& eInv, std::ostream&os) -> void {
  };
  std::function<EquivTest<MyMatrix<Tint>>(Delaunay<T,Tint> const&,Delaunay<T,Tint> const&)> fEquiv=[&](Delaunay<T,Tint> const& x, Delaunay<T,Tint> const& y) -> EquivTest<MyMatrix<Tint>> {
    Tint PreIndex=Int_IndexLattice(x.EXT);
    Tint eIndex=T_abs(PreIndex);
    return Delaunay_TestEquivalence<T,Tint,Tgroup>(eData, x, y, eIndex);
  };
  NewEnumerationWork<Delaunay<T,Tint>> ListOrbit(AllArr.Saving, AllArr.eMemory, eData.PrefixDelaunay, CompFCT, UpgradeBalinskiStat, fEquiv, MProc.GetO(TheId));
  auto FuncInsert=[&](Delaunay<T,Tint> const& x, std::ostream&os) -> int {
    DelaunayInv<T,Tint> eInv=ComputeInvariantDelaunay<T,Tint>(eData.GramMat, x);
    return ListOrbit.InsertEntry({x, eInv}, os);
  };
  int nbPresentOrbit=ListOrbit.GetNbEntry();
  if (nbPresentOrbit == 0) {
    MyMatrix<Tint> EXT_I=FindDelaunayPolytope<T,Tint>(eData.GramMat, eData.CVPmethod, MProc.GetO(TheId));
    int RetVal=FuncInsert({EXT_I}, MProc.GetO(TheId));
    if (RetVal != -1) {
      std::cerr << "RetVal=" << RetVal << "\n";
      std::cerr << "The first orbit should definitely be new\n";
      std::cerr << "Otherwise, we have a clear bug\n";
      throw TerminalException{1};
    }
  }
  MProc.GetO(TheId) << "We have inserted eInc\n";
  int nbSpannThread=0;
  std::condition_variable cv;
  std::mutex mtx_cv;
  auto WaitStuck=[&](int const& MyId) -> void {
    std::unique_lock<std::mutex> lk(mtx_cv);
    cv.wait(lk, [&]{return ListOrbit.IsStuck() == false;});
  };
  auto WaitComplete=[&](int const& MyId) -> void {
    std::unique_lock<std::mutex> lk(mtx_cv);
    cv.wait(lk, [&]{return nbSpannThread == 0;});
  };
  auto CompDelaunayAdjacency=[&](int const& MyId, std::ostream& os) -> void {
    int nbWork=0;
    nbSpannThread++;
    while(true) {
      bool testStuck=ListOrbit.IsStuck();
      if (testStuck) {
	WaitStuck(MyId);
      }
      int eEntry=ListOrbit.GetNonTreatedOrbit(os);
      bool IsComplete=ListOrbit.GetCompleteStatus();
      if (eEntry == -1)
	break;
      if (IsComplete)
	break;
      Delaunay<T,Tint> eDEL=ListOrbit.GetRepresentative(eEntry);
      DelaunayInv<T,Tint> eInv=ListOrbit.GetInvariant(eEntry);
      Tgroup GRPlatt=Delaunay_Stabilizer<T,Tint,Tgroup>(eData, eDEL.EXT, eInv.eIndex);
      MyMatrix<T> EXT_T=ConvertMatrixUniversal<T,Tint>(eDEL.EXT);
      CondTempDirectory eDir(AllArr.Saving, eData.PrefixPolyhedral + "ADM" + IntToString(eEntry) + "/");
      std::vector<Face> TheOutput=DUALDESC_THR_AdjacencyDecomposition(MProc, MyId, TheBank, EXT_T, GRPlatt, AllArr, eDir.str(), 0);
      for (auto& eOrbB : TheOutput) {
	MyMatrix<Tint> EXTadj=FindAdjacentDelaunayPolytope<T,Tint>(eData.GramMat, EXT_T, eOrbB, eData.CVPmethod);
	int eVal=FuncInsert({EXTadj}, os);
	if (eVal == -1)
	  cv.notify_one();
      }
      ListOrbit.SetEntryAsDone(eEntry, os);
      nbWork++;
    }
    if (MyId != TheId)
      MProc.MPU_Terminate(MyId);
    nbSpannThread--;
    if (nbWork > 0) {
      cv.notify_all();
    }
  };
  int NbThr=MProc.MPU_NumberFree();
  std::vector<std::thread> ListThreads(NbThr);
  MProc.GetO(TheId) << "We have NbThr=" << NbThr << "\n";
  while(true) {
    bool IsCompleteSpann=ListOrbit.GetCompleteStatus();
    MProc.GetO(TheId) << "IsCompleteSpann=" << IsCompleteSpann << "\n";
    if (IsCompleteSpann)
      WaitComplete(TheId);
    else {
      for (int iThr=0; iThr<NbThr; iThr++) {
	int NewId=MProc.MPU_GetId();
	MProc.GetO(TheId) << "NewId=" << NewId << "\n";
	if (NewId != -1) {
	  MProc.GetO(NewId) << "Before call to my_thread\n";
	  ListThreads[iThr]=std::thread(CompDelaunayAdjacency, NewId, std::ref(MProc.GetO(NewId)));
	  MProc.GetO(NewId) << "After call to my_thread\n";
	  ListThreads[iThr].detach();
	}
      }
      MProc.GetO(TheId) << "After the spanning of threads\n";
      MProc.GetO(TheId) << "Before CompDelaunayAdjacency\n";
      CompDelaunayAdjacency(TheId, MProc.GetO(TheId));
      MProc.GetO(TheId) << "After my own call to CompDelaunayAdjacency\n";
    }
    if (nbSpannThread == 0)
      break;
  }
  bool IsComplete=ListOrbit.GetCompleteStatus();
  MProc.GetO(TheId) << "IsComplete=" << IsComplete << "\n";
  if (!IsComplete) {
    std::cerr << "Major error in the code. We should be complete now\n";
    throw TerminalException{1};
  }
  std::vector<Delaunay<T,Tint>> ReturnData;
  if (eData.ReturnAll) {
    ReturnData=ListOrbit.GetListRepr();
    ListOrbit.FuncClear();
  }
  return ReturnData;
}



FullNamelist NAMELIST_GetStandard_COMPUTE_DELAUNAY()
{
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string> > ListListStringValues1;
  ListStringValues1["GRAMfile"]="unset.gram";
  ListStringValues1["SVRfile"]="unset.svr";
  ListStringValues1["OUTfile"]="unset.out";
  ListBoolValues1["SavingDelaunay"]=false;
  ListBoolValues1["FullDataInMemory"]=true;
  ListBoolValues1["ReturnAll"]=false;
  ListStringValues1["PrefixDelaunay"]="/irrelevant/";
  SingleBlock BlockDATA;
  BlockDATA.ListIntValues=ListIntValues1;
  BlockDATA.ListBoolValues=ListBoolValues1;
  BlockDATA.ListDoubleValues=ListDoubleValues1;
  BlockDATA.ListStringValues=ListStringValues1;
  BlockDATA.ListListStringValues=ListListStringValues1;
  ListBlock["DATA"]=BlockDATA;
  // METHOD
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string> > ListListStringValues2;
  ListIntValues2["NPROC"]=1;
  ListIntValues2["UpperLimitMethod4"]=10000;
  ListBoolValues2["Saving"]=false;
  ListStringValues2["PrefixPolyhedral"]="/irrelevant/";
  ListBoolValues2["FullDataInMemory"]=true;
  ListStringValues2["SplittingHeuristicFile"]="unset.heu";
  ListStringValues2["AdditionalSymmetryHeuristicFile"]="unset.heu";
  ListStringValues2["DualDescriptionHeuristicFile"]="unset.heu";
  ListStringValues2["StabEquivFacetHeuristicFile"]="unset.heu";
  ListStringValues2["MethodInitialFacetSetFile"]="unset.heu";
  ListStringValues2["CVPmethod"]="SVexact";
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
  ListStringValues3["BankSaveHeuristicFile"]="unset.heu";
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




template<typename T, typename Tint, typename Tgroup>
void TreatDelaunayEntry(FullNamelist const& eFull)
{
  SingleBlock BlockBANK=eFull.ListBlock.at("BANK");
  SingleBlock BlockDATA=eFull.ListBlock.at("DATA");
  SingleBlock BlockMETHOD=eFull.ListBlock.at("METHOD");
  //
  bool BANK_IsSaving=BlockBANK.ListBoolValues.at("Saving");
  bool BANK_Memory=BlockBANK.ListBoolValues.at("FullDataInMemory");
  std::string BANK_Prefix=BlockBANK.ListStringValues.at("Prefix");
  CreateDirectory(BANK_Prefix);
  FctsDataBank<PolyhedralEntry<T,Tgroup>> recFct=GetRec_FctsDataBank<T,Tgroup>();
  DataBank<PolyhedralEntry<T,Tgroup>> TheBank(BANK_IsSaving, BANK_Memory, BANK_Prefix, recFct);
  //
  std::cerr << "Reading DATA\n";
  std::string GRAMfile=BlockDATA.ListStringValues.at("GRAMfile");
  IsExistingFileDie(GRAMfile);
  std::cerr << "GRAMfile=" << GRAMfile << "\n";
  std::string SVRfile=BlockDATA.ListStringValues.at("SVRfile");
  IsExistingFileDie(SVRfile);
  std::cerr << "SVRfile=" << SVRfile << "\n";
  std::string OUTfile=BlockDATA.ListStringValues.at("OUTfile");
  std::cerr << "OUTfile=" << OUTfile << "\n";
  std::ifstream GRAMfs(GRAMfile);
  MyMatrix<T> GramMat=ReadMatrix<T>(GRAMfs);
  std::cerr << "dim(Gram)=" << GramMat.rows() << "\n";
  std::ifstream SVRfs(SVRfile);
  MyMatrix<Tint> SVR=ReadMatrix<Tint>(SVRfs);
  std::cerr << "|SVR|=" << SVR.rows() << "\n";
  bool Delaunay_Saving=BlockDATA.ListBoolValues.at("SavingDelaunay");
  bool Delaunay_Memory=BlockDATA.ListBoolValues.at("FullDataInMemory");
  bool Delaunay_ReturnAll=BlockDATA.ListBoolValues.at("ReturnAll");
  std::string PrefixDelaunay=BlockDATA.ListStringValues.at("PrefixDelaunay");
  CreateDirectory(PrefixDelaunay);
  std::string PrefixPolyhedral=BlockMETHOD.ListStringValues.at("PrefixPolyhedral");
  CreateDirectory(PrefixPolyhedral);
  int UppLimit=BlockMETHOD.ListIntValues.at("UpperLimitMethod4");
  std::string CVPmethod=BlockMETHOD.ListStringValues.at("CVPmethod");
  mpz_class UppLimit_z=UppLimit;
  int n=GramMat.rows();
  DataLattice<T,Tint> eData{n, GramMat, SVR, Delaunay_Saving, Delaunay_Memory, PrefixDelaunay, PrefixPolyhedral, Delaunay_ReturnAll, UppLimit_z, CVPmethod};
  //
  std::cerr << "Creating MPROC\n";
  int NbThr=BlockMETHOD.ListIntValues.at("NPROC");
  MainProcessor MProc(NbThr);
  int TheId=MProc.MPU_GetId();
  //
  PolyHeuristic<mpz_class> AllArr=AllStandardHeuristic<mpz_class>();
  //
  std::string HeuSplitFile=BlockMETHOD.ListStringValues.at("SplittingHeuristicFile");
  if (HeuSplitFile != "unset.heu") {
    IsExistingFileDie(HeuSplitFile);
    std::ifstream SPLITfs(HeuSplitFile);
    AllArr.Splitting=ReadHeuristic<mpz_class>(SPLITfs);
  }
  //
  std::string AddiSymmFile=BlockMETHOD.ListStringValues.at("AdditionalSymmetryHeuristicFile");
  if (AddiSymmFile != "unset.heu") {
    IsExistingFileDie(AddiSymmFile);
    std::ifstream SYMMfs(AddiSymmFile);
    AllArr.AdditionalSymmetry=ReadHeuristic<mpz_class>(SYMMfs);
  }
  //
  std::string DualDescFile=BlockMETHOD.ListStringValues.at("DualDescriptionHeuristicFile");
  if (DualDescFile != "unset.heu") {
    IsExistingFileDie(DualDescFile);
    std::ifstream DDfs(DualDescFile);
    AllArr.DualDescriptionProgram=ReadHeuristic<mpz_class>(DDfs);
  }
  //
  std::string StabEquivFile=BlockMETHOD.ListStringValues.at("StabEquivFacetHeuristicFile");
  if (StabEquivFile != "unset.heu") {
    IsExistingFileDie(StabEquivFile);
    std::ifstream STEQfs(StabEquivFile);
    AllArr.StabEquivFacet=ReadHeuristic<mpz_class>(STEQfs);
  }
  //
  std::string MethodFacetFile=BlockMETHOD.ListStringValues.at("MethodInitialFacetSetFile");
  if (MethodFacetFile != "unset.heu") {
    IsExistingFileDie(MethodFacetFile);
    std::ifstream MIFSfs(MethodFacetFile);
    AllArr.InitialFacetSet=ReadHeuristic<mpz_class>(MIFSfs);
  }
  //
  std::string BankSaveFile=BlockBANK.ListStringValues.at("BankSaveHeuristicFile");
  if (BankSaveFile != "unset.heu") {
    IsExistingFileDie(BankSaveFile);
    std::ifstream BANKfs(BankSaveFile);
    AllArr.BankSave=ReadHeuristic<mpz_class>(BANKfs);
  }
  //
  bool DD_Saving=BlockMETHOD.ListBoolValues.at("Saving");
  bool DD_Memory=BlockMETHOD.ListBoolValues.at("FullDataInMemory");
  AllArr.Saving=DD_Saving;
  AllArr.eMemory=DD_Memory;
  //
  std::vector<Delaunay<T,Tint>> ListDel=EnumerationDelaunayPolytopes<T>(MProc, TheId,
								   TheBank, eData, AllArr);
  std::cerr << "We now have ListDel\n";
  //
  if (Delaunay_ReturnAll) {
    std::ofstream OUTfs(OUTfile);
    int nbDel=ListDel.size();
    OUTfs << "nbDel=" << nbDel << "\n";
    for (int iDel=0; iDel<nbDel; iDel++) {
      OUTfs << "iDel=" << iDel << "/" << nbDel << "\n";
      WriteMatrix(OUTfs, ListDel[iDel].EXT);
    }
  }
}




#endif

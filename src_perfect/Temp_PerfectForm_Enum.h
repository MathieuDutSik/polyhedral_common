// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_TEMP_PERFECTFORM_ENUM_H_
#define SRC_PERFECT_TEMP_PERFECTFORM_ENUM_H_

// clang-format off
#include "POLY_ThreadDualDescription.h"
#include "Parallel_Classes.h"
#include "Temp_PerfectForm.h"
#include <map>
#include <string>
#include <vector>
// clang-format on

template <typename T> struct DataLinSpa {
  LinSpaceMatrix<T> LinSpa;
  bool SavingPerfect;
  bool FullDataInMemory;
  std::string PrefixPerfect;
  std::string PrefixPolyhedral;
  bool ReturnAll;
  mpz_class UpperLimitMethod4;
  bool NeedCheckStabilization;
};

template <typename T, typename Tint, typename Tgroup> struct ListPerfectForm {
public:
  QueryEquivInfo IsPresentNoLock(int const &iOrbitStart,
                                 LinSpaceMatrix<T> const &LinSpa,
                                 MyMatrix<T> const &eGram,
                                 MyMatrix<int> const &SHV) {
    int nbOrbit = ListPerfect.size();
    Face eTrivFace;
    for (int iOrbit = iOrbitStart; iOrbit < nbOrbit; iOrbit++) {
      MyMatrix<T> oldGram = GetGram(iOrbit);
      MyMatrix<int> oldSHV = GetSHV(iOrbit);
      std::optional<MyMatrix<int>> eLinSpaRes =
          PERF_TestEquivalence<T, Tint, Tgroup>(LinSpa, eGram, oldGram, SHV,
                                                oldSHV);
      if (eLinSpaRes) {
        PerfEquivInfo eEquiv{iOrbit, *eLinSpaRes, eTrivFace};
        return {true, nbOrbit, eEquiv};
      }
    }
    return {false, nbOrbit, {-1, {}, eTrivFace}};
  }
  PerfEquivInfo InsertForm(LinSpaceMatrix<T> const &LinSpa,
                           MyMatrix<T> const &eGram, std::ostream &os) {
    std::vector<PerfEquivInfo> eListEquivInfo;
    Face eTrivFace;
    int n = LinSpa.n;
    Tshortest<T, int> RecSHV = T_ShortestVector<T, int>(eGram, os);
    QueryEquivInfo eQuery = IsPresentNoLock(0, LinSpa, eGram, RecSHV.SHV);
    if (eQuery.result)
      return eQuery.eEquiv;
    std::lock_guard<std::mutex> lk(mul);
    QueryEquivInfo fQuery =
        IsPresentNoLock(eQuery.nbOrbit, LinSpa, eGram, RecSHV.SHV);
    if (fQuery.result)
      return fQuery.eEquiv;
    SinglePerfect<T, Tint, Tgroup> eSing =
        GetPerfectCone(LinSpa, eGram, RecSHV);
    ListPerfect.push_back(eSing);
    int len = ListPerfect.size();
    MyMatrix<int> eMatEquiv = IdentityMat<int>(n);
    return {len, eMatEquiv, eTrivFace};
  }
  std::vector<int> GetListOrbitWorkEff(std::vector<int> const &ListOrbWork) {
    std::vector<int> ListOrbitWorkEff;
    std::lock_guard<std::mutex> lk(mul);
    for (auto eOrb : ListOrbWork) {
      if (ListPerfect[eOrb].eCons == 0) {
        ListPerfect[eOrb].eCons = 1;
        ListOrbitWorkEff.push_back(eOrb);
      }
    }
    return ListOrbitWorkEff;
  }
  void SetStatus(int const &eOrb, int const &eVal) {
    ListPerfect[eOrb].eStatus = eVal;
  }
  MyMatrix<int> GetSHV(int eOrb) const { return ListPerfect[eOrb].SHV; }
  MyMatrix<T> GetPerfDomEXT(int const &eOrb) const {
    return ListPerfect[eOrb].PerfDomEXT;
  }
  Tgroup GetPerfDomGRP(int const &eOrb) const {
    return ListPerfect[eOrb].PerfDomGRP;
  }
  MyMatrix<T> GetGram(int eOrb) const { return ListPerfect[eOrb].eGram; }
  void InsertEquivInfo(int eOrb, PerfEquivInfo const &eEquiv) {
    ListPerfect[eOrb].ListEquivInfo.push_back(eEquiv);
  }
  int GetNbPerf() const {
    int len = ListPerfect.size();
    return len;
  }
  std::vector<std::vector<int>>
  VoronoiAlgo_THR_GetPartition(int NbThr, int &nbConsTodo, int &IsFinished) {
    int iOrbit, iThr;
    int nbOrbit = ListPerfect.size();
    std::vector<std::vector<int>> ThePartition;
    nbConsTodo = 0;
    if (NbThr == 0) {
      std::cerr << "Error in NbThr\n";
      std::cerr << "NbThr=" << NbThr << "\n";
      throw TerminalException{1};
    }
    std::vector<int> eNewList;
    for (iThr = 0; iThr < NbThr; iThr++)
      ThePartition.push_back(eNewList);
    iThr = 0;
    IsFinished = 1;
    for (iOrbit = 0; iOrbit < nbOrbit; iOrbit++) {
      if (ListPerfect[iOrbit].eCons == 0) {
        nbConsTodo++;
        ThePartition[iThr].push_back(iOrbit);
        iThr++;
        if (iThr == NbThr)
          iThr = 0;
      }
      if (ListPerfect[iOrbit].eStatus == 0)
        IsFinished = 0;
    }
    return ThePartition;
  }

private:
  std::mutex mul;
  std::vector<SinglePerfect<T, Tint, Tgroup>> ListPerfect;
};

template <typename T, typename Tint, typename Tgroup>
void VoronoiAlgo_PrintListMatrix(
    std::ostream &os, [[maybe_unused]] LinSpaceMatrix<T> const &LinSpa,
    ListPerfectForm<T, Tint, Tgroup> const &ListPerf) {
  int nbGram = ListPerf.GetNbPerf();
  os << nbGram << "\n";
  for (int iGram = 0; iGram < nbGram; iGram++) {
    os << iGram << " " << nbGram << "\n";
    MyMatrix<T> eGram = ListPerf.GetGram(iGram);
    MyMatrix<int> eSHV = ListPerf.GetSHV(iGram);
    WriteMatrix(os, eGram);
    WriteMatrix(os, eSHV);
  }
}

template <typename T, typename Tint, typename Tgroup>
void VoronoiAlgo_THR_BlockTreatment(
    MainProcessor &MProc, int TheId,
    DataBank<PolyhedralEntry<T, Tgroup>> &TheBank,
    LinSpaceMatrix<T> const &LinSpa, ListPerfectForm<T, Tint, Tgroup> &ListPerf,
    PolyHeuristic<mpz_class> const &AllArr,
    std::vector<int> const &ListOrbWork) {
  std::vector<int> ListOrbitWorkEff = ListPerf.GetListOrbitWorkEff(ListOrbWork);
  for (auto &eOrb : ListOrbitWorkEff) {
    MyMatrix<T> eGram = ListPerf.GetGram(eOrb);
    MyMatrix<T> PerfDomEXT = ListPerf.GetPerfDomEXT(eOrb);
    Tgroup PerfDomGRP = ListPerf.GetPerfDomGRP(eOrb);
    vectface TheOutput = DUALDESC_THR_AdjacencyDecomposition(
        MProc, TheId, TheBank, PerfDomEXT, PerfDomGRP, AllArr);
    for (auto &eOrbB : TheOutput) {
      MyVector<T> eVectOrb = FindFacetInequality(PerfDomEXT, eOrbB);
      MyMatrix<T> DirMat = LINSPA_GetMatrixInTspace(LinSpa, eVectOrb);
      MyMatrix<T> NewPerf = Flipping_Perfect<T, int>(eGram, DirMat).first;
      PerfEquivInfo eEquiv = ListPerf.InsertForm(LinSpa, NewPerf);
      eEquiv.eInc = eOrbB;
      ListPerf.InsertEquivInfo(eOrb, eEquiv);
    }
    ListPerf.SetStatus(eOrb, 1);
  }
}

template <typename T, typename Tint, typename Tgroup>
void VoronoiAlgo_THR_EnumeratePerfectForm(
    MainProcessor &MProc, int TheId,
    DataBank<PolyhedralEntry<T, Tgroup>> &TheBank,
    LinSpaceMatrix<T> const &LinSpa, PolyHeuristic<mpz_class> const &AllArr,
    ListPerfectForm<T, Tint, Tgroup> &ListPerf) {
  int IsFinished, NewId;
  int nbConsTodo, iThr, NbThr;
  MyMatrix<T> ePerf = GetOnePerfectForm<T>(LinSpa);
  PerfEquivInfo eEquiv = ListPerf.InsertForm(LinSpa, ePerf);
  while (true) {
    NbThr = MProc.MPU_NumberFree();
    std::vector<std::vector<int>> ThePartition =
        ListPerf.VoronoiAlgo_THR_GetPartition(NbThr + 1, nbConsTodo,
                                              IsFinished);
    if (IsFinished == 1)
      break;
    for (iThr = 1; iThr < NbThr; iThr++)
      if (ThePartition[iThr].size() > 0) {
        NewId = MProc.MPU_GetId();
        if (NewId != -1) {
          std::thread my_thread(VoronoiAlgo_THR_BlockTreatment<T>,
                                std::ref(MProc), NewId, std::ref(TheBank),
                                std::ref(LinSpa), std::ref(ListPerf),
                                std::ref(AllArr), std::ref(ThePartition[iThr]));
          my_thread.detach();
          MProc.MPU_Terminate(NewId);
        }
      }
    if (ThePartition[0].size() > 0)
      VoronoiAlgo_THR_BlockTreatment<T>(MProc, TheId, TheBank, LinSpa, ListPerf,
                                        AllArr, ThePartition[0]);
  }
}

template <typename T, typename Tint, typename Tgroup>
std::vector<SimplePerfect<T, Tint>>
EnumerationPerfectMatrices(MainProcessor &MProc, int const &TheId,
                           DataBank<PolyhedralEntry<T, Tgroup>> &TheBank,
                           DataLinSpa<T> const &eData,
                           PolyHeuristic<mpz_class> const &AllArr) {
  std::function<bool(SimplePerfectInv<T> const &, SimplePerfectInv<T> const &)>
      CompFCT = [](SimplePerfectInv<T> const &x,
                   SimplePerfectInv<T> const &y) -> bool { return x < y; };
  std::function<void(TrivialBalinski &, SimplePerfect<T, Tint> const &,
                     SimplePerfectInv<T> const &, std::ostream &)>
      UpgradeBalinskiStat =
          []([[maybe_unused]] TrivialBalinski const &eStat,
             [[maybe_unused]] SimplePerfect<T, Tint> const &eEnt,
             [[maybe_unused]] SimplePerfectInv<T> const &eInv,
             [[maybe_unused]] std::ostream &os) -> void {};
  std::function<std::optional<MyMatrix<Tint>>(SimplePerfect<T, Tint> const &,
                                              SimplePerfect<T, Tint> const &)>
      fEquiv =
          [&](SimplePerfect<T, Tint> const &x, SimplePerfect<T, Tint> const &y)
      -> std::optional<MyMatrix<Tint>> {
    return SimplePerfect_TestEquivalence<T, Tint, Tgroup>(eData, x.Gram,
                                                          y.Gram);
  };
  NewEnumerationWork<SimplePerfect<T, Tint>> ListOrbit(
      AllArr.DD_Saving, AllArr.DD_Memory, eData.PrefixPerfect, CompFCT,
      UpgradeBalinskiStat, fEquiv, MProc.GetO(TheId));
  auto FuncInsert = [&](SimplePerfect<T, Tint> const &x,
                        std::ostream &os) -> int {
    SimplePerfectInv<T> eInv = ComputeInvariantSimplePerfect<T, Tint>(x.Gram);
    return ListOrbit.InsertEntry({x, eInv}, os);
  };
  int nbPresentOrbit = ListOrbit.GetNbEntry();
  if (nbPresentOrbit == 0) {
    MyMatrix<T> eGram = GetOnePerfectForm(eData.LinSpa);
    int RetVal = FuncInsert({eGram}, MProc.GetO(TheId));
    if (RetVal != -1) {
      std::cerr << "RetVal=" << RetVal << "\n";
      std::cerr << "The first orbit should definitely be new\n";
      std::cerr << "Otherwise, we have a clear bug\n";
      throw TerminalException{1};
    }
  }
  MProc.GetO(TheId) << "We have inserted eInc\n";
  int nbSpannThread = 0;
  std::condition_variable cv;
  std::mutex mtx_cv;
  auto WaitStuck = [&]([[maybe_unused]] int const &MyId) -> void {
    std::unique_lock<std::mutex> lk(mtx_cv);
    cv.wait(lk, [&] { return ListOrbit.IsStuck() == false; });
  };
  auto WaitComplete = [&]([[maybe_unused]] int const &MyId) -> void {
    std::unique_lock<std::mutex> lk(mtx_cv);
    cv.wait(lk, [&] { return nbSpannThread == 0; });
  };
  auto CompSimplePerfectAdjacency = [&](int const &MyId,
                                        std::ostream &os) -> void {
    int nbWork = 0;
    nbSpannThread++;
    while (true) {
      bool testStuck = ListOrbit.IsStuck();
      if (testStuck) {
        WaitStuck(MyId);
      }
      int eEntry = ListOrbit.GetNonTreatedOrbit(os);
      bool IsComplete = ListOrbit.GetCompleteStatus();
      if (eEntry == -1)
        break;
      if (IsComplete)
        break;
      SimplePerfect<T, Tint> ePERF = ListOrbit.GetRepresentative(eEntry);
      Tshortest<T, Tint> RecSHV = T_ShortestVector<T, Tint>(ePERF.Gram, os);
      NakedPerfect<T, Tint> eNaked =
          GetNakedPerfectCone(eData.LinSpa, ePERF.Gram, RecSHV);
      Tgroup GRPshv =
          SimplePerfect_Stabilizer<T, Tint, Tgroup>(eData, ePERF.Gram, RecSHV);
      Tgroup PerfDomGRP = MapLatticeGroupToConeGroup(eNaked, GRPshv);
      CondTempDirectory eDir(AllArr.DD_Saving, eData.PrefixPolyhedral + "ADM" +
                                                   IntToString(eEntry) + "/");
      vectface TheOutput = DUALDESC_THR_AdjacencyDecomposition(
          MProc, MyId, TheBank, eNaked.PerfDomEXT, PerfDomGRP, AllArr,
          eDir.str(), 0);
      for (auto &eOrbB : TheOutput) {
        MyVector<T> eVectOrb = FindFacetInequality(eNaked.PerfDomEXT, eOrbB);
        MyMatrix<T> DirMat = LINSPA_GetMatrixInTspace(eData.LinSpa, eVectOrb);
        MyMatrix<T> NewPerf =
            Flipping_Perfect<T, int>(ePERF.Gram, DirMat).first;
        int eVal = FuncInsert({NewPerf}, os);
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
  int NbThr = MProc.MPU_NumberFree();
  std::vector<std::thread> ListThreads(NbThr);
  MProc.GetO(TheId) << "We have NbThr=" << NbThr << "\n";
  while (true) {
    bool IsCompleteSpann = ListOrbit.GetCompleteStatus();
    MProc.GetO(TheId) << "IsCompleteSpann=" << IsCompleteSpann << "\n";
    if (IsCompleteSpann) {
      WaitComplete(TheId);
    } else {
      for (int iThr = 0; iThr < NbThr; iThr++) {
        int NewId = MProc.MPU_GetId();
        MProc.GetO(TheId) << "NewId=" << NewId << "\n";
        if (NewId != -1) {
          MProc.GetO(NewId) << "Before call to my_thread\n";
          ListThreads[iThr] = std::thread(CompSimplePerfectAdjacency, NewId,
                                          std::ref(MProc.GetO(NewId)));
          MProc.GetO(NewId) << "After call to my_thread\n";
          ListThreads[iThr].detach();
        }
      }
      CompSimplePerfectAdjacency(TheId, MProc.GetO(TheId));
    }
    if (nbSpannThread == 0)
      break;
  }
  bool IsComplete = ListOrbit.GetCompleteStatus();
  MProc.GetO(TheId) << "IsComplete=" << IsComplete << "\n";
  if (!IsComplete) {
    std::cerr << "Major error in the code. We should be complete now\n";
    throw TerminalException{1};
  }
  std::vector<SimplePerfect<T, Tint>> ReturnData;
  if (eData.ReturnAll) {
    ReturnData = ListOrbit.GetListRepr();
    ListOrbit.FuncClear();
  }
  return ReturnData;
}

FullNamelist NAMELIST_GetStandard_COMPUTE_PERFECT() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListStringValues1["LinSpaFile"] = "unset.gram";
  ListStringValues1["OUTfile"] = "unset.out";
  ListBoolValues1["SavingPerfect"] = false;
  ListBoolValues1["NeedCheckStabilization"] = false;
  ListBoolValues1["FullDataInMemory"] = true;
  ListBoolValues1["ReturnAll"] = false;
  ListStringValues1["PrefixPerfect"] = "/irrelevant/";
  SingleBlock BlockDATA;
  BlockDATA.setListIntValues(ListIntValues1);
  BlockDATA.setListBoolValues(ListBoolValues1);
  BlockDATA.setListStringValues(ListStringValues1);
  ListBlock["DATA"] = BlockDATA;
  // METHOD
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, std::string> ListStringValues2;
  ListIntValues2["NPROC"] = 1;
  ListIntValues2["UpperLimitMethod4"] = 10000;
  ListBoolValues2["Saving"] = false;
  ListBoolValues2["FullDataInMemory"] = true;
  ListStringValues2["PrefixPolyhedral"] = "/irrelevant/";
  ListStringValues2["SplittingHeuristicFile"] = "unset.heu";
  ListStringValues2["AdditionalSymmetryHeuristicFile"] = "unset.heu";
  ListStringValues2["DualDescriptionHeuristicFile"] = "unset.heu";
  ListStringValues2["StabEquivFacetHeuristicFile"] = "unset.heu";
  ListStringValues2["MethodInitialFacetSetFile"] = "unset.heu";
  SingleBlock BlockMETHOD;
  BlockMETHOD.setListIntValues(ListIntValues2);
  BlockMETHOD.setListBoolValues(ListBoolValues2);
  BlockMETHOD.setListStringValues(ListStringValues2);
  ListBlock["METHOD"] = BlockMETHOD;
  // BANK
  std::map<std::string, bool> ListBoolValues3;
  std::map<std::string, std::string> ListStringValues3;
  ListStringValues3["BankSaveHeuristicFile"] = "unset.heu";
  ListStringValues3["Prefix"] = "./unset/";
  ListBoolValues3["Saving"] = false;
  ListBoolValues3["FullDataInMemory"] = true;
  SingleBlock BlockBANK;
  BlockBANK.setListBoolValues(ListBoolValues3);
  BlockBANK.setListStringValues(ListStringValues3);
  ListBlock["BANK"] = BlockBANK;
  // Merging all data
  return FullNamelist(ListBlock);
}

template <typename T, typename Tint, typename Tgroup>
void TreatPerfectLatticesEntry(FullNamelist const &eFull) {
  SingleBlock const& BlockBANK = eFull.get_block("BANK");
  SingleBlock const& BlockDATA = eFull.get_block("DATA");
  SingleBlock const& BlockMETHOD = eFull.get_block("METHOD");
  //
  bool BANK_Saving = BlockBANK.get_bool("Saving");
  bool BANK_Memory = BlockBANK.get_bool("FullDataInMemory");
  std::string BANK_Prefix = BlockBANK.get_string("Prefix");
  CreateDirectory(BANK_Prefix);
  FctsDataBank<PolyhedralEntry<T, Tgroup>> recFct =
      GetRec_FctsDataBank<T, Tgroup>();
  DataBank<PolyhedralEntry<T, Tgroup>> TheBank(BANK_Saving, BANK_Memory,
                                               BANK_Prefix, recFct);
  //
  std::cerr << "Reading DATA\n";
  std::string LINSPAfile = BlockDATA.get_string("LinSpaFile");
  IsExistingFileDie(LINSPAfile);
  std::cerr << "LINSPAfile=" << LINSPAfile << "\n";
  std::string OUTfile = BlockDATA.get_string("OUTfile");
  std::cerr << "OUTfile=" << OUTfile << "\n";
  //
  std::ifstream LINSPAfs(LINSPAfile);
  LinSpaceMatrix<T> LinSpa;
  LINSPAfs >> LinSpa;
  bool NeedCheckStabilization =
    BlockDATA.get_bool("NeedCheckStabilization");
  bool SavingPerfect = BlockDATA.get_bool("SavingPerfect");
  bool FullDataInMemory = BlockDATA.get_bool("FullDataInMemory");
  bool ReturnAll = BlockDATA.get_bool("ReturnAll");
  std::string PrefixPerfect = BlockDATA.get_string("PrefixPerfect");
  CreateDirectory(PrefixPerfect);
  std::string PrefixPolyhedral =
    BlockMETHOD.get_string("PrefixPolyhedral");
  CreateDirectory(PrefixPolyhedral);
  int UpperLimitMethod4 = BlockMETHOD.get_int("UpperLimitMethod4");
  std::string CVPmethod = BlockMETHOD.get_string("CVPmethod");
  mpz_class UpperLimitMethod4_z = UpperLimitMethod4;
  DataLinSpa<T> eData;
  eData.LinSpa = LinSpa;
  eData.SavingPerfect = SavingPerfect;
  eData.FullDataInMemory = FullDataInMemory;
  eData.PrefixPerfect = PrefixPerfect;
  eData.PrefixPolyhedral = PrefixPolyhedral;
  eData.ReturnAll = ReturnAll;
  eData.UpperLimitMethod4 = UpperLimitMethod4_z;
  eData.NeedCheckStabilization = NeedCheckStabilization;
  //
  std::cerr << "Creating MPROC\n";
  int NbThr = BlockMETHOD.get_int("NPROC");
  MainProcessor MProc(NbThr);
  int TheId = MProc.MPU_GetId();
  //
  PolyHeuristic<mpz_class> AllArr = AllStandardHeuristic<mpz_class>();
  //
  std::string HeuSplitFile =
    BlockMETHOD.get_string("SplittingHeuristicFile");
  if (HeuSplitFile != "unset.heu") {
    IsExistingFileDie(HeuSplitFile);
    std::ifstream SPLITfs(HeuSplitFile);
    AllArr.Splitting = ReadHeuristic<mpz_class>(SPLITfs);
  }
  //
  std::string AddiSymmFile =
    BlockMETHOD.get_string("AdditionalSymmetryHeuristicFile");
  if (AddiSymmFile != "unset.heu") {
    IsExistingFileDie(AddiSymmFile);
    std::ifstream SYMMfs(AddiSymmFile);
    AllArr.AdditionalSymmetry = ReadHeuristic<mpz_class>(SYMMfs);
  }
  //
  std::string DualDescFile =
    BlockMETHOD.get_string("DualDescriptionHeuristicFile");
  if (DualDescFile != "unset.heu") {
    IsExistingFileDie(DualDescFile);
    std::ifstream DDfs(DualDescFile);
    AllArr.DualDescriptionProgram = ReadHeuristic<mpz_class>(DDfs);
  }
  //
  std::string StabEquivFile =
    BlockMETHOD.get_string("StabEquivFacetHeuristicFile");
  if (StabEquivFile != "unset.heu") {
    IsExistingFileDie(StabEquivFile);
    std::ifstream STEQfs(StabEquivFile);
    AllArr.StabEquivFacet = ReadHeuristic<mpz_class>(STEQfs);
  }
  //
  std::string MethodFacetFile =
    BlockMETHOD.get_string("MethodInitialFacetSetFile");
  if (MethodFacetFile != "unset.heu") {
    IsExistingFileDie(MethodFacetFile);
    std::ifstream MIFSfs(MethodFacetFile);
    AllArr.InitialFacetSet = ReadHeuristic<mpz_class>(MIFSfs);
  }
  //
  std::string BankSaveFile =
    BlockBANK.get_string("BankSaveHeuristicFile");
  if (BankSaveFile != "unset.heu") {
    IsExistingFileDie(BankSaveFile);
    std::ifstream BANKfs(BankSaveFile);
    AllArr.BankSave = ReadHeuristic<mpz_class>(BANKfs);
  }
  //
  bool DD_Saving = BlockMETHOD.get_bool("Saving");
  bool DD_Memory = BlockMETHOD.get_bool("FullDataInMemory");
  AllArr.DD_Saving = DD_Saving;
  AllArr.DD_Memory = DD_Memory;
  //
  std::vector<SimplePerfect<T, Tint>> LPerf =
      EnumerationPerfectMatrices<T, Tint>(MProc, TheId, TheBank, eData, AllArr);
  std::cerr << "We now have ListDel\n";
  //
  if (ReturnAll) {
    std::ofstream OUTfs(OUTfile);
    int nbPerf = LPerf.size();
    OUTfs << "nbPerf=" << nbPerf << "\n";
    for (int iPerf = 0; iPerf < nbPerf; iPerf++) {
      OUTfs << "iPerf=" << iPerf << "/" << nbPerf << "\n";
      WriteMatrix(OUTfs, LPerf[iPerf].Gram);
    }
  }
}

// clang-format off
#endif  // SRC_PERFECT_TEMP_PERFECTFORM_ENUM_H_
// clang-format on

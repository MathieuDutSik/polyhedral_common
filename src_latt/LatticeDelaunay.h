// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEDELAUNAY_H_
#define SRC_LATT_LATTICEDELAUNAY_H_

// clang-format off
#include "FundamentalDelaunay.h"
#include "GRP_DoubleCoset.h"
#include "MatrixGroup.h"
#include "Namelist.h"
#include "POLY_ThreadDualDescription.h"
#include <map>
#include <string>
#include <vector>
// clang-format on

template <typename T, typename Tint> struct DataLattice {
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

template <typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value>
GetWeightMatrixFromGramEXT(MyMatrix<T> const &EXT, MyMatrix<T> const &GramMat,
                           MyMatrix<T> const &SHV) {
  int n = GramMat.rows();
  CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral(GramMat, EXT);
  MyVector<T> TheCenter = eCP.eCent;
  std::vector<MyMatrix<T>> CharPair = CharacterizingPair(GramMat, TheCenter);
  MyMatrix<T> Qmat = CharPair[0];
  int nbVect = SHV.rows();
  int nbVert = EXT.rows();
  MyMatrix<T> EXText(nbVect + nbVert, n + 1);
  for (int iVert = 0; iVert < nbVert; iVert++)
    EXText.row(iVert) = EXT.row(iVert);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    EXText(nbVert + iVect, 0) = 0;
    for (int i = 0; i < n; i++)
      EXText(nbVert + iVect, i + 1) = SHV(iVect, i);
  }
  return GetSimpleWeightMatrix<T, Tidx_value>(EXText, Qmat);
}

template <typename T, typename Tgroup>
bool IsGroupCorrect(MyMatrix<T> const &EXT_T, Tgroup const &eGRP) {
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> LGen = eGRP.GeneratorsOfGroup();
  for (auto &eGen : LGen) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    if (!IsIntegralMatrix(eMat))
      return false;
  }
  return true;
}

template <typename T, typename Tint, typename Tgroup>
Tgroup Delaunay_Stabilizer(DataLattice<T, Tint> const &eData,
                           MyMatrix<Tint> const &EXT, Tint const &eIndex) {
  using Tidx_value = int16_t;
  using Tgr = GraphListAdj;
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(EXT);
  WeightMatrix<true, T, Tidx_value> WMatRed =
      GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT_T, eData.GramMat, {});
  Tgroup GRPisomRed =
      GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMatRed);
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
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(eData.SHV);
  WeightMatrix<true, T, Tidx_value> WMat =
      GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT_T, eData.GramMat, SHV_T);
  int nbVert = EXT_T.rows();
  int nbSHV = SHV_T.rows();
  Face eFace(nbVert + nbSHV);
  for (int iVert = 0; iVert < nbVert; iVert++)
    eFace[iVert] = 1;
  for (int iSHV = 0; iSHV < nbSHV; iSHV++)
    eFace[nbVert + iSHV] = 0;
  Tgroup PreGRPisom =
      GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat);
  Tgroup GRPisom = ReducedGroupAction(PreGRPisom, eFace);
  if (IsGroupCorrect(EXT_T, GRPisom)) {
    return GRPisom;
  }
  if (GRPisom.size() < eData.UpperLimitMethod4) {
    auto f_correct = [](MyMatrix<T> const &M) -> bool {
      return IsIntegralMatrix(M);
    };
    return LinPolytopeIntegral_Stabilizer_Method4(EXT_T, GRPisom, f_correct);
  } else {
    return LinPolytopeIntegral_Stabilizer_Method8(EXT_T, GRPisom);
  }
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
Delaunay_TestEquivalence(DataLattice<T, Tint> const &eData,
                         Delaunay<T, Tint> const &RecEXT1,
                         Delaunay<T, Tint> const &RecEXT2, Tint const &eIndex) {
  using Telt = typename Tgroup::Telt;
  using Tgr = GraphListAdj;
  using Tidx_value = int16_t;
  std::cerr << "Begin Delaunay_TestEquivalence\n";
  auto ConvertEquiv = [](std::optional<MyMatrix<T>> const &eEq)
      -> std::optional<MyMatrix<Tint>> {
    if (!eEq)
      return {};
    MyMatrix<Tint> eMat_I = UniversalMatrixConversion<Tint, T>(eEq.TheEquiv);
    return eMat_I;
  };
  MyMatrix<T> EXT1_T = UniversalMatrixConversion<T, Tint>(RecEXT1.EXT);
  MyMatrix<T> EXT2_T = UniversalMatrixConversion<T, Tint>(RecEXT2.EXT);
  WeightMatrix<true, T, Tidx_value> WMatRed1 =
      GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT1_T, eData.GramMat, {});
  WeightMatrix<true, T, Tidx_value> WMatRed2 =
      GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT2_T, eData.GramMat, {});
  std::optional<Telt> eResRed =
      TestEquivalenceWeightMatrix<T, Telt, Tidx_value>(WMatRed1, WMatRed2);
  if (!eResRed) {
    std::cerr << "Leaving Delaunay_TestEquivalence with false\n";
    return {};
  }
  MyMatrix<T> MatEquivRed_T =
      FindTransformation(EXT1_T, EXT2_T, eResRed.TheEquiv);
  if (IsIntegralMatrix(MatEquivRed_T)) {
    MyMatrix<Tint> MatEquiv_I =
        UniversalMatrixConversion<Tint, T>(MatEquivRed_T);
    std::cerr << "Leaving Delaunay_TestEquivalence with true\n";
    return MatEquiv_I;
  }
  if (eIndex == 1) {
    std::cerr << "We should not reach that stage. If eIndex=1 then\n";
    std::cerr << "the basic algorithm should work\n";
    throw TerminalException{1};
  }
  //
  // Now extending by adding more vectors.
  //
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(eData.SHV);
  WeightMatrix<true, T, Tidx_value> WMat1 =
      GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT1_T, eData.GramMat, SHV_T);
  WeightMatrix<true, T, Tidx_value> WMat2 =
      GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT2_T, eData.GramMat, SHV_T);
  std::optional<Telt> eRes =
      TestEquivalenceWeightMatrix<T, Telt, Tidx_value>(WMat1, WMat2);
  if (!eRes) {
    std::cerr << "Leaving Delaunay_TestEquivalence with false\n";
    return {};
  }
  MyMatrix<T> MatEquiv_T = FindTransformation(EXT1_T, EXT2_T, eRes.TheEquiv);
  if (IsIntegralMatrix(MatEquiv_T)) {
    MyMatrix<Tint> MatEquiv_I = UniversalMatrixConversion<Tint, T>(MatEquiv_T);
    std::cerr << "Leaving Delaunay_TestEquivalence with true\n";
    return {true, MatEquiv_I};
  }
  std::cerr << "Trying other strategies\n";
  Tgroup GRP1 = GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat1);
  if (GRP1.size() < eData.UpperLimitMethod4) {
    auto f_correct = [](MyMatrix<T> const &M) -> bool {
      return IsIntegralMatrix(M);
    };
    return ConvertEquiv(LinPolytopeIntegral_Isomorphism_Method4(
        EXT1_T, EXT2_T, GRP1, eRes.TheEquiv, f_correct));
  }
  return ConvertEquiv(LinPolytopeIntegral_Isomorphism_Method8(
      EXT1_T, EXT2_T, GRP1, eRes.TheEquiv));
}

template <typename T, typename Tint>
size_t ComputeInvariantDelaunay(DataLattice<T, Tint> const &eData,
                                size_t const& seed,
                                MyMatrix<Tint> const& EXT) {
  using Tidx_value = int16_t;
  int nbVert = EXT.rows();
  int n = EXT.cols() - 1;
  Tint PreIndex = Int_IndexLattice(EXT);
  Tint eIndex = T_abs(PreIndex);
  MyVector<T> V(n);
  CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral(eData.GramMat, EXT);
  MyMatrix<T> EXT_T(nbVert, n);
  for (int iVert=0; iVert<nbVert; iVert++) {
    for (int i=0; i<n; i++) {
      T val = UniversalScalarConversion<T,Tint>(EXT(iVert, i+1));
      EXT_T(iVert, i) = val - eCP.eCent(i+1);
    }
  }
  std::map<T, size_t> ListDiagNorm;
  std::map<T, size_t> ListOffDiagNorm;
  for (int iVert = 0; iVert < nbVert; iVert++) {
    MyVector<T> V(n);
    for (int i=0; i<n; i++) {
      T eSum = 0;
      for (int j=0; j<n; j++) {
        eSum += eData.GramMat(i,j) * EXT_T(iVert, j);
      }
      V(i) = eSum;
    }
    T scal = 0;
    for (int i=0; i<n; i++) {
      scal += V(i) * EXT_T(iVert, i);
    }
    ListDiagNorm[scal] += 1;
    for (int jVert=iVert+1; jVert<nbVert; jVert++) {
      T scal = 0;
      for (int i=0; i<n; i++) {
        scal += V(i) * EXT_T(jVert, i);
      }
      ListOffDiagNorm[scal] += 1;
    }
  }
  auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
    seed ^= new_hash + 0x9e3779b8 + (seed << 6) + (seed >> 2);
  };
  size_t hash = seed;
  auto update_from_map=[&](std::map<T, size_t> const& map) -> void {
    for (auto & kv : map) {
      size_t hash1 = std::hash<T>()(kv.first);
      size_t hash2 = kv.second;
      combine_hash(hash, hash1);
      combine_hash(hash, hash2);
    }
  };
  update_from_map(ListDiagNorm);
  update_from_map(ListOffDiagNorm);
  return hash;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>> EnumerationDelaunayPolytopes(boost::mpi::communicator &comm,
                                                         DataBank<PolyhedralEntry<T, Tgroup>> &TheBank,
                                                         DataLattice<T, Tint> const &eData,
                                                         PolyHeuristic<mpz_class> const &AllArr) {
  int i_rank = comm.rank();
  int n_proc = comm.size();
  std::string FileLog = "log_" + std::to_string(n_proc) + "_" + std::sto_string(i_rank);
  std::ofstream os(FileLog);
  using Tobj = MyMatrix<Tint>;
  auto f_init=[&]() -> Tobj {
    Tobj EXT = FindDelaunayPolytope<T, Tint>(
       eData.GramMat, eData.CVPmethod, MProc.GetO(TheId));
    os << "Creation of a Delaunay with |V|=" << EXT.rows() << " vertices\n";
  };
  auto f_adj=[&](Tobj const& x) -> std::vector<Tobj> {
    Tgroup GRPlatt = Delaunay_Stabilizer<T, Tint, Tgroup>(eData, eDEL.EXT, eInv.eIndex);
    vectface TheOutput = DualDescriptionStandard(x, GRPlatt);
    std::vector<Tobj> l_obj;
    for (auto &eOrbB : TheOutput) {
      MyMatrix<Tint> EXTadj = FindAdjacentDelaunayPolytope<T, Tint>(eData.GramMat, EXT_T, eOrbB, eData.CVPmethod);
      l_obj.push_back(EXTadj);
    }
    return l_obj;
  };
  auto f_equiv=[&](Tobj const& x, Tobj const& y) -> bool {
    Tint PreIndex = Int_IndexLattice(x.EXT);
    Tint eIndex = T_abs(PreIndex);
    return Delaunay_TestEquivalence<T, Tint, Tgroup>(eData, x, y, eIndex);
  };
  std::vector<Tobj> l_obj;
  std::vector<int> l_status;
  auto f_hash=[&](size_t const& seed, Tobj const& x) -> size_t {
    return ComputeInvariantDelaunay(eData, seed, x);
  };
  auto f_exists=[&](int const& n_obj) -> bool {
    return false;
  };
  auto f_insert=[&](Tobj const& x) -> void {
    l_obj.push(x);
  };
  auto f_load=[&](size_t const& pos) -> Tobj {
    return l_obj[pos];
  };
  auto f_save_status=[&](size_t const& pos, bool const& val) -> void {
    int val_i = static_cast<int>(val);
    if (l_status.size() <= pos) {
      l_status
    } else {
    }
  };
  auto f_load_status=[&](size_t const& pos) -> bool {
    static_cast<bool>(l_status[pos]);
  };
  return l_obj;
}

FullNamelist NAMELIST_GetStandard_COMPUTE_DELAUNAY() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["GRAMfile"] = "unset.gram";
  ListStringValues1["SVRfile"] = "unset.svr";
  ListStringValues1["OUTfile"] = "unset.out";
  ListBoolValues1["SavingDelaunay"] = false;
  ListBoolValues1["FullDataInMemory"] = true;
  ListBoolValues1["ReturnAll"] = false;
  ListStringValues1["PrefixDelaunay"] = "/irrelevant/";
  SingleBlock BlockDATA;
  BlockDATA.ListIntValues = ListIntValues1;
  BlockDATA.ListBoolValues = ListBoolValues1;
  BlockDATA.ListDoubleValues = ListDoubleValues1;
  BlockDATA.ListStringValues = ListStringValues1;
  BlockDATA.ListListStringValues = ListListStringValues1;
  ListBlock["DATA"] = BlockDATA;
  // METHOD
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListIntValues2["NPROC"] = 1;
  ListIntValues2["UpperLimitMethod4"] = 10000;
  ListBoolValues2["Saving"] = false;
  ListStringValues2["PrefixPolyhedral"] = "/irrelevant/";
  ListBoolValues2["FullDataInMemory"] = true;
  ListStringValues2["SplittingHeuristicFile"] = "unset.heu";
  ListStringValues2["AdditionalSymmetryHeuristicFile"] = "unset.heu";
  ListStringValues2["DualDescriptionHeuristicFile"] = "unset.heu";
  ListStringValues2["StabEquivFacetHeuristicFile"] = "unset.heu";
  ListStringValues2["MethodInitialFacetSetFile"] = "unset.heu";
  ListStringValues2["CVPmethod"] = "SVexact";
  SingleBlock BlockMETHOD;
  BlockMETHOD.ListIntValues = ListIntValues2;
  BlockMETHOD.ListBoolValues = ListBoolValues2;
  BlockMETHOD.ListDoubleValues = ListDoubleValues2;
  BlockMETHOD.ListStringValues = ListStringValues2;
  BlockMETHOD.ListListStringValues = ListListStringValues2;
  ListBlock["METHOD"] = BlockMETHOD;
  // BANK
  std::map<std::string, int> ListIntValues3;
  std::map<std::string, bool> ListBoolValues3;
  std::map<std::string, double> ListDoubleValues3;
  std::map<std::string, std::string> ListStringValues3;
  std::map<std::string, std::vector<std::string>> ListListStringValues3;
  ListStringValues3["BankSaveHeuristicFile"] = "unset.heu";
  ListStringValues3["Prefix"] = "./unset/";
  ListBoolValues3["Saving"] = false;
  ListBoolValues3["FullDataInMemory"] = true;
  SingleBlock BlockBANK;
  BlockBANK.ListIntValues = ListIntValues3;
  BlockBANK.ListBoolValues = ListBoolValues3;
  BlockBANK.ListDoubleValues = ListDoubleValues3;
  BlockBANK.ListStringValues = ListStringValues3;
  BlockBANK.ListListStringValues = ListListStringValues3;
  ListBlock["BANK"] = BlockBANK;
  // Merging all data
  return {ListBlock, "undefined"};
}

template <typename T, typename Tint, typename Tgroup>
void ComputeDelaunayPolytope(FullNamelist const &eFull) {
  SingleBlock BlockBANK = eFull.ListBlock.at("BANK");
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  SingleBlock BlockMETHOD = eFull.ListBlock.at("METHOD");
  //
  bool BANK_Saving = BlockBANK.ListBoolValues.at("Saving");
  bool BANK_Memory = BlockBANK.ListBoolValues.at("FullDataInMemory");
  std::string BANK_Prefix = BlockBANK.ListStringValues.at("Prefix");
  CreateDirectory(BANK_Prefix);
  FctsDataBank<PolyhedralEntry<T, Tgroup>> recFct =
      GetRec_FctsDataBank<T, Tgroup>();
  DataBank<PolyhedralEntry<T, Tgroup>> TheBank(BANK_Saving, BANK_Memory,
                                               BANK_Prefix, recFct);
  //
  std::cerr << "Reading DATA\n";
  std::string GRAMfile = BlockDATA.ListStringValues.at("GRAMfile");
  MyMatrix<T> GramMat = ReadMatrixFile<T>(GRAMfile);
  //
  std::string SVRfile = BlockDATA.ListStringValues.at("SVRfile");
  MyMatrix<Tint> SVR = ReadMatrixFile<Tint>(SVRfile);
  std::cerr << "|SVR|=" << SVR.rows() << "\n";
  //
  std::string OUTfile = BlockDATA.ListStringValues.at("OUTfile");
  std::cerr << "OUTfile=" << OUTfile << "\n";

  bool Delaunay_Saving = BlockDATA.ListBoolValues.at("SavingDelaunay");
  bool Delaunay_Memory = BlockDATA.ListBoolValues.at("FullDataInMemory");
  bool Delaunay_ReturnAll = BlockDATA.ListBoolValues.at("ReturnAll");
  std::string PrefixDelaunay = BlockDATA.ListStringValues.at("PrefixDelaunay");
  CreateDirectory(PrefixDelaunay);
  std::string PrefixPolyhedral =
      BlockMETHOD.ListStringValues.at("PrefixPolyhedral");
  CreateDirectory(PrefixPolyhedral);
  int UppLimit = BlockMETHOD.ListIntValues.at("UpperLimitMethod4");
  std::string CVPmethod = BlockMETHOD.ListStringValues.at("CVPmethod");
  mpz_class UppLimit_z = UppLimit;
  int n = GramMat.rows();
  DataLattice<T, Tint> eData{n,
                             GramMat,
                             SVR,
                             Delaunay_Saving,
                             Delaunay_Memory,
                             PrefixDelaunay,
                             PrefixPolyhedral,
                             Delaunay_ReturnAll,
                             UppLimit_z,
                             CVPmethod};
  //
  std::cerr << "Creating MPROC\n";
  int NbThr = BlockMETHOD.ListIntValues.at("NPROC");
  MainProcessor MProc(NbThr);
  int TheId = MProc.MPU_GetId();
  //
  PolyHeuristic<mpz_class> AllArr = AllStandardHeuristic<mpz_class>();
  //
  std::string HeuSplitFile =
      BlockMETHOD.ListStringValues.at("SplittingHeuristicFile");
  if (HeuSplitFile != "unset.heu") {
    IsExistingFileDie(HeuSplitFile);
    std::ifstream SPLITfs(HeuSplitFile);
    AllArr.Splitting = ReadHeuristic<mpz_class>(SPLITfs);
  }
  //
  std::string AddiSymmFile =
      BlockMETHOD.ListStringValues.at("AdditionalSymmetryHeuristicFile");
  if (AddiSymmFile != "unset.heu") {
    IsExistingFileDie(AddiSymmFile);
    std::ifstream SYMMfs(AddiSymmFile);
    AllArr.AdditionalSymmetry = ReadHeuristic<mpz_class>(SYMMfs);
  }
  //
  std::string DualDescFile =
      BlockMETHOD.ListStringValues.at("DualDescriptionHeuristicFile");
  if (DualDescFile != "unset.heu") {
    IsExistingFileDie(DualDescFile);
    std::ifstream DDfs(DualDescFile);
    AllArr.DualDescriptionProgram = ReadHeuristic<mpz_class>(DDfs);
  }
  //
  std::string StabEquivFile =
      BlockMETHOD.ListStringValues.at("StabEquivFacetHeuristicFile");
  if (StabEquivFile != "unset.heu") {
    IsExistingFileDie(StabEquivFile);
    std::ifstream STEQfs(StabEquivFile);
    AllArr.StabEquivFacet = ReadHeuristic<mpz_class>(STEQfs);
  }
  //
  std::string MethodFacetFile =
      BlockMETHOD.ListStringValues.at("MethodInitialFacetSetFile");
  if (MethodFacetFile != "unset.heu") {
    IsExistingFileDie(MethodFacetFile);
    std::ifstream MIFSfs(MethodFacetFile);
    AllArr.InitialFacetSet = ReadHeuristic<mpz_class>(MIFSfs);
  }
  //
  std::string BankSaveFile =
      BlockBANK.ListStringValues.at("BankSaveHeuristicFile");
  if (BankSaveFile != "unset.heu") {
    IsExistingFileDie(BankSaveFile);
    std::ifstream BANKfs(BankSaveFile);
    AllArr.BankSave = ReadHeuristic<mpz_class>(BANKfs);
  }
  //
  bool DD_Saving = BlockMETHOD.ListBoolValues.at("Saving");
  bool DD_Memory = BlockMETHOD.ListBoolValues.at("FullDataInMemory");
  AllArr.Saving = DD_Saving;
  AllArr.eMemory = DD_Memory;
  //
  std::vector<Delaunay<T, Tint>> ListDel =
      EnumerationDelaunayPolytopes<T>(MProc, TheId, TheBank, eData, AllArr);
  std::cerr << "We now have ListDel\n";
  //
  if (Delaunay_ReturnAll) {
    std::ofstream OUTfs(OUTfile);
    int nbDel = ListDel.size();
    OUTfs << "nbDel=" << nbDel << "\n";
    for (int iDel = 0; iDel < nbDel; iDel++) {
      OUTfs << "iDel=" << iDel << "/" << nbDel << "\n";
      WriteMatrix(OUTfs, ListDel[iDel].EXT);
    }
  }
}

// clang-format off
#endif  // SRC_LATT_LATTICEDELAUNAY_H_
// clang-format on

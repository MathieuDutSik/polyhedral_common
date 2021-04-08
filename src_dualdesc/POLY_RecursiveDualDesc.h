#ifndef INCLUDE_POLY_RECURSIVE_DUAL_DESC_H
#define INCLUDE_POLY_RECURSIVE_DUAL_DESC_H

#include "Namelist.h"
#include "POLY_Heuristics.h"
#include "POLY_DirectDualDesc.h"
#include "POLY_SamplingFacet.h"
#include "Temp_PolytopeEquiStab.h"
#include "GRP_GroupFct.h"

#include "POLY_GAP.h"
#include "POLY_netcdf.h"



template<typename T, typename Tgroup>
struct EquivariantDualDescription {
  MyMatrix<T> EXT;
  Tgroup GRP;
  std::vector<Face> ListFace;
};



template<typename T, typename Tgroup>
EquivariantDualDescription<T,Tgroup> ConvertGAPread_EquivDualDesc(datagap::DataGAP<T,typename Tgroup::Telt> const& dataEXT, datagap::DataGAP<T,typename Tgroup::Telt> const& dataFAC)
{
  if (dataEXT.Nature != datagap::int_record) {
    std::cerr << "For EquivDualDesc, we need to have a record as entry\n";
    throw TerminalException{1};
  }
  int pos_EXT = -1;
  int pos_GRP = -1;
  int n_pos = dataEXT.ListRec.size();
  for (int pos=0; pos<n_pos; pos++) {
    if (dataEXT.ListRec[pos].first == "EXT")
      pos_EXT = pos;
    if (dataEXT.ListRec[pos].first == "Group")
      pos_GRP = pos;
  }
  if (pos_EXT == -1) {
    std::cerr << "Failed to find entry EXT in the record\n";
    throw TerminalException{1};
  }
  if (pos_GRP == -1) {
    std::cerr << "Failed to find entry Group in the record\n";
    throw TerminalException{1};
  }
  MyMatrix<T> EXT = datagap::ConvertGAPread_MyMatrixT(dataEXT.ListRec[pos_EXT].second);
  int n_rows = EXT.rows();
  Tgroup GRP = datagap::ConvertGAPread_PermutationGroup<T, Tgroup>(dataEXT.ListRec[pos_GRP].second, n_rows);
  //
  std::vector<Face> ListFace = ConvertGAPread_ListFace(dataFAC, n_rows);
  //
  return {std::move(EXT), std::move(GRP), std::move(ListFace)};
}




template<typename T, typename Tgroup>
void Write_EquivDualDesc(EquivariantDualDescription<T,Tgroup> const& eRec, std::string const& eFile)
{
  if (!FILE_IsFileMakeable(eFile)) {
    std::cerr << "Error in Write_EquivDualDesc: File eFile=" << eFile << " is not makeable\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  POLY_NC_WritePolytope(dataFile, eRec.EXT);
  bool orbit_setup = true;
  bool orbit_status = false;
  POLY_NC_WriteGroup(dataFile, eRec.GRP, orbit_setup, orbit_status);
  //
  size_t n_orbit = eRec.ListFace.size();
  for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++)
    POLY_NC_WriteFace(dataFile, i_orbit, eRec.ListFace[i_orbit]);
}



template<typename T>
std::pair<MyMatrix<T>, std::vector<Face>> Read_BankEntry(std::string const& eFile)
{
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  MyMatrix<T> EXT = POLY_NC_ReadPolytope<T>(dataFile);
  //
  std::vector<Face> ListFace = POLY_NC_ReadAllFaces(dataFile);
  return {std::move(EXT), std::move(ListFace)};
}

template<typename T>
void Write_BankEntry(std::string const& eFile, MyMatrix<T> const& EXT, std::vector<Face> const& ListFace)
{
  if (!FILE_IsFileMakeable(eFile)) {
    std::cerr << "Error in Write_BankEntry: File eFile=" << eFile << " is not makeable\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  POLY_NC_WritePolytope(dataFile, EXT);
  int n_act = EXT.rows();
  POLY_NC_WriteOrbitDimVars(dataFile, n_act);
  //
  size_t n_orbit = ListFace.size();
  for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++)
    POLY_NC_WriteFace(dataFile, i_orbit, ListFace[i_orbit]);
}



template<typename T, typename Tidx>
std::pair<MyMatrix<T>, std::vector<Tidx>> CanonicalizationPolytopePair(MyMatrix<T> const& EXT, WeightMatrix<T,T> const& WMat)
{
  std::pair<std::vector<Tidx>, std::vector<Tidx>> PairCanonic = GetCanonicalizationVector<T,T,GraphBitset,Tidx>(WMat);
  int n_row=EXT.rows();
  int n_col=EXT.cols();
  MyMatrix<T> EXTcan(n_row, n_col);
  for (int i_row=0; i_row<n_row; i_row++) {
    int j_row=PairCanonic.second[i_row];
    EXTcan.row(i_row) = EXT.row(j_row);
  }
  MyMatrix<T> RowRed = RowReduction(EXTcan);
  MyMatrix<T> EXTret = EXTcan * Inverse(RowRed);
  MyMatrix<T> EXTretB = RemoveFractionMatrix(EXTret);
  //
  return {std::move(EXTretB), std::move(PairCanonic.second)};
}



template<typename T>
MyMatrix<T> CanonicalizationPolytope(MyMatrix<T> const& EXT)
{
  WeightMatrix<T, T> WMat=GetWeightMatrix(EXT);
  ReorderingSetWeight(WMat);
  return CanonicalizationPolytopePair<T,int>(EXT, WMat).first;
}



template<typename T, typename Telt>
struct DataBank {
private:
  using Tidx = typename Telt::Tidx;
  int MinSize;
  std::unordered_map<MyMatrix<T>, std::vector<Face>> ListEnt;
  bool Saving;
  std::string SavingPrefix;
public:
  DataBank(bool const& _Saving, std::string const& _SavingPrefix) : Saving(_Saving), SavingPrefix(_SavingPrefix)
  {
    MinSize = std::numeric_limits<int>::max();
    if (Saving) {
      size_t iOrbit=0;
      while(true) {
        std::string eFileBank = SavingPrefix + "DualDesc" + std::to_string(iOrbit) + ".nc";
        if (!IsExistingFile(eFileBank))
          break;
        std::pair<MyMatrix<T>, std::vector<Face>> ePair = Read_BankEntry<T>(eFileBank);
        int e_size = ePair.first.rows();
        std::cerr << "Read iOrbit=" << iOrbit << " FileBank=" << eFileBank << " |EXT|=" << e_size << " |ListFace|=" << ePair.second.size() << "\n";
        ListEnt[ePair.first] = ePair.second;
        MinSize = std::min(MinSize, e_size);
        iOrbit++;
      }
    }
  }
  void InsertEntry(MyMatrix<T> const& EXT, WeightMatrix<T,T> const& WMat, std::vector<Face> const& ListFace)
  {
    std::pair<MyMatrix<T>, std::vector<Tidx>> ePair = CanonicalizationPolytopePair<T,Tidx>(EXT, WMat);
    std::vector<Face> ListFaceO;
    Telt perm1 = Telt(ePair.second);
    Telt ePerm = ~perm1;
    for (auto & eFace : ListFace) {
      Face eInc = OnFace(eFace, ePerm);
      ListFaceO.push_back(eInc);
    }
    if (Saving) {
      size_t n_orbit = ListEnt.size();
      std::string eFile = SavingPrefix + "DualDesc" + std::to_string(n_orbit) + ".nc";
      std::cerr << "Insert entry to file eFile=" << eFile << "\n";
      Write_BankEntry(eFile, ePair.first, ListFaceO);
    }
    ListEnt[ePair.first] = ListFaceO;
    int e_size = ePair.first.rows();
    MinSize = std::min(MinSize, e_size);
  }
  std::vector<Face> GetDualDesc(MyMatrix<T> const& EXT, WeightMatrix<T,T> const& WMat) const
  {
    std::cerr << "Passing by GetDualDesc |ListEnt|=" << ListEnt.size() << "\n";
    std::pair<MyMatrix<T>, std::vector<Tidx>> ePair = CanonicalizationPolytopePair<T,Tidx>(EXT, WMat);
    auto iter = ListEnt.find(ePair.first);
    std::vector<Face> ListReprTrans;
    if (iter != ListEnt.end()) {
      std::cerr << "Finding a matching entry\n";
      Telt ePerm = Telt(ePair.second);
      for (auto const& eOrbit : iter->second) {
	Face eListJ=OnFace(eOrbit, ePair.second);
	ListReprTrans.push_back(eListJ);
      }
    }
    return ListReprTrans; // If returning empty then it means nothing has been found.
  }
  int get_minsize() const
  {
    return MinSize;
  }
};



template<typename T, typename Tint, typename Tgroup>
struct DatabaseOrbits {
private:
  using Telt = typename Tgroup::Telt;
  MyMatrix<T> EXT;
  Tgroup GRP;
  Tint groupOrder;
  std::string eFile;
  netCDF::NcFile dataFile;
  bool SavingTrigger;
  struct SingEnt {
    size_t pos;
    Tint orbSize;
  };
  std::unordered_map<Face, SingEnt> DictOrbit;
  std::map<size_t, std::unordered_set<Face>> CompleteList_SetUndone;
  std::vector<Face> ListOrbit;
  Tint TotalNumber;
  size_t nbOrbitDone;
  Tint nbUndone;
  size_t nbOrbit;
  size_t n_grpsize;

public:
  void InsertEntryDatabase(Face const& face, bool const& status, Tint const& orbSize, size_t const& pos)
  {
    //    std::cerr << "status=" << status << " orbSize=" << orbSize << " pos=" << pos << "\n";
    DictOrbit[face] = {pos, orbSize};
    if (!status) {
      size_t len = face.count();
      CompleteList_SetUndone[len].insert(face);
    }
    ListOrbit.push_back(face);
    TotalNumber += orbSize;
    if (status) {
      nbOrbitDone++;
    } else {
      nbUndone += orbSize;
    }
    nbOrbit++;
  }
  DatabaseOrbits(MyMatrix<T> const& _EXT, Tgroup const& _GRP, std::string const& _eFile, bool const& _SavingTrigger) : EXT(_EXT), GRP(_GRP), eFile(_eFile), SavingTrigger(_SavingTrigger)
  {
    TotalNumber = 0;
    nbOrbitDone = 0;
    nbUndone = 0;
    nbOrbit = 0;
    std::vector<Telt> LGen = GRP.GeneratorsOfGroup();
    groupOrder = GRP.size();
    if (SavingTrigger) {
      std::cerr << "eFile=" << eFile << "\n";
      if (IsExistingFile(eFile)) {
        dataFile.open(eFile, netCDF::NcFile::write);
      } else {
        if (!FILE_IsFileMakeable(eFile)) {
          std::cerr << "Error in DatabaseOrbits: File eFile=" << eFile << " is not makeable\n";
          throw TerminalException{1};
        }
        dataFile.open(eFile, netCDF::NcFile::replace, netCDF::NcFile::nc4);
        POLY_NC_WritePolytope(dataFile, EXT);
        bool orbit_setup = true;
        bool orbit_status = true;
        POLY_NC_WriteGroup(dataFile, GRP, orbit_setup, orbit_status);
      }
      netCDF::NcDim eDim = dataFile.getDim("n_orbit");
      size_t n_orbit = eDim.getSize();
      netCDF::NcDim fDim = dataFile.getDim("n_grpsize");
      n_grpsize = fDim.getSize();
      //
      for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++) {
        SingleEntryStatus<Tint> eEnt = POLY_NC_ReadSingleEntryStatus<Tint>(dataFile, i_orbit);
        InsertEntryDatabase(eEnt.face, eEnt.status, eEnt.OrbSize, i_orbit);
      }
    }
  }
  ~DatabaseOrbits()
  {
  }
  void FuncInsert(Face const& face)
  {
    Face face_can = GRP.CanonicalImage(face);
    if (DictOrbit.count(face_can) == 1)
      return;
    //
    Tint ordStab = GRP.Stabilizer_OnSets(face).size();
    Tint orbSize = groupOrder / ordStab;
    InsertEntryDatabase(face_can, false, orbSize, nbOrbit);
    //
    if (SavingTrigger) {
      SingleEntryStatus<Tint> eEnt{false, face_can, orbSize};
      POLY_NC_WriteSingleEntryStatus(dataFile, nbOrbit - 1, eEnt, n_grpsize);
    }
  }
  void FuncPutOrbitAsDone(size_t const& iOrb)
  {
    Face face = ListOrbit[iOrb];
    Tint orbSize = DictOrbit[face].orbSize;
    if (SavingTrigger) {
      POLY_NC_SetBit(dataFile, iOrb, true);
    }
    size_t len = face.count();
    CompleteList_SetUndone[len].erase(face);
    nbUndone -= orbSize;
    nbOrbitDone++;
  }
  Face ComputeIntersectionUndone() const
  {
    size_t n_row = EXT.rows();
    Face eSetReturn(n_row);
    for (size_t i_row=0; i_row<n_row; i_row++)
      eSetReturn[i_row] = 1;
    for (auto & eEnt : CompleteList_SetUndone) {
      for (auto & eFace : eEnt.second) {
        eSetReturn &= OrbitIntersection(GRP, eFace);
        if (eSetReturn.count() == 0)
          return eSetReturn;
      }
    }
    return eSetReturn;
  }
  std::vector<Face> FuncListOrbitIncidence()
  {
    if (SavingTrigger) {
      RemoveFile(eFile);
    }
    return ListOrbit;
  }
  Tint FuncNumber() const
  {
    return TotalNumber;
  }
  Tint FuncNumberUndone() const
  {
    return nbUndone;
  }
  size_t FuncNumberOrbit() const
  {
    return nbOrbit;
  }
  size_t FuncNumberOrbitDone() const
  {
    return nbOrbitDone;
  }
  std::pair<size_t,Face> FuncGetMinimalUndoneOrbit()
  {
    for (auto & eEnt : CompleteList_SetUndone) {
      if (eEnt.second.size() > 0) {
        auto iter = eEnt.second.begin();
        Face face = *iter;
        return {DictOrbit[face].pos, face};
      }
    }
    return {-1,{}};
  }
};





//
// A number of appoximations are done in this code:
// ---In the bank we assume that the full symmetry is used.
//    This means less things to store
// ---We use the canonicalization approach which allows to treat smaller cases.
// ---Serial mode. Should be faster indeed.
// ---
//
template<typename T,typename Tgroup>
std::vector<Face> DUALDESC_AdjacencyDecomposition(
         DataBank<T,typename Tgroup::Telt> & TheBank,
	 MyMatrix<T> const& EXT,
	 Tgroup const& GRP,
	 PolyHeuristicSerial<typename Tgroup::Tint> const& AllArr,
	 std::string const& ePrefix,
	 int const& TheLevel)
{
  using Tint=typename Tgroup::Tint;
  MyMatrix<T> EXTred=ColumnReduction(EXT);
  int nbRow=EXTred.rows();
  int eRank=EXTred.cols();
  WeightMatrix<T, T> WMat;
  bool HaveWMat=false;
  auto ComputeWMat=[&]() -> void {
    if (HaveWMat)
      return;
    WMat=GetWeightMatrix(EXTred);
    ReorderingSetWeight(WMat);
  };
  //
  // Checking if the entry is present in the map.
  //
  if (nbRow >= TheBank.get_minsize()) {
    ComputeWMat();
    std::vector<Face> ListFace = TheBank.GetDualDesc(EXTred, WMat);
    if (ListFace.size() > 0)
      return ListFace;
  }
  //
  // Now computing the groups
  //
  Tgroup TheGRPrelevant;
  std::map<std::string, Tint> TheMap;
  int nbVert=EXTred.rows();
  int delta=nbVert - eRank;
  TheMap["groupsize"]=GRP.size();
  TheMap["incidence"]=nbVert;
  TheMap["level"]=TheLevel;
  TheMap["rank"]=eRank;
  TheMap["delta"]=delta;
  std::string ansSplit=HeuristicEvaluation(TheMap, AllArr.Splitting);
  //
  // The computations themselves
  //
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  std::vector<Face> ListOrbitFaces;
  std::string ansSymm;
  if (ansSplit != "split") {
    TheGRPrelevant = GRP;
    std::string ansProg=HeuristicEvaluation(TheMap, AllArr.DualDescriptionProgram);
    ListOrbitFaces = DirectFacetOrbitComputation(EXTred, GRP, ansProg, std::cerr);
    ansSymm="no";
  } else {
    ansSymm = HeuristicEvaluation(TheMap, AllArr.AdditionalSymmetry);
    std::cerr << "ansSymm=" << ansSymm << "\n";
    if (ansSymm == "yes") {
      ComputeWMat();
      TheGRPrelevant = GetStabilizerWeightMatrix<T,T,Tgroup>(WMat);
    } else {
      TheGRPrelevant = GRP;
    }
    if (TheGRPrelevant.size() == GRP.size())
      ansSymm="no";
    TheMap["groupsizerelevant"]=TheGRPrelevant.size();
    std::string eFile = ePrefix + "Database_" + std::to_string(TheLevel) + "_" + std::to_string(nbVert) + "_" + std::to_string(eRank) + ".nc";
    bool SavingTrigger=AllArr.Saving;
    DatabaseOrbits<T,Tint,Tgroup> RPL(EXTred, GRP, eFile, SavingTrigger);
    int NewLevel = TheLevel + 1;
    if (RPL.FuncNumberOrbit() == 0) {
      std::string ansSamp=HeuristicEvaluation(TheMap, AllArr.InitialFacetSet);
      std::vector<Face> ListFace=DirectComputationInitialFacetSet(EXTred, ansSamp);
      std::cerr << "After DirectComputationInitialFacetSet |ListFace|=" << ListFace.size() << "\n";
      for (auto & eInc : ListFace)
	RPL.FuncInsert(eInc);
    }
    Tint TheDim = eRank-1;
    while(true) {
      Tint nbUndone=RPL.FuncNumberUndone();
      if (RPL.FuncNumberOrbitDone() > 0) {
        Face eSetUndone=RPL.ComputeIntersectionUndone();
        if (nbUndone <= TheDim-1 || eSetUndone.count() > 0) {
          std::cerr << "End of computation, nbObj=" << RPL.FuncNumber() << " nbUndone=" << nbUndone << " |eSetUndone|=" << eSetUndone.count() << " Depth=" << TheLevel << " |EXT|=" << nbRow << "\n";
          break;
        }
      }
      std::pair<size_t,Face> ePair = RPL.FuncGetMinimalUndoneOrbit();
      size_t SelectedOrbit = ePair.first;
      Face eInc = ePair.second;
      MyMatrix<T> EXTredFace=SelectRow(EXTred, eInc);
      FlippingFramework<T> FF(EXTred, eInc);
      Tgroup TheStab=TheGRPrelevant.Stabilizer_OnSets(eInc);
      Tint OrbSize=TheGRPrelevant.size() / TheStab.size();
      Tgroup GRPred=ReducedGroupAction(TheStab, eInc);
      std::cerr << "Considering orbit " << SelectedOrbit << " |EXT|=" << eInc.size() << " |inc|=" << eInc.count() << " Level=" << TheLevel << " |stab|=" << GRPred.size() << " dim=" << TheDim << "\n";
      std::string eDir = ePrefix + "ADM" + std::to_string(SelectedOrbit) + "_";
      std::vector<Face> TheOutput=DUALDESC_AdjacencyDecomposition(TheBank, EXTredFace, GRPred, AllArr, eDir, NewLevel);
      for (auto& eOrbB : TheOutput) {
        Face eFlip=FF.Flip(eOrbB);
        RPL.FuncInsert(eFlip);
      }
      RPL.FuncPutOrbitAsDone(SelectedOrbit);
      nbUndone = RPL.FuncNumberUndone();
      std::cerr << "We have " << RPL.FuncNumberOrbit() << " orbits  Nr treated=" << RPL.FuncNumberOrbitDone() << " orbits  nbUndone=" << nbUndone << "\n";
      std::cerr << "\n";
    };
    ListOrbitFaces = RPL.FuncListOrbitIncidence();
  }
  end = std::chrono::system_clock::now();
  int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
  TheMap["time"]=elapsed_seconds;
  std::string ansBank=HeuristicEvaluation(TheMap, AllArr.BankSave);
  std::cerr << "elapsed_seconds=" << elapsed_seconds << " ansBank=" << ansBank << " ansSymm=" << ansSymm << "\n";
  if (ansBank == "yes") {
    ComputeWMat();
    TheBank.InsertEntry(EXTred, WMat, ListOrbitFaces);
  }
  if (ansSymm == "yes") {
    return OrbitSplittingListOrbit(TheGRPrelevant, GRP, ListOrbitFaces, std::cerr);
  } else {
    return ListOrbitFaces;
  }
}



FullNamelist NAMELIST_GetStandard_RecursiveDualDescription()
{
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, std::string> ListStringValues1;
  ListStringValues1["EXTfile"]="unset.ext";
  ListStringValues1["GRPfile"]="unset.grp";
  ListStringValues1["OUTfile"]="unset.out";
  SingleBlock BlockDATA;
  BlockDATA.ListStringValues=ListStringValues1;
  ListBlock["DATA"]=BlockDATA;
  // HEURISTIC
  std::map<std::string, std::string> ListStringValuesH;
  ListStringValuesH["SplittingHeuristicFile"]="unset.heu";
  ListStringValuesH["AdditionalSymmetryHeuristicFile"]="unset.heu";
  ListStringValuesH["DualDescriptionHeuristicFile"]="unset.heu";
  ListStringValuesH["MethodInitialFacetSetFile"]="unset.heu";
  ListStringValuesH["BankSaveHeuristicFile"]="unset.heu";
  SingleBlock BlockHEURIS;
  BlockHEURIS.ListStringValues=ListStringValuesH;
  ListBlock["HEURISTIC"]=BlockHEURIS;
  // METHOD
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  ListBoolValues2["Saving"]=false;
  ListStringValues2["Prefix"]="/irrelevant/";
  SingleBlock BlockMETHOD;
  BlockMETHOD.ListIntValues=ListIntValues2;
  BlockMETHOD.ListBoolValues=ListBoolValues2;
  BlockMETHOD.ListDoubleValues=ListDoubleValues2;
  BlockMETHOD.ListStringValues=ListStringValues2;
  ListBlock["METHOD"]=BlockMETHOD;
  // BANK
  std::map<std::string, int> ListIntValues3;
  std::map<std::string, bool> ListBoolValues3;
  std::map<std::string, double> ListDoubleValues3;
  std::map<std::string, std::string> ListStringValues3;
  ListStringValues3["Prefix"]="./unset/";
  ListBoolValues3["Saving"]=false;
  SingleBlock BlockBANK;
  BlockBANK.ListBoolValues=ListBoolValues3;
  BlockBANK.ListStringValues=ListStringValues3;
  ListBlock["BANK"]=BlockBANK;
  // Merging all data
  return {std::move(ListBlock), "undefined"};
}



template<typename T, typename Tgroup>
void MainFunctionSerialDualDesc(FullNamelist const& eFull)
{
  using Telt=typename Tgroup::Telt;
  using Tint=typename Tgroup::Tint;
  SingleBlock BlockBANK=eFull.ListBlock.at("BANK");
  bool BANK_IsSaving=BlockBANK.ListBoolValues.at("Saving");
  std::string BANK_Prefix=BlockBANK.ListStringValues.at("Prefix");
  DataBank<T,Telt> TheBank(BANK_IsSaving, BANK_Prefix);
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
  SingleBlock BlockMETHOD=eFull.ListBlock.at("METHOD");
  //
  PolyHeuristicSerial<Tint> AllArr=AllStandardHeuristicSerial<Tint>();
  //
  SetHeuristic(eFull, "SplittingHeuristicFile", AllArr.Splitting);
  SetHeuristic(eFull, "AdditionalSymmetryHeuristicFile", AllArr.AdditionalSymmetry);
  SetHeuristic(eFull, "DualDescriptionHeuristicFile", AllArr.DualDescriptionProgram);
  SetHeuristic(eFull, "MethodInitialFacetSetFile", AllArr.InitialFacetSet);
  SetHeuristic(eFull, "BankSaveHeuristicFile", AllArr.BankSave);
  std::cerr << "SplittingHeuristicFile\n" << AllArr.Splitting << "\n";
  std::cerr << "AdditionalSymmetryHeuristicFile\n" << AllArr.AdditionalSymmetry << "\n";
  std::cerr << "DualDescriptionHeuristicFile\n" << AllArr.DualDescriptionProgram << "\n";
  std::cerr << "MethodInitialFacetSetFile\n" << AllArr.InitialFacetSet << "\n";
  std::cerr << "BankSaveHeuristicFile\n" << AllArr.BankSave << "\n";
  //
  bool DD_Saving=BlockMETHOD.ListBoolValues.at("Saving");
  std::string DD_Prefix=BlockMETHOD.ListStringValues.at("Prefix");
  AllArr.Saving=DD_Saving;
  //
  int TheLevel=0;
  std::vector<Face> TheOutput=DUALDESC_AdjacencyDecomposition(TheBank, EXT, GRP, AllArr, DD_Prefix, TheLevel);
  std::cerr << "|TheOutput|=" << TheOutput.size() << "\n";
  //
  std::ofstream OUTfs(OUTfile);
  VectVectInt_Magma_Print(OUTfs, TheOutput);
  //
}












#endif

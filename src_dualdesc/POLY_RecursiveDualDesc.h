#ifndef INCLUDE_POLY_RECURSIVE_DUAL_DESC_H
#define INCLUDE_POLY_RECURSIVE_DUAL_DESC_H

#include "POLY_ThreadDualDescription.h"
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
  return {EXT, GRP, ListFace};
}




template<typename T, typename Tgroup>
void Write_EquivDualDesc(EquivariantDualDescription<T,Tgroup> const& eRec, std::string const& eFile)
{
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



template<typename Tgroup>
Face OrbitIntersection(Tgroup const& GRP, Face const& gList)
{
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> LGen = GRP.GeneratorsAsGroup();
  int n = GRP.n_act();
  Face rList = gList;
  while(true) {
    size_t eSumPrev = rList.count();
    for (auto & eGen : LGen) {
      for (int i=0; i<n; i++) {
        std::size_t j = OnPoints(i, eGen);
        if (rList[i] == 0)
          rList[j]=0;
      }
    }
    size_t eSum = rList.count();
    if (eSum == eSumPrev)
      break;
  }
  return rList;
}


template<typename Tgroup>
Face OrbitUnion(Tgroup const& GRP, Face const& gList)
{
  int n = GRP.n_act();
  Face gListB(n);
  for (int i=0; i<n; i++)
    gListB[i] = 1 - gList[i];
  Face rListB=OrbitIntersection(GRP, gListB);
  for (int i=0; i<n; i++)
    rListB[i] = 1 - rListB[i];
  return rListB;
}







template<typename T, typename Tint, typename Tgroup>
struct DatabaseOrbits {
private:
  using Telt = typename Tgroup::Telt;
  Tgroup GRP;
  Tint groupOrder;
  MyMatrix<T> EXT;
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
  void InsertEntryDatavase(Face const& face, bool const& status, Tint const& orbSize, size_t const& pos)
  {
    DictOrbit[face] = {pos, orbSize};
    if (!status) {
      size_t len = face.count();
      CompleteList_SetUndone[len].insert(face);
    }
    ListOrbit.push_back(face);
    TotalNumber += orbSize;
    if (!status) {
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
    int n_act = EXT.rows();
    std::vector<Telt> LGen = GeneratorsOfGroup(GRP);
    groupOrder = GRP.size();
    if (SavingTrigger) {
      if (IsExistingFile(eFile)) {
        dataFile.open(eFile, netCDF::NcFile::write);
      } else {
        dataFile.open(eFile, netCDF::NcFile::replace, netCDF::NcFile::nc4);
        POLY_NC_WritePolytope(dataFile, EXT);
        bool orbit_setup = true;
        bool orbit_status = true;
        POLY_NC_WriteGroup(dataFile, GRP, orbit_setup, orbit_status);
      }
      netCDF::NcDim eDim = dataFile.getDim("n_orbit");
      size_t n_orbit = eDim.getSize();
      netCDF::NcDim fDim = dataFile.getDim("n_grpsize");
      size_t n_grpsize = fDim.getSize();
      //
      for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++) {
        SingleEntryStatus<Tint> eEnt = POLY_NC_ReadSingleEntryStatus<Tint>(dataFile, i_orbit);
        InsertEntryDatavase(eEnt.face, eEnt.status, eEnt.orbSize, i_orbit);
      }
    }
  }
  ~DatabaseOrbits()
  {
    if (SavingTrigger) {
      RemoveFile(eFile);
    }
  }
  void FuncInsert(Face const& face)
  {
    Face face_can = GRP.CanonicalImage(face);
    if (DictOrbit.count(face_can) == 1)
      return;
    //
    Tint ordStab = GRP.Stabilizer_OnSets(face).size();
    Tint orbSize = groupOrder / ordStab;
    InsertEntryDatabase(face_can, false, orbSize, nbOrbit+1);
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
    for (size_t i=0; i<n_row; i++)
      eSetReturn[i] = 1;
    for (auto & eEnt : CompleteList_SetUndone) {
      for (auto & eFace : eEnt.second) {
        eSetReturn = Intersection(eSetReturn, OrbitIntersection(GRP, eFace));
        if (eSetReturn.size() == 0)
          return eSetReturn;
      }
    }
    return eSetReturn;
  }
  std::vector<Face> FuncListOrbitIncidence()
  {
    return ListOrbit;
  }
  Face FuncRecord(size_t const& iOrb) const
  {
    return ListOrbit[iOrb];
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
  size_t FuncGetMinimalUndoneOrbit()
  {
    for (auto & eEnt : CompleteList_SetUndone) {
      if (eEnt.second.size() > 0) {
        auto iter = eEnt.second.begin();
        Face face = *iter;
        return DictOrbit[face].pos;
      }
    }
  }
};












#endif

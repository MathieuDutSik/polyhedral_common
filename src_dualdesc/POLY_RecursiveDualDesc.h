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
EquivariantDualDescription<T,Tgroup> ConvertGAPread_EquivDualDesc(datagap::DataGAP<T,Tgroup::Telt> const& dataEXT, datagap::DataGAP<T,Tgroup::Telt> const& dataFAC)
{
  if (dataEXT.Nature != datagap::int_record) {
    std::cerr << "For EquivDualDesc, we need to have a record as entry\n";
    throw TerminalException{1};
  }
  int pos_EXT = -1;
  int pos_GRP = -1;
  for (int pos=0; pos<dataEXT.ListRec.size(); pos++) {
    if (dataEXT.ListRec[pos].first == "EXT")
      pos_EXT = pos;
    if (dataEXT.ListRec[pos].first == "Group")
      pos_GRP = pos;
  }
  MyMatrix<T> EXT = ConvertGAPread_MyMatrixT(dataEXT.ListRec[pos_EXT].second);
  int n_rows = EXT.rows();
  Tgroup GRP = ConvertGAPread_PermutationGroup(dataEXT.ListRec[pos_GRP].second, n_rows);
  //
  std::vector<Face> ListFace = ConvertGAPread_ListFace(dataFAC);
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
  int n_orbit = eRec.ListFace.size();
  for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++)
    POLY_NC_WriteFace(dataFile, i_orbit, eRec.ListFace[i_orbit]);
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
      nbOrbitdone++;
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
      if (IsExisingFile(eFile)) {
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
  ~DatabaseOrbit()
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
      for (auto & eFace : eEnt) {
        eSetReturn = Intersection(eSetReturn, GetOrbitIntersection(eFace));
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
      if (eEnt.size() > 0) {
        auto iter = eEnt.begin();
        Face face = *iter;
        return DictOrbit[face].pos;
      }
    }
  }
};



#endif

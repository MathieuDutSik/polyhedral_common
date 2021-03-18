#ifndef INCLUDE_POLY_RECURSIVE_DUAL_DESC_H
#define INCLUDE_POLY_RECURSIVE_DUAL_DESC_H


template<typename T, typename Tgroup>
struct EquivariantDualDescription {
  MyMatrix<T> EXT;
  Tgroup GRP;
  std::vector<Face> ListFace;
};




template<typename T, typename Tgroup>
EquivariantDualDescription<T,Tgroup> ConvertGAPread_EquivDualDesc(datagap::DataGAP<T,Tgroup::Telt> const& data)
{
  if (data.Nature != datagap::int_record) {
    std::cerr << "For EquivDualDesc, we need to have a record as entry\n";
  }
  
}









template<typename T, typename Tint, typename Tgroup>
struct DatabaseOrbits {
private:
  using Telt = Tgroup::element_type;
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
  std::vector<int> ListOrbit;
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
  DatabaseOrbit(MyMatrix<T> const& _EXT, Tgroup const& _GRP, std::string const& _eFile, bool const& _SavingTrigger) : EXT(_EXT), GRP(_GRP), eFile(_eFile), SavingTrigger(_SavingTrigger)
  {
    TotalNumber = 0;
    nbOrbitDone = 0;
    nbUndone = 0;
    nbOrbit = 0;
    int n_act = EXT.rows();
    std::vector<Telt> LGen = GeneratorsOfGroup(GRP);
    groupOrder = Size<Tint>(GRP);
    GRPinfo eGRP{n_act, LGen, groupOrder};
    if (SavingTrigger) {
      if (IsExisingFile(eFile)) {
        dataFile.open(eFile, netCDF::NcFile::write);
      } else {
        dataFile.open(eFile, netCDF::NcFile::replace, netCDF::NcFile::nc4);
        POLY_NC_WritePolytope(dataFile, EXT);
        bool orbit_setup = true;
        bool orbit_status = true;
        POLY_NC_WriteGroup(dataFile, eGRP, orbit_setup, orbit_status);
      }
      netCDF::NcDim eDim = dataFile.getDim("n_orbit");
      size_t n_orbit = eDim.getSize();
      netCDF::NcDim fDim = dataFile.getDim("n_grpsize");
      size_t n_grpsize = fDim.getSize();
      //
      for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++) {
        SingleEntryStatus<Tint> eEnt = POLY_NC_ReadSingleEntryStatus(dataFile, i_orbit);
        InsertEntryDatavase(eEnt.face, eEnt.status, eEnt.orbSize, i_orbit);
      }
    }
  };
  ~DatabaseOrbit()
  {
    if (SavingTrigger) {
      RemoveFile(eFile);
    }
  }
  void FuncInsert(face)
  {
    Face face_can = CanonicalImage_OnSets(GRP, face);
    if (DictOrbit.count(face_can) == 1)
      return;
    //
    Tint ordStab = Size<Tint>(Stabilizer_OnSets(GRP, face));
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
}



#endif

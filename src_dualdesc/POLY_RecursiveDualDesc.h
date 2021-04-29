#ifndef INCLUDE_POLY_RECURSIVE_DUAL_DESC_H
#define INCLUDE_POLY_RECURSIVE_DUAL_DESC_H

#include "Namelist.h"
#include "POLY_Heuristics.h"
#include "POLY_DirectDualDesc.h"
#include "POLY_SamplingFacet.h"
#include "Temp_PolytopeEquiStab.h"
#include "GRP_GroupFct.h"

#include "POLY_GAP.h"
#include "POLY_netcdf_file.h"
#include "MatrixGroupBasic.h"
#include <signal.h>


//#define UNORDERED_MAP
#define TSL_SPARSE_MAP
//#define TSL_ROBIN_MAP
//#define TSL_HOPSCOTCH_MAP

#ifdef UNORDERED_MAP
# include <unordered_map>
# include <unordered_set>
# define UNORD_MAP std::unordered_map
# define UNORD_SET std::unordered_set
#endif

#ifdef TSL_SPARSE_MAP
# include "sparse_map.h"
# include "sparse_set.h"
# define UNORD_MAP tsl::sparse_map
# define UNORD_SET tsl::sparse_set
#endif

#ifdef TSL_ROBIN_MAP
# include "robin_map.h"
# include "robin_set.h"
# define UNORD_MAP tsl::robin_map
# define UNORD_SET tsl::robin_set
#endif

#ifdef TSL_HOPSCOTCH_MAP
# include "hopscotch_map.h"
# include "hopscotch_set.h"
# define UNORD_MAP tsl::hopscotch_map
# define UNORD_SET tsl::hopscotch_set
#endif




std::atomic<bool> ExitEvent;

void signal_callback_handler(int signum) {
  std::cout << "Caught signal " << signum << "\n";
  std::cout << "We are going to exit hopefully\n";
  ExitEvent = true;
}





template<typename T, typename Tgroup>
struct EquivariantDualDescription {
  MyMatrix<T> EXT;
  Tgroup GRP;
  vectface ListFace;
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
  vectface ListFace = ConvertGAPread_ListFace(dataFAC, n_rows);
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



template<typename T, typename Tgroup>
EquivariantDualDescription<T,Tgroup> Read_BankEntry(std::string const& eFile)
{
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  MyMatrix<T> EXT = POLY_NC_ReadPolytope<T>(dataFile);
  Tgroup GRP = POLY_NC_ReadGroup<Tgroup>(dataFile);
  //
  vectface ListFace = POLY_NC_ReadAllFaces(dataFile);
  return {std::move(EXT), std::move(GRP), std::move(ListFace)};
}

template<typename T, typename Tgroup>
void Write_BankEntry(std::string const& eFile, MyMatrix<T> const& EXT, Tgroup const& GRP, vectface const& ListFace)
{
  if (!FILE_IsFileMakeable(eFile)) {
    std::cerr << "Error in Write_BankEntry: File eFile=" << eFile << " is not makeable\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  POLY_NC_WritePolytope(dataFile, EXT);
  bool orbit_setup = true;
  bool orbit_status = false;
  POLY_NC_WriteGroup(dataFile, GRP, orbit_setup, orbit_status);
  //
  size_t n_orbit = ListFace.size();
  for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++)
    POLY_NC_WriteFace(dataFile, i_orbit, ListFace[i_orbit]);
}





template<typename T, typename Tgroup>
void CheckGroupPolytope(MyMatrix<T> const& EXT, Tgroup const& GRP, std::string const& step)
{
  for (auto & eGen : GRP.GeneratorsOfGroup()) {
    resultFT<T> eRes = FindTransformationGeneral(EXT, EXT, eGen);
    if (!eRes.test) {
      std::cerr << "Error in CheckGroupPolytope at step " << step << "\n";
      throw TerminalException{1};
    }
  }
}



template<typename T, typename Tgroup>
struct TripleCanonic {
  MyMatrix<T> EXT;
  Tgroup GRP;
  std::vector<typename Tgroup::Telt::Tidx> ListIdx;
};



template<typename T, typename Tidx>
std::pair<MyMatrix<T>, std::vector<Tidx>> CanonicalizationPolytopePair(MyMatrix<T> const& EXT, WeightMatrix<T,T> const& WMat)
{
  std::pair<std::vector<Tidx>, std::vector<Tidx>> PairCanonic = GetCanonicalizationVector<T,T,GraphBitset,Tidx>(WMat);
  Tidx n_row=EXT.rows();
  Tidx n_col=EXT.cols();
  MyMatrix<T> EXTcan(n_row, n_col);
  for (int i_row=0; i_row<n_row; i_row++) {
    Tidx j_row=PairCanonic.second[i_row];
    EXTcan.row(i_row) = EXT.row(j_row);
  }
  MyMatrix<T> RowRed = RowReduction(EXTcan);
  MyMatrix<T> EXTret = EXTcan * Inverse(RowRed);
  MyMatrix<T> EXTretB = RemoveFractionMatrix(EXTret);
  //
  return {std::move(EXTretB), std::move(PairCanonic.second)};
}




template<typename T, typename Tgroup>
TripleCanonic<T,Tgroup> CanonicalizationPolytopeTriple(MyMatrix<T> const& EXT, WeightMatrix<T,T> const& WMat)
{
  using Telt=typename Tgroup::Telt;
  using Tidx=typename Telt::Tidx;
  std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> PairCanGrp = GetGroupCanonicalizationVector<T,T,GraphBitset,Tidx>(WMat);
  Tidx n_row=EXT.rows();
  Tidx n_col=EXT.cols();
  MyMatrix<T> EXTcan(n_row, n_col);
  std::vector<Tidx> RevMap(n_row);
  for (Tidx i_row=0; i_row<n_row; i_row++) {
    Tidx j_row=PairCanGrp.first[i_row];
    RevMap[j_row] = i_row;
    EXTcan.row(i_row) = EXT.row(j_row);
  }
  MyMatrix<T> RowRed = RowReduction(EXTcan);
  MyMatrix<T> EXTret = EXTcan * Inverse(RowRed);
  MyMatrix<T> EXTretB = RemoveFractionMatrix(EXTret);
  //
  std::vector<Telt> LGen;
  for (auto & eGen : PairCanGrp.second) {
    std::vector<Tidx> eList(n_row);
    for (Tidx i_row=0; i_row<n_row; i_row++) {
      Tidx i_row2 = PairCanGrp.first[i_row];
      Tidx i_row3 = eGen[i_row2];
      Tidx i_row4 = RevMap[i_row3];
      eList[i_row] = i_row4;
    }
    Telt nGen(eList);
    LGen.push_back(nGen);
  }
  Tgroup GRP(LGen, n_row);
  //
  return {std::move(EXTretB), std::move(GRP), std::move(PairCanGrp.first)};
}


template<typename T>
MyMatrix<T> CanonicalizationPolytope(MyMatrix<T> const& EXT)
{
  WeightMatrix<T, T> WMat=GetWeightMatrix(EXT);
  ReorderingSetWeight(WMat);
  return CanonicalizationPolytopePair<T,int>(EXT, WMat).first;
}



template<typename T, typename Tgroup>
struct DataBank {
private:
  struct PairStore {
    Tgroup GRP;
    vectface ListFace;
  };
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int MinSize;
  UNORD_MAP<MyMatrix<T>, PairStore> ListEnt;
  bool Saving;
  std::string SavingPrefix;
public:
  DataBank(bool const& _Saving, std::string const& _SavingPrefix) : Saving(_Saving), SavingPrefix(_SavingPrefix)
  {
    MinSize = std::numeric_limits<int>::max();
    if (Saving) {
      size_t iOrbit=0;
      while(true) {
        std::string eFileBank = SavingPrefix + "Dual_Desc_" + std::to_string(iOrbit) + ".nc";
        if (!IsExistingFile(eFileBank))
          break;
        EquivariantDualDescription<T,Tgroup> eTriple = Read_BankEntry<T,Tgroup>(eFileBank);
        int e_size = eTriple.EXT.rows();
        std::cerr << "Read iOrbit=" << iOrbit << " FileBank=" << eFileBank << " |EXT|=" << e_size << " |ListFace|=" << eTriple.ListFace.size() << "\n";
        ListEnt[eTriple.EXT] = {eTriple.GRP, eTriple.ListFace};
        MinSize = std::min(MinSize, e_size);
        iOrbit++;
      }
    }
  }
  void InsertEntry(MyMatrix<T> const& EXT, WeightMatrix<T,T> const& WMat, Tgroup const& TheGRPrelevant, bool const& BankSymmCheck, vectface const& ListFace)
  {
    if (!BankSymmCheck) {
      // The computation was already done for the full symmetry group. Only canonic form is needed.
      std::pair<MyMatrix<T>, std::vector<Tidx>> ePair = CanonicalizationPolytopePair<T,Tidx>(EXT, WMat);
      vectface ListFaceO(EXT.rows());
      Telt perm1 = Telt(ePair.second);
      Telt ePerm = ~perm1;
      for (auto & eFace : ListFace) {
        Face eInc = OnFace(eFace, ePerm);
        ListFaceO.push_back(eInc);
      }
      Tgroup GrpConj = TheGRPrelevant.GroupConjugate(ePerm);
      if (Saving) {
        size_t n_orbit = ListEnt.size();
        std::string eFile = SavingPrefix + "DualDesc" + std::to_string(n_orbit) + ".nc";
        std::cerr << "Insert entry to file eFile=" << eFile << "\n";
        Write_BankEntry(eFile, ePair.first, GrpConj, ListFaceO);
      }
      ListEnt[ePair.first] = {GrpConj, ListFaceO};
      int e_size = ePair.first.rows();
      MinSize = std::min(MinSize, e_size);
    } else {
      TripleCanonic<T,Tgroup> eTriple = CanonicalizationPolytopeTriple<T,Tgroup>(EXT, WMat);
      bool NeedRemapOrbit = eTriple.GRP.size() == TheGRPrelevant.size();
      vectface ListFaceO(EXT.rows());
      Telt perm1 = Telt(eTriple.ListIdx);
      Telt ePerm = ~perm1;
      if (!NeedRemapOrbit) {
        // We needed to compute the full group, but it turned out to be the same as the input group.
        for (auto & eFace : ListFace) {
          Face eInc = OnFace(eFace, ePerm);
          ListFaceO.push_back(eInc);
        }
      } else {
        // The full group is bigger than the input group. So we need to reduce.
        UNORD_SET<Face> SetFace;
        for (auto & eFace : ListFace) {
          Face eInc = OnFace(eFace, ePerm);
          Face eIncCan = eTriple.GRP.CanonicalImage(eInc);
          SetFace.insert(eIncCan);
        }
        for (auto & eInc : SetFace) {
          ListFaceO.push_back(eInc);
        }
      }
      if (Saving) {
        size_t n_orbit = ListEnt.size();
        std::string eFile = SavingPrefix + "DualDesc" + std::to_string(n_orbit) + ".nc";
        std::cerr << "Insert entry to file eFile=" << eFile << "\n";
        Write_BankEntry(eFile, eTriple.EXT, eTriple.GRP, ListFaceO);
      }
      ListEnt[eTriple.EXT] = {eTriple.GRP, ListFaceO};
      int e_size = eTriple.EXT.rows();
      MinSize = std::min(MinSize, e_size);
    }
  }
  vectface GetDualDesc(MyMatrix<T> const& EXT, WeightMatrix<T,T> const& WMat, Tgroup const& GRP) const
  {
    std::cerr << "Passing by GetDualDesc |ListEnt|=" << ListEnt.size() << "\n";
    std::pair<MyMatrix<T>, std::vector<Tidx>> ePair = CanonicalizationPolytopePair<T,Tidx>(EXT, WMat);
    auto iter = ListEnt.find(ePair.first);
    if (iter == ListEnt.end())
      return vectface(0); // If returning empty then it means nothing has been found.
    std::cerr << "Finding a matching entry\n";
    vectface ListReprTrans(EXT.rows());
    Telt ePerm = Telt(ePair.second);
    for (auto const& eOrbit : iter->second.ListFace) {
      Face eListJ=OnFace(eOrbit, ePerm);
      ListReprTrans.push_back(eListJ);
    }
    if (GRP.size() == iter->second.GRP.size())
      return ListReprTrans;
    Tgroup GrpConj = iter->second.GRP.GroupConjugate(ePerm);
    return OrbitSplittingListOrbit(GrpConj, GRP, ListReprTrans, std::cerr);
  }
  int get_minsize() const
  {
    return MinSize;
  }
};


size_t get_matching_power(size_t const& val)
{
  size_t pow = 1;
  size_t pos = 0;
  while(true) {
    if (pow >= val) {
      return pos;
    }
    pow *= 2;
    pos++;
  }
}



template<typename Tidx, typename Tint>
std::vector<Tint> GetAllPossibilities(std::map<Tidx,int> const& eMap)
{
  std::vector<Tint> LVal = {1};
  for (auto & kv : eMap) {
    std::vector<Tint> NewVal;
    Tint ePow = 1;
    for (int i=0; i<=kv.second; i++) {
      for (auto & eVal : LVal)
        NewVal.push_back(ePow * eVal);
      ePow *= kv.first;
    }
    LVal = NewVal;
  }
  return LVal;
}





template<typename T, typename Tint, typename Tgroup>
struct DatabaseOrbits {
private:
  using Torbsize=uint16_t;
  using Tidx = typename Tgroup::Telt::Tidx;
  MyMatrix<T> EXT;
  Tgroup GRP;
  Tint groupOrder;
  std::string MainPrefix;
  netCDF::NcFile dataFile;
  FileBool* fb;
  FileFace* ff;
  bool SavingTrigger;
  /* TRICK2: We keep the list of orbit and the map. We could in principle have built the map
     from the start since we know the occurring orders. However, since some orbitsize never occur
     this would have populated it with entries that never occur and so slow it down. */
  UNORD_MAP<Tint,Torbsize> OrbSize_Map;
  std::vector<Tint> ListPossOrbsize;
  struct SingEnt {
    Face face;
    Torbsize idx_orb;
  };
  UNORD_SET<size_t,std::function<size_t(size_t)>,std::function<bool(size_t,size_t)>> DictOrbit;
  std::map<size_t, std::vector<size_t>> CompleteList_SetUndone;
  /* TRICK3: Encoding the pair of face and idx_orb as bits allow us to save memory */
  std::vector<uint8_t> ListOrbit; // This CANNOT be replaced by vectface as we hack our way and so
                                  // making a vectface will not allow
  Tint TotalNumber;
  size_t nbOrbitDone;
  Tint nbUndone;
  size_t nbOrbit;
  size_t n_grpsize;
  /* TRICK3: Knowing the factorization of the order of the group allow us to know exactly
     what are the possible orbitsize occurring and so the number of bits needed to encode them */
  size_t n_bit_orbsize;
  size_t n_act;
  size_t delta;
  size_t n_act_div8;
  std::vector<uint8_t> V_hash;
public:
  // conversion functions that depend only on n_act and n_bit_orbsize.
  SingEnt FaceToSingEnt(Face const& f_in) const
  {
    Face f(n_act);
    for (size_t i=0; i<n_act; i++)
      f[i] = f_in[i];
    size_t i_acc = n_act;
    Torbsize idx_orb=0;
    Torbsize pow = 1;
    for (size_t i=0; i<n_bit_orbsize; i++) {
      idx_orb += Torbsize(f_in[i_acc]) * pow;
      i_acc++;
      pow *= 2;
    }
    return {f,idx_orb};
  }
  Face SingEntToFace(SingEnt const& eEnt) const
  {
    Face f(n_act + n_bit_orbsize);
    for (size_t i=0; i<n_act; i++)
      f[i] = eEnt.face[i];
    size_t work_idx = eEnt.idx_orb;
    size_t i_acc = n_act;
    for (size_t i=0; i<n_bit_orbsize; i++) {
      bool val = work_idx % 2;
      f[i_acc] = val;
      i_acc++;
      work_idx = work_idx / 2;
    }
    return f;
  }
  // Database code that uses ListOrbit;
  SingEnt RetrieveListOrbitEntry(size_t const& i_orb) const
  {
    size_t i_acc = delta * i_orb;
    Face f(n_act);
    for (size_t i=0; i<n_act; i++) {
      f[i] = getbit(ListOrbit, i_acc);
      i_acc++;
    }
    Torbsize idx_orb=0;
    Torbsize pow = 1;
    for (size_t i=0; i<n_bit_orbsize; i++) {
      idx_orb += Torbsize(getbit(ListOrbit, i_acc)) * pow;
      i_acc++;
      pow *= 2;
    }
    return {f,idx_orb};
  }
  void InsertListOrbitEntry(SingEnt const& eEnt)
  {
    size_t curr_len = ListOrbit.size();
    size_t needed_bits = (nbOrbit + 1) * delta;
    size_t needed_len = (needed_bits + 7) / 8;
    for (size_t i=curr_len; i<needed_len; i++)
      ListOrbit.push_back(0);
    // Now setting up the bits,
    size_t i_acc = nbOrbit * delta;
    for (size_t i=0; i<n_act; i++) {
      bool val = eEnt.face[i];
      setbit(ListOrbit, i_acc, val);
      i_acc++;
    }
    size_t work_idx = eEnt.idx_orb;
    for (size_t i=0; i<n_bit_orbsize; i++) {
      bool val = work_idx % 2;
      setbit(ListOrbit, i_acc, val);
      i_acc++;
      work_idx = work_idx / 2;
    }
  }
  void InsertListOrbitFace(Face const& face)
  {
    size_t curr_len = ListOrbit.size();
    size_t needed_bits = (nbOrbit + 1) * delta;
    size_t needed_len = (needed_bits + 7) / 8;
    for (size_t i=curr_len; i<needed_len; i++)
      ListOrbit.push_back(0);
    // Now setting up the bits,
    size_t i_acc = nbOrbit * delta;
    for (size_t i=0; i<n_act; i++) {
      bool val = face[i];
      setbit(ListOrbit, i_acc, val);
      i_acc++;
    }
  }
  void InsertListOrbitIdxOrb(Torbsize const& idx_orb)
  {
    // The position have been
    size_t i_acc = nbOrbit * delta + n_act;
    Torbsize work_idx = idx_orb;
    for (size_t i=0; i<n_bit_orbsize; i++) {
      bool val = work_idx % 2;
      setbit(ListOrbit, i_acc, val);
      i_acc++;
      work_idx = work_idx / 2;
    }
  }
  // Group functionqlities.
  uint16_t GetOrbSizeIndex(Tint const& orbSize)
  {
    /* TRICK4: value 0 is the default constructed one and so using it we can find if the entry is new or not
       in only one call */
    Torbsize & idx = OrbSize_Map[orbSize];
    if (idx == 0) { // A rare case. The linear loop should be totally ok
      auto set=[&]() -> int {
        for (size_t u=0; u<ListPossOrbsize.size(); u++)
          if (ListPossOrbsize[u] == orbSize) {
            return u + 1;
          }
        return 0;
      };
      idx = set();
    }
    return idx - 1;
  }
  void InsertEntryDatabase(Face const& face, bool const& status, size_t const& idx_orb, size_t const& pos)
  {
    if (!status) {
      size_t len = face.count();
      CompleteList_SetUndone[len].push_back(pos);
    }
    Tint orbSize = ListPossOrbsize[idx_orb];
    TotalNumber += orbSize;
    if (status) {
      nbOrbitDone++;
    } else {
      nbUndone += orbSize;
    }
    nbOrbit++;
  }
  DatabaseOrbits(MyMatrix<T> const& _EXT, Tgroup const& _GRP, std::string const& _MainPrefix, bool const& _SavingTrigger) : EXT(_EXT), GRP(_GRP), MainPrefix(_MainPrefix), SavingTrigger(_SavingTrigger)
  {
    TotalNumber = 0;
    nbOrbitDone = 0;
    nbUndone = 0;
    nbOrbit = 0;
    groupOrder = GRP.size();
    std::map<Tidx, int> LFact = GRP.factor_size();
    size_t n_factor = 1;
    for (auto & kv : LFact) {
      n_factor *= (1 + kv.second);
    }
    /* TRICK4: We need to add 1 because of shift by 1 in the OrbSize_Map */
    n_bit_orbsize = get_matching_power(n_factor + 1);
    ListPossOrbsize = GetAllPossibilities<Tidx,Tint>(LFact);

    /* TRICK6: The UNORD_SET only the index and this saves in memory usage. */
    n_act = GRP.n_act();
    delta = n_bit_orbsize + n_act;
    n_act_div8 = (n_act + 7) / 8;
    V_hash = std::vector<uint8_t>(n_act_div8);
    std::function<size_t(size_t)> fctHash=[&](size_t idx) -> size_t {
      for (size_t i=0; i<n_act_div8; i++)
        V_hash[i] = 0;
      size_t pos = delta * idx;
      for (size_t i=0; i<n_act; i++) {
        bool val = getbit(ListOrbit, pos);
        setbit(V_hash, i, val);
        pos++;
      }
      uint32_t seed= 0x1b873560;
      return murmur3_32(V_hash.data(), n_act_div8, seed);
    };
    std::function<bool(size_t,size_t)> fctEqual=[&](size_t idx1, size_t idx2) -> bool {
      size_t pos1 = delta * idx1;
      size_t pos2 = delta * idx2;
      for (size_t i=0; i<n_act; i++) {
        bool val1 = getbit(ListOrbit, pos1);
        bool val2 = getbit(ListOrbit, pos2);
        if (val1 != val2)
          return false;
        pos1++;
        pos2++;
      }
      return true;
    };
    DictOrbit = UNORD_SET<size_t,std::function<size_t(size_t)>,std::function<bool(size_t,size_t)>>({}, fctHash, fctEqual);
    fb = nullptr;
    ff = nullptr;
    if (SavingTrigger) {
      std::string eFileNC = MainPrefix + ".nc";
      std::string eFileFB = MainPrefix + ".fb";
      std::string eFileFF = MainPrefix + ".ff";
      std::cerr << "MainPrefix=" << MainPrefix << "\n";
      size_t n_orbit;
      if (IsExistingFile(eFileNC)) {
        std::cerr << "Opening existing files (NC, FB, FF)\n";
        dataFile.open(eFileNC, netCDF::NcFile::write);
        n_orbit = POLY_NC_ReadNbOrbit(dataFile);
        fb = new FileBool(eFileFB, n_orbit);
        ff = new FileFace(eFileFF, n_act + n_bit_orbsize, n_orbit);
      } else {
        if (!FILE_IsFileMakeable(eFileNC)) {
          std::cerr << "Error in DatabaseOrbits: File eFileNC=" << eFileNC << " is not makeable\n";
          throw TerminalException{1};
        }
        std::cerr << "Creating the files (NC, FB, FF)\n";
        dataFile.open(eFileNC, netCDF::NcFile::replace, netCDF::NcFile::nc4);
        POLY_NC_WritePolytope(dataFile, EXT);
        bool orbit_setup = false;
        bool orbit_status = false;
        POLY_NC_WriteGroup(dataFile, GRP, orbit_setup, orbit_status);
        POLY_NC_SetNbOrbit(dataFile);
        POLY_NC_WriteNbOrbit(dataFile, 0);
        n_orbit = 0;
        fb = new FileBool(eFileFB);
        ff = new FileFace(eFileFF, n_act + n_bit_orbsize);
      }
      netCDF::NcDim fDim = dataFile.getDim("n_grpsize");
      n_grpsize = fDim.getSize();
      //
      for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++) {
        Face f = ff->getface(i_orbit);
        SingEnt eEnt = FaceToSingEnt(f);
        bool status = fb->getbit(i_orbit);
        // The DictOrbit
        InsertListOrbitEntry(eEnt);
        DictOrbit.insert(i_orbit);
        // The other fields
        InsertEntryDatabase(eEnt.face, status, eEnt.idx_orb, i_orbit);
      }
      std::cerr << "Starting with nbOrbitDone=" << nbOrbitDone << " nbUndone=" << nbUndone << " TotalNumber=" << TotalNumber << "\n";
    }
  }
  ~DatabaseOrbits()
  {
    // TRICK5: The destructor does NOT destroy the database! This is because it can be used in another call.
    POLY_NC_WriteNbOrbit(dataFile, nbOrbit);
    if (fb != nullptr)
      delete fb;
    if (ff != nullptr)
      delete ff;
    std::cerr << "Clean closing of the DatabaseOrbits\n";
  }
  void FuncInsert(Face const& face)
  {
    Face face_can = GRP.CanonicalImage(face);
    InsertListOrbitFace(face_can);
    DictOrbit.insert(nbOrbit);
    if (DictOrbit.size() == nbOrbit) // Insertion did not raise the count and so it was already present
      return;
    // Setting up the orbSize
    Tint ordStab = GRP.Stabilizer_OnSets(face).size();
    Tint orbSize = groupOrder / ordStab;
    Torbsize idx_orb = GetOrbSizeIndex(orbSize);
    InsertListOrbitIdxOrb(idx_orb);
    InsertEntryDatabase(face_can, false, idx_orb, nbOrbit);
    //
    if (SavingTrigger) {
      Face f = SingEntToFace({face_can, idx_orb});
      fb->setbit(nbOrbit - 1, false);
      ff->setface(nbOrbit - 1, f);
    }
  }
  void FuncPutOrbitAsDone(size_t const& iOrb)
  {
    const SingEnt & eEnt = RetrieveListOrbitEntry(iOrb);
    if (SavingTrigger) {
      fb->setbit(iOrb, true);
    }
    size_t len = eEnt.face.count();
    /* TRICK1: We copy the last element in first position to erase it and then pop_back the vector */
    std::vector<size_t> & V = CompleteList_SetUndone[len];
    V[0] = V[V.size()-1];
    V.pop_back();
    nbUndone -= ListPossOrbsize[eEnt.idx_orb];
    nbOrbitDone++;
  }
  Face ComputeIntersectionUndone() const
  {
    size_t n_row = EXT.rows();
    Face eSetReturn(n_row);
    for (size_t i_row=0; i_row<n_row; i_row++)
      eSetReturn[i_row] = 1;
    for (auto & eEnt : CompleteList_SetUndone) {
      for (auto & pos : eEnt.second) {
        const Face & eFace = RetrieveListOrbitEntry(pos).face;
        eSetReturn &= OrbitIntersection(GRP, eFace);
        if (eSetReturn.count() == 0)
          return eSetReturn;
      }
    }
    return eSetReturn;
  }
  vectface FuncListOrbitIncidence()
  {
    if (SavingTrigger) {
      std::string eFileNC = MainPrefix + ".nc";
      dataFile.close();
      RemoveFile(eFileNC);
      //
      std::string eFileFB = MainPrefix + ".fb";
      delete fb;
      fb = nullptr;
      RemoveFile(eFileFB);
      //
      std::string eFileFF = MainPrefix + ".ff";
      delete ff;
      ff = nullptr;
      RemoveFile(eFileFF);
    }
    DictOrbit.clear();
    CompleteList_SetUndone.clear();
    vectface retListOrbit(n_act);
    for (size_t i_orbit=0; i_orbit<nbOrbit; i_orbit++) {
      Face f = RetrieveListOrbitEntry(i_orbit).face;
      retListOrbit.push_back(f);
    }
    return retListOrbit;
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
      size_t len = eEnt.second.size();
      if (len > 0) {
        /* TRICK1: Take the first element in the vector. This first element will remain
           in place but the vector will be extended. */
        size_t pos = eEnt.second[0];
        Face f = RetrieveListOrbitEntry(pos).face;
        return {pos, f};
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
vectface DUALDESC_AdjacencyDecomposition(
         DataBank<T,Tgroup> & TheBank,
	 MyMatrix<T> const& EXT,
	 Tgroup const& GRP,
	 PolyHeuristicSerial<typename Tgroup::Tint> const& AllArr,
	 std::string const& ePrefix,
	 int const& TheLevel)
{
  if (ExitEvent) {
    std::cerr << "Terminating the program by Ctrl-C\n";
    throw TerminalException{1};
  }
  using Tint=typename Tgroup::Tint;
  int nbRow=EXT.rows();
  int eRank=EXT.cols();
  WeightMatrix<T, T> WMat;
  bool HaveWMat=false;
  auto ComputeWMat=[&]() -> void {
    if (HaveWMat)
      return;
    WMat=GetWeightMatrix(EXT);
    ReorderingSetWeight(WMat);
  };
  //
  // Checking if the entry is present in the map.
  //
  if (nbRow >= TheBank.get_minsize()) {
    ComputeWMat();
    vectface ListFace = TheBank.GetDualDesc(EXT, WMat, GRP);
    if (ListFace.size() > 0)
      return ListFace;
  }
  //
  // Now computing the groups
  //
  Tgroup TheGRPrelevant;
  std::map<std::string, Tint> TheMap;
  int nbVert=EXT.rows();
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
  vectface ListOrbitFaces(nbRow);
  std::string ansSymm;
  // 3 scenario
  // --- 1 : We have the full symmetry group and the computation was done with respect to it.
  // --- 2 : We have computed for a subgroup which actually is the full group.
  // --- 3 : We have computed for a subgroup which actually is a strict subgroup.
  bool BankSymmCheck;
  if (ansSplit != "split") {
    TheGRPrelevant = GRP;
    std::string ansProg=HeuristicEvaluation(TheMap, AllArr.DualDescriptionProgram);
    ListOrbitFaces = DirectFacetOrbitComputation(EXT, GRP, ansProg, std::cerr);
    BankSymmCheck = true;
    ansSymm = "no";
  } else {
    ansSymm = HeuristicEvaluation(TheMap, AllArr.AdditionalSymmetry);
    std::cerr << "ansSymm=" << ansSymm << "\n";
    if (ansSymm == "yes") {
      ComputeWMat();
      TheGRPrelevant = GetStabilizerWeightMatrix<T,T,Tgroup>(WMat);
      BankSymmCheck = false;
    } else {
      TheGRPrelevant = GRP;
      BankSymmCheck = true;
    }
    Tint GroupSizeComp = TheGRPrelevant.size();
    std::cerr << "RESPAWN a new ADM computation |GRP|=" << GroupSizeComp << " TheDim=" << (eRank-1) << " |EXT|=" << nbRow << "\n";
    TheMap["groupsizerelevant"] = GroupSizeComp;
    std::string MainPrefix = ePrefix + "Database_" + std::to_string(TheLevel) + "_" + std::to_string(nbVert) + "_" + std::to_string(eRank);
    std::cerr << "MainPrefix=" << MainPrefix << "\n";
    bool SavingTrigger=AllArr.Saving;
    DatabaseOrbits<T,Tint,Tgroup> RPL(EXT, TheGRPrelevant, MainPrefix, SavingTrigger);
    int NewLevel = TheLevel + 1;
    if (RPL.FuncNumberOrbit() == 0) {
      std::string ansSamp=HeuristicEvaluation(TheMap, AllArr.InitialFacetSet);
      vectface ListFace=DirectComputationInitialFacetSet(EXT, ansSamp);
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
      FlippingFramework<T> FF(EXT, eInc);
      Tgroup TheStab=TheGRPrelevant.Stabilizer_OnSets(eInc);
      Tint OrbSize=GroupSizeComp / TheStab.size();
      Tgroup GRPred=ReducedGroupAction(TheStab, eInc);
      std::cerr << "Considering orbit " << SelectedOrbit << " |EXT|=" << eInc.size() << " |inc|=" << eInc.count() << " Level=" << TheLevel << " |stab|=" << GRPred.size() << " dim=" << TheDim << "\n";
      std::string NewPrefix = ePrefix + "ADM" + std::to_string(SelectedOrbit) + "_";
      std::cerr << "NewPrefix=" << NewPrefix << "\n";
      vectface TheOutput=DUALDESC_AdjacencyDecomposition(TheBank, FF.Get_EXT_face(), GRPred, AllArr, NewPrefix, NewLevel);
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
    TheBank.InsertEntry(EXT, WMat, TheGRPrelevant, BankSymmCheck, ListOrbitFaces);
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
  // Setting up the Control C event.
  ExitEvent = false;
  signal(SIGINT, signal_callback_handler);
  //
  using Tint=typename Tgroup::Tint;
  SingleBlock BlockBANK=eFull.ListBlock.at("BANK");
  bool BANK_IsSaving=BlockBANK.ListBoolValues.at("Saving");
  std::string BANK_Prefix=BlockBANK.ListStringValues.at("Prefix");
  DataBank<T,Tgroup> TheBank(BANK_IsSaving, BANK_Prefix);
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
  MyMatrix<T> EXTred=ColumnReduction(EXT);
  vectface TheOutput=DUALDESC_AdjacencyDecomposition(TheBank, EXTred, GRP, AllArr, DD_Prefix, TheLevel);
  std::cerr << "|TheOutput|=" << TheOutput.size() << "\n";
  //
  std::ofstream OUTfs(OUTfile);
  VectVectInt_Magma_Print(OUTfs, TheOutput);
  //
}



#endif

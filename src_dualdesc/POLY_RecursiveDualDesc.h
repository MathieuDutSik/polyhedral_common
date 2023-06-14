// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_RECURSIVEDUALDESC_H_
#define SRC_DUALDESC_POLY_RECURSIVEDUALDESC_H_

// clang-format off
#include "GRP_DoubleCoset.h"
#include "MAT_MatrixInt.h"
#include "Namelist.h"
#include "POLY_DirectDualDesc.h"
#include "POLY_Heuristics.h"
#include "POLY_Kskeletton.h"
#include "POLY_SamplingFacet.h"
#include "Temp_PolytopeEquiStab.h"
#include "Timings.h"

#include "POLY_GAP.h"
#include "Balinski_basic.h"
#include "Databank.h"
#include "MatrixGroupBasic.h"
#include "basic_datafile.h"
#include <limits>
#include <set>
#include <map>
#include <signal.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

// #define MURMUR_HASH
// #define ROBIN_HOOD_HASH
#define SUBSET_HASH

// #define UNORDERED_MAP
#define TSL_SPARSE_MAP
// #define TSL_ROBIN_MAP
// #define TSL_HOPSCOTCH_MAP

#ifdef UNORDERED_MAP
#define UNORD_MAP std::unordered_map
#define UNORD_SET std::unordered_set
#endif

#ifdef TSL_SPARSE_MAP
#include "sparse_map.h"
#include "sparse_set.h"
#define UNORD_MAP tsl::sparse_map
#define UNORD_SET tsl::sparse_set
#endif

#ifdef TSL_ROBIN_MAP
#include "robin_map.h"
#include "robin_set.h"
#define UNORD_MAP tsl::robin_map
#define UNORD_SET tsl::robin_set
#endif

#ifdef TSL_HOPSCOTCH_MAP
#include "hopscotch_map.h"
#include "hopscotch_set.h"
#define UNORD_MAP tsl::hopscotch_map
#define UNORD_SET tsl::hopscotch_set
#endif

static const int CANONIC_STRATEGY__CANONICAL_IMAGE = 0;
static const int CANONIC_STRATEGY__STORE = 1;
static const int CANONIC_STRATEGY__INITIAL_TRIV = 2;

static const int CANONIC_STRATEGY__DEFAULT = CANONIC_STRATEGY__CANONICAL_IMAGE;
static const int REPR_STRATEGY__DEFAULT = 0; // More or less irrelevant here


std::atomic<bool> ExitEvent;

void signal_callback_handler(int signum) {
  std::cout << "Caught signal " << signum << "\n";
  std::cout << "We are going to exit hopefully\n";
  ExitEvent = true;
}

template <typename T, typename Tgroup> struct EquivariantDualDescription {
  MyMatrix<T> EXT;
  Tgroup GRP;
  vectface ListFace;
};

template <typename T, typename Tgroup>
EquivariantDualDescription<T, Tgroup> ConvertGAPread_EquivDualDesc(
    datagap::DataGAP<T, typename Tgroup::Telt> const &dataEXT,
    datagap::DataGAP<T, typename Tgroup::Telt> const &dataFAC) {
  if (dataEXT.Nature != datagap::int_record) {
    std::cerr << "For EquivDualDesc, we need to have a record as entry\n";
    throw TerminalException{1};
  }
  int pos_EXT = -1;
  int pos_GRP = -1;
  int n_pos = dataEXT.ListRec.size();
  for (int pos = 0; pos < n_pos; pos++) {
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
  MyMatrix<T> EXT =
      datagap::ConvertGAPread_MyMatrixT(dataEXT.ListRec[pos_EXT].second);
  int n_rows = EXT.rows();
  Tgroup GRP = datagap::ConvertGAPread_PermutationGroup<T, Tgroup>(
      dataEXT.ListRec[pos_GRP].second, n_rows);
  //
  vectface ListFace = ConvertGAPread_ListFace(dataFAC, n_rows);
  //
  return {std::move(EXT), std::move(GRP), std::move(ListFace)};
}

template <typename T, typename Tgroup>
void CheckGroupPolytope(MyMatrix<T> const &EXT, Tgroup const &GRP,
                        std::string const &step) {
  for (auto &eGen : GRP.GeneratorsOfGroup()) {
    std::optional<MyMatrix<T>> opt = FindTransformationGeneral(EXT, EXT, eGen);
    if (!opt) {
      std::cerr << "Error in CheckGroupPolytope at step " << step << "\n";
      throw TerminalException{1};
    }
  }
}

template <typename T, typename Tgroup> struct TripleCanonic {
  MyMatrix<T> EXT;
  Tgroup GRP;
  std::vector<typename Tgroup::Telt::Tidx> ListIdx;
};

template <typename T, typename Tidx, typename Tidx_value>
std::pair<MyMatrix<T>, std::vector<Tidx>>
CanonicalizationPolytopePair(MyMatrix<T> const &EXT,
                             WeightMatrix<true, T, Tidx_value> const &WMat) {
  std::vector<Tidx> CanonicOrd =
      GetCanonicalizationVector_Kernel<T, GraphBitset, Tidx>(WMat);
  Tidx n_row = EXT.rows();
  Tidx n_col = EXT.cols();
  MyMatrix<T> EXTcan(n_row, n_col);
  for (Tidx i_row = 0; i_row < n_row; i_row++) {
    Tidx j_row = CanonicOrd[i_row];
    EXTcan.row(i_row) = EXT.row(j_row);
  }
  MyMatrix<T> EXTret = CanonicalizeOrderedMatrix(EXTcan);
  return {std::move(EXTret), std::move(CanonicOrd)};
}

template <typename T, typename Tgroup, typename Tidx_value>
TripleCanonic<T, Tgroup>
CanonicalizationPolytopeTriple(MyMatrix<T> const &EXT,
                               WeightMatrix<true, T, Tidx_value> const &WMat) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> PairCanGrp =
      GetGroupCanonicalizationVector_Kernel<T, GraphBitset, Tidx>(WMat);
  Tidx n_row = EXT.rows();
  Tidx n_col = EXT.cols();
  MyMatrix<T> EXTcan(n_row, n_col);
  std::vector<Tidx> RevMap(n_row);
  for (Tidx i_row = 0; i_row < n_row; i_row++) {
    Tidx j_row = PairCanGrp.first[i_row];
    RevMap[j_row] = i_row;
    EXTcan.row(i_row) = EXT.row(j_row);
  }
  MyMatrix<T> RowRed = RowReduction(EXTcan);
  MyMatrix<T> EXTret = EXTcan * Inverse(RowRed);
  MyMatrix<T> EXTretB = ScalarCanonicalizationMatrix(EXTret);
  //
  std::vector<Telt> LGen;
  for (auto &eGen : PairCanGrp.second) {
    std::vector<Tidx> eList(n_row);
    for (Tidx i_row = 0; i_row < n_row; i_row++) {
      Tidx i_row2 = PairCanGrp.first[i_row];
      Tidx i_row3 = eGen[i_row2];
      Tidx i_row4 = RevMap[i_row3];
      eList[i_row] = i_row4;
    }
    Telt nGen(eList);
    LGen.emplace_back(std::move(nGen));
  }
  Tgroup GRP(LGen, n_row);
  //
  return {std::move(EXTretB), std::move(GRP), std::move(PairCanGrp.first)};
}

template <typename T>
MyMatrix<T> CanonicalizationPolytope(MyMatrix<T> const &EXT) {
  using Tidx_value = uint16_t;
  WeightMatrix<true, T, Tidx_value> WMat = GetWeightMatrix<T, Tidx_value>(EXT);
  WMat.ReorderingSetWeight();
  return CanonicalizationPolytopePair<T, int, Tidx_value>(EXT, WMat).first;
}

template <typename Tidx>
std::pair<size_t, size_t> get_delta(const std::map<Tidx, int> &LFact,
                                    const size_t &n_act) {
  size_t n_factor = 1;
  for (auto &kv : LFact) {
    n_factor *= (1 + kv.second);
  }
  /* TRICK 4: We need to add 1 because of shift by 1 in the OrbSize_Map */
  size_t n_bit_orbsize = get_matching_power(n_factor + 1);
  size_t delta = n_bit_orbsize + n_act;
  return {n_bit_orbsize, delta};
}

template <typename Tidx, typename Tint>
std::vector<Tint> GetAllPossibilities(std::map<Tidx, int> const &eMap) {
  std::vector<Tint> LVal = {1};
  for (auto &kv : eMap) {
    std::vector<Tint> NewVal;
    Tint ePow = 1;
    // Needed conversion below because of Mac typing issues.
    Tint kv_first_tint = size_t(kv.first);
    for (int i = 0; i <= kv.second; i++) {
      for (auto &eVal : LVal)
        NewVal.push_back(ePow * eVal);
      ePow *= kv_first_tint;
    }
    LVal = NewVal;
  }
  return LVal;
}

/*
  This is a simple canonicalization function that does not return the Orbitsize
 */
template <typename Tgroup>
Face CanonicalImageDualDesc(int const& method_choice, Tgroup const& GRP, Face const& f) {
  if (method_choice == CANONIC_STRATEGY__CANONICAL_IMAGE)
    return GRP.CanonicalImage(f);
  if (method_choice == CANONIC_STRATEGY__STORE)
    return GRP.StoreCanonicalImage(f);
  if (method_choice == CANONIC_STRATEGY__INITIAL_TRIV)
    return GRP.CanonicalImageInitialTriv(f);
  std::cerr << "Error in CanonicalImageDualDesc, no method found\n";
  std::cerr << "method_choice=" << method_choice << "\n";
  throw TerminalException{1};
}

template<typename Torbsize, typename Tgroup>
struct DataFaceOrbitSize {
  using Tint = typename Tgroup::Tint;
  /* TRICK 3: Knowing the factorization of the order of the group allow us to
     know exactly
     what are the possible orbitsize occurring and so the number of bits needed
     to encode them */
  std::vector<Tint> ListPossOrbsize;
  /* TRICK 2: We keep the list of orbit and the map. We could in principle have
     built the map from the start since we know the occurring orders. However,
     since some orbitsize never occur
     this would have populated it with entries that never occur and so slow it
     down. */
  UNORD_MAP<Tint, Torbsize> OrbSize_Map;
  size_t n;
  size_t n_bit_orbsize;
  size_t delta;
  DataFaceOrbitSize(Tgroup const& GRP) {
    using Tidx = typename Tgroup::Telt::Tidx;
    std::map<Tidx, int> LFact = GRP.factor_size();
    n = GRP.n_act();
    std::pair<size_t, size_t> ep = get_delta(LFact, n);
    ListPossOrbsize = GetAllPossibilities<Tidx, Tint>(LFact);
    n_bit_orbsize = ep.first;
    delta = ep.second;
  }
  Torbsize GetOrbSizeIndex(Tint const &orbSize) {
    /* TRICK 4: value 0 is the default constructed one and so using it we can
       find if the entry is new or not in only one call */
    Torbsize &idx = OrbSize_Map[orbSize];
    if (idx == 0) {
      // A rare case. The linear loop should be totally ok
      auto set = [&]() -> int {
        for (size_t u = 0; u < ListPossOrbsize.size(); u++)
          if (ListPossOrbsize[u] == orbSize) {
            return u + 1;
          }
        return 0;
      };
      idx = set();
    }
    return idx - 1;
  }
  Face ConvertFaceOrbitSize(std::pair<Face,Tint> const& pair) {
    Face const& f = pair.first;
    Tint const& orbitSize = pair.second;
    Torbsize idx_orb = GetOrbSizeIndex(orbitSize);
    //
    Face f_ret(delta);
    for (size_t i=0; i<n; i++)
      f_ret[i] = f[i];
    size_t work_idx = idx_orb;
    size_t i_acc = n;
    for (size_t i = 0; i < n_bit_orbsize; i++) {
      bool val = work_idx % 2;
      f_ret[i_acc] = val;
      i_acc++;
      work_idx = work_idx / 2;
    }
    return f_ret;
  }
  std::pair<Face,Tint> ConvertFace(Face const& f) const {
    Face f_ret(n);
    for (size_t i=0; i<n; i++)
      f_ret[i] = f[i];
    size_t idx_orb = 0;
    size_t pow = 1;
    size_t i_acc = n;
    for (size_t i = 0; i < n_bit_orbsize; i++) {
      if (f[i_acc] == 1)
        idx_orb += pow;
      i_acc++;
      pow *= 2;
    }
    return {f_ret, ListPossOrbsize[idx_orb]};
  }
};

/*
  Return the canonical form and the orbit stabilizer if available and encoded as
  a face (which creates a lot of possibility of errors)
 */
template <typename Torbsize, typename Tgroup>
Face CanonicalImageGeneralDualDesc(int const& method_choice, Tgroup const& GRP, DataFaceOrbitSize<Torbsize, Tgroup> & recConvert, Face const& f) {
  using Tint = typename Tgroup::Tint;
  if (method_choice == CANONIC_STRATEGY__CANONICAL_IMAGE) {
    std::pair<Face,Tint> pair = GRP.CanonicalImageOrbitSize(f);
    return recConvert.ConvertFaceOrbitSize(pair);
  }
  if (method_choice == CANONIC_STRATEGY__STORE) {
    std::pair<Face,Tint> pair = GRP.StoreCanonicalImageOrbitSize(f);
    return recConvert.ConvertFaceOrbitSize(pair);
  }
  if (method_choice == CANONIC_STRATEGY__INITIAL_TRIV)
    return GRP.CanonicalImageInitialTriv(f);
  std::cerr << "Error in CanonicalImageOrbitSizeDualDesc, no method found\n";
  std::cerr << "method_choice=" << method_choice << "\n";
  throw TerminalException{1};
}

vectface vectface_reduction(vectface const& vf, size_t n_red) {
  vectface vf_red(n_red);
  Face f_red(n_red);
  for (auto & f : vf) {
    for (size_t i=0; i<n_red; i++)
      f_red[i] = f[i];
    vf_red.push_back(f_red);
  }
  return vf_red;
}

template<typename Tint>
struct FaceOrbitsizeTableContainer {
public:
  std::vector<Tint> ListPossOrbsize;
  size_t n;
  vectface vfo;
  FaceOrbitsizeTableContainer(std::vector<Tint> const& _ListPossOrbsize, size_t _n, vectface && _vfo) : ListPossOrbsize(std::move(_ListPossOrbsize)), n(_n), vfo(std::move(_vfo)) {
  }
  template<typename Tgroup>
  FaceOrbitsizeTableContainer(vectface const& vf, Tgroup const& GRP) {
    n = vf.get_n();
    using Tidx=typename Tgroup::Telt::Tidx;
    std::map<Tidx, int> LFact = GRP.factor_size();
    std::pair<size_t, size_t> ep = get_delta(LFact, n);
    size_t n_bit_orbsize = ep.first;
    size_t delta = ep.second;
    ListPossOrbsize = GetAllPossibilities<Tidx, Tint>(LFact);
    UNORD_MAP<Tint, size_t> OrbSize_Map;
    for (size_t i=0; i<ListPossOrbsize.size(); i++) {
      OrbSize_Map[ListPossOrbsize[i]] = i;
    }
    vectface vf_ins(delta);
    vfo = std::move(vf_ins);
    for (auto & eFace : vf) {
      Tint orbitSize = GRP.OrbitSize_OnSets(eFace);
      size_t idx_orb = OrbSize_Map[orbitSize];
      Face f(delta);
      for (size_t i=0; i<n; i++)
        f[i] = eFace[i];
      size_t work_idx = idx_orb;
      size_t i_acc = n;
      for (size_t i = 0; i < n_bit_orbsize; i++) {
        bool val = work_idx % 2;
        f[i_acc] = val;
        i_acc++;
        work_idx = work_idx / 2;
      }
      vfo.push_back(f);
    }
  }
  std::pair<Face,Tint> GetPair(size_t const& idx) const {
    Face f = vfo[idx];
    Face f_red(n);
    for (size_t i=0; i<n; i++)
      f_red[i] = f[i];
    size_t pow = 1;
    size_t idx_orb = 0;
    for (size_t i = n; i < vfo.get_n(); i++) {
      if (f[i] == 1)
        idx_orb += pow;
      pow *= 2;
    }
    Tint orbSize = ListPossOrbsize[idx_orb];
    return {f_red, orbSize};
  }
  size_t size() const {
    return vfo.size();
  }
  vectface GetListFaces() const {
    return vectface_reduction(vfo, n);
  }
};

template<typename Telt>
Telt trivial_extension(Telt const& ePerm, size_t const& delta) {
  using Tidx = typename Telt::Tidx;
  Tidx delta_tidx = delta;
  std::vector<Tidx> V(delta);
  Tidx siz = ePerm.size();
  for (Tidx i=0; i<siz; i++) {
    Tidx ePt = ePerm.at(i);
    V[i] = ePt;
  }
  for (Tidx i=siz; i<delta_tidx; i++) {
    V[i] = i;
  }
  return Telt(V);
}

template<typename Tgroup>
Tgroup trivial_extension_group(Tgroup const& eGroup, size_t const& delta) {
  using Telt=typename Tgroup::Telt;
  using Tidx=typename Telt::Tidx;
  std::vector<Telt> LGen = eGroup.GeneratorsOfGroup();
  std::vector<Telt> LGenExt;
  for (auto & eGen : LGen) {
    Telt eGenExt = trivial_extension(eGen, delta);
    LGenExt.push_back(eGenExt);
  }
  std::vector<Tidx> V(delta);
  for (size_t i=0; i<delta; i++) {
    V[i] = i;
  }
  Telt id(V);
  return Tgroup(LGenExt, id);
}

template <typename Tgroup>
std::vector<int> GetPossibleCanonicalizationMethod(Tgroup const &GRP) {
  // We put first the CANONIC_STRATEGY__CANONICAL_IMAGE as it is an all around reasonable method
  // on which other methods have to compete with.
  std::vector<int> list_considered = {CANONIC_STRATEGY__CANONICAL_IMAGE, CANONIC_STRATEGY__INITIAL_TRIV};
  if (GRP.size() < 20000) { // We need to exclude that strategy if too large as that strategy has no chance.
    list_considered.push_back(CANONIC_STRATEGY__STORE);
  }
  return list_considered;
}

template <typename Tgroup>
int64_t time_evaluation_can_method(int const& method, vectface const& vf, Tgroup const &GRP, int64_t upper_limit) {
  NanosecondTime time;
  int64_t duration = 0;
  int64_t miss_val = std::numeric_limits<int64_t>::max();
  for (auto & f : vf) {
    (void)CanonicalImageDualDesc(method, GRP, f);
    duration = time.const_eval_int64();
    if (duration > upper_limit)
      return miss_val;
  }
  return duration;
}

template <typename Tgroup>
int GetCanonicalizationMethod_Serial(vectface const& vf, Tgroup const &GRP) {
  std::vector<int> list_considered = GetPossibleCanonicalizationMethod(GRP);
  int64_t upper_limit = std::numeric_limits<int64_t>::max();
  int chosen_method = -1;
  for (auto& method : list_considered) {
    int64_t runtime = time_evaluation_can_method(method, vf, GRP, upper_limit);
    if (runtime < upper_limit) {
      chosen_method = method;
      upper_limit = runtime;
    }
  }
  return chosen_method;
}

template <typename T, typename Tgroup>
int GetCanonicalizationMethodRandom(MyMatrix<T> const &EXT, Tgroup const &GRP, size_t size) {
  if (size < 10000) {
    if (GRP.size() < 200)
      return CANONIC_STRATEGY__STORE;
    return CANONIC_STRATEGY__INITIAL_TRIV;
  }
  int n = EXT.rows();
  vectface vf(n);
  for (int iter=0; iter<100; iter++) {
    Face f = RandomFace(n);
    vf.push_back(f);
  }
  return GetCanonicalizationMethod_Serial(vf, GRP);
}


template <typename T, typename Tgroup, typename Tidx_value>
std::pair<MyMatrix<T>, TripleStore<Tgroup>> GetCanonicalInformation(
    MyMatrix<T> const &EXT, WeightMatrix<true, T, Tidx_value> const &WMat,
    Tgroup const &TheGRPrelevant,
    FaceOrbitsizeTableContainer<typename Tgroup::Tint> const& ListOrbitFaceOrbitsize) {
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
  std::vector<Tint> ListPossOrbSize = ListOrbitFaceOrbitsize.ListPossOrbsize;
  TripleCanonic<T, Tgroup> eTriple =
      CanonicalizationPolytopeTriple<T, Tgroup>(EXT, WMat);
  bool NeedRemapOrbit = eTriple.GRP.size() == TheGRPrelevant.size();
  size_t delta = ListOrbitFaceOrbitsize.vfo.size();
  Telt perm1 = Telt(eTriple.ListIdx);
  Telt ePerm = ~perm1;
  Telt ePermExt = trivial_extension(ePerm, delta);
  vectface ListFaceO(delta);
  Face eFaceImg(delta);
  if (!NeedRemapOrbit) {
    // We needed to compute the full group, but it turned out to be the same
    // as the input group.
    for (auto &eFace : ListOrbitFaceOrbitsize.vfo) {
      OnFace_inplace(eFaceImg, eFace, ePermExt);
      ListFaceO.push_back(eFaceImg);
    }
  } else {
    // The full group is bigger than the input group. So we need to reduce.
    // The used method for canonicalization does not matter, so everything
    // is correct.
    size_t size = ListOrbitFaceOrbitsize.size();
    int can_method = GetCanonicalizationMethodRandom(EXT, TheGRPrelevant, size);
    UNORD_SET<Face> SetFace;
    Tgroup GRPext = trivial_extension_group(eTriple.GRP, delta);
    for (auto &eFace : ListOrbitFaceOrbitsize.vfo) {
      OnFace_inplace(eFaceImg, eFace, ePermExt);
      Face eIncCan = CanonicalImageDualDesc(can_method, GRPext, eFaceImg);
      SetFace.insert(eIncCan);
    }
    for (auto &eInc : SetFace) {
      ListFaceO.push_back(eInc);
    }
  }
  TripleStore<Tgroup> ePair{eTriple.GRP, std::move(ListPossOrbSize), std::move(ListFaceO)};
  return {std::move(eTriple.EXT), std::move(ePair)};
}

template <typename Tbank, typename T, typename Tgroup, typename Tidx_value>
void insert_entry_in_bank(Tbank &bank, MyMatrix<T> const &EXT,
                          WeightMatrix<true, T, Tidx_value> const &WMat,
                          Tgroup const &TheGRPrelevant,
                          bool const &BankSymmCheck,
                          FaceOrbitsizeTableContainer<typename Tgroup::Tint> const& ListOrbitFaceOrbitsize) {
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
  using Tidx = typename Telt::Tidx;
  size_t delta = ListOrbitFaceOrbitsize.vfo.get_n();
  if (!BankSymmCheck) {
    // The computation was already done for the full symmetry group. Only
    // canonic form is needed.
    std::pair<MyMatrix<T>, std::vector<Tidx>> ePair =
        CanonicalizationPolytopePair<T, Tidx, Tidx_value>(EXT, WMat);
    vectface ListFaceO(delta);
    Telt perm1 = Telt(ePair.second);
    Telt ePerm = ~perm1;
    Telt ePermExt = trivial_extension(ePerm, delta);
    Face eFaceImg(delta);
    for (auto &eFace : ListOrbitFaceOrbitsize.vfo) {
      OnFace_inplace(eFaceImg, eFace, ePermExt);
      ListFaceO.push_back(eFaceImg);
    }
    Tgroup GrpConj = TheGRPrelevant.GroupConjugate(ePerm);
    std::vector<Tint> ListPossOrbSize = ListOrbitFaceOrbitsize.ListPossOrbsize;
    bank.InsertEntry(std::move(ePair.first),
                     {std::move(GrpConj), std::move(ListPossOrbSize), std::move(ListFaceO)});
  } else {
    std::pair<MyMatrix<T>, TripleStore<Tgroup>> eP =
        GetCanonicalInformation(EXT, WMat, TheGRPrelevant, ListOrbitFaceOrbitsize);
    bank.InsertEntry(std::move(eP.first), std::move(eP.second));
  }
}

template<typename Tgroup, typename Tface_orbitsize>
vectface OrbitSplittingListOrbitGen(const Tgroup& GRPbig, const Tgroup& GRPsma, Tface_orbitsize const& ListFaceOrbitsize,
                                    PolyHeuristicSerial<typename Tgroup::Tint> &AllArr, std::ostream & os) {
  using Tint = typename Tgroup::Tint;
  std::map<std::string, Tint> TheMap;
  Tint ordGRPbig = GRPbig.size();
  Tint ordGRPsma = GRPsma.size();
  if (ordGRPbig == ordGRPsma) {
    return ListFaceOrbitsize.GetListFaces();
  }
  Tint index = ordGRPbig / ordGRPsma;
  TheMap["groupsize_big"] = ordGRPbig;
  TheMap["groupsize_sma"] = ordGRPsma;
  TheMap["index"] = index;
  TheMap["n_orbit"] = ListFaceOrbitsize.size();
  std::string method_split = HeuristicEvaluation(TheMap, AllArr.OrbitSplitTechnique);
  return OrbitSplittingListOrbit_spec(GRPbig, GRPsma, ListFaceOrbitsize, method_split, os);
}

template <typename Tbank, typename T, typename Tgroup, typename Tidx_value>
vectface getdualdesc_in_bank(Tbank &bank, MyMatrix<T> const &EXT,
                             WeightMatrix<true, T, Tidx_value> const &WMat,
                             Tgroup const &GRP,
                             PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
                             std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
  using Tidx = typename Telt::Tidx;
  std::pair<MyMatrix<T>, std::vector<Tidx>> ePair =
      CanonicalizationPolytopePair<T, Tidx, Tidx_value>(EXT, WMat);
  const TripleStore<Tgroup> &RecAns = bank.GetDualDesc(ePair.first);
  if (RecAns.ListFace.size() == 0) {
    return vectface(0);
  }
  os << "Finding a matching entry\n";
  Telt ePerm = Telt(ePair.second);
  size_t n = EXT.rows();
  if (GRP.size() == RecAns.GRP.size()) {
    vectface ListReprTrans(n);
    Face eFaceImg(n);
    vectface ListFace = vectface_reduction(RecAns.ListFace, n);
    for (auto const &eFace : ListFace) {
      OnFace_inplace(eFaceImg, eFace, ePerm);
      ListReprTrans.push_back(eFaceImg);
    }
    return ListReprTrans;
  }
  Tgroup GrpConj = RecAns.GRP.GroupConjugate(ePerm);
  size_t delta = RecAns.ListFace.get_n();
  Telt ePermExt = trivial_extension(ePerm, delta);
  vectface ListReprTrans(delta);
  Face eFaceImg(delta);
  for (auto const &eFace : RecAns.ListFace) {
    OnFace_inplace(eFaceImg, eFace, ePermExt);
    ListReprTrans.push_back(eFaceImg);
  }
  FaceOrbitsizeTableContainer<Tint> fotc(RecAns.ListPossOrbsize, n, std::move(ListReprTrans));
  return OrbitSplittingListOrbitGen(GrpConj, GRP, fotc, AllArr, os);
}

template <typename T, typename Torbsize, typename Tgroup> struct DataFacetCan {
  using Tint = typename Tgroup::Tint;
  size_t SelectedOrbit;
  Face eInc;
  FlippingFramework<T> FF;
  const Tgroup &GRP;
  DataFaceOrbitSize<Torbsize, Tgroup> & recConvert;
  Tgroup Stab;
  int can_method;
  Face FlipFace(const Face &f, [[maybe_unused]] std::ostream & os) const {
#ifdef TIMINGS
    MicrosecondTime time;
#endif
    Face eFlip = FF.FlipFace(f);
#ifdef TIMINGS
    os << "|FlipFace|=" << time << "\n";
#endif
    Face result = CanonicalImageGeneralDualDesc(can_method, GRP, recConvert, eFlip);
#ifdef TIMINGS
    os << "|canonicalization|=" << time << "\n";
#endif
    return result;
  }
  Face FlipFaceIneq(std::pair<Face,MyVector<T>> const& pair, [[maybe_unused]] std::ostream & os) const {
#ifdef TIMINGS
    MicrosecondTime time;
#endif
    Face eFlip = FF.FlipFaceIneq(pair);
#ifdef TIMINGS
    os << "|FlipFaceIneq|=" << time << "\n";
#endif
    Face result = CanonicalImageGeneralDualDesc(can_method, GRP, recConvert, eFlip);
#ifdef TIMINGS
    os << "|canonicalization|=" << time << "\n";
#endif
    return result;
  }
};

template <typename T, typename Tgroup> struct DataFacetRepr {
  using Tint = typename Tgroup::Tint;
  size_t SelectedOrbit;
  Face eInc;
  FlippingFramework<T> FF;
  const Tgroup &GRP;
  Tgroup Stab;
  // There is only one method for Repr and it does not create a stabilizer.
  Face FlipFace(const Face &f, [[maybe_unused]] std::ostream & os) const {
#ifdef TIMINGS
    MicrosecondTime time;
#endif
    Face result = FF.FlipFace(f);
#ifdef TIMINGS
    os << "|FlipFace|=" << time << "\n";
#endif
    return result;
  }
  Face FlipFaceIneq(std::pair<Face,MyVector<T>> const& pair, [[maybe_unused]] std::ostream & os) const {
#ifdef TIMINGS
    MicrosecondTime time;
#endif
    Face result = FF.FlipFaceIneq(pair);
#ifdef TIMINGS
    os << "|FlipFaceIneq|=" << time << "\n";
#endif
    return result;
  }
};



template <typename Tgroup, typename Torbsize, typename Tidx>
struct FaceOrbsizeContainer {
public:
  using Tint = typename Tgroup::Tint;
  size_t n_act;
  size_t n_bit_orbsize;
  size_t delta;
  Tint TotalNumber;
  size_t nbOrbitDone;
  Tint nbUndone;
  size_t nbOrbit;
  std::vector<uint8_t> Vappend;
  DataFaceOrbitSize<Torbsize, Tgroup> recConvert;
  // We CANNOT replace ListOrbit by vectface as we use a number of hacks that
  // would not be available with a vectface.
  std::vector<uint8_t> ListOrbit;
  FaceOrbsizeContainer() = delete;
  FaceOrbsizeContainer(const FaceOrbsizeContainer &) = delete;
  FaceOrbsizeContainer &operator=(const FaceOrbsizeContainer &) = delete;
  FaceOrbsizeContainer(FaceOrbsizeContainer &&) = delete;
  FaceOrbsizeContainer(const Tgroup &GRP) : recConvert(GRP) {
    std::map<Tidx, int> LFact = GRP.factor_size();
    n_act = GRP.n_act();
    TotalNumber = 0;
    nbOrbitDone = 0;
    nbUndone = 0;
    nbOrbit = 0;
    std::pair<size_t, size_t> ep = get_delta(LFact, n_act);
    n_bit_orbsize = ep.first;
    delta = ep.second;
    Vappend = std::vector<uint8_t>((delta + 7) / 8, 0);
  }
  void clear() {
    TotalNumber = 0;
    nbOrbitDone = 0;
    nbUndone = 0;
    nbOrbit = 0;
    ListOrbit.clear();
  }
  // conversion functions that depend only on n_act and n_bit_orbsize.
  std::pair<Face,Tint> FaceToPair(Face const &f_in) const {
    Face f(n_act);
    Torbsize idx_orb = 0;
    for (size_t i = 0; i < n_act; i++)
      f[i] = f_in[i];
    size_t i_acc = n_act;
    Torbsize pow = 1;
    for (size_t i = 0; i < n_bit_orbsize; i++) {
      idx_orb += Torbsize(f_in[i_acc]) * pow;
      i_acc++;
      pow *= 2;
    }
    return {std::move(f), recConvert.ListPossOrbsize[idx_orb]};
  }
  // Extracting a block of faces for test cases
  vectface ExtractFirstNFace(size_t const& siz) const {
    vectface vf(n_act);
    Face f(n_act);
    for (size_t u=0; u<siz; u++) {
      size_t i_acc = delta * u;
      for (size_t i = 0; i < n_act; i++) {
        f[i] = getbit_vector(ListOrbit, i_acc);
        i_acc++;
      }
      vf.push_back(f);
    }
    return vf;
  }
  // Obtain a block of test faces that can be used for runtime estimation
  vectface get_initial_testface(bool const& use_database, size_t const& siz) const {
    if (use_database) {
      size_t siz_red = T_min(siz, nbOrbit);
      return ExtractFirstNFace(siz_red);
    } else {
      vectface vf(n_act);
      for (size_t iter=0; iter<siz; iter++) {
        Face f = RandomFace(n_act);
        vf.push_back(f);
      }
      return vf;
    }
  }
  // Database code that uses ListOrbit;
  std::pair<Face,Tint> RetrieveListOrbitEntry(size_t const &i_orb) const {
    Face f(n_act);
    Torbsize idx_orb = 0;
    size_t i_acc = delta * i_orb;
    for (size_t i = 0; i < n_act; i++) {
      f[i] = getbit_vector(ListOrbit, i_acc);
      i_acc++;
    }
    Torbsize pow = 1;
    for (size_t i = 0; i < n_bit_orbsize; i++) {
      idx_orb += Torbsize(getbit_vector(ListOrbit, i_acc)) * pow;
      i_acc++;
      pow *= 2;
    }
    return {std::move(f), recConvert.ListPossOrbsize[idx_orb]};
  }
  Face RetrieveListOrbitFace(size_t const &i_orb) const {
    Face face(n_act);
    size_t i_acc = delta * i_orb;
    for (size_t i = 0; i < n_act; i++) {
      face[i] = getbit_vector(ListOrbit, i_acc);
      i_acc++;
    }
    return face;
  }
  void InsertListOrbitEntry(Face const &f) {
    // Insert bytes to avoid a memory segfault.
    size_t curr_len = ListOrbit.size();
    size_t needed_bits = (nbOrbit + 1) * delta;
    size_t needed_len = (needed_bits + 7) / 8;
    size_t incr = needed_len - curr_len;
    if (incr > 0)
      ListOrbit.insert(ListOrbit.end(), Vappend.begin(),
                       Vappend.begin() + incr);
    // Now setting up the bits for face and idx_orb.
    size_t i_acc = nbOrbit * delta;
    for (size_t i = 0; i < n_act; i++) {
      bool val = f[i];
      setbit_vector(ListOrbit, i_acc, val);
      i_acc++;
    }
  }
  void InsertListOrbitFace(Face const &face) {
    size_t curr_len = ListOrbit.size();
    size_t needed_bits = (nbOrbit + 1) * delta;
    size_t needed_len = (needed_bits + 7) / 8;
    size_t incr = needed_len - curr_len;
    if (incr > 0)
      ListOrbit.insert(ListOrbit.end(), Vappend.begin(),
                       Vappend.begin() + incr);
    // Now setting up the bits but only for the faces as this suffices for the
    // comparison of novelty.
    size_t i_acc = nbOrbit * delta;
    for (size_t i = 0; i < n_act; i++) {
      bool val = face[i];
      setbit_vector(ListOrbit, i_acc, val);
      i_acc++;
    }
  }
  void InsertListOrbitFaceComplete(Face const &face) {
    size_t curr_len = ListOrbit.size();
    size_t needed_bits = (nbOrbit + 1) * delta;
    size_t needed_len = (needed_bits + 7) / 8;
    size_t incr = needed_len - curr_len;
    if (incr > 0)
      ListOrbit.insert(ListOrbit.end(), Vappend.begin(),
                       Vappend.begin() + incr);
    size_t i_acc = nbOrbit * delta;
    for (size_t i = 0; i < delta; i++) {
      bool val = face[i];
      setbit_vector(ListOrbit, i_acc, val);
      i_acc++;
    }
  }
  void InsertListOrbitIdxOrb(Tint const &orbSize) {
    Torbsize idx_orb = recConvert.GetOrbSizeIndex(orbSize);
    /* TRICK 8: The computation of the stabilizer is needed for getting the
       orbitsize but this is expensive to do. Therefore we first insert the list
       of faces and if found to be new then we insert afterwards the idx_orb */
    size_t i_acc = nbOrbit * delta + n_act;
    Torbsize work_idx = idx_orb;
    for (size_t i = 0; i < n_bit_orbsize; i++) {
      bool val = work_idx % 2;
      setbit_vector(ListOrbit, i_acc, val);
      i_acc++;
      work_idx = work_idx / 2;
    }
  }
  void Counts_InsertOrbit(const bool &status, const Tint &orbSize) {
    TotalNumber += orbSize;
    if (status) {
      nbOrbitDone++;
    } else {
      nbUndone += orbSize;
    }
    nbOrbit++;
  }
  void Counts_SetOrbitDone(const Tint &orbSize) {
    nbUndone -= orbSize;
    nbOrbitDone++;
  }
  FaceOrbitsizeTableContainer<Tint> GetListFaceOrbitsize() {
    vectface vfo;
    vfo.build_vectface(delta, nbOrbit, std::move(ListOrbit));
    std::vector<Tint> ListPoss = recConvert.ListPossOrbsize;
    return FaceOrbitsizeTableContainer(std::move(ListPoss), n_act, std::move(vfo));
  }
};

template <typename T, typename Tgroup>
vectface DirectComputationInitialFacetSet_Group(const MyMatrix<T> &EXT,
                                                const Tgroup &GRP,
                                                int const& can_method,
                                                const std::string &ansSamp,
                                                std::ostream &os) {
  // We can do a little better by passing a lambda to the
  // DirectComputationInitialFacetSet but that is a little overkill right now
  size_t nbRow = EXT.rows();
  vectface list_face(nbRow);
  for (auto &eFace : DirectComputationInitialFacetSet(EXT, ansSamp, os))
    list_face.push_back(CanonicalImageDualDesc(can_method, GRP, eFace));
  return list_face;
}

template <typename T_inp, typename Tint_inp, typename Tgroup_inp>
struct DatabaseCanonic {
public:
  // The number of orbits with 32 bits, is limited to 4294967296
  // which ought to be enough.
  using Tidx_orbit = uint32_t;
  using T = T_inp;
  using Tint = Tint_inp;
  using Tgroup = Tgroup_inp;
  using Telt = typename Tgroup::Telt;
  const MyMatrix<T> &EXT;
  const Tgroup &GRP;
  using Torbsize = uint16_t;
  using DataFacet = DataFacetCan<T, Torbsize, Tgroup>;
  using Tidx = typename Telt::Tidx;
  int nbRow;
  int nbCol;
  size_t delta;
  FaceOrbsizeContainer<Tgroup, Torbsize, Tidx> foc;
  int the_method; // we use an indifferent name because we also have it for DatabaseRepr
private:
  UNORD_SET<size_t, std::function<size_t(Tidx_orbit)>,
            std::function<bool(Tidx_orbit, Tidx_orbit)>>
      DictOrbit;
  std::map<size_t, std::vector<Tidx_orbit>> CompleteList_SetUndone;
  /* TRICK 3: Encoding the pair of face and idx_orb as bits allow us to save
   * memory */
#ifdef SUBSET_HASH
  std::vector<Tidx> subset_index;
  size_t n_bit_hash;
#endif
  size_t n_act;
  size_t n_act_div8;
#if defined MURMUR_HASH || defined ROBIN_HOOD_HASH
  std::vector<uint8_t> V_hash;
#endif
public:
  DatabaseCanonic() = delete;
  DatabaseCanonic(const DatabaseCanonic<T, Tint, Tgroup> &) = delete;
  DatabaseCanonic(DatabaseCanonic<T, Tint, Tgroup> &&) = delete;
  DatabaseCanonic &operator=(const DatabaseCanonic<T, Tint, Tgroup> &) = delete;

  void InsertEntryDatabase(std::pair<Face,Tint> const &face_pair, bool const &status,
                           size_t const &pos) {
    Face const& face = face_pair.first;
    Tint const& orbSize = face_pair.second;
    if (!status) {
      size_t len = face.count();
      CompleteList_SetUndone[len].push_back(pos);
    }
    foc.Counts_InsertOrbit(status, orbSize);
  }
  DatabaseCanonic(MyMatrix<T> const &_EXT, Tgroup const &_GRP)
      : EXT(_EXT), GRP(_GRP), foc(GRP) {
    the_method = std::numeric_limits<int>::max();

    /* TRICK 6: The UNORD_SET only the index and this saves in memory usage. */
    n_act = GRP.n_act();
    delta = foc.delta;
    n_act_div8 = (n_act + 7) / 8;
    nbRow = EXT.rows();
    nbCol = EXT.cols();
#if defined MURMUR_HASH || defined ROBIN_HOOD_HASH
    V_hash = std::vector<uint8_t>(n_act_div8, 0);
#endif
#ifdef SUBSET_HASH
    // The selection of indices for the hash
    size_t n_ent_bit = 8 * sizeof(size_t);
    if (n_act <= n_ent_bit) {
      n_bit_hash = n_act;
      for (size_t i = 0; i < n_ent_bit; i++)
        subset_index.push_back(Tidx(i));
    } else {
      n_bit_hash = n_ent_bit;
      double frac =
          static_cast<double>(n_act - 1) / static_cast<double>(n_ent_bit - 1);
      for (size_t i = 0; i < n_ent_bit; i++) {
        Tidx pos = Tidx(round(frac * static_cast<double>(i)));
        if (pos < 0)
          pos = 0;
        if (pos >= n_act)
          pos = n_act - 1;
        subset_index.push_back(pos);
      }
    }
#endif
    std::function<size_t(Tidx_orbit)> fctHash = [&](Tidx_orbit idx) -> size_t {
      size_t pos = delta * idx;
#if defined MURMUR_HASH || defined ROBIN_HOOD_HASH
      for (size_t i = 0; i < n_act; i++) {
        bool val = getbit_vector(foc.ListOrbit, pos);
        setbit_vector(V_hash, i, val);
        pos++;
      }
#ifdef MURMUR_HASH
      const uint32_t seed = 0x1b873560;
      return murmur3_32(V_hash.data(), n_act_div8, seed);
#endif
#ifdef ROBIN_HOOD_HASH
      const uint64_t seed = UINT64_C(0xe17a1465);
      size_t hash = robin_hood_hash_bytes(V_hash.data(), n_act_div8, seed);
      os << "hash=" << hash << "\n";
      return hash;
#endif
#endif
#ifdef SUBSET_HASH
      size_t hash = 0;
      size_t *ptr1 = &hash;
      uint8_t *ptr2 = (uint8_t *)ptr1;
      for (size_t i = 0; i < n_bit_hash; i++) {
        double idx = pos + size_t(subset_index[i]);
        bool val = getbit_vector(foc.ListOrbit, idx);
        setbit_ptr(ptr2, i, val);
      }
      return hash;
#endif
    };
    std::function<bool(Tidx_orbit, Tidx_orbit)> fctEqual = [&](Tidx_orbit idx1,
                                                               Tidx_orbit idx2) -> bool {
      size_t pos1 = delta * idx1;
      size_t pos2 = delta * idx2;
      for (size_t i = 1; i < n_act; i++) {
        // TRICK 9: Two faces will differ by at least 2 bits
        bool val1 = getbit_vector(foc.ListOrbit, pos1);
        bool val2 = getbit_vector(foc.ListOrbit, pos2);
        if (val1 != val2)
          return false;
        pos1++;
        pos2++;
      }
      return true;
    };
    DictOrbit =
        UNORD_SET<size_t, std::function<size_t(Tidx_orbit)>,
                  std::function<bool(Tidx_orbit, Tidx_orbit)>>({}, fctHash, fctEqual);
  }
  ~DatabaseCanonic() {}
  void clear() {
    foc.clear();
    DictOrbit.clear();
    CompleteList_SetUndone.clear();
  }
  bool use_f_insert_pair() {
    if (the_method == CANONIC_STRATEGY__CANONICAL_IMAGE) {
      return true;
    }
    return CANONIC_STRATEGY__INITIAL_TRIV;
    if (the_method == CANONIC_STRATEGY__STORE) {
      return true;
    }
    if (the_method == CANONIC_STRATEGY__INITIAL_TRIV) {
      return false;
    }
    std::cerr << "The value of the_method was not correctly set\n";
    throw TerminalException{1};
  }
  int get_default_strategy() {
    return CANONIC_STRATEGY__DEFAULT;
  }
  FaceOrbitsizeTableContainer<Tint> GetListFaceOrbitsize() {
    DictOrbit.clear();
    CompleteList_SetUndone.clear();
    return foc.GetListFaceOrbitsize();
  }
  void FuncInsert(Face const &face_can) {
    // The face should have been canonicalized beforehand.
    foc.InsertListOrbitFace(face_can);
    DictOrbit.insert(foc.nbOrbit);
    if (DictOrbit.size() == foc.nbOrbit) {
      // Insertion did not raise the count
      // and so it was already present
      return;
    }
    /* TRICK 8: The insertion yield something new. So now we compute the
     * expensive stabilizer */
    Tint orbSize = GRP.OrbitSize_OnSets(face_can);
    foc.InsertListOrbitIdxOrb(orbSize);
    InsertEntryDatabase({face_can, orbSize}, false, foc.nbOrbit);
  }
  void FuncInsertPair(Face const &face_orbsize) {
    // The face should have been canonicalized beforehand and also contains the orbits
    foc.InsertListOrbitFace(face_orbsize);
    DictOrbit.insert(foc.nbOrbit);
    if (DictOrbit.size() == foc.nbOrbit) {
      // Insertion did not raise the count
      // and so it was already present
      return;
    }
    std::pair<Face,Tint> pair = foc.FaceToPair(face_orbsize);
    InsertEntryDatabase(pair, false, foc.nbOrbit);
  }
  vectface ComputeInitialSet(const std::string &ansSamp, std::ostream &os) {
    return DirectComputationInitialFacetSet_Group(EXT, GRP, the_method, ansSamp, os);
  }
  void FuncPutOrbitAsDone(size_t const &i_orb) {
    std::pair<Face,Tint> eEnt = foc.RetrieveListOrbitEntry(i_orb);
    size_t len = eEnt.first.count();
    /* TRICK 1: We copy the last element in first position to erase it and then
     * pop_back the vector. */
    std::vector<Tidx_orbit> &V = CompleteList_SetUndone[len];
    if (V.size() == 1) {
      CompleteList_SetUndone.erase(len);
    } else {
      V[0] = V[V.size() - 1];
      V.pop_back();
      if (2*V.size() < V.capacity()) {
        V.shrink_to_fit();
      }
    }
    foc.Counts_SetOrbitDone(eEnt.second);
  }
  DataFacetCan<T, Torbsize, Tgroup> FuncGetMinimalUndoneOrbit() {
    for (auto &eEnt : CompleteList_SetUndone) {
      size_t len = eEnt.second.size();
      if (len > 0) {
        /* TRICK 1: Take the first element in the vector. This first element
           will remain
           in place but the vector may be extended without impacting this first
           entry. */
        size_t pos = eEnt.second[0];
        Face f = foc.RetrieveListOrbitFace(pos);
        Tgroup Stab = GRP.Stabilizer_OnSets(f);
        return {pos, f, FlippingFramework<T>(EXT, f), GRP, foc.recConvert,
                ReducedGroupAction(Stab, f), the_method};
      }
    }
    std::cerr << "Failed to find an undone orbit\n";
    throw TerminalException{1};
  }
  void InsertListOrbitEntry(Face const &eEnt, const size_t &i_orbit) {
    foc.InsertListOrbitEntry(eEnt);
    DictOrbit.insert(i_orbit);
  }

private:
  struct IteratorIndexType {
  private:
    std::map<size_t, std::vector<Tidx_orbit>>::const_iterator iter;
    size_t pos;

  public:
    IteratorIndexType(std::map<size_t, std::vector<Tidx_orbit>>::const_iterator iter,
                      size_t pos)
        : iter(iter), pos(pos) {}
    size_t operator*() const {
      return iter->second[pos];
    }
    IteratorIndexType &operator++() {
      pos++;
      if (pos == iter->second.size()) {
        iter++;
        pos = 0;
      }
      return *this;
    }
    IteratorIndexType operator++(int) {
      IteratorIndexType tmp = *this;
      pos++;
      if (pos == iter->second.size()) {
        iter++;
        pos = 0;
      }
      return tmp;
    }
    IteratorIndexType &operator--() {
      if (pos == 0) {
        iter--;
        pos = iter->second.size();
      }
      pos--;
      return *this;
    }
    IteratorIndexType operator--(int) {
      IteratorIndexType tmp = *this;
      if (pos == 0) {
        iter--;
        pos = iter->second.size();
      }
      pos--;
      return tmp;
    }

    bool operator!=(const IteratorIndexType &x) const {
      return pos != x.pos || iter != x.iter;
    }
    bool operator==(const IteratorIndexType &x) const {
      return pos == x.pos && iter == x.iter;
    }
  };

  struct IteratorFaceType {
  private:
    const FaceOrbsizeContainer<Tgroup, Torbsize, Tidx> &foc;
    IteratorIndexType iter;

  public:
    IteratorFaceType(const FaceOrbsizeContainer<Tgroup, Torbsize, Tidx> &foc,
                     IteratorIndexType iter)
        : foc(foc), iter(iter) {}
    Face operator*() const {
      return foc.RetrieveListOrbitFace(*iter);
    }
    IteratorFaceType &operator++() {
      iter++;
      return *this;
    }
    IteratorFaceType operator++(int) {
      IteratorFaceType tmp = *this;
      iter++;
      return tmp;
    }
    IteratorFaceType &operator--() {
      iter--;
      return *this;
    }
    IteratorFaceType operator--(int) {
      IteratorFaceType tmp = *this;
      tmp--;
      return tmp;
    }

    bool operator!=(const IteratorFaceType &x) const {
      return iter != x.iter;
    }
    bool operator==(const IteratorFaceType &x) const {
      return iter == x.iter;
    }
  };

public:
  using iterator_index = IteratorIndexType;
  using iterator_face = IteratorFaceType;
  iterator_index begin_index_undone() const {
    return IteratorIndexType(CompleteList_SetUndone.begin(), 0);
  }
  iterator_index end_index_undone() const {
    return IteratorIndexType(CompleteList_SetUndone.end(), 0);
  }
  iterator_face begin_face_undone() const {
    return IteratorFaceType(foc, begin_index_undone());
  }
  iterator_face end_face_undone() const {
    return IteratorFaceType(foc, end_index_undone());
  }
};

template <typename T_inp, typename Tint_inp, typename Tgroup_inp,
          typename Frepr, typename Forbitsize, typename Finv>
struct DatabaseRepr {
public:
  using T = T_inp;
  using Tint = Tint_inp;
  using Tgroup = Tgroup_inp;
  using Telt = typename Tgroup::Telt;
  const MyMatrix<T> &EXT;
  const Tgroup &GRP;
  using DataFacet = DataFacetRepr<T, Tgroup>;
  using Torbsize = uint16_t;
  using Tidx = typename Telt::Tidx;
  int nbRow;
  int nbCol;
  size_t delta;
  FaceOrbsizeContainer<Tgroup, Torbsize, Tidx> foc;
  int the_method;

private:
  std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>>
      CompleteList_SetUndone;
  std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>> CompleteList_SetDone;
  size_t n_act;
  Frepr f_repr;
  Forbitsize f_orbitsize;
  Finv f_inv;

public:
  DatabaseRepr() = delete;
  DatabaseRepr(const DatabaseRepr<T, Tint, Tgroup, Frepr, Forbitsize, Finv> &) =
      delete;
  DatabaseRepr(DatabaseRepr<T, Tint, Tgroup, Frepr, Forbitsize, Finv> &&) = delete;
  DatabaseRepr &
  operator=(const DatabaseRepr<T, Tint, Tgroup, Frepr, Forbitsize, Finv> &) = delete;

  void InsertEntryDatabase(std::pair<Face,Tint> const &face_pair, bool const &status,
                           size_t const &pos) {
    Face const& face = face_pair.first;
    Tint const& orbSize = face_pair.second;
    size_t len = face.count();
    size_t eInv = f_inv(face);
    if (status) {
      CompleteList_SetDone[len][eInv].push_back(pos);
    } else {
      CompleteList_SetUndone[len][eInv].push_back(pos);
    }
    foc.Counts_InsertOrbit(status, orbSize);
  }
  DatabaseRepr(MyMatrix<T> const &_EXT, Tgroup const &_GRP, Frepr f_repr,
               Forbitsize f_orbitsize, Finv f_inv)
      : EXT(_EXT), GRP(_GRP), foc(GRP),
        f_repr(f_repr), f_orbitsize(f_orbitsize), f_inv(f_inv) {
    /* TRICK 6: The UNORD_SET only the index and this saves in memory usage. */
    n_act = GRP.n_act();
    delta = foc.delta;
    nbRow = EXT.rows();
    nbCol = EXT.cols();
  }
  ~DatabaseRepr() {}
  void clear() {
    foc.clear();
    CompleteList_SetUndone.clear();
    CompleteList_SetDone.clear();
  }
  bool use_f_insert_pair() {
    return false;
  }
  int get_default_strategy() {
    return REPR_STRATEGY__DEFAULT;
  }
  FaceOrbitsizeTableContainer<Tint> GetListFaceOrbitsize() {
    CompleteList_SetUndone.clear();
    CompleteList_SetDone.clear();
    return foc.GetListFaceOrbitsize();
  }
  void FuncInsert(Face const &face_i) {
    size_t len = face_i.count();
    size_t eInv = f_inv(face_i);
    if (CompleteList_SetDone.count(len) == 1) {
      if (CompleteList_SetDone[len].count(eInv) == 1) {
        for (size_t &i_orb : CompleteList_SetDone[len][eInv]) {
          Face face_e = foc.RetrieveListOrbitFace(i_orb);
          bool test = f_repr(face_i, face_e);
          if (test)
            return;
        }
      }
    }
    if (CompleteList_SetUndone.count(len) == 1) {
      if (CompleteList_SetUndone[len].count(eInv) == 1) {
        for (size_t &i_orb : CompleteList_SetUndone[len][eInv]) {
          Face face_e = foc.RetrieveListOrbitFace(i_orb);
          bool test = f_repr(face_i, face_e);
          if (test)
            return;
        }
      }
    }
    foc.InsertListOrbitFace(face_i);
    // We need to recompute
    Tint orbSize = f_orbitsize(face_i);
    foc.InsertListOrbitIdxOrb(orbSize);
    InsertEntryDatabase({face_i,orbSize}, false, foc.nbOrbit);
  }
  void FuncInsertPair(Face const &face) {
    Face f_red(nbRow);
    for (int i=0; i<nbRow; i++) {
      f_red[i] = face[i];
    }
    FuncInsert(f_red);
  }
  vectface ComputeInitialSet(const std::string &ansSamp, std::ostream &os) {
    return DirectComputationInitialFacetSet(EXT, ansSamp, os);
  }
  void FuncPutOrbitAsDone(size_t const &iOrb) {
    std::pair<Face,Tint> eEnt = foc.RetrieveListOrbitEntry(iOrb);
    size_t len = eEnt.first.count();
    /* TRICK 1: We copy the last element in first position to erase it and then
     * pop_back the vector. */
    size_t eInv = f_inv(eEnt.first);
    std::vector<size_t> &V = CompleteList_SetUndone[len][eInv];
    size_t spec_orbit = V[0];
    V[0] = V[V.size() - 1];
    V.pop_back();
    if (V.size() == 0) {
      CompleteList_SetUndone[len].erase(eInv);
    }
    if (CompleteList_SetUndone[len].size() == 0) {
      CompleteList_SetUndone.erase(len);
    }
    CompleteList_SetDone[len][eInv].push_back(spec_orbit);
    foc.Counts_SetOrbitDone(eEnt.second);
  }
  DataFacetRepr<T, Tgroup> FuncGetMinimalUndoneOrbit() {
    auto iter1 = CompleteList_SetUndone.begin();
    if (iter1 == CompleteList_SetUndone.end()) {
      std::cerr << "We have an empty map (1) We should not reach that stage\n";
      throw TerminalException{1};
    }
    UNORD_MAP<size_t, std::vector<size_t>> &map_inv = iter1->second;
    auto iter2 = map_inv.begin();
    if (iter2 == map_inv.end()) {
      std::cerr << "We have an empty map (2). We should not reach that stage\n";
      throw TerminalException{1};
    }
    const std::vector<size_t> &V = iter2->second;
    /* TRICK 1: Take the first element in the vector. This first element will
       remain in place but the vector may be extended without impacting this
       first entry. */
    size_t pos = V[0];
    Face f = foc.RetrieveListOrbitFace(pos);
    Tgroup Stab = GRP.Stabilizer_OnSets(f);
    return {pos, f, FlippingFramework<T>(EXT, f), GRP,
            ReducedGroupAction(Stab, f)};
  }
  void InsertListOrbitEntry(Face const &eEnt,
                            [[maybe_unused]] const size_t &i_orbit) {
    foc.InsertListOrbitEntry(eEnt);
  }

private:
  struct IteratorIndexType {
  private:
    std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>>::const_iterator
        iter1;
    std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>>::const_iterator
        iter1_end;
    UNORD_MAP<size_t, std::vector<size_t>>::const_iterator iter2;
    size_t pos;

  public:
    IteratorIndexType(
        std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>>::const_iterator
            iter1,
        std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>>::const_iterator
            iter1_end,
        UNORD_MAP<size_t, std::vector<size_t>>::const_iterator iter2,
        size_t pos)
        : iter1(iter1), iter1_end(iter1_end), iter2(iter2), pos(pos) {
    }
    size_t operator*() const {
      return iter2->second[pos];
    }
    void PtrIncrease() {
      pos++;
      if (pos == iter2->second.size()) {
        pos = 0;
        iter2++;
        if (iter2 == iter1->second.end()) {
          iter1++;
          if (iter1 != iter1_end) {
            iter2 = iter1->second.begin();
          }
        }
      }
    }
    IteratorIndexType &operator++() {
      PtrIncrease();
      return *this;
    }
    IteratorIndexType operator++(int) {
      IteratorIndexType tmp = *this;
      PtrIncrease();
      return tmp;
    }
    bool operator!=(const IteratorIndexType &x) const {
      // If one of iter1 operator is the end then the comparison of the first
      // one suffices to conclude
      if (x.iter1 == iter1_end)
        return iter1 != iter1_end;
      if (iter1 == iter1_end)
        return x.iter1 != iter1_end;
      //
      return x.iter1 != iter1 || x.iter2 != iter2 || x.pos != pos;
    }
  };

  struct IteratorFaceType {
  private:
    const FaceOrbsizeContainer<Tgroup, Torbsize, Tidx> &foc;
    IteratorIndexType iter;

  public:
    IteratorFaceType(const FaceOrbsizeContainer<Tgroup, Torbsize, Tidx> &foc,
                     IteratorIndexType iter)
        : foc(foc), iter(iter) {}
    Face operator*() const {
      return foc.RetrieveListOrbitFace(*iter);
    }
    IteratorFaceType &operator++() {
      iter++;
      return *this;
    }
    IteratorFaceType operator++(int) {
      IteratorFaceType tmp = *this;
      iter++;
      return tmp;
    }
    IteratorFaceType &operator--() {
      iter--;
      return *this;
    }
    IteratorFaceType operator--(int) {
      IteratorFaceType tmp = *this;
      tmp--;
      return tmp;
    }
    bool operator!=(const IteratorFaceType &x) const {
      return iter != x.iter;
    }
    bool operator==(const IteratorFaceType &x) const {
      return iter == x.iter;
    }
  };

public:
  using iterator_index = IteratorIndexType;
  using iterator_face = IteratorFaceType;
  iterator_index begin_index_undone() const {
    std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>>::const_iterator
        iter1 = CompleteList_SetUndone.begin();
    std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>>::const_iterator
        iter1_end = CompleteList_SetUndone.end();
    if (iter1 == iter1_end)
      return IteratorIndexType(iter1, iter1_end, {}, 0);
    UNORD_MAP<size_t, std::vector<size_t>>::const_iterator iter2 =
        CompleteList_SetUndone.at(iter1->first).begin();
    return IteratorIndexType(iter1, iter1_end, iter2, 0);
  }
  iterator_index end_index_undone() const {
    return IteratorIndexType(CompleteList_SetUndone.end(),
                             CompleteList_SetUndone.end(), {}, 0);
  }
  iterator_face begin_face_undone() const {
    return IteratorFaceType(foc, begin_index_undone());
  }
  iterator_face end_face_undone() const {
    return IteratorFaceType(foc, end_index_undone());
  }
};

template <typename TbasicBank> struct DatabaseOrbits {
public:
  using Tgroup = typename TbasicBank::Tgroup;
  using T = typename TbasicBank::T;
  using Telt = typename Tgroup::Telt;
  using Tint = typename TbasicBank::Tint;
  Tint CritSiz;
  TbasicBank &bb;

private:
  std::string MainPrefix;
  std::string eFileEXT, eFileGRP, eFileNB, eFileFB, eFileFF, eFileMethod1, eFileMethod2;
  /* TRICK 7: Using separate files for faces and status allow us to gain
     locality. The faces are written one by one while the access to status is
     random */
  bool SavingTrigger;
  bool NeedToFlush;
  bool AdvancedTerminationCriterion;
  std::ostream &os;
  size_t delta;
  std::string strPresChar;

public:
  // method1 is the preceding method computed, method2 another one.
  // ---If method1 and method2 are missing then we compute a best guess and
  //    then we update the database. We can accelerate the method to identify
  //    if at some point we identify that the guess is as good as what we have
  //    then this can accelerate.
  // ---If method1 is present but not method2 then we sample and get the best.
  //    If method1 = method2 then nothing to do. Otherwise, we recompute
  //    everything.
  // ---If method1 and method2 are present then nothing to be done.
  int the_method1, the_method2;
  int read_method(std::string const& eFileMethod) const {
    if (IsExistingFile(eFileMethod)) {
      return std::numeric_limits<int>::max();
    } else {
      std::ifstream is(eFileMethod);
      int method;
      is >> method;
      return method;
    }
  }
  void write_method(std::string const& eFileMethod, int const& method) const {
    std::ofstream os(eFileMethod);
    os << method;
  }
  bool need_recompute_representative() const {
    if (the_method2 == std::numeric_limits<int>::max())
      return true;
    return false;
  }
  bool is_database_present() const {
    return IsExistingFile(eFileEXT);
  }
  DatabaseOrbits() = delete;
  DatabaseOrbits(const DatabaseOrbits<TbasicBank> &) = delete;
  DatabaseOrbits(DatabaseOrbits<TbasicBank> &&) = delete;
  DatabaseOrbits &operator=(const DatabaseOrbits<TbasicBank> &) = delete;
  void print_status() const {
    os << "Status : orbit=(" << bb.foc.nbOrbit << "," << bb.foc.nbOrbitDone
       << "," << (bb.foc.nbOrbit - bb.foc.nbOrbitDone) << ") facet=("
       << bb.foc.TotalNumber << "," << (bb.foc.TotalNumber - bb.foc.nbUndone)
       << "," << bb.foc.nbUndone << ")\n\n";
  }
  DatabaseOrbits(TbasicBank &bb, const std::string &MainPrefix,
                 const bool &_SavingTrigger,
                 const bool &_AdvancedTerminationCriterion, std::ostream &os)
      : CritSiz(bb.EXT.cols() - 2), bb(bb), SavingTrigger(_SavingTrigger),
        AdvancedTerminationCriterion(_AdvancedTerminationCriterion), os(os) {
    os << "MainPrefix=" << MainPrefix << "\n";
    eFileEXT = MainPrefix + ".ext";
    eFileGRP = MainPrefix + ".grp";
    eFileNB = MainPrefix + ".nb";
    eFileFB = MainPrefix + ".fb";
    eFileFF = MainPrefix + ".ff";
    eFileMethod1 = MainPrefix + ".method1";
    eFileMethod2 = MainPrefix + ".method2";
    the_method1 = read_method(eFileMethod1);
    if (the_method1 == std::numeric_limits<int>::max()) {
      // Put the default if missing
      the_method1 = bb.get_default_strategy();
    }
    the_method2 = read_method(eFileMethod2);
    strPresChar = "|EXT|=" + std::to_string(bb.nbRow) + "/" +
                  std::to_string(bb.nbCol) +
                  " |GRP|=" + std::to_string(bb.GRP.size());
    delta = bb.delta;
    NeedToFlush = true;
    // At the start we are working with a default strategy.
    bb.the_method = bb.get_default_strategy();
    if (SavingTrigger) {
      if (is_database_present()) {
        LoadDatabase();
      } else {
        if (!FILE_IsFileMakeable(eFileEXT)) {
          os << "Error in DatabaseOrbits: File eFileEXT=" << eFileEXT
             << " is not makeable\n";
          throw TerminalException{1};
        }
        os << "Creating the files (NB, FB, FF)\n";
        // Writing Group
        os << "eFileGRP=" << eFileGRP << "\n";
        std::ofstream os_grp(eFileGRP);
        os_grp << bb.GRP;
        WriteMatrixFile(eFileEXT, bb.EXT);
      }
      //
      print_status();
    }
  }
  size_t prelaod_nb_orbit() {
    if (SavingTrigger) {
      if (is_database_present()) {
        FileNumber fn(eFileNB, false);
        return fn.getval();
      }
    }
    return 0;
  }
  void LoadDatabase() {
    os << "Opening existing files (NB, FB, FF)\n";
#ifdef TIMINGS
    MicrosecondTime time;
#endif
    FileNumber fn(eFileNB, false);
    size_t n_orbit = fn.getval();
    os << "Loading database with n_orbit=" << n_orbit << "\n";
    FileBool fb(eFileFB, n_orbit);
    FileFace ff(eFileFF, bb.delta, n_orbit);
    for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
      Face f = ff.getface(i_orbit);
      std::pair<Face,Tint> eEnt = bb.foc.FaceToPair(f);
      bool status = fb.getbit(i_orbit);
      bb.InsertListOrbitEntry(f, i_orbit);
      bb.InsertEntryDatabase(eEnt, status, i_orbit);
    }
#ifdef TIMINGS
    os << "|Loading Database|=" << time << "\n";
#endif
  }
  vectface ReadDatabase() {
    os << "Opening existing files (NB, FB, FF)\n";
#ifdef TIMINGS
    MicrosecondTime time;
#endif
    FileNumber fn(eFileNB, false);
    size_t n_orbit = fn.getval();
    os << "Reading database with n_orbit=" << n_orbit << "\n";
    FileBool fb(eFileFB, n_orbit);
    FileFace ff(eFileFF, bb.delta, n_orbit);
    vectface vf(bb.delta + 1);
    for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
      Face f = ff.getface(i_orbit);
      bool status = fb.getbit(i_orbit);
      Face f_insert(bb.delta + 1);
      for (size_t u=0; u<bb.delta; u++) {
        f_insert[u] = f[u];
      }
      f_insert[bb.delta] = status;
      vf.push_back(f_insert);
    }
#ifdef TIMINGS
    os << "|Reading Database|=" << time << "\n";
#endif
    return vf;
  }
  void DirectAppendDatabase(vectface && vf) {
    bb.clear();
    size_t n_orbit = vf.size();
    size_t len_ff = 0;
    size_t len_fb = 0;
    if (SavingTrigger) {
      len_ff = (n_orbit * bb.delta + 7) / 8;
      len_fb = (n_orbit + 7) / 8;
    }
    std::vector<uint8_t> ListOrbit_ff(len_ff);
    std::vector<uint8_t> V_status(len_fb);
    size_t pos_ff = 0;
    for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++) {
      Face f = vf[i_orbit];
      Face f_red(bb.delta);
      for (size_t u=0; u<bb.delta; u++) {
        bool val = f[u];
        f_red[u] = val;
        if (SavingTrigger) {
          setbit_vector(ListOrbit_ff, pos_ff, val);
        }
        pos_ff++;
      }
      bool status = f[bb.delta];
      if (SavingTrigger) {
        setbit_vector(V_status, i_orbit, status);
      }
      std::pair<Face,Tint> eEnt = bb.foc.FaceToPair(f);
      bb.InsertListOrbitEntry(eEnt, i_orbit);
      bb.InsertEntryDatabase(eEnt, status, i_orbit);
    }
    if (SavingTrigger) {
      FileNumber fn(eFileNB, true);
      FileBool fb(eFileFB);
      FileFace ff(eFileFF, bb.delta);
      fn.setval(n_orbit);
      ff.direct_write(ListOrbit_ff);
      fb.direct_write(V_status);
    }
  }
  void DatabaseMethodUpgrade(int const& the_method) {
    FileNumber fn(eFileNB, false);
    size_t n_orbit = fn.getval();
    FileBool fb(eFileFB, n_orbit);
    FileFace ff(eFileFF, bb.delta, n_orbit);
    for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
      Face f = ff.getface(i_orbit);
      std::pair<Face,Tint> eEnt = bb.foc.FaceToPair(f);
      bool status = fb.getbit(i_orbit);
      bb.InsertListOrbitEntry(eEnt, i_orbit);
      bb.InsertEntryDatabase(eEnt, status, i_orbit);
    }
  }
  ~DatabaseOrbits() {
    /* TRICK 5: The destructor does NOT destroy the database! This is because it
       can be used in another call. Note that the returning of the list of orbit
       does destroy the database and this gives a small window in which bad
       stuff can happen.
     */
    if (SavingTrigger && NeedToFlush) {
      flush();
    }
    os << "Clean closing of the DatabaseOrbits\n";
  }
  void flush() const {
    std::cerr << "Doing the flushing operation\n";
#ifdef TIMINGS
    MicrosecondTime time;
#endif
    FileNumber fn(eFileNB, true);
    FileBool fb(eFileFB);
    FileFace ff(eFileFF, bb.delta);
    ff.direct_write(bb.foc.ListOrbit);
    size_t nbOrbit = bb.foc.nbOrbit;
    fn.setval(nbOrbit);
    size_t len = (nbOrbit + 7) / 8;
    std::vector<uint8_t> V_status(len, 255);
    auto iter = bb.begin_index_undone();
    while (iter != bb.end_index_undone()) {
      size_t pos = *iter;
      setbit_vector(V_status, pos, false);
      iter++;
    }
    fb.direct_write(V_status);
#ifdef TIMINGS
    os << "|flush|=" << time << "\n";
#endif
  }
  // FuncListOrbitIncidence() {
  FaceOrbitsizeTableContainer<Tint> GetListFaceOrbitsize() {
    NeedToFlush = false;
    if (SavingTrigger) {
      RemoveFileIfExist(eFileNB);
      RemoveFileIfExist(eFileFB);
      RemoveFileIfExist(eFileFF);
      RemoveFileIfExist(eFileEXT);
      RemoveFileIfExist(eFileGRP);
      RemoveFileIfExist(eFileMethod1);
      RemoveFileIfExist(eFileMethod2);
    }
    return bb.GetListFaceOrbitsize();
  }
  void FuncInsert(Face const &face) {
    bb.FuncInsert(face);
  }
  void FuncInsertPair(Face const &face) {
    bb.FuncInsertPair(face);
  }
  vectface ComputeInitialSet(const std::string &ansSamp, std::ostream &os) {
    return bb.ComputeInitialSet(ansSamp, os);
  }
  void FuncPutOrbitAsDone(size_t const &i_orb) {
    bb.FuncPutOrbitAsDone(i_orb);
    print_status();
  }
  Face ComputeIntersectionUndone() const {
    size_t n_row = bb.EXT.rows();
    Face eSetReturn(n_row);

    // don't do full computation if many orbit remaining
    // for some polytopes only the last orbit sets eSetReturn = 0
    // resulting in large slowdowns here
    // alternative fix: enumerate in decending order
    if (bb.foc.nbOrbit - bb.foc.nbOrbitDone > 1000)
      return eSetReturn;

    for (size_t i_row = 0; i_row < n_row; i_row++)
      eSetReturn[i_row] = 1;
    typename TbasicBank::iterator_face iter = bb.begin_face_undone();
    while (iter != bb.end_face_undone()) {
      eSetReturn &= OrbitIntersection(bb.GRP, *iter);
      if (eSetReturn.count() == 0) {
        return eSetReturn;
      }
      iter++;
    }
    return eSetReturn;
  }
  size_t FuncNumberOrbit() const { return bb.foc.nbOrbit; }
  bool IsFinished() const { return bb.foc.nbOrbit == bb.foc.nbOrbitDone; }
  typename TbasicBank::DataFacet FuncGetMinimalUndoneOrbit() {
    typename TbasicBank::DataFacet data = bb.FuncGetMinimalUndoneOrbit();
    os << strPresChar << " Considering orbit " << data.SelectedOrbit
       << " |inc|=" << data.eInc.count() << " |stab|=" << data.Stab.size()
       << "\n";
    return data;
  }
  bool GetTerminationStatus() const {
    if (bb.foc.nbOrbitDone > 0) {
      Face eSetUndone = ComputeIntersectionUndone();
      if (bb.foc.nbUndone <= CritSiz || eSetUndone.count() > 0) {
        os << "End of computation, nbObj=" << bb.foc.TotalNumber
           << " nbUndone=" << bb.foc.nbUndone
           << " |eSetUndone|=" << eSetUndone.count() << " |EXT|=" << bb.nbRow
           << "\n";
        return true;
      }
    }
    if (AdvancedTerminationCriterion)
      return EvaluationConnectednessCriterion_Serial(bb, os);
    return false;
  }
  UndoneOrbitInfo<Tint> GetTerminationInfo() const {
    return {bb.foc.nbOrbitDone, bb.foc.nbUndone, ComputeIntersectionUndone()};
  }
};

template <typename T, typename Tidx_value> struct LazyWMat {
public:
  LazyWMat(const MyMatrix<T> &EXT) : EXT(EXT), HaveWMat(false) {}
  WeightMatrix<true, T, Tidx_value> &GetWMat() {
    if (HaveWMat)
      return WMat;
    WMat = GetWeightMatrix<T, Tidx_value>(EXT);
    WMat.ReorderingSetWeight();
    HaveWMat = true;
    return WMat;
  }

private:
  const MyMatrix<T> &EXT;
  WeightMatrix<true, T, Tidx_value> WMat;
  bool HaveWMat;
};

template <typename Tint, typename T, typename Tgroup>
std::map<std::string, Tint>
ComputeInitialMap(const MyMatrix<T> &EXT, const Tgroup &GRP,
                  PolyHeuristicSerial<typename Tgroup::Tint> &AllArr) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  std::map<std::string, Tint> TheMap;
  int delta = nbRow - nbCol;
  int level = AllArr.dimEXT - nbCol;
  TheMap["groupsize"] = GRP.size();
  TheMap["incidence"] = nbRow;
  TheMap["rank"] = nbCol;
  TheMap["delta"] = delta;
  TheMap["level"] = level;
  return TheMap;
}

template<typename Tgroup>
void CheckTermination(PolyHeuristicSerial<typename Tgroup::Tint> &AllArr) {
  if (ExitEvent) {
    std::cerr << "Terminating the program by Ctrl-C\n";
    throw TerminalException{1};
  }
  if (AllArr.max_runtime > 0) {
    int runtime = si(AllArr.start);
    if (runtime > AllArr.max_runtime) {
      std::cerr << "The maximum runtime has been elapsed. max_runtime = "
                << AllArr.max_runtime << "\n";
      throw RuntimeException{1};
    }
  }
}


// Needs initial definition due to template crazyness
template<typename Tbank, typename T, typename Tgroup, typename Tidx_value>
vectface DUALDESC_AdjacencyDecomposition(
    Tbank &TheBank, MyMatrix<T> const &EXT, Tgroup const &GRP,
    std::map<std::string, typename Tgroup::Tint> & TheMap,
    PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
    std::string const &ePrefix);




template<typename Tbank, typename T, typename Tgroup, typename Tidx_value,
         typename TbasicBank, typename Finsert>
void DUALDESC_AdjacencyDecomposition_and_insert(
    Tbank &TheBank, typename TbasicBank::DataFacet const& df,
    PolyHeuristicSerial<typename Tgroup::Tint> &AllArr, Finsert f_insert,
    std::string const &ePrefix, std::ostream& os) {
  using Tint = typename Tgroup::Tint;
  CheckTermination<Tgroup>(AllArr);
  std::map<std::string, Tint> TheMap =
    ComputeInitialMap<Tint>(df.FF.EXT_face, df.Stab, AllArr);
  std::string ansSplit = HeuristicEvaluation(TheMap, AllArr.Splitting);
  if (ansSplit != "split") {
    std::string ansProg = AllArr.DualDescriptionProgram.get_eval(TheMap);
    vectface TheOutput = DirectFacetOrbitComputation(df.FF.EXT_face, df.Stab, ansProg, os);
    AllArr.DualDescriptionProgram.pop(os);
#ifdef TIMINGS
    MicrosecondTime time_full;
    os << "|outputsize|=" << TheOutput.size() << "\n";
#endif
    for (auto &eOrb : TheOutput) {
      Face eFlip = df.FlipFace(eOrb, os);
#ifdef TIMINGS
      MicrosecondTime time;
#endif
      f_insert(eFlip);
#ifdef TIMINGS
      os << "|insert1|=" << time << "\n";
#endif
    }
#ifdef TIMINGS
    os << "|outputtime|=" << time_full << "\n";
#endif
  } else {
    vectface TheOutput =
      DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value>(
        TheBank, df.FF.EXT_face, df.Stab, TheMap, AllArr, ePrefix, os);
#ifdef TIMINGS
    MicrosecondTime time_full;
    os << "|outputsize|=" << TheOutput.size() << "\n";
#endif
    for (auto &eOrb : TheOutput) {
      Face eFlip = df.FlipFace(eOrb, os);
#ifdef TIMINGS
      MicrosecondTime time;
#endif
      f_insert(eFlip);
#ifdef TIMINGS
      os << "|insert2|=" << time << "\n";
#endif
    }
#ifdef TIMINGS
    os << "|outputtime|=" << time_full << "\n";
#endif
  }
}


template <typename Tbank, typename T, typename Tgroup, typename Tidx_value,
          typename TbasicBank>
FaceOrbitsizeTableContainer<typename Tgroup::Tint> Kernel_DUALDESC_AdjacencyDecomposition(
    Tbank &TheBank, TbasicBank &bb,
    PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
    std::string const &ePrefix,
    std::map<std::string, typename Tgroup::Tint> const &TheMap,
    std::ostream &os) {
  using DataFacet = typename TbasicBank::DataFacet;
  DatabaseOrbits<TbasicBank> RPL(bb, ePrefix, AllArr.Saving,
                                 AllArr.AdvancedTerminationCriterion, os);
  if (RPL.FuncNumberOrbit() == 0) {
    std::string ansSamp = HeuristicEvaluation(TheMap, AllArr.InitialFacetSet);
    for (auto &face : RPL.ComputeInitialSet(ansSamp, os))
      RPL.FuncInsert(face);
  }
  bool use_f_insert_pair = bb.use_f_insert_pair();
  while (true) {
    if (RPL.GetTerminationStatus())
      break;
    DataFacet df = RPL.FuncGetMinimalUndoneOrbit();
    size_t SelectedOrbit = df.SelectedOrbit;
    // Alternative way is to use CondTempDirectory. BUT
    // For many we actually do not need to have such a construction.
    // Need to think.
    std::string NewPrefix =
        ePrefix + "ADM" + std::to_string(SelectedOrbit) + "_";
    if (use_f_insert_pair) {
      auto f_insert=[&](Face const& eFlip) -> void {
        RPL.FuncInsertPair(eFlip);
      };
      DUALDESC_AdjacencyDecomposition_and_insert<Tbank,T,Tgroup,Tidx_value,TbasicBank,decltype(f_insert)>(TheBank, df, AllArr, f_insert, NewPrefix, os);
    } else {
      auto f_insert=[&](Face const& eFlip) -> void {
        RPL.FuncInsert(eFlip);
      };
      DUALDESC_AdjacencyDecomposition_and_insert<Tbank,T,Tgroup,Tidx_value,TbasicBank,decltype(f_insert)>(TheBank, df, AllArr, f_insert, NewPrefix, os);
    }
    RPL.FuncPutOrbitAsDone(SelectedOrbit);
  }
  return RPL.GetListFaceOrbitsize();
}




//
// A number of appoximations are done in this code:
// ---In the bank we assume that the full symmetry is used.
//    This means less things to store
// ---We use the canonicalization approach which allows to treat smaller cases.
// ---Serial mode. Should be faster indeed.
//
template <typename Tbank, typename T, typename Tgroup, typename Tidx_value>
vectface DUALDESC_AdjacencyDecomposition(Tbank &TheBank,
    MyMatrix<T> const &EXT, Tgroup const &GRP,
    std::map<std::string, typename Tgroup::Tint> & TheMap,
    PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
    std::string const &ePrefix, std::ostream &os) {
  using Tgr = GraphListAdj;
  using Tint = typename Tgroup::Tint;
  os << "Beginning of DUALDESC_AdjacencyDecomposition\n";
  CheckTermination<Tgroup>(AllArr);
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  LazyWMat<T, Tidx_value> lwm(EXT);
  //
  // Checking if the entry is present in the map.
  //
  std::string ansBankCheck =
      HeuristicEvaluation(TheMap, AllArr.CheckDatabaseBank);
  if (ansBankCheck == "yes") {
    vectface ListFace =
      getdualdesc_in_bank(TheBank, EXT, lwm.GetWMat(), GRP, AllArr, os);
    if (ListFace.size() > 0)
      return ListFace;
  }
  //
  // Now computing the groups
  //
  Tgroup TheGRPrelevant;
  //
  // The computations themselves
  //
  SingletonTime start;
  bool NeedSplit = false;
  // 3 scenarii
  // --- 1 : We have the full symmetry group and the computation was done with
  // respect to it.
  // --- 2 : We have computed for a subgroup which actually is the full group.
  // --- 3 : We have computed for a subgroup which actually is a strict
  // subgroup.
  bool BankSymmCheck;
  auto compute_decomposition = [&]() -> FaceOrbitsizeTableContainer<Tint> {
    std::string ansSymm =
      HeuristicEvaluation(TheMap, AllArr.AdditionalSymmetry);
    os << "ansSymm=" << ansSymm << "\n";
    if (ansSymm == "yes") {
      TheGRPrelevant = GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(lwm.GetWMat());
      NeedSplit = TheGRPrelevant.size() != GRP.size();
      BankSymmCheck = false;
    } else {
      TheGRPrelevant = GRP;
      BankSymmCheck = true;
    }
    Tint GroupSizeComp = TheGRPrelevant.size();
    os << "RESPAWN a new ADM computation |GRP|=" << GroupSizeComp
       << " TheDim=" << nbCol << " |EXT|=" << nbRow << "\n";
    std::string MainPrefix = ePrefix + "D_" + std::to_string(nbRow);
    std::string ansChosenDatabase =
      HeuristicEvaluation(TheMap, AllArr.ChosenDatabase);
    os << "DUALDESC_ChosenDatabase : ChosenDatabase = " << ansChosenDatabase
       << "\n";
    if (ansChosenDatabase == "canonic") {
      using TbasicBank = DatabaseCanonic<T, Tint, Tgroup>;
      TbasicBank bb(EXT, TheGRPrelevant);
      return Kernel_DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup,Tidx_value, TbasicBank>(
          TheBank, bb, AllArr, MainPrefix, TheMap, os);
    }
    if (ansChosenDatabase == "repr") {
      WeightMatrix<true, int, Tidx_value> WMat =
        WeightMatrixFromPairOrbits<Tgroup, Tidx_value>(TheGRPrelevant);
      auto f_repr = [&](const Face &f1, const Face &f2) -> bool {
        auto test = TheGRPrelevant.RepresentativeAction_OnSets(f1, f2);
        if (test)
          return true;
        return false;
      };
      auto f_orbitsize = [&](const Face &f) -> Tint {
        return TheGRPrelevant.OrbitSize_OnSets(f);
      };
      auto f_inv = [&](const Face &f) -> size_t {
        return GetLocalInvariantWeightMatrix(WMat, f);
      };
      using TbasicBank = DatabaseRepr<T, Tint, Tgroup, decltype(f_repr),
                                      decltype(f_orbitsize), decltype(f_inv)>;
      TbasicBank bb(EXT, TheGRPrelevant, f_repr, f_orbitsize, f_inv);
      return Kernel_DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value, TbasicBank>(
        TheBank, bb, AllArr, MainPrefix, TheMap, os);
    }
    std::cerr << "compute_split_or_not: Failed to find a matching entry\n";
    std::cerr << "Authorized values: canonic, repr\n";
    throw TerminalException{1};
  };
  FaceOrbitsizeTableContainer<Tint> ListOrbitFaceOrbitsize = compute_decomposition();
  SingletonTime end;
  TheMap["time"] = s(start, end);
  std::string ansBank = HeuristicEvaluation(TheMap, AllArr.BankSave);
  os << "elapsed_seconds=" << s(start, end) << " ansBank=" << ansBank
     << " NeedSplit=" << NeedSplit << "\n";
  if (ansBank == "yes") {
    os << "Before insert_entry_in_bank\n";
    insert_entry_in_bank(TheBank, EXT, lwm.GetWMat(), TheGRPrelevant,
                         BankSymmCheck, ListOrbitFaceOrbitsize);
  }
  os << "Before return section\n";
  if (NeedSplit) {
    os << "Before OrbitSplittingListOrbit\n";
    return OrbitSplittingListOrbitGen(TheGRPrelevant, GRP, ListOrbitFaceOrbitsize, AllArr, os);
  } else {
    os << "Returning ListOrbitFaces\n";
    return ListOrbitFaceOrbitsize.GetListFaces();
  }
}

std::vector<size_t> get_subset_index_rev(const size_t &n_act) {
  // The size of the selection done
  size_t n_ent_bit = 8 * sizeof(size_t);
  size_t n_bit_hash = n_ent_bit;
  if (n_act <= n_ent_bit)
    n_bit_hash = n_act;
  std::vector<size_t> subset_index(n_bit_hash);
  size_t pos_wrt = n_bit_hash;
  ;
  if (n_act <= n_ent_bit) {
    for (size_t i = 0; i < n_ent_bit; i++) {
      pos_wrt--;
      subset_index[pos_wrt] = i;
    }
  } else {
    double frac =
        static_cast<double>(n_act - 1) / static_cast<double>(n_ent_bit - 1);
    for (size_t i = 0; i < n_ent_bit; i++) {
      int pos_i = lround(frac * static_cast<double>(i));
      size_t pos;
      if (pos_i < 0)
        pos = 0;
      else
        pos = pos_i;
      if (pos >= n_act)
        pos = n_act - 1;
      pos_wrt--;
      subset_index[pos_wrt] = pos;
    }
  }
  return subset_index;
}

template <typename Tidx>
size_t evaluate_subset_hash(const std::vector<Tidx> &subset_index,
                            const Face &f) {
  size_t hash = 0;
  size_t *ptr1 = &hash;
  uint8_t *ptr2 = reinterpret_cast<uint8_t *>(ptr1);
  size_t n_bit_hash = subset_index.size();
  for (size_t i = 0; i < n_bit_hash; i++) {
    bool val = f[subset_index[i]];
    setbit_ptr(ptr2, i, val);
  }
  return hash;
}

FullNamelist NAMELIST_GetStandard_RecursiveDualDescription() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, std::string> ListStringValues1_doc;
  std::map<std::string, std::string> ListBoolValues1_doc;
  std::map<std::string, std::string> ListIntValues1_doc;
  ListStringValues1_doc["NumericalType"] = "Default: rational\n\
The numerical type being used for the computation. Possible values:\n\
rational: the rational type, what you want in 99.999\% of cases\n\
Qsqrt5: coordinates in the field Q(sqrt(5))\n\
Qsqrt2: coordinates in the field Q(sqrt(2))\n\
RealAlgebraic: coordinate in a real algebraic field";
  ListStringValues1_doc["FileAlgebraicField"] = "Default: unset\n\
The file containing the description of the real algebraic field.\n\
This is needed of RealAlgebraic is selected";
  ListStringValues1_doc["EXTfile"] =
      "The file containing the coordinate of the output file";
  ListStringValues1_doc["GRPfile"] =
      "The file containing the symmetry group used in the computation";
  ListStringValues1_doc["OUTfile"] =
      "The file containing the output of the result";
  ListStringValues1_doc["OutFormat"] = "Default: GAP\n\
The formatting used for the output. Possible values:\n\
Magma: a file to be read in magma\n\
GAP: a file encoding the incidence as list and made a file openable in GAP\n\
SetInt: a file encoding the incidence as a single integer\n\
BankEntry: a file encoding the dual description as a bank entry that can be used for the bank system";
  ListBoolValues1_doc["DeterministicRuntime"] = "Default: F\n\
There is some randomness in several algorithms. With DeterministicRuntime:\n\
T: If you run again the program you will get exactly the same result which is good for debugging\n\
F: Running again the program will get you something different";
  ListBoolValues1_doc["ApplyStdUnitbuf"] = "Default: F\n\
There is some logging being done in the running of the program. With AppluStdUnit\n\
T: the output is done character by character which is slower but useful for debugging\n\
F: the output is buffered which is typically faster";
  ListStringValues1_doc["bank_parallelization_method"] = "Default: serial\n\
The method used for parallelizing the banking system\n\
serial: Every thread has its own banking system, which may be suboptimal\n\
  since other thread may have the dual description you computed\n\
bank_asio: a parallel bank used by several process\n\
bank_mpi: a bank shared by all the mpi threads";
  ListIntValues1_doc["port"] = "Default: 1234\n\
The port used for the bank_asio";
  ListIntValues1_doc["max_runtime"] = "Default: -1\n\
The maximum runtime of the run in seconds.\n\
If data is saved then you can rerun with the saved state\n\
if max_runtime is negative then there is no maximum runtime";
  ListBoolValues1_doc["AdvancedTerminationCriterion"] = "Default: F\n\
This is about whether to used the advanced Balinski termination criterion";
  ListBoolValues1_doc["SimpleExchangeScheme"] = "Default: F\n\
If selected then a message sent to another node can be sent only after the previously sent is marked as finished";
  SingleBlock BlockDATA;
  BlockDATA.setListStringValues(ListStringValues1_doc);
  BlockDATA.setListBoolValues(ListBoolValues1_doc);
  BlockDATA.setListIntValues(ListIntValues1_doc);
  ListBlock["DATA"] = BlockDATA;
  // HEURISTIC
  std::map<std::string, std::string> ListStringValuesH_doc;
  ListStringValuesH_doc["SplittingHeuristicFile"] = "Default: unset.heu\n\
The splitting heuristic file.\n\
If set to unset.heu then basic heuristics are applied which should be fine for small case";
  ListStringValuesH_doc["AdditionalSymmetryHeuristicFile"] =
      "Default: unset.heu\n\
The additional symmetry heuristic file\n\
If set to unset.heu then basic heuristics are applied which should be fine for small case";
  ListStringValuesH_doc["DualDescriptionThompsonFile"] = "Default: unset.ts\n\
The Thompson samspling heuristic file for choosing the dual description program.\n\
If set to unset.ts then basic heuristics are applied which should be fine for small case";
  ListStringValuesH_doc["MethodInitialFacetSetFile"] = "Default: unset.heu\n\
The heuristic for computing the initial set of facets.\n\
If set to unset.heu then basic heuristics are applied which should be fine for small case";
  ListStringValuesH_doc["BankSaveHeuristicFile"] = "Default: unset.heu\n\
The heuristic file whether to save computed data to the bank or not.\n\
If set to unset.heu then basic heuristics are applied which should be fine for small case";
  ListStringValuesH_doc["CheckDatabaseBankFile"] = "Default: unset.heu\n\
The heuristic file file for checking if entries are present in the bank.\n\
If set to unset.heu then basic heuristics are applied which should be fine for small case";
  ListStringValuesH_doc["ChosenDatabaseFile"] = "Default: unset.heu\n\
The heuristic for choosing between canonic or repr.\n\
If set to unset.heu then basic heuristics are applied which should be fine for small case";
  ListStringValuesH_doc["OrbitSplitTechniqueFile"] = "Default: unset.heu\n\
The heuristic for choosing the orbit splitting technique.\n\
If set to unset.heu then basic heuristics are applied which should be fine";
  SingleBlock BlockHEURIS;
  BlockHEURIS.setListStringValues(ListStringValuesH_doc);
  ListBlock["HEURISTIC"] = BlockHEURIS;
  // METHOD
  std::map<std::string, std::string> ListBoolValues2_doc;
  std::map<std::string, std::string> ListStringValues2_doc;
  ListBoolValues2_doc["Saving"] = "Default: F\n\
Whether to save the bank information to a disk for further reuse";
  ListStringValues2_doc["Prefix"] = "Default: /irrelevant/\n\
The directory in which the bank is saved. Put something significant if Saving = T";
  SingleBlock BlockMETHOD;
  BlockMETHOD.setListBoolValues(ListBoolValues2_doc);
  BlockMETHOD.setListStringValues(ListStringValues2_doc);
  ListBlock["METHOD"] = BlockMETHOD;
  // BANK
  std::map<std::string, std::string> ListBoolValues3_doc;
  std::map<std::string, std::string> ListStringValues3_doc;
  ListBoolValues3_doc["Saving"] = "Default: F\n\
Whether to track the computation on file or not";
  ListStringValues3_doc["Prefix"] = "Default: /irrelevant/\n\
The prefix in which data is saved. Put something significant if Saving = T";
  SingleBlock BlockBANK;
  BlockBANK.setListBoolValues(ListBoolValues3_doc);
  BlockBANK.setListStringValues(ListStringValues3_doc);
  ListBlock["BANK"] = BlockBANK;
  // Merging all data
  return {std::move(ListBlock), "undefined"};
}

FullNamelist NAMELIST_GetStandard_BankingSystem() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListBoolValues1["Saving"] = false;
  ListStringValues1["Prefix"] = "/irrelevant/";
  ListIntValues1["port"] = 1234;
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues = ListIntValues1;
  BlockPROC.ListBoolValues = ListBoolValues1;
  BlockPROC.ListStringValues = ListStringValues1;
  ListBlock["PROC"] = BlockPROC;
  // Merging all data
  return {std::move(ListBlock), "undefined"};
}

template <typename T, typename Tgroup>
void OutputFacets(const MyMatrix<T> &EXT, Tgroup const &GRP,
                  const vectface &TheOutput, const std::string &OUTfile,
                  const std::string &OutFormat) {
  if (OutFormat == "Magma") {
    return VectVectInt_Magma_PrintFile(OUTfile, TheOutput);
  }
  if (OutFormat == "GAP") {
    return VectVectInt_Gap_PrintFile(OUTfile, TheOutput);
  }
  if (OutFormat == "SetInt") {
    return VectVectInt_SetInt_PrintFile<mpz_class>(OUTfile, TheOutput);
  }
  if (OutFormat == "BankEntry") {
    // We are creating a bank entry for further works.
    using Tidx_value = uint16_t;
    using Tint=typename Tgroup::Tint;
    WeightMatrix<true, T, Tidx_value> WMat =
        GetWeightMatrix<T, Tidx_value>(EXT);
    WMat.ReorderingSetWeight();
    FaceOrbitsizeTableContainer<Tint> fotc(TheOutput, GRP);
    std::pair<MyMatrix<T>, TripleStore<Tgroup>> eP =
        GetCanonicalInformation(EXT, WMat, GRP, fotc);
    Write_BankEntry(OUTfile, eP.first, eP.second);
  }
  std::cerr << "No option has been chosen\n";
  throw TerminalException{1};
}

template <typename T> MyMatrix<T> GetEXT_from_efull(FullNamelist const &eFull) {
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  std::string EXTfile = BlockDATA.ListStringValues.at("EXTfile");
  IsExistingFileDie(EXTfile);
  std::ifstream EXTfs(EXTfile);
  return ReadMatrix<T>(EXTfs);
}

std::string GetNumericalType(FullNamelist const &eFull) {
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  std::string NumericalType = BlockDATA.ListStringValues.at("NumericalType");
  std::vector<std::string> Ltype{"rational", "Qsqrt5", "RealAlgebraic"};
  if (PositionVect(Ltype, NumericalType) == -1) {
    std::cerr << "NumericalType=" << NumericalType << "\n";
    std::cerr << "Ltype =";
    for (auto &e_type : Ltype)
      std::cerr << " " << e_type;
    std::cerr << "\n";
    throw TerminalException{1};
  }
  return NumericalType;
}

template <typename T, typename Tidx>
MyMatrix<T> Get_EXT_DualDesc(FullNamelist const &eFull, std::ostream &os) {
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  std::string EXTfile = BlockDATA.ListStringValues.at("EXTfile");
  IsExistingFileDie(EXTfile);
  os << "EXTfile=" << EXTfile << "\n";
  std::ifstream EXTfs(EXTfile);
  MyMatrix<T> EXT = ReadMatrix<T>(EXTfs);
  if (size_t(EXT.rows()) > size_t(std::numeric_limits<Tidx>::max())) {
    std::cerr << "We have |EXT|=" << EXT.rows() << "\n";
    std::cerr << "But <Tidx>::max()="
              << size_t(std::numeric_limits<Tidx>::max()) << "\n";
    throw TerminalException{1};
  }
  return EXT;
}

bool ApplyStdUnitbuf(FullNamelist const &eFull) {
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  bool result = BlockDATA.ListBoolValues.at("ApplyStdUnitbuf");
  return result;
}

template <typename Tgroup>
Tgroup Get_GRP_DualDesc(FullNamelist const &eFull, std::ostream &os) {
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  std::string GRPfile = BlockDATA.ListStringValues.at("GRPfile");
  IsExistingFileDie(GRPfile);
  os << "GRPfile=" << GRPfile << "\n";
  std::ifstream GRPfs(GRPfile);
  Tgroup GRP = ReadGroup<Tgroup>(GRPfs);
  return GRP;
}

template <typename T, typename Tint>
PolyHeuristicSerial<Tint>
Read_AllStandardHeuristicSerial(FullNamelist const &eFull,
                                MyMatrix<T> const &EXTred, std::ostream &os) {
  PolyHeuristicSerial<Tint> AllArr = AllStandardHeuristicSerial<Tint>(os);
  os << "We have AllArr\n";
  //
  SingleBlock BlockMETHOD = eFull.ListBlock.at("METHOD");
  SingleBlock BlockBANK = eFull.ListBlock.at("BANK");
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  //
  bool BANK_IsSaving = BlockBANK.ListBoolValues.at("Saving");
  AllArr.BANK_IsSaving = BANK_IsSaving;
  //
  std::string BANK_Prefix = BlockBANK.ListStringValues.at("Prefix");
  AllArr.BANK_Prefix = BANK_Prefix;
  //
  std::string OUTfile = BlockDATA.ListStringValues.at("OUTfile");
  AllArr.OUTfile = OUTfile;
  os << "OUTfile=" << OUTfile << "\n";
  //
  bool DeterministicRuntime =
      BlockDATA.ListBoolValues.at("DeterministicRuntime");
  os << "DeterministicRuntime=" << DeterministicRuntime << "\n";
  if (!DeterministicRuntime)
    srand_random_set();
  //
  std::string OutFormat = BlockDATA.ListStringValues.at("OutFormat");
  AllArr.OutFormat = OutFormat;
  os << "OutFormat=" << OutFormat << "\n";
  //
  int port_i = BlockDATA.ListIntValues.at("port");
  os << "port_i=" << port_i << "\n";
  short unsigned int port = port_i;
  AllArr.port = port;
  //
  std::string bank_parallelization_method =
      BlockDATA.ListStringValues.at("bank_parallelization_method");
  AllArr.bank_parallelization_method = bank_parallelization_method;
  os << "bank_parallelization_method=" << bank_parallelization_method << "\n";
  //
  SetHeuristic(eFull, "SplittingHeuristicFile", AllArr.Splitting, os);
  os << "SplittingHeuristicFile\n" << AllArr.Splitting << "\n";
  //
  SetHeuristic(eFull, "AdditionalSymmetryHeuristicFile",
               AllArr.AdditionalSymmetry, os);
  os << "AdditionalSymmetryHeuristicFile\n"
     << AllArr.AdditionalSymmetry << "\n";
  //
  SetThompsonSampling(eFull, "DualDescriptionThompsonFile",
                      AllArr.DualDescriptionProgram, os);
  os << "DualDescriptionThompsonFile\n"
     << AllArr.DualDescriptionProgram << "\n";
  //
  SetHeuristic(eFull, "MethodInitialFacetSetFile", AllArr.InitialFacetSet, os);
  os << "MethodInitialFacetSetFile\n" << AllArr.InitialFacetSet << "\n";
  //
  SetHeuristic(eFull, "BankSaveHeuristicFile", AllArr.BankSave, os);
  os << "BankSaveHeuristicFile\n" << AllArr.BankSave << "\n";
  //
  SetHeuristic(eFull, "CheckDatabaseBankFile", AllArr.CheckDatabaseBank, os);
  os << "CheckDatabaseBank\n" << AllArr.CheckDatabaseBank << "\n";
  //
  SetHeuristic(eFull, "ChosenDatabaseFile", AllArr.ChosenDatabase, os);
  os << "ChosenDatabase\n" << AllArr.ChosenDatabase << "\n";
  //
  SetHeuristic(eFull, "OrbitSplitTechniqueFile", AllArr.OrbitSplitTechnique, os);
  os << "OrbitSplitTechnique\n" << AllArr.OrbitSplitTechnique << "\n";
  //
  bool DD_Saving = BlockMETHOD.ListBoolValues.at("Saving");
  AllArr.Saving = DD_Saving;
  //
  std::string DD_Prefix = BlockMETHOD.ListStringValues.at("Prefix");
  AllArr.DD_Prefix = DD_Prefix;
  //
  int max_runtime = BlockDATA.ListIntValues.at("max_runtime");
  AllArr.max_runtime = max_runtime;
  //
  bool AdvancedTerminationCriterion =
      BlockDATA.ListBoolValues.at("AdvancedTerminationCriterion");
  AllArr.AdvancedTerminationCriterion = AdvancedTerminationCriterion;
  //
  bool SimpleExchangeScheme = BlockDATA.ListBoolValues.at("SimpleExchangeScheme");
  AllArr.SimpleExchangeScheme = SimpleExchangeScheme;
  //
  AllArr.dimEXT = EXTred.cols();
  //
  return AllArr;
}

template <typename T, typename Tgroup, typename Tidx_value>
void MainFunctionSerialDualDesc(FullNamelist const &eFull) {
  // Setting up the Control C event.
  ExitEvent = false;
  std::cerr << "Before submission of signal_callback_handler\n";
  signal(SIGINT, signal_callback_handler);
  //
  using Tint = typename Tgroup::Tint;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tkey = MyMatrix<T>;
  using Tval = TripleStore<Tgroup>;
  MyMatrix<T> EXT = Get_EXT_DualDesc<T, Tidx>(eFull, std::cerr);
  Tgroup GRP = Get_GRP_DualDesc<Tgroup>(eFull, std::cerr);
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  PolyHeuristicSerial<Tint> AllArr =
      Read_AllStandardHeuristicSerial<T, Tint>(eFull, EXTred, std::cerr);
  //
  std::map<std::string, Tint> TheMap =
    ComputeInitialMap<Tint>(EXTred, GRP, AllArr);
  auto get_vectface = [&]() -> vectface {
    if (AllArr.bank_parallelization_method == "serial") {
      using Tbank = DataBank<Tkey, Tval>;
      Tbank TheBank(AllArr.BANK_IsSaving, AllArr.BANK_Prefix, std::cerr);
      return DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value>(
        TheBank, EXTred, GRP, TheMap, AllArr, AllArr.DD_Prefix, std::cerr);
    }
    if (AllArr.bank_parallelization_method == "bank_asio") {
      using Tbank = DataBankAsioClient<Tkey, Tval>;
      Tbank TheBank(AllArr.port);
      return DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value>(
        TheBank, EXTred, GRP, TheMap, AllArr, AllArr.DD_Prefix, std::cerr);
    }
    std::cerr
        << "Failed to find a matching entry for bank_parallelization_method\n";
    std::cerr << "Allowed methods are serial, bank_asio\n";
    throw TerminalException{1};
  };
  vectface TheOutput = get_vectface();
  std::cerr << "|TheOutput|=" << TheOutput.size() << "\n";
  //
  OutputFacets(EXT, GRP, TheOutput, AllArr.OUTfile, AllArr.OutFormat);
}

template <typename T, typename Tgroup>
vectface DualDescriptionStandard(const MyMatrix<T> &EXT, const Tgroup &GRP) {
  using Tint = typename Tgroup::Tint;
  using Tkey = MyMatrix<T>;
  using Tval = TripleStore<Tgroup>;
  using Tidx_value = int32_t;
  bool BANK_IsSaving = false;
  std::string BANK_Prefix = "totally_irrelevant_first";
  //
  PolyHeuristicSerial<Tint> AllArr =
      AllStandardHeuristicSerial<Tint>(std::cerr);
  std::cerr << "SplittingHeuristicFile\n" << AllArr.Splitting << "\n";
  std::cerr << "AdditionalSymmetryHeuristicFile\n"
            << AllArr.AdditionalSymmetry << "\n";
  std::cerr << "DualDescriptionHeuristicFile\n"
            << AllArr.DualDescriptionProgram << "\n";
  std::cerr << "MethodInitialFacetSetFile\n" << AllArr.InitialFacetSet << "\n";
  std::cerr << "BankSaveHeuristicFile\n" << AllArr.BankSave << "\n";
  std::cerr << "CheckDatabaseBank\n" << AllArr.CheckDatabaseBank << "\n";
  std::cerr << "ChosenDatabase\n" << AllArr.ChosenDatabase << "\n";
  std::cerr << "OrbitSplitTechique\n" << AllArr.OrbitSplitTechnique << "\n";
  //
  bool DD_Saving = false;
  std::string DD_Prefix = "totally_irrelevant_second";
  AllArr.Saving = DD_Saving;
  //
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  using Tbank = DataBank<Tkey, Tval>;
  Tbank TheBank(BANK_IsSaving, BANK_Prefix, std::cerr);
  std::map<std::string, Tint> TheMap =
    ComputeInitialMap<Tint>(EXTred, GRP, AllArr);
  return DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value>(
    TheBank, EXTred, GRP, TheMap, AllArr, DD_Prefix, std::cerr);
}

// clang-format off
#endif  // SRC_DUALDESC_POLY_RECURSIVEDUALDESC_H_
// clang-format on

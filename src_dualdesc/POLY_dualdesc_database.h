// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_DUALDESC_DATABASE_H_
#define SRC_DUALDESC_POLY_DUALDESC_DATABASE_H_

// clang-format off
#include "GRP_DoubleCoset.h"
#include "POLY_Fundamental.h"
// clang-format on


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

// Those constants are for the canonic strategy since we have 3 different
// methods available
static const int CANONIC_STRATEGY__CANONICAL_IMAGE = 0;
static const int CANONIC_STRATEGY__STORE = 1;
static const int CANONIC_STRATEGY__INITIAL_TRIV = 2;
// This is for the constant 5000
static const int CANONIC_STRATEGY__INITIAL_TRIV_LIMITED1 = 4;

// Those constants are for the default strategy
static const int CANONIC_STRATEGY__DEFAULT = CANONIC_STRATEGY__CANONICAL_IMAGE;
static const int REPR_STRATEGY__DEFAULT = 0; // More or less irrelevant here

// Constant indicating failure of finding the miss
static const int STRATEGY_MISS = -1;

// Those constants express the choice for the database
static const int DATABASE_ACTION__SIMPLE_LOAD = 0;
static const int DATABASE_ACTION__GUESS = 1;
static const int DATABASE_ACTION__RECOMPUTE_AND_SHUFFLE = 2;

template <typename Tgroup_impl> struct TripleStore {
  using Tgroup = Tgroup_impl;
  using Tint = typename Tgroup::Tint;
  Tgroup GRP;
  std::vector<Tint> ListPossOrbsize;
  vectface ListFace;
};

namespace boost::serialization {

template <class Archive, typename Tgroup>
inline void serialize(Archive &ar, TripleStore<Tgroup> &triple,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("GRP", triple.GRP);
  ar &make_nvp("ListPossOrbsize", triple.ListPossOrbsize);
  ar &make_nvp("ListFace", triple.ListFace);
}

// clang-format off
}  // namespace boost::serialization
// clang-format on

size_t get_matching_power(size_t const &val) {
  size_t pow = 1;
  size_t pos = 0;
  while (true) {
    if (pow >= val)
      return pos;
    pow *= 2;
    pos++;
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
                             WeightMatrix<true, T, Tidx_value> const &WMat,
                             std::ostream &os) {
  std::vector<Tidx> CanonicOrd =
      GetCanonicalizationVector_Kernel<T, GraphBitset, Tidx>(WMat, os);
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
                               WeightMatrix<true, T, Tidx_value> const &WMat,
                               std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> PairCanGrp =
      GetGroupCanonicalizationVector_Kernel<T, GraphBitset, Tidx>(WMat, os);
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
MyMatrix<T> CanonicalizationPolytope(MyMatrix<T> const &EXT, std::ostream &os) {
  using Tidx_value = uint16_t;
  WeightMatrix<true, T, Tidx_value> WMat =
      GetWeightMatrix<T, Tidx_value>(EXT, os);
  WMat.ReorderingSetWeight();
  return CanonicalizationPolytopePair<T, int, Tidx_value>(EXT, WMat, os).first;
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
      for (auto &eVal : LVal) {
        NewVal.push_back(ePow * eVal);
      }
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
Face CanonicalImageDualDesc(int const &method_choice, Tgroup const &GRP,
                            Face const &f, [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_CANONICAL_LIMITED_V2
  os << "RDD: Entry " << StringGroup(GRP) << " " << StringFace(f) << "\n";
#endif
#ifdef DEBUG_CANONICAL_LIMITED
  os << "RDD: CAN_LIM: Beginning of CanonicalImageDualDesc method_choice="
     << method_choice << "\n";
  WriteGroup(os, GRP);
  os << "RDD: f=" << StringFace(f) << "\n";
#endif
  if (method_choice == CANONIC_STRATEGY__CANONICAL_IMAGE) {
    Face f_red = GRP.CanonicalImage(f);
#ifdef DEBUG_CANONICAL_LIMITED
    os << "RDD: CAN_LIM: After CanonicalImage\n";
#endif
    return f_red;
  }
  if (method_choice == CANONIC_STRATEGY__STORE) {
    Face f_red = GRP.StoreCanonicalImage(f);
#ifdef DEBUG_CANONICAL_LIMITED
    os << "RDD: CAN_LIM: After StoreCanonicalImage\n";
#endif
    return f_red;
  }
  if (method_choice == CANONIC_STRATEGY__INITIAL_TRIV) {
    Face f_red = GRP.CanonicalImageInitialTriv(f);
#ifdef DEBUG_CANONICAL_LIMITED
    os << "RDD: CAN_LIM: After CanonicalImageInitialTriv\n";
#endif
    return f_red;
  }
  if (method_choice == CANONIC_STRATEGY__INITIAL_TRIV_LIMITED1) {
    try {
      Face f_red = GRP.CanonicalImageInitialTrivLimited(f, LIMIT_INITIAL_TRIV);
#ifdef DEBUG_CANONICAL_LIMITED
      os << "RDD: CAN_LIM: After CanonicalImageInitialTrivLimited\n";
#endif
      return f_red;
    } catch (...) {
      std::cerr
          << "Catching some exception from CanonicalImageInitialTrivLimited\n";
      throw TerminalException{1};
    }
  }
  std::cerr << "Error in CanonicalImageDualDesc, no method found\n";
  std::cerr << "method_choice=" << method_choice << "\n";
  throw TerminalException{1};
}

/*
  The pair Face, OrbSize can be encoded as a single Face.
 */
template <typename Torbsize, typename Tgroup> struct DataFaceOrbitSize {
  using TintGroup = typename Tgroup::Tint;
  /* TRICK 3: Knowing the factorization of the order of the group allow us to
     know exactly
     what are the possible orbitsize occurring and so the number of bits needed
     to encode them */
  std::vector<TintGroup> ListPossOrbsize;
  /* TRICK 2: We keep the list of orbit and the map. We could in principle have
     built the map from the start since we know the occurring orders. However,
     since some orbitsize never occur
     this would have populated it with entries that never occur and so slow it
     down. */
  UNORD_MAP<TintGroup, Torbsize> OrbSize_Map;
  size_t n;
  size_t n_bit_orbsize;
  size_t delta;
  DataFaceOrbitSize(Tgroup const &GRP) {
    using Tidx = typename Tgroup::Telt::Tidx;
    std::map<Tidx, int> LFact = GRP.factor_size();
    n = GRP.n_act();
    std::pair<size_t, size_t> ep = get_delta(LFact, n);
    ListPossOrbsize = GetAllPossibilities<Tidx, TintGroup>(LFact);
    n_bit_orbsize = ep.first;
    delta = ep.second;
  }
  Torbsize GetOrbSizeIndex(TintGroup const &orbSize) {
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
  Face ConvertFaceOrbitSize(std::pair<Face, TintGroup> const &pair) {
    Face const &f = pair.first;
    TintGroup const &orbitSize = pair.second;
    Torbsize idx_orb = GetOrbSizeIndex(orbitSize);
    //
    Face f_ret(delta);
    for (size_t i = 0; i < n; i++)
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
  std::pair<Face, TintGroup> ConvertFace(Face const &f) const {
    Face f_ret(n);
    for (size_t i = 0; i < n; i++)
      f_ret[i] = f[i];
    size_t idx_orb = 0;
    size_t pow = 1;
    size_t i_acc = n;
    for (size_t i = 0; i < n_bit_orbsize; i++) {
      if (f[i_acc] == 1) {
        idx_orb += pow;
      }
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
Face CanonicalImageGeneralDualDesc(
    int const &method_choice, Tgroup const &GRP,
    DataFaceOrbitSize<Torbsize, Tgroup> &recConvert, Face const &f,
    [[maybe_unused]] std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
#ifdef DEBUG_CANONICAL_LIMITED_V2
  os << "RDD: Entry " << StringGroup(GRP) << " " << StringFace(f) << "\n";
#endif
#ifdef DEBUG_CANONICAL_LIMITED
  os << "RDD: CAN_LIM: Beginning of CanonicalImageGeneralDualDesc "
        "method_choice="
     << method_choice << "\n";
  WriteGroup(os, GRP);
  os << "RDD: CAN_LIM: f=" << StringFace(f) << "\n";
#endif
  if (method_choice == CANONIC_STRATEGY__CANONICAL_IMAGE) {
    std::pair<Face, TintGroup> pair = GRP.CanonicalImageOrbitSize(f);
#ifdef DEBUG_CANONICAL_LIMITED
    os << "RDD: CAN_LIM: After CanonicalImageOrbitSize\n";
#endif
    return recConvert.ConvertFaceOrbitSize(pair);
  }
  if (method_choice == CANONIC_STRATEGY__STORE) {
    std::pair<Face, TintGroup> pair = GRP.StoreCanonicalImageOrbitSize(f);
#ifdef DEBUG_CANONICAL_LIMITED
    os << "RDD: CAN_LIM: After StoreCanonicalImageOrbitSize\n";
#endif
    return recConvert.ConvertFaceOrbitSize(pair);
  }
  if (method_choice == CANONIC_STRATEGY__INITIAL_TRIV) {
    Face f_red = GRP.CanonicalImageInitialTriv(f);
#ifdef DEBUG_CANONICAL_LIMITED
    os << "RDD: CAN_LIM: After CanonicalImageInitialTriv\n";
#endif
    return f_red;
  }
  if (method_choice == CANONIC_STRATEGY__INITIAL_TRIV_LIMITED1) {
    try {
      Face f_red = GRP.CanonicalImageInitialTrivLimited(f, LIMIT_INITIAL_TRIV);
#ifdef DEBUG_CANONICAL_LIMITED
      os << "RDD: CAN_LIM: After CanonicalImageInitialTrivLimited\n";
#endif
      return f_red;
    } catch (...) {
      std::cerr
          << "Catching some exception from CanonicalImageInitialTrivLimited\n";
      throw TerminalException{1};
    }
  }
  std::cerr << "Error in CanonicalImageOrbitSizeDualDesc, no method found\n";
  std::cerr << "method_choice=" << method_choice << "\n";
  throw TerminalException{1};
}

vectface vectface_reduction(vectface const &vf, size_t n_red) {
  vectface vf_red(n_red);
  Face f_red(n_red);
  for (auto &f : vf) {
    for (size_t i = 0; i < n_red; i++) {
      f_red[i] = f[i];
    }
    vf_red.push_back(f_red);
  }
  return vf_red;
}

Face face_reduction(Face const &f, size_t n_red) {
  Face f_red(n_red);
  for (size_t i = 0; i < n_red; i++) {
    f_red[i] = f[i];
  }
  return f_red;
}

void set_face_partial(Face &f_out, Face const &f_in, size_t const &n_red) {
  for (size_t i = 0; i < n_red; i++) {
    f_out[i] = f_in[i];
  }
}

template <typename TintGroup> struct FaceOrbitsizeTableContainer {
public:
  std::vector<TintGroup> ListPossOrbsize;
  size_t n;
  vectface vfo;
  FaceOrbitsizeTableContainer(std::vector<TintGroup> const &_ListPossOrbsize,
                              size_t _n, vectface &&_vfo)
      : ListPossOrbsize(std::move(_ListPossOrbsize)), n(_n),
        vfo(std::move(_vfo)) {}
  template <typename Tgroup>
  FaceOrbitsizeTableContainer(vectface const &vf, Tgroup const &GRP) {
    n = vf.get_n();
    using Tidx = typename Tgroup::Telt::Tidx;
    std::map<Tidx, int> LFact = GRP.factor_size();
    std::pair<size_t, size_t> ep = get_delta(LFact, n);
    size_t n_bit_orbsize = ep.first;
    size_t delta = ep.second;
    ListPossOrbsize = GetAllPossibilities<Tidx, TintGroup>(LFact);
    UNORD_MAP<TintGroup, size_t> OrbSize_Map;
    for (size_t i = 0; i < ListPossOrbsize.size(); i++) {
      OrbSize_Map[ListPossOrbsize[i]] = i;
    }
    vectface vf_ins(delta);
    vfo = std::move(vf_ins);
    for (auto &eFace : vf) {
      TintGroup orbitSize = GRP.OrbitSize_OnSets(eFace);
      size_t idx_orb = OrbSize_Map[orbitSize];
      Face f(delta);
      for (size_t i = 0; i < n; i++)
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
  std::pair<Face, TintGroup> GetPair(size_t const &idx) const {
    Face f = vfo[idx];
    Face f_red(n);
    for (size_t i = 0; i < n; i++)
      f_red[i] = f[i];
    size_t pow = 1;
    size_t idx_orb = 0;
    for (size_t i = n; i < vfo.get_n(); i++) {
      if (f[i] == 1)
        idx_orb += pow;
      pow *= 2;
    }
    TintGroup orbSize = ListPossOrbsize[idx_orb];
    return {f_red, orbSize};
  }
  size_t size() const { return vfo.size(); }
  vectface GetListFaces() const { return vectface_reduction(vfo, n); }
};

template <typename Telt>
Telt trivial_extension(Telt const &ePerm, size_t const &delta) {
  using Tidx = typename Telt::Tidx;
  Tidx delta_tidx = delta;
  std::vector<Tidx> V(delta);
  Tidx siz = ePerm.size();
  for (Tidx i = 0; i < siz; i++) {
    Tidx ePt = ePerm.at(i);
    V[i] = ePt;
  }
  for (Tidx i = siz; i < delta_tidx; i++) {
    V[i] = i;
  }
  return Telt(V);
}

template <typename Tgroup>
Tgroup trivial_extension_group(Tgroup const &eGroup, size_t const &delta) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<Telt> LGen = eGroup.GeneratorsOfGroup();
  std::vector<Telt> LGenExt;
  for (auto &eGen : LGen) {
    Telt eGenExt = trivial_extension(eGen, delta);
    LGenExt.emplace_back(std::move(eGenExt));
  }
  std::vector<Tidx> V(delta);
  for (size_t i = 0; i < delta; i++) {
    V[i] = i;
  }
  Telt id(V);
  return Tgroup(LGenExt, id);
}

template <typename Tgroup>
std::vector<int> GetPossibleCanonicalizationMethod(Tgroup const &GRP) {
  // We put first the CANONIC_STRATEGY__CANONICAL_IMAGE as it is an all around
  // reasonable method on which other methods have to compete with.
  std::vector<int> list_considered = {CANONIC_STRATEGY__CANONICAL_IMAGE,
                                      CANONIC_STRATEGY__INITIAL_TRIV_LIMITED1};
  if (GRP.size() < 20000) { // We need to exclude that strategy if too large as
                            // that strategy has no chance.
    list_considered.push_back(CANONIC_STRATEGY__STORE);
  }
  return list_considered;
}

template <typename Tgroup>
int64_t time_evaluation_can_method(int const &method, vectface const &vf,
                                   Tgroup const &GRP, int64_t upper_limit,
                                   std::ostream &os) {
  NanosecondTime time;
  int64_t duration = 0;
  int64_t miss_val = std::numeric_limits<int64_t>::max();
  for (auto &f : vf) {
    (void)CanonicalImageDualDesc(method, GRP, f, os);
    duration = time.const_eval_int64();
    if (duration > upper_limit)
      return miss_val;
  }
  return duration;
}

template <typename Tgroup>
int GetCanonicalizationMethod_Serial(vectface const &vf, Tgroup const &GRP,
                                     std::ostream &os) {
  std::vector<int> list_considered = GetPossibleCanonicalizationMethod(GRP);
  int64_t upper_limit = std::numeric_limits<int64_t>::max();
  int chosen_method = CANONIC_STRATEGY__DEFAULT;
  for (auto &method : list_considered) {
    int64_t runtime =
        time_evaluation_can_method(method, vf, GRP, upper_limit, os);
    if (runtime < upper_limit) {
      chosen_method = method;
      upper_limit = runtime;
    }
  }
  return chosen_method;
}

template <typename T, typename Tgroup>
int GetCanonicalizationMethodRandom(MyMatrix<T> const &EXT, Tgroup const &GRP,
                                    size_t size, std::ostream &os) {
  if (size < 10000) {
    if (GRP.size() < 200)
      return CANONIC_STRATEGY__STORE;
    return CANONIC_STRATEGY__INITIAL_TRIV_LIMITED1;
  }
  int n = EXT.rows();
  vectface vf(n);
  for (int iter = 0; iter < 100; iter++) {
    Face f = RandomFace(n);
    vf.push_back(f);
  }
  return GetCanonicalizationMethod_Serial(vf, GRP, os);
}

template <typename T, typename Tgroup, typename Tidx_value>
std::pair<MyMatrix<T>, TripleStore<Tgroup>>
GetCanonicalInformation(MyMatrix<T> const &EXT,
                        WeightMatrix<true, T, Tidx_value> const &WMat,
                        Tgroup const &TheGRPrelevant,
                        FaceOrbitsizeTableContainer<typename Tgroup::Tint> const
                            &ListOrbitFaceOrbitsize,
                        std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  std::vector<TintGroup> ListPossOrbSize =
      ListOrbitFaceOrbitsize.ListPossOrbsize;
  TripleCanonic<T, Tgroup> eTriple =
      CanonicalizationPolytopeTriple<T, Tgroup>(EXT, WMat, os);
  bool NeedRemapOrbit = eTriple.GRP.size() == TheGRPrelevant.size();
  size_t delta = ListOrbitFaceOrbitsize.vfo.get_n();
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
    int can_method =
        GetCanonicalizationMethodRandom(EXT, TheGRPrelevant, size, os);
    UNORD_SET<Face> SetFace;
    Tgroup GRPext = trivial_extension_group(eTriple.GRP, delta);
    for (auto &eFace : ListOrbitFaceOrbitsize.vfo) {
      OnFace_inplace(eFaceImg, eFace, ePermExt);
      Face eIncCan = CanonicalImageDualDesc(can_method, GRPext, eFaceImg, os);
      SetFace.insert(eIncCan);
    }
    for (auto &eInc : SetFace) {
      ListFaceO.push_back(eInc);
    }
  }
  TripleStore<Tgroup> ePair{eTriple.GRP, std::move(ListPossOrbSize),
                            std::move(ListFaceO)};
  return {std::move(eTriple.EXT), std::move(ePair)};
}

template <typename Tbank, typename T, typename Tgroup, typename Tidx_value>
void insert_entry_in_bank(
    Tbank &bank, MyMatrix<T> const &EXT,
    WeightMatrix<true, T, Tidx_value> const &WMat, Tgroup const &TheGRPrelevant,
    bool const &BankSymmCheck,
    FaceOrbitsizeTableContainer<typename Tgroup::Tint> const
        &ListOrbitFaceOrbitsize,
    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  using Tidx = typename Telt::Tidx;
  size_t delta = ListOrbitFaceOrbitsize.vfo.get_n();
  if (!BankSymmCheck) {
    // The computation was already done for the full symmetry group. Only
    // canonic form is needed.
    std::pair<MyMatrix<T>, std::vector<Tidx>> ePair =
        CanonicalizationPolytopePair<T, Tidx, Tidx_value>(EXT, WMat, os);
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
    std::vector<TintGroup> ListPossOrbSize =
        ListOrbitFaceOrbitsize.ListPossOrbsize;
    bank.InsertEntry(
        std::move(ePair.first),
        {std::move(GrpConj), std::move(ListPossOrbSize), std::move(ListFaceO)});
  } else {
    std::pair<MyMatrix<T>, TripleStore<Tgroup>> eP = GetCanonicalInformation(
        EXT, WMat, TheGRPrelevant, ListOrbitFaceOrbitsize, os);
    bank.InsertEntry(std::move(eP.first), std::move(eP.second));
  }
}

template <typename T, typename Tgroup> struct DataFacet {
  using Tint = typename Tgroup::Tint;
  size_t SelectedOrbit;
  Face eInc;
  FlippingFramework<T> FF;
  const Tgroup &GRP;
  Tgroup Stab;
  Face FlipFace(const Face &f) { return FF.FlipFace(f); }
  Face FlipFaceIneq(std::pair<Face, MyVector<T>> const &pair) {
    return FF.FlipFaceIneq(pair);
  }
};

template <typename Tgroup, typename Torbsize, typename Tidx>
struct FaceOrbsizeContainer {
public:
  using Tint = typename Tgroup::Tint;
  // We CANNOT replace ListOrbit by vectface as we use a number of hacks that
  // would not be available with a vectface.
  std::vector<uint8_t> ListOrbit;
  Tint TotalNumber;
  size_t nbOrbitDone;
  Tint nbUndone;
  size_t nbOrbit;
  size_t n_act;
  size_t n_bit_orbsize;
  size_t delta;
  std::vector<uint8_t> Vappend;
  DataFaceOrbitSize<Torbsize, Tgroup> recConvert;
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
  // Database code that uses ListOrbit;
  std::pair<Face, Tint> RetrieveListOrbitEntry(size_t const &i_orb) const {
    Face f(n_act);
    size_t i_acc = delta * i_orb;
    for (size_t i = 0; i < n_act; i++) {
      f[i] = getbit_vector(ListOrbit, i_acc);
      i_acc++;
    }
    Torbsize idx_orb = 0;
    Torbsize pow = 1;
    for (size_t i = 0; i < n_bit_orbsize; i++) {
      if (getbit_vector(ListOrbit, i_acc))
        idx_orb += pow;
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
  void InsertListOrbitFace_size(Face const &face, size_t const &rel_siz) {
    size_t curr_len = ListOrbit.size();
    size_t needed_bits = (nbOrbit + 1) * delta;
    size_t needed_len = (needed_bits + 7) / 8;
    size_t incr = needed_len - curr_len;
    if (incr > 0)
      ListOrbit.insert(ListOrbit.end(), Vappend.begin(),
                       Vappend.begin() + incr);
    size_t i_acc = nbOrbit * delta;
    for (size_t i = 0; i < rel_siz; i++) {
      bool val = face[i];
      setbit_vector(ListOrbit, i_acc, val);
      i_acc++;
    }
  }
  void InsertListOrbitFace(Face const &face) {
    // Now setting up the bits but only for the faces as this suffices for the
    // comparison of novelty.
    InsertListOrbitFace_size(face, n_act);
  }
  void InsertListOrbitFaceComplete(Face const &face) {
    InsertListOrbitFace_size(face, delta);
  }
  void FinishWithOrbSizeAssignation(Tint const &orbSize) {
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
    return FaceOrbitsizeTableContainer(std::move(ListPoss), n_act,
                                       std::move(vfo));
  }
};







// clang-format off
#endif  // SRC_DUALDESC_POLY_DUALDESC_DATABASE_H_
// clang-format on

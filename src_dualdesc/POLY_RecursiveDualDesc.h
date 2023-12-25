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
#include "MPI_basic.h"
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
#include <algorithm>
#include <vector>
// clang-format on

#ifdef TIMINGS
#define TIMINGS_RECURSIVE_DUAL_DESC
#endif

#ifdef DEBUG
#define DEBUG_CANONICAL_LIMITED
#define DEBUG_RECURSIVE_DUAL_DESC
#endif

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

// Those constants express the choice for the database
static const int DATABASE_ACTION__SIMPLE_LOAD = 0;
static const int DATABASE_ACTION__GUESS = 1;
static const int DATABASE_ACTION__RECOMPUTE_AND_SHUFFLE = 2;

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
Face CanonicalImageDualDesc(int const &method_choice, Tgroup const &GRP,
                            Face const &f, [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_CANONICAL_LIMITED_V2
  os << "Entry " << StringGroup(GRP) << " " << StringFace(f) << "\n";
#endif
#ifdef DEBUG_CANONICAL_LIMITED
  os << "CAN_LIM: Beginning of CanonicalImageDualDesc method_choice=" << method_choice
     << "\n";
  WriteGroup(os, GRP);
  os << "f=" << StringFace(f) << "\n";
#endif
  if (method_choice == CANONIC_STRATEGY__CANONICAL_IMAGE) {
    Face f_red = GRP.CanonicalImage(f);
#ifdef DEBUG_CANONICAL_LIMITED
    os << "CAN_LIM: After CanonicalImage\n";
#endif
    return f_red;
  }
  if (method_choice == CANONIC_STRATEGY__STORE) {
    Face f_red = GRP.StoreCanonicalImage(f);
#ifdef DEBUG_CANONICAL_LIMITED
    os << "CAN_LUM: After StoreCanonicalImage\n";
#endif
    return f_red;
  }
  if (method_choice == CANONIC_STRATEGY__INITIAL_TRIV) {
    Face f_red = GRP.CanonicalImageInitialTriv(f);
#ifdef DEBUG_CANONICAL_LIMITED
    os << "CAN_LIM: After CanonicalImageInitialTriv\n";
#endif
    return f_red;
  }
  if (method_choice == CANONIC_STRATEGY__INITIAL_TRIV_LIMITED1) {
    try {
      Face f_red = GRP.CanonicalImageInitialTrivLimited(f, LIMIT_INITIAL_TRIV);
#ifdef DEBUG_CANONICAL_LIMITED
      os << "CAN_LIM: After CanonicalImageInitialTrivLimited\n";
#endif
      return f_red;
    } catch (...) {
      os << "Catching some exception from CanonicalImageInitialTrivLimited\n";
      throw TerminalException{1};
    }
  }
  std::cerr << "Error in CanonicalImageDualDesc, no method found\n";
  std::cerr << "method_choice=" << method_choice << "\n";
  throw TerminalException{1};
}

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
Face CanonicalImageGeneralDualDesc(
    int const &method_choice, Tgroup const &GRP,
    DataFaceOrbitSize<Torbsize, Tgroup> &recConvert, Face const &f,
    [[maybe_unused]] std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
#ifdef DEBUG_CANONICAL_LIMITED_V2
  os << "Entry " << StringGroup(GRP) << " " << StringFace(f) << "\n";
#endif
#ifdef DEBUG_CANONICAL_LIMITED
  os << "CAN_LIM: Beginning of CanonicalImageGeneralDualDesc method_choice="
     << method_choice << "\n";
  WriteGroup(os, GRP);
  os << "CAN_LIM: f=" << StringFace(f) << "\n";
#endif
  if (method_choice == CANONIC_STRATEGY__CANONICAL_IMAGE) {
    std::pair<Face, TintGroup> pair = GRP.CanonicalImageOrbitSize(f);
#ifdef DEBUG_CANONICAL_LIMITED
    os << "CAN_LIM: After CanonicalImageOrbitSize\n";
#endif
    return recConvert.ConvertFaceOrbitSize(pair);
  }
  if (method_choice == CANONIC_STRATEGY__STORE) {
    std::pair<Face, TintGroup> pair = GRP.StoreCanonicalImageOrbitSize(f);
#ifdef DEBUG_CANONICAL_LIMITED
    os << "CAN_LIM: After StoreCanonicalImageOrbitSize\n";
#endif
    return recConvert.ConvertFaceOrbitSize(pair);
  }
  if (method_choice == CANONIC_STRATEGY__INITIAL_TRIV) {
    Face f_red = GRP.CanonicalImageInitialTriv(f);
#ifdef DEBUG_CANONICAL_LIMITED
    os << "CAN_LIM: After CanonicalImageInitialTriv\n";
#endif
    return f_red;
  }
  if (method_choice == CANONIC_STRATEGY__INITIAL_TRIV_LIMITED1) {
    try {
      Face f_red = GRP.CanonicalImageInitialTrivLimited(f, LIMIT_INITIAL_TRIV);
#ifdef DEBUG_CANONICAL_LIMITED
      os << "CAN_LIM: After CanonicalImageInitialTrivLimited\n";
#endif
      return f_red;
    } catch (...) {
      std::cerr << "Catching some exception from CanonicalImageInitialTrivLimited\n";
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
    for (size_t i = 0; i < n_red; i++)
      f_red[i] = f[i];
    vf_red.push_back(f_red);
  }
  return vf_red;
}

Face face_reduction(Face const &f, size_t n_red) {
  Face f_red(n_red);
  for (size_t i = 0; i < n_red; i++)
    f_red[i] = f[i];
  return f_red;
}

void set_face_partial(Face &f_out, Face const &f_in, size_t const &n_red) {
  for (size_t i = 0; i < n_red; i++)
    f_out[i] = f_in[i];
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
    LGenExt.push_back(eGenExt);
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

template <typename Tgroup>
int GetCanonicalizationMethod_MPI(boost::mpi::communicator &comm,
                                  vectface const &vf, Tgroup const &GRP,
                                  std::ostream &os) {
  int n_proc = comm.size();
  std::vector<int> list_considered = GetPossibleCanonicalizationMethod(GRP);
  int64_t miss_val = std::numeric_limits<int64_t>::max();
  int64_t upper_limit_local = miss_val;
  int64_t upper_limit_global = miss_val;
  int chosen_method = CANONIC_STRATEGY__DEFAULT;
  int64_t effective_upper_limit = miss_val;
  std::vector<int64_t> V_runtime;
  for (auto &method : list_considered) {
    if (upper_limit_local != miss_val) {
      // That is a tolerance. If a method gives 5 times worse locally, then that
      // is enough to discard it.
      effective_upper_limit = 5 * upper_limit_local;
    }
    int64_t runtime_local =
        time_evaluation_can_method(method, vf, GRP, effective_upper_limit, os);
    if (runtime_local < upper_limit_local) {
      upper_limit_local = runtime_local;
    }
    boost::mpi::all_gather<int64_t>(comm, runtime_local, V_runtime);
    int64_t runtime_global = 0;
    for (int i_proc = 0; i_proc < n_proc; i_proc++) {
      if (runtime_global != miss_val) {
        int64_t runtime = V_runtime[i_proc];
        if (runtime == miss_val) {
          runtime_global = miss_val;
        } else {
          runtime_global += runtime;
        }
      }
    }
    if (runtime_global < upper_limit_global) {
      chosen_method = method;
      upper_limit_global = runtime_global;
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
  std::vector<TintGroup> ListPossOrbSize = ListOrbitFaceOrbitsize.ListPossOrbsize;
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
    std::vector<TintGroup> ListPossOrbSize = ListOrbitFaceOrbitsize.ListPossOrbsize;
    bank.InsertEntry(
        std::move(ePair.first),
        {std::move(GrpConj), std::move(ListPossOrbSize), std::move(ListFaceO)});
  } else {
    std::pair<MyMatrix<T>, TripleStore<Tgroup>> eP = GetCanonicalInformation(
        EXT, WMat, TheGRPrelevant, ListOrbitFaceOrbitsize, os);
    bank.InsertEntry(std::move(eP.first), std::move(eP.second));
  }
}

template <typename Tgroup, typename Tface_orbitsize>
vectface
OrbitSplittingListOrbitGen(const Tgroup &GRPbig, const Tgroup &GRPsma,
                           Tface_orbitsize const &ListFaceOrbitsize,
                           PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
                           std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
  TintGroup ordGRPbig = GRPbig.size();
  TintGroup ordGRPsma = GRPsma.size();
  if (ordGRPbig == ordGRPsma) {
    return ListFaceOrbitsize.GetListFaces();
  }
  TintGroup index = ordGRPbig / ordGRPsma;
  std::map<std::string, TintGroup> TheMap;
  TheMap["groupsize_big"] = ordGRPbig;
  TheMap["groupsize_sma"] = ordGRPsma;
  TheMap["index"] = index;
  TheMap["n_orbit"] = ListFaceOrbitsize.size();
  std::string method_split =
      HeuristicEvaluation(TheMap, AllArr.OrbitSplitTechnique);
  return OrbitSplittingListOrbit_spec(GRPbig, GRPsma, ListFaceOrbitsize,
                                      method_split, os);
}

template <typename Tbank, typename T, typename Tgroup, typename Tidx_value>
vectface getdualdesc_in_bank(Tbank &bank, MyMatrix<T> const &EXT,
                             WeightMatrix<true, T, Tidx_value> const &WMat,
                             Tgroup const &GRP,
                             PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
                             std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  using Tidx = typename Telt::Tidx;
  std::pair<MyMatrix<T>, std::vector<Tidx>> ePair =
      CanonicalizationPolytopePair<T, Tidx, Tidx_value>(EXT, WMat, os);
  const TripleStore<Tgroup> &RecAns = bank.GetDualDesc(ePair.first);
  if (RecAns.ListFace.size() == 0) {
    return vectface(0);
  }
#ifdef DEBUG_RECURSIVE_DUAL_DESC
  os << "REC_DD: Finding a matching entry in the bank\n";
#endif
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
  FaceOrbitsizeTableContainer<TintGroup> fotc(RecAns.ListPossOrbsize, n,
                                              std::move(ListReprTrans));
  return OrbitSplittingListOrbitGen(GrpConj, GRP, fotc, AllArr, os);
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

template <typename Tgroup>
Tgroup StabilizerUsingOrbSize_OnSets(
    Tgroup const &GRP, std::pair<Face, typename Tgroup::Tint> const &pair) {
  using Tint = typename Tgroup::Tint;
  Face const &f = pair.first;
  Tint orbSize = pair.second;
  // If the orbit is equal to the full group then the stabilizer is trivial
  if (GRP.size() == orbSize) {
    size_t len = f.count();
    return Tgroup(len);
  }
  // If the orbit is of size 1 then the stabilizer is the full group. However,
  // we do not test for this because the case is incredibly rare and the
  // Stabilizer_OnSets would anyway short circuit that.
  //
  Tgroup Stab = GRP.Stabilizer_OnSets(f);
  return ReducedGroupAction(Stab, f);
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
  using Text_int = typename SubsetRankOneSolver<T>::Tint;
  using Telt = typename Tgroup::Telt;
  const MyMatrix<T> &EXT;
  const MyMatrix<Text_int> &EXT_int;
  const Tgroup &GRP;
  using Torbsize = uint32_t;
  using Tidx = typename Telt::Tidx;
  int nbRow;
  int nbCol;
  size_t delta;
  FaceOrbsizeContainer<Tgroup, Torbsize, Tidx> foc;
  int the_method; // we use an indifferent name because we also have it for
                  // DatabaseRepr
private:
  UNORD_SET<size_t, std::function<size_t(Tidx_orbit)>,
            std::function<bool(Tidx_orbit, Tidx_orbit)>>
      DictOrbit;
  std::map<size_t, std::vector<Tidx_orbit>> CompleteList_SetUndone;
  std::ostream &os;
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

  void InsertEntryDatabase(std::pair<Face, Tint> const &face_pair,
                           bool const &status, size_t const &pos) {
    Face const &face = face_pair.first;
    Tint const &orbSize = face_pair.second;
#ifdef TRACK_DATABASE
    os << "REC_DD: InsertEntryDatabase |EXT|=" << nbRow << "/" << nbCol
       << " status=" << status << " face.size=" << face.size()
       << " face.count=" << face.count() << " orbSize=" << orbSize
       << " pos=" << pos << "\n";
#endif
    if (!status) {
      size_t len = face.count();
      CompleteList_SetUndone[len].push_back(pos);
    }
    foc.Counts_InsertOrbit(status, orbSize);
  }
  DatabaseCanonic(MyMatrix<T> const &_EXT, MyMatrix<Text_int> const &_EXT_int,
                  Tgroup const &_GRP, std::ostream &_os)
      : EXT(_EXT), EXT_int(_EXT_int), GRP(_GRP), foc(GRP), os(_os) {
    the_method = std::numeric_limits<int>::max();

    /* TRICK 6: The UNORD_SET only the index and this saves in memory usage. */
    n_act = GRP.n_act();
    delta = foc.delta;
    n_act_div8 = (n_act + 7) / 8;
    nbRow = EXT.rows();
    nbCol = EXT.cols();
    //    os << "DatabaseCanonic constructor |EXT|=" << nbRow << "/" << nbCol <<
    //    " delta=" << delta << " n_act=" << n_act << "\n";
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
      os << "REC_DD: hash=" << hash << "\n";
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
    std::function<bool(Tidx_orbit, Tidx_orbit)> fctEqual =
        [&](Tidx_orbit idx1, Tidx_orbit idx2) -> bool {
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
    DictOrbit = UNORD_SET<size_t, std::function<size_t(Tidx_orbit)>,
                          std::function<bool(Tidx_orbit, Tidx_orbit)>>(
        {}, fctHash, fctEqual);
  }
  ~DatabaseCanonic() {}
  void clear() {
    foc.clear();
    DictOrbit.clear();
    CompleteList_SetUndone.clear();
  }
  int evaluate_method_serial(vectface const &vf) const {
    return GetCanonicalizationMethod_Serial(vf, GRP, os);
  }
  int evaluate_method_mpi(boost::mpi::communicator &comm,
                          vectface const &vf) const {
    return GetCanonicalizationMethod_MPI(comm, vf, GRP, os);
  }
  Face operation_face(Face const &face) {
    return CanonicalImageGeneralDualDesc(the_method, GRP, foc.recConvert, face,
                                         os);
  }
  int convert_string_method(std::string const &choice) const {
    if (choice == "canonic")
      return CANONIC_STRATEGY__CANONICAL_IMAGE;
    if (choice == "canonic_initial_triv")
      return CANONIC_STRATEGY__INITIAL_TRIV;
    if (choice == "canonic_initial_triv_limited1")
      return CANONIC_STRATEGY__INITIAL_TRIV_LIMITED1;
    if (choice == "store")
      return CANONIC_STRATEGY__STORE;
    std::cerr << "The value of choice is not an allowed one\n";
    std::cerr
        << "Only possibilities are canonic, canonic_initial_triv, store\n";
    throw TerminalException{1};
  }
  bool use_f_insert_pair() {
    if (the_method == CANONIC_STRATEGY__CANONICAL_IMAGE) {
      return true;
    }
    if (the_method == CANONIC_STRATEGY__STORE) {
      return true;
    }
    if (the_method == CANONIC_STRATEGY__INITIAL_TRIV) {
      return false;
    }
    if (the_method == CANONIC_STRATEGY__INITIAL_TRIV_LIMITED1) {
      return false;
    }
    std::cerr << "The value of the_method was not correctly set\n";
    throw TerminalException{1};
  }
  int get_default_strategy() { return CANONIC_STRATEGY__DEFAULT; }
  FaceOrbitsizeTableContainer<Tint> GetListFaceOrbitsize() {
    DictOrbit.clear();
    CompleteList_SetUndone.clear();
    return foc.GetListFaceOrbitsize();
  }
  void FuncInsert(Face const &face_can) {
#ifdef TRACK_DATABASE
    if (face_can.size() != static_cast<size_t>(nbRow)) {
      std::cerr << "We have |face_can|=" << face_can.size()
                << " but nbRow=" << nbRow << "\n";
      os << "REC_DD: We have |face_can|=" << face_can.size() << " but nbRow=" << nbRow
         << "\n";
      throw TerminalException{1};
    }
#endif
    // The face should have been canonicalized beforehand.
#ifdef TRACK_DATABASE
    os << "REC_DD: FuncInsert face_can.size=" << face_can.size()
       << " face_can.count=" << face_can.count() << "\n";
#endif
    foc.InsertListOrbitFace(face_can);
    DictOrbit.insert(foc.nbOrbit);
#ifdef TRACK_DATABASE
    os << "REC_DD: DictOrbit.size=" << DictOrbit.size()
       << " foc.nbOrbit=" << foc.nbOrbit << "\n";
#endif
    if (DictOrbit.size() == foc.nbOrbit) {
#ifdef TRACK_DATABASE
      os << "REC_DD: FuncInsert : Already present exit\n";
#endif
      // Insertion did not raise the count
      // and so it was already present
      return;
    }
    /* TRICK 8: The insertion yield something new. So now we compute the
     * expensive stabilizer */
    Tint orbSize = GRP.OrbitSize_OnSets(face_can);
#ifdef TRACK_DATABASE
    os << "REC_DD: FuncInsert : New orbSize=" << orbSize << "\n";
#endif
    foc.FinishWithOrbSizeAssignation(orbSize);
#ifdef CHECK_INSERT
    Tgroup eStab = GRP.Stabilizer_OnSets(face_can);
    os << "REC_DD: FuncInsert: Inserting a face of size |face_can|=" << face_can.count()
       << " |eStab|=" << eStab.size() << " f=" << StringFace(face_can) << "\n";
#endif
    InsertEntryDatabase({face_can, orbSize}, false, foc.nbOrbit);
  }
  void FuncInsertPair(Face const &face_orbsize) {
#ifdef TRACK_DATABASE
    if (face_orbsize.size() != delta) {
      std::cerr << "We have |face_orbsize|=" << face_orbsize.size()
                << " but delta=" << delta << "\n";
      os << "REC_DD: We have |face_orbsize|=" << face_orbsize.size()
         << " but delta=" << delta << "\n";
      throw TerminalException{1};
    }
#endif
#ifdef TRACK_DATABASE
    os << "REC_DD: FuncInsertPair face_orbsize.size=" << face_orbsize.size()
       << " face_orbsize.count=" << face_orbsize.count() << "\n";
#endif
    // The face should have been canonicalized beforehand and also contains the
    // orbits
    foc.InsertListOrbitFaceComplete(face_orbsize);
    DictOrbit.insert(foc.nbOrbit);
    if (DictOrbit.size() == foc.nbOrbit) {
#ifdef TRACK_DATABASE
      os << "REC_DD: FuncInsertPair : Already present exit\n";
#endif
      // Insertion did not raise the count
      // and so it was already present
      return;
    }
    std::pair<Face, Tint> pair = foc.recConvert.ConvertFace(face_orbsize);
#ifdef TRACK_DATABASE
    os << "REC_DD: FuncInsertPair : New orbSize(pair.second)=" << pair.second << "\n";
#endif
#ifdef CHECK_INSERT
    os << "REC_DD: FuncInsertPair: Inserting a face of size |face_can|="
       << pair.first.count() << " |O|=" << pair.second << "\n";
#endif
    InsertEntryDatabase(pair, false, foc.nbOrbit);
  }
  void FuncPutOrbitAsDone(size_t const &i_orb) {
#ifdef TRACK_DATABASE
    os << "REC_DD: FuncPutOrbitAsDone : nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    os << "REC_DD: FuncPutOrbitAsDone : CompleteList_SetUndone\n";
    for (auto &kv : CompleteList_SetUndone) {
      os << "kv.first=" << kv.first << " |kv.second|=" << kv.second.size()
         << "\n";
    }
    os << "REC_DD: FuncPutOrbitAsDone : i_orb=" << i_orb << "\n";
#endif
    std::pair<Face, Tint> eEnt = foc.RetrieveListOrbitEntry(i_orb);
    size_t len = eEnt.first.count();
#ifdef TRACK_DATABASE
    os << "REC_DD: FuncPutOrbitAsDone : len=" << len << "\n";
#endif
    /* TRICK 1: We copy the last element in first position to erase it and then
     * pop_back the vector. */
    std::vector<Tidx_orbit> &V = CompleteList_SetUndone[len];
#ifdef TRACK_DATABASE
    os << "REC_DD: FuncPutOrbitAsDone : |V|=" << V.size() << "\n";
#endif
    if (V.size() == 1) {
      CompleteList_SetUndone.erase(len);
    } else {
#ifdef TRACK_DATABASE
      if (V.size() == 0) {
        std::cerr << "We have V.size = 0 which is not allowed here\n";
        throw TerminalException{1};
      }
#endif
      V[0] = V[V.size() - 1];
      V.pop_back();
      if (2 * V.size() < V.capacity()) {
        V.shrink_to_fit();
      }
    }
    foc.Counts_SetOrbitDone(eEnt.second);
  }
  DataFacet<T, Tgroup> FuncGetMinimalUndoneOrbit() {
    for (auto &eEnt : CompleteList_SetUndone) {
      size_t len = eEnt.second.size();
      if (len > 0) {
        /* TRICK 1: Take the first element in the vector. This first element
           will remain
           in place but the vector may be extended without impacting this first
           entry. */
        size_t pos = eEnt.second[0];
        std::pair<Face, Tint> pair = foc.RetrieveListOrbitEntry(pos);
        Face const &f = pair.first;
        Tgroup StabRed = StabilizerUsingOrbSize_OnSets(GRP, pair);
        return {pos, f, FlippingFramework<T>(EXT, EXT_int, f, os), GRP,
                StabRed};
      }
    }
    std::cerr << "Failed to find an undone orbit\n";
    throw TerminalException{1};
  }
  void InsertListOrbitEntry(Face const &f, const size_t &i_orbit) {
    foc.InsertListOrbitFaceComplete(f);
    DictOrbit.insert(i_orbit);
  }

private:
  struct IteratorIndexType {
  private:
    std::map<size_t, std::vector<Tidx_orbit>>::const_iterator iter;
    size_t pos;

  public:
    IteratorIndexType(
        std::map<size_t, std::vector<Tidx_orbit>>::const_iterator iter,
        size_t pos)
        : iter(iter), pos(pos) {}
    size_t operator*() const { return iter->second[pos]; }
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
    Face operator*() const { return foc.RetrieveListOrbitFace(*iter); }
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

    bool operator!=(const IteratorFaceType &x) const { return iter != x.iter; }
    bool operator==(const IteratorFaceType &x) const { return iter == x.iter; }
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
  using Text_int = typename SubsetRankOneSolver<T>::Tint;
  using Telt = typename Tgroup::Telt;
  const MyMatrix<T> &EXT;
  const MyMatrix<Text_int> &EXT_int;
  const Tgroup &GRP;
  using Torbsize = uint32_t;
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
  std::ostream &os;

public:
  DatabaseRepr() = delete;
  DatabaseRepr(const DatabaseRepr<T, Tint, Tgroup, Frepr, Forbitsize, Finv> &) =
      delete;
  DatabaseRepr(DatabaseRepr<T, Tint, Tgroup, Frepr, Forbitsize, Finv> &&) =
      delete;
  DatabaseRepr &operator=(
      const DatabaseRepr<T, Tint, Tgroup, Frepr, Forbitsize, Finv> &) = delete;

  void InsertEntryDatabase(std::pair<Face, Tint> const &face_pair,
                           bool const &status, size_t const &pos) {
    Face const &face = face_pair.first;
    Tint const &orbSize = face_pair.second;
    size_t len = face.count();
    size_t eInv = f_inv(face);
    if (status) {
      CompleteList_SetDone[len][eInv].push_back(pos);
    } else {
      CompleteList_SetUndone[len][eInv].push_back(pos);
    }
    foc.Counts_InsertOrbit(status, orbSize);
  }
  DatabaseRepr(MyMatrix<T> const &_EXT, MyMatrix<Text_int> const &_EXT_int,
               Tgroup const &_GRP, Frepr f_repr, Forbitsize f_orbitsize,
               Finv f_inv, std::ostream &_os)
      : EXT(_EXT), EXT_int(_EXT_int), GRP(_GRP), foc(GRP), f_repr(f_repr),
        f_orbitsize(f_orbitsize), f_inv(f_inv), os(_os) {
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
  int evaluate_method_serial([[maybe_unused]] vectface const &vf) const {
    return REPR_STRATEGY__DEFAULT;
  }
  int evaluate_method_mpi([[maybe_unused]] boost::mpi::communicator &comm,
                          [[maybe_unused]] vectface const &vf) const {
    return REPR_STRATEGY__DEFAULT;
  }
  Face operation_face(Face const &eFlip) { return eFlip; }
  int convert_string_method(std::string const &choice) const {
    if (choice == "default")
      return REPR_STRATEGY__DEFAULT;
    std::cerr << "The value of choice is not an allowed one\n";
    std::cerr << "Only possibility is default\n";
    throw TerminalException{1};
  }
  bool use_f_insert_pair() { return false; }
  int get_default_strategy() { return REPR_STRATEGY__DEFAULT; }
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
    foc.FinishWithOrbSizeAssignation(orbSize);
    InsertEntryDatabase({face_i, orbSize}, false, foc.nbOrbit);
  }
  void FuncInsertPair(Face const &face) {
    Face f_red(nbRow);
    for (int i = 0; i < nbRow; i++) {
      f_red[i] = face[i];
    }
    FuncInsert(f_red);
  }
  void FuncPutOrbitAsDone(size_t const &iOrb) {
    std::pair<Face, Tint> eEnt = foc.RetrieveListOrbitEntry(iOrb);
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
  DataFacet<T, Tgroup> FuncGetMinimalUndoneOrbit() {
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
    std::pair<Face, Tint> pair = foc.RetrieveListOrbitEntry(pos);
    Face const &f = pair.first;
    Tgroup StabRed = StabilizerUsingOrbSize_OnSets(GRP, pair);
    return {pos, f, FlippingFramework<T>(EXT, EXT_int, f, os), GRP, StabRed};
  }
  void InsertListOrbitEntry(Face const &f,
                            [[maybe_unused]] const size_t &i_orbit) {
    foc.InsertListOrbitFaceComplete(f);
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
        : iter1(iter1), iter1_end(iter1_end), iter2(iter2), pos(pos) {}
    size_t operator*() const { return iter2->second[pos]; }
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
    Face operator*() const { return foc.RetrieveListOrbitFace(*iter); }
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
    bool operator!=(const IteratorFaceType &x) const { return iter != x.iter; }
    bool operator==(const IteratorFaceType &x) const { return iter == x.iter; }
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
  std::string eFileEXT, eFileGRP, eFileNB, eFileFB, eFileFF, eFileMethod;
  /* TRICK 7: Using separate files for faces and status allow us to gain
     locality. The faces are written one by one while the access to status is
     random */
  bool SavingTrigger;
  bool NeedToFlush;
  bool AdvancedTerminationCriterion;
  std::ostream &os;
  size_t delta;
  std::string strPresChar;
  HumanTime time;

public:
  // method encodes the algorithm used for the database and essentially applies
  // only to the canonic.
  // ---If method file is absent then we assume it was computed with the
  // default.
  // ---Otherwise we read it.
  void write_method(std::string const &eFileMethod, int const &method) const {
    std::ofstream os_file(eFileMethod);
    os_file << method;
  }
  int read_method(std::string const &eFileMethod) const {
#ifdef TRACK_DATABASE
    os << "SavingTrigger=" << SavingTrigger << " eFileMethod=" << eFileMethod
       << "\n";
#endif
    if (SavingTrigger) {
#ifdef TRACK_DATABASE
      os << "Running the save system\n";
#endif
      if (!IsExistingFile(eFileMethod)) {
#ifdef TRACK_DATABASE
        os << "The file does not exists\n";
#endif
        int the_method = bb.get_default_strategy();
        write_method(eFileMethod, the_method);
        return the_method;
      } else {
#ifdef TRACK_DATABASE
        os << "The file exists\n";
#endif
        std::ifstream is_file(eFileMethod);
        int method;
        is_file >> method;
        return method;
      }
    } else {
      return bb.get_default_strategy();
    }
  }
  bool is_database_present() const {
    if (IsExistingFile(eFileEXT) == false) {
      return false;
    }
    // verify that EXT file is same as bb.EXT
    MyMatrix<T> EXT = ReadMatrixFile<T>(eFileEXT);
    if (EXT == bb.EXT) {
      return true;
    }
    // else database got changed e.g. due to method change
    // remove it
    // optional future function: check for equivalence and convert
    if (SavingTrigger) {
#ifdef TRACK_DATABASE
      os << "Database got changed, removing old one\n";
#endif
      RemoveFileIfExist(eFileNB);
      RemoveFileIfExist(eFileFB);
      RemoveFileIfExist(eFileFF);
      RemoveFileIfExist(eFileEXT);
      RemoveFileIfExist(eFileGRP);
      RemoveFileIfExist(eFileMethod);
    }
    return false;
  }
  DatabaseOrbits() = delete;
  DatabaseOrbits(const DatabaseOrbits<TbasicBank> &) = delete;
  DatabaseOrbits(DatabaseOrbits<TbasicBank> &&) = delete;
  DatabaseOrbits &operator=(const DatabaseOrbits<TbasicBank> &) = delete;
  void print_status() const {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "REC_DD: Status : orbit=(" << bb.foc.nbOrbit << "," << bb.foc.nbOrbitDone
       << "," << (bb.foc.nbOrbit - bb.foc.nbOrbitDone) << ") facet=("
       << bb.foc.TotalNumber << "," << (bb.foc.TotalNumber - bb.foc.nbUndone)
       << "," << bb.foc.nbUndone << ")\n\n";
#endif
  }
  DatabaseOrbits(TbasicBank &bb, const std::string &MainPrefix,
                 const bool &_SavingTrigger,
                 const bool &_AdvancedTerminationCriterion, std::ostream &os)
      : CritSiz(bb.EXT.cols() - 2), bb(bb), SavingTrigger(_SavingTrigger),
        AdvancedTerminationCriterion(_AdvancedTerminationCriterion), os(os) {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "REC_DD: MainPrefix=" << MainPrefix << "\n";
#endif
    eFileEXT = MainPrefix + ".ext";
    eFileGRP = MainPrefix + ".grp";
    eFileNB = MainPrefix + ".nb";
    eFileFB = MainPrefix + ".fb";
    eFileFF = MainPrefix + ".ff";
    eFileMethod = MainPrefix + ".method";
    strPresChar = "|EXT|=" + std::to_string(bb.nbRow) + "/" +
                  std::to_string(bb.nbCol) +
                  " |GRP|=" + std::to_string(bb.GRP.size());
    delta = bb.delta;
    NeedToFlush = true;
    int val = read_method(eFileMethod);
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "REC_DD: read_method val=" << val << "\n";
#endif
    bb.the_method = val;
    if (SavingTrigger && !is_database_present()) {
      if (!FILE_IsFileMakeable(eFileEXT)) {
        std::cerr << "Error in DatabaseOrbits: File eFileEXT=" << eFileEXT
                  << " is not makeable\n";
        throw TerminalException{1};
      }
      initial_writes();
    }
  }
  void initial_writes() {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "REC_DD: Creating the initial files (NB, FB, FF) with zero state\n";
#endif
    FileNumber fn(eFileNB, true);
    FileBool fb(eFileFB);
    FileFace ff(eFileFF, bb.delta);
    std::vector<uint8_t> V_empty; // empty write, maybe useless.
    fn.setval(0);
    ff.direct_write(V_empty);
    fb.direct_write(V_empty);
    write_method(eFileMethod, bb.the_method);
    std::ofstream os_grp(eFileGRP);
    os_grp << bb.GRP;
    WriteMatrixFile(eFileEXT, bb.EXT);
  }
  size_t preload_nb_orbit() const {
    if (SavingTrigger && is_database_present()) {
      FileNumber fn(eFileNB, false);
      return fn.getval();
    }
    return 0;
  }
  void LoadDatabase() {
    if (is_database_present()) {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
      os << "REC_DD: Opening existing files (NB, FB, FF)\n";
#endif
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      MicrosecondTime time;
#endif
      FileNumber fn(eFileNB, false);
      size_t n_orbit = fn.getval();
#ifdef DEBUG_RECURSIVE_DUAL_DESC
      os << "REC_DD: Loading database with n_orbit=" << n_orbit << "\n";
#endif
      FileBool fb(eFileFB, n_orbit);
      FileFace ff(eFileFF, bb.delta, n_orbit);
      for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
        Face f = ff.getface(i_orbit);
        std::pair<Face, Tint> eEnt = bb.foc.recConvert.ConvertFace(f);
        bool status = fb.getbit(i_orbit);
        bb.InsertListOrbitEntry(f, i_orbit);
        bb.InsertEntryDatabase(eEnt, status, i_orbit);
      }
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|Loading Database|=" << time << "\n";
#endif
    } else {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
      os << "REC_DD: No database present\n";
#endif
    }
    print_status();
  }
  vectface ReadDatabase(size_t const &n_read) const {
    vectface vfo(bb.delta + 1);
    if (is_database_present()) {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
      os << "REC_DD: Opening existing files (NB, FB, FF)\n";
#endif
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      MicrosecondTime time;
#endif
      FileNumber fn(eFileNB, false);
      size_t n_orbit = fn.getval();
#ifdef DEBUG_RECURSIVE_DUAL_DESC
      os << "REC_DD: Reading database with n_orbit=" << n_orbit << "\n";
#endif
      FileBool fb(eFileFB, n_orbit);
      FileFace ff(eFileFF, bb.delta, n_orbit);
      size_t n_read_eff = std::min(n_read, n_orbit);
      for (size_t i_orbit = 0; i_orbit < n_read_eff; i_orbit++) {
        Face f = ff.getface(i_orbit);
        bool status = fb.getbit(i_orbit);
        Face f_insert(bb.delta + 1);
        for (size_t u = 0; u < bb.delta; u++) {
          f_insert[u] = f[u];
        }
        f_insert[bb.delta] = status;
        vfo.push_back(f_insert);
      }
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|Reading Database|=" << time << "\n";
#endif
    } else {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
      os << "REC_DD: No database present\n";
#endif
    }
    return vfo;
  }
  vectface get_runtime_testcase() const {
    size_t n_orbit = preload_nb_orbit();
    size_t n_target = 100;
    int nbRow = bb.nbRow;
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "REC_DD: get_runtime_testcase n_orbit=" << n_orbit
       << " n_target=" << n_target << " nbRow=" << nbRow << "\n";
#endif
    if (n_orbit == 0) {
      vectface vf(nbRow);
      for (size_t i = 0; i < n_target; i++) {
        Face f = RandomFace(nbRow);
        vf.push_back(f);
      }
      return vf;
    } else {
      vectface vfo = ReadDatabase(n_target);
      return vectface_reduction(vfo, nbRow);
    }
  }
  int determine_action_database(std::string const &choice) {
    if (choice == "load")
      return DATABASE_ACTION__SIMPLE_LOAD;
    if (choice == "guess")
      return DATABASE_ACTION__GUESS;
    int choice_i = bb.convert_string_method(choice);
    if (bb.the_method == choice_i)
      return DATABASE_ACTION__SIMPLE_LOAD;
    return DATABASE_ACTION__RECOMPUTE_AND_SHUFFLE;
  }
  void set_method(int const &the_method) {
    bb.the_method = the_method;
    if (SavingTrigger) {
      write_method(eFileMethod, the_method);
    }
  }
  void DirectAppendDatabase(vectface &&vf) {
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
    for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
      Face f = vf[i_orbit];
      Face f_red(bb.delta);
      for (size_t u = 0; u < bb.delta; u++) {
        bool val = f[u];
        f_red[u] = val;
        if (SavingTrigger) {
          setbit_vector(ListOrbit_ff, pos_ff, val);
          pos_ff++;
        }
      }
      bool status = f[bb.delta];
      if (SavingTrigger) {
        setbit_vector(V_status, i_orbit, status);
      }
      std::pair<Face, Tint> eEnt = bb.foc.recConvert.ConvertFace(f_red);
      bb.InsertListOrbitEntry(f_red, i_orbit);
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
    print_status();
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
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "REC_DD: Clean closing of the DatabaseOrbits\n";
#endif
  }
  void flush() const {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "REC_DD: Doing the flushing operation\n";
#endif
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
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
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
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
      RemoveFileIfExist(eFileMethod);
    }
    return bb.GetListFaceOrbitsize();
  }
  void FuncInsert(Face const &face) { bb.FuncInsert(face); }
  void FuncInsertPair(Face const &face) { bb.FuncInsertPair(face); }
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
  DataFacet<T, Tgroup> FuncGetMinimalUndoneOrbit() {
    DataFacet<T, Tgroup> data = bb.FuncGetMinimalUndoneOrbit();
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "REC_DD: " << strPresChar << " Considering orbit " << data.SelectedOrbit
       << " |inc|=" << data.eInc.count() << " |stab|=" << data.Stab.size()
       << "\n";
#endif
    return data;
  }
  bool GetTerminationStatus() const {
    auto get_val = [&]() -> bool {
      if (bb.foc.nbOrbitDone > 0) {
        if (bb.foc.nbUndone <= CritSiz) {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
          os << "REC_DD: Termination by classic Balinski criterion nbUndone="
             << bb.foc.nbUndone << "\n";
#endif
          return true;
        }
        Face eSetUndone = ComputeIntersectionUndone();
        if (eSetUndone.count() > 0) {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
          os << "REC_DD: Termination by linear programming criterion |eSetUndone|="
             << eSetUndone.count() << "\n";
#endif
          return true;
        }
      }
      if (AdvancedTerminationCriterion)
        return EvaluationConnectednessCriterion_Serial(bb, os);
      return false;
    };
    if (get_val()) {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
      os << "REC_DD: End of computation, nbObj=" << bb.foc.TotalNumber
         << " |EXT|=" << bb.nbRow << "/" << bb.nbCol
         << " time=" << time.const_eval() << "\n";
#endif
      return true;
    }
    return false;
  }
  UndoneOrbitInfo<Tint> GetTerminationInfo() const {
    return {bb.foc.nbOrbitDone, bb.foc.nbUndone, ComputeIntersectionUndone()};
  }
};

template <typename T, typename Tidx_value> struct LazyWMat {
public:
  LazyWMat(const MyMatrix<T> &EXT, std::ostream &os)
      : EXT(EXT), WMat(os), HaveWMat(false), os(os) {}
  WeightMatrix<true, T, Tidx_value> &GetWMat() {
    if (HaveWMat)
      return WMat;
    WMat = GetWeightMatrix<T, Tidx_value>(EXT, os);
    WMat.ReorderingSetWeight();
    HaveWMat = true;
    return WMat;
  }

private:
  const MyMatrix<T> &EXT;
  WeightMatrix<true, T, Tidx_value> WMat;
  bool HaveWMat;
  std::ostream &os;
};

template <typename Tint, typename T, typename Tgroup>
std::map<std::string, Tint>
ComputeInitialMap(const MyMatrix<T> &EXT, const Tgroup &GRP,
                  PolyHeuristicSerial<typename Tgroup::Tint> const& AllArr) {
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

template <typename Tgroup>
void CheckTermination(PolyHeuristicSerial<typename Tgroup::Tint> const& AllArr) {
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

template <typename Tbank, typename T, typename Tgroup, typename Tidx_value,
          typename TbasicBank, typename Finsert, typename Fdd>
void DUALDESC_AdjacencyDecomposition_and_insert(
    Tbank &TheBank, TbasicBank &bb, DataFacet<T, Tgroup> &df,
    PolyHeuristicSerial<typename Tgroup::Tint> & AllArr, Finsert f_insert,
    Fdd & f_dd, std::string const &ePrefix, std::ostream &os) {
  using Tint = typename Tgroup::Tint;
  CheckTermination<Tgroup>(AllArr);
  std::map<std::string, Tint> TheMap =
    ComputeInitialMap<Tint,T,Tgroup>(df.FF.EXT_face, df.Stab, AllArr);
  std::string ansSplit = HeuristicEvaluation(TheMap, AllArr.Splitting);
  if (ansSplit != "split") {
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
    MicrosecondTime time_complete;
    MicrosecondTime time_step;
#endif
    std::string ansProg = AllArr.DualDescriptionProgram.get_eval(TheMap);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
    os << "|ansProg|=" << time_step << "\n";
#endif
    if (df.Stab.size() == 1) {
      auto f_process =
          [&](std::pair<Face, MyVector<T>> const &pair_face) -> void {
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
        MicrosecondTime time_loc;
#endif
        Face eFlipPre = df.FlipFaceIneq(pair_face);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
        os << "|FlipFaceIneq1|=" << time_loc << "\n";
#endif
        Face eFlip = bb.operation_face(eFlipPre);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
        os << "|operation_face1|=" << time_loc << "\n";
#endif
        f_insert(eFlip);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
        os << "|insert1|=" << time_loc << "\n";
#endif
      };
      DirectFacetComputationFaceIneq(df.FF.EXT_face, ansProg, f_process, os);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|DirectFacetComputationFaceIneq|=" << time_step << "\n";
#endif
      AllArr.DualDescriptionProgram.pop(os);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|pop|=" << time_step << "\n";
#endif
    } else {
      vectface TheOutput =
          DirectFacetOrbitComputation(df.FF.EXT_face, df.Stab, ansProg, os);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|TheOutput|=" << time_step << "\n";
      os << "Number of facets being generated=" << TheOutput.size() << "\n";
#endif
      AllArr.DualDescriptionProgram.pop(os);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|pop|=" << time_step << "\n";
#endif
      for (auto &eOrb : TheOutput) {
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
        MicrosecondTime time_loc;
#endif
        Face eFlipPre = df.FlipFace(eOrb);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
        os << "|FlipFace1|=" << time_loc << "\n";
#endif
        Face eFlip = bb.operation_face(eFlipPre);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
        os << "|operation_face1|=" << time_loc << "\n";
#endif
        f_insert(eFlip);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
        os << "|insert1|=" << time_loc << "\n";
#endif
      }
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|Adjacency processing|=" << time_step << "\n";
#endif
    }
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
    os << "|DualDesc+flip+insertion|=" << time_complete << "\n";
#endif
  } else {
    vectface TheOutput =
      f_dd(
            TheBank, df.FF.EXT_face, df.FF.EXT_face_int, df.Stab, TheMap,
            AllArr, ePrefix, os);
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "REC_DD: |outputsize|=" << TheOutput.size() << "\n";
#endif
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
    MicrosecondTime time_full;
#endif
    for (auto &eOrb : TheOutput) {
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      MicrosecondTime time;
#endif
      Face eFlipPre = df.FlipFace(eOrb);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|FlipFace2|=" << time << "\n";
#endif
      Face eFlip = bb.operation_face(eFlipPre);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|operation_face2|=" << time << "\n";
#endif
      f_insert(eFlip);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|insert2|=" << time << "\n";
#endif
    }
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
    os << "|outputtime|=" << time_full << "\n";
#endif
  }
}

template <typename TbasicBank>
void vectface_update_method(vectface &vfo, TbasicBank &bb, [[maybe_unused]] std::ostream &os) {
  size_t n_orbit = vfo.size();
  int nbRow = bb.nbRow;
#ifdef DEBUG_RECURSIVE_DUAL_DESC
  os << "REC_DD: vectface_update_method n_orbit=" << n_orbit << " nbRow=" << nbRow
     << "\n";
  os << "REC_DD: vectface_update_method bb.the_method=" << bb.the_method << "\n";
#endif
  for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
    Face fo = vfo[i_orbit];
    Face f = face_reduction(fo, nbRow);
    Face f_new = bb.operation_face(f);
    set_face_partial(fo, f_new, nbRow);
    vfo.AssignEntry(fo, i_orbit);
  }
}

template <typename Tbank, typename T, typename Tgroup, typename Tidx_value,
          typename TbasicBank, typename Fdd>
FaceOrbitsizeTableContainer<typename Tgroup::Tint>
Kernel_DUALDESC_AdjacencyDecomposition(
    Tbank &TheBank, TbasicBank &bb,
    PolyHeuristicSerial<typename Tgroup::Tint> & AllArr,
    std::string const &ePrefix,
    std::map<std::string, typename Tgroup::Tint> & TheMap,
    Fdd & f_dd, std::ostream &os) {
  DatabaseOrbits<TbasicBank> RPL(bb, ePrefix, AllArr.DD_Saving,
                                 AllArr.AdvancedTerminationCriterion, os);
  // The choice only really makes sense for the canonic, for repr no choice is
  // implied.
  auto set_up = [&]() -> void {
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
    MicrosecondTime time;
#endif
    std::string ansChoiceCanonic =
        HeuristicEvaluation(TheMap, AllArr.ChoiceCanonicalization);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
    os << "|HeuristicEvaluation|=" << time
       << " ansChoiceCanonic=" << ansChoiceCanonic << "\n";
#endif
    int action = RPL.determine_action_database(ansChoiceCanonic);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
    os << "|determine_action_database|=" << time << " action=" << action
       << "\n";
#endif
    auto f_recompute = [&](int const &method) -> void {
      size_t n_orbit = RPL.preload_nb_orbit();
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "n_orbit=" << n_orbit << " |n_orbit|=" << time << "\n";
#endif
      vectface vfo = RPL.ReadDatabase(n_orbit);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "delta=" << bb.delta << " nbRow=" << bb.nbRow << "\n";
      os << "|vfo|=" << vfo.size() << " / " << vfo.get_n()
         << " |ReadDatabase|=" << time << "\n";
#endif
      RPL.set_method(method);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|set_method|=" << time << "\n";
#endif
      vectface_update_method(vfo, bb, os);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "bb.the_method=" << bb.the_method << " |method update|=" << time
         << "\n";
#endif
      RPL.DirectAppendDatabase(std::move(vfo));
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|vfo|=" << vfo.size() << " / " << vfo.get_n()
         << " |DirectAppendDatabase|=" << time << "\n";
#endif
    };
    if (action == DATABASE_ACTION__SIMPLE_LOAD) {
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "Before RPL.LoadDatabase()\n";
#endif
      return RPL.LoadDatabase();
    }
    if (action == DATABASE_ACTION__RECOMPUTE_AND_SHUFFLE) {
      int method = bb.convert_string_method(ansChoiceCanonic);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "Before f_recompute, method=" << method
         << " ansChoiceCanonic=" << ansChoiceCanonic << "\n";
#endif
      return f_recompute(method);
    }
    if (action == DATABASE_ACTION__GUESS) {
      vectface vf = RPL.get_runtime_testcase();
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|get_runtime_testcase|=" << time << "\n";
#endif
      int method = RPL.bb.evaluate_method_serial(vf);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "method=" << method << " |evaluate_method_serial|=" << time << "\n";
#endif
      if (method == bb.the_method) {
        return RPL.LoadDatabase();
      } else {
        return f_recompute(method);
      }
    }
    std::cerr << "Failed to find a matching entry for action=" << action
              << "\n";
    throw TerminalException{1};
  };
  set_up();
  bool use_f_insert_pair = bb.use_f_insert_pair();
#ifdef DEBUG_RECURSIVE_DUAL_DESC
  os << "REC_DD: use_f_insert_pair=" << use_f_insert_pair << "\n";
#endif
  auto f_insert = [&](Face const &face) -> void {
    if (use_f_insert_pair) {
      RPL.FuncInsertPair(face);
    } else {
      RPL.FuncInsert(face);
    }
  };
  if (RPL.FuncNumberOrbit() == 0) {
    std::string ansSamp = HeuristicEvaluation(TheMap, AllArr.InitialFacetSet);
    for (auto &face : DirectComputationInitialFacetSet(bb.EXT, ansSamp, os)) {
      Face face_can = bb.operation_face(face);
      f_insert(face_can);
    }
  }
  while (true) {
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
    MicrosecondTime time;
#endif
    bool test_final = RPL.GetTerminationStatus();
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
    os << "test_final=" << test_final << "\n";
    os << "|GetTerminationStatus|=" << time << "\n";
#endif
    if (test_final)
      break;
    DataFacet<T, Tgroup> df = RPL.FuncGetMinimalUndoneOrbit();
    size_t SelectedOrbit = df.SelectedOrbit;
    // Alternative way is to use CondTempDirectory. BUT
    // For many we actually do not need to have such a construction.
    // Need to think.
    std::string NewPrefix =
        ePrefix + "ADM" + std::to_string(SelectedOrbit) + "_";
    DUALDESC_AdjacencyDecomposition_and_insert<Tbank, T, Tgroup, Tidx_value,
                                               TbasicBank, decltype(f_insert), decltype(f_dd)>(
        TheBank, bb, df, AllArr, f_insert, f_dd, NewPrefix, os);
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
vectface DUALDESC_AdjacencyDecomposition(
    Tbank &TheBank, MyMatrix<T> const &EXT,
    MyMatrix<typename SubsetRankOneSolver<T>::Tint> const &EXT_int,
    Tgroup const &GRP, std::map<std::string, typename Tgroup::Tint> & TheMap,
    PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
    std::string const &ePrefix, std::ostream &os) {
  auto f_dd=[](Tbank & TheBank_i, MyMatrix<T> const& EXT_i,
               MyMatrix<typename SubsetRankOneSolver<T>::Tint> const &EXT_int_i,
               Tgroup const& GRP_i, std::map<std::string, typename Tgroup::Tint> & TheMap_i,
               PolyHeuristicSerial<typename Tgroup::Tint> &AllArr_i,
               std::string const &ePrefix_i, std::ostream &os_i) -> vectface {
    return DUALDESC_AdjacencyDecomposition<Tbank,T,Tgroup,Tidx_value>(TheBank_i, EXT_i, EXT_int_i, GRP_i, TheMap_i, AllArr_i, ePrefix_i, os_i);
  };
  using Tgr = GraphListAdj;
  using Tint = typename Tgroup::Tint;
#ifdef DEBUG_RECURSIVE_DUAL_DESC
  os << "REC_DD: Beginning of DUALDESC_AdjacencyDecomposition\n";
#endif
  CheckTermination<Tgroup>(AllArr);
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  LazyWMat<T, Tidx_value> lwm(EXT, os);
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
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "REC_DD: ansSymm=" << ansSymm << "\n";
#endif
    if (ansSymm == "yes") {
      TheGRPrelevant = GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(
          lwm.GetWMat(), os);
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
      TbasicBank bb(EXT, EXT_int, TheGRPrelevant, os);
      return Kernel_DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup,
                                                    Tidx_value, TbasicBank>(
                                                                            TheBank, bb, AllArr, MainPrefix, TheMap, f_dd, os);
    }
    if (ansChosenDatabase == "repr") {
      WeightMatrix<true, int, Tidx_value> WMat =
          WeightMatrixFromPairOrbits<Tgroup, Tidx_value>(TheGRPrelevant, os);
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
      TbasicBank bb(EXT, EXT_int, TheGRPrelevant, f_repr, f_orbitsize, f_inv,
                    os);
      return Kernel_DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup,
                                                    Tidx_value, TbasicBank>(
                                                                            TheBank, bb, AllArr, MainPrefix, TheMap, f_dd, os);
    }
    std::cerr << "compute_split_or_not: Failed to find a matching entry\n";
    std::cerr << "Authorized values: canonic, repr\n";
    throw TerminalException{1};
  };
  FaceOrbitsizeTableContainer<Tint> ListOrbitFaceOrbitsize =
      compute_decomposition();
  TheMap["time"] = si(start);
  std::string ansBank = HeuristicEvaluation(TheMap, AllArr.BankSave);
#ifdef DEBUG_RECURSIVE_DUAL_DESC
  os << "REC_DD: elapsed_seconds=" << s(start) << " ansBank=" << ansBank
     << " NeedSplit=" << NeedSplit << "\n";
#endif
  if (ansBank == "yes") {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "REC_DD: Before insert_entry_in_bank\n";
#endif
    insert_entry_in_bank(TheBank, EXT, lwm.GetWMat(), TheGRPrelevant,
                         BankSymmCheck, ListOrbitFaceOrbitsize, os);
  }
#ifdef DEBUG_RECURSIVE_DUAL_DESC
  os << "REC_DD: Before return section\n";
#endif
  if (NeedSplit) {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "REC_DD: Before OrbitSplittingListOrbit\n";
#endif
    return OrbitSplittingListOrbitGen(TheGRPrelevant, GRP,
                                      ListOrbitFaceOrbitsize, AllArr, os);
  } else {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "REC_DD: Returning ListOrbitFaces\n";
#endif
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
safe_rational: The safe rational type. Based on int64_t but failing gracefully\n\
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
  ListBoolValues1_doc["InterceptCtrlC"] = "Default: T\n\
If a CtrlC command is thrown then the program will handle it and stop and\n\
leave a usable database that can be rerun afterwards (it seems not to work anymore)\n\
T: Activate the CtrlC mechanism\n\
F: Do not activate the interception of CtrlC";
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
  ListStringValuesH_doc["CommThreadHeuristicFile"] = "Default: unset.heu\n\
The heuristic for choosing when the communication thread is launched.\n\
If set to unset.heu then disabled";
  ListStringValuesH_doc["ChoiceCanonicalizationFile"] = "Default: unset.heu\n\
The heuristic for choosing the canonicalization method used.\n\
If set to unset.heu then disabled";
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
                  const std::string &OutFormat, std::ostream &os) {
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
    using Tint = typename Tgroup::Tint;
    WeightMatrix<true, T, Tidx_value> WMat =
        GetWeightMatrix<T, Tidx_value>(EXT, os);
    WMat.ReorderingSetWeight();
    FaceOrbitsizeTableContainer<Tint> fotc(TheOutput, GRP);
    std::pair<MyMatrix<T>, TripleStore<Tgroup>> eP =
        GetCanonicalInformation(EXT, WMat, GRP, fotc, os);
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
  std::vector<std::string> Ltype{"safe_rational", "rational", "cpp_rational",
                                 "mpq_rational",  "Qsqrt2",   "Qsqrt5",
                                 "RealAlgebraic"};
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
MyMatrix<T> Get_EXT_DualDesc(FullNamelist const &eFull, [[maybe_unused]] std::ostream &os) {
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  std::string EXTfile = BlockDATA.ListStringValues.at("EXTfile");
  IsExistingFileDie(EXTfile);
#ifdef DEBUG_RECURSIVE_DUAL_DESC
  os << "REC_DD: EXTfile=" << EXTfile << "\n";
#endif
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
Tgroup Get_GRP_DualDesc(FullNamelist const &eFull, [[maybe_unused]] std::ostream &os) {
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  std::string GRPfile = BlockDATA.ListStringValues.at("GRPfile");
  IsExistingFileDie(GRPfile);
#ifdef DEBUG_RECURSIVE_DUAL_DESC
  os << "REC_DD: GRPfile=" << GRPfile << "\n";
#endif
  std::ifstream GRPfs(GRPfile);
  Tgroup GRP = ReadGroup<Tgroup>(GRPfs);
  return GRP;
}

bool Get_InterceptCtrlC_statuc(FullNamelist const &eFull, [[maybe_unused]] std::ostream &os) {
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  bool intercept_ctrl_c = BlockDATA.ListBoolValues.at("InterceptCtrlC");
#ifdef DEBUG_RECURSIVE_DUAL_DESC
  os << "REC_DD: intercept_ctrl_c=" << intercept_ctrl_c << "\n";
#endif
  return intercept_ctrl_c;
}

template<typename Tint>
void PrintPolyHeuristicSerial(PolyHeuristicSerial<Tint> const& AllArr, std::ostream & os) {
  os << "REC_DD: OUTfile=" << AllArr.OUTfile << "\n";
  os << "REC_DD: OutFormat=" << AllArr.OutFormat << "\n";
  os << "REC_DD: DeterministicRuntime=" << AllArr.DeterministicRuntime << "\n";
  os << "REC_DD: port=" << AllArr.port << "\n";
  os << "REC_DD: bank_parallelization_method=" << AllArr.bank_parallelization_method << "\n";
  os << "REC_DD: SplittingHeuristicFile\n" << AllArr.Splitting << "\n";
  os << "REC_DD: AdditionalSymmetryHeuristicFile\n" << AllArr.AdditionalSymmetry << "\n";
  os << "REC_DD: DualDescriptionThompsonFile\n" << AllArr.DualDescriptionProgram << "\n";
  os << "REC_DD: MethodInitialFacetSetFile\n" << AllArr.InitialFacetSet << "\n";
  os << "REC_DD: BankSaveHeuristicFile\n" << AllArr.BankSave << "\n";
  os << "REC_DD: CheckDatabaseBank\n" << AllArr.CheckDatabaseBank << "\n";
  os << "REC_DD: ChosenDatabase\n" << AllArr.ChosenDatabase << "\n";
  os << "REC_DD: OrbitSplitTechnique\n" << AllArr.OrbitSplitTechnique << "\n";
  os << "REC_DD: CommThreadHeuristicFile\n" << AllArr.CommThread << "\n";
  os << "REC_DD: ChoiceCanonicalizationFile\n" << AllArr.ChoiceCanonicalization << "\n";
  os << "REC_DD: max_runtime=" << AllArr.max_runtime << "\n";
  os << "REC_DD: DD_Saving=" << AllArr.DD_Saving << "\n";
  os << "REC_DD: DD_Prefix=" << AllArr.DD_Prefix << "\n";
  os << "REC_DD: BANK_Saving=" << AllArr.BANK_Saving << "\n";
  os << "REC_DD: BANK_Prefix=" << AllArr.BANK_Prefix << "\n";
  os << "REC_DD: AdvancedTerminationCriterion=" << AllArr.AdvancedTerminationCriterion << "\n";
  os << "REC_DD: SimpleExchangeScheme=" << AllArr.SimpleExchangeScheme << "\n";
}


template <typename T, typename Tint>
PolyHeuristicSerial<Tint>
Read_AllStandardHeuristicSerial(FullNamelist const &eFull,
                                MyMatrix<T> const &EXTred, std::ostream &os) {
  PolyHeuristicSerial<Tint> AllArr = AllStandardHeuristicSerial<Tint>(os);
  //
  SingleBlock BlockMETHOD = eFull.ListBlock.at("METHOD");
  SingleBlock BlockBANK = eFull.ListBlock.at("BANK");
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  //
  bool BANK_Saving = BlockBANK.ListBoolValues.at("Saving");
  AllArr.BANK_Saving = BANK_Saving;
  //
  std::string BANK_Prefix = BlockBANK.ListStringValues.at("Prefix");
  AllArr.BANK_Prefix = BANK_Prefix;
  //
  std::string OUTfile = BlockDATA.ListStringValues.at("OUTfile");
  AllArr.OUTfile = OUTfile;
  //
  bool DeterministicRuntime =
      BlockDATA.ListBoolValues.at("DeterministicRuntime");
  if (!DeterministicRuntime) {
    unsigned seed = get_random_seed();
    srand(seed);
  }
  //
  std::string OutFormat = BlockDATA.ListStringValues.at("OutFormat");
  AllArr.OutFormat = OutFormat;
  //
  int port_i = BlockDATA.ListIntValues.at("port");
  uint16_t port = port_i;
  AllArr.port = port;
  //
  std::string bank_parallelization_method =
      BlockDATA.ListStringValues.at("bank_parallelization_method");
  AllArr.bank_parallelization_method = bank_parallelization_method;
  //
  SetHeuristic(eFull, "SplittingHeuristicFile", AllArr.Splitting, os);
  SetHeuristic(eFull, "AdditionalSymmetryHeuristicFile", AllArr.AdditionalSymmetry, os);
  SetThompsonSampling(eFull, "DualDescriptionThompsonFile",
                      AllArr.DualDescriptionProgram, os);
  SetHeuristic(eFull, "MethodInitialFacetSetFile", AllArr.InitialFacetSet, os);
  SetHeuristic(eFull, "BankSaveHeuristicFile", AllArr.BankSave, os);
  SetHeuristic(eFull, "CheckDatabaseBankFile", AllArr.CheckDatabaseBank, os);
  SetHeuristic(eFull, "ChosenDatabaseFile", AllArr.ChosenDatabase, os);
  SetHeuristic(eFull, "OrbitSplitTechniqueFile", AllArr.OrbitSplitTechnique, os);
  SetHeuristic(eFull, "CommThreadHeuristicFile", AllArr.CommThread, os);
  SetHeuristic(eFull, "ChoiceCanonicalizationFile",
               AllArr.ChoiceCanonicalization, os);
  //
  bool DD_Saving = BlockMETHOD.ListBoolValues.at("Saving");
  AllArr.DD_Saving = DD_Saving;
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
  bool SimpleExchangeScheme =
      BlockDATA.ListBoolValues.at("SimpleExchangeScheme");
  AllArr.SimpleExchangeScheme = SimpleExchangeScheme;
  //
  AllArr.dimEXT = EXTred.cols();
  //
#ifdef DEBUG_RECURSIVE_DUAL_DESC
  PrintPolyHeuristicSerial(AllArr, os);
#endif
  return AllArr;
}

template <typename T, typename Tgroup, typename Tidx_value>
void MainFunctionSerialDualDesc(FullNamelist const &eFull) {
  // Setting up the Control C event.
  ExitEvent = false;
  if (Get_InterceptCtrlC_statuc(eFull, std::cerr)) {
    std::cerr << "Before submission of signal_callback_handler\n";
    signal(SIGINT, signal_callback_handler);
  } else {
    std::cerr << "Do not capture the CtrlC event\n";
  }
  //
  using Tint = typename Tgroup::Tint;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tkey = MyMatrix<T>;
  using Tval = TripleStore<Tgroup>;
  using Text_int = typename SubsetRankOneSolver<T>::Tint;
  MyMatrix<T> EXT = Get_EXT_DualDesc<T, Tidx>(eFull, std::cerr);
  Tgroup GRP = Get_GRP_DualDesc<Tgroup>(eFull, std::cerr);
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  MyMatrix<Text_int> EXTred_int = Get_EXT_int(EXTred);
  PolyHeuristicSerial<Tint> AllArr =
      Read_AllStandardHeuristicSerial<T, Tint>(eFull, EXTred, std::cerr);
  //
  std::map<std::string, Tint> TheMap =
    ComputeInitialMap<Tint,T,Tgroup>(EXTred, GRP, AllArr);
  auto get_vectface = [&]() -> vectface {
    if (AllArr.bank_parallelization_method == "serial") {
      using Tbank = DataBank<Tkey, Tval>;
      Tbank TheBank(AllArr.BANK_Saving, AllArr.BANK_Prefix, std::cerr);
      return DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value>(
          TheBank, EXTred, EXTred_int, GRP, TheMap, AllArr, AllArr.DD_Prefix,
          std::cerr);
    }
    if (AllArr.bank_parallelization_method == "bank_asio") {
      using Tbank = DataBankAsioClient<Tkey, Tval>;
      Tbank TheBank(AllArr.port);
      return DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value>(
          TheBank, EXTred, EXTred_int, GRP, TheMap, AllArr, AllArr.DD_Prefix,
          std::cerr);
    }
    std::cerr
        << "Failed to find a matching entry for bank_parallelization_method\n";
    std::cerr << "Allowed methods are serial, bank_asio\n";
    throw TerminalException{1};
  };
  vectface TheOutput = get_vectface();
  std::cerr << "|TheOutput|=" << TheOutput.size() << "\n";
  //
  OutputFacets(EXT, GRP, TheOutput, AllArr.OUTfile, AllArr.OutFormat,
               std::cerr);
}

template <typename T, typename Tgroup>
vectface DualDescriptionStandard(const MyMatrix<T> &EXT, const Tgroup &GRP, PolyHeuristicSerial<typename Tgroup::Tint> & AllArr, std::ostream& os) {
  using TintGroup = typename Tgroup::Tint;
  using Tkey = MyMatrix<T>;
  using Tval = TripleStore<Tgroup>;
  using Tidx_value = int32_t;
  using Text_int = typename SubsetRankOneSolver<T>::Tint;
  bool BANK_Saving = false;
  std::string BANK_Prefix = "totally_irrelevant_first";
  //
#ifdef DEBUG_RECURSIVE_DUAL_DESC
  PrintPolyHeuristicSerial(AllArr, os);
#endif
  //
  std::string DD_Prefix = "totally_irrelevant_second";
  //
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  MyMatrix<Text_int> EXTred_int = Get_EXT_int(EXTred);
  using Tbank = DataBank<Tkey, Tval>;
  Tbank TheBank(BANK_Saving, BANK_Prefix, os);
  std::map<std::string, TintGroup> TheMap =
    ComputeInitialMap<TintGroup,T,Tgroup>(EXTred, GRP, AllArr);
  return DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value>(
      TheBank, EXTred, EXTred_int, GRP, TheMap, AllArr, DD_Prefix, os);
}

// clang-format off
#endif  // SRC_DUALDESC_POLY_RECURSIVEDUALDESC_H_
// clang-format on

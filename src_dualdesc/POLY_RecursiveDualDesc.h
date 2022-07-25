// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_RECURSIVEDUALDESC_H_
#define SRC_DUALDESC_POLY_RECURSIVEDUALDESC_H_

#include "GRP_DoubleCoset.h"
#include "MAT_MatrixInt.h"
#include "Namelist.h"
#include "POLY_DirectDualDesc.h"
#include "POLY_Heuristics.h"
#include "POLY_SamplingFacet.h"
#include "Temp_PolytopeEquiStab.h"
#include "Timings.h"
#include "POLY_Kskeletton.h"

#include "POLY_GAP.h"
// #include "POLY_netcdf_file.h"
#include "Databank.h"
#include "MatrixGroupBasic.h"
#include "basic_datafile.h"
#include <limits>
#include <map>
#include <signal.h>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// #define MURMUR_HASH
// #define ROBIN_HOOD_HASH
#define SUBSET_HASH

// #define UNORDERED_MAP
#define TSL_SPARSE_MAP
// #define TSL_ROBIN_MAP
// #define TSL_HOPSCOTCH_MAP

#ifdef UNORDERED_MAP
#include <unordered_map>
#include <unordered_set>
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
  MyMatrix<T> EXTretB = RemoveFractionMatrix(EXTret);
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

template <typename Tbank, typename T, typename Tgroup, typename Tidx_value>
void insert_entry_in_bank(Tbank &bank, MyMatrix<T> const &EXT,
                          WeightMatrix<true, T, Tidx_value> const &WMat,
                          Tgroup const &TheGRPrelevant,
                          bool const &BankSymmCheck, vectface const &ListFace) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  if (!BankSymmCheck) {
    // The computation was already done for the full symmetry group. Only
    // canonic form is needed.
    std::pair<MyMatrix<T>, std::vector<Tidx>> ePair =
        CanonicalizationPolytopePair<T, Tidx, Tidx_value>(EXT, WMat);
    vectface ListFaceO(EXT.rows());
    Telt perm1 = Telt(ePair.second);
    Telt ePerm = ~perm1;
    Face eFaceImg(EXT.rows());
    for (auto &eFace : ListFace) {
      OnFace_inplace(eFaceImg, eFace, ePerm);
      ListFaceO.push_back(eFaceImg);
    }
    Tgroup GrpConj = TheGRPrelevant.GroupConjugate(ePerm);
    bank.InsertEntry(std::move(ePair.first),
                     {std::move(GrpConj), std::move(ListFaceO)});
  } else {
    TripleCanonic<T, Tgroup> eTriple =
        CanonicalizationPolytopeTriple<T, Tgroup>(EXT, WMat);
    bool NeedRemapOrbit = eTriple.GRP.size() == TheGRPrelevant.size();
    vectface ListFaceO(EXT.rows());
    Telt perm1 = Telt(eTriple.ListIdx);
    Telt ePerm = ~perm1;
    if (!NeedRemapOrbit) {
      // We needed to compute the full group, but it turned out to be the same
      // as the input group.
      Face eFaceImg(EXT.rows());
      for (auto &eFace : ListFace) {
        OnFace_inplace(eFaceImg, eFace, ePerm);
        ListFaceO.push_back(eFaceImg);
      }
    } else {
      // The full group is bigger than the input group. So we need to reduce.
      UNORD_SET<Face> SetFace;
      Face eFaceImg(EXT.rows());
      for (auto &eFace : ListFace) {
        OnFace_inplace(eFaceImg, eFace, ePerm);
        Face eIncCan = eTriple.GRP.CanonicalImage(eFaceImg);
        SetFace.insert(eIncCan);
      }
      for (auto &eInc : SetFace) {
        ListFaceO.push_back(eInc);
      }
    }
    bank.InsertEntry(std::move(eTriple.EXT),
                     {std::move(eTriple.GRP), std::move(ListFaceO)});
  }
}

template <typename Tbank, typename T, typename Tgroup, typename Tidx_value>
vectface getdualdesc_in_bank(Tbank &bank, MyMatrix<T> const &EXT,
                             WeightMatrix<true, T, Tidx_value> const &WMat,
                             Tgroup const &GRP, std::ostream & os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::pair<MyMatrix<T>, std::vector<Tidx>> ePair =
      CanonicalizationPolytopePair<T, Tidx, Tidx_value>(EXT, WMat);
  const PairStore<Tgroup> &RecAns = bank.GetDualDesc(ePair.first);
  if (RecAns.ListFace.size() == 0) {
    return vectface(0);
  }
  os << "Finding a matching entry\n";
  vectface ListReprTrans(EXT.rows());
  Telt ePerm = Telt(ePair.second);
  Face eFaceImg(EXT.rows());
  for (auto const &eFace : RecAns.ListFace) {
    OnFace_inplace(eFaceImg, eFace, ePerm);
    ListReprTrans.push_back(eFaceImg);
  }
  if (GRP.size() == RecAns.GRP.size())
    return ListReprTrans;
  Tgroup GrpConj = RecAns.GRP.GroupConjugate(ePerm);
  return OrbitSplittingListOrbit(GrpConj, GRP, ListReprTrans, os);
}

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

template <typename Tidx, typename Tint>
std::vector<Tint> GetAllPossibilities(std::map<Tidx, int> const &eMap) {
  std::vector<Tint> LVal = {1};
  for (auto &kv : eMap) {
    std::vector<Tint> NewVal;
    Tint ePow = 1;
    Tint kv_first_tint =
        size_t(kv.first); // Needed conversion because of Mac typing issues.
    for (int i = 0; i <= kv.second; i++) {
      for (auto &eVal : LVal)
        NewVal.push_back(ePow * eVal);
      ePow *= kv_first_tint;
    }
    LVal = NewVal;
  }
  return LVal;
}

template <typename T, typename Tgroup> struct DataFacetCan {
  size_t SelectedOrbit;
  Face eInc;
  FlippingFramework<T> FF;
  const Tgroup &GRP;
  Tgroup Stab;
  Face flip(const Face &f) const {
    Face eFlip = FF.Flip(f);
    return GRP.CanonicalImage(eFlip);
  }
};

template <typename T, typename Tgroup> struct DataFacetRepr {
  size_t SelectedOrbit;
  Face eInc;
  FlippingFramework<T> FF;
  const Tgroup &GRP;
  Tgroup Stab;
  Face flip(const Face &f) const { return FF.Flip(f); }
};

template <typename Tint> struct UndoneOrbitInfo {
  size_t nbOrbitDone;
  Tint nbUndone;
  Face eSetUndone;
};

template <typename Tint>
UndoneOrbitInfo<Tint> get_default_undoneinfo(int n_rows) {
  Face f(n_rows);
  return {0, 0, f};
}

template <typename Tint>
UndoneOrbitInfo<Tint>
CombineUndoneOrbitInfo(const std::vector<UndoneOrbitInfo<Tint>> &LComb) {
  size_t nbOrbitDone = LComb[0].nbOrbitDone;
  Tint nbUndone = LComb[0].nbUndone;
  Face f = LComb[0].eSetUndone;
  for (size_t i = 1; i < LComb.size(); i++) {
    nbOrbitDone += LComb[i].nbOrbitDone;
    nbUndone += LComb[i].nbUndone;
    f &= LComb[i].eSetUndone;
  }
  return {nbOrbitDone, nbUndone, f};
}

template <typename Tint>
bool ComputeStatusUndone(const UndoneOrbitInfo<Tint> &eComb,
                         const Tint &CritSiz) {
  if (eComb.nbOrbitDone > 0)
    if (eComb.nbUndone <= CritSiz || eComb.eSetUndone.count() > 0)
      return true;
  return false;
}

// The condition on nbOrbitDone make the check more complex.
// For parallel, we use this monotonic partial check as heuristic
// about whether to do the major checks or not.
template <typename Tint>
bool MonotonicCheckStatusUndone(const UndoneOrbitInfo<Tint> &eComb,
                                const Tint &CritSiz) {
  if (eComb.nbUndone <= CritSiz || eComb.eSetUndone.count() > 0)
    return true;
  return false;
}

/*
  This is the advanced termination criterion.
  
 */
template <typename T, typename Tgroup, typename Teval_recur>
bool EvaluationConnectednessCriterion(const MyMatrix<T> &FAC, const Tgroup &GRP,
                                      const vectface &vf, Teval_recur f_recur,
                                      std::ostream & os) {
  using Tint = typename Tgroup::Tint;
  size_t n_rows = FAC.rows();
  size_t n_cols = FAC.cols();
  size_t n_vert = vf.size();
  MyMatrix<T> EXT(n_vert, n_cols);
  os << "  EvaluationConnectednessCriterion : n_cols=" << n_cols << " n_rows=" << n_rows << " n_vert=" << n_vert << "\n";
  for (size_t i_vert = 0; i_vert < n_vert; i_vert++) {
    Face f = vf[i_vert];
    MyVector<T> eEXT = FindFacetInequality(FAC, f);
    for (size_t i_col = 0; i_col < n_cols; i_col++)
      EXT(i_vert, i_col) = eEXT(i_col);
  }
  os << "  We have EXT\n";
  auto rank_vertset = [&](const std::vector<size_t> &elist) -> size_t {
    auto f = [&](MyMatrix<T> &M, size_t eRank, size_t iRow) -> void {
      size_t pos = elist[iRow];
      for (size_t i_col = 0; i_col < n_cols; i_col++)
        M(eRank, i_col) = EXT(pos, i_col);
    };
    SelectionRowCol<T> eSelect =
        TMat_SelectRowCol_Kernel<T>(elist.size(), n_cols, f);
    return eSelect.TheRank;
  };
  using pfr = std::pair<size_t, Face>;
  auto evaluate_single_entry = [&](const pfr &x) -> bool {
    os << "  evaluate_single_entry pfr.first=" << x.first << " |pfr.second|=" << x.second.size() << " / " << x.second.count() << "\n";
    std::vector<size_t> f_v;
    for (size_t i = 0; i < n_rows; i++)
      if (x.second[i] == 1)
        f_v.push_back(i);
    auto is_vert_in_face = [&](const Face &g) -> bool {
      for (auto &idx : f_v)
        if (g[idx] == 0)
          return false;
      return true;
    };
    std::vector<size_t> list_vert;
    Face fint(n_rows);
    for (size_t i = 0; i < n_rows; i++)
      fint[i] = 1;
    for (size_t i_vert = 0; i_vert < n_vert; i_vert++) {
      Face e_vert = vf[i_vert];
      if (is_vert_in_face(e_vert)) {
        list_vert.push_back(i_vert);
        fint &= e_vert;
      }
    }
    if (true) {
      os << "  x=(" << x.first << ",[";
      bool IsFirst=true;
      for (auto &eVal : f_v) {
        if (!IsFirst)
          os << ",";
        IsFirst=false;
        os << eVal;
      }
      os << "]) |list_vert|=" << list_vert.size()
         << " |fint|=" << fint.count() << " |f|=" << f_v.size() << "\n";
    }
    if (fint.count() > f_v.size()) {
      os << "  Exit 1: linear programming check\n";
      return true; // This is the linear programming check
    }
    size_t n_cols_rel = n_cols - x.first;
    if (list_vert.size() <= n_cols_rel - 2) {
      os << "  Exit 2: pure Balinski case\n";
      return true; // This is the pure Balinski case
    }
    if (rank_vertset(list_vert) <= n_cols_rel - 2) {
      os << "  |list_vert|=" << list_vert.size()
         << " |rank_vertset|=" << rank_vertset(list_vert)
         << " n_cols_rel=" << n_cols_rel
         << "\n";
      os << "  Exit 3: rank computation, a little subtler Balinski computation\n";
      return true; // This is the rank computation. A little advanced Balinski,
                   // see the paper
    }
    os << "  Exit 4: nothing works, exiting\n";
    return false;
  };
  std::unordered_map<Face, bool> map_face_status;
  auto get_opt_face_status = [&](const pfr &x) -> std::optional<bool> {
    Face f_can = GRP.CanonicalImage(x.second);
    auto iter = map_face_status.find(f_can);
    if (iter == map_face_status.end()) {
      return {};
    } else {
      return iter->second;
    }
  };
  auto insert_pfr = [&](const pfr &x, const bool &val) -> bool {
    Face f_can = GRP.CanonicalImage(x.second);
    map_face_status[f_can] = val;
    return val;
  };
  std::function<bool(const pfr &)> get_face_status = [&](const pfr &x) -> bool {
    std::optional<bool> val_opt = get_opt_face_status(x);
    if (val_opt) {
      return *val_opt;
    }
    bool val = evaluate_single_entry(x);
    if (val) {
      return insert_pfr(x, val);
    } else {
      if (!f_recur(x))
        return insert_pfr(x, false);
      // Looking at the facets and maybe we can so conclude
      Tgroup eStab = GRP.Stabilizer_OnSets(x.second);
      vectface vf = SPAN_face_LinearProgramming(x.second, eStab, FAC, GRP);
      auto get_value = [&]() -> bool {
        Tint siz_false = 0;
        for (auto &eFace : vf) {
          bool val_f = get_face_status({x.first + 1, eFace});
          if (!val_f) {
            Tgroup eStab_B = eStab.Stabilizer_OnSets(eFace);
            Tint orb_size = eStab.size() / eStab_B.size();
            siz_false += orb_size;
            // If we cannot prove connectivity for just 1 facet, then
            // connectivity holds.
            if (siz_false > 1)
              return false;
          }
        }
        return true;
      };
      return insert_pfr(x, get_value());
    }
  };
  pfr init_pfr{0, Face(n_rows)};
  return get_face_status(init_pfr);
}

template <typename Tint, typename Torbsize, typename Tidx>
struct FaceOrbsizeContainer {
public:
  struct _SingEnt {
    Face face;
    Torbsize idx_orb;
  };
  using SingEnt = _SingEnt;
  /* TRICK 2: We keep the list of orbit and the map. We could in principle have
     built the map from the start since we know the occurring orders. However,
     since some orbitsize never occur
     this would have populated it with entries that never occur and so slow it
     down. */
  UNORD_MAP<Tint, Torbsize> OrbSize_Map;
  std::vector<Tint>
      ListPossOrbsize; // Canonically computed from the list of factors
  /* TRICK 3: Knowing the factorization of the order of the group allow us to
     know exactly
     what are the possible orbitsize occurring and so the number of bits needed
     to encode them */
  size_t n_act;
  size_t n_bit_orbsize;
  size_t delta;
  Tint TotalNumber;
  size_t nbOrbitDone;
  Tint nbUndone;
  size_t nbOrbit;
  std::vector<uint8_t> Vappend;
  std::vector<uint8_t>
      ListOrbit; // This CANNOT be replaced by vectface as we hack our way and
                 // so making a vectface will not allow
  // conversion functions that depend only on n_act and n_bit_orbsize.
  SingEnt FaceToSingEnt(Face const &f_in) const {
    SingEnt se{Face(n_act), 0};
    for (size_t i = 0; i < n_act; i++)
      se.face[i] = f_in[i];
    size_t i_acc = n_act;
    Torbsize pow = 1;
    for (size_t i = 0; i < n_bit_orbsize; i++) {
      se.idx_orb += Torbsize(f_in[i_acc]) * pow;
      i_acc++;
      pow *= 2;
    }
    return se;
  }
  Face SingEntToFace(const Face &face, const size_t &idx_orb) const {
    Face f(delta);
    for (size_t i = 0; i < n_act; i++)
      f[i] = face[i];
    size_t work_idx = idx_orb;
    size_t i_acc = n_act;
    for (size_t i = 0; i < n_bit_orbsize; i++) {
      bool val = work_idx % 2;
      f[i_acc] = val;
      i_acc++;
      work_idx = work_idx / 2;
    }
    return f;
  }
  // Database code that uses ListOrbit;
  SingEnt RetrieveListOrbitEntry(size_t const &i_orb) const {
    SingEnt se{Face(n_act), 0};
    size_t i_acc = delta * i_orb;
    for (size_t i = 0; i < n_act; i++) {
      se.face[i] = getbit(ListOrbit, i_acc);
      i_acc++;
    }
    Torbsize pow = 1;
    for (size_t i = 0; i < n_bit_orbsize; i++) {
      se.idx_orb += Torbsize(getbit(ListOrbit, i_acc)) * pow;
      i_acc++;
      pow *= 2;
    }
    return se;
  }
  Face RetrieveListOrbitFace(size_t const &i_orb) const {
    Face face(n_act);
    size_t i_acc = delta * i_orb;
    for (size_t i = 0; i < n_act; i++) {
      face[i] = getbit(ListOrbit, i_acc);
      i_acc++;
    }
    return face;
  }
  void InsertListOrbitEntry(SingEnt const &eEnt) {
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
      bool val = eEnt.face[i];
      setbit(ListOrbit, i_acc, val);
      i_acc++;
    }
    size_t work_idx = eEnt.idx_orb;
    for (size_t i = 0; i < n_bit_orbsize; i++) {
      bool val = work_idx % 2;
      setbit(ListOrbit, i_acc, val);
      i_acc++;
      work_idx = work_idx / 2;
    }
  }
  void InsertListOrbitFace(Face const &face) {
    // Insert bytes to avoid a memory segfault.
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
      setbit(ListOrbit, i_acc, val);
      i_acc++;
    }
  }
  void InsertListOrbitIdxOrb(Torbsize const &idx_orb) {
    /* TRICK 8: The computation of the stabilizer is needed for getting the
       orbitsize but this is expensive to do. Therefore we first insert the list
       of faces and if found to be new then we insert afterwards the idx_orb */
    size_t i_acc = nbOrbit * delta + n_act;
    Torbsize work_idx = idx_orb;
    for (size_t i = 0; i < n_bit_orbsize; i++) {
      bool val = work_idx % 2;
      setbit(ListOrbit, i_acc, val);
      i_acc++;
      work_idx = work_idx / 2;
    }
  }
  // Group functionalities.
  Torbsize GetOrbSizeIndex(Tint const &orbSize) {
    /* TRICK 4: value 0 is the default constructed one and so using it we can
       find if the entry is new or not in only one call */
    Torbsize &idx = OrbSize_Map[orbSize];
    if (idx == 0) { // A rare case. The linear loop should be totally ok
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
  void Counts_InsertOrbit(const bool &status, const size_t &idx_orb) {
    Tint orbSize = ListPossOrbsize[idx_orb];
    TotalNumber += orbSize;
    if (status) {
      nbOrbitDone++;
    } else {
      nbUndone += orbSize;
    }
    nbOrbit++;
  }
  void Counts_SetOrbitDone(const size_t &idx_orb) {
    nbUndone -= ListPossOrbsize[idx_orb];
    nbOrbitDone++;
  }
  FaceOrbsizeContainer(const std::map<Tidx, int> &LFact, const size_t &n_act)
      : n_act(n_act) {
    TotalNumber = 0;
    nbOrbitDone = 0;
    nbUndone = 0;
    nbOrbit = 0;
    size_t n_factor = 1;
    for (auto &kv : LFact) {
      n_factor *= (1 + kv.second);
    }
    /* TRICK 4: We need to add 1 because of shift by 1 in the OrbSize_Map */
    n_bit_orbsize = get_matching_power(n_factor + 1);
    delta = n_bit_orbsize + n_act;
    ListPossOrbsize = GetAllPossibilities<Tidx, Tint>(LFact);
    Vappend = std::vector<uint8_t>((delta + 7) / 8, 0);
  }
};

template <typename T, typename Tgroup>
vectface DirectComputationInitialFacetSet_Group(const MyMatrix<T> &EXT,
                                                const Tgroup &GRP,
                                                const std::string &ansSamp) {
  // We can do a little better by passing a lambda to the
  // DirectComputationInitialFacetSet but that is a little overkill right now
  size_t nbRow = EXT.rows();
  vectface list_face(nbRow);
  for (auto &eFace : DirectComputationInitialFacetSet(EXT, ansSamp))
    list_face.push_back(GRP.CanonicalImage(eFace));
  return list_face;
}

template <typename T_inp, typename Tint_inp, typename Tgroup_inp>
struct DatabaseCanonic {
public:
  using T = T_inp;
  using Tint = Tint_inp;
  using Tgroup = Tgroup_inp;
  using Telt = typename Tgroup::Telt;
  const MyMatrix<T> &EXT;
  const Tgroup &GRP;
  using DataFacet = DataFacetCan<T, Tgroup>;
  using Torbsize = uint16_t;
  using Tidx = typename Telt::Tidx;
  using SingEnt = typename FaceOrbsizeContainer<Tint, Torbsize, Tidx>::SingEnt;
  int nbRow;
  int nbCol;
  size_t delta;
  FaceOrbsizeContainer<Tint, Torbsize, Tidx> foc;

private:
  Tint groupOrder;
  UNORD_SET<size_t, std::function<size_t(size_t)>,
            std::function<bool(size_t, size_t)>>
      DictOrbit;
  std::map<size_t, std::vector<size_t>> CompleteList_SetUndone;
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

  void InsertEntryDatabase(Face const &face, bool const &status,
                           size_t const &idx_orb, size_t const &pos) {
    if (!status) {
      size_t len = face.count();
      CompleteList_SetUndone[len].push_back(pos);
    }
    foc.Counts_InsertOrbit(status, idx_orb);
  }
  DatabaseCanonic(MyMatrix<T> const &_EXT, Tgroup const &_GRP)
      : EXT(_EXT), GRP(_GRP), foc(GRP.factor_size(), EXT.rows()) {
    groupOrder = GRP.size();

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
    size_t n_ent_bit = 8 * sizeof(size_t); // The selection
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
    std::function<size_t(size_t)> fctHash = [&](size_t idx) -> size_t {
      size_t pos = delta * idx;
#if defined MURMUR_HASH || defined ROBIN_HOOD_HASH
      for (size_t i = 0; i < n_act; i++) {
        bool val = getbit(foc.ListOrbit, pos);
        setbit(V_hash, i, val);
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
        bool val = getbit(foc.ListOrbit, idx);
        setbit_ptr(ptr2, i, val);
      }
      return hash;
#endif
    };
    std::function<bool(size_t, size_t)> fctEqual = [&](size_t idx1,
                                                       size_t idx2) -> bool {
      size_t pos1 = delta * idx1;
      size_t pos2 = delta * idx2;
      for (size_t i = 1; i < n_act;
           i++) { // TRICK 9: Two faces will differ by at least 2 bits
        bool val1 = getbit(foc.ListOrbit, pos1);
        bool val2 = getbit(foc.ListOrbit, pos2);
        if (val1 != val2)
          return false;
        pos1++;
        pos2++;
      }
      return true;
    };
    DictOrbit =
        UNORD_SET<size_t, std::function<size_t(size_t)>,
                  std::function<bool(size_t, size_t)>>({}, fctHash, fctEqual);
  }
  ~DatabaseCanonic() {}
  vectface FuncListOrbitIncidence() {
    DictOrbit.clear();
    CompleteList_SetUndone.clear();
    vectface retListOrbit(n_act);
    for (size_t i_orbit = 0; i_orbit < foc.nbOrbit; i_orbit++) {
      Face f = foc.RetrieveListOrbitFace(i_orbit);
      retListOrbit.push_back(f);
    }
    return retListOrbit;
  }
  std::optional<Face> FuncInsert(
      Face const
          &face_can) { // The face should have been canonicalized beforehand.
    foc.InsertListOrbitFace(face_can);
    DictOrbit.insert(foc.nbOrbit);
    if (DictOrbit.size() == foc.nbOrbit) // Insertion did not raise the count
                                         // and so it was already present
      return {};
    /* TRICK 8: The insertion yield something new. So now we compute the
     * expensive stabilizer */
    Tint ordStab = GRP.Stabilizer_OnSets(face_can).size();
    Tint orbSize = groupOrder / ordStab;
    Torbsize idx_orb = foc.GetOrbSizeIndex(orbSize);
    foc.InsertListOrbitIdxOrb(idx_orb);
    InsertEntryDatabase(face_can, false, idx_orb, foc.nbOrbit);
    return foc.SingEntToFace(face_can, idx_orb);
  }
  vectface ComputeInitialSet(const std::string &ansSamp) {
    return DirectComputationInitialFacetSet_Group(EXT, GRP, ansSamp);
  }
  void FuncPutOrbitAsDone(size_t const &i_orb) {
    SingEnt eEnt = foc.RetrieveListOrbitEntry(i_orb);
    size_t len = eEnt.face.count();
    /* TRICK 1: We copy the last element in first position to erase it and then
     * pop_back the vector. */
    std::vector<size_t> &V = CompleteList_SetUndone[len];
    if (V.size() == 1) {
      CompleteList_SetUndone.erase(len);
    } else {
      V[0] = V[V.size() - 1];
      V.pop_back();
    }
    foc.Counts_SetOrbitDone(eEnt.idx_orb);
  }
  DataFacetCan<T, Tgroup> FuncGetMinimalUndoneOrbit() {
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
        return {pos, f, FlippingFramework<T>(EXT, f), GRP,
                ReducedGroupAction(Stab, f)};
      }
    }
    std::cerr << "Failed to find an undone orbit\n";
    throw TerminalException{1};
  }
  void InsertListOrbitEntry(SingEnt const &eEnt, const size_t &i_orbit) {
    foc.InsertListOrbitEntry(eEnt);
    DictOrbit.insert(i_orbit);
  }

private:
  struct IteratorType {
  private:
    const FaceOrbsizeContainer<Tint, Torbsize, Tidx> &foc;
    std::map<size_t, std::vector<size_t>>::const_iterator iter;
    size_t pos;

  public:
    IteratorType(const FaceOrbsizeContainer<Tint, Torbsize, Tidx> &foc,
                 std::map<size_t, std::vector<size_t>>::const_iterator iter,
                 size_t pos)
        : foc(foc), iter(iter), pos(pos) {}
    Face operator*() const {
      return foc.RetrieveListOrbitFace(iter->second[pos]);
    }
    IteratorType &operator++() {
      pos++;
      if (pos == iter->second.size()) {
        iter++;
        pos = 0;
      }
      return *this;
    }
    IteratorType operator++(int) {
      IteratorType tmp = *this;
      pos++;
      if (pos == iter->second.size()) {
        iter++;
        pos = 0;
      }
      return tmp;
    }
    bool operator!=(const IteratorType &x) const {
      return pos != x.pos || iter != x.iter;
    }
    bool operator==(const IteratorType &x) const {
      return pos == x.pos && iter == x.iter;
    }
  };

public:
  using iterator = IteratorType;
  iterator begin_undone() const {
    return IteratorType(foc, CompleteList_SetUndone.begin(), 0);
  }
  iterator end_undone() const {
    return IteratorType(foc, CompleteList_SetUndone.end(), 0);
  }
};

template <typename T_inp, typename Tint_inp, typename Tgroup_inp,
          typename Frepr, typename Fstab, typename Finv>
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
  using SingEnt = typename FaceOrbsizeContainer<Tint, Torbsize, Tidx>::SingEnt;
  int nbRow;
  int nbCol;
  size_t delta;
  FaceOrbsizeContainer<Tint, Torbsize, Tidx> foc;

private:
  Tint groupOrder;
  std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>>
      CompleteList_SetUndone;
  std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>> CompleteList_SetDone;
  size_t n_act;
  Frepr f_repr;
  Fstab f_stab;
  Finv f_inv;

public:
  DatabaseRepr() = delete;
  DatabaseRepr(const DatabaseRepr<T, Tint, Tgroup, Frepr, Fstab, Finv> &) =
      delete;
  DatabaseRepr(DatabaseRepr<T, Tint, Tgroup, Frepr, Fstab, Finv> &&) = delete;
  DatabaseRepr &
  operator=(const DatabaseRepr<T, Tint, Tgroup, Frepr, Fstab, Finv> &) = delete;

  void InsertEntryDatabase(Face const &face, bool const &status,
                           size_t const &idx_orb, size_t const &pos) {
    size_t len = face.count();
    size_t eInv = f_inv(face);
    if (status) {
      CompleteList_SetDone[len][eInv].push_back(pos);
    } else {
      CompleteList_SetUndone[len][eInv].push_back(pos);
    }
    foc.Counts_InsertOrbit(status, idx_orb);
  }
  DatabaseRepr(MyMatrix<T> const &_EXT, Tgroup const &_GRP, Frepr f_repr,
               Fstab f_stab, Finv f_inv)
      : EXT(_EXT), GRP(_GRP), foc(GRP.factor_size(), EXT.rows()),
        f_repr(f_repr), f_stab(f_stab), f_inv(f_inv) {
    groupOrder = GRP.size();

    /* TRICK 6: The UNORD_SET only the index and this saves in memory usage. */
    n_act = GRP.n_act();
    delta = foc.delta;
    nbRow = EXT.rows();
    nbCol = EXT.cols();
  }
  ~DatabaseRepr() {}
  vectface FuncListOrbitIncidence() {
    CompleteList_SetUndone.clear();
    CompleteList_SetDone.clear();
    vectface retListOrbit(n_act);
    for (size_t i_orbit = 0; i_orbit < foc.nbOrbit; i_orbit++) {
      Face f = foc.RetrieveListOrbitFace(i_orbit);
      retListOrbit.push_back(f);
    }
    return retListOrbit;
  }
  std::optional<Face> FuncInsert(Face const &face_i) {
    size_t len = face_i.count();
    size_t eInv = f_inv(face_i);
    if (CompleteList_SetDone.count(len) == 1) {
      if (CompleteList_SetDone[len].count(eInv) == 1) {
        for (size_t &i_orb : CompleteList_SetDone[len][eInv]) {
          Face face_e = foc.RetrieveListOrbitFace(i_orb);
          bool test = f_repr(face_i, face_e);
          if (test)
            return {};
        }
      }
    }
    if (CompleteList_SetUndone.count(len) == 1) {
      if (CompleteList_SetUndone[len].count(eInv) == 1) {
        for (size_t &i_orb : CompleteList_SetUndone[len][eInv]) {
          Face face_e = foc.RetrieveListOrbitFace(i_orb);
          bool test = f_repr(face_i, face_e);
          if (test)
            return {};
        }
      }
    }
    foc.InsertListOrbitFace(face_i);
    Tint ordStab = f_stab(face_i).size();
    Tint orbSize = groupOrder / ordStab;
    Torbsize idx_orb = foc.GetOrbSizeIndex(orbSize);
    foc.InsertListOrbitIdxOrb(idx_orb);
    InsertEntryDatabase(face_i, false, idx_orb, foc.nbOrbit);
    //
    return foc.SingEntToFace(face_i, idx_orb);
  }
  vectface ComputeInitialSet(const std::string &ansSamp) {
    return DirectComputationInitialFacetSet(EXT, ansSamp);
  }
  void FuncPutOrbitAsDone(size_t const &iOrb) {
    SingEnt eEnt = foc.RetrieveListOrbitEntry(iOrb);
    size_t len = eEnt.face.count();
    /* TRICK 1: We copy the last element in first position to erase it and then
     * pop_back the vector. */
    size_t eInv = f_inv(eEnt.face);
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
    foc.Counts_SetOrbitDone(eEnt.idx_orb);
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
  void InsertListOrbitEntry(SingEnt const &eEnt,
                            [[maybe_unused]] const size_t &i_orbit) {
    foc.InsertListOrbitEntry(eEnt);
  }

private:
  struct IteratorType {
  private:
    const FaceOrbsizeContainer<Tint, Torbsize, Tidx> &foc;
    std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>>::const_iterator
        iter1;
    std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>>::const_iterator
        iter1_end;
    UNORD_MAP<size_t, std::vector<size_t>>::const_iterator iter2;
    size_t pos;

  public:
    IteratorType(
        const FaceOrbsizeContainer<Tint, Torbsize, Tidx> &foc,
        std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>>::const_iterator
            iter1,
        std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>>::const_iterator
            iter1_end,
        UNORD_MAP<size_t, std::vector<size_t>>::const_iterator iter2,
        size_t pos)
        : foc(foc), iter1(iter1), iter1_end(iter1_end), iter2(iter2), pos(pos) {
    }
    Face operator*() const {
      return foc.RetrieveListOrbitFace(iter2->second[pos]);
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
    IteratorType &operator++() {
      PtrIncrease();
      return *this;
    }
    IteratorType operator++(int) {
      IteratorType tmp = *this;
      PtrIncrease();
      return tmp;
    }
    bool operator!=(const IteratorType &x) const {
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

public:
  using iterator = IteratorType;
  iterator begin_undone() const {
    std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>>::const_iterator
        iter1 = CompleteList_SetUndone.begin();
    std::map<size_t, UNORD_MAP<size_t, std::vector<size_t>>>::const_iterator
        iter1_end = CompleteList_SetUndone.end();
    if (iter1 == iter1_end)
      return IteratorType(foc, iter1, iter1_end, {}, 0);
    UNORD_MAP<size_t, std::vector<size_t>>::const_iterator iter2 =
        CompleteList_SetUndone.at(iter1->first).begin();
    return IteratorType(foc, iter1, iter1_end, iter2, 0);
  }
  iterator end_undone() const {
    return IteratorType(foc, CompleteList_SetUndone.end(),
                        CompleteList_SetUndone.end(), {}, 0);
  }
};

template <typename TbasicBank> struct DatabaseOrbits {
public:
  using Tgroup = typename TbasicBank::Tgroup;
  using T = typename TbasicBank::T;
  using Telt = typename Tgroup::Telt;
  using Tint = typename TbasicBank::Tint;
  using SingEnt = typename TbasicBank::SingEnt;
  Tint CritSiz;
  TbasicBank &bb;

private:
  std::string MainPrefix;
  std::string eFileEXT, eFileGRP, eFileNB, eFileFB, eFileFF;
  /* TRICK 7: Using separate files for faces and status allow us to gain
     locality. The faces are written one by one while the access to status is
     random */
  // This is for storing the number of orbits of the polytope
  FileNumber *fn;
  // This is for storing the status of the orbits
  FileBool *fb;
  // This is for storing the faces and the index of orbit
  FileFace *ff;
  bool is_opened;
  bool SavingTrigger;
  bool AdvancedTerminationCriterion;
  std::ostream &os;
  size_t delta;
  std::string strPresChar;

public:
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
                 const bool &_SavingTrigger, const bool& _AdvancedTerminationCriterion, std::ostream &os)
      : CritSiz(bb.EXT.cols() - 2), bb(bb), SavingTrigger(_SavingTrigger),
        AdvancedTerminationCriterion(_AdvancedTerminationCriterion),
        os(os) {
    os << "MainPrefix=" << MainPrefix << "\n";
    eFileEXT = MainPrefix + ".ext";
    eFileGRP = MainPrefix + ".grp";
    eFileNB = MainPrefix + ".nb";
    eFileFB = MainPrefix + ".fb";
    eFileFF = MainPrefix + ".ff";
    strPresChar = "|EXT|=" + std::to_string(bb.nbRow) + "/" +
                  std::to_string(bb.nbCol) +
                  " |GRP|=" + std::to_string(bb.GRP.size());
    delta = bb.delta;
    fn = nullptr;
    fb = nullptr;
    ff = nullptr;
    is_opened = false;
    if (SavingTrigger) {
      size_t n_orbit;
      if (IsExistingFile(eFileEXT)) {
        os << "Opening existing files (NB, FB, FF)\n";
        fn = new FileNumber(eFileNB, false);
        n_orbit = fn->getval();
        fb = new FileBool(eFileFB, n_orbit);
        ff = new FileFace(eFileFF, bb.delta, n_orbit);
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
        // Writing polytope
        std::ofstream os_ext(eFileEXT);
        WriteMatrix(os_ext, bb.EXT);
        // Opening the files
        fn = new FileNumber(eFileNB, true);
        fb = new FileBool(eFileFB);
        ff = new FileFace(eFileFF, bb.delta);
        n_orbit = 0;
        fn->setval(n_orbit);
      }
      is_opened = true;
      //
      os << "Inserting orbits, n_orbit=" << n_orbit << "\n";
      for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
        Face f = ff->getface(i_orbit);
        SingEnt eEnt = bb.foc.FaceToSingEnt(f);
        bool status = fb->getbit(i_orbit);
        // The DictOrbit
        bb.InsertListOrbitEntry(eEnt, i_orbit);
        // The other fields
        bb.InsertEntryDatabase(eEnt.face, status, eEnt.idx_orb, i_orbit);
      }
      print_status();
    }
  }
  ~DatabaseOrbits() {
    /* TRICK 5: The destructor does NOT destroy the database! This is because it
       can be used in another call. Note that the returning of the list of orbit
       does destroy the database and this gives a small window in which bad
       stuff can happen.
     */
    if (is_opened) {
      delete fb;
      delete fn;
      delete ff;
    }
    os << "Clean closing of the DatabaseOrbits\n";
  }
  vectface FuncListOrbitIncidence() {
    if (SavingTrigger) {
      delete fb;
      delete fn;
      delete ff;
      is_opened = false;
      //
      RemoveFile(eFileNB);
      RemoveFile(eFileFB);
      RemoveFile(eFileFF);
      RemoveFile(eFileEXT);
      RemoveFile(eFileGRP);
    }
    return bb.FuncListOrbitIncidence();
  }
  void FuncInsert(Face const &face_can) {
    std::optional<Face> test = bb.FuncInsert(face_can);
    if (test && SavingTrigger) {
      fb->setbit(bb.foc.nbOrbit - 1, false);
      ff->setface(bb.foc.nbOrbit - 1, *test);
      fn->setval(bb.foc.nbOrbit);
    }
  }
  vectface ComputeInitialSet(const std::string &ansSamp) {
    return bb.ComputeInitialSet(ansSamp);
  }
  void FuncPutOrbitAsDone(size_t const &i_orb) {
    bb.FuncPutOrbitAsDone(i_orb);
    if (SavingTrigger) {
      fb->setbit(i_orb, true);
    }
    print_status();
  }
  Face ComputeIntersectionUndone() const {
    size_t n_row = bb.EXT.rows();
    Face eSetReturn(n_row);
    for (size_t i_row = 0; i_row < n_row; i_row++)
      eSetReturn[i_row] = 1;
    typename TbasicBank::iterator iter = bb.begin_undone();
    while (iter != bb.end_undone()) {
      eSetReturn &= OrbitIntersection(bb.GRP, *iter);
      if (eSetReturn.count() == 0)
        return eSetReturn;
      iter++;
    }
    return eSetReturn;
  }
  size_t FuncNumberOrbit() const { return bb.foc.nbOrbit; }
  typename TbasicBank::DataFacet FuncGetMinimalUndoneOrbit() {
    typename TbasicBank::DataFacet data = bb.FuncGetMinimalUndoneOrbit();
    os << strPresChar << " Considering orbit " << data.SelectedOrbit
       << " |inc|=" << data.eInc.count()
       << " |stab|=" << data.Stab.size() << "\n";
    return data;
  }
  bool attempt_connectedness_scheme() const {
    vectface vf(bb.nbRow);
    std::vector<Telt> LGen = bb.GRP.GeneratorsOfGroup();
    // We need an heuristic to avoid building too large orbits.
    // A better system would have to balance out the cost of
    // doing that check with respect to the dual description itsef.
    size_t max_siz = 1000;
    if (bb.foc.nbUndone > max_siz)
      return false;
    // Now explicit building of the set of vertices
    typename TbasicBank::iterator iterator = bb.begin_undone();
    while (iterator != bb.end_undone()) {
      vectface vfo = OrbitFace(*iterator, LGen);
      for (auto &uFace : vfo)
        vf.push_back(uFace);
      iterator++;
    }
    size_t max_iter = 100;
    size_t n_iter = 0;
    auto f_recur = [&](const std::pair<size_t, Face> &pfr) -> bool {
      n_iter++;
      os << "  f_recur n_iter=" << n_iter << "\n";
      if (n_iter == max_iter)
        return false;
      if (pfr.first > 1)
        return false;
      return true;
    };
    return EvaluationConnectednessCriterion(bb.EXT, bb.GRP, vf, f_recur, os);
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
      return attempt_connectedness_scheme();
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
std::map<std::string, Tint> ComputeInitialMap(const MyMatrix<T> &EXT,
                                              const Tgroup &GRP) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  std::map<std::string, Tint> TheMap;
  int delta = nbRow - nbCol;
  TheMap["groupsize"] = GRP.size();
  TheMap["incidence"] = nbRow;
  TheMap["rank"] = nbCol;
  TheMap["delta"] = delta;
  return TheMap;
}


template <typename Tint, typename T, typename Tgroup>
vectface DirectFacetOrbitComputation_ts(MyMatrix<T> const &EXT, Tgroup const &GRP,
                                        std::map<std::string, Tint> const& TheMap,
                                        ThompsonSamplingHeuristic<T> & ansprog_ts) {
  std::string ansprog = ansprog_ts.get_eval();
  vectface TheOutput = DirectFacetOrbitComputation(EXT, GRP, ansprog);
  ansprog_ts.pop();
  return TheOutput;
}



// Needs initial definition due to template crazyness
template <typename Tbank, typename T, typename Tgroup, typename Tidx_value>
vectface DUALDESC_AdjacencyDecomposition(
    Tbank &TheBank, MyMatrix<T> const &EXT, Tgroup const &GRP,
    PolyHeuristicSerial<typename Tgroup::Tint> const &AllArr,
    std::string const &ePrefix);

template <typename Tbank, typename T, typename Tgroup, typename Tidx_value,
          typename TbasicBank>
vectface Kernel_DUALDESC_AdjacencyDecomposition(
    Tbank &TheBank, TbasicBank &bb,
    PolyHeuristicSerial<typename Tgroup::Tint> const &AllArr,
    std::string const &ePrefix,
    std::map<std::string, typename Tgroup::Tint> const &TheMap, std::ostream& os) {
  using DataFacet = typename TbasicBank::DataFacet;
  DatabaseOrbits<TbasicBank> RPL(bb, ePrefix, AllArr.Saving, AllArr.AdvancedTerminationCriterion, os);
  if (RPL.FuncNumberOrbit() == 0) {
    std::string ansSamp = HeuristicEvaluation(TheMap, AllArr.InitialFacetSet);
    for (auto &face : RPL.ComputeInitialSet(ansSamp))
      RPL.FuncInsert(face);
  }
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
    vectface TheOutput =
        DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value>(
            TheBank, df.FF.EXT_face, df.Stab, AllArr, NewPrefix, os);
    for (auto &eOrbB : TheOutput) {
      Face eFlip = df.flip(eOrbB);
      RPL.FuncInsert(eFlip);
    }
    RPL.FuncPutOrbitAsDone(SelectedOrbit);
  };
  return RPL.FuncListOrbitIncidence();
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
    Tbank &TheBank, MyMatrix<T> const &EXT, Tgroup const &GRP,
    PolyHeuristicSerial<typename Tgroup::Tint> const &AllArr,
    std::string const &ePrefix, std::ostream & os) {
  using Tgr = GraphListAdj;
  using Tint = typename Tgroup::Tint;
  if (ExitEvent) {
    std::cerr << "Terminating the program by Ctrl-C\n";
    throw TerminalException{1};
  }
  if (AllArr.max_runtime > 0) {
    int runtime = si(AllArr.start);
    if (runtime > AllArr.max_runtime) {
      std::cerr << "The maximum runtime has been elapsed. max_runtime = " << AllArr.max_runtime << "\n";
      throw TerminalException{1};
    }
  }
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  LazyWMat<T, Tidx_value> lwm(EXT);
  //
  // Now computing the groups
  //
  std::map<std::string, Tint> TheMap = ComputeInitialMap<Tint>(EXT, GRP);
  //
  // Checking if the entry is present in the map.
  //
  std::string ansBankCheck =
      HeuristicEvaluation(TheMap, AllArr.CheckDatabaseBank);
  if (ansBankCheck == "yes") {
    vectface ListFace = getdualdesc_in_bank(TheBank, EXT, lwm.GetWMat(), GRP, os);
    if (ListFace.size() > 0)
      return ListFace;
  }
  //
  // Now computing the groups
  //
  std::string ansSplit = HeuristicEvaluation(TheMap, AllArr.Splitting);
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
  auto compute_split_or_not = [&]() -> vectface {
    if (ansSplit != "split") {
      TheGRPrelevant = GRP;
      std::string ansProg =
          HeuristicEvaluation(TheMap, AllArr.DualDescriptionProgram);
      BankSymmCheck = true;
      return DirectFacetOrbitComputation(EXT, GRP, ansProg);
    } else {
      std::string ansSymm =
          HeuristicEvaluation(TheMap, AllArr.AdditionalSymmetry);
      os << "ansSymm=" << ansSymm << "\n";
      if (ansSymm == "yes") {
        TheGRPrelevant = GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(
            lwm.GetWMat());
        NeedSplit = TheGRPrelevant.size() != GRP.size();
        BankSymmCheck = false;
      } else {
        TheGRPrelevant = GRP;
        BankSymmCheck = true;
      }
      Tint GroupSizeComp = TheGRPrelevant.size();
      os << "RESPAWN a new ADM computation |GRP|=" << GroupSizeComp
         << " TheDim=" << nbCol << " |EXT|=" << nbRow << "\n";
      //      std::string MainPrefix = ePrefix + "Database_" +
      //      std::to_string(nbRow) + "_" + std::to_string(nbCol);
      std::string MainPrefix = ePrefix + "D_" + std::to_string(nbRow);
      std::string ansChosenDatabase =
          HeuristicEvaluation(TheMap, AllArr.ChosenDatabase);
      os << "DUALDESC_ChosenDatabase : ChosenDatabase = "
         << ansChosenDatabase << "\n";
      if (ansChosenDatabase == "canonic") {
        using TbasicBank = DatabaseCanonic<T, Tint, Tgroup>;
        TbasicBank bb(EXT, TheGRPrelevant);
        return Kernel_DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup,
                                                      Tidx_value, TbasicBank>(
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
        auto f_stab = [&](const Face &f) -> Tgroup {
          return TheGRPrelevant.Stabilizer_OnSets(f);
        };
        auto f_inv = [&](const Face &f) -> size_t {
          return GetLocalInvariantWeightMatrix(WMat, f);
        };
        using TbasicBank = DatabaseRepr<T, Tint, Tgroup, decltype(f_repr),
                                        decltype(f_stab), decltype(f_inv)>;
        TbasicBank bb(EXT, TheGRPrelevant, f_repr, f_stab, f_inv);
        return Kernel_DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup,
                                                      Tidx_value, TbasicBank>(
            TheBank, bb, AllArr, MainPrefix, TheMap, os);
      }
      std::cerr << "Failed to find a matching entry\n";
      throw TerminalException{1};
    }
  };
  vectface ListOrbitFaces = compute_split_or_not();
  SingletonTime end;
  TheMap["time"] = s(start, end);
  std::string ansBank = HeuristicEvaluation(TheMap, AllArr.BankSave);
  os << "elapsed_seconds=" << s(start, end) << " ansBank=" << ansBank
            << " NeedSplit=" << NeedSplit << "\n";
  if (ansBank == "yes") {
    insert_entry_in_bank(TheBank, EXT, lwm.GetWMat(), TheGRPrelevant,
                         BankSymmCheck, ListOrbitFaces);
  }
  if (NeedSplit) {
    return OrbitSplittingListOrbit(TheGRPrelevant, GRP, ListOrbitFaces, os);
  } else {
    return ListOrbitFaces;
  }
}

std::vector<size_t> get_subset_index_rev(const size_t &n_act) {
  size_t n_ent_bit = 8 * sizeof(size_t); // The size of the selection
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
      size_t pos = size_t(round(frac * static_cast<double>(i)));
      if (pos < 0)
        pos = 0;
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
  uint8_t *ptr2 = (uint8_t *)ptr1;
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
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, int> ListIntValues1;
  ListStringValues1["EXTfile"] = "unset.ext";
  ListStringValues1["GRPfile"] = "unset.grp";
  ListStringValues1["OUTfile"] = "unset.out";
  ListStringValues1["OutFormat"] = "GAP";
  ListBoolValues1["DeterministicRuntime"] = true;
  ListStringValues1["parallelization_method"] = "serial";
  ListIntValues1["port"] = 1234;
  ListIntValues1["max_runtime"] = -1;
  ListBoolValues1["AdvancedTerminationCriterion"] = false;
  SingleBlock BlockDATA;
  BlockDATA.ListStringValues = ListStringValues1;
  BlockDATA.ListBoolValues = ListBoolValues1;
  BlockDATA.ListIntValues = ListIntValues1;
  ListBlock["DATA"] = BlockDATA;
  // HEURISTIC
  std::map<std::string, std::string> ListStringValuesH;
  ListStringValuesH["SplittingHeuristicFile"] = "unset.heu";
  ListStringValuesH["AdditionalSymmetryHeuristicFile"] = "unset.heu";
  ListStringValuesH["DualDescriptionHeuristicFile"] = "unset.heu";
  ListStringValuesH["MethodInitialFacetSetFile"] = "unset.heu";
  ListStringValuesH["BankSaveHeuristicFile"] = "unset.heu";
  ListStringValuesH["CheckDatabaseBankFile"] = "unset.heu";
  ListStringValuesH["ChosenDatabaseFile"] = "unset.heu";
  SingleBlock BlockHEURIS;
  BlockHEURIS.ListStringValues = ListStringValuesH;
  ListBlock["HEURISTIC"] = BlockHEURIS;
  // METHOD
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  ListBoolValues2["Saving"] = false;
  ListStringValues2["Prefix"] = "/irrelevant/";
  SingleBlock BlockMETHOD;
  BlockMETHOD.ListIntValues = ListIntValues2;
  BlockMETHOD.ListBoolValues = ListBoolValues2;
  BlockMETHOD.ListDoubleValues = ListDoubleValues2;
  BlockMETHOD.ListStringValues = ListStringValues2;
  ListBlock["METHOD"] = BlockMETHOD;
  // BANK
  std::map<std::string, int> ListIntValues3;
  std::map<std::string, bool> ListBoolValues3;
  std::map<std::string, double> ListDoubleValues3;
  std::map<std::string, std::string> ListStringValues3;
  ListStringValues3["Prefix"] = "./unset/";
  ListBoolValues3["Saving"] = false;
  SingleBlock BlockBANK;
  BlockBANK.ListBoolValues = ListBoolValues3;
  BlockBANK.ListStringValues = ListStringValues3;
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

void OutputFacets(const vectface &TheOutput, const std::string &OUTfile,
                  const std::string &OutFormat) {
  if (OutFormat == "Magma") {
    std::ofstream os(OUTfile);
    os << "return ";
    VectVectInt_Magma_Print(os, TheOutput);
    os << ";\n";
    return;
  }
  if (OutFormat == "GAP") {
    std::ofstream os(OUTfile);
    os << "return ";
    VectVectInt_Gap_Print(os, TheOutput);
    os << ";\n";
    return;
  }
  if (OutFormat == "SetInt") {
    std::ofstream os(OUTfile);
    os << TheOutput.size() << "\n";
    for (const Face &face : TheOutput) {
      mpz_class res = getsetasint<mpz_class>(face);
      os << res << "\n";
    }
    return;
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



template<typename T, typename Tidx>
MyMatrix<T> Get_EXT_DualDesc(FullNamelist const &eFull) {
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  std::string EXTfile = BlockDATA.ListStringValues.at("EXTfile");
  IsExistingFileDie(EXTfile);
  std::cerr << "EXTfile=" << EXTfile << "\n";
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

template<typename Tgroup>
Tgroup Get_GRP_DualDesc(FullNamelist const &eFull) {
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  std::string GRPfile = BlockDATA.ListStringValues.at("GRPfile");
  IsExistingFileDie(GRPfile);
  std::cerr << "GRPfile=" << GRPfile << "\n";
  std::ifstream GRPfs(GRPfile);
  Tgroup GRP = ReadGroup<Tgroup>(GRPfs);
  return GRP;
}



template<typename Tint>
PolyHeuristicSerial<Tint> Read_AllStandardHeuristicSerial(FullNamelist const &eFull) {
  PolyHeuristicSerial<Tint> AllArr = AllStandardHeuristicSerial<Tint>();
  std::cerr << "We have AllArr\n";
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
  std::cerr << "OUTfile=" << OUTfile << "\n";
  //
  bool DeterministicRuntime =
      BlockDATA.ListBoolValues.at("DeterministicRuntime");
  std::cerr << "DeterministicRuntime=" << DeterministicRuntime << "\n";
  if (!DeterministicRuntime)
    srand_random_set();
  //
  std::string OutFormat = BlockDATA.ListStringValues.at("OutFormat");
  AllArr.OutFormat = OutFormat;
  std::cerr << "OutFormat=" << OutFormat << "\n";
  //
  int port_i = BlockDATA.ListIntValues.at("port");
  std::cerr << "port_i=" << port_i << "\n";
  short unsigned int port = port_i;
  AllArr.port = port;
  //
  std::string parallelization_method =
      BlockDATA.ListStringValues.at("parallelization_method");
  AllArr.parallelization_method = parallelization_method;
  std::cerr << "parallelization_method=" << parallelization_method << "\n";
  //
  SetHeuristic(eFull, "SplittingHeuristicFile", AllArr.Splitting);
  std::cerr << "SplittingHeuristicFile\n" << AllArr.Splitting << "\n";
  //
  SetHeuristic(eFull, "AdditionalSymmetryHeuristicFile",
               AllArr.AdditionalSymmetry);
  std::cerr << "AdditionalSymmetryHeuristicFile\n"
            << AllArr.AdditionalSymmetry << "\n";
  //
  SetHeuristic(eFull, "DualDescriptionHeuristicFile",
               AllArr.DualDescriptionProgram);
  std::cerr << "DualDescriptionHeuristicFile\n"
            << AllArr.DualDescriptionProgram << "\n";
  //
  SetHeuristic(eFull, "MethodInitialFacetSetFile", AllArr.InitialFacetSet);
  std::cerr << "MethodInitialFacetSetFile\n" << AllArr.InitialFacetSet << "\n";
  //
  SetHeuristic(eFull, "BankSaveHeuristicFile", AllArr.BankSave);
  std::cerr << "BankSaveHeuristicFile\n" << AllArr.BankSave << "\n";
  //
  SetHeuristic(eFull, "CheckDatabaseBankFile", AllArr.CheckDatabaseBank);
  std::cerr << "CheckDatabaseBank\n" << AllArr.CheckDatabaseBank << "\n";
  //
  SetHeuristic(eFull, "ChosenDatabaseFile", AllArr.ChosenDatabase);
  std::cerr << "ChosenDatabase\n" << AllArr.ChosenDatabase << "\n";
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
  bool AdvancedTerminationCriterion = BlockDATA.ListBoolValues.at("AdvancedTerminationCriterion");
  AllArr.AdvancedTerminationCriterion = AdvancedTerminationCriterion;
  //
  return AllArr;
}





template <typename T, typename Tgroup, typename Tidx_value>
void MainFunctionSerialDualDesc(FullNamelist const &eFull) {
  // Setting up the Control C event.
  ExitEvent = false;
  signal(SIGINT, signal_callback_handler);
  //
  using Tint = typename Tgroup::Tint;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tkey = MyMatrix<T>;
  using Tval = PairStore<Tgroup>;
  MyMatrix<T> EXT = Get_EXT_DualDesc<T,Tidx>(eFull);
  Tgroup GRP = Get_GRP_DualDesc<Tgroup>(eFull);
  PolyHeuristicSerial<Tint> AllArr = Read_AllStandardHeuristicSerial<Tint>(eFull);
  //
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  auto get_vectface = [&]() -> vectface {
    if (AllArr.parallelization_method == "serial") {
      using Tbank = DataBank<Tkey, Tval>;
      Tbank TheBank(AllArr.BANK_IsSaving, AllArr.BANK_Prefix);
      return DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value>(
          TheBank, EXTred, GRP, AllArr, AllArr.DD_Prefix, std::cerr);
    }
    if (AllArr.parallelization_method == "bank_asio") {
      using Tbank = DataBankClient<Tkey, Tval>;
      Tbank TheBank(AllArr.port);
      return DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value>(
          TheBank, EXTred, GRP, AllArr, AllArr.DD_Prefix, std::cerr);
    }
    std::cerr << "Failed to find a matching entry for parallelization_method\n";
    throw TerminalException{1};
  };
  vectface TheOutput = get_vectface();
  std::cerr << "|TheOutput|=" << TheOutput.size() << "\n";
  //
  OutputFacets(TheOutput, AllArr.OUTfile, AllArr.OutFormat);
}

template <typename T, typename Tgroup>
vectface DualDescriptionStandard(const MyMatrix<T> &EXT, const Tgroup &GRP) {
  using Tint = typename Tgroup::Tint;
  using Tkey = MyMatrix<T>;
  using Tval = PairStore<Tgroup>;
  using Tidx_value = int32_t;
  bool BANK_IsSaving = false;
  std::string BANK_Prefix = "totally_irrelevant_first";
  //
  PolyHeuristicSerial<Tint> AllArr = AllStandardHeuristicSerial<Tint>();
  std::cerr << "SplittingHeuristicFile\n" << AllArr.Splitting << "\n";
  std::cerr << "AdditionalSymmetryHeuristicFile\n"
            << AllArr.AdditionalSymmetry << "\n";
  std::cerr << "DualDescriptionHeuristicFile\n"
            << AllArr.DualDescriptionProgram << "\n";
  std::cerr << "MethodInitialFacetSetFile\n" << AllArr.InitialFacetSet << "\n";
  std::cerr << "BankSaveHeuristicFile\n" << AllArr.BankSave << "\n";
  std::cerr << "CheckDatabaseBank\n" << AllArr.CheckDatabaseBank << "\n";
  std::cerr << "ChosenDatabase\n" << AllArr.ChosenDatabase << "\n";
  //
  bool DD_Saving = false;
  std::string DD_Prefix = "totally_irrelevant_second";
  AllArr.Saving = DD_Saving;
  //
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  using Tbank = DataBank<Tkey, Tval>;
  Tbank TheBank(BANK_IsSaving, BANK_Prefix);
  return DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value>(
      TheBank, EXTred, GRP, AllArr, DD_Prefix, std::cerr);
}

// clang-format off
#endif  // SRC_DUALDESC_POLY_RECURSIVEDUALDESC_H_
// clang-format on

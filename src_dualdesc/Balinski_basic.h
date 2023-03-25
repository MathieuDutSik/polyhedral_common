// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_BALINSKI_BASIC_H_
#define SRC_DUALDESC_BALINSKI_BASIC_H_

// clang-format off
#include "Boost_bitset_kernel.h"
#include "MAT_Matrix.h"
#include <vector>
#include <unordered_map>
#include <utility>
// clang-format on

template <typename Tint> struct UndoneOrbitInfo {
  size_t nbOrbitDone;
  Tint nbUndone;
  Face eSetUndone;
};

template <typename Tint>
std::ostream &operator<<(std::ostream &os, UndoneOrbitInfo<Tint> const &obj) {
  os << "(nbOrbitDone=" << obj.nbOrbitDone;
  os << ", nbUndone=" << obj.nbUndone;
  os << ", eSetUndone=" << obj.eSetUndone.count() << ","
     << obj.eSetUndone.size() << ")";
  return os;
}

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

template <typename Tint> struct StatusUndoneOrbitInfo {
  bool status;
  UndoneOrbitInfo<Tint> erec;
};

namespace boost::serialization {

template <class Archive, typename Tint>
inline void serialize(Archive &ar, StatusUndoneOrbitInfo<Tint> &mesg,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("status", mesg.status);
  ar &make_nvp("nborbitdone", mesg.erec.nbOrbitDone);
  ar &make_nvp("nbundone", mesg.erec.nbUndone);
  ar &make_nvp("setundone", mesg.erec.eSetUndone);
}

template <class Archive, typename Tint>
inline void serialize(Archive &ar, UndoneOrbitInfo<Tint> &erec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("nborbitdone", erec.nbOrbitDone);
  ar &make_nvp("nbundone", erec.nbUndone);
  ar &make_nvp("setundone", erec.eSetUndone);
}
// clang-format off
}  // namespace boost::serialization
// clang-format on

template <typename TbasicBank> vectface ComputeSetUndone(TbasicBank const &bb) {
  vectface vf_undone(bb.nbRow);
  using Tgroup = typename TbasicBank::Tgroup;
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> LGen = bb.GRP.GeneratorsOfGroup();
  typename TbasicBank::iterator_face iter = bb.begin_face_undone();
  while (iter != bb.end_face_undone()) {
    vectface vf_orbit = OrbitFace(*iter, LGen);
    vf_undone.append(vf_orbit);
    iter++;
  }
  return vf_undone;
}

template <typename T>
MyMatrix<T> GetVertexSet_from_vectface(MyMatrix<T> const &FAC,
                                       vectface const &vf) {
  size_t n_cols = FAC.cols();
  size_t n_vert = vf.size();
  MyMatrix<T> EXT(n_vert, n_cols);
  for (size_t i_vert = 0; i_vert < n_vert; i_vert++) {
    Face f = vf[i_vert];
    MyVector<T> eEXT = FindFacetInequality(FAC, f);
    for (size_t i_col = 0; i_col < n_cols; i_col++)
      EXT(i_vert, i_col) = eEXT(i_col);
  }
  return EXT;
}

/*
  This is the advanced termination criterion.
  This combines a number of approaches:
  ---Balinski theorem itself.
  ---The linear programming check
  ---The rank check by computing the rank of the set of missing vertices.
  If too low then the same argument as Balinski theorem applies

  If those fails, then we compute the faces and if we can prove the
  connectedness for all the facets except at most 1 then the connectedness
  holds.

  The technique can be applied recursively. This cause of course a potential
  way to explode the runtime. Thus in order to avoid that we need a function
  f_recur that does the check and returns false when this goes too deep.

  We also use the canonical storing of the computed results in order to save
  runtime.

  References:
  ---Balinski, M. L. (1961), "On the graph structure of convex polyhedra in
  n-space", Pacific Journal of Mathematics, 11 (2): 431–434,
  ---Michel Deza, Mathieu Dutour Sikirić, Enumeration of the facets of cut
  polytopes over some highly symmetric graphs, preprint at arxiv:1501.05407,
  International Transactions in Operational Research 23-5 (2016) 853--860

  The EXT_undone is precomputed because it can be done in parallel.
 */
template <typename T, typename Tgroup, typename Teval_recur>
bool EvaluationConnectednessCriterion_Kernel(
    const MyMatrix<T> &FAC, const Tgroup &GRP, const MyMatrix<T> &EXT_undone,
    const vectface &vf_undone, Teval_recur f_recur, std::ostream &os) {
  using Tint = typename Tgroup::Tint;
  size_t n_rows = FAC.rows();
  size_t n_cols = FAC.cols();
  size_t n_vert = vf_undone.size();
  os << "  We have EXT\n";
  auto rank_vertset = [&](const std::vector<size_t> &elist) -> size_t {
    // Computation of the rank of the set without computing the full set.
    auto f = [&](MyMatrix<T> &M, size_t eRank, size_t iRow) -> void {
      size_t pos = elist[iRow];
      for (size_t i_col = 0; i_col < n_cols; i_col++)
        M(eRank, i_col) = EXT_undone(pos, i_col);
    };
    SelectionRowCol<T> eSelect =
        TMat_SelectRowCol_Kernel<T>(elist.size(), n_cols, f);
    return eSelect.TheRank;
  };
  using pfr = std::pair<size_t, Face>;
  auto evaluate_single_entry = [&](const pfr &x) -> bool {
    // We test the connectedness using the known criterions:
    // ---Balinski criterion
    // ---Linear programming check
    // ---Balinski with rank check
    os << "  evaluate_single_entry pfr.first=" << x.first
       << " |pfr.second|=" << x.second.size() << " / " << x.second.count()
       << "\n";
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
      Face e_vert = vf_undone[i_vert];
      if (is_vert_in_face(e_vert)) {
        list_vert.push_back(i_vert);
        fint &= e_vert;
      }
    }
    if (true) {
      os << "  x=(" << x.first << ",[";
      bool IsFirst = true;
      for (auto &eVal : f_v) {
        if (!IsFirst)
          os << ",";
        IsFirst = false;
        os << eVal;
      }
      os << "]) |list_vert|=" << list_vert.size() << " |fint|=" << fint.count()
         << " |f|=" << f_v.size() << "\n";
    }
    if (fint.count() > f_v.size()) {
      os << "  Exit 1: linear programming check\n";
      // This is the linear programming check, See DD 2016.
      return true;
    }
    size_t n_cols_rel = n_cols - x.first;
    if (list_vert.size() <= n_cols_rel - 2) {
      os << "  Exit 2: pure Balinski case\n";
      // This is the pure Balinski case
      return true;
    }
    if (rank_vertset(list_vert) <= n_cols_rel - 2) {
      os << "  |list_vert|=" << list_vert.size()
         << " |rank_vertset|=" << rank_vertset(list_vert)
         << " n_cols_rel=" << n_cols_rel << "\n";
      os << "  Exit 3: rank computation, a little subtler Balinski "
            "computation\n";
      // This is the rank computation. A little advanced Balinski,
      // see the Balinski paper and adapt the proof.
      return true;
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
      vectface vf_span = SPAN_face_LinearProgramming(x.second, eStab, FAC, GRP);
      auto get_value = [&]() -> bool {
        Tint siz_false = 0;
        for (auto &eFace : vf_span) {
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

template <typename T, typename Tgroup>
bool EvaluationConnectednessCriterion_PreKernel_field(const MyMatrix<T> &FAC,
                                                      const Tgroup &GRP,
                                                      const vectface &vf_undone,
                                                      std::ostream &os) {
  MyMatrix<T> EXT_undone = GetVertexSet_from_vectface(FAC, vf_undone);
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
  return EvaluationConnectednessCriterion_Kernel(FAC, GRP, EXT_undone,
                                                 vf_undone, f_recur, os);
}

template <typename T, typename Tgroup>
inline typename std::enable_if<is_ring_field<T>::value, bool>::type
EvaluationConnectednessCriterion_PreKernel(const MyMatrix<T> &FAC,
                                           const Tgroup &GRP,
                                           const vectface &vf_undone,
                                           std::ostream &os) {
  return EvaluationConnectednessCriterion_PreKernel(FAC, GRP, vf_undone, os);
}

template <typename T, typename Tgroup>
inline typename std::enable_if<!is_ring_field<T>::value, bool>::type
EvaluationConnectednessCriterion_PreKernel(const MyMatrix<T> &FAC,
                                           const Tgroup &GRP,
                                           const vectface &vf_undone,
                                           std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> FACfield = UniversalMatrixConversion<Tfield, T>(FAC);
  return EvaluationConnectednessCriterion_PreKernel(FACfield, GRP, vf_undone, os);
}

template <typename TbasicBank>
bool EvaluationConnectednessCriterion_Serial(TbasicBank const &bb,
                                             std::ostream &os) {
  using T = typename TbasicBank::T;
  using Tint = typename TbasicBank::Tint;
  // We need an heuristic to avoid building too large orbits.
  // A better system would have to balance out the cost of
  // doing that check with respect to the dual description itsef.
  Tint max_siz = 1000;
  //  os << "nbUndone=" << bb.foc.nbUndone << " nbOrbit=" << bb.foc.nbOrbit <<
  //  "\n"; os << "nbOrbitDone=" << bb.foc.nbOrbitDone << " TotalNumber=" <<
  //  bb.foc.TotalNumber << "\n";
  if (bb.foc.nbOrbitDone == 0 || bb.foc.nbUndone > max_siz)
    return false;
  // Now explicit building of the set of vertices
  vectface vf_undone = ComputeSetUndone(bb);
  //
  return EvaluationConnectednessCriterion_PreKernel(bb.EXT, bb.GRP,
                                                    vf_undone, os);
}

// clang-format off
#endif  // SRC_DUALDESC_BALINSKI_BASIC_H_
// clang-format on

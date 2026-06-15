// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_ROBUST_COVERING_TYPES_H_
#define SRC_ROBUST_COVERING_ROBUST_COVERING_TYPES_H_

// clang-format off
#include "boost_serialization.h"
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_ROBUST_COVERING_TYPES
#endif


// A family of vectors with the index which is the farthest.
template <typename Tint> struct GenericRobustM {
  int index;
  MyMatrix<Tint> M;
  MyVector<Tint> v_long() const { return GetMatrixRow(M, index); }
  std::vector<MyVector<Tint>> get_short_vectors() const {
    std::vector<MyVector<Tint>> l_v;
    int n_row = M.rows();
    for (int i_row=0; i_row<n_row; i_row++) {
      if (i_row != index) {
        MyVector<Tint> V = GetMatrixRow(M, i_row);
        l_v.emplace_back(std::move(V));
      }
    }
    return l_v;
  }
  template<typename T>
  MyMatrix<T> get_ext_t() const {
    int n_row = M.rows();
    int n_vect = n_row - 1;
    int dim = M.cols();
    int i_vect = 0;
    MyMatrix<T> EXTret(n_vect, dim + 1);
    for (int i_row=0; i_row<n_row; i_row++) {
      if (i_row != index) {
        EXTret(i_vect, 0) = 1;
        for (int i=0; i<dim; i++) {
          EXTret(i_vect, i+1) = UniversalScalarConversion<T,Tint>(M(i_row, i));
        }
        i_vect += 1;
      }
    }
    return EXTret;
  }
};

template <typename Tint>
std::ostream &operator<<(std::ostream &os, GenericRobustM<Tint> const &grm) {
  os << "GenericRobustM(index=" << grm.index << " M=\n";
  WriteMatrix(os, grm.M);
  os << ")";
  return os;
}

template <typename Tint>
void WriteEntryGAP(std::ostream &os_out, GenericRobustM<Tint> const &grm) {
  os_out << "rec(index:=" << grm.index << ", M:=";
  WriteMatrixGAP(os_out, grm.M);
  os_out << ")";
}

template <typename Tint>
void WriteEntryCPP(std::ostream &os, GenericRobustM<Tint> const &grm) {
  os << grm.index << "\n";
  WriteMatrix(os, grm.M);
}

template <typename Tint>
GenericRobustM<Tint> ReadEntryCPP_GenericRobustM(std::istream &is) {
  int index;
  is >> index;
  MyMatrix<Tint> M = ReadMatrix<Tint>(is);
  return {index, M};
}

namespace boost::serialization {

template <class Archive, typename Tint>
inline void serialize(Archive &ar, GenericRobustM<Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("index", val.index);
  ar &make_nvp("M", val.M);
}

} // namespace boost::serialization

// The convex block of the construction
template <typename T, typename Tint>
struct ConvexBlock {
  std::vector<GenericRobustM<Tint>> list_robust_m;
  SinglePolytope<T> sp;
};

template <typename T, typename Tint>
std::ostream &operator<<(std::ostream &os, ConvexBlock<T, Tint> const &cb) {
  os << "ConvexBlock(|list_robust_m|=" << cb.list_robust_m.size() << " FAC=\n";
  WriteMatrix(os, cb.sp.FAC);
  os << "EXT=\n";
  WriteMatrix(os, cb.sp.EXT);
  os << ")";
  return os;
}

template <typename T, typename Tint>
void WriteEntryGAP(std::ostream &os_out, ConvexBlock<T, Tint> const &cb) {
  os_out << "rec(list_robust_m:=[";
  bool IsFirst = true;
  for (auto &grm : cb.list_robust_m) {
    if (!IsFirst) {
      os_out << ",";
    }
    IsFirst = false;
    WriteEntryGAP(os_out, grm);
  }
  os_out << "], sp:=";
  WriteEntryGAP(os_out, cb.sp);
  os_out << ")";
}

template <typename T, typename Tint>
void WriteEntryCPP(std::ostream &os, ConvexBlock<T, Tint> const &cb) {
  size_t n = cb.list_robust_m.size();
  os << n << "\n";
  for (size_t i = 0; i < n; i++) {
    WriteEntryCPP(os, cb.list_robust_m[i]);
  }
  WriteEntryCPP(os, cb.sp);
}

template <typename T, typename Tint>
ConvexBlock<T, Tint> ReadEntryCPP_ConvexBlock(std::istream &is) {
  size_t n;
  is >> n;
  std::vector<GenericRobustM<Tint>> list_robust_m;
  for (size_t i = 0; i < n; i++) {
    list_robust_m.push_back(ReadEntryCPP_GenericRobustM<Tint>(is));
  }
  SinglePolytope<T> sp = ReadEntryCPP_SinglePolytope<T>(is);
  return {list_robust_m, sp};
}

namespace boost::serialization {

template <class Archive, typename T, typename Tint>
inline void serialize(Archive &ar, ConvexBlock<T, Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("list_robust_m", val.list_robust_m);
  ar &make_nvp("sp", val.sp);
}

} // namespace boost::serialization

template<typename T>
struct HardConvexBoundary {
  ConvexBoundary<T> sp;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, HardConvexBoundary<T> const &hcb) {
  os << "HardConvexBoundary(V=" << StringVectorGAP(hcb.sp.V) << ")";
  return os;
}

template <typename T>
void WriteEntryGAP(std::ostream &os_out, HardConvexBoundary<T> const &hcb) {
  os_out << "rec(sp:=";
  WriteEntryGAP(os_out, hcb.sp);
  os_out << ")";
}

template <typename T>
void WriteEntryCPP(std::ostream &os, HardConvexBoundary<T> const &hcb) {
  WriteEntryCPP(os, hcb.sp);
}

template <typename T>
HardConvexBoundary<T> ReadEntryCPP_HardConvexBoundary(std::istream &is) {
  ConvexBoundary<T> sp = ReadEntryCPP_ConvexBoundary<T>(is);
  return {sp};
}

namespace boost::serialization {

template <class Archive, typename T>
inline void serialize(Archive &ar, HardConvexBoundary<T> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("sp", val.sp);
}

} // namespace boost::serialization

template<typename T, typename Tint>
struct SoftConvexBoundary {
  ConvexBoundary<T> cb;
  std::vector<MyVector<Tint>> l_excluded_max; // The excluded vectors. Cannot use the ConvexBlock since they can vary.
  std::vector<GenericRobustM<Tint>> l_robust_m;
  MyVector<Tint> v_long() const {
#ifdef SANITY_CHECK_ROBUST_COVERING_TYPES
    if (l_robust_m.empty()) {
      std::cerr << "ROBUST: l_robust_m should be non-empty for getting the v_long\n";
      throw TerminalException{1};
    }
#endif
    MyVector<Tint> v_long = l_robust_m[0].v_long();
#ifdef SANITY_CHECK_ROBUST_COVERING_TYPES
    size_t len = l_robust_m.size();
    for (size_t i=1; i<len; i++) {
      MyVector<Tint> v = l_robust_m[i].v_long();
      if (v != v_long) {
        std::cerr << "ROBUST: Incoherent v_long in the structure\n";
        throw TerminalException{1};
      }
    }
#endif
    return v_long;
  }
};

template <typename T, typename Tint>
std::ostream &operator<<(std::ostream &os, SoftConvexBoundary<T, Tint> const &scb) {
  os << "SoftConvexBoundary(V=" << StringVectorGAP(scb.cb.V)
     << " |l_excluded_max|=" << scb.l_excluded_max.size()
     << " |l_robust_m|=" << scb.l_robust_m.size() << ")";
  return os;
}

template <typename T, typename Tint>
void WriteEntryGAP(std::ostream &os_out, SoftConvexBoundary<T, Tint> const &scb) {
  os_out << "rec(cb:=";
  WriteEntryGAP(os_out, scb.cb);
  os_out << ", l_excluded_max:=[";
  for (size_t i = 0; i < scb.l_excluded_max.size(); i++) {
    if (i > 0) {
      os_out << ",";
    }
    os_out << StringVectorGAP(scb.l_excluded_max[i]);
  }
  os_out << "], l_robust_m:=[";
  for (size_t i = 0; i < scb.l_robust_m.size(); i++) {
    if (i > 0) {
      os_out << ",";
    }
    WriteEntryGAP(os_out, scb.l_robust_m[i]);
  }
  os_out << "])";
}

template <typename T, typename Tint>
void WriteEntryCPP(std::ostream &os, SoftConvexBoundary<T, Tint> const &scb) {
  WriteEntryCPP(os, scb.cb);
  size_t n_excl = scb.l_excluded_max.size();
  os << n_excl << "\n";
  for (size_t i = 0; i < n_excl; i++) {
    WriteVector(os, scb.l_excluded_max[i]);
  }
  size_t n_robust = scb.l_robust_m.size();
  os << n_robust << "\n";
  for (size_t i = 0; i < n_robust; i++) {
    WriteEntryCPP(os, scb.l_robust_m[i]);
  }
}

template <typename T, typename Tint>
SoftConvexBoundary<T, Tint> ReadEntryCPP_SoftConvexBoundary(std::istream &is) {
  ConvexBoundary<T> cb = ReadEntryCPP_ConvexBoundary<T>(is);
  size_t n_excl;
  is >> n_excl;
  std::vector<MyVector<Tint>> l_excluded_max;
  for (size_t i = 0; i < n_excl; i++) {
    l_excluded_max.push_back(ReadVector<Tint>(is));
  }
  size_t n_robust;
  is >> n_robust;
  std::vector<GenericRobustM<Tint>> l_robust_m;
  for (size_t i = 0; i < n_robust; i++) {
    l_robust_m.push_back(ReadEntryCPP_GenericRobustM<Tint>(is));
  }
  return {cb, l_excluded_max, l_robust_m};
}

namespace boost::serialization {

template <class Archive, typename T, typename Tint>
inline void serialize(Archive &ar, SoftConvexBoundary<T, Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("cb", val.cb);
  ar &make_nvp("l_excluded_max", val.l_excluded_max);
  ar &make_nvp("l_robust_m", val.l_robust_m);
}

} // namespace boost::serialization

/*
  The robust_m_min is defining the P-polytope.
  This is what we are after in the end.
  ----
  It is a full enumeration result describing a PVoronoi
  * v_long is the vector realizing the maximum.
  * l_robust_m_min is the set of parallelepipeds realizing
    the minimum. Could be more than 1.
  * 
 */
template <typename T, typename Tint> struct PVoronoi {
  MyVector<Tint> v_long;
  std::vector<GenericRobustM<Tint>> l_robust_m_min;
  std::vector<ConvexBlock<T,Tint>> l_cb; // The list of convex blocks.
  std::vector<HardConvexBoundary<T>> l_hcb;
  GeneralizedPolytope<T> gp;
  MyMatrix<T> EXT; // The list of vertices as defined from the generalized polytope.
};

template <typename T, typename Tint>
std::ostream &operator<<(std::ostream &os, PVoronoi<T, Tint> const &pv) {
  os << "PVoronoi(\n  v_long=" << pv.v_long << "\n";
  os << "  |l_robust_m_min|=" << pv.l_robust_m_min.size() << "\n";
  for (size_t i = 0; i < pv.l_robust_m_min.size(); i++) {
    os << "  l_robust_m_min[" << i << "]=" << pv.l_robust_m_min[i] << "\n";
  }
  os << "  |l_cb|=" << pv.l_cb.size() << "\n";
  for (size_t i = 0; i < pv.l_cb.size(); i++) {
    os << "  l_cb[" << i << "]=" << pv.l_cb[i] << "\n";
  }
  os << "  |l_hcb|=" << pv.l_hcb.size() << "\n";
  for (size_t i = 0; i < pv.l_hcb.size(); i++) {
    os << "  l_hcb[" << i << "]=" << pv.l_hcb[i] << "\n";
  }
  os << "  EXT=\n";
  WriteMatrix(os, pv.EXT);
  os << ")";
  return os;
}

template <typename T, typename Tint>
void WriteEntryGAP(std::ostream &os_out, PVoronoi<T, Tint> const &pv) {
  os_out << "rec(v_long:=" << StringVectorGAP(pv.v_long) << ", l_robust_m_min:=[";
  for (size_t i = 0; i < pv.l_robust_m_min.size(); i++) {
    if (i > 0) {
      os_out << ",";
    }
    WriteEntryGAP(os_out, pv.l_robust_m_min[i]);
  }
  os_out << "], l_cb:=[";
  for (size_t i = 0; i < pv.l_cb.size(); i++) {
    if (i > 0) {
      os_out << ",";
    }
    WriteEntryGAP(os_out, pv.l_cb[i]);
  }
  os_out << "], l_hcb:=[";
  for (size_t i = 0; i < pv.l_hcb.size(); i++) {
    if (i > 0) {
      os_out << ",";
    }
    WriteEntryGAP(os_out, pv.l_hcb[i]);
  }
  os_out << "], gp:=";
  WriteEntryGAP(os_out, pv.gp);
  os_out << ", EXT:=";
  WriteMatrixGAP(os_out, pv.EXT);
  os_out << ")";
}

template <typename T, typename Tint>
void WriteEntryCPP(std::ostream &os, PVoronoi<T, Tint> const &pv) {
  WriteVector(os, pv.v_long);
  size_t n_robust = pv.l_robust_m_min.size();
  os << n_robust << "\n";
  for (size_t i = 0; i < n_robust; i++) {
    WriteEntryCPP(os, pv.l_robust_m_min[i]);
  }
  size_t n_cb = pv.l_cb.size();
  os << n_cb << "\n";
  for (size_t i = 0; i < n_cb; i++) {
    WriteEntryCPP(os, pv.l_cb[i]);
  }
  size_t n_hcb = pv.l_hcb.size();
  os << n_hcb << "\n";
  for (size_t i = 0; i < n_hcb; i++) {
    WriteEntryCPP(os, pv.l_hcb[i]);
  }
  WriteEntryCPP(os, pv.gp);
  WriteMatrix(os, pv.EXT);
}

template <typename T, typename Tint>
PVoronoi<T, Tint> ReadEntryCPP_PVoronoi(std::istream &is) {
  MyVector<Tint> v_long = ReadVector<Tint>(is);
  size_t n_robust;
  is >> n_robust;
  std::vector<GenericRobustM<Tint>> l_robust_m_min;
  for (size_t i = 0; i < n_robust; i++) {
    l_robust_m_min.push_back(ReadEntryCPP_GenericRobustM<Tint>(is));
  }
  size_t n_cb;
  is >> n_cb;
  std::vector<ConvexBlock<T, Tint>> l_cb;
  for (size_t i = 0; i < n_cb; i++) {
    l_cb.push_back(ReadEntryCPP_ConvexBlock<T, Tint>(is));
  }
  size_t n_hcb;
  is >> n_hcb;
  std::vector<HardConvexBoundary<T>> l_hcb;
  for (size_t i = 0; i < n_hcb; i++) {
    l_hcb.push_back(ReadEntryCPP_HardConvexBoundary<T>(is));
  }
  GeneralizedPolytope<T> gp = ReadEntryCPP_GeneralizedPolytope<T>(is);
  MyMatrix<T> EXT = ReadMatrix<T>(is);
  return {v_long, l_robust_m_min, l_cb, l_hcb, gp, EXT};
}

template <typename T, typename Tint>
PVoronoi<T, Tint> ReadEntryCPP_PVoronoi_File(std::string const& file_name) {
  if (!IsExistingFile(file_name)) {
    std::cerr << "Error in ReadMatrixFile\n";
    std::cerr << "file_name=" << file_name << " does not appear to exist\n";
    throw TerminalException{1};
  }
  std::ifstream is(file_name);
  return ReadEntryCPP_PVoronoi<T,Tint>(is);
}

namespace boost::serialization {

template <class Archive, typename T, typename Tint>
inline void serialize(Archive &ar, PVoronoi<T, Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("l_robust_m_min", val.l_robust_m_min);
  ar &make_nvp("l_cb", val.l_cb);
  ar &make_nvp("l_hcb", val.l_hcb);
  ar &make_nvp("gp", val.gp);
  ar &make_nvp("EXT", val.EXT);
}

} // namespace boost::serialization


/*
  The robust_m_min is defining the P-polytope.
  This is what we are after in the end.
  ----
  It is a partial enumeration result.
 */
template <typename T, typename Tint> struct PVoronoiPart {
  MyVector<Tint> v_long;
  std::vector<GenericRobustM<Tint>> l_robust_m_min;
  std::vector<ConvexBlock<T,Tint>> l_cb; // The list of convex blocks.
  std::vector<HardConvexBoundary<T>> l_hcb;
  std::map<MyVector<T>, std::vector<SoftConvexBoundary<T,Tint>>> map_scb;
  size_t n_scb() const {
    size_t n_val = 0;
    for (auto & kv: map_scb) {
      n_val += kv.second.size();
    }
    return n_val;
  }
  SoftConvexBoundary<T,Tint> get_one_scb() const {
    auto iter = map_scb.begin();
#ifdef SANITY_CHECK_ROBUST_COVERING_TYPES
    if (iter == map_scb.end()) {
      std::cerr << "ROBUST_TYPES: map_sb should not be empty\n";
      throw TerminalException{1};
    }
#endif
    std::vector<SoftConvexBoundary<T,Tint>> const& l_scb = iter->second;
#ifdef SANITY_CHECK_ROBUST_COVERING_TYPES
    if (l_scb.empty()) {
      std::cerr << "ROBUST_TYPES: l_scb should not be empty\n";
      throw TerminalException{1};
    }
#endif
    return l_scb[0];
  }
  void insert_scb(SoftConvexBoundary<T,Tint> const& scb, std::ostream& os) {
    // The insertion of this soft convex boundary can increase but also decrease it.
    // We need to keep track of what we already have.
    MyVector<T> V = - scb.cb.V;
    MyVector<T> V_opp = - V;
    if (!map_scb.contains(V_opp)) {
      map_scb[V].push_back(scb);
      // Not present, needs to insert
      return;
    }
    std::vector<SoftConvexBoundary<T,Tint>> l_scb = map_scb.at(V);
    std::vector<SoftConvexBoundary<T,Tint>> l_scb_new;
    std::vector<ConvexBoundary<T>> l_cb_ins{scb.cb};
    for (auto & scb_old: l_scb) {
      std::vector<MyVector<Tint>> const& l_excluded_max = scb_old.l_excluded_max;
      std::vector<GenericRobustM<Tint>> const& l_robust_m = scb_old.l_robust_m;
      std::vector<ConvexBoundary<T>> l_cb_ins_new;
      for (auto & cb_ins: l_cb_ins) {
        std::vector<ConvexBoundary<T>> l_cb_a = convex_boundary_minus_cb(scb_old.cb, cb_ins, os);
        for (auto & cb_a: l_cb_a) {
          SoftConvexBoundary<T,Tint> cbN{cb_a, l_excluded_max, l_robust_m};
          l_scb_new.push_back(cbN);
        }
        std::vector<ConvexBoundary<T>> l_cb_b = convex_boundary_minus_cb(cb_ins, scb_old.cb, os);
        for (auto& cb_b: l_cb_b) {
          l_cb_ins_new.push_back(cb_b);
        }
      }
      l_cb_ins = l_cb_ins_new;
    }
    if (l_cb_ins.size() > 0) {
      for (auto & cb_ins: l_cb_ins) {
        SoftConvexBoundary<T,Tint> scb_ins{cb_ins, scb.l_excluded_max, scb.l_robust_m};
        map_scb[V].push_back(scb_ins);
      }
    }
    if (l_scb_new.size() == 0) {
      map_scb.erase(V_opp);
    } else {
      map_scb[V_opp] = l_scb_new;
    }
  }
  std::vector<SoftConvexBoundary<T,Tint>> get_l_scb() const {
    std::vector<SoftConvexBoundary<T,Tint>> l_scb;
    for (auto & kv: map_scb) {
      for (auto & entry: kv.second) {
        l_scb.push_back(entry);
      }
    }
    return l_scb;
  }
};

template <typename T, typename Tint>
std::ostream &operator<<(std::ostream &os, PVoronoiPart<T, Tint> const &pvp) {
  os << "PVoronoiPart(\n  v_long=" << pvp.v_long << "\n";
  os << "  |l_robust_m_min|=" << pvp.l_robust_m_min.size() << "\n";
  for (size_t i = 0; i < pvp.l_robust_m_min.size(); i++) {
    os << "  l_robust_m_min[" << i << "]=" << pvp.l_robust_m_min[i] << "\n";
  }
  os << "  |l_cb|=" << pvp.l_cb.size() << "\n";
  for (size_t i = 0; i < pvp.l_cb.size(); i++) {
    os << "  l_cb[" << i << "]=" << pvp.l_cb[i] << "\n";
  }
  os << "  |l_hcb|=" << pvp.l_hcb.size() << "\n";
  for (size_t i = 0; i < pvp.l_hcb.size(); i++) {
    os << "  l_hcb[" << i << "]=" << pvp.l_hcb[i] << "\n";
  }
  os << "  |map_scb|=" << pvp.map_scb.size() << "\n";
  for (size_t i = 0; i < pvp.l_scb.size(); i++) {
    os << "  l_scb[" << i << "]=" << pvp.l_scb[i] << "\n";
  }
  os << ")";
  return os;
}

// clang-format off
#endif  // SRC_ROBUST_COVERING_ROBUST_COVERING_TYPES_H_
// clang-format on


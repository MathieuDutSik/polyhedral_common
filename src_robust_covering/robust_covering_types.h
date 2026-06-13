// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_ROBUST_COVERING_TYPES_H_
#define SRC_ROBUST_COVERING_ROBUST_COVERING_TYPES_H_

// clang-format off
#include "boost_serialization.h"
// clang-format on

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
  int index_cb; // The corresponding face;
  ConvexBoundary<T> sp;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, HardConvexBoundary<T> const &hcb) {
  os << "HardConvexBoundary(index_cb=" << hcb.index_cb << " V=" << StringVectorGAP(hcb.sp.V) << ")";
  return os;
}

template <typename T>
void WriteEntryGAP(std::ostream &os_out, HardConvexBoundary<T> const &hcb) {
  os_out << "rec(index_cb:=" << hcb.index_cb << ", sp:=";
  WriteEntryGAP(os_out, hcb.sp);
  os_out << ")";
}

namespace boost::serialization {

template <class Archive, typename T>
inline void serialize(Archive &ar, HardConvexBoundary<T> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("index_cb", val.index_cb);
  ar &make_nvp("sp", val.sp);
}

} // namespace boost::serialization

template<typename T, typename Tint>
struct SoftConvexBoundary {
  int index_cb; // The corresponding face;
  ConvexBoundary<T> cb;
  std::vector<MyVector<Tint>> l_excluded_max; // The excluded vectors. Cannot use the ConvexBlock since they can vary.
  std::vector<GenericRobustM<Tint>> l_robust_m;
  MyVector<Tint> v_long() const {
#ifdef SANITY_CHECK_ENUM_P_POLYTOPES
    if (l_robust_m.empty()) {
      std::cerr << "ROBUST: l_robust_m should be non-empty for getting the v_long\n";
      throw TerminalException{1};
    }
#endif
    MyVector<Tint> v_long = l_robust_m[0].v_long();
#ifdef SANITY_CHECK_ENUM_P_POLYTOPES
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
  os << "SoftConvexBoundary(index_cb=" << scb.index_cb
     << " V=" << StringVectorGAP(scb.cb.V)
     << " |l_excluded_max|=" << scb.l_excluded_max.size()
     << " |l_robust_m|=" << scb.l_robust_m.size() << ")";
  return os;
}

template <typename T, typename Tint>
void WriteEntryGAP(std::ostream &os_out, SoftConvexBoundary<T, Tint> const &scb) {
  os_out << "rec(index_cb:=" << scb.index_cb << ", cb:=";
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

namespace boost::serialization {

template <class Archive, typename T, typename Tint>
inline void serialize(Archive &ar, SoftConvexBoundary<T, Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("index_cb", val.index_cb);
  ar &make_nvp("cb", val.cb);
  ar &make_nvp("l_excluded_max", val.l_excluded_max);
  ar &make_nvp("l_robust_m", val.l_robust_m);
}

} // namespace boost::serialization

// clang-format off
#endif  // SRC_ROBUST_COVERING_ROBUST_COVERING_TYPES_H_
// clang-format on


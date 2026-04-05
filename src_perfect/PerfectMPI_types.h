// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_MPI_PERFECTMPI_TYPES_H_
#define SRC_PERFECT_MPI_PERFECTMPI_TYPES_H_

#include "MAT_Matrix.h"
#include <functional>
#include <string>
#include <vector>
// #include <boost/container_hash/hash.hpp>

struct TypeIndex {
  int iProc;
  int idxMatrix;
  int iAdj;
};

template <typename T> struct TypePerfectExch {
  // the number of shortest vectors divided by 2
  int incd;
  MyMatrix<T> eMat;
};

template <typename T> struct PairExch {
  TypePerfectExch<T> ePerfect;
  TypeIndex eIndex;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, TypePerfectExch<T> const &obj) {
  os << obj.incd;
  int nbRow = obj.eMat.rows();
  os << " " << nbRow;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    for (int iCol = 0; iCol <= iRow; iCol++)
      os << " " << obj.eMat(iRow, iCol);
  }
  return os;
}

template <typename T>
bool operator==(TypePerfectExch<T> const &obj1,
                TypePerfectExch<T> const &obj2) {
  if (obj1.incd != obj2.incd)
    return false;
  int nbRow1 = obj1.eMat.rows();
  int nbRow2 = obj2.eMat.rows();
  if (nbRow1 != nbRow2)
    return false;
  for (int iRow = 0; iRow < nbRow1; iRow++)
    for (int iCol = 0; iCol < nbRow1; iCol++)
      if (obj1.eMat(iRow, iCol) != obj2.eMat(iRow, iCol))
        return false;
  return true;
}

std::ostream &operator<<(std::ostream &os, TypeIndex const &obj) {
  os << obj.iProc << " " << obj.idxMatrix << " " << obj.iAdj;
  return os;
}

namespace std {
template <typename T> struct less<TypePerfectExch<T>> {
  bool operator()(TypePerfectExch<T> const &eTPE1,
                  TypePerfectExch<T> const &eTPE2) const {
    if (eTPE1.incd < eTPE2.incd)
      return true;
    if (eTPE1.incd > eTPE2.incd)
      return false;
    //
    int nbRow = eTPE1.eMat.rows();
    for (int iRow = 0; iRow < nbRow; iRow++)
      for (int iCol = 0; iCol < nbRow; iCol++) {
        if (eTPE1.eMat(iRow, iCol) < eTPE2.eMat(iRow, iCol))
          return true;
        if (eTPE1.eMat(iRow, iCol) > eTPE2.eMat(iRow, iCol))
          return false;
      }
    return false;
  }
};
// clang-format off
}  // namespace std
// clang-format on

namespace boost {
namespace serialization {
// TypePerfectExch
template <class Archive, typename T>
inline void serialize(Archive &ar, TypePerfectExch<T> &eRecMat,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("incd", eRecMat.incd);
  int rows = eRecMat.eMat.rows();
  int cols = eRecMat.eMat.cols();
  ar &make_nvp("rows", rows);
  ar &make_nvp("cols", cols);
  eRecMat.eMat.resize(rows, cols);
  for (int r = 0; r < rows; ++r)
    for (int c = 0; c < cols; ++c)
      ar &make_nvp("val", eRecMat.eMat(r, c));
}

// TypePerfectExch
template <class Archive>
inline void serialize(Archive &ar, TypeIndex &eTypIdx,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("iProc", eTypIdx.iProc);
  ar &make_nvp("idxMatrix", eTypIdx.idxMatrix);
  ar &make_nvp("iAdj", eTypIdx.iAdj);
}

// PairExch
template <class Archive, typename T>
inline void serialize(Archive &ar, PairExch<T> &ePair,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("perfect", ePair.ePerfect);
  ar &make_nvp("index", ePair.eIndex);
}

// clang-format off
}  // namespace serialization
}  // namespace boost
// clang-format on

namespace std {
template <typename Tint> struct hash<TypePerfectExch<Tint>> {
  std::size_t operator()(const TypePerfectExch<Tint> &k) const {
    std::size_t h1 = std::hash<int>()(k.incd);
    int nbRow = k.eMat.rows();
    for (int iRow = 0; iRow < nbRow; iRow++)
      for (int iCol = 0; iCol < nbRow; iCol++) {
        Tint eVal = k.eMat(iRow, iCol);
        std::size_t h2 = std::hash<Tint>()(eVal);
        h1 = h2 ^ (h1 << 1);
      }
    return h1;
  }
};

// clang-format off
}  // namespace std
// clang-format on

// clang-format off
#endif  // SRC_PERFECT_MPI_PERFECTMPI_TYPES_H_
// clang-format on

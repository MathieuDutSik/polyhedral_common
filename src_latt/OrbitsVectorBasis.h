// Copyright (C) 2024 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ORBITSVECTORBASIS_H_
#define SRC_LATT_ORBITSVECTORBASIS_H_

// clang-format off
#include "OrbitEnumeration.h"
#include <unordered_map>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ORBITS_VECTOR_BASIS
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ORBITS_VECTOR_BASIS
#endif

#ifdef TIMINGS
#define TIMINGS_ORBITS_VECTOR_BASIS
#endif

template<typename T>
bool is_canonical(MyVector<T> const& V) {
  int n = V.size();
  for (int i=0; i<n; i++) {
    T const& val = V(i);
    if (val != 0) {
      return val > 0;
    }
  }
  std::cerr << "The vector V is zero\n";
  throw TerminalException{1};
}

template<typename T>
void canonicalize_sign_vector(MyVector<T> & V) {
  int n = V.size();
  for (int i=0; i<n; i++) {
    T const& val = V(i);
    if (val != 0) {
      if (val < 0) {
        for (int j=i; j<n; j++) {
          V(j) = - V(j);
        }
      }
      return;
    }
  }
}


template <typename Tgroup, typename Tint>
vectface EnumerateOrbitBasis(MyMatrix<Tint> const &SHV,
                             std::vector<MyMatrix<Tint>> const &ListGen,
                             [[maybe_unused]] std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<MyVector<Tint>> ListVert;
  std::unordered_map<MyVector<Tint>, int> MapVert;
  int pos = 0;
  for (int i = 0; i < SHV.rows(); i++) {
    MyVector<Tint> V = GetMatrixRow(SHV, i);
    if (is_canonical(V)) {
      ListVert.push_back(V);
      pos += 1;
      MapVert[V] = pos;
    }
  }
  Tidx nbVert = ListVert.size();
  int dim = SHV.cols();
  std::vector<Telt> ListPermGen;
  for (auto &eGen : ListGen) {
    std::vector<Tidx> eList(nbVert);
    for (Tidx i = 0; i < nbVert; i++) {
      MyVector<Tint> const &V = ListVert[i];
      MyVector<Tint> Vimg = eGen.transpose() * V;
      canonicalize_sign_vector(Vimg);
      int pos = MapVert[Vimg];
#ifdef SANITY_CHECK_ORBITS_VECTOR_BASIS
      if (pos == 0) {
        std::cerr << "pos should be strictly positive\n";
        throw TerminalException{1};
      }
#endif
      eList[i] = pos - 1;
    }
    Telt ePermGen(eList);
#ifdef DEBUG_ORBITS_VECTOR_BASIS
    os << "ePermGen=" << ePermGen << "\n";
#endif
    ListPermGen.push_back(ePermGen);
  }
  Tgroup GRP(ListPermGen, nbVert);
  vectface vf(nbVert);
  auto f_extensible = [&](std::vector<Tidx> const &v) -> bool {
    int n_vert = v.size();
    MyMatrix<Tint> Mtest(n_vert, dim);
    for (int iRow = 0; iRow < n_vert; iRow++) {
      int idx = v[iRow];
      for (int iCol = 0; iCol < dim; iCol++) {
        Mtest(iRow, iCol) = ListVert[idx](iCol);
      }
    }
    if (RankMat(Mtest) != n_vert) {
      return false;
    }
    if (n_vert == dim) {
      Face f(nbVert);
      for (auto &val : v) {
        f[val] = 1;
      }
      vf.push_back(f);
      return false;
    }
    return true;
  };
  SubsetOrbitEnumeration(GRP, f_extensible, os);
  return vf;
}

// clang-format off
#endif  // SRC_LATT_ORBITSVECTORBASIS_H_
// clang-format on

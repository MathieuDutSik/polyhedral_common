// Copyright (C) 2024 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ORBITSVECTORBASIS_H_
#define SRC_LATT_ORBITSVECTORBASIS_H_

// clang-format off
#include "OrbitEnumeration.h"
// clang-format on

#ifdef DEBUG
#define DEBUG_ORBITS_VECTOR_BASIS
#endif

template <typename Tgroup, typename Tint>
vectface EnumerateOrbitBasis(MyMatrix<Tint> const &SHV,
                             std::vector<MyMatrix<Tint>> const &ListGen,
                             [[maybe_unused]] std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx nbVert = SHV.rows();
  std::vector<MyVector<Tint>> ListVert;
  std::unordered_map<MyVector<Tint>, int> MapVert;
  for (int i = 0; i < SHV.rows(); i++) {
    MyVector<Tint> V = GetMatrixRow(SHV, i);
    ListVert.push_back(V);
    MapVert[V] = i + 1;
  }
  int dim = SHV.cols();
  std::vector<Telt> ListPermGen;
  for (auto &eGen : ListGen) {
    std::vector<Tidx> eList(nbVert);
    for (Tidx i = 0; i < nbVert; i++) {
      MyVector<Tint> const &V = ListVert[i];
      MyVector<Tint> Vimg = eGen.transpose() * V;
      int pos = MapVert[Vimg];
#ifdef DEBUG_ORBITS_VECTOR_BASIS
      if (pos == 0) {
        std::cerr << "pos should be strictly positive\n";
        throw TerminalException{1};
      }
#endif
      eList[i] = pos - 1;
    }
    Telt ePermGen(eList);
    ListPermGen.push_back(ePermGen);
  }
  Tgroup GRP(ListPermGen, nbVert);
  vectface vf(nbVert);
  auto f_extensible = [&](std::vector<Tidx> const &v) -> bool {
    int n_vert = v.size();
    MyMatrix<Tint> Mtest(n_vert, dim);
    for (int iRow = 0; iRow < n_vert; iRow++) {
      for (int iCol = 0; iCol < dim; iCol++) {
        Mtest(iRow, iCol) = SHV(v[iRow], iCol);
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
  SubsetOrbitEnumeration(GRP, f_extensible);
  return vf;
}

// clang-format off
#endif  // SRC_LATT_ORBITSVECTORBASIS_H_
// clang-format on

// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_TESTGROUP_H_
#define SRC_GROUP_TESTGROUP_H_

// clang-format off
#include "COMB_Combinatorics_buildset.h"
#include <algorithm>
#include <unordered_map>
#include <vector>
// clang-format on

// Computes the residue modulo N.
// This is used for debugging and allows getting finite
// groups that can test membership.
template <typename T, typename Tgroup>
Tgroup GenerateGroupModuloAction(std::vector<MyMatrix<T>> const& ListM, int const& N) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  if (ListM.size() == 0) {
    std::cerr << "COMB: Not possible to work if zero vectors are available\n";
    throw TerminalException{1};
  }
  int n = ListM[0].rows();
  MyMatrix<int> Mat_cos = BuildSet(n, N);
  std::vector<MyVector<int>> l_cos;
  std::unordered_map<MyVector<int>,size_t> map_cos;
  int n_row = Mat_cos.rows();
  for (int i_row=0; i_row<n_row; i_row++) {
    MyVector<int> eRow = GetMatrixRow(Mat_cos, i_row);
    l_cos.push_back(eRow);
    map_cos[eRow] = i_row;
  }
  Tidx n_act = n_row;
  T N_T = UniversalScalarConversion<T,int>(N);
  auto get_perm=[&](MyMatrix<T> const& eM) -> Telt {
    MyMatrix<int> eMred(n, n);
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        T res_T = ResInt(eM(i,j), N_T);
        int res = UniversalScalarConversion<int, T>(res_T);
        eMred(i,j) = res;
      }
    }
    MyVector<int> eV3(n);
    std::vector<Tidx> eList;
    for (int i_row=0; i_row<n_row; i_row++) {
      MyVector<int> eV1 = l_cos[i_row];
      MyVector<int> eV2 = eMred.transpose() * eV1;
      for (int i=0; i<n; i++) {
        int res = ResInt(eV2(i), N);
        eV3(i) = res;
      }
      size_t pos = map_cos.at(eV3);
      eList.push_back(pos);
    }
    Telt ePerm(eList);
    return ePerm;
  };
  std::vector<Telt> ListPerm;
  for (auto & eM : ListM) {
    Telt ePerm = get_perm(eM);
    ListPerm.push_back(ePerm);
  }
  return Tgroup(ListPerm, n_act);
}


template <typename T, typename Tgroup>
void CheckSubgroupInclusion(std::vector<MyMatrix<T>> const& ListGRP, std::vector<MyMatrix<T>> const& ListSubGRP, [[maybe_unused]] std::ostream& os) {
  int n = ListGRP[0].rows();
  T limit(10000);
  for (int N=2; N<=20; N++) {
    T N_T = UniversalScalarConversion<T,int>(N);
    T Npow = MyPow(N_T, n);
    if (Npow > limit) {
      break;
    }
    Tgroup GRP = GenerateGroupModuloAction<T,Tgroup>(ListGRP, N);
    Tgroup SubGRP = GenerateGroupModuloAction<T,Tgroup>(ListSubGRP, N);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: N=" << N << " Npow=" << Npow << " |GRP|=" << GRP.size() << " |SubGRP|=" << SubGRP.size() << "\n";
#endif
    bool test = GRP.IsSubgroup(SubGRP);
    if (!test) {
      std::cerr << "COMB: Found SubGRP not to be a subgroup of GRP\n";
      throw TerminalException{1};
    }
  }
}

template <typename T, typename Tgroup>
void CheckGroupEquality(std::vector<MyMatrix<T>> const& ListGens1, std::vector<MyMatrix<T>> const& ListGens2, [[maybe_unused]] std::ostream& os) {
  if (ListGens1.size() == 0) {
    return;
  }
  int n = ListGens1[0].rows();
  T limit(10000);
  for (int N=2; N<=20; N++) {
    T N_T = UniversalScalarConversion<T,int>(N);
    T Npow = MyPow(N_T, n);
    if (Npow > limit) {
      break;
    }
    Tgroup GRP1 = GenerateGroupModuloAction<T,Tgroup>(ListGens1, N);
    Tgroup GRP2 = GenerateGroupModuloAction<T,Tgroup>(ListGens2, N);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: N=" << N << " Npow=" << Npow << " |GRP1|=" << GRP1.size() << " |GRP2|=" << GRP2.size() << "\n";
#endif
    bool test1 = GRP1.IsSubgroup(GRP2);
    bool test2 = GRP2.IsSubgroup(GRP1);
    if (!test1 || !test2) {
      std::cerr << "COMB: Found test1=" << test1 << " test2=" << test2 << "\n";
      throw TerminalException{1};
    }
  }
}

// clang-format off
#endif  // SRC_GROUP_TESTGROUP_H_
// clang-format on

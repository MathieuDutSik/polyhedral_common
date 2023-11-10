// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_TESTGROUP_H_
#define SRC_GROUP_TESTGROUP_H_

// clang-format off
#include "Group.h"
#include "MatrixGroup.h"
#include "Permutation.h"
#include "Temp_PolytopeEquiStab.h"
#include <algorithm>
#include <unordered_map>
#include <vector>
// clang-format on

// Off the normal track but ok. Just used for testing purposes.

using Tidx = int32_t;
using Tint = mpz_class;
const bool use_scheme1 = true;
using Telt1 = permutalib::SingleSidedPerm<Tidx>;
using Tgroup1 = permutalib::Group<Telt1, Tint>;

template <typename T>
void TestPolytopeFace_Stabilizer(MyMatrix<T> const &M, Face const &f) {
  //
  std::cerr << "---------------------------------------- TestPolytopeFace "
               "-------------------------------------\n";
  int n_vert = M.rows();
  MyMatrix<T> id_matr = IdentityMat<T>(M.cols());
  std::cerr << "GRP_ComputeAut_ListMat_Subset_EXT : Reading input\n";
  Tgroup1 GRP1 =
      LinPolytope_Automorphism<T, use_scheme1, Tgroup1>(M, std::cerr);
  std::cerr << "We have |GRP1|=" << GRP1.size() << "  n_vert=" << n_vert
            << "\n";
  Face fo = OrbitUnion(GRP1, f);
  std::cerr << "|f|=" << f.size() << " / " << f.count()
            << "    |fo|=" << fo.size() << " / " << fo.count() << "\n";
  std::vector<int> o_v;
  for (int i = 0; i < n_vert; i++)
    if (fo[i] == 1)
      o_v.push_back(i);
  std::vector<int> o_v_rev(n_vert, -1);
  for (size_t pos = 0; pos < o_v.size(); pos++) {
    int val = o_v[pos];
    //    std::cerr << "pos=" << pos << " val=" << val << "\n";
    o_v_rev[val] = pos;
  }
  /*
  std::cerr << "o_v =";
  for (auto & val : o_v)
    std::cerr << " " << val;
  std::cerr << "\n";
  std::cerr << "o_v_rev =";
  for (auto & val : o_v_rev)
    std::cerr << " " << val;
  std::cerr << "\n";
  std::cerr << "We have fo, o_v_rev\n";
  */
  std::unordered_map<MyVector<T>, int> map;
  for (int i_vert = 0; i_vert < n_vert; i_vert++) {
    MyVector<T> V = GetMatrixRow(M, i_vert);
    map[V] = i_vert;
  }
  //  std::cerr << "We have map\n";
  //
  auto get_group_A = [&](std::vector<MyMatrix<T>> const LMat,
                         std::vector<Telt1> const &LElt,
                         Face const &f) -> Tgroup1 {
    //    std::cerr << "Before computation of LStabMatrGen\n";
    std::vector<MyMatrix<T>> LStabMatrGen =
        permutalib::StabilizerMatrixPermSubset<Telt1, MyMatrix<T>, mpz_class>(
            LMat, LElt, id_matr, f);
    //    std::cerr << " After computation of LStabMatrGen\n";
    std::vector<Telt1> LStabPermGen;
    for (auto &eGen : LStabMatrGen) {
      std::vector<Tidx> eList(n_vert);
      for (int i_vert = 0; i_vert < n_vert; i_vert++) {
        MyVector<T> V = GetMatrixRow(M, i_vert);
        MyVector<T> Vimg = eGen.transpose() * V;
        eList[i_vert] = map[Vimg];
      }
      Telt1 elt1(eList);
      LStabPermGen.push_back(elt1);
    }
    return Tgroup1(LStabPermGen, n_vert);
  };
  //
  //
  std::vector<MyMatrix<T>> ListMatrGens;
  std::vector<Telt1> ListPermGens_A = GRP1.GeneratorsOfGroup();
  std::vector<Telt1> ListPermGens_B;
  for (auto &ePerm : ListPermGens_A) {
    MyMatrix<T> eMatr = FindTransformation(M, M, ePerm);
    ListMatrGens.push_back(eMatr);
    std::vector<Tidx> eList(o_v.size());
    for (size_t i = 0; i < o_v.size(); i++) {
      Tidx val1 = o_v[i];
      Tidx val2 = OnPoints(val1, ePerm);
      Tidx val3 = o_v_rev[val2];
      eList[i] = val3;
    }
    Telt1 elt_B(eList);
    ListPermGens_B.push_back(elt_B);
  }
  std::cerr << "We have ListMatrGens, ListPermGens_B\n";
  Face f_B(o_v.size());
  for (size_t pos = 0; pos < o_v.size(); pos++) {
    f_B[pos] = f[o_v[pos]];
  }

  Tgroup1 stab1 = get_group_A(ListMatrGens, ListPermGens_A, f);
  std::cerr << "We have stab1 |stab1|=" << stab1.size() << "\n";
  Tgroup1 stab2 = get_group_A(ListMatrGens, ListPermGens_B, f_B);
  std::cerr << "We have stab2 |stab2|=" << stab2.size() << "\n";
  Tgroup1 stab3 = GRP1.Stabilizer_OnSets(f);
  std::cerr << "We have stab3 |stab3|=" << stab3.size() << "\n";
  if (stab1 != stab2 || stab1 != stab3) {
    std::cerr << "The groups are not equal. Please debug\n";
    throw TerminalException{1};
  }
}

bool IsEqual(Face const &f1, Face const &f2) {
  for (size_t i = 0; i < f1.size(); i++)
    if (f1[i] != f2[i])
      return false;
  return true;
}

template <typename T>
void TestPolytopeFace_Equivalence(MyMatrix<T> const &M, Face const &f1) {
  std::cerr << "---------------------------------------- TestPolytopeFace "
               "-------------------------------------\n";
  int n_vert = M.rows();
  MyMatrix<T> id_matr = IdentityMat<T>(M.cols());
  std::cerr << "GRP_ComputeAut_ListMat_Subset_EXT : Reading input\n";
  Tgroup1 GRP1 =
      LinPolytope_Automorphism<T, use_scheme1, Tgroup1>(M, std::cerr);
  Telt1 elt1 = GRP1.rand();
  Face f2 = permutalib::OnSets(f1, elt1);
  std::cerr << "We have |GRP1|=" << GRP1.size() << "  n_vert=" << n_vert
            << "\n";
  Face fo = OrbitUnion(GRP1, f1);
  std::cerr << "|f1|=" << f1.size() << " / " << f1.count()
            << "    |fo|=" << fo.size() << " / " << fo.count() << "\n";
  std::vector<int> o_v;
  for (int i = 0; i < n_vert; i++)
    if (fo[i] == 1)
      o_v.push_back(i);
  std::vector<int> o_v_rev(n_vert, -1);
  for (size_t pos = 0; pos < o_v.size(); pos++) {
    int val = o_v[pos];
    o_v_rev[val] = pos;
  }
  std::unordered_map<MyVector<T>, int> map;
  for (int i_vert = 0; i_vert < n_vert; i_vert++) {
    MyVector<T> V = GetMatrixRow(M, i_vert);
    map[V] = i_vert;
  }
  //  std::cerr << "We have map\n";
  //
  std::vector<Tidx> ListPoint_A(n_vert);
  for (Tidx i = 0; i < n_vert; i++)
    ListPoint_A[i] = i;
  std::vector<Tidx> ListPoint_B;
  for (auto &val : o_v)
    ListPoint_B.push_back(val);
  //
  auto test_representation = [&](std::vector<MyMatrix<T>> const LMat,
                                 std::vector<Telt1> const &LElt,
                                 std::vector<Tidx> const &ListPoint,
                                 Face const &f1_i, Face const &f2_i) -> void {
    //    std::cerr << "Before computation of LStabMatrGen\n";
    Tgroup1 GRPrepr(LElt, f1_i.size());
    std::cerr << "|GRPrepr|=" << GRPrepr.size() << "\n";
    std::optional<MyMatrix<T>> opt =
        permutalib::RepresentativeActionMatrixPermSubset<Telt1, MyMatrix<T>,
                                                         mpz_class>(
            LMat, LElt, id_matr, f1_i, f2_i);
    if (!opt) {
      std::cerr << "Failed to find an equivalence\n";
      throw TerminalException{1};
    }
    MyMatrix<T> const &eEquiv = *opt;
    //    std::cerr << "We have eEquiv\n";
    //
    std::vector<Tidx> ListPoint_Rev(n_vert);
    for (size_t i = 0; i < ListPoint.size(); i++) {
      ListPoint_Rev[ListPoint[i]] = i;
    }
    //    std::cerr << "We have ListPoint_Rev\n";
    //    std::cerr << "eEquiv=\n";
    //    WriteMatrix(std::cerr, eEquiv);
    //
    int len = ListPoint.size();
    std::vector<Tidx> eList(len);
    for (int i_vert = 0; i_vert < len; i_vert++) {
      MyVector<T> V = GetMatrixRow(M, ListPoint[i_vert]);
      //      std::cerr << "i_vert=" << i_vert << " V=" << V << "\n";
      MyVector<T> Vimg = eEquiv.transpose() * V;
      eList[i_vert] = ListPoint_Rev[map[Vimg]];
    }
    //    std::cerr << "We have eList\n";
    Telt1 elt1(eList);
    Face f1_img = permutalib::OnSets(f1_i, elt1);
    if (f1_img != f2_i) {
      std::cerr << "|f1_i|=" << f1_i.size() << "\n";
      std::cerr << "  f1_i=" << f1_i << "\n";
      std::cerr << "  f2_i=" << f2_i << "\n";
      std::cerr << "f1_img=" << f1_img << "\n";
      std::cerr << "eList=" << eList << "\n";
      std::vector<Tidx> fList = eList;
      std::sort(fList.begin(), fList.end());
      std::cerr << "fList=" << fList << "\n";
      std::cerr << "elt1=" << elt1 << "\n";
      std::cerr << "IsEqual=" << IsEqual(f1_img, f2_i) << "\n";
      std::cerr << "The obtained permutation is not an equivalence\n";
      throw TerminalException{1};
    }
  };
  //
  //
  std::vector<MyMatrix<T>> ListMatrGens;
  std::vector<Telt1> ListPermGens_A = GRP1.GeneratorsOfGroup();
  std::vector<Telt1> ListPermGens_B;
  for (auto &ePerm : ListPermGens_A) {
    MyMatrix<T> eMatr = FindTransformation(M, M, ePerm);
    ListMatrGens.push_back(eMatr);
    std::vector<Tidx> eList(o_v.size());
    for (size_t i = 0; i < o_v.size(); i++) {
      Tidx val1 = o_v[i];
      Tidx val2 = OnPoints(val1, ePerm);
      Tidx val3 = o_v_rev[val2];
      eList[i] = val3;
    }
    Telt1 elt_B(eList);
    ListPermGens_B.push_back(elt_B);
  }
  std::cerr << "We have ListMatrGens, ListPermGens_B\n";
  Face f1_B(o_v.size());
  Face f2_B(o_v.size());
  for (size_t pos = 0; pos < o_v.size(); pos++) {
    f1_B[pos] = f1[o_v[pos]];
    f2_B[pos] = f2[o_v[pos]];
  }

  test_representation(ListMatrGens, ListPermGens_A, ListPoint_A, f1, f2);
  test_representation(ListMatrGens, ListPermGens_B, ListPoint_B, f1_B, f2_B);
}

// clang-format off
#endif  // SRC_GROUP_TESTGROUP_H_
// clang-format on

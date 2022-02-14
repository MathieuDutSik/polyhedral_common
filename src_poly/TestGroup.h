#ifndef INCLUDE_TEST_GROUP_H
#define INCLUDE_TEST_GROUP_H

#include "Permutation.h"
#include "Group.h"
#include "MatrixGroup.h"
#include "Temp_PolytopeEquiStab.h"


template<typename T>
void TestPolytopeFace(MyMatrix<T> const& M, Face const& f)
{
  using Tidx = int32_t;
  using Tint = mpz_class;
  const bool use_scheme1 = true;
  using Telt1   = permutalib::SingleSidedPerm<Tidx>;
  using Tgroup1 = permutalib::Group<Telt1,Tint>;
  //
  int n_vert = M.rows();
  MyMatrix<T> id_matr = IdentityMat<T>(n_vert);
  std::cerr << "GRP_ComputeAut_ListMat_Subset_EXT : Reading input\n";
  std::vector<int> f_v;
  for (int i=0; i<n_vert; i++)
    if (f[i] == 1)
      f_v.push_back(i);
  Tgroup1 GRP1 = LinPolytope_Automorphism<T,use_scheme1,Tgroup1>(M);
  std::vector<int> o_v = OrbitUnion(GRP1, f_v);
  Face fo(n_vert);
  std::vector<int> o_v_rev(n_vert,n_vert + 1000);
  for (size_t pos=0; pos<o_v.size(); pos++) {
    int val = o_v[pos];
    fo[val] = 1;
    o_v_rev[val] = pos;
  }
  std::unordered_map<MyVector<T>, int> map;
  for (int i_vert=0; i_vert<n_vert; i_vert++) {
    MyVector<T> V = GetMatrixRow(M, i_vert);
    map[V] = i_vert;
  }
  //
  auto get_group_A=[&](std::vector<MyMatrix<T>> const LMat, std::vector<Telt1> const& LElt, Face const& f) -> Tgroup1 {
    std::vector<MyMatrix<T>> LStabMatrGen = permutalib::StabilizerMatrixPermSubset<Telt1,MyMatrix<T>,mpz_class>(LMat, LElt, id_matr, f);
    std::vector<Telt1> LStabPermGen;
    for (auto & eGen : LStabMatrGen) {
      std::vector<Tidx> eList(n_vert);
      for (int i_vert=0; i_vert<n_vert; i_vert++) {
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
  for (auto & ePerm : ListPermGens_A) {
    MyMatrix<T> eMatr = FindTransformation(M, M, ePerm);
    ListMatrGens.push_back(eMatr);
    std::vector<Tidx> eList(o_v.size());
    for (size_t i=0; i<o_v.size(); i++) {
      Tidx val1 = o_v[i];
      Tidx val2 = OnPoints(val1, ePerm);
      Tidx val3 = o_v_rev[val2];
      eList[i] = val3;
    }
    Telt1 elt_B(eList);
    ListPermGens_B.push_back(elt_B);
  }
  Face f_B(o_v.size());
  for (auto & val1 : f_v) {
    Tidx val2 = o_v_rev[val1];
    f_B[val2] = 1;
  }

  Tgroup1 stab1 = get_group_A(ListMatrGens, ListPermGens_A, f);
  Tgroup1 stab2 = get_group_A(ListMatrGens, ListPermGens_B, f_B);
  Tgroup1 stab3 = GRP1.Stabilizer_OnSets(f);
  if (stab1 != stab2 || stab1 != stab3) {
    std::cerr << "The groups are not equal. Please debug\n";
    throw TerminalException{1};
  }
}


#endif

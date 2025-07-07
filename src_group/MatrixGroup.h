// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_MATRIXGROUP_H_
#define SRC_GROUP_MATRIXGROUP_H_

// clang-format off
#include "GRP_GroupFct.h"
#include "Group.h"
#include "InvariantVectorFamily.h"
#include "MAT_MatrixInt.h"
#include "MAT_MatrixMod.h"
#include "ClassicLLL.h"
#include "PERM_Fct.h"
#include "Timings.h"
#include "factorizations.h"
#include "two_dim_lorentzian.h"
#include "MatrixGroupSimplification.h"
#include "MatrixGroupBasic.h"
#include "TestGroup.h"
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <set>
#include <map>
// clang-format on

#ifdef DEBUG
#define DEBUG_MATRIX_GROUP
#define DEBUG_DOUBLE_COSET_ENUM
#endif

#ifdef TIMINGS
#define TIMINGS_MATRIX_GROUP
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_MATRIX_GROUP
#define SANITY_CHECK_DOUBLE_COSET_ENUM
#endif

#ifdef TRACK_INFO
#define TRACK_INFO_MATRIX_GROUP
#endif

//

template <typename T, typename Telt>
struct PreImager_Finite {
  MyMatrix<T> EXTfaithful;
  MyMatrix<T> pre_image_elt(Telt const& elt) const {
    using Tidx = typename Telt::Tidx;
    Tidx nbRow_tidx = EXTfaithful.rows();
    std::vector<Tidx> v(nbRow_tidx);
    for (Tidx i = 0; i < nbRow_tidx; i++) {
      v[i] = OnPoints(i, elt);
    }
    Telt ePermB(std::move(v));
    std::optional<MyMatrix<T>> opt = FindTransformationRing(EXTfaithful, EXTfaithful, ePermB);
    if (!opt) {
      std::cerr << "MAT_GRP: matrix is rational but not integral, so no conversion possible\n";
      std::cerr << "MAT_GRP: likely a bug since this error cannot occur for other code paths\n";
      throw TerminalException{1};
    }
    return *opt;
  }
};

template <typename T, typename Telt, typename TintGroup> struct FiniteMatrixGroupHelper {
  using PreImager = PreImager_Finite<T,Telt>;
  using Tint = typename underlying_ring<T>::ring_type;
  using ThelperInt = FiniteMatrixGroupHelper<typename underlying_ring<T>::ring_type, Telt, TintGroup>;
  int n;
  MyMatrix<T> EXTfaithful;
  std::vector<MyVector<T>> ListV;
  std::unordered_map<MyVector<T>, int> MapV;
  int nbRow() const {
    return EXTfaithful.rows();
  }
  PreImager constant_pre_imager() const {
    return PreImager { EXTfaithful };
  }
  PreImager pre_imager([[maybe_unused]] std::vector<MyMatrix<T>> const& l_matr, [[maybe_unused]] std::vector<Telt> const& l_perm) const {
    return constant_pre_imager();
  }
};

//

template <typename T, typename Telt>
struct PreImager_FiniteIsotropic {
  MyMatrix<T> G;
  MyMatrix<T> EXTfaithful;
  PreImager_FiniteIsotropic(MyMatrix<T> const& _G, MyMatrix<T> const& _EXTfaithful) : G(_G), EXTfaithful(_EXTfaithful) {
  }
  MyMatrix<T> pre_image_elt(Telt const& elt) const {
    MyMatrix<T> const &Subspace1 = EXTfaithful;
    int n_rows = Subspace1.rows();
    int n_cols = Subspace1.cols();
    MyMatrix<T> Subspace2(n_rows, n_cols);
    for (int i_row = 0; i_row < n_rows; i_row++) {
      int j_row = OnPoints(i_row, elt);
      for (int i_col=0; i_col<n_cols; i_col++) {
        Subspace2(i_row, i_col) = Subspace1(j_row, i_col);
      }
    }
    std::optional<MyMatrix<T>> opt =
      LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1(G, Subspace1, G, Subspace2);
    return unfold_opt(opt, "The isotropic extension should work");
  }
};

template <typename T, typename Telt, typename TintGroup> struct FiniteIsotropicMatrixGroupHelper {
  using PreImager = PreImager_FiniteIsotropic<T,Telt>;
  using Tint = typename underlying_ring<T>::ring_type;
  using ThelperInt = FiniteIsotropicMatrixGroupHelper<typename underlying_ring<T>::ring_type, Telt, TintGroup>;
  int n;
  MyMatrix<T> G;
  MyMatrix<T> EXTfaithful;
  MyVector<T> Visotrop;
  std::vector<MyVector<T>> ListV;
  std::unordered_map<MyVector<T>, int> MapV;
  int nbRow() const {
    return EXTfaithful.rows();
  }
  PreImager constant_pre_imager() const {
    return PreImager(G, EXTfaithful);
  }
  PreImager pre_imager([[maybe_unused]] std::vector<MyMatrix<T>> const& l_matr, [[maybe_unused]] std::vector<Telt> const& l_perm) const {
    return constant_pre_imager();
  }
};

//

template <typename T, typename Telt, typename TintGroup>
struct PreImager_General {
private:
  permutalib::PreImagerElement<Telt,MyMatrix<T>, TintGroup> inner;
public:
  PreImager_General(std::vector<MyMatrix<T>> const& l_matr, std::vector<Telt> const& l_perm, int const& dim) : inner(l_matr, l_perm, IdentityMat<T>(dim)) {
#ifdef DEBUG_MATRIX_GROUP
    std::cerr << "MAT_GRP: After building inner\n";
#endif
  }
  MyMatrix<T> pre_image_elt(Telt const& elt) const {
    std::optional<MyMatrix<T>> opt = inner.get_preimage(elt);
    return unfold_opt(opt, "The element elt should belong to the group");
  }
};

template <typename T, typename Telt, typename TintGroup> struct GeneralMatrixGroupHelper {
  using PreImager = PreImager_General<T, Telt, TintGroup>;
  using Tint = typename underlying_ring<T>::ring_type;
  using ThelperInt = GeneralMatrixGroupHelper<typename underlying_ring<T>::ring_type, Telt, TintGroup>;
  int n;
  int nbRow() const {
    return 0;
  }
  PreImager pre_imager(std::vector<MyMatrix<T>> const& l_matr, std::vector<Telt> const& l_perm) const {
#ifdef DEBUG_MATRIX_GROUP
    std::cerr << "MAT_GRP: Before PreImager constructor\n";
#endif
    return PreImager(l_matr, l_perm, n);
  }
};

//

template <typename Thelper> struct has_determining_ext {
  static const bool value = false;
};

template <typename T, typename Telt, typename TintGroup>
struct has_determining_ext<FiniteMatrixGroupHelper<T, Telt, TintGroup>> {
  static const bool value = true;
};

template <typename T, typename Telt, typename TintGroup>
struct has_determining_ext<FiniteIsotropicMatrixGroupHelper<T, Telt, TintGroup>> {
  static const bool value = true;
};

//

template <typename T, typename Telt, typename TintGroup>
GeneralMatrixGroupHelper<T, Telt, TintGroup>
TransformHelper(GeneralMatrixGroupHelper<T, Telt, TintGroup> const &helper,
                [[maybe_unused]] MyMatrix<T> const &Pmat) {
  return helper;
}

template <typename T, typename Telt, typename TintGroup>
GeneralMatrixGroupHelper<typename underlying_ring<T>::ring_type, Telt, TintGroup>
ToInteger(GeneralMatrixGroupHelper<T, Telt, TintGroup> const &helper) {
  using Tint = typename underlying_ring<T>::ring_type;
  return GeneralMatrixGroupHelper<Tint, Telt, TintGroup>{helper.n};
}

template <typename T, typename Telt, typename TintGroup>
FiniteIsotropicMatrixGroupHelper<T, Telt, TintGroup>
TransformHelper(FiniteIsotropicMatrixGroupHelper<T, Telt, TintGroup> const &helper,
                MyMatrix<T> const &Pmat) {
  MyMatrix<T> PmatInv = Inverse(Pmat);
  MyMatrix<T> G_new = Pmat * helper.G * Pmat.transpose();
  MyMatrix<T> EXTfaithful_new = helper.EXTfaithful * PmatInv;
  MyVector<T> Visotrop_new = PmatInv.transpose() * helper.Visotrop;
  std::vector<MyVector<T>> ListV_new;
  std::unordered_map<MyVector<T>, int> MapV_new;
  int len = EXTfaithful_new.rows();
  for (int i = 0; i < len; i++) {
    MyVector<T> eV = GetMatrixRow(EXTfaithful_new, i);
    MapV_new[eV] = i;
    ListV_new.emplace_back(std::move(eV));
  }
  return {helper.n,
          std::move(G_new),
          std::move(EXTfaithful_new),
          std::move(Visotrop_new),
          std::move(ListV_new),
          std::move(MapV_new)};
}


template <typename T, typename Telt, typename TintGroup>
FiniteIsotropicMatrixGroupHelper<typename underlying_ring<T>::ring_type, Telt, TintGroup>
ToInteger(FiniteIsotropicMatrixGroupHelper<T, Telt, TintGroup> const &helper) {
  using Tint = typename underlying_ring<T>::ring_type;
  //
  MyMatrix<T> G2 = RemoveFractionMatrix(helper.G);
  MyMatrix<Tint> G3 = UniversalMatrixConversion<Tint,T>(G2);
  //
  MyMatrix<T> EXTfaithful2 = RemoveFractionMatrix(helper.EXTfaithful);
  MyMatrix<Tint> EXTfaithful3 = UniversalMatrixConversion<Tint,T>(EXTfaithful2);
  //
  MyVector<T> Visotrop2 = RemoveFractionVector(helper.Visotrop);
  MyVector<Tint> Visotrop3 = UniversalVectorConversion<Tint,T>(Visotrop2);
  //
  std::vector<MyVector<Tint>> ListV;
  std::unordered_map<MyVector<Tint>, int> MapV;
  int len = EXTfaithful3.rows();
  for (int i = 0; i < len; i++) {
    MyVector<Tint> eV = GetMatrixRow(EXTfaithful3, i);
    MapV[eV] = i;
    ListV.emplace_back(std::move(eV));
  }
  return {helper.n,
          std::move(G3),
          std::move(EXTfaithful3),
          std::move(Visotrop3),
          std::move(ListV),
          std::move(MapV)};
}


template <typename T, typename Telt, typename TintGroup>
FiniteMatrixGroupHelper<T, Telt, TintGroup>
TransformHelper(FiniteMatrixGroupHelper<T, Telt, TintGroup> const &helper,
                MyMatrix<T> const &Pmat) {
  MyMatrix<T> PmatInv = Inverse(Pmat);
  MyMatrix<T> EXTfaithful_new = helper.EXTfaithful * PmatInv;
  std::vector<MyVector<T>> ListV_new;
  std::unordered_map<MyVector<T>, int> MapV_new;
  int len = EXTfaithful_new.rows();
  for (int i = 0; i < len; i++) {
    MyVector<T> eV = GetMatrixRow(EXTfaithful_new, i);
    ListV_new.push_back(eV);
    MapV_new[eV] = i;
  }
  return {helper.n, std::move(EXTfaithful_new), std::move(ListV_new),
          std::move(MapV_new)};
}

template <typename T, typename Telt, typename TintGroup>
FiniteMatrixGroupHelper<typename underlying_ring<T>::ring_type, Telt, TintGroup>
ToInteger(FiniteMatrixGroupHelper<T, Telt, TintGroup> const &helper) {
  using Tint = typename underlying_ring<T>::ring_type;
  //
  MyMatrix<T> EXTfaithful2 = RemoveFractionMatrix(helper.EXTfaithful);
  MyMatrix<Tint> EXTfaithful3 = UniversalMatrixConversion<Tint,T>(EXTfaithful2);
  //
  std::vector<MyVector<Tint>> ListV;
  std::unordered_map<MyVector<Tint>, int> MapV;
  int len = EXTfaithful3.rows();
  for (int i = 0; i < len; i++) {
    MyVector<Tint> eV = GetMatrixRow(EXTfaithful3, i);
    MapV[eV] = i;
    ListV.emplace_back(std::move(eV));
  }
  return {helper.n, std::move(EXTfaithful3), std::move(ListV),
          std::move(MapV)};
}

// We need a procedure that not just transform the helper, but also
// changes the type T to an integer type when at all possible.
// ----
// What could be done:
// ---We are working with a type underlying_ring already
// ---So, the transform_helper could map to that derived
//    type.
// ---The derived type could part of the Thelper, like
//    typename Thelper::ThelperInt;
// ---The conversion from T to Tint have to be done in a
//    separate function.
// ---The commands to encapsulate are:
//    + LinearSpace_Equivalence_Kernel
//    + LinearSpace_Stabilizer_Kernel
//    + LinearSpace_Stabilizer_RightCoset_Kernel
//    + LinearSpace_Stabilizer_DoubleCoset_Kernel
//    + LinearSpace_Stabilizer_DoubleCosetStabilizer_Kernel
// ---We need independent functions like
//    ToInteger(.....)

template <typename T, typename Telt, typename TintGroup>
FiniteMatrixGroupHelper<T, Telt, TintGroup>
ComputeFiniteMatrixGroupHelper(MyMatrix<T> const &EXT) {
  std::vector<MyVector<T>> ListV;
  std::unordered_map<MyVector<T>, int> MapV;
  for (int i = 0; i < EXT.rows(); i++) {
    MyVector<T> V = GetMatrixRow(EXT, i);
    ListV.push_back(V);
    MapV[V] = i;
  }
  int n_col = EXT.cols();
  return {n_col, EXT, std::move(ListV), std::move(MapV)};
}

template <typename T, typename Telt, typename TintGroup>
FiniteIsotropicMatrixGroupHelper<T, Telt, TintGroup>
ComputeFiniteIsotropicMatrixGroupHelper(MyMatrix<T> const &G,
                                        MyMatrix<T> const &EXT,
                                        MyVector<T> const &Visotrop) {
  std::vector<MyVector<T>> ListV;
  std::unordered_map<MyVector<T>, int> MapV;
  for (int i = 0; i < EXT.rows(); i++) {
    MyVector<T> V = GetMatrixRow(EXT, i);
    ListV.push_back(V);
    MapV[V] = i;
  }
  int n_col = EXT.cols();
  return {n_col, G, EXT, Visotrop, std::move(ListV), std::move(MapV)};
}

template <typename T, typename Telt, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value, Telt>::type
GetPermutationForFiniteMatrixGroup(Thelper const &helper,
                                   MyMatrix<T> const &eMatr,
                                   [[maybe_unused]] std::ostream &os) {
  using Tidx = typename Telt::Tidx;
  Tidx len = helper.EXTfaithful.rows();
  std::vector<Tidx> V(len);
  for (Tidx i = 0; i < len; i++) {
    MyVector<T> const &eV = helper.ListV[i];
    MyVector<T> Vimg = eMatr.transpose() * eV;
    V[i] = helper.MapV.at(Vimg);
  }
  return Telt(std::move(V));
}

template <typename T, typename Tmod>
Face GetFace(std::vector<MyVector<Tmod>> const &O,
             MyMatrix<T> const &TheSpace) {
  size_t Osiz = O.size();
  Face eFace(Osiz);
  RecSolutionIntMat<T> eCan(TheSpace);
  for (size_t iO = 0; iO < Osiz; iO++) {
    MyVector<T> const &eVect = UniversalVectorConversion<T, Tmod>(O[iO]);
    bool test = eCan.has_solution_v(eVect);
    if (test) {
      eFace[iO] = 1;
    }
  }
  return eFace;
}

Face TranslateFace(int const& nbRow, Face const& face) {
  if (nbRow == 0) {
    return face;
  }
  int siz = face.size();
  Face ftrans(siz + nbRow);
  for (int i=0; i<siz; i++) {
    ftrans[nbRow + i] = face[i];
  }
  return ftrans;
}


// The subsets can be grouped with a block decomposition
//
template<typename T, typename Telt>
struct PartitionReduction {
private:
  using Tidx = typename Telt::Tidx;
public:
  std::function<Telt(const MyMatrix<T>&)> f_get_perm;
  Face face;
private:
  PartitionStorage<Telt> partition;
public:
  PartitionReduction(std::vector<MyMatrix<T>> const& ListMat, std::function<Telt(const MyMatrix<T>&)> const& f_get_perm_inp, Face const& face_inp, std::ostream& os) : partition(face_inp, os) {
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: PartitionReduction, |face_inp|=" << face_inp.size() << " / " << face_inp.count() << "\n";
#endif
    std::vector<Telt> list_gens;
    for (auto& eMat: ListMat) {
      Telt egen = f_get_perm_inp(eMat);
      list_gens.push_back(egen);
    }
    partition.RefinePartitionByListElt(list_gens);
    face = partition.map_face(face_inp);
    f_get_perm = [f_get_perm_inp, this](MyMatrix<T> const& eMat) -> Telt {
      Telt g = f_get_perm_inp(eMat);
      return partition.map_permutation(g);
    };
  }
  std::optional<Face> map_face_opt(Face const& f) const {
    return partition.map_face_opt(f);
  }
};






template <typename T, typename Telt, typename Thelper, typename Fgetperm>
inline typename std::enable_if<has_determining_ext<Thelper>::value, Telt>::type
MatrixIntegral_MapMatrix(Thelper const &helper, Fgetperm f_get_perm, MyMatrix<T> const& eMatr, std::ostream &os) {
  using Tidx = typename Telt::Tidx;
  int nbRow = helper.EXTfaithful.rows();
  Tidx nbRow_tidx = nbRow;
  Telt ePermGen = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(helper, eMatr, os);
  Telt ePermOrbit = f_get_perm(eMatr);
  Tidx Osiz = ePermOrbit.size();
  Tidx siz = nbRow_tidx + Osiz;
  std::vector<Tidx> v(siz);
  for (Tidx i = 0; i < nbRow_tidx; i++) {
    v[i] = ePermGen.at(i);
  }
  for (Tidx iO = 0; iO < Osiz; iO++) {
    Tidx jO = ePermOrbit.at(iO);
    v[nbRow_tidx + iO] = nbRow_tidx + jO;
  }
  return Telt(std::move(v));
}

template <typename T, typename Telt, typename Thelper, typename Fgetperm>
inline typename std::enable_if<!has_determining_ext<Thelper>::value, Telt>::type
MatrixIntegral_MapMatrix([[maybe_unused]] Thelper const &helper, Fgetperm f_get_perm, MyMatrix<T> const& eMatr, [[maybe_unused]] std::ostream &os) {
  return f_get_perm(eMatr);
}

template <typename T, typename Telt, typename Thelper, typename Fgetperm>
std::vector<Telt> MatrixIntegral_GeneratePermutationGroupA(
    std::vector<MyMatrix<T>> const &ListMatrGens,
    Thelper const &helper, Fgetperm f_get_perm, std::ostream &os) {
#ifdef TIMINGS_MATRIX_GROUP
  MicrosecondTime time;
#endif
  std::vector<Telt> ListPermGenProv;
  size_t nbGen = ListMatrGens.size();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_GeneratePermutationGroupA, Processing nbGen=" << nbGen << " generators\n";
#endif
  for (size_t iGen = 0; iGen < nbGen; iGen++) {
    MyMatrix<T> const &eMatrGen = ListMatrGens[iGen];
    Telt eNewPerm = MatrixIntegral_MapMatrix<T,Telt,Thelper,Fgetperm>(helper, f_get_perm, eMatrGen, os);
    ListPermGenProv.emplace_back(std::move(eNewPerm));
  }
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: MatrixIntegral_GeneratePermutationGroupA|=" << time << "\n";
#endif
  return ListPermGenProv;
}

template<typename T, typename Tgroup>
struct RetMI_S {
  typename Tgroup::Tint index;
  std::vector<MyMatrix<T>> LGen;
};

// We have a finite set on which the group is acting. Therefore, we can apply
// the partition backtrack algorithms
template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value,
                               RetMI_S<T,Tgroup>>::type
MatrixIntegral_Stabilizer(std::vector<typename Tgroup::Telt> const &ListPermGens,
                          std::vector<MyMatrix<T>> const& ListMatr,
                          [[maybe_unused]] std::function<typename Tgroup::Telt(MyMatrix<T> const&)> f_get_perm,
                          Tgroup const &GRPperm, Thelper const &helper,
                          Face const &eFace, [[maybe_unused]] std::ostream &os) {
  using PreImager = typename Thelper::PreImager;
  using TintGroup = typename Tgroup::Tint;
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Beginning of MatrixIntegral_Stabilizer(has)\n";
#endif
#ifdef TIMINGS_MATRIX_GROUP
  MicrosecondTime time;
#endif
  Tgroup eStab = GRPperm.Stabilizer_OnSets(eFace);
  TintGroup index = GRPperm.size() / eStab.size();
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: MatrixIntegral_Stabilizer(has), Stabilizer_OnSets|=" << time << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: |eStab|=" << eStab.size() << " |eFace|=" << eFace.count() << "\n";
  os << "MAT_GRP: MatrixIntegral_Stabilizer(has), index=" << index << "\n";
#endif
  PreImager pre_imager = helper.pre_imager(ListMatr, ListPermGens);
  std::vector<MyMatrix<T>> LGen1;
  for (auto &eGen : eStab.SmallGeneratingSet()) {
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: MatrixIntegral_Stabilizer(has), before pre_image_elt\n";
#endif
    MyMatrix<T> eMatr = pre_imager.pre_image_elt(eGen);
    LGen1.emplace_back(std::move(eMatr));
  }
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: MatrixIntegral_Stabilizer(has), PreImager|=" << time << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_Stabilizer(has), comp(LGen1)=" << compute_complexity_listmat(LGen1) << "\n";
#endif
  std::vector<MyMatrix<T>> LGen2 = ExhaustiveReductionComplexityGroupMatrix<T>(LGen1, os);
#ifdef SANITY_CHECK_MATRIX_GROUP
  CheckGroupEquality<T,Tgroup>(LGen1, LGen2, os);
#endif
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_Stabilizer(has), comp(LGen2)=" << compute_complexity_listmat(LGen2) << "\n";
#endif
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: MatrixIntegral_Stabilizer(has), ExhaustiveReductionComplexityGroupMatrix|=" << time << "\n";
#endif
  return {index, LGen2};
}

// We have a finite set on which the group is acting. Therefore, we can apply
// the partition backtrack algorithms
template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<
    has_determining_ext<Thelper>::value,
    std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>>>::type
MatrixIntegral_Stabilizer_RightCoset(std::vector<typename Tgroup::Telt> const &ListPermGens,
                                     std::vector<MyMatrix<T>> const& ListMatr,
                                     Tgroup const &GRPperm,
                                     Thelper const &helper, Face const &eFace,
                                     [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin MatrixIntegral_Stabilizer_RightCoset(has)\n";
#endif
  using RightCosets = typename Tgroup::RightCosets;
  using PreImager = typename Thelper::PreImager;
  Tgroup eStab = GRPperm.Stabilizer_OnSets(eFace);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: |eStab|=" << eStab.size() << " |eFace|=" << eFace.count() << "\n";
#endif
  std::vector<MyMatrix<T>> ListMatrGen;
  PreImager pre_imager = helper.pre_imager(ListMatr, ListPermGens);
  for (auto &eGen : eStab.SmallGeneratingSet()) {
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: MatrixIntegral_Stabilizer_RightCoset(has), before pre_image_elt for eGen\n";
#endif
    MyMatrix<T> eMatr = pre_imager.pre_image_elt(eGen);
    ListMatrGen.emplace_back(std::move(eMatr));
  }
#ifdef TRACK_INFO_MATRIX_GROUP
  write_matrix_group(ListMatrGen, "MatrixIntegral_Stabilizer_RightCoset_has");
#endif
  std::vector<MyMatrix<T>> ListRightCoset;
  RightCosets rc = GRPperm.right_cosets(eStab);
  for (auto &eCos : rc) {
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: MatrixIntegral_Stabilizer_RightCoset(has), before pre_image_elt for eCos\n";
#endif
    MyMatrix<T> eMatr = pre_imager.pre_image_elt(eCos);
    ListRightCoset.emplace_back(std::move(eMatr));
  }
  return {std::move(ListMatrGen), std::move(ListRightCoset)};
}

template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value,
                               std::optional<MyMatrix<T>>>::type
MatrixIntegral_RepresentativeAction(std::vector<typename Tgroup::Telt> const &ListPermGens,
                                    std::vector<MyMatrix<T>> const& ListMatr,
                                    Tgroup const &GRPperm,
                                    Thelper const &helper, Face const &eFace1,
                                    Face const &eFace2, [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Beginning of MatrixIntegral_RepresentativeAction(has)\n";
#endif
  using Telt = typename Tgroup::Telt;
  using PreImager = typename Thelper::PreImager;
  std::optional<Telt> opt = GRPperm.RepresentativeAction_OnSets(eFace1, eFace2);
  if (!opt) {
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: Exit while loop with proof that no equivalence exists\n";
#endif
    return {};
  }
  PreImager pre_imager = helper.pre_imager(ListMatr, ListPermGens);
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: MatrixIntegral_RepresentativeAction(has), before pre_image_elt\n";
#endif
  return pre_imager.pre_image_elt(*opt);
}

template<typename T, typename Tmod, typename Telt>
Telt get_permutation_from_orbit(MyMatrix<T> const& eGen, std::vector<MyVector<Tmod>> const& O, T const& TheMod, Telt const& ePermS) {
  using Tidx = typename Telt::Tidx;
  Tmod TheMod_mod = UniversalScalarConversion<Tmod, T>(TheMod);
  size_t Osiz = O.size();
  std::vector<MyVector<Tmod>> ListImage(Osiz);
  MyMatrix<Tmod> eGenMod = ModuloReductionMatrix<T,Tmod>(eGen, TheMod);
  for (size_t iV = 0; iV < Osiz; iV++) {
    MyVector<Tmod> eVect = eGenMod.transpose() * O[iV];
    ListImage[iV] = VectorMod(eVect, TheMod_mod);
  }
  Telt ePermB = Telt(SortingPerm<MyVector<Tmod>, Tidx>(ListImage));
  Telt ePermBinv = ~ePermB;
  // By the construction and above check we have
  // V1reord[i] = V1[g1.at(i)]
  // V2reord[i] = V2[g2.at(i)]
  // We have V1reord = V2reord which gets us
  // V2[i] = V1[g1 * g2^{-1}(i)]
  Telt ePermGen = ePermBinv * ePermS;
  return ePermGen;
}

template <typename T, typename Tmod, typename Telt, typename Thelper>
std::vector<Telt> MatrixIntegral_GeneratePermutationGroup(
    std::vector<MyMatrix<T>> const &ListMatrGens,
    Thelper const &helper,
    std::vector<MyVector<Tmod>> const &O, T const &TheMod,
    std::ostream &os) {
  using Tidx = typename Telt::Tidx;
  Telt ePermS = Telt(SortingPerm<MyVector<Tmod>, Tidx>(O));
  auto f_get_perm=[&](MyMatrix<T> const& eGen) -> Telt {
    return get_permutation_from_orbit(eGen, O, TheMod, ePermS);
  };
  return MatrixIntegral_GeneratePermutationGroupA<T, Telt, Thelper, decltype(f_get_perm)>(ListMatrGens, helper, f_get_perm, os);
}

// We compute the stabilizer by applying the Schreier algorithm
template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<!has_determining_ext<Thelper>::value,
                               RetMI_S<T,Tgroup>>::type
MatrixIntegral_Stabilizer(std::vector<typename Tgroup::Telt> const &ListPermGens,
                          std::vector<MyMatrix<T>> const& ListMatrGens,
                          std::function<typename Tgroup::Telt(MyMatrix<T> const&)> f_get_perm,
                          [[maybe_unused]] Tgroup const &GRPperm,
                          Thelper const &helper, Face const &f,
                          [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_MATRIX_GROUP
  MicrosecondTime time;
#endif
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin MatrixIntegral_Stabilizer(!has)\n";
#endif
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  using Tidx = typename Telt::Tidx;
  MyMatrix<T> id_matr = IdentityMat<T>(helper.n);
  Tidx len = f.size();
  Tgroup GRP(ListPermGens, len);
  Tgroup stab = GRP.Stabilizer_OnSets(f);
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: MatrixIntegral_Stabilizer, Stabilizer_OnSets|=" << time << "\n";
#endif

  TintGroup index = GRP.size() / stab.size();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_Stabilizer(!has), index=" << index << "\n";
#endif
  Telt id_perm = stab.get_identity();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_Stabilizer(!has), comp(ListMatrGens)=" << compute_complexity_listmat(ListMatrGens) << "\n";
#endif
  std::vector<MyMatrix<T>> ListGen =
    PreImageSubgroup<T, Tgroup>(ListMatrGens, ListPermGens, f_get_perm, id_matr, stab, os);
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: MatrixIntegral_Stabilizer, PreImageSubgroupAction|=" << time << "\n";
#endif
  return {index, ListGen};
}

// We compute the stabilizer and right cosets by applying the Schreier algorithm
template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<
    !has_determining_ext<Thelper>::value,
    std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>>>::type
MatrixIntegral_Stabilizer_RightCoset(std::vector<typename Tgroup::Telt> const &ListPermGens,
                                     std::vector<MyMatrix<T>> const& ListMatr,
                                     [[maybe_unused]] Tgroup const &GRPperm,
                                     Thelper const &helper, Face const &eFace,
                                     [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin MatrixIntegral_Stabilizer(!has)\n";
#endif
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> id_matr = IdentityMat<T>(helper.n);
  std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>> pair =
      permutalib::StabilizerRightCosetMatrixPermSubset<Telt, MyMatrix<T>, TintGroup>(
          ListMatr, ListPermGens, id_matr, eFace);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: After StabilizerRightCosetMatrixPermSubset\n";
#endif
#ifdef TRACK_INFO_MATRIX_GROUP
  write_matrix_group(pair.first, "MatrixIntegral_Stabilizer_RightCoset_has_not_first");
  write_matrix_group(pair.second, "MatrixIntegral_Stabilizer_RightCoset_has_not_second");
#endif
  return pair;
}

template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<!has_determining_ext<Thelper>::value,
                               std::optional<MyMatrix<T>>>::type
MatrixIntegral_RepresentativeAction(std::vector<typename Tgroup::Telt> const &ListPermGens,
                                    std::vector<MyMatrix<T>> const& ListMatr,
                                    [[maybe_unused]] Tgroup const &GRPperm,
                                    Thelper const &helper, Face const &eFace1,
                                    Face const &eFace2,
                                    [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Beginning of MatrixIntegral_RepresentativeAction(!has)\n";
#endif
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> id_matr = IdentityMat<T>(helper.n);
  std::optional<MyMatrix<T>> opt =
      permutalib::RepresentativeActionMatrixPermSubset<Telt, MyMatrix<T>, TintGroup>(
          ListMatr, ListPermGens, id_matr, eFace1, eFace2);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Ending of MatrixIntegral_RepresentativeAction 2\n";
#endif
  return opt;
}

/*
  Direct computation of orbits.
  First level of optional is for termination or not.
  Second level is for whether we find an equivalence or not.
 */
template <typename T, typename Fterminate>
std::optional<std::optional<MyMatrix<T>>> DirectSpaceOrbit_Equivalence(
    std::vector<MyMatrix<T>> const &ListMatrGen, MyMatrix<T> const &eSpace1,
    MyMatrix<T> const &eSpace2, T const &TheMod, Fterminate const &f_terminate,
    [[maybe_unused]] std::ostream &os) {
  int n = eSpace1.rows();
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
  // Here Tpair is <Space,Repr>
  using Tpair = std::pair<MyMatrix<T>, MyMatrix<T>>;
  std::vector<Tpair> ListPair;
  ListPair.push_back({eSpace1, IdentityMat<T>(n)});
  if (f_terminate(eSpace1))
    return {};
  size_t start = 0;
  while (true) {
    size_t len = ListPair.size();
    if (start == len)
      break;
    for (size_t idx = start; idx < len; idx++) {
      Tpair const &ePair = ListPair[idx];
      for (auto &eMatrGen : ListMatrGen) {
        MyMatrix<T> eSpaceImg = ePair.first * eMatrGen;
        MyMatrix<T> eReprImg = ePair.second * eMatrGen;
        //
        MyMatrix<T> eSpaceMod = Concatenate(ePair.first, ModSpace);
        RecSolutionIntMat<T> eCan(eSpaceMod);
        if (eCan.is_containing_m(eSpace2)) {
          std::optional<MyMatrix<T>> opt = eReprImg;
          return opt;
        }
        auto fInsert = [&](Tpair const &ePair) -> bool {
          for (auto &fPair : ListPair)
            if (eCan.is_containing_m(fPair.first))
              return false;
          ListPair.push_back(ePair);
          return f_terminate(ePair.first);
        };
        if (fInsert(eSpaceImg))
          return {};
      }
    }
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: start=" << start << " len=" << len << "\n";
#endif
    start = len;
  }
  std::optional<MyMatrix<T>> opt;
  return opt;
}

template <typename T, typename Fterminate>
std::optional<std::vector<MyMatrix<T>>>
DirectSpaceOrbit_Stabilizer(std::vector<MyMatrix<T>> const &ListMatrGen,
                            MyMatrix<T> const &eSpace, T const &TheMod,
                            Fterminate const &f_terminate,
                            [[maybe_unused]] std::ostream &os) {
  int n = eSpace.rows();
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: DirectSpaceOrbit_Stabilizer, We have ModSpace\n";
#endif
  // Tpair is a pair <Space, Repr>
  using Tpair = std::pair<MyMatrix<T>, MyMatrix<T>>;
  std::vector<Tpair> ListPair;
  ListPair.push_back({eSpace, IdentityMat<T>(n)});
  if (f_terminate(eSpace)) {
    return {};
  }
  size_t start = 0;
  while (true) {
    size_t len = ListPair.size();
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: DirectSpaceOrbit_Stabilizer, start=" << start << " len=" << len << "\n";
#endif
    if (start == len) {
      break;
    }
    for (size_t idx = start; idx < len; idx++) {
      Tpair ePair = ListPair[idx]; // Copy is needed since ListPair is extended
#ifdef DEBUG_MATRIX_GROUP
      size_t jdx = 0;
#endif
      for (auto &eMatrGen : ListMatrGen) {
#ifdef DEBUG_MATRIX_GROUP
        os << "MAT_GRP: DirectSpaceOrbit_Stabilizer, jdx=" << jdx << " / " << ListMatrGen.size() << "\n";
        jdx += 1;
#endif
        MyMatrix<T> eSpaceImg = ePair.first * eMatrGen;
        if (f_terminate(eSpaceImg)) {
          return {};
        }
        //
        MyMatrix<T> eSpaceMod = Concatenate(eSpaceImg, ModSpace);
        RecSolutionIntMat<T> eCan(eSpaceMod);
        auto need_insert = [&]() -> bool {
          for (auto &fPair : ListPair) {
            if (eCan.is_containing_m(fPair.first)) {
              return false;
            }
          }
          return true;
        };
        if (need_insert()) {
          MyMatrix<T> eReprImg = ePair.second * eMatrGen;
          Tpair ePairNew{std::move(eSpaceImg), std::move(eReprImg)};
          ListPair.emplace_back(std::move(ePairNew));
        }
      }
    }
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: DirectSpaceOrbit_Stabilizer, end of loop |ListPair|=" << ListPair.size() << "\n";
#endif
    start = len;
  }
  //
  // Orbit is fine, now computing the stabilizer by using the Schreier lemma.
  //
  std::unordered_set<MyMatrix<T>> SetGen;
  size_t nPair = ListPair.size();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: DirectSpaceOrbit_Stabilizer, nPair=" << nPair << "\n";
#endif
  for (size_t iPair = 0; iPair < nPair; iPair++) {
    Tpair const &ePair = ListPair[iPair];
    for (auto &eMatrGen : ListMatrGen) {
      MyMatrix<T> eSpaceImg = ePair.first * eMatrGen;
      MyMatrix<T> eSpaceImgMod = Concatenate(eSpaceImg, ModSpace);
      RecSolutionIntMat<T> eCan(eSpaceImgMod);
      auto f_insert = [&]() -> void {
        for (size_t jPair = 0; jPair < nPair; jPair++) {
          if (eCan.is_containing_m(ListPair[jPair].first)) {
            MyMatrix<T> eGenMatr_new =
                ePair.second * eMatrGen * Inverse(ListPair[jPair].second);
            if (!IsIdentity(eGenMatr_new)) {
              SetGen.insert(eGenMatr_new);
            }
          }
        }
      };
      f_insert();
    }
  }
  std::vector<MyMatrix<T>> ListGen(SetGen.begin(), SetGen.end());
#ifdef TRACK_INFO_MATRIX_GROUP
  write_matrix_group(ListGen, "DirectSpaceOrbit_Stabilizer");
#endif
  return ListGen;
}

template <typename T, typename Tmod, typename Tgroup, typename Thelper>
inline typename std::enable_if<!has_determining_ext<Thelper>::value,
                               std::optional<std::vector<MyVector<Tmod>>>>::type
FindingSmallOrbit([[maybe_unused]] std::vector<MyMatrix<T>> const &ListMatrGen,
                  [[maybe_unused]] MyMatrix<T> const &TheSpace, T const &TheMod,
                  MyVector<T> const &x, [[maybe_unused]] Thelper const &helper,
                  std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: FindingSmallOrbit (!has), start\n";
#endif
  // No determining EXT, hard to find clever ideas.
  MyVector<Tmod> x_mod = ModuloReductionVector<T, Tmod>(x, TheMod);
  Tmod TheMod_mod = UniversalScalarConversion<Tmod, T>(TheMod);
  auto f_prod = [&](MyVector<Tmod> const &eClass,
                    MyMatrix<Tmod> const &eElt) -> MyVector<Tmod> {
    MyVector<Tmod> eVect = eElt.transpose() * eClass;
    return VectorMod(eVect, TheMod_mod);
  };
  std::vector<MyMatrix<Tmod>> ListMatrGenMod =
    ModuloReductionStdVectorMatrix<T, Tmod>(ListMatrGen, TheMod);
  return OrbitComputation(ListMatrGenMod, x_mod, f_prod, os);
}

template <typename T, typename Tmod, typename Tgroup, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value,
                               std::optional<std::vector<MyVector<Tmod>>>>::type
FindingSmallOrbit(std::vector<MyMatrix<T>> const &ListMatrGen,
                  MyMatrix<T> const &TheSpace, T const &TheMod,
                  MyVector<T> const &a, Thelper const &helper,
                  std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: FindingSmallOrbit(has), start\n";
#endif
  using Telt = typename Tgroup::Telt;
  using PreImager = typename Thelper::PreImager;
  int n = TheSpace.rows();
  // The critical number for the computation
  size_t n_limit = 60000;
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: FindingSmallOrbit, n_limit=" << n_limit << "\n";
#endif
  std::vector<MyMatrix<Tmod>> ListMatrGenMod =
        ModuloReductionStdVectorMatrix<T, Tmod>(ListMatrGen, TheMod);
  auto test_adequateness =
      [&](MyVector<T> const &x) -> std::optional<std::vector<MyVector<Tmod>>> {
    MyVector<Tmod> x_mod = ModuloReductionVector<T, Tmod>(x, TheMod);
    size_t pos = 0;
    auto f_terminate = [&]([[maybe_unused]] MyVector<Tmod> const &a) -> bool {
      pos++;
      if (pos == n_limit)
        return true;
      return false;
    };
    Tmod TheMod_mod = UniversalScalarConversion<Tmod, T>(TheMod);
    auto f_prod = [&](MyVector<Tmod> const &eClass,
                      MyMatrix<Tmod> const &eElt) -> MyVector<Tmod> {
      MyVector<Tmod> eVect = eElt.transpose() * eClass;
      return VectorMod(eVect, TheMod_mod);
    };
    return OrbitComputation_limit(ListMatrGenMod, x_mod, f_prod, f_terminate,
                                  os);
  };
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
  MyMatrix<T> TheSpaceMod = Concatenate(TheSpace, ModSpace);
  RecSolutionIntMat<T> eCan(TheSpaceMod);
  auto IsStabilized = [&](MyVector<T> const &V) -> bool {
    for (auto &eMatrGen : ListMatrGen) {
      MyVector<T> Vimg = eMatrGen.transpose() * V;
      bool test = eCan.has_solution_v(Vimg);
      if (!test) {
        return false;
      }
    }
    return true;
  };

  std::optional<std::vector<MyVector<Tmod>>> opt = test_adequateness(a);
  if (opt) {
    return *opt;
  }
  std::vector<Telt> ListPermGen;
  for (auto &eMatrGen : ListMatrGen) {
    Telt ePermGen = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(
        helper, eMatrGen, os);
    ListPermGen.emplace_back(std::move(ePermGen));
  }
  size_t len = helper.EXTfaithful.rows();
  Telt id_perm(len);
  Tgroup GRP(ListPermGen, id_perm);
  std::vector<Tgroup> ListGroup = GRP.GetAscendingChain();
  size_t len_group = ListGroup.size();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: FindingSmallOrbit, len_group=" << len_group
     << " |GRP|=" << GRP.size() << "\n";
#endif
  auto try_basis = [&](MyMatrix<T> const &TheBasis)
      -> std::optional<std::vector<MyVector<Tmod>>> {
    for (int i_row = 0; i_row < TheBasis.rows(); i_row++) {
      MyVector<T> V = GetMatrixRow(TheBasis, i_row);
      if (!IsStabilized(V)) {
        std::optional<std::vector<MyVector<Tmod>>> opt = test_adequateness(V);
        if (opt) {
          std::vector<MyVector<Tmod>> const &ListV = *opt;
#ifdef DEBUG_MATRIX_GROUP
          os << "MAT_GRP: FindingSmallOrbit, |ListV|=" << ListV.size() << "\n";
#endif
          return ListV;
        }
#ifdef DEBUG_MATRIX_GROUP
        os << "MAT_GRP: FindingSmallOrbit, Too large size at i_row=" << i_row
           << "\n";
#endif
      }
    }
    return {};
  };
  PreImager pre_imager = helper.constant_pre_imager();
  for (size_t iGroup = 0; iGroup < len_group; iGroup++) {
    size_t jGroup = len_group - 1 - iGroup;
    Tgroup const &fGRP = ListGroup[jGroup];
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: FindingSmallOrbit, iGroup=" << iGroup
       << " |fGRP|=" << fGRP.size() << "\n";
#endif
    std::vector<MyMatrix<T>> LMatr;
    for (auto &eGen : fGRP.SmallGeneratingSet()) {
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: FindingSmallOrbit(has)), before pre_image_elt for eGen\n";
#endif
      MyMatrix<T> eMat = pre_imager.pre_image_elt(eGen);
      LMatr.emplace_back(std::move(eMat));
    }
    MyMatrix<T> InvBasis = ComputeBasisInvariantSpace(LMatr, TheSpace, TheMod);
    std::optional<std::vector<MyVector<Tmod>>> opt = try_basis(InvBasis);
    if (opt) {
      return *opt;
    }
  }
  bool AllowLargeOrbit = true;
  if (AllowLargeOrbit) {
    n_limit = std::numeric_limits<size_t>::max();
    MyMatrix<T> IdMat = IdentityMat<T>(n);
    std::optional<std::vector<MyVector<Tmod>>> opt = try_basis(IdMat);
    if (opt) {
      return *opt;
    }
  }
  return {};
}

/*
  Does not compile for a mysterious reason.

template <typename T, typename Tmod, typename Tgroup, typename Thelper>
std::optional<std::vector<MyVector<Tmod>>>
FindingSmallOrbitSort(std::vector<MyMatrix<T>> const &ListMatrGen,
                      MyMatrix<T> const &TheSpace, T const &TheMod,
                      MyVector<T> const &a, Thelper const &helper,
                      std::ostream &os) {
  std::optional<std::vector<MyVector<Tmod>>> opt = FindingSmallOrbit(ListMatrGen, TheSpace, TheMod, a, helper, os);
  if (!opt) {
    return {};
  }
  std::vector<MyVector<Tmod>> V = *opt;
  std::sort(V.begin(), V.end());
  return V;
}
*/



// The space must be defining a finite index subgroup of T^n
template <typename T, typename Tmod, typename Tgroup, typename Thelper,
          typename Fstab>
std::vector<MyMatrix<T>>
LinearSpace_ModStabilizer_Tmod(std::vector<MyMatrix<T>> const &ListMatr,
                               Thelper const &helper,
                               MyMatrix<T> const &TheSpace, T const &TheMod,
                               Fstab f_stab, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n = helper.n;
#ifdef DEBUG_MATRIX_GROUP
  T TotSize(1);
  for (int i = 0; i < n; i++)
    TotSize *= TheMod;
  os << "MAT_GRP: LinearSpace_ModStabilizer_Tmod, TheMod=" << TheMod
     << "  n=" << n << " TotSize=" << TotSize << "\n";
#endif
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
  MyMatrix<T> TheSpaceMod = Concatenate(TheSpace, ModSpace);
  // This is the part of the enumeration where we have problems.
  // We have too many vectors to consider which sinks the algorithm.
  // The difficulty of the work is that we have to deal with globally
  // invariant sets of vectors. This is difficult to manipulate in
  // practice.
  //
  // Possible strategies:
  // ---Using the subchains of subgroups. This could force finding
  //    some subgroup. The strategy sometimes works. But it is
  //    ultimately not sufficient for resolving the problem.
  //
  // ---Reposition the lattices so as to create an environment
  //    where we can apply Plesken-Souvignier algorithm in order to
  //    reduce the computational size.
  //    It requires coding Intersection of subgroups, but that seems
  //    feasible.
  //
  // ---Directly computing the orbit of sublattice does not appear
  //    feasible since the equality test of sublattice is very high
  //    and the orbit itself is very large.
  //
  // ---The dream would be to identify the critical orbits in order
  //    to lift the computation. But that seems out of range.
  //
  // ---The subspaces themselves correspond to rather small sets.
  //    Could this be used? It does not seem so since the thing that
  //    is bad is the size of the group itself.
  //
  // ---Could we have some combination of strategies like in good old
  //    time of GAP and the Delaunay polytopes? Yes, but we could do better
  //    with each strategy registering prograssive improvement on the
  //    problem until completely solved.
  //
  // For the iterative improvement strategy, what we would need is some
  // encoding of the partial solution.
  // ---For the stabilizer, we have two things:
  //    ---Problem statement: helper + TheSpace
  //    ---Partial oversolution: ListMatr
  // ---For the equivalence, we have two things
  //    ---TheMod, helper1, helper2, TheSpace2
  //    ---ListMatr1, TheSpace1, MatrEquiv
  // We could look at the quotient. (Z_d)^n / TheSpace and look for point
  // stabilizers Maybe we can translate to classes easily and
  RecSolutionIntMat<T> eCan(TheSpaceMod);
  auto IsNotStabilizing = [&](std::vector<MyMatrix<T>> const &ListMatrInp)
      -> std::optional<MyVector<T>> {
    for (auto &eGen : ListMatrInp) {
      MyMatrix<T> TheSpace_img = TheSpace * eGen;
      for (int i = 0; i < n; i++) {
        MyVector<T> eVectG = GetMatrixRow(TheSpace_img, i);
        bool test = eCan.has_solution_v(eVectG);
        if (!test) {
          return GetMatrixRow(TheSpace, i);
        }
      }
    }
    return {};
  };
  std::vector<MyMatrix<T>> ListMatrRet = ListMatr;
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinearSpace_ModStabilizer_Tmod(A), comp(ListMatrRet)=" << compute_complexity_listmat(ListMatrRet) << "\n";
#endif
  ListMatrRet = ExhaustiveReductionComplexityGroupMatrix(ListMatrRet, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinearSpace_ModStabilizer_Tmod(B), comp(ListMatrRet)=" << compute_complexity_listmat(ListMatrRet) << "\n";
#endif
  while (true) {
#ifdef TIMINGS_MATRIX_GROUP
    MicrosecondTime time;
#endif
    std::optional<MyVector<T>> opt = IsNotStabilizing(ListMatrRet);
#ifdef TIMINGS_MATRIX_GROUP
    os << "|MAT_GRP: IsNotStabilizing|=" << time << "\n";
#endif
    if (!opt) {
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModStabilizer_Tmod, Exiting the loop\n";
#endif
      break;
    }
    const MyVector<T> &V = *opt;
    std::optional<std::vector<MyVector<Tmod>>> opt_fso =
        FindingSmallOrbit<T, Tmod, Tgroup, Thelper>(
            ListMatrRet, TheSpace, TheMod, V, helper, os);
#ifdef TIMINGS_MATRIX_GROUP
    os << "|MAT_GRP: FindingSmallOrbit|=" << time << "\n";
#endif
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!opt_fso) {
      std::cerr << "Failed to find some entry\n";
      throw TerminalException{1};
    }
#endif
    std::vector<MyVector<Tmod>> const &O = *opt_fso;
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: LinearSpace_ModStabilizer_Tmod, |O|=" << O.size() << "\n";
#endif
    Telt ePermS = Telt(SortingPerm<MyVector<Tmod>, Tidx>(O));
    std::function<Telt(MyMatrix<T> const&)> f_get_perm=[&](MyMatrix<T> const& eGen) -> Telt {
      return get_permutation_from_orbit(eGen, O, TheMod, ePermS);
    };
    int nbRow = helper.nbRow();
    Face eFace_pre = GetFace<T, Tmod>(O, TheSpaceMod);
    PartitionReduction<T, Telt> pr(ListMatrRet, f_get_perm, eFace_pre, os);
    Face eFace = TranslateFace(nbRow, pr.face);
    std::vector<Telt> ListPermGens =
      MatrixIntegral_GeneratePermutationGroupA<T, Telt, Thelper, decltype(f_get_perm)>(ListMatrRet, helper, pr.f_get_perm, os);
    Tidx siz_act = eFace.size();
    Tgroup GRPwork(ListPermGens, siz_act);
    //
    // As it turns out, the eFace tend to be disjoint accross different embeddings.
    // That does not seem to be guaranteed by theoretical reasons. Just something
    // common.
    //
    // What we can do is collect the intersections and from that decompose the space
    // as a disjoint union
    //
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: LinearSpace_ModStabilizer_Tmod TheMod=" << TheMod
       << " |O|=" << O.size() << " |GRPwork|=" << GRPwork.size()
       << " |eFace|=" << eFace.count() << "\n";
#endif
    ListMatrRet = f_stab(ListPermGens, GRPwork, eFace, pr.f_get_perm, ListMatrRet);
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: LinearSpace_ModStabilizer_Tmod(C), comp(ListMatrRet)=" << compute_complexity_listmat(ListMatrRet) << "\n";
#endif
    ListMatrRet = ExhaustiveReductionComplexityGroupMatrix(ListMatrRet, os);
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: LinearSpace_ModStabilizer_Tmod(D), comp(ListMatrRet)=" << compute_complexity_listmat(ListMatrRet) << "\n";
#endif
  }
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinearSpace_ModStabilizer_Tmod |ListMatrRet|=" << ListMatrRet.size() << "\n";
#endif
  return ListMatrRet;
}

template <typename T, typename Tgroup, typename Thelper, typename Fstab>
std::vector<MyMatrix<T>>
LinearSpace_ModStabilizer(std::vector<MyMatrix<T>> const &ListMatr,
                          Thelper const &helper, MyMatrix<T> const &TheSpace,
                          T const &TheMod, Fstab f_stab, std::ostream &os) {
  T max_size = (TheMod - 1) * (TheMod - 1) * TheSpace.rows();
  if (max_size < T(std::numeric_limits<uint8_t>::max())) {
    return LinearSpace_ModStabilizer_Tmod<T, uint8_t, Tgroup, Thelper>(
        ListMatr, helper, TheSpace, TheMod, f_stab, os);
  }
  if (max_size < T(std::numeric_limits<uint16_t>::max())) {
    return LinearSpace_ModStabilizer_Tmod<T, uint16_t, Tgroup, Thelper>(
        ListMatr, helper, TheSpace, TheMod, f_stab, os);
  }
  if (max_size < T(std::numeric_limits<uint32_t>::max())) {
    return LinearSpace_ModStabilizer_Tmod<T, uint32_t, Tgroup, Thelper>(
        ListMatr, helper, TheSpace, TheMod, f_stab, os);
  }
  std::cerr << "Failed to find a matching arithmetic type. Quite unlikely "
               "objectively\n";
  throw TerminalException{1};
}

template <typename T, typename Tgroup, typename Thelper, typename Fstab>
std::vector<MyMatrix<T>> LinearSpace_StabilizerGen_Kernel(
    std::vector<MyMatrix<T>> const &ListGen, Thelper const &helper,
    MyMatrix<T> const &TheSpace, Fstab f_stab, std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: det(TheSpace)=" << DeterminantMat(TheSpace) << "\n";
#endif
  RecSolutionIntMat<T> eCan(TheSpace);
  auto IsStabilizing = [&](std::vector<MyMatrix<T>> const &ListGen) -> bool {
    for (auto &eGen : ListGen) {
      MyMatrix<T> TheSpace_img = TheSpace * eGen;
      if (!eCan.is_containing_m(TheSpace_img)) {
        return false;
      }
    }
    return true;
  };
  if (IsStabilizing(ListGen)) {
    return ListGen;
  }
  T LFact = LinearSpace_GetDivisor(TheSpace);
  std::vector<T> eList = FactorsInt(LFact);
  int siz = eList.size();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinearSpace_StabilizerGen_Kernel LFact=" << LFact
     << " siz=" << siz << "\n";
#endif
  std::vector<MyMatrix<T>> ListGenRet = ListGen;
  for (int i = 1; i <= siz; i++) {
    T TheMod(1);
    for (int j = 0; j < i; j++) {
      TheMod *= eList[j];
    }
    ListGenRet = LinearSpace_ModStabilizer<T, Tgroup, Thelper, Fstab>(
        ListGenRet, helper, TheSpace, TheMod, f_stab, os);
    if (IsStabilizing(ListGenRet)) {
      return ListGenRet;
    }
  }
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsStabilizing(ListGenRet)) {
    std::cerr << "MAT_GRP: Error in LinearSpace_Stabilizer_Kernel\n";
    throw TerminalException{1};
  }
#endif
  return ListGenRet;
}

template <typename T, typename Tgroup, typename Thelper>
RetMI_S<T, Tgroup> LinearSpace_Stabilizer_KernelRing(std::vector<MyMatrix<T>> const &ListGen,
                              Thelper const &helper,
                              MyMatrix<T> const &TheSpace, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  TintGroup total_index(1);
  auto f_stab = [&](std::vector<Telt> const &ListPermGens, Tgroup const &GRP,
                    Face const &eFace,
                    [[maybe_unused]] std::function<Telt(MyMatrix<T> const&)> f_get_perm,
                    std::vector<MyMatrix<T>> const& ListMatr) -> std::vector<MyMatrix<T>> {
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: f_stab(LinearSpace_Stabilizer_Kernel), before MatrixIntegral_Stabilizer\n";
#endif
    RetMI_S<T,Tgroup> ret = MatrixIntegral_Stabilizer<T, Tgroup, Thelper>(ListPermGens, ListMatr, f_get_perm, GRP, helper, eFace, os);
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: f_stab(LinearSpace_Stabilizer_Kernel), after MatrixIntegral_Stabilizer\n";
#endif
    total_index *= ret.index;
    return ret.LGen;
  };
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Before LinearSpace_StabilizerGen_Kernel\n";
#endif
  std::vector<MyMatrix<T>> ListGenRet1 =
      LinearSpace_StabilizerGen_Kernel<T, Tgroup, Thelper, decltype(f_stab)>(
          ListGen, helper, TheSpace, f_stab, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: After LinearSpace_StabilizerGen_Kernel\n";
#endif
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinearSpace_Stabilizer_Kernel, comp(ListGenRet1)=" << compute_complexity_listmat(ListGenRet1) << "\n";
#endif
  std::vector<MyMatrix<T>> ListGenRet2 = ExhaustiveReductionComplexityGroupMatrix<T>(ListGenRet1, os);
#ifdef SANITY_CHECK_MATRIX_GROUP
  CheckGroupEquality<T,Tgroup>(ListGenRet1, ListGenRet2, os);
#endif
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinearSpace_Stabilizer_Kernel, comp(ListGenRet2)=" << compute_complexity_listmat(ListGenRet2) << "\n";
#endif
  return {total_index, ListGenRet2};
}


template <typename T, typename Tgroup, typename Thelper>
RetMI_S<T, Tgroup> LinearSpace_Stabilizer_Kernel(std::vector<MyMatrix<T>> const &ListGen,
                              Thelper const &helper,
                              MyMatrix<T> const &TheSpace, std::ostream &os) {
  using Tint = typename Thelper::Tint;
  using ThelperInt = typename Thelper::ThelperInt;
  std::vector<MyMatrix<Tint>> ListGen_int = UniversalStdVectorMatrixConversion<Tint,T>(ListGen);
  MyMatrix<Tint> TheSpace_int = UniversalMatrixConversion<Tint,T>(TheSpace);
  ThelperInt helper_int = ToInteger(helper);
  RetMI_S<Tint, Tgroup> ret = LinearSpace_Stabilizer_KernelRing<Tint,Tgroup,ThelperInt>(ListGen_int, helper_int, TheSpace_int, os);
  std::vector<MyMatrix<T>> LGen = UniversalStdVectorMatrixConversion<T,Tint>(ret.LGen);
  return {ret.index, LGen};
}


template <typename T> struct Stab_RightCoset {
  std::vector<MyMatrix<T>> list_gen;
  CosetDescription<T> coset_desc;
};

template <typename T, typename Tgroup, typename Thelper>
Stab_RightCoset<T> LinearSpace_Stabilizer_RightCoset_KernelRing(
    std::vector<MyMatrix<T>> const &l_gens, Thelper const &helper,
    MyMatrix<T> const &TheSpace, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  int n = helper.n;
  CosetDescription<T> coset(n);
  auto f_stab = [&](std::vector<Telt> const &ListPermGens, Tgroup const &GRP,
                    Face const &eFace,
                    [[maybe_unused]] std::function<Telt(MyMatrix<T> const&)> f_get_perm,
                    std::vector<MyMatrix<T>> const& ListMatr) -> std::vector<MyMatrix<T>> {
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: f_stab(LinearSpace_Stabilizer_RightCoset_Kernel), beginning\n";
#endif
    std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>> pair =
    MatrixIntegral_Stabilizer_RightCoset<T, Tgroup, Thelper>(ListPermGens, ListMatr, GRP, helper, eFace, os);
    coset.insert(pair.second);
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: f_stab(LinearSpace_Stabilizer_RightCoset_Kernel), after MatrixIntegral_Stabilizer_RightCoset\n";
#endif
    return pair.first;
  };
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Before LinearSpace_StabilizerGen_Kernel\n";
#endif
  std::vector<MyMatrix<T>> l_gens_ret =
      LinearSpace_StabilizerGen_Kernel<T, Tgroup, Thelper, decltype(f_stab)>(
          l_gens, helper, TheSpace, f_stab, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: After LinearSpace_StabilizerGen_Kernel\n";
#endif
  return {std::move(l_gens_ret), coset};
}


template <typename T, typename Tgroup, typename Thelper>
Stab_RightCoset<T> LinearSpace_Stabilizer_RightCoset_Kernel(
    std::vector<MyMatrix<T>> const &l_gens, Thelper const &helper,
    MyMatrix<T> const &TheSpace, std::ostream &os) {
  using Tint = typename Thelper::Tint;
  using ThelperInt = typename Thelper::ThelperInt;
  std::vector<MyMatrix<Tint>> l_gens_int = UniversalStdVectorMatrixConversion<Tint,T>(l_gens);
  ThelperInt helper_int = ToInteger(helper);
  MyMatrix<Tint> TheSpace_int = UniversalMatrixConversion<Tint,T>(TheSpace);
  Stab_RightCoset<Tint> src =
    LinearSpace_Stabilizer_RightCoset_KernelRing<Tint,Tgroup,ThelperInt>(l_gens_int, helper_int, TheSpace_int, os);
  std::vector<MyMatrix<T>> list_gen = UniversalStdVectorMatrixConversion<T,Tint>(src.list_gen);
  int n = src.coset_desc.n;
  CosetDescription<T> coset_desc(n);
  for (auto & l_coset: src.coset_desc.ListListCoset) {
    std::vector<MyMatrix<T>> ListCoset = UniversalStdVectorMatrixConversion<T,Tint>(l_coset);
    coset_desc.insert(ListCoset);
  }
  return {list_gen, coset_desc};
}

template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value,
                               std::vector<MyMatrix<T>>>::type
MatrixIntegral_PreImageSubgroup(std::vector<typename Tgroup::Telt> const &ListPermGens,
                                std::vector<MyMatrix<T>> const& ListMatr,
                                Tgroup const &eGRP, Thelper const &helper,
                                [[maybe_unused]] std::function<typename Tgroup::Telt(MyMatrix<T> const&)> f_get_perm,
                                [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin MatrixIntegral_PreImageSubgroup(has) |ListPermGens|=" << ListPermGens.size() << " |gen(eGRP)|=" << eGRP.GeneratorsOfGroup().size() << "\n";
#endif
  using Telt = typename Tgroup::Telt;
  using PreImager = typename Thelper::PreImager;
  PreImager pre_imager = helper.pre_imager(ListMatr, ListPermGens);
  std::vector<Telt> LGenSmall = eGRP.SmallGeneratingSet();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin MatrixIntegral_PreImageSubgroup(has) |LGenSmall|=" << LGenSmall.size() << "\n";
#endif
  std::vector<MyMatrix<T>> ListMatrGen1;
  for (auto &eGen : LGenSmall) {
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: MatrixIntegral_PreImageSubgroup(has), before pre_image_elt for eGen\n";
#endif
    MyMatrix<T> eMatr = pre_imager.pre_image_elt(eGen);
    ListMatrGen1.emplace_back(std::move(eMatr));
  }
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: |ListMatrGen1|=" << ListMatrGen1.size() << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_PreImageSubgroup(has), comp(ListMatrGen1)=" << compute_complexity_listmat(ListMatrGen1) << "\n";
#endif
  std::vector<MyMatrix<T>> ListMatrGen2 = ExhaustiveReductionComplexityGroupMatrix<T>(ListMatrGen1, os);
#ifdef SANITY_CHECK_MATRIX_GROUP
  CheckGroupEquality<T,Tgroup>(ListMatrGen1, ListMatrGen2, os);
#endif
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_PreImageSubgroup(has), comp(ListMatrGen2)=" << compute_complexity_listmat(ListMatrGen2) << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: |ListMatrGen2|=" << ListMatrGen2.size() << "\n";
#endif
  return ListMatrGen2;
}

template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<!has_determining_ext<Thelper>::value,
                               std::vector<MyMatrix<T>>>::type
MatrixIntegral_PreImageSubgroup(std::vector<typename Tgroup::Telt> const &ListPermGens,
                                std::vector<MyMatrix<T>> const& ListMatrGens,
                                Tgroup const &eGRP, Thelper const &helper,
                                std::function<typename Tgroup::Telt(MyMatrix<T> const&)> f_get_perm,
                                std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin MatrixIntegral_PreImageSubgroup(!has) |ListPermGens|=" << ListPermGens.size() << " |gen(eGRP)|=" << eGRP.GeneratorsOfGroup().size() << "\n";
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n_act = ListPermGens[0].size();
  Tgroup GRP_build(ListPermGens, n_act);
  os << "MAT_GRP: Begin MatrixIntegral_PreImageSubgroup(!has) |GRP_build|=" << GRP_build.size() << " |eGRP|=" << eGRP.size() << "\n";
  bool test = GRP_build.IsSubgroup(eGRP);
  os << "MAT_GRP: Begin MatrixIntegral_PreImageSubgroup(!has) IsSubgroup=" << test << "\n";
  if (false) {
    WriteGroupFile("GRP_build", GRP_build);
    WriteGroupFile("eGRP", eGRP);
    std::cerr << "Now debugging from here\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS_MATRIX_GROUP
  MicrosecondTime time;
#endif
  MyMatrix<T> id_matr = IdentityMat<T>(helper.n);
  std::vector<MyMatrix<T>> ListGen =
    PreImageSubgroup<T,Tgroup>(ListMatrGens, ListPermGens, f_get_perm, id_matr, eGRP, os);
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: MatrixIntegral_PreImageSubgroup(!has), PreImageSubgroup|=" << time << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: After PreImageSubgroup comp(ListGen1)=" << compute_complexity_listmat(ListGen) << "\n";
#endif
  return ListGen;
}

/*
  The dcs is the representative.
  The list stab_gens is a set of generators "x"
  of a group of transformations V
  such that
  U dcs x = U.dcs
 */
template<typename T>
struct DoubleCosetEntry {
  MyMatrix<T> cos;
  std::vector<MyMatrix<T>> stab_gens;
};

template <typename T, typename Tgroup, typename Thelper>
void TestPreImageSubgroup(Thelper const &helper,
                          std::vector<typename Tgroup::Telt> const& ListPermGen,
                          std::vector<MyMatrix<T>> const& ListMatrGen,
                          std::function<typename Tgroup::Telt(MyMatrix<T> const&)> const& f_get_perm,
                          std::vector<MyMatrix<T>> const& MatrPreImage,
                          Tgroup const& OrigGRP, std::string const& context,
                          std::ostream& os) {
  using PreImager = typename Thelper::PreImager;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  PreImager pre_imager = helper.pre_imager(ListMatrGen, ListPermGen);
  if (MatrPreImage.size() == 0) {
    if (OrigGRP.size() != 1) {
      std::cerr << "MAT_GRP: Error, in context=" << context << "\n";
      std::cerr << "MAT_GRP: If MatrPreImage is empty, then necessarily OrigGRP has to be trivial\n";
      std::cerr << "MAT_GRP: Though that may not be sufficient in itself if the mapping is not an isomorphism\n";
      throw TerminalException{1};
    }
  }
  auto f_map_matr=[&](MyMatrix<T> const& eMatr) -> Telt {
    return MatrixIntegral_MapMatrix<T,Telt,Thelper,decltype(f_get_perm)>(helper, f_get_perm, eMatr, os);
  };
  size_t n_gen = ListPermGen.size();
  for (size_t i_gen=0; i_gen<n_gen; i_gen++) {
    Telt ePerm = f_map_matr(ListMatrGen[i_gen]);
    if (ePerm != ListPermGen[i_gen]) {
      std::cerr << "MAT_GRP: Error detected in TestPreImageSubgroup, context=" << context << "\n";
      std::cerr << "MAT_GRP: i_gen=" << i_gen << "\n";
      std::cerr << "MAT_GRP: Inconsistency with ListMatrGen / ListPermGen\n";
      throw TerminalException{1};
    }
  }
  Tidx n_act = OrigGRP.n_act();
  std::vector<Telt> LGen;
  for (auto & eMatrGen : MatrPreImage) {
    Telt ePermGen = f_map_matr(eMatrGen);
    if (ePermGen.size() != n_act) {
      std::cerr << "MAT_GRP: Error detected in TestPreImageSubgroup, context=" << context << "\n";
      std::cerr << "MAT_GRP: The obtained permutation is not of the right size\n";
      throw TerminalException{1};
    }
    LGen.push_back(ePermGen);
  }
  std::vector<Telt> LGenB;
  size_t n_orig_gen = 0;
  size_t n_orig_gen_error = 0;
  std::vector<std::pair<Telt, Telt>> l_pair_err;
  for (auto & eGen: OrigGRP.SmallGeneratingSet()) {
    n_orig_gen += 1;
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: TestPreImageSubgroup, before pre_image_elt for eGen\n";
#endif
    MyMatrix<T> eMatr = pre_imager.pre_image_elt(eGen);
    Telt eGenB = f_map_matr(eMatr);
    if (eGen != eGenB) {
      n_orig_gen_error += 1;
    }
    LGenB.push_back(eGenB);
  }
  Tgroup GRPfull(ListPermGen, n_act);
  Tgroup GRPimg(LGen, n_act);
  Tgroup GRPimgB(LGenB, n_act);
  if (GRPimg != OrigGRP) {
    Tgroup IntGRP = GRPimg.Intersection(OrigGRP);
    std::cerr << "MAT_GRP: Error detected in TestPreImageSubgroup, context=" << context << "\n";
    std::cerr << "MAT_GRP: n_orig_gen=" << n_orig_gen << " n_orig_gen_error=" << n_orig_gen_error << "\n";
    std::cerr << "MAT_GRP: |GRPfull|=" << GRPfull.size() << "\n";
    std::cerr << "MAT_GRP: |GRPimg|=" << GRPimg.size() << " |OrigGRP|=" << OrigGRP.size() << " |IntGRP|=" << IntGRP.size() << "\n";
    std::cerr << "MAT_GRP: The image of the PreImage is not equal to the original group\n";
    std::cerr << "MAT_GRP: This failing one basic check of correctness (this is not a complete\n";
    std::cerr << "MAT_GRP: check if the morphism is not an isomorphism)\n";
    std::cerr << "MAT_GRP: |GRPimgB|=" << GRPimgB.size() << "\n";
    if (GRPimgB == OrigGRP) {
      std::cerr << "MAT_GRP: GRPimgB is EQUAL to OrigGRP\n";
    } else {
      std::cerr << "MAT_GRP: GRPimgB is NOT EQUAL to OrigGRP\n";
    }
    if (GRPimgB == GRPimg) {
      std::cerr << "MAT_GRP: GRPimgB is EQUAL to GRPimg\n";
    } else {
      std::cerr << "MAT_GRP: GRPimgB is NOT EQUAL to GRPimg\n";
    }
    std::string FileGRP = "AllGRP";
    std::ofstream os_grp(FileGRP);
    os_grp << "return rec(GRPfull:=" << GRPfull.GapString() << ",\n"
           << " OrigGRP:=" << OrigGRP.GapString() << ",\n"
           << " GRPimg:=" << GRPimg.GapString() << ",\n"
           << " GRPimgB:=" << GRPimgB.GapString() << ");\n";
    throw TerminalException{1};
  }
}


template<typename T, typename Telt>
struct ResultSimplificationDoubleCosets {
  MyMatrix<T> cos_matr;
  Telt cos_perm;
  Telt elt_u;
  Telt elt_v;
};


template<typename T, typename Tgroup, typename Thelper>
ResultSimplificationDoubleCosets<T, typename Tgroup::Telt> IterativeSimplificationDoubleCoset(Thelper const& helper,
                                                                                              typename Thelper::PreImager pre_imager,
                                                                                              Tgroup GRP_U, Tgroup GRP_V,
                                                                                              std::vector<MyMatrix<T>> const& ListMatr_U,
                                                                                              std::vector<MyMatrix<T>> const& ListMatr_V,
                                                                                              std::function<typename Tgroup::Telt(MyMatrix<T> const&)> const& f_get_perm,
                                                                                              typename Tgroup::Telt cos_perm,
                                                                                              std::ostream& os) {
  using Telt = typename Tgroup::Telt;
  struct IntermediateState {
    MyMatrix<T> cos_matr;
    Telt cos_perm;
    DoubleCosetSimplification<T> udv;
  };
  using Tresult = std::pair<IntermediateState, T>;
  int n = helper.n;
  T absolute_minimum = UniversalScalarConversion<T,int>(n);
  auto f_norm=[&](MyMatrix<T> const& H) -> T {
    T norm(0);
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        norm += T_abs(H(i,j));
      }
    }
    return norm;
  };
  Telt id = GRP_U.get_identity();
  auto get_result_reduction=[&](Telt const& x_u, Telt const& x_v) -> Tresult {
    Telt cos_perm_cand = x_u * cos_perm * x_v;
    MyMatrix<T> cos_matr_cand = pre_imager.pre_image_elt(cos_perm_cand);
    size_t max_iter = 1000;
    DoubleCosetSimplification<T> udv = ExhaustiveMatrixDoubleCosetSimplifications(cos_matr_cand, ListMatr_U, ListMatr_V, max_iter, os);
    T norm_cand = f_norm(udv.d_cos_red);
    IntermediateState is{cos_matr_cand, cos_perm_cand, udv};
    return {is, norm_cand};
  };
  Telt elt_u = id;
  Telt elt_v = id;
  Tresult res_work = get_result_reduction(elt_u, elt_v);
  auto get_final=[&]() -> ResultSimplificationDoubleCosets<T,Telt> {
    Telt udv_u_img = f_get_perm(res_work.first.udv.u_red);
    Telt udv_v_img = f_get_perm(res_work.first.udv.v_red);
    Telt elt_u_ret = udv_u_img * elt_u;
    Telt elt_v_ret = elt_v * udv_v_img;
    Telt cos_perm_ret = elt_u_ret * cos_perm * elt_v_ret;
    MyMatrix<T> const& cos_matr = res_work.first.udv.d_cos_red;
#ifdef SANITY_CHECK_DOUBLE_COSET_ENUM
    if (f_get_perm(cos_matr) != cos_perm_ret) {
      std::cerr << "MAT_GRP: cos_matr is not as we expect\n";
      throw TerminalException{1};
    }
#endif
    return {cos_matr, cos_perm_ret, elt_u_ret, elt_v_ret};
  };
  if (res_work.second == absolute_minimum) {
    return get_final();
  }
  int n_iter = 100;
  while(true) {
    size_t n_improv = 0;
    for (int i=0; i<n_iter; i++) {
      os << "MAT_GRP: i=" << i << " / " << n_iter << " res_work.second=" << res_work.second << "\n";
      Telt rand_u = GRP_U.rand();
      Telt rand_v = GRP_V.rand();
      Telt elt_u_cand = rand_u * elt_u;
      Telt elt_v_cand = elt_v * rand_v;
      Tresult res_cand = get_result_reduction(elt_u_cand, elt_v_cand);
      if (res_cand.second < res_work.second) {
        n_improv += 1;
        res_work = res_cand;
        elt_u = elt_u_cand;
        elt_v = elt_v_cand;
        os << "MAT_GRP:   Now norm_work=" << res_work.second << " cos_matr_work=\n";
        WriteMatrix(os, res_work.first.udv.d_cos_red);
        if (res_work.second == absolute_minimum) {
          return get_final();
        }
      }
    }
    if (n_improv == 0) {
      return get_final();
    }
  }
}

/*
  If there is a single double coset, then we can set it to identity.
 */
template<typename Tgroup>
std::vector<typename Tgroup::DccEntry> simplify_span_de(Tgroup const& grp, std::vector<typename Tgroup::DccEntry> const& span_de) {
  using DccEntry = typename Tgroup::DccEntry;
  using Telt = typename Tgroup::Telt;
  if (span_de.size() > 1) {
    return span_de;
  }
  Telt id = grp.get_identity();
  DccEntry const& ent = span_de[0];
  Telt cos = ent.cos;
  Telt cos_inv = Inverse(cos);
  std::vector<Telt> stab_gens;
  for (auto & e_gen: ent.stab_gens) {
    Telt f_gen = cos * e_gen * cos_inv;
    stab_gens.push_back(f_gen);
  }
  DccEntry f_ent{id, stab_gens};
  return {f_ent};
}


/*
  Computes the Double cosets between the stabilizer of the subspace
  and the group V in argument.
  --
  PB1: In contrast to the classical algorithm for finite groups,
  we cannot know when we are at the last step. Therefore, we have to
  compute the stabilizers all the time.
  --
  PB2: The orbit being generated can be on the big side of things.
  Therefore, we need to use the double cosets from the finite
  group algorithm. Likely, we also want to have the corresponding
  stabilizing subgroup.
  ---
  PB3: Computing in permutation groups is a winning strategy. This
  works fine for the has_determining_ext = true case.
  We need to have a single function that does the job. That is in
  order to avoid recomputing the stabilizer of eFace.
  For the has_determining_ext = false case, I need to code that
  in the permutalib code. That will be a single loop. But that is
  fine because it works that way for the cosets. Will see if further
  work is needed if failing.
 */
template <typename T, typename Tgroup, typename Thelper>
std::pair<std::vector<MyMatrix<T>>, std::vector<DoubleCosetEntry<T>>>
LinearSpace_Stabilizer_DoubleCosetStabilizer_KernelRing(
    std::vector<MyMatrix<T>> const l_gens, Thelper const &helper,
    MyMatrix<T> const &TheSpace, std::vector<MyMatrix<T>> const& Vmatr_gens, std::ostream &os) {
  using PreImager = typename Thelper::PreImager;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using DoubleCosetComputer = typename Tgroup::DoubleCosetComputer;
  using DccEntry = typename Tgroup::DccEntry;
  int n = helper.n;
  std::vector<DoubleCosetEntry<T>> entries;
  entries.push_back({IdentityMat<T>(n), Vmatr_gens});
  auto f_stab = [&](std::vector<Telt> const &ListPermGens, Tgroup const &GRP,
                    Face const &eFace,
                    std::function<Telt(MyMatrix<T> const&)> f_get_perm,
                    std::vector<MyMatrix<T>> const& ListMatrGens) -> std::vector<MyMatrix<T>> {
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: f_stab(LinearSpace_Stabilizer_DoubleCosetStabilizer_Kernel), beginning |GRP|=" << GRP.size() << " n_act=" << static_cast<size_t>(GRP.n_act()) << "\n";
#endif
    Tgroup eStab_perm = GRP.Stabilizer_OnSets(eFace);
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: f_stab(LinearSpace_Stabilizer_DoubleCosetStabilizer_Kernel), |eStab_perm|=" << eStab_perm.size() << " |ListPermGens|=" << ListPermGens.size() << "\n";
    os << "MAT_GRP: f_stab(LinearSpace_Stabilizer_DoubleCosetStabilizer_Kernel), comp(ListMatrGens)=" << compute_complexity_listmat(ListMatrGens) << "\n";
#endif
#ifdef TRACK_INFO_MATRIX_GROUP
    write_matrix_group(ListMatrGens, "LinearSpace_Stabilizer_DoubleCosetStabilizer_Kernel_f_stab");
#endif
    PreImager pre_imager = helper.pre_imager(ListMatrGens, ListPermGens);
#ifdef DEBUG_DOUBLE_COSET_ENUM
    os << "MAT_GRP: We have pre_imager\n";
#endif
    std::vector<MyMatrix<T>> eStab_matr = MatrixIntegral_PreImageSubgroup<T,Tgroup,Thelper>(ListPermGens, ListMatrGens, eStab_perm, helper, f_get_perm, os);
#ifdef SANITY_CHECK_DOUBLE_COSET_ENUM_DISABLE
    TestPreImageSubgroup(helper, ListPermGens, ListMatrGens, f_get_perm, eStab_matr, eStab_perm, "eStab", os);
#endif
    std::vector<MyMatrix<T>> eStab_matr_tot = Exhaust_get_total_generators(eStab_matr);
    DoubleCosetComputer dcc_v = GRP.double_coset_computer_v(eStab_perm);
#ifdef DEBUG_DOUBLE_COSET_ENUM
    os << "MAT_GRP: We have dcc_v\n";
#endif
    Tidx siz_act = eFace.size();
    std::vector<DoubleCosetEntry<T>> new_entries;
    for (auto &entry: entries) {
      MyMatrix<T> cos_inv = Inverse(entry.cos);
      std::vector<MyMatrix<T>> const& Vmatr = entry.stab_gens;
      std::vector<MyMatrix<T>> Vmatr_conj;
      for (auto & eGen : Vmatr) {
        MyMatrix<T> NewGen = entry.cos * eGen * cos_inv;
        Vmatr_conj.emplace_back(std::move(NewGen));
      }
#ifdef DEBUG_DOUBLE_COSET_ENUM
      os << "MAT_GRP: We have Vmatr_conj\n";
#endif
      std::vector<Telt> Vperm_conj =
        MatrixIntegral_GeneratePermutationGroupA<T, Telt, Thelper, std::function<Telt(MyMatrix<T> const&)>>(Vmatr_conj, helper, f_get_perm, os);
#ifdef DEBUG_DOUBLE_COSET_ENUM
      os << "MAT_GRP: We have Vperm_conj\n";
#endif
      Tgroup Vgroup_conj = Tgroup(Vperm_conj, siz_act);
#ifdef SANITY_CHECK_DOUBLE_COSET_ENUM
      bool test_is_sub = GRP.IsSubgroup(Vgroup_conj);
      if (!test_is_sub) {
        std::cerr << "MAT_GRP: Vgroup_conj should be a subgroup of GRP\n";
        std::cerr << "MAT_GRP: |GRP|=" << GRP.size() << " |eStab_perm|=" << eStab_perm.size() << " |Vgroup_conj|=" << Vgroup_conj.size() << "\n";
        throw TerminalException{1};
      }
#endif
      // The double cosets being computed are of the form U x V.
      // The returned entry U x V has an associated group S.
      // That group satisfies U x S = U x.
      // This is equivalent to U x S x^{-1} = U or W = x S x^{-1} \subset U.
      // The new stabilizer T satisfies W = y T y^{-1} for y the new stabilizer.
      // So, T = y^{-1}x S x^{-1}y therefore we write z = y^{-1}x and zinv = x^{-1}y.
      // ---
      // In our context we write x_new = u x v    and x = u^{-1} x_new v^{-1}
      // That gets us U x_new v^{-1} S = U x_new v^{-1}
      // So, S_new = v^{-1} S v
      std::vector<DccEntry> pre_span_de = dcc_v.double_cosets_and_stabilizers(Vgroup_conj);
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_Stabilizer_DoubleCosetStabilizer_Kernel, |pre_span_de|=" << pre_span_de.size() << "\n";
#endif
      std::vector<DccEntry> span_de = pre_span_de;
      //      std::vector<DccEntry> span_de = simplify_span_de(Vgroup_conj, pre_span_de);
      for (auto & e_de: span_de) {
#ifdef DEBUG_MATRIX_GROUP
        os << "MAT_GRP: LinearSpace_Stabilizer_DoubleCosetStabilizer_Kernel, before pre_image_elt for e_de.cos\n";
#endif
        ResultSimplificationDoubleCosets<T, Telt> result = IterativeSimplificationDoubleCoset<T,Tgroup,Thelper>(helper,
                                                                                                                pre_imager,
                                                                                                                eStab_perm, Vgroup_conj,
                                                                                                                eStab_matr, Vmatr_conj,
                                                                                                                f_get_perm,
                                                                                                                e_de.cos, os);
        Telt const& v = result.elt_v;
        Telt vinv = Inverse(v);
        std::vector<Telt> stab_gens;
        for (auto & e_gen: e_de.stab_gens) {
          Telt f_gen = vinv * e_gen * v;
          stab_gens.push_back(f_gen);
        }
        MyMatrix<T> const& eCos = result.cos_matr;
#ifdef DEBUG_MATRIX_GROUP
        os << "MAT_GRP: LinearSpace_Stabilizer_DoubleCosetStabilizer_Kernel, comp(eCos)=" << compute_complexity_matrix(eCos) << " eCos=\n";
        WriteMatrix(os, eCos);
#endif
        Tgroup Stab_perm(stab_gens, siz_act);
#ifdef SANITY_CHECK_DOUBLE_COSET_ENUM
        if (!Vgroup_conj.IsSubgroup(Stab_perm)) {
          std::cerr << "MAT_GRP: Vgroup_conj should contain Stab_perm as subgroup\n";
          throw TerminalException{1};
        }
#endif
        std::vector<MyMatrix<T>> Stab_matr = MatrixIntegral_PreImageSubgroup<T,Tgroup,Thelper>(Vperm_conj, Vmatr_conj, Stab_perm, helper, f_get_perm, os);
#ifdef SANITY_CHECK_DOUBLE_COSET_ENUM_DISABLE
        TestPreImageSubgroup(helper, Vperm_conj, Vmatr_conj, f_get_perm, Stab_matr, Stab_perm, "Stab_perm", os);
#endif
        std::vector<MyMatrix<T>> Stab_matr_conj;
        for (auto & eGen : Stab_matr) {
          MyMatrix<T> NewGen = cos_inv * eGen * entry.cos;
          Stab_matr_conj.emplace_back(std::move(NewGen));
        }
        std::vector<MyMatrix<T>> Stab_matr_conj_red = ExhaustiveReductionComplexityGroupMatrix(Stab_matr_conj, os);
#ifdef SANITY_CHECK_MATRIX_GROUP
        CheckGroupEquality<T,Tgroup>(Stab_matr_conj_red, Stab_matr_conj, os);
#endif
        MyMatrix<T> new_cos = eCos * entry.cos;
        DoubleCosetEntry<T> new_de{std::move(new_cos), std::move(Stab_matr_conj)};
        new_entries.emplace_back(std::move(new_de));
      }
    }
#ifdef DEBUG_DOUBLE_COSET_ENUM
    os << "MAT_GRP: We found |new_entries|=" << new_entries.size() << "\n";
#endif
    entries = new_entries;
    return eStab_matr;
  };
  std::vector<MyMatrix<T>> l_gens_ret =
      LinearSpace_StabilizerGen_Kernel<T, Tgroup, Thelper, decltype(f_stab)>(
          l_gens, helper, TheSpace, f_stab, os);
  return {std::move(l_gens_ret), std::move(entries)};
}

template <typename T, typename Tgroup, typename Thelper>
std::pair<std::vector<MyMatrix<T>>, std::vector<DoubleCosetEntry<T>>>
LinearSpace_Stabilizer_DoubleCosetStabilizer_Kernel(
    std::vector<MyMatrix<T>> const l_gens, Thelper const &helper,
    MyMatrix<T> const &TheSpace, std::vector<MyMatrix<T>> const& Vmatr_gens, std::ostream &os) {
  using Tint = typename Thelper::Tint;
  using ThelperInt = typename Thelper::ThelperInt;
  std::vector<MyMatrix<Tint>> l_gens_int = UniversalStdVectorMatrixConversion<Tint,T>(l_gens);
  std::vector<MyMatrix<Tint>> Vmatr_gens_int = UniversalStdVectorMatrixConversion<Tint,T>(Vmatr_gens);
  MyMatrix<Tint> TheSpace_int = UniversalMatrixConversion<Tint,T>(TheSpace);
  ThelperInt helper_int = ToInteger(helper);
  std::pair<std::vector<MyMatrix<Tint>>, std::vector<DoubleCosetEntry<Tint>>> pair_int =
    LinearSpace_Stabilizer_DoubleCosetStabilizer_KernelRing<Tint,Tgroup,ThelperInt>(l_gens_int, helper_int, TheSpace_int, Vmatr_gens_int, os);
  std::vector<MyMatrix<T>> l_gens_new = UniversalStdVectorMatrixConversion<T,Tint>(pair_int.first);
  std::vector<DoubleCosetEntry<T>> l_dcs;
  for (auto & dcs: pair_int.second) {
    MyMatrix<T> cos = UniversalMatrixConversion<T,Tint>(dcs.cos);
    std::vector<MyMatrix<T>> stab_gens = UniversalStdVectorMatrixConversion<T,Tint>(dcs.stab_gens);
    DoubleCosetEntry<T> new_dcs{cos, stab_gens};
    l_dcs.push_back(new_dcs);
  }
  return {l_gens_new, l_dcs};
}

template <typename T, typename Tgroup, typename Thelper>
std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>>
LinearSpace_Stabilizer_DoubleCoset_Kernel(
    std::vector<MyMatrix<T>> const l_gens, Thelper const &helper,
    MyMatrix<T> const &TheSpace, std::vector<MyMatrix<T>> const& Vmatr_gens, std::ostream &os) {
  std::pair<std::vector<MyMatrix<T>>, std::vector<DoubleCosetEntry<T>>> pair =
    LinearSpace_Stabilizer_DoubleCosetStabilizer_Kernel<T, Tgroup, Thelper>(l_gens, helper, TheSpace, Vmatr_gens, os);
  std::vector<MyMatrix<T>> l_dcs;
  for (auto & entry : pair.second) {
    l_dcs.push_back(entry.cos);
  }
  return {std::move(pair.first), std::move(l_dcs)};
}

template <typename T, typename Tgroup, typename Thelper>
RetMI_S<T,Tgroup> LinearSpace_Stabilizer(std::vector<MyMatrix<T>> const &ListMatr,
                       Thelper const &helper, MyMatrix<T> const &TheSpace,
                       std::ostream &os) {
  using Tint = typename underlying_ring<T>::ring_type;
  std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>> pair =
    LLLMatrixGroupReduction<T, Tint, Thelper>(helper, ListMatr, os);
  std::vector<MyMatrix<T>> const &ListMatrNew = pair.first;
  MyMatrix<Tint> const &Pmat = pair.second;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  MyMatrix<T> TheSpace_B = TheSpace * PmatInv_T;
  MyMatrix<T> TheSpace_C = SublatticeBasisReduction(TheSpace_B, os);
  Thelper helper_new = TransformHelper(helper, Pmat_T);
  RetMI_S<T,Tgroup> ret =
      LinearSpace_Stabilizer_Kernel<T, Tgroup, Thelper>(ListMatrNew, helper_new,
                                                        TheSpace_C, os);
  std::vector<MyMatrix<T>> const& ListMatr_B = ret.LGen;
  std::vector<MyMatrix<T>> ListMatr_C;
  for (auto &eMatr_B : ListMatr_B) {
    MyMatrix<T> eMatr_C = PmatInv_T * eMatr_B * Pmat_T;
    ListMatr_C.push_back(eMatr_C);
  }
  if (ListMatr_C.size() == 0) {
    ListMatr_C.push_back(IdentityMat<T>(helper.n));
  }
  return {ret.index, ListMatr_C};
}

template <typename T, typename Tgroup, typename Thelper>
Stab_RightCoset<T> LinearSpace_Stabilizer_RightCoset(
    std::vector<MyMatrix<T>> const &ListMatr, Thelper const &helper,
    MyMatrix<T> const &TheSpace, std::ostream &os) {
  using Tint = typename underlying_ring<T>::ring_type;
  std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>> pair =
    LLLMatrixGroupReduction<T, Tint, Thelper>(helper, ListMatr, os);
  std::vector<MyMatrix<T>> const &ListMatrNew = pair.first;
  MyMatrix<Tint> const &Pmat = pair.second;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  MyMatrix<T> TheSpace_B = TheSpace * PmatInv_T;
  MyMatrix<T> TheSpace_C = SublatticeBasisReduction(TheSpace_B, os);
  Thelper helper_new = TransformHelper(helper, Pmat_T);
  Stab_RightCoset<T> pairB =
      LinearSpace_Stabilizer_RightCoset_Kernel<T, Tgroup, Thelper>(
          ListMatrNew, helper_new, TheSpace_C, os);
  std::vector<MyMatrix<T>> ListMatr_C;
  for (auto &eMatr_B : pairB.list_gen) {
    MyMatrix<T> eMatr_C = PmatInv_T * eMatr_B * Pmat_T;
    ListMatr_C.push_back(eMatr_C);
  }
  if (ListMatr_C.size() == 0) {
    ListMatr_C.push_back(IdentityMat<T>(helper.n));
  }
  CosetDescription<T> coset = pairB.coset_desc;
  coset.conjugate(Pmat_T);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Returning from LinearSpace_Stabilizer_RightCoset\n";
#endif
  return {std::move(ListMatr_C), coset};
}

template <typename T, typename Tgroup, typename Thelper>
std::pair<std::vector<MyMatrix<T>>,std::vector<MyMatrix<T>>> LinearSpace_Stabilizer_DoubleCoset(
    std::vector<MyMatrix<T>> const &ListMatr, Thelper const &helper,
    MyMatrix<T> const &TheSpace, std::vector<MyMatrix<T>> const& V_gens, std::ostream &os) {
  using Tint = typename underlying_ring<T>::ring_type;
  std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>> pair =
    LLLMatrixGroupReduction<T, Tint, Thelper>(helper, ListMatr, os);
  std::vector<MyMatrix<T>> const &ListMatrNew = pair.first;
  MyMatrix<Tint> const &Pmat = pair.second;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  MyMatrix<T> TheSpace_B = TheSpace * PmatInv_T;
  MyMatrix<T> TheSpace_C = SublatticeBasisReduction(TheSpace_B, os);
  Thelper helper_new = TransformHelper(helper, Pmat_T);
  std::vector<MyMatrix<T>> V_gens_B;
  for (auto &eMatr_B : V_gens) {
    MyMatrix<T> eMatr_C = Pmat_T * eMatr_B * PmatInv_T;
    V_gens_B.emplace_back(std::move(eMatr_C));
  }
  std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>> pairB =
      LinearSpace_Stabilizer_DoubleCoset_Kernel<T, Tgroup, Thelper>(ListMatrNew, helper_new, TheSpace_C, V_gens_B, os);
  auto convert=[&](std::vector<MyMatrix<T>> const& l_mat) -> std::vector<MyMatrix<T>> {
    std::vector<MyMatrix<T>> l_mat_tr;
    for (auto &eMatr_B : l_mat) {
      MyMatrix<T> eMatr_C = PmatInv_T * eMatr_B * Pmat_T;
      l_mat_tr.emplace_back(std::move(eMatr_C));
    }
    return l_mat_tr;
  };
  std::vector<MyMatrix<T>> l_gen_C = convert(pairB.first);
  std::vector<MyMatrix<T>> l_cos_C = convert(pairB.second);
  return {std::move(l_gen_C), std::move(l_cos_C)};
}

template <typename T, typename Tgroup, typename Thelper>
std::pair<std::vector<MyMatrix<T>>,std::vector<DoubleCosetEntry<T>>> LinearSpace_Stabilizer_DoubleCosetStabilizer(
    std::vector<MyMatrix<T>> const &ListMatr, Thelper const &helper,
    MyMatrix<T> const &TheSpace, std::vector<MyMatrix<T>> const& V_gens, std::ostream &os) {
  using Tint = typename underlying_ring<T>::ring_type;
  std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>> pair =
    LLLMatrixGroupReduction<T, Tint, Thelper>(helper, ListMatr, os);
  std::vector<MyMatrix<T>> const &ListMatrNew = pair.first;
  MyMatrix<Tint> const &Pmat = pair.second;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  MyMatrix<T> TheSpace_B = TheSpace * PmatInv_T;
  MyMatrix<T> TheSpace_C = SublatticeBasisReduction(TheSpace_B, os);
  Thelper helper_new = TransformHelper(helper, Pmat_T);
  std::vector<MyMatrix<T>> V_gens_B;
  for (auto &eMatr_B : V_gens) {
    MyMatrix<T> eMatr_C = Pmat_T * eMatr_B * PmatInv_T;
    V_gens_B.emplace_back(std::move(eMatr_C));
  }
  std::pair<std::vector<MyMatrix<T>>, std::vector<DoubleCosetEntry<T>>> pairB =
      LinearSpace_Stabilizer_DoubleCosetStabilizer_Kernel<T, Tgroup, Thelper>(ListMatrNew, helper_new, TheSpace_C, V_gens_B, os);
  auto convert=[&](std::vector<MyMatrix<T>> const& l_mat) -> std::vector<MyMatrix<T>> {
    std::vector<MyMatrix<T>> l_mat_tr;
    for (auto &eMatr_B : l_mat) {
      MyMatrix<T> eMatr_C = PmatInv_T * eMatr_B * Pmat_T;
      l_mat_tr.emplace_back(std::move(eMatr_C));
    }
    return l_mat_tr;
  };
  std::vector<MyMatrix<T>> l_gen_C = convert(pairB.first);
  std::vector<DoubleCosetEntry<T>> l_dcs_C;
  for (auto & dcs : pairB.second) {
    MyMatrix<T> cos_C = PmatInv_T * dcs.cos * Pmat_T;
    std::vector<MyMatrix<T>> stab_gens_C = convert(dcs.stab_gens);
    DoubleCosetEntry<T> dcs_C{std::move(cos_C), std::move(stab_gens_C)};
    l_dcs_C.emplace_back(std::move(dcs_C));
  }
  return {std::move(l_gen_C), std::move(l_dcs_C)};
}


template <typename T>
using ResultTestModEquivalence =
    std::pair<std::vector<MyMatrix<T>>, MyMatrix<T>>;

/*
  We need a number of separate functions:
  ---The list of matrices has to be separated from the helper data like the
  EXTfaithful
  ---Function that creates the permutation representation given the ListMatrMat
  and the helper
  ---A function for computing the
 */
template <typename T, typename Tmod, typename Tgroup, typename Thelper>
std::optional<ResultTestModEquivalence<T>> LinearSpace_ModEquivalence_Tmod(
    std::vector<MyMatrix<T>> const &ListMatr, Thelper const &helper,
    bool const &NeedStabilizer, MyMatrix<T> const &TheSpace1,
    MyMatrix<T> const &TheSpace2, T const &TheMod, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n = TheSpace1.rows();
#ifdef DEBUG_MATRIX_GROUP
  os << "------------------------------------------------------\n";
  os << "MAT_GRP: NeedStabilizer=" << NeedStabilizer << "\n";
  os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod, TheMod=" << TheMod << "\n";
  os << "MAT_GRP: det(TheSpace1)=" << DeterminantMat(TheSpace1)
     << " det(TheSpace2)=" << DeterminantMat(TheSpace2) << "\n";
#endif
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
  MyMatrix<T> TheSpace1Mod = Concatenate(TheSpace1, ModSpace);
  MyMatrix<T> TheSpace2Mod = Concatenate(TheSpace2, ModSpace);
  std::vector<MyMatrix<T>> ListMatrRet = ListMatr;
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod(A), comp(ListMatrRet)=" << compute_complexity_listmat(ListMatrRet) << "\n";
#endif
  ListMatrRet = ExhaustiveReductionComplexityGroupMatrix(ListMatrRet, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod(B), comp(ListMatrRet)=" << compute_complexity_listmat(ListMatrRet) << "\n";
#endif
  MyMatrix<T> eElt = IdentityMat<T>(n);
  Tmod TheMod_mod = UniversalScalarConversion<Tmod, T>(TheMod);
  auto TheAction = [&](MyVector<Tmod> const &eClass,
                       MyMatrix<Tmod> const &eElt) -> MyVector<Tmod> {
    MyVector<Tmod> eVect = eElt.transpose() * eClass;
    return VectorMod(eVect, TheMod_mod);
  };
  RecSolutionIntMat<T> eCan(TheSpace2Mod);
  auto IsEquiv =
      [&](MyMatrix<T> const &eEquiv) -> std::optional<MyVector<Tmod>> {
    MyMatrix<T> TheSpace1img = TheSpace1 * eEquiv;
    for (int i = 0; i < n; i++) {
      MyVector<T> eVect = GetMatrixRow(TheSpace1img, i);
      bool test = eCan.has_solution_v(eVect);
      if (!test) {
        return ModuloReductionVector<T, Tmod>(eVect, TheMod);
      }
    }
    return {};
  };
  auto IsStabilizing = [&](std::vector<MyMatrix<T>> const &ListGen)
      -> std::optional<MyVector<Tmod>> {
    if (!NeedStabilizer)
      return {};
    for (auto &eGen : ListGen) {
      MyMatrix<T> TheSpace2img = TheSpace2 * eGen;
      for (int i = 0; i < n; i++) {
        MyVector<T> eVect = GetMatrixRow(TheSpace2img, i);
        bool test = eCan.has_solution_v(eVect);
        if (!test) {
          return ModuloReductionVector<T, Tmod>(eVect, TheMod);
        }
      }
    }
    return {};
  };
  while (true) {
    std::optional<MyVector<Tmod>> test1 = IsEquiv(eElt);
    std::optional<MyVector<Tmod>> test2 = IsStabilizing(ListMatrRet);
    if (!test1 && !test2) {
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod, eElt and GRPwork are "
            "correct. Exiting\n";
#endif
      ResultTestModEquivalence<T> res{ListMatrRet, eElt};
      return res;
    }
    if (test1) {
      MyVector<Tmod> const &V = *test1;
      std::vector<MyMatrix<Tmod>> ListMatrRetMod =
        ModuloReductionStdVectorMatrix<T, Tmod>(ListMatrRet, TheMod);
      std::vector<MyVector<Tmod>> O =
        OrbitComputation(ListMatrRetMod, V, TheAction, os);
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod, |O|=" << O.size()
         << "\n";
#endif
      Telt ePermS = Telt(SortingPerm<MyVector<Tmod>, Tidx>(O));
      std::function<Telt(MyMatrix<T> const&)> f_get_perm=[&](MyMatrix<T> const& eGen) -> Telt {
        return get_permutation_from_orbit(eGen, O, TheMod, ePermS);
      };
      int nbRow = helper.nbRow();
      MyMatrix<T> TheSpace1work = TheSpace1 * eElt;
      MyMatrix<T> TheSpace1workMod = Concatenate(TheSpace1work, ModSpace);
      Face eFace1_pre = GetFace<T, Tmod>(O, TheSpace1workMod);
      Face eFace2_pre = GetFace<T, Tmod>(O, TheSpace2Mod);
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod, |eFace1_pre|=" << eFace1_pre.size() << " / " << eFace1_pre.count() << "\n";
      os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod, |eFace2_pre|=" << eFace2_pre.size() << " / " << eFace2_pre.count() << "\n";
#endif
      PartitionReduction<T, Telt> pr(ListMatrRet, f_get_perm, eFace1_pre, os);
      Face eFace1 = TranslateFace(nbRow, pr.face);
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod, |eFace1|=" << eFace1.size() << " / " << eFace1.count() << "\n";
#endif
      std::optional<Face> opt_face2 = pr.map_face_opt(eFace2_pre);
      if (!opt_face2) {
#ifdef DEBUG_MATRIX_GROUP
        os << "MAT_GRP: Exit as no eFace2 does not map correctly to the partition\n";
#endif
        return {};
      }
      Face const& eFace2 = TranslateFace(nbRow, *opt_face2);
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod, |eFace2|=" << eFace2.size() << " / " << eFace2.count() << "\n";
#endif
      std::vector<Telt> ListPermGens =
        MatrixIntegral_GeneratePermutationGroupA<T, Telt, Thelper, decltype(f_get_perm)>(ListMatrRet, helper, pr.f_get_perm, os);
      size_t siz_act = eFace1.size();
      Tgroup GRPperm(ListPermGens, siz_act);
#ifdef SANITY_CHECK_MATRIX_GROUP
      if (eFace1.count() == 0 && eFace2.count() == 0) {
        std::cerr << "Error in LinearSpace_ModEquivalence_Tmod. |eFace1| = "
                     "|eFace2| = 0\n";
        std::cerr << "Clear bug\n";
        throw TerminalException{1};
      }
#endif
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: ModEquivalence 1 TheMod=" << TheMod << " |O|=" << O.size()
         << " |GRPperm|=" << GRPperm.size() << " |eFace1|=" << eFace1.count()
         << " |eFace2|=" << eFace2.count() << "\n";
#endif
      std::optional<MyMatrix<T>> opt =
          MatrixIntegral_RepresentativeAction<T, Tgroup, Thelper>(ListPermGens, ListMatrRet, GRPperm, helper, eFace1, eFace2, os);
      if (!opt) {
#ifdef DEBUG_MATRIX_GROUP
        os << "MAT_GRP: Exit as no equivalence exixts\n";
#endif
        return {};
      }
      MyMatrix<T> const &M = *opt;
      eElt = eElt * M;
      if (!NeedStabilizer) {
        std::optional<MyVector<Tmod>> test1 = IsEquiv(eElt);
        if (!test1) {
#ifdef DEBUG_MATRIX_GROUP
          os << "MAT_GRP: eElt and GRPwork are correct. Exiting\n";
#endif
          ResultTestModEquivalence<T> res{std::move(ListMatrRet), std::move(eElt)};
          return res;
        }
      }
      RetMI_S<T,Tgroup> ret = MatrixIntegral_Stabilizer<T, Tgroup, Thelper>(ListPermGens, ListMatrRet, pr.f_get_perm, GRPperm, helper, eFace2, os);
      ListMatrRet = ret.LGen;
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod(C), comp(ListMatrRet)=" << compute_complexity_listmat(ListMatrRet) << "\n";
#endif
      ListMatrRet = ExhaustiveReductionComplexityGroupMatrix(ListMatrRet, os);
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod(D), comp(ListMatrRet)=" << compute_complexity_listmat(ListMatrRet) << "\n";
#endif
    } else {
      MyVector<Tmod> const &V = *test2;
      std::vector<MyMatrix<Tmod>> ListMatrRetMod =
        ModuloReductionStdVectorMatrix<T, Tmod>(ListMatrRet, TheMod);
      std::vector<MyVector<Tmod>> O =
          OrbitComputation(ListMatrRetMod, V, TheAction, os);
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: |O|=" << O.size() << "\n";
#endif
      Telt ePermS = Telt(SortingPerm<MyVector<Tmod>, Tidx>(O));
      std::function<Telt(MyMatrix<T> const&)> f_get_perm=[&](MyMatrix<T> const& eGen) -> Telt {
        return get_permutation_from_orbit(eGen, O, TheMod, ePermS);
      };
      int nbRow = helper.nbRow();
      Face eFace2_pre = GetFace<T, Tmod>(O, TheSpace2Mod);
      PartitionReduction<T, Telt> pr(ListMatrRet, f_get_perm, eFace2_pre, os);
      Face eFace2 = TranslateFace(nbRow, pr.face);
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: We have eFace2\n";
#endif
      std::vector<Telt> ListPermGens =
        MatrixIntegral_GeneratePermutationGroupA<T, Telt, Thelper, decltype(f_get_perm)>(ListMatrRet, helper, pr.f_get_perm, os);
      size_t siz_act = eFace2.size();
      Tgroup GRPperm(ListPermGens, siz_act);
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: ModEquivalence 2 TheMod=" << TheMod << " |O|=" << O.size()
         << " |GRPperm|=" << GRPperm.size() << " |eFace2|=" << eFace2.count()
         << "\n";
#endif
      RetMI_S<T,Tgroup> ret = MatrixIntegral_Stabilizer<T, Tgroup, Thelper>(ListPermGens, ListMatrRet, pr.f_get_perm, GRPperm, helper, eFace2, os);
      ListMatrRet = ret.LGen;
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod(E), comp(ListMatrRet)=" << compute_complexity_listmat(ListMatrRet) << "\n";
#endif
      ListMatrRet = ExhaustiveReductionComplexityGroupMatrix(ListMatrRet, os);
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod(F), comp(ListMatrRet)=" << compute_complexity_listmat(ListMatrRet) << "\n";
#endif
    }
  }
}

template <typename T, typename Tgroup, typename Thelper>
std::optional<ResultTestModEquivalence<T>> LinearSpace_ModEquivalence(
    std::vector<MyMatrix<T>> const &ListMatr, Thelper const &helper,
    bool const &NeedStabilizer, MyMatrix<T> const &TheSpace1,
    MyMatrix<T> const &TheSpace2, T const &TheMod, std::ostream &os) {
  T max_size = (TheMod - 1) * (TheMod - 1) * TheSpace1.rows();
  if (max_size < T(std::numeric_limits<uint8_t>::max())) {
    return LinearSpace_ModEquivalence_Tmod<T, uint8_t, Tgroup, Thelper>(
        ListMatr, helper, NeedStabilizer, TheSpace1, TheSpace2, TheMod, os);
  }
  if (max_size < T(std::numeric_limits<uint16_t>::max())) {
    return LinearSpace_ModEquivalence_Tmod<T, uint16_t, Tgroup, Thelper>(
        ListMatr, helper, NeedStabilizer, TheSpace1, TheSpace2, TheMod, os);
  }
  if (max_size < T(std::numeric_limits<uint32_t>::max())) {
    return LinearSpace_ModEquivalence_Tmod<T, uint32_t, Tgroup, Thelper>(
        ListMatr, helper, NeedStabilizer, TheSpace1, TheSpace2, TheMod, os);
  }
  std::cerr << "Failed to find a matching arithmetic type. Quite unlikely "
               "objectively\n";
  throw TerminalException{1};
}

template <typename T, typename Tgroup, typename Thelper>
std::optional<MyMatrix<T>>
LinearSpace_Equivalence_KernelRing(std::vector<MyMatrix<T>> const &ListMatr,
                                   Thelper const &helper,
                                   MyMatrix<T> const &TheSpace1,
                                   MyMatrix<T> const &TheSpace2, std::ostream &os) {
  int n = helper.n;
  T LFact1 = LinearSpace_GetDivisor(TheSpace1);
  T LFact2 = LinearSpace_GetDivisor(TheSpace2);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LFact1 = " << LFact1 << "\n";
  os << "MAT_GRP: LFact2 = " << LFact2 << "\n";
#endif
  if (LFact1 != LFact2) {
    return {};
  }
  std::vector<T> eList = FactorsInt(LFact1);
  RecSolutionIntMat<T> eCan(TheSpace2);
  auto IsEquivalence = [&](MyMatrix<T> const &eEquiv) -> bool {
    MyMatrix<T> TheSpace1img = TheSpace1 * eEquiv;
    return eCan.is_containing_m(TheSpace1img);
  };
  std::vector<MyMatrix<T>> ListMatrWork = ListMatr;
  int siz = eList.size();
  std::vector<MyMatrix<T>> ListMatrRet = ListMatr;
  MyMatrix<T> eElt = IdentityMat<T>(n);
  for (int i = 1; i <= siz; i++) {
    if (IsEquivalence(eElt))
      return eElt;
    T TheMod(1);
    for (int j = 0; j < i; j++)
      TheMod *= eList[j];
    MyMatrix<T> TheSpace1Img = TheSpace1 * eElt;
    bool NeedStabilizer = true;
    if (i == siz)
      NeedStabilizer = false;
    std::optional<ResultTestModEquivalence<T>> opt =
        LinearSpace_ModEquivalence<T, Tgroup, Thelper>(
            ListMatrWork, helper, NeedStabilizer, TheSpace1Img, TheSpace2,
            TheMod, os);
    if (!opt) {
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModEquivalence failed so we exit here\n";
#endif
      return {};
    }
    eElt = eElt * (opt->second);
    if (NeedStabilizer) {
      ListMatrWork = opt->first;
    }
  }
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsEquivalence(eElt)) {
    std::cerr << "Error in LinearSpace_Equivalence_Kernel\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Before returning from LinearSpace_Equivalence_Kernel, retuning eElt\n";
#endif
  return eElt;
}

template <typename T, typename Tgroup, typename Thelper>
std::optional<MyMatrix<T>>
LinearSpace_Equivalence_Kernel(std::vector<MyMatrix<T>> const &ListMatr,
                               Thelper const &helper,
                               MyMatrix<T> const &InSpace1,
                               MyMatrix<T> const &InSpace2, std::ostream &os) {
  using Tint = typename Thelper::Tint;
  using ThelperInt = typename Thelper::ThelperInt;
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Beginning of LinearSpace_Equivalence_Kernel\n";
  os << "MAT_GRP: |ListMatr|=" << ListMatr.size() << "\n";
  os << "MAT_GRP: Det(InSpace1)=" << DeterminantMat(InSpace1)
     << " Det(InSpace2)=" << DeterminantMat(InSpace2) << "\n";
#endif
  FractionMatrix<T> eRec1 = RemoveFractionMatrixPlusCoeff(InSpace1);
  FractionMatrix<T> eRec2 = RemoveFractionMatrixPlusCoeff(InSpace2);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: eRec1.TheMult=" << eRec1.TheMult
     << " eRec2.TheMult=" << eRec2.TheMult << "\n";
#endif
  if (eRec1.TheMult != eRec2.TheMult)
    return {};
  MyMatrix<T> const &TheSpace1 = eRec1.TheMat;
  MyMatrix<T> const &TheSpace2 = eRec2.TheMat;
  MyMatrix<Tint> TheSpace1_int = UniversalMatrixConversion<Tint,T>(TheSpace1);
  MyMatrix<Tint> TheSpace2_int = UniversalMatrixConversion<Tint,T>(TheSpace2);
  //
  ThelperInt helper_int = ToInteger(helper);
  std::vector<MyMatrix<Tint>> ListMatr_int = UniversalStdVectorMatrixConversion<Tint,T>(ListMatr);
  //
  std::optional<MyMatrix<Tint>> opt =
    LinearSpace_Equivalence_KernelRing<Tint,Tgroup,ThelperInt>(ListMatr_int, helper_int, TheSpace1_int, TheSpace2_int, os);
  if (!opt) {
    return {};
  }
  MyMatrix<Tint> const& M = *opt;
  MyMatrix<T> M_T = UniversalMatrixConversion<T,Tint>(M);
  return M_T;
}



template <typename T, typename Tgroup, typename Thelper>
std::optional<MyMatrix<T>>
LinearSpace_Equivalence(std::vector<MyMatrix<T>> const &ListMatr,
                        Thelper const &helper, MyMatrix<T> const &InSpace1,
                        MyMatrix<T> const &InSpace2, std::ostream &os) {
  using Tint = typename underlying_ring<T>::ring_type;
  std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>> pair =
    LLLMatrixGroupReduction<T, Tint, Thelper>(helper, ListMatr, os);
  std::vector<MyMatrix<T>> const &ListMatrNew = pair.first;
  MyMatrix<Tint> const &Pmat = pair.second;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  MyMatrix<T> InSpace1_B = InSpace1 * PmatInv_T;
  MyMatrix<T> InSpace2_B = InSpace2 * PmatInv_T;
  MyMatrix<T> InSpace1_C = SublatticeBasisReduction(InSpace1_B, os);
  MyMatrix<T> InSpace2_C = SublatticeBasisReduction(InSpace2_B, os);
  Thelper helper_new = TransformHelper(helper, Pmat_T);
  std::optional<MyMatrix<T>> opt =
      LinearSpace_Equivalence_Kernel<T, Tgroup, Thelper>(
          ListMatrNew, helper_new, InSpace1_C, InSpace2_C, os);
  if (!opt)
    return {};
  MyMatrix<T> RetMat = PmatInv_T * (*opt) * Pmat_T;
  return RetMat;
}

template <typename T, typename Tgroup>
RetMI_S<T,Tgroup> LinPolytopeIntegral_Automorphism_Subspaces(
    std::vector<MyMatrix<T>> const &ListMatr, MyMatrix<T> const &EXTfaithful,
    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> eBasis = GetZbasis(EXTfaithful);
  MyMatrix<T> InvBasis = Inverse(eBasis);
  MyMatrix<T> EXTbas = EXTfaithful * InvBasis;
  std::vector<MyMatrix<T>> ListMatrGens;
  for (auto &eGen : ListMatr) {
    MyMatrix<T> NewGen = eBasis * eGen * InvBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrGens.emplace_back(std::move(NewGen));
  }
  FiniteMatrixGroupHelper<T, Telt, TintGroup> helper =
    ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXTbas);
  MyMatrix<T> LattToStab = RemoveFractionMatrix(Inverse(eBasis));

  RetMI_S<T,Tgroup> ret =
      LinearSpace_Stabilizer<T, Tgroup>(ListMatrGens, helper, LattToStab, os);
  std::vector<MyMatrix<T>> ListMatrGensB;
  for (auto &eGen : ret.LGen) {
    MyMatrix<T> NewGen = InvBasis * eGen * eBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrGensB.push_back(NewGen);
  }
  return {ret.index, ListMatrGensB};
}

template <typename T, typename Tgroup>
Stab_RightCoset<T> LinPolytopeIntegral_Automorphism_RightCoset_Subspaces(
    std::vector<MyMatrix<T>> const &ListMatr, MyMatrix<T> const &EXTfaithful,
    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> eBasis = GetZbasis(EXTfaithful);
  MyMatrix<T> InvBasis = Inverse(eBasis);
  MyMatrix<T> EXTbas = EXTfaithful * InvBasis;
  std::vector<MyMatrix<T>> ListMatrGens;
  for (auto &eGen : ListMatr) {
    MyMatrix<T> NewGen = eBasis * eGen * InvBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrGens.emplace_back(std::move(NewGen));
  }
  FiniteMatrixGroupHelper<T, Telt, TintGroup> helper =
    ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXTbas);
  MyMatrix<T> LattToStab = RemoveFractionMatrix(Inverse(eBasis));
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinPolytopeIntegral_Automorphism_RightCoset_Subspaces, "
        "before LinearSpace_Stabilizer_RightCoset\n";
#endif
  Stab_RightCoset<T> pair = LinearSpace_Stabilizer_RightCoset<T, Tgroup>(
      ListMatrGens, helper, LattToStab, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinPolytopeIntegral_Automorphism_RightCoset_Subspaces, after "
        "LinearSpace_Stabilizer_RightCoset\n";
#endif
  std::vector<MyMatrix<T>> ListMatrGensB;
  for (auto &eGen : pair.list_gen) {
    MyMatrix<T> NewGen = InvBasis * eGen * eBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrGensB.emplace_back(std::move(NewGen));
  }
  pair.coset_desc.conjugate(eBasis);
  return {std::move(ListMatrGensB), pair.coset_desc};
}

template <typename T, typename Tgroup>
std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>>
LinPolytopeIntegral_Automorphism_DoubleCoset_Subspaces(std::vector<MyMatrix<T>> const &ListMatrFull, std::vector<MyMatrix<T>> const& ListMatrV, MyMatrix<T> const &EXTfaithful, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> eBasis = GetZbasis(EXTfaithful);
  MyMatrix<T> InvBasis = Inverse(eBasis);
  MyMatrix<T> EXTbas = EXTfaithful * InvBasis;
  std::vector<MyMatrix<T>> ListMatrFullGens;
  for (auto &eGen : ListMatrFull) {
    MyMatrix<T> NewGen = eBasis * eGen * InvBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrFullGens.emplace_back(std::move(NewGen));
  }
  std::vector<MyMatrix<T>> ListMatrVGens;
  for (auto &eGen : ListMatrV) {
    MyMatrix<T> NewGen = eBasis * eGen * InvBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrVGens.emplace_back(std::move(NewGen));
  }
  FiniteMatrixGroupHelper<T, Telt, TintGroup> helper =
    ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXTbas);
  MyMatrix<T> LattToStab = RemoveFractionMatrix(Inverse(eBasis));
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinPolytopeIntegral_Automorphism_RightCoset_Subspaces, "
        "before LinearSpace_Stabilizer_RightCoset\n";
#endif
  std::pair<std::vector<MyMatrix<T>>,std::vector<MyMatrix<T>>> pair =
    LinearSpace_Stabilizer_DoubleCoset<T, Tgroup>(ListMatrFullGens, helper, LattToStab, ListMatrVGens, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinPolytopeIntegral_Automorphism_DoubleCoset_Subspaces, after "
        "LinearSpace_Stabilizer_DoubleCoset\n";
#endif
  std::vector<MyMatrix<T>> ListMatrStabGens;
  for (auto &eGen : pair.first) {
    MyMatrix<T> NewGen = InvBasis * eGen * eBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrStabGens.emplace_back(std::move(NewGen));
  }
  std::vector<MyMatrix<T>> ListDoubleCosets;
  for (auto &eCos : pair.second) {
    MyMatrix<T> NewCos = InvBasis * eCos * eBasis;
    ListDoubleCosets.emplace_back(std::move(NewCos));
  }
  return {std::move(ListMatrStabGens), std::move(ListDoubleCosets)};
}

template <typename T, typename Tgroup>
std::pair<std::vector<MyMatrix<T>>, std::vector<DoubleCosetEntry<T>>>
LinPolytopeIntegral_Automorphism_DoubleCosetStabilizer_Subspaces(std::vector<MyMatrix<T>> const &ListMatrFull, std::vector<MyMatrix<T>> const& ListMatrV, MyMatrix<T> const &EXTfaithful, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> eBasis = GetZbasis(EXTfaithful);
  MyMatrix<T> InvBasis = Inverse(eBasis);
  MyMatrix<T> EXTbas = EXTfaithful * InvBasis;
  std::vector<MyMatrix<T>> ListMatrFullGens;
  for (auto &eGen : ListMatrFull) {
    MyMatrix<T> NewGen = eBasis * eGen * InvBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrFullGens.emplace_back(std::move(NewGen));
  }
  std::vector<MyMatrix<T>> ListMatrVGens;
  for (auto &eGen : ListMatrV) {
    MyMatrix<T> NewGen = eBasis * eGen * InvBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrVGens.emplace_back(std::move(NewGen));
  }
  FiniteMatrixGroupHelper<T, Telt, TintGroup> helper =
    ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXTbas);
  MyMatrix<T> LattToStab = RemoveFractionMatrix(Inverse(eBasis));
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinPolytopeIntegral_Automorphism_RightCoset_Subspaces, "
        "before LinearSpace_Stabilizer_RightCoset\n";
#endif
  std::pair<std::vector<MyMatrix<T>>,std::vector<DoubleCosetEntry<T>>> pair =
    LinearSpace_Stabilizer_DoubleCosetStabilizer<T, Tgroup>(ListMatrFullGens, helper, LattToStab, ListMatrVGens, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinPolytopeIntegral_Automorphism_DoubleCoset_Subspaces, after "
        "LinearSpace_Stabilizer_DoubleCoset\n";
#endif
  std::vector<MyMatrix<T>> ListMatrStabGens;
  for (auto &eGen : pair.first) {
    MyMatrix<T> NewGen = InvBasis * eGen * eBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrStabGens.emplace_back(std::move(NewGen));
  }
  std::vector<DoubleCosetEntry<T>> ListDoubleCosetStabilizer;
  for (auto &dcs : pair.second) {
    MyMatrix<T> NewCos = InvBasis * dcs.cos * eBasis;
    std::vector<MyMatrix<T>> new_stab_gens;
    for (auto & eGen : dcs.stab_gens) {
      MyMatrix<T> NewGen = InvBasis * eGen * eBasis;
      new_stab_gens.emplace_back(std::move(NewGen));
    }
    DoubleCosetEntry<T> new_dcs{std::move(NewCos), std::move(new_stab_gens)};
    ListDoubleCosetStabilizer.emplace_back(std::move(new_dcs));
  }
  return {std::move(ListMatrStabGens), std::move(ListDoubleCosetStabilizer)};
}

template <typename T>
std::vector<MyMatrix<T>>
ConjugateListGeneratorsTestInt(MyMatrix<T> const &Pmat,
                               std::vector<MyMatrix<T>> const &LGen) {
  std::vector<MyMatrix<T>> LGen2;
  MyMatrix<T> PmatInv = Inverse(Pmat);
  for (auto &eGen1 : LGen) {
    MyMatrix<T> eGen2 = PmatInv * eGen1 * Pmat;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(eGen2)) {
      std::cerr << "MAT_GRP: Error in ConjugateListGeneratorsTestInt\n";
      std::cerr << "MAT_GRP: The matrix eGen2 should be integral\n";
      std::cerr << "MAT_GRP: eGen2=\n";
      WriteMatrix(std::cerr, eGen2);
      throw TerminalException{1};
    }
#endif
    LGen2.emplace_back(std::move(eGen2));
  }
  return LGen2;
}

/*
  Compute the intersection of G \cap GL_n(Z)
 */
template <typename T, typename Tint, typename Tgroup>
RetMI_S<Tint,Tgroup> MatrixIntegral_Stabilizer_General(
    int const &n, std::vector<MyMatrix<T>> const &LGen1, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  using Thelper = GeneralMatrixGroupHelper<T, Telt, TintGroup>;
  MyMatrix<T> InvariantSpace = MatrixIntegral_GetInvariantSpace(n, LGen1, os);
  MyMatrix<T> InvInvariantSpace = Inverse(InvariantSpace);
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsIntegralMatrix(InvInvariantSpace)) {
    std::cerr << "The matrix InvInvariantSpace should be integral\n";
    throw TerminalException{1};
  }
#endif
  std::vector<MyMatrix<T>> LGen2 =
      ConjugateListGeneratorsTestInt(InvInvariantSpace, LGen1);
  Thelper helper{n};
  RetMI_S<T,Tgroup> ret =
    LinearSpace_Stabilizer<T, Tgroup, Thelper>(
          LGen2, helper, InvInvariantSpace, os);
  std::vector<MyMatrix<Tint>> LGen4;
  for (auto &eGen3 : ret.LGen) {
    MyMatrix<T> eGen4_T = InvInvariantSpace * eGen3 * InvariantSpace;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(eGen4_T)) {
      std::cerr << "The matrix eGen4_T should be integral\n";
      throw TerminalException{1};
    }
#endif
    MyMatrix<Tint> eGen4 = UniversalMatrixConversion<Tint, T>(eGen4_T);
    LGen4.emplace_back(std::move(eGen4));
  }
  return {ret.index, LGen4};
}

// Instead of Z^n, we now want "Sublattice" to be preserved.
template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<T>> MatrixIntegral_Stabilizer_Sublattice(
    MyMatrix<T> const& Sublattice, std::vector<MyMatrix<T>> const &LGen1, std::ostream &os) {
  MyMatrix<T> SublatticeInv = Inverse(Sublattice);
  int n = SublatticeInv.rows();
  std::vector<MyMatrix<T>> LGen2;
  for (auto & eGen1 : LGen1) {
    MyMatrix<T> eGen2 = SublatticeInv * eGen1 * Sublattice;
    LGen2.push_back(eGen2);
  }
  std::vector<MyMatrix<Tint>> LGen3 = MatrixIntegral_Stabilizer_General<T,Tint,Tgroup>(n, LGen2, os);
  std::vector<MyMatrix<T>> LGen4;
  for (auto & eGen3 : LGen3) {
    MyMatrix<T> eGen3_T = UniversalMatrixConversion<T,Tint>(eGen3);
    MyMatrix<T> eGen4 = Sublattice * eGen3_T * SublatticeInv;
    LGen4.push_back(eGen4);
  }
  return LGen4;
}


// Returns an element g in GRPrat such that   g * EquivRat   in GL(n,Z)
template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
MatrixIntegral_Equivalence_General(std::vector<MyMatrix<T>> const &LGen1,
                                   MyMatrix<T> const &EquivRat,
                                   std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  int n = EquivRat.rows();
  MyMatrix<T> TheSpace = MatrixIntegral_GetInvariantSpace(n, LGen1, os);
  MyMatrix<T> TheSpaceInv = Inverse(TheSpace);
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsIntegralMatrix(TheSpaceInv)) {
    std::cerr << "The matrix InvInvariantSpace should be integral\n";
    throw TerminalException{1};
  }
#endif
  std::vector<MyMatrix<T>> LGen2 =
      ConjugateListGeneratorsTestInt(TheSpaceInv, LGen1);
  // We search for g in GRPrat s.t. g * EquivRat in GL_n(Z).
  // So, we search g in GRPrat s.t. Z^n * g * EquivRat = Z^n
  // Writing g = TheSpaceInv g_int TheSpace we get
  // TheSpaceInv g TheSpace EquivRat = Z^n
  // Or TheSpaceInv g = Inverse(TheSpace * EquivRat)
  MyMatrix<T> TheSpaceImg = TheSpace * EquivRat;
  MyMatrix<T> TheSpaceImgInv = Inverse(TheSpaceImg);
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsIntegralMatrix(TheSpaceImgInv)) {
    std::cerr << "The matrix TheSpaceImgInv should be integral\n";
    throw TerminalException{1};
  }
#endif
  using Thelper = GeneralMatrixGroupHelper<T, Telt, TintGroup>;
  Thelper helper{n};
  std::optional<MyMatrix<T>> opt =
      LinearSpace_Equivalence_Kernel<T, Tgroup, Thelper>(
          LGen2, helper, TheSpaceInv, TheSpaceImgInv, os);
  if (!opt) {
    return {};
  }
  MyMatrix<T> const &eSpaceEquiv = *opt;
  MyMatrix<T> eMatFinal = TheSpaceInv * eSpaceEquiv * TheSpace;
  MyMatrix<T> eProd_T = eMatFinal * EquivRat;
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsIntegralMatrix(eProd_T)) {
    std::cerr << "The matrix should be integral\n";
    throw TerminalException{1};
  }
#endif
  return UniversalMatrixConversion<Tint, T>(eProd_T);
}

// Find a matrix g in GRPrat such that   EquivRat * g   in   GL(n,Z)
template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
MatrixIntegral_Equivalence_Bis_General(std::vector<MyMatrix<T>> const &GRPrat,
                                       MyMatrix<T> const &EquivRat,
                                       std::ostream &os) {
  MyMatrix<T> EquivRatInv = Inverse(EquivRat);
  std::optional<MyMatrix<Tint>> opt =
      MatrixIntegral_Equivalence_General<T, Tint, Tgroup>(GRPrat, EquivRatInv,
                                                          os);
  if (!opt) {
    return {};
  }
  MyMatrix<Tint> const &TheSol = *opt;
  // So we have TheSol = g * Inverse(EquivRat) in GL(n,Z)
  // Inverse(TheSol) = EquivRat * g in GL(n,Z)
  return Inverse(TheSol);
}

template <typename T, typename Tgroup>
std::vector<MyMatrix<T>> MatrixIntegral_RightCosets_General(
    int const &n, std::vector<MyMatrix<T>> const &LGen1, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  using Thelper = GeneralMatrixGroupHelper<T, Telt, TintGroup>;
  MyMatrix<T> InvariantSpace = MatrixIntegral_GetInvariantSpace(n, LGen1, os);
  MyMatrix<T> InvInvariantSpace = Inverse(InvariantSpace);
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsIntegralMatrix(InvInvariantSpace)) {
    std::cerr << "MAT_GRP: The matrix InvInvariantSpace should be integral\n";
    throw TerminalException{1};
  }
#endif
  std::vector<MyMatrix<T>> LGen2 =
      ConjugateListGeneratorsTestInt(InvInvariantSpace, LGen1);
  Thelper helper{n};
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_RightCosets_General, before "
        "LinearSpace_Stabilizer_RightCoset\n";
#endif
  Stab_RightCoset<T> stab_rightcoset =
      LinearSpace_Stabilizer_RightCoset<T, Tgroup, Thelper>(
          LGen2, helper, InvInvariantSpace, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_RightCosets_General, after "
        "LinearSpace_Stabilizer_RightCoset\n";
#endif
  using Iter = typename CosetDescription<T>::const_iterator;
  Iter iter = stab_rightcoset.coset_desc.begin();
  Iter end = stab_rightcoset.coset_desc.end();
  std::vector<MyMatrix<T>> LCoset2;
  while (iter != end) {
    MyMatrix<T> const &eCos1 = *iter;
    MyMatrix<T> eCos2 = InvInvariantSpace * eCos1 * InvariantSpace;
    LCoset2.push_back(eCos2);
    iter++;
  }
  return LCoset2;
}

template <typename T, typename Tgroup>
std::pair<std::vector<MyMatrix<T>>,std::vector<MyMatrix<T>>>
MatrixIntegral_DoubleCosets_General(
    int const &n, std::vector<MyMatrix<T>> const &LGenG1,
    std::vector<MyMatrix<T>> const &LGenV1,
    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  using Thelper = GeneralMatrixGroupHelper<T, Telt, TintGroup>;
#ifdef DEBUG_MATRIX_GROUP_DISABLE
  os << "MAT_GRP: MatrixIntegral_DoubleCosets_General, LGenG1=\n";
  WriteListMatrix(os, LGenG1);
  os << "MAT_GRP: MatrixIntegral_DoubleCosets_General, LGenV1=\n";
  WriteListMatrix(os, LGenV1);
#endif
  MyMatrix<T> InvariantSpace = MatrixIntegral_GetInvariantSpace(n, LGenG1, os);
  MyMatrix<T> InvInvariantSpace = Inverse(InvariantSpace);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: We have |InvariantSpace|=" << DeterminantMat(InvariantSpace) << "\n";
  os << "MAT_GRP: We have InvariantSpace=\n";
  WriteMatrix(os, InvariantSpace);
  os << "MAT_GRP: We have InvInvariantSpace=\n";
  WriteMatrix(os, InvInvariantSpace);
#endif
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsIntegralMatrix(InvInvariantSpace)) {
    std::cerr << "MAT_GRP: The matrix InvInvariantSpace should be integral\n";
    throw TerminalException{1};
  }
#endif
  std::vector<MyMatrix<T>> LGenG2 =
      ConjugateListGeneratorsTestInt(InvInvariantSpace, LGenG1);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: We have LGenG2\n";
#endif
  std::vector<MyMatrix<T>> LGenV2 =
      ConjugateListGeneratorsTestInt(InvInvariantSpace, LGenV1);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: We have LGenV2\n";
#endif
  Thelper helper{n};
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_DoubleCosets_General, before "
        "LinearSpace_Stabilizer_DoubleCoset\n";
#endif
  std::pair<std::vector<MyMatrix<T>>,std::vector<MyMatrix<T>>> pair =
      LinearSpace_Stabilizer_DoubleCoset<T, Tgroup, Thelper>(LGenG2, helper, InvInvariantSpace, LGenV2, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_DoubleCosets_General, after "
        "LinearSpace_Stabilizer_DoubleCoset\n";
#endif
  std::vector<MyMatrix<T>> LGenRet;
  for (auto & eGen1 : pair.first) {
    MyMatrix<T> eGen2 = InvInvariantSpace * eGen1 * InvariantSpace;
    LGenRet.push_back(eGen2);
  }
  std::vector<MyMatrix<T>> LCosRet;
  for (auto & eCos1 : pair.second) {
    MyMatrix<T> eCos2 = InvInvariantSpace * eCos1 * InvariantSpace;
    LCosRet.push_back(eCos2);
  }
  return {LGenRet, LCosRet};
}



template <typename T, typename Tgroup>
Tgroup LinPolytopeIntegral_Stabilizer_Method8(MyMatrix<T> const &EXT_T,
                                              Tgroup const &GRPisom,
                                              std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  int nbVert = EXT_T.rows();
  std::vector<MyMatrix<T>> ListMatrGen;
  for (auto &eGen : GRPisom.SmallGeneratingSet()) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    ListMatrGen.push_back(eMat);
  }
  using Thelper = FiniteMatrixGroupHelper<T, Telt, TintGroup>;
  Thelper helper = ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXT_T);
  RetMI_S<T,Tgroup> ret =
      LinPolytopeIntegral_Automorphism_Subspaces<T, Tgroup>(ListMatrGen, EXT_T,
                                                            os);
  std::vector<Telt> ListPermGens;
  for (auto &eMatr : ret.LGen) {
    Telt ePermGen = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(
        helper, eMatr, os);
    ListPermGens.emplace_back(std::move(ePermGen));
  }
  return Tgroup(ListPermGens, nbVert);
}

template <typename T, typename Tgroup>
std::pair<Tgroup, std::vector<typename Tgroup::Telt>>
LinPolytopeIntegral_Stabilizer_RightCoset_Method8(MyMatrix<T> const &EXT_T,
                                                  Tgroup const &GRPisom,
                                                  std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  int nbVert = EXT_T.rows();
  std::vector<MyMatrix<T>> ListMatrGen;
  for (auto &eGen : GRPisom.SmallGeneratingSet()) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    ListMatrGen.push_back(eMat);
  }
  using Thelper = FiniteMatrixGroupHelper<T, Telt, TintGroup>;
  Thelper helper = ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXT_T);
  Stab_RightCoset<T> pair =
      LinPolytopeIntegral_Automorphism_RightCoset_Subspaces<T, Tgroup>(
          ListMatrGen, EXT_T, os);
  std::vector<Telt> ListPermGens;
  for (auto &eMatr : pair.list_gen) {
    Telt ePermGen = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(
        helper, eMatr, os);
    ListPermGens.emplace_back(std::move(ePermGen));
  }
  Tgroup GRPret(ListPermGens, nbVert);
  std::vector<Telt> RightCoset;
  for (auto eMatr : pair.coset_desc) {
    Telt eCos = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(helper, eMatr, os);
    RightCoset.emplace_back(std::move(eCos));
  }
  return {std::move(GRPret), std::move(RightCoset)};
}

template <typename T, typename Tgroup>
std::pair<Tgroup, std::vector<typename Tgroup::Telt>>
LinPolytopeIntegral_Stabilizer_DoubleCoset_Method8(MyMatrix<T> const &EXT_T,
                                                   Tgroup const &GRPfull,
                                                   Tgroup const &GrpV,
                                                   std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  int nbVert = EXT_T.rows();
  std::vector<MyMatrix<T>> ListMatrGenFull;
  for (auto &eGen : GRPfull.SmallGeneratingSet()) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    ListMatrGenFull.push_back(eMat);
  }
  std::vector<MyMatrix<T>> ListMatrGenV;
  for (auto &eGen : GrpV.SmallGeneratingSet()) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    ListMatrGenV.push_back(eMat);
  }
  using Thelper = FiniteMatrixGroupHelper<T, Telt, TintGroup>;
  Thelper helper = ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXT_T);
  std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>> pair =
      LinPolytopeIntegral_Automorphism_DoubleCoset_Subspaces<T, Tgroup>(ListMatrGenFull, ListMatrGenV, EXT_T, os);
  std::vector<Telt> ListPermGens;
  for (auto &eMatr : pair.first) {
    Telt ePermGen = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(
        helper, eMatr, os);
    ListPermGens.emplace_back(std::move(ePermGen));
  }
  Tgroup GRPret(ListPermGens, nbVert);
  std::vector<Telt> DoubleCosets;
  for (auto eMatr : pair.second) {
    Telt eCos = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(helper, eMatr, os);
    DoubleCosets.emplace_back(std::move(eCos));
  }
  return {std::move(GRPret), std::move(DoubleCosets)};
}

template<typename Telt>
struct PairCosetStabGens {
  Telt cos;
  std::vector<Telt> stab_gens;
};

template <typename T, typename Tgroup>
std::pair<Tgroup, std::vector<PairCosetStabGens<typename Tgroup::Telt>>>
LinPolytopeIntegral_Stabilizer_DoubleCosetStabilizer_Method8(MyMatrix<T> const &EXT_T,
                                                             Tgroup const &GRPfull,
                                                             Tgroup const &GrpV,
                                                             std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  int nbVert = EXT_T.rows();
  std::vector<MyMatrix<T>> ListMatrGenFull;
  for (auto &eGen : GRPfull.SmallGeneratingSet()) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    ListMatrGenFull.push_back(eMat);
  }
  std::vector<MyMatrix<T>> ListMatrGenV;
  for (auto &eGen : GrpV.SmallGeneratingSet()) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    ListMatrGenV.push_back(eMat);
  }
  using Thelper = FiniteMatrixGroupHelper<T, Telt, TintGroup>;
  Thelper helper = ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXT_T);
  std::pair<std::vector<MyMatrix<T>>, std::vector<DoubleCosetEntry<T>>> pair =
      LinPolytopeIntegral_Automorphism_DoubleCosetStabilizer_Subspaces<T, Tgroup>(ListMatrGenFull, ListMatrGenV, EXT_T, os);
  std::vector<Telt> ListPermGens;
  for (auto &eMatr : pair.first) {
    Telt ePermGen = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(
        helper, eMatr, os);
    ListPermGens.emplace_back(std::move(ePermGen));
  }
  Tgroup GRPret(ListPermGens, nbVert);
  std::vector<PairCosetStabGens<Telt>> DoubleCosetStabilizer;
  for (auto eDCS : pair.second) {
    Telt NewCos = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(helper, eDCS.cos, os);
    std::vector<Telt> new_stab_gens;
    for (auto & eMatr : eDCS.stab_gens) {
      Telt NewGen = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(helper, eMatr, os);
      new_stab_gens.emplace_back(std::move(NewGen));
    }
    PairCosetStabGens<Telt> NewDCS{std::move(NewCos), std::move(new_stab_gens)};
    DoubleCosetStabilizer.emplace_back(std::move(NewDCS));
  }
  return {std::move(GRPret), std::move(DoubleCosetStabilizer)};
}

template <typename T, typename Tgroup>
std::optional<MyMatrix<T>> LinPolytopeIntegral_Isomorphism_Subspaces(
    MyMatrix<T> const &EXT1_T, MyMatrix<T> const &EXT2_T,
    std::vector<MyMatrix<T>> const &ListMatrGens2,
    typename Tgroup::Telt const &eEquiv, std::ostream &os) {
#ifdef TIMINGS_MATRIX_GROUP
  MicrosecondTime time;
#endif
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Beginning of LinPolytopeIntegral_Isomorphism_Subspaces\n";
  os << "MAT_GRP: |EXT1_T|=" << EXT1_T.rows() << " |EXT2_T|=" << EXT2_T.rows() << "\n";
#endif
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> eBasis1 = GetZbasis(EXT1_T);
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: GetZbasis1|=" << time << "\n";
#endif
  MyMatrix<T> eBasis2 = GetZbasis(EXT2_T);
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: GetZbasis2|=" << time << "\n";
#endif
  MyMatrix<T> InvBasis1 = Inverse(eBasis1);
  MyMatrix<T> InvBasis2 = Inverse(eBasis2);
  MyMatrix<T> EXTbas1 = EXT1_T * InvBasis1;
  MyMatrix<T> EXTbas2 = EXT2_T * InvBasis2;
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: InvBasisX / EXTbasX|=" << time << "\n";
#endif
#ifdef SANITY_CHECK_MATRIX_GROUP
  using Tidx = typename Telt::Tidx;
  for (auto &eMatGen2 : ListMatrGens2) {
    std::optional<std::vector<Tidx>> opt_eList =
        RepresentVertexPermutationTest<T, T, Tidx>(EXT2_T, EXT2_T, eMatGen2);
    if (!opt_eList) {
      std::cerr << "LinPolytopeIntegral_Isomorphism_Subspaces: We fail to "
                   "represent the matrix as a permutation of the rows\n";
      throw TerminalException{1};
    }
  }
# ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: RepresentVertexPermutationTest|=" << time << "\n";
# endif
#endif
  //
  MyMatrix<T> TheMatEquiv = FindTransformation(EXTbas1, EXTbas2, eEquiv);
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: FindTransformation|=" << time << "\n";
#endif
  std::vector<MyMatrix<T>> ListMatrGen;
  for (auto &eGen : ListMatrGens2) {
    MyMatrix<T> NewGen = eBasis2 * eGen * InvBasis2;
    ListMatrGen.push_back(NewGen);
  }
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: ListMatrGen|=" << time << "\n";
#endif
  FiniteMatrixGroupHelper<T, Telt, TintGroup> helper =
    ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXTbas2);
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: helper|=" << time << "\n";
#endif
  MyMatrix<T> eLatt1 = Inverse(eBasis1) * TheMatEquiv;
  MyMatrix<T> eLatt2 = Inverse(eBasis2);
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: eLattX|=" << time << "\n";
#endif
  std::optional<MyMatrix<T>> opt = LinearSpace_Equivalence<T, Tgroup>(
      ListMatrGen, helper, eLatt1, eLatt2, os);
#ifdef TIMINGS_MATRIX_GROUP
  os << "|MAT_GRP: LinearSpace_Equivalence|=" << time << "\n";
#endif
  if (!opt)
    return {};
  MyMatrix<T> const &eSpaceEquiv = *opt;
  MyMatrix<T> eMatFinal = InvBasis1 * TheMatEquiv * eSpaceEquiv * eBasis2;
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsIntegralMatrix(eMatFinal)) {
    std::cerr << "LinPolytopeIntegral_Isomorphism_Subspaces: eMatFinal should "
                 "be integral\n";
    throw TerminalException{1};
  }
#endif
  return eMatFinal;
}

// GRP1 is a group of automorphism preserving EXT1_T
// ePerm is a transformation mapping EXT1 to EXT2.
// We are searching for a transformation h in GRP1 such that
// h * ePerm is an integral transformation.
template <typename T, typename Tgroup>
std::optional<MyMatrix<T>> LinPolytopeIntegral_Isomorphism_Method8(
    MyMatrix<T> const &EXT1_T, MyMatrix<T> const &EXT2_T, Tgroup const &GRP1,
    typename Tgroup::Telt const &ePerm, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  std::vector<MyMatrix<T>> ListMatrGens;
  std::vector<Telt> LGen = GRP1.SmallGeneratingSet();
  for (auto &eGen : LGen) {
    Telt ePermGen = (~ePerm) * eGen * ePerm;
    MyMatrix<T> eMatr = FindTransformation(EXT2_T, EXT2_T, ePermGen);
    ListMatrGens.push_back(eMatr);
  }
  return LinPolytopeIntegral_Isomorphism_Subspaces<T, Tgroup>(
      EXT1_T, EXT2_T, ListMatrGens, ePerm, os);
}

template<typename Telt, typename T>
std::string get_matrs_as_string(MyMatrix<T> const& EXT, std::vector<Telt> const &l_elt) {
  std::string strGAPmatr = "[";
  bool IsFirst = true;
  for (auto &eElt : l_elt) {
    MyMatrix<T> M = RepresentVertexPermutation(EXT, EXT, eElt);
    if (!IsFirst)
      strGAPmatr += ",";
    IsFirst = false;
    strGAPmatr += StringMatrixGAP(M);
  }
  strGAPmatr += "]";
  return strGAPmatr;
}

// clang-format off
#endif  // SRC_GROUP_MATRIXGROUP_H_
// clang-format on

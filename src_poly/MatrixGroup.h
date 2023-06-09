// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_MATRIXGROUP_H_
#define SRC_POLY_MATRIXGROUP_H_

#include "GRP_GroupFct.h"
#include "Group.h"
#include "InvariantVectorFamily.h"
#include "MAT_MatrixInt.h"
#include "MAT_MatrixMod.h"
#include "MatrixGroupBasic.h"
#include "PERM_Fct.h"
#include "Timings.h"
#include "factorizations.h"
#include "two_dim_lorentzian.h"
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

//

template <typename T, typename Telt>
struct ResultGeneratePermutationGroup_Finite {
  int nbRow;
  int siz;
  std::vector<Telt> ListPermGens;
};

template <typename T, typename Telt> struct FiniteMatrixGroupHelper {
  using Treturn = ResultGeneratePermutationGroup_Finite<T, Telt>;
  int n;
  MyMatrix<T> EXTfaithful;
  std::vector<MyVector<T>> ListV;
  std::unordered_map<MyVector<T>, int> MapV;
};

//

template <typename T, typename Telt>
struct ResultGeneratePermutationGroup_FiniteIsotropic {
  int nbRow;
  int siz;
  std::vector<Telt> ListPermGens;
};

template <typename T, typename Telt> struct FiniteIsotropicMatrixGroupHelper {
  using Treturn = ResultGeneratePermutationGroup_FiniteIsotropic<T, Telt>;
  int n;
  MyMatrix<T> G;
  MyMatrix<T> EXTfaithful;
  MyVector<T> Visotrop;
  std::vector<MyVector<T>> ListV;
  std::unordered_map<MyVector<T>, int> MapV;
};

//

template <typename T, typename Telt>
struct ResultGeneratePermutationGroup_General {
  int nbRow;
  int siz;
  std::vector<MyMatrix<T>> ListMatrGens;
  std::vector<Telt> ListPermGens;
};

template <typename T, typename Telt> struct GeneralMatrixGroupHelper {
  using Treturn = ResultGeneratePermutationGroup_General<T, Telt>;
  int n;
};

//

template <typename Thelper> struct has_determining_ext {
  static const bool value = false;
};

template <typename T, typename Telt>
struct has_determining_ext<FiniteMatrixGroupHelper<T, Telt>> {
  static const bool value = true;
};

template <typename T, typename Telt>
struct has_determining_ext<FiniteIsotropicMatrixGroupHelper<T, Telt>> {
  static const bool value = true;
};

//

template <typename T, typename Telt>
GeneralMatrixGroupHelper<T, Telt>
TransformHelper(GeneralMatrixGroupHelper<T, Telt> const &helper,
                [[maybe_unused]] MyMatrix<T> const &Pmat) {
  return helper;
}

template <typename T, typename Telt>
FiniteIsotropicMatrixGroupHelper<T, Telt>
TransformHelper(FiniteIsotropicMatrixGroupHelper<T, Telt> const &helper,
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
    ListV_new.push_back(eV);
    MapV_new[eV] = i;
  }
  return {helper.n,
          std::move(G_new),
          std::move(EXTfaithful_new),
          std::move(Visotrop_new),
          std::move(ListV_new),
          std::move(MapV_new)};
}

template <typename T, typename Telt>
FiniteMatrixGroupHelper<T, Telt>
TransformHelper(FiniteMatrixGroupHelper<T, Telt> const &helper,
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

//

template <typename T, typename Thelper>
T L1normMatrixGroup(Thelper const &helper,
                    std::vector<MyMatrix<T>> const &ListMatr) {
  int n = helper.n;
  T sum = 0;
  for (auto &eMat : ListMatr) {
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        sum += T_abs(eMat(i, j));
  }
  return sum;
}

template <typename T, typename Tint, typename Thelper>
std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>>
LLLMatrixGroupReduction(Thelper const &helper,
                        std::vector<MyMatrix<T>> const &ListMatr) {
  int n = helper.n;
  MyMatrix<T> PosDefMat = IdentityMat<T>(n);
  for (auto &eMat : ListMatr) {
    if (!IsIdentity(eMat)) {
      MyMatrix<T> eProd = eMat * eMat.transpose();
      PosDefMat += eProd;
    }
  }
  LLLreduction<T, Tint> pair = LLLreducedBasis<T, Tint>(PosDefMat);
  MyMatrix<Tint> const &Pmat = pair.Pmat;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  std::vector<MyMatrix<T>> ListMatrNew;
  for (auto &eMat : ListMatr) {
    MyMatrix<T> eMatNew = Pmat_T * eMat * PmatInv_T;
    ListMatrNew.emplace_back(std::move(eMatNew));
  }
  return {std::move(ListMatrNew), Pmat};
}

template <typename T> T LinearSpace_GetDivisor(MyMatrix<T> const &TheSpace) {
  T TheDet = T_abs(DeterminantMat(TheSpace));
  T eDiv = 1;
  int n = TheSpace.rows();
  CanSolIntMat<T> eCan = ComputeCanonicalFormFastReduction(TheSpace);
  while (true) {
    bool IsOK = true;
    for (int i = 0; i < n; i++)
      if (IsOK) {
        MyVector<T> eVect = ZeroVector<T>(n);
        eVect(i) = eDiv;
        bool test = CanTestSolutionIntMat(eCan, eVect);
        if (!test)
          IsOK = false;
      }
    if (IsOK)
      return eDiv;
#ifdef SANITY_CHECK
    if (eDiv > TheDet) {
      std::cerr << "eDiv=" << eDiv << " TheDet=" << TheDet << "\n";
      std::cerr << "TheSpace=\n";
      WriteMatrix(std::cerr, TheSpace);
      std::cerr << "Clear error in LinearSpace_GetDivisor\n";
      throw TerminalException{1};
    }
#endif
    eDiv += 1;
  }
}

template <typename T>
MyMatrix<T>
MatrixIntegral_GetInvariantSpace(int const &n,
                                 std::vector<MyMatrix<T>> const &LGen) {
  std::vector<MyMatrix<T>> LGenTot;
  for (auto &eGen : LGen) {
    LGenTot.push_back(eGen);
    LGenTot.push_back(Inverse(eGen));
  }
  LGenTot.push_back(IdentityMat<T>(n));
  MyMatrix<T> TheSpace = IdentityMat<T>(n);
  T TheDet = 1;
  while (true) {
    std::vector<MyVector<T>> ConcatSpace;
    for (auto &eGen : LGenTot) {
      MyMatrix<T> TheSpaceImg = TheSpace * eGen;
      for (int i = 0; i < n; i++)
        ConcatSpace.push_back(GetMatrixRow(TheSpaceImg, i));
    }
    MyMatrix<T> NewSpace = GetZbasis(MatrixFromVectorFamily(ConcatSpace));
    T NewDet = T_abs(DeterminantMat(NewSpace));
    if (NewDet == TheDet)
      return TheSpace;
    TheSpace = NewSpace;
    TheDet = NewDet;
  }
}

// Compute Orbit of an object of type T2 under
// a group generated by elements of type T1
template <typename T1, typename T2, typename Fprod, typename Fterminate>
std::optional<std::vector<T2>>
OrbitComputation_limit(std::vector<T1> const &ListGen, T2 const &a,
                       const Fprod &f_prod, const Fterminate &f_terminate) {
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Begin of OrbitComputation_limit\n";
#endif
  std::vector<T2> TheOrbit;
  std::unordered_map<T2, uint8_t> map;
  auto fInsert = [&](T2 const &u) -> bool {
    uint8_t &pos = map[u];
    if (pos == 0) {
      pos = 1;
      TheOrbit.push_back(u);
      return f_terminate(u);
    }
    return false;
  };
  if (fInsert(a))
    return {};
  size_t pos = 0;
  while (true) {
    size_t len = TheOrbit.size();
    if (pos == len)
      break;
    for (size_t i = pos; i < len; i++)
      for (auto &eGen : ListGen) {
        T2 u = f_prod(TheOrbit[i], eGen);
        if (fInsert(u))
          return {};
      }
    pos = len;
  }
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "End of OrbitComputation_limit\n";
#endif
  return TheOrbit;
}

template <typename T1, typename T2, typename Fprod>
std::vector<T2> OrbitComputation(std::vector<T1> const &ListGen, T2 const &a,
                                 const Fprod &f_prod) {
  auto f_terminate = [&]([[maybe_unused]] T2 const &a) -> bool {
    return false;
  };
  std::optional<std::vector<T2>> opt =
      OrbitComputation_limit(ListGen, a, f_prod, f_terminate);
#ifdef SANITY_CHECK
  if (!opt) {
    std::cerr << "The opt should have been assigned\n";
    throw TerminalException{1};
  }
#endif
  return *opt;
}

template <typename T, typename Telt>
FiniteMatrixGroupHelper<T, Telt>
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

template <typename T, typename Telt>
FiniteIsotropicMatrixGroupHelper<T, Telt>
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
                                   MyMatrix<T> const &eMatr) {
#ifdef DEBUG_MATRIX_GROUP
  //  std::cerr << "Beginning of GetPermutationForFiniteMatrixGroup\n";
  //  std::cerr << "|ListV|=" << helper.ListV.size() << "\n";
  //  for (auto & eV : helper.ListV)
  //    std::cerr << "V=" << StringVectorGAP(eV) << "\n";
  //  std::cerr << "eMat=" << StringMatrixGAP(eMatr) << "\n";
  //  MyMatrix<T> eProd = eMatr * helper.G * eMatr.transpose();
  //  std::cerr << "eProd=" << StringMatrixGAP(eProd) << "\n";
  std::cerr << "eMatr=\n";
  WriteMatrix(std::cerr, eMatr);
#endif
  using Tidx = typename Telt::Tidx;
  Tidx len = helper.EXTfaithful.rows();
  std::vector<Tidx> V(len);
  for (Tidx i = 0; i < len; i++) {
    MyVector<T> const &eV = helper.ListV[i];
    MyVector<T> Vimg = eMatr.transpose() * eV;
#ifdef DEBUG_MATRIX_GROUP
    std::cerr << "i=" << i << " V=" << StringVectorGAP(eV)
              << " Vimg=" << StringVectorGAP(Vimg) << "\n";
#endif
    V[i] = helper.MapV.at(Vimg);
  }
#ifdef DEBUG_MATRIX_GROUP
  //  std::cerr << "Beginning of GetPermutationForFiniteMatrixGroup\n";
#endif
  return Telt(std::move(V));
}

template <typename T, typename Tmod>
Face GetFace(int const &nbRow, std::vector<MyVector<Tmod>> const &O,
             MyMatrix<T> const &TheSpace) {
  size_t Osiz = O.size();
  size_t siz = nbRow + Osiz;
  //  std::cerr << "GetFace : nbRow=" << nbRow << "\n";
  Face eFace(siz);
  CanSolIntMat<T> eCan = ComputeCanonicalFormFastReduction(TheSpace);
  for (size_t iO = 0; iO < Osiz; iO++) {
    MyVector<T> const &eVect = UniversalVectorConversion<T, Tmod>(O[iO]);
    bool test = CanTestSolutionIntMat(eCan, eVect);
    if (test) {
      eFace[nbRow + iO] = 1;
    }
  }
  return eFace;
}

template <typename T, typename Telt>
MyMatrix<T>
RepresentPermutationAsMatrix(FiniteMatrixGroupHelper<T, Telt> const &helper,
                             Telt const &ePerm) {
#ifdef DEBUG_MATRIX_GROUP
  //  std::cerr << "Beginning of RepresentPermutationAsMatrix for "
  //               "FiniteMatrixGroupHelper\n";
#endif
  return FindTransformation(helper.EXTfaithful, helper.EXTfaithful, ePerm);
}

template <typename T, typename Telt>
MyMatrix<T> RepresentPermutationAsMatrix(
    FiniteIsotropicMatrixGroupHelper<T, Telt> const &helper,
    Telt const &ePerm) {
#ifdef DEBUG_MATRIX_GROUP
  //  std::cerr << "Beginning of RepresentPermutationAsMatrix for "
  //               "FiniteIsotropicMatrixGroupHelper\n";
#endif
  MyMatrix<T> const &Subspace1 = helper.EXTfaithful;
  int n_rows = Subspace1.rows();
  int n_cols = Subspace1.cols();
  MyMatrix<T> Subspace2(n_rows, n_cols);
  for (int i_row = 0; i_row < n_rows; i_row++) {
    int j_row = OnPoints(i_row, ePerm);
    MyVector<T> V = GetMatrixRow(Subspace1, j_row);
    AssignMatrixRow(Subspace2, i_row, V);
  }
  std::optional<MyMatrix<T>> opt = ExtendOrthogonalIsotropicIsomorphism(
      helper.G, Subspace1, helper.G, Subspace2);
#ifdef SANITY_CHECK
  if (!opt) {
    std::cerr << "We should have opt well defined\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> const &M = *opt;
  return M;
}

template <typename T, typename Tmod, typename Telt, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value,
                               typename Thelper::Treturn>::type
MatrixIntegral_GeneratePermutationGroup(
    std::vector<MyMatrix<T>> const &ListMatrGens,
    std::vector<MyMatrix<Tmod>> const &ListMatrGensMod, Thelper const &helper,
    std::vector<MyVector<Tmod>> const &O, T const &TheMod) {
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Beginning of MatrixIntegral_GeneratePermutationGroup\n";
#endif
#ifdef TIMINGS
  MicrosecondTime time;
#endif
  using Tidx = typename Telt::Tidx;
  int Osiz = O.size();
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Osiz=" << Osiz << "\n";
#endif
  int nbRow = helper.EXTfaithful.rows();
  Tidx nbRow_tidx = nbRow;
  int siz = nbRow + Osiz;
  Telt ePermS = Telt(SortingPerm<MyVector<Tmod>, Tidx>(O));
  Tmod TheMod_mod = UniversalScalarConversion<Tmod, T>(TheMod);
  auto TheAction = [&](MyVector<Tmod> const &eClass,
                       MyMatrix<Tmod> const &eElt) -> MyVector<Tmod> {
    MyVector<Tmod> eVect = eElt.transpose() * eClass;
    return VectorMod(eVect, TheMod_mod);
  };
#ifdef TIMINGS
  std::cerr << "Timing |SortingPerm|=" << time << "\n";
#endif
  Telt ePermSinv = ~ePermS;
#ifdef TIMINGS
  std::cerr << "Timing |ePermSinv|=" << time << "\n";
#endif
  std::vector<Telt> ListPermGenProv;
  size_t nbGen = ListMatrGens.size();
  for (size_t iGen = 0; iGen < nbGen; iGen++) {
#ifdef TIMINGS
    MicrosecondTime timeB;
#endif
    MyMatrix<T> const &eMatrGen = ListMatrGens[iGen];
    MyMatrix<Tmod> const &eMatrGenMod = ListMatrGensMod[iGen];
    Telt ePermGen =
        GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(helper, eMatrGen);
    /*
    std::cerr << "siz=" << siz << " nbRow_tidx=" << int(nbRow_tidx) << "\n";
    std::cerr << "ePermGen.size()=" << int(ePermGen.size()) << "\n";
    for (Tidx i = 0; i < nbRow_tidx; i++)
      std::cerr << "i=" << int(i) << " ePermGen.at(i)=" << int(ePermGen.at(i))
    << "\n";
    */
#ifdef DEBUG_MATRIX_GROUP
    std::cerr << "iGen=" << iGen << "/" << nbGen << " ePermGen=" << ePermGen
              << "\n";
#endif
    std::vector<Tidx> v(siz);
    for (Tidx i = 0; i < nbRow_tidx; i++)
      v[i] = ePermGen.at(i);
#ifdef TIMINGS
    std::cerr << "Timing |v 1|=" << timeB << "\n";
#endif
    std::vector<MyVector<Tmod>> ListImage(Osiz);
    // That code below is shorter and it has the same speed as the above.
    // We keep the more complicate because it shows where most of the runtime
    // is: In computing Oprod.
    for (int iV = 0; iV < Osiz; iV++)
      ListImage[iV] = TheAction(O[iV], eMatrGenMod);
#ifdef TIMINGS
    std::cerr << "Timing |ListImage|=" << timeB << "\n";
#endif
    Telt ePermB = Telt(SortingPerm<MyVector<Tmod>, Tidx>(ListImage));
#ifdef TIMINGS
    std::cerr << "Timing |SortingPerm|=" << timeB << "\n";
#endif
    Telt ePermBinv = ~ePermB;
#ifdef TIMINGS
    std::cerr << "Timing |ePermBinv|=" << timeB << "\n";
#endif
    //      std::cerr << "  ePermS=" << ePermS << " ePermB=" << ePermB << "\n";
    // By the construction and above check we have
    // V1reord[i] = V1[g1.at(i)]
    // V2reord[i] = V2[g2.at(i)]
    // We have V1reord = V2reord which gets us
    // V2[i] = V1[g1 * g2^{-1}(i)]
    Telt ePermGenSelect = ePermBinv * ePermS;
#ifdef TIMINGS
    std::cerr << "Timing |ePermGenSelect|=" << timeB << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP
    //    std::cerr << "  ePermGenSelect=" << ePermGenSelect << "\n";
#endif
    for (int iO = 0; iO < Osiz; iO++) {
      int jO = ePermGenSelect.at(iO);
      v[nbRow + iO] = nbRow + jO;
    }
#ifdef TIMINGS
    std::cerr << "Timing |v 2|=" << timeB << "\n";
#endif
    Telt eNewPerm(std::move(v));
#ifdef SANITY_CHECK
    for (int iO = 0; iO < Osiz; iO++) {
      MyVector<Tmod> eVect = O[iO];
      MyVector<Tmod> eVectImg1 = TheAction(eVect, eMatrGenMod);
      size_t pos = eNewPerm.at(iO + nbRow_tidx) - nbRow_tidx;
      MyVector<Tmod> eVectImg2 = O[pos];
      if (eVectImg1 != eVectImg2) {
        std::cerr << "  Inconsistency\n";
        std::cerr << "  iGen=" << iGen << " iO=" << iO << "\n";
        std::cerr << "  eVectImg1=" << StringVectorGAP(eVectImg1) << "\n";
        std::cerr << "  eVectImg2=" << StringVectorGAP(eVectImg2) << "\n";
        throw TerminalException{1};
      }
    }
#endif
    ListPermGenProv.emplace_back(std::move(eNewPerm));
#ifdef TIMINGS
    std::cerr << "Timing |insert|=" << timeB << "\n";
#endif
  }
#ifdef TIMINGS
  std::cerr << "Timing |ListPermGenProv|=" << time << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP
  permutalib::Group<Telt, mpz_class> GRPprov(ListPermGenProv, siz);
  std::cerr << "|GRPprov|=" << GRPprov.size() << "\n";
#endif
  return {nbRow, siz, ListPermGenProv};
}

template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value,
                               std::vector<MyMatrix<T>>>::type
MatrixIntegral_PreImageSubgroup([[maybe_unused]]
                                typename Thelper::Treturn const &eret,
                                Tgroup const &eGRP, Thelper const &helper) {
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Beginning of MatrixIntegral_PreImageSubgroup\n";
#endif
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<MyMatrix<T>> ListMatrGen;
  Tidx nbRow_tidx = helper.EXTfaithful.rows();
  for (auto &eGen : eGRP.GeneratorsOfGroup()) {
    MyMatrix<T> eMatr = RepresentPermutationAsMatrix(helper, eGen);
    ListMatrGen.emplace_back(std::move(eMatr));
  }
  return ListMatrGen;
}

template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value,
                               std::vector<MyMatrix<T>>>::type
MatrixIntegral_Stabilizer([[maybe_unused]]
                          typename Thelper::Treturn const &eret,
                          Tgroup const &GRPperm, Thelper const &helper,
                          Face const &eFace) {
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Beginning of MatrixIntegral_Stabilizer\n";
#endif
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tgroup eStab = GRPperm.Stabilizer_OnSets(eFace);
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "  |eStab|=" << eStab.size() << " |eFace|=" << eFace.count()
            << "\n";
#endif
  std::vector<MyMatrix<T>> ListMatrGen;
  Tidx nbRow_tidx = helper.EXTfaithful.rows();
  for (auto &eGen : eStab.GeneratorsOfGroup()) {
    std::vector<Tidx> v(nbRow_tidx);
    for (Tidx i = 0; i < nbRow_tidx; i++)
      v[i] = OnPoints(i, eGen);
    Telt ePerm(std::move(v));
    MyMatrix<T> eMatr = RepresentPermutationAsMatrix(helper, ePerm);
    ListMatrGen.emplace_back(std::move(eMatr));
  }
  return ListMatrGen;
}

template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value,
                               std::optional<MyMatrix<T>>>::type
MatrixIntegral_RepresentativeAction([[maybe_unused]]
                                    typename Thelper::Treturn const &eret,
                                    Tgroup const &GRPperm,
                                    Thelper const &helper, Face const &eFace1,
                                    Face const &eFace2) {
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Beginning of MatrixIntegral_RepresentativeAction\n";
#endif
  using Telt = typename Tgroup::Telt;
  std::optional<Telt> opt = GRPperm.RepresentativeAction_OnSets(eFace1, eFace2);
  if (!opt) {
#ifdef DEBUG_MATRIX_GROUP
    std::cerr << "Exit while loop with proof that no equivalence exists\n";
#endif
    return {};
  }
  MyMatrix<T> eMat = RepresentPermutationAsMatrix(helper, *opt);
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "eMat=\n";
  WriteMatrix(std::cerr, eMat);
#endif
  return eMat;
}

template <typename T, typename Tmod, typename Telt, typename Thelper>
inline typename std::enable_if<
    !has_determining_ext<Thelper>::value,
    ResultGeneratePermutationGroup_General<T, Telt>>::type
MatrixIntegral_GeneratePermutationGroup(
    std::vector<MyMatrix<T>> const &ListMatrGens,
    std::vector<MyMatrix<Tmod>> const &ListMatrGensMod,
    [[maybe_unused]] Thelper const &helper,
    std::vector<MyVector<Tmod>> const &O, T const &TheMod) {
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Beginning of MatrixIntegral_GeneratePermutationGroup 2\n";
#endif
#ifdef TIMINGS
  MicrosecondTime time;
#endif
  using Tidx = typename Telt::Tidx;
  int Osiz = O.size();
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Osiz=" << Osiz << "\n";
#endif
  int siz = Osiz;
  Telt ePermS = Telt(SortingPerm<MyVector<Tmod>, Tidx>(O));
  Tmod TheMod_mod = UniversalScalarConversion<Tmod, T>(TheMod);
  auto TheAction = [&](MyVector<Tmod> const &eClass,
                       MyMatrix<Tmod> const &eElt) -> MyVector<Tmod> {
    MyVector<Tmod> eVect = eElt.transpose() * eClass;
    return VectorMod(eVect, TheMod_mod);
  };
#ifdef TIMINGS
  std::cerr << "Timing |SortingPerm|=" << time << "\n";
#endif
  Telt ePermSinv = ~ePermS;
#ifdef TIMINGS
  std::cerr << "Timing |ePermSinv|=" << time << "\n";
#endif
  std::vector<Telt> ListPermGenProv;
  size_t nbGen = ListMatrGens.size();
  for (size_t iGen = 0; iGen < nbGen; iGen++) {
#ifdef TIMINGS
    MicrosecondTime timeB;
#endif
    //    MyMatrix<T> const& eMatrGen=ListMatrGens[iGen];
    MyMatrix<Tmod> const &eMatrGenMod = ListMatrGensMod[iGen];
#ifdef DEBUG_MATRIX_GROUP
    std::cerr << "iGen=" << iGen << "/" << nbGen << "\n";
#endif
    std::vector<Tidx> v(siz);
#ifdef TIMINGS
    std::cerr << "Timing |v 1|=" << timeB << "\n";
#endif
    std::vector<MyVector<Tmod>> ListImage(Osiz);
    // That code below is shorter and it has the same speed as the above.
    // We keep the more complicate because it shows where most of the runtime
    // is: In computing Oprod.
    for (int iV = 0; iV < Osiz; iV++)
      ListImage[iV] = TheAction(O[iV], eMatrGenMod);
#ifdef TIMINGS
    std::cerr << "Timing |ListImage|=" << timeB << "\n";
#endif
    Telt ePermB = Telt(SortingPerm<MyVector<Tmod>, Tidx>(ListImage));
#ifdef TIMINGS
    std::cerr << "Timing |SortingPerm|=" << timeB << "\n";
#endif
    Telt ePermBinv = ~ePermB;
#ifdef TIMINGS
    std::cerr << "Timing |ePermBinv|=" << timeB << "\n";
#endif
    //      std::cerr << "  ePermS=" << ePermS << " ePermB=" << ePermB << "\n";
    // By the construction and above check we have
    // V1reord[i] = V1[g1.at(i)]
    // V2reord[i] = V2[g2.at(i)]
    // We have V1reord = V2reord which gets us
    // V2[i] = V1[g1 * g2^{-1}(i)]
    Telt ePermGenSelect = ePermBinv * ePermS;
#ifdef TIMINGS
    std::cerr << "Timing |ePermGenSelect|=" << timeB << "\n";
#endif
    ListPermGenProv.emplace_back(std::move(ePermGenSelect));
#ifdef TIMINGS
    std::cerr << "Timing |insert|=" << timeB << "\n";
#endif
  }
#ifdef TIMINGS
  std::cerr << "Timing |ListPermGenProv|=" << time << "\n";
#endif
  return {0, siz, ListMatrGens, std::move(ListPermGenProv)};
}

template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<!has_determining_ext<Thelper>::value,
                               std::vector<MyMatrix<T>>>::type
MatrixIntegral_PreImageSubgroup(typename Thelper::Treturn const &eret,
                                Tgroup const &eGRP, Thelper const &helper) {
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Beginning of MatrixIntegral_PreImageSubgroup\n";
#endif
  MyMatrix<T> id_matr = IdentityMat<T>(helper.n);
  return permutalib::PreImageSubgroup<Tgroup, MyMatrix<T>>(
      eret.ListMatrGens, eret.ListPermGens, id_matr, eGRP);
}

template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<!has_determining_ext<Thelper>::value,
                               std::vector<MyMatrix<T>>>::type
MatrixIntegral_Stabilizer(typename Thelper::Treturn const &eret,
                          [[maybe_unused]] Tgroup const &GRPperm,
                          Thelper const &helper, Face const &eFace) {
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Beginning of MatrixIntegral_Stabilizer 2\n";
#endif
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
  MyMatrix<T> id_matr = IdentityMat<T>(helper.n);
  return permutalib::StabilizerMatrixPermSubset<Telt, MyMatrix<T>, Tint>(
      eret.ListMatrGens, eret.ListPermGens, id_matr, eFace);
}

template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<!has_determining_ext<Thelper>::value,
                               std::optional<MyMatrix<T>>>::type
MatrixIntegral_RepresentativeAction(typename Thelper::Treturn const &eret,
                                    [[maybe_unused]] Tgroup const &GRPperm,
                                    Thelper const &helper, Face const &eFace1,
                                    Face const &eFace2) {
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Beginning of MatrixIntegral_RepresentativeAction 2\n";
#endif
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
  MyMatrix<T> id_matr = IdentityMat<T>(helper.n);
  std::optional<MyMatrix<T>> opt =
      permutalib::RepresentativeActionMatrixPermSubset<Telt, MyMatrix<T>, Tint>(
          eret.ListMatrGens, eret.ListPermGens, id_matr, eFace1, eFace2);
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Ending of MatrixIntegral_RepresentativeAction 2\n";
#endif
  return opt;
}

/*
  Direct computation of orbits.
  First level of optional is for termination or not.
  Second level is for whether we find an equivalence or not.
 */
template <typename T, typename Fterminate>
std::optional<std::optional<MyMatrix<T>>>
DirectSpaceOrbit_Equivalence(std::vector<MyMatrix<T>> const &ListMatrGen,
                             MyMatrix<T> const &eSpace1,
                             MyMatrix<T> const &eSpace2, T const &TheMod,
                             Fterminate const &f_terminate) {
  int n = eSpace1.rows();
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
  // Here Tpair is <Space,Repr>
  using Tpair = std::pair<MyMatrix<T>, MyMatrix<T>>;
  std::vector<Tpair> ListPair;
  ListPair.push_back({eSpace1, IdentityMat<T>(n)});
  if (f_terminate(eSpace1))
    return {};
  size_t pos = 0;
  while (true) {
    size_t len = ListPair.size();
    if (pos == len)
      break;
    for (size_t idx = pos; idx < len; idx++) {
      Tpair const &ePair = ListPair[idx];
      for (auto &eMatrGen : ListMatrGen) {
        MyMatrix<T> eSpaceImg = ePair.first * eMatrGen;
        MyMatrix<T> eReprImg = ePair.second * eMatrGen;
        //
        MyMatrix<T> eSpaceMod = Concatenate(ePair.first, ModSpace);
        CanSolIntMat<T> eCan = ComputeCanonicalFormFastReduction(eSpaceMod);
        auto IsEqual = [&](MyMatrix<T> const &fSpace) -> bool {
          for (int i = 0; i < n; i++) {
            MyVector<T> V = GetMatrixRow(fSpace, i);
            bool test = CanTestSolutionIntMat(eCan, V);
            if (!test)
              return false;
          }
          return true;
        };
        if (IsEqual(eSpace2)) {
          std::optional<MyMatrix<T>> opt = eReprImg;
          return opt;
        }
        auto fInsert = [&](Tpair const &ePair) -> bool {
          for (auto &fPair : ListPair)
            if (IsEqual(fPair.first))
              return false;
          ListPair.push_back(ePair);
          return f_terminate(ePair.first);
        };
        if (fInsert(eSpaceImg))
          return {};
      }
    }
    std::cerr << "pos=" << pos << " len=" << len << "\n";
    pos = len;
  }
  std::optional<MyMatrix<T>> opt;
  return opt;
}

template <typename T, typename Fterminate>
std::optional<std::vector<MyMatrix<T>>>
DirectSpaceOrbit_Stabilizer(std::vector<MyMatrix<T>> const &ListMatrGen,
                            MyMatrix<T> const &eSpace, T const &TheMod,
                            Fterminate const &f_terminate) {
  int n = eSpace.rows();
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
  // Tpair is a pair <Space, Repr>
  using Tpair = std::pair<MyMatrix<T>, MyMatrix<T>>;
  std::vector<Tpair> ListPair;
  ListPair.push_back({eSpace, IdentityMat<T>(n)});
  if (f_terminate(eSpace))
    return {};
  size_t pos = 0;
  while (true) {
    size_t len = ListPair.size();
    if (pos == len)
      break;
    for (size_t idx = pos; idx < len; idx++) {
      Tpair const &ePair = ListPair[idx];
      for (auto &eMatrGen : ListMatrGen) {
        MyMatrix<T> eSpaceImg = ePair.first * eMatrGen;
        MyMatrix<T> eReprImg = ePair.second * eMatrGen;
        //
        MyMatrix<T> eSpaceMod = Concatenate(ePair.first, ModSpace);
        CanSolIntMat<T> eCan = ComputeCanonicalFormFastReduction(eSpaceMod);
        auto IsEqual = [&](MyMatrix<T> const &fSpace) -> bool {
          for (int i = 0; i < n; i++) {
            MyVector<T> V = GetMatrixRow(fSpace, i);
            bool test = CanTestSolutionIntMat(eCan, V);
            if (!test)
              return false;
          }
          return true;
        };
        auto fInsert = [&](Tpair const &ePair) -> bool {
          for (auto &fPair : ListPair)
            if (IsEqual(fPair.first))
              return false;
          ListPair.push_back(ePair);
          return f_terminate(ePair.first);
        };
        if (fInsert(eSpaceImg))
          return {};
      }
    }
    std::cerr << "pos=" << pos << " len=" << len << "\n";
    pos = len;
  }
  //
  // Orbit is fine, now computing the stabilizer by using the Schreier lemma.
  //
  std::unordered_set<MyMatrix<T>> SetGen;
  size_t nPair = ListPair.size();
  for (size_t iPair = 0; iPair < nPair; iPair++) {
    Tpair const &ePair = ListPair[iPair];
    for (auto &eMatrGen : ListMatrGen) {
      MyMatrix<T> eSpaceImg = ePair.first * eMatrGen;
      MyMatrix<T> eSpaceImgMod = Concatenate(eSpaceImg, ModSpace);
      CanSolIntMat<T> eCan = ComputeCanonicalFormFastReduction(eSpaceImgMod);
      auto f_equal = [&](MyMatrix<T> const &fSpace) -> bool {
        for (int i = 0; i < n; i++) {
          MyVector<T> V = GetMatrixRow(fSpace, i);
          bool test = CanTestSolutionIntMat(eCan, V);
          if (!test)
            return false;
        }
        return true;
      };
      auto f_insert = [&]() -> void {
        for (size_t jPair = 0; jPair < nPair; jPair++) {
          if (f_equal(ListPair[jPair].first)) {
            MyMatrix<T> eGenMatr_new =
                ePair.second * eMatrGen * Inverse(ListPair[jPair].second);
            if (!IsIdentity(eGenMatr_new))
              SetGen.insert(eGenMatr_new);
          }
        }
      };
      f_insert();
    }
  }
  std::vector<MyMatrix<T>> ListGen;
  for (auto &eGen : SetGen)
    ListGen.push_back(eGen);
  return ListGen;
}

/*
There are several challenges for this implementation.
---We need an invariant vector family which we may obtain
   from a positive definite forms.
---We have already the one from the Lorentzian family or
   from the Q_inv. But it is unlikely to be of full rank.
---If the configuration is of full rank we win. If not, it is
   coming from an isotropic vector.
---From the fact that the group of interest is finite, we know
   that there exist a positive definite quadratic form that need
   to be found.
---From the Qinv, we can get a quadratic form preserving the
   n-1 dimensional space. That partial quadratic form, then has to be
   extended.
   In essence we ask that any isometry preserving the n-1 dimensional
   space also would preserve that one.
   We need to obtain that form from computation not involving the
   group.
---If we limit ourselves to the case of just preserving v^{perp} then
   we do not have finiteness. Let us take the form -x0^2 + x1^2 + x2^2
   The space v^{perp} is of dimension 2 and the restricted quadratic
   form is a_u^2, that is positive semidefinite but not definite. As a
   conclusion the space of isometries is
   (a_u,a_v) -> (a_u + C a_v, a_v).
   It is an isometry for all C because v is isotropic. All those isometries
   are extendible to full isometries of the 3-dim space by the known lemma.
---Under this infinite group there is a single orbit of isotropic vectors
   different from v.
---Now we have the additional finiteness requirement. That additional
   requirement is encapsulated into the preservation of the Qinv form.
   The vector v of v^{perp} is preserved. The orthogonal H of v in v^{perp}
   for Qinv is preserved as well.
---Now, we need to find another isotropic vector. The previous argument
   suggests that there is a single vector to be found.
   We can work with H^{perp} for the Lorentzian scalar product. It is
   a two dimensional space of signature (1,1). It has two isotropic
   vectors in it.
   One vector would be v, the other the one we are looking for.
---From this we can build a full dimensional vector system and nice
   quadratic form. And so with the Shortest vector problem, we can
   get vectors.
---From ListMatr, we can compute the space of invariant quadratic forms
   and this can get us good efficient signatures.
 */

template <typename T, typename Telt>
std::vector<MyMatrix<T>>
GetListQuadraticForms(FiniteIsotropicMatrixGroupHelper<T, Telt> const &helper) {
  int n = helper.n;
  int n_vect = helper.EXTfaithful.rows();
  MyMatrix<T> BasisSp = RowReduction(helper.EXTfaithful);
  std::optional<MyMatrix<T>> opt1 =
      ListSolutionMat(BasisSp, helper.EXTfaithful);
#ifdef SANITY_CHECK
  if (!opt1) {
    std::cerr << "opt1 : Failed to solution the linear systems\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> const &EXTfaithful_red = *opt1;
  std::optional<MyVector<T>> opt2 = SolutionMat(BasisSp, helper.Visotrop);
#ifdef SANITY_CHECK
  if (!opt2) {
    std::cerr << "opt2 : Failed to solution the linear systems\n";
    throw TerminalException{1};
  }
#endif
  MyVector<T> const &Visotrop_red = *opt2;
  MyMatrix<T> Qinv = GetQmatrix(EXTfaithful_red);
  MyMatrix<T> eProd1 = MatrixFromVectorFamily(Qinv * Visotrop_red);
  MyMatrix<T> PosDefSpace = NullspaceMat(eProd1);
  MyMatrix<T> eProd2 = helper.G * PosDefSpace.transpose();
  MyMatrix<T> Sign11space = NullspaceMat(eProd2);
  MyMatrix<T> G11 = Sign11space * helper.G * Sign11space.transpose();
  auto get_second_isotropic = [&]() -> MyVector<T> {
    std::vector<MyVector<T>> LVect = GetBasisIsotropicVectors(G11);
    for (auto &V : LVect) {
      if (!IsVectorMultiple(V, helper.Visotrop)) {
        MyVector<T> Vret = Sign11space.transpose() * V;
        return Vret;
      }
    }
    std::cerr << "Find the other vector\n";
    throw TerminalException{1};
  };
  MyVector<T> VisotropSecond = get_second_isotropic();

  MyMatrix<T> FullBasis(n, n);
  for (int i_vect = 0; i_vect < n_vect; i_vect++) {
    for (int i = 0; i < n; i++)
      FullBasis(i_vect, i) = BasisSp(i_vect, i);
  }
  for (int i = 0; i < n; i++)
    FullBasis(n - 1, i) = VisotropSecond(i);
  //
  MyMatrix<T> prePosDefMat1 = ZeroMatrix<T>(n, n);
  for (int i = 0; i < n - 1; i++)
    for (int j = 0; j < n - 1; j++)
      prePosDefMat1(i, j) = Qinv(i, j);
  MyMatrix<T> prePosDefMat2 = ZeroMatrix<T>(n, n);
  prePosDefMat2(n - 1, n - 1) = 1;
  //
  MyMatrix<T> InvFullBasis = Inverse(FullBasis);
  MyMatrix<T> PosDefMat1 =
      InvFullBasis * prePosDefMat1 * InvFullBasis.transpose();
  MyMatrix<T> PosDefMat2 =
      InvFullBasis * prePosDefMat2 * InvFullBasis.transpose();

  return {PosDefMat1, PosDefMat2};
}

template <typename T, typename Telt>
std::vector<MyMatrix<T>>
GetListQuadraticForms(FiniteMatrixGroupHelper<T, Telt> const &helper) {
  MyMatrix<T> Qinv = GetQmatrix(helper.EXTfaithful);
  return {Qinv};
}

/*
template <typename T, typename Tmod, typename Tgroup, typename Thelper>
std::optional<std::vector<MyMatrix<T>>>
PleskenSouvignier_Subspace_Stabilizer(std::vector<MyMatrix<T>> const &ListMatr,
                                      Thelper const &helper) {
  using Tint = typename underlying_ring<T>::ring_type;
  using Tidx = typename Tgroup::Telt::Tidx;
  int n = helper.n;
  std::vector<MyMatrix<T>> BasisSymmMat = BasisInvariantForm(n, ListMatr);
  std::vector<MyMatrix<T>> ListPosDef = GetListQuadraticForms(helper);
  std::vector<MyVector<T>> ListVect;
  std::vector<T> Vdiag;
  T idx = 0;
  for (auto &ePosDef : ListPosDef) {
    if (ePosDef.rows() > 1) {
      MyMatrix<T> eNSP = NullspaceIntMat(ePosDef);
      MyMatrix<Tint> eNSP_i = UniversalMatrixConversion<Tint, T>(eNSP);
      CanSolIntMat<Tint> eCan = ComputeCanonicalFormFastReduction(eNSP);
      MyMatrix<T> eGnsp = eNSP * ePosDef * eNSP.transpose();
      std::vector<MyMatrix<Tint>> ListMatrRed;
      for (auto &eMatr : ListMatr) {
        MyMatrix<Tint> eMatr_i = UniversalMatrixConversion<Tint, T>(eMatr);
        MyMatrix<Tint> eProd = eNSP_i * eMatr_i;
        std::optional<MyMatrix<Tint>> opt = CanSolutionIntMatMat(eCan, eProd);
#ifdef SANITY_CHECK
        if (!opt) {
          std::cerr << "The NSP space is not preserved\n";
          throw TerminalException{1};
        }
#endif
        ListMatrRed.push_back(*opt);
      }
      MyMatrix<Tint> SHV_e =
          ExtractInvariantBreakingVectorFamily(eGnsp, ListMatrRed);
      MyMatrix<Tint> SHVbreak = SHV_e * eNSP_i;
      std::vector<T> Vdiag(SHV_e.rows(), idx);
      // Getting a full dimensional family from the other vectors.
      T jdx = 0;
      for (auto &fPosDef : ListPosDef) {
        if (idx != jdx) {
          MyMatrix<T> fNSP = NullspaceIntMat(fPosDef);
          MyMatrix<Tint> fNSP_i = UniversalMatrixConversion<Tint, T>(fNSP);
          MyMatrix<T> fGnsp = fNSP * fPosDef * fNSP.transpose();
          MyMatrix<Tint> SHV_f =
              ExtractInvariantVectorFamilyFullRank<T, Tint>(fGnsp);
          MyMatrix<Tint> fProd = SHV_f * fNSP_i;
          SHVbreak = Concatenation(SHVbreak, fProd);
          std::vector<T> V(SHV_f.rows(), jdx);
          Vdiag.insert(Vdiag.end(), V.begin(), V.end());
        }
        jdx++;
      }
      //
#ifdef SANITY_CHECK
      if (RankMat(SHVbreak) != n) {
        std::cerr << "The matrix SHVbreak is not of full rank\n";
        throw TerminalException{1};
      }
#endif
      //
      const bool use_scheme = true;
      using Tfield = typename overlying_field<T>::field_type;
      std::vector<std::vector<Tidx>> ListListIdx =
          GetListGenAutomorphism_ListMat_Vdiag<T,Tfield,Tidx,use_scheme>(SHVbreak, BasisSymmMat, Vdiag);
      std::vector<MyMatrix<T>> NewListMatr;
      for (auto &eListIdx : ListListIdx) {
        std::optional<MyMatrix<T>> opt =
            FindMatrixTransformationTest(SHVbreak, SHVbreak, eListIdx);
#ifdef SANITY_CHECK
        if (!opt) {
          std::cerr << "Failed to represent the permutation\n";
          throw TerminalException{1};
        }
#endif
        NewListMatr.push_back(*opt);
      }
    }
    idx++;
  }
  std::cerr << "Failed to find a breaking invariant family\n";
  throw TerminalException{1};
}
*/



template <typename T, typename Tmod, typename Tgroup, typename Thelper>
inline typename std::enable_if<!has_determining_ext<Thelper>::value,
                               std::optional<std::vector<MyVector<Tmod>>>>::type
FindingSmallOrbit([[maybe_unused]] std::vector<MyMatrix<T>> const &ListMatrGen,
                  std::vector<MyMatrix<Tmod>> const &ListMatrGenMod,
                  [[maybe_unused]] MyMatrix<T> const &TheSpace, T const &TheMod,
                  MyVector<T> const &x,
                  [[maybe_unused]] Thelper const &helper) {
  // No determining EXT, hard to find clever ideas.
  MyVector<Tmod> x_mod = ModuloReductionVector<T, Tmod>(x, TheMod);
  Tmod TheMod_mod = UniversalScalarConversion<Tmod, T>(TheMod);
  auto f_prod = [&](MyVector<Tmod> const &eClass,
                    MyMatrix<Tmod> const &eElt) -> MyVector<Tmod> {
    MyVector<Tmod> eVect = eElt.transpose() * eClass;
    return VectorMod(eVect, TheMod_mod);
  };
  return OrbitComputation(ListMatrGenMod, x_mod, f_prod);
}

template <typename T, typename Tmod, typename Tgroup, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value,
                               std::optional<std::vector<MyVector<Tmod>>>>::type
FindingSmallOrbit(std::vector<MyMatrix<T>> const &ListMatrGen,
                  std::vector<MyMatrix<Tmod>> const &ListMatrGenMod,
                  MyMatrix<T> const &TheSpace, T const &TheMod,
                  MyVector<T> const &a, Thelper const &helper) {
  using Telt = typename Tgroup::Telt;
  int n = TheSpace.rows();
  // The critical number for the computation
  size_t n_limit = 60000;
  auto test_adequateness =
      [&](MyVector<T> const &x) -> std::optional<std::vector<MyVector<Tmod>>> {
    MyVector<Tmod> x_mod = ModuloReductionVector<T, Tmod>(x, TheMod);
    MyVector<T> x_modT = UniversalVectorConversion<T, Tmod>(x_mod);
    //    std::cerr << "x_mod = " << StringVectorGAP(x_modT) << "\n";
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
    return OrbitComputation_limit(ListMatrGenMod, x_mod, f_prod, f_terminate);
  };
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
  MyMatrix<T> TheSpaceMod = Concatenate(TheSpace, ModSpace);
  CanSolIntMat<T> eCan = ComputeCanonicalFormFastReduction(TheSpaceMod);
  auto IsStabilized = [&](MyVector<T> const &V) -> bool {
    for (auto &eMatrGen : ListMatrGen) {
      MyVector<T> Vimg = eMatrGen.transpose() * V;
      bool test = CanTestSolutionIntMat(eCan, Vimg);
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
    Telt ePermGen =
        GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(helper, eMatrGen);
    ListPermGen.push_back(ePermGen);
  }
  size_t len = helper.EXTfaithful.rows();
  Telt id_perm(len);
  Tgroup GRP(ListPermGen, id_perm);
  std::vector<Tgroup> ListGroup = GRP.GetAscendingChain();
  size_t len_group = ListGroup.size();
  std::cerr << "len_group=" << len_group << " |GRP|=" << GRP.size() << "\n";
  for (size_t iGroup = 0; iGroup < len_group; iGroup++) {
    //    std::cerr << "iGroup=" << iGroup << " |eGRP|=" <<
    //    ListGroup[iGroup].size() << "\n";
  }
  auto try_basis = [&](MyMatrix<T> const &TheBasis)
      -> std::optional<std::vector<MyVector<Tmod>>> {
    for (int i_row = 0; i_row < TheBasis.rows(); i_row++) {
      MyVector<T> V = GetMatrixRow(TheBasis, i_row);
      if (!IsStabilized(V)) {
        std::optional<std::vector<MyVector<Tmod>>> opt = test_adequateness(V);
        if (opt) {
          std::cerr << "|*opt|=" << opt->size() << "\n";
          return *opt;
        }
        std::cerr << "Too large size at i_row=" << i_row << "\n";
      }
    }
    return {};
  };
  for (size_t iGroup = 0; iGroup < len_group; iGroup++) {
    size_t jGroup = len_group - 1 - iGroup;
    Tgroup const &fGRP = ListGroup[jGroup];
    std::cerr << "iGroup=" << iGroup << " |fGRP|=" << fGRP.size() << "\n";
    std::vector<MyMatrix<T>> LMatr;
    for (auto &eGen : fGRP.GeneratorsOfGroup()) {
      MyMatrix<T> eMat = RepresentPermutationAsMatrix(helper, eGen);
      LMatr.push_back(eMat);
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

// The space must be defining a finite index subgroup of T^n
template <typename T, typename Tmod, typename Tgroup, typename Thelper>
std::vector<MyMatrix<T>>
LinearSpace_ModStabilizer_Tmod(std::vector<MyMatrix<T>> const &ListMatr,
                               Thelper const &helper,
                               MyMatrix<T> const &TheSpace, T const &TheMod) {
  using Telt = typename Tgroup::Telt;
  using Treturn = typename Thelper::Treturn;
  int n = helper.n;
#ifdef DEBUG_MATRIX_GROUP
  T TotSize = 1;
  for (int i = 0; i < n; i++)
    TotSize *= TheMod;
  std::cerr << "TheMod=" << TheMod << "  n=" << n << " TotSize=" << TotSize
            << "\n";

  //  std::cerr << "TheSpace=\n";
  //  WriteMatrix(std::cerr, TheSpace);
#endif
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
  MyMatrix<T> TheSpaceMod = Concatenate(TheSpace, ModSpace);
  CanSolIntMat<T> eCan = ComputeCanonicalFormFastReduction(TheSpaceMod);
  //  Tmod TheMod_mod = UniversalScalarConversion<Tmod, T>(TheMod);
  /*
  auto TheAction = [&](MyVector<Tmod> const &eClass,
                       MyMatrix<Tmod> const &eElt) -> MyVector<Tmod> {
    MyVector<Tmod> eVect = eElt.transpose() * eClass;
    return VectorMod(eVect, TheMod_mod);
  };
  */
  // This is the part of the enumeration where we have problems.
  // We have too many vectors to consider whih sinks the algorithm.
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
  auto IsStabilizing = [&](std::vector<MyMatrix<T>> const &ListMatrInp)
      -> std::optional<MyVector<T>> {
    for (int i = 0; i < n; i++) {
      MyVector<T> eVect = GetMatrixRow(TheSpace, i);
      for (auto &eGen : ListMatrInp) {
        MyVector<T> eVectG = eGen.transpose() * eVect;
        bool test = CanTestSolutionIntMat(eCan, eVectG);
        if (!test) {
#ifdef DEBUG_MATRIX_GROUP
          std::cerr << "i=" << i << "  eVect=" << StringVectorGAP(eVect)
                    << "\n";
#endif
          return eVect;
        }
      }
    }
    return {};
  };
  std::vector<MyMatrix<T>> ListMatrRet = ListMatr;
  std::vector<MyMatrix<Tmod>> ListMatrRetMod =
      ModuloReductionStdVectorMatrix<T, Tmod>(ListMatrRet, TheMod);
  while (true) {
#ifdef TIMINGS
    MicrosecondTime time;
#endif
    std::optional<MyVector<T>> opt = IsStabilizing(ListMatrRet);
#ifdef TIMINGS
    std::cerr << "Timing |IsStabilizing|=" << time << "\n";
#endif
    if (!opt) {
#ifdef DEBUG_MATRIX_GROUP
      std::cerr << "Exiting the loop\n";
#endif
      break;
    }
    const MyVector<T> &V = *opt;
    std::optional<std::vector<MyVector<Tmod>>> opt_fso =
        FindingSmallOrbit<T, Tmod, Tgroup, Thelper>(
            ListMatrRet, ListMatrRetMod, TheSpace, TheMod, V, helper);
#ifdef SANITY_CHECK
    if (!opt_fso) {
      std::cerr << "Failed to find some entry\n";
      throw TerminalException{1};
    }
#endif
    std::vector<MyVector<Tmod>> const &O = *opt_fso;
    //    MyVector<Tmod> V_mod = ModuloReductionVector<T,Tmod>(V, TheMod);
    //    std::vector<MyVector<Tmod>> O =
    //      OrbitComputation(ListMatrRetMod, V_mod, TheAction);
#ifdef DEBUG_MATRIX_GROUP
    std::cerr << "Orbit size |O|=" << O.size() << "\n";
#endif
#ifdef TIMINGS
    std::cerr << "Timing |OrbitComputation|=" << time << "\n";
#endif
    Treturn eret =
        MatrixIntegral_GeneratePermutationGroup<T, Tmod, Telt, Thelper>(
            ListMatrRet, ListMatrRetMod, helper, O, TheMod);

    Tgroup GRPwork(eret.ListPermGens, eret.siz);
    Face eFace = GetFace<T, Tmod>(eret.nbRow, O, TheSpaceMod);
#ifdef MATRIX_GROUP_DIAGNOSTICS
    std::cerr << "ModStabilizer TheMod=" << TheMod << " |O|=" << O.size()
              << " |GRPwork|=" << GRPwork.size() << " |eFace|=" << eFace.count()
              << "\n";
#endif

    ListMatrRet = MatrixIntegral_Stabilizer<T, Tgroup, Thelper>(eret, GRPwork,
                                                                helper, eFace);
    ListMatrRetMod =
        ModuloReductionStdVectorMatrix<T, Tmod>(ListMatrRet, TheMod);
  }
  return ListMatrRet;
}

template <typename T, typename Tgroup, typename Thelper>
std::vector<MyMatrix<T>>
LinearSpace_ModStabilizer(std::vector<MyMatrix<T>> const &ListMatr,
                          Thelper const &helper, MyMatrix<T> const &TheSpace,
                          T const &TheMod) {
  T max_size = (TheMod - 1) * (TheMod - 1) * TheSpace.rows();
  if (max_size < T(std::numeric_limits<uint8_t>::max())) {
    return LinearSpace_ModStabilizer_Tmod<T, uint8_t, Tgroup, Thelper>(
        ListMatr, helper, TheSpace, TheMod);
  }
  if (max_size < T(std::numeric_limits<uint16_t>::max())) {
    return LinearSpace_ModStabilizer_Tmod<T, uint16_t, Tgroup, Thelper>(
        ListMatr, helper, TheSpace, TheMod);
  }
  if (max_size < T(std::numeric_limits<uint32_t>::max())) {
    return LinearSpace_ModStabilizer_Tmod<T, uint32_t, Tgroup, Thelper>(
        ListMatr, helper, TheSpace, TheMod);
  }
  std::cerr << "Failed to find a matching arithmetic type. Quite unlikely "
               "objectively\n";
  throw TerminalException{1};
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
    MyMatrix<T> const &TheSpace2, T const &TheMod) {
  using Telt = typename Tgroup::Telt;
  using Treturn = typename Thelper::Treturn;
  int n = TheSpace1.rows();
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "------------------------------------------------------\n";
  std::cerr << "NeedStabilizer=" << NeedStabilizer << "\n";
  std::cerr << "LinearSpace_ModEquivalence_Tmod, TheMod=" << TheMod << "\n";
  std::cerr << "TheSpace1=\n";
  WriteMatrix(std::cerr, TheSpace1);
  std::cerr << "TheSpace2=\n";
  WriteMatrix(std::cerr, TheSpace2);
  std::cerr << "det(TheSpace1)=" << DeterminantMat(TheSpace1)
            << " det(TheSpace2)=" << DeterminantMat(TheSpace2) << "\n";
#endif
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
  MyMatrix<T> TheSpace1Mod = Concatenate(TheSpace1, ModSpace);
  MyMatrix<T> TheSpace2Mod = Concatenate(TheSpace2, ModSpace);
  std::vector<MyMatrix<T>> ListMatrRet = ListMatr;
  std::vector<MyMatrix<Tmod>> ListMatrRetMod =
      ModuloReductionStdVectorMatrix<T, Tmod>(ListMatrRet, TheMod);
  MyMatrix<T> eElt = IdentityMat<T>(n);
  Tmod TheMod_mod = UniversalScalarConversion<Tmod, T>(TheMod);
  auto TheAction = [&](MyVector<Tmod> const &eClass,
                       MyMatrix<Tmod> const &eElt) -> MyVector<Tmod> {
    MyVector<Tmod> eVect = eElt.transpose() * eClass;
    return VectorMod(eVect, TheMod_mod);
  };
  CanSolIntMat<T> eCan = ComputeCanonicalFormFastReduction(TheSpace2Mod);
  auto IsEquiv =
      [&](MyMatrix<T> const &eEquiv) -> std::optional<MyVector<Tmod>> {
    MyMatrix<T> TheSpace1img = TheSpace1 * eEquiv;
    for (int i = 0; i < n; i++) {
      MyVector<T> eVect = GetMatrixRow(TheSpace1img, i);
      bool test = CanTestSolutionIntMat(eCan, eVect);
      if (!test) {
#ifdef DEBUG_MATRIX_GROUP
        std::cerr << "   i=" << i << " eVect=" << StringVectorGAP(eVect)
                  << "\n";
        std::cerr << "   eEquiv=\n";
        WriteMatrix(std::cerr, eEquiv);
#endif
        return ModuloReductionVector<T, Tmod>(eVect, TheMod);
      }
    }
    return {};
  };
  auto IsStabilizing = [&](std::vector<MyMatrix<T>> const &ListMat)
      -> std::optional<MyVector<Tmod>> {
    if (!NeedStabilizer)
      return {};
    for (auto &eGen : ListMat) {
      MyMatrix<T> TheSpace2img = TheSpace2 * eGen;
      for (int i = 0; i < n; i++) {
        MyVector<T> eVect = GetMatrixRow(TheSpace2img, i);
        bool test = CanTestSolutionIntMat(eCan, eVect);
        if (!test)
          return ModuloReductionVector<T, Tmod>(eVect, TheMod);
      }
    }
    return {};
  };
  while (true) {
    std::optional<MyVector<Tmod>> test1 = IsEquiv(eElt);
    std::optional<MyVector<Tmod>> test2 = IsStabilizing(ListMatrRet);
    if (!test1 && !test2) {
#ifdef DEBUG_MATRIX_GROUP
      std::cerr << "eElt and GRPwork are correct. Exiting\n";
#endif
      ResultTestModEquivalence<T> res{ListMatrRet, eElt};
      return res;
    }
    if (test1) {
      MyVector<Tmod> const &V = *test1;
      std::vector<MyVector<Tmod>> O =
          OrbitComputation(ListMatrRetMod, V, TheAction);
#ifdef DEBUG_MATRIX_GROUP
      std::cerr << "|O|=" << O.size() << "\n";
#endif

      Treturn eret =
          MatrixIntegral_GeneratePermutationGroup<T, Tmod, Telt, Thelper>(
              ListMatrRet, ListMatrRetMod, helper, O, TheMod);
#ifdef DEBUG_MATRIX_GROUP
      if constexpr (has_determining_ext<Thelper>::value) {
        for (size_t iGen = 0; iGen < ListMatrRet.size(); iGen++) {
          Telt ePerm = eret.ListPermGens[iGen];
          std::cerr << "ePerm=" << ePerm;
          MyMatrix<T> eMatr = ListMatrRet[iGen];
          MyMatrix<T> M2 = RepresentPermutationAsMatrix(helper, ePerm);
          if (eMatr != M2) {
            std::cerr << " INCORRECT\n";
          } else {
            std::cerr << " correct\n";
          }
        }
      }
#endif
      Tgroup GRPperm(eret.ListPermGens, eret.siz);

      MyMatrix<T> TheSpace1work = TheSpace1 * eElt;
      MyMatrix<T> TheSpace1workMod = Concatenate(TheSpace1work, ModSpace);
#ifdef DEBUG_MATRIX_GROUP
      std::cerr << "eElt=\n";
      WriteMatrix(std::cerr, eElt);
#endif
      Face eFace1 = GetFace<T, Tmod>(eret.nbRow, O, TheSpace1workMod);
      Face eFace2 = GetFace<T, Tmod>(eret.nbRow, O, TheSpace2Mod);
#ifdef DEBUG_MATRIX_GROUP
      std::cerr << "nbRow=" << eret.nbRow << " eFace1=" << StringFace(eFace1)
                << " eFace2=" << StringFace(eFace2) << "\n";
#endif
#ifdef SANITY_CHECK
      if (eFace1.count() == 0 && eFace2.count() == 0) {
        std::cerr << "Error in LinearSpace_ModEquivalence_Tmod. |eFace1| = "
                     "|eFace2| = 0\n";
        std::cerr << "Clear bug\n";
        throw TerminalException{1};
      }
#endif
#ifdef MATRIX_GROUP_DIAGNOSTICS
      std::cerr << "ModEquivalence 1 TheMod=" << TheMod << " |O|=" << O.size()
                << " |GRPperm|=" << GRPperm.size()
                << " |eFace1|=" << eFace1.count()
                << " |eFace2|=" << eFace2.count() << "\n";
#endif
      std::optional<MyMatrix<T>> opt =
          MatrixIntegral_RepresentativeAction<T, Tgroup, Thelper>(
              eret, GRPperm, helper, eFace1, eFace2);
      if (!opt) {
#ifdef DEBUG_MATRIX_GROUP
        std::cerr << "Exit while loop with proof that no equivalence exists\n";
#endif
        return {};
      }
      MyMatrix<T> const &M = *opt;
#ifdef DEBUG_MATRIX_GROUP
      MyMatrix<Tmod> Mmod = ModuloReductionMatrix<T, Tmod>(M, TheMod);
      Treturn fret =
          MatrixIntegral_GeneratePermutationGroup<T, Tmod, Telt, Thelper>(
              {M}, {Mmod}, helper, O, TheMod);
      if (fret.ListPermGens.size() != 1) {
        std::cerr << "ListPermGens does not have the right length\n";
        throw TerminalException{1};
      }
      Telt ePerm = fret.ListPermGens[0];
      std::cerr << "ePerm=" << ePerm << "\n";
      if constexpr (has_determining_ext<Thelper>::value) {
        MyMatrix<T> M2 = RepresentPermutationAsMatrix(helper, ePerm);
        if (M != M2) {
          std::cerr << "The matrix is not the original one\n";
          throw TerminalException{1};
        }
      }
      Face eFace1_img = OnFace(eFace1, ePerm);
      if (eFace1_img != eFace2) {
        std::cerr << "eFace1 not mapped to eFace2\n";
        std::cerr << "|eFace1|=" << eFace1.size() << "\n";
        std::cerr << "nbRow=" << fret.nbRow << " siz=" << fret.siz << "\n";
        std::cerr << "eFace1_img=" << StringFace(eFace1_img) << "\n";
        throw TerminalException{1};
      }
#endif
      eElt = eElt * M;
      if (!NeedStabilizer) {
        std::optional<MyVector<Tmod>> test1 = IsEquiv(eElt);
        if (!test1) {
#ifdef DEBUG_MATRIX_GROUP
          std::cerr << "eElt and GRPwork are correct. Exiting\n";
#endif
          ResultTestModEquivalence<T> res{ListMatrRet, eElt};
          return res;
        }
      }
      ListMatrRet = MatrixIntegral_Stabilizer<T, Tgroup, Thelper>(
          eret, GRPperm, helper, eFace2);
      ListMatrRetMod =
          ModuloReductionStdVectorMatrix<T, Tmod>(ListMatrRet, TheMod);
    } else {
      MyVector<Tmod> const &V = *test2;
      std::vector<MyVector<Tmod>> O =
          OrbitComputation(ListMatrRetMod, V, TheAction);
#ifdef DEBUG_MATRIX_GROUP
      std::cerr << "|O|=" << O.size() << "\n";
#endif
      Treturn eret =
          MatrixIntegral_GeneratePermutationGroup<T, Tmod, Telt, Thelper>(
              ListMatrRet, ListMatrRetMod, helper, O, TheMod);
      Tgroup GRPperm(eret.ListPermGens, eret.siz);
      Face eFace2 = GetFace<T, Tmod>(eret.nbRow, O, TheSpace2Mod);
#ifdef MATRIX_GROUP_DIAGNOSTICS
      std::cerr << "ModEquivalence 2 TheMod=" << TheMod << " |O|=" << O.size()
                << " |GRPperm|=" << GRPperm.size()
                << " |eFace2|=" << eFace2.count() << "\n";
#endif
      ListMatrRet = MatrixIntegral_Stabilizer<T, Tgroup, Thelper>(
          eret, GRPperm, helper, eFace2);
      ListMatrRetMod =
          ModuloReductionStdVectorMatrix<T, Tmod>(ListMatrRet, TheMod);
    }
  }
}

template <typename T, typename Tgroup, typename Thelper>
std::optional<ResultTestModEquivalence<T>>
LinearSpace_ModEquivalence(std::vector<MyMatrix<T>> const &ListMatr,
                           Thelper const &helper, bool const &NeedStabilizer,
                           MyMatrix<T> const &TheSpace1,
                           MyMatrix<T> const &TheSpace2, T const &TheMod) {
  T max_size = (TheMod - 1) * (TheMod - 1) * TheSpace1.rows();
  if (max_size < T(std::numeric_limits<uint8_t>::max())) {
    return LinearSpace_ModEquivalence_Tmod<T, uint8_t, Tgroup, Thelper>(
        ListMatr, helper, NeedStabilizer, TheSpace1, TheSpace2, TheMod);
  }
  if (max_size < T(std::numeric_limits<uint16_t>::max())) {
    return LinearSpace_ModEquivalence_Tmod<T, uint16_t, Tgroup, Thelper>(
        ListMatr, helper, NeedStabilizer, TheSpace1, TheSpace2, TheMod);
  }
  if (max_size < T(std::numeric_limits<uint32_t>::max())) {
    return LinearSpace_ModEquivalence_Tmod<T, uint32_t, Tgroup, Thelper>(
        ListMatr, helper, NeedStabilizer, TheSpace1, TheSpace2, TheMod);
  }
  std::cerr << "Failed to find a matching arithmetic type. Quite unlikely "
               "objectively\n";
  throw TerminalException{1};
}

template <typename T, typename Tgroup, typename Thelper>
std::vector<MyMatrix<T>>
LinearSpace_Stabilizer_Kernel(std::vector<MyMatrix<T>> const &ListMatr,
                              Thelper const &helper,
                              MyMatrix<T> const &TheSpace) {
  int n = helper.n;
#ifdef DEBUG_MATRIX_GROUP
  //  std::cerr << "TheSpace=\n";
  //  WriteMatrixGAP(std::cerr, TheSpace);
  //  std::cerr << "\n";
  std::cerr << "det(TheSpace)=" << DeterminantMat(TheSpace) << "\n";
#endif
  CanSolIntMat<T> eCan = ComputeCanonicalFormFastReduction(TheSpace);
  auto IsStabilizing = [&](std::vector<MyMatrix<T>> const &LMat) -> bool {
    for (int i = 0; i < n; i++) {
      MyVector<T> eVect = GetMatrixRow(TheSpace, i);
      for (auto &eGen : LMat) {
        MyVector<T> eVectG = eGen.transpose() * eVect;
        bool test = CanTestSolutionIntMat(eCan, eVectG);
        if (!test) {
          return false;
        }
      }
    }
#ifdef DEBUG_MATRIX_GROUP
    std::cerr << "Leaving IsStabilzing: true\n";
#endif
    return true;
  };
  if (IsStabilizing(ListMatr))
    return ListMatr;
  T LFact = LinearSpace_GetDivisor(TheSpace);
  std::vector<T> eList = FactorsInt(LFact);
  int siz = eList.size();
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "LFact=" << LFact << " siz=" << siz << "\n";
#endif
  std::vector<MyMatrix<T>> ListMatrRet = ListMatr;
  for (int i = 1; i <= siz; i++) {
    T TheMod = 1;
    for (int j = 0; j < i; j++)
      TheMod *= eList[j];
    ListMatrRet = LinearSpace_ModStabilizer<T, Tgroup>(ListMatrRet, helper,
                                                       TheSpace, TheMod);
    if (IsStabilizing(ListMatrRet))
      return ListMatrRet;
  }
#ifdef SANITY_CHECK
  if (!IsStabilizing(ListMatrRet)) {
    std::cerr << "Error in LinearSpace_Stabilizer_Kernel\n";
    throw TerminalException{1};
  }
#endif
  return ListMatrRet;
}

template <typename T, typename Tgroup, typename Thelper>
std::vector<MyMatrix<T>>
LinearSpace_Stabilizer(std::vector<MyMatrix<T>> const &ListMatr,
                       Thelper const &helper, MyMatrix<T> const &TheSpace) {
  using Tint = typename underlying_ring<T>::ring_type;
  std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>> pair =
      LLLMatrixGroupReduction<T, Tint, Thelper>(helper, ListMatr);
  std::vector<MyMatrix<T>> const &ListMatrNew = pair.first;
  MyMatrix<Tint> const &Pmat = pair.second;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  MyMatrix<T> TheSpace_B = TheSpace * PmatInv_T;
  MyMatrix<T> TheSpace_C = LLLbasisReduction<T, Tint>(TheSpace_B).LattRed;
  Thelper helper_new = TransformHelper(helper, Pmat_T);
  std::vector<MyMatrix<T>> ListMatr_B =
      LinearSpace_Stabilizer_Kernel<T, Tgroup, Thelper>(ListMatrNew, helper_new,
                                                        TheSpace_C);
  std::vector<MyMatrix<T>> ListMatr_C;
  for (auto &eMatr_B : ListMatr_B) {
    MyMatrix<T> eMatr_C = PmatInv_T * eMatr_B * Pmat_T;
    ListMatr_C.push_back(eMatr_C);
  }
  return ListMatr_C;
}

template <typename T, typename Tgroup, typename Thelper>
std::optional<MyMatrix<T>> LinearSpace_Equivalence_Kernel(
    std::vector<MyMatrix<T>> const &ListMatr, Thelper const &helper,
    MyMatrix<T> const &InSpace1, MyMatrix<T> const &InSpace2) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in LinearSpace_Equivalence_Kernel");
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Beginning of LinearSpace_Equivalence_Kernel\n";
  std::cerr << "InSpace1=\n";
  WriteMatrix(std::cerr, InSpace1);
  std::cerr << "InSpace2=\n";
  WriteMatrix(std::cerr, InSpace2);
  std::cerr << "|ListMatr|=" << ListMatr.size() << "\n";
  for (auto &eMatr : ListMatr) {
    std::cerr << "eMatr=\n";
    WriteMatrix(std::cerr, eMatr);
  }
  std::cerr << "Det(InSpace1)=" << DeterminantMat(InSpace1)
            << " Det(InSpace2)=" << DeterminantMat(InSpace2) << "\n";
#endif
  FractionMatrix<T> eRec1 = RemoveFractionMatrixPlusCoeff(InSpace1);
  FractionMatrix<T> eRec2 = RemoveFractionMatrixPlusCoeff(InSpace2);
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "eRec1.TheMult=" << eRec1.TheMult
            << " eRec2.TheMult=" << eRec2.TheMult << "\n";
#endif
  if (eRec1.TheMult != eRec2.TheMult)
    return {};
  MyMatrix<T> const &TheSpace1 = eRec1.TheMat;
  MyMatrix<T> const &TheSpace2 = eRec2.TheMat;
  //
  int n = TheSpace1.rows();
  T LFact1 = LinearSpace_GetDivisor(TheSpace1);
  T LFact2 = LinearSpace_GetDivisor(TheSpace2);
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "LFact1 = " << LFact1 << "\n";
  std::cerr << "LFact2 = " << LFact2 << "\n";
#endif
  if (LFact1 != LFact2) {
    return {};
  }
  std::vector<T> eList = FactorsInt(LFact1);
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "|eList|=" << eList.size() << " eList =";
  for (auto &eVal : eList)
    std::cerr << " " << eVal;
  std::cerr << "\n";
#endif
  CanSolIntMat<T> eCan = ComputeCanonicalFormFastReduction(TheSpace2);
  auto IsEquivalence = [&](MyMatrix<T> const &eEquiv) -> bool {
    for (int i = 0; i < n; i++) {
      MyVector<T> eVect = GetMatrixRow(TheSpace1, i);
      MyVector<T> eVectG = eEquiv.transpose() * eVect;
      bool test = CanTestSolutionIntMat(eCan, eVectG);
      if (!test)
        return false;
    }
    return true;
  };
  std::vector<MyMatrix<T>> ListMatrWork = ListMatr;
  int siz = eList.size();
  std::vector<MyMatrix<T>> ListMatrRet = ListMatr;
  MyMatrix<T> eElt = IdentityMat<T>(n);
  for (int i = 1; i <= siz; i++) {
    if (IsEquivalence(eElt))
      return eElt;
    T TheMod = 1;
    for (int j = 0; j < i; j++)
      TheMod *= eList[j];
    MyMatrix<T> TheSpace1Img = TheSpace1 * eElt;
    bool NeedStabilizer = true;
    if (i == siz)
      NeedStabilizer = false;
    std::optional<ResultTestModEquivalence<T>> opt =
        LinearSpace_ModEquivalence<T, Tgroup, Thelper>(
            ListMatrWork, helper, NeedStabilizer, TheSpace1Img, TheSpace2,
            TheMod);
    if (!opt) {
#ifdef DEBUG_MATRIX_GROUP
      std::cerr << "LinearSpace_ModEquivalence failed so we exit here\n";
#endif
      return {};
    }
    eElt = eElt * (opt->second);
    if (NeedStabilizer)
      ListMatrWork = opt->first;
  }
#ifdef SANITY_CHECK
  if (!IsEquivalence(eElt)) {
    std::cerr << "Error in LinearSpace_Equivalence_Kernel\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Before returning from LinearSpace_Equivalence_Kernel, retuning "
               "eElt\n";
#endif
  return eElt;
}

template <typename T, typename Tgroup, typename Thelper>
std::optional<MyMatrix<T>>
LinearSpace_Equivalence(std::vector<MyMatrix<T>> const &ListMatr,
                        Thelper const &helper, MyMatrix<T> const &InSpace1,
                        MyMatrix<T> const &InSpace2) {
  //  return LinearSpace_Equivalence_Kernel<T,Tgroup,Thelper>(ListMatr, helper,
  //  InSpace1, InSpace2);
  using Tint = typename underlying_ring<T>::ring_type;
  std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>> pair =
      LLLMatrixGroupReduction<T, Tint, Thelper>(helper, ListMatr);
  std::vector<MyMatrix<T>> const &ListMatrNew = pair.first;
  MyMatrix<Tint> const &Pmat = pair.second;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  MyMatrix<T> InSpace1_B = InSpace1 * PmatInv_T;
  MyMatrix<T> InSpace2_B = InSpace2 * PmatInv_T;
  MyMatrix<T> InSpace1_C = LLLbasisReduction<T, Tint>(InSpace1_B).LattRed;
  MyMatrix<T> InSpace2_C = LLLbasisReduction<T, Tint>(InSpace2_B).LattRed;
  Thelper helper_new = TransformHelper(helper, Pmat_T);
  std::optional<MyMatrix<T>> opt =
      LinearSpace_Equivalence_Kernel<T, Tgroup, Thelper>(
          ListMatrNew, helper_new, InSpace1_C, InSpace2_C);
  if (!opt)
    return {};
  MyMatrix<T> RetMat = PmatInv_T * (*opt) * Pmat_T;
  return RetMat;
}

template <typename T, typename Tgroup>
std::vector<MyMatrix<T>> LinPolytopeIntegral_Automorphism_Subspaces(
    std::vector<MyMatrix<T>> const &ListMatr, MyMatrix<T> const &EXTfaithful) {
  static_assert(
      is_ring_field<T>::value,
      "Requires T to be a field in LinPolytopeIntegral_Automorphism_Subspaces");
  using Telt = typename Tgroup::Telt;
  MyMatrix<T> eBasis = GetZbasis(EXTfaithful);
  MyMatrix<T> InvBasis = Inverse(eBasis);
  MyMatrix<T> EXTbas = EXTfaithful * InvBasis;
  std::vector<MyMatrix<T>> ListMatrGens;
  for (auto &eGen : ListMatr) {
    MyMatrix<T> NewGen = eBasis * eGen * InvBasis;
#ifdef SANITY_CHECK
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrGens.push_back(NewGen);
  }
  FiniteMatrixGroupHelper<T, Telt> helper =
      ComputeFiniteMatrixGroupHelper<T, Telt>(EXTbas);
  MyMatrix<T> LattToStab = RemoveFractionMatrix(Inverse(eBasis));

  std::vector<MyMatrix<T>> LMat =
      LinearSpace_Stabilizer<T, Tgroup>(ListMatrGens, helper, LattToStab);
  std::vector<MyMatrix<T>> ListMatrGensB;
  for (auto &eGen : LMat) {
    MyMatrix<T> NewGen = InvBasis * eGen * eBasis;
#ifdef SANITY_CHECK
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrGensB.push_back(NewGen);
  }
  return ListMatrGensB;
}

template <typename T, typename Tgroup, typename Fcorrect>
std::optional<MyMatrix<T>> LinPolytopeIntegral_Isomorphism_Method4(
    MyMatrix<T> const &EXT1_T, MyMatrix<T> const &EXT2_T, Tgroup const &GRP1,
    typename Tgroup::Telt const &ePerm, Fcorrect f_correct) {
  using Telt = typename Tgroup::Telt;
  for (auto &fPerm : GRP1) {
    Telt eEquivCand = fPerm * ePerm;
    MyMatrix<T> eBigMat = FindTransformation(EXT1_T, EXT2_T, eEquivCand);
    if (f_correct(eBigMat))
      return eBigMat;
  }
  return {};
}

template <typename T, typename Tgroup, typename Fcorrect>
Tgroup LinPolytopeIntegral_Stabilizer_Method4(MyMatrix<T> const &EXT_T,
                                              Tgroup const &GRPisom,
                                              Fcorrect f_correct) {
  static_assert(
      is_ring_field<T>::value,
      "Requires T to be a field in LinPolytopeIntegral_Stabilizer_Method4");
  using Telt = typename Tgroup::Telt;
  int nbVert = EXT_T.rows();
  std::vector<Telt> generatorList;
  Tgroup GRPret(nbVert);
  auto fInsert = [&](Telt const &ePerm) -> void {
    bool test = GRPret.isin(ePerm);
    if (!test) {
      generatorList.push_back(ePerm);
      GRPret = Tgroup(generatorList, nbVert);
    }
  };
  for (auto &ePerm : GRPisom) {
    MyMatrix<T> eBigMat = FindTransformation(EXT_T, EXT_T, ePerm);
    if (f_correct(eBigMat))
      fInsert(ePerm);
  }
  return GRPret;
}

template <typename T, typename Tgroup>
Tgroup LinPolytopeIntegral_Stabilizer_Method8(MyMatrix<T> const &EXT_T,
                                              Tgroup const &GRPisom) {
  static_assert(
      is_ring_field<T>::value,
      "Requires T to be a field in LinPolytopeIntegral_Stabilizer_Method8");
  using Telt = typename Tgroup::Telt;
  int nbVert = EXT_T.rows();
  std::vector<MyMatrix<T>> ListMatrGen;
  for (auto &eGen : GRPisom.GeneratorsOfGroup()) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    ListMatrGen.push_back(eMat);
  }
  using Thelper = FiniteMatrixGroupHelper<T, Telt>;
  Thelper helper = ComputeFiniteMatrixGroupHelper<T, Telt>(EXT_T);
  std::vector<MyMatrix<T>> ListMatr =
      LinPolytopeIntegral_Automorphism_Subspaces<T, Tgroup>(ListMatrGen, EXT_T);
  std::vector<Telt> ListPermGens;
  for (auto &eMatr : ListMatr)
    ListPermGens.push_back(
        GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(helper, eMatr));
  return Tgroup(ListPermGens, nbVert);
}

template <typename T, typename Tgroup>
std::optional<MyMatrix<T>> LinPolytopeIntegral_Isomorphism_Subspaces(
    MyMatrix<T> const &EXT1_T, MyMatrix<T> const &EXT2_T,
    std::vector<MyMatrix<T>> const &ListMatrGens2,
    typename Tgroup::Telt const &eEquiv) {
  static_assert(
      is_ring_field<T>::value,
      "Requires T to be a field in LinPolytopeIntegral_Isomorphism_Subspaces");
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "Beginning of LinPolytopeIntegral_Isomorphism_Subspaces\n";
#endif
  using Telt = typename Tgroup::Telt;
  MyMatrix<T> eBasis1 = GetZbasis(EXT1_T);
  MyMatrix<T> eBasis2 = GetZbasis(EXT2_T);
  MyMatrix<T> InvBasis1 = Inverse(eBasis1);
  MyMatrix<T> InvBasis2 = Inverse(eBasis2);
  MyMatrix<T> EXTbas1 = EXT1_T * InvBasis1;
  MyMatrix<T> EXTbas2 = EXT2_T * InvBasis2;
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "EXT1_T=\n";
  WriteMatrix(std::cerr, EXT1_T);
  std::cerr << "EXT2_T=\n";
  WriteMatrix(std::cerr, EXT2_T);
  if (EXT1_T.rows() == EXT1_T.cols()) {
    std::cerr << "det(EXT1_T)=" << DeterminantMat(EXT1_T) << "\n";
  }
  if (EXT2_T.rows() == EXT2_T.cols()) {
    std::cerr << "det(EXT2_T)=" << DeterminantMat(EXT2_T) << "\n";
  }
  std::cerr << "eBasis1=\n";
  WriteMatrix(std::cerr, eBasis1);
  std::cerr << "eBasis2=\n";
  WriteMatrix(std::cerr, eBasis2);
  std::cerr << "EXTbas1=\n";
  WriteMatrix(std::cerr, EXTbas1);
  std::cerr << "EXTbas2=\n";
  WriteMatrix(std::cerr, EXTbas2);
#endif
#ifdef SANITY_CHECK
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
#endif
  //
  MyMatrix<T> TheMatEquiv = FindTransformation(EXTbas1, EXTbas2, eEquiv);
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "TheMatEquiv=\n";
  WriteMatrix(std::cerr, TheMatEquiv);
#endif
  std::vector<MyMatrix<T>> ListMatrGen;
  for (auto &eGen : ListMatrGens2) {
    MyMatrix<T> NewGen = eBasis2 * eGen * InvBasis2;
    ListMatrGen.push_back(NewGen);
  }
  FiniteMatrixGroupHelper<T, Telt> helper =
      ComputeFiniteMatrixGroupHelper<T, Telt>(EXTbas2);
  MyMatrix<T> eLatt1 = Inverse(eBasis1) * TheMatEquiv;
  MyMatrix<T> eLatt2 = Inverse(eBasis2);
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "eLatt1=\n";
  WriteMatrix(std::cerr, eLatt1);
  std::cerr << "eLatt2=\n";
  WriteMatrix(std::cerr, eLatt2);
#endif
  std::optional<MyMatrix<T>> opt =
      LinearSpace_Equivalence<T, Tgroup>(ListMatrGen, helper, eLatt1, eLatt2);
  if (!opt)
    return {};
  MyMatrix<T> const &eSpaceEquiv = *opt;
  MyMatrix<T> eMatFinal = InvBasis1 * TheMatEquiv * eSpaceEquiv * eBasis2;
#ifdef DEBUG_MATRIX_GROUP
  std::cerr << "We have eMatFinal=\n";
  WriteMatrix(std::cerr, eMatFinal);
#endif
#ifdef SANITY_CHECK
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
    typename Tgroup::Telt const &ePerm) {
  using Telt = typename Tgroup::Telt;
  std::vector<MyMatrix<T>> ListMatrGens;
  std::vector<Telt> LGen = GRP1.GeneratorsOfGroup();
  for (auto &eGen : LGen) {
    Telt ePermGen = (~ePerm) * eGen * ePerm;
    MyMatrix<T> eMatr = FindTransformation(EXT2_T, EXT2_T, ePermGen);
    ListMatrGens.push_back(eMatr);
  }
  return LinPolytopeIntegral_Isomorphism_Subspaces<T, Tgroup>(
      EXT1_T, EXT2_T, ListMatrGens, ePerm);
}

// clang-format off
#endif  // SRC_POLY_MATRIXGROUP_H_
// clang-format on

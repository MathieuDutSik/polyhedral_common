// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_MATRIXCANONICALFORM_H_
#define SRC_LATT_MATRIXCANONICALFORM_H_

// clang-format off
#include "InvariantVectorFamily.h"
#include "MAT_MatrixInt.h"
#include "MatrixGroup.h"
#include "Shvec_exact.h"
#include "PolytopeEquiStabInt.h"
#include <utility>
#include <string>
#include <vector>
// clang-format on

/*
  This is the set of algorithms for computing:
  * Arithmetic equivalence under GL_n(Z) of positive definite quadratic forms.
  * Stabilizer under GL_n(Z) action of positive definite quadratic forms.
  * Canonical form under GL_n(Z) action of positive definite form.

  The classic algorithm for iterm 1 and 2 is the Plesken-Souvignier algorithm:
  W. Plesken, B. Souvignier, Computing isometries of lattices, Journal of
     Symbolic Computation (1997) 24, 327-344.
  It works by using the standard basis (e_1, ...., e_n) with e_i the vector
  having
  (e_i)_j = | 1 is i=j
            | 0 if i<>j
  That standard basis is then looked for possible images in the list of short
  vectors. This is very efficient. The algorithm of this section are different.
  Instead we compute the set of pairwise scalar products and use partition
  backtrack for testing isomorphism and stabilizer:
  * This approach is not scalable for large lattice such as Leech lattice.
    Here PS algorithm works better.
  * But it may be more efficient for small dimensional lattices.

  For the canonical form computation, we have
  Mathieu Dutour Sikirić, Anna Haensch, John Voight, Wessel Van Woerden,
  A canonical form for positive definite matrices, Proceedings of the
  Fourteenth Algorithmic Number Theory Symposium (ANTS-XIV), edited by
  Steven Galbraith, Open Book Series 4, Mathematical Sciences Publishers,
  Berkeley, 2020, prepring at https://arxiv.org/abs/2004.14022

  The approach of Plesken-Souvignier does not appear to be feasible for
  computing canonical form.

  So, the approach of this file is to compute configuration of short
  vectors and then compute automorphism, isomorphism and canonical forms.

  Computing those configurations of short vectors is a difficult
  business. Two tricks are available to accelerate it:
  * Computing also for the dual. The inverse matrix A^(-1) might be
    easier to compute with.
  * We do not necessarily need a spanning configuration. Having a full
    rank one (but not neceessarily spanning) that can then be used for
    rational stabilizer and equivalence is good enough. Then we can
    use the finite index algorithm for concluding.

 */


#ifdef DEBUG
#define DEBUG_LATTICE_STAB_EQUI_CAN
#endif

#ifdef DISABLE_DEBUG_LATTICE_STAB_EQUI_CAN
#undef DEBUG_LATTICE_STAB_EQUI_CAN
#endif

#ifdef TIMINGS
#define TIMINGS_LATTICE_STAB_EQUI_CAN
#endif

template <typename T, typename Tint> struct Canonic_PosDef {
  MyMatrix<Tint> Basis;
  MyMatrix<Tint> SHV;
  MyMatrix<T> Mat;
};

//
// The canonic form
//

template <typename T, typename Tint>
MyMatrix<Tint> CanonicallyReorder_SHV(std::vector<MyMatrix<T>> const &ListMat,
                                      MyMatrix<Tint> const &SHV,
                                      std::ostream &os) {
  using Tgr = GraphListAdj;
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: Begining of ComputeCanonicalForm\n";
#endif
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  int nbRow = SHV.rows();
  int n = SHV.cols();
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: nbRow=" << nbRow << " n=" << n << "\n";
  if (!CheckCentralSymmetry(SHV)) {
    std::cerr
        << "LSEC: The set of vector does not respect the central symmetry "
           "condition\n";
    throw TerminalException{1};
  }
#endif
  //
  // Computing the scalar product matrix
  //
  using Tidx_value = int16_t;
  WeightMatrix<false, std::vector<T>, Tidx_value> WMat =
      T_TranslateToMatrix_ListMat_SHV<T, Tint, Tidx_value>(ListMat, SHV, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: WMat|=" << time << "\n";
#endif
  WMat.ReorderingSetWeight();
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ReorderingSetWeight|=" << time << "\n";
#endif
  //
  // Computing the canonicalization of the scalar product matrix
  //
  std::vector<int> CanonicOrd =
    GetCanonicalizationVector_Kernel<std::vector<T>, Tgr, int>(WMat, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CNA
  os << "|LSEC: GetCanonicalizationVector_Kernel|=" << time << "\n";
#endif
  //
  // Building the canonical basis
  //
  MyMatrix<Tint> SHVcan(n, nbRow);
  for (int iRowCan = 0; iRowCan < nbRow; iRowCan++) {
    int iRowNative = CanonicOrd[iRowCan];
    MyVector<Tint> eRow_Tint = GetMatrixRow(SHV, iRowNative);
    AssignMatrixCol(SHVcan, iRowCan, eRow_Tint);
  }
  return SHVcan;
}

template<typename Tint>
MyMatrix<Tint> get_canonicallization_matrix(MyMatrix<Tint> const& SHVcan, [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  MyMatrix<Tint> BasisCan_pre = ComputeRowHermiteNormalForm_first(SHVcan);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: get_canonicallization_matrix, BasisCan_pre|=" << time << "\n";
#endif
  MyMatrix<Tint> BasisCan = TransposedMat(Inverse(BasisCan_pre));
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: get_canonicallization_matrix, BasisCan|=" << time << "\n";
#endif
  return BasisCan;
}



template <typename T, typename Tint>
MyMatrix<Tint> ComputeCanonicalForm_inner(std::vector<MyMatrix<T>> const &ListMat,
                                          MyMatrix<Tint> const &SHV,
                                          std::ostream &os) {
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  MyMatrix<Tint> SHVcan = CanonicallyReorder_SHV<T,Tint>(ListMat, SHV, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: SHVcan|=" << time << "\n";
#endif
  MyMatrix<Tint> BasisCan = get_canonicallization_matrix(SHVcan, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: BasisCan|=" << time << "\n";
#endif
  return BasisCan;
}

template <typename T, typename Tint>
MyMatrix<Tint> ComputeCanonicalForm(MyMatrix<T> const &inpMat, std::ostream &os) {
  //
  // Computing the Z-basis on which the computation relies.
  //
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(inpMat, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ExtractInvariantVectorFamilyZbasis|=" << time << "\n";
#endif
  std::vector<MyMatrix<T>> ListMat{inpMat};
  return ComputeCanonicalForm_inner(ListMat, SHV, os);
}

template <typename T, typename Tint>
MyMatrix<Tint> ComputeCanonicalFormMultiple(std::vector<MyMatrix<T>> const &ListMat,
                                            std::ostream &os) {
  //
  // Computing the Z-basis on which the computation relies.
  //
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: Begining of ComputeCanonicalForm\n";
#endif
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  MyMatrix<T> inpMat = ListMat[0];
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(inpMat, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ExtractInvariantVectorFamilyZbasis|=" << time << "\n";
#endif
  return ComputeCanonicalForm_inner<T, Tint>(ListMat, SHV, os);
}

template <typename T, typename Tint>
MyMatrix<Tint> ComputeCanonicalFormSymplectic(MyMatrix<T> const &inpMat, std::ostream &os) {
  int n_tot = inpMat.rows();
  if (n_tot % 2 == 1) {
    std::cerr << "LSEC: The dimension is odd\n";
    throw TerminalException{1};
  }
  int n = n_tot / 2;
  MyMatrix<T> SympFormMat = ZeroMatrix<T>(2 * n, 2 * n);
  for (int i = 0; i < n; i++) {
    SympFormMat(i, n + i) = 1;
    SympFormMat(n + i, i) = -1;
  }
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(inpMat, os);
  std::vector<MyMatrix<T>> ListMat{inpMat, SympFormMat};
  MyMatrix<Tint> SHVvan = CanonicallyReorder_SHV<T,Tint>(ListMat, SHV, os);
  MyMatrix<Tint> BasisSymp = SYMPL_ComputeSymplecticBasis(SHVvan);
  return BasisSymp;
}

//
// Automorphism code
//

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>> ArithmeticAutomorphismGroupMultiple_inner(
    std::vector<MyMatrix<T>> const &ListMat, MyMatrix<Tint> const &SHV,
    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using TintGroup = typename Tgroup::Tint;
  using Thelper = FiniteMatrixGroupHelper<T, Telt, TintGroup>;
  
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  int n_row = SHV_T.rows();
  std::vector<T> Vdiag(n_row, 0);

  std::vector<std::vector<Tidx>> ListGen =
      GetListGenAutomorphism_ListMat_Vdiag<T, T, Tgroup>(SHV_T, ListMat, Vdiag,
                                                         os);
  bool all_matrices_integral = true;
  std::vector<MyMatrix<T>> ListGenRet_T;
  for (auto &eList : ListGen) {
    std::optional<MyMatrix<T>> opt =
        FindMatrixTransformationTest(SHV_T, SHV_T, eList);
    if (!opt) {
      std::cerr << "LSEC: Failed to find the matrix\n";
      throw TerminalException{1};
    }
    MyMatrix<T> const &M_T = *opt;
    if (all_matrices_integral) {
      if (!IsIntegralMatrix(M_T)) {
        all_matrices_integral = false;
      }
    }
    ListGenRet_T.push_back(M_T);
  }
  std::vector<MyMatrix<Tint>> ListGenRet;
  if (all_matrices_integral) {
    for (auto &M_T : ListGenRet_T) {
      MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(M_T);
      ListGenRet.push_back(M);
    }
    return ListGenRet;
  }
  //
  // Not integral, computing the integral stabilizer
  //
  Thelper helper = ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(SHV_T);
  RetMI_S<T, Tgroup> ret =
    LinPolytopeIntegral_Automorphism_Subspaces<T, Tgroup>(ListGenRet_T, SHV_T, os);
  for (auto &M_T : ret.LGen) {
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
    if (!IsIntegralMatrix(M_T)) {
      std::cerr << "LSEC: Bug: The matrix should be integral\n";
      throw TerminalException{1};
    }
#endif
    MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(M_T);
    ListGenRet.push_back(M);
  }
  return ListGenRet;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
ArithmeticAutomorphismGroup_inner(MyMatrix<T> const &inpMat,
                                  MyMatrix<Tint> const &SHV, std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat{inpMat};
  return ArithmeticAutomorphismGroupMultiple_inner<T, Tint, Tgroup>(ListMat,
                                                                    SHV, os);
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
ArithmeticAutomorphismGroupMultiple(std::vector<MyMatrix<T>> const &ListMat,
                                    std::ostream &os) {
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  MyMatrix<Tint> SHV =
      ExtractInvariantVectorFamilyFullRank<T, Tint>(ListMat[0], os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ExtractInvariantVectorFamilyFullRank|=" << time << "\n";
#endif
  return ArithmeticAutomorphismGroupMultiple_inner<T, Tint, Tgroup>(ListMat,
                                                                    SHV, os);
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
ArithmeticAutomorphismGroup(MyMatrix<T> const &inpMat, std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat{inpMat};
  return ArithmeticAutomorphismGroupMultiple<T, Tint, Tgroup>(ListMat, os);
}

//
// Equivalence code
//

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>> ArithmeticEquivalenceMultiple_inner(
    std::vector<MyMatrix<T>> const &ListMat1, MyMatrix<Tint> const &SHV1,
    std::vector<MyMatrix<T>> const &ListMat2, MyMatrix<Tint> const &SHV2,
    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  auto f_check_mat=[&](MyMatrix<T> const& mat) -> void {
    for (size_t i = 0; i < ListMat1.size(); i++) {
      MyMatrix<T> eMat1 = ListMat1[i];
      MyMatrix<T> eMat2 = ListMat2[i];
      MyMatrix<T> eProd = mat * eMat1 * mat.transpose();
      if (eProd != eMat2) {
        std::cerr << "LSEC: Inconsistency error in the code\n";
        throw TerminalException{1};
      }
    }
  };
#endif
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: |SHV1|=" << SHV1.rows() << " |SHV2|=" << SHV2.rows() << "\n";
#endif

  if (SHV1.rows() != SHV2.rows())
    return {};
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: After the test\n";
#endif

  MyMatrix<T> SHV1_T = UniversalMatrixConversion<T, Tint>(SHV1);
  MyMatrix<T> SHV2_T = UniversalMatrixConversion<T, Tint>(SHV2);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: SHV1_T / SHV2_T|=" << time << "\n";
#endif

  int n_rows = SHV1_T.rows();
  std::vector<T> Vdiag(n_rows, 0);
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: Before the TestEquivalence_ListMat_Vdiag\n";
#endif
  std::optional<std::vector<Tidx>> opt1 =
      TestEquivalence_ListMat_Vdiag<T, T, Tidx>(SHV2_T, ListMat2, Vdiag,
                                                SHV1_T, ListMat1, Vdiag, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: TestEquivalence_ListMat_Vdiag|=" << time << "\n";
#endif

  if (!opt1)
    return {};
  std::optional<MyMatrix<T>> opt2 =
      FindMatrixTransformationTest(SHV2_T, SHV1_T, *opt1);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: FindMatrixTransformationTest|=" << time << "\n";
#endif
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  if (!opt2) {
    std::cerr << "LSEC: We have a matrix bug\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> const &MA_T = *opt2;
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  f_check_mat(MA_T);
#endif
  if (IsIntegralMatrix(MA_T)) {
    // Early termination
    MyMatrix<Tint> MA = UniversalMatrixConversion<Tint, T>(MA_T);
    return MA;
  }
  std::vector<std::vector<Tidx>> ListGen =
      GetListGenAutomorphism_ListMat_Vdiag<T, T, Tgroup>(SHV2_T, ListMat2, Vdiag, os);
  std::vector<MyMatrix<T>> ListMatrGens2;
  for (auto & eGen: ListGen) {
    Telt elt(eGen);
    MyMatrix<T> eMatrGen2 = FindTransformation(SHV2_T, SHV2_T, elt);
    ListMatrGens2.emplace_back(std::move(eMatrGen2));
  }
  Telt eEquiv(*opt1);
  Telt eEquivInv = Inverse(eEquiv);
  std::optional<MyMatrix<T>> opt3 = LinPolytopeIntegral_Isomorphism_Subspaces<T,Tgroup>(SHV1_T, SHV2_T, ListMatrGens2, eEquivInv, os);
  if (!opt3) {
    return {};
  }
  MyMatrix<T> MB_T = Inverse(*opt3);
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  if (!IsIntegralMatrix(MB_T)) {
    std::cerr << "LSEC: Matrix should be integral\n";
    throw TerminalException{1};
  }
  f_check_mat(MB_T);
#endif
  MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(MB_T);
  return M;
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>> ArithmeticEquivalence_inner(
    MyMatrix<T> const &inpMat1, MyMatrix<Tint> const &SHV1,
    MyMatrix<T> const &inpMat2, MyMatrix<Tint> const &SHV2, std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat1{inpMat1};
  std::vector<MyMatrix<T>> ListMat2{inpMat2};
  return ArithmeticEquivalenceMultiple_inner<T,Tint,Tgroup>(ListMat1, SHV1, ListMat2, SHV2, os);
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
ArithmeticEquivalenceMultiple(std::vector<MyMatrix<T>> const &ListMat1,
                              std::vector<MyMatrix<T>> const &ListMat2,
                              std::ostream &os) {
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  MyMatrix<Tint> SHV1 =
      ExtractInvariantVectorFamilyFullRank<T, Tint>(ListMat1[0], os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ExtractInvariantVectorFamilyFullRank1|=" << time << "\n";
#endif
  MyMatrix<Tint> SHV2 =
      ExtractInvariantVectorFamilyFullRank<T, Tint>(ListMat2[0], os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ExtractInvariantVectorFamilyFullRank2|=" << time << "\n";
#endif
  return ArithmeticEquivalenceMultiple_inner<T,Tint,Tgroup>(ListMat1, SHV1, ListMat2, SHV2, os);
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>> ArithmeticEquivalence(MyMatrix<T> const &inpMat1,
                                                    MyMatrix<T> const &inpMat2,
                                                    std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat1{inpMat1};
  std::vector<MyMatrix<T>> ListMat2{inpMat2};
  return ArithmeticEquivalenceMultiple<T, Tint, Tgroup>(ListMat1, ListMat2, os);
}

// clang-format off
#endif  // SRC_LATT_MATRIXCANONICALFORM_H_
// clang-format on

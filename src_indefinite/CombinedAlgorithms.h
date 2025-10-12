// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_
#define SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_

// clang-format off
#include "ApproximateModels.h"
#include "lorentzian_perfect.h"
#include "lorentzian_linalg.h"
#include "EquiStabMemoization.h"
#include "IndefApproxCanonical.h"
#include "MatrixGroup.h"
#include <vector>
// clang-format on

// This is a reimplementation of the GAP code and implements the algorithm for
// indefinite forms.
//
// The result of the paper were published in
// Mathieu Dutour Sikirić, Klaus Hulek, Moduli of polarized Enriques surfaces --
// computational aspects, Journal of the London Mathematical Society (2023)
// preprint at https://arxiv.org/abs/2302.01679

#ifdef DEBUG
#define DEBUG_INDEFINITE_COMBINED_ALGORITHMS
#endif

#ifdef TIMINGS
#define TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
#endif

#ifdef METHOD_COMPARISON
#define METHOD_COMPARISON_INDEFINITE_COMBINED_ALGORITHMS_ISOTROPIC
#endif

#ifdef TRACK_INFO
#define TRACK_INFO_INDEFINITE_COMBINED_ALGORITHMS
#endif

static const int INDEFINITE_FORM_PLANE = 32;
static const int INDEFINITE_FORM_FLAG = 92;

static const int METHOD_GENERATION_RIGHT_COSETS = 49;
static const int METHOD_GENERATION_DOUBLE_COSETS = 81;


template <typename T, typename Tint>
void check_equivalence(MyMatrix<T> const& Q1, MyMatrix<T> const& Q2, MyMatrix<Tint> const& equiv, std::string const& context) {
  MyMatrix<T> equiv_T = UniversalMatrixConversion<T, Tint>(equiv);
  MyMatrix<T> prod = equiv_T * Q1 * equiv_T.transpose();
  if (prod != Q2) {
    std::cerr << "COMB: Q1 should map to Q2 via equiv. context=" << context << "\n";
    throw TerminalException{1};
  }
}

template <typename T, typename Tint> struct INDEF_FORM_GetVectorStructure {
public:
  T eNorm;
  MyMatrix<T> Qmat;
  MyVector<Tint> v;
  MyVector<T> v_T;
  MyMatrix<T> Pmat;
  MyMatrix<T> PmatInv;
  MyMatrix<Tint> NSP;
  MyMatrix<T> NSP_T;
  MyMatrix<T> GramMatRed;
  std::ostream &os;
  INDEF_FORM_GetVectorStructure(MyMatrix<T> const &_Qmat,
                                MyVector<Tint> const &_v, std::ostream &_os)
      : Qmat(_Qmat), v(_v), os(_os) {
    int n = Qmat.rows();
    eNorm = EvaluationQuadForm<T, Tint>(Qmat, v);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_GetVectorStructure n=" << n << " eNorm=" << eNorm
       << "\n";
#endif
    v_T = UniversalVectorConversion<T, Tint>(v);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_GetVectorStructure v_T=" << StringVectorGAP(v_T)
       << "\n";
#endif
    MyVector<T> eProd = Qmat * v_T;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_GetVectorStructure eProd=" << StringVectorGAP(eProd)
       << "\n";
#endif
    NSP_T = NullspaceIntVect(eProd);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_GetVectorStructure NSP_T=\n";
    WriteMatrix(os, NSP_T);
#endif
    NSP = UniversalMatrixConversion<Tint, T>(NSP_T);
    GramMatRed = NSP_T * Qmat * NSP_T.transpose();
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_GetVectorStructure GramMatRed=\n";
    WriteMatrix(os, GramMatRed);
#endif
    if (eNorm != 0) {
      Pmat = ZeroMatrix<T>(n, n);
      for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n; j++) {
          Pmat(i, j) = NSP_T(i, j);
        }
      }
      for (int i = 0; i < n; i++) {
        Pmat(n - 1, i) = v_T(i);
      }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      os << "COMB: INDEF_FORM_GetVectorStructure, We have Pmat\n";
#endif
      PmatInv = Inverse(Pmat);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      os << "COMB: INDEF_FORM_GetVectorStructure, We have PmatInv\n";
#endif
    }
  }
  MyMatrix<T>
  MapOrthogonalSublatticeEndomorphism(MyMatrix<Tint> const &eEndoRed) {
    MyMatrix<T> eEndoRed_T = UniversalMatrixConversion<T, Tint>(eEndoRed);
    if (eNorm != 0) {
      MyMatrix<T> TheBigMat = ExpandMatrix(eEndoRed_T);
      MyMatrix<T> RetMat = PmatInv * TheBigMat * Pmat;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
      MyVector<T> prodV = RetMat.transpose() * v_T;
      if (prodV != v_T) {
        std::cerr << "COMB: RetMat is not preserving the vector v\n";
        throw TerminalException{1};
      }
#endif
      return RetMat;
    } else {
      MyMatrix<T> Subspace1 = eEndoRed_T * NSP_T;
      MyMatrix<T> const &Subspace2 = NSP_T;
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
      MicrosecondTime time;
#endif
      std::optional<MyMatrix<T>> opt =
          LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1(
              Qmat, Subspace1, Qmat, Subspace2);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
      os << "|COMB: LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1|=" << time << "\n";
#endif
      MyMatrix<T> RetMat = unfold_opt(
          opt, "opt should be something because NSP.rows = RankMat(NSP)");
      MyVector<T> vImg = RetMat.transpose() * v_T;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
      if (vImg != v_T && vImg != -v_T) {
        std::cerr << "COMB: RetMat should map v to v or -v\n";
        throw TerminalException{1};
      }
#endif
      if (vImg == v_T) {
        return RetMat;
      } else {
        return -RetMat;
      }
    }
  }
  std::vector<MyMatrix<T>>
  MapOrthogonalSublatticeGroup(std::vector<MyMatrix<Tint>> const &GRPmatr) {
    std::vector<MyMatrix<T>> NewListGen;
    for (auto &eGen : GRPmatr) {
      NewListGen.push_back(MapOrthogonalSublatticeEndomorphism(eGen));
    }
    return NewListGen;
  }
};

template <typename T, typename Tint> struct INDEF_FORM_Rec_IsotropicKplane {
public:
  MyMatrix<T> Qmat;
  MyMatrix<Tint> Plane;
  MyMatrix<T> Plane_T;
  MyMatrix<Tint> PlaneExpr;
  MyMatrix<T> PlaneExpr_T;
  int dimSpace;
  int dim;
  MyMatrix<Tint> NSP;
  MyMatrix<T> NSP_T;
  MyMatrix<T> GramMatRed;
  MyMatrix<T> QmatRed;
  MyMatrix<T> FullBasis_T;
  MyMatrix<Tint> FullBasis;
  MyMatrix<Tint> FullBasisInv;
  MyMatrix<T> TheCompl_T;
  MyMatrix<Tint> TheCompl;
  int the_dim;
  int dimCompl;
  std::ostream &os;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
  void check_generator(MyMatrix<Tint> const &eEndoRed,
                       MyMatrix<T> const &RetMat) const {
    MyMatrix<T> eEndoRed_T = UniversalMatrixConversion<T, Tint>(eEndoRed);
    MyMatrix<T> PlaneImg = Plane_T * RetMat;
    MyMatrix<T> TransRed(dim, dim);
    for (int u = 0; u < dim; u++) {
      MyVector<T> eV = GetMatrixRow(PlaneImg, u);
      std::optional<MyVector<T>> opt = SolutionMat(Plane_T, eV);
      MyVector<T> fV = unfold_opt(opt, "Get the vector fV");
      AssignMatrixRow(TransRed, u, fV);
    }
    T eDet = DeterminantMat(TransRed);
    if (T_abs(eDet) != 1) {
      std::cerr << "COMB: TransRed should have absolute determinant 1\n";
      throw TerminalException{1};
    }
    if (!TestEqualitySpaces(PlaneImg, Plane_T)) {
      std::cerr << "COMB: Plane should be invariant (isotropic case)\n";
      throw TerminalException{1};
    }
    int dimNSP = NSP_T.rows();
    MyMatrix<T> RetMat_red(dimNSP, dimNSP);
    for (int u = 0; u < dimNSP; u++) {
      MyVector<T> eV = GetMatrixRow(NSP_T, u);
      MyVector<T> fV = RetMat.transpose() * eV;
      std::optional<MyVector<T>> opt = SolutionMat(NSP_T, fV);
      MyVector<T> gV = unfold_opt(opt, "getting gV");
      AssignMatrixRow(RetMat_red, u, gV);
    }
    if (RetMat_red != eEndoRed_T) {
      std::cerr << "COMB: RetMat_red restricted to the nullspace is not the original eEndoRed\n";
      throw TerminalException{1};
    }
  }
#endif
  INDEF_FORM_Rec_IsotropicKplane(MyMatrix<T> const &_Qmat,
                                 MyMatrix<Tint> const &_Plane,
                                 std::ostream &_os)
      : Qmat(_Qmat), Plane(_Plane), dimSpace(Qmat.rows()), dim(Plane.rows()),
        os(_os) {
    Plane_T = UniversalMatrixConversion<T, Tint>(Plane);
    MyMatrix<T> eProd = Plane_T * Qmat;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    MyMatrix<T> WitnessIsotropy = Plane_T * Qmat * Plane_T.transpose();
    if (!IsZeroMatrix(WitnessIsotropy)) {
      std::cerr
          << "The matrix Plane does not define a totally isotropic space\n";
      throw TerminalException{1};
    }
#endif
    MyMatrix<T> PreNSP_T = SublatticeBasisReduction(NullspaceIntTrMat(eProd), os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_Rec_IsotropicKplane, eProd=\n";
    WriteMatrix(os, eProd);
    os << "COMB: INDEF_FORM_Rec_IsotropicKplane, PreNSP_T=\n";
    WriteMatrix(os, PreNSP_T);
#endif
    NSP_T = SublatticeBasisReduction(PreNSP_T, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_Rec_IsotropicKplane, NSP_T=\n";
    WriteMatrix(os, NSP_T);
#endif
    NSP = UniversalMatrixConversion<Tint, T>(NSP_T);
    GramMatRed = NSP_T * Qmat * NSP_T.transpose();
    the_dim = NSP_T.rows();
    PlaneExpr_T = MyMatrix<T>(dim, the_dim);
    for (int u = 0; u < dim; u++) {
      MyVector<T> eV = GetMatrixRow(Plane_T, u);
      std::optional<MyVector<T>> opt = SolutionIntMat(NSP_T, eV);
      MyVector<T> fV = unfold_opt(opt, "getting fV");
      AssignMatrixRow(PlaneExpr_T, u, fV);
    }
    PlaneExpr = UniversalMatrixConversion<Tint, T>(PlaneExpr_T);
    TheCompl_T = SubspaceCompletionInt(PlaneExpr_T, the_dim);
    TheCompl = UniversalMatrixConversion<Tint, T>(TheCompl_T);
    dimCompl = TheCompl_T.rows();
    FullBasis_T = Concatenate(TheCompl_T, PlaneExpr_T);
    FullBasis = UniversalMatrixConversion<Tint, T>(FullBasis_T);
    FullBasisInv = Inverse(FullBasis);
    QmatRed = TheCompl_T * GramMatRed * TheCompl_T.transpose();
  }
  // We map automorphism group of a sublattice to the full group.
  // The difficulty is that the rational kernel is non-trivial
  // if dim(Plane) > 1.
  // The extension maps an integral group of automorphism to a rational
  // group of automorphism of a higher dimensional space.
  //
  // This poses some painful problems:
  // * The rational stuff can be defined. But it will be infinite
  //   dimensional.
  // * From an element, we can compute the lift. We can minimize
  //   the denominators. The solution is never unique. But can we find
  //   something better.
  // * Given a lattice, can we test if a lift exists or not? Yes,
  //   definitely.
  // * The idea is that we can find some subspace being stabilized
  //   which makes it integral. This is the big ask. We have no reason
  // * When having the stabilized subspace, we can compute the
  //   relevant kernel.
  // * We can also compute the effective lifting method that stabilizes
  //   the subspace in question.

  //
  // Description of the function being used.
  //

  // Lift mapping. The sublattice in argument is supposed to help doing that.
  MyMatrix<T> LiftToFullAutomorphism(MyMatrix<Tint> const& eGenRed, MyMatrix<T> const& HelpingSublattice) const {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS_DISABLE
    os << "COMB: LiftToFullAutomorphism, step 1\n";
#endif
    MyMatrix<T> HelpingSublatticeInv = Inverse(HelpingSublattice);
    MyMatrix<T> Subspace1 = NSP_T * HelpingSublatticeInv;
    MyMatrix<T> eGenRed_T = UniversalMatrixConversion<T, Tint>(eGenRed);
    MyMatrix<T> Subspace2 = eGenRed_T * Subspace1;
    MyMatrix<T> QmatRed = HelpingSublattice * Qmat * HelpingSublattice.transpose();
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS_DISABLE
    os << "COMB: LiftToFullAutomorphism, step 2\n";
#endif
    LORENTZ_ExtendOrthogonalIsotropicIsomorphism<T> TheRec(QmatRed, Subspace1, QmatRed, Subspace2, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS_DISABLE
    os << "COMB: LiftToFullAutomorphism, step 3\n";
#endif
    MyMatrix<T> eGen1 = TheRec.get_one_transformation();
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS_DISABLE
    os << "COMB: LiftToFullAutomorphism, step 4\n";
#endif
    MyMatrix<T> eGen2 = HelpingSublatticeInv * eGen1 * HelpingSublattice;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    //    os << "COMB: HelpingSublattice=\n";
    //    WriteMatrix(os, HelpingSublattice);
    //    os << "COMB: |HelpingSublattice|=" << DeterminantMat(HelpingSublattice) << "\n";
    check_generator(eGenRed, eGen2);
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS_DISABLE
    os << "COMB: LiftToFullAutomorphism, step 5\n";
#endif
    return eGen2;
  }

  // The following function attempts to compute an invariant sublattice
  // For the computed lift of the matrix group. There is no theoretical
  // guarantee that the lift exists, but so far it has always happened.
  // Of course what we want is the smallest sublattice
  MyMatrix<T>
  ComputeInvariantSublattice_method1(std::vector<MyMatrix<Tint>> const &GRPmatr) const {
    size_t max_possval = std::numeric_limits<size_t>::max();
    int n = Qmat.rows();
    size_t n_gen = GRPmatr.size();
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: ComputeInvariantSublattice_method1, n_gen=" << n_gen << "\n";
#endif
    Face AlreadySolved(n_gen);
    MyMatrix<T> Sublattice = IdentityMat<T>(n);
    MyMatrix<T> SublatticeInv = IdentityMat<T>(n);
    while(true) {
      std::vector<size_t> list_pos;
      std::vector<size_t> list_idx;
      std::vector<MyMatrix<T>> list_lattice;
      std::vector<T> list_det;
      T max_det(0);
      size_t found_pos = max_possval;
      //
      std::vector<MyMatrix<T>> already_done;
      for (size_t i_gen=0; i_gen<n_gen; i_gen++) {
        if (AlreadySolved[i_gen] == 1) {
          MyMatrix<Tint> eGen1 = GRPmatr[i_gen];
          MyMatrix<T> eGen2 = LiftToFullAutomorphism(eGen1, Sublattice);
          MyMatrix<T> eGen3 = SublatticeInv * eGen2 * Sublattice;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
          if (!IsIntegralMatrix(eGen3)) {
            std::cerr << "COMB: eGen3 should be integral\n";
            throw TerminalException{1};
          }
#endif
          already_done.push_back(eGen3);
        }
      }
      size_t pos = 0;
      for (size_t i_gen=0; i_gen<n_gen; i_gen++) {
        if (AlreadySolved[i_gen] == 0) {
          MyMatrix<Tint> eGen1 = GRPmatr[i_gen];
          MyMatrix<T> eGen2 = LiftToFullAutomorphism(eGen1, Sublattice);
          MyMatrix<T> eGen3 = SublatticeInv * eGen2 * Sublattice;
          std::vector<MyMatrix<T>> LGen = already_done;
          LGen.push_back(eGen3);
          MyMatrix<T> SubSublatt = MatrixIntegral_GetInvariantSpace(n, LGen, os);
          T det = T_abs(DeterminantMat(SubSublatt));
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
          os << "COMB: ComputeInvariantSublattice_method1, i_gen=" << i_gen << " |LGen|=" << LGen.size() << " det=" << det << "\n";
#endif
          list_pos.push_back(pos);
          list_idx.push_back(i_gen);
          list_lattice.push_back(SubSublatt);
          list_det.push_back(det);
          pos += 1;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
          if (det > 1) {
            std::cerr << "COMB: The determinant should be lower or equal to 1\n";
            throw TerminalException{1};
          }
#endif
          if (det < 1) {
            // We want to reduce the determinant by as little as possible
            if (det > max_det) {
              max_det = det;
              found_pos = pos;
            }
          }
        }
      }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      os << "COMB: ComputeInvariantSublattice_method1, max_det=" << max_det << " found_pos=" << found_pos << "\n";
#endif
      if (found_pos == max_possval) {
        // No new lattice found. This means that the lattice is actually preserving everything.
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
        RecSolutionIntMat<T> eCan(Sublattice);
        for (auto & eGenRed: GRPmatr) {
          MyMatrix<T> eGen = LiftToFullAutomorphism(eGenRed, Sublattice);
          MyMatrix<T> Sublattice_img = Sublattice * eGen;
          if (!eCan.is_containing_m(Sublattice_img)) {
            std::cerr << "COMB: The sublattice should be preserved\n";
            throw TerminalException{1};
          }
        }
#endif
        return Sublattice;
      } else {
        Sublattice = list_lattice[found_pos] * Sublattice;
        SublatticeInv = Inverse(Sublattice);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        os << "COMB: ComputeInvariantSublattice_method1, Det(Sublattice)=" << DeterminantMat(Sublattice) << "\n";
#endif
        Face NewAlreadySolved(n_gen);
        for (size_t i_gen=0; i_gen<n_gen; i_gen++) {
          MyMatrix<Tint> eGen1 = GRPmatr[i_gen];
          MyMatrix<T> eGen2 = LiftToFullAutomorphism(eGen1, Sublattice);
          MyMatrix<T> eGen3 = SublatticeInv * eGen2 * Sublattice;
          if (IsIntegralMatrix(eGen3)) {
            NewAlreadySolved[i_gen] = 1;
          } else {
            if (AlreadySolved[i_gen] == 1) {
              std::cerr << "COMB: We have !NewAlreadySolved but we have AlreaadySolve so that is a bug\n";
              throw TerminalException{1};
            }
          }
        }
        AlreadySolved = NewAlreadySolved;
      }
    }
  }

  // The direct method since method1 appears to require more work.
  MyMatrix<T>
  ComputeInvariantSublattice_method2(std::vector<MyMatrix<Tint>> const &GRPmatr) const {
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    MicrosecondTime time;
#endif
    int n = Qmat.rows();
    MyMatrix<T> Sublattice = IdentityMat<T>(n);
    std::vector<MyMatrix<T>> LGen;
    for (auto & eGen1 : GRPmatr) {
      MyMatrix<T> eGen2 = LiftToFullAutomorphism(eGen1, Sublattice);
      LGen.push_back(eGen2);
    }
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: ComputeInvariantSublattice_method2 / LGen|=" << time << "\n";
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: ComputeInvariantSublattice_method2, we have |LGen|=" << LGen.size() << "\n";
#endif
    MyMatrix<T> SublatticeRet = MatrixIntegral_GetInvariantSpace(n, LGen, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: ComputeInvariantSublattice_method2, Det(SublatticeRet)=" << DeterminantMat(SublatticeRet) << "\n";
#endif
    return SublatticeRet;
  }

  MyMatrix<T>
  ComputeInvariantSublattice(std::vector<MyMatrix<Tint>> const &GRPmatr) const {
    return ComputeInvariantSublattice_method2(GRPmatr);
  }

  // Computes integral kernel relevant to the subspace
  // Nice to know would be to understand the structure of the group.
  std::vector<MyMatrix<T>> ComputeRelevantKernel(MyMatrix<T> const& Sublattice) const {
    return GetOrthogonalTotallyIsotropicKernelSubspace<T,Tint>(Qmat, NSP_T, Sublattice, os);
  }

  // We use the sublattice to map the group GRPmatr to higher dimension.
  std::vector<MyMatrix<T>>
  MapOrthogonalSublatticeGroupUsingSublattice(std::vector<MyMatrix<Tint>> const &GRPmatr, MyMatrix<T> const& Sublattice) const {
    std::vector<MyMatrix<T>> ListGens;
    //    int n = Qmat.rows();
    //    MyMatrix<T> IdMat = IdentityMat
    for (auto & eGenRed: GRPmatr) {
      MyMatrix<T> eGen = LiftToFullAutomorphism(eGenRed, Sublattice);
      ListGens.push_back(eGen);
    }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: MapOrthogonalSublatticeGroupUsingSublattice, after lifts\n";
#endif
    for (auto & eGen: ComputeRelevantKernel(Sublattice)) {
      ListGens.push_back(eGen);
    }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: MapOrthogonalSublatticeGroupUsingSublattice, after kernel gens\n";
#endif
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    RecSolutionIntMat<T> eCan(Sublattice);
    size_t i_gen = 0;
    for (auto & eGen : ListGens) {
      // We expect the lifted automorphism to preserve the lattice.
      // This is because first it should exist and second the reduction of the coefficients
      // should make the algorithm work.
      MyMatrix<T> Sublattice_img = Sublattice * eGen;
      if (!eCan.is_containing_m(Sublattice_img)) {
        std::cerr << "COMB: i_gen=" << i_gen << " |GRPmatr|=" << GRPmatr.size() << " |ListGens|=" << ListGens.size() << "\n";
        std::cerr << "COMB: The sublattice should be preserved\n";
        throw TerminalException{1};
      }
      i_gen += 1;
    }
#endif
    return ListGens;
  }


  std::vector<MyMatrix<T>>
  MapOrthogonalSublatticeGroup(std::vector<MyMatrix<Tint>> const &GRPmatr) const {
    MyMatrix<T> Sublattice = ComputeInvariantSublattice(GRPmatr);
    return MapOrthogonalSublatticeGroupUsingSublattice(GRPmatr, Sublattice);
  }




};

/*
  As obtained by GeneratorsOfGroup(GeneralLinearGroup(6,Integers))
  in GAP.
 */
template <typename T>
std::vector<MyMatrix<T>> GeneralLinearGroup(int const &n) {
  std::vector<MyMatrix<T>> ListGens;
  if (n > 1) {
    MyMatrix<T> mat1 = ZeroMatrix<T>(n, n);
    for (int i = 0; i < n; i++) {
      int iNext = 0;
      if (i < n - 1) {
        iNext = i + 1;
      }
      mat1(i, iNext) = 1;
    }
    ListGens.push_back(mat1);
  }
  //
  if (n > 2) {
    MyMatrix<T> mat2 = ZeroMatrix<T>(n, n);
    mat2(1, 0) = 1;
    mat2(0, 1) = 1;
    for (int i = 2; i < n; i++) {
      mat2(i, i) = 1;
    }
    ListGens.push_back(mat2);
  }
  //
  MyMatrix<T> mat3 = IdentityMat<T>(n);
  mat3(0, 0) = -1;
  ListGens.push_back(mat3);
  //
  if (n > 1) {
    MyMatrix<T> mat4 = IdentityMat<T>(n);
    mat4(0, 1) = 1;
    ListGens.push_back(mat4);
  }
  //
  return ListGens;
}

/*
  The successive dimensions of the flag.
  For the full plane of dimension k, dims={k}
  For the full sequence of flags, we have dims={1,1,.... , 1}
 */
struct SeqDims {
  std::vector<size_t> dims;
};

SeqDims seq_dims_reduced(SeqDims const& sd, int const& k) {
  size_t k_s = static_cast<size_t>(k);
  std::vector<size_t> dims;
  size_t sum_dim = 0;
  for (auto & edim: sd.dims) {
    if (sum_dim + edim > k_s) {
      size_t n_dim = k_s - sum_dim;
      dims.push_back(n_dim);
      return SeqDims{dims};
    }
    if (sum_dim + edim == k_s) {
      dims.push_back(edim);
      return SeqDims{dims};
    }
    dims.push_back(edim);
    sum_dim += edim;
  }
  std::cerr << "COMB: We should not reach that stage\n";
  throw TerminalException{1};
}

SeqDims seq_dims_append_one(SeqDims const& sd) {
  std::vector<size_t> dims = sd.dims;
  dims.push_back(1);
  SeqDims sd_ext{dims};
  return sd_ext;
}

void write_seq_dims(SeqDims const& sd, std::string const& name, std::ostream& os) {
  os << "COMB: SeqDims " << name << "=[ ";
  bool IsFirst = true;
  for (auto & edim: sd.dims) {
    if (!IsFirst) {
      os << " | ";
    }
    IsFirst = false;
    os << edim;
  }
  os << " ]\n";
}

SeqDims seq_dims_plane(int const& k) {
  std::vector<size_t> dims{static_cast<size_t>(k)};
  return SeqDims{dims};
}

SeqDims seq_dims_flag(int const& k) {
  std::vector<size_t> dims;
  for (int i=0; i<k; i++) {
    dims.push_back(1);
  }
  return SeqDims{dims};
}

int k_dim(SeqDims const& sd) {
  size_t sum = 0;
  for (auto & val : sd.dims) {
    sum += val;
  }
  return static_cast<int>(sum);
}

template<typename Tint>
std::vector<MyMatrix<Tint>> f_get_list_spaces(MyMatrix<Tint> const &ListVect, SeqDims const& sd, [[maybe_unused]] std::ostream& os) {
  size_t n_case = sd.dims.size();
  int dim = ListVect.cols();
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS_DISABLE
  os << "COMB: f_get_list_spaces, n_case=" << n_case << " dim=" << dim << " n_rows=" << ListVect.rows() << "\n";
  os << "COMB: f_get_list_spaces, ListVect=\n";
  WriteMatrix(os, ListVect);
#endif
  std::vector<MyMatrix<Tint>> ListSpaces;
  for (size_t i_case=0; i_case<n_case; i_case++) {
    int sum_dim = 0;
    for (size_t u=0; u<=i_case; u++) {
      sum_dim += sd.dims[u];
    }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS_DISABLE
    os << "COMB: f_get_list_spaces, i_case=" << i_case << " sum_dim=" << sum_dim << "\n";
#endif
    MyMatrix<Tint> eSpace(sum_dim, dim);
    for (int u = 0; u < sum_dim; u++) {
      for (int iCol = 0; iCol < dim; iCol++) {
        eSpace(u, iCol) = ListVect(u, iCol);
      }
    }
    ListSpaces.push_back(eSpace);
  }
  return ListSpaces;
}

template <typename T>
std::vector<MyMatrix<T>> GetAutomorphismOfFlag(SeqDims const& sd) {
  size_t n_case = sd.dims.size();
  std::vector<int> starts;
  int pos = 0;
  for (auto& edim : sd.dims) {
    starts.push_back(pos);
    pos += edim;
  }
  int n = pos;
  std::vector<MyMatrix<T>> LGen;
  for (size_t i_case=0; i_case<n_case; i_case++) {
    int start = starts[i_case];
    int dim = sd.dims[i_case];
    for (auto & eGen: GeneralLinearGroup<T>(dim)) {
      MyMatrix<T> TheMat = IdentityMat<T>(n);
      for (int i=0; i<dim; i++) {
        for (int j=0; j<dim; j++) {
          TheMat(i+start, j+start) = eGen(i,j);
        }
      }
      LGen.push_back(TheMat);
    }
  }
  for (size_t i_case = 0; i_case < n_case; i_case++) {
    for (size_t j_case = 0; j_case < i_case; j_case++) {
      int i = starts[i_case];
      int j = starts[j_case];
      MyMatrix<T> TheMat = IdentityMat<T>(n);
      TheMat(i, j) = 1;
      LGen.push_back(TheMat);
    }
  }
  return LGen;
}

/*
  Isometry group defined on a p dimensional space for a quadratic form Qp.
  We extend the quadratic form to dimension n with
  Qn = | Qp 0 |
       | 0  0 |
  If the original matrices satisfy Pp Qp Pp^T = Qp
  the the extended matrices must satisfy Pn Qn Pn^T = Qn
  So, in block formulation
  Pn = | A B |
       | C D |
  and so
  Pn Qn Pn^T = | A B |     | Qp 0 |     | A^T C^T |
               | C D |  x  | 0  0 |  x  | B^T D^T |
             = | A Qp 0 |     | A^T C^T |
               | C Qp 0 |  x  | B^T D^T |
             = | A Qp A^T  A Qp C^T |
               | C Qp A^T  C Qp C^T |
             = | Qp 0 |
               | 0  0 |
  And so we get A Qp A^T = Qp , C Qp A^T = 0 , C Qp C^T = 0
  The equation A Qp A^T forces A to be an isometry of Qp.
  The equation C Qp A^T forces C to be 0 from which the rest follows.
*/
template <typename T>
std::vector<MyMatrix<T>>
ExtendIsometryGroup(std::vector<MyMatrix<T>> const &GRPmatr, int const &p,
                    int const &n) {
  std::vector<MyMatrix<T>> ListGens;
  for (auto &eGen : GRPmatr) {
    MyMatrix<T> NewMat = IdentityMat<T>(n);
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        NewMat(i, j) = eGen(i, j);
      }
    }
    ListGens.push_back(NewMat);
  }
  if (n > p) {
    for (auto &eGen : GeneralLinearGroup<T>(n - p)) {
      MyMatrix<T> NewMat = IdentityMat<T>(n);
      for (int i = 0; i < n - p; i++) {
        for (int j = 0; j < n - p; j++) {
          NewMat(i + p, j + p) = eGen(i, j);
        }
      }
      ListGens.push_back(NewMat);
    }
    for (int i = 0; i < p; i++) {
      MyMatrix<T> NewMat = IdentityMat<T>(n);
      NewMat(i, p) = 1;
      ListGens.push_back(NewMat);
    }
  }
  return ListGens;
}

/*
  p is the dimension of the group GRPmatr.
  n is the full dimension of the space,
 */
template <typename T>
std::vector<MyMatrix<T>>
ExtendIsometryGroup_Triangular(std::vector<MyMatrix<T>> const &GRPmatr,
                               int const &p, int const &n, SeqDims const& sd) {
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
  if (p + k_dim(sd) != n) {
    std::cerr << "COMB: ExtendIsometryGroup_Triangular, the dimension of sd\n";
    std::cerr << "is not adequate\n";
    throw TerminalException{1};
  }
#endif
  std::vector<MyMatrix<T>> ListGens;
  for (auto &eGen : GRPmatr) {
    MyMatrix<T> NewMat = IdentityMat<T>(n);
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        NewMat(i, j) = eGen(i, j);
      }
    }
    ListGens.push_back(NewMat);
  }
  std::vector<MyMatrix<T>> SubNPgroup = GetAutomorphismOfFlag<T>(sd);
  for (auto &eGen : SubNPgroup) {
    MyMatrix<T> NewMat = IdentityMat<T>(n);
    for (int i = 0; i < n - p; i++) {
      for (int j = 0; j < n - p; j++) {
        NewMat(i + p, j + p) = eGen(i, j);
      }
    }
    ListGens.push_back(NewMat);
  }
  for (int i = 0; i < p; i++) {
    int idx = p;
    for (auto& edim : sd.dims) {
      MyMatrix<T> NewMat = IdentityMat<T>(n);
      NewMat(i, idx) = 1;
      ListGens.push_back(NewMat);
      idx += edim;
    }
  }
  return ListGens;
}

// Extend the isometry group which should preserve
// the eRec.QmatRed
template <typename T, typename Tint>
std::vector<MyMatrix<Tint>> ExtendIsometryGroup_IsotropicOrth(std::vector<MyMatrix<Tint>> const& GRPred,
                                                              INDEF_FORM_Rec_IsotropicKplane<T, Tint> const &eRec,
                                                              SeqDims const& sd,
                                                              [[maybe_unused]] std::ostream& os) {
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
  for (auto & eGen : GRPred) {
    check_equivalence(eRec.QmatRed, eRec.QmatRed, eGen, "eGen / eRec.QmatRed");
  }
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
  os << "COMB: EIG_IO, We have |GRPred|=" << GRPred.size() << " |QmatRed|=" << eRec.QmatRed.rows() << "\n";
  os << "COMB: EIG_IO, We have dimCompl=" << eRec.dimCompl << " the_dim=" << eRec.the_dim << " k=" << eRec.PlaneExpr.rows() << "\n";
  //  os << "COMB: EIG_IO, we have GRPred=\n";
  //  WriteListMatrix(os, GRPred);
#endif
  std::vector<MyMatrix<Tint>> GRPfull =
    ExtendIsometryGroup_Triangular(GRPred, eRec.dimCompl, eRec.the_dim, sd);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
  os << "COMB: EIG_IO, We have |GRPfull|=" << GRPfull.size() << " eRec.the_dim=" << eRec.the_dim << "\n";
  //  os << "COMB: EIG_IO, we have GRPfull=\n";
  //  WriteListMatrix(os, GRPfull);
#endif
  std::vector<MyMatrix<Tint>> ListGenTot;
  ListGenTot.push_back(-IdentityMat<Tint>(eRec.the_dim));
  for (auto &eGen : GRPfull) {
    MyMatrix<Tint> eGenB = eRec.FullBasisInv * eGen * eRec.FullBasis;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    check_equivalence(eRec.GramMatRed, eRec.GramMatRed, eGenB, "eGenB / eRec.GramMatRed");
    std::vector<MyMatrix<Tint>> ListSpaces =
      f_get_list_spaces(eRec.PlaneExpr, sd, os);
    for (auto &eSpace : ListSpaces) {
      MyMatrix<Tint> eSpaceImg = eSpace * eGenB;
      if (!TestEqualitySpaces(eSpaceImg, eSpace)) {
        std::cerr << "COMB: The space eSpace is not correctly preserved\n";
        throw TerminalException{1};
      }
    }
#endif
    ListGenTot.push_back(eGenB);
  }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
  os << "COMB: EIG_IO, we have ListGensTot, comp=" << compute_complexity_listmat(ListGenTot) << "\n";
#endif
  std::vector<MyMatrix<Tint>>  ListGenRet = ExhaustiveReductionComplexityGroupMatrix<Tint>(ListGenTot, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
  os << "COMB: EIG_IO, we have ListGensRet, comp=" << compute_complexity_listmat(ListGenRet) << "\n";
#endif
  return ListGenRet;
}

template <typename T, typename Tint, typename Tgroup>
struct IndefiniteCombinedAlgo {
private:
  std::ostream &os;
  DatabaseResultEquiStab<MyMatrix<T>, MyMatrix<Tint>> database;

private:
  std::vector<MyMatrix<Tint>>
  INDEF_FORM_AutomorphismGroup_Kernel(MyMatrix<T> const &Qmat) {
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    MicrosecondTime time;
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: Beginning of INDEF_FORM_AutomorphismGroup_Kernel with Qmat=\n";
    WriteMatrix(os, Qmat);
#endif
    AttackScheme<T> eBlock = INDEF_FORM_GetAttackScheme(Qmat, os);
    MyMatrix<T> const &mat = eBlock.mat;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: AttackScheme, eBlock.h=" << eBlock.h << "\n";
#endif
    if (eBlock.h == 0) {
      std::vector<MyMatrix<Tint>> LGen = INDEF_FORM_AutomorphismGroup_PosNeg<T, Tint, Tgroup>(mat, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
      os << "|COMB: INDEF_FORM_AutomorphismGroup_PosNeg|=" << time << "\n";
#endif
      return LGen;
    }
    if (eBlock.h == 1) {
      std::vector<MyMatrix<Tint>> LGen = LORENTZ_GetGeneratorsAutom<T, Tint, Tgroup>(mat, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
      os << "|COMB: LORENTZ_GetGeneratorsAutom|=" << time << "\n";
#endif
      return LGen;
    }
    std::vector<MyMatrix<Tint>> ListGenerators;
    auto f_insert = [&](MyMatrix<Tint> const &eGen) -> void {
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
      check_equivalence(Qmat, Qmat, eGen, "Qmat / eGen");
#endif
      ListGenerators.push_back(eGen);
    };
    ApproximateModel<T, Tint> approx =
        INDEF_FORM_GetApproximateModel<T, Tint, Tgroup>(mat, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: approx|=" << time << "\n";
#endif
    FirstNorm<T, Tint> first_norm = GetFirstNorm(approx, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: first_norm|=" << time << "\n";
#endif
    T const &X = first_norm.X;
    MyVector<Tint> const &v1 = first_norm.eVect;
    std::vector<MyMatrix<Tint>> list_gen1 = approx.GetApproximateGroup(os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: GetApproximateGroup|=" << time << "\n";
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: |list_gen1|=" << list_gen1.size() << "\n";
#endif
    for (auto &eGen : list_gen1) {
      f_insert(eGen);
    }
    std::vector<MyMatrix<Tint>> list_gen2 =
        INDEF_FORM_StabilizerVector(mat, v1);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: INDEF_FORM_StabilizerVector|=" << time << "\n";
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: |list_gen2|=" << list_gen2.size() << "\n";
#endif
    for (auto &eGen : list_gen2) {
      f_insert(eGen);
    }
    std::vector<MyVector<Tint>> the_cover = approx.GetCoveringOrbitRepresentatives(X, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: GetCoveringOrbitRepresentatives|=" << time << "\n";
#endif
    for (auto &v2 : the_cover) {
      std::optional<MyMatrix<Tint>> opt =
          INDEF_FORM_EquivalenceVector(mat, mat, v1, v2);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
      os << "|COMB: INDEF_FORM_EquivalenceVector|=" << time << "\n";
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      os << "COMB: We have opt\n";
#endif
      if (opt) {
        f_insert(*opt);
      }
    }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: |ListGenerators|=" << ListGenerators.size() << "\n";
#endif
    return ListGenerators;
  }
  std::optional<MyMatrix<Tint>>
  INDEF_FORM_TestEquivalence_Kernel(MyMatrix<T> const &Qmat1,
                                    MyMatrix<T> const &Qmat2) {
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    MicrosecondTime time;
#endif
    AttackScheme<T> eBlock1 = INDEF_FORM_GetAttackScheme(Qmat1, os);
    AttackScheme<T> eBlock2 = INDEF_FORM_GetAttackScheme(Qmat2, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_TestEquivalence_Kernel, eBlock1.h=" << eBlock1.h << " eBlock2.h=" << eBlock2.h << "\n";
#endif
    if (eBlock1.h != eBlock2.h) {
      return {};
    }
    MyMatrix<T> const &mat1 = eBlock1.mat;
    MyMatrix<T> const &mat2 = eBlock2.mat;
    if (eBlock1.h == 0) {
      std::optional<MyMatrix<Tint>> opt = INDEF_FORM_TestEquivalence_PosNeg<T, Tint>(mat1, mat2, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
      os << "|COMB: INDEF_FORM_TestEquivalence_PosNeg|=" << time << "\n";
#endif
      return opt;
    }
    if (eBlock1.h == 1) {
#ifdef TRACK_INFO_INDEFINITE_COMBINED_ALGORITHMS
      std::vector<MyMatrix<T>> pair_mat{mat1, mat2};
      write_matrix_group(pair_mat, "input_lorentz_test_equivalence_matrices");
#endif
      std::optional<MyMatrix<Tint>> opt = LORENTZ_TestEquivalenceMatrices<T, Tint, Tgroup>(mat1, mat2, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
      os << "|COMB: LORENTZ_TestEquivalenceMatrices|=" << time << "\n";
#endif
      return opt;
    }
    ApproximateModel<T, Tint> approx1 =
        INDEF_FORM_GetApproximateModel<T, Tint, Tgroup>(mat1, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: INDEF_FORM_GetApproximateModel, approx1|=" << time << "\n";
#endif
    FirstNorm<T, Tint> first_norm1 = GetFirstNorm(approx1, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: GetFirstNorm|=" << time << "\n";
#endif
    T const &X = first_norm1.X;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: first_norm1.X=" << first_norm1.X << "\n";
#endif
    MyVector<Tint> const &v1 = first_norm1.eVect;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: mat2=\n";
    WriteMatrix(os, mat2);
#endif
    ApproximateModel<T, Tint> approx2 =
        INDEF_FORM_GetApproximateModel<T, Tint, Tgroup>(mat2, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: INDEF_FORM_GetApproximateModel, approx2|=" << time << "\n";
#endif
    std::vector<MyVector<Tint>> ListCand2 =
        approx2.GetCoveringOrbitRepresentatives(X, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: GetCoveringOrbitRepresentatives|=" << time << "\n";
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_TestEquivalence_Kernel |ListCand2|=" << ListCand2.size() << "\n";
#endif
    for (auto &v2 : ListCand2) {
      std::optional<MyMatrix<Tint>> opt =
          INDEF_FORM_EquivalenceVector(mat1, mat2, v1, v2);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
      os << "|COMB: INDEF_FORM_EquivalenceVector|=" << time << "\n";
#endif
      if (opt) {
        MyMatrix<Tint> const& equiv = *opt;
        std::vector<MyMatrix<Tint>> l_gen1 = approx1.GetApproximateGroup(os);
        std::vector<MyMatrix<Tint>> l_gen1_red = ExhaustiveReductionComplexityGroupMatrix(l_gen1, os);
        std::vector<MyMatrix<Tint>> l_gen2 = approx2.GetApproximateGroup(os);
        std::vector<MyMatrix<Tint>> l_gen2_red = ExhaustiveReductionComplexityGroupMatrix(l_gen2, os);
        size_t max_iter = 1000;
        DoubleCosetSimplification<Tint> dcs = ExhaustiveMatrixDoubleCosetSimplifications(equiv, l_gen2_red, l_gen1_red, max_iter, os);
        MyMatrix<Tint> const& equiv_red = dcs.d_cos_red;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        os << "COMB: INDEF_FORM_TestEquivalence_Kernel equiv=\n";
        WriteMatrix(os, equiv);
        os << "COMB: INDEF_FORM_TestEquivalence_Kernel equiv_red=\n";
        WriteMatrix(os, equiv_red);
#endif
        return equiv_red;
      }
    }
    return {};
  }
  std::vector<MyMatrix<Tint>>
  INDEF_FORM_AutomorphismGroup_FullDim(MyMatrix<T> const &Qmat) {
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    MicrosecondTime time;
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: Beginning of INDEF_FORM_AutomorphismGroup_FullDim with "
          "Qmat=\n";
    WriteMatrix(os, Qmat);
#endif
    int dim = Qmat.rows();
    if (dim == 0) {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      os << "COMB: dimension 0, producing 0 generator\n";
#endif
      return {};
    }
    if (dim == 1) {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      os << "COMB: dimension 1, producing 1 generator\n";
#endif
      MyMatrix<Tint> eGen = - IdentityMat<Tint>(dim);
      return {eGen};
    }
    ResultReduction<T, Tint> ResRed =
      ApproxCanonicalIndefiniteForm<T, Tint>(Qmat, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: ApproxCanonicalIndefiniteForm(ResRed)|=" << time << "\n";
#endif
    MyMatrix<T> const &QmatRed = ResRed.Mred;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_AutomorphismGroup_FullDim, QmatRed=\n";
    WriteMatrix(os, QmatRed);
#endif
    MyMatrix<Tint> const &B = ResRed.B;
    MyMatrix<Tint> Binv = Inverse(B);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_AutomorphismGroup_FullDim, comp(Qmat)=" << compute_complexity_matrix(Qmat) << "\n";
    os << "COMB: INDEF_FORM_AutomorphismGroup_FullDim, comp(ResRed.B)=" << compute_complexity_matrix(ResRed.B) << "\n";
    os << "COMB: INDEF_FORM_AutomorphismGroup_FullDim, comp(Binv)=" << compute_complexity_matrix(Binv) << "\n";
    os << "COMB: INDEF_FORM_AutomorphismGroup_FullDim, comp(QmatRed)=" << compute_complexity_matrix(QmatRed) << "\n";
#endif
    auto get_stab = [&]() -> std::vector<MyMatrix<Tint>> {
      std::optional<std::vector<MyMatrix<Tint>>> opt =
          database.attempt_stabilizer(QmatRed);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
      os << "|COMB: database.attempt_stabilizer|=" << time << "\n";
#endif
      if (opt) {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        os << "COMB: INDEF_FORM_AutomorphismGroup_FullDim, successful database "
              "call\n";
#endif
        return *opt;
      } else {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        os << "COMB: INDEF_FORM_AutomorphismGroup_FullDim, not in database, "
              "recomputing\n";
#endif
        std::vector<MyMatrix<Tint>> ListGen =
            INDEF_FORM_AutomorphismGroup_Kernel(QmatRed);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
        os << "|COMB: INDEF_FORM_AutomorphismGroup_Kernel|=" << time << "\n";
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        os << "COMB: INDEF_FORM_AutomorphismGroup_FullDim, We have ListGen\n";
#endif
        database.insert_stab(QmatRed, ListGen);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        os << "COMB: INDEF_FORM_AutomorphismGroup_FullDim, After database "
              "insert\n";
#endif
        return ListGen;
      }
    };
    std::vector<MyMatrix<Tint>> LGenFinal;
    for (auto &eGen : get_stab()) {
      MyMatrix<Tint> NewGen = Binv * eGen * B;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
      check_equivalence(Qmat, Qmat, NewGen, "Qmat / NewGen");
#endif
      LGenFinal.push_back(NewGen);
    }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_AutomorphismGroup_FullDim, We have comp(LGenFinal)=" << compute_complexity_listmat(LGenFinal) << "\n";
#endif
    return LGenFinal;
  }
  std::optional<MyMatrix<Tint>>
  INDEF_FORM_TestEquivalence_FullDim(MyMatrix<T> const &Qmat1,
                                     MyMatrix<T> const &Qmat2) {
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    MicrosecondTime time;
#endif
    if (Qmat1.rows() != Qmat2.rows()) {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      os << "COMB: INDEF_FORM_TestEquivalance_FullDim, different dimensions\n";
#endif
      return {};
    }
    int dim = Qmat1.rows();
    if (dim == 0) {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      os << "COMB: INDEF_FORM_TestEquivalence_FullDim, 0 dimensional, all isomorphic\n";
#endif
      MyMatrix<Tint> eGen = IdentityMat<Tint>(0);
      return eGen;
    }
    if (dim == 1) {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      os << "COMB: INDEF_FORM_TestEquivalence_FullDim, 1 dimensional, trivial to resolve\n";
#endif
      if (Qmat1(0,0) == Qmat2(0,0)) {
        MyMatrix<Tint> eGen = IdentityMat<Tint>(1);
        return eGen;
      } else {
        return {};
      }
    }
    ResultReduction<T, Tint> res1 =
      ApproxCanonicalIndefiniteForm<T, Tint>(Qmat1, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: ApproxCanonicalIndefiniteForm(Qmat1)|=" << time << "\n";
#endif
    ResultReduction<T, Tint> res2 =
      ApproxCanonicalIndefiniteForm<T, Tint>(Qmat2, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: ApproxCanonicalIndefiniteForm(Qmat2)|=" << time << "\n";
#endif
    MyMatrix<T> const &QmatRed1 = res1.Mred;
    MyMatrix<T> const &QmatRed2 = res2.Mred;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_TestEquivalence_FullDim, QmatRed1=\n";
    WriteMatrix(os, QmatRed1);
    os << "COMB: INDEF_FORM_TestEquivalence_FullDim, QmatRed2=\n";
    WriteMatrix(os, QmatRed2);
#endif
    auto get_equi = [&]() -> std::optional<MyMatrix<Tint>> {
      std::optional<std::optional<MyMatrix<Tint>>> opt =
          database.attempt_equiv(QmatRed1, QmatRed2);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
      os << "|COMB: database.attempt_equiv|=" << time << "\n";
#endif
      if (opt) {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        os << "COMB: INDEF_FORM_TestEquivalence_FullDim, successful database "
              "call\n";
#endif
        return *opt;
      } else {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        os << "COMB: INDEF_FORM_TestEquivalence_FullDim, not in database\n";
#endif
        std::optional<MyMatrix<Tint>> optB =
            INDEF_FORM_TestEquivalence_Kernel(QmatRed1, QmatRed2);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
        os << "|COMB: INDEF_FORM_TestEquivalence_Kernel|=" << time << "\n";
#endif
        database.insert_equi(QmatRed1, QmatRed2, optB);
        return optB;
      }
    };
    std::optional<MyMatrix<Tint>> opt = get_equi();
    if (opt) {
      MyMatrix<Tint> const &eEquiv = *opt;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
      check_equivalence(QmatRed1, QmatRed2, eEquiv, "QmatRed1 / QmatRed2 / eEquiv");
#endif
      MyMatrix<Tint> NewEquiv = Inverse(res2.B) * eEquiv * res1.B;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
      check_equivalence(Qmat1, Qmat2, NewEquiv, "Qmat1 / Qmat2 / NewEquiv");
#endif
      return NewEquiv;
    } else {
      return {};
    }
  }
  //
  // Specific functions f_stab, f_equiv
  //
  std::vector<MyMatrix<Tint>>
  f_stab(INDEF_FORM_Rec_IsotropicKplane<T, Tint> const &eRec, SeqDims const& sd) {
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    MicrosecondTime time;
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_stab, start\n";
#endif
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    if (k_dim(sd) != eRec.PlaneExpr.rows()) {
      std::cerr << "The input of sd is not of the correct dimension\n";
      throw TerminalException{1};
    }
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    write_seq_dims(sd, "sd(f_stab)", os);
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: Before, GRPred = INDEF_FORM_AutomorphismGroup, comp(eRec.QmatRed)=" << compute_complexity_matrix(eRec.QmatRed) << "\n";
#endif
    std::vector<MyMatrix<Tint>> GRPred =
        INDEF_FORM_AutomorphismGroup(eRec.QmatRed);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_stab, We have |GRPred|=" << GRPred.size() << " |QmatRed|=" << eRec.QmatRed.rows() << "\n";
    os << "COMB: f_stab, We have dimCompl=" << eRec.dimCompl << " the_dim=" << eRec.the_dim << " k=" << eRec.PlaneExpr.rows() << "\n";
    //    os << "COMB: f_stab, we have GRPred=\n";
    //    WriteListMatrix(os, GRPred);
#endif
    std::vector<MyMatrix<Tint>> ListGenTot =
      ExtendIsometryGroup_IsotropicOrth(GRPred, eRec, sd, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_stab, We have ListGenTot\n";
#endif
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: f_stab|=" << time << "\n";
#endif
    return ListGenTot;
  }
  std::optional<MyMatrix<Tint>>
  f_equiv(INDEF_FORM_Rec_IsotropicKplane<T, Tint> const &eRec1,
          INDEF_FORM_Rec_IsotropicKplane<T, Tint> const &eRec2,
          [[maybe_unused]] SeqDims const& sd) {
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    MicrosecondTime time;
#endif
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    if (k_dim(sd) != eRec1.PlaneExpr.rows()) {
      std::cerr << "The input of sd is not of the correct dimension\n";
      throw TerminalException{1};
    }
#endif
    std::optional<MyMatrix<Tint>> opt =
        INDEF_FORM_TestEquivalence(eRec1.QmatRed, eRec2.QmatRed);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const &test = *opt;
    MyMatrix<Tint> TheEquivTest = IdentityMat<Tint>(eRec2.the_dim);
    int p = eRec1.QmatRed.rows();
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        TheEquivTest(i, j) = test(i, j);
      }
    }
    MyMatrix<Tint> TheEquiv =
        eRec2.FullBasisInv * TheEquivTest * eRec1.FullBasis;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    check_equivalence(eRec1.GramMatRed, eRec2.GramMatRed, TheEquiv, "eRec1/eRec2.GramMatRed / TheEquiv");
    std::vector<MyMatrix<Tint>> ListSpaces1 =
      f_get_list_spaces(eRec1.PlaneExpr, sd, os);
    std::vector<MyMatrix<Tint>> ListSpaces2 =
      f_get_list_spaces(eRec2.PlaneExpr, sd, os);
    MyMatrix<Tint> TheEquivInv = Inverse(TheEquiv);
    for (size_t iSpace = 0; iSpace < ListSpaces1.size(); iSpace++) {
      MyMatrix<Tint> eSpace1_img = ListSpaces1[iSpace] * TheEquivInv;
      MyMatrix<Tint> const &eSpace2 = ListSpaces2[iSpace];
      if (!TestEqualitySpaces(eSpace1_img, eSpace2)) {
        std::cerr << "COMB: The space are not correctly mapped\n";
        throw TerminalException{1};
      }
    }
#endif
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: f_equiv|=" << time << "\n";
#endif
    return TheEquiv;
  }
  //
  // The function using choice integer input
  //
  size_t INDEF_FORM_Invariant_IsotropicKstuff_Reduced(
      MyMatrix<T> const &Qmat, MyMatrix<Tint> const &Plane, SeqDims const &sd) {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_Invariant_IsotropicKstuff_Kernel, start\n";
#endif
    INDEF_FORM_Rec_IsotropicKplane<T, Tint> eRec(Qmat, Plane, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_Invariant_IsotropicKstuff_Kernel, we have eRec\n";
#endif
    std::vector<MyMatrix<Tint>> GRP1 = f_stab(eRec, sd);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_Invariant_IsotropicKstuff_Kernel, we have GRP1\n";
#endif
    std::vector<MyMatrix<T>> GRP2 = eRec.MapOrthogonalSublatticeGroup(GRP1);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_Invariant_IsotropicKstuff_Kernel, we have GRP2\n";
#endif
    size_t eInvRed = INDEF_FORM_Invariant_IsotropicKplane_Raw(Qmat, Plane, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_Invariant_IsotropicKstuff_Kernel, we have "
          "eInvRed\n";
#endif
    size_t GRP_inv = GetRationalInvariant(GRP2);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_Invariant_IsotropicKstuff_Kernel, we have "
          "GRP_inv\n";
#endif
    return eInvRed + GRP_inv;
  }
  size_t INDEF_FORM_Invariant_IsotropicKstuff_Kernel(
      MyMatrix<T> const &Qmat, MyMatrix<Tint> const &Plane, SeqDims const &sd) {
    ResultReduction<T, Tint> res = IndefiniteReduction<T,Tint>(Qmat, os);
    MyMatrix<Tint> Binv = Inverse(res.B);
    MyMatrix<Tint> PlaneRed = Plane * Binv;
    return INDEF_FORM_Invariant_IsotropicKstuff_Reduced(res.Mred, PlaneRed, sd);
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_Equivalence_IsotropicKstuff_Reduced(
      MyMatrix<T> const &Qmat1, MyMatrix<T> const &Qmat2,
      MyMatrix<Tint> const &Plane1, MyMatrix<Tint> const &Plane2,
      SeqDims const &sd) {
    if (INDEF_FORM_Invariant_IsotropicKstuff_Reduced(Qmat1, Plane1, sd) !=
        INDEF_FORM_Invariant_IsotropicKstuff_Reduced(Qmat2, Plane2, sd)) {
      return {};
    }
    INDEF_FORM_Rec_IsotropicKplane<T, Tint> eRec1(Qmat1, Plane1, os);
    INDEF_FORM_Rec_IsotropicKplane<T, Tint> eRec2(Qmat2, Plane2, os);
    std::optional<MyMatrix<Tint>> opt = f_equiv(eRec1, eRec2, sd);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const &test = *opt;
    MyMatrix<Tint> Subspace1 = Inverse(test) * eRec2.NSP;
    MyMatrix<Tint> Subspace2 = eRec1.NSP;
    MyMatrix<T> Subspace1_T = UniversalMatrixConversion<T, Tint>(Subspace1);
    MyMatrix<T> Subspace2_T = UniversalMatrixConversion<T, Tint>(Subspace2);
    LORENTZ_ExtendOrthogonalIsotropicIsomorphism<T> TheRec(
        Qmat1, Subspace1_T, Qmat2, Subspace2_T, os);
    MyMatrix<T> EquivRat = TheRec.get_one_transformation();
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    MyMatrix<T> eProd = EquivRat * Qmat1 * EquivRat.transpose();
    if (eProd != Qmat2) {
      std::cerr << "COMB: The matrix EquivRat is not mapping Qmat1 to Qmat2\n";
      throw TerminalException{1};
    }
#endif
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    MyMatrix<T> EquivRatInv = Inverse(EquivRat);
    MyMatrix<T> Plane1_T = UniversalMatrixConversion<T, Tint>(Plane1);
    MyMatrix<T> Plane2_T = UniversalMatrixConversion<T, Tint>(Plane2);
    MyMatrix<T> Plane1_img = Plane1_T * EquivRatInv;
    if (!TestEqualitySpaces(Plane1_img, Plane2_T)) {
      std::cerr << "COMB: Plane1 and Plane2 should be mapped\n";
      throw TerminalException{1};
    }
#endif
    if (IsIntegralMatrix(EquivRat)) {
      // This is an Ansatz. If the matrix is integral, no more work to be done.
      MyMatrix<Tint> EquivRat_tint =
          UniversalMatrixConversion<Tint, T>(EquivRat);
      return EquivRat_tint;
    }
    std::vector<MyMatrix<Tint>> GRP1_A = f_stab(eRec1, sd);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_Equivalence_IsotropicKstuff_Kernel, We have "
          "GRP1_A\n";
#endif
    std::vector<MyMatrix<T>> GRP1_B =
        eRec1.MapOrthogonalSublatticeGroup(GRP1_A);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_Equivalence_IsotropicKstuff_Kernel, We have "
          "GRP1_B\n";
#endif
    // Find a g1 in GRP1_B such that EquivRat * g1 in GL(n,Z)
    std::optional<MyMatrix<Tint>> optB =
        MatrixIntegral_Equivalence_Bis_General<T, Tint, Tgroup>(GRP1_B,
                                                                EquivRat, os);
    if (!optB) {
      return {};
    }
    MyMatrix<Tint> const &TheRet = *optB;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    check_equivalence(Qmat1, Qmat2, TheRet, "Qmat1 / Qmat2 / TheRet");
    MyMatrix<Tint> TheRetInv = Inverse(TheRet);
    MyMatrix<Tint> Plane1_imgB = Plane1 * TheRetInv;
    if (!TestEqualitySpaces(Plane1_imgB, Plane2)) {
      std::cerr << "COMB: Plane1 should be mapped on Plane2\n";
      throw TerminalException{1};
    }
#endif
    return TheRet;
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_Equivalence_IsotropicKstuff_Kernel(
      MyMatrix<T> const &Qmat1, MyMatrix<T> const &Qmat2,
      MyMatrix<Tint> const &Plane1, MyMatrix<Tint> const &Plane2,
      SeqDims const &sd) {
    ResultReduction<T, Tint> res1 = IndefiniteReduction<T,Tint>(Qmat1, os);
    MyMatrix<Tint> B1_inv = Inverse(res1.B);
    MyMatrix<Tint> Plane1_red = Plane1 * B1_inv;
    //
    ResultReduction<T, Tint> res2 = IndefiniteReduction<T,Tint>(Qmat2, os);
    MyMatrix<Tint> B2_inv = Inverse(res2.B);
    MyMatrix<Tint> Plane2_red = Plane2 * B2_inv;
    //
    std::optional<MyMatrix<Tint>> opt =
      INDEF_FORM_Equivalence_IsotropicKstuff_Reduced(res1.Mred, res2.Mred, Plane1_red, Plane2_red, sd);
    if (opt) {
      MyMatrix<Tint> const& eEquiv = *opt;
      MyMatrix<Tint> eEquivRet = B2_inv * eEquiv * res1.B;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
      check_equivalence(Qmat1, Qmat2, eEquivRet, "Qmat1 / Qmat2 / eEquivRet");
      std::vector<MyMatrix<Tint>> ListSpaces1 =
        f_get_list_spaces(Plane1, sd, os);
      std::vector<MyMatrix<Tint>> ListSpaces2 =
        f_get_list_spaces(Plane2, sd, os);
      MyMatrix<Tint> TheEquivInv = Inverse(eEquivRet);
      for (size_t iSpace = 0; iSpace < ListSpaces1.size(); iSpace++) {
        MyMatrix<Tint> eSpace1_img = ListSpaces1[iSpace] * TheEquivInv;
        MyMatrix<Tint> const &eSpace2 = ListSpaces2[iSpace];
        if (!TestEqualitySpaces(eSpace1_img, eSpace2)) {
          std::cerr << "COMB: The space are not correctly mapped\n";
          throw TerminalException{1};
        }
      }
#endif
      return eEquivRet;
    }
    return {};
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_Stabilizer_IsotropicKstuff_Reduced(
      MyMatrix<T> const &Qmat, MyMatrix<Tint> const &Plane, SeqDims const &sd) {
    INDEF_FORM_Rec_IsotropicKplane<T, Tint> eRec(Qmat, Plane, os);
    std::vector<MyMatrix<Tint>> GRP1 = f_stab(eRec, sd);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_Stabilizer_IsotropicKstuff_Reduced, We have GRP1\n";
#endif
    std::vector<MyMatrix<T>> GRP2 = eRec.MapOrthogonalSublatticeGroup(GRP1);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_Stabilizer_IsotropicKstuff_Reduced, We have GRP2\n";
#endif
    int n = Qmat.rows();
    RetMI_S<Tint,Tgroup> ret = MatrixIntegral_Stabilizer_General<T, Tint, Tgroup>(n, GRP2, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_Stabilizer_IsotropicKstuff_Reduced, index=" << ret.index << "\n";
#endif
    return ret.LGen;
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel(
      MyMatrix<T> const &Qmat, MyMatrix<Tint> const &Plane, SeqDims const &sd) {
    int n = Qmat.rows();
    if (RankMat(Qmat) != n) {
      std::cerr << "COMB: Right now INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel requires Qmat to be full dimensional\n";
      throw TerminalException{1};
    }
    ResultReduction<T, Tint> res = IndefiniteReduction<T,Tint>(Qmat, os);
    MyMatrix<Tint> B_inv = Inverse(res.B);
    MyMatrix<Tint> Plane_red = Plane * B_inv;
    //
    std::vector<MyMatrix<Tint>> LGen =
      INDEF_FORM_Stabilizer_IsotropicKstuff_Reduced(res.Mred, Plane_red, sd);
    std::vector<MyMatrix<Tint>> LGenRet;
    for (auto & eGen : LGen) {
      MyMatrix<Tint> eGenRet = B_inv * eGen * res.B;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
      check_equivalence(Qmat, Qmat, eGenRet, "Qmat / eGenRet");
      std::vector<MyMatrix<Tint>> ListSpaces =
        f_get_list_spaces(Plane, sd, os);
      for (size_t iSpace = 0; iSpace < ListSpaces.size(); iSpace++) {
        MyMatrix<Tint> const &eSpace = ListSpaces[iSpace];
        MyMatrix<Tint> eSpace_img = ListSpaces[iSpace] * eGenRet;
        if (!TestEqualitySpaces(eSpace_img, eSpace)) {
          std::cerr << "COMB: The space are not correctly mapped\n";
          throw TerminalException{1};
        }
      }
#endif
      LGenRet.push_back(eGenRet);
    }
    return LGenRet;
  }
  std::vector<MyMatrix<T>>
  INDEF_FORM_RightCosets_IsotropicKstuff_Reduced(MyMatrix<T> const &Qmat,
                                                 MyMatrix<Tint> const &ePlane,
                                                 SeqDims const &sd) {
    // We have the following groups:
    // -- H1 The group stabilizing ePlane (which is integral)
    // -- G the group stabilizing ePlane^{perp} and its mapping to the full group
    //    (which is rational)
    // -- The intersection H = G \cap GL_n(Z)
    // The right cosets that are computed are the ones of G / H.
    INDEF_FORM_Rec_IsotropicKplane<T, Tint> eRec(Qmat, ePlane, os);
    std::vector<MyMatrix<Tint>> GRP1 = f_stab(eRec, sd);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_RightCosets_IsotropicKstuff_Kernel, We have GRP1\n";
#endif
    std::vector<MyMatrix<T>> GRP2 = eRec.MapOrthogonalSublatticeGroup(GRP1);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_RightCosets_IsotropicKstuff_Kernel, We have GRP2\n";
#endif
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    for (auto &eGen : GRP2) {
      check_equivalence(Qmat, Qmat, eGen, "Qmat / GRP2");
      MyMatrix<T> ePlane_T = UniversalMatrixConversion<T, Tint>(ePlane);
      MyMatrix<T> ePlaneImg = ePlane_T * eGen;
      int k = ePlane.rows();
      for (int u = 0; u < k; u++) {
        MyVector<T> eLine = GetMatrixRow(ePlaneImg, u);
        std::optional<MyVector<T>> opt = SolutionIntMat(ePlane_T, eLine);
        if (!opt) {
          std::cerr << "COMB: ePlane should be left invariant by eGen\n";
          throw TerminalException{1};
        }
      }
    }
#endif
    int n = Qmat.rows();
    std::vector<MyMatrix<T>> ListRightCoset =
        MatrixIntegral_RightCosets_General<T, Tgroup>(n, GRP2, os);
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    for (auto &eCos : ListRightCoset) {
      check_equivalence(Qmat, Qmat, eCos, "Qmat / eCos");
    }
#endif
    return ListRightCoset;
  }
  std::vector<MyMatrix<T>>
  INDEF_FORM_RightCosets_IsotropicKstuff_Kernel(MyMatrix<T> const &Qmat,
                                                MyMatrix<Tint> const &ePlane,
                                                SeqDims const &sd) {
    int n = Qmat.rows();
    if (RankMat(Qmat) != n) {
      std::cerr << "COMB: Right now INDEF_FORM_StabilizerVector requires Qmat to be full dimensional\n";
      throw TerminalException{1};
    }
    // Apply the indefinite reduction algorithm
    ResultReduction<T,Tint> res = IndefiniteReduction<T,Tint>(Qmat, os);
    MyMatrix<Tint> Binv = Inverse(res.B);
    MyMatrix<T> B_T = UniversalMatrixConversion<T,Tint>(res.B);
    MyMatrix<T> B_Tinv = Inverse(B_T);
    //
    MyMatrix<Tint> ePlaneRed = ePlane * Binv;
    std::vector<MyMatrix<T>> LCos =
      INDEF_FORM_RightCosets_IsotropicKstuff_Reduced(res.Mred, ePlaneRed, sd);
    std::vector<MyMatrix<T>> LCosRet;
    for (auto & eCos : LCos) {
      MyMatrix<T> eCosRet = B_Tinv * eCos * B_T;
      LCosRet.push_back(eCosRet);
    }
    return LCosRet;
  }
  std::vector<MyVector<Tint>>
  INDEF_FORM_GetOrbitRepresentative_Reduced(MyMatrix<T> const &Q, T const &X) {
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    MicrosecondTime time;
#endif
    AttackScheme<T> eBlock = INDEF_FORM_GetAttackScheme(Q, os);
    int n = Q.rows();
    auto orbit_decomposition = [&](std::vector<MyVector<Tint>> const &l_vect)
        -> std::vector<MyVector<Tint>> {
      std::vector<std::vector<MyVector<Tint>>> ListGroupedRepr;
      auto f_insert = [&](MyVector<Tint> fRepr) -> void {
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
        T fNorm = EvaluationQuadForm<T, Tint>(Q, fRepr);
        if (fNorm != X) {
          std::cerr << "COMB: fNorm=" << fNorm << " X=" << X << "\n";
          std::cerr << "COMB: The norm is inconsistent\n";
          throw TerminalException{1};
        }
#endif
        for (size_t u=0; u<ListGroupedRepr.size(); u++) {
          MyVector<Tint> const& eRepr = ListGroupedRepr[u][0];
          std::optional<MyMatrix<Tint>> opt =
              INDEF_FORM_EquivalenceVector(Q, Q, eRepr, fRepr);
          if (opt) {
            ListGroupedRepr[u].push_back(fRepr);
            return;
          }
        }
        ListGroupedRepr.push_back({fRepr});
      };
      for (auto &eVect : l_vect) {
        f_insert(eVect);
      }
      std::vector<MyVector<Tint>> ListRepr;
      for (auto & eGroupedRepr : ListGroupedRepr) {
        T min_norm(0);
        MyVector<Tint> ChosenRepr;
        bool IsFirst = true;
        for (auto & eRepr : eGroupedRepr) {
          T norm(0);
          for (int i=0; i<n; i++) {
            norm += T_abs(eRepr(i));
          }
          if (IsFirst) {
            ChosenRepr = eRepr;
            min_norm = norm;
            IsFirst = false;
          } else {
            if (norm < min_norm) {
              ChosenRepr = eRepr;
              min_norm = norm;
            }
          }
        }
        ListRepr.push_back(ChosenRepr);
      }
      return ListRepr;
    };
    if (eBlock.h == 0) {
      std::vector<MyVector<Tint>> LRepr = INDEF_FORM_GetOrbitRepresentative_PosNeg<T, Tint, Tgroup>(Q, X, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
      os << "|COMB: INDEF_FORM_GetOrbitRepresentative_PosNeg|=" << time << "\n";
#endif
      return LRepr;
    }
    if (eBlock.h == 1) {
      MyMatrix<T> const &mat = eBlock.mat;
      T Xcall = X * eBlock.sign;
      std::vector<MyVector<Tint>> ListRepr =
          LORENTZ_GetOrbitRepresentative<T, Tint, Tgroup>(mat, Xcall, os);
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
      std::vector<MyVector<Tint>> ListReprRed = orbit_decomposition(ListRepr);
      if (ListReprRed.size() != ListRepr.size()) {
        std::cerr << "COMB: |ListRepr|=" << ListRepr.size()
                  << " |ListReprRed|=" << ListReprRed.size() << "\n";
        std::cerr << "COMB: The ListRepr should have been reduced\n";
        throw TerminalException{1};
      }
#endif
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
      os << "|COMB: LORENTZ_GetOrbitRepresentative|=" << time << "\n";
#endif
      return ListRepr;
    }
    ApproximateModel<T, Tint> approx =
        INDEF_FORM_GetApproximateModel<T, Tint, Tgroup>(Q, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: INDEF_FORM_GetApproximateModel|=" << time << "\n";
#endif
    std::vector<MyVector<Tint>> ListCand =
        approx.GetCoveringOrbitRepresentatives(X, os);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: approx.GetCoveringOrbitRepresentatives|=" << time << "\n";
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: GetCoveringOrbitRepresentatives, |ListCand|="
       << ListCand.size() << "\n";
#endif
    std::vector<MyMatrix<Tint>> ListGenApprox = approx.GetApproximateGroup(os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: GetCoveringOrbitRepresentatives, |ListGenApprox|="
       << ListGenApprox.size() << "\n";
#endif
    std::vector<MyVector<Tint>> ListCandSimp = ExhaustiveVectorSimplifications(ListCand, ListGenApprox);

    std::vector<MyVector<Tint>> ListRepr = orbit_decomposition(ListCandSimp);
#ifdef TIMINGS_INDEFINITE_COMBINED_ALGORITHMS
    os << "|COMB: orbit_decomposition|=" << time << "\n";
#endif
    return ListRepr;
  }
  std::vector<MyMatrix<Tint>>
  INDEF_FORM_StabilizerVector_Reduced(MyMatrix<T> const &Qmat,
                                      MyVector<Tint> const &v) {
    if (RankMat(Qmat) != Qmat.rows()) {
      std::cerr << "COMB: Right now the matrix Qmat should be full dimensional\n";
      throw TerminalException{1};
    }
    int n = Qmat.rows();
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: Beginning of INDEF_FORM_StabilizerVector_Reduced\n";
#endif
    INDEF_FORM_GetVectorStructure<T, Tint> eRec(Qmat, v, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: Before, GRP1 = INDEF_FORM_AutomorphismGroup, comp(eRec.GramMatRed)=" << compute_complexity_matrix(eRec.GramMatRed) << "\n";
#endif
    std::vector<MyMatrix<Tint>> GRP1 =
        INDEF_FORM_AutomorphismGroup(eRec.GramMatRed);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_StabilizerVector_Reduced, We have GRP1\n";
#endif
    std::vector<MyMatrix<T>> GRP2_T = eRec.MapOrthogonalSublatticeGroup(GRP1);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_StabilizerVector_Reduced, We have GRP2_T\n";
#endif
    RetMI_S<Tint,Tgroup> ret =
        MatrixIntegral_Stabilizer_General<T, Tint, Tgroup>(n, GRP2_T, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_StabilizerVector_Reduced, index=" << ret.index << "\n";
#endif
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    for (auto &eMat : ret.LGen) {
      check_equivalence(Qmat, Qmat, eMat, "Qmat / eMat");
      MyVector<Tint> v_eMat = eMat.transpose() * v;
      if (v_eMat != v) {
        std::cerr << "COMB: The matrix eMat does not preserves v\n";
        throw TerminalException{1};
      }
    }
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_StabilizerVector_Reduced, returning\n";
#endif
    return ret.LGen;
  }
  std::vector<MyMatrix<Tint>>
  INDEF_FORM_StabilizerVector_NotReduced(MyMatrix<T> const &Q,
                              MyVector<Tint> const &v) {
    ResultReduction<T,Tint> res = IndefiniteReduction<T,Tint>(Q, os);
    MyMatrix<Tint> Binv = Inverse(res.B);
    //
    MyVector<Tint> vRed = Binv.transpose() * v;
    std::vector<MyMatrix<Tint>> LGen =
      INDEF_FORM_StabilizerVector_Reduced(res.Mred, vRed);
    std::vector<MyMatrix<Tint>> LGenRet;
    for (auto & eGen : LGen) {
      MyMatrix<Tint> eGenRet = Binv * eGen * res.B;
      LGenRet.push_back(eGenRet);
    }
    return LGenRet;
  }
  std::optional<MyMatrix<Tint>>
  INDEF_FORM_EquivalenceVector_Reduced(MyMatrix<T> const &Q1, MyMatrix<T> const &Q2,
                                       MyVector<Tint> const &v1,
                                       MyVector<Tint> const &v2) {
    if (INDEF_FORM_InvariantVector(Q1, v1, os) !=
        INDEF_FORM_InvariantVector(Q2, v2, os)) {
      return {};
    }
    if (RankMat(Q1) != Q1.rows()) {
      std::cerr << "COMB: We need Q1 to be a square matrix\n";
      throw TerminalException{1};
    }
    T eNorm = EvaluationQuadForm<T, Tint>(Q1, v1);
    INDEF_FORM_GetVectorStructure<T, Tint> eRec1(Q1, v1, os);
    INDEF_FORM_GetVectorStructure<T, Tint> eRec2(Q2, v2, os);
    std::optional<MyMatrix<Tint>> opt =
        INDEF_FORM_TestEquivalence(eRec1.GramMatRed, eRec2.GramMatRed);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const &test = *opt;
    MyMatrix<T> test_T = UniversalMatrixConversion<T, Tint>(test);
    auto iife_equiv_rat = [&]() -> MyMatrix<T> {
      if (eNorm != 0) {
        return eRec2.PmatInv * ExpandMatrix(test_T) * eRec1.Pmat;
      } else {
        MyMatrix<T> Subspace1 = Inverse(test_T) * eRec2.NSP_T;
        MyMatrix<T> Subspace2 = eRec1.NSP_T;
        std::optional<MyMatrix<T>> opt =
            LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1(Q1, Subspace1, Q2,
                                                              Subspace2);
        MyMatrix<T> EquivRat = unfold_opt(opt, "failed to find EquivRat");
        MyMatrix<T> EquivRatInv = Inverse(EquivRat);
        MyVector<T> v1_inv_equivrat = EquivRatInv.transpose() * eRec1.v_T;
        if (v1_inv_equivrat == eRec2.v_T) {
          return EquivRat;
        } else {
          return -EquivRat;
        }
      }
    };
    MyMatrix<T> EquivRat = iife_equiv_rat();
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    MyMatrix<T> EquivRatInv = Inverse(EquivRat);
    MyVector<T> v1_T = UniversalVectorConversion<T, Tint>(v1);
    MyVector<T> v2_T = UniversalVectorConversion<T, Tint>(v2);
    MyVector<T> v_EquivRatInv = EquivRatInv.transpose() * v1_T;
    if (v_EquivRatInv != v2_T) {
      std::cerr << "COMB: The vector v1 is not mapped correctly\n";
      throw TerminalException{1};
    }
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: Before, GRP1_A = INDEF_FORM_AutomorphismGroup, comp(eRec1.GramMatRed)=" << compute_complexity_matrix(eRec1.GramMatRed) << "\n";
#endif
    std::vector<MyMatrix<Tint>> GRP1_A =
        INDEF_FORM_AutomorphismGroup(eRec1.GramMatRed);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_EquivalenceVector, We have GRP1_A\n";
#endif
    std::vector<MyMatrix<T>> GRP1_B =
        eRec1.MapOrthogonalSublatticeGroup(GRP1_A);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_EquivalenceVector, We have GRP1_B\n";
#endif
    // Find a g1 in GRP1_B such that EquivRat * g1 in GL(n,Z)
    std::optional<MyMatrix<Tint>> optB =
        MatrixIntegral_Equivalence_Bis_General<T, Tint, Tgroup>(GRP1_B,
                                                                EquivRat, os);
    if (!optB) {
      return {};
    }
    MyMatrix<Tint> const &TheRet = *optB;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    check_equivalence(Q1, Q2, TheRet, "Q1 / Q2 / TheRet");
    MyMatrix<Tint> TheRetInv = Inverse(TheRet);
    MyVector<Tint> v_TheRetInv = TheRetInv.transpose() * v1;
    if (v_TheRetInv != v2) {
      std::cerr << "COMB: The vector v1 is not mapped correctly\n";
      throw TerminalException{1};
    }
#endif
    return TheRet;
  }
  std::vector<MyMatrix<Tint>>
  INDEF_FORM_AutomorphismGroup_Reduced(MyMatrix<T> const &Q) {
    int n = Q.rows();
    MyMatrix<T> NSP_T = SublatticeBasisReduction(NullspaceIntMat(Q), os);
    MyMatrix<Tint> NSP = UniversalMatrixConversion<Tint, T>(NSP_T);
    MyMatrix<Tint> TheCompl = SubspaceCompletionInt(NSP, n);
    MyMatrix<Tint> FullBasis = Concatenate(TheCompl, NSP);
    MyMatrix<Tint> FullBasisInv = Inverse(FullBasis);
    int p = TheCompl.rows();
    MyMatrix<T> TheCompl_T = UniversalMatrixConversion<T, Tint>(TheCompl);
    MyMatrix<T> QmatRed = TheCompl_T * Q * TheCompl_T.transpose();
    std::vector<MyMatrix<Tint>> GRPred =
        INDEF_FORM_AutomorphismGroup_FullDim(QmatRed);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_AutomorphismGroup_Reduced, We have GRPred p=" << p
       << " n=" << n << "\n";
#endif
    std::vector<MyMatrix<Tint>> GRPfull = ExtendIsometryGroup(GRPred, p, n);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_AutomorphismGroup_Reduced, We have GRPfull\n";
#endif
    std::vector<MyMatrix<Tint>> ListGenTot;
    for (auto &eGen : GRPfull) {
      MyMatrix<Tint> eGenB = FullBasisInv * eGen * FullBasis;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
      check_equivalence(Q, Q, eGenB, "Q / eGenB");
#endif
      ListGenTot.push_back(eGenB);
    }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_AutomorphismGroup_Reduced, returning comp(ListGenTot)=" << compute_complexity_listmat(ListGenTot) << "\n";
#endif
    return ListGenTot;
  }
  std::optional<MyMatrix<Tint>>
  INDEF_FORM_TestEquivalence_Reduced(MyMatrix<T> const &Q1, MyMatrix<T> const &Q2) {
    int n = Q1.rows();
    MyMatrix<T> NSP1_T = SublatticeBasisReduction(NullspaceIntMat(Q1), os);
    MyMatrix<Tint> NSP1 = UniversalMatrixConversion<Tint, T>(NSP1_T);
    MyMatrix<Tint> TheCompl1 = SubspaceCompletionInt(NSP1, n);
    MyMatrix<T> TheCompl1_T = UniversalMatrixConversion<T, Tint>(TheCompl1);
    MyMatrix<Tint> FullBasis1 = Concatenate(TheCompl1, NSP1);
    MyMatrix<T> QmatRed1 = TheCompl1_T * Q1 * TheCompl1_T.transpose();
    MyMatrix<T> NSP2_T = SublatticeBasisReduction(NullspaceIntMat(Q2), os);
    MyMatrix<Tint> NSP2 = UniversalMatrixConversion<Tint, T>(NSP2_T);
    MyMatrix<Tint> TheCompl2 = SubspaceCompletionInt(NSP2, n);
    MyMatrix<T> TheCompl2_T = UniversalMatrixConversion<T, Tint>(TheCompl2);
    MyMatrix<Tint> FullBasis2 = Concatenate(TheCompl2, NSP2);
    MyMatrix<T> QmatRed2 = TheCompl2_T * Q2 * TheCompl2_T.transpose();
    int p = TheCompl1.rows();
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: QmatRed1=\n";
    WriteMatrix(os, QmatRed1);
    os << "COMB: QmatRed2=\n";
    WriteMatrix(os, QmatRed2);
#endif
    std::optional<MyMatrix<Tint>> opt =
        INDEF_FORM_TestEquivalence_FullDim(QmatRed1, QmatRed2);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const &test = *opt;
    MyMatrix<Tint> TheEquivTest = IdentityMat<Tint>(n);
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        TheEquivTest(i, j) = test(i, j);
      }
    }
    MyMatrix<Tint> TheEquiv = Inverse(FullBasis2) * TheEquivTest * FullBasis1;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    check_equivalence(Q1, Q2, TheEquiv, "Q1 / Q2 / TheEquiv");
#endif
    return TheEquiv;
  }
  // v should be an isotropic vector of eRec.QmatRed,
  // Returns a subgroup stabilizing the isotropic vector v in eRec.
  // The code below returns the full stabilizer, but it is not really
  std::vector<MyMatrix<Tint>> f_stab_plane_v(INDEF_FORM_Rec_IsotropicKplane<T, Tint> const& eRec,
                                             MyVector<Tint> const& v,
                                             SeqDims const& sd) {
    int n = eRec.GramMatRed.rows();
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    T eNorm = EvaluationQuadForm<T, Tint>(eRec.QmatRed, v);
    if (eNorm != 0) {
      std::cerr << "COMB: v should be an isotropic vector for Q\n";
      throw TerminalException{1};
    }
#endif
    std::vector<MyMatrix<Tint>> LGenRet;
    // First part: Mapping v to -v
    LGenRet.push_back(-IdentityMat<Tint>(n));
    // Second part: the mapping to v to v.
    std::vector<MyMatrix<Tint>> LGen2 = INDEF_FORM_StabilizerVector_NotReduced(eRec.QmatRed, v);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_stab_plane_v, We have LGen2\n";
#endif
    for (auto & eGen2 : ExtendIsometryGroup_IsotropicOrth(LGen2, eRec, sd, os)) {
      LGenRet.push_back(eGen2);
    }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_stab_plane_v, step 1\n";
    os << "COMB: f_stab_plane_v, |v|=" << v.size() << "\n";
    os << "COMB: f_stab_plane_v, |eRec.TheCompl|=" << eRec.TheCompl.rows() << " / " << eRec.TheCompl.cols() << "\n";
#endif
    // Third part: The mappings from v to v + v_iso
    MyVector<Tint> v_full = eRec.TheCompl.transpose() * v;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_stab_plane_v, step 2\n";
#endif
    MyMatrix<Tint> fPlane = ConcatenateMatVec(eRec.PlaneExpr, v_full);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_stab_plane_v, step 3\n";
    os << "COMB: f_stab_plane_v, n=" << n << "\n";
    os << "COMB: f_stab_plane_v, fPlane=";
    WriteMatrix(os, fPlane);
#endif
    MyMatrix<Tint> TheCompl = SubspaceCompletionInt(fPlane, n);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_stab_plane_v, step 4\n";
#endif
    MyMatrix<Tint> FullBasis = Concatenate(fPlane, TheCompl);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_stab_plane_v, step 5\n";
    os << "COMB: f_stab_plane_v, FullBasis=";
    WriteMatrix(os, FullBasis);
    MyMatrix<T> FullBasis_T = UniversalMatrixConversion<T,Tint>(FullBasis);
    MyMatrix<T> eProd = FullBasis_T * eRec.GramMatRed * FullBasis_T.transpose();
    os << "COMB: f_stab_plane_v, eProd=";
    WriteMatrix(os, eProd);
#endif
    MyMatrix<Tint> FullBasisInv = Inverse(FullBasis);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_stab_plane_v, step 6\n";
    os << "COMB: f_stab_plane_v, FullBasisInv=";
    WriteMatrix(os, FullBasisInv);
#endif
    int k = eRec.Plane.rows();
    for (int iK=0; iK<k; iK++) {
      MyMatrix<Tint> eGenBasis = IdentityMat<Tint>(n);
      eGenBasis(k, iK) = 1;
      MyMatrix<Tint> eGen = FullBasisInv * eGenBasis * FullBasis;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      os << "COMB: f_stab_plane_v, iK=" << iK << "\n";
      os << "COMB: f_stab_plane_v, eGenBasis=\n";
      WriteMatrix(os, eGenBasis);
#endif
      LGenRet.push_back(eGen);
    }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_stab_plane_v, step 7\n";
#endif
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    SeqDims sd_ext = seq_dims_append_one(sd);
    for (auto & eGenRet: LGenRet) {
      check_equivalence(eRec.GramMatRed, eRec.GramMatRed, eGenRet, "eRec.GramMatRed / eGenRet");
      std::vector<MyMatrix<Tint>> ListSpaces =
        f_get_list_spaces(fPlane, sd_ext, os);
      for (size_t iSpace = 0; iSpace < ListSpaces.size(); iSpace++) {
        MyMatrix<Tint> const &eSpace = ListSpaces[iSpace];
        MyMatrix<Tint> eSpace_img = ListSpaces[iSpace] * eGenRet;
        if (!TestEqualitySpaces(eSpace_img, eSpace)) {
          std::cerr << "COMB: The space are not correctly mapped\n";
          throw TerminalException{1};
        }
      }
    }
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_stab_plane_v, returning LGenRet\n";
#endif
    return LGenRet;
  }
  std::vector<MyMatrix<T>> f_double_cosets(INDEF_FORM_Rec_IsotropicKplane<T, Tint> const& eRec,
                                           MyVector<Tint> const& v,
                                           SeqDims const& sd) {
    // We have the following groups:
    // -- H1 The group stabilizing ePlane (which is integral)
    // -- G the group stabilizing ePlane^{perp} and its mapping to the full group
    //    (which is rational)
    // -- The intersection H = G \cap GL_n(Z)
    // The right cosets that are computed are the ones of G / H.
    //
    // The construction should go into the following way:
    // ---We have the ePlane as an isotropic space:
    //    * We can compute the automorphism group G1 of ePlane^T.
    //    * G1 can be extended to an automorphism group G1^{ext}
    //      of R^n preserving a finite index lattice L.
    //    * We can compute the integral stabilizer G2 that will be
    //      a finite index subgroup of G1^{ext}.
    // ---We compute the orbits of vectors v for the group G1.
    //    which is the same as G1^{ext}.
    // ---What we want to obtain is the orbits for the group G2.
    //    In other words, we want an orbit splitting to occur.
    // ---That gets us the orbits for the group {ePlane,v}.
    //    * The vector v belongs to ePlane^T. Therefore the
    //      problem is within the space ePlane^T.
    //    * The stabilizer of {ePlane,v} in {ePlane,v}^T is H1.
    //    * Computing the extension H1^{ext} to R^n, the
    //      sublattice and the integral stabilization H2.
    //    * The double cosets decomposition that we are looking
    //      for is
    //      G1^{ext} = \cup_g H1^{ext} g G2
    // ---With that choice of construction, it is forced upon
    //    us that H1^{ext} be a subgroup of G1^{ext}. Therefore
    //    It is put upon us that the lattice being stabilized
    //    by both will be the same.
    // ---But how we build all that is another story.
    //    * On the one hand, H1^{ext} is a strict subgroup of G1^{ext}.
    //    * On the other hand the extension is somewhat bigger because
    //      for G1 we extend from dimension n - k and for H1 we
    //      extend from dimension n - k - 1.
    // ---What are the implications:
    //    * Computing the kernel should use the structure of the
    //      preserved spaces.
    //    * The computation of orbits of vectors could very well
    //      be followed.
    //    * The question of whether an element belongs to the
    //      group can be decided by the preservation of the lattice
    //    * We can compute the stabilizer of v explicitly in the
    //      same space as it is computed.
    //      Then map it and ob
    // ---Maybe we could reduce everything to ePlane^T. Idea
    //    goes as follows
    //    * What we can do in the iterating algorithm that gets
    //      us the integral part the following steps:
    //      std::vector<ReductionStep>
    //      with struct ReductionStep {
    //           std::vector<MyMatrix<T>> ListMatr;
    //           std::vector<Telt> ListPerm;
    //           Face eFace;
    //      }
    //    * We could build the vector std::vector<ReductionStep>
    //      and process them.
    //    * Problem: Do we have that H1 stabilize ePlane^T?
    //      If done correctly, yes, see above.
    //    * But that seems all too hackish. If we modify the
    //      ListMatr, then we migh break the homomorphisms in
    //      question.
    //
    std::vector<MyMatrix<Tint>> GRP_G_plane = f_stab(eRec, sd);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_double_cosets, we have comp(GRP_G_plane)=" << compute_complexity_listmat(GRP_G_plane) << "\n";
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    write_seq_dims(sd, "sd(f_double_cosets)", os);
    //    os << "COMB: f_double_cosets, we have GRP_G_plane\n";
    //    WriteListMatrix(os, GRP_G_plane);
#endif
    std::vector<MyMatrix<Tint>> GRP_V_plane = f_stab_plane_v(eRec, v, sd);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_double_cosets, we have comp(GRP_V_plane)=" << compute_complexity_listmat(GRP_V_plane) << "\n";
#endif
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    CheckSubgroupInclusion<Tint,Tgroup>(GRP_G_plane, GRP_V_plane, os);
#endif
    MyMatrix<T> Sublattice = eRec.ComputeInvariantSublattice(GRP_G_plane);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_double_cosets, Sublattice=\n";
    WriteMatrix(os, Sublattice);
#endif
    std::vector<MyMatrix<T>> GRP_G = eRec.MapOrthogonalSublatticeGroupUsingSublattice(GRP_G_plane, Sublattice);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_double_cosets, we have comp(GRP_G)=" << compute_complexity_listmat(GRP_G) << "\n";
#endif
    // Computation below should succeed because GRP_V_plane is a subgroup of GRP_G_plane
    std::vector<MyMatrix<T>> GRP_V = eRec.MapOrthogonalSublatticeGroupUsingSublattice(GRP_V_plane, Sublattice);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_double_cosets, we have comp(GRP_V)=" << compute_complexity_listmat(GRP_V) << "\n";
#endif
    int n = eRec.Qmat.rows();
    std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>> pair =
      MatrixIntegral_DoubleCosets_General<T, Tgroup>(n, GRP_G, GRP_V, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: f_double_cosets, we have pair\n";
#endif
    return pair.second;
  }
  std::vector<MyMatrix<Tint>>
  INDEF_FORM_GetOrbit_IsotropicKstuff_Method(MyMatrix<T> const &Qmat, int k,
                                             int const& method_generation,
                                             SeqDims const &sd) {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_GetOrbit_IsotropicKstuff_Method, beginning\n";
#endif
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    if (method_generation != METHOD_GENERATION_RIGHT_COSETS &&
        method_generation != METHOD_GENERATION_DOUBLE_COSETS) {
      std::cerr << "No valid method being provided with method_generation\n";
      throw TerminalException{1};
    }
#endif
    T ZeroNorm(0);
    std::vector<MyMatrix<Tint>> ListOrbit;
    for (auto &eVect : INDEF_FORM_GetOrbitRepresentative(Qmat, ZeroNorm)) {
      MyMatrix<Tint> ePlane = MatrixFromVector(eVect);
      ListOrbit.push_back(ePlane);
    }
    for (int iK = 2; iK <= k; iK++) {
      SeqDims sd1 = seq_dims_reduced(sd, iK - 1);
      SeqDims sd2 = seq_dims_reduced(sd, iK);
      struct PlaneInv {
        MyMatrix<Tint> ePlane;
        size_t eInv;
      };
      auto get_planeinv = [&](MyMatrix<Tint> const &Plane) -> PlaneInv {
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
        int dim = Plane.rows();
        for (int i=0; i<dim; i++) {
          MyVector<Tint> V = GetMatrixRow(Plane, i);
          if (EvaluationQuadForm<T, Tint>(Qmat, V) != 0) {
            std::cerr << "COMB: V is not isotropic\n";
            throw TerminalException{1};
          }
        }
#endif
        size_t eInv = INDEF_FORM_Invariant_IsotropicKstuff_Reduced(Qmat, Plane, sd2);
        return {Plane, eInv};
      };
      std::vector<PlaneInv> ListRecReprKplane;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      size_t n_over_generation = 0;
      size_t n_insert = 0;
#endif
      auto fInsert = [&](MyMatrix<Tint> const& fPlane) -> void {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        n_insert += 1;
#endif
        PlaneInv fRecReprKplane = get_planeinv(fPlane);
        for (auto &eRecReprKplane : ListRecReprKplane) {
          if (eRecReprKplane.eInv == fRecReprKplane.eInv) {
            std::optional<MyMatrix<Tint>> opt =
              INDEF_FORM_Equivalence_IsotropicKstuff_Reduced(Qmat, Qmat, eRecReprKplane.ePlane,
                                                             fRecReprKplane.ePlane, sd2);
            if (opt) {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
              n_over_generation += 1;
#endif
              return;
            }
          }
        }
        ListRecReprKplane.push_back(fRecReprKplane);
      };
      auto SpanRepresentatives = [&](MyMatrix<Tint> const &ePlane) {
        // Some possible improvement. Use the double cosets
        // The double coset consists in splitting an orbit x G as
        // y1 H \cup ..... yN H
        // Or in other words G = \cup_i Stab(x) y_i H
        // This is the canonical stuff.
        //
        // The groups being used.
        // G = group of rational transformation preserving L = S^{perp}
        //     and acting integrally on it.
        // Stab({S, x}) = Integral transformations preserving {L, x}
        //     considered as plane or flag and extended to acting integrally
        //     on L in the same way as G.
        // H = group of integral transformations preserving the big lattice
        //     and preserving L.
        //
        // The first method is to compute the cosets for the group G in
        // G_int. This can be computed quite efficiently with the stabilizer
        // algorithm. The spanned objects are then tested for equivalence.
        // That method is reasonable, but slow.
        //
        // What could be done?
        // * Get a correct computation of the initial set to consider.
        //   Right now, it works by kind of chance since we compute
        //   only for k=1 and k=2.
        // * So, we need to have a good strategy for computing the
        //   The problem is to compute the orbits of the vectors.
        // * The computation can be done in the space S^{perp}.
        // * From the computation in the group of the space S^{perp}
        //   We can get the orbits for the space L and the group
        //   that is an extension G_{ext} of this one.
        // * We can compute the stabilizer in the space S^{perp} which
        //   we need just above anyway. We extend it and then we have the
        //   Stab_{ext} group by the same process as above.
        // * This is what the double coset is supposed to enter into the
        //   picture.
        // * The group G is G_ext. The group U is Stab_{ext} and V is
        //   G_{ext} \cap GL(L).
        // * So, that fits with the scheme of the Stabilizer_RightCoset.
        //   When the computation is started, the group V is not known
        //   and V is computed at the same time as the double cosets are.
        // * So, the designs seems relatively clear, but what we would
        //   need is sensible examples. Those have to be finite groups
        //   - The lattice simplices are good examples, since they have
        //     some sequence of stabilizers.
        //   - The fundamental domains of reflective lorentzian groups.
        //   - The full group and the arithmetic group are always clear.
        //   - For the other group we can take a cyclic group generated
        //     by a random permutation.
        //  All, in all, that looks like a reasonable stretgy for
        //  debugging the scheme.
        // * Now our mystery: How to compute the list of orbits up to the
        //   relevant isometry group. This is we compute the isometries for
        //   QmatRed
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        os << "COMB: SpanRepresentatives, beginning\n";
#endif
        INDEF_FORM_Rec_IsotropicKplane<T, Tint> eRec(Qmat, ePlane, os);
        std::vector<MyVector<Tint>> ListOrbitF =
            INDEF_FORM_GetOrbitRepresentative(eRec.QmatRed, ZeroNorm);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        os << "COMB: |ListOrbitF|=" << ListOrbitF.size() << "\n";
#endif
        auto get_right_cosets=[&]() -> std::vector<MyMatrix<T>> {
          if (method_generation == METHOD_GENERATION_RIGHT_COSETS) {
            return INDEF_FORM_RightCosets_IsotropicKstuff_Reduced(Qmat, ePlane, sd1);
          } else {
            return {};
          }
        };
        std::vector<MyMatrix<T>> ListRightCosets = get_right_cosets();
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        os << "COMB: |ListRightCosets|=" << ListRightCosets.size() << "\n";
#endif
        std::vector<MyMatrix<Tint>> LGenStab = INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel(Qmat, ePlane, sd1);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        os << "COMB: |LGenStab|=" << LGenStab.size() << "\n";
#endif
        std::vector<MyMatrix<T>> LGenStab_T = UniversalStdVectorMatrixConversion<T,Tint>(LGenStab);
        MyMatrix<Tint> Embed = eRec.TheCompl * eRec.NSP;
        for (auto &eVect : ListOrbitF) {
          MyVector<Tint> eVectB = Embed.transpose() * eVect;
          MyVector<T> eVectB_T = UniversalVectorConversion<T, Tint>(eVectB);
          auto get_cosets=[&]() -> std::vector<MyMatrix<T>> {
            if (method_generation == METHOD_GENERATION_RIGHT_COSETS) {
              return ListRightCosets;
            } else {
              return f_double_cosets(eRec, eVect, sd1);
            }
          };
          std::vector<MyMatrix<T>> list_cosets = get_cosets();
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
          os << "COMB: |list_cosets|=" << list_cosets.size() << "\n";
#endif
          std::vector<MyMatrix<T>> list_cosetsB = ExhaustiveMatrixRightCosetSimplifications(list_cosets, LGenStab_T);
          for (auto &eCos : list_cosetsB) {
            MyVector<T> eVectC_T = eCos.transpose() * eVectB_T;
            MyVector<Tint> eVectC =
              UniversalVectorConversion<Tint, T>(eVectC_T);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
            //            os << "COMB: LGenStab=\n";
            //            WriteListMatrix(os, LGenStab);
            os << "COMB: eVectC=" << StringVectorGAP(eVectC) << "\n";
#endif
            MyVector<Tint> eVectD = ExhaustiveVectorSimplification(eVectC, LGenStab);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
            os << "COMB: eVectD=" << StringVectorGAP(eVectD) << "\n";
#endif
            MyMatrix<Tint> ePlaneB = ConcatenateMatVec(ePlane, eVectD);
            fInsert(ePlaneB);
          }
        }
      };
      for (auto &eRepr : ListOrbit) {
        SpanRepresentatives(eRepr);
      }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      os << "COMB: n_over_generation = " << n_over_generation << " n_insert=" << n_insert << "\n";
#endif
      ListOrbit.clear();
      for (auto &eRec : ListRecReprKplane) {
        ListOrbit.push_back(eRec.ePlane);
      }
    }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: |ListOrbit|=" << ListOrbit.size() << "\n";
#endif
    return ListOrbit;
  }
  std::vector<MyMatrix<Tint>>
  INDEF_FORM_GetOrbit_IsotropicKstuff_Kernel(MyMatrix<T> const &Qmat, int k,
                                             SeqDims const &sd) {
    ResultReduction<T, Tint> res = IndefiniteReduction<T,Tint>(Qmat, os);
    // Double cosets method is more efficient in principle than right cosets.
    std::vector<MyMatrix<Tint>> LOrb =
      INDEF_FORM_GetOrbit_IsotropicKstuff_Method(res.Mred, k, METHOD_GENERATION_DOUBLE_COSETS, sd);
#ifdef METHOD_COMPARISON_INDEFINITE_COMBINED_ALGORITHMS_ISOTROPIC
    std::vector<MyMatrix<Tint>> LOrbB =
      INDEF_FORM_GetOrbit_IsotropicKstuff_Method(res.Mred, k, METHOD_GENERATION_RIGHT_COSETS, sd);
    if (LOrbB.size() != LOrb.size()) {
      std::cerr << "COMB: Inconsistencies in the computation methods |LOrb|=" << LOrb.size() << " |LOrbB|=" << LOrbB.size() << "\n";
      throw TerminalException{1};
    }
#endif
    std::vector<MyMatrix<Tint>> LOrbRet;
    for (auto & eOrb : LOrb) {
      MyMatrix<Tint> eOrbRet = eOrb * res.B;
      LOrbRet.push_back(eOrbRet);
    }
    return LOrbRet;
  }




public:
  IndefiniteCombinedAlgo(std::ostream &_os) : os(_os) {}
  // Now the specific implementations
  size_t INDEF_FORM_Invariant_IsotropicKplane(MyMatrix<T> const &Q,
                                              MyMatrix<Tint> const &Plane) {
    SeqDims sd = seq_dims_plane(Plane.rows());
    return INDEF_FORM_Invariant_IsotropicKstuff_Kernel(Q, Plane, sd);
  }
  size_t INDEF_FORM_Invariant_IsotropicKflag(MyMatrix<T> const &Q,
                                             MyMatrix<Tint> const &Plane) {
    SeqDims sd = seq_dims_flag(Plane.rows());
    return INDEF_FORM_Invariant_IsotropicKstuff_Kernel(Q, Plane, sd);
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_Equivalence_IsotropicKplane(
      MyMatrix<T> const &Qmat1, MyMatrix<T> const &Qmat2,
      MyMatrix<Tint> const &Plane1, MyMatrix<Tint> const &Plane2) {
    SeqDims sd = seq_dims_plane(Plane1.rows());
    return INDEF_FORM_Equivalence_IsotropicKstuff_Kernel(
        Qmat1, Qmat2, Plane1, Plane2, sd);
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_Equivalence_IsotropicKflag(
      MyMatrix<T> const &Qmat1, MyMatrix<T> const &Qmat2,
      MyMatrix<Tint> const &Plane1, MyMatrix<Tint> const &Plane2) {
    SeqDims sd = seq_dims_flag(Plane1.rows());
    return INDEF_FORM_Equivalence_IsotropicKstuff_Kernel(
        Qmat1, Qmat2, Plane1, Plane2, sd);
  }
  std::vector<MyMatrix<Tint>>
  INDEF_FORM_Stabilizer_IsotropicKplane(MyMatrix<T> const &Q,
                                        MyMatrix<Tint> const &Plane) {
    SeqDims sd = seq_dims_plane(Plane.rows());
    return INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel(Q, Plane, sd);
  }
  std::vector<MyMatrix<Tint>>
  INDEF_FORM_Stabilizer_IsotropicKflag(MyMatrix<T> const &Q,
                                       MyMatrix<Tint> const &Plane) {
    SeqDims sd = seq_dims_flag(Plane.rows());
    return INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel(Q, Plane, sd);
  }
  std::vector<MyMatrix<T>>
  INDEF_FORM_RightCosets_IsotropicKplane(MyMatrix<T> const &Q,
                                         MyMatrix<Tint> const &Plane) {
    SeqDims sd = seq_dims_plane(Plane.rows());
    return INDEF_FORM_RightCosets_IsotropicKstuff_Kernel(Q, Plane, sd);
  }
  std::vector<MyMatrix<T>>
  INDEF_FORM_RightCosets_IsotropicKflag(MyMatrix<T> const &Q,
                                        MyMatrix<Tint> const &Plane) {
    SeqDims sd = seq_dims_flag(Plane.rows());
    return INDEF_FORM_RightCosets_IsotropicKstuff_Kernel(Q, Plane, sd);
  }
  std::vector<MyMatrix<Tint>>
  INDEF_FORM_GetOrbit_IsotropicKplane(MyMatrix<T> const &Q, int k) {
    SeqDims sd = seq_dims_plane(k);
    return INDEF_FORM_GetOrbit_IsotropicKstuff_Kernel(Q, k, sd);
  }
  std::vector<MyMatrix<Tint>>
  INDEF_FORM_GetOrbit_IsotropicKflag(MyMatrix<T> const &Q, int k) {
    SeqDims sd = seq_dims_flag(k);
    return INDEF_FORM_GetOrbit_IsotropicKstuff_Kernel(Q, k, sd);
  }
  std::vector<MyVector<Tint>>
  INDEF_FORM_GetOrbitRepresentative(MyMatrix<T> const &Q, T const &X) {
    ResultReduction<T,Tint> res = IndefiniteReduction<T,Tint>(Q, os);
    std::vector<MyVector<Tint>> LRepr =
      INDEF_FORM_GetOrbitRepresentative_Reduced(res.Mred, X);
    std::vector<MyVector<Tint>> LReprRet;
    for (auto & eRepr : LRepr) {
      MyVector<Tint> eReprRet = res.B.transpose() * eRepr;
      LReprRet.push_back(eReprRet);
    }
    return LReprRet;
  }
  std::vector<MyMatrix<Tint>>
  INDEF_FORM_StabilizerVector(MyMatrix<T> const &Q,
                              MyVector<Tint> const &v) {
    std::vector<MyMatrix<Tint>> LGen =
      INDEF_FORM_StabilizerVector_NotReduced(Q, v);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_StabilizerVector, We have LGen\n";
#endif
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
    for (auto & eGen : LGen) {
      check_equivalence(Q, Q, eGen, "Q / eGen");
      MyVector<Tint> v_img = eGen.transpose() * v;
      if (v_img != v) {
        std::cerr << "COMB: eGen does not preserve v\n";
        throw TerminalException{1};
      }
    }
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_StabilizerVector, after the checks\n";
#endif
    return LGen;
  }
  std::optional<MyMatrix<Tint>>
  INDEF_FORM_EquivalenceVector(MyMatrix<T> const &Q1, MyMatrix<T> const &Q2,
                               MyVector<Tint> const &v1,
                               MyVector<Tint> const &v2) {
    ResultReduction<T,Tint> res1 = IndefiniteReduction<T,Tint>(Q1, os);
    MyMatrix<Tint> B1_inv = Inverse(res1.B);
    ResultReduction<T,Tint> res2 = IndefiniteReduction<T,Tint>(Q2, os);
    MyMatrix<Tint> B2_inv = Inverse(res2.B);
    //
    MyVector<Tint> v1_red = B1_inv.transpose() * v1;
    MyVector<Tint> v2_red = B2_inv.transpose() * v2;
    std::optional<MyMatrix<Tint>> opt =
      INDEF_FORM_EquivalenceVector_Reduced(res1.Mred, res2.Mred, v1_red, v2_red);
    if (opt) {
      MyMatrix<Tint> const& eEquiv = *opt;
      MyMatrix<Tint> eEquivRet = B2_inv * eEquiv * res1.B;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
      check_equivalence(res1.Mred, res2.Mred, eEquiv, "res1.Mred / res2.Mred / eEquiv");
      check_equivalence(Q1, Q2, eEquivRet, "Q1 / Q2 / eEquivRet");
#endif
      return eEquivRet;
    }
    return {};
  }
  std::vector<MyMatrix<Tint>>
  INDEF_FORM_AutomorphismGroup(MyMatrix<T> const &Q) {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_AutomorphismGroup, start, comp(Q)=" << compute_complexity_matrix(Q) << "\n";
#endif
    if (Q.rows() == 0) {
      MyMatrix<Tint> eGen = IdentityMat<Tint>(0);
      return {eGen};
    }
    if (Q.rows() == 1) {
      MyMatrix<Tint> eGen = - IdentityMat<Tint>(1);
      return {eGen};
    }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_AutomorphismGroup, start |Q|=" << Q.rows() << "\n";
#endif
    ResultReduction<T,Tint> res = IndefiniteReduction<T,Tint>(Q, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_AutomorphismGroup, We have res\n";
#endif
    MyMatrix<Tint> B_inv = Inverse(res.B);
    //
    std::vector<MyMatrix<Tint>> LGen = INDEF_FORM_AutomorphismGroup_Reduced(res.Mred);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_AutomorphismGroup, comp(res.B)=" << compute_complexity_matrix(res.B) << "\n";
    os << "COMB: INDEF_FORM_AutomorphismGroup, comp(res.Mred)=" << compute_complexity_matrix(res.Mred) << "\n";
    os << "COMB: INDEF_FORM_AutomorphismGroup, comp(Q)=" << compute_complexity_matrix(Q) << "\n";
    os << "COMB: INDEF_FORM_AutomorphismGroup, comp(B_inv)=" << compute_complexity_matrix(B_inv) << "\n";
    os << "COMB: INDEF_FORM_AutomorphismGroup, We have comp(LGen)=" << compute_complexity_listmat(LGen) << "\n";
#endif
    std::vector<MyMatrix<Tint>> LGenRet;
    for (auto & eGen : LGen) {
      MyMatrix<Tint> eGenRet = B_inv * eGen * res.B;
#ifdef SANITY_CHECK_INDEFINITE_COMBINED_ALGORITHMS
      check_equivalence(Q, Q, eGenRet, "Q / eGenRet");
#endif
      LGenRet.push_back(eGenRet);
    }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_AutomorphismGroup, comp(LGenRet)=" << compute_complexity_listmat(LGenRet) << "\n";
#endif
    std::vector<MyMatrix<Tint>>  LGenRetB = ExhaustiveReductionComplexityGroupMatrix<Tint>(LGenRet, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: INDEF_FORM_AutomorphismGroup, comp(LGenRetB)=" << compute_complexity_listmat(LGenRetB) << "\n";
#endif
    return LGenRetB;
  }
  std::optional<MyMatrix<Tint>>
  INDEF_FORM_TestEquivalence(MyMatrix<T> const &Q1, MyMatrix<T> const &Q2) {
    if (Q1.rows() == 0) {
      MyMatrix<Tint> eGen = IdentityMat<Tint>(0);
      return eGen;
    }
    if (Q1.rows() == 1) {
      if (Q1(0,0) != Q2(0,0)) {
        return {};
      }
      MyMatrix<Tint> eGen = IdentityMat<Tint>(1);
      return eGen;
    }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: Q1=\n";
    WriteMatrix(os, Q1);
    os << "COMB: Q2=\n";
    WriteMatrix(os, Q2);
#endif
    ResultReduction<T,Tint> res1 = IndefiniteReduction<T,Tint>(Q1, os);
    ResultReduction<T,Tint> res2 = IndefiniteReduction<T,Tint>(Q2, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "COMB: res1.Mred=\n";
    WriteMatrix(os, res1.Mred);
    os << "COMB: res2.Mred=\n";
    WriteMatrix(os, res2.Mred);
#endif
    MyMatrix<Tint> B2_inv = Inverse(res2.B);
    //
    std::optional<MyMatrix<Tint>> opt =
      INDEF_FORM_TestEquivalence_Reduced(res1.Mred, res2.Mred);
    if (opt) {
      MyMatrix<Tint> const& eEquiv = *opt;
      MyMatrix<Tint> eEquivRet = B2_inv * eEquiv * res1.B;
      return eEquivRet;
    }
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    os << "Finally returning non-equivalence conclusion\n";
#endif
    return {};
  }
};

// clang-format off
#endif  // SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_
// clang-format on

// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_MODELS_FUNDAMENTALS_H_
#define SRC_INDEFINITE_MODELS_FUNDAMENTALS_H_

template<typename T>
struct AttackScheme {
  int h;
  MyMatrix<T> mat;
  int sign;
};

template<typename T>
AttackScheme<T> INDEF_FORM_GetAttackScheme(MyMatrix<T> const& Qmat) {
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(Qmat);
  int nbPlus=DiagInfo.nbPlus;
  int nbMinus=DiagInfo.nbMinus;
  if (nbMinus < nbPlus) {
    MyMatrix<T> RetMat = -Qmat;
    return {nbMinus, RetMat, -1};
  } else {
    return {nbPlus, Qmat, 1};
  }
}

template<typename T>
bool INDEF_FORM_IsEven(MyMatrix<T> const& Qmat) {
  if (!IsIntegralMatrix(Qmat)) {
    return false;
  }
  T two(2);
  for (int u=0; u<Qmat.rows(); u++) {
    T res = ResInt(Qmat(u,u), two);
    if (res != 0) {
      return false;
    }
  }
  return true;
}

template<typename T>
struct INDEF_InvariantQ {
  int n;
  T eDet;
  int nbPlus;
  int nbMinus;
  int nbZero;
  bool IsEven;
};

namespace std {
  template <typename T> struct hash<INDEF_InvariantQ<T>> {
    std::size_t operator()(const INDEF_InvariantQ<T> &x) const {
      auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
        seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      };
      size_t n_s = x.n;
      size_t nbPlus_s = x.nbPlus;
      size_t nbMinus_s = x.nbMinus;
      size_t nbZero_s = x.nbZero;
      size_t IsEven_s = x.IsEven;
      size_t hash_det = std::hash<T>()(x.eDet);
      //
      size_t seed = 1234;
      combine_hash(seed, n_s);
      combine_hash(seed, nbPlus_s);
      combine_hash(seed, nbMinus_s);
      combine_hash(seed, nbZero_s);
      combine_hash(seed, IsEven_s);
      combine_hash(seed, hash_det);
      return seed;
    }
  };
  // clang-format off
}  // namespace std
// clang-format on


template<typename T>
size_t INDEF_FORM_Invariant(MyMatrix<T> const& Qmat) {
  int n = Qmat.rows();
  MyMatrix<T> NSP = NullspaceIntMat(Qmat);
  MyMatrix<T> TheCompl = SubspaceCompletionInt(NSP, n);
  MyMatrix<T> GramRed = TheCompl * Qmat * TheCompl.transpose();
  T eDet = DeterminantMat(GramRed);
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(Qmat);
  int nbPlus = DiagInfo.nbPlus;
  int nbMinus = DiagInfo.nbMinus;
  int nbZero = DiagInfo.nbZero;
  bool IsEven = INDEF_FORM_IsEven(Qmat);
  INDEF_InvariantQ<T> eInv{n, eDet, nbPlus, nbMinus, nbZero, IsEven};
  return std::hash<INDEF_InvariantQ<T>>()(eInv);
}


template<typename T>
struct INDEF_InvariantQV {
  int eRank;
  T eNorm;
  T index;
  T divisor;
  size_t GramInv;
};

namespace std {
  template <typename T> struct hash<INDEF_InvariantQV<T>> {
    std::size_t operator()(const INDEF_InvariantQV<T> &x) const {
      auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
        seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      };
      size_t eRank_s = x.eRank;
      size_t eNorm_s = std::hash<T>()(x.eNorm);
      size_t index_s = std::hash<T>()(x.index);
      size_t divisor_s = std::hash<T>()(x.divisor);
      //
      size_t seed = 1234;
      combine_hash(seed, eRank_s);
      combine_hash(seed, eNorm_s);
      combine_hash(seed, index_s);
      combine_hash(seed, divisor_s);
      combine_hash(seed, x.GramInv);
      return seed;
    }
  };
  // clang-format off
}  // namespace std
// clang-format on





template<typename T, typename Tint>
size_t INDEF_FORM_InvariantVector(MyMatrix<T> const& Qmat, MyVector<Tint> const& v) {
  int eRank = RankMat(Qmat);
  T eNorm = EvaluationQuadForm<T,Tint>(Qmat, v);
  MyVector<T> v_T = UniversalVectorConversion<T,Tint>(v);
  MyVector<T> eProd = Qmat * v_T;
  T divisor = RemoveFractionVectorPlusCoeff(v_T).TheMult;
  T index = RemoveFractionVectorPlusCoeff(eProd).TheMult;
  MyMatrix<T> NSP = NullspaceIntVect(eProd);
  MyMatrix<T> GramRed = NSP * Qmat * NSP.transpose();
  size_t typeInv = INDEF_FORM_Invariant(GramRed);
  //
  INDEF_InvariantQV<T> eInv{eRank, eNorm, index, divisor, typeInv};
  return std::hash<INDEF_InvariantQV<T>>()(eInv);
}

struct InvariantIsotropic {
  int k;
  size_t eInv1;
  size_t eInv2;
};

namespace std {
  template <> struct hash<InvariantIsotropic> {
    std::size_t operator()(const InvariantIsotropic &x) const {
      auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
        seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      };
      //
      size_t seed = x.k;
      combine_hash(seed, x.eInv1);
      combine_hash(seed, x.eInv2);
      return seed;
    }
  };
  // clang-format off
}  // namespace std
// clang-format on




template<typename T, typename Tint>
size_t INDEF_FORM_Invariant_IsotropicKplane_Raw(MyMatrix<T> const& Qmat, MyMatrix<Tint> const& ePlane) {
  int k = ePlane.rows();
  MyMatrix<T> ePlane_T = UniversalMatrixConversion<T,Tint>(ePlane);
  MyMatrix<T> eProd = ePlane_T * Qmat;
  MyMatrix<T> NSP_T = NullspaceIntTrMat(eProd);
  MyMatrix<Tint> NSP = UniversalMatrixConversion<Tint,T>(NSP_T);
  int dimNSP = NSP.rows();
  MyMatrix<T> ePlaneB(k, dimNSP);
  for (int u=0; u<k; u++) {
    MyVector<Tint> eV = GetMatrixRow(ePlane, u);
    std::optional<MyVector<Tint>> opt = SolutionIntMat(NSP, eV);
    if (opt) {
      MyVector<Tint> const& fV = *opt;
      MyVector<T> fV_T = UniversalVectorConversion<T,Tint>(fV);
      AssignMatrixRow(ePlaneB, u, fV_T);
    } else {
      std::cerr << "eV should belong to the space by the virtue of being isotropic\n";
      throw TerminalException{1};
    }
  }
  MyMatrix<T> ComplBasisInNSP = SubspaceCompletionInt(ePlaneB, dimNSP);
  MyMatrix<T> NSP_sub = ComplBasisInNSP * NSP_T;
  MyMatrix<T> QmatRed = NSP_sub * Qmat * NSP_sub.transpose();
  size_t eInv1 = INDEF_FORM_Invariant(Qmat);
  size_t eInv2 = INDEF_FORM_Invariant(QmatRed);
  InvariantIsotropic eInv{k, eInv1, eInv2};
  return std::hash<InvariantIsotropic>()(eInv);
}

template<typename T>
MyMatrix<T> ExpandMatrix(MyMatrix<T> const& M) {
  int n = M.rows();
  MyMatrix<T> TheBigMat = IdentityMat<T>(n+1);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      TheBigMat(i,j) = M(i,j);
    }
  }
  return TheBigMat;
}

// clang-format off
#endif  // SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_
// clang-format on


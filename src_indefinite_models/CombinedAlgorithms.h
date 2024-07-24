// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_
#define SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_

// clang-format off
#include "ApproximateModels.h"
// clang-format on


// Used for the memoization process for isomorphism
template<typename T, typename Tint>
struct ResultIsomorphism {
  MyMatrix<T> Q1;
  MyMatrix<T> Q2;
  std::optional<MyMatrix<Tint>> result;
};

template<typename T, typename Tint>
struct ResultStabilizer {
  MyMatrix<T> Q;
  std::vector<MyMatrix<Tint>> ListGens;
};

template<typename T, typename Tint>
struct InvariantIsotropic {
};




template<typename T, typename Tint>
struct CombinedAlgo {
private:
  std::vector<ResultIsomorphism<T,Tint>> ListResultIsomorphism;
  std::vector<ResultStabilizer<T,Tint>> ListResultStabilizer;

private:
  std::vector<MyMatrix<Tint>> INDEF_FORM_AutomorphismGroup_Kernel(MyMatrix<T> const& Qmat) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_AutomorphismGroup_Memoize(MyMatrix<T> const& Qmat) {
  }

  InvariantIsotropic<T,Tint> INDEF_FORM_Invariant_IsotropicKplane(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
  }
  InvariantIsotropic<T,Tint> INDEF_FORM_Invariant_IsotropicKflag(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
  }

public:
  std::vector<MyVector<Tint>> INDEF_FORM_GetOrbitRepresentative(MyMatrix<T> const& Qmat) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_StabilizerVector(MyMatrix<T> const& Qmat) {
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_EquivalenceVector(MyMatrix<T> const& Q1, MyMatrix<T> const& Q2, MyVector<Tint> const& v1, MyVector<Tint> const& v2) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_StabilizerVector(MyMatrix<T> const& Q, MyVector<Tint> const& v) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_AutomorphismGroup(MyMatrix<T> const& Q) {
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_TestEquivalence(MyMatrix<T> const& Q1, MyMatrix<T> const& Q2) {
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_Equivalence_IsotropicKplane(MyMatrix<T> const& Q1, MyMatrix<T> const& Q2, MyMatrix<Tint> const& Plane1, MyMatrix<Tint> const& Plane2) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_Stabilizer_IsotropicKplane(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_RightCosets_IsotropicKplane(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_GetOrbit_IsotropicKplane(MyMatrix<T> const& Q) {
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_Equivalence_IsotropicKflag(MyMatrix<T> const& Q1, MyMatrix<T> const& Q2, MyMatrix<Tint> const& Plane1, MyMatrix<Tint> const& Plane2) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_Stabilizer_IsotropicKflag(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_RightCosets_IsotropicKflag(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_GetOrbit_IsotropicKflag(MyMatrix<T> const& Q) {
  }

};

// clang-format off
#endif  // SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_
// clang-format on

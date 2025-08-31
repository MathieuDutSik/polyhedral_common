// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_EUTACTICITY_H_
#define SRC_PERFECT_EUTACTICITY_H_

// clang-format off
#include "MatrixGroup.h"
#include "Positivity.h"
#include <optional>
// clang-format on

#ifdef DEBUG
#define DEBUG_EUTACTICITY
#endif

// Test if a quadratic form is eutactic
// A quadratic form is eutactic if the inverse of eGram can be expressed
// as a positive linear combination of the outer products v ⊗ v of shortest vectors
template <typename T>
bool IsEutactic(MyMatrix<T> const &eGram, MyMatrix<T> const& SHV_T, std::ostream& os) {
#ifdef DEBUG_EUTACTICITY
  os << "EUTACTIC: IsEutactic, beginning\n";
#endif
  int n = eGram.rows();
  if (eGram.cols() != n) {
    std::cerr << "Error: eGram must be square\n";
    return false;
  }
  if (SHV_T.cols() != n) {
    std::cerr << "Error: SHV_T must have n columns\n";
    return false;
  }

  int nbSHV = SHV_T.rows();
  if (nbSHV == 0) {
#ifdef DEBUG_EUTACTICITY
    os << "EUTACTIC: No shortest vectors provided\n";
#endif
    return false;
  }

  // Dimension of the space of symmetric matrices
  int dimSymm = n * (n + 1) / 2;

#ifdef DEBUG_EUTACTICITY
  os << "EUTACTIC: n=" << n << " nbSHV=" << nbSHV << " dimSymm=" << dimSymm << "\n";
#endif

  // Compute the inverse of eGram
  MyMatrix<T> eGramInv;
  try {
    eGramInv = Inverse(eGram);
  } catch (...) {
    std::cerr << "Error: Cannot compute inverse of eGram\n";
    return false;
  }

  // Convert inverse to vector form
  MyVector<T> eGramInvVec = SymmetricMatrixToVector(eGramInv);

#ifdef DEBUG_EUTACTICITY
  os << "EUTACTIC: eGram inverse computed and vectorized\n";
#endif

  // Create matrix of outer products v ⊗ v for each shortest vector v
  MyMatrix<T> OuterProductMat(nbSHV, dimSymm);

  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<T> v = GetMatrixRow(SHV_T, iSHV);

    // Compute outer product v ⊗ v as a symmetric matrix
    MyMatrix<T> OuterProd(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        OuterProd(i, j) = v(i) * v(j);
      }
    }

    // Convert symmetric matrix to vector form
    MyVector<T> OuterProdVec = SymmetricMatrixToVector(OuterProd);
    AssignMatrixRow(OuterProductMat, iSHV, OuterProdVec);
  }

#ifdef DEBUG_EUTACTICITY
  os << "EUTACTIC: OuterProductMat constructed, size: " << OuterProductMat.rows() << "x" << OuterProductMat.cols() << "\n";
#endif

  // Solve the linear system: OuterProductMat^T * coeffs = eGramInvVec
  // We want to find coefficients such that eGramInv = sum(coeffs[i] * v_i ⊗ v_i)
  MyMatrix<T> OuterProductMatT = TransposedMat(OuterProductMat);

  // Check if the system is solvable
  MyMatrix<T> AugmentedMat = ConcatenateMatVec(OuterProductMatT, eGramInvVec);
  int rankA = RankMat(OuterProductMatT);
  int rankAug = RankMat(AugmentedMat);

#ifdef DEBUG_EUTACTICITY
  os << "EUTACTIC: rankA=" << rankA << " rankAug=" << rankAug << "\n";
#endif

  if (rankA != rankAug) {
    // System is inconsistent - eGramInv cannot be expressed as linear combination
#ifdef DEBUG_EUTACTICITY
    os << "EUTACTIC: System is inconsistent - not eutactic\n";
#endif
    return false;
  }

  // Try to solve the system to get coefficients
  std::optional<MyVector<T>> opt_coeffs = SolutionMatNonnegative(OuterProductMatT, eGramInvVec);

  if (!opt_coeffs) {
#ifdef DEBUG_EUTACTICITY
    os << "EUTACTIC: No nonnegative solution found - not eutactic\n";
#endif
    return false;
  }

  MyVector<T> coeffs = *opt_coeffs;

  // Check if all coefficients are positive (not just nonnegative)
  bool all_positive = true;
  T tolerance = T(1)/T(1000000); // Small tolerance for numerical errors

  for (int i = 0; i < coeffs.size(); i++) {
    if (coeffs(i) <= tolerance) {
      all_positive = false;
      break;
    }
  }

#ifdef DEBUG_EUTACTICITY
  if (all_positive) {
    os << "EUTACTIC: All coefficients are positive - the form is eutactic\n";
    os << "EUTACTIC: Coefficients: ";
    for (int i = 0; i < coeffs.size(); i++) {
      os << coeffs(i) << " ";
    }
    os << "\n";
  } else {
    os << "EUTACTIC: Some coefficients are zero or negative - not eutactic\n";
  }
#endif

  return all_positive;
}

// Test if a quadratic form is perfect
// A quadratic form is perfect if it is determined by its shortest vectors
template <typename T>
bool IsPerfect(MyMatrix<T> const &eGram, MyMatrix<T> const& SHV_T, std::ostream& os) {
#ifdef DEBUG_EUTACTICITY
  os << "EUTACTIC: IsPerfect, beginning\n";
#endif
  int n = eGram.rows();
  if (eGram.cols() != n) {
    std::cerr << "Error: eGram must be square\n";
    return false;
  }
  if (SHV_T.cols() != n) {
    std::cerr << "Error: SHV_T must have n columns\n";
    return false;
  }

  int nbSHV = SHV_T.rows();
  if (nbSHV == 0) {
#ifdef DEBUG_EUTACTICITY
    os << "EUTACTIC: No shortest vectors provided\n";
#endif
    return false;
  }

  // Dimension of the space of symmetric matrices  
  int dimSymm = n * (n + 1) / 2;

  // Number of degrees of freedom for a positive definite quadratic form
  // This is dimSymm minus 1 (for scaling)
  int dimQuadForm = dimSymm - 1;

#ifdef DEBUG_EUTACTICITY
  os << "EUTACTIC: n=" << n << " nbSHV=" << nbSHV << " dimQuadForm=" << dimQuadForm << "\n";
#endif

  // Create the system of equations from the constraint that
  // all shortest vectors have the same length
  // For each pair of shortest vectors v_i, v_j: v_i^T G v_i = v_j^T G v_j

  // We need at least dimQuadForm independent constraints
  int numConstraints = nbSHV - 1; // We can form nbSHV-1 independent constraints

  if (numConstraints < dimQuadForm) {
#ifdef DEBUG_EUTACTICITY
    os << "EUTACTIC: Not enough shortest vectors for perfect form\n";
#endif
    return false;
  }

  MyMatrix<T> ConstraintMat(numConstraints, dimSymm);

  // Reference vector (first shortest vector)
  MyVector<T> v0 = GetMatrixRow(SHV_T, 0);
  MyMatrix<T> M0(n, n);
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < n; v++) {
      M0(u, v) = v0(u) * v0(v);
    }
  }
  MyVector<T> V0 = SymmetricMatrixToVector(M0);

  // Create constraint matrix
  for (int iSHV = 1; iSHV < nbSHV; iSHV++) {
    MyVector<T> vi = GetMatrixRow(SHV_T, iSHV);
    MyMatrix<T> Mi(n, n);
    for (int u = 0; u < n; u++) {
      for (int v = 0; v < n; v++) {
        Mi(u, v) = vi(u) * vi(v);
      }
    }
    MyVector<T> Vi = SymmetricMatrixToVector(Mi);

    // Constraint: V0 - Vi = 0 (same quadratic form value)
    MyVector<T> constraint = V0 - Vi;
    AssignMatrixRow(ConstraintMat, iSHV - 1, constraint);
  }

#ifdef DEBUG_EUTACTICITY
  os << "EUTACTIC: ConstraintMat constructed, size: " << ConstraintMat.rows() << "x" << ConstraintMat.cols() << "\n";
#endif

  // Check if we have exactly the right number of independent constraints
  int rank = RankMat(ConstraintMat);

#ifdef DEBUG_EUTACTICITY
  os << "EUTACTIC: constraint rank=" << rank << " dimQuadForm=" << dimQuadForm << "\n";
#endif

  bool is_perfect = (rank == dimQuadForm);

#ifdef DEBUG_EUTACTICITY
  if (is_perfect) {
    os << "EUTACTIC: The form is perfect\n";
  } else {
    os << "EUTACTIC: The form is NOT perfect\n";
  }
#endif

  return is_perfect;
}

// Test if a quadratic form is extreme (both perfect and eutactic)
template <typename T>
bool IsExtreme(MyMatrix<T> const &eGram, MyMatrix<T> const& SHV_T, std::ostream& os) {
#ifdef DEBUG_EUTACTICITY
  os << "EUTACTIC: IsExtreme, beginning\n";
#endif

  bool perfect = IsPerfect(eGram, SHV_T, os);
  bool eutactic = IsEutactic(eGram, SHV_T, os);

  bool extreme = perfect && eutactic;

#ifdef DEBUG_EUTACTICITY
  os << "EUTACTIC: perfect=" << perfect << " eutactic=" << eutactic << " extreme=" << extreme << "\n";
#endif

  return extreme;
}

// clang-format off
#endif  // SRC_PERFECT_EUTACTICITY_H_
// clang-format on
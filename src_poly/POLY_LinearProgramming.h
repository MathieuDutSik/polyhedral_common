// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_LINEARPROGRAMMING_H_
#define SRC_POLY_POLY_LINEARPROGRAMMING_H_

// clang-format off
#include "Basic_file.h"
#include "Basic_string.h"
#include "MAT_MatrixInt.h"
#include "POLY_Fundamental.h"
#include "POLY_LinearProgrammingFund.h"
#include "POLY_cddlib.h"
#include "basic_datafile.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
// clang-format on

/*
  Code for computing by using linear programming functionalities.

  The linear programming code is typically not in that place
  but in another. Typically we use CDD.
  But there are potential alternatives, for example GLPK or LRS
  though they do not yet offer the full functionalities as CDD.
 */

#ifdef DEBUG
#define DEBUG_LINEAR_PROGRAM
#define DEBUG_POLYTOPIZATION
#define DEBUG_SEARCH_POSITIVE_RELATION
#define DEBUG_GEOMETRICALLY_UNIQUE
#define DEBUG_FULL_RANK_FACET_SET
#define DEBUG_FIND_SINGLE_VERTEX
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_SEARCH_POSITIVE_RELATION
#define SANITY_CHECK_FULL_RANK_FACET_SET
#define SANITY_CHECK_FIND_SINGLE_VERTEX
#endif

#ifdef TIMINGS
#define TIMINGS_FULL_RANK_FACET_SET
#endif

template <typename T>
std::optional<MyMatrix<T>> AffinizeSubspace(MyMatrix<T> const &NSP) {
  int TheDim = NSP.rows();
  int nbCol = NSP.cols();
  if (TheDim == 0)
    return {};
  int FirstVertex = -1;
  for (int iRow = 0; iRow < TheDim; iRow++) {
    if (NSP(iRow, 0) != 0)
      FirstVertex = iRow;
  }
  if (FirstVertex == -1)
    return {};
  T ePivot = NSP(FirstVertex, 0);
  MyMatrix<T> NSPred(TheDim, nbCol);
  MyVector<T> eVertex(nbCol);
  for (int i = 0; i < nbCol; i++)
    NSPred(0, i) = NSP(FirstVertex, i) / ePivot;
  int pos = 0;
  for (int iRow = 0; iRow < TheDim; iRow++)
    if (iRow != FirstVertex) {
      pos++;
      for (int i = 0; i < nbCol; i++)
        NSPred(pos, i) = NSP(iRow, i) - NSP(iRow, 0) * NSPred(0, i);
    }
  return NSPred;
}

template <typename T>
std::optional<MyMatrix<T>>
ComputeStandardAffineSubspace(MyMatrix<T> const &ListEqua) {
  MyMatrix<T> NSP = NullspaceTrMat(ListEqua);
  return AffinizeSubspace(NSP);
}

template <typename T>
bool TestRealizabilityInequalitiesEqualities(MyMatrix<T> const &ListIneq,
                                             MyMatrix<T> const &ListEqua,
                                             MyVector<T> const &ToBeMinimized,
                                             std::ostream &os) {
  std::optional<MyMatrix<T>> EquivArr = ComputeStandardAffineSubspace(ListEqua);
  if (!EquivArr)
    return false;
  MyMatrix<T> NSPred = EquivArr.TheEquiv;
  MyMatrix<T> ListIneqRed = ListIneq * (NSPred.transpose());
  MyVector<T> ToBeMinimizedRed = NSPred * ToBeMinimized;
  LpSolution<T> eSol = CDD_LinearProgramming(ListIneqRed, ToBeMinimizedRed, os);
  if (eSol.PrimalDefined && eSol.DualDefined)
    return true;
  return false;
}

template <typename T> bool IsPolytopal(MyMatrix<T> const &EXT) {
  int nbRow = EXT.rows();
  for (int iRow = 0; iRow < nbRow; iRow++) {
    if (EXT(iRow, 0) != 1)
      return false;
  }
  return true;
}

template <typename T>
MyMatrix<T> Polytopization(MyMatrix<T> const &EXT, std::ostream &os) {
#ifdef DEBUG_POLYTOPIZATION
  os << "LP: Polytopization, beginning\n";
#endif
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  if (IsPolytopal(EXT)) {
#ifdef DEBUG_POLYTOPIZATION
    os << "LP: Early termination of the polytopization stuff\n";
#endif
    return EXT;
  }
  //
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  T eVal, prov;
  MyMatrix<T> nMat(nbRow, nbCol + 1);
  MyMatrix<T> eBasis(nbCol, nbCol);
  MyVector<T> eVect(nbCol + 1);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    eVal = T(-1);
    nMat(iRow, 0) = eVal;
    for (int iCol = 0; iCol < nbCol; iCol++)
      nMat(iRow, iCol + 1) = EXT(iRow, iCol);
  }
  for (int iCol = 0; iCol <= nbCol; iCol++) {
    T eSum(0);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      eSum += nMat(iRow, iCol);
    }
    eVect(iCol) = eSum;
  }
  //
#ifdef DEBUG_POLYTOPIZATION
  os << "LP: Polytopization, Before CDD_LinearProgramming\n";
#endif
  LpSolution<T> eSol = CDD_LinearProgramming(nMat, eVect, os);
#ifdef DEBUG_POLYTOPIZATION
  os << "LP: Polytopization, After CDD_LinearProgramming\n";
#endif
  if (!eSol.PrimalDefined) {
    std::cerr << "The optimization resulted in a result whose primal is not "
                 "defined\n";
    std::cerr << "This likely means that the set of inequalities on input does "
                 "not define\n";
    std::cerr << "A full dimensional polytope, i.e. that the inequalitity "
                 "together imply some equalities\n";
    throw TerminalException{1};
  }
  MyVector<T> SolDir = eSol.DirectSolution;
  ZeroAssignation(eBasis);
  for (int iCol = 0; iCol < nbCol; iCol++)
    eBasis(iCol, 0) = SolDir(iCol);
  int iColSelect = -1;
  for (int iCol = 0; iCol < nbCol; iCol++)
    if (iColSelect == -1) {
      if (eBasis(iCol, 0) != 0)
        iColSelect = iCol;
    }
  if (iColSelect == -1) {
    std::cerr << "Apparently, we did not find the column\n";
    std::cerr << "Please solve error in Polytopization\n";
    throw TerminalException{1};
  }
  int iRowWrite = 0;
  eVal = 1;
  for (int iCol = 0; iCol < nbCol; iCol++)
    if (iCol != iColSelect) {
      iRowWrite++;
      eBasis(iCol, iRowWrite) = eVal;
    }
  MyMatrix<T> EXTret = EXT * eBasis;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    T eMult = 1 / EXTret(iRow, 0);
    for (int iCol = 0; iCol < nbCol; iCol++)
      EXTret(iRow, iCol) = EXTret(iRow, iCol) * eMult;
  }
  return EXTret;
}

template <typename T> struct SolVertex {
  Face eInc;
  MyVector<T> eVert;
};

template <typename T> MyMatrix<T> SetIsobarycenter(MyMatrix<T> const &EXT) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  MyVector<T> eVect = MyVector<T>(nbCol);
  for (int iCol = 0; iCol < nbCol; iCol++) {
    T eSum(0);
    for (int iRow = 0; iRow < nbRow; iRow++)
      eSum += EXT(iRow, iCol);
    eSum = eSum / T(nbRow);
    eVect(iCol) = eSum;
  }
  MyMatrix<T> nMat(nbRow, nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    nMat(iRow, 0) = EXT(iRow, 0);
    for (int iCol = 1; iCol < nbCol; iCol++)
      nMat(iRow, iCol) = EXT(iRow, iCol) - eVect(iCol);
  }
  return nMat;
}

template <typename T>
Face Kernel_FindSingleVertex(MyMatrix<T> const &EXT, std::ostream &os) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  MyVector<T> eVect = MyVector<T>(nbCol);
  eVect(0) = 0;
#ifdef DEBUG_FIND_SINGLE_VERTEX
  size_t n_iter = 0;
#endif
#ifdef DEBUG_FIND_SINGLE_VERTEX_DISABLE
  os << "LP: Kernel_FindSingleVertex, EXT=\n";
  WriteMatrix(os, EXT);
#endif
  while (true) {
    for (int iCol = 1; iCol < nbCol; iCol++) {
      int a = random();
      int b = random();
      T eVal(a - b);
      eVect(iCol) = eVal;
    }
    LpSolution<T> eSol = CDD_LinearProgramming(EXT, eVect, os);
    MyVector<T> const& SolDir = eSol.DirectSolution;
    T const& OptimalValue = eSol.OptimalValue;
#ifdef DEBUG_FIND_SINGLE_VERTEX
    os << "LP: nbCol=" << nbCol << " |SolDir|=" << SolDir.size() << "\n";
    os << "LP: optimal_value=" << OptimalValue << "\n";
    os << "LP: SolDir=" << StringVector(SolDir) << "\n";
#endif
    Face Inc_eSum(nbRow);
    Face Inc_contrib(nbRow);
    T min_contrib(0);
    for (int iRow = 0; iRow < nbRow; iRow++) {
#ifdef SANITY_CHECK_FIND_SINGLE_VERTEX
      if (EXT(iRow,0) == 1) {
        std::cerr << "The first entry should be equal to 1\n";
        throw TerminalException{1};
      }
#endif
      T contrib(0);
      for (int iCol = 0; iCol < nbCol-1; iCol++) {
        contrib += EXT(iRow, iCol + 1) * SolDir(iCol);
      }
      T eSum = EXT(iRow, 0) + contrib;
      if (eSum == 0) {
        Inc_eSum[iRow] = 1;
      }
      if (contrib == 0) {
        Inc_contrib[iRow] = 1;
      }
      if (min_contrib > contrib) {
        min_contrib = contrib;
      }
#ifdef DEBUG_FIND_SINGLE_VERTEX
      os << "LP: iRow=" << iRow << " contrib=" << contrib << " eSum=" << eSum << "\n";
#endif
#ifdef SANITY_CHECK_FIND_SINGLE_VERTEX
      if (eSum < 0) {
        std::cerr << "The obtained vertex violates facet iRow=" << iRow << "\n";
        throw TerminalException{1};
      }
#endif
    }
#ifdef DEBUG_FIND_SINGLE_VERTEX
    n_iter += 1;
    os << "LP: Kernel_FindSingleVertex, min_contrib=" << min_contrib << "\n";
    os << "LP: Kernel_FindSingleVertex, n_iter=" << n_iter << "\n";
#endif
    if (min_contrib == 0 && OptimalValue == 0) {
      if (TestFacetInequality(EXT, Inc_contrib)) {
#ifdef DEBUG_FIND_SINGLE_VERTEX
        os << "LP: Kernel_FindSingleVertex, returning Inc_contrib at n_iter=" << n_iter << "\n";
#endif
        return Inc_contrib;
      }
    }
    if (TestFacetInequality(EXT, Inc_eSum)) {
#ifdef DEBUG_FIND_SINGLE_VERTEX
      os << "LP: Kernel_FindSingleVertex, returning Inc_eSum at n_iter=" << n_iter << "\n";
#endif
      return Inc_eSum;
    }
  }
}

template <typename T>
vectface Kernel_FindVertices(MyMatrix<T> const &EXT, size_t const &nb,
                             std::ostream &os) {
  vectface vf(EXT.rows());
  for (size_t i = 0; i < nb; i++) {
    Face f = Kernel_FindSingleVertex(EXT, os);
    vf.push_back(f);
  }
  return vf;
}

template <typename T>
vectface FindVertices(MyMatrix<T> const &EXT, int const &nb, std::ostream &os) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  MyMatrix<T> EXT_B = ColumnReduction(EXT);
  MyMatrix<T> EXT_C = Polytopization(EXT_B, os);
  MyMatrix<T> EXT_D = SetIsobarycenter(EXT_C);
  return Kernel_FindVertices(EXT_D, nb, os);
}

template <typename T>
Face FindOneInitialVertex(MyMatrix<T> const &TheEXT, std::ostream &os) {
  return FindVertices(TheEXT, 1, os)[0];
}

template <typename T>
MyMatrix<T> LP_GetExpressionForLP(MyMatrix<T> const &EXT) {
  MyVector<T> OneVert = EXT.row(0);
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  MyMatrix<T> ListDiff(nbRow, nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<T> eDiff = GetMatrixRow(EXT, iRow) - OneVert;
    AssignMatrixRow(ListDiff, iRow, eDiff);
  }
  MyMatrix<T> TheBasis = RowReduction(ListDiff);
  int TheDim = TheBasis.rows();
  MyMatrix<T> NewListCoord(nbRow, TheDim + 1);
  MyVector<T> eEntOne(1);
  eEntOne(0) = 1;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<T> eVect1 = ListDiff.row(iRow);
    std::optional<MyVector<T>> opt = SolutionMat(TheBasis, eVect1);
    const MyVector<T> &eVect2 = *opt;
    MyVector<T> eVect3 = Concatenation(eEntOne, eVect2);
    NewListCoord.row(iRow) = eVect3;
  }
  MyVector<T> eIso = Isobarycenter(NewListCoord);
  eIso(0) = 0;
  MyMatrix<T> RetMat(nbRow, TheDim + 1);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<T> V = GetMatrixRow(NewListCoord, iRow) - eIso;
    AssignMatrixRow(RetMat, iRow, V);
  }
  return RetMat;
}

template <typename T>
Face FindViolatedFace(MyMatrix<T> const &EXT, MyVector<T> const &eVect,
                      std::ostream &os) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  MyMatrix<T> EXT_ext(nbRow, nbCol + 1);
  MyVector<T> ToMinimize(nbCol + 1);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    EXT_ext(iRow, 0);
    for (int iCol = 0; iCol < nbCol; iCol++)
      EXT_ext(iRow, iCol + 1) = EXT(iRow, iCol);
  }
  ToMinimize(0) = 0;
  for (int iCol = 0; iCol < nbCol; iCol++)
    ToMinimize(iCol + 1) = eVect(iCol);
  LpSolution<T> eSol = CDD_LinearProgramming(EXT_ext, ToMinimize, os);
  Face eFace(nbRow);
  for (int i_row = 0; i_row < nbRow; i_row++) {
    T sum(0);
    for (int i_col = 0; i_col < nbCol; i_col++) {
      sum += EXT(i_row, i_col) * eSol.DirectSolution(i_col);
    }
    if (sum == 0)
      eFace[i_row] = 1;
  }
#ifdef DEBUG_LINEAR_PROGRAM
  MyVector<T> eFAC = FindFacetInequality(EXT, eFace);
  T scal = eVect.dot(eFAC);
  if (scal >= 0) {
    std::cerr << "Faied to find a correct solution\n";
    throw TerminalException{1};
  }
#endif
  return eFace;
}

// Given a configuration of vectors v_i, we are looking for the lambda_i >= 0
// such that sum_i lambda_i v_i = 0 with not all lambda_i zero.
// The existence of such a system of lambda is equivalent to the absence of a
// vector x such that x.v_i > 0 for all i.
//
// The following data type is the output of making such checks.
// eTestExist: If true then there exist a configuration of such beta. If false
// then no. InternalVector: The internal vector that satisfes x.v_i > 0 if
// existing. TheRelat: The relation between the vectors.
template <typename T> struct PosRelRes {
  bool eTestExist;
  bool has_relation;
  bool has_internal_vector;
  std::optional<MyVector<T>> InternalVector;
  std::optional<MyVector<T>> TheRelat;
};

template <typename T>
void CheckResult_PositiveRelationSimple(MyMatrix<T> const &ListVect,
                                        PosRelRes<T> const &pos) {
  int nbRow = ListVect.rows();
  int nbCol = ListVect.cols();
  // Checking the internal vector (if existing of course)
  if (pos.has_internal_vector) {
    if (!pos.eTestExist) {
      if (!pos.InternalVector.has_value()) {
        std::cerr << "There is an internal vector while eTestExist is true\n";
        throw TerminalException{1};
      }
      MyVector<T> const &V = *pos.InternalVector;
      for (int iRow = 0; iRow < nbRow; iRow++) {
        T eScal(0);
        for (int iCol = 0; iCol < nbCol; iCol++)
          eScal += V(iCol) * ListVect(iRow, iCol);
        if (eScal <= 0) {
          std::cerr << "Error in SearchPositiveRelationSimple_DualMethod 1\n";
          std::cerr << "We have eScal=" << eScal << "\n";
          throw TerminalException{1};
        }
      }
    } else {
      if (pos.InternalVector.has_value()) {
        std::cerr << "There is no internal vector while eTestExist is false\n";
        throw TerminalException{1};
      }
    }
  }
  // Checking the relation now.
  if (pos.has_relation) {
    if (!pos.eTestExist) {
      if (pos.TheRelat.has_value()) {
        std::cerr
            << "There is no relation for TheRelat while eTestExist is true\n";
        throw TerminalException{1};
      }
    } else {
      if (!pos.TheRelat.has_value()) {
        std::cerr
            << "There is a relation for TheRelat while eTestExist is false\n";
        throw TerminalException{1};
      }
      MyVector<T> const &V = *pos.TheRelat;
      T sumV(0);
      for (int iRow = 0; iRow < nbRow; iRow++) {
        if (V(iRow) < 0) {
          std::cerr << "iRow=" << iRow << "\n";
          std::cerr << "We have coeff=" << V(iRow) << " < 0\n";
          throw TerminalException{1};
        }
        sumV += V(iRow);
      }
      if (sumV == 0) {
        std::cerr << "All the coefficient are 0, which is not allowed\n";
        throw TerminalException{1};
      }
      for (int iCol = 0; iCol < nbCol; iCol++) {
        T eSum(0);
        for (int iRow = 0; iRow < nbRow; iRow++)
          eSum += V(iRow) * ListVect(iRow, iCol);
        if (eSum != 0) {
          std::cerr << "CheckResult_PositiveRelationSimple : error checking "
                       "TheRelat\n";
          std::cerr << "CheckResult_PositiveRelationSimple : iCol=" << iCol
                    << "\n";
          std::cerr << "CheckResult_PositiveRelationSimple : We have eSum="
                    << eSum << "\n";
          WriteMatrixFile("CheckResult_PositiveRelationSimple", ListVect);
          throw TerminalException{1};
        }
      }
    }
  }
}

// Given a configuration of vectors v_i, we are looking for the lambda_i >= 0
// such that sum_i lambda_i v_i = 0.
// The existence of such a system of lambda is equivalent to the absence of a
// vector x such that x.v_i > 0 for all i.
//
// I need to finish understanding what is the idea below.
template <typename T>
PosRelRes<T>
SearchPositiveRelationSimple_DualMethod(MyMatrix<T> const &ListVect,
                                        std::ostream &os) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int nbRow = ListVect.rows();
  int nbCol = ListVect.cols();
  MyMatrix<T> ListVectExt(nbRow, nbCol + 1);
  MyVector<T> eMinimize(nbCol + 1);
  eMinimize(0) = 0;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    ListVectExt(iRow, 0) = -1;
    for (int iCol = 0; iCol < nbCol; iCol++) {
      ListVectExt(iRow, iCol + 1) = ListVect(iRow, iCol);
      eMinimize(iCol + 1) += ListVect(iRow, iCol);
    }
  }
  LpSolution<T> eSol = CDD_LinearProgramming(ListVectExt, eMinimize, os);
  PosRelRes<T> eResult;
  eResult.has_relation = false;
  eResult.has_internal_vector = true;
  bool IsDone = false;
  if (eSol.PrimalDefined && eSol.DualDefined) {
    IsDone = true;
    eResult.eTestExist = false;
    eResult.InternalVector = eSol.DirectSolution;
  }
  if (!eSol.PrimalDefined && eSol.DualDefined) {
    IsDone = true;
    eResult.eTestExist = true;
#ifdef PRINT_LINEAR_PROGRAM_RESULT
    MyVector<T> const &V = eSol.DualSolution;
    // That seems actually less easy to obtain than we expected.
    // More work is needed here.
    T max_V = V(0);
    T min_V = V(0);
    for (int iRow = 1; iRow < nbRow; iRow++) {
      if (V(iRow) < min_V) {
        min_V = V(iRow);
      }
      if (V(iRow) > max_V) {
        max_V = V(iRow);
      }
    }
    double max_V_d = UniversalScalarConversion<double, T>(max_V);
    double min_V_d = UniversalScalarConversion<double, T>(min_V);
    std::cerr << "max_V_d=" << max_V_d << " min_V_d=" << min_V_d << "\n";
    for (int iCol = 0; iCol < nbCol; iCol++) {
      T eSum(0);
      for (int iRow = 0; iRow < nbRow; iRow++) {
        eSum += V(iRow) * ListVect(iRow, iCol);
      }
      std::cerr << "iCol=" << iCol << " eSum=" << eSum
                << " eMin(iRow)=" << eMinimize(iCol + 1) << "\n";
    }
#endif
  }
  if (!IsDone) {
    std::cerr << "Error. No value assigned\n";
    throw TerminalException{1};
  }
#ifdef DEBUG_SEARCH_POSITIVE_RELATION
  CheckResult_PositiveRelationSimple(ListVect, eResult);
#endif
  return eResult;
}

struct Constraint {
  std::vector<int> ListStrictlyPositive;
  std::vector<int> ListPositive;
  std::vector<std::vector<int>> ListSetStrictPositive;
};

// Given a configuration of vectors v_i, we are looking for the lambda_i >= 0
// such that sum_i lambda_i v_i = 0.
// We compute the kernel of the space of the equation space of the system.
// So, if we have a configuration Lambda_1, ...., Lambda_N
// So we write lambda = y_1 Lambda_1 + ... + y_N Lambda_N.
// We have the constraint on the y_i.
//
// The dimensionality is the following:
// If we have m vectors in dimension n, then the dimension of NSP is m-n.
// and the number of constraints is at least m.
template <typename T>
PosRelRes<T> SearchPositiveRelation(MyMatrix<T> const &ListVect,
                                    Constraint const &eConstraint,
                                    std::ostream &os) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  MyMatrix<T> NSP = NullspaceMat(ListVect);
  int nbVect = ListVect.rows();
  int nbRelation = NSP.rows();
  std::vector<MyVector<T>> ListInequalities;
  for (auto &eVect : eConstraint.ListStrictlyPositive) {
    MyVector<T> eV(nbRelation + 1);
    eV(0) = -1;
    for (int iRel = 0; iRel < nbRelation; iRel++)
      eV(iRel + 1) = NSP(iRel, eVect);
    ListInequalities.push_back(eV);
  }
  for (auto &eVect : eConstraint.ListPositive) {
    MyVector<T> eV(nbRelation + 1);
    eV(0) = 0;
    for (int iRel = 0; iRel < nbRelation; iRel++)
      eV(iRel + 1) = NSP(iRel, eVect);
    ListInequalities.push_back(eV);
  }
  for (auto &eSet : eConstraint.ListSetStrictPositive) {
    MyVector<T> eV = ZeroVector<T>(nbRelation + 1);
    eV(0) = -1;
    for (auto &eVect : eSet) {
      for (int iRel = 0; iRel < nbRelation; iRel++)
        eV(iRel + 1) += NSP(iRel, eVect);
    }
    ListInequalities.push_back(eV);
  }
  MyVector<T> ToBeMinimized(nbRelation + 1);
  ToBeMinimized(0) = 0;
  for (int iRel = 0; iRel < nbRelation; iRel++) {
    T eSum(0);
    for (int iVect = 0; iVect < nbVect; iVect++)
      eSum += NSP(iRel, iVect);
    ToBeMinimized(iRel + 1) = eSum;
  }
#ifdef DEBUG_SEARCH_POSITIVE_RELATION
  os << "LP: |ListInequalities|=" << ListInequalities.size() << "\n";
#endif
  MyMatrix<T> MatInequalities = MatrixFromVectorFamily(ListInequalities);
  LpSolution<T> eSol =
      CDD_LinearProgramming(MatInequalities, ToBeMinimized, os);
  PosRelRes<T> eResult;
  eResult.has_relation = true;
  eResult.has_internal_vector = false;
  if (eSol.PrimalDefined && eSol.DualDefined) {
    MyVector<T> DirSol = eSol.DirectSolution;
    eResult.eTestExist = true;
    MyVector<T> eVectRel(nbRelation);
    for (int iRel = 0; iRel < nbRelation; iRel++)
      eVectRel(iRel) = DirSol(iRel);
    MyVector<T> TheRelat = ProductVectorMatrix(eVectRel, NSP);
    for (int iVect = 0; iVect < nbVect; iVect++) {
      if (TheRelat(iVect) < 0) {
        std::cerr << "We have a clear and present error\n";
        throw TerminalException{1};
      }
    }
    eResult.TheRelat = TheRelat;
  } else {
    eResult.eTestExist = false;
    // Finding the x vector if there is no relation seems a little difficult.
    // So, we pass it out for now.
  }
  return eResult;
}

template <typename T>
PosRelRes<T> SearchPositiveRelationSimple_Direct(MyMatrix<T> const &ListVect,
                                                 std::ostream &os) {
  int nbVect = ListVect.rows();
  std::vector<int> ListStrictlyPositive;
  std::vector<int> ListPositive(nbVect);
  for (int iVect = 0; iVect < nbVect; iVect++)
    ListPositive[iVect] = iVect;
  std::vector<std::vector<int>> ListSetStrictPositive = {ListPositive};
  Constraint eConstraint{ListStrictlyPositive, ListPositive,
                         ListSetStrictPositive};
  return SearchPositiveRelation(ListVect, eConstraint, os);
}

template <typename T>
bool TestExistPositiveRelation(MyMatrix<T> const &ListVect, std::ostream &os) {
  int nbRow = ListVect.rows();
  int nbCol = ListVect.cols();
  int dim_direct = nbRow - nbCol;
  int dim_dual = nbCol;
#ifdef SANITY_CHECK_SEARCH_POSITIVE_RELATION
  // Relatively expensive checks to do
  PosRelRes<T> sol1 = SearchPositiveRelationSimple_Direct(ListVect, os);
  PosRelRes<T> sol2 = SearchPositiveRelationSimple_DualMethod(ListVect, os);
  if (sol1.eTestExist != sol2.eTestExist) {
    std::cerr << "We get different result for Direct and DualMethod\n";
    std::cerr << "Clearly not what we expected\n";
  }
#endif
  auto get_solution = [&]() -> PosRelRes<T> {
    // We take the dimensionality as the driving factor for the complexity of
    // the computation.
    if (dim_direct < dim_dual) {
      return SearchPositiveRelationSimple_Direct(ListVect, os);
    }
    return SearchPositiveRelationSimple_DualMethod(ListVect, os);
  };
  PosRelRes<T> the_sol = get_solution();
#ifdef DEBUG_SEARCH_POSITIVE_RELATION
  CheckResult_PositiveRelationSimple(ListVect, the_sol);
#endif
  return the_sol.eTestExist;
}

template <typename T>
std::optional<MyVector<T>>
SolutionMatNonnegative_Version1(MyMatrix<T> const &ListVect,
                                MyVector<T> const &eVect, std::ostream &os) {
  int nbVect = ListVect.rows();
  int nbCol = ListVect.cols();
  MyMatrix<T> InputListVect(nbVect + 1, nbCol);
  for (int iVect = 0; iVect < nbVect; iVect++)
    InputListVect.row(iVect) = ListVect.row(iVect);
  for (int iCol = 0; iCol < nbCol; iCol++)
    InputListVect(nbVect, iCol) = -eVect(iCol);
  //
  std::vector<int> ListStrictlyPositive;
  std::vector<int> ListPositive(nbVect + 1);
  for (int iVect = 0; iVect <= nbVect; iVect++)
    ListPositive[iVect] = iVect;
  std::vector<int> Part1(nbVect);
  for (int iVect = 0; iVect < nbVect; iVect++)
    Part1[iVect] = iVect;
  std::vector<int> Part2 = {nbVect};
  std::vector<std::vector<int>> ListSetStrictPositive = {Part1, Part2};
  Constraint eConstraint{ListStrictlyPositive, ListPositive,
                         ListSetStrictPositive};
  //
  PosRelRes<T> PRR = SearchPositiveRelation(InputListVect, eConstraint, os);
  if (!PRR.eTestExist)
    return {};
  MyVector<T> TheSol(nbVect);
  MyVector<T> const &V = *PRR.TheRelat;
  for (int iVect = 0; iVect < nbVect; iVect++)
    TheSol(iVect) = V(iVect) / V(nbVect);
  return TheSol;
}

template <typename T>
std::optional<MyVector<T>>
SolutionMatNonnegative_LP(MyMatrix<T> const &ListVect, MyVector<T> const &eVect,
                          std::ostream &os) {
  int nbVect = ListVect.rows();
  int nbCol = ListVect.cols();
  MyMatrix<T> ListIneq(nbVect, nbCol + 1);
  MyVector<T> eIneq(nbCol + 1);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    ListIneq(iVect, 0) = 0;
    for (int iCol = 0; iCol < nbCol; iCol++)
      ListIneq(iVect, 1 + iCol) = ListVect(iVect, iCol);
  }
  eIneq(0) = 0;
  for (int iCol = 0; iCol < nbCol; iCol++)
    eIneq(1 + iCol) = eVect(iCol);
  //
  LpSolution<T> eSol = CDD_LinearProgramming(ListIneq, eIneq, os);
  if (!eSol.DualDefined) {
    return {};
  }
  MyVector<T> TheRet = -eSol.DualSolution;
  return TheRet;
}

template <typename T> struct SolutionMatNonnegativeComplete {
  std::optional<MyVector<T>> ExtremeRay;
  std::optional<MyVector<T>> SolNonnegative;
};

template <typename T>
SolutionMatNonnegativeComplete<T>
GetSolutionMatNonnegativeComplete(MyMatrix<T> const &ListVect,
                                  MyVector<T> const &eVect, std::ostream &os) {
  int nbVect = ListVect.rows();
  int nbCol = ListVect.cols();
  MyMatrix<T> ListIneq(nbVect, nbCol + 1);
  MyVector<T> eIneq(nbCol + 1);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    ListIneq(iVect, 0) = 0;
    for (int iCol = 0; iCol < nbCol; iCol++)
      ListIneq(iVect, 1 + iCol) = ListVect(iVect, iCol);
  }
  eIneq(0) = 0;
  for (int iCol = 0; iCol < nbCol; iCol++)
    eIneq(1 + iCol) = eVect(iCol);
  //
  LpSolution<T> eSol = CDD_LinearProgramming(ListIneq, eIneq, os);
  auto GetSolNonnegative = [&]() -> std::optional<MyVector<T>> {
    if (!eSol.DualDefined) {
      return {};
    }
    MyVector<T> TheRet = -eSol.DualSolution;
    return TheRet;
  };
  std::optional<MyVector<T>> SolNonnegative = GetSolNonnegative();
  auto GetExtremeRay = [&]() -> std::optional<MyVector<T>> {
    if (eSol.PrimalDefined) {
      MyVector<T> V = eSol.DirectSolution;
      MyVector<T> LScal = ListVect * V;
      for (int iVect = 0; iVect < nbVect; iVect++) {
        if (LScal(iVect) < 0) {
          std::cerr << "LScal=" << StringVectorGAP(LScal) << "\n";
          std::cerr << "Negative value at iVect=" << iVect
                    << " LScal(iVect)=" << LScal(iVect) << "\n";
          throw TerminalException{1};
        }
      }
      T eScal = V.dot(eVect);
      if (eScal > 0) {
        std::cerr << "eScal=" << eScal << "\n";
        std::cerr << "The direction is not a counter example\n";
        throw TerminalException{1};
      }
      if (eScal == 0 && !SolNonnegative) {
        std::cerr << "eScal = 0 and no Non-negative solution were found\n";
        std::cerr << "Possibly a bug here\n";
        throw TerminalException{1};
      }
      return V;
    }
    return {};
  };
  std::optional<MyVector<T>> ExtremeRay = GetExtremeRay();
  return {ExtremeRay, SolNonnegative};
}

template <typename T>
std::optional<MyVector<T>>
SolutionMatNonnegative_Check(MyMatrix<T> const &ListVect,
                             MyVector<T> const &eVect, std::ostream &os) {
  std::optional<MyVector<T>> opt1 =
      SolutionMatNonnegative_Version1(ListVect, eVect, os);
  std::optional<MyVector<T>> opt2 =
      SolutionMatNonnegative_LP(ListVect, eVect, os);
  if (opt1 && !opt2) {
    std::cerr << "opt1 is defined but not opt2, incoherent\n";
    throw TerminalException{1};
  }
  if (!opt1 && opt2) {
    std::cerr << "opt1 is not defined but opt2 is, incoherent\n";
    throw TerminalException{1};
  }
  if (!opt1 && !opt2)
    return {};
  MyVector<T> const &V = *opt1;
  int nbRow = V.size();
  if (nbRow != ListVect.rows()) {
    std::cerr << "Dimension incoherency\n";
    throw TerminalException{1};
  }
  for (int iRow = 0; iRow < nbRow; iRow++) {
    if (V(iRow) < 0) {
      std::cerr << "V should be non-negative\n";
      throw TerminalException{1};
    }
  }
  MyVector<T> prod = ListVect.transpose() * V;
  if (prod != eVect) {
    std::cerr << "prod does not matcheVect\n";
    throw TerminalException{1};
  }
  return V;
}

template <typename T>
std::optional<MyVector<T>>
SolutionMatNonnegative_FailSafe(MyMatrix<T> const &ListVect,
                                MyVector<T> const &eVect, std::ostream &os) {
  auto local_terminate = [&]() -> void {
    std::string FILE_FAC = "DEBUG_Nonnegative_FAC";
    WriteMatrixFile(FILE_FAC, ListVect);
    std::string FILE_INEQ = "DEBUG_Nonnegative_INEQ";
    WriteVectorFile(FILE_INEQ, eVect);
    throw TerminalException{1};
  };
  std::optional<MyVector<T>> opt_LP =
      SolutionMatNonnegative_LP(ListVect, eVect, os);
  if (opt_LP) {
    MyVector<T> const &V = *opt_LP;
    int nbRow = V.size();
    if (nbRow != ListVect.rows()) {
      std::cerr << "Dimension incoherency\n";
      local_terminate();
    }
    bool HasError = false;
    for (int iRow = 0; iRow < nbRow; iRow++) {
      T val = V(iRow);
      if (val < 0) {
        double val_d = UniversalScalarConversion<double, T>(val);
        std::cerr << "iRow=" << iRow << " val=" << val << " val_d=" << val_d
                  << "\n";
        HasError = true;
      }
    }
    MyVector<T> prod = ListVect.transpose() * V;
    if (prod != eVect) {
      std::cerr << "prod does not matcheVect\n";
      std::cerr << "|prod|=" << prod.rows() << " / " << prod.cols() << "\n";
      std::cerr << "|eVect|=" << eVect.rows() << " / " << eVect.cols() << "\n";
      for (int iRow = 0; iRow < prod.rows(); iRow++) {
        double prod_d = UniversalScalarConversion<double, T>(prod(iRow, 0));
        double eVect_d = UniversalScalarConversion<double, T>(eVect(iRow, 0));
        std::cerr << "iRow=" << iRow << " prod_d=" << prod_d
                  << " eVect_d=" << eVect_d << "\n";
      }
      HasError = true;
    }
    if (HasError) {
      local_terminate();
    }
    return V;
  }
  std::optional<MyVector<T>> opt_V1 =
      SolutionMatNonnegative_Version1(ListVect, eVect, os);
  if (opt_V1) {
    std::cerr << "opt_V1 is defined but not opt_LP, incoherent\n";
    local_terminate();
  }
  return {};
}

template <typename T>
std::optional<MyVector<T>> SolutionMatNonnegative(MyMatrix<T> const &ListVect,
                                                  MyVector<T> const &eVect,
                                                  std::ostream &os) {
  //  return SolutionMatNonnegative_Check(ListVect, eVect);
  return SolutionMatNonnegative_FailSafe(ListVect, eVect, os);
}

template <typename T>
Face ComputeSkeletonClarkson(MyMatrix<T> const &FACinp, std::ostream &os) {
  MyMatrix<T> FAC = ColumnReduction(FACinp);
  int n_fac = FAC.rows();
  int dim = FAC.cols();
  Face f_adj(n_fac * n_fac);
  for (int i_fac = 0; i_fac < n_fac; i_fac++) {
    MyMatrix<T> Equa(dim, 1);
    for (int i = 0; i < dim; i++)
      Equa(i, 0) = FAC(i_fac, i);
    MyMatrix<T> NSP = NullspaceMat(Equa);
    std::vector<int> ListIdx;
    for (int j_fac = 0; j_fac < i_fac; j_fac++) {
      if (f_adj[i_fac + j_fac * n_fac] == 1) {
        ListIdx.push_back(j_fac);
      }
    }
    for (int j_fac = i_fac + 1; j_fac < n_fac; j_fac++)
      ListIdx.push_back(j_fac);
    int n_ineq = ListIdx.size();
    MyMatrix<T> FAC_local(n_ineq, dim);
    for (int i_ineq = 0; i_ineq < n_ineq; i_ineq++) {
      int idx = ListIdx[i_ineq];
      MyVector<T> eFAC = GetMatrixRow(FAC, idx);
      MyVector<T> eFACred = NSP * eFAC;
      FAC_local(i_ineq, 0) = 0;
      for (int i = 0; i < dim - 1; i++)
        FAC_local(i_ineq, i + 1) = eFACred(i);
    }
    std::vector<int> ListIrred =
        cdd::RedundancyReductionClarkson(FAC_local, os);
    for (auto &eIrred : ListIrred) {
      int pos = ListIdx[eIrred];
      f_adj[pos + i_fac * n_fac] = 1;
    }
  }
  return f_adj;
}

template <typename T>
std::pair<MyMatrix<T>, Face>
KernelLinearDeterminedByInequalitiesAndIndices_DualMeth(MyMatrix<T> const &FAC,
                                                        std::ostream &os) {
  int nbCol = FAC.cols();
  int nbRow = FAC.rows();
  PosRelRes<T> eRes = SearchPositiveRelationSimple_Direct(FAC, os);
  //  std::cerr << "eRes.eTestExist=" << eRes.eTestExist << "\n";
  if (!eRes.eTestExist) {
    MyMatrix<T> Spa = IdentityMat<T>(nbCol);
    Face f(nbRow);
    return {std::move(Spa), std::move(f)};
  } else {
    std::vector<int> ListIdx;
    MyVector<T> const &TheRelat = *eRes.TheRelat;
    for (int iRow = 0; iRow < nbRow; iRow++)
      if (TheRelat(iRow) > 0)
        ListIdx.push_back(iRow);
    MyMatrix<T> FACred = SelectRow(FAC, ListIdx);
    SelectionRowCol<T> eSelect = TMat_SelectRowCol(FACred);
    if (eSelect.NSP.rows() == 0) {
      MyMatrix<T> Spa = MyMatrix<T>(0, nbCol);
      Face f(nbRow);
      for (int u = 0; u < nbRow; u++)
        f[u] = 1;
      return {std::move(Spa), std::move(f)};
    }
    MyMatrix<T> FACproj = FAC * eSelect.NSP.transpose();
    int FACproj_col = FACproj.cols();
    auto IsZeroLine = [&](int const &iLine) -> bool {
      for (int u = 0; u < FACproj_col; u++)
        if (FACproj(iLine, u) != 0)
          return false;
      return true;
    };
    std::vector<int> ListIdxZ, ListIdxNZ;
    for (int u = 0; u < FACproj.rows(); u++) {
      if (IsZeroLine(u)) {
        ListIdxZ.push_back(u);
      } else {
        ListIdxNZ.push_back(u);
      }
    }
    MyMatrix<T> FACprojCor = SelectRow(FACproj, ListIdxNZ);
    auto f_recursive = [&]() -> std::pair<MyMatrix<T>, Face> {
      if (FACprojCor.rows() == 0) {
        MyMatrix<T> TheSpann = IdentityMat<T>(eSelect.NSP.rows());
        Face f(0);
        return {std::move(TheSpann), f};
      } else {
        return KernelLinearDeterminedByInequalitiesAndIndices_DualMeth(
            FACprojCor, os);
      }
    };
    std::pair<MyMatrix<T>, Face> pair = f_recursive();
    if (pair.first.rows() == 0) {
      MyMatrix<T> Spann = MyMatrix<T>(0, nbCol);
      Face f(nbRow);
      for (int u = 0; u < nbRow; u++) {
        f[u] = 1;
      }
      return {std::move(Spann), std::move(f)};
    } else {
      MyMatrix<T> TheSpann = pair.first * eSelect.NSP;
      Face f(nbRow);
      for (auto &eVal : ListIdxZ)
        f[eVal] = 1;
      int len = ListIdxNZ.size();
      for (int u = 0; u < len; u++) {
        f[ListIdxNZ[u]] = pair.second[u];
      }
      return {std::move(TheSpann), std::move(f)};
    }
  }
}

template <typename T>
std::pair<MyMatrix<T>, Face>
KernelLinearDeterminedByInequalitiesAndIndices(MyMatrix<T> const &FAC,
                                               std::ostream &os) {
  int n_row = FAC.rows();
  int n_col = FAC.cols();
  int dim_direct = n_col;
  int dim_dual = n_row - n_col;
#ifdef TIMINGS_LINEAR_PROGRAM
  HumanTime time;
  (void)KernelLinearDeterminedByInequalitiesAndIndices_DualMeth(FAC, os);
  os << "|LP: KernelLinearDeterminedByInequalitiesAndIndices_DualMeth|=" << time
     << "\n";
  size_t maxiter1 = 0; // An abstracting leakage for sure
  (void)cdd::KernelLinearDeterminedByInequalitiesAndIndices_DirectLP(FAC, maxiter1, os);
  os << "|LP: KernelLinearDeterminedByInequalitiesAndIndices_DirectLP|=" << time
     << "\n";
  (void)cdd::KernelLinearDeterminedByInequalitiesAndIndices_LPandNullspace(FAC,
                                                                           os);
  os << "|LP: KernelLinearDeterminedByInequalitiesAndIndices_LPandNullspace|="
     << time << "\n";
#endif
#ifdef DEBUG_LINEAR_PROGRAM
  auto test_equa = [](std::pair<MyMatrix<T>, Face> const &res1,
                      std::pair<MyMatrix<T>, Face> const &res2) -> bool {
    if (res1.second != res2.second) {
      return false;
    }
    return TestEqualitySpannedSpaces(res1.first, res2.first);
  };
  auto f_print = [&os](std::string const &origin,
                       std::pair<MyMatrix<T>, Face> const &res) -> void {
    os << "LP: f_print origin=" << origin << "\n";
    os << "LP: Face=" << StringFace(res.second) << "\n";
    os << "LP: Space=\n";
    WriteMatrix(os, res.first);
  };
  os << "LP: Before KernelLinearDeterminedByInequalitiesAndIndices_DualMeth\n";
  std::pair<MyMatrix<T>, Face> res_dual =
      KernelLinearDeterminedByInequalitiesAndIndices_DualMeth(FAC, os);
  os << "LP: |res_dual|=" << res_dual.first.rows() << " / "
     << res_dual.first.cols() << "\n";
  os << "LP: res_dual.second=" << StringFace(res_dual.second) << "\n";
  os << "Before KernelLinearDeterminedByInequalitiesAndIndices_DirectLP\n";
  os << "res_dual.first=\n";
  WriteMatrix(os, res_dual.first);
  size_t maxiter2 = 0; // An abstracting leakage for sure
  std::pair<MyMatrix<T>, Face> res_dir_lp =
    cdd::KernelLinearDeterminedByInequalitiesAndIndices_DirectLP(FAC, maxiter2, os);
  os << "|res_dir_lp|=" << res_dir_lp.first.rows() << " / "
     << res_dir_lp.first.cols() << "\n";
  os << "LP: res_dir_lp.second=" << StringFace(res_dir_lp.second) << "\n";
  os << "res_dir_lp.first=\n";
  WriteMatrix(os, res_dir_lp.first);
  os << "Before "
        "KernelLinearDeterminedByInequalitiesAndIndices_LPandNullspace\n";
  std::pair<MyMatrix<T>, Face> res_dir_lp_nsp =
      cdd::KernelLinearDeterminedByInequalitiesAndIndices_LPandNullspace(FAC,
                                                                         os);
  os << "|res_dir_lp_nsp|=" << res_dir_lp_nsp.first.rows() << " / "
     << res_dir_lp_nsp.first.cols() << "\n";
  os << "LP: res_dir_lp_nsp.second=" << StringFace(res_dir_lp_nsp.second)
     << "\n";
  os << "res_dir_lp_nsp.first=\n";
  WriteMatrix(os, res_dir_lp_nsp.first);
  auto f_print_all = [&]() -> void {
    os << "LP: f_print_all. FAC=\n";
    WriteMatrix(os, FAC);
    f_print("res_dual", res_dual);
    f_print("res_dir_lp", res_dir_lp);
    f_print("res_dir_lp_nsp", res_dir_lp_nsp);
  };
  os << "After KernelLinearDeterminedByInequalitiesAndIndices_LPandNullspace\n";
  if (!test_equa(res_dual, res_dir_lp)) {
    std::cerr << "res_dual and res_dir_lp are not equal\n";
    f_print_all();
    throw TerminalException{1};
  }
  os << "After test_equa1\n";
  if (!test_equa(res_dual, res_dir_lp_nsp)) {
    std::cerr << "res_dual and res_dir_lp_nsp are not equal\n";
    f_print_all();
    throw TerminalException{1};
  }
  os << "After test_equa2\n";
  auto test_correct = [&](std::pair<MyMatrix<T>, Face> const &res) -> void {
    for (int i_row = 0; i_row < n_row; i_row++) {
      if (res.second[i_row] == 1) {
        MyVector<T> eRow = GetMatrixRow(FAC, i_row);
        MyVector<T> Vscal = res.first * eRow;
        if (!IsZeroVector(Vscal)) {
          std::cerr << "The vector Vscal is not zero\n";
          throw TerminalException{1};
        }
      }
    }
  };
  test_correct(res_dual);
  test_correct(res_dir_lp);
  test_correct(res_dir_lp_nsp);
#endif
  if (dim_dual < dim_direct) {
    return KernelLinearDeterminedByInequalitiesAndIndices_DualMeth(FAC, os);
  } else {
    return cdd::KernelLinearDeterminedByInequalitiesAndIndices_LPandNullspace(
        FAC, os);
  }
}

/*
  Return the corresponding subspace and the indices of the inequalities that are
  matching with equality.
 */
template <typename T>
std::pair<MyMatrix<T>, Face>
LinearDeterminedByInequalitiesAndIndices(MyMatrix<T> const &FAC,
                                         std::ostream &os) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int n_row = FAC.rows();
  int n_col = FAC.cols();
  std::vector<int> ListIdxZ;
  std::unordered_map<MyVector<T>, std::vector<int>> map;
  for (int u = 0; u < n_row; u++) {
    MyVector<T> eRow = GetMatrixRow(FAC, u);
    if (IsZeroVector(eRow)) {
      ListIdxZ.push_back(u);
    } else {
      MyVector<T> eRowCan = CanonicalizeVector(eRow);
      map[eRowCan].push_back(u);
    }
  }
  int siz = map.size();
#ifdef DEBUG_LINEAR_PROGRAM
  os << "LP: n_row=" << n_row << " n_col=" << n_col << " siz=" << siz
     << " |ListIdxZ|=" << ListIdxZ.size() << "\n";
#endif
  MyMatrix<T> FACred(siz, n_col);
  std::vector<std::vector<int>> ll_idx;
  int pos = 0;
  for (auto &kv : map) {
    for (int u = 0; u < n_col; u++) {
      FACred(pos, u) = kv.first(u);
    }
    ll_idx.push_back(kv.second);
    pos++;
  }
  std::pair<MyMatrix<T>, Face> pair =
      KernelLinearDeterminedByInequalitiesAndIndices(FACred, os);
#ifdef DEBUG_LINEAR_PROGRAM
  os << "LP: LinearDeterminedByInequalitiesAndIndices pair.second="
     << StringFace(pair.second) << "\n";
#endif
  Face f(n_row);
  for (auto &eVal : ListIdxZ) {
    f[eVal] = 1;
  }
  for (int pos = 0; pos < siz; pos++) {
#ifdef DEBUG_LINEAR_PROGRAM
    os << "LP: LinearDeterminedByInequalitiesAndIndices pos=" << pos
       << " |ll_idx|=" << ll_idx[pos].size() << "\n";
#endif
    if (pair.second[pos] == 1) {
      for (auto &&idx : ll_idx[pos]) {
        f[idx] = 1;
      }
    }
  }
#ifdef DEBUG_LINEAR_PROGRAM
  os << "LP: LinearDeterminedByInequalitiesAndIndices f=" << StringFace(f)
     << "\n";
#endif
  return {std::move(pair.first), std::move(f)};
}

template <typename T>
MyMatrix<T> KernelLinearDeterminedByInequalities(MyMatrix<T> const &FAC,
                                                 std::ostream &os) {
  int nbCol = FAC.cols();
  int nbRow = FAC.rows();
  PosRelRes<T> eRes = SearchPositiveRelationSimple_Direct(FAC, os);
  //  std::cerr << "eRes.eTestExist=" << eRes.eTestExist << "\n";
  if (!eRes.eTestExist) {
    return IdentityMat<T>(nbCol);
  } else {
    std::vector<int> ListIdx;
    MyVector<T> const &TheRelat = *eRes.TheRelat;
    for (int iRow = 0; iRow < nbRow; iRow++)
      if (TheRelat(iRow) > 0)
        ListIdx.push_back(iRow);
    MyMatrix<T> FACred = SelectRow(FAC, ListIdx);
    SelectionRowCol<T> eSelect = TMat_SelectRowCol(FACred);
    if (eSelect.NSP.rows() == 0)
      return MyMatrix<T>(0, nbCol);
    MyMatrix<T> FACproj = FAC * eSelect.NSP.transpose();
    MyMatrix<T> FACprojCor = SelectNonZeroRows(FACproj);
    MyMatrix<T> TheSpann;
    if (FACprojCor.rows() == 0)
      TheSpann = IdentityMat<T>(eSelect.NSP.rows());
    else
      TheSpann = KernelLinearDeterminedByInequalities(FACprojCor, os);
    if (TheSpann.rows() == 0)
      return MyMatrix<T>(0, nbCol);
    else
      return TheSpann * eSelect.NSP;
  }
}

template <typename T>
bool IsFullDimensional_V1(MyMatrix<T> const &FAC, std::ostream &os) {
  return !TestExistPositiveRelation(FAC, os);
}

template <typename T>
MyMatrix<T> LinearDeterminedByInequalities(MyMatrix<T> const &FAC,
                                           std::ostream &os) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  return KernelLinearDeterminedByInequalities(SelectNonZeroRows(FAC), os);
}

template <typename T>
bool IsFullDimensional(MyMatrix<T> const &FAC, std::ostream &os) {
  int nbRow = FAC.rows();
  int nbCol = FAC.cols();
  MyMatrix<T> FACexp(nbRow, 1 + nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    FACexp(iRow, 0) = T(-1);
    for (int iCol = 0; iCol < nbCol; iCol++)
      FACexp(iRow, 1 + iCol) = FAC(iRow, iCol);
  }
  MyVector<T> SumIneq(1 + nbCol);
  for (int iCol = 0; iCol <= nbCol; iCol++) {
    T sum(0);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      sum += FACexp(iRow, iCol);
    }
    SumIneq(iCol) = sum;
  }
  LpSolution<T> eSol = CDD_LinearProgramming(FACexp, SumIneq, os);
  if (eSol.PrimalDefined && eSol.DualDefined) {
    return true;
  }
  return false;
}

/* Finds an interior point in the cone determined by the inequalities.
   No group used here nor any equalities.
 */
template <typename T>
MyVector<T> GetSpaceInteriorPoint_Basic(MyMatrix<T> const &FAC,
                                        std::ostream &os) {
  int n_rows = FAC.rows();
  int n_cols = FAC.cols();
  MyMatrix<T> ListInequalities = ZeroMatrix<T>(n_rows, n_cols + 1);
  MyVector<T> ToBeMinimized = ZeroVector<T>(n_cols + 1);
  for (int i_row = 0; i_row < n_rows; i_row++) {
    ListInequalities(i_row, 0) = -1;
    for (int i_col = 0; i_col < n_cols; i_col++)
      ListInequalities(i_row, i_col + 1) = FAC(i_row, i_col);
    for (int i_col = 0; i_col <= n_cols; i_col++)
      ToBeMinimized(i_col) += ListInequalities(i_row, i_col);
  }
  LpSolution<T> eSol =
      CDD_LinearProgramming(ListInequalities, ToBeMinimized, os);
  if (!eSol.PrimalDefined || !eSol.DualDefined) {
    std::cerr << "Failed to find an interior point by linear programming\n";
    std::cerr << "Maybe the cone is actually not full dimensional\n";
    throw TerminalException{1};
  }
  MyVector<T> eVect = eSol.DirectSolution;
#ifdef DEBUG_LINEAR_PROGRAM
  MyVector<T> ListScal = FAC * eVect;
  for (int i_row = 0; i_row < n_rows; i_row++) {
    T eScal = ListScal(i_row);
    if (eScal <= 0) {
      std::cerr << "We have eScal=" << eScal << " which should be positive\n";
      throw TerminalException{1};
    }
  }
#endif
  return eVect;
}

// We check whether the linear program is non-degenerate or not.
// By this, we mean, whether there is a unique solution to the
// linear program or not.
//
// The argument should be to take a solution of the linear program.
// eSol should a valid solution to the linear program.
//
// If the code returns true, then the vertex is indeed non-degenerate.
// If not, then the vertex may or may not be non-degenerate.
template <typename T>
bool TestCriterionNonDegenerate(
    MyMatrix<T> const &ListIneq,
    [[maybe_unused]] MyVector<T> const &ToBeMinimized,
    LpSolution<T> const &eSol, [[maybe_unused]] std::ostream &os) {
  int n_row = ListIneq.rows();
  int n_col = ListIneq.cols();
  MyVector<T> eSolDual = -eSol.DualSolution;
#ifdef DEBUG_LINEAR_PROGRAM
  os << "LP: OptimalValue=" << eSol.OptimalValue << "\n";
  os << "LP: eSolDual=" << StringVectorGAP(eSolDual) << "\n";
  MyVector<T> TheSum = ToBeMinimized;
  TheSum(0) -= eSol.OptimalValue;
  os << "LP: 1: TheSum=" << StringVectorGAP(TheSum) << "\n";
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyVector<T> eRow = GetMatrixRow(ListIneq, i_row);
    TheSum -= eSolDual(i_row) * eRow;
  }
  os << "LP: 2: TheSum=" << StringVectorGAP(TheSum) << "\n";
  if (!IsZeroVector(TheSum)) {
    std::cerr << "The vector TheSum should be zero\n";
    throw TerminalException{1};
  }
  std::vector<MyVector<T>> vect_incd1;
  for (int i_row = 0; i_row < n_row; i_row++) {
    if (eSolDual(i_row) > 0) {
      MyVector<T> eRow = GetMatrixRow(ListIneq, i_row);
      vect_incd1.push_back(eRow);
    }
  }
  os << "LP: |vect_incd1|=" << vect_incd1.size() << "\n";
  MyMatrix<T> MatIncd1 = MatrixFromVectorFamily(vect_incd1);
  int rnk_incd = RankMat(MatIncd1);
  os << "LP: rnk_incd=" << rnk_incd << "\n";
  std::vector<MyVector<T>> vect_incd2;
  MyVector<T> DirectSolutionExt = GetDirectSolutionExt(eSol);
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyVector<T> eRow = GetMatrixRow(ListIneq, i_row);
    T scal = eRow.dot(DirectSolutionExt);
    if (scal < 0) {
      std::cerr << "The scalar product is negative, which is not allowed\n";
      throw TerminalException{1};
    }
    if (scal == 0) {
      vect_incd2.push_back(eRow);
    }
  }
  os << "LP: |vect_incd2|=" << vect_incd2.size() << "\n";
  MyMatrix<T> MatIncd2 = MatrixFromVectorFamily(vect_incd2);
  int rnk_incd2 = RankMat(MatIncd2);
  os << "LP: rnk_incd2=" << rnk_incd2 << "\n";
#endif
  int n_non_zero = 0;
  for (int i_row = 0; i_row < n_row; i_row++) {
    if (eSolDual(i_row) > 0) {
      n_non_zero += 1;
    }
  }
  int target_dim = n_col - 1;
#ifdef DEBUG_LINEAR_PROGRAM
  os << "LP: n_non_zero=" << n_non_zero << " n_col=" << n_col
     << " target_dim=" << target_dim << " rnk_incd=" << rnk_incd << "\n";
#endif
  if (n_non_zero == target_dim) {
    return true;
  }
  return false;
}

/*
  The Interior point obtained by a regular linear program corresponds
  essentially to just one vertex solution of the linear program. This code gives
  a solution with is geometrically unique:
  * If you permute the rows then you get the same point.
  * If you apply a matrix transformation to the rows of FAC then the interior
  point is changed by the contragredient action. Note that rescaling each rows
  by some scalar which varies from row to row will get you a different interior
  point.
  ---
  The whole computation is done with only linear programming though lots of it
  and does not involve full dual description.
 */
template <typename T>
MyVector<T> GetGeometricallyUniqueInteriorPoint(MyMatrix<T> const &FAC,
                                                std::ostream &os) {
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  auto check_return = [&](MyVector<T> const &Vret) -> void {
    int len = FAC.rows();
    for (int i = 0; i < len; i++) {
      MyVector<T> eFAC = GetMatrixRow(FAC, i);
      os << "i=" << i << " |eFAC|=" << eFAC.size() << " |retV|=" << Vret.size()
         << "\n";
      T scal = eFAC.dot(Vret);
      if (scal <= 0) {
        std::cerr << "We should have strictly positive scalar product\n";
        throw TerminalException{1};
      }
    }
  };
  os << "LP: GetGeometricallyUniqueInteriorPoint(GGUIP) : beginning\n";
#endif
  int n_rows = FAC.rows();
  int n_cols = FAC.cols();
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  os << "LP: GGUIP, n_rows=" << n_rows << " n_cols=" << n_cols << "\n";
#endif
  MyMatrix<T> ListInequalities = ZeroMatrix<T>(n_rows, n_cols + 1);
  MyVector<T> ToBeMinimized = ZeroVector<T>(n_cols + 1);
  for (int i_row = 0; i_row < n_rows; i_row++) {
    ListInequalities(i_row, 0) = -1;
    for (int i_col = 0; i_col < n_cols; i_col++)
      ListInequalities(i_row, i_col + 1) = FAC(i_row, i_col);
    for (int i_col = 0; i_col <= n_cols; i_col++)
      ToBeMinimized(i_col) += ListInequalities(i_row, i_col);
  }
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  os << "LP: GGUIP, CDD_LinearProgramming, before\n";
#endif
  LpSolution<T> eSol =
      CDD_LinearProgramming(ListInequalities, ToBeMinimized, os);
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  os << "LP: GGUIP, CDD_LinearProgramming, after\n";
#endif
  if (!eSol.PrimalDefined || !eSol.DualDefined) {
    std::cerr << "Failed to find an interior point by linear programming\n";
    std::cerr << "Maybe the cone is actually not full dimensional\n";
    throw TerminalException{1};
  }
  bool is_non_degenerate =
      TestCriterionNonDegenerate(ListInequalities, ToBeMinimized, eSol, os);
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  os << "LP: GGUIP, is_non_degenerate=" << is_non_degenerate << "\n";
#endif
  if (is_non_degenerate) {
    MyVector<T> V = eSol.DirectSolution;
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
    os << "LP: GGUIP, returns a solution\n";
    os << "LP: GGUIP V=" << StringVectorGAP(V) << " |V|=" << V.size()
       << " case 1\n";
    check_return(V);
#endif
    // The happy scenario where the linear programming directly gets us
    // something unique
    return V;
  }
  // Now computing the kernel of the attained kernel.
  ToBeMinimized(0) -= eSol.OptimalValue;
  MyMatrix<T> NSP = NullspaceMatSingleVector(ToBeMinimized);
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  os << "LP: GGUIP, |NSP|=" << NSP.rows() << " / " << NSP.cols() << "\n";
#endif
  MyMatrix<T> ListIneqRed = ListInequalities * NSP.transpose();
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  os << "LP: GGUIP, |ListIneqRed|=" << ListIneqRed.rows() << " / "
     << ListIneqRed.cols() << "\n";
#endif
  // Now the subspace built by the matching inequalities. It is of strictly
  // smaller dimension since otherwise the solution is not optimal.
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  os << "LP: GGUIP, LinearDeterminedByInequalitiesAndIndices, before\n";
  os << "LP: ListIneqRed=\n";
  WriteMatrix(os, ListIneqRed);
#endif
  std::pair<MyMatrix<T>, Face> pair =
      LinearDeterminedByInequalitiesAndIndices(ListIneqRed, os);
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  os << "LP: GGUIP, LinearDeterminedByInequalitiesAndIndices, after\n";
  os << "LP: pair.second=" << StringFace(pair.second)
     << " |pair.first|=" << pair.first.rows() << " / " << pair.first.cols()
     << " RankMat(pair.first)=" << RankMat(pair.first) << "\n";
  if (pair.first.rows() == 0) {
    std::cerr << "The dimension is 0 which is impossible\n";
    throw TerminalException{1};
  }
#endif
  if (pair.first.rows() == 1) {
    MyVector<T> V = eSol.DirectSolution;
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
    os << "LP: GGUIP, returns a solution\n";
    os << "LP: GGUIP V=" << StringVectorGAP(V) << " |V|=" << V.size()
       << " case 2\n";
    check_return(V);
#endif
    // The happy scenario where the linear programming directly gets us
    // something unique
    return V;
  }
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  os << "LP: GGUIP, Recurring\n";
#endif
  std::vector<int> ListIdxNZ;
  int len = ListIneqRed.rows();
  for (int u = 0; u < len; u++) {
    if (pair.second[u] == 0) {
      ListIdxNZ.push_back(u);
    }
  }
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  os << "LP: GGUIP, |ListIdxNZ|=" << ListIdxNZ.size() << "\n";
#endif
  MyMatrix<T> ListIneqRedB =
      SelectRow(ListIneqRed, ListIdxNZ) * pair.first.transpose();
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  os << "LP: GGUIP, |ListIneqRedB|=" << ListIneqRedB.rows() << " / "
     << ListIneqRedB.cols() << "\n";
  os << "LP: ListIneqRedB=\n";
  WriteMatrix(os, ListIneqRedB);
#endif
  MyVector<T> SolSubspaceB =
      GetGeometricallyUniqueInteriorPoint(ListIneqRedB, os);
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  os << "LP: GGUIP, We have SolSubspaceB=" << StringVectorGAP(SolSubspaceB)
     << "\n";
#endif
  MyVector<T> SolSubspace = pair.first.transpose() * SolSubspaceB;
  MyVector<T> V = NSP.transpose() * SolSubspace;
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  if (V(0) <= 0) {
    std::cerr << "The first coordinate should be strictly positive\n";
    throw TerminalException{1};
  }
#endif
  MyVector<T> Vret(n_cols);
  for (int i = 0; i < n_cols; i++) {
    Vret(i) = V(i + 1) / V(0);
  }
#ifdef DEBUG_GEOMETRICALLY_UNIQUE
  os << "LP: GGUIP, eSolFinal found\n";
  os << "LP: GGUIP Vret=" << StringVectorGAP(Vret) << " |V|=" << Vret.size()
     << " case 3\n";
  check_return(Vret);
#endif
  return Vret;
}

template <typename T>
MyVector<T> GetSpaceInteriorPoint(MyMatrix<T> const &FAC,
                                  MyMatrix<T> const &Equa, std::ostream &os) {
  MyMatrix<T> NSP = NullspaceTrMat(Equa);
  MyMatrix<T> FACred = FAC * NSP.transpose();
  MyVector<T> eVectInt = GetSpaceInteriorPoint_Basic(FACred, os);
  MyVector<T> TheSol = NSP.transpose() * eVectInt;
#ifdef DEBUG_LINEAR_PROGRAM
  MyVector<T> ListScal_FAC = FAC * TheSol;
  for (int i_row = 0; i_row < FAC.rows(); i_row++) {
    T eScal = ListScal_FAC(i_row);
    if (eScal <= 0) {
      std::cerr << "We have eScal=" << eScal << " which should be positive\n";
      throw TerminalException{1};
    }
  }
  MyVector<T> ListScal_Equa = Equa * TheSol;
  for (int i_row = 0; i_row < Equa.rows(); i_row++) {
    T eScal = ListScal_Equa(i_row);
    if (eScal != 0) {
      std::cerr << "We have eScal=" << eScal << " which should be zero\n";
      throw TerminalException{1};
    }
  }
#endif
  return TheSol;
}

template <typename T>
MyVector<T> GetSpaceInteriorPointFace(MyMatrix<T> const &FAC, Face const &f,
                                      std::ostream &os) {
  int n_row = FAC.rows();
  int n = FAC.cols();
  int n_f = f.size();
  if (n_row != n_f) {
    std::cerr << "FAC and f have incoherent lengths\n";
    throw TerminalException{1};
  }
  int cnt = f.count();
  MyMatrix<T> FACred(n_row - cnt, n);
  MyMatrix<T> Equa(cnt, n);
  int pos_fac = 0;
  int pos_equa = 0;
  for (int i_row = 0; i_row < n_row; i_row++) {
    if (f[i_row] == 0) {
      FACred.row(pos_fac) = FAC.row(i_row);
      pos_fac++;
    } else {
      Equa.row(pos_equa) = FAC.row(i_row);
      pos_equa++;
    }
  }
  return GetSpaceInteriorPoint(FACred, Equa, os);
}

template <typename T>
MyVector<T> GetSpaceInteriorPointFacet(MyMatrix<T> const &FAC,
                                       int const &iFacet, std::ostream &os) {
  int nbRow = FAC.rows();
  Face f(nbRow);
  f[iFacet] = 1;
  return GetSpaceInteriorPointFace(FAC, f, os);
}

template <typename T> struct EmbeddedPolytope {
  MyMatrix<T> LinSpace;
  MyMatrix<T> ListIneq;
};

template <typename T>
EmbeddedPolytope<T> ComputeEmbeddedPolytope(MyMatrix<T> const &ListIneq,
                                            MyMatrix<T> const &ListEqua) {
  MyMatrix<T> NSP = NullspaceTrMat(ListEqua);
  MyMatrix<T> ListIneqRed = ListIneq * (NSP.transpose());
  MyMatrix<T> LinSpa = LinearDeterminedByInequalities(ListIneqRed);
  MyMatrix<T> FinalSpace = (LinSpa.transpose()) * NSP;
  std::optional<MyMatrix<T>> eRes = AffinizeSubspace(FinalSpace);
  if (!eRes) {
    std::cerr << "Call to ComputeEmbeddedPolytope failed\n";
    std::cerr << "Because the affinization operation failed\n";
    throw TerminalException{1};
  }
  MyMatrix<T> FinalSpace_Aff = *eRes;
  MyMatrix<T> ListIneqFinal = ListIneq * (FinalSpace_Aff.transpose());
  return {FinalSpace_Aff, ListIneqFinal};
}

template <typename T>
vectface Kernel_GetFullRankFacetSet(const MyMatrix<T> &EXT, std::ostream &os) {
  size_t dim = EXT.cols();
  size_t n_rows = EXT.rows();
  if (dim == 2) {
    if (n_rows != 2) {
      std::cerr << "SAMP: In dimension 2, the cone should have exactly\n";
      std::cerr << "SAMP: two extreme rays\n";
      throw TerminalException{1};
    }
    vectface vf_ret(2);
    Face f1(2), f2(2);
    f1[0] = 1;
    f2[1] = 1;
    vf_ret.push_back(f1);
    vf_ret.push_back(f2);
    return vf_ret;
  }
#ifdef TIMINGS_FULL_RANK_FACET_SET
  MicrosecondTime time;
#endif
#ifdef DEBUG_FULL_RANK_FACET_SET
  os << "SAMP: Before Kernel_FindSingleVertex\n";
#endif
  Face eSet = Kernel_FindSingleVertex(EXT, os);
#ifdef TIMINGS_FULL_RANK_FACET_SET
  os << "|SAMP: Kernel_FindSingleVertex|=" << time << "\n";
#endif
  // Here we use a trick that the ColumnReduction will select the first column
  // and so will return a matrix that is polytopal
  MyMatrix<T> EXTsel_pre1 = SelectRow(EXT, eSet);
  std::vector<int> l_cols = ColumnReductionSet(EXTsel_pre1);
  MyMatrix<T> EXTsel_pre2 = SelectColumn(EXTsel_pre1, l_cols);
  MyMatrix<T> EXTsel = SetIsobarycenter(EXTsel_pre2); // Essential step, selecting a subset offsets the isobarycenter.
#ifdef TIMINGS_FULL_RANK_FACET_SET
  os << "|SAMP: SelectRow / SelectColumn / ColumnReductionSet|=" << time << "\n";
#endif
#ifdef DEBUG_FULL_RANK_FACET_SET
  os << "SAMP: |EXTsel|=" << EXTsel.rows() << " / " << EXTsel.cols()
     << " rnk=" << RankMat(EXTsel) << "\n";
#endif
#ifdef SANITY_CHECK_FULL_RANK_FACET_SET
  if (!IsPolytopal(EXTsel)) {
    std::cerr << "SAMP: The configuration EXTsel is not polytopal\n";
    throw TerminalException{1};
  }
#endif
  vectface ListRidge = Kernel_GetFullRankFacetSet(EXTsel, os);
#ifdef TIMINGS_FULL_RANK_FACET_SET
  os << "|SAMP: ListRidge|=" << time << "\n";
#endif
#ifdef DEBUG_FULL_RANK_FACET_SET
  os << "SAMP: We have ListRidge\n";
#endif
  SimplifiedFlippingFramework<T> RPLlift(EXT, eSet, os);
#ifdef TIMINGS_FULL_RANK_FACET_SET
  os << "|SAMP: FlippingFramework|=" << time << "\n";
#endif
#ifdef DEBUG_FULL_RANK_FACET_SET
  os << "SAMP: We have FlippingFramework\n";
#endif
  vectface vf_ret(n_rows);
  vf_ret.push_back(eSet);
  for (auto &eRidge : ListRidge) {
    Face eFace = RPLlift.FlipFace(eRidge);
    vf_ret.push_back(eFace);
  }
#ifdef TIMINGS_FULL_RANK_FACET_SET
  os << "|SAMP: vf_ret|=" << time << "\n";
#endif
#ifdef DEBUG_FULL_RANK_FACET_SET
  os << "SAMP: We have vf_ret\n";
#endif
#ifdef SANITY_CHECK_FULL_RANK_FACET_SET
  MyMatrix<T> FACsamp(vf_ret.size(), dim);
  int pos = 0;
  for (auto &face : vf_ret) {
    MyVector<T> eFAC = FindFacetInequality(EXT, face);
    AssignMatrixRow(FACsamp, pos, eFAC);
    pos++;
  }
  if (RankMat(FACsamp) != static_cast<int>(dim)) {
    std::cerr << "SAMP: Failed to find a ful rank vector configuration\n";
    throw TerminalException{1};
  }
#endif
  return vf_ret;
}

template <typename T>
vectface GetFullRankFacetSet(const MyMatrix<T> &EXT, std::ostream &os) {
#ifdef TIMINGS_FULL_RANK_FACET_SET
  MicrosecondTime time;
#endif
  MyMatrix<T> EXT_B = ColumnReduction(EXT);
#ifdef TIMINGS_FULL_RANK_FACET_SET
  os << "|SAMP: ColumnReduction|=" << time << "\n";
#endif
#ifdef DEBUG_FULL_RANK_FACET_SET
  os << "SAMP: Before Polytopization\n";
#endif
  MyMatrix<T> EXT_C = Polytopization(EXT_B, os);
#ifdef TIMINGS_FULL_RANK_FACET_SET
  os << "|SAMP: Polytopization|=" << time << "\n";
#endif
  MyMatrix<T> EXT_D = SetIsobarycenter(EXT_C);
#ifdef TIMINGS_FULL_RANK_FACET_SET
  os << "|SAMP: SetIsobarycenter|=" << time << "\n";
#endif
#ifdef DEBUG_FULL_RANK_FACET_SET
  os << "SAMP: Before Kernel_GetFullRankFacetSet\n";
#endif
  vectface vf = Kernel_GetFullRankFacetSet(EXT_D, os);
#ifdef TIMINGS_FULL_RANK_FACET_SET
  os << "|SAMP: Kernel_GetFullRankFacetSet|=" << time << "\n";
#endif
  return vf;
}

// clang-format off
#endif  // SRC_POLY_POLY_LINEARPROGRAMMING_H_
// clang-format on

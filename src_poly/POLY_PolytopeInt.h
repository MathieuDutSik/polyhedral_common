// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_POLYTOPEINT_H_
#define SRC_POLY_POLY_POLYTOPEINT_H_

// clang-format off
#include "COMB_Combinatorics.h"
#include "POLY_LinearProgramming.h"
#include <limits>
#include <vector>
// clang-format on

template <typename T, typename Tint, typename Finsert>
void Kernel_GetListIntegralPoint(MyMatrix<T> const &FAC, MyMatrix<T> const &EXT,
                                 Finsert f_insert, std::ostream &os) {
  int nbVert = EXT.rows();
  int n = EXT.cols();
  int dim = n - 1;
  MyMatrix<T> EXTnorm(nbVert, n);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    T eVal = EXT(iVert, 0);
    if (eVal <= 0) {
      std::cerr << "We have eVal=" << eVal << "\n";
      std::cerr << "which means the enumeration of integer points does not "
                   "make sense\n";
      throw TerminalException{1};
    }
    T eValInv = 1 / eVal;
    for (int i = 0; i < n; i++)
      EXTnorm(iVert, i) = eValInv * EXT(iVert, i);
  }
  std::vector<Tint> ListLow(dim);
  std::vector<Tint> ListUpp(dim);
  std::vector<int> ListSize(dim);
  Tint eProd = 1;
  for (int iDim = 0; iDim < dim; iDim++) {
    T val = EXTnorm(0, iDim + 1);
    Tint eLow = UniversalCeilScalarInteger<Tint, T>(val);
    Tint eUpp = UniversalFloorScalarInteger<Tint, T>(val);
    for (int iVert = 1; iVert < nbVert; iVert++) {
      T val_B = EXTnorm(iVert, iDim + 1);
      Tint eLow_B = UniversalCeilScalarInteger<Tint, T>(val_B);
      Tint eUpp_B = UniversalFloorScalarInteger<Tint, T>(val_B);
      if (eLow_B < eLow)
        eLow = eLow_B;
      if (eUpp_B > eUpp)
        eUpp = eUpp_B;
    }
    ListLow[iDim] = eLow;
    ListUpp[iDim] = eUpp;
    Tint eSiz = 1 + eUpp - eLow;
    ListSize[iDim] = UniversalScalarConversion<int, Tint>(eSiz);
    eProd *= eSiz;
  }
  os << "ListBound =";
  for (int iDim = 0; iDim < dim; iDim++)
    os << " [" << ListLow[iDim] << "," << ListUpp[iDim] << "]";
  os << " dim=" << dim << " eProd=" << eProd << "\n";
  int nbFac = FAC.rows();
  auto IsCorrect = [&](MyVector<Tint> const &eVect) -> bool {
    for (int iFac = 0; iFac < nbFac; iFac++) {
      T eScal = FAC(iFac, 0);
      for (int i = 0; i < dim; i++)
        eScal += FAC(iFac, i + 1) * eVect(i);
      if (eScal < 0)
        return false;
    }
    return true;
  };
  MyVector<Tint> ePoint(dim);
  BlockIterationMultiple BlIter(ListSize);
  for (auto const &eVect : BlIter) {
    for (int iDim = 0; iDim < dim; iDim++)
      ePoint(iDim) = ListLow[iDim] + eVect[iDim];
    bool test = IsCorrect(ePoint);
    if (test) {
      bool retval = f_insert(ePoint);
      if (!retval)
        return;
    }
  }
}

template <typename T, typename Tint>
std::vector<MyVector<Tint>> GetListIntegralPoint(MyMatrix<T> const &FAC,
                                                 MyMatrix<T> const &EXT,
                                                 std::ostream &os) {
  std::vector<MyVector<Tint>> ListPoint;
  auto f_insert = [&](const MyVector<Tint> &ePoint) -> bool {
    ListPoint.push_back(ePoint);
    return true;
  };
  Kernel_GetListIntegralPoint<T, Tint, decltype(f_insert)>(FAC, EXT, f_insert,
                                                           os);
  return ListPoint;
}

template <typename T, typename Tint, typename Finsert>
void Kernel_GetListIntegralPoint_LP(MyMatrix<T> const &FAC, Finsert f_insert,
                                    std::ostream &os) {
  //
  // Basic functionality
  //
  size_t dim = FAC.cols() - 1;
  size_t n_row = FAC.rows();
  std::vector<Tint> ListLow(dim);
  std::vector<Tint> ListUpp(dim);
  // Getting the bounds with the coordinates i<pos1 being fixed with values
  // set by eVert. The values of the coordinate pos2 is then computed.
  auto set_bound = [&](const MyVector<Tint> &ePoint,
                       const size_t &pos) -> void {
    MyMatrix<T> FACred(n_row, 1 + dim - pos);
    for (size_t i_row = 0; i_row < n_row; i_row++) {
      T val = FAC(i_row, 0);
      for (size_t i = 0; i < pos; i++)
        val += FAC(i_row, 1 + i) * ePoint(i);
      FACred(i_row, 0) = val;
      for (size_t i = 0; i < dim - pos; i++)
        FACred(i_row, 1 + i) = FAC(i_row, 1 + pos + i);
    }
    MyVector<T> Vminimize = ZeroVector<T>(1 + dim - pos);
    LpSolution<T> eSol;
    for (size_t i = 0; i < dim - pos; i++) {
      Vminimize(1 + i) = 1;
      eSol = CDD_LinearProgramming(FACred, Vminimize, os);
      if (!eSol.DualDefined || !eSol.PrimalDefined) {
        std::cerr << "eSol.DualDefined=" << eSol.DualDefined
                  << " eSol.PrimalDefined=" << eSol.PrimalDefined << "\n";
        std::cerr << "Failed computation of ListLow at i=" << i
                  << " which means that\n";
        std::cerr << "the polytope is unbounded and thus the integer point "
                     "enumeration will not work\n";
        throw TerminalException{1};
      }
      ListLow[i + pos] = UniversalCeilScalarInteger<Tint, T>(eSol.OptimalValue);
      //
      Vminimize(1 + i) = -1;
      eSol = CDD_LinearProgramming(FACred, Vminimize, os);
      if (!eSol.DualDefined || !eSol.PrimalDefined) {
        std::cerr << "eSol.DualDefined=" << eSol.DualDefined
                  << " eSol.PrimalDefined=" << eSol.PrimalDefined << "\n";
        std::cerr << "Failed computation of ListUpp at i=" << i
                  << " which means that\n";
        std::cerr << "the polytope is unbounded and thus the integer point "
                     "enumeration will not work\n";
        throw TerminalException{1};
      }
      ListUpp[i + pos] =
          UniversalFloorScalarInteger<Tint, T>(-eSol.OptimalValue);
      //
      Vminimize(1 + i) = 0;
    }
  };
  auto get_number_poss = [&](const size_t &pos) -> size_t {
    size_t nb_pos = 1;
    for (size_t i = pos; i < dim; i++) {
      size_t len = size_t(
          UniversalScalarConversion<int, Tint>(1 + ListUpp[i] - ListLow[i]));
      size_t new_nb_pos = nb_pos * len;
      if (new_nb_pos < nb_pos) {
        // Case of going overflow
        return std::numeric_limits<size_t>::max();
      }
      nb_pos = new_nb_pos;
    }
    return nb_pos;
  };
  auto IsCorrect = [&](MyVector<Tint> const &eVect) -> bool {
    for (size_t i_row = 0; i_row < n_row; i_row++) {
      T eScal = FAC(i_row, 0);
      for (size_t i = 0; i < dim; i++)
        eScal += FAC(i_row, i + 1) * eVect(i);
      if (eScal < 0)
        return false;
    }
    return true;
  };
  //
  // Setting up the initial entries
  //
  // size_crit : This empirical value is obtained from an analysis
  // of a number of cases.
  size_t crit_siz = 10000;
  MyVector<Tint> ePoint(dim);
  set_bound(ePoint, 0);
  std::cerr << "ListBound =";
  for (size_t iDim = 0; iDim < dim; iDim++)
    std::cerr << " [" << ListLow[iDim] << "," << ListUpp[iDim] << "]";
  std::cerr << "\n";
  size_t pos = 0;
  //
  // While loop for iterating
  //
  while (true) {
    size_t nb_poss = get_number_poss(pos);
    if (nb_poss < crit_siz) {
      std::vector<int> ListSize(dim - pos);
      size_t len = dim - pos;
      for (size_t i = 0; i < len; i++) {
        ListSize[i] = UniversalScalarConversion<int, Tint>(
            1 + ListUpp[i + pos] - ListLow[i + pos]);
      }
      BlockIterationMultiple BlIter(ListSize);
      for (auto const &eVect : BlIter) {
        for (size_t i = 0; i < len; i++)
          ePoint(pos + i) = ListLow[pos + i] + eVect[i];
        bool test = IsCorrect(ePoint);
        if (test) {
          bool retval = f_insert(ePoint);
          if (!retval)
            return;
        }
      }
      if (pos == 0)
        break;
      pos--;
      // Now we need to increase the positions
      while (true) {
        if (ePoint(pos) < ListUpp[pos]) {
          ePoint(pos) += 1;
          pos++;
          set_bound(ePoint, pos);
          break;
        }
        if (pos == 0)
          break;
        pos--;
      }
      if (pos == 0)
        break;
    } else { // If the number of cases is large, then implicitly, we can go
             // deeper
      ePoint(pos) = ListLow[pos];
      pos++;
      set_bound(ePoint, pos);
    }
  }
}

template <typename T, typename Tint>
std::vector<MyVector<Tint>> GetListIntegralPoint_LP(MyMatrix<T> const &FAC,
                                                    std::ostream &os) {
  std::vector<MyVector<Tint>> ListPoint;
  auto f_insert = [&](const MyVector<Tint> &ePoint) -> bool {
    ListPoint.push_back(ePoint);
    return true;
  };
  Kernel_GetListIntegralPoint_LP<T, Tint, decltype(f_insert)>(FAC, f_insert,
                                                              os);
  return ListPoint;
}

/*
  We are looking for the 01 solutions X of the equation XA = B.
  We have m equations and n unknowns.

  -

  What could be a polyhedral approach to the problem:
  * Compute the Q-vector space of possible solutions of dimension r.
  * Define the polytope of possible solutions with 0 <= x_i <= 1.
    It will have dimension r and have 2n defining inequalities.
  * Backtrack on the possible partial choices for x_i, that is for
    sure. We have a C-style data set containing three possible values
    for the indices:
    - 0 for 0
    - 1 for 1
    - 2 for unknown
  * Compute by linear programming the upper bound of function x_i over
    the corresponding restricted polytopes. Both 0 <= x_i  and x_i <= 1
    have to be considered. If only one extreme value is attained, then
    we can set it to be equal and iterate the refinement strategy.
    If none are attained then no solution are possible and we backtrack.
  * Could we use redundancy informations?
    It would seem yes. If some redundancy is like f = C + sum C_i f_i
    with C > 0 then we can detect it right away by the linear programming.
    If C = 0, then the choice of some C would imply many other choices.
    So, said variable becomes pivotal. And the cost does not appear to
    be higher than a linear program. We simply use the dual solution.
  * Can we do a backtracking on pow(2,r) choices?
    It should be doable. We select r independent coordinate {i_1, ..., i_r}.
    We iterate over their possible values, and find bz matrix product their
    possible values and keep the one that are integral.
 */

/*
  The matrix A(m,n) and B(n)   and we look for x A = B
  We search for solution x in {0,1}^m.

 */
/*
template <typename T, typename Tint, typename Finsert>
void Kernel_Enumerate01_solutions_Iterate(MyMatrix<T> const &A, MyVector<T>
const &B, Finsert f_insert) { int m = A.rows(); int n = A.cols(); MyMatrix<T>
Bresidual(m+1,n); for (int i=0; i<n; i++) { Bresidual(0,i) = B(i);
  }
  MyVector<T> TheSol;
  int level_select = 0;
  auto set_residual=[&](int level, int coef) -> void {
    if (coef == 0) {
      for (int i=0; i<n; i++)
        Bresidual(level+1,i) = Bresidual(level,i);
    } else {
      for (int i=0; i<n; i++)
        Bresidual(level+1,i) = Bresidual(level,i) - A(level,i);
    }
  };
  auto is_feasible=[&](int level) -> bool {
    // Filtering by the GCD
    for (int i=0; i<n; i++) {
      T TheGcd = A(level,i);
      for (int j=0; j<m-level; j++)
        TheGcd = GcdPair(TheGcd, A(j+level,i));
      T TheRes = ResInt(Bresidual(level+1,i), TheGcd);
      if (TheRes != 0)
        return false;
    }
    return true;
  };
  while(true) {



  }

}
*/

template <typename T, typename Tint, typename Finsert>
void Kernel_Enumerate01_solutions_Iterate(MyMatrix<T> const &A,
                                          MyVector<T> const &B,
                                          [[maybe_unused]] Finsert f_insert) {
  std::optional<MyVector<T>> opt = SolutionIntMat(A, B);
  if (!opt)
    return;
  MyVector<T> const &ePt = *opt;
  MyMatrix<T> NSP = NullspaceIntMat(A);
}

// clang-format off
#endif  // SRC_POLY_POLY_POLYTOPEINT_H_
// clang-format on

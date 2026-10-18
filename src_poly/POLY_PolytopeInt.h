// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_POLYTOPEINT_H_
#define SRC_POLY_POLY_POLYTOPEINT_H_

// clang-format off
#include "COMB_Combinatorics.h"
#include "POLY_LinearProgramming.h"
#include "ClassicLLL.h"
#include <limits>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_POLYTOPE_INT
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_POLYTOPE_INT
#endif

template<typename T, typename Tint>
MyMatrix<Tint> lll_reduction_matrix(MyMatrix<T> const& FAC, std::ostream& os) {
  int dim = FAC.cols() - 1;
  int nbRow = FAC.rows();
  MyMatrix<T> TheGram = ZeroMatrix<T>(dim, dim);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        AddMul(TheGram(i,j), FAC(iRow, i+1), FAC(iRow, j+1));
      }
    }
  }
  LLLreduction<T, Tint> res = LLLreducedBasis<T, Tint>(TheGram, os);
  return TransposedMat(res.Pmat);
}

template<typename T, typename Tint>
MyMatrix<T> apply_lll_reduction(MyMatrix<T> const&FAC, MyMatrix<Tint> const& Pmat, [[maybe_unused]] std::ostream& os) {
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T,Tint>(Pmat);
  int nbRow = FAC.rows();
  int dim = FAC.cols() - 1;
  MyMatrix<T> FACred(nbRow, dim + 1);
  MyVector<T> V1(dim), V2(dim);
  for (int iRow=0; iRow<nbRow; iRow++) {
    for (int i=0; i<dim; i++) {
      V1(i) = FAC(iRow, i+1);
    }
    FACred(iRow, 0) = FAC(iRow, 0);
    V2 = Pmat_T.transpose() * V1;
    for (int i=0; i<dim; i++) {
      FACred(iRow, i+1) = V2(i);
    }
  }
#ifdef DEBUG_POLYTOPE_INT
  os << "POLYINT: L1(FAC)    = " << L1_norm_mat(FAC) << " L1(FACred) = " << L1_norm_mat(FACred) << "\n";
#endif
  return FACred;
}


// The bounds of the coordinates pos..dim-1 over the polytope FACin. Returns
// false when the polytope is empty, and throws when it is unbounded (the
// integer point enumeration needs a polytope).
template<typename T, typename Tint>
bool set_bound_lp(MyMatrix<T> const& FACin, int const& dim, int const& pos, std::vector<Tint> & ListUpp, std::vector<Tint> & ListLow, std::ostream& os) {
#ifdef DEBUG_POLYTOPE_INT
  os << "POLYINT: set_bound_lp, step 1\n";
#endif
  MyVector<T> Vminimize = ZeroVector<T>(1 + dim - pos);
  LpSolution<T> eSol;
#ifdef DEBUG_POLYTOPE_INT
  os << "POLYINT: set_bound_lp, step 2\n";
#endif
  for (int i = 0; i < dim - pos; i++) {
    Vminimize(1 + i) = 1;
    eSol = SIMPLEX_LinearProgramming(FACin, Vminimize, os);
    if (!eSol.DirectSolution) {
      // Infeasible: the polytope is empty.
      return false;
    }
    if (!eSol.DualSolution) {
      std::cerr << "eSol.DualDefined=" << eSol.DualSolution.has_value()
                << " eSol.PrimalDefined=" << eSol.DirectSolution.has_value() << "\n";
      std::cerr << "Failed computation of ListLow at i=" << i
                << " which means that\n";
      std::cerr << "the polytope is unbounded and thus the integer point "
        "enumeration will not work\n";
      throw TerminalException{1};
    }
    ListLow[i + pos] = UniversalCeilScalarInteger<Tint, T>(eSol.OptimalValue);
    //
    Vminimize(1 + i) = -1;
    eSol = SIMPLEX_LinearProgramming(FACin, Vminimize, os);
    if (!eSol.DualSolution || !eSol.DirectSolution) {
      std::cerr << "eSol.DualDefined=" << eSol.DualSolution.has_value()
                << " eSol.PrimalDefined=" << eSol.DirectSolution.has_value() << "\n";
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
#ifdef DEBUG_POLYTOPE_INT
  os << "POLYINT: set_bound_lp, step 3\n";
#endif
  return true;
}


template <typename T, typename Tint, typename Finsert>
void Kernel_GetListIntegralPoint_ITER_no_LLL(MyMatrix<T> const &FAC,
                                             Finsert f_insert,
                                             std::ostream &os) {
  int n = FAC.cols();
  int dim = n - 1;
  std::vector<Tint> ListLow(dim);
  std::vector<Tint> ListUpp(dim);
  if (!set_bound_lp(FAC, dim, 0, ListUpp, ListLow, os)) {
    // An empty polytope has no integer point.
    return;
  }
  std::vector<int> ListSize(dim);
  for (int iDim = 0; iDim < dim; iDim++) {
    Tint const& eLow = ListLow[iDim];
    Tint const& eUpp = ListUpp[iDim];
    Tint eSiz = 1 + eUpp - eLow;
    ListSize[iDim] = UniversalScalarConversion<int, Tint>(eSiz);
  }
#ifdef DEBUG_POLYTOPE_INT
  os << "ListBound =";
  for (int iDim = 0; iDim < dim; iDim++) {
    os << " [" << ListLow[iDim] << "," << ListUpp[iDim] << "]";
  }
  Tint eProd = 1;
  for (int iDim = 0; iDim < dim; iDim++) {
    Tint eSiz = ListSize[iDim];
    eProd *= eSiz;
  }
  os << " dim=" << dim << " eProd=" << eProd << "\n";
#endif
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
    for (int iDim = 0; iDim < dim; iDim++) {
      ePoint(iDim) = ListLow[iDim] + eVect[iDim];
    }
    bool test = IsCorrect(ePoint);
    if (test) {
      bool retval = f_insert(ePoint);
      if (!retval) {
        return;
      }
    }
  }
}

template <typename T, typename Tint, typename Finsert>
void Kernel_GetListIntegralPoint_ITER(MyMatrix<T> const &FAC,
                                      Finsert f_insert,
                                      std::ostream &os) {
  int dim = FAC.cols() - 1;
  MyMatrix<Tint> Pmat = lll_reduction_matrix<T,Tint>(FAC, os);
  MyMatrix<T> FACred = apply_lll_reduction(FAC, Pmat, os);
  MyVector<Tint> V(dim);
  auto f_insert_bis=[&](MyVector<Tint> const& Vred) -> bool {
    V = Pmat * Vred;
    return f_insert(V);
  };
  Kernel_GetListIntegralPoint_ITER_no_LLL<T,Tint,decltype(f_insert_bis)>(FACred, f_insert_bis, os);
}


template <typename T, typename Tint>
std::vector<MyVector<Tint>> GetListIntegralPoint_ITER_no_LLL(MyMatrix<T> const &FAC,
                                                             std::ostream &os) {
  std::vector<MyVector<Tint>> ListPoint;
  auto f_insert = [&](const MyVector<Tint> &ePoint) -> bool {
    ListPoint.push_back(ePoint);
    return true;
  };
  Kernel_GetListIntegralPoint_ITER_no_LLL<T, Tint, decltype(f_insert)>(FAC, f_insert,
                                                                       os);
  return ListPoint;
}

template <typename T, typename Tint>
std::vector<MyVector<Tint>> GetListIntegralPoint_ITER(MyMatrix<T> const &FAC,
                                                      std::ostream &os) {
  std::vector<MyVector<Tint>> ListPoint;
  auto f_insert = [&](const MyVector<Tint> &ePoint) -> bool {
    ListPoint.push_back(ePoint);
    return true;
  };
  Kernel_GetListIntegralPoint_ITER<T, Tint, decltype(f_insert)>(FAC, f_insert,
                                                                os);
  return ListPoint;
}



template <typename T, typename Tint, typename Finsert>
void Kernel_GetListIntegralPoint_LP_no_LLL(MyMatrix<T> const &FAC, Finsert f_insert,
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
                       const size_t &pos) -> bool {
    MyMatrix<T> FACred(n_row, 1 + dim - pos);
    for (size_t i_row = 0; i_row < n_row; i_row++) {
      T val = FAC(i_row, 0);
      for (size_t i = 0; i < pos; i++)
        val += FAC(i_row, 1 + i) * ePoint(i);
      FACred(i_row, 0) = val;
      for (size_t i = 0; i < dim - pos; i++)
        FACred(i_row, 1 + i) = FAC(i_row, 1 + pos + i);
    }
    return set_bound_lp(FACred, dim, pos, ListUpp, ListLow, os);
  };
  // Below the top level a slice cannot be empty: the fixed coordinate lies
  // between the bounds of the projection of the previous slice.
  auto set_bound_slice = [&](const MyVector<Tint> &ePoint,
                             const size_t &pos) -> void {
    [[maybe_unused]] bool feasible = set_bound(ePoint, pos);
#ifdef SANITY_CHECK_POLYTOPE_INT
    if (!feasible) {
      std::cerr << "POLYINT: the slice at pos=" << pos << " is empty\n";
      throw TerminalException{1};
    }
#endif
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
  if (!set_bound(ePoint, 0)) {
    // An empty polytope has no integer point.
    return;
  }
#ifdef DEBUG_POLYTOPE_INT
  os << "ListBound =";
  for (size_t iDim = 0; iDim < dim; iDim++)
    os << " [" << ListLow[iDim] << "," << ListUpp[iDim] << "]";
  os << "\n";
#endif
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
          set_bound_slice(ePoint, pos);
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
      set_bound_slice(ePoint, pos);
    }
  }
}

template <typename T, typename Tint, typename Finsert>
void Kernel_GetListIntegralPoint_LP(MyMatrix<T> const &FAC,
                                    Finsert f_insert,
                                    std::ostream &os) {
  int dim = FAC.cols() - 1;
  MyMatrix<Tint> Pmat = lll_reduction_matrix<T,Tint>(FAC, os);
  MyMatrix<T> FACred = apply_lll_reduction(FAC, Pmat, os);
  MyVector<Tint> V(dim);
  auto f_insert_bis=[&](MyVector<Tint> const& Vred) -> bool {
    V = Pmat * Vred;
    return f_insert(V);
  };
  Kernel_GetListIntegralPoint_LP_no_LLL<T,Tint,decltype(f_insert_bis)>(FACred, f_insert_bis, os);
}

template <typename T, typename Tint>
std::vector<MyVector<Tint>> GetListIntegralPoint_LP_no_LLL(MyMatrix<T> const &FAC,
                                                           std::ostream &os) {
  std::vector<MyVector<Tint>> ListPoint;
  auto f_insert = [&](const MyVector<Tint> &ePoint) -> bool {
    ListPoint.push_back(ePoint);
    return true;
  };
  Kernel_GetListIntegralPoint_LP_no_LLL<T, Tint, decltype(f_insert)>(FAC, f_insert,
                                                                     os);
  return ListPoint;
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


template <typename T, typename Tint>
std::vector<MyVector<Tint>> GetListIntegralPoint(MyMatrix<T> const &FAC, std::ostream &os) {
  // Use the iteration systematically. Combined with LLL, that is quite powerful.
  return GetListIntegralPoint_ITER<T,Tint>(FAC, os);
}

// clang-format off
#endif  // SRC_POLY_POLY_POLYTOPEINT_H_
// clang-format on

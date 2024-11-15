// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_REDUNDANCYELIMINATION_H_
#define SRC_POLY_POLY_REDUNDANCYELIMINATION_H_

// clang-format off
#include "GRP_GroupFct.h"
#include "POLY_LinearProgramming.h"
#include "POLY_cddlib.h"
#include <algorithm>
#include <set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ELIMINATION_REDUNDANCY
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ELIMINATION_REDUNDANCY
#endif

#ifdef PRINT
#define PRINT_ELIMINATION_REDUNDANCY
#endif

// Fairly expensive function. But useful for debugging
template <typename T>
std::vector<int> Kernel_GetNonRedundant_CDD(const MyMatrix<T> &M, std::ostream & os) {
  MyMatrix<T> Mred = ColumnReduction(M);
  MyMatrix<T> EXT = cdd::DualDescription(Mred, os);
  int n_row = Mred.rows();
  int n_vert = EXT.rows();
  int n_col = Mred.cols();
  std::vector<int> TheSel;
  for (int i_row = 0; i_row < n_row; i_row++) {
    std::vector<int> eIncd;
    for (int i_vert = 0; i_vert < n_vert; i_vert++) {
      T scal = 0;
      for (int i_col = 0; i_col < n_col; i_col++)
        scal += EXT(i_vert, i_col) * Mred(i_row, i_col);
      if (scal == 0)
        eIncd.push_back(i_vert);
    }
    MyMatrix<T> EXT_face = SelectRow(EXT, eIncd);
    if (RankMat(EXT_face) == n_col - 1)
      TheSel.push_back(i_row);
  }
  return TheSel;
}

template <typename T> MyMatrix<T> GetNonRedundant_CDD(const MyMatrix<T> &M, std::ostream& os) {
  return SelectRow(M, Kernel_GetNonRedundant_CDD(M, os));
}

template <typename T> struct FacetizationInfo {
  MyMatrix<T> EXT;
  MyMatrix<T> BoundingFAC;
};

template <typename T>
FacetizationInfo<T> FacetizationCone(MyMatrix<T> const &EXT,
                                     MyMatrix<T> const &BoundingFac, [[maybe_unused]] std::ostream &os) {
  int n_rows = EXT.rows();
  int n_cols = EXT.cols();
#ifdef DEBUG_ELIMINATION_REDUNDANCY
  os << "FacetizationCone : n_rows=" << n_rows << " n_cols=" << n_cols << "\n";
#endif
  MyMatrix<T> TheSum = ZeroMatrix<T>(1, n_cols);
  for (int i_row = 0; i_row < n_rows; i_row++)
    for (int i_col = 0; i_col < n_cols; i_col++)
      TheSum(0, i_col) += EXT(i_row, i_col);
#ifdef DEBUG_ELIMINATION_REDUNDANCY
  os << "FacetizationCone : TheSum\n";
#endif
  MyMatrix<T> EXTtot = Concatenate(TheSum, EXT);
#ifdef DEBUG_ELIMINATION_REDUNDANCY
  os << "FacetizationCone : EXTtot\n";
  os << "|EXTtot|=" << EXTtot.rows() << " / " << EXTtot.cols() << "\n";
#endif
  MyMatrix<T> EXTbas = RowReduction(EXTtot);
#ifdef DEBUG_ELIMINATION_REDUNDANCY
  os << "FacetizationCone : EXTbas\n";
  os << "|EXTbas|=" << EXTbas.rows() << " / " << EXTbas.cols() << "\n";
#endif
  MyMatrix<T> eInvMat = Inverse(EXTbas);
#ifdef DEBUG_ELIMINATION_REDUNDANCY
  os << "FacetizationCone : eInvMat\n";
#endif
  MyMatrix<T> EXT_ret = EXT * eInvMat;
#ifdef DEBUG_ELIMINATION_REDUNDANCY
  os << "FacetizationCone : EXT_ret\n";
#endif
  //
  MyMatrix<T> BoundingFac_ret = BoundingFac * TransposedMat(EXTbas);
#ifdef DEBUG_ELIMINATION_REDUNDANCY
  os << "FacetizationCone : BoundingFac_ret\n";
#endif
  return {std::move(EXT_ret), std::move(BoundingFac_ret)};
}

template <typename T>
std::vector<int> EliminationByRedundance_HitAndRun(MyMatrix<T> const &EXTin,
                                                   std::ostream &os) {
  MyMatrix<T> EXT = ColumnReduction(EXTin);
  int n_rows = EXT.rows();
  int n_cols = EXT.cols();
#ifdef DEBUG_ELIMINATION_REDUNDANCY
  os << "n_rows=" << n_rows << " n_cols=" << n_cols << "\n";
#endif
  MyMatrix<T> BoundingFac(0, n_cols);
  FacetizationInfo<T> RegCone = FacetizationCone(EXT, BoundingFac, os);
#ifdef DEBUG_ELIMINATION_REDUNDANCY
  os << "We have regCone\n";
#endif
  MyMatrix<T> EXT_work = RegCone.EXT;
  // return true if it is redundant. False if irredundant
  // If true, also returns a list of indices showing up positively in the
  // decomposition
  auto GetRedundancyInfo =
      [&](std::vector<int> const &ListIRow,
          int const &iRow) -> std::pair<bool, std::vector<int>> {
    MyMatrix<T> EXT_sel = SelectRow(EXT, ListIRow);
    MyVector<T> eRow = GetMatrixRow(EXT, iRow);
#ifdef PRINT_ELIMINATION_REDUNDANCY
    os << "-------------------\n";
    os << "H-representation\n";
    os << "begin\n";
    os << " " << EXT_sel.rows() << " " << EXT_sel.cols() << " integer\n";
    for (int i_r = 0; i_r < EXT_sel.rows(); i_r++) {
      for (int i_c = 0; i_c < EXT_sel.cols(); i_c++)
        os << " " << EXT_sel(i_r, i_c);
      os << "\n";
    }
    os << "end\n";
    os << "minimize\n";
    for (int i_c = 0; i_c < EXT_sel.cols(); i_c++)
      os << " " << eRow(i_c);
    os << "\n";
    os << "-------------------\n";
#endif
    LpSolution<T> eSol = CDD_LinearProgramming(EXT_sel, eRow, os);
    if (!eSol.DualDefined || !eSol.PrimalDefined) {
      return {false, {}};
    }
    if (eSol.OptimalValue < 0)
      return {false, {}};
    std::vector<int> ListIdx;
    int n_rows_ineq = ListIRow.size();
    for (int i_row = 0; i_row < n_rows_ineq; i_row++) {
      T e_val = -eSol.DualSolution(i_row);
#ifdef SANITY_CHECK_ELIMINATION_REDUNDANCY
      if (e_val < 0) {
        std::cerr << "The coefficient should be non-negative\n";
        for (int j_row = 0; j_row < n_rows_ineq; j_row++) {
          std::cerr << " " << eSol.DualSolution(j_row);
        }
        std::cerr << "\n";
        throw TerminalException{1};
      }
#endif
      if (e_val > 0)
        ListIdx.push_back(ListIRow[i_row]);
    }
    return {true, ListIdx};
  };
  //
  // The variables showing our current state
  //
  // -1: Still unknown 0: redundant 1: irredundant
  std::vector<int> RedundancyStatus(n_rows, -1);
  std::vector<int> ListKnownIrred;
  int nbFoundIrred = 0;
  RankTool<T> tool(n_cols);
  auto SetIRowIrredundant = [&](int const &i_row) -> void {
    RedundancyStatus[i_row] = 1;
    ListKnownIrred.push_back(i_row);
    nbFoundIrred++;
    MyVector<T> V = GetMatrixRow(EXT, i_row);
    tool.insert_if_indep(V);
  };
#ifdef DEBUG_ELIMINATION_REDUNDANCY
  os << "We have RankTool\n";
#endif
  //
  // Now computing one interior point.
  //
  MyVector<T> eVectInterior = GetSpaceInteriorPoint_Basic(EXT_work, os);
  MyVector<T> ListScalInterior = EXT_work * eVectInterior;
#ifdef DEBUG_ELIMINATION_REDUNDANCY
  for (int i_row = 0; i_row < n_rows; i_row++) {
    T eScal = ListScalInterior(i_row);
    if (eScal <= 0) {
      std::cerr << "Inconsistent scalar product\n";
      throw TerminalException{1};
    }
  }
#endif
  //
  // Now running the hit and run operations.
  //
  int N = 10;
  MyVector<T> eVect(n_cols);
  auto SetRandomVector = [&]() -> void {
    for (int i_col = 0; i_col < n_cols; i_col++) {
      int val = -N + random() % (2 * N + 1);
      eVect(i_col) = val;
    }
  };
  // Return the found facet if successful. Returns -1 otherwise.
  MyVector<T> ListScal;
  auto GetRandomOutsideVector_and_HitAndRun = [&]() -> int {
    ListScal.noalias() = EXT_work * eVect;
    auto HasOneViolatedFacet = [&](int const &h) -> bool {
      for (int i_row = 0; i_row < n_rows; i_row++) {
        T eScal = ListScal(i_row) - h * ListScalInterior(i_row);
        if (eScal < 0)
          return true;
      }
      return false;
    };
    auto GetSmallestValue = [&](int const &h) -> int {
      bool FoundUpper = false;
      // CurrUpper : specific value that will stay unused
      T CurrUpper(-1);
      int nbMatch = 0;
      int idxFound = -1;
      for (int i_row = 0; i_row < n_rows; i_row++) {
        T eScalDir = ListScal(i_row) - h * ListScalInterior(i_row);
        T eScal = ListScalInterior(i_row);
        if (eScalDir < 0) {
          // eScal + eScalDir * alpha = 0
          T UpperValue = -eScal / eScalDir;
          if (FoundUpper) {
            if (CurrUpper == UpperValue) {
              nbMatch++;
            } else {
              if (UpperValue < CurrUpper) {
                CurrUpper = UpperValue;
                nbMatch = 1;
                idxFound = i_row;
              }
            }
          } else {
            FoundUpper = true;
            CurrUpper = UpperValue;
            nbMatch = 1;
            idxFound = i_row;
          }
        }
      }
      if (nbMatch > 1)
        return -1;
      return idxFound;
    };
    int h = 0;
    while (true) {
      bool test = HasOneViolatedFacet(h);
#ifdef DEBUG_ELIMINATION_REDUNDANCY
      os << "HasOneViolatedFacet h=" << h << " test=" << test << "\n";
#endif
      if (test)
        return GetSmallestValue(h);
      h++;
    }
  };
  int nbRuns = 0;
#ifdef UNSET_TRIES
  for (int i_col = 0; i_col < n_cols; i_col++) {
    eVect(i_col) = 0;
  }
  for (int i_col = 0; i_col < n_cols; i_col++) {
    for (int i = 0; i < 2; i++) {
      eVect(i_col) = -1 + 2 * i;
      int idxFound = GetRandomOutsideVector_and_HitAndRun();
      if (idxFound != -1)
        if (RedundancyStatus[idxFound] == -1)
          SetIRowIrredundant(idxFound);
    }
    eVect(i_col) = 0;
  }
  for (int i_row = 0; i_row < n_rows; i_row++) {
    for (int i = 0; i < 2; i++) {
      int alpha = -1 + 2 * i;
      for (int i_col = 0; i_col < n_cols; i_col++)
        eVect(i_col) = EXT_work(i_row, i_col) * alpha;
      int idxFound = GetRandomOutsideVector_and_HitAndRun();
      if (idxFound != -1)
        if (RedundancyStatus[idxFound] == -1)
          SetIRowIrredundant(idxFound);
    }
  }
#endif
  while (true) {
#ifdef DEBUG_ELIMINATION_REDUNDANCY
    os << "nbRuns=" << nbRuns << " nbFoundIrred=" << nbFoundIrred << "\n";
#endif
    SetRandomVector();
    int idxFound = GetRandomOutsideVector_and_HitAndRun();
    if (idxFound != -1)
      if (RedundancyStatus[idxFound] == -1)
        SetIRowIrredundant(idxFound);
    nbRuns++;

    // Updating indexes
    // If we found less than n_cols then we know more can be found
    if (nbFoundIrred >= n_cols) {
      if (nbRuns > 100 * nbFoundIrred) {
        break;
      }
    }
  }
  //
  // Above we have determined an initial list of facets.
  // Now we use it to do the full enumeration
  //
  // Determine the status with absolute certainty.
  // ---If false then it is irredundant.
  // ---If true then it is redundant and retuns a list of candidates
  //    for finding a new facet.
  auto DetermineStatusSurely =
      [&](int const &i_row) -> std::pair<bool, std::vector<int>> {
    std::vector<int> ListIRow;
    for (int j_row = 0; j_row < n_rows; j_row++)
      if (RedundancyStatus[j_row] != 0 && j_row != i_row)
        ListIRow.push_back(j_row);
    std::pair<bool, std::vector<int>> ePair =
        GetRedundancyInfo(ListIRow, i_row);
    if (!ePair.first)
      return {false, {}};
    std::vector<int> LIdx;
    for (auto &eIdx : ePair.second)
      if (RedundancyStatus[eIdx] == -1)
        LIdx.push_back(eIdx);
    return {true, LIdx};
  };
  // If it is true then it is redundant.
  // If this function return false then the facet is redundant or not.
  auto FastStatusDetermination = [&](int const &i_row) -> bool {
    MyVector<T> V = GetMatrixRow(EXT, i_row);
    if (!tool.is_in_span(V))
      // We cannot conclude because the vector does not belong to the span.
      return false;
    std::vector<int> ListIRow;
    for (int j_row = 0; j_row < n_rows; j_row++)
      if (RedundancyStatus[j_row] == 1)
        ListIRow.push_back(j_row);
    return GetRedundancyInfo(ListIRow, i_row).first;
  };
  auto FindOneMoreFacet = [&](int const &i_start) -> void {
    std::vector<int> WorkLIdx{i_start};
    while (true) {
      std::set<int> NewCand;
      for (auto &eIdx : WorkLIdx) {
        std::pair<bool, std::vector<int>> ePair = DetermineStatusSurely(eIdx);
        if (ePair.first) {
          // The facet is redundant
          RedundancyStatus[eIdx] = 0;
          for (auto &fIdx : ePair.second)
            NewCand.insert(fIdx);
        } else {
          // The facet is irredundant. Mission accomplished.
          SetIRowIrredundant(eIdx);
          return;
        }
      }
      WorkLIdx.clear();
      for (auto &eIdx : NewCand)
        // Only those unconcluded need to be considered.
        if (RedundancyStatus[eIdx] == -1)
          WorkLIdx.push_back(eIdx);
#ifdef SANITY_CHECK_ELIMINATION_REDUNDANCY
      if (WorkLIdx.size() == 0) {
        std::cerr << "WorkLIdx is empty. Not allowed\n";
        throw TerminalException{1};
      }
#endif
    }
  };
  auto ProcessOnePoint = [&](int const &i_row) -> void {
    bool test = FastStatusDetermination(i_row);
    if (test) {
      // The heuristic works. Facet is redundant. We do conclude.
      RedundancyStatus[i_row] = 0;
      return;
    }
    // If the heuristic fails then it means that we miss one irreducible facet
    // and we should find at least one
    FindOneMoreFacet(i_row);
  };
  // Now processing all the data
  for (int i_row = 0; i_row < n_rows; i_row++)
    if (RedundancyStatus[i_row] == -1)
      ProcessOnePoint(i_row);
#ifdef SANITY_CHECK_ELIMINATION_REDUNDANCY
  for (int i_row = 0; i_row < n_rows; i_row++)
    if (RedundancyStatus[i_row] == -1) {
      std::cerr << "The algorithm failed to treat all the points\n";
      throw TerminalException{1};
    }
#endif
  std::sort(ListKnownIrred.begin(), ListKnownIrred.end());
  return ListKnownIrred;
}

template <typename T>
MyMatrix<T> SelectColumnAddZero(MyMatrix<T> const &TheMat,
                                std::vector<int> const &eList) {
  size_t nbRow = TheMat.rows();
  size_t nbColRed = eList.size();
  MyMatrix<T> TheProv(nbRow, 1 + nbColRed);
  for (size_t iRow = 0; iRow < nbRow; iRow++)
    TheProv(iRow, 0) = 0;
  for (size_t iCol = 0; iCol < nbColRed; iCol++) {
    size_t jCol = eList[iCol];
    TheProv.col(1 + iCol) = TheMat.col(jCol);
  }
  return TheProv;
}

template <typename T>
MyVector<T> SelectColumnVectorAddZero(MyVector<T> const &TheV,
                                      std::vector<int> const &eList) {
  int nbColRed = eList.size();
  MyVector<T> TheProv(1 + nbColRed);
  TheProv(0) = 0;
  for (int iCol = 0; iCol < nbColRed; iCol++) {
    int jCol = eList[iCol];
    TheProv(1 + iCol) = TheV(jCol);
  }
  return TheProv;
}

template <typename T, typename Tgroup>
std::vector<int> GetNonRedundant_Equivariant(const MyMatrix<T> &EXT,
                                             const Tgroup &GRP,
                                             std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  size_t n_rows = EXT.rows();
  size_t n_cols = EXT.cols();
  Face status_ret(n_rows);
  //
  Face work(n_rows);
  for (size_t i_row = 0; i_row < n_rows; i_row++)
    work[i_row] = 1;
  vectface vf = DecomposeOrbitPoint(GRP, work);
  size_t n_orbit = vf.size();
#ifdef DEBUG_ELIMINATION_REDUNDANCY
  os << "n_rows=" << n_rows << " n_cols=" << n_cols
     << " |GRP|=" << GRP.size() << " n_orbit=" << n_orbit << "\n";
#endif
  Face status_orbit(n_orbit);
  for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++)
    status_orbit[i_orbit] = 1;
  for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
    Face e_orbit = vf[i_orbit];
#ifdef DEBUG_ELIMINATION_REDUNDANCY
    os << "i_orbit=" << i_orbit << "/" << n_orbit
       << " |e_orbit|=" << e_orbit.count() << "\n";
    os << "O =";
    for (size_t i_row = 0; i_row < n_rows; i_row++) {
      if (e_orbit[i_row] == 1) {
        int val = i_row + 1;
        os << " " << val;
      }
    }
    os << "\n";
#endif
    //
    // Selecting the relevant rows to test against
    //
    Face select(n_rows);
#ifdef DEBUG_ELIMINATION_REDUNDANCY
    size_t sel_siz = 0;
#endif
    for (size_t j_orbit = 0; j_orbit < n_orbit; j_orbit++) {
      if (i_orbit != j_orbit && status_orbit[j_orbit] == 1) {
        Face sing_orbit = vf[j_orbit];
        boost::dynamic_bitset<>::size_type i_row = sing_orbit.find_first();
        while (i_row != boost::dynamic_bitset<>::npos) {
          select[i_row] = 1;
#ifdef DEBUG_ELIMINATION_REDUNDANCY
          sel_siz++;
#endif
          i_row = sing_orbit.find_next(i_row);
        }
      }
    }
#ifdef DEBUG_ELIMINATION_REDUNDANCY
    os << "    sel_siz=" << sel_siz << "\n";
#endif
    //
    // The single vertex, stabilizer and orbit breakdown
    //
    boost::dynamic_bitset<>::size_type i_row = e_orbit.find_first();
    Tgroup Stab = GRP.Stabilizer_OnPoints(Tidx(i_row));
    vectface vf_stab = DecomposeOrbitPoint(Stab, select);
    size_t n_row_sel = vf_stab.size();
#ifdef DEBUG_ELIMINATION_REDUNDANCY
    os << "    |Stab|=" << Stab.size() << " n_row_sel=" << n_row_sel << "\n";
#endif
    MyMatrix<T> M(n_row_sel, n_cols);
    size_t i_row_sel = 0;
    for (auto &e_orb_stab : vf_stab) {
      for (size_t i_col = 0; i_col < n_cols; i_col++)
        M(i_row_sel, i_col) = 0;
      boost::dynamic_bitset<>::size_type j_row = e_orb_stab.find_first();
      while (j_row != boost::dynamic_bitset<>::npos) {
        for (size_t i_col = 0; i_col < n_cols; i_col++)
          M(i_row_sel, i_col) += EXT(j_row, i_col);
        j_row = e_orb_stab.find_next(j_row);
      }
      i_row_sel++;
    }
    //
    // The computation itself
    //
    SelectionRowCol<T> eSelect = TMat_SelectRowCol(M);
#ifdef DEBUG_ELIMINATION_REDUNDANCY
    os << "    |eSelect.ListColSelect|=" << eSelect.ListColSelect.size()
       << " n_cols=" << n_cols << "\n";
#endif
    MyVector<T> V = GetMatrixRow(EXT, i_row);
    // return true if it is redundant. False if irredundant
    auto get_status = [&]() -> bool {
      bool test1 = IsVectorInSpace(eSelect, V);
      if (!test1)
        return false;
      MyMatrix<T> M_sel = SelectColumnAddZero(M, eSelect.ListColSelect);
      MyVector<T> V_sel = SelectColumnVectorAddZero(V, eSelect.ListColSelect);
      LpSolution<T> eSol = CDD_LinearProgramming(M_sel, V_sel, os);
      if (!eSol.DualDefined || !eSol.PrimalDefined) {
        return false;
      }
      if (eSol.OptimalValue < 0)
        return false;
      return true;
    };
    bool status = get_status();
#ifdef DEBUG_ELIMINATION_REDUNDANCY
    os << "    status=" << status << "\n";
#endif
    if (!status) {
      boost::dynamic_bitset<>::size_type j_row = e_orbit.find_first();
      while (j_row != boost::dynamic_bitset<>::npos) {
        status_ret[j_row] = 1;
        j_row = e_orbit.find_next(j_row);
      }
    } else {
      status_orbit[i_orbit] = 0;
    }
  }
  std::vector<int> ListIrred;
  for (size_t i_row = 0; i_row < n_rows; i_row++) {
    if (status_ret[i_row] == 1) {
      ListIrred.push_back(i_row);
    }
  }
  return ListIrred;
}

// clang-format off
#endif  // SRC_POLY_POLY_REDUNDANCYELIMINATION_H_
// clang-format on

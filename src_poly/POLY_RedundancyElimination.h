#ifndef INCLUDE_REDUNDANCY_CHECK
#define INCLUDE_REDUNDANCY_CHECK

#include "POLY_LinearProgramming.h"


template<typename T>
struct FacetizationInfo {
  MyMatrix<T> EXT;
  MyMatrix<T> BoundingFAC;
};



template<typename T>
FacetizationInfo<T> FacetizationCone(MyMatrix<T> const& EXT, MyMatrix<T> const& BoundingFac)
{
  int n_rows=EXT.rows();
  int n_cols=EXT.cols();
  MyMatrix<T> TheSum=ZeroMatrix<T>(1,n_cols);
  for (int i_row=0; i_row<n_rows; i_row++)
    for (int i_col=0; i_col<n_cols; i_col++)
      TheSum(0, i_col) += EXT(i_row, i_col);
  MyMatrix<T> EXTtot = Concatenate(TheSum, EXT);
  std::cerr << "FacetizationCone\n";
  MyMatrix<T> EXTbas = RowReduction(EXTtot);
  std::cerr << "EXTbas=\n";
  WriteMatrix(std::cerr, EXTbas);
  MyMatrix<T> eInvMat = Inverse(EXTbas);
  MyMatrix<T> EXT_ret = EXT * eInvMat;
  //
  MyMatrix<T> BoundingFac_ret = BoundingFac * TransposedMat(EXTbas);
  return {EXT_ret, BoundingFac_ret};
}




#define DEBUG_REDUND

template<typename T>
std::vector<int> EliminationByRedundance_HitAndRun(MyMatrix<T> const& EXT)
{
  int n_rows=EXT.rows();
  int n_cols=EXT.cols();
  std::cerr << "n_rows=" << n_rows << " n_cols=" << n_cols << "\n";
  MyMatrix<T> BoundingFac(0, n_cols);
  FacetizationInfo<T> RegCone = FacetizationCone(EXT, BoundingFac);
  MyMatrix<T> EXT_work = RegCone.EXT;
  std::cerr << "EXT_work=\n";
  WriteMatrix(std::cerr, EXT_work);
  // return true if it is redundant. False if irredundant
  // If true, also returns a list of indices showing up positively in the decomposition
  auto GetRedundancyInfo=[&](std::vector<int> const& ListIRow, int const& iRow) -> std::pair<bool,std::vector<int>> {
    std::cerr << "ListIRow =";
    for (auto& eVal : ListIRow)
      std::cerr << " " << eVal;
    std::cerr << " iRow=" << iRow << "\n";
    MyMatrix<T> EXT_sel = SelectRow(EXT, ListIRow);
    MyVector<T> eRow = GetMatrixRow(EXT, iRow);
    std::cerr << "-------------------\n";
    std::cerr << "H-representation\n";
    std::cerr << "begin\n";
    std::cerr << " " << EXT_sel.rows() << " " << EXT_sel.cols() << " integer\n";
    for (int i_r=0; i_r<EXT_sel.rows(); i_r++) {
      for (int i_c=0; i_c<EXT_sel.cols(); i_c++)
        std::cerr << " " << EXT_sel(i_r, i_c);
      std::cerr << "\n";
    }
    std::cerr << "end\n";
    std::cerr << "minimize\n";
    for (int i_c=0; i_c<EXT_sel.cols(); i_c++)
      std::cerr << " " << eRow(i_c);
    std::cerr << "\n";
    std::cerr << "-------------------\n";
    std::cerr << "Before call to CDD_LinearProgramming\n";
    LpSolution<T> eSol = CDD_LinearProgramming(EXT_sel, eRow);
    std::cerr << " After call to CDD_LinearProgramming\n";
    if (!eSol.DualDefined || !eSol.PrimalDefined) {
      return {false, {}};
    }
    if (eSol.OptimalValue < 0)
      return {false, {}};
    std::cerr << "OptimalValue=" << eSol.OptimalValue << "\n";
    std::vector<int> ListIdx;
    int n_rows_ineq=ListIRow.size();
    for (int i_row=0; i_row<n_rows_ineq; i_row++) {
      T e_val = -eSol.DualSolution(i_row);
      if (e_val < 0) {
        std::cerr << "The coefficient should be non-negative\n";
        for (int j_row=0; j_row<n_rows_ineq; j_row++) {
          std::cerr << " " << eSol.DualSolution(j_row);
        }
        std::cerr << "\n";
        throw TerminalException{1};
      }
      if (e_val > 0)
        ListIdx.push_back(ListIRow[i_row]);
    }
    std::cerr << "Exiting GetRedundancyInfo\n";
    return {true, ListIdx};
  };
  //
  // The variables showing our current state
  //
  // -1: Still unknown 0: redundant 1: irredundant
  std::vector<int> RedundancyStatus(n_rows,-1);
  std::set<int> ListKnownIrred;
  int nbFoundIrred=0;
  RankTool<T> tool(n_cols);
  auto SetIRowIrredundant=[&](int const& i_row) -> void {
      RedundancyStatus[i_row] = 1;
      ListKnownIrred.insert(i_row);
      nbFoundIrred++;
      MyVector<T> V = GetMatrixRow(EXT, i_row);
      tool.insert_if_indep(V);
  };
  //
  // Now computing one interior point.
  //
  MyVector<T> eVectInterior = GetSpaceInteriorPoint_Basic(EXT_work);
  std::cerr << "We have eVectInterior\n";
  MyVector<T> ListScalInterior = EXT_work * eVectInterior;
  for (int i_row=0; i_row<n_rows; i_row++) {
    T eScal = ListScalInterior(i_row);
    if (eScal <= 0) {
      std::cerr << "Inconsistent scalar product\n";
      throw TerminalException{1};
    }
  }
  //
  // Now running the hit and run operations.
  //
  int N=10;
  auto GetRandomVector=[&]() -> MyVector<T> {
    MyVector<T> eVect(n_cols);
    for (int i_col=0; i_col<n_cols; i_col++) {
      int val = -N + rand() % (2*N+1);
      eVect(i_col) = val;
    }
    return eVect;
  };
  // Return the found facet if successful. Returns -1 otherwise.
  auto GetRandomOutsideVector_and_HitAndRun=[&](MyVector<T> const& eVect) -> int {
    MyVector<T> ListScal = EXT_work * eVect;
    auto HasOneViolatedFacet=[&](int const& h) -> bool {
      for (int i_row=0; i_row<n_rows; i_row++) {
        T eScal = ListScal(i_row) - h * ListScalInterior(i_row);
        if (eScal < 0) {
          return true;
        }
      }
      return false;
    };
    auto GetSmallestValue=[&](int const& h) -> int {
      bool FoundUpper=false;
      T CurrUpper = -1; // specific value that will stay unused
      int nbMatch = 0;
      int idxFound = -1;
      for (int i_row=0; i_row<n_rows; i_row++) {
        T eScalDir = ListScal(i_row) - h * ListScalInterior(i_row);
        T eScal = ListScalInterior(i_row);
        if (eScalDir < 0) {
          // eScal + eScalDir * alpha = 0
          T UpperValue = - eScal / eScalDir;
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
    while(true) {
      bool test = HasOneViolatedFacet(h);
      std::cerr << "h=" << h << "\n";
      if (test)
        return GetSmallestValue(h);
      h++;
    }
  };
  int nbRuns=0;
  while(true) {
    std::cerr << "nbRuns=" << nbRuns << " nbFoundIrred=" << nbFoundIrred << "\n";
    MyVector<T> eVect = GetRandomVector();
    int idxFound = GetRandomOutsideVector_and_HitAndRun(eVect);
    if (idxFound != -1) {
      if (RedundancyStatus[idxFound] == -1)
        SetIRowIrredundant(idxFound);
    }
    nbRuns++;

    // Updating indexes
    // If we found less than n_cols then we know more can be found
    if (nbFoundIrred >= n_cols) {
      if (nbRuns > 100 * nbFoundIrred) {
        break;
      }
    }
  }
  std::cerr << "nbFoundIrred=" << nbFoundIrred << " nbRuns=" << nbRuns << "\n";
  std::cerr << "ListKnownIrred =";
  for (auto & eVal : ListKnownIrred)
    std::cerr << " " << eVal;
  std::cerr << "\n";
  //
  // Above we have determined an initial list of facets.
  // Now we use it to do the full enumeration
  //
  // Determine the status with absolute certainty.
  // ---If false then it is irredundant.
  // ---If true then it is redundant and retuns a list of candidates
  //    for finding a new facet.
  auto DetermineStatusSurely=[&](int const& i_row) -> std::pair<bool,std::vector<int>> {
    std::vector<int> ListIRow;
    for (int j_row=0; j_row<n_rows; j_row++)
      if (RedundancyStatus[j_row] != 0 && j_row != i_row)
        ListIRow.push_back(j_row);
    std::pair<bool,std::vector<int>> ePair = GetRedundancyInfo(ListIRow, i_row);
    if (!ePair.first)
      return {false, {}};
    std::vector<int> LIdx;
    for (auto & eIdx : ePair.second) {
      if (RedundancyStatus[eIdx] == -1)
        LIdx.push_back(eIdx);
    }
    return {true, LIdx};
  };
  // If it is true then it is redundant.
  // If this function return false then the facet is redundant or not.
  auto FastStatusDetermination=[&](int const& i_row) -> bool {
    MyVector<T> V = GetMatrixRow(EXT, i_row);
    if (!tool.is_in_span(V))
      return false; // We cannot conclude because the vector does not belong to the span.
    std::vector<int> ListIRow;
    for (int j_row=0; j_row<n_rows; j_row++)
      if (RedundancyStatus[j_row] == 1)
        ListIRow.push_back(j_row);
    return GetRedundancyInfo(ListIRow, i_row).first;
  };
  auto FindOneMoreFacet=[&](int const& i_start) -> void {
    std::vector<int> WorkLIdx{i_start};
    while(true) {
      std::cerr << "WorkLIdx =";
      for (auto& eVal : WorkLIdx)
        std::cerr << " " << eVal;
      std::cerr << "\n";
      std::set<int> NewCand;
      for (auto& eIdx : WorkLIdx) {
        std::pair<bool,std::vector<int>> ePair = DetermineStatusSurely(eIdx);
        std::cerr << "ePair=" << ePair.first << " V =";
        for (auto& fIdx : ePair.second)
          std::cerr << " " << fIdx;
        std::cerr << "\n";
        if (ePair.first) { // The facet is redundant
          RedundancyStatus[eIdx] = 0;
          for (auto& fIdx : ePair.second)
            NewCand.insert(fIdx);
        } else { // The facet is irredundant. Mission accomplished.
          SetIRowIrredundant(eIdx);
          return;
        }
      }
      WorkLIdx.clear();
      for (auto& eIdx : NewCand)
        if (RedundancyStatus[eIdx] == -1) // Only those unconcluded need to be considered.
          WorkLIdx.push_back(eIdx);
#ifdef DEBUG_REDUND
      if (WorkLIdx.size() == 0) {
        std::cerr << "WorkLIdx is empty. Not allowed\n";
        throw TerminalException{1};
      }
#endif
    }
  };
  auto ProcessOnePoint=[&](int const& i_row) -> void {
    bool test = FastStatusDetermination(i_row);
    std::cerr << "FastStatusDetermination : i_row=" << i_row << " test=" << test << "\n";
    if (test) { // The heuristic works. Facet is redundant. We do conclude.
      RedundancyStatus[i_row] = 0;
      return;
    }
    // If the heuristic fails then it means that we miss one irreducible facet and we should
    // find at least one
    FindOneMoreFacet(i_row);
  };
  // Now processing all the data
  for (int i_row=0; i_row<n_rows; i_row++) {
    if (RedundancyStatus[i_row] == -1) {
      std::cerr << "Before ProcessOnePoint\n";
      ProcessOnePoint(i_row);
      std::cerr << " After ProcessOnePoint\n";
    }
  }
  for (int i_row=0; i_row<n_rows; i_row++)
    if (RedundancyStatus[i_row] == -1) {
      std::cerr << "The algorithm failed to treat all the points\n";
      throw TerminalException{1};
    }
  std::vector<int> ListKnownIrred_vect;
  for (auto & eVal : ListKnownIrred)
    ListKnownIrred_vect.push_back(eVal);
  std::cerr << "|ListKnownIrred|=" << ListKnownIrred_vect.size() << "\n";
  return ListKnownIrred_vect;
}












#endif

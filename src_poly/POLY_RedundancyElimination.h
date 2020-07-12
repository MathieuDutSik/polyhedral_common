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
  MyMatrix<T> EXTbas = RowReduction(EXTtot);
  MyMatrix<T> eInvMat = Inverse(EXTbas);
  MyMatrix<T> EXT_ret = EXT * eInvMat;
  //
  MyMatrix<T> BoundingFac_ret = BoundingFac * TransposedMat(EXTbas);
  return {EXT_ret, BoundingFac_ret};
}






template<typename T>
std::vector<int> EliminationByRedundance_HitAndRun(MyMatrix<T> const& EXT)
{
  int n_rows=EXT.rows();
  int n_cols=EXT.cols();
  MyMatrix<T> BoundingFac(0, n_cols);
  FacetizationInfo<T> RegCone = FacetizationCone(EXT, BoundingFac);
  MyMatrix<T> EXT_work = RegCone.EXT;
  // return true if it is redundant. False otherwise.
  // If true, also returns a list of indices showing up positively in the decomposition
  auto GetRedundancyInfo=[&](std::vector<int> const& ListIRow, int const& iRow) -> std::pair<bool,std::vector<int>> {
    MyMatrix<T> EXT_sel = SelectRow(EXT, ListIRow);
    MyVector<T> eRow = GetMatrixRow(EXT, iRow);
    LpSolution<T> eSol = CDD_LinearProgramming(EXT_sel, eRow);
    if (!eSol.DualDefined || !eSol.PrimalDefined) {
      return {false,{}};
    }
    if (eSol.OptimalValue < 0)
      return {false,{}};
    std::vector<int> ListIdx;
    int n_rows=ListIRow.size();
    for (int i_row=0; i_row<n_rows; i_row++) {
      T e_val = eSol.DualSolution(i_row);
      if (e_val < 0) {
        std::cerr << "The coefficient should be non-negative\n";
        throw TerminalException{1};
      }
      if (e_val > 0)
        ListIdx.push_back(ListIRow[i_row]);
    }
    return {true,ListIdx};
  };
  //
  // The variables showing our current state
  //
  // -1: Still unknown 0: redundant 1: irredundant
  std::vector<int> RedundancyStatus(n_rows,-1);
  std::vector<int> ListKnownIrred;
  //
  // Now computing one interior point.
  //
  MyVector<T> eVectInterior = GetSpaceInteriorPoint_Basic(EXT_work);
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
      if (test)
        return GetSmallestValue(h);
      h++;
    }
  };
  int nbFoundIrred=0;
  int nbRuns=0;
  while(true) {
    MyVector<T> eVect = GetRandomVector();
    int idxFound = GetRandomOutsideVector_and_HitAndRun(eVect);
    if (idxFound != -1) {
      if (RedundancyStatus[idxFound] == -1) {
        RedundancyStatus[idxFound] = 1;
        ListKnownIrred.push_back(idxFound);
        nbFoundIrred++;
      }
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
      if (RedundancyStatus[j_row] != 0)
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
  // If this function return false then the facet IS irredundant
  // If it is true then it may or may not be.
  auto FastStatusDetermination=[&](int const& i_row) -> bool {
    std::vector<int> ListIRow;
    for (int j_row=0; j_row<n_rows; j_row++)
      if (RedundancyStatus[j_row] == 1)
        ListIRow.push_back(j_row);
    return GetRedundancyInfo(ListIRow, i_row).first;
  };
  auto FindOneMoreFacet=[&](std::vector<int> const& LIdx) -> void {
    std::vector<int> WorkLIdx = LIdx;
    while(true) {
      std::set<int> NewCand;
      for (auto& eIdx : WorkLIdx) {
        std::pair<bool,std::vector<int>> ePair = DetermineStatusSurely(eIdx);
        if (ePair.first) { // The facet is redundant
          RedundancyStatus[eIdx] = 0;
          for (auto& fIdx : ePair.second)
            NewCand.insert(fIdx);
        } else { // The facet is irredundant. End of story
          RedundancyStatus[eIdx] = 1;
          ListKnownIrred.push_back(eIdx);
          return;
        }
      }
      WorkLIdx.clear();
      for (auto& eIdx : NewCand)
        if (RedundancyStatus[eIdx] == -1)
          WorkLIdx.push_back(eIdx);
    }
  };
  auto ProcessOnePoint=[&](int const& i_row) -> void {
    bool test = FastStatusDetermination(i_row);
    if (!test) { // The heuristic works. We conclude
      RedundancyStatus[i_row] = 0;
      return;
    }
    // If the heuristic fails then it means that we miss one irreducible facet and we should
    // find at least one
    FindOneMoreFacet({i_row});
  };
  // Now processing all the data
  for (int i_row=0; i_row<n_rows; i_row++) {
    if (RedundancyStatus[i_row] == -1)
      ProcessOnePoint(i_row);
  }
  for (int i_row=0; i_row<n_rows; i_row++)
    if (RedundancyStatus[i_row] == -1) {
      std::cerr << "The algorithm failed to treat all the points\n";
      throw TerminalException{1};
    }
  return ListKnownIrred;
}












#endif

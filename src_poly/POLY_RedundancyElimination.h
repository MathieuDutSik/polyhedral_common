#ifndef INCLUDE_REDUNDANCY_CHECK
#define INCLUDE_REDUNDANCY_CHECK

#include "POLY_LinearProgramming.h"
#include "GRP_GroupFct.h"
#include "POLY_cddlib.h"


// Fairly expensive function. But useful for debugging
template<typename T>
std::vector<int> Kernel_GetNonRedundant_CDD(const MyMatrix<T>& M)
{
  MyMatrix<T> Mred = ColumnReduction(M);
  MyMatrix<T> EXT = cdd::DualDescription(Mred);
  size_t n_row = Mred.rows();
  size_t n_vert = EXT.rows();
  size_t n_col = Mred.cols();
  std::vector<int> TheSel;
  for (size_t i_row=0; i_row<n_row; i_row++) {
    std::vector<int> eIncd;
    for (size_t i_vert=0; i_vert<n_vert; i_vert++) {
      T scal=0;
      for (size_t i_col=0; i_col<n_col; i_col++)
        scal += EXT(i_vert,i_col) * Mred(i_row,i_col);
      if (scal == 0)
        eIncd.push_back(i_vert);
    }
    MyMatrix<T> EXT_face = SelectRow(EXT, eIncd);
    if (RankMat(EXT_face) == int(n_col-1))
      TheSel.push_back(i_row);
  }
  return TheSel;
}

template<typename T>
MyMatrix<T> GetNonRedundant_CDD(const MyMatrix<T>& M)
{
  return SelectRow(M, Kernel_GetNonRedundant_CDD(M));
}







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
  //  std::cerr << "FacetizationCone\n";
  MyMatrix<T> EXTbas = RowReduction(EXTtot);
  //  std::cerr << "EXTbas=\n";
  //  WriteMatrix(std::cerr, EXTbas);
  MyMatrix<T> eInvMat = Inverse(EXTbas);
  MyMatrix<T> EXT_ret = EXT * eInvMat;
  //
  MyMatrix<T> BoundingFac_ret = BoundingFac * TransposedMat(EXTbas);
  return {EXT_ret, BoundingFac_ret};
}




#undef DEBUG_REDUND
#undef UNSET_OUTPUT

template<typename T>
std::vector<int> EliminationByRedundance_HitAndRun(MyMatrix<T> const& EXT)
{
  int n_rows=EXT.rows();
  int n_cols=EXT.cols();
  //  std::cerr << "n_rows=" << n_rows << " n_cols=" << n_cols << "\n";
  MyMatrix<T> BoundingFac(0, n_cols);
  FacetizationInfo<T> RegCone = FacetizationCone(EXT, BoundingFac);
  MyMatrix<T> EXT_work = RegCone.EXT;
  //  std::cerr << "EXT_work=\n";
  //  WriteMatrix(std::cerr, EXT_work);
  // return true if it is redundant. False if irredundant
  // If true, also returns a list of indices showing up positively in the decomposition
  auto GetRedundancyInfo=[&](std::vector<int> const& ListIRow, int const& iRow) -> std::pair<bool,std::vector<int>> {
#ifdef UNSET_OUTPUT
    std::cerr << "ListIRow =";
    for (auto& eVal : ListIRow)
      std::cerr << " " << eVal;
    std::cerr << " iRow=" << iRow << "\n";
#endif
    //    std::cerr << "|ListIRow|=" << ListIRow.size() << " iRow=" << iRow << "\n";
    MyMatrix<T> EXT_sel = SelectRow(EXT, ListIRow);
    MyVector<T> eRow = GetMatrixRow(EXT, iRow);
#ifdef UNSET_OUTPUT
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
#endif
    //    std::cerr << "Before call to CDD_LinearProgramming\n";
    LpSolution<T> eSol = CDD_LinearProgramming(EXT_sel, eRow);
    //    std::cerr << " After call to CDD_LinearProgramming\n";
    if (!eSol.DualDefined || !eSol.PrimalDefined) {
      return {false, {}};
    }
    if (eSol.OptimalValue < 0)
      return {false, {}};
    std::vector<int> ListIdx;
    int n_rows_ineq=ListIRow.size();
    //    std::cerr << "n_rows_ineq=" << n_rows_ineq << " OptimalValue=" << eSol.OptimalValue << "\n";
    for (int i_row=0; i_row<n_rows_ineq; i_row++) {
      T e_val = -eSol.DualSolution(i_row);
#ifdef UNSET_OUTPUT
      if (e_val < 0) {
        std::cerr << "The coefficient should be non-negative\n";
        for (int j_row=0; j_row<n_rows_ineq; j_row++) {
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
  std::vector<int> RedundancyStatus(n_rows,-1);
  std::vector<int> ListKnownIrred;
  int nbFoundIrred=0;
  RankTool<T> tool(n_cols);
  auto SetIRowIrredundant=[&](int const& i_row) -> void {
    RedundancyStatus[i_row] = 1;
    ListKnownIrred.push_back(i_row);
    nbFoundIrred++;
    MyVector<T> V = GetMatrixRow(EXT, i_row);
    tool.insert_if_indep(V);
  };
  //
  // Now computing one interior point.
  //
  MyVector<T> eVectInterior = GetSpaceInteriorPoint_Basic(EXT_work);
  MyVector<T> ListScalInterior = EXT_work * eVectInterior;
#ifdef DEBUG_REDUND
  for (int i_row=0; i_row<n_rows; i_row++) {
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
  int N=10;
  MyVector<T> eVect(n_cols);
  auto SetRandomVector=[&]() -> void {
    for (int i_col=0; i_col<n_cols; i_col++) {
      int val = -N + rand() % (2*N+1);
      eVect(i_col) = val;
    }
  };
  // Return the found facet if successful. Returns -1 otherwise.
  MyVector<T> ListScal;
  auto GetRandomOutsideVector_and_HitAndRun=[&]() -> int {
    ListScal.noalias() = EXT_work * eVect;
    auto HasOneViolatedFacet=[&](int const& h) -> bool {
      for (int i_row=0; i_row<n_rows; i_row++) {
        T eScal = ListScal(i_row) - h * ListScalInterior(i_row);
        if (eScal < 0)
          return true;
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
      //      std::cerr << "h=" << h << "\n";
      if (test)
        return GetSmallestValue(h);
      h++;
    }
  };
  int nbRuns=0;
#ifdef UNSET_TRIES
  for (int i_col=0; i_col<n_cols; i_col++) {
    eVect(i_col) = 0;
  }
  for (int i_col=0; i_col<n_cols; i_col++) {
    for (int i=0; i<2; i++) {
      eVect(i_col) = -1 + 2*i;
      int idxFound = GetRandomOutsideVector_and_HitAndRun();
      if (idxFound != -1)
        if (RedundancyStatus[idxFound] == -1)
          SetIRowIrredundant(idxFound);
    }
    eVect(i_col) = 0;
  }
  for (int i_row=0; i_row<n_rows; i_row++) {
    for (int i=0; i<2; i++) {
      int alpha = -1 + 2*i;
      for (int i_col=0; i_col<n_cols; i_col++)
        eVect(i_col) = EXT_work(i_row,i_col) * alpha;
      int idxFound = GetRandomOutsideVector_and_HitAndRun();
      if (idxFound != -1)
        if (RedundancyStatus[idxFound] == -1)
          SetIRowIrredundant(idxFound);
    }
  }
#endif
  while(true) {
    //    std::cerr << "nbRuns=" << nbRuns << " nbFoundIrred=" << nbFoundIrred << "\n";
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
  //  std::cerr << "nbFoundIrred=" << nbFoundIrred << " nbRuns=" << nbRuns << "\n";
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
    for (auto & eIdx : ePair.second)
      if (RedundancyStatus[eIdx] == -1)
        LIdx.push_back(eIdx);
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
#ifdef UNSET_OUTPUT
      std::cerr << "WorkLIdx =";
      for (auto& eVal : WorkLIdx)
        std::cerr << " " << eVal;
      std::cerr << "\n";
#endif
      std::set<int> NewCand;
      for (auto& eIdx : WorkLIdx) {
        std::pair<bool,std::vector<int>> ePair = DetermineStatusSurely(eIdx);
#ifdef UNSET_OUTPUT
        std::cerr << "ePair=" << ePair.first << " V =";
        for (auto& fIdx : ePair.second)
          std::cerr << " " << fIdx;
        std::cerr << "\n";
#endif
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
    //    std::cerr << "FastStatusDetermination : i_row=" << i_row << " test=" << test << "\n";
    if (test) { // The heuristic works. Facet is redundant. We do conclude.
      RedundancyStatus[i_row] = 0;
      return;
    }
    // If the heuristic fails then it means that we miss one irreducible facet and we should
    // find at least one
    FindOneMoreFacet(i_row);
  };
  // Now processing all the data
  for (int i_row=0; i_row<n_rows; i_row++)
    if (RedundancyStatus[i_row] == -1)
      ProcessOnePoint(i_row);
#ifdef DEBUG_REDUND
  for (int i_row=0; i_row<n_rows; i_row++)
    if (RedundancyStatus[i_row] == -1) {
      std::cerr << "The algorithm failed to treat all the points\n";
      throw TerminalException{1};
    }
#endif
  std::sort(ListKnownIrred.begin(), ListKnownIrred.end());
  //  std::cerr << "|ListKnownIrred|=" << ListKnownIrred.size() << "\n";
  return ListKnownIrred;
}









template<typename T,typename Tgroup>
Face GetNonRedundant_Equivariant(const MyMatrix<T>& EXT, const Tgroup& GRP)
{
  using Telt=typename Tgroup::Telt;
  using Tidx=typename Telt::Tidx;
  size_t n_rows=EXT.rows();
  size_t n_cols=EXT.cols();
  Face status_ret(n_rows);
  //
  Face work(n_rows);
  for (size_t i_row=0; i_row<n_rows; i_row++)
    work[i_row] = 1;
  vectface vf = DecomposeOrbitPoint(GRP, work);
  size_t n_orbit = vf.size();
  Face status_orbit(n_orbit);
  for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++)
    status_orbit[i_orbit] = 1;
  for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++) {
    Face e_orbit = vf[i_orbit];
    //
    // Selecting the relevant rows to test against
    //
    Face select(n_rows);
    for (size_t j_orbit=0; j_orbit<n_orbit; j_orbit++) {
      if (i_orbit != j_orbit && status_orbit[j_orbit] == 1) {
        Face sing_orbit = vf[j_orbit];
        boost::dynamic_bitset<>::size_type i_row = sing_orbit.find_first();
        while (i_row != boost::dynamic_bitset<>::npos) {
          select[i_row] = 1;
          i_row = sing_orbit.find_next(i_row);
        }
      }
    }
    //
    // The single vertex, stabilizer and orbit breakdown
    //
    boost::dynamic_bitset<>::size_type i_row = e_orbit.find_first();
    Tgroup Stab = GRP.Stabilizer_OnPoints(Tidx(i_row));
    vectface vf_stab = DecomposeOrbitPoint(Stab, select);
    size_t n_row_sel = vf.size();
    MyMatrix<T> M(n_row_sel, n_cols);
    size_t i_row_sel = 0;
    for (auto & e_orb_stab : vf_stab) {
      for (size_t i_col=0; i_col<n_cols; i_col++)
        M(i_row_sel, i_col) = 0;
      boost::dynamic_bitset<>::size_type j_row = e_orb_stab.find_first();
      while (j_row != boost::dynamic_bitset<>::npos) {
        for (size_t i_col=0; i_col<n_cols; i_col++)
          M(i_row_sel, i_col) += EXT(j_row,i_col);
        j_row = e_orb_stab.find_next(j_row);
      }
      i_row_sel++;
    }
    //
    // The computation itself
    //
    SelectionRowCol<T> eSelect = TMat_SelectRowCol(M);
    MyVector<T> V = GetMatrixRow(EXT, i_row);
    // return true if it is redundant. False if irredundant
    auto get_status=[&]() -> bool {
      bool test1 = IsVectorInSpace(eSelect, V);
      if (!test1)
        return false;
      MyMatrix<T> M_sel = SelectColumn(M, eSelect.ListColSelect);
      MyVector<T> V_sel = SelectColumnVector(V, eSelect.ListColSelect);
      LpSolution<T> eSol = CDD_LinearProgramming(M_sel, V_sel);
      if (!eSol.DualDefined || !eSol.PrimalDefined) {
        return false;
      }
      if (eSol.OptimalValue < 0)
        return false;
      return true;
    };
    bool status = get_status();
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
  return status_ret;
}







#endif

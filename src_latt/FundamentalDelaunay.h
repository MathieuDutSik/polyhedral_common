#ifndef DELAUNAY_FUNCTION_INCLUDE
#define DELAUNAY_FUNCTION_INCLUDE

#include "CVP_NiemeierAlgo.h"
#include "POLY_LinearProgramming.h"
#include "ShortestUniversal.h"
#include "Temp_Positivity.h"
#include <string>
#include <vector>

template <typename T, typename Tint>
resultCVP<T, Tint> CVPVallentinProgram(MyMatrix<T> const &GramMat,
                                       MyVector<T> const &eV,
                                       std::string const &NameMeth) {
  bool DoCheck = false;
  if (DoCheck) {
    resultCVP<T, Tint> res1 = CVPVallentinProgram_exact<T, Tint>(GramMat, eV);
    //  resultCVP<T> res2=CVPVallentinProgram_double(GramMat, eV);
    resultCVP<T, Tint> res2 = CVP_N23_24A1<T, Tint>(eV);
    if (res1 != res2) {
      std::cerr << "res1.TheNorm=" << res1.TheNorm << "\n";
      std::cerr << "res2.TheNorm=" << res2.TheNorm << "\n";
      std::cerr << "|res1.ListVect|=" << res1.ListVect.rows()
                << " |res2.ListVect|=" << res2.ListVect.rows() << "\n";
      std::cerr << "res1.ListVect=\n";
      WriteMatrixGAP(std::cerr, res1.ListVect);
      std::cerr << "res2.ListVect=\n";
      WriteMatrixGAP(std::cerr, res2.ListVect);
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
    std::cerr << "All correct\n";
    return res1;
  }
  //
  if (NameMeth == "SVexact")
    return CVPVallentinProgram_exact<T, Tint>(GramMat, eV);
  //
  if (NameMeth == "SVdouble")
    return CVPVallentinProgram_double<T, Tint>(GramMat, eV);
  //
  if (NameMeth == "CVP_N23_24A1")
    return CVP_N23_24A1<T, Tint>(eV);
  //
  std::cerr << "No matching method found\n";
  throw TerminalException{1};
}

template <typename T>
MyVector<T> FuncRandomDirection(int const &n, int const &siz) {
  MyVector<T> eVect(n);
  int siz2 = 2 * siz + 1;
  for (int i = 0; i < n; i++) {
    int eVal = rand() % siz2;
    eVect(i) = eVal - siz;
  }
  return eVect;
}

template <typename T, typename Tint>
MyMatrix<Tint> FindDelaunayPolytope(MyMatrix<T> const &GramMat,
                                    std::string const &CVPmethod,
                                    std::ostream &os) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int dim = GramMat.rows();
  std::vector<MyVector<T>> ListRelevantPoint;
  for (int i = 0; i < dim; i++) {
    MyVector<T> V1 = ZeroVector<T>(dim);
    V1(i) = 1;
    ListRelevantPoint.push_back(V1);
    MyVector<T> V2 = ZeroVector<T>(dim);
    V2(i) = -1;
    ListRelevantPoint.push_back(V2);
  }
  auto DefiningInequality = [&](MyVector<T> const &eVect) -> MyVector<T> {
    MyVector<T> eIneq(1 + dim);
    T eSum = EvaluationQuadForm<T, T>(GramMat, eVect);
    eIneq(0) = eSum;
    for (int iCol = 0; iCol < dim; iCol++) {
      T eSumB = 0;
      for (int iRow = 0; iRow < dim; iRow++)
        eSumB += eVect(iRow) * GramMat(iRow, iCol);
      eIneq(iCol + 1) = -2 * eSumB;
    }
    return eIneq;
  };
  MyVector<T> TheRandomDirection = FuncRandomDirection<T>(dim + 1, 10);
  while (true) {
    int nbVect = ListRelevantPoint.size();
    MyMatrix<T> ListIneq(nbVect, dim + 1);
    for (int iVect = 0; iVect < nbVect; iVect++) {
      MyVector<T> fVect = DefiningInequality(ListRelevantPoint[iVect]);
      for (int i = 0; i <= dim; i++)
        ListIneq(iVect, i) = fVect(i);
    }
    LpSolution<T> eSol = CDD_LinearProgramming(ListIneq, TheRandomDirection);
    assert(eSol.PrimalDefined);
    MyVector<T> eVect = eSol.DirectSolution;
    T TheNorm = EvaluationQuadForm<T, T>(GramMat, eVect);
    //    std::cerr << "Before CVPVallentinProgram\n";
    resultCVP<T, Tint> TheCVP =
        CVPVallentinProgram<T, Tint>(GramMat, eVect, CVPmethod);
    //    std::cerr << "After CVPVallentinProgram\n";
    int nbVectTot = TheCVP.ListVect.rows();
    if (TheCVP.TheNorm == TheNorm) {
      MyMatrix<Tint> RetEXT(nbVectTot, dim + 1);
      for (int iVect = 0; iVect < nbVectTot; iVect++) {
        RetEXT(iVect, 0) = 1;
        for (int i = 0; i < dim; i++)
          RetEXT(iVect, 1 + i) = TheCVP.ListVect(iVect, i);
      }
      return RetEXT;
    } else {
      for (int iVect = 0; iVect < nbVectTot; iVect++) {
        MyVector<Tint> fVect = GetMatrixRow(TheCVP.ListVect, iVect);
        MyVector<T> fVect_T = UniversalVectorConversion<T, Tint>(fVect);
        ListRelevantPoint.push_back(fVect_T);
      }
    }
    os << "|ListRelevantPoint|=" << ListRelevantPoint.size() << "\n";
  }
}

template <typename T> struct CP {
  T SquareRadius;
  MyVector<T> eCent;
};

template <typename T>
CP<T> CenterRadiusDelaunayPolytopeGeneral(MyMatrix<T> const &GramMat,
                                          MyMatrix<T> const &EXT) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int n = GramMat.rows();
  /*  std::cerr << "GramMat=\n";
  WriteMatrixGAP(std::cerr, GramMat);
  std::cerr << "EXT=\n";
  WriteMatrixGAP(std::cerr, EXT);*/

  SelectionRowCol<T> eSelect = TMat_SelectRowCol(EXT);
  int dimNSP = eSelect.NSP.rows();
  //  std::cerr << "dimNSP=" << dimNSP << "\n";
  int nbEqua = dimNSP + EXT.rows() - 1;
  //  std::cerr << "nbEqua=" << nbEqua << "\n";
  int nbVert = EXT.rows();
  MyMatrix<T> ListEquation(n, nbEqua);
  MyVector<T> ListB(nbEqua);
  int idx = 0;
  for (int iVert = 1; iVert < nbVert; iVert++) {
    MyVector<T> eV(n);
    MyVector<T> fV(n);
    for (int i = 0; i < n; i++) {
      eV(i) = EXT(0, i + 1);
      fV(i) = EXT(iVert, i + 1);
    }
    T Sum = EvaluationQuadForm<T, T>(GramMat, eV) -
            EvaluationQuadForm<T, T>(GramMat, fV);
    ListB(idx) = Sum;
    for (int i = 0; i < n; i++) {
      T vSum = 0;
      for (int j = 0; j < n; j++)
        vSum += GramMat(i, j) * (eV(j) - fV(j));
      ListEquation(i, idx) = 2 * vSum;
    }
    idx++;
  }
  for (int iNSP = 0; iNSP < dimNSP; iNSP++) {
    ListB(idx) = eSelect.NSP(idx, 0);
    for (int i = 0; i < n; i++)
      ListEquation(i, idx) = -eSelect.NSP(idx, i + 1);
    idx++;
  }
  SolMatResult<T> Solu = SolutionMat(ListEquation, ListB);
  assert(Solu.result);
  MyVector<T> eCent(1 + n);
  eCent(0) = 1;
  for (int i = 0; i < n; i++)
    eCent(i + 1) = Solu.eSol(i);
  T SquareRadius;
  for (int iVert = 0; iVert < nbVert; iVert++) {
    MyVector<T> eW(n);
    for (int i = 0; i < n; i++)
      eW(i) = Solu.eSol(i) - EXT(iVert, i + 1);
    T eNorm = EvaluationQuadForm<T, T>(GramMat, eW);
    if (iVert == 0) {
      SquareRadius = eNorm;
    } else {
      if (SquareRadius != eNorm) {
        std::cerr << "SquareRadius=" << SquareRadius << "\n";
        std::cerr << "eNorm=" << eNorm << "\n";
        std::cerr << "Error in the Solutioning\n";
        throw TerminalException{1};
      }
    }
  }
  return {SquareRadius, eCent};
}

template <typename T, typename Tint>
MyMatrix<Tint>
FindAdjacentDelaunayPolytope(MyMatrix<T> const &GramMat, MyMatrix<T> const &EXT,
                             Face const &eInc, std::string const &CVPmethod) {
  //  std::cerr << "FindAdjacentDelaunayPolytope, step 1\n";
  int dim = GramMat.rows();
  //  std::cerr << "FindAdjacentDelaunayPolytope, step 1.1\n";
  MyMatrix<T> EXTred = SelectRow(EXT, eInc);
  //  std::cerr << "FindAdjacentDelaunayPolytope, step 1.2\n";
  MyMatrix<T> IndependentBasis = RowReduction(EXTred);
  //  std::cerr << "FindAdjacentDelaunayPolytope, step 1.3\n";
  MyVector<T> TheFac = FindFacetInequality(EXT, eInc);
  //  std::cerr << "FindAdjacentDelaunayPolytope, step 2\n";
  int iColFind = -1;
  for (int iCol = 0; iCol < dim; iCol++)
    if (iColFind == -1 && TheFac(1 + iCol) != 0)
      iColFind = iCol;
  MyVector<T> eVect = ZeroVector<T>(1 + dim);
  eVect(1 + iColFind) = -T_sign(TheFac(1 + iColFind));
  int jRow = eInc.find_first();
  MyVector<T> SelectedVertex(1 + dim);
  //  std::cerr << "FindAdjacentDelaunayPolytope, step 3\n";
  for (int iCol = 0; iCol <= dim; iCol++)
    SelectedVertex(iCol) = EXT(jRow, iCol) + eVect(iCol);
  auto GetCenterRadius = [&](MyVector<T> const &TheVert) -> CP<T> {
    MyMatrix<T> VSet(dim + 1, dim + 1);
    for (int i = 0; i < dim; i++)
      VSet.row(i) = IndependentBasis.row(i);
    //    std::cerr << "IndependentBasis part of VSet assigned\n";
    for (int j = 0; j <= dim; j++)
      VSet(dim, j) = TheVert(j);
    //    std::cerr << "TheVert put in VSet\n";
    return CenterRadiusDelaunayPolytopeGeneral(GramMat, VSet);
  };
  T MinRadius = GetCenterRadius(SelectedVertex).SquareRadius;
  std::vector<MyVector<T>> ListGraverOptions;
  //  std::cerr << "FindAdjacentDelaunayPolytope, step 4\n";
  for (int i = 0; i < dim; i++) {
    MyVector<T> V1 = ZeroVector<T>(dim + 1);
    V1(i + 1) = 1;
    ListGraverOptions.push_back(V1);
    MyVector<T> V2 = ZeroVector<T>(dim + 1);
    V2(i + 1) = -1;
    ListGraverOptions.push_back(V2);
  }
  //  std::cerr << "FindAdjacentDelaunayPolytope, step 5\n";
  auto fGraverUpdate = [&]() -> void {
    while (true) {
      bool IsImprovement = false;
      for (auto &fVect : ListGraverOptions) {
        MyVector<T> NewTestVert = SelectedVertex + fVect;
        T eScal = 0;
        for (int i = 0; i <= dim; i++)
          eScal += NewTestVert(i) * TheFac(i);
        if (eScal < 0) {
          T TheRadius = GetCenterRadius(NewTestVert).SquareRadius;
          if (TheRadius < MinRadius) {
            IsImprovement = true;
            SelectedVertex = NewTestVert;
            MinRadius = TheRadius;
          }
        }
      }
      if (!IsImprovement)
        break;
    }
  };
  fGraverUpdate();
  resultCVP<T, Tint> reply;
  //  std::cerr << "FindAdjacentDelaunayPolytope, step 6\n";
  while (true) {
    CP<T> eCP = GetCenterRadius(SelectedVertex);
    MyVector<T> eCenter(dim);
    for (int i = 0; i < dim; i++)
      eCenter(i) = eCP.eCent(i + 1);
    //    std::cerr << "Before CVPVallentinProgram\n";
    reply = CVPVallentinProgram<T, Tint>(GramMat, eCenter, CVPmethod);
    //    std::cerr << "After CVPVallentinProgram\n";
    if (reply.TheNorm == eCP.SquareRadius)
      break;
    for (int i = 0; i < dim; i++)
      SelectedVertex(i + 1) = reply.ListVect(0, i);
    fGraverUpdate();
  }
  //  std::cerr << "FindAdjacentDelaunayPolytope, step 7\n";
  int nbVect = reply.ListVect.rows();
  MyMatrix<Tint> RetEXT(nbVect, dim + 1);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    RetEXT(iVect, 0) = 1;
    for (int i = 0; i < dim; i++)
      RetEXT(iVect, i + 1) = reply.ListVect(iVect, i);
  }
  //  std::cerr << "FindAdjacentDelaunayPolytope, step 8\n";
  return RetEXT;
}

template <typename T>
std::vector<MyMatrix<T>> CharacterizingPair(MyMatrix<T> const &GramMat,
                                            MyVector<T> const &TheCenter) {
  int n = GramMat.rows();
  MyVector<T> TheCenterRed(n);
  for (int i = 0; i < n; i++)
    TheCenterRed(i) = TheCenter(i + 1);
  MyMatrix<T> Mat1(n + 1, n + 1);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      Mat1(i + 1, j + 1) = GramMat(i, j);
  MyMatrix<T> eVect = ProductVectorMatrix(TheCenterRed, GramMat);
  for (int i = 0; i < n; i++) {
    Mat1(0, i + 1) = -eVect(i);
    Mat1(i + 1, 0) = -eVect(i);
  }
  T eScal = 1;
  for (int i = 0; i < n; i++)
    eScal += eVect(i) * TheCenterRed(i);
  Mat1(0, 0) = eScal;
  MyMatrix<T> Mat2 = ZeroMatrix<T>(n + 1, n + 1);
  Mat2(0, 0) = 1;
  while (true) {
    bool test = IsPositiveDefinite(Mat1);
    if (test)
      break;
    Mat1(0, 0)++;
  }
  return {Mat1, Mat2};
}

#endif

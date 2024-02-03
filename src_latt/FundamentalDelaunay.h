// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_FUNDAMENTALDELAUNAY_H_
#define SRC_LATT_FUNDAMENTALDELAUNAY_H_

// clang-format off
#include "POLY_LinearProgramming.h"
#include "ShortestUniversal.h"
#include "Temp_Positivity.h"
#include <string>
#include <vector>
// clang-format on

#ifdef TIMINGS
#define TIMINGS_DELAUNAY_ENUMERATION
#endif

#ifdef DEBUG
#define DEBUG_DELAUNAY_ENUMERATION
#endif

template <typename T, typename Tint>
resultCVP<T, Tint> CVPVallentinProgram(MyMatrix<T> const &GramMat,
                                       MyVector<T> const &eV,
                                       std::string const &NameMeth, std::ostream& os) {
  if (NameMeth == "SVexact")
    return CVPVallentinProgram_exact<T, Tint>(GramMat, eV, os);
  //  if (NameMeth == "SVdouble")
  //    return CVPVallentinProgram_double<T, Tint>(GramMat, eV);
  std::cerr << "No matching method found\n";
  throw TerminalException{1};
}

template <typename T>
MyVector<T> FuncRandomDirection(int const &n, int const &siz) {
  MyVector<T> eVect(n);
  int siz2 = 2 * siz + 1;
  for (int i = 0; i < n; i++) {
    int eVal = random() % siz2;
    eVect(i) = eVal - siz;
  }
  return eVect;
}

template <typename T, typename Tint>
MyMatrix<Tint> FindDelaunayPolytope(MyMatrix<T> const &GramMat,
                                    CVPSolver<T,Tint> & solver,
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
    resultCVP<T, Tint> TheCVP = solver.SingleSolver(eVect);
    int nbVectTot = TheCVP.ListVect.rows();
    if (TheCVP.TheNorm == TheNorm) {
      MyMatrix<Tint> RetEXT(nbVectTot, dim + 1);
      for (int iVect = 0; iVect < nbVectTot; iVect++) {
        RetEXT(iVect, 0) = 1;
        for (int i = 0; i < dim; i++)
          RetEXT(iVect, 1 + i) = TheCVP.ListVect(iVect, i);
      }
#ifdef DEBUG_DELAUNAY_ENUMERATION
      os << "DEL_ENUM: f_init, Creation of a Delaunay with |V|=" << RetEXT.rows() << " vertices\n";
#endif
      return RetEXT;
    } else {
      for (int iVect = 0; iVect < nbVectTot; iVect++) {
        MyVector<Tint> fVect = GetMatrixRow(TheCVP.ListVect, iVect);
        MyVector<T> fVect_T = UniversalVectorConversion<T, Tint>(fVect);
        ListRelevantPoint.push_back(fVect_T);
      }
    }
#ifdef DEBUG_DELAUNAY_ENUMERATION
    os << "|ListRelevantPoint|=" << ListRelevantPoint.size() << "\n";
#endif
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
  std::optional<MyVector<T>> opt = SolutionMat(ListEquation, ListB);
  MyVector<T> const& eSol = *opt;
  MyVector<T> eCent(1 + n);
  eCent(0) = 1;
  for (int i = 0; i < n; i++)
    eCent(i + 1) = eSol(i);
  auto get_sqr_dist=[&](int iVert) -> T {
    MyVector<T> eW(n);
    for (int i = 0; i < n; i++)
      eW(i) = eSol(i) - EXT(iVert, i + 1);
    return EvaluationQuadForm<T, T>(GramMat, eW);
  };
  T SquareRadius = get_sqr_dist(0);
#ifdef DEBUG_DELAUNAY
  for (int iVert = 1; iVert < nbVert; iVert++) {
    T eNorm = get_sqr_dist(iVert);
    if (SquareRadius != eNorm) {
      std::cerr << "SquareRadius=" << SquareRadius << "\n";
      std::cerr << "eNorm=" << eNorm << "\n";
      std::cerr << "Error in the Solutioning\n";
      throw TerminalException{1};
    }
  }
#endif
  return {std::move(SquareRadius), std::move(eCent)};
}


// This is for finding cirsumscribing spheres and their radius
//
// In the constructor we compute the point at the center of the facet
// and the direction. It is basically an affine space.
//
// G[v1 - v] = G[v2 - v]
// G[v1] - G[v2] = 2 v1 G v - 2 v2 G v
//               = (2 v1 G - 2 v2 G) v
// So, this is a linear system for which we need the solution
// Expressed as an homogeneous system we obtain
// 0 = h (G[v1] - G[v2]) + (2 v2 G - 2 v1 G) v
//
// Next we want to get the center and radius.
// The equation to resolve is
// G[ePt + t eDir - v1] = G[ePt + t eDir - v2]
// with v1 in the facet and v2 outside of it.
// That reduces to
// G[v1] - G[v2] = (2 v1 G - 2 v2 G) (ePt + t eDir)
// The linear system is then written as
// C = t D
template<typename T>
struct AdjacentDelaunayPointSolver {
private:
  int n;
  MyMatrix<T> const& GramMat;
  MyMatrix<T> const& EXT;
  Face const& eInc;
  std::vector<T> LineGramMat;
  MyVector<T> ePt;
  MyVector<T> eDir;
  MyVector<T> v1;
  T G_v1;
  MyVector<T> two_v1_G;
  std::ostream& os;
public:
  AdjacentDelaunayPointSolver(MyMatrix<T> const &_GramMat, MyMatrix<T> const &_EXT, Face const &_eInc, std::ostream& _os) : n(_GramMat.rows()), GramMat(_GramMat), EXT(_EXT), eInc(_eInc), LineGramMat(GetLineVector(GramMat)), ePt(GramMat.rows()), eDir(GramMat.rows()), v1(n), os(_os) {
    std::vector<size_t> V;
    size_t n_vert = EXT.rows();
    for (size_t i=0; i<n_vert; i++) {
      if (eInc[i] == 1) {
        V.push_back(i);
      }
    }
    size_t cnt = V.size();
    for (int i=0; i<n; i++) {
      v1(i) = EXT(V[0], i+1);
    }
    G_v1 = EvaluateLineVector(LineGramMat, v1);
    two_v1_G = 2 * GramMat * v1;
    MyMatrix<T> Equa(cnt-1, n+1);
    MyVector<T> v2(n);
    for (size_t i_vert=0; i_vert<cnt-1; i_vert++) {
      for (int i=0; i<n; i++) {
        v2(i) = EXT(V[i_vert + 1], i+1);
      }
      T G_v2 = EvaluateLineVector(LineGramMat, v2);
      MyVector<T> two_v2_G = 2 * GramMat * v2;
      Equa(i_vert, 0) = G_v1 - G_v2;
      for (int i=0; i<n; i++) {
        Equa(i_vert, i + 1) = two_v2_G(i) - two_v1_G(i);
      }
    }
    MyMatrix<T> NSP = NullspaceTrMat(Equa);
#ifdef DEBUG_DELAUNAY_ENUMERATION
    if (NSP.rows() != 2) {
      std::cerr << "Error in the computation of the kernel\n";
      throw TerminalException{1};
    }
#endif
    auto get_point_direction=[&]() -> void {
      for (int pos=0; pos<2; pos++) {
        if (NSP(pos,0) != 0) {
          for (int i=0; i<n; i++) {
            ePt(i) = NSP(pos,i+1) / NSP(pos,0);
          }
          for (int i=0; i<n; i++) {
            eDir(i) = NSP(1-pos, i+1) - ePt(i) * NSP(1-pos,0);
          }
          return;
        }
      }
      std::cerr << "Failed to find a non-zero entry\n";
      throw TerminalException{1};
    };
    get_point_direction();
  }
  // v2 should be a point outside of the plane
  CP<T> GetCenterRadius(MyVector<T> const& v2) {
    T G_v2 = EvaluateLineVector(LineGramMat, v2);
    MyVector<T> two_v1_m_v2_G = two_v1_G - 2 * GramMat * v2;
    T C_cst = G_v1 - G_v2 - two_v1_m_v2_G.dot(ePt);
    T D_cst = two_v1_m_v2_G.dot(eDir);
    T t = C_cst / D_cst;
    MyVector<T> eCent = ePt + t * eDir;
    MyVector<T> delta = eCent - v1;
    T eSqrDist = EvaluateLineVector(LineGramMat, delta);
    MyVector<T> eCentRet(n+1);
    eCentRet(0) = 1;
    for (int i=0; i<n; i++)
      eCentRet(i+1) = eCent(i);
    return {eSqrDist, eCentRet};
  }
  T GetRadius(MyVector<T> const& v2) {
    T G_v2 = EvaluateLineVector(LineGramMat, v2);
    MyVector<T> two_v1_m_v2_G = two_v1_G - 2 * GramMat * v2;
    T C_cst = G_v1 - G_v2 - two_v1_m_v2_G.dot(ePt);
    T D_cst = two_v1_m_v2_G.dot(eDir);
    T t = C_cst / D_cst;
    MyVector<T> eCent = ePt + t * eDir;
    MyVector<T> delta = eCent - v1;
    return EvaluateLineVector(LineGramMat, delta);
  }
};



template <typename T, typename Tint>
MyMatrix<Tint>
FindAdjacentDelaunayPolytope(MyMatrix<T> const &GramMat, CVPSolver<T,Tint> & solver, MyMatrix<Tint> const& ShvGraverBasis,
                             MyMatrix<T> const &EXT, Face const &eInc, std::ostream& os) {
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  MicrosecondTime time;
#endif
  int dim = GramMat.rows();
  MyVector<T> TheFac = FindFacetInequality(EXT, eInc);
  auto get_iColFind=[&]() -> int {
    for (int iCol = 0; iCol < dim; iCol++)
      if (TheFac(1 + iCol) != 0)
        return iCol;
    std::cerr << "Failed to find the matching iCol\n";
    throw TerminalException{1};
  };
  int iColFind = get_iColFind();
  T delta = -T_sign(TheFac(1 + iColFind));
  int jRow = eInc.find_first();
  MyVector<T> SelectedVertex(dim);
  for (int i=0; i<dim; i++)
    SelectedVertex(i) = EXT(jRow, i + 1);
  SelectedVertex(iColFind) += delta;
  AdjacentDelaunayPointSolver<T> adps(GramMat, EXT, eInc, os);
  T MinRadius = adps.GetRadius(SelectedVertex);
  auto fGraverUpdate = [&]() -> void {
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    MicrosecondTime time_graver;
#endif
#ifdef DEBUG_DELAUNAY_ENUMERATION
    size_t n_update = 0, n_loop = 0;
#endif
    while (true) {
      bool IsImprovement = false;
      MyVector<T> NewTestVert(dim);
      int n_graver = ShvGraverBasis.rows();
      for (int i_graver=0; i_graver<n_graver; i_graver++) {
        for (int i=0; i<dim; i++) {
          NewTestVert(i) = SelectedVertex(i) + ShvGraverBasis(i_graver, i);
        }
        T eScal = TheFac(0);
        for (int i = 0; i < dim; i++)
          eScal += NewTestVert(i) * TheFac(i + 1);
        if (eScal < 0) {
          T TheRadius = adps.GetRadius(NewTestVert);
          if (TheRadius < MinRadius) {
            IsImprovement = true;
            SelectedVertex = NewTestVert;
            MinRadius = TheRadius;
#ifdef DEBUG_DELAUNAY_ENUMERATION
            n_update++;
#endif
          }
        }
      }
#ifdef DEBUG_DELAUNAY_ENUMERATION
      n_loop++;
#endif
      if (!IsImprovement) {
#ifdef DEBUG_DELAUNAY_ENUMERATION
        os << "DEL_ENUM: n_loop=" << n_loop << " n_update=" << n_update << " |ShvGraverBasis|=" << n_graver << "\n";
#endif
#ifdef TIMINGS_DELAUNAY_ENUMERATION
        os << "|fGraverUpdate|=" << time_graver << "\n";
#endif
        break;
      }
    }
  };
  fGraverUpdate();
  auto get_reply=[&]() -> resultCVP<T, Tint> {
    while (true) {
      CP<T> eCP = adps.GetCenterRadius(SelectedVertex);
      MyVector<T> eCenter(dim);
      for (int i = 0; i < dim; i++)
        eCenter(i) = eCP.eCent(i + 1);
      resultCVP<T, Tint> reply = solver.SingleSolver(eCenter);
      if (reply.TheNorm == eCP.SquareRadius)
        return reply;
      for (int i = 0; i < dim; i++)
        SelectedVertex(i) = UniversalScalarConversion<T,Tint>(reply.ListVect(0, i));
      fGraverUpdate();
    }
  };
  resultCVP<T, Tint> reply = get_reply();
  int nbVect = reply.ListVect.rows();
  MyMatrix<Tint> RetEXT(nbVect, dim + 1);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    RetEXT(iVect, 0) = 1;
    for (int i = 0; i < dim; i++)
      RetEXT(iVect, i + 1) = reply.ListVect(iVect, i);
  }
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  os << "|FindAdjacentDelaunayPolytope|=" << time << "\n";
#endif
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

// clang-format off
#endif  // SRC_LATT_FUNDAMENTALDELAUNAY_H_
// clang-format on

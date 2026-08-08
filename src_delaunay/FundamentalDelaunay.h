// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_FUNDAMENTALDELAUNAY_H_
#define SRC_LATT_FUNDAMENTALDELAUNAY_H_

// clang-format off
#include "POLY_SamplingFacet.h"
#include "POLY_LinearProgramming.h"
#include "Shvec_exact.h"
#include "Positivity.h"
#include <limits>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef TIMINGS
#define TIMINGS_FUNDAMENTAL_DELAUNAY
#endif

#ifdef DEBUG
#define DEBUG_FUNDAMENTAL_DELAUNAY
#define DEBUG_DELAUNAY
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_FUNDAMENTAL_DELAUNAY
#endif

#ifdef DISABLE_DEBUG_FUNDAMENTAL_DELAUNAY
#undef DEBUG_FUNDAMENTAL_DELAUNAY
#endif

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

/*
  The Delaunay geometry below is written against a solver rather than
  against a fixed one: the only things it asks of Tsolver are the member
  GramMat and the method nearest_vectors, returning the closest points of
  the point set to a given point together with their common distance.

  That is the whole of what separates the lattice case from the periodic
  one, so both use this code: the lattice solver enumerates the closest
  points of Z^n, the periodic one the closest points of the periodic point
  set expressed in the coordinates scaled by its denominator, in which it
  is a union of cosets of a sublattice of Z^n. The other operation that
  has to move between points of the set, the local improvement by the
  Graver moves, is already a parameter (ShvGraverBasis) and only needs to
  be given the moves of the relevant lattice.
 */
template <typename T, typename Tint, typename Tsolver>
std::optional<MyMatrix<Tint>>
FindDelaunayPolytope_direction(Tsolver const &solver,
                               MyVector<T> const &TheRandomDirection,
                               std::ostream &os) {
  MyMatrix<T> const& GramMat = solver.GramMat;
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int dim = GramMat.rows();
  std::vector<MyVector<T>> ListIneq_vect;
  auto DefiningInequality = [&](MyVector<T> const &eVect) -> MyVector<T> {
    MyVector<T> eIneq(1 + dim);
    T eSum = EvaluationQuadForm<T, T>(GramMat, eVect);
    eIneq(0) = eSum;
    for (int iCol = 0; iCol < dim; iCol++) {
      T eSumB(0);
      for (int iRow = 0; iRow < dim; iRow++)
        AddMul(eSumB, eVect(iRow), GramMat(iRow, iCol));
      eIneq(iCol + 1) = -2 * eSumB;
    }
    return eIneq;
  };
  // The seeds are differences between points of the point set, which the
  // solver provides: writing them as the unit vectors would assume that
  // the point set is the lattice Z^n itself.
  for (auto &eVect : solver.get_seed_differences()) {
    ListIneq_vect.push_back(DefiningInequality(eVect));
  }
  while (true) {
#ifdef TIMINGS_FUNDAMENTAL_DELAUNAY
    MicrosecondTime time;
#endif
    MyMatrix<T> ListIneq = MatrixFromVectorFamily(ListIneq_vect);
#ifdef TIMINGS_FUNDAMENTAL_DELAUNAY
    os << "|DEL: Building ListIneq|=" << time << "\n";
#endif
#ifdef DEBUG_FUNDAMENTAL_DELAUNAY
    os << "DEL: nbVect=" << ListIneq.rows() << "\n";
#endif
    LpSolution<T> eSol =
        SIMPLEX_LinearProgramming(ListIneq, TheRandomDirection, os);
#ifdef TIMINGS_FUNDAMENTAL_DELAUNAY
    os << "|DEL: Computing LP eSol|=" << time << "\n";
#endif
#ifdef SANITY_CHECK_FUNDAMENTAL_DELAUNAY
    if (!eSol.DirectSolution) {
      std::cerr << "DEL: We should have a primal solution\n";
      throw TerminalException{1};
    }
#endif
    MyVector<T> const &eVect = *eSol.DirectSolution;
    if (!IsVertexOfPolyhedron(ListIneq, eVect)) {
      return {};
    }
    T TheNorm = EvaluationQuadForm<T, T>(GramMat, eVect);
    // There has been an attempt to accelerate the computation by stopping
    // when we found a vector of norm lower than TheNorm. This turn out
    // badly as the number of iterations grow immensely.
    resultCVP<T, Tint> TheCVP = solver.nearest_vectors(eVect);
#ifdef TIMINGS_FUNDAMENTAL_DELAUNAY
    os << "|DEL: Computing TheCVP|=" << time << "\n";
#endif
    int nbVectTot = TheCVP.ListVect.rows();
    if (TheCVP.TheNorm == TheNorm) {
      MyMatrix<Tint> RetEXT(nbVectTot, dim + 1);
      for (int iVect = 0; iVect < nbVectTot; iVect++) {
        RetEXT(iVect, 0) = 1;
        for (int i = 0; i < dim; i++) {
          RetEXT(iVect, 1 + i) = TheCVP.ListVect(iVect, i);
        }
      }
#ifdef TIMINGS_FUNDAMENTAL_DELAUNAY
      os << "|DEL: Processing TheCVP 1|=" << time << "\n";
#endif
#ifdef DEBUG_FUNDAMENTAL_DELAUNAY
      os << "DEL_ENUM: f_init, Creation of a Delaunay with |V|="
         << RetEXT.rows() << " vertices\n";
#endif
#ifdef SANITY_CHECK_FUNDAMENTAL_DELAUNAY
      if (RankMat(RetEXT) != dim + 1) {
        std::cerr << "DEL: FindDelaunayPolytope (initial) RetEXT rank="
                  << RankMat(RetEXT) << " expected " << (dim + 1) << "\n";
        throw TerminalException{1};
      }
#endif
      return RetEXT;
    } else {
      for (int iVect = 0; iVect < nbVectTot; iVect++) {
        MyVector<Tint> fVect = GetMatrixRow(TheCVP.ListVect, iVect);
        MyVector<T> fVect_T = UniversalVectorConversion<T, Tint>(fVect);
        MyVector<T> nVect = DefiningInequality(fVect_T);
        ListIneq_vect.push_back(nVect);
      }
    }
#ifdef TIMINGS_FUNDAMENTAL_DELAUNAY
    os << "|DEL: Processing TheCVP 2|=" << time << "\n";
#endif
#ifdef DEBUG_FUNDAMENTAL_DELAUNAY
    os << "DEL: |ListIneq_vect|=" << ListIneq_vect.size() << "\n";
#endif
  }


}




template <typename T, typename Tint, typename Tsolver>
MyMatrix<Tint> FindDelaunayPolytope(Tsolver const &solver, std::ostream &os) {
  int dim = solver.GramMat.rows();
  int N = 3;
  while (true) {
    MyVector<T> TheRandomDirection = FuncRandomDirection<T>(dim + 1, N);
    std::optional<MyMatrix<Tint>> opt = FindDelaunayPolytope_direction<T,Tint>(solver, TheRandomDirection, os);
    if (opt) {
#ifdef DEBUG_FUNDAMENTAL_DELAUNAY
      os << "DEL: FindDelaunayPolytope, find a Delaunay at N=" << N << "\n";
#endif
      return *opt;
    }
    N += 1;
  }
}

template <typename T> struct CP {
  T SquareRadius;
  MyVector<T> eCent;
};

// Squared circumscribing radius from a KNOWN center (homogeneous eCent): the
// squared distance (in GramMat) from the center to the first vertex.
template <typename T>
T SquareRadiusFromCenter(MyMatrix<T> const &GramMat, MyMatrix<T> const &EXT,
                         MyVector<T> const &eCent) {
  int n = GramMat.rows();
  MyVector<T> eW(n);
  for (int i = 0; i < n; i++)
    eW(i) = eCent(i + 1) - EXT(0, i + 1);
  return EvaluationQuadForm<T, T>(GramMat, eW);
}

// Center (homogeneous, length n+1, eCent(0) = 1) of the circumscribing sphere of
// a Delaunay polytope, general version: it determines the affine span of EXT (a
// rank/nullspace computation, TMat_SelectRowCol) and adds the corresponding
// equations, so the center is well defined even when EXT is lower-dimensional.
// Many callers only need the center; the radius is one extra quadratic form,
// added by CenterRadiusDelaunayPolytopeGeneral.
template <typename T>
MyVector<T> CenterDelaunayPolytopeGeneral(MyMatrix<T> const &GramMat,
                                          MyMatrix<T> const &EXT) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int n = GramMat.rows();
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(EXT);
  int dimNSP = eSelect.NSP.rows();
  int nbVert = EXT.rows();
  int nbEqua = dimNSP + nbVert - 1;
  MyMatrix<T> ListEquation(n, nbEqua);
  MyVector<T> ListB(nbEqua);
  int idx = 0;
  MyVector<T> eV(n);
  for (int i = 0; i < n; i++)
    eV(i) = EXT(0, i + 1);
  T normV = EvaluationQuadForm<T, T>(GramMat, eV);
  for (int iVert = 1; iVert < nbVert; iVert++) {
    MyVector<T> fV(n);
    for (int i = 0; i < n; i++)
      fV(i) = EXT(iVert, i + 1);
    ListB(idx) = normV - EvaluationQuadForm<T, T>(GramMat, fV);
    for (int i = 0; i < n; i++) {
      T vSum(0);
      for (int j = 0; j < n; j++)
        vSum += GramMat(i, j) * (eV(j) - fV(j));
      ListEquation(i, idx) = 2 * vSum;
    }
    idx++;
  }
  for (int iNSP = 0; iNSP < dimNSP; iNSP++) {
    ListB(idx) = eSelect.NSP(iNSP, 0);
    for (int i = 0; i < n; i++)
      ListEquation(i, idx) = -eSelect.NSP(iNSP, i + 1);
    idx++;
  }
  std::optional<MyVector<T>> opt = SolutionMat(ListEquation, ListB);
#ifdef SANITY_CHECK_FUNDAMENTAL_DELAUNAY
  if (!opt) {
    std::cerr << "DEL: Error in CenterDelaunayPolytopeGeneral\n";
    std::cerr << "DEL: The linear system has no solution\n";
    throw TerminalException{1};
  }
#endif
  MyVector<T> const &eSol = *opt;
  MyVector<T> eCent(1 + n);
  eCent(0) = 1;
  for (int i = 0; i < n; i++)
    eCent(i + 1) = eSol(i);
  return eCent;
}

// Center and squared radius, general version: the center (the expensive part)
// plus one quadratic form for the radius.
template <typename T>
CP<T> CenterRadiusDelaunayPolytopeGeneral(MyMatrix<T> const &GramMat,
                                          MyMatrix<T> const &EXT) {
  MyVector<T> eCent = CenterDelaunayPolytopeGeneral(GramMat, EXT);
  T SquareRadius = SquareRadiusFromCenter(GramMat, EXT, eCent);
#ifdef DEBUG_DELAUNAY
  int n = GramMat.rows();
  int nbVert = EXT.rows();
  for (int iVert = 1; iVert < nbVert; iVert++) {
    MyVector<T> eW(n);
    for (int i = 0; i < n; i++)
      eW(i) = eCent(i + 1) - EXT(iVert, i + 1);
    T eNorm = EvaluationQuadForm<T, T>(GramMat, eW);
    if (SquareRadius != eNorm) {
      std::cerr << "DEL: SquareRadius=" << SquareRadius << " eNorm=" << eNorm
                << "\n";
      std::cerr << "DEL: Error in the Solutioning\n";
      throw TerminalException{1};
    }
  }
#endif
  return {std::move(SquareRadius), std::move(eCent)};
}

// Center of the circumscribing sphere of a FULL-DIMENSIONAL Delaunay polytope
// (its vertices affinely span the whole space). The general routine must first
// determine the affine span of EXT (TMat_SelectRowCol) so that the center is
// defined for a lower-dimensional polytope; when EXT is known to be
// full-dimensional the equidistance equations alone determine the center, so we
// skip the span computation entirely -- one fewer elimination and fewer
// equations. Same center as the general routine (checked under SANITY).
template <typename T>
MyVector<T> CenterDelaunayPolytope(MyMatrix<T> const &GramMat,
                                   MyMatrix<T> const &EXT) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int n = GramMat.rows();
  int nbVert = EXT.rows();
#ifdef SANITY_CHECK_FUNDAMENTAL_DELAUNAY
  int rnk = RankMat(EXT);
  if (rnk != n + 1) {
    std::cerr << "DEL: CenterDelaunayPolytope called on a NON-full-dimensional "
                 "polytope: rank(EXT)="
              << rnk << " != n+1=" << (n + 1) << "\n";
    throw TerminalException{1};
  }
#endif
  int nbEqua = nbVert - 1;
  MyMatrix<T> ListEquation(n, nbEqua);
  MyVector<T> ListB(nbEqua);
  MyVector<T> eV(n);
  for (int i = 0; i < n; i++)
    eV(i) = EXT(0, i + 1);
  T normV = EvaluationQuadForm<T, T>(GramMat, eV);
  for (int iVert = 1; iVert < nbVert; iVert++) {
    MyVector<T> fV(n);
    for (int i = 0; i < n; i++)
      fV(i) = EXT(iVert, i + 1);
    ListB(iVert - 1) = normV - EvaluationQuadForm<T, T>(GramMat, fV);
    for (int i = 0; i < n; i++) {
      T vSum(0);
      for (int j = 0; j < n; j++)
        vSum += GramMat(i, j) * (eV(j) - fV(j));
      ListEquation(i, iVert - 1) = 2 * vSum;
    }
  }
  std::optional<MyVector<T>> opt = SolutionMat(ListEquation, ListB);
#ifdef SANITY_CHECK_FUNDAMENTAL_DELAUNAY
  if (!opt) {
    std::cerr << "DEL: CenterDelaunayPolytope: the equidistance system has no "
                 "solution (polytope not co-spherical?)\n";
    throw TerminalException{1};
  }
#endif
  MyVector<T> const &eSol = *opt;
  MyVector<T> eCent(1 + n);
  eCent(0) = 1;
  for (int i = 0; i < n; i++)
    eCent(i + 1) = eSol(i);
#ifdef SANITY_CHECK_FUNDAMENTAL_DELAUNAY
  MyVector<T> ref = CenterDelaunayPolytopeGeneral<T>(GramMat, EXT);
  bool ok = (ref.size() == eCent.size());
  for (int i = 0; ok && i < eCent.size(); i++)
    if (ref(i) != eCent(i))
      ok = false;
  if (!ok) {
    std::cerr << "DEL: CenterDelaunayPolytope disagrees with the general "
                 "routine\n";
    throw TerminalException{1};
  }
#endif
  return eCent;
}

// Center and squared radius, full-dimensional version.
template <typename T>
CP<T> CenterRadiusDelaunayPolytope(MyMatrix<T> const &GramMat,
                                   MyMatrix<T> const &EXT) {
  MyVector<T> eCent = CenterDelaunayPolytope(GramMat, EXT);
  T SquareRadius = SquareRadiusFromCenter(GramMat, EXT, eCent);
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
template <typename T> struct AdjacentDelaunayPointSolver {
private:
  int n;
  MyMatrix<T> const &GramMat;
  MyMatrix<T> const &EXT;
  Face const &eInc;
  std::vector<T> LineGramMat;
  MyVector<T> ePt;
  MyVector<T> eDir;
  MyVector<T> v1;
  T G_v1;
  MyVector<T> two_v1_G;
  std::ostream &os;

public:
  AdjacentDelaunayPointSolver(MyMatrix<T> const &_GramMat,
                              MyMatrix<T> const &_EXT, Face const &_eInc,
                              std::ostream &_os)
      : n(_GramMat.rows()), GramMat(_GramMat), EXT(_EXT), eInc(_eInc),
        LineGramMat(GetLineVector(GramMat)), ePt(GramMat.rows()),
        eDir(GramMat.rows()), v1(n), os(_os) {
    std::vector<size_t> V;
    size_t n_vert = EXT.rows();
    for (size_t i = 0; i < n_vert; i++) {
      if (eInc[i] == 1) {
        V.push_back(i);
      }
    }
    size_t cnt = V.size();
    for (int i = 0; i < n; i++) {
      v1(i) = EXT(V[0], i + 1);
    }
    G_v1 = EvaluateLineVector(LineGramMat, v1);
    two_v1_G = T(2) * GramMat * v1;
    MyMatrix<T> Equa(cnt - 1, n + 1);
    MyVector<T> v2(n);
    for (size_t i_vert = 0; i_vert < cnt - 1; i_vert++) {
      for (int i = 0; i < n; i++) {
        v2(i) = EXT(V[i_vert + 1], i + 1);
      }
      T G_v2 = EvaluateLineVector(LineGramMat, v2);
      MyVector<T> two_v2_G = T(2) * GramMat * v2;
      Equa(i_vert, 0) = G_v1 - G_v2;
      for (int i = 0; i < n; i++) {
        Equa(i_vert, i + 1) = two_v2_G(i) - two_v1_G(i);
      }
    }
    MyMatrix<T> NSP = NullspaceTrMat(Equa);
#ifdef DEBUG_FUNDAMENTAL_DELAUNAY
    if (NSP.rows() != 2) {
      std::cerr << "DEL: Error in the computation of the kernel\n";
      throw TerminalException{1};
    }
#endif
    auto get_point_direction = [&]() -> void {
      for (int pos = 0; pos < 2; pos++) {
        if (NSP(pos, 0) != 0) {
          for (int i = 0; i < n; i++) {
            ePt(i) = NSP(pos, i + 1) / NSP(pos, 0);
          }
          for (int i = 0; i < n; i++) {
            eDir(i) = NSP(1 - pos, i + 1) - ePt(i) * NSP(1 - pos, 0);
          }
          return;
        }
      }
      std::cerr << "DEL: Failed to find a non-zero entry\n";
      throw TerminalException{1};
    };
    get_point_direction();
  }
  // v2 should be a point outside of the plane
  CP<T> GetCenterRadius(MyVector<T> const &v2) {
    T G_v2 = EvaluateLineVector(LineGramMat, v2);
    MyVector<T> two_v1_m_v2_G = two_v1_G - T(2) * GramMat * v2;
    T C_cst = G_v1 - G_v2 - two_v1_m_v2_G.dot(ePt);
    T D_cst = two_v1_m_v2_G.dot(eDir);
    T t = C_cst / D_cst;
    MyVector<T> eCent = ePt + t * eDir;
    MyVector<T> delta = eCent - v1;
    T eSqrDist = EvaluateLineVector(LineGramMat, delta);
    MyVector<T> eCentRet(n + 1);
    eCentRet(0) = T(1);
    for (int i = 0; i < n; i++)
      eCentRet(i + 1) = eCent(i);
    return {eSqrDist, eCentRet};
  }
  T GetRadius(MyVector<T> const &v2) {
    T G_v2 = EvaluateLineVector(LineGramMat, v2);
    MyVector<T> two_v1_m_v2_G = two_v1_G - T(2) * GramMat * v2;
    T C_cst = G_v1 - G_v2 - two_v1_m_v2_G.dot(ePt);
    T D_cst = two_v1_m_v2_G.dot(eDir);
    T t = C_cst / D_cst;
    MyVector<T> eCent = ePt + t * eDir;
    MyVector<T> delta = eCent - v1;
    return EvaluateLineVector(LineGramMat, delta);
  }
};

template <typename T, typename Tint, typename Tsolver>
MyMatrix<Tint> FindAdjacentDelaunayPolytope(
    Tsolver const &solver,
    MyMatrix<Tint> const &ShvGraverBasis, MyMatrix<T> const &EXT,
    SubsetRankOneSolver<T> &ext_solver, Face const &eInc, std::ostream &os) {
#ifdef TIMINGS_FUNDAMENTAL_DELAUNAY
  MicrosecondTime time;
#endif
  MyMatrix<T> const& GramMat = solver.GramMat;
  int dim = GramMat.rows();
  MyVector<T> TheFac = ext_solver.GetPositiveKernelVector(eInc);
  auto get_iColFind = [&]() -> int {
    for (int iCol = 0; iCol < dim; iCol++)
      if (TheFac(1 + iCol) != 0)
        return iCol;
    std::cerr << "Failed to find the matching iCol\n";
    throw TerminalException{1};
  };
  int iColFind = get_iColFind();
  T delta(-T_sign(TheFac(1 + iColFind)));
  int jRow = eInc.find_first();
  MyVector<T> SelectedVertex(dim);
  for (int i = 0; i < dim; i++)
    SelectedVertex(i) = EXT(jRow, i + 1);
  SelectedVertex(iColFind) += delta;
  AdjacentDelaunayPointSolver<T> adps(GramMat, EXT, eInc, os);
  T MinRadius = adps.GetRadius(SelectedVertex);
  auto fGraverUpdate = [&]() -> void {
#ifdef TIMINGS_FUNDAMENTAL_DELAUNAY
    MicrosecondTime time_graver;
#endif
#ifdef DEBUG_FUNDAMENTAL_DELAUNAY
    size_t n_update = 0, n_loop = 0;
#endif
    while (true) {
      bool IsImprovement = false;
      MyVector<T> NewTestVert(dim);
      int n_graver = ShvGraverBasis.rows();
      for (int i_graver = 0; i_graver < n_graver; i_graver++) {
        for (int i = 0; i < dim; i++) {
          NewTestVert(i) = SelectedVertex(i) + ShvGraverBasis(i_graver, i);
        }
        T eScal = TheFac(0);
        for (int i = 0; i < dim; i++)
          AddMul(eScal, NewTestVert(i), TheFac(i + 1));
        if (eScal < 0) {
          T TheRadius = adps.GetRadius(NewTestVert);
          if (TheRadius < MinRadius) {
            IsImprovement = true;
            SelectedVertex = NewTestVert;
            MinRadius = TheRadius;
#ifdef DEBUG_FUNDAMENTAL_DELAUNAY
            n_update++;
#endif
          }
        }
      }
#ifdef DEBUG_FUNDAMENTAL_DELAUNAY
      n_loop++;
#endif
      if (!IsImprovement) {
#ifdef DEBUG_FUNDAMENTAL_DELAUNAY
        os << "DEL_ENUM: n_loop=" << n_loop << " n_update=" << n_update
           << " |ShvGraverBasis|=" << n_graver << "\n";
#endif
#ifdef TIMINGS_FUNDAMENTAL_DELAUNAY
        os << "|DEL: fGraverUpdate|=" << time_graver << "\n";
#endif
        break;
      }
    }
  };
  fGraverUpdate();
  auto get_reply = [&]() -> resultCVP<T, Tint> {
    while (true) {
      CP<T> eCP = adps.GetCenterRadius(SelectedVertex);
      MyVector<T> eCenter(dim);
      for (int i = 0; i < dim; i++)
        eCenter(i) = eCP.eCent(i + 1);
      resultCVP<T, Tint> reply = solver.nearest_vectors(eCenter);
      if (reply.TheNorm == eCP.SquareRadius) {
#ifdef SANITY_CHECK_FUNDAMENTAL_DELAUNAY
        size_t n_error = 0;
        std::unordered_set<MyVector<T>> setEXT_adj;
        int n_row = reply.ListVect.rows();
        MyMatrix<T> EXT_adj_aff(n_row, dim + 1);
        for (int iRow=0; iRow<n_row; iRow++) {
          MyVector<Tint> v_i = GetMatrixRow(reply.ListVect, iRow);
          MyVector<T> v = UniversalVectorConversion<T,Tint>(v_i);
          EXT_adj_aff(iRow, 0) = 1;
          for (int i=0; i<dim; i++) {
            EXT_adj_aff(iRow, i+1) = v(i);
          }
          setEXT_adj.insert(v);
        }
        if (RankMat(EXT_adj_aff) != dim + 1) {
          n_error += 1;
          std::cerr << "DEL: Incorrect rank of output EXT_adj_aff\n";
        }
        for (int i=0; i<EXT.rows(); i++) {
          MyVector<T> V(dim);
          for (int iCol=0; iCol<dim; iCol++) {
            V(iCol) = EXT(i, iCol+1);
          }
          if (eInc[i] == 1) {
            if (!setEXT_adj.contains(V)) {
              n_error += 1;
              std::cerr << "DEL: Failed to contain a vertex from eInc\n";
            }
          } else {
            if (setEXT_adj.contains(V)) {
              n_error += 1;
              std::cerr << "DEL: Containing a vertex out of eInc\n";
            }
          }
        }
        if (n_error > 0) {
          std::cerr << "DEL: |EXT|=" << EXT.rows() << " / " << EXT.cols()
                    << " |eInc|=" << eInc.count() << " / " << eInc.size() << "\n";
          std::cerr << "DEL: There has been n_error=" << n_error << "\n";
          throw TerminalException{1};
        }
#endif
        return reply;
      }
      for (int i = 0; i < dim; i++)
        SelectedVertex(i) =
            UniversalScalarConversion<T, Tint>(reply.ListVect(0, i));
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
#ifdef SANITY_CHECK_FUNDAMENTAL_DELAUNAY
  if (RankMat(RetEXT) != dim+1) {
    std::cerr << "DEL: RetEXT should be of the right rank\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS_FUNDAMENTAL_DELAUNAY
  os << "|DEL: FindAdjacentDelaunayPolytope|=" << time << "\n";
#endif
  return RetEXT;
}

template <typename T, typename Tint, typename Tsolver>
MyMatrix<Tint> FindDelaunayPolytope_random(Tsolver const &solver,
                                           MyMatrix<Tint> const &ShvGraverBasis,
                                           int const &target_ext, int max_iter,
                                           std::string method,
                                           std::ostream &os) {
  MyMatrix<Tint> EXTret = FindDelaunayPolytope<T, Tint>(solver, os);
  for (int iter = 0; iter < max_iter; iter++) {
    MyMatrix<T> EXTret_T = UniversalMatrixConversion<T, Tint>(EXTret);
#ifdef DEBUG_FUNDAMENTAL_DELAUNAY
    os << "iter=" << iter << "/" << max_iter << " method=" << method
       << " |EXTret|=" << EXTret.rows() << "\n";
#endif
    vectface vf = DirectComputationInitialFacetSet<T>(EXTret_T, method, os);
#ifdef DEBUG_FUNDAMENTAL_DELAUNAY
    os << "  |vf|=" << vf.size() << "\n";
#endif
    SubsetRankOneSolver<T> ext_solver(EXTret_T);
    auto get_adj_cand = [&]() -> MyMatrix<Tint> {
      MyMatrix<Tint> EXTwork;
      int n_max = std::numeric_limits<int>::max();
      for (auto &eFace : vf) {
        MyMatrix<Tint> EXTadj = FindAdjacentDelaunayPolytope<T, Tint>(
            solver, ShvGraverBasis, EXTret_T, ext_solver, eFace, os);
        int nbRow = EXTadj.rows();
        if (nbRow < n_max) {
          EXTwork = EXTadj;
          n_max = nbRow;
        }
      }
      return EXTwork;
    };
    MyMatrix<Tint> EXTcand = get_adj_cand();
    if (target_ext > 0) {
      if (EXTcand.rows() <= target_ext) {
        return EXTcand;
      }
    } else {
      // If target_ext = 0, we terminate if the sampling does not give us
      // anything better.
      if (EXTcand.rows() >= EXTret.rows()) {
        return EXTret;
      }
    }
    EXTret = EXTcand;
  }
  return EXTret;
}

template <typename T>
std::vector<MyMatrix<T>> CharacterizingPair(MyMatrix<T> const &GramMat,
                                            MyVector<T> const &TheCenter,
                                            std::ostream &os) {
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
  T eScal(1);
  for (int i = 0; i < n; i++)
    AddMul(eScal, eVect(i), TheCenterRed(i));
  Mat1(0, 0) = eScal;
  MyMatrix<T> Mat2 = ZeroMatrix<T>(n + 1, n + 1);
  Mat2(0, 0) = 1;
  while (true) {
    bool test = IsPositiveDefinite(Mat1, os);
    if (test)
      break;
    Mat1(0, 0)++;
  }
  return {Mat1, Mat2};
}

template <typename T> struct ResultCov {
  double CovDensityNormalized;
  double CovDensity;
  double CoveringRadius;
  double DeterminantLattice;
  double VolumeBall;
  double TheCov_d;
  int TheDim;
  T TheDet;
  T TheCov;
};

template <typename T>
ResultCov<T> ComputeCoveringDensityFromDimDetCov(int TheDim, T TheDet,
                                                 T TheCov) {

  double TheCov_d = UniversalScalarConversion<double, T>(TheCov);
  double TheCovRad = std::sqrt(TheCov_d);
  double TheDet_d = UniversalScalarConversion<double, T>(TheDet);
  double TheDetLattice = std::sqrt(TheDet_d);
  //
  double val = 1.0;
  for (int i = 0; i < TheDim; i++) {
    val *= TheCov_d;
  }
  double quot = val / TheDet_d;
  double red = std::sqrt(quot);
  //
  auto get_volume_ball = [&]() -> double {
    double pi = 2 * std::acos(0);
    int res = TheDim % 2;
    if (res == 0) {
      int half_d = TheDim / 2;
      double powpi = 1.0;
      double fact = 1.0;
      for (int i = 1; i <= half_d; i++) {
        powpi *= pi;
        fact *= i;
      }
      return powpi / fact;
    } else {
      int half_d = (TheDim - 1) / 2;
      double powpi = 1.0;
      for (int i = 1; i <= half_d; i++) {
        powpi *= pi;
      }
      double fact = 1.0;
      for (int i = 0; i <= half_d; i++) {
        fact *= (0.5 + i);
      }
      return powpi / fact;
    }
  };
  double VolumeBall = get_volume_ball();
  double CovDens = red * VolumeBall;
  return {red,      CovDens, TheCovRad, TheDetLattice, VolumeBall,
          TheCov_d, TheDim,  TheDet,    TheCov};
}

template <typename T> std::string to_stringGAP(ResultCov<T> const &x) {
  return "rec(CovDensityNormalized:=" + std::to_string(x.CovDensityNormalized) +
         ", CovDensity:=" + std::to_string(x.CovDensity) +
         ", CoveringRadius:=" + std::to_string(x.CoveringRadius) +
         ", DeterminantLattice:=" + std::to_string(x.DeterminantLattice) +
         ", VolumeBall:=" + std::to_string(x.VolumeBall) +
         ", TheCov_d:=" + std::to_string(x.TheCov_d) +
         ", TheDim:=" + std::to_string(x.TheDim) +
         ", TheDet:=" + std::to_string(x.TheDet) +
         ", TheCov:=" + std::to_string(x.TheCov) + ")";
}

// clang-format off
#endif  // SRC_LATT_FUNDAMENTALDELAUNAY_H_
// clang-format on

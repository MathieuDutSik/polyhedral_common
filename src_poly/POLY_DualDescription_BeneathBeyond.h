// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DUALDESCRIPTION_BENEATHBEYOND_H_
#define SRC_POLY_POLY_DUALDESCRIPTION_BENEATHBEYOND_H_

// clang-format off
#include "POLY_Fundamental.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include <map>
#include <set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_BENEATH_BEYOND
#endif

#ifdef DISABLE_DEBUG_BENEATH_BEYOND
#undef DEBUG_BENEATH_BEYOND
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_BENEATH_BEYOND
#endif

#ifdef TIMINGS
#define TIMINGS_BENEATH_BEYOND
#endif

// Beneath-and-beyond convex hull method (a.k.a. the incremental / placing
// algorithm; see Joswig "Beneath-and-Beyond revisited" and the polymake
// implementation beneath_beyond_impl.h).
//
// This is a native, in-idiom port for the case that is used throughout the
// rest of src_poly: EXT is the set of extreme rays of a *full-dimensional
// pointed cone*, one ray per row in homogeneous coordinates. A facet is a
// linear inequality f with f . x >= 0 on every row and f . x == 0 exactly on
// the incident rows.
//
// Precondition: RankMat(EXT) == EXT.cols() (full dimensional). The lineality
// space / low-dimensional / redundancy-detecting generality of the polymake
// original is deliberately not handled here; feed such inputs through the
// existing preprocessing (or use cdd/lrs) instead.
//
// The facet enumeration requires T to be a field (facet normals are obtained
// through FindFacetInequality / a nullspace computation). The triangulation
// helper below only uses determinants and a row basis, so it also accepts ring
// types.
//
// Facet enumeration: the rays are inserted one at a time. After the first
// d = EXT.cols() linearly independent rays (a simplicial cone) are set up with
// their d facets, each further ray p is classified against every current facet
// by sign(normal . p):
//   < 0  beneath   : the facet is visible from p and dies;
//   > 0  beyond    : the facet survives untouched;
//   == 0 incident  : p lies on the facet.
// Every ridge shared by a dying (visible) facet and a surviving *strictly
// positive* facet is a horizon ridge; ridge u {p} spans a new facet. Interior
// rays (no visible facet) produce nothing.
//
// A key simplification versus the polymake code: each facet stores its full
// incidence over *all* rows, computed once from its normal at creation. Because
// a facet's supporting hyperplane is fixed for its whole lifetime, its
// incidence {r : normal . r == 0} is time invariant, so incident rays are
// absorbed automatically and no incremental incidence bookkeeping is needed.

template <typename T> struct BeneathBeyondFacet {
  MyVector<T> normal; // normal . row >= 0 on every row, == 0 on the incidence
  Face incd;          // incidence over all EXT rows
};

namespace beneath_beyond {

// Scalar product normal . EXT.row(iRow).
template <typename T>
inline T facet_scal(MyMatrix<T> const &EXT, MyVector<T> const &normal,
                    int iRow) {
  int nbCol = EXT.cols();
  T sum(0);
  for (int iCol = 0; iCol < nbCol; iCol++)
    sum += normal(iCol) * EXT(iRow, iCol);
  return sum;
}

// Full incidence {r : normal . r == 0} of a facet normal over all rows.
template <typename T>
Face facet_incidence(MyMatrix<T> const &EXT, MyVector<T> const &normal) {
  int nbRow = EXT.rows();
  MyVector<T> prod = EXT * normal;
  Face f(nbRow);
  for (int iRow = 0; iRow < nbRow; iRow++)
    if (prod(iRow) == 0)
      f[iRow] = 1;
  return f;
}

// The sorted list of set positions of a face.
inline std::vector<int> face_to_vector(Face const &f) {
  std::vector<int> ListRow;
  boost::dynamic_bitset<>::size_type iRow = f.find_first();
  while (iRow != boost::dynamic_bitset<>::npos) {
    ListRow.push_back(static_cast<int>(iRow));
    iRow = f.find_next(iRow);
  }
  return ListRow;
}

// Rank of the sub-family of rows selected by the face.
template <typename T> int face_rank(MyMatrix<T> const &EXT, Face const &f) {
  return RankMat(SelectRow(EXT, face_to_vector(f)));
}

} // namespace beneath_beyond

// Core kernel: returns the facets (normal + full incidence) of the cone spanned
// by the rows of EXT.
template <typename T>
std::vector<BeneathBeyondFacet<T>>
BeneathBeyond_Kernel(MyMatrix<T> const &EXT,
                     [[maybe_unused]] std::ostream &os) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
#ifdef TIMINGS_BENEATH_BEYOND
  MicrosecondTime time;
#endif
#ifdef SANITY_CHECK_BENEATH_BEYOND
  if (nbCol < 2) {
    std::cerr << "BeneathBeyond: requires EXT.cols() >= 2, nbCol=" << nbCol
              << "\n";
    throw TerminalException{1};
  }
  int rank = RankMat(EXT);
  if (rank != nbCol) {
    std::cerr << "BeneathBeyond: EXT must span a full-dimensional cone, "
              << "RankMat(EXT)=" << rank << " but nbCol=" << nbCol << "\n";
    throw TerminalException{1};
  }
#endif

  // The initial simplicial cone: d = nbCol linearly independent rays.
  std::vector<int> basis = TMat_ListRowSelect(EXT);
#ifdef SANITY_CHECK_BENEATH_BEYOND
  if (static_cast<int>(basis.size()) != nbCol) {
    std::cerr << "BeneathBeyond: could not extract a full-rank basis, "
              << "|basis|=" << basis.size() << " nbCol=" << nbCol << "\n";
    throw TerminalException{1};
  }
#endif
  Face in_basis(nbRow);
  for (auto &eRow : basis)
    in_basis[eRow] = 1;

  auto make_facet = [&](Face const &seed) -> BeneathBeyondFacet<T> {
    MyVector<T> normal = FindFacetInequality(EXT, seed);
    normal = ScalarCanonicalizationVector(normal);
    Face incd = beneath_beyond::facet_incidence(EXT, normal);
    return {std::move(normal), std::move(incd)};
  };

  std::vector<BeneathBeyondFacet<T>> facets;
  // The nbCol facets of the simplicial cone: drop one basis ray at a time.
  for (int iBas = 0; iBas < nbCol; iBas++) {
    Face seed(nbRow);
    for (int jBas = 0; jBas < nbCol; jBas++)
      if (jBas != iBas)
        seed[basis[jBas]] = 1;
    facets.push_back(make_facet(seed));
  }

  // Insert the remaining rays one at a time.
  for (int p = 0; p < nbRow; p++) {
    if (in_basis[p])
      continue;
#ifdef DEBUG_BENEATH_BEYOND
    os << "BENEATH_BEYOND: inserting ray p=" << p
       << " |facets|=" << facets.size() << "\n";
#endif
    // Classify the current facets w.r.t. p.
    std::vector<size_t> visible, positive;
    for (size_t i = 0; i < facets.size(); i++) {
      T scal = beneath_beyond::facet_scal(EXT, facets[i].normal, p);
      if (scal < 0)
        visible.push_back(i);
      else if (scal > 0)
        positive.push_back(i);
      // scal == 0: incident, facet already contains p via its full incidence.
    }
    if (visible.empty())
      continue; // interior / boundary ray, nothing changes

    // Mark the visible facets; the others (strictly positive and incident)
    // survive this step.
    Face to_delete(facets.size());
    for (size_t iv : visible)
      to_delete[iv] = 1;

    // Horizon ridges: a ridge shared by a visible and a strictly positive facet
    // (rank nbCol-2). ridge u {p} spans a new facet. New facets are deduplicated
    // by normal against each other AND against the surviving facets: in
    // degenerate position an incident facet created earlier may already coincide
    // with a horizon facet of this step (its full-incidence scan already listed
    // p), and must not be duplicated.
    std::vector<BeneathBeyondFacet<T>> new_facets;
    std::set<MyVector<T>> seen_normals;
    for (size_t i = 0; i < facets.size(); i++)
      if (!to_delete[i])
        seen_normals.insert(facets[i].normal);
    for (size_t iv : visible) {
      for (size_t ip : positive) {
        Face ridge = facets[iv].incd & facets[ip].incd;
        if (static_cast<int>(ridge.count()) < nbCol - 2)
          continue;
        if (beneath_beyond::face_rank(EXT, ridge) != nbCol - 2)
          continue;
        Face seed = ridge;
        seed[p] = 1;
        BeneathBeyondFacet<T> nf = make_facet(seed);
        if (seen_normals.insert(nf.normal).second)
          new_facets.push_back(std::move(nf));
      }
    }

    // Remove the visible facets (compact) and append the new ones.
    std::vector<BeneathBeyondFacet<T>> kept;
    kept.reserve(facets.size() - visible.size() + new_facets.size());
    for (size_t i = 0; i < facets.size(); i++)
      if (!to_delete[i])
        kept.push_back(std::move(facets[i]));
    for (auto &nf : new_facets)
      kept.push_back(std::move(nf));
    facets = std::move(kept);
  }

#ifdef TIMINGS_BENEATH_BEYOND
  os << "BENEATH_BEYOND: |EXT|=" << nbRow << "/" << nbCol
     << " |facets|=" << facets.size() << " time=" << time << "\n";
#endif
  return facets;
}

// Facet incidences, matching DirectFacetComputationIncidence.
template <typename T>
vectface POLY_DualDescription_BeneathBeyondIncidence(MyMatrix<T> const &EXT,
                                                     std::ostream &os) {
  std::vector<BeneathBeyondFacet<T>> facets = BeneathBeyond_Kernel(EXT, os);
  vectface vf(EXT.rows());
  for (auto &facet : facets)
    vf.push_back(facet.incd);
  return vf;
}

// Facet normals, matching DirectFacetComputationInequalities.
template <typename T>
MyMatrix<T> POLY_DualDescription_BeneathBeyondInequalities(MyMatrix<T> const &EXT,
                                                           std::ostream &os) {
  std::vector<BeneathBeyondFacet<T>> facets = BeneathBeyond_Kernel(EXT, os);
  int n_facet = facets.size();
  int nbCol = EXT.cols();
  MyMatrix<T> FAC(n_facet, nbCol);
  for (int i = 0; i < n_facet; i++)
    AssignMatrixRow(FAC, i, facets[i].normal);
  return FAC;
}

// (Face, inequality) pairs, matching DirectFacetComputationFaceIneq.
template <typename T, typename Fprocess>
void POLY_DualDescription_BeneathBeyondFaceIneq(MyMatrix<T> const &EXT,
                                                Fprocess f_process,
                                                std::ostream &os) {
  std::vector<BeneathBeyondFacet<T>> facets = BeneathBeyond_Kernel(EXT, os);
  for (auto &facet : facets) {
    std::pair<Face, MyVector<T>> pair{facet.incd, facet.normal};
    f_process(pair);
  }
}

// Placing triangulation produced by the beneath-and-beyond insertion order:
// a list of simplices, each the sorted list of its nbCol vertex row indices
// (same format as lrs::GetTriangulation).
//
// This is computed by the direct/robust definition of the placing
// triangulation (Lee): start from a simplex on d = nbCol independent rays; for
// each further ray p, cone p over every boundary (d-1)-cell of the current
// triangulation that is strictly visible from p (p on the opposite side of the
// cell's hyperplane from the cell's apex). Redundant rays see no boundary cell
// and are skipped automatically. The boundary cells are recomputed each step
// (a (d-1)-face is on the boundary iff it belongs to exactly one simplex),
// which is simpler and less error-prone than threading the triangulation
// through the facet update, at the cost of some speed.
template <typename T>
std::vector<std::vector<int>>
POLY_DualDescription_BeneathBeyondTriangulation(MyMatrix<T> const &EXT,
                                                [[maybe_unused]] std::ostream &os) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
#ifdef TIMINGS_BENEATH_BEYOND
  MicrosecondTime time;
#endif
#ifdef SANITY_CHECK_BENEATH_BEYOND
  int rank = RankMat(EXT);
  if (rank != nbCol) {
    std::cerr << "BeneathBeyond(triang): EXT must span a full-dimensional cone, "
              << "RankMat(EXT)=" << rank << " but nbCol=" << nbCol << "\n";
    throw TerminalException{1};
  }
#endif

  // Signed d x d determinant of the rows (face vertices, then extra), where the
  // face keeps its stored order so signs are comparable across calls.
  auto det_face_point = [&](std::vector<int> const &face, int extra) -> T {
    MyMatrix<T> M(nbCol, nbCol);
    for (int i = 0; i < nbCol - 1; i++)
      for (int j = 0; j < nbCol; j++)
        M(i, j) = EXT(face[i], j);
    for (int j = 0; j < nbCol; j++)
      M(nbCol - 1, j) = EXT(extra, j);
    return DeterminantMat(M);
  };

  std::vector<int> basis = TMat_ListRowSelect(EXT);
  Face in_basis(nbRow);
  for (auto &eRow : basis)
    in_basis[eRow] = 1;

  std::vector<std::vector<int>> simplices;
  {
    std::vector<int> simplex = basis;
    std::sort(simplex.begin(), simplex.end());
    simplices.push_back(std::move(simplex));
  }

  for (int p = 0; p < nbRow; p++) {
    if (in_basis[p])
      continue;
    // Boundary (d-1)-cells of the current triangulation: a face is on the
    // boundary iff exactly one simplex contains it. Keep one apex (the dropped
    // vertex) to orient the visibility test.
    std::map<std::vector<int>, std::pair<int, int>> faces; // face -> (count, apex)
    for (auto &simplex : simplices) {
      for (int i = 0; i < nbCol; i++) {
        std::vector<int> face;
        face.reserve(nbCol - 1);
        for (int j = 0; j < nbCol; j++)
          if (j != i)
            face.push_back(simplex[j]);
        auto &ref = faces[face];
        ref.first++;
        ref.second = simplex[i];
      }
    }
    std::vector<std::vector<int>> new_simplices;
    for (auto &kv : faces) {
      if (kv.second.first != 1)
        continue; // interior face
      std::vector<int> const &face = kv.first;
      int apex = kv.second.second;
      T dp = det_face_point(face, p);
      if (dp == 0)
        continue; // p coplanar with the face: not strictly visible
      T da = det_face_point(face, apex);
      // p is visible from outside iff on the opposite side from the apex.
      if ((dp > 0) == (da > 0))
        continue;
      std::vector<int> simplex = face;
      simplex.push_back(p);
      std::sort(simplex.begin(), simplex.end());
      new_simplices.push_back(std::move(simplex));
    }
    for (auto &simplex : new_simplices)
      simplices.push_back(std::move(simplex));
  }

#ifdef TIMINGS_BENEATH_BEYOND
  os << "BENEATH_BEYOND: triangulation |EXT|=" << nbRow << "/" << nbCol
     << " |simplices|=" << simplices.size() << " time=" << time << "\n";
#endif
  return simplices;
}

// clang-format off
#endif  // SRC_POLY_POLY_DUALDESCRIPTION_BENEATHBEYOND_H_
// clang-format on

// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DUALDESCRIPTION_BENEATHBEYOND_H_
#define SRC_POLY_POLY_DUALDESCRIPTION_BENEATHBEYOND_H_

// clang-format off
#include "POLY_Fundamental.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include <algorithm>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
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
// The new facet normal is NOT recomputed from a nullspace. Its supporting
// hyperplane must contain the ridge (= the intersection of the two parent
// hyperplanes) and the new ray p, so it is the unique combination of the two
// parent normals killing p:
//   normal_new = (normal_ip . p) * normal_iv - (normal_iv . p) * normal_ip,
// with iv the visible parent (normal_iv . p < 0) and ip the positive parent
// (normal_ip . p > 0). This is a positive combination of two inward normals, so
// it is inward-oriented for free. Only the nbCol initial simplex facets need a
// genuine nullspace. (This mirrors the ratio-step of the FlippingFramework in
// POLY_Fundamental.h.)
//
// Each facet stores its full incidence over *all* rows, computed once from its
// normal at creation. Because a facet's supporting hyperplane is fixed for its
// whole lifetime, its incidence {r : normal . r == 0} is time invariant, so
// incident rays are absorbed automatically and no incremental incidence
// bookkeeping is needed. The incidence cannot be replaced by ridge u {p}: for
// non-simplicial facets (e.g. cut / metric polytopes) the true facet is
// strictly larger than the seed. Deduplication is by this incidence Face.

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

// Build a facet record from an (unnormalized) inward normal: canonicalize it
// (only to keep coordinates small over repeated flip combinations; sign is
// restored to the orientation of `raw`, so the inward direction is untouched)
// and scan its full incidence.
template <typename T>
BeneathBeyondFacet<T> facet_from_normal(MyMatrix<T> const &EXT,
                                        MyVector<T> const &raw) {
  int nbCol = EXT.cols();
  MyVector<T> normal = ScalarCanonicalizationVector(raw);
  for (int iCol = 0; iCol < nbCol; iCol++) {
    if (raw(iCol) != 0) {
      if ((normal(iCol) > 0) != (raw(iCol) > 0))
        normal = -normal;
      break;
    }
  }
  Face incd = facet_incidence(EXT, normal);
  return {std::move(normal), std::move(incd)};
}

// Same, but the orientation is fixed by an "interior" row instead of by the raw
// sign: the single EXT * normal pass is reused to both orient the facet (so it
// is positive on `interior_row`) and read off its incidence -- no separate
// scalar product is computed for the sign. Used for the initial simplex facets,
// oriented on the dropped basis ray.
template <typename T>
BeneathBeyondFacet<T> facet_from_normal_interior(MyMatrix<T> const &EXT,
                                                 MyVector<T> const &raw,
                                                 int interior_row) {
  int nbRow = EXT.rows();
  MyVector<T> normal = ScalarCanonicalizationVector(raw);
  MyVector<T> prod = EXT * normal;
  if (prod(interior_row) < 0) {
    normal = -normal;
    prod = -prod;
  }
  Face incd(nbRow);
  for (int iRow = 0; iRow < nbRow; iRow++)
    if (prod(iRow) == 0)
      incd[iRow] = 1;
  return {std::move(normal), std::move(incd)};
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

// Containment relation of two faces (polymake's incl): 0 if equal, -1 if
// A is a proper subset of B, 1 if A is a proper superset of B, 2 if the two
// are incomparable.
inline int face_incl(Face const &A, Face const &B) {
  Face AB = A & B;
  bool a_in_b = (AB == A);
  bool b_in_a = (AB == B);
  if (a_in_b && b_in_a)
    return 0;
  if (a_in_b)
    return -1;
  if (b_in_a)
    return 1;
  return 2;
}

// Test whether the rows selected by the face have rank exactly `target`, used
// with target = nbCol-2 to decide whether an intersection of two facets is a
// ridge. If the face carries fewer than `target` rays its rank cannot reach
// `target`, so we reject it before building anything; otherwise the rank is
// computed with the well-pivoted RankMat kernel (a bad pivot choice blows up
// the rational coefficients and dominates the cost, so this matters).
template <typename T>
bool face_rank_is(MyMatrix<T> const &EXT, Face const &f, int target) {
  if (static_cast<int>(f.count()) < target)
    return false;
  return face_rank(EXT, f) == target;
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

  auto facet_from_normal = [&](MyVector<T> const &raw) -> BeneathBeyondFacet<T> {
    return beneath_beyond::facet_from_normal(EXT, raw);
  };

  std::vector<BeneathBeyondFacet<T>> facets;
  // The nbCol facets of the simplicial cone: drop one basis ray at a time. These
  // are the only normals obtained from a nullspace; every later facet is a cheap
  // linear combination of two existing ones (see below). Orient each inward, so
  // it is positive on the dropped basis ray.
  for (int iBas = 0; iBas < nbCol; iBas++) {
    Face seed(nbRow);
    for (int jBas = 0; jBas < nbCol; jBas++)
      if (jBas != iBas)
        seed[basis[jBas]] = 1;
    MyVector<T> normal = FindFacetInequality(EXT, seed);
    facets.push_back(
        beneath_beyond::facet_from_normal_interior(EXT, normal, basis[iBas]));
  }

  // Insert the remaining rays one at a time.
  for (int p = 0; p < nbRow; p++) {
    if (in_basis[p])
      continue;
#ifdef DEBUG_BENEATH_BEYOND
    os << "BENEATH_BEYOND: inserting ray p=" << p
       << " |facets|=" << facets.size() << "\n";
#endif
    // Classify the current facets w.r.t. p, keeping the scalar products so the
    // new facet normals can be assembled from them without a nullspace.
    std::vector<T> scal_p(facets.size());
    std::vector<size_t> visible, positive;
    for (size_t i = 0; i < facets.size(); i++) {
      T scal = beneath_beyond::facet_scal(EXT, facets[i].normal, p);
      scal_p[i] = scal;
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

    // Horizon ridges: a ridge shared by a visible facet (scal < 0) and a
    // strictly positive facet (scal > 0), of rank nbCol-2. The new facet spans
    // ridge u {p}; its supporting hyperplane is the unique combination of the
    // two parent hyperplanes that also passes through p:
    //   normal = scal_p[ip] * normal[iv] - scal_p[iv] * normal[ip].
    // Since scal_p[ip] > 0 and -scal_p[iv] > 0 and both parents are >= 0 on the
    // whole cone, this combination is automatically inward-oriented, is zero on
    // the ridge, and is zero on p. No nullspace is needed.
    // Facets are deduplicated by their incidence Face, which merges coplanar
    // horizon facets (non-generic position) and rejects a horizon facet that
    // coincides with a surviving (incident) facet.
    std::vector<BeneathBeyondFacet<T>> new_facets;
    std::unordered_set<Face> seen_faces;
    for (size_t i = 0; i < facets.size(); i++)
      if (!to_delete[i])
        seen_faces.insert(facets[i].incd);
    for (size_t iv : visible) {
      for (size_t ip : positive) {
        Face ridge = facets[iv].incd & facets[ip].incd;
        if (!beneath_beyond::face_rank_is(EXT, ridge, nbCol - 2))
          continue;
        MyVector<T> normal =
            scal_p[ip] * facets[iv].normal - scal_p[iv] * facets[ip].normal;
        BeneathBeyondFacet<T> nf = facet_from_normal(normal);
        if (seen_faces.insert(nf.incd).second)
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

namespace beneath_beyond {
// A facet node of the dual graph: its supporting inequality, its incidence, the
// ids of the facets it shares a ridge with, and a liveness flag (deleted facets
// are tombstoned rather than erased, so ids stay stable).
template <typename T> struct DGFacet {
  MyVector<T> normal;
  Face incd;
  std::vector<int> adj;
  bool alive;
};
} // namespace beneath_beyond

// Dual-graph beneath-and-beyond. Same result as BeneathBeyond_Kernel, but the
// horizon is read off an incrementally-maintained facet adjacency graph instead
// of being rediscovered by testing every visible x positive pair -- which was
// 60-80% of the plain kernel's time on non-generic polytopes. Per inserted ray:
//   * classify the live facets by sign(normal . p)  (O(#facets));
//   * every graph edge from a visible facet to a strictly positive facet is a
//     horizon ridge -> a new facet (normal by the flip formula), with no rank
//     test, since adjacency already certifies the ridge;
//   * new/incident facets are rewired among themselves (a small local set) with
//     the rank-based ridge test;
//   * visible facets are tombstoned and unlinked.
template <typename T>
std::vector<BeneathBeyondFacet<T>>
BeneathBeyond_Kernel_DualGraph(MyMatrix<T> const &EXT,
                               [[maybe_unused]] std::ostream &os) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
#ifdef TIMINGS_BENEATH_BEYOND
  MicrosecondTime time;
#endif
#ifdef SANITY_CHECK_BENEATH_BEYOND
  if (nbCol < 2) {
    std::cerr << "BeneathBeyond(DG): requires EXT.cols() >= 2, nbCol=" << nbCol
              << "\n";
    throw TerminalException{1};
  }
  int rank = RankMat(EXT);
  if (rank != nbCol) {
    std::cerr << "BeneathBeyond(DG): EXT must span a full-dimensional cone, "
              << "RankMat(EXT)=" << rank << " but nbCol=" << nbCol << "\n";
    throw TerminalException{1};
  }
#endif

  std::vector<int> basis = TMat_ListRowSelect(EXT);
  Face in_basis(nbRow);
  for (auto &eRow : basis)
    in_basis[eRow] = 1;

  using DGF = beneath_beyond::DGFacet<T>;
  std::vector<DGF> F;
  auto add_facet = [&](BeneathBeyondFacet<T> &&bf) -> int {
    int id = F.size();
    F.push_back(DGF{std::move(bf.normal), std::move(bf.incd), {}, true});
    return id;
  };
  auto edge_exists = [&](int a, int b) -> bool {
    for (int x : F[a].adj)
      if (x == b)
        return true;
    return false;
  };
  auto add_edge = [&](int a, int b) {
    F[a].adj.push_back(b);
    F[b].adj.push_back(a);
  };

  // Initial simplicial cone: nbCol facets, every pair sharing a ridge, so the
  // dual graph starts as the complete graph.
  std::vector<int> init;
  for (int iBas = 0; iBas < nbCol; iBas++) {
    Face seed(nbRow);
    for (int jBas = 0; jBas < nbCol; jBas++)
      if (jBas != iBas)
        seed[basis[jBas]] = 1;
    MyVector<T> normal = FindFacetInequality(EXT, seed);
    init.push_back(add_facet(
        beneath_beyond::facet_from_normal_interior(EXT, normal, basis[iBas])));
  }
  for (size_t i = 0; i < init.size(); i++)
    for (size_t j = i + 1; j < init.size(); j++)
      add_edge(init[i], init[j]);

  for (int p = 0; p < nbRow; p++) {
    if (in_basis[p])
      continue;
    // Classify the live facets: the visible ones (to be replaced) and the
    // incident ones (facets p lies on, which survive and can gain adjacencies).
    std::vector<T> scal(F.size());
    std::vector<int> visible, incident;
    for (size_t i = 0; i < F.size(); i++) {
      if (!F[i].alive)
        continue;
      T s = beneath_beyond::facet_scal(EXT, F[i].normal, p);
      scal[i] = s;
      if (s < 0)
        visible.push_back(static_cast<int>(i));
      else if (s == 0)
        incident.push_back(static_cast<int>(i));
    }
    if (visible.empty())
      continue;
#ifdef DEBUG_BENEATH_BEYOND
    os << "BENEATH_BEYOND(DG): ray p=" << p << " |facets|=" << F.size()
       << " |visible|=" << visible.size() << " |incident|=" << incident.size()
       << "\n";
#endif

    // Create the new facets from the horizon ridges (visible -> strictly positive
    // edges), with the normal from the flip formula. Deduplicate by incidence
    // against each other and against the incident facets: in non-generic position
    // a horizon facet can coincide with an incident facet (whose full incidence
    // already lists p), which is then reused instead of duplicated.
    std::unordered_map<Face, int> incd_to_facet;
    for (int fi : incident)
      incd_to_facet.emplace(F[fi].incd, fi);
    std::vector<int> new_facets;
    std::vector<std::pair<int, int>> nf_parent; // (new facet, positive parent)
    for (int iv : visible) {
      std::vector<int> nbrs = F[iv].adj; // copy: add_facet may reallocate F
      for (int g : nbrs) {
        if (!F[g].alive || scal[g] <= 0)
          continue;
        MyVector<T> raw = scal[g] * F[iv].normal - scal[iv] * F[g].normal;
        BeneathBeyondFacet<T> bf = beneath_beyond::facet_from_normal(EXT, raw);
        auto it = incd_to_facet.find(bf.incd);
        int nf;
        if (it == incd_to_facet.end()) {
          Face key = bf.incd;
          nf = add_facet(std::move(bf));
          incd_to_facet.emplace(std::move(key), nf);
          new_facets.push_back(nf);
        } else {
          nf = it->second;
        }
        nf_parent.push_back({nf, g});
      }
    }
    for (auto &pr : nf_parent)
      if (!edge_exists(pr.first, pr.second))
        add_edge(pr.first, pr.second);

    // Tombstone the visible facets and unlink them BEFORE rewiring, so the
    // maximality test below only ever compares against surviving facets.
    for (int iv : visible) {
      for (int g : F[iv].adj) {
        auto &ga = F[g].adj;
        ga.erase(std::remove(ga.begin(), ga.end(), iv), ga.end());
      }
      F[iv].adj.clear();
      F[iv].alive = false;
    }

    // Rewire the facets that now touch p (new + incident) among themselves. A
    // pair (fa, fb) is joined iff R = fa ∩ fb is a ridge (rank nbCol-2). Two
    // cheap filters cut down the number of rank tests: R must carry at least
    // nbCol-2 rays, and -- since every current edge of fa is a genuine ridge and
    // all ridges share the same dimension -- if R is comparable (subset or
    // superset) to any existing ridge of fa it cannot be a new ridge, so it is
    // skipped. Only candidates incomparable to all of fa's ridges are verified.
    std::vector<int> touch = new_facets;
    touch.insert(touch.end(), incident.begin(), incident.end());
    for (size_t a = 0; a < touch.size(); a++) {
      int fa = touch[a];
      for (size_t b = a + 1; b < touch.size(); b++) {
        int fb = touch[b];
        if (edge_exists(fa, fb))
          continue;
        Face R = F[fa].incd & F[fb].incd;
        if (static_cast<int>(R.count()) < nbCol - 2)
          continue;
        bool comparable = false;
        for (int nbr : F[fa].adj) {
          if (beneath_beyond::face_incl(F[fa].incd & F[nbr].incd, R) != 2) {
            comparable = true;
            break;
          }
        }
        if (comparable)
          continue;
        if (beneath_beyond::face_rank(EXT, R) == nbCol - 2)
          add_edge(fa, fb);
      }
    }
  }

  std::vector<BeneathBeyondFacet<T>> result;
  for (auto &f : F)
    if (f.alive)
      result.push_back({std::move(f.normal), std::move(f.incd)});
#ifdef TIMINGS_BENEATH_BEYOND
  os << "BENEATH_BEYOND(DG): |EXT|=" << nbRow << "/" << nbCol
     << " |facets|=" << result.size() << " time=" << time << "\n";
#endif
  return result;
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

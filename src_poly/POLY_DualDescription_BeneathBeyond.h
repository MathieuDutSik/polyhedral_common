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
// Two kernels are provided, producing identical output (cross-checked):
//  * BeneathBeyond_Kernel: the straightforward version described just below --
//    it rediscovers the horizon by testing visible x positive facet pairs and
//    stores each facet's incidence by a full scan. Kept as a simple correctness
//    oracle and used by the tests.
//  * BeneathBeyond_Kernel_DualGraph: a port of polymake's update_facets that
//    maintains a facet adjacency graph and carries incidence forward
//    incrementally, avoiding the per-facet scan and the pairwise ridge search.
//    It is the default behind the public entry points below because it is one
//    to three orders of magnitude faster on the degenerate polytopes
//    (cut / metric / TSP) that dominate this codebase.
//
// The algorithm uses no field division: the flip formula and the gcd
// canonicalization are ring operations, and the only nullspace (for the nbCol
// initial facets) is taken over the overlying field and scaled back. So the
// kernel runs on a ring as well as a field, and for a field input the public
// entry points run it on the underlying integer ring (typically ~2x faster than
// the rational arithmetic) and convert the integer normals back.
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
// A facet node of the dual graph. `incd` is the incidence over the points
// inserted so far, maintained incrementally (a new facet starts as ridge u {p},
// an incident facet absorbs p), so no EXT scan is ever performed. Once every
// point has been inserted, `incd` is the full incidence.
template <typename T> struct DGFacet {
  MyVector<T> normal;
  Face incd;
  std::vector<int> adj;
  bool alive;
};

// Dual-graph beneath-and-beyond, ported from polymake's update_facets. It carries
// the incidence forward incrementally and reconnects new facets by a purely
// combinatorial (set-containment) ridge test, so -- unlike BeneathBeyond_Kernel --
// it does no per-facet incidence scan and no rank test in the main loop. It
// returns the same facets as BeneathBeyond_Kernel (cross-checked). Field type,
// full-dimensional pointed cone only.
//
// Per inserted ray p:
//   * classify the live facets by sign(normal . p): visible (<0), incident (=0);
//   * absorb p into the incident facets (incd += p);
//   * for every dual-graph edge from a visible facet to a strictly positive one
//     (a horizon ridge) create a new facet with incidence ridge u {p} and normal
//     from the flip formula; coplanar new facets (equal normal) are merged;
//   * delete the visible facets;
//   * reconnect the touched facets (new + incident) with the incl-maximality test
//     on their incidence sets -- the crucial point being that these sets contain
//     only already-processed points, which is what makes the combinatorial test
//     exact.
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
  auto add_facet = [&](MyVector<T> normal, Face incd) -> int {
    int id = F.size();
    F.push_back(DGF{std::move(normal), std::move(incd), {}, true});
    return id;
  };
  auto edge_exists = [&](int a, int b) -> bool {
    for (int x : F[a].adj)
      if (x == b)
        return true;
    return false;
  };
  auto add_edge = [&](int a, int b) {
    if (!edge_exists(a, b)) {
      F[a].adj.push_back(b);
      F[b].adj.push_back(a);
    }
  };
  auto remove_edge = [&](int a, int b) {
    auto rm = [&](int x, int y) {
      auto &v = F[x].adj;
      v.erase(std::remove(v.begin(), v.end(), y), v.end());
    };
    rm(a, b);
    rm(b, a);
  };
  // Scratch used by the reconnection to answer "is fb already a neighbour of the
  // current fa?" in O(1): nbr_tok[x] holds the token stamped when x was last a
  // neighbour of the fa being processed. Monotonic tokens avoid clearing it.
  std::vector<long> nbr_tok;
  long cur_tok = 0;
  // Reused scratch bitsets so the reconnection's set intersections are done
  // in place (boost::dynamic_bitset's operator& allocates a temporary each call,
  // which dominated the reconnection over hundreds of millions of pairs).
  Face sR(nbRow), sE(nbRow);
  // Canonicalize an inward normal and restore the orientation of `raw` (which
  // canonicalization may reverse for non-rational fields).
  auto canon = [&](MyVector<T> const &raw) -> MyVector<T> {
    MyVector<T> normal = ScalarCanonicalizationVector(raw);
    for (int iCol = 0; iCol < nbCol; iCol++)
      if (raw(iCol) != 0) {
        if ((normal(iCol) > 0) != (raw(iCol) > 0))
          normal = -normal;
        break;
      }
    return normal;
  };

  // Initial simplicial cone: nbCol facets (incidence = the nbCol-1 basis rays
  // spanning each), forming a complete dual graph. The initial facet normals are
  // the only place a nullspace (a field operation) is needed; when T is a ring we
  // compute them over the overlying field and scale the result back to T, so the
  // whole rest of the algorithm runs in ring arithmetic.
  using Tfield = typename overlying_field<T>::field_type;
  [[maybe_unused]] MyMatrix<Tfield> EXTfield;
  if constexpr (!is_ring_field<T>::value)
    EXTfield = UniversalMatrixConversion<Tfield, T>(EXT);
  std::vector<int> init;
  for (int iBas = 0; iBas < nbCol; iBas++) {
    Face seed(nbRow);
    for (int jBas = 0; jBas < nbCol; jBas++)
      if (jBas != iBas)
        seed[basis[jBas]] = 1;
    MyVector<T> normal;
    if constexpr (is_ring_field<T>::value) {
      normal = ScalarCanonicalizationVector(FindFacetInequality(EXT, seed));
    } else {
      MyVector<Tfield> nf = FindFacetInequality(EXTfield, seed);
      normal = ScalarCanonicalizationVector(
          UniversalVectorConversion<T, Tfield>(NonUniqueScaleToIntegerVector(nf)));
    }
    if (beneath_beyond::facet_scal(EXT, normal, basis[iBas]) < 0)
      normal = -normal;
    init.push_back(add_facet(std::move(normal), seed));
  }
  for (size_t i = 0; i < init.size(); i++)
    for (size_t j = i + 1; j < init.size(); j++)
      add_edge(init[i], init[j]);

  for (int p = 0; p < nbRow; p++) {
    if (in_basis[p])
      continue;
    // Classify the live facets: sign of normal . p.
    std::vector<T> sval(F.size());
    std::vector<int> orient(F.size(), 2);
    std::vector<int> visible, incident;
    for (size_t i = 0; i < F.size(); i++) {
      if (!F[i].alive)
        continue;
      T s = beneath_beyond::facet_scal(EXT, F[i].normal, p);
      sval[i] = s;
      if (s < 0) {
        orient[i] = -1;
        visible.push_back(static_cast<int>(i));
      } else if (s > 0) {
        orient[i] = 1;
      } else {
        orient[i] = 0;
        incident.push_back(static_cast<int>(i));
      }
    }
    if (visible.empty())
      continue;

    // p lies on every incident facet: absorb it.
    for (int fi : incident)
      F[fi].incd[p] = 1;

    // Horizon: each edge from a visible facet to a strictly positive one is a
    // ridge; ridge u {p} is a new facet, its normal the flip combination of the
    // two parents. Deduplicate by normal against the incident facets and among
    // the new facets, merging incidences of coplanar facets.
    std::map<MyVector<T>, int> norm_to_facet;
    for (int fi : incident)
      norm_to_facet.emplace(F[fi].normal, fi);
    std::vector<int> new_facets;
    std::vector<std::pair<int, int>> nf_parent;
    for (int iv : visible) {
      std::vector<int> nbrs = F[iv].adj; // copy: add_facet may reallocate F
      for (int g : nbrs) {
        if (!F[g].alive || orient[g] <= 0)
          continue;
        Face incd = F[iv].incd & F[g].incd;
        incd[p] = 1;
        MyVector<T> raw = sval[g] * F[iv].normal - sval[iv] * F[g].normal;
        MyVector<T> normal = canon(raw);
        auto it = norm_to_facet.find(normal);
        int nf;
        if (it == norm_to_facet.end()) {
          nf = add_facet(normal, incd);
          norm_to_facet.emplace(std::move(normal), nf);
          new_facets.push_back(nf);
        } else {
          nf = it->second;
          F[nf].incd |= incd; // coplanar facet: merge incidences
        }
        nf_parent.push_back({nf, g});
      }
    }
    for (auto &pr : nf_parent)
      add_edge(pr.first, pr.second);

    // Delete the visible facets.
    for (int iv : visible) {
      for (int g : F[iv].adj) {
        auto &ga = F[g].adj;
        ga.erase(std::remove(ga.begin(), ga.end(), iv), ga.end());
      }
      F[iv].adj.clear();
      F[iv].alive = false;
    }

    // Reconnect the touched facets (new + incident). A pair (fa, fb) becomes an
    // edge iff R = fa ∩ fb is a maximal proper common face (a ridge). Maximality
    // is decided combinatorially (polymake's incl) against fa's current ridges:
    // an existing ridge contained in R is not maximal and is dropped; if R is
    // contained in an existing ridge it is not maximal and rejected.
    std::vector<int> touch = new_facets;
    touch.insert(touch.end(), incident.begin(), incident.end());
    int min_ridge = nbCol - 2;
    if (nbr_tok.size() < F.size())
      nbr_tok.resize(F.size(), 0);
    for (size_t a = 0; a < touch.size(); a++) {
      int fa = touch[a];
      // Stamp fa's current neighbours so "does edge (fa,fb) exist?" is O(1).
      ++cur_tok;
      for (int x : F[fa].adj)
        nbr_tok[x] = cur_tok;
      Face const &fa_incd = F[fa].incd;
      for (size_t b = a + 1; b < touch.size(); b++) {
        int fb = touch[b];
        if (nbr_tok[fb] == cur_tok)
          continue;
        sR = fa_incd;
        sR &= F[fb].incd; // sR = R = fa ∩ fb
        long cR = sR.count();
        if (cR < min_ridge)
          continue;
        bool add = true;
        std::vector<int> nbrs = F[fa].adj; // copy: edges may be removed below
        for (int nbr : nbrs) {
          // Relation of the existing ridge E = fa ∩ nbr to the candidate R,
          // from set sizes (E ⊆ R iff |E∩R|==|E|, R ⊆ E iff |E∩R|==|R|).
          sE = fa_incd;
          sE &= F[nbr].incd; // sE = E
          long cE = sE.count();
          sE &= sR;              // sE = E ∩ R
          long cER = sE.count();
          bool E_sub_R = (cER == cE);
          bool R_sub_E = (cER == cR);
          if (!E_sub_R && !R_sub_E)
            continue; // incomparable
          if (E_sub_R) {
            remove_edge(fa, nbr);
            nbr_tok[nbr] = 0; // fa lost this neighbour
          }
          if (R_sub_E) {
            add = false;
            break;
          }
        }
        if (add) {
          // fb is known to be absent from fa's neighbours (nbr_tok check above).
          F[fa].adj.push_back(fb);
          F[fb].adj.push_back(fa);
          nbr_tok[fb] = cur_tok; // fa gained this neighbour
        }
      }
    }
  }

  // The incrementally-carried incidence lists the extreme rays on each facet,
  // which is exactly the facet incidence under this header's extreme-ray
  // precondition (verified against the full-scan reference kernel). Redundant
  // rays -- which also lie on a facet -- are outside that precondition and would
  // need BeneathBeyond_Kernel's full scan; here we keep the incremental incidence
  // and avoid rescanning every facet.
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
} // namespace beneath_beyond

// Runs the dual-graph kernel and returns facets with T-typed normals. The
// beneath-and-beyond arithmetic uses no field division (the flip formula and the
// gcd canonicalization are ring operations), so when T is a field the whole
// computation is done on the underlying integer ring -- typically ~3x faster
// than the rational arithmetic -- with the resulting integer normals converted
// back to T. Facet incidences are combinatorial and type independent.
template <typename T>
std::vector<BeneathBeyondFacet<T>>
BeneathBeyond_run(MyMatrix<T> const &EXT, std::ostream &os) {
  if constexpr (is_ring_field<T>::value) {
    using Tring = typename underlying_ring<T>::ring_type;
    int nbRow = EXT.rows();
    int nbCol = EXT.cols();
    MyMatrix<Tring> EXTring(nbRow, nbCol);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      MyVector<T> eRow = NonUniqueScaleToIntegerVector(GetMatrixRow(EXT, iRow));
      AssignMatrixRow(EXTring, iRow, UniversalVectorConversion<Tring, T>(eRow));
    }
    std::vector<BeneathBeyondFacet<Tring>> ring_facets =
        beneath_beyond::BeneathBeyond_Kernel_DualGraph(EXTring, os);
    std::vector<BeneathBeyondFacet<T>> result;
    result.reserve(ring_facets.size());
    for (auto &f : ring_facets)
      result.push_back(
          {UniversalVectorConversion<T, Tring>(f.normal), std::move(f.incd)});
    return result;
  } else {
    return beneath_beyond::BeneathBeyond_Kernel_DualGraph(EXT, os);
  }
}

// Facet incidences, matching DirectFacetComputationIncidence.
template <typename T>
vectface POLY_DualDescription_BeneathBeyondIncidence(MyMatrix<T> const &EXT,
                                                     std::ostream &os) {
  std::vector<BeneathBeyondFacet<T>> facets = BeneathBeyond_run(EXT, os);
  vectface vf(EXT.rows());
  for (auto &facet : facets)
    vf.push_back(facet.incd);
  return vf;
}

// Facet normals, matching DirectFacetComputationInequalities.
template <typename T>
MyMatrix<T> POLY_DualDescription_BeneathBeyondInequalities(MyMatrix<T> const &EXT,
                                                           std::ostream &os) {
  std::vector<BeneathBeyondFacet<T>> facets = BeneathBeyond_run(EXT, os);
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
  std::vector<BeneathBeyondFacet<T>> facets = BeneathBeyond_run(EXT, os);
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

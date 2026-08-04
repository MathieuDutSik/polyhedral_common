// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DUALDESC_BENEATH_AND_BEYOND_H_
#define SRC_POLY_POLY_DUALDESC_BENEATH_AND_BEYOND_H_

// clang-format off
#include "POLY_Fundamental.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include <algorithm>
#include <bit>
#include <cstdint>
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

// A row-set packed into 64-bit blocks. It exists because boost::dynamic_bitset
// (our Face) offers no block access, so its operator& / operator&= / count()
// carry per-call overhead that dominates when a hot loop does millions of
// intersection-size and subset tests on small, fixed-width row sets -- exactly
// the beneath-and-beyond reconnection. The block loops below inline to a couple
// of `and`+`popcount` instructions. This is deliberately generic (nothing
// beneath-and-beyond specific) so it can be reused wherever such inner loops
// show up.
struct PackedSet {
  std::vector<uint64_t> w;
  PackedSet() = default;
  explicit PackedSet(int nbBit) : w((nbBit + 63) / 64, 0) {}
  void set(int i) { w[i >> 6] |= uint64_t(1) << (i & 63); }
  void or_with(PackedSet const &o) {
    for (size_t i = 0; i < w.size(); i++)
      w[i] |= o.w[i];
  }
  long count() const {
    long c = 0;
    for (uint64_t x : w)
      c += std::popcount(x);
    return c;
  }
};

// |a intersect b|, without materializing the intersection.
inline long inter_count(PackedSet const &a, PackedSet const &b) {
  long c = 0;
  for (size_t i = 0; i < a.w.size(); i++)
    c += std::popcount(a.w[i] & b.w[i]);
  return c;
}

// out = a intersect b (out already sized like a and b).
inline void inter_into(PackedSet const &a, PackedSet const &b, PackedSet &out) {
  for (size_t i = 0; i < a.w.size(); i++)
    out.w[i] = a.w[i] & b.w[i];
}

// Conversions to/from the boost Face used at the algorithm's boundaries.
inline PackedSet packed_from_face(Face const &f) {
  PackedSet s(f.size());
  boost::dynamic_bitset<>::size_type i = f.find_first();
  while (i != boost::dynamic_bitset<>::npos) {
    s.set(static_cast<int>(i));
    i = f.find_next(i);
  }
  return s;
}
inline Face face_from_packed(PackedSet const &s, int nbBit) {
  Face f(nbBit);
  for (int i = 0; i < nbBit; i++)
    if ((s.w[i >> 6] >> (i & 63)) & 1)
      f[i] = 1;
  return f;
}

// Scalar product normal . EXT.row(iRow).
template <typename T>
inline T facet_scal(MyMatrix<T> const &EXT, MyVector<T> const &normal,
                    int iRow) {
  int nbCol = EXT.cols();
  T sum(0);
  for (int iCol = 0; iCol < nbCol; iCol++)
    AddMul(sum, normal(iCol), EXT(iRow, iCol));
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
  // linear combination of two existing ones (see below). The subset solver is
  // built once for the nbCol seeds and delivers the kernel vector together with
  // its incidence in one pass. Orient each inward, so it is positive on the
  // dropped basis ray.
  {
    SubsetRankOneSolver<T> solver(EXT);
    for (int iBas = 0; iBas < nbCol; iBas++) {
      Face seed(nbRow);
      for (int jBas = 0; jBas < nbCol; jBas++)
        if (jBas != iBas)
          seed[basis[jBas]] = 1;
      std::pair<MyVector<T>, Face> pair =
          solver.GetPositiveKernelVectorAndFace(seed);
      MyVector<T> normal = ScalarCanonicalizationVector(pair.first);
      if (beneath_beyond::facet_scal(EXT, normal, basis[iBas]) < 0)
        normal = -normal;
      facets.push_back({std::move(normal), std::move(pair.second)});
    }
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
  PackedSet incd; // incidence over the inserted-so-far rays, packed for speed
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
  auto add_facet = [&](MyVector<T> normal, PackedSet incd) -> int {
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
  // Reused scratch sets so the reconnection's intersections are materialized in
  // place, only when actually needed (never per rejected pair).
  PackedSet sR(nbRow), sE(nbRow);
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
  // the only place a kernel computation is needed; the subset solver is built
  // once for the nbCol seeds and handles the arithmetic dispatch, so the whole
  // rest of the algorithm runs in ring arithmetic.
  SubsetRankOneSolver<T> solver(EXT);
  std::vector<int> init;
  for (int iBas = 0; iBas < nbCol; iBas++) {
    Face seed(nbRow);
    for (int jBas = 0; jBas < nbCol; jBas++)
      if (jBas != iBas)
        seed[basis[jBas]] = 1;
    MyVector<T> normal =
        ScalarCanonicalizationVector(solver.GetPositiveKernelVector(seed));
    if (beneath_beyond::facet_scal(EXT, normal, basis[iBas]) < 0)
      normal = -normal;
    init.push_back(
        add_facet(std::move(normal), beneath_beyond::packed_from_face(seed)));
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
      F[fi].incd.set(p);

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
        PackedSet incd(nbRow);
        beneath_beyond::inter_into(F[iv].incd, F[g].incd, incd);
        incd.set(p);
        MyVector<T> raw = sval[g] * F[iv].normal - sval[iv] * F[g].normal;
        MyVector<T> normal = canon(raw);
        auto it = norm_to_facet.find(normal);
        int nf;
        if (it == norm_to_facet.end()) {
          nf = add_facet(normal, std::move(incd));
          norm_to_facet.emplace(std::move(normal), nf);
          new_facets.push_back(nf);
        } else {
          nf = it->second;
          F[nf].incd.or_with(incd); // coplanar facet: merge incidences
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
      PackedSet const &fa_incd = F[fa].incd;
      for (size_t b = a + 1; b < touch.size(); b++) {
        int fb = touch[b];
        if (nbr_tok[fb] == cur_tok)
          continue;
        // Size filter first, without materializing R -- most pairs fail here.
        long cR = beneath_beyond::inter_count(fa_incd, F[fb].incd);
        if (cR < min_ridge)
          continue;
        beneath_beyond::inter_into(fa_incd, F[fb].incd, sR); // sR = R
        bool add = true;
        std::vector<int> nbrs = F[fa].adj; // copy: edges may be removed below
        for (int nbr : nbrs) {
          // Relation of the existing ridge E = fa ∩ nbr to the candidate R,
          // from set sizes (E ⊆ R iff |E∩R|==|E|, R ⊆ E iff |E∩R|==|R|).
          beneath_beyond::inter_into(fa_incd, F[nbr].incd, sE); // sE = E
          long cE = sE.count();
          long cER = beneath_beyond::inter_count(sE, sR); // |E ∩ R|
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
      result.push_back({std::move(f.normal),
                        beneath_beyond::face_from_packed(f.incd, nbRow)});
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

// Core placing-triangulation enumeration, running entirely in the working type
// Twork on EXTwork. For each simplex, as soon as it is created, it calls
//   f_core(simplex, det_work)
// where simplex is the sorted list of its nbCol vertex row indices and det_work
// is the SIGNED determinant of SelectRow(EXTwork, simplex) in Twork. Twork may
// be an integer-scaled copy of a rational input (see the field wrapper below),
// in which case det_work is off from the true determinant by the product of the
// per-row scales; the caller corrects for that. All arithmetic here is ring
// arithmetic (determinants, sign tests, gcd), so Twork can be the underlying
// ring even when the original type is a field.
//
// The determinant comes essentially for free: the visibility test that decides
// whether ray p cones over a boundary cell already computes dp, the signed
// determinant of the rows [face..., p], and that set of rows *is* the new
// simplex. The only determinant not produced by the enumeration is that of the
// initial simplex, fixed once by a single DeterminantMat (negligible).
//
// The triangulation itself is the direct/robust placing triangulation (Lee):
// start from a simplex on d = nbCol independent rays; for each further ray p,
// cone p over every boundary (d-1)-cell of the current triangulation that is
// strictly visible from p (p on the opposite side of the cell's hyperplane from
// the cell's apex). Redundant rays see no boundary cell and are skipped
// automatically. The placing triangulation is append-only (existing simplices
// are never subdivided or removed), which is exactly what makes emit-on-creation
// valid. The boundary cells are recomputed each step (a (d-1)-face is on the
// boundary iff it belongs to exactly one simplex), which is simpler and less
// error-prone than threading the triangulation through the facet update, at the
// cost of some speed.
template <typename Twork, typename Fcore>
void BeneathBeyond_TriangulationDet_core(MyMatrix<Twork> const &EXTwork,
                                         [[maybe_unused]] std::ostream &os,
                                         Fcore f_core) {
  int nbRow = EXTwork.rows();
  int nbCol = EXTwork.cols();
#ifdef TIMINGS_BENEATH_BEYOND
  MicrosecondTime time;
#endif
#ifdef SANITY_CHECK_BENEATH_BEYOND
  int rank = RankMat(EXTwork);
  if (rank != nbCol) {
    std::cerr << "BeneathBeyond(triang): EXT must span a full-dimensional cone, "
              << "RankMat(EXT)=" << rank << " but nbCol=" << nbCol << "\n";
    throw TerminalException{1};
  }
#endif

  // Signed d x d determinant of the rows (face vertices, then extra), where the
  // face keeps its stored order so signs are comparable across calls.
  auto det_face_point = [&](std::vector<int> const &face, int extra) -> Twork {
    MyMatrix<Twork> M(nbCol, nbCol);
    for (int i = 0; i < nbCol - 1; i++)
      for (int j = 0; j < nbCol; j++)
        M(i, j) = EXTwork(face[i], j);
    for (int j = 0; j < nbCol; j++)
      M(nbCol - 1, j) = EXTwork(extra, j);
    return DeterminantMat(M);
  };

  std::vector<int> basis = TMat_ListRowSelect(EXTwork);
  Face in_basis(nbRow);
  for (auto &eRow : basis)
    in_basis[eRow] = 1;

  // Producer working state: only the vertex lists are needed to recompute the
  // boundary; determinants are streamed to f_core, never retained here.
  std::vector<std::vector<int>> simplices;
  [[maybe_unused]] size_t n_simplices = 0;
  {
    std::vector<int> simplex = basis;
    std::sort(simplex.begin(), simplex.end());
    // The initial simplex is the only determinant the visibility test does not
    // yield; one DeterminantMat fixes it, in sorted-row order like the rest.
    Twork det = DeterminantMat(SelectRow(EXTwork, simplex));
    f_core(simplex, det);
    n_simplices++;
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
      Twork dp = det_face_point(face, p);
      if (dp == 0)
        continue; // p coplanar with the face: not strictly visible
      Twork da = det_face_point(face, apex);
      // p is visible from outside iff on the opposite side from the apex. The
      // sign of a determinant is invariant under the positive per-row scaling
      // that produces EXTwork, so this test is identical to the field one.
      if ((dp > 0) == (da > 0))
        continue;
      // dp is the signed determinant of the rows [face..., p]. The stored
      // simplex is the same rows in sorted index order; since face is already
      // ascending (it is a std::map key of ascending indices), sorting only
      // moves p left past every face vertex greater than it, so the sign of the
      // sort permutation is (-1)^{#{f in face : f > p}}.
      std::vector<int> simplex = face;
      simplex.push_back(p);
      std::sort(simplex.begin(), simplex.end());
      int n_inv = 0;
      for (int f : face)
        if (f > p)
          n_inv++;
      Twork det = ((n_inv % 2) == 0) ? dp : -dp;
#ifdef SANITY_CHECK_BENEATH_BEYOND
      Twork det_check = DeterminantMat(SelectRow(EXTwork, simplex));
      if (det_check != det) {
        std::cerr << "BeneathBeyond(triang): signed-determinant mismatch, "
                  << "recovered=" << det << " actual=" << det_check << "\n";
        throw TerminalException{1};
      }
#endif
      f_core(simplex, det);
      n_simplices++;
      new_simplices.push_back(std::move(simplex));
    }
    for (auto &e : new_simplices)
      simplices.push_back(std::move(e));
  }

#ifdef TIMINGS_BENEATH_BEYOND
  os << "BENEATH_BEYOND: triangulationDet |EXT|=" << nbRow << "/" << nbCol
     << " |simplices|=" << n_simplices << " time=" << time << "\n";
#endif
}

// Placing triangulation produced by the beneath-and-beyond insertion order.
// Streaming variant: for each simplex, as soon as it is created, this calls
//   f_trig_det(simplex, det)
// where
//  * simplex is the sorted list of its nbCol vertex row indices, and
//  * det is the SIGNED determinant of SelectRow(EXT, simplex), in type T.
// This is the callback analogue of lrs::GetTriangulationDet_f: the caller's
// consumer never has to hold the full simplex list. Note, however, that unlike
// lrs' reverse search the *producer* still keeps the growing triangulation as
// working state (needed to recompute the boundary at each step), so this bounds
// the consumer's memory, not the producer's -- see the discussion in
// POLY_DualDesc_lrslib.h on why beneath-and-beyond is not a bounded-memory tree search.
//
// Ring/field split (same idiom as lrs::DualDescription): the enumeration uses
// only ring operations, so when T is a field each ray is scaled to an integer
// vector and the whole placing triangulation runs on the underlying ring
// (typically ~3x faster). The determinant then needs the care the field version
// hides: row i of the ring matrix is scale[i] * (row i of EXT), so the ring
// determinant of a simplex is (prod_i scale[i]) times the true one. We keep the
// per-row scales (the cleared denominators) and divide them back out to recover
// the exact rational determinant. The scales are positive, so simplex identity
// and orientation are untouched -- only the determinant magnitude is corrected.
template <typename T, typename Ftrig_det>
void POLY_DualDescription_BeneathBeyondTriangulationDet_f(
    MyMatrix<T> const &EXT, std::ostream &os, Ftrig_det f_trig_det) {
  // The ring path applies to exact fields only; floating fields (double/float)
  // have no exact underlying ring and no denominators to clear, so they run
  // directly.
  if constexpr (is_ring_field<T>::value && !std::is_floating_point_v<T>) {
    using Tring = typename underlying_ring<T>::ring_type;
    int nbRow = EXT.rows();
    int nbCol = EXT.cols();
    MyMatrix<Tring> EXTring(nbRow, nbCol);
    std::vector<T> scale(nbRow); // per-row cleared denominator, positive
    for (int iRow = 0; iRow < nbRow; iRow++) {
      FractionVector<T> fr =
          NonUniqueScaleToIntegerVectorPlusCoeff(GetMatrixRow(EXT, iRow));
      scale[iRow] = fr.TheMult; // EXTring row = fr.TheMult * (EXT row)
      AssignMatrixRow(EXTring, iRow,
                      UniversalVectorConversion<Tring, T>(fr.TheVect));
    }
    auto f_core = [&](std::vector<int> const &simplex,
                      Tring const &det_ring) -> void {
      // det(EXT rows) = det(EXTring rows) / prod scale[i] over the simplex.
      T denom(1);
      for (int idx : simplex)
        denom *= scale[idx];
      T det = UniversalScalarConversion<T, Tring>(det_ring) / denom;
      f_trig_det(simplex, det);
    };
    BeneathBeyond_TriangulationDet_core(EXTring, os, f_core);
  } else {
    BeneathBeyond_TriangulationDet_core(EXT, os, f_trig_det);
  }
}

// Placing triangulation with, for each simplex, the SIGNED determinant of its
// vertex matrix SelectRow(EXT, simplex) -- the same (simplex, det) format as
// lrs::GetTriangulationDet, so the two are interchangeable. Thin wrapper that
// collects the streamed simplices into a vector.
template <typename T>
std::vector<std::pair<std::vector<int>, T>>
POLY_DualDescription_BeneathBeyondTriangulationDet(MyMatrix<T> const &EXT,
                                                   std::ostream &os) {
  std::vector<std::pair<std::vector<int>, T>> l_trig;
  auto f_trig_det = [&](std::vector<int> const &trig, T const &det) -> void {
    l_trig.push_back({trig, det});
  };
  POLY_DualDescription_BeneathBeyondTriangulationDet_f<T>(EXT, os, f_trig_det);
  return l_trig;
}

// Placing triangulation as a bare list of simplices, each the sorted list of
// its nbCol vertex row indices (same format as lrs::GetTriangulation). Thin
// wrapper that collects the streamed simplices, dropping the determinants.
template <typename T>
std::vector<std::vector<int>>
POLY_DualDescription_BeneathBeyondTriangulation(MyMatrix<T> const &EXT,
                                                std::ostream &os) {
  std::vector<std::vector<int>> simplices;
  auto f_trig_det = [&](std::vector<int> const &trig,
                        [[maybe_unused]] T const &det) -> void {
    simplices.push_back(trig);
  };
  POLY_DualDescription_BeneathBeyondTriangulationDet_f<T>(EXT, os, f_trig_det);
  return simplices;
}

// clang-format off
#endif  // SRC_POLY_POLY_DUALDESC_BENEATH_AND_BEYOND_H_
// clang-format on

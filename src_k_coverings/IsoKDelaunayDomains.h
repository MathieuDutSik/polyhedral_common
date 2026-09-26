// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_K_COVERINGS_ISOKDELAUNAYDOMAINS_H_
#define SRC_K_COVERINGS_ISOKDELAUNAYDOMAINS_H_

// clang-format off
#include "LatticeKDelaunay.h"
#include "IsoDelaunayDomains.h"
#include "CoveringMaxdet.h"
#include <algorithm>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ISO_K_DELAUNAY
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ISO_K_DELAUNAY
#endif

#ifdef TIMINGS
#define TIMINGS_ISO_K_DELAUNAY
#endif

/*
  Iso-k-Delaunay domains, or (L,k)-types.

  For a T-space of quadratic forms, the (L,k)-type of a positive definite
  Q of the space is the combinatorial type of its order-k Delaunay tiling
  (LatticeKDelaunay.h): the list of the tiles (P_-, P_0) up to the affine
  lattice isometries. The set of forms of the space with a given (L,k)-type
  is a polyhedral cone, whose top-dimensional cells are enumerated here.

  The lifting heights of the point set Q_k of Aurenhammer,
  sum_{p in S} Q[p], are linear in Q, so a tiling stays the lower hull of
  the lifted set on a polyhedral cone of forms, exactly as for the
  Delaunay subdivisions (k = 1). The cone of a generic tiling is given by
  two families of linear inequalities, one per tile (P_-, P_0) with
  circumcenter v and radius r:

  (a) For each tile adjacent across a facet, the lifted vertices of the
      adjacent tile are above the hyperplane carrying the lift of the tile.
      With phi(p) = |p - v|_Q^2 - r^2 expressed through the affine
      coordinates of p in an affine basis of P_0 (the Voronoi regulator
      phi(p) = Q[p] - sum lambda_i Q[b_i], linear in Q), the condition for
      a k-subset S' of the adjacent tile off the shared facet reads

          sum_{p in S'} phi(p) - sum_{p in P_-} phi(p) >= 0.

      All the vertices off the facet give proportional inequalities.

  (b) For each p in P_-, the point stays inside the sphere: phi(p) <= 0.

  The walls of the cone are the forms where the tiling changes. A change
  where some tile grows by the merging of adjacent tiles is caught by (a).
  The remaining changes are those where the sphere of a tile acquires a
  new point without the tile growing: an outside point cannot do it (the
  hypersimplex would have to contain the sum of that point with j-1 points
  of P_0 for every choice of the j-1 points, impossible for a point on the
  sphere), but a point of P_- can, when j = n and the point lies in the
  simplex P_0 (the stellar subdivision of the tile appears on the other
  side of the wall). Those are caught by (b). For k = 1 the family (b) is
  empty and the cone is the usual L-type domain.

  The enumeration mirrors the one of the iso-Delaunay domains: a generic
  form of the space seeds the enumeration, the facets of the cone are
  computed up to the stabilizer of the domain, and the domain across a
  facet is obtained by recomputing the tiling at a form slightly beyond an
  interior point of the facet and checking that the facet point lies in
  the closure of the resulting domain. The equivalence and the invariants
  of domains are those of their canonical interior forms, through the
  invariant vector families as for the iso-Delaunay domains.
 */

template <typename T, typename Tint, typename Tgroup> struct IsoKDelaunayDomain {
  KDelaunayTesselation<Tint, Tgroup> DT;
  // The defining inequalities (a) and (b), as primitive integral vectors in
  // the coordinates of the T-space, made unique.
  MyMatrix<Tint> ListIneq;
  // The canonical interior form (primitive integral) and its invariant
  // vector family.
  MyMatrix<Tint> GramMat;
  MyMatrix<Tint> SHV;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, IsoKDelaunayDomain<T, Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("DT", eRec.DT);
  ar &make_nvp("ListIneq", eRec.ListIneq);
  ar &make_nvp("GramMat", eRec.GramMat);
  ar &make_nvp("SHV", eRec.SHV);
}
} // namespace boost::serialization

// Whether some tile of the tiling is co-spherical only on a wall of the
// T-space, which makes the domain lower dimensional.
template <typename Tint, typename Tgroup>
bool IsKDelaunayTesselationInducingEqualities(
    KDelaunayTesselation<Tint, Tgroup> const &DT,
    std::vector<std::vector<Tint>> const &ListGramRing, std::ostream &os) {
  for (auto &eEnt : DT.l_tiles) {
    if (IsDelaunayPolytopeInducingEqualities(eEnt.tile.EXT, ListGramRing, os)) {
      return true;
    }
  }
  return false;
}

/*
  The defining inequalities (a) and (b) of the domain of a generic tiling,
  over the ring: the Voronoi regulators are computed by
  VoronoiLinearInequality (multiplied by the positive determinant of the
  affine basis, which the canonicalization absorbs) and combined
  linearly.
 */
template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint> FillKDelaunayInequalities(
    KDelaunayTesselation<Tint, Tgroup> &DT,
    std::vector<std::vector<Tint>> const &ListGramRing, std::ostream &os) {
  int k = DT.k;
  int dimSpace = ListGramRing.size();
  std::unordered_set<MyVector<Tint>> set_ineq;
  std::vector<MyVector<Tint>> l_ineq;
  auto insert = [&](MyVector<Tint> const &V) -> MyVector<Tint> {
#ifdef SANITY_CHECK_ISO_K_DELAUNAY
    if (IsZeroVector(V)) {
      std::cerr << "ISO_K_DELAUNAY: a defining inequality is zero, the form "
                   "is not generic in the T-space\n";
      throw TerminalException{1};
    }
#endif
    MyVector<Tint> Vcan = ScalarCanonicalizationVector(V);
    if (set_ineq.insert(Vcan).second) {
      l_ineq.push_back(Vcan);
    }
    return Vcan;
  };
  for (auto &eEnt : DT.l_tiles) {
    KDelaunayTile<Tint> const &tile = eEnt.tile;
    int n = tile.EXT.cols() - 1;
    int i = tile.INT.rows();
    // The inequalities already carried by the tiling (the untouched tiles
    // of a flip) are reused: an inequality of a tile depends only on the
    // tile and on its neighbour.
    bool need_int = (static_cast<int>(eEnt.ListIntIneq.size()) != i);
    bool need_adj = false;
    for (auto &eAdj : eEnt.ListAdj) {
      if (eAdj.eIneq.size() == 0) {
        need_adj = true;
      }
    }
    if (!need_int && !need_adj) {
      for (auto &V : eEnt.ListIntIneq) {
        insert(V);
      }
      for (auto &eAdj : eEnt.ListAdj) {
        insert(eAdj.eIneq);
      }
      continue;
    }
    VoronoiInequalityPreComput<Tint> vipc =
        BuildVoronoiIneqPreComputeChecked<Tint>(tile.EXT, ListGramRing, os);
    auto get_phi = [&](MyVector<Tint> const &p_hom) -> MyVector<Tint> {
      return VoronoiLinearInequality(vipc, p_hom, ListGramRing, os);
    };
    MyVector<Tint> sumInt = ZeroVector<Tint>(dimSpace);
    if (need_int) {
      eEnt.ListIntIneq.resize(i);
    }
    for (int i_row = 0; i_row < i; i_row++) {
      MyVector<Tint> p_hom = GetMatrixRow(tile.INT, i_row);
      MyVector<Tint> phi = get_phi(p_hom);
      sumInt += phi;
      // (b): the point stays inside, -phi(p) >= 0.
      MyVector<Tint> V = -phi;
      if (need_int) {
        eEnt.ListIntIneq[i_row] = insert(V);
      } else {
        insert(eEnt.ListIntIneq[i_row]);
      }
    }
    if (!need_adj) {
      for (auto &eAdj : eEnt.ListAdj) {
        insert(eAdj.eIneq);
      }
      continue;
    }
    // (a): one inequality per adjacency, from a k-subset of the adjacent
    // tile whose sum is off the shared facet.
    KDelaunayTileVertices<Tint> vert = GetKDelaunayTileVertices(tile, k);
    MyMatrix<T> EXTsum_T = UniversalMatrixConversion<T, Tint>(vert.EXTsum);
    SubsetRankOneSolver<T> ext_solver(EXTsum_T);
    for (auto &eAdj : eEnt.ListAdj) {
      if (eAdj.eIneq.size() > 0) {
        insert(eAdj.eIneq);
        continue;
      }
      KDelaunayTile<Tint> const &tile2 = DT.l_tiles[eAdj.iOrb].tile;
      MyMatrix<Tint> EXTadj = tile2.EXT * eAdj.eBigMat;
      MyMatrix<Tint> INTadj = tile2.INT * eAdj.eBigMat;
      MyVector<T> ell_hom = ext_solver.GetPositiveKernelVector(eAdj.eInc);
      int i2 = INTadj.rows();
      int m2 = EXTadj.rows();
      int j2 = k - i2;
      MyVector<Tint> sumIntAdj = ZeroVector<Tint>(dimSpace);
      MyVector<Tint> baseAdj = ZeroVector<Tint>(n);
      for (int i_row = 0; i_row < i2; i_row++) {
        MyVector<Tint> p_hom = GetMatrixRow(INTadj, i_row);
        sumIntAdj += get_phi(p_hom);
        for (int u = 0; u < n; u++) {
          baseAdj(u) += INTadj(i_row, u + 1);
        }
      }
      std::vector<MyVector<Tint>> l_phi_ext;
      for (int i_row = 0; i_row < m2; i_row++) {
        MyVector<Tint> p_hom = GetMatrixRow(EXTadj, i_row);
        l_phi_ext.push_back(get_phi(p_hom));
      }
      // The j2-subsets of the adjacent P_0 whose sum is off the facet.
      std::vector<int> subset(j2);
      for (int u = 0; u < j2; u++) {
        subset[u] = u;
      }
      std::optional<MyVector<Tint>> found;
#ifdef SANITY_CHECK_ISO_K_DELAUNAY
      std::unordered_set<MyVector<Tint>> set_can;
#endif
      while (true) {
        MyVector<Tint> sumV = baseAdj;
        for (int u = 0; u < j2; u++) {
          for (int w = 0; w < n; w++) {
            sumV(w) += EXTadj(subset[u], w + 1);
          }
        }
        T scal = ell_hom(0);
        for (int w = 0; w < n; w++) {
          AddMul(scal, ell_hom(w + 1), UniversalScalarConversion<T, Tint>(sumV(w)));
        }
        if (scal != 0) {
          MyVector<Tint> V = sumIntAdj - sumInt;
          for (int u = 0; u < j2; u++) {
            V += l_phi_ext[subset[u]];
          }
#ifdef SANITY_CHECK_ISO_K_DELAUNAY
          // The adjacent tile lies on the far side of the facet, where the
          // positive kernel vector of the facet is negative.
          if (scal > 0) {
            std::cerr << "ISO_K_DELAUNAY: a vertex of the adjacent tile is on "
                         "the wrong side of the facet\n";
            std::cerr << "ISO_K_DELAUNAY: tile.EXT=\n";
            WriteMatrix(std::cerr, tile.EXT);
            std::cerr << "ISO_K_DELAUNAY: tile.INT=\n";
            WriteMatrix(std::cerr, tile.INT);
            std::cerr << "ISO_K_DELAUNAY: vertices=\n";
            WriteMatrix(std::cerr, vert.EXTsum);
            std::cerr << "ISO_K_DELAUNAY: eInc=" << eAdj.eInc
                      << " ell_hom=" << StringVectorGAP(ell_hom) << "\n";
            std::cerr << "ISO_K_DELAUNAY: EXTadj=\n";
            WriteMatrix(std::cerr, EXTadj);
            std::cerr << "ISO_K_DELAUNAY: INTadj=\n";
            WriteMatrix(std::cerr, INTadj);
            std::cerr << "ISO_K_DELAUNAY: sumV=" << StringVectorGAP(sumV)
                      << " scal=" << scal << "\n";
            throw TerminalException{1};
          }
          set_can.insert(ScalarCanonicalizationVector(V));
#endif
          if (!found) {
            found = V;
#ifndef SANITY_CHECK_ISO_K_DELAUNAY
            break;
#endif
          }
        }
        int pos = j2 - 1;
        while (pos >= 0 && subset[pos] == m2 - j2 + pos) {
          pos--;
        }
        if (pos < 0) {
          break;
        }
        subset[pos]++;
        for (int u = pos + 1; u < j2; u++) {
          subset[u] = subset[u - 1] + 1;
        }
      }
      if (!found) {
        std::cerr << "ISO_K_DELAUNAY: no vertex of the adjacent tile is off "
                     "the shared facet\n";
        throw TerminalException{1};
      }
#ifdef SANITY_CHECK_ISO_K_DELAUNAY
      // The vertices of the adjacent tile off the facet all express the
      // same wall, so their inequalities are proportional.
      if (set_can.size() != 1) {
        std::cerr << "ISO_K_DELAUNAY: the vertices of an adjacent tile off "
                     "the facet give " << set_can.size()
                  << " different inequalities instead of one\n";
        throw TerminalException{1};
      }
#endif
      eAdj.eIneq = insert(*found);
    }
  }
  int n_ineq = l_ineq.size();
  MyMatrix<Tint> ListIneq(n_ineq, dimSpace);
  for (int i_ineq = 0; i_ineq < n_ineq; i_ineq++) {
    for (int u = 0; u < dimSpace; u++) {
      ListIneq(i_ineq, u) = l_ineq[i_ineq](u);
    }
  }
  return ListIneq;
}

template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint> ComputeIsoKDelaunayInequalities(
    KDelaunayTesselation<Tint, Tgroup> const &DT,
    std::vector<std::vector<Tint>> const &ListGramRing, std::ostream &os) {
  KDelaunayTesselation<Tint, Tgroup> DTcopy = DT;
  return FillKDelaunayInequalities<T, Tint, Tgroup>(DTcopy, ListGramRing, os);
}

// The order-k tiling of a form of the T-space, or nothing when a tile
// induces equalities on the space (the form is not generic).
template <typename T, typename Tint, typename Tgroup>
std::optional<KDelaunayTesselation<Tint, Tgroup>>
ComputeGenericKDelaunayTesselation(
    MyMatrix<T> const &GramMat, int const &k,
    std::vector<std::vector<Tint>> const &ListGramRing,
    PolyHeuristicSerial<typename Tgroup::Tint> &AllArr, std::ostream &os) {
  DataLattice<T, Tint, Tgroup> data =
      GetDataLattice<T, Tint, Tgroup>(GramMat, AllArr, os);
  auto f_incorrect = [&](KDelaunay_Obj<Tint, Tgroup> const &x) -> bool {
    bool test = IsDelaunayPolytopeInducingEqualities(x.tile.EXT, ListGramRing, os);
#ifdef DEBUG_ISO_K_DELAUNAY
    if (test) {
      os << "ISO_K_DELAUNAY: non-generic tile |INT|=" << x.tile.INT.rows()
         << " |EXT|=" << x.tile.EXT.rows() << " EXT=\n";
      WriteMatrix(os, x.tile.EXT);
      os << "ISO_K_DELAUNAY: INT=\n";
      WriteMatrix(os, x.tile.INT);
    }
#endif
    return test;
  };
  int max_runtime_second = 0;
  return EnumerationKDelaunayTiles<T, Tint, Tgroup, decltype(f_incorrect)>(
      data, k, f_incorrect, max_runtime_second);
}

// The domain of a generic tiling: its inequalities, its canonical interior
// form and the invariant vector family of that form.
template <typename T, typename Tint, typename Tgroup>
IsoKDelaunayDomain<T, Tint, Tgroup>
BuildIsoKDelaunayDomain(KDelaunayTesselation<Tint, Tgroup> DT,
                        [[maybe_unused]] MyMatrix<T> const &GramMat,
                        LinSpaceMatrix<T> const &LinSpa,
                        std::vector<std::vector<Tint>> const &ListGramRing,
                        std::ostream &os,
                        [[maybe_unused]] bool const &check_generating = true) {
#ifdef TIMINGS_ISO_K_DELAUNAY
  MicrosecondTime time_build;
#endif
  MyMatrix<Tint> ListIneq =
      FillKDelaunayInequalities<T, Tint, Tgroup>(DT, ListGramRing, os);
#ifdef TIMINGS_ISO_K_DELAUNAY
  os << "|ISO_K_DELAUNAY: build, inequalities|=" << time_build << "\n";
#endif
#ifdef SANITY_CHECK_ISO_K_DELAUNAY
  // The generating form is interior to its own domain: every inequality is
  // strictly positive at it. Not applicable to a flipped tiling, whose
  // generating form is on the other side of the wall.
  if (check_generating) {
    MyVector<T> c_vec = LINSPA_GetVectorOfMatrixExpression(LinSpa, GramMat);
    int n_ineq = ListIneq.rows();
    int dimSpace = ListIneq.cols();
    for (int i_ineq = 0; i_ineq < n_ineq; i_ineq++) {
      T val(0);
      for (int u = 0; u < dimSpace; u++) {
        AddMul(val, UniversalScalarConversion<T, Tint>(ListIneq(i_ineq, u)),
               c_vec(u));
      }
      if (val <= 0) {
        std::cerr << "ISO_K_DELAUNAY: SANITY_CHECK failed: the inequality "
                  << StringVectorGAP(GetMatrixRow(ListIneq, i_ineq))
                  << " has value " << val << " at the generating form\n";
        std::cerr << "ISO_K_DELAUNAY: GramMat=\n";
        WriteMatrix(std::cerr, GramMat);
        throw TerminalException{1};
      }
    }
  }
#endif
  MyMatrix<T> FAC_T = UniversalMatrixConversion<T, Tint>(ListIneq);
  MyMatrix<T> M = get_interior_gram_matrix_lp(LinSpa, FAC_T, os);
#ifdef TIMINGS_ISO_K_DELAUNAY
  os << "|ISO_K_DELAUNAY: build, interior point|=" << time_build << "\n";
#endif
  MyMatrix<Tint> M_ring = RemoveFractionMatrixPlusCoeffRing(M).TheMat;
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyFullRank<Tint, Tint>(M_ring, os);
#ifdef TIMINGS_ISO_K_DELAUNAY
  os << "|ISO_K_DELAUNAY: build, invariant family|=" << time_build << "\n";
#endif
  return {std::move(DT), std::move(ListIneq), std::move(M_ring), std::move(SHV)};
}

template <typename T, typename Tint, typename Tgroup>
struct DataIsoKDelaunayDomains {
  DataIsoDelaunayDomains<T, Tint, Tgroup> data;
  int k;
  std::vector<std::vector<Tint>> ListGramRing;
  // "incremental" (the flip of the tiling across the wall) or "recompute"
  // (the tiling of a trial form beyond the wall).
  std::string FlipMethod = "incremental";
  // Whether to drop the stabilizers and the adjacency records of the tiling
  // of a domain once its neighbours have been computed. They are the bulk
  // of the memory of a domain and are not needed afterwards: the tiles
  // themselves are kept for the covering optimization and the outputs.
  bool SlimDomains = true;
  // Whether f_adj computes the k-covering optimum of each domain (needed
  // when the tilings are dropped).
  bool ComputeCovering = false;
  // The lattice data of a generic reference form of the T-space, used by
  // the incremental flip for the stabilizers and equivalences of point
  // sets; set with the seed form.
  std::unique_ptr<DataLattice<T, Tint, Tgroup>> data_ref;
};

template <typename T, typename Tint, typename Tgroup>
void SetReferenceForm(DataIsoKDelaunayDomains<T, Tint, Tgroup> &data,
                      MyMatrix<T> const &GramMat) {
  std::ostream &os = data.data.rddo.os;
  int n = GramMat.rows();
  MyMatrix<T> SHV(0, n);
  CVPSolver<T, Tint> solver(GramMat, os);
  std::string choice_initial = "direct";
  MyMatrix<Tint> ShvGraverBasis = GetGraverBasis<T, Tint>(GramMat);
  data.data_ref.reset(new DataLattice<T, Tint, Tgroup>{
      n, SHV, solver, ShvGraverBasis, choice_initial,
      RecordDualDescOperation<T, Tgroup>(data.data.rddo.AllArr, os)});
}

template <typename T, typename Tint, typename Tgroup>
IsoKDelaunayDomain<T, Tint, Tgroup>
GetInitialIsoKDelaunayDomain(DataIsoKDelaunayDomains<T, Tint, Tgroup> &data) {
  std::ostream &os = data.data.rddo.os;
  LinSpaceMatrix<T> const &LinSpa = data.data.LinSpa;
  if (!is_integrally_saturated_matrix_space(LinSpa.ListMat)) {
    std::cerr << "ISO_K_DELAUNAY: The space should be integral and fully "
                 "span it\n";
    throw TerminalException{1};
  }
  // The order-k tiling of a form has many more tiles than its Delaunay
  // tessellation, and an integer form with small entries almost always has
  // an accidentally co-spherical tile in dimension 4 and above. The
  // coefficient bound therefore grows geometrically over the attempts, a
  // failed attempt being detected at the first non-generic tile.
  int N = 16;
  [[maybe_unused]] size_t n_iter = 0;
  while (true) {
    MyMatrix<T> GramMat =
        GetRandomPositiveDefiniteNoNontrivialSymm<T, Tint, Tgroup>(LinSpa, N, os);
#ifdef DEBUG_ISO_K_DELAUNAY
    os << "ISO_K_DELAUNAY: GetInitialIsoKDelaunayDomain n_iter=" << n_iter
       << " N=" << N << " GramMat=\n";
    WriteMatrix(os, GramMat);
#endif
    std::optional<KDelaunayTesselation<Tint, Tgroup>> opt =
        ComputeGenericKDelaunayTesselation<T, Tint, Tgroup>(
            GramMat, data.k, data.ListGramRing, data.data.rddo.AllArr, os);
    if (opt) {
      SetReferenceForm(data, GramMat);
      return BuildIsoKDelaunayDomain<T, Tint, Tgroup>(*opt, GramMat, LinSpa,
                                                      data.ListGramRing, os);
    }
    n_iter += 1;
    N *= 3;
  }
}

/*
  The domain across the facet of the domain x whose interior point (in
  T-space coordinates) is TestPt and whose inequality is eIneq: the tiling
  is recomputed at TestPt - eps eIneq + eps^2 r with r a random integral
  vector, and accepted when the form is positive definite, generic in the
  T-space and TestPt lies in the closure of its domain (exactly one wall was
  crossed). The initial step is a fraction of the size of TestPt. When the
  domain found does not contain TestPt, its violated inequalities V tell
  where the segment from TestPt to the trial point left the adjacent
  domain: the crossing parameter t = V.TestPt / (V.TestPt - V.trial) is in
  (0, 1] and the next step is half the smallest crossing. No symmetry test
  is done on the trial form: the tiling of any generic form of the adjacent
  domain gives its inequalities, whatever the automorphisms of the form.
 */
template <typename T, typename Tint, typename Tgroup>
IsoKDelaunayDomain<T, Tint, Tgroup>
FlipIsoKDelaunayDomain(DataIsoKDelaunayDomains<T, Tint, Tgroup> &data,
                       MyVector<T> const &TestPt, MyVector<T> const &eIneq) {
  std::ostream &os = data.data.rddo.os;
  LinSpaceMatrix<T> const &LinSpa = data.data.LinSpa;
  int dimSpace = LinSpa.ListMat.size();
  T norm_pt(0), norm_ineq(0);
  for (int u = 0; u < dimSpace; u++) {
    norm_pt += T_abs(TestPt(u));
    norm_ineq += T_abs(eIneq(u));
  }
  T eps = norm_pt / (T(8) * norm_ineq);
  [[maybe_unused]] size_t n_iter = 0;
  while (true) {
#ifdef DEBUG_ISO_K_DELAUNAY
    os << "ISO_K_DELAUNAY: FlipIsoKDelaunayDomain n_iter=" << n_iter
       << " eps=" << eps << " TestPt=" << StringVectorGAP(TestPt)
       << " eIneq=" << StringVectorGAP(eIneq) << "\n";
#endif
    n_iter += 1;
    MyVector<T> rnd = FuncRandomDirection<T>(dimSpace, 2);
    MyVector<T> cand = TestPt - eps * eIneq + eps * eps * rnd;
    T scal = eIneq.dot(cand);
    if (scal >= 0) {
#ifdef DEBUG_ISO_K_DELAUNAY
      os << "ISO_K_DELAUNAY: FlipIsoKDelaunayDomain: wrong side, halving\n";
#endif
      eps /= 2;
      continue;
    }
    MyMatrix<T> GramMat = LINSPA_GetMatrixInTspace(LinSpa, cand);
    if (!IsPositiveDefinite(GramMat, os)) {
#ifdef DEBUG_ISO_K_DELAUNAY
      os << "ISO_K_DELAUNAY: FlipIsoKDelaunayDomain: not positive definite, halving\n";
#endif
      eps /= 2;
      continue;
    }
#ifdef DEBUG_ISO_K_DELAUNAY
    MicrosecondTime time_flip;
#endif
    std::optional<KDelaunayTesselation<Tint, Tgroup>> opt =
        ComputeGenericKDelaunayTesselation<T, Tint, Tgroup>(
            GramMat, data.k, data.ListGramRing, data.data.rddo.AllArr, os);
#ifdef DEBUG_ISO_K_DELAUNAY
    os << "|ISO_K_DELAUNAY: FlipIsoKDelaunayDomain: tiling|=" << time_flip
       << " n_tile=" << (opt ? opt->l_tiles.size() : 0) << "\n";
#endif
    if (!opt) {
#ifdef DEBUG_ISO_K_DELAUNAY
      os << "ISO_K_DELAUNAY: FlipIsoKDelaunayDomain: non-generic tiling, halving\n";
#endif
      eps /= 2;
      continue;
    }
    IsoKDelaunayDomain<T, Tint, Tgroup> dom = BuildIsoKDelaunayDomain<T, Tint, Tgroup>(
        *opt, GramMat, LinSpa, data.ListGramRing, os);
    // TestPt has to be in the closure of the new domain, otherwise the step
    // went too far and crossed more than one wall: the violated
    // inequalities bound the admissible step.
    MyMatrix<T> FAC_T = UniversalMatrixConversion<T, Tint>(dom.ListIneq);
    int n_ineq = FAC_T.rows();
    std::optional<T> t_min;
    for (int i_ineq = 0; i_ineq < n_ineq; i_ineq++) {
      T val_pt(0), val_cand(0);
      for (int u = 0; u < dimSpace; u++) {
        AddMul(val_pt, FAC_T(i_ineq, u), TestPt(u));
        AddMul(val_cand, FAC_T(i_ineq, u), cand(u));
      }
      if (val_pt < 0) {
        T t = val_pt / (val_pt - val_cand);
#ifdef DEBUG_ISO_K_DELAUNAY
        os << "ISO_K_DELAUNAY: violated V=" << StringVectorGAP(GetMatrixRow(FAC_T, i_ineq))
           << " val_pt=" << val_pt << " val_cand=" << val_cand << " t=" << t << "\n";
#endif
        if (!t_min || t < *t_min) {
          t_min = t;
        }
      }
    }
    if (t_min) {
#ifdef DEBUG_ISO_K_DELAUNAY
      os << "ISO_K_DELAUNAY: FlipIsoKDelaunayDomain: the facet point is not in "
            "the closure of the new domain, crossing at t=" << *t_min << "\n";
#endif
      eps = eps * (*t_min) / 2;
      continue;
    }
    return dom;
  }
}

/*
  Incremental flip.

  Crossing the wall V of a domain changes only the tiles carrying V as an
  inequality: the pairs of adjacent tiles whose wall inequality (a) is V and
  the tiles having a point of P_- whose inequality (b) is V. Those tiles are
  grouped into clusters, the connected components of the graph of the
  V-adjacencies, each cluster being one merged tile at the wall with merged
  sphere data (P_-^w, P_0^w): P_0^w is the union of the sets P_0 of the
  members together with the points of P_- leaving through V, and P_-^w is
  the common remainder of the sets P_-.

  The k-subsets P_-^w cup T, T a (k - |P_-^w|)-subset of P_0^w, lifted with
  the heights of a form of the current domain, have the current tiles of the
  cluster as the facets of their lower hull. The heights being linear in the
  form and coplanar on the wall, the facets of the UPPER hull of the same
  lifted points are the tiles of the cluster on the other side of the wall.
  The new tiles are thus read off one small convex hull per cluster, with
  no trial form and no search for adjacent tiles by sweeping spheres.

  The adjacencies of the new tiles are of two kinds. Two new tiles of a
  cluster sharing a facet have the same facet points (the sums of the
  lifted points on their common ridge), so the internal adjacencies are
  matched by facet point sets. A boundary facet of the cluster is a facet of
  an old member facing a tile outside the cluster, and it is a facet of
  exactly one new tile: that tile is found geometrically, as the new tile
  containing the centroid of the old facet with the facet inequality tight
  (the lifted points over a boundary facet need not be coplanar, so the
  vertex sets of the old and new facets can differ by non-extreme points,
  which rules out a matching by vertex sets). The records of the untouched
  tiles pointing into a cluster are re-resolved the same way.

  Everything is done up to the group G x Z^n with G the automorphisms of a
  fixed generic reference form of the T-space, through the stabilizer and
  equivalence tests of the Delaunay code applied to point sets: a tile is
  determined by P_0 and a cluster by P_0^w. The representatives of the new
  tiles are the first occurrences within a cluster class; their stabilizers
  are computed once. No metric enters except through the reference form,
  whose choice does not matter for generic forms.
 */
namespace iso_k_flip {

template <typename Tint> using Tkey = std::vector<MyVector<Tint>>;

template <typename Tint>
bool SameVector(MyVector<Tint> const &a, MyVector<Tint> const &b) {
  if (a.size() != b.size()) {
    return false;
  }
  for (int u = 0; u < a.size(); u++) {
    if (a(u) != b(u)) {
      return false;
    }
  }
  return true;
}

template <typename Tint> Tkey<Tint> GetKey(MyMatrix<Tint> const &M) {
  Tkey<Tint> key;
  int n_row = M.rows();
  for (int i_row = 0; i_row < n_row; i_row++) {
    key.push_back(GetMatrixRow(M, i_row));
  }
  std::sort(key.begin(), key.end());
  return key;
}

template <typename Tint>
MyMatrix<Tint> MatrixFromRows(std::vector<MyVector<Tint>> const &l_row,
                              int const &n_col) {
  int n_row = l_row.size();
  MyMatrix<Tint> M(n_row, n_col);
  for (int i_row = 0; i_row < n_row; i_row++) {
    for (int u = 0; u < n_col; u++) {
      M(i_row, u) = l_row[i_row](u);
    }
  }
  return M;
}

template <typename Tint>
KDelaunayTile<Tint> PlaceTile(KDelaunayTile<Tint> const &tile,
                              MyMatrix<Tint> const &M) {
  return {tile.EXT * M, tile.INT * M};
}

// The tiles are in sum coordinates: an affine map x -> x A + b of the
// lattice moves a vertex, the sum of k points, by k b. This is the matrix
// acting on the vertices for the matrix M acting on the points.
template <typename Tint>
MyMatrix<Tint> SumMatrix(MyMatrix<Tint> const &M, int const &k) {
  MyMatrix<Tint> Ms = M;
  int n1 = M.cols();
  for (int u = 1; u < n1; u++) {
    Ms(0, u) *= k;
  }
  return Ms;
}

// A facet of a representative with the tile across it (as rep * eBigMat).
template <typename Tint> struct FacetRecord {
  Face eInc;
  MyMatrix<Tint> eBigMat;
  int iOrb;
  MyVector<Tint> eIneq;
};

template <typename Tint> struct RepInfo {
  KDelaunayTileVertices<Tint> vert;
  // All the facets, obtained from the stored orbit representatives by the
  // action of the stabilizer.
  std::vector<FacetRecord<Tint>> l_facet;
};

// The permutation of the rows of EXT induced on the sums of the tile
// vertices, together with the affine matrix of the isometry.
template <typename Tint, typename Tgroup>
RepInfo<Tint> GetRepInfo(KDelaunay_Entry<Tint, Tgroup> const &ent,
                         int const &k) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  KDelaunayTile<Tint> const &tile = ent.tile;
  int n = tile.EXT.cols() - 1;
  KDelaunayTileVertices<Tint> vert = GetKDelaunayTileVertices(tile, k);
  int n_vert = vert.EXTsum.rows();
  std::unordered_map<MyVector<Tint>, size_t> map;
  for (int i_vert = 0; i_vert < n_vert; i_vert++) {
    MyVector<Tint> V(n);
    for (int u = 0; u < n; u++) {
      V(u) = vert.EXTsum(i_vert, u + 1);
    }
    map[V] = i_vert;
  }
  MyVector<Tint> base = ZeroVector<Tint>(n);
  for (int i_row = 0; i_row < tile.INT.rows(); i_row++) {
    for (int u = 0; u < n; u++) {
      base(u) += tile.INT(i_row, u + 1);
    }
  }
  std::vector<std::vector<Tidx>> l_perm;
  std::vector<MyMatrix<Tint>> l_mat;
  for (auto &eGen : ent.GRP.GeneratorsOfGroup()) {
    std::vector<Tidx> l_pos(n_vert);
    for (int i_vert = 0; i_vert < n_vert; i_vert++) {
      MyVector<Tint> V = base;
      for (auto &i_row : vert.ListSubset[i_vert]) {
        int i_img = OnPoints(i_row, eGen);
        for (int u = 0; u < n; u++) {
          V(u) += tile.EXT(i_img, u + 1);
        }
      }
      auto iter = map.find(V);
      if (iter == map.end()) {
        std::cerr << "ISO_K_FLIP: GetRepInfo: the image of a vertex is not a "
                     "vertex\n";
        throw TerminalException{1};
      }
      l_pos[i_vert] = iter->second;
    }
    l_perm.push_back(l_pos);
    l_mat.push_back(RepresentVertexPermutation(tile.EXT, tile.EXT, eGen));
  }
  std::vector<FacetRecord<Tint>> l_facet;
  std::set<Face> seen;
  MyMatrix<Tint> Id = IdentityMat<Tint>(n + 1);
  for (auto &rec : ent.ListAdj) {
    if (seen.count(rec.eInc) > 0) {
      continue;
    }
    std::vector<std::pair<Face, MyMatrix<Tint>>> queue{{rec.eInc, Id}};
    seen.insert(rec.eInc);
    size_t head = 0;
    while (head < queue.size()) {
      Face f = queue[head].first;
      MyMatrix<Tint> A = queue[head].second;
      head++;
      l_facet.push_back({f, rec.eBigMat * A, rec.iOrb, rec.eIneq});
      for (size_t i_gen = 0; i_gen < l_perm.size(); i_gen++) {
        Face f_img(n_vert);
        for (int i_vert = 0; i_vert < n_vert; i_vert++) {
          if (f[i_vert] == 1) {
            f_img[l_perm[i_gen][i_vert]] = 1;
          }
        }
        if (seen.count(f_img) == 0) {
          seen.insert(f_img);
          queue.push_back({f_img, A * l_mat[i_gen]});
        }
      }
    }
  }
  return {std::move(vert), std::move(l_facet)};
}

// The (homogeneous) points of the facet f of a tile placed by M.
template <typename Tint>
Tkey<Tint> GetFacetPoints(KDelaunayTileVertices<Tint> const &vert,
                          Face const &f, MyMatrix<Tint> const &M) {
  std::vector<MyVector<Tint>> l_row;
  int n_vert = vert.EXTsum.rows();
  for (int i_vert = 0; i_vert < n_vert; i_vert++) {
    if (f[i_vert] == 1) {
      MyVector<Tint> V = GetMatrixRow(vert.EXTsum, i_vert);
      l_row.push_back(M.transpose() * V);
    }
  }
  std::sort(l_row.begin(), l_row.end());
  return l_row;
}

template <typename Tint> struct Cluster {
  std::vector<int> l_orb;
  std::vector<MyMatrix<Tint>> l_M;
  std::vector<KDelaunayTile<Tint>> l_tile;
  std::vector<Tkey<Tint>> l_key;
  MyMatrix<Tint> P0w;
  MyMatrix<Tint> Pmw;
  Tkey<Tint> P0w_key;
};

template <typename Tint, typename Tgroup>
Cluster<Tint> BuildCluster(KDelaunayTesselation<Tint, Tgroup> const &DT,
                           std::vector<RepInfo<Tint>> const &l_info,
                           MyVector<Tint> const &V, int const &iOrb0,
                           MyMatrix<Tint> const &M0) {
  int n = DT.l_tiles[0].tile.EXT.cols() - 1;
  Cluster<Tint> C;
  std::set<Tkey<Tint>> seen;
  std::vector<MyVector<Tint>> l_ein;
  auto insert = [&](int iOrb, MyMatrix<Tint> const &M) -> bool {
    KDelaunayTile<Tint> tile = PlaceTile(DT.l_tiles[iOrb].tile, M);
    Tkey<Tint> key = GetKey(tile.EXT);
    if (!seen.insert(key).second) {
      return false;
    }
    C.l_orb.push_back(iOrb);
    C.l_M.push_back(M);
    C.l_tile.push_back(std::move(tile));
    C.l_key.push_back(std::move(key));
    return true;
  };
  insert(iOrb0, M0);
  size_t head = 0;
  while (head < C.l_orb.size()) {
    int iOrb = C.l_orb[head];
    MyMatrix<Tint> M = C.l_M[head];
    head++;
    for (auto &rec : l_info[iOrb].l_facet) {
      if (SameVector(rec.eIneq, V)) {
        insert(rec.iOrb, rec.eBigMat * M);
      }
    }
    KDelaunay_Entry<Tint, Tgroup> const &ent = DT.l_tiles[iOrb];
    for (int i_row = 0; i_row < ent.tile.INT.rows(); i_row++) {
      if (SameVector(ent.ListIntIneq[i_row], V)) {
        MyVector<Tint> p = GetMatrixRow(ent.tile.INT, i_row);
        l_ein.push_back(M.transpose() * p);
      }
    }
  }
  std::set<MyVector<Tint>> set_p0w;
  for (auto &tile : C.l_tile) {
    for (int i_row = 0; i_row < tile.EXT.rows(); i_row++) {
      set_p0w.insert(GetMatrixRow(tile.EXT, i_row));
    }
  }
  for (auto &p : l_ein) {
    set_p0w.insert(p);
  }
  C.P0w_key = Tkey<Tint>(set_p0w.begin(), set_p0w.end());
  C.P0w = MatrixFromRows(C.P0w_key, n + 1);
  std::vector<MyVector<Tint>> l_pmw;
  KDelaunayTile<Tint> const &tile0 = C.l_tile[0];
  for (int i_row = 0; i_row < tile0.INT.rows(); i_row++) {
    MyVector<Tint> p = GetMatrixRow(tile0.INT, i_row);
    if (set_p0w.count(p) == 0) {
      l_pmw.push_back(p);
    }
  }
  std::sort(l_pmw.begin(), l_pmw.end());
  C.Pmw = MatrixFromRows(l_pmw, n + 1);
#ifdef SANITY_CHECK_ISO_K_DELAUNAY
  // The remainder of P_- is the same for all the members.
  for (auto &tile : C.l_tile) {
    std::vector<MyVector<Tint>> l_test;
    for (int i_row = 0; i_row < tile.INT.rows(); i_row++) {
      MyVector<Tint> p = GetMatrixRow(tile.INT, i_row);
      if (set_p0w.count(p) == 0) {
        l_test.push_back(p);
      }
    }
    std::sort(l_test.begin(), l_test.end());
    if (l_test != l_pmw) {
      std::cerr << "ISO_K_FLIP: BuildCluster: the members of a cluster do not "
                   "share the same inner points\n";
      throw TerminalException{1};
    }
  }
#endif
  return C;
}

// A new tile of a cluster with its facets. A facet is either internal (shared
// with another new tile of the cluster) or on the boundary of the cluster,
// where it faces an untouched tile or a new tile of an adjacent cluster.
// The boundary facets are matched by their point sets: an untouched tile
// keeps its facet, and two adjacent clusters re-partition their common
// boundary identically, so the facet points coincide on both sides.
template <typename Tint> struct NewFacet {
  Face eInc;
  Tkey<Tint> pts;
  bool internal = false;
  int other_tile = -1; // internal: the other new tile of the cluster
  bool assigned = false;
  // The record: the representative in the new tiling and the matrix
  // placing it, filled by the assignment phase.
  int rec_orb = -1;
  MyMatrix<Tint> rec_M;
};

template <typename Tint> struct NewTile {
  KDelaunayTile<Tint> tile;
  KDelaunayTileVertices<Tint> vert;
  std::vector<NewFacet<Tint>> l_facet;
  // The representative among the new tiles of the class and the matrix
  // with rep * M = this tile.
  int i_rep = -1;
  MyMatrix<Tint> M;
};

template <typename Tint> struct ClusterClass {
  Cluster<Tint> C;
  std::vector<NewTile<Tint>> l_new;
  std::vector<int> l_rep; // indices in l_new of the representatives
};

template <typename Tint>
Tkey<Tint> TransformKey(Tkey<Tint> const &pts, MyMatrix<Tint> const &Ms) {
  Tkey<Tint> ret;
  for (auto &p : pts) {
    ret.push_back(Ms.transpose() * p);
  }
  std::sort(ret.begin(), ret.end());
  return ret;
}

/*
  The new tiles of a cluster from the upper hull of its lifted points, with
  their facets and the internal ones matched.
 */
template <typename T, typename Tint, typename Tgroup>
ClusterClass<Tint> ComputeClusterClass(Cluster<Tint> const &C,
                                       MyMatrix<T> const &GramOld,
                                       int const &k, std::ostream &os) {
  int n = GramOld.rows();
  int m = C.P0w.rows();
  int iw = C.Pmw.rows();
  int jw = k - iw;
  if (jw <= 0 || jw >= m) {
    std::cerr << "ISO_K_FLIP: ComputeClusterClass: jw=" << jw << " m=" << m
              << " is not a valid merged tile\n";
    throw TerminalException{1};
  }
  // The lifted points: the sums of the k-subsets and their heights.
  MyVector<Tint> base = ZeroVector<Tint>(n);
  T h_base(0);
  for (int i_row = 0; i_row < iw; i_row++) {
    MyVector<T> p(n);
    for (int u = 0; u < n; u++) {
      base(u) += C.Pmw(i_row, u + 1);
      p(u) = UniversalScalarConversion<T, Tint>(C.Pmw(i_row, u + 1));
    }
    h_base += EvaluationQuadForm<T, T>(GramOld, p);
  }
  std::vector<T> l_height(m);
  for (int i_row = 0; i_row < m; i_row++) {
    MyVector<T> p(n);
    for (int u = 0; u < n; u++) {
      p(u) = UniversalScalarConversion<T, Tint>(C.P0w(i_row, u + 1));
    }
    l_height[i_row] = EvaluationQuadForm<T, T>(GramOld, p);
  }
  std::vector<MyVector<T>> l_lift;
  std::vector<std::vector<std::vector<int>>> l_subsets;
  std::map<MyVector<T>, size_t> map_lift;
  std::vector<int> subset(jw);
  for (int u = 0; u < jw; u++) {
    subset[u] = u;
  }
  while (true) {
    MyVector<T> lift(n + 2);
    lift(0) = 1;
    MyVector<Tint> sigma = base;
    T h = h_base;
    for (int u = 0; u < jw; u++) {
      for (int w = 0; w < n; w++) {
        sigma(w) += C.P0w(subset[u], w + 1);
      }
      h += l_height[subset[u]];
    }
    for (int w = 0; w < n; w++) {
      lift(w + 1) = UniversalScalarConversion<T, Tint>(sigma(w));
    }
    lift(n + 1) = h;
    auto iter = map_lift.find(lift);
    if (iter == map_lift.end()) {
      map_lift[lift] = l_lift.size();
      l_lift.push_back(lift);
      l_subsets.push_back({subset});
    } else {
      l_subsets[iter->second].push_back(subset);
    }
    int pos = jw - 1;
    while (pos >= 0 && subset[pos] == m - jw + pos) {
      pos--;
    }
    if (pos < 0) {
      break;
    }
    subset[pos]++;
    for (int u = pos + 1; u < jw; u++) {
      subset[u] = subset[u - 1] + 1;
    }
  }
  MyMatrix<T> LIFT = MatrixFromVectorFamily(l_lift);
  if (RankMat(LIFT) != n + 2) {
    std::cerr << "ISO_K_FLIP: the lifted points of a cluster do not span, "
                 "the form is on the wall\n";
    throw TerminalException{1};
  }
  MyMatrix<T> FAC = DirectDualDescription_mat(LIFT, os);
  int n_fac = FAC.rows();
  int n_lift = l_lift.size();
  ClusterClass<Tint> cl;
  cl.C = C;
  for (int i_fac = 0; i_fac < n_fac; i_fac++) {
    if (FAC(i_fac, n + 1) >= 0) {
      continue; // lower or vertical facet
    }
    // The subsets of the points on the facet.
    std::vector<int> count(m, 0);
    int n_sub = 0;
    for (int i_lift = 0; i_lift < n_lift; i_lift++) {
      T val(0);
      for (int u = 0; u < n + 2; u++) {
        AddMul(val, FAC(i_fac, u), l_lift[i_lift](u));
      }
      if (val == 0) {
        for (auto &sub : l_subsets[i_lift]) {
          n_sub++;
          for (auto &i_row : sub) {
            count[i_row]++;
          }
        }
      }
    }
    std::vector<MyVector<Tint>> l_int, l_ext;
    for (int i_row = 0; i_row < iw; i_row++) {
      l_int.push_back(GetMatrixRow(C.Pmw, i_row));
    }
    for (int i_row = 0; i_row < m; i_row++) {
      if (count[i_row] == n_sub) {
        l_int.push_back(GetMatrixRow(C.P0w, i_row));
      } else if (count[i_row] > 0) {
        l_ext.push_back(GetMatrixRow(C.P0w, i_row));
      }
    }
    NewTile<Tint> nt;
    nt.tile = {MatrixFromRows(l_ext, n + 1), MatrixFromRows(l_int, n + 1)};
    if (!IsValidKDelaunayTile(nt.tile, k)) {
      std::cerr << "ISO_K_FLIP: a new tile is not valid |INT|="
                << nt.tile.INT.rows() << " |EXT|=" << nt.tile.EXT.rows()
                << "\n";
      throw TerminalException{1};
    }
    nt.vert = GetKDelaunayTileVertices(nt.tile, k);
    MyMatrix<T> EXTsum_T = UniversalMatrixConversion<T, Tint>(nt.vert.EXTsum);
    MyMatrix<T> FACtile_T = DirectDualDescription_mat(EXTsum_T, os);
    int n_facet = FACtile_T.rows();
    int n_vert = nt.vert.EXTsum.rows();
    for (int i_facet = 0; i_facet < n_facet; i_facet++) {
      NewFacet<Tint> nf;
      nf.eInc = Face(n_vert);
      for (int i_vert = 0; i_vert < n_vert; i_vert++) {
        T val(0);
        for (int u = 0; u <= n; u++) {
          AddMul(val, FACtile_T(i_facet, u), EXTsum_T(i_vert, u));
        }
        if (val == 0) {
          nf.eInc[i_vert] = 1;
        }
      }
      nf.pts = GetFacetPoints(nt.vert, nf.eInc, IdentityMat<Tint>(n + 1));
      nt.l_facet.push_back(std::move(nf));
    }
    cl.l_new.push_back(std::move(nt));
  }
  int n_new = cl.l_new.size();
  // Internal facets: matched by their point sets.
  std::map<Tkey<Tint>, std::pair<int, int>> map_facet;
  for (int i_new = 0; i_new < n_new; i_new++) {
    int n_facet = cl.l_new[i_new].l_facet.size();
    for (int i_facet = 0; i_facet < n_facet; i_facet++) {
      NewFacet<Tint> &nf = cl.l_new[i_new].l_facet[i_facet];
      auto iter = map_facet.find(nf.pts);
      if (iter == map_facet.end()) {
        map_facet[nf.pts] = {i_new, i_facet};
      } else {
        NewFacet<Tint> &nf2 =
            cl.l_new[iter->second.first].l_facet[iter->second.second];
        nf.internal = true;
        nf.other_tile = iter->second.first;
        nf2.internal = true;
        nf2.other_tile = i_new;
      }
    }
  }
  return cl;
}

} // namespace iso_k_flip

/*
  The domain across the wall V of the domain x, by the incremental flip of
  its tiling. TestPt (an interior point of the facet) certifies the result:
  it has to lie in the closure of the new domain, which has to have the
  reversed wall among its inequalities.
 */
template <typename T, typename Tint, typename Tgroup>
IsoKDelaunayDomain<T, Tint, Tgroup>
FlipIsoKDelaunayDomainIncremental(DataIsoKDelaunayDomains<T, Tint, Tgroup> &data,
                                  IsoKDelaunayDomain<T, Tint, Tgroup> const &x,
                                  std::vector<iso_k_flip::RepInfo<Tint>> const &l_info,
                                  MyVector<Tint> const &V,
                                  MyVector<T> const &TestPt) {
  using namespace iso_k_flip;
  std::ostream &os = data.data.rddo.os;
  LinSpaceMatrix<T> const &LinSpa = data.data.LinSpa;
  DataLattice<T, Tint, Tgroup> &data_ref = *data.data_ref;
  KDelaunayTesselation<Tint, Tgroup> const &DT = x.DT;
  int k = DT.k;
  int n = LinSpa.n;
  int n_orb = DT.l_tiles.size();
  MyMatrix<Tint> Id = IdentityMat<Tint>(n + 1);
  MyMatrix<T> GramOld = UniversalMatrixConversion<T, Tint>(x.GramMat);
#ifdef TIMINGS_ISO_K_DELAUNAY
  MicrosecondTime time;
#endif
  // The changed representatives.
  std::vector<int> l_changed(n_orb, 0);
  for (int iOrb = 0; iOrb < n_orb; iOrb++) {
    for (auto &rec : l_info[iOrb].l_facet) {
      if (SameVector(rec.eIneq, V)) {
        l_changed[iOrb] = 1;
      }
    }
    for (auto &eIneq : DT.l_tiles[iOrb].ListIntIneq) {
      if (SameVector(eIneq, V)) {
        l_changed[iOrb] = 1;
      }
    }
  }
  // Phase A: the cluster classes, from the changed representatives, with
  // the new tiles and their representatives. Every cluster of the tiling
  // contains a changed tile, hence is equivalent to one of those classes.
  std::vector<ClusterClass<Tint>> l_class;
  std::map<Tkey<Tint>, std::pair<int, MyMatrix<Tint>>> map_cluster;
  auto get_class = [&](int iOrb, MyMatrix<Tint> const &M)
      -> std::pair<int, MyMatrix<Tint>> {
    Cluster<Tint> C = BuildCluster(DT, l_info, V, iOrb, M);
    auto iter = map_cluster.find(C.P0w_key);
    if (iter != map_cluster.end()) {
      return iter->second;
    }
    MyMatrix<T> P0w_T = UniversalMatrixConversion<T, Tint>(C.P0w);
    int n_class = l_class.size();
    for (int i_class = 0; i_class < n_class; i_class++) {
      MyMatrix<T> P0w_ref_T =
          UniversalMatrixConversion<T, Tint>(l_class[i_class].C.P0w);
      std::optional<MyMatrix<T>> opt =
          Polytope_TestEquivalence<T, Tint, Tgroup>(data_ref, P0w_ref_T, P0w_T);
      if (opt) {
        std::pair<int, MyMatrix<Tint>> ret{
            i_class, UniversalMatrixConversion<Tint, T>(*opt)};
        map_cluster[C.P0w_key] = ret;
        return ret;
      }
    }
    ClusterClass<Tint> cl = ComputeClusterClass<T, Tint, Tgroup>(C, GramOld, k, os);
    int n_new = cl.l_new.size();
    for (int i_new = 0; i_new < n_new; i_new++) {
      NewTile<Tint> &nt = cl.l_new[i_new];
      MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(nt.tile.EXT);
      for (auto &i_rep : cl.l_rep) {
        MyMatrix<T> EXTrep_T =
            UniversalMatrixConversion<T, Tint>(cl.l_new[i_rep].tile.EXT);
        std::optional<MyMatrix<T>> opt =
            Polytope_TestEquivalence<T, Tint, Tgroup>(data_ref, EXTrep_T, EXT_T);
        if (opt) {
          nt.i_rep = i_rep;
          nt.M = UniversalMatrixConversion<Tint, T>(*opt);
          break;
        }
      }
      if (nt.i_rep < 0) {
        nt.i_rep = i_new;
        nt.M = Id;
        cl.l_rep.push_back(i_new);
      }
    }
    std::pair<int, MyMatrix<Tint>> ret{n_class, Id};
    map_cluster[C.P0w_key] = ret;
    l_class.push_back(std::move(cl));
    return ret;
  };
  for (int iOrb = 0; iOrb < n_orb; iOrb++) {
    if (l_changed[iOrb] == 1) {
      get_class(iOrb, Id);
    }
  }
  int n_class = l_class.size();
#ifdef TIMINGS_ISO_K_DELAUNAY
  os << "|ISO_K_FLIP: clusters|=" << time << "\n";
#endif
#ifdef DEBUG_ISO_K_DELAUNAY
  {
    int n_changed = 0;
    for (auto &c : l_changed) {
      n_changed += c;
    }
    os << "ISO_K_FLIP: n_orb=" << n_orb << " n_changed=" << n_changed
       << " n_class=" << n_class << "\n";
  }
#endif
  // The indices of the representatives in the new tiling: the untouched
  // ones first, then the representatives of the new tiles per class.
  std::vector<int> map_old(n_orb, -1);
  int n_untouched = 0;
  for (int iOrb = 0; iOrb < n_orb; iOrb++) {
    if (l_changed[iOrb] == 0) {
      map_old[iOrb] = n_untouched;
      n_untouched++;
    }
  }
  std::vector<std::vector<int>> map_new(n_class);
  int pos = n_untouched;
  for (int i_class = 0; i_class < n_class; i_class++) {
    map_new[i_class].resize(l_class[i_class].l_new.size(), -1);
    for (auto &i_rep : l_class[i_class].l_rep) {
      map_new[i_class][i_rep] = pos;
      pos++;
    }
  }
  int n_orb_new = pos;
  // The record of the new tile owning the facet with points pts (in the
  // coordinates of the tiling) on the side of the placed changed tile
  // (iOrb_t, M_t): the facet of a new tile of that cluster with the same
  // points.
  auto resolve = [&](int iOrb_t, MyMatrix<Tint> const &M_t,
                     Tkey<Tint> const &pts) -> std::pair<int, MyMatrix<Tint>> {
    std::pair<int, MyMatrix<Tint>> pair = get_class(iOrb_t, M_t);
    int i_class = pair.first;
    MyMatrix<Tint> const &M_c = pair.second;
    MyMatrix<Tint> Ms = SumMatrix(M_c, k);
    ClusterClass<Tint> const &cl = l_class[i_class];
    int n_new = cl.l_new.size();
    for (int i_new = 0; i_new < n_new; i_new++) {
      NewTile<Tint> const &nt = cl.l_new[i_new];
      for (auto &nf : nt.l_facet) {
        if (!nf.internal && TransformKey(nf.pts, Ms) == pts) {
          return {map_new[i_class][nt.i_rep], nt.M * M_c};
        }
      }
    }
    std::cerr << "ISO_K_FLIP: resolve: no facet of the new tiles of the "
                 "cluster has the requested points\n";
    throw TerminalException{1};
  };
  // Phase B: the records of the facets of the new tiles, in the coordinates
  // of their class. The boundary facets are reached from the facets of the
  // members facing the outside: for an untouched neighbour the facet is
  // unchanged, for a changed neighbour the new facets of its cluster are
  // matched.
  for (int i_class = 0; i_class < n_class; i_class++) {
    ClusterClass<Tint> &cl = l_class[i_class];
    int n_new = cl.l_new.size();
    std::map<Tkey<Tint>, std::pair<int, int>> map_boundary;
    for (int i_new = 0; i_new < n_new; i_new++) {
      NewTile<Tint> &nt = cl.l_new[i_new];
      int n_facet = nt.l_facet.size();
      for (int i_facet = 0; i_facet < n_facet; i_facet++) {
        NewFacet<Tint> &nf = nt.l_facet[i_facet];
        if (nf.internal) {
          NewTile<Tint> const &nt2 = cl.l_new[nf.other_tile];
          nf.rec_orb = map_new[i_class][nt2.i_rep];
          nf.rec_M = nt2.M;
          nf.assigned = true;
        } else {
          map_boundary[nf.pts] = {i_new, i_facet};
        }
      }
    }
    auto assign = [&](Tkey<Tint> const &pts, int rec_orb,
                      MyMatrix<Tint> const &rec_M) -> void {
      auto iter = map_boundary.find(pts);
      if (iter == map_boundary.end()) {
        std::cerr << "ISO_K_FLIP: a facet reached from the outside is not a "
                     "boundary facet of the new tiles\n";
        throw TerminalException{1};
      }
      NewFacet<Tint> &nf =
          cl.l_new[iter->second.first].l_facet[iter->second.second];
      if (nf.assigned) {
        if (nf.rec_orb != rec_orb || nf.rec_M != rec_M) {
          std::cerr << "ISO_K_FLIP: a boundary facet is assigned twice with "
                       "different neighbours\n";
          throw TerminalException{1};
        }
        return;
      }
      nf.assigned = true;
      nf.rec_orb = rec_orb;
      nf.rec_M = rec_M;
    };
    std::set<Tkey<Tint>> set_member(cl.C.l_key.begin(), cl.C.l_key.end());
    int n_mem = cl.C.l_orb.size();
    for (int i_mem = 0; i_mem < n_mem; i_mem++) {
      int iOrb = cl.C.l_orb[i_mem];
      MyMatrix<Tint> const &M = cl.C.l_M[i_mem];
      MyMatrix<Tint> Ms = SumMatrix(M, k);
      for (auto &rec : l_info[iOrb].l_facet) {
        MyMatrix<Tint> M_out = rec.eBigMat * M;
        if (l_changed[rec.iOrb] == 0) {
          Tkey<Tint> pts = GetFacetPoints(l_info[iOrb].vert, rec.eInc, Ms);
          assign(pts, map_old[rec.iOrb], M_out);
          continue;
        }
        KDelaunayTile<Tint> tile_out = PlaceTile(DT.l_tiles[rec.iOrb].tile, M_out);
        if (set_member.count(GetKey(tile_out.EXT)) > 0) {
          continue; // an internal facet of the old subdivision
        }
        // A facet facing an adjacent cluster: its new facets on our side.
        std::pair<int, MyMatrix<Tint>> pair = get_class(rec.iOrb, M_out);
        int i_class2 = pair.first;
        MyMatrix<Tint> const &M_c2 = pair.second;
        MyMatrix<Tint> Ms2 = SumMatrix(M_c2, k);
        ClusterClass<Tint> const &cl2 = l_class[i_class2];
        for (auto &nt2 : cl2.l_new) {
          for (auto &nf2 : nt2.l_facet) {
            if (nf2.internal) {
              continue;
            }
            Tkey<Tint> pts2 = TransformKey(nf2.pts, Ms2);
            if (map_boundary.count(pts2) > 0) {
              assign(pts2, map_new[i_class2][nt2.i_rep], nt2.M * M_c2);
            }
          }
        }
      }
    }
    for (auto &nt : cl.l_new) {
      for (auto &nf : nt.l_facet) {
        if (!nf.assigned) {
          std::cerr << "ISO_K_FLIP: a facet of a new tile is neither internal "
                       "nor reached from the outside\n";
          throw TerminalException{1};
        }
      }
    }
  }
#ifdef TIMINGS_ISO_K_DELAUNAY
  os << "|ISO_K_FLIP: facet records|=" << time << "\n";
#endif
  // Phase C: the new tiling.
  std::vector<KDelaunay_Entry<Tint, Tgroup>> l_tiles_new(n_orb_new);
  for (int iOrb = 0; iOrb < n_orb; iOrb++) {
    if (l_changed[iOrb] == 1) {
      continue;
    }
    KDelaunay_Entry<Tint, Tgroup> const &ent = DT.l_tiles[iOrb];
    std::vector<KDelaunay_AdjO<Tint>> ListAdj;
    for (auto &rec : ent.ListAdj) {
      if (l_changed[rec.iOrb] == 0) {
        // The neighbour is unchanged, and so is the wall inequality.
        ListAdj.push_back({rec.eInc, rec.eBigMat, map_old[rec.iOrb], rec.eIneq});
      } else {
        Tkey<Tint> pts = GetFacetPoints(l_info[iOrb].vert, rec.eInc, Id);
        std::pair<int, MyMatrix<Tint>> pair = resolve(rec.iOrb, rec.eBigMat, pts);
        ListAdj.push_back({rec.eInc, pair.second, pair.first, {}});
      }
    }
    l_tiles_new[map_old[iOrb]] = {ent.tile, ent.GRP, std::move(ListAdj),
                                  ent.ListIntIneq};
  }
  for (int i_class = 0; i_class < n_class; i_class++) {
    ClusterClass<Tint> const &cl = l_class[i_class];
    for (auto &i_rep : cl.l_rep) {
      NewTile<Tint> const &nt = cl.l_new[i_rep];
      std::vector<KDelaunay_AdjO<Tint>> ListAdj;
      for (auto &nf : nt.l_facet) {
        ListAdj.push_back({nf.eInc, nf.rec_M, nf.rec_orb, {}});
      }
      MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(nt.tile.EXT);
      Tgroup GRP = Polytope_StabilizerKernel<T, Tint, Tgroup>(data_ref, EXT_T);
      l_tiles_new[map_new[i_class][i_rep]] = {nt.tile, std::move(GRP),
                                              std::move(ListAdj), {}};
    }
  }
  KDelaunayTesselation<Tint, Tgroup> DT_new{k, std::move(l_tiles_new)};
#ifdef TIMINGS_ISO_K_DELAUNAY
  os << "|ISO_K_FLIP: assembly|=" << time << "\n";
#endif
#ifdef SANITY_CHECK_ISO_K_DELAUNAY
  check_k_delaunay_tessellation(DT_new, os);
#endif
  IsoKDelaunayDomain<T, Tint, Tgroup> dom = BuildIsoKDelaunayDomain<T, Tint, Tgroup>(
      std::move(DT_new), GramOld, LinSpa, data.ListGramRing, os, false);
#ifdef TIMINGS_ISO_K_DELAUNAY
  os << "|ISO_K_FLIP: domain|=" << time << "\n";
#endif
  // Certification: the facet point lies in the closure of the new domain
  // and the crossed wall is one of its inequalities, reversed.
  {
    MyMatrix<T> FAC_T = UniversalMatrixConversion<T, Tint>(dom.ListIneq);
    int n_ineq = FAC_T.rows();
    int dimSpace = FAC_T.cols();
    bool has_reverse = false;
    MyVector<Tint> Vrev = ScalarCanonicalizationVector(MyVector<Tint>(-V));
    for (int i_ineq = 0; i_ineq < n_ineq; i_ineq++) {
      T val(0);
      for (int u = 0; u < dimSpace; u++) {
        AddMul(val, FAC_T(i_ineq, u), TestPt(u));
      }
      if (val < 0) {
        std::cerr << "ISO_K_FLIP: the facet point is not in the closure of "
                     "the flipped domain, inequality "
                  << StringVectorGAP(GetMatrixRow(dom.ListIneq, i_ineq))
                  << " has value " << val << "\n";
        throw TerminalException{1};
      }
      if (SameVector(MyVector<Tint>(GetMatrixRow(dom.ListIneq, i_ineq)), Vrev)) {
        has_reverse = true;
      }
    }
    if (!has_reverse) {
      std::cerr << "ISO_K_FLIP: the reversed wall is not an inequality of the "
                   "flipped domain\n";
      throw TerminalException{1};
    }
  }
  return dom;
}

template <typename T, typename Tint, typename Tgroup>
struct IsoKDelaunayDomain_AdjI {
  MyVector<Tint> V;
  IsoKDelaunayDomain<T, Tint, Tgroup> dom;
  // The interior point of the crossed facet, which certifies a flip.
  MyVector<T> TestPt;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, IsoKDelaunayDomain_AdjI<T, Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("V", eRec.V);
  ar &make_nvp("dom", eRec.dom);
  ar &make_nvp("TestPt", eRec.TestPt);
}
} // namespace boost::serialization


/*
  Compact storage of a domain for the enumeration. The enumeration keeps
  every domain found until it is processed, and a domain with its tiling
  (tiles, stabilizer chains, adjacency records and inequalities, all over
  the ring) costs hundreds of kilobytes, so the frontier of found domains
  dominates the memory. The compact form stores the coordinates and the
  matrices as machine integers, the stabilizers by their generators, and
  keeps only the Gram matrix and the invariant vector family exact, which
  are what the equivalence tests and the hash read. A domain is expanded
  when it is processed and, once processed, stored without its tiling.
 */
template <typename Tint>
int64_t CompactScalar(Tint const &x) {
  int64_t v = UniversalScalarConversion<int64_t, Tint>(x);
  if (UniversalScalarConversion<Tint, int64_t>(v) != x) {
    std::cerr << "ISO_K_DELAUNAY: CompactScalar: the value " << x
              << " does not fit in 64 bits\n";
    throw TerminalException{1};
  }
  return v;
}

template <typename Tint>
std::vector<int64_t> CompactMatrix(MyMatrix<Tint> const &M) {
  std::vector<int64_t> v;
  v.reserve(M.rows() * M.cols());
  for (int i = 0; i < M.rows(); i++) {
    for (int j = 0; j < M.cols(); j++) {
      v.push_back(CompactScalar(M(i, j)));
    }
  }
  return v;
}

template <typename Tint>
MyMatrix<Tint> ExpandMatrix(std::vector<int64_t> const &v, int const &n_row,
                            int const &n_col) {
  MyMatrix<Tint> M(n_row, n_col);
  size_t pos = 0;
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n_col; j++) {
      M(i, j) = UniversalScalarConversion<Tint, int64_t>(v[pos]);
      pos++;
    }
  }
  return M;
}

template <typename Tint>
std::vector<int64_t> CompactVector(MyVector<Tint> const &V) {
  std::vector<int64_t> v(V.size());
  for (int i = 0; i < V.size(); i++) {
    v[i] = CompactScalar(V(i));
  }
  return v;
}

template <typename Tint>
MyVector<Tint> ExpandVector(std::vector<int64_t> const &v) {
  MyVector<Tint> V(v.size());
  for (size_t i = 0; i < v.size(); i++) {
    V(i) = UniversalScalarConversion<Tint, int64_t>(v[i]);
  }
  return V;
}

struct KTileRecordCompact {
  Face eInc;
  std::vector<int64_t> eBigMat;
  int iOrb;
  std::vector<int64_t> eIneq;
};

struct KTileCompact {
  int n_ext;
  int n_int;
  std::vector<int64_t> ext;
  std::vector<int64_t> intp;
  std::vector<std::vector<uint32_t>> gens;
  std::vector<KTileRecordCompact> l_rec;
  std::vector<std::vector<int64_t>> l_int_ineq;
};

template <typename Tint> struct IsoKDelaunayDomainCompact {
  int k = 0;
  int n = 0;
  std::vector<KTileCompact> l_tile;
  std::vector<std::vector<int64_t>> ListIneq;
  MyMatrix<Tint> GramMat;
  MyMatrix<Tint> SHV;
};

namespace boost::serialization {
template <class Archive>
inline void serialize(Archive &ar, KTileRecordCompact &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("eInc", eRec.eInc);
  ar &make_nvp("eBigMat", eRec.eBigMat);
  ar &make_nvp("iOrb", eRec.iOrb);
  ar &make_nvp("eIneq", eRec.eIneq);
}
template <class Archive>
inline void serialize(Archive &ar, KTileCompact &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("n_ext", eRec.n_ext);
  ar &make_nvp("n_int", eRec.n_int);
  ar &make_nvp("ext", eRec.ext);
  ar &make_nvp("intp", eRec.intp);
  ar &make_nvp("gens", eRec.gens);
  ar &make_nvp("l_rec", eRec.l_rec);
  ar &make_nvp("l_int_ineq", eRec.l_int_ineq);
}
template <class Archive, typename Tint>
inline void serialize(Archive &ar, IsoKDelaunayDomainCompact<Tint> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("k", eRec.k);
  ar &make_nvp("n", eRec.n);
  ar &make_nvp("l_tile", eRec.l_tile);
  ar &make_nvp("ListIneq", eRec.ListIneq);
  ar &make_nvp("GramMat", eRec.GramMat);
  ar &make_nvp("SHV", eRec.SHV);
}
} // namespace boost::serialization

template <typename T, typename Tint, typename Tgroup>
IsoKDelaunayDomainCompact<Tint>
CompactDomain(IsoKDelaunayDomain<T, Tint, Tgroup> const &dom) {
  IsoKDelaunayDomainCompact<Tint> c;
  c.k = dom.DT.k;
  c.n = dom.GramMat.rows();
  for (auto &ent : dom.DT.l_tiles) {
    KTileCompact t;
    t.n_ext = ent.tile.EXT.rows();
    t.n_int = ent.tile.INT.rows();
    t.ext = CompactMatrix(ent.tile.EXT);
    t.intp = CompactMatrix(ent.tile.INT);
    for (auto &eGen : ent.GRP.GeneratorsOfGroup()) {
      std::vector<uint32_t> perm(t.n_ext);
      for (int i = 0; i < t.n_ext; i++) {
        perm[i] = OnPoints(i, eGen);
      }
      t.gens.push_back(std::move(perm));
    }
    for (auto &rec : ent.ListAdj) {
      t.l_rec.push_back({rec.eInc, CompactMatrix(rec.eBigMat), rec.iOrb,
                         CompactVector(rec.eIneq)});
    }
    for (auto &V : ent.ListIntIneq) {
      t.l_int_ineq.push_back(CompactVector(V));
    }
    c.l_tile.push_back(std::move(t));
  }
  for (int i = 0; i < dom.ListIneq.rows(); i++) {
    c.ListIneq.push_back(CompactVector(MyVector<Tint>(GetMatrixRow(dom.ListIneq, i))));
  }
  c.GramMat = dom.GramMat;
  c.SHV = dom.SHV;
  return c;
}

template <typename T, typename Tint, typename Tgroup>
IsoKDelaunayDomain<T, Tint, Tgroup>
ExpandDomain(IsoKDelaunayDomainCompact<Tint> const &c) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  IsoKDelaunayDomain<T, Tint, Tgroup> dom;
  dom.DT.k = c.k;
  int n1 = c.n + 1;
  for (auto &t : c.l_tile) {
    KDelaunay_Entry<Tint, Tgroup> ent;
    ent.tile = {ExpandMatrix<Tint>(t.ext, t.n_ext, n1),
                ExpandMatrix<Tint>(t.intp, t.n_int, n1)};
    std::vector<Telt> l_gen;
    for (auto &perm : t.gens) {
      std::vector<Tidx> v(perm.begin(), perm.end());
      l_gen.push_back(Telt(v));
    }
    ent.GRP = Tgroup(l_gen, static_cast<Tidx>(t.n_ext));
    for (auto &rec : t.l_rec) {
      ent.ListAdj.push_back({rec.eInc, ExpandMatrix<Tint>(rec.eBigMat, n1, n1),
                             rec.iOrb, ExpandVector<Tint>(rec.eIneq)});
    }
    for (auto &v : t.l_int_ineq) {
      ent.ListIntIneq.push_back(ExpandVector<Tint>(v));
    }
    dom.DT.l_tiles.push_back(std::move(ent));
  }
  int n_ineq = c.ListIneq.size();
  int dimSpace = n_ineq > 0 ? c.ListIneq[0].size() : 0;
  dom.ListIneq = MyMatrix<Tint>(n_ineq, dimSpace);
  for (int i = 0; i < n_ineq; i++) {
    for (int u = 0; u < dimSpace; u++) {
      dom.ListIneq(i, u) = UniversalScalarConversion<Tint, int64_t>(c.ListIneq[i][u]);
    }
  }
  dom.GramMat = c.GramMat;
  dom.SHV = c.SHV;
  return dom;
}

template <typename T, typename Tint, typename Tgroup>
struct IsoKDelaunayDomain_Obj {
  IsoKDelaunayDomainCompact<Tint> dom;
  // The irredundant inequalities and the permutation group induced on them
  // by the stabilizer of the domain, filled by f_adj.
  MyMatrix<Tint> ListIneqRed;
  Tgroup GRPperm;
  // The summary of the tiling (number of tile orbits and their types by
  // (|P_-|, |P_0|)), filled by f_adj, available after the tiling is dropped.
  int n_tile = 0;
  std::map<std::pair<int, int>, size_t> map_tile_type;
  // The k-covering optimum of the domain, computed by f_adj when requested.
  std::optional<covering_maxdet::MaxdetResult<double>> cov_opt;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, IsoKDelaunayDomain_Obj<T, Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("dom", eRec.dom);
  ar &make_nvp("ListIneqRed", eRec.ListIneqRed);
  ar &make_nvp("GRPperm", eRec.GRPperm);
}
} // namespace boost::serialization

template <typename T, typename Tint, typename Tgroup>
void WriteEntryGAP(std::ostream &os_out,
                   IsoKDelaunayDomain<T, Tint, Tgroup> const &ent) {
  os_out << "rec(DT:=";
  WriteKDelaunayTesselationGAP(os_out, ent.DT);
  os_out << ", ListIneq:=" << StringMatrixGAP(ent.ListIneq)
         << ", GramMat:=" << StringMatrixGAP(ent.GramMat) << ")";
}

template <typename T, typename Tint, typename Tgroup>
void WriteEntryPYTHON(std::ostream &os_out,
                      IsoKDelaunayDomain<T, Tint, Tgroup> const &ent) {
  os_out << "{\"DT\":";
  WriteKDelaunayTesselationPYTHON(os_out, ent.DT);
  os_out << ", \"ListIneq\":" << StringMatrixPYTHON(ent.ListIneq)
         << ", \"GramMat\":" << StringMatrixPYTHON(ent.GramMat) << "}";
}

template <typename T, typename Tint, typename Tgroup>
void WriteEntryGAP(std::ostream &os_out,
                   IsoKDelaunayDomain_Obj<T, Tint, Tgroup> const &ent) {
  os_out << "rec(dom:=";
  WriteEntryGAP(os_out, ExpandDomain<T, Tint, Tgroup>(ent.dom));
  os_out << ", ListIneqRed:=" << StringMatrixGAP(ent.ListIneqRed)
         << ", GRPperm:=" << ent.GRPperm.GapString() << ")";
}

template <typename T, typename Tint, typename Tgroup>
void WriteEntryPYTHON(std::ostream &os_out,
                      IsoKDelaunayDomain_Obj<T, Tint, Tgroup> const &ent) {
  os_out << "{\"dom\":";
  WriteEntryPYTHON(os_out, ExpandDomain<T, Tint, Tgroup>(ent.dom));
  os_out << ", \"ListIneqRed\":" << StringMatrixPYTHON(ent.ListIneqRed)
         << ", \"GRPperm\":" << ent.GRPperm.PythonString() << "}";
}

template <typename T, typename Tint, typename Tgroup>
void WriteDetailedEntryGAP(std::ostream &os_out,
                           DataIsoKDelaunayDomains<T, Tint, Tgroup> const &data,
                           IsoKDelaunayDomain_Obj<T, Tint, Tgroup> const &ent,
                           std::ostream &os) {
  LinSpaceMatrix<T> const &LinSpa = data.data.LinSpa;
  int dimSpace = LinSpa.ListMat.size();
  int n = LinSpa.n;
  os_out << "rec(k:=" << data.k;
  os_out << ", GRPpermSize:=" << ent.GRPperm.size();
  os_out << ", n_ineq:=" << ent.dom.ListIneq.size();
  os_out << ", n_ineq_red:=" << ent.ListIneqRed.rows();
  os_out << ", det:=" << DeterminantMat(ent.dom.GramMat);
  os_out << ", n_shv:=" << ent.dom.SHV.rows();
  os_out << ", n_tile:=" << ent.n_tile;
  MyMatrix<T> FAC = UniversalMatrixConversion<T, Tint>(ent.ListIneqRed);
  MyMatrix<T> EXT = DirectDualDescription_mat(FAC, os);
  int n_row = EXT.rows();
  std::map<int, size_t> map_rank;
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyMatrix<T> RayMat = ZeroMatrix<T>(n, n);
    for (int u = 0; u < dimSpace; u++) {
      MatAddMul(RayMat, EXT(i_row, u), LinSpa.ListMat[u]);
    }
    map_rank[RankMat(RayMat)] += 1;
  }
  auto write_map_val = [&](std::map<int, size_t> const &map) -> void {
    bool IsFirst = true;
    os_out << "[";
    for (auto &kv : map) {
      if (!IsFirst) {
        os_out << ",";
      }
      IsFirst = false;
      os_out << "[" << kv.first << "," << kv.second << "]";
    }
    os_out << "]";
  };
  os_out << ", n_ray:=" << n_row << ", ListRank:=";
  write_map_val(map_rank);
  // The tiles by their sphere data (|P_-|, |P_0|).
  os_out << ", ListTileType:=[";
  bool IsFirst = true;
  for (auto &kv : ent.map_tile_type) {
    if (!IsFirst) {
      os_out << ",";
    }
    IsFirst = false;
    os_out << "rec(nbInt:=" << kv.first.first << ", nbExt:=" << kv.first.second
           << ", count:=" << kv.second << ")";
  }
  os_out << "])";
}

template <typename Tint> struct IsoKDelaunayDomain_AdjO {
  MyVector<Tint> V;
  MyMatrix<Tint> eBigMat;
};

template <typename Tint>
void WriteEntryGAP(std::ostream &os_out, IsoKDelaunayDomain_AdjO<Tint> const &ent) {
  os_out << "rec(V:=" << StringVectorGAP(ent.V)
         << ", eBigMat:=" << StringMatrixGAP(ent.eBigMat) << ")";
}

template <typename Tint>
void WriteEntryPYTHON(std::ostream &os_out,
                      IsoKDelaunayDomain_AdjO<Tint> const &ent) {
  os_out << "{\"V\":" << StringVectorPYTHON(ent.V)
         << ", \"eBigMat\":" << StringMatrixPYTHON(ent.eBigMat) << "}";
}

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, IsoKDelaunayDomain_AdjO<Tint> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("V", eRec.V);
  ar &make_nvp("eBigMat", eRec.eBigMat);
}
} // namespace boost::serialization

template <typename T, typename Tint, typename Tgroup>
struct ResultIsoKDelaunayAdj {
  MyMatrix<Tint> ListIneqRed;
  Tgroup GRPperm;
  std::vector<IsoKDelaunayDomain_AdjI<T, Tint, Tgroup>> l_adj;
};

/*
  The adjacent domains: the irredundant inequalities are computed, the
  stabilizer of the interior form permutes them, and one flip is done per
  orbit of facets whose interior contains positive definite forms.
 */
template <typename T, typename Tint, typename Tgroup>
ResultIsoKDelaunayAdj<T, Tint, Tgroup>
get_result_iso_k_delaunay_adj(IsoKDelaunayDomain<T, Tint, Tgroup> const &x,
                              DataIsoKDelaunayDomains<T, Tint, Tgroup> &data) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::ostream &os = data.data.rddo.os;
  LinSpaceMatrix<T> const &LinSpa = data.data.LinSpa;
#ifdef TIMINGS_ISO_K_DELAUNAY
  MicrosecondTime time;
#endif
  Result_ComputeStabilizer_SHV<Tint, Tgroup> result_stab =
      LINSPA_ComputeStabilizer_SHV<Tint, Tgroup>(data.data.LinSpaRing, x.GramMat,
                                                 x.SHV, {}, os);
  std::vector<MyMatrix<Tint>> ListGenTot =
      result_stab.get_list_matrix(x.SHV, x.GramMat, data.data.LinSpaRing, os);
#ifdef TIMINGS_ISO_K_DELAUNAY
  os << "|ISO_K_DELAUNAY: f_adj, stabilizer|=" << time << "\n";
#endif
  MyMatrix<T> FAC_T = UniversalMatrixConversion<T, Tint>(x.ListIneq);
  std::vector<int> ListIrred = get_non_redundant_indices(FAC_T, os);
  size_t nbIrred = ListIrred.size();
  MyMatrix<Tint> FACred = SelectRow(x.ListIneq, ListIrred);
  MyMatrix<T> FACred_T = SelectRow(FAC_T, ListIrred);
#ifdef TIMINGS_ISO_K_DELAUNAY
  os << "|ISO_K_DELAUNAY: f_adj, get_non_redundant_indices|=" << time << "\n";
#endif
#ifdef DEBUG_ISO_K_DELAUNAY
  os << "ISO_K_DELAUNAY: f_adj: |FAC|=" << x.ListIneq.rows()
     << " nbIrred=" << nbIrred << " |ListGenTot|=" << ListGenTot.size() << "\n";
#endif
  std::vector<MyVector<T>> l_ineq;
  std::unordered_map<MyVector<Tint>, size_t> map_ineq;
  for (size_t i = 0; i < nbIrred; i++) {
    MyVector<Tint> eV = GetMatrixRow(FACred, i);
    l_ineq.push_back(UniversalVectorConversion<T, Tint>(eV));
    map_ineq[eV] = i;
  }
  std::vector<Telt> ListPermGens;
  for (auto &eGenTot : ListGenTot) {
    MyMatrix<T> eGenTot_T = UniversalMatrixConversion<T, Tint>(eGenTot);
    MyMatrix<T> MatSpace = matrix_in_t_space(eGenTot_T, LinSpa);
    std::vector<Tidx> l_pos(nbIrred);
    for (size_t i = 0; i < nbIrred; i++) {
      MyVector<T> eVimg = MatSpace * l_ineq[i];
      MyVector<Tint> eVimg_red = RemoveFractionVectorPlusCoeffRing(eVimg).TheVect;
      auto iter = map_ineq.find(eVimg_red);
      if (iter == map_ineq.end()) {
        std::cerr << "ISO_K_DELAUNAY: the stabilizer does not permute the "
                     "irredundant inequalities\n";
        throw TerminalException{1};
      }
      l_pos[i] = iter->second;
    }
    ListPermGens.push_back(Telt(l_pos));
  }
  Tgroup GRPperm = Tgroup(ListPermGens, nbIrred);
  std::vector<size_t> l_idx = DecomposeOrbitPoint_FullRepr(GRPperm);
#ifdef DEBUG_ISO_K_DELAUNAY
  os << "ISO_K_DELAUNAY: f_adj: |GRPperm|=" << GRPperm.size()
     << " |l_idx|=" << l_idx.size() << "\n";
#endif
  // The extreme rays of the cone: the interior point of a facet used for the
  // flip is the sum of its rays, each normalized, which is as far from the
  // boundary of the positive definite cone as the facet allows. The
  // geometrically unique interior point of the facet, an LP by-product, is
  // often positive definite only by a hair, and a step across the wall from
  // it has to be tiny, which makes the trial tiling large and slow.
  MyMatrix<T> EXTrays = DirectDualDescription_mat(FACred_T, os);
  int n_ray = EXTrays.rows();
  int dimSpace = LinSpa.ListMat.size();
#ifdef DEBUG_ISO_K_DELAUNAY
  os << "ISO_K_DELAUNAY: f_adj: n_ray=" << n_ray << "\n";
#endif
  auto get_facet_point = [&](size_t i) -> MyVector<T> {
    MyVector<T> sum = ZeroVector<T>(dimSpace);
    for (int i_ray = 0; i_ray < n_ray; i_ray++) {
      T scal(0);
      T norm(0);
      for (int u = 0; u < dimSpace; u++) {
        AddMul(scal, FACred_T(i, u), EXTrays(i_ray, u));
        norm += T_abs(EXTrays(i_ray, u));
      }
      if (scal == 0) {
        for (int u = 0; u < dimSpace; u++) {
          sum(u) += EXTrays(i_ray, u) / norm;
        }
      }
    }
    return sum;
  };
  // The tile information used by the incremental flip, shared by the flips
  // across all the facets of the domain.
  std::vector<iso_k_flip::RepInfo<Tint>> l_info;
  if (data.FlipMethod == "incremental") {
    for (auto &ent : x.DT.l_tiles) {
      l_info.push_back(iso_k_flip::GetRepInfo(ent, x.DT.k));
    }
  }
#ifdef TIMINGS_ISO_K_DELAUNAY
  os << "|ISO_K_DELAUNAY: f_adj, rep info|=" << time << "\n";
#endif
  std::vector<IsoKDelaunayDomain_AdjI<T, Tint, Tgroup>> l_adj;
  for (auto &i : l_idx) {
    MyVector<T> TestPt = get_facet_point(i);
    MyMatrix<T> TestMat = LINSPA_GetMatrixInTspace(LinSpa, TestPt);
    if (!IsPositiveDefinite(TestMat, os)) {
      // The sum of the rays of the facet is not positive definite. When the
      // rays are all positive semidefinite this says the facet lies on the
      // boundary of the positive definite cone; the LP interior point of the
      // facet settles the remaining cases.
      TestPt = GetSpaceInteriorPointFacet(FACred_T, i, os);
      TestMat = LINSPA_GetMatrixInTspace(LinSpa, TestPt);
      if (!IsPositiveDefinite(TestMat, os)) {
#ifdef DEBUG_ISO_K_DELAUNAY
        os << "ISO_K_DELAUNAY: f_adj: facet i=" << i
           << " is on the boundary of the positive definite cone\n";
#endif
        continue;
      }
    }
    MyVector<Tint> eIneq = GetMatrixRow(FACred, i);
    IsoKDelaunayDomain<T, Tint, Tgroup> dom = [&]() {
      if (data.FlipMethod == "recompute") {
        return FlipIsoKDelaunayDomain<T, Tint, Tgroup>(data, TestPt, l_ineq[i]);
      }
      if (data.FlipMethod == "incremental") {
        return FlipIsoKDelaunayDomainIncremental<T, Tint, Tgroup>(
            data, x, l_info, eIneq, TestPt);
      }
      std::cerr << "ISO_K_DELAUNAY: FlipMethod=" << data.FlipMethod
                << " should be incremental or recompute\n";
      throw TerminalException{1};
    }();
#ifdef TIMINGS_ISO_K_DELAUNAY
    os << "|ISO_K_DELAUNAY: f_adj, flip|=" << time << "\n";
#endif
    l_adj.push_back({eIneq, std::move(dom), TestPt});
  }
  return {std::move(FACred), std::move(GRPperm), std::move(l_adj)};
}

/*
  The extreme rays of a domain, from the irredundant inequalities filled by
  f_adj, as forms of the T-space. The rays of full rank are the positive
  definite forms among them; the rank distribution of all the rays is
  reported alongside.
 */
template <typename T, typename Tint, typename Tgroup>
void WriteFullRankRaysGAP(std::ostream &os_out,
                          IsoKDelaunayDomain_Obj<T, Tint, Tgroup> const &ent,
                          LinSpaceMatrix<T> const &LinSpa, std::ostream &os) {
  int dimSpace = LinSpa.ListMat.size();
  int n = LinSpa.n;
  MyMatrix<T> FAC = UniversalMatrixConversion<T, Tint>(ent.ListIneqRed);
  MyMatrix<T> EXT = DirectDualDescription_mat(FAC, os);
  int n_row = EXT.rows();
  std::map<int, size_t> map_rank;
  std::vector<MyMatrix<Tint>> l_full;
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyMatrix<T> RayMat = ZeroMatrix<T>(n, n);
    for (int u = 0; u < dimSpace; u++) {
      MatAddMul(RayMat, EXT(i_row, u), LinSpa.ListMat[u]);
    }
    int rnk = RankMat(RayMat);
    map_rank[rnk] += 1;
    if (rnk == n) {
      l_full.push_back(RemoveFractionMatrixPlusCoeffRing(RayMat).TheMat);
    }
  }
  os_out << "rec(n_ray:=" << n_row << ", ListRank:=[";
  bool IsFirst = true;
  for (auto &kv : map_rank) {
    if (!IsFirst) {
      os_out << ",";
    }
    IsFirst = false;
    os_out << "[" << kv.first << "," << kv.second << "]";
  }
  os_out << "], n_full_rank:=" << l_full.size() << ", ListFullRankRay:=[";
  IsFirst = true;
  for (auto &M : l_full) {
    if (!IsFirst) {
      os_out << ",\n";
    }
    IsFirst = false;
    os_out << StringMatrixGAP(M);
  }
  os_out << "])";
}

/*
  The k-covering optimization over a domain: the same determinant
  maximization as for the sphere covering (CoveringMaxdet.h), with the
  circumradius constraint imposed on the simplex of P_0 of every tile,
  since the k-covering radius is the largest sphere radius of the tiles.
 */
template <typename T, typename Tint, typename Tgroup>
covering_maxdet::CoveringData<T>
BuildKCoveringData(IsoKDelaunayDomain<T, Tint, Tgroup> const &x,
                   LinSpaceMatrix<T> const &LinSpa, MyMatrix<T> const &FAC,
                   [[maybe_unused]] std::ostream &os) {
  int n = LinSpa.n;
  int dim = LinSpa.ListMat.size();
  if (FAC.cols() != dim) {
    std::cerr << "ISO_K_DELAUNAY: BuildKCoveringData: the inequalities have "
              << FAC.cols() << " columns while the T-space has dimension "
              << dim << "\n";
    throw TerminalException{1};
  }
  std::vector<MyMatrix<T>> ListSimplex;
  for (auto &eEnt : x.DT.l_tiles) {
    ListSimplex.push_back(
        covering_maxdet::GetSimplexFromDelaunay<T, Tint>(eEnt.tile.EXT));
  }
  T point_density(1);
  return {n, dim, LinSpa.ListMat, FAC, std::move(ListSimplex), point_density};
}

template <typename T, typename Tint, typename Tgroup>
covering_maxdet::MaxdetResult<double>
OptimizeKCovering(IsoKDelaunayDomain<T, Tint, Tgroup> const &x,
                  LinSpaceMatrix<T> const &LinSpa, std::ostream &os) {
  using Tfloat = double;
  MyMatrix<T> FAC_T = UniversalMatrixConversion<T, Tint>(x.ListIneq);
  std::vector<int> ListIrred = get_non_redundant_indices(FAC_T, os);
  MyMatrix<T> FACred = SelectRow(FAC_T, ListIrred);
  // (a slimmed domain already carries only its irredundant inequalities,
  // in which case the elimination above is a no-op)
  covering_maxdet::CoveringData<T> cd =
      BuildKCoveringData<T, Tint, Tgroup>(x, LinSpa, FACred, os);
  covering_maxdet::CoveringData<Tfloat> cd_f =
      covering_maxdet::ConvertCoveringData<Tfloat, T>(cd);
  std::optional<MyVector<Tfloat>> opt_start =
      covering_maxdet::GetStartingPoint<T, Tfloat>(cd, cd_f, os);
  covering_maxdet::MaxdetResult<Tfloat> res;
  if (opt_start) {
    covering_maxdet::MaxdetOptions<Tfloat> opts =
        covering_maxdet::GetDefaultMaxdetOptions<Tfloat>();
    res = covering_maxdet::SolveCoveringMaxdet(cd_f, *opt_start, opts, os);
  } else {
    res.success = false;
    res.has_point = false;
    res.message = "the domain is too thin to be optimized in floating point";
  }
  return res;
}

// The optimization result as a GAP record (without the return), with the
// k-normalized density.
template <typename Tfloat>
void WriteKCoveringOptimumRecordGAP(std::ostream &os_out, int const &k,
                                    covering_maxdet::MaxdetResult<Tfloat> const &res) {
  std::streamsize prec = os_out.precision();
  os_out << std::setprecision(17) << std::showpoint;
  os_out << "rec(k:=" << k << ", success:=" << (res.success ? "true" : "false");
  os_out << ", has_point:=" << (res.has_point ? "true" : "false");
  os_out << ", message:=\"" << res.message << "\"";
  os_out << ", k_covering_density:=" << res.cov_density;
  os_out << ", k_covering_density_normalized:="
         << res.cov_density / static_cast<Tfloat>(k);
  os_out << ", k_covering_radius_sq:=" << res.cov_radius_sq;
  os_out << ", det:=" << res.det;
  os_out << ", gap_bound:=" << res.gap_bound;
  os_out << ", n_newton:=" << res.n_newton;
  os_out << ", GramMat:=[";
  for (int i = 0; i < res.Q.rows(); i++) {
    if (i > 0) {
      os_out << ", ";
    }
    os_out << "[";
    for (int j = 0; j < res.Q.cols(); j++) {
      if (j > 0) {
        os_out << ", ";
      }
      os_out << res.Q(i, j);
    }
    os_out << "]";
  }
  os_out << "])";
  os_out << std::setprecision(prec) << std::noshowpoint;
}

template <typename T, typename Tint, typename Tgroup>
struct DataIsoKDelaunayDomainsFunc {
  DataIsoKDelaunayDomains<T, Tint, Tgroup> data;
  using Tobj = IsoKDelaunayDomain_Obj<T, Tint, Tgroup>;
  using TadjI = IsoKDelaunayDomain_AdjI<T, Tint, Tgroup>;
  using TadjO = IsoKDelaunayDomain_AdjO<Tint>;
  std::ostream &get_os() { return data.data.rddo.os; }
  static Tobj make_obj(IsoKDelaunayDomain<T, Tint, Tgroup> const &dom) {
    Tobj x;
    x.dom = CompactDomain(dom);
    return x;
  }
  Tobj f_init() {
    IsoKDelaunayDomain<T, Tint, Tgroup> dom = GetInitialIsoKDelaunayDomain(data);
    return make_obj(dom);
  }
  size_t f_hash(size_t const &seed, Tobj const &x) {
    return LINSPA_Invariant_SHV<Tint>(seed, data.data.LinSpaRing, x.dom.GramMat,
                                      x.dom.SHV, {}, data.data.rddo.os);
  }
  std::optional<TadjO> f_repr(Tobj const &x, TadjI const &y) {
    std::optional<MyMatrix<Tint>> opt =
        LINSPA_TestEquivalenceGramMatrix_SHV<Tint, Tgroup>(
            data.data.LinSpaRing, x.dom.GramMat, y.dom.GramMat, x.dom.SHV,
            y.dom.SHV, {}, data.data.rddo.os);
    if (!opt) {
      return {};
    }
    return TadjO{y.V, *opt};
  }
  std::pair<Tobj, TadjO> f_spann(TadjI const &x) {
    n_found++;
    Tobj x_ret = make_obj(x.dom);
    MyMatrix<Tint> eBigMat = IdentityMat<Tint>(data.data.LinSpa.n);
    TadjO ret{x.V, eBigMat};
    return {std::move(x_ret), std::move(ret)};
  }
  size_t n_processed = 0;
  size_t n_found = 0;
  std::optional<std::vector<TadjI>> f_adj(Tobj &x_in) {
    IsoKDelaunayDomain<T, Tint, Tgroup> dom = ExpandDomain<T, Tint, Tgroup>(x_in.dom);
    ResultIsoKDelaunayAdj<T, Tint, Tgroup> result =
        get_result_iso_k_delaunay_adj(dom, data);
    x_in.ListIneqRed = result.ListIneqRed;
    x_in.GRPperm = result.GRPperm;
    x_in.n_tile = dom.DT.l_tiles.size();
    for (auto &eEnt : dom.DT.l_tiles) {
      x_in.map_tile_type[{eEnt.tile.INT.rows(), eEnt.tile.EXT.rows()}] += 1;
    }
    if (data.ComputeCovering) {
#ifdef TIMINGS_ISO_K_DELAUNAY
      MicrosecondTime time_cov;
#endif
      x_in.cov_opt = OptimizeKCovering<T, Tint, Tgroup>(dom, data.data.LinSpa,
                                                        data.data.rddo.os);
#ifdef TIMINGS_ISO_K_DELAUNAY
      data.data.rddo.os << "|ISO_K_DELAUNAY: f_adj, covering n_newton="
                        << x_in.cov_opt->n_newton << " success="
                        << x_in.cov_opt->success << "|=" << time_cov << "\n";
#endif
    }
    if (data.SlimDomains) {
      // The irredundant inequalities describe the same cone and are what the
      // extreme rays use; the tiling has served (neighbours, summary,
      // covering optimum) and is dropped.
      dom.ListIneq = result.ListIneqRed;
      dom.DT.l_tiles.clear();
    }
    x_in.dom = CompactDomain(dom);
    n_processed++;
    if (n_processed % 1000 == 0) {
      data.data.rddo.os << "ISO_K_DELAUNAY: progress n_processed=" << n_processed
                        << " n_found=" << n_found << "\n";
    }
    return result.l_adj;
  }
  Tobj f_adji_obj(TadjI const &x) { return make_obj(x.dom); }
  size_t f_complexity([[maybe_unused]] Tobj const &x) { return 0; }
};

/*
  The enumeration loop of the (L,k)-types.

  The generic adjacency scheme keeps every domain found until it is
  processed, with its tiling, and processes them depth first, so the memory
  is that of the frontier of found domains. Here a found domain is kept as
  a stub: its canonical interior form and invariant vector family (what the
  hash and the equivalence test read), and a recipe, the index of the domain
  it was found from and the wall crossed. When the stub is processed its
  tiling is rebuilt by one incremental flip of the tiling of its parent,
  which costs one more flip per domain on top of the flips finding its
  neighbours. The domains are processed first in first out, so that the
  parents whose children are pending are those of the last layer; a
  parent's tiling is freed when its last child has been processed.
 */
template <typename T, typename Tint, typename Tgroup> struct IsoKDomainStub {
  MyMatrix<Tint> GramMat;
  MyMatrix<Tint> SHV;
  int parent = -1;
  MyVector<Tint> V;
  MyVector<T> TestPt;
  // The compact tiling, present while some child is pending (and for the
  // root until it is processed).
  std::optional<IsoKDelaunayDomainCompact<Tint>> tiling;
  int n_pending = 0;
  bool processed = false;
  // Filled when processed.
  MyMatrix<Tint> ListIneqRed;
  Tgroup GRPperm;
  int n_tile = 0;
  std::map<std::pair<int, int>, size_t> map_tile_type;
  std::optional<covering_maxdet::MaxdetResult<double>> cov_opt;
};

template <typename T, typename Tint, typename Tgroup>
std::vector<DatabaseEntry_Serial<IsoKDelaunayDomain_Obj<T, Tint, Tgroup>,
                                 IsoKDelaunayDomain_AdjO<Tint>>>
EnumerateIsoKDelaunayDomains(DataIsoKDelaunayDomains<T, Tint, Tgroup> &data,
                             std::ostream &os) {
  using Tobj = IsoKDelaunayDomain_Obj<T, Tint, Tgroup>;
  using TadjO = IsoKDelaunayDomain_AdjO<Tint>;
  using Tstub = IsoKDomainStub<T, Tint, Tgroup>;
  using Tdom = IsoKDelaunayDomain<T, Tint, Tgroup>;
  LinSpaceMatrix<T> const &LinSpa = data.data.LinSpa;
  int n = LinSpa.n;
  MyMatrix<Tint> Id = IdentityMat<Tint>(n);
  std::vector<Tstub> l_obj;
  std::vector<std::vector<AdjO_Serial<TadjO>>> l_adj;
  std::unordered_map<size_t, std::vector<size_t>> indices_by_hash;
  size_t seed = 1234;
  auto get_hash = [&](MyMatrix<Tint> const &GramMat,
                      MyMatrix<Tint> const &SHV) -> size_t {
    return LINSPA_Invariant_SHV<Tint>(seed, data.data.LinSpaRing, GramMat, SHV,
                                      {}, os);
  };
  // The root.
  {
    Tdom dom0 = GetInitialIsoKDelaunayDomain(data);
    Tstub stub;
    stub.GramMat = dom0.GramMat;
    stub.SHV = dom0.SHV;
    stub.tiling = CompactDomain(dom0);
    indices_by_hash[get_hash(stub.GramMat, stub.SHV)].push_back(0);
    l_obj.push_back(std::move(stub));
    l_adj.push_back({});
  }
  size_t head = 0;
  while (head < l_obj.size()) {
    size_t i = head;
    head++;
    // The tiling of the domain: stored for the root, rebuilt from the
    // parent otherwise.
    Tdom dom;
    if (l_obj[i].tiling) {
      dom = ExpandDomain<T, Tint, Tgroup>(*l_obj[i].tiling);
      l_obj[i].tiling.reset();
    } else {
      int p = l_obj[i].parent;
      if (!l_obj[p].tiling) {
        std::cerr << "ISO_K_DELAUNAY: the tiling of the parent is missing\n";
        throw TerminalException{1};
      }
      Tdom dom_p = ExpandDomain<T, Tint, Tgroup>(*l_obj[p].tiling);
      std::vector<iso_k_flip::RepInfo<Tint>> l_info;
      for (auto &ent : dom_p.DT.l_tiles) {
        l_info.push_back(iso_k_flip::GetRepInfo(ent, dom_p.DT.k));
      }
      dom = FlipIsoKDelaunayDomainIncremental<T, Tint, Tgroup>(
          data, dom_p, l_info, l_obj[i].V, l_obj[i].TestPt);
#ifdef SANITY_CHECK_ISO_K_DELAUNAY
      if (dom.GramMat != l_obj[i].GramMat) {
        std::cerr << "ISO_K_DELAUNAY: the rebuilt domain does not have the "
                     "stored interior form\n";
        throw TerminalException{1};
      }
#endif
      l_obj[p].n_pending--;
      if (l_obj[p].n_pending == 0) {
        l_obj[p].tiling.reset();
      }
    }
    ResultIsoKDelaunayAdj<T, Tint, Tgroup> result =
        get_result_iso_k_delaunay_adj(dom, data);
    Tstub &x = l_obj[i];
    x.processed = true;
    x.ListIneqRed = result.ListIneqRed;
    x.GRPperm = result.GRPperm;
    x.n_tile = dom.DT.l_tiles.size();
    for (auto &eEnt : dom.DT.l_tiles) {
      x.map_tile_type[{eEnt.tile.INT.rows(), eEnt.tile.EXT.rows()}] += 1;
    }
    if (data.ComputeCovering) {
      x.cov_opt = OptimizeKCovering<T, Tint, Tgroup>(dom, LinSpa, os);
    }
    // The neighbours: known ones get an adjacency record, new ones a stub.
    std::vector<AdjO_Serial<TadjO>> ListAdj;
    for (auto &adj : result.l_adj) {
      size_t hash = get_hash(adj.dom.GramMat, adj.dom.SHV);
      std::vector<size_t> &vect = indices_by_hash[hash];
      std::optional<size_t> found;
      MyMatrix<Tint> eBigMat;
      for (auto &idx : vect) {
        std::optional<MyMatrix<Tint>> opt =
            LINSPA_TestEquivalenceGramMatrix_SHV<Tint, Tgroup>(
                data.data.LinSpaRing, l_obj[idx].GramMat, adj.dom.GramMat,
                l_obj[idx].SHV, adj.dom.SHV, {}, os);
        if (opt) {
          found = idx;
          eBigMat = *opt;
          break;
        }
      }
      if (found) {
        ListAdj.push_back({{adj.V, eBigMat}, static_cast<int>(*found)});
      } else {
        size_t idx_new = l_obj.size();
        Tstub stub;
        stub.GramMat = adj.dom.GramMat;
        stub.SHV = adj.dom.SHV;
        stub.parent = i;
        stub.V = adj.V;
        stub.TestPt = adj.TestPt;
        vect.push_back(idx_new);
        l_obj.push_back(std::move(stub));
        l_adj.push_back({});
        l_obj[i].n_pending++;
        ListAdj.push_back({{adj.V, Id}, static_cast<int>(idx_new)});
      }
    }
    l_adj[i] = std::move(ListAdj);
    if (l_obj[i].n_pending > 0) {
      // Records with all their facets, as the flips of the children need.
      l_obj[i].tiling = CompactDomain(dom);
    }
    if (head % 1000 == 0) {
      size_t n_active = 0;
      for (auto &stub : l_obj) {
        if (stub.tiling) {
          n_active++;
        }
      }
      os << "ISO_K_DELAUNAY: progress n_processed=" << head
         << " n_found=" << l_obj.size() << " n_active_tilings=" << n_active
         << "\n";
    }
  }
  // The database entries, in the form the writers read.
  std::vector<DatabaseEntry_Serial<Tobj, TadjO>> l_ret;
  size_t n_obj = l_obj.size();
  for (size_t i = 0; i < n_obj; i++) {
    Tstub &stub = l_obj[i];
    Tobj x;
    x.dom.k = data.k;
    x.dom.n = n;
    for (int i_row = 0; i_row < stub.ListIneqRed.rows(); i_row++) {
      x.dom.ListIneq.push_back(
          CompactVector(MyVector<Tint>(GetMatrixRow(stub.ListIneqRed, i_row))));
    }
    x.dom.GramMat = stub.GramMat;
    x.dom.SHV = stub.SHV;
    x.ListIneqRed = stub.ListIneqRed;
    x.GRPperm = stub.GRPperm;
    x.n_tile = stub.n_tile;
    x.map_tile_type = stub.map_tile_type;
    x.cov_opt = stub.cov_opt;
    l_ret.push_back({std::move(x), std::move(l_adj[i])});
  }
  return l_ret;
}

FullNamelist NAMELIST_GetStandard_COMPUTE_LATTICE_IsoKDelaunayDomains() {
  std::map<std::string, SingleBlock> ListBlock;
  // SYSTEM
  ListBlock["SYSTEM"] = SINGLEBLOCK_Get_System();
  // DATA
  {
    std::map<std::string, std::string> ListStringValues;
    ListStringValues["arithmetic"] = "gmp";
    ListStringValues["FileDualDescription"] = "unset";
    // Read by get_data_isodelaunay_domains; a common Gram matrix is not
    // supported by the order-k enumeration and has to stay unset.
    ListStringValues["CommonGramMat"] = "unset";
    // When set to a prefix P, each enumerated domain is written as a boost
    // text-archive of an IsoKDelaunayDomain to the file P<i>.
    ListStringValues["PrefixIsoKDelaunayDomains"] = "unset";
    // When set, the k-covering density is minimized over each enumerated
    // domain by the determinant maximization of CoveringMaxdet.h and the
    // results are written to that file as a GAP list.
    ListStringValues["FileCoveringOptimum"] = "null";
    // When set, the extreme rays of each enumerated domain are computed and
    // the ones of full rank (positive definite forms) are written to that
    // file as a GAP list, together with the rank distribution of all rays.
    ListStringValues["FileFullRankRays"] = "null";
    // The flip across a wall: "incremental" (the tiling is flipped, the
    // default) or "recompute" (the tiling of a trial form beyond the wall
    // is computed from scratch, slower, kept as a cross-check).
    ListStringValues["FlipMethod"] = "incremental";
    std::map<std::string, bool> ListBoolValues;
    // Drop the stabilizers and adjacencies of the tiling of a domain once
    // its neighbours are computed (they dominate the memory and are not
    // needed afterwards). With T the ObjectGAP output has trivial tile
    // stabilizers and no tile adjacencies.
    ListBoolValues["SlimDomains"] = true;
    std::map<std::string, int> ListIntValues;
    ListIntValues["k"] = 2;
    SingleBlock BlockDATA;
    BlockDATA.setListStringValues(ListStringValues);
    BlockDATA.setListIntValues(ListIntValues);
    BlockDATA.setListBoolValues(ListBoolValues);
    ListBlock["DATA"] = BlockDATA;
  }
  // TSPACE
  ListBlock["TSPACE"] = SINGLEBLOCK_Get_Tspace_Description();
  return FullNamelist(ListBlock);
}

template <typename T, typename Tint, typename Tgroup>
DataIsoKDelaunayDomains<T, Tint, Tgroup>
get_data_iso_k_delaunay_domains(FullNamelist const &eFull,
                                PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
                                std::ostream &os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  int k = BlockDATA.get_int("k");
  if (k < 1) {
    std::cerr << "ISO_K_DELAUNAY: k=" << k << " should be at least 1\n";
    throw TerminalException{1};
  }
  DataIsoDelaunayDomains<T, Tint, Tgroup> data =
      get_data_isodelaunay_domains<T, Tint, Tgroup>(eFull, AllArr, os);
  if (data.CommonGramMat) {
    std::cerr << "ISO_K_DELAUNAY: CommonGramMat is not supported by the "
                 "order-k enumeration\n";
    throw TerminalException{1};
  }
  std::vector<std::vector<Tint>> ListGramRing =
      GetListGramRing(data.LinSpa.ListLineMat);
  DataIsoKDelaunayDomains<T, Tint, Tgroup> data_k{
      std::move(data), k, std::move(ListGramRing),
      BlockDATA.get_string("FlipMethod"), BlockDATA.get_bool("SlimDomains"),
      false, nullptr};
  return data_k;
}

// clang-format off
#endif  // SRC_K_COVERINGS_ISOKDELAUNAYDOMAINS_H_
// clang-format on

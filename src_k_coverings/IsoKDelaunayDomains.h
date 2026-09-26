// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_K_COVERINGS_ISOKDELAUNAYDOMAINS_H_
#define SRC_K_COVERINGS_ISOKDELAUNAYDOMAINS_H_

// clang-format off
#include "LatticeKDelaunay.h"
#include "IsoDelaunayDomains.h"
#include "CoveringMaxdet.h"
#include <iomanip>
#include <map>
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
MyMatrix<Tint> ComputeIsoKDelaunayInequalities(
    KDelaunayTesselation<Tint, Tgroup> const &DT,
    std::vector<std::vector<Tint>> const &ListGramRing, std::ostream &os) {
  int k = DT.k;
  int dimSpace = ListGramRing.size();
  std::unordered_set<MyVector<Tint>> set_ineq;
  std::vector<MyVector<Tint>> l_ineq;
  auto insert = [&](MyVector<Tint> const &V) -> void {
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
  };
  for (auto &eEnt : DT.l_tiles) {
    KDelaunayTile<Tint> const &tile = eEnt.tile;
    int n = tile.EXT.cols() - 1;
    VoronoiInequalityPreComput<Tint> vipc =
        BuildVoronoiIneqPreComputeChecked<Tint>(tile.EXT, ListGramRing, os);
    auto get_phi = [&](MyVector<Tint> const &p_hom) -> MyVector<Tint> {
      return VoronoiLinearInequality(vipc, p_hom, ListGramRing, os);
    };
    int i = tile.INT.rows();
    MyVector<Tint> sumInt = ZeroVector<Tint>(dimSpace);
    for (int i_row = 0; i_row < i; i_row++) {
      MyVector<Tint> p_hom = GetMatrixRow(tile.INT, i_row);
      MyVector<Tint> phi = get_phi(p_hom);
      sumInt += phi;
      // (b): the point stays inside, -phi(p) >= 0.
      MyVector<Tint> V = -phi;
      insert(V);
    }
    // (a): one inequality per adjacency, from a k-subset of the adjacent
    // tile whose sum is off the shared facet.
    KDelaunayTileVertices<Tint> vert = GetKDelaunayTileVertices(tile, k);
    MyMatrix<T> EXTsum_T = UniversalMatrixConversion<T, Tint>(vert.EXTsum);
    SubsetRankOneSolver<T> ext_solver(EXTsum_T);
    for (auto &eAdj : eEnt.ListAdj) {
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
      insert(*found);
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
BuildIsoKDelaunayDomain(KDelaunayTesselation<Tint, Tgroup> const &DT,
                        [[maybe_unused]] MyMatrix<T> const &GramMat,
                        LinSpaceMatrix<T> const &LinSpa,
                        std::vector<std::vector<Tint>> const &ListGramRing,
                        std::ostream &os) {
  MyMatrix<Tint> ListIneq =
      ComputeIsoKDelaunayInequalities<T, Tint, Tgroup>(DT, ListGramRing, os);
#ifdef SANITY_CHECK_ISO_K_DELAUNAY
  // The generating form is interior to its own domain: every inequality is
  // strictly positive at it.
  {
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
  MyMatrix<Tint> M_ring = RemoveFractionMatrixPlusCoeffRing(M).TheMat;
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyFullRank<Tint, Tint>(M_ring, os);
  return {DT, std::move(ListIneq), std::move(M_ring), std::move(SHV)};
}

template <typename T, typename Tint, typename Tgroup>
struct DataIsoKDelaunayDomains {
  DataIsoDelaunayDomains<T, Tint, Tgroup> data;
  int k;
  std::vector<std::vector<Tint>> ListGramRing;
};

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

template <typename T, typename Tint, typename Tgroup>
struct IsoKDelaunayDomain_AdjI {
  MyVector<Tint> V;
  IsoKDelaunayDomain<T, Tint, Tgroup> dom;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, IsoKDelaunayDomain_AdjI<T, Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("V", eRec.V);
  ar &make_nvp("dom", eRec.dom);
}
} // namespace boost::serialization

template <typename T, typename Tint, typename Tgroup>
struct IsoKDelaunayDomain_Obj {
  IsoKDelaunayDomain<T, Tint, Tgroup> dom;
  // The irredundant inequalities and the permutation group induced on them
  // by the stabilizer of the domain, filled by f_adj.
  MyMatrix<Tint> ListIneqRed;
  Tgroup GRPperm;
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
  WriteEntryGAP(os_out, ent.dom);
  os_out << ", ListIneqRed:=" << StringMatrixGAP(ent.ListIneqRed)
         << ", GRPperm:=" << ent.GRPperm.GapString() << ")";
}

template <typename T, typename Tint, typename Tgroup>
void WriteEntryPYTHON(std::ostream &os_out,
                      IsoKDelaunayDomain_Obj<T, Tint, Tgroup> const &ent) {
  os_out << "{\"dom\":";
  WriteEntryPYTHON(os_out, ent.dom);
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
  os_out << ", n_ineq:=" << ent.dom.ListIneq.rows();
  os_out << ", n_ineq_red:=" << ent.ListIneqRed.rows();
  os_out << ", det:=" << DeterminantMat(ent.dom.GramMat);
  os_out << ", n_shv:=" << ent.dom.SHV.rows();
  os_out << ", n_tile:=" << ent.dom.DT.l_tiles.size();
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
  std::map<std::pair<int, int>, size_t> map_tile;
  for (auto &eEnt : ent.dom.DT.l_tiles) {
    map_tile[{eEnt.tile.INT.rows(), eEnt.tile.EXT.rows()}] += 1;
  }
  os_out << ", ListTileType:=[";
  bool IsFirst = true;
  for (auto &kv : map_tile) {
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
    IsoKDelaunayDomain<T, Tint, Tgroup> dom =
        FlipIsoKDelaunayDomain<T, Tint, Tgroup>(data, TestPt, l_ineq[i]);
#ifdef TIMINGS_ISO_K_DELAUNAY
    os << "|ISO_K_DELAUNAY: f_adj, flip|=" << time << "\n";
#endif
    l_adj.push_back({eIneq, std::move(dom)});
  }
  return {std::move(FACred), std::move(GRPperm), std::move(l_adj)};
}

template <typename T, typename Tint, typename Tgroup>
struct DataIsoKDelaunayDomainsFunc {
  DataIsoKDelaunayDomains<T, Tint, Tgroup> data;
  using Tobj = IsoKDelaunayDomain_Obj<T, Tint, Tgroup>;
  using TadjI = IsoKDelaunayDomain_AdjI<T, Tint, Tgroup>;
  using TadjO = IsoKDelaunayDomain_AdjO<Tint>;
  std::ostream &get_os() { return data.data.rddo.os; }
  Tobj f_init() {
    IsoKDelaunayDomain<T, Tint, Tgroup> dom = GetInitialIsoKDelaunayDomain(data);
    return {std::move(dom), {}, {}};
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
    Tobj x_ret{x.dom, {}, {}};
    MyMatrix<Tint> eBigMat = IdentityMat<Tint>(data.data.LinSpa.n);
    TadjO ret{x.V, eBigMat};
    return {std::move(x_ret), std::move(ret)};
  }
  std::optional<std::vector<TadjI>> f_adj(Tobj &x_in) {
    ResultIsoKDelaunayAdj<T, Tint, Tgroup> result =
        get_result_iso_k_delaunay_adj(x_in.dom, data);
    x_in.ListIneqRed = result.ListIneqRed;
    x_in.GRPperm = result.GRPperm;
    return result.l_adj;
  }
  Tobj f_adji_obj(TadjI const &x) { return {x.dom, {}, {}}; }
  size_t f_complexity([[maybe_unused]] Tobj const &x) { return 0; }
};

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
    std::map<std::string, int> ListIntValues;
    ListIntValues["k"] = 2;
    SingleBlock BlockDATA;
    BlockDATA.setListStringValues(ListStringValues);
    BlockDATA.setListIntValues(ListIntValues);
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
  return {std::move(data), k, std::move(ListGramRing)};
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

// clang-format off
#endif  // SRC_K_COVERINGS_ISOKDELAUNAYDOMAINS_H_
// clang-format on

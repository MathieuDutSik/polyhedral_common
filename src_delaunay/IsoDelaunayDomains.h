// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ISODELAUNAYDOMAINS_H_
#define SRC_LATT_ISODELAUNAYDOMAINS_H_

// clang-format off
#include "POLY_LinearProgramming.h"
#include "POLY_RedundancyElimination.h"
#include "POLY_DirectDualDesc_TS.h"
#include "Shvec_exact.h"
#include "Positivity.h"
#include "LatticeDelaunay.h"
#include "Tspace_General.h"
#include "SystemNamelist.h"
#include <format>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <limits>
// clang-format on

#ifdef DEBUG
#define DEBUG_ISO_DELAUNAY_DOMAIN
#define DEBUG_ISO_DELAUNAY_DOMAIN_PRESENTATION
#endif

#ifdef DISABLE_DEBUG_ISO_DELAUNAY_DOMAIN
#undef DEBUG_ISO_DELAUNAY_DOMAIN
#endif

#ifdef TIMINGS
#define TIMINGS_ISO_DELAUNAY_DOMAIN
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ISO_DELAUNAY_DOMAIN
#define SANITY_CHECK_MULTIPLE_INEQUALITY
#define SANITY_CHECK_DELAUNAY_INTERSECTION
#define SANITY_CHECK_DELAUNAY_NO_EQUALITY
#define SANITY_CHECK_FLIP_CASE1_INEQ
#define SANITY_CHECK_FLIP_INEQ_CONSISTENCY
#endif

template <typename T, typename Tint, typename Tgroup>
struct DataIsoDelaunayDomains {
  // The Linear space of matrices
  LinSpaceMatrix<T> LinSpa;
  // The same space over the ring, built once. The stabilizer, the
  // equivalence and the invariant of a Gram matrix are computed from it:
  // their inner loop is the pairwise scalar products of an integral vector
  // family, which costs far less over the ring.
  LinSpaceMatrix<Tint> LinSpaRing;
  // The database of dual descriptions.
  RecordDualDescOperation<T, Tgroup> rddo;
  // The matrices common to all the IsoDelaunay domains if one is chosen
  // to be common.
  std::optional<MyMatrix<T>> CommonGramMat;
  // The same over the ring, for the stabilizer / equivalence / invariant.
  // Those read it as the condition P G P^T = G, which a positive scaling of
  // G leaves unchanged.
  std::optional<MyMatrix<Tint>> CommonGramMatRing;
};

// The primitive integral representative of an imposed common Gram matrix.
template <typename Tint, typename T>
std::optional<MyMatrix<Tint>>
GetCommonGramMatRing(std::optional<MyMatrix<T>> const &CommonGramMat) {
  if (!CommonGramMat) {
    return {};
  }
  return RemoveFractionMatrixPlusCoeffRing(*CommonGramMat).TheMat;
}

/*
  Code for the L-type domains.

  Two main use case:
  ---Lattice case: Then the Tvert is actually a Tint and can be mpz_class,
  int32_t, etc.
  ---Periodic structure case: Then the coordinates are no longer integral.
    Also the equivalence are no longer integral. Sure the matrix transformation
  is integral, but the translation vector is not necessarily so.

  We use the definitions from LatticeDelaunay.h
  They are for the lattice case but can be generalized for the periodic case.
 */

/*
  Precompute for the Voronoi inequalities of a Delaunay polytope with
  vertex matrix EXT. The basis submatrix VertBasis is inverted through its
  adjugate and determinant so that the whole machinery runs over a ring
  (the intended instantiation: the vertices are integral, possibly after
  the periodic scaling): with the sign flipped to make VertBasisDet
  positive, the affine coordinates of a vertex v with respect to the rows
  of VertBasis are (VertBasisAdjTr * v) / VertBasisDet.
  VoronoiLinearInequality clears that denominator, so its output is the
  Voronoi inequality multiplied by VertBasisDet > 0, which none of its
  uses see (zero test, sign test, scalar canonicalization).
 */
template <typename T> struct VoronoiInequalityPreComput {
  int n;
  Face f_basis;
  MyMatrix<T> VertBasis;
  MyMatrix<T> VertBasisAdjTr;
  T VertBasisDet;
  std::vector<MyVector<T>> VertBasisRed;
};

/*
  We implement the hash of a Delaunay tessellationn. The constraint is that two
  different but equivalent tessellations, must have the same hash. Hopefully,
  this is not a problem since we have many invariants:
  * The number of vertices of the orbit representatives of Delaunay polytopes.
  * The size of their automorphism groups
  * The number of vertices of the orbit representative of their facets.
 */
template <typename T, typename Tgroup>
size_t ComputeInvariantDelaunayTessellation(
    DelaunayTesselation<T, Tgroup> const &DT, size_t const &seed,
    [[maybe_unused]] std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
  std::map<size_t, size_t> map;
  auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
    seed ^= new_hash + 0x9e3779b8 + (seed << 6) + (seed >> 2);
  };
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_PRESENTATION
  std::unordered_map<int, size_t> map_vert;
  std::unordered_map<TintGroup, size_t> map_ord_grp;
#endif

  for (auto &e_del : DT.l_dels) {
    std::map<size_t, size_t> map_siz;
    for (auto &eAdj : e_del.ListAdj) {
      size_t siz = eAdj.eInc.count();
      map_siz[siz] += 1;
    }
    size_t hash_del = 123;
    int n_ext = e_del.EXT.rows();
    TintGroup ord_grp = e_del.GRP.order();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_PRESENTATION
    map_vert[n_ext] += 1;
    map_ord_grp[ord_grp] += 1;
#endif
    size_t hash1 = std::hash<int>()(n_ext);
    size_t hash2 = std::hash<TintGroup>()(ord_grp);
    combine_hash(hash_del, hash1);
    combine_hash(hash_del, hash2);
    for (auto &kv : map_siz) {
      combine_hash(hash_del, kv.first);
      combine_hash(hash_del, kv.second);
    }
    map[hash_del] += 1;
  }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_PRESENTATION
  std::string str_vert = "ISODEL: map_vert =";
  for (auto &kv : map_vert) {
    str_vert +=
        " [" + std::to_string(kv.first) + "," + std::to_string(kv.second) + "]";
  }
  std::string str_ord_grp = "ISODEL: map_ord_grp =";
  for (auto &kv : map_ord_grp) {
    str_ord_grp +=
        " [" + std::format("{}", kv.first) + "," + std::to_string(kv.second) +
        "]";
  }
  os << str_vert << "\n";
  os << str_ord_grp << "\n";
#endif
  size_t hash_ret = seed;
  for (auto &kv : map) {
    combine_hash(hash_ret, kv.first);
    combine_hash(hash_ret, kv.second);
  }
  return hash_ret;
}

template <typename T>
VoronoiInequalityPreComput<T>
BuildVoronoiIneqPreCompute(MyMatrix<T> const &EXT,
                           [[maybe_unused]] std::ostream &os) {
  int n = EXT.cols() - 1;
  SelectionRowCol<T> eSelect = TMat_SelectRowColRing(EXT);
  std::vector<int> const &ListRowSelect = eSelect.ListRowSelect;
  Face f_basis(EXT.rows());
  for (auto &eIdx : ListRowSelect) {
    f_basis[eIdx] = 1;
  }
  MyMatrix<T> VertBasis = SelectRow(EXT, ListRowSelect);
  std::pair<MyMatrix<T>, T> pair = AdjugateDeterminant(VertBasis);
  MyMatrix<T> VertBasisAdjTr = TransposedMat(pair.first);
  T VertBasisDet = pair.second;
  if (VertBasisDet < 0) {
    VertBasisAdjTr = -VertBasisAdjTr;
    VertBasisDet = -VertBasisDet;
  }
  std::vector<MyVector<T>> VertBasisRed;
  for (int i = 0; i <= n; i++) {
    MyVector<T> V(n);
    for (int u = 0; u < n; u++) {
      V(u) = VertBasis(i, u + 1);
    }
    VertBasisRed.push_back(V);
  }
  return {n,
          f_basis,
          std::move(VertBasis),
          std::move(VertBasisAdjTr),
          std::move(VertBasisDet),
          std::move(VertBasisRed)};
}

template <typename T>
MyVector<T> VoronoiLinearInequality(VoronoiInequalityPreComput<T> const &vipc,
                                    MyVector<T> const &TheVert,
                                    std::vector<std::vector<T>> const &ListGram,
                                    [[maybe_unused]] std::ostream &os) {
  int n = vipc.n;
  int dimSpace = ListGram.size();
  // B is VertBasisDet times the affine coordinates of TheVert in the rows
  // of VertBasis, so the inequality comes out multiplied by the positive
  // VertBasisDet and no division occurs.
  MyVector<T> B = vipc.VertBasisAdjTr * TheVert;
  MyVector<T> Ineq(dimSpace);
  MyVector<T> TheVertRed(n);
  for (int u = 0; u < n; u++) {
    TheVertRed(u) = TheVert(u + 1);
  }
  int iGram = 0;
  for (auto &eLineMat : ListGram) {
    T val = vipc.VertBasisDet * EvaluateLineVector(eLineMat, TheVertRed);
    for (int k = 0; k <= n; k++) {
      SubMul(val, B(k), EvaluateLineVector(eLineMat, vipc.VertBasisRed[k]));
    }
    Ineq(iGram) = val;
    iGram++;
  }
  return Ineq;
}

/*
  The list of Gram matrices (as line vectors) over the ring. The whole
  family is scaled by one common positive factor clearing all the
  denominators at once: a uniform positive scaling multiplies every
  Voronoi inequality by that factor, which none of their uses see (zero
  tests, sign tests, scalar canonicalization).
 */
template <typename T>
std::vector<std::vector<typename underlying_ring<T>::ring_type>>
GetListGramRing(std::vector<std::vector<T>> const &ListGram) {
  using Tring = typename underlying_ring<T>::ring_type;
  size_t len = 0;
  for (auto &eLine : ListGram) {
    len += eLine.size();
  }
  MyVector<T> V(len);
  size_t pos = 0;
  for (auto &eLine : ListGram) {
    for (auto &val : eLine) {
      V(pos) = val;
      pos++;
    }
  }
  MyVector<Tring> Vring = RemoveFractionVectorPlusCoeffRing(V).TheVect;
  size_t n_gram = ListGram.size();
  std::vector<std::vector<Tring>> ListGramRing(n_gram);
  pos = 0;
  for (size_t i = 0; i < n_gram; i++) {
    size_t e_len = ListGram[i].size();
    std::vector<Tring> eLine(e_len);
    for (size_t u = 0; u < e_len; u++) {
      eLine[u] = Vring(pos);
      pos++;
    }
    ListGramRing[i] = std::move(eLine);
  }
  return ListGramRing;
}

/*
  Whether some vertex outside the affine basis has a non-zero Voronoi
  regulator, i.e. the polytope forces an equality on the T-space. The
  vertices are over the ring (in the periodic case, in the scaled frame,
  which multiplies the regulators by a positive factor and so changes
  nothing to the zero test); the T-space enters through its ring line
  vectors (GetListGramRing), converted once by the caller.
 */
template <typename Tint>
bool IsDelaunayPolytopeInducingEqualities(
    MyMatrix<Tint> const &EXT,
    std::vector<std::vector<Tint>> const &ListGramRing, std::ostream &os) {
  int n_row = EXT.rows();
  VoronoiInequalityPreComput<Tint> vipc =
      BuildVoronoiIneqPreCompute<Tint>(EXT, os);
  for (int i_row = 0; i_row < n_row; i_row++) {
    if (vipc.f_basis[i_row] == 0) {
      MyVector<Tint> TheVert = GetMatrixRow(EXT, i_row);
      MyVector<Tint> V =
          VoronoiLinearInequality(vipc, TheVert, ListGramRing, os);
      if (!IsZeroVector(V)) {
        return true;
      }
    }
  }
  return false;
}

/*
  Whether every Voronoi inequality of the Delaunay polytope EXT accepts the
  Gram matrix whose T-space coordinate vector is TestV (i.e. the matrix is
  in the closed iso-Delaunay domain on the side of EXT). The ring
  inequalities are converted to the field only for the sign test against
  the rational TestV.
 */
template <typename Tint, typename T>
bool IsDelaunayAcceptableForGramMat(
    MyMatrix<Tint> const &EXT,
    std::vector<std::vector<Tint>> const &ListGramRing,
    MyVector<T> const &TestV, std::ostream &os) {
  int n_row = EXT.rows();
  VoronoiInequalityPreComput<Tint> vipc =
      BuildVoronoiIneqPreCompute<Tint>(EXT, os);
  for (int i_row = 0; i_row < n_row; i_row++) {
    if (vipc.f_basis[i_row] == 0) {
      MyVector<Tint> TheVert = GetMatrixRow(EXT, i_row);
      MyVector<Tint> V =
          VoronoiLinearInequality(vipc, TheVert, ListGramRing, os);
      MyVector<T> V_T = UniversalVectorConversion<T, Tint>(V);
      T TheScal = V_T.dot(TestV);
      if (TheScal < 0) {
        return false;
      }
    }
  }
  return true;
}

struct AdjInfo {
  int iOrb;
  int i_adj;
};

void WriteEntryGAP(std::ostream &os_out, AdjInfo const &eai) {
  os_out << "rec(iOrb:=" << (eai.iOrb + 1) << ", i_adj:=" << (eai.i_adj + 1)
         << ")";
}

void WriteEntryPYTHON(std::ostream &os_out, AdjInfo const &eai) {
  os_out << "{\"iOrb\":" << (eai.iOrb + 1) << ", \"i_adj\":" << (eai.i_adj + 1)
         << "}";
}

namespace boost::serialization {
template <class Archive>
inline void serialize(Archive &ar, AdjInfo &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("iOrb", eRec.iOrb);
  ar &make_nvp("i_adj", eRec.i_adj);
}
} // namespace boost::serialization

// The defining inequality (a primitive integral vector, hence over the
// ring) together with the adjacencies of the tessellation realizing it.
template <typename Tint> struct FullAdjInfo {
  MyVector<Tint> eIneq;
  std::vector<AdjInfo> ListAdjInfo;
};

template <typename T>
void WriteEntryGAP(std::ostream &os_out, FullAdjInfo<T> const &ent) {
  os_out << "rec(eIneq:=" << StringVectorGAP(ent.eIneq) << ", ListAdjInfo:=[";
  bool IsFirst = true;
  for (auto &eAdj : ent.ListAdjInfo) {
    if (!IsFirst) {
      os_out << ",";
    }
    IsFirst = false;
    WriteEntryGAP(os_out, eAdj);
  }
  os_out << "])";
}

template <typename T>
void WriteEntryPYTHON(std::ostream &os_out, FullAdjInfo<T> const &ent) {
  os_out << "{\"eIneq\":" << StringVectorPYTHON(ent.eIneq)
         << ", \"ListAdjInfo\":[";
  bool IsFirst = true;
  for (auto &eAdj : ent.ListAdjInfo) {
    if (!IsFirst) {
      os_out << ",";
    }
    IsFirst = false;
    WriteEntryPYTHON(os_out, eAdj);
  }
  os_out << "]}";
}

namespace boost::serialization {
template <class Archive, typename T>
inline void serialize(Archive &ar, FullAdjInfo<T> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("eIneq", eRec.eIneq);
  ar &make_nvp("ListAdjInfo", eRec.ListAdjInfo);
}
} // namespace boost::serialization

/*
  A DelaunayTesselation annotated with, on top of each Delaunay_AdjO, the
  (ScalarCanonicalizationVector-reduced) Voronoi inequality of that
  adjacency. Carrying the inequalities alongside the tessellation lets a
  freshly flipped tessellation reuse the inequalities of the edges that
  survive the flip unchanged instead of recomputing every inequality of
  the whole tessellation from scratch (see the Case 1 handling in
  FlippingLtype).
 */
template <typename Tint> struct Delaunay_AdjIneqO {
  Face eInc;
  MyMatrix<Tint> eBigMat;
  int iOrb;
  // The (ScalarCanonicalizationVector-reduced) Voronoi inequality; its
  // canonical representative is primitive, hence stored over the ring.
  MyVector<Tint> eIneq;
};

template <typename Tint, typename Tgroup>
struct Delaunay_EntryIneq {
  MyMatrix<Tint> EXT;
  Tgroup GRP;
  std::vector<Delaunay_AdjIneqO<Tint>> ListAdj;
};

template <typename Tint, typename Tgroup>
struct DelaunayTesselationIneq {
  std::vector<Delaunay_EntryIneq<Tint, Tgroup>> l_dels;
};

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, Delaunay_AdjIneqO<Tint> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("eInc", eRec.eInc);
  ar &make_nvp("eBigMat", eRec.eBigMat);
  ar &make_nvp("iOrb", eRec.iOrb);
  ar &make_nvp("eIneq", eRec.eIneq);
}
template <class Archive, typename Tint, typename Tgroup>
inline void serialize(Archive &ar,
                      Delaunay_EntryIneq<Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("EXT", eRec.EXT);
  ar &make_nvp("GRP", eRec.GRP);
  ar &make_nvp("ListAdj", eRec.ListAdj);
}
template <class Archive, typename Tint, typename Tgroup>
inline void serialize(Archive &ar,
                      DelaunayTesselationIneq<Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("l_dels", eRec.l_dels);
}
} // namespace boost::serialization

// The GAP / PYTHON records do not carry the per-adjacency inequality, so the
// Ineq tessellation is stripped to the plain tessellation and the shared
// writers are used on that.
template <typename Tint, typename Tgroup>
void WriteEntryGAP(std::ostream &os_out,
                   DelaunayTesselationIneq<Tint, Tgroup> const &DT) {
  DelaunayTesselation<Tint, Tgroup> DT_strip = StripDelaunayTesselationIneq(DT);
  WriteEntryGAP_l_dels<Tint, Tgroup>(os_out, DT_strip.l_dels);
}

template <typename Tint, typename Tgroup>
void WriteEntryPYTHON(std::ostream &os_out,
                      DelaunayTesselationIneq<Tint, Tgroup> const &DT) {
  DelaunayTesselation<Tint, Tgroup> DT_strip = StripDelaunayTesselationIneq(DT);
  WriteEntryPYTHON_l_dels<Tint, Tgroup>(os_out, DT_strip.l_dels);
}

template <typename Tint, typename Tgroup>
DelaunayTesselation<Tint, Tgroup>
StripDelaunayTesselationIneq(
    DelaunayTesselationIneq<Tint, Tgroup> const &DTI) {
  std::vector<Delaunay_Entry<Tint, Tgroup>> l_dels;
  for (auto &eDel : DTI.l_dels) {
    std::vector<Delaunay_AdjO<Tint>> ListAdj;
    for (auto &eAdj : eDel.ListAdj) {
      ListAdj.push_back({eAdj.eInc, eAdj.eBigMat, eAdj.iOrb});
    }
    l_dels.push_back({eDel.EXT, eDel.GRP, std::move(ListAdj)});
  }
  return {std::move(l_dels)};
}

/*
  Build the Voronoi inequality precompute for a Delaunay polytope EXT.
 */
template <typename T>
VoronoiInequalityPreComput<T> BuildVoronoiIneqPreComputeChecked(
    MyMatrix<T> const &EXT,
    [[maybe_unused]] std::vector<std::vector<T>> const &ListGram,
    std::ostream &os) {
  VoronoiInequalityPreComput<T> vipc = BuildVoronoiIneqPreCompute<T>(EXT, os);
#ifdef SANITY_CHECK_DELAUNAY_NO_EQUALITY
  // A Delaunay polytope of a generic form of the T-space must be co-spherical
  // throughout the T-space: every non-basis vertex must give the zero Voronoi
  // regulator. A non-zero one means the polytope is co-spherical only on a
  // wall, i.e. the form lies on an iso-Delaunay wall (it is not generic) --
  // exactly what makes the defining inequalities ill-defined.
  {
    int n_row_cur = EXT.rows();
    for (int i_row = 0; i_row < n_row_cur; i_row++) {
      if (vipc.f_basis[i_row] == 0) {
        MyVector<T> TheVert = GetMatrixRow(EXT, i_row);
        MyVector<T> V = VoronoiLinearInequality(vipc, TheVert, ListGram, os);
        if (!IsZeroVector(V)) {
          std::cerr << "DELAUNAY_NO_EQUALITY: SANITY_CHECK failed: Delaunay "
                       "polytope induces an equality in the T-space "
                       "(non-basis vertex regulator V="
                    << StringVectorGAP(V) << " is non-zero); the form is on "
                       "an iso-Delaunay wall, not generic\n";
          throw TerminalException{1};
        }
      }
    }
  }
#endif
  return vipc;
}

/*
  Compute the (ScalarCanonicalizationVector-reduced) Voronoi inequality of
  the adjacency taking the polytope described by (vipc, cont) to EXT_target
  via eBigMat. Shared by BuildDelaunayTesselationIneq (whole tessellation,
  from scratch) and FlippingLtype (single adjacency, after a flip).
 */
template <typename T>
MyVector<T> ComputeDelaunayAdjIneq(VoronoiInequalityPreComput<T> const &vipc,
                                   ContainerMatrix<T> const &cont,
                                   MyMatrix<T> const &EXT_target,
                                   MyMatrix<T> const &eBigMat,
                                   std::vector<std::vector<T>> const &ListGram,
                                   std::ostream &os) {
  MyMatrix<T> EXTadj = EXT_target * eBigMat;
  int len = EXTadj.rows();
#ifdef SANITY_CHECK_DELAUNAY_INTERSECTION
  // The current polytope and the adjacent one EXTadj must share a facet, so
  // the vertices common to both span a hyperplane: as homogeneous vectors
  // they must have rank exactly n. A smaller rank means the recorded
  // adjacency is not a genuine facet-neighbour (the Delaunay polytopes /
  // their adjacencies were not computed correctly).
  {
    std::vector<MyVector<T>> common;
    for (int u = 0; u < len; u++) {
      MyVector<T> TheVert = GetMatrixRow(EXTadj, u);
      std::optional<size_t> opt = cont.GetIdx_v(TheVert);
      if (opt) {
        common.push_back(TheVert);
      }
    }
    int n = EXTadj.cols() - 1;
    MyMatrix<T> Mcommon = MatrixFromVectorFamily(common);
    int rnk = RankMat(Mcommon);
    if (rnk != n) {
      std::cerr << "DELAUNAY_INTERSECTION: SANITY_CHECK failed: the vertices "
                   "common to a Delaunay polytope and its recorded "
                   "neighbour have rank " << rnk << " but a facet requires "
                   "rank n=" << n << " (|common|=" << common.size() << ")\n";
      throw TerminalException{1};
    }
  }
#endif
#ifdef SANITY_CHECK_MULTIPLE_INEQUALITY
  // The vertices of the adjacent Delaunay polytope that are not in the
  // current polytope all lie on the far side of the shared facet, so the
  // Voronoi inequalities they induce must all be scalar multiples of each
  // other (they express the single condition that the shared facet stays
  // Delaunay). After ScalarCanonicalizationVector they must therefore all
  // be equal: the set below must have size 1.
  {
    std::set<MyVector<T>> set_ineq;
    for (int u = 0; u < len; u++) {
      MyVector<T> TheVert = GetMatrixRow(EXTadj, u);
      std::optional<size_t> opt = cont.GetIdx_v(TheVert);
      if (!opt) {
        MyVector<T> V = VoronoiLinearInequality(vipc, TheVert, ListGram, os);
        set_ineq.insert(ScalarCanonicalizationVector(V));
      }
    }
    if (set_ineq.size() != 1) {
      std::cerr << "MULTIPLE_INEQUALITY: SANITY_CHECK failed: the Voronoi "
                   "inequalities from the vertices of an adjacent Delaunay "
                   "polytope are not all multiples of each other (canonical "
                   "set size=" << set_ineq.size() << ")\n";
      for (auto &eV : set_ineq) {
        std::cerr << "MULTIPLE_INEQUALITY:   canon=" << StringVectorGAP(eV)
                  << "\n";
      }
      throw TerminalException{1};
    }
  }
#endif
  for (int u = 0; u < len; u++) {
    MyVector<T> TheVert = GetMatrixRow(EXTadj, u);
    std::optional<size_t> opt = cont.GetIdx_v(TheVert);
    if (!opt) {
      MyVector<T> V = VoronoiLinearInequality(vipc, TheVert, ListGram, os);
      return ScalarCanonicalizationVector(V);
    }
  }
  std::cerr << "Failed to find a matching entry in ComputeDelaunayAdjIneq\n";
  throw TerminalException{1};
}

/*
  Build a DelaunayTesselationIneq from a DelaunayTesselation: same orbits
  and adjacencies, each adjacency additionally carrying its inequality.
 */
template <typename Tint, typename Tgroup>
DelaunayTesselationIneq<Tint, Tgroup> BuildDelaunayTesselationIneq(
    DelaunayTesselation<Tint, Tgroup> const &DT,
    std::vector<std::vector<Tint>> const &ListGramRing, std::ostream &os) {
  int n_del = DT.l_dels.size();
  std::vector<Delaunay_EntryIneq<Tint, Tgroup>> l_dels(n_del);
  for (int i_del = 0; i_del < n_del; i_del++) {
    Delaunay_Entry<Tint, Tgroup> const &eDel = DT.l_dels[i_del];
    VoronoiInequalityPreComput<Tint> vipc =
        BuildVoronoiIneqPreComputeChecked<Tint>(eDel.EXT, ListGramRing, os);
    ContainerMatrix<Tint> cont(eDel.EXT);
    int n_adj = eDel.ListAdj.size();
    std::vector<Delaunay_AdjIneqO<Tint>> ListAdjIneq(n_adj);
    for (int i_adj = 0; i_adj < n_adj; i_adj++) {
      Delaunay_AdjO<Tint> const &adj = eDel.ListAdj[i_adj];
      MyVector<Tint> eIneq = ComputeDelaunayAdjIneq(
          vipc, cont, DT.l_dels[adj.iOrb].EXT, adj.eBigMat, ListGramRing, os);
      ListAdjIneq[i_adj] = {adj.eInc, adj.eBigMat, adj.iOrb,
                            std::move(eIneq)};
    }
    l_dels[i_del] = {eDel.EXT, eDel.GRP, std::move(ListAdjIneq)};
  }
  return {std::move(l_dels)};
}

/*
  Collect the (already computed) per-adjacency inequalities of a
  DelaunayTesselationIneq into the grouped-by-canonical-inequality form
  used to build the defining system of the iso-Delaunay domain.
 */
template <typename Tint, typename Tgroup>
std::vector<FullAdjInfo<Tint>> ComputeListIneqFromTesselationIneq(
    DelaunayTesselationIneq<Tint, Tgroup> const &DTI) {
  std::unordered_map<MyVector<Tint>, std::vector<AdjInfo>> map;
  int n_del = DTI.l_dels.size();
  for (int i_del = 0; i_del < n_del; i_del++) {
    int n_adj = DTI.l_dels[i_del].ListAdj.size();
    for (int i_adj = 0; i_adj < n_adj; i_adj++) {
      MyVector<Tint> const &V_red = DTI.l_dels[i_del].ListAdj[i_adj].eIneq;
      AdjInfo eAdj{i_del, i_adj};
      map[V_red].push_back(eAdj);
    }
  }
  std::vector<FullAdjInfo<Tint>> l_ret;
  for (auto &kv : map) {
    l_ret.push_back({kv.first, kv.second});
  }
  return l_ret;
}

/*
  Compute the defining inequalities of an iso-Delaunay domain
 */
template <typename Tint, typename Tgroup>
std::vector<FullAdjInfo<Tint>> ComputeDefiningIneqIsoDelaunayDomain(
    DelaunayTesselation<Tint, Tgroup> const &DT,
    std::vector<std::vector<Tint>> const &ListGramRing, std::ostream &os) {
  DelaunayTesselationIneq<Tint, Tgroup> DTI =
      BuildDelaunayTesselationIneq(DT, ListGramRing, os);
  return ComputeListIneqFromTesselationIneq(DTI);
}

template <typename T>
MyMatrix<T> GetFACineq(std::vector<FullAdjInfo<T>> const &ListIneq) {
  int dim_spa = ListIneq[0].eIneq.size();
  int n_ineq = ListIneq.size();
  MyMatrix<T> FAC(n_ineq, dim_spa);
  for (int i_ineq = 0; i_ineq < n_ineq; i_ineq++) {
    for (int u = 0; u < dim_spa; u++) {
      FAC(i_ineq, u) = ListIneq[i_ineq].eIneq(u);
    }
  }
  return FAC;
}

template <typename Tint, typename Tgroup>
MyMatrix<Tint> get_list_ineq(
    DelaunayTesselationIneq<Tint, Tgroup> const &DTI,
    [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  MicrosecondTime time_interior;
#endif
  std::vector<FullAdjInfo<Tint>> ListIneq =
      ComputeListIneqFromTesselationIneq(DTI);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: get_list_ineq, ComputeListIneqFromTesselationIneq|="
     << time_interior << "\n";
#endif
  MyMatrix<Tint> FAC = GetFACineq(ListIneq);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: get_list_ineq, GetFACineq|=" << time_interior << "\n";
#endif
  return FAC;
}

template <typename T>
MyMatrix<T> get_interior_gram_matrix_lp(LinSpaceMatrix<T> const &LinSpa,
                                        MyMatrix<T> const &FAC,
                                        std::ostream &os) {
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  MicrosecondTime time;
#endif
  int n = LinSpa.n;
  int dimSpace = LinSpa.ListMat.size();
  MyVector<T> ThePt = GetGeometricallyUniqueInteriorPoint(FAC, os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: get_interior_gram_matrix_lp, "
        "GetGeometricallyUniqueInteriorPoint|="
     << time << "\n";
#endif
  MyMatrix<T> RetMat = ZeroMatrix<T>(n, n);
  for (int u = 0; u < dimSpace; u++) {
    MatAddMul(RetMat, ThePt(u), LinSpa.ListMat[u]);
  }
  return RetMat;
}

template <typename T>
MyMatrix<T> get_interior_gram_matrix_lrs(LinSpaceMatrix<T> const &LinSpa,
                                         MyMatrix<T> const &FAC,
                                         std::ostream &os) {
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  MicrosecondTime time;
#endif
  int n = LinSpa.n;
  int dimSpace = LinSpa.ListMat.size();
  MyMatrix<T> EXT = DirectDualdescription_mat(FAC, os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: get_interior_gram_matrix_lrs, DirectDualdescription|="
     << time << "\n";
#endif
  int n_row = EXT.rows();
  MyMatrix<T> SumMatExtRay = ZeroMatrix<T>(n, n);
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyMatrix<T> RayMat = ZeroMatrix<T>(n, n);
    for (int u = 0; u < dimSpace; u++) {
      MatAddMul(RayMat, EXT(i_row, u), LinSpa.ListMat[u]);
    }
    SumMatExtRay += RemoveFractionMatrix(RayMat);
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN
    DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(RayMat, os);
    if (DiagInfo.nbMinus > 0) {
      std::cerr
          << "One ray of the IsoDelaunay domain is not positive semidefinite\n";
      std::cerr << "This is not allowed\n";
      throw TerminalException{1};
    }
#endif
  }
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: get_interior_gram_matrix_lrs, SumMatExtRay|=" << time << "\n";
#endif
  return SumMatExtRay;
}

template <typename T, typename Tint, typename Tgroup>
MyMatrix<T> GetInteriorGramMatrix(
    LinSpaceMatrix<T> const &LinSpa,
    DelaunayTesselationIneq<Tint, Tgroup> const &DTI,
    std::ostream &os) {
  MyMatrix<Tint> FAC = get_list_ineq(DTI, os);
  MyMatrix<T> FAC_T = UniversalMatrixConversion<T, Tint>(FAC);
  return get_interior_gram_matrix_lp(LinSpa, FAC_T, os);
}

template <typename T, typename Tint, typename Tgroup>
std::pair<MyMatrix<T>, MyMatrix<T>>
GetInteriorGramMatrices(
    LinSpaceMatrix<T> const &LinSpa,
    DelaunayTesselationIneq<Tint, Tgroup> const &DTI,
    std::ostream &os) {
  MyMatrix<Tint> FAC = get_list_ineq(DTI, os);
  MyMatrix<T> FAC_T = UniversalMatrixConversion<T, Tint>(FAC);
  MyMatrix<T> Mat1 = get_interior_gram_matrix_lp(LinSpa, FAC_T, os);
  MyMatrix<T> Mat2 = get_interior_gram_matrix_lrs(LinSpa, FAC_T, os);
  return {std::move(Mat1), std::move(Mat2)};
}

template <typename T, typename Tint, typename Tgroup>
DelaunayTesselation<Tint, Tgroup> GetInitialGenericDelaunayTesselation(
    DataIsoDelaunayDomains<T, Tint, Tgroup> &data) {
  std::ostream &os = data.rddo.os;
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: GetInitialGenericDelaunayTesselation, beginning\n";
#endif
  std::vector<std::vector<Tint>> ListGramRing =
      GetListGramRing(data.LinSpa.ListLineMat);
  auto f_incorrect = [&](Delaunay_Obj<Tint, Tgroup> const &x) -> bool {
    MyMatrix<Tint> const &EXT = x.EXT;
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: Before IsDelaunayPolytopeInducingEqualities\n";
#endif
    bool test1 = IsDelaunayPolytopeInducingEqualities(EXT, ListGramRing, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: After IsDelaunayPolytopeInducingEqualities test1=" << test1
       << "\n";
#endif
    if (test1) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: GetInitialGenericDelaunayTesselation, true by "
            "IsDelaunayPolytopeInducingEqualities\n";
#endif
      return true;
    }
    // Note: when CommonGramMat is imposed, the requirement that the seed domain
    // contains G in its closure is a whole-domain condition, so it is enforced
    // by contains_common below rather than by this per-cell predicate.
    return false;
  };
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: GetInitialGenericDelaunayTesselation, we have f_incorrect\n";
#endif
  auto test_matrix = [&](MyMatrix<T> const &GramMat)
      -> std::optional<DelaunayTesselation<Tint, Tgroup>> {
    bool test =
        IsSymmetryGroupCorrect<T, Tint, Tgroup>(GramMat, data.LinSpa, os);
    if (!test) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: GetInitialGenericDelaunayTesselation, failing by "
            "symmetry\n";
#endif
      return {};
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: GetInitialGenericDelaunayTesselation, symmetry is ok\n";
#endif
    using TintGroup = typename Tgroup::Tint;
    int dimEXT = GramMat.rows() + 1;
    PolyHeuristicSerial<TintGroup> AllArr =
        AllStandardHeuristicSerial<T, TintGroup>(dimEXT, os);
    DataLattice<T, Tint, Tgroup> data_lattice =
        GetDataLattice<T, Tint, Tgroup>(GramMat, AllArr, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: After GetInitialGenericDelaunayTesselation "
          "data_lattice.rddo.AllArr.OutFile="
       << data_lattice.rddo.AllArr.OutFile << "\n";
#endif
    int max_runtime_second = 0;
    std::optional<DelaunayTesselation<Tint, Tgroup>> result =
        EnumerationDelaunayPolytopes<T, Tint, Tgroup, decltype(f_incorrect)>(
            data_lattice, f_incorrect, max_runtime_second);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: GetInitialGenericDelaunayTesselation, return something\n";
#endif
    return result;
  };
  // When a common Gram matrix G is imposed, the seed must be a generic
  // (top-dimensional) iso-Delaunay domain whose closure contains G. G generally
  // lies on the walls shared by several such domains, so we seed near its ray:
  // GramMat = P + scale*G with P a generic random form. Since iso-Delaunay
  // domains are scale-invariant cones, P + scale*G tends to the ray of G as
  // scale grows, hence lands in a domain whose closure contains G. We double
  // scale until the containment test succeeds. Without CommonGramMat this is
  // the ordinary random seed (scale unused, contains_common vacuously true) and
  // we instead grow the random amplitude N to escape a box of finite
  // possibilities.
  std::optional<MyVector<T>> g_vec;
  if (data.CommonGramMat) {
    g_vec = LINSPA_GetVectorOfMatrixExpression(data.LinSpa, *data.CommonGramMat);
  }
  auto contains_common = [&](DelaunayTesselation<Tint, Tgroup> const &DT,
                             MyMatrix<T> const &GramMat) -> bool {
    if (!g_vec) {
      return true;
    }
    std::vector<FullAdjInfo<Tint>> ListIneq =
        ComputeDefiningIneqIsoDelaunayDomain(DT, ListGramRing, os);
    MyVector<T> c_vec = LINSPA_GetVectorOfMatrixExpression(data.LinSpa, GramMat);
    for (auto &eRec : ListIneq) {
      // GramMat is strictly interior, so eIneq.dot(c_vec) != 0 fixes the
      // accepted orientation of each wall; G is in the closure iff it is never
      // on the strictly opposite side.
      MyVector<T> eIneq_T = UniversalVectorConversion<T, Tint>(eRec.eIneq);
      T s_g = eIneq_T.dot(*g_vec);
      T s_c = eIneq_T.dot(c_vec);
      if (s_g * s_c < 0) {
        return false;
      }
    }
    return true;
  };
  size_t n_iter = 0;
  int N = 2;
  T scale(1);
  while (true) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: Before GetRandomPositiveDefiniteNoNontrivialSymm, n_iter="
       << n_iter << " N=" << N << " scale=" << scale << "\n";
#endif
    MyMatrix<T> GramMat =
        GetRandomPositiveDefiniteNoNontrivialSymm<T, Tint, Tgroup>(data.LinSpa,
                                                                   N, os);
    if (data.CommonGramMat) {
      GramMat += scale * (*data.CommonGramMat);
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: After GetRandomPositiveDefiniteNoNontrivialSymm, GramMat=\n";
    WriteMatrix(os, GramMat);
#endif
    std::optional<DelaunayTesselation<Tint, Tgroup>> opt =
        test_matrix(GramMat);
    if (opt && contains_common(*opt, GramMat)) {
      return *opt;
    }
    n_iter += 1;
    scale *= T(2);
    N += n_iter;
  }
  std::cerr << "ISODEL: Failed to find a matching entry in "
               "GetInitialGenericDelaunayTesselation\n";
  throw TerminalException{1};
}

template <typename Tint, typename Tgroup> struct RepartEntry {
  MyMatrix<Tint> EXT;
  Tgroup TheStab;
  int8_t Position; // -1: lower, 0: barrel, 1: higher
  int iDelaunayOrigin;
  std::vector<Delaunay_AdjO<Tint>> ListAdj;
  MyMatrix<Tint> eBigMat;
  // Cache of the inverse of eBigMat, filled at most once by get_bigmat_inv
  // and read through it only. The transformations are unimodular, so the
  // inverse stays over the ring. The adjacency loop of FlippingLtype needs
  // it once per adjacency while it depends only on the cell. The entry stays
  // an aggregate (it is brace-initialized without this field), so the cache
  // is a mutable member rather than a private one.
  mutable std::optional<MyMatrix<Tint>> eBigMatInv;
  MyMatrix<Tint> const &get_bigmat_inv() const {
    if (!eBigMatInv) {
      eBigMatInv = Inverse(eBigMat);
    }
    return *eBigMatInv;
  }
};

template <typename T>
std::vector<MyVector<T>>
Orbit_MatrixGroup(std::vector<MyMatrix<T>> const &ListGen,
                  MyVector<T> const &eV, std::ostream &os) {
  auto f_prod = [](MyVector<T> const &v,
                   MyMatrix<T> const &M) -> MyVector<T> {
    return M.transpose() * v;
  };
  return OrbitComputation<MyMatrix<T>, MyVector<T>, decltype(f_prod)>(
      ListGen, eV, f_prod, os);
}

// eIso is the isobarycenter of the repartitioning polytope scaled by the
// number of vertices, (nVert, sum of the vertex coordinates), which keeps it
// over the ring. Its only use is as a point fixed by the transformations of
// the barrel facets, for which the scaling is immaterial (a row and its
// image scale alike in the solve).
template <typename Tint, typename Tgroup> struct FullRepart {
  std::vector<RepartEntry<Tint, Tgroup>> cells;
  MyVector<Tint> eIso;
};

/*
  The function is named "NextGeneration" since there was a slower version in GAP
  which was superseded by a newer one and the name is inherited.
  What the function do is merge Delaunay polytopes and form the repartitioning
  polytopes by lifting the height.
  The lower facets correspond to the old Delaunay tessellation, the upper ones
  to the new one. The lateral ones on the side are named "barrel" and are used
  when switching groups of repartitioning polytopes simultaneously.
 */
template <typename Tring, typename Tgroup>
FullRepart<Tring, Tgroup> FindRepartitionningInfoNextGeneration(
    size_t eIdx,
    DelaunayTesselationIneq<Tring, Tgroup> const &ListOrbitDelaunay,
    std::vector<AdjInfo> const &ListInformationsOneFlipping,
    MyMatrix<Tring> const &InteriorElement,
    PolyHeuristicSerial<typename Tgroup::Tint> &AllArr, std::ostream &os) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FindRepartitionningInfoNextGeneration, begin FRING\n";
  os << "ISODEL: |ListInformationsOneFlipping|="
     << ListInformationsOneFlipping.size() << "\n";
  for (auto &eEnt : ListInformationsOneFlipping) {
    int iOrbAdj = ListOrbitDelaunay.l_dels[eEnt.iOrb].ListAdj[eEnt.i_adj].iOrb;
    os << "ISODEL:    iOrb=" << eEnt.iOrb << " i_adj=" << eEnt.i_adj
       << " iOrbAdj=" << iOrbAdj << "\n";
  }
#endif
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  // The lifted repartitionning geometry is integral, so its dual
  // descriptions run over the ring. The bank of this record only lives for
  // this repartitionning; the facet polytopes of one repartitionning are
  // the candidates for isomorphy anyway.
  RecordDualDescOperation<Tring, Tgroup> rddo_ring(AllArr, os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  MicrosecondTime time_fring_full;
  int64_t time_fring_group = 0, time_fring_stab = 0, time_fring_dualdesc = 0,
          time_fring_flip = 0, time_fring_insert_facet = 0,
          time_fring_onsets = 0, time_fring_reprvert = 0,
          time_fring_centers = 0, time_fring_vertmat = 0,
          time_fring_kernel = 0, time_fring_ext2 = 0, time_fring_frame = 0,
          time_fring_lev = 0;
#endif
  int n = InteriorElement.rows();
  std::vector<std::vector<Tidx>> ListPermGenList;
  std::vector<MyMatrix<Tring>> ListMatGens;
  std::vector<MyVector<Tring>> ListVertices;
  std::unordered_map<MyVector<Tring>, size_t> ListVertices_rev;
  size_t idx_vertices = 1;
  Tgroup PermGRP;
  // Building PermGRP (a stabilizer-chain computation) from scratch is not
  // cheap, and the vertex/generator discovery loop below can trigger many
  // insertion events in a row. Only the LAST rebuild before PermGRP is
  // actually read (the isin() check below, or the facet-orbit computations
  // in the second half of this function) matters, so the rebuild is done
  // lazily: insertions just mark the group stale via group_dirty.
  bool group_dirty = true;
  auto StandardGroupUpdate = [&]() -> void {
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
    MicrosecondTime time_g;
#endif
    std::vector<Telt> ListGen;
    Tidx n_act = idx_vertices - 1;
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: StandardGroupUpdate n_act="
       << static_cast<size_t>(n_act) << "\n";
#endif
    for (auto &eList : ListPermGenList) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: StandardGroupUpdate |eList|=" << eList.size()
         << "\n";
#endif
      Telt x(eList);
      ListGen.push_back(x);
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: StandardGroupUpdate Before Tgroup |PermGRP|="
       << PermGRP.size() << "\n";
#endif
    PermGRP = Tgroup(ListGen, n_act);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: StandardGroupUpdate After Tgroup |PermGRP|="
       << PermGRP.size() << "\n";
#endif
    group_dirty = false;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
    time_fring_group += time_g.eval_int64();
#endif
  };
  auto ensure_group_updated = [&]() -> void {
    if (group_dirty) {
      StandardGroupUpdate();
    }
  };
  struct TypeOrbitCenter {
    int iDelaunay;
    MyMatrix<Tring> eBigMat;
    std::vector<Tidx> Linc;
    MyMatrix<Tring> EXT;
  };
  struct TypeOrbitCenterMin {
    int iDelaunay;
    MyMatrix<Tring> eBigMat;
  };
  std::vector<TypeOrbitCenter> ListOrbitCenter;
  auto FuncInsertVertex = [&](MyVector<Tring> const &eVert) -> void {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertVertex, step 1\n";
#endif
    if (ListVertices_rev.contains(eVert)) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: is_present eVert=" << StringVector(eVert) << "\n";
#endif
      return;
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertVertex, step 2\n";
#endif
    std::vector<MyVector<Tring>> O = Orbit_MatrixGroup(ListMatGens, eVert, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertVertex, step 3 |O|=" << O.size() << "\n";
#endif
    for (auto &eV : O) {
      ListVertices.push_back(eV);
      ListVertices_rev[eV] = idx_vertices;
      idx_vertices++;
    }
    group_dirty = true;
    size_t n_gen = ListPermGenList.size();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertVertex, step 4, n_gen=" << n_gen << "\n";
#endif
    for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
      std::vector<Tidx> &ePermGen = ListPermGenList[i_gen];
      MyMatrix<Tring> const &eMatrGen = ListMatGens[i_gen];
      for (auto &eVert : O) {
        MyVector<Tring> eVertImg = eMatrGen.transpose() * eVert;
        size_t pos = ListVertices_rev.at(eVertImg);
        Tidx pos_idx = static_cast<Tidx>(pos - 1);
        ePermGen.push_back(pos_idx);
      }
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertVertex, step 5\n";
#endif
  };
  auto get_v_positions =
      [&](MyMatrix<Tring> const &eMat) -> std::optional<std::vector<Tidx>> {
    size_t len = idx_vertices - 1;
    std::vector<Tidx> v_ret(len);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: get_v_positions, len=" << len << "\n";
#endif
    for (size_t i = 0; i < len; i++) {
      MyVector<Tring> eVimg = eMat.transpose() * ListVertices[i];
      if (!ListVertices_rev.contains(eVimg)) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: FRING: get_v_positions, i=" << i
           << " eV=" << StringVector(ListVertices[i])
           << " eVimg=" << StringVector(eVimg) << "\n";
#endif
        return {};
      }
      size_t pos = ListVertices_rev.at(eVimg);
      Tidx pos_idx = static_cast<Tidx>(pos - 1);
      v_ret[i] = pos_idx;
    }
    return v_ret;
  };
  auto FuncInsertGenerator = [&](MyMatrix<Tring> const &eMat) -> void {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertGenerator, step 1\n";
#endif
    std::optional<std::vector<Tidx>> opt = get_v_positions(eMat);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertGenerator, step 2\n";
#endif
    auto get_test_belong = [&]() -> bool {
      if (opt) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: FRING: opt, case exist\n";
#endif
        Telt elt(*opt);
        ensure_group_updated();
        return PermGRP.isin(elt);
      } else {
        std::vector<MyMatrix<Tring>> LGen = ListMatGens;
        LGen.push_back(eMat);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: FRING: opt, case non-exist\n";
        os << "ISODEL: FRING: Before |ListVertices|=" << ListVertices.size()
           << "\n";
#endif
        size_t nVert = ListVertices.size();
        Face f_att(nVert);
        size_t miss_val = std::numeric_limits<size_t>::max();
        auto get_idx = [&](MyVector<Tring> const &eVert) -> size_t {
          if (!ListVertices_rev.contains(eVert)) {
            return miss_val;
          }
          size_t u = ListVertices_rev.at(eVert) - 1;
          if (u >= nVert) {
            return miss_val;
          }
          return u;
        };
        for (size_t pos = 0; pos < nVert; pos++) {
          if (f_att[pos] == 0) {
            MyVector<Tring> const &eVert = ListVertices[pos];
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
            os << "ISODEL: FRING: Treating one vertex pos=" << pos << "\n";
#endif
            std::vector<MyVector<Tring>> TheOrb =
                Orbit_MatrixGroup(LGen, eVert, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
            os << "ISODEL: FRING: We have pos=" << pos
               << " |TheOrb|=" << TheOrb.size() << "\n";
#endif
            for (auto &eVertB : TheOrb) {
              size_t idx = get_idx(eVertB);
              if (idx == miss_val) {
                FuncInsertVertex(eVertB);
              } else {
                f_att[idx] = 1;
              }
            }
          }
        }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: FRING: After |ListVertices|=" << ListVertices.size()
           << "\n";
#endif
        return false;
      }
    };
    bool test = get_test_belong();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertGenerator, step 3, test=" << test << "\n";
#endif
    if (!test) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      for (size_t pos = 0; pos < ListVertices.size(); pos++) {
        os << "ISODEL: FRING: pos=" << pos
           << " eV=" << StringVector(ListVertices[pos]) << "\n";
      }
#endif
      std::optional<std::vector<Tidx>> opt = get_v_positions(eMat);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      if (!opt) {
        std::cerr << "The get_v_positions returned a null\n";
        throw TerminalException{1};
      }
#endif
      ListPermGenList.push_back(*opt);
      ListMatGens.push_back(eMat);
      group_dirty = true;
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertGenerator, step 4\n";
#endif
  };
  auto FuncInsertCenter = [&](TypeOrbitCenterMin const &TheRec) -> void {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertCenter, step 1\n";
#endif
    MyMatrix<Tring> LVert =
        ListOrbitDelaunay.l_dels[TheRec.iDelaunay].EXT * TheRec.eBigMat;
    for (int u = 0; u < LVert.rows(); u++) {
      MyVector<Tring> eVert = GetMatrixRow(LVert, u);
      FuncInsertVertex(eVert);
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertCenter, step 2\n";
#endif
    // No group rebuild here: FuncInsertVertex above already marked the group
    // dirty if it grew, and the only functional read within this call is
    // FuncInsertGenerator's isin() check, which rebuilds lazily itself.
    auto get_iOrbFound = [&]() -> std::optional<size_t> {
      for (size_t iOrb = 0; iOrb < ListOrbitCenter.size(); iOrb++) {
        if (ListOrbitCenter[iOrb].iDelaunay == TheRec.iDelaunay) {
          return iOrb;
        }
      }
      return {};
    };
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertCenter, step 3\n";
#endif
    std::optional<size_t> opt = get_iOrbFound();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertCenter, step 4\n";
#endif
    if (opt) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: FuncInsertCenter, step 5, A\n";
#endif
      MyMatrix<Tring> eGen =
          Inverse(TheRec.eBigMat) * ListOrbitCenter[*opt].eBigMat;
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: FuncInsertCenter, step 5, B\n";
#endif
      FuncInsertGenerator(eGen);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: FuncInsertCenter, step 5, C\n";
#endif
    } else {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: FuncInsertCenter, step 6, A\n";
#endif
      MyMatrix<Tring> const &BigMatR = TheRec.eBigMat;
      MyMatrix<Tring> BigMatI = Inverse(BigMatR);
      MyMatrix<Tring> const &EXT =
          ListOrbitDelaunay.l_dels[TheRec.iDelaunay].EXT;
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: FuncInsertCenter, step 6, B\n";
#endif
      RepresentVertexPermutationPreComput<Tring> precomp_ext(EXT);
      for (auto &ePermGen : ListOrbitDelaunay.l_dels[TheRec.iDelaunay]
                                .GRP.SmallGeneratingSet()) {
        MyMatrix<Tring> eBigMat = precomp_ext.represent(ePermGen);
        MyMatrix<Tring> eBigMat_new = BigMatI * eBigMat * BigMatR;
        FuncInsertGenerator(eBigMat_new);
      }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: FuncInsertCenter, step 6, C\n";
#endif
      std::vector<Tidx> Linc;
      for (int u = 0; u < LVert.rows(); u++) {
        MyVector<Tring> eV = GetMatrixRow(LVert, u);
        size_t pos = ListVertices_rev.at(eV);
        Tidx pos_idx = static_cast<Tidx>(pos - 1);
        Linc.push_back(pos_idx);
      }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: FuncInsertCenter, step 6, D\n";
#endif
      TypeOrbitCenter OrbCent{TheRec.iDelaunay, TheRec.eBigMat, Linc, LVert};
      ListOrbitCenter.push_back(OrbCent);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: FuncInsertCenter, step 6, E\n";
#endif
    }
  };
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FRING: before insert\n";
#endif
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  MicrosecondTime time_c;
#endif
  TypeOrbitCenterMin TheRec{static_cast<int>(eIdx), IdentityMat<Tring>(n + 1)};
  FuncInsertCenter(TheRec);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FRING: after single FuncInsertCenter, before \n";
#endif
  size_t nbCent, nbCentStart = 0;
  while (true) {
    nbCent = ListOrbitCenter.size();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: first while loop nbCentStart=" << nbCentStart
       << " nbCent=" << nbCent << "\n";
#endif
    if (nbCentStart == nbCent) {
      break;
    }
    for (size_t iCent = nbCentStart; iCent < nbCent; iCent++) {
      // Copy is needed since we are inserting so if we used a reference, it
      // might become invalid
      TypeOrbitCenter eEnt = ListOrbitCenter[iCent];
      for (auto &eCase : ListInformationsOneFlipping) {
        if (eEnt.iDelaunay == eCase.iOrb) {
          MyMatrix<Tring> const &eBigMat =
              ListOrbitDelaunay.l_dels[eCase.iOrb].ListAdj[eCase.i_adj].eBigMat;
          int iOrbAdj =
              ListOrbitDelaunay.l_dels[eCase.iOrb].ListAdj[eCase.i_adj].iOrb;
          MyMatrix<Tring> eBigMatNew = eBigMat * eEnt.eBigMat;
          TypeOrbitCenterMin TheRec{iOrbAdj, std::move(eBigMatNew)};
          FuncInsertCenter(TheRec);
        }
      }
    }
    nbCentStart = nbCent;
  }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FRING: after first big loop\n";
#endif
  // Vertex/generator discovery is now complete: PermGRP must reflect all of
  // it before the facet-orbit code below (FuncInsertFacet, Stabilizer_OnSets)
  // starts reading it.
  ensure_group_updated();
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  time_fring_centers += time_c.eval_int64();
  MicrosecondTime time_vm;
#endif
  // second part, the convex decomposition
  int nVert = ListVertices.size();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FRING: nVert=" << nVert << "\n";
#endif
  // The lifted vertex matrix is over the ring: the vertex coordinates are
  // lattice points and the heights are integral since the interior element
  // is handed in fraction-free (a positive rescaling of the lift does not
  // change the regular subdivision).
  MyMatrix<Tring> TotalListVertices(nVert, n + 2);
  MyMatrix<Tring> TotalListVerticesRed(nVert, n + 1);
  std::vector<Tring> LineInterior = GetLineVector(InteriorElement);
  MyVector<Tring> eV(n);
  // The scaled isobarycenter (nVert, sum of the vertex coordinates), see
  // FullRepart.
  MyVector<Tring> eIso = ZeroVector<Tring>(n + 1);
  eIso(0) = nVert;
  for (int iVert = 0; iVert < nVert; iVert++) {
    MyVector<Tring> const &eVert = ListVertices[iVert];
    TotalListVertices(iVert, 0) = 1;
    TotalListVerticesRed(iVert, 0) = 1;
    for (int iCol = 0; iCol < n; iCol++) {
      Tring const &val = eVert(iCol + 1);
      TotalListVerticesRed(iVert, iCol + 1) = val;
      TotalListVertices(iVert, iCol + 1) = val;
      eIso(iCol + 1) += val;
      eV(iCol) = val;
    }
    TotalListVertices(iVert, n + 1) = EvaluateLineVector(LineInterior, eV);
  }
  // TotalListVerticesRed is fixed from here on and FuncInsertFacet represents
  // many permutations of it, so the basis inversion is done once.
  RepresentVertexPermutationPreComput<Tring> precomp_tlvr(TotalListVerticesRed);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  time_fring_vertmat += time_vm.eval_int64();
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FRING: we have TotalListVertices\n";
  os << "ISODEL: FRING: TotalListVertices=\n";
  WriteMatrix(os, TotalListVertices);
  os << "ISODEL: FRING: Before testing symmetry for TotalListVertices\n";
  if (!IsSymmetryGroupOfPolytope(TotalListVertices, PermGRP)) {
    std::cerr << "ISODEL: FRING: PermGRP is not a symmetry group of "
                 "TotalListVertices\n";
    throw TerminalException{1};
  }
  os << "ISODEL: FRING: Before testing symmetry for TotalListVerticesRed_T\n";
  if (!IsSymmetryGroupOfPolytope(TotalListVerticesRed, PermGRP)) {
    std::cerr << "ISODEL: FRING: PermGRP is not a symmetry group of "
                 "TotalListVertices\n";
    throw TerminalException{1};
  }
#endif
  auto get_incd_status = [&](int iVert, MyVector<Tring> const &eFac) -> bool {
    Tring eSum(0);
    for (int u = 0; u <= n + 1; u++) {
      AddMul(eSum, TotalListVertices(iVert, u), eFac(u));
    }
    return eSum == Tring(0);
  };
  // The Linc is not ordered while the Linc_face by virtue of being built as a
  // Face has an ordered intrinsic to it.
  // We need this design for the following reason:
  // * If we reorder the vertices of EXT to match the ListVertices, then the
  // incidence are broken or need to be recomputed.
  // * Converting from one ordering to another is not too expensive overall.
  struct RepartEntryProv {
    MyVector<Tring> eFac;
    std::vector<Tidx> Linc;
    Face Linc_face;
  };
  std::vector<RepartEntry<Tring, Tgroup>> ListOrbitFacet;
  std::vector<RepartEntryProv> ListOrbitFacet_prov;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  MicrosecondTime time_kv;
#endif
  SubsetRankOneSolver<Tring> solver(TotalListVertices);
  for (auto &eRec : ListOrbitCenter) {
    Face Linc_face = VectorToFace(eRec.Linc, nVert);
    MyVector<Tring> eFac = solver.GetPositiveKernelVector(Linc_face);
    Tgroup TheStab;
    int8_t Position = -1;
    std::vector<Delaunay_AdjO<Tring>> ListAdj;
    int iDelaunayOrigin = eRec.iDelaunay;
    MyMatrix<Tring> const &eBigMat = eRec.eBigMat;
    MyMatrix<Tring> const &EXT = eRec.EXT;
    RepartEntry<Tring, Tgroup> re{EXT,     TheStab, Position, iDelaunayOrigin,
                                  ListAdj, eBigMat, {}};
    RepartEntryProv rep{eFac, eRec.Linc, Linc_face};
    ListOrbitFacet.push_back(re);
    ListOrbitFacet_prov.push_back(rep);
  }
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  time_fring_kernel += time_kv.eval_int64();
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FRING: we have ListOrbitFacet\n";
#endif
  auto FuncInsertFacet =
      [&](MyVector<Tring> const &eFac) -> Delaunay_AdjO<Tring> {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertFacet, step 1\n";
#endif
    std::vector<Tidx> Linc;
    Face Linc_face(nVert);
    for (int iVert = 0; iVert < nVert; iVert++) {
      if (get_incd_status(iVert, eFac)) {
        Linc.push_back(iVert);
        Linc_face[iVert] = 1;
      }
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertFacet, step 2 |Linc|=" << Linc.size()
       << "\n";
#endif
    int nOrb = ListOrbitFacet.size();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertFacet, step 2.1 nOrb=" << nOrb << "\n";
#endif
    for (int iOrb = 0; iOrb < nOrb; iOrb++) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: FuncInsertFacet iOrb=" << iOrb
         << " |PermGRP|=" << PermGRP.size() << "\n";
      os << "ISODEL: FRING: FuncInsertFacet Linc_face1="
         << ListOrbitFacet_prov[iOrb].Linc_face.size() << "\n";
      os << "ISODEL: FRING: FuncInsertFacet Linc_face2=" << Linc_face.size()
         << "\n";
#endif
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
      MicrosecondTime time_o;
#endif
      std::optional<Telt> opt = PermGRP.RepresentativeAction_OnSets(
          ListOrbitFacet_prov[iOrb].Linc_face, Linc_face);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
      time_fring_onsets += time_o.eval_int64();
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: FuncInsertFacet, after "
            "RepresentativeAction_OnSets\n";
#endif
      if (opt) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: FRING: FuncInsertFacet, finding an equivalence\n";
#endif
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
        MicrosecondTime time_rv;
#endif
        MyMatrix<Tring> eBigMat = precomp_tlvr.represent(*opt);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
        time_fring_reprvert += time_rv.eval_int64();
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: FRING: FuncInsertFacet, We have eBigMat\n";
#endif
        return {{}, std::move(eBigMat), iOrb};
      }
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertFacet, step 3\n";
#endif
    // A new facet is either barrel or higher because we put all the loiwer ones
    // already
    int8_t Position = 1;
    if (eFac(n + 1) == 0) {
      Position = 0;
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertFacet, step 4\n";
#endif
    Tgroup TheStab;
    int iDelaunayOrigin = -1;
    std::vector<Delaunay_AdjO<Tring>> ListAdj;
    MyMatrix<Tring> eMatUnused; // That matrix should never be used
    MyMatrix<Tring> EXT = SelectRow(TotalListVerticesRed, Linc_face);
    RepartEntry<Tring, Tgroup> re{EXT,     TheStab,    Position, iDelaunayOrigin,
                                  ListAdj, eMatUnused, {}};
    RepartEntryProv rep{eFac, Linc, Linc_face};
    ListOrbitFacet.push_back(re);
    ListOrbitFacet_prov.push_back(rep);
    int iOrb = nOrb;
    MyMatrix<Tring> eBigMat = IdentityMat<Tring>(n + 1);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertFacet, step 5\n";
#endif
    return {{}, std::move(eBigMat), iOrb};
  };
  size_t nOrbStart = 0;
  size_t nOrb;
  while (true) {
    nOrb = ListOrbitFacet.size();
    if (nOrbStart == nOrb) {
      break;
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: second while loop nOrbStart=" << nOrbStart
       << " nOrb=" << nOrb << "\n";
#endif
    for (size_t iOrb = nOrbStart; iOrb < nOrb; iOrb++) {
      Face const &Linc_face = ListOrbitFacet_prov[iOrb].Linc_face;
      std::vector<Tidx> const &Linc = ListOrbitFacet_prov[iOrb].Linc;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
      MicrosecondTime time_s;
#endif
      Tgroup Stab = PermGRP.Stabilizer_OnSets(Linc_face);
      Tgroup TheStabFace = ReducedGroupActionFace(Stab, Linc_face);
      Tgroup TheStabVect = ReducedGroupActionVect(Stab, Linc);
      ListOrbitFacet[iOrb].TheStab = TheStabVect;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
      time_fring_stab += time_s.eval_int64();
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: iOrb=" << iOrb << "\n";
      os << "ISODEL: FRING: |Linc_face|=" << Linc_face.size() << " / "
         << Linc_face.count() << "\n";
      os << "ISODEL: FRING: Position="
         << static_cast<int>(ListOrbitFacet[iOrb].Position) << "\n";
      os << "ISODEL: FRING: Before IsSymmetryGroupOfPolytope after "
            "computation\n";
      if (!IsSymmetryGroupOfPolytope(ListOrbitFacet[iOrb].EXT, TheStabVect)) {
        std::cerr << "ISODEL: FRING: The initial computation went wrong\n";
        throw TerminalException{1};
      }
#endif
      std::vector<Delaunay_AdjO<Tring>> ListAdj;
      // Depending on the nature of the facet (low, barrel, top), the idx_drop
      // can very much vary. We cannot set it to
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
      MicrosecondTime time_e2;
#endif
      int idx_drop = get_idx_drop(ListOrbitFacet_prov[iOrb].eFac);
      MyMatrix<Tring> EXT2 =
          SelectRowDropColumnFace(TotalListVertices, Linc_face, idx_drop);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
      time_fring_ext2 += time_e2.eval_int64();
      MicrosecondTime time_d;
#endif
      vectface vf = DualDescriptionRecordFullDim(EXT2, TheStabFace, rddo_ring);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
      time_fring_dualdesc += time_d.eval_int64();
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: Second while |vf|=" << vf.size()
         << " |TheStabFace|=" << TheStabFace.size() << " |EXT2|=" << EXT2.rows()
         << "/" << EXT2.cols() << " rnk=" << RankMat(EXT2) << "\n";
      CheckFacetInequality(TotalListVertices, Linc_face,
                           "FuncInsertFace TotalListVertices Linc_face");
#endif
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
      MicrosecondTime time_fr;
#endif
      FlippingFramework<Tring> frame(TotalListVertices, Linc_face, os);
      SubsetRankOneSolver<Tring> solver_tlv(TotalListVertices);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
      time_fring_frame += time_fr.eval_int64();
#endif
      for (auto &eFace : vf) {
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN
        CheckFacetInequality(EXT2, eFace, "FuncInsertFace EXT2 eFace");
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: Before FlipFace |EXT2|=" << EXT2.rows() << " / "
           << EXT2.cols() << " |eFace|=" << eFace.size() << " / "
           << eFace.count() << "\n";
#endif
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
        MicrosecondTime time_fl;
#endif
        Face eInc = frame.FlipFace(eFace);
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN
        CheckFacetInequality(TotalListVertices, eInc,
                             "FuncInsertFace TotalListVertices eInc");
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: After FlipFace |eInc|=" << eInc.size() << " / "
           << eInc.count() << "\n";
#endif
        MyVector<Tring> eFac = solver_tlv.GetPositiveKernelVector(eInc);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
        time_fring_flip += time_fl.eval_int64();
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: FRING: We have eFac\n";
#endif
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
        MicrosecondTime time_if;
#endif
        Delaunay_AdjO<Tring> eAdj = FuncInsertFacet(eFac);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
        time_fring_insert_facet += time_if.eval_int64();
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: FRING: We have eAdj\n";
#endif
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
        MicrosecondTime time_lev;
#endif
        size_t nVert = ListOrbitFacet_prov[iOrb].Linc.size();
        Face LEV(nVert);
        for (size_t iInc = 0; iInc < nVert; iInc++) {
          Tidx jInc = ListOrbitFacet_prov[iOrb].Linc[iInc];
          if (get_incd_status(jInc, eFac)) {
            LEV[iInc] = 1;
          }
        }
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
        time_fring_lev += time_lev.eval_int64();
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: FRING: We have |LEV|=" << LEV.size() << " / "
           << LEV.count() << "\n";
#endif
        eAdj.eInc = LEV;
        ListAdj.push_back(eAdj);
      }
      ListOrbitFacet[iOrb].ListAdj = ListAdj;
    }
    nOrbStart = nOrb;
  }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FRING: we have ListOrbitFacet\n";
#endif
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: FRING: group_update|=" << time_fring_group << "\n";
  os << "|ISODEL: FRING: stabilizer_onsets|=" << time_fring_stab << "\n";
  os << "|ISODEL: FRING: dualdesc|=" << time_fring_dualdesc << "\n";
  os << "|ISODEL: FRING: flip_kernel|=" << time_fring_flip << "\n";
  os << "|ISODEL: FRING: insert_facet|=" << time_fring_insert_facet << "\n";
  os << "|ISODEL: FRING: insert_facet_onsets|=" << time_fring_onsets << "\n";
  os << "|ISODEL: FRING: insert_facet_reprvert|=" << time_fring_reprvert
     << "\n";
  os << "|ISODEL: FRING: centers|=" << time_fring_centers << "\n";
  os << "|ISODEL: FRING: vertex_matrix|=" << time_fring_vertmat << "\n";
  os << "|ISODEL: FRING: kernel_vectors|=" << time_fring_kernel << "\n";
  os << "|ISODEL: FRING: select_ext2|=" << time_fring_ext2 << "\n";
  os << "|ISODEL: FRING: frame_setup|=" << time_fring_frame << "\n";
  os << "|ISODEL: FRING: lev|=" << time_fring_lev << "\n";
  os << "|ISODEL: FRING: full|=" << time_fring_full << "\n";
#endif
  return {ListOrbitFacet, eIso};
}

/*
  Effectively do the flipping. Note that the whole computation is combinatorial
  and all the polyhedral computations are done in the
  FindRepartitionningInfoNextGeneration code.
 */
template <typename Tring, typename Tgroup>
DelaunayTesselationIneq<Tring, Tgroup>
FlippingLtype(
    DelaunayTesselationIneq<Tring, Tgroup> const &ListOrbitDelaunay,
    MyMatrix<Tring> const &InteriorElement,
    std::vector<AdjInfo> const &ListInformationsOneFlipping,
    std::vector<std::vector<Tring>> const &ListGramRing,
    PolyHeuristicSerial<typename Tgroup::Tint> &AllArr, std::ostream &os) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FLT: FlippingLtype, begin\n";
#endif
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  MicrosecondTime time_flt_full;
  int64_t time_flt_repart = 0, time_flt_onsets = 0, time_flt_bigmat = 0,
          time_flt_vipc = 0, time_flt_ineq = 0, time_flt_setup = 0,
          time_flt_symbols = 0, time_flt_adjacencies = 0,
          time_flt_matching = 0, time_flt_matching_old = 0,
          time_flt_transport = 0, time_flt_symbolpos = 0,
          time_flt_chain = 0;
  MicrosecondTime time_su;
#endif
  using Tgr = GraphListAdj;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n_dels = ListOrbitDelaunay.l_dels.size();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FLT: n_dels=" << n_dels << "\n";
#endif
  Tgr Gra(n_dels);
  Face ListMatched(n_dels);
  for (auto &eAI : ListInformationsOneFlipping) {
    ListMatched[eAI.iOrb] = 1;
    int iOrbAdj = ListOrbitDelaunay.l_dels[eAI.iOrb].ListAdj[eAI.i_adj].iOrb;
    if (eAI.iOrb != iOrbAdj) {
      Gra.AddAdjacent(eAI.iOrb, iOrbAdj);
    }
  }
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN
  if (!IsSymmetricGraph(Gra)) {
    std::cerr << "ISODEL: The graph is not symmetric\n";
    throw TerminalException{1};
  }
#endif
  std::vector<std::vector<size_t>> ListConn = ConnectedComponents_set(Gra);
  std::vector<std::vector<size_t>> ListGroupMelt, ListGroupUnMelt;
  auto is_melt = [&](std::vector<size_t> const &eConn) -> bool {
    size_t n_matched = 0;
    for (auto &eVal : eConn) {
      n_matched += ListMatched[eVal];
    }
    if (n_matched == eConn.size()) {
      return true;
    }
    if (n_matched == 0) {
      return false;
    }
    std::cerr << "ISODEL: FLT: The melt should be total or none at all\n";
    throw TerminalException{1};
  };
  for (auto &eConn : ListConn) {
    if (is_melt(eConn)) {
      ListGroupMelt.push_back(eConn);
    } else {
      ListGroupUnMelt.push_back(eConn);
    }
  }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FLT: |ListGroupMelt|=" << ListGroupMelt.size()
     << " |ListGroupUnMelt|=" << ListGroupUnMelt.size() << "\n";
#endif
  std::vector<std::vector<RepartEntry<Tring, Tgroup>>> ListInfo;
  std::vector<MyVector<Tring>> ListIso;
  std::vector<int> vect_iInfo(n_dels, -1);
  std::vector<int> vect_lower_iFacet(n_dels, -1);
  int iInfo = 0;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  time_flt_setup += time_su.eval_int64();
#endif
  for (auto &eConn : ListGroupMelt) {
    size_t eIdx = eConn[0];
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
    MicrosecondTime time_r;
#endif
    FullRepart<Tring, Tgroup> fr = FindRepartitionningInfoNextGeneration(
        eIdx, ListOrbitDelaunay, ListInformationsOneFlipping, InteriorElement,
        AllArr, os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
    time_flt_repart += time_r.eval_int64();
#endif
    int n_facet = fr.cells.size();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FLT: iInfo=" << iInfo << " |eConn|=" << eConn.size()
       << " n_facet=" << n_facet << "\n";
    size_t n_zero = 0, n_plus = 0, n_minus = 0;
#endif
    for (int iFacet = 0; iFacet < n_facet; iFacet++) {
      RepartEntry<Tring, Tgroup> const &eFacet = fr.cells[iFacet];
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FLT: iFacet=" << iFacet
         << " Position=" << static_cast<int>(eFacet.Position) << "\n";
      if (eFacet.Position == -1) {
        n_minus += 1;
      }
      if (eFacet.Position == 0) {
        n_zero += 1;
      }
      if (eFacet.Position == 1) {
        n_plus += 1;
      }
#endif
      if (eFacet.Position == -1) {
        vect_lower_iFacet[eFacet.iDelaunayOrigin] = iFacet;
      }
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FLT: n_minus=" << n_minus << " n_zero=" << n_zero
       << " n_plus=" << n_plus << "\n";
#endif
    ListInfo.push_back(fr.cells);
    ListIso.push_back(fr.eIso);
    for (auto &pos : eConn) {
      vect_iInfo[pos] = iInfo;
    }
    iInfo++;
  }
  struct MatchedFacet {
    Delaunay_AdjO<Tring> adj;
    MyMatrix<Tring> eBigMat;
  };
  // Lazy per-orbit basis-inversion caches for representing vertex
  // permutations: entry iDelaunay for the old tessellation orbits, entry
  // (iInfo, iFacet) for the repartitioning facets. A barrel facet caches
  // the precompute of its isobarycenter-extended matrix (see get_bigmat),
  // the other facets that of their plain vertex matrix; each slot is only
  // ever used by the branch matching its Position.
  std::vector<std::optional<RepresentVertexPermutationPreComput<Tring>>>
      l_precomp_old(n_dels);
  std::vector<
      std::vector<std::optional<RepresentVertexPermutationPreComput<Tring>>>>
      ll_precomp_info;
  for (auto &eInfo : ListInfo) {
    ll_precomp_info.emplace_back(eInfo.size());
  }
  auto get_matching_listinfo = [&](int const &iInfo, int const &iFacet,
                                   Face const &eInc) -> MatchedFacet {
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
    MicrosecondTime time_gm;
    struct AccumOnExit {
      MicrosecondTime &t;
      int64_t &acc;
      ~AccumOnExit() { acc += t.eval_int64(); }
    } accum_gm{time_gm, time_flt_matching};
#endif
    MyMatrix<Tring> const &EXT = ListInfo[iInfo][iFacet].EXT;
    Tgroup const &TheStab = ListInfo[iInfo][iFacet].TheStab;
    int dim = InteriorElement.rows();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FLT: iInfo=" << iInfo << " iFacet=" << iFacet
       << " |eInc|=" << eInc.size() << " / " << eInc.count() << "\n";
    os << "ISODEL: FLT: Position="
       << static_cast<int>(ListInfo[iInfo][iFacet].Position) << "\n";
    os << "ISODEL: FLT: |EXT|=" << EXT.rows() << " / " << EXT.cols()
       << " rnk=" << RankMat(EXT) << "\n";
    os << "ISODEL: FLT: |TheStab|=" << TheStab.size() << "\n";
    os << "ISODEL: FLT: |ListAdj|=" << ListInfo[iInfo][iFacet].ListAdj.size()
       << "\n";
#endif
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN
    if (RankMat(EXT) == dim + 1) {
      CheckFacetInequality(EXT, eInc, "get_matching_listinfo EXT eInc");
    }
    if (!IsSymmetryGroupOfPolytope(EXT, TheStab)) {
      std::cerr << "The group TheStab is not a symmetry group\n";
      throw TerminalException{1};
    }
#endif
    auto get_bigmat = [&](Telt const &ePerm) -> MyMatrix<Tring> {
      int8_t Position = ListInfo[iInfo][iFacet].Position;
      std::optional<RepresentVertexPermutationPreComput<Tring>> &precomp =
          ll_precomp_info[iInfo][iFacet];
      if (Position == 0) {
        // For barrel case, we need to extend by one point in order to get
        // a full dimensional set. That point has to be preserved by the
        // corresponding group and so we take the (scaled, hence integral)
        // isobarycenter of the repartitioning polytope: the transformation
        // fixes it whatever the scaling, since a row and its image scale
        // alike in the solve.
        int nbRow = EXT.rows();
        if (!precomp) {
          MyMatrix<Tring> EXT_ext(nbRow + 1, dim + 1);
          for (int iRow = 0; iRow < nbRow; iRow++) {
            for (int i = 0; i <= dim; i++) {
              EXT_ext(iRow, i) = EXT(iRow, i);
            }
          }
          for (int i = 0; i <= dim; i++) {
            EXT_ext(nbRow, i) = ListIso[iInfo](i);
          }
          precomp.emplace(std::move(EXT_ext));
        }
        std::vector<Tidx> eList(nbRow + 1);
        for (int u = 0; u < nbRow; u++) {
          eList[u] = ePerm.at(u);
        }
        eList[nbRow] = nbRow;
        Telt ePermExt(std::move(eList));
        return precomp->represent(ePermExt);
      } else {
        if (!precomp) {
          precomp.emplace(EXT);
        }
        return precomp->represent(ePerm);
      }
    };
    for (auto &eAdj : ListInfo[iInfo][iFacet].ListAdj) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FLT: |eAdj.eInc|=" << eAdj.eInc.size() << " / "
         << eAdj.eInc.count() << " |eInc|=" << eInc.size() << " / "
         << eInc.count() << "\n";
      os << "ISODEL: FLT: TheStab, n_act=" << TheStab.n_act()
         << " order=" << TheStab.size() << "\n";
#endif
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN
      if (RankMat(EXT) == dim + 1) {
        CheckFacetInequality(EXT, eAdj.eInc,
                             "get_matching_listinfo EXT eAdj.eInc");
      }
#endif
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
      MicrosecondTime time_o;
#endif
      std::optional<Telt> opt =
          TheStab.RepresentativeAction_OnSets(eAdj.eInc, eInc);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
      time_flt_onsets += time_o.eval_int64();
#endif
      if (opt) {
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
        MicrosecondTime time_b;
#endif
        MyMatrix<Tring> eBigMat = get_bigmat(*opt);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
        time_flt_bigmat += time_b.eval_int64();
#endif
        return {eAdj, eBigMat};
      }
    }
    std::cerr << "ISO_EDGE: FLT: iInfo=" << iInfo << " iFacet=" << iFacet
              << " |eInc|=" << eInc.size() << " / " << eInc.count() << "\n";
    std::cerr << "ISO_EDGE: FLT: Failed to find a matching entry in "
                 "get_matching_listinfo\n";
    throw TerminalException{1};
  };
  auto get_matching_old_tessel = [&](int const &iDelaunayOld,
                                     Face const &eInc) -> MatchedFacet {
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
    MicrosecondTime time_go;
    struct AccumOnExit {
      MicrosecondTime &t;
      int64_t &acc;
      ~AccumOnExit() { acc += t.eval_int64(); }
    } accum_go{time_go, time_flt_matching_old};
#endif
    MyMatrix<Tring> const &EXT = ListOrbitDelaunay.l_dels[iDelaunayOld].EXT;
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN
    CheckFacetInequality(EXT, eInc, "get_matching_old_tessel EXT eInc");
#endif
    for (auto &eAdj : ListOrbitDelaunay.l_dels[iDelaunayOld].ListAdj) {
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
      MicrosecondTime time_o;
#endif
      std::optional<Telt> opt =
          ListOrbitDelaunay.l_dels[iDelaunayOld]
              .GRP.RepresentativeAction_OnSets(eAdj.eInc, eInc);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
      time_flt_onsets += time_o.eval_int64();
#endif
      if (opt) {
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
        MicrosecondTime time_b;
#endif
        std::optional<RepresentVertexPermutationPreComput<Tring>> &precomp =
            l_precomp_old[iDelaunayOld];
        if (!precomp) {
          precomp.emplace(EXT);
        }
        MyMatrix<Tring> eBigMat = precomp->represent(*opt);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
        time_flt_bigmat += time_b.eval_int64();
#endif
        return {{eAdj.eInc, eAdj.eBigMat, eAdj.iOrb}, eBigMat};
      }
    }
    std::cerr << "ISO_EDGE: FLT: iDelaunayOld=" << iDelaunayOld
              << " |eInc|=" << eInc.size() << " / " << eInc.count() << "\n";
    std::cerr << "ISO_EDGE: FLT: Failed to find a matching entry in "
                 "get_matching_old_tessel\n";
    throw TerminalException{1};
  };
  auto get_lower_adjacency = [&](int const &iInfo,
                                 int const &iFacet) -> Delaunay_AdjO<Tring> {
    for (auto &fAdj : ListInfo[iInfo][iFacet].ListAdj) {
      if (ListInfo[iInfo][fAdj.iOrb].Position == -1) {
        return fAdj;
      }
    }
    std::cerr << "ISO_EDGE: FLT: iInfo=" << iInfo << " iFacet=" << iFacet
              << "\n";
    std::cerr
        << "ISO_EDGE: FLT: Failed to find an entry in get_lower_adjacency\n";
    throw TerminalException{1};
  };
  auto get_face_m_m = [&](MyMatrix<Tring> const &M1,
                          MyMatrix<Tring> const &M2) -> Face {
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
    MicrosecondTime time_tr;
#endif
    int nVert1 = M1.rows();
    int nVert2 = M2.rows();
    Face f2(nVert2);
    ContainerMatrix<Tring> cont(M2);
    for (int i1 = 0; i1 < nVert1; i1++) {
      MyVector<Tring> V = GetMatrixRow(M1, i1);
      std::optional<size_t> opt = cont.GetIdx_v(V);
      size_t pos2 = unfold_opt(opt, "Error in get_face_m_m");
      f2[pos2] = 1;
    }
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
    time_flt_transport += time_tr.eval_int64();
#endif
    return f2;
  };
  auto get_face_msub_m = [&](MyMatrix<Tring> const &M1, Face const &f1,
                             MyMatrix<Tring> const &M2) -> Face {
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
    MicrosecondTime time_tr;
#endif
    int nVert1 = M1.rows();
    int nVert2 = M2.rows();
    Face f2(nVert2);
    ContainerMatrix<Tring> cont(M2);
    for (int i1 = 0; i1 < nVert1; i1++) {
      if (f1[i1] == 1) {
        MyVector<Tring> V = GetMatrixRow(M1, i1);
        std::optional<size_t> opt = cont.GetIdx_v(V);
        int pos2 = unfold_opt(opt, "Error in get_face_msub_m");
        f2[pos2] = 1;
      }
    }
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
    time_flt_transport += time_tr.eval_int64();
#endif
    return f2;
  };
  int n_info = ListInfo.size();
  int8_t Position_old = 4, Position_new = 5;
  struct DelaunaySymb {
    int8_t Position;
    int iDelaunay;
    int iInfo;
    int iFacet;
  };
  std::vector<DelaunaySymb> NewListOrbitDelaunay;
  auto get_symbol_position =
      [&](DelaunaySymb const &ds) -> std::optional<size_t> {
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
    MicrosecondTime time_sp;
    struct AccumOnExit {
      MicrosecondTime &t;
      int64_t &acc;
      ~AccumOnExit() { acc += t.eval_int64(); }
    } accum_sp{time_sp, time_flt_symbolpos};
#endif
    if (ds.Position == Position_old) {
      for (size_t i = 0; i < NewListOrbitDelaunay.size(); i++) {
        if (NewListOrbitDelaunay[i].Position == Position_old) {
          if (NewListOrbitDelaunay[i].iDelaunay == ds.iDelaunay) {
            return i;
          }
        }
      }
    }
    if (ds.Position == Position_new) {
      for (size_t i = 0; i < NewListOrbitDelaunay.size(); i++) {
        if (NewListOrbitDelaunay[i].Position == Position_new) {
          if (NewListOrbitDelaunay[i].iInfo == ds.iInfo &&
              NewListOrbitDelaunay[i].iFacet == ds.iFacet) {
            return i;
          }
        }
      }
    }
    return {};
  };
  std::vector<Delaunay_EntryIneq<Tring, Tgroup>> l_dels;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  MicrosecondTime time_sy;
#endif
  for (auto &eConn : ListGroupUnMelt) {
    if (eConn.size() > 1) {
      std::cerr << "Error of connected component computation\n";
      throw TerminalException{1};
    }
    int iDelaunay = eConn[0];
    DelaunaySymb ds{Position_old, iDelaunay, -1, -1};
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: |l_dels|=" << l_dels.size()
       << " Position_old=" << static_cast<int>(Position_old)
       << " iDelaunay=" << iDelaunay << "\n";
#endif
    NewListOrbitDelaunay.push_back(ds);
    Delaunay_EntryIneq<Tring, Tgroup> del{
      ListOrbitDelaunay.l_dels[iDelaunay].EXT,
      ListOrbitDelaunay.l_dels[iDelaunay].GRP,
      {}};
    l_dels.push_back(del);
  }
  for (int iInfo = 0; iInfo < n_info; iInfo++) {
    int n_facet = ListInfo[iInfo].size();
    for (int iFacet = 0; iFacet < n_facet; iFacet++) {
      int8_t Position = ListInfo[iInfo][iFacet].Position;
      if (Position == 1) {
        DelaunaySymb ds{Position_new, -1, iInfo, iFacet};
        NewListOrbitDelaunay.push_back(ds);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: |l_dels|=" << l_dels.size() << " iInfo=" << iInfo
           << " iFacet=" << iFacet << "\n";
#endif
        Tgroup GRP = ListInfo[iInfo][iFacet].TheStab;
        Delaunay_EntryIneq<Tring, Tgroup> del{ListInfo[iInfo][iFacet].EXT, GRP,
                                          {}};
        l_dels.push_back(del);
      }
    }
  }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FLT: |l_dels|=" << l_dels.size() << "\n";
  for (size_t i_del = 0; i_del < l_dels.size(); i_del++) {
    DelaunaySymb ds = NewListOrbitDelaunay[i_del];
    int Position = ds.Position;
    int iDelaunay = ds.iDelaunay;
    int iInfo = ds.iInfo;
    int iFacet = ds.iFacet;
    os << "ISODEL: FLT: i_del=" << i_del << " Position=" << Position
       << " iDelaunay=" << iDelaunay << " iInfo=" << iInfo
       << " iFacet=" << iFacet << "\n";
  }
  auto check_adj = [&](int iOrb, Delaunay_AdjIneqO<Tring> const &NewAdj,
                       std::string const &context) -> void {
    MyMatrix<Tring> const &EXT = l_dels[iOrb].EXT;
    ContainerMatrix<Tring> cont(EXT);
    Face f_att(EXT.rows());
    MyMatrix<Tring> EXTadj = l_dels[NewAdj.iOrb].EXT * NewAdj.eBigMat;
    os << "ISODEL: FLT: check_adj iOrb=" << iOrb
       << " NewAdj.iOrb=" << NewAdj.iOrb << "\n";
    os << "ISODEL: FLT:      |EXT|=" << EXT.rows() << " / " << EXT.cols()
       << "\n";
    os << "ISODEL: FLT:   |EXTadj|=" << EXTadj.rows() << " / " << EXTadj.cols()
       << "\n";
    int len = EXTadj.rows();
    for (int iVert = 0; iVert < len; iVert++) {
      MyVector<Tring> V = GetMatrixRow(EXTadj, iVert);
      std::optional<size_t> opt = cont.GetIdx_v(V);
      if (opt) {
        size_t pos = *opt;
        os << "ISODEL: FLT: iVert=" << iVert << " pos=" << pos << "\n";
        f_att[*opt] = 1;
      }
    }
    os << "ISODEL: FLT:     |f_att|=" << f_att.size() << " / " << f_att.count()
       << "\n";
    os << "ISODEL: FLT: NewAdj.eInc=" << NewAdj.eInc.size() << " / "
       << NewAdj.eInc.count() << "\n";
    if (f_att != NewAdj.eInc) {
      std::cerr << "f_att should be equal to NewAdj.eInc\n";
      std::cerr << "Consistency error in context=" << context << "\n";
      throw TerminalException{1};
    }
  };
#endif
  int n_del_ret = l_dels.size();
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  time_flt_symbols += time_sy.eval_int64();
  MicrosecondTime time_ad;
#endif
  for (int iOrb = 0; iOrb < n_del_ret; iOrb++) {
    DelaunaySymb ds = NewListOrbitDelaunay[iOrb];
    std::vector<Delaunay_AdjIneqO<Tring>> ListAdj;
    // Most orbits surviving a flip have every adjacency in Case 1 (no
    // recomputation needed at all), so building the Voronoi precompute for
    // this orbit is deferred until an adjacency actually requires it.
    std::optional<VoronoiInequalityPreComput<Tring>> vipc_opt;
    std::optional<ContainerMatrix<Tring>> cont_opt;
    auto ensure_vipc_cont = [&]() -> void {
      if (!vipc_opt) {
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
        MicrosecondTime time_v;
#endif
        vipc_opt.emplace(BuildVoronoiIneqPreComputeChecked<Tring>(
            l_dels[iOrb].EXT, ListGramRing, os));
        cont_opt.emplace(l_dels[iOrb].EXT);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
        time_flt_vipc += time_v.eval_int64();
#endif
      }
    };
    if (ds.Position == Position_old) {
      int iDelaunay = ds.iDelaunay;
      int n_adj_old = ListOrbitDelaunay.l_dels[iDelaunay].ListAdj.size();
      for (int i_adj_old = 0; i_adj_old < n_adj_old; i_adj_old++) {
        Delaunay_AdjIneqO<Tring> const &eAdj =
            ListOrbitDelaunay.l_dels[iDelaunay].ListAdj[i_adj_old];
        int iDelaunayOld = eAdj.iOrb;
        DelaunaySymb dss{Position_old, iDelaunayOld, -1, -1};
        std::optional<size_t> opt = get_symbol_position(dss);
        if (opt) {
          int iOrbAdj = *opt;
          // Case 1: neither this orbit nor its neighbour were touched by the
          // flip (both stay Position_old), so the adjacency (eInc, eBigMat)
          // is byte-identical to the old one and so is its inequality.
          MyVector<Tring> eIneq = eAdj.eIneq;
#ifdef SANITY_CHECK_FLIP_CASE1_INEQ
          {
            ensure_vipc_cont();
            MyVector<Tring> eIneqCheck = ComputeDelaunayAdjIneq(
                *vipc_opt, *cont_opt, l_dels[iOrbAdj].EXT, eAdj.eBigMat,
                ListGramRing, os);
            if (eIneqCheck != eIneq) {
              std::cerr << "FLIP_CASE1_INEQ: SANITY_CHECK failed: the "
                           "inequality of an unmelted (Case 1) adjacency "
                           "changed after the flip, but it should be "
                           "identical to the old one\n";
              throw TerminalException{1};
            }
          }
#endif
          Delaunay_AdjIneqO<Tring> NAdj{eAdj.eInc, eAdj.eBigMat, iOrbAdj,
                                        eIneq};
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
          os << "ISODEL: eAdj.eBigMat=\n";
          WriteMatrix(os, eAdj.eBigMat);
          os << "ISODEL: iOrb=" << iOrb << " iOrbAdj=" << iOrbAdj << "\n";
          check_adj(iOrb, NAdj, "Case 1");
#endif
          ListAdj.push_back(NAdj);
        } else {
          int iInfo = vect_iInfo[iDelaunayOld];
          int iFacet = vect_lower_iFacet[iDelaunayOld];
          RepartEntry<Tring, Tgroup> const &eFacet = ListInfo[iInfo][iFacet];
          MyMatrix<Tring> ImageEXT =
              ListOrbitDelaunay.l_dels[iDelaunayOld].EXT * eAdj.eBigMat;
          Face Linc = get_face_msub_m(ListOrbitDelaunay.l_dels[iDelaunay].EXT,
                                      eAdj.eInc, ImageEXT);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
          os << "ISODEL: FLT: iInfo=" << iInfo << " iFacet=" << iFacet
             << " |Linc|=" << Linc.size() << " / " << Linc.count() << "\n";
          os << "ISODEL: FLT: Case 2\n";
#endif
          MatchedFacet RecMatch = get_matching_listinfo(iInfo, iFacet, Linc);
          int iOrbFound = RecMatch.adj.iOrb;
          if (ListInfo[iInfo][iFacet].Position == 0) {
            std::cerr << "Illogic error concerning the structure of "
                         "repartitionning polytope\n";
            throw TerminalException{1};
          }
          DelaunaySymb dss{Position_new, -1, iInfo, iOrbFound};
          std::optional<int> optN = get_symbol_position(dss);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
          os << "Position_new=" << static_cast<int>(Position_new)
             << " iInfo=" << iInfo << " iOrbFound=" << iOrbFound << "\n";
#endif
          int Pos = unfold_opt(optN, "Failed to find entry for Case 2");
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
          MicrosecondTime time_ch1;
#endif
          MyMatrix<Tring> BigMat1 = RecMatch.adj.eBigMat * RecMatch.eBigMat *
                                    eFacet.get_bigmat_inv() * eAdj.eBigMat;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
          time_flt_chain += time_ch1.eval_int64();
#endif
          // Case 2: this orbit is unmelted but its neighbour was melted into
          // a new polytope, so the shared facet's inequality must be
          // recomputed against the new neighbour.
          ensure_vipc_cont();
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
          MicrosecondTime time_i;
#endif
          MyVector<Tring> eIneq = ComputeDelaunayAdjIneq(
              *vipc_opt, *cont_opt, l_dels[Pos].EXT, BigMat1, ListGramRing,
              os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
          time_flt_ineq += time_i.eval_int64();
#endif
          Delaunay_AdjIneqO<Tring> NAdj{eAdj.eInc, std::move(BigMat1), Pos,
                                    std::move(eIneq)};
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
          check_adj(iOrb, NAdj, "Case 2");
#endif
          ListAdj.push_back(NAdj);
        }
      }
    } else {
      int iInfo = ds.iInfo;
      int iFacet = ds.iFacet;
      for (auto &eAdj : ListInfo[iInfo][iFacet].ListAdj) {
        int jFacet = eAdj.iOrb;
        if (ListInfo[iInfo][jFacet].Position == 0) {
          MyMatrix<Tring> const &eMat1 = eAdj.eBigMat;
          MyMatrix<Tring> LincEXT =
              SelectRow(ListInfo[iInfo][iFacet].EXT, eAdj.eInc);
          MyMatrix<Tring> ImageEXTbarrel = ListInfo[iInfo][jFacet].EXT * eMat1;
          Face LLinc = get_face_m_m(LincEXT, ImageEXTbarrel);
          Delaunay_AdjO<Tring> TheFoundAdj = get_lower_adjacency(iInfo, jFacet);
          int kFacet = TheFoundAdj.iOrb;
          MyMatrix<Tring> const &eMat2 = TheFoundAdj.eBigMat;
          int iDelaunayOrigin = ListInfo[iInfo][kFacet].iDelaunayOrigin;
          MyMatrix<Tring> ImageEXT =
              ListInfo[iInfo][kFacet].EXT * TheFoundAdj.eBigMat;
          Face LLinc2 = get_face_msub_m(ListInfo[iInfo][jFacet].EXT,
                                        TheFoundAdj.eInc, ImageEXT);
          MatchedFacet match2 =
              get_matching_old_tessel(iDelaunayOrigin, LLinc2);
          Delaunay_AdjO<Tring> const &TheFoundAdj2 = match2.adj;
          MyMatrix<Tring> const &TheMat2 = match2.eBigMat;
          int jDelaunayOrigin = TheFoundAdj2.iOrb;
          MyMatrix<Tring> const &eMat3 = TheFoundAdj2.eBigMat;
          ImageEXT = ListOrbitDelaunay.l_dels[jDelaunayOrigin].EXT * eMat3;
          Face LLinc3 =
              get_face_msub_m(ListOrbitDelaunay.l_dels[iDelaunayOrigin].EXT,
                              TheFoundAdj2.eInc, ImageEXT);
          int jInfo = vect_iInfo[jDelaunayOrigin];
          int iFacet2 = vect_lower_iFacet[jDelaunayOrigin];
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
          os << "ISODEL: FLT: Case 3, step 1\n";
#endif
          MatchedFacet match3 = get_matching_listinfo(jInfo, iFacet2, LLinc3);
          Delaunay_AdjO<Tring> const &TheFoundAdj3 = match3.adj;
          MyMatrix<Tring> const &TheMat3 = match3.eBigMat;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
          MicrosecondTime time_ch2;
#endif
          MyMatrix<Tring> BigMat1 =
              TheFoundAdj3.eBigMat * TheMat3 *
              ListInfo[jInfo][iFacet2].get_bigmat_inv() * eMat3 * TheMat2 *
              ListInfo[iInfo][kFacet].eBigMat * eMat2 * eMat1;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
          time_flt_chain += time_ch2.eval_int64();
#endif
          MyMatrix<Tring> EXT7 =
              ListInfo[jInfo][TheFoundAdj3.iOrb].EXT * BigMat1;
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN
          if (SortMatrix(EXT7) != SortMatrix(ImageEXTbarrel)) {
            std::cerr << "We fail an important test with the barrel images\n";
            throw TerminalException{1};
          }
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
          os << "ISODEL: FLT: |EXT7|=" << EXT7.rows() << " / " << EXT7.cols()
             << " rnk=" << RankMat(EXT7) << "\n";
          os << "ISODEL: FLT: |LincEXT|=" << LincEXT.rows() << " / "
             << LincEXT.cols() << " rnk=" << RankMat(LincEXT) << "\n";
#endif
          Face LLinc4 = get_face_m_m(LincEXT, EXT7);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
          os << "ISODEL: FLT: Case 3, step 2\n";
#endif
          MatchedFacet match4 =
              get_matching_listinfo(jInfo, TheFoundAdj3.iOrb, LLinc4);
          Delaunay_AdjO<Tring> const &TheFoundAdj4 = match4.adj;
          MyMatrix<Tring> const &TheMat4 = match4.eBigMat;
          DelaunaySymb dss{Position_new, -1, jInfo, TheFoundAdj4.iOrb};
          std::optional<size_t> optN = get_symbol_position(dss);
          int Pos = unfold_opt(optN, "Failed to find entry for Case 3");
          MyMatrix<Tring> BigMat2 = TheFoundAdj4.eBigMat * TheMat4 * BigMat1;
          ensure_vipc_cont();
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
          MicrosecondTime time_i;
#endif
          MyVector<Tring> eIneq = ComputeDelaunayAdjIneq(
              *vipc_opt, *cont_opt, l_dels[Pos].EXT, BigMat2, ListGramRing,
              os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
          time_flt_ineq += time_i.eval_int64();
#endif
          Delaunay_AdjIneqO<Tring> NAdj{eAdj.eInc, std::move(BigMat2), Pos,
                                    std::move(eIneq)};
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
          check_adj(iOrb, NAdj, "Case 3");
#endif
          ListAdj.push_back(NAdj);
        }
        if (ListInfo[iInfo][jFacet].Position == -1) {
          int iDelaunayOrigin = ListInfo[iInfo][jFacet].iDelaunayOrigin;
          MyMatrix<Tring> const &eMat1 = eAdj.eBigMat;
          MyMatrix<Tring> ImageEXT = ListInfo[iInfo][jFacet].EXT * eMat1;
          Face LLinc =
              get_face_msub_m(ListInfo[iInfo][iFacet].EXT, eAdj.eInc, ImageEXT);
          MatchedFacet match1 = get_matching_old_tessel(iDelaunayOrigin, LLinc);
          Delaunay_AdjO<Tring> const &TheFoundAdj1 = match1.adj;
          MyMatrix<Tring> const &TheMat1 = match1.eBigMat;
          int jDelaunayOld = TheFoundAdj1.iOrb;
          if (vect_iInfo[jDelaunayOld] == -1) {
            DelaunaySymb dss{Position_old, jDelaunayOld, -1, -1};
            std::optional<size_t> opt = get_symbol_position(dss);
            int Pos2 = unfold_opt(opt, "Case 4");
            MyMatrix<Tring> BigMat1 = TheFoundAdj1.eBigMat * TheMat1 *
                                      ListInfo[iInfo][jFacet].eBigMat *
                                      eAdj.eBigMat;
            ensure_vipc_cont();
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
            MicrosecondTime time_i;
#endif
            MyVector<Tring> eIneq = ComputeDelaunayAdjIneq(
                *vipc_opt, *cont_opt, l_dels[Pos2].EXT, BigMat1, ListGramRing,
                os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
            time_flt_ineq += time_i.eval_int64();
#endif
            Delaunay_AdjIneqO<Tring> NAdj{eAdj.eInc, std::move(BigMat1), Pos2,
                                      std::move(eIneq)};
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
            check_adj(iOrb, NAdj, "Case 4");
#endif
            ListAdj.push_back(NAdj);
          } else {
            int jInfo = vect_iInfo[jDelaunayOld];
            int iFacet2 = vect_lower_iFacet[jDelaunayOld];
            MyMatrix<Tring> ImageEXT =
                ListOrbitDelaunay.l_dels[jDelaunayOld].EXT *
                TheFoundAdj1.eBigMat;
            Face LLinc2 =
                get_face_msub_m(ListOrbitDelaunay.l_dels[iDelaunayOrigin].EXT,
                                TheFoundAdj1.eInc, ImageEXT);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
            os << "ISODEL: FLT: Case 5\n";
#endif
            MatchedFacet match2 = get_matching_listinfo(jInfo, iFacet2, LLinc2);
            Delaunay_AdjO<Tring> const &TheFoundAdj2 = match2.adj;
            MyMatrix<Tring> const &TheMat2 = match2.eBigMat;
            DelaunaySymb dss{Position_new, -1, jInfo, TheFoundAdj2.iOrb};
            std::optional<size_t> opt = get_symbol_position(dss);
            int Pos = unfold_opt(opt, "Case 5");
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
            MicrosecondTime time_ch3;
#endif
            MyMatrix<Tring> BigMat1 =
                TheFoundAdj2.eBigMat * TheMat2 *
                ListInfo[jInfo][iFacet2].get_bigmat_inv() *
                TheFoundAdj1.eBigMat * TheMat1 *
                ListInfo[iInfo][jFacet].eBigMat * eAdj.eBigMat;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
            time_flt_chain += time_ch3.eval_int64();
#endif
            ensure_vipc_cont();
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
            MicrosecondTime time_i;
#endif
            MyVector<Tring> eIneq = ComputeDelaunayAdjIneq(
                *vipc_opt, *cont_opt, l_dels[Pos].EXT, BigMat1, ListGramRing,
                os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
            time_flt_ineq += time_i.eval_int64();
#endif
            Delaunay_AdjIneqO<Tring> NAdj{eAdj.eInc, std::move(BigMat1), Pos,
                                      std::move(eIneq)};
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
            check_adj(iOrb, NAdj, "Case 5");
#endif
            ListAdj.push_back(NAdj);
          }
        }
        if (ListInfo[iInfo][jFacet].Position == 1) {
          DelaunaySymb dss{Position_new, -1, iInfo, jFacet};
          std::optional<size_t> opt = get_symbol_position(dss);
          int Pos = unfold_opt(opt, "Case 5");
          ensure_vipc_cont();
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
          MicrosecondTime time_i;
#endif
          MyVector<Tring> eIneq = ComputeDelaunayAdjIneq(
              *vipc_opt, *cont_opt, l_dels[Pos].EXT, eAdj.eBigMat,
              ListGramRing, os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
          time_flt_ineq += time_i.eval_int64();
#endif
          Delaunay_AdjIneqO<Tring> NAdj{eAdj.eInc, eAdj.eBigMat, Pos,
                                    std::move(eIneq)};
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
          check_adj(iOrb, NAdj, "Case 6");
#endif
          ListAdj.push_back(NAdj);
        }
      }
    }
    l_dels[iOrb].ListAdj = ListAdj;
  }
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  time_flt_adjacencies += time_ad.eval_int64();
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FLT: Exiting\n";
#endif
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: FLT: repartitionning|=" << time_flt_repart << "\n";
  os << "|ISODEL: FLT: representative_onsets|=" << time_flt_onsets << "\n";
  os << "|ISODEL: FLT: represent_vertex_perm|=" << time_flt_bigmat << "\n";
  os << "|ISODEL: FLT: voronoi_precompute|=" << time_flt_vipc << "\n";
  os << "|ISODEL: FLT: adj_inequality|=" << time_flt_ineq << "\n";
  os << "|ISODEL: FLT: setup|=" << time_flt_setup << "\n";
  os << "|ISODEL: FLT: symbols|=" << time_flt_symbols << "\n";
  os << "|ISODEL: FLT: adjacencies|=" << time_flt_adjacencies << "\n";
  os << "|ISODEL: FLT: get_matching|=" << time_flt_matching << "\n";
  os << "|ISODEL: FLT: get_matching_old|=" << time_flt_matching_old << "\n";
  os << "|ISODEL: FLT: symbol_position|=" << time_flt_symbolpos << "\n";
  os << "|ISODEL: FLT: bigmat_chain|=" << time_flt_chain << "\n";
  os << "|ISODEL: FLT: face_transport|=" << time_flt_transport << "\n";
  os << "|ISODEL: FLT: full|=" << time_flt_full << "\n";
#endif
  return {std::move(l_dels)};
}

FullNamelist NAMELIST_GetStandard_COMPUTE_LATTICE_IsoDelaunayDomains() {
  std::map<std::string, SingleBlock> ListBlock;
  // SYSTEM
  ListBlock["SYSTEM"] = SINGLEBLOCK_Get_System();
  // DATA
  {
    std::map<std::string, std::string> ListStringValues;
    ListStringValues["arithmetic"] = "gmp";
    ListStringValues["FileDualDescription"] = "unset";
    ListStringValues["CommonGramMat"] = "unset";
    // When set to a directory/prefix P (not "unset"), each enumerated
    // iso-Delaunay domain is written as a boost text-archive of an
    // IsoDelaunayDomain to the file P<i> (i = 0, 1, ...), consumable by
    // LATT_AnalysisIsoDelaunay.
    ListStringValues["PrefixIsoDelaunayDomains"] = "unset";
    ListStringValues["CVPmethod"] = "SVexact";
    std::map<std::string, bool> ListBoolValues;
    // Default false keeps the historical behaviour of enumerating only
    // the top-dimensional iso-Delaunay domains. Set to true to run the
    // full L-type complex enumeration (top-dim + lower-dim cells, with
    // the positive-definite filter applied per cell).
    ListBoolValues["compute_full_dimensional"] = false;
    SingleBlock BlockDATA;
    BlockDATA.setListStringValues(ListStringValues);
    BlockDATA.setListBoolValues(ListBoolValues);
    ListBlock["DATA"] = BlockDATA;
  }
  // TSPACE
  ListBlock["TSPACE"] = SINGLEBLOCK_Get_Tspace_Description();
  // Merging all data
  return FullNamelist(ListBlock);
}

/*
  An iso-Delaunay domain is a cone, so its interior point is only defined up
  to a positive factor: the primitive integral representative is taken, which
  makes the stabilizer, the equivalence and the invariant of the domain run
  over the ring. SHV is the invariant vector family of that Gram matrix, which
  is integral as well.
 */
template <typename T, typename Tint, typename Tgroup> struct IsoDelaunayDomain {
  DelaunayTesselationIneq<Tint, Tgroup> DT;
  MyMatrix<Tint> GramMat;
  MyMatrix<Tint> SHV;
};

template <typename T, typename Tint, typename Tgroup>
void WriteEntryGAP(std::ostream &os_out,
                   IsoDelaunayDomain<T, Tint, Tgroup> const &ent) {
  os_out << "rec(DT:=";
  WriteEntryGAP(os_out, ent.DT);
  os_out << ", GramMat:=" << StringMatrixGAP(ent.GramMat) << ")";
}

template <typename T, typename Tint, typename Tgroup>
void WriteEntryPYTHON(std::ostream &os_out,
                      IsoDelaunayDomain<T, Tint, Tgroup> const &ent) {
  os_out << "{\"DT\":";
  WriteEntryPYTHON(os_out, ent.DT);
  os_out << ", \"GramMat\":" << StringMatrixPYTHON(ent.GramMat) << "}";
}

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, IsoDelaunayDomain<T, Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("DT", eRec.DT);
  ar &make_nvp("GramMat", eRec.GramMat);
}
} // namespace boost::serialization

template <typename T, typename Tint, typename Tgroup>
struct IsoDelaunayDomain_AdjI {
  // The defining inequality of the crossed wall, primitive integral. The
  // second scalar T enters through DT_gram, whose interior Gram matrix and
  // SHV are genuinely rational.
  MyVector<Tint> V;
  IsoDelaunayDomain<T, Tint, Tgroup> DT_gram;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar,
                      IsoDelaunayDomain_AdjI<T, Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("V", eRec.V);
  ar &make_nvp("DT_gram", eRec.DT_gram);
}
} // namespace boost::serialization

template <typename Tint> struct IsoDelaunayDomain_AdjO {
  MyVector<Tint> V;
  MyMatrix<Tint> eBigMat;
};

template <typename Tint>
void WriteEntryGAP(std::ostream &os_out,
                   IsoDelaunayDomain_AdjO<Tint> const &ent) {
  os_out << "rec(V:=" << StringVectorGAP(ent.V)
         << " eBigMat:=" << StringMatrixGAP(ent.eBigMat) << ")";
}

template <typename Tint>
void WriteEntryPYTHON(std::ostream &os_out,
                      IsoDelaunayDomain_AdjO<Tint> const &ent) {
  os_out << "{\"V\":" << StringVectorPYTHON(ent.V)
         << " \"eBigMat\":" << StringMatrixPYTHON(ent.eBigMat) << "}";
}

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, IsoDelaunayDomain_AdjO<Tint> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("V", eRec.V);
  ar &make_nvp("eBigMat", eRec.eBigMat);
}
} // namespace boost::serialization

/*
  We do not use the ComputeInvariantDelaunay function since it requires having
  the GramMatrix, which would be an additional computation.
  ---
  BUT we should have some invariants coming from the Tspace and right now we
  have none.
 */
template <typename T, typename Tint, typename Tgroup>
size_t ComputeInvariantIsoDelaunayDomain(
    [[maybe_unused]] DataIsoDelaunayDomains<T, Tint, Tgroup> &data,
    size_t const &seed,
    DelaunayTesselationIneq<Tint, Tgroup> const &DT,
    [[maybe_unused]] std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
  std::map<size_t, size_t> map_delaunays;
  auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
    seed ^= new_hash + 0x9e3779b8 + (seed << 6) + (seed >> 2);
  };
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_PRESENTATION
  std::unordered_map<size_t, size_t> map_vert;
  std::unordered_map<TintGroup, size_t> map_ord_grp;
#endif
  size_t seed_delaunay = 10;
  for (auto &eDel : DT.l_dels) {
    size_t hash = seed_delaunay;
    size_t hash_nVert = eDel.EXT.rows();
    TintGroup ord_grp = eDel.GRP.size();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_PRESENTATION
    map_vert[hash_nVert] += 1;
    map_ord_grp[ord_grp] += 1;
#endif
    size_t hash_grp = std::hash<TintGroup>()(ord_grp);
    combine_hash(hash, hash_nVert);
    combine_hash(hash, hash_grp);
    std::map<size_t, size_t> map_facesiz;
    for (auto &eAdj : eDel.ListAdj) {
      size_t cnt = eAdj.eInc.count();
      map_facesiz[cnt] += 1;
    }
    for (auto &kv : map_facesiz) {
      combine_hash(hash, kv.first);
      combine_hash(hash, kv.second);
    }
    map_delaunays[hash] += 1;
  }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_PRESENTATION
  std::string str_vert = "ISODEL: map_vert =";
  for (auto &kv : map_vert) {
    str_vert +=
        " [" + std::to_string(kv.first) + "," + std::to_string(kv.second) + "]";
  }
  std::string str_ord_grp = "ISODEL: map_ord_grp =";
  for (auto &kv : map_ord_grp) {
    str_ord_grp +=
        " [" + std::format("{}", kv.first) + "," + std::to_string(kv.second) +
        "]";
  }
  os << str_vert << "\n";
  os << str_ord_grp << "\n";
#endif
  size_t hash_ret = seed;
  for (auto &kv : map_delaunays) {
    combine_hash(hash_ret, kv.first);
    combine_hash(hash_ret, kv.second);
  }
  return hash_ret;
}

template <typename T, typename Tint, typename Tgroup>
IsoDelaunayDomain<T, Tint, Tgroup>
GetInitialIsoDelaunayDomain(DataIsoDelaunayDomains<T, Tint, Tgroup> &data) {
  if (!is_integrally_saturated_matrix_space(data.LinSpa.ListMat)) {
    std::cerr << "ISODEL: The space should be integral and fully span it\n";
    throw TerminalException{1};
  }
  std::ostream &os = data.rddo.os;
  DelaunayTesselation<Tint, Tgroup> DT =
      GetInitialGenericDelaunayTesselation(data);
  std::vector<std::vector<Tint>> ListGramRing =
      GetListGramRing(data.LinSpa.ListLineMat);
  DelaunayTesselationIneq<Tint, Tgroup> DTI =
      BuildDelaunayTesselationIneq(DT, ListGramRing, os);
  MyMatrix<T> M = GetInteriorGramMatrix(data.LinSpa, DTI, os);
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(M, os);
  MyMatrix<Tint> M_ring = RemoveFractionMatrixPlusCoeffRing(M).TheMat;
  return {std::move(DTI), std::move(M_ring), std::move(SHV)};
}

template <typename T, typename Tint, typename Tgroup>
struct IsoDelaunayDomain_Obj {
  IsoDelaunayDomain<T, Tint, Tgroup> DT_gram;
  std::vector<FullAdjInfo<Tint>> ListIneqRed;
  Tgroup GRPperm;
};

template <typename T, typename Tint, typename Tgroup>
void WriteBasicEntryGAP(std::ostream &os_out,
                        IsoDelaunayDomain_Obj<T, Tint, Tgroup> const &ent) {
  os_out << "DT_gram:=";
  WriteEntryGAP(os_out, ent.DT_gram);
  //
  os_out << ", ListIneqRed:=[";
  bool IsFirst = true;
  for (auto &eFullAI : ent.ListIneqRed) {
    if (!IsFirst) {
      os_out << ",";
    }
    IsFirst = false;
    WriteEntryGAP(os_out, eFullAI);
  }
  os_out << "]";
  //
  os_out << ", GRPperm:=" << ent.GRPperm.GapString();
}

template <typename T, typename Tint, typename Tgroup>
void WriteEntryGAP(std::ostream &os_out,
                   IsoDelaunayDomain_Obj<T, Tint, Tgroup> const &ent) {
  os_out << "rec(";
  WriteBasicEntryGAP(os_out, ent);
  os_out << ")";
}

template <typename T, typename Tint, typename Tgroup>
void WriteDetailedEntryGAP(std::ostream &os_out,
                           DataIsoDelaunayDomains<T, Tint, Tgroup> const &data,
                           IsoDelaunayDomain_Obj<T, Tint, Tgroup> const &ent,
                           std::ostream &os) {
  os_out << "rec(";
  //  WriteBasicEntryGAP(os_out, ent);
  os_out << "GRPpermSize:=" << ent.GRPperm.size();
  int dimSpace = data.LinSpa.ListMat.size();
  int n = ent.DT_gram.GramMat.rows();
  MyMatrix<T> FAC =
      UniversalMatrixConversion<T, Tint>(GetFACineq(ent.ListIneqRed));
  int nbIneq = FAC.rows();
  os_out << ", n_ineq_red:=" << nbIneq;
  os_out << ", det:=" << DeterminantMat(ent.DT_gram.GramMat);
  os_out << ", n_shv:=" << ent.DT_gram.SHV.rows();
  //
  MyMatrix<T> EXT = DirectDualDescription_mat(FAC, os);
  int n_row = EXT.rows();
  std::map<int, size_t> map_rank;
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyMatrix<T> RayMat = ZeroMatrix<T>(n, n);
    for (int u = 0; u < dimSpace; u++) {
      MatAddMul(RayMat, EXT(i_row, u), data.LinSpa.ListMat[u]);
    }
    int rnk = RankMat(RayMat);
    map_rank[rnk] += 1;
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
  os_out << ", ListRank:=";
  write_map_val(map_rank);
  //
  std::map<int, size_t> map_nb_vert;
  for (auto &eDel : ent.DT_gram.DT.l_dels) {
    int nb_vert = eDel.EXT.rows();
    map_nb_vert[nb_vert] += 1;
  }
  os_out << ", ListNbVert:=";
  write_map_val(map_nb_vert);
  //
  os_out << ")";
}

template <typename T, typename Tint, typename Tgroup>
void WriteEntryPYTHON(std::ostream &os_out,
                      IsoDelaunayDomain_Obj<T, Tint, Tgroup> const &ent) {
  os_out << "{\"DT_gram\":";
  WriteEntryPYTHON(os_out, ent.DT_gram);
  //
  os_out << ", \"ListIneqRed\":[";
  bool IsFirst = true;
  for (auto &eFullAI : ent.ListIneqRed) {
    if (!IsFirst) {
      os_out << ",";
    }
    IsFirst = false;
    WriteEntryPYTHON(os_out, eFullAI);
  }
  os_out << "]";
  //
  os_out << ", \"GRPperm\":" << ent.GRPperm.PythonString();
  //
  os_out << "}";
}

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, IsoDelaunayDomain_Obj<T, Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("DT_gram", eRec.DT_gram);
  ar &make_nvp("ListIneqRed", eRec.ListIneqRed);
  ar &make_nvp("GRPperm", eRec.GRPperm);
}
} // namespace boost::serialization




template <typename T, typename Tint, typename Tgroup>
struct ResultDelaunayAdj {
  std::vector<FullAdjInfo<Tint>> ListIneqRed;
  Tgroup GRPperm;
  std::vector<IsoDelaunayDomain_AdjI<T, Tint, Tgroup>> l_adj;
};

/*
  The combinatorial part of the adjacency computation: the reduced defining
  inequalities, the permutation group acting on them and, for each orbit
  representative, whether the facet leads to a genuine (positive-definite)
  adjacent iso-Delaunay domain. This deliberately stops short of the
  expensive FlippingLtype: l_test[k] tells whether ListIneqRed[k] is
  flippable, but the flip itself is left to get_adjacent. Callers that only
  need one neighbour (e.g. RandomWalkIsoDelaunay) thus avoid flipping every
  facet.
 */
template <typename T, typename Tint, typename Tgroup>
struct PreResultDelaunayAdj {
  std::vector<FullAdjInfo<Tint>> ListIneqRed;
  Tgroup GRPperm;
  std::vector<bool> l_test;
};

/*
  Number of extreme rays r of the L-type domain of x (interpreted as
  vectors in the t-space) such that the Gram matrix
    sum_u r_u * LinSpa.ListMat[u]
  is NOT of full rank. This is the objective function minimized by the
  random walk in LookForFullRankRayDomain: 0 means every ray of the
  domain corresponds to a full-rank Gram matrix.
 */
template <typename T, typename Tint, typename Tgroup>
int CountNonFullRankRays(IsoDelaunayDomain<T, Tint, Tgroup> const &x,
                         DataIsoDelaunayDomains<T, Tint, Tgroup> &data,
                         ThompsonSamplingHeuristic<typename Tgroup::Tint> &eTS,
                         std::ostream &os) {
  int n = data.LinSpa.n;
  int dimSpace = data.LinSpa.ListMat.size();
  MicrosecondTime time_non;
  std::vector<FullAdjInfo<Tint>> ListIneq =
      ComputeListIneqFromTesselationIneq(x.DT);
  os << "|ISODEL: CountNonFullRankRays, ComputeListIneqFromTesselationIneq|=" << time_non << "\n";
  MyMatrix<T> FAC =
      UniversalMatrixConversion<T, Tint>(GetFACineq(ListIneq));
  os << "|ISODEL: CountNonFullRankRays, GetFACineq|=" << time_non << "\n";
  std::vector<int> ListIrred = get_non_redundant_indices(FAC, os);
  os << "|ISODEL: CountNonFullRankRays, get_non_redundant_indices|=" << time_non << "\n";
  MyMatrix<T> FACred = SelectRow(FAC, ListIrred);
  MyMatrix<T> EXT = DirectDualDescription_mat_ts(FACred, eTS, os);
  os << "|ISODEL: CountNonFullRankRays, DirectDualDescription_mat|=" << time_non << "\n";
  int n_row = EXT.rows();
  int count = 0;
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyMatrix<T> RayMat = ZeroMatrix<T>(n, n);
    for (int u = 0; u < dimSpace; u++) {
      MatAddMul(RayMat, EXT(i_row, u), data.LinSpa.ListMat[u]);
    }
    int delta = n - RankMat(RayMat);
    count += delta;
  }
  os << "|ISODEL: CountNonFullRankRays, count|=" << time_non << "\n";
  return count;
}

/*
  The stabilizer of the domain's Gram matrix in the T-space, computed over the
  ring from the invariant vector family already carried by x.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
get_stabilizer_gens(IsoDelaunayDomain<T, Tint, Tgroup> const &x,
                    DataIsoDelaunayDomains<T, Tint, Tgroup> const &data,
                    std::ostream &os) {
  Result_ComputeStabilizer_SHV<Tint, Tgroup> result =
      LINSPA_ComputeStabilizer_SHV<Tint, Tgroup>(
          data.LinSpaRing, x.GramMat, x.SHV, data.CommonGramMatRing, os);
  return result.get_list_matrix(x.SHV, x.GramMat, data.LinSpaRing, os);
}

/*
  Combinatorial pre-computation of the adjacencies: everything except the
  expensive per-facet FlippingLtype. See PreResultDelaunayAdj.

  The kernel takes the generators of the stabilizer of x.GramMat as an
  argument instead of computing them, so that the enumeration of the
  periodic iso-Delaunay domains can pass the subgroup that also preserves
  its point set: the wall orbits are the orbits under that group.
 */
template <typename T, typename Tint, typename Tgroup>
PreResultDelaunayAdj<T, Tint, Tgroup>
get_pre_result_delaunay_adj_kernel(IsoDelaunayDomain<T, Tint, Tgroup> const &x,
                            DataIsoDelaunayDomains<T, Tint, Tgroup> &data,
                            std::vector<MyMatrix<Tint>> const &ListGenTot) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::ostream &os = data.rddo.os;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  MicrosecondTime time_f_adj;
#endif
  int n = data.LinSpa.n;
  int dimSpace = data.LinSpa.ListMat.size();
  // When a common Gram matrix G is imposed, the enumeration is restricted to
  // the "star of G": the iso-Delaunay domains whose closure contains G. G is a
  // point of the current domain's closure (the seed passes the acceptability
  // check and this invariant is preserved on every accepted adjacency below).
  // Crossing a wall reaches an adjacent domain that still contains G if and
  // only if G lies exactly on that wall, i.e. the wall's defining inequality
  // vanishes at G. The stored inequality is ScalarCanonicalizationVector(V), a
  // nonzero multiple of the Voronoi inequality V, so V.dot(g_vec)==0 is tested
  // sign-independently as eIneq.dot(g_vec)==0. g_vec is the T-space coordinate
  // vector of G. The wall selection itself is applied to l_test below; the
  // actual flip is performed later by get_adjacent.
  std::optional<MyVector<T>> g_vec;
  if (data.CommonGramMat) {
    g_vec = LINSPA_GetVectorOfMatrixExpression(data.LinSpa, *data.CommonGramMat);
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN
    // Invariant: G must lie in the closure of the current domain. Check it with
    // the same acceptability predicate that defines domain membership for the
    // seed, which is orientation-free (unlike the per-wall inequality sign).
    // The T-space coordinate vector of G is g_vec itself.
    std::vector<std::vector<Tint>> ListGramRing =
        GetListGramRing(data.LinSpa.ListLineMat);
    for (auto &eDel : x.DT.l_dels) {
      if (!IsDelaunayAcceptableForGramMat(eDel.EXT, ListGramRing, *g_vec,
                                          os)) {
        std::cerr << "ISODEL: f_adj, CommonGramMat G is not in the closure of "
                     "the current domain. The seed domain must contain G.\n";
        throw TerminalException{1};
      }
    }
#endif
  }
  std::vector<FullAdjInfo<Tint>> ListIneq =
    ComputeListIneqFromTesselationIneq(x.DT);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: f_adj, ComputeListIneqFromTesselationIneq|=" << time_f_adj
       << "\n";
#endif
  // Compute the irredundant ones as well as the l_ineq / map_ineq. The
  // inequalities are canonical primitive integral vectors; the field enters
  // only in the redundancy elimination, the interior points of the walls and
  // the action of the rational stabilizer generators.
  MyMatrix<Tint> FAC = GetFACineq(ListIneq);
  MyMatrix<T> FAC_T = UniversalMatrixConversion<T, Tint>(FAC);
  std::vector<int> ListIrred = get_non_redundant_indices(FAC_T, os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: f_adj, get_non_redundant_indices|=" << time_f_adj << "\n";
#endif
  size_t nbIrred = ListIrred.size();
  MyMatrix<Tint> FACred = SelectRow(FAC, ListIrred);
  MyMatrix<T> FACred_T = SelectRow(FAC_T, ListIrred);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: f_adj: |FAC|=" << FAC.rows() << " / " << FAC.cols()
     << " nbIrred=" << nbIrred << "\n";
  os << "ISODEL: f_adj, x.GramMat=\n";
  WriteMatrix(os, x.GramMat);
#endif
  std::vector<MyVector<T>> l_ineq;
  std::unordered_map<MyVector<Tint>, size_t> map_ineq;
  for (size_t i = 0; i < nbIrred; i++) {
    MyVector<Tint> eV = GetMatrixRow(FACred, i);
    l_ineq.push_back(UniversalVectorConversion<T, Tint>(eV));
    map_ineq[eV] = i;
  }
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: f_adj, l_ineq - map_ineq|=" << time_f_adj << "\n";
#endif
  // The permutations of the facets induced by the stabilizer generators.
  std::vector<Telt> ListPermGens;
  for (auto &eGenTot : ListGenTot) {
    // The generators are integral, but their expression in the basis of the
    // T-space is rational in general, so the action is computed over T.
    MyMatrix<T> eGenTot_T = UniversalMatrixConversion<T, Tint>(eGenTot);
    MyMatrix<T> MatSpace = matrix_in_t_space(eGenTot_T, data.LinSpa);
    std::vector<Tidx> l_pos(nbIrred);
    for (size_t i = 0; i < nbIrred; i++) {
      MyVector<T> const &eV = l_ineq[i];
      // There are actually two transpose here:
      // * One from going from action in the space to action on the dual
      // * One for going from row to column action
      MyVector<T> eVimg = MatSpace * eV;
      MyVector<Tint> eVimg_red =
          RemoveFractionVectorPlusCoeffRing(eVimg).TheVect;
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN
      if (map_ineq.find(eVimg_red) == map_ineq.end()) {
        std::cerr << "ISODEL: eGenTot=\n";
        WriteMatrix(std::cerr, eGenTot);
        std::cerr << "ISODEL: i=" << i << " eVimg=" << StringVector(eVimg)
                  << "\n";
        std::cerr << "ISODEL: i=" << i << " eVimg=" << StringVector(eVimg_red)
                  << "\n";
        std::cerr << "ISODEL: FACred=\n";
        WriteMatrix(std::cerr, FACred);
        std::cerr << "ISODEL: The entry eVimg is missing from the map\n";
        throw TerminalException{1};
      }
#endif
      size_t pos = map_ineq.at(eVimg_red);
      l_pos[i] = pos;
    }
    Telt ePermGen(l_pos);
    ListPermGens.push_back(ePermGen);
  }
  Tgroup GRPperm = Tgroup(ListPermGens, nbIrred);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: f_adj, GRPperm|=" << time_f_adj << "\n";
#endif

  std::vector<size_t> l_idx = DecomposeOrbitPoint_FullRepr(GRPperm);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: f_adj: |GRPperm|=" << GRPperm.size()
     << " nbIrred=" << nbIrred << " |l_idx|=" << l_idx.size() << "\n";
#endif
  std::vector<FullAdjInfo<Tint>> ListIneqRed;
  std::vector<bool> l_test;
  for (auto &i : l_idx) {
    int idxIrred = ListIrred[i];
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: f_adj, i=" << i << " idxIrred=" << idxIrred << "\n";
#endif
    MyVector<T> TestPt = GetSpaceInteriorPointFacet(FACred_T, i, os);
    MyMatrix<T> TestMat = ZeroMatrix<T>(n, n);
    for (int u = 0; u < dimSpace; u++) {
      MatAddMul(TestMat, TestPt(u), data.LinSpa.ListMat[u]);
    }
    bool test = IsPositiveDefinite(TestMat, os);
    FullAdjInfo<Tint> eRecIneq = ListIneq[idxIrred];
    ListIneqRed.push_back(eRecIneq);
    // Star-of-G restriction: when a common Gram matrix G is imposed, a wall is
    // flippable only if G lies on it, i.e. the adjacent domain still contains G
    // in its closure. G lies on this wall iff its Voronoi inequality vanishes at
    // G (sign-independent, see above). Folding this into l_test lets master's
    // downstream (get_result_delaunay_adj / RandomWalkIsoDelaunay) skip the
    // expensive flip for walls that leave the star of G.
    bool flippable = test;
    if (flippable && g_vec) {
      MyVector<T> eIneq_T = UniversalVectorConversion<T, Tint>(eRecIneq.eIneq);
      T scal = eIneq_T.dot(*g_vec);
      if (scal != 0) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: f_adj, skipping wall i=" << i
           << " not containing CommonGramMat (scal=" << scal << ")\n";
#endif
        flippable = false;
      }
    }
    l_test.push_back(flippable);
  }
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: f_adj, pre-result|=" << time_f_adj << "\n";
#endif
  return {std::move(ListIneqRed), std::move(GRPperm), std::move(l_test)};
}

template <typename T, typename Tint, typename Tgroup>
PreResultDelaunayAdj<T, Tint, Tgroup>
get_pre_result_delaunay_adj(IsoDelaunayDomain<T, Tint, Tgroup> const &x,
                            DataIsoDelaunayDomains<T, Tint, Tgroup> &data) {
  std::ostream &os = data.rddo.os;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  MicrosecondTime time_f_adj;
#endif
  std::vector<MyMatrix<Tint>> ListGenTot = get_stabilizer_gens(x, data, os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: f_adj, LINSPA_ComputeStabilizer|=" << time_f_adj << "\n";
#endif
  return get_pre_result_delaunay_adj_kernel<T, Tint, Tgroup>(x, data,
                                                             ListGenTot);
}

/*
  Perform the actual flip across the facet eRecIneq and build the resulting
  adjacent iso-Delaunay domain. This is the expensive step, isolated so that
  it can be run only for the neighbour(s) actually needed.
 */
template <typename T, typename Tint, typename Tgroup>
IsoDelaunayDomain_AdjI<T, Tint, Tgroup>
get_adjacent(IsoDelaunayDomain<T, Tint, Tgroup> const &x,
             DataIsoDelaunayDomains<T, Tint, Tgroup> &data,
             FullAdjInfo<Tint> const &eRecIneq) {
  std::ostream &os = data.rddo.os;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  MicrosecondTime time_s_adj;
#endif
  std::vector<std::vector<Tint>> ListGramRing =
      GetListGramRing(data.LinSpa.ListLineMat);
  // The flipping runs entirely over the ring: the interior element enters
  // only through the lifting heights, for which a positive fraction-free
  // rescaling changes nothing.
  DelaunayTesselationIneq<Tint, Tgroup> DTIadj =
    FlippingLtype(x.DT, x.GramMat, eRecIneq.ListAdjInfo, ListGramRing,
                  data.rddo.AllArr, os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: s_adj, FlippingLtype|=" << time_s_adj << "\n";
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: After FlippingLtype\n";
#endif
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN
  if (data.CommonGramMat) {
    // Star-of-G restriction: get_adjacent is only called for a wall on which G
    // lies (get_pre_result_delaunay_adj marks it flippable only when the wall's
    // inequality vanishes at G), so the flipped domain must still contain G in
    // its closure. Verify it with the acceptability predicate.
    MyVector<T> g_vec =
        LINSPA_GetVectorOfMatrixExpression(data.LinSpa, *data.CommonGramMat);
    for (auto &eDel : DTIadj.l_dels) {
      if (!IsDelaunayAcceptableForGramMat(eDel.EXT, ListGramRing, g_vec, os)) {
        std::cerr << "ISODEL: get_adjacent, star-of-G: crossed a wall on which "
                     "G lies but the adjacent domain does not contain G\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  std::vector<FullAdjInfo<Tint>> ListIneqAdj =
    ComputeListIneqFromTesselationIneq(DTIadj);
#ifdef SANITY_CHECK_FLIP_INEQ_CONSISTENCY
  // End-to-end check that the incrementally-updated inequalities
  // (Case 1 reused, other cases recomputed) match a full from-scratch
  // recomputation on the flipped tessellation.
  {
    std::vector<FullAdjInfo<Tint>> ListIneqAdjSlow =
      ComputeDefiningIneqIsoDelaunayDomain(
          StripDelaunayTesselationIneq(DTIadj), ListGramRing, os);
    std::set<MyVector<Tint>> setFast, setSlow;
    for (auto &eFullAI : ListIneqAdj) {
      setFast.insert(eFullAI.eIneq);
    }
    for (auto &eFullAI : ListIneqAdjSlow) {
      setSlow.insert(eFullAI.eIneq);
    }
    if (setFast != setSlow) {
      std::cerr << "FLIP_INEQ_CONSISTENCY: SANITY_CHECK failed: the "
                   "incrementally flipped inequalities do not match a "
                   "full recomputation on the flipped tessellation\n";
      throw TerminalException{1};
    }
  }
#endif
  MyMatrix<T> FACadj =
      UniversalMatrixConversion<T, Tint>(GetFACineq(ListIneqAdj));
  MyMatrix<T> M = get_interior_gram_matrix_lp(data.LinSpa, FACadj, os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: s_adj, GetInteriorGramMatrix|=" << time_s_adj << "\n";
#endif
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(M, os);
  MyMatrix<Tint> M_ring = RemoveFractionMatrixPlusCoeffRing(M).TheMat;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
  os << "|ISODEL: s_adj, shv_extract|=" << time_s_adj << "\n";
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: After GetInteriorGramMatrix\n";
#endif
  IsoDelaunayDomain<T, Tint, Tgroup> IsoDelAdj{
    std::move(DTIadj), std::move(M_ring), std::move(SHV)};
  return {eRecIneq.eIneq, std::move(IsoDelAdj)};
}

/*
  Full adjacency computation: the combinatorial pre-result followed by an
  actual flip for every flippable facet. Same kernel / wrapper split as
  get_pre_result_delaunay_adj, and for the same reason.
 */
template <typename T, typename Tint, typename Tgroup>
ResultDelaunayAdj<T, Tint, Tgroup>
get_result_delaunay_adj_kernel(IsoDelaunayDomain<T, Tint, Tgroup> const &x,
                        DataIsoDelaunayDomains<T, Tint, Tgroup> &data,
                        std::vector<MyMatrix<Tint>> const &ListGenTot) {
  PreResultDelaunayAdj<T, Tint, Tgroup> pre =
    get_pre_result_delaunay_adj_kernel<T, Tint, Tgroup>(x, data, ListGenTot);
  std::vector<IsoDelaunayDomain_AdjI<T, Tint, Tgroup>> l_adj;
  for (size_t k = 0; k < pre.l_test.size(); k++) {
    if (pre.l_test[k]) {
      l_adj.push_back(get_adjacent<T, Tint, Tgroup>(x, data, pre.ListIneqRed[k]));
    }
  }
  return {std::move(pre.ListIneqRed), std::move(pre.GRPperm), std::move(l_adj)};
}

template <typename T, typename Tint, typename Tgroup>
ResultDelaunayAdj<T, Tint, Tgroup>
get_result_delaunay_adj(IsoDelaunayDomain<T, Tint, Tgroup> const &x,
                        DataIsoDelaunayDomains<T, Tint, Tgroup> &data) {
  std::vector<MyMatrix<Tint>> ListGenTot =
      get_stabilizer_gens(x, data, data.rddo.os);
  return get_result_delaunay_adj_kernel<T, Tint, Tgroup>(x, data, ListGenTot);
}

/*
  Take n_iter unbiased random steps in the adjacency graph of iso-Delaunay
  domains. Used to escape a local minimum of the ray-rank objective.
 */
template <typename T, typename Tint, typename Tgroup>
IsoDelaunayDomain<T, Tint, Tgroup>
RandomWalkIsoDelaunay(IsoDelaunayDomain<T, Tint, Tgroup> const &x,
                      DataIsoDelaunayDomains<T, Tint, Tgroup> &data,
                      int n_iter, std::ostream &os) {
  IsoDelaunayDomain<T, Tint, Tgroup> Work = x;
  for (int iter = 0; iter < n_iter; iter++) {
    os << "ISODEL: RandomWalkIsoDelaunay iter=" << iter << "/" << n_iter << "\n";
    // Only the combinatorial part is needed to know which facets are
    // flippable; the actual flip is done for the single chosen neighbour.
    PreResultDelaunayAdj<T, Tint, Tgroup> pre =
        get_pre_result_delaunay_adj(Work, data);
    std::vector<size_t> l_flippable;
    for (size_t k = 0; k < pre.l_test.size(); k++) {
      if (pre.l_test[k]) {
        l_flippable.push_back(k);
      }
    }
    int n_adj = l_flippable.size();
    os << "ISODEL: RandomWalkIsoDelaunay n_adj=" << n_adj << "\n";
    if (n_adj == 0) {
      os << "ISODEL: RandomWalkIsoDelaunay, no adjacent domain at iter="
         << iter << ", stopping early\n";
      break;
    }
    int pos = random() % n_adj;
    Work = get_adjacent(Work, data, pre.ListIneqRed[l_flippable[pos]]).DT_gram;
  }
  return Work;
}

/*
  Greedy descent on the ray-rank objective, modelled on
  src_ctype/CTYP_LookForNoFreeVector.cpp::GetLocalFreenessMinimum.

  Strategy:
  * Compute the number of non-full-rank rays of the current domain
    (CountNonFullRankRays).
  * Enumerate adjacent domains and evaluate the same objective on each.
  * If at least one adjacent strictly improves, jump to a (uniformly
    random) one among those tying for the minimum.
  * Otherwise, perform n_walk_steps random adjacency jumps to leave
    the basin.
  * Return as soon as the objective reaches zero, writing the
    corresponding Gram matrix to a fresh file under Prefix. If the
    objective is at most max_s but still positive we also dump the
    Gram matrix (useful for collecting near-misses).
  * Bail out after max_iter outer-loop iterations even if the
    objective never reaches zero. This matters because in low
    dimension (n = 2..5 of the Classic t-space) every iso-Delaunay
    domain provably has some rank-1 extreme ray, so the target is
    unreachable and a strict "return on zero" would loop forever.
 */
template <typename T, typename Tint, typename Tgroup>
void LookForFullRankRayDomain(DataIsoDelaunayDomains<T, Tint, Tgroup> &data,
                              std::string const &Prefix, int const &max_s,
                              int const &n_walk_steps, int const &max_iter,
                              std::ostream &os) {
  IsoDelaunayDomain<T, Tint, Tgroup> Work = GetInitialIsoDelaunayDomain(data);
  // One stateful lrs/cdd/bb Thompson sampler shared across every
  // CountNonFullRankRays call so the posterior keeps accumulating over
  // the whole random walk.
  ThompsonSamplingHeuristic<typename Tgroup::Tint> eTS =
      MakeDualDescThompsonSampler<typename Tgroup::Tint>(os);
  os << "ISODEL: LookForFullRankRayDomain, before CountNonFullRankRays\n";
  int curr_count = CountNonFullRankRays(Work, data, eTS, os);
  os << "ISODEL: LookForFullRankRayDomain, initial curr_count=" << curr_count
     << "\n";
  int iter1 = 0;
  int iter2 = 0;
  int total_iter = 0;
  while (true) {
    if (total_iter >= max_iter) {
      os << "ISODEL: LookForFullRankRayDomain, reached max_iter=" << max_iter
         << " without hitting curr_count=0 (current curr_count=" << curr_count
         << "), returning\n";
      return;
    }
    total_iter++;
    ResultDelaunayAdj<T, Tint, Tgroup> result =
        get_result_delaunay_adj(Work, data);
    int n_adj = result.l_adj.size();
    if (n_adj == 0) {
      std::cerr << "ISODEL: We should have at least one adjacent domain\n";
      throw TerminalException{1};
    }
    std::vector<int> ListCount;
    for (int i = 0; i < n_adj; i++) {
      int c = CountNonFullRankRays(result.l_adj[i].DT_gram, data, eTS, os);
      os << "ISODEL: LookForFullRankRayDomain, i=" << i << "/" << n_adj << " c=" << c << "\n";
      ListCount.push_back(c);
    }
    int the_min = VectorMin(ListCount);
    os << "ISODEL: LookForFullRankRayDomain, the_min=" << the_min << " curr_count=" << curr_count << "\n";
    if (the_min >= curr_count && curr_count > 0) {
      os << "ISODEL: Before RandomWalkIsoDelaunay\n";
      Work = RandomWalkIsoDelaunay(Work, data, n_walk_steps, os);
      curr_count = CountNonFullRankRays(Work, data, eTS, os);
      os << "ISODEL: LookForFullRankRayDomain, iter1=" << iter1
         << " iter2=" << iter2
         << " After RandomWalk curr_count=" << curr_count << "\n";
      iter1++;
      iter2 = 0;
    } else {
      curr_count = the_min;
      std::vector<size_t> ListIdx;
      for (size_t u = 0; u < ListCount.size(); u++) {
        if (ListCount[u] == the_min) {
          ListIdx.push_back(u);
        }
      }
      int n_min = ListIdx.size();
      int pos = random() % n_min;
      os << "ISODEL: LookForFullRankRayDomain, iter1=" << iter1
         << " iter2=" << iter2 << " curr_count=" << curr_count
         << " n_min=" << n_min << " pos=" << pos << "\n";
      Work = result.l_adj[ListIdx[pos]].DT_gram;
      if (curr_count <= max_s) {
        std::string FileOut = FILE_FindAvailableFileFromPrefix(Prefix);
        WriteMatrixFile(FileOut, Work.GramMat);
        os << "ISODEL: LookForFullRankRayDomain, wrote candidate Gram matrix "
              "to "
           << FileOut << " (curr_count=" << curr_count << ")\n";
      }
      if (curr_count == 0) {
        return;
      }
      iter2++;
    }
  }
}


template <typename T, typename Tint, typename Tgroup>
struct DataIsoDelaunayDomainsFunc {
  DataIsoDelaunayDomains<T, Tint, Tgroup> data;
  using Tobj = IsoDelaunayDomain_Obj<T, Tint, Tgroup>;
  using TadjI = IsoDelaunayDomain_AdjI<T, Tint, Tgroup>;
  using TadjO = IsoDelaunayDomain_AdjO<Tint>;
  std::ostream &get_os() { return data.rddo.os; }
  Tobj f_init() {
    IsoDelaunayDomain<T, Tint, Tgroup> IsoDel =
        GetInitialIsoDelaunayDomain(data);
    Tobj x{IsoDel, {}, {}};
    return x;
  }
  size_t f_hash(size_t const &seed, Tobj const &x) {
    std::ostream &os = data.rddo.os;
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
    MicrosecondTime time_hash;
#endif
    int method = 2;
    if (method == 1) {
      // That method does not seem to work very well
      // because it can happen that many triangulations have all Delaunays
      // being simplices with trivial stabilizer. Those would have exactly
      // the same hash invariant.
      return ComputeInvariantIsoDelaunayDomain<T, Tint, Tgroup>(
          data, seed, x.DT_gram.DT, os);
    }
    if (method == 2) {
      size_t hash =
          LINSPA_Invariant_SHV<Tint>(seed, data.LinSpaRing, x.DT_gram.GramMat,
                                     x.DT_gram.SHV, data.CommonGramMatRing, os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
      os << "|ISODEL: f_hash|=" << time_hash << "\n";
#endif
      return hash;
    }
    std::cerr << "No valid method for computing the hash\n";
    throw TerminalException{1};
  }
  std::optional<TadjO> f_repr(Tobj const &x, TadjI const &y) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    data.rddo.os
        << "ISODEL: f_repr, before LINSPA_TestEquivalenceGramMatrix_SHV\n";
#endif
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
    MicrosecondTime time_repr;
#endif
    std::optional<MyMatrix<Tint>> opt =
        LINSPA_TestEquivalenceGramMatrix_SHV<Tint, Tgroup>(
            data.LinSpaRing, x.DT_gram.GramMat, y.DT_gram.GramMat,
            x.DT_gram.SHV, y.DT_gram.SHV, data.CommonGramMatRing,
            data.rddo.os);
#ifdef TIMINGS_ISO_DELAUNAY_DOMAIN
    // The outcome is logged with the time: an invariant that separates well
    // makes almost every call return an equivalence, a weak one wastes the
    // (expensive) test on pairs that are not equivalent.
    data.rddo.os << "|ISODEL: f_repr equiv=" << opt.has_value()
                 << "|=" << time_repr << "\n";
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    data.rddo.os
        << "ISODEL: f_repr, after LINSPA_TestEquivalenceGramMatrix_SHV\n";
#endif
    if (!opt) {
      return {};
    }
    TadjO ret{y.V, *opt};
    return ret;
  }
  std::pair<Tobj, TadjO> f_spann(TadjI const &x) {
    IsoDelaunayDomain<T, Tint, Tgroup> IsoDel = x.DT_gram;
    Tobj x_ret{IsoDel, {}, {}};
    MyMatrix<Tint> eBigMat = IdentityMat<Tint>(data.LinSpa.n);
    TadjO ret{x.V, eBigMat};
    return {std::move(x_ret), std::move(ret)};
  }
  std::optional<std::vector<TadjI>> f_adj(Tobj &x_in) {
    IsoDelaunayDomain<T, Tint, Tgroup> const& x = x_in.DT_gram;
    ResultDelaunayAdj<T,Tint,Tgroup> result = get_result_delaunay_adj(x, data);
    x_in.ListIneqRed = result.ListIneqRed;
    x_in.GRPperm = result.GRPperm;
    return result.l_adj;
  }
  Tobj f_adji_obj(TadjI const &x) { return {x.DT_gram, {}, {}}; }
  size_t f_complexity([[maybe_unused]] Tobj const &x) {
    // It is not clear what we could put. The total number of
    // inequality could be one. But actually the number after reduction
    // by redundancy would be the true metric.
    // But on the other hand, that number is actually not really a complexity.
    // It just tells that the computation last longer but there is no
    // exponential increase in complexity. It is not a dual description.
    return 0;
  }
};

template <typename T, typename Tint, typename Tgroup>
DataIsoDelaunayDomains<T, Tint, Tgroup>
get_data_isodelaunay_domains(FullNamelist const &eFull,
                             PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
                             std::ostream &os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
  auto get_common = [&]() -> std::optional<MyMatrix<T>> {
    std::string const &CommonGramMat = BlockDATA.get_string("CommonGramMat");
    if (CommonGramMat == "unset") {
      return {};
    }
    MyMatrix<T> eMat = ReadMatrixFile<T>(CommonGramMat);
    return eMat;
  };
  std::optional<MyMatrix<T>> CommonGramMat = get_common();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "We have CommonGramMat\n";
#endif
  //
  LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint, Tgroup>(BlockTSPACE, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "We have LinSpa\n";
#endif
  //
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "We have rddo\n";
#endif
  //
  LinSpaceMatrix<Tint> LinSpaRing = LINSPA_GetRingVersion(LinSpa);
  std::optional<MyMatrix<Tint>> CommonGramMatRing =
      GetCommonGramMatRing<Tint, T>(CommonGramMat);
  DataIsoDelaunayDomains<T, Tint, Tgroup> data{
      LinSpa, std::move(LinSpaRing), std::move(rddo), CommonGramMat,
      CommonGramMatRing};
  return data;
}

// clang-format off
#endif  // SRC_LATT_ISODELAUNAYDOMAINS_H_
// clang-format on

// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ISODELAUNAYDOMAINS_H_
#define SRC_LATT_ISODELAUNAYDOMAINS_H_

// clang-format off
#include "POLY_LinearProgramming.h"
#include "ShortestUniversal.h"
#include "Positivity.h"
#include "LatticeDelaunay.h"
#include "Tspace_General.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <memory>
#include <utility>
#include <limits>
// clang-format on

#ifdef DEBUG
#define DEBUG_ISO_DELAUNAY_DOMAIN
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ISO_DELAUNAY_DOMAIN
#endif

template <typename T, typename Tint, typename Tgroup>
struct DataIsoDelaunayDomains {
  // The Linear space of matrices
  LinSpaceMatrix<T> LinSpa;
  // The database of dual descriptions.
  RecordDualDescOperation<T, Tgroup> rddo;
  // The matrices common to all the IsoDelaunay domains if one is chosen
  // to be common.
  std::optional<MyMatrix<T>> CommonGramMat;
};

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

template <typename T> struct VoronoiInequalityPreComput {
  //  MyMatrix<Tvert> VertBasis;
  int n;
  Face f_basis;
  MyMatrix<T> VertBasis_T;
  MyMatrix<T> VertBasisInv_T;
  std::vector<MyVector<T>> VertBasisRed_T;
};

/*
  We implement the hash of a Delaunay tessellationn. The constraint is that two
  different but equivalent tessellations, must have the same hash. Hopefully,
  this is not a problem since we have many invariants:
  * The number of vertices of the orbit representatives of Delaunay polytopes.
  * The size of their automorphism groups
  * The number of vertices of the orbit representative of their facets.
 */
template <typename Tvert, typename Tgroup>
size_t ComputeInvariantDelaunayTessellation(
    DelaunayTesselation<Tvert, Tgroup> const &DT, size_t const &seed,
    [[maybe_unused]] std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
  std::map<size_t, size_t> map;
  auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
    seed ^= new_hash + 0x9e3779b8 + (seed << 6) + (seed >> 2);
  };
  for (auto &e_del : DT.l_dels) {
    std::map<size_t, size_t> map_siz;
    for (auto &eAdj : e_del.ListAdj) {
      size_t siz = eAdj.eInc.count();
      map_siz[siz] += 1;
    }
    size_t hash_del = 123;
    size_t hash1 = std::hash<int>()(e_del.EXT.rows());
    size_t hash2 = std::hash<TintGroup>()(e_del.GRP.order());
    combine_hash(hash_del, hash1);
    combine_hash(hash_del, hash2);
    for (auto &kv : map_siz) {
      combine_hash(hash_del, kv.first);
      combine_hash(hash_del, kv.second);
    }
    map[hash_del] += 1;
  }
  size_t hash_ret = seed;
  for (auto &kv : map) {
    combine_hash(hash_ret, kv.first);
    combine_hash(hash_ret, kv.second);
  }
  return hash_ret;
}

template <typename T, typename Tvert>
VoronoiInequalityPreComput<T>
BuildVoronoiIneqPreCompute(MyMatrix<Tvert> const &EXT,
                           [[maybe_unused]] std::ostream &os) {
  int n = EXT.cols() - 1;
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tvert>(EXT);
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(EXT_T);
  std::vector<int> ListRowSelect = eSelect.ListRowSelect;
  Face f_basis(EXT_T.rows());
  for (auto &eIdx : ListRowSelect) {
    f_basis[eIdx] = 1;
  }
  MyMatrix<T> VertBasis_T = SelectRow(EXT_T, ListRowSelect);
  MyMatrix<T> VertBasisInv_T = TransposedMat(Inverse(VertBasis_T));
  std::vector<MyVector<T>> VertBasisRed_T;
  for (int i = 0; i <= n; i++) {
    MyVector<T> V(n);
    for (int u = 0; u < n; u++) {
      V(u) = VertBasis_T(i, u + 1);
    }
    VertBasisRed_T.push_back(V);
  }
  return {n, f_basis, std::move(VertBasis_T), std::move(VertBasisInv_T),
          std::move(VertBasisRed_T)};
}

template <typename T, typename Tvert>
MyVector<T> VoronoiLinearInequality(VoronoiInequalityPreComput<T> const &vipc,
                                    MyVector<Tvert> const &TheVert,
                                    std::vector<std::vector<T>> const &ListGram,
                                    [[maybe_unused]] std::ostream &os) {
  int n = vipc.n;
  int dimSpace = ListGram.size();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_DISABLE
  os << "ISODEL: VoronoiLinearInequality n=" << n << " dimSpace=" << dimSpace
     << "\n";
#endif
  MyVector<T> TheVert_T = UniversalVectorConversion<T, Tvert>(TheVert);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_DISABLE
  os << "ISODEL: VoronoiLinearInequality |TheVert_T|=" << TheVert_T.size()
     << "\n";
  os << "ISODEL: VoronoiLinearInequality |Tvipc.VertBasisInv_T}="
     << vipc.VertBasisInv_T.rows() << " / " << vipc.VertBasisInv_T.cols()
     << "\n";
#endif
  MyVector<T> B = vipc.VertBasisInv_T * TheVert_T;
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_DISABLE
  os << "ISODEL: VoronoiLinearInequality We have B\n";
#endif
  MyVector<T> Ineq(dimSpace);
  MyVector<T> TheVertRed(n);
  for (int u = 0; u < n; u++) {
    TheVertRed(u) = TheVert_T(u + 1);
  }
  int iGram = 0;
  for (auto &eLineMat : ListGram) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_DISABLE
    os << "ISODEL: VoronoiLinearInequality Before EvaluateLineVector "
          "(first)\n";
#endif
    T val = EvaluateLineVector(eLineMat, TheVertRed);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_DISABLE
    os << "ISODEL: VoronoiLinearInequality After EvaluateLineVector (first)\n";
#endif
    for (int k = 0; k <= n; k++) {
      val -= B(k) * EvaluateLineVector(eLineMat, vipc.VertBasisRed_T[k]);
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_DISABLE
    os << "ISODEL: VoronoiLinearInequality Summation of the entries\n";
#endif
    Ineq(iGram) = val;
    iGram++;
  }
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN_DISABLE
  for (int iVert = 0; iVert <= n; iVert++) {
    for (int u = 0; u < n; u++) {
      TheVertRed(u) = vipc.VertBasis_T(iVert, u + 1);
    }
    for (int u = 0; u <= n; u++) {
      TheVert_T(u) = vipc.VertBasis_T(iVert, u);
    }
    MyVector<T> B = vipc.VertBasisInv_T * TheVert_T;
    for (auto &eLineMat : ListGram) {
      T val = EvaluateLineVector(eLineMat, TheVertRed);
      for (int k = 0; k <= n; k++) {
        val -= B(k) * EvaluateLineVector(eLineMat, vipc.VertBasisRed_T[k]);
      }
      if (val != 0) {
        std::cerr << "For the vertex, the condition should be void and so we "
                     "should get 0\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  return Ineq;
}

template <typename T, typename Tvert>
bool IsDelaunayPolytopeInducingEqualities(MyMatrix<Tvert> const &EXT,
                                          LinSpaceMatrix<T> const &LinSpa,
                                          std::ostream &os) {
  int n_row = EXT.rows();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_DISABLE
  os << "ISODEL: IsDelaunayPolytopeInducingEqualities, before "
        "BuildVoronoiIneqPreCompute\n";
#endif
  VoronoiInequalityPreComput<T> vipc =
      BuildVoronoiIneqPreCompute<T, Tvert>(EXT, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_DISABLE
  os << "ISODEL: IsDelaunayPolytopeInducingEqualities, we have vipc\n";
#endif
  for (int i_row = 0; i_row < n_row; i_row++) {
    if (vipc.f_basis[i_row] == 0) {
      MyVector<Tvert> TheVert = GetMatrixRow(EXT, i_row);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_DISABLE
      os << "ISODEL: IsDelaunayPolytopeInducingEqualities, Before "
            "VoronoiLinearInequality\n";
#endif
      MyVector<T> V =
          VoronoiLinearInequality(vipc, TheVert, LinSpa.ListLineMat, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN_DISABLE
      os << "ISODEL: IsDelaunayPolytopeInducingEqualities, After "
            "VoronoiLinearInequality\n";
#endif
      if (!IsZeroVector(V)) {
        return true;
      }
    }
  }
  return false;
}

template <typename T, typename Tvert>
bool IsDelaunayAcceptableForGramMat(MyMatrix<Tvert> const &EXT,
                                    LinSpaceMatrix<T> const &LinSpa,
                                    MyMatrix<T> const &TestGram,
                                    std::ostream &os) {
  int n_row = EXT.rows();
  MyVector<T> TestV = LINSPA_GetVectorOfMatrixExpression(LinSpa, TestGram);
  VoronoiInequalityPreComput<T> vipc =
      BuildVoronoiIneqPreCompute<T, Tvert>(EXT, os);
  for (int i_row = 0; i_row < n_row; i_row++) {
    if (vipc.f_basis[i_row] == 0) {
      MyVector<Tvert> TheVert = GetMatrixRow(EXT, i_row);
      MyVector<T> V =
          VoronoiLinearInequality(vipc, TheVert, LinSpa.ListLineMat, os);
      T TheScal = V.dot(TestV);
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
  os_out << "{\"iOrb\":" << (eai.iOrb + 1) << ", \"i_adj\":" << (eai.i_adj + 1) << "}";
}

namespace boost::serialization {
template <class Archive>
inline void serialize(Archive &ar, AdjInfo &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("iOrb", eRec.iOrb);
  ar &make_nvp("i_adj", eRec.i_adj);
}
} // namespace boost::serialization

template <typename T> struct FullAdjInfo {
  MyVector<T> eIneq;
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
  os_out << "{\"eIneq\":" << StringVectorPYTHON(ent.eIneq) << ", \"ListAdjInfo\":[";
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
  Compute the defining inequalities of an iso-Delaunay domain
 */
template <typename T, typename Tvert, typename Tgroup>
std::vector<FullAdjInfo<T>> ComputeDefiningIneqIsoDelaunayDomain(
    DelaunayTesselation<Tvert, Tgroup> const &DT,
    std::vector<std::vector<T>> const &ListGram, std::ostream &os) {
  std::unordered_map<MyVector<T>, std::vector<AdjInfo>> map;
  int n_del = DT.l_dels.size();
  for (int i_del = 0; i_del < n_del; i_del++) {
    int n_adj = DT.l_dels[i_del].ListAdj.size();
    VoronoiInequalityPreComput<T> vipc =
        BuildVoronoiIneqPreCompute<T, Tvert>(DT.l_dels[i_del].EXT, os);
    ContainerMatrix<Tvert> cont(DT.l_dels[i_del].EXT);
    auto get_ineq = [&](int const &i_adj) -> MyVector<T> {
      Delaunay_AdjO<Tvert> adj = DT.l_dels[i_del].ListAdj[i_adj];
      int j_del = adj.iOrb;
      MyMatrix<Tvert> EXTadj = DT.l_dels[j_del].EXT * adj.eBigMat;
      int len = EXTadj.rows();
      for (int u = 0; u < len; u++) {
        MyVector<Tvert> TheVert = GetMatrixRow(EXTadj, u);
        std::optional<size_t> opt = cont.GetIdx_v(TheVert);
        if (!opt) {
          return VoronoiLinearInequality(vipc, TheVert, ListGram, os);
        }
      }
      std::cerr << "Failed to find a matching entry in IsoDelaunay::get_ineq\n";
      throw TerminalException{1};
    };
    for (int i_adj = 0; i_adj < n_adj; i_adj++) {
      MyVector<T> V = get_ineq(i_adj);
      // That canonicalization is incorrect because it is not invariant under
      // the group of transformation. The right group transformations are the
      // ones that The list of matrices ListGram must be an integral basis of
      // the T-space. This forces the transformations to be integral in that
      // basis and so the canonicalization by the integer
      MyVector<T> V_red = ScalarCanonicalizationVector(V);
      AdjInfo eAdj{i_del, i_adj};
      map[V_red].push_back(eAdj);
    }
  }
  std::vector<FullAdjInfo<T>> l_ret;
  for (auto &kv : map) {
    FullAdjInfo<T> fai{kv.first, kv.second};
    l_ret.push_back(fai);
  }
  return l_ret;
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

template <typename T, typename Tvert, typename Tgroup>
std::pair<int, MyMatrix<T>>
GetInteriorGramMatrix(LinSpaceMatrix<T> const &LinSpa,
                      DelaunayTesselation<Tvert, Tgroup> const &DT,
                      std::ostream &os) {
  int n = LinSpa.n;
  int dimSpace = LinSpa.ListMat.size();
  std::vector<FullAdjInfo<T>> ListIneq =
      ComputeDefiningIneqIsoDelaunayDomain<T, Tvert, Tgroup>(
          DT, LinSpa.ListLineMat, os);
  MyMatrix<T> FAC = GetFACineq(ListIneq);
  int nbIneq = FAC.rows();
  MyVector<T> ThePt = GetGeometricallyUniqueInteriorPoint(FAC, os);
  MyMatrix<T> RetMat = ZeroMatrix<T>(n, n);
  for (int u = 0; u < dimSpace; u++) {
    RetMat += ThePt(u) * LinSpa.ListMat[u];
  }
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN
  MyMatrix<T> EXT = DirectFacetComputationInequalities(FAC, "lrs", os);
  int n_row = EXT.rows();
  MyMatrix<T> SumMatExtRay = ZeroMatrix<T>(n, n);
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyMatrix<T> RayMat = ZeroMatrix<T>(n, n);
    for (int u = 0; u < dimSpace; u++) {
      RayMat += EXT(i_row, u) * LinSpa.ListMat[u];
    }
    SumMatExtRay += RemoveFractionMatrix(RayMat);
    DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(RayMat);
    if (DiagInfo.nbMinus > 0) {
      std::cerr
          << "One ray of the IsoDelaunay domain is not positive semidefinite\n";
      std::cerr << "This is not allowed\n";
      throw TerminalException{1};
    }
  }
  RetMat = SumMatExtRay;
#endif
  return {nbIneq, std::move(RetMat)};
}

template <typename T, typename Tint, typename Tgroup>
DelaunayTesselation<Tint, Tgroup> GetInitialGenericDelaunayTesselation(
    DataIsoDelaunayDomains<T, Tint, Tgroup> &data) {
  std::ostream &os = data.rddo.os;
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: GetInitialGenericDelaunayTesselation, beginning\n";
#endif
  auto f_incorrect = [&](Delaunay_Obj<Tint, Tgroup> const &x) -> bool {
    MyMatrix<Tint> const &EXT = x.EXT;
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: Before IsDelaunayPolytopeInducingEqualities\n";
#endif
    bool test1 = IsDelaunayPolytopeInducingEqualities(EXT, data.LinSpa, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: After IsDelaunayPolytopeInducingEqualities test1=" << test1
       << "\n";
#endif
    if (test1) {
      return true;
    }
    if (data.CommonGramMat) {
      MyMatrix<T> const &TestGram = *data.CommonGramMat;
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: Before IsDelaunayAcceptableForGramMat\n";
#endif
      bool test2 =
          IsDelaunayAcceptableForGramMat(EXT, data.LinSpa, TestGram, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: After IsDelaunayAcceptableForGramMat test2=" << test2
         << "\n";
#endif
      if (!test2) {
        return true;
      }
    }
    return false;
  };
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: GetInitialGenericDelaunayTesselation, we have f_incorrect\n";
#endif
  auto test_matrix = [&](MyMatrix<T> const &GramMat)
      -> std::optional<DelaunayTesselation<Tint, Tgroup>> {
    bool test = IsSymmetryGroupCorrect<T, Tint, Tgroup>(GramMat, data.LinSpa, os);
    if (!test) {
      return {};
    }
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
    return EnumerationDelaunayPolytopes<T, Tint, Tgroup, decltype(f_incorrect)>(
        data_lattice, f_incorrect, max_runtime_second);
  };
  while (true) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: Before GetRandomPositiveDefiniteNoNontrivialSymm\n";
#endif
    MyMatrix<T> GramMat =
      GetRandomPositiveDefiniteNoNontrivialSymm<T, Tint, Tgroup>(data.LinSpa, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: After GetRandomPositiveDefiniteNoNontrivialSymm\n";
#endif
    std::optional<DelaunayTesselation<Tint, Tgroup>> opt = test_matrix(GramMat);
    if (opt) {
      return *opt;
    }
  }
  std::cerr << "ISODEL: Failed to find a matching entry in "
               "GetInitialGenericDelaunayTesselation\n";
  throw TerminalException{1};
}

template <typename Tvert, typename Tgroup> struct RepartEntry {
  MyMatrix<Tvert> EXT;
  Tgroup TheStab;
  int8_t Position; // -1: lower, 0: barrel, 1: higher
  int iDelaunayOrigin;
  std::vector<Delaunay_AdjO<Tvert>> ListAdj;
  MyMatrix<Tvert> eBigMat;
};

template <typename Tvert>
std::vector<MyVector<Tvert>>
Orbit_MatrixGroup(std::vector<MyMatrix<Tvert>> const &ListGen,
                  MyVector<Tvert> const &eV, std::ostream &os) {
  auto f_prod = [](MyVector<Tvert> const &v,
                   MyMatrix<Tvert> const &M) -> MyVector<Tvert> {
    return M.transpose() * v;
  };
  return OrbitComputation<MyMatrix<Tvert>, MyVector<Tvert>, decltype(f_prod)>(
      ListGen, eV, f_prod, os);
}

template <typename T, typename Tvert, typename Tgroup> struct FullRepart {
  std::vector<RepartEntry<Tvert, Tgroup>> cells;
  MyVector<T> eIso;
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
template <typename T, typename Tvert, typename Tgroup>
FullRepart<T, Tvert, Tgroup> FindRepartitionningInfoNextGeneration(
    size_t eIdx, DelaunayTesselation<Tvert, Tgroup> const &ListOrbitDelaunay,
    std::vector<AdjInfo> const &ListInformationsOneFlipping,
    MyMatrix<T> const &InteriorElement,
    RecordDualDescOperation<T, Tgroup> &rddo) {
  std::ostream &os = rddo.os;
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FindRepartitionningInfoNextGeneration, begin FRING\n";
  os << "ISODEL: |ListInformationsOneFlipping|=" << ListInformationsOneFlipping.size()
     << "\n";
  for (auto &eEnt : ListInformationsOneFlipping) {
    int iOrbAdj = ListOrbitDelaunay.l_dels[eEnt.iOrb].ListAdj[eEnt.i_adj].iOrb;
    os << "ISODEL:    iOrb=" << eEnt.iOrb << " i_adj=" << eEnt.i_adj
       << " iOrbAdj=" << iOrbAdj << "\n";
  }
#endif
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n = InteriorElement.rows();
  std::vector<std::vector<Tidx>> ListPermGenList;
  std::vector<MyMatrix<Tvert>> ListMatGens;
  std::vector<MyVector<Tvert>> ListVertices;
  std::unordered_map<MyVector<Tvert>, size_t> ListVertices_rev;
  size_t idx_vertices = 1;
  Tgroup PermGRP;
  auto StandardGroupUpdate = [&]() -> void {
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
  };
  StandardGroupUpdate();
  struct TypeOrbitCenter {
    int iDelaunay;
    MyMatrix<Tvert> eBigMat;
    std::vector<Tidx> Linc;
    MyMatrix<Tvert> EXT;
  };
  struct TypeOrbitCenterMin {
    int iDelaunay;
    MyMatrix<Tvert> eBigMat;
  };
  struct FacetEntryEquiv {
    int iOrb;
    MyMatrix<Tvert> eBigMat;
  };
  std::vector<TypeOrbitCenter> ListOrbitCenter;
  auto FuncInsertVertex = [&](MyVector<Tvert> const &eVert) -> void {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertVertex, step 1\n";
#endif
    if (ListVertices_rev.count(eVert) == 1) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: is_present eVert=" << StringVector(eVert) << "\n";
#endif
      return;
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertVertex, step 2\n";
#endif
    std::vector<MyVector<Tvert>> O = Orbit_MatrixGroup(ListMatGens, eVert, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertVertex, step 3 |O|=" << O.size() << "\n";
#endif
    for (auto &eV : O) {
      ListVertices.push_back(eV);
      ListVertices_rev[eV] = idx_vertices;
      idx_vertices++;
    }
    size_t n_gen = ListPermGenList.size();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertVertex, step 4, n_gen=" << n_gen << "\n";
#endif
    for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
      std::vector<Tidx> &ePermGen = ListPermGenList[i_gen];
      MyMatrix<Tvert> eMatrGen = ListMatGens[i_gen];
      for (auto &eVert : O) {
        MyVector<Tvert> eVertImg = eMatrGen.transpose() * eVert;
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
      [&](MyMatrix<Tvert> const &eMat) -> std::optional<std::vector<Tidx>> {
    size_t len = idx_vertices - 1;
    std::vector<Tidx> v_ret(len);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: get_v_positions, len=" << len << "\n";
#endif
    for (size_t i = 0; i < len; i++) {
      MyVector<Tvert> eVimg = eMat.transpose() * ListVertices[i];
      if (ListVertices_rev.count(eVimg) == 0) {
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
  auto FuncInsertGenerator = [&](MyMatrix<Tvert> const &eMat) -> void {
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
        return PermGRP.isin(elt);
      } else {
        std::vector<MyMatrix<Tvert>> LGen = ListMatGens;
        LGen.push_back(eMat);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: FRING: opt, case non-exist\n";
        os << "ISODEL: FRING: Before |ListVertices|=" << ListVertices.size()
           << "\n";
#endif
        size_t nVert = ListVertices.size();
        Face f_att(nVert);
        size_t miss_val = std::numeric_limits<size_t>::max();
        auto get_idx = [&](MyVector<Tvert> const &eVert) -> size_t {
          if (ListVertices_rev.count(eVert) == 0) {
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
            MyVector<Tvert> const &eVert = ListVertices[pos];
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
            os << "ISODEL: FRING: Treating one vertex pos=" << pos << "\n";
#endif
            std::vector<MyVector<Tvert>> TheOrb =
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
      StandardGroupUpdate();
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertGenerator, step 4\n";
#endif
  };
  auto FuncInsertCenter = [&](TypeOrbitCenterMin const &TheRec) -> void {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertCenter, step 1\n";
#endif
    MyMatrix<Tvert> LVert =
        ListOrbitDelaunay.l_dels[TheRec.iDelaunay].EXT * TheRec.eBigMat;
    for (int u = 0; u < LVert.rows(); u++) {
      MyVector<Tvert> eVert = GetMatrixRow(LVert, u);
      FuncInsertVertex(eVert);
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertCenter, step 2\n";
#endif
    StandardGroupUpdate();
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
      MyMatrix<Tvert> eGen =
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
      MyMatrix<Tvert> const &BigMatR = TheRec.eBigMat;
      MyMatrix<Tvert> BigMatI = Inverse(BigMatR);
      MyMatrix<Tvert> const &EXT =
          ListOrbitDelaunay.l_dels[TheRec.iDelaunay].EXT;
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: FuncInsertCenter, step 6, B\n";
#endif
      for (auto &ePermGen : ListOrbitDelaunay.l_dels[TheRec.iDelaunay]
                                .GRP.SmallGeneratingSet()) {
        MyMatrix<Tvert> eBigMat =
            RepresentVertexPermutation(EXT, EXT, ePermGen);
        MyMatrix<Tvert> eBigMat_new = BigMatI * eBigMat * BigMatR;
        FuncInsertGenerator(eBigMat_new);
      }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: FuncInsertCenter, step 6, C\n";
#endif
      std::vector<Tidx> Linc;
      for (int u = 0; u < LVert.rows(); u++) {
        MyVector<Tvert> eV = GetMatrixRow(LVert, u);
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
  TypeOrbitCenterMin TheRec{static_cast<int>(eIdx), IdentityMat<Tvert>(n + 1)};
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
          MyMatrix<Tvert> const &eBigMat =
              ListOrbitDelaunay.l_dels[eCase.iOrb].ListAdj[eCase.i_adj].eBigMat;
          int iOrbAdj =
              ListOrbitDelaunay.l_dels[eCase.iOrb].ListAdj[eCase.i_adj].iOrb;
          MyMatrix<Tvert> eBigMatNew = eBigMat * eEnt.eBigMat;
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
  // second part, the convex decomposition
  int nVert = ListVertices.size();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FRING: nVert=" << nVert << "\n";
#endif
  MyMatrix<T> TotalListVertices(nVert, n + 2);
  MyMatrix<Tvert> TotalListVerticesRed(nVert, n + 1);
  std::vector<T> LineInterior = GetLineVector(InteriorElement);
  MyVector<T> eV(n);
  MyVector<T> eIso(n + 1);
  eIso(0) = 1;
  for (int iVert = 0; iVert < nVert; iVert++) {
    MyVector<Tvert> const &eVert = ListVertices[iVert];
    TotalListVertices(iVert, 0) = 1;
    TotalListVerticesRed(iVert, 0) = 1;
    for (int iCol = 0; iCol < n; iCol++) {
      TotalListVerticesRed(iVert, iCol + 1) = eVert(iCol + 1);
      T val = UniversalScalarConversion<T, Tvert>(eVert(iCol + 1));
      eIso(iCol + 1) += val;
      TotalListVertices(iVert, iCol + 1) = val;
      eV(iCol) = val;
    }
    T Height = EvaluateLineVector(LineInterior, eV);
    TotalListVertices(iVert, n + 1) = Height;
  }
  for (int iCol = 0; iCol < n; iCol++) {
    eIso(iCol + 1) /= nVert;
  }
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
  MyMatrix<T> TotalListVerticesRed_T =
      UniversalMatrixConversion<T, Tvert>(TotalListVerticesRed);
  os << "ISODEL: FRING: Before testing symmetry for TotalListVerticesRed_T\n";
  if (!IsSymmetryGroupOfPolytope(TotalListVerticesRed_T, PermGRP)) {
    std::cerr << "ISODEL: FRING: PermGRP is not a symmetry group of "
                 "TotalListVertices\n";
    throw TerminalException{1};
  }
#endif
  auto get_incd_status = [&](int iVert, MyVector<T> const &eFac) -> bool {
    T eSum = 0;
    for (int u = 0; u <= n + 1; u++) {
      eSum += TotalListVertices(iVert, u) * eFac(u);
    }
    return eSum == 0;
  };
  // The Linc is not ordered while the Linc_face by virtue of being built as a
  // Face has an ordered intrinsic to it.
  // We need this design for the following reason:
  // * If we reorder the vertices of EXT to match the ListVertices, then the
  // incidence are broken or need to be recomputed.
  // * Converting from one ordering to another is not too expensive overall.
  struct RepartEntryProv {
    MyVector<T> eFac;
    std::vector<Tidx> Linc;
    Face Linc_face;
  };
  std::vector<RepartEntry<Tvert, Tgroup>> ListOrbitFacet;
  std::vector<RepartEntryProv> ListOrbitFacet_prov;
  for (auto &eRec : ListOrbitCenter) {
    Face Linc_face = VectorToFace(eRec.Linc, nVert);
    MyVector<T> eFac = FindFacetInequality(TotalListVertices, Linc_face);
    Tgroup TheStab;
    int8_t Position = -1;
    std::vector<Delaunay_AdjO<Tvert>> ListAdj;
    int iDelaunayOrigin = eRec.iDelaunay;
    MyMatrix<Tvert> const &eBigMat = eRec.eBigMat;
    MyMatrix<Tvert> const &EXT = eRec.EXT;
    RepartEntry<Tvert, Tgroup> re{EXT,     TheStab, Position, iDelaunayOrigin,
                                  ListAdj, eBigMat};
    RepartEntryProv rep{eFac, eRec.Linc, Linc_face};
    ListOrbitFacet.push_back(re);
    ListOrbitFacet_prov.push_back(rep);
  }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FRING: we have ListOrbitFacet\n";
#endif
  auto FuncInsertFacet = [&](MyVector<T> const &eFac) -> Delaunay_AdjO<Tvert> {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertFacet, step 1\n";
#endif
    std::vector<Tidx> Linc;
    Face Linc_face(nVert);
    std::vector<MyVector<Tvert>> EXT_list;
    for (int iVert = 0; iVert < nVert; iVert++) {
      if (get_incd_status(iVert, eFac)) {
        Linc.push_back(iVert);
        Linc_face[iVert] = 1;
        EXT_list.push_back(ListVertices[iVert]);
      }
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertFacet, step 2 |EXT_list|="
       << EXT_list.size() << "\n";
#endif
    MyMatrix<Tvert> EXT = MatrixFromVectorFamily(EXT_list);
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
      std::optional<Telt> opt = PermGRP.RepresentativeAction_OnSets(
          ListOrbitFacet_prov[iOrb].Linc_face, Linc_face);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: FuncInsertFacet, after "
            "RepresentativeAction_OnSets\n";
#endif
      if (opt) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: FRING: FuncInsertFacet, finding an equivalence\n";
#endif
        MyMatrix<Tvert> eBigMat = RepresentVertexPermutation(
            TotalListVerticesRed, TotalListVerticesRed, *opt);
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
    std::vector<Delaunay_AdjO<Tvert>> ListAdj;
    MyMatrix<Tvert> eMatUnused; // That matrix should never be used
    RepartEntry<Tvert, Tgroup> re{EXT,     TheStab,   Position, iDelaunayOrigin,
                                  ListAdj, eMatUnused};
    RepartEntryProv rep{eFac, Linc, Linc_face};
    ListOrbitFacet.push_back(re);
    ListOrbitFacet_prov.push_back(rep);
    int iOrb = nOrb;
    MyMatrix<Tvert> eBigMat = IdentityMat<Tvert>(n + 1);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FRING: FuncInsertFacet, step 5\n";
#endif
    return {{}, std::move(eBigMat), iOrb};
  };
  using Text_int = typename SubsetRankOneSolver<T>::Tint;
  MyMatrix<Text_int> TotalListVertices_int = Get_EXT_int(TotalListVertices);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FRING: we have TotalListVertices_int\n";
#endif
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
      Tgroup Stab = PermGRP.Stabilizer_OnSets(Linc_face);
      Tgroup TheStabFace = ReducedGroupActionFace(Stab, Linc_face);
      Tgroup TheStabVect = ReducedGroupActionVect(Stab, Linc);
      ListOrbitFacet[iOrb].TheStab = TheStabVect;
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: FRING: iOrb=" << iOrb << "\n";
      os << "ISODEL: FRING: |Linc_face|=" << Linc_face.size() << " / "
         << Linc_face.count() << "\n";
      os << "ISODEL: FRING: Position="
         << static_cast<int>(ListOrbitFacet[iOrb].Position) << "\n";
      MyMatrix<T> EXT_T =
          UniversalMatrixConversion<T, Tvert>(ListOrbitFacet[iOrb].EXT);
      os << "ISODEL: FRING: Before IsSymmetryGroupOfPolytope after "
            "computation\n";
      if (!IsSymmetryGroupOfPolytope(EXT_T, TheStabVect)) {
        std::cerr << "ISODEL: FRING: The initial computation went wrong\n";
        throw TerminalException{1};
      }
#endif
      std::vector<Delaunay_AdjO<Tvert>> ListAdj;
      // Depending on the nature of the facet (low, barrel, top), the idx_drop
      // can very much vary. We cannot set it to
      int idx_drop = get_idx_drop(ListOrbitFacet_prov[iOrb].eFac);
      MyMatrix<T> EXT2 =
          SelectRowDropColumnFace(TotalListVertices, Linc_face, idx_drop);
      vectface vf = DualDescriptionRecordFullDim(EXT2, TheStabFace, rddo);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: Second while |vf|=" << vf.size()
         << " |TheStabFace|=" << TheStabFace.size() << " |EXT2|=" << EXT2.rows()
         << "/" << EXT2.cols() << " rnk=" << RankMat(EXT2) << "\n";
      CheckFacetInequality(TotalListVertices, Linc_face,
                           "FuncInsertFace TotalListVertices Linc_face");
#endif
      FlippingFramework<T> frame(TotalListVertices, TotalListVertices_int,
                                 Linc_face, os);
      for (auto &eFace : vf) {
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN
        CheckFacetInequality(EXT2, eFace, "FuncInsertFace EXT2 eFace");
#endif
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: Before FlipFace |EXT2|=" << EXT2.rows() << " / "
           << EXT2.cols() << " |eFace|=" << eFace.size() << " / "
           << eFace.count() << "\n";
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
        MyVector<T> eFac = FindFacetInequality(TotalListVertices, eInc);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: FRING: We have eFac\n";
#endif
        Delaunay_AdjO<Tvert> eAdj = FuncInsertFacet(eFac);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: FRING: We have eAdj\n";
#endif
        size_t nVert = ListOrbitFacet_prov[iOrb].Linc.size();
        Face LEV(nVert);
        for (size_t iInc = 0; iInc < nVert; iInc++) {
          Tidx jInc = ListOrbitFacet_prov[iOrb].Linc[iInc];
          if (get_incd_status(jInc, eFac)) {
            LEV[iInc] = 1;
          }
        }
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
  return {ListOrbitFacet, eIso};
}

/*
  Effectively do the flipping. Note that the whole computation is combinatorial
  and all the polyhedral computations are done in the
  FindRepartitionningInfoNextGeneration code.
 */
template <typename T, typename Tvert, typename Tgroup>
DelaunayTesselation<Tvert, Tgroup>
FlippingLtype(DelaunayTesselation<Tvert, Tgroup> const &ListOrbitDelaunay,
              MyMatrix<T> const &InteriorElement,
              std::vector<AdjInfo> const &ListInformationsOneFlipping,
              RecordDualDescOperation<T, Tgroup> &rddo) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  std::ostream &os = rddo.os;
  os << "ISODEL: FLT: FlippingLtype, begin\n";
#endif
  using Tgr = GraphListAdj;
  using Telt = typename Tgroup::Telt;
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
  std::vector<std::vector<RepartEntry<Tvert, Tgroup>>> ListInfo;
  std::vector<MyVector<T>> ListIso;
  std::vector<int> vect_iInfo(n_dels, -1);
  std::vector<int> vect_lower_iFacet(n_dels, -1);
  int iInfo = 0;
  for (auto &eConn : ListGroupMelt) {
    size_t eIdx = eConn[0];
    FullRepart<T, Tvert, Tgroup> fr = FindRepartitionningInfoNextGeneration(
        eIdx, ListOrbitDelaunay, ListInformationsOneFlipping, InteriorElement,
        rddo);
    int n_facet = fr.cells.size();
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: FLT: iInfo=" << iInfo << " |eConn|=" << eConn.size()
       << " n_facet=" << n_facet << "\n";
    size_t n_zero = 0, n_plus = 0, n_minus = 0;
#endif
    for (int iFacet = 0; iFacet < n_facet; iFacet++) {
      RepartEntry<Tvert, Tgroup> const &eFacet = fr.cells[iFacet];
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
    Delaunay_AdjO<Tvert> adj;
    MyMatrix<Tvert> eBigMat;
  };
  auto get_matching_listinfo = [&](int const &iInfo, int const &iFacet,
                                   Face const &eInc) -> MatchedFacet {
    MyMatrix<Tvert> const &EXT = ListInfo[iInfo][iFacet].EXT;
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
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tvert>(EXT);
    if (!IsSymmetryGroupOfPolytope(EXT_T, TheStab)) {
      std::cerr << "The group TheStab is not a symmetry group\n";
      throw TerminalException{1};
    }
#endif
    auto get_bigmat = [&](Telt const &ePerm) -> MyMatrix<Tvert> {
      int8_t Position = ListInfo[iInfo][iFacet].Position;
      if (Position == 0) {
        // For barrel case, we need to extend by one point in order to get
        // a full dimensional set. That point has to be preserved by the
        // corresponding group and so we take the isobarycenter of the
        // repartitioning polytope.
        int nbRow = EXT.rows();
        MyMatrix<T> EXT_ext(nbRow + 1, dim + 1);
        for (int iRow = 0; iRow < nbRow; iRow++) {
          EXT_ext(iRow, 0) = 1;
          for (int i = 0; i < dim; i++) {
            T val = UniversalScalarConversion<T, Tvert>(EXT(iRow, i + 1));
            EXT_ext(iRow, i + 1) = val;
          }
        }
        for (int i = 0; i <= dim; i++) {
          EXT_ext(nbRow, i) = ListIso[iInfo](i);
        }
        auto f = [&](int const &u) -> int {
          if (u == nbRow) {
            return nbRow;
          }
          return ePerm.at(u);
        };
        MyMatrix<T> BigMat_T = FindTransformation_f(EXT_ext, EXT_ext, f);
        return UniversalMatrixConversion<Tvert, T>(BigMat_T);
      } else {
        return RepresentVertexPermutation(EXT, EXT, ePerm);
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
      std::optional<Telt> opt =
          TheStab.RepresentativeAction_OnSets(eAdj.eInc, eInc);
      if (opt) {
        MyMatrix<Tvert> eBigMat = get_bigmat(*opt);
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
    MyMatrix<Tvert> const &EXT = ListOrbitDelaunay.l_dels[iDelaunayOld].EXT;
#ifdef SANITY_CHECK_ISO_DELAUNAY_DOMAIN
    CheckFacetInequality(EXT, eInc, "get_matching_old_tessel EXT eInc");
#endif
    for (auto &eAdj : ListOrbitDelaunay.l_dels[iDelaunayOld].ListAdj) {
#ifdef PRINT_ISO_DELAUNAY_DOMAIN_REPRESENTATIVE_ACTION
      PrintRepresentativeAction_OnSets_GRP_f1_f2(ListOrbitDelaunay.l_dels[iDelaunayOld].GRP, eAdj.eInc, eInc);
#endif
      std::optional<Telt> opt =
          ListOrbitDelaunay.l_dels[iDelaunayOld]
              .GRP.RepresentativeAction_OnSets(eAdj.eInc, eInc);
      if (opt) {
        MyMatrix<Tvert> eBigMat = RepresentVertexPermutation(EXT, EXT, *opt);
        return {eAdj, eBigMat};
      }
    }
    std::cerr << "ISO_EDGE: FLT: iDelaunayOld=" << iDelaunayOld
              << " |eInc|=" << eInc.size() << " / " << eInc.count() << "\n";
    std::cerr << "ISO_EDGE: FLT: Failed to find a matching entry in "
                 "get_matching_old_tessel\n";
    throw TerminalException{1};
  };
  auto get_lower_adjacency = [&](int const &iInfo,
                                 int const &iFacet) -> Delaunay_AdjO<Tvert> {
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
  auto get_face_m_m = [](MyMatrix<Tvert> const &M1,
                         MyMatrix<Tvert> const &M2) -> Face {
    int nVert1 = M1.rows();
    int nVert2 = M2.rows();
    Face f2(nVert2);
    ContainerMatrix<Tvert> cont(M2);
    for (int i1 = 0; i1 < nVert1; i1++) {
      MyVector<Tvert> V = GetMatrixRow(M1, i1);
      std::optional<size_t> opt = cont.GetIdx_v(V);
      size_t pos2 = unfold_opt(opt, "Error in get_face_m_m");
      f2[pos2] = 1;
    }
    return f2;
  };
  auto get_face_msub_m = [](MyMatrix<Tvert> const &M1, Face const &f1,
                            MyMatrix<Tvert> const &M2) -> Face {
    int nVert1 = M1.rows();
    int nVert2 = M2.rows();
    Face f2(nVert2);
    ContainerMatrix<Tvert> cont(M2);
    for (int i1 = 0; i1 < nVert1; i1++) {
      if (f1[i1] == 1) {
        MyVector<Tvert> V = GetMatrixRow(M1, i1);
        std::optional<size_t> opt = cont.GetIdx_v(V);
        int pos2 = unfold_opt(opt, "Error in get_face_msub_m");
        f2[pos2] = 1;
      }
    }
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
  std::vector<Delaunay_Entry<Tvert, Tgroup>> l_dels;
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
    Delaunay_Entry<Tvert, Tgroup> del{ListOrbitDelaunay.l_dels[iDelaunay].EXT,
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
        MyMatrix<Tvert> EXT = ListInfo[iInfo][iFacet].EXT;
        Tgroup GRP = ListInfo[iInfo][iFacet].TheStab;
        Delaunay_Entry<Tvert, Tgroup> del{EXT, GRP, {}};
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
  auto check_adj = [&](int iOrb, Delaunay_AdjO<Tvert> const &NewAdj,
                       std::string const &context) -> void {
    MyMatrix<Tvert> const &EXT = l_dels[iOrb].EXT;
    ContainerMatrix<Tvert> cont(EXT);
    Face f_att(EXT.rows());
    MyMatrix<Tvert> EXTadj = l_dels[NewAdj.iOrb].EXT * NewAdj.eBigMat;
    os << "ISODEL: FLT: check_adj iOrb=" << iOrb
       << " NewAdj.iOrb=" << NewAdj.iOrb << "\n";
    os << "ISODEL: FLT:      |EXT|=" << EXT.rows() << " / " << EXT.cols()
       << "\n";
    os << "ISODEL: FLT:   |EXTadj|=" << EXTadj.rows() << " / " << EXTadj.cols()
       << "\n";
    int len = EXTadj.rows();
    for (int iVert = 0; iVert < len; iVert++) {
      MyVector<Tvert> V = GetMatrixRow(EXTadj, iVert);
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
  for (int iOrb = 0; iOrb < n_del_ret; iOrb++) {
    DelaunaySymb ds = NewListOrbitDelaunay[iOrb];
    std::vector<Delaunay_AdjO<Tvert>> ListAdj;
    if (ds.Position == Position_old) {
      int iDelaunay = ds.iDelaunay;
      for (auto &eAdj : ListOrbitDelaunay.l_dels[iDelaunay].ListAdj) {
        int iDelaunayOld = eAdj.iOrb;
        DelaunaySymb dss{Position_old, iDelaunayOld, -1, -1};
        std::optional<size_t> opt = get_symbol_position(dss);
        if (opt) {
          int iOrbAdj = *opt;
          Delaunay_AdjO<Tvert> NAdj{eAdj.eInc, eAdj.eBigMat, iOrbAdj};
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
          RepartEntry<Tvert, Tgroup> const &eFacet = ListInfo[iInfo][iFacet];
          MyMatrix<Tvert> const &BigMat2 = eFacet.eBigMat;
          MyMatrix<Tvert> ImageEXT =
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
          MyMatrix<Tvert> BigMat1 = RecMatch.adj.eBigMat * RecMatch.eBigMat *
                                    Inverse(BigMat2) * eAdj.eBigMat;
          Delaunay_AdjO<Tvert> NAdj{eAdj.eInc, BigMat1, Pos};
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
          MyMatrix<Tvert> const &eMat1 = eAdj.eBigMat;
          MyMatrix<Tvert> LincEXT =
              SelectRow(ListInfo[iInfo][iFacet].EXT, eAdj.eInc);
          MyMatrix<Tvert> ImageEXTbarrel = ListInfo[iInfo][jFacet].EXT * eMat1;
          Face LLinc = get_face_m_m(LincEXT, ImageEXTbarrel);
          Delaunay_AdjO<Tvert> TheFoundAdj = get_lower_adjacency(iInfo, jFacet);
          int kFacet = TheFoundAdj.iOrb;
          MyMatrix<Tvert> const &eMat2 = TheFoundAdj.eBigMat;
          int iDelaunayOrigin = ListInfo[iInfo][kFacet].iDelaunayOrigin;
          MyMatrix<Tvert> ImageEXT =
              ListInfo[iInfo][kFacet].EXT * TheFoundAdj.eBigMat;
          Face LLinc2 = get_face_msub_m(ListInfo[iInfo][jFacet].EXT,
                                        TheFoundAdj.eInc, ImageEXT);
          MatchedFacet match2 =
              get_matching_old_tessel(iDelaunayOrigin, LLinc2);
          Delaunay_AdjO<Tvert> const &TheFoundAdj2 = match2.adj;
          MyMatrix<Tvert> const &TheMat2 = match2.eBigMat;
          int jDelaunayOrigin = TheFoundAdj2.iOrb;
          MyMatrix<Tvert> const &eMat3 = TheFoundAdj2.eBigMat;
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
          Delaunay_AdjO<Tvert> const &TheFoundAdj3 = match3.adj;
          MyMatrix<Tvert> const &TheMat3 = match3.eBigMat;
          MyMatrix<Tvert> BigMat1 =
              TheFoundAdj3.eBigMat * TheMat3 *
              Inverse(ListInfo[jInfo][iFacet2].eBigMat) * eMat3 * TheMat2 *
              ListInfo[iInfo][kFacet].eBigMat * eMat2 * eMat1;
          MyMatrix<Tvert> EXT7 =
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
          Delaunay_AdjO<Tvert> const &TheFoundAdj4 = match4.adj;
          MyMatrix<Tvert> const &TheMat4 = match4.eBigMat;
          DelaunaySymb dss{Position_new, -1, jInfo, TheFoundAdj4.iOrb};
          std::optional<size_t> optN = get_symbol_position(dss);
          int Pos = unfold_opt(optN, "Failed to find entry for Case 3");
          MyMatrix<Tvert> BigMat2 = TheFoundAdj4.eBigMat * TheMat4 * BigMat1;
          Delaunay_AdjO<Tvert> NAdj{eAdj.eInc, BigMat2, Pos};
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
          check_adj(iOrb, NAdj, "Case 3");
#endif
          ListAdj.push_back(NAdj);
        }
        if (ListInfo[iInfo][jFacet].Position == -1) {
          int iDelaunayOrigin = ListInfo[iInfo][jFacet].iDelaunayOrigin;
          MyMatrix<Tvert> const &eMat1 = eAdj.eBigMat;
          MyMatrix<Tvert> ImageEXT = ListInfo[iInfo][jFacet].EXT * eMat1;
          Face LLinc =
              get_face_msub_m(ListInfo[iInfo][iFacet].EXT, eAdj.eInc, ImageEXT);
          MatchedFacet match1 = get_matching_old_tessel(iDelaunayOrigin, LLinc);
          Delaunay_AdjO<Tvert> const &TheFoundAdj1 = match1.adj;
          MyMatrix<Tvert> const &TheMat1 = match1.eBigMat;
          int jDelaunayOld = TheFoundAdj1.iOrb;
          if (vect_iInfo[jDelaunayOld] == -1) {
            DelaunaySymb dss{Position_old, jDelaunayOld, -1, -1};
            std::optional<size_t> opt = get_symbol_position(dss);
            int Pos2 = unfold_opt(opt, "Case 4");
            MyMatrix<Tvert> BigMat1 = TheFoundAdj1.eBigMat * TheMat1 *
                                      ListInfo[iInfo][jFacet].eBigMat *
                                      eAdj.eBigMat;
            Delaunay_AdjO<Tvert> NAdj{eAdj.eInc, BigMat1, Pos2};
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
            check_adj(iOrb, NAdj, "Case 4");
#endif
            ListAdj.push_back(NAdj);
          } else {
            int jInfo = vect_iInfo[jDelaunayOld];
            int iFacet2 = vect_lower_iFacet[jDelaunayOld];
            MyMatrix<Tvert> ImageEXT =
                ListOrbitDelaunay.l_dels[jDelaunayOld].EXT *
                TheFoundAdj1.eBigMat;
            Face LLinc2 =
                get_face_msub_m(ListOrbitDelaunay.l_dels[iDelaunayOrigin].EXT,
                                TheFoundAdj1.eInc, ImageEXT);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
            os << "ISODEL: FLT: Case 5\n";
#endif
            MatchedFacet match2 = get_matching_listinfo(jInfo, iFacet2, LLinc2);
            Delaunay_AdjO<Tvert> const &TheFoundAdj2 = match2.adj;
            MyMatrix<Tvert> const &TheMat2 = match2.eBigMat;
            DelaunaySymb dss{Position_new, -1, jInfo, TheFoundAdj2.iOrb};
            std::optional<size_t> opt = get_symbol_position(dss);
            int Pos = unfold_opt(opt, "Case 5");
            MyMatrix<Tvert> BigMat1 =
                TheFoundAdj2.eBigMat * TheMat2 *
                Inverse(ListInfo[jInfo][iFacet2].eBigMat) *
                TheFoundAdj1.eBigMat * TheMat1 *
                ListInfo[iInfo][jFacet].eBigMat * eAdj.eBigMat;
            Delaunay_AdjO<Tvert> NAdj{eAdj.eInc, BigMat1, Pos};
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
          Delaunay_AdjO<Tvert> NAdj{eAdj.eInc, eAdj.eBigMat, Pos};
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
          check_adj(iOrb, NAdj, "Case 6");
#endif
          ListAdj.push_back(NAdj);
        }
      }
    }
    l_dels[iOrb].ListAdj = ListAdj;
  }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  os << "ISODEL: FLT: Exiting\n";
#endif
  return {l_dels};
}

FullNamelist NAMELIST_GetStandard_COMPUTE_LATTICE_IsoDelaunayDomains() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["arithmetic_T"] = "gmp_rational";
  ListStringValues1["arithmetic_Tint"] = "gmp_integer";
  ListStringValues1["OutFormat"] = "nothing";
  ListStringValues1["OutFile"] = "unset.out";
  ListStringValues1["FileDualDescription"] = "unset";
  ListStringValues1["CommonGramMat"] = "unset";
  ListIntValues1["max_runtime_second"] = 0;
  ListBoolValues1["ApplyStdUnitbuf"] = false;
  ListBoolValues1["Saving"] = false;
  ListStringValues1["Prefix"] = "/irrelevant/";
  ListStringValues1["CVPmethod"] = "SVexact";
  SingleBlock BlockDATA;
  BlockDATA.ListIntValues = ListIntValues1;
  BlockDATA.ListBoolValues = ListBoolValues1;
  BlockDATA.ListDoubleValues = ListDoubleValues1;
  BlockDATA.ListStringValues = ListStringValues1;
  BlockDATA.ListListStringValues = ListListStringValues1;
  ListBlock["DATA"] = BlockDATA;
  // TSPACE
  ListBlock["TSPACE"] = SINGLEBLOCK_Get_Tspace_Description();
  // Merging all data
  return {ListBlock, "undefined"};
}

template <typename T, typename Tint, typename Tgroup> struct IsoDelaunayDomain {
  DelaunayTesselation<Tint, Tgroup> DT;
  int nbIneq;
  MyMatrix<T> GramMat;
};

template <typename T, typename Tint, typename Tgroup>
void WriteEntryGAP(std::ostream &os_out,
                   IsoDelaunayDomain<T, Tint, Tgroup> const &ent) {
  os_out << "rec(DT:=";
  WriteEntryGAP(os_out, ent.DT);
  os_out << ", nbIneq:=" << ent.nbIneq;
  os_out << ", GramMat:=" << StringMatrixGAP(ent.GramMat) << ")";
}

template <typename T, typename Tint, typename Tgroup>
void WriteEntryPYTHON(std::ostream &os_out,
                      IsoDelaunayDomain<T, Tint, Tgroup> const &ent) {
  os_out << "{\"DT\":";
  WriteEntryPYTHON(os_out, ent.DT);
  os_out << ", \"nbIneq\":" << ent.nbIneq;
  os_out << ", \"GramMat\":" << StringMatrixPYTHON(ent.GramMat) << "}";
}

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, IsoDelaunayDomain<T, Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("DT", eRec.DT);
  ar &make_nvp("nbIneq", eRec.nbIneq);
  ar &make_nvp("GramMat", eRec.GramMat);
}
} // namespace boost::serialization

template <typename T, typename Tint, typename Tgroup>
struct IsoDelaunayDomain_AdjI {
  MyVector<T> V;
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

template <typename T, typename Tint> struct IsoDelaunayDomain_AdjO {
  MyVector<T> V;
  MyMatrix<Tint> eBigMat;
};

template <typename T, typename Tint>
void WriteEntryGAP(std::ostream &os_out,
                   IsoDelaunayDomain_AdjO<T, Tint> const &ent) {
  os_out << "rec(V:=" << StringVectorGAP(ent.V)
         << " eBigMat:=" << StringMatrixGAP(ent.eBigMat) << ")";
}

template <typename T, typename Tint>
void WriteEntryPYTHON(std::ostream &os_out,
                      IsoDelaunayDomain_AdjO<T, Tint> const &ent) {
  os_out << "{\"V\":" << StringVectorPYTHON(ent.V)
         << " \"eBigMat\":" << StringMatrixPYTHON(ent.eBigMat) << "}";
}

namespace boost::serialization {
template <class Archive, typename T, typename Tint>
inline void serialize(Archive &ar, IsoDelaunayDomain_AdjO<T, Tint> &eRec,
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
    size_t const &seed, DelaunayTesselation<Tint, Tgroup> const &DT) {
  using TintGroup = typename Tgroup::Tint;
  std::map<size_t, size_t> map_delaunays;
  auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
    seed ^= new_hash + 0x9e3779b8 + (seed << 6) + (seed >> 2);
  };
  size_t seed_delaunay = 10;
  for (auto &eDel : DT.l_dels) {
    size_t hash = seed_delaunay;
    size_t hash_nVert = eDel.EXT.rows();
    TintGroup order = eDel.GRP.size();
    size_t hash_grp = std::hash<TintGroup>()(order);
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
  DelaunayTesselation<Tint, Tgroup> DT =
      GetInitialGenericDelaunayTesselation(data);
  std::pair<int, MyMatrix<T>> pair =
      GetInteriorGramMatrix(data.LinSpa, DT, data.rddo.os);
  return {std::move(DT), pair.first, std::move(pair.second)};
}

template <typename T, typename Tint, typename Tgroup>
struct IsoDelaunayDomain_Obj {
  IsoDelaunayDomain<T, Tint, Tgroup> DT_gram;
  std::vector<FullAdjInfo<T>> ListIneq;
  Tgroup GRPperm;
};

template <typename T, typename Tint, typename Tgroup>
void WriteBasicEntryGAP(std::ostream &os_out,
                        IsoDelaunayDomain_Obj<T, Tint, Tgroup> const &ent) {
  os_out << "DT_gram:=";
  WriteEntryGAP(os_out, ent.DT_gram);
  //
  os_out << ", ListIneq:=[";
  bool IsFirst = true;
  for (auto &eFullAI : ent.ListIneq) {
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
                           DataIsoDelaunayDomains<T, Tint, Tgroup> const& data,
                           IsoDelaunayDomain_Obj<T, Tint, Tgroup> const &ent, std::ostream& os) {
  os_out << "rec(";
  //  WriteBasicEntryGAP(os_out, ent);
  os_out << "GRPpermSize:=" << ent.GRPperm.size();
  int dimSpace = data.LinSpa.ListMat.size();
  int n = ent.DT_gram.GramMat.rows();
  MyMatrix<T> FAC = GetFACineq(ent.ListIneq);
  MyMatrix<T> FAC_extend = AddFirstZeroColumn(FAC);
  std::vector<int> ListIrred =
    cdd::RedundancyReductionClarkson(FAC_extend, os);
  int nbIneq = FAC.rows();
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(ent.DT_gram.GramMat, os);
  os_out << ", n_ineq:=" << nbIneq;
  os_out << ", n_irred:=" << ListIrred.size();
  os_out << ", det:=" << DeterminantMat(ent.DT_gram.GramMat);
  os_out << ", n_shv:=" << SHV.rows();
  //
  MyMatrix<T> EXT = DirectFacetComputationInequalities(FAC, "lrs", os);
  int n_row = EXT.rows();
  std::map<int, size_t> map_rank;
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyMatrix<T> RayMat = ZeroMatrix<T>(n, n);
    for (int u = 0; u < dimSpace; u++) {
      RayMat += EXT(i_row, u) * data.LinSpa.ListMat[u];
    }
    int rnk = RankMat(RayMat);
    map_rank[rnk] += 1;
  }
  auto write_map_val=[&](std::map<int,size_t> const& map) -> void {
    bool IsFirst = true;
    os_out << "[";
    for (auto & kv: map) {
      if (!IsFirst) {
        os_out << ",";
      }
      IsFirst=false;
      os_out << "[" << kv.first << "," << kv.second << "]";
    }
    os_out << "]";
  };
  os_out << ", ListRank:=";
  write_map_val(map_rank);
  //
  std::map<int, size_t> map_nb_vert;
  for (auto & eDel: ent.DT_gram.DT.l_dels) {
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
  os_out << ", \"ListIneq\":[";
  bool IsFirst = true;
  for (auto &eFullAI : ent.ListIneq) {
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
  ar &make_nvp("ListIneq", eRec.ListIneq);
  ar &make_nvp("GRPperm", eRec.GRPperm);
}
} // namespace boost::serialization

template <typename T, typename Tint, typename Tgroup>
struct DataIsoDelaunayDomainsFunc {
  DataIsoDelaunayDomains<T, Tint, Tgroup> data;
  using Tobj = IsoDelaunayDomain_Obj<T, Tint, Tgroup>;
  using TadjI = IsoDelaunayDomain_AdjI<T, Tint, Tgroup>;
  using TadjO = IsoDelaunayDomain_AdjO<T, Tint>;
  std::ostream &get_os() { return data.rddo.os; }
  Tobj f_init() {
    IsoDelaunayDomain<T, Tint, Tgroup> IsoDel =
        GetInitialIsoDelaunayDomain(data);
    Tobj x{IsoDel, {}, {}};
    return x;
  }
  size_t f_hash(size_t const &seed, Tobj const &x) {
    return ComputeInvariantIsoDelaunayDomain<T, Tint, Tgroup>(data, seed,
                                                              x.DT_gram.DT);
  }
  std::optional<TadjO> f_repr(Tobj const &x, TadjI const &y) {
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    data.rddo.os << "ISODEL: f_repr, before LINSPA_TestEquivalenceGramMatrix\n";
#endif
    std::optional<MyMatrix<Tint>> opt =
        LINSPA_TestEquivalenceGramMatrix<T, Tint, Tgroup>(
            data.LinSpa, x.DT_gram.GramMat, y.DT_gram.GramMat, data.rddo.os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    data.rddo.os << "ISODEL: f_repr, after LINSPA_TestEquivalenceGramMatrix\n";
#endif
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const &eBigMat = *opt;
    TadjO ret{y.V, eBigMat};
    return ret;
  }
  std::pair<Tobj, TadjO> f_spann(TadjI const &x) {
    IsoDelaunayDomain<T, Tint, Tgroup> IsoDel = x.DT_gram;
    Tobj x_ret{IsoDel, {}, {}};
    MyMatrix<Tint> eBigMat = IdentityMat<Tint>(data.LinSpa.n);
    TadjO ret{x.V, eBigMat};
    return {std::move(x_ret), std::move(ret)};
  }
  std::vector<TadjI> f_adj(Tobj &x_in) {
    using Telt = typename Tgroup::Telt;
    using Tidx = typename Telt::Tidx;
    std::ostream &os = data.rddo.os;
    int n = data.LinSpa.n;
    int dimSpace = data.LinSpa.ListMat.size();
    IsoDelaunayDomain<T, Tint, Tgroup> &x = x_in.DT_gram;
    // compute the inequalities
    std::vector<FullAdjInfo<T>> ListIneq =
        ComputeDefiningIneqIsoDelaunayDomain<T, Tint, Tgroup>(
            x.DT, data.LinSpa.ListLineMat, os);
    x_in.ListIneq = ListIneq;
    // Compute the irredundant ones as well as the l_ineq / map_ineq
    MyMatrix<T> FAC = GetFACineq(ListIneq);
    MyMatrix<T> FAC_extend = AddFirstZeroColumn(FAC);
    std::vector<int> ListIrred =
      cdd::RedundancyReductionClarkson(FAC_extend, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    std::vector<int> ListIrred_DD = RedundancyReductionDualDescription(FAC_extend, "lrs", os);
    Face f_irred = VectorToFace(ListIrred, FAC_extend.rows());
    Face f_irred_DD = VectorToFace(ListIrred_DD, FAC_extend.rows());
    if (f_irred != f_irred_DD) {
      std::cerr << "Inconsistency in the computation of redundancy\n";
      throw TerminalException{1};
    } else {
      os << "ISODEL: Coherency of the redundancy computation\n";
    }
#endif
    size_t nbIrred = ListIrred.size();
    MyMatrix<T> FACred = SelectRow(FAC, ListIrred);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: f_adj: |FAC|=" << FAC.rows() << " / " << FAC.cols()
       << " nbIrred=" << nbIrred << "\n";
    os << "ISODEL: x.GramMat=\n";
    WriteMatrix(os, x.GramMat);
#endif
    std::vector<MyVector<T>> l_ineq;
    std::unordered_map<MyVector<T>, size_t> map_ineq;
    for (size_t i = 0; i < nbIrred; i++) {
      MyVector<T> eV = GetMatrixRow(FACred, i);
      l_ineq.push_back(eV);
      map_ineq[eV] = i;
    }
    // Compute the automorphism group on the central gram and then the facets
    std::vector<MyMatrix<T>> ListGenTot =
        LINSPA_ComputeStabilizer<T, Tint, Tgroup>(data.LinSpa, x.GramMat,
                                                  os);
    std::vector<Telt> ListPermGens;
    for (auto &eGenTot : ListGenTot) {
      MyMatrix<T> MatSpace = matrix_in_t_space(eGenTot, data.LinSpa);
      std::vector<Tidx> l_pos(nbIrred);
      for (size_t i = 0; i < nbIrred; i++) {
        MyVector<T> const &eV = l_ineq[i];
        // There are actually two transpose here:
        // * One from going from action in the space to action on the dual
        // * One for going from row to column action
        MyVector<T> eVimg = MatSpace * eV;
        size_t pos = map_ineq.at(eVimg);
        l_pos[i] = pos;
      }
      Telt ePermGen(l_pos);
      ListPermGens.push_back(ePermGen);
    }
    Tgroup GRPperm = Tgroup(ListPermGens, nbIrred);
    x_in.GRPperm = GRPperm;


    std::vector<size_t> l_idx = DecomposeOrbitPoint_FullRepr(GRPperm);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: f_adj: |GRPperm|=" << GRPperm.size()
       << " nbIrred=" << nbIrred << " |l_idx|=" << l_idx.size() << "\n";
#endif
    std::vector<TadjI> l_adj;
    for (auto &i : l_idx) {
      int idxIrred = ListIrred[i];
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: f_adj, i=" << i << " idxIrred=" << idxIrred << "\n";
#endif
      MyVector<T> TestPt = GetSpaceInteriorPointFacet(FACred, i, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: f_adj, We have TestPt\n";
#endif
      MyMatrix<T> TestMat = ZeroMatrix<T>(n, n);
      for (int u = 0; u < dimSpace; u++) {
        TestMat += TestPt(u) * data.LinSpa.ListMat[u];
      }
      bool test = IsPositiveDefinite(TestMat);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
      os << "ISODEL: f_adj, We have test for TestMat\n";
#endif
      if (test) {
        FullAdjInfo<T> eRecIneq = ListIneq[idxIrred];
        DelaunayTesselation<Tint, Tgroup> DTadj =
            FlippingLtype<T, Tint, Tgroup>(x.DT, x.GramMat,
                                           eRecIneq.ListAdjInfo, data.rddo);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: After FlippingLtype\n";
#endif
        std::pair<int, MyMatrix<T>> pair =
            GetInteriorGramMatrix(data.LinSpa, DTadj, os);
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
        os << "ISODEL: After GetInteriorGramMatrix\n";
#endif
        IsoDelaunayDomain<T, Tint, Tgroup> IsoDelAdj{
            std::move(DTadj), pair.first, std::move(pair.second)};
        TadjI eAdj{eRecIneq.eIneq, IsoDelAdj};
        l_adj.push_back(eAdj);
      }
    }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
    os << "ISODEL: Before returning l_adj\n";
#endif
    return l_adj;
  }
  Tobj f_adji_obj(TadjI const &x) { return {x.DT_gram, {}, {}}; }
  size_t f_complexity(Tobj const &x) { return x.nbIneq; }
};

template <typename T, typename Tint, typename Tgroup>
DataIsoDelaunayDomains<T, Tint, Tgroup> get_data_isodelaunay_domains(FullNamelist const &eFull,
                                                                     PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
                                                                     std::ostream& os) {
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  SingleBlock BlockTSPACE = eFull.ListBlock.at("TSPACE");
  auto get_common = [&]() -> std::optional<MyMatrix<T>> {
    std::string CommonGramMat = BlockDATA.ListStringValues.at("CommonGramMat");
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
  DataIsoDelaunayDomains<T, Tint, Tgroup> data{LinSpa, std::move(rddo),
                                               CommonGramMat};
  return data;
}


// clang-format off
#endif  // SRC_LATT_ISODELAUNAYDOMAINS_H_
// clang-format on

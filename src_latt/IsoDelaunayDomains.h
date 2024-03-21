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
// clang-format on

template <typename T, typename Tint, typename Tgroup>
struct DataIsoDelaunayDomains {
  LinSpaceMatrix<T> LinSpa;
  RecordDualDescOperation<T,Tgroup> rddo;
  std::optional<MyMatrix<T>> CommonGramMat;
};

/*
  Code for the L-type domains.

  Two main use case:
  ---Lattice case: Then the Tvert is actually a Tint and can be mpz_class, int32_t, etc.
  ---Periodic structure case: Then the coordinates are no longer integral.
    Also the equivalence are no longer integral. Sure the matrix transformation is
    integral, but the translation vector is not necessarily so.

  We use the definitions from LatticeDelaunay.h
  They are for the lattice case but can be generalized for the periodic case.
 */


template<typename T>
struct VoronoiInequalityPreComput {
  //  MyMatrix<Tvert> VertBasis;
  int n;
  Face f_basis;
  MyMatrix<T> VertBasis_T;
  MyMatrix<T> VertBasisInv_T;
  std::vector<MyVector<T>> VertBasisRed_T;
};

/*
  We implement the hash of a Delaunay tessellationn. The constraint is that two different
  but equivalent tessellations, must have the same hash. Hopefully, this is not a problem
  since we have many invariants:
  * The number of vertices of the orbit representatives of Delaunay polytopes.
  * The size of their automorphism groups
  * The number of vertices of the orbit representative of their facets.
 */
template<typename Tvert, typename Tgroup>
size_t ComputeInvariantDelaunayTessellation(DelaunayTesselation<Tvert,Tgroup> const& DT,
                                            size_t const& seed,
                                            [[maybe_unused]] std::ostream & os) {
  using TintGroup = typename Tgroup::Tint;
  std::map<size_t,size_t> map;
  auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
    seed ^= new_hash + 0x9e3779b8 + (seed << 6) + (seed >> 2);
  };
  for (auto & e_del : DT.l_dels) {
    std::map<size_t,size_t> map_siz;
    for (auto & eAdj : e_del.ListAdj) {
      size_t siz = eAdj.eInc.count();
      map_siz[siz] += 1;
    }
    size_t hash_del = 123;
    size_t hash1 = std::hash<int>()(e_del.EXT.rows());
    size_t hash2 = std::hash<TintGroup>()(e_del.GRP.order());
    combine_hash(hash_del, hash1);
    combine_hash(hash_del, hash2);
    for (auto & kv : map_siz) {
      combine_hash(hash_del, kv.first);
      combine_hash(hash_del, kv.second);
    }
    map[hash_del] += 1;
  }
  size_t hash_ret = seed;
  for (auto & kv : map) {
    combine_hash(hash_ret, kv.first);
    combine_hash(hash_ret, kv.second);
  }
  return hash_ret;
}


template<typename T, typename Tvert>
VoronoiInequalityPreComput<T> BuildVoronoiIneqPreCompute(MyMatrix<Tvert> const& EXT) {
  int n = EXT.cols() - 1;
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T,Tvert>(EXT);
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(EXT_T);
  std::vector<int> ListRowSelect = eSelect.ListRowSelect;
  Face f_basis(EXT_T.rows());
  for (auto & eIdx : ListRowSelect) {
    f_basis[eIdx] = 1;
  }
  MyMatrix<T> VertBasis_T = SelectRow(EXT_T, ListRowSelect);
  MyMatrix<T> VertBasisInv_T = TransposedMat(Inverse(VertBasis_T));
  std::vector<MyVector<T>> VertBasisRed_T;
  for (int i=0; i<=n; i++) {
    MyVector<T> V(n);
    for (int u=0; u<n; u++) {
      V(u) = VertBasis_T(i,u+1);
    }
    VertBasisRed_T.push_back(V);
  }
  return {n, f_basis, std::move(VertBasis_T), std::move(VertBasisInv_T), std::move(VertBasisRed_T)};
}


template<typename T, typename Tvert>
MyVector<T> VoronoiLinearInequality(VoronoiInequalityPreComput<T> const& vipc, MyVector<Tvert> const& TheVert, std::vector<std::vector<T>> const& ListGram) {
  int n = vipc.n;
  int dimSpace = ListGram.size();
  MyVector<T> TheVert_T = UniversalVectorConversion<T,Tvert>(TheVert);
  MyVector<T> B = vipc.VertBasisInv_T * TheVert_T.transpose();
  MyVector<T> Ineq(dimSpace);
  MyVector<T> TheVertRed(n);
  for (int u=0; u<n; u++) {
    TheVertRed(u) = TheVert_T(u+1);
  }
  int iGram = 0;
  for (auto & eLineMat : ListGram) {
    T val = EvaluateLineVector(eLineMat, TheVertRed);
    for (int k=0; k<=n; k++) {
      val += B(k) * EvaluateLineVector(eLineMat, vipc.VertBasisRed_T[k]);
    }
    Ineq(iGram) = val;
    iGram++;
  }
  return Ineq;
}

template<typename T, typename Tvert>
bool IsDelaunayPolytopeInducingEqualities(MyMatrix<Tvert> const& EXT, LinSpaceMatrix<T> const &LinSpa) {
  int n_row = EXT.rows();
  VoronoiInequalityPreComput<T> vipc = BuildVoronoiIneqPreCompute<T,Tvert>(EXT);
  for (int i_row=0; i_row<n_row; i_row++) {
    if (vipc.f_basis[i_row] == 0) {
      MyVector<Tvert> TheVert = GetMatrixRow(EXT, i_row);
      MyVector<T> V = VoronoiLinearInequality(vipc, TheVert, LinSpa.ListLineMat);
      if (!IsZeroVector(V)) {
        return true;
      }
    }
  }
  return false;
}

template<typename T, typename Tvert>
bool IsDelaunayAcceptableForGramMat(MyMatrix<Tvert> const& EXT, LinSpaceMatrix<T> const &LinSpa, MyMatrix<T> const& TestGram) {
  int n_row = EXT.rows();
  MyVector<T> TestV = LINSPA_GetVectorOfMatrixExpression(LinSpa, TestGram);
  VoronoiInequalityPreComput<T> vipc = BuildVoronoiIneqPreCompute<T,Tvert>(EXT);
  for (int i_row=0; i_row<n_row; i_row++) {
    if (vipc.f_basis[i_row] == 0) {
      MyVector<Tvert> TheVert = GetMatrixRow(EXT, i_row);
      MyVector<T> V = VoronoiLinearInequality(vipc, TheVert, LinSpa.ListLineMat);
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

namespace boost::serialization {
  template <class Archive>
  inline void serialize(Archive &ar, AdjInfo &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("iOrb", eRec.iOrb);
    ar &make_nvp("i_adj", eRec.i_adj);
  }
}

template<typename T>
struct FullAdjInfo {
  MyVector<T> eIneq;
  std::vector<AdjInfo> ListAdjInfo;
};

namespace boost::serialization {
  template <class Archive, typename T>
  inline void serialize(Archive &ar, FullAdjInfo<T> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("eIneq", eRec.eIneq);
    ar &make_nvp("ListAdjInfo", eRec.ListAdjInfo);
  }
}

/*
  Compute the defining inequalities of an iso-Delaunay domain
 */
template<typename T, typename Tvert, typename Tgroup>
std::vector<FullAdjInfo<T>> ComputeDefiningIneqIsoDelaunayDomain(DelaunayTesselation<Tvert, Tgroup> const& DT, std::vector<std::vector<T>> const& ListGram) {
  std::unordered_map<MyVector<T>,std::vector<AdjInfo>> map;
  int n_del = DT.l_dels.size();
  for (int i_del=0; i_del<n_del; i_del++) {
    int n_adj = DT.l_dels[i_del].ListAdj.size();
    VoronoiInequalityPreComput<T> vipc = BuildVoronoiIneqPreCompute<T,Tvert>(DT.l_dels[i_del].EXT);
    ContainerMatrix<Tvert> cont(DT.l_dels[i_del].EXT);
    auto get_ineq=[&](int const& i_adj) -> MyVector<T> {
      Delaunay_AdjO<Tvert> adj = DT.l_dels[i_del].ListAdj[i_adj];
      int j_del = adj.iOrb;
      MyMatrix<Tvert> EXTadj = DT.l_dels[j_del].EXT * adj.eBigMat;
      int len = EXTadj.rows();
      for (int u=0; u<len; u++) {
        MyVector<Tvert> TheVert = GetMatrixRow(EXTadj, u);
        std::optional<size_t> opt = cont.GetIdx_v(TheVert);
        if (!opt) {
          return VoronoiLinearInequality(vipc, TheVert, ListGram);
        }
      }
      std::cerr << "Failed to find a matching entry\n";
      throw TerminalException{1};
    };
    for (int i_adj=0; i_adj<n_adj; i_adj++) {
      MyVector<T> V = get_ineq(i_adj);
      // That canonicalization is incorrect because it is not invariant under the group of transformation.
      // The right group transformations are the ones that
      // The list of matrices ListGram must be an integral basis of the T-space. This forces the transformations
      // to be integral in that basis and so the canonicalization by the integer
      MyVector<T> V_red = ScalarCanonicalizationVector(V);
      AdjInfo eAdj{i_del, i_adj};
      map[V_red].push_back(eAdj);
    }
  }
  std::vector<FullAdjInfo<T>> l_ret;
  for (auto & kv : map) {
    FullAdjInfo<T> fai{kv.first, kv.second};
    l_ret.push_back(fai);
  }
  return l_ret;
}

template<typename T>
MyMatrix<T> GetFACineq(std::vector<FullAdjInfo<T>> const& ListIneq) {
  int dim_spa = ListIneq[0].eIneq.size();
  int n_ineq = ListIneq.size();
  MyMatrix<T> FAC(n_ineq, dim_spa);
  for (int i_ineq=0; i_ineq<n_ineq; i_ineq++) {
    for (int u=0; u<dim_spa; u++) {
      FAC(i_ineq, u) = ListIneq[i_ineq].eIneq(u);
    }
  }
  return FAC;
}

template<typename T, typename Tvert, typename Tgroup>
MyMatrix<T> GetInteriorGramMatrix(LinSpaceMatrix<T> const &LinSpa, DelaunayTesselation<Tvert, Tgroup> const& DT, std::ostream& os) {
  int n = LinSpa.n;
  std::vector<FullAdjInfo<T>> ListIneq = ComputeDefiningIneqIsoDelaunayDomain<T,Tvert,Tgroup>(DT, LinSpa.ListLineMat);
  MyMatrix<T> FAC = GetFACineq(ListIneq);
  MyVector<T> ThePt = GetGeometricallyUniqueInteriorPoint(FAC, os);
  MyMatrix<T> RetMat = ZeroMatrix<T>(n, n);
  for (int u=0; u<n; u++) {
    RetMat += ThePt(u) * LinSpa.ListMat[u];
  }
  return RetMat;
}






template<typename T, typename Tint, typename Tgroup>
DelaunayTesselation<Tint, Tgroup> GetInitialGenericDelaunayTesselation(DataIsoDelaunayDomains<T,Tint,Tgroup> const& eData) {
  std::ostream & os = eData.rddo.os;
  auto f_incorrect=[&](Delaunay_Obj<Tint,Tgroup> const& x) -> bool {
    MyMatrix<Tint> const& EXT = x.EXT;
    bool test1 = IsDelaunayPolytopeInducingEqualities(EXT, eData.LinSpa);
    if (test1) {
      return true;
    }
    if (eData.CommonGramMat) {
      MyMatrix<T> const& TestGram = *eData.CommonGramMat;
      bool test2 = IsDelaunayAcceptableForGramMat(EXT, eData.LinSpa, TestGram);
      if (!test2) {
        return true;
      }
    }
    return false;
  };
  auto test_matrix=[&](MyMatrix<T> const& GramMat) -> std::optional<DelaunayTesselation<Tint, Tgroup>> {
    bool test = IsSymmetryGroupCorrect<T,Tint>(GramMat, eData.LinSpa, os);
    if (!test) {
      return {};
    }
    DataLattice<T,Tint,Tgroup> eData = GetDataLattice<T,Tint,Tgroup>(GramMat, os);
    int max_runtime_second = 0;
    return EnumerationDelaunayPolytopes<T, Tint, Tgroup, decltype(f_incorrect)>(eData, f_incorrect, max_runtime_second);
  };
  while(true) {
    MyMatrix<T> GramMat = GetRandomPositiveDefiniteNoNontrialSymm<T,Tint>(eData.LinSpa, os);
    std::optional<DelaunayTesselation<Tint, Tgroup>> opt = test_matrix(GramMat);
    if (opt) {
      return *opt;
    }
  }
  std::cerr << "Failed to find a matching entry\n";
  throw TerminalException{1};
}

template<typename Tvert, typename Tgroup>
struct RepartEntry {
  MyMatrix<Tvert> EXT;
  Tgroup TheStab;
  int8_t Position; // -1: lower, 0: barrel, 1: higher
  int iDelaunayOrigin;
  std::vector<Delaunay_AdjO<Tvert>> ListAdj;
  MyMatrix<Tvert> eBigMat;
};

template<typename Tvert>
std::vector<MyVector<Tvert>> Orbit_MatrixGroup(std::vector<MyMatrix<Tvert>> const& ListGen, MyVector<Tvert> const& eV, std::ostream & os) {
  auto f_prod=[](MyVector<Tvert> const& v, MyMatrix<Tvert> const& M) -> MyVector<Tvert> {
    return M.transpose() * v;
  };
  return OrbitComputation<MyMatrix<Tvert>,MyVector<Tvert>,decltype(f_prod)>(ListGen, eV, f_prod, os);
}


template<typename T, typename Tvert, typename Tgroup>
std::vector<RepartEntry<Tvert, Tgroup>> FindRepartitionningInfoNextGeneration(size_t eIdx, DelaunayTesselation<Tvert, Tgroup> const& ListOrbitDelaunay, std::vector<AdjInfo> const& ListInformationsOneFlipping, MyMatrix<T> const& InteriorElement, RecordDualDescOperation<T, Tgroup> & rddo) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n = InteriorElement.rows();
  std::vector<std::vector<Tidx>> ListPermGenList;
  std::vector<MyMatrix<Tvert>> ListMatGens;
  std::vector<MyVector<Tvert>> ListVertices;
  std::unordered_map<MyVector<Tvert>, size_t> ListVertices_rev;
  size_t n_vertices = 1;
  Tgroup PermGRP;
  auto StandardGroupUpdate=[&]() -> void {
    std::vector<Telt> ListGen;
    Tidx n_act = n_vertices;
    for (auto & eList : ListPermGenList) {
      Telt x(eList);
      ListGen.push_back(x);
    }
    PermGRP = Tgroup(ListGen, n_act);
  };
  StandardGroupUpdate();
  struct TypeOrbitCenter {
    int iDelaunay;
    MyMatrix<Tvert> eBigMat;
    bool status;
    std::vector<Tidx> Linc;
    MyVector<Tvert> EXT;
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
  auto FuncInsertVertex=[&](MyVector<Tvert> const& eVert) -> void {
    if (ListVertices_rev.count(eVert) == 1) {
      return;
    }
    std::vector<MyVector<Tvert>> O = Orbit_MatrixGroup(ListMatGens, eVert, rddo.os);
    for (auto &eV : O) {
      ListVertices.push_back(eV);
      ListVertices_rev[eV] = n_vertices;
      n_vertices++;
    }
    size_t n_gen = ListPermGenList.size();
    for (size_t i_gen=0; i_gen<n_gen; i_gen++) {
      std::vector<Tidx> & ePermGen = ListPermGenList[i_gen];
      MyMatrix<Tvert> eMatrGen = ListMatGens[i_gen];
      for (auto & eVert : O) {
        MyVector<Tvert> eVertImg = eMatrGen.transpose() * eVert;
        size_t pos = ListVertices_rev.at(eVertImg);
        Tidx pos_idx = static_cast<Tidx>(pos - 1);
        ePermGen.push_back(pos_idx);
      }
    }
  };
  auto get_sub_vertices=[&](Face const& f) -> MyMatrix<Tvert> {
    int siz = f.count();
    MyMatrix<Tvert> EXT(siz, n+1);
    int nVert=f.size();
    int pos = 0;
    for (int iVert=0; iVert<nVert; iVert++) {
      if (f[iVert] == 1) {
        for (int i=0; i<=n; i++) {
          EXT(pos,i) = ListVertices[iVert](i);
        }
        pos++;
      }
    }
    return EXT;
  };
  auto get_v_positions=[&](MyMatrix<Tvert> const& eMat) -> std::optional<std::vector<Tidx>> {
    size_t len = n_vertices - 1;
    std::vector<Tidx> v_ret(len);
    for (size_t i=0; i<len; i++) {
      MyVector<Tvert> eVimg = eMat.transpose() * ListVertices[i];
      size_t pos = ListVertices_rev.at(eVimg);
      if (pos == 0) {
        return {};
      }
      Tidx pos_idx = static_cast<Tidx>(pos - 1);
      v_ret[i] = pos_idx;
    }
    return v_ret;
  };
  auto FuncInsertGenerator=[&](MyMatrix<Tvert> const& eMat) -> void {
    std::optional<std::vector<Tidx>> opt = get_v_positions(eMat);
    auto get_test_belong=[&]() -> bool {
      if (opt) {
        Telt elt(*opt);
        return PermGRP.isin(elt);
      } else {
        std::vector<MyMatrix<Tvert>> LGen = ListMatGens;
        LGen.push_back(eMat);
        for (auto & eVert : ListVertices) {
          for (auto & eVertB : Orbit_MatrixGroup(LGen, eVert, rddo.os)) {
            FuncInsertVertex(eVertB);
          }
        }
        return false;
      }
    };
    if (!get_test_belong()) {
      std::optional<std::vector<Tidx>> opt = get_v_positions(eMat);
      ListPermGenList.push_back(*opt);
      ListMatGens.push_back(eMat);
      StandardGroupUpdate();
    }
  };
  auto FuncInsertCenter=[&](TypeOrbitCenterMin const& TheRec) -> void {
    MyMatrix<Tvert> LVert = ListOrbitDelaunay.l_dels[TheRec.iDelaunay].EXT * TheRec.eBigMat;
    for (int u=0; u<LVert.rows(); u++) {
      MyVector<Tvert> eVert = GetMatrixRow(LVert, u);
      FuncInsertVertex(eVert);
    }
    StandardGroupUpdate();
    auto get_iOrbFound=[&]() -> std::optional<size_t> {
      for (size_t iOrb=0; iOrb<ListOrbitCenter.size(); iOrb++) {
        if (ListOrbitCenter[iOrb].iDelaunay == TheRec.iDelaunay) {
          return iOrb;
        }
      }
      return {};
    };
    std::optional<size_t> opt = get_iOrbFound();
    if (opt) {
      MyMatrix<Tvert> eGen = Inverse(TheRec.eBigMat) * ListOrbitCenter[*opt].eBigMat;
      FuncInsertGenerator(eGen);
    } else {
      MyMatrix<Tvert> const& BigMatR = TheRec.eBigMat;
      MyMatrix<Tvert> BigMatI = Inverse(BigMatR);
      MyMatrix<Tvert> const& EXT = ListOrbitDelaunay.l_dels[TheRec.iDelaunay].EXT;
      for (auto & ePermGen : ListOrbitDelaunay.l_dels[TheRec.iDelaunay].GRP.SmallGeneratingSet()) {
        MyMatrix<Tvert> eBigMat = RepresentVertexPermutation(EXT, EXT, ePermGen);
        MyMatrix<Tvert> eBigMat_new = BigMatI * eBigMat * BigMatR;
        FuncInsertGenerator(eBigMat_new);
      }
      std::vector<Tidx> Linc;
      for (int u=0; u<LVert.rows(); u++) {
        MyVector<Tvert> eV = GetMatrixRow(LVert, u);
        size_t pos = ListVertices_rev.at(eV);
        Tidx pos_idx = static_cast<Tidx>(pos - 1);
        Linc.push_back(pos_idx);
      }
      TypeOrbitCenter OrbCent{TheRec.iDelaunay, TheRec.eBigMat, false, Linc, LVert};
      ListOrbitCenter.push_back(OrbCent);
    }
  };
  TypeOrbitCenterMin TheRec{static_cast<int>(eIdx), IdentityMat<Tvert>(n+1)};
  FuncInsertCenter(TheRec);
  while(true) {
    bool IsFinished = true;
    size_t nbCent = ListOrbitCenter.size();
    for (size_t iCent=0; iCent<nbCent; iCent++) {
      TypeOrbitCenter & eEnt = ListOrbitCenter[iCent];
      if (!eEnt.status) {
        IsFinished = false;
        eEnt.status = true;
        for (auto & eCase : ListInformationsOneFlipping) {
          if (eEnt.iDelaunay == eCase.iOrb) {
            MyMatrix<Tvert> const& eBigMat = ListOrbitDelaunay.l_dels[eCase.iOrb].ListAdj[eCase.i_adj].eBigMat;
            MyMatrix<Tvert> eBigMatNew = eBigMat * eEnt.eBigMat;
            TypeOrbitCenterMin TheRec{eCase.iOrb, std::move(eBigMatNew)};
            FuncInsertCenter(TheRec);
          }
        }
      }
    }
    if (IsFinished) {
      break;
    }
  }
  // second part, the convex decomposition
  int nVert = ListVertices.size();
  MyMatrix<T> TotalListVertices(nVert, 2+n);
  std::vector<T> LineInterior = GetLineVector(InteriorElement);
  MyVector<T> eV(n);
  for (int iVert=0; iVert<nVert; iVert++) {
    MyVector<Tvert> const& eVert = ListVertices[iVert];
    TotalListVertices(iVert, 0) = 1;
    for (int iCol=1; iCol<n; iCol++) {
      T val = UniversalScalarConversion<T,Tvert>(eVert(iCol+1));
      TotalListVertices(iVert, iCol+1) = val;
      eV(iCol) = val;
    }
    T Height = EvaluateLineVector(LineInterior, eV);
    TotalListVertices(iVert, n+1) = Height;
  }
  auto get_incd_status=[&](int iVert, MyVector<T> const& eFac) -> bool {
    T eSum = 0;
    for (int u=0; u<=n+1; u++) {
      eSum += TotalListVertices(iVert,u) * eFac(u);
    }
    return eSum == 0;
  };
  struct RepartEntryProv {
    MyVector<T> eFac;
    std::vector<Tidx> Linc;
    Face Linc_face;
    bool Status; // true: YES, false: NO
  };
  std::vector<RepartEntry<Tvert, Tgroup>> ListOrbitFacet;
  std::vector<RepartEntryProv> ListOrbitFacet_prov;
  for (auto & eRec : ListOrbitCenter) {
    Face Linc_face = VectorToFace(eRec.Linc, nVert);
    MyVector<T> eFac = FindFacetInequality(TotalListVertices, Linc_face);
    Tgroup TheStab;
    bool Status = false;
    int8_t Position = -1;
    std::vector<Delaunay_AdjO<Tvert>> ListAdj;
    int iDelaunayOrigin = eRec.iDelaunay;
    MyMatrix<Tvert> const& eBigMat = eRec.eBigMat;
    MyMatrix<Tvert> const& EXT = eRec.EXT;
    RepartEntry<Tvert, Tgroup> re{EXT, TheStab, Position, iDelaunayOrigin, ListAdj, eBigMat};
    RepartEntryProv rep{eFac, eRec.Linc, Linc_face, Status};
    ListOrbitFacet.push_back(re);
    ListOrbitFacet_prov.push_back(rep);
  }
  auto FuncInsertFacet=[&](MyVector<T> const& eFac) -> Delaunay_AdjO<Tvert> {
    std::vector<Tidx> Linc;
    Face Linc_face(nVert);
    std::vector<MyVector<Tvert>> EXT_list;
    for (int iVert=0; iVert<nVert; iVert++) {
      if (get_incd_status(iVert, eFac)) {
        Linc.push_back(iVert);
        Linc_face[iVert] = 1;
        EXT_list.push_back(ListVertices[iVert]);
      }
    }
    MyMatrix<Tvert> EXT = MatrixFromVectorFamily(EXT_list);
    int nOrb = ListOrbitFacet.size();
    for (int iOrb=0; iOrb<nOrb; iOrb++) {
      std::optional<Telt> opt = PermGRP.RepresentativeAction_OnSets(ListOrbitFacet_prov[iOrb].Linc_face, Linc_face);
      if (opt) {
        MyMatrix<Tvert> const& EXT = ListOrbitFacet[iOrb].EXT;
        MyMatrix<Tvert> eBigMat = RepresentVertexPermutation(EXT, EXT, *opt);
        return {{}, std::move(eBigMat), iOrb};
      }
    }
    // A new facet is either barrel or higher because we put all the loiwer ones already
    int8_t Position = 1;
    if (eFac(n+1) == 0) {
      Position = 0;
    }
    bool Status = false;
    Tgroup TheStab;
    int iDelaunayOrigin = -1;
    std::vector<Delaunay_AdjO<Tvert>> ListAdj;
    MyMatrix<Tvert> eMatUnused; // That matrix should never be used
    RepartEntry<Tvert, Tgroup> re{EXT, TheStab, Position, iDelaunayOrigin, ListAdj, eMatUnused};
    RepartEntryProv rep{eFac, Linc, Linc_face, Status};
    ListOrbitFacet.push_back(re);
    ListOrbitFacet_prov.push_back(rep);
    int iOrb = nOrb;
    MyMatrix<Tvert> eBigMat = IdentityMat<Tvert>(n+1);
    return { {}, std::move(eBigMat), iOrb};
  };
  using Text_int = typename SubsetRankOneSolver<T>::Tint;
  MyMatrix<Text_int> TotalListVertices_int = Get_EXT_int(TotalListVertices);
  while(true) {
    bool IsFinished = true;
    int nOrb = ListOrbitFacet.size();
    for (int iOrb=0; iOrb<nOrb; iOrb++) {
      if (!ListOrbitFacet_prov[iOrb].Status) {
        Face const& Linc_face = ListOrbitFacet_prov[iOrb].Linc_face;
        Tgroup Stab = PermGRP.Stabilizer_OnSets(Linc_face);
        Tgroup TheStab = ReducedGroupAction(Stab, Linc_face);
        ListOrbitFacet[iOrb].TheStab = TheStab;
        std::vector<Delaunay_AdjO<Tvert>> ListAdj;
        MyMatrix<Tvert> EXT1 = get_sub_vertices(Linc_face);
        MyMatrix<T> EXT2 = UniversalMatrixConversion<T,Tvert>(EXT1);
        vectface vf = DualDescriptionRecord(EXT2, TheStab, rddo);
        FlippingFramework<T> frame(TotalListVertices, TotalListVertices_int, Linc_face, rddo.os);
        for (auto & eFace : vf) {
          Face eInc = frame.FlipFace(eFace);
          MyVector<T> eFac = FindFacetInequality(TotalListVertices, eInc);
          Delaunay_AdjO<Tvert> eAdj = FuncInsertFacet(eFac);
          size_t nVert = ListOrbitFacet_prov[iOrb].Linc.size();
          Face LEV(nVert);
          for (size_t iInc=0; iInc<nVert; iInc++) {
            Tidx jInc = ListOrbitFacet_prov[iOrb].Linc[iInc];
            if (get_incd_status(jInc, eFac)) {
              LEV[iInc] = 1;
            }
          }
          eAdj.eInc = LEV;
          ListAdj.push_back(eAdj);
        }
        ListOrbitFacet[iOrb].ListAdj = ListAdj;
        ListOrbitFacet_prov[iOrb].Status = true;
        IsFinished = false;
      }
    }
    if (IsFinished) {
      break;
    }
  }
  return ListOrbitFacet;
}


template<typename T, typename Tvert, typename Tgroup>
DelaunayTesselation<Tvert, Tgroup> FlippingLtype(DelaunayTesselation<Tvert, Tgroup> const& ListOrbitDelaunay, MyMatrix<T> const& InteriorElement, std::vector<AdjInfo> const& ListInformationsOneFlipping, RecordDualDescOperation<T, Tgroup> & rddo) {
  using Tgr = GraphListAdj;
  using Telt = typename Tgroup::Telt;
  int n_dels = ListOrbitDelaunay.l_dels.size();
  Tgr Gra(n_dels);
  Face ListMatched(n_dels);
  for (auto & eAI : ListInformationsOneFlipping) {
    ListMatched[eAI.iOrb] = 1;
    int iOrbAdj = ListOrbitDelaunay.l_dels[eAI.iOrb].ListAdj[eAI.i_adj].iOrb;
    if (eAI.iOrb != iOrbAdj) {
      Gra.AddAdjacent(eAI.iOrb, iOrbAdj);
    }
  }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  if (!IsSymmetricGraph(Gra)) {
    std::cerr << "ISO_DEL: The graph is not symmetric\n";
    throw TerminalException{1};
  }
#endif
  std::vector<std::vector<size_t>> ListConn = ConnectedComponents_set(Gra);
  std::vector<std::vector<size_t>> ListGroupMelt, ListGroupUnMelt;
  auto is_melt=[&](std::vector<size_t> const& eConn) -> bool {
    size_t n_matched = 0;
    for (auto & eVal : eConn) {
      n_matched += ListMatched[eVal];
    }
    if (n_matched == eConn.size()) {
      return true;
    }
    if (n_matched == 0) {
      return false;
    }
    std::cerr << "ISO_DEL: The melt should be total or none at all\n";
    throw TerminalException{1};
  };
  for (auto & eConn : ListConn) {
    if (is_melt(eConn)) {
      ListGroupMelt.push_back(eConn);
    } else {
      ListGroupUnMelt.push_back(eConn);
    }
  }
  std::vector<std::vector<RepartEntry<Tvert, Tgroup>>> ListInfo;
  std::vector<int> vect_iInfo(n_dels, -1);
  std::vector<int> vect_lower_iFacet(n_dels, -1);
  int iInfo = 0;
  for (auto &eConn : ListGroupMelt) {
    size_t eIdx = eConn[0];
    std::vector<RepartEntry<Tvert, Tgroup>> LORB2 = FindRepartitionningInfoNextGeneration(eIdx, ListOrbitDelaunay, ListInformationsOneFlipping, InteriorElement, rddo);
    int n_facet = LORB2.size();
    for (int iFacet=0; iFacet<n_facet; iFacet++) {
      RepartEntry<Tvert, Tgroup> const& eFacet = LORB2[iFacet];
      if (eFacet.Position == -1) {
        vect_lower_iFacet[eFacet.iDelaunayOrigin] = iFacet;
      }
    }
    ListInfo.push_back(LORB2);
    for (auto &pos : eConn) {
      vect_iInfo[pos] = iInfo;
    }
    iInfo++;
  }
  struct MatchedFacet {
    Delaunay_AdjO<Tvert> adj;
    MyMatrix<Tvert> eBigMat;
  };
  auto get_matching_listinfo=[&](int iInfo, int iFacet, Face eInc) -> MatchedFacet {
    MyMatrix<Tvert> const& EXT = ListInfo[iInfo][iFacet].EXT;
    for (auto &eAdj : ListInfo[iInfo][iFacet].ListAdj) {
      std::optional<Telt> opt = ListInfo[iInfo][iFacet].TheStab.RepresentativeAction_OnSets(eAdj.eInc, eInc);
      if (opt) {
        MyMatrix<Tvert> eBigMat = RepresentVertexPermutation(EXT, EXT, *opt);
        return {eAdj, eBigMat};
      }
    }
    std::cerr << "Failed to find a matching entry in get_matching_listinfo\n";
    throw TerminalException{1};
  };
  auto get_matching_old_tessel=[&](int iDelaunayOld, Face eInc) -> MatchedFacet {
    MyMatrix<Tvert> const& EXT = ListOrbitDelaunay.l_dels[iDelaunayOld].EXT;
    for (auto &eAdj : ListOrbitDelaunay.l_dels[iDelaunayOld].ListAdj) {
      std::optional<Telt> opt = ListOrbitDelaunay.l_dels[iDelaunayOld].GRP.RepresentativeAction_OnSets(eAdj.eInc, eInc);
      if (opt) {
        MyMatrix<Tvert> eBigMat = RepresentVertexPermutation(EXT, EXT, *opt);
        return {eAdj, eBigMat};
      }
    }
    std::cerr << "Failed to find a matching entry in get_matching_listinfo\n";
    throw TerminalException{1};
  };
  auto get_lower_adjacency=[&](int iInfo, int iFacet) -> Delaunay_AdjO<Tvert> {
    for (auto & fAdj : ListInfo[iInfo][iFacet].ListAdj) {
      if (ListInfo[iInfo][fAdj.iOrb].Position == -1) {
        return fAdj;
      }
    }
    std::cerr << "Failed to find a lower facet\n";
    throw TerminalException{1};
  };
  auto get_face_m_m=[&](MyMatrix<Tvert> const M1, MyMatrix<Tvert> const& M2) -> Face {
    int nVert1=M1.rows();
    int nVert2=M2.rows();
    Face f2(nVert2);
    ContainerMatrix<Tvert> cont(M2);
    for (int i1=0; i1<nVert1; i1++) {
      MyVector<Tvert> V = GetMatrixRow(M1, i1);
      std::optional<size_t> opt = cont.GetIdx_v(V);
      size_t pos2 = unfold_opt(opt, "Error in get_face_m_m");
      f2[pos2] = 1;
    }
    return f2;
  };
  auto get_face_msub_m=[&](MyMatrix<Tvert> const M1, Face const& f1, MyMatrix<Tvert> const& M2) -> Face {
    int nVert1=M1.rows();
    int nVert2=M2.rows();
    Face f2(nVert2);
    ContainerMatrix<Tvert> cont(M2);
    for (int i1=0; i1<nVert1; i1++) {
      if (f1[i1] == 1) {
        MyVector<Tvert> V = GetMatrixRow(M1, i1);
        std::optional<size_t> opt = cont.GetIdx_v(V);
        int pos2 = unfold_opt(opt, "Error in get_face_m_m");
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
  auto get_symbol_position=[&](DelaunaySymb const& ds) -> std::optional<size_t> {
    if (ds.Position == Position_old) {
      for (size_t i=0; i<NewListOrbitDelaunay.size(); i++) {
        if (NewListOrbitDelaunay[i].Position == Position_old) {
          if (NewListOrbitDelaunay[i].iDelaunay == ds.iDelaunay) {
            return i;
          }
        }
      }
    }
    if (ds.Position == Position_new) {
      for (size_t i=0; i<NewListOrbitDelaunay.size(); i++) {
        if (NewListOrbitDelaunay[i].Position == Position_new) {
          if (NewListOrbitDelaunay[i].iInfo == ds.iInfo && NewListOrbitDelaunay[i].iFacet == ds.iFacet) {
            return i;
          }
        }
      }
    }
    return {};
  };
  std::vector<Delaunay_Entry<Tvert,Tgroup>> l_dels;
  for (auto & eConn : ListGroupUnMelt) {
    if (eConn.size() > 1) {
      std::cerr << "Error of connected component computation\n";
      throw TerminalException{1};
    }
    int iDelaunay = eConn[0];
    DelaunaySymb ds{Position_old, iDelaunay, -1, -1};
    NewListOrbitDelaunay.push_back(ds);
    Delaunay_Entry<Tvert, Tgroup> del{ListOrbitDelaunay.l_dels[iDelaunay].EXT, ListOrbitDelaunay.l_dels[iDelaunay].GRP, {}};
    l_dels.push_back(del);
  }
  for (int iInfo=0; iInfo<n_info; iInfo++) {
    int n_facet = ListInfo[iInfo].size();
    for (int iFacet=0; iFacet<n_facet; iFacet++) {
      int8_t Position = ListInfo[iInfo][iFacet].Position;
      if (Position == 1) {
        DelaunaySymb ds{Position_old, -1, iInfo, iFacet};
        NewListOrbitDelaunay.push_back(ds);
        MyMatrix<Tvert> EXT = ListInfo[iInfo][iFacet].EXT;
        Tgroup GRP = ListInfo[iInfo][iFacet].TheStab;
        Delaunay_Entry<Tvert, Tgroup> del{EXT, GRP, {} };
        l_dels.push_back(del);
      }
    }
  }
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
  auto check_adj=[&](int iOrb, Delaunay_AdjO<Tvert> const& NewAdj, std::string const& context) -> void {
    MyMatrix<Tvert> const& EXT = l_dels[iOrb].EXT;
    ContainerMatrix<Tint> cont(EXT);
    Face f_att(EXT.rows());
    MyMatrix<Tvert> EXTadj = l_dels[NewAdj.iOrb] * NewAdj.eBigMat;
    int len = EXTadj.rows();
    for (int iVert=0; iVert<len; iVert++) {
      MyVector<Tvert> V = GetMatrixRow(EXTadj, u);
      std::optional<size_t> opt = cont.GetIdx_v(V);
      if (opt) {
        f_att[*opt] = 1;
      }
    }
    if (f_att != NewAdj.f) {
      std::cerr << "Consistency error in context=" << context << "\n";
      throw TerminalException{1};
    }
  };
#endif
  int n_del_ret = l_dels.size();
  for (int iOrb=0; iOrb<n_del_ret; iOrb++) {
    DelaunaySymb ds = NewListOrbitDelaunay[iOrb];
    std::vector<Delaunay_AdjO<Tvert>> ListAdj;
    if (ds.Position == Position_old) {
      int iDelaunay = ds.iDelaunay;
      for (auto & eAdj : ListOrbitDelaunay.l_dels[iDelaunay].ListAdj) {
        int iDelaunayOld = eAdj.iOrb;
        DelaunaySymb dss{Position_old, iDelaunayOld, -1, -1};
        std::optional<size_t> opt = get_symbol_position(dss);
        if (opt) {
          int iOrb = *opt;
          Delaunay_AdjO<Tvert> NAdj{eAdj.eInc, eAdj.eBigMat, iOrb};
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
          check_adj(iOrb, NAdj, "Case 1");
#endif
          ListAdj.push_back(NAdj);
        } else {
          int iInfo = vect_iInfo[iDelaunayOld];
          int iFacet = vect_lower_iFacet[iDelaunayOld];
          RepartEntry<Tvert, Tgroup> const& eFacet = ListInfo[iInfo][iFacet];
          MyMatrix<Tvert> const& BigMat2 = eFacet.eBigMat;
          MyMatrix<Tvert> ImageEXT = ListOrbitDelaunay.l_dels[iDelaunayOld].EXT * eAdj.eBigMat;
          Face Linc = get_face_msub_m(ListOrbitDelaunay.l_dels[iDelaunay].EXT, eAdj.eInc, ImageEXT);
          MatchedFacet RecMatch = get_matching_listinfo(iInfo, iFacet, Linc);
          int iOrbFound = RecMatch.adj.iOrb;
          if (ListInfo[iInfo][iFacet].Position == 0) {
            std::cerr << "Illogic error concerning the structure of repartitionning polytope\n";
            throw TerminalException{1};
          }
          DelaunaySymb dss{Position_new, -1, iInfo, iOrbFound};
          std::optional<int> optN = get_symbol_position(dss);
          int Pos = unfold_opt(optN, "Failed to find entry for Case2");
          MyMatrix<Tvert> BigMat1 = RecMatch.adj.eBigMat * RecMatch.eBigMat * Inverse(BigMat2) * eAdj.eBigMat;
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
      for (auto & eAdj : ListInfo[iInfo][iFacet].ListAdj) {
        int jFacet = eAdj.iOrb;
        if (ListInfo[iInfo][jFacet].Position == 0) {
          MyMatrix<Tvert> const& eMat1 = eAdj.eBigMat;
          MyMatrix<Tvert> LincEXT = SelectRow(ListInfo[iInfo][iFacet].EXT, eAdj.eInc);
          MyMatrix<Tvert> ImageEXTbarrel = ListInfo[iInfo][jFacet].EXT * eMat1;
          Face LLinc = get_face_m_m(LincEXT, ImageEXTbarrel);
          Delaunay_AdjO<Tvert> TheFoundAdj = get_lower_adjacency(iInfo, jFacet);
          int kFacet = TheFoundAdj.iOrb;
          MyMatrix<Tvert> const& eMat2 = TheFoundAdj.eBigMat;
          int iDelaunayOrigin = ListInfo[iInfo][kFacet].iDelaunayOrigin;
          MyMatrix<Tvert> ImageEXT = ListInfo[iInfo][kFacet].EXT * TheFoundAdj.eBigMat;
          Face LLinc2 = get_face_msub_m(ListInfo[iInfo][jFacet].EXT, TheFoundAdj.eInc, ImageEXT);
          MatchedFacet match2 = get_matching_old_tessel(iDelaunayOrigin, LLinc2);
          Delaunay_AdjO<Tvert> const& TheFoundAdj2 = match2.adj;
          MyMatrix<Tvert> const& TheMat2 = match2.eBigMat;
          int jDelaunayOrigin = TheFoundAdj2.iOrb;
          MyMatrix<Tvert> const& eMat3 = TheFoundAdj2.eBigMat;
          ImageEXT = ListOrbitDelaunay.l_dels[jDelaunayOrigin].EXT * eMat3;
          Face LLinc3 = get_face_msub_m(ListOrbitDelaunay.l_dels[iDelaunayOrigin].EXT, TheFoundAdj2.eInc, ImageEXT);
          int jInfo = vect_iInfo[jDelaunayOrigin];
          int iFacet2 = vect_lower_iFacet[jDelaunayOrigin];
          MatchedFacet match3 = get_matching_listinfo(jInfo, iFacet2, LLinc3);
          Delaunay_AdjO<Tvert> const& TheFoundAdj3 = match3.adj;
          MyMatrix<Tvert> const& TheMat3 = match3.eBigMat;
          MyMatrix<Tvert> BigMat1 = TheFoundAdj3.eBigMat*TheMat3*Inverse(ListInfo[jInfo][iFacet2].eBigMat)*eMat3*TheMat2*ListInfo[iInfo][kFacet].eBigMat*eMat2*eMat1;
          MyMatrix<Tvert> EXT7 = ListInfo[jInfo][TheFoundAdj3.iOrb].EXT * BigMat1;
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
          if (SortMatrix(EXT7) != SortMatrix(ImageEXTbarrel)) {
            std::cerr << "We fail an important test with the barrel images\n";
            throw TerminalException{1};
          }
#endif
          Face LLinc4 = get_face_m_m(LincEXT, EXT7);
          MatchedFacet match4 = get_matching_listinfo(jInfo, TheFoundAdj3.iOrb, LLinc4);
          Delaunay_AdjO<Tvert> const& TheFoundAdj4 = match4.adj;
          MyMatrix<Tvert> const& TheMat4 = match4.eBigMat;
          DelaunaySymb dss{Position_new, -1, jInfo, TheFoundAdj4.iOrb};
          std::optional<size_t> optN = get_symbol_position(dss);
          int Pos = unfold_opt(optN, "Failed to find entry for Case3");
          MyMatrix<Tvert> BigMat2 = TheFoundAdj4.eBigMat * TheMat4 * BigMat1;
          Delaunay_AdjO<Tvert> NAdj{eAdj.eInc, BigMat2, Pos};
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
          check_adj(iOrb, NAdj, "Case 3");
#endif
          ListAdj.push_back(NAdj);
        }
        if (ListInfo[iInfo][jFacet].Position == -1) {
          int iDelaunayOrigin = ListInfo[iInfo][jFacet].iDelaunayOrigin;
          MyMatrix<Tvert> const& eMat1 = eAdj.eBigMat;
          MyMatrix<Tvert> ImageEXT = ListInfo[iInfo][jFacet].EXT * eMat1;
          Face LLinc = get_face_msub_m(ListInfo[iInfo][iFacet].EXT, eAdj.eInc, ImageEXT);
          MatchedFacet match1 = get_matching_old_tessel(iDelaunayOrigin, LLinc);
          Delaunay_AdjO<Tvert> const& TheFoundAdj1 = match1.adj;
          MyMatrix<Tvert> const& TheMat1 = match1.eBigMat;
          int jDelaunayOld = TheFoundAdj1.iOrb;
          if (vect_iInfo[jDelaunayOld] == -1) {
            DelaunaySymb dss{Position_old, jDelaunayOld, -1, -1};
            std::optional<size_t> opt = get_symbol_position(dss);
            int Pos2 = unfold_opt(opt, "Case 4");
            MyMatrix<Tvert> BigMat1 = TheFoundAdj1.eBigMat*TheMat1*ListInfo[iInfo][jFacet].eBigMat*eAdj.eBigMat;
            Delaunay_AdjO<Tvert> NAdj{eAdj.eInc, BigMat1, Pos2};
#ifdef DEBUG_ISO_DELAUNAY_DOMAIN
            check_adj(iOrb, NAdj, "Case 4");
#endif
            ListAdj.push_back(NAdj);
          } else {
            int jInfo = vect_iInfo[jDelaunayOld];
            int iFacet2 = vect_lower_iFacet[jDelaunayOld];
            MyMatrix<Tvert> ImageEXT = ListOrbitDelaunay.l_dels[jDelaunayOld].EXT * TheFoundAdj1.eBigMat;
            Face LLinc2 = get_face_msub_m(ListOrbitDelaunay.l_dels[iDelaunayOrigin].EXT, TheFoundAdj1.eInc, ImageEXT);
            MatchedFacet match2 = get_matching_listinfo(jInfo, iFacet2, LLinc2);
            Delaunay_AdjO<Tvert> const& TheFoundAdj2 = match2.adj;
            MyMatrix<Tvert> const& TheMat2 = match2.eBigMat;
            DelaunaySymb dss{Position_new, -1, jInfo, TheFoundAdj2.iOrb};
            std::optional<size_t> opt = get_symbol_position(dss);
            int Pos = unfold_opt(opt, "Case 5");
            MyMatrix<Tvert> BigMat1 = TheFoundAdj2.eBigMat*TheMat2*Inverse(ListInfo[jInfo][iFacet2].eBigMat)*TheFoundAdj1.eBigMat*TheMat1*ListInfo[iInfo][jFacet].eBigMat*eAdj.eBigMat;
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

template<typename T, typename Tint>
struct IsoDelaunayDomain_MPI_AdjO {
  int iProc;
  int iOrb;
  MyVector<T> V;
  MyMatrix<Tint> eBigMat;
};

namespace boost::serialization {
  template <class Archive, typename T, typename Tint>
  inline void serialize(Archive &ar, IsoDelaunayDomain_MPI_AdjO<T, Tint> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("iProc", eRec.iProc);
    ar &make_nvp("iOrb", eRec.iOrb);
    ar &make_nvp("V", eRec.V);
    ar &make_nvp("eBigMat", eRec.eBigMat);
  }
}

template<typename T, typename Tint, typename Tgroup>
struct IsoDelaunayDomain {
  DelaunayTesselation<Tint, Tgroup> DT;
  MyMatrix<T> GramMat;
};

namespace boost::serialization {
  template <class Archive, typename T, typename Tint, typename Tgroup>
  inline void serialize(Archive &ar, IsoDelaunayDomain<T, Tint, Tgroup> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("DT", eRec.DT);
    ar &make_nvp("GramMat", eRec.GramMat);
  }
}

template<typename T, typename Tint, typename Tgroup>
struct IsoDelaunayDomain_AdjI {
  MyVector<T> V;
  IsoDelaunayDomain<T, Tint, Tgroup> DT_gram;
};

namespace boost::serialization {
  template <class Archive, typename T, typename Tint, typename Tgroup>
  inline void serialize(Archive &ar, IsoDelaunayDomain_AdjI<T, Tint, Tgroup> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("V", eRec.V);
    ar &make_nvp("DT_gram", eRec.DT_gram);
  }
}

template<typename T, typename Tint>
struct IsoDelaunayDomain_AdjO {
  MyVector<T> V;
  MyMatrix<Tint> eBigMat;
};

namespace boost::serialization {
  template <class Archive, typename T, typename Tint>
  inline void serialize(Archive &ar, IsoDelaunayDomain_AdjO<T, Tint> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("V", eRec.V);
    ar &make_nvp("eBigMat", eRec.eBigMat);
  }
}

template<typename T, typename Tint, typename Tgroup>
struct IsoDelaunayDomain_MPI_Entry {
  IsoDelaunayDomain<T, Tint, Tgroup> DT_gram;
  std::vector<FullAdjInfo<T>> ListIneq;
  std::vector<IsoDelaunayDomain_MPI_AdjO<T, Tint>> ListAdj;
};

namespace boost::serialization {
  template <class Archive, typename T, typename Tint, typename Tgroup>
  inline void serialize(Archive &ar, IsoDelaunayDomain_MPI_Entry<T, Tint, Tgroup> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("DT_gram", eRec.DT_gram);
    ar &make_nvp("ListIneq", eRec.ListIneq);
    ar &make_nvp("ListAdj", eRec.ListAdj);
  }
}

/*
  We do not use the ComputeInvariantDelaunay function since it requires having the
  GramMatrix, which would be an additional computation.
  ---
  BUT we should have some invariants coming from the Tspace and right now we have none.
 */
template<typename T, typename Tint, typename Tgroup>
size_t ComputeInvariantIsoDelaunayDomain(DataIsoDelaunayDomains<T,Tint,Tgroup> const& eData, size_t const& seed, DelaunayTesselation<Tint, Tgroup> const& DT) {
  using TintGroup = typename Tgroup::Tint;
  std::map<size_t, size_t> map_delaunays;
  auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
    seed ^= new_hash + 0x9e3779b8 + (seed << 6) + (seed >> 2);
  };
  size_t seed_delaunay = 10;
  for (auto & eDel : DT.l_dels) {
    size_t hash = seed_delaunay;
    size_t hash_nVert = eDel.EXT.rows();
    TintGroup order = eDel.GRP.size();
    size_t hash_grp = std::hash<TintGroup>()(order);
    combine_hash(hash, hash_nVert);
    combine_hash(hash, hash_grp);
    std::map<size_t, size_t> map_facesiz;
    for (auto& eAdj : eDel.ListAdj) {
      size_t cnt = eAdj.eInc.count();
      map_facesiz[cnt] += 1;
    }
    for (auto & kv : map_facesiz) {
      combine_hash(hash, kv.first);
      combine_hash(hash, kv.second);
    }
    map_delaunays[hash] += 1;
  }
  size_t hash_ret = seed;
  for (auto & kv : map_delaunays) {
    combine_hash(hash_ret, kv.first);
    combine_hash(hash_ret, kv.second);
  }
  return hash_ret;
}

template<typename T, typename Tint, typename Tgroup>
IsoDelaunayDomain<T, Tint, Tgroup> GetInitialIsoDelaunayDomain(DataIsoDelaunayDomains<T,Tint,Tgroup> const& data) {
  DelaunayTesselation<Tint, Tgroup> DT = GetInitialGenericDelaunayTesselation(data);
  MyMatrix<T> GramMat = GetInteriorGramMatrix(data.LinSpa, DT, data.rddo.os);
  return {DT, GramMat};
}

template<typename T, typename Tint, typename Tgroup>
struct IsoDelaunayDomain_Obj {
  IsoDelaunayDomain<T, Tint, Tgroup> DT_gram;
  std::vector<FullAdjInfo<T>> ListIneq;
};

namespace boost::serialization {
  template <class Archive, typename T, typename Tint, typename Tgroup>
  inline void serialize(Archive &ar, IsoDelaunayDomain_Obj<T, Tint, Tgroup> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("DT_gram", eRec.DT_gram);
    ar &make_nvp("ListIneq", eRec.ListIneq);
  }
}






template<typename T, typename Tint, typename Tgroup>
struct DataIsoDelaunayDomainsFunc {
  DataIsoDelaunayDomains<T,Tint,Tgroup> data;
  using Tobj = IsoDelaunayDomain_Obj<T, Tint, Tgroup>;
  using TadjI = IsoDelaunayDomain_AdjI<T, Tint, Tgroup>;
  using TadjO = IsoDelaunayDomain_AdjO<T, Tint>;
  std::ostream& get_os() {
    return data.rddo.os;
  }
  Tobj f_init() {
    IsoDelaunayDomain<T, Tint, Tgroup> IsoDel = GetInitialIsoDelaunayDomain(data);
    Tobj x{IsoDel, {} };
    return x;
  }
  size_t f_hash(size_t const& seed, Tobj const& x) {
    return ComputeInvariantIsoDelaunayDomain<T,Tint,Tgroup>(data, seed, x.DT_gram.DT);
  }
  std::optional<TadjO> f_repr(Tobj const& x, TadjI const& y) {
    std::optional<MyMatrix<Tint>> opt = LINSPA_TestEquivalenceGramMatrix<T,Tint,Tgroup>(data.LinSpa, x.DT_gram.GramMat, y.DT_gram.GramMat, data.rddo.os);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const& eBigMat = *opt;
    TadjO ret{y.V, eBigMat};
    return ret;
  }
  std::pair<Tobj, TadjO> f_spann(TadjI const& x) {
    IsoDelaunayDomain<T, Tint, Tgroup> IsoDel = x.DT_gram;
    Tobj x_ret{IsoDel, {} };
    MyMatrix<Tint> eBigMat = IdentityMat<Tint>(data.LinSpa.n);
    TadjO ret{x.V, eBigMat};
    return {std::move(x_ret), std::move(ret)};
  }
  std::vector<TadjI> f_adj(Tobj & x_in) {
    IsoDelaunayDomain<T, Tint, Tgroup> & x = x_in.DT_gram;
    std::vector<FullAdjInfo<T>> ListIneq = ComputeDefiningIneqIsoDelaunayDomain<T,Tint,Tgroup>(x.DT, data.LinSpa.ListLineMat);
    x_in.ListIneq = ListIneq;
    MyMatrix<T> FAC = GetFACineq(ListIneq);
    std::vector<int> ListIrred = cdd::RedundancyReductionClarkson(FAC);
    std::vector<TadjI> l_adj;
    for (auto & idxIrred : ListIrred) {
      FullAdjInfo<T> eRecIneq = ListIneq[idxIrred];
      DelaunayTesselation<Tint, Tgroup> DTadj = FlippingLtype<T,Tint,Tgroup>(x.DT, x.GramMat, eRecIneq.ListAdjInfo, data.rddo);
      MyMatrix<T> GramMatAdj = GetInteriorGramMatrix(data.LinSpa, DTadj, data.rddo.os);
      IsoDelaunayDomain<T, Tint, Tgroup> IsoDelAdj{DTadj, GramMatAdj};
      TadjI eAdj{eRecIneq.eIneq, IsoDelAdj};
      l_adj.push_back(eAdj);
    }
    return l_adj;
  }
  Tobj f_adji_obj(TadjI const& x) {
    return {x.DT_gram, {} };
  }
};

template<typename T, typename Tint, typename Tgroup>
void WriteFamilyIsoDelaunayDomain(boost::mpi::communicator &comm, std::string const& OutFormat, std::string const& OutFile, std::vector<DatabaseEntry_MPI<typename DataIsoDelaunayDomainsFunc<T,Tint,Tgroup>::Tobj, typename DataIsoDelaunayDomainsFunc<T,Tint,Tgroup>::TadjO>> const& ListIDD, [[maybe_unused]] std::ostream & os) {
  //  int i_rank = comm.rank();
  if (OutFormat == "nothing") {
    std::cerr << "No output\n";
    return;
  }
  if (OutFormat == "GAP") {
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
  throw TerminalException{1};
}


template<typename T, typename Tint, typename Tgroup>
void ComputeLatticeIsoDelaunayDomains(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  int i_rank = comm.rank();
  int n_proc = comm.size();
  std::string FileLog = "log_" + std::to_string(n_proc) + "_" + std::to_string(i_rank);
  std::ofstream os(FileLog);
  if (ApplyStdUnitbuf(eFull)) {
    os << std::unitbuf;
    os << "Apply UnitBuf\n";
  } else {
    os << "Do not apply UnitBuf\n";
  }
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  SingleBlock BlockTSPACE = eFull.ListBlock.at("TSPACE");
  //
  bool STORAGE_Saving = BlockDATA.ListBoolValues.at("Saving");
  std::string STORAGE_Prefix = BlockDATA.ListStringValues.at("Prefix");
  CreateDirectory(STORAGE_Prefix);
  //
  int max_runtime_second = BlockDATA.ListIntValues.at("max_runtime_second");
  std::cerr << "max_runtime_second=" <<	max_runtime_second << "\n";
  std::string OutFormat = BlockDATA.ListStringValues.at("OutFormat");
  std::string OutFile = BlockDATA.ListStringValues.at("OutFile");
  std::cerr << "OutFormat=" << OutFormat << " OutFile=" << OutFile << "\n";
  auto get_common=[&]() -> std::optional<MyMatrix<T>> {
    std::string CommonGramMat = BlockDATA.ListStringValues.at("CommonGramMat");
    if (CommonGramMat == "unset") {
      return {};
    }
    MyMatrix<T> eMat = ReadMatrixFile<T>(CommonGramMat);
    return eMat;
  };
  std::optional<MyMatrix<T>> CommonGramMat = get_common();
  //
  using TintGroup = typename Tgroup::Tint;
  PolyHeuristicSerial<TintGroup> AllArr = AllStandardHeuristicSerial<TintGroup>(os);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);

  LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint>(BlockTSPACE, os);

  DataIsoDelaunayDomains<T,Tint,Tgroup> data{LinSpa,
    std::move(rddo),
    CommonGramMat};
  using Tdata = DataIsoDelaunayDomainsFunc<T,Tint,Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_MPI<Tobj, TadjO>;
  std::pair<bool, std::vector<Tout>> pair = EnumerateAndStore_MPI<Tdata>(comm, data_func, STORAGE_Prefix, STORAGE_Saving, max_runtime_second);
  if (pair.first) {
    WriteFamilyIsoDelaunayDomain<T, Tint, Tgroup>(comm, OutFormat, OutFile, pair.second, os);
  }
}


// clang-format off
#endif  // SRC_LATT_ISODELAUNAYDOMAINS_H_
// clang-format on

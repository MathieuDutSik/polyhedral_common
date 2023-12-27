// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ISODELAUNAYDOMAINS_H_
#define SRC_LATT_ISODELAUNAYDOMAINS_H_

// clang-format off
#include "POLY_LinearProgramming.h"
#include "ShortestUniversal.h"
#include "Temp_Positivity.h"
#include "LatticeDelaunay.h"
#include <string>
#include <vector>
// clang-format on

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
  Face f_basis;
  MyMatrix<T> VertBasis_T;
  MyMatrix<T> VertBasisInv_T;
  std::vector<MyVector<T>> VertBasisRed_T;
};

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
  return {f_basis, std::move(VertBasis_T), std::move(VertBasisInv_T), std::move(VertBasisRed_T)};
}


template<typename T, typename Tvert>
MyVector<T> VoronoiLinearInequality(VoronoiInequalityPreComput<T> const& vipc, MyVector<Tvert> const& TheVert, std::vector<std::vector<T>> const& ListGram) {
  int n = ListGram[0].rows();
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
bool IsDelaunayPolytopeInducingEqualities(MyMatrix<Tvert> const& EXT, std::vector<std::vector<T>> const& ListGram) {
  int n_row = EXT.rows();
  VoronoiInequalityPreComput<T> vipc = BuildVoronoiIneqPreCompute<T,Tvert>(EXT);
  for (int i_row=0; i_row<n_row; i_row++) {
    if (vipc.f_basis[i_row] == 0) {
      MyVector<Tvert> TheVert = GetMatrixRow(EXT, i_row);
      MyVector<T> V = VoronoiLinearInequality(vipc, TheVert, ListGram);
      if (!IsZeroVector(V)) {
        return true;
      }
    }
  }
  return false;
}



struct AdjInfo {
  int iOrb;
  int i_adj;
};


/*
  Compute the defining inequalities of an iso-Delaunay domain
 */
template<typename T, typename Tvert, typename Tgroup>
std::unordered_map<MyVector<T>,std::vector<AdjInfo>> ComputeDefiningIneqIsoDelaunayDomain(DelaunayTesselation<Tvert, Tgroup> const& DT, std::vector<std::vector<T>> const& ListGram) {
  std::unordered_map<MyVector<T>,std::vector<AdjInfo>> map;
  int n_del = DT.l_dels.size();
  for (int i_del=0; i_del<n_del; i_del++) {
    int n_adj = DT.l_dels.l_adj.size();
    VoronoiInequalityPreComput<T> vipc = BuildVoronoiIneqPreCompute(DT.l_dels[i_del].obj);
    ContainerMatrix<Tvert> cont(DT.l_dels[i_del].obj);
    auto get_ineq=[&](int const& i_adj) -> MyVector<T> {
      Delaunay_AdjO<Tvert> adj = DT.l_dels[i_del].l_adj[i_adj];
      int j_del = adj.iOrb;
      MyMatrix<Tvert> EXTadj = DT.l_dels[j_del].obj * adj.P;
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
  return map;
}


template<typename T, typename Tint>
bool IsSymmetryGroupCorrect(MyMatrix<T> const& GramMat, LinSpaceMatrix<T> const& LinSpa, std::ostream & os) {
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(GramMat);
  int n_row = SHV.rows();
  std::vector<T> Vdiag(n_row,0);
  std::vector<MyMatrix<T>> ListMat = {GramMat};
  const bool use_scheme = true;
  std::vector<std::vector<Tidx>> ListGen =
    GetListGenAutomorphism_ListMat_Vdiag<T, Tfield, Tidx, use_scheme>(SHV_T, ListMat, Vdiag, os);
  for (auto &eList : ListGen) {
    std::optional<MyMatrix<T>> opt =
      FindMatrixTransformationTest(SHV_T, SHV_T, eList);
    if (!opt) {
      std::cerr << "Failed to find the matrix\n";
      throw TerminalException{1};
    }
    MyMatrix<T> const &M_T = *opt;
    if (!IsIntegralMatrix(M_T)) {
      std::cerr << "Bug: The matrix should be integral\n";
      throw TerminalException{1};
    }
    for (auto & eMat : LinSpa.ListMat) {
      MyMatrix<T> eMatImg = M_T * eMat * M_T.transpose();
      if (eMatImg != eMat) {
        return false;
      }
    }
  }
  return true;
}



template<typename T, typename Tint, typename Tgroup>
DelaunayTesselation<Tint, Tgroup> GetInitialGenericDelaunayTesselation(LinSpaceMatrix<T> const& LinSpa, std::ostream & os) {
  auto f_incorrect=[&](MyMatrix<Tint> const& EXT) -> bool {
    return IsDelaunayPolytopeInducingEqualities(EXT, LinSpa.ListLineMat);
  };
  auto test_matrix=[&](MyMatrix<T> const& GramMat) -> std::optional<DelaunayTesselation<Tint, Tgroup>> {
    bool test = IsSymmetryGroupCorrect(GramMat, LinSpa, os);
    if (!test) {
      return {};
    }
    DataLattice<T,Tint,Tgroup> eData = GetDataLattice<T,Tint,Tgroup>(GramMat);
    return EnumerationDelaunayPolytopes<T, Tint, Tgroup, decltype(f_incorrect)>(eData, f_incorrect, os);
  };
  while(true) {
    MyMatrix<T> GramMat = GetRandomPositiveDefinite(LinSpa);
    std::optional<DelaunayTesselation<Tint, Tgroup>> opt = test_matrix();
    if (opt) {
      return *opt;
    }
  }
  std::cerr << "Failed to find a matching entry\n";
  throw TerminalException{1};
}

template<typename Tvert, typename Tgroup>
struct AdjRepart {
  int iOrb;
  std::vector<typename Tgroup::Telt::Tidx> eInc;
  MyMatrix<Tvert> eBigMat;
};

template<typename Tvert, typename Tgroup>
struct RepartEntry {
  MyMatrix<Tvert> EXT;
  Tgroup TheStab;
  int8_t Position; // -1: lower, 0: barrel, 1: higher
  int iDelaunayOrigin;
  std::vector<AdjRepart<Tvert, Tgroup>> ListAdj;
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
std::vector<RepartEntry<Tvert>> FindRepartitionningInfoNextGeneration(size_t eIdx, DelaunayTesselation<Tvert, Tgroup> const& ListOrbitDelaunay, std::vector<AdjInfo> const& ListInformationsOneFlipping, MyMatrix<T> const& InteriorElement, RecordDualDescOperation<T, Tgroup> & rddo) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n = InteriorElement.rows();
  std::vector<std::vector<Tidx>> ListPermGenList;
  std::vector<MyMatrix<Tvert>> ListMatGens;
  Tgroup PermGRP;
  auto StandardGroupUpdate=[&]() -> void {
    std::vector<Telt> ListGen;
    for (auto & eList : ListPermGenList) {
      Telt x(eList);
      ListGen.push_back(x);
    }
    PermGRP = Tgroup(ListGen);
  };
  StandardGroupUpdate();
  std::vector<MyVector<Tvert>> ListVertices;
  std::unordered_map<MyVector<Tvert>, size_t> ListVertices_rev;
  size_t n_vertices = 1;
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
      break;
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
    MyMatrix<Tvert> LVert = ListOrbitDelaunay.l_dels[TheRec.iDelaunay] * TheRec.eBigMat;
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
      MyMatrix<Tvert> const& EXT = ListOrbitDelaunay.l_dels[TheRec.iDelaunay].obj;
      for (auto & ePermGen : ListOrbitDelaunay.l_dels[TheRec.iDelaunay].GRPlatt.GeneratorsOfGroup()) {
        MyMatrix<Tvert> eBigMat = RepresentVertexPermutation(EXT, EXT, ePermGen);
        MyMatrix<Tvert> eBigMat_new = BigMatI * eBigMat * BigMatR;
        FuncInsertGenerator(eBigMat_new);
      }
      std::vector<Tidx> Linc;
      for (auto & eV : LVert) {
        size_t pos = ListVertices_rev.at(eV);
        Tidx pos_idx = static_cast<Tidx>(pos - 1);
        Linc.push_back(pos_idx);
      }
      TypeOrbitCenter OrbCent{TheRec.iDelaunay, TheRec.eBigMat, false, Linc, LVert};
      ListOrbitCenter.push_back(OrbCent);
    }
  };
  TypeOrbitCenterMin TheRec{eIdx, IdentityMat<Tvert>(n+1)};
  FuncInsertCenter(TheRec);
  while(true) {
    bool IsFinished = true;
    size_t nbCent = ListOrbitCenter.size();
    for (size_t iCent=0; iCent<nbCent; iCent++) {
      TypeOrbitCenter & eEnr = ListOrbitCenter[iCent];
      if (!eEnt.Status) {
        IsFinished = false;
        eEnr.Status = true;
        for (auto & eCase : ListInformationsOneFlipping) {
          if (eEnr.iDelaunay == eCase.iOrb) {
            MyMatrix<Tvert> const& eBigMat = ListOrbitDelaunay.l_dels[eCase.iOrb].l_adj[eCase.i_adj].P;
            MyMatrix<Tvert> eBigMatNew = eBigMat * eEnr.eBigMat;
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
  int nVert = ListVertices.rows();
  MyMatrix<T> TotalListVertices(len, 2+n);
  std::vector<T> LineInterior = GetLineVector(InteriorElement);
  MyVector<Tvert> eV(n);
  for (int iVert=0; iVert<nVert; iVert++) {
    for (int iCol=1; iCol<n; iCol++) {
      T val = UniversalScalarConversion<T,Tvert>(ListVertices(iVert, iCol+1));
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
    Face Linc_face = VectorToFace(eRec.Linc);
    MyVector<T> eFac = FindFacetInequality(TotalListVertices, Linc_face);
    Tgroup TheStab;
    bool Status = false;
    int8_t Position = -1;
    std::vector<AdjRepart<Tvert, Tgroup>> ListAdj;
    int iDelaunayOrigin = eRec.iDelaunay;
    MyMatrix<Tvert> const& eBigMat = eRec.eBigMat;
    MyMatrix<Tvert> const& EXT = eRec.EXT;
    RepartEntry<Tvert, Tgroup> re{EXT, TheStab, Position, iDelaunayOrigin, ListAdj, eBigMat};
    RepartEntryProv rep{eFac, eRec.Linc, Linc_face, Status};
    ListOrbitFacet.push_back(re);
    ListOrbitFacet_prov.push_back(rep);
  }
  auto FuncInsertFacet=[&](MyVector<T> const& eFac) -> AdjReport<Tvert, Tgroup> {
    std::vector<Tidx> Linc;
    Face Linc_face(nVert);
    for (int iVert=0; iVert<nVert; iVert++) {
      if (get_incd_status(iVert, eFac)) {
        Linc.push_back(iVert);
        Linc_face[iVert] = 1;
      }
    }
    int nOrb = ListOrbitFacet.size();
    for (int iOrb=0; iOrb<nOrb; iOrb++) {
      std::optional<Telt> opt = PermGRP.RepresentativeAction(ListOrbitFacet_prov[iOrb].Linc_face, Linc_face);
      if (opt) {
        MyMatrix<Tvert> const& EXT = ListOrbitFacet[iOrb].EXT;
        MyMatrix<Tvert> eBigMat = RepresentVertexPermutation(EXT, EXT, *opt);
        return {iOrb, {}, eBigMat};
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
    std::vector<AdjRepart<Tvert, Tgroup>> ListAdj;
    MyMatrix<Tvert> eBigMat;
    RepartEntry<Tvert, Tgroup> re{EXT, TheStab, Position, iDelaunayOrigin, ListAdj, eBigMat};
    RepartEntryProv rep{eFac, eLinc, Linc_face, Status};
    ListOrbitFacet.push_back(re);
    ListOrbitFacet_prov.push_back(rep);
    int iOrb = nOrb;
    MyMatrix<Tvert> eBigMat = IdentityMat<Tvert>(n+1);
    return {iOrb, {}, eBigMat};
  };
  using Text_int = typename SubsetRankOneSolver<T>::Tint;
  MyMatrix<Text_int> TotalListVertices_int = Get_EXT_int(TotalListVertices);
  while(true) {
    bool IsFinished = true;
    int nOrb = ListOrbitFacet.size();
    for (int iOrb=0; iOrb<nOrb; iOrb++) {
      if (!ListOrbitFacet_prov[iOrb].Status) {
        Face const& Linc_face = ListOrbitFacet_prov[iOrb].Linc_face;
        Tgroup Stab = PermGRP.Stabiliser_OnSets(Linc_face);
        Tgroup TheStab = RenormStabilizer(Stab);
        ListOrbitFacet[iOrb].TheStab = TheStab;
        std::vector<AdjRepart<Tvert, Tgroup>> ListAdj;
        MyMatrix<Tvert> EXT1 = SelectRow(ListVertices, Linc_face);
        MyMatrix<T> EXT2 = UniversalMatrixConversion<T,Tvert>(EXT1);
        vectface vf = DualDescriptionRecord(EXT2, TheStab, rddo);
        FlippingFramework<T> frame(TotalListVertices, TotalListVertices_int, Linc_face, rddo.os);
        for (auto & eFace : vf) {
          Face eInc = frame.FlipFace(eFace);
          MyVector<T> eFac = FindFacetInequality(TotalListVertices, eInc);
          AdjReport<Tvert, Tgroup> eAdj = FuncInsertFacet(eFac);
          std::vector<Tidx> eInc;
          for (size_t iInc=0; iInc<ListOrbitFacet_prov[iOrb].Linc.size(); iInc++) {
            Tidx jInc = ListOrbitFacet_prov[iOrb].Linc[iInc];
            if (get_incd_status(jInc, eFac)) {
              eInc.push_back(iInc);
            }
          }
          eAdj.eInc = eInc;
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
DelaunayTesselation<Tint, Tgroup> FlippingLtype(DelaunayTesselation<Tvert, Tgroup> const& ListOrbitDelaunay, MyMatrix<T> const& InteriorElement, std::vector<AdjInfo> const& ListInformationsOneFlipping, RecordDualDescOperation<T, Tgroup> & rddo) {
  using Tgr = GraphListAdj;
  int n_dels = ListOrbitDelaunay.l_dels.size();
  Tgr Gra(n_dels);
  Face ListMatched(n_dels);
  for (auto & eAI : ListInformationsOneFlipping) {
    ListMatched[eAI.iOrb] = 1;
    int iOrbAdj = ListOrbitDelaunay.l_dels[eAI.iOrb].l_adj[eAI.i_adj].iOrb;
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
  std::vector<std::vector<size_t>> ListConn = ConnectedComponents_set(GR);
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
  std::vector<std::vector<RepartEntry<Tvert>>> ListInfo;
  for (auto &eConn : ListGroupMelt) {
    Tidx eIdx = eConn[0];
    std::vector<RepartEntry<Tvert>> LORB2 = FindRepartitionningInfoNextGeneration(eIdx, ListOrbitDelaunay, ListInformationsOneFlipping, InteriorElement, rddo);
    ListInfo.push_back(LORB2);
  }
  int n_info = ListInfo.size();
  int8_t Position_old = 4, Position_old = 5;
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
        if (NewListOrbitDelaunay[i].Position == ds.Position_old) {
          if (NewListOrbitDelaunay[i].iDelaunay == ds.iDelaunay) {
            return i;
          }
        }
      }
    }
    if (ds.Position == Position_new) {
      for (size_t i=0; i<NewListOrbitDelaunay.size(); i++) {
        if (NewListOrbitDelaunay[i].Position == ds.Position_new) {
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
    Delaunay_Entry<Tvert, Tgroup> del{ListOrbitDelaunay[iDelaunay].obj, ListOrbitDelaunay[iDelaunay].GRP, {}};
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
    MyMatrix<Tvert> const& EXT = l_dels[iOrb].obj;
    ContainerMatrix<Tint> cont(EXT);
    Face f_att(EXT.rows());
    MyMatrix<Tvert> EXTadj = l_dels[NewAdj.iOrb] * NewAdj.P;
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
      for (auto & eAdj : ListOrbitDelaunay.l_dels[iDelaunay].l_adj) {
        int iDelaunayOld = eAdj.iOrb;
        if (ListMatched[iDelaunayOld] == 0) {
          
        } else {
        }
      }
    } else {
    }
  }
  return {l_dels};
}









// clang-format off
#endif  // SRC_LATT_ISODELAUNAYDOMAINS_H_
// clang-format on

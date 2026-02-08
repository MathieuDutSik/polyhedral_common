// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLYDECOMP_DECOMPOSITIONS_H_
#define SRC_POLYDECOMP_DECOMPOSITIONS_H_

// clang-format off
#include "Group.h"
#include "MAT_Matrix.h"
#include "LatticeStabEquiCan.h"
#include "POLY_Kskeletton.h"
#include "POLY_RecursiveDualDesc.h"
#include "Permutation.h"
#include "PolytopeEquiStabInt.h"
#include "triples.h"
#include "Namelist.h"
#include <string>
#include <set>
#include <map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_POLYHEDRAL_DECOMPOSITIONS
#endif

#ifdef TIMINGS
#define TIMINGS_POLYHEDRAL_DECOMPOSITIONS
#endif

// possible strategies for computing isomorphsim

template <typename T, typename Tint_inp, typename Tgroup_inp> struct ConeDesc {
  using Tint = Tint_inp;
  using Tgroup = Tgroup_inp;
  MyMatrix<T> EXT_T;
  MyMatrix<Tint> EXT;
  MyMatrix<T> FAC;
  Face extfac_incd;
  Face facext_incd;
  Tgroup GRP_ext;
  Tgroup GRP_fac;
  std::vector<sing_adj<Tint>> l_sing_adj;
  MyMatrix<Tint> find_matrix(typename Tgroup::Telt const& x, [[maybe_unused]] std::ostream &os) const {
    return FindTransformation(EXT, EXT, x);
  }
};

struct FaceDesc {
  size_t iCone;
  Face f_fac;
};

template <typename T, typename Tint, typename Tidx_value> struct Tent {
  MyMatrix<Tint> M;
  MyMatrix<Tint> Spann;
  MyMatrix<Tint> Qmat;
  WeightMatrix<true, std::vector<Tint>, Tidx_value> WMat;
  FaceDesc fd;
};

template <typename T, typename Tint, typename Tidx_value>
void PrintTent(std::ostream &os, Tent<T, Tint, Tidx_value> const &e_ent) {
  os << "M=\n";
  WriteMatrix(os, e_ent.M);
  os << "Spann=\n";
  WriteMatrix(os, e_ent.Spann);
  os << "Qmat=\n";
  WriteMatrix(os, e_ent.Qmat);
  os << "WMat=\n";
  PrintWeightedMatrix(os, e_ent.WMat);
}

template <typename T>
std::vector<T> f_vsub(const size_t &n_row, const size_t &len) {
  std::vector<T> Vsub(n_row, 23);
  for (size_t i = 0; i < len; i++)
    Vsub[i] = 37;
  return Vsub;
}

Face f_subset(const size_t &n_row, const size_t &len) {
  Face subset(n_row);
  for (size_t i = 0; i < len; i++)
    subset[i] = 1;
  return subset;
}

template <typename Tgroup, typename Tint> struct stab_info {
  Tgroup GRPfull;
  Tgroup GRPres;
  std::vector<std::pair<typename Tgroup::Telt, MyMatrix<Tint>>> ListGenMat;
};

template <typename T, typename Tint, typename Tgroup, typename Tidx_value>
stab_info<Tgroup, Tint> f_stab(const Tent<T, Tint, Tidx_value> &eEnt,
                               std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tgr = GraphBitset;
  Tgroup GRP1 =
      GetStabilizerWeightMatrix<std::vector<Tint>, Tgr, Tgroup, Tidx_value>(
          eEnt.WMat, os);
  MyMatrix<T> Concat_T =
      UniversalMatrixConversion<T, Tint>(Concatenate(eEnt.M, eEnt.Spann));
  Tgroup GRPfull = LinPolytopeIntegral_Stabilizer(Concat_T, GRP1, os);
  Face subset = f_subset(Concat_T.rows(), eEnt.M.rows());
  std::vector<Telt> LGenRed;
  std::vector<std::pair<typename Tgroup::Telt, MyMatrix<Tint>>> ListGenMat;
  for (auto &eGen : GRPfull.GeneratorsOfGroup()) {
    Telt eGenRed = ReduceElementActionFace(eGen, subset);
    LGenRed.push_back(eGenRed);
    MyMatrix<T> eMat_T = FindTransformation(Concat_T, Concat_T, eGen);
    MyMatrix<Tint> eMat = UniversalMatrixConversion<Tint, T>(eMat_T);
    ListGenMat.push_back({eGenRed, eMat});
  }
  Tgroup GRPres(LGenRed, eEnt.M.rows());
  return {GRPfull, GRPres, ListGenMat};
}

template <typename T, typename Tint, typename Tgroup, typename Tidx_value>
std::optional<MyMatrix<Tint>> f_equiv(const Tent<T, Tint, Tidx_value> &eEnt,
                                      const Tent<T, Tint, Tidx_value> &fEnt,
                                      std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tgr = GraphBitset;
  if ((eEnt.M.rows() != fEnt.M.rows()) ||
      (eEnt.Spann.rows() != fEnt.Spann.rows()))
    return {};
  MyMatrix<T> eConcat_T =
      UniversalMatrixConversion<T, Tint>(Concatenate(eEnt.M, eEnt.Spann));
  MyMatrix<T> fConcat_T =
      UniversalMatrixConversion<T, Tint>(Concatenate(fEnt.M, fEnt.Spann));
  std::vector<Tidx> eCanonicReord =
      GetGroupCanonicalizationVector_Kernel<std::vector<Tint>, Tgr, Tidx,
                                            Tidx_value>(eEnt.WMat, os)
          .first;
  std::vector<Tidx> fCanonicReord =
      GetGroupCanonicalizationVector_Kernel<std::vector<Tint>, Tgr, Tidx,
                                            Tidx_value>(fEnt.WMat, os)
          .first;
  // Computing the isomorphism
  using Tfield = typename overlying_field<Tint>::field_type;
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo =
      IsomorphismFromCanonicReord<T, Tfield, Tidx>(
          eConcat_T, fConcat_T, eCanonicReord, fCanonicReord, os);
  if (!IsoInfo) {
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITION
    std::cerr << "Failed isomorphism at the graph level\n";
#endif
    return {};
  }
  Telt ePerm(IsoInfo->first);
  Tgroup GRP1 =
      GetStabilizerWeightMatrix<std::vector<Tint>, Tgr, Tgroup, Tidx_value>(
          eEnt.WMat, os);
  std::optional<MyMatrix<T>> eRes = LinPolytopeIntegral_Isomorphism(
      eConcat_T, fConcat_T, GRP1, ePerm, os);
  if (eRes)
    return UniversalMatrixConversion<Tint, T>(*eRes);
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITION
  std::cerr << "Failed isomorphism at the integral level\n";
#endif
  return {};
}

template <typename T, typename Tint, typename Tgroup, typename Tidx_value>
Tent<T, Tint, Tidx_value>
f_ent(std::vector<ConeDesc<T, Tint, Tgroup>> const &ListCones,
      const MyMatrix<Tint> &G, const FaceDesc &fd, std::ostream &os) {
  using Tidx = typename Tgroup::Telt::Tidx;
  int nbFac = ListCones[fd.iCone].FAC.rows();
  int nbExt = ListCones[fd.iCone].EXT_T.rows();
  Face face_ext = Compute_faceEXT_from_faceFAC(ListCones[fd.iCone].extfac_incd,
                                               nbFac, nbExt, fd.f_fac);
  MyMatrix<Tint> M = SelectRow(ListCones[fd.iCone].EXT, face_ext);
  MyMatrix<Tint> P = M * G;
  MyMatrix<Tint> Spann;
  MyMatrix<Tint> Concat = M;
  MyMatrix<Tint> NSP =
      SublatticeBasisReduction(NullspaceIntMat(TransposedMat(P)), os);
  if (NSP.rows() > 0) {
    MyMatrix<Tint> Gres = -NSP * G * NSP.transpose();
    MyMatrix<T> Gres_T = UniversalMatrixConversion<T, Tint>(Gres);
    MyMatrix<Tint> SHV =
        ExtractInvariantVectorFamilyZbasis<T, Tint>(Gres_T, os);
    Spann = SHV * NSP;
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITION
    std::cerr << "Gres_T=\n";
    WriteMatrix(std::cerr, Gres_T);
    std::cerr << "NSP=\n";
    WriteMatrix(std::cerr, NSP);
    std::cerr << "SHV=\n";
    WriteMatrix(std::cerr, SHV);
    std::cerr << "Spann=\n";
    WriteMatrix(std::cerr, Spann);
#endif
    Concat = Concatenate(M, Spann);
  }
  MyMatrix<Tint> Qmat = GetQmatrix(Concat, os);
  std::vector<Tint> Vsubset = f_vsub<Tint>(Concat.rows(), M.rows());
  std::vector<MyMatrix<Tint>> ListMat{Qmat, G};
  using Tfield = typename overlying_field<Tint>::field_type;
  WeightMatrix<true, std::vector<Tint>, Tidx_value> WMat =
      GetWeightMatrix_ListMat_Vdiag<Tint, Tfield, Tidx, Tidx_value>(
          Concat, ListMat, Vsubset, os);
  WMat.ReorderingSetWeight();
  return {M, Spann, Qmat, std::move(WMat), fd};
}

template <typename T, typename Tint, typename Tgroup, typename Tidx_value>
size_t f_inv(const Tent<T, Tint, Tidx_value> &eEnt) {
  static_assert(std::is_integral<Tidx_value>::value,
                "Tidx_value should be integral");
  const WeightMatrix<true, std::vector<Tint>, Tidx_value> &WMat = eEnt.WMat;
  return std::hash<WeightMatrix<true, std::vector<Tint>, Tidx_value>>()(WMat);
}

template <typename T> MyVector<T> GetFirstNonZeroVector(const MyMatrix<T> &M) {
  size_t n_rows = M.rows();
  for (size_t i_row = 0; i_row < n_rows; i_row++) {
    MyVector<T> V = GetMatrixRow(M, i_row);
    if (!IsZeroVector(V))
      return V;
  }
  std::cerr << "Failed to find a non-zero line\n";
  throw TerminalException{1};
}

template <typename T, typename Tint, typename Tgroup, typename Tidx_value>
void compute_adjacency_structure(
    std::vector<ConeDesc<T, Tint, Tgroup>> & ListCones,
    MyMatrix<Tint> const &G, std::ostream &os) {
  static_assert(std::is_integral<Tidx_value>::value,
                "Tidx_value should be integral");
  std::vector<std::vector<sing_adj<Tint>>> adjacency_information;
  size_t n_domain = ListCones.size();
  struct ent_info {
    size_t i_domain;
    size_t i_adj;
    Face f_ext;
    Tent<T, Tint, Tidx_value> eEnt;
    size_t hash;
  };
  std::vector<ent_info> l_ent_info;
  std::vector<size_t> l_n_orb_adj;
  for (size_t i_domain = 0; i_domain < n_domain; i_domain++) {
    size_t n_fac = ListCones[i_domain].FAC.rows();
    size_t n_ext = ListCones[i_domain].EXT_T.rows();
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITION
    std::cerr << "i_domain=" << i_domain << " n_ext=" << n_ext
              << " n_fac=" << n_fac << "\n";
#endif
    vectface vf = DecomposeOrbitPoint_Full(ListCones[i_domain].GRP_fac);
    size_t i_adj = 0;
    for (auto &eOrb : vf) {
      boost::dynamic_bitset<>::size_type MinVal = eOrb.find_first();
      size_t i_fac = MinVal;
      Face f_ext(n_ext);
      for (size_t i_ext = 0; i_ext < n_ext; i_ext++)
        if (ListCones[i_domain].extfac_incd[i_fac * n_ext + i_ext] == 1)
          f_ext[i_ext] = 1;
      Face f_fac(n_fac);
      f_fac[i_fac] = 1;
      FaceDesc fd{i_domain, f_fac};
      Tent<T, Tint, Tidx_value> eEnt =
          f_ent<T, Tint, Tgroup, Tidx_value>(ListCones, G, fd, os);
      size_t hash = f_inv<T, Tint, Tgroup, Tidx_value>(eEnt);
      ent_info e_ent_info{i_domain, i_adj, f_ext, std::move(eEnt), hash};
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITION
      std::cerr << "  i_domain=" << i_domain << " i_adj=" << i_adj
                << " i_fac=" << i_fac << " |f_ext|=" << f_ext.count()
                << " hash=" << hash << "\n";
      if (i_domain == 5 && f_ext.count() == 201) {
        std::string FileO =
            "DEBUG_" + std::to_string(i_domain) + "_" + std::to_string(i_fac);
        WriteMatrixGAPfile(FileO, e_ent_info.eEnt.M);
      }
#endif
      l_ent_info.emplace_back(std::move(e_ent_info));
      i_adj++;
    }
    l_n_orb_adj.push_back(i_adj);
  }
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITION
  std::cerr << "First block of data built\n";
#endif
  auto get_ent_info = [&](size_t const &i_domain,
                          size_t const &i_adj) -> const ent_info & {
    for (auto &e_ent : l_ent_info)
      if (e_ent.i_domain == i_domain && e_ent.i_adj == i_adj)
        return e_ent;
    std::cerr << "Failed to find the matching entry\n";
    throw TerminalException{1};
  };
  auto get_reverting_transformation = [&](const Tent<T, Tint, Tidx_value> &eEnt)
      -> std::optional<MyMatrix<Tint>> {
    stab_info<Tgroup, Tint> e_stab_info =
        f_stab<T, Tint, Tgroup, Tidx_value>(eEnt, os);
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITION
    std::cerr << "After f_stab call |GRPfull|=" << e_stab_info.GRPfull.size()
              << " |GRPres|=" << e_stab_info.GRPres.size() << "\n";
#endif
    MyVector<Tint> eSpann = GetFirstNonZeroVector(eEnt.Spann);
    for (auto &ePair : e_stab_info.ListGenMat) {
      MyVector<Tint> eSpannImg = ePair.second.transpose() * eSpann;
      if (eSpannImg == -eSpann) {
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITION
        std::cerr << "It is opposite\n";
#endif
        return ePair.second;
      }
      if (eSpannImg != eSpann) {
        std::cerr << "Some inconsistency in the matrix transformation\n";
        throw TerminalException{1};
      }
    }
    return {};
  };
  auto get_mapped = [&](const ent_info &a_ent) -> sing_adj<Tint> {
    for (auto &b_ent : l_ent_info) {
      if (b_ent.hash == a_ent.hash &&
          (a_ent.i_domain != b_ent.i_domain || a_ent.i_adj != b_ent.i_adj)) {
        std::optional<MyMatrix<Tint>> e_equiv =
            f_equiv<T, Tint, Tgroup, Tidx_value>(b_ent.eEnt, a_ent.eEnt, os);
        if (e_equiv) {
          return {b_ent.i_domain, a_ent.f_ext, *e_equiv};
        }
      }
    }
    std::cerr << "a_ent: i_domain=" << a_ent.i_domain
              << " i_adj=" << a_ent.i_adj << "\n";
    std::cerr << "f_ext=" << a_ent.f_ext << " |f_ext|=" << a_ent.f_ext.count()
              << " / " << a_ent.f_ext.size() << "\n";
    std::cerr << "f_fac=" << StringFace(a_ent.eEnt.fd.f_fac)
              << " |f_act|=" << a_ent.eEnt.fd.f_fac.count() << " / "
              << a_ent.eEnt.fd.f_fac.size() << "\n";
    std::cerr << "M=\n";
    WriteMatrix(std::cerr, a_ent.eEnt.M);
    std::cerr << "Spann=\n";
    WriteMatrix(std::cerr, a_ent.eEnt.Spann);
    std::cerr << "Qmat=\n";
    WriteMatrix(std::cerr, a_ent.eEnt.Qmat);
    std::cerr << "Failed to find a matching entry in get_mapped\n";
    throw TerminalException{1};
  };
  auto get_sing_adj = [&](size_t const &i_domain,
                          size_t const &i_adj) -> sing_adj<Tint> {
    const ent_info &e_ent = get_ent_info(i_domain, i_adj);
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITION
    std::cerr << "get_sing_adj i_domain=" << i_domain << " i_adj=" << i_adj
              << "\n";
#endif
    std::optional<MyMatrix<Tint>> trans_opt =
        get_reverting_transformation(e_ent.eEnt);
    if (trans_opt) {
      return {i_domain, e_ent.f_ext, *trans_opt};
    }
    return get_mapped(e_ent);
  };
  for (size_t i_domain = 0; i_domain < n_domain; i_domain++) {
    std::vector<sing_adj<Tint>> l_sing_adj;
    for (size_t i_adj = 0; i_adj < l_n_orb_adj[i_domain]; i_adj++) {
      l_sing_adj.push_back(get_sing_adj(i_domain, i_adj));
    }
    ListCones[i_domain].l_sing_adj = l_sing_adj;
  }
}

template <typename T, typename Tint, typename Tgroup, typename Tidx_value>
std::vector<std::vector<FaceDesc>> Compute_ListListDomain_strategy2(
    std::vector<ConeDesc<T, Tint, Tgroup>> const &ListCones,
    MyMatrix<Tint> const &G, int TheLev, std::ostream &os) {
  static_assert(std::is_integral<Tidx_value>::value,
                "Tidx_value should be integral");
  std::vector<FaceDesc> ListDomain;
  for (size_t i = 0; i < ListCones.size(); i++) {
    size_t len = ListCones[i].FAC.rows();
    Face f_fac(len);
    FaceDesc fd = {i, f_fac};
    ListDomain.push_back(fd);
  }
  std::vector<std::vector<FaceDesc>> ListListDomain;
  ListListDomain.push_back(ListDomain);
  //
  for (int i = 1; i < TheLev; i++) {
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITION
    os << "i=" << i << "\n";
#endif
    std::vector<std::pair<size_t, Tent<T, Tint, Tidx_value>>> NewListCand;
#ifdef SANITY_CHECK_POLYEDRAL_DECOMPOSITION
    using Tface =
      std::pair<std::vector<triple<Tint>>, std::vector<MyMatrix<Tint>>>;
    std::vector<Tface> list_face;
    size_t n_equiv_found = 0;
    size_t dim = ListCones[0].FAC.cols();
    auto f_insert_face = [&](const FaceDesc &fd_A) -> void {
      const ConeDesc<T, Tint, Tgroup> &eC = ListCones[fd_A.iCone];
      Face f_ext = Compute_faceEXT_from_faceFAC(eC.extfac_incd, eC.FAC.rows(),
                                                eC.EXT_T.rows(), fd_A.f_fac);
      triple<Tint> ef_A_pre{fd_A.iCone, f_ext, IdentityMat<Tint>(dim)};
      triple<Tint> ef_A = canonicalize_triple(ListCones, ef_A_pre, os);
      for (auto &eOrbit : list_face) {
        for (auto &ef_B : eOrbit.first) {
          std::optional<MyMatrix<Tint>> equiv_opt =
            test_equiv_triple(ListCones, ef_A, ef_B, os);
          if (equiv_opt) {
            n_equiv_found++;
            Tent<T, Tint, Tidx_value> ent_A =
                f_ent<T, Tint, Tgroup, Tidx_value>(ListCones, G, fd_A, os);
            const ConeDesc<T, Tint, Tgroup> &eC_B = ListCones[ef_B.iCone];
            Face f_fac = Compute_faceFAC_from_faceEXT(
                eC_B.extfac_incd, eC_B.FAC.rows(), eC_B.EXT_T.rows(), ef_B.f_ext);
            FaceDesc fd_B{ef_B.iCone, f_fac};
            Tent<T, Tint, Tidx_value> ent_B =
                f_ent<T, Tint, Tgroup, Tidx_value>(ListCones, G, fd_B, os);
            bool test = false;
            if (f_equiv<T, Tint, Tgroup, Tidx_value>(ent_A, ent_B, os))
              test = true;
            std::cerr << "------------------ FOUND A BUG ---------------\n";
            std::cerr << "ent_A=\n";
            PrintTent(std::cerr, ent_A);
            std::cerr << "ent_B=\n";
            PrintTent(std::cerr, ent_B);
            std::cerr << "test=" << test << "\n";
            return;
          }
        }
      }
      list_face.push_back(get_spanning_list_triple(ListCones, ef_A, os));
    };
#endif
    auto f_insert = [&](Tent<T, Tint, Tidx_value> &&eEnt) -> void {
      size_t e_inv = f_inv<T, Tint, Tgroup, Tidx_value>(eEnt);
      for (auto &eP : NewListCand) {
        if (eP.first == e_inv) {
          std::optional<MyMatrix<Tint>> eEquiv =
              f_equiv<T, Tint, Tgroup, Tidx_value>(eP.second, eEnt, os);
          if (eEquiv)
            return;
        }
      }
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITION
      f_insert_face(eEnt.fd);
#endif
      std::pair<size_t, Tent<T, Tint, Tidx_value>> e_pair{e_inv,
                                                          std::move(eEnt)};
      NewListCand.emplace_back(std::move(e_pair));
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITION
      std::cerr << "Now |NewListCand|=" << NewListCand.size() << "\n";
#endif
    };
    for (auto &eDomain : ListListDomain[i - 1]) {
      Tent<T, Tint, Tidx_value> eEnt =
          f_ent<T, Tint, Tgroup, Tidx_value>(ListCones, G, eDomain, os);
      size_t iCone = eDomain.iCone;
      Tgroup StabFace_fac =
          ListCones[iCone].GRP_fac.Stabilizer_OnSets(eDomain.f_fac);
      int RankFace = i - 1;
      vectface ListFace = SPAN_face_ExtremeRays(
          eDomain.f_fac, StabFace_fac, RankFace, ListCones[iCone].extfac_incd,
          ListCones[iCone].FAC, ListCones[iCone].EXT_T, os);
      for (auto &eFace_fac : ListFace) {
        FaceDesc fdn{iCone, eFace_fac};
        Tent<T, Tint, Tidx_value> fEnt =
            f_ent<T, Tint, Tgroup, Tidx_value>(ListCones, G, fdn, os);
        f_insert(std::move(fEnt));
      }
    }
    std::vector<FaceDesc> NewListDomain;
    for (auto &eEnt : NewListCand) {
      NewListDomain.push_back(eEnt.second.fd);
    }
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITION
    if (list_face.size() != NewListDomain.size()) {
      std::cerr << "Inconsistent sizes\n";
      std::cerr << "|list_face|=" << list_face.size()
                << " |NewListDomain|=" << NewListDomain.size()
                << " n_equiv_found=" << n_equiv_found << "\n";
      throw TerminalException{1};
    }
    std::cerr << "i=" << i << " |NewListDomain|=" << NewListDomain.size()
              << "\n";
#endif
    ListListDomain.emplace_back(std::move(NewListDomain));
  }
  return ListListDomain;
}

template <typename T, typename Tint, typename Tgroup, typename Tidx_value>
std::vector<std::vector<FaceDesc>> Compute_ListListDomain_strategy1(
    std::vector<ConeDesc<T, Tint, Tgroup>> const &ListCones,
    MyMatrix<Tint> const &G, int TheLev, std::ostream &os) {
  std::vector<std::vector<FaceDesc>> ListListDomain;
  size_t n_col = G.rows();
  using Tface =
      std::pair<std::vector<triple<Tint>>, std::vector<MyMatrix<Tint>>>;
  std::vector<std::vector<Tface>> list_list_face;
  for (int iLev = 1; iLev <= TheLev; iLev++) {
    std::cerr << "iLev=" << iLev << "\n";
    std::vector<Tface> list_face;
    auto f_insert = [&](const triple<Tint> &ef_A) -> void {
      for (auto &eOrbit : list_face) {
        for (auto &ef_B : eOrbit.first) {
          std::optional<MyMatrix<Tint>> equiv_opt =
            test_equiv_triple(ListCones, ef_A, ef_B, os);
          if (equiv_opt) {
            return;
          }
        }
      }
      list_face.push_back(get_spanning_list_triple(ListCones, ef_A, os));
    };
    if (iLev == 1) {
      for (size_t iCone = 0; iCone < ListCones.size(); iCone++) {
        vectface vf = DecomposeOrbitPoint_Full(ListCones[iCone].GRP_ext);
        size_t nExt = ListCones[iCone].EXT_T.rows();
        for (auto &eOrb : vf) {
          Face f_ext(nExt);
          boost::dynamic_bitset<>::size_type eVal = eOrb.find_first();
          size_t iExt = eVal;
          f_ext[iExt] = 1;
          triple<Tint> e_ent_pre{iCone, f_ext, IdentityMat<Tint>(n_col)};
          triple<Tint> e_ent = canonicalize_triple(ListCones, e_ent_pre, os);
          f_insert(e_ent);
        }
      }
    } else {
      for (auto &eOrbit : list_list_face[iLev - 2]) {
        for (auto &eRepr : eOrbit.first) {
          const ConeDesc<T, Tint, Tgroup> &eC = ListCones[eRepr.iCone];
          Tgroup StabFace_ext = eC.GRP_ext.Stabilizer_OnSets(eRepr.f_ext);
          int RankFace_ext = iLev - 1;
          vectface vf =
              SPAN_face_ExtremeRays(eRepr.f_ext, StabFace_ext, RankFace_ext,
                                    eC.facext_incd, eC.EXT_T, eC.FAC, os);
          for (auto &f_ext_new : vf) {
            triple<Tint> e_ent_pre{eRepr.iCone, f_ext_new,
                                   IdentityMat<Tint>(n_col)};
            triple<Tint> e_ent = canonicalize_triple(ListCones, e_ent_pre, os);
            f_insert(e_ent);
          }
        }
      }
    }
    std::vector<FaceDesc> ListDomain;
    for (auto &eOrbit : list_face) {
      triple<Tint> e_ent = eOrbit.first[0];
      const ConeDesc<T, Tint, Tgroup> &eC = ListCones[e_ent.iCone];
      Face f_fac = Compute_faceFAC_from_faceEXT(eC.extfac_incd, eC.FAC.rows(),
                                                eC.EXT_T.rows(), e_ent.f_ext);
      FaceDesc fd{e_ent.iCone, f_fac};
      ListDomain.push_back(fd);
    }
    ListListDomain.push_back(ListDomain);
    list_list_face.push_back(list_face);
  }
  return ListListDomain;
}

FullNamelist NAMELIST_TestUnionCones() {
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, std::string> ListStringValues1_doc;
  std::map<std::string, std::string> ListIntValues1_doc;
  std::map<std::string, std::string> ListBoolValues1_doc;
  ListStringValues1_doc["FileI"] = "Default: unset.ext\n\
The input file for the list of polyhedral cone (EXT first, then FAC)";
  ListStringValues1_doc["FileO"] = "Default: stderr\n\
The output file. stderr means writing to std::cerr, stdout means for std::cout and otherwise to FileO";
  ListBoolValues1_doc["TestPairwiseIntersection"] = "Default: true\n\
Whether to test for pairwise intersection of cones. That test can be expensive and so it can be disabled";
  ListBoolValues1_doc["BreakConnectedComponents"] = "Default: false\n\
It can be that the cones are union of disjoint cones. In that first process this";
  SingleBlock BlockPROC;
  BlockPROC.setListIntValues_doc(ListIntValues1_doc);
  BlockPROC.setListBoolValues_doc(ListBoolValues1_doc);
  BlockPROC.setListStringValues_doc(ListStringValues1_doc);
  ListBlock["PROC"] = BlockPROC;
  // Merging all data
  return FullNamelist(ListBlock);
}

template <typename T> struct ConeSimpDesc {
  MyMatrix<T> EXT;
  MyMatrix<T> FAC;
};

template <typename T>
std::optional<ConeSimpDesc<T>>
TestPolyhedralPartition(bool const &TestPairwiseIntersection,
                        std::vector<ConeSimpDesc<T>> const &l_cone,
                        std::ostream &os) {
  size_t n_cone = l_cone.size();
  int dim = l_cone[0].FAC.cols();
  HumanTime time;
  if (TestPairwiseIntersection) {
    for (size_t i_cone = 0; i_cone < n_cone; i_cone++) {
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITIONS
      std::cerr << "i_cone=" << i_cone << " / " << n_cone << "\n";
#endif
#ifdef TIMINGS_POLYHEDRAL_DECOMPOSITIONS
      std::cerr << "|DECOMP: IsFullDimensional tests|=" << time << "\n";
#endif
      for (size_t j_cone = i_cone + 1; j_cone < n_cone; j_cone++) {
        MyMatrix<T> FACtot =
            Concatenate(l_cone[i_cone].FAC, l_cone[j_cone].FAC);
        if (IsFullDimensional_V1(FACtot, os)) {
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITIONS
          std::cerr << "Cone i_cone=" << i_cone << " and j_cone=" << j_cone
                    << " are overlapping\n";
#endif
          return {};
        }
      }
    }
#ifdef TIMINGS_POLYHEDRAL_DECOMPOSITIONS
    std::cerr << "|DECOMP: pairwise intersection tests|=" << time << "\n";
#endif
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITIONS
    std::cerr << "Passing the pairwise intersection test\n";
#endif
  } else {
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITIONS
    std::cerr << "Not doing the pairwise intersection test\n";
#endif
  }
  std::map<MyVector<T>, int> MapEXT;
  using Tent = std::pair<size_t, int>;
  for (size_t i_cone = 0; i_cone < n_cone; i_cone++) {
    MyMatrix<T> const &EXT = l_cone[i_cone].EXT;
    int n_ext = EXT.rows();
    for (int i_ext = 0; i_ext < n_ext; i_ext++) {
      MyVector<T> V = GetMatrixRow(EXT, i_ext);
      MapEXT[V] = 0;
    }
  }
  int n_ext_tot = MapEXT.size();
#ifdef TIMINGS_POLYHEDRAL_DECOMPOSITIONS
  std::cerr << "|DECOMP: MapEXT 1|=" << time << "\n";
#endif
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITIONS
  std::cerr << "DECOMP: n_ext_tot=" << n_ext_tot << "\n";
#endif
  MyMatrix<T> EXTtot(n_ext_tot, dim);
  int pos = 0;
  for (auto &kv : MapEXT) {
    MyVector<T> const &eEXT = kv.first;
    for (int i = 0; i < dim; i++) {
      EXTtot(pos, i) = eEXT(i);
    }
    pos++;
  }
#ifdef TIMINGS_POLYHEDRAL_DECOMPOSITIONS
  std::cerr << "|DECOMP: EXTtot|=" << time << "\n";
#endif
  for (int i_ext_tot = 0; i_ext_tot < n_ext_tot; i_ext_tot++) {
    MyVector<T> eEXT = GetMatrixRow(EXTtot, i_ext_tot);
    MapEXT[eEXT] = i_ext_tot;
  }
#ifdef TIMINGS_POLYHEDRAL_DECOMPOSITIONS
  std::cerr << "|DECOMP: MapEXT 2|=" << time << "\n";
#endif
  // We match the facets in order to find which ones are
  // contained into just one and so get the facet of our cones.
  std::map<Face, std::vector<Tent>> FACmult;
  for (size_t i_cone = 0; i_cone < n_cone; i_cone++) {
    MyMatrix<T> const &FAC = l_cone[i_cone].FAC;
    MyMatrix<T> const &EXT = l_cone[i_cone].EXT;
    int n_fac = FAC.rows();
    int n_ext = EXT.rows();
    for (int i_fac = 0; i_fac < n_fac; i_fac++) {
      Face f(n_ext_tot);
      for (int i_ext = 0; i_ext < n_ext; i_ext++) {
        MyVector<T> eEXT = GetMatrixRow(EXT, i_ext);
        int i_ext_tot = MapEXT[eEXT];
        T sum = 0;
        for (int i = 0; i < dim; i++) {
          sum += eEXT(i) * FAC(i_fac, i);
        }
        if (sum == 0) {
          f[i_ext_tot] = 1;
        }
      }
      Tent eEnt{i_cone, i_fac};
      FACmult[f].push_back(eEnt);
    }
  }
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITIONS
  std::cerr << "DECOMP: |FACmult|=" << FACmult.size() << "\n";
#endif
#ifdef TIMINGS_POLYHEDRAL_DECOMPOSITIONS
  std::cerr << "|DECOMP: FACmult|=" << time << "\n";
#endif
  // We check that the facets matched only once have all
  // the vertices on the right side.
  // At the same time, we build the set of facets (with
  // no duplication)
  std::set<MyVector<T>> SetFAC;
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITIONS
  size_t n_size1 = 0;
#endif
  auto is_matching_facet = [&](MyVector<T> const &eFAC) -> bool {
    for (int i_ext_tot = 0; i_ext_tot < n_ext_tot; i_ext_tot++) {
      T sum(0);
      for (int i = 0; i < dim; i++) {
        sum += eFAC(i) * EXTtot(i_ext_tot, i);
      }
      if (sum < 0)
        return false;
    }
    return true;
  };
  for (auto &kv : FACmult) {
    std::vector<Tent> const &eList = kv.second;
    if (eList.size() != 1 && eList.size() != 2) {
      std::cerr << "The length of the list is not 1 or 2\n";
      return {};
    }
    if (eList.size() == 1) {
      size_t i_cone = eList[0].first;
      int i_fac = eList[0].second;
      MyVector<T> eFAC = GetMatrixRow(l_cone[i_cone].FAC, i_fac);
      if (!is_matching_facet(eFAC)) {
        std::cerr << "The facet has some violating vertices\n";
        return {};
      }
      SetFAC.insert(eFAC);
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITIONS
      n_size1++;
#endif
    }
  }
#ifdef TIMINGS_POLYHEDRAL_DECOMPOSITIONS
  std::cerr << "|DECOMP: SetFAC|=" << time << "\n";
#endif
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITIONS
  std::cerr << "DECOMP: SetFAC n_size1=" << n_size1
            << " |SetFAC|=" << SetFAC.size() << "\n";
#endif
  // Building the FAC matrix
  int n_fac_final = SetFAC.size();
  std::cerr << "n_fac_final=" << n_fac_final << "\n";
  MyMatrix<T> FACfinal(n_fac_final, dim);
  std::vector<MyVector<T>> ListFACfinal;
  int pos_fac = 0;
  for (auto &eFAC : SetFAC) {
    for (int i = 0; i < dim; i++)
      FACfinal(pos_fac, i) = eFAC(i);
    ListFACfinal.push_back(eFAC);
    pos_fac++;
  }
#ifdef TIMINGS_POLYHEDRAL_DECOMPOSITIONS
  std::cerr << "|DECOMP: FACfinal|=" << time << "\n";
#endif
  // Selecting the vertices of the right rank.
  std::vector<MyVector<T>> ListEXT;
  size_t min_len = dim - 1;
  for (int i_ext_tot = 0; i_ext_tot < n_ext_tot; i_ext_tot++) {
    std::vector<MyVector<T>> ListIncd;
    MyVector<T> eEXT = GetMatrixRow(EXTtot, i_ext_tot);
    for (int i_fac_final = 0; i_fac_final < n_fac_final; i_fac_final++) {
      MyVector<T> const &eFAC = ListFACfinal[i_fac_final];
      T sum = 0;
      for (int i = 0; i < dim; i++) {
        sum += eFAC(i) * eEXT(i);
      }
      if (sum == 0) {
        ListIncd.push_back(eFAC);
      }
    }
    if (ListIncd.size() >= min_len) {
      MyMatrix<T> FACincd = MatrixFromVectorFamily(ListIncd);
      if (RankMat(FACincd) == dim - 1) {
        ListEXT.push_back(eEXT);
      }
    }
  }
  MyMatrix<T> EXTfinal = MatrixFromVectorFamily(ListEXT);
#ifdef TIMINGS_POLYHEDRAL_DECOMPOSITIONS
  std::cerr << "|DECOMP: EXTfinal|=" << time << "\n";
#endif
#ifdef DEBUG_POLYHEDRAL_DECOMPOSITIONS
  std::cerr << "DECOMP: EXTfinal |ListEXT|=" << ListEXT.size() << "\n";
#endif
  ConeSimpDesc<T> cone{EXTfinal, FACfinal};
  return cone;
}

template <typename T>
std::vector<std::vector<size_t>>
ConnectedComponentsPolyhedral(std::vector<ConeSimpDesc<T>> const &l_cone) {
  std::map<MyVector<T>, size_t> MapEXT;
  int dim = 0;
  for (auto &e_cone : l_cone) {
    int n_ext = e_cone.EXT.rows();
    dim = e_cone.EXT.cols();
    for (int i_ext = 0; i_ext < n_ext; i_ext++) {
      MyVector<T> eEXT = GetMatrixRow(e_cone.EXT, i_ext);
      MapEXT[eEXT] = 0;
    }
  }
  size_t pos = 0;
  for (auto &kv : MapEXT) {
    kv.second = pos;
    pos++;
  }
  size_t n_ext_tot = pos;
  std::cerr << "Building the full set of vertices n_ext_tot=" << n_ext_tot
            << "\n";
  //
  std::map<Face, std::vector<size_t>> MapFace_idx;
  size_t n_cone = l_cone.size();
  for (size_t i_cone = 0; i_cone < n_cone; i_cone++) {
    ConeSimpDesc<T> const &e_cone = l_cone.at(i_cone);
    int n_fac = e_cone.FAC.rows();
    int n_ext = e_cone.EXT.rows();
    for (int i_fac = 0; i_fac < n_fac; i_fac++) {
      MyVector<T> eFAC = GetMatrixRow(e_cone.FAC, i_fac);
      Face f(n_ext_tot);
      for (int i_ext = 0; i_ext < n_ext; i_ext++) {
        MyVector<T> eEXT = GetMatrixRow(e_cone.EXT, i_ext);
        T scal = 0;
        for (int i = 0; i < dim; i++)
          scal += eEXT(i) * eFAC(i);
        if (scal == 0) {
          size_t i_ext_tot = MapEXT.at(eEXT);
          f[i_ext_tot] = 1;
        }
      }
      MapFace_idx[f].push_back(i_cone);
    }
  }
  std::cerr << "The map has been built\n";
  //
  using Tgr = GraphListAdj;
  std::vector<std::pair<size_t, size_t>> l_pair;
  for (auto &kv : MapFace_idx) {
    std::vector<size_t> const &eList = kv.second;
    if (eList.size() != 1 && eList.size() != 2) {
      std::cerr << "Wrong size for eList\n";
      throw TerminalException{1};
    }
    if (eList.size() == 2) {
      size_t a = eList[0];
      size_t b = eList[1];
      std::pair<size_t, size_t> epair{a, b};
      l_pair.push_back(epair);
    }
  }
  Tgr eGR(l_pair, n_cone);
  std::vector<std::vector<size_t>> vect_cone = ConnectedComponents_set(eGR);
  return vect_cone;
}

// clang-format off
#endif  // SRC_POLYDECOMP_DECOMPOSITIONS_H_
// clang-format on

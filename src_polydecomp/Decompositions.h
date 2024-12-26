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

template <typename T, typename Tint, typename Tgroup> struct ConeDesc {
  MyMatrix<T> EXT;
  MyMatrix<Tint> EXT_i;
  MyMatrix<T> FAC;
  Face extfac_incd;
  Face facext_incd;
  Tgroup GRP_ext;
  Tgroup GRP_fac;
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
  Tgroup GRPfull = LinPolytopeIntegral_Stabilizer_Method8(Concat_T, GRP1, os);
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
          eConcat_T, fConcat_T, eCanonicReord, fCanonicReord);
  if (!IsoInfo) {
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
    std::cerr << "Failed isomorphism at the graph level\n";
#endif
    return {};
  }
  Telt ePerm(IsoInfo->first);
  Tgroup GRP1 =
      GetStabilizerWeightMatrix<std::vector<Tint>, Tgr, Tgroup, Tidx_value>(
          eEnt.WMat, os);
  std::optional<MyMatrix<T>> eRes = LinPolytopeIntegral_Isomorphism_Method8(
      eConcat_T, fConcat_T, GRP1, ePerm, os);
  if (eRes)
    return UniversalMatrixConversion<Tint, T>(*eRes);
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
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
  int nbExt = ListCones[fd.iCone].EXT.rows();
  Face face_ext = Compute_faceEXT_from_faceFAC(ListCones[fd.iCone].extfac_incd,
                                               nbFac, nbExt, fd.f_fac);
  MyMatrix<Tint> M = SelectRow(ListCones[fd.iCone].EXT_i, face_ext);
  MyMatrix<Tint> P = M * G;
  MyMatrix<Tint> Spann;
  MyMatrix<Tint> Concat = M;
  MyMatrix<Tint> NSP = NullspaceIntMat(TransposedMat(P));
  if (NSP.rows() > 0) {
    MyMatrix<Tint> Gres = -NSP * G * NSP.transpose();
    MyMatrix<T> Gres_T = UniversalMatrixConversion<T, Tint>(Gres);
    MyMatrix<Tint> SHV =
        ExtractInvariantVectorFamilyZbasis<T, Tint>(Gres_T, os);
    Spann = SHV * NSP;
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
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

template <typename Tint> struct sing_adj {
  size_t jCone;
  Face f_ext;
  MyMatrix<Tint> eMat;
};

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
std::vector<std::vector<sing_adj<Tint>>> compute_adjacency_structure(
    std::vector<ConeDesc<T, Tint, Tgroup>> const &ListCones,
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
    size_t n_ext = ListCones[i_domain].EXT.rows();
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
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
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
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
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
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
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
    std::cerr << "After f_stab call |GRPfull|=" << e_stab_info.GRPfull.size()
              << " |GRPres|=" << e_stab_info.GRPres.size() << "\n";
#endif
    MyVector<Tint> eSpann = GetFirstNonZeroVector(eEnt.Spann);
    for (auto &ePair : e_stab_info.ListGenMat) {
      MyVector<Tint> eSpannImg = ePair.second.transpose() * eSpann;
      if (eSpannImg == -eSpann) {
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
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
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
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
  std::vector<std::vector<sing_adj<Tint>>> ll_adj_struct;
  for (size_t i_domain = 0; i_domain < n_domain; i_domain++) {
    std::vector<sing_adj<Tint>> l_adj_struct;
    for (size_t i_adj = 0; i_adj < l_n_orb_adj[i_domain]; i_adj++) {
      l_adj_struct.push_back(get_sing_adj(i_domain, i_adj));
    }
    ll_adj_struct.push_back(l_adj_struct);
  }
  return ll_adj_struct;
}

// A face of the cellular complex is determined by all the ways in which
// f_ext is the subset of the corresponding face
template <typename Tint> struct ent_face {
  size_t iCone;
  Face f_ext;
  MyMatrix<Tint> eMat;
};

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
test_equiv_ent_face(std::vector<ConeDesc<T, Tint, Tgroup>> const &ListCones,
                    ent_face<Tint> const &ef1, ent_face<Tint> const &ef2) {
  using Telt = typename Tgroup::Telt;
  if (ef1.iCone != ef2.iCone)
    return {};
  size_t iC = ef1.iCone;
  const ConeDesc<T, Tint, Tgroup> &eC = ListCones[iC];
  std::optional<Telt> test =
      eC.GRP_ext.RepresentativeAction_OnSets(ef1.f_ext, ef2.f_ext);
  if (!test)
    return {};
  MyMatrix<T> eMat_T = FindTransformation(eC.EXT, eC.EXT, *test);
  MyMatrix<Tint> eMat = UniversalMatrixConversion<Tint, T>(eMat_T);
  return Inverse(ef1.eMat) * eMat * ef2.eMat;
}

/*
  Generate the list of entries in the face and the list of stabilizer generators
 */
template <typename T, typename Tint, typename Tgroup>
std::pair<std::vector<ent_face<Tint>>, std::vector<MyMatrix<Tint>>>
get_spanning_list_ent_face(
    std::vector<ConeDesc<T, Tint, Tgroup>> const &ListCones,
    std::vector<std::vector<sing_adj<Tint>>> const &ll_sing_adj,
    const ent_face<Tint> &ef_input) {
  //  std::cerr << "Beginning of get_spanning_list_ent_face\n";
  using Telt = typename Tgroup::Telt;
  std::vector<MyMatrix<Tint>> ListMatrGen;
  std::set<MyVector<Tint>> set_EXT;
  // That value of dim should be overwritten later
  for (auto &ePt : FaceToVector<int>(ef_input.f_ext)) {
    MyVector<Tint> V = GetMatrixRow(ListCones[ef_input.iCone].EXT_i, ePt);
    // In GAP Vimg = V A and in transpose we get Vimg^T = A^T V^T
    MyVector<Tint> Vimg = ef_input.eMat.transpose() * V;
    set_EXT.insert(Vimg);
  }
  auto f_insert_generator = [&](const MyMatrix<Tint> &eMatrGen) -> void {
    ListMatrGen.push_back(eMatrGen);
    for (auto &eV : set_EXT) {
      MyVector<Tint> Vimg = eMatrGen.transpose() * eV;
      if (set_EXT.count(Vimg) != 1) {
        std::cerr
            << "Error: The generator does not preserve the face globally\n";
        throw TerminalException{1};
      }
    }
  };
  std::vector<ent_face<Tint>> l_ent_face;
  auto f_insert = [&](const ent_face<Tint> &ef_A) -> void {
    for (const auto &ef_B : l_ent_face) {
      std::optional<MyMatrix<Tint>> equiv_opt =
          test_equiv_ent_face(ListCones, ef_A, ef_B);
      if (equiv_opt) {
        f_insert_generator(*equiv_opt);
        return;
      }
    }
    l_ent_face.push_back(ef_A);
    const ConeDesc<T, Tint, Tgroup> &uC = ListCones[ef_A.iCone];
    Tgroup stab = uC.GRP_ext.Stabilizer_OnSets(ef_A.f_ext);
    MyMatrix<Tint> eInv = Inverse(ef_A.eMat);
    for (auto &eGen : stab.GeneratorsOfGroup()) {
      MyMatrix<T> eMatGen_T = FindTransformation(uC.EXT, uC.EXT, eGen);
      MyMatrix<Tint> eMatGen = UniversalMatrixConversion<Tint, T>(eMatGen_T);
      MyMatrix<Tint> TransGen = eInv * eMatGen * ef_A.eMat;
      f_insert_generator(TransGen);
    }
  };
  f_insert(ef_input);
  size_t curr_pos = 0;
  while (true) {
    size_t len = l_ent_face.size();
    if (curr_pos == len)
      break;
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
    std::cerr << "curr_pos=" << curr_pos << " len=" << len << "\n";
#endif
    for (size_t i = curr_pos; i < len; i++) {
      // We cannot use const& for ent_face for mysterious
      // reasons (3 hours debugging)
      ent_face<Tint> ef = l_ent_face[i];
      const ConeDesc<T, Tint, Tgroup> &eC = ListCones[ef.iCone];
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
      std::cerr << "i=" << i << " iCone=" << ef.iCone
                << " |ll_sing_adj[ef.iCone]|=" << ll_sing_adj[ef.iCone].size()
                << "\n";
      size_t n_facet = 0;
#endif
      for (auto &e_sing_adj : ll_sing_adj[ef.iCone]) {
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
        std::cerr << "|ef.f_ext|=" << SignatureFace(ef.f_ext) << "\n";
#endif
        std::vector<std::pair<Face, Telt>> l_pair =
            FindContainingOrbit(eC.GRP_ext, e_sing_adj.f_ext, ef.f_ext);
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
        vectface vfo =
            OrbitFace(e_sing_adj.f_ext, eC.GRP_ext.GeneratorsOfGroup());
        n_facet += vfo.size();
        Tgroup stab = eC.GRP_ext.Stabilizer_OnSets(ef.f_ext);
        std::cerr << "|vfo|=" << vfo.size() << " |stab|=" << stab.size()
                  << "\n";
        size_t n_match = 0;
        vectface vfcont(vfo.get_n());
        for (auto &eFace : vfo) {
          if (is_subset(ef.f_ext, eFace)) {
            std::cerr << "V=" << StringFace(eFace) << "\n";
            vfcont.push_back(eFace);
            n_match++;
          }
        }
        vectface vfs = OrbitSplittingSet(vfcont, stab);
        std::cerr << "|l_pair|=" << l_pair.size()
                  << " |e_sing_adj.f_ext|=" << SignatureFace(e_sing_adj.f_ext)
                  << " |ef.f_ext|=" << SignatureFace(ef.f_ext)
                  << " n_match=" << n_match << " |vfs|=" << vfs.size() << "\n";
        if (vfs.size() != l_pair.size()) {
          std::cerr << "We have a size error\n";
          throw TerminalException{1};
        }
#endif
        for (auto &e_pair : l_pair) {
          MyMatrix<T> eMat1_T =
              FindTransformation(eC.EXT, eC.EXT, e_pair.second);
          MyMatrix<Tint> eMat1 = UniversalMatrixConversion<Tint, T>(eMat1_T);
          size_t jCone = e_sing_adj.jCone;
          const ConeDesc<T, Tint, Tgroup> &fC = ListCones[jCone];
          MyMatrix<Tint> eMatAdj = e_sing_adj.eMat * eMat1 * ef.eMat;
          MyMatrix<Tint> EXTimg = fC.EXT_i * eMatAdj;
          ContainerMatrix<Tint> Cont(EXTimg);
          Face faceNew(EXTimg.rows());
          for (auto &e_line : set_EXT) {
            std::optional<size_t> opt = Cont.GetIdx_v(e_line);
            if (!opt) {
              std::cerr << "The vector is not in the image. Clear bug\n";
              throw TerminalException{1};
            }
            size_t idx = *opt;
            faceNew[idx] = 1;
          }
          ent_face<Tint> efNew{jCone, faceNew, eMatAdj};
          f_insert(efNew);
        }
      }
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
      std::cerr << "iCone=" << ef.iCone
                << " |ef.f_ext|=" << SignatureFace(ef.f_ext)
                << " n_facet=" << n_facet << "\n";
#endif
    }
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
    std::cerr << "Now |l_ent_face|=" << l_ent_face.size() << "\n";
#endif
    curr_pos = len;
  }
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
  std::cerr << "|l_ent_face|=" << l_ent_face.size()
            << " |ListMatrGen|=" << ListMatrGen.size() << "\n";
  std::cerr << "l_ent_face.iCon =";
  for (auto &e_ent : l_ent_face)
    std::cerr << " " << e_ent.iCone;
  std::cerr << "\n";
#endif
  return {l_ent_face, ListMatrGen};
}

template <typename T, typename Tint, typename Tgroup, typename Tidx_value>
std::vector<std::vector<FaceDesc>> Compute_ListListDomain_strategy2(
    std::vector<ConeDesc<T, Tint, Tgroup>> const &ListCones,
    MyMatrix<Tint> const &G,
    std::vector<std::vector<sing_adj<Tint>>> const &ll_sing_adj, int TheLev,
    std::ostream &os) {
  static_assert(std::is_integral<Tidx_value>::value,
                "Tidx_value should be integral");
  std::vector<FaceDesc> ListDomain;
  size_t dim;
  for (size_t i = 0; i < ListCones.size(); i++) {
    size_t len = ListCones[i].FAC.rows();
    dim = ListCones[i].FAC.cols();
    Face f_fac(len);
    FaceDesc fd = {i, f_fac};
    ListDomain.push_back(fd);
  }
  std::vector<std::vector<FaceDesc>> ListListDomain;
  ListListDomain.push_back(ListDomain);
  //
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
  using Tface =
      std::pair<std::vector<ent_face<Tint>>, std::vector<MyMatrix<Tint>>>;
#endif
  for (int i = 1; i < TheLev; i++) {
    std::cerr << "i=" << i << "\n";
    std::vector<std::pair<size_t, Tent<T, Tint, Tidx_value>>> NewListCand;
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
    std::vector<Tface> list_face;
    size_t n_equiv_found = 0;
    auto f_insert_face = [&](const FaceDesc &fd_A) -> void {
      const ConeDesc<T, Tint, Tgroup> &eC = ListCones[fd_A.iCone];
      Face f_ext = Compute_faceEXT_from_faceFAC(eC.extfac_incd, eC.FAC.rows(),
                                                eC.EXT.rows(), fd_A.f_fac);
      ent_face<Tint> ef_A{fd_A.iCone, f_ext, IdentityMat<Tint>(dim)};
      for (auto &eOrbit : list_face) {
        for (auto &ef_B : eOrbit.first) {
          std::optional<MyMatrix<Tint>> equiv_opt =
              test_equiv_ent_face(ListCones, ef_A, ef_B);
          if (equiv_opt) {
            n_equiv_found++;
            Tent<T, Tint, Tidx_value> ent_A =
                f_ent<T, Tint, Tgroup, Tidx_value>(ListCones, G, fd_A, os);
            const ConeDesc<T, Tint, Tgroup> &eC_B = ListCones[ef_B.iCone];
            Face f_fac = Compute_faceFAC_from_faceEXT(
                eC_B.extfac_incd, eC_B.FAC.rows(), eC_B.EXT.rows(), ef_B.f_ext);
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
      list_face.push_back(
          get_spanning_list_ent_face(ListCones, ll_sing_adj, ef_A));
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
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
      f_insert_face(eEnt.fd);
#endif
      std::pair<size_t, Tent<T, Tint, Tidx_value>> e_pair{e_inv,
                                                          std::move(eEnt)};
      NewListCand.emplace_back(std::move(e_pair));
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
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
          ListCones[iCone].FAC, ListCones[iCone].EXT, os);
      for (auto &eFace_fac : ListFace) {
        FaceDesc fdn{iCone, eFace_fac};
        Tent<T, Tint, Tidx_value> fEnt =
            f_ent<T, Tint, Tgroup, Tidx_value>(ListCones, G, fdn, os);
        f_insert(std::move(fEnt));
      }
    }
    std::vector<FaceDesc> NewListDomain;
    for (auto &eEnt : NewListCand)
      NewListDomain.push_back(eEnt.second.fd);
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
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
    MyMatrix<Tint> const &G,
    const std::vector<std::vector<sing_adj<Tint>>> &ll_sing_adj, int TheLev,
    std::ostream &os) {
  std::vector<std::vector<FaceDesc>> ListListDomain;
  size_t n_col = G.rows();
  using Tface =
      std::pair<std::vector<ent_face<Tint>>, std::vector<MyMatrix<Tint>>>;
  std::vector<std::vector<Tface>> list_list_face;
  for (int iLev = 1; iLev <= TheLev; iLev++) {
    std::cerr << "iLev=" << iLev << "\n";
    std::vector<Tface> list_face;
    auto f_insert = [&](const ent_face<Tint> &ef_A) -> void {
      for (auto &eOrbit : list_face) {
        for (auto &ef_B : eOrbit.first) {
          std::optional<MyMatrix<Tint>> equiv_opt =
              test_equiv_ent_face(ListCones, ef_A, ef_B);
          if (equiv_opt)
            return;
        }
      }
      list_face.push_back(
          get_spanning_list_ent_face(ListCones, ll_sing_adj, ef_A));
    };
    if (iLev == 1) {
      for (size_t iCone = 0; iCone < ListCones.size(); iCone++) {
        vectface vf = DecomposeOrbitPoint_Full(ListCones[iCone].GRP_ext);
        size_t nExt = ListCones[iCone].EXT.rows();
        for (auto &eOrb : vf) {
          Face f_ext(nExt);
          boost::dynamic_bitset<>::size_type eVal = eOrb.find_first();
          size_t iExt = eVal;
          f_ext[iExt] = 1;
          ent_face<Tint> e_ent{iCone, f_ext, IdentityMat<Tint>(n_col)};
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
                                    eC.facext_incd, eC.EXT, eC.FAC, os);
          for (auto &f_ext_new : vf) {
            ent_face<Tint> e_ent{eRepr.iCone, f_ext_new,
                                 IdentityMat<Tint>(n_col)};
            f_insert(e_ent);
          }
        }
      }
    }
    std::vector<FaceDesc> ListDomain;
    for (auto &eOrbit : list_face) {
      ent_face<Tint> e_ent = eOrbit.first[0];
      const ConeDesc<T, Tint, Tgroup> &eC = ListCones[e_ent.iCone];
      Face f_fac = Compute_faceFAC_from_faceEXT(eC.extfac_incd, eC.FAC.rows(),
                                                eC.EXT.rows(), e_ent.f_ext);
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
  BlockPROC.setListIntValues(ListIntValues1_doc);
  BlockPROC.setListBoolValues(ListBoolValues1_doc);
  BlockPROC.setListStringValues(ListStringValues1_doc);
  ListBlock["PROC"] = BlockPROC;
  // Merging all data
  return {std::move(ListBlock), "undefined"};
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
    for (int i = 0; i < dim; i++)
      EXTtot(pos, i) = eEXT(i);
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
      T sum = 0;
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
  std::cerr << "DECOMP: SetFAC n_size1=" << n_size1 << " |SetFAC|=" << SetFAC.size() << "\n";
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

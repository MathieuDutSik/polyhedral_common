// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LORENTZIAN_EDGEWALK_H_
#define SRC_LORENTZIAN_EDGEWALK_H_

#include "Heuristic_ThompsonSampling.h"
#include "LatticeDefinitions.h"
#include "MatrixCanonicalForm.h"
#include "Namelist.h"
#include "POLY_RecursiveDualDesc.h"
#include "Temp_PolytopeEquiStab.h"
#include "Temp_Positivity.h"
#include "Timings.h"
#include "coxeter_dynkin.h"
#include "fund_domain_vertices.h"
#include "lorentzian_linalg.h"
#include "two_dim_lorentzian.h"
#include "vinberg.h"
#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#define ALLOW_VINBERG_ALGORITHM_FOR_INITIAL_VERTEX

template <typename F>
void print_stderr_stdout_file(std::string const &FileOut, F f) {
  if (FileOut == "stderr")
    return f(std::cerr);
  if (FileOut == "stdout")
    return f(std::cout);
  std::ofstream os(FileOut);
  return f(os);
}

FullNamelist NAMELIST_GetStandard_EDGEWALK() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListStringValues1["FileLorMat"] = "the lorentzian matrix used";
  ListStringValues1["OptionInitialVertex"] =
      "vinberg or FileVertex or FileVertexRoots and if FileVertex or "
      "FileVertexRoots selected use FileVertDomain as initial vertex";
  ListStringValues1["FileInitialVertex"] =
      "unset put the name of the file used for the initial vertex";
  ListStringValues1["OptionNorms"] =
      "possible option K3 (then just 2) or all where all norms are considered";
  ListStringValues1["DualDescProg"] = "lrs_iterate";
  ListStringValues1["OutFormat"] = "GAP for gap use or TXT for text output";
  ListStringValues1["FileOut"] =
      "stdout, or stderr or the filename of the file you want to write to";
  // Sometimes we can terminate by proving that it is not reflective
  ListBoolValues1["EarlyTerminationIfNotReflective"] = false;
  // Normally, we want to ApplyReduction, this is for debug only
  ListBoolValues1["ApplyReduction"] = true;
  ListBoolValues1["ComputeAllSimpleRoots"] = true;
  ListStringValues1["FileHeuristicIdealStabEquiv"] = "unset.heu";
  ListStringValues1["FileHeuristicTryTerminateDualDescription"] = "unset.heu";
  SingleBlock BlockPROC;
  BlockPROC.ListStringValues = ListStringValues1;
  BlockPROC.ListBoolValues = ListBoolValues1;
  ListBlock["PROC"] = BlockPROC;
  // Merging all data
  return {ListBlock, "undefined"};
}

FullNamelist NAMELIST_GetStandard_EDGEWALK_Isomorphism() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListStringValues1["FileLorMat1"] = "the lorentzian matrix used";
  ListStringValues1["FileLorMat2"] = "the lorentzian matrix used";
  ListStringValues1["OptionNorms"] =
      "possible option K3 (then just 2) or all where all norms are considered";
  ListStringValues1["DualDescProg"] = "lrs_iterate";
  ListStringValues1["OutFormat"] = "GAP for gap use or TXT for text output";
  ListStringValues1["FileOut"] =
      "stdout, or stderr or the filename of the file you want to write to";
  // Normally, we want to ApplyReduction, this is for debug only
  ListBoolValues1["ApplyReduction"] = true;
  SingleBlock BlockPROC;
  BlockPROC.ListStringValues = ListStringValues1;
  BlockPROC.ListBoolValues = ListBoolValues1;
  ListBlock["PROC"] = BlockPROC;
  // Merging all data
  return {ListBlock, "undefined"};
}

template <typename T>
MyMatrix<T> ComputeLattice_LN(MyMatrix<T> const &G, T const &N) {
  int n = G.rows();
  MyMatrix<T> M1 = IdentityMat<T>(n);
  MyMatrix<T> M2 = (N / 2) * Inverse(G);
  return IntersectionLattice(M1, M2);
}

template <typename T> struct SublattInfos {
  MyMatrix<T> G;
  std::vector<T> l_norms;
  std::unordered_map<T, MyMatrix<T>> map_norm_latt;
};

template <typename T, typename Tint>
SublattInfos<T> ComputeSublatticeInfos(MyMatrix<T> const &G,
                                       std::vector<T> const &l_norms) {
  std::unordered_map<T, MyMatrix<T>> map_norm_latt;
  for (auto &e_norm : l_norms) {
    MyMatrix<T> Latt_pre = ComputeLattice_LN(G, e_norm);
    MyMatrix<T> Latt = LLLbasisReduction<T, Tint>(Latt_pre).LattRed;
    map_norm_latt[e_norm] = Latt;
  }
  return {G, l_norms, std::move(map_norm_latt)};
}


// Computation related to the enumeration algorithm.
// Since it has several square roots, we need to keep track.
//
// We have
// sign: 0 for 0, 1 for positive, -1 for negative
// quant1: this is (k.alpha_{N,\Delta'})^2 / R_{N,\Delta'}
// quant2: this is (k.alpha_{N,\Delta'})^2 / N
template <typename T, typename Tint> struct RootCandidate {
  int sign;
  T quant1;
  T quant2;
  T e_norm;
  MyVector<Tint> alpha;
  FundDomainVertex<T, Tint> fund_v;
};

template <typename T> int get_sign_sing(T const &val) {
  if (val > 0)
    return 1;
  if (val < 0)
    return -1;
  return 0;
}

template <typename T> int get_sign_pair_t(T const &p1, T const &p2) {
  if (p1 < p2)
    return 1;
  if (p1 > p2)
    return -1;
  return 0;
}

// return 0 is p1 == p2 :  1 if p1 < p2 : -1 if p1 > p2
template <typename T>
int get_sign_pair_stdpair(std::pair<int, T> const &p1,
                          std::pair<int, T> const &p2) {
  if (p1.first != p2.first)
    return get_sign_pair_t(p1.first, p2.first);
  if (p1.first == 1)
    return get_sign_pair_t(p1.second, p2.second);
  if (p1.first == -1)
    return -get_sign_pair_t(p1.second, p2.second);
  return 0;
}

template <typename T, typename Tint>
RootCandidate<T, Tint>
gen_possible_extension(MyMatrix<T> const &G, MyVector<T> const &k,
                       MyVector<Tint> const &alpha, T const &res_norm,
                       T const &e_norm,
                       FundDomainVertex<T, Tint> const &fund_v) {
  MyVector<T> alpha_T = UniversalVectorConversion<T, Tint>(alpha);
  T scal = -alpha_T.dot(G * k);
  T quant1 = (scal * scal) / res_norm;
  T quant2 = (scal * scal) / e_norm;
  return {get_sign_sing(scal), quant1, quant2, e_norm, alpha, fund_v};
}

// sign as before + means preferable according to page 27.
// return 1 if poss1 is preferable to poss2
template <typename T, typename Tint>
int get_sign_cand(RootCandidate<T, Tint> const &poss1,
                  RootCandidate<T, Tint> const &poss2) {
  int sign1 = get_sign_pair_stdpair<T>({poss1.sign, poss1.quant1},
                                       {poss2.sign, poss2.quant1});
  if (sign1 != 0) {
    // because -k.alpha1 / sqrt(R1)    <     -k.alpha2 / sqrt(R2)
    // correspond to 1 in the above.
    return sign1;
  }
  int sign2 = get_sign_pair_stdpair<T>({poss1.sign, poss1.quant2},
                                       {poss2.sign, poss2.quant2});
  if (sign2 != 0) {
    // because -k.alpha1 / sqrt(N1)    <     -k.alpha2 / sqrt(N2)
    // correspond to 1 in the above.
    return sign2;
  }
  // because N1 < N2 corresponds to 1
  return get_sign_pair_t(poss1.e_norm, poss2.e_norm);
}

template <typename T, typename Tint>
RootCandidate<T, Tint>
get_best_candidate(std::vector<RootCandidate<T, Tint>> const &l_cand) {
  if (l_cand.size() == 0) {
    std::cerr << "We have zero candidates. Abort\n";
    throw TerminalException{1};
  }
  RootCandidate<T, Tint> best_cand = l_cand[0];
  for (size_t i = 1; i < l_cand.size(); i++)
    if (get_sign_cand(l_cand[i], best_cand) == 1)
      best_cand = l_cand[i];
  return best_cand;
}

template <typename T, typename Tint> struct CuspidalRequest {
  std::vector<MyVector<Tint>> l_ui;
  MyVector<T> k;
  MyVector<T> kP;
};

template <typename T, typename Tint> struct CuspidalRequest_FullInfo {
  CuspidalRequest<T, Tint> eRequest;
  pair_char<T> e_pair;
  size_t hash;
};

template <typename T, typename Tint> struct CuspidalBank {
  std::vector<CuspidalRequest_FullInfo<T, Tint>> l_request;
  std::vector<std::vector<MyVector<Tint>>> l_answer;
};

template <typename T, typename Tint>
CuspidalRequest_FullInfo<T, Tint>
gen_cuspidal_request_full_info(MyMatrix<T> const &G,
                               CuspidalRequest<T, Tint> const &eReq) {
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "gen_cuspidal_request_full_info, step 1\n";
#endif
  std::unordered_map<MyVector<Tint>, uint8_t> map_v;
  std::vector<MyVector<Tint>> l_vect;
  std::vector<T> Vdiag;
  for (auto &eV : eReq.l_ui) {
    l_vect.push_back(eV);
    Vdiag.push_back(1);
  }
  MyVector<Tint> k_tint =
      UniversalVectorConversion<Tint, T>(RemoveFractionVector(eReq.k));
  l_vect.push_back(k_tint);
  Vdiag.push_back(2);
  MyVector<Tint> kp_tint =
      UniversalVectorConversion<Tint, T>(RemoveFractionVector(eReq.kP));
  l_vect.push_back(kp_tint);
  Vdiag.push_back(3);
  size_t n_row = Vdiag.size();
  //
  using Tidx = uint32_t;
  using Tidx_value = uint16_t;
  using Tgr = GraphListAdj;
  std::vector<MyMatrix<T>> ListMat{G};
  //
  MyMatrix<T> MatV =
      UniversalMatrixConversion<T, Tint>(MatrixFromVectorFamily(l_vect));
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "gen_cuspidal_request_full_info, step 2\n";
#endif
  WeightMatrix<true, std::vector<T>, Tidx_value> WMat =
      GetWeightMatrix_ListMat_Vdiag<T, Tidx, Tidx_value>(MatV, ListMat, Vdiag);
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "gen_cuspidal_request_full_info, step 3\n";
#endif
  WMat.ReorderingSetWeight();
  std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> epair =
      GetGroupCanonicalizationVector_Kernel<std::vector<T>, Tgr, Tidx,
                                            Tidx_value>(WMat);
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "gen_cuspidal_request_full_info, step 4\n";
#endif
  const std::vector<Tidx> &ListIdx = epair.first;
  WMat.RowColumnReordering(ListIdx);
  //
  std::vector<MyVector<Tint>> l_vect_reord(n_row);
  for (size_t i = 0; i < n_row; i++) {
    size_t j = ListIdx[i];
    l_vect_reord[i] = l_vect[j];
  }
  MyMatrix<T> MatV_reord =
      UniversalMatrixConversion<T, Tint>(MatrixFromVectorFamily(l_vect_reord));
  //
  size_t seed = 1440;
  size_t hash = ComputeHashWeightMatrix_raw(WMat, seed);
  pair_char<T> e_pair{std::move(MatV_reord), std::move(WMat)};
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "gen_cuspidal_request_full_info, step 5\n";
#endif
  return {eReq, std::move(e_pair), hash};
}

/*
  We use solution of Problem 7.1 in the edgewalk text.
  It is indeed true that we are searching roots that contain k.
  The dimension of that space is n.
  But since k is isotropic that space is actually the span of k and the l_ui.
  --
  Now if we write a vector as u + ck then we get N(u + ck) = N(u)
 */
template <typename T, typename Tint>
std::vector<MyVector<Tint>>
DetermineRootsCuspidalCase(SublattInfos<T> const &si,
                           CuspidalRequest<T, Tint> const &eReq) {
  MyMatrix<T> const &G = si.G;
  std::vector<T> const &l_norms = si.l_norms;
#ifdef TIMINGS
  MicrosecondTime time;
#endif
  std::vector<MyVector<Tint>> const &l_ui = eReq.l_ui;
  MyVector<T> const &k = eReq.k;
  MyVector<T> const &kP = eReq.kP;
  // 0 for 0, 1 for positive, -1 for negative
  // this is (kP.v_{N,\Delta'})^2 / N
  struct RootCandidateCuspidal {
    int sign;
    T quant;
    T e_norm;
    MyVector<Tint> v;
  };
  auto gen_possible_cuspidalextension =
      [&](MyVector<T> const &v_T, T const &e_norm) -> RootCandidateCuspidal {
    MyVector<Tint> v = UniversalVectorConversion<Tint, T>(v_T);
    T scal = -kP.dot(G * v_T);
    T quant = (scal * scal) / e_norm;
    return {get_sign_sing(scal), quant, e_norm, v};
  };
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "DetermineRootsCuspidalCase, beginning\n";
#endif
  bool only_spherical = false;
  std::vector<Possible_Extension<T>> l_extension =
      ComputePossibleExtensions(G, l_ui, l_norms, only_spherical);
#ifdef TIMINGS
  std::cerr << "Timing |ComputePossibleExtensions|=" << time << "\n";
#endif
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "DetermineRootsCuspidalCase : |l_extension|="
            << l_extension.size() << "\n";
#endif
  std::vector<RootCandidateCuspidal> l_candidates;
  for (auto &e_extension : l_extension) {
    if (e_extension.res_norm == 0) {
      MyMatrix<T> const &Latt = si.map_norm_latt.at(e_extension.e_norm);
      std::optional<MyVector<T>> opt_v =
          ResolveLattEquation(Latt, e_extension.u_component, k);
      if (opt_v) {
        const MyVector<T> &v_T = *opt_v;
        RootCandidateCuspidal e_cand =
            gen_possible_cuspidalextension(v_T, e_extension.e_norm);
        l_candidates.push_back(e_cand);
      }
    }
  }
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "DetermineRootsCuspidalCase : |l_candidates|="
            << l_candidates.size() << "\n";
#endif
#ifdef TIMINGS
  std::cerr << "Timing |l_candidates|=" << time << "\n";
#endif
  /* std::sort is sorting from the highest to the smallest
   */
  std::sort(l_candidates.begin(), l_candidates.end(),
            [&](RootCandidateCuspidal const &x,
                RootCandidateCuspidal const &y) -> bool {
              // We want x > y if
              // -k.alpha(x) / sqrt(Nx) < -k.alpha(y) / sqrt(Ny)
              // or if equality if Nx < Ny
              int sign = get_sign_pair_stdpair<T>({x.sign, x.quant},
                                                  {y.sign, y.quant});
              if (sign != 0) {
                // because -k.alpha1 / sqrt(R1)    < -k.alpha2 / sqrt(R2)
                // correspond to 1 in the above.
                return sign > 0;
              }
              return x.e_norm < y.e_norm;
            });
#ifdef TIMINGS
  std::cerr << "Timing |sort|=" << time << "\n";
#endif
#ifdef DEBUG_EDGEWALK_GENERIC
  for (auto &x : l_candidates) {
    std::cerr << "x : sign=" << x.sign << " quant=" << x.quant
              << " norm=" << x.e_norm << " v=" << StringVectorGAP(x.v) << "\n";
  }
#endif
  std::vector<MyVector<Tint>> l_ui_ret = l_ui;
  auto is_approved = [&](MyVector<Tint> const &cand) -> bool {
    MyVector<T> G_cand_T = G * UniversalVectorConversion<T, Tint>(cand);
    for (auto &v : l_ui_ret) {
      MyVector<T> v_T = UniversalVectorConversion<T, Tint>(v);
      T scal = v_T.dot(G_cand_T);
      if (scal > 0)
        return false;
    }
    return true;
  };
  for (auto &eV : l_candidates) {
    MyVector<Tint> eV_i = eV.v;
    if (is_approved(eV_i)) {
#ifdef DEBUG_EDGEWALK_GENERIC
      std::cerr << "Inserting eV_i=";
      WriteVector(std::cerr, eV_i);
#endif
      l_ui_ret.push_back(eV_i);
    }
  }
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "DetermineRootsCuspidalCase, exiting |l_ui_ret|="
            << l_ui_ret.size() << "\n";
#endif
#ifdef TIMINGS
  std::cerr << "Timing |l_ui_ret|=" << time << "\n";
#endif
  return l_ui_ret;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyVector<Tint>>
DetermineRootsCuspidalCase_Memoized(CuspidalBank<T, Tint> &cusp_bank,
                                    SublattInfos<T> const &si,
                                    CuspidalRequest<T, Tint> const &eReq) {
#ifdef TIMINGS
  MicrosecondTime time;
#endif
  MyMatrix<T> const &G = si.G;
  CuspidalRequest_FullInfo<T, Tint> eReq_full =
      gen_cuspidal_request_full_info(G, eReq);
#ifdef TIMINGS
  std::cerr << "Timing |gen_cuspidal_request_full_info|=" << time << "\n";
#endif
  size_t len = cusp_bank.l_request.size();
  for (size_t i = 0; i < len; i++) {
    const CuspidalRequest_FullInfo<T, Tint> &fReq_full = cusp_bank.l_request[i];
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "i=" << i << "/" << len << " hash: eReq=" << eReq_full.hash
              << " fReq=" << fReq_full.hash << "\n";
#endif
    if (fReq_full.hash == eReq_full.hash) {
      std::optional<MyMatrix<T>> equiv_opt =
          LinPolytopeIntegralWMat_Isomorphism<T, Tgroup, std::vector<T>,
                                              uint16_t>(fReq_full.e_pair,
                                                        eReq_full.e_pair);
      if (equiv_opt) {
        MyMatrix<Tint> eEquiv = UniversalMatrixConversion<Tint, T>(*equiv_opt);
        std::vector<MyVector<Tint>> l_ui_ret;
        for (auto &eV : cusp_bank.l_answer[i]) {
          MyVector<Tint> Vret = eEquiv.transpose() * eV;
          l_ui_ret.push_back(Vret);
        }
#ifdef DEBUG_EDGEWALK_GENERIC
        std::cerr
            << "DetermineRootsCuspidalCase_Memoized, find some isomorphism\n";
#endif
#ifdef TIMINGS
        std::cerr << "Timing |query(succ)|=" << time << "\n";
#endif
        return l_ui_ret;
      }
    }
  }
#ifdef TIMINGS
  std::cerr << "Timing |query(fail)|=" << time << "\n";
#endif
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "DetermineRootsCuspidalCase_Memoized, failed to find some "
               "isomorphism\n";
#endif
  std::vector<MyVector<Tint>> l_ui_ret = DetermineRootsCuspidalCase(si, eReq);
  cusp_bank.l_request.emplace_back(std::move(eReq_full));
  cusp_bank.l_answer.push_back(l_ui_ret);
#ifdef TIMINGS
  std::cerr << "Timing |l_ui_ret|=" << time << "\n";
#endif
  return l_ui_ret;
}

template <typename Tint> struct AdjacencyDirection {
  std::vector<MyVector<Tint>> l_ui;
  MyVector<Tint> v_disc;
};

template <typename Tint>
void PrintAdjacencyDirection(std::ostream &os,
                             AdjacencyDirection<Tint> const &ad) {
  os << "l_ui =";
  for (auto &root : ad.l_ui)
    os << " " << StringVectorGAP(root);
  os << "\n";
  os << "v_disc =" << StringVectorGAP(ad.v_disc) << "\n";
}

template <typename Tint>
std::string StringAdjacencyDirectionGAP(AdjacencyDirection<Tint> const &ad) {
  std::string ret = "rec(v_disc:=" + StringVectorGAP(ad.v_disc) + ", l_ui:=" +
                    StringMatrixGAP(MatrixFromVectorFamily(ad.l_ui)) + ")";
  return ret;
}

template <typename Tint>
AdjacencyDirection<Tint> GetAdjacencyDirection(MyMatrix<Tint> const &MatRoot,
                                               Face const &f) {
  size_t n_root = MatRoot.rows();
  size_t i_disc = std::numeric_limits<size_t>::max();
  std::vector<MyVector<Tint>> l_ui;
  for (size_t i_root = 0; i_root < n_root; i_root++) {
    if (f[i_root] == 1) {
      MyVector<Tint> root = GetMatrixRow(MatRoot, i_root);
      l_ui.push_back(root);
    } else {
      i_disc = i_root;
    }
  }
  MyVector<Tint> v_disc = GetMatrixRow(MatRoot, i_disc);
  return {l_ui, v_disc};
}

/*
  We take the notations as in EDGEWALK paper.
  ---The dimension is n+1
  ---We have an edge e between two rays k and k'.
  We have u1, .... u(n-1) roots that are orthogonal and U the real span of those
  vectors.
  ---P is the real plane orthogonal to U.
  ---pi_U and pi_P are the corresponding projectors
  ---How is (1/2) P defined and correspond to k (typo correction)
 */
template <typename T, typename Tint, typename Tgroup>
FundDomainVertex<T, Tint>
EdgewalkProcedure(CuspidalBank<T, Tint> &cusp_bank, SublattInfos<T> const &si,
                  MyVector<T> const &k, AdjacencyDirection<Tint> const ad) {
  MyMatrix<T> const &G = si.G;
  std::vector<T> const &l_norms = si.l_norms;
#ifdef TIMINGS
  Microsecond time;
#endif
  const std::vector<MyVector<Tint>> &l_ui = ad.l_ui;
  const MyVector<Tint> &v_disc = ad.v_disc;
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "-------------------------------------- EDGEWALK PROCEDURE "
               "---------------------------------------------\n";
  std::cerr << "k=" << StringVectorGAP(k) << "\n";
  std::cerr << "l_norms =";
  for (auto &eN : l_norms)
    std::cerr << " " << eN;
  std::cerr << "\n";
  std::cerr << "l_ui =";
  for (auto &eV : l_ui)
    std::cerr << " " << StringVectorGAP(eV);
  std::cerr << "\n";
  std::cerr << "v_disc=" << StringVectorGAP(v_disc) << "\n";
  std::cerr << "    Real work starts now\n";
#endif
  //
  // Initial computation of linear algebra nature:
  // Find a basis (k,r0) of the plane P
  //
  MyVector<T> v_disc_t = UniversalVectorConversion<T, Tint>(v_disc);
  int n = G.rows();
  size_t n_root = l_ui.size();
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "n_root=" << n_root << "\n";
#endif
  MyMatrix<T> EquaRvect(n_root + 1, n);
  for (size_t i_root = 0; i_root < n_root; i_root++) {
    MyVector<T> eV = UniversalVectorConversion<T, Tint>(l_ui[i_root]);
    MyVector<T> eP = G * eV;
    T eScal = k.dot(eP);
    if (eScal != 0) {
      std::cerr << "The scalar product should be 0\n";
      throw TerminalException{1};
    }
    AssignMatrixRow(EquaRvect, i_root, eP);
  }
  MyMatrix<T> Pplane = Get_Pplane(G, l_ui);
  MyVector<T> eP = G * k;
  T norm = k.dot(eP);
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "k=" << StringVectorGAP(k) << " norm=" << norm << "\n";
#endif
  AssignMatrixRow(EquaRvect, n_root, eP);
  MyMatrix<T> NSP = NullspaceTrMat(EquaRvect);
  if (NSP.rows() != 1) {
    std::cerr << "|NSP|=" << NSP.rows() << "/" << NSP.cols() << "\n";
    std::cerr << "The dimension should be exactly 1\n";
    throw TerminalException{1};
  }
  /*
    The vector r0 is orthogonal to k and is well defined up to a sign.
    The half plane (1/2)P is shown on Figure 8.1 as being orthogonal to
    k. Note that the mention on page 26 "(1/)2P is the open half plane
    corresponding to e" is not correct. It should be corresponding to k. But how
    to build it? How to select the sign?
   */
  MyVector<T> r0 = GetMatrixRow(NSP, 0);
  T scal_r0 = r0.dot(G * v_disc_t);
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "scal_r0=" << scal_r0 << "\n";
#endif
  if (scal_r0 > 0)
    r0 = -r0;
  // We follow here the convention on oriented basis of Section 8:
  // "First member lies in the interior of (1/2)P and whose second member is k"
  MyMatrix<T> OrientedBasis(2, n);
  if (norm < 0) { // the point is inner, the oriented basis is clear.
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "Builging OrientedBasis, ordinary case\n";
#endif
    AssignMatrixRow(OrientedBasis, 0, r0);
    AssignMatrixRow(OrientedBasis, 1, k);
  } else { // Now the point is ideal
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "Builging OrientedBasis, ideal case\n";
#endif
    /*
      The roots to be found are of positive norms.
      In general, we cannot think of which zone of positive norms to select from
      since those are connected. However, in the specific case of dimension 2,
      yes the vectors of positive norms are in two connected components. How to
      select which connected component? The scalar product with k allows us to
      select which one
     */
    auto get_positive_norm_vector = [&]() -> MyVector<T> {
      T two = 2;
      T alpha = 1;
      while (true) {
        for (int i_bas = 0; i_bas < 2; i_bas++) {
          MyVector<T> v_bas = GetMatrixRow(Pplane, i_bas);
          for (int u = -1; u < 2; u += 2) {
            MyVector<T> v_pos_cand = k + u * alpha * v_bas;
            T norm = v_pos_cand.dot(G * v_pos_cand);
#ifdef DEBUG_EDGEWALK_GENERIC
            std::cerr << "u=" << u << " alpha=" << alpha
                      << " v_pos_cand=" << v_pos_cand << " norm=" << norm
                      << "\n";
#endif
            if (norm > 0)
              return v_pos_cand;
          }
        }
        alpha /= two;
      }
    };
    MyVector<T> v_pos = get_positive_norm_vector();
    T scal = v_pos.dot(G * k);
    if (scal >
        0) { // The convention is that negative scalar product is for facets.
      v_pos = -v_pos;
    }
    AssignMatrixRow(OrientedBasis, 0, v_pos);
    AssignMatrixRow(OrientedBasis, 1, k);
    r0 = -k; // Follows Right part of Figure 8.1
  }
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "r0=" << StringVectorGAP(r0) << "\n";
#endif
#ifdef TIMINGS
  std::cerr << "Timing |paperwork|=" << time << "\n";
#endif
  //
  // Computing the extension and the maximum norms from that.
  //
  std::vector<RootCandidate<T, Tint>> l_candidates;
  bool only_spherical = true;
  std::vector<Possible_Extension<T>> l_extension =
      ComputePossibleExtensions(G, l_ui, l_norms, only_spherical);
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "EdgewalkProcedure : |l_extension|=" << l_extension.size()
            << "\n";
#endif
#ifdef TIMINGS
  std::cerr << "Timing |l_extension|=" << time << "\n";
#endif

  //
  // Determine if the plane P is isotropic and if not compute the set of test
  // vectors
  //
  MyMatrix<T> G_Pplane = Pplane * G * Pplane.transpose();
  struct SingCompAnisotropic {
    MyMatrix<T> Latt;
    MyVector<Tint> r0_work;
    MyMatrix<T> Basis_ProjP_LN;
    MyMatrix<T> Basis_P_inter_LN;
    MyMatrix<T> Gwork;
    std::vector<MyVector<Tint>> l_vect;
  };
  struct SingCompIsotropic {
    MyMatrix<T> Latt;
    MyMatrix<T> Basis_ProjP_LN;
    MyMatrix<T> GP_LN;
    MyMatrix<T> Factor_GP_LN;
    MyVector<Tint> r0_work;
    std::map<T, std::optional<std::vector<MyVector<Tint>>>> map_res_norm;
  };
  auto get_basis_projp_ln = [&](MyMatrix<T> const &Latt) -> MyMatrix<T> {
    LatticeProjectionFramework<T> ProjFram(G, Pplane, Latt);
    MyMatrix<T> BasisProj = ProjFram.BasisProj;
    if (BasisProj.rows() != 2) {
      std::cerr << "The BasisProj should be of rank 2\n";
      throw TerminalException{1};
    }
    MyMatrix<T> Expr =
        ExpressVectorsInIndependentFamilt(BasisProj, OrientedBasis);
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "Det(Expr)=" << DeterminantMat(Expr) << "\n";
#endif
    if (DeterminantMat(Expr) < 0) { // Change to get positive determinant
      for (int i = 0; i < n; i++)
        BasisProj(0, i) = -BasisProj(0, i);
    }
    return BasisProj;
  };
  auto get_r0work = [&](MyMatrix<T> const &Basis_ProjP_LN,
                        MyVector<T> const &r0) -> MyVector<Tint> {
    std::optional<MyVector<T>> opt_r0 = SolutionMat(Basis_ProjP_LN, r0);
    if (!opt_r0) {
      std::cerr << "Failed to resolve the SolutionMat problem\n";
      throw TerminalException{1};
    }
    MyVector<T> r0_NSP = *opt_r0;
    MyVector<Tint> r0_work =
        UniversalVectorConversion<Tint, T>(RemoveFractionVector(r0_NSP));
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "r0_work=" << StringVectorGAP(r0_work) << "\n";
#endif
    return r0_work;
  };
  auto get_sing_comp_anisotropic = [&](T const &e_norm) -> SingCompAnisotropic {
#ifdef TIMINGS
    MicrosecondTime timeA;
#endif
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << " -------------- get_sing_comp_anisotropic, e_norm=" << e_norm
              << " ------------------------\n";
#endif
    MyMatrix<T> const &Latt = si.map_norm_latt.at(e_norm);
    MyMatrix<T> Basis_ProjP_LN = get_basis_projp_ln(Latt);
    MyMatrix<T> Basis_P_inter_LN =
        IntersectionLattice_VectorSpace(Latt, Pplane);
    MyMatrix<T> Gwork = Basis_ProjP_LN * G * Basis_ProjP_LN.transpose();
    // The residual norm is res_norm = N - u_{N,Delta}^2
    // u_{N,Delta} belongs to a positive definite lattice.
    // Therefore res_norm <= N
    // Thus we compute all the vectors up to norm res_norm because
    // res_norm is always realizable with the vector 2 2 2 ..... 2
    T res_norm = e_norm;
    MyVector<Tint> r0_work = get_r0work(Basis_ProjP_LN, r0);
    T r0_norm = eval_quad(Gwork, r0_work);
    MyVector<Tint> l_A = GetTwoComplement(r0_work);
    MyVector<Tint> l_B = Canonical(Gwork, r0_norm, r0_work, l_A);
    std::optional<std::pair<MyMatrix<Tint>, std::vector<MyVector<Tint>>>> opt =
        Anisotropic<T, Tint>(Gwork, res_norm, r0_work, l_B);
    if (!opt) { // No solution, this definitely can happen
      return {Latt, r0_work, Basis_ProjP_LN, Basis_P_inter_LN, Gwork, {}};
    }
    const std::vector<MyVector<Tint>> &l_vect1 = opt->second;
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "|l_vect1|=" << l_vect1.size() << "\n";
#endif
    const MyMatrix<Tint> &TransformSma = opt->first;
    MyMatrix<T> TransformSma_T =
        UniversalMatrixConversion<T, Tint>(TransformSma);
    /*
      Transformation rule
     */
    MyMatrix<T> TransformBig_red = IdentityMat<T>(n);
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
        TransformBig_red(i, j) =
            UniversalScalarConversion<T, Tint>(TransformSma(i, j));
    MyMatrix<T> l_ui_MatT =
        UniversalMatrixConversion<T, Tint>(MatrixFromVectorFamily(l_ui));
    MyMatrix<T> ReductionBasis = Concatenate(Basis_ProjP_LN, l_ui_MatT);
    MyMatrix<T> TransformBig =
        Inverse(ReductionBasis) * TransformBig_red * ReductionBasis;
    MyMatrix<T> Expr_t =
        ExpressVectorsInIndependentFamilt(Basis_P_inter_LN, Basis_ProjP_LN);
    if (!IsIntegralMatrix(Expr_t)) {
      std::cerr << "The matrix should be integral\n";
      throw TerminalException{1};
    }
    size_t order =
        GetMatrixExponentSublattice_TrivClass(TransformSma_T, Expr_t);
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "order=" << order << "\n";
#endif
    std::vector<MyMatrix<Tint>> l_vect2;
    for (auto &e_vect1 : l_vect1) {
      T norm1 = eval_quad(Gwork, e_vect1);
      size_t ord = 1;
      while (true) {
        T norm2 = ord * ord * norm1;
        if (norm2 > res_norm)
          break;
        MyVector<Tint> e_vect2 = ord * e_vect1;
        l_vect2.push_back(e_vect2);
        ord++;
      }
    }
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "|l_vect2|=" << l_vect2.size() << "\n";
#endif
    std::vector<MyVector<Tint>> l_vect3;
    MyMatrix<Tint> TheMat = IdentityMat<Tint>(2);
    for (size_t i = 0; i < order; i++) {
      for (auto &e_vect2 : l_vect2) {
        MyVector<Tint> e_vect3 = TheMat.transpose() * e_vect2;
        l_vect3.push_back(e_vect3);
      }
      TheMat = TheMat * TransformSma;
    }
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "|l_vect3|=" << l_vect3.size() << "\n";
#endif
#ifdef TIMINGS
    std::cerr << "Timing |get_sing_comp_anisotropic|=" << timeA << "\n";
#endif
    return {Latt, r0_work, Basis_ProjP_LN, Basis_P_inter_LN, Gwork, l_vect3};
  };
  auto get_sing_comp_isotropic = [&](T const &e_norm) -> SingCompIsotropic {
#ifdef TIMINGS
    MicrosecondTime timeB;
#endif
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "get_sing_comp_isotropic, e_norm=" << e_norm << "\n";
#endif
    MyMatrix<T> const &Latt = si.map_norm_latt.at(e_norm);
    MyMatrix<T> Basis_ProjP_LN = get_basis_projp_ln(Latt);
    MyMatrix<T> GP_LN = Basis_ProjP_LN * G * Basis_ProjP_LN.transpose();
    MyVector<Tint> r0_work = get_r0work(Basis_ProjP_LN, r0);
    std::optional<MyMatrix<T>> opt_factor = GetIsotropicFactorization(GP_LN);
    if (!opt_factor) {
      std::cerr << "Computation of isotropy factorization failed\n";
      throw TerminalException{1};
    }
    MyMatrix<T> Factor_GP_LN = *opt_factor;
#ifdef TIMINGS
    std::cerr << "Timing |get_sing_comp_isotropic|=" << timeB << "\n";
#endif
    return {Latt, Basis_ProjP_LN, GP_LN, Factor_GP_LN, r0_work, {}};
  };
  bool is_isotropic;
  std::map<T, SingCompAnisotropic> map_anisotropic;
  std::map<T, SingCompIsotropic> map_isotropic;
  {
    MyMatrix<T> Gtest = Pplane * G * Pplane.transpose();
    std::optional<MyMatrix<T>> opt = GetIsotropicFactorization(Gtest);
    if (opt) {
      is_isotropic = true;
#ifdef DEBUG_EDGEWALK_GENERIC
      std::cerr << "Case is_isotropic = true\n";
#endif
      for (auto &u_norm : l_norms)
        map_isotropic[u_norm] = get_sing_comp_isotropic(u_norm);
    } else {
      is_isotropic = false;
#ifdef DEBUG_EDGEWALK_GENERIC
      std::cerr << "Case is_isotropic = false\n";
#endif
      for (auto &u_norm : l_norms)
        map_anisotropic[u_norm] = get_sing_comp_anisotropic(u_norm);
    }
  }
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "Edgewalk Procedure, step 7\n";
#endif
  // Evaluation of fun
  auto get_next_anisotropic =
      [&](Possible_Extension<T> const &poss) -> std::optional<MyVector<Tint>> {
    T const &e_norm = poss.e_norm;
    SingCompAnisotropic const &e_comp = map_anisotropic[e_norm];
    for (auto &e_vect : e_comp.l_vect) {
      T val = eval_quad(e_comp.Gwork, e_vect);
      if (val == poss.res_norm) {
        MyVector<T> v_T =
            poss.u_component + e_comp.Basis_ProjP_LN.transpose() *
                                   UniversalVectorConversion<T, Tint>(e_vect);
        if (IsIntegerVector(v_T)) {
          std::optional<MyVector<T>> eSol = SolutionIntMat(e_comp.Latt, v_T);
          if (eSol) {
            MyVector<Tint> v_i = UniversalVectorConversion<Tint, T>(v_T);
#ifdef DEBUG_EDGEWALK_GENERIC
            std::cerr << "Returning v_i=" << StringVectorGAP(v_i) << "\n";
#endif
            return v_i;
          }
        }
      }
    }
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "No good vector found\n";
#endif
    return {};
  };
  auto get_successive_list_cand =
      [&](SingCompIsotropic &e_comp,
          T const &res_norm) -> std::vector<MyVector<Tint>> {
    std::optional<std::vector<MyVector<Tint>>> &opt =
        e_comp.map_res_norm[res_norm];
    if (opt)
      return *opt;
    std::vector<MyVector<Tint>> l_vect1 =
        EnumerateVectorFixedNorm_Factorization<T, Tint>(e_comp.Factor_GP_LN,
                                                        res_norm);
    std::vector<MyVector<Tint>> l_vect2;
    auto det = [&](MyVector<Tint> const &v1, MyVector<Tint> const &v2) -> Tint {
      return v1(0) * v2(1) - v1(1) * v2(0);
    };
    auto is_corr = [&](MyVector<Tint> const &x) -> bool {
      // In the isotropic case, there is no relevant test
      // See Right picture on Figure 8.1 of the edgewalk paper
      if (norm < 0) {
        T scal = eval_scal(e_comp.GP_LN, e_comp.r0_work, x);
        if (scal <= 0) {
          return false;
        }
      }
      Tint eDet = det(e_comp.r0_work, x);
      return eDet > 0;
    };
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "|l_vect1|=" << l_vect1.size() << "\n";
#endif
    for (auto &e_vect1 : l_vect1)
      if (is_corr(e_vect1))
        l_vect2.push_back(e_vect1);
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "We found |l_vect2|=" << l_vect2.size() << "\n";
#endif
    std::sort(l_vect2.begin(), l_vect2.end(),
              [&](MyVector<Tint> const &x, MyVector<Tint> const &y) -> bool {
                return det(x, y) > 0;
              });
    opt = l_vect2;
    return l_vect2;
  };
  auto get_next_isotropic =
      [&](Possible_Extension<T> const &poss) -> std::optional<MyVector<Tint>> {
    T const &e_norm = poss.e_norm;
    SingCompIsotropic &e_comp = map_isotropic[e_norm];
    T const &res_norm = poss.res_norm;
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "get_next_isotropic with e_norm=" << e_norm
              << " res_norm=" << res_norm << "\n";
#endif
    std::vector<MyVector<Tint>> l_vect =
        get_successive_list_cand(e_comp, res_norm);
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "|l_vect|=" << l_vect.size() << "\n";
#endif
    for (auto &e_vect : l_vect) {
#ifdef DEBUG_EDGEWALK_GENERIC
      std::cerr << "e_vect=" << StringVectorGAP(e_vect) << "\n";
#endif
      MyVector<T> v_T =
          poss.u_component + e_comp.Basis_ProjP_LN.transpose() *
                                 UniversalVectorConversion<T, Tint>(e_vect);
#ifdef DEBUG_EDGEWALK_GENERIC
      std::cerr << "After v_T computation v_T=" << StringVectorGAP(v_T) << "\n";
#endif
      if (IsIntegerVector(v_T)) {
        std::optional<MyVector<T>> eSol = SolutionIntMat(e_comp.Latt, v_T);
        if (eSol) {
          MyVector<Tint> v_i = UniversalVectorConversion<Tint, T>(v_T);
#ifdef DEBUG_EDGEWALK_GENERIC
          std::cerr << "get_next_isotropic. Returning v_i="
                    << StringVectorGAP(v_i) << "\n";
#endif
          return v_i;
        }
      }
    }
    return {};
  };
  auto get_next =
      [&](Possible_Extension<T> const &poss) -> std::optional<MyVector<Tint>> {
    if (is_isotropic) {
      return get_next_isotropic(poss);
    } else {
      return get_next_anisotropic(poss);
    }
  };
  //
  //
  //
  for (auto &e_extension : l_extension) {
    T e_norm = e_extension.e_norm;
    T res_norm = e_extension.res_norm;
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "------ u_component="
              << StringVectorGAP(e_extension.u_component) << " norm=" << e_norm
              << " res_norm=" << res_norm << " ----------\n";
#endif

    std::optional<MyVector<Tint>> opt_v = get_next(e_extension);
    if (opt_v) {
      MyVector<Tint> const &alpha = *opt_v;
#ifdef DEBUG_EDGEWALK_GENERIC
      std::cerr << "alpha=" << StringVectorGAP(alpha) << "\n";
#endif
      std::vector<MyVector<Tint>> l_roots = l_ui;
      l_roots.push_back(alpha);
      MyMatrix<T> Mat_root =
          UniversalMatrixConversion<T, Tint>(MatrixFromVectorFamily(l_roots));
      MyMatrix<T> EquaMat = Mat_root * G;
      MyMatrix<T> NSP = NullspaceTrMat(EquaMat);
      if (NSP.rows() != 1) {
        std::cerr << "We should have exactly one row\n";
        throw TerminalException{1};
      }
      MyVector<T> gen = GetMatrixRow(NSP, 0);
#ifdef DEBUG_EDGEWALK_GENERIC
      std::cerr << "gen=" << StringVectorGAP(gen) << " k=" << StringVectorGAP(k)
                << "\n";
#endif
      T scal = gen.dot(G * k);
      auto get_gen = [&]() -> std::optional<MyVector<T>> {
        if (scal < 0) { // The sign convention means that two vectors in the
                        // same cone have negative scalar product.
          return gen;
        }
        if (scal > 0) {
          return -gen;
        }
        return {};
      };
      std::optional<MyVector<T>> opt_k_new = get_gen();
      if (opt_k_new) {
        MyVector<T> k_new = RemoveFractionVector(*opt_k_new);
        T scal = v_disc_t.dot(G * k_new);
        if (scal < 0) { // The convention in Lorentzian is negative scalar (see
                        // end of Sect 2 of edgewalk paper)
          MyMatrix<Tint> MatRoot = MatrixFromVectorFamily(l_roots);
          FundDomainVertex<T, Tint> fund_v{k_new, MatRoot};
          RootCandidate<T, Tint> eCand =
              gen_possible_extension(G, k, alpha, res_norm, e_norm, fund_v);
          l_candidates.push_back(eCand);
        }
      }
    }
  }
#ifdef TIMINGS
  std::cerr << "Timing |l_candidates|=" << time << "\n";
#endif
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "EdgewalkProcedure : |l_candidates|=" << l_candidates.size()
            << "\n";
#endif
  if (l_candidates.size() > 0) {
    RootCandidate<T, Tint> best_cand = get_best_candidate(l_candidates);
    return best_cand.fund_v;
  }
#ifdef DEBUG_EDGEWALK_GENERIC
  std::cerr << "         --------------- Looking for an isotropic vector "
               "------------\n";
#endif
  // So, no candidates were found. We need to find isotropic vectors.
  const MyMatrix<T> Gred = Pplane * G * Pplane.transpose();
  std::vector<MyVector<T>> BasisIsotrop = GetBasisIsotropicVectors(Gred);
#ifdef TIMINGS
  std::cerr << "Timing |Factor_opt|=" << time << "\n";
#endif
  // We want a vector inside of the cone (there are two: C and -C)
  auto get_can_gen = [&](MyVector<T> const &v) -> MyVector<T> {
    T scal = k.dot(G * v);
#ifdef DEBUG_EDGEWALK_GENERIC
    std::cerr << "  scal=" << scal << "\n";
#endif
    if (scal < 0) {
      // The value should be negative because with the chosen convention,
      // interior vectors have negative pairwise scalar products
      return v;
    }
    if (scal > 0)
      return -v;
    std::cerr << "k=" << StringVectorGAP(k) << " v=" << StringVectorGAP(v)
              << "\n";
    std::cerr << "We should have scal != 0 to be able to conclude\n";
    throw TerminalException{1};
  };
  std::vector<MyVector<T>> l_gens;
  for (auto &U : BasisIsotrop) {
    T sum = U.dot(Gred * U);
    MyVector<T> gen = Pplane.transpose() * U;
    T sum_B = gen.dot(G * gen);
    if (!IsVectorMultiple(gen, k)) {
      MyVector<T> can_gen = get_can_gen(gen);
      T scal = v_disc_t.dot(G * can_gen);
      if (scal < 0) {
        // Convention is negative scalar in Lorentzian theory (see
        // end of sect 2 of edgewalk paper)
        l_gens.push_back(can_gen);
      }
    }
  }
  if (l_gens.size() != 1) {
    std::cerr << "We should have just one vector in order to conclude. Rethink "
                 "needed\n";
    throw TerminalException{1};
  }
  const MyVector<T> &k_new = l_gens[0];
  CuspidalRequest<T, Tint> eReq{l_ui, k_new, k};
#ifdef TIMINGS
  std::cerr << "Timing |CuspidalRequest|=" << time << "\n";
#endif
  std::vector<MyVector<Tint>> l_roots_ret =
      DetermineRootsCuspidalCase_Memoized<T, Tint, Tgroup>(cusp_bank, si, eReq);
  return {RemoveFractionVector(k_new), MatrixFromVectorFamily(l_roots_ret)};
}

template <typename Tint> TheHeuristic<Tint> GetHeuristicIdealStabEquiv() {
  // Allowed returned values: "linalg", "orbmin"
  // Allowed input values "groupsize", "size"
  std::vector<std::string> ListString = {
      "1", "1 groupsize > 500000 size > 100 linalg", "linalg"};
  return HeuristicFromListString<Tint>(ListString);
}

template <typename Tint>
TheHeuristic<Tint> GetHeuristicTryTerminateDualDescription() {
  // Allowed returned values: "trydualdesc", "notry"
  // Allowed input values: "increase_treat_nothingnew"
  std::vector<std::string> ListString = {
      "1", "1 increase_treat_nothingnew > 10 notry", "notry"};
  return HeuristicFromListString<Tint>(ListString);
}

template <typename T, typename Tint, typename Tgroup>
FundDomainVertex_FullInfo<T, Tint, Tgroup>
gen_fund_domain_fund_info(CuspidalBank<T, Tint> &cusp_bank,
                          SublattInfos<T> const &si,
                          FundDomainVertex<T, Tint> const &vert,
                          TheHeuristic<Tint> const &HeuristicIdealStabEquiv) {
  MyMatrix<T> const &G = si.G;
#ifdef TIMINGS
  Microsecond time;
#endif
  //
  // Put the stuff that can help for invariant first
  InitialComputation<T, Tint> ic = GetInitialComputation(G, vert);
  std::string method = "extendedvectfamily";
  if (ic.norm == 0) {
    // In isotropic case, we are unfortunately forced to do more complex stuff
    // Those needs
    ret_type<T, Tint, Tgroup> erec =
        get_canonicalized_record<T, Tint, Tgroup>(ic.ListMat, ic.map_v);
    // Add new vertices
    MyMatrix<T> FAC = UniversalMatrixConversion<T, Tint>(erec.MatRoot);
    MyMatrix<T> FACred = ColumnReduction(FAC);
    std::map<std::string, Tint> mapV;
    mapV["groupsize"] = UniversalScalarConversion<Tint, typename Tgroup::Tint>(
        erec.GRP1.size());
    mapV["size"] = vert.MatRoot.rows();
    std::string choice = HeuristicEvaluation(mapV, HeuristicIdealStabEquiv);
    if (choice == "orbmin") {
      vectface vf = lrs::DualDescription_temp_incd(FACred);
      vectface vf_min = OrbitSplittingSet_GetMinimalOrbit(vf, erec.GRP1);
      for (auto &eFAC : vf_min) {
        AdjacencyDirection<Tint> ad = GetAdjacencyDirection(erec.MatRoot, eFAC);
        FundDomainVertex<T, Tint> fVert =
            EdgewalkProcedure<T, Tint, Tgroup>(cusp_bank, si, vert.gen, ad);
        MyVector<Tint> fVert_tint =
            UniversalVectorConversion<Tint, T>(RemoveFractionVector(fVert.gen));
        ic.map_v[fVert_tint] = 3;
      }
    } else {
      method = "isotropstabequiv";
    }
  }
  ret_type<T, Tint, Tgroup> frec =
      get_canonicalized_record<T, Tint, Tgroup>(ic.ListMat, ic.map_v);
#ifdef TIMINGS
  std::cerr << "Timing |gen_fund_domain_fund_info|=" << time << "\n";
#endif
  return get_full_info(vert, frec, method);
}

template <typename T, typename Tint> struct ResultEdgewalk {
  std::vector<MyMatrix<Tint>> l_gen_isom_cox;
  std::vector<FundDomainVertex<T, Tint>> l_orbit_vertices;
  LorentzianFinitenessGroupTester<T, Tint> group_tester;
  std::optional<bool> is_reflective;
};

template <typename Tint, typename Titer_root, typename Titer_isom>
std::vector<MyVector<Tint>> compute_full_root_orbit_iter(
    Titer_root const &iter_root_begin, Titer_root const &iter_root_end,
    Titer_isom const &iter_isom_begin, Titer_isom const &iter_isom_end) {
  std::unordered_set<MyVector<Tint>> TotalList;
  auto f_insert = [&](MyVector<Tint> const &v) -> void {
    if (TotalList.count(v) != 0)
      return;
    std::unordered_set<MyVector<Tint>> s_v;
    std::vector<MyVector<Tint>> l_v;
    auto f_ins = [&](MyVector<Tint> const &w) -> void {
      if (s_v.count(w) != 0)
        return;
      s_v.insert(w);
      l_v.push_back(w);
      TotalList.insert(w);
    };
    f_ins(v);
    size_t pos = 0;
    while (true) {
      size_t len = l_v.size();
      if (pos == len)
        break;
      Titer_isom iter_isom = iter_isom_begin;
      while (iter_isom != iter_isom_end) {
        for (size_t i = pos; i < len; i++) {
          MyVector<Tint> w_img = iter_isom->transpose() * l_v[i];
          f_ins(w_img);
        }
        iter_isom++;
      }
      pos = len;
    }
  };
  Titer_root iter_root = iter_root_begin;
  while (iter_root != iter_root_end) {
    f_insert(*iter_root);
    iter_root++;
  }
  std::vector<MyVector<Tint>> l_root;
  for (auto &v : TotalList)
    l_root.push_back(v);
  return l_root;
}

template <typename T, typename Tint>
std::vector<MyVector<Tint>>
compute_full_root_orbit(ResultEdgewalk<T, Tint> const &re) {
  std::unordered_set<MyVector<Tint>> s_root;
  for (auto &fdv : re.l_orbit_vertices) {
    size_t len = fdv.MatRoot.rows();
    for (size_t i = 0; i < len; i++) {
      MyVector<Tint> e_root = GetMatrixRow(fdv.MatRoot, i);
      s_root.insert(e_root);
    }
  }
  return compute_full_root_orbit_iter<Tint>(s_root.begin(), s_root.end(),
                                            re.l_gen_isom_cox.begin(),
                                            re.l_gen_isom_cox.end());
}

template <typename T, typename Tint>
std::vector<T> get_list_norms(MyMatrix<T> const &G,
                              ResultEdgewalk<T, Tint> const &re) {
  std::set<T> set_norms;
  auto proc_vertex = [&](FundDomainVertex<T, Tint> const &vert) -> void {
    size_t len = vert.MatRoot.rows();
    for (size_t i = 0; i < len; i++) {
      MyVector<Tint> root = GetMatrixRow(vert.MatRoot, i);
      MyVector<T> root_t = UniversalVectorConversion<T, Tint>(root);
      T norm = root_t.dot(G * root_t);
      set_norms.insert(norm);
    }
  };
  for (auto &evert : re.l_orbit_vertices)
    proc_vertex(evert);
  std::vector<T> l_norms;
  for (auto &v : set_norms)
    l_norms.push_back(v);
  return l_norms;
}

template <typename T, typename Tint>
void PrintResultEdgewalk(MyMatrix<T> const &G,
                         ResultEdgewalk<T, Tint> const &re, std::ostream &os,
                         const std::string &OutFormat,
                         bool const &ComputeAllSimpleRoots) {
  std::vector<T> l_norms = get_list_norms(G, re);
  std::vector<MyVector<Tint>> l_simple_root;
  if (ComputeAllSimpleRoots)
    l_simple_root = compute_full_root_orbit(re);
  std::cerr << "We write G\n";
  std::cerr << "We write l_norms\n";
  if (re.is_reflective) {
    bool val = *re.is_reflective;
    if (val) {
      std::cerr << "lattice found to be reflective\n";
    } else {
      std::cerr << "lattice found NOT to be reflective\n";
    }
  } else {
    std::cerr << "No reflectivity computation\n";
  }
  std::cerr << "We have |l_gen_isom_cox|=" << re.l_gen_isom_cox.size() << "\n";
  size_t n_orbit_vertices = re.l_orbit_vertices.size();
  std::cerr << "We have |l_orbit_vertices|=" << n_orbit_vertices << "\n";
  size_t n_simple = l_simple_root.size();
  std::cerr << "ComputeAllSimpleRoots=" << ComputeAllSimpleRoots
            << " n_simple=" << n_simple << "\n";
  if (OutFormat == "GAP") {
    os << "return rec(LorMat:=";
    WriteMatrixGAP(os, G);
    os << ", l_norms:=";
    WriteStdVectorGAP(os, l_norms);
    os << ", ListIsomCox:=";
    WriteVectorMatrixGAP(os, re.l_gen_isom_cox);
    if (re.is_reflective) {
      bool val = *re.is_reflective;
      if (val) {
        os << ", is_reflective:=true";
      } else {
        os << ", is_reflective:=false";
      }
    }
    os << ", group_tester:=rec(InvariantBasis:=";
    WriteMatrixGAP(os, re.group_tester.get_invariant_basis());
    os << ", max_finite_order:=" << re.group_tester.get_max_finite_order();
    DiagSymMat<T> DiagInfo = re.group_tester.get_diag_info();
    os << ", nbZero:=" << DiagInfo.nbZero;
    os << ", nbPlus:=" << DiagInfo.nbPlus;
    os << ", nbMinus:=" << DiagInfo.nbMinus;
    os << ", is_finite:=" << re.group_tester.get_finiteness_status() << ")";
    os << ", ListVertices:=[";
    bool IsFirst = true;
    for (size_t i = 0; i < n_orbit_vertices; i++) {
      if (!IsFirst)
        os << ",\n";
      IsFirst = false;
      const FundDomainVertex<T, Tint> &evert = re.l_orbit_vertices[i];
      WriteFundDomainVertex(G, evert, os, OutFormat);
    }
    os << "], n_orbit_vertices:=" << n_orbit_vertices;
    if (ComputeAllSimpleRoots) {
      os << ", ListSimpleRoots:=[";
      for (size_t i = 0; i < n_simple; i++) {
        if (i > 0)
          os << ",";
        os << StringVectorGAP(l_simple_root[i]);
      }
      os << "], n_simple:=" << n_simple;
    }
    os << ");\n";
    return;
  }
  if (OutFormat == "TXT") {
    os << "List of found generators of Isom / Cox\n";
    return;
  }
  std::cerr << "Failed to find a matching entry in PrintResultEdgewalk\n";
  std::cerr << "OutFormat=" << OutFormat
            << " but allowed values are GAP and TXT\n";
  throw TerminalException{1};
}

template <typename T> std::string StringStdVectorGAP(std::vector<T> const &V) {
  std::ostringstream os;
  os << "[";
  bool IsFirst = true;
  for (auto &val : V) {
    if (!IsFirst)
      os << ",";
    IsFirst = false;
    os << val;
  }
  os << "]";
  return os.str();
}

template <typename T, typename Tint, typename Tgroup, typename Fvertex,
          typename Fisom, typename Fincrease>
void LORENTZ_RunEdgewalkAlgorithm_Kernel(
    SublattInfos<T> const &si, FundDomainVertex<T, Tint> const &eVert,
    Fvertex f_vertex, Fisom f_isom, Fincrease f_increase_nbdone,
    TheHeuristic<Tint> const &HeuristicIdealStabEquiv) {
  MyMatrix<T> const &G = si.G;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  CuspidalBank<T, Tint> cusp_bank;
  std::vector<int> l_status;
  std::vector<FundDomainVertex_FullInfo<T, Tint, Tgroup>> l_orbit_vertices;
#ifdef TRACK_INFOS_LOG
  std::cout << "return [\n";
#endif
  size_t nbDone = 0;
  auto func_insert_vertex =
      [&](FundDomainVertex_FullInfo<T, Tint, Tgroup> &vertFull1) -> bool {
    size_t len = l_orbit_vertices.size();
#ifdef TIMINGS
    MicrosecondTime time;
#endif
    for (size_t i = 0; i < len; i++) {
      const FundDomainVertex_FullInfo<T, Tint, Tgroup> &vertFull2 =
          l_orbit_vertices[i];
      if (vertFull1.hash == vertFull2.hash) {
        std::optional<MyMatrix<T>> equiv_opt =
            LORENTZ_TestEquivalence<T, Tint, Tgroup>(G, vertFull1, G,
                                                     vertFull2);
        if (equiv_opt) {
#ifdef DEBUG_ENUM_PROCESS
          std::cerr << "Find some isomorphism\n";
#endif
#ifdef TIMINGS
          std::cerr << "Timing |func_insert_vertex(find iso)|=" << time << "\n";
#endif
          bool test = f_isom(UniversalMatrixConversion<Tint, T>(*equiv_opt));
          if (test) {
#ifdef DEBUG_ENUM_PROCESS
            std::cerr << "Exiting at f_isom in LORENTZ_TestEquivalence, return "
                         "true\n";
#endif
            return true;
          } else {
#ifdef DEBUG_ENUM_PROCESS
            std::cerr << "Exiting at f_isom in LORENTZ_TestEquivalence, return "
                         "false\n";
#endif
            return false;
          }
        }
      }
    }
#ifdef TIMINGS
    std::cerr << "Timing |func_insert_vertex(no iso)|=" << time << "\n";
#endif
#ifdef DEBUG_ENUM_PROCESS
    std::cerr << "               Failed to find some isomorphism\n";
    std::cerr << "Before the LORENTZ_GetStabilizerGenerator nbDone=" << nbDone
              << " |l_orbit_vertices|=" << l_orbit_vertices.size() << "\n";
#endif
    std::vector<Telt> LGenIntegral;
    for (auto &eGen_Mat :
         LORENTZ_GetStabilizerGenerator<T, Tint, Tgroup>(G, vertFull1)) {
      bool test = f_isom(UniversalMatrixConversion<Tint, T>(eGen_Mat));
      if (test) {
#ifdef DEBUG_ENUM_PROCESS
        std::cerr << "Exiting at f_isom in func_insert_vertex\n";
#endif
        return true;
      }
      std::optional<std::vector<Tidx>> opt =
          RepresentVertexPermutationTest<Tint, T, Tidx>(
              vertFull1.vert.MatRoot, vertFull1.vert.MatRoot, eGen_Mat);
      if (!opt) {
#ifdef DEBUG_ENUM_PROCESS
        std::cerr << "Failed to find the representation\n";
#endif
        throw TerminalException{1};
      }
      LGenIntegral.push_back(Telt(*opt));
    }
    vertFull1.GRP1_integral =
        Tgroup(LGenIntegral, vertFull1.vert.MatRoot.rows());
#ifdef TIMINGS
    std::cerr << "Timing |Automorphism|=" << time << "\n";
#endif
    bool test = f_vertex(vertFull1);
    if (test) {
#ifdef DEBUG_ENUM_PROCESS
      std::cerr << "Exiting at f_vertex in func_insert_vertex\n";
#endif
      return true;
    }
    l_status.push_back(1);
    l_orbit_vertices.emplace_back(std::move(vertFull1));
#ifdef DEBUG_ENUM_PROCESS
    std::cerr << "Exiting the func_insert_vertex\n";
#endif
#ifdef TIMINGS
    std::cerr << "Timing |func_insert_vertex(end)|=" << time << "\n";
#endif
    return false;
  };
  // We have to do a copy of the Vert since when the vector is extended the
  // previous entries are desttroyed when a new array is built. This would then
  // invalidates a const& theVert reference. Took 1 week to fully debug that
  // problem.
  auto insert_adjacent_vertices =
      [&](FundDomainVertex_FullInfo<T, Tint, Tgroup> const &vertFull) -> bool {
#ifdef TIMINGS
    MicrosecondTime time;
#endif
    const FundDomainVertex<T, Tint> &theVert = vertFull.vert;
#ifdef DEBUG_ENUM_PROCESS
    std::cerr << "insert_edges_from_vertex theVert="
              << StringVectorGAP(RemoveFractionVector(theVert.gen)) << "\n";
#endif
    MyMatrix<T> FAC = UniversalMatrixConversion<T, Tint>(theVert.MatRoot);
    MyMatrix<T> FACred = ColumnReduction(FAC);
    vectface vf = lrs::DualDescription_temp_incd(FACred);
    vectface vf_orb = OrbitSplittingSet(vf, vertFull.GRP1_integral);
#ifdef TIMINGS
    std::cerr << "Timing |vf_orb|=" << time << "\n";
#endif
    //
#ifdef DEBUG_ENUM_PROCESS
    std::cerr << "nbDone=" << nbDone << " |vf_orb|=" << vf_orb.size()
              << " |GRP1|=" << vertFull.GRP1.size()
              << " |GRP1_int|=" << vertFull.GRP1_integral.size() << "\n";
#endif
    for (auto &eFAC : vf_orb) {
      AdjacencyDirection<Tint> ad =
          GetAdjacencyDirection(theVert.MatRoot, eFAC);
      FundDomainVertex<T, Tint> fVert =
          EdgewalkProcedure<T, Tint, Tgroup>(cusp_bank, si, theVert.gen, ad);
      { // Output. Fairly important to see what is happening
#ifdef DEBUG_ENUM_PROCESS
        T norm = fVert.gen.dot(G * fVert.gen);
        std::cerr << "Result of EdgewalkProcedure\n";
        std::cerr << "k=" << StringVectorGAP(theVert.gen) << " l_ui=";
        PrintAdjacencyDirection(std::cerr, ad);
        std::cerr << " fVert=" << StringVectorGAP(fVert.gen) << " norm=" << norm
                  << "\n";
#endif
#ifdef TRACK_INFOS_LOG
        std::vector<T> const &l_norms = si.l_norms;
        std::cout << "rec(k1:=" << StringFundDomainVertexGAP(theVert)
                  << ", k2:=" << StringFundDomainVertexGAP(fVert)
                  << ", ad:=" << StringAdjacencyDirectionGAP(ad)
                  << ", G:=" << StringMatrixGAP(G)
                  << ", l_norms:=" << StringStdVectorGAP(l_norms) << "),\n";
#endif
      }
      FundDomainVertex_FullInfo<T, Tint, Tgroup> fVertFull =
          gen_fund_domain_fund_info<T, Tint, Tgroup>(cusp_bank, si, fVert,
                                                     HeuristicIdealStabEquiv);
      bool test = func_insert_vertex(fVertFull);
      if (test) {
#ifdef DEBUG_ENUM_PROCESS
        std::cerr
            << "Exiting at func_insert_vertex in insert_adjacent_vertices\n";
#endif
        return true;
      }
    }
#ifdef TIMINGS
    std::cerr << "Timing |process vf_orb|=" << time << "\n";
#endif
#ifdef DEBUG_ENUM_PROCESS
    std::cerr << "Exiting from the insert_edges_from_vertex\n";
#endif
    return false;
  };
  FundDomainVertex_FullInfo<T, Tint, Tgroup> eVertFull =
      gen_fund_domain_fund_info<T, Tint, Tgroup>(cusp_bank, si, eVert,
                                                 HeuristicIdealStabEquiv);
  bool test = func_insert_vertex(eVertFull);
  if (test) {
#ifdef DEBUG_ENUM_PROCESS
    std::cerr << "Exiting at initial func_insert_vertex\n";
#endif
    return;
  }
  while (true) {
    bool IsFinished = true;
    size_t len = l_status.size();
    for (size_t i = 0; i < len; i++) {
      if (l_status[i] == 1) {
        nbDone++;
        IsFinished = false;
        l_status[i] = 0;
        // We need to do a direct copy in that case.
        // This is because we are working with a std::vector, which contains a
        // T* array. That array got deallocated and reallocated. This
        // invalidates the const& reference.
        //
        // See for details
        // https://stackoverflow.com/questions/6438086/iterator-invalidation-rules-for-c-containers
        // Which writes: "vector: all iterators and references before the point
        // of insertion are unaffected, unless the new container size is greater
        // than the previous capacity (in which case all iterators and
        // references are invalidated) [23.2.4.3/1]
        //
        // The original problem originally took one week to debug.
        FundDomainVertex_FullInfo<T, Tint, Tgroup> VertFullCp =
            DirectCopy(l_orbit_vertices[i]);
        bool test1 = insert_adjacent_vertices(VertFullCp);
        if (test1) {
#ifdef DEBUG_ENUM_PROCESS
          std::cerr << "Exiting after insert_adjacent_vertices\n";
#endif
          return;
        }
        bool test2 = f_increase_nbdone();
        if (test2) {
#ifdef DEBUG_ENUM_PROCESS
          std::cerr << "Exiting after f_increase_nbdone\n";
#endif
          return;
        }
      }
    }
    if (IsFinished) {
#ifdef DEBUG_ENUM_PROCESS
      std::cerr << "Exiting because all orbits have been treated\n";
#endif
      break;
    }
  }
#ifdef DEBUG_ENUM_PROCESS
  std::cerr
      << "Exiting from the infinite loop of enumeration of vertex pairs\n";
#endif
}

template <typename T, typename Tint, typename Tgroup>
ResultEdgewalk<T, Tint> LORENTZ_RunEdgewalkAlgorithm(
    SublattInfos<T> const &si, FundDomainVertex<T, Tint> const &eVert,
    bool EarlyTerminationIfNotReflective,
    TheHeuristic<Tint> const &HeuristicIdealStabEquiv,
    TheHeuristic<Tint> const &HeuristicTryTerminateDualDescription) {
  MyMatrix<T> const &G = si.G;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<FundDomainVertex<T, Tint>> l_orbit_vertices_ret;
  int dim = G.rows();
  std::unordered_set<MyMatrix<Tint>> s_gen_isom_cox;
  std::unordered_set<MyVector<Tint>> s_simple_roots;

  MyMatrix<Tint> IdMat = IdentityMat<Tint>(dim);
  std::optional<bool> is_reflective;
  LorentzianFinitenessGroupTester<T, Tint> group_tester =
      LorentzianFinitenessGroupTester<T, Tint>(G);
  if (EarlyTerminationIfNotReflective) {
    is_reflective = true;
  }
  int nonew_nbdone = 0;
  size_t n_simple_roots = 0;
  auto f_try_terminate = [&]() -> bool {
    std::vector<MyVector<Tint>> LVect = compute_full_root_orbit_iter<Tint>(
        s_simple_roots.begin(), s_simple_roots.end(), s_gen_isom_cox.begin(),
        s_gen_isom_cox.end());
    size_t n_simple_roots_new = LVect.size();
    if (n_simple_roots_new > n_simple_roots) {
      n_simple_roots = n_simple_roots_new;
      std::unordered_map<MyVector<Tint>, Tidx> MapVectRev;
      for (Tidx i = 0; i < Tidx(n_simple_roots); i++)
        MapVectRev[LVect[i]] = i;
      std::vector<Telt> LGenPerm;
      for (auto &eGen : s_gen_isom_cox) {
        std::vector<Tidx> v(n_simple_roots);
        for (size_t i = 0; i < n_simple_roots; i++) {
          MyVector<Tint> Vimg = eGen.transpose() * LVect[i];
          v[i] = MapVectRev[Vimg];
        }
        Telt eGenPerm(std::move(v));
        LGenPerm.push_back(eGenPerm);
      }
      Tgroup GRP(LGenPerm, n_simple_roots);
      MyMatrix<T> ListIneq =
          -UniversalMatrixConversion<T, Tint>(MatrixFromVectorFamily(LVect)) *
          G;
      vectface vf = DualDescriptionStandard(ListIneq, GRP);
      bool AllRaysInside = true;
      for (auto &eFace : vf) {
        MyVector<T> V = FindFacetInequality(ListIneq, eFace);
        T scal = V.dot(G * V);
        if (scal > 0)
          AllRaysInside = false;
      }
      return AllRaysInside;
    }
    return false;
  };
  auto f_maybe_terminate = [&]() -> bool {
    std::map<std::string, Tint> mapV;
    mapV["increase_treat_nothingnew"] =
        UniversalScalarConversion<Tint, int>(nonew_nbdone);
    std::string choice =
        HeuristicEvaluation(mapV, HeuristicTryTerminateDualDescription);
    if (choice == "") {
      return f_try_terminate();
    }
    return false;
  };
  auto f_increase_nbdone = [&]() -> bool {
    nonew_nbdone++;
    return f_maybe_terminate();
  };
  auto f_vertex =
      [&](FundDomainVertex_FullInfo<T, Tint, Tgroup> const &vertFull) -> bool {
    nonew_nbdone = 0;
    l_orbit_vertices_ret.push_back(vertFull.vert);
    int n_rows = vertFull.vert.MatRoot.rows();
    for (int i_row = 0; i_row < n_rows; i_row++) {
      MyVector<Tint> V = GetMatrixRow(vertFull.vert.MatRoot, i_row);
      s_simple_roots.insert(V);
    }
    return f_maybe_terminate();
  };
  auto f_isom = [&](MyMatrix<Tint> const &eP) -> bool {
    if (eP == IdMat)
      return false;
    if (s_gen_isom_cox.count(eP) > 0)
      return false;
    s_gen_isom_cox.insert(eP);
#ifdef TRACK_INFOS_LOG
    std::cout << "rec(isom:=" << StringMatrixGAP(eP) << "),\n";
#endif
    group_tester.GeneratorUpdate(eP);
    if (is_reflective) {
      if (!group_tester.get_finiteness_status()) {
        is_reflective = false;
        return true;
      }
    }
    return false;
  };
  LORENTZ_RunEdgewalkAlgorithm_Kernel<T, Tint, Tgroup, decltype(f_vertex),
                                      decltype(f_isom),
                                      decltype(f_increase_nbdone)>(
      si, eVert, f_vertex, f_isom, f_increase_nbdone, HeuristicIdealStabEquiv);
  std::vector<MyMatrix<Tint>> l_gen_isom_cox;
  for (auto &e_gen : s_gen_isom_cox)
    l_gen_isom_cox.push_back(e_gen);
  return {std::move(l_gen_isom_cox), std::move(l_orbit_vertices_ret),
          group_tester, is_reflective};
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>> LORENTZ_RunEdgewalkAlgorithm_Isomorphism(
    SublattInfos<T> const &si1, MyMatrix<T> const &G2,
    FundDomainVertex<T, Tint> const &eVert1,
    FundDomainVertex<T, Tint> const &eVert2,
    TheHeuristic<Tint> const &HeuristicIdealStabEquiv) {
  CuspidalBank<T, Tint> cusp_bank;
  std::optional<MyMatrix<Tint>> answer;
  //
  FundDomainVertex_FullInfo<T, Tint, Tgroup> vertFull2 =
      gen_fund_domain_fund_info<T, Tint, Tgroup>(cusp_bank, si1, eVert2,
                                                 HeuristicIdealStabEquiv);
  auto f_increase_nbdone = [&]() -> bool { return false; };
  auto f_vertex =
      [&](FundDomainVertex_FullInfo<T, Tint, Tgroup> const &vertFull1) -> bool {
    if (vertFull1.hash == vertFull2.hash) {
      std::optional<MyMatrix<T>> equiv_opt =
          LORENTZ_TestEquivalence(si1.G, vertFull1, G2, vertFull2);
      if (equiv_opt) {
        answer = UniversalMatrixConversion<Tint, T>(*equiv_opt);
        return true;
      }
    }
    return false;
  };
  auto f_isom = [&](MyMatrix<Tint> const &eP) -> bool { return false; };
  LORENTZ_RunEdgewalkAlgorithm_Kernel<T, Tint, Tgroup, decltype(f_vertex),
                                      decltype(f_isom),
                                      decltype(f_increase_nbdone)>(
      si1, eVert1, f_vertex, f_isom, f_increase_nbdone,
      HeuristicIdealStabEquiv);
  return answer;
}

template <typename T, typename Tint>
std::vector<MyVector<Tint>>
get_simple_cone_from_lattice(SublattInfos<T> const &si,
                             MyMatrix<Tint> const &NSP_tint) {
  MyMatrix<T> const &G = si.G;
  std::vector<T> const &l_norms = si.l_norms;
  std::cerr << "Beginning of get_simple_cone\n";
  int dimSpace = NSP_tint.rows();
  MyMatrix<T> NSP = UniversalMatrixConversion<T, Tint>(NSP_tint);
  MyMatrix<Tint> G_int = UniversalMatrixConversion<Tint, T>(G);
  std::vector<MyVector<Tint>> l_roots;
  MyVector<T> zeroVect = ZeroVector<T>(dimSpace);
  for (auto &e_norm : l_norms) {
    std::cerr << "-------- e_norm=" << e_norm << " --------------\n";
    MyMatrix<T> const &Latt = si.map_norm_latt.at(e_norm);
    MyMatrix<T> Latt_i_Orth = IntersectionLattice(NSP, Latt);
    MyMatrix<T> G_P = Latt_i_Orth * G * Latt_i_Orth.transpose();
    CheckPositiveDefinite(G_P);
    std::vector<MyVector<Tint>> l_v =
        FindFixedNormVectors<T, Tint>(G_P, zeroVect, e_norm);
    for (auto &e_v : l_v) {
      MyVector<Tint> e_root = Latt_i_Orth_tint.transpose() * e_v;
      std::optional<MyVector<Tint>> opt = SolutionIntMat(NSP_tint, e_root);
      if (opt) {
        l_roots.push_back(*opt);
      } else {
        std::cerr << "Failed to find the solution in the subspace\n";
        throw TerminalException{1};
      }
    }
    std::cerr << "e_norm=" << e_norm << " |l_v|=" << l_v.size() << "\n";
  }
  std::vector<MyVector<Tint>> facet_one_cone = GetFacetOneDomain(l_roots);
  if (facet_one_cone.size() != size_t(dimSpace)) {
    std::cerr << "dimSpace =" << dimSpace
              << " |facet_one_cone|=" << facet_one_cone.size() << "\n";
    std::cerr << "and they should be equal\n";
    throw TerminalException{1};
  }
  std::vector<MyVector<Tint>> l_ui;
  for (auto &e_root : facet_one_cone) {
    MyVector<Tint> v = NSP_tint.transpose() * e_root;
    l_ui.push_back(v);
  }
  return l_ui;
}

template <typename T, typename Tint>
MyMatrix<Tint> get_simple_cone(SublattInfos<T> const &si,
                               MyVector<T> const &V) {
  std::cerr << "------------------------------ get_simple_cone "
               "--------------------------\n";
  MyMatrix<T> const &G = si.G;
  std::cerr << "G=\n";
  WriteMatrixGAP(std::cerr, G);
  std::cerr << "\n";
  T norm = V.dot(G * V);
  std::cerr << "V=" << StringVectorGAP(V) << " norm=" << norm << "\n";
  if (norm > 0) {
    std::cerr << "We need a vector of negative norm or zero norm in order to "
                 "build a system of simple roots\n";
    throw TerminalException{1};
  }
  int dim = G.rows();
  MyVector<T> eProd = G * V;
  MyMatrix<T> eProdB(1, dim);
  AssignMatrixRow(eProdB, 0, eProd);
  MyMatrix<T> NSP = NullspaceIntTrMat(eProdB);
  MyMatrix<Tint> NSP_tint = UniversalMatrixConversion<Tint, T>(NSP);
  if (norm < 0) {
    std::cerr << "Ordinary case\n";
    // ordinary point case
    std::vector<MyVector<Tint>> l_vect =
        get_simple_cone_from_lattice(si, NSP_tint);
    return MatrixFromVectorFamily(l_vect);
  } else {
    std::cerr << "Ideal case\n";
    // ideal point case
    MyVector<Tint> V_i =
        UniversalVectorConversion<Tint, T>(RemoveFractionVector(V));
    std::optional<MyVector<Tint>> opt = SolutionIntMat(NSP_tint, V_i);
    if (!opt) {
      std::cerr << "The vector V does not below to NSP which contradicts it "
                   "being isotrop\n";
      throw TerminalException{1};
    }
    MyVector<Tint> const &Vnsp = *opt;
    std::cerr << "We have Vnsp\n";
    /*
      We need a more general code for finding complement of subspace, possibly
      using HermiteNormalForm
     */
    MyMatrix<Tint> Basis = ComplementToBasis(Vnsp);
    MyMatrix<Tint> Basis_NSP = Basis * NSP_tint;
    MyMatrix<T> Subspace = UniversalMatrixConversion<T, Tint>(Basis_NSP);
    std::map<T, size_t> MapIdxFr;
    std::vector<LatticeProjectionFramework<T>> ListFr;
    MyVector<T> zeroVect = ZeroVector<T>(Subspace.rows());
    std::vector<MyVector<T>> list_vect;
    std::vector<MyVector<T>> list_vect_big;
    std::vector<T> list_norm;
    size_t pos = 0;
    for (auto &e_norm : si.l_norms) {
      std::cerr << "e_norm=" << e_norm << "\n";
      MyMatrix<T> const &Latt = si.map_norm_latt.at(e_norm);
      MyMatrix<T> Latt_inter_NSP_pre = IntersectionLattice(Latt, NSP);
      MyMatrix<T> Latt_inter_NSP =
          LLLbasisReduction<T, Tint>(Latt_inter_NSP_pre).LattRed;
      LatticeProjectionFramework<T> fr(G, Subspace, Latt_inter_NSP);
      MapIdxFr[e_norm] = pos;
      ListFr.push_back(fr);
      //
      MyMatrix<T> const &RelBasis = fr.BasisProj;
      MyMatrix<T> G_P = RelBasis * G * RelBasis.transpose();
      //      CheckPositiveDefinite(G_P);
      std::vector<MyVector<Tint>> l_v =
          FindFixedNormVectors<T, Tint>(G_P, zeroVect, e_norm);
      std::cerr << "|l_v|=" << l_v.size() << "\n";
      for (auto &e_v : l_v) {
        MyVector<T> e_vt = UniversalVectorConversion<T, Tint>(e_v);
        MyVector<T> e_vect = RelBasis.transpose() * e_vt;
        std::optional<MyVector<T>> opt = SolutionMat(Subspace, e_vect);
        if (opt) {
          list_vect.push_back(*opt);
          list_vect_big.push_back(e_vect);
          list_norm.push_back(e_norm);
        } else {
          std::cerr << "Failed to find the solution in the subspace\n";
          throw TerminalException{1};
        }
      }
      std::cerr << "Insertion done for all vectors\n";
      pos++;
    }
    if (list_vect.size() == 0) {
      std::cerr << "The list of vectors is empty. Cannot be reflective. In any "
                   "case, we cannot continue\n";
      throw NonReflectivityException{};
    }
    int rnk = RankMat(MatrixFromVectorFamily(list_vect));
    std::cerr << "|list_vect|=" << list_vect.size()
              << " Rank(list_vect)=" << rnk << "\n";
    if (rnk < G.rows() - 2) {
      std::cerr << "The list of roots is not of correct rank. Cannot be "
                   "reflective. In any case, we cannot continue\n";
      throw NonReflectivityException{};
    }
    auto get_one_root = [&](MyVector<T> const &e_vect) -> MyVector<Tint> {
      size_t len = list_vect.size();
      for (size_t i = 0; i < len; i++) {
        MyVector<T> const &f_vect = list_vect[i];
        if (f_vect == e_vect) {
          T const &e_norm = list_norm[i];
          MyVector<T> const &e_vect_big = list_vect_big[i];
          size_t idx = MapIdxFr[e_norm];
          LatticeProjectionFramework<T> const &fr = ListFr[idx];
          std::optional<MyVector<T>> opt = fr.GetOnePreimage(e_vect_big);
          if (!opt) {
            std::cerr << "Failed to find the Preimage\n";
            throw TerminalException{1};
          }
          MyVector<T> const &V = *opt;
          MyVector<Tint> V_i = UniversalVectorConversion<Tint, T>(V);
          return V_i;
        }
      }
      std::cerr << "Failed to find the vector in the list\n";
      throw TerminalException{1};
    };
    std::vector<MyVector<T>> facet_one_cone = GetFacetOneDomain(list_vect);
    std::cerr << "We have facet_one_cone\n";
    std::vector<MyVector<Tint>> l_ui;
    for (auto &e_vt : facet_one_cone) {
      MyVector<Tint> e_vi = get_one_root(e_vt);
      l_ui.push_back(e_vi);
    }
    std::cerr << "We have l_ui\n";
    MyMatrix<T> Pplane = Get_Pplane(G, l_ui);
    std::cerr << "We have Pplane\n";
    auto get_kP = [&]() -> MyVector<T> {
      MyMatrix<T> Gprod = Pplane * G * Pplane.transpose();
      T CritNorm = 0;
      std::cerr << "Gprod=\n";
      WriteMatrix(std::cerr, Gprod);
      MyVector<T> eVect_A = GetNegativeNormVector(Gprod);
      MyVector<T> eVect_B = Pplane.transpose() * eVect_A;
      T scal = V.dot(G * eVect_B);
      if (scal < 0) { // This is because of the sign convention
        return eVect_B;
      } else {
        return -eVect_B;
      }
    };
    MyVector<T> kP = get_kP();
    std::cerr << "we have kP\n";
    CuspidalRequest<T, Tint> eReq{l_ui, V, kP};
    std::vector<MyVector<Tint>> l_vect = DetermineRootsCuspidalCase(si, eReq);
    std::cerr << "get_simple_cone, step 8\n";
    return MatrixFromVectorFamily(l_vect);
  }
}

template <typename T, typename Tint, typename Tgroup>
MyVector<T> GetOneVertex(SublattInfos<T> const &si, bool const &ApplyReduction,
                         std::string const &DualDescProg,
                         bool const &EarlyTerminationIfNotReflective) {
  MyMatrix<T> const &G = si.G;
  std::vector<T> const &l_norms = si.l_norms;
  ResultReduction<T, Tint> ResRed =
      ComputeReductionIndefinite_opt<T, Tint>(G, ApplyReduction);
  /*
    We have ResRed.B and ResRed.Mred    with Mred = B * G * B^T
  */
  VinbergTot<T, Tint> Vtot = GetVinbergFromG<T, Tint>(
      ResRed.Mred, l_norms, DualDescProg, EarlyTerminationIfNotReflective);
  MyVector<Tint> V1 = FindOneInitialRay<T, Tint, Tgroup>(Vtot);
  MyVector<Tint> V2 = ResRed.B.transpose() * V1;
  MyVector<Tint> V3 = RemoveFractionVector(V2);
  MyVector<T> V4 = UniversalVectorConversion<T, Tint>(V3);
  return V4;
}

template <typename T, typename Tint, typename Tgroup>
FundDomainVertex<T, Tint>
get_initial_vertex(SublattInfos<T> const &si, bool const &ApplyReduction,
                   std::string const &DualDescProg,
                   bool const &EarlyTerminationIfNotReflective,
                   std::string const &OptionInitialVertex,
                   std::string const &FileInitialVertex) {
  std::cerr << "Beginning of get_initial_vertex\n";
  if (OptionInitialVertex == "FileVertex") {
    if (!IsExistingFile(FileInitialVertex)) {
      std::cerr << "The file FileInitialVertex=" << FileInitialVertex
                << " is missing\n";
      throw TerminalException{1};
    }
    std::ifstream is(FileInitialVertex);
    MyVector<T> gen = ReadVector<T>(is);
    MyMatrix<Tint> MatRoot = get_simple_cone<T, Tint>(si, gen);
    return {RemoveFractionVector(gen), MatRoot};
  }
  if (OptionInitialVertex == "FileVertexRoots") {
    if (!IsExistingFile(FileInitialVertex)) {
      std::cerr << "The file FileInitialVertex=" << FileInitialVertex
                << " is missing\n";
      throw TerminalException{1};
    }
    std::ifstream is(FileInitialVertex);
    MyVector<T> gen = ReadVector<T>(is);
    MyMatrix<Tint> MatRoot = ReadMatrix<Tint>(is);
    return {RemoveFractionVector(gen), MatRoot};
  }
#ifdef ALLOW_VINBERG_ALGORITHM_FOR_INITIAL_VERTEX
  if (OptionInitialVertex == "vinberg") {
    MyVector<T> V = GetOneVertex<T, Tint, Tgroup>(
        si, ApplyReduction, DualDescProg, EarlyTerminationIfNotReflective);
    MyMatrix<Tint> MatRoot = get_simple_cone<T, Tint>(si, V);
    return {RemoveFractionVector(V), MatRoot};
  }
#endif
  std::cerr << "Failed to find a matching entry in get_initial_vertex\n";
  std::cerr << "OptionInitialVertex=" << OptionInitialVertex
            << " but allowed values are FileVertex or FileVertexRoots\n";
#ifdef ALLOW_VINBERG_ALGORITHM_FOR_INITIAL_VERTEX
  std::cerr << "and vinberg has also been allowed\n";
#else
  std::cerr << "option vinberg has not been allowed\n";
#endif
  throw TerminalException{1};
}

template <typename T> void TestLorentzianity(MyMatrix<T> const &G) {
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(G);
  if (DiagInfo.nbZero != 0 || DiagInfo.nbMinus != 1) {
    std::cerr << "G=\n";
    WriteMatrix(std::cerr, G);
    std::cerr << "We have nbZero=" << DiagInfo.nbZero
              << " nbPlus=" << DiagInfo.nbPlus
              << " nbMinus=" << DiagInfo.nbMinus << "\n";
    std::cerr
        << "In the hyperbolic geometry we should have nbZero=0 and nbMinus=1\n";
    throw TerminalException{1};
  }
}

template <typename T, typename Tint>
void PrintVertexInformation(MyMatrix<T> const &G,
                            FundDomainVertex<T, Tint> const &eVert) {
  T norm = eVert.gen.dot(G * eVert.gen);
  std::cerr << "Initial vertex is eVert=" << StringVectorGAP(eVert.gen)
            << " norm=" << norm << "\n";
  std::cerr << "|MatRoot|=" << eVert.MatRoot.rows() << "\n";
  std::vector<MyVector<Tint>> l_root;
  for (int i = 0; i < eVert.MatRoot.rows(); i++) {
    MyVector<Tint> eLine = GetMatrixRow(eVert.MatRoot, i);
    std::cerr << StringVectorGAP(eLine) << "\n";
    l_root.push_back(eLine);
  }
  std::pair<MyMatrix<T>, MyMatrix<T>> ep = ComputeCoxeterMatrix(G, l_root);
  const MyMatrix<T> &CoxMat = ep.first;
  const MyMatrix<T> &ScalMat = ep.second;
  std::cerr << "ScalMat=\n";
  WriteMatrix(std::cerr, ScalMat);
  std::cerr << "CoxMat=\n";
  WriteMatrix(std::cerr, CoxMat);
  std::cerr << "We have CoxMat\n";
  std::string symb = coxdyn_matrix_to_string(CoxMat);
  std::cerr << "symb=" << symb << "\n";
  std::cerr << "l_roots=\n";
  WriteMatrix(std::cerr, eVert.MatRoot);
}

template <typename T, typename Tint, typename Tgroup>
void MainFunctionEdgewalk(FullNamelist const &eFull) {
  SingleBlock BlockPROC = eFull.ListBlock.at("PROC");
  std::string FileLorMat = BlockPROC.ListStringValues.at("FileLorMat");
  MyMatrix<T> G = ReadMatrixFile<T>(FileLorMat);
  TestLorentzianity(G);
  //
  std::string OptionNorms = BlockPROC.ListStringValues.at("OptionNorms");
  std::string DualDescProg = BlockPROC.ListStringValues.at("DualDescProg");
  bool EarlyTerminationIfNotReflective =
      BlockPROC.ListBoolValues.at("EarlyTerminationIfNotReflective");
  bool ApplyReduction = BlockPROC.ListBoolValues.at("ApplyReduction");
  std::vector<T> l_norms = get_initial_list_norms<T, Tint>(G, OptionNorms);
  SublattInfos<T> si = ComputeSublatticeInfos<T, Tint>(G, l_norms);
  std::cerr << "We have l_norms\n";
  //
  std::string FileHeuristicIdealStabEquiv =
      BlockPROC.ListStringValues.at("FileHeuristicIdealStabEquiv");
  TheHeuristic<Tint> HeuristicIdealStabEquiv =
      GetHeuristicIdealStabEquiv<Tint>();
  ReadHeuristicFileCond(FileHeuristicIdealStabEquiv, HeuristicIdealStabEquiv);
  //
  std::string FileHeuristicTryTerminateDualDescription =
      BlockPROC.ListStringValues.at("FileHeuristicTryTerminateDualDescription");
  TheHeuristic<Tint> HeuristicTryTerminateDualDescription =
      GetHeuristicTryTerminateDualDescription<Tint>();
  ReadHeuristicFileCond(FileHeuristicTryTerminateDualDescription,
                        HeuristicTryTerminateDualDescription);
  //
  auto print_result_edgewalk = [&](ResultEdgewalk<T, Tint> const &re) -> void {
    std::string OutFormat = BlockPROC.ListStringValues.at("OutFormat");
    std::string FileOut = BlockPROC.ListStringValues.at("FileOut");
    bool ComputeAllSimpleRoots =
        BlockPROC.ListBoolValues.at("ComputeAllSimpleRoots");
    std::cerr << "OutFormat=" << OutFormat << " FileOut=" << FileOut
              << " ComputeAllSimpleRoots=" << ComputeAllSimpleRoots << "\n";
    auto f_print = [&](std::ostream &os) -> void {
      PrintResultEdgewalk(G, re, os, OutFormat, ComputeAllSimpleRoots);
    };
    print_stderr_stdout_file(FileOut, f_print);
  };
  //
  std::string OptionInitialVertex =
      BlockPROC.ListStringValues.at("OptionInitialVertex");
  std::string FileInitialVertex =
      BlockPROC.ListStringValues.at("FileInitialVertex");
  try {
    FundDomainVertex<T, Tint> eVert = get_initial_vertex<T, Tint, Tgroup>(
        si, ApplyReduction, DualDescProg, EarlyTerminationIfNotReflective,
        OptionInitialVertex, FileInitialVertex);
#ifdef PRINT_SYMBOL_INFORMATION
    PrintVertexInformation(G, eVert);
#endif
    //
    ResultEdgewalk<T, Tint> re = LORENTZ_RunEdgewalkAlgorithm<T, Tint, Tgroup>(
        si, eVert, EarlyTerminationIfNotReflective, HeuristicIdealStabEquiv,
        HeuristicTryTerminateDualDescription);
    print_result_edgewalk(re);
  } catch (NonReflectivityException const &e) {
    if (!EarlyTerminationIfNotReflective) {
      std::cerr << "The program cannot go forward. Since we have "
                   "EarlyTerminationIfNotReflective = F\n";
      std::cerr << "this is actually a runtime error\n";
    }
    ResultEdgewalk<T, Tint> re{
        {}, {}, LorentzianFinitenessGroupTester<T, Tint>(G), false};
    print_result_edgewalk(re);
  }
  std::cerr << "We are after the PrintResultEdgewalk\n";
}

template <typename T, typename Tint, typename Tgroup>
void MainFunctionEdgewalk_Isomorphism(FullNamelist const &eFull) {
  SingleBlock BlockPROC = eFull.ListBlock.at("PROC");
  std::string FileLorMat1 = BlockPROC.ListStringValues.at("FileLorMat1");
  std::string FileLorMat2 = BlockPROC.ListStringValues.at("FileLorMat2");
  MyMatrix<T> G1 = ReadMatrixFile<T>(FileLorMat1);
  MyMatrix<T> G2 = ReadMatrixFile<T>(FileLorMat2);
  TestLorentzianity(G1);
  TestLorentzianity(G2);
  std::string FileHeuristicIdealStabEquiv =
      BlockPROC.ListStringValues.at("FileHeuristicIdealStabEquiv");
  TheHeuristic<Tint> HeuristicIdealStabEquiv =
      GetHeuristicIdealStabEquiv<Tint>();
  ReadHeuristicFileCond(FileHeuristicIdealStabEquiv, HeuristicIdealStabEquiv);
  //
  auto print_result = [&](std::optional<MyMatrix<Tint>> const &opt) -> void {
    std::string OutFormat = BlockPROC.ListStringValues.at("OutFormat");
    auto print_result_isomorphism = [&](std::ostream &os) -> void {
      if (OutFormat == "GAP") {
        if (opt) {
          os << "return ";
          WriteMatrixGAP(os, *opt);
          os << ";\n";
        } else {
          os << "return fail;\n";
        }
      }
      std::cerr << "We fil to have a matching format. OutFormat=" << OutFormat
                << "\n";
      throw TerminalException{1};
    };
    std::string FileOut = BlockPROC.ListStringValues.at("FileOut");
    print_stderr_stdout_file(FileOut, print_result_isomorphism);
  };
  //
  std::string OptionNorms = "all";
  bool ApplyReduction = BlockPROC.ListBoolValues.at("ApplyReduction");
  std::string DualDescProg = BlockPROC.ListStringValues.at("DualDescProg");
  std::vector<T> l_norms1 = get_initial_list_norms<T, Tint>(G1, OptionNorms);
  std::vector<T> l_norms2 = get_initial_list_norms<T, Tint>(G2, OptionNorms);
  if (l_norms1 != l_norms2) {
    print_result({});
    return;
  }
  std::vector<T> l_norms = l_norms1;
  SublattInfos<T> si1 = ComputeSublatticeInfos<T, Tint>(G1, l_norms);
  SublattInfos<T> si2 = ComputeSublatticeInfos<T, Tint>(G1, l_norms);
  std::cerr << "We have l_norms\n";
  //
  std::string OptionInitialVertex = "vinberg";
  std::string FileInitialVertex = "irrelevant";
  bool EarlyTerminationIfNotReflective = false;
  FundDomainVertex<T, Tint> eVert1 = get_initial_vertex<T, Tint, Tgroup>(
      si1, ApplyReduction, DualDescProg, EarlyTerminationIfNotReflective,
      OptionInitialVertex, FileInitialVertex);
  FundDomainVertex<T, Tint> eVert2 = get_initial_vertex<T, Tint, Tgroup>(
      si2, ApplyReduction, DualDescProg, EarlyTerminationIfNotReflective,
      OptionInitialVertex, FileInitialVertex);
  //
  std::optional<MyMatrix<Tint>> opt =
      LORENTZ_RunEdgewalkAlgorithm_Isomorphism<T, Tint, Tgroup>(
          si1, G2, eVert1, eVert2, HeuristicIdealStabEquiv);
  print_result(opt);
}

// clang-format off
#endif  // SRC_LORENTZIAN_EDGEWALK_H_
// clang-format on

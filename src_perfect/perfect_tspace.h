// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_PERFECT_TSPACE_H_
#define SRC_PERFECT_PERFECT_TSPACE_H_

// clang-format off
#include "LatticeStabEquiCan.h"
#include "OnlineExhaustiveReduction.h"
#include "Temp_PerfectForm.h"
#include "POLY_lrslib.h"
#include "Positivity.h"
#include "POLY_RecursiveDualDesc.h"
#include "POLY_AdjacencyScheme.h"
#include "hash_functions.h"
#include "Tspace_Namelist.h"
// clang-format on

#ifdef DEBUG
#define DEBUG_PERFECT_TSPACE
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_PERFECT_TSPACE
#endif

#ifdef TIMINGS
#define TIMINGS_PERFECT_TSPACE
#endif


void WriteEntryGAP(std::ostream &os, Face const &eFace) {
  os << "[";
  size_t len = eFace.size();
  for (size_t i = 0; i < len; i++) {
    if (i > 0)
      os << ",";
    os << (eFace[i] ? 1 : 0);
  }
  os << "]";
}

void WriteEntryPYTHON(std::ostream &os, Face const &eFace) {
  os << "[";
  size_t len = eFace.size();
  for (size_t i = 0; i < len; i++) {
    if (i > 0)
      os << ",";
    os << (eFace[i] ? 1 : 0);
  }
  os << "]";
}

template <typename T, typename Tint, typename Tgroup>
struct PerfectTspace_Obj {
  MyMatrix<T> Gram;
  Tshortest<T, Tint> RecSHV;
  Tgroup GRP;
};

template <typename T, typename Tint, typename Tgroup>
void check_correctness(PerfectTspace_Obj<T,Tint,Tgroup> const& po) {
  int n_row = po.RecSHV.SHV.rows();
  T min = po.RecSHV.min;
  for (int i_row=0; i_row<n_row; i_row++) {
    MyVector<Tint> V = GetMatrixRow(po.RecSHV.SHV, i_row);
    T norm = EvaluationQuadForm<T,Tint>(po.Gram, V);
    if (min != norm) {
      std::cerr << "PERF: inconsistency at i_row=" << i_row << "\n";
      std::cerr << "PERF: norm=" << norm << " min=" << min << "\n";
      throw TerminalException{1};
    }
  }
}




template <typename T, typename Tint> struct PerfectTspace_AdjI {
  Face eInc;
  MyMatrix<T> Gram;
  Tshortest<T, Tint> RecSHV;
};

template <typename Tint> struct PerfectTspace_AdjO {
  Face eInc;
  MyMatrix<Tint> eBigMat;
};

template <typename T, typename Tint, typename Tgroup>
struct DataPerfectTspace {
  LinSpaceMatrix<T> LinSpa;
  OnlineHierarchicalMatrixReduction<Tint> onl_gens;
  RecordDualDescOperation<T, Tgroup> rddo;
  std::ostream &get_os() { return rddo.os; }
};


template <typename T, typename Tint>
std::pair<MyMatrix<Tint>, MyMatrix<T>> get_reduced_gram(OnlineHierarchicalMatrixReduction<Tint> const& onl_gens, MyMatrix<T> const& g) {
  std::vector<MyMatrix<Tint>> l_gens = onl_gens.get_current_matrix_t();
  std::vector<MyMatrix<T>> l_gens_T = UniversalStdVectorMatrixConversion<T,Tint>(l_gens);
  size_t n_gen = l_gens.size();
  int n = g.rows();
  auto norm=[&](MyMatrix<T> const& M) -> T {
    T e_norm(0);
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        e_norm += T_abs(M(i,j));
      }
    }
    return e_norm;
  };
  MyMatrix<Tint> P = IdentityMat<Tint>(n);
  MyMatrix<T> g_ret = g;
  T norm_ret = norm(g);
  auto atp_improve=[&]() -> bool {
    for (size_t i_gen=0; i_gen<n_gen; i_gen++) {
      MyMatrix<T> const& egen = l_gens_T[i_gen];
      MyMatrix<T> prod = egen * g_ret * egen.transpose();
      T norm_atp = norm(prod);
      if (norm_atp < norm_ret) {
        norm_ret = norm_atp;
        g_ret = prod;
        P = l_gens[i_gen] * P;
        return true;
      }
    }
    return false;
  };
  while(1) {
    bool test = atp_improve();
    if (!test) {
      break;
    }
  }
  return {P, g_ret};
}



template <typename T, typename Tint, typename Tgroup>
std::vector<PerfectTspace_AdjI<T, Tint>>
TSPACE_GetAdjacencies(LinSpaceMatrix<T> const &LinSpa,
                      MyMatrix<T> const &eGram,
                      Tshortest<T, Tint> const &RecSHV,
                      Tgroup const& GRP,
                      std::ostream &os) {
#ifdef TIMINGS_PERFECT_TSPACE
  HumanTime time;
#endif
#ifdef DEBUG_PERFECT_TSPACE
  os << "PERF_TSPACE: TSPACE_GetAdjacencies, begin\n";
#endif
  MyMatrix<T> SHV_T = conversion_and_duplication<T,Tint>(RecSHV.SHV);
#ifdef TIMINGS_PERFECT_TSPACE
  os << "|PERF_TSPACE: GetAdj_shv_t|=" << time << "\n";
#endif
#ifdef SANITY_CHECK_PERFECT_TSPACE
  if (!is_perfect_in_space(LinSpa, RecSHV)) {
    std::cerr << "Error in input of TSPACE_GetAdjacencies, the input matrix is not perfect\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_PERFECT_TSPACE
  MyMatrix<T> ScalMat = get_scal_mat<T,Tint>(LinSpa, RecSHV);
  os << "PERF_TSPACE: The ScalMat is the following\n";
  WriteMatrix(os, ScalMat);
  bool test = is_perfect_in_space<T,Tint>(LinSpa, RecSHV);
  os << "PERF_TSPACE: RankMat(ScalMat)=" << RankMat(ScalMat) << " test=" << test << "\n";
#endif

  RyshkovGRP<T,Tgroup> eCone = GetNakedPerfectCone_GRP<T,Tgroup>(LinSpa,
                                                                 eGram,
                                                                 SHV_T,
                                                                 GRP, os);
#ifdef TIMINGS_PERFECT_TSPACE
  os << "|PERF_TSPACE: GetAdj_eCone|=" << time << "\n";
#endif
  vectface ListIncd = DualDescriptionStandard<T,Tgroup>(eCone.PerfDomEXT, eCone.GRPsub);
#ifdef TIMINGS_PERFECT_TSPACE
  os << "|PERF_TSPACE: GetAdj_ListIncd|=" << time << "\n";
#endif
#ifdef DEBUG_PERFECT_TSPACE_DISABLE
  os << "PERF_TSPACE: The eCone.PerfDomEXT is the following\n";
  WriteMatrix(os, eCone.PerfDomEXT);
  os << "PERF_TSPACE: RankMat(eCone.PerfDomEXT)=" << RankMat(eCone.PerfDomEXT) << "\n";
#endif
  std::vector<PerfectTspace_AdjI<T, Tint>> ListAdj;
#ifdef DEBUG_PERFECT_TSPACE
  os << "PERF_TSPACE: |ListIncd|=" << ListIncd.size() << "\n";
  os << "PERF_TSPACE: |eCone.PerfDomEXT|=" << eCone.PerfDomEXT.cols() << " / " << eCone.PerfDomEXT.rows() << "\n";
  os << "PERF_TSPACE: rk(eCone.PerfDomEXT)=" << RankMat(eCone.PerfDomEXT) << "\n";
  size_t pos = 0;
#endif

  for (auto &eIncd : ListIncd) {
#ifdef DEBUG_PERFECT_TSPACE
    os << "PERF_TSPACE: pos=" << pos << " eIncd=" << eIncd << "\n";
    pos += 1;
#endif
    MyVector<T> eFacet = FindFacetInequality(eCone.PerfDomEXT, eIncd);
    MyMatrix<T> eMatDir = LINSPA_GetMatrixInTspace(LinSpa, eFacet);
    // That seems to slow down things
    //    MyMatrix<T> eMatDirScal = ScalarCanonicalizationMatrix(eMatDir);
#ifdef TIMINGS_PERFECT_TSPACE
    HumanTime time_flip;
#endif
    std::pair<MyMatrix<T>, Tshortest<T, Tint>> pair =
        Flipping_Perfect<T, Tint>(eGram, eMatDir, os);
#ifdef TIMINGS_PERFECT_TSPACE
    os << "|PERF_TSPACE: Flipping_Perfect|=" << time_flip << "\n";
#endif
#ifdef DEBUG_PERFECT_TSPACE
    os << "PERF_TSPACE: Result flipping pos=" << pos << " min=" << pair.second.min << " |SHV|=" << pair.second.SHV.rows() << " det=" << DeterminantMat(pair.first) << "\n";
#endif
#ifdef SANITY_CHECK_PERFECT_TSPACE
    if (!is_perfect_in_space(LinSpa, pair.second)) {
      std::cerr << "The matrix created by Flipping_Perfect is not perfect\n";
      throw TerminalException{1};
    }
#endif

    PerfectTspace_AdjI<T, Tint> eAdj{eIncd, pair.first, pair.second};
    ListAdj.push_back(eAdj);
  }
#ifdef TIMINGS_PERFECT_TSPACE
  os << "|PERF_TSPACE: GetAdj_ListAdj|=" << time << "\n";
#endif
  return ListAdj;
}

template <typename T, typename Tint, typename Tgroup>
struct DataPerfectTspaceFunc {
  DataPerfectTspace<T, Tint, Tgroup> data;
  using Tobj = PerfectTspace_Obj<T, Tint, Tgroup>;
  using TadjI = PerfectTspace_AdjI<T, Tint>;
  using TadjO = PerfectTspace_AdjO<Tint>;
  std::ostream &get_os() { return data.rddo.os; }

  Tobj f_init() {
    std::ostream &os = get_os();
    std::pair<MyMatrix<T>, Tshortest<T, Tint>> pair = GetOnePerfectForm<T,Tint>(data.LinSpa, os);
    Tobj x{pair.first, pair.second, {}};
    return x;
  }

  size_t f_hash(size_t const &seed, Tobj const &x) {
    std::ostream &os = get_os();
#ifdef DEBUG_PERFECT_TSPACE
    os << "PERF_TSPACE: Before f_hash\n";
#endif
    return SimplePerfect_Invariant<T,Tint>(seed, data.LinSpa, x.Gram, x.RecSHV, os);
  }

  std::optional<TadjO> f_repr(Tobj const &x, TadjI const &y) {
    std::ostream &os = get_os();
    std::optional<MyMatrix<Tint>> opt = SimplePerfect_TestEquivalence<T, Tint, Tgroup>(data.LinSpa, x.Gram, y.Gram, x.RecSHV, y.RecSHV, os);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> eBigMat = *opt;
    data.onl_gens.insert_generator(eBigMat);
#ifdef DEBUG_PERFECT_TSPACE
    os << "PERF_TSPACE: A, |onl_gens|=" << data.onl_gens.size() << "\n";
#endif
    TadjO ret{y.eInc, eBigMat};
    return ret;
  }

  std::pair<Tobj, TadjO> f_spann(TadjI const &y) {
    Tobj x_ret{y.Gram, y.RecSHV, {}};
#ifdef DEBUG_PERFECT_TSPACE
    check_correctness(x_ret);
#endif
    std::pair<MyMatrix<Tint>, MyMatrix<T>> pair = get_reduced_gram(data.onl_gens, y.Gram);
    MyMatrix<Tint> Pinv = Inverse(pair.first);
    Tshortest<T,Tint> RecSHV = apply_transformation(y.RecSHV, Pinv);
    Tobj x_reduced{pair.second, RecSHV, {}};
#ifdef DEBUG_PERFECT_TSPACE
    check_correctness(x_reduced);
#endif
    MyMatrix<Tint> eBigMat = Pinv;
    TadjO ret{y.eInc, eBigMat};
    return {x_reduced, ret};
  }

  std::vector<TadjI> f_adj(Tobj &x) {
    std::ostream &os = get_os();
    std::pair<Tgroup, std::vector<MyMatrix<Tint>>> pair = SimplePerfect_Stabilizer<T, Tint, Tgroup>(data.LinSpa, x.Gram, x.RecSHV, os);
    x.GRP = pair.first;
#ifdef DEBUG_PERFECT_TSPACE
    os << "PERF_TSPACE: After x.GRP set\n";
#endif
    for (auto & eGenMat: pair.second) {
      data.onl_gens.insert_generator(eGenMat);
    }
#ifdef DEBUG_PERFECT_TSPACE
    os << "PERF_TSPACE: B, |onl_gens|=" << data.onl_gens.size() << "\n";
#endif
    return TSPACE_GetAdjacencies<T, Tint>(data.LinSpa, x.Gram, x.RecSHV, x.GRP, os);
  }

  Tobj f_adji_obj(TadjI const &x) {
    Tobj x_ret{x.Gram, x.RecSHV, {}};
    return x_ret;
  }
};

template <typename T, typename Tint, typename Tgroup>
void WriteEntryGAP(std::ostream &os, PerfectTspace_Obj<T, Tint, Tgroup> const &obj) {
  os << "rec(Gram:=";
  WriteMatrixGAP(os, obj.Gram);
  os << ")";
}

template <typename T, typename Tint, typename Tgroup>
void WriteEntryPYTHON(std::ostream &os, PerfectTspace_Obj<T, Tint, Tgroup> const &obj) {
  os << "{\"Gram\":" << StringMatrixPYTHON(obj.Gram) << "}";
}

template <typename Tint>
void WriteEntryGAP(std::ostream &os, PerfectTspace_AdjO<Tint> const &adj) {
  os << "rec(eInc:=";
  WriteEntryGAP(os, adj.eInc);
  os << ", eBigMat:=";
  WriteMatrixGAP(os, adj.eBigMat);
  os << ")";
}

template <typename Tint>
void WriteEntryPYTHON(std::ostream &os, PerfectTspace_AdjO<Tint> const &adj) {
  os << "{\"eInc\":";
  WriteEntryPYTHON(os, adj.eInc);
  os << ", \"eBigMat\":" << StringMatrixPYTHON(adj.eBigMat) << "}";
}

template <typename T, typename Tint, typename Tgroup>
void WriteDetailedEntryGAP(std::ostream &os_out,
                           DataPerfectTspace<T, Tint, Tgroup> const& data,
                           PerfectTspace_Obj<T, Tint, Tgroup> const &obj, 
                           [[maybe_unused]] std::ostream& os) {
  os_out << "rec(";
  WriteEntryGAP(os_out, obj);
  os_out << ", GRPsize:=" << obj.GRP.size();
  os_out << ", n:=" << data.LinSpa.n;
  os_out << ")";
}


namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, PerfectTspace_Obj<T, Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("Gram", eRec.Gram);
  ar &make_nvp("GRP", eRec.GRP);
}
template <class Archive, typename T, typename Tint>
inline void serialize(Archive &ar, PerfectTspace_AdjI<T, Tint> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("eInc", eRec.eInc);
  ar &make_nvp("Gram", eRec.Gram);
}
template <class Archive, typename Tint>
inline void serialize(Archive &ar, PerfectTspace_AdjO<Tint> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("eInc", eRec.eInc);
  ar &make_nvp("eBigMat", eRec.eBigMat);
}
} // namespace boost::serialization

FullNamelist NAMELIST_GetStandard_ENUMERATE_PERFECT_TSPACE() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListStringValues1["arithmetic_T"] = "gmp_rational";
  ListStringValues1["arithmetic_Tint"] = "gmp_integer";
  ListStringValues1["OutFormat"] = "nothing";
  ListStringValues1["OutFile"] = "unset.out";
  ListStringValues1["FileDualDescription"] = "unset";
  ListStringValues1["Prefix"] = "unset";
  ListIntValues1["max_runtime_second"] = 0;
  ListBoolValues1["ApplyStdUnitbuf"] = false;
  SingleBlock BlockDATA;
  BlockDATA.setListIntValues(ListIntValues1);
  BlockDATA.setListBoolValues(ListBoolValues1);
  BlockDATA.setListStringValues(ListStringValues1);
  ListBlock["DATA"] = BlockDATA;
  // STORAGE
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, std::string> ListStringValues2;
  ListBoolValues2["Saving"] = false;
  ListStringValues2["Prefix"] = "/irrelevant/";
  SingleBlock BlockSTORAGE;
  BlockSTORAGE.setListBoolValues(ListBoolValues2);
  BlockSTORAGE.setListStringValues(ListStringValues2);
  ListBlock["STORAGE"] = BlockSTORAGE;
  // TSPACE
  ListBlock["TSPACE"] = SINGLEBLOCK_Get_Tspace_Description();
  // Merging all data
  return FullNamelist(ListBlock);
}

// clang-format off
#endif  // SRC_PERFECT_PERFECT_TSPACE_H_
// clang-format on

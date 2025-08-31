// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_PERFECT_TSPACE_H_
#define SRC_PERFECT_PERFECT_TSPACE_H_

// clang-format off
#include "LatticeStabEquiCan.h"
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
  int n;
  LinSpaceMatrix<T> LinSpa;
  RecordDualDescOperation<T, Tgroup> rddo;
  std::ostream &get_os() { return rddo.os; }
};


template <typename T, typename Tint, typename Tgroup>
std::vector<PerfectTspace_AdjI<T, Tint>>
TSPACE_GetAdjacencies(LinSpaceMatrix<T> const &LinSpa,
                      MyMatrix<T> const &eGram,
                      Tshortest<T, Tint> const &RecSHV,
                      Tgroup const& GRP,
                      std::ostream &os) {
#ifdef DEBUG_PERFECT_TSPACE
  os << "PERF_TSPACE: TSPACE_GetAdjacencies, begin\n";
#endif
  MyMatrix<T> SHV_T = conversion_and_duplication<T,Tint>(RecSHV.SHV);

  RyshkovGRP<T,Tgroup> eCone = GetNakedPerfectCone_GRP<T,Tgroup>(LinSpa,
                                                                 eGram,
                                                                 SHV_T,
                                                                 GRP, os);
  vectface ListIncd = DualDescriptionStandard<T,Tgroup>(eCone.PerfDomEXT, eCone.GRPsub);
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
    std::pair<MyMatrix<T>, Tshortest<T, Tint>> pair =
        Flipping_Perfect<T, Tint>(eGram, eMatDir, os);

    PerfectTspace_AdjI<T, Tint> eAdj{eIncd, pair.first, pair.second};
    ListAdj.push_back(eAdj);
  }
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
    return ComputeInvariantPerfectTspace<T, Tint>(seed, x.Gram, x.RecSHV, data.rddo.os);
  }

  std::optional<TadjO> f_repr(Tobj const &x, TadjI const &y) {
    std::ostream &os = get_os();
    std::optional<MyMatrix<Tint>> opt = SimplePerfect_TestEquivalence<T, Tint, Tgroup>(data.LinSpa, x.Gram, y.Gram, x.RecSHV, x.RecSHV, os);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> eBigMat = *opt;
    TadjO ret{y.eInc, eBigMat};
    return ret;
  }

  std::pair<Tobj, TadjO> f_spann(TadjI const &y) {
    Tobj x_ret{y.Gram, y.RecSHV, {}};
    MyMatrix<Tint> eBigMat = IdentityMat<Tint>(data.n);
    TadjO ret{y.eInc, eBigMat};
    return {x_ret, ret};
  }

  std::vector<TadjI> f_adj(Tobj &x) {
    std::ostream &os = get_os();
    x.GRP = SimplePerfect_Stabilizer<T, Tint, Tgroup>(data.LinSpa, x.Gram, x.RecSHV, os);
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
  os_out << ", n:=" << data.n;
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
  std::map<std::string, int> ListIntValues3;
  std::map<std::string, std::string> ListStringValues3;
  ListStringValues3["TypeTspace"] = "File";
  ListStringValues3["FileLinSpa"] = "unset.linspa";
  ListStringValues3["SuperMatMethod"] = "NotNeeded";
  ListStringValues3["ListComm"] = "Trivial";
  ListStringValues3["PtGroupMethod"] = "Trivial";
  ListStringValues3["FileListSubspaces"] = "unset";
  ListIntValues3["RealImagDim"] = 0;
  ListIntValues3["RealImagSum"] = 0;
  ListIntValues3["RealImagProd"] = 0;
  SingleBlock BlockTSPACE;
  BlockTSPACE.setListIntValues(ListIntValues3);
  BlockTSPACE.setListStringValues(ListStringValues3);
  ListBlock["TSPACE"] = BlockTSPACE;
  // Merging all data
  return FullNamelist(ListBlock);
}

// clang-format off
#endif  // SRC_PERFECT_PERFECT_TSPACE_H_
// clang-format on

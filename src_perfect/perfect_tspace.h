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

template <typename T, typename Tint>
size_t ComputeInvariantPerfectTspace(size_t const &seed,
                                     MyMatrix<T> const &eGram,
                                     Tshortest<T, Tint> const &RecSHV,
                                     std::ostream &os) {
#ifdef DEBUG_PERFECT_TSPACE
  os << "PERF_TSPACE: ComputeInvariantPerfectTspace, begin\n";
#endif
  using Tidx_value = int16_t;
  int n = eGram.rows();
  MyMatrix<T> eG = eGram / RecSHV.eMin;
  int nbSHV = RecSHV.SHV.rows();

  MyVector<T> V(n);
  auto f1 = [&](size_t iRow) -> void {
    for (int i=0; i<n; i++) {
      T eSum(0);
      for (int j=0; j<n; j++)
        eSum += eG(i, j) * RecSHV.SHV(iRow, i);
      V(i) = eSum;
    }
  };
  auto f2 = [&](size_t jRow) -> T {
    T eSum(0);
    for (int i=0; i<n; i++)
      eSum += V(i) * RecSHV.SHV(jRow, i);
    return eSum;
  };
  WeightMatrix<true, T, Tidx_value> WMat(nbSHV, f1, f2, os);
#ifdef DEBUG_PERFECT_TSPACE
  os << "PERF_TSPACE: ComputeInvariantPerfectTspace, We have WMat\n";
#endif
  size_t hash = GetInvariantWeightMatrix(seed, WMat);
#ifdef DEBUG_PERFECT_TSPACE
  os << "PERF_TSPACE: ComputeInvariantPerfectTspace, We have hash\n";
#endif
  return hash;
}

template <typename T, typename Tint>
std::vector<PerfectTspace_AdjI<T, Tint>>
TSPACE_GetAdjacencies(LinSpaceMatrix<T> const &LinSpa, MyMatrix<T> const &Gram,
                     std::ostream &os) {
#ifdef DEBUG_PERFECT_TSPACE
  os << "PERF_TSPACE: TSPACE_GetAdjacencies, begin\n";
#endif
  Tshortest<T, Tint> eRec = T_ShortestVector<T, Tint>(Gram, os);
  int n = eRec.SHV.cols();
  int nbShort = eRec.SHV.rows() / 2;
  MyMatrix<Tint> SHVred(nbShort, n);
  for (int iShort = 0; iShort < nbShort; iShort++)
    for (int i = 0; i < n; i++)
      SHVred(iShort, i) = eRec.SHV(2 * iShort, i);

  MyMatrix<T> ConeClassical = GetNakedPerfectConeClassical<T, Tint>(SHVred);
  vectface ListIncd = lrs::DualDescription_incd(ConeClassical);

  std::vector<PerfectTspace_AdjI<T, Tint>> ListAdj;
  for (auto &eIncd : ListIncd) {
    MyVector<T> eFacet = FindFacetInequality(ConeClassical, eIncd);
    MyMatrix<T> eMatDir = LINSPA_GetMatrixInTspace(LinSpa, eFacet);
    std::pair<MyMatrix<T>, Tshortest<T, Tint>> pair =
        Flipping_Perfect<T, Tint>(Gram, eMatDir, os);

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
    return TSPACE_GetAdjacencies<T, Tint>(data.LinSpa, x.Gram, os);
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

// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_PERFECT_TSPACE_H_
#define SRC_PERFECT_PERFECT_TSPACE_H_

// clang-format off
#include "LatticeStabEquiCan.h"
#include "Temp_PerfectForm.h"
#include "POLY_lrslib.h"
#include "Positivity.h"
#include "POLY_RecursiveDualDesc.h"
#include "hash_functions.h"
// clang-format on

template <typename T, typename Tint> struct PerfectTspaceEntry {
  MyMatrix<Tint> Gram;
  int incd;
};

template <typename T, typename Tint, typename Tgroup>
struct PerfectTspace_Obj {
  MyMatrix<Tint> Gram;
  Tgroup GRP;
};

template <typename Tint> struct PerfectTspace_AdjI {
  Face eInc;
  MyMatrix<Tint> Gram;
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
size_t ComputeInvariantPerfectTspace(size_t const &seed, MyMatrix<Tint> const &Gram,
                                    std::ostream &os) {
  std::ostringstream ostringstream;
  ostringstream << "Perfect_Tspace Gram=" << StringMatrixGAP(Gram);
  std::string converted(ostringstream.str());
  size_t h = std::hash<std::string>()(converted);
  return robin_hood_hash_bytes(&h, sizeof(h), seed);
}

template <typename T, typename Tint>
PerfectTspaceEntry<T, Tint> TSPACE_GetOnePerfect(LinSpaceMatrix<T> const &LinSpa,
                                                std::ostream &os) {
  MyMatrix<T> eGram = GetOnePerfectForm<T>(LinSpa, os);
  Tshortest<T, Tint> eRec = T_ShortestVector<T, Tint>(eGram, os);
  int incd = eRec.SHV.rows() / 2;
  MyMatrix<Tint> Gram_int = UniversalMatrixConversion<Tint, T>(eGram);
  return {Gram_int, incd};
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
TSPACE_TestEquivalence(LinSpaceMatrix<T> const &LinSpa,
                      MyMatrix<Tint> const &Gram1,
                      MyMatrix<Tint> const &Gram2,
                      std::ostream &os) {
  MyMatrix<T> Gram1_T = UniversalMatrixConversion<T, Tint>(Gram1);
  MyMatrix<T> Gram2_T = UniversalMatrixConversion<T, Tint>(Gram2);
  
  Tshortest<T, Tint> eRec1 = T_ShortestVector<T, Tint>(Gram1_T, os);
  Tshortest<T, Tint> eRec2 = T_ShortestVector<T, Tint>(Gram2_T, os);
  
  return PERF_TestEquivalence<T, Tint, Tgroup>(LinSpa, Gram1_T, Gram2_T,
                                              eRec1.SHV, eRec2.SHV);
}

template <typename T, typename Tint, typename Tgroup>
Tgroup TSPACE_ComputeStabilizer(LinSpaceMatrix<T> const &LinSpa,
                               MyMatrix<Tint> const &Gram,
                               std::ostream &os) {
  MyMatrix<T> Gram_T = UniversalMatrixConversion<T, Tint>(Gram);
  Tshortest<T, Tint> eRec = T_ShortestVector<T, Tint>(Gram_T, os);
  return SimplePerfect_Stabilizer<T, Tint, Tgroup>(LinSpa, Gram_T, eRec);
}

template <typename T, typename Tint>
std::vector<PerfectTspace_AdjI<Tint>>
TSPACE_GetAdjacencies(LinSpaceMatrix<T> const &LinSpa, MyMatrix<Tint> const &Gram,
                     std::ostream &os) {
  MyMatrix<T> Gram_T = UniversalMatrixConversion<T, Tint>(Gram);
  Tshortest<T, Tint> eRec = T_ShortestVector<T, Tint>(Gram_T, os);
  int n = eRec.SHV.cols();
  int nbShort = eRec.SHV.rows() / 2;
  MyMatrix<Tint> SHVred(nbShort, n);
  for (int iShort = 0; iShort < nbShort; iShort++)
    for (int i = 0; i < n; i++)
      SHVred(iShort, i) = eRec.SHV(2 * iShort, i);
      
  MyMatrix<T> ConeClassical = GetNakedPerfectConeClassical<T, Tint>(SHVred);
  vectface ListIncd = lrs::DualDescription_incd(ConeClassical);
  
  std::vector<PerfectTspace_AdjI<Tint>> ListAdj;
  for (auto &eIncd : ListIncd) {
    MyVector<T> eFacet = FindFacetInequality(ConeClassical, eIncd);
    MyMatrix<T> eMatDir = LINSPA_GetMatrixInTspace(LinSpa, eFacet);
    std::pair<MyMatrix<T>, Tshortest<T, Tint>> ePairAdj =
        Flipping_Perfect<T, Tint>(Gram_T, eMatDir, os);
    
    MyMatrix<T> eMat2 = ComputeCanonicalForm<T, Tint>(ePairAdj.first, std::cerr).Mat;
    MyMatrix<T> eMat3 = RemoveFractionMatrix(eMat2);
    MyMatrix<Tint> Gram_adj = UniversalMatrixConversion<Tint, T>(eMat3);
    
    PerfectTspace_AdjI<Tint> eAdj{eIncd, Gram_adj};
    ListAdj.push_back(eAdj);
  }
  return ListAdj;
}

template <typename T, typename Tint, typename Tgroup>
struct DataPerfectTspaceFunc {
  DataPerfectTspace<T, Tint, Tgroup> data;
  using Tobj = PerfectTspace_Obj<T, Tint, Tgroup>;
  using TadjI = PerfectTspace_AdjI<Tint>;
  using TadjO = PerfectTspace_AdjO<Tint>;
  std::ostream &get_os() { return data.rddo.os; }
  
  Tobj f_init() {
    std::ostream &os = get_os();
    PerfectTspaceEntry<T, Tint> eRec = TSPACE_GetOnePerfect<T, Tint>(data.LinSpa, os);
    Tobj x{eRec.Gram, {}};
    return x;
  }
  
  size_t f_hash(size_t const &seed, Tobj const &x) {
    return ComputeInvariantPerfectTspace<T, Tint>(seed, x.Gram, data.rddo.os);
  }
  
  std::optional<TadjO> f_repr(Tobj const &x, TadjI const &y) {
    std::ostream &os = get_os();
    std::optional<MyMatrix<Tint>> opt = TSPACE_TestEquivalence<T, Tint, Tgroup>(
        data.LinSpa, x.Gram, y.Gram, os);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> eBigMat = *opt;
    TadjO ret{y.eInc, eBigMat};
    return ret;
  }
  
  std::pair<Tobj, TadjO> f_spann(TadjI const &y) {
    Tobj x_ret{y.Gram, {}};
    MyMatrix<Tint> eBigMat = IdentityMat<Tint>(data.n);
    TadjO ret{y.eInc, eBigMat};
    return {x_ret, ret};
  }
  
  std::vector<TadjI> f_adj(Tobj &x) {
    std::ostream &os = get_os();
    MyMatrix<T> Gram_T = UniversalMatrixConversion<T, Tint>(x.Gram);
    x.GRP = TSPACE_ComputeStabilizer<T, Tint, Tgroup>(data.LinSpa, x.Gram, os);
    return TSPACE_GetAdjacencies<T, Tint>(data.LinSpa, x.Gram, os);
  }
  
  Tobj f_adji_obj(TadjI const &x) {
    Tobj x_ret{x.Gram, {}};
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


namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, PerfectTspace_Obj<T, Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("Gram", eRec.Gram);
  ar &make_nvp("GRP", eRec.GRP);
}
template <class Archive, typename Tint>
inline void serialize(Archive &ar, PerfectTspace_AdjI<Tint> &eRec,
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

// clang-format off
#endif  // SRC_PERFECT_PERFECT_TSPACE_H_
// clang-format on
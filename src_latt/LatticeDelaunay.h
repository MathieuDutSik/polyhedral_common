// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEDELAUNAY_H_
#define SRC_LATT_LATTICEDELAUNAY_H_

// clang-format off
#include "boost_serialization.h"
#include "FundamentalDelaunay.h"
#include "GRP_DoubleCoset.h"
#include "MatrixGroup.h"
#include "POLY_RecursiveDualDesc.h"
#include "POLY_AdjacencyScheme.h"
#include "GraverBasis.h"
#include <map>
#include <string>
#include <vector>
// clang-format on

#ifdef TIMINGS
#define TIMINGS_DELAUNAY_ENUMERATION
#endif

#ifdef DEBUG
#define DEBUG_DELAUNAY_ENUMERATION
#endif

template <typename T, typename Tint, typename Tgroup>
struct DataLattice {
  int n;
  MyMatrix<T> GramMat;
  MyMatrix<T> SHV;
  CVPSolver<T,Tint> solver;
  MyMatrix<Tint> ShvGraverBasis;
  RecordDualDescOperation<T,Tgroup> rddo;
};

template <typename T, typename Tint, typename Tgroup>
DataLattice<T,Tint,Tgroup> GetDataLattice(MyMatrix<T> const& GramMat, std::ostream& os) {
  using TintGroup = typename Tgroup::Tint;
  int n = GramMat.rows();
  int dimEXT = n + 1;
  MyMatrix<T> SHV(0,n);
  CVPSolver<T,Tint> solver(GramMat, os);
  MyMatrix<Tint> ShvGraverBasis = GetGraverBasis<T,Tint>(GramMat);
  PolyHeuristicSerial<TintGroup> AllArr = AllStandardHeuristicSerial<TintGroup>(dimEXT, os);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);
  return {n, GramMat, SHV, solver, ShvGraverBasis, std::move(rddo)};
}



template <typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value>
GetWeightMatrixFromGramEXT(MyMatrix<T> const &EXT, MyMatrix<T> const &GramMat,
                           MyMatrix<T> const &SHV, std::ostream& os) {
  int n = GramMat.rows();
  CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral(GramMat, EXT);
  MyVector<T> TheCenter = eCP.eCent;
  std::vector<MyMatrix<T>> CharPair = CharacterizingPair(GramMat, TheCenter);
  MyMatrix<T> Qmat = CharPair[0];
  int nbVect = SHV.rows();
  int nbVert = EXT.rows();
  MyMatrix<T> EXText(nbVect + nbVert, n + 1);
  for (int iVert = 0; iVert < nbVert; iVert++)
    EXText.row(iVert) = EXT.row(iVert);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    EXText(nbVert + iVect, 0) = 0;
    for (int i = 0; i < n; i++)
      EXText(nbVert + iVect, i + 1) = SHV(iVect, i);
  }
  return GetSimpleWeightMatrix<T, Tidx_value>(EXText, Qmat, os);
}

template <typename T, typename Tgroup>
bool IsGroupCorrect(MyMatrix<T> const &EXT_T, Tgroup const &eGRP) {
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> LGen = eGRP.GeneratorsOfGroup();
  for (auto &eGen : LGen) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    if (!IsIntegralMatrix(eMat))
      return false;
  }
  return true;
}

template <typename T, typename Tint, typename Tgroup>
Tgroup Delaunay_Stabilizer(DataLattice<T, Tint, Tgroup> const &eData,
                           MyMatrix<Tint> const &EXT, std::ostream& os) {
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  MicrosecondTime time;
#endif
  using Tidx_value = int16_t;
  using Tgr = GraphListAdj;
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(EXT);
  //
  // Now extending with the SHV vector set
  //
  WeightMatrix<true, T, Tidx_value> WMat =
    GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT_T, eData.GramMat, eData.SHV, os);
  int nbVert = EXT_T.rows();
  int nbSHV = eData.SHV.rows();
  Face eFace(nbVert + nbSHV);
  for (int iVert = 0; iVert < nbVert; iVert++)
    eFace[iVert] = 1;
  for (int iSHV = 0; iSHV < nbSHV; iSHV++)
    eFace[nbVert + iSHV] = 0;
  Tgroup PreGRPisom =
    GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat, os);
  Tgroup GRPisom = ReducedGroupAction(PreGRPisom, eFace);
  Tgroup GRPlatt = LinPolytopeIntegral_Stabilizer_Method8(EXT_T, GRPisom, os);
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  os << "|Delaunay_Stabilizer|=" << time << "\n";
#endif
  return GRPlatt;
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
Delaunay_TestEquivalence(DataLattice<T, Tint, Tgroup> const &eData,
                         MyMatrix<Tint> const &EXT1,
                         MyMatrix<Tint> const &EXT2,
                         std::ostream & os) {
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  MicrosecondTime time;
#endif
  using Telt = typename Tgroup::Telt;
  using Tgr = GraphListAdj;
  using Tidx_value = int16_t;
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: Begin Delaunay_TestEquivalence\n";
#endif
  MyMatrix<T> EXT1_T = UniversalMatrixConversion<T, Tint>(EXT1);
  MyMatrix<T> EXT2_T = UniversalMatrixConversion<T, Tint>(EXT2);
  //
  // Now extending by adding more vectors.
  //
  WeightMatrix<true, T, Tidx_value> WMat1 =
    GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT1_T, eData.GramMat, eData.SHV, os);
  WeightMatrix<true, T, Tidx_value> WMat2 =
    GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT2_T, eData.GramMat, eData.SHV, os);
  std::optional<Telt> eRes =
    TestEquivalenceWeightMatrix<T, Telt, Tidx_value>(WMat1, WMat2, os);
  if (!eRes) {
#ifdef DEBUG_DELAUNAY_ENUMERATION
    os << "DEL_ENUM: Leaving Delaunay_TestEquivalence 1 with false\n";
#endif
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    os << "|Delaunay_TestEquivalence|=" << time << "\n";
#endif
    return {};
  }
  Telt const& eElt = *eRes;
  MyMatrix<T> MatEquiv_T = FindTransformation(EXT1_T, EXT2_T, eElt);
  if (IsIntegralMatrix(MatEquiv_T)) {
    MyMatrix<Tint> MatEquiv_I = UniversalMatrixConversion<Tint, T>(MatEquiv_T);
#ifdef DEBUG_DELAUNAY_ENUMERATION
    os << "DEL_ENUM: Leaving Delaunay_TestEquivalence 2 with true\n";
#endif
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    os << "|Delaunay_TestEquivalence|=" << time << "\n";
#endif
    return MatEquiv_I;
  }
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: Trying other strategies\n";
#endif
  Tgroup GRP1 = GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat1, os);
  std::optional<MyMatrix<T>> eEq = LinPolytopeIntegral_Isomorphism_Method8(EXT1_T, EXT2_T, GRP1, eElt, os);
  if (!eEq) {
#ifdef DEBUG_DELAUNAY_ENUMERATION
    os << "DEL_ENUM: Leaving Delaunay_TestEquivalence 3 with false\n";
#endif
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    os << "|Delaunay_TestEquivalence|=" << time << "\n";
#endif
    return {};
  }
  MyMatrix<Tint> eMat_I = UniversalMatrixConversion<Tint, T>(*eEq);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: Leaving Delaunay_TestEquivalence 4 with true\n";
#endif
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  os << "|Delaunay_TestEquivalence|=" << time << "\n";
#endif
  return eMat_I;
}

template <typename T, typename Tint, typename Tgroup>
size_t ComputeInvariantDelaunay(DataLattice<T, Tint, Tgroup> const &eData,
                                size_t const& seed,
                                MyMatrix<Tint> const& EXT, [[maybe_unused]] std::ostream & os) {
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  MicrosecondTime time;
#endif
  int nbVert = EXT.rows();
  int n = EXT.cols() - 1;
  Tint PreIndex = Int_IndexLattice(EXT);
  Tint eIndex = T_abs(PreIndex);
  MyMatrix<T> EXT_fullT = UniversalMatrixConversion<T,Tint>(EXT);
  CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral(eData.GramMat, EXT_fullT);
  MyMatrix<T> EXT_T(nbVert, n);
  for (int iVert=0; iVert<nbVert; iVert++) {
    for (int i=0; i<n; i++) {
      T val = UniversalScalarConversion<T,Tint>(EXT(iVert, i+1));
      EXT_T(iVert, i) = val - eCP.eCent(i+1);
    }
  }
  std::map<T, size_t> ListDiagNorm;
  std::map<T, size_t> ListOffDiagNorm;
  MyVector<T> V(n);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    for (int i=0; i<n; i++) {
      T eSum = 0;
      for (int j=0; j<n; j++) {
        eSum += eData.GramMat(i,j) * EXT_T(iVert, j);
      }
      V(i) = eSum;
    }
    T scal = 0;
    for (int i=0; i<n; i++) {
      scal += V(i) * EXT_T(iVert, i);
    }
    ListDiagNorm[scal] += 1;
    for (int jVert=iVert+1; jVert<nbVert; jVert++) {
      T scal = 0;
      for (int i=0; i<n; i++) {
        scal += V(i) * EXT_T(jVert, i);
      }
      ListOffDiagNorm[scal] += 1;
    }
  }
  size_t hash = ComputeHashTwoMap(seed, ListDiagNorm, ListOffDiagNorm);
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  os << "|ComputeInvariantDelaunay|=" << time << "\n";
#endif
  return hash;
}

template<typename Tint>
struct Delaunay_AdjI {
  Face eInc;
  MyMatrix<Tint> EXT;
};

namespace boost::serialization {
  template <class Archive, typename Tint>
  inline void serialize(Archive &ar, Delaunay_AdjI<Tint> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("eInc", eRec.eInc);
    ar &make_nvp("EXT", eRec.EXT);
  }
}

template<typename Tint>
struct Delaunay_AdjO {
  Face eInc;
  MyMatrix<Tint> eBigMat;
  int iOrb;
};

namespace boost::serialization {
  template <class Archive, typename Tint>
  inline void serialize(Archive &ar, Delaunay_AdjO<Tint> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("eInc", eRec.eInc);
    ar &make_nvp("eBigMat", eRec.eBigMat);
    ar &make_nvp("iOrb", eRec.iOrb);
  }
}

template<typename Tint>
struct Delaunay_AdjO_spec {
  Face eInc;
  MyMatrix<Tint> eBigMat;
};

namespace boost::serialization {
  template <class Archive, typename Tint>
  inline void serialize(Archive &ar, Delaunay_AdjO_spec<Tint> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("eInc", eRec.eInc);
    ar &make_nvp("eBigMat", eRec.eBigMat);
  }
}

template<typename Tint, typename Tgroup>
struct Delaunay_Entry {
  MyMatrix<Tint> EXT;
  Tgroup GRP;
  std::vector<Delaunay_AdjO<Tint>> ListAdj;
};

namespace boost::serialization {
  template <class Archive, typename Tint, typename Tgroup>
  inline void serialize(Archive &ar, Delaunay_Entry<Tint,Tgroup> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("EXT", eRec.EXT);
    ar &make_nvp("GRP", eRec.GRP);
    ar &make_nvp("ListAdj", eRec.ListAdj);
  }
}

template<typename Tvert, typename Tgroup>
struct DelaunayTesselation {
  std::vector<Delaunay_Entry<Tvert,Tgroup>> l_dels;
};

namespace boost::serialization {
  template <class Archive, typename Tvert, typename Tgroup>
  inline void serialize(Archive &ar, DelaunayTesselation<Tvert,Tgroup> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("l_dels", eRec.l_dels);
  }
}

template<typename Tvert, typename Tgroup>
void check_delaunay_tessellation(DelaunayTesselation<Tvert,Tgroup> const& DT, [[maybe_unused]] std::ostream& os) {
  for (auto & eDel : DT.l_dels) {
    MyMatrix<Tvert> const& EXT = eDel.EXT;
    ContainerMatrix<Tvert> cont(EXT);
    for (auto & eAdj : eDel.ListAdj) {
      Face const& eInc = eAdj.eInc;
      Face eIncEff(EXT.rows());
      MyMatrix<Tvert> const& EXT2 = DT.l_dels[eAdj.iOrb].EXT;
      MyMatrix<Tvert> const& eBigMat = eAdj.eBigMat;
#ifdef DEBUG_DELAUNAY_ENUMERATION
      os << "DEL_ENUM: check_delaunay_tessellation |EXT2}=" << EXT2.rows() << "/" << EXT2.cols() << " |eBigMat|=" << eBigMat.rows() << "/" << eBigMat.cols() << "\n";
#endif
      MyMatrix<Tvert> EXTadj = EXT2 * eBigMat;
      int len = EXTadj.rows();
      for (int u=0; u<len; u++) {
        MyVector<Tvert> V = GetMatrixRow(EXTadj, u);
        std::optional<size_t> opt = cont.GetIdx_v(V);
        if (opt) {
          eIncEff[*opt] = 1;
        }
      }
      if (eIncEff != eInc) {
        std::cerr << "Inconsistency in the adjacency\n";
        throw TerminalException{1};
      }
    }
  }
}


template<typename Tvert, typename Tgroup>
void WriteEntryGAP(std::ostream& os_out, DelaunayTesselation<Tvert,Tgroup> const& DT) {
  using Telt = typename Tgroup::Telt;
  os_out << "[";
  size_t n_del = DT.l_dels.size();
  for (size_t i_del=0; i_del<n_del; i_del++) {
    Delaunay_Entry<Tvert,Tgroup> const& eDel = DT.l_dels[i_del];
    MyMatrix<Tvert> const& EXT = eDel.EXT;
    if (i_del > 0)
      os_out << ",";
    os_out << "rec(EXT:=" << StringMatrixGAP(EXT) << ",\n";
    std::vector<Telt> LGen = eDel.GRP.SmallGeneratingSet();
    auto get_gap_string=[&]() -> std::string {
      if (LGen.size() == 0) {
        return "Group(())";
      } else {
        std::string str_ret = "Group([";
        bool IsFirst = true;
        for (auto & eElt : LGen) {
          if (!IsFirst)
            str_ret += ",";
          IsFirst = false;
          str_ret += GapStyleString(eElt);
        }
        str_ret += "])";
        return str_ret;
      }
    };
    std::string str_perm = "[", str_matr = "[";
    bool IsFirst = true;
    for (auto & eElt : LGen) {
      if (!IsFirst) {
        str_perm += ",";
        str_matr += ",";
      }
      IsFirst = false;
      MyMatrix<Tvert> M = RepresentVertexPermutation(EXT, EXT, eElt);
      str_perm += GapStyleString(eElt);
      str_matr += StringMatrixGAP(M);
    }
    str_perm += "]";
    str_matr += "]";
    std::string str_phi = "GroupHomomorphismByImagesNC(Group(" + str_perm + "), Group(" + str_matr + "), " + str_perm + ", " + str_matr + ")";
    os_out << "TheStab:=rec(PermutationStabilizer:=" << get_gap_string() << ", PhiPermMat:=" << str_phi << "), ";
    os_out << "Adjacencies:=[";
    IsFirst = true;
    for (auto & eAdj : eDel.ListAdj) {
      if (!IsFirst)
        os_out << ",";
      IsFirst = false;
      os_out << "rec(iDelaunay:=" << (eAdj.iOrb + 1) << ", ";
      os_out << "eInc:=[";
      std::vector<int> V = FaceToVector<int>(eAdj.eInc);
      for (size_t u=0; u<V.size(); u++) {
        if (u > 0)
          os_out << ",";
        os_out << (V[u] + 1);
      }
      os_out << "],\n";
      os_out << "eBigMat:=" << StringMatrixGAP(eAdj.eBigMat) << ")";
    }
    os_out << "])";
  }
  os_out << "]";
}

template<typename Tvert, typename Tgroup>
void WriteGAPformat(DelaunayTesselation<Tvert,Tgroup> const& DT, std::string const& OutFile) {
  //  using T = typename overlying_field<Tvert>::field_type;
  std::ofstream OUTfs(OutFile);
  OUTfs << "return ";
  WriteEntryGAP(OUTfs, DT);
  OUTfs << ";\n";
}

template<typename T, typename Tint, typename Tgroup>
std::pair<Tgroup, std::vector<Delaunay_AdjI<Tint>>> ComputeGroupAndAdjacencies(DataLattice<T, Tint, Tgroup> & eData, MyMatrix<Tint> const& x) {
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T,Tint>(x);
  std::ostream& os = eData.rddo.os;
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: |EXT_T|=" << EXT_T.rows() << " / " << EXT_T.cols() << "\n";
#endif
  Tgroup GRPlatt = Delaunay_Stabilizer<T, Tint, Tgroup>(eData, x, os);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: |GRPlatt|=" << GRPlatt.size() << "\n";
#endif
  vectface TheOutput = DualDescriptionRecord(EXT_T, GRPlatt, eData.rddo);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: |TheOutput|=" << TheOutput.size() << "\n";
#endif
  std::vector<Delaunay_AdjI<Tint>> ListAdj;
  for (auto &eOrbB : TheOutput) {
    MyMatrix<Tint> EXTadj = FindAdjacentDelaunayPolytope<T, Tint>(eData.GramMat, eData.solver, eData.ShvGraverBasis, EXT_T, eOrbB, os);
    Delaunay_AdjI<Tint> eAdj{eOrbB, EXTadj};
    ListAdj.push_back(eAdj);
  }
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: |ListAdj|=" << ListAdj.size() << "\n";
#endif
  return {GRPlatt, std::move(ListAdj)};
}

template<typename Tint, typename Tgroup>
struct Delaunay_Obj {
  MyMatrix<Tint> EXT;
  Tgroup GRP;
};

namespace boost::serialization {
  template <class Archive, typename Tint, typename Tgroup>
  inline void serialize(Archive &ar, Delaunay_Obj<Tint, Tgroup> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("EXT", eRec.EXT);
    ar &make_nvp("GRP", eRec.GRP);
  }
}


template <typename T, typename Tint, typename Tgroup>
struct DataLatticeFunc {
  DataLattice<T, Tint, Tgroup> data;
  using Tobj = Delaunay_Obj<Tint, Tgroup>;
  using TadjI = Delaunay_AdjI<Tint>;
  using TadjO = Delaunay_AdjO_spec<Tint>;
  std::ostream& get_os() {
    return data.rddo.os;
  }
  Tobj f_init() {
    MyMatrix<Tint> EXT = FindDelaunayPolytope<T, Tint>(data.GramMat, data.solver, data.rddo.os);
    Tobj x{std::move(EXT), {} };
    return x;
  }
  size_t f_hash(size_t const& seed, Tobj const& x) {
    return ComputeInvariantDelaunay(data, seed, x.EXT, data.rddo.os);
  }
  std::optional<TadjO> f_repr(Tobj const& x, TadjI const& y) {
    std::optional<MyMatrix<Tint>> opt = Delaunay_TestEquivalence<T, Tint, Tgroup>(data, x.EXT, y.EXT, data.rddo.os);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const& eBigMat = *opt;
    TadjO ret{y.eInc, eBigMat};
    return ret;
  }
  std::pair<Tobj,TadjO> f_spann(TadjI const& x) {
    MyMatrix<Tint> EXT = x.EXT;
    Tobj x_ret{EXT, {} };
    MyMatrix<Tint> eBigMat = IdentityMat<Tint>(data.n + 1);
    TadjO ret{x.eInc, eBigMat};
    return {x_ret, ret};
  }
  std::vector<TadjI> f_adj(Tobj & x) {
    std::pair<Tgroup, std::vector<TadjI>> pair = ComputeGroupAndAdjacencies<T,Tint,Tgroup>(data, x.EXT);
    x.GRP = pair.first;
    return pair.second;
  }
  Tobj f_adji_obj(TadjI const& x) {
    MyMatrix<Tint> EXT = x.EXT;
    Tobj x_ret{EXT, {} };
    return x_ret;
  };
  size_t f_complexity(Tobj const& x) {
    return x.EXT.rows();
  }
};

template<typename T, typename Tvert, typename Tgroup>
DelaunayTesselation<Tvert, Tgroup> DelaunayTesselation_From_DatabaseEntries_MPI(std::vector<DatabaseEntry_Serial<typename DataLatticeFunc<T, Tvert, Tgroup>::Tobj, typename DataLatticeFunc<T, Tvert, Tgroup>::TadjO>> const& l_ent) {
  std::vector<Delaunay_Entry<Tvert,Tgroup>> l_dels;
  for (auto & eDel : l_ent) {
    std::vector<Delaunay_AdjO<Tvert>> ListAdj;
    for (auto & eAdj : eDel.ListAdj) {
      Delaunay_AdjO<Tvert> fAdj{eAdj.x.eInc, eAdj.x.eBigMat, eAdj.iOrb};
      ListAdj.push_back(fAdj);
    }
    Delaunay_Entry<Tvert,Tgroup> fDel{eDel.x.EXT, eDel.x.GRP, ListAdj};
    l_dels.push_back(fDel);
  }
  return {l_dels};
}








template <typename T, typename Tint, typename Tgroup, typename Fincorrect>
std::optional<DelaunayTesselation<Tint,Tgroup>> EnumerationDelaunayPolytopes(DataLattice<T, Tint, Tgroup> & data,
                                                                             Fincorrect f_incorrect,
                                                                             int const& max_runtime_second) {
#ifdef DEBUG_DELAUNAY_ENUMERATION
  std::ostream& os = data.rddo.os;
  os << "DEL_ENUM: EnumerationDelaunayPolytopes, begin\n";
#endif
  using Tdata = DataLatticeFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = std::vector<DatabaseEntry_Serial<Tobj, TadjO>>;
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: EnumerationDelaunayPolytopes, before EnumerateAndStore_Serial\n";
#endif
  std::optional<Tout> opt = EnumerateAndStore_Serial<Tdata>(data_func, f_incorrect, max_runtime_second);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: EnumerationDelaunayPolytopes, after EnumerateAndStore_Serial\n";
#endif
  if (opt) {
    DelaunayTesselation<Tint, Tgroup> DT = DelaunayTesselation_From_DatabaseEntries_MPI<T,Tint,Tgroup>(*opt);
    return DT;
  } else {
    return {};
  }
}


FullNamelist NAMELIST_GetStandard_COMPUTE_DELAUNAY() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["arithmetic_T"] = "gmp_rational";
  ListStringValues1["arithmetic_Tint"] = "gmp_integer";
  ListStringValues1["GRAMfile"] = "unset.gram";
  ListStringValues1["SVRfile"] = "unset.svr";
  ListStringValues1["OutFormat"] = "nothing";
  ListStringValues1["OutFile"] = "unset.out";
  ListStringValues1["FileDualDescription"] = "unset";
  ListIntValues1["max_runtime_second"] = 0;
  ListBoolValues1["ApplyStdUnitbuf"] = false;
  SingleBlock BlockDATA;
  BlockDATA.ListIntValues = ListIntValues1;
  BlockDATA.ListBoolValues = ListBoolValues1;
  BlockDATA.ListDoubleValues = ListDoubleValues1;
  BlockDATA.ListStringValues = ListStringValues1;
  BlockDATA.ListListStringValues = ListListStringValues1;
  ListBlock["DATA"] = BlockDATA;
  // STORAGE
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListBoolValues2["Saving"] = false;
  ListStringValues2["Prefix"] = "/irrelevant/";
  SingleBlock BlockSTORAGE;
  BlockSTORAGE.ListIntValues = ListIntValues2;
  BlockSTORAGE.ListBoolValues = ListBoolValues2;
  BlockSTORAGE.ListDoubleValues = ListDoubleValues2;
  BlockSTORAGE.ListStringValues = ListStringValues2;
  BlockSTORAGE.ListListStringValues = ListListStringValues2;
  ListBlock["STORAGE"] = BlockSTORAGE;
  // Merging all data
  return {ListBlock, "undefined"};
}

template<typename T, typename Tvert, typename Tgroup>
void WriteFamilyDelaunay(boost::mpi::communicator &comm, std::string const& OutFormat, std::string const& OutFile, std::vector<DatabaseEntry_MPI<typename DataLatticeFunc<T, Tvert, Tgroup>::Tobj, typename DataLatticeFunc<T, Tvert, Tgroup>::TadjO>> const& ListDel, std::ostream & os) {
  int i_rank = comm.rank();
  if (OutFormat == "nothing") {
    std::cerr << "No output\n";
    return;
  }
  using Tout = DatabaseEntry_Serial<typename DataLatticeFunc<T, Tvert, Tgroup>::Tobj, typename DataLatticeFunc<T, Tvert, Tgroup>::TadjO>;
  if (OutFormat == "CheckMergedOutput") {
    int i_proc_out = 0;
    std::vector<Tout> l_ent = my_mpi_gather(comm, ListDel, i_proc_out);
    if (i_proc_out == i_rank) {
      DelaunayTesselation<Tvert, Tgroup> DT = DelaunayTesselation_From_DatabaseEntries_MPI<T,Tvert,Tgroup>(l_ent);
      check_delaunay_tessellation(DT, os);
    }
    std::cerr << "The Delaunay tesselation passed the adjacency check\n";
    return;
  }
  if (OutFormat == "GAP") {
    int i_proc_out = 0;
    std::vector<Tout> l_ent = my_mpi_gather(comm, ListDel, i_proc_out);
    if (i_proc_out == i_rank) {
      DelaunayTesselation<Tvert, Tgroup> DT = DelaunayTesselation_From_DatabaseEntries_MPI<T,Tvert,Tgroup>(l_ent);
      WriteGAPformat(DT, OutFile);
    }
    std::cerr << "The Delaunay tesselation has been written to file\n";
    return;
  }
  if (OutFormat == "RAW") {
    std::ofstream OUTfs(OutFile);
    int nbDel = ListDel.size();
    OUTfs << "nbDel=" << nbDel << "\n";
    for (int iDel = 0; iDel < nbDel; iDel++) {
      OUTfs << "iDel=" << iDel << "/" << nbDel << "\n";
      WriteMatrix(OUTfs, ListDel[iDel].x.EXT);
    }
  }
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
  throw TerminalException{1};
}


template<typename T, typename Tint, typename Tgroup>
void ComputeDelaunayPolytope(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  std::unique_ptr<std::ofstream> os_ptr = get_mpi_log_stream(comm, eFull);
  std::ostream& os = *os_ptr;
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  SingleBlock BlockSTORAGE = eFull.ListBlock.at("STORAGE");
  //
  bool STORAGE_Saving = BlockSTORAGE.ListBoolValues.at("Saving");
  std::string STORAGE_Prefix = BlockSTORAGE.ListStringValues.at("Prefix");
  CreateDirectory(STORAGE_Prefix);
  //
  int max_runtime_second = BlockDATA.ListIntValues.at("max_runtime_second");
  std::cerr << "max_runtime_second=" << max_runtime_second << "\n";
  std::string GRAMfile = BlockDATA.ListStringValues.at("GRAMfile");
  MyMatrix<T> GramMat = ReadMatrixFile<T>(GRAMfile);
  //
  std::string SVRfile = BlockDATA.ListStringValues.at("SVRfile");
  auto get_SVR=[&]() -> MyMatrix<T> {
    if (IsExistingFile(SVRfile)) {
      return ReadMatrixFile<T>(SVRfile);
    }
    int n = GramMat.rows();
    return ZeroMatrix<T>(0, n);
  };
  MyMatrix<T> SVR = get_SVR();
  std::cerr << "|SVR|=" << SVR.rows() << "\n";
  //
  std::string OutFormat = BlockDATA.ListStringValues.at("OutFormat");
  std::string OutFile = BlockDATA.ListStringValues.at("OutFile");
  std::cerr << "OutFormat=" << OutFormat << " OutFile=" << OutFile << "\n";
  //
  int n = GramMat.rows();
  int dimEXT = n + 1;
  using TintGroup = typename Tgroup::Tint;
  std::string FileDualDesc = BlockDATA.ListStringValues.at("FileDualDescription");
  PolyHeuristicSerial<TintGroup> AllArr =
    Read_AllStandardHeuristicSerial_File<TintGroup>(FileDualDesc, dimEXT, os);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);
  //
  CVPSolver<T,Tint> solver(GramMat, os);
  MyMatrix<Tint> ShvGraverBasis = GetGraverBasis<T,Tint>(GramMat);
  DataLattice<T, Tint, Tgroup> data{n,
                                    GramMat,
                                    SVR,
                                    solver,
                                    ShvGraverBasis,
                                    std::move(rddo)};
  using Tdata = DataLatticeFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_MPI<Tobj, TadjO>;
  //
  std::pair<bool, std::vector<Tout>> pair = EnumerateAndStore_MPI<Tdata>(comm, data_func, STORAGE_Prefix, STORAGE_Saving, max_runtime_second);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: We now have IsFinished=" << pair.first << "\n";
  os << "DEL_ENUM: We now have ListDel |ListDel|=" << pair.second.size() << "\n";
#endif
  //
  if (pair.first) {
    WriteFamilyDelaunay<T, Tint, Tgroup>(comm, OutFormat, OutFile, pair.second, os);
  }
}




// clang-format off
#endif  // SRC_LATT_LATTICEDELAUNAY_H_
// clang-format on

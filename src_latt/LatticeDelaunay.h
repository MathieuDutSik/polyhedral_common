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
#include "SystemNamelist.h"
#include <map>
#include <string>
#include <vector>
#include <memory>
#include <utility>
// clang-format on

#ifdef TIMINGS
#define TIMINGS_DELAUNAY_ENUMERATION
#endif

#ifdef DEBUG
#define DEBUG_DELAUNAY_ENUMERATION
#endif

#ifdef DISABLE_DEBUG_DELAUNAY_ENUMERATION
#undef DEBUG_DELAUNAY_ENUMERATION
#endif

template <typename T, typename Tint, typename Tgroup> struct DataLattice {
  int n;
  MyMatrix<T> GramMat;
  MyMatrix<T> SHV;
  CVPSolver<T, Tint> solver;
  MyMatrix<Tint> ShvGraverBasis;
  std::string choice_initial;
  RecordDualDescOperation<T, Tgroup> rddo;
};

template <typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value>
GetWeightMatrixFromGramEXT(MyMatrix<T> const &EXT, MyMatrix<T> const &GramMat,
                           MyMatrix<T> const &SHV, std::ostream &os) {
  int n = GramMat.rows();
  CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral(GramMat, EXT);
  MyVector<T> TheCenter = eCP.eCent;
  std::vector<MyMatrix<T>> CharPair =
      CharacterizingPair(GramMat, TheCenter, os);
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
    MyMatrix<T> eMat = FindTransformation<T, Telt>(EXT_T, EXT_T, eGen);
    if (!IsIntegralMatrix(eMat))
      return false;
  }
  return true;
}

template <typename T, typename Tint, typename Tgroup>
Tgroup Delaunay_StabilizerKernel(MyMatrix<T> const &GramMat,
                                 MyMatrix<T> const &SHV,
                                 MyMatrix<Tint> const &EXT, std::ostream &os) {
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
      GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT_T, GramMat, SHV, os);
  int nbVert = EXT_T.rows();
  int nbSHV = SHV.rows();
  Face eFace(nbVert + nbSHV);
  for (int iVert = 0; iVert < nbVert; iVert++)
    eFace[iVert] = 1;
  for (int iSHV = 0; iSHV < nbSHV; iSHV++)
    eFace[nbVert + iSHV] = 0;
  Tgroup PreGRPisom =
      GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat, os);
  Tgroup GRPisom = ReducedGroupActionFace(PreGRPisom, eFace);
  Tgroup GRPlatt = LinPolytopeIntegral_Stabilizer_Method8(EXT_T, GRPisom, os);
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  os << "|DEL_ENUM: Delaunay_Stabilizer|=" << time << "\n";
#endif
  return GRPlatt;
}

template <typename T, typename Tint, typename Tgroup>
Tgroup Delaunay_Stabilizer(DataLattice<T, Tint, Tgroup> const &eData,
                           MyMatrix<Tint> const &EXT, std::ostream &os) {
  return Delaunay_StabilizerKernel<T, Tint, Tgroup>(eData.GramMat, eData.SHV,
                                                    EXT, os);
}

template <typename T> bool is_affine_integral(MyMatrix<T> const &M) {
  int dim = M.rows() - 1;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      T val = M(i + 1, j + 1);
      if (!IsInteger(val)) {
        return false;
      }
    }
  }
  for (int i = 0; i < dim; i++) {
    T val = M(i + 1, 0);
    if (val != 0) {
      return false;
    }
  }
  if (M(0, 0) != 1) {
    return false;
  }
  return true;
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<T>>
Polytope_TestEquivalence(DataLattice<T, Tint, Tgroup> &eData,
                         MyMatrix<T> const &EXT1_T, MyMatrix<T> const &EXT2_T) {
  std::ostream &os = eData.rddo.os;
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  MicrosecondTime time;
#endif
  using Telt = typename Tgroup::Telt;
  using Tgr = GraphListAdj;
  using Tidx_value = int16_t;
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: Begin Delaunay_TestEquivalence\n";
#endif
  //
  // Now extending by adding more vectors.
  //
  WeightMatrix<true, T, Tidx_value> WMat1 =
      GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT1_T, eData.GramMat,
                                                eData.SHV, os);
  WeightMatrix<true, T, Tidx_value> WMat2 =
      GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT2_T, eData.GramMat,
                                                eData.SHV, os);
  std::optional<Telt> eRes =
      TestEquivalenceWeightMatrix<T, Telt, Tidx_value>(WMat1, WMat2, os);
  if (!eRes) {
#ifdef DEBUG_DELAUNAY_ENUMERATION
    os << "DEL_ENUM: Leaving Delaunay_TestEquivalence 1 with false\n";
#endif
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    os << "|DEL_ENUM: Delaunay_TestEquivalence|=" << time << "\n";
#endif
    return {};
  }
  Telt const &eElt = *eRes;
  MyMatrix<T> MatEquiv_T = FindTransformation<T, Telt>(EXT1_T, EXT2_T, eElt);
  if (is_affine_integral(MatEquiv_T)) {
#ifdef DEBUG_DELAUNAY_ENUMERATION
    os << "DEL_ENUM: Leaving Delaunay_TestEquivalence 2 with true\n";
#endif
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    os << "|DEL_ENUM: Delaunay_TestEquivalence|=" << time << "\n";
#endif
    return MatEquiv_T;
  }
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: Trying other strategies\n";
#endif
  Tgroup GRP1 =
      GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat1, os);
  // We iterate over the group elements since we do not have yet have a
  // way to encode that the want integrality only on the linear part.
  for (auto &eElt1 : GRP1) {
    Telt fElt = eElt1 * eElt;
    MyMatrix<T> MatEquiv_T = FindTransformation<T, Telt>(EXT1_T, EXT2_T, fElt);
    if (is_affine_integral(MatEquiv_T)) {
      return MatEquiv_T;
    }
  }
  return {};
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
Delaunay_TestEquivalence(DataLattice<T, Tint, Tgroup> &eData,
                         MyMatrix<Tint> const &EXT1,
                         MyMatrix<Tint> const &EXT2) {
  std::ostream &os = eData.rddo.os;
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
      GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT1_T, eData.GramMat,
                                                eData.SHV, os);
  WeightMatrix<true, T, Tidx_value> WMat2 =
      GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT2_T, eData.GramMat,
                                                eData.SHV, os);
  std::optional<Telt> eRes =
      TestEquivalenceWeightMatrix<T, Telt, Tidx_value>(WMat1, WMat2, os);
  if (!eRes) {
#ifdef DEBUG_DELAUNAY_ENUMERATION
    os << "DEL_ENUM: Leaving Delaunay_TestEquivalence 1 with false\n";
#endif
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    os << "|DEL_ENUM: Delaunay_TestEquivalence|=" << time << "\n";
#endif
    return {};
  }
  Telt const &eElt = *eRes;
  MyMatrix<T> MatEquiv_T = FindTransformation<T, Telt>(EXT1_T, EXT2_T, eElt);
  if (IsIntegralMatrix(MatEquiv_T)) {
    MyMatrix<Tint> MatEquiv_I = UniversalMatrixConversion<Tint, T>(MatEquiv_T);
#ifdef DEBUG_DELAUNAY_ENUMERATION
    os << "DEL_ENUM: Leaving Delaunay_TestEquivalence 2 with true\n";
#endif
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    os << "|DEL_ENUM: Delaunay_TestEquivalence|=" << time << "\n";
#endif
    return MatEquiv_I;
  }
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: Trying other strategies\n";
#endif
  Tgroup GRP1 =
      GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat1, os);
  std::optional<MyMatrix<T>> eEq =
      LinPolytopeIntegral_Isomorphism_Method8(EXT1_T, EXT2_T, GRP1, eElt, os);
  if (!eEq) {
#ifdef DEBUG_DELAUNAY_ENUMERATION
    os << "DEL_ENUM: Leaving Delaunay_TestEquivalence 3 with false\n";
#endif
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    os << "|DEL_ENUM: Delaunay_TestEquivalence|=" << time << "\n";
#endif
    return {};
  }
  MyMatrix<Tint> eMat_I = UniversalMatrixConversion<Tint, T>(*eEq);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: Leaving Delaunay_TestEquivalence 4 with true\n";
#endif
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  os << "|DEL_ENUM: Delaunay_TestEquivalence|=" << time << "\n";
#endif
  return eMat_I;
}

template <typename T, typename Tint, typename Tgroup>
size_t ComputeInvariantDelaunay(DataLattice<T, Tint, Tgroup> const &eData,
                                size_t const &seed, MyMatrix<Tint> const &EXT,
                                [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  MicrosecondTime time;
#endif
  int nbVert = EXT.rows();
  int n = EXT.cols() - 1;
  Tint PreIndex = Int_IndexLattice(EXT);
  Tint eIndex = T_abs(PreIndex);
  MyMatrix<T> EXT_fullT = UniversalMatrixConversion<T, Tint>(EXT);
  CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral(eData.GramMat, EXT_fullT);
  MyMatrix<T> EXT_T(nbVert, n);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    for (int i = 0; i < n; i++) {
      T val = UniversalScalarConversion<T, Tint>(EXT(iVert, i + 1));
      EXT_T(iVert, i) = val - eCP.eCent(i + 1);
    }
  }
  std::map<T, size_t> ListDiagNorm;
  std::map<T, size_t> ListOffDiagNorm;
  MyVector<T> V(n);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    for (int i = 0; i < n; i++) {
      T eSum = 0;
      for (int j = 0; j < n; j++) {
        eSum += eData.GramMat(i, j) * EXT_T(iVert, j);
      }
      V(i) = eSum;
    }
    T scal = 0;
    for (int i = 0; i < n; i++) {
      scal += V(i) * EXT_T(iVert, i);
    }
    ListDiagNorm[scal] += 1;
    for (int jVert = iVert + 1; jVert < nbVert; jVert++) {
      T scal = 0;
      for (int i = 0; i < n; i++) {
        scal += V(i) * EXT_T(jVert, i);
      }
      ListOffDiagNorm[scal] += 1;
    }
  }
  size_t hash = ComputeHashTwoMap(seed, ListDiagNorm, ListOffDiagNorm);
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  os << "|DEL_ENUM: ComputeInvariantDelaunay|=" << time << "\n";
#endif
  return hash;
}

template <typename Tint> struct Delaunay_AdjI {
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
} // namespace boost::serialization

template <typename Tint> struct Delaunay_AdjO {
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
} // namespace boost::serialization

template <typename Tint> struct Delaunay_AdjO_spec {
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
} // namespace boost::serialization

template <typename Tint, typename Tgroup> struct Delaunay_Entry {
  MyMatrix<Tint> EXT;
  Tgroup GRP;
  std::vector<Delaunay_AdjO<Tint>> ListAdj;
};

namespace boost::serialization {
template <class Archive, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, Delaunay_Entry<Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("EXT", eRec.EXT);
  ar &make_nvp("GRP", eRec.GRP);
  ar &make_nvp("ListAdj", eRec.ListAdj);
}
} // namespace boost::serialization

template <typename Tvert, typename Tgroup> struct DelaunayTesselation {
  std::vector<Delaunay_Entry<Tvert, Tgroup>> l_dels;
};

namespace boost::serialization {
template <class Archive, typename Tvert, typename Tgroup>
inline void serialize(Archive &ar, DelaunayTesselation<Tvert, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("l_dels", eRec.l_dels);
}
} // namespace boost::serialization

template <typename Tvert, typename Tgroup>
void check_delaunay_tessellation(DelaunayTesselation<Tvert, Tgroup> const &DT,
                                 [[maybe_unused]] std::ostream &os) {
  for (auto &eDel : DT.l_dels) {
    MyMatrix<Tvert> const &EXT = eDel.EXT;
    ContainerMatrix<Tvert> cont(EXT);
    for (auto &eAdj : eDel.ListAdj) {
      Face const &eInc = eAdj.eInc;
      Face eIncEff(EXT.rows());
      MyMatrix<Tvert> const &EXT2 = DT.l_dels[eAdj.iOrb].EXT;
      MyMatrix<Tvert> const &eBigMat = eAdj.eBigMat;
#ifdef DEBUG_DELAUNAY_ENUMERATION
      os << "DEL_ENUM: check_delaunay_tessellation |EXT2}=" << EXT2.rows()
         << "/" << EXT2.cols() << " |eBigMat|=" << eBigMat.rows() << "/"
         << eBigMat.cols() << "\n";
#endif
      MyMatrix<Tvert> EXTadj = EXT2 * eBigMat;
      int len = EXTadj.rows();
      for (int u = 0; u < len; u++) {
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

template <typename Tvert, typename Tgroup>
void WriteEntryGAP(std::ostream &os_out,
                   DelaunayTesselation<Tvert, Tgroup> const &DT) {
  using Telt = typename Tgroup::Telt;
  os_out << "[";
  size_t n_del = DT.l_dels.size();
  for (size_t i_del = 0; i_del < n_del; i_del++) {
    Delaunay_Entry<Tvert, Tgroup> const &eDel = DT.l_dels[i_del];
    MyMatrix<Tvert> const &EXT = eDel.EXT;
    if (i_del > 0)
      os_out << ",";
    os_out << "rec(EXT:=" << StringMatrixGAP(EXT) << ",\n";
    std::vector<Telt> LGen = eDel.GRP.SmallGeneratingSet();
    if (LGen.size() == 0) {
      LGen.push_back(eDel.GRP.get_identity());
    }
    auto get_gap_string = [&]() -> std::string {
      if (LGen.size() == 0) {
        return "Group(())";
      } else {
        std::string str_ret = "Group([";
        bool IsFirst = true;
        for (auto &eElt : LGen) {
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
    for (auto &eElt : LGen) {
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
    std::string str_phi = "GroupHomomorphismByImagesNC(Group(" + str_perm +
                          "), Group(" + str_matr + "), " + str_perm + ", " +
                          str_matr + ")";
    os_out << "TheStab:=rec(PermutationStabilizer:=" << get_gap_string()
           << ", PhiPermMat:=" << str_phi << "), ";
    os_out << "Adjacencies:=[";
    IsFirst = true;
    for (auto &eAdj : eDel.ListAdj) {
      if (!IsFirst)
        os_out << ",";
      IsFirst = false;
      os_out << "rec(iDelaunay:=" << (eAdj.iOrb + 1) << ", ";
      os_out << "eInc:=[";
      std::vector<int> V = FaceToVector<int>(eAdj.eInc);
      for (size_t u = 0; u < V.size(); u++) {
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

template <typename Tvert, typename Tint, typename Tgroup>
void WriteDetailedEntryGAP(
    std::ostream &os_out,
    [[maybe_unused]] DataLattice<Tvert, Tint, Tgroup> const &data,
    DelaunayTesselation<Tvert, Tgroup> const &DT,
    [[maybe_unused]] std::ostream &os) {
  WriteEntryGAP(os_out, DT);
}

template <typename Tvert, typename Tgroup>
void WriteEntryPYTHON(std::ostream &os_out,
                      DelaunayTesselation<Tvert, Tgroup> const &DT) {
  using Telt = typename Tgroup::Telt;
  os_out << "[";
  size_t n_del = DT.l_dels.size();
  for (size_t i_del = 0; i_del < n_del; i_del++) {
    Delaunay_Entry<Tvert, Tgroup> const &eDel = DT.l_dels[i_del];
    MyMatrix<Tvert> const &EXT = eDel.EXT;
    int n_vert = EXT.rows();
    if (i_del > 0)
      os_out << ",";
    os_out << "{\"EXT\":" << StringMatrixPYTHON(EXT);
    std::vector<Telt> LGen = eDel.GRP.SmallGeneratingSet();
    auto get_gap_string = [&]() -> std::string {
      std::string str_ret = "[";
      bool IsFirst = true;
      for (auto &eElt : LGen) {
        if (!IsFirst) {
          str_ret += ",";
        }
        IsFirst = false;
        str_ret += "[";
        for (int i_vert = 0; i_vert < n_vert; i_vert++) {
          if (i_vert > 0) {
            str_ret += ",";
          }
          str_ret += std::to_string(OnPoints(i_vert, eElt));
        }
        str_ret += "]";
      }
      str_ret += "]";
      return str_ret;
    };
    os_out << ", \"TheStab\":" << get_gap_string();
    os_out << ", \"Adjacencies\":[";
    bool IsFirstAdj = true;
    for (auto &eAdj : eDel.ListAdj) {
      if (!IsFirstAdj)
        os_out << ",";
      IsFirstAdj = false;
      os_out << "{\"iDelaunay\":" << eAdj.iOrb << ", \"eInc\":[";
      std::vector<int> V = FaceToVector<int>(eAdj.eInc);
      for (size_t u = 0; u < V.size(); u++) {
        if (u > 0) {
          os_out << ",";
        }
        os_out << V[u];
      }
      os_out << "], \"eBigMat\":" << StringMatrixPYTHON(eAdj.eBigMat) << "}";
    }
    os_out << "]}";
  }
  os_out << "]";
}

template <typename T, typename Tint, typename Tgroup>
std::pair<Tgroup, std::vector<Delaunay_AdjI<Tint>>>
ComputeGroupAndAdjacencies(DataLattice<T, Tint, Tgroup> &eData,
                           MyMatrix<Tint> const &x) {
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(x);
  std::ostream &os = eData.rddo.os;
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: |EXT_T|=" << EXT_T.rows() << " / " << EXT_T.cols() << "\n";
#endif
  Tgroup GRPlatt = Delaunay_Stabilizer<T, Tint, Tgroup>(eData, x, os);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: |GRPlatt|=" << GRPlatt.size() << "\n";
#endif
  vectface TheOutput = DualDescriptionRecordFullDim(EXT_T, GRPlatt, eData.rddo);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: |TheOutput|=" << TheOutput.size() << "\n";
#endif
  std::vector<Delaunay_AdjI<Tint>> ListAdj;
  for (auto &eOrbB : TheOutput) {
    MyMatrix<Tint> EXTadj = FindAdjacentDelaunayPolytope<T, Tint>(
        eData.GramMat, eData.solver, eData.ShvGraverBasis, EXT_T, eOrbB, os);
    Delaunay_AdjI<Tint> eAdj{eOrbB, EXTadj};
    ListAdj.push_back(eAdj);
  }
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: |ListAdj|=" << ListAdj.size() << "\n";
#endif
  return {GRPlatt, std::move(ListAdj)};
}

template <typename Tint, typename Tgroup> struct Delaunay_Obj {
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
} // namespace boost::serialization

template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint>
FindDelaunayPolytopeExtended(DataLattice<T, Tint, Tgroup> &data) {
  MyMatrix<T> const &GramMat = data.GramMat;
  CVPSolver<T, Tint> &solver = data.solver;
  MyMatrix<Tint> const &ShvGraverBasis = data.ShvGraverBasis;
  std::string const &choice = data.choice_initial;
  std::ostream &os = data.rddo.os;
  if (choice == "direct") {
    return FindDelaunayPolytope<T, Tint>(GramMat, solver, os);
  }
  if (choice == "simple_iter") {
    int target_ext = 2 * GramMat.rows();
    int max_iter = 10;
    std::string method = "lp_cdd";
    return FindDelaunayPolytope_random<T, Tint>(
        GramMat, solver, ShvGraverBasis, target_ext, max_iter, method, os);
  }
  if (choice == "sampling") {
    int target_ext = 0;
    int max_iter = std::numeric_limits<int>::max();
    std::string method = "sampling";
    return FindDelaunayPolytope_random<T, Tint>(
        GramMat, solver, ShvGraverBasis, target_ext, max_iter, method, os);
  }
  std::cerr
      << "Failed to find a relevant method for FindDelaunayPolytopeExtended\n";
  std::cerr << "Allowed methods: direct, simple_iter, sampling\n";
  throw TerminalException{1};
}

template <typename T, typename Tint, typename Tgroup> struct DataLatticeFunc {
  DataLattice<T, Tint, Tgroup> data;
  using Tobj = Delaunay_Obj<Tint, Tgroup>;
  using TadjI = Delaunay_AdjI<Tint>;
  using TadjO = Delaunay_AdjO_spec<Tint>;
  std::ostream &get_os() { return data.rddo.os; }
  Tobj f_init() {
#ifdef DEBUG_DELAUNAY_ENUMERATION
    data.rddo.os << "DEL_ENUM: DataLatticeFunc : f_init, OutFile="
                 << data.rddo.AllArr.OutFile << "\n";
#endif
    MyMatrix<Tint> EXT = FindDelaunayPolytopeExtended<T, Tint>(data);
    Tobj x{std::move(EXT), {}};
    return x;
  }
  size_t f_hash(size_t const &seed, Tobj const &x) {
    return ComputeInvariantDelaunay(data, seed, x.EXT, data.rddo.os);
  }
  std::optional<TadjO> f_repr(Tobj const &x, TadjI const &y) {
    std::optional<MyMatrix<Tint>> opt =
        Delaunay_TestEquivalence<T, Tint, Tgroup>(data, x.EXT, y.EXT);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const &eBigMat = *opt;
    TadjO ret{y.eInc, eBigMat};
    return ret;
  }
  std::pair<Tobj, TadjO> f_spann(TadjI const &x) {
    MyMatrix<Tint> EXT = x.EXT;
    Tobj x_ret{EXT, {}};
    MyMatrix<Tint> eBigMat = IdentityMat<Tint>(data.n + 1);
    TadjO ret{x.eInc, eBigMat};
    return {x_ret, ret};
  }
  std::vector<TadjI> f_adj(Tobj &x) {
    std::pair<Tgroup, std::vector<TadjI>> pair =
        ComputeGroupAndAdjacencies<T, Tint, Tgroup>(data, x.EXT);
    x.GRP = pair.first;
    return pair.second;
  }
  Tobj f_adji_obj(TadjI const &x) {
    MyMatrix<Tint> EXT = x.EXT;
    Tobj x_ret{EXT, {}};
    return x_ret;
  }
  size_t f_complexity(Tobj const &x) { return x.EXT.rows(); }
};

template <typename T, typename Tvert, typename Tgroup>
DelaunayTesselation<Tvert, Tgroup>
DelaunayTesselation_From_DatabaseEntries_Serial(
    std::vector<DatabaseEntry_Serial<
        typename DataLatticeFunc<T, Tvert, Tgroup>::Tobj,
        typename DataLatticeFunc<T, Tvert, Tgroup>::TadjO>> const &l_ent) {
  std::vector<Delaunay_Entry<Tvert, Tgroup>> l_dels;
  for (auto &eDel : l_ent) {
    std::vector<Delaunay_AdjO<Tvert>> ListAdj;
    for (auto &eAdj : eDel.ListAdj) {
      Delaunay_AdjO<Tvert> fAdj{eAdj.x.eInc, eAdj.x.eBigMat, eAdj.iOrb};
      ListAdj.push_back(fAdj);
    }
    Delaunay_Entry<Tvert, Tgroup> fDel{eDel.x.EXT, eDel.x.GRP, ListAdj};
    l_dels.push_back(fDel);
  }
  return {l_dels};
}

template <typename T, typename Tvert, typename Tgroup>
void WriteDelaunayTesselation(std::string const &OutFormat,
                              std::ostream &os_out, MyMatrix<T> const &GramMat,
                              DelaunayTesselation<Tvert, Tgroup> const &DT) {
  if (OutFormat == "nothing") {
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return ";
    WriteEntryGAP(os_out, DT);
    os_out << ";\n";
    return;
  }
  if (OutFormat == "PYTHON") {
    return WriteEntryPYTHON(os_out, DT);
  }
  if (OutFormat == "GAP_covering") {
    T TheCovSqr(0);
    for (auto &eDel : DT.l_dels) {
      MyMatrix<Tvert> const &EXT = eDel.EXT;
      MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tvert>(EXT);
      CP<T> cp = CenterRadiusDelaunayPolytopeGeneral<T>(GramMat, EXT_T);
      T SquareRadius = cp.SquareRadius;
      if (SquareRadius > TheCovSqr) {
        TheCovSqr = SquareRadius;
      }
    }
    T TheDet = DeterminantMat(GramMat);
    int TheDim = GramMat.rows();
    ResultCov<T> x =
        ComputeCoveringDensityFromDimDetCov<T>(TheDim, TheDet, TheCovSqr);
    os_out << "return " << to_stringGAP(x) << ";\n";
    return;
  }
  std::cerr << "DEL_ENUM: WriteDelaunayTesselation failed for OutFormat="
            << OutFormat << "\n";
  throw TerminalException{1};
}

template <typename T, typename Tint, typename Tgroup, typename Fincorrect>
std::optional<DelaunayTesselation<Tint, Tgroup>>
EnumerationDelaunayPolytopes(DataLattice<T, Tint, Tgroup> &data,
                             Fincorrect f_incorrect,
                             int const &max_runtime_second) {
#ifdef DEBUG_DELAUNAY_ENUMERATION
  std::ostream &os = data.rddo.os;
  os << "DEL_ENUM: EnumerationDelaunayPolytopes, begin\n";
  os << "DEL_ENUM: EnumerationDelaunayPolytopes, OutFile="
     << data.rddo.AllArr.OutFile << "\n";
#endif
  using Tdata = DataLatticeFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = std::vector<DatabaseEntry_Serial<Tobj, TadjO>>;
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: EnumerationDelaunayPolytopes, before "
        "EnumerateAndStore_Serial\n";
#endif
  bool is_incorrect = false;
  auto f_incorrect_bis = [&](Tobj const &x) -> bool {
    bool test = f_incorrect(x);
    if (test) {
      is_incorrect = true;
    }
    return test;
  };
  Tout result = EnumerateAndStore_Serial<Tdata>(data_func, f_incorrect_bis,
                                                max_runtime_second);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: EnumerationDelaunayPolytopes, after "
        "EnumerateAndStore_Serial\n";
#endif
  if (is_incorrect) {
#ifdef DEBUG_DELAUNAY_ENUMERATION
    os << "DEL_ENUM: EnumerationDelaunayPolytopes: opt not matching\n";
#endif
    return {};
  }
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: EnumerationDelaunayPolytopes: opt match\n";
#endif
  DelaunayTesselation<Tint, Tgroup> DT =
      DelaunayTesselation_From_DatabaseEntries_Serial<T, Tint, Tgroup>(result);
  return DT;
}

FullNamelist NAMELIST_GetStandard_COMPUTE_DELAUNAY() {
  std::map<std::string, SingleBlock> ListBlock;
  // SYSTEM
  ListBlock["SYSTEM"] = SINGLEBLOCK_Get_System();
  // DATA
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["arithmetic_T"] = "gmp_rational";
  ListStringValues1["arithmetic_Tint"] = "gmp_integer";
  ListStringValues1["GRAMfile"] = "unset.gram";
  ListStringValues1["SVRfile"] = "unset.svr";
  ListStringValues1["choice_initial"] = "direct";
  ListStringValues1["FileDualDescription"] = "unset";
  SingleBlock BlockDATA;
  BlockDATA.setListDoubleValues(ListDoubleValues1);
  BlockDATA.setListStringValues(ListStringValues1);
  BlockDATA.setListListStringValues(ListListStringValues1);
  ListBlock["DATA"] = BlockDATA;
  // Merging all data
  return FullNamelist(ListBlock);
}

template <typename T, typename Tint, typename Tgroup>
DataLattice<T, Tint, Tgroup>
GetDataLattice(MyMatrix<T> const &GramMat,
               PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
               std::ostream &os) {
  int n = GramMat.rows();
  MyMatrix<T> SHV(0, n);
  CVPSolver<T, Tint> solver(GramMat, os);
  std::string choice_initial = "direct";
  MyMatrix<Tint> ShvGraverBasis = GetGraverBasis<T, Tint>(GramMat);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: GetDataLattice, AllArr.OutFile=" << AllArr.OutFile << "\n";
#endif
  return {n,
          GramMat,
          SHV,
          solver,
          ShvGraverBasis,
          choice_initial,
          RecordDualDescOperation<T, Tgroup>(AllArr, os)};
}

template <typename T, typename Tint, typename Tgroup>
DataLattice<T, Tint, Tgroup>
get_data_lattice(FullNamelist const &eFull,
                 PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
                 std::ostream &os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");

  std::string GRAMfile = BlockDATA.get_string("GRAMfile");
  MyMatrix<T> GramMat = ReadMatrixFile<T>(GRAMfile);
  //
  std::string SVRfile = BlockDATA.get_string("SVRfile");
  auto get_SVR = [&]() -> MyMatrix<T> {
    if (IsExistingFile(SVRfile)) {
      return ReadMatrixFile<T>(SVRfile);
    }
    int n = GramMat.rows();
    return ZeroMatrix<T>(0, n);
  };
  MyMatrix<T> SVR = get_SVR();
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: |SVR|=" << SVR.rows() << "\n";
#endif
  //
  std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: OutFormat=" << OutFormat << " OutFile=" << OutFile << "\n";
#endif
  //
  int n = GramMat.rows();
  //
  CVPSolver<T, Tint> solver(GramMat, os);
  MyMatrix<Tint> ShvGraverBasis = GetGraverBasis<T, Tint>(GramMat);
  std::string choice_initial = BlockDATA.get_string("choice_initial");
  DataLattice<T, Tint, Tgroup> data{
      n,
      GramMat,
      SVR,
      solver,
      ShvGraverBasis,
      choice_initial,
      RecordDualDescOperation<T, Tgroup>(AllArr, os)};
  return data;
}

// clang-format off
#endif  // SRC_LATT_LATTICEDELAUNAY_H_
// clang-format on

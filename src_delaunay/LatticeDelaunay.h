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
  MyMatrix<T> SHV;
  CVPSolver<T, Tint> solver;
  MyMatrix<Tint> ShvGraverBasis;
  std::string choice_initial;
  RecordDualDescOperation<T, Tgroup> rddo;
};

template <typename T>
MyMatrix<T> get_reduced_delaunay_shv(MyMatrix<T> const &EXT, MyMatrix<T> const &GramMat,
                                     MyMatrix<T> const &SHV, MyVector<T> const& TheCenter) {
  int n = GramMat.rows();
  int nbVect = SHV.rows();
  int nbVert = EXT.rows();
  MyMatrix<T> EXText(nbVect + nbVert, n);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    for (int i=0; i<n; i++) {
      EXText(iVert, i) = EXT(iVert, i+1) - TheCenter(i);
    }
  }
  for (int iVect = 0; iVect < nbVect; iVect++) {
    for (int i = 0; i < n; i++) {
      EXText(nbVert + iVect, i) = SHV(iVect, i);
    }
  }
  return EXText;
}

template <typename T>
Face get_face_delaunay_shv(MyMatrix<T> const &EXT, MyMatrix<T> const &SHV) {
  int nbVect = SHV.rows();
  int nbVert = EXT.rows();
  Face eFace(nbVert + nbVect);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    eFace[iVert] = 1;
  }
  for (int iVect = 0; iVect < nbVect; iVect++) {
    eFace[nbVert + iVect] = 0;
  }
  return eFace;
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

template <typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value>
GetWeightMatrixFromGramEXT(MyMatrix<T> const &EXT, MyMatrix<T> const &GramMat,
                           MyMatrix<T> const &SHV, std::ostream &os) {
  CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral(GramMat, EXT);
  MyMatrix<T> EXText = get_reduced_delaunay_shv(EXT, GramMat, SHV, eCP.eCent);
  return GetSimpleWeightMatrix<T, Tidx_value>(EXText, GramMat, os);
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

template <typename T, typename Tint, typename Tgroup, typename Fcent>
Tgroup PolytopeGen_StabilizerKernel(MyMatrix<T> const &GramMat,
                                    Fcent f_cent,
                                    MyMatrix<T> const &SHV,
                                    MyMatrix<T> const &EXT_T, std::ostream &os) {
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  MicrosecondTime time;
#endif
  using Tidx_value = int16_t;
  using Tgr = GraphListAdj;
  //
  // Now extending with the SHV vector set
  //
  MyVector<T> Cent = f_cent(EXT_T);
  MyMatrix<T> EXText = get_reduced_delaunay_shv(EXT_T, GramMat, SHV, Cent);
  WeightMatrix<true, T, Tidx_value> WMat =
    GetSimpleWeightMatrix<T, Tidx_value>(EXText, GramMat, os);
  Face eFace = get_face_delaunay_shv(EXT_T, SHV);
  Tgroup PreGRPisom =
      GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat, os);
  Tgroup GRPisom = ReducedGroupActionFace(PreGRPisom, eFace);
  MyMatrix<T> EXTextInt = RemoveFractionMatrix(EXText);
  Tgroup GRPlatt = LinPolytopeIntegral_Stabilizer(EXTextInt, GRPisom, os);
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  os << "|DEL_ENUM: Delaunay_Stabilizer|=" << time << "\n";
#endif
  return GRPlatt;
}

template <typename T, typename Tint, typename Tgroup>
Tgroup Polytope_StabilizerKernel(MyMatrix<T> const &GramMat,
                                 MyMatrix<T> const &SHV,
                                 MyMatrix<T> const &EXT_T, std::ostream &os) {
  auto f_cent=[&](MyMatrix<T> const& EXT_T) -> MyVector<T> {
    CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral(GramMat, EXT_T);
    return eCP.eCent;
  };
  return PolytopeGen_StabilizerKernel<T,Tint,Tgroup, decltype(f_cent)>(GramMat, f_cent, SHV, EXT_T, os);
}



template <typename T, typename Tint, typename Tgroup>
Tgroup Delaunay_Stabilizer(DataLattice<T, Tint, Tgroup> const &eData,
                           MyMatrix<T> const &EXT_T, std::ostream &os) {
  auto f_cent=[&](MyMatrix<T> const& EXT_T) -> MyVector<T> {
    CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral(eData.solver.GramMat, EXT_T);
    return eCP.eCent;
  };
  return PolytopeGen_StabilizerKernel<T, Tint, Tgroup, decltype(f_cent)>(eData.solver.GramMat, f_cent, eData.SHV, EXT_T, os);
}

template <typename T, typename Tint, typename Tgroup, typename Fcent>
std::optional<MyMatrix<T>>
PolytopeGen_TestEquivalence(DataLattice<T, Tint, Tgroup> &eData,
                            Fcent f_cent,
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
  MyMatrix<T> const& GramMat = eData.solver.GramMat;
  //
  // Now extending by adding more vectors.
  //
  MyVector<T> Cent1 = f_cent(EXT1_T);
  MyVector<T> Cent2 = f_cent(EXT2_T);
  auto extend_linear_transform=[&](MyMatrix<T> const& M) -> MyMatrix<T> {
    int n = GramMat.rows();
    MyMatrix<T> eMatRet = ZeroMatrix<T>(n+1, n+1);
    eMatRet(0,0) = 1;
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        eMatRet(i+1, j+1) = M(i,j);
      }
    }
    MyVector<T> delta = Cent2 - M.transpose() * Cent1;
    for (int i = 0; i<n; i++) {
      eMatRet(0, i+1) = delta(i);
    }
    return eMatRet;
  };
  MyMatrix<T> EXText1 = get_reduced_delaunay_shv(EXT1_T, GramMat, eData.SHV, Cent1);
  MyMatrix<T> EXText2 = get_reduced_delaunay_shv(EXT2_T, GramMat, eData.SHV, Cent2);
  WeightMatrix<true, T, Tidx_value> WMat1 =
    GetSimpleWeightMatrix<T, Tidx_value>(EXText1, GramMat, os);
  WeightMatrix<true, T, Tidx_value> WMat2 =
    GetSimpleWeightMatrix<T, Tidx_value>(EXText2, GramMat, os);
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
  Face eFace = get_face_delaunay_shv(EXT1_T, eData.SHV);
  Telt const &eElt = *eRes;
  Telt eEltRed = ReduceElementActionFace(eElt, eFace);
  MyMatrix<T> MatEquiv_T = FindTransformation<T, Telt>(EXText1, EXText2, eElt);
  if (IsIntegralMatrix(MatEquiv_T)) {
#ifdef DEBUG_DELAUNAY_ENUMERATION
    os << "DEL_ENUM: Leaving Delaunay_TestEquivalence 2 with true\n";
#endif
    MyMatrix<T> eMatRet = extend_linear_transform(MatEquiv_T);
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    os << "|DEL_ENUM: Delaunay_TestEquivalence|=" << time << "\n";
#endif
    return eMatRet;
  }
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: Trying other strategy\n";
#endif
  Tgroup GRP1 =
      GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat1, os);
  Tgroup GRPisom1 =
      GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat1, os);
  MyMatrix<T> EXTextInt1 = RemoveFractionMatrix(EXText1);
  MyMatrix<T> EXTextInt2 = RemoveFractionMatrix(EXText2);
  std::optional<MyMatrix<T>> opt =
    LinPolytopeIntegral_Isomorphism<T,Tgroup>(EXTextInt1, EXTextInt2, GRPisom1, eElt, os);
  if (!opt) {
#ifdef DEBUG_DELAUNAY_ENUMERATION
    os << "DEL_ENUM: Leaving Delaunay_TestEquivalence 3 with false\n";
#endif
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    os << "|DEL_ENUM: Delaunay_TestEquivalence|=" << time << "\n";
#endif
    return {};
  }
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: Leaving Delaunay_TestEquivalence 4 with true\n";
#endif
  MyMatrix<T> eMatRet = extend_linear_transform(*opt);
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  os << "|DEL_ENUM: Delaunay_TestEquivalence|=" << time << "\n";
#endif
  return eMatRet;
}

template<typename T>
MyVector<T> get_reduced_isobarycenter(MyMatrix<T> const& EXT) {
  int n_col = EXT.cols();
  int n_row = EXT.rows();
  MyVector<T> eIso(n_col-1);
  for (int i_col=0; i_col<n_col-1; i_col++) {
    T scal(0);
    for (int i_row=0; i_row<n_row; i_row++) {
      scal += EXT(i_row, i_col+1);
    }
    eIso(i_col) = scal / n_row;
  }
  return eIso;
}


template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<T>>
Polytope_TestEquivalence(DataLattice<T, Tint, Tgroup> &eData,
                         MyMatrix<T> const &EXT1_T, MyMatrix<T> const &EXT2_T) {
  auto f_cent=[&](MyMatrix<T> const& EXT) -> MyVector<T> {
    return get_reduced_isobarycenter(EXT);
  };
  return PolytopeGen_TestEquivalence<T,Tint,Tgroup,decltype(f_cent)>(eData, f_cent, EXT1_T, EXT2_T);
}



template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<T>>
Delaunay_TestEquivalence(DataLattice<T, Tint, Tgroup> &eData,
                         MyMatrix<T> const &EXT1_T,
                         MyMatrix<T> const &EXT2_T) {
  auto f_cent=[&](MyMatrix<T> const& EXT_T) -> MyVector<T> {
    CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral(eData.solver.GramMat, EXT_T);
    return eCP.eCent;
  };
  return PolytopeGen_TestEquivalence<T,Tint,Tgroup,decltype(f_cent)>(eData, f_cent, EXT1_T, EXT2_T);
}

template <typename T, typename Tint, typename Tgroup>
size_t ComputeInvariantDelaunay(DataLattice<T, Tint, Tgroup> const &eData,
                                size_t const &seed, MyMatrix<T> const &EXT_T,
                                [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  MicrosecondTime time;
#endif
  int nbVert = EXT_T.rows();
  int n = EXT_T.cols() - 1;
  MyMatrix<T> const& GramMat = eData.solver.GramMat;
  CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral(GramMat, EXT_T);
  MyMatrix<T> EXTtrans_T(nbVert, n);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    for (int i = 0; i < n; i++) {
      T val = EXT_T(iVert, i + 1);
      EXTtrans_T(iVert, i) = val - eCP.eCent(i + 1);
    }
  }
  std::map<T, size_t> ListDiagNorm;
  std::map<T, size_t> ListOffDiagNorm;
  MyVector<T> V(n);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    for (int i = 0; i < n; i++) {
      T eSum(0);
      for (int j = 0; j < n; j++) {
        eSum += GramMat(i, j) * EXTtrans_T(iVert, j);
      }
      V(i) = eSum;
    }
    T scal(0);
    for (int i = 0; i < n; i++) {
      scal += V(i) * EXTtrans_T(iVert, i);
    }
    ListDiagNorm[scal] += 1;
    for (int jVert = iVert + 1; jVert < nbVert; jVert++) {
      T scal(0);
      for (int i = 0; i < n; i++) {
        scal += V(i) * EXTtrans_T(jVert, i);
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

template <typename T> struct Delaunay_AdjI {
  Face eInc;
  MyMatrix<T> EXT;
};

namespace boost::serialization {
template <class Archive, typename T>
inline void serialize(Archive &ar, Delaunay_AdjI<T> &eRec,
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

template <typename T, typename Tgroup> struct DelaunayTesselation {
  std::vector<Delaunay_Entry<T, Tgroup>> l_dels;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tgroup>
inline void serialize(Archive &ar, DelaunayTesselation<T, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("l_dels", eRec.l_dels);
}
} // namespace boost::serialization

template <typename T, typename Tgroup>
void check_delaunay_tessellation(DelaunayTesselation<T, Tgroup> const &DT,
                                 [[maybe_unused]] std::ostream &os) {
  for (auto &eDel : DT.l_dels) {
    MyMatrix<T> const &EXT = eDel.EXT;
    ContainerMatrix<T> cont(EXT);
    for (auto &eAdj : eDel.ListAdj) {
      Face const &eInc = eAdj.eInc;
      Face eIncEff(EXT.rows());
      MyMatrix<T> const &EXT2 = DT.l_dels[eAdj.iOrb].EXT;
      MyMatrix<T> const &eBigMat = eAdj.eBigMat;
#ifdef DEBUG_DELAUNAY_ENUMERATION
      os << "DEL_ENUM: check_delaunay_tessellation |EXT2}=" << EXT2.rows()
         << "/" << EXT2.cols() << " |eBigMat|=" << eBigMat.rows() << "/"
         << eBigMat.cols() << "\n";
#endif
      MyMatrix<T> EXTadj = EXT2 * eBigMat;
      int len = EXTadj.rows();
      for (int u = 0; u < len; u++) {
        MyVector<T> V = GetMatrixRow(EXTadj, u);
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

template <typename T, typename Tgroup>
void WriteEntryGAP(std::ostream &os_out,
                   DelaunayTesselation<T, Tgroup> const &DT) {
  using Telt = typename Tgroup::Telt;
  os_out << "[";
  size_t n_del = DT.l_dels.size();
  for (size_t i_del = 0; i_del < n_del; i_del++) {
    Delaunay_Entry<T, Tgroup> const &eDel = DT.l_dels[i_del];
    MyMatrix<T> const &EXT = eDel.EXT;
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
      MyMatrix<T> M = RepresentVertexPermutation(EXT, EXT, eElt);
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

template <typename T, typename Tint, typename Tgroup>
void WriteDetailedEntryGAP(
    std::ostream &os_out,
    [[maybe_unused]] DataLattice<T, Tint, Tgroup> const &data,
    DelaunayTesselation<T, Tgroup> const &DT,
    [[maybe_unused]] std::ostream &os) {
  WriteEntryGAP(os_out, DT);
}

template <typename T, typename Tgroup>
void WriteEntryPYTHON(std::ostream &os_out,
                      DelaunayTesselation<T, Tgroup> const &DT) {
  using Telt = typename Tgroup::Telt;
  os_out << "[";
  size_t n_del = DT.l_dels.size();
  for (size_t i_del = 0; i_del < n_del; i_del++) {
    Delaunay_Entry<T, Tgroup> const &eDel = DT.l_dels[i_del];
    MyMatrix<T> const &EXT = eDel.EXT;
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
std::pair<Tgroup, std::vector<Delaunay_AdjI<T>>>
ComputeGroupAndAdjacencies(DataLattice<T, Tint, Tgroup> &eData,
                           MyMatrix<T> const &EXT_T) {
  std::ostream &os = eData.rddo.os;
  MyMatrix<T> const& GramMat = eData.solver.GramMat;
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: |EXT_T|=" << EXT_T.rows() << " / " << EXT_T.cols() << "\n";
#endif
  Tgroup GRPlatt = Delaunay_Stabilizer<T, Tint, Tgroup>(eData, EXT_T, os);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: |GRPlatt|=" << GRPlatt.size() << "\n";
#endif
  vectface TheOutput = DualDescriptionRecordFullDim(EXT_T, GRPlatt, eData.rddo);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: |TheOutput|=" << TheOutput.size() << "\n";
#endif
  std::vector<Delaunay_AdjI<T>> ListAdj;
  for (auto &eOrbB : TheOutput) {
    MyMatrix<Tint> EXTadj = FindAdjacentDelaunayPolytope<T, Tint>(
        GramMat, eData.solver, eData.ShvGraverBasis, EXT_T, eOrbB, os);
    MyMatrix<T> EXTadj_T = UniversalMatrixConversion<T,Tint>(EXTadj);
    Delaunay_AdjI<T> eAdj{eOrbB, EXTadj_T};
    ListAdj.push_back(eAdj);
  }
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "DEL_ENUM: |ListAdj|=" << ListAdj.size() << "\n";
#endif
  return {GRPlatt, std::move(ListAdj)};
}

template <typename T, typename Tgroup> struct Delaunay_Obj {
  MyMatrix<T> EXT;
  Tgroup GRP;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tgroup>
inline void serialize(Archive &ar, Delaunay_Obj<T, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("EXT", eRec.EXT);
  ar &make_nvp("GRP", eRec.GRP);
}
} // namespace boost::serialization

template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint>
FindDelaunayPolytopeExtended(DataLattice<T, Tint, Tgroup> &data) {
  MyMatrix<T> const &GramMat = data.solver.GramMat;
  CVPSolver<T, Tint> const& solver = data.solver;
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
  using Tobj = Delaunay_Obj<T, Tgroup>;
  using TadjI = Delaunay_AdjI<T>;
  using TadjO = Delaunay_AdjO_spec<T>;
  std::ostream &get_os() { return data.rddo.os; }
  Tobj f_init() {
#ifdef DEBUG_DELAUNAY_ENUMERATION
    data.rddo.os << "DEL_ENUM: DataLatticeFunc : f_init, OutFile="
                 << data.rddo.AllArr.OutFile << "\n";
#endif
    MyMatrix<Tint> EXT = FindDelaunayPolytopeExtended<T, Tint>(data);
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T,Tint>(EXT);
    Tobj x{std::move(EXT_T), {}};
    return x;
  }
  size_t f_hash(size_t const &seed, Tobj const &x) {
    return ComputeInvariantDelaunay(data, seed, x.EXT, data.rddo.os);
  }
  std::optional<TadjO> f_repr(Tobj const &x, TadjI const &y) {
    std::optional<MyMatrix<T>> opt =
        Delaunay_TestEquivalence<T, Tint, Tgroup>(data, x.EXT, y.EXT);
    if (!opt) {
      return {};
    }
    MyMatrix<T> const &eBigMat = *opt;
    TadjO ret{y.eInc, eBigMat};
    return ret;
  }
  std::pair<Tobj, TadjO> f_spann(TadjI const &x) {
    MyMatrix<T> EXT = x.EXT;
    Tobj x_ret{EXT, {}};
    MyMatrix<T> eBigMat = IdentityMat<T>(data.n + 1);
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
    MyMatrix<T> EXT = x.EXT;
    Tobj x_ret{EXT, {}};
    return x_ret;
  }
  size_t f_complexity(Tobj const &x) { return x.EXT.rows(); }
};

template <typename T, typename Tint, typename Tgroup>
DelaunayTesselation<T, Tgroup>
DelaunayTesselation_From_DatabaseEntries_Serial(
    std::vector<DatabaseEntry_Serial<
    typename DataLatticeFunc<T, Tint, Tgroup>::Tobj,
    typename DataLatticeFunc<T, Tint, Tgroup>::TadjO>> const &l_ent) {
  std::vector<Delaunay_Entry<T, Tgroup>> l_dels;
  for (auto &eDel : l_ent) {
    std::vector<Delaunay_AdjO<T>> ListAdj;
    for (auto &eAdj : eDel.ListAdj) {
      Delaunay_AdjO<T> fAdj{eAdj.x.eInc, eAdj.x.eBigMat, eAdj.iOrb};
      ListAdj.push_back(fAdj);
    }
    Delaunay_Entry<T, Tgroup> fDel{eDel.x.EXT, eDel.x.GRP, ListAdj};
    l_dels.push_back(fDel);
  }
  return {l_dels};
}

template <typename T, typename Tgroup>
void WriteDelaunayTesselation(std::string const &OutFormat,
                              std::ostream &os_out, MyMatrix<T> const &GramMat,
                              DelaunayTesselation<T, Tgroup> const &DT) {
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
      MyMatrix<T> const &EXT = eDel.EXT;
      CP<T> cp = CenterRadiusDelaunayPolytopeGeneral<T>(GramMat, EXT);
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
std::optional<DelaunayTesselation<T, Tgroup>>
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
  DelaunayTesselation<T, Tgroup> DT =
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

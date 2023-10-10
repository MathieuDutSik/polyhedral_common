// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_SAMPLINGFACET_H_
#define SRC_POLY_POLY_SAMPLINGFACET_H_

// clang-format off
#include "POLY_LinearProgramming.h"
#include "POLY_PolytopeFct.h"
#include "POLY_lrslib.h"
#include <map>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_SAMPLING_FACET
#endif

#ifdef TIMINGS
#define TIMINGS_SAMPLING_FACET
#endif

struct recSamplingOption {
  int critlevel;
  int maxnbcall;
  int maxnbsize;
  std::string prog;
};

template <typename T>
vectface Kernel_DUALDESC_SamplingFacetProcedure(
    MyMatrix<T> const &EXT, recSamplingOption const &eOption, int &nbCall) {
  int dim = RankMat(EXT);
  int len = EXT.rows();
  std::string prog = eOption.prog;
  int critlevel = eOption.critlevel;
  int maxnbcall = eOption.maxnbcall;
  int maxnbsize = eOption.maxnbsize;
  std::cerr << "critlevel=" << critlevel << " prog=" << prog
            << " maxnbcall=" << maxnbcall << "\n";
  auto IsRecursive = [&]() -> bool {
    if (len < critlevel)
      return false;
    if (dim < 15)
      return false;
    return true;
  };
  bool DoRecur = IsRecursive();
  vectface ListFace(EXT.rows());
  std::vector<int> ListStatus;
  auto FuncInsert = [&](Face const &eFace) -> void {
    for (auto &fFace : ListFace) {
      if (fFace.count() == eFace.count())
        return;
    }
    ListFace.push_back(eFace);
    ListStatus.push_back(0);
  };
  std::cerr << "dim=" << dim << "  len=" << len << "\n";
  if (!DoRecur) {
    auto comp_dd = [&]() -> vectface {
      if (prog == "lrs")
        return lrs::DualDescription_incd(EXT);
      if (prog == "cdd")
        return cdd::DualDescription_incd(EXT);
      std::cerr << "Failed to find a matching program\n";
      throw TerminalException{1};
    };
    vectface ListIncd = comp_dd();
    for (auto &eFace : ListIncd)
      FuncInsert(eFace);
    std::cerr << "DirectDualDesc |ListFace|=" << ListFace.size() << "\n";
    nbCall++;
    return ListFace;
  }
  Face eInc = FindOneInitialVertex(EXT);
  FuncInsert(eInc);
  while (true) {
    int nbCases = ListFace.size();
    bool IsFinished = true;
    for (int iC = 0; iC < nbCases; iC++)
      if (ListStatus[iC] == 0) {
        // we liberally increase the nbCall value
        nbCall++;
        IsFinished = false;
        ListStatus[iC] = 1;
        Face eFace = ListFace[iC];
        MyMatrix<T> EXTred = SelectRow(EXT, eFace);
        vectface ListRidge =
            Kernel_DUALDESC_SamplingFacetProcedure(EXTred, eOption, nbCall);
        for (auto &eRidge : ListRidge) {
          Face eFlip = ComputeFlipping(EXT, eFace, eRidge);
          FuncInsert(eFlip);
        }
        if (maxnbsize != -1) {
          int siz = ListFace.size();
          if (maxnbsize > siz) {
            std::cerr << "Ending by maxsize criterion\n";
            std::cerr << "siz=" << siz << " maxnbsize=" << maxnbsize << "\n";
            return ListFace;
          }
        }
        if (maxnbcall != -1) {
          if (nbCall > maxnbcall) {
            std::cerr << "Ending by maxnbcall\n";
            return ListFace;
          }
        }
      }
    if (IsFinished)
      break;
  }
  std::cerr << "RecursiveDualDesc |ListFace|=" << ListFace.size() << "\n";
  return ListFace;
}

template <typename T>
vectface
DUALDESC_SamplingFacetProcedure(MyMatrix<T> const &EXT,
                                std::vector<std::string> const &ListOpt) {
  std::string prog = "lrs";
  int critlevel = 50;
  int maxnbcall = -1;
  int maxnbsize = 20;
  for (auto &eOpt : ListOpt) {
    std::vector<std::string> ListStrB = STRING_Split(eOpt, "_");
    if (ListStrB.size() == 2) {
      if (ListStrB[0] == "prog")
        prog = ListStrB[1];
      if (ListStrB[0] == "critlevel")
        std::istringstream(ListStrB[1]) >> critlevel;
      if (ListStrB[0] == "maxnbcall")
        std::istringstream(ListStrB[1]) >> maxnbcall;
      if (ListStrB[0] == "maxnbsize")
        std::istringstream(ListStrB[1]) >> maxnbsize;
    }
  }
  if (prog != "lrs" && prog != "cdd") {
    std::cerr << "We have prog=" << prog << "\n";
    std::cerr << "but the only allowed input formats are lrs and cdd\n";
    throw TerminalException{1};
  }
  recSamplingOption eOption;
  eOption.maxnbcall = maxnbcall;
  eOption.prog = prog;
  eOption.critlevel = critlevel;
  eOption.maxnbsize = maxnbsize;
  int nbcall = 0;
  return Kernel_DUALDESC_SamplingFacetProcedure(EXT, eOption, nbcall);
}

template <typename T>
vectface
Kernel_DirectComputationInitialFacetSet(MyMatrix<T> const &EXT,
                                        std::string const &ansSamp,
                                        [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_SAMPLING_FACET
  MicrosecondTime time;
#endif
  int n_rows = EXT.rows();
  std::vector<std::string> ListStr = STRING_Split(ansSamp, ":");
  std::string ansOpt = ListStr[0];
  auto get_iter = [&]() -> int {
    int iter = 10;
    if (ListStr.size() > 1) {
      std::vector<std::string> ListStrB = STRING_Split(ListStr[1], "_");
      if (ListStrB.size() == 2 && ListStrB[0] == "iter")
        std::istringstream(ListStrB[1]) >> iter;
    }
    return iter;
  };
  auto get_face = [&]() -> Face {
    if (ListStr.size() != 2) {
      std::cerr << "We should have exactly two entries\n";
      throw TerminalException{1};
    }
    Face f(n_rows);
    std::string face_str = ListStr[1];
    for (int i_row = 0; i_row < n_rows; i_row++) {
      if (face_str.substr(i_row, 1) == "1") {
        f[i_row] = 1;
      }
    }
    return f;
  };
  auto compute_samp = [&]() -> vectface {
    if (ansOpt == "lp_cdd") {
      // So possible format is lp_cdd:iter_100
      int iter = get_iter();
      return FindVertices(EXT, iter);
    }
    if (ansOpt == "lp_cdd_min") {
      // So possible format is lp_cdd_min:iter_100
      int iter = get_iter();
      vectface vf = FindVertices(EXT, iter);
      return select_minimum_count(vf);
    }
    if (ansOpt == "sampling") {
      std::vector<std::string> ListOpt;
      int n_ent = ListStr.size();
      for (int i_ent = 1; i_ent < n_ent; i_ent++)
        ListOpt.push_back(ListStr[i_ent]);
      return DUALDESC_SamplingFacetProcedure(EXT, ListOpt);
    }
    if (ansOpt == "lrs_limited") {
      int upperlimit = 100;
      // So possible format is lrs_limited:upperlimit_1000
      if (ListStr.size() > 1) {
        std::vector<std::string> ListStrB = STRING_Split(ListStr[1], "_");
        if (ListStrB.size() == 2 && ListStrB[0] == "upperlimit")
          std::istringstream(ListStrB[1]) >> upperlimit;
      }
      return lrs::DualDescription_incd_limited(EXT, upperlimit);
    }
    if (ansOpt == "specific") {
      Face f = get_face();
      vectface vf(n_rows);
      vf.push_back(f);
      return vf;
    }
    std::cerr << "No right program found\n";
    std::cerr << "Let us die\n";
    throw TerminalException{1};
  };
  vectface ListIncd = compute_samp();
  if (ListIncd.size() == 0) {
    std::cerr << "We found 0 facet and that is not good\n";
    throw TerminalException{1};
  }
  std::map<size_t, size_t> map_incd;
  for (auto &eFace : ListIncd) {
    map_incd[eFace.count()] += 1;
  }
  os << "Found incidences =";
  for (auto &kv : map_incd) {
    os << " (" << kv.first << "," << kv.second << ")";
  }
  os << "\n";
#ifdef TIMINGS_SAMPLING_FACET
  os << "|DirectComputationInitialFacetSet|=" << time << "\n";
#endif
  return ListIncd;
}

template <typename T>
inline typename std::enable_if<is_ring_field<T>::value, vectface>::type
DirectComputationInitialFacetSet(MyMatrix<T> const &EXT,
                                 std::string const &ansSamp, std::ostream &os) {
  return Kernel_DirectComputationInitialFacetSet(EXT, ansSamp, os);
}

template <typename T>
inline typename std::enable_if<!is_ring_field<T>::value, vectface>::type
DirectComputationInitialFacetSet(MyMatrix<T> const &EXT,
                                 std::string const &ansSamp, std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> EXTfield = UniversalMatrixConversion<Tfield, T>(EXT);
  return Kernel_DirectComputationInitialFacetSet(EXTfield, ansSamp, os);
}

template <typename T>
vectface Kernel_GetFullRankFacetSet(
    const MyMatrix<T> &EXT,
    const MyMatrix<typename underlying_ring<T>::ring_type> &EXT_int,
    std::ostream &os) {
  using Tint = typename underlying_ring<T>::ring_type;
  size_t dim = EXT.cols();
  size_t n_rows = EXT.rows();
  if (dim == 2) {
    if (n_rows != 2) {
      std::cerr
          << "In dimension 2, the cone should have exactly two extreme rays\n";
      throw TerminalException{1};
    }
    vectface vf_ret(2);
    Face f1(2), f2(2);
    f1[0] = 1;
    f2[1] = 1;
    vf_ret.push_back(f1);
    vf_ret.push_back(f2);
    return vf_ret;
  }
  Face eSet = Kernel_FindSingleVertex(EXT);
  // Here we use a trick that the ColumnReduction will select the first column
  // and so will return a matrix that is polytopal
  MyMatrix<T> EXTsel_pre = SelectRow(EXT, eSet);
  MyMatrix<Tint> EXTsel_pre_int = SelectRow(EXT_int, eSet);
  std::vector<int> l_cols = ColumnReductionSet(EXTsel_pre);
  MyMatrix<T> EXTsel = SelectColumn(EXTsel_pre, l_cols);
  MyMatrix<Tint> EXTsel_int = SelectColumn(EXTsel_pre_int, l_cols);
  os << "|EXTsel|=" << EXTsel.rows() << " / " << EXTsel.cols()
     << " rnk=" << RankMat(EXTsel) << "\n";
#ifdef DEBUG_SAMPLING_FACET
  if (!IsPolytopal(EXTsel)) {
    std::cerr << "The configuration EXTsel is not polytopal\n";
    throw TerminalException{1};
  }
#endif
  vectface ListRidge = Kernel_GetFullRankFacetSet(EXTsel, EXTsel_int, os);
  os << "We have ListRidge\n";
  FlippingFramework<T> RPLlift(EXT, EXT_int, eSet);
  os << "We have FlippingFramework\n";
  vectface vf_ret(n_rows);
  vf_ret.push_back(eSet);
  for (auto &eRidge : ListRidge) {
    Face eFace = RPLlift.FlipFace(eRidge);
    vf_ret.push_back(eFace);
  }
  os << "We have vf_ret\n";
#ifdef DEBUG_SAMPLING_FACET
  MyMatrix<T> FACsamp(vf_ret.size(), dim);
  int pos = 0;
  for (auto &face : vf_ret) {
    MyVector<T> eFAC = FindFacetInequality(EXT, face);
    AssignMatrixRow(FACsamp, pos, eFAC);
    pos++;
  }
  if (RankMat(FACsamp) != static_cast<int>(dim)) {
    std::cerr << "Failed to find a ful rank vector configuration\n";
    throw TerminalException{1};
  }
#endif
  return vf_ret;
}

template <typename T>
vectface GetFullRankFacetSet(const MyMatrix<T> &EXT, std::ostream &os) {
#ifdef TIMINGS_SAMPLING_FACET
  MicrosecondTime time;
#endif
  MyMatrix<T> EXT_B = ColumnReduction(EXT);
#ifdef TIMINGS_SAMPLING_FACET
  os << "|ColumnReduction|=" << time << "\n";
#endif
  MyMatrix<T> EXT_C = Polytopization(EXT_B);
#ifdef TIMINGS_SAMPLING_FACET
  os << "|Polytopization|=" << time << "\n";
#endif
  MyMatrix<T> EXT_D = SetIsobarycenter(EXT_C);
#ifdef TIMINGS_SAMPLING_FACET
  os << "|SetIsobarycenter|=" << time << "\n";
#endif
  using Tint = typename underlying_ring<T>::ring_type;
  MyMatrix<Tint> EXT_D_int = Get_EXT_int(EXT_D);
  vectface vf = Kernel_GetFullRankFacetSet(EXT_D, EXT_D_int, os);
#ifdef TIMINGS_SAMPLING_FACET
  os << "|Kernel_GetFullRankFacetSet|=" << time << "\n";
#endif
  return vf;
}

// clang-format off
#endif  // SRC_POLY_POLY_SAMPLINGFACET_H_
// clang-format on

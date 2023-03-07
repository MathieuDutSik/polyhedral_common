// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_SAMPLINGFACET_H_
#define SRC_POLY_POLY_SAMPLINGFACET_H_

#include "POLY_LinearProgramming.h"
#include "POLY_PolytopeFct.h"
#include "POLY_lrslib.h"
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

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
vectface DirectComputationInitialFacetSet(MyMatrix<T> const &EXT,
                                          std::string const &ansSamp,
                                          std::ostream &os) {
  os << "DirectComputationInitialFacetSet ansSamp=" << ansSamp << "\n";
  std::vector<std::string> ListStr = STRING_Split(ansSamp, ":");
  std::string ansOpt = ListStr[0];
  auto compute_samp = [&]() -> vectface {
    if (ansOpt == "lp_cdd") {
      // So possible format is lp_cdd:iter_100
      int iter = 10;
      if (ListStr.size() > 1) {
        std::vector<std::string> ListStrB = STRING_Split(ListStr[1], "_");
        if (ListStrB.size() == 2 && ListStrB[0] == "iter")
          std::istringstream(ListStrB[1]) >> iter;
      }
      return FindVertices(EXT, iter);
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
    std::cerr << "No right program found\n";
    std::cerr << "Let us die\n";
    throw TerminalException{1};
  };
  vectface ListIncd = compute_samp();
  if (ListIncd.size() == 0) {
    std::cerr << "We found 0 facet and that is not good\n";
    throw TerminalException{1};
  }
  return ListIncd;
}

template <typename T> vectface GetFullRankFacetSet(const MyMatrix<T> &EXT) {
  // Heuristic first, should work in many cases
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  size_t dim = EXTred.cols();
  size_t n_rows = EXT.rows();
  size_t nb = 10 * dim;
  vectface ListSets = FindVertices(EXTred, nb);
  std::unordered_set<Face> set_face;
  for (auto &eFace : ListSets)
    set_face.insert(eFace);
  size_t n_face = set_face.size();
  MyMatrix<T> FAC(n_face, dim);
  size_t i_face = 0;
  for (auto &eFace : set_face) {
    MyVector<T> V = FindFacetInequality(EXT, eFace);
    AssignMatrixRow(FAC, i_face, V);
    i_face++;
  }
  size_t rnk = RankMat(FAC);
  std::cerr << "GetFullRankFacetSet |FAC|=" << FAC.rows() << " |EXT|=" << n_rows
            << " rnk=" << rnk << " dim=" << dim << "\n";
  if (rnk == dim) {
    vectface vf(n_rows);
    for (auto &eFace : set_face)
      vf.push_back(eFace);
    return vf;
  }
  // Otherwise we call recursively
  std::cerr << "Failing, so calling the recursive algo\n";
  auto get_minincd = [&]() -> Face {
    size_t min_incd = std::numeric_limits<size_t>::max();
    for (auto &eFace : set_face) {
      size_t incd = eFace.count();
      if (incd < min_incd)
        min_incd = incd;
    }
    for (auto &eFace : set_face) {
      size_t incd = eFace.count();
      if (incd == min_incd)
        return eFace;
    }
    std::cerr << "We should not reach that stage\n";
    throw TerminalException{1};
  };
  Face eSet = get_minincd();
  MyMatrix<T> EXTsel = ColumnReduction(SelectRow(EXTred, eSet));
  std::cerr << "|EXTsel|=" << EXTsel.rows() << " / " << EXTsel.cols()
            << " rnk=" << RankMat(EXTsel) << "\n";
  vectface ListRidge = GetFullRankFacetSet(EXTsel);
  std::cerr << "We have ListRidge\n";
  FlippingFramework<T> RPLlift(EXTred, eSet);
  std::cerr << "We have FlippingFramework\n";
  vectface vf_ret(n_rows);
  vf_ret.push_back(eSet);
  for (auto &eRidge : ListRidge) {
    Face eFace = RPLlift.Flip(eRidge);
    vf_ret.push_back(eFace);
  }
  std::cerr << "We have vf_ret\n";
  return vf_ret;
}

// clang-format off
#endif  // SRC_POLY_POLY_SAMPLINGFACET_H_
// clang-format on

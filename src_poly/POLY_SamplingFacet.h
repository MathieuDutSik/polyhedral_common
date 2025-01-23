// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_SAMPLINGFACET_H_
#define SRC_POLY_POLY_SAMPLINGFACET_H_

// clang-format off
#include "POLY_LinearProgramming.h"
#include "POLY_Fundamental.h"
#include "POLY_DirectDualDesc.h"
#include <map>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_SAMPLING_FACET
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_SAMPLING_FACET
#endif

#ifdef TIMINGS
#define TIMINGS_SAMPLING_FACET
#endif

struct recSamplingOption {
  int critlevel;
  int critdim;
  int maxnboper;
  int maxnbsize;
  bool greedy_termination;
  std::string prog;
};

template <typename T>
vectface
Kernel_DUALDESC_SamplingFacetProcedure(MyMatrix<T> const &EXT,
                                       recSamplingOption const &eOption,
                                       int &nbOper, std::ostream &os) {
  int dim = RankMat(EXT);
  int len = EXT.rows();
  std::string prog = eOption.prog;
  int critlevel = eOption.critlevel;
  int critdim = eOption.critdim;
  int maxnboper = eOption.maxnboper;
  int maxnbsize = eOption.maxnbsize;
  bool greedy_termination = eOption.greedy_termination;
#ifdef DEBUG_SAMPLING_FACET
  os << "SAMP: critlevel=" << critlevel << " critdim=" << critdim
     << " prog=" << prog << " maxnboper=" << maxnboper
     << " maxnbsize=" << maxnbsize << "\n";
#endif
  auto IsRecursive = [&]() -> bool {
    if (len < critlevel) {
      return false;
    }
    if (dim < critdim) {
      return false;
    }
    return true;
  };
  bool DoRecur = IsRecursive();
#ifdef DEBUG_SAMPLING_FACET
  os << "SAMP: dim=" << dim << "  len=" << len << " DoRecur=" << DoRecur << " nbOper=" << nbOper << "\n";
#endif
  vectface ListFace(EXT.rows());
  std::unordered_set<size_t> set_incd;
  auto FuncInsert = [&](Face const &eFace) -> void {
#ifdef DEBUG_SAMPLING_FACET_DISABLE
    os << "SAMP: eFace.count()=" << eFace.count() << "\n";
    os << "SAMP: ListFace=";
    for (auto &fFace : ListFace) {
      os << fFace.count() << " ";
    }
    os << "\n";
#endif
    size_t incd = eFace.count();
    if (set_incd.count(incd) == 1) {
      return;
    }
    ListFace.push_back(eFace);
    set_incd.insert(incd);
  };
  if (!DoRecur) {
    auto comp_dd = [&]() -> vectface {
      if (is_method_supported<T>(prog)) {
        return DirectFacetComputationIncidence(EXT, prog, os);
      }
      if (prog == "fullrankfacetset") {
        return GetFullRankFacetSet(EXT, os);
      }
      std::cerr << "SAMP: Failed to find a matching method\n";
      throw TerminalException{1};
    };
#ifdef TIMINGS_SAMPLING_FACET
    MicrosecondTime time;
#endif
    vectface ListIncd = comp_dd();
#ifdef TIMINGS_SAMPLING_FACET
    os << "|SAMP: ListIncd, comp_dd|=" << time << "\n";
#endif
    for (auto &eFace : ListIncd) {
      FuncInsert(eFace);
    }
#ifdef DEBUG_SAMPLING_FACET
    os << "SAMP: DirectDualDesc |ListFace|=" << ListFace.size() << "\n";
#endif
    nbOper++;
    return ListFace;
  }
  Face eInc = FindOneInitialVertex(EXT, os);
  FuncInsert(eInc);
  size_t start = 0;
  while (true) {
    size_t nbCases = ListFace.size();
#ifdef DEBUG_SAMPLING_FACET
    os << "SAMP: while loop start=" << start << " nbCases=" << nbCases << "\n";
#endif
    for (size_t iC=start; iC<nbCases; iC++) {
      size_t nbCaseBefore = ListFace.size();
      nbOper++;
      Face eFace = ListFace[iC];
#ifdef DEBUG_SAMPLING_FACET_DISABLE
      os << "SAMP: len=" << len << " dim=" << dim << " treating iC=" << iC << " |eFace|=" << eFace.size() << "/" << eFace.count() << "\n";
#endif
      MyMatrix<T> EXTred = SelectRow(EXT, eFace);
      vectface ListRidge =
        Kernel_DUALDESC_SamplingFacetProcedure(EXTred, eOption, nbOper, os);
      SimplifiedFlippingFramework<T> sff(EXT, eFace, os);
      for (auto &eRidge : ListRidge) {
        Face eFlip = sff.FlipFace(eRidge);
        FuncInsert(eFlip);
      }
      size_t nbCaseAfter = ListFace.size();
      // We have enough, ending the stuff.
      if (maxnbsize != -1) {
        int siz = ListFace.size();
        if (siz > maxnbsize) {
#ifdef DEBUG_SAMPLING_FACET
          os << "SAMP: Ending by maxsize criterion\n";
          os << "SAMP: siz=" << siz << " maxnbsize=" << maxnbsize << "\n";
#endif
          return ListFace;
        }
      }
      // Too many calls, we drop out
      if (maxnboper != -1) {
        if (nbOper > maxnboper) {
#ifdef DEBUG_SAMPLING_FACET
          os << "SAMP: Ending by maxnboper\n";
#endif
          return ListFace;
        }
      }
      // Aggressive termination if selected, we 
      if (greedy_termination) {
        if (nbCaseBefore == nbCaseAfter) {
#ifdef DEBUG_SAMPLING_FACET
          os << "SAMP: Ending due to lack of progress\n";
#endif
          return ListFace;
        }
      }
    }
    start = nbCases;
    // All orbits have been treated.
    if (nbCases == ListFace.size()) {
#ifdef DEBUG_SAMPLING_FACET
      os << "SAMP: Ending since all orbits are treated\n";
#endif
      break;
    }
  }
#ifdef DEBUG_SAMPLING_FACET
  os << "SAMP: RecursiveDualDesc |ListFace|=" << ListFace.size() << "\n";
#endif
  return ListFace;
}

template <typename T>
vectface
DUALDESC_SamplingFacetProcedure(MyMatrix<T> const &EXT,
                                std::vector<std::string> const &ListOpt,
                                std::ostream &os) {
  std::string prog = "lrs";
  int critlevel = 50;
  int critdim = 15;
  int maxnboper = 1000;
  int maxnbsize = 20;
  bool greedy_termination = true;
  for (auto &eOpt : ListOpt) {
    std::vector<std::string> ListStrB = STRING_Split(eOpt, "_");
    if (ListStrB.size() == 2) {
      if (ListStrB[0] == "prog") {
        prog = ListStrB[1];
      }
      if (ListStrB[0] == "critlevel") {
        std::istringstream(ListStrB[1]) >> critlevel;
      }
      if (ListStrB[0] == "critdim") {
        std::istringstream(ListStrB[1]) >> critdim;
      }
      if (ListStrB[0] == "maxnboper") {
        std::istringstream(ListStrB[1]) >> maxnboper;
      }
      if (ListStrB[0] == "maxnbsize") {
        std::istringstream(ListStrB[1]) >> maxnbsize;
      }
      if (ListStrB[0] == "greedy_termination") {
        if (ListStrB[1] == "true") {
          greedy_termination = true;
        }
        if (ListStrB[1] == "false") {
          greedy_termination = false;
        }
      }
    }
  }
  if (prog != "lrs" && prog != "cdd") {
    std::cerr << "SAMP: We have prog=" << prog << "\n";
    std::cerr << "SAMP: but the only allowed input formats are lrs and cdd\n";
    throw TerminalException{1};
  }
  recSamplingOption eOption;
  eOption.maxnboper = maxnboper;
  eOption.prog = prog;
  eOption.critlevel = critlevel;
  eOption.critdim = critdim;
  eOption.maxnbsize = maxnbsize;
  eOption.greedy_termination = greedy_termination;
  int nbOper = 0;
  return Kernel_DUALDESC_SamplingFacetProcedure(EXT, eOption, nbOper, os);
}

template <typename T>
vectface Kernel_DirectComputationInitialFacetSet(MyMatrix<T> const &EXT,
                                                 std::string const &ansSamp,
                                                 std::ostream &os) {
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
      std::cerr << "SAMP: We should have exactly two entries\n";
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
      return FindVertices(EXT, iter, os);
    }
    if (ansOpt == "lp_cdd_min") {
      // So possible format is lp_cdd_min:iter_100
      int iter = get_iter();
      vectface vf = FindVertices(EXT, iter, os);
      return select_minimum_count(vf);
    }
    if (ansOpt == "sampling") {
      std::vector<std::string> ListOpt;
      int n_ent = ListStr.size();
      for (int i_ent = 1; i_ent < n_ent; i_ent++)
        ListOpt.push_back(ListStr[i_ent]);
      return DUALDESC_SamplingFacetProcedure(EXT, ListOpt, os);
    }
    if (ansOpt == "fullrankfacetset") {
      return GetFullRankFacetSet(EXT, os);
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
    std::cerr << "SAMP: No right program found\n";
    std::cerr << "SAMP: Let us die\n";
    throw TerminalException{1};
  };
  vectface ListIncd = compute_samp();
  if (ListIncd.size() == 0) {
    std::cerr << "SAMP: We found 0 facet and that is not good\n";
    throw TerminalException{1};
  }
  std::map<size_t, size_t> map_incd;
  for (auto &eFace : ListIncd) {
    map_incd[eFace.count()] += 1;
  }
#ifdef DEBUG_SAMPLING_FACET
  os << "SAMP: Found incidences =";
  for (auto &kv : map_incd) {
    os << " (" << kv.first << "," << kv.second << ")";
  }
  os << "\n";
#endif
#ifdef TIMINGS_SAMPLING_FACET
  os << "|SAMP: DirectComputationInitialFacetSet|=" << time << "\n";
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

// clang-format off
#endif  // SRC_POLY_POLY_SAMPLINGFACET_H_
// clang-format on

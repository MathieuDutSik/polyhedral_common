// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_WEIGHTMATRIXSPECIFIED_H_
#define SRC_GROUP_WEIGHTMATRIXSPECIFIED_H_

/*
  The two main applications of weight matrices are:
  ---Computing the canonical ordering of the vertices.
  ---Computing the stabilizer of weight matrix.

  This can be slow at times so we need to apply a
  number of tricks to accelerate. In this section
  we apply the idea of computing on a subset of
  the vertices and then checking from that that
  everything is correct.

  The check of correctness is done in the following way:
  ---For the stabilizer we check that the obtained
  automorphism are actually ok.
  ---For the canonical form, we also compute the
  automorphism group. If the automorphism happens to be
  correct, then this implies that the canonical form
  is actually also ok.

  The code uses a number of template functions provided
  as input: f1, f2, f3, f4, f5.
 */

// The hash map do not seem to make much difference in the overall
// performance.

// #define UNORDERED_MAP_SPECIFIC
// #define TSL_SPARSE_MAP_SPECIFIC
// #define TSL_ROBIN_MAP_SPECIFIC
#define TSL_HOPSCOTCH_MAP_SPECIFIC

#ifdef UNORDERED_MAP_SPECIFIC
#include <unordered_map>
#include <unordered_set>
#define UNORD_MAP_SPECIFIC std::unordered_map
#define UNORD_SET_SPECIFIC std::unordered_set
#endif

#ifdef TSL_SPARSE_MAP_SPECIFIC
#include "sparse_map.h"
#include "sparse_set.h"
#define UNORD_MAP_SPECIFIC tsl::sparse_map
#define UNORD_SET_SPECIFIC tsl::sparse_set
#endif

#ifdef TSL_ROBIN_MAP_SPECIFIC
#include "robin_map.h"
#include "robin_set.h"
#define UNORD_MAP_SPECIFIC tsl::robin_map
#define UNORD_SET_SPECIFIC tsl::robin_set
#endif

#ifdef TSL_HOPSCOTCH_MAP_SPECIFIC
#include "hopscotch_map.h"
#include "hopscotch_set.h"
#define UNORD_MAP_SPECIFIC tsl::hopscotch_map
#define UNORD_SET_SPECIFIC tsl::hopscotch_set
#endif

// clang-format off
#include "GRAPH_Bindings.h"
#include "WeightMatrix.h"
#include <algorithm>
#include <limits>
#include <map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_WEIGHT_MATRIX_SPECIFIED
#define DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
#endif

#ifdef DISABLE_DEBUG_WEIGHT_MATRIX_SPECIFIED
#undef DEBUG_WEIGHT_MATRIX_SPECIFIED
#endif

#ifdef DISABLE_DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
#undef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
#endif

#ifdef TIMINGS
#define TIMINGS_WEIGHT_MATRIX_SPECIFIED
#endif

template <typename T>
std::pair<std::vector<T>, std::vector<int>>
GetReorderingInfoWeight(const std::vector<T> &ListWeight) {
  size_t n_Wei = ListWeight.size();
  std::map<T, size_t> map;
  for (size_t i = 0; i < n_Wei; i++)
    map[ListWeight[i]] = i;
  std::vector<int> g(n_Wei);
  int idx = 0;
  for (auto &kv : map) {
    int pos = kv.second;
    g[pos] = idx;
    idx++;
  }
  std::vector<T> NewListWeight(n_Wei);
  for (size_t iW = 0; iW < n_Wei; iW++) {
    size_t jW = g[iW];
    NewListWeight[jW] = ListWeight[iW];
  }
  return {std::move(NewListWeight), std::move(g)};
}

/*
  Data structure for the partitionning of the vertex set.
  --- MapVertexBlock maps the vertices to the blocks.
  --- The ListBlocks is the corresponding data structure that
  is the list of blocks (has to be coherent with MapVertexBlock
  of course).
 */
template <typename Tidx> struct VertexPartition {
  size_t nbRow;
  std::vector<Tidx> MapVertexBlock;
  std::vector<std::vector<Tidx>> ListBlocks;
};

template <typename T, typename Tidx, typename F1, typename F2>
VertexPartition<Tidx> ComputeInitialVertexPartition(size_t nbRow, F1 f1, F2 f2,
                                                    bool canonically,
                                                    [[maybe_unused]] std::ostream &os) {
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "ComputeInitialVertexPartition : We have nbRow=" << nbRow
              << " which is larger than the possible values of Tidx : "
              << max_poss_rows << "\n";
    throw TerminalException{1};
  }
  std::vector<T> ListWeight;
  UNORD_MAP_SPECIFIC<T, size_t> ValueMap_T;
  size_t idxWeight = 0;
  auto get_T_idx = [&](T eval) -> Tidx {
    size_t &idx = ValueMap_T[eval];
    if (idx == 0) {
      // value is missing
      idxWeight++;
      idx = idxWeight;
      ListWeight.push_back(eval);
    }
    return static_cast<Tidx>(idx - 1);
  };
  std::vector<Tidx> MapVertexBlock(nbRow);
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    f1(iRow);
    T val = f2(iRow);
    MapVertexBlock[iRow] = get_T_idx(val);
  }
#ifdef SANITY_CHECK_WEIGHT_MATRIX_SPECIFIED
  if (idxWeight != ListWeight.size()) {
    std::cerr << "WMS: Incoherent lengths\n";
    throw TerminalException{1};
  }
#endif
  if (canonically) {
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
    os << "WMS: canonically reordering\n";
#endif
    std::pair<std::vector<T>, std::vector<int>> rec_pair =
        GetReorderingInfoWeight(ListWeight);
    const std::vector<int> &g = rec_pair.second;
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
    os << "WMS: |g|=" << g.size() << "\n";
    os << "WMS: |MapVertexBlock|=" << MapVertexBlock.size() << " nbRow=" << nbRow << "\n";
    os << "WMS: ComputeInitialVertexPartition, ListWeight=" << rec_pair.first << "\n";
#endif
    for (size_t iRow = 0; iRow < nbRow; iRow++) {
      int NewIdx = g[MapVertexBlock[iRow]];
      MapVertexBlock[iRow] = NewIdx;
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
      os << "WMS:   iRow=" << iRow << " MapVertexBlock[iRow]=" << MapVertexBlock[iRow] << " NewIdx=" << NewIdx << "\n";
#endif
    }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
    const std::vector<T> & NewListWeight = rec_pair.first;
    for (size_t idx=1; idx<idxWeight; idx++) {
      size_t idx1 = idx - 1;
      size_t idx2 = idx;
      if (NewListWeight[idx1] > NewListWeight[idx2]) {
        std::cerr << "WMS: idx1=" << idx1 << " weight1=" << NewListWeight[idx1] << "\n";
        std::cerr << "WMS: idx2=" << idx2 << " weight2=" << NewListWeight[idx1] << "\n";
        throw TerminalException{1};
      }
    }
#endif
  }
  size_t n_weight = ListWeight.size();
  std::vector<std::vector<Tidx>> ListBlocks(n_weight);
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    Tidx pos = MapVertexBlock[iRow];
    ListBlocks[pos].push_back(iRow);
  }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
  os << "WMS: n_weight=" << n_weight << "\n";
  for (size_t i_weight=0; i_weight<n_weight; i_weight++) {
    os << "WMS: i_weight=" << i_weight << " |Block|=" << ListBlocks[i_weight].size() << "\n";
  }
#endif
  return {nbRow, std::move(MapVertexBlock), std::move(ListBlocks)};
}

// jBlock : is the block for which we look at breaking down
template <typename T, typename Tidx, typename F1, typename F2>
bool RefineSpecificVertexPartition(VertexPartition<Tidx> &VP, const int &jBlock,
                                   const int &iBlock, F1 f1, F2 f2,
                                   bool canonically,
                                   [[maybe_unused]] std::ostream &os) {
  size_t nbRow = VP.nbRow;
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "RefineSpecificVertexPartition : We have nbRow=" << nbRow
              << " which is larger than the possible values of Tidx : "
              << max_poss_rows << "\n";
    throw TerminalException{1};
  }
  // The code for finding the indices
  UNORD_MAP_SPECIFIC<T, int> ValueMap_T;
  std::vector<T> ListWeight;
  int idxWeight = 0;
  auto get_T_idx = [&](const T &eval) -> int {
    int &idx = ValueMap_T[eval];
    if (idx == 0) {
      // value is missing
      idxWeight++;
      idx = idxWeight;
      ListWeight.push_back(eval);
    }
    return idx - 1;
  };
  UNORD_MAP_SPECIFIC<std::vector<int>, int> ValueMap_Tvs;
  std::vector<std::vector<int>> ListPossibleSignatures;
  int idxSign = 0;
  auto get_Tvs_idx = [&](std::vector<int> const &esign) -> int {
    int &idx = ValueMap_Tvs[esign];
    if (idx == 0) {
      // value is missing
      idxSign++;
      idx = idxSign;
      ListPossibleSignatures.push_back(esign);
    }
    return idx - 1;
  };
  // Now computing the indices
  // We have to do a copy operation for eBlockRead
  std::vector<Tidx> eBlockBreak = VP.ListBlocks[jBlock];
  const std::vector<Tidx> &eBlockSpec = VP.ListBlocks[iBlock];
  size_t siz_block_break = eBlockBreak.size();
  size_t siz_block_spec = eBlockSpec.size();
  std::vector<int> V(siz_block_break);
  int len = 0;
  std::vector<int> list_mult;
  std::vector<int> esign;
  for (size_t iBreak = 0; iBreak < siz_block_break; iBreak++) {
    size_t iRow = eBlockBreak[iBreak];
    f1(iRow);
    for (size_t iSpec = 0; iSpec < siz_block_spec; iSpec++) {
      size_t jRow = eBlockSpec[iSpec];
      T val = f2(jRow);
      int idx = get_T_idx(val);
      if (idx == len) {
        list_mult.push_back(0);
        len++;
      }
      list_mult[idx]++;
    }
    for (int u = 0; u < len; u++) {
      if (list_mult[u] > 0) {
        esign.push_back(u);
        esign.push_back(list_mult[u]);
        list_mult[u] = 0;
      }
    }
    int idx_sign = get_Tvs_idx(esign);
    esign.clear();
    V[iBreak] = idx_sign;
  }
  size_t n_block = ListPossibleSignatures.size();
  if (n_block == 1) {
    return false;
  }
  if (canonically) {
    // First reordering the weights
    std::pair<std::vector<T>, std::vector<int>> rec_pair_A =
        GetReorderingInfoWeight(ListWeight);
    ListWeight = rec_pair_A.first;
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
    os << "WMS: RefineSpecificVertexPartition, ListWeight=" << ListWeight << "\n";
#endif
    size_t n_Wei = ListWeight.size();
    std::vector<int> list_mult(n_Wei, 0);
    const std::vector<int> &g = rec_pair_A.second;
    for (auto &esign : ListPossibleSignatures) {
      size_t n_ent = esign.size() / 2;
      for (size_t i = 0; i < n_ent; i++) {
        int NewVal = g[esign[2 * i]];
        int eMult = esign[1 + 2 * i];
        list_mult[NewVal] = eMult;
      }
      //
      std::vector<int> newsign;
      for (size_t i = 0; i < n_Wei; i++) {
        if (list_mult[i] > 0) {
          newsign.push_back(i);
          newsign.push_back(list_mult[i]);
          list_mult[i] = 0;
        }
      }
      esign = std::move(newsign);
    }
    // Now reordering the signatures
    std::pair<std::vector<std::vector<int>>, std::vector<int>> rec_pair_B =
        GetReorderingInfoWeight(ListPossibleSignatures);
    ListPossibleSignatures = rec_pair_B.first;
    const std::vector<int> &h = rec_pair_B.second;
    for (size_t iBreak = 0; iBreak < siz_block_break; iBreak++) {
      int NewIdx = h[V[iBreak]];
      V[iBreak] = NewIdx;
    }
  }
  int n_block_orig = VP.ListBlocks.size();
  std::vector<int> Vmap(n_block);
  Vmap[0] = jBlock;
  for (size_t i = 1; i < n_block; i++) {
    int pos = n_block_orig + i - 1;
    Vmap[i] = pos;
  }
  VP.ListBlocks[jBlock].clear();
  for (size_t i = 1; i < n_block; i++) {
    VP.ListBlocks.push_back(std::vector<Tidx>());
  }
  for (size_t iBreak = 0; iBreak < siz_block_break; iBreak++) {
    int iRow = eBlockBreak[iBreak];
    int iBlockLoc = V[iBreak];
    int iBlockGlob = Vmap[iBlockLoc];
    VP.MapVertexBlock[iRow] = iBlockGlob;
    VP.ListBlocks[iBlockGlob].push_back(iRow);
  }
  return true;
}

#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
template <typename Tidx>
void PrintVertexPartitionInfo(const VertexPartition<Tidx> &VP,
                              const std::vector<uint8_t> &status,
                              std::ostream &os) {
  size_t n_block = VP.ListBlocks.size();
  os << "WMS: VP, nbRow=" << VP.nbRow << "  ListBlocks=";
  for (size_t i = 0; i < n_block; i++) {
    os << VP.ListBlocks[i].size() << "," << int(status[i]) << " ";
  }
  os << "\n";
}
#endif

template <typename T, typename Tidx, typename F1, typename F2>
VertexPartition<Tidx>
ComputeVertexPartition(size_t nbRow, F1 f1, F2 f2, bool canonically,
                       size_t max_globiter, std::ostream &os) {
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "ComputeVertexPartition : We have nbRow=" << nbRow
              << " which is larger than the possible values of Tidx : "
              << max_poss_rows << "\n";
    throw TerminalException{1};
  }
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  MicrosecondTime time;
#endif
  VertexPartition<Tidx> VP =
      ComputeInitialVertexPartition<T, Tidx>(nbRow, f1, f2, canonically, os);
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  os << "|WMS: ComputeInitialVertexPartition|=" << time << "\n";
#endif
  std::vector<uint8_t> status(VP.ListBlocks.size(), 0);
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
  os << "WMS: After ComputeInitialVertexPartition\n";
  PrintVertexPartitionInfo(VP, status, os);
#endif
  auto GetPreferable_iBlock = [&]() -> int {
    size_t min_size = nbRow + 1;
    int iBlockSel = -1;
    int n_block = VP.ListBlocks.size();
    for (int iBlock = 0; iBlock < n_block; iBlock++) {
      size_t siz = VP.ListBlocks[iBlock].size();
      if (status[iBlock] == 0 && siz < max_globiter) {
        if (iBlockSel == -1) {
          iBlockSel = iBlock;
          min_size = siz;
        } else {
          if (siz < min_size) {
            iBlockSel = iBlock;
            min_size = siz;
          }
        }
      }
    }
    return iBlockSel;
  };
  auto DoRefinement = [&](const int &iBlock) -> bool {
    int n_block = VP.ListBlocks.size();
    // First looking at the diagonal
    bool test1 = RefineSpecificVertexPartition<T, Tidx>(VP, iBlock, iBlock, f1,
                                                        f2, canonically, os);
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
    os << "WMS:   iBlock/iBlock=" << iBlock << " test1=" << test1 << "\n";
#endif
    if (test1) {
      size_t len1 = VP.ListBlocks.size();
      size_t len2 = status.size();
      for (size_t i = len2; i < len1; i++) {
        status.push_back(0);
      }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
      os << "WMS:   DoRefinement early, len1=" << len1 << " len2=" << len2 << "\n";
#endif
      return true;
    }
    status[iBlock] = 1;
    bool DidSomething = false;
    for (int jBlock = 0; jBlock < n_block; jBlock++) {
      if (iBlock != jBlock) {
        bool test2 = RefineSpecificVertexPartition<T, Tidx>(VP, jBlock, iBlock, f1, f2, canonically, os);
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
        os << "WMS:   iBlock=" << iBlock << " jBlock=" << jBlock
           << " test2=" << test2
           << " len1=" << VP.ListBlocks.size() << " len2=" << status.size() << "\n";
#endif
        if (test2) {
          status[jBlock] = 0;
          size_t len1 = VP.ListBlocks.size();
          size_t len2 = status.size();
          for (size_t i = len2; i < len1; i++) {
            status.push_back(0);
          }
          DidSomething = true;
        }
      }
    }
    return DidSomething;
  };
  while (true) {
    int iBlock = GetPreferable_iBlock();
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
    os << "WMS: GetPreferable_iBlock, iBlock=" << iBlock << "\n";
#endif
    if (iBlock == -1) {
      break;
    }
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
    MicrosecondTime time;
#endif
    bool test = DoRefinement(iBlock);
#ifdef TIMINGS
    os << "|WMS: DoRefinement|=" << time << "\n";
#endif
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
    os << "WMS: After Dorefinement\n";
    PrintVertexPartitionInfo(VP, status, os);
    os << "WMS: DoRefinement, test=" << test << "\n";
#endif
    if (!test) {
      break;
    }
  }
  return VP;
}

// We use stable_sort to make sure that from the canonical ordering of the
// blocks we get a canonical ordering for the block sizes.
template <typename Tidx>
std::vector<size_t> GetOrdering_ListIdx(const VertexPartition<Tidx> &VP) {
  size_t nbCase = VP.ListBlocks.size();
  std::vector<size_t> ListIdx(nbCase);
  for (size_t iCase = 0; iCase < nbCase; iCase++)
    ListIdx[iCase] = iCase;
  std::stable_sort(
      ListIdx.begin(), ListIdx.end(), [&](int idx1, int idx2) -> bool {
        return VP.ListBlocks[idx1].size() < VP.ListBlocks[idx2].size();
      });
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
  for (size_t iCase=1; iCase<nbCase; iCase++) {
    size_t idx1 = ListIdx[iCase-1];
    size_t idx2 = ListIdx[iCase];
    size_t siz1 = VP.ListBlocks[idx1].size();
    size_t siz2 = VP.ListBlocks[idx2].size();
    if (siz1 > siz2) {
      std::cerr << "WMS: iCase=" << iCase << " siz1=" << siz1 << " siz2=" << siz2 << "\n";
      throw TerminalException{1};
    }
  }
#endif
  return ListIdx;
}

//
// The computation of vertex signatures is needed for the computation of vertex
// degrees which are needed for TRACES graph construction. Therefore said
// computation cannot be avoided when building our objects. However, for
// breaking down the vertex set into blocks, they are too expensive
//
// Explanations:
// ---nbRow: The number of vertices of the corresponding graph corresponding
//    graph.
// ---ListWeight: The list of weights occurring in that specified subset.
// ---ListPossibleSignatures: The list of possible signatures.
//    The signature is the list of possible values.
// ---ListSignatureByVertex: The signature by the vertices.
// ---ListNbCase: For each possible signature, the number of matching entries.
template <typename T> struct WeightMatrixVertexSignatures {
  size_t nbRow;
  size_t nbWeight;
  std::vector<T> ListWeight;
  std::vector<std::vector<int>> ListPossibleSignatures;
  std::vector<int> ListSignatureByVertex;
  std::vector<int> ListNbCase;
};

template<typename T>
WeightMatrixVertexSignatures<T> EmptyWeightMatrixVertexSignatures() {
  size_t nbRow = 0;
  size_t nbWeight = 0;
  std::vector<T> ListWeight;
  std::vector<std::vector<int>> ListPossibleSignatures;
  std::vector<int> ListSignatureByVertex;
  std::vector<int> ListNbCase;
  return {nbRow,
          nbWeight,
          std::move(ListWeight),
          std::move(ListPossibleSignatures),
          std::move(ListSignatureByVertex),
          std::move(ListNbCase)};
}

template <typename T> struct PairWeightMatrixVertexSignatures {
  WeightMatrixVertexSignatures<T> WMVS_direct;
  WeightMatrixVertexSignatures<T> WMVS_dual;
};

/*
template <typename T>
std::vector<int>
GetOrdering_ListIdx(WeightMatrixVertexSignatures<T> const &WMVS,
                    [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  MicrosecondTime time;
#endif
  size_t nbCase = WMVS.ListNbCase.size();
  std::vector<int> ListIdx(nbCase);
  for (size_t iCase = 0; iCase < nbCase; iCase++)
    ListIdx[iCase] = iCase;
  auto fctComp = [](auto val1, auto val2) -> std::pair<bool, bool> {
    if (val1 < val2)
      return {true, true};
    if (val1 > val2)
      return {true, false};
    return {false, true};
  };
  std::sort(ListIdx.begin(), ListIdx.end(), [&](int idx1, int idx2) -> bool {
    // First selection by the number of cases.
    std::pair<bool, bool> test1 =
        fctComp(WMVS.ListNbCase[idx1], WMVS.ListNbCase[idx2]);
    if (test1.first) {
      return test1.second;
    }
    // The cases with high number of cases are preferable.
    size_t len1 = WMVS.ListPossibleSignatures[idx1].size();
    size_t len2 = WMVS.ListPossibleSignatures[idx2].size();
    std::pair<bool, bool> test2 = fctComp(len2, len1);
    if (test2.first)
      return test2.second;
    // Now the order does not really matter for speed but it has to be
    // fully unicized for the canonical form.
    // Now going after the diagonal value
    const std::vector<int> &list_pair1 = WMVS.ListPossibleSignatures[idx1];
    const std::vector<int> &list_pair2 = WMVS.ListPossibleSignatures[idx2];
    for (size_t i = 0; i < len1; i++) {
      std::pair<bool, bool> test4 = fctComp(list_pair1[i], list_pair2[i]);
      if (test4.first) {
        return test4.second;
      }
    }
    return false;
  });
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  os << "|WMS: GetOrdering_ListIdx|=" << time << "\n";
#endif
  return ListIdx;
}
*/

template <typename T>
void PrintWMVS(WeightMatrixVertexSignatures<T> const &WMVS, std::ostream &os) {
  size_t nbCase = WMVS.ListNbCase.size();
  os << "WMS:   nbRow=" << WMVS.nbRow << " nbWeight=" << WMVS.nbWeight
     << " nbCase=" << nbCase << "\n";
  os << "WMS:   ListWeight =";
  for (size_t iWei = 0; iWei < WMVS.nbWeight; iWei++) {
    os << " (" << iWei << "," << WMVS.ListWeight[iWei] << ")";
  }
  os << "\n";
  //
  for (size_t iCase = 0; iCase < nbCase; iCase++) {
    std::vector<int> eCase = WMVS.ListPossibleSignatures[iCase];
    os << "WMS:   iCase=" << iCase << "/" << nbCase
       << " nb=" << WMVS.ListNbCase[iCase] << " eCase=" << eCase[0];
    size_t len = eCase.size() / 2;
    os << " LV=";
    for (size_t i = 0; i < len; i++) {
      int first = eCase[1 + 2 * i];
      int second = eCase[2 + 2 * i];
      os << " [" << first << "," << second << "]";
    }
    os << "\n";
  }
}

template <typename T, typename F1, typename F2>
WeightMatrixVertexSignatures<T>
ComputeVertexSignatures(size_t nbRow, F1 f1, F2 f2,
                        [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  MicrosecondTime time;
#endif
  UNORD_MAP_SPECIFIC<T, int> ValueMap_T;
  std::vector<T> ListWeight;
  int idxWeight = 0;
  auto get_T_idx = [&](const T &eval) -> int {
    int &idx = ValueMap_T[eval];
    if (idx == 0) {
      // value is missing
      idxWeight++;
      idx = idxWeight;
      ListWeight.push_back(eval);
    }
    return idx - 1;
  };
  UNORD_MAP_SPECIFIC<std::vector<int>, int> ValueMap_Tvs;
  std::vector<std::vector<int>> ListPossibleSignatures;
  int idxSign = 0;
  auto get_Tvs_idx = [&](std::vector<int> const &esign) -> int {
    int &idx = ValueMap_Tvs[esign];
    if (idx == 0) {
      // value is missing
      idxSign++;
      idx = idxSign;
      ListPossibleSignatures.push_back(esign);
    }
    return idx - 1;
  };
  std::vector<int> ListSignatureByVertex(nbRow);
  int len = 0;
  std::vector<int> list_mult;
  std::vector<int> esign;
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    f1(iRow);
    int idx_diagonal;
    for (size_t jRow = 0; jRow < nbRow; jRow++) {
      T val = f2(jRow);
      int idx = get_T_idx(val);
      if (idx == len) {
        list_mult.push_back(0);
        len++;
      }
      if (iRow != jRow) {
        list_mult[idx]++;
      } else {
        idx_diagonal = idx;
      }
    }
    esign.push_back(idx_diagonal);
    for (int u = 0; u < len; u++) {
      if (list_mult[u] > 0) {
        esign.push_back(u);
        esign.push_back(list_mult[u]);
        list_mult[u] = 0;
      }
    }
    int idx_sign = get_Tvs_idx(esign);
    esign.clear();
    ListSignatureByVertex[iRow] = idx_sign;
  }
  size_t nbWeight = ListWeight.size();
  //
  size_t nbCase = ListPossibleSignatures.size();
  std::vector<int> ListNbCase(nbCase, 0);
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    int iCase = ListSignatureByVertex[iRow];
    ListNbCase[iCase]++;
  }
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  os << "|WMS: ComputeVertexSignature|=" << time << "\n";
#endif
  return {nbRow,
          nbWeight,
          std::move(ListWeight),
          std::move(ListPossibleSignatures),
          std::move(ListSignatureByVertex),
          std::move(ListNbCase)};
}

template <typename T, bool is_symm, typename F1, typename F2, typename F1tr, typename F2tr>
inline typename std::enable_if<is_symm, PairWeightMatrixVertexSignatures<T>>::type
ComputePairVertexSignatures(size_t nbRow, bool canonically, F1 f1, F2 f2,
                            [[maybe_unused]] F1tr f1tr, [[maybe_unused]] F2tr f2tr,
                            std::ostream &os) {
  WeightMatrixVertexSignatures<T> WMVS_direct = ComputeVertexSignatures<T>(nbRow, f1, f2, os);
  if (canonically) {
    RenormalizeWMVS(WMVS_direct, os);
  }
  WeightMatrixVertexSignatures<T> WMVS_dual = EmptyWeightMatrixVertexSignatures<T>();
  return {std::move(WMVS_direct), std::move(WMVS_dual)};
}

template <typename T>
void RenormalizeWMVS(WeightMatrixVertexSignatures<T> &WMVS,
                     [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  MicrosecondTime time;
#endif
  // Building the permutation on the weights
  std::pair<std::vector<T>, std::vector<int>> rec_pair =
      GetReorderingInfoWeight(WMVS.ListWeight);
  size_t n_Wei = WMVS.ListWeight.size();
  const std::vector<int> &g = rec_pair.second;
  WMVS.ListWeight = rec_pair.first;
  // Changing the list of signatures
  std::vector<std::vector<int>> NewListPossibleSignatures;
  for (auto &ePossSignature : WMVS.ListPossibleSignatures) {
    int NewDiag = g[ePossSignature[0]];
    size_t len = ePossSignature.size() / 2;
    std::vector<int> list_mult(n_Wei, 0);
    for (size_t i = 0; i < len; i++) {
      int NewVal = g[ePossSignature[1 + 2 * i]];
      int eMult = ePossSignature[2 + 2 * i];
      list_mult[NewVal] = eMult;
    }
    //
    std::vector<int> newsign;
    newsign.push_back(NewDiag);
    for (size_t i = 0; i < n_Wei; i++) {
      if (list_mult[i] > 0) {
        newsign.push_back(i);
        newsign.push_back(list_mult[i]);
      }
    }
    NewListPossibleSignatures.push_back(newsign);
  }
  WMVS.ListPossibleSignatures = NewListPossibleSignatures;
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  os << "|WMS: RenormalizeWMVS|=" << time << "\n";
#endif
}

template <typename T, bool is_symm, typename F1, typename F2, typename F1tr, typename F2tr>
inline typename std::enable_if<!is_symm, PairWeightMatrixVertexSignatures<T>>::type
ComputePairVertexSignatures(size_t nbRow, [[maybe_unused]] bool canonically, F1 f1, F2 f2,
                            F1tr f1tr, F2tr f2tr,
                            std::ostream &os) {
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
  std::vector<std::vector<T>> Mat_direct;
  for (size_t i=0; i<nbRow; i++) {
    std::vector<T> eV;
    f1(i);
    for (size_t j=0; j<nbRow; j++) {
      T val = f2(j);
      eV.push_back(val);
    }
    Mat_direct.push_back(eV);
  }
  std::vector<std::vector<T>> Mat_dual;
  for (size_t i=0; i<nbRow; i++) {
    std::vector<T> eV;
    f1tr(i);
    for (size_t j=0; j<nbRow; j++) {
      T val = f2tr(j);
      eV.push_back(val);
    }
    Mat_dual.push_back(eV);
  }
  size_t n_error_direct_dual = 0;
  size_t n_nonsymm_direct = 0;
  size_t n_nonsymm_dual = 0;
  for (size_t i=0; i<nbRow; i++) {
    for (size_t j=0; j<nbRow; j++) {
      T val1 = Mat_direct[i][j];
      T val2 = Mat_dual[j][i];
      if (val1 != val2) {
        n_error_direct_dual += 1;
      }
      if (Mat_direct[i][j] != Mat_direct[j][i]) {
        n_nonsymm_direct += 1;
      }
      if (Mat_dual[i][j] != Mat_dual[j][i]) {
        n_nonsymm_dual += 1;
      }
    }
  }
  os << "WMS: n_error_direct_dual=" << n_error_direct_dual << "\n";
  os << "WMS: n_nonsymm_direct=" << n_nonsymm_direct << "\n";
  os << "WMS: n_nonsymm_dual=" << n_nonsymm_dual << "\n";
  if (n_error_direct_dual > 0) {
    std::cerr << "We have n_error_direct_dual=" << n_error_direct_dual << "\n";
    throw TerminalException{1};
  }
#endif
  WeightMatrixVertexSignatures<T> WMVS_direct = ComputeVertexSignatures<T>(nbRow, f1, f2, os);
  RenormalizeWMVS(WMVS_direct, os);
  WeightMatrixVertexSignatures<T> WMVS_dual = ComputeVertexSignatures<T>(nbRow, f1tr, f2tr, os);
  RenormalizeWMVS(WMVS_dual, os);
#ifdef SANITY_CHECK_WEIGHT_MATRIX_SPECIFIED
  if (WMVS_direct.ListWeight != WMVS_dual.ListWeight) {
    std::cerr << "WMS: We have WMVS_direct.ListWeight != WMVS_dual.ListWeight which is not what we want\n";
    throw TerminalException{1};
  }
#endif
  return {std::move(WMVS_direct), std::move(WMVS_dual)};
}

/*
  See WeightMatrix.h : get_effective_weight_index for underlying logic
 */
template <typename T, bool is_symm, typename F2>
inline typename std::enable_if<is_symm, size_t>::type
evaluate_f2(size_t nbRow, size_t nWei, std::unordered_map<T, size_t> const &map,
            size_t iVert, size_t jVert, F2 f2) {
  if (jVert == nbRow + 1) {
    if (iVert == nbRow) {
      return nWei;
    } else {
      return nWei + 1;
    }
  } else {
    if (jVert == nbRow) {
      return map.at(f2(iVert));
    } else {
      return map.at(f2(jVert));
    }
  }
}

template <typename T, bool is_symm, typename F2>
inline typename std::enable_if<!is_symm, size_t>::type
evaluate_f2(size_t nbRow, size_t nWei, std::unordered_map<T, size_t> const &map,
            size_t iVert, size_t jVert, F2 f2) {
#ifdef SANITY_CHECK_WEIGHT_MATRIX_SPECIFIED
  if (iVert >= jVert) {
    std::cerr << "WMS: iVert=" << iVert << " jVert=" << jVert << " but we should have iVert < jVert\n";
    throw TerminalException{1};
  }
#endif
  size_t iVertRed = iVert % nbRow;
  size_t jVertRed = jVert % nbRow;
  if (jVert == 2 * nbRow) {
    // This is the diagonal case of the symmetries
    if (iVert < nbRow) {
      // Case C
      return map.at(f2(iVert));
    } else {
      // Case D
      return nWei + 1;
    }
  }
#ifdef SANITY_CHECK_WEIGHT_MATRIX_SPECIFIED
  if (iVert >= 2 * nbRow) {
    std::cerr << "WMS: We should have iVert < 2 * nbRow\n";
    throw TerminalException{1};
  }
  if (jVert >= 2 * nbRow) {
    std::cerr << "WMS: We should have jVert < 2 * nbRow\n";
    throw TerminalException{1};
  }
#endif
  if (iVertRed == jVertRed) {
#ifdef SANITY_CHECK_WEIGHT_MATRIX_SPECIFIED
    if (iVert + nbRow != jVert) {
      std::cerr << "WMS: iVert=" << iVert << " jVert=" << jVert << " but we should have iVert + nbRow = jVert\n";
      throw TerminalException{1};
    }
#endif
    // Because iVert < jVert, that case can only occurs if iVert = v1 and jVert
    // = v2 for some v. Case B
    return nWei;
  }
  if (jVert < nbRow) {
    // Because iVert < jVert, that case only occurs if iVert = v1, jVert = w1
    // with v != w Case A1
    return nWei + 1;
  }
  if (nbRow <= iVert) {
    // Because iVert < jVert, that case only occurs if iVert = v2, jVert = w2
    // with v != w Case A2
    return nWei + 2;
  }
#ifdef SANITY_CHECK_WEIGHT_MATRIX_SPECIFIED
  if (iVert >= nbRow || nbRow > jVert) {
    std::cerr << "WMS: iVert=" << iVert << " jVert=" << jVert << " for the Case E entry\n";
    throw TerminalException{1};
  }
#endif
  // From previous check, we are now in the situation where iVert = v1 and jVert
  // = w2 with v != w. Case E
  return map.at(f2(jVertRed));
}

/*
  When computing the big graph, ze need to have the symbolic information.
  Those information allow to determine the number of adjacency for each vertex.
  -- For each vertex, add the possible contribution of other vertices.
  -- Expand the list of vertex_to_signature to cover everything.
 */
struct ExpandedSymbolic {
  std::vector<std::vector<std::pair<int, int>>> list_signature;
  std::vector<int> vertex_to_signature;
};

template <typename T, bool is_symm>
inline typename std::enable_if<is_symm, ExpandedSymbolic>::type
get_expanded_symbolic(size_t nbWeight,
                      PairWeightMatrixVertexSignatures<T> const &PairWMVS,
                      [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  MicrosecondTime time;
#endif
  WeightMatrixVertexSignatures<T> const &WMVS = PairWMVS.WMVS_direct;
  size_t nbRow = WMVS.nbRow;
  // The old vertices, but there ar now two more vertices being added
  std::vector<std::vector<std::pair<int, int>>> list_signature;
  for (auto &esign : WMVS.ListPossibleSignatures) {
    std::vector<std::pair<int, int>> e_vect;
    size_t len = esign.size() / 2;
    for (size_t u = 0; u < len; u++)
      e_vect.push_back({esign[1 + 2 * u], esign[2 + 2 * u]});
    // We use a O(n) algorithm since we insert just one entry.
    auto fInsert = [&](int e_val) -> void {
      for (auto &e_pair : e_vect) {
        if (e_pair.first == e_val) {
          e_pair.second++;
          return;
        }
      }
      e_vect.push_back({e_val, 1});
    };
    int diag_val = esign[0]; // The diagonal value
    fInsert(diag_val);
    fInsert(nbWeight + 1);
    list_signature.push_back(e_vect);
  }
  // nbRow : Adding the column corresponding to the diagonal values
  std::unordered_map<int, int> map_vert_nRow;
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    int isign = WMVS.ListSignatureByVertex[iRow];
    int iweight_diag = WMVS.ListPossibleSignatures[isign][0];
    map_vert_nRow[iweight_diag]++;
  }
  map_vert_nRow[nbWeight]++;
  std::vector<std::pair<int, int>> f_vect;
  for (auto &kv : map_vert_nRow)
    f_vect.push_back({kv.first, kv.second});
  list_signature.push_back(f_vect);
  // nbRow+1 : Adding the column corresponding to the external vertex
  std::vector<std::pair<int, int>> g_vect{{nbWeight, 1}, {nbWeight + 1, nbRow}};
  list_signature.push_back(g_vect);
  size_t nbCase = list_signature.size();
  // Now adding the list_of index
  std::vector<int> vertex_to_signature = WMVS.ListSignatureByVertex;
  // nbCase - 2 corresponds to for nbRow
  // nbCase - 1 corresponds to for nbRow + 1
  vertex_to_signature.push_back(nbCase - 2);
  vertex_to_signature.push_back(nbCase - 1);
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  os << "|WMS: get_expanded_symbolic(is_symm=true)|=" << time << "\n";
#endif
  return {std::move(list_signature), std::move(vertex_to_signature)};
}

// The computation for the non-symmetric;
// It corresponds to the code in get_effective_weight_index in WeightMatrix.h
// The matrix is built according to following ideas:
// | A B C |
// |   E F |
// |     I |
// * C is a column matrix formed by WMat.GetValue(iVert, iVert) (Case C)
// * F is nWei + 1 (Case D)
// * A is nWei + 1 all around (Case A1)
// * E is nWei + 2 all around (Case A2)
// * B is nWei on the diagonal
// *      and WMat.GetValue(iVert, jVert) off the diagonal.
//
// The diagonal value is available via WMVS in the first value.
template <typename T, bool is_symm>
inline typename std::enable_if<!is_symm, ExpandedSymbolic>::type
get_expanded_symbolic(size_t nWei, PairWeightMatrixVertexSignatures<T> const &PairWMVS,
                      [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  MicrosecondTime time;
#endif
  WeightMatrixVertexSignatures<T> const &WMVS_direct = PairWMVS.WMVS_direct;
  WeightMatrixVertexSignatures<T> const &WMVS_dual = PairWMVS.WMVS_dual;
  size_t nbRow = WMVS_direct.nbRow;
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
  os << "WMS: nWei=" << nWei << " nbRow=" << nbRow << "\n";
  os << "WMS: WMVS_direct=\n";
  PrintWMVS(WMVS_direct, os);
  os << "WMS: WMVS_dual=\n";
  PrintWMVS(WMVS_dual, os);
#endif
  std::vector<std::vector<std::pair<int, int>>> list_signature;
  std::vector<int> vertex_to_signature;
  auto fUpdate = [](std::vector<std::pair<int, int>> &e_vect, int i_weight,
                    int shift) -> void {
    // We use a O(n) algorithm since we update only few entries. shift can be
    // negative.
    for (auto &e_pair : e_vect) {
      if (e_pair.first == i_weight) {
        e_pair.second += shift;
        return;
      }
    }
    e_vect.push_back({i_weight, shift}); // Thus shift shouldn be positive here
  };
  size_t n_sign_direct = WMVS_direct.ListPossibleSignatures.size();
  size_t n_sign_dual = WMVS_dual.ListPossibleSignatures.size();
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
  os << "WMS: n_sign_direct=" << n_sign_direct << " n_sign_dual=" << n_sign_dual << "\n";
  auto f_check = [&](std::vector<std::pair<int, int>> const &e_vect,
                     int the_case) -> void {
    size_t sum_vert = 0;
    for (auto &ent : e_vect) {
      sum_vert += ent.second;
    }
    if (sum_vert != 2 * nbRow) {
      std::cerr << "case=" << the_case << "\n";
      std::cerr << "sum_vert=" << sum_vert
                << " but should be 2*nbRow=" << (2 * nbRow) << "\n";
      throw TerminalException{1};
    }
  };
#endif
  // The vertices are written first with a block of v1 then a block of v2 and
  // then spec = 2 * nbRow First block of vertices; updates:
  // -- Add the weight nWei+1 for nbRow-1 times (the (v1, w1) adjacencies for v
  // != w)
  // -- Add the weight nWei for just one time (the (v1,v2) adjacency)
  // -- Add the diagonal weight just one time (on
  for (auto &esign : WMVS_direct.ListPossibleSignatures) {
    // Our vertex is named V and is of the form v1
    std::vector<std::pair<int, int>> e_vect;
    size_t len = esign.size() / 2;
    // The off diagonal values of the old matrix. So adjacencies (V = v1, w2)
    // with v <> w.
    for (size_t u = 0; u < len; u++) {
      e_vect.push_back({esign[1 + 2 * u], esign[2 + 2 * u]});
    }
    // Adjacency (V = v1, v2) that we missed from the above
    fUpdate(e_vect, nWei, 1);
    // Adjacencies (V = v1, w1) with v  <> w.
    fUpdate(e_vect, nWei + 1, nbRow - 1);
    // Adjacency (V, spec = 2*nbRow)
    int diag_val = esign[0]; // The diagonal value
    fUpdate(e_vect, diag_val, 1);
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
    f_check(e_vect, 1);
#endif
    list_signature.push_back(e_vect);
  }
  // Second block of vertices; updates:
  // -- Add the value nWei + 2 for nbRow-1 times (the (v2,w2) adjacencies for v
  // != w)
  // -- Subtract the diagonal value by 1 (the diagonal value is replaced by nWei
  // and not moved)
  // -- Add the value mWei for just one time (the (v1,v2) adjacency)
  // -- Add the value nWei + 1 for two times (PURE GUESS at this point)
  for (auto &esign : WMVS_dual.ListPossibleSignatures) {
    // Our vertex is of the form V = v2
    std::vector<std::pair<int, int>> e_vect;
    size_t len = esign.size() / 2;
    // The adjacencies (V = v2, w1) with v <> w
    for (size_t u = 0; u < len; u++) {
      e_vect.push_back({esign[1 + 2 * u], esign[2 + 2 * u]});
    }
    // The adjacency (V = v2, w1) that we missed from above
    fUpdate(e_vect, nWei, 1);
    // The adjacency (V = v2, w2) for v<> w
    fUpdate(e_vect, nWei + 2, nbRow - 1);
    // The adjacency (V = v2, spec)
    fUpdate(e_vect, nWei + 1, 1);
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
    f_check(e_vect, 2);
#endif
    list_signature.push_back(e_vect);
  }
  // Third vertex added:
  // -- We take the diagonal values
  // -- The value nWei + 1 is for nbRow times
  std::unordered_map<int, int> map_vert_nRow;
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    int isign = WMVS_direct.ListSignatureByVertex[iRow];
    int iweight_diag = WMVS_direct.ListPossibleSignatures[isign][0];
    map_vert_nRow[iweight_diag]++;
  }
  map_vert_nRow[nWei + 1] = nbRow;
  std::vector<std::pair<int, int>> f_vect;
  for (auto &kv : map_vert_nRow) {
    f_vect.push_back({kv.first, kv.second});
  }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
  f_check(f_vect, 3);
#endif
  list_signature.push_back(f_vect);
  // Now the signatures.
  for (size_t i_row = 0; i_row < nbRow; i_row++) {
    int iCaseDirect = WMVS_direct.ListSignatureByVertex[i_row];
    vertex_to_signature.push_back(iCaseDirect);
  }
  for (size_t i_row = 0; i_row < nbRow; i_row++) {
    int iCaseDual = WMVS_dual.ListSignatureByVertex[i_row];
    int iCaseNew = iCaseDual + n_sign_direct;
    vertex_to_signature.push_back(iCaseNew);
  }
  vertex_to_signature.push_back(n_sign_direct + n_sign_dual);
  // Moving the results.
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  os << "|WMS: get_expanded_symbolic(is_symm=false)|=" << time << "\n";
#endif
  return {std::move(list_signature), std::move(vertex_to_signature)};
}

template <typename T, bool is_symm, typename F1, typename F2>
SimplifiedVertexColoredGraph
GetSimplifiedVCG(F1 f1, F2 f2, PairWeightMatrixVertexSignatures<T> const &PairWMVS,
                 std::ostream &os) {
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  MicrosecondTime time;
#endif
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
  os << "WMS: Beginning of GetSimplifiedVCG\n";
#endif
  WeightMatrixVertexSignatures<T> const &WMVS = PairWMVS.WMVS_direct;
  size_t nbRow = WMVS.nbRow;
  size_t nbWeight = WMVS.nbWeight;
  std::vector<T> const &ListWeight = WMVS.ListWeight;
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
  os << "WMS: GetSimplifiedVCG, ListWeight=" << ListWeight << "\n";
#endif
  std::unordered_map<T, size_t> map;
  for (size_t iWei = 0; iWei < nbWeight; iWei++) {
    map[ListWeight[iWei]] = iWei;
  }
  //
  size_t nbMult = get_effective_nb_weight<is_symm>(nbWeight);
  size_t hS = Pairs_GetNeededN(nbMult);
  size_t nbVert = get_effective_nb_vert<is_symm>(nbRow);
  size_t nbVertTot = nbVert * hS;
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
  os << "WMS: nbMult=" << nbMult << " hS=" << hS << " nbRow=" << nbRow
     << " nbVert=" << nbVert << " nbVertTot=" << nbVertTot << "\n";
#endif
  //
  std::vector<int> V = Pairs_GetListPair(hS, nbMult);
  size_t e_pow = V.size() / 2;
  ExpandedSymbolic expand =
      get_expanded_symbolic<T, is_symm>(nbWeight, PairWMVS, os);
  size_t nbCase = expand.list_signature.size();
  //
  // Determining the number of cases
  //
  std::vector<int> ListNbCase(nbCase, 0);
  for (size_t i = 0; i < nbVert; i++) {
    int iCase = expand.vertex_to_signature[i];
    ListNbCase[iCase]++;
  }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
  os << "WMS: |vertex_to_signature|=" << expand.vertex_to_signature.size()
     << "\n";
  os << "WMS: |List_signature|=" << expand.list_signature.size() << "\n";
  for (size_t iCase = 0; iCase < nbCase; iCase++) {
    os << "WMS:   iCase=" << iCase << "/" << nbCase << " ListNbCase=" << ListNbCase[iCase] << " V =";
    for (auto &ent : expand.list_signature[iCase]) {
      os << " (" << ent.first << "," << ent.second << ")";
    }
    os << "\n";
  }
#endif
  //
  // Now computing the weight by multiplier
  //
  std::vector<size_t> WeightByMult(nbMult, 0);
  for (size_t iCase = 0; iCase < nbCase; iCase++) {
    int sizCase = ListNbCase[iCase];
    std::vector<std::pair<int, int>> const &e_vect =
        expand.list_signature[iCase];
    for (auto &e_pair : e_vect) {
      int iWeight = e_pair.first;
      int eMult = e_pair.second;
      WeightByMult[iWeight] += sizCase * eMult;
    }
  }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
  os << "WMS: nbMult=" << nbMult << " WeightByMult(A)=";
  for (size_t iMult = 0; iMult < nbMult; iMult++) {
    os << WeightByMult[iMult] << " ";
  }
  os << "\n";
#endif
  //
  // Reordering the multiplier to maximize the sparsity
  // (This is not perfect since if iH1 != iH2 we need to add more edges. See
  // below) (There are more strange things happening. We should have gotten some
  // improvement accross the board in memory usage, but we do not)
  //
  std::vector<size_t> ListIdx(nbMult);
  for (size_t iMult = 0; iMult < nbMult; iMult++) {
    ListIdx[iMult] = iMult;
  }
  std::stable_sort(ListIdx.begin(), ListIdx.end(),
                   [&](const size_t &idx1, const size_t &idx2) -> bool {
                     return WeightByMult[idx1] > WeightByMult[idx2];
                   });
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
  os << "WMS: After sort\n";
  os << "WMS: nbMult=" << nbMult << " WeightByMult(B)=";
  for (size_t iMult = 0; iMult < nbMult; iMult++) {
    os << WeightByMult[ListIdx[iMult]] << " ";
  }
  os << "\n";
#endif
  //
  // Now computing the list of signature
  //
  MyMatrix<size_t> MatrixAdj = ZeroMatrix<size_t>(hS, nbCase);
  // Adjacencies (iVert + nbVert*iH , iVert + nbVert*jH)
  size_t nb_adj1 = hS - 1;
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
  os << "WMS: hS=" << hS << " nb_adj1=" << nb_adj1 << "\n";
#endif
  for (size_t iH = 0; iH < hS; iH++) {
    for (size_t iCase = 0; iCase < nbCase; iCase++) {
      MatrixAdj(iH, iCase) += nb_adj1;
    }
  }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
  auto f_print = [&](std::ostream& os_o) -> void {
    os_o << "WMS: hS=" << hS << " nbCase=" << nbCase << " MatrixAdj=\n";
    for (size_t iH = 0; iH < hS; iH++) {
      os_o << "WMS:   iH=" << iH << " MatrixAdj =";
      for (size_t iCase = 0; iCase < nbCase; iCase++) {
        os_o << " " << MatrixAdj(iH, iCase);
      }
      os_o << "\n";
    }
    os_o << "WMS: nbMult=" << nbMult << "\n";
    os_o << "WMS: ListIdx =";
    for (size_t iMult = 0; iMult < nbMult; iMult++) {
      os_o << " " << ListIdx[iMult];
    }
    os_o << "\n";
    os_o << "WMS: WeightByMult =";
    for (size_t iMult = 0; iMult < nbMult; iMult++) {
      os_o << " " << WeightByMult[iMult];
    }
    os_o << "\n";
  };
#endif
  //
  // Now computing the contribution of the adjacencies
  //
  Face f_total = GetAllBinaryExpressionsByWeight(e_pow, nbMult);
  for (size_t iCase = 0; iCase < nbCase; iCase++) {
    std::vector<std::pair<int, int>> const &e_vect =
        expand.list_signature[iCase];
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
    for (auto &epair : e_vect) {
      if (epair.second < 0) {
        std::cerr << "WMS: iCase=" << iCase << "\n";
        f_print(std::cerr);
        std::cerr
            << "A negative multiplicity is inserted which is not allowed\n";
        throw TerminalException{1};
      }
    }
#endif
    for (auto &e_pair : e_vect) {
      int iWeight = ListIdx[e_pair.first];
      size_t eMult = e_pair.second;
      size_t shift = iWeight * e_pow;
      for (size_t i_pow = 0; i_pow < e_pow; i_pow++) {
        if (f_total[shift + i_pow] == 1) {
          int iH1 = V[2 * i_pow];
          int iH2 = V[2 * i_pow + 1];
          // Adjacency (iVert + nbVert*iH1 , jVert + nbVert*iH2)
          MatrixAdj(iH1, iCase) += eMult;
          if (iH1 != iH2) {
            // Adjacency (iVert + nbVert*iH2 , jVert + nbVert*iH1)
            MatrixAdj(iH2, iCase) += eMult;
          }
        }
      }
    }
  }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
  os << "WMS: MatrixAdj built\n";
#endif
  //
  // Now building the adjacencies for traces
  //
  size_t nbAdjacent = 0;
  for (size_t iCase = 0; iCase < nbCase; iCase++) {
    for (size_t iH = 0; iH < hS; iH++) {
      nbAdjacent += ListNbCase[iCase] * MatrixAdj(iH, iCase);
    }
  }

  SimplifiedVertexColoredGraph s =
      GetSimplifiedVertexColoredGraph(nbVertTot, nbAdjacent, hS);
  // Now setting up the adjacencies
  size_t pos = 0;
  std::vector<size_t> ListShift(nbVertTot);
  for (size_t i = 0; i < nbVertTot; i++) {
    size_t iVert = i % nbVert;
    size_t iH = i / nbVert;
    int iCase = expand.vertex_to_signature[iVert];
    int nbAdj = MatrixAdj(iH, iCase);
    s.d[i] = nbAdj;
    ListShift[i] = pos;
    pos += nbAdj;
  }
  // Assigning the color stuff
  for (size_t iH = 0; iH < hS; iH++) {
    s.ListBlockSize[iH] = nbVert;
  }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
  os << "WMS: SimplifiedVCG stuff set up\n";
#endif
  //
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
  std::vector<int> ListDegExpe(nbVertTot, 0);
#endif
  auto f_adj = [&](size_t iVert, size_t jVert) -> void {
    s.e[ListShift[iVert]] = jVert;
    ListShift[iVert]++;
    s.e[ListShift[jVert]] = iVert;
    ListShift[jVert]++;
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
    ListDegExpe[iVert]++;
    ListDegExpe[jVert]++;
#endif
  };
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    for (size_t iH = 0; iH < hS - 1; iH++) {
      for (size_t jH = iH + 1; jH < hS; jH++) {
        size_t aVert = iVert + nbVert * iH;
        size_t bVert = iVert + nbVert * jH;
        f_adj(aVert, bVert);
      }
    }
  }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
  MyMatrix<size_t> Mval(nbVert, nbVert);
  for (size_t iVert = 0; iVert < nbVert - 1; iVert++) {
    if (iVert < nbRow) {
      // Both for is_symm = true/false, only the iVert < nbRow needs to be
      // computed.
      f1(iVert);
    }
    for (size_t jVert = iVert + 1; jVert < nbVert; jVert++) {
      size_t val =
          evaluate_f2<T, is_symm>(nbRow, nbWeight, map, iVert, jVert, f2);
      Mval(iVert, jVert) = val;
      Mval(jVert, iVert) = val;
    }
  }
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    std::set<size_t> indices;
    std::map<size_t, size_t> map_eval;
    for (size_t jVert = 0; jVert < nbVert; jVert++) {
      if (iVert != jVert) {
        size_t val = Mval(iVert, jVert);
        map_eval[val] += 1;
        indices.insert(val);
      }
    }
    std::map<size_t, size_t> map_matrix;
    int iCase = expand.vertex_to_signature[iVert];
    std::vector<std::pair<int, int>> const &e_vect =
        expand.list_signature[iCase];
    for (auto & epair : e_vect) {
      size_t val = epair.first;
      map_matrix[val] = epair.second;
      indices.insert(val);
    }
    size_t n_error = 0;
    for (auto & val : indices) {
      size_t mult_eval = map_eval[val];
      size_t mult_matrix = map_matrix[val];
      if (mult_eval != mult_matrix) {
        n_error += 1;
      }
    }
    if (n_error > 0) {
      os << "WMS: iVert=" << iVert << " mult_disc =";
      for (auto & val : indices) {
        size_t mult_eval = map_eval[val];
        size_t mult_matrix = map_matrix[val];
        if (mult_eval != mult_matrix) {
          os << " (" << val << "," << mult_eval << "/" << mult_matrix << ")";
        }
      }
      os << "\n";
      os << "WMS:     iCase=" << iCase;
      for (auto & epair : e_vect) {
        os << " (" << epair.first << "," << epair.second << ")";
      }
      os << "\n";
      os << "WMS:     eff_mult=";
      for (auto & kv : map_eval) {
        if (kv.second > 0) {
          os << " (" << kv.first << "," << kv.second << ")";
        }
      }
      os << "\n";
      std::cerr << "WMS: That is not matching\n";
      throw TerminalException{1};
    }
  }
#endif
  for (size_t iVert = 0; iVert < nbVert - 1; iVert++) {
    if (iVert < nbRow) {
      // Both for is_symm = true/false, only the iVert < nbRow needs to be
      // computed.
      f1(iVert);
    }
    for (size_t jVert = iVert + 1; jVert < nbVert; jVert++) {
      size_t val =
          evaluate_f2<T, is_symm>(nbRow, nbWeight, map, iVert, jVert, f2);
      size_t shift = e_pow * ListIdx[val];
      for (size_t i_pow = 0; i_pow < e_pow; i_pow++) {
        if (f_total[shift + i_pow] == 1) {
          int iH1 = V[2 * i_pow];
          int iH2 = V[2 * i_pow + 1];
          size_t aVert = iVert + nbVert * iH1;
          size_t bVert = jVert + nbVert * iH2;
          f_adj(aVert, bVert);
          if (iH1 != iH2) {
            aVert = iVert + nbVert * iH2;
            bVert = jVert + nbVert * iH1;
            f_adj(aVert, bVert);
          }
        }
      }
    }
  }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
  int sum_adj = 0;
  int nb_error01 = 0;
  size_t sum_deg0 = 0;
  size_t sum_deg1 = 0;
  std::unordered_set<std::pair<size_t, size_t>> set_case_error;
  for (size_t i = 0; i < nbVertTot; i++) {
    size_t iVert = i % nbVert;
    size_t iH = i / nbVert;
    int deg0 = s.d[i];
    int deg1 = ListDegExpe[i];
    sum_deg0 += static_cast<size_t>(deg0);
    sum_deg1 += static_cast<size_t>(deg1);
    if (deg0 != deg1) {
      int delta = deg0 - deg1;
      std::cerr << "i=" << i << " iVert=" << iVert << " iH=" << iH
                << " deg0=" << deg0
                << " deg1=" << deg1;
      if (deg0 > deg1) {
        std::cerr << " A";
      } else {
        std::cerr << " B";
      }
      std::cerr << " delta=" << delta << "\n";
      std::pair<size_t, size_t> pair{iVert, iH};
      set_case_error.insert(pair);
      nb_error01++;
    }
    sum_adj += deg1;
  }
  double frac_adj =
      static_cast<double>(sum_adj) /
      (static_cast<double>(nbVertTot) * static_cast<double>(nbVertTot - 1));
  os << "WMS: sum_adj=" << sum_adj << " nbAdjacent=" << nbAdjacent
     << " frac_adj=" << frac_adj << "\n";
  std::cerr << "WMS: nb_error01=" << nb_error01 << "\n";
  std::cerr << "WMS: sum_deg0=" << sum_deg0 << "  sum_deg1=" << sum_deg1 << " nbAdjacent=" << nbAdjacent << "\n";
  /*
  for (auto & case_error : set_case_error) {
    std::cerr << "case_error: iVert=" << case_error.first << " iH=" << case_error.second << "\n";
  }
  */
  if (nb_error01 > 0) {
    std::cerr << "some error occurred\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  os << "|WMS: GetSimplifiedVCG|=" << time << "\n";
#endif
  return s;
}

template <typename TidxC, typename Tidx, bool is_symm>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>>
GetGroupCanonicalization_KnownSignature_TidxC(SimplifiedVertexColoredGraph const& s, size_t const& nbRow, std::ostream& os) {
  std::pair<std::vector<TidxC>, std::vector<std::vector<Tidx>>> ePair =
    GRAPH_GetCanonicalOrdering_ListGenerators_Simp<TidxC, Tidx>(s, nbRow, os);
  std::vector<Tidx> MapVectRev2 =
    GetCanonicalizationVector_KernelBis<Tidx, TidxC, is_symm>(nbRow, ePair.first, os);
  return {std::move(MapVectRev2), std::move(ePair.second)};
}

template <typename T, typename Tidx, bool is_symm, typename F1, typename F2>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>>
GetGroupCanonicalization_KnownSignature(
    PairWeightMatrixVertexSignatures<T> const &PairWMVS, F1 f1, F2 f2,
    std::ostream &os) {
  size_t nbRow = PairWMVS.WMVS_direct.nbRow;
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "GetGroupCanonicalization_KnownSignature : We have nbRow="
              << nbRow << " which is larger than the possible values of Tidx : "
              << max_poss_rows << "\n";
    throw TerminalException{1};
  }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
  os << "WMS: Before calling GetSimplifiedVCG from "
        "GetGroupCanonicalization_KnownSignature\n";
#endif
  SimplifiedVertexColoredGraph s =
      GetSimplifiedVCG<T, is_symm, F1, F2>(f1, f2, PairWMVS, os);
  if (s.nbVert < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return GetGroupCanonicalization_KnownSignature_TidxC<uint8_t,Tidx,is_symm>(s, nbRow, os);
  }
  if (s.nbVert < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return GetGroupCanonicalization_KnownSignature_TidxC<uint16_t,Tidx,is_symm>(s, nbRow, os);
  }
  if (s.nbVert < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return GetGroupCanonicalization_KnownSignature_TidxC<uint32_t,Tidx,is_symm>(s, nbRow, os);
  }
#if !defined __APPLE__
  if (s.nbVert < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return GetGroupCanonicalization_KnownSignature_TidxC<uint64_t,Tidx,is_symm>(s, nbRow, os);
  }
#endif
  std::cerr << "Failed to find matching numeric in "
               "GetGroupCanonicalization_KnownSignature\n";
  throw TerminalException{1};
}

template <typename T, typename Tidx, bool is_symm, typename F1, typename F2>
std::vector<std::vector<Tidx>> GetStabilizerWeightMatrix_KnownSignature(
    PairWeightMatrixVertexSignatures<T> const &PairWMVS, F1 f1, F2 f2,
    std::ostream &os) {
  size_t nbRow = PairWMVS.WMVS_direct.nbRow;
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "GetStabilizerWeightMatrix_KnownSignature : We have nbRow="
              << nbRow << " which is larger than the possible values of Tidx : "
              << max_poss_rows << "\n";
    throw TerminalException{1};
  }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
  os << "WMS: Before calling GetSimplifiedVCG from "
        "GetStabilizerWeightMatrix_KnownSignature\n";
#endif
  SimplifiedVertexColoredGraph s =
      GetSimplifiedVCG<T, is_symm, F1, F2>(f1, f2, PairWMVS, os);
  return GRAPH_GetListGenerators_Simp<Tidx>(s, nbRow, os);
}

/*
  ---F1/F2 : The first/second template function for creating the Weight matrix
  ---F3    : The function for testing acceptability of sets for consideration
  (e.g. rank function)
  ---F4    : The function for testing acceptability of small generators and
  returning the big generators. It also takes blocks in input and determine
  their nature, that is preserved or not
  ---Fproc1 : function processing the graph (and using traces for this)
  ---Fproc2 : function taking the output and returning the list of generators
  ---Fproc3 : function taking all of it and returning the output.
*/
template <bool canonically, bool is_symm, typename T, typename Tidx, typename Tret1,
          typename Tret2, typename Tret3, typename F1, typename F2, typename F1tr, typename F2tr,
          typename F3, typename F4, typename Fproc1, typename Fproc2, typename Fproc3>
Tret3 BlockBreakdown_Heuristic(size_t nbRow, F1 f1, F2 f2, F1tr f1tr, F2tr f2tr,
                               F3 f3, F4 f4,
                               Fproc1 fproc1, Fproc2 fproc2, Fproc3 fproc3,
                               std::ostream &os) {
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "BlockBreakdown_Heuristic : We have nbRow=" << nbRow
              << " which is larger than the possible values of Tidx : "
              << max_poss_rows << "\n";
    throw TerminalException{1};
  }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
  os << "WMS: Beginning of BlockBreakdown_Heuristic\n";
#endif
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  MicrosecondTime time;
#endif
  size_t max_globiter = 1000;
  VertexPartition<Tidx> VP = ComputeVertexPartition<T, Tidx>(
      nbRow, f1, f2, canonically, max_globiter, os);
#ifdef TIMINGS_WEIGHT_MATRIX_SPECIFIED
  os << "|WMS: ComputeVertexPartition|=" << time << "\n";
#endif
  size_t nbCase = VP.ListBlocks.size();
  std::vector<size_t> ListIdx = GetOrdering_ListIdx(VP);
  Face f_covered(nbCase);
  std::vector<int> StatusCase(nbCase);
  size_t idx = 0;
  auto set_status_case = [&]() -> void {
    for (size_t iCase = 0; iCase < nbCase; iCase++) {
      StatusCase[iCase] = 0;
    }
    for (size_t u = 0; u < idx; u++) {
      size_t pos = ListIdx[u];
      if (f_covered[pos] == 0) {
        StatusCase[pos] = 1;
      }
    }
  };
  auto sum_status_case = [&]() -> int {
    int sum = 0;
    for (size_t iCase = 0; iCase < nbCase; iCase++)
      sum += StatusCase[iCase];
    return sum;
  };
  auto increase_idx = [&]() -> bool {
    int sum_prev = sum_status_case();
    while (true) {
      idx++;
      if (idx == nbCase + 1) {
        return false;
      }
      set_status_case();
      int sum_new = sum_status_case();
      if (sum_new > sum_prev) {
        return true;
      }
    }
  };
  while (true) {
    bool test = increase_idx();
    if (!test) {
      break;
    }
    std::vector<Tidx> CurrentListIdx;
    for (size_t iRow = 0; iRow < nbRow; iRow++) {
      int iCase = VP.MapVertexBlock[iRow];
      if (StatusCase[iCase] == 1) {
        CurrentListIdx.push_back(iRow);
      }
    }
    size_t nbRow_res = CurrentListIdx.size();
    //
    bool test_f3 = f3(CurrentListIdx);
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
    os << "WMS: idx=" << idx << " |CurrentListIdx|=" << nbRow_res << " test_f3=" << test_f3 << "\n";
#endif
    if (test_f3) {
      std::vector<std::vector<Tidx>> ListBlocks;
      std::vector<size_t> ListEnt;
      for (size_t u = 0; u < nbCase; u++)
        if (StatusCase[u] == 0) {
          ListBlocks.push_back(VP.ListBlocks[u]);
          ListEnt.push_back(u);
        }
      auto f1_res = [&](size_t iRow) -> void { f1(CurrentListIdx[iRow]); };
      auto f2_res = [&](size_t jRow) -> T { return f2(CurrentListIdx[jRow]); };
      auto f1tr_res = [&](size_t iRow) -> void { f1tr(CurrentListIdx[iRow]); };
      auto f2tr_res = [&](size_t jRow) -> T { return f2tr(CurrentListIdx[jRow]); };
      PairWeightMatrixVertexSignatures<T> PairWMVS_res =
        ComputePairVertexSignatures<T, is_symm>(nbRow_res, canonically, f1_res, f2_res, f1tr_res, f2tr_res, os);
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
      os << "WMS: PairWMVS_res.WMVS_direct=\n";
      PrintWMVS(PairWMVS_res.WMVS_direct, os);
      os << "WMS: PairWMVS_res.WMVS_dual=\n";
      PrintWMVS(PairWMVS_res.WMVS_dual, os);
#endif
      Tret1 ret1 = fproc1(PairWMVS_res, f1_res, f2_res);
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
      os << "WMS: After fproc1\n";
#endif
      bool IsCorrect = true;
      std::vector<std::vector<Tidx>> LGen;
      Face f_incorr(ListEnt.size());
      const Tret2& ret_fproc2 = fproc2(ret1);
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
      os << "WMS: After fproc2 |ret_fproc2|=" << ret_fproc2.size() << "\n";
#endif
      for (auto &eList : ret_fproc2) {
        DataMapping<Tidx> test = f4(CurrentListIdx, eList, ListBlocks);
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
        os << "WMS: After f4 test.correct=" << test.correct << "\n";
#endif
        if (test.correct) {
          LGen.push_back(test.eGen);
        } else {
          IsCorrect = false;
        }
        size_t n_incorr = 0;
        for (size_t u = 0; u < ListEnt.size(); u++) {
          if (!test.block_status[u])
            f_incorr[u] = 1;
          n_incorr += f_incorr[u];
        }
        if (n_incorr == ListEnt.size() && n_incorr > 0) {
          // All entries are incorrect. So, no need to continue
          break;
        }
      }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
      os << "WMS: f_incorr=";
      for (size_t u = 0; u < ListEnt.size(); u++)
        os << int(f_incorr[u]);
      os << "\n";
      os << "WMS: |ListEnt|=" << ListEnt.size() << "\n";
#endif
      for (size_t u = 0; u < ListEnt.size(); u++) {
        size_t pos = ListEnt[u];
        if (f_incorr[u] == 0) {
          f_covered[pos] = 1;
        }
      }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
      os << "WMS: Before IsCorrect test IsCorrect=" << IsCorrect << "\n";
#endif
      if (IsCorrect) {
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
        os << "WMS: Before fproc3\n";
#endif
        return fproc3(CurrentListIdx, ret1, LGen);
      }
    }
  }
  std::cerr << "BlockBreakdown_Heuristic : We should never reach that stage\n";
  throw TerminalException{1};
}

template <typename T, typename Tidx, bool is_symm, typename F1, typename F2,
          typename F1tr, typename F2tr,
          typename F3, typename F4>
std::vector<std::vector<Tidx>>
GetStabilizerWeightMatrix_Heuristic(size_t nbRow, F1 f1, F2 f2, F1tr f1tr, F2tr f2tr,
                                    F3 f3, F4 f4, std::ostream &os) {
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "GetStabilizerWeightMatrix_Heuristic : We have nbRow=" << nbRow
              << " which is larger than the possible values of Tidx : "
              << max_poss_rows << "\n";
    throw TerminalException{1};
  }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
  os << "WMS: Beginning of GetStabilizerWeightMatrix_Heuristic\n";
#endif
  using Tret1 = std::vector<std::vector<Tidx>>;
  using Tret2 = std::vector<std::vector<Tidx>>;
  using Tret3 = std::vector<std::vector<Tidx>>;
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED_EXTENSIVE
  os << "WMS: Thematrix(" << nbRow << "|" << nbRow << ")=\n";
  std::unordered_map<T, size_t> map_T;
  size_t n_val = 0;
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    f1(iRow);
    for (size_t iCol = 0; iCol < nbRow; iCol++) {
      T val = f2(iCol);
      size_t &pos = map_T[val];
      if (pos == 0) {
        n_val++;
        pos = n_val;
      }
    }
  }
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    f1(iRow);
    for (size_t iCol = 0; iCol < nbRow; iCol++) {
      T val = f2(iCol);
      size_t pos = map_T[val] - 1;
      os << " " << pos;
    }
    os << "\n";
  }
  os << "WMS: Beginning of GetStabilizerWeightMatrix_Heuristic\n";
#endif
  auto fproc1 = [&](const PairWeightMatrixVertexSignatures<T> &PairWMVS_res,
                    auto f1_res, auto f2_res) -> Tret1 {
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
    os << "WMS: GetStabilizerWeightMatrix_Heuristic(fproc1) |WMVS_res|="
       << PairWMVS_res.WMVS_direct.nbRow << " nbWeight=" << PairWMVS_res.WMVS_direct.nbWeight << "\n";
#endif
    return GetStabilizerWeightMatrix_KnownSignature<T, Tidx, is_symm>(
        PairWMVS_res, f1_res, f2_res, os);
  };
  auto fproc2 = [&](const Tret1 &ListGen) -> const Tret2 & {
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
    os << "WMS: GetStabilizerWeightMatrix_Heuristic(fproc2) : |ListGen|="
       << ListGen.size() << "\n";
#endif
    return ListGen;
  };
  auto fproc3 = [&]([[maybe_unused]] const std::vector<Tidx> &Vsubset,
                    [[maybe_unused]] const Tret1 &ret1,
                    const std::vector<std::vector<Tidx>> &LGenFinal) -> Tret3 {
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
    os << "WMS: GetStabilizerWeightMatrix_Heuristic(fproc3) : |LGenFinal|="
       << LGenFinal.size() << "\n";
#endif
    return LGenFinal;
  };
  const bool canonically = false;
  return BlockBreakdown_Heuristic<canonically, is_symm, T, Tidx, Tret1, Tret2, Tret3>(nbRow, f1, f2, f1tr, f2tr, f3, f4, fproc1, fproc2, fproc3, os);
}

/*
  ---F1/F2 : The first/second template function for creating the Weight matrix
  ---F1tr/F2tr : The first/second template function for creating the transposed
     Weight matrix
  ---F3    : The function for testing acceptability of sets for consideration
     (e.g. rank function)
  ---F4    : The function for testing acceptability of small generators and
     returning the big generators.
  ---F5    : The lifting of canonicalization from a subset to the full set.
*/
template <typename T, typename Tidx, bool is_symm, typename F1, typename F2,
          typename F1tr, typename F2tr,
          typename F3, typename F4, typename F5>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>>
GetGroupCanonicalizationVector_Heuristic(size_t nbRow, F1 f1, F2 f2, F1tr f1tr, F2tr f2tr,
                                         F3 f3, F4 f4, F5 f5, std::ostream &os) {
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "GetGroupCanonicalizationVector_Heuristic : We have nbRow="
              << nbRow << " which is larger than the possible values of Tidx : "
              << max_poss_rows << "\n";
    throw TerminalException{1};
  }
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
  os << "WMS: Beginning of GetGroupCanonicalizationVector_Heuristic\n";
#endif
  using Tret1 = std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>>;
  using Tret2 = std::vector<std::vector<Tidx>>;
  using Tret3 = std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>>;
  auto fproc1 = [&](const PairWeightMatrixVertexSignatures<T> &PairWMVS_res,
                    auto f1_res, auto f2_res) -> Tret1 {
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
    os << "WMS: GetGroupCanonicalizationVector_Heuristic, before GetGroupCanonicalization_KnownSignature\n";
#endif
    return GetGroupCanonicalization_KnownSignature<T, Tidx, is_symm>(
        PairWMVS_res, f1_res, f2_res, os);
  };
  auto fproc2 = [&](const Tret1 &ePair) -> const Tret2 & {
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
    os << "WMS: GetGroupCanonicalizationVector_Heuristic : |ListGen|="
       << ePair.second.size() << "\n";
#endif
    return ePair.second;
  };
  auto fproc3 = [&](const std::vector<Tidx> &Vsubset, const Tret1 &ePair,
                    const std::vector<std::vector<Tidx>> &LGenFinal) -> Tret3 {
#ifdef DEBUG_WEIGHT_MATRIX_SPECIFIED
    os << "WMS: GetGroupCanonicalizationVector_Heuristic : |LGenFinal|="
       << LGenFinal.size() << "\n";
#endif
    return {f5(Vsubset, ePair.first), LGenFinal};
  };
  const bool canonically = true;
  return BlockBreakdown_Heuristic<canonically, is_symm, T, Tidx, Tret1, Tret2, Tret3>(nbRow, f1, f2, f1tr, f2tr, f3, f4, fproc1, fproc2, fproc3, os);
}

// clang-format off
#endif  // SRC_GROUP_WEIGHTMATRIXSPECIFIED_H_
// clang-format on

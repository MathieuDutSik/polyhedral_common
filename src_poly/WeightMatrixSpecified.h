#ifndef INCLUDE_WEIGHT_MATRIX_SPECIFIED_H
#define INCLUDE_WEIGHT_MATRIX_SPECIFIED_H


// The hash map do not seem to make much difference in the overall
// performance.

//#define UNORDERED_MAP_SPECIFIC
//#define TSL_SPARSE_MAP_SPECIFIC
//#define TSL_ROBIN_MAP_SPECIFIC
#define TSL_HOPSCOTCH_MAP_SPECIFIC

#ifdef UNORDERED_MAP_SPECIFIC
# include <unordered_map>
# include <unordered_set>
# define UNORD_MAP_SPECIFIC std::unordered_map
# define UNORD_SET_SPECIFIC std::unordered_set
#endif

#ifdef TSL_SPARSE_MAP_SPECIFIC
# include "sparse_map.h"
# include "sparse_set.h"
# define UNORD_MAP_SPECIFIC tsl::sparse_map
# define UNORD_SET_SPECIFIC tsl::sparse_set
#endif

#ifdef TSL_ROBIN_MAP_SPECIFIC
# include "robin_map.h"
# include "robin_set.h"
# define UNORD_MAP_SPECIFIC tsl::robin_map
# define UNORD_SET_SPECIFIC tsl::robin_set
#endif

#ifdef TSL_HOPSCOTCH_MAP_SPECIFIC
# include "hopscotch_map.h"
# include "hopscotch_set.h"
# define UNORD_MAP_SPECIFIC tsl::hopscotch_map
# define UNORD_SET_SPECIFIC tsl::hopscotch_set
#endif



using SignVertex = std::vector<int>;


template<typename T>
std::pair<std::vector<T>, std::vector<int>> GetReorderingInfoWeight(const std::vector<T>& ListWeight)
{
  size_t n_Wei = ListWeight.size();
  std::map<T, size_t> map;
  for (size_t i=0; i<n_Wei; i++)
    map[ListWeight[i]] = i;
  std::vector<int> g(n_Wei);
  int idx=0;
  for (auto & kv : map) {
    int pos = kv.second;
    g[pos] = idx;
    idx++;
  }
  std::vector<T> NewListWeight(n_Wei);
  for (size_t iW=0; iW<n_Wei; iW++) {
    size_t jW = g[iW];
    NewListWeight[jW] = ListWeight[iW];
  }
  return {std::move(NewListWeight), std::move(g)};
}



template<typename Tidx>
struct VertexPartition {
  size_t nbRow;
  std::vector<Tidx> MapVertexBlock;
  std::vector<std::vector<Tidx>> ListBlocks;
};


template<typename T, typename Tidx, typename F1, typename F2>
VertexPartition<Tidx> ComputeInitialVertexPartition(size_t nbRow, F1 f1, F2 f2, bool canonically)
{
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "ComputeInitialVertexPartition : We have nbRow=" << nbRow << " which is larger than the possible values of Tidx : " << max_poss_rows << "\n";
    throw TerminalException{1};
  }
  std::vector<T> ListWeight;
  UNORD_MAP_SPECIFIC<T,size_t> ValueMap_T;
  size_t idxWeight = 0;
  auto get_T_idx=[&](T eval) -> Tidx {
    size_t& idx = ValueMap_T[eval];
    if (idx == 0) { // value is missing
      idxWeight++;
      idx = idxWeight;
      ListWeight.push_back(eval);
    }
    return Tidx(idx - 1);
  };
  std::vector<Tidx> MapVertexBlock(nbRow);
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    f1(iRow);
    T val = f2(iRow);
    MapVertexBlock[iRow] = get_T_idx(val);
  }
  if (canonically) {
    std::pair<std::vector<T>, std::vector<int>> rec_pair = GetReorderingInfoWeight(ListWeight);
    const std::vector<int>& g = rec_pair.second;
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      int NewIdx = g[MapVertexBlock[iRow]];
      MapVertexBlock[iRow] = NewIdx;
    }
  }
  std::vector<std::vector<Tidx>> ListBlocks(ListWeight.size());
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    Tidx pos = MapVertexBlock[iRow];
    ListBlocks[pos].push_back(iRow);
  }
  return {nbRow, std::move(MapVertexBlock), std::move(ListBlocks)};
}

// jBlock : is the block for which we look at breaking down
template<typename T, typename Tidx, typename F1, typename F2>
bool RefineSpecificVertexPartition(VertexPartition<Tidx> & VP, const int& jBlock, const int& iBlock, F1 f1, F2 f2, bool canonically)
{
  size_t nbRow = VP.nbRow;
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "RefineSpecificVertexPartition : We have nbRow=" << nbRow << " which is larger than the possible values of Tidx : " << max_poss_rows << "\n";
    throw TerminalException{1};
  }
  // The code for finding the indices
  UNORD_MAP_SPECIFIC<T,int> ValueMap_T;
  std::vector<T> ListWeight;
  int idxWeight = 0;
  auto get_T_idx=[&](const T& eval) -> int {
    int& idx = ValueMap_T[eval];
    if (idx == 0) { // value is missing
      idxWeight++;
      idx = idxWeight;
      ListWeight.push_back(eval);
    }
    return idx - 1;
  };
  UNORD_MAP_SPECIFIC<SignVertex,int> ValueMap_Tvs;
  std::vector<SignVertex> ListPossibleSignatures;
  int idxSign = 0;
  auto get_Tvs_idx=[&](SignVertex const& esign) -> int {
    int& idx = ValueMap_Tvs[esign];
    if (idx == 0) { // value is missing
      idxSign++;
      idx = idxSign;
      ListPossibleSignatures.push_back(esign);
    }
    return idx - 1;
  };
  // Now computing the indices
  std::vector<Tidx> eBlockBreak = VP.ListBlocks[jBlock]; // We do need a copy operation here
  const std::vector<Tidx>& eBlockSpec = VP.ListBlocks[iBlock];
  size_t siz_block_break = eBlockBreak.size();
  size_t siz_block_spec = eBlockSpec.size();
  std::vector<int> V(siz_block_break);
  int len = 0;
  std::vector<int> list_mult;
  std::vector<int> esign;
  for (size_t iBreak=0; iBreak<siz_block_break; iBreak++) {
    size_t iRow = eBlockBreak[iBreak];
    f1(iRow);
    for (size_t iSpec=0; iSpec<siz_block_spec; iSpec++) {
      size_t jRow = eBlockSpec[iSpec];
      T val = f2(jRow);
      int idx = get_T_idx(val);
      if (idx == len) {
        list_mult.push_back(0);
        len++;
      }
      list_mult[idx]++;
    }
    for (int u=0; u<len; u++)
      if (list_mult[u] > 0) {
        esign.push_back(u);
        esign.push_back(list_mult[u]);
        list_mult[u] = 0;
      }
    int idx_sign = get_Tvs_idx(esign);
    esign.clear();
    V[iBreak] = idx_sign;
  }
  size_t n_block = ListPossibleSignatures.size();
  if (n_block == 1)
    return false;
  if (canonically) {
    // First reordering the weights
    std::pair<std::vector<T>, std::vector<int>> rec_pair_A = GetReorderingInfoWeight(ListWeight);
    ListWeight = rec_pair_A.first;
    size_t n_Wei = ListWeight.size();
    std::vector<int> list_mult(n_Wei,0);
    const std::vector<int>& g = rec_pair_A.second;
    for (auto & esign : ListPossibleSignatures) {
      size_t n_ent = esign.size() / 2;
      for (size_t i=0; i<n_ent; i++) {
        int NewVal = g[esign[2*i]];
        int eMult = esign[1 + 2*i];
        list_mult[NewVal] = eMult;
      }
      //
      std::vector<int> newsign;
      for (size_t i=0; i<n_Wei; i++)
        if (list_mult[i] > 0) {
          newsign.push_back(i);
          newsign.push_back(list_mult[i]);
        }
      esign = std::move(newsign);
    }
    // Now reordering the signatures
    std::pair<std::vector<std::vector<int>>, std::vector<int>> rec_pair_B = GetReorderingInfoWeight(ListPossibleSignatures);
    ListPossibleSignatures = rec_pair_B.first;
    const std::vector<int>& h = rec_pair_B.second;
    for (size_t iBreak=0; iBreak<siz_block_break; iBreak++) {
      int NewIdx = h[V[iBreak]];
      V[iBreak] = NewIdx;
    }
  }
  int n_block_orig = VP.ListBlocks.size();
  std::vector<int> Vmap(n_block);
  Vmap[0] = jBlock;
  for (size_t i=1; i<n_block; i++) {
    int pos = n_block_orig + i - 1;
    Vmap[i] = pos;
  }
  VP.ListBlocks[jBlock].clear();
  for (size_t i=1; i<n_block; i++) {
    VP.ListBlocks.push_back(std::vector<Tidx>());
  }
  for (size_t iBreak=0; iBreak<siz_block_break; iBreak++) {
    int iRow = eBlockBreak[iBreak];
    int iBlockLoc = V[iBreak];
    int iBlockGlob = Vmap[iBlockLoc];
    VP.MapVertexBlock[iRow] = iBlockGlob;
    VP.ListBlocks[iBlockGlob].push_back(iRow);
  }
  return true;
}

template<typename Tidx>
void PrintVertexParttionInfo(const VertexPartition<Tidx>& VP, const std::vector<uint8_t>& status)
{
  std::cerr << "nbRow=" << VP.nbRow << "\n";
  size_t n_block = VP.ListBlocks.size();
  for (size_t i=0; i<n_block; i++) {
    std::cerr << "i=" << i << "/" << n_block << " siz=" << VP.ListBlocks[i].size() << " status=" << int(status[i]) << "\n";
  }
}

template<typename T, typename Tidx, typename F1, typename F2>
VertexPartition<Tidx> ComputeVertexPartition(size_t nbRow, F1 f1, F2 f2, bool canonically, size_t max_globiter)
{
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "ComputeVertexPartition : We have nbRow=" << nbRow << " which is larger than the possible values of Tidx : " << max_poss_rows << "\n";
    throw TerminalException{1};
  }
  using Ttime = std::chrono::time_point<std::chrono::system_clock>;
  Ttime time_vp1 = std::chrono::system_clock::now();
  VertexPartition<Tidx> VP = ComputeInitialVertexPartition<T,Tidx>(nbRow, f1, f2, canonically);
  Ttime time_vp2 = std::chrono::system_clock::now();
  uint64_t dur_vp = std::chrono::duration_cast<std::chrono::microseconds>(time_vp2 - time_vp1).count();
#ifdef TIMINGS
  std::cerr << "|ComputeInitialVertexPartition|=" << dur_vp << "\n";
#endif
  std::vector<uint8_t> status(VP.ListBlocks.size(), 0);
#ifdef DEBUG_SPECIFIC
  std::cerr << "After ComputeInitialVertexPartition\n";
  PrintVertexParttionInfo(VP, status);
#endif
  auto GetPreferable_iBlock=[&]() -> int {
    size_t min_size = nbRow + 1;
    int iBlockSel = -1;
    int n_block = VP.ListBlocks.size();
    for (int iBlock=0; iBlock<n_block; iBlock++) {
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
  auto DoRefinement=[&](const int& iBlock) -> bool {
    int n_block = VP.ListBlocks.size();
    // First looking at the diagonal
    bool test1 = RefineSpecificVertexPartition<T,Tidx>(VP, iBlock, iBlock, f1, f2, canonically);
#ifdef DEBUG_SPECIFIC
    std::cerr << "iBlock=" << iBlock << " test1=" << test1 << "\n";
#endif
    if (test1) {
      size_t len1 = VP.ListBlocks.size();
      size_t len2 = status.size();
      for (size_t i=len2; i<len1; i++)
        status.push_back(0);
#ifdef DEBUG_SPECIFIC
      std::cerr << "len1=" << len1 << " len2=" << len2 << "\n";
#endif
      return true;
    }
    status[iBlock] =1;
    bool DidSomething = false;
    for (int jBlock=0; jBlock<n_block; jBlock++) {
      if (iBlock != jBlock) {
        bool test2 = RefineSpecificVertexPartition<T,Tidx>(VP, jBlock, iBlock, f1, f2, canonically);
#ifdef DEBUG_SPECIFIC
        std::cerr << "iBlock=" << iBlock << " jBlock=" << jBlock << " test2=" << test2 << "\n";
#endif
        if (test2) {
          status[jBlock] = 0;
          size_t len1 = VP.ListBlocks.size();
          size_t len2 = status.size();
          for (size_t i=len2; i<len1; i++)
            status.push_back(0);
          DidSomething = true;
        }
      }
    }
    return DidSomething;
  };
  while(true) {
    int iBlock = GetPreferable_iBlock();
    if (iBlock == -1)
      break;
#ifdef DEBUG_SPECIFIC
    std::cerr << "iBlock=" << iBlock << "\n";
#endif
    Ttime time_ref1 = std::chrono::system_clock::now();
    bool test = DoRefinement(iBlock);
    Ttime time_ref2 = std::chrono::system_clock::now();
    uint64_t dur_ref = std::chrono::duration_cast<std::chrono::microseconds>(time_ref2 - time_ref1).count();
#ifdef TIMINGS
    std::cerr << "|DoRefinement|=" << dur_ref << "\n";
#endif
#ifdef DEBUG_SPECIFIC
    std::cerr << "After Dorefinement\n";
    PrintVertexParttionInfo(VP, status);
    std::cerr << "test=" << test << "\n";
#endif
    if (!test && dur_ref > 1000 && dur_ref > 10 * dur_vp)
      break;
  }
  return VP;
}





// We use stable_sort to make sure that from the canonical ordering of the blocks
// we get a canonical ordering for the block sizes.
template<typename Tidx>
std::vector<int> GetOrdering_ListIdx(const VertexPartition<Tidx>& VP)
{
  size_t nbCase = VP.ListBlocks.size();
  std::vector<int> ListIdx(nbCase);
  for (size_t iCase=0; iCase<nbCase; iCase++)
    ListIdx[iCase] = iCase;
  std::stable_sort(ListIdx.begin(), ListIdx.end(),
                   [&](int idx1, int idx2) -> bool {
                     return VP.ListBlocks[idx1].size() < VP.ListBlocks[idx2].size();});
  return ListIdx;
}










//
// The computation of vertex signatures is needed for the computation of vertex degrees
// which are needed for TRACES graph construction. Therefore said computation cannot be
// avoided when building our objects. However, for breaking down the vertex set into
// blocks, they are too expensive
//


template<typename T>
struct WeightMatrixVertexSignatures {
  size_t nbRow;
  size_t nbWeight;
  std::vector<T> ListWeight;
  std::vector<SignVertex> ListPossibleSignatures;
  std::vector<int> ListSignatureByVertex;
  std::vector<int> ListNbCase;
};


template<typename T>
std::vector<int> GetOrdering_ListIdx(WeightMatrixVertexSignatures<T> const& WMVS)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  size_t nbCase = WMVS.ListNbCase.size();
  std::vector<int> ListIdx(nbCase);
  for (size_t iCase=0; iCase<nbCase; iCase++)
    ListIdx[iCase] = iCase;
  auto fctComp=[](auto val1, auto val2) -> std::pair<bool,bool> {
    if (val1 < val2)
      return {true, true};
    if (val1 > val2)
      return {true, false};
    return {false,true};
  };
  std::sort(ListIdx.begin(), ListIdx.end(),
            [&](int idx1, int idx2) -> bool {
              // First selection by the number of cases.
              std::pair<bool,bool> test1 = fctComp(WMVS.ListNbCase[idx1], WMVS.ListNbCase[idx2]);
              if (test1.first)
                return test1.second;
              // The cases with high number of cases are preferable.
              size_t len1 = WMVS.ListPossibleSignatures[idx1].size();
              size_t len2 = WMVS.ListPossibleSignatures[idx2].size();
              std::pair<bool,bool> test2 = fctComp(len2, len1);
              if (test2.first)
                return test2.second;
              // Now the order does not really matter for speed but it has to be
              // fully unicized for the canonical form.
              // Now going after the diagonal value
              const std::vector<int>& list_pair1 = WMVS.ListPossibleSignatures[idx1];
              const std::vector<int>& list_pair2 = WMVS.ListPossibleSignatures[idx2];
              for (size_t i=0; i<len1; i++) {
                std::pair<bool,bool> test4 = fctComp(list_pair1[i], list_pair2[i]);
                if (test4.first)
                  return test4.second;
              }
              return false;
            });
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetOrdering_ListIdx|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  return ListIdx;
}



template<typename T>
void PrintWMVS(std::ostream & os, WeightMatrixVertexSignatures<T> const& WMVS)
{
  size_t nbCase = WMVS.ListNbCase.size();
  os << "nbRow=" << WMVS.nbRow << " nbWeight=" << WMVS.nbWeight << " nbCase=" << nbCase << "\n";
  os << "ListWeight =";
  for (size_t iWei=0; iWei<WMVS.nbWeight; iWei++) {
    os << " (" << iWei << "," << WMVS.ListWeight[iWei] << ")";
  }
  os << "\n";
  //
  for (size_t iCase=0; iCase<nbCase; iCase++) {
    SignVertex eCase = WMVS.ListPossibleSignatures[iCase];
    os << "iCase=" << iCase << "/" << nbCase << " nb=" << WMVS.ListNbCase[iCase] << " eCase=" << eCase[0];
    size_t len = eCase.size() / 2;
    os << " LV=";
    for (size_t i=0; i<len; i++) {
      int first = eCase[1 + 2*i];
      int second = eCase[2 + 2*i];
      os << " [" << first << "," << second << "]";
    }
    os << "\n";
  }
}



template<typename T, typename F1, typename F2>
WeightMatrixVertexSignatures<T> ComputeVertexSignatures(size_t nbRow, F1 f1, F2 f2)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  UNORD_MAP_SPECIFIC<T,int> ValueMap_T;
  std::vector<T> ListWeight;
  int idxWeight = 0;
  auto get_T_idx=[&](const T& eval) -> int {
    int& idx = ValueMap_T[eval];
    if (idx == 0) { // value is missing
      idxWeight++;
      idx = idxWeight;
      ListWeight.push_back(eval);
    }
    return idx - 1;
  };
  UNORD_MAP_SPECIFIC<SignVertex,int> ValueMap_Tvs;
  std::vector<SignVertex> ListPossibleSignatures;
  int idxSign = 0;
  auto get_Tvs_idx=[&](SignVertex const& esign) -> int {
    int& idx = ValueMap_Tvs[esign];
    if (idx == 0) { // value is missing
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
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    f1(iRow);
    int idx_specific;
    for (size_t jRow=0; jRow<nbRow; jRow++) {
      T val = f2(jRow);
      int idx = get_T_idx(val);
      if (idx == len) {
        list_mult.push_back(0);
        len++;
      }
      if (iRow != jRow) {
        list_mult[idx]++;
      } else {
        idx_specific = idx;
      }
    }
    esign.push_back(idx_specific);
    for (int u=0; u<len; u++)
      if (list_mult[u] > 0) {
        esign.push_back(u);
        esign.push_back(list_mult[u]);
        list_mult[u] = 0;
      }
    int idx_sign = get_Tvs_idx(esign);
    esign.clear();
    ListSignatureByVertex[iRow] = idx_sign;
  }
  size_t nbWeight = ListWeight.size();
  //
  size_t nbCase = ListPossibleSignatures.size();
  std::vector<int> ListNbCase(nbCase, 0);
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    int iCase = ListSignatureByVertex[iRow];
    ListNbCase[iCase]++;
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|ComputeVertexSignature|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  return {nbRow, nbWeight, std::move(ListWeight), std::move(ListPossibleSignatures), std::move(ListSignatureByVertex), std::move(ListNbCase)};
}


template<typename T>
void RenormalizeWMVS(WeightMatrixVertexSignatures<T>& WMVS)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  // Building the permutation on the weights
  std::pair<std::vector<T>, std::vector<int>> rec_pair = GetReorderingInfoWeight(WMVS.ListWeight);
  size_t n_Wei = WMVS.ListWeight.size();
  const std::vector<int>& g = rec_pair.second;
  WMVS.ListWeight = rec_pair.first;
  // Changing the list of signatures
  std::vector<SignVertex> NewListPossibleSignatures;
  for (auto & ePossSignature : WMVS.ListPossibleSignatures) {
    int NewDiag = g[ePossSignature[0]];
    size_t len = ePossSignature.size() / 2;
    std::vector<int> list_mult(n_Wei, 0);
    for (size_t i=0; i<len; i++) {
      int NewVal = g[ePossSignature[1 + 2*i]];
      int eMult = ePossSignature[2 + 2*i];
      list_mult[NewVal] = eMult;
    }
    //
    std::vector<int> newsign;
    newsign.push_back(NewDiag);
    for (size_t i=0; i<n_Wei; i++)
      if (list_mult[i] > 0) {
        newsign.push_back(i);
        newsign.push_back(list_mult[i]);
      }
    NewListPossibleSignatures.push_back(newsign);
  }
  WMVS.ListPossibleSignatures=NewListPossibleSignatures;
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|RenormalizeWMVS|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
}



template<typename T, typename F1, typename F2>
DataTraces GetDataTraces(F1 f1, F2 f2, WeightMatrixVertexSignatures<T> const& WMVS)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  size_t nbRow = WMVS.nbRow;
  size_t nbWeight = WMVS.nbWeight;
  std::vector<T> const& ListWeight = WMVS.ListWeight;
  std::unordered_map<T,int> map;
  for (size_t iWei=0; iWei<nbWeight; iWei++) {
    map[ListWeight[iWei]] = iWei;
  }
  //
  size_t nbMult=nbWeight + 2;
  size_t hS=Pairs_GetNeededN(nbMult);
  size_t nbVert=nbRow + 2;
  size_t nbVertTot=nbVert * hS;
#ifdef DEBUG_SPECIFIC
  std::cerr << "nbMult=" << nbMult << " hS=" << hS << " nbRow=" << nbRow << " nbVert=" << nbVert << " nbVertTot=" << nbVertTot << "\n";
#endif
  //
  std::vector<int> V = Pairs_GetListPair(hS, nbMult);
  size_t e_pow = V.size() / 2;
  //
  // Now building the list of possibilities for the 2 added vertices
  //
  std::vector<std::vector<std::pair<int,int>>> new_list_signature;
  for (auto & esign : WMVS.ListPossibleSignatures) {
    std::vector<std::pair<int,int>> e_vect;
    size_t len = esign.size() / 2;
    for (size_t u=0; u<len; u++)
      e_vect.push_back({esign[1 + 2*u], esign[2 + 2*u]});
    auto fInsert=[&](int e_val) -> void {
      for (auto & e_pair : e_vect) {
        if (e_pair.first == e_val) {
          e_pair.second++;
          return;
        }
      }
      e_vect.push_back({e_val,1});
    };
    fInsert(esign[0]);
    fInsert(nbWeight + 1);
    new_list_signature.push_back(e_vect);
  }
  // nbRow : Adding the column corresponding to the diagonal values
  std::map<int,int> map_vert_nRow;
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    int isign = WMVS.ListSignatureByVertex[iRow];
    int iweight_specific = WMVS.ListPossibleSignatures[isign][0];
    map_vert_nRow[iweight_specific]++;
  }
  map_vert_nRow[nbWeight]++;
  std::vector<std::pair<int,int>> f_vect;
  for (auto & kv : map_vert_nRow)
    f_vect.push_back({kv.first, kv.second});
  new_list_signature.push_back(f_vect);
  // nbRow+1 : Adding the column corresponding to the external vertex
  std::vector<std::pair<int,int>> g_vect { {nbWeight,1} , {nbWeight+1, nbRow}};
  new_list_signature.push_back(g_vect);
  size_t nbCase = new_list_signature.size();
  // Now adding the list_of index
  std::vector<int> list_signature = WMVS.ListSignatureByVertex;
  list_signature.push_back(nbCase-2); // for nbRow
  list_signature.push_back(nbCase-1); // for nbRow+1
  //
  // Determining the number of cases
  //
  std::vector<int> ListNbCase(nbCase,0);
  for (size_t i=0; i<nbVert; i++) {
    int iCase = list_signature[i];
    ListNbCase[iCase]++;
  }
  //
  // Now computing the weight by multiplier
  //
  std::vector<size_t> WeightByMult(nbMult, 0);
  for (size_t iCase=0; iCase<nbCase; iCase++) {
    int sizCase = ListNbCase[iCase];
    std::vector<std::pair<int,int>> e_vect = new_list_signature[iCase];
    for (auto & e_pair : e_vect) {
      int iWeight = e_pair.first;
      int eMult = e_pair.second;
      WeightByMult[iWeight] += sizCase * eMult;
    }
  }
  //
  // Reordering the multiplier to maximize the sparsity
  // (This is not perfect since if iH1 != iH2 we need to add more edges. See below)
  // (There are more strange things happening. We should have gotten some improvement accross
  // the board in memory usage, but we do not)
  //
  std::vector<size_t> ListIdx(nbMult);
  for (size_t iMult=0; iMult<nbMult; iMult++)
    ListIdx[iMult] = iMult;
  std::stable_sort(ListIdx.begin(), ListIdx.end(), [&](const size_t& idx1, const size_t& idx2) -> bool {
    return WeightByMult[idx1] > WeightByMult[idx2];
  });
  //
  // Now computing the list of signature
  //
  MyMatrix<int> MatrixAdj = ZeroMatrix<int>(hS, nbCase);
  // Adjacencies (iVert + nbVert*iH , iVert + nbVert*jH)
  size_t nb_adj1 = hS - 1;
  for (size_t iH=0; iH<hS; iH++)
    for (size_t iCase=0; iCase<nbCase; iCase++)
      MatrixAdj(iH, iCase) += nb_adj1;
  //
  // Now computing the contribution of the adjacencies
  //
  Face f_total = GetAllBinaryExpressionsByWeight(e_pow, nbMult);
  for (size_t iCase=0; iCase<nbCase; iCase++) {
    std::vector<std::pair<int,int>> e_vect = new_list_signature[iCase];
    for (auto & e_pair : e_vect) {
      int iWeight = ListIdx[e_pair.first];
      int eMult = e_pair.second;
      size_t shift = iWeight * e_pow;
      for (size_t i_pow=0; i_pow<e_pow; i_pow++) {
        if (f_total[shift + i_pow] == 1) {
          int iH1 = V[2*i_pow];
          int iH2 = V[2*i_pow+1];
          // Adjacency (iVert + nbVert*iH , jVert + nbVert*iH2)
          MatrixAdj(iH1, iCase) += eMult;
          if (iH1 != iH2) {
            // Adjacency (iVert + nbVert*iH2 , jVert + nbVert*iH1)
            MatrixAdj(iH2, iCase) += eMult;
          }
        }
      }
    }
  }
  //
  // Now building the adjacencies for traces
  //
  int nbAdjacent = 0;
  for (size_t iCase=0; iCase<nbCase; iCase++)
    for (size_t iH=0; iH<hS; iH++)
      nbAdjacent += ListNbCase[iCase] * MatrixAdj(iH, iCase);
  DataTraces DT(nbVertTot, nbAdjacent);
  // Now setting up the adjacencies
  int pos=0;
  std::vector<int> ListShift(nbVertTot);
  for (size_t i=0; i<nbVertTot; i++) {
    size_t iVert = i % nbVert;
    size_t iH = i / nbVert;
    int iCase = list_signature[iVert];
    int nbAdj = MatrixAdj(iH, iCase);
    DT.sg1.d[i] = nbAdj;
    DT.sg1.v[i] = pos;
    ListShift[i] = pos;
    pos += nbAdj;
  }
  // Assigning the color stuff
  for (size_t i=0; i<nbVertTot; i++) {
    DT.lab1[i] = i;
  }
  for (size_t i=0; i<nbVertTot; i++) {
    DT.ptn[i] = NAUTY_INFINITY;
  }
  for (size_t iH=0; iH<hS; iH++) {
    size_t pos = iH * nbVert + nbVert - 1;
    DT.ptn[pos] = 0;
  }
  //
#ifdef DEBUG_SPECIFIC
  std::vector<int> ListDegExpe1(nbVertTot,0);
  std::vector<int> ListDegExpe2(nbVertTot,0);
#endif
  auto f_adj=[&](size_t iVert, size_t jVert) -> void {
    DT.sg1.e[ListShift[iVert]] = jVert;
    ListShift[iVert]++;
#ifdef DEBUG_SPECIFIC
    ListDegExpe1[iVert]++;
    ListDegExpe2[jVert]++;
#endif
  };
  for (size_t iVert=0; iVert<nbVert; iVert++)
    for (size_t iH=0; iH<hS-1; iH++)
      for (size_t jH=iH+1; jH<hS; jH++) {
        size_t aVert=iVert + nbVert*iH;
        size_t bVert=iVert + nbVert*jH;
        f_adj(aVert, bVert);
        f_adj(bVert, aVert);
      }
  for (size_t iVert=0; iVert<nbVert-1; iVert++) {
    if (iVert < nbRow)
      f1(iVert);
    for (size_t jVert=iVert+1; jVert<nbVert; jVert++) {
      int eVal;
      if (jVert == nbRow+1) {
        if (iVert == nbRow)
          eVal = nbWeight;
        else
          eVal = nbWeight + 1;
      } else {
        if (jVert == nbRow) {
          eVal = map.at(f2(iVert));
        } else {
          eVal = map.at(f2(jVert));
        }
      }
      size_t shift = e_pow * ListIdx[eVal];
      for (size_t i_pow=0; i_pow<e_pow; i_pow++)
        if (f_total[shift + i_pow] == 1) {
          int iH1 = V[2*i_pow];
          int iH2 = V[2*i_pow+1];
          size_t aVert=iVert + nbVert*iH1;
          size_t bVert=jVert + nbVert*iH2;
          f_adj(aVert, bVert);
          f_adj(bVert, aVert);
          if (iH1 != iH2) {
            aVert=iVert + nbVert*iH2;
            bVert=jVert + nbVert*iH1;
            f_adj(aVert, bVert);
            f_adj(bVert, aVert);
          }
        }
    }
  }
#ifdef DEBUG_SPECIFIC
  int sum_adj = 0;
  int nb_error = 0;
  for (size_t iVert=0; iVert<nbVertTot; iVert++) {
    int deg0 = DT.sg1.d[iVert];
    int deg1 = ListDegExpe1[iVert];
    int deg2 = ListDegExpe2[iVert];
    if (deg0 != deg1 || deg0 != deg2 || deg1 != deg2) {
      std::cerr << "iVert=" << iVert << " deg0=" << DT.sg1.d[iVert] << " deg1=" << ListDegExpe1[iVert] << " deg2=" << ListDegExpe2[iVert] << "\n";
      nb_error++;
    }
    sum_adj += ListDegExpe2[iVert];
  }
  double frac_adj = double(sum_adj) / (double(nbVertTot) * double(nbVertTot-1));
  std::cerr << "sum_adj=" << sum_adj << " nbAdjacent=" << nbAdjacent << " frac_adj=" << frac_adj << "\n";
  if (nb_error > 0) {
    std::cerr << "nb_error=" << nb_error << "\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetDataTraces|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  return DT;
}






template<typename T, typename Tidx, typename F1, typename F2>
std::vector<Tidx> GetCanonicalizationVector_KnownSignature(WeightMatrixVertexSignatures<T> const& WMVS, F1 f1, F2 f2)
{
  size_t nbRow = WMVS.nbRow;
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "GetCanonicalizationVector_KnownSignature : We have nbRow=" << nbRow << " which is larger than the possible values of Tidx : " << max_poss_rows << "\n";
    throw TerminalException{1};
  }
  DataTraces DT = GetDataTraces<T,F1,F2>(f1, f2, WMVS);
  if (DT.n < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    using TidxIn = uint8_t;
    std::vector<TidxIn> cl = TRACES_GetCanonicalOrdering_Arr<TidxIn>(DT);
    return GetCanonicalizationVector_KernelBis<Tidx,TidxIn>(nbRow, cl);
  }
  if (DT.n < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    using TidxIn = uint16_t;
    std::vector<TidxIn> cl = TRACES_GetCanonicalOrdering_Arr<TidxIn>(DT);
    return GetCanonicalizationVector_KernelBis<Tidx,TidxIn>(nbRow, cl);
  }
  if (DT.n < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    using TidxIn = uint32_t;
    std::vector<TidxIn> cl = TRACES_GetCanonicalOrdering_Arr<TidxIn>(DT);
    return GetCanonicalizationVector_KernelBis<Tidx,TidxIn>(nbRow, cl);
  }
#if !defined __APPLE__
  if (DT.n < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    using TidxIn = uint64_t;
    std::vector<TidxIn> cl = TRACES_GetCanonicalOrdering_Arr<TidxIn>(DT);
    return GetCanonicalizationVector_KernelBis<Tidx,TidxIn>(nbRow, cl);
  }
#endif
  std::cerr << "Failed to find matching type in GetCanonicalizationVector_KnownSignature\n";
  throw TerminalException{1};
}



template<typename T, typename Tidx, typename F1, typename F2>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> GetGroupCanonicalization_KnownSignature(WeightMatrixVertexSignatures<T> const& WMVS, F1 f1, F2 f2)
{
  size_t nbRow = WMVS.nbRow;
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "GetGroupCanonicalization_KnownSignature : We have nbRow=" << nbRow << " which is larger than the possible values of Tidx : " << max_poss_rows << "\n";
    throw TerminalException{1};
  }
  DataTraces DT = GetDataTraces<T,F1,F2>(f1, f2, WMVS);
  if (DT.n < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    using TidxC = uint8_t;
    std::pair<std::vector<TidxC>, std::vector<std::vector<Tidx>>> ePair = TRACES_GetCanonicalOrdering_ListGenerators_Arr<TidxC,Tidx>(DT, nbRow);
    std::vector<Tidx> MapVectRev2 = GetCanonicalizationVector_KernelBis<Tidx,TidxC>(nbRow, ePair.first);
    return {std::move(MapVectRev2), std::move(ePair.second)};
  }
  if (DT.n < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    using TidxC = uint16_t;
    std::pair<std::vector<TidxC>, std::vector<std::vector<Tidx>>> ePair = TRACES_GetCanonicalOrdering_ListGenerators_Arr<TidxC,Tidx>(DT, nbRow);
    std::vector<Tidx> MapVectRev2 = GetCanonicalizationVector_KernelBis<Tidx,TidxC>(nbRow, ePair.first);
    return {std::move(MapVectRev2), std::move(ePair.second)};
  }
  if (DT.n < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    using TidxC = uint32_t;
    std::pair<std::vector<TidxC>, std::vector<std::vector<Tidx>>> ePair = TRACES_GetCanonicalOrdering_ListGenerators_Arr<TidxC,Tidx>(DT, nbRow);
    std::vector<Tidx> MapVectRev2 = GetCanonicalizationVector_KernelBis<Tidx,TidxC>(nbRow, ePair.first);
    return {std::move(MapVectRev2), std::move(ePair.second)};
  }
#if !defined __APPLE__
  if (DT.n < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    using TidxC = uint64_t;
    std::pair<std::vector<TidxC>, std::vector<std::vector<Tidx>>> ePair = TRACES_GetCanonicalOrdering_ListGenerators_Arr<TidxC,Tidx>(DT, nbRow);
    std::vector<Tidx> MapVectRev2 = GetCanonicalizationVector_KernelBis<Tidx,TidxC>(nbRow, ePair.first);
    return {std::move(MapVectRev2), std::move(ePair.second)};
  }
#endif
  std::cerr << "Failed to find matching numeric in GetGroupCanonicalization_KnownSignature\n";
  throw TerminalException{1};
}



template<typename T, typename Tidx, typename F1, typename F2>
std::vector<std::vector<Tidx>> GetStabilizerWeightMatrix_KnownSignature(WeightMatrixVertexSignatures<T> const& WMVS, F1 f1, F2 f2)
{
  size_t nbRow = WMVS.nbRow;
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "GetStabilizerWeightMatrix_KnownSignature : We have nbRow=" << nbRow << " which is larger than the possible values of Tidx : " << max_poss_rows << "\n";
    throw TerminalException{1};
  }
  DataTraces DT = GetDataTraces<T,F1,F2>(f1, f2, WMVS);
  return TRACES_GetListGenerators_Arr<Tidx>(DT, nbRow);
}


/*
  ---F1/F2 : The first/second template function for creating the Weight matrix
  ---F3    : The function for testing acceptability of sets for consideration (e.g. rank function)
  ---F4    : The function for testing acceptability of small generators and returning the big
             generators. It also takes blocks in input and determine their nature, that is preserved
             or not
  ---Fproc1 : function processing the graph (and using traces for this)
  ---Fproc2 : function taking the output and returning the list of generators
  ---Fproc3 : function taking all of it and returning the output.
*/
template<bool canonically, typename T, typename Tidx, typename Tret1, typename Tret2, typename Tret3, typename F1, typename F2, typename F3, typename F4, typename Fproc1, typename Fproc2, typename Fproc3>
Tret3 BlockBreakdown_Heuristic(size_t nbRow, F1 f1, F2 f2, F3 f3, F4 f4, Fproc1 fproc1, Fproc2 fproc2, Fproc3 fproc3)
{
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "BlockBreakdown_Heuristic : We have nbRow=" << nbRow << " which is larger than the possible values of Tidx : " << max_poss_rows << "\n";
    throw TerminalException{1};
  }
#ifdef DEBUG_SPECIFIC
  std::cerr << "Beginning of BlockBreakdown_Heuristic\n";
#endif
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  size_t max_globiter = 1000;
  VertexPartition<Tidx> VP = ComputeVertexPartition<T,Tidx>(nbRow, f1, f2, canonically, max_globiter);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|ComputeVertexPartition|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  size_t nbCase = VP.ListBlocks.size();
  std::vector<int> ListIdx = GetOrdering_ListIdx(VP);
#ifdef DEBUG_SPECIFIC
  for (size_t iCase=0; iCase<nbCase; iCase++) {
    int idx = ListIdx[iCase];
    std::cerr << "iCase=" << iCase << " ListIdx[iCase]=" << idx << " nbCase=" << VP.ListBlocks[idx].size() << "\n";
  }
#endif
  Face f_covered(nbCase);
  std::vector<int> StatusCase(nbCase);
  size_t idx=0;
  auto set_status_case=[&]() -> void {
    for (size_t iCase=0; iCase<nbCase; iCase++)
      StatusCase[iCase] = 0;
    for (size_t u=0; u<idx; u++) {
      size_t pos = ListIdx[u];
      if (f_covered[pos] == 0)
        StatusCase[ListIdx[u]] = 1;
    }
  };
  auto sum_status_case=[&]() -> int {
    int sum=0;
    for (size_t iCase=0; iCase<nbCase; iCase++)
      sum += StatusCase[iCase];
    return sum;
  };
  auto increase_idx=[&]() -> bool {
    int sum_prev = sum_status_case();
    while(true) {
      idx++;
      set_status_case();
      int sum_new = sum_status_case();
      if (sum_new > sum_prev)
        return true;
      if (idx == nbCase)
        return false;
    }
  };
  while(true) {
    bool test = increase_idx();
    if (!test)
      break;
    std::vector<Tidx> CurrentListIdx;
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      int iCase = VP.MapVertexBlock[iRow];
      if (StatusCase[iCase] == 1)
        CurrentListIdx.push_back(iRow);
    }
    size_t nbRow_res = CurrentListIdx.size();
#ifdef DEBUG_SPECIFIC
    std::cerr << "idx=" << idx << " |CurrentListIdx|=" << nbRow_res << "\n";
#endif
    //
    if (f3(CurrentListIdx)) {
      std::vector<std::vector<Tidx>> ListBlocks;
      std::vector<size_t> ListEnt;
      for (size_t u=0; u<nbCase; u++)
        if (StatusCase[u] == 0) {
          ListBlocks.push_back(VP.ListBlocks[u]);
          ListEnt.push_back(u);
        }
      auto f1_res = [&](size_t iRow) -> void {
        f1(CurrentListIdx[iRow]);
      };
      auto f2_res = [&](size_t jRow) -> T {
        return f2(CurrentListIdx[jRow]);
      };
      WeightMatrixVertexSignatures<T> WMVS_res = ComputeVertexSignatures<T>(nbRow_res, f1_res, f2_res);
#ifdef DEBUG_SPECIFIC
      std::cerr << "WMVS_res=\n";
      PrintWMVS(std::cerr, WMVS_res);
#endif
      Tret1 ret1 = fproc1(WMVS_res, f1_res, f2_res);
      bool IsCorrect = true;
      std::vector<std::vector<Tidx>> LGen;
      Face f_incorr(ListEnt.size());
      for (auto & eList : fproc2(ret1)) {
        DataMapping<Tidx> test = f4(CurrentListIdx, eList, ListBlocks);
        if (test.correct) {
          LGen.push_back(test.eGen);
        } else {
          IsCorrect = false;
        }
        size_t n_incorr=0;
        for (size_t u=0; u<ListEnt.size(); u++) {
          if (!test.block_status[u])
            f_incorr[u] = 1;
          n_incorr += f_incorr[u];
        }
        if (n_incorr == ListEnt.size()) // All entries are incorrect. So, no need to continue
          break;
      }
      std::cerr << "f_incorr=";
      for (size_t u=0; u<ListEnt.size(); u++)
        std::cerr << int(f_incorr[u]);
      std::cerr << "\n";

      for (size_t u=0; u<ListEnt.size(); u++) {
        size_t pos = ListEnt[u];
        if (f_incorr[u] == 0)
          f_covered[pos] = 1;
      }
      if (IsCorrect) {
        return fproc3(CurrentListIdx, ret1, LGen);
      }
    }
  }
  std::cerr << "BlockBreakdown_Heuristic : We should never reach that stage\n";
  throw TerminalException{1};
}


template<typename T, typename Tidx, typename F1, typename F2, typename F3, typename F4>
std::vector<std::vector<Tidx>> GetStabilizerWeightMatrix_Heuristic(size_t nbRow, F1 f1, F2 f2, F3 f3, F4 f4)
{
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "GetStabilizerWeightMatrix_Heuristic : We have nbRow=" << nbRow << " which is larger than the possible values of Tidx : " << max_poss_rows << "\n";
    throw TerminalException{1};
  }
#ifdef DEBUG_SPECIFIC
  std::cerr << "Beginning of GetStabilizerWeightMatrix_Heuristic\n";
#endif
  using Tret1 = std::vector<std::vector<Tidx>>;
  using Tret2 = std::vector<std::vector<Tidx>>;
  using Tret3 = std::vector<std::vector<Tidx>>;
  auto fproc1=[&](const WeightMatrixVertexSignatures<T>& WMVS_res, auto f1_res, auto f2_res) -> Tret1 {
    return GetStabilizerWeightMatrix_KnownSignature<T,Tidx>(WMVS_res, f1_res, f2_res);
  };
  auto fproc2=[&](const Tret1& ListGen) -> const Tret2& {
#ifdef DEBUG_SPECIFIC
    std::cerr << "|ListGen|=" << ListGen.size() << "\n";
#endif
    return ListGen;
  };
  auto fproc3=[&](const std::vector<Tidx>& Vsubset, const Tret1& ret1, const std::vector<std::vector<Tidx>>& LGenFinal) -> Tret3 {
    return LGenFinal;
  };
  const bool canonically=false;
  return BlockBreakdown_Heuristic<canonically,T,Tidx,Tret1,Tret2,Tret3>(nbRow, f1, f2, f3, f4, fproc1, fproc2, fproc3);
}



/*
  ---F1/F2 : The first/second template function for creating the Weight matrix
  ---F3    : The function for testing acceptability of sets for consideration (e.g. rank function)
  ---F4    : The function for testing acceptability of small generators and returning the big
             generators.
  ---F5    : The lifting of canonicalization from a subset to the full set.
*/
template<typename T, typename Tidx, typename F1, typename F2, typename F3, typename F4, typename F5>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> GetGroupCanonicalizationVector_Heuristic(size_t nbRow, F1 f1, F2 f2, F3 f3, F4 f4, F5 f5)
{
  size_t max_poss_rows = size_t(std::numeric_limits<Tidx>::max());
  if (nbRow >= max_poss_rows) {
    std::cerr << "GetGroupCanonicalizationVector_Heuristic : We have nbRow=" << nbRow << " which is larger than the possible values of Tidx : " << max_poss_rows << "\n";
    throw TerminalException{1};
  }
#ifdef DEBUG_SPECIFIC
  std::cerr << "Beginning of GetGroupCanonicalizationVector_Heuristic\n";
#endif
  using Tret1 = std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>>;
  using Tret2 = std::vector<std::vector<Tidx>>;
  using Tret3 = std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>>;
  auto fproc1=[&](const WeightMatrixVertexSignatures<T>& WMVS_res, auto f1_res, auto f2_res) -> Tret1 {
    return GetGroupCanonicalization_KnownSignature<T,Tidx>(WMVS_res, f1_res, f2_res);
  };
  auto fproc2=[&](const Tret1& ePair) -> const Tret2& {
#ifdef DEBUG_SPECIFIC
    std::cerr << "|ListGen|=" << ePair.second.size() << "\n";
#endif
    return ePair.second;
  };
  auto fproc3=[&](const std::vector<Tidx>& Vsubset, const Tret1& ePair, const std::vector<std::vector<Tidx>>& LGenFinal) -> Tret3 {
    return {f5(Vsubset, ePair.first), LGenFinal};
  };
  const bool canonically=true;
  return BlockBreakdown_Heuristic<canonically,T,Tidx,Tret1,Tret2,Tret3>(nbRow, f1, f2, f3, f4, fproc1, fproc2, fproc3);
}


#endif

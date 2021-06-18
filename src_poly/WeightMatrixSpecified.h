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



struct VertexPartition {
  std::vector<int> MapVertexBlock;
  std::vector<std::vector<int>> ListBlocks;
};


template<typename T, typename F1, typename F2>
VertexPartition ComputeInitialVertexPartition(size_t nbRow, F1 f1, F2 f2, bool canonically)
{
  std::vector<T> ListWeight;
  UNORD_MAP_SPECIFIC<T,int> ValueMap_T;
  int idxWeight = 0;
  auto get_T_idx=[&](T eval) -> int {
    int& idx = ValueMap_T[eval];
    if (idx == 0) { // value is missing
      idxWeight++;
      idx = idxWeight;
      ListWeight.push_back(eval);
    }
    return idx - 1;
  };
  std::vector<int> MapVertexBlock;
  std::vector<std::vector<int>> ListBlocks;
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    f1(iRow);
    T val = f2(iRow);
    int idx = get_T_idx(val);
    MapVertexBlock[iRow] = idx;
    if (idx < ListBlocks.size()) {
      ListBlocks[idx].push_back(iRow);
    } else {
      std::vector<int> blk{int(iRow)};
      ListBlocks.push_back(blk);
    }
  }
  if (canonically) {
    std::pair<std::vector<T>, std::vector<int>> rec_pair = GetReorderingInfoWeight(ListWeight);
    const std::vector<int>& g = rec_pair.second;
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      int NewIdx = g[MapVertexBlock[iRow]];
      MapVertexBlock[iRow] = NewIdx;
    }
    size_t n_block = ListBlocks.size();
    std::vector<std::vector<int>> NewListBlocks(n_block);
    for (size_t iBlk=0; iBlk<n_block; iBlk++) {
      size_t jBlk = g[iBlk];
      NewListBlocks[jBlk] = std::move(ListBlocks[iBlk]);
    }
    ListBlocks = std::move(NewListBlocks);
  }
  return {std::move(MapVertexBlock), std::move(ListBlocks)};
}

/*
template<typename T, typename F1, typename F2>
VertexPartition ComputeVertexPartition(size_t nbRow, F1 f1, F2 f2, bool canonically, int max_globiter)
{
  
}
*/


//
// The computation of vertex signatures is needed for the computation of vertex degrees
// which are needed for TRACES graph construction. Therefore said computation cannot be
// avoided when building our objects. However, for breaking down the vertex set into
// blocks, they are too expensive
//


using SignVertex = std::vector<int>;


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
  auto get_T_idx=[&](T eval) -> int {
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
    for (size_t i=0; i<len; i++)
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



template<typename T, typename Tidx, typename F1, typename F2>
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
  std::vector<int> eVect(e_pow);
  for (size_t iCase=0; iCase<nbCase; iCase++) {
    std::vector<std::pair<int,int>> e_vect = new_list_signature[iCase];
    for (auto & e_pair : e_vect) {
      int iWeight = e_pair.first;
      int eMult = e_pair.second;
      GetBinaryExpression(iWeight, e_pow, eVect);
      for (size_t i_pow=0; i_pow<e_pow; i_pow++) {
        if (eVect[i_pow] == 1) {
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
  std::vector<int> ListNbCase(nbCase,0);
  for (size_t i=0; i<nbVert; i++) {
    int iCase = list_signature[i];
    ListNbCase[iCase]++;
  }
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
          eVal = map[f2(iVert)];
        } else {
          eVal = map[f2(jVert)];
        }
      }
      GetBinaryExpression(eVal, e_pow, eVect);
      for (size_t i_pow=0; i_pow<e_pow; i_pow++)
        if (eVect[i_pow] == 1) {
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
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  int nbRow = WMVS.nbRow;
  DataTraces DT = GetDataTraces<T,Tidx>(f1, f2, WMVS);
  std::vector<Tidx> cl = TRACES_GetCanonicalOrdering_Arr<Tidx>(DT);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetCanonicalizationVector_KnownSignature|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  return GetCanonicalizationVector_KernelBis<Tidx>(nbRow, cl);
}



template<typename T, typename Tidx, typename F1, typename F2>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> GetGroupCanonicalization_KnownSignature(WeightMatrixVertexSignatures<T> const& WMVS, F1 f1, F2 f2)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  int nbRow = WMVS.nbRow;
  DataTraces DT = GetDataTraces<T,Tidx>(f1, f2, WMVS);
  std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> ePair = TRACES_GetCanonicalOrdering_ListGenerators_Arr<Tidx>(DT, nbRow);
  std::vector<std::vector<Tidx>> LGen;
  for (auto& eListI : ePair.second) {
    std::vector<Tidx> eListO(nbRow);
    for (Tidx i=0; i<Tidx(nbRow); i++)
      eListO[i] = eListI[i];
    LGen.push_back(eListO);
  }
  //
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
#endif
  std::vector<Tidx> MapVectRev2 = GetCanonicalizationVector_KernelBis<Tidx>(nbRow, ePair.first);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
  std::cerr << "|GetDataTraces + CanonicalListGen|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
  std::cerr << "|GetCanonicalizationVector_KernelBis|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
#endif
  return {std::move(MapVectRev2), std::move(LGen)};
}



template<typename T, typename Tidx, typename F1, typename F2>
std::vector<std::vector<Tidx>> GetStabilizerWeightMatrix_KnownSignature(WeightMatrixVertexSignatures<T> const& WMVS, F1 f1, F2 f2)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  int nbRow = WMVS.nbRow;
  DataTraces DT = GetDataTraces<T,Tidx>(f1, f2, WMVS);
  std::vector<std::vector<Tidx>> ListGen = TRACES_GetListGenerators_Arr<Tidx>(DT, nbRow);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetDataTraces + GetListGenerators|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  return ListGen;
}


/*
  ---F1/F2 : The first/second template function for creating the Weight matrix
  ---F3    : The function for testing acceptability of sets for consideration (e.g. rank function)
  ---F4    : The function for testing acceptability of small generators and returning the big
             generators.
*/
template<typename T, typename Tidx, typename F1, typename F2, typename F3, typename F4>
std::vector<std::vector<Tidx>> GetStabilizerWeightMatrix_Heuristic(size_t nbRow, F1 f1, F2 f2, F3 f3, F4 f4)
{
#ifdef DEBUG_SPECIFIC
  std::cerr << "Beginning of GetStabilizerWeightMatrix_Heuristic\n";
#endif

  WeightMatrixVertexSignatures<T> WMVS = ComputeVertexSignatures<T>(nbRow, f1, f2);
#ifdef DEBUG_SPECIFIC
  PrintWMVS(std::cerr, WMVS);
#endif
  size_t nbCase = WMVS.ListPossibleSignatures.size();
  std::vector<int> ListIdx = GetOrdering_ListIdx(WMVS);
#ifdef DEBUG_SPECIFIC
  for (size_t iCase=0; iCase<nbCase; iCase++) {
    int idx = ListIdx[iCase];
    std::cerr << "iCase=" << iCase << " ListIdx[iCase]=" << idx << " nbCase=" << WMVS.ListNbCase[idx] << "\n";
  }
#endif

  for (size_t idx=1; idx<=nbCase; idx++) {
    std::vector<int> StatusCase(nbCase,0);
    for (size_t u=0; u<idx; u++)
      StatusCase[ListIdx[u]] = 1;
    std::vector<Tidx> CurrentListIdx;
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      int iCase = WMVS.ListSignatureByVertex[iRow];
      if (StatusCase[iCase] == 1)
        CurrentListIdx.push_back(iRow);
    }
    size_t nbRow_res = CurrentListIdx.size();
#ifdef DEBUG_SPECIFIC
    std::cerr << "idx=" << idx << " |CurrentListIdx|=" << nbRow_res << "\n";
#endif
    //
    if (f3(CurrentListIdx)) {
      auto f1_res = [&](size_t iRow) -> void {
        f1(CurrentListIdx[iRow]);
      };
      auto f2_res = [&](size_t jRow) -> T {
        return f2(CurrentListIdx[jRow]);
      };
      WeightMatrixVertexSignatures<T> WMVS_res = ComputeVertexSignatures<T>(nbRow_res, f1_res, f2_res);
      std::vector<std::vector<Tidx>> ListGen = GetStabilizerWeightMatrix_KnownSignature<T,Tidx>(WMVS_res, f1_res, f2_res);
#ifdef DEBUG_SPECIFIC
      std::cerr << "|ListGen|=" << ListGen.size() << "\n";
#endif
      bool IsCorrect = true;
      std::vector<std::vector<Tidx>> LGen;
      for (auto & eList : ListGen) {
        if (IsCorrect) {
          EquivTest<std::vector<Tidx>> test = f4(CurrentListIdx, eList);
          if (test.TheReply) {
            LGen.push_back(test.TheEquiv);
          } else {
            IsCorrect = false;
          }
        }
      }
      if (IsCorrect) {
        return LGen;
      }
    }
  }
  std::cerr << "GetStabilizerWeightMatrix_Heuristic : We should never reach that stage\n";
  throw TerminalException{1};
}






/*
  ---F1/F2 : The first/second template function for creating the Weight matrix
  ---F3    : The function for testing acceptability of sets for consideration (e.g. rank function)
  ---F4    : The function for testing acceptability of small generators and returning the big
             generators.
  ---F5    : The lifting of canonicalization from a subset to the full set.
*/
template<typename T, typename Tidx, typename F1, typename F2, typename F3, typename F4, typename F5>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>>  GetGroupCanonicalizationVector_Heuristic(size_t nbRow, F1 f1, F2 f2, F3 f3, F4 f4, F5 f5)
{
  WeightMatrixVertexSignatures<T> WMVS = ComputeVertexSignatures<T>(nbRow, f1, f2);
  RenormalizeWMVS(WMVS);
#ifdef DEBUG_SPECIFIC
  std::cerr << "WMVS=\n";
  PrintWMVS(std::cerr, WMVS);
#endif
  size_t nbCase = WMVS.ListPossibleSignatures.size();
  std::vector<int> ListIdx = GetOrdering_ListIdx(WMVS);
  for (size_t idx=1; idx<=nbCase; idx++) {
    std::vector<int> StatusCase(nbCase,0);
    for (size_t u=0; u<idx; u++)
      StatusCase[ListIdx[u]] = 1;
    std::vector<Tidx> CurrentListIdx;
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      int iCase = WMVS.ListSignatureByVertex[iRow];
      if (StatusCase[iCase] == 1)
        CurrentListIdx.push_back(iRow);
    }
    size_t nbRow_res = CurrentListIdx.size();
#ifdef DEBUG_SPECIFIC
    std::cerr << "|CurrentListIdx|=" << nbRow_res << "\n";
#endif
    //
    if (f3(CurrentListIdx)) {
      auto f1_res = [&](size_t iRow) -> void {
        f1(CurrentListIdx[iRow]);
      };
      auto f2_res = [&](size_t jRow) -> T {
        return f2(CurrentListIdx[jRow]);
      };
      WeightMatrixVertexSignatures<T> WMVS_res = ComputeVertexSignatures<T>(nbRow_res, f1_res, f2_res);
      RenormalizeWMVS(WMVS_res);
#ifdef DEBUG_SPECIFIC
      std::cerr << "WMVS_res=\n";
      PrintWMVS(std::cerr, WMVS_res);
#endif
      std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> ePair = GetGroupCanonicalization_KnownSignature<T,Tidx>(WMVS_res, f1_res, f2_res);
      bool IsCorrect = true;
      std::vector<std::vector<Tidx>> LGen;
      for (auto & eList : ePair.second) {
        if (IsCorrect) {
          EquivTest<std::vector<Tidx>> test = f4(CurrentListIdx, eList);
          if (test.TheReply) {
            LGen.push_back(test.TheEquiv);
          } else {
            IsCorrect = false;
          }
        }
      }
      if (IsCorrect) {
        return {f5(CurrentListIdx, ePair.first), LGen};
      }
    }
  }
  std::cerr << "GetGroupCanonicalizationVector_Heuristic : We should never reach that stage\n";
  throw TerminalException{1};
}


#endif

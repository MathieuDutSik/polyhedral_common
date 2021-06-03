#ifndef INCLUDE_WEIGHT_MATRIX_SPECIFIED_H
#define INCLUDE_WEIGHT_MATRIX_SPECIFIED_H




template<typename T>
struct WeightMatrixVertexSignatures {
  size_t nbRow;
  size_t nbWeight;
  std::vector<T> ListWeight;
  std::vector<std::pair<int, std::vector<std::pair<int,int>>>> ListPossibleSignatures;
  std::vector<int> ListSignatureByVertex;
  std::vector<int> ListNbCase;
};



template<typename T, typename F1, typename F2>
WeightMatrixVertexSignatures<T> ComputeVertexSignatures(size_t nbRow, F1 f1, F2 f2)
{
  std::unordered_map<T,int> ValueMap_T;
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

  using Tvertex_signature = std::pair<int, std::vector<std::pair<int,int>>>;
  std::unordered_map<Tvertex_signature,int> ValueMap_Tvs;
  std::vector<Tvertex_signature> ListPossibleSignatures;
  int idxSign = 0;
  auto get_Tvs_idx=[&](Tvertex_signature const& esign) -> int {
    int& idx = ValueMap_Tvs[esign];
    if (idx == 0) { // value is missing
      idxSign++;
      idx = idxSign;
      ListPossibleSignatures.push_back(esign);
    }
    return idxSign - 1;
  };
  std::vector<int> ListSignatureByVertex(nbRow);
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    f1(iRow);
    std::map<int,int> map_mult;
    int idx_specific;
    for (size_t jRow=0; jRow<nbRow; jRow++) {
      T val = f2(jRow);
      int idx = get_T_idx(val);
      if (iRow != jRow) {
        map_mult[idx]++;
      } else {
        idx_specific = idx;
      }
    }
    std::vector<std::pair<int,int>> list_pair;
    for (auto & kv : map_mult)
      list_pair.push_back({kv.first, kv.second});
    std::pair<int, std::vector<std::pair<int,int>>> e_pair{idx_specific, list_pair};
    int idx_sign = get_Tvs_idx(e_pair);
    ListSignatureByVertex[iRow] = idx_sign;
  }
  size_t nbWeight = ListWeight.size();
  return {nbRow, nbWeight, std::move(ListWeight), std::move(ListPossibleSignatures), std::move(ListSignatureByVertex)};
}


template<typename T>
void RenormalizeWMVS(WeightMatrixVertexSignatures<T>& WMVS)
{
  // Building the permutation on the weights
  std::map<T, int> map;
  size_t n_Wei = WMVS.nbWeight;
  for (size_t i=0; i<n_Wei; i++)
    map[WMVS.ListWeight[i]] = i;
  std::vector<int> g(n_Wei);
  int idx=0;
  for (auto & kv : map) {
    int pos = kv.second;
    g[pos] = idx;
    idx++;
  }
  // Changing the weight
  std::vector<T> NewListWeight(n_Wei);
  for (size_t iW=0; iW<n_Wei; iW++) {
    size_t jW = g[iW];
    NewListWeight[jW] = WMVS.ListWeight[iW];
  }
  WMVS.ListWeight = NewListWeight;
  // Changing the list of signatures
  std::vector<std::pair<int, std::vector<std::pair<int,int>>>> NewListPossibleSignatures;
  for (auto & ePossSignature : WMVS.ListPossibleSignatures) {
    int NewDiag = g[ePossSignature.first];
    std::map<int,int> NewMap_mult;
    for (int i=0; i<int(ePossSignature.second.size()); i++) {
      std::pair<int,int> ePair = ePossSignature.second[i];
      int NewVal = g[ePair.first];
      int eMult = ePair.second;
      NewMap_mult[NewVal] = eMult;
    }
    std::vector<std::pair<int,int>> NewList_pair;
    for (auto & kv : NewMap_mult)
      NewList_pair.push_back({kv.first, kv.second});
    std::pair<int, std::vector<std::pair<int,int>>> e_pair{NewDiag, NewList_pair};
    NewListPossibleSignatures.push_back(e_pair);
  }
  WMVS.ListPossibleSignatures=NewListPossibleSignatures;
}



template<typename T, typename Tidx, typename F1, typename F2>
DataTraces GetDataTraces(F1 f1, F2 f2, WeightMatrixVertexSignatures<T> const& WMVS)
{
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
  //
  std::vector<int> V = Pairs_GetListPair(hS, nbMult);
  size_t e_pow = V.size() / 2;
  //
  // Now building the list of possibilities for the 2 added vertices
  //
  std::vector<std::vector<std::pair<int,int>>> new_list_signature;
  for (auto & e_pair : WMVS.ListPossibleSignatures) {
    std::vector<std::pair<int,int>> e_vect = e_pair.second;
    auto fInsert=[&](int e_val) -> void {
      for (auto & e_pair : e_vect) {
        if (e_pair.first == e_val) {
          e_pair.second++;
          return;
        }
      }
      e_vect.push_back({e_val,1});
    };
    fInsert(e_pair.first);
    fInsert(nbWeight + 1);
    new_list_signature.push_back(e_vect);
  }
  // nbRow : Adding the column corresponding to the diagonal values
  std::map<int,int> map_vert_nRow;
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    int isign = WMVS.ListSignatureByVertex[iRow];
    int iweight_specific = WMVS.ListPossibleSignatures[isign].first;
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
  for (size_t iCase=0; iCase<nbCase; iCase++) {
    for (size_t iH=0; iH<hS; iH++) {
      nbAdjacent += ListNbCase[iCase] * MatrixAdj(iH, iCase);
    }
  }
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
  auto f_adj=[&](size_t iVert, size_t jVert) -> void {
    DT.sg1.e[ListShift[iVert]] = jVert;
    ListShift[iVert]++;
  };
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
  std::cerr << "|GetCanonicalizationVector_KnownSignature|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time2).count() << "\n";
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
  std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> ePair = TRACES_GetCanonicalOrdering_ListGenerators_Arr<Tidx>(DT);
  std::vector<std::vector<Tidx>> LGen;
  for (auto& eListI : ePair.second) {
    std::vector<Tidx> eListO(nbRow);
    for (Tidx i=0; i<Tidx(nbRow); i++)
      eListO[i] = eListI[i];
    LGen.push_back(eListO);
  }
  //
  std::vector<Tidx> MapVectRev2 = GetCanonicalizationVector_KernelBis<Tidx>(nbRow, ePair.first);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetGroupCanonicalization_KnownSignature|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time2).count() << "\n";
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
  std::cerr << "|GetStabilizerWeightMatrix_KnownSignature|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
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
  WeightMatrixVertexSignatures<T> WMVS = ComputeVertexSignatures<T>(nbRow, f1, f2);
  size_t nbCase = WMVS.ListPossibleSignatures.size();
  std::vector<int> ListNbCase(nbCase, 0);
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    int iCase = WMVS.ListSignatureByVertex[iRow];
    ListNbCase[iCase]++;
  }
  std::vector<int> ListIdx;
  for (size_t iCase=0; iCase<nbCase; iCase++)
    ListIdx.push_back(iCase);
  std::sort(ListIdx.begin(), ListIdx.end(),
            [&](int idx1, int idx2) -> bool {
              return ListNbCase[idx1] < ListNbCase[idx2];});

  for (size_t idx=1; idx<nbCase; idx++) {
    std::vector<int> StatusCase(nbCase,0);
    for (size_t u=0; u<=idx; u++)
      StatusCase[ListIdx[u]] = 1;
    std::vector<Tidx> CurrentListIdx;
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      int iCase = WMVS.ListSignatureByVertex[iRow];
      if (StatusCase[iCase] == 1)
        CurrentListIdx.push_back(iRow);
    }
    size_t nbRow_res = CurrentListIdx.size();
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
  std::cerr << "We should never reach that stage\n";
  throw TerminalException{1};
}






/*
  ---F1/F2 : The first/second template function for creating the Weight matrix
  ---F3    : The function for testing acceptability of sets for consideration (e.g. rank function)
  ---F4    : The function for testing acceptability of small generators and returning the big
             generators.
*/
template<typename T, typename Tidx, typename F1, typename F2, typename F3, typename F4, typename F5>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>>  GetGroupCanonicalizationVector_Heuristic(size_t nbRow, F1 f1, F2 f2, F3 f3, F4 f4, F5 f5)
{
  WeightMatrixVertexSignatures<T> WMVS = ComputeVertexSignatures<T>(nbRow, f1, f2);
  RenormalizeWMVS(WMVS);
  size_t nbCase = WMVS.ListPossibleSignatures.size();
  std::vector<int> ListNbCase(nbCase, 0);
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    int iCase = WMVS.ListSignatureByVertex[iRow];
    ListNbCase[iCase]++;
  }
  std::vector<int> ListIdx;
  for (size_t iCase=0; iCase<nbCase; iCase++)
    ListIdx.push_back(iCase);
  std::sort(ListIdx.begin(), ListIdx.end(),
            [&](int idx1, int idx2) -> bool {
              return ListNbCase[idx1] < ListNbCase[idx2];});

  for (size_t idx=1; idx<nbCase; idx++) {
    std::vector<int> StatusCase(nbCase,0);
    for (size_t u=0; u<=idx; u++)
      StatusCase[ListIdx[u]] = 1;
    std::vector<Tidx> CurrentListIdx;
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      int iCase = WMVS.ListSignatureByVertex[iRow];
      if (StatusCase[iCase] == 1)
        CurrentListIdx.push_back(iRow);
    }
    size_t nbRow_res = CurrentListIdx.size();
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
        return {f5(ePair.first), LGen};
      }
    }
  }
  std::cerr << "We should never reach that stage\n";
  throw TerminalException{1};
}









#endif

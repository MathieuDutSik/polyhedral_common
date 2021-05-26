#ifndef INCLUDE_WEIGHT_MATRIX_SPECIFIED_H
#define INCLUDE_WEIGHT_MATRIX_SPECIFIED_H




template<typename T>
struct WeightMatrixVertexSignatures {
  int nbRow;
  int nbWeight;
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
  int nbWeight = ListWeight.size();
  return {nbRow, nbWeight, std::move(ListWeight), std::move(ListPossibleSignatures), std::move(ListSignatureByVertex)};
}


template<typename T, typename F1, typename F2, typename Tidx>
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




template<typename T, typename F1, typename F2, typename Tidx>
std::vector<Tidx> GetCanonicalizationVector_KnownSignature(F1 f1, F2 f2, WeightMatrixVertexSignatures<T> const& WMVS)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  int nbRow = WMVS.nbRow;
  DataTraces DT = GetDataTraces(f1, f2, WMVS);
  std::vector<unsigned int> cl = TRACES_GetCanonicalOrdering_Arr(DT);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetCanonicalizationVector_KnownSignature|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time2).count() << "\n";
#endif
  return GetCanonicalizationVector_KernelBis<Tidx>(nbRow, cl);
}



template<typename T, typename F1, typename F2, typename Tidx>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> GetGroupCanonicalization_KnownSignature(F1 f1, F2 f2, WeightMatrixVertexSignatures<T> const& WMVS)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  int nbRow = WMVS.nbRow;
  DataTraces DT = GetDataTraces(f1, f2, WMVS);
  std::pair<std::vector<unsigned int>, std::vector<std::vector<unsigned int>>> ePair = TRACES_GetCanonicalOrdering_ListGenerators_Arr(DT);
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



template<typename T, typename F1, typename F2, typename Tgroup>
Tgroup GetStabilizerWeightMatrix_KnownSignature(F1 f1, F2 f2, WeightMatrixVertexSignatures<T> const& WMVS)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int nbRow = WMVS.nbRow;
  DataTraces DT = GetDataTraces(f1, f2, WMVS);
  std::vector<std::vector<unsigned int>> ListGen = TRACES_GetListGenerators_Arr(DT);
  std::vector<std::vector<Tidx>> LGen;
  for (auto& eListI : ListGen) {
    std::vector<Tidx> eListO(nbRow);
    for (Tidx i=0; i<Tidx(nbRow); i++)
      eListO[i] = eListI[i];
    LGen.push_back(eListO);
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetStabilizerWeightMatrix_KnownSignature|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time2).count() << "\n";
#endif
  return Tgroup(LGen, nbRow);
}


/*
  ---F1/F2 : The first/second template function for creating the Weight matrix
  ---F3    : The function for testing acceptability of sets for consideration (e.g. rank function)
  ---F4    : The function for testing acceptability of small generators
*/
template<typename T, typename F1, typename F2, typename F3, typename F4, typename Tgroup>
Tgroup GetStabilizerWeightMatrix_Heuristic(int nbRow, F1 f1, F2 f2, F3 f3, F4 f4)
{
  using Telt = typename Tgroup::Telt;
  WeightMatrixVertexSignatures<T> WMVS = ComputeVertexSignatures(nbRow, f1, f2);
  size_t nbCase = WMVS.ListPossibleSignatures.size();
  std::vector<int> ListNbCase(nbCase, 0);
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    int iCase = WMVS.ListSignatureByVertex[iRow];
    ListNbCase[iCase]++;
  }
  std::vector<int> ListIdx;
  for (size_t iCase=0; iCase<nbCase; iCase++)
    ListIdx.push_back(iCase);
  std::sort(ListIdx.begin(), ListIdx.end(), [&](int idx1, int idx2) -> bool {
                                              return ListNbCase[idx1] < ListNbCase[idx2];});

  for (int idx=1; idx<nbCase; idx++) {
    std::vector<int> StatusCase(nbCase,0);
    for (int u=0; u<=idx; u++)
      StatusCase[ListIdx[u]] = 1;
    std::vector<int> CurrentListIdx;
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
      auto f2_res = [&](size_t jRow) -> void {
        f2(CurrentListIdx[jRow]);
      };
      WeightMatrixVertexSignatures<T> WMVS_res = ComputeVertexSignatures(nbRow_res, f1_res, f2_res);
      Tgroup GRP = GetStabilizerWeightMatrix_KnownSignature(f1_res, f2_res, WMVS_res);
      bool IsCorrect = true;
      std::vector<Telt> LGen;
      for (auto & eGen : GRP.GeneratorsOfGroup()) {
        if (IsCorrect) {
          EquivTest<Telt> test = f4(eGen);
          if (test.TheReply) {
            LGen.push_back(test.TheEquiv);
          } else {
            IsCorrect = false;
          }
        }
      }
      if (IsCorrect) {
        return Group(LGen, nbRow);
      }
    }
  }
}









#endif

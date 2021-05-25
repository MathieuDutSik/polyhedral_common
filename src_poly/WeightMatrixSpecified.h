#ifndef INCLUDE_WEIGHT_MATRIX_SPECIFIED_H
#define INCLUDE_WEIGHT_MATRIX_SPECIFIED_H




template<typename T>
struct WeightMatrixVertexSignatures {
  int nbRow;
  int nbWeight;
  std::vector<T> ListWeight;
  using Tvertex_signature = std::vector<std::pair<int,int>>;
  std::vector<Tvertex_signature> ListPossibleSignatures;
  std::vector<int> ListSignatureByVertex;
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

  using Tvertex_signature = std::vector<std::pair<int,int>>;
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
    for (size_t jRow=0; jRow<nbRow; jRow++) {
      T val = f2(jRow);
      int idx = get_T_idx(val);
      map_mult[idx]++;
    }
    Tvertex_signature the_sign;
    for (auto & kv : map_mult)
      the_sign.push_back({kv.first, kv.second});
    int idx_sign = get_Tvs_idx(the_sign);
    ListSignatureByVertex[iRow] = idx_sign;
  }
  int nbWeight = ListWeight.size();
  return {nbRow, nbWeight, std::move(ListWeight), std::move(ListPossibleSignatures), std::move(ListSignatureByVertex)};
}


template<typename T, typename F1, typename F2, typename Tidx>
DataTraces GetDataTraces(F1 f1, F2 f2, WeightMatrixVertexSignatures<T> const& DT)
{
  int nbRow = DT.nbRow;
  
}




template<typename T, typename F1, typename F2, typename Tidx>
std::vector<Tidx> GetCanonicalizationVector_KnownSignature(F1 f1, F2 f2, WeightMatrixVertexSignatures<T> const& DT)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  int nbRow = DT.nbRow;
  DataTraces data = GetDataTraces(f1, f2, DataSign);
  std::vector<unsigned int> cl = TRACES_GetCanonicalOrdering_Arr(DT);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetCanonicalizationVector_KnownSignature|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time2).count() << "\n";
#endif
  return GetCanonicalizationVector_KernelBis<Tidx>(nbRow, cl);
}



template<typename T, typename F1, typename F2, typename Tidx>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> GetGroupCanonicalization_KnownSignature(F1 f1, F2 f2, WeightMatrixVertexSignatures<T> const& DT)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  int nbRow = DT.nbRow;
  DataTraces data = GetDataTraces(f1, f2, DataSign);
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
Tgroup GetStabilizerWeightMatrix_KnownSignature(F1 f1, F2 f2, WeightMatrixVertexSignatures<T> const& DT)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int nbRow = DT.nbRow;
  DataTraces data = GetDataTraces(f1, f2, DataSign);
  std::pair<std::vector<unsigned int>, std::vector<std::vector<unsigned int>>> ePair = TRACES_GetCanonicalOrdering_ListGenerators_Arr(DT);
  std::vector<std::vector<Tidx>> LGen;
  for (auto& eListI : ePair.second) {
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




#endif

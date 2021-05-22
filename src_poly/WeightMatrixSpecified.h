#ifndef INCLUDE_WEIGHT_MATRIX_SPECIFIED_H
#define INCLUDE_WEIGHT_MATRIX_SPECIFIED_H




template<typename T>
struct WeightMatrixVertexSignatures {
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
  return {std::move(ListWeight), std::move(ListPossibleSignatures), std::move(ListSignatureByVertex)};
}









#endif

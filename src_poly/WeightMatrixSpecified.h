#ifndef INCLUDE_WEIGHT_MATRIX_SPECIFIED_H
#define INCLUDE_WEIGHT_MATRIX_SPECIFIED_H




template<typename T>
struct SummaryWeightMatrix {
  std::unordered_map<T,int> map_value;
  using Tvertex_signature = std::vector<std::pair<int,int>>;
  std::vector<Tvertex_signature> ListPossibleSignatures;
  std::vector<int> ListSignatureByVertex;
};





#endif

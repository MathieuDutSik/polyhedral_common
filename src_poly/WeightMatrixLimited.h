#ifndef INCLUDE_WEIGHT_MATRIX_LIMITED_H
#define INCLUDE_WEIGHT_MATRIX_LIMITED_H






template<bool is_symmetric_impl, typename T_impl>
struct WeightMatrixLimited {
public:
  static const bool is_symmetric = is_symmetric_impl;
  using T = T_impl;
private:
  void reorder_entries(std::vector<T> & list_weight, std::vector<size_t> & list_idx)
  {
    std::map<T, size_t> map_value;
    size_t n_ent = list_weight.size();
    for (size_t i=0; i<n_ent; i++)
      map_value[list_diag_weight_tmp[i]] = i;
    std::vector<size_t> g(n_ent);
    size_t idx=0;
    for (auto & kv : map_value) {
      g[kv.second] = idx;
      idx++;
    }
    size_t len = list_idx.size();
    for (size_t i=0; i<len; i++)
      list_idx[i] = g[list_idx[i]];
    std::vector<T> list_weight_new(n_ent);
    for (size_t i=0; i<n_ent; i++)
      list_weight_new[g[i]] = list_weight[i];
    list_weight = std::move(list_weight_new);
  }
  void compute_shift_size(std::vector<size_t> & list_sizes, std::vector<size_t> & list_shift, std::vector<size_t> & list_element,
                          const std::vector<size_t> & list_idx, const std::vector<T> & list_weight)
  {
    size_t n_weight = list_weight.size();
    size_t len = list_idx.size();
    list_sizes.resize(n_weight);
    list_shift.resize(n_weight+1);
    for (size_t i=0; i<n_weight; i++)
      list_sizes[i] = 0;
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      size_t iWeight = list_diag_idx[iRow];
      list_diag_sizes[iWeight]++;
    }
    list_shift[0] = 0;
    for (size_t i=0; i<n_weight; i++)
      list_shift[i+1] = list_shift[i] + list_sizes[i];
    std::vector<size_t> list_pos = list_shift;
    list_element.resize(len);
    for (size_t i=0; i<len; i++) {
      size_t iWeight = list_idx[i];
      size_t pos = list_pos[iWeight];
      list_element[pos] = iRow;
      list_pos[iWeight] = pos + 1;
    }
  }
public:
  template<typename F1, typename F2>
  WeightMatrix(size_t const& _nbRow, F1 f1, F2 f2, size_t max_offdiag) : nbRow(_nbRow)
  {
    // Computing the diagional values
    list_diag_idx.resize(nbRow);
    std::unordered_map<T, size_t> unordmap_value;
    size_t idxWeight=0;
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      f1(iRow);
      T val = f2(iRow);
      size_t & idx = unordmap_value[val];
      if (idx == 0) {
        idxWeight++;
        idx = idxWeight;
        list_diag_weight.push_back(val);
      }
      size_t pos_val = idx - 1;
      list_diag_idx[iRow] = pos_val;
    }
    size_t n_weight = unordmap_value.size();
    // Coding the reordering of the list_diag_weight
    reorder_entries(list_diag_weight, list_diag_idx);
    // Computing the shift and sizes
    compute_shift_size(list_diag_sizes, list_diag_shift, list_diag_element, list_diag_idx, list_diag_weight);
    // Computing the set of off diagonal entries that we are interested in.
    using Tpair = std::pair<size_t, size_t>;
    std::map<size_t,std::vector<Tpair>> map_pair;
    for (size_t i_weight=0; i_weight<n_weight; i_weight++) {
      size_t start=0;
      if (constexpr is_symmetric)
        start = i_weight;
      size_t siz1 = list_diag_sizes[i_weight];
      for (size_t j_weight=0; j_weight<n_weight; j_weight++) {
        if (i_weight == j_weight) {
          if (constexpr is_symmetric) {
            size_t n_poss = siz1 * (siz1 - 1) / 2;
            map_pair[n_poss].push_back({i_weight, i_weight});
          } else {
            size_t n_poss = siz1 * (siz1 - 1);
            map_pair[n_poss].push_back({i_weight, i_weight});
          }
        } else {
          size_t siz2 = list_diag_sizes[j_weight];
          size_t n_poss = siz1 * siz2;
          map_pair[n_poss].push_back({i_weight, j_weight});
        }
      }
    }
    std::vector<Tpair> list_selected;
    auto iter = map_pair.cbegin();
    size_t sel_size = 0;
    while(true) {
      if (iter == map_pair.cend())
        break;
      const std::vector<Tpair> & e_pair = iter->second;
      size_t tot_siz = iter->first * (e_pair.size());
      if (sel_size + tot_siz > max_offdiag)
        break;
      sel_size += tot_siz;
      list_selected.insert(list_selected.end(), e_pair.begin(), e_pair.end());
      iter++;
    }
    // Computing the arrays themselves
    list_offdiag_idx.resize(sel_size);
    unordmap_value.clear();
    idxWeight=0;
    size_t pos = 0;
    auto insertval=[&](size_t const& jRow) -> void {
      T val = f2(jRow);
      size_t & idx = unordmap_value[val];
      if (idx == 0) {
        idxWeight++;
        idx = idxWeight;
        list_offdiag_weight.push_back(val);
      }
      size_t pos_val = idx - 1;
      list_offdiag_idx[pos] = pos_val;
      pos++;
    };
    for (auto& e_ent : list_selected) {
      size_t i_wei = e_ent.first;
      size_t j_wei = e_ent.second;
      if (i_wei == j_wei) {
        size_t siz = list_sizes[i_wei];
        size_t shift = list_shift[i_wei];
        if (constexpr is_symmetric) {
          for (size_t i=0; i<siz; i++) {
            size_t ipos = shift + i;
            size_t iRow = list_diag_element[ipos];
            f1(iRow);
            for (size_t j=i+1; j<siz; j++) {
              size_t jpos = shift + j;
              size_t jRow = list_diag_element[jpos];
              insertval(jRow);
            }
          }
        } else {
          for (size_t i=0; i<siz; i++) {
            size_t ipos = shift + i;
            size_t iRow = list_diag_element[ipos];
            f1(iRow);
            for (size_t j=0; j<siz; j++) {
              if (i != j) {
                size_t jpos = shift + j;
                size_t jRow = list_diag_element[jpos];
                insertval(jRow);
              }
            }
          }
        }
      } else {
        size_t i_siz = list_diag_sizes[i_wei];
        size_t i_shift = list_diag_shift[i_wei];
        size_t j_siz = list_diag_sizes[j_wei];
        size_t j_shift = list_diag_shift[j_wei];
        for (size_t i=0; i<i_siz; i++) {
          size_t ipos = i_shift + i;
          size_t iRow = list_diag_element[ipos];
          f1(iRow);
          for (size_t j=0; j<j_siz; j++) {
            size_t jpos = j_shift + j;
            size_t jRow = list_diag_element[jpos];
            insertval(jRow);
          }
        }
      }
    }
    // Reordering of the unordmap_value.
    reorder_entries(list_offdiag_weight, list_offdiag_idx);
    // Now computing shifting information
    compute_shift_size(list_offdiag_sizes, list_offdiag_shift, list_offdiag_element, list_offdiag_idx, list_offdiag_weight);
  }
  template<typename Tgroup>
  WeightMatrix(const Tgroup& GRP) : nbRow(GRP.n_act())
  {
    
  }
private:
  size_t nbRow;
  //
  std::vector<T> list_diag_weight;
  std::vector<size_t> list_diag_idx;
  std::vector<size_t> list_diag_sizes;
  std::vector<size_t> list_diag_shift;
  std::vector<size_t> list_diag_elements;
  //
  std::vector<T> list_offdiag_weight;
  std::vector<size_t> list_offdiag_idx;
  std::vector<size_t> list_offdiag_sizes;
  std::vector<size_t> list_offdiag_shifts;
  std::vector<size_t> list_offdiag_elements;
};



















#endif

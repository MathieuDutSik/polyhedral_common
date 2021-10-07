#ifndef INCLUDE_WEIGHT_MATRIX_LIMITED_H
#define INCLUDE_WEIGHT_MATRIX_LIMITED_H






template<bool is_symmetric_impl, typename T_impl, typename Tidx_value_impl>
struct WeightMatrixLimited {
public:
  static const bool is_symmetric = is_symmetric_impl;
  using T = T_impl;
  using Tidx_value = Tidx_value_impl;
  template<typename F1, typename F2>
  WeightMatrix(size_t const& _nbRow, F1 f1, F2 f2, size_t max_offdiag) : nbRow(_nbRow)
  {
    // Computing the diagional values
    list_idx_diag.resize(nbRow);
    std::unordered_map<T, Tidx_value> ValueMap;
    Tidx_value idxWeight=0;
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      f1(iRow);
      T val = f2(iRow);
      Tidx_value & idx = ValueMap[val];
      if (idx == 0) {
        idxWeight++;
        idx = idxWeight;
        list_weight_diag.push_back(val);
      }
      Tidx_value pos_val = idx - 1;
      list_idx_diag[iRow] = pos_val;
    }
    size_t n_weight = ValueMap.size();
    if (n_weight > size_t(std::numeric_limits<Tidx_value>::max() - 1)) {
      std::cerr << "We have |ValueMap|=" << n_weight << "\n";
      std::cerr << "std::numeric_limits<Tidx_value>::max()=" << size_t(std::numeric_limits<Tidx_value>::max()) << "\n";
      std::cerr << "We also need some space for the missing values\n";
      throw TerminalException{1};
    }
    // Coding the reordering of the list_weight_diag

    // Computing the shift and sizes
    list_sizes.resize(n_weight);
    list_shift.resize(n_weight+1);
    for (size_t i=0; i<n_weight; i++)
      list_sizes[i] = 0;
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      size_t iWeight = list_idx_diag[iRow];
      list_sizes[iWeight]++;
    }
    list_shift[0] = 0;
    for (size_t i=0; i<n_weight; i++)
      list_shift[i+1] = list_shift[i] + list_sizes[i];
    std::vector<size_t> list_pos = list_shift;
    list_element.resize(nbRow);
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      size_t iWeight = list_idx_diag[iRow];
      size_t pos = list_shift[iWeight];
      list_element[pos] = iRow;
      list_shift[iWeight] = pos + 1;
    }
    // Computing the set of off diagonal entries that we are interested in.
    using Tsize = std::pair<size_t, std::pair<size_t, size_t>>;
    auto f_ker_comp=[](const size_t& x, const size_t& y) -> std::pair<bool,bool> {
      if (x < y)
        return {true, true};
      if (x > y)
        return {true, false};
      return {false, true};
    };
    auto f_comp = [](const Tsize& x, const Tsize& y) -> bool {
      auto test1 = f_ker_comp(x.first, y.first);
      if (test1.first) return test1.second;
      //
      auto test2 = f_ker_comp(x.second.first, y.second.first);
      if (test2.first) return test2.second;
      //
      return x.second.second < y.second.second;
    };
    std::unordered_set<T,decltype(f_comp)> set_pair(f_comp);
    for (size_t i_weight=0; i_weight<n_weight; i_weight++) {
      size_t start=0;
      if (constexpr is_symmetric)
        start = i_weight;
      size_t siz1 = list_sizes[i_weight];
      for (size_t j_weight=0; j_weight<n_weight; j_weight++) {
        if (i_weight == j_weight) {
          if (constexpr is_symmetric) {
            size_t n_poss = siz1 * (siz1 - 1) / 2;
            set_pair.push_back({n_poss, {i_weight, i_weight}});
          } else {
            size_t n_poss = siz1 * (siz1 - 1);
            set_pair.push_back({n_poss, {i_weight, i_weight}});
          }
        } else {
          size_t siz2 = list_sizes[j_weight];
          size_t n_poss = siz1 * siz2;
          set_pair.push_back({n_poss, {i_weight, j_weight}});
        }
      }
    }
    std::vector<T> list_selected;
    auto iter = set_pair.cbegin();
    size_t sel_size = 0;
    while(true) {
      if (iter == set_pair.cend())
        break;
      const Tsize& ent = *iter;
      if (sel_size + ent.first > max_offdiag)
        break;
      list_selected.push_back(ent);
      iter++;
    }
    // Computing the arrays themselves
    for (auto& e_ent : list_selected) {
      
    }
  }

private:
  size_t nbRow;
  std::vector<T> list_weight_diag;
  //
  std::vector<Tidx_value> list_idx_diag;
  std::vector<size_t> list_sizes, list_shift;
  std::vector<size_t> list_elements;
  //
  std::vector<T> list_weight_offdiag;
  std::vector<size_t> list_offdiag_sizes;
  std::vector<size_t> list_offdiag_shifts;
  std::vectior<size_t
};



















#endif

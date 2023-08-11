// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_WEIGHTMATRIXLIMITED_H_
#define SRC_POLY_WEIGHTMATRIXLIMITED_H_

// clang-format off
#include "GRP_GroupFct.h"
#include <limits>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

template <bool is_symmetric_impl, typename T_impl> struct WeightMatrixLimited {
public:
  static const bool is_symmetric = is_symmetric_impl;
  using T = T_impl;
  using Tpair = std::pair<size_t, size_t>;

private:
  void reorder_entries(std::vector<T> &list_weight,
                       std::vector<size_t> &list_idx) const {
    std::map<T, size_t> map_value;
    size_t n_ent = list_weight.size();
    for (size_t i = 0; i < n_ent; i++)
      map_value[list_weight[i]] = i;
    std::vector<size_t> g(n_ent);
    size_t idx = 0;
    for (auto &kv : map_value) {
      g[kv.second] = idx;
      idx++;
    }
    size_t len = list_idx.size();
    for (size_t i = 0; i < len; i++)
      list_idx[i] = g[list_idx[i]];
    std::vector<T> list_weight_new(n_ent);
    for (size_t i = 0; i < n_ent; i++)
      list_weight_new[g[i]] = list_weight[i];
    list_weight = std::move(list_weight_new);
  }
  void compute_shift_size(std::vector<size_t> &list_sizes,
                          std::vector<size_t> &list_shift,
                          std::vector<size_t> &list_element,
                          const std::vector<size_t> &list_idx,
                          const size_t &n_weight) const {
    size_t len = list_idx.size();
    std::cerr << "compute_shift_size n_weight=" << n_weight << " len=" << len
              << "\n";
    list_sizes.resize(n_weight);
    list_shift.resize(n_weight + 1);
    for (size_t i = 0; i < n_weight; i++)
      list_sizes[i] = 0;
    for (size_t i = 0; i < len; i++) {
      size_t iWeight = list_idx[i];
      list_sizes[iWeight]++;
    }
    list_shift[0] = 0;
    for (size_t i = 0; i < n_weight; i++)
      list_shift[i + 1] = list_shift[i] + list_sizes[i];
    std::vector<size_t> list_pos = list_shift;
    list_element.resize(len);
    for (size_t i = 0; i < len; i++) {
      size_t iWeight = list_idx[i];
      size_t pos = list_pos[iWeight];
      list_element[pos] = i;
      list_pos[iWeight] = pos + 1;
    }
  }
  void set_mat_select_pair(const size_t &n_weight,
                           const std::vector<Tpair> &list_selected) {
    size_t eProd = n_weight * n_weight;
    mat_select_pair.resize(eProd);
    for (size_t i = 0; i < eProd; i++)
      mat_select_pair[i] = std::numeric_limits<size_t>::max();
    size_t len = list_selected.size();
    for (size_t u = 0; u < len; u++) {
      size_t i = list_selected[u].first;
      size_t j = list_selected[u].second;
      std::cerr << "u=" << u << " i=" << i << " j=" << j
                << " n_weight=" << n_weight << "\n";
      if constexpr (is_symmetric) {
        mat_select_pair[i + n_weight * j] = u;
        mat_select_pair[j + n_weight * i] = u;
      } else {
        mat_select_pair[i + n_weight * j] = u;
      }
    }
    std::cerr << "mat_select_pair=\n";
    for (size_t i = 0; i < n_weight; i++) {
      for (size_t j = 0; j < n_weight; j++) {
        std::cerr << " " << mat_select_pair[i + n_weight * j];
      }
      std::cerr << "\n";
    }
  }
  std::pair<size_t, std::vector<Tpair>> get_list_selected(size_t max_offdiag) {
    std::map<size_t, std::vector<Tpair>> map_pair;
    size_t n_weight = list_diag_weight.size();
    for (size_t i_weight = 0; i_weight < n_weight; i_weight++) {
      size_t start = 0;
      if constexpr (is_symmetric)
        start = i_weight;
      size_t siz1 = list_diag_sizes[i_weight];
      for (size_t j_weight = start; j_weight < n_weight; j_weight++) {
        if (i_weight == j_weight) {
          if constexpr (is_symmetric) {
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
    while (true) {
      if (iter == map_pair.cend())
        break;
      const std::vector<Tpair> &e_pair = iter->second;
      size_t tot_siz = iter->first * (e_pair.size());
      if (sel_size + tot_siz > max_offdiag)
        break;
      sel_size += tot_siz;
      list_selected.insert(list_selected.end(), e_pair.begin(), e_pair.end());
      iter++;
    }
    set_mat_select_pair(n_weight, list_selected);
    return {sel_size, std::move(list_selected)};
  }
  void compute_list_revdiag_elements() {
    list_revdiag_elements.resize(nbRow);
    size_t n_weight = list_diag_weight.size();
    std::cerr << "|list_diag_shift|=" << list_diag_shift.size()
              << " |list_diag_sizes|=" << list_diag_sizes.size() << "\n";
    for (size_t i_w = 0; i_w < n_weight; i_w++) {
      std::cerr << "i_w=" << i_w << "\n";
      size_t shift = list_diag_shift[i_w];
      size_t siz = list_diag_sizes[i_w];
      for (size_t u = 0; u < siz; u++) {
        size_t pt = list_diag_element[shift + u];
        list_revdiag_elements[pt] = u;
      }
    }
  }
  void print_variables() const {
    std::cerr << "list_diag_sizes =";
    for (auto &val : list_diag_sizes)
      std::cerr << " " << val;
    std::cerr << "\n";
    //
    std::cerr << "list_offdiag_sizes =";
    for (auto &val : list_offdiag_sizes)
      std::cerr << " " << val;
    std::cerr << "\n";
    //
    std::cerr << "list_offdiag_shifts =";
    for (auto &val : list_offdiag_shifts)
      std::cerr << " " << val;
    std::cerr << "\n";
    //
    std::cerr << "|list_offdiag_idx|=" << list_offdiag_idx.size() << "\n";
  }

public:
  template <typename F1, typename F2>
  WeightMatrixLimited(size_t const &_nbRow, F1 f1, F2 f2, size_t max_offdiag)
      : nbRow(_nbRow) {
    // Computing the diagional values
    list_diag_idx.resize(nbRow);
    std::unordered_map<T, size_t> unordmap_value;
    size_t idxWeight = 0;
    for (size_t iRow = 0; iRow < nbRow; iRow++) {
      f1(iRow);
      T val = f2(iRow);
      size_t &idx = unordmap_value[val];
      if (idx == 0) {
        idxWeight++;
        idx = idxWeight;
        list_diag_weight.push_back(val);
      }
      size_t pos_val = idx - 1;
      list_diag_idx[iRow] = pos_val;
    }
    // Coding the reordering of the list_diag_weight
    reorder_entries(list_diag_weight, list_diag_idx);
    // Computing the shift and sizes
    compute_shift_size(list_diag_sizes, list_diag_shift, list_diag_element,
                       list_diag_idx, list_diag_weight.size());
    // Computing the shift and sizes
    compute_list_revdiag_elements();
    // Computing the set of off diagonal entries that we are interested in.
    std::pair<size_t, std::vector<Tpair>> pair = get_list_selected(max_offdiag);
    const std::vector<Tpair> &list_selected = pair.second;
    // Computing the arrays themselves
    size_t n_offdiag_weight = 0;
    size_t shift_weight = 0;
    for (auto &e_ent : list_selected) {
      unordmap_value.clear();
      size_t jdxWeight = 0;
      std::vector<size_t> list_part_idx;
      std::vector<T> list_part_weight;
      auto insertval = [&](size_t const &jRow) -> void {
        T val = f2(jRow);
        size_t &idx = unordmap_value[val];
        if (idx == 0) {
          jdxWeight++;
          idx = jdxWeight;
          list_part_weight.push_back(val);
        }
        size_t pos_val = idx - 1;
        list_part_idx.push_back(pos_val);
      };
      size_t i_wei = e_ent.first;
      size_t j_wei = e_ent.second;
      if (i_wei == j_wei) {
        size_t siz = list_diag_sizes[i_wei];
        size_t shift = list_diag_shift[i_wei];
        if constexpr (is_symmetric) {
          for (size_t i = 0; i < siz; i++) {
            size_t ipos = shift + i;
            size_t iRow = list_diag_element[ipos];
            f1(iRow);
            for (size_t j = i + 1; j < siz; j++) {
              size_t jpos = shift + j;
              size_t jRow = list_diag_element[jpos];
              insertval(jRow);
            }
          }
        } else {
          for (size_t i = 0; i < siz; i++) {
            size_t ipos = shift + i;
            size_t iRow = list_diag_element[ipos];
            f1(iRow);
            for (size_t j = 0; j < siz; j++) {
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
        for (size_t i = 0; i < i_siz; i++) {
          size_t ipos = i_shift + i;
          size_t iRow = list_diag_element[ipos];
          f1(iRow);
          for (size_t j = 0; j < j_siz; j++) {
            size_t jpos = j_shift + j;
            size_t jRow = list_diag_element[jpos];
            insertval(jRow);
          }
        }
      }
      // reordering the weight
      reorder_entries(list_part_weight, list_part_idx);
      for (auto &val : list_part_idx)
        list_offdiag_idx.push_back(shift_weight + val);
      shift_weight += list_part_weight.size();
      n_offdiag_weight += list_part_weight.size();
    }
    // Now computing shifting information
    compute_shift_size(list_offdiag_sizes, list_offdiag_shifts,
                       list_offdiag_elements, list_offdiag_idx,
                       n_offdiag_weight);
    print_variables();
  }
  template <typename Tgroup>
  WeightMatrixLimited(const Tgroup &GRP, const size_t &max_offdiag)
      : nbRow(GRP.n_act()) {
    using Telt = typename Tgroup::Telt;
    using Tidx = typename Tgroup::Telt::Tidx;
    std::vector<Telt> LGen = GRP.GeneratorsOfGroup();
    vectface vf = DecomposeOrbitPoint_KernelFull<Telt>(nbRow, LGen);
    size_t len = vf.size();
    list_diag_weight.resize(len);
    for (size_t i = 0; i < len; i++)
      list_diag_weight[i] = i;
    list_diag_idx.resize(nbRow);
    for (size_t i = 0; i < len; i++) {
      Face f = vf[i];
      boost::dynamic_bitset<>::size_type pos = f.find_first();
      while (pos != boost::dynamic_bitset<>::npos) {
        list_diag_idx[pos] = i;
        pos = f.find_next(pos);
      }
    }
    compute_shift_size(list_diag_sizes, list_diag_shift, list_diag_element,
                       list_diag_idx, list_diag_weight.size());
    compute_list_revdiag_elements();
    // Selecting the groups.
    std::pair<size_t, std::vector<Tpair>> pair = get_list_selected(max_offdiag);
    const size_t &sel_size = pair.first;
    const std::vector<Tpair> &list_selected = pair.second;
    //
    size_t idx_weight = 0;
    list_offdiag_idx.resize(sel_size);
    size_t pos_offdiag = 0;
    auto f_insert = [&](const size_t &i_wei, const size_t &j_wei) -> void {
      size_t siz = list_diag_sizes[i_wei];
      size_t shift = list_diag_shift[i_wei];
      size_t siz_B, shift_B;
      size_t len;
      std::vector<size_t> map_rev;
      std::vector<size_t> map;
      size_t pos = 0;
      if (i_wei == j_wei) {
        siz_B = siz;
        shift_B = shift;
        if constexpr (is_symmetric) {
          map_rev.resize(siz * siz);
          len = siz * (siz - 1) / 2;
          map.resize(2 * len);
          for (size_t i = 0; i < siz; i++) {
            for (size_t j = i + 1; j < siz; j++) {
              map_rev[i + siz * j] = pos;
              map_rev[j + siz * i] = pos;
              map[2 * pos] = i;
              map[2 * pos + 1] = j;
              pos++;
            }
          }
        } else {
          map_rev.resize(siz * siz);
          len = siz * (siz - 1);
          map.resize(2 * len);
          for (size_t i = 0; i < siz; i++) {
            for (size_t j = 0; j < siz; j++) {
              if (i != j) {
                map_rev[i + siz * j] = pos;
                map[2 * pos] = i;
                map[2 * pos + 1] = j;
                pos++;
              }
            }
          }
        }
      } else {
        siz_B = list_diag_sizes[j_wei];
        shift_B = list_diag_shift[j_wei];
        map_rev.resize(siz * siz_B);
        len = siz * siz_B;
        map.resize(2 * len);
        for (size_t i = 0; i < siz; i++) {
          for (size_t j = 0; j < siz_B; j++) {
            map_rev[i + siz * j] = pos;
            map[2 * pos] = i;
            map[2 * pos + 1] = j;
            pos++;
          }
        }
      }
      std::vector<Telt> LGenMap;
      for (auto &eGen : LGen) {
        std::vector<Tidx> eList(len);
        for (size_t u = 0; u < len; u++) {
          size_t i = map[2 * u];
          size_t i_B = map[2 * u + 1];
          size_t ePt1 = list_diag_element[shift + i];
          size_t fPt1 = list_diag_element[shift_B + i_B];
          size_t ePt2 = OnPoints(ePt1, eGen);
          size_t fPt2 = OnPoints(fPt1, eGen);
          size_t j = list_revdiag_elements[ePt2];
          size_t j_B = list_revdiag_elements[fPt2];
          size_t u_img = map_rev[j + siz * j_B];
          eList[u] = u_img;
        }
        Telt eGB(eList);
        LGenMap.emplace_back(std::move(eGB));
      }
      vectface vf_b = DecomposeOrbitPoint_KernelFull(len, LGenMap);
      size_t len_b = vf_b.size();
      for (size_t i_b = 0; i_b < len_b; i_b++) {
        Face f_b = vf_b[i_b];
        boost::dynamic_bitset<>::size_type pos = f_b.find_first();
        while (pos != boost::dynamic_bitset<>::npos) {
          list_offdiag_idx[pos_offdiag + pos] = idx_weight;
          pos = f_b.find_next(pos);
        }
        list_offdiag_weight.push_back(idx_weight);
        idx_weight++;
      }
      pos_offdiag += len;
    };
    for (auto &e_ent : list_selected) {
      size_t i_wei = e_ent.first;
      size_t j_wei = e_ent.second;
      f_insert(i_wei, j_wei);
    }
    compute_shift_size(list_offdiag_sizes, list_offdiag_shifts,
                       list_offdiag_elements, list_offdiag_idx,
                       list_offdiag_weight.size());
    print_variables();
  }
  size_t get_hash(const Face &f) const {
    std::pair<int, std::vector<size_t>> pair = get_smallest_set(f);
    const std::vector<size_t> &eList = pair.second;
    const int &set_val = pair.first;
    //
    size_t nbVert = eList.size();
    size_t nw_diag = list_diag_weight.size();
    size_t nw_offdiag = list_offdiag_weight.size();
    std::vector<size_t> eInv(nw_diag + 2 * nw_offdiag, 0);
    for (auto &eVal : eList) {
      size_t iWeight = list_diag_idx[eVal];
      eInv[iWeight]++;
    }
    auto insert_offdiag_pair = [&](size_t eVal, size_t fVal,
                                   const size_t &shift_index) -> void {
      size_t iWeight = list_diag_idx[eVal];
      size_t jWeight = list_diag_idx[fVal];
      size_t pos = mat_select_pair[iWeight + nw_diag * jWeight];
      if (pos != std::numeric_limits<size_t>::max()) {
        size_t i = list_revdiag_elements[eVal];
        size_t j = list_revdiag_elements[fVal];
        size_t siz1 = list_diag_sizes[iWeight];
        size_t siz2 = list_diag_sizes[jWeight];
        size_t pos_B;
        if (iWeight != jWeight) {
          pos_B = i + siz1 * j;
        } else {
          if constexpr (is_symmetric) {
            if (i < j) {
              pos_B = i * (i - 1) / 2 + j - 1;
            } else {
              pos_B = j * (j - 1) / 2 + i - 1;
            }
          } else {
            std::cerr << "Using the wrong formula\n";
            // false
            pos_B = i + siz2 * j;
          }
        }
        size_t pos_C = list_offdiag_shifts[pos] + pos_B;
        if (pos_C >= list_offdiag_idx.size()) {
          throw TerminalException{1};
        }
        size_t kWeight = list_offdiag_idx[pos_C];
        eInv[nw_diag + shift_index + kWeight]++;
      }
    };
    auto insert_offdiag_pair_gen = [&](const size_t &eVal, const size_t &fVal,
                                       const size_t &shift_index) -> void {
      if constexpr (is_symmetric) {
        size_t iWeight = list_diag_idx[eVal];
        size_t jWeight = list_diag_idx[fVal];
        if (iWeight < jWeight) {
          return insert_offdiag_pair(eVal, fVal, shift_index);
        } else {
          return insert_offdiag_pair(fVal, eVal, shift_index);
        }
      } else {
        return insert_offdiag_pair(eVal, fVal, shift_index);
      }
    };
    if (nbVert < 50) {
      // In that case, we look at the full structure
      for (size_t u = 0; u < nbVert; u++) {
        size_t eVal = eList[u];
        for (size_t iRow = 0; iRow < nbRow; iRow++) {
          if (iRow != eVal) {
            size_t shift_index = 0;
            if (f[iRow] == set_val)
              shift_index = nw_offdiag;
            insert_offdiag_pair_gen(eVal, iRow, shift_index);
          }
        }
      }
    } else {
      for (size_t u = 0; u < nbVert; u++) {
        size_t eVal = eList[u];
        for (size_t v = 0; v < nbVert; v++) {
          if (u != v) {
            size_t fVal = eList[u];
            insert_offdiag_pair_gen(eVal, fVal, 0);
          }
        }
      }
    }
    return std::hash<std::vector<size_t>>()(eInv);
  }

private:
  size_t nbRow;
  //
  std::vector<T> list_diag_weight;
  std::vector<size_t> list_diag_idx;
  std::vector<size_t> list_diag_sizes;
  std::vector<size_t> list_diag_shift;
  std::vector<size_t> list_diag_element;
  std::vector<size_t> list_revdiag_elements;
  //
  std::vector<size_t> mat_select_pair;
  std::vector<T> list_offdiag_weight;
  std::vector<size_t> list_offdiag_idx;
  std::vector<size_t> list_offdiag_sizes;
  std::vector<size_t> list_offdiag_shifts;
  std::vector<size_t> list_offdiag_elements;
};

// clang-format off
#endif  // SRC_POLY_WEIGHTMATRIXLIMITED_H_
// clang-format on

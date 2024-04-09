// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POINCARE_POLYHEDRON_FINITE_MATRIX_GROUP_H_
#define SRC_POINCARE_POLYHEDRON_FINITE_MATRIX_GROUP_H_

template <typename T>
std::vector<CombElt<T>>
GroupGeneration(std::vector<CombElt<T>> const &input_l_ent) {
  std::vector<CombElt<T>> l_ent = input_l_ent;
  while (true) {
    std::unordered_set<CombElt<T>> s_GenElt;
    std::vector<CombElt<T>> l_GenElt;
    auto f_insert = [&](CombElt<T> const &eElt) -> void {
      for (auto &fElt : l_GenElt)
        if (fElt == eElt)
          return;
      l_GenElt.push_back(eElt);
    };
    size_t n_ent = l_ent.size();
    for (size_t i_ent = 0; i_ent < n_ent; i_ent++) {
      CombElt<T> const &pe1 = l_ent[i_ent];
      for (size_t j_ent = 0; j_ent < n_ent; j_ent++) {
        CombElt<T> const &pe2 = l_ent[j_ent];
        CombElt<T> prod = ProductComb(pe1, pe2);
        s_GenElt.insert(prod);
        f_insert(prod);
      }
    }
    if (s_GenElt.size() != l_GenElt.size()) {
      std::cerr << "|s_GenElt|=" << s_GenElt.size()
                << " |l_GenElt|=" << l_GenElt.size() << "\n";
      std::cerr << "sizes are different\n";
      throw TerminalException{1};
    }
    if (s_GenElt.size() == n_ent) {
      return l_ent;
    }
    l_ent.clear();
    for (auto &e_ent : s_GenElt)
      l_ent.push_back(e_ent);
  }
}

// a right coset is of the form Ug
template <typename T>
std::vector<CombElt<T>>
IdentifyRightCosets(std::vector<CombElt<T>> const &l_ent,
                    std::vector<CombElt<T>> const &list_grp_elt) {
  std::unordered_set<CombElt<T>> s_coset;
  auto f_insert = [&](CombElt<T> const &pe) -> void {
    for (auto &e_grp_elt : list_grp_elt) {
      CombElt<T> prod = ProductComb(e_grp_elt, pe);
      if (s_coset.count(prod) == 1)
        break;
    }
    s_coset.insert(pe);
  };
  for (auto &pe : l_ent)
    f_insert(pe);
  std::vector<CombElt<T>> l_coset;
  for (auto &e_coset : s_coset)
    l_coset.push_back(e_coset);
  return l_coset;
}

// a left coset is of the form gU
template <typename T>
std::vector<CombElt<T>>
IdentifyLeftCosets(std::vector<CombElt<T>> const &l_ent,
                   std::vector<CombElt<T>> const &list_grp_elt) {
  std::unordered_set<CombElt<T>> s_coset;
  auto f_insert = [&](CombElt<T> const &pe) -> void {
    for (auto &e_grp_elt : list_grp_elt) {
      CombElt<T> prod = ProductComb(pe, e_grp_elt);
      if (s_coset.find(prod) != s_coset.end()) {
        std::cerr << "find matching\n";
        return;
      }
    }
    s_coset.insert(pe);
  };
  for (auto &pe : l_ent)
    f_insert(pe);
  std::vector<CombElt<T>> l_coset(s_coset.begin(), s_coset.end());
  std::cerr << "|l_ent|=" << l_ent.size()
            << " |list_grp_elt|=" << list_grp_elt.size()
            << " |l_coset|=" << l_coset.size() << "\n";
  return l_coset;
}

template <typename T>
std::vector<CombElt<T>>
IdentifyDoubleCosets(MyVector<T> const &x, std::vector<CombElt<T>> const &l_ent,
                     std::vector<CombElt<T>> const &list_grp_elt) {
  std::unordered_map<MyVector<T>, CombElt<T>> map;
  auto f_insert = [&](CombElt<T> const &pe) -> void {
    MyVector<T> x2 = pe.mat.transpose() * x;
    for (auto &e_grp_elt : list_grp_elt) {
      MyVector<T> x3 = e_grp_elt.mat.transpose() * x2;
      auto iter = map.find(x3);
      if (iter != map.end()) {
        std::cerr << "find matching\n";
        return;
      }
    }
    map[x2] = pe;
  };
  for (auto &pe : l_ent)
    f_insert(pe);
  std::vector<CombElt<T>> l_coset;
  for (auto &kv : map)
    l_coset.push_back(kv.second);
  std::cerr << "|l_ent|=" << l_ent.size()
            << " |list_grp_elt|=" << list_grp_elt.size()
            << " |l_coset|=" << l_coset.size() << "\n";
  return l_coset;
}

// clang-format off
#endif  // SRC_POINCARE_POLYHEDRON_FINITE_MATRIX_GROUP_H_
// clang-format on

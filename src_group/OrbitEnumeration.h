// Copyright (C) 2024 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_GRP_ORBITENUMERATION_H_
#define SRC_GROUP_GRP_ORBITENUMERATION_H_

// clang-format off
#include "GRP_GroupFct.h"
// clang-format on

/*
  We are iterating by finding the minimal orbit.
  This relies on the group being relatively small.
  This is a tree search. 
 */
template<typename Tgroup, typename Fextensible>
void SubsetOrbitEnumeration(Tgroup const& GRP, Fextensible f_extensible) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx miss_val = std::numeric_limits<Tidx>::max();
  Tidx nbVert = GRP.n_act();
  std::vector<Telt> l_elt;
  for (auto & elt : GRP) {
    l_elt.push_back(elt);
  }
  Face f1(nbVert), f2(nbVert), f3(nbVert);
  auto compare_f1_f2=[&]() -> bool {
    // returns true if f1 <= f2
    for (Tidx u=0; u<nbVert; u++) {
      if (f1[u] == 1 && f2[u] == 0) {
        return true;
      }
      if (f2[u] == 1 && f1[u] == 0) {
        return false;
      }
    }
    return true;
  };
  auto is_representative_minimal=[&](std::vector<Tidx> const& V) -> bool {
    for (Tidx u=0; u<nbVert; u++) {
      f1[u] = 0;
    }
    for (auto & val : V) {
      f1[val] = 1;
    }
    for (auto & elt : l_elt) {
      OnFace_inplace(f2, f1, elt);
      if (!compare_f1_f2()) {
        return false;
      }
    }
    return true;
  };
  struct Level {
    std::vector<Tidx> vect;
    std::vector<std::vector<Tidx>> l_extensions;
    size_t choice;
  };
  std::vector<Level> ListLevel;
  int level = -1;
  std::vector<std::vector<Tidx>> l_extensions;
  auto compute_extensions=[&](std::vector<Tidx> const& v) {
    for (Tidx u=0; u<nbVert; u++) {
      f3[u] = 0;
    }
    for (auto & eVert : v) {
      f3[eVert] = 1;
    }
    size_t len = v.size();
    auto iife_first_element=[&]() -> Tidx {
      if (len == 0) {
        return 0;
      }
      return v[len-1] + 1;
    };
    Tidx first_elt = iife_first_element();
    std::vector<Tidx> w = v;
    w.push_back(miss_val);
    l_extensions.clear();
    for (Tidx u=first_elt; u<nbVert; u++) {
      w[len] = u;
      if (is_representative_minimal(w)) {
        l_extensions.push_back(w);
      }
    }
  };
  auto GoUpNextInTree=[&]() -> bool {
    while(true) {
      size_t choice = ListLevel[level].choice;
      size_t len = ListLevel[level].l_extensions.size();
      if (choice + 1 < len) {
        ListLevel[level].vect = ListLevel[level].l_extensions[choice + 1];
        ListLevel[level].choice = choice + 1;
        return true;
      } else {
        if (level == 0) {
          return false;
        }
        level -= 1;
        // Iterating
      }
    }
  };
  auto NextInTree=[&]() -> bool {
    if (level == -1) {
      // Start from 0
      std::vector<Tidx> v;
      compute_extensions(v);
      ListLevel[0] = {l_extensions[0], l_extensions, 0};
      level = 0;
      return true;
    } else {
      std::vector<Tidx> const& vect = ListLevel[level].vect;
      if (!f_extensible(vect)) {
        return GoUpNextInTree();
      }
      compute_extensions(vect);
      if (l_extensions.size() == 0) {
        return GoUpNextInTree();
      }
      level += 1;
      int len = ListLevel.size();
      if (len <= level) {
        ListLevel.push_back({l_extensions[0], l_extensions, 0});
      } else {
        ListLevel[level] = {l_extensions[0], l_extensions, 0};
      }
      return true;
    }
  };
  while(true) {
    bool test = NextInTree();
    if (!test) {
      break;
    }
  }
}

// clang-format off
#endif  // SRC_GROUP_GRP_ORBITENUMERATION_H_
// clang-format on

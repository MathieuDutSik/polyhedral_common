// Copyright (C) 2024 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_GRP_ORBITENUMERATION_H_
#define SRC_GROUP_GRP_ORBITENUMERATION_H_

// clang-format off
#include "GRP_GroupFct.h"
#include <limits>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ORBIT_ENUMERATION
#endif

/*
  Algorithm for generating orbits of subsets.
  ---The minimal algorithm is a tree search that passes through all lexicographically minimal elements.
  ---The canonical form uses the canonical form but does a storing of generated data,
  and so we do not do an ordered enumeration.

  Input:
  ---GRP: The permutation group acting on the elements.
  ---f_extensible: This function should be passed just once for each orbit representative:
     ---It receives a vector as input (and it can anything it wants with it like storing it).
     ---It returns true if it is extensible to a bigger cell.
 */



/*
  We are iterating by finding the minimal orbit.
  This relies on the group being relatively small.
  This is a tree search.
 */
template <typename Tgroup, typename Fextensible>
void SubsetOrbitEnumeration_minimal(Tgroup const &GRP, Fextensible f_extensible, [[maybe_unused]] std::ostream& os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx miss_val = std::numeric_limits<Tidx>::max();
#ifdef DEBUG_ORBIT_ENUMERATION
  os << "OE: SubsetOrbitEnumeration_minimal, begin\n";
#endif
  Tidx nbVert = GRP.n_act();
  std::vector<Telt> l_elt;
  for (auto &elt : GRP) {
    l_elt.push_back(elt);
  }
#ifdef DEBUG_ORBIT_ENUMERATION
  os << "OE: |l_elt|=" << l_elt << "\n";
#endif
  Face f1(nbVert), f2(nbVert), f3(nbVert);
  auto compare_f1_f2 = [&]() -> bool {
    // returns true if f1 <= f2
    for (Tidx u = 0; u < nbVert; u++) {
      if (f1[u] == 1 && f2[u] == 0) {
        return true;
      }
      if (f2[u] == 1 && f1[u] == 0) {
        return false;
      }
    }
    return true;
  };
  auto is_representative_minimal = [&](std::vector<Tidx> const &V) -> bool {
    for (Tidx u = 0; u < nbVert; u++) {
      f1[u] = 0;
    }
    for (auto &val : V) {
      f1[val] = 1;
    }
    for (auto &elt : l_elt) {
      OnFace_inplace(f2, f1, elt);
      if (!compare_f1_f2()) {
#ifdef DEBUG_ORBIT_ENUMERATION
        os << "OE: is_representative_minimal, return false\n";
#endif
        return false;
      }
    }
#ifdef DEBUG_ORBIT_ENUMERATION
    os << "OE: is_representative_minimal, return true\n";
#endif
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
  auto compute_extensions = [&](std::vector<Tidx> const &v) {
    for (Tidx u = 0; u < nbVert; u++) {
      f3[u] = 0;
    }
    for (auto &eVert : v) {
      f3[eVert] = 1;
    }
    size_t len = v.size();
#ifdef DEBUG_ORBIT_ENUMERATION
    os << "OE: len=" << len << " level=" << level << "\n";
#endif
    auto iife_first_element = [&]() -> Tidx {
      if (len == 0) {
        return 0;
      }
      return v[len - 1] + 1;
    };
    Tidx first_elt = iife_first_element();
#ifdef DEBUG_ORBIT_ENUMERATION
    os << "OE: first_elt=" << static_cast<int>(first_elt) << "\n";
#endif
    std::vector<Tidx> w = v;
    w.push_back(miss_val);
    l_extensions.clear();
    for (Tidx u = first_elt; u < nbVert; u++) {
#ifdef DEBUG_ORBIT_ENUMERATION
      os << "OE: u=" << static_cast<int>(u) << "\n";
#endif
      w[len] = u;
      if (is_representative_minimal(w)) {
#ifdef DEBUG_ORBIT_ENUMERATION
        os << "OE: w =";
        for (auto & val : w) {
          os << " " << static_cast<int>(val);
        }
        os << "\n";
#endif
        l_extensions.push_back(w);
      }
    }
  };
  auto GoUpNextInTree = [&]() -> bool {
    while (true) {
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
  auto NextInTree = [&]() -> bool {
#ifdef DEBUG_ORBIT_ENUMERATION
    os << "OE: NextInTree, level=" << level << "\n";
#endif
    if (level == -1) {
      // Start from 0
      std::vector<Tidx> v;
      compute_extensions(v);
#ifdef DEBUG_ORBIT_ENUMERATION
      os << "OE: We have 1, |l_extensions|=" << l_extensions.size() << "\n";
#endif
      if (ListLevel.size() == 0) {
        ListLevel.push_back({l_extensions[0], l_extensions, 0});
      } else {
        ListLevel[0] = {l_extensions[0], l_extensions, 0};
      }
      level = 0;
      return true;
    } else {
      std::vector<Tidx> const &vect = ListLevel[level].vect;
      if (!f_extensible(vect)) {
        return GoUpNextInTree();
      }
      compute_extensions(vect);
#ifdef DEBUG_ORBIT_ENUMERATION
      os << "OE: We have 2, |l_extensions|=" << l_extensions.size() << "\n";
#endif
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
  while (true) {
    bool test = NextInTree();
    if (!test) {
      break;
    }
  }
}

template <typename Tgroup, typename Fextensible>
void SubsetOrbitEnumeration_canform(Tgroup const &GRP, Fextensible f_extensible, [[maybe_unused]] std::ostream& os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n = GRP.n_act();
  {
    // A limit case that has to be covered.
    std::vector<Tidx> v_test;
    if (!f_extensible(v_test)) {
      return;
    }
  }
  // Generating all the elements
  std::vector<std::vector<Tidx>> list_prev;
  bool is_first = true;
  while(true) {
    if (is_first) {
      vectface vf = DecomposeOrbitPoint_Full(GRP);
      for (auto & eFace : vf) {
        std::vector<Tidx> v;
        boost::dynamic_bitset<>::size_type aRow = eFace.find_first();
        v.push_back(aRow);
        if (f_extensible(v)) {
          list_prev.push_back(v);
        }
      }
      is_first = false;
    } else {
      std::unordered_set<std::vector<Tidx>> set_new;
      for (auto & v : list_prev) {
        Face f_face(n);
        for (auto & val : v) {
          f_face[val] = 1;
        }
        for (Tidx u=0; u<n; u++) {
          if (f_face[u] == 0) {
            Face fins = f_face;
            fins[u] = 1;
            Face f2 = GRP.CanonicalImage(fins);
            std::vector<Tidx> v2 = FaceToVector<Tidx>(f2);
            set_new.insert(v2);
          }
        }
      }
      list_prev.clear();
      for (auto & v : set_new) {
        if (f_extensible(v)) {
          list_prev.push_back(v);
        }
      }
    }
    if (list_prev.size() == 0) {
      return;
    }
  }

}

template <typename Tgroup, typename Fextensible>
void SubsetOrbitEnumeration(Tgroup const &GRP, Fextensible f_extensible, [[maybe_unused]] std::ostream& os) {
  if (GRP.size() > 80000) {
    return SubsetOrbitEnumeration_canform<Tgroup,Fextensible>(GRP, f_extensible, os);
  } else {
    return SubsetOrbitEnumeration_minimal<Tgroup,Fextensible>(GRP, f_extensible, os);
  }
}


// clang-format off
#endif  // SRC_GROUP_GRP_ORBITENUMERATION_H_
// clang-format on

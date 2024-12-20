// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POINCARE_POLYHEDRON_MATRIX_GROUP_WORD_H_
#define SRC_POINCARE_POLYHEDRON_MATRIX_GROUP_WORD_H_

// clang-format off
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

// The data structure for keeping track of the group elements:
// ---Positive value (so element 0 correspond to 1, X to X+1, ...) are the
// elements themselves.
// ---Negative values correspond to their inverse (that is -1 correspond to
// inverse of generator 0)
// ---0 should never show up in the list.
// DI stands for "Direct or Inverse"
struct TrackGroup {
  std::vector<int> ListDI;
};

// We do operations but we can keep track of what is happening.
template <typename T> struct CombElt {
  TrackGroup tg;
  MyMatrix<T> mat;
  MyMatrix<double> mat_d;
};

namespace boost::serialization {
template <class Archive>
inline void serialize(Archive &ar, TrackGroup &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("ListDI", eRec.ListDI);
}

template <class Archive, typename T>
inline void serialize(Archive &ar, CombElt<T> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("tg", eRec.tg);
  ar &make_nvp("mat", eRec.mat);
  ar &make_nvp("mat_d", eRec.mat_d);
}
} // namespace boost::serialization

TrackGroup ProductTrack(TrackGroup const &tg1, TrackGroup const &tg2) {
  std::vector<int> ListDI = tg1.ListDI;
  size_t len = ListDI.size();
  for (auto &eVal : tg2.ListDI) {
    if (len > 0) {
      if (ListDI[len - 1] == -eVal) {
        ListDI.pop_back();
        len--;
      } else {
        ListDI.push_back(eVal);
        len++;
      }
    } else {
      ListDI.push_back(eVal);
      len++;
    }
  }
  return {ListDI};
}

TrackGroup InverseTrack(TrackGroup const &tg) {
  std::vector<int> ListDI_ret;
  size_t len = tg.ListDI.size();
  for (size_t u = 0; u < len; u++) {
    size_t v = len - 1 - u;
    ListDI_ret.push_back(-tg.ListDI[v]);
  }
  return {ListDI_ret};
}

template <typename T>
CombElt<T> ProductComb(CombElt<T> const &p1, CombElt<T> const &p2) {
  TrackGroup tg = ProductTrack(p1.tg, p2.tg);
  MyMatrix<T> mat = p1.mat * p2.mat;
  MyMatrix<double> mat_d = p1.mat_d * p2.mat_d;
  return {tg, mat, mat_d};
}

template <typename T> CombElt<T> InverseComb(CombElt<T> const &p) {
  MyMatrix<T> eInv = Inverse(p.mat);
  MyMatrix<double> eInv_d = UniversalMatrixConversion<double, T>(eInv);
  return {InverseTrack(p.tg), std::move(eInv), std::move(eInv_d)};
}

template <typename T> CombElt<T> GenerateIdentity(int const &n) {
  TrackGroup tg;
  MyMatrix<T> mat = IdentityMat<T>(n);
  MyMatrix<double> mat_d = IdentityMat<double>(n);
  return {tg, mat, mat_d};
}

void WriteTrackGroup(std::ofstream &os, TrackGroup const &tg) {
  size_t n_elt = tg.ListDI.size();
  os << n_elt;
  for (size_t i_elt = 0; i_elt < n_elt; i_elt++) {
    os << " " << tg.ListDI[i_elt];
  }
  os << "\n";
}

TrackGroup ReadTrackGroup(std::istream &is) {
  int n_elt;
  is >> n_elt;
  std::vector<int> ListDI;
  for (int i = 0; i < n_elt; i++) {
    int val;
    is >> val;
    ListDI.push_back(val);
  }
  return {ListDI};
}

template <typename T>
void WriteComb(std::ofstream &os, CombElt<T> const &eElt) {
  WriteTrackGroup(os, eElt.tg);
  WriteMatrix(os, eElt.mat);
}

template <typename T> CombElt<T> ReadComb(std::istream &is) {
  TrackGroup tg = ReadTrackGroup(is);
  MyMatrix<T> mat = ReadMatrix<T>(is);
  MyMatrix<double> mat_d = UniversalMatrixConversion<double, T>(mat);
  return {tg, mat, mat_d};
}

template <typename T>
bool operator==(CombElt<T> const &pe1, CombElt<T> const &pe2) {
  return pe1.mat == pe2.mat;
}

namespace std {
template <typename T> struct hash<CombElt<T>> {
  std::size_t operator()(CombElt<T> const &pe) const {
    return std::hash<MyMatrix<T>>()(pe.mat);
  }
};
} // namespace std

template <typename T> T NormCombElt(CombElt<T> const &e_elt) {
  T val = 0;
  int n = e_elt.mat.rows();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      T delta = 0;
      if (i == j)
        delta = 1;
      T u = e_elt.mat(i, j) - delta;
      val += T_abs(u);
    }
  }
  return val;
}

template <typename T>
std::vector<CombElt<T>>
InverseSaturation(std::vector<CombElt<T>> const &l_ent) {
  std::unordered_set<CombElt<T>> s_sat;
  for (auto &eElt : l_ent) {
    s_sat.insert(eElt);
    CombElt<T> eEltInv = InverseComb(eElt);
    s_sat.insert(eEltInv);
  }
  std::vector<CombElt<T>> l_ret;
  for (auto &eElt : s_sat)
    l_ret.push_back(eElt);
  return l_ret;
}

template <typename T>
std::vector<CombElt<T>> ListExpansion(std::vector<CombElt<T>> const &l_previous,
                                      std::vector<CombElt<T>> const &l_gen) {
  std::unordered_set<CombElt<T>> s_expand;
  for (auto &eElt : l_previous) {
    for (auto &fElt : l_gen) {
      CombElt<T> newElt = ProductComb(eElt, fElt);
      s_expand.insert(newElt);
    }
  }
  std::vector<CombElt<T>> l_ret;
  for (auto &eElt : s_expand)
    l_ret.push_back(eElt);
  return l_ret;
}

template <typename T> struct ShortVectorGroup {
  MyVector<T> x;
  std::vector<CombElt<T>> ListGen;
  ShortVectorGroup(MyVector<T> const &_x,
                   std::vector<CombElt<T>> const &_ListGen)
      : x(_x), ListGen(_ListGen) {}

  bool IsSolution(MyVector<T> const &v, T const &target_scal,
                  CombElt<T> const &eElt) const {
    MyVector<T> xImg = eElt.mat.transpose() * x;
    T scal = v.dot(xImg);
    return scal < target_scal;
  }

  CombElt<T> GetShortVectorNoDuplication(MyVector<T> const &y,
                                         T const &target_scal) const {
    if (ListGen.size() == 0) {
      std::cerr << "The number of generators is 0\n";
      std::cerr << "Therefore calling the function does not make sense 1\n";
      throw TerminalException{1};
    }
    HumanTime time;
    std::cerr << "Beginning of GetShortVectorNoDuplication\n";
    std::unordered_set<MyVector<T>> set_done;
    std::unordered_map<MyVector<T>, std::vector<size_t>> list_active;
    list_active[x] = std::vector<size_t>();
    set_done.insert(x);
    size_t nGen = ListGen.size();
    int n = y.size();
    int iter = 0;
    int n_cons = 0;
    while (true) {
      std::cerr << "iter=" << iter << " |list_active|=" << list_active.size()
                << "\n";
      std::unordered_map<MyVector<T>, std::vector<size_t>> list_curr =
          std::move(list_active);
      for (auto &kv : list_curr) {
        for (size_t iGen = 0; iGen < nGen; iGen++) {
          CombElt<T> const &eGen = ListGen[iGen];
          MyVector<T> xNew = eGen.mat.transpose() * kv.first;
          std::vector<size_t> eList = kv.second;
          eList.push_back(iGen);
          T scal = xNew.dot(y);
          if (scal < target_scal) {
            CombElt<T> RetElt = GenerateIdentity<T>(n);
            for (auto &pos : eList) {
              RetElt = ProductComb(RetElt, ListGen[pos]);
            }
            std::cerr << "MGW: GetShortVectorNoDuplication n_cons=" << n_cons << "\n";
            std::cerr << "|MGW: GetShortVectorNoDuplication|=" << time << "\n";
            return RetElt;
          }
          n_cons++;
          if (set_done.count(xNew) == 0) {
            list_active[xNew] = eList;
            set_done.insert(xNew);
          }
        }
      }
      iter += 1;
    }
  }

  CombElt<T> GetShortVectorIteration(MyVector<T> const &y,
                                     T const &target_scal) const {
    if (ListGen.size() == 0) {
      std::cerr << "The number of generators is 0\n";
      std::cerr << "Therefore calling the function does not make sense 2\n";
      throw TerminalException{1};
    }
    HumanTime time;
    std::cerr << "Beginning of GetShortVectorIteration\n";
    size_t nGen = ListGen.size();
    int n = y.size();
    int n_cons = 0;
    int n_iter = 1;
    while (true) {
      std::cerr << "Passing by the loop for n_iter=" << n_iter << "\n";
      std::vector<MyVector<T>> ListX(n_iter + 1);
      std::vector<size_t> eList(n_iter, 0);
      ListX[0] = x;
      for (int i = 0; i < n_iter; i++) {
        CombElt<T> const &eGen = ListGen[0];
        MyVector<T> xNew = eGen.mat.transpose() * ListX[i];
        ListX[i + 1] = xNew;
      }
      auto get_iter = [&]() -> int {
        for (int i = n_iter - 1; i >= 0; i--) {
          if (eList[i] < nGen - 1) {
            return i;
          }
        }
        return -1;
      };
      while (true) {
        int pos = get_iter();
        if (pos == -1) {
          break;
        }
        eList[pos]++;
        for (int i = pos + 1; i < n_iter; i++) {
          eList[i] = 0;
        }
        for (int i = pos; i < n_iter; i++) {
          CombElt<T> const &eGen = ListGen[eList[i]];
          MyVector<T> xNew = eGen.mat.transpose() * ListX[i];
          ListX[i + 1] = xNew;
        }
        MyVector<T> const &xTest = ListX[n_iter];
        T scal = xTest.dot(y);
        if (scal < target_scal) {
          CombElt<T> RetElt = GenerateIdentity<T>(n);
          for (auto &pos : eList) {
            RetElt = ProductComb(RetElt, ListGen[pos]);
          }
          std::cerr << "MGW: Exiting GetShortVectorIteration n_cons=" << n_cons << "\n";
          std::cerr << "|MGW: GetShortVectorIteration|=" << time << "\n";
          return RetElt;
        }
        n_cons++;
      }
      n_iter++;
    }
  }

  CombElt<T> GetShortVector(MyVector<T> const &y, T const &target_scal) const {
    //    CombElt<T> eElt1 = GetShortVectorNoDuplication(y, target_scal);
    CombElt<T> eElt2 = GetShortVectorIteration(y, target_scal);
    return eElt2;
  }
};

// As it happens, when we find some element we are likely to find
// the same one later. So testing the ones we already have is a good
// idea
template <typename T> struct ShortVectorGroupMemoize {
  ShortVectorGroup<T> const &svg;
  std::vector<CombElt<T>> ListMiss;
  ShortVectorGroupMemoize(ShortVectorGroup<T> const &_svg) : svg(_svg) {}

  void ComputeInsertSolution(MyVector<T> const &y, T const &target_scal) {
    for (auto &eElt : ListMiss) {
      if (svg.IsSolution(y, target_scal, eElt)) {
        return;
      }
    }
    CombElt<T> eElt = svg.GetShortVector(y, target_scal);
    ListMiss.push_back(eElt);
  }

  std::vector<CombElt<T>> const &GetListMiss() const { return ListMiss; }
};

// clang-format off
#endif  // SRC_POINCARE_POLYHEDRON_POINCARE_POLYHEDRON_H_
// clang-format on

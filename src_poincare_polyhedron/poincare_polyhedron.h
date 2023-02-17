// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POINCARE_POLYHEDRON_POINCARE_POLYHEDRON_H_
#define SRC_POINCARE_POLYHEDRON_POINCARE_POLYHEDRON_H_
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "QuadField.h"
#include "Namelist.h"
#include "POLY_DirectDualDesc.h"
#include "POLY_PolytopeFct.h"
#include "POLY_LinearProgramming.h"
// clang-format on

//
// First part, the data set used
//

// The data structure for keeîng track of the group elements:
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
template <typename T> struct PairElt {
  TrackGroup tg;
  MyMatrix<T> mat;
};

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
PairElt<T> ProductPair(PairElt<T> const &p1, PairElt<T> const &p2) {
  TrackGroup tg = ProductTrack(p1.tg, p2.tg);
  MyMatrix<T> mat = p1.mat * p2.mat;
  return {tg, mat};
}

template <typename T> PairElt<T> InversePair(PairElt<T> const &p) {
  return {InverseTrack(p.tg), Inverse(p.mat)};
}

template <typename T> PairElt<T> GenerateIdentity(int const &n) {
  TrackGroup tg;
  MyMatrix<T> mat = IdentityMat<T>(n);
  return {tg, mat};
}

void WriteTrackGroup(std::ofstream &os, TrackGroup const &tg) {
  size_t n_elt = tg.ListDI.size();
  os << n_elt;
  for (size_t i_elt = 0; i_elt < n_elt; i_elt++) {
    os << ":" << tg.ListDI[i_elt];
  }
}

template <typename T>
bool operator==(PairElt<T> const &pe1, PairElt<T> const &pe2) {
  return pe1.mat == pe2.mat;
}

namespace std {
template <typename T> struct hash<PairElt<T>> {
  std::size_t operator()(PairElt<T> const &pe) const {
    return std::hash<MyMatrix<T>>()(pe.mat);
  }
};
} // namespace std

//
// Second part, the finite group
//

template <typename T>
std::vector<PairElt<T>>
GroupGeneration(std::vector<PairElt<T>> const &input_l_ent) {
  std::vector<PairElt<T>> l_ent = input_l_ent;
  while (true) {
    std::unordered_set<PairElt<T>> GenElt;
    size_t n_ent = l_ent.size();
    for (size_t i_ent = 0; i_ent < n_ent; i_ent++) {
      PairElt<T> const &pe1 = l_ent[i_ent];
      for (size_t j_ent = 0; j_ent < n_ent; j_ent++) {
        PairElt<T> const &pe2 = l_ent[j_ent];
        PairElt<T> prod = ProductPair(pe1, pe2);
        GenElt.insert(prod);
      }
    }
    if (GenElt.size() == n_ent) {
      return input_l_ent;
    }
    l_ent.clear();
    for (auto &e_ent : GenElt)
      l_ent.push_back(e_ent);
  }
}

// a right coset is of the form Ug
template <typename T>
std::vector<PairElt<T>>
IdentifyRightCosets(std::vector<PairElt<T>> const &l_ent,
                    std::vector<PairElt<T>> const &list_grp_elt) {
  std::unordered_set<PairElt<T>> s_coset;
  auto f_insert = [&](PairElt<T> const &pe) -> void {
    for (auto &e_grp_elt : list_grp_elt) {
      PairElt<T> prod = ProductPair(e_grp_elt, pe);
      if (s_coset.count(prod) == 1)
        break;
    }
    s_coset.insert(pe);
  };
  for (auto &pe : l_ent)
    f_insert(pe);
  std::vector<PairElt<T>> l_coset;
  for (auto &e_coset : s_coset)
    l_coset.push_back(e_coset);
  return l_coset;
}

// a left coset is of the form gU
template <typename T>
std::vector<PairElt<T>>
IdentifyLeftCosets(std::vector<PairElt<T>> const &l_ent,
                   std::vector<PairElt<T>> const &list_grp_elt) {
  std::unordered_set<PairElt<T>> s_coset;
  auto f_insert = [&](PairElt<T> const &pe) -> void {
    for (auto &e_grp_elt : list_grp_elt) {
      PairElt<T> prod = ProductPair(pe, e_grp_elt);
      if (s_coset.count(prod) == 1)
        break;
    }
    s_coset.insert(pe);
  };
  for (auto &pe : l_ent)
    f_insert(pe);
  std::vector<PairElt<T>> l_coset;
  for (auto &e_coset : s_coset)
    l_coset.push_back(e_coset);
  return l_coset;
}

template <typename T>
std::vector<PairElt<T>>
InverseSaturation(std::vector<PairElt<T>> const &l_ent) {
  std::unordered_set<PairElt<T>> s_sat;
  int i_ent = 0;
  for (auto &eElt : l_ent) {
    s_sat.insert(eElt);
    PairElt<T> eEltInv = InversePair(eElt);
    s_sat.insert(eEltInv);
    i_ent++;
  }
  std::vector<PairElt<T>> l_ret;
  for (auto &eElt : s_sat)
    l_ret.push_back(eElt);
  return l_ret;
}

template <typename T>
std::vector<PairElt<T>> ListExpansion(std::vector<PairElt<T>> const &l_previous,
                                      std::vector<PairElt<T>> const &l_gen) {
  std::unordered_set<PairElt<T>> s_expand;
  for (auto &eElt : l_previous) {
    for (auto &fElt : l_gen) {
      PairElt<T> newElt = ProductPair(eElt, fElt);
      s_expand.insert(newElt);
    }
  }
  std::vector<PairElt<T>> l_ret;
  for (auto &eElt : s_expand)
    l_ret.push_back(eElt);
  return l_ret;
}

//
// The common function for paperwork
//

template <typename T> struct DataEXT {
  MyMatrix<T> EXT;
  std::vector<Face> v_red;
};

template <typename T>
DataEXT<T> GetTransposedDualDesc(vectface const &vf, MyMatrix<T> const &FAC) {
  int n = FAC.cols();
  int n_fac = FAC.rows();
  int n_ext = vf.size();
  MyMatrix<T> EXT(n_ext, n);
  int pos = 0;
  std::vector<Face> v_red(n_fac, Face(n_ext));
  for (auto &eInc : vf) {
    MyVector<T> eEXT = FindFacetInequality(FAC, eInc);
    AssignMatrixRow(EXT, pos, eEXT);
    for (int i_fac = 0; i_fac < n_fac; i_fac++) {
      v_red[i_fac][pos] = eInc[i_fac];
    }
    pos++;
  }
  return {EXT, v_red};
}

std::unordered_map<Face, size_t>
get_map_face(std::vector<Face> const &l_facet) {
  size_t n_facet = l_facet.size();
  std::unordered_map<Face, size_t> s_facet;
  for (size_t i_facet = 0; i_facet < n_facet; i_facet++) {
    Face const &f = l_facet[i_facet];
    s_facet[f] = i_facet + 1;
  }
  if (l_facet.size() != s_facet.size()) {
    std::cerr << "l_facet contains some duplicate and that is illegal\n";
    throw TerminalException{1};
  }
  return s_facet;
}

//
// The elementary data structures for the Poincare stuff
//

// The initial data for the Poincare Polyhedron Theorem
// ---a point x
// ---a list of group element which is of course assumed to generate the group
template <typename T> struct DataPoincare {
  MyVector<T> x;
  std::vector<PairElt<T>> ListGroupElt;
};

template <typename T>
DataPoincare<T> ReadDataPoincare(std::string const &FileI,
                                 int const &n_expand) {
  IsExistingFileDie(FileI);
  std::ifstream is(FileI);
  MyVector<T> x = ReadVector<T>(is);
  int n = x.size();
  std::cerr << "ReadDataPoincare : |x|=" << n << "\n";
  int n_elt;
  is >> n_elt;
  std::cerr << "ReadDataPoincare : n_elt=" << n_elt << "\n";
  std::unordered_set<PairElt<T>> s_elt;
  s_elt.insert(GenerateIdentity<T>(n));
  for (int i_elt = 0; i_elt < n_elt; i_elt++) {
    int pos = i_elt + 1;
    MyMatrix<T> eElt = ReadMatrix<T>(is);
    T TheDet = DeterminantMat(eElt);
    std::cerr << "i_elt=" << i_elt << " TheDet=" << TheDet << "\n";
    TrackGroup tg{{pos}};
    PairElt<T> pe{tg, eElt};
    s_elt.insert(pe);
  }
  std::vector<PairElt<T>> l_elt;
  for (auto &eElt : s_elt)
    l_elt.push_back(eElt);
  std::cerr << "|Initial set|=" << l_elt.size() << "\n";
  l_elt = InverseSaturation(l_elt);
  std::vector<PairElt<T>> l_gen = l_elt;
  std::cerr << "|Inverse saturation|=" << l_elt.size() << "\n";
  std::cerr << "n_expand=" << n_expand << "\n";
  for (int i_expand = 0; i_expand < n_expand; i_expand++) {
    l_elt = ListExpansion(l_elt, l_gen);
    std::cerr << "i_expand=" << i_expand << " |l_elt|=" << l_elt.size() << "\n";
  }
  return {x, l_elt};
}

// This is for the single adjacency in the polyhedron
// * iFaceAdj is the index of the facet in the polyhedron
//   which is adjacent to it in the ridge.
// * iPolyAdj is the index in the adjacent iFaceAdj face
//   in the ridge.
// * iFaceOpp is the index of the facet in the mapped polyhedron
// * iPolyAdj is similarly the corresponding index of the ridge
// * EXTadj is the record of the vertices of the codimension 2 face
struct TsingAdj {
  size_t iFaceAdj;
  size_t iPolyAdj;
  size_t iFaceOpp;
  size_t iPolyOpp;
  Face IncdRidge;
};

struct Tfacet {
  Face IncdFacet;
  std::vector<TsingAdj> l_sing_adj;
};

template <typename T> struct ResultAdjacencyInfo {
  MyMatrix<T> FAC;
  MyMatrix<T> EXT;
  std::unordered_map<Face, size_t> s_facet;
};

template <typename T> struct AdjacencyInfo {
  MyMatrix<T> EXT;
  std::vector<Tfacet> ll_adj;
};

template <typename T> MyMatrix<T> Contragredient(MyMatrix<T> const &M) {
  return Inverse(TransposedMat(M));
}

/*
  The program works for A3 (where actually we used B3)
  But it has some problems for the case of the cocompact group of interest to us.
  ---
  Problems:
  - Bad numerics where we are forced to have evaluation of number very near 0 and the continued
    fraction algorithm fails us.
    => Look at the absolute values of the coefficients in double precision.
    => How hard it is to compute the symmetric function P(x_1) ... P(x_k)?
       This ought to give us a value that is rational. If P(x_1) is very small, likely the others are not.
    => See exactly where things go haywire.
  - Non-convergence of finding the dual description. Our hope is definitely that there are a lot of
    facets that are redundant. This is the foundation of our approach. Another element is that for an element,
    we can easily identify if it is new or already present.
    So, we have already accepted that insertion done step by step is the business of the day.
  - If the checks of non-redundancy are so so expensive, then we should keep track of them
    in a database of known redundant entries.
    This should probably be just a std::unordered_map<Myvector<T>,std::pair<size_t,size_t>>
    We would have a statbilizer type that is passed as a reference to the insertion process.
    and then when an insertion is done.
    But all that seems overkill right now because the stabilizer are so far trivial.
  - What should be reasonably accessible for us:
    - redundancy checks using Clarkson method (based on Linear programming)
    - Linear programming over those special fields has to be doable.
      Maybe we could use shifting to doule precision to help solve those linear problems.
    - dual description for a set of inequalities that are all defining facets.

 */
template <typename T> struct StepEnum {
public:
  MyVector<T> x;
  std::vector<PairElt<T>> stabilizerElt;
  std::unordered_set<PairElt<T>> stabilizerElt_map;
  std::vector<PairElt<T>> ListNeighborCoset;
  std::vector<MyVector<T>> ListNeighborX;
  std::vector<std::pair<size_t, size_t>> ListNeighborData;
  std::unordered_map<MyVector<T>, std::pair<size_t, size_t>> map;
  void print_statistics(std::ostream& os) const {
    os << "|stabilizerElt|=" << stabilizerElt.size() << "\n";
    os << "|ListNeighborCoset|=" << ListNeighborCoset.size() << " |ListNeighborX|=" << ListNeighborX.size() << "\n";
    std::map<size_t,size_t> map_by_cos;
    for (auto & ePair : ListNeighborData) {
      map_by_cos[ePair.first]++;
    }
    std::map<size_t,size_t> map_by_size;
    for (auto & kv: map_by_cos) {
      map_by_size[kv.second]++;
    }
    os << "Sizes :";
    for (auto & kv : map_by_size) {
      os << " [" << kv.first << "," << kv.second << "]";
    }
    os << "\n";
  }
  bool InsertStabilizerGenerator(PairElt<T> const &eElt) {
    if (IsPresentInStabilizer(eElt))
      return false;
    std::vector<PairElt<T>> ExtListGen = stabilizerElt;
    ExtListGen.push_back(eElt);
    stabilizerElt = GroupGeneration(ExtListGen);
    stabilizerElt_map.clear();
    for (auto & eElt : stabilizerElt)
      stabilizerElt_map.insert(eElt);
    return true;
  }
  PairElt<T> GetElement(std::pair<size_t, size_t> const &val) const {
    size_t i_coset = val.first;
    size_t i_elt = val.second;
    return ProductPair(ListNeighborCoset[i_coset], stabilizerElt[i_elt]);
  }
  void InsertCoset(PairElt<T> const &eCoset) {
    size_t n_coset = ListNeighborCoset.size();
    ListNeighborCoset.push_back(eCoset);
    std::unordered_map<MyVector<T>, std::pair<size_t, size_t>> map_local;
    MyVector<T> x_img = eCoset.mat.transpose() * x;
    size_t n_elt = stabilizerElt.size();
    for (size_t i_elt = 0; i_elt < n_elt; i_elt++) {
      MyVector<T> x2 = stabilizerElt[i_elt].mat.transpose() * x_img;
      std::pair<size_t, size_t> &val = map_local[x2];
      val.first = n_coset;
      val.second = i_elt;
    }
    for (auto &kv : map_local) {
      ListNeighborX.push_back(kv.first);
      ListNeighborData.push_back(kv.second);
      map[kv.first] = kv.second;
    }
  }
  void ComputeCosets(std::vector<PairElt<T>> const &l_elt) {
    ListNeighborX.clear();
    ListNeighborData.clear();
    ListNeighborCoset.clear();
    map.clear();
    std::vector<PairElt<T>> l_cos = IdentifyLeftCosets(l_elt, stabilizerElt);
    for (auto &eCoset : l_cos) {
      InsertCoset(eCoset);
    }
  }
  bool IsPresentInStabilizer(PairElt<T> const& eElt) const {
    return stabilizerElt_map.find(eElt) != stabilizerElt_map.end();
  }
  bool IsPresentInCosetOrStabilizer(PairElt<T> const& eElt) const {
    MyVector<T> x_img = eElt.mat.transpose() * x;
    if (x_img == x) {
      return IsPresentInStabilizer(eElt);
    } else {
      auto iter = map.find(x_img);
      if (iter == map.end()) {
        return false;
      } else {
        PairElt<T> TheProd = GetElement(iter->second);
        PairElt<T> eElt_inv = InversePair(eElt);
        PairElt<T> stab_elt = ProductPair(TheProd, eElt_inv);
        MyVector<T> x2 = stab_elt.mat.transpose() * x;
        if (x2 != x) {
          std::cerr << "x is not stabilized\n";
          throw TerminalException{1};
        }
        return IsPresentInStabilizer(stab_elt);
      }
    }
  }
  bool InsertGenerators(std::vector<PairElt<T>> const &ListGen) {
    bool DidSomething = false;
    auto generator_upgrade = [&](PairElt<T> const &e_elt) -> void {
      bool test = InsertStabilizerGenerator(e_elt);
      if (test) {
        // Copy needed of the old data then recompute
        std::vector<PairElt<T>> OldListCos = ListNeighborCoset;
        ComputeCosets(OldListCos);
        DidSomething = true;
      }
    };
    for (auto &e_elt : ListGen) {
      MyVector<T> x_img = e_elt.mat.transpose() * x;
      if (x_img == x) {
        generator_upgrade(e_elt);
      } else {
        auto iter = map.find(x_img);
        if (iter == map.end()) {
          InsertCoset(e_elt);
          DidSomething = true;
        } else {
          PairElt<T> TheProd = GetElement(iter->second);
          PairElt<T> e_elt_inv = InversePair(e_elt);
          PairElt<T> stab_elt = ProductPair(TheProd, e_elt_inv);
          MyVector<T> x2 = stab_elt.mat.transpose() * x;
          if (x2 != x) {
            std::cerr << "x is not stabilized\n";
            throw TerminalException{1};
          }
          generator_upgrade(stab_elt);
        }
      }
    }
    print_statistics(std::cerr);
    return DidSomething;
  }
  StepEnum(MyVector<T> const &_x) {
    x = _x;
    int n = x.size();
    PairElt<T> IdMat = GenerateIdentity<T>(n);
    stabilizerElt = {IdMat};
    stabilizerElt_map.insert(IdMat);
  }
  MyVector<T> GetIneq(PairElt<T> const e_elt) const {
    MyMatrix<T> const &eMat = e_elt.mat;
    MyVector<T> x_img = eMat.transpose() * x;
    MyVector<T> x_diff = x_img - x;
    return x_diff;
  }
  MyMatrix<T> GetFAC() const {
    int n = x.size();
    int n_mat = ListNeighborX.size();
    MyMatrix<T> FAC(n_mat, n);
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      MyVector<T> x_diff = ListNeighborX[i_mat] - x;
      AssignMatrixRow(FAC, i_mat, x_diff);
    }
    return FAC;
  }
  void RemoveRedundancy(std::string const &eCommand) {
    MyMatrix<T> FAC = GetFAC();
    int n = x.size();
    int n_mat = FAC.rows();
    int rnk = RankMat(FAC);
    std::cerr << "RemoveRedundancy : n=" << n << " n_mat=" << n_mat << " rnk=" << rnk << "\n";
    if (rnk != n) {
      std::cerr << "Error in RemoveRedundancy\n";
      std::cerr << "n=" << n << " n_mat=" << n_mat << " rnk=" << rnk << "\n";
      throw TerminalException{1};
    }
    vectface vf = DirectFacetOrbitComputation_nogroup(FAC, eCommand, std::cerr);
    std::cerr << "We have vf\n";
    DataEXT<T> dataext = GetTransposedDualDesc(vf, FAC);
    std::cerr << "RemoveRedundancy : n_ext=" << dataext.EXT.rows() << "\n";
    Face f_status_keep(n_mat);
    std::set<size_t> l_keep;
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      MyMatrix<T> EXT_red = SelectRow(dataext.EXT, dataext.v_red[i_mat]);
      std::pair<size_t, size_t> epair = ListNeighborData[i_mat];
      size_t i_coset = epair.first;
      int rnk = RankMat(EXT_red);
      if (rnk == n - 1) {
        l_keep.insert(i_coset);
        f_status_keep[i_mat] = 1;
      }
    }
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      std::pair<size_t, size_t> epair = ListNeighborData[i_mat];
      size_t i_coset = epair.first;
      if (l_keep.count(i_coset) != f_status_keep[i_mat]) {
        std::cerr << "There is incoherency in the orbit nature\n";
        throw TerminalException{1};
      }
    }
    std::cerr << "RemoveRedundancy : |l_keep|=" << l_keep.size() << "\n";
    std::vector<PairElt<T>> ListNeighborCosetRed;
    for (auto &i_coset : l_keep) {
      ListNeighborCosetRed.push_back(ListNeighborCoset[i_coset]);
    }
    ComputeCosets(ListNeighborCosetRed);
  }
  // The domain is defined originally as
  // Tr(AX) <= Tr(PAP^T X)
  // which we rewrite
  // a.x <= phi(a).x        0 <= (phi(a) - a).x
  //
  // Matrixwise the scalar product a.x is rewritten as
  // A X^T
  // phi(a) is expressed as AQ
  // 0 <= (AQ - A) X^T
  // The mapping A ----> AQ maps to the adjacent domain.
  // The mapping of the X ----> X c(Q)^T with c(Q) the
  // contragredient representation.
  AdjacencyInfo<T> ComputeAdjacencyInfo(std::string const &eCommand) {
    size_t miss_val = std::numeric_limits<size_t>::max();
    MyMatrix<T> FAC = GetFAC();
    int n = x.size();
    int n_mat = ListNeighborData.size();
    int rnk = RankMat(FAC);
    if (rnk != n) {
      std::cerr << "Error in ComputeAdjacencyInfo\n";
      std::cerr << "n=" << n << " n_mat=" << n_mat << " rnk=" << rnk << "\n";
      throw TerminalException{1};
    }
    vectface vf = DirectFacetOrbitComputation_nogroup(FAC, eCommand, std::cerr);
    DataEXT<T> dataext = GetTransposedDualDesc(vf, FAC);
    int n_ext = dataext.EXT.rows();
    std::unordered_map<Face, size_t> s_facet = get_map_face(dataext.v_red);

    std::cerr << "First part: adjacency structure within the polyhedron\n";
    std::vector<Tfacet> ll_adj;
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      Face const &f1 = dataext.v_red[i_mat];
      std::vector<TsingAdj> l_adj;
      for (int j_mat = 0; j_mat < n_mat; j_mat++) {
        if (i_mat != j_mat) {
          Face const &f2 = dataext.v_red[j_mat];
          Face f(n_ext);
          for (int i_ext = 0; i_ext < n_ext; i_ext++)
            if (f1[i_ext] == 1 && f2[i_ext] == 1)
              f[i_ext] = 1;
          MyMatrix<T> EXT_red = SelectRow(dataext.EXT, f);
          int rnk = RankMat(EXT_red);
          if (rnk == n - 2) {
            size_t j_mat_s = static_cast<size_t>(j_mat);
            l_adj.push_back({j_mat_s, miss_val, miss_val, miss_val, f});
          }
        }
      }
      Face IncdFacet = dataext.v_red[i_mat];
      ll_adj.push_back({IncdFacet, l_adj});
    }
    auto get_iPoly = [&](size_t iFace, Face const &f1) -> size_t {
      size_t n_adjB = ll_adj[iFace].l_sing_adj.size();
      for (size_t i_adjB = 0; i_adjB < n_adjB; i_adjB++) {
        Face f2 = ll_adj[iFace].l_sing_adj[i_adjB].IncdRidge;
        if (f1 == f2)
          return i_adjB;
      }
      std::cerr << "Failed to find a matching for f1\n";
      throw TerminalException{1};
    };
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      size_t n_adj = ll_adj[i_mat].l_sing_adj.size();
      for (size_t i_adj = 0; i_adj < n_adj; i_adj++) {
        size_t iFaceAdj = ll_adj[i_mat].l_sing_adj[i_adj].iFaceAdj;
        Face f1 = ll_adj[i_mat].l_sing_adj[i_adj].IncdRidge;
        ll_adj[i_mat].l_sing_adj[i_adj].iPolyAdj = get_iPoly(iFaceAdj, f1);
      }
    }
    std::cerr << "Second part: computing the opposite facets\n";
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      size_t n_adj = ll_adj[i_mat].l_sing_adj.size();
      MyMatrix<T> Q = GetElement(ListNeighborData[i_mat]).mat;
      MyMatrix<T> cQ = Contragredient(Q);
      MyMatrix<T> EXTimg = dataext.EXT * cQ;
      ContainerMatrix<T> Cont(EXTimg);
      Face MapFace(n_ext);
      std::vector<size_t> l_idx(n_ext, 0);
      std::map<int, size_t> map_index;
      for (int i_ext = 0; i_ext < n_ext; i_ext++) {
        if (ll_adj[i_mat].IncdFacet[i_ext] == 1) {
          std::optional<size_t> opt = Cont.GetIdx_f([&](int i) -> T {
            return dataext.EXT(i_ext, i);
          });
          if (opt) {
            size_t idx = *opt;
            MapFace[idx] = 1;
            map_index[i_ext] = idx;
          }
        }
      }
      size_t pos = s_facet[MapFace];
      if (pos == 0) {
        std::cerr << "Failed to find the MapFace in s_facet\n";
        throw TerminalException{1};
      }
      size_t iFaceOpp = pos - 1;
      for (size_t i_adj = 0; i_adj < n_adj; i_adj++) {
        Face f_map(n_ext);
        for (int i_ext = 0; i_ext < n_ext; i_ext++) {
          //          std::cerr << "i_mat=" << i_mat << " i_adj=" << i_adj << " i_ext=" << i_ext << "\n";
          if (ll_adj[i_mat].l_sing_adj[i_adj].IncdRidge[i_ext] == 1) {
            size_t idx = map_index[i_ext];
            f_map[idx] = 1;
          }
        }
        size_t iPolyOpp = get_iPoly(iFaceOpp, f_map);
        ll_adj[i_mat].l_sing_adj[i_adj].iFaceOpp = iFaceOpp;
        ll_adj[i_mat].l_sing_adj[i_adj].iPolyOpp = iPolyOpp;
      }
    }
    bool print_ai = false;
    if (print_ai) {
      std::cerr << "ai: n_mat=" << n_mat << "\n";
      for (int i_mat=0; i_mat<n_mat; i_mat++) {
        int n_adj = ll_adj[i_mat].l_sing_adj.size();
        std::cerr << "  i_mat=" << i_mat << " |l_sing_adj|=" << n_adj << "\n";
        for (int i_adj=0; i_adj<n_adj; i_adj++) {
          TsingAdj const& singAdj = ll_adj[i_mat].l_sing_adj[i_adj];
          std::cerr << "    i_adj=" << i_adj << " iFaceAdj=" << singAdj.iFaceAdj << " iPolyAdj=" << singAdj.iPolyAdj << " iFaceOpp=" << singAdj.iFaceOpp << " iPolyOpp=" << singAdj.iPolyOpp << "\n";
        }
      }
    }
    return {dataext.EXT, ll_adj};
  }
  std::vector<PairElt<T>>
  GenerateTypeIneighbors(std::vector<PairElt<T>> const &l_elt) {
    std::cerr << "Beginning of GenerateTypeIneighbors\n";
    int n_mat = ListNeighborX.size();
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      MyVector<T> u = ListNeighborX[i_mat];
      std::cerr << "i_mat=" << i_mat << " neighbor=" << StringVector(u) << "\n";
    }
    MyMatrix<T> FAC = GetFAC();
    std::cerr << "FAC=\n";
    WriteMatrix(std::cerr, FAC);
    auto f_insert = [&](PairElt<T> const &TestElt) -> std::optional<PairElt<T>> {
      PairElt<T> WorkElt = TestElt;
      int n_iter = 0;
      std::cerr << "Beginning of f_insert\n";
      while(true) {
        std::cerr << "  n_iter=" << n_iter << "\n";
        MyVector<T> x_ineq = GetIneq(WorkElt);
        if (IsZeroVector(x_ineq)) {
          bool test_stab = IsPresentInStabilizer(WorkElt);
          if (test_stab) {
            return {};
          } else {
            return WorkElt;
          }
        }
        std::optional<MyVector<T>> opt = SolutionMatNonnegative(FAC, x_ineq);
        if (opt) {
          MyVector<T> V = *opt;
          std::cerr << "  SolMatNonNeg V=" << StringVector(V) << "\n";
          for (int u = 0; u < V.size(); u++) {
            if (V(u) > 0) {
              PairElt<T> uElt = GetElement(ListNeighborData[u]);
              PairElt<T> uEltInv = InversePair(uElt);
              WorkElt = ProductPair(WorkElt, uEltInv);
            }
          }
        } else {
          std::cerr << "  SolMatNonNeg : no solution found\n";
          return WorkElt;
        }
        n_iter++;
      }
    };
    std::vector<PairElt<T>> ListMiss;
    for (auto &e_elt : l_elt) {
      std::optional<PairElt<T>> opt = f_insert(e_elt);
      if (opt) {
        PairElt<T> const& RedElt = *opt;
        ListMiss.push_back(RedElt);
      }
    }
    return ListMiss;
  }
  bool TestIntersection(MyMatrix<T> const &FAC, PairElt<T> const &eElt) {
    MyMatrix<T> FACimg = FAC * eElt.mat;
    MyMatrix<T> FACtot = Concatenate(FAC, FACimg);
    return IsFullDimensional(FACtot);
  }
  std::vector<PairElt<T>> GenerateTypeIIneighbors(AdjacencyInfo<T> const &ai) {
    int n = x.size();
    std::vector<PairElt<T>> ListMiss;
    MyMatrix<T> FAC = GetFAC();
    int n_mat = FAC.rows();
    std::vector<PairElt<T>> ListAdj;
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      ListAdj.push_back(GetElement(ListNeighborData[i_mat]));
    }
    auto GetMissedGenerator = [&](int i_mat,
                                  int i_facet) -> std::optional<PairElt<T>> {
      PairElt<T> TheMat = GenerateIdentity<T>(n);
      int i_mat_work = i_mat;
      int i_facet_work = i_facet;
      //      std::cerr << "i_mat_work=" << i_mat_work << " i_facet_work=" << i_facet_work << "\n";
      while (true) {
        TheMat = ProductPair(TheMat, ListAdj[i_mat_work]);
        int iFaceOpp = ai.ll_adj[i_mat_work].l_sing_adj[i_facet_work].iFaceOpp;
        //        std::cerr << "We have iFaceOpp=" << iFaceOpp << "\n";
        int iPolyOpp = ai.ll_adj[i_mat_work].l_sing_adj[i_facet_work].iPolyOpp;
        //        std::cerr << "We have iPolyOpp=" << iPolyOpp << "\n";
        i_mat_work = ai.ll_adj[iFaceOpp].l_sing_adj[iPolyOpp].iFaceAdj;
        //        std::cerr << "Now i_mat_work=" << i_mat_work << "\n";
        i_facet_work = ai.ll_adj[iFaceOpp].l_sing_adj[iPolyOpp].iPolyAdj;
        //        std::cerr << "Now i_facet_work=" << i_facet_work << "\n";
        MyVector<T> x_img = TheMat.mat.transpose() * x;
        if (x_img == x) {
          return {};
        }
        bool test = TestIntersection(FAC, TheMat);
        //        std::cerr << "test=" << test << "\n";
        if (test) {
          return TheMat;
        }
      }
    };
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      int n_facet = ai.ll_adj[i_mat].l_sing_adj.size();
      for (int i_facet = 0; i_facet < n_facet; i_facet++) {
        std::optional<PairElt<T>> opt = GetMissedGenerator(i_mat, i_facet);
        if (opt) {
          ListMiss.push_back(*opt);
        }
      }
    }
    return ListMiss;
  }
  template <typename Tgroup> Tgroup GetPermutationGroup() const {
    using Telt = typename Tgroup::Telt;
    using Tidx = typename Telt::Tidx;
    int n_neigh = ListNeighborX.size();
    std::unordered_map<MyVector<T>, int> map_rev;
    for (int i_neigh = 0; i_neigh < n_neigh; i_neigh++)
      map_rev[ListNeighborX[i_neigh]] = i_neigh + 1;
    auto GetPermutation = [&](PairElt<T> const &eElt) -> Telt {
      std::vector<Tidx> eList(n_neigh);
      for (int i_neigh = 0; i_neigh < n_neigh; i_neigh++) {
        MyVector<T> V = ListNeighborX[i_neigh];
        MyVector<T> Vimg = eElt.mat.transpose() * V;
        int pos = map_rev[Vimg];
        if (pos == 0) {
          std::cerr << "Vimg should belong to the list of entries\n";
          throw TerminalException{1};
        }
        eList[i_neigh] = pos - 1;
      }
      return Telt(eList);
    };
    std::vector<Tidx> idList(n_neigh);
    for (int i_neigh = 0; i_neigh < n_neigh; i_neigh++)
      idList[i_neigh] = i_neigh;
    Telt id(idList);
    std::vector<Telt> LGen;
    Tgroup GRP(LGen, id);
    auto f_insert = [&](PairElt<T> const &eElt) -> void {
      Telt ePerm = GetPermutation(eElt);
      if (GRP.isin(ePerm))
        return;
      LGen.push_back(ePerm);
      GRP = Tgroup(LGen, idList);
    };
    for (auto &eElt : stabilizerElt)
      f_insert(eElt);
    return GRP;
  }
  std::pair<int, std::vector<std::vector<int>>> GetGroupPresentation(AdjacencyInfo<T> const &ai) {
    int n = x.size();
    std::vector<std::vector<int>> ListWord;
    MyMatrix<T> FAC = GetFAC();
    int n_mat = FAC.rows();
    std::vector<PairElt<T>> ListAdj;
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      ListAdj.push_back(GetElement(ListNeighborData[i_mat]));
    }
    auto InsertWordRidge = [&](int i_mat, int i_facet) -> void {
      PairElt<T> TheMat = GenerateIdentity<T>(n);
      int i_mat_work = i_mat;
      int i_facet_work = i_facet;
      std::vector<int> TheWord;
      while (true) {
        TheWord.push_back(i_mat_work);
        TheMat = ProductPair(TheMat, ListAdj[i_mat_work]);
        int iFaceOpp = ai.ll_adj[i_mat_work].l_sing_adj[i_facet_work].iFaceOpp;
        int iPolyOpp = ai.ll_adj[i_mat_work].l_sing_adj[i_facet_work].iPolyOpp;
        i_mat_work = ai.ll_adj[iFaceOpp].l_sing_adj[iPolyOpp].iFaceAdj;
        i_facet_work = ai.ll_adj[iFaceOpp].l_sing_adj[iPolyOpp].iPolyAdj;
        MyVector<T> x_img = TheMat.mat.transpose() * x;
        if (x_img == x) {
          ListWord.push_back(TheWord);
          return;
        }
      }
    };
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      int iFaceOpp = ai.ll_adj[i_mat].l_sing_adj[0].iFaceOpp;
      std::vector<int> TheWord{i_mat, iFaceOpp};
      ListWord.push_back(TheWord);
    }
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      int n_facet = ai.ll_adj[i_mat].l_sing_adj.size();
      for (int i_facet = 0; i_facet < n_facet; i_facet++) {
        InsertWordRidge(i_mat, i_facet);
      }
    }
    return {n_mat, ListWord};
  }
};

//
// Now the polyhedral stuff
//

template <typename T>
std::optional<StepEnum<T>> ComputeMissingNeighbors(AdjacencyInfo<T> const &ai) {
  return {};
}

struct RecOption {
  std::string eCommand;
  std::string FileI;
  std::string FileO;
  std::string Arithmetic;
  int n_iter_max;
  int n_expand;
  bool ComputeStabilizerPermutation;
  bool ComputeGroupPresentation;
};

template <typename T>
StepEnum<T> IterativePoincareRefinement(DataPoincare<T> const &dp,
                                        RecOption const &rec_option) {
  std::string eCommand = rec_option.eCommand;
  StepEnum<T> se(dp.x);
  se.InsertGenerators(dp.ListGroupElt);
  se.RemoveRedundancy(eCommand);
  bool DidSomething = false;
  auto insert_block = [&](std::vector<PairElt<T>> const &ListMiss) -> void {
    if (ListMiss.size() > 0) {
      bool test = se.InsertGenerators(ListMiss);
      if (test) {
        se.RemoveRedundancy(eCommand);
        DidSomething = true;
      }
    }
  };
  int n_iter = 0;
  while (true) {
    DidSomething = false;
    //
    // Iteration Type I
    //
    std::cerr << "IterativePoincareRefinement n_iter=" << n_iter << "\n";
    std::vector<PairElt<T>> ListMissI =
        se.GenerateTypeIneighbors(dp.ListGroupElt);
    std::cerr << "|ListMissI|=" << ListMissI.size() << "\n";
    insert_block(ListMissI);
    //
    // Iteration Type II
    //
    AdjacencyInfo<T> ai = se.ComputeAdjacencyInfo(eCommand);
    std::vector<PairElt<T>> ListMissII = se.GenerateTypeIIneighbors(ai);
    std::cerr << "|ListMissII|=" << ListMissII.size() << "\n";
    insert_block(ListMissII);
    //
    // Terminating if ok.
    //
    if (!DidSomething) {
      return se;
    }
    //
    // Iteration checks
    //
    if (rec_option.n_iter_max > 0) {
      if (n_iter > rec_option.n_iter_max) {
        std::cerr << "Reached the maximum number of iterations for type I\n";
        throw TerminalException{1};
      }
    }
    n_iter++;
  }
}

//
// The code for calling it.
//

FullNamelist NAMELIST_GetPoincareInput() {
  std::map<std::string, SingleBlock> ListBlock;
  // METHOD
  std::map<std::string, std::string> ListBoolValues_doc;
  std::map<std::string, std::string> ListIntValues_doc;
  std::map<std::string, std::string> ListStringValues_doc;
  ListBoolValues_doc["ComputeStabilizerPermutation"] = "Default: F\n\
Compute the action of its stabilizer on the facets as permutation group";
  ListBoolValues_doc["ComputeGroupPresentation"] = "Default: F\n\
Compute the presentation of the greoup from facets and ridges";
  ListIntValues_doc["n_iter_max"] = "Default: -1\n\
The maximum number of iteration. If negative then infinite";
  ListIntValues_doc["n_expand"] = "Default: 0\n\
The number of iteration to expand the initial set of group elements";
  ListStringValues_doc["eCommand"] = "eCommand: lrs\n\
The serial program for computing the dual description. Possibilities: lrs, cdd";
  ListStringValues_doc["FileI"] = "The input file of the computation";
  ListStringValues_doc["FileO"] = "The output file of the computation";
  ListStringValues_doc["Arithmetic"] = "Default: rational\n\
Other possibilities are Qsqrt2, Qsqrt5 and RealAlgebraic=FileDesc where FileDesc is the description";
  SingleBlock BlockPROC;
  BlockPROC.setListBoolValues(ListBoolValues_doc);
  BlockPROC.setListIntValues(ListIntValues_doc);
  BlockPROC.setListStringValues(ListStringValues_doc);
  ListBlock["PROC"] = BlockPROC;
  // Merging all data
  return {std::move(ListBlock), "undefined"};
}

RecOption ReadInitialData(FullNamelist const &eFull) {
  SingleBlock BlockPROC = eFull.ListBlock.at("PROC");
  std::string eCommand = BlockPROC.ListStringValues.at("eCommand");
  std::string FileI = BlockPROC.ListStringValues.at("FileI");
  std::string FileO = BlockPROC.ListStringValues.at("FileO");
  std::string Arithmetic = BlockPROC.ListStringValues.at("Arithmetic");
  int n_iter_max = BlockPROC.ListIntValues.at("n_iter_max");
  int n_expand = BlockPROC.ListIntValues.at("n_expand");
  bool ComputeStabilizerPermutation = BlockPROC.ListBoolValues.at("ComputeStabilizerPermutation");
  bool ComputeGroupPresentation = BlockPROC.ListBoolValues.at("ComputeGroupPresentation");
  return {eCommand, FileI, FileO, Arithmetic, n_iter_max, n_expand, ComputeStabilizerPermutation, ComputeGroupPresentation};
}

template <typename T>
void PrintAdjacencyInfo(StepEnum<T> const &se, std::string const &FileO) {
  std::ofstream os(FileO);
  size_t n_neigh = se.ListNeighborX.size();
  os << n_neigh << "\n";
  for (size_t i_neigh = 0; i_neigh < n_neigh; i_neigh++) {
    std::pair<size_t, size_t> epair = se.ListNeighborData[i_neigh];
    PairElt<T> eElt = se.GetElement(epair);
    WriteTrackGroup(os, eElt.tg);
    os << "\n";
  }
}

void PrintGroupPresentation(std::ostream& os, std::pair<int, std::vector<std::vector<int>>> const& ThePres) {
  os << "The number of generators is " << ThePres.first << "\n";
  size_t n_word = ThePres.second.size();
  os << "The number of words is " << n_word << "\n";
  for (size_t i_word=0; i_word<n_word; i_word++) {
    os << i_word << " : {";
    std::vector<int> const& eWord = ThePres.second[i_word];
    size_t len = eWord.size();
    for (size_t u=0; u<len; u++) {
      if (u>0)
        os << ",";
      os << eWord[u];
    }
    os << "}\n";
  }
}

template <typename T,typename Tgroup> void full_process_type(RecOption const &rec_option) {
  std::cerr << "Beginning of full_process_type\n";
  DataPoincare<T> dp =
      ReadDataPoincare<T>(rec_option.FileI, rec_option.n_expand);
  std::cerr << "We have dp\n";
  StepEnum<T> se = IterativePoincareRefinement(dp, rec_option);
  std::cerr << "We have se\n";
  PrintAdjacencyInfo(se, rec_option.FileO);
  std::cerr << "se has been written to file\n";
  bool ComputeStabilizerPermutation = rec_option.ComputeStabilizerPermutation;
  bool ComputeGroupPresentation = rec_option.ComputeGroupPresentation;
  std::cerr << "ComputeStabilizerPermutation = " << ComputeStabilizerPermutation << "\n";
  std::cerr << "ComputeGroupPresentation = " << ComputeGroupPresentation << "\n";
  if (ComputeStabilizerPermutation) {
    std::cerr << "Writing the permutation group stabilizer\n";
    Tgroup GRP = se.template GetPermutationGroup<Tgroup>();
    std::cerr << "GRP=" << GRP.GapString() << "\n";
  }
  if (ComputeGroupPresentation) {
    AdjacencyInfo<T> ai = se.ComputeAdjacencyInfo(rec_option.eCommand);
    std::cerr << "Writing the group presentation\n";
    std::pair<int, std::vector<std::vector<int>>> ThePres = se.GetGroupPresentation(ai);
    PrintGroupPresentation(std::cerr, ThePres);
  }
}

template <typename Tgroup>
void Process_rec_option(RecOption const &rec_option) {
  std::string arith = rec_option.Arithmetic;
  if (arith == "rational") {
    using T = mpq_class;
    return full_process_type<T,Tgroup>(rec_option);
  }
  if (arith == "Qsqrt5") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return full_process_type<T,Tgroup>(rec_option);
  }
  if (arith == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return full_process_type<T,Tgroup>(rec_option);
  }
  std::optional<std::string> opt_realalgebraic =
      get_postfix(arith, "RealAlgebraic=");
  if (opt_realalgebraic) {
    std::string const &FileAlgebraicField = *opt_realalgebraic;
    if (!IsExistingFile(FileAlgebraicField)) {
      std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                << " is missing\n";
      throw TerminalException{1};
    }
    using T_rat = mpq_class;
    HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
    int const idx_real_algebraic_field = 1;
    insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
    using T = RealField<idx_real_algebraic_field>;
    return full_process_type<T,Tgroup>(rec_option);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

// clang-format off
#endif  // SRC_POINCARE_POLYHEDRON_POINCARE_POLYHEDRON_H_
// clang-format on

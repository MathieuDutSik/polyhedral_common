// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POINCARE_POLYHEDRON_POINCARE_POLYHEDRON_H_
#define SRC_POINCARE_POLYHEDRON_POINCARE_POLYHEDRON_H_

// clang-format off
#include "Namelist.h"
#include "POLY_DirectDualDesc.h"
#include "POLY_Fundamental.h"
#include "POLY_LinearProgramming.h"
#include "matrix_group_word.h"
#include "finite_matrix_group.h"
#include <set>
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <string>
#include <utility>
#include <vector>
// clang-format on

template <typename T> struct DataEXT {
  MyMatrix<T> EXT;
  vectface vf;
  std::vector<Face> v_red;
};

template<typename T>
DataEXT<T> DirectDataExtComputation(MyMatrix<T> const& FAC, std::string const& eCommand_DD, std::ostream& os) {
  int n_fac = FAC.rows();
  vectface vf(n_fac);
  std::vector<MyVector<T>> ListVert;
  auto f_process=[&](std::pair<Face, MyVector<T>> const& pair_face) -> void {
    vf.push_back(pair_face.first);
    ListVert.push_back(pair_face.second);
  };
  DirectFacetComputationFaceIneq(FAC, eCommand_DD, f_process, os);
  MyVector<T> EXT = MatrixFromVectorFamily(ListVert);
  int n_ext = EXT.rows();
  std::vector<Face> v_red(n_fac, Face(n_ext));
  size_t pos = 0;
  for (auto &eInc : vf) {
    for (int i_fac = 0; i_fac < n_fac; i_fac++) {
      v_red[i_fac][pos] = eInc[i_fac];
    }
    pos++;
  }
  return {std::move(EXT), std::move(vf), std::move(v_red)};
}



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
  return {std::move(EXT), std::move(v_red)};
}

//
// The elementary data structures for the Poincare stuff
//

// The initial data for the Poincare Polyhedron Theorem
// ---a point x
// ---a list of group element which is of course assumed to generate the group
template <typename T> struct DataPoincare {
  MyVector<T> x;
  std::vector<CombElt<T>> ListGroupElt;
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
  std::unordered_set<CombElt<T>> s_elt;
  s_elt.insert(GenerateIdentity<T>(n));
  for (int i_elt = 0; i_elt < n_elt; i_elt++) {
    int pos = i_elt + 1;
    MyMatrix<T> eElt = ReadMatrix<T>(is);
    MyMatrix<double> eElt_d = UniversalMatrixConversion<double, T>(eElt);
    T TheDet = DeterminantMat(eElt);
    std::cerr << "i_elt=" << i_elt << "/" << n_elt << " TheDet=" << TheDet
              << "\n";
    TrackGroup tg{{pos}};
    CombElt<T> pe{tg, eElt, eElt_d};
    s_elt.insert(pe);
  }
  std::vector<CombElt<T>> l_elt;
  for (auto &eElt : s_elt)
    l_elt.push_back(eElt);
  std::cerr << "|Initial set|=" << l_elt.size() << "\n";
  l_elt = InverseSaturation(l_elt);
  std::vector<CombElt<T>> l_gen = l_elt;
  std::cerr << "|Inverse saturation|=" << l_elt.size() << "\n";
  std::cerr << "n_expand=" << n_expand << "\n";
  for (int i_expand = 0; i_expand < n_expand; i_expand++) {
    l_elt = ListExpansion(l_elt, l_gen);
    std::cerr << "i_expand=" << i_expand << " |l_elt|=" << l_elt.size() << "\n";
  }
  return {x, l_elt};
}

template <typename T> MyMatrix<T> Contragredient(MyMatrix<T> const &M) {
  return Inverse(TransposedMat(M));
}

template <typename T>
MyMatrix<T> SmallestCanonicalization(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  MyMatrix<T> Mret(nbRow, nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<T> V = GetMatrixRow(M, iRow);
    MyVector<T> Vcan =
        CanonicalizationSmallestCoefficientVectorPlusCoeff(V).TheVect;
    AssignMatrixRow(Mret, iRow, Vcan);
  }
  return Mret;
}

template <typename T> MyMatrix<T> AddZeroColumn(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  MyMatrix<T> Mret(nbRow, nbCol + 1);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    Mret(iRow, 0) = 0;
    for (int iCol = 0; iCol < nbCol; iCol++)
      Mret(iRow, 1 + iCol) = M(iRow, iCol);
  }
  return Mret;
}

template <typename T> struct DataFAC {
  int n_mat;
  int rnk;
  MyMatrix<T> FAC;
  std::optional<MyVector<T>> eVectInt;
  std::vector<CombElt<T>> ListAdj;
  std::vector<CombElt<T>> ListAdjInv;
};

namespace boost::serialization {
  template <class Archive, typename T>
  inline void serialize(Archive &ar, DataFAC<T> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("n_mat", eRec.n_mat);
    ar &make_nvp("rnk", eRec.rnk);
    ar &make_nvp("FAC", eRec.FAC);
    ar &make_nvp("eVectInt", eRec.eVectInt);
    ar &make_nvp("ListAdjInv", eRec.ListAdjInv);
  }
}

template<typename T>
void write_description_to_file(std::string const &eFile,
                               DataFAC<T> const &datafac) {
  std::ofstream osF(eFile);
  size_t n_mat = datafac.ListAdj.size();
  osF << n_mat << "\n";
  for (size_t i_mat = 0; i_mat < n_mat; i_mat++) {
    WriteTrackGroup(osF, datafac.ListAdj[i_mat].tg);
    WriteMatrix(osF, datafac.ListAdj[i_mat].mat);
  }
  WriteMatrix(osF, datafac.FAC);
}

//
// Now the polyhedral stuff
//

struct RecOption {
  std::string method_adjacent;
  std::string eCommand_DD;
  std::string PrefixSave;
  std::string FileDataPoincare;
  std::string FileO;
  std::string Arithmetic;
  std::string Approach;
  std::string MethodMissingI;
  std::string MethodVertexMatching;
  int n_iter_max;
  int n_expand;
  bool ComputeStabilizerPermutation;
  bool ComputeGroupPresentation;
};

/*
  The program works for A3 (where actually we used B3)
  But it has some problems for the case of the cocompact group of interest to
  us.
  ---
  Problems:
  - Bad numerics where we are forced to have evaluation of number very near 0
  and the continued fraction algorithm fails us.
    => Look at the absolute values of the coefficients in double precision.
    => How hard it is to compute the symmetric function P(x_1) ... P(x_k)?
       This ought to give us a value that is rational. If P(x_1) is very small,
  likely the others are not.
    => See exactly where things go haywire.
  - Non-convergence of finding the dual description. Our hope is definitely that
  there are a lot of facets that are redundant. This is the foundation of our
  approach. Another element is that for an element, we can easily identify if it
  is new or already present. So, we have already accepted that insertion done
  step by step is the business of the day.
  - If the checks of non-redundancy are so so expensive, then we should keep
  track of them in a database of known redundant entries. This should probably
  be just a std::unordered_map<Myvector<T>,std::pair<size_t,size_t>> We would
  have a statbilizer type that is passed as a reference to the insertion
  process. and then when an insertion is done. But all that seems overkill right
  now because the stabilizer are so far trivial.
  - What should be reasonably accessible for us:
    - redundancy checks using Clarkson method (based on Linear programming)
    - Linear programming over those special fields has to be doable.
      Maybe we could use shifting to doule precision to help solve those linear
  problems.
    - dual description for a set of inequalities that are all defining facets.
  -
 */
template <typename T> struct StepEnum {
public:
  MyVector<T> x;
  std::vector<CombElt<T>> stabilizerElt;
  std::unordered_set<CombElt<T>> stabilizerElt_set;
  std::vector<CombElt<T>> ListNeighborCoset;
  std::vector<MyVector<T>> ListNeighborX;
  std::vector<std::pair<size_t, size_t>> ListNeighborData;
  std::unordered_map<MyVector<T>, std::pair<size_t, size_t>> map;
  std::unordered_set<CombElt<T>> known_redundant;
  // For a facet of the cone, there should be a matching element in the adjacent
  // facet. Otherwise return -1 and proceed from that.
  std::vector<int> ComputeMatchingVector(std::ostream& os) const {
    HumanTime time;
    int n_mat = ListNeighborData.size();
    std::vector<MyMatrix<T>> l_mat;
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      MyMatrix<T> Q = GetElement(ListNeighborData[i_mat]).mat;
      l_mat.push_back(Q);
    }
    auto get_j_mat = [&](int const &i_mat) -> int {
      MyVector<T> x_img = l_mat[i_mat].transpose() * x;
      for (int j_mat = 0; j_mat < n_mat; j_mat++) {
        MyVector<T> x2 = l_mat[j_mat].transpose() * x_img;
        if (x2 == x)
          return j_mat;
      }
      return -1;
    };
    std::vector<int> MatchVector;
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      int j_mat = get_j_mat(i_mat);
      MatchVector.push_back(j_mat);
    }
    os << "ComputingMatchingVector time=" << time << "\n";
    return MatchVector;
  }
  std::vector<int> ComputeMatchingVectorCheck(std::ostream& os) const {
    std::vector<int> V = ComputeMatchingVector(os);
    for (auto & eVal : V) {
      if (eVal == -1) {
        std::cerr << "We have a missing value which is not allowed here\n";
        throw TerminalException{1};
      }
    }
    return V;
  }
  void print_statistics(std::ostream &os_out, std::ostream& os) const {
    os_out << "|stabilizerElt|=" << stabilizerElt.size() << "\n";
    os_out << "|ListNeighborCoset|=" << ListNeighborCoset.size()
           << " |ListNeighborX|=" << ListNeighborX.size() << "\n";
    std::map<size_t, size_t> map_by_cos;
    for (auto &ePair : ListNeighborData) {
      map_by_cos[ePair.first]++;
    }
    std::map<size_t, size_t> map_by_size;
    for (auto &kv : map_by_cos) {
      map_by_size[kv.second]++;
    }
    os_out << "Sizes :";
    for (auto &kv : map_by_size) {
      os_out << " [" << kv.first << "," << kv.second << "]";
    }
    os_out << "\n";
    std::vector<int> V = ComputeMatchingVector(os);
    int n_missing = 0;
    for (size_t i_mat = 0; i_mat < V.size(); i_mat++) {
      int val = V[i_mat];
      //      os_out << "i_mat=" << i_mat << " j=" << val << "\n";
      if (val == -1) {
        n_missing++;
      }
    }
    if (n_missing > 0) {
      os_out << "ERROR: We have some matching missing n_missing=" << n_missing
         << "\n";
    } else {
      os_out << "OK: All facets have matching on the other side\n";
    }
  }
  bool IsPresentInStabilizer(CombElt<T> const &eElt) const {
    return stabilizerElt_set.count(eElt) == 1;
  }
  bool InsertStabilizerGenerator(CombElt<T> const &eElt, std::ostream & os) {
    if (IsPresentInStabilizer(eElt))
      return false;
    os << "InsertStabilizerGenerator 1 |stabilizerElt|="
       << stabilizerElt.size() << "\n";
    std::vector<CombElt<T>> ExtListGen = stabilizerElt;
    ExtListGen.push_back(eElt);
    os << "InsertStabilizerGenerator 2 |ExtListGen|="
       << ExtListGen.size() << "\n";
    stabilizerElt = GroupGeneration(ExtListGen);
    os << "InsertStabilizerGenerator 3 |stabilizerElt|="
       << stabilizerElt.size() << "\n";
    stabilizerElt_set.clear();
    for (auto &eElt : stabilizerElt)
      stabilizerElt_set.insert(eElt);
    return true;
  }
  CombElt<T> GetElement(std::pair<size_t, size_t> const &val) const {
    size_t i_coset = val.first;
    size_t i_elt = val.second;
    if (i_coset > ListNeighborCoset.size()) {
      std::cerr << "Accessing something over the index\n";
      throw TerminalException{1};
    }
    CombElt<T> prod =
        ProductComb(ListNeighborCoset[i_coset], stabilizerElt[i_elt]);
    CombElt<T> eInv = InverseComb(stabilizerElt[i_elt]);
    return ProductComb(eInv, prod);
  }
  void InsertCoset(CombElt<T> const &eCoset) {
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
      auto iter = map.find(kv.first);
      if (iter != map.end()) {
        std::cerr << "find overlap V=" << StringVector(kv.first)
                  << " pair=" << kv.second.first << "/" << kv.second.second
                  << "\n";
        std::cerr << "          iter=" << StringVector(iter->first)
                  << " pair=" << iter->second.first << "/"
                  << iter->second.second << "\n";
        throw TerminalException{1};
      }
      map[kv.first] = kv.second;
    }
    if (ListNeighborX.size() != map.size()) {
      std::cerr << "|map_local|=" << map_local.size() << "\n";
      std::cerr << "|ListNeighborX|=" << ListNeighborX.size()
                << " |map|=" << map.size() << "\n";
      throw TerminalException{1};
    }
  }
  void ComputeCosets(std::vector<CombElt<T>> const &l_elt) {
    ListNeighborX.clear();
    ListNeighborData.clear();
    ListNeighborCoset.clear();
    map.clear();
    std::vector<CombElt<T>> l_cos =
        IdentifyDoubleCosets(x, l_elt, stabilizerElt);
    for (auto &eCoset : l_cos) {
      InsertCoset(eCoset);
    }
  }
  bool IsPresentInCosetOrStabilizer(CombElt<T> const &eElt) const {
    MyVector<T> x_img = eElt.mat.transpose() * x;
    if (x_img == x) {
      return IsPresentInStabilizer(eElt);
    } else {
      auto iter = map.find(x_img);
      if (iter == map.end()) {
        return false;
      } else {
        CombElt<T> TheProd = GetElement(iter->second);
        CombElt<T> eElt_inv = InverseComb(eElt);
        CombElt<T> stab_elt = ProductComb(TheProd, eElt_inv);
        MyVector<T> x2 = stab_elt.mat.transpose() * x;
        if (x2 != x) {
          std::cerr << "x is not stabilized\n";
          throw TerminalException{1};
        }
        return IsPresentInStabilizer(stab_elt);
      }
    }
  }
  bool InsertGenerators(std::vector<CombElt<T>> const &ListGen, std::ostream& os) {
    bool DidSomething = false;
    auto generator_upgrade = [&](CombElt<T> const &e_elt) -> void {
      bool test = InsertStabilizerGenerator(e_elt, os);
      if (test) {
        // Copy needed of the old data then recompute
        std::vector<CombElt<T>> OldListCos = ListNeighborCoset;
        ComputeCosets(OldListCos);
        DidSomething = true;
      }
    };
    auto f_insert = [&](CombElt<T> const &e_elt) -> void {
      if (known_redundant.count(e_elt) == 1) {
        os << "Exiting f_insert because the element is already known to be redundant\n";
        return;
      }
      MyVector<T> x_img = e_elt.mat.transpose() * x;
      if (x_img == x) {
        generator_upgrade(e_elt);
      } else {
        auto iter = map.find(x_img);
        if (iter == map.end()) {
          InsertCoset(e_elt);
          DidSomething = true;
        } else {
          CombElt<T> TheProd = GetElement(iter->second);
          CombElt<T> e_elt_inv = InverseComb(e_elt);
          CombElt<T> stab_elt = ProductComb(TheProd, e_elt_inv);
          MyVector<T> x2 = stab_elt.mat.transpose() * x;
          if (x2 != x) {
            std::cerr << "x is not stabilized\n";
            throw TerminalException{1};
          }
          generator_upgrade(stab_elt);
        }
      }
    };
    for (auto &e_elt : ListGen)
      f_insert(e_elt);
    print_statistics(os, os);
    return DidSomething;
  }
  void initial_set(MyVector<T> const &_x) {
    x = _x;
    int n = x.size();
    CombElt<T> IdMat = GenerateIdentity<T>(n);
    stabilizerElt = {IdMat};
    stabilizerElt_set.insert(IdMat);
  }
  MyVector<T> GetIneq(CombElt<T> const& e_elt) const {
    MyMatrix<T> const &eMat = e_elt.mat;
    MyVector<T> x_img = eMat.transpose() * x;
    MyVector<T> x_diff = x_img - x;
    return x_diff;
  }
  MyMatrix<T> GetFAC(std::ostream& os) const {
    int n = x.size();
    int n_mat = ListNeighborX.size();
    MyMatrix<T> FAC(n_mat, n);
    std::unordered_map<MyVector<T>, size_t> map_test;
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      MyVector<T> x_diff = ListNeighborX[i_mat] - x;
      size_t &val = map_test[x_diff];
      if (val != 0) {
        size_t j_mat = val - 1;
        os << "Collision found between i_mat=" << i_mat
           << " j_mat=" << j_mat << "\n";
        std::pair<size_t, size_t> i_pair = ListNeighborData[i_mat];
        std::pair<size_t, size_t> j_pair = ListNeighborData[j_mat];
        std::cerr << "i_pair=" << i_pair.first << "/" << i_pair.second
                  << " j_pair=" << j_pair.first << "/" << j_pair.second << "\n";
        CombElt<T> ePair = GetElement(i_pair);
        CombElt<T> fPair = GetElement(j_pair);
        MyVector<T> eV = ePair.mat.transpose() * x;
        MyVector<T> fV = fPair.mat.transpose() * x;
        CombElt<T> ePairInv = InverseComb(ePair);
        CombElt<T> eStabElt = ProductComb(fPair, ePairInv);
        MyVector<T> xImg = eStabElt.mat.transpose() * x;
        std::cerr << "eV=" << StringVector(eV) << " fV=" << StringVector(fV)
                  << "\n";
        std::cerr << "x=" << StringVector(x) << " xImg=" << StringVector(xImg)
                  << "\n";
        bool test = IsPresentInStabilizer(eStabElt);
        std::cerr << "test=" << test << "\n";
        throw TerminalException{1};
      }
      val = i_mat + 1;
      AssignMatrixRow(FAC, i_mat, x_diff);
    }
    return FAC;
  }
  DataFAC<T> GetDataCone(std::ostream& os) const {
    HumanTime time;
    MyMatrix<T> FAC = GetFAC(os);
    int n_mat = FAC.rows();
    int rnk = RankMat(FAC);
    std::optional<MyVector<T>> eVectInt;
    if (rnk == FAC.cols()) {
      eVectInt = GetSpaceInteriorPoint_Basic(FAC, os);
    }
    std::vector<CombElt<T>> ListAdj;
    std::vector<CombElt<T>> ListAdjInv;
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      CombElt<T> uElt = GetElement(ListNeighborData[i_mat]);
      CombElt<T> uEltInv = InverseComb(uElt);
      ListAdj.push_back(uElt);
      ListAdjInv.push_back(uEltInv);
    }
    os << "New DataFAC computed time=" << time << "\n";
    return {n_mat, rnk, FAC, eVectInt, ListAdj, ListAdjInv};
  }
  int RemoveRedundancy(std::ostream& os) {
    HumanTime time;
    MyMatrix<T> FAC = GetFAC(os);
    int n = x.size();
    int n_mat = FAC.rows();
    int rnk = RankMat(FAC);
    os << "RemoveRedundancy : n=" << n << " n_mat=" << n_mat << " rnk=" << rnk << "\n";
    if (rnk != n) {
      std::cerr << "Error in RemoveRedundancy\n";
      std::cerr << "n=" << n << " n_mat=" << n_mat << " rnk=" << rnk << "\n";
      throw TerminalException{1};
    }
    MyMatrix<T> FACexp = AddZeroColumn(FAC);
    //
    // Doing the redundancy computation
    //
    os << "Before RedundancyReductionClarkson n_mat=" << n_mat << "\n";
    std::vector<int> ListIrred = cdd::RedundancyReductionClarkson(FACexp, os);
    os << "|ListIrred|=" << ListIrred.size() << " n_mat=" << n_mat << " time=" << time << "\n";
    //
    // Paperwork
    //
    Face f_status_keep(n_mat);
    std::set<size_t> l_keep;
    for (auto &i_mat : ListIrred) {
      std::pair<size_t, size_t> epair = ListNeighborData[i_mat];
      size_t i_coset = epair.first;
      l_keep.insert(i_coset);
      f_status_keep[i_mat] = 1;
    }
    //
    // Updating the list of known redundants
    //
    int n_remove = 0;
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      if (f_status_keep[i_mat] == 0) {
        CombElt<T> uElt = GetElement(ListNeighborData[i_mat]);
        known_redundant.insert(uElt);
        n_remove++;
      }
    }
    //
    // Check
    //
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      std::pair<size_t, size_t> epair = ListNeighborData[i_mat];
      size_t i_coset = epair.first;
      if (l_keep.count(i_coset) != f_status_keep[i_mat]) {
        std::cerr << "There is incoherency in the orbit nature\n";
        throw TerminalException{1};
      }
    }
    os << "RemoveRedundancy : |l_keep|=" << l_keep.size()
       << " n_remove=" << n_remove << " time=" << time << "\n";
    if (n_remove > 0) {
      std::vector<CombElt<T>> ListNeighborCosetRed;
      for (auto &i_coset : l_keep) {
        ListNeighborCosetRed.push_back(ListNeighborCoset[i_coset]);
      }
      ComputeCosets(ListNeighborCosetRed);
      print_statistics(os, os);
    }
    os << "ComputeCosets time=" << time << "\n";
    return n_remove;
  }
  template <typename Tgroup>
  Tgroup GetPermutationGroup() const {
    using Telt = typename Tgroup::Telt;
    using Tidx = typename Telt::Tidx;
    int n_neigh = ListNeighborX.size();
    std::unordered_map<MyVector<T>, int> map_rev;
    for (int i_neigh = 0; i_neigh < n_neigh; i_neigh++)
      map_rev[ListNeighborX[i_neigh]] = i_neigh + 1;
    auto GetPermutation = [&](CombElt<T> const &eElt) -> Telt {
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
    auto f_insert = [&](CombElt<T> const &eElt) -> void {
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
};

namespace boost::serialization {
  template <class Archive, typename T>
  inline void serialize(Archive &ar, StepEnum<T> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("x", eRec.x);
    ar &make_nvp("stabilizerElt", eRec.stabilizerElt);
    ar &make_nvp("stabilizerElt_set", eRec.stabilizerElt_set);
    ar &make_nvp("ListNeighborCoset", eRec.ListNeighborCoset);
    ar &make_nvp("ListNeighborX", eRec.ListNeighborX);
    ar &make_nvp("ListNeighborData", eRec.ListNeighborData);
    ar &make_nvp("map", eRec.map);
    ar &make_nvp("known_redundant", eRec.known_redundant);
  }
}

// If there is an element g defining a facet then there should be
// another matching element on the other side. That element has
// to be g^{-1} if the stabilizer is trivial. But that is not
// necessarily the case. And of course g^{-1} = g in the Coxeter
// case.
//
// In the process of building the set of group elements that
// define inequalities, we often find ourself in the situation
// where an element g defines an irredundant inequality but g^{-1}
// defines a redundant inequality. That means that g is actually
// made redundant by elements that are yet to be discovered.
//
// The idea is then to find an interior point to the facet
// and from that we are using the SGE to find the missing
// inequality.
//
// Using the inverse and other strategies do not seem to work well.
//
// BELOW is an attempt to find the missing inequality by solving
// with non-negative coefficients. It failed, but here goes:
// Inequalities is
//   x.a <= x.aw
// Inequality is f_w(x) = x.aw - x.a
// If there are redundancy for w^(-1) then we can write
// f_{w^(-1)}(x) = a_1 f_w1(x) + .... + a_N f_wN(x)   with a_j >= 0
// The transformation rule that we have is
// x.a = xc(w).aw  with c(w) = (w^(-1))^T
// Expanding this gets us
// x.aw^(-1) - x.a = sum a_j (x.aw_i - x.a)
// xc(w).a - xc(w).aw = sum a_j (xc(w).aw_iw - xc(w).aw)
// So maybe inserting the w_i w is a good idea.
template<typename T>
std::vector<CombElt<T>>
GetMissingInverseElement(StepEnum<T> const& se,
                         DataFAC<T> const &datafac,
                         ShortVectorGroup<T> const &svg, std::ostream& os) {
  HumanTime time;
  std::vector<int> V = se.ComputeMatchingVector(os);
  ShortVectorGroupMemoize<T> svg_mem(svg);
  int n_mat = se.ListNeighborData.size();
  std::vector<CombElt<T>> ListMiss;
  bool attempt_inverse = false;
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    if (V[i_mat] == -1) {
      auto get_relevant_inverse=[&]() -> std::optional<CombElt<T>> {
        if (!attempt_inverse) {
          return {};
        }
        CombElt<T> w = se.GetElement(se.ListNeighborData[i_mat]);
        CombElt<T> wInv = InverseComb(w);
        MyVector<T> x_ineq = se.GetIneq(wInv);
        HumanTime time1;
        std::optional<MyVector<T>> opt =
          SolutionMatNonnegative(datafac.FAC, x_ineq, os);
        os << "|SolutionMatNonnegative|=" << time1 << "\n";
        if (opt) {
          return {};
        } else {
          return wInv;
        }
      };
      std::optional<CombElt<T>> opt = get_relevant_inverse();
      if (opt) {
        ListMiss.push_back(*opt);
        os << "wInv actually define a new inequality\n";
      } else {
        // Finding by nearest group point.
        Face f(n_mat);
        f[i_mat] = 1;
        HumanTime time2;
        MyVector<T> eVectInt = GetSpaceInteriorPointFace(datafac.FAC, f, os);
        os << "|GetSpaceInteriorPointFace|=" << time2 << "\n";
        T target_scal = eVectInt.dot(se.x);
        svg_mem.ComputeInsertSolution(eVectInt, target_scal);
        os << "Found new elements by Short Group Element\n";
      }
    }
  }
  for (auto &eElt : svg_mem.GetListMiss()) {
    ListMiss.push_back(eElt);
    ListMiss.push_back(InverseComb(eElt));
  }
  os << "GetMissingInverseElement time=" << time << "\n";
  return ListMiss;
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
};

struct Tfacet {
  std::vector<TsingAdj> l_sing_adj;
};

template <typename T> struct AdjacencyInfo {
  std::vector<Tfacet> ll_adj;
};

// For each element g defining a facet, we have a corresponding
// element that defines the corresponding facet. That element can
// be g^{-1} or it can be something else.
//
// But it the facet are actually matching, that does not mean
// that all is fine. The facet could have the same defining inequality
// but they could not be geometrically the same.
//
// If the facet are not matching then we can find offending vertices
// by linear programming and then with SGE we can find the corresponding
// missed element.
template<typename T>
std::vector<CombElt<T>>
GetMissingFacetMatchingElement_LP(StepEnum<T> const& se,
                                  DataFAC<T> const &datafac,
                                  ShortVectorGroup<T> const &svg,
                                  std::ostream& os) {
  HumanTime time;
  //
  // Preprocessing information
  //
  int dim = datafac.FAC.cols();
  int n_fac = datafac.FAC.rows();
  os << "GetMissingFacetMatchingElement_LP, n_fac=" << n_fac << " dim=" << dim << "\n";
  ShortVectorGroupMemoize<T> svg_mem(svg);
  Face f_adj = ComputeSkeletonClarkson(datafac.FAC, os);
  os << "We have f_adj, n_fac=" << n_fac << " time=" << time << "\n";
  //
  auto get_nsp = [&](int const &i_fac) -> MyMatrix<T> {
    MyMatrix<T> Equa(dim, 1);
    for (int i = 0; i < dim; i++)
      Equa(i, 0) = datafac.FAC(i_fac, i);
    return NullspaceMat(Equa);
  };
  auto get_localfac = [&](int const &i_fac, MyMatrix<T> const &TheFAC,
                          MyMatrix<T> const &NSP) -> MyMatrix<T> {
    int count = 0;
    for (int j_fac = 0; j_fac < n_fac; j_fac++) {
      if (f_adj[i_fac + j_fac * n_fac] == 1) {
        count++;
      }
    }
    MyMatrix<T> FAC_local(count, dim - 1);
    int pos = 0;
    for (int k_fac = 0; k_fac < n_fac; k_fac++) {
      if (f_adj[i_fac + k_fac * n_fac] == 1) {
        MyVector<T> eFAC = GetMatrixRow(TheFAC, k_fac);
        MyVector<T> eFACred = NSP * eFAC;
        AssignMatrixRow(FAC_local, pos, eFACred);
        pos++;
      }
    }
    return FAC_local;
  };
  std::vector<int> V = se.ComputeMatchingVectorCheck(os);
  for (int i_fac = 0; i_fac < n_fac; i_fac++) {
    MyMatrix<T> const &Q = datafac.ListAdj[i_fac].mat;
    int j_fac = V[i_fac];
    MyMatrix<T> NSP = get_nsp(i_fac);
    MyMatrix<T> FACimg = datafac.FAC * Q;
    MyMatrix<T> FAC_local = get_localfac(i_fac, datafac.FAC, NSP);
    int n_adj = FAC_local.rows();
    MyMatrix<T> FAC_localImg = get_localfac(j_fac, FACimg, NSP);
    int n_adj_img = FAC_localImg.rows();
    int n_found = 0;
    for (int i_adj_img = 0; i_adj_img < n_adj_img; i_adj_img++) {
      MyVector<T> eVect = GetMatrixRow(FAC_localImg, i_adj_img);
      SolutionMatNonnegativeComplete<T> SolCompl =
        GetSolutionMatNonnegativeComplete(FAC_local, eVect, os);
      if (!SolCompl.SolNonnegative) {
        if (!SolCompl.ExtremeRay) {
          std::cerr << "Failed to have an extreme ray\n";
          throw TerminalException{1};
        }
        MyVector<T> const &eEXTred = *SolCompl.ExtremeRay;
        MyVector<T> eEXT = NSP.transpose() * eEXTred;
        T target_scal = eEXT.dot(se.x);
        svg_mem.ComputeInsertSolution(eEXT, target_scal);
        n_found++;
      }
    }
    os << "i_fac=" << i_fac << "/" << n_fac << "  n_adj=" << n_adj
       << " n_adj_img=" << n_adj_img << " n_found=" << n_found << "\n";
  }
  //
  // Now calling the SGE code
  //
  os << "We have computed svg_mem, time=" << time << "\n";
  std::vector<CombElt<T>> ListMiss;
  for (auto &eElt : svg_mem.GetListMiss()) {
    ListMiss.push_back(eElt);
    ListMiss.push_back(InverseComb(eElt));
  }
  os << "Returning |ListMiss|=" << ListMiss.size() << " time=" << time << "\n";
  return ListMiss;
}

// The same as above but using a dual description.
template<typename T>
std::vector<CombElt<T>>
GetMissingFacetMatchingElement_DD(StepEnum<T> const& se,
                                  DataFAC<T> const &datafac,
                                  std::string const &eCommand_DD,
                                  ShortVectorGroup<T> const &svg,
                                  std::ostream& os) {
  HumanTime time;
  //
  // Preprocessing information
  //
  os << "GetMissingFacetMatchingElement_DD, beginning\n";
  DataEXT<T> dataext = DirectDataExtComputation(datafac.FAC, eCommand_DD, os);
  os << "|EXT|=" << dataext.EXT.rows() << " / " << dataext.EXT.cols() << " time=" << time << "\n";
  (void)se.ComputeMatchingVectorCheck(os);
  //
  // Building incidence informations
  //
  int n_fac = datafac.FAC.rows();
  int n_ext = dataext.EXT.rows();
  int n = dataext.EXT.cols();
  os << "n_fac=" << n_fac << " n_ext=" << n_ext << " n=" << n << "\n";
  Face f_adjacency(n_fac * n_fac);
  for (int i_fac = 0; i_fac < n_fac; i_fac++) {
    for (int j_fac = i_fac + 1; j_fac < n_fac; j_fac++) {
      Face f(n_ext);
      Face f1 = dataext.vf[i_fac];
      Face f2 = dataext.vf[j_fac];
      for (int i_ext = 0; i_ext < n_ext; i_ext++) {
        if (f1[i_ext] && f2[i_ext])
          f[i_ext] = 1;
      }
      if (static_cast<int>(f.count()) >= n - 2) {
        MyMatrix<T> EXT_red = SelectRow(dataext.EXT, f);
        int rnk = RankMat(EXT_red);
        if (rnk == n - 2) {
          f_adjacency[i_fac + n_fac * j_fac] = 1;
          f_adjacency[j_fac + n_fac * i_fac] = 1;
        }
      }
    }
  }
  for (int i_fac = 0; i_fac < n_fac; i_fac++) {
    int n_adj = 0;
    for (int j_fac = 0; j_fac < n_fac; j_fac++) {
      if (f_adjacency[i_fac + n_fac * j_fac] == 1) {
        n_adj++;
      }
    }
    os << "i_fac=" << i_fac << " n_adj=" << n_adj << "\n";
  }
  os << "We have f_adjacency, time=" << time << "\n";
  //
  // Determining the vertices which are
  //
  Face f_insert_svg(n_ext);
  for (int i_fac = 0; i_fac < n_fac; i_fac++) {
    MyMatrix<T> const &Q = datafac.ListAdj[i_fac].mat;
    MyMatrix<T> cQ = Contragredient(Q);
    MyMatrix<T> EXTimg = dataext.EXT * cQ;
    MyMatrix<T> FACimg = datafac.FAC * Q;
    Face f = dataext.vf[i_fac];
    auto f_set = [&](int i_ext) {
      if (f[i_ext] == 0) {
        // Not in face, so no insertion to do
        return;
      }
      if (f_insert_svg[i_ext] == 1) {
        // Already selected, so no need to compute
        return;
      }
      for (int i_fac = 0; i_fac < n_fac; i_fac++) {
        T scal = 0;
        for (int i = 0; i < n; i++) {
          scal += dataext.EXT(i_ext, i) * FACimg(i_fac, i);
        }
        if (scal < 0) {
          f_insert_svg[i_ext] = 1;
          return;
        }
      }
    };
    for (int i_ext = 0; i_ext < n_ext; i_ext++)
      f_set(i_ext);
  }
  size_t count = f_insert_svg.count();
  os << "We have f_insert_svg |f_insert_svg|=" << count
     << " n_ext=" << n_ext << " time=" << time << "\n";
  //
  // Now calling the SGE code
  //
  ShortVectorGroupMemoize<T> svg_mem(svg);
  size_t pos = 0;
  for (int i_ext = 0; i_ext < n_ext; i_ext++) {
    if (f_insert_svg[i_ext] == 1) {
      os << "pos=" << pos << " / " << count << "      i_ext=" << i_ext << " / " << n_ext << "\n";
      MyVector<T> eEXT = GetMatrixRow(dataext.EXT, i_ext);
      T target_scal = eEXT.dot(se.x);
      svg_mem.ComputeInsertSolution(eEXT, target_scal);
      pos++;
    }
  }
  os << "We have computed svg_mem, time=" << time << "\n";
  std::vector<CombElt<T>> ListMiss;
  for (auto &eElt : svg_mem.GetListMiss()) {
    ListMiss.push_back(eElt);
    ListMiss.push_back(InverseComb(eElt));
  }
  os << "Returning |ListMiss|=" << ListMiss.size() << " time=" << time << "\n";
  return ListMiss;
}

template<typename T>
std::vector<CombElt<T>>
GetMissingFacetMatchingElement(StepEnum<T> const& se,
                               DataFAC<T> const &datafac,
                               std::string const& method_adjacent,
                               std::string const &eCommand_DD,
                               ShortVectorGroup<T> const &svg,
                               std::ostream& os) {
  if (method_adjacent == "linear_programming") {
    return GetMissingFacetMatchingElement_LP(se, datafac, svg, os);
  }
  if (method_adjacent == "dual_description") {
    return GetMissingFacetMatchingElement_DD(se, datafac, eCommand_DD, svg, os);
  }
  std::cerr << "Failed to find a matching for\n";
  std::cerr << "method_adjacent=" << method_adjacent << "\n";
  throw TerminalException{1};
}

// The ComputeAdjacencyInfo gives the corresponding
// matching between facets and the adjacent domain.
// That computation makes sense only if the
// GetMissingFacetMatchingElement functions return
// no element.
template<typename T>
AdjacencyInfo<T> ComputeAdjacencyInfo(StepEnum<T> & se,
                                      std::string const &eCommand_DD, std::ostream& os) {
  size_t miss_val = std::numeric_limits<size_t>::max();
  MyMatrix<T> FAC = se.GetFAC(os);
  int n = se.x.size();
  int n_mat = se.ListNeighborData.size();
  int rnk = RankMat(FAC);
  if (rnk != n) {
    std::cerr << "Error in ComputeAdjacencyInfo\n";
    std::cerr << "n=" << n << " n_mat=" << n_mat << " rnk=" << rnk << "\n";
    throw TerminalException{1};
  }
  os << "ComputeAdjacencyInfo FAC.rows=" << FAC.rows()
     << " FAC.cols=" << FAC.cols() << " n_mat=" << n_mat << " n=" << n
     << " nk=" << rnk << "\n";
  DataEXT<T> dataext = DirectDataExtComputation(FAC, eCommand_DD, os);
  int n_ext = dataext.EXT.rows();
  os << "n_ext=" << n_ext << "\n";
  os << "First part: adjacency structure within the polyhedron\n";
  std::vector<Tfacet> ll_adj;
  std::vector<std::vector<Face>> ll_ridges;
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    Face const &f1 = dataext.v_red[i_mat];
    std::vector<TsingAdj> l_adj;
    std::vector<Face> l_ridges;
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
          l_adj.push_back({j_mat_s, miss_val, miss_val, miss_val});
          l_ridges.push_back(f);
        }
      }
    }
    ll_adj.push_back({l_adj});
    ll_ridges.push_back(l_ridges);
  }
  auto get_iPoly = [&](size_t iFace, Face const &f1) -> size_t {
    size_t n_adjB = ll_adj[iFace].l_sing_adj.size();
    for (size_t i_adjB = 0; i_adjB < n_adjB; i_adjB++) {
      Face f2 = ll_ridges[iFace][i_adjB];
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
      Face f1 = ll_ridges[i_mat][i_adj];
      ll_adj[i_mat].l_sing_adj[i_adj].iPolyAdj = get_iPoly(iFaceAdj, f1);
    }
  }
  os << "Second part: computing the opposite facets\n";
  MyMatrix<T> EXTcan = SmallestCanonicalization(dataext.EXT);
  std::vector<int> V = se.ComputeMatchingVectorCheck(os);
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    int j_mat = V[i_mat];
    size_t j_mat_s = j_mat;
    size_t n_adj = ll_adj[i_mat].l_sing_adj.size();
    MyMatrix<T> Q = se.GetElement(se.ListNeighborData[i_mat]).mat;
    MyMatrix<T> cQ = Contragredient(Q);
    MyMatrix<T> EXTimg = dataext.EXT * cQ;
    MyMatrix<T> EXTimgCan = SmallestCanonicalization(EXTimg);
    ContainerMatrix<T> Cont(EXTimgCan);
    std::map<int, size_t> map_index;
    for (int i_ext = 0; i_ext < n_ext; i_ext++) {
      if (dataext.v_red[i_mat][i_ext] == 1) {
        std::optional<size_t> opt =
          Cont.GetIdx_f([&](int i) -> T { return EXTcan(i_ext, i); });
        if (opt) {
          size_t idx = *opt;
          map_index[i_ext] = idx;
        }
      }
    }
    size_t iFaceOpp = j_mat_s;
    for (size_t i_adj = 0; i_adj < n_adj; i_adj++) {
      Face f_map(n_ext);
      Face const& f = ll_ridges[i_mat][i_adj];
      for (int i_ext = 0; i_ext < n_ext; i_ext++) {
        if (f[i_ext] == 1) {
          size_t idx = map_index[i_ext];
          f_map[idx] = 1;
        }
      }
      size_t iPolyOpp = get_iPoly(iFaceOpp, f_map);
      ll_adj[i_mat].l_sing_adj[i_adj].iFaceOpp = iFaceOpp;
      ll_adj[i_mat].l_sing_adj[i_adj].iPolyOpp = iPolyOpp;
    }
  }
  return {ll_adj};
}

//
// Now the code for advancing the optimization procedure.
//
template<typename T>
std::optional<CombElt<T>> GetMissing_TypeI_Gen2(StepEnum<T> const& se,
                                                DataFAC<T> const &datafac,
                                                CombElt<T> const &TestElt,
                                                int const &max_iter,
                                                std::ostream& os) {
  struct ResultOptim {
    T the_scal;
    MyVector<T> the_x;
    std::vector<int> eList;
  };
  MyVector<T> eVectInt = unfold_opt(datafac.eVectInt, "eVectInt not assigned");
  int n_mat = datafac.FAC.rows();
  auto f_start = [&](MyVector<T> const &start_x) -> ResultOptim {
    T scal = start_x.dot(eVectInt);
    return {scal, start_x, {}};
  };
  auto f_get_ineq = [&](ResultOptim const &ro) -> MyVector<T> {
    MyVector<T> x_ret = ro.the_x - se.x;
    return x_ret;
  };
  auto f_get_combelt = [&](int i_mat) -> CombElt<T> {
    if (i_mat > 0) {
      int pos = i_mat - 1;
      return datafac.ListAdj[pos];
    } else {
      int pos = -i_mat - 1;
      return datafac.ListAdjInv[pos];
    }
  };
  auto f_increment = [&](ResultOptim const &ro, int i_mat) -> ResultOptim {
    CombElt<T> eAdj = f_get_combelt(i_mat);
    MyVector<T> test_x = eAdj.mat.transpose() * ro.the_x;
    T test_scal = test_x.dot(eVectInt);
    std::vector<int> eList;
    eList.push_back(i_mat);
    return {std::move(test_scal), std::move(test_x), std::move(eList)};
  };
  auto f_all_decrease =
    [&](ResultOptim const &ro) -> std::vector<ResultOptim> {
    std::unordered_set<MyVector<T>> l_done;
    std::vector<ResultOptim> l_active{ro};
    std::vector<ResultOptim> l_total;
    int iter = 0;
    while (true) {
      std::vector<ResultOptim> l_result;
      os << "iter=" << iter << " f_all_decrease |l_active|=" << l_active.size() << "\n";
      // It is uncertain whether looking at negative values is a good idea at
      // all. What is certain is that it ends the infinite loop.
      bool has_negative = false;
      for (auto &ro : l_active) {
        for (int i = 0; i < datafac.n_mat; i++) {
          ResultOptim ro_new = f_increment(ro, i + 1);
          if (ro_new.the_scal < ro.the_scal) {
            if (l_done.count(ro_new.the_x) == 0) {
              l_result.push_back(ro_new);
              l_total.push_back(ro_new);
              l_done.insert(ro_new.the_x);
              if (ro_new.the_scal < 0) {
                has_negative = true;
              }
            }
          }
        }
      }
      os << "   f_all_decrease |l_result|=" << l_result.size() << " has_negative=" << has_negative << "\n";
      if (l_result.size() == 0 || has_negative) {
        return l_total;
      }
      l_active = std::move(l_result);
      iter++;
    }
  };
  auto f_decrease_best = [&](ResultOptim const &ro) -> ResultOptim {
    ResultOptim ro_ret = ro;
    std::vector<ResultOptim> l_result = f_all_decrease(ro);
    os << "f_decrease_best: |l_result|=" << l_result.size() << "\n";
    for (auto &ro_cand : l_result) {
      if (ro_cand.the_scal < ro_ret.the_scal) {
        ro_ret = ro_cand;
      }
    }
    return ro_ret;
  };
  auto f_evaluate = [&](ResultOptim const &ro) -> CombElt<T> {
    CombElt<T> RetElt = TestElt;
    for (auto &i_mat : ro.eList) {
      RetElt = ProductComb(RetElt, f_get_combelt(i_mat));
    }
    return RetElt;
  };
  std::unordered_set<MyVector<T>> l_new_point;
  struct ResultLP {
    // 0 for finding new inequality
    // 1 for finding the zero vector
    // 2 for managed to find a new position
    // 3 Failed to find a new position by linear programming
    int result;
    ResultOptim ro_new;
  };
  auto f_linear_programming = [&](ResultOptim const &ro) -> ResultLP {
    MyVector<T> x_ineq = f_get_ineq(ro);
    if (IsZeroVector(x_ineq)) {
      return {1, {}};
    }
    std::optional<MyVector<T>> opt =
      SolutionMatNonnegative(datafac.FAC, x_ineq, os);
    if (opt) {
      MyVector<T> V = *opt;
      std::vector<ResultOptim> l_cand;
      for (int i = 0; i < n_mat; i++) {
        if (V(i) > 0) {
          ResultOptim ro_new = f_increment(ro, -1 - i);
          if (l_new_point.count(ro_new.the_x) == 0) {
            l_cand.push_back(ro_new);
          }
        }
      }
      size_t len = l_cand.size();
      if (len == 0) {
        os << "Case 2, failed to find a new position\n";
        return {3, {}};
      }
      size_t pos = rand() % len;
      return {2, l_cand[pos]};
    } else {
      os << "  SolMatNonNeg : no solution found\n";
      return {0, {}};
    }
  };
  auto f_random_move = [&](ResultOptim const &ro) -> ResultOptim {
    ResultOptim ro_new = ro;
    while (true) {
      int sign = 2 * (rand() % 2) - 1;
      size_t i_mat = rand() % n_mat;
      int pos = sign * (i_mat + 1);
      ro_new = f_increment(ro_new, pos);
      if (l_new_point.count(ro_new.the_x) == 0) {
        return ro_new;
      }
    }
  };
  MyVector<T> start_x = TestElt.mat.transpose() * se.x;
  ResultOptim ro = f_start(start_x);
  int n_iter = 0;
  while (true) {
    os << "  n_iter=" << n_iter << "\n";
    l_new_point.insert(ro.the_x);
    ro = f_decrease_best(ro);
    os << "    We have ro\n";
    ResultLP res = f_linear_programming(ro);
    os << "    res.result=" << res.result << "\n";
    if (res.result == 0) {
      CombElt<T> ResidualElt = f_evaluate(ro);
      os << "  SolMatNonNeg : no solution found\n";
      return ResidualElt;
    }
    if (res.result == 1) {
      CombElt<T> ResidualElt = f_evaluate(ro);
      bool test_stab = se.IsPresentInStabilizer(ResidualElt);
      if (test_stab) {
        os << "It is already in stabilizer so nothing to be done\n";
        return {};
      } else {
        os << "It stabilizes but is not already known so interesting by itself\n";
        return ResidualElt;
      }
    }
    if (res.result == 2) {
      ro = res.ro_new;
    }
    if (res.result == 3) {
      ro = f_random_move(ro);
    }
    n_iter++;
    if (max_iter > 0) {
      os << "    Going to iteration termination check\n";
      if (max_iter < n_iter) {
        os << "    Exiting the enumeration\n";
        return {};
      }
    }
  }
}

template<typename T>
std::optional<CombElt<T>> GetMissing_TypeI_Gen1(StepEnum<T> const& se,
                                                DataFAC<T> const &datafac,
                                                CombElt<T> const &TestElt,
                                                int const &max_iter,
                                                std::ostream& os) {
  std::string strategy = "strategy2";
  CombElt<T> WorkElt = TestElt;
  int n_iter = 0;
  os << "Beginning of f_insert\n";
  MyVector<T> eVectInt = unfold_opt(datafac.eVectInt, "eVectInt not assigned");
  MyVector<T> curr_x = WorkElt.mat.transpose() * se.x;
  T curr_scal = curr_x.dot(eVectInt);
  T target_scal = se.x.dot(eVectInt);
  double target_scal_d = UniversalScalarConversion<double, T>(target_scal);
  os << "target_scal=" << target_scal_d << "\n";
  while (true) {
    os << "  n_iter=" << n_iter << "\n";
    MyVector<T> x_ineq = se.GetIneq(WorkElt);
    if (IsZeroVector(x_ineq)) {
      bool test_stab = se.IsPresentInStabilizer(WorkElt);
      if (test_stab) {
        os << "It is already in stabilizer so nothing to be done\n";
        return {};
      } else {
        os << "It stabilizes but is not already known so interesting by itself\n";
          return WorkElt;
      }
    }
    if (strategy == "strategy1") {
      std::optional<MyVector<T>> opt =
        SolutionMatNonnegative(datafac.FAC, x_ineq, os);
      if (opt) {
        MyVector<T> V = *opt;
        // If we take just 1 then we go into infinite loops.
        for (int u = 0; u < V.size(); u++) {
          if (V(u) > 0) {
            CombElt<T> uElt = se.GetElement(se.ListNeighborData[u]);
            CombElt<T> uEltInv = InverseComb(uElt);
            WorkElt = ProductComb(WorkElt, uEltInv);
          }
        }
      } else {
        os << "  SolMatNonNeg : no solution found\n";
        return WorkElt;
      }
    }
    if (strategy == "strategy2") {
      int n_oper = 0;
      if (curr_scal < target_scal) {
        os << "The decreasing process has found some contradiction\n";
        return WorkElt;
      }
      for (auto &eAdj : datafac.ListAdj) {
        MyVector<T> test_x = eAdj.mat.transpose() * curr_x;
        T test_scal = test_x.dot(eVectInt);
        if (test_scal < curr_scal) {
          double curr_scal_d =
            UniversalScalarConversion<double, T>(curr_scal);
          double test_scal_d =
            UniversalScalarConversion<double, T>(test_scal);
          os << "    Improving scalar product from curr_scal_d="
             << curr_scal_d << " to test_scal_d=" << test_scal_d
             << "\n";
          curr_x = test_x;
          curr_scal = test_scal;
          WorkElt = ProductComb(WorkElt, eAdj);
          n_oper++;
          if (curr_scal < target_scal) {
            os << "    The decreasing process has found some contradiction\n";
            return WorkElt;
          }
        }
      }
      if (n_oper == 0) {
        os << "  Switching to strategy1\n";
        strategy = "strategy1";
      }
    }
    n_iter++;
    if (max_iter > 0) {
      os << "    Going to iteration termination check\n";
      if (max_iter < n_iter) {
        os << "    Exiting the enumeration\n";
        return {};
      }
    }
  }
}

template<typename T>
std::optional<CombElt<T>>
GetMissing_TypeI(StepEnum<T> const& se,
                 DataFAC<T> const &datafac, CombElt<T> const &TestElt,
                 int const &max_iter,
                 std::string const &MethodMissingI, std::ostream& os) {
  if (MethodMissingI == "Gen1") {
    return GetMissing_TypeI_Gen1(se, datafac, TestElt, max_iter, os);
  }
  if (MethodMissingI == "Gen2") {
    return GetMissing_TypeI_Gen2(se, datafac, TestElt, max_iter, os);
  }
  std::cerr << "Failed to find a matching entry. MethodMissingI="
            << MethodMissingI << "\n";
  throw TerminalException{1};
}


template<typename T>
bool TestIntersection(MyMatrix<T> const &FAC, CombElt<T> const &eElt, std::ostream& os) {
  MyMatrix<T> FACimg = FAC * eElt.mat;
  MyMatrix<T> FACtot = Concatenate(FAC, FACimg);
  return IsFullDimensional(FACtot, os);
}

// Generate the missing neighbors that were discovered
// by the ridges
template<typename T>
std::vector<CombElt<T>> GenerateTypeIIneighbors(StepEnum<T> const& se,
                                                AdjacencyInfo<T> const &ai, std::ostream& os) {
  int n = se.x.size();
  std::vector<CombElt<T>> ListMiss;
  MyMatrix<T> FAC = se.GetFAC(os);
  int n_mat = FAC.rows();
  std::vector<CombElt<T>> ListAdj;
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    ListAdj.push_back(se.GetElement(se.ListNeighborData[i_mat]));
  }
  auto GetMissedGenerator = [&](int i_mat,
                                int i_facet) -> std::optional<CombElt<T>> {
    CombElt<T> TheMat = GenerateIdentity<T>(n);
    int i_mat_work = i_mat;
    int i_facet_work = i_facet;
    while (true) {
      TheMat = ProductComb(TheMat, ListAdj[i_mat_work]);
      int iFaceOpp = ai.ll_adj[i_mat_work].l_sing_adj[i_facet_work].iFaceOpp;
      int iPolyOpp = ai.ll_adj[i_mat_work].l_sing_adj[i_facet_work].iPolyOpp;
      i_mat_work = ai.ll_adj[iFaceOpp].l_sing_adj[iPolyOpp].iFaceAdj;
      i_facet_work = ai.ll_adj[iFaceOpp].l_sing_adj[iPolyOpp].iPolyAdj;
      MyVector<T> x_img = TheMat.mat.transpose() * se.x;
      if (x_img == se.x) {
        return {};
      }
      bool test = TestIntersection(FAC, TheMat, os);
      if (test) {
        return TheMat;
      }
    }
  };
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    int n_facet = ai.ll_adj[i_mat].l_sing_adj.size();
    for (int i_facet = 0; i_facet < n_facet; i_facet++) {
      std::optional<CombElt<T>> opt = GetMissedGenerator(i_mat, i_facet);
      if (opt) {
        ListMiss.push_back(*opt);
      }
    }
  }
  return ListMiss;
}

template<typename T>
void InsertAndCheckRedundancy(StepEnum<T> & se,
                              std::vector<CombElt<T>> const &l_elt_pre,
                              RecOption const &rec_option,
                              std::ostream & os) {
  std::string PrefixSave = rec_option.PrefixSave;
  std::string MethodMissingI = rec_option.MethodMissingI;
  std::string method_adjacent = rec_option.method_adjacent;
  std::string eCommand_DD = rec_option.eCommand_DD;
  std::vector<CombElt<T>> l_elt = l_elt_pre;
  os << "InsertAndCheckRedundancy before std::sort\n";
  std::sort(l_elt.begin(), l_elt.end(),
            [](CombElt<T> const &x, CombElt<T> const &y) -> bool {
              T norm_x = NormCombElt(x);
              T norm_y = NormCombElt(y);
              if (norm_x < norm_y)
                return true;
              if (norm_x > norm_y)
                return false;
              return x.mat < y.mat;
            });
  int n_elt = l_elt.size();
  for (int i_elt = 0; i_elt < n_elt; i_elt++) {
    double norm =
      UniversalScalarConversion<double, T>(NormCombElt(l_elt[i_elt]));
    os << "i_elt=" << i_elt << " norm=" << norm << "\n";
  }
  DataFAC<T> datafac;
  auto write_files = [&]() -> void {
    std::string PrefixDatafac = PrefixSave + "DATAFAC_";
    SingleData_IncrementalWrite(PrefixDatafac, datafac);
    std::string PrefixStepenum = PrefixSave + "STEPENUM_";
    SingleData_IncrementalWrite(PrefixStepenum, se);
  };
  auto insert_generator = [&](std::vector<CombElt<T>> const f_list) -> bool {
    HumanTime time;
    bool test = se.InsertGenerators(f_list, os);
    if (test) {
      datafac = se.GetDataCone(os);
    }
    if (test && datafac.eVectInt) {
      int n_remove = se.RemoveRedundancy(os);
      if (n_remove > 0) {
        datafac = se.GetDataCone(os);
      }
    }
    os << "insert_generator, time=" << time << "\n";
    write_files();
    return test;
  };
  std::unordered_set<CombElt<T>> ListTried;
  ShortVectorGroup<T> svg(se.x, l_elt);
  auto f_inverses_clear = [&]() -> bool {
    HumanTime time;
    bool DidSomething = false;
    while (true) {
      std::vector<CombElt<T>> ListMiss =
        GetMissingInverseElement(se, datafac, svg, os);
      std::vector<CombElt<T>> ListMissB;
      for (auto &eElt : ListMiss) {
        if (ListTried.count(eElt) == 0) {
          ListMissB.push_back(eElt);
          ListTried.insert(eElt);
        }
      }
      os << "|ListMiss|=" << ListMiss.size()
         << " |ListMissB|=" << ListMissB.size()
         << " |ListTried|=" << ListTried.size() << " time=" << time
         << "\n";
      if (ListMissB.size() == 0) {
        os << "Exiting the f_inverses_clear loop time=" << time << "\n";
        return DidSomething;
      }
      bool test = insert_generator(ListMissB);
      if (test) {
        DidSomething = true;
      }
    }
  };
  auto f_facet_matching = [&]() -> bool {
    std::vector<CombElt<T>> ListMiss = GetMissingFacetMatchingElement(se, datafac, method_adjacent, eCommand_DD, svg, os);
    return insert_generator(ListMiss);
  };
  auto f_coherency_update = [&]() -> bool {
    HumanTime time;
    bool DidSomething = false;
    if (datafac.eVectInt) {
      while (true) {
        bool result1 = f_inverses_clear();
        if (result1)
          DidSomething = true;
        os << "f_coherency_update, f_inverses_clear : time=" << time << "\n";
        bool result2 = f_facet_matching();
        if (result1)
          DidSomething = true;
        os << "f_coherency_update, f_facet_matching : time=" << time << "\n";
        if (!result1 && !result2)
          return DidSomething;
      }
    }
    return false;
  };
  auto f_get_candidates =
    [&](CombElt<T> const &e_elt) -> std::vector<CombElt<T>> {
    CombElt<T> e_eltInv = InverseComb(e_elt);
    std::vector<CombElt<T>> e_pair{e_elt, e_eltInv};
    if (datafac.eVectInt) {
      std::vector<CombElt<T>> n_pair;
      for (auto &TestElt : e_pair) {
        std::optional<CombElt<T>> opt =
          GetMissing_TypeI(se, datafac, TestElt, 0, MethodMissingI, os);
        if (opt) {
          CombElt<T> const &uElt1 = *opt;
          CombElt<T> uElt2 = InverseComb(uElt1);
          n_pair.push_back(uElt1);
          n_pair.push_back(uElt2);
        }
      }
      os << "|n_pair|=" << n_pair.size() << "\n";
      return n_pair;
    } else {
      return e_pair;
    }
  };
  int pos = 0;
  for (auto &e_elt : l_elt) {
    os << "       pos = " << pos << "\n";
    os << "       |known_redundant| = " << se.known_redundant.size() << "\n";
    std::vector<CombElt<T>> l_cand = f_get_candidates(e_elt);
    bool test = insert_generator(l_cand);
    if (test) {
      os << "We did something therefore we need to do a coherency update\n";
      f_coherency_update();
    }
    pos++;
  }
}

template<typename T>
std::pair<int, std::vector<std::vector<int>>>
GetGroupPresentation(StepEnum<T> const& se,
                     AdjacencyInfo<T> const &ai, std::ostream& os) {
  int n = se.x.size();
  std::vector<std::vector<int>> ListWord;
  MyMatrix<T> FAC = se.GetFAC(os);
  int n_mat = FAC.rows();
  std::vector<CombElt<T>> ListAdj;
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    ListAdj.push_back(se.GetElement(se.ListNeighborData[i_mat]));
  }
  auto InsertWordRidge = [&](int i_mat, int i_facet) -> void {
    CombElt<T> TheMat = GenerateIdentity<T>(n);
    int i_mat_work = i_mat;
    int i_facet_work = i_facet;
    std::vector<int> TheWord;
    while (true) {
      TheWord.push_back(i_mat_work);
      TheMat = ProductComb(TheMat, ListAdj[i_mat_work]);
      int iFaceOpp = ai.ll_adj[i_mat_work].l_sing_adj[i_facet_work].iFaceOpp;
      int iPolyOpp = ai.ll_adj[i_mat_work].l_sing_adj[i_facet_work].iPolyOpp;
      i_mat_work = ai.ll_adj[iFaceOpp].l_sing_adj[iPolyOpp].iFaceAdj;
      i_facet_work = ai.ll_adj[iFaceOpp].l_sing_adj[iPolyOpp].iPolyAdj;
      MyVector<T> x_img = TheMat.mat.transpose() * se.x;
      if (x_img == se.x) {
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

template <typename T>
StepEnum<T> IterativePoincareRefinement(StepEnum<T> se,
                                        RecOption const &rec_option,
                                        std::ostream& os) {
  std::string method_adjacent = rec_option.method_adjacent;
  std::string eCommand_DD = rec_option.eCommand_DD;
  bool DidSomething = false;
  auto insert_block = [&](std::vector<CombElt<T>> const &ListMiss) -> void {
    if (ListMiss.size() > 0) {
      bool test = se.InsertGenerators(ListMiss, os);
      if (test) {
        se.RemoveRedundancy(os);
        DidSomething = true;
      }
    }
  };
  int n_iter = 0;
  while (true) {
    DidSomething = false;
    //
    // Iteration Type II
    //
    AdjacencyInfo<T> ai = ComputeAdjacencyInfo(se, eCommand_DD, os);
    std::vector<CombElt<T>> ListMissII = GenerateTypeIIneighbors(se, ai, os);
    os << "|ListMissII|=" << ListMissII.size() << "\n";
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
  ListStringValues_doc["method_adjacent"] = "method_adjacent: linear_programming\n\
The available methods are linear_programming and dual_description";
  ListStringValues_doc["eCommand_DD"] = "eCommand_DD: lrs\n\
The serial program for computing the dual description. Possibilities: lrs, cdd";
  ListStringValues_doc["PrefixSave"] = "Default: unset\n\
The step enum current state";
  ListStringValues_doc["FileDataPoincare"] = "Default: unset\n\
The input file of the computation";
  ListStringValues_doc["FileO"] = "The output file of the computation";
  ListStringValues_doc["Arithmetic"] = "Default: rational\n\
Other possibilities are Qsqrt2, Qsqrt5 and RealAlgebraic=FileDesc where FileDesc is the description";
  ListStringValues_doc["Approach"] = "IncrementallyAdd or FacetAdjacencies";
  ListStringValues_doc["MethodMissingI"] = "Default: Gen1\n\
Method used for computing TypeI neighbors";
  ListStringValues_doc["MethodVertexMatching"] = "Default: none\n\
Whether to generate new elements from vertex matchings";
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
  std::string method_adjacent = BlockPROC.ListStringValues.at("method_adjacent");
  std::string eCommand_DD = BlockPROC.ListStringValues.at("eCommand_DD");
  std::string PrefixSave = BlockPROC.ListStringValues.at("PrefixSave");
  std::string FileDataPoincare =
      BlockPROC.ListStringValues.at("FileDataPoincare");
  std::string FileO = BlockPROC.ListStringValues.at("FileO");
  std::string Arithmetic = BlockPROC.ListStringValues.at("Arithmetic");
  std::string Approach = BlockPROC.ListStringValues.at("Approach");
  std::string MethodMissingI = BlockPROC.ListStringValues.at("MethodMissingI");
  std::string MethodVertexMatching =
      BlockPROC.ListStringValues.at("MethodVertexMatching");
  int n_iter_max = BlockPROC.ListIntValues.at("n_iter_max");
  int n_expand = BlockPROC.ListIntValues.at("n_expand");
  bool ComputeStabilizerPermutation =
      BlockPROC.ListBoolValues.at("ComputeStabilizerPermutation");
  bool ComputeGroupPresentation =
      BlockPROC.ListBoolValues.at("ComputeGroupPresentation");
  return {method_adjacent,
          eCommand_DD,
          PrefixSave,
          FileDataPoincare,
          FileO,
          Arithmetic,
          Approach,
          MethodMissingI,
          MethodVertexMatching,
          n_iter_max,
          n_expand,
          ComputeStabilizerPermutation,
          ComputeGroupPresentation};
}

template <typename T>
void PrintAdjacencyInfo(StepEnum<T> const &se, std::string const &FileO) {
  std::ofstream os(FileO);
  size_t n_neigh = se.ListNeighborX.size();
  os << n_neigh << "\n";
  for (size_t i_neigh = 0; i_neigh < n_neigh; i_neigh++) {
    std::pair<size_t, size_t> epair = se.ListNeighborData[i_neigh];
    CombElt<T> eElt = se.GetElement(epair);
    WriteTrackGroup(os, eElt.tg);
    os << "\n";
  }
}

void PrintGroupPresentation(
    std::ostream &os,
    std::pair<int, std::vector<std::vector<int>>> const &ThePres) {
  os << "The number of generators is " << ThePres.first << "\n";
  size_t n_word = ThePres.second.size();
  os << "The number of words is " << n_word << "\n";
  for (size_t i_word = 0; i_word < n_word; i_word++) {
    os << i_word << " : {";
    std::vector<int> const &eWord = ThePres.second[i_word];
    size_t len = eWord.size();
    for (size_t u = 0; u < len; u++) {
      if (u > 0)
        os << ",";
      os << eWord[u];
    }
    os << "}\n";
  }
}

template <typename T>
StepEnum<T> compute_step_enum(RecOption const &rec_option, std::ostream& os) {
  DataPoincare<T> dp =
      ReadDataPoincare<T>(rec_option.FileDataPoincare, rec_option.n_expand);
  os << "We have dp\n";
  auto f_init = [&]() -> StepEnum<T> {
    if (rec_option.Approach == "IncrementallyAdd") {
      StepEnum<T> se;
      se.initial_set(dp.x);
      return se;
    }
    if (rec_option.Approach == "Restart") {
      std::string PrefixStepenum = rec_option.PrefixSave + "STEPENUM_";
      std::optional<StepEnum<T>> opt = SingleData_LoadLast<StepEnum<T>>(PrefixStepenum);
      if (opt) {
        return *opt;
      } else {
        std::cerr << "Failed to find a file \n";
        throw TerminalException{1};
      }
    }
    std::cerr << "Failed to find a matching entry in compute_step_enum\n";
    throw TerminalException{1};
  };
  StepEnum<T> se = f_init();
  InsertAndCheckRedundancy(se, dp.ListGroupElt, rec_option, os);
  return IterativePoincareRefinement(se, rec_option, os);
}

template <typename T, typename Tgroup>
void full_process_type(RecOption const &rec_option, std::ostream& os) {
  os << "Beginning of full_process_type\n";
  StepEnum<T> se = compute_step_enum<T>(rec_option, os);
  os << "We have se\n";
  PrintAdjacencyInfo(se, rec_option.FileO);
  os << "se has been written to file\n";
  bool ComputeStabilizerPermutation = rec_option.ComputeStabilizerPermutation;
  bool ComputeGroupPresentation = rec_option.ComputeGroupPresentation;
  os << "ComputeStabilizerPermutation = " << ComputeStabilizerPermutation
            << "\n";
  os << "ComputeGroupPresentation = " << ComputeGroupPresentation
            << "\n";
  if (ComputeStabilizerPermutation) {
    os << "Writing the permutation group stabilizer\n";
    Tgroup GRP = se.template GetPermutationGroup<Tgroup>();
    os << "GRP=" << GRP.GapString() << "\n";
  }
  if (ComputeGroupPresentation) {
    AdjacencyInfo<T> ai = ComputeAdjacencyInfo(se, rec_option.eCommand_DD, os);
    os << "Writing the group presentation\n";
    std::pair<int, std::vector<std::vector<int>>> ThePres =
      GetGroupPresentation(se, ai, os);
    PrintGroupPresentation(os, ThePres);
  }
}

// clang-format off
#endif  // SRC_POINCARE_POLYHEDRON_POINCARE_POLYHEDRON_H_
// clang-format on

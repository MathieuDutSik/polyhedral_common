// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_MATRIXGROUPBASIC_H_
#define SRC_GROUP_MATRIXGROUPBASIC_H_

// clang-format off
#include "GRP_GroupFct.h"
#include "MAT_MatrixInt.h"
#include "MAT_MatrixMod.h"
#include "ClassicLLL.h"
#include "Timings.h"
#include "TestGroup.h"
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <set>
#include <map>
// clang-format on

#ifdef DEBUG
#define DEBUG_MATRIX_GROUP_BASIC
#endif

#ifdef TIMINGS
#define TIMINGS_MATRIX_GROUP_BASIC
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_MATRIX_GROUP_BASIC
#endif

#ifdef TRACK_INFO
#define TRACK_INFO_MATRIX_GROUP_BASIC
#endif

template <typename T>
void write_matrix_group(std::vector<MyMatrix<T>> const& list_mat, std::string const& context) {
  if (list_mat.size() == 0) {
    return;
  }
  int dim = list_mat[0].rows();
  std::string Prefix = "MatrixGroup_" + context + "_dim" + std::to_string(dim) + "_idx";
  std::string FileOut = FindAvailableFileFromPrefix(Prefix);
  WriteListMatrixFile(FileOut, list_mat);
}


template <typename T>
size_t GetRationalInvariant(std::vector<MyMatrix<T>> const &ListGen) {
  std::set<T> set_den;
  for (auto &eGen : ListGen) {
    if (!IsIntegralMatrix(eGen)) {
      int n = eGen.rows();
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          T eDen = GetDenominator(eGen(i, j));
          set_den.insert(eDen);
        }
      }
    }
  }
  std::set<T> primes;
  for (auto &eDen : set_den) {
    std::map<T, size_t> l_primes = FactorsIntMap(eDen);
    for (auto &kv : l_primes) {
      primes.insert(kv.first);
    }
  }
  T prod(1);
  for (auto &p : primes) {
    prod *= p;
  }
  return std::hash<T>()(prod);
}






template <typename T, typename Thelper>
T L1normMatrixGroup(Thelper const &helper,
                    std::vector<MyMatrix<T>> const &ListMatr) {
  int n = helper.n;
  T sum(0);
  for (auto &eMat : ListMatr) {
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        sum += T_abs(eMat(i, j));
  }
  return sum;
}

template <typename T, typename Tint, typename Thelper>
std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>>
LLLMatrixGroupReduction(Thelper const &helper,
                        std::vector<MyMatrix<T>> const &ListMatr, std::ostream& os) {
  int n = helper.n;
  MyMatrix<T> PosDefMat = IdentityMat<T>(n);
  for (auto &eMat : ListMatr) {
    if (!IsIdentity(eMat)) {
      MyMatrix<T> eProd = eMat * eMat.transpose();
      PosDefMat += eProd;
    }
  }
  LLLreduction<T, Tint> pair = LLLreducedBasis<T, Tint>(PosDefMat, os);
  MyMatrix<Tint> const &Pmat = pair.Pmat;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  std::vector<MyMatrix<T>> ListMatrNew;
  for (auto &eMat : ListMatr) {
    MyMatrix<T> eMatNew = Pmat_T * eMat * PmatInv_T;
    ListMatrNew.emplace_back(std::move(eMatNew));
  }
  return {std::move(ListMatrNew), Pmat};
}

template <typename T> T LinearSpace_GetDivisor(MyMatrix<T> const &TheSpace) {
  T TheDet = T_abs(DeterminantMat(TheSpace));
  T eDiv(1);
  int n = TheSpace.rows();
  RecSolutionIntMat<T> eCan(TheSpace);
  while (true) {
    MyMatrix<T> M = eDiv * IdentityMat<T>(n);
    bool test = eCan.is_containing_m(M);
    if (test) {
      return eDiv;
    }
#ifdef SANITY_CHECK_MATRIX_GROUP_BASIC
    if (eDiv > TheDet) {
      std::cerr << "eDiv=" << eDiv << " TheDet=" << TheDet << "\n";
      std::cerr << "TheSpace=\n";
      WriteMatrix(std::cerr, TheSpace);
      std::cerr << "Clear error in LinearSpace_GetDivisor\n";
      throw TerminalException{1};
    }
#endif
    eDiv += 1;
  }
}

template <typename T>
MyMatrix<T>
MatrixIntegral_GetInvariantSpace(int const &n,
                                 std::vector<MyMatrix<T>> const &LGen,
                                 [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_MATRIX_GROUP_BASIC
  MicrosecondTime time;
#endif
  std::vector<MyMatrix<T>> LGenTot;
  for (auto &eGen : LGen) {
    LGenTot.push_back(eGen);
    LGenTot.push_back(Inverse(eGen));
  }
  LGenTot.push_back(IdentityMat<T>(n));
  MyMatrix<T> TheSpace = IdentityMat<T>(n);
  T TheDet(1);
#ifdef DEBUG_MATRIX_GROUP_BASIC
  size_t iter = 0;
#endif
  while (true) {
    std::vector<MyVector<T>> ConcatSpace;
    for (auto &eGen : LGenTot) {
      MyMatrix<T> TheSpaceImg = TheSpace * eGen;
      for (int i = 0; i < n; i++) {
        ConcatSpace.push_back(GetMatrixRow(TheSpaceImg, i));
      }
    }
    MyMatrix<T> NewSpace1 = GetZbasis(MatrixFromVectorFamily(ConcatSpace));
    // The LLL reduction appears quite efficient
    MyMatrix<T> NewSpace = SublatticeBasisReduction(NewSpace1, os);
    T NewDet = T_abs(DeterminantMat(NewSpace));
    if (NewDet == TheDet) {
#ifdef DEBUG_MATRIX_GROUP_BASIC
      os << "MATGRPBAS: MatrixIntegral_GetInvariantSpace, NewSpace=\n";
      WriteMatrix(os, NewSpace);
      os << "MATGRPBAS: MatrixIntegral_GetInvariantSpace, returning after n_iter=" << iter << " TheDet=" << TheDet << "\n";
#endif
#ifdef TIMINGS_MATRIX_GROUP_BASIC
      os << "|MATGRPBAS: MatrixIntegral_GetInvariantSpace|=" << time << "\n";
#endif
      return NewSpace;
    }
    TheSpace = NewSpace;
    TheDet = NewDet;
#ifdef DEBUG_MATRIX_GROUP_BASIC
    os << "MATGRPBAS: MatrixIntegral_GetInvariantSpace, iter=" << iter << " TheDet=" << TheDet << "\n";
    iter += 1;
#endif
  }
}


template <typename Telt> struct CosetIterator {
  Telt id;
  size_t len;
  std::vector<std::vector<Telt>> ListListCoset;
  bool is_end;
  std::vector<Telt> ListRes;
  std::vector<size_t> ListLen;
  std::vector<size_t> ListPos;
  CosetIterator(Telt const &_id,
                std::vector<std::vector<Telt>> const &_ListListCoset,
                bool const &_is_end)
      : id(_id), len(_ListListCoset.size()), ListListCoset(_ListListCoset),
        is_end(_is_end) {
    for (auto &eListCoset : ListListCoset) {
      ListLen.push_back(eListCoset.size());
      ListPos.push_back(0);
      ListRes.push_back(id);
    }
  }
  Telt evaluate() const {
    Telt x = id;
    for (size_t u = 0; u < len; u++) {
      size_t v = len - 1 - u;
      x *= ListListCoset[v][ListPos[v]];
    }
    return x;
  }
  bool increment() {
    for (size_t u = 0; u < len; u++) {
      if (ListPos[u] < ListLen[u] - 1) {
        ListPos[u]++;
        for (size_t i = 0; i < u; i++)
          ListPos[i] = 0;
        return true;
      }
    }
    is_end = true;
    return false;
  }
  CosetIterator &operator++() {
    (void)increment();
    return *this;
  }
  CosetIterator operator++(int) {
    CosetIterator tmp = *this;
    (void)increment();
    return tmp;
  }
  Telt operator*() const { return evaluate(); }
  bool operator==(const CosetIterator<Telt> &x) const {
    if (is_end != x.is_end)
      return false;
    if (!is_end) {
      for (size_t u = 0; u < len; u++) {
        if (ListPos[u] != x.ListPos[u])
          return false;
      }
    }
    return true;
  }
  bool operator!=(const CosetIterator<Telt> &x) const {
    if (is_end != x.is_end)
      return true;
    if (!is_end) {
      for (size_t u = 0; u < len; u++) {
        if (ListPos[u] != x.ListPos[u])
          return true;
      }
    }
    return false;
  }
};

template <typename T> struct CosetDescription {
  using iterator = CosetIterator<MyMatrix<T>>;
  using const_iterator = iterator;
  int n;
  std::vector<std::vector<MyMatrix<T>>> ListListCoset;
  CosetDescription(int const &_n) : n(_n) {}
  void insert(std::vector<MyMatrix<T>> const &ListCoset) {
    if (ListCoset.size() > 1) {
      ListListCoset.push_back(ListCoset);
    }
  }
  size_t total_size() const {
    size_t total_size = 1;
    for (auto & eList : ListListCoset) {
      total_size *= eList.size();
    }
    return total_size;
  }
  void conjugate(MyMatrix<T> const &P) {
    MyMatrix<T> Pinv = Inverse(P);
    size_t len = ListListCoset.size();
    for (size_t u = 0; u < len; u++) {
      size_t n_cos = ListListCoset[u].size();
      for (size_t i_cos = 0; i_cos < n_cos; i_cos++) {
        MyMatrix<T> A = ListListCoset[u][i_cos];
        MyMatrix<T> A_cj = Pinv * A * P;
        ListListCoset[u][i_cos] = A_cj;
      }
    }
  }
  iterator get_begin() const {
    return CosetIterator(IdentityMat<T>(n), ListListCoset, false);
  }
  iterator get_end() const {
    return CosetIterator(IdentityMat<T>(n), {}, true);
  }
  const_iterator cbegin() const { return get_begin(); }
  const_iterator cend() const { return get_end(); }
  const_iterator begin() const { return get_begin(); }
  const_iterator end() const { return get_end(); }
  template <typename Tout> std::vector<MyMatrix<Tout>> expand() const {
    std::vector<MyMatrix<Tout>> ListCos;
    const_iterator iter = get_begin();
    while (iter != get_end()) {
      MyMatrix<T> eCos = *iter;
      MyMatrix<Tout> eCos_out = UniversalMatrixConversion<Tout, T>(eCos);
      ListCos.emplace_back(std::move(eCos_out));
      iter++;
    }
    return ListCos;
  }
};


// Compute Orbit of an object of type T2 under
// a group generated by elements of type T1
template <typename T1, typename T2, typename Fprod, typename Fterminate>
std::optional<std::vector<T2>>
OrbitComputation_limit(std::vector<T1> const &ListGen, T2 const &a,
                       const Fprod &f_prod, const Fterminate &f_terminate,
                       [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MATGRPBAS: Begin of OrbitComputation_limit\n";
#endif
  std::vector<T2> TheOrbit;
  std::unordered_map<T2, uint8_t> map;
  auto fInsert = [&](T2 const &u) -> bool {
    uint8_t &pos = map[u];
    if (pos == 0) {
      pos = 1;
      TheOrbit.push_back(u);
      return f_terminate(u);
    }
    return false;
  };
  if (fInsert(a)) {
    return {};
  }
  size_t start = 0;
  while (true) {
    size_t len = TheOrbit.size();
#ifdef DEBUG_MATRIX_GROUP
    os << "MATGRPBAS: OrbitComputation_limit start=" << start << " len=" << len
       << "\n";
#endif
    if (start == len) {
      break;
    }
    for (size_t i = start; i < len; i++) {
      // Doing a copy to avoid memory problem.
      T2 x = TheOrbit[i];
      for (auto &eGen : ListGen) {
        T2 u = f_prod(x, eGen);
        if (fInsert(u)) {
          return {};
        }
      }
    }
    start = len;
  }
#ifdef DEBUG_MATRIX_GROUP
  os << "MATGRPBAS: End of OrbitComputation_limit\n";
#endif
  return TheOrbit;
}


template <typename T1, typename T2, typename Fprod>
std::vector<T2> OrbitComputation(std::vector<T1> const &ListGen, T2 const &a,
                                 const Fprod &f_prod, std::ostream &os) {
  auto f_terminate = [&]([[maybe_unused]] T2 const &a) -> bool {
    return false;
  };
  std::optional<std::vector<T2>> opt =
      OrbitComputation_limit(ListGen, a, f_prod, f_terminate, os);
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!opt) {
    std::cerr << "The opt should have been assigned\n";
    throw TerminalException{1};
  }
#endif
  return *opt;
}

/*
template <typename T1, typename T2, typename Fprod>
std::vector<T2> OrbitComputationSort(std::vector<T1> const &ListGen, T2 const &a,
                                     const Fprod &f_prod, std::ostream &os) {
  std::vector<T2> orb = OrbitComputation(ListGen, a, f_prod, os);
  std::sort(orb.begin(), orb.end());
  return orb;
}
*/




template<typename T, typename Tgroup>
std::vector<MyMatrix<T>> PreImageSubgroupOneStep(std::vector<MyMatrix<T>> const& ListMatr, std::vector<typename Tgroup::Telt> const& ListPerm, MyMatrix<T> const& id_matr, Tgroup const& eGRP, std::ostream& os) {
#ifdef TIMINGS_MATRIX_GROUP_BASIC
  MicrosecondTime time;
#endif
#ifdef DEBUG_MATRIX_GROUP_BASIC
  os << "MATGRPBAS: PreImageSubgroupOneStep, |eGRP|=" << eGRP.size() << " |ListPerm|=" << ListPerm.size() << "\n";
#endif
  std::vector<MyMatrix<T>> ListMatr1 =
    permutalib::PreImageSubgroup<Tgroup, MyMatrix<T>>(ListMatr, ListPerm, id_matr, eGRP);
#ifdef TIMINGS_MATRIX_GROUP_BASIC
  os << "|MATGRPBAS: PreImageSubgroupOneStep, permutalib::PreImageSubgroup|=" << time << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP_BASIC
  os << "MATGRPBAS: PreImageSubgroupOneStep, comp(ListMatr1)=" << compute_complexity_listmat(ListMatr1) << "\n";
#endif
  std::vector<MyMatrix<T>> ListMatr2 = ExhaustiveReductionComplexityGroupMatrix<T>(ListMatr1, os);
#ifdef TIMINGS_MATRIX_GROUP_BASIC
  os << "|MATGRPBAS: PreImageSubgroupOneStep, ExhaustiveReductionComplexityGroupMatrix|=" << time << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP_BASIC
  os << "MATGRPBAS: PreImageSubgroupOneStep, comp(ListMatr2)=" << compute_complexity_listmat(ListMatr2) << "\n";
#endif
#ifdef SANITY_CHECK_MATRIX_GROUP_BASIC_DISABLE
  CheckGroupEquality<T,Tgroup>(ListMatr1, ListMatr2, os);
#endif
  return ListMatr2;
}





template<typename T, typename Tgroup>
std::vector<MyMatrix<T>> PreImageSubgroup(std::vector<MyMatrix<T>> const& ListMatr, std::vector<typename Tgroup::Telt> const& ListPerm, std::function<typename Tgroup::Telt(MyMatrix<T> const&)> f_get_perm, MyMatrix<T> const& id_matr, Tgroup const& eGRP, std::ostream& os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Telt id_perm = eGRP.get_identity();
  Tidx len = id_perm.size();
  Tgroup GRPbig(ListPerm, len);
  if (GRPbig.size() == eGRP.size()) {
    return ListMatr;
  }
#ifdef TRACK_INFO_MATRIX_GROUP_BASIC
  WriteGroupFile("GRPbig", GRPbig);
  WriteGroupFile("GRPsub", eGRP);
  WriteGroupFileGAP("GRPbig_gap", GRPbig);
  WriteGroupFileGAP("GRPsub_gap", eGRP);
#endif
#ifdef DEBUG_MATRIX_GROUP_BASIC
  os << "MATGRPBAS: PreImageSubgroup, beginning\n";
#endif
  std::vector<Tgroup> l_grp = GRPbig.GetAscendingChainSubgroup(eGRP);
  size_t len_stab = l_grp.size() - 1;
#ifdef DEBUG_MATRIX_GROUP_BASIC
  os << "MATGRPBAS: PreImageSubgroup, len_stab=" << len_stab << "\n";
  for (size_t iGRP=0; iGRP<=len_stab; iGRP++) {
    os << "MATGRPBAS: PreImageSubgroup, iGRP=" << iGRP << "/" << len_stab << " |eGRP|=" << l_grp[iGRP].size() << "\n";
  }
#endif
  std::vector<MyMatrix<T>> LGenMatr = ListMatr;
  std::vector<Telt> LGenPerm = ListPerm;
  for (size_t u=0; u<len_stab; u++) {
    size_t idx = len_stab - 1 - u;
#ifdef DEBUG_MATRIX_GROUP_BASIC
    os << "MATGRPBAS: PreImageSubgroup, len_stab=" << len_stab << " u=" << u << " idx=" << idx << "\n";
#endif
    LGenMatr = PreImageSubgroupOneStep<T,Tgroup>(LGenMatr, LGenPerm, id_matr, l_grp[idx], os);
    if (idx > 0) {
      LGenPerm.clear();
      for (auto & eMatr: LGenMatr) {
        Telt ePerm = f_get_perm(eMatr);
        LGenPerm.push_back(ePerm);
      }
    }
  }
  return LGenMatr;
}



// clang-format off
#endif  // SRC_GROUP_MATRIXGROUPBASIC_H_
// clang-format on

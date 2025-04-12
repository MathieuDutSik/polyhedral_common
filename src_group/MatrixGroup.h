// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_MATRIXGROUP_H_
#define SRC_GROUP_MATRIXGROUP_H_

// clang-format off
#include "GRP_GroupFct.h"
#include "Group.h"
#include "InvariantVectorFamily.h"
#include "MAT_MatrixInt.h"
#include "MAT_MatrixMod.h"
#include "PERM_Fct.h"
#include "Timings.h"
#include "factorizations.h"
#include "two_dim_lorentzian.h"
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <set>
#include <map>
// clang-format on

#ifdef DEBUG
#define DEBUG_MATRIX_GROUP
#endif

#ifdef TIMINGS
#define TIMINGS_MATRIX_GROUP
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_MATRIX_GROUP
#endif

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

//

template <typename T, typename Telt>
struct PreImager_Finite {
  MyMatrix<T> EXTfaithful;
  MyMatrix<T> pre_image_elt(Telt const& elt) const {
    using Tidx = typename Telt::Tidx;
    Tidx nbRow_tidx = EXTfaithful.rows();
    std::vector<Tidx> v(nbRow_tidx);
    for (Tidx i = 0; i < nbRow_tidx; i++) {
      v[i] = OnPoints(i, elt);
    }
    Telt ePermB(std::move(v));
    return FindTransformation(EXTfaithful, EXTfaithful, ePermB);
  }
};

template <typename T, typename Telt, typename Tint> struct FiniteMatrixGroupHelper {
  using PreImager = PreImager_Finite<T,Telt>;
  int n;
  MyMatrix<T> EXTfaithful;
  std::vector<MyVector<T>> ListV;
  std::unordered_map<MyVector<T>, int> MapV;
  int nbRow() const {
    return EXTfaithful.rows();
  }
  PreImager constant_pre_imager() const {
    return PreImager { EXTfaithful };
  }
  PreImager pre_imager([[maybe_unused]] std::vector<MyMatrix<T>> const& l_matr, [[maybe_unused]] std::vector<Telt> const& l_perm) const {
    return constant_pre_imager();
  }
};

//

template <typename T, typename Telt>
struct PreImager_FiniteIsotropic {
  MyMatrix<T> G;
  MyMatrix<T> EXTfaithful;
  PreImager_FiniteIsotropic(MyMatrix<T> const& _G, MyMatrix<T> const& _EXTfaithful) : G(_G), EXTfaithful(_EXTfaithful) {
  }
  MyMatrix<T> pre_image_elt(Telt const& elt) const {
    MyMatrix<T> const &Subspace1 = EXTfaithful;
    int n_rows = Subspace1.rows();
    int n_cols = Subspace1.cols();
    MyMatrix<T> Subspace2(n_rows, n_cols);
    for (int i_row = 0; i_row < n_rows; i_row++) {
      int j_row = OnPoints(i_row, elt);
      for (int i_col=0; i_col<n_cols; i_col++) {
        Subspace2(i_row, i_col) = Subspace1(j_row, i_col);
      }
    }
    std::optional<MyMatrix<T>> opt =
      LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1(G, Subspace1, G, Subspace2);
    return unfold_opt(opt, "The isotropic extension should work");
  }
};

template <typename T, typename Telt, typename Tint> struct FiniteIsotropicMatrixGroupHelper {
  using PreImager = PreImager_FiniteIsotropic<T,Telt>;
  int n;
  MyMatrix<T> G;
  MyMatrix<T> EXTfaithful;
  MyVector<T> Visotrop;
  std::vector<MyVector<T>> ListV;
  std::unordered_map<MyVector<T>, int> MapV;
  int nbRow() const {
    return EXTfaithful.rows();
  }
  PreImager constant_pre_imager() const {
    return PreImager(G, EXTfaithful);
  }
  PreImager pre_imager([[maybe_unused]] std::vector<MyMatrix<T>> const& l_matr, [[maybe_unused]] std::vector<Telt> const& l_perm) const {
    return constant_pre_imager();
  }
};

//

template <typename T, typename Telt, typename Tint>
struct PreImager_General {
private:
  permutalib::PreImagerElement<Telt,MyMatrix<T>, Tint> inner;
public:
  PreImager_General(std::vector<MyMatrix<T>> const& l_matr, std::vector<Telt> const& l_perm, int const& dim) : inner(l_matr, l_perm, IdentityMat<T>(dim)) {
  }
  MyMatrix<T> pre_image_elt(Telt const& elt) const {
    std::optional<MyMatrix<T>> opt = inner.get_preimage(elt);
    return unfold_opt(opt, "The element elt should belong to the group");
  }
};


template <typename T, typename Telt, typename Tint> struct GeneralMatrixGroupHelper {
  using PreImager = PreImager_General<T, Telt, Tint>;
  int n;
  int nbRow() const {
    return 0;
  }
  PreImager pre_imager(std::vector<MyMatrix<T>> const& l_matr, std::vector<Telt> const& l_perm) const {
    return PreImager(l_matr, l_perm, n);
  }
};

//

template <typename Thelper> struct has_determining_ext {
  static const bool value = false;
};

template <typename T, typename Telt, typename Tint>
struct has_determining_ext<FiniteMatrixGroupHelper<T, Telt, Tint>> {
  static const bool value = true;
};

template <typename T, typename Telt, typename Tint>
struct has_determining_ext<FiniteIsotropicMatrixGroupHelper<T, Telt, Tint>> {
  static const bool value = true;
};

//

template <typename T, typename Telt, typename Tint>
GeneralMatrixGroupHelper<T, Telt, Tint>
TransformHelper(GeneralMatrixGroupHelper<T, Telt, Tint> const &helper,
                [[maybe_unused]] MyMatrix<T> const &Pmat) {
  return helper;
}

template <typename T, typename Telt, typename Tint>
FiniteIsotropicMatrixGroupHelper<T, Telt, Tint>
TransformHelper(FiniteIsotropicMatrixGroupHelper<T, Telt, Tint> const &helper,
                MyMatrix<T> const &Pmat) {
  MyMatrix<T> PmatInv = Inverse(Pmat);
  MyMatrix<T> G_new = Pmat * helper.G * Pmat.transpose();
  MyMatrix<T> EXTfaithful_new = helper.EXTfaithful * PmatInv;
  MyVector<T> Visotrop_new = PmatInv.transpose() * helper.Visotrop;
  std::vector<MyVector<T>> ListV_new;
  std::unordered_map<MyVector<T>, int> MapV_new;
  int len = EXTfaithful_new.rows();
  for (int i = 0; i < len; i++) {
    MyVector<T> eV = GetMatrixRow(EXTfaithful_new, i);
    MapV_new[eV] = i;
    ListV_new.emplace_back(std::move(eV));
  }
  return {helper.n,
          std::move(G_new),
          std::move(EXTfaithful_new),
          std::move(Visotrop_new),
          std::move(ListV_new),
          std::move(MapV_new)};
}

template <typename T, typename Telt, typename Tint>
FiniteMatrixGroupHelper<T, Telt, Tint>
TransformHelper(FiniteMatrixGroupHelper<T, Telt, Tint> const &helper,
                MyMatrix<T> const &Pmat) {
  MyMatrix<T> PmatInv = Inverse(Pmat);
  MyMatrix<T> EXTfaithful_new = helper.EXTfaithful * PmatInv;
  std::vector<MyVector<T>> ListV_new;
  std::unordered_map<MyVector<T>, int> MapV_new;
  int len = EXTfaithful_new.rows();
  for (int i = 0; i < len; i++) {
    MyVector<T> eV = GetMatrixRow(EXTfaithful_new, i);
    ListV_new.push_back(eV);
    MapV_new[eV] = i;
  }
  return {helper.n, std::move(EXTfaithful_new), std::move(ListV_new),
          std::move(MapV_new)};
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

//

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
                        std::vector<MyMatrix<T>> const &ListMatr) {
  int n = helper.n;
  MyMatrix<T> PosDefMat = IdentityMat<T>(n);
  for (auto &eMat : ListMatr) {
    if (!IsIdentity(eMat)) {
      MyMatrix<T> eProd = eMat * eMat.transpose();
      PosDefMat += eProd;
    }
  }
  LLLreduction<T, Tint> pair = LLLreducedBasis<T, Tint>(PosDefMat);
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
#ifdef SANITY_CHECK_MATRIX_GROUP
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
  std::vector<MyMatrix<T>> LGenTot;
  for (auto &eGen : LGen) {
    LGenTot.push_back(eGen);
    LGenTot.push_back(Inverse(eGen));
  }
  LGenTot.push_back(IdentityMat<T>(n));
  MyMatrix<T> TheSpace = IdentityMat<T>(n);
  T TheDet(1);
#ifdef DEBUG_MATRIX_GROUP
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
    MyMatrix<T> NewSpace = GetZbasis(MatrixFromVectorFamily(ConcatSpace));
    T NewDet = T_abs(DeterminantMat(NewSpace));
    if (NewDet == TheDet) {
      return TheSpace;
    }
    TheSpace = NewSpace;
    TheDet = NewDet;
#ifdef DEBUG_MATRIX_GROUP
    std::cerr << "MAT_GRP: MatrixIntegral_GetInvariantSpace, iter=" << iter << " TheDet=" << TheDet << "\n";
    iter += 1;
#endif
  }
}

// Compute Orbit of an object of type T2 under
// a group generated by elements of type T1
template <typename T1, typename T2, typename Fprod, typename Fterminate>
std::optional<std::vector<T2>>
OrbitComputation_limit(std::vector<T1> const &ListGen, T2 const &a,
                       const Fprod &f_prod, const Fterminate &f_terminate,
                       [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin of OrbitComputation_limit\n";
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
    os << "MAT_GRP: OrbitComputation_limit start=" << start << " len=" << len
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
  os << "MAT_GRP: End of OrbitComputation_limit\n";
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

template <typename T, typename Telt, typename Tint>
FiniteMatrixGroupHelper<T, Telt, Tint>
ComputeFiniteMatrixGroupHelper(MyMatrix<T> const &EXT) {
  std::vector<MyVector<T>> ListV;
  std::unordered_map<MyVector<T>, int> MapV;
  for (int i = 0; i < EXT.rows(); i++) {
    MyVector<T> V = GetMatrixRow(EXT, i);
    ListV.push_back(V);
    MapV[V] = i;
  }
  int n_col = EXT.cols();
  return {n_col, EXT, std::move(ListV), std::move(MapV)};
}

template <typename T, typename Telt, typename Tint>
FiniteIsotropicMatrixGroupHelper<T, Telt, Tint>
ComputeFiniteIsotropicMatrixGroupHelper(MyMatrix<T> const &G,
                                        MyMatrix<T> const &EXT,
                                        MyVector<T> const &Visotrop) {
  std::vector<MyVector<T>> ListV;
  std::unordered_map<MyVector<T>, int> MapV;
  for (int i = 0; i < EXT.rows(); i++) {
    MyVector<T> V = GetMatrixRow(EXT, i);
    ListV.push_back(V);
    MapV[V] = i;
  }
  int n_col = EXT.cols();
  return {n_col, G, EXT, Visotrop, std::move(ListV), std::move(MapV)};
}

template <typename T, typename Telt, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value, Telt>::type
GetPermutationForFiniteMatrixGroup(Thelper const &helper,
                                   MyMatrix<T> const &eMatr,
                                   [[maybe_unused]] std::ostream &os) {
  using Tidx = typename Telt::Tidx;
  Tidx len = helper.EXTfaithful.rows();
  std::vector<Tidx> V(len);
  for (Tidx i = 0; i < len; i++) {
    MyVector<T> const &eV = helper.ListV[i];
    MyVector<T> Vimg = eMatr.transpose() * eV;
    V[i] = helper.MapV.at(Vimg);
  }
  return Telt(std::move(V));
}

template <typename T, typename Tmod>
Face GetFace(int const &nbRow, std::vector<MyVector<Tmod>> const &O,
             MyMatrix<T> const &TheSpace) {
  size_t Osiz = O.size();
  size_t siz = nbRow + Osiz;
  Face eFace(siz);
  RecSolutionIntMat<T> eCan(TheSpace);
  for (size_t iO = 0; iO < Osiz; iO++) {
    MyVector<T> const &eVect = UniversalVectorConversion<T, Tmod>(O[iO]);
    bool test = eCan.has_solution_v(eVect);
    if (test) {
      eFace[nbRow + iO] = 1;
    }
  }
  return eFace;
}

template <typename T, typename Telt, typename Thelper, typename Fgetperm>
inline typename std::enable_if<has_determining_ext<Thelper>::value, std::vector<Telt>>::type
MatrixIntegral_GeneratePermutationGroupA(
    std::vector<MyMatrix<T>> const &ListMatrGens,
    Thelper const &helper, Fgetperm f_get_perm, std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin MatrixIntegral_GeneratePermutationGroupA(has)\n";
#endif
  using Tidx = typename Telt::Tidx;
  int nbRow = helper.EXTfaithful.rows();
  Tidx nbRow_tidx = nbRow;
  std::vector<Telt> ListPermGenProv;
  size_t nbGen = ListMatrGens.size();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Processing nbGen=" << nbGen << " generators\n";
#endif
  for (size_t iGen = 0; iGen < nbGen; iGen++) {
    MyMatrix<T> const &eMatrGen = ListMatrGens[iGen];
    Telt ePermGen = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(
        helper, eMatrGen, os);
    Telt ePermOrbit = f_get_perm(eMatrGen);
    Tidx Osiz = ePermOrbit.size();
    Tidx siz = nbRow_tidx + Osiz;
    std::vector<Tidx> v(siz);
    for (Tidx i = 0; i < nbRow_tidx; i++) {
      v[i] = ePermGen.at(i);
    }
    for (Tidx iO = 0; iO < Osiz; iO++) {
      Tidx jO = ePermOrbit.at(iO);
      v[nbRow_tidx + iO] = nbRow_tidx + jO;
    }
    Telt eNewPerm(std::move(v));
    ListPermGenProv.emplace_back(std::move(eNewPerm));
  }
#ifdef DEBUG_MATRIX_GROUP
  if (ListPermGenProv.size() > 0) {
    Tidx siz = ListPermGenProv[0].size();
    permutalib::Group<Telt, mpz_class> GRPprov(ListPermGenProv, siz);
    os << "MAT_GRP: |GRPprov|=" << GRPprov.size() << "\n";
  } else {
    os << "MAT_GRP: |GRPprov|=1 (no generators)\n";
  }
#endif
  return ListPermGenProv;
}

template<typename T, typename Tgroup>
struct RetMI_S {
  typename Tgroup::Tint index;
  std::vector<MyMatrix<T>> LGen;
};

// We have a finite set on which the group is acting. Therefore, we can apply
// the partition backtrack algorithms
template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value,
                               RetMI_S<T,Tgroup>>::type
MatrixIntegral_Stabilizer(std::vector<typename Tgroup::Telt> const &ListPermGens,
                          std::vector<MyMatrix<T>> const& ListMatr,
                          Tgroup const &GRPperm, Thelper const &helper,
                          Face const &eFace, [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Beginning of MatrixIntegral_Stabilizer(has)\n";
#endif
  using PreImager = typename Thelper::PreImager;
  using TintGroup = typename Tgroup::Tint;
  Tgroup eStab = GRPperm.Stabilizer_OnSets(eFace);
  TintGroup index = GRPperm.size() / eStab.size();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: |eStab|=" << eStab.size() << " |eFace|=" << eFace.count() << "\n";
  os << "MAT_GRP: MatrixIntegral_Stabilizer(has), index=" << index << "\n";
#endif
  std::vector<MyMatrix<T>> LGen;
  PreImager pre_imager = helper.pre_imager(ListMatr, ListPermGens);
  for (auto &eGen : eStab.GeneratorsOfGroup()) {
    MyMatrix<T> eMatr = pre_imager.pre_image_elt(eGen);
    LGen.emplace_back(std::move(eMatr));
  }
  return {index, LGen};
}

// We have a finite set on which the group is acting. Therefore, we can apply
// the partition backtrack algorithms
template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<
    has_determining_ext<Thelper>::value,
    std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>>>::type
MatrixIntegral_Stabilizer_RightCoset(std::vector<typename Tgroup::Telt> const &ListPermGens,
                                     std::vector<MyMatrix<T>> const& ListMatr,
                                     Tgroup const &GRPperm,
                                     Thelper const &helper, Face const &eFace,
                                     [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin MatrixIntegral_Stabilizer_RightCoset(has)\n";
#endif
  using RightCosets = typename Tgroup::RightCosets;
  using PreImager = typename Thelper::PreImager;
  Tgroup eStab = GRPperm.Stabilizer_OnSets(eFace);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: |eStab|=" << eStab.size() << " |eFace|=" << eFace.count() << "\n";
#endif
  std::vector<MyMatrix<T>> ListMatrGen;
  PreImager pre_imager = helper.pre_imager(ListMatr, ListPermGens);
  for (auto &eGen : eStab.GeneratorsOfGroup()) {
    MyMatrix<T> eMatr = pre_imager.pre_image_elt(eGen);
    ListMatrGen.emplace_back(std::move(eMatr));
  }
  std::vector<MyMatrix<T>> ListRightCoset;
  RightCosets rc = GRPperm.right_cosets(eStab);
  for (auto &eCos : rc) {
    MyMatrix<T> eMatr = pre_imager.pre_image_elt(eCos);
    ListRightCoset.emplace_back(std::move(eMatr));
  }
  return {std::move(ListMatrGen), std::move(ListRightCoset)};
}

template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value,
                               std::optional<MyMatrix<T>>>::type
MatrixIntegral_RepresentativeAction(std::vector<typename Tgroup::Telt> const &ListPermGens,
                                    std::vector<MyMatrix<T>> const& ListMatr,
                                    Tgroup const &GRPperm,
                                    Thelper const &helper, Face const &eFace1,
                                    Face const &eFace2, [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Beginning of MatrixIntegral_RepresentativeAction(has)\n";
#endif
  using Telt = typename Tgroup::Telt;
  using PreImager = typename Thelper::PreImager;
  std::optional<Telt> opt = GRPperm.RepresentativeAction_OnSets(eFace1, eFace2);
  if (!opt) {
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: Exit while loop with proof that no equivalence exists\n";
#endif
    return {};
  }
  PreImager pre_imager = helper.pre_imager(ListMatr, ListPermGens);
  return pre_imager.pre_image_elt(*opt);
}

template <typename T, typename Telt, typename Thelper, typename Fgetperm>
inline typename std::enable_if<!has_determining_ext<Thelper>::value, std::vector<Telt>>::type
MatrixIntegral_GeneratePermutationGroupA(
    std::vector<MyMatrix<T>> const &ListMatrGens,
    [[maybe_unused]] Thelper const &helper,
    Fgetperm f_get_perm,
    [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin MatrixIntegral_GeneratePermutationGroupA(!has)\n";
#endif
  std::vector<Telt> ListPermGenProv;
  size_t nbGen = ListMatrGens.size();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Processing nbGen=" << nbGen << " generators\n";
#endif
  for (size_t iGen = 0; iGen < nbGen; iGen++) {
    Telt ePermGenSelect = f_get_perm(ListMatrGens[iGen]);
    ListPermGenProv.emplace_back(std::move(ePermGenSelect));
  }
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: After generation\n";
#endif
  return ListPermGenProv;
}

template<typename T, typename Tmod, typename Telt>
Telt get_permutation_from_orbit(MyMatrix<T> const& eGen, std::vector<MyVector<Tmod>> const& O, T const& TheMod, Telt const& ePermS) {
  using Tidx = typename Telt::Tidx;
  Tmod TheMod_mod = UniversalScalarConversion<Tmod, T>(TheMod);
  size_t Osiz = O.size();
  std::vector<MyVector<Tmod>> ListImage(Osiz);
  MyMatrix<Tmod> eGenMod = ModuloReductionMatrix<T,Tmod>(eGen, TheMod);
  for (size_t iV = 0; iV < Osiz; iV++) {
    MyVector<Tmod> eVect = eGenMod.transpose() * O[iV];
    ListImage[iV] = VectorMod(eVect, TheMod_mod);
  }
  Telt ePermB = Telt(SortingPerm<MyVector<Tmod>, Tidx>(ListImage));
  Telt ePermBinv = ~ePermB;
  // By the construction and above check we have
  // V1reord[i] = V1[g1.at(i)]
  // V2reord[i] = V2[g2.at(i)]
  // We have V1reord = V2reord which gets us
  // V2[i] = V1[g1 * g2^{-1}(i)]
  Telt ePermGen = ePermBinv * ePermS;
  return ePermGen;
}

template <typename T, typename Tmod, typename Telt, typename Thelper>
std::vector<Telt> MatrixIntegral_GeneratePermutationGroup(
    std::vector<MyMatrix<T>> const &ListMatrGens,
    Thelper const &helper,
    std::vector<MyVector<Tmod>> const &O, T const &TheMod,
    std::ostream &os) {
  using Tidx = typename Telt::Tidx;
  Telt ePermS = Telt(SortingPerm<MyVector<Tmod>, Tidx>(O));
  auto f_get_perm=[&](MyMatrix<T> const& eGen) -> Telt {
    return get_permutation_from_orbit(eGen, O, TheMod, ePermS);
  };
  return MatrixIntegral_GeneratePermutationGroupA<T, Telt, Thelper, decltype(f_get_perm)>(ListMatrGens, helper, f_get_perm, os);
}

// We compute the stabilizer by applying the Schreier algorithm
template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<!has_determining_ext<Thelper>::value,
                               RetMI_S<T,Tgroup>>::type
MatrixIntegral_Stabilizer(std::vector<typename Tgroup::Telt> const &ListPermGens,
                          std::vector<MyMatrix<T>> const& ListMatrGens,
                          [[maybe_unused]] Tgroup const &GRPperm,
                          Thelper const &helper, Face const &f,
                          [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin MatrixIntegral_Stabilizer(!has)\n";
#endif
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  using Tidx = typename Telt::Tidx;
  using Tobj = Face;
  using TeltMatr = MyMatrix<T>;
  MyMatrix<T> id_matr = IdentityMat<T>(helper.n);
  Tidx len = f.size();
  Tgroup GRP(ListPermGens, len);
  Tgroup stab = GRP.Stabilizer_OnSets(f);

  auto f_op = [&](Tobj const &x, Telt const &u) -> Tobj {
    return OnSets(x, u);
  };

  TintGroup index = GRP.size() / stab.size();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_Stabilizer(!has), index=" << index << "\n";
#endif
  std::vector<MyMatrix<T>> LGen =
    permutalib::PreImageSubgroupAction<Tgroup, TeltMatr, Tobj, decltype(f_op)>(
      ListMatrGens, ListPermGens, id_matr, stab, f, f_op);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: After StabilizerMatrixPermSubset\n";
#endif
  return {index, LGen};
}

// We compute the stabilizer and right cosets by applying the Schreier algorithm
template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<
    !has_determining_ext<Thelper>::value,
    std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>>>::type
MatrixIntegral_Stabilizer_RightCoset(std::vector<typename Tgroup::Telt> const &ListPermGens,
                                     std::vector<MyMatrix<T>> const& ListMatr,
                                     [[maybe_unused]] Tgroup const &GRPperm,
                                     Thelper const &helper, Face const &eFace,
                                     [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin MatrixIntegral_Stabilizer(!has)\n";
#endif
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> id_matr = IdentityMat<T>(helper.n);
  std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>> pair =
      permutalib::StabilizerRightCosetMatrixPermSubset<Telt, MyMatrix<T>, TintGroup>(
          ListMatr, ListPermGens, id_matr, eFace);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: After StabilizerRightCosetMatrixPermSubset\n";
#endif
  return pair;
}

template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<!has_determining_ext<Thelper>::value,
                               std::optional<MyMatrix<T>>>::type
MatrixIntegral_RepresentativeAction(std::vector<typename Tgroup::Telt> const &ListPermGens,
                                    std::vector<MyMatrix<T>> const& ListMatr,
                                    [[maybe_unused]] Tgroup const &GRPperm,
                                    Thelper const &helper, Face const &eFace1,
                                    Face const &eFace2,
                                    [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Beginning of MatrixIntegral_RepresentativeAction(!has)\n";
#endif
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> id_matr = IdentityMat<T>(helper.n);
  std::optional<MyMatrix<T>> opt =
      permutalib::RepresentativeActionMatrixPermSubset<Telt, MyMatrix<T>, TintGroup>(
          ListMatr, ListPermGens, id_matr, eFace1, eFace2);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Ending of MatrixIntegral_RepresentativeAction 2\n";
#endif
  return opt;
}

/*
  Direct computation of orbits.
  First level of optional is for termination or not.
  Second level is for whether we find an equivalence or not.
 */
template <typename T, typename Fterminate>
std::optional<std::optional<MyMatrix<T>>> DirectSpaceOrbit_Equivalence(
    std::vector<MyMatrix<T>> const &ListMatrGen, MyMatrix<T> const &eSpace1,
    MyMatrix<T> const &eSpace2, T const &TheMod, Fterminate const &f_terminate,
    [[maybe_unused]] std::ostream &os) {
  int n = eSpace1.rows();
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
  // Here Tpair is <Space,Repr>
  using Tpair = std::pair<MyMatrix<T>, MyMatrix<T>>;
  std::vector<Tpair> ListPair;
  ListPair.push_back({eSpace1, IdentityMat<T>(n)});
  if (f_terminate(eSpace1))
    return {};
  size_t start = 0;
  while (true) {
    size_t len = ListPair.size();
    if (start == len)
      break;
    for (size_t idx = start; idx < len; idx++) {
      Tpair const &ePair = ListPair[idx];
      for (auto &eMatrGen : ListMatrGen) {
        MyMatrix<T> eSpaceImg = ePair.first * eMatrGen;
        MyMatrix<T> eReprImg = ePair.second * eMatrGen;
        //
        MyMatrix<T> eSpaceMod = Concatenate(ePair.first, ModSpace);
        RecSolutionIntMat<T> eCan(eSpaceMod);
        if (eCan.is_containing_m(eSpace2)) {
          std::optional<MyMatrix<T>> opt = eReprImg;
          return opt;
        }
        auto fInsert = [&](Tpair const &ePair) -> bool {
          for (auto &fPair : ListPair)
            if (eCan.is_containing_m(fPair.first))
              return false;
          ListPair.push_back(ePair);
          return f_terminate(ePair.first);
        };
        if (fInsert(eSpaceImg))
          return {};
      }
    }
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: start=" << start << " len=" << len << "\n";
#endif
    start = len;
  }
  std::optional<MyMatrix<T>> opt;
  return opt;
}

template <typename T, typename Fterminate>
std::optional<std::vector<MyMatrix<T>>>
DirectSpaceOrbit_Stabilizer(std::vector<MyMatrix<T>> const &ListMatrGen,
                            MyMatrix<T> const &eSpace, T const &TheMod,
                            Fterminate const &f_terminate,
                            [[maybe_unused]] std::ostream &os) {
  int n = eSpace.rows();
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: DirectSpaceOrbit_Stabilizer, We have ModSpace\n";
#endif
  // Tpair is a pair <Space, Repr>
  using Tpair = std::pair<MyMatrix<T>, MyMatrix<T>>;
  std::vector<Tpair> ListPair;
  ListPair.push_back({eSpace, IdentityMat<T>(n)});
  if (f_terminate(eSpace)) {
    return {};
  }
  size_t start = 0;
  while (true) {
    size_t len = ListPair.size();
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: DirectSpaceOrbit_Stabilizer, start=" << start << " len=" << len << "\n";
#endif
    if (start == len) {
      break;
    }
    for (size_t idx = start; idx < len; idx++) {
      Tpair ePair = ListPair[idx]; // Copy is needed since ListPair is extended
#ifdef DEBUG_MATRIX_GROUP
      size_t jdx = 0;
#endif
      for (auto &eMatrGen : ListMatrGen) {
#ifdef DEBUG_MATRIX_GROUP
        os << "MAT_GRP: DirectSpaceOrbit_Stabilizer, jdx=" << jdx << " / " << ListMatrGen.size() << "\n";
        jdx += 1;
#endif
        MyMatrix<T> eSpaceImg = ePair.first * eMatrGen;
        if (f_terminate(eSpaceImg)) {
          return {};
        }
        //
        MyMatrix<T> eSpaceMod = Concatenate(eSpaceImg, ModSpace);
        RecSolutionIntMat<T> eCan(eSpaceMod);
        auto need_insert = [&]() -> bool {
          for (auto &fPair : ListPair) {
            if (eCan.is_containing_m(fPair.first)) {
              return false;
            }
          }
          return true;
        };
        if (need_insert()) {
          MyMatrix<T> eReprImg = ePair.second * eMatrGen;
          Tpair ePairNew{std::move(eSpaceImg), std::move(eReprImg)};
          ListPair.emplace_back(std::move(ePairNew));
        }
      }
    }
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: DirectSpaceOrbit_Stabilizer, end of loop |ListPair|=" << ListPair.size() << "\n";
#endif
    start = len;
  }
  //
  // Orbit is fine, now computing the stabilizer by using the Schreier lemma.
  //
  std::unordered_set<MyMatrix<T>> SetGen;
  size_t nPair = ListPair.size();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: DirectSpaceOrbit_Stabilizer, nPair=" << nPair << "\n";
#endif
  for (size_t iPair = 0; iPair < nPair; iPair++) {
    Tpair const &ePair = ListPair[iPair];
    for (auto &eMatrGen : ListMatrGen) {
      MyMatrix<T> eSpaceImg = ePair.first * eMatrGen;
      MyMatrix<T> eSpaceImgMod = Concatenate(eSpaceImg, ModSpace);
      RecSolutionIntMat<T> eCan(eSpaceImgMod);
      auto f_insert = [&]() -> void {
        for (size_t jPair = 0; jPair < nPair; jPair++) {
          if (eCan.is_containing_m(ListPair[jPair].first)) {
            MyMatrix<T> eGenMatr_new =
                ePair.second * eMatrGen * Inverse(ListPair[jPair].second);
            if (!IsIdentity(eGenMatr_new)) {
              SetGen.insert(eGenMatr_new);
            }
          }
        }
      };
      f_insert();
    }
  }
  std::vector<MyMatrix<T>> ListGen(SetGen.begin(), SetGen.end());
  return ListGen;
}

template <typename T, typename Tmod, typename Tgroup, typename Thelper>
inline typename std::enable_if<!has_determining_ext<Thelper>::value,
                               std::optional<std::vector<MyVector<Tmod>>>>::type
FindingSmallOrbit([[maybe_unused]] std::vector<MyMatrix<T>> const &ListMatrGen,
                  [[maybe_unused]] MyMatrix<T> const &TheSpace, T const &TheMod,
                  MyVector<T> const &x, [[maybe_unused]] Thelper const &helper,
                  std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: FindingSmallOrbit (!has), start\n";
#endif
  // No determining EXT, hard to find clever ideas.
  MyVector<Tmod> x_mod = ModuloReductionVector<T, Tmod>(x, TheMod);
  Tmod TheMod_mod = UniversalScalarConversion<Tmod, T>(TheMod);
  auto f_prod = [&](MyVector<Tmod> const &eClass,
                    MyMatrix<Tmod> const &eElt) -> MyVector<Tmod> {
    MyVector<Tmod> eVect = eElt.transpose() * eClass;
    return VectorMod(eVect, TheMod_mod);
  };
  std::vector<MyMatrix<Tmod>> ListMatrGenMod =
    ModuloReductionStdVectorMatrix<T, Tmod>(ListMatrGen, TheMod);
  return OrbitComputation(ListMatrGenMod, x_mod, f_prod, os);
}

template <typename T, typename Tmod, typename Tgroup, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value,
                               std::optional<std::vector<MyVector<Tmod>>>>::type
FindingSmallOrbit(std::vector<MyMatrix<T>> const &ListMatrGen,
                  MyMatrix<T> const &TheSpace, T const &TheMod,
                  MyVector<T> const &a, Thelper const &helper,
                  std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: FindingSmallOrbit(has), start\n";
#endif
  using Telt = typename Tgroup::Telt;
  using PreImager = typename Thelper::PreImager;
  int n = TheSpace.rows();
  // The critical number for the computation
  size_t n_limit = 60000;
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: FindingSmallOrbit, n_limit=" << n_limit << "\n";
#endif
  std::vector<MyMatrix<Tmod>> ListMatrGenMod =
        ModuloReductionStdVectorMatrix<T, Tmod>(ListMatrGen, TheMod);
  auto test_adequateness =
      [&](MyVector<T> const &x) -> std::optional<std::vector<MyVector<Tmod>>> {
    MyVector<Tmod> x_mod = ModuloReductionVector<T, Tmod>(x, TheMod);
    size_t pos = 0;
    auto f_terminate = [&]([[maybe_unused]] MyVector<Tmod> const &a) -> bool {
      pos++;
      if (pos == n_limit)
        return true;
      return false;
    };
    Tmod TheMod_mod = UniversalScalarConversion<Tmod, T>(TheMod);
    auto f_prod = [&](MyVector<Tmod> const &eClass,
                      MyMatrix<Tmod> const &eElt) -> MyVector<Tmod> {
      MyVector<Tmod> eVect = eElt.transpose() * eClass;
      return VectorMod(eVect, TheMod_mod);
    };
    return OrbitComputation_limit(ListMatrGenMod, x_mod, f_prod, f_terminate,
                                  os);
  };
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
  MyMatrix<T> TheSpaceMod = Concatenate(TheSpace, ModSpace);
  RecSolutionIntMat<T> eCan(TheSpaceMod);
  auto IsStabilized = [&](MyVector<T> const &V) -> bool {
    for (auto &eMatrGen : ListMatrGen) {
      MyVector<T> Vimg = eMatrGen.transpose() * V;
      bool test = eCan.has_solution_v(Vimg);
      if (!test) {
        return false;
      }
    }
    return true;
  };

  std::optional<std::vector<MyVector<Tmod>>> opt = test_adequateness(a);
  if (opt) {
    return *opt;
  }
  std::vector<Telt> ListPermGen;
  for (auto &eMatrGen : ListMatrGen) {
    Telt ePermGen = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(
        helper, eMatrGen, os);
    ListPermGen.emplace_back(std::move(ePermGen));
  }
  size_t len = helper.EXTfaithful.rows();
  Telt id_perm(len);
  Tgroup GRP(ListPermGen, id_perm);
  std::vector<Tgroup> ListGroup = GRP.GetAscendingChain();
  size_t len_group = ListGroup.size();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: FindingSmallOrbit, len_group=" << len_group
     << " |GRP|=" << GRP.size() << "\n";
#endif
  auto try_basis = [&](MyMatrix<T> const &TheBasis)
      -> std::optional<std::vector<MyVector<Tmod>>> {
    for (int i_row = 0; i_row < TheBasis.rows(); i_row++) {
      MyVector<T> V = GetMatrixRow(TheBasis, i_row);
      if (!IsStabilized(V)) {
        std::optional<std::vector<MyVector<Tmod>>> opt = test_adequateness(V);
        if (opt) {
          std::vector<MyVector<Tmod>> const &ListV = *opt;
#ifdef DEBUG_MATRIX_GROUP
          os << "MAT_GRP: FindingSmallOrbit, |ListV|=" << ListV.size() << "\n";
#endif
          return ListV;
        }
#ifdef DEBUG_MATRIX_GROUP
        os << "MAT_GRP: FindingSmallOrbit, Too large size at i_row=" << i_row
           << "\n";
#endif
      }
    }
    return {};
  };
  PreImager pre_imager = helper.constant_pre_imager();
  for (size_t iGroup = 0; iGroup < len_group; iGroup++) {
    size_t jGroup = len_group - 1 - iGroup;
    Tgroup const &fGRP = ListGroup[jGroup];
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: FindingSmallOrbit, iGroup=" << iGroup
       << " |fGRP|=" << fGRP.size() << "\n";
#endif
    std::vector<MyMatrix<T>> LMatr;
    for (auto &eGen : fGRP.GeneratorsOfGroup()) {
      MyMatrix<T> eMat = pre_imager.pre_image_elt(eGen);
      LMatr.emplace_back(std::move(eMat));
    }
    MyMatrix<T> InvBasis = ComputeBasisInvariantSpace(LMatr, TheSpace, TheMod);
    std::optional<std::vector<MyVector<Tmod>>> opt = try_basis(InvBasis);
    if (opt) {
      return *opt;
    }
  }
  bool AllowLargeOrbit = true;
  if (AllowLargeOrbit) {
    n_limit = std::numeric_limits<size_t>::max();
    MyMatrix<T> IdMat = IdentityMat<T>(n);
    std::optional<std::vector<MyVector<Tmod>>> opt = try_basis(IdMat);
    if (opt) {
      return *opt;
    }
  }
  return {};
}


// The space must be defining a finite index subgroup of T^n
template <typename T, typename Tmod, typename Tgroup, typename Thelper,
          typename Fstab>
std::vector<MyMatrix<T>>
LinearSpace_ModStabilizer_Tmod(std::vector<MyMatrix<T>> const &ListMatr,
                               Thelper const &helper,
                               MyMatrix<T> const &TheSpace, T const &TheMod,
                               Fstab f_stab, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n = helper.n;
#ifdef DEBUG_MATRIX_GROUP
  T TotSize(1);
  for (int i = 0; i < n; i++)
    TotSize *= TheMod;
  os << "MAT_GRP: LinearSpace_ModStabilizer_Tmod, TheMod=" << TheMod
     << "  n=" << n << " TotSize=" << TotSize << "\n";
#endif
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
  MyMatrix<T> TheSpaceMod = Concatenate(TheSpace, ModSpace);
  // This is the part of the enumeration where we have problems.
  // We have too many vectors to consider which sinks the algorithm.
  // The difficulty of the work is that we have to deal with globally
  // invariant sets of vectors. This is difficult to manipulate in
  // practice.
  //
  // Possible strategies:
  // ---Using the subchains of subgroups. This could force finding
  //    some subgroup. The strategy sometimes works. But it is
  //    ultimately not sufficient for resolving the problem.
  //
  // ---Reposition the lattices so as to create an environment
  //    where we can apply Plesken-Souvignier algorithm in order to
  //    reduce the computational size.
  //    It requires coding Intersection of subgroups, but that seems
  //    feasible.
  //
  // ---Directly computing the orbit of sublattice does not appear
  //    feasible since the equality test of sublattice is very high
  //    and the orbit itself is very large.
  //
  // ---The dream would be to identify the critical orbits in order
  //    to lift the computation. But that seems out of range.
  //
  // ---The subspaces themselves correspond to rather small sets.
  //    Could this be used? It does not seem so since the thing that
  //    is bad is the size of the group itself.
  //
  // ---Could we have some combination of strategies like in good old
  //    time of GAP and the Delaunay polytopes? Yes, but we could do better
  //    with each strategy registering prograssive improvement on the
  //    problem until completely solved.
  //
  // For the iterative improvement strategy, what we would need is some
  // encoding of the partial solution.
  // ---For the stabilizer, we have two things:
  //    ---Problem statement: helper + TheSpace
  //    ---Partial oversolution: ListMatr
  // ---For the equivalence, we have two things
  //    ---TheMod, helper1, helper2, TheSpace2
  //    ---ListMatr1, TheSpace1, MatrEquiv
  // We could look at the quotient. (Z_d)^n / TheSpace and look for point
  // stabilizers Maybe we can translate to classes easily and
  RecSolutionIntMat<T> eCan(TheSpaceMod);
  auto IsNotStabilizing = [&](std::vector<MyMatrix<T>> const &ListMatrInp)
      -> std::optional<MyVector<T>> {
    for (auto &eGen : ListMatrInp) {
      MyMatrix<T> TheSpace_img = TheSpace * eGen;
      for (int i = 0; i < n; i++) {
        MyVector<T> eVectG = GetMatrixRow(TheSpace_img, i);
        bool test = eCan.has_solution_v(eVectG);
        if (!test) {
          return GetMatrixRow(TheSpace, i);
        }
      }
    }
    return {};
  };
  std::vector<MyMatrix<T>> ListMatrRet = ListMatr;
  while (true) {
#ifdef TIMINGS_MATRIX_GROUP
    MicrosecondTime time;
#endif
    std::optional<MyVector<T>> opt = IsNotStabilizing(ListMatrRet);
#ifdef TIMINGS_MATRIX_GROUP
    os << "|MAT_GRP: IsNotStabilizing|=" << time << "\n";
#endif
    if (!opt) {
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModStabilizer_Tmod, Exiting the loop\n";
#endif
      break;
    }
    const MyVector<T> &V = *opt;
    std::optional<std::vector<MyVector<Tmod>>> opt_fso =
        FindingSmallOrbit<T, Tmod, Tgroup, Thelper>(
            ListMatrRet, TheSpace, TheMod, V, helper, os);
#ifdef TIMINGS_MATRIX_GROUP
    os << "|MAT_GRP: FindingSmallOrbit|=" << time << "\n";
#endif
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!opt_fso) {
      std::cerr << "Failed to find some entry\n";
      throw TerminalException{1};
    }
#endif
    std::vector<MyVector<Tmod>> const &O = *opt_fso;
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: LinearSpace_ModStabilizer_Tmod, |O|=" << O.size() << "\n";
#endif
    Telt ePermS = Telt(SortingPerm<MyVector<Tmod>, Tidx>(O));
    std::function<Telt(MyMatrix<T> const&)> f_get_perm=[&](MyMatrix<T> const& eGen) -> Telt {
      return get_permutation_from_orbit(eGen, O, TheMod, ePermS);
    };
    std::vector<Telt> ListPermGens =
      MatrixIntegral_GeneratePermutationGroupA<T, Telt, Thelper, decltype(f_get_perm)>(ListMatrRet, helper, f_get_perm, os);
    int nbRow = helper.nbRow();
    Tidx siz_act = nbRow + O.size();
    Tgroup GRPwork(ListPermGens, siz_act);
    Face eFace = GetFace<T, Tmod>(nbRow, O, TheSpaceMod);
#ifdef DEBUG_MATRIX_GROUP
    os << "MAT_GRP: LinearSpace_ModStabilizer_Tmod TheMod=" << TheMod
       << " |O|=" << O.size() << " |GRPwork|=" << GRPwork.size()
       << " |eFace|=" << eFace.count() << "\n";
#endif
    ListMatrRet = f_stab(ListPermGens, GRPwork, eFace, f_get_perm, ListMatrRet);
  }
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinearSpace_ModStabilizer_Tmod |ListMatrRet|=" << ListMatrRet.size() << "\n";
#endif
  return ListMatrRet;
}

template <typename T, typename Tgroup, typename Thelper, typename Fstab>
std::vector<MyMatrix<T>>
LinearSpace_ModStabilizer(std::vector<MyMatrix<T>> const &ListMatr,
                          Thelper const &helper, MyMatrix<T> const &TheSpace,
                          T const &TheMod, Fstab f_stab, std::ostream &os) {
  T max_size = (TheMod - 1) * (TheMod - 1) * TheSpace.rows();
  if (max_size < T(std::numeric_limits<uint8_t>::max())) {
    return LinearSpace_ModStabilizer_Tmod<T, uint8_t, Tgroup, Thelper>(
        ListMatr, helper, TheSpace, TheMod, f_stab, os);
  }
  if (max_size < T(std::numeric_limits<uint16_t>::max())) {
    return LinearSpace_ModStabilizer_Tmod<T, uint16_t, Tgroup, Thelper>(
        ListMatr, helper, TheSpace, TheMod, f_stab, os);
  }
  if (max_size < T(std::numeric_limits<uint32_t>::max())) {
    return LinearSpace_ModStabilizer_Tmod<T, uint32_t, Tgroup, Thelper>(
        ListMatr, helper, TheSpace, TheMod, f_stab, os);
  }
  std::cerr << "Failed to find a matching arithmetic type. Quite unlikely "
               "objectively\n";
  throw TerminalException{1};
}

template <typename T, typename Tgroup, typename Thelper, typename Fstab>
std::vector<MyMatrix<T>> LinearSpace_StabilizerGen_Kernel(
    std::vector<MyMatrix<T>> const &ListGen, Thelper const &helper,
    MyMatrix<T> const &TheSpace, Fstab f_stab, std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: det(TheSpace)=" << DeterminantMat(TheSpace) << "\n";
#endif
  RecSolutionIntMat<T> eCan(TheSpace);
  auto IsStabilizing = [&](std::vector<MyMatrix<T>> const &ListGen) -> bool {
    for (auto &eGen : ListGen) {
      MyMatrix<T> TheSpace_img = TheSpace * eGen;
      if (!eCan.is_containing_m(TheSpace_img)) {
        return false;
      }
    }
    return true;
  };
  if (IsStabilizing(ListGen)) {
    return ListGen;
  }
  T LFact = LinearSpace_GetDivisor(TheSpace);
  std::vector<T> eList = FactorsInt(LFact);
  int siz = eList.size();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinearSpace_StabilizerGen_Kernel LFact=" << LFact
     << " siz=" << siz << "\n";
#endif
  std::vector<MyMatrix<T>> ListGenRet = ListGen;
  for (int i = 1; i <= siz; i++) {
    T TheMod(1);
    for (int j = 0; j < i; j++) {
      TheMod *= eList[j];
    }
    ListGenRet = LinearSpace_ModStabilizer<T, Tgroup, Thelper, Fstab>(
        ListGenRet, helper, TheSpace, TheMod, f_stab, os);
    if (IsStabilizing(ListGenRet)) {
      return ListGenRet;
    }
  }
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsStabilizing(ListGenRet)) {
    std::cerr << "MAT_GRP: Error in LinearSpace_Stabilizer_Kernel\n";
    throw TerminalException{1};
  }
#endif
  return ListGenRet;
}

template <typename T, typename Tgroup, typename Thelper>
RetMI_S<T, Tgroup> LinearSpace_Stabilizer_Kernel(std::vector<MyMatrix<T>> const &ListGen,
                              Thelper const &helper,
                              MyMatrix<T> const &TheSpace, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  TintGroup total_index(1);
  auto f_stab = [&](std::vector<Telt> const &ListPermGens, Tgroup const &GRP,
                    Face const &eFace,
                    [[maybe_unused]] std::function<Telt(MyMatrix<T> const&)> f_get_perm,
                    std::vector<MyMatrix<T>> const& ListMatr) -> std::vector<MyMatrix<T>> {
    RetMI_S<T,Tgroup> ret = MatrixIntegral_Stabilizer<T, Tgroup, Thelper>(ListPermGens, ListMatr, GRP, helper, eFace, os);
    total_index *= ret.index;
    return ret.LGen;
  };
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Before LinearSpace_StabilizerGen_Kernel\n";
#endif
  std::vector<MyMatrix<T>> ListGenRet =
      LinearSpace_StabilizerGen_Kernel<T, Tgroup, Thelper, decltype(f_stab)>(
          ListGen, helper, TheSpace, f_stab, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: After LinearSpace_StabilizerGen_Kernel\n";
#endif
  return {total_index, ListGenRet};
}

template <typename T> struct Stab_RightCoset {
  std::vector<MyMatrix<T>> list_gen;
  CosetDescription<T> coset_desc;
};

template <typename T, typename Tgroup, typename Thelper>
Stab_RightCoset<T> LinearSpace_Stabilizer_RightCoset_Kernel(
    std::vector<MyMatrix<T>> const &l_gens, Thelper const &helper,
    MyMatrix<T> const &TheSpace, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  int n = helper.n;
  CosetDescription<T> coset(n);
  auto f_stab = [&](std::vector<Telt> const &ListPermGens, Tgroup const &GRP,
                    Face const &eFace,
                    [[maybe_unused]] std::function<Telt(MyMatrix<T> const&)> f_get_perm,
                    std::vector<MyMatrix<T>> const& ListMatr) -> std::vector<MyMatrix<T>> {
    std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>> pair =
    MatrixIntegral_Stabilizer_RightCoset<T, Tgroup, Thelper>(ListPermGens, ListMatr, GRP, helper, eFace, os);
    coset.insert(pair.second);
    return pair.first;
  };
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Before LinearSpace_StabilizerGen_Kernel\n";
#endif
  std::vector<MyMatrix<T>> l_gens_ret =
      LinearSpace_StabilizerGen_Kernel<T, Tgroup, Thelper, decltype(f_stab)>(
          l_gens, helper, TheSpace, f_stab, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: After LinearSpace_StabilizerGen_Kernel\n";
#endif
  return {std::move(l_gens_ret), coset};
}

template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<has_determining_ext<Thelper>::value,
                               std::vector<MyMatrix<T>>>::type
MatrixIntegral_PreImageSubgroup(std::vector<typename Tgroup::Telt> const &ListPermGens,
                                std::vector<MyMatrix<T>> const& ListMatr,
                                Tgroup const &eGRP, Thelper const &helper,
                                [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin MatrixIntegral_PreImageSubgroup(has) |ListPermGens|=" << ListPermGens.size() << " |gen(eGRP)|=" << eGRP.GeneratorsOfGroup().size() << "\n";
#endif
  using Telt = typename Tgroup::Telt;
  using PreImager = typename Thelper::PreImager;
  std::vector<MyMatrix<T>> ListMatrGen;
  PreImager pre_imager = helper.pre_imager(ListMatr, ListPermGens);
  std::vector<Telt> LGenSmall = eGRP.SmallGeneratingSet();
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin MatrixIntegral_PreImageSubgroup(has) |LGenSmall|=" << LGenSmall.size() << "\n";
#endif
  for (auto &eGen : eGRP.GeneratorsOfGroup()) {
    MyMatrix<T> eMatr = pre_imager.pre_image_elt(eGen);
    ListMatrGen.emplace_back(std::move(eMatr));
  }
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: |ListMatrGen|=" << ListMatrGen.size() << "\n";
#endif
  return ListMatrGen;
}

template <typename T, typename Tgroup, typename Thelper>
inline typename std::enable_if<!has_determining_ext<Thelper>::value,
                               std::vector<MyMatrix<T>>>::type
MatrixIntegral_PreImageSubgroup(std::vector<typename Tgroup::Telt> const &ListPermGens,
                                std::vector<MyMatrix<T>> const& ListMatr,
                                Tgroup const &eGRP, Thelper const &helper,
                                [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Begin MatrixIntegral_PreImageSubgroup(!has) |ListPermGens|=" << ListPermGens.size() << " |gen(eGRP)|=" << eGRP.GeneratorsOfGroup().size() << "\n";
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n_act = ListPermGens[0].size();
  Tgroup GRP_build(ListPermGens, n_act);
  os << "MAT_GRP: Begin MatrixIntegral_PreImageSubgroup(!has) |GRP_build|=" << GRP_build.size() << " |eGRP|=" << eGRP.size() << "\n";
  bool test = GRP_build.IsSubgroup(eGRP);
  os << "MAT_GRP: Begin MatrixIntegral_PreImageSubgroup(!has) IsSubgroup=" << test << "\n";
  if (false) {
    WriteGroupFile("GRP_build", GRP_build);
    WriteGroupFile("eGRP", eGRP);
    std::cerr << "Now debugging from here\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> id_matr = IdentityMat<T>(helper.n);
  std::vector<MyMatrix<T>> ListGen =
      permutalib::PreImageSubgroup<Tgroup, MyMatrix<T>>(
          ListMatr, ListPermGens, id_matr, eGRP);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: After PreImageSubgroup |ListGen|=" << ListGen.size() << "\n";
#endif
  return ListGen;
}

/*
  The dcs is the representative.
  The list stab_gens is a set of generators "x"
  of a group of transformations V
  such that
  U dcs x = U.dcs
 */
template<typename T>
struct DoubleCosetEntry {
  MyMatrix<T> cos;
  std::vector<MyMatrix<T>> stab_gens;
};

/*
  Computes the Double cosets between the stabilizer of the subspace
  and the group V in argument.
  --
  PB1: In contrast to the classical algorithm for finite groups,
  we cannot know when we are at the last step. Therefore, we have to
  compute the stabilizers all the time.
  --
  PB2: The orbit being generated can be on the big side of things.
  Therefore, we need to use the double cosets from the finite
  group algorithm. Likely, we also want to have the corresponding
  stabilizing subgroup.
  ---
  PB3: Computing in permutation groups is a winning strategy. This
  works fine for the has_determining_ext = true case.
  We need to have a single function that does the job. That is in
  order to avoid recomputing the stabilizer of eFace.
  For the has_determining_ext = false case, I need to code that
  in the permutalib code. That will be a single loop. But that is
  fine because it works that way for the cosets. Will see if further
  work is needed if failing.
 */
template <typename T, typename Tgroup, typename Thelper>
std::pair<std::vector<MyMatrix<T>>, std::vector<DoubleCosetEntry<T>>>
LinearSpace_Stabilizer_DoubleCosetStabilizer_Kernel(
    std::vector<MyMatrix<T>> const l_gens, Thelper const &helper,
    MyMatrix<T> const &TheSpace, std::vector<MyMatrix<T>> const& Vmatr_gens, std::ostream &os) {
  using PreImager = typename Thelper::PreImager;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using DoubleCosetComputer = typename Tgroup::DoubleCosetComputer;
  using DccEntry = typename Tgroup::DccEntry;
  int n = helper.n;
  std::vector<DoubleCosetEntry<T>> entries;
  entries.push_back({IdentityMat<T>(n), Vmatr_gens});
  auto f_stab = [&](std::vector<Telt> const &ListPermGens, Tgroup const &GRP,
                    Face const &eFace,
                    std::function<Telt(MyMatrix<T> const&)> f_get_perm,
                    std::vector<MyMatrix<T>> const& ListMatr) -> std::vector<MyMatrix<T>> {
    Tgroup eStab = GRP.Stabilizer_OnSets(eFace);
    PreImager pre_imager = helper.pre_imager(ListMatr, ListPermGens);
#ifdef DEBUG_DOUBLE_COSET_ENUM
    os << "MAT_GRP: Before GRP.double_coset_computer_v |GRP|=" << GRP.size() << " |eStab|=" << eStab.size() << " n_act=" << static_cast<size_t>(eStab.n_act()) << "\n";
#endif
    DoubleCosetComputer dcc_v = GRP.double_coset_computer_v(eStab);
    Tidx siz_act = eFace.size();
    std::vector<DoubleCosetEntry<T>> new_entries;
    for (auto &entry: entries) {
      MyMatrix<T> cos_inv = Inverse(entry.cos);
      std::vector<MyMatrix<T>> const& V_gens = entry.stab_gens;
      std::vector<MyMatrix<T>> V_gens_conj;
      for (auto & eGen : V_gens) {
        MyMatrix<T> NewGen = entry.cos * eGen * cos_inv;
        V_gens_conj.emplace_back(std::move(NewGen));
      }
      std::vector<Telt> ListPermGens_B =
        MatrixIntegral_GeneratePermutationGroupA<T, Telt, Thelper, std::function<Telt(MyMatrix<T> const&)>>(V_gens_conj, helper, f_get_perm, os);
      Tgroup Vperm_gens = Tgroup(ListPermGens_B, siz_act);
      std::vector<DccEntry> span_de = dcc_v.double_cosets_and_stabilizers(Vperm_gens);
      for (auto & e_de: span_de) {
        MyMatrix<T> eCos = pre_imager.pre_image_elt(e_de.cos);
        Tgroup Vred_perm(e_de.stab_gens, siz_act);
        std::vector<MyMatrix<T>> Vred_matr = MatrixIntegral_PreImageSubgroup<T,Tgroup,Thelper>(ListPermGens_B, V_gens, Vred_perm, helper, os);
        std::vector<MyMatrix<T>> Vred_matr_conj;
        for (auto & eGen : Vred_matr) {
          MyMatrix<T> NewGen = cos_inv * eGen * entry.cos;
          Vred_matr_conj.emplace_back(std::move(NewGen));
        }
        MyMatrix<T> NewCos = eCos * entry.cos;
        DoubleCosetEntry<T> new_de{std::move(NewCos), std::move(Vred_matr_conj)};
        new_entries.emplace_back(std::move(new_de));
      }
    }
#ifdef DEBUG_DOUBLE_COSET_ENUM
    os << "MAT_GRP: We found |new_entries|=" << new_entries.size() << "\n";
#endif
    entries = new_entries;
    return MatrixIntegral_PreImageSubgroup<T,Tgroup,Thelper>(ListPermGens, ListMatr, eStab, helper, os);
  };
  std::vector<MyMatrix<T>> l_gens_ret =
      LinearSpace_StabilizerGen_Kernel<T, Tgroup, Thelper, decltype(f_stab)>(
          l_gens, helper, TheSpace, f_stab, os);
  return {std::move(l_gens_ret), std::move(entries)};
}

template <typename T, typename Tgroup, typename Thelper>
std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>>
LinearSpace_Stabilizer_DoubleCoset_Kernel(
    std::vector<MyMatrix<T>> const l_gens, Thelper const &helper,
    MyMatrix<T> const &TheSpace, std::vector<MyMatrix<T>> const& Vmatr_gens, std::ostream &os) {
  std::pair<std::vector<MyMatrix<T>>, std::vector<DoubleCosetEntry<T>>> pair =
    LinearSpace_Stabilizer_DoubleCosetStabilizer_Kernel<T, Tgroup, Thelper>(l_gens, helper, TheSpace, Vmatr_gens, os);
  std::vector<MyMatrix<T>> l_dcs;
  for (auto & entry : pair.second) {
    l_dcs.push_back(entry.cos);
  }
  return {std::move(pair.first), std::move(l_dcs)};
}

template <typename T, typename Tgroup, typename Thelper>
RetMI_S<T,Tgroup> LinearSpace_Stabilizer(std::vector<MyMatrix<T>> const &ListMatr,
                       Thelper const &helper, MyMatrix<T> const &TheSpace,
                       std::ostream &os) {
  using Tint = typename underlying_ring<T>::ring_type;
  std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>> pair =
      LLLMatrixGroupReduction<T, Tint, Thelper>(helper, ListMatr);
  std::vector<MyMatrix<T>> const &ListMatrNew = pair.first;
  MyMatrix<Tint> const &Pmat = pair.second;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  MyMatrix<T> TheSpace_B = TheSpace * PmatInv_T;
  MyMatrix<T> TheSpace_C = LLLbasisReduction<T, Tint>(TheSpace_B).LattRed;
  Thelper helper_new = TransformHelper(helper, Pmat_T);
  RetMI_S<T,Tgroup> ret =
      LinearSpace_Stabilizer_Kernel<T, Tgroup, Thelper>(ListMatrNew, helper_new,
                                                        TheSpace_C, os);
  std::vector<MyMatrix<T>> const& ListMatr_B = ret.LGen;
  std::vector<MyMatrix<T>> ListMatr_C;
  for (auto &eMatr_B : ListMatr_B) {
    MyMatrix<T> eMatr_C = PmatInv_T * eMatr_B * Pmat_T;
    ListMatr_C.push_back(eMatr_C);
  }
  if (ListMatr_C.size() == 0) {
    ListMatr_C.push_back(IdentityMat<T>(helper.n));
  }
  return {ret.index, ListMatr_C};
}

template <typename T, typename Tgroup, typename Thelper>
Stab_RightCoset<T> LinearSpace_Stabilizer_RightCoset(
    std::vector<MyMatrix<T>> const &ListMatr, Thelper const &helper,
    MyMatrix<T> const &TheSpace, std::ostream &os) {
  using Tint = typename underlying_ring<T>::ring_type;
  std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>> pair =
      LLLMatrixGroupReduction<T, Tint, Thelper>(helper, ListMatr);
  std::vector<MyMatrix<T>> const &ListMatrNew = pair.first;
  MyMatrix<Tint> const &Pmat = pair.second;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  MyMatrix<T> TheSpace_B = TheSpace * PmatInv_T;
  MyMatrix<T> TheSpace_C = LLLbasisReduction<T, Tint>(TheSpace_B).LattRed;
  Thelper helper_new = TransformHelper(helper, Pmat_T);
  Stab_RightCoset<T> pairB =
      LinearSpace_Stabilizer_RightCoset_Kernel<T, Tgroup, Thelper>(
          ListMatrNew, helper_new, TheSpace_C, os);
  std::vector<MyMatrix<T>> ListMatr_C;
  for (auto &eMatr_B : pairB.list_gen) {
    MyMatrix<T> eMatr_C = PmatInv_T * eMatr_B * Pmat_T;
    ListMatr_C.push_back(eMatr_C);
  }
  if (ListMatr_C.size() == 0) {
    ListMatr_C.push_back(IdentityMat<T>(helper.n));
  }
  CosetDescription<T> coset = pairB.coset_desc;
  coset.conjugate(Pmat_T);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Returning from LinearSpace_Stabilizer_RightCoset\n";
#endif
  return {std::move(ListMatr_C), coset};
}

template <typename T, typename Tgroup, typename Thelper>
std::pair<std::vector<MyMatrix<T>>,std::vector<MyMatrix<T>>> LinearSpace_Stabilizer_DoubleCoset(
    std::vector<MyMatrix<T>> const &ListMatr, Thelper const &helper,
    MyMatrix<T> const &TheSpace, std::vector<MyMatrix<T>> const& V_gens, std::ostream &os) {
  using Tint = typename underlying_ring<T>::ring_type;
  std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>> pair =
      LLLMatrixGroupReduction<T, Tint, Thelper>(helper, ListMatr);
  std::vector<MyMatrix<T>> const &ListMatrNew = pair.first;
  MyMatrix<Tint> const &Pmat = pair.second;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  MyMatrix<T> TheSpace_B = TheSpace * PmatInv_T;
  MyMatrix<T> TheSpace_C = LLLbasisReduction<T, Tint>(TheSpace_B).LattRed;
  Thelper helper_new = TransformHelper(helper, Pmat_T);
  std::vector<MyMatrix<T>> V_gens_B;
  for (auto &eMatr_B : V_gens) {
    MyMatrix<T> eMatr_C = Pmat_T * eMatr_B * PmatInv_T;
    V_gens_B.emplace_back(std::move(eMatr_C));
  }
  std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>> pairB =
      LinearSpace_Stabilizer_DoubleCoset_Kernel<T, Tgroup, Thelper>(ListMatrNew, helper_new, TheSpace_C, V_gens_B, os);
  auto convert=[&](std::vector<MyMatrix<T>> const& l_mat) -> std::vector<MyMatrix<T>> {
    std::vector<MyMatrix<T>> l_mat_tr;
    for (auto &eMatr_B : l_mat) {
      MyMatrix<T> eMatr_C = PmatInv_T * eMatr_B * Pmat_T;
      l_mat_tr.emplace_back(std::move(eMatr_C));
    }
    return l_mat_tr;
  };
  std::vector<MyMatrix<T>> l_gen_C = convert(pairB.first);
  std::vector<MyMatrix<T>> l_cos_C = convert(pairB.second);
  return {std::move(l_gen_C), std::move(l_cos_C)};
}

template <typename T, typename Tgroup, typename Thelper>
std::pair<std::vector<MyMatrix<T>>,std::vector<DoubleCosetEntry<T>>> LinearSpace_Stabilizer_DoubleCosetStabilizer(
    std::vector<MyMatrix<T>> const &ListMatr, Thelper const &helper,
    MyMatrix<T> const &TheSpace, std::vector<MyMatrix<T>> const& V_gens, std::ostream &os) {
  using Tint = typename underlying_ring<T>::ring_type;
  std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>> pair =
      LLLMatrixGroupReduction<T, Tint, Thelper>(helper, ListMatr);
  std::vector<MyMatrix<T>> const &ListMatrNew = pair.first;
  MyMatrix<Tint> const &Pmat = pair.second;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  MyMatrix<T> TheSpace_B = TheSpace * PmatInv_T;
  MyMatrix<T> TheSpace_C = LLLbasisReduction<T, Tint>(TheSpace_B).LattRed;
  Thelper helper_new = TransformHelper(helper, Pmat_T);
  std::vector<MyMatrix<T>> V_gens_B;
  for (auto &eMatr_B : V_gens) {
    MyMatrix<T> eMatr_C = Pmat_T * eMatr_B * PmatInv_T;
    V_gens_B.emplace_back(std::move(eMatr_C));
  }
  std::pair<std::vector<MyMatrix<T>>, std::vector<DoubleCosetEntry<T>>> pairB =
      LinearSpace_Stabilizer_DoubleCosetStabilizer_Kernel<T, Tgroup, Thelper>(ListMatrNew, helper_new, TheSpace_C, V_gens_B, os);
  auto convert=[&](std::vector<MyMatrix<T>> const& l_mat) -> std::vector<MyMatrix<T>> {
    std::vector<MyMatrix<T>> l_mat_tr;
    for (auto &eMatr_B : l_mat) {
      MyMatrix<T> eMatr_C = PmatInv_T * eMatr_B * Pmat_T;
      l_mat_tr.emplace_back(std::move(eMatr_C));
    }
    return l_mat_tr;
  };
  std::vector<MyMatrix<T>> l_gen_C = convert(pairB.first);
  std::vector<DoubleCosetEntry<T>> l_dcs_C;
  for (auto & dcs : pairB.second) {
    MyMatrix<T> cos_C = PmatInv_T * dcs.cos * Pmat_T;
    std::vector<MyMatrix<T>> stab_gens_C = convert(dcs.stab_gens);
    DoubleCosetEntry<T> dcs_C{std::move(cos_C), std::move(stab_gens_C)};
    l_dcs_C.emplace_back(std::move(dcs_C));
  }
  return {std::move(l_gen_C), std::move(l_dcs_C)};
}


template <typename T>
using ResultTestModEquivalence =
    std::pair<std::vector<MyMatrix<T>>, MyMatrix<T>>;

/*
  We need a number of separate functions:
  ---The list of matrices has to be separated from the helper data like the
  EXTfaithful
  ---Function that creates the permutation representation given the ListMatrMat
  and the helper
  ---A function for computing the
 */
template <typename T, typename Tmod, typename Tgroup, typename Thelper>
std::optional<ResultTestModEquivalence<T>> LinearSpace_ModEquivalence_Tmod(
    std::vector<MyMatrix<T>> const &ListMatr, Thelper const &helper,
    bool const &NeedStabilizer, MyMatrix<T> const &TheSpace1,
    MyMatrix<T> const &TheSpace2, T const &TheMod, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  int n = TheSpace1.rows();
#ifdef DEBUG_MATRIX_GROUP
  os << "------------------------------------------------------\n";
  os << "MAT_GRP: NeedStabilizer=" << NeedStabilizer << "\n";
  os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod, TheMod=" << TheMod << "\n";
  os << "MAT_GRP: det(TheSpace1)=" << DeterminantMat(TheSpace1)
     << " det(TheSpace2)=" << DeterminantMat(TheSpace2) << "\n";
#endif
  MyMatrix<T> ModSpace = TheMod * IdentityMat<T>(n);
  MyMatrix<T> TheSpace1Mod = Concatenate(TheSpace1, ModSpace);
  MyMatrix<T> TheSpace2Mod = Concatenate(TheSpace2, ModSpace);
  std::vector<MyMatrix<T>> ListMatrRet = ListMatr;
  MyMatrix<T> eElt = IdentityMat<T>(n);
  Tmod TheMod_mod = UniversalScalarConversion<Tmod, T>(TheMod);
  auto TheAction = [&](MyVector<Tmod> const &eClass,
                       MyMatrix<Tmod> const &eElt) -> MyVector<Tmod> {
    MyVector<Tmod> eVect = eElt.transpose() * eClass;
    return VectorMod(eVect, TheMod_mod);
  };
  RecSolutionIntMat<T> eCan(TheSpace2Mod);
  auto IsEquiv =
      [&](MyMatrix<T> const &eEquiv) -> std::optional<MyVector<Tmod>> {
    MyMatrix<T> TheSpace1img = TheSpace1 * eEquiv;
    for (int i = 0; i < n; i++) {
      MyVector<T> eVect = GetMatrixRow(TheSpace1img, i);
      bool test = eCan.has_solution_v(eVect);
      if (!test) {
        return ModuloReductionVector<T, Tmod>(eVect, TheMod);
      }
    }
    return {};
  };
  auto IsStabilizing = [&](std::vector<MyMatrix<T>> const &ListGen)
      -> std::optional<MyVector<Tmod>> {
    if (!NeedStabilizer)
      return {};
    for (auto &eGen : ListGen) {
      MyMatrix<T> TheSpace2img = TheSpace2 * eGen;
      for (int i = 0; i < n; i++) {
        MyVector<T> eVect = GetMatrixRow(TheSpace2img, i);
        bool test = eCan.has_solution_v(eVect);
        if (!test) {
          return ModuloReductionVector<T, Tmod>(eVect, TheMod);
        }
      }
    }
    return {};
  };
  while (true) {
    std::optional<MyVector<Tmod>> test1 = IsEquiv(eElt);
    std::optional<MyVector<Tmod>> test2 = IsStabilizing(ListMatrRet);
    if (!test1 && !test2) {
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod, eElt and GRPwork are "
            "correct. Exiting\n";
#endif
      ResultTestModEquivalence<T> res{ListMatrRet, eElt};
      return res;
    }
    if (test1) {
      MyVector<Tmod> const &V = *test1;
      std::vector<MyMatrix<Tmod>> ListMatrRetMod =
        ModuloReductionStdVectorMatrix<T, Tmod>(ListMatrRet, TheMod);
      std::vector<MyVector<Tmod>> O =
        OrbitComputation(ListMatrRetMod, V, TheAction, os);
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod, |O|=" << O.size()
         << "\n";
#endif

      std::vector<Telt> ListPermGens =
          MatrixIntegral_GeneratePermutationGroup<T, Tmod, Telt, Thelper>(
              ListMatrRet, helper, O, TheMod, os);
      size_t nbRow = helper.nbRow();
      size_t siz_act = nbRow + O.size();
      Tgroup GRPperm(ListPermGens, siz_act);
      MyMatrix<T> TheSpace1work = TheSpace1 * eElt;
      MyMatrix<T> TheSpace1workMod = Concatenate(TheSpace1work, ModSpace);
      Face eFace1 = GetFace<T, Tmod>(nbRow, O, TheSpace1workMod);
      Face eFace2 = GetFace<T, Tmod>(nbRow, O, TheSpace2Mod);
#ifdef SANITY_CHECK_MATRIX_GROUP
      if (eFace1.count() == 0 && eFace2.count() == 0) {
        std::cerr << "Error in LinearSpace_ModEquivalence_Tmod. |eFace1| = "
                     "|eFace2| = 0\n";
        std::cerr << "Clear bug\n";
        throw TerminalException{1};
      }
#endif
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: ModEquivalence 1 TheMod=" << TheMod << " |O|=" << O.size()
         << " |GRPperm|=" << GRPperm.size() << " |eFace1|=" << eFace1.count()
         << " |eFace2|=" << eFace2.count() << "\n";
#endif
      std::optional<MyMatrix<T>> opt =
          MatrixIntegral_RepresentativeAction<T, Tgroup, Thelper>(ListPermGens, ListMatrRet, GRPperm, helper, eFace1, eFace2, os);
      if (!opt) {
#ifdef DEBUG_MATRIX_GROUP
        os << "MAT_GRP: Exit as no equivalence exixts\n";
#endif
        return {};
      }
      MyMatrix<T> const &M = *opt;
      eElt = eElt * M;
      if (!NeedStabilizer) {
        std::optional<MyVector<Tmod>> test1 = IsEquiv(eElt);
        if (!test1) {
#ifdef DEBUG_MATRIX_GROUP
          os << "MAT_GRP: eElt and GRPwork are correct. Exiting\n";
#endif
          ResultTestModEquivalence<T> res{std::move(ListMatrRet), std::move(eElt)};
          return res;
        }
      }
      RetMI_S<T,Tgroup> ret = MatrixIntegral_Stabilizer<T, Tgroup, Thelper>(ListPermGens, ListMatrRet, GRPperm, helper, eFace2, os);
      ListMatrRet = ret.LGen;
    } else {
      MyVector<Tmod> const &V = *test2;
      std::vector<MyMatrix<Tmod>> ListMatrRetMod =
        ModuloReductionStdVectorMatrix<T, Tmod>(ListMatrRet, TheMod);
      std::vector<MyVector<Tmod>> O =
          OrbitComputation(ListMatrRetMod, V, TheAction, os);
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: |O|=" << O.size() << "\n";
#endif
      std::vector<Telt> ListPermGens =
          MatrixIntegral_GeneratePermutationGroup<T, Tmod, Telt, Thelper>(
              ListMatrRet, helper, O, TheMod, os);
      int nbRow = helper.nbRow();
      size_t siz_act = nbRow + O.size();
      Tgroup GRPperm(ListPermGens, siz_act);
      Face eFace2 = GetFace<T, Tmod>(nbRow, O, TheSpace2Mod);
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: ModEquivalence 2 TheMod=" << TheMod << " |O|=" << O.size()
         << " |GRPperm|=" << GRPperm.size() << " |eFace2|=" << eFace2.count()
         << "\n";
#endif
      RetMI_S<T,Tgroup> ret = MatrixIntegral_Stabilizer<T, Tgroup, Thelper>(ListPermGens, ListMatrRet, GRPperm, helper, eFace2, os);
      ListMatrRet = ret.LGen;
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModEquivalence_Tmod, We have ListMatrRet\n";
#endif
    }
  }
}

template <typename T, typename Tgroup, typename Thelper>
std::optional<ResultTestModEquivalence<T>> LinearSpace_ModEquivalence(
    std::vector<MyMatrix<T>> const &ListMatr, Thelper const &helper,
    bool const &NeedStabilizer, MyMatrix<T> const &TheSpace1,
    MyMatrix<T> const &TheSpace2, T const &TheMod, std::ostream &os) {
  T max_size = (TheMod - 1) * (TheMod - 1) * TheSpace1.rows();
  if (max_size < T(std::numeric_limits<uint8_t>::max())) {
    return LinearSpace_ModEquivalence_Tmod<T, uint8_t, Tgroup, Thelper>(
        ListMatr, helper, NeedStabilizer, TheSpace1, TheSpace2, TheMod, os);
  }
  if (max_size < T(std::numeric_limits<uint16_t>::max())) {
    return LinearSpace_ModEquivalence_Tmod<T, uint16_t, Tgroup, Thelper>(
        ListMatr, helper, NeedStabilizer, TheSpace1, TheSpace2, TheMod, os);
  }
  if (max_size < T(std::numeric_limits<uint32_t>::max())) {
    return LinearSpace_ModEquivalence_Tmod<T, uint32_t, Tgroup, Thelper>(
        ListMatr, helper, NeedStabilizer, TheSpace1, TheSpace2, TheMod, os);
  }
  std::cerr << "Failed to find a matching arithmetic type. Quite unlikely "
               "objectively\n";
  throw TerminalException{1};
}

template <typename T, typename Tgroup, typename Thelper>
std::optional<MyMatrix<T>>
LinearSpace_Equivalence_Kernel(std::vector<MyMatrix<T>> const &ListMatr,
                               Thelper const &helper,
                               MyMatrix<T> const &InSpace1,
                               MyMatrix<T> const &InSpace2, std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Beginning of LinearSpace_Equivalence_Kernel\n";
  os << "MAT_GRP: |ListMatr|=" << ListMatr.size() << "\n";
  os << "MAT_GRP: Det(InSpace1)=" << DeterminantMat(InSpace1)
     << " Det(InSpace2)=" << DeterminantMat(InSpace2) << "\n";
#endif
  FractionMatrix<T> eRec1 = RemoveFractionMatrixPlusCoeff(InSpace1);
  FractionMatrix<T> eRec2 = RemoveFractionMatrixPlusCoeff(InSpace2);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: eRec1.TheMult=" << eRec1.TheMult
     << " eRec2.TheMult=" << eRec2.TheMult << "\n";
#endif
  if (eRec1.TheMult != eRec2.TheMult)
    return {};
  MyMatrix<T> const &TheSpace1 = eRec1.TheMat;
  MyMatrix<T> const &TheSpace2 = eRec2.TheMat;
  //
  int n = TheSpace1.rows();
  T LFact1 = LinearSpace_GetDivisor(TheSpace1);
  T LFact2 = LinearSpace_GetDivisor(TheSpace2);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LFact1 = " << LFact1 << "\n";
  os << "MAT_GRP: LFact2 = " << LFact2 << "\n";
#endif
  if (LFact1 != LFact2) {
    return {};
  }
  std::vector<T> eList = FactorsInt(LFact1);
  RecSolutionIntMat<T> eCan(TheSpace2);
  auto IsEquivalence = [&](MyMatrix<T> const &eEquiv) -> bool {
    MyMatrix<T> TheSpace1img = TheSpace1 * eEquiv;
    return eCan.is_containing_m(TheSpace1img);
  };
  std::vector<MyMatrix<T>> ListMatrWork = ListMatr;
  int siz = eList.size();
  std::vector<MyMatrix<T>> ListMatrRet = ListMatr;
  MyMatrix<T> eElt = IdentityMat<T>(n);
  for (int i = 1; i <= siz; i++) {
    if (IsEquivalence(eElt))
      return eElt;
    T TheMod(1);
    for (int j = 0; j < i; j++)
      TheMod *= eList[j];
    MyMatrix<T> TheSpace1Img = TheSpace1 * eElt;
    bool NeedStabilizer = true;
    if (i == siz)
      NeedStabilizer = false;
    std::optional<ResultTestModEquivalence<T>> opt =
        LinearSpace_ModEquivalence<T, Tgroup, Thelper>(
            ListMatrWork, helper, NeedStabilizer, TheSpace1Img, TheSpace2,
            TheMod, os);
    if (!opt) {
#ifdef DEBUG_MATRIX_GROUP
      os << "MAT_GRP: LinearSpace_ModEquivalence failed so we exit here\n";
#endif
      return {};
    }
    eElt = eElt * (opt->second);
    if (NeedStabilizer) {
      ListMatrWork = opt->first;
    }
  }
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsEquivalence(eElt)) {
    std::cerr << "Error in LinearSpace_Equivalence_Kernel\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Before returning from LinearSpace_Equivalence_Kernel, retuning eElt\n";
#endif
  return eElt;
}

template <typename T, typename Tgroup, typename Thelper>
std::optional<MyMatrix<T>>
LinearSpace_Equivalence(std::vector<MyMatrix<T>> const &ListMatr,
                        Thelper const &helper, MyMatrix<T> const &InSpace1,
                        MyMatrix<T> const &InSpace2, std::ostream &os) {
  using Tint = typename underlying_ring<T>::ring_type;
  std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>> pair =
      LLLMatrixGroupReduction<T, Tint, Thelper>(helper, ListMatr);
  std::vector<MyMatrix<T>> const &ListMatrNew = pair.first;
  MyMatrix<Tint> const &Pmat = pair.second;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  MyMatrix<T> InSpace1_B = InSpace1 * PmatInv_T;
  MyMatrix<T> InSpace2_B = InSpace2 * PmatInv_T;
  MyMatrix<T> InSpace1_C = LLLbasisReduction<T, Tint>(InSpace1_B).LattRed;
  MyMatrix<T> InSpace2_C = LLLbasisReduction<T, Tint>(InSpace2_B).LattRed;
  Thelper helper_new = TransformHelper(helper, Pmat_T);
  std::optional<MyMatrix<T>> opt =
      LinearSpace_Equivalence_Kernel<T, Tgroup, Thelper>(
          ListMatrNew, helper_new, InSpace1_C, InSpace2_C, os);
  if (!opt)
    return {};
  MyMatrix<T> RetMat = PmatInv_T * (*opt) * Pmat_T;
  return RetMat;
}

template <typename T, typename Tgroup>
RetMI_S<T,Tgroup> LinPolytopeIntegral_Automorphism_Subspaces(
    std::vector<MyMatrix<T>> const &ListMatr, MyMatrix<T> const &EXTfaithful,
    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> eBasis = GetZbasis(EXTfaithful);
  MyMatrix<T> InvBasis = Inverse(eBasis);
  MyMatrix<T> EXTbas = EXTfaithful * InvBasis;
  std::vector<MyMatrix<T>> ListMatrGens;
  for (auto &eGen : ListMatr) {
    MyMatrix<T> NewGen = eBasis * eGen * InvBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrGens.emplace_back(std::move(NewGen));
  }
  FiniteMatrixGroupHelper<T, Telt, TintGroup> helper =
    ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXTbas);
  MyMatrix<T> LattToStab = RemoveFractionMatrix(Inverse(eBasis));

  RetMI_S<T,Tgroup> ret =
      LinearSpace_Stabilizer<T, Tgroup>(ListMatrGens, helper, LattToStab, os);
  std::vector<MyMatrix<T>> ListMatrGensB;
  for (auto &eGen : ret.LGen) {
    MyMatrix<T> NewGen = InvBasis * eGen * eBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrGensB.push_back(NewGen);
  }
  return {ret.index, ListMatrGensB};
}

template <typename T, typename Tgroup>
Stab_RightCoset<T> LinPolytopeIntegral_Automorphism_RightCoset_Subspaces(
    std::vector<MyMatrix<T>> const &ListMatr, MyMatrix<T> const &EXTfaithful,
    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> eBasis = GetZbasis(EXTfaithful);
  MyMatrix<T> InvBasis = Inverse(eBasis);
  MyMatrix<T> EXTbas = EXTfaithful * InvBasis;
  std::vector<MyMatrix<T>> ListMatrGens;
  for (auto &eGen : ListMatr) {
    MyMatrix<T> NewGen = eBasis * eGen * InvBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrGens.emplace_back(std::move(NewGen));
  }
  FiniteMatrixGroupHelper<T, Telt, TintGroup> helper =
    ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXTbas);
  MyMatrix<T> LattToStab = RemoveFractionMatrix(Inverse(eBasis));
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinPolytopeIntegral_Automorphism_RightCoset_Subspaces, "
        "before LinearSpace_Stabilizer_RightCoset\n";
#endif
  Stab_RightCoset<T> pair = LinearSpace_Stabilizer_RightCoset<T, Tgroup>(
      ListMatrGens, helper, LattToStab, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinPolytopeIntegral_Automorphism_RightCoset_Subspaces, after "
        "LinearSpace_Stabilizer_RightCoset\n";
#endif
  std::vector<MyMatrix<T>> ListMatrGensB;
  for (auto &eGen : pair.list_gen) {
    MyMatrix<T> NewGen = InvBasis * eGen * eBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrGensB.emplace_back(std::move(NewGen));
  }
  pair.coset_desc.conjugate(eBasis);
  return {std::move(ListMatrGensB), pair.coset_desc};
}

template <typename T, typename Tgroup>
std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>>
LinPolytopeIntegral_Automorphism_DoubleCoset_Subspaces(std::vector<MyMatrix<T>> const &ListMatrFull, std::vector<MyMatrix<T>> const& ListMatrV, MyMatrix<T> const &EXTfaithful, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> eBasis = GetZbasis(EXTfaithful);
  MyMatrix<T> InvBasis = Inverse(eBasis);
  MyMatrix<T> EXTbas = EXTfaithful * InvBasis;
  std::vector<MyMatrix<T>> ListMatrFullGens;
  for (auto &eGen : ListMatrFull) {
    MyMatrix<T> NewGen = eBasis * eGen * InvBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrFullGens.emplace_back(std::move(NewGen));
  }
  std::vector<MyMatrix<T>> ListMatrVGens;
  for (auto &eGen : ListMatrV) {
    MyMatrix<T> NewGen = eBasis * eGen * InvBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrVGens.emplace_back(std::move(NewGen));
  }
  FiniteMatrixGroupHelper<T, Telt, TintGroup> helper =
    ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXTbas);
  MyMatrix<T> LattToStab = RemoveFractionMatrix(Inverse(eBasis));
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinPolytopeIntegral_Automorphism_RightCoset_Subspaces, "
        "before LinearSpace_Stabilizer_RightCoset\n";
#endif
  std::pair<std::vector<MyMatrix<T>>,std::vector<MyMatrix<T>>> pair =
    LinearSpace_Stabilizer_DoubleCoset<T, Tgroup>(ListMatrFullGens, helper, LattToStab, ListMatrVGens, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinPolytopeIntegral_Automorphism_DoubleCoset_Subspaces, after "
        "LinearSpace_Stabilizer_DoubleCoset\n";
#endif
  std::vector<MyMatrix<T>> ListMatrStabGens;
  for (auto &eGen : pair.first) {
    MyMatrix<T> NewGen = InvBasis * eGen * eBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrStabGens.emplace_back(std::move(NewGen));
  }
  std::vector<MyMatrix<T>> ListDoubleCosets;
  for (auto &eCos : pair.second) {
    MyMatrix<T> NewCos = InvBasis * eCos * eBasis;
    ListDoubleCosets.emplace_back(std::move(NewCos));
  }
  return {std::move(ListMatrStabGens), std::move(ListDoubleCosets)};
}

template <typename T, typename Tgroup>
std::pair<std::vector<MyMatrix<T>>, std::vector<DoubleCosetEntry<T>>>
LinPolytopeIntegral_Automorphism_DoubleCosetStabilizer_Subspaces(std::vector<MyMatrix<T>> const &ListMatrFull, std::vector<MyMatrix<T>> const& ListMatrV, MyMatrix<T> const &EXTfaithful, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> eBasis = GetZbasis(EXTfaithful);
  MyMatrix<T> InvBasis = Inverse(eBasis);
  MyMatrix<T> EXTbas = EXTfaithful * InvBasis;
  std::vector<MyMatrix<T>> ListMatrFullGens;
  for (auto &eGen : ListMatrFull) {
    MyMatrix<T> NewGen = eBasis * eGen * InvBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrFullGens.emplace_back(std::move(NewGen));
  }
  std::vector<MyMatrix<T>> ListMatrVGens;
  for (auto &eGen : ListMatrV) {
    MyMatrix<T> NewGen = eBasis * eGen * InvBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrVGens.emplace_back(std::move(NewGen));
  }
  FiniteMatrixGroupHelper<T, Telt, TintGroup> helper =
    ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXTbas);
  MyMatrix<T> LattToStab = RemoveFractionMatrix(Inverse(eBasis));
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinPolytopeIntegral_Automorphism_RightCoset_Subspaces, "
        "before LinearSpace_Stabilizer_RightCoset\n";
#endif
  std::pair<std::vector<MyMatrix<T>>,std::vector<DoubleCosetEntry<T>>> pair =
    LinearSpace_Stabilizer_DoubleCosetStabilizer<T, Tgroup>(ListMatrFullGens, helper, LattToStab, ListMatrVGens, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: LinPolytopeIntegral_Automorphism_DoubleCoset_Subspaces, after "
        "LinearSpace_Stabilizer_DoubleCoset\n";
#endif
  std::vector<MyMatrix<T>> ListMatrStabGens;
  for (auto &eGen : pair.first) {
    MyMatrix<T> NewGen = InvBasis * eGen * eBasis;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(NewGen)) {
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
#endif
    ListMatrStabGens.emplace_back(std::move(NewGen));
  }
  std::vector<DoubleCosetEntry<T>> ListDoubleCosetStabilizer;
  for (auto &dcs : pair.second) {
    MyMatrix<T> NewCos = InvBasis * dcs.cos * eBasis;
    std::vector<MyMatrix<T>> new_stab_gens;
    for (auto & eGen : dcs.stab_gens) {
      MyMatrix<T> NewGen = InvBasis * eGen * eBasis;
      new_stab_gens.emplace_back(std::move(NewGen));
    }
    DoubleCosetEntry<T> new_dcs{std::move(NewCos), std::move(new_stab_gens)};
    ListDoubleCosetStabilizer.emplace_back(std::move(new_dcs));
  }
  return {std::move(ListMatrStabGens), std::move(ListDoubleCosetStabilizer)};
}

template <typename T>
std::vector<MyMatrix<T>>
ConjugateListGeneratorsTestInt(MyMatrix<T> const &Pmat,
                               std::vector<MyMatrix<T>> const &LGen) {
  std::vector<MyMatrix<T>> LGen2;
  MyMatrix<T> PmatInv = Inverse(Pmat);
  for (auto &eGen1 : LGen) {
    MyMatrix<T> eGen2 = PmatInv * eGen1 * Pmat;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(eGen2)) {
      std::cerr << "MAT_GRP: Error in ConjugateListGeneratorsTestInt\n";
      std::cerr << "MAT_GRP: The matrix eGen2 should be integral\n";
      std::cerr << "MAT_GRP: eGen2=\n";
      WriteMatrix(std::cerr, eGen2);
      throw TerminalException{1};
    }
#endif
    LGen2.emplace_back(std::move(eGen2));
  }
  return LGen2;
}

/*
  Compute the intersection of G \cap GL_n(Z)
 */
template <typename T, typename Tint, typename Tgroup>
RetMI_S<Tint,Tgroup> MatrixIntegral_Stabilizer_General(
    int const &n, std::vector<MyMatrix<T>> const &LGen1, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  using Thelper = GeneralMatrixGroupHelper<T, Telt, TintGroup>;
  MyMatrix<T> InvariantSpace = MatrixIntegral_GetInvariantSpace(n, LGen1, os);
  MyMatrix<T> InvInvariantSpace = Inverse(InvariantSpace);
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsIntegralMatrix(InvInvariantSpace)) {
    std::cerr << "The matrix InvInvariantSpace should be integral\n";
    throw TerminalException{1};
  }
#endif
  std::vector<MyMatrix<T>> LGen2 =
      ConjugateListGeneratorsTestInt(InvInvariantSpace, LGen1);
  Thelper helper{n};
  RetMI_S<T,Tgroup> ret =
    LinearSpace_Stabilizer<T, Tgroup, Thelper>(
          LGen2, helper, InvInvariantSpace, os);
  std::vector<MyMatrix<Tint>> LGen4;
  for (auto &eGen3 : ret.LGen) {
    MyMatrix<T> eGen4_T = InvInvariantSpace * eGen3 * InvariantSpace;
#ifdef SANITY_CHECK_MATRIX_GROUP
    if (!IsIntegralMatrix(eGen4_T)) {
      std::cerr << "The matrix eGen4_T should be integral\n";
      throw TerminalException{1};
    }
#endif
    MyMatrix<Tint> eGen4 = UniversalMatrixConversion<Tint, T>(eGen4_T);
    LGen4.emplace_back(std::move(eGen4));
  }
  return {ret.index, LGen4};
}

// Instead of Z^n, we now want "Sublattice" to be preserved.
template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<T>> MatrixIntegral_Stabilizer_Sublattice(
    MyMatrix<T> const& Sublattice, std::vector<MyMatrix<T>> const &LGen1, std::ostream &os) {
  MyMatrix<T> SublatticeInv = Inverse(Sublattice);
  int n = SublatticeInv.rows();
  std::vector<MyMatrix<T>> LGen2;
  for (auto & eGen1 : LGen1) {
    MyMatrix<T> eGen2 = SublatticeInv * eGen1 * Sublattice;
    LGen2.push_back(eGen2);
  }
  std::vector<MyMatrix<Tint>> LGen3 = MatrixIntegral_Stabilizer_General<T,Tint,Tgroup>(n, LGen2, os);
  std::vector<MyMatrix<T>> LGen4;
  for (auto & eGen3 : LGen3) {
    MyMatrix<T> eGen3_T = UniversalMatrixConversion<T,Tint>(eGen3);
    MyMatrix<T> eGen4 = Sublattice * eGen3_T * SublatticeInv;
    LGen4.push_back(eGen4);
  }
  return LGen4;
}


// Returns an element g in GRPrat such that   g * EquivRat   in GL(n,Z)
template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
MatrixIntegral_Equivalence_General(std::vector<MyMatrix<T>> const &LGen1,
                                   MyMatrix<T> const &EquivRat,
                                   std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  int n = EquivRat.rows();
  MyMatrix<T> TheSpace = MatrixIntegral_GetInvariantSpace(n, LGen1, os);
  MyMatrix<T> TheSpaceInv = Inverse(TheSpace);
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsIntegralMatrix(TheSpaceInv)) {
    std::cerr << "The matrix InvInvariantSpace should be integral\n";
    throw TerminalException{1};
  }
#endif
  std::vector<MyMatrix<T>> LGen2 =
      ConjugateListGeneratorsTestInt(TheSpaceInv, LGen1);
  // We search for g in GRPrat s.t. g * EquivRat in GL_n(Z).
  // So, we search g in GRPrat s.t. Z^n * g * EquivRat = Z^n
  // Writing g = TheSpaceInv g_int TheSpace we get
  // TheSpaceInv g TheSpace EquivRat = Z^n
  // Or TheSpaceInv g = Inverse(TheSpace * EquivRat)
  MyMatrix<T> TheSpaceImg = TheSpace * EquivRat;
  MyMatrix<T> TheSpaceImgInv = Inverse(TheSpaceImg);
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsIntegralMatrix(TheSpaceImgInv)) {
    std::cerr << "The matrix TheSpaceImgInv should be integral\n";
    throw TerminalException{1};
  }
#endif
  using Thelper = GeneralMatrixGroupHelper<T, Telt, TintGroup>;
  Thelper helper{n};
  std::optional<MyMatrix<T>> opt =
      LinearSpace_Equivalence_Kernel<T, Tgroup, Thelper>(
          LGen2, helper, TheSpaceInv, TheSpaceImgInv, os);
  if (!opt) {
    return {};
  }
  MyMatrix<T> const &eSpaceEquiv = *opt;
  MyMatrix<T> eMatFinal = TheSpaceInv * eSpaceEquiv * TheSpace;
  MyMatrix<T> eProd_T = eMatFinal * EquivRat;
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsIntegralMatrix(eProd_T)) {
    std::cerr << "The matrix should be integral\n";
    throw TerminalException{1};
  }
#endif
  return UniversalMatrixConversion<Tint, T>(eProd_T);
}

// Find a matrix g in GRPrat such that   EquivRat * g   in   GL(n,Z)
template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
MatrixIntegral_Equivalence_Bis_General(std::vector<MyMatrix<T>> const &GRPrat,
                                       MyMatrix<T> const &EquivRat,
                                       std::ostream &os) {
  MyMatrix<T> EquivRatInv = Inverse(EquivRat);
  std::optional<MyMatrix<Tint>> opt =
      MatrixIntegral_Equivalence_General<T, Tint, Tgroup>(GRPrat, EquivRatInv,
                                                          os);
  if (!opt) {
    return {};
  }
  MyMatrix<Tint> const &TheSol = *opt;
  // So we have TheSol = g * Inverse(EquivRat) in GL(n,Z)
  // Inverse(TheSol) = EquivRat * g in GL(n,Z)
  return Inverse(TheSol);
}

template <typename T, typename Tgroup>
std::vector<MyMatrix<T>> MatrixIntegral_RightCosets_General(
    int const &n, std::vector<MyMatrix<T>> const &LGen1, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  using Thelper = GeneralMatrixGroupHelper<T, Telt, TintGroup>;
  MyMatrix<T> InvariantSpace = MatrixIntegral_GetInvariantSpace(n, LGen1, os);
  MyMatrix<T> InvInvariantSpace = Inverse(InvariantSpace);
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsIntegralMatrix(InvInvariantSpace)) {
    std::cerr << "MAT_GRP: The matrix InvInvariantSpace should be integral\n";
    throw TerminalException{1};
  }
#endif
  std::vector<MyMatrix<T>> LGen2 =
      ConjugateListGeneratorsTestInt(InvInvariantSpace, LGen1);
  Thelper helper{n};
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_RightCosets_General, before "
        "LinearSpace_Stabilizer_RightCoset\n";
#endif
  Stab_RightCoset<T> stab_rightcoset =
      LinearSpace_Stabilizer_RightCoset<T, Tgroup, Thelper>(
          LGen2, helper, InvInvariantSpace, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_RightCosets_General, after "
        "LinearSpace_Stabilizer_RightCoset\n";
#endif
  using Iter = typename CosetDescription<T>::const_iterator;
  Iter iter = stab_rightcoset.coset_desc.begin();
  Iter end = stab_rightcoset.coset_desc.end();
  std::vector<MyMatrix<T>> LCoset2;
  while (iter != end) {
    MyMatrix<T> const &eCos1 = *iter;
    MyMatrix<T> eCos2 = InvInvariantSpace * eCos1 * InvariantSpace;
    LCoset2.push_back(eCos2);
    iter++;
  }
  return LCoset2;
}

template <typename T, typename Tgroup>
std::pair<std::vector<MyMatrix<T>>,std::vector<MyMatrix<T>>>
MatrixIntegral_DoubleCosets_General(
    int const &n, std::vector<MyMatrix<T>> const &LGenG1,
    std::vector<MyMatrix<T>> const &LGenV1,
    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  using Thelper = GeneralMatrixGroupHelper<T, Telt, TintGroup>;
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: We have LGenG1=\n";
  WriteListMatrix(os, LGenG1);
  os << "MAT_GRP: We have LGenV1=\n";
  WriteListMatrix(os, LGenV1);
#endif
  MyMatrix<T> InvariantSpace = MatrixIntegral_GetInvariantSpace(n, LGenG1, os);
  MyMatrix<T> InvInvariantSpace = Inverse(InvariantSpace);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: We have |InvariantSpace|=" << DeterminantMat(InvariantSpace) << "\n";
  os << "MAT_GRP: We have InvariantSpace=\n";
  WriteMatrix(os, InvariantSpace);
  os << "MAT_GRP: We have InvInvariantSpace=\n";
  WriteMatrix(os, InvInvariantSpace);
#endif
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsIntegralMatrix(InvInvariantSpace)) {
    std::cerr << "MAT_GRP: The matrix InvInvariantSpace should be integral\n";
    throw TerminalException{1};
  }
#endif
  std::vector<MyMatrix<T>> LGenG2 =
      ConjugateListGeneratorsTestInt(InvInvariantSpace, LGenG1);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: We have LGenG2\n";
#endif
  std::vector<MyMatrix<T>> LGenV2 =
      ConjugateListGeneratorsTestInt(InvInvariantSpace, LGenV1);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: We have LGenV2\n";
#endif
  Thelper helper{n};
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_DoubleCosets_General, before "
        "LinearSpace_Stabilizer_DoubleCoset\n";
#endif
  std::pair<std::vector<MyMatrix<T>>,std::vector<MyMatrix<T>>> pair =
      LinearSpace_Stabilizer_DoubleCoset<T, Tgroup, Thelper>(LGenG2, helper, InvInvariantSpace, LGenV2, os);
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: MatrixIntegral_DoubleCosets_General, after "
        "LinearSpace_Stabilizer_DoubleCoset\n";
#endif
  std::vector<MyMatrix<T>> LGenRet;
  for (auto & eGen1 : pair.first) {
    MyMatrix<T> eGen2 = InvInvariantSpace * eGen1 * InvariantSpace;
    LGenRet.push_back(eGen2);
  }
  std::vector<MyMatrix<T>> LCosRet;
  for (auto & eCos1 : pair.second) {
    MyMatrix<T> eCos2 = InvInvariantSpace * eCos1 * InvariantSpace;
    LCosRet.push_back(eCos2);
  }
  return {LGenRet, LCosRet};
}



template <typename T, typename Tgroup>
Tgroup LinPolytopeIntegral_Stabilizer_Method8(MyMatrix<T> const &EXT_T,
                                              Tgroup const &GRPisom,
                                              std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  int nbVert = EXT_T.rows();
  std::vector<MyMatrix<T>> ListMatrGen;
  for (auto &eGen : GRPisom.GeneratorsOfGroup()) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    ListMatrGen.push_back(eMat);
  }
  using Thelper = FiniteMatrixGroupHelper<T, Telt, TintGroup>;
  Thelper helper = ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXT_T);
  RetMI_S<T,Tgroup> ret =
      LinPolytopeIntegral_Automorphism_Subspaces<T, Tgroup>(ListMatrGen, EXT_T,
                                                            os);
  std::vector<Telt> ListPermGens;
  for (auto &eMatr : ret.LGen) {
    Telt ePermGen = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(
        helper, eMatr, os);
    ListPermGens.emplace_back(std::move(ePermGen));
  }
  return Tgroup(ListPermGens, nbVert);
}

template <typename T, typename Tgroup>
std::pair<Tgroup, std::vector<typename Tgroup::Telt>>
LinPolytopeIntegral_Stabilizer_RightCoset_Method8(MyMatrix<T> const &EXT_T,
                                                  Tgroup const &GRPisom,
                                                  std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  int nbVert = EXT_T.rows();
  std::vector<MyMatrix<T>> ListMatrGen;
  for (auto &eGen : GRPisom.GeneratorsOfGroup()) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    ListMatrGen.push_back(eMat);
  }
  using Thelper = FiniteMatrixGroupHelper<T, Telt, TintGroup>;
  Thelper helper = ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXT_T);
  Stab_RightCoset<T> pair =
      LinPolytopeIntegral_Automorphism_RightCoset_Subspaces<T, Tgroup>(
          ListMatrGen, EXT_T, os);
  std::vector<Telt> ListPermGens;
  for (auto &eMatr : pair.list_gen) {
    Telt ePermGen = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(
        helper, eMatr, os);
    ListPermGens.emplace_back(std::move(ePermGen));
  }
  Tgroup GRPret(ListPermGens, nbVert);
  std::vector<Telt> RightCoset;
  for (auto eMatr : pair.coset_desc) {
    Telt eCos = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(helper, eMatr, os);
    RightCoset.emplace_back(std::move(eCos));
  }
  return {std::move(GRPret), std::move(RightCoset)};
}

template <typename T, typename Tgroup>
std::pair<Tgroup, std::vector<typename Tgroup::Telt>>
LinPolytopeIntegral_Stabilizer_DoubleCoset_Method8(MyMatrix<T> const &EXT_T,
                                                   Tgroup const &GRPfull,
                                                   Tgroup const &GrpV,
                                                   std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  int nbVert = EXT_T.rows();
  std::vector<MyMatrix<T>> ListMatrGenFull;
  for (auto &eGen : GRPfull.GeneratorsOfGroup()) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    ListMatrGenFull.push_back(eMat);
  }
  std::vector<MyMatrix<T>> ListMatrGenV;
  for (auto &eGen : GrpV.GeneratorsOfGroup()) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    ListMatrGenV.push_back(eMat);
  }
  using Thelper = FiniteMatrixGroupHelper<T, Telt, TintGroup>;
  Thelper helper = ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXT_T);
  std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>> pair =
      LinPolytopeIntegral_Automorphism_DoubleCoset_Subspaces<T, Tgroup>(ListMatrGenFull, ListMatrGenV, EXT_T, os);
  std::vector<Telt> ListPermGens;
  for (auto &eMatr : pair.first) {
    Telt ePermGen = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(
        helper, eMatr, os);
    ListPermGens.emplace_back(std::move(ePermGen));
  }
  Tgroup GRPret(ListPermGens, nbVert);
  std::vector<Telt> DoubleCosets;
  for (auto eMatr : pair.second) {
    Telt eCos = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(helper, eMatr, os);
    DoubleCosets.emplace_back(std::move(eCos));
  }
  return {std::move(GRPret), std::move(DoubleCosets)};
}

template<typename Telt>
struct PairCosetStabGens {
  Telt cos;
  std::vector<Telt> stab_gens;
};

template <typename T, typename Tgroup>
std::pair<Tgroup, std::vector<PairCosetStabGens<typename Tgroup::Telt>>>
LinPolytopeIntegral_Stabilizer_DoubleCosetStabilizer_Method8(MyMatrix<T> const &EXT_T,
                                                             Tgroup const &GRPfull,
                                                             Tgroup const &GrpV,
                                                             std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  int nbVert = EXT_T.rows();
  std::vector<MyMatrix<T>> ListMatrGenFull;
  for (auto &eGen : GRPfull.GeneratorsOfGroup()) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    ListMatrGenFull.push_back(eMat);
  }
  std::vector<MyMatrix<T>> ListMatrGenV;
  for (auto &eGen : GrpV.GeneratorsOfGroup()) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    ListMatrGenV.push_back(eMat);
  }
  using Thelper = FiniteMatrixGroupHelper<T, Telt, TintGroup>;
  Thelper helper = ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXT_T);
  std::pair<std::vector<MyMatrix<T>>, std::vector<DoubleCosetEntry<T>>> pair =
      LinPolytopeIntegral_Automorphism_DoubleCosetStabilizer_Subspaces<T, Tgroup>(ListMatrGenFull, ListMatrGenV, EXT_T, os);
  std::vector<Telt> ListPermGens;
  for (auto &eMatr : pair.first) {
    Telt ePermGen = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(
        helper, eMatr, os);
    ListPermGens.emplace_back(std::move(ePermGen));
  }
  Tgroup GRPret(ListPermGens, nbVert);
  std::vector<PairCosetStabGens<Telt>> DoubleCosetStabilizer;
  for (auto eDCS : pair.second) {
    Telt NewCos = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(helper, eDCS.cos, os);
    std::vector<Telt> new_stab_gens;
    for (auto & eMatr : eDCS.stab_gens) {
      Telt NewGen = GetPermutationForFiniteMatrixGroup<T, Telt, Thelper>(helper, eMatr, os);
      new_stab_gens.emplace_back(std::move(NewGen));
    }
    PairCosetStabGens<Telt> NewDCS{std::move(NewCos), std::move(new_stab_gens)};
    DoubleCosetStabilizer.emplace_back(std::move(NewDCS));
  }
  return {std::move(GRPret), std::move(DoubleCosetStabilizer)};
}

template <typename T, typename Tgroup>
std::optional<MyMatrix<T>> LinPolytopeIntegral_Isomorphism_Subspaces(
    MyMatrix<T> const &EXT1_T, MyMatrix<T> const &EXT2_T,
    std::vector<MyMatrix<T>> const &ListMatrGens2,
    typename Tgroup::Telt const &eEquiv, std::ostream &os) {
#ifdef DEBUG_MATRIX_GROUP
  os << "MAT_GRP: Beginning of LinPolytopeIntegral_Isomorphism_Subspaces\n";
#endif
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> eBasis1 = GetZbasis(EXT1_T);
  MyMatrix<T> eBasis2 = GetZbasis(EXT2_T);
  MyMatrix<T> InvBasis1 = Inverse(eBasis1);
  MyMatrix<T> InvBasis2 = Inverse(eBasis2);
  MyMatrix<T> EXTbas1 = EXT1_T * InvBasis1;
  MyMatrix<T> EXTbas2 = EXT2_T * InvBasis2;
#ifdef SANITY_CHECK_MATRIX_GROUP
  using Tidx = typename Telt::Tidx;
  for (auto &eMatGen2 : ListMatrGens2) {
    std::optional<std::vector<Tidx>> opt_eList =
        RepresentVertexPermutationTest<T, T, Tidx>(EXT2_T, EXT2_T, eMatGen2);
    if (!opt_eList) {
      std::cerr << "LinPolytopeIntegral_Isomorphism_Subspaces: We fail to "
                   "represent the matrix as a permutation of the rows\n";
      throw TerminalException{1};
    }
  }
#endif
  //
  MyMatrix<T> TheMatEquiv = FindTransformation(EXTbas1, EXTbas2, eEquiv);
  std::vector<MyMatrix<T>> ListMatrGen;
  for (auto &eGen : ListMatrGens2) {
    MyMatrix<T> NewGen = eBasis2 * eGen * InvBasis2;
    ListMatrGen.push_back(NewGen);
  }
  FiniteMatrixGroupHelper<T, Telt, TintGroup> helper =
    ComputeFiniteMatrixGroupHelper<T, Telt, TintGroup>(EXTbas2);
  MyMatrix<T> eLatt1 = Inverse(eBasis1) * TheMatEquiv;
  MyMatrix<T> eLatt2 = Inverse(eBasis2);
  std::optional<MyMatrix<T>> opt = LinearSpace_Equivalence<T, Tgroup>(
      ListMatrGen, helper, eLatt1, eLatt2, os);
  if (!opt)
    return {};
  MyMatrix<T> const &eSpaceEquiv = *opt;
  MyMatrix<T> eMatFinal = InvBasis1 * TheMatEquiv * eSpaceEquiv * eBasis2;
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!IsIntegralMatrix(eMatFinal)) {
    std::cerr << "LinPolytopeIntegral_Isomorphism_Subspaces: eMatFinal should "
                 "be integral\n";
    throw TerminalException{1};
  }
#endif
  return eMatFinal;
}

// GRP1 is a group of automorphism preserving EXT1_T
// ePerm is a transformation mapping EXT1 to EXT2.
// We are searching for a transformation h in GRP1 such that
// h * ePerm is an integral transformation.
template <typename T, typename Tgroup>
std::optional<MyMatrix<T>> LinPolytopeIntegral_Isomorphism_Method8(
    MyMatrix<T> const &EXT1_T, MyMatrix<T> const &EXT2_T, Tgroup const &GRP1,
    typename Tgroup::Telt const &ePerm, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  std::vector<MyMatrix<T>> ListMatrGens;
  std::vector<Telt> LGen = GRP1.GeneratorsOfGroup();
  for (auto &eGen : LGen) {
    Telt ePermGen = (~ePerm) * eGen * ePerm;
    MyMatrix<T> eMatr = FindTransformation(EXT2_T, EXT2_T, ePermGen);
    ListMatrGens.push_back(eMatr);
  }
  return LinPolytopeIntegral_Isomorphism_Subspaces<T, Tgroup>(
      EXT1_T, EXT2_T, ListMatrGens, ePerm, os);
}

template<typename Telt, typename T>
std::string get_matrs_as_string(MyMatrix<T> const& EXT, std::vector<Telt> const &l_elt) {
  std::string strGAPmatr = "[";
  bool IsFirst = true;
  for (auto &eElt : l_elt) {
    MyMatrix<T> M = RepresentVertexPermutation(EXT, EXT, eElt);
    if (!IsFirst)
      strGAPmatr += ",";
    IsFirst = false;
    strGAPmatr += StringMatrixGAP(M);
  }
  strGAPmatr += "]";
  return strGAPmatr;
}

// clang-format off
#endif  // SRC_GROUP_MATRIXGROUP_H_
// clang-format on

// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_INVARIANTVECTORFAMILY_H_
#define SRC_LATT_INVARIANTVECTORFAMILY_H_

// clang-format off
#include "Shvec_exact.h"
#include "Positivity.h"
#include "Timings.h"
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <utility>
#include <map>
// clang-format on

#ifdef TIMINGS
#define TIMINGS_INVARIANT_VECTOR_FAMILY
#endif

#ifdef DEBUG
#define DEBUG_INVARIANT_VECTOR_FAMILY
#endif

#ifdef DISABLE_DEBUG_INVARIANT_VECTOR_FAMILY
#undef DEBUG_INVARIANT_VECTOR_FAMILY
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_INVARIANT_VECTOR_FAMILY
#endif

/*
  We are considering here the enumeration of configurations of vectors for
  positive definite quadratic forms. There are many possible optimization
  and choices.

  The basis assumption of the code are the following:
  --- If enumerating vectors at norm N, then the cost for vectors of norms
  smaller than N is very small compared to the cost for the vectors of norm
  exactly N. Therefore, it is ok, to redo many small enumeration until we get
  the configuration of vector that we want.
  --- The increment are usually monotonous like norm 2,3,4. So, going with
  the increment given by GCD is probably a good heuristic.

 */

template <typename T, typename Tint>
T GetMaxNorm(MyMatrix<T> const &eMat, std::ostream &os) {
  LLLreduction<T, Tint> recLLL = LLLreducedBasis<T, Tint>(eMat, os);
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(recLLL.Pmat);
  MyMatrix<T> eMatRed = Pmat_T * eMat * TransposedMat(Pmat_T);
  return MaximumDiagonal(eMatRed);
}

template <typename T> T GetSmallestIncrement(MyMatrix<T> const &eMat) {
  int n = eMat.rows();
  T eGcd = eMat(0, 0);
  // No early exit at eGcd = 1: over the rationals the gcd can decrease
  // below 1 (gcd(1, 2/3) = 1/3), so 1 is not a stopping point. A rational
  // Gram matrix does occur, e.g. the realizability matrices of the
  // short-vector configurations.
  for (int i = 1; i < n; i++) {
    eGcd = GcdPair(eGcd, eMat(i, i));
  }
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      T val = 2 * eMat(i, j);
      eGcd = GcdPair(eGcd, val);
    }
  }
  return eGcd;
}

template <typename T, typename Tint>
std::set<T> GetSetNormConsider(MyMatrix<T> const &eMat, std::ostream &os) {
  T incr = GetSmallestIncrement(eMat);
  T MaxNorm = GetMaxNorm<T, Tint>(eMat, os);
  std::set<T> AllowedNorms;
  T norm = incr;
  while (true) {
    if (norm > MaxNorm)
      break;
    AllowedNorms.insert(norm);
    norm += incr;
  }
  return AllowedNorms;
}

/*
  The rank ought to be self explanatory.
  The index is the index in the lattice obtaine by the saturation
  of the spanned lattce
 */
template <typename Tint> struct FundInvariantVectorFamily {
  int rank;
  Tint index;
};

template <typename Tint>
FundInvariantVectorFamily<Tint> TrivFundamentalInvariant() {
  return {0, 1};
}

template <typename Tint>
bool IsCompleteSystem(FundInvariantVectorFamily<Tint> const &fi, int n) {
  if (fi.rank != n) {
    return false;
  }
  return fi.index == 1;
}

template <typename T> bool is_antipodal(MyMatrix<T> const &SHV) {
  int n_row = SHV.rows();
  std::unordered_set<MyVector<T>> set;
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyVector<T> V = GetMatrixRow(SHV, i_row);
    set.insert(V);
  }
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyVector<T> V = -GetMatrixRow(SHV, i_row);
    if (!set.contains(V)) {
      return false;
    }
  }
  return true;
}

template <typename T> MyMatrix<T> matrix_duplication(MyMatrix<T> const &SHV) {
  int dim = SHV.cols();
  int nbSHV = SHV.rows();
  MyMatrix<T> SHV_T(2 * nbSHV, dim);
  for (int i_row = 0; i_row < nbSHV; i_row++) {
    for (int i = 0; i < dim; i++) {
      T val = SHV(i_row, i);
      SHV_T(2 * i_row, i) = val;
      SHV_T(2 * i_row + 1, i) = -val;
    }
  }
  return SHV_T;
}

template <typename T> bool has_duplication(MyMatrix<T> const &SHV) {
  int n_row = SHV.rows();
  std::unordered_set<MyVector<T>> set;
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyVector<T> V = GetMatrixRow(SHV, i_row);
    if (!set.insert(V).second) {
      return true;
    }
  }
  return false;
}

template <typename T> void check_antipodality_mymatrix(MyMatrix<T> const &SHV) {
  if (!is_antipodal(SHV)) {
    std::cerr << "TSPACE: The family SHV is not antipodal\n";
    throw TerminalException{1};
  }
}

template <typename Tint> bool IsFullDimZbasis(MyMatrix<Tint> const &M, [[maybe_unused]] std::ostream &os) {
  int n = M.cols();
  int rnk = RankMat(M);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: IsFullDimZbasis rnk=" << rnk << " n=" << n << "\n";
#endif
  if (rnk < n) {
    return false;
  }
  Tint indx = Int_IndexLattice(M);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: IsFullDimZbasis indx=" << indx << "\n";
#endif
  if (T_NormGen(indx) == 1) {
    return true;
  }
  return false;
}

template <typename Tint> bool IsFullDim(MyMatrix<Tint> const &M, [[maybe_unused]] std::ostream &os) {
  int n = M.cols();
  int rnk = RankMat(M);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: IsFullDimZbasis rnk=" << rnk << " n=" << n << "\n";
#endif
  if (rnk < n) {
    return false;
  }
  return true;
}

template <typename Tint>
FundInvariantVectorFamily<Tint>
ComputeFundamentalInvariant(MyMatrix<Tint> const &M,
                            [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: ComputeFundamentalInvariant, beginning\n";
  WriteMatrix(os, M);
#endif
  MyVector<Tint> eVect = SmithNormalFormInvariant(M);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: ComputeFundamentalInvariant, We have eVect\n";
#endif
  int dim = eVect.size();
  int rank = 0;
  Tint index = 1;
  for (int i = 0; i < dim; i++) {
    Tint val = eVect(i);
    if (val != 0) {
      rank++;
      index *= val;
    }
  }
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
  int n = M.cols();
  if (rank != RankMat(M)) {
    std::cerr << "Something is inconsistent here\n";
    throw TerminalException{1};
  }
  if (rank == n) {
    Tint index_B = T_abs(Int_IndexLattice(M));
    if (index != index_B) {
      std::cerr << "index=" << index << " but index_B=" << index_B << "\n";
      throw TerminalException{1};
    }
  }
#endif
  return {rank, index};
}

template <typename Tint>
FundInvariantVectorFamily<Tint>
ComputeRankInvariant(MyMatrix<Tint> const &M,
                     [[maybe_unused]] std::ostream &os) {
  int rank = RankMat(M);
  Tint index(1);
  return {rank, index};
}

template <typename Tint>
bool operator>(FundInvariantVectorFamily<Tint> const &x,
               FundInvariantVectorFamily<Tint> const &y) {
  if (x.rank != y.rank) {
    return x.rank > y.rank;
  }
  if (x.index != y.index) {
    return x.index < y.index;
  }
  return false;
}

template <typename Tint>
bool operator<(FundInvariantVectorFamily<Tint> const &x,
               FundInvariantVectorFamily<Tint> const &y) {
  if (x.rank != y.rank) {
    return x.rank < y.rank;
  }
  if (x.index != y.index) {
    return x.index > y.index;
  }
  return false;
}

// Core search loop. Returns the family before antipodal duplication (one
// representative per +/-v pair); ExtractInvariantVectorFamily below
// duplicates it to restore the historical (fully duplicated) public
// return convention. Callers that only need one representative per pair
// (e.g. to feed conversion_and_duplication themselves later) should call
// this directly instead of duplicating then immediately halving again.
/*
  The shell based family, built one shell at a time so that its cost can be
  measured and compared with that of another construction while it runs.

  A shell has to be taken whole, so the family grows in jumps that can be
  enormous: on a genus of determinant 351 and rank 14 it reaches 76448 pairs
  where V_cv gives 1576. Which of the two is the cheaper is not known in
  advance -- V_cv pays 2^k closest-vector computations when the minimal
  vectors span a sublattice of index 2^k -- so neither may be run to the end
  before the other is looked at.
 */
template <typename T, typename Tint, typename Ffinal, typename Finvariant>
struct ShellFamilyBuilder {
  ShellFamilyBuilder(MyMatrix<T> const &GramMat, Ffinal f_final,
                     Finvariant f_invariant, std::ostream &os)
      : f_final_(f_final), f_invariant_(f_invariant), n_(GramMat.rows()),
        incr_(GetSmallestIncrement(GramMat)),
        max_norm_(GetMaxNorm<T, Tint>(GramMat, os)), norm_(incr_),
        SHVret_(0, GramMat.rows()),
        fi_ret_(TrivFundamentalInvariant<Tint>()), solver_(GramMat, os),
        done_(false) {}

  // Add the next shell. After it the family may be final, see is_done.
  void one_shell(std::ostream &os) {
    if (done_) {
      return;
    }
    if (norm_ > max_norm_) {
      std::cerr << "IVF: failed to find a relevant vector configuration\n";
      throw TerminalException{1};
    }
    std::vector<MyVector<Tint>> ListVect = solver_.fixed_norm_vectors(norm_);
    int n_prev = SHVret_.rows();
    int n_new = ListVect.size();
    MyMatrix<Tint> SHVtest(n_prev + n_new, n_);
    for (int i_row = 0; i_row < n_prev; i_row++) {
      for (int i_col = 0; i_col < n_; i_col++) {
        SHVtest(i_row, i_col) = SHVret_(i_row, i_col);
      }
    }
    for (int i_row = 0; i_row < n_new; i_row++) {
      for (int i_col = 0; i_col < n_; i_col++) {
        SHVtest(n_prev + i_row, i_col) = ListVect[i_row](i_col);
      }
    }
    FundInvariantVectorFamily<Tint> fi_test = f_invariant_(SHVtest);
    if (fi_ret_ < fi_test) {
      SHVret_ = SHVtest;
      fi_ret_ = fi_test;
      if (f_final_(SHVret_, fi_ret_)) {
        done_ = true;
      }
    }
    norm_ += incr_;
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
    os << "IVF: shell up to norm " << norm_ << ", " << SHVret_.rows()
       << " vectors, done=" << done_ << "\n";
#endif
  }
  bool is_done() const { return done_; }
  size_t size() const { return static_cast<size_t>(SHVret_.rows()); }
  MyMatrix<Tint> const &get() const { return SHVret_; }

private:
  Ffinal f_final_;
  Finvariant f_invariant_;
  int n_;
  T incr_, max_norm_, norm_;
  MyMatrix<Tint> SHVret_;
  FundInvariantVectorFamily<Tint> fi_ret_;
  CVPSolver<T, Tint> solver_;
  bool done_;
};

template <typename T, typename Tint, typename Ffinal, typename Finvariant>
MyMatrix<Tint> ExtractInvariantVectorFamilyHalf(MyMatrix<T> const &GramMat,
                                                Ffinal f_final,
                                                Finvariant f_invariant,
                                                std::ostream &os) {
  ShellFamilyBuilder<T, Tint, Ffinal, Finvariant> builder(GramMat, f_final,
                                                          f_invariant, os);
  while (!builder.is_done()) {
    builder.one_shell(os);
  }
  return builder.get();
}

template <typename T, typename Tint, typename Ffinal, typename Finvariant>
MyMatrix<Tint> ExtractInvariantVectorFamily(MyMatrix<T> const &GramMat,
                                            Ffinal f_final,
                                            Finvariant f_invariant,
                                            std::ostream &os) {
  MyMatrix<Tint> SHVhalf =
      ExtractInvariantVectorFamilyHalf<T, Tint, Ffinal, Finvariant>(
          GramMat, f_final, f_invariant, os);
  return matrix_duplication(SHVhalf);
}

// The half family, one vector per antipodal pair, before the duplication of
// ExtractInvariantVectorFamilyFullRank above; for callers that want the pairs.
template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyFullRankHalf(MyMatrix<T> const &eMat,
                                                        std::ostream &os) {
  int n = eMat.rows();
  auto f_final = [&]([[maybe_unused]] MyMatrix<Tint> const &M,
                     FundInvariantVectorFamily<Tint> const &fi) -> bool {
    return fi.rank == n;
  };
  auto f_invariant = [&](MyMatrix<Tint> const &M)
      -> FundInvariantVectorFamily<Tint> {
    return ComputeRankInvariant(M, os);
  };
  return ExtractInvariantVectorFamilyHalf<T, Tint, decltype(f_final),
                                          decltype(f_invariant)>(
      eMat, f_final, f_invariant, os);
}

// Representatives of the first K orbits of nonzero integer vectors under the
// group generated by ListGen (acting by v -> g v), each paired with its orbit
// size, ordered by increasing GramMat-norm v^T GramMat v. Vectors are enumerated
// shell by shell with the norm increasing by GetSmallestIncrement(GramMat), so
// there is no bound parameter. (Antipodes -v and v fall in the same orbit
// exactly when -Id is in the group, which holds e.g. for an arithmetic
// automorphism group.)
template <typename T, typename Tint>
std::vector<std::pair<MyVector<Tint>, size_t>>
get_k_short_orbit_vectors(MyMatrix<T> const &GramMat,
                          std::vector<MyMatrix<Tint>> const &ListGen, size_t K,
                          std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  T incr = GetSmallestIncrement(GramMat);
  T norm = incr;
  std::vector<std::pair<MyVector<Tint>, size_t>> reps;
  while (reps.size() < K) {
    std::vector<MyVector<Tint>> ListVect = solver.fixed_norm_vectors(norm);
    // The orbit of v under v -> g v stays within this norm shell (the ListGen
    // preserve the GramMat), so the dedup set only needs the current shell.
    // Vectors are sign-canonicalized, so antipodes -v and v count as one (the
    // orbit size is then the number of distinct lines, i.e. of distinct v v^T
    // directions).
    std::unordered_set<MyVector<Tint>> assigned;
    for (auto &v : ListVect) {
      if (reps.size() >= K) {
        break;
      }
      MyVector<Tint> cv = SignCanonicalizeVector(v);
      // insert().second is false when cv was already assigned to an earlier
      // orbit of this shell.
      if (!assigned.insert(cv).second) {
        continue;
      }
      std::vector<MyVector<Tint>> todo{cv};
      size_t orbit_size = 1;
      while (!todo.empty()) {
        MyVector<Tint> w = todo.back();
        todo.pop_back();
        for (auto &g : ListGen) {
          MyVector<Tint> gw = g * w;
          MyVector<Tint> cgw = SignCanonicalizeVector(gw);
          if (assigned.insert(cgw).second) {
            todo.push_back(cgw);
            orbit_size++;
          }
        }
      }
      reps.push_back({v, orbit_size});
    }
    norm += incr;
  }
  return reps;
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyFullRank(MyMatrix<T> const &eMat,
                                                    std::ostream &os) {
  int n = eMat.rows();
  auto f_final = [&]([[maybe_unused]] MyMatrix<Tint> const &M, FundInvariantVectorFamily<Tint> const& fi) -> bool {
    return fi.rank == n;
  };
  auto f_invariant=[&](MyMatrix<Tint> const &M) -> FundInvariantVectorFamily<Tint> {
    return ComputeRankInvariant(M, os);
  };
  return ExtractInvariantVectorFamily<T, Tint, decltype(f_final), decltype(f_invariant)>(eMat, f_final, f_invariant, os);
}

// Half version (one representative per +/-v pair, no antipodal
// duplication) of ExtractInvariantVectorFamilyZbasis below, for callers
// that would otherwise immediately undo the duplication themselves.
template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyZbasisHalf(MyMatrix<T> const &eMat,
                                                      std::ostream &os) {
  int n = eMat.cols();
  auto f_final = [&]([[maybe_unused]] MyMatrix<Tint> const &M, FundInvariantVectorFamily<Tint> const& fi) -> bool {
    return IsCompleteSystem(fi, n);
  };
  auto f_invariant=[&](MyMatrix<Tint> const &M) -> FundInvariantVectorFamily<Tint> {
    return ComputeFundamentalInvariant(M, os);
  };
  return ExtractInvariantVectorFamilyHalf<T, Tint, decltype(f_final), decltype(f_invariant)>(eMat, f_final, f_invariant, os);
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyZbasis(MyMatrix<T> const &eMat,
                                                  std::ostream &os) {
  return matrix_duplication(ExtractInvariantVectorFamilyZbasisHalf<T, Tint>(eMat, os));
}

/*
  The base block of one level of an iterated construction, both signs. With
  use_roots it is the set of roots -- the vectors of norm 2 -- when the
  lattice has any, and the shortest vectors otherwise; without it, always the
  shortest vectors. The two agree on an even lattice of minimum 2, where the
  shortest vectors are exactly the roots; they differ only when the minimum
  is below 2, in which case the root variant skips the sub-root vectors and
  starts the tower from the root sublattice, matching the root decomposition
  of LatticeRootDecomposition.h. The set is always closed under v -> -v.
 */
template <typename T, typename Tint>
MyMatrix<Tint> get_shv_block(MyMatrix<T> const &GramMat, bool use_roots,
                             std::ostream &os) {
  int n = GramMat.rows();
  if (use_roots) {
    MyMatrix<Tint> roots = T_ShortVector_fixed<T, Tint>(GramMat, T(2), os);
    if (roots.rows() > 0) {
      // Symmetrise to both signs, whatever the enumerator's convention.
      std::unordered_set<MyVector<Tint>> seen;
      std::vector<MyVector<Tint>> rows;
      for (int i = 0; i < roots.rows(); i++) {
        MyVector<Tint> v = GetMatrixRow(roots, i);
        if (seen.insert(v).second) {
          rows.push_back(v);
        }
        MyVector<Tint> mv = -v;
        if (seen.insert(mv).second) {
          rows.push_back(mv);
        }
      }
      return MatrixFromVectorFamilyDim(n, rows);
    }
  }
  return T_ShortestVector<T, Tint>(GramMat, os).SHV;
}

/*
  The iterated-shortest vector family, one vector per antipodal pair: the
  shortest vectors of the lattice, and, when they are not of full rank, the
  same construction applied to the integral orthogonal complement of their
  span, lifted back to L and concatenated. The result is of full rank but
  need not span L: it is the direct sum of the shortest-vector families of
  an orthogonal tower of sublattices, which has finite index in L. Unlike
  the norm-ball or the well-rounded families it never enumerates a shell
  above the minimum, so it stays small even when one lattice direction is
  far longer than the others -- which is exactly the case that makes the
  norm-ball family explode. A full-rank but not necessarily Z-spanning
  family, which is all the Plesken-Souvignier engine requires.

  With use_roots the tower starts from the roots rather than the shortest
  vectors (see get_shv_block); this is the "root_iterated_shortest" family.
 */
template <typename T, typename Tint>
MyMatrix<Tint> IteratedShortestVectorFamilyHalf_gen(MyMatrix<T> const &GramMat,
                                                    bool use_roots,
                                                    std::ostream &os) {
  int n = GramMat.rows();
  MyMatrix<Tint> SHV = get_shv_block<T, Tint>(GramMat, use_roots, os);
  int r = RankMat(SHV);
  // One representative per +/-v pair of the shortest vectors.
  std::unordered_set<MyVector<Tint>> seen;
  std::vector<MyVector<Tint>> half_rows;
  for (int i = 0; i < SHV.rows(); i++) {
    MyVector<Tint> c = SignCanonicalizeVector(GetMatrixRow(SHV, i));
    if (seen.insert(c).second) {
      half_rows.push_back(c);
    }
  }
  MyMatrix<Tint> SHVhalf = MatrixFromVectorFamilyDim(n, half_rows);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: iterated_shortest_half n=" << n << " |SHVhalf|="
     << SHVhalf.rows() << " r=" << r << "\n";
#endif
  if (r == n) {
    return SHVhalf;
  }
  // The integral orthogonal complement of the span of the shortest vectors,
  // as a Z-basis (rows) in the coordinates of L: the x of L with <x, s> = 0
  // for every shortest vector s, i.e. x (G B^T) = 0 over Z.
  MyMatrix<Tint> B = GetZbasis(SHV);
  MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
  MyMatrix<T> Prod_T = GramMat * B_T.transpose();
  MyMatrix<Tint> Prod = UniversalMatrixConversion<Tint, T>(Prod_T);
  MyMatrix<Tint> Perp = NullspaceIntMat(Prod);
  MyMatrix<T> Perp_T = UniversalMatrixConversion<T, Tint>(Perp);
  MyMatrix<T> PerpGram = Perp_T * GramMat * Perp_T.transpose();
  // Recurse on the complement, then lift: a row w in Perp coordinates is
  // the vector w * Perp of L. The lift keeps one representative per pair,
  // and the two blocks lie in orthogonal subspaces, so the concatenation
  // is still one vector per antipodal pair.
  MyMatrix<Tint> PerpFamHalf =
      IteratedShortestVectorFamilyHalf_gen<T, Tint>(PerpGram, use_roots, os);
  MyMatrix<Tint> PerpFamHalf_L = PerpFamHalf * Perp;
  MyMatrix<Tint> result = Concatenate(SHVhalf, PerpFamHalf_L);
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
  if (RankMat(result) != n) {
    std::cerr << "IVF: iterated_shortest did not reach full rank\n";
    throw TerminalException{1};
  }
#endif
  return result;
}

// The iterated-shortest family, one vector per antipodal pair.
template <typename T, typename Tint>
MyMatrix<Tint> IteratedShortestVectorFamilyHalf(MyMatrix<T> const &GramMat,
                                                std::ostream &os) {
  return IteratedShortestVectorFamilyHalf_gen<T, Tint>(GramMat, false, os);
}

// The root-iterated-shortest family (roots first, then the iterated-shortest
// tower on the orthogonal complement), one vector per antipodal pair.
template <typename T, typename Tint>
MyMatrix<Tint>
RootIteratedShortestVectorFamilyHalf(MyMatrix<T> const &GramMat,
                                     std::ostream &os) {
  return IteratedShortestVectorFamilyHalf_gen<T, Tint>(GramMat, true, os);
}

// The full (both signs) root-iterated-shortest family.
template <typename T, typename Tint>
MyMatrix<Tint> root_iterated_shortest(MyMatrix<T> const &GramMat,
                                      std::ostream &os) {
  return matrix_duplication(
      RootIteratedShortestVectorFamilyHalf<T, Tint>(GramMat, os));
}

// The full (both signs) iterated-shortest family.
template <typename T, typename Tint>
MyMatrix<Tint> IteratedShortestVectorFamily(MyMatrix<T> const &GramMat,
                                            std::ostream &os) {
  return matrix_duplication(IteratedShortestVectorFamilyHalf<T, Tint>(GramMat, os));
}

/*
  Complete a full-rank vector family to a Z-spanning one, generically -- it
  does not matter how the family was built. The family spans a finite-index
  sublattice K = <V> of L = Z^n; for every nonzero coset of Z^n / K this
  adds ALL its minimal-norm representatives (a closest-vector problem per
  coset: the representatives are the alpha - k for k a closest point of K to
  the coset rep alpha). Adding every minimal representative, not one, keeps
  the result an invariant of the lattice -- an isometry permutes the cosets
  and preserves norms, so it maps a coset minimum-set onto a coset
  minimum-set, but not a single chosen representative onto a chosen one.

  This is the saturation completion of Hecke's characteristic vectors
  (_characteristic_vectors) specialised to an already-full-rank family: the
  saturation of a full-rank sublattice is the whole lattice, so there is no
  recursion, only the one completion over Z^n / K. The result Z-spans, so a
  searched basis drawn from it is unimodular and the Plesken-Souvignier
  backtrack has no rational leaves.
 */
template <typename T, typename Tint>
MyMatrix<Tint> ivf_z_spanning(MyMatrix<T> const &GramMat,
                              MyMatrix<Tint> const &V, std::ostream &os) {
  int n = GramMat.rows();
  MyMatrix<Tint> B0 = GetZbasis(V);
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
  if (B0.rows() != n) {
    std::cerr << "IVF: ivf_z_spanning requires a full-rank family\n";
    throw TerminalException{1};
  }
#endif
  std::unordered_set<MyVector<Tint>> present;
  std::vector<MyVector<Tint>> rows;
  for (int i = 0; i < V.rows(); i++) {
    MyVector<Tint> v = GetMatrixRow(V, i);
    if (present.insert(v).second) {
      rows.push_back(v);
    }
  }
  // K = <V> already equals Z^n: nothing to complete.
  if (T_abs(DeterminantMatBareiss(B0)) == Tint(1)) {
    return MatrixFromVectorFamilyDim(n, rows);
  }
  // A reduced basis of K. The closest-vector problems are on the form that
  // K carries in its own basis, and a skewed basis (which GetZbasis of the
  // iterated-shortest family typically is) makes them very slow -- so
  // LLL-reduce it first, as Hecke does before its coset CVPs.
  MyMatrix<T> B0_T = UniversalMatrixConversion<T, Tint>(B0);
  MyMatrix<T> G_K0 = B0_T * GramMat * B0_T.transpose();
  LLLreduction<T, Tint> lll = LLLreducedBasis<T, Tint>(G_K0, os);
  MyMatrix<Tint> B = lll.Pmat * B0;
  MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
  MyMatrix<T> Binv = Inverse(B_T);
  MyMatrix<T> G_K = lll.GramMatRed;
  CVPSolver<T, Tint> solver(G_K, os);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  MicrosecondTime time_tc;
#endif
  std::vector<MyVector<Tint>> cosets = ComputeTranslationClasses<Tint, Tint>(B);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: ivf_z_spanning index " << cosets.size() << " translation classes "
     << "in " << time_tc << "\n";
  MicrosecondTime time_cvp;
#endif
  for (auto &alpha : cosets) {
    if (IsZeroVector(alpha)) {
      continue;
    }
    // alpha in the coordinates of the basis B: a = alpha B^{-1}.
    MyVector<T> a(n);
    for (int j = 0; j < n; j++) {
      T s(0);
      for (int i = 0; i < n; i++) {
        s += UniversalScalarConversion<T, Tint>(alpha(i)) * Binv(i, j);
      }
      a(j) = s;
    }
    resultCVP<T, Tint> cvp = solver.nearest_vectors(a);
    for (int r = 0; r < cvp.ListVect.rows(); r++) {
      MyVector<Tint> c = GetMatrixRow(cvp.ListVect, r);
      // w = alpha - c B, a vector of Z^n minimal in the coset alpha + K.
      MyVector<Tint> w(n);
      for (int k = 0; k < n; k++) {
        Tint s = alpha(k);
        for (int i = 0; i < n; i++) {
          s -= c(i) * B(i, k);
        }
        w(k) = s;
      }
      if (present.insert(w).second) {
        rows.push_back(w);
      }
    }
  }
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: ivf_z_spanning " << (cosets.size() - 1) << " CVPs in " << time_cvp
     << ", family " << V.rows() << " -> " << rows.size() << "\n";
#endif
  return MatrixFromVectorFamilyDim(n, rows);
}

/*
  The Z-spanning characteristic vector family of a Gram matrix, the port of
  Hecke's _characteristic_vectors (following Sikiric-Haensch-Voight-van
  Woerden). It returns the full set (both signs) of vectors in the
  coordinates of L = Z^n. Unlike the generic ivf_z_spanning, which completes
  an already-flat full-rank family in one shot and pays for a closest vector
  problem on the whole mixed-scale sublattice, this builds the family through
  the same orthogonal recursion as the iterated-shortest family, so every
  closest vector problem lives on the shortest-vector sublattice S1 or its
  saturation P1 -- lattices whose vectors are all near the minimum, hence
  well conditioned -- and never on a lattice mixing a long direction with the
  short ones. That is what keeps the completion tractable where the one-shot
  completion explodes.

  At each level, with SHV the base block (get_shv_block: the shortest vectors,
  or the roots when use_roots is set), S1 = <SHV> and P1 its saturation in
  Z^n:
   - the block vectors are characteristic vectors;
   - for every non-zero coset of P1 / S1, the minimal vectors of that coset
     (a closest vector problem in S1) are characteristic vectors -- this
     completes the span of S1 up to its saturation P1;
   - if S1 is not of full rank, we recurse on the PROJECTION L2 of Z^n onto
     the orthogonal complement of P1 (not the integral intersection, which
     would drop the glue between the blocks). A characteristic vector a of L2
     that lifts integrally into Z^n is one directly; otherwise its minimal
     integral lifts are a + (minimal P1-glue), a closest vector problem in P1.
  The result Z-spans L. With use_roots this is the Z-spanning
  root-iterated-shortest family (Hecke's _reduced_characteristic_vectors).
 */
template <typename T, typename Tint>
std::vector<MyVector<Tint>>
inner_span_iterated_shortest_rec(MyMatrix<T> const &G, bool use_roots,
                                 std::ostream &os) {
  int n = G.rows();
  MyMatrix<Tint> SHV = get_shv_block<T, Tint>(G, use_roots, os);
  int r = RankMat(SHV);
  std::vector<MyVector<Tint>> cvL;
  for (int i = 0; i < SHV.rows(); i++) {
    cvL.push_back(GetMatrixRow(SHV, i));
  }
  // S1 = <SHV>, P1 = saturation of S1 in Z^n, both as row bases in Z^n.
  MyMatrix<Tint> S1 = GetZbasis(SHV);
  MyMatrix<Tint> P1 = IntegralSpaceSaturation(S1);
  // C: the coordinates of S1 in P1, i.e. the integer matrix with S1 = C P1.
  // |det C| = [P1 : S1] is the index whose cosets the completion runs over.
  MyMatrix<Tint> C(r, r);
  for (int i = 0; i < r; i++) {
    MyVector<Tint> s = GetMatrixRow(S1, i);
    std::optional<MyVector<Tint>> opt = SolutionIntMat(P1, s);
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
    if (!opt) {
      std::cerr << "IVF: inner_span, a shortest-vector basis vector is not in "
                << "the saturation\n";
      throw TerminalException{1};
    }
#endif
    for (int j = 0; j < r; j++) {
      C(i, j) = (*opt)(j);
    }
  }
  MyMatrix<T> S1_T = UniversalMatrixConversion<T, Tint>(S1);
  // Saturation completion S1 -> P1: only needed when S1 is a proper
  // sublattice of its saturation.
  if (T_abs(DeterminantMatBareiss(C)) != Tint(1)) {
    MyMatrix<T> G_S1 = S1_T * G * S1_T.transpose();
    CVPSolver<T, Tint> solverS1(G_S1, os);
    std::vector<MyVector<Tint>> cosets =
        ComputeTranslationClasses<Tint, Tint>(C);
    for (auto &cr : cosets) {
      if (IsZeroVector(cr)) {
        continue;
      }
      // p = cr P1, a representative of the coset in Z^n, and a_S1 its
      // coordinates in the basis of S1 (rational, since p is not in S1).
      MyVector<Tint> p = P1.transpose() * cr;
      MyVector<T> p_T = UniversalVectorConversion<T, Tint>(p);
      std::optional<MyVector<T>> a_opt = SolutionMat(S1_T, p_T);
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
      if (!a_opt) {
        std::cerr << "IVF: inner_span, a coset rep is not in the span of S1\n";
        throw TerminalException{1};
      }
#endif
      resultCVP<T, Tint> cvp = solverS1.nearest_vectors(*a_opt);
      for (int rr = 0; rr < cvp.ListVect.rows(); rr++) {
        MyVector<Tint> j = GetMatrixRow(cvp.ListVect, rr);
        // w = p - j S1, a minimal vector of the coset p + S1 in Z^n.
        MyVector<Tint> w = p - S1.transpose() * j;
        cvL.push_back(w);
      }
    }
  }
  if (r == n) {
    return cvL;
  }
  // The G-orthogonal projection onto the span of P1 and its complement:
  // pr1 = G P1^T (P1 G P1^T)^{-1} P1, proj2 = I - pr1.
  MyMatrix<T> P1_T = UniversalMatrixConversion<T, Tint>(P1);
  MyMatrix<T> GP1t = G * P1_T.transpose();
  MyMatrix<T> P1GP1t = P1_T * GP1t;
  MyMatrix<T> P1GP1t_inv = Inverse(P1GP1t);
  MyMatrix<T> pr1 = GP1t * P1GP1t_inv * P1_T;
  MyMatrix<T> proj2 = IdentityMat<T>(n) - pr1;
  // L2 = projection of Z^n onto the complement, as a rational row basis:
  // clear denominators, take an integer Z-basis, divide back.
  FractionMatrix<T> fr = RemoveFractionMatrixPlusCoeff(proj2);
  MyMatrix<Tint> proj2_scaled = UniversalMatrixConversion<Tint, T>(fr.TheMat);
  MyMatrix<Tint> B2int = GetZbasis(proj2_scaled);
  MyMatrix<T> L2basis = UniversalMatrixConversion<T, Tint>(B2int) / fr.TheMult;
  MyMatrix<T> G_L2 = L2basis * G * L2basis.transpose();
  // P_Z: coordinates of proj2(e_i) in the basis of L2, an integer n x (n-r)
  // matrix mapping x in Z^n to its projection in L2 coordinates.
  int nc = L2basis.rows();
  MyMatrix<Tint> P_Z(n, nc);
  for (int i = 0; i < n; i++) {
    MyVector<T> row = GetMatrixRow(proj2, i);
    std::optional<MyVector<T>> opt = SolutionMat(L2basis, row);
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
    if (!opt || !IsIntegralVector(*opt)) {
      std::cerr << "IVF: inner_span, a projection is not in the projection "
                << "lattice\n";
      throw TerminalException{1};
    }
#endif
    for (int j = 0; j < nc; j++) {
      P_Z(i, j) = UniversalScalarConversion<Tint, T>((*opt)(j));
    }
  }
  CVPSolver<T, Tint> solverP1(P1GP1t, os);
  std::vector<MyVector<Tint>> rec2 =
      inner_span_iterated_shortest_rec<T, Tint>(G_L2, use_roots, os);
  for (auto &a : rec2) {
    MyVector<T> a_T = UniversalVectorConversion<T, Tint>(a);
    MyVector<T> aL = L2basis.transpose() * a_T;
    if (IsIntegralVector(aL)) {
      // The lift already lies in Z^n: it is a characteristic vector as is.
      cvL.push_back(UniversalVectorConversion<Tint, T>(aL));
      continue;
    }
    // vL: an actual element of Z^n projecting to a; its minimal-norm coset
    // representatives modulo P1 are the lifts we add.
    std::optional<MyVector<Tint>> vL_opt = SolutionIntMat(P_Z, a);
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
    if (!vL_opt) {
      std::cerr << "IVF: inner_span, no integral preimage of a projection "
                << "characteristic vector\n";
      throw TerminalException{1};
    }
#endif
    MyVector<Tint> const &vL = *vL_opt;
    MyVector<T> vL_T = UniversalVectorConversion<T, Tint>(vL);
    MyVector<T> w_amb = pr1.transpose() * vL_T;
    std::optional<MyVector<T>> w_opt = SolutionMat(P1_T, w_amb);
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
    if (!w_opt) {
      std::cerr << "IVF: inner_span, the P1 part is not in the span of P1\n";
      throw TerminalException{1};
    }
#endif
    resultCVP<T, Tint> cvp = solverP1.nearest_vectors(*w_opt);
    for (int rr = 0; rr < cvp.ListVect.rows(); rr++) {
      MyVector<Tint> jp = GetMatrixRow(cvp.ListVect, rr);
      // vL - jp P1, a minimal integral lift of a into Z^n.
      MyVector<Tint> glue = vL - P1.transpose() * jp;
      cvL.push_back(glue);
    }
  }
  return cvL;
}

// The Z-spanning characteristic vector family as a matrix, deduplicated.
// use_roots selects the root variant (see get_shv_block).
template <typename T, typename Tint>
MyMatrix<Tint> inner_span_iterated_shortest_gen(MyMatrix<T> const &GramMat,
                                                bool use_roots,
                                                std::ostream &os) {
  int n = GramMat.rows();
  std::vector<MyVector<Tint>> cvL =
      inner_span_iterated_shortest_rec<T, Tint>(GramMat, use_roots, os);
  std::unordered_set<MyVector<Tint>> seen;
  std::vector<MyVector<Tint>> rows;
  for (auto &v : cvL) {
    if (seen.insert(v).second) {
      rows.push_back(v);
    }
  }
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: inner_span_iterated_shortest use_roots=" << use_roots << " n=" << n
     << " |raw|=" << cvL.size() << " |dedup|=" << rows.size() << "\n";
#endif
  return MatrixFromVectorFamilyDim(n, rows);
}

// The Z-spanning iterated-shortest characteristic vector family.
template <typename T, typename Tint>
MyMatrix<Tint> inner_span_iterated_shortest(MyMatrix<T> const &GramMat,
                                            std::ostream &os) {
  return inner_span_iterated_shortest_gen<T, Tint>(GramMat, false, os);
}

// The Z-spanning root-iterated-shortest characteristic vector family.
template <typename T, typename Tint>
MyMatrix<Tint> inner_span_root_iterated_shortest(MyMatrix<T> const &GramMat,
                                                 std::ostream &os) {
  return inner_span_iterated_shortest_gen<T, Tint>(GramMat, true, os);
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantBreakingVectorFamily(
    MyMatrix<T> const &eMat, std::vector<MyMatrix<Tint>> const &ListMatr,
    std::ostream &os) {
  auto f_check = [&](MyMatrix<Tint> const &M) -> bool {
    int n_row = M.rows();
    std::unordered_set<MyVector<Tint>> set;
    for (int i = 0; i < n_row; i++) {
      MyVector<Tint> V = GetMatrixRow(M, i);
      set.insert(V);
    }
    for (auto &eMatr : ListMatr) {
      MyMatrix<Tint> Mprod = M * eMatr;
      for (int i = 0; i < n_row; i++) {
        MyVector<Tint> V = GetMatrixRow(Mprod, i);
        if (!set.contains(V))
          return true;
      }
    }
    return false;
  };
  auto f_invariant = [&](MyMatrix<Tint> const &M) -> FundInvariantVectorFamily<Tint> {
    int artificial_rank = 0;
    if (f_check(M)) {
      artificial_rank = 1;
    }
    Tint index(1);
    return {artificial_rank, index};
  };
  auto f_final = [&]([[maybe_unused]] MyMatrix<Tint> const &M, FundInvariantVectorFamily<Tint> const& fi) -> bool {
    return fi.rank == 1;
  };
  return ExtractInvariantVectorFamily<T, Tint, decltype(f_final), decltype(f_invariant)>(eMat, f_final, f_invariant, os);
}

template <typename Tint> bool CheckCentralSymmetry(MyMatrix<Tint> const &M) {
  int nbRow = M.rows();
  std::unordered_map<MyVector<Tint>, int> map;
  for (int i = 0; i < nbRow; i++) {
    MyVector<Tint> V = GetMatrixRow(M, i);
    if (!IsZeroVector(V)) {
      MyVector<Tint> Vcan = SignCanonicalizeVector(V);
      map[Vcan]++;
    }
  }
  for (auto &[vcan, multiplicity] : map) {
    if (multiplicity != 2)
      return false;
  }
  return true;
}

template <typename T, typename Tint>
MyMatrix<Tint> ComputeVoronoiRelevantVector(MyMatrix<T> const &GramMat,
                                            std::ostream &os) {
  int n = GramMat.rows();
  std::vector<MyVector<Tint>> ListVect;
  BlockCppIterator blk(n, 2);
  CVPSolver<T, Tint> solver(GramMat, os);
  for (auto &eVect : blk) {
    int sum = 0;
    for (auto &eVal : eVect) {
      sum += eVal;
    }
    if (sum > 0) {
      MyVector<T> eV(n);
      for (int u = 0; u < n; u++) {
        T val = UniversalScalarConversion<T, int>(eVect[u]);
        eV(u) = val / T(2);
      }
      resultCVP<T, Tint> result = solver.nearest_vectors(eV);
      if (result.ListVect.rows() == 2) {
        MyVector<Tint> Vins(n);
        for (int u = 0; u < n; u++) {
          Vins(u) = result.ListVect(0, u) - result.ListVect(1, u);
        }
        ListVect.push_back(Vins);
        ListVect.push_back(-Vins);
      }
    }
  }
  return MatrixFromVectorFamily(ListVect);
}

template <typename T, typename Tint>
MyMatrix<Tint> FilterByNorm(MyMatrix<T> const &GramMat,
                            MyMatrix<Tint> const &ListVect, std::ostream &os) {
  int n = GramMat.rows();
  std::map<T, std::vector<MyVector<Tint>>> map;
  std::vector<T> LineMat = GetLineVector(GramMat);
  int n_vect = ListVect.rows();
  for (int i_vect = 0; i_vect < n_vect; i_vect++) {
    MyVector<Tint> V = GetMatrixRow(ListVect, i_vect);
    MyVector<T> V_T = UniversalVectorConversion<T, Tint>(V);
    T norm = EvaluateLineVector(LineMat, V_T);
    map[norm].push_back(V);
  }
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: FilterByNorm, map built\n";
#endif
  MyMatrix<Tint> SHV_ret(0, n);
  FundInvariantVectorFamily<Tint> fi_ret = TrivFundamentalInvariant<Tint>();
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  size_t pos = 0;
#endif
  for (auto &[norm, vectors] : map) {
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
    os << "IVF: FilterByNorm, pos=" << pos
       << " |vectors|=" << vectors.size() << "\n";
#endif
    MyMatrix<Tint> BlkMat = MatrixFromVectorFamily(vectors);
    MyMatrix<Tint> SHV_new = Concatenate(SHV_ret, BlkMat);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
    os << "IVF: FilterByNorm, We have SHV_new\n";
#endif
    FundInvariantVectorFamily<Tint> fi_new =
        ComputeFundamentalInvariant(SHV_new, os);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
    os << "IVF: FilterByNorm, We have fi_new\n";
#endif
    if (fi_new > fi_ret) {
      SHV_ret = SHV_new;
      fi_ret = fi_new;
      if (IsCompleteSystem(fi_ret, n)) {
        return SHV_ret;
      }
    }
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
    pos += 1;
#endif
  }
  std::cerr
      << "Failed to terminate and find a correct family for the filtration\n";
  std::cerr
      << "This indicates that the original family was not a spanning one\n";
  throw TerminalException{1};
}

/*
  The characteristic vector set V_cv of Section 2.2 of the paper
  "A canonical form for positive definite matrices" (Dutour Sikiric,
  Haensch, Voight, van Woerden).

  A characteristic vector set function (Definition 1.2.1) is a map A -> V(A)
  from positive definite n x n matrices to finite subsets of Z^n such that
  --- V(A) generates Z^n as a group.
  --- U^{-1} V(A) = V(U^T A U) for every U in GL_n(Z).
  The second property is what makes it usable for canonicalization: the
  canonical form of A can be computed from V(A) alone.

  The construction of V_cv is inductive.  Write Min(A) for the set of
  vectors realizing the minimum of A and cvd(A,v) = min_{x in Z^n} A[x - v]
  for the squared distance from v to the lattice, and
      CV(A, v) = { x in Z^n : A[x - v] = cvd(A, v) }
  for the set of closest vectors.  Both are canonical: CV(U^T A U, U^{-1} v)
  = U^{-1} CV(A, v), which is what makes the whole construction work.

  Step 1: the well-rounded case (2.2.5).  If A is well-rounded, that is if
  Min(A) spans R^n, let L_min be the lattice generated by Min(A) and B a
  Z-basis of it.  Then
      V_wr-cv(A) = Min(A) union
                   union_{c in Z^n / L_min} (c - B CV(B^T A B, B^{-1} c))
  The deep-hole-like vectors c - B x pick up the part of Z^n that Min(A)
  misses, so V_wr-cv(A) generates Z^n even when Min(A) does not.

  Step 2: the general case (2.2.6).  Let L_1 = satspan(Min(A)) = QMin(A) cap
  Z^n with Z-basis B_1 (r x n), and A_1 = B_1 A B_1^T.  By construction A_1
  is well-rounded.  Let proj be the orthogonal projection onto the A-
  orthogonal complement of L_1, L_2 = proj(Z^n) with basis B_2
  ((n-r) x n, rational), and A_2 = B_2 A B_2^T.  Then
      V_cv(A) = B_1 V_wr-cv(A_1) union
                union_{v in B_2 V_cv(A_2)} CV(A, v)
  with V_cv of the empty matrix set to the empty set.  Since Z^n / L_1 is
  isomorphic to L_2 by proj, the closest vectors to the lifts of a
  generating set of L_2 complete B_1 V_wr-cv(A_1) to a generating set of
  Z^n.

  Both functions take a `spanning` flag. With spanning = true the family
  generates Z^n, as Definition 1.2.1 asks. With spanning = false the coset
  unions of (2.2.5) are dropped at every level of the recursion and the
  family is only of full rank, generating a finite index subgroup of Z^n.
  The second is usually much smaller and is the better choice whenever the
  consumer can work from a full rank family, since those coset unions are
  the expensive part of the construction: one closest-vector computation per
  coset, and there are [Z^n : L_min] of them at each level. What remains is
  still canonical, the decision to stop at full rank being uniform rather
  than read off the data.

  Two deviations from the text of the paper, both needed for correctness:
  --- The coset c = 0 in (2.2.5) contributes 0 - B CV(A_1, 0) = {0}.  The
      zero vector is dropped, it carries no information and a vector family
      containing it is not usable downstream.
  --- CV(A, v) in (2.2.6) is not guaranteed to consist of vectors x with
      proj(x) = v: a lattice point in a neighbouring L_1-coset can be
      closer to v than every point of the coset of v.  When that happens
      the union misses the coset of v altogether and V_cv(A) does not
      generate Z^n, so (2.2.6) as written is not a characteristic vector
      set function.  A family of counterexamples: for k >= 3 let A(k, N) be
      the Gram matrix of e_1, ..., e_k, f with the e_i pairwise orthogonal
      of norm 4, <f, e_i> = 2 and <f, f> = N.  Then Min(A) = {+- e_i} spans
      only L_1 = <e_1, ..., e_k>, the covering radius squared of L_1 is k,
      and proj(f) has norm N - k.  For N - k < k the origin is strictly
      closer to v = proj(f) than every point of f + L_1, so CV(A, v) = {0}
      and the formula of the paper returns a family generating a subgroup
      of index 2; A(3,5), A(8,11) and A(8,15) are worked out in
      results/charvect_cv.md of the strongly perfect lattices notes, where
      the measurement is recorded.  We therefore add, for each v, the
      closest vectors *within the coset determined by v*, which is the set
      the generation argument actually needs.  This is canonical for the
      same reason CV is, the coset being read off from v.  Both sets are
      kept so that no information is lost.

  Everything is in the row convention: vectors are rows, a basis matrix B
  has the basis vectors as rows, and the Gram matrix of the sublattice
  spanned by the rows of B is B A B^T.  A row vector y in the B-coordinates
  corresponds to the ambient vector y B, that is B^T y in column form.
 */

// Min(A), with both signs, as rows of the returned matrix.
template <typename T, typename Tint>
MyMatrix<Tint> CharVectSet_MinimalVectors(MyMatrix<T> const &GramMat,
                                          std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  Tshortest<T, Tint> RecSHV = solver.shortest_vectors();
  // Already antipodal, unlike fixed_norm_vectors which returns one vector per
  // pair and is the one the shell based families have to duplicate.
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
  if (!CheckCentralSymmetry(RecSHV.SHV)) {
    std::cerr << "IVF: shortest_vectors did not return an antipodal family\n";
    throw TerminalException{1};
  }
#endif
  return RecSHV.SHV;
}

// Insertion into a list of rows with duplicate removal, the sets built below
// overlap a lot (CV(A,v) for nearby v, Min(A) and the c = 0 coset, ...).
template <typename Tint> struct CharVectSet_Accumulator {
  std::set<MyVector<Tint>> TheSet;
  // When set, the sign of every vector is normalized on its first nonzero
  // coordinate before insertion, so the set holds one vector per antipodal
  // pair and never the pair itself. That is what the consumers want: the
  // weight matrix and the graph are quadratic in the number of vectors, so
  // building the family already halved is worth more than halving it later.
  bool antipodal_half = false;
  void insert(MyVector<Tint> const &V) {
    if (IsZeroVector(V)) {
      return;
    }
    if (antipodal_half) {
      int n = V.size();
      for (int i = 0; i < n; i++) {
        if (V(i) != 0) {
          if (V(i) < 0) {
            TheSet.insert(MyVector<Tint>(-V));
            return;
          }
          break;
        }
      }
    }
    TheSet.insert(V);
  }
  void insert_rows(MyMatrix<Tint> const &M) {
    for (int i_row = 0; i_row < M.rows(); i_row++) {
      insert(GetMatrixRow(M, i_row));
    }
  }
  size_t size() const { return TheSet.size(); }
  MyMatrix<Tint> get_matrix(int const &n) const {
    MyMatrix<Tint> M(TheSet.size(), n);
    int pos = 0;
    for (auto &V : TheSet) {
      for (int i = 0; i < n; i++) {
        M(pos, i) = V(i);
      }
      pos++;
    }
    return M;
  }
};

/*
  The two ways the V_cv construction can be cut short when it is raced
  against the shells. Nothing is returned in either case. A zero means no
  limit.

  budget_ns is a time limit, used only to decide which of the two
  constructions to advance next. It must never decide WHICH family is kept:
  the canonical form depends on the family, so a choice made on timing would
  give two isometric lattices different canonical forms from one run to the
  next, and the genus enumeration would count them as distinct classes and
  overshoot its mass.

  size_bound is what the choice is made on. It is a property of the lattice,
  so a decision taken with it is the same in every run.

  The bound is checked against the accumulator of every level of the
  recursion, not only the outermost one, and it has to abort exactly when
  the final size of the outermost one would pass it, or the two ways of
  entering the race would not agree. They do, because the size of a level
  never exceeds the size of the level above it:

  - the level above inserts every row of the V_wr-cv it received, and B_1 is
    injective, so distinct pairs stay distinct;
  - it inserts, for each row y of the V_cv of the level below, the vectors
    of the coset of xLift, whose coordinates along C are exactly y, so
    distinct rows again give distinct vectors.

  An accumulator only grows, so a level passing the bound at any moment has
  a final size past it, hence so has the outermost one; and conversely an
  outermost size within the bound keeps every level within it.

  The construction is restarted rather than resumed when the budget grows.
  That wastes at most the work of the previous, shorter, attempt, so the
  total stays within twice what the winning attempt costs.
 */
struct CharVectSet_Budget {
  int64_t budget_ns;
  size_t size_bound;
  NanosecondTime time;
  CharVectSet_Budget() : budget_ns(0), size_bound(0), time() {}
  CharVectSet_Budget(int64_t const &ns, size_t const &sz)
      : budget_ns(ns), size_bound(sz), time() {}
  bool expired() const {
    return budget_ns > 0 && time.const_eval_int64() > budget_ns;
  }
  bool too_big(size_t const &siz) const {
    return size_bound > 0 && siz > size_bound;
  }
};

/*
  V_wr-cv(A) of (2.2.5), for a well-rounded A.

  With spanning = false only Min(A) is returned. That is enough to be of full
  rank, A being well-rounded, and it drops the union over Z^n / L_min, which
  is the only expensive part: it costs one closest-vector computation per
  coset and there are [Z^n : L_min] of them. What is lost is that the family
  then generates L_min instead of Z^n.
 */
template <typename T, typename Tint>
std::optional<MyMatrix<Tint>>
CharacteristicVectorSetWellRoundedCV(MyMatrix<T> const &GramMat,
                                     bool const &spanning,
                                     bool const &antipodal_half,
                                     CharVectSet_Budget const &budget,
                                     std::ostream &os) {
#ifdef TIMINGS_INVARIANT_VECTOR_FAMILY
  MicrosecondTime time;
#endif
  int n = GramMat.rows();
  MyMatrix<Tint> SHV = CharVectSet_MinimalVectors<T, Tint>(GramMat, os);
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
  if (RankMat(SHV) != n) {
    std::cerr << "IVF: CharacteristicVectorSetWellRoundedCV called on a "
              << "matrix that is not well-rounded\n";
    throw TerminalException{1};
  }
#endif
  CharVectSet_Accumulator<Tint> acc;
  acc.antipodal_half = antipodal_half;
  acc.insert_rows(SHV);
  if (!spanning) {
    // Min(A) alone, which is of full rank since A is well-rounded. The union
    // over Z^n / L_min below is what makes the family span Z^n and it is also
    // the whole cost of this function.
#ifdef TIMINGS_INVARIANT_VECTOR_FAMILY
    os << "|IVF: CharacteristicVectorSetWellRoundedCV(fullrank)|=" << time
       << "\n";
#endif
    return acc.get_matrix(n);
  }
  // B is a Z-basis of the lattice L_min generated by Min(A).
  MyMatrix<Tint> B = GetZbasis(SHV);
  MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
  MyMatrix<T> GramB = B_T * GramMat * B_T.transpose();
  // c = c_B B in row convention, that is c = B^T c_B in column form.
  MyMatrix<T> B_Ttr = B_T.transpose();
  MyMatrix<T> BinvTr = Inverse(B_Ttr);
  std::vector<MyVector<Tint>> ListClasses =
      ComputeTranslationClasses<Tint, Tint>(B);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: WellRoundedCV n=" << n << " |Min|=" << SHV.rows()
     << " |Z^n/L_min|=" << ListClasses.size() << "\n";
#endif
  CVPSolver<T, Tint> solver(GramB, os);
  for (auto &eClass : ListClasses) {
    // One closest-vector computation per coset, and there are
    // [Z^n : L_min] of them, so this is where the limits have to be honoured.
    if (budget.expired() || budget.too_big(acc.size())) {
      return {};
    }
    MyVector<T> eClass_T = UniversalVectorConversion<T, Tint>(eClass);
    MyVector<T> eClassB = BinvTr * eClass_T;
    resultCVP<T, Tint> res = solver.nearest_vectors(eClassB);
    for (int i_row = 0; i_row < res.ListVect.rows(); i_row++) {
      MyVector<Tint> y = GetMatrixRow(res.ListVect, i_row);
      // The vector is c - y B.
      MyVector<Tint> eVect = eClass - B.transpose() * y;
      acc.insert(eVect);
    }
  }
  MyMatrix<Tint> RetMat = acc.get_matrix(n);
#ifdef TIMINGS_INVARIANT_VECTOR_FAMILY
  os << "|IVF: CharacteristicVectorSetWellRoundedCV|=" << time << "\n";
#endif
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: WellRoundedCV |V_wr-cv|=" << RetMat.rows() << "\n";
#endif
  return RetMat;
}

/*
  V_cv(A) of (2.2.6).  Recursive over the successive saturated spans of the
  minimal vectors.
 */
template <typename T, typename Tint>
std::optional<MyMatrix<Tint>>
CharacteristicVectorSetCV(MyMatrix<T> const &GramMat, bool const &spanning,
                                         bool const &antipodal_half,
                                         CharVectSet_Budget const &budget,
                                         std::ostream &os) {
#ifdef TIMINGS_INVARIANT_VECTOR_FAMILY
  MicrosecondTime time;
#endif
  int n = GramMat.rows();
  MyMatrix<Tint> SHV = CharVectSet_MinimalVectors<T, Tint>(GramMat, os);
  int r = RankMat(SHV);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: CV n=" << n << " |Min|=" << SHV.rows() << " r=" << r << "\n";
#endif
  if (r == n) {
    // A is well-rounded, the recursion stops.
    return CharacteristicVectorSetWellRoundedCV<T, Tint>(
        GramMat, spanning, antipodal_half, budget, os);
  }
  CharVectSet_Accumulator<Tint> acc;
  acc.antipodal_half = antipodal_half;
  // B_1 is a Z-basis of L_1 = satspan(Min(A)), r x n.
  MyMatrix<Tint> B1 = IntegralSpaceSaturation(GetZbasis(SHV));
  MyMatrix<T> B1_T = UniversalMatrixConversion<T, Tint>(B1);
  MyMatrix<T> Gram1 = B1_T * GramMat * B1_T.transpose();
  std::optional<MyMatrix<Tint>> opt_wr =
      CharacteristicVectorSetWellRoundedCV<T, Tint>(Gram1, spanning,
                                                    antipodal_half, budget, os);
  if (!opt_wr) {
    return {};
  }
  // B_1 V_wr-cv(A_1), that is the rows y of Vwr sent to y B_1.
  acc.insert_rows((*opt_wr) * B1);
  // C completes B_1 to a Z-basis of Z^n, and B_2 = proj(C) is a basis of
  // L_2 = proj(Z^n) since proj has kernel L_1 on Z^n, L_1 being saturated.
  MyMatrix<Tint> C = SubspaceCompletionInt(B1, n);
  MyMatrix<T> C_T = UniversalMatrixConversion<T, Tint>(C);
  MyMatrix<T> CrossGram = C_T * GramMat * B1_T.transpose();
  MyMatrix<T> Gram1inv = Inverse(Gram1);
  MyMatrix<T> Coeff = CrossGram * Gram1inv;
  MyMatrix<T> B2 = C_T - Coeff * B1_T;
  MyMatrix<T> Gram2 = B2 * GramMat * B2.transpose();
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
  MyMatrix<T> Ortho = B2 * GramMat * B1_T.transpose();
  if (!IsZeroMatrix(Ortho)) {
    std::cerr << "IVF: the projection B2 is not A-orthogonal to L_1\n";
    throw TerminalException{1};
  }
#endif
  // The recursion is on a rational Gram matrix.  Since V_cv only depends on
  // argmin sets it is invariant under scaling of the Gram matrix, but the
  // rational entries are kept as they are, no rescaling is needed for
  // correctness.
  /*
    The recursion may return one vector per pair. That loses nothing: the
    closest vectors to -v are the negatives of those to v, so the two
    contribute the same set once the signs are normalized.
   */
  std::optional<MyMatrix<Tint>> opt_V2 = CharacteristicVectorSetCV<T, Tint>(
      Gram2, spanning, antipodal_half, budget, os);
  if (!opt_V2) {
    return {};
  }
  MyMatrix<Tint> const &V2 = *opt_V2;
  CVPSolver<T, Tint> solver(GramMat, os);
  CVPSolver<T, Tint> solver1(Gram1, os);
  // Projection of Z^n onto L_2 written in the B_2 coordinates, that is
  // x |-> x - (x A B_1^T) Gram1^{-1} B_1 followed by the B_2 coordinates.
  // In the basis (B_1, C) of Z^n the projection of B_1 is 0 and that of C
  // is B_2, so the B_2 coordinate of proj(x) is the C part of x.
  MyMatrix<Tint> FullBasis(n, n);
  for (int i = 0; i < r; i++) {
    FullBasis.row(i) = B1.row(i);
  }
  for (int i = 0; i < n - r; i++) {
    FullBasis.row(r + i) = C.row(i);
  }
  for (int i_row = 0; i_row < V2.rows(); i_row++) {
    // Two closest-vector computations per element of V2.
    if (budget.expired() || budget.too_big(acc.size())) {
      return {};
    }
    MyVector<Tint> y = GetMatrixRow(V2, i_row);
    MyVector<T> y_T = UniversalVectorConversion<T, Tint>(y);
    // v = y B_2 in R^n.
    MyVector<T> v = B2.transpose() * y_T;
    resultCVP<T, Tint> res = solver.nearest_vectors(v);
    acc.insert_rows(res.ListVect);
    // The closest vectors of the L_1-coset of v itself.  See the header
    // comment: CV(A, v) alone can miss that coset.
    MyVector<Tint> x0(n);
    for (int i = 0; i < r; i++) {
      x0(i) = 0;
    }
    for (int i = 0; i < n - r; i++) {
      x0(r + i) = y(i);
    }
    MyVector<Tint> xLift = FullBasis.transpose() * x0;
    // The coset is xLift + L_1, so we look for the point of L_1 closest to
    // v - xLift, in the B_1 coordinates.
    MyVector<T> xLift_T = UniversalVectorConversion<T, Tint>(xLift);
    MyVector<T> resid = v - xLift_T;
    // The B_1 coordinates of the L_1 component of resid.
    MyVector<T> residB1 = Gram1inv * (B1_T * (GramMat * resid));
    resultCVP<T, Tint> res1 = solver1.nearest_vectors(residB1);
    for (int i_r = 0; i_r < res1.ListVect.rows(); i_r++) {
      MyVector<Tint> z = GetMatrixRow(res1.ListVect, i_r);
      acc.insert(xLift + B1.transpose() * z);
    }
  }
  MyMatrix<Tint> RetMat = acc.get_matrix(n);
#ifdef TIMINGS_INVARIANT_VECTOR_FAMILY
  os << "|IVF: CharacteristicVectorSetCV|=" << time << "\n";
#endif
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: CV |V_cv|=" << RetMat.rows() << "\n";
#endif
  return RetMat;
}

/*
  The invariant vector family a canonicalization works from, together with
  whether it generates Z^n.

  Which family is best is not decidable in advance. The shell based
  ExtractInvariantVectorFamily is cheap when the successive minima are close
  together and catastrophic when they are not, a whole shell having to be
  taken at once: 16942 vectors on TestData/SlowCanonic/slow_canonic_2, 16816
  of them in the single shell of norm 14. V_cv does not care how far apart
  the minima are but pays 2^k closest-vector computations when the minimal
  vectors span a sublattice of index 2^k. Neither dominates: on
  slow_canonic_1 the shells give 2358 vectors against 8516 for the spanning
  V_cv, on slow_canonic_2 it is 16942 against 238.

  Building a family costs a fraction of a percent of what is done with it,
  the weight matrix being quadratic in its size and the canonical labelling
  worse than quadratic, so the cheapest way to avoid the bad case of either
  is to build both and keep the smaller. Both sizes are invariants of the
  isometry class, so the choice is itself canonical and the canonical forms
  obtained from it remain comparable.
 */
template <typename Tint> struct CanonicVectorFamily {
  /*
    One vector per antipodal pair. The families are closed under v -> -v, and
    everything downstream is quadratic in the number of vectors, so they are
    built halved rather than built whole and halved afterwards. Consumers
    that need both signs, the ones computing the automorphism group of the
    configuration as a permutation group, call get_full.
   */
  MyMatrix<Tint> SHVhalf;
  // Measured by FamilySpansLattice, never inferred from how the family was
  // built: a family built without asking to span usually spans anyway, and
  // every family that does can take the cheap canonicalization. The half
  // family spans what the whole one does, -v being in the span of v.
  bool spans_lattice;
  /*
    The whole family, the chosen representatives first and their negatives
    after, in the same order. That is the indexing the generators lifted by
    AbsTrick_LiftGenerators use, [0, nbPair) for +v and [nbPair, 2 nbPair)
    for -v, so the two agree without a translation table.
   */
  MyMatrix<Tint> get_full() const {
    int nbPair = SHVhalf.rows();
    int n = SHVhalf.cols();
    MyMatrix<Tint> SHV(2 * nbPair, n);
    for (int iPair = 0; iPair < nbPair; iPair++) {
      for (int i = 0; i < n; i++) {
        SHV(iPair, i) = SHVhalf(iPair, i);
        SHV(nbPair + iPair, i) = -SHVhalf(iPair, i);
      }
    }
    return SHV;
  }
  int n_pair() const { return SHVhalf.rows(); }
};

/*
  Whether the rows of SHV generate Z^n rather than a finite index subgroup.
  A Z-basis of the span followed by its determinant, which costs nothing
  next to what the answer saves.
 */
template <typename Tint>
bool FamilySpansLattice(MyMatrix<Tint> const &SHV) {
  int n = SHV.cols();
  if (RankMat(SHV) != n) {
    return false;
  }
  MyMatrix<Tint> Basis = GetZbasis(SHV);
  Tint det = DeterminantMat(Basis);
  return T_abs(det) == Tint(1);
}

template <typename Tint>
CanonicVectorFamily<Tint> get_canonic_vector_family(MyMatrix<Tint> &&SHVhalf) {
  bool spans = FamilySpansLattice<Tint>(SHVhalf);
  return {std::move(SHVhalf), spans};
}

/*
  The smaller of the two families, both stopped at full rank, found by
  running them against each other rather than by building either to the end.

  Neither may be built blindly. The shells grow in jumps, a shell having to
  be taken whole, and reach 76448 pairs on a genus of determinant 351 and
  rank 14 where V_cv gives 1576. But V_cv is not the cheap one either: on the
  lattices of that same genus it takes 3.0 s against 1.0 s for the shells,
  and it pays 2^k closest-vector computations when the minimal vectors span a
  sublattice of index 2^k. Which is cheaper is a property of the lattice and
  is not known in advance.

  So the two are interleaved, as test_finiteness_group of
  src_latt/FiniteMatrixGroupTest.h interleaves its two methods: one shell is
  added, the time it took is handed to V_cv as a budget, and so on until one
  of them finishes. Then the other is continued with the SIZE of the winner
  as a bound, and abandoned if it passes it.

  The timing decides only which construction to advance next. The answer is
  always "the smaller family, V_cv on a tie", which depends on the lattice
  alone. That distinction matters: the canonical form depends on which family
  is used, so a choice made on timing would give two isometric lattices
  different canonical forms from one run to the next, and the genus
  enumeration would count them as distinct classes and never reach its mass.
 */
template <typename Tfield, typename Tint>
CanonicVectorFamily<Tint>
GetCanonicVectorFamily_kernel(MyMatrix<Tfield> const &GramMat,
                              std::ostream &os) {
  using T = Tfield;
#ifdef TIMINGS_INVARIANT_VECTOR_FAMILY
  MicrosecondTime time_tot;
#endif
  const bool antipodal_half = true;
  int n = GramMat.rows();
  auto f_final = [&]([[maybe_unused]] MyMatrix<Tint> const &M,
                     FundInvariantVectorFamily<Tint> const &fi) -> bool {
    return fi.rank == n;
  };
  auto f_invariant =
      [&](MyMatrix<Tint> const &M) -> FundInvariantVectorFamily<Tint> {
    return ComputeRankInvariant(M, os);
  };
  auto f_cv = [&](int64_t const &ns,
                  size_t const &sz) -> std::optional<MyMatrix<Tint>> {
    CharVectSet_Budget budget(ns, sz);
    return CharacteristicVectorSetCV<T, Tint>(GramMat, false, antipodal_half,
                                              budget, os);
  };
  auto f_ret = [&](MyMatrix<Tint> &&SHV,
                   [[maybe_unused]] const char *who) -> CanonicVectorFamily<Tint> {
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
    os << "IVF: keeping " << who << ", " << SHV.rows() << " pairs\n";
#endif
#ifdef TIMINGS_INVARIANT_VECTOR_FAMILY
    os << "|IVF: GetCanonicVectorFamily|=" << time_tot << "\n";
#endif
    return get_canonic_vector_family<Tint>(std::move(SHV));
  };
  ShellFamilyBuilder<T, Tint, decltype(f_final), decltype(f_invariant)> shells(
      GramMat, f_final, f_invariant, os);
  int64_t shell_ns = 0;
  while (true) {
    NanosecondTime time_shell;
    shells.one_shell(os);
    shell_ns += time_shell.const_eval_int64();
    if (shells.is_done()) {
      // The shells won the race. V_cv is kept only if it is smaller, and it
      // is bounded by their size so that a hopeless attempt is cut short.
      std::optional<MyMatrix<Tint>> opt = f_cv(0, shells.size());
      if (opt && static_cast<size_t>(opt->rows()) <= shells.size()) {
        return f_ret(std::move(*opt), "V_cv");
      }
      MyMatrix<Tint> SHV_shell = shells.get();
      return f_ret(std::move(SHV_shell), "the shells");
    }
    // V_cv is given the time the shells have spent so far. It is restarted
    // rather than resumed, which wastes at most the previous attempt.
    std::optional<MyMatrix<Tint>> opt = f_cv(shell_ns, 0);
    if (opt) {
      // V_cv won the race. The shells are kept only if they are smaller, and
      // they are abandoned as soon as they pass its size.
      size_t cv_size = opt->rows();
      while (!shells.is_done() && shells.size() <= cv_size) {
        shells.one_shell(os);
      }
      if (shells.is_done() && shells.size() < cv_size) {
        MyMatrix<Tint> SHV_shell = shells.get();
        return f_ret(std::move(SHV_shell), "the shells");
      }
      return f_ret(std::move(*opt), "V_cv");
    }
  }
}

/*
  The same for a caller whose Gram matrix is integral.

  The construction has to be done over a field whatever the caller holds:
  V_cv projects onto the orthogonal of the span of the minimal vectors, and
  that projection inverts the Gram matrix of that span. Over a ring the two
  inversions in CharacteristicVectorSetCV and
  CharacteristicVectorSetWellRoundedCV would divide with truncation and
  return a family that is not the one asked for, silently, and only on the
  lattices whose minimal vectors do not have full rank.

  The family itself is integral, so nothing rational reaches the caller, and
  the canonicalization built on it can be done over the ring.
 */
template <typename T, typename Tint>
CanonicVectorFamily<Tint> GetCanonicVectorFamily(MyMatrix<T> const &GramMat,
                                                 std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  if constexpr (std::is_same_v<T, Tfield>) {
    return GetCanonicVectorFamily_kernel<T, Tint>(GramMat, os);
  } else {
    MyMatrix<Tfield> GramMat_F = UniversalMatrixConversion<Tfield, T>(GramMat);
    return GetCanonicVectorFamily_kernel<Tfield, Tint>(GramMat_F, os);
  }
}

// clang-format off
#endif  // SRC_LATT_INVARIANTVECTORFAMILY_H_
// clang-format on

// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_PERFECTFORM_H_
#define SRC_PERFECT_PERFECTFORM_H_

// clang-format off>
#include "MatrixGroup.h"
#include "PolytopeEquiStab.h"
#include "EquiStabMemoization.h"
#include "shortest_flipping.h"
#include "Positivity.h"
#include "Tspace_General.h"
#include "fractions.h"
#include "POLY_RecursiveDualDesc.h"
#include "MatrixGroupAverage.h"
#include <map>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_PERFECT_FORM
#define DEBUG_INITIAL_PERFECT
#define DEBUG_PERFECT_REPR
#define DEBUG_BOUNDED_FACE
#endif

#ifdef DISABLE_DEBUG_PERFECT_FORM
#undef DEBUG_PERFECT_FORM
#endif

#ifdef DISABLE_DEBUG_INITIAL_PERFECT
#undef DEBUG_INITIAL_PERFECT
#endif

#ifdef DISABLE_DEBUG_PERFECT_REPR
#undef DEBUG_PERFECT_REPR
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_INITIAL_PERFECT
#define SANITY_CHECK_PERFECT_REPR
#define SANITY_CHECK_PERFECT_FORM
#define SANITY_CHECK_BOUNDED_FACE
#define SANITY_CHECK_POSITIVE_SEMIDEFINITE_SELF_DUAL
#endif

#ifdef TIMINGS
#define TIMINGS_PERFECT_FORM
#endif

template <typename T, typename Tgroup> struct RyshkovGRP {
  MyMatrix<T> PerfDomEXT;
  Tgroup GRPsub;
  vectface ListIncd;
  std::vector<int> ListPos;
  std::vector<typename Tgroup::Telt::Tidx> ListBlockFirst;
  typename Tgroup::Telt map_elt(typename Tgroup::Telt const& x) const {
    using Telt = typename Tgroup::Telt;
    using Tidx = typename Telt::Tidx;
    size_t nbBlock = ListBlockFirst.size();
    std::vector<Tidx> eList(nbBlock);
    for (size_t iBlock = 0; iBlock < nbBlock; iBlock++) {
      Tidx iSHV = ListBlockFirst[iBlock];
      Tidx jSHV = x.at(iSHV);
      Tidx jBlock = ListPos[jSHV];
      eList[iBlock] = jBlock;
    }
    Telt elt(eList);
    return elt;
  }
};

template <typename T, typename Tgroup>
RyshkovGRP<T, Tgroup>
GetNakedPerfectCone_GRP(LinSpaceMatrix<T> const &LinSpa,
                        MyMatrix<T> const &SHV_T,
                        Tgroup const &GRP, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int nbSHV = SHV_T.rows();
  std::vector<int> ListPos(nbSHV);
  int nbMat = LinSpa.ListMat.size();
  int n = LinSpa.n;
  MyMatrix<T> RyshkovLoc(nbSHV, nbMat);
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<T> eVect = GetMatrixRow(SHV_T, iSHV);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      T eSum = EvaluationQuadForm<T, T>(LinSpa.ListMat[iMat], eVect);
      RyshkovLoc(iSHV, iMat) = eSum;
    }
  }
  //  os << "RyshkovLoc=\n";
  //  WriteMatrix(os, RyshkovLoc);
  std::map<MyVector<T>, std::vector<int>> map;
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<T> V = GetMatrixRow(RyshkovLoc, iSHV);
    map[V].push_back(iSHV);
  }
  std::vector<std::vector<int>> ListBlock;
  int i_block = 0;
  for (auto &kv : map) {
    ListBlock.push_back(kv.second);
    for (auto &iSHV : kv.second) {
      ListPos[iSHV] = i_block;
    }
    i_block += 1;
  }
  int nbBlock = i_block;
#ifdef DEBUG_PERFECT_FORM
  os << "m=" << n << " nbBlock=" << nbBlock << " nbSHV=" << nbSHV
     << " nbMat=" << nbMat << "\n";
#endif
  MyMatrix<T> PerfDomEXT(nbBlock, nbMat);
  MyMatrix<T> SHVred(nbBlock, n);
  for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
    int iSHV = ListBlock[iBlock][0];
    MyVector<T> eVect = GetMatrixRow(SHV_T, iSHV);
    AssignMatrixRow(SHVred, iBlock, eVect);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      PerfDomEXT(iBlock, iMat) = RyshkovLoc(iSHV, iMat);
    }
  }
  std::vector<Telt> l_gens;
  std::vector<Tidx> ListBlockFirst;
  for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
    Tidx iSHV = ListBlock[iBlock][0];
    ListBlockFirst.push_back(iSHV);
  }
  for (auto &eGen : GRP.GeneratorsOfGroup()) {
    std::vector<Tidx> eList(nbBlock);
    for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
      Tidx iSHV = ListBlockFirst[iBlock];
      Tidx jSHV = eGen.at(iSHV);
      Tidx jBlock = ListPos[jSHV];
      eList[iBlock] = jBlock;
    }
    Telt eElt(eList);
    l_gens.emplace_back(std::move(eElt));
  }
  Tgroup GRPsub(l_gens, nbBlock);
  vectface ListIncd =
    DualDescriptionStandard<T, Tgroup>(PerfDomEXT, GRPsub, os);
  return {std::move(PerfDomEXT), std::move(GRPsub), std::move(ListIncd), std::move(ListPos), std::move(ListBlockFirst)};
}

template <typename T, typename Tgroup>
Face get_big_incd(RyshkovGRP<T,Tgroup> const& ryshk, Face const& f_sma) {
  int n_shv = ryshk.ListPos.size();
  Face f_big(n_shv);
  for (int i=0; i<n_shv; i++) {
    int i_block = ryshk.ListPos[i];
    if (f_sma[i_block] == 1) {
      f_big[i] = 1;
    }
  }
  return f_big;
}

template <typename T, typename Tint>
std::pair<MyMatrix<T>, Tshortest<T, Tint>>
Flipping_Perfect(MyMatrix<T> const &eMatIn, MyMatrix<T> const &eMatDir,
                 std::ostream &os) {
  auto f_admissible = [&os](MyMatrix<T> const &eMat) -> bool {
    return IsPositiveDefinite<T>(eMat, os);
  };
  auto f_shortest = [&os](MyMatrix<T> const &eMat) -> Tshortest<T, Tint> {
    return T_ShortestVectorHalf<T, Tint>(eMat, os);
  };
  return Kernel_Flipping_Perfect<T, Tint, decltype(f_admissible),
                                 decltype(f_shortest)>(f_admissible, f_shortest,
                                                       eMatIn, eMatDir, os);
}

template <typename T, typename Tint>
MyMatrix<T> get_scal_mat(std::vector<MyMatrix<T>> const &ListMat,
                         Tshortest<T, Tint> const &rec_shv) {
  int nbMat = ListMat.size();
  int nbShort = rec_shv.SHV.rows();
  MyMatrix<T> ScalMat(nbShort, nbMat);
  for (int iShort = 0; iShort < nbShort; iShort++) {
    MyVector<Tint> eVectShort = rec_shv.SHV.row(iShort);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      T eNorm = EvaluationQuadForm<T, Tint>(ListMat[iMat], eVectShort);
      ScalMat(iShort, iMat) = eNorm;
    }
  }
  return ScalMat;
}

template <typename T, typename Tint>
bool is_perfect_in_space(LinSpaceMatrix<T> const &LinSpa,
                         Tshortest<T, Tint> const &rec_shv) {
  int nbMat = LinSpa.ListMat.size();
  MyMatrix<T> ScalMat = get_scal_mat<T, Tint>(LinSpa.ListMat, rec_shv);
  return RankMat(ScalMat) == nbMat;
}

/*
  Look for a positive definite matrix in the T-space such that
  lambda = A[v] for v in SHV
  ---
  That function does not address the issue for which it was created
  (finding bounded
  but it has independent interest.
 */
template<typename T, typename Tint>
bool find_positive_definite_shv_equal(LinSpaceMatrix<T> const &LinSpa, MyMatrix<Tint> const& SHV, std::ostream& os) {
  int n = LinSpa.n;
  std::vector<MyVector<Tint>> ListVect = VectorFamilyFromMatrix(SHV);
  std::vector<MyVector<Tint>> TestV_init = get_initial_vector_test_v<Tint>(n, ListVect, os);
  std::vector<MyMatrix<T>> BasisSpace = get_spec_basis_shv_equal(LinSpa.ListMat, ListVect);
  std::optional<MyMatrix<T>> opt = GetOnePositiveDefiniteMatrix_ListV<T,Tint>(BasisSpace, TestV_init, os);
  if (opt) {
    return false;
  } else {
    return true;
  }
}

/*
  Look for a positive semi-definite matrix in the T-space such that
  0 = A[v] for v in SHV
  ---
  That function does not address the issue for which it was created
  but it has independent interest.

 */
template<typename T, typename Tint>
bool find_positive_semidefinite_shv_zero(LinSpaceMatrix<T> const &LinSpa, MyMatrix<Tint> const& SHV, std::ostream& os) {
#ifdef DEBUG_BOUNDED_FACE
  os << "PERF: is_bounded_face_iterative_bis, step 1\n";
#endif
  int n = LinSpa.n;
  std::vector<MyVector<Tint>> ListVect = VectorFamilyFromMatrix(SHV);
  std::vector<MyVector<Tint>> TestV_init = get_initial_vector_test_v<Tint>(n, ListVect, os);
  std::vector<MyMatrix<T>> BasisSpace = get_spec_basis_shv_zero<T,Tint>(LinSpa.ListMat, ListVect);
#ifdef DEBUG_BOUNDED_FACE
  os << "PERF: is_bounded_face_iterative_bis, step 5\n";
#endif
  std::optional<MyMatrix<T>> opt = GetOnePositiveSemiDefiniteMatrix_ListV<T,Tint>(BasisSpace, TestV_init, os);
  if (opt) {
    return false;
  } else {
    return true;
  }
}

/*
  Look for a positive semi-definite matrix in the T-space such that
  Min(A) = SHV
  ---
  The difficulty is that the set SHV is not necessarily full dimensional.
  * This makes computing the symmetry group extremely difficult.
  * One possible improvement could be to find a quadratic form which
    is positive semidefinite and whose kernel is equal to the SHV.
 */
template<typename T, typename Tint>
std::optional<MyMatrix<T>> is_bounded_face_iterative(LinSpaceMatrix<T> const &LinSpa, MyMatrix<Tint> const& SHV, std::ostream& os) {

  std::vector<MyVector<Tint>> ListVectWork = VectorFamilyFromMatrix(SHV);
  std::vector<MyMatrix<T>> const& ListMat = LinSpa.ListMat;
  bool NoExtension = true;
  ReplyRealizability<T, Tint> RecTest = SHORT_TestRealizabilityShortestFamily_Raw<T, Tint>(ListVectWork, ListMat, NoExtension, os);
  if (!RecTest.reply) {
    return {};
  }
  return RecTest.eMat;
}

/*
  We consider here the question of well roundedness of
  the faces of the cellular complex.

  Known properties of those vector configurations:
  * Any cell is the set of minimum vectors of a positive
    definite quadratic form. This can be proved by taking
    the average of the quadratic forms corresponding to
    the perfect domains containing it.

  * The key property in question has to be descending.
    That is if X does not satisfy P, then faces of X do
    not satisfy it as well.
    + The infinite stabilizer. If a face F has an
      infinite stabilizer and a subface G\subset F has
      a finite stabilizer. That seems impossible.
    + The condition that there exist a positive semidefinite
      form in the T-space which is zero on F. It is of
      course true for the subface trivially. [It does
      not appear sensible to consider the form being
      zero on the whole linear span of the subset]
      . That property appears very hard to prove on its
        own.
      . However, if the cell in question might be maximal
        for the property, then for each cell containing G
        and which does not have this property, we can provide
        some vectors which should be strictly positive. This
        might help the linear program if this is something we
        want to pursue
    + Whether the set of shortest vectors are of full rank.
      That property can be tested relatively easily. It is
      not an adequate notion, but it exists.

  Specific notions:
  * The self-duality of some cones. Some cones are self-dual
    and some are not.
    + The theory seems to align very well:
      . The real places, we embed into S^n, the complex in
        H^n.
      . The expression of the quadratic form as a sum of
        individual quadratic forms appear sensible and not
        break the self-duality.
      . For the case of groups
    + The imaginary quadratic case might be the typical
      example of this kind. It does appear to work.
    + We have to understand how this works out because the
      use of the inverse matrix A0 is specific.
    + Meanwhile the check that the corresponding matrix
      is positive semidefinite is a good test of
      correctness of the algorithm.
    + Then the face is bounded if and only if the
      corresponding form is positive definite.

  * The span of the vectors is full dimensional.
    This property on page 462 of
    A. Ash, "Small dimensional classifying spaces for
    arithmetic subgroups of general linear groups"
    is M(L) S = S^n.
    Points:
    + That property can be interpreted in our setting.
    + We provide a number of generators of the ring.
      If the rank of the span of the shortest vectors
      under it is full dimensional, then we are in
      the bounded case.

  So, we have 3 properties that should be equivalent
  for the bounded face:
  * finite stabilizer.
  * self-dual property.
  * spanning property.
  If we pass those tests, then we should be reasonably
  confident about our understanding of the situation.

  ----

  Design:
  * Self-dual property and spanning property are easy
    to resolve.
  * But the finite stabilizer requires computing the
    stabilizer. That is expensive.
  * So, the self-dual and spanning properties are
    determined if they are available. If they cannot
    be computed, then marked as None.
  * If full dimensional directly, then mark everything
    as {true}
  * The bounded finite stabilizer requires more work.
    But may be needed if we have neither spanning nor
    self-duality.
 */
struct PerfectBoundednessProperty {
  std::optional<bool> bounded_self_dual;
  std::optional<bool> bounded_spanning;
  std::optional<bool> bounded_finite_stabilizer;
};

void check_pbp(PerfectBoundednessProperty const& pbp) {
  std::unordered_map<std::string, bool> map_result;
  std::set<bool> results;
  if (pbp.bounded_self_dual) {
    bool val = *pbp.bounded_self_dual;
    map_result["self_dual"] = val;
    results.insert(val);
  }
  if (pbp.bounded_spanning) {
    bool val = *pbp.bounded_spanning;
    map_result["spanning"] = val;
    results.insert(val);
  }
  if (pbp.bounded_finite_stabilizer) {
    bool val = *pbp.bounded_finite_stabilizer;
    map_result["finite_stabilizer"] = val;
    results.insert(val);
  }
  if (results.size() != 1) {
    std::cerr << "Inconsistent results\n";
    for (auto & kv: map_result) {
      std::cerr << "key=" << kv.first << " result=" << kv.second << "\n";
    }
    throw TerminalException{1};
  }
}

bool get_result(PerfectBoundednessProperty const& pbp) {
  if (pbp.bounded_self_dual) {
    return *pbp.bounded_self_dual;
  }
  if (pbp.bounded_spanning) {
    return *pbp.bounded_spanning;
  }
  if (pbp.bounded_finite_stabilizer) {
    return *pbp.bounded_finite_stabilizer;
  }
  std::cerr << "One of the methods should have worked\n";
  throw TerminalException{1};
}

std::string to_string(PerfectBoundednessProperty const& pbp) {
  std::string str_out;
  //
  auto str_opt_bool=[&](std::optional<bool> const& opt) -> std::string {
    if (opt) {
      bool val = *opt;
      return "Some(" + GAP_logical(val) + ")";
    } else {
      return "None";
    }
  };
  //
  str_out += "bounded_self_dual=";
  str_out += str_opt_bool(pbp.bounded_self_dual);
  str_out += ", bounded_spanning=";
  str_out += str_opt_bool(pbp.bounded_spanning);
  str_out += ", bounded_finite_stabilizer=";
  str_out += str_opt_bool(pbp.bounded_finite_stabilizer);
  //
  return str_out;
}


template<typename T, typename Tint>
bool is_bounded_self_dual(PairwiseScalarInfo<T> const& psi, std::vector<MyMatrix<T>> const& ListMat, MyMatrix<Tint> const& SHV, [[maybe_unused]] std::ostream &os) {
  int n_vect = SHV.rows();
  int n = SHV.cols();
  int n_mat = ListMat.size();
  MyVector<T> SumRay = ZeroVector<T>(n_mat);
  for (int i_vect=0; i_vect<n_vect; i_vect++) {
    MyVector<Tint> V = GetMatrixRow(SHV, i_vect);
    for (int i_mat=0; i_mat<n_mat; i_mat++) {
      T val = EvaluationQuadForm<T, Tint>(ListMat[i_mat], V);
      SumRay(i_mat) += val;
    }
  }
  MyVector<T> ExprBasis = psi.PairwiseScalarInv * SumRay;
  MyMatrix<T> SumMat = ZeroMatrix<T>(n,n);
  for (int i_mat=0; i_mat<n_mat; i_mat++) {
    SumMat += ExprBasis(i_mat) * ListMat[i_mat];
  }
#ifdef SANITY_CHECK_POSITIVE_SEMIDEFINITE_SELF_DUAL
  if (!IsPositiveSemiDefinite(SumMat, os)) {
    std::cerr << "The matrix should be positive semi-definite\n";
    std::cerr << "So, maybe the cone is not self-dual\n";
    throw TerminalException{1};
  }
#endif
  int rnk = RankMat(SumMat);
  if (rnk == n) {
    return true;
  } else {
    return false;
  }
}

template<typename T, typename Tint>
bool is_bounded_spanning_elements(std::vector<MyMatrix<T>> const& l_spanning_elements, MyMatrix<Tint> const& SHV, [[maybe_unused]] std::ostream &os) {
  int n = SHV.cols();
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T,Tint>(SHV);
  MyMatrix<T> SHV_spann = DirectSpannEquivariantSpace(SHV_T, l_spanning_elements);
  int rnk = RankMat(SHV_spann);
  if (rnk == n) {
    return true;
  } else {
    return false;
  }
}


template<typename T, typename Tint>
PerfectBoundednessProperty initial_bounded_property(LinSpaceMatrix<T> const &LinSpa, MyMatrix<Tint> const& SHV, std::ostream& os) {
  int n = LinSpa.n;
  int rnk = RankMat(SHV);
  std::optional<bool> bounded_self_dual;
  std::optional<bool> bounded_spanning;
  std::optional<bool> bounded_finite_stabilizer;
  PerfectBoundednessProperty pbp{bounded_self_dual, bounded_spanning, bounded_finite_stabilizer};

  if (rnk == n) {
    // For the case of classic Voronoi, that would be a necessary and
    // sufficient condition.
    // That case should be very frequent and so this accelerate.
    pbp.bounded_self_dual = {true};
    pbp.bounded_spanning = {true};
    pbp.bounded_finite_stabilizer = {true};
    return pbp;
  }
  if (LinSpa.is_self_dual) {
    bool val = is_bounded_self_dual(LinSpa.pairwise_scalar_info, LinSpa.ListMat, SHV, os);
    pbp.bounded_self_dual = val;
  }
  if (LinSpa.l_spanning_elements.size() > 0) {
    bool val = is_bounded_spanning_elements(LinSpa.l_spanning_elements, SHV, os);
    pbp.bounded_spanning = val;
  }
  return pbp;
}


/*
  Finds a perfect form by starting from the super matrix of the space
  and then doing some moves if that is not perfect.
  ---
  The direction of the movements are obtained by computing the
  shortest vectors, getting the kernel.
  The vector in the kernel gives direction of change.
  However, if that direction is positive definite then we will not obtain
  a new form. This scenario is explained in "Enumerating Perfect Forms"
  arXiv:0901.1587
  ---
  During the running of this test, the minimal norm of the matrix remains the
  same, but the number of vectors attaining it is increasing.
 */
template <typename T, typename Tint>
std::pair<MyMatrix<T>, Tshortest<T, Tint>>
GetOnePerfectForm_Kernel(std::vector<MyMatrix<T>> const& ListMat, MyMatrix<T> const& InitialMat, std::ostream &os) {
  int nbMat = ListMat.size();
  MyMatrix<T> ThePerfMat = InitialMat;
  Tshortest<T, Tint> rec_shv = T_ShortestVectorHalf<T, Tint>(ThePerfMat, os);
#ifdef SANITY_CHECK_INITIAL_PERFECT
  T TheMin = rec_shv.min;
#endif
#ifdef DEBUG_INITIAL_PERFECT
  int iter = 0;
#endif
  while (true) {
    MyMatrix<T> ScalMat = get_scal_mat<T, Tint>(ListMat, rec_shv);
    SelectionRowCol<T> eSelect = TMat_SelectRowCol(ScalMat);
    int TheRank = eSelect.TheRank;
#ifdef DEBUG_INITIAL_PERFECT
    os << "PERF: GetOnePerfectForm, iter=" << iter << " min=" << rec_shv.min
       << " |SHV|=" << rec_shv.SHV.rows() << "\n";
#endif
    if (TheRank == nbMat) {
#ifdef DEBUG_INITIAL_PERFECT
      os << "PERF: GetOnePerfectForm, returning at iter=" << iter << "\n";
#endif
      return {std::move(ThePerfMat), std::move(rec_shv)};
    }
    MyVector<T> V = eSelect.NSP.row(0);
    auto iife_get_dir = [&]() -> MyMatrix<T> {
      MyMatrix<T> M = GetMatrixFromBasis(ListMat, V);
      if (IsPositiveDefinite<T>(M, os)) {
        // For a positive definite matrix, we need to take the opposite
        // because we need a direction outside of the cone.
        return -M;
      }
      return M;
    };
    MyMatrix<T> DirMat = iife_get_dir();
#ifdef DEBUG_INITIAL_PERFECT
    os << "PERF: GetOnePerfectForm, iter=" << iter << " DirMat=\n";
    WriteMatrix(os, DirMat);
#endif
    auto pair = Flipping_Perfect<T, Tint>(ThePerfMat, DirMat, os);
    ThePerfMat = pair.first;
    rec_shv = pair.second;
#ifdef SANITY_CHECK_INITIAL_PERFECT
    if (rec_shv.min != TheMin) {
      std::cerr << "PERF: The rec_shv minimum should remain invariant\n";
      throw TerminalException{1};
    }
#endif
#ifdef DEBUG_INITIAL_PERFECT
    iter += 1;
#endif
  }
}

template <typename T, typename Tint>
std::pair<MyMatrix<T>, Tshortest<T, Tint>>
GetOnePerfectForm(LinSpaceMatrix<T> const &LinSpa, std::ostream &os) {
  return GetOnePerfectForm_Kernel<T,Tint>(LinSpa.ListMat, LinSpa.SuperMat, os);
}



template <typename T, typename Tint> struct SimplePerfect {
  MyMatrix<T> Gram;
  Tshortest<T, Tint> rec_shv;
};

template <typename T, typename Tint>
std::istream &operator>>(std::istream &is, SimplePerfect<T, Tint> &obj) {
  MyMatrix<T> eG = ReadMatrix<T>(is);
  MyMatrix<Tint> SHV = ReadMatrix<Tint>(is);
  MyVector<Tint> V = GetMatrixRow(SHV, 0);
  T min = EvaluationQuadForm<T, Tint>(eG, V);
  Tshortest<T, Tint> rec_shv{min, SHV};
  obj = {eG, rec_shv};
  return is;
}

template <typename T, typename Tint>
std::ostream &operator<<(std::ostream &os, SimplePerfect<T, Tint> const &obj) {
  WriteMatrix(os, obj.Gram);
  WriteMatrix(os, obj.rec_shv.SHV);
  return os;
}

template <typename T, typename Tint>
MyMatrix<T> conversion_and_duplication(MyMatrix<Tint> const &SHV) {
  int dim = SHV.cols();
  int nbSHV = SHV.rows();
  MyMatrix<T> SHV_T(2 * nbSHV, dim);
  for (int i_row = 0; i_row < nbSHV; i_row++) {
    for (int i = 0; i < dim; i++) {
      T val = UniversalScalarConversion<T, Tint>(SHV(i_row, i));
      SHV_T(2 * i_row, i) = val;
      SHV_T(2 * i_row + 1, i) = -val;
    }
  }
  return SHV_T;
}

template <typename T, typename Tint>
MyMatrix<T> get_fulldim_zbasis_shv_t(MyMatrix<T> const &eMat, Tshortest<T, Tint> const &rec_shv,
                      std::ostream &os) {
  MyMatrix<T> SHVorig_T = UniversalMatrixConversion<T, Tint>(rec_shv.SHV);
  if (IsFullDimZbasis(rec_shv.SHV, os)) {
    return conversion_and_duplication<T, Tint>(rec_shv.SHV);
  }
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(eMat, os);
  return UniversalMatrixConversion<T, Tint>(SHV);
}

template <typename T, typename Tint>
MyMatrix<Tint> get_fulldim_shv_tint(MyMatrix<T> const &eMat, Tshortest<T, Tint> const &rec_shv,
                              std::ostream &os) {
  MyMatrix<T> SHVorig_T = UniversalMatrixConversion<T, Tint>(rec_shv.SHV);
  if (IsFullDim(rec_shv.SHV, os)) {
    return conversion_and_duplication<Tint, Tint>(rec_shv.SHV);
  }
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyFullRank<T, Tint>(eMat, os);
  return SHV;
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>> SimplePerfect_TestEquivalence(
    LinSpaceMatrix<T> const &LinSpa, MyMatrix<T> const &eMat1,
    MyMatrix<T> const &eMat2, Tshortest<T, Tint> const &rec_shv1,
    Tshortest<T, Tint> const &rec_shv2, std::ostream &os) {
  MyMatrix<T> SHV1_T = get_fulldim_zbasis_shv_t(eMat1, rec_shv1, os);
  MyMatrix<T> SHV2_T = get_fulldim_zbasis_shv_t(eMat2, rec_shv2, os);
#ifdef SANITY_CHECK_PERFECT_REPR
  if (has_duplication(SHV1_T)) {
    std::cerr << "PERF: SHV1 has duplication\n";
    throw TerminalException{1};
  }
  if (has_duplication(SHV2_T)) {
    std::cerr << "PERF: SHV2 has duplication\n";
    throw TerminalException{1};
  }
  if (!is_antipodal(SHV1_T)) {
    std::cerr << "PERF: SHV1 is not antipodal\n";
    throw TerminalException{1};
  }
  if (!is_antipodal(SHV2_T)) {
    std::cerr << "PERF: SHV2 is not antipodal\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_PERFECT_REPR
  os << "PERF: TestEquivalence, before LINSPA_TestEquivalenceGramMatrix_SHV\n";
  os << "PERF: TestEquivalence, eMat1=\n";
  WriteMatrix(os, RemoveFractionMatrix(eMat1));
  os << "PERF: TestEquivalence, eMat2=\n";
  WriteMatrix(os, RemoveFractionMatrix(eMat2));
#endif
  std::optional<MyMatrix<T>> opt =
      LINSPA_TestEquivalenceGramMatrix_SHV<T, Tgroup>(LinSpa, eMat1, eMat2,
                                                      SHV1_T, SHV2_T, os);
#ifdef DEBUG_PERFECT_REPR
  os << "PERF: TestEquivalence, det1=" << DeterminantMat(eMat1)
     << " det2=" << DeterminantMat(eMat2)
     << " opt.has_value()=" << opt.has_value() << "\n";
#endif
  if (!opt) {
    return {};
  }
  MyMatrix<T> const &M_T = *opt;
  MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(M_T);
  return M;
}

template <typename T, typename Tint>
size_t
SimplePerfect_Invariant(size_t const &seed, LinSpaceMatrix<T> const &LinSpa,
                        MyMatrix<T> const &eMat,
                        Tshortest<T, Tint> const &rec_shv, std::ostream &os) {
  MyMatrix<T> SHV_T = get_fulldim_zbasis_shv_t(eMat, rec_shv, os);
  return LINSPA_Invariant_SHV<T>(seed, LinSpa, eMat, SHV_T, os);
}

template <typename T, typename Tint, typename Tgroup>
std::pair<Tgroup, std::vector<MyMatrix<Tint>>>
SimplePerfect_Stabilizer(LinSpaceMatrix<T> const &LinSpa,
                         MyMatrix<T> const &eMat,
                         Tshortest<T, Tint> const &rec_shv, std::ostream &os) {
  //
  // Functionality for checking quality of equivalences
  //
  MyMatrix<T> SHVorig_T = conversion_and_duplication<T, Tint>(rec_shv.SHV);
  MyMatrix<T> SHV_T = get_fulldim_zbasis_shv_t(eMat, rec_shv, os);
  Result_ComputeStabilizer_SHV<T, Tgroup> result =
      LINSPA_ComputeStabilizer_SHV<T, Tgroup>(LinSpa, eMat, SHV_T, os);
#ifdef DEBUG_PERFECT_FORM
  os << "PERFECT: SimplePerfect_Stabilizer, we have result\n";
  os << "PERFECT: SimplePerfect_Stabilizer |SHVorig_T|=" << SHVorig_T.rows()
     << " |SHV_T|=" << SHV_T.rows() << "\n";
#endif
  if (TestEqualityMatrix(SHVorig_T, SHV_T)) {
    // This is the most likely scenario: The original
    // SHVorig is adequate
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(A), Case SHVorig_T == SHV_T\n";
#endif
    std::vector<MyMatrix<T>> l_matr_t =
        result.get_list_matrix(SHV_T, eMat, LinSpa, os);
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(A), We have l_matr_t\n";
#endif
    std::vector<MyMatrix<Tint>> l_matr =
        UniversalStdVectorMatrixConversion<Tint, T>(l_matr_t);
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(A), We have l_matr\n";
#endif
    Tgroup GRP = result.get_perm_group(SHV_T, os);
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(A), We have GRP\n";
#endif
    return {std::move(GRP), std::move(l_matr)};
  } else {
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(B), Case SHVorig_T != SHV_T\n";
#endif
    std::vector<MyMatrix<T>> l_matr_t =
        result.get_list_matrix(SHV_T, eMat, LinSpa, os);
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(B), We have l_matr_t\n";
#endif
    std::vector<MyMatrix<Tint>> l_matr =
        UniversalStdVectorMatrixConversion<Tint, T>(l_matr_t);
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(B), We have l_matr\n";
#endif
    Tgroup GRP =
        get_perm_group_from_list_matrices<T, Tgroup>(l_matr_t, SHVorig_T, os);
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(B), We have GRP\n";
    os << "PERFECT: SimplePerfect_Stabilizer(B), |SHVorig_T|=" << SHVorig_T.rows() << " / " << SHVorig_T.cols() << "\n";
    WriteMatrix(os, SHVorig_T);
    os << "PERFECT: SimplePerfect_Stabilizer(B),     |SHV_T|=" << SHV_T.rows() << " / " << SHV_T.cols() << "\n";
    WriteMatrix(os, SHV_T);
#endif
    return {std::move(GRP), std::move(l_matr)};
  }
}

// clang-format off
#endif  // SRC_PERFECT_PERFECTFORM_H_
// clang-format on

// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_PERFECT_LOW_DIM_H_
#define SRC_PERFECT_PERFECT_LOW_DIM_H_

// clang-format off
#include "perfect_complex.h"
#include "POLY_Kskeletton.h"
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_PERFECT_LOW_DIM
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_PERFECT_LOW_DIM
#endif

/*
  The cells of small rank of the perfect form complex, without the
  descent from the top.

  Every cell of the complex is a face of a perfect domain. So the cells of
  rank at most r_max are obtained by:
  + Enumerating, in every perfect domain, the orbits of faces of rank at
    most r_max under the stabilizer of the domain. This is done from the
    bottom, starting from the extreme rays, by the k-skeleton enumeration
    of POLY_Kskeletton.h with the linear programming method, so that the
    facets of the perfect domains are not needed.
  + Reducing the well rounded faces of the requested ranks up to the
    action of the group.

  The rank of a face is the dimension of the cone, that is the rank of its
  rays q(v) = v v^T in the T-space. A face of rank r is at the level
  dimEXT - r of the complex computed by perfect_complex.h.

  The rays of a perfect domain: the minimal vectors come in pairs +-v
  (EXT has the rows v, -v, v', -v', ...). Several pairs can give the same
  ray q(v): in a Hermitian space the multiples u v by the units u of the
  ring all do (the units i of Z[i] and the sixth roots of unity of the
  Eisenstein ring). The rays are therefore the distinct rows of the
  evaluations of the basis of the T-space on the vectors, and a face is
  made of all the vectors of its rays.
 */

template <typename T, typename Tint, typename Tgroup>
struct PerfectConeLowRankFaces {
  // One minimal vector per pair +-v; row i is row 2i of the EXT of the cone.
  MyMatrix<Tint> SHV;
  // The distinct rays q(v) in the coordinates of the T-space.
  MyMatrix<T> FAC;
  // ray_vectors[i_ray]: the rows of SHV giving the ray i_ray.
  std::vector<std::vector<int>> ray_vectors;
  // The action of the stabilizer of the cone on the rays.
  Tgroup GRP;
  // levels[r-1]: one representative per orbit of faces of rank r, as sets
  // of rays.
  std::vector<vectface> levels;
  // The vectors (rows of SHV) of a face given by its rays.
  MyMatrix<Tint> get_face_vectors(Face const &f) const {
    std::vector<int> l_row;
    for (size_t i_ray = f.find_first(); i_ray != Face::npos;
         i_ray = f.find_next(i_ray)) {
      for (int i_row : ray_vectors[i_ray]) {
        l_row.push_back(i_row);
      }
    }
    return SelectRow(SHV, l_row);
  }
};

/*
  The action of GRP_ext (on the rows of EXT) induced on the rays.
  ray_of_vect[i] is the ray of the row i of SHV, that is of the rows 2i and
  2i+1 of EXT.
 */
template <typename T, typename Tint, typename Tgroup>
Tgroup get_ray_group(PerfectFormInfoForComplex<T, Tint, Tgroup> const &cone,
                     std::vector<int> const &ray_of_vect,
                     std::vector<std::vector<int>> const &ray_vectors) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n_ray = ray_vectors.size();
  std::vector<Telt> l_gens;
  for (auto &ePerm : cone.GRP_ext.GeneratorsOfGroup()) {
    std::vector<Tidx> v(n_ray);
    for (int i_ray = 0; i_ray < n_ray; i_ray++) {
      int i_vect = ray_vectors[i_ray][0];
      v[i_ray] = ray_of_vect[ePerm.at(2 * i_vect) / 2];
#ifdef SANITY_CHECK_PERFECT_LOW_DIM
      for (int j_vect : ray_vectors[i_ray]) {
        for (int i_sign = 0; i_sign < 2; i_sign++) {
          int j_img = ePerm.at(2 * j_vect + i_sign) / 2;
          if (ray_of_vect[j_img] != static_cast<int>(v[i_ray])) {
            std::cerr << "PERFLOWDIM: The element does not map a ray to a ray\n";
            throw TerminalException{1};
          }
        }
      }
#endif
    }
    l_gens.push_back(Telt(v));
  }
  return Tgroup(l_gens, n_ray);
}

/*
  The orbits of faces of rank at most max_rank of one perfect domain.
 */
template <typename T, typename Tint, typename Tgroup>
PerfectConeLowRankFaces<T, Tint, Tgroup>
get_cone_low_rank_faces(LinSpaceMatrix<T> const &LinSpa,
                        PerfectFormInfoForComplex<T, Tint, Tgroup> const &cone,
                        int const &max_rank, std::ostream &os) {
  MyMatrix<Tint> const &EXT = cone.EXT;
  int n_ext = EXT.rows();
  int n_vect = n_ext / 2;
  int dim = EXT.cols();
#ifdef SANITY_CHECK_PERFECT_LOW_DIM
  if (n_ext != 2 * n_vect) {
    std::cerr << "PERFLOWDIM: EXT should have an even number of rows\n";
    throw TerminalException{1};
  }
  for (int i_vect = 0; i_vect < n_vect; i_vect++) {
    for (int i = 0; i < dim; i++) {
      if (EXT(2 * i_vect + 1, i) != -EXT(2 * i_vect, i)) {
        std::cerr << "PERFLOWDIM: The rows 2i and 2i+1 of EXT should be v and -v\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  MyMatrix<Tint> SHV(n_vect, dim);
  for (int i_vect = 0; i_vect < n_vect; i_vect++) {
    for (int i = 0; i < dim; i++) {
      SHV(i_vect, i) = EXT(2 * i_vect, i);
    }
  }
  MyMatrix<T> ScalMat = get_scal_mat<T, Tint>(LinSpa.ListMat, SHV);
  std::map<MyVector<T>, int> map_ray;
  std::vector<int> ray_of_vect(n_vect);
  std::vector<std::vector<int>> ray_vectors;
  std::vector<MyVector<T>> l_ray;
  for (int i_vect = 0; i_vect < n_vect; i_vect++) {
    MyVector<T> V = GetMatrixRow(ScalMat, i_vect);
    auto iter = map_ray.find(V);
    if (iter == map_ray.end()) {
      int i_ray = l_ray.size();
      map_ray[V] = i_ray;
      l_ray.push_back(V);
      ray_vectors.push_back({i_vect});
      ray_of_vect[i_vect] = i_ray;
    } else {
      ray_vectors[iter->second].push_back(i_vect);
      ray_of_vect[i_vect] = iter->second;
    }
  }
  MyMatrix<T> FAC = MatrixFromVectorFamily(l_ray);
  int dimEXT = LinSpa.ListMat.size();
  if (RankMat(FAC) != dimEXT) {
    std::cerr << "PERFLOWDIM: The rays of a perfect domain should span the T-space\n";
    throw TerminalException{1};
  }
  if (max_rank < 1 || max_rank > dimEXT) {
    std::cerr << "PERFLOWDIM: max_rank=" << max_rank << " should be between 1 and dimEXT=" << dimEXT << "\n";
    throw TerminalException{1};
  }
  Tgroup GRP = get_ray_group(cone, ray_of_vect, ray_vectors);
  // The vertices are only used by the extreme ray methods.
  MyMatrix<T> EXT_unused(0, dimEXT);
  int LevSearch = max_rank - 1;
  std::string method_spann = "LinearProgramming";
  std::string method_final = "all";
  bool ComputeTotalNumberFaces = false;
  std::vector<vectface> levels =
      EnumerationFaces(GRP, FAC, EXT_unused, LevSearch, method_spann,
                       method_final, ComputeTotalNumberFaces, os);
#ifdef DEBUG_PERFECT_LOW_DIM
  os << "PERFLOWDIM: n_vect=" << n_vect << " n_ray=" << l_ray.size() << " |GRP|=" << GRP.size() << " n_orbit=[";
  for (size_t i_lev = 0; i_lev < levels.size(); i_lev++) {
    if (i_lev > 0) {
      os << ",";
    }
    os << levels[i_lev].size();
  }
  os << "]\n";
#endif
  return {std::move(SHV), std::move(FAC), std::move(ray_vectors), std::move(GRP), std::move(levels)};
}

/*
  The equivalence of the well rounded cells.

  The group acts on the vectors by v -> v x^{-1} when it acts on the forms
  by A -> x A x^T (the convention of LINSPA_TestEquivalenceGramMatrix_SHV).
  To a well rounded cell M we attach the form
    R_M = ( sum_{c in C} sum_{v in M} (v c)^T (v c) )^{-1}
  scaled to an integral matrix, which satisfies R_{M x^{-1}} = x R_M x^T for every element x of the group
  that preserves M. C is the set {Id} and ListComm, closed under the
  conjugation c -> o^{-1} c o by the outer elements o. The elements of the
  group commute with ListComm, and an outer equivalence g o conjugates
  ListComm by o, so C is mapped to itself. For the classical space C = {Id}
  and R_M is the inverse of sum v v^T. For a Hermitian space the vectors of
  a cell span K^m but not always Q^{2m}, and the multiplication by the
  matrices of ListComm makes R_M positive definite. If R_M is not positive
  definite for a well rounded cell, the program stops.

  The equivalence is then that of (R_M, SHVdisc) with SHVdisc made of an
  invariant family of R_M that spans the lattice Z^n (the characteristic
  vector set V_cv of InvariantVectorFamily.h), and of the vectors
  +-k v for v in M. The family spans Z^n, so the transformation found is
  integral (with a family that only has full rank it is integral only on
  the lattice spanned by the family, and two cells of different index were
  identified). The multiplier k is the smallest integer k >= 2 that is not
  the gcd of the coordinates of a vector of the family. An integral
  transformation preserves the gcd of the coordinates, so it maps the
  vectors k v of one cell to those of the other: the equivalences of
  (R_M, SHVdisc) are those of the cells. The family is invariant, so k is
  an invariant of the cell.
 */
template <typename T>
std::vector<MyMatrix<T>> get_cell_form_multipliers(LinSpaceMatrix<T> const &LinSpa) {
  std::vector<MyMatrix<T>> l_mult{IdentityMat<T>(LinSpa.n)};
  auto f_insert = [&](MyMatrix<T> const &c) -> void {
    for (auto &c_old : l_mult) {
      if (c_old == c) {
        return;
      }
    }
    l_mult.push_back(c);
  };
  for (auto &c : LinSpa.ListComm) {
    f_insert(c);
  }
  std::vector<MyMatrix<T>> l_outer_inv;
  for (auto &o : LinSpa.l_outer_elements) {
    l_outer_inv.push_back(Inverse(o));
  }
  size_t pos = 0;
  while (pos < l_mult.size()) {
    MyMatrix<T> c = l_mult[pos];
    pos++;
    for (size_t i_outer = 0; i_outer < l_outer_inv.size(); i_outer++) {
      f_insert(l_outer_inv[i_outer] * c * LinSpa.l_outer_elements[i_outer]);
    }
  }
  return l_mult;
}

template <typename T, typename Tint>
struct LowRankCell {
  // The rank of the cone, that is the rank of the q(v) in the T-space.
  int rank;
  // One vector per pair +-v.
  MyMatrix<Tint> SHV;
  // The form R_M.
  MyMatrix<T> R;
  // The multiplier k of the vectors of the cell in SHVdisc_T.
  int mult;
  // The vectors for the equivalence: an invariant family of R_M spanning
  // Z^n, then the +-k v.
  MyMatrix<T> SHVdisc_T;
  // The invariant of (R_M, SHVdisc_T), for the splitting into buckets.
  size_t invariant;
  // Generators of the stabilizer of the cell, acting on the vectors by
  // v -> v g (set by set_low_rank_cell_stabilizer).
  std::vector<MyMatrix<Tint>> l_gens;
  // Whether every element of the stabilizer preserves the orientation of
  // the cone (set by set_low_rank_cell_stabilizer).
  bool is_orientable;
};

template <typename T, typename Tint>
bool is_well_rounded_cell(LinSpaceMatrix<T> const &LinSpa,
                          MyMatrix<Tint> const &SHV, std::ostream &os) {
  MyMatrix<Tint> EXT = conversion_and_duplication<Tint, Tint>(SHV);
  PerfectBoundednessProperty pbp = initial_bounded_property(LinSpa, EXT, os);
  // get_result stops the program if none of the methods applies (the
  // finite stabilizer method needs the descent from the top).
  return get_result(pbp);
}

template <typename T, typename Tint>
LowRankCell<T, Tint> get_low_rank_cell(LinSpaceMatrix<T> const &LinSpa,
                                       std::vector<MyMatrix<T>> const &l_mult,
                                       MyMatrix<Tint> const &SHV,
                                       int const &rank, std::ostream &os) {
  int n = LinSpa.n;
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  MyMatrix<T> S = ZeroMatrix<T>(n, n);
  for (auto &c : l_mult) {
    MyMatrix<T> W = SHV_T * c;
    S += W.transpose() * W;
  }
  if (!IsPositiveDefinite(S, os)) {
    std::cerr << "PERFLOWDIM: The form of a well rounded cell is not positive definite\n";
    std::cerr << "PERFLOWDIM: This T-space is not supported by the equivalence of the cells\n";
    throw TerminalException{1};
  }
  // The inverse scaled to an integral matrix: the scaling does not depend
  // on the cell in its orbit (for x in GL_n(Z) the entries of x R x^T are
  // integral combinations of those of R and conversely).
  MyMatrix<T> R = RemoveFractionMatrix(Inverse(S));
  // The characteristic vector set V_cv: it generates Z^n and is canonical.
  // The family made of shells (ExtractInvariantVectorFamilyZbasis) is not
  // used: the minimum of R is large, of the order of |M|^{n-1}, and the
  // shells are enumerated one norm at a time from the smallest step.
  bool spanning = true;
  bool antipodal_half = false;
  CharVectSet_Budget budget;
  MyMatrix<Tint> FAM = unfold_opt(
      CharacteristicVectorSetCV<T, Tint>(R, spanning, antipodal_half, budget, os),
      "CharacteristicVectorSetCV has no budget here");
  std::set<Tint> set_gcd;
  for (int i_row = 0; i_row < FAM.rows(); i_row++) {
    Tint g = 0;
    for (int i = 0; i < n; i++) {
      g = GcdPair(g, FAM(i_row, i));
    }
    set_gcd.insert(T_abs(g));
  }
  int mult = 2;
  while (set_gcd.count(Tint(mult)) > 0) {
    mult++;
  }
  std::vector<MyVector<T>> l_disc;
  for (int i_row = 0; i_row < FAM.rows(); i_row++) {
    MyVector<Tint> V = GetMatrixRow(FAM, i_row);
    l_disc.push_back(UniversalVectorConversion<T, Tint>(V));
  }
  for (int i_row = 0; i_row < SHV_T.rows(); i_row++) {
    MyVector<T> V = GetMatrixRow(SHV_T, i_row);
    l_disc.push_back(T(mult) * V);
    l_disc.push_back(T(-mult) * V);
  }
  MyMatrix<T> SHVdisc_T = MatrixFromVectorFamily(l_disc);
  size_t seed = 1729;
  std::optional<MyMatrix<T>> CommonGramMat;
  size_t invariant = LINSPA_Invariant_SHV<T>(seed, LinSpa, R, SHVdisc_T, CommonGramMat, os);
  return {rank, SHV, std::move(R), mult, std::move(SHVdisc_T), invariant, {}, false};
}

/*
  The set of lines of the rows of a matrix: each row multiplied by the
  sign of its first nonzero entry.
 */
template <typename T>
std::set<MyVector<T>> get_line_set(MyMatrix<T> const &M) {
  std::set<MyVector<T>> set;
  for (int i_row = 0; i_row < M.rows(); i_row++) {
    MyVector<T> V = GetMatrixRow(M, i_row);
    for (int i = 0; i < V.size(); i++) {
      if (V(i) != 0) {
        if (V(i) < 0) {
          V = -V;
        }
        break;
      }
    }
    set.insert(V);
  }
  return set;
}

/*
  An element x of the group with SHV1 x^{-1} = SHV2 up to the signs, if
  one exists.
 */
template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<T>>
test_equivalence_low_rank_cell(LinSpaceMatrix<T> const &LinSpa,
                               LowRankCell<T, Tint> const &cell1,
                               LowRankCell<T, Tint> const &cell2,
                               std::ostream &os) {
  if (cell1.rank != cell2.rank || cell1.SHV.rows() != cell2.SHV.rows() ||
      cell1.mult != cell2.mult || cell1.invariant != cell2.invariant) {
    return {};
  }
  std::optional<MyMatrix<T>> CommonGramMat;
  std::optional<MyMatrix<T>> opt =
      LINSPA_TestEquivalenceGramMatrix_SHV<T, Tgroup>(
          LinSpa, cell1.R, cell2.R, cell1.SHVdisc_T, cell2.SHVdisc_T,
          CommonGramMat, os);
#ifdef SANITY_CHECK_PERFECT_LOW_DIM
  if (opt) {
    if (!IsIntegralMatrix(*opt) || T_abs(DeterminantMat(*opt)) != 1) {
      std::cerr << "PERFLOWDIM: The equivalence should be an element of GL_n(Z)\n";
      throw TerminalException{1};
    }
    MyMatrix<T> x_inv = Inverse(*opt);
    MyMatrix<T> SHV1_T = UniversalMatrixConversion<T, Tint>(cell1.SHV);
    MyMatrix<T> SHV2_T = UniversalMatrixConversion<T, Tint>(cell2.SHV);
    MyMatrix<T> SHV1_img = SHV1_T * x_inv;
    if (get_line_set(SHV1_img) != get_line_set(SHV2_T)) {
      std::cerr << "PERFLOWDIM: The equivalence does not map the cell to the other\n";
      throw TerminalException{1};
    }
  }
#endif
  return opt;
}

/*
  The stabilizer of a cell and its orientability.

  The stabilizer of (R_M, SHVdisc) is that of the cell: an element x with
  x R_M x^T = R_M acts on the vectors by v -> v x^{-1}, so the generators
  kept are the g = x^{-1}, in the convention v -> v g of perfect_complex.h.
  The cell is orientable if every generator preserves the orientation of
  the span of the q(v) in the T-space (is_face_orientable of
  perfect_complex.h, which needs the vectors with both signs).
 */
template <typename T, typename Tint, typename Tgroup>
void set_low_rank_cell_stabilizer(LinSpaceMatrix<T> const &LinSpa,
                                  LowRankCell<T, Tint> &cell,
                                  std::ostream &os) {
  std::optional<MyMatrix<T>> CommonGramMat;
  Result_ComputeStabilizer_SHV<T, Tgroup> result =
      LINSPA_ComputeStabilizer_SHV<T, Tgroup>(LinSpa, cell.R, cell.SHVdisc_T,
                                              CommonGramMat, os);
  std::vector<MyMatrix<T>> l_x =
      result.get_list_matrix(cell.SHVdisc_T, cell.R, LinSpa, os);
  MyMatrix<Tint> EXT = conversion_and_duplication<Tint, Tint>(cell.SHV);
  std::vector<MyMatrix<Tint>> l_gens;
  for (auto &x : l_x) {
    MyMatrix<T> g_T = Inverse(x);
#ifdef SANITY_CHECK_PERFECT_LOW_DIM
    if (!IsIntegralMatrix(g_T)) {
      std::cerr << "PERFLOWDIM: A stabilizer element should be integral\n";
      throw TerminalException{1};
    }
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(EXT);
    MyMatrix<T> EXT_img = EXT_T * g_T;
    if (get_line_set(EXT_img) != get_line_set(EXT_T)) {
      std::cerr << "PERFLOWDIM: A stabilizer element should preserve the cell\n";
      throw TerminalException{1};
    }
#endif
    l_gens.push_back(UniversalMatrixConversion<Tint, T>(g_T));
  }
  OrientationInfo or_info = get_orientation_info(EXT, LinSpa.ListMat);
  cell.is_orientable = is_face_orientable(EXT, LinSpa.ListMat, or_info, l_gens);
  cell.l_gens = std::move(l_gens);
}

/*
  The well rounded cells of rank min_rank to max_rank up to the group: the
  well rounded faces of every perfect domain, reduced by equivalence. The
  candidates are split into buckets by (rank, number of vectors,
  multiplier, invariant), so the equivalence is only tested within a
  bucket.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<LowRankCell<T, Tint>> get_low_rank_well_rounded_cells(
    LinSpaceMatrix<T> const &LinSpa,
    std::vector<PerfectFormInfoForComplex<T, Tint, Tgroup>> const &l_perfect,
    int const &min_rank, int const &max_rank, std::ostream &os) {
  std::vector<MyMatrix<T>> l_mult = get_cell_form_multipliers(LinSpa);
  std::vector<LowRankCell<T, Tint>> l_cell;
  std::map<std::tuple<int, int, int, size_t>, std::vector<size_t>> buckets;
  for (size_t i_cone = 0; i_cone < l_perfect.size(); i_cone++) {
    PerfectConeLowRankFaces<T, Tint, Tgroup> cf =
        get_cone_low_rank_faces(LinSpa, l_perfect[i_cone], max_rank, os);
    for (int rank = min_rank; rank <= max_rank; rank++) {
      for (Face const &f : cf.levels[rank - 1]) {
        MyMatrix<Tint> SHV = cf.get_face_vectors(f);
        if (!is_well_rounded_cell(LinSpa, SHV, os)) {
          continue;
        }
        LowRankCell<T, Tint> cell = get_low_rank_cell(LinSpa, l_mult, SHV, rank, os);
        std::tuple<int, int, int, size_t> key{rank, static_cast<int>(SHV.rows()), cell.mult, cell.invariant};
        std::vector<size_t> &bucket = buckets[key];
        bool is_new = true;
        for (size_t idx : bucket) {
          std::optional<MyMatrix<T>> opt =
              test_equivalence_low_rank_cell<T, Tint, Tgroup>(LinSpa, cell, l_cell[idx], os);
          if (opt) {
            is_new = false;
            break;
          }
        }
        if (is_new) {
          bucket.push_back(l_cell.size());
          l_cell.emplace_back(std::move(cell));
        }
      }
    }
#ifdef DEBUG_PERFECT_LOW_DIM
    os << "PERFLOWDIM: i_cone=" << i_cone << " |l_cell|=" << l_cell.size() << "\n";
#endif
  }
  return l_cell;
}

// clang-format off
#endif  // SRC_PERFECT_PERFECT_LOW_DIM_H_
// clang-format on

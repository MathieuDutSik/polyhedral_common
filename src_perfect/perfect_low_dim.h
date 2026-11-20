// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_PERFECT_LOW_DIM_H_
#define SRC_PERFECT_PERFECT_LOW_DIM_H_

// clang-format off
#include "perfect_complex.h"
#include "POLY_Kskeletton.h"
#include <map>
#include <string>
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
  + Reducing the faces of the requested ranks up to the action of the
    group (to be added).

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

// clang-format off
#endif  // SRC_PERFECT_PERFECT_LOW_DIM_H_
// clang-format on

// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "siegel_voronoi.h"
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>
// clang-format on

/*
  The description file (see README.md) of the projective system of the
  group Sp(J_N, Z), J_N = [[0, D], [-D, 0]], D = diag(1, N), which is the
  paramodular group K(N) (N = 1 gives Sp_4(Z)): the representative
  configurations of the cells of the Siegel space cut out by the Voronoi
  decomposition of the positive definite quaternary forms.

  The top cells are found by a walk: from the perfect cone containing a
  random rational point of the Siegel space, across the facets of the
  perfect cones that meet the Siegel space in a ridge (a cell of dimension
  5), the cones being taken up to Sp(J_N, Z). The cells of the lower
  levels are the faces of the top cells meeting the Siegel space, and the
  faces at infinity: the faces whose vectors span an isotropic sublattice
  (a line or a Lagrangian plane), of the dimension of their cone minus
  one. The level index of a cell of dimension d is 6 - d.
 */

namespace {

using namespace psys;

// The ring Z^{2g} with the form J = [[0, D], [-D, 0]], D = diag(1, ..., 1, N).
ProjectiveRing paramodular_ring(int const &g, int const &N) {
  int n = 2 * g;
  ProjectiveRing ring;
  ring.d = 1;
  ring.m = n;
  MyMatrix<Tint> B = ZeroMatrix<Tint>(n, n);
  for (int i = 0; i < g; i++) {
    int d_i = (i == g - 1) ? N : 1;
    B(i, g + i) = d_i;
    B(g + i, i) = -d_i;
  }
  ring.B = B;
  MyMatrix<Tint> id = IdentityMat<Tint>(n);
  ring.l_ring_basis.push_back(id);
  ring.l_unit = {id, MyMatrix<Tint>(-id)};
  return ring;
}

// A rational point of the Siegel space of J: g^T A_0 g for g the Cayley
// transform of a random element of the Lie algebra of Sp(J), A_0 the
// diagonal form of entries N / D_i, N / D_i.
MyMatrix<T> random_siegel_point(int const &N, MyMatrix<Tint> const &J, std::mt19937 &gen) {
  int n = J.rows();
  int g = n / 2;
  std::uniform_int_distribution<int> dist_num(-4, 4);
  std::uniform_int_distribution<int> dist_den(1, 4);
  MyMatrix<T> K(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      T x = T(dist_num(gen)) / T(dist_den(gen));
      K(i, j) = x;
      K(j, i) = x;
    }
  }
  MyMatrix<T> J_T = UniversalMatrixConversion<T, Tint>(J);
  // X = K J^{-1} satisfies X J + J X^T = 0.
  MyMatrix<T> X = K * Inverse(J_T);
  MyMatrix<T> I = IdentityMat<T>(n);
  MyMatrix<T> ImX = I - X;
  MyMatrix<T> gmat = (I + X) * Inverse(ImX);
  MyMatrix<T> A0 = ZeroMatrix<T>(n, n);
  for (int i = 0; i < g; i++) {
    T val = T(N) / T(J(i, g + i));
    A0(i, i) = val;
    A0(g + i, g + i) = val;
  }
  MyMatrix<T> A = gmat.transpose() * A0 * gmat;
  MyMatrix<T> AJ = A * J_T;
  MyMatrix<T> C = AJ * AJ;
  if (C != C(0, 0) * I || C(0, 0) >= 0) {
    std::cerr << "SIEGEL: the point is not in the Siegel space\n";
    throw TerminalException{1};
  }
  return A;
}

/*
  Whether the form A is in the interior of the perfect cone of the rows of
  EXT_cone, of facets l_facet: the functional Tr(L X) of each facet (L
  symmetric, vanishing on the v^T v of the facet, positive on the other
  rows) is positive on A.
 */
bool in_cone_interior(MyMatrix<T> const &A, MyMatrix<Tint> const &EXT_cone, std::vector<Face> const &l_facet) {
  int n = A.rows();
  int dim_sym = n * (n + 1) / 2;
  auto f_coeff = [&](MyVector<Tint> const &v) -> MyVector<T> {
    MyVector<T> c(dim_sym);
    int pos = 0;
    for (int a = 0; a < n; a++) {
      for (int b = a; b < n; b++) {
        c(pos) = T(v(a) * v(b)) * (a == b ? 1 : 2);
        pos++;
      }
    }
    return c;
  };
  int n_row = EXT_cone.rows();
  for (auto &f : l_facet) {
    MyMatrix<T> Mrow(f.count(), dim_sym);
    int i_f = 0;
    for (int i = 0; i < n_row; i++) {
      if (f[i] == 1) {
        AssignMatrixRow(Mrow, i_f, f_coeff(GetMatrixRow(EXT_cone, i)));
        i_f++;
      }
    }
    MyMatrix<T> NSP = NullspaceTrMat(Mrow);
    if (NSP.rows() != 1) {
      std::cerr << "SIEGEL: the facet should have a nullspace of dimension 1\n";
      throw TerminalException{1};
    }
    MyVector<T> L = GetMatrixRow(NSP, 0);
    for (int i = 0; i < n_row; i++) {
      if (f[i] == 0) {
        if (L.dot(f_coeff(GetMatrixRow(EXT_cone, i))) < 0) {
          L = -L;
        }
        break;
      }
    }
    T val(0);
    int pos = 0;
    for (int a = 0; a < n; a++) {
      for (int b = a; b < n; b++) {
        val += L(pos) * A(a, b) * (a == b ? 1 : 2);
        pos++;
      }
    }
    if (val <= 0) {
      return false;
    }
  }
  return true;
}

/*
  The rank of the forms v^T v of the rows of EXT selected by f (the
  dimension of the span of the face of the cone), modulo the prime p.
 */
int rank_forms_face_mod(MyMatrix<Tint> const &EXT, Face const &f, int64_t const &p) {
  int n = EXT.cols();
  int dim_sym = n * (n + 1) / 2;
  std::vector<std::vector<int64_t>> rows;
  for (int i = 0; i < EXT.rows(); i++) {
    if (f[i] == 0) {
      continue;
    }
    std::vector<int64_t> row;
    for (int a = 0; a < n; a++) {
      for (int b = a; b < n; b++) {
        int64_t x = UniversalScalarConversion<int64_t, Tint>(EXT(i, a) * EXT(i, b)) % p;
        row.push_back(x < 0 ? x + p : x);
      }
    }
    rows.push_back(row);
  }
  int rnk = 0;
  for (int col = 0; col < dim_sym && rnk < static_cast<int>(rows.size()); col++) {
    int piv = -1;
    for (int r = rnk; r < static_cast<int>(rows.size()); r++) {
      if (rows[r][col] != 0) {
        piv = r;
        break;
      }
    }
    if (piv == -1) {
      continue;
    }
    std::swap(rows[rnk], rows[piv]);
    // the inverse of the pivot by Fermat
    int64_t inv = 1, base = rows[rnk][col], e = p - 2;
    while (e > 0) {
      if (e & 1) {
        inv = static_cast<int64_t>((static_cast<__int128>(inv) * base) % p);
      }
      base = static_cast<int64_t>((static_cast<__int128>(base) * base) % p);
      e >>= 1;
    }
    for (int r = rnk + 1; r < static_cast<int>(rows.size()); r++) {
      if (rows[r][col] == 0) {
        continue;
      }
      int64_t c = static_cast<int64_t>((static_cast<__int128>(rows[r][col]) * inv) % p);
      for (int k = col; k < dim_sym; k++) {
        int64_t v = (rows[r][k] - static_cast<int64_t>((static_cast<__int128>(c) * rows[rnk][k]) % p)) % p;
        rows[r][k] = v < 0 ? v + p : v;
      }
    }
    rnk++;
  }
  return rnk;
}

/*
  The rank of the forms of the face f: modulo the primes 2^31 - 1 and
  2^31 - 19 (the rank modulo a prime is at most the rank over Q, and
  smaller only if the prime divides all the minors of that size), exactly
  when the two disagree or are below the expected rank only for one.
 */
int rank_forms_face(MyMatrix<Tint> const &EXT, Face const &f, [[maybe_unused]] int const &expected) {
  int n = EXT.cols();
  int dim_sym = n * (n + 1) / 2;
  int r1 = rank_forms_face_mod(EXT, f, 2147483647);
  int r2 = rank_forms_face_mod(EXT, f, 2147483629);
  if (r1 == r2) {
    return r1;
  }
  MyMatrix<T> M(f.count(), dim_sym);
  int pos = 0;
  for (int i = 0; i < EXT.rows(); i++) {
    if (f[i] == 0) {
      continue;
    }
    int k = 0;
    for (int a = 0; a < n; a++) {
      for (int b = a; b < n; b++) {
        M(pos, k++) = T(EXT(i, a) * EXT(i, b));
      }
    }
    pos++;
  }
  return RankMat(M);
}

/*
  The stabilizer of a top cell (the elements l_elt of G) as permutations of
  the rows of its cone EXT_cone: the rows may contain v and -v (the minimal
  vectors of the perfect form), so the images are matched exactly, and up
  to sign only for a row whose negative is not a row.
 */
std::vector<std::vector<int>> cone_row_permutations(ProjectiveRing const &ring, MyMatrix<Tint> const &EXT_cone,
                                                    std::vector<MyMatrix<Tint>> const &l_elt) {
  int n_ray = EXT_cone.rows();
  std::map<MyVector<Tint>, int> map_row, map_row_unit;
  for (int i = 0; i < n_ray; i++) {
    map_row[GetMatrixRow(EXT_cone, i)] = i;
    map_row_unit[unit_normalize(ring, GetMatrixRow(EXT_cone, i))] = i;
  }
  std::vector<std::vector<int>> l_perm;
  for (auto &elt : l_elt) {
    std::vector<int> perm(n_ray);
    std::vector<int> used(n_ray, 0);
    for (int i = 0; i < n_ray; i++) {
      MyVector<Tint> img = elt.transpose() * GetMatrixRow(EXT_cone, i);
      auto iter = map_row.find(img);
      if (iter == map_row.end()) {
        iter = map_row_unit.find(unit_normalize(ring, img));
        if (iter == map_row_unit.end()) {
          std::cerr << "SIEGEL: the stabilizer does not preserve the cone\n";
          throw TerminalException{1};
        }
      }
      perm[i] = iter->second;
      if (used[perm[i]]++ > 0) {
        std::cerr << "SIEGEL: the action of the stabilizer on the rows is not a permutation\n";
        throw TerminalException{1};
      }
    }
    l_perm.push_back(perm);
  }
  return l_perm;
}

Face permute_face(Face const &f, std::vector<int> const &perm) {
  Face img(f.size());
  for (size_t i = 0; i < f.size(); i++) {
    if (f[i] == 1) {
      img[perm[i]] = 1;
    }
  }
  return img;
}

// Whether f is the smallest face of its orbit under the permutations.
bool is_orbit_rep(Face const &f, std::vector<std::vector<int>> const &l_perm) {
  for (auto &perm : l_perm) {
    if (permute_face(f, perm) < f) {
      return false;
    }
  }
  return true;
}

/*
  The faces of dimension dim_face of the cone of EXT_cone contained in the
  faces of l_level (of dimension dim_face + 1): the f cap F of that
  dimension for F in l_cand. With l_cand all the facets of the cone these
  are all the facets of the faces of l_level; with a subset, the faces that
  are intersections with it.
 */
std::set<Face> next_level_faces(std::vector<Face> const &l_level, std::vector<Face> const &l_cand,
                                MyMatrix<Tint> const &EXT_cone, int const &dim_face) {
  std::set<Face> set_next;
  for (auto &f : l_level) {
    std::set<Face> set_cand;
    for (auto &F : l_cand) {
      Face h = f & F;
      if (h != f && static_cast<int>(h.count()) >= dim_face) {
        set_cand.insert(h);
      }
    }
    for (auto &h : set_cand) {
      if (set_next.count(h) > 0) {
        continue;
      }
      if (rank_forms_face(EXT_cone, h, dim_face) == dim_face) {
        set_next.insert(h);
      }
    }
  }
  return set_next;
}

/*
  An invariant of the configurations up to G (J the form), with an integer
  (the dimension of the cell, or -1). With Q the sum of the v^T v of the
  rows (Q becomes g^T Q g with the configuration), the values v Q^{-1} w
  are invariant under GL_n(Z) and the v J w under G: each vector gets the
  sorted list of the pairs (|v Q^{-1} w|, |v J w|) with the others, and
  the invariant is the sorted list of these profiles. (With the pairings
  v J w alone, many inequivalent cones of Sp_6(Z) have the same invariant,
  and the equivalence tests dominate.)
 */
std::vector<T> config_invariant(MyMatrix<Tint> const &EXT, MyMatrix<Tint> const &J, int const &dim) {
  int k = EXT.rows();
  int n = EXT.cols();
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(EXT);
  MyMatrix<T> Q = ZeroMatrix<T>(n, n);
  for (int i = 0; i < k; i++) {
    MyVector<T> v = GetMatrixRow(EXT_T, i);
    Q += v * v.transpose();
  }
  std::vector<T> V{T(dim), T(k)};
  if (RankMat(Q) < n) {
    // At infinity: the pairings only.
    MyMatrix<Tint> G = EXT * J * EXT.transpose();
    std::vector<std::vector<T>> l_prof(k);
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < k; j++) {
        if (j != i) {
          l_prof[i].push_back(T(abs(G(i, j))));
        }
      }
      std::sort(l_prof[i].begin(), l_prof[i].end());
    }
    std::sort(l_prof.begin(), l_prof.end());
    for (auto &prof : l_prof) {
      V.insert(V.end(), prof.begin(), prof.end());
    }
    return V;
  }
  MyMatrix<T> Gq = EXT_T * Inverse(Q) * EXT_T.transpose();
  MyMatrix<Tint> Gj = EXT * J * EXT.transpose();
  std::vector<std::vector<std::pair<T, T>>> l_prof(k);
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      if (j != i) {
        T q = Gq(i, j);
        l_prof[i].push_back({q < 0 ? T(-q) : q, T(abs(Gj(i, j)))});
      }
    }
    std::sort(l_prof[i].begin(), l_prof[i].end());
  }
  std::sort(l_prof.begin(), l_prof.end());
  for (int i = 0; i < k; i++) {
    V.push_back(Gq(i, i));
  }
  std::sort(V.begin() + 2, V.end());
  for (auto &prof : l_prof) {
    for (auto &pair : prof) {
      V.push_back(pair.first);
      V.push_back(pair.second);
    }
  }
  return V;
}

/*
  An element g of G = Sp(J, Z) with EXT1 g = EXT2 (rows up to order and
  sign), if any, for configurations spanning R^n. The isomorphism of the
  weighted configurations is found by the canonical forms of
  LinPolytope_Isomorphism_GramMat with the Gram matrix Q^{-1} + J (Q the sum
  of the v^T v): an element of G preserves both Q^{-1} and J on the
  vectors, and from v (Q^{-1} + J) w and w (Q^{-1} + J) v both are
  recovered, so no isomorphism means no element of G. The linear map of an
  isomorphism found is checked to be integral and in G; otherwise (an
  isomorphism over Q only) list_embedding decides. The backtracking of
  list_embedding alone, with the pairings v J w only (mostly 0), is far too
  slow on the cones of Sp_6(Z).
 */
std::optional<MyMatrix<Tint>> sp_equivalence(ProjectiveRing const &ring, MyMatrix<Tint> const &EXT1,
                                             MyMatrix<Tint> const &EXT2, std::ostream &os) {
  int n = EXT1.cols();
  if (EXT1.rows() != EXT2.rows()) {
    return {};
  }
  auto f_fallback = [&]() -> std::optional<MyMatrix<Tint>> {
    std::vector<MyMatrix<Tint>> l_M = list_embedding(ring, EXT1, EXT2);
    if (l_M.empty()) {
      return {};
    }
    return l_M[0];
  };
  if (RankMat(EXT1) < n || RankMat(EXT2) < n) {
    return f_fallback();
  }
  MyMatrix<Tint> tot1 = totality(ring, EXT1);
  MyMatrix<Tint> tot2 = totality(ring, EXT2);
  MyMatrix<T> tot1_T = UniversalMatrixConversion<T, Tint>(tot1);
  MyMatrix<T> tot2_T = UniversalMatrixConversion<T, Tint>(tot2);
  MyMatrix<T> J_T = UniversalMatrixConversion<T, Tint>(ring.B);
  auto f_gram = [&](MyMatrix<T> const &tot) -> MyMatrix<T> {
    MyMatrix<T> Q = tot.transpose() * tot;
    return MyMatrix<T>(Inverse(Q) + J_T);
  };
  std::optional<std::vector<Tidx>> opt =
      LinPolytope_Isomorphism_GramMat<T, Tidx>(tot1_T, f_gram(tot1_T), tot2_T, f_gram(tot2_T), os);
  if (!opt) {
    return {};
  }
  std::vector<int> l_basis = basis_rows(tot1);
  MyMatrix<T> B1(n, n), B2(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      B1(i, j) = tot1_T(l_basis[i], j);
      B2(i, j) = tot2_T((*opt)[l_basis[i]], j);
    }
  }
  MyMatrix<T> g_T = Inverse(B1) * B2;
  if (IsIntegralMatrix(g_T) && g_T * J_T * g_T.transpose() == J_T) {
    MyMatrix<Tint> g = UniversalMatrixConversion<Tint, T>(g_T);
    MyMatrix<Tint> img = tot1 * g;
    bool is_ok = true;
    for (int i = 0; i < img.rows(); i++) {
      if (!position_unit(ring, EXT2, GetMatrixRow(img, i))) {
        is_ok = false;
        break;
      }
    }
    if (is_ok) {
      return g;
    }
  }
  return f_fallback();
}

struct TopRep {
  int i_cone;
  MyMatrix<Tint> M;
  MyMatrix<Tint> EXT;
  // An exact point of the Siegel space in the open cone, when the cell was
  // reached across a certified ridge.
  std::optional<MyMatrix<T>> A_interior;
};

}  // namespace

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 7) {
      std::cerr << "PROJ_SiegelVoronoiCells [g] [N] [mode] [n_start] [seed] [output]\n";
      std::cerr << "  g: the genus, the group acts on Z^{2g}\n";
      std::cerr << "  N: J = [[0, D], [-D, 0]], D = diag(1, ..., 1, N) (N = 1: Sp_{2g}(Z), g = 2: the paramodular group K(N))\n";
      std::cerr << "  mode: full (the description file, g = 2 only), top (the top cells only) or rays (the top cells,\n";
      std::cerr << "        the ridges found by shooting n_start rays from an interior point instead of testing all facets)\n";
      std::cerr << "        or fast (as top, the points on the facets searched with the margin 0.01 / |X| only)\n";
      std::cerr << "        or probeC, C >= 1 (the walk as fast, then the faces of the top cells meeting the Siegel space\n";
      std::cerr << "        down to the codimension C, their orbits, the times and the orbifold Euler characteristic)\n";
      std::cerr << "        or ridgesC (the ridge facets of the top cells, then without numerical test the faces of the\n";
      std::cerr << "        cones down to the codimension C that are intersections of ridge facets: bounds on the\n";
      std::cerr << "        number of transverse cells)\n";
      std::cerr << "        or cellsC (the transverse cells down to the codimension C, by orbits under G: the faces of\n";
      std::cerr << "        each orbit of cells up to its stabilizer in G, the orbifold Euler characteristic)\n";
      std::cerr << "        In the modes fast, probeC, ridgesC and cellsC the top cells are saved in output.top after the walk, and\n";
      std::cerr << "        read from it, without walking, when it exists.\n";
      std::cerr << "  n_start: the number of random starts of the numerical search of a point (rays: also the number of rays)\n";
      std::cerr << "  seed: the random seed\n";
      std::cerr << "  output: the file written\n";
      return -1;
    }
    int g = ParseScalar<int>(argv[1]);
    int N = ParseScalar<int>(argv[2]);
    std::string mode = argv[3];
    int n_start = ParseScalar<int>(argv[4]);
    unsigned seed = ParseScalar<unsigned>(argv[5]);
    std::string file_out = argv[6];
    // The modes probeC and ridgesC, C >= 1.
    auto f_mode_codim = [&](std::string const &prefix) -> int {
      size_t len = prefix.size();
      if (mode.size() > len && mode.substr(0, len) == prefix &&
          mode.find_first_not_of("0123456789", len) == std::string::npos) {
        return std::stoi(mode.substr(len));
      }
      return 0;
    };
    int probe_codim = f_mode_codim("probe");
    int ridges_codim = f_mode_codim("ridges");
    int cells_codim = f_mode_codim("cells");
    if (mode != "full" && mode != "top" && mode != "fast" && mode != "rays" && probe_codim == 0 &&
        ridges_codim == 0 && cells_codim == 0) {
      std::cerr << "SIEGEL: the mode should be full, top, fast, rays, probeC, ridgesC or cellsC\n";
      throw TerminalException{1};
    }
    bool use_checkpoint = (mode == "fast" || probe_codim > 0 || ridges_codim > 0 || cells_codim > 0);
    std::string file_top = file_out + ".top";
    if (mode == "full" && g != 2) {
      std::cerr << "SIEGEL: the mode full is for g = 2 only\n";
      throw TerminalException{1};
    }
    int n = 2 * g;
    int dim_top = g * (g + 1);
    std::ostream &os = std::cerr;
    ProjectiveRing ring = paramodular_ring(g, N);
    MyMatrix<Tint> const &J = ring.B;
    // The perfect cones of GL_n(Z) with their adjacencies.
    LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(n);
    PerfectComplexOptions pco{false, false, false, false};
    FullComplexEnumeration<T, Tint, Tgroup> fce = get_full_complex_enumeration_kernel<T, Tint, Tgroup>(LinSpa, pco, os);
    auto const &l_perfect = fce.pctdi.l_perfect;
    std::cerr << "|perfect cones|=" << l_perfect.size() << "\n";
    ConeWalker<T, Tint, Tgroup> walker(fce, os);
    SiegelOracle oracle(J, n_start, seed);
    oracle.set_check_doubtful(mode == "full");
    // The margins of the search on the facets (on the top cells, always
    // all of them: a top cell reached across a certified ridge exists).
    std::vector<double> l_margin_all{0.2, 0.05, 0.01};
    std::vector<double> l_margin_facet = l_margin_all;
    if (mode == "fast" || probe_codim > 0 || ridges_codim > 0 || cells_codim > 0) {
      l_margin_facet = {0.01};
    }
    std::mt19937 gen(seed);
    // The intersections already computed, by configuration.
    std::map<ConfigKey, SiegelIntersection> map_inter;
    auto f_intersect = [&](MyMatrix<Tint> const &EXT_half) -> SiegelIntersection {
      ConfigKey key = config_key(ring, EXT_half);
      auto iter = map_inter.find(key);
      if (iter != map_inter.end()) {
        return iter->second;
      }
      SiegelIntersection si = oracle.intersect(EXT_half);
      map_inter[key] = si;
      return si;
    };
    // The seed: the perfect cone containing a random point in its interior.
    std::vector<TopRep> l_top;
    // The top cells of a previous walk: i_cone, M and EXT for each.
    bool from_checkpoint = false;
    if (use_checkpoint && FILE_IsExistingFile(file_top)) {
      std::ifstream is(file_top);
      int n_top;
      is >> n_top;
      for (int i = 0; i < n_top; i++) {
        int i_cone;
        is >> i_cone;
        MyMatrix<Tint> M = ReadMatrix<Tint>(is);
        MyMatrix<Tint> EXT = ReadMatrix<Tint>(is);
        l_top.push_back({i_cone, M, EXT, {}});
      }
      from_checkpoint = true;
      std::cerr << "SIEGEL: " << n_top << " top cells read from " << file_top << ", no walk\n";
    }
    while (l_top.empty()) {
      MyMatrix<T> A = random_siegel_point(N, J, gen);
      TopPerfectCone<Tint> tpc = walker.find_containing_cone_form(A);
      MyMatrix<Tint> EXT_cone = l_perfect[tpc.i_perfect].EXT * tpc.M;
      MyMatrix<Tint> EXT = rows_up_to_sign(ring, EXT_cone);
      std::vector<Face> l_facet;
      for (auto &adj : l_perfect[tpc.i_perfect].full_adjacencies(os)) {
        l_facet.push_back(adj.f_ext);
      }
      if (in_cone_interior(A, EXT_cone, l_facet)) {
        l_top.push_back({tpc.i_perfect, tpc.M, EXT, {}});
      }
    }
    // The walk through the top cells, up to G, with the exact certificates
    // of the top cells and of the ridges crossed.
    SiegelParametrization sp(J);
    auto f_numeric_point = [&](MyMatrix<Tint> const &EXT_half, std::vector<double> const &weights) -> MyMatrix<double> {
      MyMatrix<double> A = ZeroMatrix<double>(n, n);
      for (int i = 0; i < EXT_half.rows(); i++) {
        for (int a = 0; a < n; a++) {
          for (int b = 0; b < n; b++) {
            A(a, b) += weights[i] * EXT_half(i, a).get_d() * EXT_half(i, b).get_d();
          }
        }
      }
      return A;
    };
    /*
      The facet functionals of the cone l_perfect[i_cone].EXT M: with the
      functional L of a facet of the reference cone as the symmetric matrix
      Lambda (L(A) = Tr(Lambda A)), the points of the cone of EXT M are the
      M^T A M, so the functional becomes M^{-1} Lambda M^{-T}. The reference
      ones are computed once.
     */
    std::map<int, std::vector<MyMatrix<T>>> map_ref_fct;
    auto f_cone_functionals = [&](int const &i_cone, MyMatrix<Tint> const &M) -> std::vector<MyVector<T>> {
      auto iter = map_ref_fct.find(i_cone);
      if (iter == map_ref_fct.end()) {
        std::vector<Face> l_facet;
        for (auto &adj : l_perfect[i_cone].full_adjacencies(os)) {
          l_facet.push_back(adj.f_ext);
        }
        std::vector<MyMatrix<T>> l_mat;
        for (auto &L : facet_functionals(l_perfect[i_cone].EXT, l_facet)) {
          MyMatrix<T> Lambda(n, n);
          int pos = 0;
          for (int a = 0; a < n; a++) {
            for (int b = a; b < n; b++) {
              Lambda(a, b) = L(pos);
              Lambda(b, a) = L(pos);
              pos++;
            }
          }
          l_mat.push_back(Lambda);
        }
        iter = map_ref_fct.emplace(i_cone, std::move(l_mat)).first;
      }
      MyMatrix<T> Minv = Inverse(UniversalMatrixConversion<T, Tint>(M));
      MyMatrix<T> MinvT = Minv.transpose();
      std::vector<MyVector<T>> l_fct;
      for (auto &Lambda : iter->second) {
        MyMatrix<T> Lp = Minv * Lambda * MinvT;
        MyVector<T> L(n * (n + 1) / 2);
        int pos = 0;
        for (int a = 0; a < n; a++) {
          for (int b = a; b < n; b++) {
            L(pos) = Lp(a, b);
            pos++;
          }
        }
        l_fct.push_back(L);
      }
      return l_fct;
    };
    int n_top_certified = 0;
    // The configurations of the certified ridges: in the mode full the cells
    // of dimension g (g + 1) - 1 are exactly the faces equivalent to one of
    // them.
    std::vector<MyMatrix<Tint>> l_ridge_certified;
    int n_ridge = 0;
    int n_ridge_certified = 0;
    size_t pos = from_checkpoint ? l_top.size() : 0;
    while (pos < l_top.size()) {
      TopRep top = l_top[pos];
      std::optional<MyMatrix<T>> opt_A;
      if (top.A_interior) {
        opt_A = top.A_interior;
      } else {
        oracle.set_margins(l_margin_all);
        SiegelIntersection si = f_intersect(top.EXT);
        oracle.set_margins(l_margin_facet);
        if (!si.nonempty || si.dim != dim_top) {
          std::cerr << "SIEGEL: the top cell " << pos << " should meet the Siegel space in dimension " << dim_top << ", dim=" << si.dim << "\n";
          throw TerminalException{1};
        }
        opt_A = certify_top_cell(sp, f_numeric_point(top.EXT, si.weights), f_cone_functionals(top.i_cone, top.M));
      }
      MyMatrix<Tint> EXT_cone = l_perfect[top.i_cone].EXT * top.M;
      std::vector<sing_adj<Tint>> const &l_adj = l_perfect[top.i_cone].full_adjacencies(os);
      std::vector<MyVector<T>> l_fct = f_cone_functionals(top.i_cone, top.M);
      if (opt_A && !in_cone_interior_exact(*opt_A, l_fct)) {
        std::cerr << "SIEGEL: the point of the top cell " << pos << " is not in the open cone\n";
        throw TerminalException{1};
      }
      if (opt_A) {
        n_top_certified++;
      } else {
        std::cerr << "SIEGEL: the top cell " << pos << " is not certified\n";
      }
      HumanTime time_top;
      // The candidate ridges: a facet with a numerical point of the Siegel
      // space on it.
      std::vector<std::pair<int, MyMatrix<double>>> l_cand;
      std::vector<MyMatrix<Tint>> l_stab;
      int n_facet_rep = 0;
      if (mode == "rays") {
        if (!opt_A) {
          std::cerr << "SIEGEL: the rays need a certified interior point\n";
          throw TerminalException{1};
        }
        Eigen::VectorXd v0 = sp.to_vector(sp.param_of(UniversalMatrixConversion<double, T>(*opt_A)));
        for (auto &kv : shoot_rays(sp, v0, l_fct, n_start, gen)) {
          l_cand.push_back(kv);
        }
      } else {
        // One facet per orbit of the stabilizer of the top cell in G: the
        // facet whose key is the smallest of its orbit.
        l_stab = projective_stabilizer(ring, top.EXT, os).l_elt;
        for (size_t i_F = 0; i_F < l_adj.size(); i_F++) {
          MyMatrix<Tint> EXT_facet = rows_up_to_sign(ring, SelectRow(EXT_cone, l_adj[i_F].f_ext));
          ConfigKey key_facet = config_key(ring, EXT_facet);
          bool is_rep = true;
          for (auto &elt : l_stab) {
            if (config_key(ring, EXT_facet * elt) < key_facet) {
              is_rep = false;
              break;
            }
          }
          if (!is_rep) {
            continue;
          }
          n_facet_rep++;
          SiegelIntersection si_f = f_intersect(EXT_facet);
          if (!si_f.nonempty || si_f.dim != dim_top - 1) {
            continue;
          }
          l_cand.push_back({static_cast<int>(i_F), f_numeric_point(EXT_facet, si_f.weights)});
        }
      }
      for (auto &cand : l_cand) {
        int i_F = cand.first;
        sing_adj<Tint> const &adj = l_adj[i_F];
        MyMatrix<Tint> EXT_facet = rows_up_to_sign(ring, SelectRow(EXT_cone, adj.f_ext));
        ConfigKey key_facet = config_key(ring, EXT_facet);
        MyMatrix<Tint> M_adj = adj.eMat * top.M;
        MyMatrix<Tint> EXT_adj = rows_up_to_sign(ring, l_perfect[adj.jCone].EXT * M_adj);
        // The certificate of the ridge: the same facet in the adjacent cone.
        n_ridge++;
        MyMatrix<Tint> EXT_cone_adj = l_perfect[adj.jCone].EXT * M_adj;
        std::vector<sing_adj<Tint>> const &l_adj1 = l_perfect[adj.jCone].full_adjacencies(os);
        int i_F1 = -1;
        for (size_t j = 0; j < l_adj1.size(); j++) {
          if (config_key(ring, rows_up_to_sign(ring, SelectRow(EXT_cone_adj, l_adj1[j].f_ext))) == key_facet) {
            i_F1 = j;
          }
        }
        if (i_F1 == -1) {
          std::cerr << "SIEGEL: the facet is not a facet of the adjacent cone\n";
          throw TerminalException{1};
        }
        std::vector<MyVector<T>> l_fct1 = f_cone_functionals(adj.jCone, M_adj);
        std::optional<MyMatrix<T>> opt_A_adj = certify_ridge(sp, cand.second, l_fct, i_F, l_fct1, i_F1);
        if (opt_A_adj) {
          n_ridge_certified++;
          l_ridge_certified.push_back(EXT_facet);
        } else {
          // Only the certified ridges are crossed, so that every top cell
          // found is certified.
          std::cerr << "SIEGEL: ridge not certified, top cell " << pos << " EXT_facet=" << StringMatrixGAP(EXT_facet) << "\n";
          continue;
        }
        bool is_new = true;
        for (auto &top2 : l_top) {
          if (are_equivalent_configurations(ring, EXT_adj, top2.EXT)) {
            is_new = false;
            break;
          }
        }
        if (is_new) {
          l_top.push_back({static_cast<int>(adj.jCone), M_adj, EXT_adj, opt_A_adj});
          std::cerr << "SIEGEL: new top cell |EXT|=" << EXT_adj.rows() << " |l_top|=" << l_top.size() << "\n";
        }
      }
      std::cerr << "SIEGEL: top cell " << pos << " |EXT|=" << top.EXT.rows() << " |Stab|=" << l_stab.size()
                << " |facets|=" << l_adj.size() << " |facet orbits|=" << n_facet_rep << " |candidate ridges|=" << l_cand.size()
                << " time=" << time_top << "\n";
      pos++;
    }
    std::cerr << "|top cells|=" << l_top.size() << " certified=" << n_top_certified << " |ridges crossed|=" << n_ridge
              << " certified=" << n_ridge_certified << "\n";
    if (use_checkpoint && !from_checkpoint) {
      std::string file_tmp = file_top + ".tmp";
      {
        std::ofstream os_top(file_tmp);
        os_top << l_top.size() << "\n";
        for (auto &top : l_top) {
          os_top << top.i_cone << "\n";
          WriteMatrix(os_top, top.M);
          WriteMatrix(os_top, top.EXT);
        }
      }
      std::rename(file_tmp.c_str(), file_top.c_str());
      std::cerr << "SIEGEL: the top cells written to " << file_top << "\n";
    }
    if (cells_codim > 0) {
      /*
        The cells by orbits under G, level by level (the codimension in the
        cone). An orbit of cells Y of codimension k is represented in one top
        cell, as a face of its cone; its faces of codimension k + 1 are
        taken up to the stabilizer of Y in G (not the one of the top cell),
        and a face Z is a candidate only if the faces of codimension k of the
        cone containing Z are all cells (by their invariant, a necessary
        condition): near a transverse cell Z the Siegel space meets all of
        them. The candidates are taken up to G (each orbit is tested once)
        and tested numerically from the point of Y, then by the full
        search. The cells found are the transverse cells (the non
        transverse ones, deeper in the cone than their dimension, are
        generally missed).
       */
      SiegelOracle oracle_warm(J, 2, seed);
      oracle_warm.set_check_doubtful(false);
      oracle_warm.set_margins({0.01});
      // For the facets of the top cells the filter is empty (their only
      // cover is the top cell) and most of them miss the Siegel space: no
      // full search (0.5 s each), more random starts instead, and a negative
      // result is not cached (a ridge missed from one top cell is found from
      // the other one).
      SiegelOracle oracle_warm_facet(J, 6, seed + 1);
      oracle_warm_facet.set_check_doubtful(false);
      oracle_warm_facet.set_margins({0.01});
      int dim_sym = n * (n + 1) / 2;
      struct CellRep {
        int i_top;
        Face f;  // the rows of the cone of the top cell
        MyMatrix<Tint> EXT;  // rows_up_to_sign of these rows
        std::vector<double> weights;  // a point of the cell, on the rows of EXT
        int dim;
        int64_t stab_order;
      };
      struct TopData {
        MyMatrix<Tint> EXT_cone;
        std::vector<Face> l_facet;
        // For each facet, the configuration of the perfect cone on the other
        // side of it.
        std::map<Face, MyMatrix<Tint>> map_neighbor;
      };
      std::vector<TopData> l_top_data;
      for (auto &top : l_top) {
        TopData td;
        td.EXT_cone = l_perfect[top.i_cone].EXT * top.M;
        for (auto &adj : l_perfect[top.i_cone].full_adjacencies(os)) {
          td.l_facet.push_back(adj.f_ext);
          MyMatrix<Tint> M_adj = adj.eMat * top.M;
          td.map_neighbor[adj.f_ext] = rows_up_to_sign(ring, l_perfect[adj.jCone].EXT * M_adj);
        }
        l_top_data.push_back(std::move(td));
      }
      // The top cells by invariant: a facet is a ridge only if the Siegel
      // space crosses it, so only if the cone on the other side is a top
      // cell.
      std::map<std::vector<T>, std::vector<int>> map_top_inv;
      for (size_t i_top = 0; i_top < l_top.size(); i_top++) {
        map_top_inv[config_invariant(l_top[i_top].EXT, J, -1)].push_back(static_cast<int>(i_top));
      }
      // The top cell equivalent to EXT, with g such that EXT g is its
      // configuration (up to sign), if any.
      auto f_top_match = [&](MyMatrix<Tint> const &EXT) -> std::optional<std::pair<int, MyMatrix<Tint>>> {
        auto iter = map_top_inv.find(config_invariant(EXT, J, -1));
        if (iter == map_top_inv.end()) {
          return {};
        }
        for (auto &i_top : iter->second) {
          if (EXT.rows() != l_top[i_top].EXT.rows()) {
            continue;
          }
          std::optional<MyMatrix<Tint>> opt_g = sp_equivalence(ring, EXT, l_top[i_top].EXT, os);
          if (opt_g) {
            return std::pair<int, MyMatrix<Tint>>{i_top, *opt_g};
          }
        }
        return {};
      };
      // The orbits of the configurations tested, with the result: by level
      // and invariant, the configurations and their intersections.
      std::vector<std::map<std::vector<T>, std::vector<std::pair<MyMatrix<Tint>, SiegelIntersection>>>>
          l_map_tested(cells_codim + 1);
      // The invariants of the cells of each level (for the filter).
      std::vector<std::set<std::vector<T>>> l_set_inv_cell(cells_codim + 1);
      std::vector<std::vector<CellRep>> ll_cell(cells_codim + 1);
      std::map<int, std::vector<int64_t>> map_dim_stab;
      auto f_weights_from = [&](CellRep const &Y, MyMatrix<Tint> const &EXT) -> std::vector<double> {
        std::map<MyVector<Tint>, double> map_w;
        for (int i = 0; i < Y.EXT.rows(); i++) {
          map_w[unit_normalize(ring, GetMatrixRow(Y.EXT, i))] = Y.weights[i];
        }
        std::vector<double> w;
        for (int i = 0; i < EXT.rows(); i++) {
          w.push_back(map_w[unit_normalize(ring, GetMatrixRow(EXT, i))]);
        }
        return w;
      };
      // Level 0: the top cells, with a point of each when one is found (a
      // thin top cell may have none: its facets then start from the point
      // of the top cell on the other side).
      SiegelOracle oracle_top(J, 50, seed + 2);
      for (size_t i_top = 0; i_top < l_top.size(); i_top++) {
        std::optional<std::vector<double>> opt_w = oracle_top.point_of_existing_cell(l_top[i_top].EXT);
        std::vector<double> w;
        if (opt_w) {
          w = *opt_w;
        } else {
          std::cerr << "CELLS: no point found in the top cell " << i_top << " (|EXT|=" << l_top[i_top].EXT.rows()
                    << "), its facets start from their other side\n";
        }
        Face f_all(l_top_data[i_top].EXT_cone.rows());
        f_all.set();
        int64_t order = projective_stabilizer(ring, l_top[i_top].EXT, os).l_elt.size();
        ll_cell[0].push_back({static_cast<int>(i_top), f_all, l_top[i_top].EXT, w, dim_top, order});
        l_set_inv_cell[0].insert(config_invariant(l_top[i_top].EXT, J, -1));
        map_dim_stab[dim_top].push_back(order);
      }
      for (int codim = 1; codim <= cells_codim && !ll_cell[codim - 1].empty(); codim++) {
        HumanTime time_level;
        int dim_face = dim_sym - codim;
        size_t n_face = 0, n_rep = 0, n_filtered = 0, n_tested = 0, n_cached = 0, n_warm = 0, n_full = 0;
        size_t n_equiv = 0;
        for (auto &Y : ll_cell[codim - 1]) {
          TopData const &td = l_top_data[Y.i_top];
          // The faces of Y of codimension one, and the action of its
          // stabilizer in G on the rows of its cone.
          std::set<Face> set_face = next_level_faces({Y.f}, td.l_facet, td.EXT_cone, dim_face);
          n_face += set_face.size();
          std::vector<int> l_pos;
          for (int i = 0; i < td.EXT_cone.rows(); i++) {
            if (Y.f[i] == 1) {
              l_pos.push_back(i);
            }
          }
          MyMatrix<Tint> EXT_Y = SelectRow(td.EXT_cone, Y.f);
          std::vector<std::vector<int>> l_perm_Y =
              cone_row_permutations(ring, EXT_Y, projective_stabilizer(ring, Y.EXT, os).l_elt);
          auto f_restrict = [&](Face const &h) -> Face {
            Face g(l_pos.size());
            for (size_t u = 0; u < l_pos.size(); u++) {
              g[u] = h[l_pos[u]];
            }
            return g;
          };
          // The facets of the cone containing h, for its upper covers.
          for (auto &h : set_face) {
            if (!is_orbit_rep(f_restrict(h), l_perm_Y)) {
              continue;
            }
            if (n_rep % 5000 == 4999) {
              std::cerr << "CELLS: codim=" << codim << " progress: reps=" << n_rep + 1 << " cells="
                        << ll_cell[codim].size() << " tested=" << n_tested << " time=" << time_level << "\n";
            }
            n_rep++;
            MyMatrix<Tint> EXT_h = rows_up_to_sign(ring, SelectRow(td.EXT_cone, h));
            if (RankMat(EXT_h) < n) {
              continue;
            }
            // At the codimension 1: the cone on the other side of the facet
            // must be a top cell.
            std::optional<std::pair<int, MyMatrix<Tint>>> opt_match;
            if (codim == 1) {
              opt_match = f_top_match(td.map_neighbor.at(h));
              if (!opt_match) {
                n_filtered++;
                continue;
              }
            }
            // The filter: the faces of the cone of codimension codim - 1
            // containing h are cells.
            std::vector<Face> l_fac_h;
            for (auto &F : td.l_facet) {
              if (h.is_subset_of(F)) {
                l_fac_h.push_back(F);
              }
            }
            bool pass = true;
            std::set<Face> set_cover;
            for (size_t r = 0; r < h.size() && pass; r++) {
              if (h[r] == 1) {
                continue;
              }
              Face cl(h.size());
              cl.set();
              bool any = false;
              for (auto &F : l_fac_h) {
                if (F[r] == 1) {
                  cl &= F;
                  any = true;
                }
              }
              if (!any || set_cover.count(cl) > 0) {
                continue;
              }
              set_cover.insert(cl);
              if (rank_forms_face(td.EXT_cone, cl, dim_face + 1) != dim_face + 1) {
                continue;
              }
              MyMatrix<Tint> EXT_cl = rows_up_to_sign(ring, SelectRow(td.EXT_cone, cl));
              if (RankMat(EXT_cl) == n && l_set_inv_cell[codim - 1].count(config_invariant(EXT_cl, J, -1)) == 0) {
                pass = false;
              }
            }
            if (!pass) {
              n_filtered++;
              continue;
            }
            // The orbit of h tested once.
            std::vector<T> inv = config_invariant(EXT_h, J, -1);
            std::vector<std::pair<MyMatrix<Tint>, SiegelIntersection>> &l_bucket = l_map_tested[codim][inv];
            bool known = false;
            for (auto &pair : l_bucket) {
              n_equiv++;
              if (sp_equivalence(ring, EXT_h, pair.first, os)) {
                known = true;
                break;
              }
            }
            if (known) {
              n_cached++;
              continue;
            }
            n_tested++;
            SiegelIntersection si;
            if (codim == 1) {
              // The start: the point of Y, else the one of the top cell on
              // the other side, mapped by the equivalence.
              std::optional<std::vector<double>> init;
              if (!Y.weights.empty()) {
                init = f_weights_from(Y, EXT_h);
              } else {
                CellRep const &Y2 = ll_cell[0][opt_match->first];
                if (!Y2.weights.empty()) {
                  // The rows of EXT_h are rows of the cone on the other side,
                  // so their images are rows of that top cell (same order).
                  init = f_weights_from(Y2, EXT_h * opt_match->second);
                }
              }
              si = oracle_warm_facet.intersect(EXT_h, init);
              if (!si.nonempty) {
                continue;
              }
              n_warm++;
            } else {
              si = oracle_warm.intersect(EXT_h, f_weights_from(Y, EXT_h));
              if (si.nonempty) {
                n_warm++;
              } else {
                si = f_intersect(EXT_h);
                n_full++;
              }
            }
            l_bucket.push_back({EXT_h, si});
            if (!si.nonempty) {
              continue;
            }
            int64_t order = projective_stabilizer(ring, EXT_h, os).l_elt.size();
            ll_cell[codim].push_back({Y.i_top, h, EXT_h, si.weights, si.dim, order});
            l_set_inv_cell[codim].insert(inv);
            map_dim_stab[si.dim].push_back(order);
          }
        }
        std::map<int, size_t> map_dim_count;
        for (auto &Z : ll_cell[codim]) {
          map_dim_count[Z.dim]++;
        }
        std::cerr << "CELLS: codim=" << codim << " |orbits of cells|=" << ll_cell[codim].size() << " by dimension:";
        for (auto &kv : map_dim_count) {
          std::cerr << " " << kv.first << ":" << kv.second;
        }
        std::cerr << " faces=" << n_face << " reps=" << n_rep << " filtered=" << n_filtered << " cached=" << n_cached
                  << " tested=" << n_tested << " (warm " << n_warm << ", full " << n_full << ") equivalence tests="
                  << n_equiv << " time=" << time_level << "\n";
      }
      T chi(0);
      for (auto &kv : map_dim_stab) {
        T sum(0);
        for (auto &order : kv.second) {
          sum += T(1) / T(static_cast<long>(order));
        }
        std::cerr << "CELLS: dim=" << kv.first << " |orbits|=" << kv.second.size() << " sum 1/|Stab|=" << sum << "\n";
        chi += (kv.first % 2 == 0) ? sum : T(-sum);
      }
      std::cerr << "CELLS: orbifold Euler characteristic of the cells found=" << chi << "\n";
      std::cerr << "Normal termination of PROJ_SiegelVoronoiCells\n";
      runtime(time);
      return 0;
    }
    if (ridges_codim > 0) {
      /*
        An estimate of the number of transverse cells. A transverse cell of
        codimension c in its top cell is cut on a face of codimension c of
        the cone, and near it the Siegel space meets every facet of the
        cone containing that face in a ridge: so the face is an
        intersection of ridge facets. The ridge facets of each top cell
        are found numerically (one solve per orbit of facets, started from
        a point of the top cell), then the faces of the cone that are
        intersections of ridge facets are enumerated down to the
        codimension ridges_codim, one per orbit of the stabilizer of the
        top cell, without numerical test. Their numbers summed over the top
        cells bound the number of orbits of transverse cells from above
        (a face may not meet the Siegel space, and a cell is counted in
        each top cell containing it); the number of distinct invariants
        (config_invariant) bounds the number of orbits of these faces under
        G from below.
       */
      SiegelOracle oracle_warm(J, 2, seed);
      oracle_warm.set_check_doubtful(false);
      oracle_warm.set_margins({0.01});
      int dim_sym = n * (n + 1) / 2;
      std::vector<size_t> l_n_rep(ridges_codim + 1, 0);
      std::vector<double> l_time(ridges_codim + 1, 0);
      std::vector<std::set<std::vector<T>>> l_set_inv(ridges_codim + 1);
      size_t n_facet_rep_tot = 0;
      double time_ridge = 0;
      for (size_t i_top = 0; i_top < l_top.size(); i_top++) {
        TopRep const &top = l_top[i_top];
        MyMatrix<Tint> EXT_cone = l_perfect[top.i_cone].EXT * top.M;
        std::vector<Face> l_facet;
        for (auto &adj : l_perfect[top.i_cone].full_adjacencies(os)) {
          l_facet.push_back(adj.f_ext);
        }
        std::vector<std::vector<int>> l_perm =
            cone_row_permutations(ring, EXT_cone, projective_stabilizer(ring, top.EXT, os).l_elt);
        // A point of the top cell, its weights indexed by the rows of top.EXT.
        MicrosecondTime time_r;
        oracle.set_margins(l_margin_all);
        SiegelIntersection si_top = f_intersect(top.EXT);
        oracle.set_margins(l_margin_facet);
        if (!si_top.nonempty || si_top.dim != dim_top) {
          std::cerr << "SIEGEL: no point found in the top cell " << i_top << "\n";
          throw TerminalException{1};
        }
        std::map<MyVector<Tint>, double> map_weight;
        for (int i = 0; i < top.EXT.rows(); i++) {
          map_weight[unit_normalize(ring, GetMatrixRow(top.EXT, i))] = si_top.weights[i];
        }
        // The ridge facets: one test per orbit, started from the weights of
        // the point of the top cell on the rows of the facet.
        std::vector<Face> l_ridge;
        int n_rep_facet = 0;
        int n_rep_ridge = 0;
        for (auto &f : l_facet) {
          if (!is_orbit_rep(f, l_perm)) {
            continue;
          }
          n_rep_facet++;
          MyMatrix<Tint> EXT_facet = rows_up_to_sign(ring, SelectRow(EXT_cone, f));
          std::vector<double> init;
          for (int i = 0; i < EXT_facet.rows(); i++) {
            init.push_back(map_weight[unit_normalize(ring, GetMatrixRow(EXT_facet, i))]);
          }
          SiegelIntersection si = oracle_warm.intersect(EXT_facet, init);
          if (si.nonempty && si.dim == dim_top - 1) {
            n_rep_ridge++;
            std::set<Face> set_orb;
            for (auto &perm : l_perm) {
              set_orb.insert(permute_face(f, perm));
            }
            l_ridge.insert(l_ridge.end(), set_orb.begin(), set_orb.end());
          }
        }
        n_facet_rep_tot += n_rep_facet;
        double t_ridge = static_cast<double>(time_r.const_eval_int64()) / 1e6;
        time_ridge += t_ridge;
        std::cerr << "RIDGES: top cell " << i_top << " |EXT|=" << top.EXT.rows() << " facet orbits=" << n_rep_facet
                  << " ridge orbits=" << n_rep_ridge << " ridge facets=" << l_ridge.size() << " time=" << t_ridge
                  << "s\n";
        // The faces that are intersections of ridge facets.
        std::vector<Face> l_level;
        for (auto &f : l_ridge) {
          if (is_orbit_rep(f, l_perm)) {
            l_level.push_back(f);
          }
        }
        for (int codim = 1; codim <= ridges_codim && !l_level.empty(); codim++) {
          MicrosecondTime time_l;
          if (codim > 1) {
            std::vector<Face> l_rep;
            for (auto &f : next_level_faces(l_level, l_ridge, EXT_cone, dim_sym - codim)) {
              if (is_orbit_rep(f, l_perm)) {
                l_rep.push_back(f);
              }
            }
            l_level = l_rep;
          }
          l_n_rep[codim] += l_level.size();
          for (auto &f : l_level) {
            l_set_inv[codim].insert(config_invariant(rows_up_to_sign(ring, SelectRow(EXT_cone, f)), J, -1));
          }
          double t_l = static_cast<double>(time_l.const_eval_int64()) / 1e6;
          l_time[codim] += t_l;
          std::cerr << "RIDGES: top cell " << i_top << " codim=" << codim << " |reps|=" << l_level.size()
                    << " time=" << t_l << "s\n";
        }
      }
      std::cerr << "RIDGES: facet orbits tested=" << n_facet_rep_tot << " time of the ridge tests=" << time_ridge
                << "s\n";
      for (int codim = 1; codim <= ridges_codim; codim++) {
        std::cerr << "RIDGES: codim=" << codim << " (cells of dimension " << dim_top - codim
                  << ") orbits of transverse cells between " << l_set_inv[codim].size() << " and " << l_n_rep[codim]
                  << ", time=" << l_time[codim] << "s\n";
      }
      std::cerr << "Normal termination of PROJ_SiegelVoronoiCells\n";
      runtime(time);
      return 0;
    }
    if (probe_codim > 0) {
      /*
        The faces of the top cells down to the codimension probe_codim (as
        faces of the cones: the facets of a face f of dimension d are the
        f cap F of dimension d - 1, F a facet of the cone), one per orbit
        of the stabilizer of the top cell (acting on the rows of the cone
        by permutations), their intersections with the Siegel space
        (numerical, as in the mode full), and their orbits under G (an
        invariant, then are_equivalent_configurations), with the orders of
        their stabilizers. Only the faces meeting the Siegel space are
        expanded: a cell is in the closure of a cell of a face one level
        up (checked against the mode full for g = 2). Down to the bottom
        this gives all the cells in the Siegel space and the orbifold Euler
        characteristic.
       */
      size_t n_equiv_test = 0;
      size_t max_bucket = 0;
      struct ProbeStat {
        size_t n_face = 0;
        size_t n_rep = 0;
        size_t n_infinity = 0;
        size_t n_empty = 0;
        size_t n_doubtful = 0;
        double time_gen = 0;
        double time_oracle = 0;
        double time_equiv = 0;
      };
      std::vector<ProbeStat> l_stat(probe_codim + 1);
      std::map<std::vector<T>, std::vector<MyMatrix<Tint>>> map_orbit;
      // The orbits of cells by dimension, with the orders of the stabilizers.
      std::map<int, std::vector<int64_t>> map_dim_stab;
      int dim_sym = n * (n + 1) / 2;
      for (auto &top : l_top) {
        map_dim_stab[dim_top].push_back(static_cast<int64_t>(projective_stabilizer(ring, top.EXT, os).l_elt.size()));
      }
      for (size_t i_top = 0; i_top < l_top.size(); i_top++) {
        TopRep const &top = l_top[i_top];
        MyMatrix<Tint> EXT_cone = l_perfect[top.i_cone].EXT * top.M;
        int n_ray = EXT_cone.rows();
        std::vector<Face> l_facet;
        for (auto &adj : l_perfect[top.i_cone].full_adjacencies(os)) {
          l_facet.push_back(adj.f_ext);
        }
        std::vector<std::vector<int>> l_perm =
            cone_row_permutations(ring, EXT_cone, projective_stabilizer(ring, top.EXT, os).l_elt);
        Face f_all(n_ray);
        f_all.set();
        std::vector<Face> l_level{f_all};
        for (int codim = 1; codim <= probe_codim && !l_level.empty(); codim++) {
          ProbeStat &st = l_stat[codim];
          MicrosecondTime time_g;
          int dim_face = dim_sym - codim;
          std::set<Face> set_next = next_level_faces(l_level, l_facet, EXT_cone, dim_face);
          std::vector<Face> l_rep;
          for (auto &f : set_next) {
            if (is_orbit_rep(f, l_perm)) {
              l_rep.push_back(f);
            }
          }
          st.n_face += set_next.size();
          st.time_gen += static_cast<double>(time_g.const_eval_int64());
          std::vector<Face> l_meet;
          for (auto &f : l_rep) {
            MyMatrix<Tint> EXT_face = rows_up_to_sign(ring, SelectRow(EXT_cone, f));
            st.n_rep++;
            if (RankMat(EXT_face) < n) {
              st.n_infinity++;
              continue;
            }
            MicrosecondTime time_o;
            SiegelIntersection si = f_intersect(EXT_face);
            st.time_oracle += static_cast<double>(time_o.const_eval_int64());
            if (si.doubtful) {
              st.n_doubtful++;
            }
            if (!si.nonempty) {
              st.n_empty++;
              continue;
            }
            l_meet.push_back(f);
            MicrosecondTime time_e;
            std::vector<MyMatrix<Tint>> &l_bucket = map_orbit[config_invariant(EXT_face, J, si.dim)];
            max_bucket = std::max(max_bucket, l_bucket.size());
            bool is_new = true;
            for (auto &EXT2 : l_bucket) {
              n_equiv_test++;
              if (are_equivalent_configurations(ring, EXT_face, EXT2)) {
                is_new = false;
                break;
              }
            }
            if (is_new) {
              l_bucket.push_back(EXT_face);
              size_t order = 0;
              try {
                order = projective_stabilizer(ring, EXT_face, os).l_elt.size();
              } catch (TerminalException const &e) {
                std::cerr << "PROBE: the stabilizer failed for EXT_face=" << StringMatrixGAP(EXT_face) << "\n";
                throw;
              }
              map_dim_stab[si.dim].push_back(static_cast<int64_t>(order));
            }
            st.time_equiv += static_cast<double>(time_e.const_eval_int64());
          }
          std::cerr << "PROBE: top cell " << i_top << " |EXT|=" << top.EXT.rows() << " codim=" << codim
                    << " |faces|=" << set_next.size() << " |reps|=" << l_rep.size() << " |meeting|=" << l_meet.size()
                    << " |orbits so far|=" << map_dim_stab.size() << " equivalence tests=" << n_equiv_test
                    << " max bucket=" << max_bucket << "\n";
          // The next level from the representatives meeting the Siegel space
          // (the faces of the other ones are their images).
          l_level = l_meet;
        }
      }
      for (int codim = 1; codim <= probe_codim; codim++) {
        ProbeStat const &st = l_stat[codim];
        std::cerr << "PROBE: codim=" << codim << " |faces|=" << st.n_face << " |reps|=" << st.n_rep
                  << " at infinity=" << st.n_infinity << " empty=" << st.n_empty << " doubtful=" << st.n_doubtful
                  << " time generation=" << st.time_gen / 1e6 << "s oracle=" << st.time_oracle / 1e6
                  << "s equivalence=" << st.time_equiv / 1e6 << "s\n";
      }
      T chi(0);
      for (auto &kv : map_dim_stab) {
        T sum(0);
        for (auto &order : kv.second) {
          sum += T(1) / T(static_cast<long>(order));
        }
        std::cerr << "PROBE: dim=" << kv.first << " |orbits|=" << kv.second.size() << " sum 1/|Stab|=" << sum << "\n";
        chi += (kv.first % 2 == 0) ? sum : T(-sum);
      }
      std::cerr << "PROBE: orbifold Euler characteristic of the cells found=" << chi << "\n";
      std::cerr << "Normal termination of PROJ_SiegelVoronoiCells\n";
      runtime(time);
      return 0;
    }
    if (mode != "full") {
      /*
        An independent check of the walk: random rational points A(X, W)
        of the Siegel space (X with entries in [-3, 3], W = L L^T for L
        lower triangular with a random diagonal in [1/8, 8]), the perfect
        cone containing each, which has to be equivalent to a top cell
        found when the point is in its interior.
       */
      int n_sample = 200;
      int n_interior = 0;
      int n_unknown = 0;
      std::uniform_int_distribution<int> dist_x(-30, 30);
      std::uniform_int_distribution<int> dist_l(-20, 20);
      std::uniform_int_distribution<int> dist_e(-3, 3);
      for (int i_sample = 0; i_sample < n_sample; i_sample++) {
        SiegelParam<T> P{MyMatrix<T>(g, g), MyMatrix<T>(g, g)};
        MyMatrix<T> L = ZeroMatrix<T>(g, g);
        for (int i = 0; i < g; i++) {
          for (int j = 0; j <= i; j++) {
            T x = T(dist_x(gen)) / T(10);
            P.X(i, j) = x;
            P.X(j, i) = x;
            L(i, j) = (i == j) ? T(1 << (dist_e(gen) + 3)) / T(8) : T(dist_l(gen)) / T(10);
          }
        }
        P.W = L * L.transpose();
        MyMatrix<T> A = sp.point(P, false);
        TopPerfectCone<Tint> tpc = walker.find_containing_cone_form(A);
        if (!in_cone_interior_exact(A, f_cone_functionals(tpc.i_perfect, tpc.M))) {
          continue;
        }
        n_interior++;
        MyMatrix<Tint> EXT = rows_up_to_sign(ring, l_perfect[tpc.i_perfect].EXT * tpc.M);
        bool is_known = false;
        for (auto &top : l_top) {
          if (are_equivalent_configurations(ring, EXT, top.EXT)) {
            is_known = true;
            break;
          }
        }
        if (!is_known) {
          n_unknown++;
          std::cerr << "SIEGEL: sample in an unknown top cell |EXT|=" << EXT.rows() << "\n";
        }
      }
      std::cerr << "SIEGEL: samples=" << n_sample << " in the interior of a cone=" << n_interior
                << " in an unknown top cell=" << n_unknown << "\n";
      // The top cells with the orders of their stabilizers in Sp(J, Z).
      std::ofstream os_out(file_out);
      os_out << "# The top cells of Sp(J, Z), J = [[0, D], [-D, 0]], D = diag(1, ..., 1, " << N << "), g = " << g << "\n";
      os_out << "# for each: |EXT|, |Stab|, the perfect cone, the vectors\n";
      os_out << l_top.size() << "\n";
      for (auto &top : l_top) {
        FiniteStabilizer fs = projective_stabilizer(ring, top.EXT, os);
        std::cerr << "top cell |EXT|=" << top.EXT.rows() << " perfect cone=" << top.i_cone << " |Stab|=" << fs.l_elt.size() << "\n";
        os_out << top.EXT.rows() << " " << fs.l_elt.size() << " " << top.i_cone << "\n";
        WriteMatrix(os_out, top.EXT);
      }
      std::cerr << "Normal termination of PROJ_SiegelVoronoiCells\n";
      runtime(time);
      return 0;
    }
    /*
      The cells: the faces of the top cells, up to G, by dimension. The
      intersection with the Siegel space is decided once per orbit of faces
      under G: the numerical oracle, run on two equivalent faces, can give
      them different dimensions (for K(3) a copy of a facet rejected by the
      walk was found of dimension 5, a ridge without a second top cell).
      The facets of dimension 5 are the ones equivalent to a certified
      ridge, the others are decided by the oracle on the first copy met.
     */
    std::vector<std::vector<MyMatrix<Tint>>> ll_cell(7);
    std::vector<std::vector<double>> ll_margin(7);
    int n_ridge_non_facet = 0;
    int n_doubtful = 0;
    auto f_invariant = [&](MyMatrix<Tint> const &EXT) -> std::vector<T> { return config_invariant(EXT, J, -1); };
    std::map<std::vector<T>, std::vector<std::pair<MyMatrix<Tint>, SiegelIntersection>>> map_orbit_si;
    auto f_find_orbit = [&](MyMatrix<Tint> const &EXT) -> std::optional<SiegelIntersection> {
      auto iter = map_orbit_si.find(f_invariant(EXT));
      if (iter == map_orbit_si.end()) {
        return {};
      }
      for (auto &pair : iter->second) {
        if (are_equivalent_configurations(ring, EXT, pair.first)) {
          return pair.second;
        }
      }
      return {};
    };
    for (auto &EXT_ridge : l_ridge_certified) {
      if (!f_find_orbit(EXT_ridge)) {
        SiegelIntersection si{true, dim_top - 1, 1, false, {}};
        map_orbit_si[f_invariant(EXT_ridge)].push_back({EXT_ridge, si});
      }
    }
    int n_facet_orbit_rejected = 0;
    for (auto &top : l_top) {
      MyMatrix<Tint> EXT_cone = l_perfect[top.i_cone].EXT * top.M;
      std::vector<Face> l_facet;
      for (auto &adj : l_perfect[top.i_cone].full_adjacencies(os)) {
        l_facet.push_back(adj.f_ext);
      }
      std::set<Face> set_facet(l_facet.begin(), l_facet.end());
      for (auto &f : all_faces_from_facets(EXT_cone.rows(), l_facet)) {
        MyMatrix<Tint> EXT_face = rows_up_to_sign(ring, SelectRow(EXT_cone, f));
        int dim = -1;
        double margin = 1;
        if (RankMat(EXT_face) == n) {
          std::optional<SiegelIntersection> opt_si = f_find_orbit(EXT_face);
          SiegelIntersection si;
          if (opt_si) {
            si = *opt_si;
          } else {
            si = f_intersect(EXT_face);
            if (si.doubtful) {
              n_doubtful++;
              std::cerr << "SIEGEL: doubtful EXT=" << StringMatrixGAP(EXT_face) << "\n";
            }
            if (set_facet.count(f) > 0 && si.nonempty && si.dim == dim_top - 1) {
              // Not equivalent to a certified ridge: the walk tested an
              // equivalent facet and found no ridge.
              n_facet_orbit_rejected++;
              std::cerr << "SIEGEL: facet of numerical dimension 5 not equivalent to a certified ridge, margin="
                        << si.margin << " EXT=" << StringMatrixGAP(EXT_face) << "\n";
              si.nonempty = false;
            }
            map_orbit_si[f_invariant(EXT_face)].push_back({EXT_face, si});
          }
          if (!si.nonempty) {
            continue;
          }
          dim = si.dim;
          margin = si.margin;
          if (dim == 5 && set_facet.count(f) == 0) {
            n_ridge_non_facet++;
            std::cerr << "SIEGEL: non facet ridge EXT=" << StringMatrixGAP(EXT_face) << " margin=" << margin << "\n";
          }
          if (margin < 1e-3 / EXT_face.rows()) {
            std::cerr << "SIEGEL: small margin dim=" << dim << " margin=" << margin << " EXT=" << StringMatrixGAP(EXT_face) << "\n";
          }
        } else {
          // At infinity: the faces spanning an isotropic sublattice.
          MyMatrix<Tint> G = EXT_face * J * EXT_face.transpose();
          if (!IsZeroMatrix(G)) {
            continue;
          }
          int pos_mat = 0;
          MyMatrix<T> ScalMat(EXT_face.rows(), 10);
          for (int i = 0; i < EXT_face.rows(); i++) {
            pos_mat = 0;
            for (int a = 0; a < 4; a++) {
              for (int b = a; b < 4; b++) {
                ScalMat(i, pos_mat) = T(EXT_face(i, a) * EXT_face(i, b));
                pos_mat++;
              }
            }
          }
          dim = RankMat(ScalMat) - 1;
        }
        if (dim < 0 || dim > 6) {
          std::cerr << "SIEGEL: wrong dimension " << dim << "\n";
          throw TerminalException{1};
        }
        bool is_new = true;
        for (auto &EXT2 : ll_cell[dim]) {
          if (are_equivalent_configurations(ring, EXT_face, EXT2)) {
            is_new = false;
            break;
          }
        }
        if (is_new) {
          ll_cell[dim].push_back(EXT_face);
          ll_margin[dim].push_back(margin);
        }
      }
    }
    if (n_ridge_non_facet > 0) {
      std::cerr << "SIEGEL: " << n_ridge_non_facet << " cells of dimension 5 are not facets of their cone, not supported\n";
      throw TerminalException{1};
    }
    std::cerr << "|doubtful faces|=" << n_doubtful << " |facets of dimension 5 not certified|=" << n_facet_orbit_rejected
              << "\n";
    /*
      The incidences: a configuration Y of the dimension d - 1 contained in
      a configuration X of the dimension d (X spanning R^4) is a face of X
      only if the cell of Y is in the closure of the cell of X. The other
      ones are written in the section non_faces.
     */
    struct NonFace {
      int index;
      int iConf;
      MyMatrix<Tint> EXT;
    };
    std::vector<NonFace> l_non_face;
    for (int dim = 6; dim >= 1; dim--) {
      for (size_t iConf = 0; iConf < ll_cell[dim].size(); iConf++) {
        MyMatrix<Tint> const &EXT = ll_cell[dim][iConf];
        if (RankMat(EXT) < 4) {
          continue;
        }
        int n_face = 0;
        int n_non_face = 0;
        for (auto &EXT_face : ll_cell[dim - 1]) {
          for (auto &l_pos : get_all_embeddings(ring, EXT_face, EXT)) {
            std::vector<int> f_in(EXT.rows(), 0);
            for (auto &i : l_pos) {
              f_in[i] = 1;
            }
            if (oracle.in_closure(EXT, f_in)) {
              n_face++;
            } else {
              n_non_face++;
              l_non_face.push_back({6 - dim, static_cast<int>(iConf), select_rows(EXT, l_pos)});
            }
          }
        }
        std::cerr << "dim=" << dim << " iConf=" << iConf << " |EXT|=" << EXT.rows() << " |faces|=" << n_face
                  << " |non faces|=" << n_non_face << "\n";
      }
    }
    std::ofstream os_out(file_out);
    os_out << "# The projective system of the paramodular group K(" << N << ") = Sp(J_N, Z),\n";
    os_out << "# J_N = [[0, D], [-D, 0]], D = diag(1, " << N << "), computed by PROJ_SiegelVoronoiCells:\n";
    os_out << "# the cells of the Siegel space cut out by the Voronoi decomposition.\n";
    os_out << "ring Z\ndimension 4\nform\n4 4\n";
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        os_out << (j > 0 ? " " : "") << J(i, j);
      }
      os_out << "\n";
    }
    // The orientations from the top cells: the complex is not regular at
    // infinity for N > 1 (see compute_boundaries_top_down).
    os_out << "at_infinity rank\nwalk voronoi\norientation top\nlevels 7\n";
    for (int index = 0; index < 7; index++) {
      int dim = 6 - index;
      os_out << "# index " << index << " (dimension " << dim << ")\n";
      os_out << "level " << ll_cell[dim].size() << "\n";
      std::cerr << "dim=" << dim << " |orbits|=" << ll_cell[dim].size() << " sizes/rank/margin:";
      for (size_t i = 0; i < ll_cell[dim].size(); i++) {
        MyMatrix<Tint> const &EXT = ll_cell[dim][i];
        os_out << EXT.rows() << " 4\n";
        for (int r = 0; r < EXT.rows(); r++) {
          for (int j = 0; j < 4; j++) {
            os_out << (j > 0 ? " " : "") << EXT(r, j);
          }
          os_out << "\n";
        }
        std::cerr << " " << EXT.rows() << "/" << RankMat(EXT) << "/" << ll_margin[dim][i];
      }
      std::cerr << "\n";
    }
    if (!l_non_face.empty()) {
      os_out << "# The sub-configurations that are not faces: index, configuration, vectors\n";
      os_out << "non_faces " << l_non_face.size() << "\n";
      for (auto &nf : l_non_face) {
        os_out << nf.index << " " << nf.iConf << " " << nf.EXT.rows() << " 4\n";
        for (int r = 0; r < nf.EXT.rows(); r++) {
          for (int j = 0; j < 4; j++) {
            os_out << (j > 0 ? " " : "") << nf.EXT(r, j);
          }
          os_out << "\n";
        }
      }
    }
    std::cerr << "|non faces|=" << l_non_face.size() << "\n";
    std::cerr << "Normal termination of PROJ_SiegelVoronoiCells\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of PROJ_SiegelVoronoiCells\n";
    exit(e.eVal);
  }
  runtime(time);
}

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
      std::cerr << "        In the modes fast and probeC the top cells are saved in output.top after the walk, and\n";
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
    int probe_codim = 0;
    if (mode.size() > 5 && mode.substr(0, 5) == "probe" &&
        mode.find_first_not_of("0123456789", 5) == std::string::npos) {
      probe_codim = std::stoi(mode.substr(5));
    }
    if (mode != "full" && mode != "top" && mode != "fast" && mode != "rays" && probe_codim == 0) {
      std::cerr << "SIEGEL: the mode should be full, top, fast, rays or probeC\n";
      throw TerminalException{1};
    }
    bool use_checkpoint = (mode == "fast" || probe_codim > 0);
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
    if (mode == "fast" || probe_codim > 0) {
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
      // An invariant of the faces up to G, with the dimension of their
      // intersection with the Siegel space: the sorted list of the sorted
      // absolute pairings of each vector with the others.
      auto f_invariant = [&](MyMatrix<Tint> const &EXT, int const &dim) -> std::vector<Tint> {
        MyMatrix<Tint> G = EXT * J * EXT.transpose();
        int k = EXT.rows();
        std::vector<std::vector<Tint>> l_prof(k);
        for (int i = 0; i < k; i++) {
          for (int j = 0; j < k; j++) {
            if (j != i) {
              l_prof[i].push_back(abs(G(i, j)));
            }
          }
          std::sort(l_prof[i].begin(), l_prof[i].end());
        }
        std::sort(l_prof.begin(), l_prof.end());
        std::vector<Tint> V{Tint(dim), Tint(k)};
        for (auto &prof : l_prof) {
          V.insert(V.end(), prof.begin(), prof.end());
        }
        return V;
      };
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
      std::map<std::vector<Tint>, std::vector<MyMatrix<Tint>>> map_orbit;
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
        // The stabilizer as permutations of the rows of EXT_cone: the rows
        // may contain v and -v (the minimal vectors of the perfect form),
        // so the images are matched exactly, and up to sign only for a row
        // whose negative is not a row.
        std::map<MyVector<Tint>, int> map_row, map_row_unit;
        for (int i = 0; i < n_ray; i++) {
          map_row[GetMatrixRow(EXT_cone, i)] = i;
          map_row_unit[unit_normalize(ring, GetMatrixRow(EXT_cone, i))] = i;
        }
        std::vector<std::vector<int>> l_perm;
        for (auto &elt : projective_stabilizer(ring, top.EXT, os).l_elt) {
          std::vector<int> perm(n_ray);
          for (int i = 0; i < n_ray; i++) {
            MyVector<Tint> img = elt.transpose() * GetMatrixRow(EXT_cone, i);
            auto iter = map_row.find(img);
            if (iter == map_row.end()) {
              iter = map_row_unit.find(unit_normalize(ring, img));
              if (iter == map_row_unit.end()) {
                std::cerr << "PROBE: the stabilizer does not preserve the cone\n";
                throw TerminalException{1};
              }
            }
            perm[i] = iter->second;
          }
          l_perm.push_back(perm);
        }
        for (auto &perm : l_perm) {
          std::vector<int> sorted = perm;
          std::sort(sorted.begin(), sorted.end());
          for (int i = 0; i < n_ray; i++) {
            if (sorted[i] != i) {
              std::cerr << "PROBE: the action of the stabilizer on the rows is not a permutation\n";
              throw TerminalException{1};
            }
          }
        }
        auto f_is_rep = [&](Face const &f) -> bool {
          for (auto &perm : l_perm) {
            Face img(n_ray);
            for (int i = 0; i < n_ray; i++) {
              if (f[i] == 1) {
                img[perm[i]] = 1;
              }
            }
            if (img < f) {
              return false;
            }
          }
          return true;
        };
        Face f_all(n_ray);
        f_all.set();
        std::vector<Face> l_level{f_all};
        for (int codim = 1; codim <= probe_codim && !l_level.empty(); codim++) {
          ProbeStat &st = l_stat[codim];
          MicrosecondTime time_g;
          int dim_face = dim_sym - codim;
          std::set<Face> set_next;
          for (auto &f : l_level) {
            std::set<Face> set_cand;
            for (auto &F : l_facet) {
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
          std::vector<Face> l_rep;
          for (auto &f : set_next) {
            if (f_is_rep(f)) {
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
            std::vector<MyMatrix<Tint>> &l_bucket = map_orbit[f_invariant(EXT_face, si.dim)];
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
    // An invariant of the configurations up to G: the size and the sorted
    // absolute values of the pairing.
    auto f_invariant = [&](MyMatrix<Tint> const &EXT) -> std::vector<Tint> {
      MyMatrix<Tint> G = EXT * J * EXT.transpose();
      std::vector<Tint> V{Tint(EXT.rows())};
      for (int i = 0; i < EXT.rows(); i++) {
        for (int j = i + 1; j < EXT.rows(); j++) {
          V.push_back(abs(G(i, j)));
        }
      }
      std::sort(V.begin() + 1, V.end());
      return V;
    };
    std::map<std::vector<Tint>, std::vector<std::pair<MyMatrix<Tint>, SiegelIntersection>>> map_orbit_si;
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

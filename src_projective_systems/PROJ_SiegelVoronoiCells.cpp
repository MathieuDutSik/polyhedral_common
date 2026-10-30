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

struct TopRep {
  int i_cone;
  MyMatrix<Tint> M;
  MyMatrix<Tint> EXT;
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
    if (mode != "full" && mode != "top" && mode != "fast" && mode != "rays") {
      std::cerr << "SIEGEL: the mode should be full, top, fast or rays\n";
      throw TerminalException{1};
    }
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
    if (mode == "fast") {
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
        l_top.push_back({tpc.i_perfect, tpc.M, EXT});
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
    auto f_cone_functionals = [&](int const &i_cone, MyMatrix<Tint> const &M) -> std::vector<MyVector<T>> {
      std::vector<Face> l_facet;
      for (auto &adj : l_perfect[i_cone].full_adjacencies(os)) {
        l_facet.push_back(adj.f_ext);
      }
      return facet_functionals(l_perfect[i_cone].EXT * M, l_facet);
    };
    int n_top_certified = 0;
    int n_ridge = 0;
    int n_ridge_certified = 0;
    size_t pos = 0;
    while (pos < l_top.size()) {
      TopRep top = l_top[pos];
      oracle.set_margins(l_margin_all);
      SiegelIntersection si = f_intersect(top.EXT);
      oracle.set_margins(l_margin_facet);
      if (!si.nonempty || si.dim != dim_top) {
        std::cerr << "SIEGEL: the top cell " << pos << " should meet the Siegel space in dimension " << dim_top << ", dim=" << si.dim << "\n";
        throw TerminalException{1};
      }
      MyMatrix<Tint> EXT_cone = l_perfect[top.i_cone].EXT * top.M;
      std::vector<sing_adj<Tint>> const &l_adj = l_perfect[top.i_cone].full_adjacencies(os);
      std::vector<MyVector<T>> l_fct = f_cone_functionals(top.i_cone, top.M);
      std::optional<MyMatrix<T>> opt_A = certify_top_cell(sp, f_numeric_point(top.EXT, si.weights), l_fct);
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
        if (certify_ridge(sp, cand.second, l_fct, i_F, l_fct1, i_F1)) {
          n_ridge_certified++;
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
          l_top.push_back({static_cast<int>(adj.jCone), M_adj, EXT_adj});
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
    // The cells: the faces of the top cells, up to G, by dimension.
    std::vector<std::vector<MyMatrix<Tint>>> ll_cell(7);
    std::vector<std::vector<double>> ll_margin(7);
    int n_ridge_non_facet = 0;
    int n_doubtful = 0;
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
          SiegelIntersection si = f_intersect(EXT_face);
          if (si.doubtful) {
            n_doubtful++;
            std::cerr << "SIEGEL: doubtful EXT=" << StringMatrixGAP(EXT_face) << "\n";
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
    std::cerr << "|doubtful faces|=" << n_doubtful << " (with repetitions)\n";
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
    os_out << "at_infinity rank\nwalk voronoi\nlevels 7\n";
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

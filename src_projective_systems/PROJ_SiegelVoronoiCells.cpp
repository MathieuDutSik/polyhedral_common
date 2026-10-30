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

ProjectiveRing paramodular_ring(int const &N) {
  ProjectiveRing ring;
  ring.d = 1;
  ring.m = 4;
  MyMatrix<Tint> B = ZeroMatrix<Tint>(4, 4);
  B(0, 2) = 1;
  B(2, 0) = -1;
  B(1, 3) = N;
  B(3, 1) = -N;
  ring.B = B;
  MyMatrix<Tint> id = IdentityMat<Tint>(4);
  ring.l_ring_basis.push_back(id);
  ring.l_unit = {id, MyMatrix<Tint>(-id)};
  return ring;
}

// A rational point of the Siegel space of J: g^T A_0 g for g the Cayley
// transform of a random element of the Lie algebra of Sp(J).
MyMatrix<T> random_siegel_point(int const &N, MyMatrix<Tint> const &J, std::mt19937 &gen) {
  std::uniform_int_distribution<int> dist_num(-4, 4);
  std::uniform_int_distribution<int> dist_den(1, 4);
  MyMatrix<T> K(4, 4);
  for (int i = 0; i < 4; i++) {
    for (int j = i; j < 4; j++) {
      T x = T(dist_num(gen)) / T(dist_den(gen));
      K(i, j) = x;
      K(j, i) = x;
    }
  }
  MyMatrix<T> J_T = UniversalMatrixConversion<T, Tint>(J);
  // X = K J^{-1} satisfies X J + J X^T = 0.
  MyMatrix<T> X = K * Inverse(J_T);
  MyMatrix<T> I = IdentityMat<T>(4);
  MyMatrix<T> ImX = I - X;
  MyMatrix<T> g = (I + X) * Inverse(ImX);
  MyMatrix<T> A0 = ZeroMatrix<T>(4, 4);
  A0(0, 0) = N;
  A0(1, 1) = 1;
  A0(2, 2) = N;
  A0(3, 3) = 1;
  MyMatrix<T> A = g.transpose() * A0 * g;
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
  auto f_coeff = [](MyVector<Tint> const &v) -> MyVector<T> {
    MyVector<T> c(10);
    int pos = 0;
    for (int a = 0; a < 4; a++) {
      for (int b = a; b < 4; b++) {
        c(pos) = T(v(a) * v(b)) * (a == b ? 1 : 2);
        pos++;
      }
    }
    return c;
  };
  int n_row = EXT_cone.rows();
  for (auto &f : l_facet) {
    MyMatrix<T> Mrow(f.count(), 10);
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
    for (int a = 0; a < 4; a++) {
      for (int b = a; b < 4; b++) {
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
    if (argc != 5) {
      std::cerr << "PROJ_SiegelVoronoiCells [N] [n_start] [seed] [output]\n";
      std::cerr << "  N: the level of the paramodular group (1 for Sp_4(Z))\n";
      std::cerr << "  n_start: the number of random starts of the numerical search of a point\n";
      std::cerr << "  seed: the random seed\n";
      std::cerr << "  output: the description file written\n";
      return -1;
    }
    int N = ParseScalar<int>(argv[1]);
    int n_start = ParseScalar<int>(argv[2]);
    unsigned seed = ParseScalar<unsigned>(argv[3]);
    std::string file_out = argv[4];
    std::ostream &os = std::cerr;
    ProjectiveRing ring = paramodular_ring(N);
    MyMatrix<Tint> const &J = ring.B;
    // The perfect cones of GL_4(Z) with their adjacencies.
    LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(4);
    PerfectComplexOptions pco{false, false, false, false};
    FullComplexEnumeration<T, Tint, Tgroup> fce = get_full_complex_enumeration_kernel<T, Tint, Tgroup>(LinSpa, pco, os);
    auto const &l_perfect = fce.pctdi.l_perfect;
    std::cerr << "|perfect cones|=" << l_perfect.size() << "\n";
    ConeWalker<T, Tint, Tgroup> walker(fce, os);
    SiegelOracle oracle(J, n_start, seed);
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
    // The walk through the top cells, up to G.
    size_t pos = 0;
    while (pos < l_top.size()) {
      TopRep top = l_top[pos];
      SiegelIntersection si = f_intersect(top.EXT);
      if (!si.nonempty || si.dim != 6) {
        std::cerr << "SIEGEL: the top cell " << pos << " should meet the Siegel space in dimension 6, dim=" << si.dim << "\n";
        throw TerminalException{1};
      }
      MyMatrix<Tint> EXT_cone = l_perfect[top.i_cone].EXT * top.M;
      std::vector<sing_adj<Tint>> const &l_adj = l_perfect[top.i_cone].full_adjacencies(os);
      for (auto &adj : l_adj) {
        MyMatrix<Tint> EXT_facet = rows_up_to_sign(ring, SelectRow(EXT_cone, adj.f_ext));
        SiegelIntersection si_f = f_intersect(EXT_facet);
        if (!si_f.nonempty || si_f.dim != 5) {
          continue;
        }
        MyMatrix<Tint> M_adj = adj.eMat * top.M;
        MyMatrix<Tint> EXT_adj = rows_up_to_sign(ring, l_perfect[adj.jCone].EXT * M_adj);
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
      pos++;
    }
    std::cerr << "|top cells|=" << l_top.size() << "\n";
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
        if (RankMat(EXT_face) == 4) {
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

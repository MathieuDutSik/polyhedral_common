// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_PYVISTA_JSON_H_
#define SRC_ROBUST_COVERING_PYVISTA_JSON_H_

// clang-format off
#include "generalized_polytopes.h"
#include "robust_covering_types.h"
#include <algorithm>
#include <iomanip>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

/*
  Output of the P-polytopes as a list of JSON entries for 3D plotting
  with PyVista (see plot_geometry_index.py). Each entry contains:
  * The v_long vertex (as a point).
  * The parallelepipeds realizing the minimum (as wireframe lines).
  * The generalized polytope (as polygonal faces obtained from
    find_generalized_polytope_boundary, with the vertices of each
    polygon put in cyclic order).
  The lattice coordinates are embedded into R^3 by vectors v_1, v_2, v_3
  with <v_i, v_j> = G_ij obtained from the Cholesky decomposition of the
  Gram matrix G. This conversion to real coordinates is done at the very
  last stage, when writing the JSON.
 */

#ifdef DEBUG
#define DEBUG_ROBUST_PYVISTA
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ROBUST_PYVISTA
#endif

// The lower triangular matrix L with L L^T = G. The rows v_i of L
// then satisfy <v_i, v_j> = G_ij and are the embedding vectors.
inline MyMatrix<double> cholesky_lower(MyMatrix<double> const &G) {
  int n = G.rows();
  MyMatrix<double> L = ZeroMatrix<double>(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      double sum = G(i, j);
      for (int k = 0; k < j; k++) {
        sum -= L(i, k) * L(j, k);
      }
      if (i == j) {
#ifdef SANITY_CHECK_ROBUST_PYVISTA
        if (sum <= 0) {
          std::cerr << "PYVISTA: cholesky_lower: the matrix is not positive definite\n";
          throw TerminalException{1};
        }
#endif
        L(i, i) = sqrt(sum);
      } else {
        L(i, j) = sum / L(j, j);
      }
    }
  }
  return L;
}

// The edges of a parallelepiped given by its 2^dim vertices (the rows
// of M, in arbitrary order). Each edge is returned as a pair of row
// indices of M.
// The method: from a base vertex p0, the differences q - p0 form the
// set {sum_i eps_i g_i, eps in {0,1}^dim} for the generators g_i.
// The generators are exactly the nonzero differences that are not the
// sum of two nonzero differences. Each vertex is then identified by
// its eps vector and two vertices are adjacent if and only if their
// eps vectors differ in exactly one position.
template <typename Tint>
std::vector<std::pair<int, int>>
get_parallelepiped_edges(MyMatrix<Tint> const &M) {
  int n_row = M.rows();
  int dim = M.cols();
#ifdef SANITY_CHECK_ROBUST_PYVISTA
  if (n_row != (1 << dim)) {
    std::cerr << "PYVISTA: get_parallelepiped_edges: n_row=" << n_row
              << " is not 2^dim for dim=" << dim << "\n";
    throw TerminalException{1};
  }
#endif
  MyVector<Tint> p0 = GetMatrixRow(M, 0);
  std::vector<MyVector<Tint>> diffs;
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyVector<Tint> diff = GetMatrixRow(M, i_row) - p0;
    diffs.push_back(diff);
  }
  auto is_zero = [](MyVector<Tint> const &V) -> bool {
    for (int i = 0; i < V.size(); i++) {
      if (V(i) != 0) {
        return false;
      }
    }
    return true;
  };
  // The generators: the nonzero differences that are not the sum of
  // two nonzero differences.
  std::vector<MyVector<Tint>> gens;
  for (int i = 0; i < n_row; i++) {
    if (is_zero(diffs[i])) {
      continue;
    }
    bool is_decomposable = false;
    for (int j = 0; j < n_row && !is_decomposable; j++) {
      if (is_zero(diffs[j])) {
        continue;
      }
      for (int k = 0; k < n_row; k++) {
        if (is_zero(diffs[k])) {
          continue;
        }
        if (diffs[j] + diffs[k] == diffs[i]) {
          is_decomposable = true;
          break;
        }
      }
    }
    if (!is_decomposable) {
      gens.push_back(diffs[i]);
    }
  }
#ifdef SANITY_CHECK_ROBUST_PYVISTA
  if (gens.size() != static_cast<size_t>(dim)) {
    std::cerr << "PYVISTA: get_parallelepiped_edges: found " << gens.size()
              << " generators instead of " << dim << "\n";
    throw TerminalException{1};
  }
#endif
  // Identify each subset of generators with the matching row of M.
  std::vector<int> subset_to_row(n_row, -1);
  for (int subset = 0; subset < n_row; subset++) {
    MyVector<Tint> pt = p0;
    for (int i = 0; i < dim; i++) {
      if (subset & (1 << i)) {
        pt += gens[i];
      }
    }
    for (int i_row = 0; i_row < n_row; i_row++) {
      if (GetMatrixRow(M, i_row) == pt) {
        subset_to_row[subset] = i_row;
        break;
      }
    }
#ifdef SANITY_CHECK_ROBUST_PYVISTA
    if (subset_to_row[subset] == -1) {
      std::cerr << "PYVISTA: get_parallelepiped_edges: subset " << subset
                << " does not match a vertex\n";
      throw TerminalException{1};
    }
#endif
  }
  std::vector<std::pair<int, int>> edges;
  for (int subset = 0; subset < n_row; subset++) {
    for (int i = 0; i < dim; i++) {
      int subset_adj = subset | (1 << i);
      if (subset_adj != subset) {
        edges.push_back({subset_to_row[subset], subset_to_row[subset_adj]});
      }
    }
  }
  return edges;
}

// The cyclic ordering of the vertices of a 2-dimensional polytope
// (homogeneous coordinates of dimension 3). Each facet of the polygon
// is an edge incident to exactly 2 vertices, which gives the adjacency
// from which the cycle is read.
template <typename T>
std::vector<int> get_polygon_cycle(SinglePolytope<T> const &sp) {
  int n_ext = sp.EXT.rows();
#ifdef SANITY_CHECK_ROBUST_PYVISTA
  if (n_ext < 3) {
    std::cerr << "PYVISTA: get_polygon_cycle: only " << n_ext << " vertices\n";
    throw TerminalException{1};
  }
#endif
  std::vector<std::vector<int>> adj(n_ext);
  for (Face const &f : sp.facets) {
    std::vector<int> incd;
    for (int i_ext = 0; i_ext < n_ext; i_ext++) {
      if (f[i_ext] == 1) {
        incd.push_back(i_ext);
      }
    }
#ifdef SANITY_CHECK_ROBUST_PYVISTA
    if (incd.size() != 2) {
      std::cerr << "PYVISTA: get_polygon_cycle: a facet of the polygon is "
                << "incident to " << incd.size() << " vertices instead of 2\n";
      throw TerminalException{1};
    }
#endif
    adj[incd[0]].push_back(incd[1]);
    adj[incd[1]].push_back(incd[0]);
  }
#ifdef SANITY_CHECK_ROBUST_PYVISTA
  for (int i_ext = 0; i_ext < n_ext; i_ext++) {
    if (adj[i_ext].size() != 2) {
      std::cerr << "PYVISTA: get_polygon_cycle: vertex " << i_ext << " has "
                << adj[i_ext].size() << " neighbors instead of 2\n";
      throw TerminalException{1};
    }
  }
#endif
  std::vector<int> cycle;
  int prev = -1;
  int cur = 0;
  for (int step = 0; step < n_ext; step++) {
    cycle.push_back(cur);
    int next = adj[cur][0] != prev ? adj[cur][0] : adj[cur][1];
    prev = cur;
    cur = next;
  }
#ifdef SANITY_CHECK_ROBUST_PYVISTA
  if (cur != 0) {
    std::cerr << "PYVISTA: get_polygon_cycle: the walk does not close up\n";
    throw TerminalException{1};
  }
#endif
  return cycle;
}

// Accumulator for one JSON entry: the deduplicated vertices together
// with the points / lines / polygons expressed by vertex index. The
// vertices are kept in exact homogeneous coordinates (first coordinate
// normalized to 1) and converted to real coordinates only when writing.
template <typename T> struct PyVistaGeometry {
  std::unordered_map<MyVector<T>, int> map_vert;
  std::vector<MyVector<T>> l_vert;
  std::vector<int> point_indices;
  std::vector<std::pair<int, int>> line_cells;
  std::vector<int> line_colors; // one palette index per line cell
  std::vector<std::vector<int>> polygons;
  std::vector<int> polygon_groups;
  int get_point_index(MyVector<T> const &hom) {
#ifdef SANITY_CHECK_ROBUST_PYVISTA
    if (hom(0) == 0) {
      std::cerr << "PYVISTA: get_point_index: the point is at infinity\n";
      throw TerminalException{1};
    }
#endif
    MyVector<T> hom_norm = hom / hom(0);
    auto iter = map_vert.find(hom_norm);
    if (iter != map_vert.end()) {
      return iter->second;
    }
    int idx = l_vert.size();
    map_vert[hom_norm] = idx;
    l_vert.push_back(hom_norm);
    return idx;
  }
};

// The palette used to color the parallelotopes. The i-th parallelotope
// realizing the minimum gets the i-th color (cycling if there are more
// parallelotopes than colors).
inline std::vector<std::string> const &pyvista_palette() {
  static std::vector<std::string> const palette = {
      "black", "green",   "red",     "violet", "blue",
      "orange", "cyan", "magenta", "brown",  "gray"};
  return palette;
}

template <typename Tint, typename T>
MyVector<T> get_homogeneous_vector(MyVector<Tint> const &V) {
  int dim = V.size();
  MyVector<T> hom(dim + 1);
  hom(0) = T(1);
  for (int i = 0; i < dim; i++) {
    hom(i + 1) = UniversalScalarConversion<T, Tint>(V(i));
  }
  return hom;
}

// One JSON entry for one P-polytope.
template <typename T, typename Tint>
void write_pyvista_entry(std::ostream &os_out, PVoronoi<T, Tint> const &pv,
                         MyMatrix<double> const &L, std::ostream &os) {
  int dim = L.rows();
  PyVistaGeometry<T> geo;
  // The v_long vertex.
  MyVector<T> v_long_hom = get_homogeneous_vector<Tint, T>(pv.v_long);
  int idx_v_long = geo.get_point_index(v_long_hom);
  geo.point_indices.push_back(idx_v_long);
  // The parallelepipeds realizing the minimum, as wireframes. Each
  // parallelotope has its own color. An edge may belong to several
  // parallelotopes: we collect for each edge the set of parallelotopes
  // that contain it, then draw the edge as 3*k consecutive segments
  // cycling through the k colors (so C1-C2-...-Ck repeated three times).
  std::map<std::pair<int, int>, std::set<int>> map_edge;
  for (size_t p = 0; p < pv.l_robust_m_min.size(); p++) {
    MyMatrix<Tint> const &M = pv.l_robust_m_min[p].M;
    std::vector<std::pair<int, int>> edges = get_parallelepiped_edges(M);
    for (auto &edge : edges) {
      MyVector<Tint> V1 = GetMatrixRow(M, edge.first);
      MyVector<Tint> V2 = GetMatrixRow(M, edge.second);
      int idx1 = geo.get_point_index(get_homogeneous_vector<Tint, T>(V1));
      int idx2 = geo.get_point_index(get_homogeneous_vector<Tint, T>(V2));
      int a = std::min(idx1, idx2);
      int b = std::max(idx1, idx2);
      map_edge[{a, b}].insert(static_cast<int>(p));
    }
  }
  for (auto &kv : map_edge) {
    MyVector<T> A = geo.l_vert[kv.first.first];
    MyVector<T> B = geo.l_vert[kv.first.second];
    std::vector<int> ps(kv.second.begin(), kv.second.end());
    int k = ps.size();
    int n_seg = 3 * k;
    MyVector<T> diff = B - A;
    T den(n_seg);
    for (int j = 0; j < n_seg; j++) {
      T t0 = T(j) / den;
      T t1 = T(j + 1) / den;
      MyVector<T> P0 = A + t0 * diff;
      MyVector<T> P1 = A + t1 * diff;
      int idx0 = geo.get_point_index(P0);
      int idx1 = geo.get_point_index(P1);
      geo.line_cells.push_back({idx0, idx1});
      geo.line_colors.push_back(ps[j % k]);
    }
  }
  // The generalized polytope, as polygonal faces. All the polygons of
  // one planar facet share one polygon_group so that they get the same
  // color in the plot.
  BoundaryGeneralizedPolytope<T> bnd =
      find_generalized_polytope_boundary(pv.gp, os);
  int i_group = 0;
  for (auto &kv : bnd.full_data_facets) {
    MyMatrix<T> const &NSP = kv.second.NSP;
    auto f_insert_gp = [&](GeneralizedPolytope<T> const &gp_face) -> void {
      for (auto &sp : gp_face.polytopes) {
        std::vector<int> cycle = get_polygon_cycle(sp);
        std::vector<int> polygon;
        for (int i_ext : cycle) {
          MyVector<T> eEXT = GetMatrixRow(sp.EXT, i_ext);
          MyVector<T> hom = NSP.transpose() * eEXT;
          polygon.push_back(geo.get_point_index(hom));
        }
        geo.polygons.push_back(polygon);
        geo.polygon_groups.push_back(i_group);
      }
    };
    f_insert_gp(kv.second.gp_plus);
    f_insert_gp(kv.second.gp_minus);
    i_group += 1;
  }
#ifdef DEBUG_ROBUST_PYVISTA
  os << "PYVISTA: v_long=" << StringVectorGAP(pv.v_long)
     << " n_vert=" << geo.l_vert.size()
     << " n_lines=" << geo.line_cells.size()
     << " n_polygons=" << geo.polygons.size() << " n_groups=" << i_group
     << "\n";
#endif
  // The conversion to real coordinates: the lattice point x maps to
  // sum_i x_i v_i with v_i the rows of L.
  auto f_coord = [&](MyVector<T> const &hom, int j) -> double {
    double coord = 0;
    for (int i = 0; i < dim; i++) {
      double x_i = UniversalScalarConversion<double, T>(hom(i + 1));
      coord += x_i * L(i, j);
    }
    return coord;
  };
  os_out << "{\n";
  os_out << "  \"coordinates\": [\n";
  for (size_t i_vert = 0; i_vert < geo.l_vert.size(); i_vert++) {
    if (i_vert > 0) {
      os_out << ",\n";
    }
    os_out << "    [";
    for (int j = 0; j < dim; j++) {
      if (j > 0) {
        os_out << ", ";
      }
      os_out << f_coord(geo.l_vert[i_vert], j);
    }
    os_out << "]";
  }
  os_out << "\n  ],\n";
  os_out << "  \"point_indices\": [";
  for (size_t i = 0; i < geo.point_indices.size(); i++) {
    if (i > 0) {
      os_out << ", ";
    }
    os_out << geo.point_indices[i];
  }
  os_out << "],\n";
  os_out << "  \"lines\": [\n";
  for (size_t i_line = 0; i_line < geo.line_cells.size(); i_line++) {
    if (i_line > 0) {
      os_out << ",\n";
    }
    os_out << "    [" << geo.line_cells[i_line].first << ", "
           << geo.line_cells[i_line].second << "]";
  }
  os_out << "\n  ],\n";
  os_out << "  \"line_colors\": [";
  for (size_t i_line = 0; i_line < geo.line_colors.size(); i_line++) {
    if (i_line > 0) {
      os_out << ", ";
    }
    os_out << geo.line_colors[i_line];
  }
  os_out << "],\n";
  os_out << "  \"polygons\": [\n";
  for (size_t i_poly = 0; i_poly < geo.polygons.size(); i_poly++) {
    if (i_poly > 0) {
      os_out << ",\n";
    }
    os_out << "    [";
    for (size_t i = 0; i < geo.polygons[i_poly].size(); i++) {
      if (i > 0) {
        os_out << ", ";
      }
      os_out << geo.polygons[i_poly][i];
    }
    os_out << "]";
  }
  os_out << "\n  ],\n";
  os_out << "  \"polygon_groups\": [";
  for (size_t i = 0; i < geo.polygon_groups.size(); i++) {
    if (i > 0) {
      os_out << ", ";
    }
    os_out << geo.polygon_groups[i];
  }
  os_out << "],\n";
  os_out << "  \"style\": {\n";
  os_out << "    \"background\": \"white\",\n";
  os_out << "    \"window_size\": [1280, 720],\n";
  os_out << "    \"axes\": true,\n";
  os_out << "    \"points\": {\n";
  os_out << "      \"color\": \"black\",\n";
  os_out << "      \"size\": 36,\n";
  os_out << "      \"render_as_spheres\": true\n";
  os_out << "    },\n";
  os_out << "    \"lines\": {\n";
  os_out << "      \"color\": \"darkred\",\n";
  os_out << "      \"width\": 5,\n";
  os_out << "      \"render_as_tubes\": true,\n";
  os_out << "      \"palette\": [";
  {
    std::vector<std::string> const &palette = pyvista_palette();
    for (size_t i = 0; i < palette.size(); i++) {
      if (i > 0) {
        os_out << ", ";
      }
      os_out << "\"" << palette[i] << "\"";
    }
  }
  os_out << "]\n";
  os_out << "    },\n";
  os_out << "    \"polygons\": {\n";
  os_out << "      \"cmap\": \"tab20\",\n";
  os_out << "      \"opacity\": 0.72,\n";
  os_out << "      \"show_edges\": true,\n";
  os_out << "      \"edge_color\": \"black\",\n";
  os_out << "      \"edge_width\": 1.5,\n";
  os_out << "      \"lighting\": true,\n";
  os_out << "      \"show_scalar_bar\": false\n";
  os_out << "    }\n";
  os_out << "  },\n";
  os_out << "  \"camera\": {\n";
  os_out << "    \"parallel_projection\": true,\n";
  os_out << "    \"initial_direction\": [1.4, -1.8, 1.1],\n";
  os_out << "    \"view_up\": [0.0, 0.0, 1.0],\n";
  os_out << "    \"start_azimuth_degrees\": -55.0,\n";
  os_out << "    \"elevation_degrees\": 24.0,\n";
  os_out << "    \"elevation_swing_degrees\": 9.0\n";
  os_out << "  }\n";
  os_out << "}";
}

// The list of P-polytopes as a JSON list of entries, each in the
// format of geometry.json read by plot_geometry_index.py.
template <typename T, typename Tint>
void write_pyvista_json(std::ostream &os_out,
                        std::vector<PVoronoi<T, Tint>> const &l_ppoly,
                        MyMatrix<T> const &G, std::ostream &os) {
#ifdef SANITY_CHECK_ROBUST_PYVISTA
  if (G.rows() != 3) {
    std::cerr << "PYVISTA: write_pyvista_json is only for dimension 3\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<double> G_d = UniversalMatrixConversion<double, T>(G);
  MyMatrix<double> L = cholesky_lower(G_d);
  os_out << std::setprecision(17);
  os_out << "[\n";
  for (size_t i_ppoly = 0; i_ppoly < l_ppoly.size(); i_ppoly++) {
    if (i_ppoly > 0) {
      os_out << ",\n";
    }
    write_pyvista_entry(os_out, l_ppoly[i_ppoly], L, os);
  }
  os_out << "\n]\n";
}

// clang-format off
#endif  // SRC_ROBUST_COVERING_PYVISTA_JSON_H_
// clang-format on

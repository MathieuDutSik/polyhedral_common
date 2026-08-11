// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// Unit tests for the geometric bookkeeping of SoftBoundaryStore: the two
// primitives "add" (same-normal, keep the stored regions interior-disjoint)
// and "reduce" (opposite-normal, subtract an inserted block facet's footprint).
// No lattice / CVP involved: the fixtures are plain 2-D polytopes and the
// boundaries are the 1-D facets on the line x = 1.

// clang-format off
#include "NumberTheory.h"
#include "generalized_polytopes.h"
#include "robust_covering_types.h"
// clang-format on

using T = mpq_class;
using Tint = mpz_class;

static int n_fail = 0;

static void check(bool cond, std::string const &name) {
  if (cond) {
    std::cerr << "PASS: " << name << "\n";
  } else {
    std::cerr << "FAIL: " << name << "\n";
    n_fail += 1;
  }
}

static MyMatrix<T> mat_from_rows(std::vector<std::vector<int>> const &rows) {
  int n = rows.size();
  int m = rows[0].size();
  MyMatrix<T> M(n, m);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      M(i, j) = T(rows[i][j]);
    }
  }
  return M;
}

// Build the polytope from the homogeneous inequalities (a0 + a1 x + a2 y >= 0)
// and return its facet whose normal is "normal".
static ConvexBoundary<T> boundary_on(std::vector<std::vector<int>> const &fac_rows,
                                     std::vector<int> const &normal,
                                     std::ostream &os) {
  MyMatrix<T> FAC = mat_from_rows(fac_rows);
  SinglePolytope<T> sp = generate_single_polytope(FAC, os);
  MyVector<T> V(normal.size());
  for (size_t i = 0; i < normal.size(); i++) {
    V(i) = T(normal[i]);
  }
  MyVector<T> cn = ScalarCanonicalizationVector(V);
  int i_fac = get_matching_face_position(sp.FAC, cn);
  if (i_fac == -1) {
    std::cerr << "boundary_on: the requested normal is not a facet\n";
    throw TerminalException{1};
  }
  return get_convex_boundary(sp, i_fac, os);
}

static SoftConvexBoundary<T, Tint> scb_of(ConvexBoundary<T> const &cb) {
  return {cb, {}, {}};
}

// Every pair of regions stored at the same normal must be interior-disjoint.
static bool all_disjoint(SoftBoundaryStore<T, Tint> const &store, std::ostream &os) {
  for (auto &kv : store.map_scb) {
    std::vector<SoftConvexBoundary<T, Tint>> const &v = kv.second;
    for (size_t i = 0; i < v.size(); i++) {
      for (size_t j = i + 1; j < v.size(); j++) {
        if (is_pairwise_intersecting(v[i].cb.sp, v[j].cb.sp, os)) {
          return false;
        }
      }
    }
  }
  return true;
}

int main() {
  std::ostream &os = std::cerr;
  try {
    // Boundaries on the line x = 1, normal (1,-1,0) (interior x <= 1):
    ConvexBoundary<T> b_01 =
        boundary_on({{0, 1, 0}, {1, -1, 0}, {0, 0, 1}, {1, 0, -1}}, {1, -1, 0}, os);
    // footprint y in [1/2, 3/2], same normal:
    ConvexBoundary<T> b_half_32 =
        boundary_on({{0, 1, 0}, {1, -1, 0}, {-1, 0, 2}, {3, 0, -2}}, {1, -1, 0}, os);
    // footprint y in [1/4, 3/4] (subset of [0,1]), same normal:
    ConvexBoundary<T> b_sub =
        boundary_on({{0, 1, 0}, {1, -1, 0}, {-1, 0, 4}, {3, 0, -4}}, {1, -1, 0}, os);
    // From the right side, normal (-1,1,0) (interior x >= 1), footprint [0,1/2]:
    ConvexBoundary<T> b_right_half =
        boundary_on({{-1, 1, 0}, {2, -1, 0}, {0, 0, 1}, {1, 0, -2}}, {-1, 1, 0}, os);
    // From the right side, footprint [0,1] (covers all of b_01):
    ConvexBoundary<T> b_right_full =
        boundary_on({{-1, 1, 0}, {2, -1, 0}, {0, 0, 1}, {1, 0, -1}}, {-1, 1, 0}, os);
    // From the right side, footprint [1/4,3/4] (strictly inside b_01):
    ConvexBoundary<T> b_right_mid =
        boundary_on({{-1, 1, 0}, {2, -1, 0}, {-1, 0, 4}, {3, 0, -4}}, {-1, 1, 0}, os);

    // 3-D ambient: boundaries on the plane x = 1 are 2-D squares.
    // Box [0,1]^3, face x = 1 (normal (1,-1,0,0)), footprint (y,z) in [0,1]^2:
    ConvexBoundary<T> f_unit =
        boundary_on({{0, 1, 0, 0}, {1, -1, 0, 0}, {0, 0, 1, 0}, {1, 0, -1, 0},
                     {0, 0, 0, 1}, {1, 0, 0, -1}}, {1, -1, 0, 0}, os);
    // Box [0,1]x[1/2,3/2]x[0,1], face x = 1, footprint y in [1/2,3/2], z in [0,1]:
    ConvexBoundary<T> f_yhalf =
        boundary_on({{0, 1, 0, 0}, {1, -1, 0, 0}, {-1, 0, 2, 0}, {3, 0, -2, 0},
                     {0, 0, 0, 1}, {1, 0, 0, -1}}, {1, -1, 0, 0}, os);
    // Box [1,2]x[0,1]x[0,1], face x = 1 (normal (-1,1,0,0)), covers f_unit:
    ConvexBoundary<T> f_right_full =
        boundary_on({{-1, 1, 0, 0}, {2, -1, 0, 0}, {0, 0, 1, 0}, {1, 0, -1, 0},
                     {0, 0, 0, 1}, {1, 0, 0, -1}}, {-1, 1, 0, 0}, os);
    // Box [1,2]x[1/4,3/4]x[1/4,3/4], face x = 1, central square strictly inside f_unit:
    ConvexBoundary<T> f_right_mid =
        boundary_on({{-1, 1, 0, 0}, {2, -1, 0, 0}, {-1, 0, 4, 0}, {3, 0, -4, 0},
                     {-1, 0, 0, 4}, {3, 0, 0, -4}}, {-1, 1, 0, 0}, os);

    // 1. add fresh.
    {
      SoftBoundaryStore<T, Tint> store;
      store.add(scb_of(b_01), os);
      check(store.n_scb() == 1, "add_fresh_count");
      check(all_disjoint(store, os), "add_fresh_disjoint");
    }
    // 2. add same-key partial overlap => split, disjoint.
    {
      SoftBoundaryStore<T, Tint> store;
      store.add(scb_of(b_01), os);
      store.add(scb_of(b_half_32), os);
      check(store.n_scb() == 2, "add_partial_count");
      check(all_disjoint(store, os), "add_partial_disjoint");
    }
    // 3a. add a region fully contained in an existing one => no change.
    {
      SoftBoundaryStore<T, Tint> store;
      store.add(scb_of(b_01), os);
      store.add(scb_of(b_sub), os);
      check(store.n_scb() == 1, "add_contained_idempotent");
    }
    // 3b. add the identical region twice => no change.
    {
      SoftBoundaryStore<T, Tint> store;
      store.add(scb_of(b_01), os);
      store.add(scb_of(b_01), os);
      check(store.n_scb() == 1, "add_identical_idempotent");
    }
    // 4. reduce with an opposite-normal facet covering half => remnant stays,
    //    reducing again changes nothing.
    {
      SoftBoundaryStore<T, Tint> store;
      store.add(scb_of(b_01), os);
      store.reduce(b_right_half, os);
      check(store.n_scb() == 1, "reduce_partial_count");
      check(all_disjoint(store, os), "reduce_partial_disjoint");
      store.reduce(b_right_half, os);
      check(store.n_scb() == 1, "reduce_partial_idempotent");
    }
    // 5. reduce with an opposite-normal facet covering everything => empty.
    {
      SoftBoundaryStore<T, Tint> store;
      store.add(scb_of(b_01), os);
      store.reduce(b_right_full, os);
      check(store.n_scb() == 0, "reduce_full_empties");
      check(store.empty(), "reduce_full_empty_flag");
    }
    // 6. reduce with a normal that has no opposite entry => no-op.
    {
      SoftBoundaryStore<T, Tint> store;
      store.add(scb_of(b_01), os);
      store.reduce(b_01, os); // key -(1,-1,0) is absent
      check(store.n_scb() == 1, "reduce_no_match_noop");
    }
    // 7. reduce, middle cut (1-D): a facet strictly inside leaves two remnants.
    {
      SoftBoundaryStore<T, Tint> store;
      store.add(scb_of(b_01), os);
      store.reduce(b_right_mid, os);
      check(store.n_scb() == 2, "reduce_middle_two_remnants");
      check(all_disjoint(store, os), "reduce_middle_disjoint");
    }
    // 8. add, incoming strictly contains the stored one => two new remnants.
    {
      SoftBoundaryStore<T, Tint> store;
      store.add(scb_of(b_sub), os);  // [1/4,3/4]
      store.add(scb_of(b_01), os);   // [0,1] minus [1/4,3/4] = 2 pieces
      check(store.n_scb() == 3, "add_around_two_remnants");
      check(all_disjoint(store, os), "add_around_disjoint");
    }
    // 9. 3-D add fresh (2-D facet).
    {
      SoftBoundaryStore<T, Tint> store;
      store.add(scb_of(f_unit), os);
      check(store.n_scb() == 1, "add3d_fresh_count");
      check(all_disjoint(store, os), "add3d_fresh_disjoint");
    }
    // 10. 3-D add same-key partial overlap => split, disjoint.
    {
      SoftBoundaryStore<T, Tint> store;
      store.add(scb_of(f_unit), os);
      store.add(scb_of(f_yhalf), os);
      check(store.n_scb() == 2, "add3d_partial_count");
      check(all_disjoint(store, os), "add3d_partial_disjoint");
    }
    // 11. 3-D reduce full cover => empty.
    {
      SoftBoundaryStore<T, Tint> store;
      store.add(scb_of(f_unit), os);
      store.reduce(f_right_full, os);
      check(store.n_scb() == 0, "reduce3d_full_empties");
    }
    // 12. 3-D reduce, central hole => the face minus a middle square is a frame
    //     (difference_p_p returns it as 4 rectangles).
    {
      SoftBoundaryStore<T, Tint> store;
      store.add(scb_of(f_unit), os);
      store.reduce(f_right_mid, os);
      check(store.n_scb() == 4, "reduce3d_frame_count");
      check(all_disjoint(store, os), "reduce3d_frame_disjoint");
    }
    // 13. mixed 3-D sequence: the disjoint invariant must hold after every step.
    {
      SoftBoundaryStore<T, Tint> store;
      bool ok = true;
      store.add(scb_of(f_unit), os);   ok = ok && all_disjoint(store, os);
      store.add(scb_of(f_yhalf), os);  ok = ok && all_disjoint(store, os);
      store.reduce(f_right_mid, os);   ok = ok && all_disjoint(store, os);
      store.add(scb_of(f_unit), os);   ok = ok && all_disjoint(store, os);
      store.reduce(f_right_full, os);  ok = ok && all_disjoint(store, os);
      check(ok, "sequence3d_disjoint_invariant");
    }
  } catch (TerminalException const &e) {
    std::cerr << "Test aborted with TerminalException\n";
    return e.eVal;
  }
  if (n_fail == 0) {
    std::cerr << "All SoftBoundaryStore tests passed\n";
    return 0;
  }
  std::cerr << n_fail << " SoftBoundaryStore test(s) failed\n";
  return 1;
}

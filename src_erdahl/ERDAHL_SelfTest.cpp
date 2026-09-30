// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "erdahl_enumeration.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

/*
  Checks of the src_erdahl functionality on cases where the answer is known
  by hand. Returns a nonzero exit code if any check fails.
 */

struct TestCounter {
  int n_pass = 0;
  int n_fail = 0;
  void check(bool test, std::string const &name) {
    if (test) {
      n_pass++;
      std::cerr << "PASS: " << name << "\n";
    } else {
      n_fail++;
      std::cerr << "FAIL: " << name << "\n";
    }
  }
};

template <typename T>
MyMatrix<T> parse_matrix(int n_row, int n_col, std::vector<T> const &l_val) {
  MyMatrix<T> M(n_row, n_col);
  int pos = 0;
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n_col; j++) {
      M(i, j) = l_val[pos];
      pos++;
    }
  }
  return M;
}

// The function sum_i x_i (x_i - 1) whose zero set is the cube {0,1}^n.
template <typename T> MyMatrix<T> cube_function(int n) {
  MyMatrix<T> F = ZeroMatrix<T>(n + 1, n + 1);
  for (int i = 0; i < n; i++) {
    F(i + 1, i + 1) = 1;
    F(0, i + 1) = T(-1) / T(2);
    F(i + 1, 0) = T(-1) / T(2);
  }
  return F;
}

template <typename Tint>
MyMatrix<Tint> lattice_basis(int n, std::vector<std::vector<int>> const &l_v) {
  MyMatrix<Tint> L(l_v.size(), n);
  for (size_t i = 0; i < l_v.size(); i++) {
    for (int j = 0; j < n; j++) {
      L(i, j) = l_v[i][j];
    }
  }
  return L;
}

template <typename T, typename Tint, typename Tgroup>
void test_zero_sets(TestCounter &tc, std::ostream &os) {
  // Indefinite quadratic form: x1^2 - x2^2.
  {
    MyMatrix<T> F = parse_matrix<T>(3, 3, {0, 0, 0, 0, 1, 0, 0, 0, -1});
    auto res = erdahl_zero_set<T, Tint>(F, os);
    bool test = res.is_err();
    if (test) {
      for (auto &e : res.get_err().points) {
        if (EvaluationQuadForm<T, Tint>(F, e) >= 0) {
          test = false;
        }
      }
    }
    tc.check(test, "zero set: indefinite form gives negative points");
  }
  // Linear part not vanishing on the kernel: x1^2 + x2.
  {
    MyMatrix<T> F = parse_matrix<T>(
        3, 3, {0, 0, T(1) / T(2), 0, 1, 0, T(1) / T(2), 0, 0});
    auto res = erdahl_zero_set<T, Tint>(F, os);
    bool test = res.is_err() && res.get_err().points.size() > 0 &&
                EvaluationQuadForm<T, Tint>(F, res.get_err().points[0]) < 0;
    tc.check(test, "zero set: linear part on the kernel gives negative points");
  }
  // Negative minimum: (x1 - 1/2)^2 + x2^2 - 1 = x1^2 - x1 + x2^2 - 3/4.
  {
    MyMatrix<T> F = parse_matrix<T>(
        3, 3, {T(-3) / T(4), T(-1) / T(2), 0, T(-1) / T(2), 1, 0, 0, 0, 1});
    auto res = erdahl_zero_set<T, Tint>(F, os);
    tc.check(res.is_err(), "zero set: negative minimum");
  }
  // Positive minimum: x1^2 + x2^2 + 1, empty zero set.
  {
    MyMatrix<T> F = parse_matrix<T>(3, 3, {1, 0, 0, 0, 1, 0, 0, 0, 1});
    auto res = erdahl_zero_set<T, Tint>(F, os);
    tc.check(res.is_ok() && res.get_ok().EXT.rows() == 0,
             "zero set: positive minimum gives the empty set");
  }
  // The zero function: the zero set is Z^n.
  {
    MyMatrix<T> F = ZeroMatrix<T>(4, 4);
    DelaunayPolyhedron<T, Tint> D =
        erdahl_polyhedron_from_function<T, Tint>(F, os);
    tc.check(D.EXT.rows() == 1 && D.L.rows() == 3,
             "zero set: the zero function gives Z^3");
  }
  // The strip {0,1} x Z.
  {
    MyMatrix<T> F = parse_matrix<T>(
        3, 3, {0, T(-1) / T(2), 0, T(-1) / T(2), 1, 0, 0, 0, 0});
    DelaunayPolyhedron<T, Tint> D =
        erdahl_polyhedron_from_function<T, Tint>(F, os);
    tc.check(D.EXT.rows() == 2 && D.L.rows() == 1,
             "zero set: the strip {0,1} x Z");
  }
}

template <typename T, typename Tint>
void test_spaces(TestCounter &tc, std::ostream &) {
  int n = 3;
  MyVector<T> c(n);
  c(0) = T(1) / T(2);
  c(1) = 0;
  c(2) = T(1) / T(2);
  ErdahlFunctionSpace<T> Wc = erdahl_centered_space(c);
  // The point reflection x -> 2c - x as a row transformation.
  MyMatrix<Tint> g = ZeroMatrix<Tint>(n + 1, n + 1);
  g(0, 0) = 1;
  for (int i = 0; i < n; i++) {
    g(0, i + 1) = UniversalScalarConversion<Tint, T>(2 * c(i));
    g(i + 1, i + 1) = -1;
  }
  ErdahlFunctionSpace<T> Winv = erdahl_invariant_space<T, Tint>(n, {g});
  bool test = Wc.basis.size() == Winv.basis.size();
  for (auto &F : Winv.basis) {
    if (!erdahl_is_in_space(Wc, F)) {
      test = false;
    }
  }
  tc.check(test, "spaces: the invariant space of the reflection is the "
                 "centered space");
  tc.check(erdahl_preserves_space(Wc, g),
           "spaces: the reflection preserves the centered space");
  ErdahlFunctionSpace<T> Wfull = erdahl_full_space<T>(n);
  tc.check(Wfull.basis.size() == 10, "spaces: dimension of the full space");
}

template <typename T, typename Tint, typename Tgroup>
size_t count_facets_no_group(MyMatrix<T> const &EXTcone, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tgroup grp_triv(std::vector<Telt>{}, Telt(static_cast<Tidx>(EXTcone.rows())));
  return DualDescriptionStandard<T, Tgroup>(EXTcone, grp_triv, os).size();
}

template <typename T, typename Tint, typename Tgroup>
void test_polytopes(TestCounter &tc, std::ostream &os) {
  // The square {0,1}^2 in the full space: rank 2, the facets are the four
  // triangles, in a single orbit.
  {
    int n = 2;
    ErdahlFunctionSpace<T> W = erdahl_full_space<T>(n);
    DelaunayPolyhedron<T, Tint> D =
        erdahl_polyhedron_from_function<T, Tint>(cube_function<T>(n), os);
    tc.check(D.EXT.rows() == 4 && D.L.rows() == 0, "square: vertices");
    tc.check(erdahl_perfection_rank(W, D) == 2, "square: rank 2");
    ErdahlPolytopeGroup<T, Tint, Tgroup> grp =
        erdahl_polytope_group<T, Tint, Tgroup>(W, D, os);
    tc.check(grp.grp.size() == 8, "square: group of order 8");
    ErdahlSubPolyhedraPolytope<Tint> sub =
        erdahl_sub_polyhedra_polytope<T, Tint, Tgroup>(W, D, "unset", os);
    bool test = sub.l_face.size() == 1 && sub.l_face_lowdim.size() == 0 &&
                sub.l_face[0].count() == 3;
    tc.check(test, "square: one orbit of triangles");
    for (auto &eFace : sub.l_face) {
      DelaunayPolyhedron<T, Tint> Dsub =
          erdahl_sub_polyhedron_of_face(W, D, eFace, os);
      tc.check(erdahl_perfection_rank(W, Dsub) == 3 &&
                   erdahl_is_subset(Dsub, D),
               "square: the triangle has rank 3");
    }
  }
  // The cube {0,1}^3 in the full space: since x_i^2 = x_i on the cube the
  // evaluation vectors span a space of dimension 7, so the rank is 3. By
  // Theorem "CUTn_and_ErdahlHn" of the paper the cone of evaluations is the
  // cone over the cut polytope CUT_4, whose 16 facets are the 12 triangle
  // inequalities and the 4 perimeter inequalities.
  {
    int n = 3;
    ErdahlFunctionSpace<T> W = erdahl_full_space<T>(n);
    DelaunayPolyhedron<T, Tint> D =
        erdahl_polyhedron_from_function<T, Tint>(cube_function<T>(n), os);
    tc.check(erdahl_perfection_rank(W, D) == 3, "cube: rank 3");
    ErdahlEvaluationClasses<T> ec = erdahl_evaluation_classes(W, D.EXT);
    size_t n_fac = count_facets_no_group<T, Tint, Tgroup>(ec.EXTcone, os);
    tc.check(n_fac == 16, "cube: 16 facets of the evaluation cone");
    ErdahlSubPolyhedraPolytope<Tint> sub =
        erdahl_sub_polyhedra_polytope<T, Tint, Tgroup>(W, D, "unset", os);
    ErdahlPolytopeGroup<T, Tint, Tgroup> grp =
        erdahl_polytope_group<T, Tint, Tgroup>(W, D, os);
    tc.check(grp.grp.size() == 48, "cube: group of order 48");
    // Sum of the orbit sizes |G| / |Stab(face)| should be 16.
    using Tint_grp = typename Tgroup::Tint;
    Tint_grp sum(0);
    for (auto &eFace : sub.l_face) {
      Tgroup stab = grp.grp.Stabilizer_OnSets(eFace);
      sum += grp.grp.size() / stab.size();
      DelaunayPolyhedron<T, Tint> Dsub =
          erdahl_sub_polyhedron_of_face(W, D, eFace, os);
      tc.check(erdahl_perfection_rank(W, Dsub) == 4,
               "cube: the sub-polyhedron has rank 4");
    }
    tc.check(sum == 16 && sub.l_face_lowdim.empty(),
             "cube: the orbits cover the 16 facets");
  }
  // The cube {0,1}^3 in the space of functions of center (1/2,1/2,1/2):
  // rank 3, the sub-polyhedra are the cube minus an antipodal pair, one
  // orbit of 4.
  {
    int n = 3;
    MyVector<T> c(n);
    for (int i = 0; i < n; i++) {
      c(i) = T(1) / T(2);
    }
    ErdahlFunctionSpace<T> W = erdahl_centered_space(c);
    DelaunayPolyhedron<T, Tint> D =
        erdahl_polyhedron_from_function<T, Tint>(cube_function<T>(n), os);
    tc.check(erdahl_is_in_space(W, D.F), "centered cube: F is centered");
    tc.check(erdahl_perfection_rank(W, D) == 3, "centered cube: rank 3");
    ErdahlSubPolyhedraPolytope<Tint> sub =
        erdahl_sub_polyhedra_polytope<T, Tint, Tgroup>(W, D, "unset", os);
    bool test = sub.l_face.size() == 1 && sub.l_face[0].count() == 6;
    tc.check(test, "centered cube: one orbit of 6-vertex sub-polytopes");
    for (auto &eFace : sub.l_face) {
      DelaunayPolyhedron<T, Tint> Dsub =
          erdahl_sub_polyhedron_of_face(W, D, eFace, os);
      tc.check(erdahl_perfection_rank(W, Dsub) == 4 &&
                   erdahl_is_in_space(W, Dsub.F),
               "centered cube: the sub-polytope has rank 4 and is centered");
    }
  }
  // Equivalence of the square with its image under an affine map.
  {
    int n = 2;
    ErdahlFunctionSpace<T> W = erdahl_full_space<T>(n);
    DelaunayPolyhedron<T, Tint> D1 =
        erdahl_polyhedron_from_function<T, Tint>(cube_function<T>(n), os);
    MyMatrix<Tint> g =
        parse_matrix<Tint>(3, 3, {1, 3, -2, 0, 1, 1, 0, 1, 2});
    MyMatrix<Tint> ginv = Inverse(g);
    // The image of D1 under g is the zero set of F o g^{-1}.
    MyMatrix<T> F2 = erdahl_transform_function(D1.F, ginv);
    DelaunayPolyhedron<T, Tint> D2 =
        erdahl_polyhedron_from_function<T, Tint>(F2, os);
    std::optional<MyMatrix<Tint>> opt =
        erdahl_polytope_equivalence<T, Tint, Tgroup>(W, D1, D2, os);
    bool test = opt.has_value();
    if (opt) {
      MyMatrix<Tint> EXTimg = D1.EXT * (*opt);
      std::vector<MyVector<Tint>> l_row;
      for (int i = 0; i < EXTimg.rows(); i++) {
        l_row.push_back(GetMatrixRow(EXTimg, i));
      }
      test = erdahl_sorted_rows(l_row, n + 1) == D2.EXT;
    }
    tc.check(test, "equivalence: the square and its image");
  }
}

template <typename T, typename Tint, typename Tgroup>
void test_flips(TestCounter &tc, std::ostream &os) {
  int n = 2;
  // Full space: square subset {0,1} x Z subset Z^2, the flip is Z x {0,1}.
  {
    ErdahlFunctionSpace<T> W = erdahl_full_space<T>(n);
    DelaunayPolyhedron<T, Tint> D3 =
        erdahl_polyhedron_from_function<T, Tint>(ZeroMatrix<T>(3, 3), os);
    MyMatrix<T> F2 = parse_matrix<T>(
        3, 3, {0, T(-1) / T(2), 0, T(-1) / T(2), 1, 0, 0, 0, 0});
    DelaunayPolyhedron<T, Tint> D2 =
        erdahl_polyhedron_from_function<T, Tint>(F2, os);
    DelaunayPolyhedron<T, Tint> D1 =
        erdahl_polyhedron_from_function<T, Tint>(cube_function<T>(n), os);
    tc.check(erdahl_perfection_rank(W, D3) == 0 &&
                 erdahl_perfection_rank(W, D2) == 1 &&
                 erdahl_perfection_rank(W, D1) == 2,
             "flip full: ranks 2, 1, 0");
    DelaunayPolyhedron<T, Tint> D2p = erdahl_flip(W, D1, D2, D3, true, os);
    MyMatrix<T> F2exp = parse_matrix<T>(
        3, 3, {0, 0, T(-1) / T(2), 0, 0, 0, T(-1) / T(2), 0, 1});
    DelaunayPolyhedron<T, Tint> D2exp =
        erdahl_polyhedron_from_function<T, Tint>(F2exp, os);
    tc.check(erdahl_is_equal(D2p, D2exp), "flip full: the flip is Z x {0,1}");
  }
  // Centered space of center (1/2, 0): parallelogram subset {0,1} x Z subset
  // Z^2, the flip is {x : x1 + x2 in {0,1}}.
  {
    MyVector<T> c(n);
    c(0) = T(1) / T(2);
    c(1) = 0;
    ErdahlFunctionSpace<T> W = erdahl_centered_space(c);
    DelaunayPolyhedron<T, Tint> D3 =
        erdahl_polyhedron_from_function<T, Tint>(ZeroMatrix<T>(3, 3), os);
    MyMatrix<T> F2 = parse_matrix<T>(
        3, 3, {0, T(-1) / T(2), 0, T(-1) / T(2), 1, 0, 0, 0, 0});
    DelaunayPolyhedron<T, Tint> D2 =
        erdahl_polyhedron_from_function<T, Tint>(F2, os);
    MyMatrix<T> F1 =
        parse_matrix<T>(3, 3, {0, -1, T(-1) / T(2), -1, 2, 1, T(-1) / T(2), 1, 1});
    DelaunayPolyhedron<T, Tint> D1 =
        erdahl_polyhedron_from_function<T, Tint>(F1, os);
    tc.check(D1.EXT.rows() == 4 && D1.L.rows() == 0,
             "flip centered: the parallelogram");
    tc.check(erdahl_is_in_space(W, F1) && erdahl_is_in_space(W, F2),
             "flip centered: the functions are centered");
    tc.check(erdahl_perfection_rank(W, D3) == 0 &&
                 erdahl_perfection_rank(W, D2) == 1 &&
                 erdahl_perfection_rank(W, D1) == 2,
             "flip centered: ranks 2, 1, 0");
    DelaunayPolyhedron<T, Tint> D2p = erdahl_flip(W, D1, D2, D3, true, os);
    tc.check(erdahl_is_in_space(W, D2p.F), "flip centered: F is centered");
    // The expected polyhedron: (x1 + x2)(x1 + x2 - 1).
    MyMatrix<T> F2exp = parse_matrix<T>(
        3, 3,
        {0, T(-1) / T(2), T(-1) / T(2), T(-1) / T(2), 1, 1, T(-1) / T(2), 1, 1});
    DelaunayPolyhedron<T, Tint> D2exp =
        erdahl_polyhedron_from_function<T, Tint>(F2exp, os);
    tc.check(erdahl_is_equal(D2p, D2exp),
             "flip centered: the flip is {x1 + x2 in {0,1}}");
  }
}

template <typename T, typename Tint, typename Tgroup>
void test_enumeration(TestCounter &tc, std::ostream &os) {
  // In the full space, for n <= 5 the only perfect Delaunay polyhedron is
  // {0,1} x Z^{n-1} (Erdahl 1992, Deza-Grishukhin-Laurent 1992).
  for (int n = 1; n <= 4; n++) {
    ErdahlFunctionSpace<T> W = erdahl_full_space<T>(n);
    std::vector<DelaunayPolyhedron<T, Tint>> l_perf =
        erdahl_enumerate_perfect<T, Tint, Tgroup>(W, "unset", os);
    bool test = l_perf.size() == 1;
    for (auto &D : l_perf) {
      if (erdahl_perfection_rank(W, D) != 1) {
        test = false;
      }
    }
    if (l_perf.size() == 1) {
      test = test && l_perf[0].EXT.rows() == 2 &&
             l_perf[0].L.rows() == n - 1;
    }
    tc.check(test, "enumeration full n=" + std::to_string(n) +
                       ": the only perfect one is {0,1} x Z^{n-1}");
  }
  // The centered spaces of center (1/2, 0, ..., 0). A perfect polytope
  // for them needs at least n(n+1)/2 antipodal pairs of vertices, which
  // does not happen for small n, so the only perfect one is a {0,1} strip.
  for (int n = 1; n <= 4; n++) {
    MyVector<T> c = ZeroVector<T>(n);
    c(0) = T(1) / T(2);
    ErdahlFunctionSpace<T> W = erdahl_centered_space(c);
    std::vector<DelaunayPolyhedron<T, Tint>> l_perf =
        erdahl_enumerate_perfect<T, Tint, Tgroup>(W, "unset", os);
    bool test = l_perf.size() == 1;
    for (auto &D : l_perf) {
      if (erdahl_perfection_rank(W, D) != 1 || !erdahl_is_in_space(W, D.F) ||
          D.EXT.rows() != 2 || D.L.rows() != n - 1) {
        test = false;
      }
    }
    tc.check(test, "enumeration centered n=" + std::to_string(n) +
                       ": the only perfect one is a {0,1} strip");
  }
}

int main() {
  using T = mpq_class;
  using Tint = mpz_class;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  HumanTime time;
  TestCounter tc;
  try {
    std::ostream &os = std::cerr;
    test_zero_sets<T, Tint, Tgroup>(tc, os);
    test_spaces<T, Tint>(tc, os);
    test_polytopes<T, Tint, Tgroup>(tc, os);
    test_flips<T, Tint, Tgroup>(tc, os);
    test_enumeration<T, Tint, Tgroup>(tc, os);
  } catch (TerminalException const &e) {
    std::cerr << "ERDAHL_SelfTest: an exception was raised\n";
    exit(e.eVal);
  }
  std::cerr << "ERDAHL_SelfTest: n_pass=" << tc.n_pass
            << " n_fail=" << tc.n_fail << " time=" << time << "\n";
  if (tc.n_fail > 0) {
    return 1;
  }
  return 0;
}

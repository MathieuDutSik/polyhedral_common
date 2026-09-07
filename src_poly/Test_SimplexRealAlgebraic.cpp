// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// Linear programming over a real algebraic field.
//
// SIMPLEX_LinearProgramming reduces to SIMPLEX_LinearProgramming_ring, which
// runs the whole enumeration over underlying_ring<T>. For a real algebraic
// field that ring is the order Z[x] spanned by the powers of the generator, a
// type that is neither a field nor a euclidean domain, and nothing else
// exercised that combination on a linear program. The rational cases of
// CI(15) cannot: their inequalities are read from GAP records, which carry no
// algebraic entries, and the GAP output format does not accept them either.
// So the program is built here rather than driven from a data file.
//
// The field is the cubic field of discriminant 49 that basic_common_cpp ships
// for its own tests: the generator is 2*cos(2*pi/7), of minimal polynomial
// X^3 + X^2 - 2X - 1, monic so that the powers of the generator span a ring.
//
// The rows of ListIneq are f_i(x) = b_i + a_i . x >= 0 and the objective is
// ToBeMinimized, in the same convention.

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "POLY_SimplexClarkson.h"
// clang-format on
#include <string>
#include <vector>

int const idx_field = 1;
using T = RealField<idx_field>;

static int n_error = 0;

static void check(bool test, std::string const &name) {
  if (test) {
    std::cerr << "PASS: " << name << "\n";
  } else {
    std::cerr << "FAIL: " << name << "\n";
    n_error++;
  }
}

// The element c0 + c1 x + c2 x^2 of the field.
static T Elt(int c0, int c1, int c2) {
  std::vector<mpq_class> V{mpq_class(c0), mpq_class(c1), mpq_class(c2)};
  return T(V);
}

// f_i(x) for the row i of the system, at the point x.
static T EvalRow(MyMatrix<T> const &M, int iRow, MyVector<T> const &x) {
  T val = M(iRow, 0);
  for (int j = 1; j < M.cols(); j++)
    val += M(iRow, j) * x(j - 1);
  return val;
}

// The box 0 <= x_i <= c_i with algebraic bounds, minimizing -sum x_i. The
// optimum is at x_i = c_i and is worth -sum c_i, which is known exactly and
// is a genuine element of the field, not a rational.
static void process_box(int n, std::vector<T> const &bounds) {
  MyMatrix<T> ListIneq(2 * n, n + 1);
  for (int i = 0; i < 2 * n; i++)
    for (int j = 0; j <= n; j++)
      ListIneq(i, j) = T(0);
  for (int i = 0; i < n; i++) {
    ListIneq(2 * i, i + 1) = T(1);       // x_i >= 0
    ListIneq(2 * i + 1, 0) = bounds[i];  // c_i - x_i >= 0
    ListIneq(2 * i + 1, i + 1) = T(-1);
  }
  MyVector<T> obj(n + 1);
  obj(0) = T(0);
  for (int i = 1; i <= n; i++)
    obj(i) = T(-1);
  T expected(0);
  for (int i = 0; i < n; i++)
    expected -= bounds[i];
  //
  LpSolution<T> sol = SIMPLEX_LinearProgramming(ListIneq, obj, std::cerr);
  if (!sol.DirectSolution) {
    check(false, "the bounded program has a primal solution");
    return;
  }
  check(true, "the bounded program has a primal solution");
  MyVector<T> const &x = *sol.DirectSolution;
  bool feasible = true;
  for (int i = 0; i < 2 * n; i++)
    if (EvalRow(ListIneq, i, x) < 0)
      feasible = false;
  check(feasible, "the primal solution is feasible");
  T value = obj(0);
  for (int i = 1; i <= n; i++)
    value += obj(i) * x(i - 1);
  check(value == sol.OptimalValue,
        "the objective at the primal solution is OptimalValue");
  check(sol.OptimalValue == expected,
        "OptimalValue is the exact algebraic optimum");
  std::cerr << "INFO: OptimalValue=" << sol.OptimalValue
            << " expected=" << expected << "\n";
}

// An infeasible system: x >= c and x <= c - 1 with c algebraic. Only the dual
// is set, and it is a Farkas certificate.
static void process_infeasible() {
  T c = Elt(0, 1, 0);
  MyMatrix<T> ListIneq(2, 2);
  ListIneq(0, 0) = -c;
  ListIneq(0, 1) = T(1);  // x - c >= 0
  ListIneq(1, 0) = c - T(1);
  ListIneq(1, 1) = T(-1); // c - 1 - x >= 0
  MyVector<T> obj(2);
  obj(0) = T(0);
  obj(1) = T(1);
  LpSolution<T> sol = SIMPLEX_LinearProgramming(ListIneq, obj, std::cerr);
  check(!sol.DirectSolution, "the infeasible program has no primal solution");
  check(bool(sol.DualSolution), "the infeasible program returns a dual certificate");
}

int main() {
  HumanTime time;
  try {
    std::string eFile =
        "basic_common_cpp/CI_tests/RealAlgebraicField/CubicFieldDisc_49";
    bool found = false;
    for (int level = 0; level <= 10; level++) {
      if (FILE_IsExistingFile(eFile)) {
        found = true;
        break;
      }
      eFile = "../" + eFile;
    }
    if (!found) {
      std::cerr << "Failed to find the cubic field description after checking "
                   "paths from basic_common_cpp/CI_tests/ up to 10 parent "
                   "levels\n";
      throw TerminalException{1};
    }
    HelperClassRealField<mpq_class> hcrf(eFile);
    insert_helper_real_algebraic_field(idx_field, hcrf);
    check(hcrf.is_monic(), "the minimal polynomial is monic");
    std::cerr << "STEP 1: the field is registered, degree " << hcrf.deg << "\n";
    //
    // x = 2*cos(2*pi/7) is about 1.247, so the bounds below are positive and
    // the box is not empty.
    process_box(3, {Elt(0, 1, 0), Elt(1, 1, 0), Elt(-1, 0, 1)});
    std::cerr << "STEP 2: the bounded program\n";
    process_box(4, {Elt(2, 0, 0), Elt(0, 1, 0), Elt(1, 0, 1), Elt(3, -1, 1)});
    std::cerr << "STEP 3: a second bounded program\n";
    process_infeasible();
    std::cerr << "STEP 4: the infeasible program\n";
    //
    std::cerr << "n_error=" << n_error << "\n";
    if (n_error > 0) {
      std::cerr << "Erroneous termination of Test_SimplexRealAlgebraic\n";
      return 1;
    }
    std::cerr << "Normal termination of Test_SimplexRealAlgebraic\n";
    std::cerr << "runtime = " << time << "\n";
    return 0;
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happened\n";
    exit(e.eVal);
  }
}

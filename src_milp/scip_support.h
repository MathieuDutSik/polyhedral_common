// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MILP_SCIP_SUPPORT_H_
#define SRC_MILP_SCIP_SUPPORT_H_

/*
  The SCIP specific code. It is included only when compiled with
  ENABLE_SCIP_SUPPORT:

  #ifdef ENABLE_SCIP_SUPPORT
  #include "scip_support.h"
  #endif

  SCIP works in floating point. The rows are scaled to integers, which must
  be exactly representable as doubles, and the returned solution is rounded
  and checked exactly against the constraints.
 */

// clang-format off
#include "ilp_fundamental.h"
#include "MAT_MatrixInt.h"
#include <cmath>
#include <string>
#include <vector>
// The SCIP headers are third party, their warnings are not ours to fix.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#pragma GCC diagnostic pop
// clang-format on

inline void ILP_ScipCheck(SCIP_RETCODE rc, char const *where) {
  if (rc != SCIP_OKAY) {
    std::cerr << "ILP: SCIP call " << where << " failed with code "
              << static_cast<int>(rc) << "\n";
    throw TerminalException{1};
  }
}

// Scale a row by a positive integer so that it becomes integral, and
// convert it to doubles, checking that the conversion is exact.
template <typename T>
std::vector<double> ILP_ScaledRowDouble(MyVector<T> const &row) {
  std::vector<double> ret(row.size(), 0.0);
  if (IsZeroVector(row)) {
    return ret;
  }
  // The scaling factor of RemoveFractionVector is positive.
  MyVector<T> row_int = RemoveFractionVector(row);
  double const max_exact = 9007199254740992.0; // 2^53
  for (int j = 0; j < row.size(); j++) {
    double val = UniversalScalarConversion<double, T>(row_int(j));
    if (std::abs(val) >= max_exact) {
      std::cerr << "ILP: A coefficient is too large to be passed exactly to "
                   "SCIP\n";
      throw TerminalException{1};
    }
    ret[j] = val;
  }
  return ret;
}

template <typename T>
IlpSolution<T> ILP_Scip(MyMatrix<T> const &ListIneq,
                        MyMatrix<T> const &ListEqua,
                        MyVector<T> const &ToBeMinimized,
                        [[maybe_unused]] std::ostream &os) {
  int n_col = ToBeMinimized.size();
  int d = n_col - 1;
  SCIP *scip = nullptr;
  ILP_ScipCheck(SCIPcreate(&scip), "SCIPcreate");
  ILP_ScipCheck(SCIPincludeDefaultPlugins(scip), "SCIPincludeDefaultPlugins");
  SCIPsetMessagehdlrQuiet(scip, TRUE);
  ILP_ScipCheck(SCIPcreateProbBasic(scip, "ilp"), "SCIPcreateProbBasic");
  ILP_ScipCheck(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE),
                "SCIPsetObjsense");
  std::vector<double> obj = ILP_ScaledRowDouble(ToBeMinimized);
  std::vector<SCIP_VAR *> vars(d, nullptr);
  double infinity = SCIPinfinity(scip);
  for (int j = 0; j < d; j++) {
    std::string name = "x" + std::to_string(j);
    ILP_ScipCheck(SCIPcreateVarBasic(scip, &vars[j], name.c_str(), -infinity,
                                     infinity, obj[j + 1],
                                     SCIP_VARTYPE_INTEGER),
                  "SCIPcreateVarBasic");
    ILP_ScipCheck(SCIPaddVar(scip, vars[j]), "SCIPaddVar");
  }
  auto add_row = [&](MyVector<T> const &row, bool is_equality,
                     std::string const &name) -> void {
    std::vector<double> row_d = ILP_ScaledRowDouble(row);
    std::vector<SCIP_VAR *> l_var;
    std::vector<double> l_val;
    for (int j = 0; j < d; j++) {
      if (row_d[j + 1] != 0) {
        l_var.push_back(vars[j]);
        l_val.push_back(row_d[j + 1]);
      }
    }
    double lhs = -row_d[0];
    double rhs = is_equality ? lhs : infinity;
    SCIP_CONS *cons = nullptr;
    ILP_ScipCheck(SCIPcreateConsBasicLinear(scip, &cons, name.c_str(),
                                            l_var.size(), l_var.data(),
                                            l_val.data(), lhs, rhs),
                  "SCIPcreateConsBasicLinear");
    ILP_ScipCheck(SCIPaddCons(scip, cons), "SCIPaddCons");
    ILP_ScipCheck(SCIPreleaseCons(scip, &cons), "SCIPreleaseCons");
  };
  for (int i = 0; i < ListIneq.rows(); i++) {
    add_row(GetMatrixRow(ListIneq, i), false, "ineq" + std::to_string(i));
  }
  for (int i = 0; i < ListEqua.rows(); i++) {
    add_row(GetMatrixRow(ListEqua, i), true, "equa" + std::to_string(i));
  }
  ILP_ScipCheck(SCIPsolve(scip), "SCIPsolve");
  SCIP_STATUS status = SCIPgetStatus(scip);
  if (status == SCIP_STATUS_INFORUNBD) {
    // The presolving could not decide, redo without it.
    ILP_ScipCheck(SCIPfreeTransform(scip), "SCIPfreeTransform");
    ILP_ScipCheck(SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE),
                  "SCIPsetPresolving");
    ILP_ScipCheck(SCIPsolve(scip), "SCIPsolve");
    status = SCIPgetStatus(scip);
  }
  auto clean = [&]() -> void {
    for (int j = 0; j < d; j++) {
      ILP_ScipCheck(SCIPreleaseVar(scip, &vars[j]), "SCIPreleaseVar");
    }
    ILP_ScipCheck(SCIPfree(&scip), "SCIPfree");
  };
  if (status == SCIP_STATUS_INFEASIBLE) {
    clean();
    return {IlpStatus::Infeasible, {}, T(0)};
  }
  if (status == SCIP_STATUS_UNBOUNDED) {
    clean();
    return {IlpStatus::Unbounded, {}, T(0)};
  }
  if (status != SCIP_STATUS_OPTIMAL) {
    std::cerr << "ILP: SCIP returned the unexpected status "
              << static_cast<int>(status) << "\n";
    clean();
    throw TerminalException{1};
  }
  SCIP_SOL *sol = SCIPgetBestSol(scip);
  MyVector<T> x(d);
  for (int j = 0; j < d; j++) {
    double val = SCIPgetSolVal(scip, sol, vars[j]);
    // The rows are integral and exactly representable, so the solution
    // entries are integers of moderate size.
    long val_i = std::lround(val);
    x(j) = T(val_i);
  }
  clean();
  if (!ILP_IsFeasible(ListIneq, ListEqua, x)) {
    std::cerr << "ILP: The rounded SCIP solution is not feasible\n";
    throw TerminalException{1};
  }
  T value = ILP_EvaluateRow(ToBeMinimized, x);
  return {IlpStatus::Optimal, x, value};
}

// clang-format off
#endif  // SRC_MILP_SCIP_SUPPORT_H_
// clang-format on

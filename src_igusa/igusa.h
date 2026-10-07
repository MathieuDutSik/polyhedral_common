// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_IGUSA_IGUSA_H_
#define SRC_IGUSA_IGUSA_H_

// clang-format off
#include "PerfectForm.h"
#include "POLY_AdjacencyScheme.h"
#include "Tspace_Namelist.h"
#include "Positivity.h"
#include "SignatureSymmetric.h"
#include "integer_linear_programming.h"
#include "GRP_GroupFile.h"
#include <cmath>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_IGUSA
#endif

#ifdef DISABLE_DEBUG_IGUSA
#undef DEBUG_IGUSA
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_IGUSA
#endif

#ifdef TIMINGS
#define TIMINGS_IGUSA
#endif

/*
  The Igusa polyhedron of a T-space.

  Let T be a T-space of symmetric n x n matrices and L the lattice of the
  matrices X of T that are integral valued, that is with X[v] in Z for all
  v in Z^n (integral diagonal and half-integral off-diagonal entries).
  The Igusa polyhedron is

      P = conv(I)   with   I = { X in L : X positive definite }.

  For T the full space of symmetric matrices this is the central cone of
  Igusa. P is contained in the Ryshkov polyhedron { X : X[v] >= 1 } and has
  the same recession cone, the positive semidefinite matrices of T.
  The code enumerates the vertices of P up to the arithmetic group of the
  T-space, and for each vertex A its local cone

      C_A = cone{ X - A : X in I } = cone{ D in L : A + D positive definite }

  with its extreme rays (the edges of P at A) and facets (the facets of P
  containing A).

  An extreme ray D, taken primitive in L, gives the edge [A, A + t D] with
  t the largest integer such that A + t D is positive definite. If D is
  positive semidefinite the edge is infinite; otherwise A + t D is the
  neighboring vertex.

  The minimization of a linear function over I is the basic tool. It is
  done by cutting planes: a linear program over finitely many inequalities
  X[v] >= 1, refined until the optimum is in the Ryshkov polyhedron, then
  an integer program over L whose solution is refined until it is positive
  definite. The inequalities X[v] >= 1 are valid on all of I, so they are
  kept in a pool shared by all the computations.

  The local cone at a vertex A is computed in the following way:
  --- The known edges: the reverse of the edges by which A was reached
      from the vertices already processed.
  --- Candidate rays are the D in L with A + D positive definite and
      tr((A^{-1} D)^2) <= NormBound. They are closed under the stabilizer of
      A and reduced to the extreme ones by linear programming.
  --- Random facets of the cone of candidates are obtained by linear
      programming and checked by integer programming: a facet f is valid
      for C_A iff min_{X in I} f(X) = f(A). A violating X gives a new ray.
      This repeats until NbCleanRound consecutive rounds of NbRandomFacet
      facets find no violation.
  --- Norm closure: the candidates up to the largest norm of the extreme
      rays are added, unless the estimated number of lattice vectors to
      enumerate exceeds MaxNormEnumeration.
  --- Only then the dual description, which is the expensive part, is
      computed, and each of its orbits of facets is checked. With the
      rays found in the previous step it is normally correct at the first
      attempt. If not, the new rays are inserted and the process repeats.
 */

template <typename T> struct IgusaParameters {
  // Candidate rays D with tr((A^{-1} D)^2) <= NormBound
  int NormBound;
  // Number of random facets tested per round
  int NbRandomFacet;
  // Number of consecutive rounds without violation before the dual
  // description is computed
  int NbCleanRound;
  // The integer programming method: default, scip, exact_bb
  std::string IlpMethod;
  // Whether to insert all the candidates up to the largest norm of the
  // extreme rays before the dual description
  bool NormClosure;
  // The norm closure is skipped when the estimated number of lattice
  // vectors to enumerate is above this
  int MaxNormEnumeration;
  // If not "unset", each dual description instance (rays and permutation
  // group) is written to files with this prefix before being computed
  std::string FileInstancePrefix;
  // Stop (by IgusaStopException) after writing the first instance
  bool StopBeforeDualDescription;
  // The heuristics of the dual description (namelist of
  // POLY_SerialDualDesc) or "unset" for igusa_dual_description_heuristic
  std::string FileDualDescription;
  // If positive, a sampling phase of that many random facets, each checked,
  // is done before the dual description
  int NbSampleFacet;
  // StopBeforeDualDescription applies only to the instances with at least
  // that many rays
  int StopMinRays;
  // Whether the candidates and the new rays are inserted with their whole
  // orbit under the stabilizer (the known edges always are). Without it,
  // the group of the cone is taken trivial: for the large stabilizers of
  // dimension 9 and more the orbits cannot be stored.
  bool OrbitClosure;
  // If positive, the bound of the enumeration by norm is not doubled beyond
  // it once the rays span the space
  int MaxEnumerationBound;
  // Write the type of the neighbor of each new ray found by a facet check
  // as soon as it is found
  bool LogNewRays;
};

// Thrown to stop after writing a dual description instance
struct IgusaStopException {};

/*
  The heuristics of the dual descriptions of the local cones. They are the
  standard ones except that for the polytopes with more than 10000
  vertices the bank is neither queried nor filled, and no additional
  symmetries are computed. All of them need a complete weighted graph on the
  vertices: at the vertex D7 (40768 rays) it ran out of memory. The bank is
  empty at the top level anyway.
 */
template <typename T, typename TintGroup>
PolyHeuristicSerial<TintGroup>
igusa_dual_description_heuristic(std::string const &FileDualDescription,
                                 int const &dimEXT, std::ostream &os) {
  if (FileDualDescription != "unset") {
    return Read_AllStandardHeuristicSerial_File<T, TintGroup>(
        FileDualDescription, dimEXT, os);
  }
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, os);
  auto get_heu = [](std::vector<std::string> const &l_str)
      -> DualDescHeuristic<TintGroup> {
    return convert_dual_desc_heuristic(
        HeuristicFromListString<TintGroup>(l_str));
  };
  AllArr.CheckDatabaseBank = get_heu(
      {"2", "1 incidence > 10000 no", "1 incidence > 50 yes", "no"});
  AllArr.BankSave = get_heu({"3", "1 incidence > 10000 no", "1 time < 300 no",
                             "1 incidence < 50 no", "yes"});
  // The additional symmetries are computed on the same kind of graph
  AllArr.AdditionalSymmetry = get_heu(
      {"2", "1 incidence > 10000 no", "1 incidence < 50 no", "yes"});
  return AllArr;
}

template <typename T, typename Tint> struct IgusaSpace {
  LinSpaceMatrix<T> LinSpa;
  int n;
  int dim;
  // A Z-basis of the lattice L of integral valued matrices of the T-space
  std::vector<MyMatrix<T>> ListMatInt;
  // For computing the coordinates in ListMatInt of a matrix of the T-space
  std::vector<int> ListColSel;
  MyMatrix<T> InvSel;
  // The Gram matrix of the trace form on ListMatInt
  MyMatrix<T> TraceGram;
  IgusaParameters<T> params;
  // The pool of the vectors v of the cuts X[v] >= 1 valid on I
  std::vector<MyVector<Tint>> ListCutVect;
  std::set<MyVector<Tint>> SetCutVect;
  MyMatrix<T> ListCutIneq;
};

// The coordinates in which the integral valued forms are the integral
// vectors: the diagonal entries and twice the off-diagonal ones.
template <typename T>
MyMatrix<T> igusa_double_off_diagonal(MyMatrix<T> const &M, T const &scal) {
  int n = M.rows();
  MyMatrix<T> Mret = M;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j) {
        Mret(i, j) = scal * M(i, j);
      }
    }
  }
  return Mret;
}

template <typename T>
std::vector<MyMatrix<T>>
igusa_integral_valued_basis(std::vector<MyMatrix<T>> const &ListMat,
                            std::ostream &os) {
  int n_mat = ListMat.size();
  int n = ListMat[0].rows();
  int sym_dim = n * (n + 1) / 2;
  std::vector<MyMatrix<T>> ListMatRet;
  if (n_mat == sym_dim) {
    for (int i = 0; i < n; i++) {
      for (int j = i; j < n; j++) {
        MyMatrix<T> M = ZeroMatrix<T>(n, n);
        if (i == j) {
          M(i, i) = T(1);
        } else {
          M(i, j) = T(1) / T(2);
          M(j, i) = T(1) / T(2);
        }
        ListMatRet.push_back(M);
      }
    }
    return ListMatRet;
  }
  std::vector<MyMatrix<T>> ListMatDouble;
  for (auto &eMat : ListMat) {
    ListMatDouble.push_back(igusa_double_off_diagonal<T>(eMat, T(2)));
  }
  std::vector<MyMatrix<T>> ListMatDoubleSat =
      IntegralSaturationSpace(ListMatDouble, os);
  for (auto &eMat : ListMatDoubleSat) {
    ListMatRet.push_back(igusa_double_off_diagonal<T>(eMat, T(1) / T(2)));
  }
  return ListMatRet;
}

template <typename T, typename Tint>
IgusaSpace<T, Tint> build_igusa_space(LinSpaceMatrix<T> const &LinSpa,
                                      IgusaParameters<T> const &params,
                                      std::ostream &os) {
  int n = LinSpa.n;
  std::vector<MyMatrix<T>> ListMatInt =
      igusa_integral_valued_basis(LinSpa.ListMat, os);
  int dim = ListMatInt.size();
  int sym_dim = n * (n + 1) / 2;
  MyMatrix<T> BigMat(dim, sym_dim);
  for (int i = 0; i < dim; i++) {
    MyVector<T> V = SymmetricMatrixToVector(ListMatInt[i]);
    AssignMatrixRow(BigMat, i, V);
  }
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(BigMat);
  if (static_cast<int>(eSelect.TheRank) != dim) {
    std::cerr << "IGUSA: The basis of the integral valued forms is not free\n";
    throw TerminalException{1};
  }
  std::vector<int> ListColSel = eSelect.ListColSelect;
  MyMatrix<T> SelMat = SelectColumn(BigMat, ListColSel);
  MyMatrix<T> InvSel = Inverse(SelMat);
  MyMatrix<T> TraceGram(dim, dim);
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      TraceGram(i, j) = frobenius_inner(ListMatInt[i], ListMatInt[j]);
    }
  }
#ifdef DEBUG_IGUSA
  os << "IGUSA: build_igusa_space, n=" << n << " dim=" << dim << "\n";
#endif
  MyMatrix<T> ListCutIneq(0, dim + 1);
  return {LinSpa,     n,         dim,    std::move(ListMatInt),
          ListColSel, InvSel,    TraceGram, params,
          {},         {},        ListCutIneq};
}

template <typename T, typename Tint>
MyVector<T> igusa_coordinates(IgusaSpace<T, Tint> const &space,
                              MyMatrix<T> const &M) {
  MyVector<T> V = SymmetricMatrixToVector(M);
  int dim = space.dim;
  MyVector<T> Vsel(dim);
  for (int i = 0; i < dim; i++) {
    Vsel(i) = V(space.ListColSel[i]);
  }
  MyVector<T> x = space.InvSel.transpose() * Vsel;
#ifdef SANITY_CHECK_IGUSA
  MyMatrix<T> Mrec = GetMatrixFromBasis(space.ListMatInt, x);
  if (Mrec != M) {
    std::cerr << "IGUSA: The matrix is not in the T-space\n";
    throw TerminalException{1};
  }
#endif
  return x;
}

template <typename T, typename Tint>
MyVector<Tint> igusa_integral_coordinates(IgusaSpace<T, Tint> const &space,
                                          MyMatrix<T> const &M) {
  MyVector<T> x = igusa_coordinates(space, M);
  if (!IsIntegralVector(x)) {
    std::cerr << "IGUSA: The matrix is not integral valued\n";
    throw TerminalException{1};
  }
  return UniversalVectorConversion<Tint, T>(x);
}

template <typename T, typename Tint>
MyMatrix<T> igusa_matrix(IgusaSpace<T, Tint> const &space,
                         MyVector<T> const &x) {
  return GetMatrixFromBasis(space.ListMatInt, x);
}

template <typename T, typename Tint>
MyMatrix<T> igusa_matrix_int(IgusaSpace<T, Tint> const &space,
                             MyVector<Tint> const &x) {
  MyVector<T> x_T = UniversalVectorConversion<T, Tint>(x);
  return GetMatrixFromBasis(space.ListMatInt, x_T);
}

// The matrix G of the T-space with tr(G X) = sum_j g(j+1) x_j.
template <typename T, typename Tint>
MyMatrix<T> igusa_functional_matrix(IgusaSpace<T, Tint> const &space,
                                    MyVector<T> const &g) {
  int dim = space.dim;
  MyVector<T> glin(dim);
  for (int j = 0; j < dim; j++) {
    glin(j) = g(j + 1);
  }
  MyVector<T> y = Inverse(space.TraceGram) * glin;
  return igusa_matrix(space, y);
}

// Insert the cut X[v] >= 1 in the pool. Returns false if already present.
template <typename T, typename Tint>
bool igusa_insert_cut(IgusaSpace<T, Tint> &space, MyVector<Tint> const &v) {
  MyVector<Tint> v_can = CanonicalizeVectorToInvertible(v);
  if (space.SetCutVect.count(v_can) > 0) {
    return false;
  }
  space.SetCutVect.insert(v_can);
  space.ListCutVect.push_back(v_can);
  int dim = space.dim;
  int n_row = space.ListCutIneq.rows();
  MyMatrix<T> NewIneq(n_row + 1, dim + 1);
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j <= dim; j++) {
      NewIneq(i, j) = space.ListCutIneq(i, j);
    }
  }
  NewIneq(n_row, 0) = T(-1);
  for (int j = 0; j < dim; j++) {
    NewIneq(n_row, j + 1) =
        EvaluationQuadForm<T, Tint>(space.ListMatInt[j], v_can);
  }
  space.ListCutIneq = NewIneq;
  return true;
}

// The cuts from the vectors of norm at most 1 of a positive definite X,
// and whether there was any vector of norm less than 1.
template <typename T, typename Tint>
bool igusa_insert_short_vector_cuts(IgusaSpace<T, Tint> &space,
                                    MyMatrix<T> const &X, std::ostream &os) {
  std::vector<MyVector<Tint>> l_v =
      computeLevel_GramMat<T, Tint>(X, T(1), os);
  bool has_short = false;
  for (auto &v : l_v) {
    if (EvaluationQuadForm<T, Tint>(X, v) < 1) {
      has_short = true;
    }
    (void)igusa_insert_cut(space, v);
  }
  return has_short;
}

/*
  For X not positive semidefinite (or not positive definite if MaxNorm > 0),
  an integral vector v with X[v] < MaxNorm, for MaxNorm = 0 or 1.
  The rational vectors w with X[w] <= 0 of the diagonalization give one
  exactly once their denominators are cleared, but it can be long. So
  we first look for a short one among the roundings of the multiples of
  the w, with a bounded number of multiples, since for an almost
  degenerate X the unbounded search can take arbitrarily long.
 */
template <typename T, typename Tint>
MyVector<Tint> igusa_non_positive_vector(MyMatrix<T> const &X,
                                         T const &MaxNorm, std::ostream &os) {
  std::vector<MyVector<T>> ListNeg = GetSetNegativeOrZeroVector(X, os);
  std::optional<int> max_mult = 100;
  std::optional<MyVector<Tint>> opt = GetShortVectorSpecifiedBounded<T, Tint>(
      X, ListNeg, MaxNorm, max_mult, os);
  if (opt) {
    return *opt;
  }
  for (auto &w : ListNeg) {
    MyVector<Tint> v =
        UniversalVectorConversion<Tint, T>(RemoveFractionVector(w));
    if (EvaluationQuadForm<T, Tint>(X, v) < MaxNorm) {
      return v;
    }
  }
  std::cerr << "IGUSA: Failed to find a vector of norm below MaxNorm="
            << MaxNorm << "\n";
  throw TerminalException{1};
}

// For X not positive definite, the cut from a vector v with X[v] < 1.
template <typename T, typename Tint>
void igusa_insert_non_positive_cut(IgusaSpace<T, Tint> &space,
                                   MyMatrix<T> const &X, std::ostream &os) {
  MyVector<Tint> v = igusa_non_positive_vector<T, Tint>(X, T(1), os);
  if (!igusa_insert_cut(space, v)) {
    std::cerr << "IGUSA: The cut from the non-positive matrix was already "
                 "present\n";
    throw TerminalException{1};
  }
}

template <typename T> struct IgusaMinimum {
  MyVector<T> x;
  T value;
};

/*
  The inequalities of the cut pool followed by ListExtraIneq, as rows (1, x).
 */
template <typename T, typename Tint>
MyMatrix<T> igusa_cut_inequalities(IgusaSpace<T, Tint> const &space,
                                   MyMatrix<T> const &ListExtraIneq) {
  int dim = space.dim;
  int n_cut = space.ListCutIneq.rows();
  int n_extra = ListExtraIneq.rows();
  MyMatrix<T> M(n_cut + n_extra, dim + 1);
  for (int i = 0; i < n_cut; i++) {
    for (int j = 0; j <= dim; j++) {
      M(i, j) = space.ListCutIneq(i, j);
    }
  }
  for (int i = 0; i < n_extra; i++) {
    for (int j = 0; j <= dim; j++) {
      M(n_cut + i, j) = ListExtraIneq(i, j);
    }
  }
  return M;
}

/*
  Minimize f(x) = f(0) + sum_j f(j) x_j over the real x such that
  X = sum_j x_j ListMatInt[j] satisfies X[v] >= 1 for all v in Z^n - {0},
  the inequalities ListExtraIneq.(1,x) >= 0 and the equalities
  ListEqua.(1,x) = 0. Without extra constraints this is the minimum of f
  over the Ryshkov polyhedron of the T-space. The cuts X[v] >= 1 needed are
  added to the pool of the space. Returns none if infeasible. The function
  f has to be bounded below on the feasible set, which holds if f is
  positive on the nonzero positive semidefinite matrices of the T-space or
  if the equalities bound the feasible set.
 */
template <typename T, typename Tint>
std::optional<IgusaMinimum<T>>
igusa_lp_minimization(IgusaSpace<T, Tint> &space, MyVector<T> const &f,
                      MyMatrix<T> const &ListExtraIneq,
                      MyMatrix<T> const &ListEqua, std::ostream &os) {
  int dim = space.dim;
  auto get_ineq = [&]() -> MyMatrix<T> {
    MyMatrix<T> Mcut = igusa_cut_inequalities(space, ListExtraIneq);
    int n_row = Mcut.rows();
    int n_equa = ListEqua.rows();
    MyMatrix<T> M(n_row + 2 * n_equa, dim + 1);
    for (int i = 0; i < n_row; i++) {
      for (int j = 0; j <= dim; j++) {
        M(i, j) = Mcut(i, j);
      }
    }
    for (int i = 0; i < n_equa; i++) {
      for (int j = 0; j <= dim; j++) {
        M(n_row + 2 * i, j) = ListEqua(i, j);
        M(n_row + 2 * i + 1, j) = -ListEqua(i, j);
      }
    }
    return M;
  };
#ifdef DEBUG_IGUSA
  size_t iter_lp = 0;
#endif
  MyVector<T> x_lp;
  while (true) {
#ifdef DEBUG_IGUSA
    iter_lp++;
#endif
    MyMatrix<T> M = get_ineq();
    LpSolution<T> eSol = SIMPLEX_LinearProgramming(M, f, os);
    if (!eSol.DirectSolution) {
#ifdef DEBUG_IGUSA
      os << "IGUSA: igusa_lp_minimization, LP infeasible at iter_lp="
         << iter_lp << "\n";
#endif
      return {};
    }
    if (!eSol.DualSolution) {
      // Unbounded: the direct solution is a ray d. Since f is positive on
      // the nonzero positive semidefinite matrices and f(d) < 0, D is not
      // positive semidefinite and there is a v with D[v] < 0.
      MyVector<T> ray = *eSol.DirectSolution;
      MyMatrix<T> D = igusa_matrix(space, ray);
      MyVector<Tint> v = igusa_non_positive_vector<T, Tint>(D, T(0), os);
      if (!igusa_insert_cut(space, v)) {
        std::cerr << "IGUSA: The cut for the ray was already present\n";
        throw TerminalException{1};
      }
      continue;
    }
    x_lp = *eSol.DirectSolution;
    MyMatrix<T> X = igusa_matrix(space, x_lp);
    if (!IsPositiveDefinite(X, os)) {
      igusa_insert_non_positive_cut(space, X, os);
      continue;
    }
    if (igusa_insert_short_vector_cuts(space, X, os)) {
      continue;
    }
    break;
  }
#ifdef DEBUG_IGUSA
  os << "IGUSA: igusa_lp_minimization, done iter_lp=" << iter_lp
     << " |cuts|=" << space.ListCutIneq.rows() << "\n";
#endif
  T value = ILP_EvaluateRow(f, x_lp);
  return IgusaMinimum<T>{x_lp, value};
}

/*
  Minimize the function f(x) = f(0) + sum_j f(j) x_j over the x in Z^dim
  such that X = sum_j x_j ListMatInt[j] is positive definite, subject to
  the inequalities ListExtraIneq.(1,x) >= 0 and the equalities
  ListEqua.(1,x) = 0. Returns none if there is no such x.
  ---
  The function f has to be positive on the nonzero positive semidefinite
  matrices of the T-space. Then the part of P where f is below any given
  value is bounded, which is what makes the cutting planes terminate and
  keeps the integer solutions of moderate size. Minimizing a function
  that vanishes on some positive semidefinite direction instead gives
  unbounded sets of optimal solutions of the integer programs, far away
  and almost degenerate.
 */
template <typename T, typename Tint>
std::optional<IgusaMinimum<T>>
igusa_integral_minimization(IgusaSpace<T, Tint> &space, MyVector<T> const &f,
                            MyMatrix<T> const &ListExtraIneq,
                            MyMatrix<T> const &ListEqua, std::ostream &os) {
#ifdef TIMINGS_IGUSA
  MicrosecondTime time;
#endif
  // Without equalities bounding the feasible set, f has to be bounded below
  // on I, that is its matrix has to be positive semidefinite: otherwise f
  // decreases along a ray X + t w w^T, no cut X[v] >= 1 removes it and the
  // cutting planes do not terminate. With equalities (a bounded face of P,
  // as for the coordinate functions minimized on such a face) this is not
  // required.
  if (ListEqua.rows() == 0 &&
      !IsPositiveSemiDefinite(igusa_functional_matrix(space, f), os)) {
    std::cerr << "IGUSA: igusa_integral_minimization, the matrix of the "
                 "function should be positive semidefinite\n";
    throw TerminalException{1};
  }
  // The linear programming phase
  std::optional<IgusaMinimum<T>> opt_lp =
      igusa_lp_minimization(space, f, ListExtraIneq, ListEqua, os);
  if (!opt_lp) {
    return {};
  }
  MyVector<T> const &x_lp = opt_lp->x;
  if (IsIntegralVector(x_lp)) {
#ifdef TIMINGS_IGUSA
    os << "|IGUSA: igusa_integral_minimization(LP)|=" << time << "\n";
#endif
    return opt_lp;
  }
  // The integer programming phase
#ifdef DEBUG_IGUSA
  size_t iter_ilp = 0;
#endif
  while (true) {
#ifdef DEBUG_IGUSA
    iter_ilp++;
#endif
    MyMatrix<T> Mineq = igusa_cut_inequalities(space, ListExtraIneq);
    IlpSolution<T> sol = ILP_IntegerLinearProgramming(
        Mineq, ListEqua, f, space.params.IlpMethod, os);
    if (sol.status == IlpStatus::Infeasible) {
#ifdef DEBUG_IGUSA
      os << "IGUSA: igusa_integral_minimization, ILP infeasible at iter_ilp="
         << iter_ilp << "\n";
#endif
#ifdef TIMINGS_IGUSA
      os << "|IGUSA: igusa_integral_minimization(ILP)|=" << time << "\n";
#endif
      return {};
    }
    if (sol.status == IlpStatus::Unbounded) {
      std::cerr << "IGUSA: The integer program should not be unbounded\n";
      throw TerminalException{1};
    }
    MyMatrix<T> X = igusa_matrix(space, sol.solution);
    if (IsPositiveDefinite(X, os)) {
#ifdef DEBUG_IGUSA
      os << "IGUSA: igusa_integral_minimization, ILP done iter_ilp="
         << iter_ilp << " value=" << sol.OptimalValue << "\n";
#endif
#ifdef TIMINGS_IGUSA
      os << "|IGUSA: igusa_integral_minimization(ILP)|=" << time << "\n";
#endif
      return IgusaMinimum<T>{sol.solution, sol.OptimalValue};
    }
    igusa_insert_non_positive_cut(space, X, os);
  }
}

// A starting set of cuts: X[v] >= 1 for v = e_i and e_i +- e_j.
template <typename T, typename Tint>
void igusa_insert_initial_cuts(IgusaSpace<T, Tint> &space) {
  int n = space.n;
  for (int i = 0; i < n; i++) {
    MyVector<Tint> v = ZeroVector<Tint>(n);
    v(i) = 1;
    (void)igusa_insert_cut(space, v);
    for (int j = i + 1; j < n; j++) {
      MyVector<Tint> w1 = v;
      w1(j) = 1;
      (void)igusa_insert_cut(space, w1);
      MyVector<Tint> w2 = v;
      w2(j) = -1;
      (void)igusa_insert_cut(space, w2);
    }
  }
}

/*
  A vertex of P as the lexicographic minimum over I of
  (tr(SuperMat X), x_1, ..., x_dim). The function tr(SuperMat X) is positive
  on the nonzero positive semidefinite matrices so the face of P where it is
  minimal is a polytope and the lexicographic minimum is a vertex.
 */
template <typename T, typename Tint>
MyMatrix<T> igusa_initial_vertex(IgusaSpace<T, Tint> &space, std::ostream &os) {
  int dim = space.dim;
  MyVector<T> f = ZeroVector<T>(dim + 1);
  for (int j = 0; j < dim; j++) {
    f(j + 1) = frobenius_inner(space.LinSpa.SuperMat, space.ListMatInt[j]);
  }
  std::vector<MyVector<T>> l_equa;
  auto get_equa = [&]() -> MyMatrix<T> {
    return MatrixFromVectorFamilyDim(dim + 1, l_equa);
  };
  MyMatrix<T> ListExtraIneq(0, dim + 1);
  MyVector<T> x;
  // After the first step, the coordinates are minimized on the face where
  // tr(SuperMat X) is minimal, which is bounded.
  for (int k = 0; k <= dim; k++) {
    std::optional<IgusaMinimum<T>> opt =
        igusa_integral_minimization(space, f, ListExtraIneq, get_equa(), os);
    IgusaMinimum<T> result =
        unfold_opt(opt, "The lexicographic minimization should succeed");
    x = result.x;
    MyVector<T> equa = f;
    equa(0) -= result.value;
    l_equa.push_back(equa);
    if (k < dim) {
      f = ZeroVector<T>(dim + 1);
      f(k + 1) = T(1);
    }
  }
  MyMatrix<T> A = igusa_matrix(space, x);
#ifdef DEBUG_IGUSA
  os << "IGUSA: igusa_initial_vertex, A=\n";
  WriteMatrix(os, A);
#endif
  return A;
}

/*
  A test of whether A in I is a vertex of P. The function s(X) =
  tr(A^{-1} X) is minimized over I. If the minimum is s(A), A lies on the
  face F of P where s is minimal, a polytope. Then random generic integral
  functions g are minimized and maximized over F: if the minimum and the
  maximum of g are both g(A), F is the point A, which is then a vertex.
  With random g this is a probabilistic test. With rigorous = true the g
  are the coordinate functions, and F = {A} is then proven.
  Returns true if A is proven to be a vertex. If the minimum of s is below
  s(A), the test is inconclusive for this s.
 */
template <typename T, typename Tint>
bool igusa_vertex_test(IgusaSpace<T, Tint> &space, MyMatrix<T> const &A,
                       int n_direction, bool rigorous, std::ostream &os) {
  int dim = space.dim;
  MyVector<T> a = igusa_coordinates(space, A);
  // The cuts X[v] >= 1 for the minimal vectors of A: without them the first
  // linear programs are far from bounded and the cutting planes are slow.
  (void)igusa_insert_short_vector_cuts(space, A, os);
  MyMatrix<T> Ainv = Inverse(A);
  MyVector<T> s = ZeroVector<T>(dim + 1);
  for (int j = 0; j < dim; j++) {
    s(j + 1) = frobenius_inner(Ainv, space.ListMatInt[j]);
  }
  T s_A = ILP_EvaluateRow(s, a);
  MyMatrix<T> ListExtraIneq(0, dim + 1);
  MyMatrix<T> ListEqua(0, dim + 1);
  std::optional<IgusaMinimum<T>> opt =
      igusa_integral_minimization(space, s, ListExtraIneq, ListEqua, os);
  IgusaMinimum<T> res = unfold_opt(opt, "The minimum of s should exist");
  os << "IGUSA: vertex test, s(A)=" << s_A << " min s over I=" << res.value
     << "\n";
  if (res.value < s_A) {
    os << "IGUSA: vertex test inconclusive: A does not minimize s\n";
    return false;
  }
  MyMatrix<T> Equa(1, dim + 1);
  Equa(0, 0) = -s_A;
  for (int j = 0; j < dim; j++) {
    Equa(0, j + 1) = s(j + 1);
  }
  int n_dir = rigorous ? dim : n_direction;
  for (int i_dir = 0; i_dir < n_dir; i_dir++) {
    MyVector<T> g = ZeroVector<T>(dim + 1);
    if (rigorous) {
      g(i_dir + 1) = T(1);
    } else {
      for (int j = 0; j < dim; j++) {
        g(j + 1) = T(random() % 201 - 100);
      }
    }
    T g_A = ILP_EvaluateRow(g, a);
    for (int sign = -1; sign <= 1; sign += 2) {
      MyVector<T> gs = T(sign) * g;
      std::optional<IgusaMinimum<T>> opt_g =
          igusa_integral_minimization(space, gs, ListExtraIneq, Equa, os);
      IgusaMinimum<T> res_g = unfold_opt(opt_g, "The face is not empty");
      if (res_g.value != T(sign) * g_A) {
        os << "IGUSA: vertex test: the face where s is minimal contains\n";
        WriteMatrix(os, igusa_matrix(space, res_g.x));
        os << "IGUSA: so A is not proven a vertex with this s\n";
        return false;
      }
    }
  }
  os << "IGUSA: vertex test: A is a vertex (exposed by tr(A^{-1} X))";
  if (rigorous) {
    os << ", proven with the coordinate functions";
  }
  os << "\n";
  return true;
}

/*
  The action of the elements g of the stabilizer on the coordinates:
  g ListMatInt[i] g^T = sum_j M(i,j) ListMatInt[j], so the matrix of
  coordinates x is mapped to M^T x.
 */
template <typename T, typename Tint>
MyMatrix<Tint> igusa_coordinate_action(IgusaSpace<T, Tint> const &space,
                                       MyMatrix<Tint> const &g) {
  int dim = space.dim;
  MyMatrix<T> g_T = UniversalMatrixConversion<T, Tint>(g);
  MyMatrix<Tint> M(dim, dim);
  for (int i = 0; i < dim; i++) {
    MyMatrix<T> eImg = g_T * space.ListMatInt[i] * g_T.transpose();
    MyVector<Tint> x = igusa_integral_coordinates(space, eImg);
    AssignMatrixRow(M, i, x);
  }
  return M.transpose();
}

// The orbit index of each point under a permutation group
template <typename Tgroup>
std::vector<int> igusa_orbit_index(Tgroup const &GRP, int n_pt) {
  using Telt = typename Tgroup::Telt;
  std::vector<int> orbit(n_pt, -1);
  std::vector<Telt> l_gens = GRP.GeneratorsOfGroup();
  int i_orb = 0;
  for (int i = 0; i < n_pt; i++) {
    if (orbit[i] == -1) {
      std::vector<int> l_pos{i};
      orbit[i] = i_orb;
      size_t pos = 0;
      while (pos < l_pos.size()) {
        int j = l_pos[pos];
        pos++;
        for (auto &elt : l_gens) {
          int k = elt.at(j);
          if (orbit[k] == -1) {
            orbit[k] = i_orb;
            l_pos.push_back(k);
          }
        }
      }
      i_orb++;
    }
  }
  return orbit;
}

/*
  A facet of the local cone C_A, which is also a facet of P:
  the inequality tr(F X) >= rhs on P, with equality at A.
 */
template <typename T> struct IgusaFacet {
  // The facet functional f in the dual coordinates: f.D >= 0 on C_A
  MyVector<T> f;
  // The matrix F of the T-space with tr(F D) = f.D, primitive
  MyMatrix<T> F;
  // The value tr(F A)
  T rhs;
  // The incidence on the extreme rays of C_A
  Face incd;
};

template <typename T, typename Tint, typename Tgroup> struct IgusaLocalCone {
  // The extreme rays of C_A, primitive in L, in coordinates.
  MyMatrix<Tint> EXT;
  // The permutation action of the stabilizer on the rays
  Tgroup GRP;
  // Representatives of the orbits of rays
  std::vector<int> ListRayRepr;
  // Representatives of the orbits of facets
  std::vector<IgusaFacet<T>> ListFacetRepr;
  // The number of dual descriptions computed, 1 when the rays found before
  // the first one were all the rays.
  int nb_dual_description;
};

template <typename T, typename Tint>
IgusaFacet<T> igusa_get_facet(IgusaSpace<T, Tint> const &space,
                              MyMatrix<T> const &A, MyVector<T> const &f,
                              Face const &incd) {
  MyVector<T> y = Inverse(space.TraceGram) * f;
  MyMatrix<T> F = RemoveFractionMatrix(igusa_matrix(space, y));
  T rhs = frobenius_inner(F, A);
  // f is nonnegative on the positive semidefinite cone and nonzero, so
  // it is positive on A. This fixes the sign of the rescaling.
  if (rhs < 0) {
    F = -F;
    rhs = -rhs;
  }
  return {f, F, rhs, incd};
}

/*
  Returns, if the facet f of the cone of the rays found so far is not a
  facet of C_A, a primitive ray of C_A on which f is negative.
  ---
  The facet is valid iff there is no X in I with f(X) < f(A), that is
  with f(X) <= f(A) - 1 once f is made integral. The function minimized
  over those X is s(X) = tr(A^{-1} X), which is positive on the nonzero
  positive semidefinite matrices (see igusa_integral_minimization). The
  violating X returned is thus the closest to A for s, which is typically
  a neighbor of A.
 */
template <typename T, typename Tint>
std::optional<MyVector<Tint>>
igusa_test_facet(IgusaSpace<T, Tint> &space, MyVector<Tint> const &a,
                 MyVector<T> const &s, MyVector<T> const &f_in,
                 std::ostream &os) {
  int dim = space.dim;
  MyVector<T> f_int = RemoveFractionVector(f_in);
  T val_a(0);
  for (int j = 0; j < dim; j++) {
    val_a += f_int(j) * UniversalScalarConversion<T, Tint>(a(j));
  }
  // The inequality f(A) - 1 - f(X) >= 0
  MyMatrix<T> ListExtraIneq(1, dim + 1);
  ListExtraIneq(0, 0) = val_a - T(1);
  for (int j = 0; j < dim; j++) {
    ListExtraIneq(0, j + 1) = -f_int(j);
  }
  MyMatrix<T> ListEqua(0, dim + 1);
  std::optional<IgusaMinimum<T>> opt =
      igusa_integral_minimization(space, s, ListExtraIneq, ListEqua, os);
  if (!opt) {
    return {};
  }
  MyVector<Tint> x = UniversalVectorConversion<Tint, T>(opt->x);
  MyVector<Tint> diff = x - a;
  return RemoveFractionVector(diff);
}

template <typename T, typename Tint, typename Tgroup>
IgusaLocalCone<T, Tint, Tgroup>
igusa_local_cone(IgusaSpace<T, Tint> &space, MyMatrix<T> const &A,
                 std::vector<MyMatrix<Tint>> const &ListGen,
                 std::vector<MyMatrix<T>> const &ListKnownEdge,
                 std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
#ifdef TIMINGS_IGUSA
  MicrosecondTime time;
#endif
  int dim = space.dim;
  MyVector<Tint> a = igusa_integral_coordinates(space, A);
  std::vector<MyMatrix<Tint>> ListGenCoord;
  for (auto &g : ListGen) {
    ListGenCoord.push_back(igusa_coordinate_action(space, g));
  }
  (void)igusa_insert_short_vector_cuts(space, A, os);
  MyMatrix<T> Ainv = Inverse(A);
  // The function s(X) = tr(A^{-1} X) used for testing the facets
  MyVector<T> s = ZeroVector<T>(dim + 1);
  for (int j = 0; j < dim; j++) {
    s(j + 1) = frobenius_inner(Ainv, space.ListMatInt[j]);
  }
  //
  // The set of candidate rays, closed under the stabilizer.
  //
  std::vector<MyVector<Tint>> ListRay;
  std::map<MyVector<Tint>, int> MapRay;
  auto insert_orbit_kernel = [&](MyVector<Tint> const &x,
                                 bool closure) -> void {
    if (MapRay.count(x) > 0) {
      return;
    }
    std::vector<MyVector<Tint>> l_new{x};
    MapRay[x] = ListRay.size();
    ListRay.push_back(x);
    if (!closure) {
      return;
    }
    size_t pos = 0;
    while (pos < l_new.size()) {
      MyVector<Tint> y = l_new[pos];
      pos++;
      for (auto &M : ListGenCoord) {
        MyVector<Tint> z = M * y;
        if (MapRay.count(z) == 0) {
          MapRay[z] = ListRay.size();
          ListRay.push_back(z);
          l_new.push_back(z);
        }
      }
    }
  };
  auto insert_orbit = [&](MyVector<Tint> const &x) -> void {
    insert_orbit_kernel(x, space.params.OrbitClosure);
  };
  auto is_candidate = [&](MyVector<Tint> const &x) -> bool {
    MyMatrix<T> X = A + igusa_matrix_int(space, x);
    return IsPositiveDefinite(X, os);
  };
  //
  // The candidates from the norm QA(D) = tr((A^{-1} D)^2) = sum lambda_i^2
  // with lambda_i the eigenvalues of A^{-1} D. The condition for D to be a
  // candidate is lambda_i > -1.
  //
  MyMatrix<T> QA(dim, dim);
  {
    std::vector<MyMatrix<T>> ListProd;
    for (int i = 0; i < dim; i++) {
      ListProd.push_back(Ainv * space.ListMatInt[i]);
    }
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        // tr(P_i P_j) for the non-symmetric P_i = A^{-1} ListMatInt[i]
        MyMatrix<T> ProdTr = ListProd[j].transpose();
        QA(i, j) = frobenius_inner(ListProd[i], ProdTr);
      }
    }
  }
  // The number of vectors of the last enumeration and its bound, for
  // estimating the size of the next one.
  size_t last_count = 0;
  T last_bound(0);
  auto enumerate_candidates = [&](T const &bound) -> void {
#ifdef TIMINGS_IGUSA
    MicrosecondTime time_enum;
#endif
    std::vector<MyVector<Tint>> l_x =
        computeLevel_GramMat<T, Tint>(QA, bound, os);
    last_count = l_x.size();
    last_bound = bound;
#ifdef TIMINGS_IGUSA
    os << "|IGUSA: computeLevel_GramMat, bound=" << bound
       << " |l_x|=" << l_x.size() << "|=" << time_enum << "\n";
#endif
    for (auto &x : l_x) {
      for (int sign = -1; sign <= 1; sign += 2) {
        MyVector<Tint> y = RemoveFractionVector(MyVector<Tint>(Tint(sign) * x));
        if (is_candidate(y)) {
          insert_orbit(y);
        }
      }
    }
#ifdef DEBUG_IGUSA
    os << "IGUSA: enumerate_candidates, bound=" << bound
       << " |l_x|=" << l_x.size() << " |ListRay|=" << ListRay.size() << "\n";
#endif
#ifdef TIMINGS_IGUSA
    os << "|IGUSA: candidate test and orbit insertion, |ListRay|="
       << ListRay.size() << "|=" << time_enum << "\n";
#endif
  };
  // The edges already known, typically the reverse of the edges by which
  // A was reached from the vertices processed before. They are genuine
  // extreme rays of C_A, of any norm, so they give the rays that the
  // enumeration by norm misses.
  for (auto &eEdge : ListKnownEdge) {
    MyVector<Tint> x = igusa_integral_coordinates(space, eEdge);
    insert_orbit_kernel(RemoveFractionVector(x), true);
  }
  // The rays of the known edges, and their orbits, are edges of P, hence
  // extreme rays of C_A and of any subcone containing them. They need no
  // redundancy test.
  std::set<MyVector<Tint>> SetCertified(ListRay.begin(), ListRay.end());
#ifdef TIMINGS_IGUSA
  os << "|IGUSA: known edges, |ListKnownEdge|=" << ListKnownEdge.size()
     << " |ListRay|=" << ListRay.size() << "|=" << time << "\n";
#endif
  // Since A is a vertex, there is no D with |lambda_i| < 1 for all i, so
  // no candidate for a bound below 1. The bound is doubled until the
  // candidates span the space.
  // The known edges can span the space alone; the enumeration is still
  // continued until it gives some vectors, since it finds other orbits.
  T current_bound(space.params.NormBound);
  while (true) {
    enumerate_candidates(current_bound);
    bool is_full_rank = false;
    if (!ListRay.empty()) {
      MyMatrix<Tint> EXT = MatrixFromVectorFamilyDim(dim, ListRay);
      is_full_rank = RankMat(EXT) == dim;
    }
    if (is_full_rank && last_count > 0) {
      break;
    }
    if (is_full_rank && space.params.MaxEnumerationBound > 0 &&
        current_bound >= T(space.params.MaxEnumerationBound)) {
      break;
    }
    current_bound *= 2;
  }
  // The neighbors found so far give the scale of the norms of the rays.
  // All the candidates up to the largest norm of the extreme rays are
  // inserted. This is what finds the orbits of rays that the random
  // facets miss, before the dual description is computed. Returns true
  // if the bound was raised.
  auto norm_closure = [&]() -> bool {
    T max_norm(0);
    for (auto &x : ListRay) {
      MyVector<T> x_T = UniversalVectorConversion<T, Tint>(x);
      T norm = EvaluationQuadForm<T, T>(QA, x_T);
      if (norm > max_norm) {
        max_norm = norm;
      }
    }
    if (!space.params.NormClosure || max_norm <= current_bound) {
      return false;
    }
    // The number of lattice vectors of norm at most b grows like
    // b^(dim/2). The closure is skipped when too large: for n = 6 it
    // reached 5 10^7 vectors at the vertex E6 and dominated everything.
    double ratio = UniversalScalarConversion<double, T>(max_norm / last_bound);
    // An empty last enumeration still means at least one vector at the
    // scale of its bound, so the count is taken to be at least 1.
    size_t count = std::max(last_count, static_cast<size_t>(1));
    double estimate =
        static_cast<double>(count) * std::pow(ratio, 0.5 * dim);
#ifdef TIMINGS_IGUSA
    os << "|IGUSA: norm closure, bound=" << max_norm
       << " estimated |l_x|=" << estimate << "|=0\n";
#endif
    if (estimate > static_cast<double>(space.params.MaxNormEnumeration)) {
      return false;
    }
    current_bound = max_norm;
    enumerate_candidates(current_bound);
    return true;
  };
  auto insert_new_ray = [&](MyVector<Tint> const &x) -> void {
#ifdef SANITY_CHECK_IGUSA
    if (!is_candidate(x)) {
      std::cerr << "IGUSA: The new ray should be a candidate\n";
      throw TerminalException{1};
    }
#endif
    insert_orbit(x);
  };
  //
  // The permutation group, the orbits and the redundancy elimination
  //
  auto get_perm_group = [&](std::vector<MyVector<Tint>> const &l_ray,
                            std::map<MyVector<Tint>, int> const &map_ray)
      -> Tgroup {
    int n_ray = l_ray.size();
    std::vector<Telt> l_gens;
    if (!space.params.OrbitClosure) {
      return Tgroup(l_gens, n_ray);
    }
    for (auto &M : ListGenCoord) {
      std::vector<Tidx> eList(n_ray);
      for (int i = 0; i < n_ray; i++) {
        MyVector<Tint> z = M * l_ray[i];
        auto iter = map_ray.find(z);
        if (iter == map_ray.end()) {
          std::cerr << "IGUSA: The set of rays is not invariant\n";
          throw TerminalException{1};
        }
        eList[i] = iter->second;
      }
      l_gens.push_back(Telt(eList));
    }
    return Tgroup(l_gens, n_ray);
  };
  // Replace ListRay by its extreme rays
  // With certified rays, only the orbits of the other rays are tested, one
  // linear program per orbit: the representative r is extreme iff some y
  // has y.g >= 0 for the other rays g and y.r < 0. A non-extreme ray is in
  // the cone of the extreme ones, so removing all the non-extreme orbits
  // keeps the cone.
  auto get_nonredundant_certified =
      [&](MyMatrix<T> const &EXT_T,
          std::vector<int> const &BlockBelong) -> std::vector<int> {
    int n_ray = EXT_T.rows();
    std::vector<uint8_t> keep(n_ray, 1);
    std::map<int, int> orbit_repr;
    for (int i = 0; i < n_ray; i++) {
      if (orbit_repr.count(BlockBelong[i]) == 0) {
        orbit_repr[BlockBelong[i]] = i;
      }
    }
    for (auto &kv : orbit_repr) {
      int i_rep = kv.second;
      if (SetCertified.count(ListRay[i_rep]) > 0) {
        continue;
      }
      // Rows: y.g >= 0 for the kept rays g != r, and 1 + y.r >= 0
      int n_row = 0;
      for (int i = 0; i < n_ray; i++) {
        if (keep[i] && i != i_rep) {
          n_row++;
        }
      }
      MyMatrix<T> ListIneq(n_row + 1, dim + 1);
      int pos = 0;
      for (int i = 0; i < n_ray; i++) {
        if (keep[i] && i != i_rep) {
          ListIneq(pos, 0) = T(0);
          for (int j = 0; j < dim; j++) {
            ListIneq(pos, j + 1) = EXT_T(i, j);
          }
          pos++;
        }
      }
      MyVector<T> obj(dim + 1);
      obj(0) = T(0);
      ListIneq(n_row, 0) = T(1);
      for (int j = 0; j < dim; j++) {
        ListIneq(n_row, j + 1) = EXT_T(i_rep, j);
        obj(j + 1) = EXT_T(i_rep, j);
      }
      LpSolution<T> eSol = SIMPLEX_LinearProgramming(ListIneq, obj, os);
      bool is_extreme = eSol.DirectSolution.has_value() &&
                        eSol.DualSolution.has_value() &&
                        eSol.OptimalValue < 0;
      if (!is_extreme) {
        for (int i = 0; i < n_ray; i++) {
          if (BlockBelong[i] == kv.first) {
            keep[i] = 0;
          }
        }
      }
    }
    std::vector<int> l_idx;
    for (int i = 0; i < n_ray; i++) {
      if (keep[i]) {
        l_idx.push_back(i);
      }
    }
    return l_idx;
  };
  auto reduce_rays = [&]() -> void {
#ifdef TIMINGS_IGUSA
    MicrosecondTime time_red;
#endif
    int n_ray = ListRay.size();
    Tgroup GRP = get_perm_group(ListRay, MapRay);
    std::vector<int> BlockBelong = igusa_orbit_index(GRP, n_ray);
    MyMatrix<Tint> EXT = MatrixFromVectorFamilyDim(dim, ListRay);
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(EXT);
    auto get_l_idx = [&]() -> std::vector<int> {
      if (!SetCertified.empty()) {
        return get_nonredundant_certified(EXT_T, BlockBelong);
      }
      MyMatrix<T> ListIneq(n_ray, dim + 1);
      for (int i = 0; i < n_ray; i++) {
        ListIneq(i, 0) = T(0);
        for (int j = 0; j < dim; j++) {
          ListIneq(i, j + 1) = EXT_T(i, j);
        }
      }
      return SIMPLEX_RedundancyReductionClarksonBlocks(ListIneq, BlockBelong,
                                                       os);
    };
    std::vector<int> l_idx = get_l_idx();
    std::vector<MyVector<Tint>> NewListRay;
    std::map<MyVector<Tint>, int> NewMapRay;
    for (auto &idx : l_idx) {
      NewMapRay[ListRay[idx]] = NewListRay.size();
      NewListRay.push_back(ListRay[idx]);
    }
#ifdef DEBUG_IGUSA
    os << "IGUSA: reduce_rays, " << n_ray << " -> " << NewListRay.size()
       << "\n";
#endif
    ListRay = NewListRay;
    MapRay = NewMapRay;
#ifdef TIMINGS_IGUSA
    os << "|IGUSA: reduce_rays " << n_ray << " -> " << ListRay.size()
       << "|=" << time_red << "\n";
#endif
  };
  auto get_ext_t = [&]() -> MyMatrix<T> {
    MyMatrix<Tint> EXT = MatrixFromVectorFamilyDim(dim, ListRay);
    return UniversalMatrixConversion<T, Tint>(EXT);
  };
  // Test a facet, and insert the violating ray. Returns true if valid.
  auto test_facet = [&](MyVector<T> const &f) -> bool {
#ifdef TIMINGS_IGUSA
    MicrosecondTime time_test;
#endif
    std::optional<MyVector<Tint>> opt = igusa_test_facet(space, a, s, f, os);
#ifdef TIMINGS_IGUSA
    os << "|IGUSA: test_facet valid=" << !opt.has_value() << "|=" << time_test
       << "\n";
#endif
    if (!opt) {
      return true;
    }
    if (space.params.LogNewRays && MapRay.count(*opt) == 0) {
      std::optional<Tint> opt_t = igusa_edge_length(space, A, *opt, os);
      os << "IGUSA: new ray";
      if (opt_t) {
        MyMatrix<T> B = A + UniversalScalarConversion<T, Tint>(*opt_t) *
                                igusa_matrix_int(space, *opt);
        Tshortest<T, Tint> rec_shv = T_ShortestVectorHalf<T, Tint>(B, os);
        os << ", neighbor det=" << DeterminantMat(B)
           << " min=" << rec_shv.min << " pairs=" << rec_shv.SHV.rows()
           << " Neighbor=" << StringMatrixGAP(B);
      } else {
        os << ", infinite edge";
      }
      os << "\n";
    }
    insert_new_ray(*opt);
    return false;
  };
  int nb_dual_description = 0;
  while (true) {
    //
    // The cheap phase with random facets
    //
    int n_clean_round = 0;
    while (true) {
      if (n_clean_round >= space.params.NbCleanRound) {
        reduce_rays();
        if (!norm_closure()) {
          break;
        }
        n_clean_round = 0;
      }
      reduce_rays();
      MyMatrix<T> EXT_T = get_ext_t();
#ifdef TIMINGS_IGUSA
      MicrosecondTime time_fv;
#endif
      vectface vf = FindVertices(EXT_T, space.params.NbRandomFacet, os);
#ifdef TIMINGS_IGUSA
      os << "|IGUSA: FindVertices|=" << time_fv << "\n";
#endif
      SubsetRankOneSolver<T> solver(EXT_T);
      bool is_clean = true;
      for (auto &face : vf) {
        MyVector<T> f = solver.GetPositiveKernelVector(face);
        if (!test_facet(f)) {
          is_clean = false;
        }
      }
#ifdef DEBUG_IGUSA
      os << "IGUSA: random facets, |ListRay|=" << ListRay.size()
         << " is_clean=" << is_clean << "\n";
#endif
      if (is_clean) {
        n_clean_round++;
      } else {
        n_clean_round = 0;
      }
    }
    //
    // The optional sampling phase: NbSampleFacet random facets in total,
    // each checked, the cone being updated after each batch.
    //
    if (space.params.NbSampleFacet > 0 && nb_dual_description == 0) {
      int n_sampled = 0;
#ifdef TIMINGS_IGUSA
      int n_violated = 0;
#endif
      while (n_sampled < space.params.NbSampleFacet) {
        reduce_rays();
        MyMatrix<T> EXT_T = get_ext_t();
        int n_batch = std::min(space.params.NbRandomFacet,
                               space.params.NbSampleFacet - n_sampled);
        vectface vf = FindVertices(EXT_T, n_batch, os);
        SubsetRankOneSolver<T> solver(EXT_T);
        for (auto &face : vf) {
          MyVector<T> f = solver.GetPositiveKernelVector(face);
          bool is_valid = test_facet(f);
#ifdef TIMINGS_IGUSA
          if (!is_valid) {
            n_violated++;
          }
#else
          (void)is_valid;
#endif
        }
        n_sampled += n_batch;
#ifdef TIMINGS_IGUSA
        os << "|IGUSA: sampling, n_sampled=" << n_sampled
           << " n_violated=" << n_violated << " |ListRay|=" << ListRay.size()
           << "|=" << time << "\n";
#endif
      }
    }
    //
    // The dual description and the check of all the facets
    //
    reduce_rays();
    MyMatrix<T> EXT_T = get_ext_t();
    Tgroup GRP = get_perm_group(ListRay, MapRay);
    if (space.params.FileInstancePrefix != "unset") {
      std::string prefix = space.params.FileInstancePrefix + "_" +
                           std::to_string(nb_dual_description + 1);
      WriteMatrixFile(prefix + ".ext", EXT_T);
      WriteGroupFile(prefix + ".grp", GRP);
      WriteGroupFileGAP(prefix + ".grp.gap", GRP);
      // One neighbor B = A + t D per orbit of rays, or the ray D if the
      // edge is infinite, for identifying the neighbors without the dual
      // description.
      {
        int n_ray = ListRay.size();
        std::vector<int> orbit = igusa_orbit_index(GRP, n_ray);
        std::vector<int> orbit_size;
        for (auto &o : orbit) {
          if (o >= static_cast<int>(orbit_size.size())) {
            orbit_size.resize(o + 1, 0);
          }
          orbit_size[o]++;
        }
        std::ofstream os_n(prefix + ".neighbors.g");
        // The types of the neighbors: (det, min, pairs of minimal vectors)
        std::map<std::tuple<T, T, int>, int> map_type;
        os_n << "return [";
        int i_orb_next = 0;
        bool is_first = true;
        for (int i = 0; i < n_ray; i++) {
          if (orbit[i] != i_orb_next) {
            continue;
          }
          if (!is_first) {
            os_n << ",\n";
          }
          is_first = false;
          MyVector<Tint> ray = ListRay[i];
          std::optional<Tint> opt = igusa_edge_length(space, A, ray, os);
          MyMatrix<T> D = igusa_matrix_int(space, ray);
          os_n << "rec(OrbitSize:=" << orbit_size[i_orb_next];
          if (opt) {
            MyMatrix<T> B = A + UniversalScalarConversion<T, Tint>(*opt) * D;
            Tshortest<T, Tint> rec_shv = T_ShortestVectorHalf<T, Tint>(B, os);
            T det = DeterminantMat(B);
            int nb_pair = rec_shv.SHV.rows();
            map_type[{det, rec_shv.min, nb_pair}] += orbit_size[i_orb_next];
            os_n << ", t:=" << *opt << ", det:=" << det
                 << ", min:=" << rec_shv.min << ", nbPairMin:=" << nb_pair
                 << ", Neighbor:=";
            WriteMatrixGAP(os_n, B);
          } else {
            os_n << ", InfiniteRay:=";
            WriteMatrixGAP(os_n, D);
          }
          os_n << ")";
          i_orb_next++;
        }
        os_n << "];\n";
        os << "IGUSA: neighbor types (det, min, pairs, number of rays):";
        for (auto &kv : map_type) {
          os << " (" << std::get<0>(kv.first) << ", "
             << std::get<1>(kv.first) << ", " << std::get<2>(kv.first) << ", "
             << kv.second << ")";
        }
        os << "\n";
      }
      os << "IGUSA: dual description instance written to " << prefix
         << ".ext / .grp / .grp.gap, |EXT|=" << EXT_T.rows()
         << " |GRP|=" << GRP.size() << "\n";
      if (space.params.StopBeforeDualDescription &&
          EXT_T.rows() >= space.params.StopMinRays) {
        throw IgusaStopException{};
      }
    }
#ifdef TIMINGS_IGUSA
    MicrosecondTime time_dd;
#endif
    using TintGroup = typename Tgroup::Tint;
    PolyHeuristicSerial<TintGroup> AllArr =
        igusa_dual_description_heuristic<T, TintGroup>(
            space.params.FileDualDescription, ColumnReduction(EXT_T).cols(),
            os);
    vectface vf = DualDescriptionStandard<T, Tgroup>(EXT_T, GRP, AllArr, os);
    nb_dual_description++;
#ifdef TIMINGS_IGUSA
    os << "|IGUSA: DualDescriptionStandard|=" << time_dd << "\n";
#endif
#ifdef DEBUG_IGUSA
    os << "IGUSA: |EXT|=" << ListRay.size() << " |GRP|=" << GRP.size()
       << " |vf|=" << vf.size() << "\n";
#endif
    SubsetRankOneSolver<T> solver(EXT_T);
    std::vector<IgusaFacet<T>> ListFacetRepr;
    bool is_correct = true;
    for (auto &face : vf) {
      MyVector<T> f = solver.GetPositiveKernelVector(face);
      if (!test_facet(f)) {
        is_correct = false;
      }
      ListFacetRepr.push_back(igusa_get_facet(space, A, f, face));
    }
    if (is_correct) {
      MyMatrix<Tint> EXT = MatrixFromVectorFamilyDim(dim, ListRay);
      int n_ray = ListRay.size();
      std::vector<int> orbit = igusa_orbit_index(GRP, n_ray);
      std::vector<int> ListRayRepr;
      int i_orb_next = 0;
      for (int i = 0; i < n_ray; i++) {
        if (orbit[i] == i_orb_next) {
          ListRayRepr.push_back(i);
          i_orb_next++;
        }
      }
#ifdef TIMINGS_IGUSA
      os << "|IGUSA: igusa_local_cone|=" << time << "\n";
#endif
      return {EXT, GRP, ListRayRepr, ListFacetRepr, nb_dual_description};
    }
#ifdef DEBUG_IGUSA
    os << "IGUSA: The dual description found some violating facets, redoing\n";
#endif
  }
}

/*
  The edge [A, A + t D] for a primitive ray D of C_A: t is the largest
  integer with A + t D positive definite, or none if D is positive
  semidefinite (infinite edge).
 */
template <typename T, typename Tint>
std::optional<Tint> igusa_edge_length(IgusaSpace<T, Tint> const &space,
                                      MyMatrix<T> const &A,
                                      MyVector<Tint> const &x,
                                      std::ostream &os) {
  MyMatrix<T> D = igusa_matrix_int(space, x);
  if (IsPositiveSemiDefinite(D, os)) {
    return {};
  }
  auto is_pd = [&](Tint const &t) -> bool {
    MyMatrix<T> M = A + UniversalScalarConversion<T, Tint>(t) * D;
    return IsPositiveDefinite(M, os);
  };
  // t_low is positive definite, t_upp is not
  Tint t_low(1);
  Tint t_upp(2);
  while (is_pd(t_upp)) {
    t_low = t_upp;
    t_upp *= 2;
  }
  while (t_upp - t_low > 1) {
    Tint t_mid = (t_low + t_upp) / 2;
    if (is_pd(t_mid)) {
      t_low = t_mid;
    } else {
      t_upp = t_mid;
    }
  }
  return t_low;
}

/*
  The enumeration of the vertices of P up to equivalence
 */

template <typename T, typename Tint, typename Tgroup> struct IgusaVertex {
  MyMatrix<T> Gram;
  TshortestPerfect<T, Tint> tsp;
  // Filled by f_adj
  std::vector<MyMatrix<Tint>> GRP_matr;
  typename Tgroup::Tint stab_size;
  MyMatrix<Tint> EXT;
  // The action of the stabilizer on the rows of EXT
  Tgroup GRP;
  std::vector<IgusaFacet<T>> ListFacetRepr;
  // The infinite edges up to the stabilizer, as matrices
  std::vector<MyMatrix<T>> ListInfiniteRay;
  int nb_dual_description;
};

template <typename T, typename Tint> struct IgusaVertex_AdjI {
  MyMatrix<T> Gram;
  TshortestPerfect<T, Tint> tsp;
  // The direction of the edge
  MyMatrix<T> Direction;
};

template <typename T, typename Tint> struct IgusaVertex_AdjO {
  MyMatrix<T> Direction;
  MyMatrix<Tint> eBigMat;
};

template <typename T, typename Tint, typename Tgroup> struct DataIgusaFunc {
  IgusaSpace<T, Tint> space;
  std::ostream &os;
  // For each vertex representative, the reverse of the edges by which it
  // was reached, in its frame.
  std::map<MyMatrix<T>, std::vector<MyMatrix<T>>> MapKnownEdge;
  // If set, the enumeration starts from that vertex instead of the
  // lexicographic minimum
  std::optional<MyMatrix<T>> InitialVertex;
  using Tobj = IgusaVertex<T, Tint, Tgroup>;
  using TadjI = IgusaVertex_AdjI<T, Tint>;
  using TadjO = IgusaVertex_AdjO<T, Tint>;
  std::ostream &get_os() { return os; }

  Tobj make_vertex(MyMatrix<T> const &Gram) {
    Tshortest<T, Tint> rec_shv = T_ShortestVectorHalf<T, Tint>(Gram, os);
    TshortestPerfect<T, Tint> tsp =
        build_tshortest_perfect<T, Tint>(Gram, rec_shv, os);
    return {Gram, std::move(tsp), {}, 0, {}, {}, {}, {}, 0};
  }

  Tobj f_init() {
    igusa_insert_initial_cuts(space);
    if (InitialVertex) {
      return make_vertex(*InitialVertex);
    }
    MyMatrix<T> A = igusa_initial_vertex(space, os);
    return make_vertex(A);
  }

  size_t f_hash(size_t const &seed, Tobj const &x) {
    return SimplePerfect_Invariant<T, Tint>(seed, space.LinSpa, x.Gram, x.tsp,
                                            os);
  }

  std::optional<TadjO> f_repr(Tobj const &x, TadjI const &y) {
    std::optional<MyMatrix<Tint>> opt =
        SimplePerfect_TestEquivalence<T, Tint, Tgroup>(
            space.LinSpa, x.Gram, y.Gram, x.tsp, y.tsp, os);
    if (!opt) {
      return {};
    }
    // M x.Gram M^T = y.Gram, and the reverse edge at y.Gram is -Direction
    MyMatrix<T> Minv =
        UniversalMatrixConversion<T, Tint>(Inverse(*opt));
    MyMatrix<T> eEdge = -Minv * y.Direction * Minv.transpose();
    MapKnownEdge[x.Gram].push_back(eEdge);
    TadjO ret{y.Direction, *opt};
    return ret;
  }

  std::pair<Tobj, TadjO> f_spann(TadjI const &y) {
    Tobj x{y.Gram, y.tsp, {}, 0, {}, {}, {}, {}, 0};
    MapKnownEdge[y.Gram].push_back(-y.Direction);
    TadjO ret{y.Direction, IdentityMat<Tint>(space.n)};
    return {x, ret};
  }

  std::optional<std::vector<TadjI>> f_adj(Tobj &x) {
    std::pair<Tgroup, std::vector<MyMatrix<Tint>>> pair =
        SimplePerfect_Stabilizer<T, Tint, Tgroup>(space.LinSpa, x.Gram, x.tsp,
                                                  os);
    x.stab_size = pair.first.size();
    x.GRP_matr = pair.second;
    IgusaLocalCone<T, Tint, Tgroup> cone =
        igusa_local_cone<T, Tint, Tgroup>(space, x.Gram, x.GRP_matr,
                                          MapKnownEdge[x.Gram], os);
    x.EXT = cone.EXT;
    x.GRP = cone.GRP;
    x.ListFacetRepr = cone.ListFacetRepr;
    x.nb_dual_description = cone.nb_dual_description;
    std::vector<TadjI> ListAdj;
    for (auto &i_ray : cone.ListRayRepr) {
      MyVector<Tint> ray = GetMatrixRow(cone.EXT, i_ray);
      MyMatrix<T> D = igusa_matrix_int(space, ray);
      std::optional<Tint> opt = igusa_edge_length(space, x.Gram, ray, os);
      if (!opt) {
        x.ListInfiniteRay.push_back(D);
        continue;
      }
      MyMatrix<T> Direction = UniversalScalarConversion<T, Tint>(*opt) * D;
      MyMatrix<T> B = x.Gram + Direction;
      Tobj y = make_vertex(B);
      ListAdj.push_back({B, y.tsp, Direction});
    }
#ifdef DEBUG_IGUSA
    os << "IGUSA: f_adj, |EXT|=" << x.EXT.rows()
       << " |ListRayRepr|=" << cone.ListRayRepr.size()
       << " |ListAdj|=" << ListAdj.size()
       << " |ListInfiniteRay|=" << x.ListInfiniteRay.size()
       << " |ListFacetRepr|=" << x.ListFacetRepr.size() << "\n";
#endif
    return ListAdj;
  }

  Tobj f_adji_obj(TadjI const &x) {
    return {x.Gram, x.tsp, {}, 0, {}, {}, {}, {}, 0};
  }
};

template <typename T, typename Tint, typename Tgroup>
void WriteEntryGAP(std::ostream &os_out,
                   IgusaVertex<T, Tint, Tgroup> const &obj) {
  os_out << "rec(Gram:=";
  WriteMatrixGAP(os_out, obj.Gram);
  os_out << ", GRPsize:=" << obj.stab_size;
  os_out << ", nbRay:=" << obj.EXT.rows();
  os_out << ", nbDualDescription:=" << obj.nb_dual_description;
  os_out << ", ListFacet:=[";
  bool is_first = true;
  for (auto &facet : obj.ListFacetRepr) {
    if (!is_first) {
      os_out << ",";
    }
    is_first = false;
    os_out << "rec(F:=";
    WriteMatrixGAP(os_out, facet.F);
    os_out << ", rhs:=" << facet.rhs;
    os_out << ", nbIncd:=" << facet.incd.count();
    os_out << ", OrbitSize:="
           << obj.GRP.size() / obj.GRP.Stabilizer_OnSets(facet.incd).size()
           << ")";
  }
  os_out << "], ListInfiniteRay:=";
  WriteListMatrixGAP(os_out, obj.ListInfiniteRay);
  os_out << ")";
}

template <typename T, typename Tint>
void WriteEntryGAP(std::ostream &os_out,
                   IgusaVertex_AdjO<T, Tint> const &adj) {
  os_out << "rec(Direction:=";
  WriteMatrixGAP(os_out, adj.Direction);
  os_out << ", eBigMat:=";
  WriteMatrixGAP(os_out, adj.eBigMat);
  os_out << ")";
}

template <typename T, typename Tint, typename Tgroup>
void WriteEntryPYTHON(std::ostream &os_out,
                      IgusaVertex<T, Tint, Tgroup> const &obj) {
  os_out << "{\"Gram\":" << StringMatrixPYTHON(obj.Gram);
  os_out << ", \"GRPsize\":" << obj.stab_size;
  os_out << ", \"nbRay\":" << obj.EXT.rows();
  os_out << ", \"nbDualDescription\":" << obj.nb_dual_description;
  os_out << ", \"ListFacet\":[";
  bool is_first = true;
  for (auto &facet : obj.ListFacetRepr) {
    if (!is_first) {
      os_out << ",";
    }
    is_first = false;
    os_out << "{\"F\":" << StringMatrixPYTHON(facet.F)
           << ", \"rhs\":\"" << facet.rhs << "\"}";
  }
  os_out << "], \"ListInfiniteRay\":[";
  is_first = true;
  for (auto &M : obj.ListInfiniteRay) {
    if (!is_first) {
      os_out << ",";
    }
    is_first = false;
    os_out << StringMatrixPYTHON(M);
  }
  os_out << "]}";
}

template <typename T, typename Tint>
void WriteEntryPYTHON(std::ostream &os_out,
                      IgusaVertex_AdjO<T, Tint> const &adj) {
  os_out << "{\"Direction\":" << StringMatrixPYTHON(adj.Direction);
  os_out << ", \"eBigMat\":" << StringMatrixPYTHON(adj.eBigMat) << "}";
}

/*
  The orbits of facets of P. A facet of P appears at each of its vertices
  as a facet of the local cone, so the same orbit can appear at several
  vertex representatives, or several times at one representative (with
  different Stab(A)-orbits). Moving along an edge [A, B] contained in the
  facet and mapping B to its representative identifies the facet with a
  facet at the representative of B. Since the graph of the vertices of a
  facet is connected, the union-find over those moves gives the orbits.
 */
template <typename T> struct IgusaFacetOrbit {
  MyMatrix<T> F;
  T rhs;
  int rank;
  // The (vertex, facet) representatives in the orbit
  std::vector<std::pair<int, int>> l_repr;
};

template <typename T, typename Tint, typename Tgroup>
std::vector<IgusaFacetOrbit<T>> igusa_facet_orbits(
    DataIgusaFunc<T, Tint, Tgroup> &data,
    std::vector<IgusaVertex<T, Tint, Tgroup>> const &l_vert) {
  std::ostream &os = data.os;
  IgusaSpace<T, Tint> const &space = data.space;
  int n_vert = l_vert.size();
  std::vector<int> l_shift(n_vert + 1, 0);
  for (int i = 0; i < n_vert; i++) {
    l_shift[i + 1] = l_shift[i] + l_vert[i].ListFacetRepr.size();
  }
  int n_pair = l_shift[n_vert];
  std::vector<int> parent(n_pair);
  for (int i = 0; i < n_pair; i++) {
    parent[i] = i;
  }
  auto find = [&](int i) -> int {
    while (parent[i] != i) {
      parent[i] = parent[parent[i]];
      i = parent[i];
    }
    return i;
  };
  std::vector<size_t> l_hash;
  size_t seed = 1234;
  for (auto &v : l_vert) {
    l_hash.push_back(data.f_hash(seed, v));
  }
  std::vector<std::vector<Face>> l_can;
  for (auto &v : l_vert) {
    std::vector<Face> l_face;
    for (auto &facet : v.ListFacetRepr) {
      l_face.push_back(v.GRP.CanonicalImage(facet.incd));
    }
    l_can.push_back(l_face);
  }
  for (int i_vert = 0; i_vert < n_vert; i_vert++) {
    IgusaVertex<T, Tint, Tgroup> const &v = l_vert[i_vert];
    int n_ray = v.EXT.rows();
    int n_facet = v.ListFacetRepr.size();
    for (int i_facet = 0; i_facet < n_facet; i_facet++) {
      IgusaFacet<T> const &facet = v.ListFacetRepr[i_facet];
      // One ray for each orbit of the stabilizer of the facet
      Tgroup GRPfacet = v.GRP.Stabilizer_OnSets(facet.incd);
      std::vector<int> orbit = igusa_orbit_index(GRPfacet, n_ray);
      std::set<int> set_orbit_done;
      for (int i_ray = 0; i_ray < n_ray; i_ray++) {
        if (facet.incd[i_ray] == 0 || set_orbit_done.count(orbit[i_ray]) > 0) {
          continue;
        }
        set_orbit_done.insert(orbit[i_ray]);
        MyVector<Tint> ray = GetMatrixRow(v.EXT, i_ray);
        std::optional<Tint> opt = igusa_edge_length(space, v.Gram, ray, os);
        if (!opt) {
          continue;
        }
        MyMatrix<T> B = v.Gram + UniversalScalarConversion<T, Tint>(*opt) *
                                     igusa_matrix_int(space, ray);
        IgusaVertex<T, Tint, Tgroup> vB = data.make_vertex(B);
        size_t hashB = data.f_hash(seed, vB);
        bool is_found = false;
        for (int j_vert = 0; j_vert < n_vert; j_vert++) {
          if (l_hash[j_vert] != hashB) {
            continue;
          }
          IgusaVertex<T, Tint, Tgroup> const &w = l_vert[j_vert];
          std::optional<MyMatrix<Tint>> optP =
              SimplePerfect_TestEquivalence<T, Tint, Tgroup>(
                  space.LinSpa, B, w.Gram, vB.tsp, w.tsp, os);
          if (!optP) {
            continue;
          }
          is_found = true;
          // P B P^T = w.Gram, so f'(X) = f(P^{-1} X P^{-T})
          MyMatrix<Tint> Pinv = Inverse(*optP);
          MyMatrix<Tint> M = igusa_coordinate_action(space, Pinv);
          MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(M);
          MyVector<T> f_img = M_T.transpose() * facet.f;
          int n_ray_w = w.EXT.rows();
          Face incd(n_ray_w);
          for (int j_ray = 0; j_ray < n_ray_w; j_ray++) {
            MyVector<T> ray_w =
                UniversalVectorConversion<T, Tint>(GetMatrixRow(w.EXT, j_ray));
            T scal = f_img.dot(ray_w);
#ifdef SANITY_CHECK_IGUSA
            if (scal < 0) {
              std::cerr << "IGUSA: The mapped facet should be valid\n";
              throw TerminalException{1};
            }
#endif
            if (scal == 0) {
              incd[j_ray] = 1;
            }
          }
          Face can = w.GRP.CanonicalImage(incd);
          int j_facet = -1;
          for (size_t u = 0; u < l_can[j_vert].size(); u++) {
            if (l_can[j_vert][u] == can) {
              j_facet = u;
            }
          }
          if (j_facet == -1) {
            std::cerr << "IGUSA: Failed to find the mapped facet\n";
            throw TerminalException{1};
          }
          int a = find(l_shift[i_vert] + i_facet);
          int b = find(l_shift[j_vert] + j_facet);
          parent[a] = b;
          break;
        }
        if (!is_found) {
          std::cerr << "IGUSA: Failed to find the neighbor among the "
                       "vertices\n";
          throw TerminalException{1};
        }
      }
    }
  }
  std::map<int, int> map_orbit;
  std::vector<IgusaFacetOrbit<T>> l_orbit;
  for (int i_vert = 0; i_vert < n_vert; i_vert++) {
    int n_facet = l_vert[i_vert].ListFacetRepr.size();
    for (int i_facet = 0; i_facet < n_facet; i_facet++) {
      int root = find(l_shift[i_vert] + i_facet);
      if (map_orbit.count(root) == 0) {
        IgusaFacet<T> const &facet = l_vert[i_vert].ListFacetRepr[i_facet];
        map_orbit[root] = l_orbit.size();
        int rank = RankMat(facet.F);
        l_orbit.push_back({facet.F, facet.rhs, rank, {}});
      }
      l_orbit[map_orbit[root]].l_repr.push_back({i_vert, i_facet});
    }
  }
  return l_orbit;
}

template <typename T>
void WriteEntryGAP(std::ostream &os_out, IgusaFacetOrbit<T> const &orb) {
  os_out << "rec(F:=";
  WriteMatrixGAP(os_out, orb.F);
  os_out << ", rhs:=" << orb.rhs << ", rank:=" << orb.rank
         << ", ListVertexFacet:=[";
  bool is_first = true;
  for (auto &pair : orb.l_repr) {
    if (!is_first) {
      os_out << ",";
    }
    is_first = false;
    os_out << "[" << (pair.first + 1) << "," << (pair.second + 1) << "]";
  }
  os_out << "])";
}

template <typename T>
void WriteEntryPYTHON(std::ostream &os_out, IgusaFacetOrbit<T> const &orb) {
  os_out << "{\"F\":" << StringMatrixPYTHON(orb.F) << ", \"rhs\":\""
         << orb.rhs << "\", \"rank\":" << orb.rank
         << ", \"ListVertexFacet\":[";
  bool is_first = true;
  for (auto &pair : orb.l_repr) {
    if (!is_first) {
      os_out << ",";
    }
    is_first = false;
    os_out << "[" << pair.first << "," << pair.second << "]";
  }
  os_out << "]}";
}

/*
  For the analysis of the complexity: the dual description in the other
  direction, from all the facets of the local cone to its extreme rays.
  Returns the number of orbits of extreme rays found, which should be the
  number of orbits of rays of the vertex.
 */
template <typename T, typename Tint, typename Tgroup>
size_t igusa_dual_side_description(IgusaVertex<T, Tint, Tgroup> const &x,
                                   std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n_ray = x.EXT.rows();
  std::vector<Telt> l_gens = x.GRP.GeneratorsOfGroup();
  std::vector<Face> l_face;
  std::map<Face, int> map_face;
  for (auto &facet : x.ListFacetRepr) {
    if (map_face.count(facet.incd) > 0) {
      continue;
    }
    size_t pos = l_face.size();
    map_face[facet.incd] = pos;
    l_face.push_back(facet.incd);
    while (pos < l_face.size()) {
      Face f = l_face[pos];
      pos++;
      for (auto &elt : l_gens) {
        Face g(n_ray);
        for (int i = 0; i < n_ray; i++) {
          if (f[i] == 1) {
            g[elt.at(i)] = 1;
          }
        }
        if (map_face.count(g) == 0) {
          map_face[g] = l_face.size();
          l_face.push_back(g);
        }
      }
    }
  }
  int n_fac = l_face.size();
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(x.EXT);
  SubsetRankOneSolver<T> solver(EXT_T);
  int dim = EXT_T.cols();
  MyMatrix<T> FAC(n_fac, dim);
  for (int i = 0; i < n_fac; i++) {
    MyVector<T> f = solver.GetPositiveKernelVector(l_face[i]);
    AssignMatrixRow(FAC, i, f);
  }
  std::vector<Telt> l_gens_fac;
  for (auto &elt : l_gens) {
    std::vector<Tidx> eList(n_fac);
    for (int i = 0; i < n_fac; i++) {
      Face g(n_ray);
      for (int j = 0; j < n_ray; j++) {
        if (l_face[i][j] == 1) {
          g[elt.at(j)] = 1;
        }
      }
      eList[i] = map_face.at(g);
    }
    l_gens_fac.push_back(Telt(eList));
  }
  Tgroup GRPfac(l_gens_fac, n_fac);
  MicrosecondTime time;
  vectface vf = DualDescriptionStandard<T, Tgroup>(FAC, GRPfac, os);
  os << "IGUSA: dual side, |FAC|=" << n_fac << " |EXT|=" << n_ray
     << " orbits of rays found=" << vf.size() << " time=" << time << "\n";
  return vf.size();
}

inline FullNamelist NAMELIST_GetStandard_ENUMERATE_IGUSA_TSPACE() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  {
    std::map<std::string, std::string> ListStringValues;
    std::map<std::string, int> ListIntValues;
    ListStringValues["arithmetic"] = "gmp";
    ListStringValues["IlpMethod"] = "default";
    ListStringValues["OutFormat"] = "GAP";
    ListStringValues["OutFile"] = "stderr";
    // If set, only the local cone of the vertex in that file is computed
    ListStringValues["FileVertex"] = "unset";
    ListIntValues["NormBound"] = 2;
    ListIntValues["NbRandomFacet"] = 10;
    ListIntValues["NbCleanRound"] = 2;
    ListIntValues["MaxNormEnumeration"] = 1000000;
    // If positive, that many random facets are sampled and checked before
    // the dual description
    ListIntValues["NbSampleFacet"] = 0;
    // StopBeforeDualDescription only for instances with that many rays
    ListIntValues["StopMinRays"] = 0;
    // Stop doubling the bound of the enumeration by norm beyond this once
    // the rays span the space (0: no limit)
    ListIntValues["MaxEnumerationBound"] = 0;
    // Write each dual description instance with this prefix
    ListStringValues["FileInstancePrefix"] = "unset";
    // With FileVertex, the edges known at that vertex (list of matrices)
    ListStringValues["FileKnownEdges"] = "unset";
    // If set, the enumeration starts from that vertex (it must be a vertex)
    ListStringValues["FileInitialVertex"] = "unset";
    // The heuristics of the dual description, a namelist file in the format
    // of POLY_SerialDualDesc, or unset for the standard ones without the
    // bank above 10000 vertices
    ListStringValues["FileDualDescription"] = "unset";
    ListIntValues["max_runtime_second"] = 0;
    std::map<std::string, bool> ListBoolValues;
    ListBoolValues["NormClosure"] = true;
    // In the FileVertex mode, also time the dual description from the
    // facets to the rays, for the analysis of the complexity
    ListBoolValues["CompareDualSide"] = false;
    // Stop after writing the first dual description instance
    ListBoolValues["StopBeforeDualDescription"] = false;
    // With FileVertex, only test whether the vertex is a vertex of P
    ListBoolValues["VertexTest"] = false;
    // Insert the candidates and the new rays with their whole orbit
    ListBoolValues["OrbitClosure"] = true;
    // Write the neighbor type of each new ray found by a facet check
    ListBoolValues["LogNewRays"] = false;
    // The face test with the coordinate functions instead of random ones
    ListBoolValues["VertexTestRigorous"] = false;
    SingleBlock BlockDATA;
    BlockDATA.setListStringValues(ListStringValues);
    BlockDATA.setListIntValues(ListIntValues);
    BlockDATA.setListBoolValues(ListBoolValues);
    ListBlock["DATA"] = BlockDATA;
  }
  // TSPACE
  ListBlock["TSPACE"] = SINGLEBLOCK_Get_Tspace_Description();
  // Merging all data
  return FullNamelist(ListBlock);
}

// clang-format off
#endif  // SRC_IGUSA_IGUSA_H_
// clang-format on

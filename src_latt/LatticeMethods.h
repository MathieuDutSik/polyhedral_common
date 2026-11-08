// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEMETHODS_H_
#define SRC_LATT_LATTICEMETHODS_H_

// clang-format off
#include "LatticeStabEquiCan.h"
#include "LatticeOrbitSplitCanonical.h"
#include "LatticeRootDecomposition.h"
#include "LatticePleskenSouvignier.h"
#include "LatticeAutomorphismVinberg.h"
#include "ClassicLLL.h"
#include <optional>
#include <string>
#include <utility>
#include <vector>
// clang-format on

/*
  The methods for the three basic computations on a positive definite form:
  its canonical form, the isometry test of two forms, and the automorphism
  group. LATT_Canonicalize, LATT_Isomorphism, LATT_Automorphism and the genus
  enumeration choose among them by name, through ListLatticeStabCanIsomMethod,
  and call the dispatchers below, so that a method is described and checked
  in one place.

  Most methods do the three computations. weyl differs from vinberg only in
  its canonical form: for the isometry and the automorphism group the
  residual family of the decomposition is already the Weyl-reduced one, so
  they are the vinberg computations. iterated_shortest gives only the
  isometry and the automorphism group. The
  root-system methods (vinberg, weyl, weyl_chamber) need an integral form,
  whose roots are its vectors of norm 2, and for the isometry and the
  automorphism group a single form. The canonical forms of two different
  methods are not comparable with each other.

  The dispatchers take a Gram matrix in any basis and LLL-reduce it once
  (the first form of the list for the isometry and the group); the methods
  behind them work on the reduced matrix and assume it is reduced. The
  results are given back in the basis of the input.
 */
enum class LatticeStabCanIsomMethod {
  Graph,
  PleskenSouvignier,
  Vinberg,
  Weyl,
  WeylChamber,
  IteratedShortest
};

// A method: its name and a description for the usage and namelist texts.
struct LatticeStabCanIsomMethodEntry {
  std::string name;
  LatticeStabCanIsomMethod method;
  std::string description;
};

inline std::vector<LatticeStabCanIsomMethodEntry>
ListLatticeStabCanIsomMethod() {
  return {
      {"graph", LatticeStabCanIsomMethod::Graph,
       "from the characteristic vector family"},
      {"plesken_souvignier", LatticeStabCanIsomMethod::PleskenSouvignier,
       "the Plesken-Souvignier backtrack on the Z-spanning family; the "
       "canonical form by orbit splitting over its automorphism group, "
       "viable on large families (Leech), can decline on some inputs"},
      {"vinberg", LatticeStabCanIsomMethod::Vinberg,
       "along the root decomposition, the automorphism group as W(R) rtimes "
       "Aut(L, rho); integral form"},
      {"weyl", LatticeStabCanIsomMethod::Weyl,
       "as vinberg, with only the simple roots in the family; the isometry "
       "and the automorphism group are those of vinberg"},
      {"weyl_chamber", LatticeStabCanIsomMethod::WeylChamber,
       "the characteristic vectors in the Weyl chamber and the simple roots, "
       "a family spanning the lattice; integral form"},
      {"iterated_shortest", LatticeStabCanIsomMethod::IteratedShortest,
       "the MatrixGroup algorithm on the full-rank iterated-shortest family; "
       "isometry and automorphism group only"}};
}

inline LatticeStabCanIsomMethod
ParseLatticeStabCanIsomMethod(std::string const &name) {
  for (auto &entry : ListLatticeStabCanIsomMethod()) {
    if (entry.name == name) {
      return entry.method;
    }
  }
  std::cerr << "LATTICE_METHODS: unknown method " << name
            << ". The methods are";
  for (auto &entry : ListLatticeStabCanIsomMethod()) {
    std::cerr << " " << entry.name;
  }
  std::cerr << "\n";
  throw TerminalException{1};
}

inline std::string NameLatticeStabCanIsomMethod(LatticeStabCanIsomMethod method) {
  for (auto &entry : ListLatticeStabCanIsomMethod()) {
    if (entry.method == method) {
      return entry.name;
    }
  }
  std::cerr << "LATTICE_METHODS: a method without a name\n";
  throw TerminalException{1};
}

// The names and descriptions, one per line, for the usage texts.
inline void PrintLatticeStabCanIsomMethods(std::ostream &os) {
  for (auto &entry : ListLatticeStabCanIsomMethod()) {
    os << "  " << entry.name << " : " << entry.description << "\n";
  }
}

// The same as one string, for a namelist documentation.
inline std::string DocLatticeStabCanIsomMethods() {
  std::string doc;
  for (auto &entry : ListLatticeStabCanIsomMethod()) {
    doc += "\n" + entry.name + ": " + entry.description;
  }
  return doc;
}

[[noreturn]] inline void
LatticeStabCanIsomMethodUnavailable(LatticeStabCanIsomMethod method,
                                    std::string const &task) {
  std::cerr << "LATTICE_METHODS: the method "
            << NameLatticeStabCanIsomMethod(method) << " does not compute "
            << task << "\n";
  throw TerminalException{1};
}

template <typename T>
void CheckLatticeMethodIntegral(std::vector<MyMatrix<T>> const &ListMat,
                                std::string const &name) {
  if (ListMat.size() != 1 || !IsIntegralMatrix(ListMat[0])) {
    std::cerr << "LATTICE_METHODS: the " << name << " method needs a single "
              << "integral Gram matrix\n";
    throw TerminalException{1};
  }
}

// The forms of ListMat in the basis Pmat, rows: Pmat M Pmat^T for each M.
template <typename T, typename Tint>
std::vector<MyMatrix<T>>
LatticeMethodsChangeBasis(std::vector<MyMatrix<T>> const &ListMat,
                          MyMatrix<Tint> const &Pmat) {
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  std::vector<MyMatrix<T>> ListMatRed;
  for (auto &M : ListMat) {
    ListMatRed.push_back(Pmat_T * M * Pmat_T.transpose());
  }
  return ListMatRed;
}

/*
  The canonical form of the form GramMat by the method: a basis B, as rows,
  with B GramMat B^T canonical. LatticeCanonicalBasisReduced assumes GramMat
  LLL-reduced.
 */
template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint> LatticeCanonicalBasisReduced(MyMatrix<T> const &GramMat,
                                            LatticeStabCanIsomMethod method,
                                            std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat{GramMat};
  std::string name = NameLatticeStabCanIsomMethod(method);
  switch (method) {
  case LatticeStabCanIsomMethod::Graph:
    return ComputeCanonicalForm<T, Tint, Tgroup>(GramMat, os);
  case LatticeStabCanIsomMethod::PleskenSouvignier: {
    // The failure is a property of the isometry class, so the forms the
    // method does canonicalize remain comparable with each other.
    std::optional<MyMatrix<Tint>> opt =
        OrbitSplitCanonicalizeGram<T, Tint, Tgroup>(GramMat, os);
    if (!opt) {
      std::cerr << "LATTICE_METHODS: the orbit splitting could not select a "
                   "canonical prefix, no canonical form for this input\n";
      throw TerminalException{1};
    }
    return *opt;
  }
  case LatticeStabCanIsomMethod::Vinberg:
    CheckLatticeMethodIntegral(ListMat, name);
    return ComputeCanonicalFormRootDecomposed<T, Tint, Tgroup>(GramMat, os);
  case LatticeStabCanIsomMethod::Weyl:
    CheckLatticeMethodIntegral(ListMat, name);
    return ComputeCanonicalFormWeylReduced<T, Tint, Tgroup>(GramMat, os);
  case LatticeStabCanIsomMethod::WeylChamber:
    CheckLatticeMethodIntegral(ListMat, name);
    return ComputeCanonicalFormWeylChamber<T, Tint, Tgroup>(GramMat, os);
  case LatticeStabCanIsomMethod::IteratedShortest:
    break;
  }
  LatticeStabCanIsomMethodUnavailable(method, "a canonical form");
}

template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint> LatticeCanonicalBasis(MyMatrix<T> const &GramMat,
                                     LatticeStabCanIsomMethod method,
                                     std::ostream &os) {
  LLLreduction<T, Tint> red = LLLreducedBasis<T, Tint>(GramMat, os);
  MyMatrix<Tint> B_red = LatticeCanonicalBasisReduced<T, Tint, Tgroup>(
      red.GramMatRed, method, os);
  // B_red GramRed B_red^T is canonical and GramRed = Pmat GramMat Pmat^T.
  return B_red * red.Pmat;
}

/*
  An isometry from the forms ListMat1 to ListMat2 by the method, if there is
  one. ListGenAut2, generators of the automorphism group of ListMat2 (possibly
  empty), prunes the Plesken-Souvignier search and is ignored otherwise.
 */
template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
LatticeIsometryReduced(std::vector<MyMatrix<T>> const &ListMat1,
                       std::vector<MyMatrix<T>> const &ListMat2,
                       std::vector<MyMatrix<Tint>> const &ListGenAut2,
                       LatticeStabCanIsomMethod method, std::ostream &os) {
  std::string name = NameLatticeStabCanIsomMethod(method);
  switch (method) {
  case LatticeStabCanIsomMethod::Graph:
    return ArithmeticEquivalenceMultiple<T, Tint, Tgroup>(ListMat1, ListMat2,
                                                          os);
  case LatticeStabCanIsomMethod::PleskenSouvignier:
    return PleskenSouvignierLatticeIsometry<T, Tint>(ListMat1, ListMat2,
                                                     ListGenAut2, os);
  case LatticeStabCanIsomMethod::Vinberg:
  case LatticeStabCanIsomMethod::Weyl:
    CheckLatticeMethodIntegral(ListMat1, name);
    CheckLatticeMethodIntegral(ListMat2, name);
    return ComputeIsometryVinberg<T, Tint, Tgroup>(ListMat1[0], ListMat2[0],
                                                   os);
  case LatticeStabCanIsomMethod::WeylChamber:
    CheckLatticeMethodIntegral(ListMat1, name);
    CheckLatticeMethodIntegral(ListMat2, name);
    return ComputeIsometryWeylChamber<T, Tint, Tgroup>(ListMat1[0],
                                                       ListMat2[0], os);
  case LatticeStabCanIsomMethod::IteratedShortest:
    return IteratedShortestEquivalenceMultiple<T, Tint, Tgroup>(ListMat1,
                                                                ListMat2, os);
  }
  std::cerr << "LATTICE_METHODS: unreachable isometry method\n";
  throw TerminalException{1};
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
LatticeIsometry(std::vector<MyMatrix<T>> const &ListMat1,
                std::vector<MyMatrix<T>> const &ListMat2,
                std::vector<MyMatrix<Tint>> const &ListGenAut2,
                LatticeStabCanIsomMethod method, std::ostream &os) {
  if (ListMat1.empty() || ListMat1.size() != ListMat2.size() ||
      ListMat1[0].rows() != ListMat2[0].rows()) {
    return {};
  }
  LLLreduction<T, Tint> red1 = LLLreducedBasis<T, Tint>(ListMat1[0], os);
  LLLreduction<T, Tint> red2 = LLLreducedBasis<T, Tint>(ListMat2[0], os);
  MyMatrix<Tint> P2inv = Inverse(red2.Pmat);
  // g preserves the forms of ListMat2, so P2 g P2^{-1} preserves the reduced
  // ones.
  std::vector<MyMatrix<Tint>> ListGenAut2Red;
  for (auto &g : ListGenAut2) {
    ListGenAut2Red.push_back(red2.Pmat * g * P2inv);
  }
  std::optional<MyMatrix<Tint>> opt = LatticeIsometryReduced<T, Tint, Tgroup>(
      LatticeMethodsChangeBasis(ListMat1, red1.Pmat),
      LatticeMethodsChangeBasis(ListMat2, red2.Pmat), ListGenAut2Red, method,
      os);
  if (!opt) {
    return {};
  }
  // P_red (P1 M1 P1^T) P_red^T = P2 M2 P2^T, so P = P2^{-1} P_red P1.
  return P2inv * (*opt) * red1.Pmat;
}

// The automorphism group of the forms ListMat by the method: generators and
// the order.
template <typename T, typename Tint, typename Tgroup>
IntAutomorphismGroup<Tint, Tgroup>
LatticeAutomorphismGroupReduced(std::vector<MyMatrix<T>> const &ListMat,
                                LatticeStabCanIsomMethod method,
                                std::ostream &os) {
  std::string name = NameLatticeStabCanIsomMethod(method);
  switch (method) {
  case LatticeStabCanIsomMethod::Graph:
    return ArithmeticAutomorphismGroupMultiple<T, Tint, Tgroup>(ListMat, os);
  case LatticeStabCanIsomMethod::PleskenSouvignier: {
    PleskenSouvignierAutomResult<Tint, Tgroup> result =
        PleskenSouvignierLatticeAutomorphism<T, Tint, Tgroup>(ListMat, os);
    return {std::move(result.ListGen), result.order};
  }
  case LatticeStabCanIsomMethod::Vinberg:
  case LatticeStabCanIsomMethod::Weyl: {
    CheckLatticeMethodIntegral(ListMat, name);
    VinbergAutom<Tint, Tgroup> vin =
        ComputeAutomorphismVinberg<T, Tint, Tgroup>(ListMat[0], os);
    return {std::move(vin.ListGen), vin.order};
  }
  case LatticeStabCanIsomMethod::WeylChamber: {
    CheckLatticeMethodIntegral(ListMat, name);
    VinbergAutom<Tint, Tgroup> vin =
        ComputeAutomorphismWeylChamber<T, Tint, Tgroup>(ListMat[0], os);
    return {std::move(vin.ListGen), vin.order};
  }
  case LatticeStabCanIsomMethod::IteratedShortest:
    return IteratedShortestAutomorphismGroupMultiple<T, Tint, Tgroup>(ListMat,
                                                                      os);
  }
  std::cerr << "LATTICE_METHODS: unreachable automorphism method\n";
  throw TerminalException{1};
}

template <typename T, typename Tint, typename Tgroup>
IntAutomorphismGroup<Tint, Tgroup>
LatticeAutomorphismGroup(std::vector<MyMatrix<T>> const &ListMat,
                         LatticeStabCanIsomMethod method, std::ostream &os) {
  LLLreduction<T, Tint> red = LLLreducedBasis<T, Tint>(ListMat[0], os);
  IntAutomorphismGroup<Tint, Tgroup> aut =
      LatticeAutomorphismGroupReduced<T, Tint, Tgroup>(
          LatticeMethodsChangeBasis(ListMat, red.Pmat), method, os);
  // g_red preserves Pmat M Pmat^T, so Pmat^{-1} g_red Pmat preserves M.
  MyMatrix<Tint> PmatInv = Inverse(red.Pmat);
  for (auto &g : aut.ListGen) {
    g = PmatInv * g * red.Pmat;
  }
  return aut;
}

// clang-format off
#endif  // SRC_LATT_LATTICEMETHODS_H_
// clang-format on

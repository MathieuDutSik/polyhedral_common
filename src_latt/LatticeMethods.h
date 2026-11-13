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
#include "SignatureSymmetric.h"
#include <optional>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_LATTICE_METHODS
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_LATTICE_METHODS
#endif

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

  heuristic is not a method of its own: it chooses one of the others for
  each input, from isometry invariants only (see ChooseHeuristicMethod), so
  that isometric lattices always get the same method and their canonical
  forms remain comparable.

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
  Heuristic
};

/*
  A method: its name, the algorithm and the vector family it stands for,
  what it computes, and a description for the usage and namelist texts.

  The algorithm and the family are separate choices. graph builds the graph
  of an invariant family of vectors weighted by the forms, and keeps the
  transformations that preserve Z^n (MatrixGroup.h), so that the family need
  only have full rank; plesken_souvignier searches the images of a basis in
  the family, and needs full rank too, with fewer rational leaves when it
  spans. Both take any kind of InvariantFamilyKind (graph:KIND,
  plesken_souvignier:KIND); the root-system methods (vinberg, weyl,
  weyl_chamber) build their own families. iterated_shortest is the name of
  graph:iterated_shortest.

  canonical_form: whether it gives a canonical form (the canonical form of
  plesken_souvignier, by orbit splitting, needs a family that spans).
  several_forms: whether it applies to a tuple of forms
  (A_1, ..., A_k), A_1 positive definite, whose isometries U A_i U^T = B_i
  and automorphisms are wanted together; the root-system methods read the
  first form only. A tuple whose forms are all multiples of the first is the
  problem of the first form alone, which every method solves.
 */
struct LatticeStabCanIsomMethodEntry {
  std::string name;
  LatticeStabCanIsomMethod method;
  std::optional<InvariantFamilyKind> family;
  bool canonical_form;
  bool several_forms;
  std::string description;
};

inline std::vector<LatticeStabCanIsomMethodEntry>
ListLatticeStabCanIsomMethod() {
  return {
      {"graph", LatticeStabCanIsomMethod::Graph, {}, true, true,
       "the graph of an invariant vector family, the characteristic one by "
       "default (graph:KIND chooses it, see InvariantFamilyKind)"},
      {"plesken_souvignier", LatticeStabCanIsomMethod::PleskenSouvignier, {},
       true, true,
       "the Plesken-Souvignier backtrack on an invariant vector family "
       "(plesken_souvignier:KIND chooses it, see InvariantFamilyKind); the "
       "canonical form by orbit splitting over its automorphism group, "
       "viable on large families (Leech), can decline on some inputs"},
      {"vinberg", LatticeStabCanIsomMethod::Vinberg, {}, true, false,
       "along the root decomposition, the automorphism group as W(R) rtimes "
       "Aut(L, rho); integral form"},
      {"weyl", LatticeStabCanIsomMethod::Weyl, {}, true, false,
       "as vinberg, with only the simple roots in the family; the isometry "
       "and the automorphism group are those of vinberg"},
      {"weyl_chamber", LatticeStabCanIsomMethod::WeylChamber, {}, true, false,
       "the characteristic vectors in the Weyl chamber and the simple roots, "
       "a family spanning the lattice; integral form"},
      {"iterated_shortest", LatticeStabCanIsomMethod::Graph,
       InvariantFamilyKind::IteratedShortest, true, true,
       "graph:iterated_shortest, the graph of the iterated shortest vectors "
       "up to full rank"},
      {"heuristic", LatticeStabCanIsomMethod::Heuristic, {}, true, true,
       "chooses from isometry invariants (ChooseHeuristicMethod): for one "
       "integral form weyl_chamber if it has vectors of norm at most 2, "
       "otherwise plesken_souvignier if it has many minimal vectors (Leech) "
       "and graph if not; for a tuple of forms plesken_souvignier on the "
       "shells_zbasis family"}};
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

inline LatticeStabCanIsomMethodEntry const &
LatticeStabCanIsomMethodGetEntry(LatticeStabCanIsomMethod method) {
  static std::vector<LatticeStabCanIsomMethodEntry> const list =
      ListLatticeStabCanIsomMethod();
  for (auto &entry : list) {
    if (entry.method == method) {
      return entry;
    }
  }
  std::cerr << "LATTICE_METHODS: a method without an entry\n";
  throw TerminalException{1};
}

// Whether the method applies to a tuple of forms (see the entries).
inline bool LatticeMethodAppliesToMultipleForms(LatticeStabCanIsomMethod method) {
  return LatticeStabCanIsomMethodGetEntry(method).several_forms;
}

// Whether the method gives a canonical form.
inline bool LatticeMethodGivesCanonicalForm(LatticeStabCanIsomMethod method) {
  return LatticeStabCanIsomMethodGetEntry(method).canonical_form;
}

/*
  A method with, for plesken_souvignier, the invariant vector family it runs
  on; written NAME or plesken_souvignier:KIND. Without a family the
  graph takes the characteristic family, plesken_souvignier
  DefaultInvariantFamilyKind, and the heuristic its own choice. The
  conversion from a method is implicit so that a method is a choice.
 */
struct LatticeMethodChoice {
  LatticeStabCanIsomMethod method;
  std::optional<InvariantFamilyKind> family;
  LatticeMethodChoice(LatticeStabCanIsomMethod method_)
      : method(method_), family() {}
  LatticeMethodChoice(LatticeStabCanIsomMethod method_,
                      InvariantFamilyKind family_)
      : method(method_), family(family_) {}
  InvariantFamilyKind family_for(size_t n_forms) const {
    if (family) {
      return *family;
    }
    if (method == LatticeStabCanIsomMethod::Graph) {
      return InvariantFamilyKind::Characteristic;
    }
    return DefaultInvariantFamilyKind(n_forms);
  }
};

inline LatticeMethodChoice ParseLatticeMethodChoice(std::string const &name) {
  size_t pos = name.find(':');
  if (pos == std::string::npos) {
    for (auto &entry : ListLatticeStabCanIsomMethod()) {
      if (entry.name == name) {
        return entry.family ? LatticeMethodChoice(entry.method, *entry.family)
                            : LatticeMethodChoice(entry.method);
      }
    }
    // the error message of the parse
    return LatticeMethodChoice(ParseLatticeStabCanIsomMethod(name));
  }
  LatticeStabCanIsomMethod method =
      ParseLatticeStabCanIsomMethod(name.substr(0, pos));
  if (method != LatticeStabCanIsomMethod::PleskenSouvignier &&
      method != LatticeStabCanIsomMethod::Graph) {
    std::cerr << "LATTICE_METHODS: only graph and plesken_souvignier take a "
              << "family, not " << name << "\n";
    throw TerminalException{1};
  }
  return LatticeMethodChoice(method,
                             ParseInvariantFamilyKind(name.substr(pos + 1)));
}

inline std::string NameLatticeMethodChoice(LatticeMethodChoice const &choice) {
  std::string name = NameLatticeStabCanIsomMethod(choice.method);
  if (choice.family) {
    name += ":" + NameInvariantFamilyKind(*choice.family);
  }
  return name;
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

// The root-system methods: a single integral form.
template <typename T>
void CheckLatticeMethodApplies(std::vector<MyMatrix<T>> const &ListMat,
                               LatticeStabCanIsomMethod method) {
  if (ListMat.size() != 1 && !LatticeMethodAppliesToMultipleForms(method)) {
    std::cerr << "LATTICE_METHODS: the method "
              << NameLatticeStabCanIsomMethod(method)
              << " applies to a single form, not to a tuple of "
              << ListMat.size() << "\n";
    throw TerminalException{1};
  }
  bool root_system = method == LatticeStabCanIsomMethod::Vinberg ||
                     method == LatticeStabCanIsomMethod::Weyl ||
                     method == LatticeStabCanIsomMethod::WeylChamber;
  if (root_system && !IsIntegralMatrix(ListMat[0])) {
    std::cerr << "LATTICE_METHODS: the method "
              << NameLatticeStabCanIsomMethod(method)
              << " needs an integral form\n";
    throw TerminalException{1};
  }
}

/*
  The ratios c_i with A_i = c_i A_1 when every form of the tuple is a
  multiple of the first, the tuple being then the problem of A_1 alone:
  a basis change preserving or mapping A_1 does the same to the c_i A_1.
 */
template <typename T>
std::optional<std::vector<T>>
LatticeMethodsProportionalRatios(std::vector<MyMatrix<T>> const &ListMat) {
  MyMatrix<T> const &A = ListMat[0];
  int n = A.rows();
  int i0 = 0, j0 = 0;
  // A is positive definite, so its diagonal is nonzero
  std::vector<T> ratios{T(1)};
  for (size_t k = 1; k < ListMat.size(); k++) {
    T c = ListMat[k](i0, j0) / A(i0, j0);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (ListMat[k](i, j) != c * A(i, j)) {
          return {};
        }
      }
    }
    ratios.push_back(c);
  }
  return ratios;
}

#ifdef SANITY_CHECK_LATTICE_METHODS
// The dispatchers below work on positive definite forms: the LLL reduction
// and the short vector enumerations of the methods do not terminate on an
// indefinite or negative definite form, so a wrong input hangs instead of
// failing.
template <typename T>
void CheckLatticeMethodPositiveDefinite(MyMatrix<T> const &GramMat,
                                        std::string const &context,
                                        std::ostream &os) {
  if (!IsSymmetricMatrix(GramMat) || !IsPositiveDefinite(GramMat, os)) {
    std::cerr << "LATTICE_METHODS: the input of " << context
              << " is not a symmetric positive definite matrix\n";
    throw TerminalException{1};
  }
}
#endif

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

// The family of the kind of the reduced Gram matrix GramRed, with whether
// it spans the lattice, for the graph methods.
template <typename T, typename Tint>
CanonicVectorFamily<Tint> LatticeMethodsFamily(MyMatrix<T> const &GramRed,
                                               InvariantFamilyKind kind,
                                               std::ostream &os) {
  return get_canonic_vector_family<Tint>(
      InvariantFamilyHalf<T, Tint>(GramRed, kind, os));
}

/*
  The method chosen by heuristic for the LLL-reduced forms ListMat. The
  choice depends only on isometry invariants of the forms (whether the
  lattice has vectors of norm at most 2, its number of minimal vectors, the
  number of forms), never on the basis, so isometric inputs get the same
  method; this is what keeps the canonical forms of heuristic comparable
  with each other. A tuple of proportional forms has been reduced to its
  first form by the caller.

  * A tuple of forms: plesken_souvignier on the shells_zbasis family, the
    canonical form by orbit splitting with the graph as fallback (the
    splitting declines for whole isometry classes). The smaller
    shells_fullrank family (the roots of a Niemeier lattice) does not span
    when the roots do not, and the search then meets many rational leaves:
    one centraliser of (E8^3, i) took minutes on it. The graph of the
    characteristic family of the tuple is slow when the common automorphism
    group is large (10^12 for W(A3) on E8+E8), and the Z-span closure of the
    other families has millions of vectors on some rank 24 lattices.
  * A non-integral form: graph for the canonical form, plesken_souvignier
    for the isometry and the group.
  * An integral lattice with vectors of norm at most 2 (roots): weyl_chamber,
    which removes the Weyl group. A diagonal entry at most 2 of the reduced
    Gram matrix shows it without any enumeration.
  * Otherwise the lattice has no roots, weyl_chamber would fall back on the
    graph of its characteristic family, and that graph is quadratic in the
    family: with more than heuristic_max_pairs pairs of minimal vectors
    (98280 for the Leech lattice) plesken_souvignier, which never builds it,
    and graph below.
 */
inline constexpr size_t heuristic_max_pairs = 5000;

template <typename T, typename Tint>
LatticeMethodChoice
ChooseHeuristicMethod(std::vector<MyMatrix<T>> const &ListMat, bool canonical,
                      std::ostream &os) {
  if (ListMat.size() != 1) {
    return LatticeMethodChoice(LatticeStabCanIsomMethod::PleskenSouvignier,
                               InvariantFamilyKind::ShellsZbasis);
  }
  if (!IsIntegralMatrix(ListMat[0])) {
    return canonical ? LatticeStabCanIsomMethod::Graph
                     : LatticeStabCanIsomMethod::PleskenSouvignier;
  }
  MyMatrix<T> const &G = ListMat[0];
  for (int i = 0; i < G.rows(); i++) {
    if (G(i, i) <= T(2)) {
      return LatticeStabCanIsomMethod::WeylChamber;
    }
  }
  Tshortest<T, Tint> rec = T_ShortestVectorHalf<T, Tint>(G, os);
  if (rec.min <= T(2)) {
    return LatticeStabCanIsomMethod::WeylChamber;
  }
  size_t n_pairs = static_cast<size_t>(rec.SHV.rows());
#ifdef DEBUG_LATTICE_METHODS
  os << "LATTICE_METHODS: heuristic, no roots, minimum " << rec.min << ", "
     << n_pairs << " pairs of minimal vectors\n";
#endif
  return n_pairs > heuristic_max_pairs
             ? LatticeStabCanIsomMethod::PleskenSouvignier
             : LatticeStabCanIsomMethod::Graph;
}

/*
  The canonical form of the tuple of forms ListMat by the method: a basis B,
  as rows, with the B A_i B^T canonical. The Reduced variant assumes ListMat
  LLL-reduced (its first form) and not proportional.
 */
template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint>
LatticeCanonicalBasisMultipleReduced(std::vector<MyMatrix<T>> const &ListMat,
                                     LatticeMethodChoice choice,
                                     std::ostream &os) {
  bool heuristic = choice.method == LatticeStabCanIsomMethod::Heuristic;
  if (heuristic) {
    choice = ChooseHeuristicMethod<T, Tint>(ListMat, true, os);
  }
  LatticeStabCanIsomMethod method = choice.method;
  if (!LatticeMethodGivesCanonicalForm(method)) {
    LatticeStabCanIsomMethodUnavailable(method, "a canonical form");
  }
  CheckLatticeMethodApplies(ListMat, method);
  switch (method) {
  case LatticeStabCanIsomMethod::Graph: {
    InvariantFamilyKind kind = choice.family_for(ListMat.size());
    if (kind == InvariantFamilyKind::DiagonalBound) {
      std::cerr << "LATTICE_METHODS: a canonical form needs an intrinsic "
                   "family, the diagonal_bound one depends on the basis\n";
      throw TerminalException{1};
    }
    if (kind == InvariantFamilyKind::Characteristic) {
      if (ListMat.size() == 1) {
        return ComputeCanonicalForm<T, Tint, Tgroup>(ListMat[0], os);
      }
      return ComputeCanonicalFormMultiple<T, Tint, Tgroup>(ListMat, os);
    }
    return ComputeCanonicalForm_family<T, Tint, Tgroup>(
        ListMat, LatticeMethodsFamily<T, Tint>(ListMat[0], kind, os), os);
  }
  case LatticeStabCanIsomMethod::PleskenSouvignier: {
    // The failure is a property of the isometry class, so the forms the
    // method does canonicalize remain comparable with each other.
    std::optional<MyMatrix<Tint>> opt =
        OrbitSplitCanonicalizeListMat<T, Tint, Tgroup>(
            ListMat, choice.family_for(ListMat.size()), os);
    if (!opt && heuristic) {
      // The same failure for every input of the isometry class, so the
      // fallback keeps the canonical forms of heuristic comparable.
      if (ListMat.size() == 1) {
        return ComputeCanonicalForm<T, Tint, Tgroup>(ListMat[0], os);
      }
      return ComputeCanonicalFormMultiple<T, Tint, Tgroup>(ListMat, os);
    }
    if (!opt) {
      std::cerr << "LATTICE_METHODS: the orbit splitting could not select a "
                   "canonical prefix, no canonical form for this input\n";
      throw TerminalException{1};
    }
    return *opt;
  }
  case LatticeStabCanIsomMethod::Vinberg:
    return ComputeCanonicalFormRootDecomposed<T, Tint, Tgroup>(ListMat[0], os);
  case LatticeStabCanIsomMethod::Weyl:
    return ComputeCanonicalFormWeylReduced<T, Tint, Tgroup>(ListMat[0], os);
  case LatticeStabCanIsomMethod::WeylChamber:
    return ComputeCanonicalFormWeylChamber<T, Tint, Tgroup>(ListMat[0], os);
  case LatticeStabCanIsomMethod::Heuristic:
    break;
  }
  LatticeStabCanIsomMethodUnavailable(method, "a canonical form");
}

template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint> LatticeCanonicalBasisReduced(MyMatrix<T> const &GramMat,
                                            LatticeMethodChoice choice,
                                            std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat{GramMat};
  return LatticeCanonicalBasisMultipleReduced<T, Tint, Tgroup>(ListMat, choice,
                                                               os);
}

template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint>
LatticeCanonicalBasisMultiple(std::vector<MyMatrix<T>> const &ListMat,
                              LatticeMethodChoice choice, std::ostream &os) {
#ifdef SANITY_CHECK_LATTICE_METHODS
  CheckLatticeMethodPositiveDefinite(ListMat[0], "LatticeCanonicalBasis", os);
#endif
  LLLreduction<T, Tint> red = LLLreducedBasis<T, Tint>(ListMat[0], os);
  std::vector<MyMatrix<T>> ListMatRed =
      LatticeMethodsChangeBasis(ListMat, red.Pmat);
  if (ListMat.size() > 1 && LatticeMethodsProportionalRatios(ListMat)) {
    ListMatRed = {ListMatRed[0]};
  }
  MyMatrix<Tint> B_red = LatticeCanonicalBasisMultipleReduced<T, Tint, Tgroup>(
      ListMatRed, choice, os);
  // B_red GramRed B_red^T is canonical and GramRed = Pmat GramMat Pmat^T.
  return B_red * red.Pmat;
}

template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint> LatticeCanonicalBasis(MyMatrix<T> const &GramMat,
                                     LatticeMethodChoice choice,
                                     std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat{GramMat};
  return LatticeCanonicalBasisMultiple<T, Tint, Tgroup>(ListMat, choice, os);
}

/*
  An isometry from the forms ListMat1 to ListMat2 by the method, if there is
  one: P with P A_i P^T = B_i. ListGenAut2, generators of the automorphism
  group of ListMat2 (possibly empty), prunes the Plesken-Souvignier search
  and is ignored otherwise.
 */
template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
LatticeIsometryReduced(std::vector<MyMatrix<T>> const &ListMat1,
                       std::vector<MyMatrix<T>> const &ListMat2,
                       std::vector<MyMatrix<Tint>> const &ListGenAut2,
                       LatticeMethodChoice choice, std::ostream &os) {
  // Lattices with different invariants are not isometric, and the method
  // chosen for the first one answers that correctly.
  if (choice.method == LatticeStabCanIsomMethod::Heuristic) {
    choice = ChooseHeuristicMethod<T, Tint>(ListMat1, false, os);
  }
  LatticeStabCanIsomMethod method = choice.method;
  CheckLatticeMethodApplies(ListMat1, method);
  CheckLatticeMethodApplies(ListMat2, method);
  switch (method) {
  case LatticeStabCanIsomMethod::Graph: {
    InvariantFamilyKind kind = choice.family_for(ListMat1.size());
    if (kind == InvariantFamilyKind::Characteristic) {
      return ArithmeticEquivalenceMultiple<T, Tint, Tgroup>(ListMat1, ListMat2,
                                                            os);
    }
    if (kind == InvariantFamilyKind::DiagonalBound) {
      std::cerr << "LATTICE_METHODS: the diagonal_bound family depends on the "
                   "basis, it cannot be used on two lattices\n";
      throw TerminalException{1};
    }
    return ArithmeticEquivalenceMultiple_family<T, Tint, Tgroup>(
        ListMat1, LatticeMethodsFamily<T, Tint>(ListMat1[0], kind, os),
        ListMat2, LatticeMethodsFamily<T, Tint>(ListMat2[0], kind, os), os);
  }
  case LatticeStabCanIsomMethod::PleskenSouvignier:
    return PleskenSouvignierLatticeIsometry<T, Tint>(
        ListMat1, ListMat2, ListGenAut2, choice.family_for(ListMat1.size()),
        os);
  case LatticeStabCanIsomMethod::Vinberg:
  case LatticeStabCanIsomMethod::Weyl:
    return ComputeIsometryVinberg<T, Tint, Tgroup>(ListMat1[0], ListMat2[0],
                                                   os);
  case LatticeStabCanIsomMethod::WeylChamber:
    return ComputeIsometryWeylChamber<T, Tint, Tgroup>(ListMat1[0],
                                                       ListMat2[0], os);
  case LatticeStabCanIsomMethod::Heuristic:
    break;
  }
  std::cerr << "LATTICE_METHODS: unreachable isometry method\n";
  throw TerminalException{1};
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
LatticeIsometry(std::vector<MyMatrix<T>> const &ListMat1,
                std::vector<MyMatrix<T>> const &ListMat2,
                std::vector<MyMatrix<Tint>> const &ListGenAut2,
                LatticeMethodChoice choice, std::ostream &os) {
  if (ListMat1.empty() || ListMat1.size() != ListMat2.size() ||
      ListMat1[0].rows() != ListMat2[0].rows()) {
    return {};
  }
  // The integrality of each form is preserved by the isometries, and the
  // method, chosen from the first tuple, may refuse a non-integral form of
  // the second (the root-system methods): a pair that differs in it is
  // answered here.
  for (size_t i = 0; i < ListMat1.size(); i++) {
    if (IsIntegralMatrix(ListMat1[i]) != IsIntegralMatrix(ListMat2[i])) {
      return {};
    }
  }
#ifdef SANITY_CHECK_LATTICE_METHODS
  CheckLatticeMethodPositiveDefinite(ListMat1[0], "LatticeIsometry", os);
  CheckLatticeMethodPositiveDefinite(ListMat2[0], "LatticeIsometry", os);
#endif
  LLLreduction<T, Tint> red1 = LLLreducedBasis<T, Tint>(ListMat1[0], os);
  LLLreduction<T, Tint> red2 = LLLreducedBasis<T, Tint>(ListMat2[0], os);
  MyMatrix<Tint> P2inv = Inverse(red2.Pmat);
  // g preserves the forms of ListMat2, so P2 g P2^{-1} preserves the reduced
  // ones.
  std::vector<MyMatrix<Tint>> ListGenAut2Red;
  for (auto &g : ListGenAut2) {
    ListGenAut2Red.push_back(red2.Pmat * g * P2inv);
  }
  std::vector<MyMatrix<T>> ListMatRed1 =
      LatticeMethodsChangeBasis(ListMat1, red1.Pmat);
  std::vector<MyMatrix<T>> ListMatRed2 =
      LatticeMethodsChangeBasis(ListMat2, red2.Pmat);
  if (ListMat1.size() > 1) {
    // proportional tuples: the same ratios, then the first forms alone
    std::optional<std::vector<T>> r1 = LatticeMethodsProportionalRatios(ListMat1);
    std::optional<std::vector<T>> r2 = LatticeMethodsProportionalRatios(ListMat2);
    if (r1.has_value() != r2.has_value()) {
      return {};
    }
    if (r1) {
      if (*r1 != *r2) {
        return {};
      }
      ListMatRed1 = {ListMatRed1[0]};
      ListMatRed2 = {ListMatRed2[0]};
    }
  }
  std::optional<MyMatrix<Tint>> opt = LatticeIsometryReduced<T, Tint, Tgroup>(
      ListMatRed1, ListMatRed2, ListGenAut2Red, choice, os);
  if (!opt) {
    return {};
  }
  // P_red (P1 M1 P1^T) P_red^T = P2 M2 P2^T, so P = P2^{-1} P_red P1.
  MyMatrix<Tint> P = P2inv * (*opt) * red1.Pmat;
#ifdef SANITY_CHECK_LATTICE_METHODS
  MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(P);
  for (size_t i = 0; i < ListMat1.size(); i++) {
    if (P_T * ListMat1[i] * P_T.transpose() != ListMat2[i]) {
      std::cerr << "LATTICE_METHODS: the isometry of "
                << NameLatticeMethodChoice(choice)
                << " does not satisfy P A P^T = B\n";
      throw TerminalException{1};
    }
  }
#endif
  return P;
}

// The automorphism group of the forms ListMat by the method: generators g
// with g A g^T = A for every form, and the order.
template <typename T, typename Tint, typename Tgroup>
IntAutomorphismGroup<Tint, Tgroup>
LatticeAutomorphismGroupReduced(std::vector<MyMatrix<T>> const &ListMat,
                                LatticeMethodChoice choice,
                                std::ostream &os) {
  if (choice.method == LatticeStabCanIsomMethod::Heuristic) {
    choice = ChooseHeuristicMethod<T, Tint>(ListMat, false, os);
  }
  LatticeStabCanIsomMethod method = choice.method;
  CheckLatticeMethodApplies(ListMat, method);
  switch (method) {
  case LatticeStabCanIsomMethod::Graph: {
    InvariantFamilyKind kind = choice.family_for(ListMat.size());
    if (kind == InvariantFamilyKind::Characteristic) {
      return ArithmeticAutomorphismGroupMultiple<T, Tint, Tgroup>(ListMat, os);
    }
    CanonicVectorFamily<Tint> fam =
        LatticeMethodsFamily<T, Tint>(ListMat[0], kind, os);
    return GetIntAutomorphism_Family<T, Tint, Tgroup>(ListMat, fam,
                                                      fam.get_full(), os);
  }
  case LatticeStabCanIsomMethod::PleskenSouvignier: {
    PleskenSouvignierAutomResult<Tint, Tgroup> result =
        PleskenSouvignierLatticeAutomorphism<T, Tint, Tgroup>(
            ListMat, choice.family_for(ListMat.size()), os);
    return {std::move(result.ListGen), result.order};
  }
  case LatticeStabCanIsomMethod::Vinberg:
  case LatticeStabCanIsomMethod::Weyl: {
    VinbergAutom<Tint, Tgroup> vin =
        ComputeAutomorphismVinberg<T, Tint, Tgroup>(ListMat[0], os);
    return {std::move(vin.ListGen), vin.order};
  }
  case LatticeStabCanIsomMethod::WeylChamber: {
    VinbergAutom<Tint, Tgroup> vin =
        ComputeAutomorphismWeylChamber<T, Tint, Tgroup>(ListMat[0], os);
    return {std::move(vin.ListGen), vin.order};
  }
  case LatticeStabCanIsomMethod::Heuristic:
    break;
  }
  std::cerr << "LATTICE_METHODS: unreachable automorphism method\n";
  throw TerminalException{1};
}

template <typename T, typename Tint, typename Tgroup>
IntAutomorphismGroup<Tint, Tgroup>
LatticeAutomorphismGroup(std::vector<MyMatrix<T>> const &ListMat,
                         LatticeMethodChoice choice, std::ostream &os) {
#ifdef SANITY_CHECK_LATTICE_METHODS
  CheckLatticeMethodPositiveDefinite(ListMat[0], "LatticeAutomorphismGroup",
                                     os);
#endif
  LLLreduction<T, Tint> red = LLLreducedBasis<T, Tint>(ListMat[0], os);
  std::vector<MyMatrix<T>> ListMatRed =
      LatticeMethodsChangeBasis(ListMat, red.Pmat);
  if (ListMat.size() > 1 && LatticeMethodsProportionalRatios(ListMat)) {
    ListMatRed = {ListMatRed[0]};
  }
  IntAutomorphismGroup<Tint, Tgroup> aut =
      LatticeAutomorphismGroupReduced<T, Tint, Tgroup>(ListMatRed, choice, os);
  // g_red preserves Pmat M Pmat^T, so Pmat^{-1} g_red Pmat preserves M.
  MyMatrix<Tint> PmatInv = Inverse(red.Pmat);
  for (auto &g : aut.ListGen) {
    g = PmatInv * g * red.Pmat;
  }
#ifdef SANITY_CHECK_LATTICE_METHODS
  for (auto &g : aut.ListGen) {
    MyMatrix<T> g_T = UniversalMatrixConversion<T, Tint>(g);
    for (auto &M : ListMat) {
      if (g_T * M * g_T.transpose() != M) {
        std::cerr << "LATTICE_METHODS: an automorphism of "
                  << NameLatticeMethodChoice(choice)
                  << " does not satisfy g A g^T = A\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  return aut;
}

// clang-format off
#endif  // SRC_LATT_LATTICEMETHODS_H_
// clang-format on

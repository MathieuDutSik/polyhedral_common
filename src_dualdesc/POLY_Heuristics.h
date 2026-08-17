// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_HEURISTICS_H_
#define SRC_POLY_POLY_HEURISTICS_H_

// clang-format off
#include "Basic_file.h"
#include "Heuristic_ThompsonSampling.h"
#include "Namelist.h"
#include <array>
#include <format>
#include <limits>
#include <map>
#include <string>
#include <string_view>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_HEURISTICS
#endif

//
// Typed dual-description heuristic input and heuristic
//
// The dual-description heuristics no longer read a std::map<std::string,
// TintGroup>. Each heuristic is typed over a set of named fields (a field enum)
// and evaluated against an input struct exposing get_field(field). The generic
// machinery is shared by the polytope-info heuristics and the orbit-splitting
// heuristic.

// --- The fields, one enum per heuristic input kind ---
enum class PolytopeField { groupsize, incidence, rank, delta, level, time };
enum class OrbitSplitField { groupsize_big, groupsize_sma, index, n_orbit };

// The full list of polytope fields, the single source of truth used both by
// the heuristic evaluation and by the Thompson-sampler named-map view below.
inline constexpr std::array<PolytopeField, 6> polytope_all_fields = {
    PolytopeField::groupsize, PolytopeField::incidence, PolytopeField::rank,
    PolytopeField::delta,     PolytopeField::level,     PolytopeField::time};

// Parse a field name into the field enum (tag-dispatched on the enum type).
inline PolytopeField parse_heuristic_field(std::string const &name,
                                           PolytopeField *) {
  if (name == "groupsize")
    return PolytopeField::groupsize;
  if (name == "incidence")
    return PolytopeField::incidence;
  if (name == "rank")
    return PolytopeField::rank;
  if (name == "delta")
    return PolytopeField::delta;
  if (name == "level")
    return PolytopeField::level;
  if (name == "time")
    return PolytopeField::time;
  std::cerr << "HEU: unknown polytope field name=" << name << "\n";
  throw TerminalException{1};
}

inline OrbitSplitField parse_heuristic_field(std::string const &name,
                                             OrbitSplitField *) {
  if (name == "groupsize_big")
    return OrbitSplitField::groupsize_big;
  if (name == "groupsize_sma")
    return OrbitSplitField::groupsize_sma;
  if (name == "index")
    return OrbitSplitField::index;
  if (name == "n_orbit")
    return OrbitSplitField::n_orbit;
  std::cerr << "HEU: unknown orbit split field name=" << name << "\n";
  throw TerminalException{1};
}

template <>
struct std::formatter<PolytopeField> : std::formatter<std::string_view> {
  static constexpr std::string_view name(PolytopeField f) {
    switch (f) {
    case PolytopeField::groupsize:
      return "groupsize";
    case PolytopeField::incidence:
      return "incidence";
    case PolytopeField::rank:
      return "rank";
    case PolytopeField::delta:
      return "delta";
    case PolytopeField::level:
      return "level";
    case PolytopeField::time:
      return "time";
    }
    return "unknown";
  }
  template <typename FormatContext>
  auto format(PolytopeField f, FormatContext &ctx) const {
    return std::formatter<std::string_view>::format(name(f), ctx);
  }
};

template <>
struct std::formatter<OrbitSplitField> : std::formatter<std::string_view> {
  static constexpr std::string_view name(OrbitSplitField f) {
    switch (f) {
    case OrbitSplitField::groupsize_big:
      return "groupsize_big";
    case OrbitSplitField::groupsize_sma:
      return "groupsize_sma";
    case OrbitSplitField::index:
      return "index";
    case OrbitSplitField::n_orbit:
      return "n_orbit";
    }
    return "unknown";
  }
  template <typename FormatContext>
  auto format(OrbitSplitField f, FormatContext &ctx) const {
    return std::formatter<std::string_view>::format(name(f), ctx);
  }
};

// --- The input structs, one per heuristic input kind ---
// The group sizes stay full precision, the combinatorial quantities are plain
// ints / size_t and are promoted to TintGroup only for comparison.
template <typename TintGroup> struct PolytopeInputInfo {
  TintGroup groupsize;
  int incidence;
  int rank;
  int delta;
  int level;
  uint64_t time;
  TintGroup get_field(PolytopeField f) const {
    switch (f) {
    case PolytopeField::groupsize:
      return groupsize;
    case PolytopeField::incidence:
      return UniversalScalarConversion<TintGroup, int>(incidence);
    case PolytopeField::rank:
      return UniversalScalarConversion<TintGroup, int>(rank);
    case PolytopeField::delta:
      return UniversalScalarConversion<TintGroup, int>(delta);
    case PolytopeField::level:
      return UniversalScalarConversion<TintGroup, int>(level);
    case PolytopeField::time:
      return UniversalScalarConversion<TintGroup, uint64_t>(time);
    }
    return TintGroup(0);
  }
  // Named-map view of the typed info, consumed by the map-based heuristic
  // layers: the shared, disk-persisted Thompson learner (intrinsically keyed by
  // field name) and the generic HeuristicEvaluation. The names come from
  // to_string, so there is a single source of truth and no
  // hand-written string literals.
  std::map<std::string, TintGroup> named_map() const {
    std::map<std::string, TintGroup> m;
    for (PolytopeField f : polytope_all_fields)
      m[std::format("{}", f)] = get_field(f);
    return m;
  }
};

template <typename TintGroup> struct OrbitSplitInputInfo {
  TintGroup groupsize_big;
  TintGroup groupsize_sma;
  TintGroup index;
  size_t n_orbit;
  TintGroup get_field(OrbitSplitField f) const {
    switch (f) {
    case OrbitSplitField::groupsize_big:
      return groupsize_big;
    case OrbitSplitField::groupsize_sma:
      return groupsize_sma;
    case OrbitSplitField::index:
      return index;
    case OrbitSplitField::n_orbit:
      return UniversalScalarConversion<TintGroup, size_t>(n_orbit);
    }
    return TintGroup(0);
  }
};

// --- The generic typed heuristic ---
// Same shape as the generic TheHeuristic, but conditions reference a typed
// field. It is obtained once, at load time, from the string-parsed generic
// heuristic and afterwards evaluated with no map.
template <typename TintGroup, typename TField> struct TypedHeuristicCondition {
  TField field;
  HeuristicOp op;
  TintGroup value;
};

template <typename TintGroup, typename TField>
struct TypedHeuristicFullCondition {
  std::vector<TypedHeuristicCondition<TintGroup, TField>> conditions;
  std::string result;
};

template <typename TintGroup, typename TField> struct TypedHeuristic {
  std::vector<TypedHeuristicFullCondition<TintGroup, TField>> tests;
  std::string default_result;
};

template <typename TintGroup>
using DualDescHeuristic = TypedHeuristic<TintGroup, PolytopeField>;
template <typename TintGroup>
using OrbitSplitHeuristic = TypedHeuristic<TintGroup, OrbitSplitField>;

template <typename TintGroup, typename TField>
TypedHeuristic<TintGroup, TField>
convert_typed_heuristic(TheHeuristic<TintGroup> const &heu) {
  TypedHeuristic<TintGroup, TField> result;
  for (auto const &eFullCond : heu.AllTests) {
    TypedHeuristicFullCondition<TintGroup, TField> new_full;
    new_full.result = eFullCond.TheResult;
    for (auto const &eSingCond : eFullCond.TheConditions)
      new_full.conditions.push_back(
          {parse_heuristic_field(eSingCond.eCond,
                                 static_cast<TField *>(nullptr)),
           eSingCond.eType, eSingCond.NumValue});
    result.tests.push_back(std::move(new_full));
  }
  result.default_result = heu.DefaultResult;
  return result;
}

template <typename TintGroup, typename TField, typename TInput>
std::string typed_heuristic_evaluation(TypedHeuristic<TintGroup, TField> const &heu,
                                       TInput const &info) {
  for (auto const &eFullCond : heu.tests) {
    bool IsOK = true;
    for (auto const &eSingCond : eFullCond.conditions) {
      TintGroup eValue = info.get_field(eSingCond.field);
      bool WeMatch = false;
      switch (eSingCond.op) {
      case HeuristicOp::Gt:
        WeMatch = (eValue > eSingCond.value);
        break;
      case HeuristicOp::Ge:
        WeMatch = (eValue >= eSingCond.value);
        break;
      case HeuristicOp::Eq:
        WeMatch = (eValue == eSingCond.value);
        break;
      case HeuristicOp::Lt:
        WeMatch = (eValue < eSingCond.value);
        break;
      case HeuristicOp::Le:
        WeMatch = (eValue <= eSingCond.value);
        break;
      }
      if (!WeMatch) {
        IsOK = false;
        break;
      }
    }
    if (IsOK)
      return eFullCond.result;
  }
  return heu.default_result;
}

template <typename TintGroup, typename TField>
std::ostream &operator<<(std::ostream &os,
                         TypedHeuristic<TintGroup, TField> const &heu) {
  size_t len = heu.tests.size();
  for (size_t i = 0; i < len; i++) {
    TypedHeuristicFullCondition<TintGroup, TField> const &eFullCond =
        heu.tests[i];
    os << "   i=" << i << "/" << len;
    for (auto const &eSingCond : eFullCond.conditions)
      os << " (" << std::format("{}", eSingCond.field) << " "
         << static_cast<int>(eSingCond.op) << " " << eSingCond.value << ")";
    os << " => " << eFullCond.result << "\n";
  }
  os << "      Default=" << heu.default_result;
  return os;
}

// Thin wrappers preserving the polytope-info heuristic call sites.
template <typename TintGroup>
DualDescHeuristic<TintGroup>
convert_dual_desc_heuristic(TheHeuristic<TintGroup> const &heu) {
  return convert_typed_heuristic<TintGroup, PolytopeField>(heu);
}

template <typename TintGroup>
std::string
dual_desc_heuristic_evaluation(DualDescHeuristic<TintGroup> const &heu,
                               PolytopeInputInfo<TintGroup> const &info) {
  return typed_heuristic_evaluation(heu, info);
}

//
// Heuristic business
//

// Description: Heuristic for whether to apply the adjacency
// decomposition method or not. Following points have to be taken
// into account:
// * If the symmetry group is very large then one can expect benefits
// from using it.
// * Applying the method recursively makes everything slower.
//
// Variable available for making the decision:
// * "incidence": High incidence means that classic algorithm will be slower.
// * "groupsize": The size of the group available for the computation
// * "level": The depth of the calls to the recursive adjacency
//    decomposition.
// * "rank": the dimension of the polytope in question
// * "delta": the difference between the incidence and the rank.
//
// Possible output:
// * "nosplit": Call directly the classical adjacency decomposition.
// * "split": Call the Adjacency decomposition method on the polytope
//    or face considered.
// The first rule guards against the recursive decomposition being entered on
// polytopes that have almost as few vertices as their rank: there the direct
// dual description costs microseconds while one level of the decomposition
// costs milliseconds (its database of orbits canonicalizes every face). It
// comes first because the group of such a polytope is often large, and the
// group-size rule below would otherwise send a simplex into the recursion.
template <typename T> TheHeuristic<T> StandardHeuristicSplitting() {
  std::vector<std::string> ListString = {"5",
                                         "1 delta < 20 nosplit",
                                         "2 groupsize > 5000 level < 3 split",
                                         "1 incidence < 40 nosplit",
                                         "1 level > 1 nosplit",
                                         "1 groupsize < 100 nosplit",
                                         "split"};
  return HeuristicFromListString<T>(ListString);
}

// Description: Whether to use a communication thread or not for the
// MPI calls. Having this thread avoids the accumulation of data
// in the buffer when a computation is done deeper.
//
// Variable available for making the decision:
// * "incidence": High incidence means that classic algorithm will be slower.
// * "groupsize": The size of the group available for the computation
// * "level": The depth of the calls to the recursive adjacency
//    decomposition.
// * "rank": the dimension of the polytope in question
// * "delta": the difference between the incidence and the rank.
//
// Possible output:
// * "no": Do not create the communication thread.
// * "yes": Create the communication thread.
template <typename T> TheHeuristic<T> StandardHeuristicCommThread() {
  std::vector<std::string> ListString = {"0", "no"};
  return HeuristicFromListString<T>(ListString);
}

// Description: When computing with the adjacency decomposition
// method, we may want ddditional symmetries. This require some
// call to partition backtrack algorithms. Things to consider:
// * The computation of additional symmetries can make the
// computation much faster. If the number of orbits is divided
// by 2 then the computation should be divided by 2.
// * There is a cost to computing the additional symmetries.
// which may yield nothing if the full symmetry group is equal
// to the stabilizer.
// * If we have additional symmetries, then we need to split the
// orbits for the big group to the orbits for the small group.
// This has costs and require algorithmic choices.
//
// Variable available for making the decision:
// * "incidence": High incidence means that classic algorithm will be slower.
// * "groupsize": The size of the group available for the computation
// * "level": The depth of the calls to the recursive adjacency
//    decomposition.
// * "rank": the dimension of the polytope in question
// * "delta": the difference between the incidence and the rank.
//
// Possible output:
// * "no": Do not compute additional symmetries
// * "yes": Compute additional symmetries
template <typename T> TheHeuristic<T> StandardHeuristicAdditionalSymmetry() {
  std::vector<std::string> ListString = {"1", "1 incidence < 50 no", "yes"};
  return HeuristicFromListString<T>(ListString);
}

// Description: When computing with the adjacency decomposition
// method, we compute the dual description of faces. Those dual
// descriptions are expensive and if by any chance they would occur
// again then it would be nice to use what has been computed.
// This is what a banking system is.
// Things to consider:
// * Storing everything that has been computed is clearly
// a bad idea as the banking system gets flooded and the disk as
// wel if saving is selected.
// * Not storing at all negates the usefulness of the method.
//
// Variable available for making the decision:
// * "incidence": High incidence means that classic algorithm will be slower.
// * "groupsize": The size of the group available for the computation
// * "level": The depth of the calls to the recursive adjacency
//    decomposition.
// * "rank": the dimension of the polytope in question
// * "delta": the difference between the incidence and the rank.
// * "time" the tuntime in second. Saving what took long time is a good
//   idea in general.
//
// Possible output:
// * "no": Do not compute additional symmetries
// * "yes": Compute additional symmetries
template <typename T> TheHeuristic<T> StandardHeuristicBankSave() {
  std::vector<std::string> ListString = {"2", "1 time < 300 no",
                                         "1 incidence < 50 no", "yes"};
  return HeuristicFromListString<T>(ListString);
}

// Description: This is the most important method choice. There are
// many algorithm to consider for computing the dual description and
// you have to choose the one that suits you best.
// Things to consider:
// * Algorithm like lrs are working for polytopes with a lot of
// vertices. Not so much for polytopes with degenerate vertices.
// * Algorithm like CDD works the reverse.
// * Not all algorithm are available for all the types. For rational
// mpq_class everything works, but for algebraic types, that will not
// be true.
// * The normaliz external function has an overhead to calling external
// programs.
// * The lrs algorithm is only using addition, difference and multiplication.
// Hence it makes sense to reduce the fractions to integer and then do
// integer computations which are faster than rational computations. This
// is not possible for all fields.
//
// Variable available for making the decision:
// * "incidence": High incidence means that classic algorithm will be slower.
// * "groupsize": The size of the group available for the computation
// * "level": The depth of the calls to the recursive adjacency
//    decomposition.
// * "rank": the dimension of the polytope in question
// * "delta": the difference between the incidence and the rank.
//
// Possible output:
// * "cdd": The templatized cdd (only for fields)
// * "lrs": The templatized lrs (for fields and rings). When T is a field,
//   the implementation internally uses the fraction-reduction variant;
//   this is transparent to callers.
// * "normaliz": the native port of the Normaliz primal algorithm
//   (Fourier-Motzkin with pyramid decomposition, fields and rings)
// The defaults come from a measurement over the 7066 polytopes that one
// enumeration of iso-Delaunay domains sends to the dual description (ring
// arithmetic, rank 8 to 10, "delta" = incidence - rank from 0 to 138). Mean
// time per polytope, in microseconds:
//
//   delta      0     3     8    20    48    96   118
//   lrs        3     7    47   432  2480 84392 247911
//   normaliz   6    11    60    71   105   857   3506
//   cdd      194   210   343   434   651  3133   6643
//   bb        36    61   271   367   657  6201  30096
//
// So lrs is unbeatable while the polytope is close to a simplex (it is
// then almost pure pivoting) and collapses as soon as the polytope has
// many more vertices than its rank, where normaliz is one to two orders of
// magnitude ahead. cdd never wins on this family. The crossover sits
// around delta = 20.
template <typename T>
TheHeuristic<T> StandardHeuristicDualDescriptionProgram() {
  std::vector<std::string> ListString = {"1", "1 delta < 20 lrs", "normaliz"};
  return HeuristicFromListString<T>(ListString);
}

template <typename T> TheHeuristic<T> StandardHeuristicStabEquiv() {
  std::vector<std::string> ListString = {
      "2", "1 groupsizerelevant < -1 exhaustive",
      "1 groupsizerelevant < 10000 classic", "partition"};
  return HeuristicFromListString<T>(ListString);
}

// Description: When computing with the adjacency method
// one needs to have some initial facets. There are many
// methods for computing an initial facet.
// Things to consider:
// * In most cases the linear programming would give you
// a good enough starting point.
// * In some cases the linear programming would give you
// the most degenerate facet as starting point. This is
// annoying since you do not want to start by the hardest
// one and instead you want to start with the simplest one
// * The hardest one can sometimes be circumvented by
// early termination criterion. This is an additional reason
// not to start with the hardest.
// * The most sophisticated method like "sampling" takes
// a lot of time to run and you typically do not want that.
// * Sometimes it just pays off to increase the nuber of
// iterations.
//
// Variable available for making the decision:
// * "incidence": High incidence means that classic algorithm will be slower.
// * "groupsize": The size of the group available for the computation
// * "level": The depth of the calls to the recursive adjacency
//    decomposition.
// * "rank": the dimension of the polytope in question
// * "delta": the difference between the incidence and the rank.
//
// Possible output:
// * "lp_cdd": The linear programming iterated 10 times.
// * "lp_cdd:iter_100": The linear programming iterated 100 times
// * "lp_cdd_min": The linear programming iterated 10 times
//    but keeping only the ones of minimal incidence (recommended choice)
// * "lp_cdd_min:iter_100": Same but iterated 100 times
// * "sampling": A sophisticated sampling approach that is a sibling
//   of the recursive adjacency decomposition but does that for sampling.
//   Expensive, but allows to find facet of low incidence on the hardest
//   stuff
// * "lrs_limited": lrs but limited to the first 100 choices
// * "lrs_limited:upper_limit_1000": same but for 1000
template <typename T> TheHeuristic<T> MethodInitialFacetSet() {
  std::vector<std::string> ListString = {"2", "1 incidence < 100 lp_cdd",
                                         "1 incidence < 1000 lp_cdd:iter_100",
                                         "lp_cdd:iter_200"};
  return HeuristicFromListString<T>(ListString);
}

// Description: When computing with the dual description the use of the
// database bank is supposed to help us. But if we check for every facet
// then it is going to be so very slow.
// Things to consider:
// * Essentially the incidence has to be used, maybe also the rank and delta
// * It has to match somehow with the banksave except for the runtime of
//   of course.
//
// Variable available for making the decision:
// * "incidence": High incidence means that classic algorithm will be slower.
// * "groupsize": The size of the group available for the computation
// * "level": The depth of the calls to the recursive adjacency
//    decomposition.
// * "rank": the dimension of the polytope in question
// * "delta": the difference between the incidence and the rank.
//
// Possible output:
// * "yes": Do query the database bank
// * "no": Do not query the database
template <typename T> TheHeuristic<T> MethodCheckDatabaseBank() {
  std::vector<std::string> ListString = {"1", "1 incidence > 50 yes", "no"};
  return HeuristicFromListString<T>(ListString);
}

// Description: When using the canonic scheme for storing the database
// for the Adjacency Decomposition Method we have to choose the canonic
// method.
// Things to consider:
// * Different canonicalization methods produce different results.
// * The standard technique is "canonic" that first computes the stabilizer
// for computing the canonical form. So there is that cost builtin.
// But it is the best method if most faces have a relatively big stabilizer.
// * We can also use a trivial stabilizer for computing the stabilizer.
// It works best if the faces have trivial or small stabilizer.
// * We can also iterate over the group elements which makes sense if
// the group is small.
//
// Variable available for making the decision:
// * "incidence": High incidence means that classic algorithm will be slower.
// * "groupsize": The size of the group available for the computation
// * "level": The depth of the calls to the recursive adjacency
//    decomposition.
// * "rank": the dimension of the polytope in question
// * "delta": the difference between the incidence and the rank.
//
// Possible output:
// * "canonic": Use the canonicalization with the full stabilizer
// * "canonic_initial_triv": Use the canonicalization with a subgroup of the
//   stabilizer that is trivial.
// * "store": Iterate over all the group elements.
// * "guess": If there is an existing database, we use the first ones
//   (at most 100) for runtime evaluation. If none, we randomly generate some
//   as test case.
// * "load": Load from the database. If info missing then use "canonic"
//   If database missing
template <typename T> TheHeuristic<T> MethodChoiceCanonicalization() {
  std::vector<std::string> ListString = {"1", "1 groupsize < 500 store",
                                         "canonic"};
  return HeuristicFromListString<T>(ListString);
}

// Description: When computing with dual description, we need to split
// from the list of orbits for one big group to the list of orbits for
// a smaller group. There are many possible methods to consider.
// Things to consider:
// * If the index is small, then the "single_cosets" is the best one
// * If the small group is really small then "exhaustive" method is
// * The "canonic" method is an all around reasonably efficient method.
// * We have not yet implemented the algorithms from GAP for doing
//   this computation. There is room for improvement here.
//
// Variable available for making the decision:
// * "groupsize_big": The big group size
// * "groupsize_sma": The small group size
// * "index": The index of the small group in the big one.
// * "n_orbit": The number of orbits used.
//
// Possible output:
// * "guess": Apply a number of guesses in order to compute the
//    the best method.
// * "repr": Using the repr for splitting the orbit using generators.
// * "canonic": Using the canonic for splitting the orbit using generators.
// * "canonic_initial_triv": Using the canonic with trivial initial
//   stabilizer for splitting the orbit using generators.
// * "exhaustive_std": exchaustive splitting using the std::unordered_map
// * "exhaustive_sparse": exchaustive splitting using the tsl::sparse_set
// * "exhaustive_robin": exchaustive splitting using the tsl::robin_set
// * "exhaustive_hopscotch": exchaustive splitting using the tsl::hopscotch_set
// * "exchaustive": exhaustive using exhaustive_robin, that is generate the
//    full orbit for the big group and then split it using the small group.
// * "guess": The beauty of orbit splitting is that the data is available
//    when we choose the method, so we can do some tests before embarking
//    one one specific method.
template <typename T> TheHeuristic<T> MethodOrbitSplitTechnique() {
  std::vector<std::string> ListString = {
      "2", "1 groupsize_sma < 100 exhaustive", "1 index < 20 single_cosets",
      "canonic"};
  return HeuristicFromListString<T>(ListString);
}

// Description: When computing with the adjacency decomposition method,
// we have different possible methods.
// Things to consider:
// * "repr" is here mostly for experimental reasons, it is almost never used.
// * "canonic" is the method to use.
//
// Variable available for making the decision:
// * "incidence": High incidence means that classic algorithm will be slower.
// * "groupsize": The size of the group available for the computation
// * "level": The depth of the calls to the recursive adjacency
//    decomposition.
// * "rank": the dimension of the polytope in question
// * "delta": the difference between the incidence and the rank.
//
// Possible output:
// * "repr": Use the RepresentativeAction for testing equivalence of faces.
// * "canonic": Use the Canonic representative for faces
template <typename T> TheHeuristic<T> MethodChosenDatabase() {
  std::vector<std::string> ListString = {
      "2", "2 groupsize > 5000000000000 incidence > 1000000 repr",
      "2 groupsize > 5000000000000 incidence > 100000 repr", "canonic"};
  return HeuristicFromListString<T>(ListString);
}

// Description: When choosing the dual description method, it can be hard
// to decide what is the best technique. Thus, the Thompson sampling
// can be used for better running of that sort of things.
// Things to consider:
// * We can take an heuristic and generate the Thompson sampling input
//   file by using the program "Convert_Heuristic_to_Thompson_sampling"
// * The logs are quite important to see what is happening.
// * It is still an experimental feature.
template <typename T>
FullNamelist StandardHeuristicDualDescriptionProgram_TS() {
  std::vector<std::string> lstr_proba = {"&PROBABILITY_DISTRIBUTIONS",
                                         " ListName = \"distri1\" ",
                                         " ListNmax = 25 ",
                                         " ListNstart = 2",
                                         " ListNature = \"dirac\"",
                                         " ListDescription = \"0.0\"",
                                         "/"};
  //
  std::vector<std::string> lstr_thompson_prior{"&THOMPSON_PRIOR"};
  std::string s1 = " ListAnswer = \"cdd\", \"lrs\"";
  std::string s2 = " ListName = \"only_cdd\", \"only_lrs\"";
  std::string s3 = " ListDescription = \"cdd:distri1\", \"lrs:distri1\"";
  lstr_thompson_prior.push_back(s1);
  lstr_thompson_prior.push_back(s2);
  lstr_thompson_prior.push_back(s3);
  lstr_thompson_prior.push_back("/");
  //
  std::vector<std::string> lstr_key = {"&KEY_COMPRESSION",
                                       " ListKey = \"delta\"",
                                       " ListDescription = \"superfine\"", "/"};
  //
  std::vector<std::string> lstr_heuristic_prior = {
      "&HEURISTIC_PRIOR", " DefaultPrior = \"noprior:10\"",
      " ListFullCond = \"delta > 30\""};
  lstr_heuristic_prior.push_back(" ListConclusion = \"only_cdd\"");
  lstr_heuristic_prior.push_back("/");
  //
  std::vector<std::string> lstr_io = {"&IO",
                                      " name = \"split\"",
                                      " WriteLog = F",
                                      " ProcessExistingDataIfExist = F",
                                      " LogFileToProcess = \"input_logfile\"",
                                      "/"};
  //
  // Putting things together
  //
  std::vector<std::string> lstr;
  for (auto &e_lstr : {lstr_proba, lstr_thompson_prior, lstr_key,
                       lstr_heuristic_prior, lstr_io}) {
    lstr.push_back("");
    for (auto &estr : e_lstr)
      lstr.push_back(estr);
  }
  FullNamelist eFull = NAMELIST_ThompsonSamplingRuntime();
  NAMELIST_ReadListString(eFull, lstr);
  return eFull;
}

template <typename TintGroup>
void SetThompsonSampling(FullNamelist const &eFull,
                         std::string const &NamelistEnt,
                         ThompsonSamplingHeuristic<TintGroup> &eTS,
                         std::ostream &os) {
  SingleBlock const &BlockHEU = eFull.get_block("HEURISTIC");
  std::string const &NamelistEntFile = BlockHEU.get_string(NamelistEnt);
  if (NamelistEntFile != "unset.ts") {
#ifdef DEBUG_HEURISTICS
    os << "NamelistEntFile for Thompson sampling=" << NamelistEntFile << "\n";
#endif
    FILE_IsExistingFileDie(NamelistEntFile);
    try {
      FullNamelist eFullN = NAMELIST_ThompsonSamplingRuntime();
      NAMELIST_ReadNamelistFile(NamelistEntFile, eFullN);
      eTS = ThompsonSamplingHeuristic<TintGroup>(eFullN, os);
    } catch (TerminalException const &e) {
      std::cerr << "Failed in reading the file NamelistEntFile=" << NamelistEnt
                << "\n";
      throw TerminalException{1};
    }
  }
}

// How the polytope bank is shared across the computation: kept in-process
// ("serial"), served over asio to a separate bank process ("bank_asio"), or
// distributed over MPI ranks ("bank_mpi", MPI runs only).
enum class BankParallelizationMethod { serial, bank_asio, bank_mpi };

inline std::string
bank_parallelization_method_to_string(BankParallelizationMethod method) {
  switch (method) {
  case BankParallelizationMethod::serial:
    return "serial";
  case BankParallelizationMethod::bank_asio:
    return "bank_asio";
  case BankParallelizationMethod::bank_mpi:
    return "bank_mpi";
  }
  return "unknown";
}

inline BankParallelizationMethod
bank_parallelization_method_from_string(std::string const &method) {
  if (method == "serial")
    return BankParallelizationMethod::serial;
  if (method == "bank_asio")
    return BankParallelizationMethod::bank_asio;
  if (method == "bank_mpi")
    return BankParallelizationMethod::bank_mpi;
  std::cerr << "HEU: unknown bank_parallelization_method=" << method << "\n";
  throw TerminalException{1};
}

template <typename TintGroup> struct PolyHeuristicSerial {
  DualDescHeuristic<TintGroup> Splitting;
  DualDescHeuristic<TintGroup> BankSave;
  DualDescHeuristic<TintGroup> AdditionalSymmetry;
  ThompsonSamplingHeuristic<TintGroup> DualDescriptionProgram;
  DualDescHeuristic<TintGroup> InitialFacetSet;
  DualDescHeuristic<TintGroup> CheckDatabaseBank;
  DualDescHeuristic<TintGroup> ChosenDatabase;
  OrbitSplitHeuristic<TintGroup> OrbitSplitTechnique;
  DualDescHeuristic<TintGroup> CommThread;
  DualDescHeuristic<TintGroup> ChoiceCanonicalization;
  bool DD_Saving;
  bool AdvancedTerminationCriterion;
  bool SimpleExchangeScheme;
  SingletonTime start;
  int max_runtime;
  uint16_t port;
  bool BANK_Saving;
  std::string BANK_Prefix;
  std::string OutFormat;
  std::string OutFile;
  BankParallelizationMethod bank_parallelization_method;
  std::string DD_Prefix;
  int dimEXT;
  bool DeterministicRuntime;
};

template <typename T, typename TintGroup>
PolyHeuristicSerial<TintGroup> AllStandardHeuristicSerial(int const &dimEXT,
                                                          std::ostream &os) {
  FullNamelist eFull = StandardHeuristicDualDescriptionProgram_TS<T>();
  bool DD_Saving = false;
  bool AdvancedTerminationCriterion = false;
  bool SimpleExchangeScheme = false;
  int max_runtime = -1;
  uint16_t port = 1234;
  bool BANK_Saving = false;
  std::string BANK_Prefix = "/unset/";
  std::string OutFormat = "GAP";
  std::string OutFile = "unset.out";
  BankParallelizationMethod bank_parallelization_method =
      BankParallelizationMethod::serial;
  std::string DD_Prefix = "/irrelevant/";
  bool DeterministicRuntime = true;
  return {convert_dual_desc_heuristic(StandardHeuristicSplitting<TintGroup>()),
          convert_dual_desc_heuristic(StandardHeuristicBankSave<TintGroup>()),
          convert_dual_desc_heuristic(
              StandardHeuristicAdditionalSymmetry<TintGroup>()),
          ThompsonSamplingHeuristic<TintGroup>(eFull, os),
          convert_dual_desc_heuristic(MethodInitialFacetSet<TintGroup>()),
          convert_dual_desc_heuristic(MethodCheckDatabaseBank<TintGroup>()),
          convert_dual_desc_heuristic(MethodChosenDatabase<TintGroup>()),
          convert_typed_heuristic<TintGroup, OrbitSplitField>(
              MethodOrbitSplitTechnique<TintGroup>()),
          convert_dual_desc_heuristic(StandardHeuristicCommThread<TintGroup>()),
          convert_dual_desc_heuristic(MethodChoiceCanonicalization<TintGroup>()),
          DD_Saving,
          AdvancedTerminationCriterion,
          SimpleExchangeScheme,
          SingletonTime(),
          max_runtime,
          port,
          BANK_Saving,
          BANK_Prefix,
          OutFormat,
          OutFile,
          bank_parallelization_method,
          DD_Prefix,
          dimEXT,
          DeterministicRuntime};
}

template <typename Tint>
void PrintPolyHeuristicSerial(PolyHeuristicSerial<Tint> const &AllArr,
                              std::ostream &os) {
  os << "RDD: PrintPolyHeuristicSerial, beginning\n";
  os << "RDD: OutFile=" << AllArr.OutFile << "\n";
  os << "RDD: OutFormat=" << AllArr.OutFormat << "\n";
  os << "RDD: DeterministicRuntime=" << AllArr.DeterministicRuntime << "\n";
  os << "RDD: port=" << AllArr.port << "\n";
  os << "RDD: bank_parallelization_method="
     << bank_parallelization_method_to_string(AllArr.bank_parallelization_method)
     << "\n";
  os << "RDD: SplittingHeuristicFile\n" << AllArr.Splitting << "\n";
  os << "RDD: AdditionalSymmetryHeuristicFile\n"
     << AllArr.AdditionalSymmetry << "\n";
  os << "RDD: DualDescriptionThompsonFile\n"
     << AllArr.DualDescriptionProgram << "\n";
  os << "RDD: MethodInitialFacetSetFile\n" << AllArr.InitialFacetSet << "\n";
  os << "RDD: BankSaveHeuristicFile\n" << AllArr.BankSave << "\n";
  os << "RDD: CheckDatabaseBank\n" << AllArr.CheckDatabaseBank << "\n";
  os << "RDD: ChosenDatabase\n" << AllArr.ChosenDatabase << "\n";
  os << "RDD: OrbitSplitTechnique\n" << AllArr.OrbitSplitTechnique << "\n";
  os << "RDD: CommThreadHeuristicFile\n" << AllArr.CommThread << "\n";
  os << "RDD: ChoiceCanonicalizationFile\n"
     << AllArr.ChoiceCanonicalization << "\n";
  os << "RDD: max_runtime=" << AllArr.max_runtime << "\n";
  os << "RDD: DD_Saving=" << AllArr.DD_Saving << "\n";
  os << "RDD: DD_Prefix=" << AllArr.DD_Prefix << "\n";
  os << "RDD: BANK_Saving=" << AllArr.BANK_Saving << "\n";
  os << "RDD: BANK_Prefix=" << AllArr.BANK_Prefix << "\n";
  os << "RDD: AdvancedTerminationCriterion="
     << AllArr.AdvancedTerminationCriterion << "\n";
  os << "RDD: SimpleExchangeScheme=" << AllArr.SimpleExchangeScheme << "\n";
}

// clang-format off
#endif  // SRC_POLY_POLY_HEURISTICS_H_
// clang-format on

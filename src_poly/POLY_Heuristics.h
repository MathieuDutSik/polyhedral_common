// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_HEURISTICS_H_
#define SRC_POLY_POLY_HEURISTICS_H_

// clang-format off
#include "Basic_file.h"
#include "Heuristic_ThompsonSampling.h"
#include "Namelist.h"
#include <limits>
#include <string>
#include <vector>
// clang-format on

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
template <typename T> TheHeuristic<T> StandardHeuristicSplitting() {
  std::vector<std::string> ListString = {"4",
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
// * Some algorithm also have external function like glrs and ppl_lcdd
// which tend to be faster but there is an overhead to calling external
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
// * "cdd_cbased": Available only for mpq_class, uses the C-based version
//  by C-linking
// * "cdd": The templatized cdd (only for fields)
// * "lrs_ring": The templatized lrs (only for fields) but fractions reduced
//   to integers or more precisely just ring operations.
// * "lrs": The templatized lrs (for fields and rings)
// * "glrs": The external program glrs (only for types implementing Q)
// * "ppl_ext": The external program ppl_lcdd (only for types implementing Q)
// * "cdd_ext": The external program cdd (only for types implementing Q)
// * "normaliz": he external program normaliz (only for types implementing Q)
template <typename T>
TheHeuristic<T> StandardHeuristicDualDescriptionProgram() {
  std::vector<std::string> ListString = {"1", "1 incidence < 50 cdd", "cdd"};
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

template <typename T>
inline typename std::enable_if<is_implementation_of_Q<T>::value, bool>::type
IsPPLpossible() {
  return IsProgramInPath("ppl_lcdd");
}

template <typename T>
inline typename std::enable_if<!is_implementation_of_Q<T>::value, bool>::type
IsPPLpossible() {
  return false;
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
  std::string s1 = " ListAnswer = \"cdd\", \"lrs_ring\"";
  std::string s2 = " ListName = \"only_cdd\", \"only_lrs\"";
  std::string s3 = " ListDescription = \"cdd:distri1\", \"lrs_ring:distri1\"";
  bool test = IsPPLpossible<T>();
  if (test) {
    s1 += ", \"ppl_ext\"";
    s2 += ", \"only_ppl\"";
    s3 += ", \"ppl_ext:distri1\"";
  }
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
  if (test) {
    lstr_heuristic_prior.push_back(" ListConclusion = \"only_ppl\"");
  } else {
    lstr_heuristic_prior.push_back(" ListConclusion = \"only_cdd\"");
  }
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
void SetHeuristic(FullNamelist const &eFull, std::string const &NamelistEnt,
                  TheHeuristic<TintGroup> &eHeu, std::ostream &os) {
  SingleBlock const& BlockHEU = eFull.get_block("HEURISTIC");
  std::string NamelistEntFile = BlockHEU.get_string(NamelistEnt);
  if (NamelistEntFile != "unset.heu") {
    os << "NamelistEntFile for heuristic=" << NamelistEntFile << "\n";
    IsExistingFileDie(NamelistEntFile);
    std::ifstream is(NamelistEntFile);
    try {
      eHeu = ReadHeuristic<TintGroup>(is);
    } catch (TerminalException const &e) {
      std::cerr << "Failed in reading the file NamelistEntFile=" << NamelistEnt
                << "\n";
      throw TerminalException{1};
    }
  }
}

template <typename TintGroup>
void SetThompsonSampling(FullNamelist const &eFull,
                         std::string const &NamelistEnt,
                         ThompsonSamplingHeuristic<TintGroup> &eTS,
                         std::ostream &os) {
  SingleBlock const& BlockHEU = eFull.get_block("HEURISTIC");
  std::string const& NamelistEntFile = BlockHEU.get_string(NamelistEnt);
  if (NamelistEntFile != "unset.ts") {
    os << "NamelistEntFile for Thompson sampling=" << NamelistEntFile << "\n";
    IsExistingFileDie(NamelistEntFile);
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

template <typename TintGroup> struct PolyHeuristic {
  TheHeuristic<TintGroup> Splitting;
  TheHeuristic<TintGroup> CommThread;
  TheHeuristic<TintGroup> BankSave;
  TheHeuristic<TintGroup> AdditionalSymmetry;
  TheHeuristic<TintGroup> DualDescriptionProgram;
  TheHeuristic<TintGroup> StabEquivFacet;
  TheHeuristic<TintGroup> InitialFacetSet;
  TheHeuristic<TintGroup> ChoiceCanonicalization;
  bool DD_Saving;
  bool DD_Memory;
};

template <typename TintGroup> PolyHeuristic<TintGroup> AllStandardHeuristic() {
  PolyHeuristic<TintGroup> AllArr;
  AllArr.Splitting = StandardHeuristicSplitting<TintGroup>();
  AllArr.CommThread = StandardHeuristicCommThread<TintGroup>();
  AllArr.BankSave = StandardHeuristicBankSave<TintGroup>();
  AllArr.AdditionalSymmetry = StandardHeuristicAdditionalSymmetry<TintGroup>();
  AllArr.DualDescriptionProgram =
      StandardHeuristicDualDescriptionProgram<TintGroup>();
  AllArr.StabEquivFacet = StandardHeuristicStabEquiv<TintGroup>();
  AllArr.InitialFacetSet = MethodInitialFacetSet<TintGroup>();
  AllArr.ChoiceCanonicalization = MethodChoiceCanonicalization<TintGroup>();
  return AllArr;
}

template <typename TintGroup> struct PolyHeuristicSerial {
  TheHeuristic<TintGroup> Splitting;
  TheHeuristic<TintGroup> BankSave;
  TheHeuristic<TintGroup> AdditionalSymmetry;
  ThompsonSamplingHeuristic<TintGroup> DualDescriptionProgram;
  TheHeuristic<TintGroup> InitialFacetSet;
  TheHeuristic<TintGroup> CheckDatabaseBank;
  TheHeuristic<TintGroup> ChosenDatabase;
  TheHeuristic<TintGroup> OrbitSplitTechnique;
  TheHeuristic<TintGroup> CommThread;
  TheHeuristic<TintGroup> ChoiceCanonicalization;
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
  std::string bank_parallelization_method;
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
  std::string bank_parallelization_method = "serial";
  std::string DD_Prefix = "/irrelevant/";
  bool DeterministicRuntime = true;
  return {StandardHeuristicSplitting<TintGroup>(),
          StandardHeuristicBankSave<TintGroup>(),
          StandardHeuristicAdditionalSymmetry<TintGroup>(),
          ThompsonSamplingHeuristic<TintGroup>(eFull, os),
          MethodInitialFacetSet<TintGroup>(),
          MethodCheckDatabaseBank<TintGroup>(),
          MethodChosenDatabase<TintGroup>(),
          MethodOrbitSplitTechnique<TintGroup>(),
          StandardHeuristicCommThread<TintGroup>(),
          MethodChoiceCanonicalization<TintGroup>(),
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

// clang-format off
#endif  // SRC_POLY_POLY_HEURISTICS_H_
// clang-format on

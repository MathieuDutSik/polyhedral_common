// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_HEURISTICS_H_
#define SRC_POLY_POLY_HEURISTICS_H_

#include "Basic_file.h"
#include "Heuristic_ThompsonSampling.h"
#include "Namelist.h"
#include <limits>
#include <string>
#include <vector>

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
  std::vector<std::string> ListString = {"0",
                                         "no"};
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

template <typename T> TheHeuristic<T> MethodInvariantQuality() {
  std::vector<std::string> ListString = {
      "1", "1 groupsizerelevant > 660602000 pairinv:use_pair_orbit", "pairinv"};
  return HeuristicFromListString<T>(ListString);
}

template <typename T> TheHeuristic<T> MethodCheckDatabaseBank() {
  std::vector<std::string> ListString = {"1", "1 incidence > 50 yes", "no"};
  return HeuristicFromListString<T>(ListString);
}

template <typename T> TheHeuristic<T> MethodOrbitSplitTechnique() {
  std::vector<std::string> ListString = {
      "2", "1 groupsize_sma < 100 exhaustive", "1 index < 20 single_cosets",
      "canonic"};
  return HeuristicFromListString<T>(ListString);
}

template <typename T> TheHeuristic<T> MethodChosenDatabase() {
  std::vector<std::string> ListString = {
      "2", "2 groupsize > 5000000000000 incidence > 1000000 traces",
      "2 groupsize > 5000000000000 incidence > 100000 repr", "canonic"};
  return HeuristicFromListString<T>(ListString);
}

FullNamelist StandardHeuristicDualDescriptionProgram_TS() {
  std::vector<std::string> lstr_proba = {"&PROBABILITY_DISTRIBUTIONS",
                                         " ListName = \"distri1\" ",
                                         " ListNmax = 100 ",
                                         " ListNstart = 100",
                                         " ListNature = \"dirac\"",
                                         " ListDescription = \"1:100\"",
                                         "/"};
  //
  std::vector<std::string> lstr_thompson_prior{"&THOMPSON_PRIOR"};
  bool test = IsProgramInPath("ppl_lcdd");
  if (test) {
    lstr_thompson_prior.push_back(
        " ListAnswer = \"cdd\", \"lrs_ring\", \"ppl_ext\"");
    lstr_thompson_prior.push_back(" ListName = \"only_cdd\", \"only_ppl\"");
    lstr_thompson_prior.push_back(
        " ListDescription = \"cdd:distri1\", \"ppl_ext:distri1\"");
  } else {
    lstr_thompson_prior.push_back(" ListAnswer = \"cdd\", \"lrs_ring\"");
    lstr_thompson_prior.push_back(" ListName = \"only_cdd\", \"only_lrs\"");
    lstr_thompson_prior.push_back(
        " ListDescription = \"cdd:distri1\", \"lrs_ring:distri1\"");
  }
  lstr_thompson_prior.push_back("/");
  //
  std::vector<std::string> lstr_key = {
      "&KEY_COMPRESSION", " ListKey = \"incidence\"",
      " ListDescription = "
      "\"1-30,31-35,36-40,41-45,46-50,51-55,56-60,61-65,66-70,71-infinity\"",
      "/"};
  //
  std::vector<std::string> lstr_heuristic_prior = {
      "&HEURISTIC_PRIOR", " DefaultPrior = \"noprior:10\"",
      " ListFullCond = \"incidence > 70\""};
  if (test) {
    lstr_heuristic_prior.push_back(" ListConclusion = \"only_ppl\"");
  } else {
    lstr_heuristic_prior.push_back(" ListConclusion = \"only_cdd\"");
  }
  lstr_heuristic_prior.push_back("/");
  //
  std::vector<std::string> lstr_io = {"&IO",
                                      " name = \"split\"",
                                      " WriteLog = T",
                                      " ProcessExistingDataIfExist = F",
                                      " LogFileToProcess = \"input_logfile\"",
                                      "/"};
  //
  // Putting things together
  //
  std::vector<std::string> lstr = lstr_proba;
  for (auto &e_lstr :
       {lstr_thompson_prior, lstr_key, lstr_heuristic_prior, lstr_io}) {
    lstr.push_back("");
    for (auto &estr : e_lstr)
      lstr.push_back(estr);
  }
  std::cerr << "lstr=\n";
  for (auto estr : lstr)
    std::cerr << estr << "\n";
  FullNamelist eFull = NAMELIST_ThompsonSamplingRuntime();
  NAMELIST_ReadListString(eFull, lstr);
  return eFull;
}

template <typename T>
void SetHeuristic(FullNamelist const &eFull, std::string const &NamelistEnt,
                  TheHeuristic<T> &eHeu, std::ostream &os) {
  SingleBlock BlockHEU = eFull.ListBlock.at("HEURISTIC");
  std::string NamelistEntFile = BlockHEU.ListStringValues.at(NamelistEnt);
  if (NamelistEntFile != "unset.heu") {
    os << "NamelistEntFile for heuristic=" << NamelistEntFile << "\n";
    IsExistingFileDie(NamelistEntFile);
    std::ifstream is(NamelistEntFile);
    try {
      eHeu = ReadHeuristic<T>(is);
    } catch (TerminalException const &e) {
      std::cerr << "Failed in reading the file NamelistEntFile=" << NamelistEnt
                << "\n";
      throw TerminalException{1};
    }
  }
}

template <typename T>
void SetThompsonSampling(FullNamelist const &eFull,
                         std::string const &NamelistEnt,
                         ThompsonSamplingHeuristic<T> &eTS, std::ostream &os) {
  SingleBlock BlockHEU = eFull.ListBlock.at("HEURISTIC");
  std::string NamelistEntFile = BlockHEU.ListStringValues.at(NamelistEnt);
  if (NamelistEntFile != "unset.ts") {
    os << "NamelistEntFile for Thompson sampling=" << NamelistEntFile << "\n";
    IsExistingFileDie(NamelistEntFile);
    try {
      FullNamelist eFullN = NAMELIST_ThompsonSamplingRuntime();
      NAMELIST_ReadNamelistFile(NamelistEntFile, eFullN);
      eTS = ThompsonSamplingHeuristic<T>(eFullN, os);
    } catch (TerminalException const &e) {
      std::cerr << "Failed in reading the file NamelistEntFile=" << NamelistEnt
                << "\n";
      throw TerminalException{1};
    }
  }
}

template <typename T> struct PolyHeuristic {
  TheHeuristic<T> Splitting;
  TheHeuristic<T> CommThread;
  TheHeuristic<T> BankSave;
  TheHeuristic<T> AdditionalSymmetry;
  TheHeuristic<T> DualDescriptionProgram;
  TheHeuristic<T> StabEquivFacet;
  TheHeuristic<T> InitialFacetSet;
  TheHeuristic<T> InvariantQuality;
  bool Saving;
  bool eMemory;
};

template <typename T> PolyHeuristic<T> AllStandardHeuristic() {
  PolyHeuristic<T> AllArr;
  AllArr.Splitting = StandardHeuristicSplitting<T>();
  AllArr.CommThread = StandardHeuristicCommThread<T>();
  AllArr.BankSave = StandardHeuristicBankSave<T>();
  AllArr.AdditionalSymmetry = StandardHeuristicAdditionalSymmetry<T>();
  AllArr.DualDescriptionProgram = StandardHeuristicDualDescriptionProgram<T>();
  AllArr.StabEquivFacet = StandardHeuristicStabEquiv<T>();
  AllArr.InitialFacetSet = MethodInitialFacetSet<T>();
  AllArr.InvariantQuality = MethodInvariantQuality<T>();
  return AllArr;
}

template <typename T> struct PolyHeuristicSerial {
  TheHeuristic<T> Splitting;
  TheHeuristic<T> BankSave;
  TheHeuristic<T> AdditionalSymmetry;
  ThompsonSamplingHeuristic<T> DualDescriptionProgram;
  TheHeuristic<T> InitialFacetSet;
  TheHeuristic<T> CheckDatabaseBank;
  TheHeuristic<T> ChosenDatabase;
  TheHeuristic<T> OrbitSplitTechnique;
  TheHeuristic<T> CommThread;
  bool Saving;
  bool AdvancedTerminationCriterion;
  bool SimpleExchangeScheme;
  SingletonTime start;
  int max_runtime;
  short unsigned int port;
  bool BANK_IsSaving;
  std::string BANK_Prefix;
  std::string OutFormat;
  std::string OUTfile;
  std::string bank_parallelization_method;
  std::string DD_Prefix;
  size_t dimEXT;
};

template <typename T>
PolyHeuristicSerial<T> AllStandardHeuristicSerial(std::ostream &os) {
  FullNamelist eFull = StandardHeuristicDualDescriptionProgram_TS();
  bool Saving = false;
  bool AdvancedTerminationCriterion = false;
  bool SimpleExchangeScheme = false;
  int max_runtime = -1;
  short unsigned int port = 1234;
  bool BANK_IsSaving = false;
  std::string BANK_Prefix = "/unset/";
  std::string OutFormat = "GAP";
  std::string OUTfile = "unset.out";
  std::string bank_parallelization_method = "serial";
  std::string DD_Prefix = "/irrelevant/";
  size_t dimEXT = std::numeric_limits<size_t>::max();
  return {StandardHeuristicSplitting<T>(),
          StandardHeuristicBankSave<T>(),
          StandardHeuristicAdditionalSymmetry<T>(),
          ThompsonSamplingHeuristic<T>(eFull, os),
          MethodInitialFacetSet<T>(),
          MethodCheckDatabaseBank<T>(),
          MethodChosenDatabase<T>(),
          MethodOrbitSplitTechnique<T>(),
          StandardHeuristicCommThread<T>(),
          Saving,
          AdvancedTerminationCriterion,
          SimpleExchangeScheme,
          SingletonTime(),
          max_runtime,
          port,
          BANK_IsSaving,
          BANK_Prefix,
          OutFormat,
          OUTfile,
          bank_parallelization_method,
          DD_Prefix,
          dimEXT};
}

// clang-format off
#endif  // SRC_POLY_POLY_HEURISTICS_H_
// clang-format on

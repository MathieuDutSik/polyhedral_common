#ifndef SRC_POLY_POLY_HEURISTICS_H_
#define SRC_POLY_POLY_HEURISTICS_H_

#include "Basic_file.h"
#include "Heuristic_fct.h"
#include "Namelist.h"
#include <vector>
#include <string>

//
// Heuristic business
//

template <typename T> TheHeuristic<T> StandardHeuristicADM() {
  std::vector<std::string> ListString = {"4",
                                         "2 groupsize > 5000 level < 3 split",
                                         "1 incidence < 40 nosplit",
                                         "1 level > 1 nosplit",
                                         "1 groupsize < 100 nosplit",
                                         "split"};
  return HeuristicFromListString<T>(ListString);
}

template <typename T> TheHeuristic<T> StandardHeuristicAdditionalSymmetry() {
  std::vector<std::string> ListString = {"1", "1 incidence < 50 no", "yes"};
  return HeuristicFromListString<T>(ListString);
}

template <typename T> TheHeuristic<T> StandardHeuristicBankSave() {
  std::vector<std::string> ListString = {"2", "1 time < 300 no",
                                         "1 incidence < 50 no", "yes"};
  return HeuristicFromListString<T>(ListString);
}

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

template <typename T> TheHeuristic<T> MethodChosenDatabase() {
  std::vector<std::string> ListString = {
      "2", "2 groupsize > 5000000000000 incidence > 1000000 traces",
      "2 groupsize > 5000000000000 incidence > 100000 repr", "canonic"};
  return HeuristicFromListString<T>(ListString);
}

template <typename T>
void SetHeuristic(FullNamelist const &eFull, std::string const &NamelistEnt,
                  TheHeuristic<T> &eHeu) {
  SingleBlock BlockHEU = eFull.ListBlock.at("HEURISTIC");
  std::string NamelistEntFile = BlockHEU.ListStringValues.at(NamelistEnt);
  if (NamelistEntFile != "unset.heu") {
    std::cerr << "NamelistEntFile=" << NamelistEntFile << "\n";
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

template <typename T> struct PolyHeuristic {
  TheHeuristic<T> Splitting;
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
  AllArr.Splitting = StandardHeuristicADM<T>();
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
  TheHeuristic<T> DualDescriptionProgram;
  TheHeuristic<T> InitialFacetSet;
  TheHeuristic<T> CheckDatabaseBank;
  TheHeuristic<T> ChosenDatabase;
  bool Saving;
};

template <typename T> PolyHeuristicSerial<T> AllStandardHeuristicSerial() {
  PolyHeuristicSerial<T> AllArr;
  AllArr.Splitting = StandardHeuristicADM<T>();
  AllArr.BankSave = StandardHeuristicBankSave<T>();
  AllArr.AdditionalSymmetry = StandardHeuristicAdditionalSymmetry<T>();
  AllArr.DualDescriptionProgram = StandardHeuristicDualDescriptionProgram<T>();
  AllArr.InitialFacetSet = MethodInitialFacetSet<T>();
  AllArr.CheckDatabaseBank = MethodCheckDatabaseBank<T>();
  AllArr.ChosenDatabase = MethodChosenDatabase<T>();
  return AllArr;
}

// clang-format off
#endif  // SRC_POLY_POLY_HEURISTICS_H_
// clang-format on

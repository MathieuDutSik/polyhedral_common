// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_HEURISTICS_H_
#define SRC_POLY_POLY_HEURISTICS_H_

#include "Basic_file.h"
#include "Heuristic_fct.h"
#include "Namelist.h"
#include <string>
#include <vector>

//
// Heuristic business
//

template <typename T> TheHeuristic<T> StandardHeuristicSplitting() {
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


FullNamelist StandardHeuristicDualDescriptionProgram_TS() {
  std::vector<std::string> ListString = {
    "&PROBABILITY_DISTRIBUTIONS",
    " ListName = \"distri1\" ",
    " ListNmax = 100 ",
    " ListNstart = 100",
    " ListNature = \"dirac\"",
    " ListDescription = \"1:100\"",
    "/",
    "",
    "&THOMPSON_PRIOR",
    " ListAnswer = \"cdd\", \"lrs_ring\", \"ppl_ext\"",
    " ListName = \"only_cdd\", \"only_ppl\"",
    " ListDescription = \"cdd:distri1\", \"ppl_ext:distri1\"",
    "/",
    "",
    "&KEY_COMPRESSION",
    " ListKey = \"incidence\"",
    " ListDescription = \"1-30,31-35,36-40,41-45,46-50,51-55,56-60,61-65,66-70,71-infinity\"",
    "/",
    "",
    "&HEURISTIC_PRIOR",
    " DefaultPrior = \"noprior:10\"",
    " ListFullCond = \"incidence > 70\"",
    " ListConclusion = \"only_ppl\"",
    "/",
    "",
    "&IO",
    " name = \"split\"",
    " WriteLog = T",
    " ProcessExistingDataIfExist = F",
    " LogFileToProcess = \"input_logfile\"",
    "/"};
  FullNamelist eFull = NAMELIST_ThompsonSamplingRuntime();
  NAMELIST_ReadListString(eFull, ListString);
  return eFull;
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

template <typename T>
void SetThompsonSampling(FullNamelist const &eFull, std::string const &NamelistEnt,
                         ThompsonSamplingHeuristic<T> & eTS) {
  SingleBlock BlockHEU = eFull.ListBlock.at("HEURISTIC");
  std::string NamelistEntFile = BlockHEU.ListStringValues.at(NamelistEnt);
  if (NamelistEntFile != "unset.ts") {
    std::cerr << "NamelistEntFile=" << NamelistEntFile << "\n";
    IsExistingFileDie(NamelistEntFile);
    try {
      FullNamelist eFullN = NAMELIST_ThompsonSamplingRuntime();
      NAMELIST_ReadNamelistFile(NamelistEntFile, eFullN);
      eTS = ThompsonSamplingHeuristic<T>(eFullN);
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
  AllArr.Splitting = StandardHeuristicSplitting<T>();
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
  bool Saving;
  bool AdvancedTerminationCriterion;
  SingletonTime start;
  int max_runtime;
  short unsigned int port;
  bool BANK_IsSaving;
  std::string BANK_Prefix;
  std::string OutFormat;
  std::string OUTfile;
  std::string parallelization_method;
  std::string DD_Prefix;
};

template <typename T> PolyHeuristicSerial<T> AllStandardHeuristicSerial() {
  FullNamelist eFull = StandardHeuristicDualDescriptionProgram_TS();
  bool Saving = false;
  bool AdvancedTerminationCriterion = false;
  int max_runtime = -1;
  short unsigned int port = 1234;
  bool BANK_IsSaving = false;
  std::string BANK_Prefix = "/unset/";
  std::string OutFormat = "GAP";
  std::string OUTfile = "unset.out";
  std::string parallelization_method = "serial";
  std::string DD_Prefix = "/irrelevant/";
  return {StandardHeuristicSplitting<T>(),
          StandardHeuristicBankSave<T>(),
          StandardHeuristicAdditionalSymmetry<T>(),
          ThompsonSamplingHeuristic<T>(eFull),
          MethodInitialFacetSet<T>(),
          MethodCheckDatabaseBank<T>(),
          MethodChosenDatabase<T>(),
          Saving, AdvancedTerminationCriterion,
          SingletonTime(), max_runtime,
          port, BANK_IsSaving, BANK_Prefix,
          OutFormat, OUTfile,
          parallelization_method, DD_Prefix};
}

// clang-format off
#endif  // SRC_POLY_POLY_HEURISTICS_H_
// clang-format on

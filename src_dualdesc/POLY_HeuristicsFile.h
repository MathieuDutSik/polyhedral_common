// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_HEURISTICSFILE_H_
#define SRC_DUALDESC_POLY_HEURISTICSFILE_H_

// File-system glue for heuristic configuration:
//   - SetHeuristic                 (moved from POLY_Heuristics.h)
//   - UpdateHeuristicSerial_eFull  (moved from POLY_RecursiveDualDesc.h)
// Both functions pull data out of a FullNamelist and either read an external
// .heu file (SetHeuristic) or thread the values from blocks BANK/DATA/METHOD
// into a PolyHeuristicSerial (UpdateHeuristicSerial_eFull). They live in
// their own header so callers that only need the in-memory heuristic types
// from POLY_Heuristics.h don't have to drag in IsExistingFileDie / ifstream.

// clang-format off
#include "Basic_file.h"
#include "Namelist.h"
#include "POLY_Heuristics.h"
#include <fstream>
#include <string>
// clang-format on

template <typename TintGroup>
void SetHeuristic(FullNamelist const &eFull, std::string const &NamelistEnt,
                  TheHeuristic<TintGroup> &eHeu,
                  [[maybe_unused]] std::ostream &os) {
  SingleBlock const &BlockHEU = eFull.get_block("HEURISTIC");
  std::string NamelistEntFile = BlockHEU.get_string(NamelistEnt);
  if (NamelistEntFile != "unset.heu") {
#ifdef DEBUG_HEURISTICS
    os << "NamelistEntFile for heuristic=" << NamelistEntFile << "\n";
#endif
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

template <typename Tint>
void UpdateHeuristicSerial_eFull(FullNamelist const &eFull,
                                 PolyHeuristicSerial<Tint> &AllArr,
                                 std::ostream &os) {
  //
  SingleBlock const &BlockMETHOD = eFull.get_block("METHOD");
  SingleBlock const &BlockBANK = eFull.get_block("BANK");
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  //
  bool BANK_Saving = BlockBANK.get_bool("Saving");
  AllArr.BANK_Saving = BANK_Saving;
  //
  std::string BANK_Prefix = BlockBANK.get_string("Prefix");
  AllArr.BANK_Prefix = BANK_Prefix;
  //
  std::string OutFile = BlockDATA.get_string("OutFile");
  AllArr.OutFile = OutFile;
  //
  bool DeterministicRuntime = BlockDATA.get_bool("DeterministicRuntime");
  if (!DeterministicRuntime) {
    unsigned seed = get_random_seed();
    srand(seed);
  }
  //
  std::string OutFormat = BlockDATA.get_string("OutFormat");
  AllArr.OutFormat = OutFormat;
  //
  int port_i = BlockDATA.get_int("port");
  uint16_t port = port_i;
  AllArr.port = port;
  //
  std::string bank_parallelization_method =
      BlockDATA.get_string("bank_parallelization_method");
  AllArr.bank_parallelization_method = bank_parallelization_method;
  //
  SetHeuristic(eFull, "SplittingHeuristicFile", AllArr.Splitting, os);
  SetHeuristic(eFull, "AdditionalSymmetryHeuristicFile",
               AllArr.AdditionalSymmetry, os);
  SetThompsonSampling(eFull, "DualDescriptionThompsonFile",
                      AllArr.DualDescriptionProgram, os);
  SetHeuristic(eFull, "MethodInitialFacetSetFile", AllArr.InitialFacetSet, os);
  SetHeuristic(eFull, "BankSaveHeuristicFile", AllArr.BankSave, os);
  SetHeuristic(eFull, "CheckDatabaseBankFile", AllArr.CheckDatabaseBank, os);
  SetHeuristic(eFull, "ChosenDatabaseFile", AllArr.ChosenDatabase, os);
  SetHeuristic(eFull, "OrbitSplitTechniqueFile", AllArr.OrbitSplitTechnique,
               os);
  SetHeuristic(eFull, "CommThreadHeuristicFile", AllArr.CommThread, os);
  SetHeuristic(eFull, "ChoiceCanonicalizationFile",
               AllArr.ChoiceCanonicalization, os);
  //
  bool DD_Saving = BlockMETHOD.get_bool("Saving");
  AllArr.DD_Saving = DD_Saving;
  //
  std::string DD_Prefix = BlockMETHOD.get_string("Prefix");
  AllArr.DD_Prefix = DD_Prefix;
  //
  int max_runtime = BlockDATA.get_int("max_runtime");
  AllArr.max_runtime = max_runtime;
  //
  bool AdvancedTerminationCriterion =
      BlockDATA.get_bool("AdvancedTerminationCriterion");
  AllArr.AdvancedTerminationCriterion = AdvancedTerminationCriterion;
  //
  bool SimpleExchangeScheme = BlockDATA.get_bool("SimpleExchangeScheme");
  AllArr.SimpleExchangeScheme = SimpleExchangeScheme;
  //
#ifdef DEBUG_RECURSIVE_DUAL_DESC
  PrintPolyHeuristicSerial(AllArr, os);
#endif
}

// clang-format off
#endif  // SRC_DUALDESC_POLY_HEURISTICSFILE_H_
// clang-format on

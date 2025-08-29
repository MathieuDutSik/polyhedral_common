// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_PERFECTENUMERATION_H_
#define SRC_PERFECT_PERFECTENUMERATION_H_

#include "Namelist.h"
#include <map>
#include <string>
#include <vector>

FullNamelist NAMELIST_GetStandard_ENUMERATE_PERFECT() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["arithmetic_T"] = "gmp_rational";
  ListStringValues1["arithmetic_Tint"] = "gmp_integer";
  ListIntValues1["n"] = 9;
  ListIntValues1["MaxNumberFlyingMessage"] = 100;
  ListIntValues1["MaxIncidenceTreating"] = 45 + 20;
  ListIntValues1["MaxStoredUnsentMatrices"] = 1000;
  //  ListStringValues1["PrefixDataSave"]="Output_";
  SingleBlock BlockDATA;
  BlockDATA.setListIntValues(ListIntValues1);
  BlockDATA.setListBoolValues(ListBoolValues1);
  BlockDATA.setListDoubleValues(ListDoubleValues1);
  BlockDATA.setListStringValues(ListStringValues1);
  BlockDATA.setListListStringValues(ListListStringValues1);
  ListBlock["DATA"] = BlockDATA;
  // Merging all data
  return FullNamelist(ListBlock);
}

// clang-format off
#endif  // SRC_PERFECT_PERFECTENUMERATION_H_
// clang-format on

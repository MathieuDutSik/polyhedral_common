#ifndef PERFECT_ENUMERATION_H
#define PERFECT_ENUMERATION_H

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
  ListIntValues1["n"] = 9;
  ListIntValues1["MaxNumberFlyingMessage"] = 100;
  ListIntValues1["MaxIncidenceTreating"] = 45 + 20;
  ListIntValues1["MaxStoredUnsentMatrices"] = 1000;
  //  ListStringValues1["PrefixDataSave"]="Output_";
  SingleBlock BlockDATA;
  BlockDATA.ListIntValues = ListIntValues1;
  BlockDATA.ListBoolValues = ListBoolValues1;
  BlockDATA.ListDoubleValues = ListDoubleValues1;
  BlockDATA.ListStringValues = ListStringValues1;
  BlockDATA.ListListStringValues = ListListStringValues1;
  ListBlock["DATA"] = BlockDATA;
  // Merging all data
  return {ListBlock, "undefined"};
}

#endif

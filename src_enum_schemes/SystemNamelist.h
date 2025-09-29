// Copyright (C) 2023 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ENUM_SCHEMES_SYSTEM_NAMELIST_H_
#define SRC_ENUM_SCHEMES_SYSTEM_NAMELIST_H_

// clang-format off
#include "Namelist.h"
#include "MatrixGroupAverage.h"
#include <map>
#include <vector>
#include <string>
// clang-format on

SingleBlock SINGLEBLOCK_Get_System() {
  std::map<std::string, std::string> ListStringValues1_doc;
  std::map<std::string, std::string> ListBoolValues1_doc;
  std::map<std::string, std::string> ListIntValues1_doc;
  ListIntValues1_doc["max_runtime_second"] = "Default: 0\n\
The maximum runtime of the process";
  ListBoolValues1_doc["ApplyStdUnitbuf"] = "Default: false\n\
Whether we apply buffering to writing. Needs to be true for debugging crashes";
  ListBoolValues1_doc["Saving"] = "Default: false\n\
Whether we save the data in the run of the enumeration";
  ListStringValues1_doc["Prefix"] = "The place where data should be written to";
  ListBoolValues1_doc["DeterministicRuntime"] = "Default: F\n\
There is some randomness in several algorithms. With DeterministicRuntime:\n\
T: If you run again the program you will get exactly the same result which is good for debugging\n\
F: Running again the program will get you something different";
  ListStringValues1_doc["OUTfile"] =
      "The file containing the output of the result";
  ListStringValues1_doc["OutFormat"] = "Default: GAP\n\
The formatting used for the output. Possible values:\n\
nothing: Do not write anything\n\
NumberGAP: Only the number of entries created as a GAP readable file\n\
ObjectGAP: The list of objects generated as a GAP readable file\n\
DetailedObjectGAP: The list of objects generated in details as a GAP readable file\n\
ObjectPYTHON: The list of objects generated as a PYTHON readable file\n\
AdjacencyGAP: The adjacencies between the object in a GAP format\n\
AdjacencyPYTHON: The adjacencies between the object in a PYTHON format";
  SingleBlock BlockSYSTEM;
  BlockSYSTEM.setListStringValues_doc(ListStringValues1_doc);
  BlockSYSTEM.setListBoolValues_doc(ListBoolValues1_doc);
  BlockSYSTEM.setListIntValues_doc(ListIntValues1_doc);
  return BlockSYSTEM;
}

// clang-format off
#endif  // SRC_ENUM_SCHEMES_SYSTEM_NAMELIST_H_
// clang-format on

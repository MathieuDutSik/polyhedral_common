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

// The SYSTEM block. The values of OutFormat depend on the function writing
// the output, so its documentation, default included, is given by the caller,
// usually from a doc_OutFormat_* function next to that writer.
SingleBlock SINGLEBLOCK_Get_System(std::string const &doc_OutFormat) {
  std::map<std::string, std::string> ListStringValues1_doc;
  std::map<std::string, std::string> ListBoolValues1_doc;
  std::map<std::string, std::string> ListIntValues1_doc;
  ListIntValues1_doc["max_runtime_second"] = "Default: 0\n\
The maximum runtime of the process in seconds (0 = unlimited)";
  ListBoolValues1_doc["ApplyStdUnitbuf"] = "Default: false\n\
Whether to apply unbuffered output for stdout/stderr. Set to true for debugging crashes";
  ListBoolValues1_doc["Saving"] = "Default: false\n\
Whether to save computation data to disk for checkpointing and recovery";
  ListStringValues1_doc["Prefix"] = "Default: /irrelevant/\n\
The directory prefix where data files should be saved. Only relevant if Saving = T";
  ListBoolValues1_doc["DeterministicRuntime"] = "Default: false\n\
There is some randomness in several algorithms. With DeterministicRuntime:\n\
T: If you run the program again you will get exactly the same result (good for debugging)\n\
F: Running the program again will get you something different";
  ListStringValues1_doc["OutFile"] = "Default: unset.out\n\
The file path for the output results";
  ListStringValues1_doc["OutFormat"] = doc_OutFormat;
  SingleBlock BlockSYSTEM;
  BlockSYSTEM.setListStringValues_doc(ListStringValues1_doc);
  BlockSYSTEM.setListBoolValues_doc(ListBoolValues1_doc);
  BlockSYSTEM.setListIntValues_doc(ListIntValues1_doc);
  return BlockSYSTEM;
}

// clang-format off
#endif  // SRC_ENUM_SCHEMES_SYSTEM_NAMELIST_H_
// clang-format on

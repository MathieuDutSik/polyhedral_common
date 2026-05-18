// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_GRP_GROUPFILE_H_
#define SRC_GROUP_GRP_GROUPFILE_H_

// File-system glue for reading and writing groups: the functions that
// existed in GRP_GroupFct.h but pulled in Basic_file.h (IsExistingFile,
// FindAvailableFileFromPrefix, print_stderr_stdout_file).  Stream-based
// counterparts (ReadGroup, WriteGroup, WriteGroupGAP, ...) remain in
// GRP_GroupFct.h so that translation units that only manipulate groups
// in memory do not have to drag in the file-system header.

// clang-format off
#include "Basic_file.h"
#include "GRP_GroupFct.h"
#include <fstream>
#include <string>
// clang-format on

template <typename Tgroup> Tgroup ReadGroupFile(std::string const &file_name) {
  if (!IsExistingFile(file_name)) {
    std::cerr << "Error in ReadGroupFile\n";
    std::cerr << "file_name=" << file_name << " does not appear to exist\n";
    throw TerminalException{1};
  }
  std::ifstream is(file_name);
  return ReadGroup<Tgroup>(is);
}

template <typename Tgroup>
void PrintRepresentativeAction_OnSets_GRP_f1_f2(Tgroup const &GRP,
                                                Face const &f1,
                                                Face const &f2) {
  std::string prefix = "RepresentativeAction_OnSets_GRP_f1_f2_idx";
  std::string FileOut = FindAvailableFileFromPrefix(prefix);
  std::ofstream os(FileOut);
  WriteGroup(os, GRP);
  auto f_print = [&](Face const &f) -> void {
    for (size_t i = 0; i < f.size(); i++) {
      os << " " << f[i];
    }
    os << "\n";
  };
  f_print(f1);
  f_print(f2);
}

template <typename Tgroup>
void WriteGroupFile(std::string const &eFile, Tgroup const &TheGRP) {
  std::ofstream os(eFile);
  WriteGroup(os, TheGRP);
}

template <typename Tgroup>
void WriteGroupFileGAP(std::string const &eFile, Tgroup const &TheGRP) {
  std::ofstream osf(eFile);
  WriteGroupGAP(osf, TheGRP);
}

template <typename Tgroup>
void WriteGroupFormat(std::string const &FileGroup,
                      std::string const &OutFormat, Tgroup const &TheGRP) {
  auto f_print = [&](std::ostream &os) -> void {
    if (OutFormat == "CPP") {
      return WriteGroup(os, TheGRP);
    }
    if (OutFormat == "GAP") {
      return WriteGroupGAP(os, TheGRP);
    }
    std::cerr
        << "Failed to find a matching entry. Allowed types are GAP and CPP\n";
    throw TerminalException{1};
  };
  print_stderr_stdout_file(FileGroup, f_print);
}

// clang-format off
#endif  // SRC_GROUP_GRP_GROUPFILE_H_
// clang-format on

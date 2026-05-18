// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_CDDLIB_EXTERNAL_H_
#define SRC_POLY_POLY_CDDLIB_EXTERNAL_H_

#include "POLY_cddlib.h"

template <typename T>
LpSolution<T>
CDD_LinearProgramming_External(MyMatrix<T> const &InequalitySet,
                               MyVector<T> const &ToBeMinimized,
                               [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_CDD
  os << "Begin CDD_LinearProgramming_External\n";
#endif
  std::string eStr = random_string(20);
  std::string FileIne = "/tmp/LP_" + eStr + ".ine";
  std::string FileLps = "/tmp/LP_" + eStr + ".lps";
  std::string FileErr = "/tmp/LP_" + eStr + ".error";
  std::string FileDdl = "/tmp/LP_" + eStr + ".ddl";
  std::string FileLog = "/tmp/LP_" + eStr + ".log";
  std::string FileCpp = "/tmp/LP_" + eStr + ".cpp";
  auto CleanFile = [&]() -> void {
    RemoveFileIfExist(FileIne);
    RemoveFileIfExist(FileLps);
    RemoveFileIfExist(FileErr);
    RemoveFileIfExist(FileDdl);
    RemoveFileIfExist(FileLog);
    RemoveFileIfExist(FileCpp);
  };
  CleanFile();
  WriteInputFileCdd(FileIne, InequalitySet, ToBeMinimized);
  std::string FileTestlp2 = "testlp2_gmp";
  int iret1 = RunExternalProgram(FileTestlp2, {FileIne}, std::nullopt, FileLog,
                                 FileErr);
  if (iret1 != 0) {
    std::cerr << "iret1=" << iret1 << "\n";
    std::cerr << "Call to " << FileTestlp2 << " failed\n";
    throw TerminalException{1};
  }
  std::string FilelpcddcleanerCpp = "lpcddcleanerCpp";
  int iret2 = RunExternalProgram(FilelpcddcleanerCpp, {}, FileLog, FileCpp,
                                 std::nullopt);
  if (iret2 != 0) {
    std::cerr << "iret2=" << iret2 << "\n";
    std::cerr << "Call to " << FilelpcddcleanerCpp << "\n";
    throw TerminalException{1};
  }
  std::ifstream is(FileCpp);
  LpSolution<T> eSol;
  bool PrimalDefined;
  is >> PrimalDefined;
  eSol.PrimalDefined = PrimalDefined;
  bool DualDefined;
  is >> DualDefined;
  eSol.DualDefined = DualDefined;
  MyVector<T> DualSolution = ReadVector<T>(is);
  eSol.DualSolution = DualSolution;
  T OptimalValue;
  is >> OptimalValue;
  eSol.OptimalValue = OptimalValue;
  MyVector<T> DirectSolution = ReadVector<T>(is);
  eSol.DirectSolution = DirectSolution;
  CleanFile();
  return eSol;
}

template <typename T>
LpSolution<T> CDD_LinearProgramming_BugSearch(MyMatrix<T> const &TheEXT,
                                              MyVector<T> const &eVect,
                                              std::ostream &os) {
  LpSolution<T> eSol1 = CDD_LinearProgramming(TheEXT, eVect, os);
  LpSolution<T> eSol2 = CDD_LinearProgramming_External(TheEXT, eVect, os);
  if (eSol1.PrimalDefined != eSol2.PrimalDefined ||
      eSol1.DualDefined != eSol2.DualDefined) {
    WriteInputFileCdd("bugSearch.ine", TheEXT, eVect);
    std::cerr << "We find the bug we were after\n";
    throw TerminalException{1};
  }
  return eSol1;
}

#endif  // SRC_POLY_POLY_CDDLIB_EXTERNAL_H_

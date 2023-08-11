// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "POLY_LinearProgramming.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_LinearProgram [DATA_LP] [RESULT]\n";
      std::cerr << "\n";
      std::cerr << "DATA_LP: The file containing the data on the problem "
                   "considered\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    std::ifstream is(argv[1]);
    using T = mpq_class;
    MyMatrix<T> FAC = ReadMatrix<T>(is);
    MyVector<T> eMinimize = ReadVector<T>(is);
    //
    LpSolution<T> eSol = CDD_LinearProgramming(FAC, eMinimize);
    //
    std::ofstream os(argv[2]);
    os << "return rec(FAC:=";
    WriteMatrixGAP(os, FAC);
    os << ",\n eMinimize:=";
    WriteVectorGAP(os, eMinimize);
    os << ",\n PrimalDefined:=" << GAP_logical(eSol.PrimalDefined);
    os << ",\n DualDefined:=" << GAP_logical(eSol.DualDefined);
    os << ", DualSolution:=";
    WriteVectorGAP(os, eSol.DualSolution);
    os << ", OptimalValue:=" << eSol.OptimalValue;
    os << ", DirectSolution:=";
    WriteVectorGAP(os, eSol.DirectSolution);
    os << ", DirectSolutionExt:=";
    WriteVectorGAP(os, eSol.DirectSolutionExt);
    os << ", face:=";
    WriteFaceGAP(os, eSol.eFace);
    os << ", rankDirectSol:=" << eSol.rankDirectSol << ");\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_LinearProgram\n";
    exit(e.eVal);
  }
  runtime(time1);
}

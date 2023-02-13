// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheory.h"
#include "POLY_LinearProgramming.h"

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_SolutionMatNonnegative [FAC] [INEQ]\n";
      std::cerr << "\n";
      std::cerr << "FAC : The list of defining inequalities\n";
      std::cerr << "INEQ: A single inequality\n";
      return -1;
    }
    using T = mpq_class;
    //
    std::string FileFAC = argv[1];
    MyMatrix<T> FAC = ReadMatrixFile<T>(FileFAC);
    //
    std::string FileINEQ = argv[2];
    MyVector<T> Ineq = ReadVectorFile<T>(FileINEQ);
    //
    std::optional<MyVector<T>> opt = SolutionMatNonnegative(FAC, Ineq);
    //
    if (opt) {
      MyVector<T> V = *opt;
      std::cerr << "Found a solution V=" << StringVector(V) << "\n";
    } else {
      std::cerr << "No solution found\n";
    }
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_SolutionMatNonnegative\n";
    exit(e.eVal);
  }
  runtime(time1);
}

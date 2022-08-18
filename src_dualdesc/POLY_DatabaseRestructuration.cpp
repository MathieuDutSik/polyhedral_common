// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
//#include "NumberTheory.h"
#include "Group.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "POLY_RecursiveDualDesc.h"
#include "Permutation.h"

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_DatabaseRestructuration [DatabaseInput] [NprocInput] [DatabaseOutput] [NprocOutput]\n";
      return -1;
    }
    std::string DatabaseI = argv[1];
    int NprocI = ParseScalar<int>(argv[2]);
    std::string DatabaseO = argv[3];
    int NprocO = ParseScalar<int>(argv[4]);

    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the computation, please debug\n";
    exit(e.eVal);
  }
  runtime(time1);
}

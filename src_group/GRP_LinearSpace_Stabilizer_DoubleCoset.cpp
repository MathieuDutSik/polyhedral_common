// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "MatrixGroup.h"
#include "Permutation.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinearSpace_Stabilizer_DoubleCoset [GRP_file] "
                   "[SPA_file] [VGRP_file] [OUT_file]\n";
      std::cerr << "\n";
      std::cerr << "GRP_file    : The file containing the group\n";
      std::cerr << "SPA_file    : The file containing the space\n";
      std::cerr << "VGRP_file   : The file containing the group\n";
      std::cerr << "OUT_file    : The file containing the result\n";
      std::cerr << "              to be read in GAP\n";
      return -1;
    }
    using T = mpq_class;
    using TintGroup = mpz_class;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    //
    std::cerr << "GRP_LinearSpace_Stabilizer_DoubleCoset : Reading input\n";
    std::string GRP_file = argv[1];
    std::string SPA_file = argv[2];
    std::string VGRP_file = argv[3];
    std::string OUT_file = argv[4];
    //
    std::vector<MyMatrix<T>> ListMatrGen = ReadListMatrixFile<T>(GRP_file);
    MyMatrix<T> eLatt = ReadMatrixFile<T>(SPA_file);
    std::vector<MyMatrix<T>> V_gens = ReadListMatrixFile<T>(VGRP_file);
    //
    int n = eLatt.rows();
    GeneralMatrixGroupHelper<T, Telt, TintGroup> helper{n};
    std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>> pair =
        LinearSpace_Stabilizer_DoubleCoset<T, Tgroup>(ListMatrGen, helper,
                                                      eLatt, V_gens, std::cerr);
    //
    std::ofstream os(OUT_file);
    os << "return rec(GRPmatr:=Group(";
    WriteListMatrixGAP(os, pair.first);
    os << "), ListCos:=";
    WriteListMatrixGAP(os, pair.second);
    os << ");\n";
    std::cerr
        << "Normal termination of GRP_LinearSpace_Stabilizer_DoubleCoset\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinearSpace_Stabilizer_DoubleCoset\n";
    exit(e.eVal);
  }
  runtime(time1);
}

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
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "GRP_LinearSpace_Equivalence [GRP_file] [SPA_file] [OUT_file]\n";
      std::cerr << "\n";
      std::cerr << "GRP_file    : The file containing the group\n";
      std::cerr << "SPA_file    : The file containing the space\n";
      std::cerr
          << "OUT_file    : The file containing the result to be read in GAP\n";
      return -1;
    }
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    std::cerr << "GRP_ComputeAut_ListMat_Subset_EXT : Reading input\n";
    std::string GRP_file = argv[1];
    std::string SPA_file = argv[2];
    std::string OUT_file = argv[3];
    //
    std::vector<MyMatrix<T>> ListMatrGen;
    {
      std::ifstream is(GRP_file);
      int nbMat;
      is >> nbMat;
      for (int iMat = 0; iMat < nbMat; iMat++) {
        MyMatrix<T> eMatrGen = ReadMatrix<T>(is);
        ListMatrGen.push_back(eMatrGen);
      }
    }
    //
    MyMatrix<T> eLatt;
    {
      std::ifstream is(SPA_file);
      eLatt = ReadMatrix<T>(is);
    }
    //
    int n = eLatt.rows();
    GeneralMatrixGroupHelper<T, Telt> helper{n};
    std::pair<std::vector<MyMatrix<T>>, CosetDescription<T>> pair =
        LinearSpace_Stabilizer_RightCoset<T, Tgroup>(ListMatrGen, helper, eLatt,
                                                     std::cerr);
    CosetDescription<T> coset = pair.second;
    CosetDescription<T>::iterator iter = coset.begin();
    std::vector<MyMatrix<T>> RightCoset;
    while (iter != coset.end()) {
      RightCoset.push_back(*iter);
      iter++;
    }
    //
    std::ofstream os(OUT_file);
    auto f_print = [&](std::vector<MyMatrix<T>> const &Lmat) -> void {
      os << "[";
      bool IsFirst = true;
      for (auto &eMat : Lmat) {
        if (!IsFirst)
          os << ",\n";
        IsFirst = false;
        WriteMatrixGAP(os, eMat);
      }
      os << "]";
    };
    os << "return rec(GRP:=Group(";
    f_print(pair.first);
    os << "), ListCos:=";
    f_print(RightCoset);
    os << ");\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinearSpace_Stabilizer\n";
    exit(e.eVal);
  }
  runtime(time1);
}

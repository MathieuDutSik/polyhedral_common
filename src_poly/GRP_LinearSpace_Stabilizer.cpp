// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "GRP_GroupFct.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinearSpace_Equivalence [GRP_file] [SPA_file] [OUT_file]\n";
      std::cerr << "\n";
      std::cerr << "GRP_file    : The file containing the group\n";
      std::cerr << "SPA_file    : The file containing the space\n";
      std::cerr << "OUT_file    : The file containing the result to be read in GAP\n";
      return -1;
    }
    using T = mpz_class;
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
    GeneralMatrixGroupHelper<T,Telt> helper(n);
    std::vector<MyMatrix<T>> LGen =
      LinearSpace_Equivalence<T, Tgroup>(ListMatrGen, helper, eLatt);
    //
    {
      if (LGen.size() == 0)
        LGen.push_back(IdentityMat<T>(n));
      std::ofstream os(OUT_file);
      os << "return Group([";
      bool IsFirst = true;
      for (auto & eGen : LGen) {
        if (!IsFirst)
          os ",\n";
        IsFirst = false;
        WriteMatrixGAP(os, eGen);
      }
      os << "]);\n";
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}

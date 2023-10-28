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
      std::cerr << "GRP_LinearSpace_Equivalence [GRP_file] [SPA1_file] "
                   "[SPA2_file] [OUT_file]\n";
      std::cerr << "\n";
      std::cerr << "GRP_file    : The file containing the group\n";
      std::cerr << "SPA1_file   : The file containing the space1\n";
      std::cerr << "SPA2_file   : The file containing the space2\n";
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
    std::string SPA1_file = argv[2];
    std::string SPA2_file = argv[3];
    std::string OUT_file = argv[4];
    //
    std::vector<MyMatrix<T>> ListMatrGen;
    {
      std::ifstream is(GRP_file);
      int nbMat;
      is >> nbMat;
      int n = -1;
      for (int iMat = 0; iMat < nbMat; iMat++) {
        MyMatrix<T> eMatrGen = ReadMatrix<T>(is);
        n = eMatrGen.rows();
        ListMatrGen.push_back(eMatrGen);
      }
      GeneralMatrixGroupHelper<T, Telt> helper{n};
      std::cerr << "|MG|=" << L1normMatrixGroup(helper, ListMatrGen) << "\n";
      std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>> pair =
          LLLMatrixGroupReduction<T, Tint>(helper, ListMatrGen);
      std::vector<MyMatrix<T>> ListGenNew = pair.first;
      std::cerr << "|MGnew|=" << L1normMatrixGroup(helper, ListGenNew) << "\n";
      std::cerr << "Pmat=\n";
      WriteMatrix(std::cerr, pair.second);
    }
    //
    MyMatrix<T> eLatt1;
    {
      std::ifstream is(SPA1_file);
      eLatt1 = ReadMatrix<T>(is);
    }
    //
    MyMatrix<T> eLatt2;
    {
      std::ifstream is(SPA2_file);
      eLatt2 = ReadMatrix<T>(is);
    }
    //
    int n = eLatt2.rows();
    GeneralMatrixGroupHelper<T, Telt> helper{n};
    std::optional<MyMatrix<T>> opt =
      LinearSpace_Equivalence<T, Tgroup>(ListMatrGen, helper, eLatt1, eLatt2, std::cerr);
    //
    {
      std::ofstream os(OUT_file);
      os << "return ";
      if (!opt) {
        os << "fail";
      } else {
        MyMatrix<T> const &eEquiv = *opt;
        WriteMatrixGAP(os, eEquiv);
      }
      os << ";\n";
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinearSpace_Equivalence\n";
    exit(e.eVal);
  }
  runtime(time1);
}

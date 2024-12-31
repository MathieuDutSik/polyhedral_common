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
          << "GRP_LinearSpace_Stabilizer [GRP_file] [SPA_file] [OUT_file]\n";
      std::cerr << "\n";
      std::cerr << "GRP_file    : The file containing the group\n";
      std::cerr << "SPA_file    : The file containing the space\n";
      std::cerr
          << "OUT_file    : The file containing the result to be read in GAP\n";
      return -1;
    }
    using T = mpq_class;
    using TintGroup = mpz_class;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    //
    std::cerr << "GRP_ComputeAut_ListMat_Subset_EXT : Reading input\n";
    std::string GRP_file = argv[1];
    std::string SPA_file = argv[2];
    std::string OUT_file = argv[3];
    //
    std::vector<MyMatrix<T>> ListMatrGen = ReadListMatrixFile<T>(GRP_file);
    //
    MyMatrix<T> eLatt = ReadMatrixFile<T>(SPA_file);
    //
    int n = eLatt.rows();
    GeneralMatrixGroupHelper<T, Telt, TintGroup> helper{n};
    std::vector<MyMatrix<T>> LGen = LinearSpace_Stabilizer<T, Tgroup>(
        ListMatrGen, helper, eLatt, std::cerr);
    //
    {
      if (LGen.size() == 0)
        LGen.push_back(IdentityMat<T>(n));
      std::ofstream os(OUT_file);
      os << "return Group([";
      bool IsFirst = true;
      for (auto &eGen : LGen) {
        if (!IsFirst)
          os << ",\n";
        IsFirst = false;
        WriteMatrixGAP(os, eGen);
      }
      os << "]);\n";
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinearSpace_Stabilizer\n";
    exit(e.eVal);
  }
  runtime(time1);
}

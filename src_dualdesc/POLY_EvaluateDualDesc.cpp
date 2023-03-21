// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "NumberTheoryRealField.h"
#include "POLY_RecursiveDualDesc.h"
#include "Permutation.h"
#include "Group.h"
#include "QuadField.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    if (argc != 6) {
      std::cerr << "POLY_EvaluateDualDesc [ListOrbit] [FileEXT] [Incidence] [n_case] [ListProg]\n";
      return -1;
    }
    std::string FileListOrbit = argv[1];
    std::string FileEXT = argv[2];
    std::string strIncidence = argv[3];
    std::string strNcase = argv[4];
    std::string strListProg = argv[5];
    //
    int TheIncidence = ParseScalar<int>(strIncidence);
    int n_case = ParseScalar<int>(strNcase);
    using T = mpq_class;
    MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
    int n_ext = EXT.rows();
    int n_col = EXT.cols();
    int rnk = RankMat(EXT);
    std::cerr << "n_ext=" << n_ext << " n_col=" << n_col << " rnk=" << rnk << "\n";
    vectface vf(n_ext);
    //
    {
      std::vector<std::string> ListLines = ReadFullFile(FileListOrbit);
      int n_orbit = ParseScalar<int>(ListLines[0]);
      std::cerr << "n_orbit=" << n_orbit << "\n";
      int pos = 0;
      int min_incd = std::numeric_limits<int>::max();
      int max_incd = std::numeric_limits<int>::min();
      for (int i_orbit=0; i_orbit<n_orbit; i_orbit++) {
        std::string eLine = ListLines[1+i_orbit];
        Face f(n_ext);
        std::vector<int> LVal = STRING_Split_Int(eLine, " ");
        for (auto pos : LVal) {
          f[pos - 1] = 1;
        }
        int eIncd = f.count();
        if (eIncd < min_incd)
          min_incd = eIncd;
        if (eIncd > max_incd)
          max_incd = eIncd;
        if (eIncd == TheIncidence && pos < n_case) {
          vf.push_back(f);
          pos++;
        }
      }
      std::cerr << "min_incd=" << min_incd << " max_incd=" << max_incd << "\n";
      std::cerr << "|vf|=" << vf.size() << "\n";
    }
    //
    std::vector<std::string> ListProg = STRING_Split(strListProg, ",");
    for (auto & eProg : ListProg) {
      size_t tot_size = 0;
      HumanTime time;
      for (auto & eFace : vf) {
        MyMatrix<T> EXTface = SelectRow(EXT, eFace);
        MyMatrix<T> EXTfaceRed = ColumnReduction(EXTface);
        vectface ListIncd = DirectFacetComputationIncidence(EXTfaceRed, eProg, std::cerr);
        tot_size += ListIncd.size();
      }
      std::cerr << "Result eProg=" << eProg << " tot_size=" << tot_size << " time=" << time << "\n";
    }
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
}

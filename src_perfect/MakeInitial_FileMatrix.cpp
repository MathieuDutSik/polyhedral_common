// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRational.h"
#include "MAT_Matrix.h"
#include "LatticeStabEquiCan.h"
#include "Namelist.h"
#include "Temp_PerfectForm_Enum.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "MakeInitial_FileMatrix [FileIn] [FileOut]\n";
      throw TerminalException{1};
    }
    //
    using T = mpq_class;
    using Tint = mpz_class;
    //
    std::string FileIn = argv[1];
    std::string FileOut = argv[2];
    //
    std::ifstream is(FileIn);
    std::ofstream os(FileOut);
    int nbPerfect;
    is >> nbPerfect;
    os << nbPerfect << "\n";
    std::cerr << "nbPerfect=" << nbPerfect << "\n";
    for (int iPerfect = 0; iPerfect < nbPerfect; iPerfect++) {
      MyMatrix<T> ePerfect_T = ReadMatrix<T>(is);
      //
      int eStatus = 0;
      //
      Tshortest<T, Tint> eRec =
          T_ShortestVector<T, Tint>(ePerfect_T, std::cerr);
      int incd = (eRec.SHV.rows()) / 2;
      //
      HumanTime time;
      MyMatrix<Tint> B = ComputeCanonicalForm<T, Tint>(ePerfect_T, std::cerr);
      MyMatrix<T> B_T = UniversalMatrixConversion<T,Tint>(B);
      MyMatrix<T> ePerfectCan_T = B_T * ePerfect_T * B_T.transpose();
      std::cerr << "iPerfect=" << iPerfect << " / " << nbPerfect
                << " elapsed_seconds=" << time << "\n";
      //
      os << eStatus << "\n";
      os << incd << "\n";
      WriteMatrix(os, ePerfectCan_T);
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}

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
    using Tmat = mpq_class;
    using Tint = mpz_class;
    //    using Tmat = mpq_class;
    //    using Tint = long;
    // using Tmat = mpq_class;
    //    using Tint = int64_t;
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
      MyMatrix<Tmat> ePerfect_Tmat = ReadMatrix<Tmat>(is);
      //
      int eStatus = 0;
      //
      Tshortest<Tmat, Tint> eRec =
          T_ShortestVector<Tmat, Tint>(ePerfect_Tmat, std::cerr);
      int incd = (eRec.SHV.rows()) / 2;
      //
      HumanTime time;
      MyMatrix<Tmat> eMatCan_Tmat =
          ComputeCanonicalForm<Tmat, Tint>(ePerfect_Tmat, std::cerr).Mat;
      std::cerr << "iPerfect=" << iPerfect << " / " << nbPerfect
                << " elapsed_seconds=" << time << "\n";

      //
      os << eStatus << "\n";
      os << incd << "\n";
      WriteMatrix(os, eMatCan_Tmat);
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}

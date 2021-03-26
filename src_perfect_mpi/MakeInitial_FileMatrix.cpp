#include "Permlib_specific.h"
#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "Namelist.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PerfectForm.h"



int main(int argc, char* argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "MakeInitial_FileMatrix [FileIn] [FileOut]\n";
      throw TerminalException{1};
    }
    //
    //    using Tmat=mpq_class;
    //    using Tint=mpz_class;
    using Tmat=mpq_class;
    using Tint=long;
    //
    std::string FileIn  = argv[1];
    std::string FileOut = argv[2];
    //
    std::ifstream is(FileIn);
    std::ofstream os(FileOut);
    int nbPerfect;
    is >> nbPerfect;
    os << nbPerfect << "\n";
    std::cerr << "nbPerfect=" << nbPerfect << "\n";
    for (int iPerfect=0; iPerfect<nbPerfect; iPerfect++) {
      MyMatrix<Tmat> ePerfect_Tmat = ReadMatrix<Tmat>(is);
      //
      int eStatus=0;
      //
      Tshortest<Tmat,Tint> eRec = T_ShortestVector<Tmat,Tint>(ePerfect_Tmat);
      int incd = (eRec.SHV.rows()) / 2;
      //
      std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
      MyMatrix<Tmat> eMatCan_Tmat = ComputeCanonicalForm<Tmat,Tint>(ePerfect_Tmat).Mat;
      std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
      int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
      std::cerr << "iPerfect=" << iPerfect << " / " << nbPerfect << " elapsed_seconds=" << elapsed_seconds << "\n";

      //
      os << eStatus << "\n";
      os << incd << "\n";
      WriteMatrix(os, eMatCan_Tmat);
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}

// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "MAT_MatrixInt.h"
#include "NumberTheory.h"

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CreateAffineBasis [eMat] [FileSave]\n";
      std::cerr << "\n";
      std::cerr << "eMat: the symmetric matrix which we want to express\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    //
    using T = mpq_class;
    std::ifstream EXTfs(argv[1]);
    MyMatrix<T> eMat = ReadMatrix<T>(EXTfs);
    std::cerr << "After read matrix\n";
    //
    AffineBasisResult eBasRes = ComputeAffineBasis<T>(eMat);
    std::cerr << "result=" << eBasRes.result << "\n";
    if (eBasRes.result == true) {
      int n = eBasRes.ListIdx.size();
      std::cerr << "ListIdx=";
      for (int i = 0; i < n; i++)
        std::cerr << eBasRes.ListIdx[i] << " ";
      std::cerr << "\n";
      //
      MyMatrix<T> eBasis = SelectRow(eMat, eBasRes.ListIdx);
      T eDet = DeterminantMat(eBasis);
      std::cerr << "eDet=" << eDet << "\n";
      std::ofstream BASfs(argv[2]);
      BASfs << "return ";
      WriteMatrixGAP(BASfs, eBasis);
      BASfs << ";\n";
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_CreateAffineBasis\n";
    exit(e.eVal);
  }
  runtime(time1);
}

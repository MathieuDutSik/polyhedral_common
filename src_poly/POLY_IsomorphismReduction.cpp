// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheory.h"
#include "POLY_lrslib.h"

int main(int argc, char *argv[]) {
  try {
    using T = mpq_class;
    //
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    if (argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_IsomorphismReduction [EXT] [FAC] [GRP] [OUT]\n";
      std::cerr << "\n";
      std::cerr << "EXT : The vertices\n";
      std::cerr << "FAC : The facets\n";
      std::cerr << "GRP : The group\n";
      std::cerr << "OUT : The gap readable file of output\n";
      return -1;
    }
    //
    std::string FileEXT = argv[1];
    MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
    //
    std::string FileFAC = argv[2];
    MyMatrix<T> FAC = ReadMatrixFile<T>(FileFAC);
    //
    std::string FileGRP = argv[3];
    Tgroup GRP = ReadGRoupFile<Tgroup>(FileGRP);
    //
    std::string FileOUT = argv[4];
    //
    int n_ext = EXT.rows();
    int n_col = EXT.cols(); 
    vectface vf(n_ext);
    int n_fac = FAC.rows();
    for (int i_fac=0; i_fac<n_fac; i_fac++) {
      Face f(n_ext);
      for (int i_ext=0; i_ext<n_ext; i_ext++) {
        T eSum = 0;
        for (int i_col=0; i_col<n_col; i_col++)
          eSum += EXT(i_ext,i_col) * FAC(i_fac,i_col);
        if (eSum == 0)
          f[i_ext] = 1;
        else
          f[i_ext] = 0;
      }
      vf.push_back(f);
    }
    //
    vectface vf_red = OrbitSplittingSet(vf, GRP);
    //
    VectVectInt_Gap_PrintFile(FileOUT, vf_red);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the computation, please debug\n";
    exit(e.eVal);
  }
}

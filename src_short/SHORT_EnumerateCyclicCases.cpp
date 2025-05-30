// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "rational.h"
#include "Group.h"
#include "NumberTheory.h"
#include "Permutation.h"
#include "SHORT_ShortestConfig.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    using Tfield = mpq_class;
    using Tint = int;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, mpz_class>;
    if (argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "SHORT_EnumerateCyclicCases n d iProc Nproc TheMethod [outfile]\n";
      std::cerr << "\n";
      std::cerr << "n          : The dimension of the lattice\n";
      std::cerr << "d          : The dimension of the code subspace\n";
      std::cerr << "iProc      : The processor index\n";
      std::cerr << "NProc      : The number of processors\n";
      std::cerr << "TheMethod  : can be cdd or glpk_secure\n";
      std::cerr << "[outfile]  : optionary argument for writing GAP readable "
                   "output\n";
      std::cerr
          << "It returns if those parameters are feasible for a LCD code\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    //
    int n;
    sscanf(argv[1], "%d", &n);
    //
    int d;
    sscanf(argv[2], "%d", &d);
    //
    int iProc;
    sscanf(argv[3], "%d", &iProc);
    //
    int NProc;
    sscanf(argv[4], "%d", &NProc);
    //
    std::string TheMethod = argv[5];
    //
    if (NProc < 1) {
      std::cerr << "NProc=" << NProc << "\n";
      std::cerr << "while it should be at least 1\n";
      throw TerminalException{1};
    }
    if (iProc >= NProc || iProc < 0) {
      std::cerr << "NProc=" << NProc << " iProc=" << iProc << "\n";
      std::cerr << "while we should have 0<= iProc < NProc \n";
      throw TerminalException{1};
    }
    std::string outfile = argv[6];
    std::ofstream os(outfile);
    os << "return [";
    bool IsFirst = true;
    //
    std::cerr << "n=" << n << " d=" << d << "\n";
    std::vector<std::vector<int>> ListCases =
        SHORT_GetCandidateCyclic_Optimized(n, d);
    int nbCase = ListCases.size();
    std::cerr << "nbCase=" << nbCase << "\n";
    std::vector<MyMatrix<int>> ListSHV;
    int nbFound = 0;
    for (int iCase = 0; iCase < nbCase; iCase++) {
      int res = iCase % NProc;
      std::cerr << "iCase=" << iCase << " / " << nbCase << " NProc=" << NProc
                << " iProc=" << iProc << " res=" << res << "\n";
      if (res == iProc) {
        std::vector<int> eCase = ListCases[iCase];
        MyMatrix<Tfield> eFrame(n + 1, n);
        for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++) {
            int eVal;
            if (i == j)
              eVal = 1;
            else
              eVal = 0;
            eFrame(i, j) = eVal;
          }
        for (int i = 0; i < n; i++)
          eFrame(n, i) = Tfield(eCase[i]) / Tfield(d);
        MyMatrix<Tfield> eBasis = GetZbasis(eFrame);
        MyMatrix<Tfield> eBasisInv = Inverse(eBasis);
        MyMatrix<int> SHV = UniversalMatrixConversion<int, Tfield>(eBasisInv);
        std::cerr << "SHV : \n";
        WriteMatrix(std::cerr, SHV);
        ReplyRealizability<Tfield, Tint> eRes =
            SHORT_TestRealizabilityShortestFamily<Tfield, Tint, Tgroup>(
                SHV, TheMethod, std::cerr);
        if (eRes.reply) {
          if (IsFirst == false)
            os << ",\n";
          IsFirst = false;
          WriteMatrixGAP(os, SHV);
          nbFound++;
        }
      }
    }
    os << "];\n";
    std::cerr << "nbFound=" << nbFound << "\n";
    std::cerr << "Normal termination of SHORT_EnumerateCyclicCases\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SHORT_EnumerateCyclicCases\n";
    exit(e.eVal);
  }
  runtime(time1);
}

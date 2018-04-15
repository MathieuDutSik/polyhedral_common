#include "StrictPositivity.h"
int main(int argc, char *argv[])
{
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
    std::ifstream EXTfs(argv[1]);
    MyMatrix<mpq_class> eMat=ReadMatrix<mpq_class>(EXTfs);
    std::cerr << "After read matrix\n";
    //
    AffineBasisResult eBasRes=ComputeAffineBasis<mpq_class>(eMat);
    std::cerr << "result=" << eBasRes.result << "\n";
    if (eBasRes.result == true) {
      int n=eBasRes.ListIdx.size();
      std::cerr << "ListIdx=";
      for (int i=0; i<n; i++)
	std::cerr << eBasRes.ListIdx[i] << " ";
      std::cerr << "\n";
      //
      MyMatrix<mpq_class> eBasis=SelectRow(eMat, eBasRes.ListIdx);
      mpq_class eDet=DeterminantMat(eBasis);
      std::cerr << "eDet=" << eDet << "\n";
      std::ofstream BASfs(argv[2]);
      BASfs << "return ";
      WriteMatrixGAP(BASfs, eBasis);
      BASfs << ";\n";
    }
    std::cerr << "Completion of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}

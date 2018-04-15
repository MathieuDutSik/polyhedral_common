#include "Temp_PolytopeEquiStab.h"
int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cerr << "Number of argument is = " << argc << "\n";
    std::cerr << "This program is used as\n";
    std::cerr << "AutomorphismGroupVectorFamily [TheInput] [FileSave]\n";
    std::cerr << "\n";
    std::cerr << "eMat: the symmetric matrix which we want to express\n";
    return -1;
  }
  //
  std::cerr << "Reading input\n";
  //
  std::ifstream INfs(argv[1]);
  int nbMat;
  INfs >> nbMat;
  std::vector<MyMatrix<mpq_class> > ListMatrix(nbMat);
  for (int iMat=0; iMat<nbMat; iMat++) {
    MyMatrix<mpq_class> eMat=ReadMatrix<mpq_class>(INfs);
    ListMatrix[iMat]=eMat;
  }
  MyMatrix<int> SHV=ReadMatrix<int>(INfs);
  INfs.close();
  //
  TheGroupFormat GRP=LinPolytopeGram_Automorphism_Nauty<mpq_class,GraphBitset>(ListMatrix, SHV);
  //
  std::ofstream OUTfs(argv[2]);
  WriteGroupGAP(OUTfs, GRP);
  OUTfs.close();
  //
  std::cerr << "End of the program\n";
}

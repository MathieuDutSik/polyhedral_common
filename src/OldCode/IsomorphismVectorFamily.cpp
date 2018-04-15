#include "Temp_PolytopeEquiStab.h"
int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cerr << "Number of argument is = " << argc << "\n";
    std::cerr << "This program is used as\n";
    std::cerr << "IsomorphismVectorFamily [TheInput] [FileSave]\n";
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
  std::vector<MyMatrix<mpq_class> > ListMatrix1(nbMat);
  for (int iMat=0; iMat<nbMat; iMat++) {
    MyMatrix<mpq_class> eMat=ReadMatrix<mpq_class>(INfs);
    ListMatrix1[iMat]=eMat;
  }
  MyMatrix<int> SHV1=ReadMatrix<int>(INfs);
  std::vector<MyMatrix<mpq_class> > ListMatrix2(nbMat);
  for (int iMat=0; iMat<nbMat; iMat++) {
    MyMatrix<mpq_class> eMat=ReadMatrix<mpq_class>(INfs);
    ListMatrix2[iMat]=eMat;
  }
  MyMatrix<int> SHV2=ReadMatrix<int>(INfs);
  INfs.close();
  //
  ResultEquivWeightMatrix eRes=LinPolytopeGram_Isomorphism_Nauty<mpq_class,GraphBitset>(ListMatrix1, SHV1, ListMatrix2, SHV2);
  //
  std::ofstream OUTfs(argv[2]);
  if (eRes.TheReply == false) {
    OUTfs << "return rec(result:=false);\n";
  }
  else {
    OUTfs << "return rec(result:=false, eList:=[";
    int nbVert=eRes.TheEquiv.size();
    for (int iVert=0; iVert<nbVert; iVert++) {
      if (iVert>0)
	OUTfs << ",";
      int eVal=eRes.TheEquiv.at(iVert)+1;
      OUTfs << eVal;
    }
    OUTfs << "]);\n";
  }
  OUTfs.close();
  //
  std::cerr << "End of the program\n";
}

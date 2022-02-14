#include "Permutation.h"
#include "Group.h"
#include "GRP_GroupFct.h"
#include "Temp_PolytopeEquiStab.h"



#include "TestGroup.h"




int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_PermMatrTest_Automorphism [MatrFile] [FaceFile]\n";
      return -1;
    }
    using T = mpq_class;
    //
    std::cerr << "GRP_ComputeAut_ListMat_Subset_EXT : Reading input\n";
    //
    std::string MatrFile = argv[1];
    std::string FaceFile = argv[2];

    MyMatrix<T> M = ReadMatrixFile<T>(MatrFile);
    Face f = ReadFaceFile(FaceFile);

    TestPolytopeFace(M, f);

  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}

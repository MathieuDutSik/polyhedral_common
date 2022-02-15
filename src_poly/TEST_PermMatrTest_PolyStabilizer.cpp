#include "NumberTheory.h"



#include "TestGroup.h"




int main(int argc, char *argv[])
{
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_PermMatrTest_PolyAutomorphism single [MatrFile] [FaceFile]\n";
      std::cerr << "GRP_PermMatrTest_PolyAutomorphism random [MatrFile] [k]\n";
      std::cerr << "GRP_PermMatrTest_PolyAutomorphism multiple [ListMatrFile] [k]\n";
      return -1;
    }
    using T = mpq_class;
    //
    std::cerr << "GRP_ComputeAut_ListMat_Subset_EXT : Reading input\n";
    //
    auto process=[&]() -> void {
      std::string option = argv[1];
      if (option == "single") {
        std::string MatrFile = argv[2];
        MyMatrix<T> M = ReadMatrixFile<T>(MatrFile);
        std::string FaceFile = argv[3];
        Face f = ReadFaceFile(FaceFile);
        TestPolytopeFace(M, f);
        return;
      }
      if (option == "random") {
        std::string MatrFile = argv[2];
        MyMatrix<T> M = ReadMatrixFile<T>(MatrFile);
        int k;
        sscanf(argv[3], "%d", &k);
        for (int iter=0; iter<100; iter++) {
          Face f = RandomKFace(M.rows(), k);
          TestPolytopeFace(M, f);
        }
        return;
      }
      if (option == "multiple") {
        std::string ListMatrFile = argv[2];
        std::ifstream is(ListMatrFile);
        int k;
        sscanf(argv[3], "%d", &k);
        int n_file;
        is >> n_file;
        for (int i_file=0; i_file<n_file; i_file++) {
          MyMatrix<T> M = ReadMatrix<T>(is);
          for (int iter=0; iter<100; iter++) {
            Face f = RandomKFace(M.rows(), k);
            TestPolytopeFace(M, f);
          }
        }
        return;
      }
      std::cerr << "Allowed options are single, random, multiple\n";
      throw TerminalException{1};
    };
    process();
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
}

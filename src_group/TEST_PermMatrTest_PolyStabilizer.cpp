// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "TestGroup.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 4) {
      std::cerr << "This program is used as\n";
      std::cerr
          << "TEST_PermMatrTest_PolyAutomorphism single [MatrFile] [FaceFile]\n";
      std::cerr << "TEST_PermMatrTest_PolyAutomorphism random [MatrFile] [k]\n";
      std::cerr
          << "TEST_PermMatrTest_PolyAutomorphism multiple [ListMatrFile] [k]\n";
      return -1;
    }
    using T = mpq_class;
    //
    auto process = [&]() -> void {
      std::string option = argv[1];
      if (option == "single") {
        std::string MatrFile = argv[2];
        MyMatrix<T> M = ReadMatrixFile<T>(MatrFile);
        std::string FaceFile = argv[3];
        Face f = ReadFaceFile(FaceFile);
        TestPolytopeFace_Stabilizer(M, f);
        return;
      }
      if (option == "random") {
        std::string MatrFile = argv[2];
        MyMatrix<T> M = ReadMatrixFile<T>(MatrFile);
        int k;
        sscanf(argv[3], "%d", &k);
        for (int iter = 0; iter < 100; iter++) {
          Face f = RandomKFace(M.rows(), k);
          TestPolytopeFace_Stabilizer(M, f);
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
        for (int i_file = 0; i_file < n_file; i_file++) {
          MyMatrix<T> M = ReadMatrix<T>(is);
          for (int iter = 0; iter < 100; iter++) {
            Face f = RandomKFace(M.rows(), k);
            TestPolytopeFace_Stabilizer(M, f);
          }
        }
        return;
      }
      std::cerr << "Allowed options are single, random, multiple\n";
      throw TerminalException{1};
    };
    process();
    std::cerr << "Normal termination of TEST_PermMatrTest_PolyStabilizer\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_PermMatrTest_PolyStabilizer\n";
    exit(e.eVal);
  }
  runtime(time1);
}

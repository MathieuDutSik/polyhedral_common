// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheory.h"
#include "POLY_lrslib.h"
#include "NumberTheoryRealField.h"

template<typename T>
void process(std::string const& eFile) {
  std::ifstream is(eFile);
  MyMatrix<T> EXT = ReadMatrixLrsCdd<T>(is);
  int nbCol = EXT.cols();
  //
  std::cout << "V-representation\n";
  std::cout << "begin\n";
  std::cout << "****** " << nbCol << " rational\n";
  long nVertices = 0;
  std::function<void(T *)> fPrint = [&](T *out) -> void {
    for (int iCol = 0; iCol < nbCol; iCol++)
      std::cout << " " << out[iCol];
    std::cout << "\n";
    nVertices++;
  };
  lrs::Kernel_DualDescription(EXT, fPrint);
  std::cout << "end\n";
  std::cout << "*Total: nvertices=" << nVertices << "\n";
}


int main(int argc, char *argv[]) {
  try {
    if (argc != 2 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Temp_lrs [DATAIN]\n";
      std::cerr << "or\n";
      std::cerr << "Temp_lrs [DATAIN] [DATA_ALGEBRAIC_FIELD]\n";
      std::cerr << "\n";
      std::cerr << "DATAIN : The polyhedral cone inequalities\n";
      std::cerr << "DATA_ALGEBRAIC_FIELD : The algebraic field used\n";
      return -1;
    }
    //
    std::string eFile = argv[1];
    if (argc == 2) {
      process<mpq_class>(eFile);
    } else {
      using T_rat = mpq_class;
      std::string FileAlgebraicField = argv[2];
      if (!IsExistingFile(FileAlgebraicField)) {
        std::cerr << "FileAlgebraicField=" << FileAlgebraicField << " is missing\n";
        throw TerminalException{1};
      }
      HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
      int const idx_real_algebraic_field = 1;
      insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
      using T = RealField<idx_real_algebraic_field>;
      process<T>(eFile);
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}

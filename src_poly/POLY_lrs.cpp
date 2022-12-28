// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheory.h"
#include "POLY_lrslib.h"
#include "QuadField.h"
#include "NumberTheoryRealField.h"

template<typename T>
void process(std::string const& eFile, std::string const& choice) {
  std::ifstream is(eFile);
  MyMatrix<T> EXT = ReadMatrixLrsCdd<T>(is);
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  //
  if (choice == "lrs") {
    std::cout << "V-representation\n";
    std::cout << "begin\n";
    std::cout << "****** " << nbCol << " rational\n";
    long nVertices = 0;
    auto fPrint = [&](T *out) -> void {
      for (int iCol = 0; iCol < nbCol; iCol++)
        std::cout << " " << out[iCol];
      std::cout << "\n";
      nVertices++;
    };
    lrs::Kernel_DualDescription(EXT, fPrint);
    std::cout << "end\n";
    std::cout << "*Total: nvertices=" << nVertices << "\n";
    return;
  }
  if (choice == "vertex_incidence") {
    std::vector<size_t> VertexIncd(nbRow,0);
    T eScal;
    auto fUpdateIncd = [&](T *out) -> void {
      for (int iRow=0; iRow<nbRow; iRow++) {
        eScal = 0;
        for (int iCol=0; iCol<nbCol; iCol++)
          eScal += out[iCol] * EXT(iRow,iCol);
        if (eScal == 0)
          VertexIncd[iRow] += 1;
      }
    };
    lrs::Kernel_DualDescription(EXT, fUpdateIncd);
    std::cout << "VertexIncd=[";
    for (int iRow=0; iRow<nbRow; iRow++) {
      if (iRow > 0)
        std::cout << ",";
      std::cout << VertexIncd[iRow];
    }
    std::cout << "]\n";
    return;
  }
  std::cerr << "Failed to find a matching entry in POLY_lrs\n";
  throw TerminalException{1};
}


int main(int argc, char *argv[]) {
  try {
    if (argc != 4 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_lrs choice rational [DATAIN]\n";
      std::cerr << "POLY_lrs choice Qsqrt5 [DATAIN]\n";
      std::cerr << "POLY_lrs choice Qsqrt2 [DATAIN]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_lrs choice RealAlgebraic [DATAIN] [DATA_ALGEBRAIC_FIELD]\n";
      std::cerr << "\n";
      std::cerr << "DATAIN : The polyhedral cone inequalities\n";
      std::cerr << "DATA_ALGEBRAIC_FIELD : The algebraic field used\n";
      std::cerr << "\n";
      std::cerr << "        -----\n";
      std::cerr << "\n";
      std::cerr << "lrs: exact behavior as in lrs\n";
      std::cerr << "vertex_incidence: the incidence of the vertices is printed\n";
      return -1;
    }
    //
    std::string choice = argv[1];
    std::string eFile = argv[3];
    std::string type = argv[2];
    auto call_lrs=[&]() -> void {
      if (type == "rational") {
        return process<mpq_class>(eFile, choice);
      }
      if (type == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat,5>;
        return process<T>(eFile, choice);
      }
      if (type == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat,2>;
        return process<T>(eFile, choice);
      }
      if (type == "RealAlgebraic") {
        using T_rat = mpq_class;
        if (argc != 5) {
          std::cerr << "Missing the file for the real algebraic field description on input\n";
          throw TerminalException{1};
        }
        std::string FileAlgebraicField = argv[4];
        if (!IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField << " is missing\n";
          throw TerminalException{1};
        }
        HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        using T = RealField<idx_real_algebraic_field>;
        return process<T>(eFile, choice);
      }
      std::cerr << "Failed to find a matching field for type=" << type << "\n";
      std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, RealAlgebraic\n";
      throw TerminalException{1};
    };
    call_lrs();
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}

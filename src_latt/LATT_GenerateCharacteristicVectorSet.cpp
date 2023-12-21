// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "InvariantVectorFamily.h"
// clang-format on

template<typename T, typename Tint>
void process_C(std::string choice, std::string MatFile, std::string OutFile) {
  MyMatrix<T> GramMat = ReadMatrixFile<T>(MatFile);
  auto f=[&]() -> MyMatrix<Tint> {
    if (choice == "shortest") {
      T_shvec_info<T, Tint> info = computeMinimum_GramMat<T,Tint>(GramMat);
      return MatrixFromVectorFamily(info.short_vectors);
    }
    if (choice == "fullrank") {
      return ExtractInvariantVectorFamilyFullRank<T,Tint>(GramMat);
    }
    if (choice == "spanning") {
      return ExtractInvariantVectorFamilyZbasis<T, Tint>(GramMat);
    }
    std::cerr << "Failed to find a matching entry for choice\n";
    throw TerminalException{1};
  };
  MyMatrix<Tint> M = f();
  int n_vect = M.rows();
  int n = GramMat.cols();
  std::cerr << "|M|=" << n_vect << " / " << n << "\n";
  std::map<T, size_t> map;
  for (int i_vect=0; i_vect<n_vect; i_vect++) {
    T eSum = 0;
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        T val1 = UniversalScalarConversion<T,Tint>(M(i_vect, i));
        T val2 = UniversalScalarConversion<T,Tint>(M(i_vect, j));
        eSum += GramMat(i, j) * val1 * val2;
      }
    }
    map[eSum] += 1;
  }
  std::cerr << "norms =";
  for (auto &kv : map) {
    std::cerr << " [" << kv.first << " : " << kv.second << " ]";
  }
  std::cerr << "\n";
  WriteMatrixFile(OutFile, M);
}


template<typename T>
void process_B(std::string const& arithmetic_vec, std::string choice, std::string MatFile, std::string OutFile) {
  if (arithmetic_vec == "mpz_integer") {
    using Tint = mpz_class;
    return process_C<T,Tint>(choice, MatFile, OutFile);
  }
  std::cerr << "process_A failure: No matching entry for arithmetic_mat\n";
  throw	TerminalException{1};  
}

void process_A(std::string const& arithmetic_mat, std::string const& arithmetic_vec, std::string choice, std::string MatFile, std::string OutFile) {
  if (arithmetic_mat = "mpq_rational") {
    using T = mpq_class;
    return process_B<T>(arithmetic_vec, choice, MatFile, OutFile);
  }
  std::cerr << "process_A failure: No matching entry for arithmetic_mat\n";
  throw TerminalException{1};
}



int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_GenerateCharacteristicVectorSet arithmetic_mat arithmetic_vect choice [MaFilet] [OutFile]\n";
      return -1;
    }
    using Tidx = uint32_t;
    //
    std::string arithmetic_mat = argv[1];
    std::string arithmetic_vec = argv[2];
    std::string choice = argv[3];
    std::string MatFile = argv[4];
    std::string OutFile = argv[5];
    std::vector<MyMatrix<T>> ListMat = ReadListMatrixFile<T>(FileListMat);

    MyMatrix<Tint> SHV =
        ExtractInvariantVectorFamilyZbasis<T, Tint>(ListMat[0]);

    MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
    int n_row = SHV_T.rows();
    std::vector<T> Vdiag(n_row, 0);

    const bool use_scheme = true;
    std::vector<std::vector<Tidx>> ListGen =
        GetListGenAutomorphism_ListMat_Vdiag<T, Tfield, Tidx, use_scheme>(
            SHV_T, ListMat, Vdiag, std::cerr);

    std::vector<MyMatrix<Tint>> ListGenEquiv;
    for (auto &eList : ListGen) {
      std::optional<MyMatrix<T>> opt =
          FindMatrixTransformationTest(SHV_T, SHV_T, eList);
      if (!opt) {
        std::cerr << "Failed to find the matrix\n";
        throw TerminalException{1};
      }
      MyMatrix<T> const &M_T = *opt;
      if (!IsIntegralMatrix(M_T)) {
        std::cerr << "Bug: The matrix should be integral\n";
        throw TerminalException{1};
      }
      MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(M_T);
      ListGenEquiv.push_back(M);
    }
    //
    auto prt = [&](std::ostream &os) -> void {
      if (OutFormat == "GAP") {
        os << "return [";
        bool IsFirst = true;
        for (auto &eMat : ListGenEquiv) {
          if (!IsFirst)
            os << ",\n";
          IsFirst = false;
          WriteMatrixGAP(os, eMat);
        }
        os << "];\n";
        return;
      }
      if (OutFormat == "Oscar") {
        os << ListGenEquiv.size() << "\n";
        for (auto &eMat : ListGenEquiv) {
          WriteMatrix(os, eMat);
        }
        return;
      }
      std::cerr << "Failed to find a matching type for OutFormat=" << OutFormat
                << "\n";
      throw TerminalException{1};
    };
    if (OutFile == "stderr") {
      prt(std::cerr);
    } else {
      if (OutFile == "stdout") {
        prt(std::cout);
      } else {
        std::ofstream os(OutFile);
        prt(os);
      }
    }
    std::cerr << "Normal termination of LATT_Automorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Raised exception led to premature end of LATT_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time);
}

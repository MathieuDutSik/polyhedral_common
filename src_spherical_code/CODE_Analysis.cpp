// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "Group.h"
#include "Permutation.h"
#include "PolytopeEquiStab.h"
#include "POLY_RecursiveDualDesc.h"
// clang-format on

size_t get_in(Face const &f) {
  size_t len = f.size();
  for (size_t i = 0; i < len; i++) {
    if (f[i] == 1) {
      return i;
    }
  }
  std::cerr << "Error at get_in\n";
  throw TerminalException{1};
}

size_t get_out(Face const &f) {
  size_t len = f.size();
  for (size_t i = 0; i < len; i++) {
    if (f[i] == 0) {
      return i;
    }
  }
  std::cerr << "Error at get_out\n";
  throw TerminalException{1};
}

template <typename T> MyMatrix<T> DropZeroColumn(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  std::vector<int> lnz;
  auto is_nz = [&](int const &iCol) -> bool {
    for (int iRow = 0; iRow < nbRow; iRow++) {
      if (M(iRow, iCol) != 0) {
        return true;
      }
    }
    return false;
  };
  for (int iCol = 0; iCol < nbCol; iCol++) {
    if (is_nz(iCol)) {
      lnz.push_back(iCol);
    }
  }
  MyMatrix<T> Mred(nbRow, lnz.size());
  int pos = 0;
  for (auto &eCol : lnz) {
    for (int iRow = 0; iRow < nbRow; iRow++) {
      Mred(iRow, pos) = M(iRow, eCol);
    }
    pos += 1;
  }
  return Mred;
}

template <typename T, typename Tgroup>
void process_entry_type(std::string const &FileCode,
                        std::string const &FileGramMat) {
  MyMatrix<T> preCODE = ReadMatrixFile<T>(FileCode);
  int nbEnt = preCODE.rows();
  int nbCol = preCODE.cols();
  std::cerr << "nbEnt=" << nbEnt << " nbCol=" << nbCol << "\n";
  MyMatrix<T> CODE = DropZeroColumn(preCODE);
  std::cerr << "CODE=\n";
  WriteMatrix(std::cerr, CODE);
  int dim = CODE.cols();
  auto get_gram_mat = [&]() -> MyMatrix<T> {
    if (FileGramMat == "identity") {
      return IdentityMat<T>(dim);
    } else {
      return ReadMatrixFile<T>(FileGramMat);
    }
  };
  MyMatrix<T> GramMat = get_gram_mat();
  std::cerr << "GramMat=\n";
  WriteMatrix(std::cerr, GramMat);
  //
  int rnk = RankMat(CODE);
  std::cerr << "nbEnt=" << nbEnt << " dim=" << dim << " rnk=" << rnk << "\n";
  if (dim != rnk) {
    std::cerr << "We have dim != rnk\n";
    throw TerminalException{1};
  }
  std::map<T, size_t> s_norm;
  std::vector<T> norm_by_vertex;
  for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
    MyVector<T> V1 = GetMatrixRow(CODE, iEnt);
    std::map<T, size_t> map;
    for (int jEnt = 0; jEnt < nbEnt; jEnt++) {
      if (iEnt != jEnt) {
        MyVector<T> V2 = GetMatrixRow(CODE, jEnt);
        T scal = ScalarProductQuadForm(GramMat, V1, V2);
        map[scal] += 1;
      }
    }
    std::cerr << "iEnt=" << iEnt;
    for (auto &kv : map) {
      std::cerr << " (" << kv.first << "/" << kv.second << ")";
    }
    std::cerr << "\n";
    //
    T norm = ScalarProductQuadForm(GramMat, V1, V1);
    std::cerr << "iEnt=" << iEnt << " norm=" << norm << "\n";
    std::cerr << "  eLine =";
    for (int i = 0; i < dim; i++) {
      std::cerr << " " << CODE(iEnt, i);
    }
    std::cerr << "\n";
    s_norm[norm] += 1;
    norm_by_vertex.push_back(norm);
  }
  std::cerr << "s_norm =";
  for (auto &kv : s_norm) {
    std::cerr << " (" << kv.first << "," << kv.second << ")";
  }
  std::cerr << "\n";
  for (auto &kv : s_norm) {
    T norm = kv.first;
    std::cerr << "norm=" << norm << " lv=";
    for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
      if (norm_by_vertex[iEnt] == norm) {
        std::cerr << " " << iEnt;
      }
    }
    std::cerr << "\n";
  }
  auto get_main_norm = [&]() -> T {
    for (auto &kv : s_norm) {
      return kv.first;
    }
    return T(0);
  };
  T main_norm = get_main_norm();
  //
  MyMatrix<T> EXT(nbEnt, dim + 1);
  for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
    EXT(iEnt, 0) = 1;
    for (int i = 0; i < dim; i++) {
      EXT(iEnt, i + 1) = CODE(iEnt, i);
    }
  }
  Tgroup GRP =
      LinPolytope_Automorphism_GramMat<T, Tgroup>(CODE, GramMat, std::cerr);
  std::cerr << "|GRP|=" << GRP.size() << "\n";
  //
  vectface vf = DualDescriptionStandard(EXT, GRP);
  std::cerr << "|vf|=" << vf.size() << "\n";
  T MaxCovRadiusSqr = 0;
  double MaxCovRadius = 0;
  bool IsFirst = true;
  for (auto &eFace : vf) {
    std::cerr << "|eFace|=" << eFace.size() << " / " << eFace.count() << "\n";
    MyVector<T> eSol = FindFacetInequality(EXT, eFace);
    MyVector<T> fRay(dim);
    for (int i = 0; i < dim; i++) {
      fRay(i) = eSol(i + 1);
    }
    MyVector<T> eRay = Inverse(GramMat) * fRay;
    T normRay = ScalarProductQuadForm(GramMat, eRay, eRay);
    std::vector<T> ListScal(nbEnt);
    for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
      MyVector<T> V = GetMatrixRow(CODE, iEnt);
      T scal = ScalarProductQuadForm(GramMat, eRay, V);
      ListScal[iEnt] = scal;
    }
    size_t idx_in = get_in(eFace);
    size_t idx_out = get_out(eFace);
    if (ListScal[idx_out] > ListScal[idx_in]) {
      for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
        ListScal[iEnt] = -ListScal[iEnt];
      }
      eRay = -eRay;
    }
    T MaxScal = ListScal[idx_in];
    std::cerr << "eRay=" << StringVector(eRay) << "\n";
    std::cerr << "incd=";
    for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
      if (eFace[iEnt] == 1) {
        std::cerr << " " << iEnt;
      }
    }
    std::cerr << "\n";
    std::cerr << "ListScal=";
    for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
      std::cerr << " " << ListScal[iEnt];
    }
    std::cerr << "\n";
    for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
      if (eFace[iEnt] == 1) {
        if (ListScal[iEnt] != MaxScal) {
          std::cerr << "Error for == MaxScal\n";
          throw TerminalException{1};
        }
      }
      if (eFace[iEnt] == 0) {
        if (ListScal[iEnt] >= MaxScal) {
          std::cerr << "Error for < MaxScal\n";
          throw TerminalException{1};
        }
      }
    }
    T CovRadiusSqr = MaxScal * MaxScal / (normRay * main_norm);
    double MaxScal_d = UniversalScalarConversion<double, T>(MaxScal);
    double normRay_d = UniversalScalarConversion<double, T>(normRay);
    double main_norm_d = UniversalScalarConversion<double, T>(main_norm);
    double CovRadius = MaxScal_d / sqrt(normRay_d * main_norm_d);
    if (IsFirst) {
      MaxCovRadiusSqr = CovRadiusSqr;
      MaxCovRadius = CovRadius;
    } else {
      if (CovRadius < MaxCovRadius) {
        MaxCovRadius = CovRadius;
      }
      if (CovRadiusSqr < MaxCovRadiusSqr) {
        MaxCovRadiusSqr = CovRadiusSqr;
      }
    }
  }
  std::cerr << "MaxCovRadius=" << MaxCovRadius << "\n";
  std::cerr << "MaxCovRadiusSqr=" << MaxCovRadiusSqr << "\n";
  double cov_angle = acos(MaxCovRadius);
  std::cerr << "cov_angle=" << cov_angle << "\n";
}

template <typename Tgroup>
void process_entry(std::string const &arith, std::string const &FileCode,
                   std::string const &FileGramMat) {
  if (arith == "rational") {
    using T = mpq_class;
    return process_entry_type<T, Tgroup>(FileCode, FileGramMat);
  }
  if (arith == "Qsqrt3") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 3>;
    return process_entry_type<T, Tgroup>(FileCode, FileGramMat);
  }
  if (arith == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return process_entry_type<T, Tgroup>(FileCode, FileGramMat);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4) {
      std::cerr << "CODE_Analysis [arith] [FileCODE] [FileGRAM]\n";
      throw TerminalException{1};
    }
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string arith = argv[1];
    std::string FileCode = argv[2];
    std::string FileGramMat = argv[3];
    process_entry<Tgroup>(arith, FileCode, FileGramMat);
    std::cerr << "Normal termination of the program time=" << time << "\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CODE_Analysis time=" << time << "\n";
    exit(e.eVal);
  }
}

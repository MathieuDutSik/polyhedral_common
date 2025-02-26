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
                        std::string const &FileGramMat,
                        std::string const &OutFormat,
                        std::ostream& osf) {
  using TintGroup = typename Tgroup::Tint;
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
    if (norm == 0) {
      std::cerr << "The norm should not be equal to zero\n";
      throw TerminalException{1};
    }
    s_norm[norm] += 1;
    norm_by_vertex.push_back(norm);
  }
  std::cerr << "s_norm =";
  for (auto &kv : s_norm) {
    std::cerr << " (" << kv.first << "," << kv.second << ")";
  }
  std::cerr << "\n";
  if (s_norm.size() > 1) {
    std::cerr << "The number of different norms is incorrect\n";
    throw TerminalException{1};
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
  T MinCosineSqr = 0;
  double MinCosine = 0;
  bool IsFirst = true;
  size_t pos = 0;
  std::vector<MyVector<T>> ListExtremeRay;
  for (auto &eFace : vf) {
    Tgroup eStab = GRP.Stabilizer_OnSets(eFace);
    TintGroup OrbSize = GRP.size() / eStab.size();
    std::cerr << "pos=" << pos << "  |eFace|=" << eFace.size() << " / " << eFace.count() << " |stab|=" << eStab.size() << " |Orb|=" << OrbSize << "\n";
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
    std::cerr << "  eRay=" << StringVector(eRay) << "\n";
    std::cerr << "  incd=";
    for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
      if (eFace[iEnt] == 1) {
        std::cerr << " " << iEnt;
      }
    }
    std::cerr << "\n";
    std::map<T, size_t> map;
    for (auto & val : ListScal) {
      map[val] += 1;
    }
    std::cerr << "  map(ListScal)=";
    for (auto & kv : map) {
      T val = kv.first;
      size_t mult = kv.second;
      double val_d = UniversalScalarConversion<double,T>(val);
      std::cerr << " (" << val << "|" << val_d << "|" << mult << ")";
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
    T CosineSqr = MaxScal * MaxScal / (normRay * main_norm);
    double Scal_d = UniversalScalarConversion<double, T>(MaxScal);
    double normRay_d = UniversalScalarConversion<double, T>(normRay);
    double main_norm_d = UniversalScalarConversion<double, T>(main_norm);
    double Cosine = Scal_d / sqrt(normRay_d * main_norm_d);
    std::cerr << "  Cosine=" << Cosine << " CosineSqr=" << CosineSqr << "\n";
    if (IsFirst) {
      MinCosineSqr = CosineSqr;
      MinCosine = Cosine;
      IsFirst = false;
      ListExtremeRay.push_back(eRay);
    } else {
      if (Cosine == MinCosine) {
        ListExtremeRay.push_back(eRay);
      } else {
        if (Cosine < MinCosine) {
          MinCosine = Cosine;
          MinCosineSqr = CosineSqr;
          ListExtremeRay.clear();
          ListExtremeRay.push_back(eRay);
        }
      }
    }
    pos += 1;
  }
  std::cerr << "MinCosine=" << MinCosine << "\n";
  std::cerr << "MinCosineSqr=" << MinCosineSqr << "\n";
  double cov_angle = acos(MinCosine);
  std::cerr << "cov_angle=" << cov_angle << "\n";
  if (OutFormat == "GAP") {
    osf << "return rec(ListExtremeRay:=[";
    MyMatrix<T> MatExtremeRay = MatrixFromVectorFamily(ListExtremeRay);
    WriteMatrixGAP(osf, MatExtremeRay);
    osf << "]);\n";
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
  throw TerminalException{1};
}

template <typename Tgroup>
void process_entry(std::string const &arith, std::string const &FileCode,
                   std::string const &FileGramMat,
                   std::string const& OutFormat, std::ostream & osf) {
  if (arith == "rational") {
    using T = mpq_class;
    return process_entry_type<T, Tgroup>(FileCode, FileGramMat, OutFormat, osf);
  }
  if (arith == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return process_entry_type<T, Tgroup>(FileCode, FileGramMat, OutFormat, osf);
  }
  if (arith == "Qsqrt3") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 3>;
    return process_entry_type<T, Tgroup>(FileCode, FileGramMat, OutFormat, osf);
  }
  if (arith == "Qsqrt5") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return process_entry_type<T, Tgroup>(FileCode, FileGramMat, OutFormat, osf);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "CODE_Analysis [arith] [FileCODE] [FileGRAM]\n";
      std::cerr << "    or\n";
      std::cerr << "CODE_Analysis [arith] [FileCODE] [FileGRAM] [OutFormat] [FileO]\n";
      std::cerr << "   ---------------\n";
      std::cerr << "arith: rational, Qsqrt2, Qsqrt3, Qsqrt5\n";
      throw TerminalException{1};
    }
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string arith = argv[1];
    std::string FileCode = argv[2];
    std::string FileGramMat = argv[3];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      FileO = argv[5];
    }
    if (FileO == "stderr") {
      process_entry<Tgroup>(arith, FileCode, FileGramMat, OutFormat, std::cerr);
    } else {
      if (FileO == "stdout") {
        process_entry<Tgroup>(arith, FileCode, FileGramMat, OutFormat, std::cout);
      } else {
        std::ofstream osf(FileO);
        process_entry<Tgroup>(arith, FileCode, FileGramMat, OutFormat, osf);
      }
    }
    std::cerr << "Normal termination of the program time=" << time << "\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CODE_Analysis time=" << time << "\n";
    exit(e.eVal);
  }
}

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

size_t get_in(Face const& f) {
  size_t len=f.size();
  for (size_t i=0; i<len; i++) {
    if (f[i] == 1) {
      return i;
    }
  }
  std::cerr << "Error at get_in\n";
  throw TerminalException{1};
}

size_t get_out(Face const& f) {
  size_t len=f.size();
  for (size_t i=0; i<len; i++) {
    if (f[i] == 0) {
      return i;
    }
  }
  std::cerr << "Error at get_out\n";
  throw TerminalException{1};
}


template<typename T, typename Tgroup>
void process_entry_type(std::string const& FileCode) {
  MyMatrix<T> CODE = ReadMatrixFile<T>(FileCode);
  int nbEnt = CODE.rows();
  int dim = CODE.cols();
  std::cerr << "nbEnt=" << nbEnt << " dim=" << dim << "\n";
  std::map<T, size_t> s_norm;
  for (int iEnt=0; iEnt<nbEnt; iEnt++) {
    std::map<T, size_t> map;
    for (int jEnt=0; jEnt<nbEnt; jEnt++) {
      if (iEnt != jEnt) {
        T scal(0);
        for (int i=0; i<dim; i++) {
          scal += CODE(iEnt, i) * CODE(jEnt, i);
        }
        map[scal] += 1;
      }
    }
    std::cerr << "iEnt=" << iEnt;
    for (auto & kv : map) {
      std::cerr << " (" << kv.first << "/" << kv.second << ")";
    }
    std::cerr << "\n";
    //
    T norm(0);
    for (int i=0; i<dim; i++) {
      norm += CODE(iEnt, i) * CODE(iEnt, i);
    }
    std::cerr << "iEnt=" << iEnt << " norm=" << norm << "\n";
    std::cerr << "  eLine =";
    for (int i=0; i<dim; i++) {
      std::cerr << " " << CODE(iEnt,i);
    }
    std::cerr << "\n";
    s_norm[norm] += 1;
  }
  std::cerr << "s_norm =";
  for (auto & kv: s_norm) {
    std::cerr << " (" << kv.first << "," << kv.second << ")";
  }
  std::cerr << "\n";
  auto get_main_norm=[&]() -> T {
    for (auto & kv: s_norm) {
      return kv.first;
    }
    return T(0);
  };
  T main_norm = get_main_norm();
  //
  MyMatrix<T> EXT(nbEnt, dim + 1);
  for (int iEnt=0; iEnt<nbEnt; iEnt++) {
    EXT(iEnt, 0) = 1;
    for (int i=0; i<dim; i++) {
      EXT(iEnt, i+1) = CODE(iEnt,i);
    }
  }
  Tgroup GRP = LinPolytope_Automorphism<T, Tgroup>(EXT, std::cerr);
  std::cerr << "|GRP|=" << GRP.size() << "\n";
  //
  vectface vf = DualDescriptionStandard(EXT, GRP);
  std::cerr << "|vf|=" << vf.size() << "\n";
  T MaxCovRadiusSqr = 0;
  double MaxCovRadius = 0;
  bool IsFirst = true;
  for (auto & eFace : vf) {
    std::cerr << "|eFace|=" << eFace.size() << " / " << eFace.count() << "\n";
    MyVector<T> eSol = FindFacetInequality(EXT, eFace);
    MyVector<T> eRay(dim);
    for (int i=0; i<dim; i++) {
      eRay(i) = eSol(i+1);
    }
    T normRay(0);
    for (int i=0; i<dim; i++) {
      normRay += eRay(i) * eRay(i);
    }
    std::vector<T> ListScal(nbEnt);
    for (int iEnt=0; iEnt<nbEnt; iEnt++) {
      T scal(0);
      for (int i=0; i<dim; i++) {
        scal += eRay(i) * CODE(iEnt,i);
      }
      ListScal[iEnt] = scal;
      std::cerr << "1: iEnt=" << iEnt << " scal=" << scal << " eFace=" << eFace[iEnt] << "\n";
    }
    size_t idx_in = get_in(eFace);
    size_t idx_out = get_out(eFace);
    if (ListScal[idx_out] > ListScal[idx_in]) {
      for (int iEnt=0; iEnt<nbEnt; iEnt++) {
        ListScal[iEnt] = - ListScal[iEnt];
      }
      eRay = - eRay;
    }
    T MaxScal = ListScal[idx_in];
    for (int iEnt=0; iEnt<nbEnt; iEnt++) {
      std::cerr << "2: iEnt=" << iEnt << " scal=" << ListScal[iEnt] << " eFace=" << eFace[iEnt] << "\n";
    }
    std::cerr << "MaxScal=" << MaxScal << "\n";
    for (int iEnt=0; iEnt<nbEnt; iEnt++) {
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
    double MaxScal_d = UniversalScalarConversion<double,T>(MaxScal);
    double normRay_d = UniversalScalarConversion<double,T>(normRay);
    double main_norm_d = UniversalScalarConversion<double,T>(main_norm);
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
void process_entry(std::string const& arith, std::string const& FileCode) {
  if (arith == "rational") {
    using T = mpq_class;
    return process_entry_type<T, Tgroup>(FileCode);
  }
  if (arith == "Qsqrt3") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 3>;
    return process_entry_type<T, Tgroup>(FileCode);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "CODE_Analysis [arith] [FileCODE]\n";
      throw TerminalException{1};
    }
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string arith = argv[1];
    std::string FileCode = argv[2];
    process_entry<Tgroup>(arith, FileCode);
    std::cerr << "Normal termination of the program time=" << time << "\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CODE_Analysis time=" << time << "\n";
    exit(e.eVal);
  }
}

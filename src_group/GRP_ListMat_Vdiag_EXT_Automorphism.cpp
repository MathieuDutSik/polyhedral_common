// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "GRP_GroupFct.h"
#include "Group.h"
#include "Permutation.h"
#include "Temp_PolytopeEquiStab.h"
// clang-format on

template <typename Tfield, typename Twork, typename Tgroup>
void process_inner1(std::string const &FileO, std::string const &OutFormat,
                    MyMatrix<Twork> const &EXT,
                    std::vector<MyMatrix<Twork>> const &ListMat,
                    std::vector<Twork> const &Vdiag) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  const bool use_scheme = true;
  std::vector<std::vector<Tidx>> ListGen =
      GetListGenAutomorphism_ListMat_Vdiag<Twork, Tfield, Tidx, use_scheme>(
          EXT, ListMat, Vdiag, std::cerr);
  //
  std::vector<Telt> LGen;
  for (auto &eList : ListGen)
    LGen.push_back(Telt(eList));
  int n_rows = EXT.rows();
  Tgroup GRP(LGen, n_rows);
  std::cerr << "|GRP|=" << GRP.size() << "\n";
  auto f_print = [&](std::ostream &os) -> void {
    if (OutFormat == "GAP") {
      os << "return " << GRP.GapString() << ";\n";
      return;
    }
    if (OutFormat == "Oscar") {
      WriteGroup(os, GRP);
      return;
    }
    std::cerr << "Failed to find a matching OutFormat=" << OutFormat << "\n";
    throw TerminalException{1};
  };
  if (FileO == "stderr") {
    f_print(std::cerr);
  } else {
    if (FileO == "stdout") {
      f_print(std::cout);
    } else {
      std::ofstream os(FileO);
      f_print(os);
    }
  }
}

template <typename Tfield, typename Twork, typename Tinput, typename Tgroup>
void process_inner2(std::string const &FileO, std::string const &OutFormat,
                    MyMatrix<Tinput> const &EXT,
                    std::vector<MyMatrix<Tinput>> const &ListMat,
                    std::vector<Tinput> const &Vdiag) {
  MyMatrix<Twork> EXT_T = UniversalMatrixConversion<Twork, Tinput>(EXT);
  std::vector<MyMatrix<Twork>> ListMat_T =
      UniversalStdVectorMatrixConversion<Twork, Tinput>(ListMat);
  std::vector<Twork> Vdiag_T =
      UniversalStdVectorScalarConversion<Twork, Tinput>(Vdiag);
  return process_inner1<Tfield, Twork, Tgroup>(FileO, OutFormat, EXT_T,
                                               ListMat_T, Vdiag_T);
}

template <typename Tfield, typename T, typename Tgroup>
void process_inner3(std::string const &FileO, std::string const &OutFormat,
                    MyMatrix<T> const &EXT,
                    std::vector<MyMatrix<T>> const &ListMat,
                    std::vector<T> const &Vdiag) {
  T Linf_EXT = Linfinity_norm_mat(EXT);
  T max_Linf_norm = 0;
  for (auto &eMat : ListMat) {
    T norm = Linfinity_norm_mat(eMat);
    max_Linf_norm = T_max(max_Linf_norm, norm);
  }
  int dim = EXT.cols();
  T max_val = Linf_EXT * Linf_EXT * max_Linf_norm * dim * dim;
  std::string type = get_matching_types(max_val);
  if (type == "int16_t")
    return process_inner2<Tfield, int16_t, T, Tgroup>(FileO, OutFormat, EXT,
                                                      ListMat, Vdiag);
  if (type == "int32_t")
    return process_inner2<Tfield, int32_t, T, Tgroup>(FileO, OutFormat, EXT,
                                                      ListMat, Vdiag);
  if (type == "int64_t")
    return process_inner2<Tfield, int64_t, T, Tgroup>(FileO, OutFormat, EXT,
                                                      ListMat, Vdiag);
  return process_inner2<Tfield, T, T, Tgroup>(FileO, OutFormat, EXT, ListMat,
                                              Vdiag);
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 2 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "GRP_ListMat_Vdiag_EXT_Automorphism [FileI] [OutFormat] [FileO]\n";
      std::cerr << "        or\n";
      std::cerr << "GRP_ListMat_Vdiag_EXT_Automorphism [FileI]\n";
      std::cerr << "\n";
      std::cerr << "FileI     : The file containing the group\n";
      std::cerr << "OutFormat : The output format that can be GAP or Oscar\n";
      std::cerr << "FileO     : The file containg the output\n";
      return -1;
    }
    using Tidx = int32_t;
    using T = mpz_class;
    using Tfield = mpq_class;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using Tint = boost::multiprecision::mpz_int;
#else
    using Tint = mpz_class;
#endif
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    std::cerr << "GRP_ListMat_Vdiag_EXT_Automorphism : Reading input\n";
    //
    std::string FileI = argv[1];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      FileO = argv[3];
    }
    std::ifstream is(FileI);
    int nbMat, len;
    is >> nbMat;
    std::vector<MyMatrix<T>> ListMat;
    for (int iMat = 0; iMat < nbMat; iMat++) {
      MyMatrix<T> eMat = ReadMatrix<T>(is);
      ListMat.push_back(eMat);
    }
    MyMatrix<T> EXT = ReadMatrix<T>(is);
    int dim = EXT.cols();
    for (auto &eMat : ListMat) {
      if (!IsSymmetricMatrix(eMat)) {
        std::cerr << "The matrix eMat should be symmetric\n";
        throw TerminalException{1};
      }
      if (eMat.cols() != dim) {
        std::cerr << "|eMat|=" << eMat.cols() << " |EXT|=" << EXT.cols()
                  << "\n";
        throw TerminalException{1};
      }
    }
    // Pessimistic maximal value
    int n_rows = EXT.rows();
    std::cerr << "n_rows=" << n_rows << "\n";
    is >> len;
    if (len != n_rows) {
      std::cerr << "We have n_rows=" << n_rows << " but len=" << len << "\n";
      throw TerminalException{1};
    }
    std::vector<T> Vdiag(n_rows);
    for (int i = 0; i < n_rows; i++) {
      T val;
      is >> val;
      Vdiag[i] = val;
    }
    //
    process_inner3<Tfield, T, Tgroup>(FileO, OutFormat, EXT, ListMat, Vdiag);
    std::cerr << "Normal termination of GRP_ListMat_Vdiag_EXT_Automorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_ListMat_Vdiag_EXT_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}

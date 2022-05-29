// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Group.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "Permutation.h"

#include "edgewalk.h"
//#include "vinberg.h"

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 2) {
      std::cerr << "Program is used as\n";
      std::cerr << "LORENTZ_SingleEdgewalk [Input]\n";
      throw TerminalException{1};
    }
    std::string FileInput = argv[1];
    //    using T = long;
    //    using Tint = long;
    using T = mpq_class;
    using Tint = mpz_class;
    //    using Tint = mpq_class;
    //    using T = boost::multiprecision::cpp_rational;
    //    using Tint = boost::multiprecision::cpp_int;
    //
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    std::ifstream is(FileInput);
    // Reading the input
    MyMatrix<T> G = ReadMatrix<T>(is);
    //
    std::vector<T> l_norms;
    int n_norms;
    is >> n_norms;
    for (int i = 0; i < n_norms; i++) {
      T val;
      is >> val;
      l_norms.push_back(val);
    }
    //
    MyVector<T> k = ReadVector<T>(is);
    //
    MyVector<Tint> v_disc = ReadVector<Tint>(is);
    MyMatrix<Tint> Mat_l_ui = ReadMatrix<Tint>(is);
    std::vector<MyVector<Tint>> l_ui;
    for (int u = 0; u < Mat_l_ui.rows(); u++) {
      MyVector<Tint> V = GetMatrixRow(Mat_l_ui, u);
      l_ui.push_back(V);
    }
    AdjacencyDirection<Tint> ad{l_ui, v_disc};
    //
    // Reading the expected result
    //
    MyVector<T> gen_expect = ReadVector<T>(is);
    MyMatrix<Tint> MatRoot_expect = ReadMatrix<Tint>(is);
    FundDomainVertex<T, Tint> vert_expect{gen_expect, MatRoot_expect};
    //
    // Doing the edgewalk
    //
    CuspidalBank<T, Tint> cusp_bank;
    FundDomainVertex<T, Tint> vert_obtain =
        EdgewalkProcedure<T, Tint, Tgroup>(cusp_bank, G, l_norms, k, ad);
    //
    // Testing equality
    //
    bool test = vert_obtain == vert_expect;
    if (test) {
      std::cerr << "We obtain the same\n";
    } else {
      std::cerr << "We obtain something different\n";
    }
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
  runtime(time1);
}

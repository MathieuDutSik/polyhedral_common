// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "Group.h"
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "Permutation.h"
#include "edgewalk.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 4) {
      std::cerr << "Program is used as\n";
      std::cerr << "LORENTZ_ComputeStabilizer_Vertex [G] [Vertex]\n";
      std::cerr << "or\n";
      std::cerr << "LORENTZ_ComputeStabilizer_Vertex [G] [Vertex] [outfile]\n";
      throw TerminalException{1};
    }
    std::string FileLorMat = argv[1];
    std::string FileVertex = argv[2];
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
    MyMatrix<T> G = ReadMatrixFile<T>(FileLorMat);
    //
    std::ifstream is2(FileVertex);
    MyVector<T> gen = ReadVector<T>(is2);
    T norm = gen.dot(G * gen);
    std::cerr << "norm=" << norm << "\n";
    MyMatrix<Tint> MatRoot = ReadMatrix<Tint>(is2);
    std::cerr << "|MatRoot| = " << MatRoot.rows() << " / " << MatRoot.cols()
              << " rnk=" << RankMat(MatRoot) << "\n";
    FundDomainVertex<T, Tint> vert{RemoveFractionVector(gen), MatRoot};
    std::cerr << "We have vert\n";
    MyMatrix<Tint> MatRootRed = ColumnReduction(MatRoot);
    std::cerr << "We have MatRootRed\n";
    Tgroup GRP = LinPolytope_Automorphism<Tint, Tgroup>(MatRootRed, std::cerr);
    std::cerr << "|GRP|=" << GRP.size() << "\n";
    //
    std::string OptionNorms = "all";
    std::vector<T> l_norms =
        get_initial_list_norms<T, Tint>(G, OptionNorms, std::cerr);
    std::cerr << "We have l_norms\n";
    SublattInfos<T> si = ComputeSublatticeInfos<T, Tint>(G, l_norms);
    std::cerr << "We have si\n";

    CuspidalBank<T, Tint> cusp_bank;
    TheHeuristic<Tint> HeuristicIdealStabEquiv =
        GetHeuristicIdealStabEquiv<Tint>();
    std::cerr << "We have HeuristicIdealStabEquiv\n";
    FundDomainVertex_FullInfo<T, Tint, Tgroup> vertFull =
        gen_fund_domain_fund_info<T, Tint, Tgroup>(
            cusp_bank, si, vert, HeuristicIdealStabEquiv, std::cerr);
    std::cerr << "We have vertFull\n";
    //
    std::cerr << "Before LORENTZ_GetStabilizerGenerator\n";
    std::vector<MyMatrix<T>> l_mat =
        LORENTZ_GetStabilizerGenerator(G, vertFull, std::cerr);
    std::cerr << "After LORENTZ_GetStabilizerGenerator\n";
    if (argc == 3) {
      for (size_t i_mat = 0; i_mat < l_mat.size(); i_mat++) {
        std::cerr << "i_mat=" << i_mat << "\n";
        WriteMatrix(std::cerr, l_mat[i_mat]);
      }
    }
    if (argc == 4) {
      std::string file = argv[3];
      std::ofstream os(file);
      os << "return [";
      for (size_t i_mat = 0; i_mat < l_mat.size(); i_mat++) {
        if (i_mat > 0)
          os << ",\n";
        os << StringMatrixGAP(l_mat[i_mat]);
      }
      os << "];\n";
    }
    std::cerr << "Normal termination of LORENTZ_ComputeStabilizer_Vertex\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LORENTZ_ComputeStabilizer_Vertex\n";
    exit(e.eVal);
  }
  runtime(time);
}

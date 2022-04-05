#include "Group.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "Permutation.h"

#include "edgewalk.h"
//#include "vinberg.h"

int main(int argc, char *argv[]) {
  SingletonTime time1;
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
    MyMatrix<Tint> MatRoot = ReadMatrix<Tint>(is2);
    std::cerr << "|MatRoot| = " << MatRoot.rows() << " / " << MatRoot.cols()
              << " rnk=" << RankMat(MatRoot) << "\n";
    FundDomainVertex<T, Tint> vert{RemoveFractionVector(gen), MatRoot};
    MyMatrix<Tint> MatRootRed = ColumnReduction(MatRoot);
    Tgroup GRP = LinPolytope_Automorphism<Tint, false, Tgroup>(MatRootRed);
    std::cerr << "|GRP|=" << GRP.size() << "\n";
    //
    std::string OptionNorms = "all";
    std::vector<T> l_norms = get_initial_list_norms<T, Tint>(G, OptionNorms);
    SublattInfos<T> si = ComputeSublatticeInfos<T, Tint>(G, l_norms);
    
    CuspidalBank<T, Tint> cusp_bank;
    TheHeuristic<Tint> HeuristicIdealStabEquiv =
        GetHeuristicIdealStabEquiv<Tint>();
    FundDomainVertex_FullInfo<T, Tint, Tgroup> vertFull =
        gen_fund_domain_fund_info<T, Tint, Tgroup>(cusp_bank, si, vert,
                                                   HeuristicIdealStabEquiv);
    //
    std::vector<MyMatrix<T>> l_mat =
        LORENTZ_GetStabilizerGenerator(G, vertFull);
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
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
  runtime(time1);
}

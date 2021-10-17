#include "NumberTheory.h"
#include "Permutation.h"
#include "Group.h"
#include "GRP_GroupFct.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TEST_WeightMatrixLimited [EXTIN]\n";
      std::cerr << "\n";
      std::cerr << "EXTIN : The list of vertices (or inequalities for that matter)\n";
      return -1;
    }
    //
    using T=mpq_class;
    using Tint = mpz_class;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt,Tint>;
    //
    std::ifstream is(argv[1]);
    MyMatrix<T> EXT=ReadMatrix<T>(is);
    size_t len = EXT.rows();
    //
    const bool use_scheme1 = true;
    Tgroup GRP = LinPolytope_Automorphism<T,use_scheme1,Tgroup>(EXT);
    //
    auto test_WMat = [&](const WeightMatrixLimited<true,T>& WMatLimited) -> void {
      for (int iter=0; iter<20; iter++) {
        Face eFace(len);
        for (Tidx i=0; i<len; i++) {
          int eVal = Tidx(rand()) % 2;
          eFace[i] = eVal;
        }
        size_t hash1 = WMatLimited.get_hash(eFace);
        for (int iter2=0; iter2<20; iter2++) {
          Telt eElt = GRP.rand();
          Face fFace = OnFace(eFace, eElt);
          size_t hash2 = WMatLimited.get_hash(eFace);
          if (hash1 != hash2) {
            std::cerr << "hash1=" << hash1 << " hash2=" << hash2 << "\n";
            throw TerminalException{1};
          }
        }
      }
    };
    //
    size_t total = len * (len - 1) / 2;
    int expo = 10;
    for (int i=1; i<=expo; i++) {
      size_t max_offdiag = size_t( double(total+1) * (double(i) / double(10)) );
      WeightMatrixLimited<true,T> WMatLimited1 = GetWeightMatrixLimited(EXT, max_offdiag);
      test_WMat(WMatLimited1);
      WeightMatrixLimited<true,T> WMatLimited2(GRP, max_offdiag);
      test_WMat(WMatLimited2);
    }
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}

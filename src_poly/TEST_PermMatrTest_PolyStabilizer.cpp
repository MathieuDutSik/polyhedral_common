#include "rational.h"
#include "Permutation.h"
#include "PermutationMatrix.h"
#include "Group.h"
#include "GRP_GroupFct.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_PermMatrTest_Automorphism [MatrFile] [FaceFile]\n";
      return -1;
    }
    using T = mpq_class;
    using Tidx = int32_t;
    using Tint = mpz_class;
    const bool use_scheme1 = true;
    using Telt1   = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup1 = permutalib::Group<Telt1,Tint>;
    using Telt2   = permutalib::PermutationMatrix<Tidx,T>;
    using Tgroup2 = permutalib::Group<Telt2,Tint>;
    //
    std::cerr << "GRP_ComputeAut_ListMat_Subset_EXT : Reading input\n";
    //
    std::string MatrFile = argv[1];
    std::string FaceFile = argv[2];

    MyMatrix<T> M = ReadMatrixFile<T>(MatrFile);
    Face f = ReadFaceFile(FaceFile);
    Tgroup1 GRP1 = LinPolytope_Automorphism<T,use_scheme1,Tgroup1>(M);
    //
    std::vector<Telt2> ListPermMatr;
    for (auto & ePerm : GRP1.GeneratorsOfGroup()) {
      MyMatrix<T> eMatr = FindTransformation(M, M, ePerm);
      Telt2 elt2(ePerm.getListVal(), eMatr);
      ListPermMatr.push_back(elt2);
    }
    Telt2 elt2(M.rows(), M.cols());
    Tgroup2 GRP2(ListPermMatr, elt2);
    Tgroup2 stab=GRP2.Stabilizer_OnSets(f);
    //
    size_t pos=0;
    for (auto & eGen : stab.GeneratorsOfGroup()) {
      std::cerr << "pos=" << pos << "M=\n";
      WriteMatrix(std::cerr, eGen.getMatrix());
      pos++;
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}

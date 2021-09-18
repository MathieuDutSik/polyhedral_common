#include "Permlib_specific.h"
#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PolytopeEquiStab.h"
#include "vinberg_code.h"


int main(int argc, char* argv[])
{
  try {
    if (argc != 3 && argc != 2) {
      std::cerr << "VIN_ComputeDomain [FileI] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "VIN_ComputeDomain [FileI]\n";
      throw TerminalException{1};
    }
    using T=mpq_class;
    using Tint=mpz_class;

    std::string FileI = argv[1];
    std::ifstream is(FileI);
    //
    MyMatrix<Tint> G = ReadMatrix<Tint>(is);
    std::cerr << "We have G\n";
    size_t n_domain;
    is >> n_domain;
    std::vector<MyMatrix<Tint>> ListDomain(n_domain);
    for (size_t i=0; i<n_domain; i++) {
      MyMatrix<Tint> Dom = ReadMatrix<Tint>(is);
      ListDomain[i] = Dom;
    }
    //
    std::vector<std::vector<MyMatrix<Tint>>> ListListDomain;
    ListListDomain.push_back(ListDomain);
    //
    int TheDim = G.rows();
    for (int i=1; i<TheDim-1; i++) {
      using Tent = std::tuple<size_t,MyMatrix<Tint>,MyMatrix<Tint>,MyMatrix<Tint>>;
      auto f_stab=[&](const Tent& eEnt) -> Tgroup {
      };
      auto f_equiv=[&](const Tent& eEnt, const Tent& fEnt) -> std::option<MyMatrix<Tint>> {
        if (std::get<0>(eEnt) == std::get<0>(fEnt))
          return {};

      };
      std::vector<Tent> NewListCand;
      auto f_ent=[&](const MyMatrix<Tint>& M) -> Tent {
        MyMatrix<Tint> P = M * G;
        MyMatrix<Tint> NSP = NullspaceIntMat(P);
        MyMatrix<Tint> Gres = NSP * G * NSP.transpose();
        MyMatrix<T> Gres_T = UniversalMatrixConversion<T,Tint>(Gres);
        MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T,Tint>(Gres_T);
        MyMatrix<Tint> Spann = SHV * NSP;
        MyMatrix<Tint> Concat = ConcatenateMatrix(M, Spann);
        size_t inv = 0;
        MyMatrix<Tint> Qmat = GetQmatrix(Concat);
        return {inv, M, Spann, Qmat},
      };
      auto f_insert=[&](const Tent& fEnt) -> void {
        for (auto & eEnt : NewListCand) {
          std::option<MyMatrix<Tint>> eEquiv = f_equiv(eEnt, fEnt);
          if (eEquiv)
            return;
        }
        NewListCand.push_back(fEnt);
      };
      for (auto & eDomain : ListListDomain[i-1]) {
        Tent eEnt = f_ent(eDomain);
        Tgroup GRP = f_stab(eEnt);
        MyMatrix<T> Mred = ColumnReduction(UniversalMatrixConversion<T,Tint>(eDomain));
        vectface ListFace = DualDescriptionStandard(Mred, GRP);
        for (auto & eFace : ListFace) {
          std::vector<int> ListIdx = FaceToVector(eFace);
          MyMatrix<Tint> eDomainFace = SelectRow(eDomain, ListIdx);
          Tent fEnt = f_ent(eDomainFace);
          f_insert(fEnt);
        }
      }
      std::vector<MyMatrix<Tint>> NewListDomain;
      for (auto & eEnt : NewListCand)
        NewListDomain.push_back(std::get<1>(eEnt));
      ListListDomain.push_back(NewListDomain);
      std::cerr << "i=" << i << " |NewListDomain|=" << NewListDomain.size() << "\n";
    }

    //
    if (argc == 2) {
      Print_DataReflectionGroup(data, std::cerr);
    } else {
      std::string FileO = argv[2];
      std::ofstream os(FileO);
      Print_DataReflectionGroup(data, os);
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}

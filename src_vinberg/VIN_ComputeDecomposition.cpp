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
      struct Tent {
        MyMatrix<Tint> M;
        MyMatrix<Tint> Spann;
        MyMatrix<Tint> Qmat;
        WeightMatrix<true,std::vector<Tint>,Tidx_value>> WMat;
      };
      auto f_subset=[&](const size_t& n_row, const size_t& len) -> Face {
        Face subset(n_row);
        for (size_t i=0; i<len; i++)
          subset[i] = 1;
        return subset;
      };
      auto f_stab=[&](const Tent& eEnt) -> Tgroup {
        Tgroup GRP1 = GetStabilizerWeightMatrix_Kernel<std::vector<Tint>,Tgr,Tidx,Tidx_value>(eEnt.WMat);
        MyMatrix<Tint> Concat = ConcatenateMatrix(eEnt.M, eEnt.Spann);
        Tgroup GRP2 = LinPolytopeIntegral_Stabilizer_Method8(Concat, GRP1);
        Face subset = f_subset(Concat.rows(), eEnt.M.rows());
        return ReducedGroupAction(GRP2, subset);
      };
      auto f_equiv=[&](const Tent& eEnt, const Tent& fEnt) -> std::option<MyMatrix<Tint>> {
        MyMatrix<Tint> eConcat = ConcatenateMatrix(eEnt.M, eEnt.Spann);
        MyMatrix<Tint> fConcat = ConcatenateMatrix(fEnt.M, fEnt.Spann);
        std::vector<Tidx> eCanonicReord = GetGroupCanonicalizationVector_Kernel<T,Tgr,Tidx,Tidx_value>(eEnt.WMat).first;
        std::vector<Tidx> fCanonicReord = GetGroupCanonicalizationVector_Kernel<T,Tgr,Tidx,Tidx_value>(fEnt.WMat).first;
        // Computing the isomorphism
        using Tfield = typename overlying_field<Tint>::field_type;
        std::optional<std::pair<std::vector<Tidx>,MyMatrix<Tfield>>> IsoInfo = IsomorphismFromCanonicReord<T,Tfield,Tidx>(eConcat, fConcat, eCanonicReord, fCanonicReord);
        if (!IsoInfo)
          return {};
        Telt ePerm(IsoInfo->first);
        return LinPolytopeIntegral_Isomorphism_Method8(eConcat, fConcat, GRP1, ePerm);
      };
      std::vector<std::pair<size_t,Tent>> NewListCand;
      auto f_ent=[&](const MyMatrix<Tint>& M) -> Tent {
        MyMatrix<Tint> P = M * G;
        MyMatrix<Tint> NSP = NullspaceIntMat(P);
        MyMatrix<Tint> Gres = NSP * G * NSP.transpose();
        MyMatrix<T> Gres_T = UniversalMatrixConversion<T,Tint>(Gres);
        MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T,Tint>(Gres_T);
        MyMatrix<Tint> Spann = SHV * NSP;
        MyMatrix<Tint> Concat = ConcatenateMatrix(M, Spann);
        MyMatrix<Tint> Qmat = GetQmatrix(Concat);
        Face subset = f_subset(Concat.rows(), M.rows());
        std::vector<MyMatrix<Tint>> ListMat{Qmat, G};
        WeightMatrix<true,std::vector<Tint>,Tidx_value> WMat = GetWeightMatrix_ListMat_Subset<Tint,Tidx,Tidx_value>(EXT1, ListMat, subset);
        WMat.ReorderingSetWeight();
        return {M, Spann, Qmat, WMat},
      };
      auto f_inv=[&](const Tent& eEnt) -> size_t {
        return std::hash<WeightMatrix<true, std::vector<Tint>, Tidx_value>>()(WMat);
      };
      auto f_insert=[&](const Tent& fEnt) -> void {
        size_t e_inv = f_inv(fEnt);
        for (auto & eP : NewListCand) {
          if (eP.first == e_inv) {
            std::option<MyMatrix<Tint>> eEquiv = f_equiv(eP.second, fEnt);
            if (eEquiv)
              return;
          }
        }
        NewListCand.push_back({e_inv,fEnt});
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
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}

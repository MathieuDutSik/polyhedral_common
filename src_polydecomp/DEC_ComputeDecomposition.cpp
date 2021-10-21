#include "Permutation.h"
#include "Group.h"
#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PolytopeEquiStab.h"
#include "POLY_RecursiveDualDesc.h"
#include "POLY_Kskeletton.h"


// possible strategies for computing isomorphsim


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
    using Tgr=GraphBitset;
    using Tint=mpz_class;
    using Tidx_value = int16_t;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt,Tint>;
    //
    std::string FileI = argv[1];
    std::ifstream is(FileI);
    //
    // The polyhedral cone.
    // So far, we encode all the EXT / FAC and the group
    struct ConeDesc {
      MyMatrix<T> EXT;
      MyMatrix<Tint> EXT_i;
      MyMatrix<T> FAC;
      Face extfac_incd;
      Tgroup GRP_ext;
    };
    std::vector<ConeDesc> ListCones;
    //
    struct FaceDesc {
      size_t iCone;
      Face f;
    };
    //
    MyMatrix<Tint> G = ReadMatrix<Tint>(is);
    std::cerr << "We have G\n";
    size_t n_domain;
    is >> n_domain;
    std::vector<FaceDesc> ListDomain(n_domain);
    for (size_t i=0; i<n_domain; i++) {
      MyMatrix<T> EXT = ReadMatrix<T>(is);
      MyMatrix<Tint> EXT_i = UniversalMatrixConversion<Tint,T>(EXT);
      MyMatrix<T> FAC = ReadMatrix<T>(is);
      Tgroup GRP_ext = ReadGroup<Tgroup>(is);
      //
      int nbFac = FAC.rows();
      int nbExt = EXT.rows();
      int nbCol = EXT.cols();
      Face extfac_incd(nbFac * nbExt);
      for (int iFac=0; iFac<nbFac; iFac++) {
        for (int iExt=0; iExt<nbExt; iExt++) {
          T sum = 0;
          for (int i=0; i<nbCol; i++)
            sum += FAC(iFac,i) * EXT(iExt,i);
          if (sum == 0)
            extfac_incd[iFac * nbExt + iExt] = 1;
        }
      }
      ListCones[i] = {EXT, EXT_i, FAC, extfac_incd, GRP_ext};
      size_t len = EXT.rows();
      Face f(len);
      for (size_t i=0; i<len; i++)
        f[i] = 1;
      ListDomain[i] = {i, f};
    }
    //
    std::vector<std::vector<FaceDesc>> ListListDomain;
    ListListDomain.push_back(ListDomain);
    //
    int TheDim = G.rows();
    for (int i=1; i<TheDim-1; i++) {
      struct Tent {
        MyMatrix<Tint> M;
        MyMatrix<Tint> Spann;
        MyMatrix<Tint> Qmat;
        WeightMatrix<true,std::vector<Tint>,Tidx_value> WMat;
        FaceDesc fd;
      };
      auto f_subset=[&](const size_t& n_row, const size_t& len) -> Face {
        Face subset(n_row);
        for (size_t i=0; i<len; i++)
          subset[i] = 1;
        return subset;
      };
      auto f_stab=[&](const Tent& eEnt) -> Tgroup {
        Tgroup GRP1 = GetStabilizerWeightMatrix<std::vector<Tint>,Tgr,Tgroup,Tidx_value>(eEnt.WMat);
        MyMatrix<T> Concat_T = UniversalMatrixConversion<T,Tint>(Concatenate(eEnt.M, eEnt.Spann));
        Tgroup GRP2 = LinPolytopeIntegral_Stabilizer_Method8(Concat_T, GRP1);
        Face subset = f_subset(Concat_T.rows(), eEnt.M.rows());
        return ReducedGroupAction(GRP2, subset);
      };
      auto f_equiv=[&](const Tent& eEnt, const Tent& fEnt) -> std::optional<MyMatrix<Tint>> {
        MyMatrix<T> eConcat_T = UniversalMatrixConversion<T,Tint>(Concatenate(eEnt.M, eEnt.Spann));
        MyMatrix<T> fConcat_T = UniversalMatrixConversion<T,Tint>(Concatenate(fEnt.M, fEnt.Spann));
        std::vector<Tidx> eCanonicReord = GetGroupCanonicalizationVector_Kernel<std::vector<Tint>,Tgr,Tidx,Tidx_value>(eEnt.WMat).first;
        std::vector<Tidx> fCanonicReord = GetGroupCanonicalizationVector_Kernel<std::vector<Tint>,Tgr,Tidx,Tidx_value>(fEnt.WMat).first;
        // Computing the isomorphism
        using Tfield = typename overlying_field<Tint>::field_type;
        std::optional<std::pair<std::vector<Tidx>,MyMatrix<Tfield>>> IsoInfo = IsomorphismFromCanonicReord<T,Tfield,Tidx>(eConcat_T, fConcat_T, eCanonicReord, fCanonicReord);
        if (!IsoInfo)
          return {};
        Telt ePerm(IsoInfo->first);
        Tgroup GRP1 = GetStabilizerWeightMatrix<std::vector<Tint>,Tgr,Tgroup,Tidx_value>(eEnt.WMat);
        std::optional<MyMatrix<T>> eRes = LinPolytopeIntegral_Isomorphism_Method8(eConcat_T, fConcat_T, GRP1, ePerm);
        if (eRes)
          return UniversalMatrixConversion<Tint,T>(*eRes);
        return {};
      };
      std::vector<std::pair<size_t,Tent>> NewListCand;
      auto f_ent=[&](const FaceDesc& fd) -> Tent {
        MyMatrix<Tint> M = SelectRow(ListCones[fd.iCone].EXT_i, fd.f);
        MyMatrix<Tint> P = M * G;
        MyMatrix<Tint> NSP = NullspaceIntMat(P);
        MyMatrix<Tint> Gres = NSP * G * NSP.transpose();
        MyMatrix<T> Gres_T = UniversalMatrixConversion<T,Tint>(Gres);
        MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T,Tint>(Gres_T);
        MyMatrix<Tint> Spann = SHV * NSP;
        MyMatrix<Tint> Concat = Concatenate(M, Spann);
        MyMatrix<Tint> Qmat = GetQmatrix(Concat);
        Face subset = f_subset(Concat.rows(), M.rows());
        std::vector<MyMatrix<Tint>> ListMat{Qmat, G};
        WeightMatrix<true,std::vector<Tint>,Tidx_value> WMat = GetWeightMatrix_ListMat_Subset<Tint,Tidx,Tidx_value>(Concat, ListMat, subset);
        WMat.ReorderingSetWeight();
        return {M, Spann, Qmat, std::move(WMat), fd};
      };
      auto f_inv=[&](const Tent& eEnt) -> size_t {
        return std::hash<WeightMatrix<true, std::vector<Tint>, Tidx_value>>()(eEnt.WMat);
      };
      auto f_insert=[&](Tent&& eEnt) -> void {
        size_t e_inv = f_inv(eEnt);
        for (auto & eP : NewListCand) {
          if (eP.first == e_inv) {
            std::optional<MyMatrix<Tint>> eEquiv = f_equiv(eP.second, eEnt);
            if (eEquiv)
              return;
          }
        }
        std::pair<size_t,Tent> e_pair{e_inv,std::move(eEnt)};
        NewListCand.emplace_back(std::move(e_pair));
      };
      for (auto & eDomain : ListListDomain[i-1]) {
        Tent eEnt = f_ent(eDomain);
        size_t iCone = eDomain.iCone;
        Tgroup StabFace = ListCones[iCone].GRP_ext.Stabilizer_OnSets(eDomain.f);
        int RankFace = TheDim - i;
        vectface ListFace = SPAN_face_ExtremeRays(eDomain.f, StabFace, RankFace,
                                                  ListCones[iCone].extfac_incd, ListCones[iCone].FAC,
                                                  ListCones[iCone].EXT, ListCones[iCone].GRP_ext);
        //        Tgroup GRP = f_stab(eEnt);
        //        vectface ListFace = DualDescriptionStandard(Mred, GRP);
        for (auto & eFace : ListFace) {
          FaceDesc fdn{iCone, eFace};
          Tent fEnt = f_ent(fdn);
          f_insert(std::move(fEnt));
        }
      }
      std::vector<FaceDesc> NewListDomain;
      for (auto & eEnt : NewListCand)
        NewListDomain.push_back(eEnt.second.fd);
      std::cerr << "i=" << i << " |NewListDomain|=" << NewListDomain.size() << "\n";
      ListListDomain.emplace_back(std::move(NewListDomain));
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}

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
    using Tidx_value = int32_t;
    using Tidx = uint32_t;
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
      Tgroup GRP_fac;
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
      std::cerr << "i=" << i << " / " << n_domain << "\n";
      MyMatrix<T> EXT = ReadMatrix<T>(is);
      //      std::cerr << "EXT=\n";
      //      WriteMatrix(std::cerr, EXT);
      MyMatrix<Tint> EXT_i = UniversalMatrixConversion<Tint,T>(EXT);
      MyMatrix<T> FAC = ReadMatrix<T>(is);
      //      std::cerr << "FAC=\n";
      //      WriteMatrix(std::cerr, FAC);
      Tgroup GRP_ext = ReadGroup<Tgroup>(is);
      Tgroup GRP_fac = ReadGroup<Tgroup>(is);
      std::cerr << "|GRP_ext|=" << GRP_ext.size() << " |GRP_fac|=" << GRP_fac.size() << "\n";
      //
      Face extfac_incd = Compute_extfac_incd(FAC, EXT);
      ConeDesc eCone{EXT, EXT_i, FAC, extfac_incd, GRP_ext, GRP_fac};
      ListCones.push_back(eCone);
      size_t len = FAC.rows();
      Face f(len);
      //      for (size_t i=0; i<len; i++)
      //        f[i] = 1;
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
        std::cerr << "f_ent, step 1\n";
        int nbFac = ListCones[fd.iCone].FAC.rows();
        int nbExt = ListCones[fd.iCone].EXT.rows();
        Face face_ext = Compute_faceEXT_from_faceFAC(ListCones[fd.iCone].extfac_incd, nbFac, nbExt, fd.f);
        MyMatrix<Tint> M = SelectRow(ListCones[fd.iCone].EXT_i, face_ext);
        std::cerr << "f_ent, step 2\n";
        //        std::cerr << "M=\n";
        //        WriteMatrix(std::cerr, M);
        MyMatrix<Tint> P = M * G;
        std::cerr << "f_ent, step 3\n";
        //        std::cerr << "P=\n";
        //        WriteMatrix(std::cerr, P);
        MyMatrix<Tint> Spann;
        std::cerr << "f_ent, step 4\n";
        MyMatrix<Tint> Concat = M;
        std::cerr << "f_ent, step 5\n";
        MyMatrix<Tint> NSP = NullspaceIntMat(TransposedMat(P));
        std::cerr << "f_ent, step 6 |NSP|=" << NSP.rows() << " / " << NSP.cols() << "\n";
        if (NSP.rows() > 0) {
          std::cerr << "f_ent, step 7\n";
          MyMatrix<Tint> Gres = NSP * G * NSP.transpose();
          std::cerr << "f_ent, step 8\n";
          MyMatrix<T> Gres_T = UniversalMatrixConversion<T,Tint>(Gres);
          std::cerr << "f_ent, step 9\n";
          std::cerr << "Gres_T=\n";
          WriteMatrix(std::cerr, Gres_T);
          MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T,Tint>(Gres_T);
          std::cerr << "f_ent, step 10\n";
          Spann = SHV * NSP;
          std::cerr << "f_ent, step 11\n";
          Concat = Concatenate(M, Spann);
          std::cerr << "f_ent, step 12\n";
        }
        std::cerr << "f_ent, step 13\n";
        MyMatrix<Tint> Qmat = GetQmatrix(Concat);
        std::cerr << "f_ent, step 14\n";
        Face subset = f_subset(Concat.rows(), M.rows());
        std::cerr << "f_ent, step 15\n";
        std::vector<MyMatrix<Tint>> ListMat{Qmat, G};
        std::cerr << "f_ent, step 16\n";
        WeightMatrix<true,std::vector<Tint>,Tidx_value> WMat = GetWeightMatrix_ListMat_Subset<Tint,Tidx,Tidx_value>(Concat, ListMat, subset);
        std::cerr << "f_ent, step 17\n";
        WMat.ReorderingSetWeight();
        std::cerr << "f_ent, step 18\n";
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
        std::cerr << "Now |NewListCand|=" << NewListCand.size() << "\n";
      };
      for (auto & eDomain : ListListDomain[i-1]) {
        std::cerr << "eDomain, step 1\n";
        Tent eEnt = f_ent(eDomain);
        std::cerr << "eDomain, step 2\n";
        size_t iCone = eDomain.iCone;
        std::cerr << "eDomain, step 3\n";
        Tgroup StabFace = ListCones[iCone].GRP_fac.Stabilizer_OnSets(eDomain.f);
        std::cerr << "eDomain, step 4\n";
        int RankFace = i;
        std::cerr << "eDomain, step 5\n";
        //        vectface ListFace = SPAN_face_ExtremeRays(eDomain.f, StabFace, RankFace,
        //                                                  ListCones[iCone].extfac_incd, ListCones[iCone].FAC,
        //                                                  ListCones[iCone].EXT, ListCones[iCone].GRP_fac);
        // For the downward road from high dimension to lower dimension, we have to work with FAC as input
        // 
        vectface ListFace = SPAN_face_ExtremeRays(eDomain.f, StabFace, RankFace,
                                                  ListCones[iCone].extfac_incd, ListCones[iCone].FAC,
                                                  ListCones[iCone].EXT, ListCones[iCone].GRP_fac);
        std::cerr << "eDomain, step 6 |ListFace|=" << ListFace.size() << "\n";
        //        Tgroup GRP = f_stab(eEnt);
        //        vectface ListFace = DualDescriptionStandard(Mred, GRP);
        for (auto & eFace : ListFace) {
          FaceDesc fdn{iCone, eFace};
          Tent fEnt = f_ent(fdn);
          f_insert(std::move(fEnt));
        }
        std::cerr << "eDomain, step 7\n";
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

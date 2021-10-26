#include "Permutation.h"
#include "Group.h"
#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PolytopeEquiStab.h"
#include "POLY_RecursiveDualDesc.h"
#include "POLY_Kskeletton.h"


// possible strategies for computing isomorphsim

template<typename T, typename Tint, typename Tgroup>
struct ConeDesc {
  MyMatrix<T> EXT;
  MyMatrix<Tint> EXT_i;
  MyMatrix<T> FAC;
  Face extfac_incd;
  Tgroup GRP_ext;
  Tgroup GRP_fac;
};

struct FaceDesc {
  size_t iCone;
  Face f;
};


template<typename Tint, typename Tidx_value>
struct Tent {
  MyMatrix<Tint> M;
  MyMatrix<Tint> Spann;
  MyMatrix<Tint> Qmat;
  WeightMatrix<true,std::vector<Tint>,Tidx_value> WMat;
  FaceDesc fd;
};


Face f_subset(const size_t& n_row, const size_t& len)
{
  Face subset(n_row);
  for (size_t i=0; i<len; i++)
    subset[i] = 1;
  return subset;
}


template<typename Tgroup, typename Tint>
struct stab_info {
  Tgroup GRPfull;
  Tgroup GRPres;
  std::vector<std::pair<typename Tgroup::Telt,MyMatrix<Tint>>> ListGenMat;
};

template<typename T, typename Tint, typename Tidx_value, typename Tgroup>
stab_info<Tgroup,Tint> f_stab(const Tent<T,Tint,Tgroup>& eEnt)
{
  using Telt=typename Tgroup::Telt;
  Tgroup GRP1 = GetStabilizerWeightMatrix<std::vector<Tint>,Tgr,Tgroup,Tidx_value>(eEnt.WMat);
  MyMatrix<T> Concat_T = UniversalMatrixConversion<T,Tint>(Concatenate(eEnt.M, eEnt.Spann));
  Tgroup GRPfull = LinPolytopeIntegral_Stabilizer_Method8(Concat_T, GRP1);
  Face subset = f_subset(Concat_T.rows(), eEnt.M.rows());
  std::vector<Telt> LGenRed;
  std::vector<std::pair<typename Tgroup::Telt,MyMatrix<Tint>>> ListGenMat;
  for (auto & eGen : GRPfull.GeneratorsOfGroup()) {
    Telt eGenRed = ReduceElementAction(eGen, subset);
    LGenRed.push_back(eGenRed);
    MyMatrix<T> eMat_T = FindTransformation(Concat_T, Concat_T, eGen);
    MyMatrix<Tint> eMat = UniversalMatrixConversion<Tint,T>(eMat_T);
    ListGenMat.push_back({eGenRed, eMat});
  }
  Tgroup GRPres(LGenRed, eEnt.M.rows());
  return {GRPfull, GRPres, ListGenMat};
}

template<typename T, typename Tint, typename Tidx_value, typename Tgroup>
std::optional<MyMatrix<Tint>> f_equiv(const Tent<T,Tint,Tgroup>& eEnt, const Tent<T,Tint,Tgroup>& fEnt)
{
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
}


Tent<T,Tint,Tgroup> f_ent(std::vector<ConeDesc<T,Tint,Tgroup>> const& ListCones, const FaceDesc& fd)
{
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
    MyMatrix<Tint> Gres = - NSP * G * NSP.transpose();
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
}

size_t f_inv(const Tent<T,Tint,Tgroup>& eEnt)
{
  return std::hash<WeightMatrix<true, std::vector<Tint>, Tidx_value>>()(eEnt.WMat);
}


template<typename T, typename Tint, typename Tgroup>
std::vector<std::vector<FaceDesc>> Compute_ListListDomain_strategy2(std::vector<ConeDesc<T,Tint,Tgroup>> const& ListCones, MyMatrix<Tint> const& G, int TheLev)
{
  std::vector<FaceDesc> ListDomain;
  for (size_t i=0; i<ListCones.size(); i++) {
    size_t len = ListCones[i].FAC.rows();
    Face f(len);
    FaceDesc fd = {i, f};
    ListDomains.push_back(fd);
  }
  std::vector<std::vector<FaceDesc>> ListListDomain;
  ListListDomain.push_back(ListDomain);
  //
  for (int i=1; i<TheLev; i++) {
    std::vector<std::pair<size_t,Tent>> NewListCand;
    auto f_insert=[&](Tent&& eEnt) -> void {
      std::cerr << "Beginning of f_insert\n";
      size_t e_inv = f_inv(eEnt);
      std::cerr << "e_inv=" << e_inv << "\n";
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
      Tent eEnt = f_ent(ListCones, eDomain);
      std::cerr << "eDomain, step 2\n";
      size_t iCone = eDomain.iCone;
      std::cerr << "eDomain, step 3\n";
      Tgroup StabFace = ListCones[iCone].GRP_fac.Stabilizer_OnSets(eDomain.f);
      std::cerr << "eDomain, step 4\n";
      int RankFace = i - 1;
      std::cerr << "eDomain, step 5\n";
      vectface ListFace = SPAN_face_ExtremeRays(eDomain.f, StabFace, RankFace,
                                                ListCones[iCone].extfac_incd, ListCones[iCone].FAC,
                                                ListCones[iCone].EXT, ListCones[iCone].GRP_fac);
      std::cerr << "eDomain, step 6 |ListFace|=" << ListFace.size() << "\n";
      //        Tgroup GRP = f_stab(eEnt);
      //        vectface ListFace = DualDescriptionStandard(Mred, GRP);
      for (auto & eFace : ListFace) {
        FaceDesc fdn{iCone, eFace};
        Tent fEnt = f_ent(ListCones, fdn);
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
  return ListListDomain;
}

template<typename Tint>
struct sing_adj {
  size_t jCone;
  Face f;
  MyMatrix<Tint> eMat;
};

template<typename T, typename Tint, typename Tgroup>
std::vector<std::vector<sing_adj>> compute_adjacency_structure(std::vector<ConeDesc<T,Tint,Tgroup>> const& ListCones, MyMatrix<Tint> const& G)
{
  std::vector<std::vector<sing_adj>> adjacency_information;
  size_t n_domain = ListCones.size();
  struct ent_info {
    size_t i_domain;
    size_t i_adj;
    Face f;
    Tent<T,Tint,Tgroup> eEnt;
    size_t hash;
  };
  std::vector<ent_info> l_ent_info;
  std::vector<size_t> l_n_orb_adj;
  for (size_t i_domain=0; i_domain<n_domain; i_domain++) {
    size_t n_fac = ListCones[i_domain].FAC.rows();
    size_t n_ext = ListCones[i_domain].EXT.rows();
    vectface vf = DecomposeOrbitPoint_Full(ListCones[i_domain].GRP_fac);
    size_t i_adj = 0;
    for (auto & eOrb : vf) {
      boost::dynamic_bitset<>::size_type MinVal=eOrb.find_first();
      size_t i_fac = MinVal;
      Face f_ext(n_ext);
      for (size_t i_ext=0; i_ext<n_ext; i_ext++)
        if (ListCone[i_domain].extfac_incd[iFac * n_ext + i_ext] == 1)
          f_ext[i_ext] = 1;
      f_fac(n_fac);
      f_fac[i_fac] = 1;
      FaceDesc fd{i_domain, f_fac};
      Tent<T,Tint,Tgroup> eEnt = f_ent(ListCones, fd);
      size_t hash = f_inv(eEnt);
      ent_info e_ent_info{i_domain, i_adj, f_ext, eEnt, hash};
      l_ent_info.push_back(e_ent_info);
      i_adj++;
    }
    l_n_orb_adj.push_back(i_adj);
  }
  auto get_ent_info=[&](size_t const& i_domain, size_t const& i_adj) -> const ent_info& {
    for (auto & e_ent : l_ent_info)
      if (e_ent.i_domain == i_domain && e_ent.i_adj == i_adj)
        return e_ent;
    std::cerr << "Failed to find the matching entry\n";
    throw TerminalException{1};
  };
  auto get_reverting_transformation=[&](const Tent<T,Tint,Tgroup> &eEnt) -> std::optional<MyMatrix<Tint>> {
    stab_info<Tgroup,Tint> e_stab_info = f_stab(eEnt);
    MyVector<Tint> eSpann = GetMatrixRow(eEnt.Spann, 0);
    for (auto& ePair : e_stab_info.ListGenMat) {
      MyVector<Tint> eSpannImg = eSpann * ePair.second;
      if (eSpannImg == - eSpann) {
        return ePair.second;
      }
      if (eSpannImg != eSpann) {
        std::cerr << "Some inconsistency in the matrix transformation\n";
        throw TerminalException{1};
      }
    }
    return {};
  };
  auto get_mapped=[&](const ent_info& e_ent) -> sing_adj {
    for (auto & f_ent : l_ent_info) {
      if (f_ent.hash == e_ent.hash && (f_ent.i_domain != e_ent.i_domain || f_ent.i_adj != e_ent.i_adj)) {
        std::optional<MyMatrix<Tint>> e_equiv = f_equiv(f_ent.eEnt, e_ent.eEnt);
        if (e_equiv) {
          return {f_ent.i_domain, e_ent.f, *e_equiv};
        }
      }
    }
    std::cerr << "Failed to find a matching entry\n";
    throw TerminalException{1};
  };
  auto get_sing_adj=[&](size_t const& i_domain, size_t const& i_adj) -> sing_adj {
    const ent_info& e_ent = get_ent_info(i_domain, i_adj);
    std::optional<MyMatrix<Tint>> trans_opt = get_reverting_transformation(e_ent.eEnt);
    if (trans_opt) {
      return {i_domain, e_ent.f, *trans_opt};
    }
    return get_mapped(e_ent);
  };

  std::vector<std::vector<sing_adj>> ll_adj_struct;
  for (size_t i_domain=0; i_domain<n_domain; i_domain++) {
    std::vector<sing_adj> l_adj_struct;
    for (size_t i_adj=0; i_adj<l_n_orb_adj[i_domain]; i_adj++)
      l_adj_struct.push_back(get_sing_adj(i_domain, i_adj));
    ll_adj_struct.push_back(l_adj);
  }
  return ll_adj_struct;
}


// A face of the cellular complex is determined by all the ways in which 
struct ent_face {
  size_t iCone;
  size_t iorb_facet;
  Telt ePerm;
};


std::optional<MyMatrix<Tint>> test_equiv_ent_face(ent_face const& ef1, ent_face const& ef2)
{
  
}

std::vector<ent_face> get_spanning_list_ent_face(std::vector<ConeDesc<T,Tint,Tgroup>> const& ListCones, const ent_face& ef)
{
  std::vector<ent_face> l_ent_face;
  std::vector<uint8_t> l_status;
  auto f_insert=[&](ent_face& ef_A) -> void {
    for (auto & ef_B : l_ent_face) {
      
    }
    l_ent_face.push_back(ef_A);
    l_status.push_back(0);
  };
  while(true) {
    bool is_finished=true;
    for (size_t i=0; i<l_ent_face.size(); i++) {
      if (l_status[i] == 0) {
        l_status[i] = 1;
        is_finished = false;
      }
    }
    if (is_finished)
      break;
  }
  return l_ent_face;
}




  
int main(int argc, char* argv[])
{
  try {
    if (argc != 4 && argc != 3) {
      std::cerr << "VIN_ComputeDomain [opt] [FileI] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "VIN_ComputeDomain [opt] [FileI]\n";
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
    std::string opt = argv[1];
    std::string FileI = argv[2];
    std::ifstream is(FileI);
    //
    // The polyhedral cone.
    // So far, we encode all the EXT / FAC and the group
    std::vector<ConeDesc<T,Tint,Tgroup>> ListCones;
    //
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
    }
    //
    std::vector<std::vector<FaceDesc>> ListListDomain;
    if (opt == "strategy2") {
      ListListDomain = Compute_ListListDomain_strategy2(ListCones, G, TheLev);
    }
    if (opt == "strategy1") {
      std::vector<std::vector<sing_adj>> ll_sing_adj = compute_adjacency_structure(ListCones, G);
      
    }
    
    if (
    
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}

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
  Face facext_incd;
  Tgroup GRP_ext;
  Tgroup GRP_fac;
};

struct FaceDesc {
  size_t iCone;
  Face f_fac;
};


template<typename T, typename Tint, typename Tidx_value>
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

template<typename T, typename Tint, typename Tgroup, typename Tidx_value>
stab_info<Tgroup,Tint> f_stab(const Tent<T,Tint,Tidx_value>& eEnt)
{
  using Telt=typename Tgroup::Telt;
  using Tgr=GraphBitset;
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

template<typename T, typename Tint, typename Tgroup, typename Tidx_value>
std::optional<MyMatrix<Tint>> f_equiv(const Tent<T,Tint,Tidx_value>& eEnt, const Tent<T,Tint,Tidx_value>& fEnt)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tgr = GraphBitset;
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


template<typename T, typename Tint, typename Tgroup, typename Tidx_value>
Tent<T,Tint,Tidx_value> f_ent(std::vector<ConeDesc<T,Tint,Tgroup>> const& ListCones, const MyMatrix<Tint>& G, const FaceDesc& fd)
{
  using Tidx = typename Tgroup::Telt::Tidx;
  std::cerr << "f_ent, step 1\n";
  int nbFac = ListCones[fd.iCone].FAC.rows();
  int nbExt = ListCones[fd.iCone].EXT.rows();
  Face face_ext = Compute_faceEXT_from_faceFAC(ListCones[fd.iCone].extfac_incd, nbFac, nbExt, fd.f_fac);
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

template<typename T, typename Tint, typename Tgroup, typename Tidx_value>
size_t f_inv(const Tent<T,Tint,Tidx_value>& eEnt)
{
  static_assert(std::is_integral<Tidx_value>::value, "Tidx_value should be integral");
  const WeightMatrix<true,std::vector<Tint>,Tidx_value>& WMat = eEnt.WMat;
  return std::hash<WeightMatrix<true, std::vector<Tint>, Tidx_value>>()(WMat);
}


template<typename T, typename Tint, typename Tgroup, typename Tidx_value>
std::vector<std::vector<FaceDesc>> Compute_ListListDomain_strategy2(std::vector<ConeDesc<T,Tint,Tgroup>> const& ListCones, MyMatrix<Tint> const& G, int TheLev)
{
  static_assert(std::is_integral<Tidx_value>::value, "Tidx_value should be integral");
  std::vector<FaceDesc> ListDomain;
  for (size_t i=0; i<ListCones.size(); i++) {
    size_t len = ListCones[i].FAC.rows();
    Face f_fac(len);
    FaceDesc fd = {i, f_fac};
    ListDomain.push_back(fd);
  }
  std::vector<std::vector<FaceDesc>> ListListDomain;
  ListListDomain.push_back(ListDomain);
  //
  for (int i=1; i<TheLev; i++) {
    std::vector<std::pair<size_t,Tent<T,Tint,Tidx_value>>> NewListCand;
    auto f_insert=[&](Tent<T,Tint,Tidx_value>&& eEnt) -> void {
      std::cerr << "Beginning of f_insert\n";
      size_t e_inv = f_inv<T,Tint,Tgroup,Tidx_value>(eEnt);
      std::cerr << "e_inv=" << e_inv << "\n";
      for (auto & eP : NewListCand) {
        if (eP.first == e_inv) {
          std::optional<MyMatrix<Tint>> eEquiv = f_equiv<T,Tint,Tgroup,Tidx_value>(eP.second, eEnt);
          if (eEquiv)
            return;
        }
      }
      std::pair<size_t,Tent<T,Tint,Tidx_value>> e_pair{e_inv,std::move(eEnt)};
      NewListCand.emplace_back(std::move(e_pair));
      std::cerr << "Now |NewListCand|=" << NewListCand.size() << "\n";
    };
    for (auto & eDomain : ListListDomain[i-1]) {
      std::cerr << "eDomain, step 1\n";
      Tent<T,Tint,Tidx_value> eEnt = f_ent<T,Tint,Tgroup,Tidx_value>(ListCones, G, eDomain);
      std::cerr << "eDomain, step 2\n";
      size_t iCone = eDomain.iCone;
      std::cerr << "eDomain, step 3\n";
      Tgroup StabFace_fac = ListCones[iCone].GRP_fac.Stabilizer_OnSets(eDomain.f_fac);
      std::cerr << "eDomain, step 4\n";
      int RankFace = i - 1;
      std::cerr << "eDomain, step 5\n";
      vectface ListFace = SPAN_face_ExtremeRays(eDomain.f_fac, StabFace_fac, RankFace,
                                                ListCones[iCone].extfac_incd, ListCones[iCone].FAC,
                                                ListCones[iCone].EXT);
      std::cerr << "eDomain, step 6 |ListFace|=" << ListFace.size() << "\n";
      //        Tgroup GRP = f_stab(eEnt);
      //        vectface ListFace = DualDescriptionStandard(Mred, GRP);
      for (auto & eFace_fac : ListFace) {
        FaceDesc fdn{iCone, eFace_fac};
        Tent<T,Tint,Tidx_value> fEnt = f_ent<T,Tint,Tgroup,Tidx_value>(ListCones, G, fdn);
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
  Face f_ext;
  MyMatrix<Tint> eMat;
};

template<typename T, typename Tint, typename Tgroup, typename Tidx_value>
std::vector<std::vector<sing_adj<Tint>>> compute_adjacency_structure(std::vector<ConeDesc<T,Tint,Tgroup>> const& ListCones, MyMatrix<Tint> const& G)
{
  static_assert(std::is_integral<Tidx_value>::value, "Tidx_value should be integral");
  std::vector<std::vector<sing_adj<Tint>>> adjacency_information;
  size_t n_domain = ListCones.size();
  struct ent_info {
    size_t i_domain;
    size_t i_adj;
    Face f_ext;
    Tent<T,Tint,Tidx_value> eEnt;
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
        if (ListCones[i_domain].extfac_incd[i_fac * n_ext + i_ext] == 1)
          f_ext[i_ext] = 1;
      Face f_fac(n_fac);
      f_fac[i_fac] = 1;
      FaceDesc fd{i_domain, f_fac};
      Tent<T,Tint,Tidx_value> eEnt = f_ent<T,Tint,Tgroup,Tidx_value>(ListCones, G, fd);
      size_t hash = f_inv<T,Tint,Tgroup,Tidx_value>(eEnt);
      ent_info e_ent_info{i_domain, i_adj, f_ext, std::move(eEnt), hash};
      l_ent_info.emplace_back(std::move(e_ent_info));
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
  auto get_reverting_transformation=[&](const Tent<T,Tint,Tidx_value> &eEnt) -> std::optional<MyMatrix<Tint>> {
    stab_info<Tgroup,Tint> e_stab_info = f_stab<T,Tint,Tgroup,Tidx_value>(eEnt);
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
  auto get_mapped=[&](const ent_info& a_ent) -> sing_adj<Tint> {
    for (auto & b_ent : l_ent_info) {
      if (b_ent.hash == a_ent.hash && (a_ent.i_domain != a_ent.i_domain || a_ent.i_adj != b_ent.i_adj)) {
        std::optional<MyMatrix<Tint>> e_equiv = f_equiv<T,Tint,Tgroup,Tidx_value>(b_ent.eEnt, a_ent.eEnt);
        if (e_equiv) {
          return {b_ent.i_domain, a_ent.f_ext, *e_equiv};
        }
      }
    }
    std::cerr << "Failed to find a matching entry\n";
    throw TerminalException{1};
  };
  auto get_sing_adj=[&](size_t const& i_domain, size_t const& i_adj) -> sing_adj<Tint> {
    const ent_info& e_ent = get_ent_info(i_domain, i_adj);
    std::optional<MyMatrix<Tint>> trans_opt = get_reverting_transformation(e_ent.eEnt);
    if (trans_opt) {
      return {i_domain, e_ent.f_ext, *trans_opt};
    }
    return get_mapped(e_ent);
  };
  std::vector<std::vector<sing_adj<Tint>>> ll_adj_struct;
  for (size_t i_domain=0; i_domain<n_domain; i_domain++) {
    std::vector<sing_adj<Tint>> l_adj_struct;
    for (size_t i_adj=0; i_adj<l_n_orb_adj[i_domain]; i_adj++)
      l_adj_struct.push_back(get_sing_adj(i_domain, i_adj));
    ll_adj_struct.push_back(l_adj_struct);
  }
  return ll_adj_struct;
}


// A face of the cellular complex is determined by all the ways in which 
template<typename Tint>
struct ent_face {
  size_t iCone;
  Face f_ext; // The subset itself 
  MyMatrix<Tint> eMat;
};


template<typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>> test_equiv_ent_face(std::vector<ConeDesc<T,Tint,Tgroup>> const& ListCones, ent_face<Tint> const& ef1, ent_face<Tint> const& ef2)
{
  using Telt = typename Tgroup::Telt;
  if (ef1.iCone != ef2.iCone)
    return {};
  size_t iC = ef1.iCone;
  const ConeDesc<T,Tint,Tgroup>& eC = ListCones[iC];
  std::pair<bool,Telt> test = eC.GRP_ext.RepresentativeAction_OnSets(ef1.f_ext, ef2.f_ext);
  if (!test.first)
    return {};
  MyMatrix<T> eMat_T = FindTransformation(eC.EXT, eC.EXT, test.second);
  MyMatrix<Tint> eMat = UniversalMatrixConversion<Tint,T>(eMat_T);
  return Inverse(ef1.eMat) * eMat * ef2.eMat;
}

/*
  Generate the list of entries in the face and the list of stabilizer generators
 */
template<typename T, typename Tint, typename Tgroup>
std::pair<std::vector<ent_face<Tint>>,std::vector<MyMatrix<Tint>>> get_spanning_list_ent_face(std::vector<ConeDesc<T,Tint,Tgroup>> const& ListCones, std::vector<std::vector<sing_adj<Tint>>> const& ll_sing_adj, const ent_face<Tint>& ef_input)
{
  using Telt=typename Tgroup::Telt;
  std::vector<MyMatrix<Tint>> ListMatrGen;
  std::set<MyVector<Tint>> EXT;
  size_t dim;
  for (auto & ePt : FaceToVector(ef_input.f_ext)) {
    MyVector<Tint> V = GetMatrixRow(ListCones[ef_input.iCone].EXT_i, ePt);
    dim = V.size();
    // In GAP Vimg = V A and in transpose we get Vimg^T = A^T V^T
    MyVector<Tint> Vimg = ef_input.eMat.transpose() * V;
    EXT.insert(Vimg);
  }
  auto f_insert_generator=[&](const MyMatrix<Tint>& eMatrGen) -> void {
    ListMatrGen.push_back(eMatrGen);
    for (auto & eV : EXT) {
      MyVector<Tint> Vimg = eMatrGen.transpose() * eV;
      if (EXT.count(Vimg) != 1) {
        std::cerr << "Error: The generator does not preserve the face globally\n";
        throw TerminalException{1};
      }
    }
  };
  std::vector<ent_face<Tint>> l_ent_face;
  auto f_insert=[&](const ent_face<Tint>& ef_A) -> void {
    for (auto & ef_B : l_ent_face) {
      std::optional<MyMatrix<Tint>> equiv_opt = test_equiv_ent_face(ListCones, ef_A, ef_B);
      if (equiv_opt) {
        f_insert_generator(*equiv_opt);
        return;
      }
    }
    l_ent_face.push_back(ef_A);
    const ConeDesc<T,Tint,Tgroup>& eC = ListCones[ef_A.iCone];
    Tgroup stab = eC.GRP_ext.Stabilizer_OnSets(ef_A.f_ext);
    for (auto & eGen : stab.GeneratorsOfGroup()) {
      MyMatrix<T> eMatGen_T = FindTransformation(eC.EXT, eC.EXT, eGen);
      MyMatrix<Tint> eMatGen = UniversalMatrixConversion<Tint,T>(eMatGen_T);
      f_insert_generator(eMatGen);
    }
  };
  f_insert(ef_input);
  size_t curr_pos = 0;
  while(true) {
    size_t len = l_ent_face.size();
    if (curr_pos == len)
      break;
    for (size_t i=curr_pos; i<len; i++) {
      const ent_face<Tint>& ef = l_ent_face[i];
      const ConeDesc<T,Tint,Tgroup>& eC = ListCones[ef.iCone];
      for (auto & e_sing_adj : ll_sing_adj[ef.iCone]) {
        std::vector<std::pair<Face, Telt>> l_pair = FindContainingOrbit(eC.GRP_ext, e_sing_adj.f_ext, ef.f_ext);
        for (auto & e_pair : l_pair) {
          MyMatrix<T> eMat1_T = FindTransformation(eC.EXT, eC.EXT, e_pair.second);
          MyMatrix<Tint> eMat1 = UniversalMatrixConversion<Tint,T>(eMat1_T);
          size_t jCone = e_sing_adj.jCone;
          MyMatrix<Tint> eMatAdj = e_sing_adj.eMat * ef.eMat * eMat1; // Needs to be cleaned up
          MyMatrix<Tint> EXTimg = eMatAdj.transpose() * eC.EXT_i;
          MyMatrix<Tint> VectorContain(1,dim);
          ContainerMatrix<Tint> Cont(EXTimg, VectorContain);
          Face faceNew(EXTimg.rows());
          for (auto & e_line : EXT) {
            for (size_t i_col=0; i_col<dim; i_col++)
              VectorContain(0,i_col) = e_line(i_col);
            std::pair<bool, size_t> epair = Cont.GetIdx();
            if (!epair.first) {
              std::cerr << "The vector is not in the image. Clear bug\n";
              throw TerminalException{1};
            }
            faceNew[epair.second] = 1;
          }
          ent_face<Tint> efNew{jCone, faceNew, eMatAdj};
          f_insert(efNew);
        }
      }
    }
  }
  return {l_ent_face, ListMatrGen};
}



template<typename T, typename Tint, typename Tgroup, typename Tidx_value>
std::vector<std::vector<FaceDesc>> Compute_ListListDomain_strategy1(std::vector<ConeDesc<T,Tint,Tgroup>> const& ListCones, MyMatrix<Tint> const& G, const std::vector<std::vector<sing_adj<Tint>>> & ll_sing_adj, int TheLev)
{
  std::vector<std::vector<FaceDesc>> ListListDomain;
  size_t n_col = G.rows();
  using Tface = std::pair<std::vector<ent_face<Tint>>,std::vector<MyMatrix<Tint>>>;
  for (int iLev=1; iLev<=TheLev; iLev++) {
    std::vector<Tface> list_face;
    auto f_insert=[&](const ent_face<Tint>& ef_A) -> void {
      for (auto & eOrbit : list_face) {
        for (auto & ef_B : eOrbit.first) {
          std::optional<MyMatrix<Tint>> equiv_opt = test_equiv_ent_face(ListCones, ef_A, ef_B);
          if (equiv_opt)
            return;
        }
      }
      list_face.push_back(get_spanning_list_ent_face(ListCones, ll_sing_adj, ef_A));
    };
    if (iLev == 1) {
      for (size_t iCone=0; iCone<ListCones.size(); iCone++) {
        vectface vf = DecomposeOrbitPoint_Full(ListCones[iCone].GRP_ext);
        size_t nExt = ListCones[iCone].EXT.rows();
        for (auto & eOrb : vf) {
          Face f_ext(nExt);
          boost::dynamic_bitset<>::size_type eVal=eOrb.find_first();
          size_t iExt = eVal;
          f_ext[iExt] = 1;
          ent_face<Tint> e_ent{iCone, f_ext, IdentityMat<Tint>(n_col)};
          f_insert(e_ent);
        }
      }
    } else {
      for (auto & eOrbit : ListListDomain[iLev - 2]) {
        const ConeDesc<T,Tint,Tgroup>& eC = ListCones[eOrbit.iCone];
        Face f_ext = Compute_faceEXT_from_faceFAC(eC.extfac_incd, eC.FAC.rows(), eC.EXT.rows(), eOrbit.f_fac);
        Tgroup StabFace_ext = eC.GRP_ext.Stabilizer_OnSets(f_ext);
        int RankFace_ext = iLev - 1;
        vectface vf = SPAN_face_ExtremeRays(f_ext, StabFace_ext, RankFace_ext,
                                            eC.facext_incd, eC.EXT, eC.FAC);
        for (auto & f_ext_new : vf) {
          ent_face<Tint> e_ent{eOrbit.iCone, f_ext_new, IdentityMat<Tint>(n_col)};
          f_insert(e_ent);
        }
      }
    }
    std::vector<FaceDesc> ListDomain;
    for (auto & eOrbit : list_face) {
      const ent_face<Tint>& e_ent = eOrbit.first[0];
      const ConeDesc<T,Tint,Tgroup>& eC = ListCones[e_ent.iCone];
      Face f_fac = Compute_faceFAC_from_faceEXT(eC.extfac_incd, eC.FAC.rows(), eC.EXT.rows(), e_ent.f_ext);
      FaceDesc fd{e_ent.iCone, f_fac};
      ListDomain.push_back(fd);
    }
    ListListDomain.push_back(ListDomain);
  }
  return ListListDomain;
}




int main(int argc, char* argv[])
{
  try {
    if (argc != 5 && argc != 4) {
      std::cerr << "VIN_ComputeDomain [opt] [TheLev] [FileI] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "VIN_ComputeDomain [opt] [TheLev] [FileI]\n";
      throw TerminalException{1};
    }
    using T=mpq_class;
    using Tint=mpz_class;
    using Tidx_value = int32_t;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt,Tint>;
    //
    std::string opt = argv[1];
    //
    std::string FileI = argv[2];
    std::ifstream is(FileI);
    //
    int TheLev;
    sscanf(argv[3], "%d", &TheLev);
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
      Face facext_incd = Compute_facext(extfac_incd, FAC.rows(), EXT.rows());
      ConeDesc<T,Tint,Tgroup> eCone{EXT, EXT_i, FAC, extfac_incd, facext_incd, GRP_ext, GRP_fac};
      ListCones.push_back(eCone);
    }
    //
    std::vector<std::vector<FaceDesc>> ListListDomain;
    if (opt == "strategy2_from_high") {
      ListListDomain = Compute_ListListDomain_strategy2<T,Tint,Tgroup,Tidx_value>(ListCones, G, TheLev);
    }
    if (opt == "strategy1_fom_low"|| opt == "strategy1_fom_high") {
      std::vector<std::vector<sing_adj<Tint>>> ll_sing_adj = compute_adjacency_structure<T,Tint,Tgroup,Tidx_value>(ListCones, G);
      if (opt == "strategy1_fom_high") {
        std::cerr << "This needs to be implemented as it does not appear needed right now\n";
        // Implementation should use above code
        throw TerminalException{1};
      }
      if (opt == "strategy1_fom_low") {
        ListListDomain = Compute_ListListDomain_strategy1<T,Tint,Tgroup,Tidx_value>(ListCones, G, ll_sing_adj, TheLev);
      }
    }
    //
    auto do_print=[&](std::ostream& os) -> void {
      size_t n_lev = ListListDomain.size();
      os << "Number of levels = " << n_lev << "\n";
      for (size_t i_lev=0; i_lev<n_lev; i_lev++) {
        size_t n_orbit = ListListDomain[i_lev].size();
        os << "Number of orbits at level " << i_lev << " = " << n_orbit << "\n";
        for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++) {
          const FaceDesc& fd = ListListDomain[i_lev][i_orbit];
          const ConeDesc<T,Tint,Tgroup>& eC = ListCones[fd.iCone];
          Face f_ext = Compute_faceEXT_from_faceFAC(eC.extfac_incd, eC.FAC.rows(), eC.EXT.rows(), fd.f_fac);
          //
          MyMatrix<Tint> EXT_sel = SelectRow(eC.EXT_i, f_ext);
          os << "i_orbit=" << i_orbit << " |EXT|=" << EXT_sel.rows() << "\n";
          WriteMatrix(os, EXT_sel);
        }
      }
    };
    //
    if (argc == 4) {
      do_print(std::cerr);
    } else {
      std::string FileO = argv[4];
      std::ofstream os(FileO);
      do_print(os);
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}

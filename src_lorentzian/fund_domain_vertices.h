#ifndef SRC_LORENTZIAN_FUND_DOMAIN_VERTICES_H_
#define SRC_LORENTZIAN_FUND_DOMAIN_VERTICES_H_

#include <limits>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

template <typename T, typename Tint> struct FundDomainVertex {
  MyVector<T> gen;
  MyMatrix<Tint> MatRoot;
};

template <typename T, typename Tint>
bool operator==(FundDomainVertex<T, Tint> const &k1,
                FundDomainVertex<T, Tint> const &k2) {
  MyVector<T> gen1 = RemoveFractionVector(k1.gen);
  MyVector<T> gen2 = RemoveFractionVector(k2.gen);
  if (gen1 != gen2)
    return false;
  if (k1.MatRoot.rows() != k2.MatRoot.rows())
    return false;
  std::set<MyVector<Tint>> set1;
  for (int i = 0; i < k1.MatRoot.rows(); i++) {
    MyVector<Tint> V = GetMatrixRow(k1.MatRoot, i);
    set1.insert(V);
  }
  std::set<MyVector<Tint>> set2;
  for (int i = 0; i < k2.MatRoot.rows(); i++) {
    MyVector<Tint> V = GetMatrixRow(k2.MatRoot, i);
    set2.insert(V);
  }
  return set1 == set2;
}

template <typename T, typename Tint>
bool operator!=(FundDomainVertex<T, Tint> const &k1,
                FundDomainVertex<T, Tint> const &k2) {
  return !(k1 == k2);
}

template <typename T, typename Tint>
std::string StringFundDomainVertexGAP(FundDomainVertex<T, Tint> const &vert) {
  std::string ret = "rec(gen:=" + StringVectorGAP(vert.gen) +
                    ", l_roots:=" + StringMatrixGAP(vert.MatRoot) + ")";
  return ret;
}

template <typename T, typename Tint>
void WriteFundDomainVertex(MyMatrix<T> const &G,
                           FundDomainVertex<T, Tint> const &vert,
                           std::ostream &os, std::string const &OutFormat) {
  const MyMatrix<Tint> &Mroot = vert.MatRoot;
  MyVector<T> gen_nf = RemoveFractionVector(vert.gen);
  T norm = gen_nf.dot(G * gen_nf);
  if (OutFormat == "GAP") {
    os << "rec(gen:=";
    WriteVectorGAP(os, gen_nf);
    os << ", norm:=" << norm << ", l_roots:=";
    WriteMatrixGAP(os, Mroot);
    os << ")";
    return;
  }
  if (OutFormat == "TXT") {
    os << "gen=";
    WriteVector(os, gen_nf);
    os << "l_roots=\n";
    WriteMatrixGAP(os, Mroot);
    return;
  }
  std::cerr << "Failed to find a matching entry for WriteFundDomainVertex\n";
  std::cerr << "OutFormat=" << OutFormat
            << " but allowed values are GAP and TXT\n";
  throw TerminalException{1};
}

template <typename T>
using pair_char =
    std::pair<MyMatrix<T>, WeightMatrix<true, std::vector<T>, uint16_t>>;

template <typename T, typename Tint, typename Tgroup>
struct FundDomainVertex_FullInfo {
  FundDomainVertex<T, Tint> vert;
  pair_char<T> e_pair_char;
  size_t hash;
  Tgroup GRP1;
  Tgroup GRP1_integral;
  std::string method;
};

template <typename T, typename Tint, typename Tgroup> struct ret_type {
  pair_char<T> e_pair_char;
  Tgroup GRP1;
  MyMatrix<Tint> MatRoot;
};

template <typename T, typename Tint> struct InitialComputation {
  T norm;
  std::unordered_map<MyVector<Tint>, uint8_t> map_v;
  std::vector<MyMatrix<T>> ListMat;
};

template <typename T, typename Tint>
InitialComputation<T, Tint>
GetInitialComputation(MyMatrix<T> const &G,
                      FundDomainVertex<T, Tint> const &vert) {
  std::unordered_map<MyVector<Tint>, uint8_t> map_v;
  size_t len = vert.MatRoot.rows();
  for (size_t i = 0; i < len; i++) {
    MyVector<Tint> eV = GetMatrixRow(vert.MatRoot, i);
    map_v[eV] = 1;
  }
  MyVector<Tint> gen_tint =
      UniversalVectorConversion<Tint, T>(RemoveFractionVector(vert.gen));
  map_v[gen_tint] = 2;
  T norm = vert.gen.dot(G * vert.gen);
  //
  std::vector<MyMatrix<T>> ListMat{G};
  if (norm == 0) {
    MyMatrix<T> MatRoot_T = UniversalMatrixConversion<T, Tint>(vert.MatRoot);
    MyMatrix<T> Qmat = GetQmatrix_NotFullRank(MatRoot_T);
    std::cerr << "Qmat=\n";
    WriteMatrix(std::cerr, Qmat);
    std::cerr << "MatRoot=\n";
    WriteMatrix(std::cerr, MatRoot_T);
    MyMatrix<T> Hmat1 = MatRoot_T * Qmat * MatRoot_T.transpose();
    std::cerr << "Hmat1=\n";
    WriteMatrix(std::cerr, Hmat1);
    MyMatrix<T> Hmat2 = MatRoot_T * G * MatRoot_T.transpose();
    std::cerr << "Hmat2=\n";
    WriteMatrix(std::cerr, Hmat2);
    ListMat.emplace_back(std::move(Qmat));
  }
  return {norm, map_v, ListMat};
}

template <typename T, typename Tint, typename Tgroup>
ret_type<T, Tint, Tgroup> get_canonicalized_record(
    std::vector<MyMatrix<T>> const &ListMat,
    std::unordered_map<MyVector<Tint>, uint8_t> const &the_map) {
  std::cerr << "Beginning of get_canonicalized_record\n";
  using Tidx_value = uint16_t;
  using Tidx = uint32_t;
  using Telt = typename Tgroup::Telt;
  using Telt_idx = typename Telt::Tidx;
  using Tgr = GraphListAdj;
  size_t n_row = the_map.size();
  std::vector<MyVector<Tint>> l_vect;
  std::vector<T> Vdiag;
  size_t idx = 0;
  std::vector<size_t> V;
  for (auto &kv : the_map) {
    std::cerr << "kv=" << kv.first << " / " << static_cast<int>(kv.second) << "\n";
    l_vect.push_back(kv.first);
    Vdiag.push_back(T(kv.second));
    if (kv.second == 1)
      V.push_back(idx);
    idx++;
  }
  size_t n1 = V.size();
  MyMatrix<T> MatV =
      UniversalMatrixConversion<T, Tint>(MatrixFromVectorFamily(l_vect));
  MyMatrix<T> MatV_red = ColumnReduction(MatV);
  WeightMatrix<true, std::vector<T>, Tidx_value> WMat =
      GetWeightMatrix_ListMat_Vdiag<T, Tidx, Tidx_value>(MatV, ListMat, Vdiag);
  WMat.ReorderingSetWeight();
  std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> epair =
      GetGroupCanonicalizationVector_Kernel<std::vector<T>, Tgr, Tidx,
                                            Tidx_value>(WMat);
  std::cerr << "We have epair\n";
  const std::vector<Tidx> &ListIdx = epair.first;
  const std::vector<std::vector<Tidx>> &ListGen = epair.second;
  WMat.RowColumnReordering(ListIdx);
  std::vector<Tidx> ListIdxRev(n_row);
  for (size_t i = 0; i < n_row; i++)
    ListIdxRev[ListIdx[i]] = i;
  std::cerr << "     ePermIdx=" << Telt(ListIdx) << "\n";
  std::cerr << "inv(ePermIdx)=" << Telt(ListIdxRev) << "\n";
  std::vector<MyVector<Tint>> l_vect_reord(n_row);
  std::vector<T> Vdiag_reord(n_row);
  std::vector<MyVector<Tint>> l_vect1;
  std::vector<size_t> Map1;
  std::vector<size_t> Map1_rev(n_row, std::numeric_limits<size_t>::max());
  size_t pos = 0;
  for (size_t i = 0; i < n_row; i++) {
    size_t j = ListIdx[i];
    MyVector<Tint> const& eV = l_vect[j];
    l_vect_reord[i] = eV;
    Vdiag_reord[i] = Vdiag[j];
    if (Vdiag[j] == 1) {
      l_vect1.push_back(eV);
      Map1.push_back(i);
      Map1_rev[i] = pos;
      pos++;
    }
    std::cerr << "i=" << i << " evect=" << StringVectorGAP(eV) << " diag=" << Vdiag[j] << "\n";
  }
  std::cerr << "Map1 =";
  for (auto & eval : Map1)
    std::cerr << " " << int(eval);
  std::cerr << "\n";
  std::cerr << "Map1_rev =";
  for (auto & eval : Map1_rev)
    std::cerr << " " << int(eval);
  std::cerr << "\n";
  MyMatrix<T> MatV_reord =
      UniversalMatrixConversion<T, Tint>(MatrixFromVectorFamily(l_vect_reord));
  MyMatrix<T> MatV_reord_red = ColumnReduction(MatV_reord);
  // There are two use case of computing the group
  // ---For the computation of minimal adjacencies that would get us a full
  // dimensional system
  // ---For the computation of orbits of adjacent vertices
  // So, in both cases, we need to reduce to the group for values 1.
  std::vector<Telt> LGen1;
  for (auto &eGen : ListGen) {
    MyMatrix<T> eEquivMat = RepresentVertexPermutation(MatV_red, MatV_red, eGen);
    std::optional<MyMatrix<T>> opt_EquivMat = FindMatrixTransformationTest(MatV_red, MatV_red, eGen);
    if (!opt_EquivMat) {
      std::cerr << "Failed to find opt_EquivMat\n";
      throw TerminalException{1};
    }
    std::cerr << "get_canonicalized_record eGen=" << Telt(eGen) << " eEquivMat=\n";
    WriteMatrix(std::cerr, eEquivMat);
    std::cerr << "*opt_EquivMat=\n";
    WriteMatrix(std::cerr, *opt_EquivMat);
    //
    std::vector<Telt_idx> Vloc(n_row);
    for (size_t i1 = 0; i1 < n_row; i1++) {
      Tidx i2 = ListIdx[i1];
      Tidx i3 = eGen[i2];
      Tidx i4 = ListIdxRev[i3];
      Vloc[i1] = i4;
    }
    Telt eGenReord(Vloc);
    //
    MyMatrix<T> eEquivMatReord = RepresentVertexPermutation(MatV_reord_red, MatV_reord_red, eGenReord);
    std::cerr << "get_canonicalized_record eGenReord=" << eGenReord << " eMat=\n";
    WriteMatrix(std::cerr, eEquivMatReord);
    if (!TestEqualityMatrix(eEquivMatReord, eEquivMat)) {
      std::cerr << "eEquivMat and eEquivMatReord should be equal\n";
      throw TerminalException{1};
    }
    //
    std::vector<Telt_idx> V1(n1);
    for (size_t i1 = 0; i1 < n1; i1++) {
      std::cerr << "i1=" << static_cast<int>(i1) << "\n";
      size_t i2 = Map1[i1];
      std::cerr << "  i2=" << static_cast<int>(i2) << "\n";
      size_t i3 = Vloc[i2];
      std::cerr << "  i3=" << static_cast<int>(i3) << "\n";
      size_t i4 = Map1_rev[i3];
      std::cerr << "  i4=" << static_cast<int>(i4) << "\n";
      V1[i1] = i4;
    }
    Telt ePerm1(V1);
    LGen1.push_back(ePerm1);
  }
  Tgroup GRP1(LGen1, n1);
  std::cerr << "GRP1 has been built\n";
  MyMatrix<Tint> MatV_ret = MatrixFromVectorFamily(l_vect1);
  std::cerr << "We have MatV_ret=\n";
  WriteMatrix(std::cerr, MatV_ret);
  MyMatrix<T> MatV_ret_red = UniversalMatrixConversion<T,Tint>(ColumnReduction(MatV_ret));
  std::cerr << "We have MatV_ret_red\n";
  for (auto & eGen : GRP1.GeneratorsOfGroup()) {
    std::cerr << "Treating generator eGen=" << eGen << "\n";
    MyMatrix<T> eMat = RepresentVertexPermutation(MatV_ret_red, MatV_ret_red, eGen);
    std::cerr << "get_canonicalized_record eMat=\n";
    WriteMatrix(std::cerr, eMat);
  }
  pair_char<T> e_pair_char{std::move(MatV_reord), std::move(WMat)};
  return {std::move(e_pair_char), std::move(GRP1), std::move(MatV_ret)};
}

template <typename T, typename Tint, typename Tgroup>
FundDomainVertex_FullInfo<T, Tint, Tgroup>
get_full_info(FundDomainVertex<T, Tint> const &vert,
              ret_type<T, Tint, Tgroup> &frec, std::string const &method) {
  auto &WMat = frec.e_pair_char.second;
  //  std::cerr << "frec.e_pair_char.first=\n";
  //  WriteMatrix(std::cerr, frec.e_pair_char.first);
  //  std::cerr << "RankMat = " << RankMat(frec.e_pair_char.first) << "\n";
  //  std::cerr << "frec.e_pair_char.second=\n";
  //  PrintWeightedMatrix(std::cerr, frec.e_pair_char.second);
  size_t seed = 1440;
  size_t hash = ComputeHashWeightMatrix_raw(WMat, seed);
  FundDomainVertex<T, Tint> new_vert{vert.gen, frec.MatRoot};
  //  std::cerr << "gen_fund_domain_fund_info gen=" << StringVectorGAP(vert.gen)
  //  << " |GRP1|=" << frec.GRP1.size() << "\n"; for (auto & eGen :
  //  frec.GRP1.GeneratorsOfGroup())
  //    std::cerr << "eGen=" << eGen << "\n";
  return {std::move(new_vert),
          std::move(frec.e_pair_char),
          hash,
          std::move(frec.GRP1),
          {},
          method};
}

/*
template <typename T, typename Tint, typename Tgroup>
FundDomainVertex_FullInfo<T, Tint, Tgroup>
get_fund_domain_full_info_noiso(MyMatrix<T> const &G,
                                FundDomainVertex<T, Tint> const &vert) {
  InitialComputation<T, Tint> ic = GetInitialComputation(G, vert);
  std::string method = "extendedvectfamily";
  ret_type<T, Tint, Tgroup> frec =
      get_canonicalized_record<T, Tint, Tgroup>(ic.ListMat, ic.map_v);
  return get_full_info(vert, frec, method);
}
*/

template <typename T, typename Tint, typename Tgroup>
FundDomainVertex_FullInfo<T, Tint, Tgroup>
DirectCopy(FundDomainVertex_FullInfo<T, Tint, Tgroup> const &fdfi) {
  return {fdfi.vert,
          {fdfi.e_pair_char.first, fdfi.e_pair_char.second.DirectCopy()},
          fdfi.hash,
          fdfi.GRP1,
          fdfi.GRP1_integral,
          fdfi.method};
}


template<typename T, typename Tint>
MyMatrix<T> ComputeSpanningSpace(MyMatrix<T> const& M)
{
  MyMatrix<T> Orth = NullspaceIntTrMat(M);
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
  std::cerr << "Orth=\n";
  WriteMatrix(std::cerr, Orth);
#endif
  MyMatrix<T> TheSpann_pre = NullspaceIntTrMat(Orth);
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
  std::cerr << "TheSpann_pre\n";
  WriteMatrix(std::cerr, TheSpann_pre);
#endif
  MyMatrix<T> TheSpann = LLLbasisReduction<T, Tint>(TheSpann_pre).LattRed;
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
  std::cerr << "TheSpann\n";
  WriteMatrix(std::cerr, TheSpann);
#endif
  CanSolIntMat<T> eCan = ComputeCanonicalFormFastReduction(TheSpann);
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
  std::cerr << "We have eCan\n";
#endif
  std::optional<MyMatrix<T>> opt = CanSolutionIntMatMat(eCan, M);
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
  std::cerr << "We have opt\n";
#endif
#ifdef CHECK_LORENTZIAN_STAB_EQUIV
  if (!opt) {
    std::cerr << "opt should have been assigned\n";
    throw TerminalException{1};
  }
#endif
  return *opt;
}


template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<T>> LORENTZ_GetStabilizerGenerator(
    MyMatrix<T> const &G,
    FundDomainVertex_FullInfo<T, Tint, Tgroup> const &vertFull) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  MyMatrix<Tint> const &MatRoot = vertFull.vert.MatRoot;
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
  std::cerr << "LORENTZ_GetStabilizerGenerator, vertFull.method="
            << vertFull.method << " gen=" << StringVectorGAP(vertFull.vert.gen)
            << "\n";
  std::cerr << "gen=" << StringVector(vertFull.vert.gen) << "\n";
  std::cerr << "LORENTZ_GetStabilizerGenerator MatRoot=\n";
  WriteMatrix(std::cerr, MatRoot);
#endif
  if (vertFull.method == "extendedvectfamily") {
    return LinPolytopeIntegralWMat_Automorphism<T, Tgroup, std::vector<T>,
                                                uint16_t>(vertFull.e_pair_char);
  }
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
  std::cerr << "|GRP1|=" << vertFull.GRP1.size() << "\n";
#endif
  if (vertFull.method == "isotropstabequiv_V1" ||
      vertFull.method == "isotropstabequiv") {
    int n = G.rows();
    std::vector<MyMatrix<T>> LGen1;
    MyMatrix<T> Subspace1 = UniversalMatrixConversion<T, Tint>(MatRoot);
    MyMatrix<T> Subspace1_proj = ComputeSpanningSpace<T, Tint>(Subspace1);
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
    std::cerr << "Subspace1_proj=\n";
    WriteMatrix(std::cerr, Subspace1_proj);
#endif
    std::vector<MyMatrix<T>> LGen1_B;
    for (auto & eGen : vertFull.GRP1.GeneratorsOfGroup()) {
      MyMatrix<T> eGen_M = RepresentVertexPermutation(Subspace1_proj, Subspace1_proj, eGen);
      std::cerr << "Gen_M=\n";
      WriteMatrix(std::cerr, eGen_M);
      LGen1_B.push_back(eGen_M);
    }
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
    std::cerr << "We have LGen1_B\n";
#endif
    std::vector<MyMatrix<T>> LGen1_C = LinPolytopeIntegral_Automorphism_Subspaces<T,Tgroup>(LGen1_B, Subspace1_proj);
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
    std::cerr << "We have LGen1_C\n";
#endif
    std::vector<Telt> LGen1_D;
    for (auto & eGen : LGen1_C) {
      MyMatrix<T> Subspace1_projImg = Subspace1_proj * eGen;
      Telt ePerm = GetPermutationOnVectors<T,Telt>(Subspace1_proj, Subspace1_projImg);
      LGen1_D.push_back(ePerm);
    }
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
    std::cerr << "We have LGen1_D\n";
#endif
    //    MyMatrix<T> Subspace1red = ColumnReduction(Subspace1);
    //    Tgroup GRPlin =
    //    LinPolytope_Automorphism<T,false,Tgroup>(Subspace1red);
    int nRow = Subspace1.rows();
    Tidx nRow_tidx = nRow;
#ifdef TRACK_INFOS_LOG
    std::string strListGen = "[";
    bool IsFirst = true;
#endif
    for (auto &eGen : LGen1_D) {
#ifdef TRACK_INFOS_LOG
      if (!IsFirst)
        strListGen += ",";
      IsFirst = false;
      strListGen += std::to_string(eGen);
#endif
      //      std::cerr << "eGen=" << eGen << "\n";
      //      bool test = GRPlin.isin(eGen);
      //      std::cerr << "test=" << test << "\n";
      //      MyMatrix<T> eGen_as_matr = FindTransformation(Subspace1red,
      //      Subspace1red, eGen); std::cerr << "We have eGen_as_matr\n";
      MyMatrix<T> Subspace2(nRow, Subspace1.cols());
      for (Tidx iRow = 0; iRow < nRow_tidx; iRow++) {
        Tidx jRow = eGen.at(iRow);
        MyVector<T> V = GetMatrixRow(Subspace1, iRow);
        AssignMatrixRow(Subspace2, jRow, V);
      }
      std::optional<MyMatrix<T>> opt =
          ExtendOrthogonalIsotropicIsomorphism(G, Subspace1, G, Subspace2);
#ifdef CHECK_LORENTZIAN_STAB_EQUIV
      if (!opt) {
        std::cerr << "opt found to be missing\n";
        throw TerminalException{1};
      }
#endif
      MyMatrix<T> const &eGen1 = *opt;
      LGen1.push_back(eGen1);
    }
#ifdef TRACK_INFOS_LOG
    strListGen += "]";
    std::string strGroup = "Group(" + strListGen + ")";
    std::cout << "rec(group:=" << strGroup << "),\n";
#endif
    MyMatrix<T> InvariantSpace = MatrixIntegral_GetInvariantSpace(n, LGen1);
    MyMatrix<T> InvInvariantSpace = Inverse(InvariantSpace);
    std::vector<MyMatrix<T>> LGen2;
    for (auto &eGen1 : LGen1) {
      MyMatrix<T> eGen2 = InvariantSpace * eGen1 * InvInvariantSpace;
#ifdef CHECK_LORENTZIAN_STAB_EQUIV
      if (!IsIntegralMatrix(eGen2)) {
        std::cerr << "The matrix eGen2 should be integral\n";
        throw TerminalException{1};
      }
#endif
      LGen2.emplace_back(std::move(eGen2));
    }
    auto get_gen3 = [&]() -> std::vector<MyMatrix<T>> {
      if (vertFull.method == "isotropstabequiv_V1") {
        GeneralMatrixGroupHelper<T, Telt> helper{n};
        return LinearSpace_Stabilizer<T, Tgroup,
                                      GeneralMatrixGroupHelper<T, Telt>>(
            LGen2, helper, InvInvariantSpace);
      } else {
        MyVector<T> const &Visotrop = vertFull.vert.gen;
        MyMatrix<T> eProd = Subspace1 * InvInvariantSpace;
        MyMatrix<T> G_new = InvariantSpace * G * InvariantSpace.transpose();
        FiniteIsotropicMatrixGroupHelper<T, Telt> helper =
            ComputeFiniteIsotropicMatrixGroupHelper<T, Telt>(G_new, eProd,
                                                             Visotrop);
        return LinearSpace_Stabilizer<
            T, Tgroup, FiniteIsotropicMatrixGroupHelper<T, Telt>>(
            LGen2, helper, InvInvariantSpace);
      }
    };
#ifdef CHECK_LORENTZIAN_STAB_EQUIV
    std::vector<MyVector<Tint>> ListV;
    std::unordered_set<MyVector<Tint>> SetV;
    for (int i = 0; i < MatRoot.rows(); i++) {
      MyVector<Tint> V = GetMatrixRow(MatRoot, i);
      ListV.push_back(V);
      SetV.insert(V);
    }
#endif
    std::vector<MyMatrix<T>> LGen4;
    for (auto &eGen3 : get_gen3()) {
      MyMatrix<T> eGen4 = InvInvariantSpace * eGen3 * InvariantSpace;
#ifdef CHECK_LORENTZIAN_STAB_EQUIV
      if (!IsIntegralMatrix(eGen4)) {
        std::cerr << "The matrix eGen4 should be integral\n";
        throw TerminalException{1};
      }
      MyMatrix<Tint> eGen4_i = UniversalMatrixConversion<Tint, T>(eGen4);
      for (int i = 0; i < MatRoot.rows(); i++) {
        MyVector<Tint> Vimg = eGen4_i.transpose() * ListV[i];
        if (SetV.count(Vimg) == 0) {
          std::cerr << "The vertor at i=" << i << " is not mapped in MatRoot\n";
          throw TerminalException{1};
        }
      }
#endif
      LGen4.emplace_back(std::move(eGen4));
    }
    return LGen4;
  }
  std::cerr << "Error in LORENTZ_GetStabilizerGenerator\n";
  throw TerminalException{1};
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<T>> LORENTZ_TestEquivalence(
    MyMatrix<T> const &G1,
    FundDomainVertex_FullInfo<T, Tint, Tgroup> const &vertFull1,
    MyMatrix<T> const &G2,
    FundDomainVertex_FullInfo<T, Tint, Tgroup> const &vertFull2) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
#ifdef PRINT_DEBUG_STAB_EQUIV_INFO
  std::cerr << "LORENTZ_TestEquivalence, vertFull1.method=" << vertFull1.method
            << "\n";
  std::cerr << "LORENTZ_TestEquivalence, gen1="
            << StringVector(vertFull1.vert.gen)
            << " gen2=" << StringVector(vertFull2.vert.gen) << "\n";
#endif
  if (vertFull1.method != vertFull2.method) {
    return {};
  }
  if (vertFull1.method == "extendedvectfamily") {
    return LinPolytopeIntegralWMat_Isomorphism<T, Tgroup, std::vector<T>,
                                               uint16_t>(vertFull1.e_pair_char,
                                                         vertFull2.e_pair_char);
  }
  if (vertFull1.method == "isotropstabequiv_V1" ||
      vertFull1.method == "isotropstabequiv") {
    if (vertFull1.e_pair_char.first.rows() !=
        vertFull2.e_pair_char.first.rows())
      return {};
    if (vertFull1.e_pair_char.second.GetWeight() !=
        vertFull2.e_pair_char.second.GetWeight())
      return {};
    MyMatrix<T> Subspace1 =
        UniversalMatrixConversion<T, Tint>(vertFull1.vert.MatRoot);
    MyMatrix<T> Subspace2 =
        UniversalMatrixConversion<T, Tint>(vertFull2.vert.MatRoot);
    MyMatrix<T> Subspace1_proj = ComputeSpanningSpace<T, Tint>(Subspace1);
    MyMatrix<T> Subspace2_proj = ComputeSpanningSpace<T, Tint>(Subspace1);
    //    std::cerr << "Subspace1=\n";
    //    WriteMatrix(std::cerr, Subspace1);
    //    std::cerr << "Subspace2=\n";
    //    WriteMatrix(std::cerr, Subspace2);
    std::optional<MyMatrix<T>> opt1 =
        ExtendOrthogonalIsotropicIsomorphism(G1, Subspace1, G1, Subspace2);
    if (!opt1) {
      std::cerr << "Failed at extending equivalence\n";
      return {};
    }
    MyMatrix<T> const &EquivRat = *opt1;
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
    std::cerr << "EquivRat=" << StringMatrixGAP_line(EquivRat) << "\n";
#endif
    //    WriteMatrix(std::cerr, EquivRat);
    //
    int n = G1.rows();
    std::vector<MyMatrix<T>> LGen1;
    int nRow = Subspace1.rows();
    Tidx nRow_tidx = nRow;
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
    std::cerr << "|GRP1|=" << vertFull1.GRP1.size() << "\n";
#endif
    std::vector<std::string> LGenStr;
    for (auto &eGen : vertFull1.GRP1.GeneratorsOfGroup()) {
      MyMatrix<T> Subspace2(nRow, Subspace1.cols());
      for (Tidx iRow = 0; iRow < nRow_tidx; iRow++) {
        Tidx jRow = eGen.at(iRow);
        MyVector<T> V = GetMatrixRow(Subspace1, iRow);
        AssignMatrixRow(Subspace2, jRow, V);
      }
      std::optional<MyMatrix<T>> opt =
          ExtendOrthogonalIsotropicIsomorphism(G1, Subspace1, G2, Subspace2);
#ifdef CHECK_LORENTZIAN_STAB_EQUIV
      if (!opt) {
        std::cerr << "opt should point to something\n";
        throw TerminalException{1};
      }
#endif
      MyMatrix<T> const &eGen1 = *opt;
      //      std::cerr << "eGen1=\n";
      //      WriteMatrix(std::cerr, eGen1);
      LGenStr.push_back(StringMatrixGAP_line(eGen1));
      LGen1.emplace_back(std::move(eGen1));
    }
    // Original question: Does there exist g in GRP(LGen1) s.t. g * EquivRat in
    // GLn(Z)
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
    std::cerr << "LGen1=" << StringVectorStringGAP(LGenStr) << "\n";
    std::cerr << "We have LGen1\n";
#endif
    MyMatrix<T> InvariantSpace = MatrixIntegral_GetInvariantSpace(n, LGen1);
    MyMatrix<T> InvariantSpaceInv = Inverse(InvariantSpace);
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
    std::cerr << "InvariantSpace=\n";
    WriteMatrix(std::cerr, InvariantSpace);
    std::cerr << "InvariantSpaceInv=\n";
    WriteMatrix(std::cerr, InvariantSpaceInv);
#endif
    std::vector<MyMatrix<T>> LGen2;
    for (auto &eGen1 : LGen1) {
      MyMatrix<T> eGen2 = InvariantSpace * eGen1 * InvariantSpaceInv;
#ifdef CHECK_LORENTZIAN_STAB_EQUIV
      if (!IsIntegralMatrix(eGen2)) {
        std::cerr << "The matrix eGen2 should be integral\n";
        throw TerminalException{1};
      }
#endif
      LGen2.emplace_back(std::move(eGen2));
    }
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
    std::cerr << "We have LGen2\n";
#endif
    //
    MyMatrix<T> InvariantSpaceImg = InvariantSpace * EquivRat;
    MyMatrix<T> InvariantSpaceImgInv = Inverse(InvariantSpaceImg);
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
    std::cerr << "We have InvariantSpaceImg\n";
#endif

    auto get_opt2 = [&]() -> std::optional<MyMatrix<T>> {
      if (vertFull1.method == "isotropstabequiv_V1") {
        GeneralMatrixGroupHelper<T, Telt> helper{n};
        return LinearSpace_Equivalence<T, Tgroup,
                                       GeneralMatrixGroupHelper<T, Telt>>(
            LGen2, helper, InvariantSpaceInv, InvariantSpaceImgInv);
      } else {
        MyMatrix<T> eProd = Subspace1 * InvariantSpaceInv;
        MyMatrix<T> G1_new = InvariantSpace * G1 * InvariantSpace.transpose();
        MyVector<T> const &Visotrop = vertFull1.vert.gen;
        FiniteIsotropicMatrixGroupHelper<T, Telt> helper =
            ComputeFiniteIsotropicMatrixGroupHelper<T, Telt>(G1_new, eProd,
                                                             Visotrop);
        return LinearSpace_Equivalence<
            T, Tgroup, FiniteIsotropicMatrixGroupHelper<T, Telt>>(
            LGen2, helper, InvariantSpaceInv, InvariantSpaceImgInv);
      }
    };
    std::optional<MyMatrix<T>> opt2 = get_opt2();
#ifdef DEBUG_LORENTZIAN_STAB_EQUIV
    std::cerr << "We have opt2\n";
#endif
    if (!opt2)
      return {};
    //
    MyMatrix<T> const &eSpaceEquiv = *opt2;
    MyMatrix<T> eMatFinal = InvariantSpaceInv * eSpaceEquiv * InvariantSpace;
    MyMatrix<T> eProd = eMatFinal * EquivRat;
#ifdef CHECK_LORENTZIAN_STAB_EQUIV
    if (!IsIntegralMatrix(eProd)) {
      std::cerr << "eProd should be integral\n";
      throw TerminalException{1};
    }
#endif
    return eProd;
  }
  std::cerr << "Error in LORENTZ_TestEquivalence\n";
  throw TerminalException{1};
}

#endif // SRC_LORENTZIAN_FUND_DOMAIN_VERTICES_H_

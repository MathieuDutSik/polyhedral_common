#ifndef INCLUDE_EDGEWALK_H
#define INCLUDE_EDGEWALK_H

#include "MatrixCanonicalForm.h"
#include "Temp_PolytopeEquiStab.h"
#include "two_dim_lorentzian.h"
#include "coxeter_dynkin.h"
#include "vinberg_code.h"
#include "Namelist.h"
#include "Temp_Positivity.h"



FullNamelist NAMELIST_GetStandard_EDGEWALK()
{
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListIntValues1["n"]=9;
  ListStringValues1["FileLorMat"]="/tmp/unset";
  ListStringValues1["OptionInitialVertex"]="vinberg";
  ListStringValues1["FileVertDomain"]="/tmp/unset";
  ListStringValues1["OptionNorms"]="all";
  ListStringValues1["OutFormat"]="unset";
  ListStringValues1["FileOut"]="unset";
  SingleBlock BlockPROC;
  BlockPROC.ListStringValues=ListStringValues1;
  ListBlock["PROC"]=BlockPROC;
  // Merging all data
  return {ListBlock, "undefined"};
}


template<typename T>
MyMatrix<T> ComputeLattice_LN(MyMatrix<T> const& G, T const& N)
{
  int n = G.rows();
  MyMatrix<T> M1 = IdentityMat<T>(n);
  MyMatrix<T> M2 = (N / 2) * Inverse(G);
  return IntersectionLattice(M1, M2);
}


template<typename T, typename Tint>
struct FundDomainVertex {
  MyVector<T> gen;
  std::vector<MyVector<Tint>> l_roots;
};

template<typename T, typename Tint>
void WriteFundDomainVertex(FundDomainVertex<T,Tint> const& vert, std::ostream & os, std::string const& OutFormat)
{
  MyMatrix<Tint> Mroot = MatrixFromVectorFamily(vert.l_roots);
  if (OutFormat == "GAP") {
    os << "rec(gen:=";
    WriteMatrixGAP(os, vert.gen);
    os << ", l_roots:=";
    WriteMatrixGAP(os, Mroot);
    os << ")";
  }
  if (OutFormat == "TXT") {
    os << "gen=";
    WriteMatrix(os, vert.gen);
    os << "l_roots=\n";
    WriteMatrixGAP(os, Mroot);
  }
  std::cerr << "Failed to find a matching entry for WritePairVertices\n";
  throw TerminalException{1};
}






template<typename T, typename Tint>
struct RootCandidate {
  int sign; // 0 for 0, 1 for positive, -1 for negative
  T quant1; // this is (k.alpha_{N,\Delta'})^2 / R_{N,\Delta'}
  T quant2; // this is (k.alpha_{N,\Delta'})^2 / N
  T e_norm;
  MyVector<Tint> alpha;
};

template<typename T>
int get_sign_sing(T const& val)
{
  if (val > 0)
    return 1;
  if (val < 0)
    return -1;
  return 0;
}


template<typename T>
int get_sign_pair_t(T const& p1, T const& p2)
{
  if (p1 < p2)
    return 1;
  if (p1 > p2)
    return -1;
  return 0;
}



// return 0 is p1 == p2 :  1 if p1 < p2 : -1 if p1 > p2
template<typename T>
int get_sign_pair_stdpair(std::pair<int,T> const& p1, std::pair<int,T> const& p2)
{
  if (p1.first == p2.first)
    return get_sign_pair_t(p1.second, p2.second);
  return get_sign_pair_t(p1.first, p2.first);
}


template<typename T, typename Tint>
RootCandidate<T,Tint> gen_possible_extension(MyMatrix<T> const& G, MyVector<T> const& k, MyVector<Tint> const& alpha, T const& res_norm, T const& e_norm)
{
  MyVector<T> alpha_T = UniversalVectorConversion<T,Tint>(alpha);
  T scal = - k.dot(G * alpha_T);
  T quant1 = (scal * scal) / res_norm;
  T quant2 = (scal * scal) / e_norm;
  return {get_sign_sing(scal), quant1, quant2, e_norm, alpha};
}


// sign as before + means preferable according to page 27.
// return 1 if poss1 is preferable to poss2
template<typename T, typename Tint>
int get_sign_cand(RootCandidate<T,Tint> const& poss1, RootCandidate<T,Tint> const& poss2)
{
  int sign1 = get_sign_pair_stdpair<T>({poss1.sign, poss1.quant1}, {poss2.sign, poss2.quant1});
  if (sign1 != 0)
    return sign1; // because -k.alpha1 / sqrt(R1)    <     -k.alpha2 / sqrt(R2)   correspond to 1 in the above.
  int sign2 = get_sign_pair_stdpair<T>({poss1.sign, poss1.quant2}, {poss2.sign, poss2.quant2});
  if (sign2 != 0)
    return sign2; // because -k.alpha1 / sqrt(N1)    <     -k.alpha2 / sqrt(N2)   correspond to 1 in the above.
  int sign3 = get_sign_pair_t(poss1.e_norm, poss2.e_norm);
  return sign3; // because N1 < N2 corresponds to 1
}



template<typename T, typename Tint>
RootCandidate<T,Tint> get_best_candidate(std::vector<RootCandidate<T,Tint>> const& l_cand)
{
  if (l_cand.size() == 0) {
    std::cerr << "We have zero candidates. Abort\n";
    throw TerminalException{1};
  }
  RootCandidate<T,Tint> best_cand = l_cand[0];
  for (auto & eCand : l_cand) {
    if (get_sign_cand(eCand, best_cand) == 1) {
      best_cand = eCand;
    }
  }
  return best_cand;
}






/*
  We take the notations as in EDGEWALK paper.
  ---The dimension is n+1
  ---We have an edge e between two rays k and k'.
  We have u1, .... u(n-1) roots that are orthogonal and U the real span of those vectors.
  ---P is the real plane orthogonal to U.
  ---pi_U and pi_P are the corresponding projectors
  ---How is (1/2) P defined and correspond to k (typo correction)
  ---


 */
template<typename T, typename Tint>
FundDomainVertex<T,Tint> EdgewalkProcedure(MyMatrix<T> const& G, MyVector<T> const& k, std::vector<MyVector<Tint>> l_ui, std::vector<T> const& l_norms, MyVector<Tint> const& v_disc)
{
  int n = G.size();
  size_t n_root = l_ui.size();
  MyMatrix<T> Space(n_root,n);
  MyMatrix<T> EquaB(n_root+1,n);
  for (size_t i_root=0; i_root<n_root; i_root++) {
    MyVector<T> eV = UniversalVectorConversion<T,Tint>(l_ui[i_root]);
    MyVector<T> eP = G * eV;
    T eScal = k.dot(eP);
    if (eScal != 0) {
      std::cerr << "The scalar product should be 0\n";
      throw TerminalException{1};
    }
    for (int i=0; i<n; i++) {
      Space(i_root,i) = eV(i);
      EquaB(i_root,i) = eP(i);
    }
  }
  MyVector<T> eP = G * k;
  for (int i=0; i<n; i++)
    EquaB(n_root,i) = eP(i);
  MyMatrix<T> NSP = NullspaceMat(EquaB);
  if (NSP.rows() != 1) {
    std::cerr << "The dimension should be exactly 2\n";
    throw TerminalException{1};
  }
  MyVector<T> r0 = GetMatrixRow(NSP,0);
  std::vector<RootCandidate<T,Tint>> l_candidates;
  bool allow_euclidean = false;
  std::vector<Possible_Extension<T>> l_extension = ComputePossibleExtensions(G, l_ui, l_norms, allow_euclidean);
  for (auto & e_extension : l_extension) {
    T e_norm = e_extension.e_norm;
    MyMatrix<T> Latt = ComputeLattice_LN(G, e_norm);
    // Now getting into the LN space
    MyMatrix<T> Space_LN = Space * Inverse(Latt);
    MyMatrix<T> G_LN = Latt * G * Latt.transpose();
    MyMatrix<T> Equas = Space_LN * G_LN;
    MyMatrix<T> NSP = NullspaceIntMat(TransposedMat(Equas));
    MyMatrix<T> GP_LN = NSP * G_LN * NSP.transpose();
    MyVector<T> r0_LN = Inverse(Latt).transpose() * r0;
    auto RecSol = SolutionMat(NSP, r0_LN);
    if (!RecSol.result) {
      std::cerr << "Failed to resolve the SolutionMat problem\n";
      throw TerminalException{1};
    }
    MyVector<T> r0_NSP = RecSol.eSol;
    MyVector<Tint> r0_work = UniversalVectorConversion<Tint,T>(RemoveFractionVector(r0_NSP));
    std::optional<MyVector<Tint>> opt_v = get_first_next_vector(GP_LN, r0_work, e_extension.res_norm);
    if (opt_v) {
      MyVector<T> v = UniversalVectorConversion<T,Tint>(*opt_v);
      MyVector<T> alpha_T = e_extension.u_component + NSP.transpose() * v;
      MyVector<Tint> alpha = UniversalVectorConversion<Tint,T>(alpha_T);
      RootCandidate<T,Tint> eCand = gen_possible_extension(G, k, alpha, e_extension.res_norm, e_norm);
      l_candidates.push_back(eCand);
    }
  }
  if (l_candidates.size() > 0) {
    RootCandidate<T,Tint> best_cand = get_best_candidate(l_candidates);
    std::vector<MyVector<Tint>> l_roots = l_ui;
    l_roots.push_back(best_cand.alpha);
    MyMatrix<T> Mat_root = UniversalMatrixConversion<T,Tint>(MatrixFromVectorFamily(l_roots));
    MyMatrix<T> EquaMat = Mat_root * G;
    MyMatrix<T> NSP = NullspaceTrMat(EquaMat);
    if (NSP.rows() != 1) {
      std::cerr << "We should have exactly one row\n";
      throw TerminalException{1};
    }
    MyVector<T> gen = GetMatrixRow(NSP, 0);
    T scal = gen.dot(G * k);
    if (scal > 0) {
      return {gen, l_roots};
    }
    if (scal < 0) {
      return {-gen, l_roots};
    }
    std::cerr << "Failed to find a matching entry\n";
    throw TerminalException{1};
  }
  // So, no candidates were found. We need to find isotropic vectors.
  MyMatrix<T> NSPbas(n,2);
  AssignMatrixRow(NSPbas, 0, k);
  AssignMatrixRow(NSPbas, 1, r0);
  MyMatrix<T> Gred = NSPbas * G * NSPbas.transpose();
  std::optional<MyMatrix<T>> Factor_opt = GetIsotropicFactorization(Gred);
  if (!Factor_opt) {
    std::cerr << "The matrix is not isotropic. Major rethink are needed\n";
    throw TerminalException{1};
  }
  MyMatrix<T> Factor = *Factor_opt;
  // We want a vector inside of the cone (there are two: C and -C)
  auto get_can_gen=[&](MyVector<T> const& v) -> MyVector<T> {
    T scal = k.dot(G * v);
    if (scal > 0)
      return v;
    if (scal < 0)
      return -v;
    std::cerr << "We should have scal != 0 to be able to conclude\n";
    throw TerminalException{1};
  };
  std::vector<MyVector<T>> l_gens;
  for (size_t i=0; i<2; i++) {
    // a x + by correspond to the ray (u0, u1) = (-b, a)
    T u0 = -Factor(i,1);
    T u1 =  Factor(i,0);
    MyVector<T> gen = u0 * k + u1 * r0;
    MyVector<T> can_gen = get_can_gen(gen);
    MyVector<T> v_disc_t = UniversalVectorConversion<T,Tint>(v_disc);
    T scal = v_disc_t.dot(G * can_gen);
    if (scal > 0)
      l_gens.push_back(can_gen);
  }
  if (l_gens.size() != 1) {
    std::cerr << "We should have just one vector in order to conclude. Rethink needed\n";
    throw TerminalException{1};
  }
  std::vector<MyVector<Tint>> l_roots = l_ui;
  MyVector<Tint> w(n);
  /// FILL OUT THE CODE
  l_roots.push_back(w);
  return {l_gens[0], l_roots};
}







template<typename T, typename Tint>
struct PairVertices {
  FundDomainVertex<T,Tint> vert1;
  FundDomainVertex<T,Tint> vert2;
  std::pair<MyMatrix<T>,WeightMatrix<true,T,uint16_t>> pair_char;
};

template<typename T, typename Tint>
void WritePairVertices(PairVertices<T,Tint> const& epair, std::ostream & os, std::string const& OutFormat)
{
  if (OutFormat == "GAP") {
    os << "rec(vert1:=";
    WriteFundDomainVertex(epair.vert1, os, OutFormat);
    os << ", vert2:=";
    WriteFundDomainVertex(epair.vert2, os, OutFormat);
    os << ")";
  }
  if (OutFormat == "TXT") {
    os << "vert&=\n";
    WriteFundDomainVertex(epair.vert1, os, OutFormat);
    os << "vert2=\n";
    WriteFundDomainVertex(epair.vert2, os, OutFormat);
  }
  std::cerr << "Failed to find a matching entry for WritePairVertices\n";
  throw TerminalException{1};
}


template<typename T, typename Tint>
PairVertices<T,Tint> gen_pair_vertices(MyMatrix<T> const& G, FundDomainVertex<T,Tint> const& vert1, FundDomainVertex<T,Tint> const& vert2)
{
  std::unordered_set<MyVector<Tint>> set_v;
  for (auto & eV : vert1.l_roots)
    set_v.insert(eV);
  for (auto & eV : vert2.l_roots)
    set_v.insert(eV);
  std::vector<MyVector<Tint>> l_roots;
  for (auto & eV : set_v)
    l_roots.push_back(eV);
  MyMatrix<T> MatV = UniversalMatrixConversion<T,Tint>(MatrixFromVectorFamily(l_roots));
  using Tidx_value = uint16_t;
  WeightMatrix<true, T, Tidx_value> WMat = GetSimpleWeightMatrix<T,Tidx_value>(MatV, G);
  std::pair<MyMatrix<T>,WeightMatrix<true,T,uint16_t>> pair_char{std::move(MatV),std::move(WMat)};
  return {vert1, vert2, std::move(pair_char)};
}



template<typename T, typename Tint>
struct ResultEdgewalk {
  std::vector<MyMatrix<Tint>> l_gen_isom_cox;
  std::vector<PairVertices<T,Tint>> l_orbit_pair_vertices;
};

template<typename T, typename Tint>
std::vector<T> get_list_norms(MyMatrix<T> const& G, ResultEdgewalk<T,Tint> const& re)
{
  std::set<T> set_norms;
  auto proc_vertex=[&](FundDomainVertex<T,Tint> const& vert) -> void {
    for (auto root : vert.l_roots) {
      MyVector<T> root_t = UniversalVectorConversion<T,Tint>(root);
      T norm = root_t.dot(G * root_t);
      set_norms.insert(norm);
    }
  };
  for (auto & e_pair : re.l_orbit_pair_vertices) {
    proc_vertex(e_pair.vert1);
    proc_vertex(e_pair.vert2);
  }
  std::vector<T> l_norms;
  for (auto & v : set_norms)
    l_norms.push_back(v);
  return l_norms;
}


template<typename T, typename Tint>
void PrintResultEdgewalk(MyMatrix<T> const& G, ResultEdgewalk<T,Tint> const& re, std::ostream& os, const std::string& OutFormat)
{
  std::vector<T> l_norms = get_list_norms(G, re);
  if (OutFormat == "GAP") {
    os << "return rec(l_norms:=";
    WriteStdVectorGAP(os, l_norms);
    os << ", ListIsomCox:=";
    WriteVectorMatrixGAP(os, re.l_gen_isom_cox);
    os << ", ListVertices:=";
    
  }
  if (OutFormat == "TXT") {
    os << "List of found generators of Isom / Cox\n";
  }
  std::cerr << "Failed to find a matching entry in PrintResultEdgewalk\n";
  throw TerminalException{1};
}





template<typename T, typename Tint, typename Tgroup>
ResultEdgewalk<T,Tint> LORENTZ_RunEdgewalkAlgorithm(MyMatrix<T> const& G, std::vector<T> const& l_norms, FundDomainVertex<T,Tint> const& eVert)
{
  std::vector<MyMatrix<Tint>> l_gen_isom_cox;
  struct EnumEntry {
    bool stat1;
    bool stat2;
    PairVertices<T,Tint> val;
  };
  std::vector<EnumEntry> l_entry;
  auto f_insert_gen=[&](MyMatrix<Tint> const& eP) -> void {
    MyMatrix<T> eP_T = UniversalMatrixConversion<T,Tint>(eP);
    MyMatrix<T> G2 = eP_T * G * eP_T.transpose();
    if (G2 != G) {
      std::cerr << "The matrix eP should leave the quadratic form invariant\n";
      throw TerminalException{1};
    }
    l_gen_isom_cox.push_back(eP);
  };
  auto func_insert_pair_vertices=[&](EnumEntry & v_pair) -> void {
    for (auto & u_pair : l_entry) {
      std::optional<MyMatrix<T>> equiv_opt = LinPolytopeWMat_Isomorphism<T,Tgroup,T,uint16_t>(u_pair.val.pair_char, v_pair.val.pair_char);
      if (equiv_opt) {
        f_insert_gen(UniversalMatrixConversion<Tint,T>(*equiv_opt));
        return;
      }
    }
    for (auto & eGen : LinPolytopeWMat_Automorphism<T,Tgroup,T,uint16_t>(v_pair.val.pair_char))
      f_insert_gen(UniversalMatrixConversion<Tint,T>(eGen));
    l_entry.emplace_back(std::move(v_pair));
  };
  size_t len = eVert.l_roots.size();
  auto insert_edges_from_vertex=[&](FundDomainVertex<T,Tint> const& theVert) -> void {
    for (size_t i=0; i<len; i++) {
      std::vector<MyVector<Tint>> l_ui;
      for (size_t j=0; j<len; j++) {
        if (i != j) {
          l_ui.push_back(theVert.l_roots[j]);
        }
      }
      MyVector<Tint> v_disc = theVert.l_roots[i];
      FundDomainVertex<T,Tint> fVert = EdgewalkProcedure(G, theVert.gen, l_ui, l_norms, v_disc);
      PairVertices<T,Tint> epair = gen_pair_vertices(G, theVert, fVert);
      EnumEntry entry{true, false, std::move(epair)};
      func_insert_pair_vertices(entry);
    }
  };
  while(true) {
    bool IsFinished = true;
    for (auto & entry : l_entry) {
      if (entry.stat1) {
        entry.stat1 = false;
        insert_edges_from_vertex(entry.val.vert1);
        IsFinished = false;
      }
      if (entry.stat2) {
        entry.stat2 = false;
        insert_edges_from_vertex(entry.val.vert2);
        IsFinished = false;
      }
    }
    if (IsFinished)
      break;
  }
  std::vector<PairVertices<T,Tint>> l_orbit_pair_vertices;
  for (auto & epair : l_entry)
    l_orbit_pair_vertices.emplace_back(std::move(epair.val));
  return {l_gen_isom_cox, std::move(l_orbit_pair_vertices)};
}




template<typename T, typename Tint>
std::vector<T> get_initial_list_norms(MyMatrix<T> const& G, std::string const& OptionNorms)
{
  if (OptionNorms == "K3")
    return {T(2)};
  if (OptionNorms == "all") {
    MyMatrix<Tint> G_Tint = UniversalMatrixConversion<Tint,T>(G);
    std::vector<Tint> l_norms_tint = Get_root_lengths(G_Tint);
    std::vector<T> l_norms;
    for (auto & eN : l_norms_tint)
      l_norms.push_back(T(eN));
    return l_norms;
  }
  std::cerr << "Failed to find a matching entry in get_initial_list_norms\n";
  throw TerminalException{1};
}



template<typename T, typename Tint>
FundDomainVertex<T,Tint> get_initial_vertex(MyMatrix<T> const& G, std::string const& OptionInitialVertex, std::string const& FileInitialVertex)
{
  if (OptionInitialVertex == "File") {
    std::ifstream is(FileInitialVertex);
    MyVector<T> gen = ReadVector<T>(is);
    MyMatrix<Tint> Mroot = ReadMatrix<Tint>(is);
    std::vector<MyVector<Tint>> l_roots;
    size_t n_root=Mroot.rows();
    for (size_t i=0; i<n_root; i++) {
      MyVector<Tint> root = GetMatrixRow(Mroot,i);
      l_roots.push_back(root);
    }
    return {gen, l_roots};
  }
  if (OptionInitialVertex == "vinberg") {
    VinbergTot<T,Tint> Vtot = GetVinbergFromG<T,Tint>(G);
    std::pair<MyVector<Tint>, std::vector<MyVector<Tint>>> epair = FindOneInitialRay(Vtot);
    return {UniversalVectorConversion<T,Tint>(epair.first), epair.second};
  }
  std::cerr << "Failed to find a matching entry in get_initial_list_norms\n";
  throw TerminalException{1};
}



template<typename T, typename Tint, typename Tgroup>
void MainFunctionEdgewalk(FullNamelist const& eFull)
{
  SingleBlock BlockPROC=eFull.ListBlock.at("PROC");
  std::string FileLorMat=BlockPROC.ListStringValues.at("FileLorMat");
  MyMatrix<T> G = ReadMatrixFile<T>(FileLorMat);
  DiagSymMat<T> DiagInfo = DiagonalizeNonDegenerateSymmetricMatrix(G);
  if (DiagInfo.nbZero != 0 || DiagInfo.nbMinus != 1) {
    std::cerr << "We have nbZero=" << DiagInfo.nbZero << " nbPlus=" << DiagInfo.nbPlus << " nbMinus=" << DiagInfo.nbMinus << "\n";
    std::cerr << "In the hyperbolic geometry we should have nbZero=0 and nbMinus=1\n";
    throw TerminalException{1};
  }
  //
  std::string OptionNorms=BlockPROC.ListStringValues.at("OptionIniti");
  std::vector<T> l_norms = get_initial_list_norms<T,Tint>(G, OptionNorms);
  //
  std::string OptionInitialVertex=BlockPROC.ListStringValues.at("OptionInitialVertex");
  std::string FileInitialVertex=BlockPROC.ListStringValues.at("FileInitialVertex");
  FundDomainVertex<T,Tint> eVert = get_initial_vertex<T,Tint>(G, OptionInitialVertex, FileInitialVertex);
  //
  ResultEdgewalk<T,Tint> re = LORENTZ_RunEdgewalkAlgorithm<T,Tint,Tgroup>(G, l_norms, eVert);
  std::string OutFormat=BlockPROC.ListStringValues.at("OutFormat");
  std::string FileOut=BlockPROC.ListStringValues.at("FileOut");
  if (FileOut == "unset") {
    PrintResultEdgewalk(G, re, std::cerr, OutFormat);
  } else {
    std::ofstream os(FileOut);
    PrintResultEdgewalk(G, re, os, OutFormat);
  }

}



#endif

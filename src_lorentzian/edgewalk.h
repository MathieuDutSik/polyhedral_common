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
  ListStringValues1["FileLorMat"]="the lorentzian matrix used";
  ListStringValues1["OptionInitialVertex"]="vinberg or File and if File selected use FileVertDomain as initial vertex";
  ListStringValues1["FileInitialVertex"]="unset put the name of the file used for the initial vertex";
  ListStringValues1["OptionNorms"]="possible option K3 (then just 2) or all where all norms are considered";
  ListStringValues1["OutFormat"]="GAP for gap use or TXT for text output";
  ListStringValues1["FileOut"]="stdout, or stderr or the filename you want to write to";
  SingleBlock BlockPROC;
  BlockPROC.ListStringValues=ListStringValues1;
  ListBlock["PROC"]=BlockPROC;
  // Merging all data
  return {ListBlock, "undefined"};
}



template<typename T>
size_t GetMatrixExponentSublattice(MyMatrix<T> const& g, MyMatrix<T> const& Latt)
{
  int n = Latt.rows();
  auto is_preserving=[&](MyMatrix<T> const& h) -> bool {
    for (int i=0; i<n; i++) {
      MyVector<T> eV = GetMatrixRow(Latt, i);
      MyVector<T> eVimg = h.transpose() * eV;
      std::optional<MyVector<T>> opt = SolutionIntMat(Latt, eVimg);
      if (!opt)
        return false;
    }
    return true;
  };
  size_t ord = 1;
  MyMatrix<T> h = g;
  while(true) {
    if (is_preserving(h))
      break;
    ord++;
    h = h * g;
  }
  return ord;
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
  FundDomainVertex<T,Tint> fund_v;
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
  if (p1.first != p2.first)
    return get_sign_pair_t(p1.first, p2.first);
  if (p1.first == 1)
    return get_sign_pair_t(p1.second, p2.second);
  if (p1.first == -1)
    return -get_sign_pair_t(p1.second, p2.second);
  return 0;
}


template<typename T, typename Tint>
RootCandidate<T,Tint> gen_possible_extension(MyMatrix<T> const& G, MyVector<T> const& k, MyVector<Tint> const& alpha, T const& res_norm, T const& e_norm, FundDomainVertex<T,Tint> const& fund_v)
{
  MyVector<T> alpha_T = UniversalVectorConversion<T,Tint>(alpha);
  T scal = - k.dot(G * alpha_T);
  T quant1 = (scal * scal) / res_norm;
  T quant2 = (scal * scal) / e_norm;
  return {get_sign_sing(scal), quant1, quant2, e_norm, alpha, fund_v};
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


template<typename T, typename Tint>
struct RootCandidateCuspidal {
  int sign; // 0 for 0, 1 for positive, -1 for negative
  T quant; // this is (kP.v_{N,\Delta'})^2 / N
  T e_norm;
  MyVector<Tint> v;
};


template<typename T, typename Tint>
RootCandidateCuspidal<T,Tint> gen_possible_cuspidalextension(MyMatrix<T> const& G, MyVector<T> const& kP, MyVector<T> const& v_T, T const& e_norm)
{
  MyVector<Tint> v = UniversalVectorConversion<Tint,T>(v_T);
  T scal = - kP.dot(G * v_T);
  T quant = (scal * scal) / e_norm;
  return {get_sign_sing(scal), quant, e_norm, v};
}



/*
  We are looking forthe smallest solution c>0 for the equation
  u + c k in Latt
  By selecting an adequate basis for Latt we can reduce the problem to
  u + ck = a1 v1 + a2 v2
  with a1, a2 in Z and v1, v2 in Latt.
  We write u = u1 v1 + x2 u2   and   k = k1 v1 + k2 v2
  This gets us
  u1 + ck1 = a1
  u2 + ck2 = a2
  The right way to solve the equation is to compute kG = gcd(k1, k2) and a basis of the kernel.
  We thus remap the equation to
  u1 + c k1 = a1
  u2        = a2
  solvability condition becomes u2 in Z.
  c0 = -u1 / k1
  cS = 1/k1
  Solution is c = c0 + h cS
  k = - c0 / cS
 */
template<typename T>
std::optional<MyVector<T>> ResolveLattEquation(MyMatrix<T> const& Latt, MyVector<T> const& u, MyVector<T> const& k)
{
  std::cerr << "ResolveLattEquation k="; WriteVector(std::cerr, RemoveFractionVector(k));
  int n = Latt.rows();
  std::vector<MyVector<T>> l_v = {u,k};
  MyMatrix<T> eIndep = MatrixFromVectorFamily(l_v);
  MyMatrix<T> eBasis = ExtendToBasis(eIndep);
  MyMatrix<T> Latt2 = Latt * Inverse(eBasis);
  std::vector<int> V(n-2);
  for (int i=0; i<n-2; i++)
    V[i] = i+2;
  MyMatrix<T> Latt3 = SelectColumn(Latt2, V);
  //  std::cerr << "Latt3=\n";
  //  WriteMatrix(std::cerr, Latt3);
  MyMatrix<T> NSP = NullspaceIntMat(Latt3);
  if (!IsIntegralMatrix(NSP)) {
    std::cerr << "NSP should be integral\n";
    throw TerminalException{1};
  }
  //  std::cerr << "NSP=\n";
  //  WriteMatrix(std::cerr, NSP);
  MyMatrix<T> IntBasis = NSP * Latt;
  std::cerr << "IntBasis=\n";
  WriteMatrix(std::cerr, IntBasis);
  std::optional<MyVector<T>> opt_u = SolutionMat(IntBasis, u);
  if (!opt_u) {
    std::cerr << "We failed to find a solution for u\n";
    throw TerminalException{1};
  }
  MyVector<T> sol_u = *opt_u;
  T u1 = sol_u(0);
  T u2 = sol_u(1);
  std::cerr << "u1=" << u1 << " u2=" << u2 << "\n";
  std::optional<MyVector<T>> opt_k = SolutionMat(IntBasis, k);
  if (!opt_k) {
    std::cerr << "We failed to find a solution for k\n";
    throw TerminalException{1};
  }
  MyVector<T> sol_k = *opt_k;
  T k1 = sol_k(0);
  T k2 = sol_k(1);
  std::cerr << "k1=" << k1 << " k2=" << k2 << "\n";
  std::cerr << "u=";
  WriteVector(std::cerr, u);
  std::cerr << "k=";
  WriteVector(std::cerr, k);
  //
  GCD_int<T> ep = ComputePairGcd(k1, k2);
  T u1_norm = ep.Pmat(0,0) * u1 + ep.Pmat(1,0) * u2;
  T u2_norm = ep.Pmat(0,1) * u1 + ep.Pmat(1,1) * u2;
  T k1_norm = ep.Pmat(0,0) * k1 + ep.Pmat(1,0) * k2;
  T k2_norm = ep.Pmat(0,1) * k1 + ep.Pmat(1,1) * k2;
  std::cerr << "norm : u1=" << u1_norm << " u2=" << u2_norm << "\n";
  std::cerr << "norm : k1=" << k1_norm << " k2=" << k2_norm << "\n";
  if (k2_norm != 0) {
    std::cerr << "We should have k2_norm = 0. Likely a bug here\n";
    throw TerminalException{1};
  }
  if (!IsInteger(u2_norm)) // No solution then
    return {};
  //
  T c0 = - u1_norm / k1_norm;
  T cS = 1 / k1_norm;
  std::cerr << "c0=" << c0 << " cS=" << cS << "\n";
  T hinp = -c0 / cS;
  T h;
  if (cS > 0) {
    h = UniversalCeilScalarInteger<T,T>(hinp);
    std::cerr << "1 : hinp=" << hinp << " h=" << h << "\n";
    if (hinp == h)
      h += 1;
  } else {
    h = UniversalFloorScalarInteger<T,T>(hinp);
    std::cerr << "2 : hinp=" << hinp << " h=" << h << "\n";
    if (hinp == h)
      h -= 1;
  }
  T c = c0 + h * cS;
  std::cerr << "h=" << h << " c=" << c << "\n";
  if (c <= 0) {
    std::cerr << "We should have c>0\n";
    throw TerminalException{1};
  }
  return u + c * k;
}





/*
  We use solution of Problem 7.1 in the edgewalk text.
  It is indeed true that we are searching roots that contain k.
  The dimension of that space is n.
  But since k is isotropic that space is actually the span of k and the l_ui.
  --
  Now if we write a vector as u + ck then we get N(u + ck) = N(u)
 */
template<typename T, typename Tint>
std::vector<MyVector<Tint>> DetermineRootsCuspidalCase(MyMatrix<T> const& G, std::vector<MyVector<Tint>> const& l_ui, std::vector<T> const& l_norms,
                                                       MyVector<T> const& k, MyVector<T> const& kP)
{
  std::cerr << "DetermineRootsCuspidalCase, step 1\n";
  bool only_spherical = false;
  std::vector<Possible_Extension<T>> l_extension = ComputePossibleExtensions(G, l_ui, l_norms, only_spherical);
  std::cerr << "DetermineRootsCuspidalCase : |l_extension|=" << l_extension.size() << "\n";
  std::cerr << "DetermineRootsCuspidalCase, step 2\n";
  std::vector<RootCandidateCuspidal<T,Tint>> l_candidates;
  for (auto & e_extension : l_extension) {
    std::cerr << "res_norm=" << e_extension.res_norm << "\n";
    if (e_extension.res_norm == 0) {
      MyMatrix<T> Latt = ComputeLattice_LN(G, e_extension.e_norm);
      std::cerr << "We have det(Latt)=" << DeterminantMat(Latt) << " Latt=\n";
      WriteMatrix(std::cerr, Latt);
      std::optional<MyVector<T>> opt_v = ResolveLattEquation(Latt, e_extension.u_component, k);
      std::cerr << "We have opt_v\n";
      if (opt_v) {
        const MyVector<T>& v_T = *opt_v;
        std::cerr << "Proposed v_T =";
        WriteVector(std::cerr, v_T);
        RootCandidateCuspidal<T,Tint> e_cand = gen_possible_cuspidalextension<T,Tint>(G, kP, v_T, e_extension.e_norm);
        std::cerr << "We have e_cand\n";
        l_candidates.push_back(e_cand);
      }
    }
  }
  std::cerr << "DetermineRootsCuspidalCase : |l_candidates|=" << l_candidates.size() << "\n";
  std::cerr << "DetermineRootsCuspidalCase, step 3\n";
  /* std::sort is sorting from the highest to the smallest
   */
  std::sort(l_candidates.begin(), l_candidates.end(),
            [&](RootCandidateCuspidal<T,Tint> const& x, RootCandidateCuspidal<T,Tint> const& y) -> bool {
              // We want x > y if -k.alpha(x) / sqrt(Nx) < -k.alpha(y) / sqrt(Ny) or if equality if
              // Nx < Ny
              int sign = get_sign_pair_stdpair<T>({x.sign, x.quant}, {y.sign, y.quant});
              //              std::cerr << "x: (" << x.sign << "," << x.quant << ") y: (" << y.sign << "," << y.quant << ") sign=" << sign << "\n";
              if (sign != 0)
                return sign > 0; // because -k.alpha1 / sqrt(R1)    <     -k.alpha2 / sqrt(R2)   correspond to 1 in the above.
              return x.e_norm < y.e_norm;
            });
 
 for (auto & x : l_candidates) {
    std::cerr << "x : sign=" << x.sign << " quant=" << x.quant << " norm=" << x.e_norm << " v="; WriteVector(std::cerr, x.v);
  }
  std::cerr << "DetermineRootsCuspidalCase, step 4\n";
  std::vector<MyVector<Tint>> l_ui_ret = l_ui;
  auto is_approved=[&](MyVector<Tint> const& cand) -> bool {
    MyVector<T> G_cand_T = G * UniversalVectorConversion<T,Tint>(cand);
    for (auto & v : l_ui_ret) {
      MyVector<T> v_T = UniversalVectorConversion<T,Tint>(v);
      T scal = v_T.dot(G_cand_T);
      if (scal > 0)
        return false;
    }
    return true;
  };
  for (auto & eV : l_candidates) {
    MyVector<Tint> eV_i = eV.v;
    if (is_approved(eV_i)) {
      std::cerr << "Inserting eV_i="; WriteVector(std::cerr, eV_i);
      l_ui_ret.push_back(eV_i);
    }
  }
  std::cerr << "DetermineRootsCuspidalCase, step 5\n";
  return l_ui_ret;
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
FundDomainVertex<T,Tint> EdgewalkProcedure(MyMatrix<T> const& G, MyVector<T> const& k, std::vector<MyVector<Tint>> const& l_ui, std::vector<T> const& l_norms, MyVector<Tint> const& v_disc)
{
  //
  // Initial computation of linear algebra nature:
  // Find a basis (k,r0) of the plane P
  //
  MyVector<T> v_disc_t = UniversalVectorConversion<T,Tint>(v_disc);
  std::cerr << "v_disc_t="; WriteVector(std::cerr, v_disc_t);
  int n = G.rows();
  size_t n_root = l_ui.size();
  std::cerr << "n_root=" << n_root << "\n";
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
  std::cerr << "l_roots(Common)=\n";
  WriteMatrix(std::cerr, Space);
  MyVector<T> eP = G * k;
  std::cerr << "k="; WriteVectorGAP(std::cerr, k); std::cerr << "\n";
  T norm = k.dot(eP);
  std::cerr << "norm=" << norm << "\n";
  for (int i=0; i<n; i++)
    EquaB(n_root,i) = eP(i);
  //  std::cerr << "EquaB=\n";
  //  WriteMatrix(std::cerr, EquaB);
  std::cerr << "RankMat(EquaB)=" << RankMat(EquaB) << "\n";
  MyMatrix<T> NSP = NullspaceTrMat(EquaB);
  std::cerr << "Edgewalk Procedure, step 1\n";
  if (NSP.rows() != 1) {
    std::cerr << "|NSP|=" << NSP.rows() << "/" << NSP.cols() << "\n";
    std::cerr << "The dimension should be exactly 1\n";
    throw TerminalException{1};
  }
  std::cerr << "Edgewalk Procedure, step 1\n";
  MyVector<T> r0 = GetMatrixRow(NSP,0);
  MyMatrix<T> NSPbas(2,n);
  AssignMatrixRow(NSPbas, 0, k);
  AssignMatrixRow(NSPbas, 1, r0);
  std::cerr << "|NSPbas|=" << NSPbas.rows() << " / " << NSPbas.cols() << "\n";
  std::cerr << "r0="; WriteVectorGAP(std::cerr, r0); std::cerr << "\n";
  //
  // Computing the extension and the maximum norms from that.
  //
  std::cerr << "Edgewalk Procedure, step 2\n";
  std::vector<RootCandidate<T,Tint>> l_candidates;
  bool only_spherical = true;
  std::cerr << "Edgewalk Procedure, step 3\n";
  std::vector<Possible_Extension<T>> l_extension = ComputePossibleExtensions(G, l_ui, l_norms, only_spherical);
  std::cerr << "EdgewalkProcedure : |l_extension|=" << l_extension.size() << "\n";
  std::cerr << "Edgewalk Procedure, step 4\n";
  std::map<T,T> map_max_resnorm;
  for (auto & e_extension : l_extension) {
    T norm = e_extension.e_norm;
    T res_norm = e_extension.res_norm;
    map_max_resnorm[norm] = std::max(map_max_resnorm[norm], res_norm);
  }
  //
  // Determine if the plane P is isotropic and if not compute the set of test vectors
  //
  MyMatrix<T> ProjP = GetProjectionMatrix(G, NSPbas);
  struct SingComp {
    MyMatrix<T> Latt;
    MyMatrix<T> Basis_ProjP_LN;
    MyMatrix<T> Basis_P_inter_LN;
    MyMatrix<T> Gwork;
    std::vector<MyVector<T>> l_vect;
  };
  auto get_sing_comp=[&](T const& e_norm) -> SingComp {
    MyMatrix<T> Latt = ComputeLattice_LN(G, e_norm);
    MyMatrix<T> ProjFamily(n,n);
    for (int i=0; i<n; i++) {
      MyVector<T> eVect = GetMatrixRow(Latt, i);
      MyVector<T> eVectProj = ProjMat * eVect;
      AssignMatrixRow(ProjFamily, i, eVectProj);
    }
    MyMatrix<T> Basis_ProjP_LN = GetZbasis(ProjFamily);
    MyMatrix<T> Basis_P_inter_LN = IntersectionLattice_VectorSpace(Latt, NSPbas);
    MyMatrix<T> Gwork = Basis_ProjP_LN * G * Basis_ProjP_LN.transpose();
    T res_norm = map_max_resnorm[e_norm];
    std::optional<MyVector<T>> opt = SolutionMat(Basis_ProjP_LN, r0);
    if (!opt) {
      std::cerr << "Failed to resolve the SolutionMat problem\n";
      throw TerminalException{1};
    }
    MyVector<T> r0_NSP = *opt;
    MyVector<Tint> r0_work = UniversalVectorConversion<Tint,T>(RemoveFractionVector(r0_NSP));
    MyVector<Tint> l_A = GetTwoComplement(r0_work);
    MyVector<Tint> l_B = Canonical(Gwork, res_norm, r0_work, l_A);
    std::optional<std::pair<MyMatrix<Tint>,std::vector<MyVector<Tint>>>> opt = Anisotropic(Gwork, res_norm, r0, l_B);
    if (!opt) {
      return {Latt, Basis_ProjP_LN, Basis_P_inter_LN, Gwork, {}};
    }
    const std::vector<MyVector<Tint>>& l_vect1 = opt->second;
    const MyMatrix<Tint>& P = opt->first;
    MyMatrix<T> Expr_t = Basis_P_inter_LN * Inverse(Basis_ProjP_LN);
    if (!IsIntegerMatrix(Expr_t)) {
      std::cerr << "The matrix should be integral\n";
      throw TerminalException{1};
    }
    MyMatrix<Tint> Expr_i = UniversalMatrixConversion<Tint,T>(Expr_t);
    size_t order = GetMatrixExponentSublattice(P, Expr_i);
    std::vector<MyMatrix<Tint>> l_vect2;
    for (auto & e_vect1 : l_vect1) {
      T norm1 = eval_norm(Gwork, e_vect1);
      size_t ord = 1;
      while(true) {
        T norm2 = ord * ord * norm1;
        if (norm2 > res_norm)
          break;
        l_vect2.push_back(ord * e_vect1);
        ord++;
      }
    }
    std::vector<MyMatrix<Tint>> l_vect3;
    MyMatrix<Tint> TheMat = IdentityMat<Tint>(2);
    for (size_t i=0; i<order; i++) {
      for (auto & e_vect2 : l_vect2) {
        MyVector<Tint> e_vect3 = TheMat.transpose() * e_vect2;
        l_vect3.push_back(e_vect3);
      }
      TheMat = TheMat * P;
    }
    return {Latt, Basis_ProjP_LN, Basis_P_inter_LN, Gwork, l_vect3};
  };
  using Tanisotropic = std::map<T, SingComp>;
  auto compute_anisotropic_data=[&]() -> std::optional<Tanisotropic> {
    MyMatrix<T> Gtest = NSPbas * G * NSPbas.transpose();
    std::optional<MyMatrix<T>> opt = GetIsotropicFactorization(Gtest);
    if (opt)
      return {};
    Tanisotropic map;
    for (auto & u_norm : l_norms)
      map[u_norm] = get_sing_comp(u_norm);
    return map;
  };
  std::optional<Tanisotropic> opt_anisotrop = compute_anisotropic_data();
  auto get_next_anisotropic=[&](Possible_extension<T> const& poss) -> std::optional<MyVector<Tint>> {
    T const& e_norm = poss.e_norm;
    SingComp const& e_comp = (*opt_anisotrop)[e_norm];
    for (auto & e_vect : e_comp.l_vect) {
      T val = eval_quad(e_comp, e_vect);
      if (val == poss.res_norm) {
        MyVector<T> v_T = poss.u_component + e_comp.Basis_ProjP_LN * e_vect;
        if (IsIntegerVector(v_T)) {
          std::optional<MyVector<T>> eSol = SolutionIntMat(e_comp.Latt, v_T);
          if (!eSol) {
            MyVector<Tint> v_i = UniversalVectorConversion<Tint,T>(v_T);
            return v_i;
          }
        }
      }
    }
    return false;
  };



  
  auto get_next_isotropic=[&](Possible_extension<T> const& poss) -> std::optional<MyVector<Tint>> {
    
  };

  
  //
  //
  //
  for (auto & e_extension : l_extension) {
    std::cerr << "------ u_component=" << StringVectorGAP(e_extension.u_component) << " ----------\n";
    T e_norm = e_extension.e_norm;
    MyMatrix<T> Latt = ComputeLattice_LN(G, e_norm);
    MyMatrix<T> LattInv = Inverse(Latt);
    MyMatrix<T> Space_LN = Space * LattInv;
    MyMatrix<T> G_LN = Latt * G * Latt.transpose();
    MyMatrix<T> Equas = Space_LN * G_LN;
    MyMatrix<T> NSP = NullspaceIntTrMat(Equas);
    MyMatrix<T> GP_LN = NSP * G_LN * NSP.transpose();
    MyVector<T> r0_LN = LattInv.transpose() * r0;
    std::optional<MyVector<T>> opt = SolutionMat(NSP, r0_LN);
    if (!opt) {
      std::cerr << "Failed to resolve the SolutionMat problem\n";
      throw TerminalException{1};
    }
    MyVector<T> r0_NSP = *opt;
    MyVector<Tint> r0_work = UniversalVectorConversion<Tint,T>(RemoveFractionVector(r0_NSP));
    std::optional<MyVector<Tint>> opt_v = get_first_next_vector(GP_LN, r0_work, e_extension.res_norm);
    std::cerr << "We have opt_v\n";
    if (opt_v) {
      MyVector<T> v = UniversalVectorConversion<T,Tint>(*opt_v);
      MyVector<T> r_component = Latt.transpose() * NSP.transpose() * v;
      std::optional<MyVector<T>> opt_s = SolutionMat(NSPbas, r_component);
      if (!opt_s) {
        std::cerr << "The r_component should belong to the space spanned by k and r0\n";
        throw TerminalException{1};
      }
      MyVector<T> alpha_T = e_extension.u_component + r_component;
      MyVector<Tint> alpha = UniversalVectorConversion<Tint,T>(alpha_T);
      std::cerr << "alpha="; WriteVector(std::cerr, alpha);
      auto f_ins=[&]() -> void {
        std::vector<MyVector<Tint>> l_roots = l_ui;
        l_roots.push_back(alpha);
        MyMatrix<T> Mat_root = UniversalMatrixConversion<T,Tint>(MatrixFromVectorFamily(l_roots));
        MyMatrix<T> EquaMat = Mat_root * G;
        MyMatrix<T> NSP = NullspaceTrMat(EquaMat);
        if (NSP.rows() != 1) {
          std::cerr << "We should have exactly one row\n";
          return;
        }
        MyVector<T> gen = GetMatrixRow(NSP, 0);
        T scal = gen.dot(G * k);
        auto get_gen=[&]() -> std::optional<MyVector<T>> {
          if (scal < 0) { // The sign convention means that two vectors in the same cone have negative scalar product.
            return gen;
          }
          if (scal > 0) {
            return -gen;
          }
          return {};
        };
        std::optional<MyVector<T>> opt_k_new = get_gen();
        if (opt_k_new) {
          const MyVector<T> & k_new = *opt_k_new;
          T scal = v_disc_t.dot(G * k_new);
          if (scal < 0) { // The convention in Lorentzian is negative scalar (see end of Sect 2 of edgewalk paper)
            FundDomainVertex<T,Tint> fund_v{k_new,l_roots};
            RootCandidate<T,Tint> eCand = gen_possible_extension(G, k, alpha, e_extension.res_norm, e_norm, fund_v);
            l_candidates.push_back(eCand);
          }
        }
      };
      f_ins();
    }
  }
  std::cerr << "EdgewalkProcedure : |l_candidates|=" << l_candidates.size() << "\n";
  if (l_candidates.size() > 0) {
    RootCandidate<T,Tint> best_cand = get_best_candidate(l_candidates);
    return best_cand.fund_v;
  }
  // So, no candidates were found. We need to find isotropic vectors.
  const MyMatrix<T> Gred = NSPbas * G * NSPbas.transpose();
  std::cerr << "We have Gred=\n";
  WriteMatrix(std::cerr, Gred);
  std::optional<MyMatrix<T>> Factor_opt = GetIsotropicFactorization(Gred);
  std::cerr << "We have Factor_opt\n";
  if (!Factor_opt) {
    std::cerr << "The matrix is not isotropic. Major rethink are needed\n";
    throw TerminalException{1};
  }
  MyMatrix<T> Factor = *Factor_opt;
  std::cerr << "We have Factor=\n";
  WriteMatrix(std::cerr, Factor);
  // We want a vector inside of the cone (there are two: C and -C)
  auto get_can_gen=[&](MyVector<T> const& v) -> MyVector<T> {
    T scal = k.dot(G * v);
    std::cerr << "  scal=" << scal << "\n";
    if (scal < 0) // The value should be negative because with the chosen convention, interior vectors have negative pairwise scalar products
      return v;
    if (scal > 0)
      return -v;
    std::cerr << "We should have scal != 0 to be able to conclude\n";
    throw TerminalException{1};
  };
  std::vector<MyVector<T>> l_gens;
  for (size_t i=0; i<2; i++) {
    std::cerr << "i=" << i << "\n";
    // a x + by correspond to the ray (u0, u1) = (-b, a)
    T u0 = -Factor(i,1);
    T u1 =  Factor(i,0);
    std::cerr << "u0=" << u0 << " u1=" << u1 << "\n";
    T sum = u0 * u0 * Gred(0,0) + 2 * u0 * u1 * Gred(0,1) + u1 * u1 * Gred(1,1);
    std::cerr << "sum=" << sum << "\n";
    MyVector<T> gen = u0 * k + u1 * r0;
    //    std::cerr << "k="; WriteVectorGAP(std::cerr, k); std::cerr << "\n";
    //    std::cerr << "r0="; WriteVectorGAP(std::cerr, r0); std::cerr << "\n";
    //    std::cerr << "gen="; WriteVectorGAP(std::cerr, gen); std::cerr << "\n";
    T sum_B = gen.dot(G * gen);
    std::cerr << "sum_B=" << sum_B << "\n";
    std::cerr << "gen="; WriteVector(std::cerr, gen);
    MyVector<T> can_gen = get_can_gen(gen);
    std::cerr << "can_gen="; WriteVector(std::cerr, can_gen);
    std::cerr << "RemoveFraction(can_gen)="; WriteVector(std::cerr, RemoveFractionVector(can_gen));
    T scal = v_disc_t.dot(G * can_gen);
    std::cerr << "scal=" << scal << "\n";
    if (scal < 0) // Convention is negative scalar in Lorentzian theory (see end of sect 2 of edgewalk paper)
      l_gens.push_back(can_gen);
  }
  std::cerr << "|l_gens|=" << l_gens.size() << "\n";
  if (l_gens.size() != 1) {
    std::cerr << "We should have just one vector in order to conclude. Rethink needed\n";
    throw TerminalException{1};
  }
  const MyVector<T> & k_new = l_gens[0];
  std::cerr << "k_new="; WriteVector(std::cerr, RemoveFractionVector(k_new));
  /* FILL OUT THE CODE
     What we are looking for is the roots satisfying say A[x] = 2
     and x.v = 0.
     This can only be infinite. Since otherwise, there would only be
     a finite set of roots incident to the vector.
     ----
     But we could add the constraint of negative scalar product with the known
     roots.
     This just could work. But is would be expensive, as the problem is N-dimensional.
     Most likely, it would just introduce Vinberg back and break our speed
     improvements.

   */
  std::vector<MyVector<Tint>> l_roots_ret = DetermineRootsCuspidalCase(G, l_ui, l_norms, k_new, k);
  return {k_new, l_roots_ret};
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
      std::cerr <<  "Before LinPolytopeWMat_Isomorphism\n";
      std::optional<MyMatrix<T>> equiv_opt = LinPolytopeWMat_Isomorphism<T,Tgroup,T,uint16_t>(u_pair.val.pair_char, v_pair.val.pair_char);
      std::cerr <<  "After  LinPolytopeWMat_Isomorphism\n";
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
  std::cerr << "|l_roots| len=" << len << "\n";
  auto insert_edges_from_vertex=[&](FundDomainVertex<T,Tint> const& theVert) -> void {
    for (size_t i=1; i<len; i++) {
      std::cerr << "ADJ i=" << i << "/" << len << "\n";
      std::vector<MyVector<Tint>> l_ui;
      for (size_t j=0; j<len; j++) {
        if (i != j) {
          l_ui.push_back(theVert.l_roots[j]);
        }
      }
      MyVector<Tint> v_disc = theVert.l_roots[i];
      FundDomainVertex<T,Tint> fVert = EdgewalkProcedure(G, theVert.gen, l_ui, l_norms, v_disc);
      T norm = fVert.gen.dot(G * fVert.gen);
      std::cerr << "Result of EdgewalkProcedure norm=" << norm << "\n";
      std::cerr << "We have fVert=";
      WriteVector(std::cerr, fVert.gen);
      std::cerr << "l_roots=\n";
      WriteMatrix(std::cerr, MatrixFromVectorFamily(fVert.l_roots));

      if (i == 1) {
        std::cerr << "ENDING FOR THAT DEBUG PART i=" << i << "\n";
        throw TerminalException{1};
      }
      PairVertices<T,Tint> epair = gen_pair_vertices(G, theVert, fVert);
      std::cerr << "We have epair\n";
      EnumEntry entry{true, false, std::move(epair)};
      std::cerr << "We have entry\n";
      func_insert_pair_vertices(entry);
    }
  };
  insert_edges_from_vertex(eVert);
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
  std::string OptionNorms=BlockPROC.ListStringValues.at("OptionNorms");
  std::vector<T> l_norms = get_initial_list_norms<T,Tint>(G, OptionNorms);
  //
  std::string OptionInitialVertex=BlockPROC.ListStringValues.at("OptionInitialVertex");
  std::string FileInitialVertex=BlockPROC.ListStringValues.at("FileInitialVertex");
  FundDomainVertex<T,Tint> eVert = get_initial_vertex<T,Tint>(G, OptionInitialVertex, FileInitialVertex);
  //
  ResultEdgewalk<T,Tint> re = LORENTZ_RunEdgewalkAlgorithm<T,Tint,Tgroup>(G, l_norms, eVert);
  std::string OutFormat=BlockPROC.ListStringValues.at("OutFormat");
  std::string FileOut=BlockPROC.ListStringValues.at("FileOut");
  if (FileOut == "stderr") {
    PrintResultEdgewalk(G, re, std::cerr, OutFormat);
  } else {
    if (FileOut == "stdout") {
      PrintResultEdgewalk(G, re, std::cout, OutFormat);
    } else {
      std::ofstream os(FileOut);
      PrintResultEdgewalk(G, re, os, OutFormat);
    }
  }

}



#endif

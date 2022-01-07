#ifndef INCLUDE_EDGEWALK_H
#define INCLUDE_EDGEWALK_H

#include "MatrixCanonicalForm.h"
#include "Temp_PolytopeEquiStab.h"
#include "two_dim_lorentzian.h"
#include "coxeter_dynkin.h"
#include "vinberg.h"
#include "Namelist.h"
#include "Temp_Positivity.h"
#include "POLY_lrslib.h"


#define ALLOW_VINBERG_ALGORITHM_FOR_INITIAL_VERTEX



FullNamelist NAMELIST_GetStandard_EDGEWALK()
{
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListStringValues1["FileLorMat"] = "the lorentzian matrix used";
  ListStringValues1["OptionInitialVertex"] = "vinberg or File and if File selected use FileVertDomain as initial vertex";
  ListStringValues1["FileInitialVertex"] = "unset put the name of the file used for the initial vertex";
  ListStringValues1["OptionNorms"] = "possible option K3 (then just 2) or all where all norms are considered";
  ListStringValues1["OutFormat"] = "GAP for gap use or TXT for text output";
  ListStringValues1["FileOut"] = "stdout, or stderr or the filename of the file you want to write to";
  ListBoolValues1["ComputeAllSimpleRoots"]=true;
  SingleBlock BlockPROC;
  BlockPROC.ListStringValues = ListStringValues1;
  BlockPROC.ListBoolValues = ListBoolValues1;
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
  MyMatrix<Tint> MatRoot;
};



/*
  Direct copy. We have some memory hack
 */
template<typename T, typename Tint>
FundDomainVertex<T,Tint> DirectExplicitCopyHack(FundDomainVertex<T,Tint> const& x)
{
  return {UniversalVectorConversion<T,T>(x.gen), UniversalMatrixConversion<Tint,Tint>(x.MatRoot)};
}

template<typename T, typename Tint>
void WriteFundDomainVertex(MyMatrix<T> const& G, FundDomainVertex<T,Tint> const& vert, std::ostream & os, std::string const& OutFormat)
{
  const MyMatrix<Tint> & Mroot = vert.MatRoot;
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
  std::cerr << "Failed to find a matching entry for WritePairVertices\n";
  std::cerr << "OutFormat=" << OutFormat << " but allowed values are GAP and TXT\n";
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
  T scal = - alpha_T.dot(G * k);
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
  //  std::cerr << "First best_cand=\n";
  //  WriteMatrix(std::cerr, best_cand.fund_v.MatRoot);
  for (size_t i=1; i<l_cand.size(); i++) {
    //    std::cerr << "i=" << i << "\n";
    //    std::cerr << "Considering l_cand[i]=\n";
    //    WriteMatrix(std::cerr, l_cand[i].fund_v.MatRoot);;
    if (get_sign_cand(l_cand[i], best_cand) == 1) {
      best_cand = l_cand[i];
      //      std::cerr << "Now best_cand=\n";
      //      WriteMatrix(std::cerr, best_cand.fund_v.MatRoot);
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
  std::vector<MyVector<T>> l_v = {u,k};
  MyMatrix<T> eIndep = MatrixFromVectorFamily(l_v);
  MyMatrix<T> IntBasis = IntersectionLattice_VectorSpace(Latt, eIndep);
  //  std::cerr << "IntBasis=\n";
  //  WriteMatrix(std::cerr, IntBasis);
  std::optional<MyVector<T>> opt_u = SolutionMat(IntBasis, u);
  if (!opt_u) {
    std::cerr << "We failed to find a solution for u\n";
    throw TerminalException{1};
  }
  MyVector<T> sol_u = *opt_u;
  T u1 = sol_u(0);
  T u2 = sol_u(1);
  //  std::cerr << "u1=" << u1 << " u2=" << u2 << "\n";
  std::optional<MyVector<T>> opt_k = SolutionMat(IntBasis, k);
  if (!opt_k) {
    std::cerr << "We failed to find a solution for k\n";
    throw TerminalException{1};
  }
  MyVector<T> sol_k = *opt_k;
  T k1 = sol_k(0);
  T k2 = sol_k(1);
  //  std::cerr << "k1=" << k1 << " k2=" << k2 << "\n";
  //  std::cerr << "u=";
  //  WriteVector(std::cerr, u);
  //  std::cerr << "k=";
  //  WriteVector(std::cerr, k);
  //
  GCD_int<T> ep = ComputePairGcd(k1, k2);
  T u1_norm = ep.Pmat(0,0) * u1 + ep.Pmat(1,0) * u2;
  T u2_norm = ep.Pmat(0,1) * u1 + ep.Pmat(1,1) * u2;
  T k1_norm = ep.Pmat(0,0) * k1 + ep.Pmat(1,0) * k2;
  T k2_norm = ep.Pmat(0,1) * k1 + ep.Pmat(1,1) * k2;
  //  std::cerr << "norm : u1=" << u1_norm << " u2=" << u2_norm << "\n";
  //  std::cerr << "norm : k1=" << k1_norm << " k2=" << k2_norm << "\n";
  if (k2_norm != 0) {
    std::cerr << "We should have k2_norm = 0. Likely a bug here\n";
    throw TerminalException{1};
  }
  if (!IsInteger(u2_norm)) // No solution then
    return {};
  //
  T c0 = - u1_norm / k1_norm;
  T cS = 1 / k1_norm;
  //  std::cerr << "c0=" << c0 << " cS=" << cS << "\n";
  T hinp = -c0 / cS;
  T h;
  if (cS > 0) {
    h = UniversalCeilScalarInteger<T,T>(hinp);
    //    std::cerr << "1 : hinp=" << hinp << " h=" << h << "\n";
    if (hinp == h)
      h += 1;
  } else {
    h = UniversalFloorScalarInteger<T,T>(hinp);
    //    std::cerr << "2 : hinp=" << hinp << " h=" << h << "\n";
    if (hinp == h)
      h -= 1;
  }
  T c = c0 + h * cS;
  //  std::cerr << "h=" << h << " c=" << c << "\n";
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
  std::cerr << "-------------------------------------- EDGEWALK PROCEDURE ---------------------------------------------\n";
  std::cerr << "k=" << StringVectorGAP(k) << "\n";
  std::cerr << "l_norms =";
  for (auto & eN : l_norms)
    std::cerr << " " << eN;
  std::cerr << "\n";
  std::cerr << "l_ui =";
  for (auto & eV : l_ui)
    std::cerr << " " << StringVectorGAP(eV);
  std::cerr << "\n";
  std::cerr << "v_disc=" << StringVectorGAP(v_disc) << "\n";
  std::cerr << "    Real work starts now\n";
  //
  // Initial computation of linear algebra nature:
  // Find a basis (k,r0) of the plane P
  //
  MyVector<T> v_disc_t = UniversalVectorConversion<T,Tint>(v_disc);
  /*
  MyMatrix<T> CoxMat = ComputeCoxeterMatrix(G, l_ui).first;
  std::string symb = coxdyn_matrix_to_string(CoxMat);
  std::cerr << "Coxeter diagram of the vertex k in u_i direction=" << symb << "\n";
  */
  int n = G.rows();
  size_t n_root = l_ui.size();
  std::cerr << "n_root=" << n_root << "\n";
  MyMatrix<T> Space(n_root,n);
  MyMatrix<T> EquaRvect(n_root+1,n);
  MyMatrix<T> EquaPplane(n_root,n);
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
      EquaRvect(i_root,i) = eP(i);
      EquaPplane(i_root,i) = eP(i);
    }
  }
  MyMatrix<T> Pplane = NullspaceTrMat(EquaPplane);
  std::cerr << "Plane P=\n";
  WriteMatrix(std::cerr, Pplane);
  //  std::cerr << "Space=\n";
  //  WriteMatrix(std::cerr, Space);
  MyVector<T> eP = G * k;
  T norm = k.dot(eP);
  std::cerr << "k=" << StringVectorGAP(k) << " norm=" << norm << "\n";
  for (int i=0; i<n; i++)
    EquaRvect(n_root,i) = eP(i);
  MyMatrix<T> NSP = NullspaceTrMat(EquaRvect);
  std::cerr << "Edgewalk Procedure, step 1\n";
  if (NSP.rows() != 1) {
    std::cerr << "|NSP|=" << NSP.rows() << "/" << NSP.cols() << "\n";
    std::cerr << "The dimension should be exactly 1\n";
    throw TerminalException{1};
  }
  std::cerr << "Edgewalk Procedure, step 1\n";
  /*
    The vector r0 is orthogonal to k and is well defined up to a sign.
    The half plane (1/2)P is shown on Figure 8.1 as being orthogonal to
    k. Note that the mention on page 26 "(1/)2P is the open half plane corresponding to e"
    is not correct. It should be corresponding to k. But how to build it?
    How to select the sign?
   */
  MyVector<T> r0 = GetMatrixRow(NSP,0);
  T scal_r0 = r0.dot(G * v_disc_t);
  std::cerr << "scal_r0=" << scal_r0 << "\n";
  if (scal_r0 > 0)
    r0 = -r0;
  // We follow here the convention on oriented basis of Section 8:
  // "First member lies in the interior of (1/2)P and whose second member is k"
  MyMatrix<T> OrientedBasis(2,n);
  if (norm < 0) { // the point is inner, the oriented basis is clear.
    std::cerr << "Builging OrientedBasis, ordinary case\n";
    AssignMatrixRow(OrientedBasis, 0, r0);
    AssignMatrixRow(OrientedBasis, 1, k);
  } else { // Now the point is ideal
    std::cerr << "Builging OrientedBasis, ideal case\n";
    /*
      The roots to be found are of positive norms.
      In general, we cannot think of which zone of positive norms to select from since
      those are connected.
      However, in the specific case of dimension 2, yes the vectors of positive norms
      are in two connected components.
      How to select which connected component?
      The scalar product with k allows us to select which one
     */
    auto get_positive_norm_vector=[&]() -> MyVector<T> {
      T two = 2;
      T alpha = 1;
      while(true) {
        for (int i_bas=0; i_bas<2; i_bas++) {
          MyVector<T> v_bas = GetMatrixRow(Pplane, i_bas);
          for (int u=-1; u<2; u += 2) {
            MyVector<T> v_pos_cand = k + u * alpha * v_bas;
            T norm = v_pos_cand.dot(G * v_pos_cand);
            std::cerr << "u=" << u << " alpha=" << alpha << " v_pos_cand=" << v_pos_cand << " norm=" << norm << "\n";
            if (norm > 0)
              return v_pos_cand;
          }
        }
        alpha /= two;
      }
    };
    MyVector<T> v_pos = get_positive_norm_vector();
    T scal = v_pos.dot(G * k);
    if (scal > 0) { // The convention is that negative scalar product is for facets.
      v_pos = -v_pos;
    }
    AssignMatrixRow(OrientedBasis, 0, v_pos);
    AssignMatrixRow(OrientedBasis, 1, k);
    r0 = -k; // Follows Right part of Figure 8.1
  }
  std::cerr << "r0=" << StringVectorGAP(r0) << "\n";
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
  // For each norm, there is a corresponding maximum possible norms and specific enumeration.
  for (auto & e_extension : l_extension) {
    T norm = e_extension.e_norm;
    T res_norm = e_extension.res_norm;
    map_max_resnorm[norm] = std::max(map_max_resnorm[norm], res_norm);
  }
  for (auto & kv : map_max_resnorm)
    std::cerr << "kv : norm=" << kv.first << " max(res_norm)=" << kv.second << "\n";
  std::cerr << "Edgewalk Procedure, step 5\n";
  //
  // Determine if the plane P is isotropic and if not compute the set of test vectors
  //
  MyMatrix<T> G_Pplane = Pplane * G * Pplane.transpose();
  std::cerr << "G_Pplane=\n";
  WriteMatrix(std::cerr, G_Pplane);
  std::cerr << "G=\n";
  WriteMatrix(std::cerr, G);
  MyMatrix<T> ProjP = GetProjectionMatrix(G, Pplane);
  std::cerr << "Edgewalk Procedure, step 6\n";
  struct SingCompAnisotropic {
    MyMatrix<T> Latt;
    MyMatrix<T> Basis_ProjP_LN;
    MyMatrix<T> Basis_P_inter_LN;
    MyMatrix<T> Gwork;
    std::vector<MyVector<Tint>> l_vect;
  };
  struct SingCompIsotropic {
    MyMatrix<T> Latt;
    MyMatrix<T> Basis_ProjP_LN;
    MyMatrix<T> GP_LN;
    MyMatrix<T> Factor_GP_LN;
    MyVector<Tint> r0_work;
    std::map<T,std::optional<std::vector<MyVector<Tint>>>> map_res_norm;
  };
  auto get_basis_projp_ln=[&](MyMatrix<T> const& Latt) -> MyMatrix<T> {
    MyMatrix<T> ProjFamily(n,n);
    for (int i=0; i<n; i++) {
      MyVector<T> eVect = GetMatrixRow(Latt, i);
      MyVector<T> eVectProj = ProjP * eVect;
      AssignMatrixRow(ProjFamily, i, eVectProj);
    }
    MyMatrix<T> BasisProj = GetZbasis(ProjFamily);
    if (BasisProj.rows() != 2) {
      std::cerr << "The BasisProj should be of rank 2\n";
      throw TerminalException{1};
    }
    MyMatrix<T> Expr = ExpressVectorsInIndependentFamilt(BasisProj, OrientedBasis);
    std::cerr << "Det(Expr)=" << DeterminantMat(Expr) << "\n";
    if (DeterminantMat(Expr) < 0) { // Change to get positive determinant
      for (int i=0; i<n; i++)
        BasisProj(0,i) = -BasisProj(0,i);
    }
    MyVector<T> Basis0 = GetMatrixRow(BasisProj, 0);
    std::cerr << "Basis0=" << StringVectorGAP(Basis0) << "\n";
    MyVector<T> Basis1 = GetMatrixRow(BasisProj, 1);
    std::cerr << "Basis1=" << StringVectorGAP(Basis1) << "\n";
    return BasisProj;
  };
  auto get_r0work=[&](MyMatrix<T> const& Basis_ProjP_LN, MyVector<T> const& r0) -> MyVector<Tint> {
    std::optional<MyVector<T>> opt_r0 = SolutionMat(Basis_ProjP_LN, r0);
    if (!opt_r0) {
      std::cerr << "Failed to resolve the SolutionMat problem\n";
      throw TerminalException{1};
    }
    MyVector<T> r0_NSP = *opt_r0;
    MyVector<Tint> r0_work = UniversalVectorConversion<Tint,T>(RemoveFractionVector(r0_NSP));
    std::cerr << "r0_work=" << StringVectorGAP(r0_work) << "\n";
    return r0_work;
  };
  auto get_sing_comp_anisotropic=[&](T const& e_norm) -> SingCompAnisotropic {
    std::cerr << "get_sing_comp_anisotropic, e_norm=" << e_norm << "\n";
    MyMatrix<T> Latt = ComputeLattice_LN(G, e_norm);
    MyMatrix<T> Basis_ProjP_LN = get_basis_projp_ln(Latt);
    std::cerr << "Basis_ProjP_LN=\n";
    WriteMatrix(std::cerr, Basis_ProjP_LN);
    MyMatrix<T> Basis_P_inter_LN = IntersectionLattice_VectorSpace(Latt, Pplane);
    MyMatrix<T> Gwork = Basis_ProjP_LN * G * Basis_ProjP_LN.transpose();
    T res_norm = map_max_resnorm[e_norm];
    MyVector<Tint> r0_work = get_r0work(Basis_ProjP_LN, r0);
    T r0_norm = eval_quad(Gwork, r0_work);
    std::cerr << "r0_norm=" << r0_norm << "\n";
    MyVector<Tint> l_A = GetTwoComplement(r0_work);
    std::cerr << "l_A=" << StringVectorGAP(l_A) << " res_norm=" << res_norm << "\n";
    MyVector<Tint> l_B = Canonical(Gwork, r0_norm, r0_work, l_A);
    std::cerr << "l_B=" << StringVectorGAP(l_B) << "\n";
    std::cerr << "get_sing_comp_anisotropic, step 2\n";
    std::optional<std::pair<MyMatrix<Tint>,std::vector<MyVector<Tint>>>> opt = Anisotropic<T,Tint>(Gwork, res_norm, r0_work, l_B);
    std::cerr << "get_sing_comp_anisotropic, step 3\n";
    if (!opt) { // No solution, this definitely can happen
      return {Latt, Basis_ProjP_LN, Basis_P_inter_LN, Gwork, {}};
    }
    const std::vector<MyVector<Tint>>& l_vect1 = opt->second;
    std::cerr << "|l_vect1|=" << l_vect1.size() << "\n";
    const MyMatrix<Tint>& P = opt->first;
    std::cerr << "Basis_ProjP_LN=\n";
    WriteMatrix(std::cerr, Basis_ProjP_LN);
    std::cerr << "Basis_P_inter_LN=\n";
    WriteMatrix(std::cerr, Basis_P_inter_LN);
    MyMatrix<T> Expr_t = ExpressVectorsInIndependentFamilt(Basis_P_inter_LN, Basis_ProjP_LN);
    std::cerr << "Expr_t=\n";
    WriteMatrix(std::cerr, Expr_t);
    std::cerr << "get_sing_comp_anisotropic, step 4\n";
    if (!IsIntegralMatrix(Expr_t)) {
      std::cerr << "The matrix should be integral\n";
      throw TerminalException{1};
    }
    MyMatrix<Tint> Expr_i = UniversalMatrixConversion<Tint,T>(Expr_t);
    size_t order = GetMatrixExponentSublattice(P, Expr_i);
    std::cerr << "order=" << order << "\n";
    std::vector<MyMatrix<Tint>> l_vect2;
    std::cerr << "get_sing_comp_anisotropic, step 5\n";
    for (auto & e_vect1 : l_vect1) {
      T norm1 = eval_quad(Gwork, e_vect1);
      size_t ord = 1;
      while(true) {
        T norm2 = ord * ord * norm1;
        if (norm2 > res_norm)
          break;
        MyVector<Tint> e_vect2 = ord * e_vect1;
        std::cerr << "norm2=" << norm2 << " e_vect2=" << StringVectorGAP(e_vect2) << "\n";
        l_vect2.push_back(e_vect2);
        ord++;
      }
    }
    std::cerr << "|l_vect2|=" << l_vect2.size() << "\n";
    std::vector<MyVector<Tint>> l_vect3;
    MyMatrix<Tint> TheMat = IdentityMat<Tint>(2);
    std::cerr << "get_sing_comp_anisotropic, step 6\n";
    for (size_t i=0; i<order; i++) {
      for (auto & e_vect2 : l_vect2) {
        MyVector<Tint> e_vect3 = TheMat.transpose() * e_vect2;
        l_vect3.push_back(e_vect3);
      }
      TheMat = TheMat * P;
    }
    std::cerr << "|l_vect3|=" << l_vect3.size() << "\n";
    std::cerr << "get_sing_comp_anisotropic, step 7\n";
    return {Latt, Basis_ProjP_LN, Basis_P_inter_LN, Gwork, l_vect3};
  };
  auto get_sing_comp_isotropic=[&](T const& e_norm) -> SingCompIsotropic {
    std::cerr << "get_sing_comp_isotropic, e_norm=" << e_norm << "\n";
    MyMatrix<T> Latt = ComputeLattice_LN(G, e_norm);
    MyMatrix<T> Basis_ProjP_LN = get_basis_projp_ln(Latt);
    MyMatrix<T> GP_LN = Basis_ProjP_LN * G * Basis_ProjP_LN.transpose();
    MyVector<Tint> r0_work = get_r0work(Basis_ProjP_LN, r0);
    std::optional<MyMatrix<T>> opt_factor = GetIsotropicFactorization(GP_LN);
    if (!opt_factor) {
      std::cerr << "Computation of isotropy factorization failed\n";
      throw TerminalException{1};
    }
    MyMatrix<T> Factor_GP_LN = *opt_factor;
    return {Latt, Basis_ProjP_LN, GP_LN, Factor_GP_LN, r0_work, {}};
  };
  bool is_isotropic;
  std::map<T, SingCompAnisotropic> map_anisotropic;
  std::map<T, SingCompIsotropic> map_isotropic;
  {
    MyMatrix<T> Gtest = Pplane * G * Pplane.transpose();
    std::optional<MyMatrix<T>> opt = GetIsotropicFactorization(Gtest);
    if (opt) {
      is_isotropic = true;
      for (auto & u_norm : l_norms)
        map_isotropic[u_norm] = get_sing_comp_isotropic(u_norm);
    } else {
      is_isotropic = false;
      for (auto & u_norm : l_norms)
        map_anisotropic[u_norm] = get_sing_comp_anisotropic(u_norm);
    }
  }
  std::cerr << "is_isotropic=" << is_isotropic << "\n";
  std::cerr << "Edgewalk Procedure, step 7\n";
  // Evaluation of fun
  auto get_next_anisotropic=[&](Possible_Extension<T> const& poss) -> std::optional<MyVector<Tint>> {
    T const& e_norm = poss.e_norm;
    SingCompAnisotropic const& e_comp = map_anisotropic[e_norm];
    //    std::cerr << "gna, step 1 res_norm=" << poss.res_norm << "\n";
    for (auto & e_vect : e_comp.l_vect) {
      T val = eval_quad(e_comp.Gwork, e_vect);
      //      std::cerr << "gna, step 2 e_vect=" << StringVectorGAP(e_vect) << " val=" << val << "\n";
      if (val == poss.res_norm) {
        //        std::cerr << "e_vect=" << StringVectorGAP(e_vect) << " Basis_ProjP_LN=\n";
        //        WriteMatrix(std::cerr, e_comp.Basis_ProjP_LN);
        MyVector<T> v_T = poss.u_component + e_comp.Basis_ProjP_LN.transpose() * UniversalVectorConversion<T,Tint>(e_vect);
        //        std::cerr << "gna, step 3 v_T=" << StringVectorGAP(v_T) << "\n";
        if (IsIntegerVector(v_T)) {
          std::optional<MyVector<T>> eSol = SolutionIntMat(e_comp.Latt, v_T);
          //          std::cerr << "get_next_anisotropic, step 4\n";
          if (eSol) {
            MyVector<Tint> v_i = UniversalVectorConversion<Tint,T>(v_T);
            std::cerr << "Returning v_i=" << StringVectorGAP(v_i) << "\n";
            return v_i;
          }
        }
      }
    }
    std::cerr << "No good vector found\n";
    return {};
  };
  auto get_successive_list_cand=[&](SingCompIsotropic & e_comp, T const& res_norm) -> std::vector<MyVector<Tint>> {
    std::optional<std::vector<MyVector<Tint>>> & opt = e_comp.map_res_norm[res_norm];
    if (opt)
      return *opt;
    std::vector<MyVector<Tint>> l_vect1 = EnumerateVectorFixedNorm_Factorization<T,Tint>(e_comp.Factor_GP_LN, res_norm);
    std::vector<MyVector<Tint>> l_vect2;
    auto det=[&](MyVector<Tint> const& v1, MyVector<Tint> const& v2) -> Tint {
      return v1(0) * v2(1) - v1(1) * v2(0);
    };
    auto is_corr=[&](MyVector<Tint> const& x) -> bool {
      T scal = eval_scal(G, e_comp.r0_work, x);
      if (scal <= 0)
        return false;
      return det(e_comp.r0_work, x) > 0;
    };
    for (auto & e_vect1 : l_vect1)
      if (is_corr(e_vect1))
        l_vect2.push_back(e_vect1);
    std::sort(l_vect2.begin(), l_vect2.end(),
              [&](MyVector<Tint> const& x, MyVector<Tint> const& y) -> bool {
                return det(x, y) > 0;
              });
    opt = l_vect2;
    return l_vect2;
  };
  auto get_next_isotropic=[&](Possible_Extension<T> const& poss) -> std::optional<MyVector<Tint>> {
    T const& e_norm = poss.e_norm;
    SingCompIsotropic & e_comp = map_isotropic[e_norm];
    T res_norm = poss.res_norm;
    std::vector<MyVector<Tint>> l_vect = get_successive_list_cand(e_comp, res_norm);
    std::cerr << "|l_vect|=" << l_vect.size() << "\n";
    for (auto & e_vect : l_vect) {
      std::cerr << "e_vect=" << StringVectorGAP(e_vect) << "\n";
      //      std::cerr << "Before v_T computation\n";
      //      std::cerr << "u_component=" << StringVectorGAP(poss.u_component) << "\n";
      //      std::cerr << "Basis_ProjP_LN=";
      //      WriteMatrix(std::cerr, e_comp.Basis_ProjP_LN);
      MyVector<T> v_T = poss.u_component + e_comp.Basis_ProjP_LN.transpose() * UniversalVectorConversion<T,Tint>(e_vect);
      std::cerr << "After v_T computation v_T=" << StringVectorGAP(v_T) << "\n";
      if (IsIntegerVector(v_T)) {
        std::optional<MyVector<T>> eSol = SolutionIntMat(e_comp.Latt, v_T);
        if (eSol) {
          MyVector<Tint> v_i = UniversalVectorConversion<Tint,T>(v_T);
          std::cerr << "get_next_isotropic. Returning v_i=" << StringVectorGAP(v_i) << "\n";
          return v_i;
        }
      }
    }
    return {};
  };
  auto get_next=[&](Possible_Extension<T> const& poss) -> std::optional<MyVector<Tint>> {
    if (is_isotropic) {
      return get_next_isotropic(poss);
    } else {
      return get_next_anisotropic(poss);
    }
  };
  //
  //
  //
  for (auto & e_extension : l_extension) {
    std::cerr << "------ u_component=" << StringVectorGAP(e_extension.u_component) << " ----------\n";
    T e_norm = e_extension.e_norm;


    std::optional<MyVector<Tint>> opt_v = get_next(e_extension);
    std::cerr << "We have opt_v\n";
    if (opt_v) {
      MyVector<Tint> alpha = *opt_v;
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
        std::cerr << "gen=" << StringVectorGAP(gen) << " k=" << StringVectorGAP(k) << "\n";
        T scal = gen.dot(G * k);
        std::cerr << "We have scal\n";
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
          MyVector<T> k_new = RemoveFractionVector(*opt_k_new);
          T scal = v_disc_t.dot(G * k_new);
          if (scal < 0) { // The convention in Lorentzian is negative scalar (see end of Sect 2 of edgewalk paper)
            MyMatrix<Tint> MatRoot = MatrixFromVectorFamily(l_roots);
            FundDomainVertex<T,Tint> fund_v{k_new,MatRoot};
            std::cerr << "k_new=" << StringVectorGAP(k_new) << "\n";
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
    for (auto e_cand : l_candidates)
      std::cerr << "e_cand sign=" << e_cand.sign << " quant1=" << e_cand.quant1 << " quant2=" << e_cand.quant2 << " e_norm=" << e_cand.e_norm << " fund_v=" << StringVectorGAP(e_cand.fund_v.gen) << " alpha=" << StringVectorGAP(e_cand.alpha) << "\n";
    RootCandidate<T,Tint> best_cand = get_best_candidate(l_candidates);
    std::cerr << "fund_v=" << StringVectorGAP(best_cand.fund_v.gen) << "\n";
    std::cerr << "MatRoot=\n";
    WriteMatrix(std::cerr, best_cand.fund_v.MatRoot);
    return best_cand.fund_v;
  }
  // So, no candidates were found. We need to find isotropic vectors.
  const MyMatrix<T> Gred = Pplane * G * Pplane.transpose();
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
    std::cerr << "k=" << StringVectorGAP(k) << " v=" << StringVectorGAP(v) << "\n";
    std::cerr << "We should have scal != 0 to be able to conclude\n";
    throw TerminalException{1};
  };
  std::vector<MyVector<T>> l_gens;
  for (size_t i=0; i<2; i++) {
    std::cerr << "i=" << i << "\n";
    // a x + by correspond to the ray (u0, u1) = (-b, a)
    MyVector<T> U(2);
    U(0) = -Factor(i,1);
    U(1) =  Factor(i,0);
    std::cerr << "U=" << StringVectorGAP(U) << "\n";
    T sum = U.dot(Gred * U);
    std::cerr << "sum=" << sum << "\n";
    MyVector<T> gen = Pplane.transpose() * U;
    //    std::cerr << "k="; WriteVectorGAP(std::cerr, k); std::cerr << "\n";
    //    std::cerr << "r0="; WriteVectorGAP(std::cerr, r0); std::cerr << "\n";
    //    std::cerr << "gen="; WriteVectorGAP(std::cerr, gen); std::cerr << "\n";
    T sum_B = gen.dot(G * gen);
    std::cerr << "sum_B=" << sum_B << "\n";
    std::cerr << "gen="; WriteVector(std::cerr, gen);
    if (!IsVectorMultiple(gen, k)) {
      MyVector<T> can_gen = get_can_gen(gen);
      std::cerr << "can_gen="; WriteVector(std::cerr, can_gen);
      std::cerr << "RemoveFraction(can_gen)="; WriteVector(std::cerr, RemoveFractionVector(can_gen));
      T scal = v_disc_t.dot(G * can_gen);
      std::cerr << "scal=" << scal << "\n";
      if (scal < 0) // Convention is negative scalar in Lorentzian theory (see end of sect 2 of edgewalk paper)
        l_gens.push_back(can_gen);
    }
  }
  std::cerr << "|l_gens|=" << l_gens.size() << "\n";
  if (l_gens.size() != 1) {
    std::cerr << "We should have just one vector in order to conclude. Rethink needed\n";
    throw TerminalException{1};
  }
  const MyVector<T> & k_new = l_gens[0];
  std::cerr << "k_new=" << StringVectorGAP(RemoveFractionVector(k_new)) << "\n";
  std::vector<MyVector<Tint>> l_roots_ret = DetermineRootsCuspidalCase(G, l_ui, l_norms, k_new, k);
  return {RemoveFractionVector(k_new), MatrixFromVectorFamily(l_roots_ret)};
}




template<typename T>
using pair_char = std::pair<MyMatrix<T>,WeightMatrix<true,std::vector<T>,uint16_t>>;



template<typename T, typename Tint>
void WritePairVertices(MyMatrix<T> const& G,
                       FundDomainVertex<T,Tint> const& evert1, FundDomainVertex<T,Tint> const& evert2,
                       std::ostream & os, std::string const& OutFormat)
{
  if (OutFormat == "GAP") {
    os << "rec(vert1:=";
    WriteFundDomainVertex(G, evert1, os, OutFormat);
    os << ", vert2:=";
    WriteFundDomainVertex(G, evert2, os, OutFormat);
    os << ")";
    return;
  }
  if (OutFormat == "TXT") {
    os << "vert1=\n";
    WriteFundDomainVertex(G, evert1, os, OutFormat);
    os << "vert2=\n";
    WriteFundDomainVertex(G, evert2, os, OutFormat);
    return;
  }
  std::cerr << "Failed to find a matching entry for WritePairVertices\n";
  std::cerr << "OutFormat=" << OutFormat << " but allowed values are GAP and TXT\n";
  throw TerminalException{1};
}






template<typename T, typename Tint>
pair_char<T> gen_pair_char(MyMatrix<T> const& G,
                           FundDomainVertex<T,Tint> const& evert1, FundDomainVertex<T,Tint> const& evert2)
{
  std::unordered_map<MyVector<Tint>,int> map_v;
  size_t len1 = evert1.MatRoot.rows();
  for (size_t i=0; i<len1; i++) {
    MyVector<Tint> eV = GetMatrixRow(evert1.MatRoot, i);
    map_v[eV]++;
  }
  size_t len2 = evert2.MatRoot.rows();
  for (size_t i=0; i<len2; i++) {
    MyVector<Tint> eV = GetMatrixRow(evert2.MatRoot, i);
    map_v[eV]++;
  }
  //  std::cerr << "|evert1.l_roots|=" << evert1.l_roots.size() << "  |evert2.l_roots|=" << evert2.l_roots.size() << "\n";
  std::vector<MyVector<Tint>> l_vect;
  std::vector<T> Vdiag;
  for (auto & kv : map_v) {
    l_vect.push_back(kv.first);
    Vdiag.push_back(T(kv.second));
  }
  l_vect.push_back(UniversalVectorConversion<Tint,T>(RemoveFractionVector(evert1.gen)));
  l_vect.push_back(UniversalVectorConversion<Tint,T>(RemoveFractionVector(evert2.gen)));
  T insVal = 10;
  Vdiag.push_back(insVal);
  Vdiag.push_back(insVal);
  //  std::cerr << "Vdiag=" << Vdiag << "\n";
  MyMatrix<T> MatV = UniversalMatrixConversion<T,Tint>(MatrixFromVectorFamily(l_vect));
  MyMatrix<T> Gred = MatV * G * MatV.transpose();
  //  std::cerr << "Gred=\n";
  //  WriteMatrix(std::cerr, Gred);
  std::vector<MyMatrix<T>> ListMat{G};
  using Tidx = uint32_t;
  using Tidx_value = uint16_t;
  WeightMatrix<true, std::vector<T>, Tidx_value> WMat = GetWeightMatrix_ListMat_Vdiag<T,Tidx,Tidx_value>(MatV, ListMat, Vdiag);
  WMat.ReorderingSetWeight();
  return {MatV,std::move(WMat)};
}



template<typename T, typename Tint>
struct ResultEdgewalk {
  std::vector<MyMatrix<Tint>> l_gen_isom_cox;
  std::vector<FundDomainVertex<T,Tint>> l_orbit_pair_vertices;
};


template<typename T, typename Tint>
std::vector<MyVector<Tint>> get_complete_finite_root_set(ResultEdgewalk<T,Tint> const& re)
{
  for (auto & eGen : re.l_gen_isom_cox) {
    std::cerr << "eGen=\n";
    WriteMatrix(std::cerr, eGen);
  }
  std::unordered_set<MyVector<Tint>> TotalList;
  auto f_insert=[&](MyVector<Tint> const& v) -> void {
    std::cerr << "f_insert call with v=" << StringVectorGAP(v) << "\n";
    if (TotalList.count(v) != 0)
      return;
    std::unordered_set<MyVector<Tint>> s_v;
    std::vector<MyVector<Tint>> l_v;
    auto f_ins=[&](MyVector<Tint> const& w) -> void {
      if (s_v.count(w) != 0)
        return;
      s_v.insert(w);
      l_v.push_back(w);
    };
    f_ins(v);
    size_t pos=0;
    while(true) {
      size_t len = l_v.size();
      std::cerr << "pos=" << pos << " len=" << len << "\n";
      if (pos == len)
        break;
      for (size_t i=pos; i<len; i++) {
        for (auto & eGen : re.l_gen_isom_cox) {
          MyVector<Tint> w_img = eGen.transpose() * l_v[i];
          f_ins(w_img);
        }
      }
      pos=len;
    }
  };
  for (auto & fdv : re.l_orbit_pair_vertices) {
    size_t len = fdv.MatRoot.rows();
    for (size_t i=0; i<len; i++) {
      MyVector<Tint> e_root = GetMatrixRow(fdv.MatRoot, i);
      f_insert(e_root);
    }
  }
  std::vector<MyVector<Tint>> l_root;
  for (auto & v : TotalList)
    l_root.push_back(v);
  return l_root;
}




template<typename T, typename Tint>
std::vector<T> get_list_norms(MyMatrix<T> const& G, ResultEdgewalk<T,Tint> const& re)
{
  std::set<T> set_norms;
  auto proc_vertex=[&](FundDomainVertex<T,Tint> const& vert) -> void {
    size_t len = vert.MatRoot.rows();
    for (size_t i=0; i<len; i++) {
      MyVector<Tint> root = GetMatrixRow(vert.MatRoot, i);
      MyVector<T> root_t = UniversalVectorConversion<T,Tint>(root);
      T norm = root_t.dot(G * root_t);
      set_norms.insert(norm);
    }
  };
  for (auto & evert : re.l_orbit_pair_vertices)
    proc_vertex(evert);
  std::vector<T> l_norms;
  for (auto & v : set_norms)
    l_norms.push_back(v);
  return l_norms;
}


template<typename T, typename Tint>
void PrintResultEdgewalk(MyMatrix<T> const& G, ResultEdgewalk<T,Tint> const& re, std::ostream& os, const std::string& OutFormat, bool const& ComputeAllSimpleRoots)
{
  std::vector<T> l_norms = get_list_norms(G, re);
  std::vector<MyVector<Tint>> l_simple_root;
  if (ComputeAllSimpleRoots)
    l_simple_root = get_complete_finite_root_set(re);
  if (OutFormat == "GAP") {
    os << "return rec(l_norms:=";
    WriteStdVectorGAP(os, l_norms);
    os << ", ListIsomCox:=";
    WriteVectorMatrixGAP(os, re.l_gen_isom_cox);
    os << ", ListVertices:=[";
    bool IsFirst = true;
    size_t len = re.l_orbit_pair_vertices.size() / 2;
    for (size_t i=0; i<len; i++) {
      os << "\n";
      if (!IsFirst)
        os << ",";
      IsFirst = false;
      const FundDomainVertex<T,Tint> & evert1 = re.l_orbit_pair_vertices[2*i];
      const FundDomainVertex<T,Tint> & evert2 = re.l_orbit_pair_vertices[2*i+1];
      WritePairVertices(G, evert1, evert2, os, OutFormat);
    }
    os << "]";
    if (ComputeAllSimpleRoots) {
      os << ", ListSimpleRoots:=[";
      for (size_t i=0; i<l_simple_root.size(); i++) {
        if (i>0)
          os << ",";
        os << StringVectorGAP(l_simple_root[i]);
      }
      os << "]";
    }
    os << ");\n";
    return;
  }
  if (OutFormat == "TXT") {
    os << "List of found generators of Isom / Cox\n";
    return;
  }
  std::cerr << "Failed to find a matching entry in PrintResultEdgewalk\n";
  std::cerr << "OutFormat=" << OutFormat << " but allowed values are GAP and TXT\n";
  throw TerminalException{1};
}


template<typename T, typename Tint, typename Tgroup>
ResultEdgewalk<T,Tint> LORENTZ_RunEdgewalkAlgorithm(MyMatrix<T> const& G, std::vector<T> const& l_norms, FundDomainVertex<T,Tint> const& eVert)
{
  std::unordered_set<MyMatrix<Tint>> s_gen_isom_cox;
  struct StatusEntry {
    bool stat1;
    bool stat2;
  };
  std::vector<StatusEntry> l_entry;
  std::vector<FundDomainVertex<T,Tint>> l_orbit_pair_vertices;
  MyMatrix<Tint> IdMat = IdentityMat<Tint>(G.rows());
  auto f_insert_gen=[&](MyMatrix<Tint> const& eP) -> void {
    if (eP == IdMat)
      return;
    MyMatrix<T> eP_T = UniversalMatrixConversion<T,Tint>(eP);
    MyMatrix<T> G_img = eP_T * G * eP_T.transpose();
    if (G_img != G) {
      std::cerr << "G="; WriteMatrix(std::cerr, G);
      std::cerr << "eP_T="; WriteMatrix(std::cerr, eP_T);
      std::cerr << "G_img="; WriteMatrix(std::cerr, G_img);
      std::cerr << "The matrix eP should leave the quadratic form invariant\n";
      throw TerminalException{1};
    }
    s_gen_isom_cox.insert(eP);
  };
  auto func_insert_pair_vertices=[&](FundDomainVertex<T,Tint> const& theVert, StatusEntry const& entry, FundDomainVertex<T,Tint> evert1, FundDomainVertex<T,Tint> evert2) -> void {
    //    theVert.l_roots.clear();
    //    std::cerr << "1 : func_insert_pair_vertices |theVert.l_roots|=" << theVert.MatRoot.rows() << "\n";
    //    std::cerr << "Before computing v_pair_char\n";
    pair_char<T> v_pair_char = gen_pair_char(G, evert1, evert2);
    //    std::cerr << "2 : func_insert_pair_vertices |theVert.l_roots|=" << theVert.MatRoot.rows() << "\n";

    size_t len = l_orbit_pair_vertices.size() / 2;
    for (size_t i=0; i<len; i++) {
      //      std::cerr <<  "Before LinPolytopeIntegralWMat_Isomorphism\n";
      //      std::cerr << "Before computing u_pair_char\n";
      const FundDomainVertex<T,Tint>& hvert1 = l_orbit_pair_vertices[2*i];
      const FundDomainVertex<T,Tint>& hvert2 = l_orbit_pair_vertices[2*i+1];
      pair_char<T> u_pair_char = gen_pair_char(G, hvert1, hvert2);
      //      std::cerr << "3 : func_insert_pair_vertices |theVert.l_roots|=" << theVert.MatRoot.rows() << "\n";
      std::optional<MyMatrix<T>> equiv_opt = LinPolytopeIntegralWMat_Isomorphism<T,Tgroup,std::vector<T>,uint16_t>(u_pair_char, v_pair_char);
      //      std::cerr << "4 : func_insert_pair_vertices |theVert.l_roots|=" << theVert.MatRoot.rows() << "\n";
      //      std::cerr <<  "After  LinPolytopeIntegralWMat_Isomorphism\n";
      if (equiv_opt) {
        std::cerr << "Find some isomorphism\n";
        /*
        std::cerr << "u : vert1=" << StringVectorGAP(hvert1.gen) << " vert2=" << StringVectorGAP(hvert2.gen) << "\n";
        std::cerr << "v : vert1=" << StringVectorGAP(evert1.gen) << " vert2=" << StringVectorGAP(evert2.gen) << "\n";
        std::cerr << "u : EXT=\n";
        WriteMatrix(std::cerr, u_pair_char.first);
        std::cerr << "u : WMat=\n";
        PrintWeightedMatrix(std::cerr, u_pair_char.second);
        //
        std::cerr << "v : EXT=\n";
        WriteMatrix(std::cerr, v_pair_char.first);
        std::cerr << "u : WMat=\n";
        PrintWeightedMatrix(std::cerr, v_pair_char.second);
        */
        f_insert_gen(UniversalMatrixConversion<Tint,T>(*equiv_opt));
        return;
      }
    }
    //    std::cerr << "ipass=" << ipass << "\n";
    //    v_pair.vert1.l_roots.clear();
    std::cerr << "Failed to find some isomorphism\n";
    //    std::cerr << "5 : func_insert_pair_vertices |theVert.l_roots|=" << theVert.MatRoot.rows() << "\n";
    l_entry.push_back(entry);
    //    std::cerr << "Before the automorphism insertions\n";
    for (auto & eGen : LinPolytopeIntegralWMat_Automorphism<T,Tgroup,std::vector<T>,uint16_t>(v_pair_char))
      f_insert_gen(UniversalMatrixConversion<Tint,T>(eGen));
    l_orbit_pair_vertices.push_back(evert1);
    l_orbit_pair_vertices.push_back(evert2);
  };
  size_t iVERT = 0;
  // We have to do a copy of the Vert since when the vector is extended the previous entries arfe desttroyed when a new
  // array is built. This would then invalidates a const& theVert reference.
  // See for details https://stackoverflow.com/questions/6438086/iterator-invalidation-rules-for-c-containers
  // Which writes: "vector: all iterators and references before the point of insertion are unaffected, unless
  // the new container size is greater than the previous capacity (in which case all iterators and references are
  // invalidated) [23.2.4.3/1]
  // Took 1 week to fully debug that problem.
  auto insert_edges_from_vertex=[&](FundDomainVertex<T,Tint> theVert) -> void {
    std::cerr << "insert_edges_from_vertex theVert=" << StringVectorGAP(RemoveFractionVector(theVert.gen)) << "\n";
    size_t n_root = theVert.MatRoot.rows();
    /*
    MyMatrix<T> CoxMat = ComputeCoxeterMatrix(G, theVert.l_roots).first;
    std::cerr << "CoxMat=\n";
    WriteMatrix(std::cerr, CoxMat);
    std::string symb = coxdyn_matrix_to_string(CoxMat);
    std::cerr << "Coxeter diagram of the vertex k=" << symb << "\n";
    */
    MyMatrix<T> FAC = UniversalMatrixConversion<T,Tint>(theVert.MatRoot);
    std::cerr << "FAC=\n";
    WriteMatrix(std::cerr, FAC);
    MyMatrix<T> FACred = ColumnReduction(FAC);
    std::cerr << "FACred=\n";
    WriteMatrix(std::cerr, FACred);
    std::cerr << "RankMat(FAC)=" << RankMat(FAC) << " RankMat(FACred)=" << RankMat(FACred) << "\n";
    vectface vf = lrs::DualDescription_temp_incd(FACred);
    std::cerr << "|vf|=" << vf.size() << "\n";
    for (auto & eIncd : vf) {
      TestFacetness(FACred, eIncd);
      MyVector<T> eFAC = RemoveFractionVector(FindFacetInequality(FACred, eIncd));
      std::cerr << "vf eIncd=" << StringFace(eIncd) << " eFAC=" << eFAC << "\n";
    }
    // Doing the std::sort(vf.begin(), vf.end()) seems science fiction right now.
    std::vector<Face> vf_ugly;
    for (auto & eFAC : vf)
      vf_ugly.push_back(eFAC);
    std::sort(vf_ugly.begin(), vf_ugly.end(),
              [&](Face const& x, Face const& y) -> bool {
                for (size_t i_root=0; i_root<n_root; i_root++) {
                  if (x[i_root] == 1 && y[i_root] == 0)
                    return true;
                  if (x[i_root] == 0 && y[i_root] == 1)
                    return false;
                }
                return false;
              });
    size_t pos=0;
    for (auto & eFAC : vf_ugly) {
      std::cerr << "pos=" << pos << " eFAC=" << StringFace(eFAC) << "\n";
      pos++;
    }

    size_t iFAC = 0;
    for (auto & eFAC : vf_ugly) {
      Face fFAC = eFAC;
      size_t i_disc = std::numeric_limits<size_t>::max();
      std::cerr << "\n";
      std::cerr << "iVERT=" << iVERT << " iFAC=" << iFAC << " n_root=" << n_root << " |eFAC|=" << eFAC.count() << "\n";
      std::vector<MyVector<Tint>> l_ui;
      for (size_t i_root=0; i_root<n_root; i_root++) {
        if (fFAC[i_root] == 1) {
          MyVector<Tint> root = GetMatrixRow(theVert.MatRoot, i_root);
          l_ui.push_back(root);
        } else {
          i_disc = i_root;
        }
      }
      MyVector<Tint> v_disc = GetMatrixRow(theVert.MatRoot, i_disc);
      std::cerr << "iVERT=" << iVERT << " iFAC=" << iFAC << " n_root=" << n_root << " |eFAC|=" << eFAC.count() << " i_disc=" << i_disc << "\n";
      FundDomainVertex<T,Tint> fVert = EdgewalkProcedure(G, theVert.gen, l_ui, l_norms, v_disc);
      T norm = fVert.gen.dot(G * fVert.gen);
      std::cerr << "k=" << StringVectorGAP(theVert.gen) << " l_ui=";
      for (auto & root : l_ui)
        std::cerr << " " << StringVectorGAP(root);
      std::cerr << "\n";
      std::cerr << "Result of EdgewalkProcedure\n";
      std::cerr << "We have fVert=" << StringVectorGAP(fVert.gen) << " norm=" << norm << "\n";
      std::cerr << "MatRoot=\n";
      WriteMatrix(std::cerr, fVert.MatRoot);
      //
      StatusEntry entry{false,true};
      func_insert_pair_vertices(theVert, entry, theVert, fVert);
      iFAC++;
    }
    iVERT++;
    std::cerr << "Exiting from the insert_edges_from_vertex\n";
    //    throw TerminalException{1};
  };
  insert_edges_from_vertex(eVert);
  while(true) {
    bool IsFinished = true;
    size_t len = l_entry.size();
    for (size_t i=0; i<len; i++) {
      if (l_entry[i].stat1) {
        l_entry[i].stat1 = false;
        insert_edges_from_vertex(l_orbit_pair_vertices[2*i]);
        IsFinished = false;
      }
      if (l_entry[i].stat2) {
        l_entry[i].stat2 = false;
        insert_edges_from_vertex(l_orbit_pair_vertices[2*i+1]);
        IsFinished = false;
      }
    }
    if (IsFinished)
      break;
  }
  std::cerr << "Exiting from the infinite loop of enumeration of vertex pairs\n";
  std::vector<MyMatrix<Tint>> l_gen_isom_cox;
  for (auto & e_gen : s_gen_isom_cox)
    l_gen_isom_cox.push_back(e_gen);
  return {l_gen_isom_cox, l_orbit_pair_vertices};
}





template<typename T, typename Tint>
FundDomainVertex<T,Tint> get_initial_vertex(MyMatrix<T> const& G, std::vector<T> const& l_norms, std::string const& OptionInitialVertex, std::string const& FileInitialVertex)
{
  if (OptionInitialVertex == "File") {
    if (!IsExistingFile(FileInitialVertex)) {
      std::cerr << "The file FileInitialVertex=" << FileInitialVertex << " is missing\n";
      throw TerminalException{1};
    }
    std::ifstream is(FileInitialVertex);
    MyVector<T> gen = ReadVector<T>(is);
    MyMatrix<Tint> Mroot = ReadMatrix<Tint>(is);
    return {gen, Mroot};
  }
#ifdef ALLOW_VINBERG_ALGORITHM_FOR_INITIAL_VERTEX
  if (OptionInitialVertex == "vinberg") {
    VinbergTot<T,Tint> Vtot = GetVinbergFromG<T,Tint>(G, l_norms);
    std::pair<MyVector<Tint>, std::vector<MyVector<Tint>>> epair = FindOneInitialRay(Vtot);
    return {UniversalVectorConversion<T,Tint>(epair.first), MatrixFromVectorFamily(epair.second)};
  }
#endif
  std::cerr << "Failed to find a matching entry in get_initial_vertex\n";
  std::cerr << "OptionInitialVertex=" << OptionInitialVertex << " but allowed values are File\n";
#ifdef ALLOW_VINBERG_ALGORITHM_FOR_INITIAL_VERTEX
  std::cerr << "and vinberg has also been allowed\n";
#else
  std::cerr << "option vinberg has not been allowed\n";
#endif
  throw TerminalException{1};
}



template<typename T, typename Tint, typename Tgroup>
void MainFunctionEdgewalk(FullNamelist const& eFull)
{
  SingleBlock BlockPROC=eFull.ListBlock.at("PROC");
  std::string FileLorMat=BlockPROC.ListStringValues.at("FileLorMat");
  MyMatrix<T> G = ReadMatrixFile<T>(FileLorMat);
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(G);
  if (DiagInfo.nbZero != 0 || DiagInfo.nbMinus != 1) {
    std::cerr << "G=\n";
    WriteMatrix(std::cerr, G);
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
  FundDomainVertex<T,Tint> eVert = get_initial_vertex<T,Tint>(G, l_norms, OptionInitialVertex, FileInitialVertex);
  std::cerr << "Initial vertex is\n";
  std::cerr << "eVert.gen=" << StringVectorGAP(eVert.gen) << "\n";
  std::cerr << "l_roots=\n";
  WriteMatrix(std::cerr, eVert.MatRoot);
  //
  ResultEdgewalk<T,Tint> re = LORENTZ_RunEdgewalkAlgorithm<T,Tint,Tgroup>(G, l_norms, eVert);
  std::string OutFormat=BlockPROC.ListStringValues.at("OutFormat");
  std::string FileOut=BlockPROC.ListStringValues.at("FileOut");
  bool ComputeAllSimpleRoots=BlockPROC.ListBoolValues.at("ComputeAllSimpleRoots");
  if (FileOut == "stderr") {
    PrintResultEdgewalk(G, re, std::cerr, OutFormat, ComputeAllSimpleRoots);
  } else {
    if (FileOut == "stdout") {
      PrintResultEdgewalk(G, re, std::cout, OutFormat, ComputeAllSimpleRoots);
    } else {
      std::ofstream os(FileOut);
      PrintResultEdgewalk(G, re, os, OutFormat, ComputeAllSimpleRoots);
    }
  }

}



#endif

#ifndef INCLUDE_EDGEWALK_H
#define INCLUDE_EDGEWALK_H

#include "MatrixCanonicalForm.h"
#include "Temp_PolytopeEquiStab.h"
#include "two_dim_lorentzian.h"
#include "coxeter_dynkin.h"
#include "vinberg.h"
#include "lorentzian_linalg.h"
#include "Namelist.h"
#include "Temp_Positivity.h"
#include "POLY_lrslib.h"
#include "POLY_cddlib.h"


#define ALLOW_VINBERG_ALGORITHM_FOR_INITIAL_VERTEX



FullNamelist NAMELIST_GetStandard_EDGEWALK()
{
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListStringValues1["FileLorMat"] = "the lorentzian matrix used";
  ListStringValues1["OptionInitialVertex"] = "vinberg or FileVertex or FileVertexRoots and if FileVertex or FileVertexRoots selected use FileVertDomain as initial vertex";
  ListStringValues1["FileInitialVertex"] = "unset put the name of the file used for the initial vertex";
  ListStringValues1["OptionNorms"] = "possible option K3 (then just 2) or all where all norms are considered";
  ListStringValues1["OutFormat"] = "GAP for gap use or TXT for text output";
  ListStringValues1["FileOut"] = "stdout, or stderr or the filename of the file you want to write to";
  ListBoolValues1["ApplyReduction"]=true; // Normally, we want to ApplyReduction, this is for debug only
  ListBoolValues1["ComputeAllSimpleRoots"]=true;
  SingleBlock BlockPROC;
  BlockPROC.ListStringValues = ListStringValues1;
  BlockPROC.ListBoolValues = ListBoolValues1;
  ListBlock["PROC"]=BlockPROC;
  // Merging all data
  return {ListBlock, "undefined"};
}


FullNamelist NAMELIST_GetStandard_EDGEWALK_Isomorphism()
{
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListStringValues1["FileLorMat1"] = "the lorentzian matrix used";
  ListStringValues1["FileLorMat2"] = "the lorentzian matrix used";
  ListStringValues1["OptionNorms"] = "possible option K3 (then just 2) or all where all norms are considered";
  ListStringValues1["OutFormat"] = "GAP for gap use or TXT for text output";
  ListStringValues1["FileOut"] = "stdout, or stderr or the filename of the file you want to write to";
  ListBoolValues1["ApplyReduction"]=true; // Normally, we want to ApplyReduction, this is for debug only
  SingleBlock BlockPROC;
  BlockPROC.ListStringValues = ListStringValues1;
  BlockPROC.ListBoolValues = ListBoolValues1;
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
  MyMatrix<Tint> MatRoot;
};


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
  std::cerr << "Failed to find a matching entry for WriteFundDomainVertex\n";
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
  std::vector<RootCandidateCuspidal<T,Tint>> l_candidates;
  for (auto & e_extension : l_extension) {
    //    std::cerr << "res_norm=" << e_extension.res_norm << "\n";
    if (e_extension.res_norm == 0) {
      MyMatrix<T> Latt = ComputeLattice_LN(G, e_extension.e_norm);
      std::optional<MyVector<T>> opt_v = ResolveLattEquation(Latt, e_extension.u_component, k);
      //      std::cerr << "We have opt_v\n";
      if (opt_v) {
        const MyVector<T>& v_T = *opt_v;
        //        std::cerr << "Proposed v_T=" << StringVectorGAP(v_T) << "\n";
        RootCandidateCuspidal<T,Tint> e_cand = gen_possible_cuspidalextension<T,Tint>(G, kP, v_T, e_extension.e_norm);
        l_candidates.push_back(e_cand);
      }
    }
  }
  std::cerr << "DetermineRootsCuspidalCase : |l_candidates|=" << l_candidates.size() << "\n";
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
    std::cerr << "x : sign=" << x.sign << " quant=" << x.quant << " norm=" << x.e_norm << " v=" << StringVectorGAP(x.v) << "\n";
  }
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
  std::cerr << "DetermineRootsCuspidalCase, exiting |l_ui_ret|=" << l_ui_ret.size() << "\n";
  return l_ui_ret;
}





template<typename Tint>
struct AdjacencyDirection {
  std::vector<MyVector<Tint>> l_ui;
  MyVector<Tint> v_disc;
};


template<typename Tint>
void PrintAdjacencyDirection(std::ostream & os, AdjacencyDirection<Tint> const& ad)
{
  os << "l_ui =";
  for (auto & root : ad.l_ui)
    os << " " << StringVectorGAP(root);
  os << "\n";
  os << "v_disc =" << StringVectorGAP(ad.v_disc) << "\n";
}



template<typename Tint>
AdjacencyDirection<Tint> GetAdjacencyDirection(MyMatrix<Tint> const& MatRoot, Face const& f)
{
  size_t n_root = MatRoot.rows();
  size_t i_disc = std::numeric_limits<size_t>::max();
  std::vector<MyVector<Tint>> l_ui;
  for (size_t i_root=0; i_root<n_root; i_root++) {
    if (f[i_root] == 1) {
      MyVector<Tint> root = GetMatrixRow(MatRoot, i_root);
      l_ui.push_back(root);
    } else {
      i_disc = i_root;
    }
  }
  MyVector<Tint> v_disc = GetMatrixRow(MatRoot, i_disc);
  return {l_ui, v_disc};
}




/*
  We take the notations as in EDGEWALK paper.
  ---The dimension is n+1
  ---We have an edge e between two rays k and k'.
  We have u1, .... u(n-1) roots that are orthogonal and U the real span of those vectors.
  ---P is the real plane orthogonal to U.
  ---pi_U and pi_P are the corresponding projectors
  ---How is (1/2) P defined and correspond to k (typo correction)
 */
template<typename T, typename Tint>
FundDomainVertex<T,Tint> EdgewalkProcedure(MyMatrix<T> const& G, std::vector<T> const& l_norms, MyVector<T> const& k, AdjacencyDirection<Tint> const ad)
{
  const std::vector<MyVector<Tint>>& l_ui = ad.l_ui;
  const MyVector<Tint>& v_disc = ad.v_disc;
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
  MyMatrix<T> EquaRvect(n_root+1,n);
  for (size_t i_root=0; i_root<n_root; i_root++) {
    MyVector<T> eV = UniversalVectorConversion<T,Tint>(l_ui[i_root]);
    MyVector<T> eP = G * eV;
    T eScal = k.dot(eP);
    if (eScal != 0) {
      std::cerr << "The scalar product should be 0\n";
      throw TerminalException{1};
    }
    AssignMatrixRow(EquaRvect, i_root, eP);
  }
  MyMatrix<T> Pplane = Get_Pplane(G, l_ui);
  //  std::cerr << "Plane P=[" << StringVectorGAP(GetMatrixRow(Pplane,0)) << ", " << StringVectorGAP(GetMatrixRow(Pplane,1)) << "]\n";
  MyVector<T> eP = G * k;
  T norm = k.dot(eP);
  std::cerr << "k=" << StringVectorGAP(k) << " norm=" << norm << "\n";
  AssignMatrixRow(EquaRvect, n_root, eP);
  MyMatrix<T> NSP = NullspaceTrMat(EquaRvect);
  if (NSP.rows() != 1) {
    std::cerr << "|NSP|=" << NSP.rows() << "/" << NSP.cols() << "\n";
    std::cerr << "The dimension should be exactly 1\n";
    throw TerminalException{1};
  }
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
  std::vector<RootCandidate<T,Tint>> l_candidates;
  bool only_spherical = true;
  std::vector<Possible_Extension<T>> l_extension = ComputePossibleExtensions(G, l_ui, l_norms, only_spherical);
  std::cerr << "EdgewalkProcedure : |l_extension|=" << l_extension.size() << "\n";

  //  for (auto & kv : map_max_resnorm)
  //    std::cerr << "kv : norm=" << kv.first << " max(res_norm)=" << kv.second << "\n";
  //
  // Determine if the plane P is isotropic and if not compute the set of test vectors
  //
  MyMatrix<T> G_Pplane = Pplane * G * Pplane.transpose();
  //  std::cerr << "G_Pplane=\n";
  //  WriteMatrix(std::cerr, G_Pplane);
  //  std::cerr << "G=\n";
  //  WriteMatrix(std::cerr, G);
  struct SingCompAnisotropic {
    MyMatrix<T> Latt;
    MyVector<Tint> r0_work;
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
    LatticeProjectionFramework<T> ProjFram(G, Pplane, Latt);
    MyMatrix<T> BasisProj = ProjFram.BasisProj;
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
    //    MyVector<T> Basis0 = GetMatrixRow(BasisProj, 0);
    //    std::cerr << "Basis0=" << StringVectorGAP(Basis0) << "\n";
    //    MyVector<T> Basis1 = GetMatrixRow(BasisProj, 1);
    //    std::cerr << "Basis1=" << StringVectorGAP(Basis1) << "\n";
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
    std::cerr << " -------------- get_sing_comp_anisotropic, e_norm=" << e_norm << " ------------------------\n";
    MyMatrix<T> Latt = ComputeLattice_LN(G, e_norm);
    //    std::cerr << "Latt=" << StringMatrixGAP(Latt) << "\n";
    MyMatrix<T> Basis_ProjP_LN = get_basis_projp_ln(Latt);
    //    std::cerr << "Basis_ProjP_LN=\n";
    //    WriteMatrix(std::cerr, Basis_ProjP_LN);
    MyMatrix<T> Basis_P_inter_LN = IntersectionLattice_VectorSpace(Latt, Pplane);
    MyMatrix<T> Gwork = Basis_ProjP_LN * G * Basis_ProjP_LN.transpose();
    //    std::cerr << "Gwork=" << StringMatrixGAP(Gwork) << "\n";
    // The residual norm is res_norm = N - u_{N,Delta}^2
    // u_{N,Delta} belongs to a positive definite lattice.
    // Therefore res_norm <= N
    // Thus we compute all the vectors up to norm res_norm because
    // res_norm is always realizable with the vector 2 2 2 ..... 2
    T res_norm = e_norm;
    MyVector<Tint> r0_work = get_r0work(Basis_ProjP_LN, r0);
    T r0_norm = eval_quad(Gwork, r0_work);
    //    std::cerr << "r0_norm=" << r0_norm << " e_norm=" << e_norm << "\n";
    MyVector<Tint> l_A = GetTwoComplement(r0_work);
    //    std::cerr << "l_A=" << StringVectorGAP(l_A) << " res_norm=" << res_norm << "\n";
    MyVector<Tint> l_B = Canonical(Gwork, r0_norm, r0_work, l_A);
    //    std::cerr << "l_B=" << StringVectorGAP(l_B) << "\n";
    //    std::cerr << "get_sing_comp_anisotropic, step 2\n";
    std::optional<std::pair<MyMatrix<Tint>,std::vector<MyVector<Tint>>>> opt = Anisotropic<T,Tint>(Gwork, res_norm, r0_work, l_B);
    //    std::cerr << "get_sing_comp_anisotropic, step 3\n";
    if (!opt) { // No solution, this definitely can happen
      return {Latt, r0_work, Basis_ProjP_LN, Basis_P_inter_LN, Gwork, {}};
    }
    const std::vector<MyVector<Tint>>& l_vect1 = opt->second;
    std::cerr << "|l_vect1|=" << l_vect1.size() << "\n";
    const MyMatrix<Tint>& TransformSma = opt->first;
    MyMatrix<T> TransformSma_T = UniversalMatrixConversion<T,Tint>(TransformSma);
    /*
      Transformation rule
     */
    MyMatrix<T> TransformBig_red = IdentityMat<T>(n);
    for (int i=0; i<2; i++)
      for (int j=0; j<2; j++)
        TransformBig_red(i,j) = UniversalScalarConversion<T,Tint>(TransformSma(i,j));
    MyMatrix<T> l_ui_MatT = UniversalMatrixConversion<T,Tint>(MatrixFromVectorFamily(l_ui));
    MyMatrix<T> ReductionBasis = Concatenate(Basis_ProjP_LN, l_ui_MatT);
    MyMatrix<T> TransformBig = Inverse(ReductionBasis) * TransformBig_red * ReductionBasis;
    //    std::cerr << "TransformBig=" << StringMatrixGAP(TransformBig) << "\n";

    //    std::cerr << "TransformSma=" << StringMatrixGAP(TransformSma) << "\n";
    //    std::cerr << "Basis_ProjP_LN=\n";
    //    WriteMatrix(std::cerr, Basis_ProjP_LN);
    //    std::cerr << "Basis_P_inter_LN=\n";
    //    WriteMatrix(std::cerr, Basis_P_inter_LN);
    MyMatrix<T> Expr_t = ExpressVectorsInIndependentFamilt(Basis_P_inter_LN, Basis_ProjP_LN);
    //    std::cerr << "Expr_t=" << StringMatrixGAP(Expr_t) << "\n";
    //    WriteMatrix(std::cerr, Expr_t);
    //    std::cerr << "get_sing_comp_anisotropic, step 4\n";
    if (!IsIntegralMatrix(Expr_t)) {
      std::cerr << "The matrix should be integral\n";
      throw TerminalException{1};
    }
    //    MyMatrix<Tint> Expr_i = UniversalMatrixConversion<Tint,T>(Expr_t);
    size_t order = GetMatrixExponentSublattice_TrivClass(TransformSma_T, Expr_t);
    std::cerr << "order=" << order << "\n";
    std::vector<MyMatrix<Tint>> l_vect2;
    //    std::cerr << "get_sing_comp_anisotropic, step 5\n";
    for (auto & e_vect1 : l_vect1) {
      T norm1 = eval_quad(Gwork, e_vect1);
      size_t ord = 1;
      while(true) {
        T norm2 = ord * ord * norm1;
        if (norm2 > res_norm)
          break;
        MyVector<Tint> e_vect2 = ord * e_vect1;
        //        std::cerr << "norm2=" << norm2 << " e_vect2=" << StringVectorGAP(e_vect2) << "\n";
        l_vect2.push_back(e_vect2);
        ord++;
      }
    }
    std::cerr << "|l_vect2|=" << l_vect2.size() << "\n";
    std::vector<MyVector<Tint>> l_vect3;
    MyMatrix<Tint> TheMat = IdentityMat<Tint>(2);
    for (size_t i=0; i<order; i++) {
      for (auto & e_vect2 : l_vect2) {
        MyVector<Tint> e_vect3 = TheMat.transpose() * e_vect2;
        l_vect3.push_back(e_vect3);
      }
      TheMat = TheMat * TransformSma;
    }
    std::cerr << "|l_vect3|=" << l_vect3.size() << "\n";
    return {Latt, r0_work, Basis_ProjP_LN, Basis_P_inter_LN, Gwork, l_vect3};
  };
  auto get_sing_comp_isotropic=[&](T const& e_norm) -> SingCompIsotropic {
    std::cerr << "get_sing_comp_isotropic, e_norm=" << e_norm << "\n";
    MyMatrix<T> Latt = ComputeLattice_LN(G, e_norm);
    //    std::cerr << "Latt=" << StringMatrixGAP(Latt) << "\n";
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
      std::cerr << "Case is_isotropic = true\n";
      for (auto & u_norm : l_norms)
        map_isotropic[u_norm] = get_sing_comp_isotropic(u_norm);
    } else {
      is_isotropic = false;
      std::cerr << "Case is_isotropic = false\n";
      for (auto & u_norm : l_norms)
        map_anisotropic[u_norm] = get_sing_comp_anisotropic(u_norm);
    }
  }
  std::cerr << "Edgewalk Procedure, step 7\n";
  // Evaluation of fun
  auto get_next_anisotropic=[&](Possible_Extension<T> const& poss) -> std::optional<MyVector<Tint>> {
    //    std::cerr << "Beginning of get_next_anisotropic\n";
    T const& e_norm = poss.e_norm;
    SingCompAnisotropic const& e_comp = map_anisotropic[e_norm];
    for (auto & e_vect : e_comp.l_vect) {
      T val = eval_quad(e_comp.Gwork, e_vect);
      //      std::cerr << "gna, step 2 e_vect=" << StringVectorGAP(e_vect) << " val=" << val << "\n";
      if (val == poss.res_norm) {
        //        std::cerr << "e_vect=" << StringVectorGAP(e_vect) << " Basis_ProjP_LN=[" << StringVectorGAP(GetMatrixRow(e_comp.Basis_ProjP_LN, 0)) << ", " << StringVectorGAP(GetMatrixRow(e_comp.Basis_ProjP_LN, 1)) << "]\n";
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
      //      std::cerr << "Testing adequateness of x=" << x << " norm=" << norm << "\n";
      // In the isotropic case, there is no relevant test
      // See Right picture on Figure 8.1 of the edgewalk paper
      if (norm < 0) {
        T scal = eval_scal(e_comp.GP_LN, e_comp.r0_work, x);
        if (scal <= 0) {
          //          std::cerr << "scal=" << scal << "\n";
          return false;
        }
      }
      Tint eDet = det(e_comp.r0_work, x);
      //      std::cerr << "eDet=" << eDet << "\n";
      return eDet > 0;
    };
    std::cerr << "|l_vect1|=" << l_vect1.size() << "\n";
    for (auto & e_vect1 : l_vect1)
      if (is_corr(e_vect1))
        l_vect2.push_back(e_vect1);
    std::cerr << "We found |l_vect2|=" << l_vect2.size() << "\n";
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
    T const& res_norm = poss.res_norm;
    std::cerr << "get_next_isotropic with e_norm=" << e_norm << " res_norm=" << res_norm << "\n";
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
    T e_norm = e_extension.e_norm;
    T res_norm = e_extension.res_norm;
    std::cerr << "------ u_component=" << StringVectorGAP(e_extension.u_component) << " norm=" << e_norm << " res_norm=" << res_norm << " ----------\n";


    std::optional<MyVector<Tint>> opt_v = get_next(e_extension);
    //    std::cerr << "We have opt_v\n";
    if (opt_v) {
      MyVector<Tint> const& alpha = *opt_v;
      std::cerr << "alpha=" << StringVectorGAP(alpha) << "\n";
      std::vector<MyVector<Tint>> l_roots = l_ui;
      l_roots.push_back(alpha);
      MyMatrix<T> Mat_root = UniversalMatrixConversion<T,Tint>(MatrixFromVectorFamily(l_roots));
      MyMatrix<T> EquaMat = Mat_root * G;
      MyMatrix<T> NSP = NullspaceTrMat(EquaMat);
      if (NSP.rows() != 1) {
        std::cerr << "We should have exactly one row\n";
        throw TerminalException{1};
      }
      MyVector<T> gen = GetMatrixRow(NSP, 0);
      std::cerr << "gen=" << StringVectorGAP(gen) << " k=" << StringVectorGAP(k) << "\n";
      T scal = gen.dot(G * k);
      //        std::cerr << "We have scal\n";
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
          //            std::cerr << "k_new=" << StringVectorGAP(k_new) << "\n";
          RootCandidate<T,Tint> eCand = gen_possible_extension(G, k, alpha, res_norm, e_norm, fund_v);
          l_candidates.push_back(eCand);
        }
      }
    }
  }
  std::cerr << "EdgewalkProcedure : |l_candidates|=" << l_candidates.size() << "\n";
  if (l_candidates.size() > 0) {
    /*
    for (auto e_cand : l_candidates)
      std::cerr << "e_cand sign=" << e_cand.sign << " quant1=" << e_cand.quant1 << " quant2=" << e_cand.quant2 << " e_norm=" << e_cand.e_norm << " fund_v=" << StringVectorGAP(e_cand.fund_v.gen) << " alpha=" << StringVectorGAP(e_cand.alpha) << "\n";
    */
    RootCandidate<T,Tint> best_cand = get_best_candidate(l_candidates);
    //    std::cerr << "fund_v=" << StringVectorGAP(best_cand.fund_v.gen) << "\n";
    //    std::cerr << "MatRoot=\n";
    //    WriteMatrix(std::cerr, best_cand.fund_v.MatRoot);
    return best_cand.fund_v;
  }
  std::cerr << "         --------------- Looking for an isotropic vector ------------\n";
  // So, no candidates were found. We need to find isotropic vectors.
  const MyMatrix<T> Gred = Pplane * G * Pplane.transpose();
  //  std::cerr << "We have Gred=\n";
  //  WriteMatrix(std::cerr, Gred);
  std::optional<MyMatrix<T>> Factor_opt = GetIsotropicFactorization(Gred);
  //  std::cerr << "We have Factor_opt\n";
  if (!Factor_opt) {
    std::cerr << "The matrix is not isotropic. Major rethink are needed\n";
    throw TerminalException{1};
  }
  MyMatrix<T> const& Factor = *Factor_opt;
  //  std::cerr << "We have Factor=\n";
  //  WriteMatrix(std::cerr, Factor);
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
    //    std::cerr << "i=" << i << "\n";
    // a x + b y correspond to the ray (u0, u1) = (-b, a)
    MyVector<T> U(2);
    U(0) = -Factor(i,1);
    U(1) =  Factor(i,0);
    //    std::cerr << "U=" << StringVectorGAP(U) << "\n";
    T sum = U.dot(Gred * U);
    //    std::cerr << "sum=" << sum << "\n";
    MyVector<T> gen = Pplane.transpose() * U;
    //    std::cerr << "k="; WriteVectorGAP(std::cerr, k); std::cerr << "\n";
    //    std::cerr << "r0="; WriteVectorGAP(std::cerr, r0); std::cerr << "\n";
    //    std::cerr << "gen="; WriteVectorGAP(std::cerr, gen); std::cerr << "\n";
    T sum_B = gen.dot(G * gen);
    //    std::cerr << "sum_B=" << sum_B << "\n";
    //    std::cerr << "gen=" << StringVectorGAP(gen) << "\n";
    if (!IsVectorMultiple(gen, k)) {
      MyVector<T> can_gen = get_can_gen(gen);
      //      std::cerr << "can_gen=" << StringVectorGAP(can_gen) << "\n";
      //      std::cerr << "RemoveFraction(can_gen)=" << StringVectorGAP(RemoveFractionVector(can_gen)) << "\n";
      T scal = v_disc_t.dot(G * can_gen);
      //      std::cerr << "scal=" << scal << "\n";
      if (scal < 0) // Convention is negative scalar in Lorentzian theory (see end of sect 2 of edgewalk paper)
        l_gens.push_back(can_gen);
    }
  }
  //  std::cerr << "|l_gens|=" << l_gens.size() << "\n";
  if (l_gens.size() != 1) {
    std::cerr << "We should have just one vector in order to conclude. Rethink needed\n";
    throw TerminalException{1};
  }
  const MyVector<T> & k_new = l_gens[0];
  //  std::cerr << "k_new=" << StringVectorGAP(RemoveFractionVector(k_new)) << "\n";
  std::vector<MyVector<Tint>> l_roots_ret = DetermineRootsCuspidalCase(G, l_ui, l_norms, k_new, k);
  return {RemoveFractionVector(k_new), MatrixFromVectorFamily(l_roots_ret)};
}




template<typename T>
using pair_char = std::pair<MyMatrix<T>,WeightMatrix<true,std::vector<T>,uint16_t>>;




template<typename T, typename Tint, typename Tgroup>
struct FundDomainVertex_FullInfo {
  FundDomainVertex<T,Tint> vert;
  pair_char<T> e_pair_char;
  size_t hash;
  Tgroup GRP1;
};


template<typename T, typename Tint, typename Tgroup>
FundDomainVertex_FullInfo<T,Tint,Tgroup> DirectCopy(FundDomainVertex_FullInfo<T,Tint,Tgroup> const& fdfi)
{
  return {fdfi.vert, {fdfi.e_pair_char.first, fdfi.e_pair_char.second.DirectCopy()}, fdfi.hash, fdfi.GRP1};
}






template<typename T, typename Tint, typename Tgroup>
FundDomainVertex_FullInfo<T,Tint,Tgroup> gen_fund_domain_fund_info(MyMatrix<T> const& G, std::vector<T> const& l_norms, FundDomainVertex<T,Tint> const& vert)
{
  std::cerr << "gen_fund_domain_fund_info, beginning\n";
  //
  // Put the stuff that can help for invariant first
  std::unordered_map<MyVector<Tint>,int> map_v;
  size_t len = vert.MatRoot.rows();
  for (size_t i=0; i<len; i++) {
    MyVector<Tint> eV = GetMatrixRow(vert.MatRoot, i);
    map_v[eV] = 1;
  }
  MyVector<Tint> gen_tint = UniversalVectorConversion<Tint,T>(RemoveFractionVector(vert.gen));
  map_v[gen_tint] = 2;
  std::cerr << "Initial map_v built\n";
  //
  using Tidx = uint32_t;
  using Tidx_value = uint16_t;
  using Telt = typename Tgroup::Telt;
  using Telt_idx = typename Telt::Tidx;
  using Tgr = GraphListAdj;
  std::vector<MyMatrix<T>> ListMat{G};
  struct ret_type {
    pair_char<T> e_pair_char;
    Tgroup GRP1;
    MyMatrix<Tint> MatRoot;
  };
  auto get_canonicalized_record=[&](std::unordered_map<MyVector<Tint>,int> const& the_map) -> ret_type {
    size_t n_row = the_map.size();
    std::vector<MyVector<Tint>> l_vect;
    std::vector<T> Vdiag;
    size_t idx = 0;
    std::vector<size_t> V;
    for (auto & kv : the_map) {
      l_vect.push_back(kv.first);
      Vdiag.push_back(T(kv.second));
      if (kv.second == 1)
        V.push_back(idx);
      idx++;
    }
    size_t n1 = V.size();
    MyMatrix<T> MatV = UniversalMatrixConversion<T,Tint>(MatrixFromVectorFamily(l_vect));
    WeightMatrix<true, std::vector<T>, Tidx_value> WMat = GetWeightMatrix_ListMat_Vdiag<T,Tidx,Tidx_value>(MatV, ListMat, Vdiag);
    WMat.ReorderingSetWeight();
    std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> epair = GetGroupCanonicalizationVector_Kernel<std::vector<T>,Tgr,Tidx,Tidx_value>(WMat);
    const std::vector<Tidx>& ListIdx = epair.first;
    const std::vector<std::vector<Tidx>>& ListGen = epair.second;
    WMat.RowColumnReordering(ListIdx);
    std::vector<Tidx> ListIdxRev(n_row);
    for (size_t i1=0; i1<n_row; i1++)
      ListIdxRev[ListIdx[i1]] = i1;
    std::vector<MyVector<Tint>> l_vect_reord(n_row);
    std::vector<T> Vdiag_reord(n_row);
    std::vector<MyVector<Tint>> l_vect1;
    std::vector<size_t> Map1;
    std::vector<size_t> Map1_rev(n_row,std::numeric_limits<size_t>::max());
    size_t pos = 0;
    for (size_t i=0; i<n_row; i++) {
      size_t j = ListIdx[i];
      l_vect_reord[i] = l_vect[j];
      Vdiag_reord[i] = Vdiag[j];
      if (Vdiag[j] == 1) {
        l_vect1.push_back(l_vect[j]);
        Map1.push_back(i);
        Map1_rev[i] = pos;
        pos++;
      }
    }
    MyMatrix<T> MatV_reord = UniversalMatrixConversion<T,Tint>(MatrixFromVectorFamily(l_vect_reord));
    // There are two use case of computing the group
    // ---For the computation of minimal adjacencies that would get us a full dimensional system
    // ---For the computation of orbits of adjacent vertices
    // So, in both cases, we need to reduce to the group for values 1.
    std::vector<Telt> LGen;
    std::vector<Telt> LGen1;
    for (auto & eGen : ListGen) {
      std::vector<Telt_idx> V(n_row);
      for (size_t i1=0; i1<n_row; i1++) {
        Tidx i2 = ListIdx[i1];
        Tidx i3 = eGen[i2];
        Tidx i4 = ListIdxRev[i3];
        V[i1] = i4;
      }
      Telt ePerm(V);
      LGen.push_back(ePerm);
      //
      std::vector<Telt_idx> V1(n1);
      for (size_t i1=0; i1<n1; i1++) {
        size_t i2 = Map1[i1];
        size_t i3 = V[i2];
        size_t i4 = Map1_rev[i3];
        V1[i1] = i4;
      }
      Telt ePerm1(V1);
      LGen1.push_back(ePerm1);
    }
    Tgroup GRP(LGen, n_row);
    Tgroup GRP1(LGen1, n1);
    MyMatrix<Tint> MatV_ret = MatrixFromVectorFamily(l_vect1);
    pair_char<T> e_pair_char{std::move(MatV_reord), std::move(WMat)};
    return {std::move(e_pair_char), std::move(GRP1), std::move(MatV_ret)};
  };
  T norm = vert.gen.dot(G * vert.gen);
  std::cerr << "norm=" << norm << "\n";
  if (norm == 0) {
    // In isotropic case, we are unfortunately forced to do more complex stuff
    // Those needs
    ret_type erec = get_canonicalized_record(map_v);
    //    std::cerr << "We have erec\n";
    // Add new vertices to
    MyMatrix<T> FAC = UniversalMatrixConversion<T,Tint>(erec.MatRoot);
    MyMatrix<T> FACred = ColumnReduction(FAC);
    vectface vf = lrs::DualDescription_temp_incd(FACred);
    //    std::cerr << "We have vf\n";
    // Finding the minimal orbit and then
    vectface vf_min = OrbitSplittingSet_GetMinimalOrbit(vf, erec.GRP1);
    std::cerr << "|vf_min|=" << vf_min.size() << " gen=" << StringVectorGAP(vert.gen) << " |GRP1|=" << erec.GRP1.size() << "\n";
    for (auto & eFAC : vf_min) {
      AdjacencyDirection<Tint> ad = GetAdjacencyDirection(erec.MatRoot, eFAC);
      FundDomainVertex<T,Tint> fVert = EdgewalkProcedure(G, l_norms, vert.gen, ad);
      MyVector<Tint> fVert_tint = UniversalVectorConversion<Tint,T>(RemoveFractionVector(fVert.gen));
      map_v[fVert_tint] = 3;
    }
    std::cerr << "map_v has been extended\n";
  }
  ret_type frec = get_canonicalized_record(map_v);
  const auto& WMat = frec.e_pair_char.second;
  size_t seed = 1440;
  size_t hash = ComputeHashWeightMatrix_raw(WMat, seed);
  FundDomainVertex<T,Tint> new_vert{vert.gen, frec.MatRoot};
  std::cerr << "gen_fund_domain_fund_info gen=" << StringVectorGAP(vert.gen) << " |GRP1|=" << frec.GRP1.size() << "\n";
  return {std::move(new_vert), std::move(frec.e_pair_char), hash, std::move(frec.GRP1)};
}





template<typename T, typename Tint>
struct ResultEdgewalk {
  std::vector<MyMatrix<Tint>> l_gen_isom_cox;
  std::vector<FundDomainVertex<T,Tint>> l_orbit_vertices;
};


template<typename T, typename Tint>
std::vector<MyVector<Tint>> compute_full_root_orbit(ResultEdgewalk<T,Tint> const& re)
{
  std::unordered_set<MyVector<Tint>> TotalList;
  auto f_insert=[&](MyVector<Tint> const& v) -> void {
    //    std::cerr << "f_insert call with v=" << StringVectorGAP(v) << "\n";
    if (TotalList.count(v) != 0)
      return;
    std::unordered_set<MyVector<Tint>> s_v;
    std::vector<MyVector<Tint>> l_v;
    auto f_ins=[&](MyVector<Tint> const& w) -> void {
      if (s_v.count(w) != 0)
        return;
      s_v.insert(w);
      l_v.push_back(w);
      TotalList.insert(w);
    };
    f_ins(v);
    size_t pos=0;
    while(true) {
      size_t len = l_v.size();
      //      std::cerr << "pos=" << pos << " len=" << len << "\n";
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
  for (auto & fdv : re.l_orbit_vertices) {
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
  for (auto & evert : re.l_orbit_vertices)
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
    l_simple_root = compute_full_root_orbit(re);
  if (OutFormat == "GAP") {
    os << "return rec(l_norms:=";
    WriteStdVectorGAP(os, l_norms);
    os << ", ListIsomCox:=";
    WriteVectorMatrixGAP(os, re.l_gen_isom_cox);
    os << ", ListVertices:=[";
    bool IsFirst = true;
    size_t len = re.l_orbit_vertices.size();
    for (size_t i=0; i<len; i++) {
      if (!IsFirst)
        os << ",\n";
      IsFirst = false;
      const FundDomainVertex<T,Tint> & evert = re.l_orbit_vertices[ i ];
      WriteFundDomainVertex(G, evert, os, OutFormat);
    }
    os << "], n_orbit_vertices:=" << len;
    if (ComputeAllSimpleRoots) {
      os << ", ListSimpleRoots:=[";
      size_t n_simple = l_simple_root.size();
      for (size_t i=0; i<n_simple; i++) {
        if (i>0)
          os << ",";
        os << StringVectorGAP(l_simple_root[i]);
      }
      os << "], n_simple:=" << n_simple;
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





template<typename T, typename Tint, typename Tgroup, typename Fvertex, typename Fisom>
void LORENTZ_RunEdgewalkAlgorithm_Kernel(MyMatrix<T> const& G, std::vector<T> const& l_norms, FundDomainVertex<T,Tint> const& eVert, Fvertex f_vertex, Fisom f_isom)
{
  std::vector<int> l_status;
  std::vector<FundDomainVertex_FullInfo<T,Tint,Tgroup>> l_orbit_vertices;
  size_t nbDone = 0;
  auto func_insert_vertex=[&](FundDomainVertex_FullInfo<T,Tint,Tgroup> & vertFull1) -> bool {
    size_t len = l_orbit_vertices.size();
    for (size_t i=0; i<len; i++) {
      const FundDomainVertex_FullInfo<T,Tint,Tgroup>& vertFull2 = l_orbit_vertices[ i ];
      std::cerr << "i=" << i << "/" << len << " vert1=" << StringVectorGAP(vertFull1.vert.gen) << " / " << StringVectorGAP(vertFull2.vert.gen) << "\n";
      std::cerr << "    hash1=" << vertFull1.hash << " hash2=" << vertFull2.hash << "\n";
      if (vertFull1.hash == vertFull2.hash) {
        std::optional<MyMatrix<T>> equiv_opt = LinPolytopeIntegralWMat_Isomorphism<T,Tgroup,std::vector<T>,uint16_t>(vertFull1.e_pair_char, vertFull2.e_pair_char);
        if (equiv_opt) {
          std::cerr << "Find some isomorphism\n";
          return f_isom(UniversalMatrixConversion<Tint,T>(*equiv_opt));
        }
      }
    }
    std::cerr << "Failed to find some isomorphism\n";
    const auto& epair = vertFull1.e_pair_char;
    std::cerr << "GAP : MatV=" << StringMatrixGAP(epair.first) << " WMat=\n";
    PrintWeightedMatrixGAP(std::cerr, epair.second);
    std::cerr << "\n";
    std::cerr << "MatV=\n";
    WriteMatrix(std::cerr, epair.first);
    std::cerr << "WMat=\n";
    PrintWeightedMatrix(std::cerr, epair.second);
    std::cerr << "Before the LinPolytopeIntegralWMat_Automorphism nbDone=" << nbDone << " |l_orbit_vertices|=" << l_orbit_vertices.size() << "\n";
    for (auto & eGen : LinPolytopeIntegralWMat_Automorphism<T,Tgroup,std::vector<T>,uint16_t>(vertFull1.e_pair_char)) {
      bool test = f_isom(UniversalMatrixConversion<Tint,T>(eGen));
      if (test)
        return true;
    }
    bool test = f_vertex(vertFull1);
    if (test)
      return true;
    l_status.push_back(1);
    l_orbit_vertices.emplace_back(std::move(vertFull1));
    std::cerr << "Exiting the func_insert_vertex\n";
    return false;
  };
  // We have to do a copy of the Vert since when the vector is extended the previous entries are desttroyed when a new
  // array is built. This would then invalidates a const& theVert reference.
  // Took 1 week to fully debug that problem.
  auto insert_adjacent_vertices=[&](FundDomainVertex_FullInfo<T,Tint,Tgroup> const& vertFull) -> bool {
    const FundDomainVertex<T,Tint>& theVert = vertFull.vert;
    std::cerr << "insert_edges_from_vertex theVert=" << StringVectorGAP(RemoveFractionVector(theVert.gen)) << "\n";
    MyMatrix<T> FAC = UniversalMatrixConversion<T,Tint>(theVert.MatRoot);
    MyMatrix<T> FACred = ColumnReduction(FAC);
    vectface vf = lrs::DualDescription_temp_incd(FACred);
    //
    for (auto & eFAC : vf) {
      AdjacencyDirection<Tint> ad = GetAdjacencyDirection(theVert.MatRoot, eFAC);
      FundDomainVertex<T,Tint> fVert = EdgewalkProcedure(G, l_norms, theVert.gen, ad);
      { // Output. Fairly important to see what is happening
        T norm = fVert.gen.dot(G * fVert.gen);
        std::cerr << "Result of EdgewalkProcedure\n";
        std::cerr << "k=" << StringVectorGAP(theVert.gen) << " l_ui=";
        PrintAdjacencyDirection(std::cerr, ad);
        std::cerr << " fVert=" << StringVectorGAP(fVert.gen) << " norm=" << norm << "\n";
        std::cerr << "rec(k1:=" << StringVectorGAP(theVert.gen) << ",  k2:=" << StringVectorGAP(fVert.gen) << "),\n";
      }
      FundDomainVertex_FullInfo<T,Tint,Tgroup> fVertFull = gen_fund_domain_fund_info<T,Tint,Tgroup>(G, l_norms, fVert);
      bool test = func_insert_vertex(fVertFull);
      if (test)
        return true;
    }
    std::cerr << "Exiting from the insert_edges_from_vertex\n";
    return false;
  };
  FundDomainVertex_FullInfo<T,Tint,Tgroup> eVertFull = gen_fund_domain_fund_info<T,Tint,Tgroup>(G, l_norms, eVert);
  bool test = func_insert_vertex(eVertFull);
  if (test)
    return;
  while(true) {
    bool IsFinished = true;
    size_t len = l_status.size();
    for (size_t i=0; i<len; i++) {
      if (l_status[i] == 1) {
        nbDone++;
        IsFinished = false;
        l_status[i] = 0;
        // We need to do a direct copy in that case.
        // This is because we are working with a std::vector, which contains a T* array.
        // That array got deallocated and reallocated. This invalidates the const& reference.
        //
        // See for details https://stackoverflow.com/questions/6438086/iterator-invalidation-rules-for-c-containers
        // Which writes: "vector: all iterators and references before the point of insertion are unaffected, unless
        // the new container size is greater than the previous capacity (in which case all iterators and references are
        // invalidated) [23.2.4.3/1]
        //
        // The original problem originally took one week to debug.
        FundDomainVertex_FullInfo<T,Tint,Tgroup> VertFullCp = DirectCopy(l_orbit_vertices[i]);
        bool test = insert_adjacent_vertices(VertFullCp);
        if (test)
          return;
      }
    }
    if (IsFinished)
      break;
  }
  std::cerr << "Exiting from the infinite loop of enumeration of vertex pairs\n";
}

template<typename T, typename Tint, typename Tgroup>
ResultEdgewalk<T,Tint> LORENTZ_RunEdgewalkAlgorithm(MyMatrix<T> const& G, std::vector<T> const& l_norms, FundDomainVertex<T,Tint> const& eVert)
{
  std::vector<FundDomainVertex<T,Tint>> l_orbit_vertices_ret;
  auto f_vertex=[&](FundDomainVertex_FullInfo<T,Tint,Tgroup> const& vertFull) -> bool {
    l_orbit_vertices_ret.push_back(vertFull.vert);
    return false;
  };
  std::unordered_set<MyMatrix<Tint>> s_gen_isom_cox;
  MyMatrix<Tint> IdMat = IdentityMat<Tint>(G.rows());
  auto f_isom=[&](MyMatrix<Tint> const& eP) -> bool {
    if (eP == IdMat)
      return false;
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
    return false;
  };
  LORENTZ_RunEdgewalkAlgorithm_Kernel<T,Tint,Tgroup,decltype(f_vertex),decltype(f_isom)>(G, l_norms, eVert, f_vertex, f_isom);
  std::vector<MyMatrix<Tint>> l_gen_isom_cox;
  for (auto & e_gen : s_gen_isom_cox)
    l_gen_isom_cox.push_back(e_gen);
  return {std::move(l_gen_isom_cox), std::move(l_orbit_vertices_ret)};
}


template<typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>> LORENTZ_RunEdgewalkAlgorithm_Isomorphism(MyMatrix<T> const& G1, MyMatrix<T> const& G2, std::vector<T> const& l_norms, FundDomainVertex<T,Tint> const& eVert1, FundDomainVertex<T,Tint> const& eVert2)
{
  std::optional<MyMatrix<Tint>> answer;
  //
  FundDomainVertex_FullInfo<T,Tint,Tgroup> vertFull2 = gen_fund_domain_fund_info<T,Tint,Tgroup>(G1, l_norms, eVert2);
  auto f_vertex=[&](FundDomainVertex_FullInfo<T,Tint,Tgroup> const& vertFull1) -> bool {
    if (vertFull1.hash == vertFull2.hash) {
      std::optional<MyMatrix<T>> equiv_opt = LinPolytopeIntegralWMat_Isomorphism<T,Tgroup,std::vector<T>,uint16_t>(vertFull1.e_pair_char, vertFull2.e_pair_char);
      if (equiv_opt) {
        answer = UniversalMatrixConversion<Tint,T>(*equiv_opt);
        return true;
      }
    }
    return false;
  };
  auto f_isom=[&](MyMatrix<Tint> const& eP) -> bool {
    return false;
  };
  LORENTZ_RunEdgewalkAlgorithm_Kernel<T,Tint,Tgroup,decltype(f_vertex),decltype(f_isom)>(G1, l_norms, eVert1, f_vertex, f_isom);
  return answer;
}






template<typename T, typename Tint>
std::vector<MyVector<Tint>> get_simple_cone_from_lattice(MyMatrix<T> const& G, std::vector<T> const& l_norms, MyMatrix<Tint> const& NSP_tint)
{
  std::cerr << "Beginning of get_simple_cone\n";
  int dimSpace = NSP_tint.rows();
  MyMatrix<T> NSP = UniversalMatrixConversion<T,Tint>(NSP_tint);
  MyMatrix<Tint> G_int = UniversalMatrixConversion<Tint,T>(G);
  std::vector<MyVector<Tint>> l_roots;
  MyVector<T> zeroVect = ZeroVector<T>(dimSpace);
  //  std::cerr << "NSP=\n";
  //  WriteMatrix(std::cerr, NSP);
  for (auto & e_norm : l_norms) {
    std::cerr << "---------------------- e_norm=" << e_norm << " ----------------------\n";
    MyMatrix<T> Latt = ComputeLattice_LN(G, e_norm);
    //    std::cerr << "Latt=\n";
    //    WriteMatrix(std::cerr, Latt);
    //    std::cerr << "|Latt|=" << Latt.rows() << " / " << Latt.cols() << "\n";
    //    std::cerr << "|NSP|=" << NSP.rows() << " / " << NSP.cols() << "\n";
    MyMatrix<T> Latt_i_Orth = IntersectionLattice(NSP, Latt);
    //    std::cerr << "Latt_i_Orth=\n";
    //    WriteMatrix(std::cerr, Latt_i_Orth);
    //    std::cerr << "We have Latt_i_Orth\n";
    MyMatrix<Tint> Latt_i_Orth_tint = UniversalMatrixConversion<Tint,T>(Latt_i_Orth);
    MyMatrix<T> G_P = Latt_i_Orth * G * Latt_i_Orth.transpose();
    //    std::cerr << "G_P=\n";
    //    WriteMatrix(std::cerr, G_P);
    CheckPositiveDefinite(G_P);
    std::vector<MyVector<Tint>> l_v = FindFixedNormVectors<T,Tint>(G_P, zeroVect, e_norm);
    for (auto & e_v : l_v) {
      MyVector<Tint> e_root = Latt_i_Orth_tint.transpose() * e_v;
      std::optional<MyVector<Tint>> opt = SolutionIntMat(NSP_tint, e_root);
      if (opt) {
        l_roots.push_back(*opt);
      } else {
        std::cerr << "Failed to find the solution in the subspace\n";
        throw TerminalException{1};
      }
    }
    std::cerr << "e_norm=" << e_norm << " |l_v|=" << l_v.size() << "\n";
  }
  std::vector<MyVector<Tint>> facet_one_cone = GetFacetOneDomain(l_roots);
  if (facet_one_cone.size() != size_t(dimSpace)) {
    std::cerr << "dimSpace =" << dimSpace << " |facet_one_cone|=" << facet_one_cone.size() << "\n";
    std::cerr << "and they should be equal\n";
    throw TerminalException{1};
  }
  std::vector<MyVector<Tint>> l_ui;
  for (auto & e_root : facet_one_cone) {
    MyVector<Tint> v = NSP_tint.transpose() * e_root;
    l_ui.push_back(v);
  }
  return l_ui;
}



template<typename T, typename Tint>
MyMatrix<Tint> get_simple_cone(MyMatrix<T> const& G, std::vector<T> const& l_norms, MyVector<T> const& V)
{
  std::cerr << "------------------------------ get_simple_cone --------------------------\n";
  std::cerr << "G=\n";
  WriteMatrixGAP(std::cerr, G);
  std::cerr << "\n";
  T norm = V.dot(G*V);
  std::cerr << "V=" << StringVectorGAP(V) << " norm=" << norm << "\n";
  if (norm > 0) {
    std::cerr << "We need a vector of negative norm or zero norm in order to build a system of simple roots\n";
    throw TerminalException{1};
  }
  int dim = G.rows();
  MyVector<T> eProd = G * V;
  MyMatrix<T> eProdB(1,dim);
  AssignMatrixRow(eProdB, 0, eProd);
  MyMatrix<T> NSP = NullspaceIntTrMat(eProdB);
  MyMatrix<Tint> NSP_tint = UniversalMatrixConversion<Tint,T>(NSP);
  if (norm < 0) {
    // ordinary point case
    std::vector<MyVector<Tint>> l_vect = get_simple_cone_from_lattice(G, l_norms, NSP_tint);
    return MatrixFromVectorFamily(l_vect);
  } else {
    std::cerr << "get_simple_cone, step 1\n";
    // ideal point case
    MyVector<Tint> V_i = UniversalVectorConversion<Tint,T>(RemoveFractionVector(V));
    std::optional<MyVector<Tint>> opt = SolutionIntMat(NSP_tint, V_i);
    if (!opt) {
      std::cerr << "The vector V does not below to NSP which contradicts it being isotrop\n";
      throw TerminalException{1};
    }
    MyVector<Tint> const& Vnsp = *opt;
    std::cerr << "get_simple_cone, step 2\n";
    /*
      We need a more general code for finding complement of subspace, possibly using HermiteNormalForm
     */
    MyMatrix<Tint> Basis = ComplementToBasis(Vnsp);
    //    MyMatrix<Tint> Basis_p_Vnsp = ConcatenateMatVec(Basis, Vnsp);
    //    std::cerr << "Det(Basis_p_Vnsp)=" << DeterminantMat(Basis_p_Vnsp) << "\n";
    MyMatrix<Tint> Basis_NSP = Basis * NSP_tint;
    MyMatrix<T> Subspace = UniversalMatrixConversion<T,Tint>(Basis_NSP);
    //    std::cerr << "Subspace=\n";
    //    WriteMatrix(std::cerr, Subspace);
    std::map<T, size_t> MapIdxFr;
    std::vector<LatticeProjectionFramework<T>> ListFr;
    MyVector<T> zeroVect = ZeroVector<T>(Subspace.rows());
    std::vector<MyVector<T>> list_vect;
    std::vector<MyVector<T>> list_vect_big;
    std::vector<T> list_norm;
    std::cerr << "get_simple_cone, step 3\n";
    size_t pos = 0;
    for (auto & e_norm : l_norms) {
      std::cerr << "e_norm=" << e_norm << "\n";
      MyMatrix<T> Latt = ComputeLattice_LN(G, e_norm);
      MyMatrix<T> Latt_inter_NSP = IntersectionLattice(Latt, NSP);
      //      std::cerr << "Latt=\n";
      //      WriteMatrix(std::cerr, Latt);
      LatticeProjectionFramework<T> fr(G, Subspace, Latt_inter_NSP);
      //      std::cerr << "We have fr\n";
      MapIdxFr[e_norm] = pos;
      ListFr.push_back(fr);
      //      std::cerr << "We have MapFr assigned\n";
      //
      MyMatrix<T> const& RelBasis = fr.BasisProj;
      MyMatrix<T> G_P = RelBasis * G * RelBasis.transpose();
      //      CheckPositiveDefinite(G_P);
      std::vector<MyVector<Tint>> l_v = FindFixedNormVectors<T,Tint>(G_P, zeroVect, e_norm);
      for (auto & e_v : l_v) {
        MyVector<T> e_vt = UniversalVectorConversion<T,Tint>(e_v);
        MyVector<T> e_vect = RelBasis.transpose() * e_vt;
        std::optional<MyVector<T>> opt = SolutionMat(Subspace, e_vect);
        if (opt) {
          list_vect.push_back(*opt);
          list_vect_big.push_back(e_vect);
          list_norm.push_back(e_norm);
        } else {
          std::cerr << "Failed to find the solution in the subspace\n";
          throw TerminalException{1};
        }
      }
      pos++;
    }
    std::cerr << "|list_vect|=" << list_vect.size() << "\n";
    for (size_t i=0; i<list_vect.size(); i++)
      std::cerr << "i=" << i << " e_vect=" << StringVectorGAP(list_vect[i]) << " norm=" << list_norm[i] << "\n";
    std::cerr << "get_simple_cone, step 4\n";
    auto get_one_root=[&](MyVector<T> const& e_vect) -> MyVector<Tint> {
      std::cerr << "Beginning of get_one_root\n";
      std::cerr << "e_vect=" << StringVectorGAP(e_vect) << "\n";
      size_t len = list_vect.size();
      for (size_t i=0; i<len; i++) {
        MyVector<T> const& f_vect = list_vect[i];
        if (f_vect == e_vect) {
          T const& e_norm = list_norm[i];
          MyVector<T> const& e_vect_big = list_vect_big[i];
          size_t idx = MapIdxFr[e_norm];
          std::cerr << "e_norm=" << e_norm << " idx=" << idx << "\n";
          LatticeProjectionFramework<T> const& fr = ListFr[idx];
          std::optional<MyVector<T>> opt = fr.GetOnePreimage(e_vect_big);
          if (!opt) {
            std::cerr << "Failed to find the Preimage\n";
            throw TerminalException{1};
          }
          MyVector<T> const& V = *opt;
          MyVector<Tint> V_i = UniversalVectorConversion<Tint,T>(V);
          return V_i;
        }
      }
      std::cerr << "Failed to find the vector in the list\n";
      throw TerminalException{1};
    };
    std::vector<MyVector<T>> facet_one_cone = GetFacetOneDomain(list_vect);
    std::cerr << "get_simple_cone, step 5\n";
    std::vector<MyVector<Tint>> l_ui;
    for (auto & e_vt : facet_one_cone) {
      MyVector<Tint> e_vi = get_one_root(e_vt);
      std::cerr << "e_vt=" << StringVectorGAP(e_vt) << " e_vi=" << StringVectorGAP(e_vi) << "\n";
      l_ui.push_back(e_vi);
    }
    std::cerr << "get_simple_cone, step 5\n";
    MyMatrix<T> Pplane = Get_Pplane(G, l_ui);
    auto get_kP=[&]() -> MyVector<T> {
      MyMatrix<T> Gprod = Pplane * G * Pplane.transpose();
      T CritNorm = 0;
      bool StrictIneq = true;
      bool NeedNonZero = true;
      MyVector<Tint> eVect_A = GetShortVector_unlimited_float<Tint,T>(Gprod, CritNorm, StrictIneq, NeedNonZero);
      MyVector<T> eVect_B = UniversalVectorConversion<T,Tint>(eVect_A);
      MyVector<T> eVect_C = Pplane.transpose() * eVect_B;
      T scal = V.dot(G * eVect_C);
      if (scal < 0) { // This is because of the sign convention
        return eVect_C;
      } else {
        return -eVect_C;
      }
    };
    MyVector<T> kP = get_kP();
    std::cerr << "get_simple_cone, step 7\n";
    std::vector<MyVector<Tint>> l_vect = DetermineRootsCuspidalCase(G, l_ui, l_norms, V, kP);
    std::cerr << "get_simple_cone, step 8\n";
    return MatrixFromVectorFamily(l_vect);
  }
}



template<typename T, typename Tint>
MyVector<T> GetOneVertex(MyMatrix<T> const& G, std::vector<T> const& l_norms, bool const& ApplyReduction)
{
  ResultReductionIndefinite<T,Tint> ResRed = ComputeReductionIndefinite_opt<T,Tint>(G, ApplyReduction);
  /*
    We have ResRed.B and ResRed.Mred    with Mred = B * G * B^T
  */
  VinbergTot<T,Tint> Vtot = GetVinbergFromG<T,Tint>(ResRed.Mred, l_norms);
  MyVector<Tint> V1 = FindOneInitialRay(Vtot);
  MyVector<Tint> V2 = ResRed.B.transpose() * V1;
  MyVector<Tint> V3 = RemoveFractionVector(V2);
  MyVector<T> V4 = UniversalVectorConversion<T,Tint>(V3);
  return V4;
}




template<typename T, typename Tint>
FundDomainVertex<T,Tint> get_initial_vertex(MyMatrix<T> const& G, std::vector<T> const& l_norms, bool const& ApplyReduction, std::string const& OptionInitialVertex, std::string const& FileInitialVertex)
{
  std::cerr << "Beginning of get_initial_vertex\n";
  if (OptionInitialVertex == "FileVertex") {
    if (!IsExistingFile(FileInitialVertex)) {
      std::cerr << "The file FileInitialVertex=" << FileInitialVertex << " is missing\n";
      throw TerminalException{1};
    }
    std::ifstream is(FileInitialVertex);
    MyVector<T> gen = ReadVector<T>(is);
    MyMatrix<Tint> MatRoot = get_simple_cone<T,Tint>(G, l_norms, gen);
    return {RemoveFractionVector(gen), MatRoot};
  }
  if (OptionInitialVertex == "FileVertexRoots") {
    if (!IsExistingFile(FileInitialVertex)) {
      std::cerr << "The file FileInitialVertex=" << FileInitialVertex << " is missing\n";
      throw TerminalException{1};
    }
    std::ifstream is(FileInitialVertex);
    MyVector<T> gen = ReadVector<T>(is);
    MyMatrix<Tint> MatRoot = ReadMatrix<Tint>(is);
    return {RemoveFractionVector(gen), MatRoot};
  }
#ifdef ALLOW_VINBERG_ALGORITHM_FOR_INITIAL_VERTEX
  if (OptionInitialVertex == "vinberg") {
    MyVector<T> V = GetOneVertex<T,Tint>(G, l_norms, ApplyReduction);
    MyMatrix<Tint> MatRoot = get_simple_cone<T,Tint>(G, l_norms, V);
    return {RemoveFractionVector(V), MatRoot};
  }
#endif
  std::cerr << "Failed to find a matching entry in get_initial_vertex\n";
  std::cerr << "OptionInitialVertex=" << OptionInitialVertex << " but allowed values are FileVertex or FileVertexRoots\n";
#ifdef ALLOW_VINBERG_ALGORITHM_FOR_INITIAL_VERTEX
  std::cerr << "and vinberg has also been allowed\n";
#else
  std::cerr << "option vinberg has not been allowed\n";
#endif
  throw TerminalException{1};
}


template<typename T>
void TestLorentzianity(MyMatrix<T> const& G)
{
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(G);
  if (DiagInfo.nbZero != 0 || DiagInfo.nbMinus != 1) {
    std::cerr << "G=\n";
    WriteMatrix(std::cerr, G);
    std::cerr << "We have nbZero=" << DiagInfo.nbZero << " nbPlus=" << DiagInfo.nbPlus << " nbMinus=" << DiagInfo.nbMinus << "\n";
    std::cerr << "In the hyperbolic geometry we should have nbZero=0 and nbMinus=1\n";
    throw TerminalException{1};
  }
}


template<typename T, typename Tint, typename Tgroup>
void MainFunctionEdgewalk(FullNamelist const& eFull)
{
  SingleBlock BlockPROC=eFull.ListBlock.at("PROC");
  std::string FileLorMat=BlockPROC.ListStringValues.at("FileLorMat");
  MyMatrix<T> G = ReadMatrixFile<T>(FileLorMat);
  TestLorentzianity(G);
  //
  std::string OptionNorms=BlockPROC.ListStringValues.at("OptionNorms");
  bool ApplyReduction=BlockPROC.ListBoolValues.at("ApplyReduction");
  std::vector<T> l_norms = get_initial_list_norms<T,Tint>(G, OptionNorms);
  std::cerr << "We have l_norms\n";
  //
  std::string OptionInitialVertex=BlockPROC.ListStringValues.at("OptionInitialVertex");
  std::string FileInitialVertex=BlockPROC.ListStringValues.at("FileInitialVertex");
  FundDomainVertex<T,Tint> eVert = get_initial_vertex<T,Tint>(G, l_norms, ApplyReduction, OptionInitialVertex, FileInitialVertex);
  T norm = eVert.gen.dot(G * eVert.gen);
  std::cerr << "Initial vertex is eVert=" << StringVectorGAP(eVert.gen) << " norm=" << norm << "\n";
  std::cerr << "|MatRoot|=" << eVert.MatRoot.rows() << "\n";
  std::vector<MyVector<Tint>> l_root;
  for (int i=0; i<eVert.MatRoot.rows(); i++) {
    MyVector<Tint> eLine = GetMatrixRow(eVert.MatRoot, i);
    l_root.push_back(eLine);
  }
  MyMatrix<T> CoxMat = ComputeCoxeterMatrix(G, l_root).first;
  std::cerr << "We have CoxMat\n";
  std::string symb = coxdyn_matrix_to_string(CoxMat);
  std::cerr << "symb=" << symb << "\n";
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





template<typename T, typename Tint, typename Tgroup>
void MainFunctionEdgewalk_Isomorphism(FullNamelist const& eFull)
{
  SingleBlock BlockPROC=eFull.ListBlock.at("PROC");
  std::string FileLorMat1=BlockPROC.ListStringValues.at("FileLorMat1");
  std::string FileLorMat2=BlockPROC.ListStringValues.at("FileLorMat2");
  MyMatrix<T> G1 = ReadMatrixFile<T>(FileLorMat1);
  MyMatrix<T> G2 = ReadMatrixFile<T>(FileLorMat2);
  TestLorentzianity(G1);
  TestLorentzianity(G2);
  //
  auto print_result=[&](std::optional<MyMatrix<Tint>> const& opt) -> void {
    std::string OutFormat=BlockPROC.ListStringValues.at("OutFormat");
    auto print_result_isomorphism=[&](std::ostream& os) -> void {
      if (OutFormat == "GAP") {
        if (opt) {
          os << "return ";
          WriteMatrixGAP(os, *opt);
          os << ";\n";
        } else {
          os << "return fail;\n";
        }
      }
      std::cerr << "We fil to have a matching format. OutFormat=" << OutFormat << "\n";
      throw TerminalException{1};
    };
    std::string FileOut=BlockPROC.ListStringValues.at("FileOut");
    if (FileOut == "stderr") {
      print_result_isomorphism(std::cerr);
    } else {
      if (FileOut == "stdout") {
        print_result_isomorphism(std::cout);
      } else {
        std::ofstream os(FileOut);
        print_result_isomorphism(os);
      }
    }
  };
  //
  std::string OptionNorms="all";
  bool ApplyReduction=BlockPROC.ListBoolValues.at("ApplyReduction");
  std::vector<T> l_norms1 = get_initial_list_norms<T,Tint>(G1, OptionNorms);
  std::vector<T> l_norms2 = get_initial_list_norms<T,Tint>(G2, OptionNorms);
  if (l_norms1 != l_norms2) {
    print_result( {} );
    return;
  }
  std::vector<T> l_norms = l_norms1;
  std::cerr << "We have l_norms\n";
  //
  std::string OptionInitialVertex="vinberg";
  std::string FileInitialVertex="irrelevant";
  FundDomainVertex<T,Tint> eVert1 = get_initial_vertex<T,Tint>(G1, l_norms, ApplyReduction, OptionInitialVertex, FileInitialVertex);
  FundDomainVertex<T,Tint> eVert2 = get_initial_vertex<T,Tint>(G2, l_norms, ApplyReduction, OptionInitialVertex, FileInitialVertex);
  //
  std::optional<MyMatrix<Tint>> opt = LORENTZ_RunEdgewalkAlgorithm_Isomorphism<T,Tint,Tgroup>(G1, G2, l_norms, eVert1, eVert2);
  print_result(opt);
}



#endif

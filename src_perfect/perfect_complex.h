// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_PERFECT_COMPLEX_H_
#define SRC_PERFECT_PERFECT_COMPLEX_H_

#include "perfect_tspace.h"
#include "triples.h"

#ifdef SANITY_CHECK
#define SANITY_CHECK_PERFECT_COMPLEX
#endif

#ifdef DEBUG
#define DEBUG_PERFECT_COMPLEX
#endif

#ifdef PRINT
#define PRINT_PERFECT_COMPLEX_DIMENSION
#endif

/*
  Construction of the perfect form complex.
  We compute the faces of the oerfect domain
  and check for equivalence.

  There are two options.
  + Computes all the faces.
  + Computes only the well rounded faces.

  Basic algorithm:
  + From the perfect cone, build the facets.
  + Check the facets for equivalence.
  + After one level is computed, go to the
    next one.
  + We also need to have the stabilizers
    . To know which are orientable for the
      homology group.
    . To compute the lower faces.

  The main difficulty is the test of equivalence.
  There are two methods:
  + For the well rounded faces, computing
    could be done with the classic algorithm.
    It has some advantages, but it is not
    necessarily the fastest algorithm since
    computing equivalences requires a partition
    backtrack call. Other inconveniences:
    . It does not work with degenerate faces.
    . The quadratic form used for the computation
      is not necessarily well defined as we
      experienced when working on that issue
      in GAP.
  + For the other faces we can apply the
    algorithm of using the top dimensional
    faces.
    That algorithm has several advantages:
    . We can work with degenerate faces.
    . We can work with stabilizers having
      infinite size.
    . No backtracking algorithm is required
      for computing equivalence.
    . The inconvenience is that we would need
      to store all the triples.

  Therefore, we use the second algorithm

  - - - - - - - - - - - - -

  
  
 */

//
// The top dimensional complex
//

template<typename T, typename Tint_impl, typename Tgroup_impl>
struct PerfectFormInfoForComplex {
  using Tint = Tint_impl;
  using Tgroup = Tgroup_impl;
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> gram;
  MyMatrix<Tint> EXT; // Not necessarily full-dimensional
  std::optional<permutalib::PreImagerElement<Telt, MyMatrix<Tint>, TintGroup>> opt_pre_imager;
  Tgroup GRP_ext; // Group acting on the shortest vectors
  std::vector<sing_adj<Tint>> l_sing_adj;
  MyMatrix<Tint> find_matrix(Telt const& x, [[maybe_unused]] std::ostream& os) const {
    if (opt_pre_imager) {
      permutalib::PreImagerElement<Telt, MyMatrix<Tint>, TintGroup> const& pre_imager = *opt_pre_imager;
      std::optional<MyMatrix<Tint>> opt = pre_imager.get_preimage(x);
      MyMatrix<Tint> M = unfold_opt(opt, "The element elt should belong to the group");
#ifdef SANITY_CHECK_PERFECT_COMPLEX
      Telt x_img = get_elt_from_matrix<Tint,Telt>(M, EXT, os);
      if (x_img != x) {
        std::cerr << "Inconsistency in the pre_image computation\n";
        throw TerminalException{1};
      }
#endif
      return M;
    }
    return FindTransformation(EXT, EXT, x);
  }
};

template<typename T, typename Tint, typename Tgroup>
struct PerfectComplexTopDimInfo {
  std::vector<PerfectFormInfoForComplex<T,Tint,Tgroup>> l_perfect;
  LinSpaceMatrix<T> LinSpa;
  bool only_well_rounded;
  bool compute_boundary;
};


template<typename T, typename Tint, typename Tgroup>
PerfectComplexTopDimInfo<T,Tint,Tgroup> generate_perfect_complex_top_dim_info(std::vector<DatabaseEntry_Serial<PerfectTspace_Obj<T,Tint,Tgroup>,PerfectTspace_AdjO<Tint>>> const& l_tot, LinSpaceMatrix<T> const& LinSpa, bool const& only_well_rounded, bool const& compute_boundary, std::ostream& os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  std::vector<PerfectFormInfoForComplex<T,Tint,Tgroup>> l_perfect;
  for (auto & ePerf: l_tot) {
    std::vector<sing_adj<Tint>> l_sing_adj;
    for (auto & eAdj: ePerf.ListAdj) {
      size_t jCone = eAdj.iOrb;
      Face f_ext = eAdj.x.eInc;
      // This inverse is because the eBigMat is acting on perfect forms
      // but for the complex, we work on vector configurations.
      MyMatrix<Tint> eMat = Inverse(eAdj.x.eBigMat);
      sing_adj<Tint> adj{jCone, f_ext, eMat};
      l_sing_adj.emplace_back(std::move(adj));
    }
    MyMatrix<Tint> EXT = conversion_and_duplication<Tint, Tint>(ePerf.x.rec_shv.SHV);
    // In some cases, EXT is not full dimensional. We need an alternate strategy for that.
    std::optional<permutalib::PreImagerElement<Telt, MyMatrix<Tint>, TintGroup>> opt_pre_imager;
    if (RankMat(EXT) < EXT.cols()) {
      std::vector<MyMatrix<Tint>> const& l_matr = ePerf.x.GRP_matr;
      std::vector<Telt> l_perm = get_list_elt_from_list_matrices<Tint,Telt>(l_matr, EXT, os);
      MyMatrix<Tint> id = IdentityMat<Tint>(EXT.cols());
      opt_pre_imager = permutalib::PreImagerElement<Telt, MyMatrix<Tint>, TintGroup>(l_matr, l_perm, id);
    }
    Tgroup const& GRP_ext = ePerf.x.GRP;
    PerfectFormInfoForComplex<T,Tint,Tgroup> perfect{ePerf.x.Gram, std::move(EXT), opt_pre_imager, GRP_ext, std::move(l_sing_adj)};
    l_perfect.emplace_back(std::move(perfect));
  }
  return {std::move(l_perfect), LinSpa, only_well_rounded, compute_boundary};
}

struct OrientationInfo {
  int sign;
  std::vector<int> ListRowSelect;
  std::vector<int> ListColSelect;
};

//
// The intermediate dimensional faces.
// This is the orientation coming from the 
//

template<typename T, typename Tint, typename Tgroup>
struct FacePerfectComplex {
  std::vector<triple<Tint>> l_triple; // The containing triples.
  std::vector<MyMatrix<Tint>> l_gens; // generating set of the stabilizer
  MyMatrix<Tint> EXT; // The list of vectors of the perfect form.
  Tgroup GRP_ext; // Group acting on the vectors (NOT on the rays of the cone)
  bool is_well_rounded;
  OrientationInfo or_info; // The info for determining the orientation.
  bool is_orientable; // Whether the cell was orientable.
  MyVector<T> get_interior_point(std::vector<MyMatrix<T>> const& ListMat) const {
    int n_row = or_info.ListRowSelect.size();
    int n_col = or_info.ListColSelect.size();
    MyVector<T> interior_pt = ZeroVector<T>(n_col);
    for (int i_row=0; i_row<n_row; i_row++) {
      int j_row = or_info.ListRowSelect[i_row];
      MyVector<Tint> V = GetMatrixRow(EXT, j_row);
      for (int i_col=0; i_col<n_col; i_col++) {
        int i_mat = or_info.ListColSelect[i_col];
        T scal = EvaluationQuadForm<T,Tint>(ListMat[i_mat], V);
        interior_pt(i_col) += scal;
      }
    }
    return interior_pt;
  }
};


template<typename T, typename Tint, typename Tgroup>
struct FacesPerfectComplex {
  std::vector<FacePerfectComplex<T,Tint,Tgroup>> l_faces;
};


template<typename T, typename Tint>
OrientationInfo get_orientation_info(MyMatrix<Tint> const& EXT, std::vector<MyMatrix<T>> const& ListMat) {
  int n_row = EXT.rows();
  int n_mat = ListMat.size();
  MyMatrix<T> MatScal(n_row, n_mat);
  for (int i_row=0; i_row<n_row; i_row++) {
    MyVector<Tint> V = GetMatrixRow(EXT, i_row);
    for (int i_mat=0; i_mat<n_mat; i_mat++) {
      T scal = EvaluationQuadForm<T,Tint>(ListMat[i_mat], V);
      MatScal(i_row, i_mat) = scal;
    }
  }
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(MatScal);
  int dim = eSelect.ListRowSelect.size();
  MyMatrix<T> M(dim, dim);
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      int i2 = eSelect.ListRowSelect[i];
      int j2 = eSelect.ListColSelect[j];
      M(i, j) = MatScal(i2, j2);
    }
  }
  T det = DeterminantMat(M);
  int sign = 1;
  if (det < 0) {
    sign = -1;
  }
  return {sign, std::move(eSelect.ListRowSelect), std::move(eSelect.ListColSelect)};
}





template<typename T, typename Tint>
int get_face_orientation(MyMatrix<Tint> const& EXT, std::vector<MyMatrix<T>> const& ListMat,
                         OrientationInfo const& or_info, MyMatrix<Tint> const& t) {
#ifdef SANITY_CHECK_PERFECT_COMPLEX
  int n_row = EXT.rows();
  std::unordered_set<MyVector<Tint>> set;
  for (int i_row=0; i_row<n_row; i_row++) {
    MyVector<Tint> V = GetMatrixRow(EXT, i_row);
    set.insert(V);
  }
  for (int i_row=0; i_row<n_row; i_row++) {
    MyVector<Tint> V = GetMatrixRow(EXT, i_row);
    MyVector<Tint> Vimg = t.transpose() * V;
    if (set.count(Vimg) != 1) {
      std::cerr << "PERFCOMP: Vimg should belong to V\n";
      throw TerminalException{1};
    }
  }
#endif
  int dim = or_info.ListRowSelect.size();
  MyMatrix<T> M(dim, dim);
  for (int i=0; i<dim; i++) {
    int i_row = or_info.ListRowSelect[i];
    MyVector<Tint> V = GetMatrixRow(EXT, i_row);
    MyVector<Tint> Vimg = t.transpose() * V;
    for (int j=0; j<dim; j++) {
      int i_mat = or_info.ListColSelect[j];
      T scal = EvaluationQuadForm<T,Tint>(ListMat[i_mat], Vimg);
      M(i, j) = scal;
    }
  }
  T det = DeterminantMat(M);
  if (det > 0) {
    return or_info.sign;
  } else {
    return -or_info.sign;
  }
}

template<typename T, typename Tint>
bool is_face_orientable(MyMatrix<Tint> const& EXT, std::vector<MyMatrix<T>> const& ListMat,
                        OrientationInfo const& or_info, std::vector<MyMatrix<Tint>> const& l_t) {
  for (auto & t: l_t) {
    int sign = get_face_orientation(EXT, ListMat, or_info, t);
    if (sign == -1) {
      return false;
    }
  }
  return true;
}



template<typename T, typename Tint, typename Tgroup>
FacesPerfectComplex<T,Tint,Tgroup> get_first_step_perfect_complex_enumeration(PerfectComplexTopDimInfo<T,Tint,Tgroup> const& pctdi, std::ostream & os) {
  using Telt = typename Tgroup::Telt;
  std::vector<FacePerfectComplex<T,Tint,Tgroup>> l_faces;
  int n = pctdi.LinSpa.n;
#ifdef DEBUG_PERFECT_COMPLEX
  os << "PERFCOMP: get_first_step_perfect_complex_enumeration n=" << n << "\n";
#endif
  auto get_l_triple=[&](size_t i_domain) -> std::vector<triple<Tint>> {
    int n_ext = pctdi.l_perfect[i_domain].EXT.rows();
    Face f_ext(n_ext);
    for (int i_ext=0; i_ext<n_ext; i_ext++) {
      f_ext[i_ext] = 1;
    }
    MyMatrix<Tint> eMat = IdentityMat<Tint>(n);
    triple<Tint> t{i_domain, f_ext, eMat};
    return {t};
  };
  auto get_l_gens=[&](size_t i_domain) -> std::vector<MyMatrix<Tint>> {
    std::vector<MyMatrix<Tint>> l_gens;
    PerfectFormInfoForComplex<T,Tint,Tgroup> const& top = pctdi.l_perfect[i_domain];
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: |top.GRP_ext|=" << top.GRP_ext.size() << "\n";
#endif
    std::vector<Telt> l_elt = top.GRP_ext.SmallGeneratingSet();
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: |l_elt|=" << l_elt.size() << "\n";
#endif
    for (auto & ePermGen: l_elt) {
      MyMatrix<Tint> eMatrGen = top.find_matrix(ePermGen, os);
      l_gens.push_back(eMatrGen);
    }
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: Returning l_gens\n";
#endif
    return l_gens;
  };
  auto get_or_info=[&](size_t i_domain) -> OrientationInfo {
    MyMatrix<Tint> const& EXT = pctdi.l_perfect[i_domain].EXT;
    std::vector<MyMatrix<T>> const& ListMat = pctdi.LinSpa.ListMat;
    return get_orientation_info(EXT, ListMat);
  };
  auto get_is_orientable=[&](size_t i_domain, std::vector<MyMatrix<Tint>> const& l_gens, OrientationInfo const& or_info) -> bool {
    MyMatrix<Tint> const& EXT = pctdi.l_perfect[i_domain].EXT;
    std::vector<MyMatrix<T>> const& ListMat = pctdi.LinSpa.ListMat;
    return is_face_orientable(EXT, ListMat, or_info, l_gens);
  };
  size_t n_domain = pctdi.l_perfect.size();
  os << "PERFCOMP: get_first_step n_domain=" << n_domain << "\n";
  for (size_t i_domain=0; i_domain<n_domain; i_domain++) {
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: get_first_step i_domain=" << i_domain << "\n";
#endif
    std::vector<triple<Tint>> l_triple = get_l_triple(i_domain);
    std::vector<MyMatrix<Tint>> l_gens = get_l_gens(i_domain);
    MyMatrix<Tint> const& EXT = pctdi.l_perfect[i_domain].EXT;
    Tgroup const& GRP_ext = pctdi.l_perfect[i_domain].GRP_ext;
    bool is_well_rounded = true; // Yes, the top dimensional cells are well rounded.
    OrientationInfo or_info = get_or_info(i_domain);
    bool is_orientable = get_is_orientable(i_domain, l_gens, or_info);
    FacePerfectComplex<T,Tint,Tgroup> face{l_triple, l_gens, EXT, GRP_ext, is_well_rounded, or_info, is_orientable};
    l_faces.push_back(face);
  }
  FacesPerfectComplex<T,Tint,Tgroup> pfc{l_faces};
  return pfc;
}

template<typename Tint>
struct BoundEntry {
  int iOrb;
  int sign;
  MyMatrix<Tint> M;
};

template<typename Tint>
struct ListBoundEntry {
  std::vector<BoundEntry<Tint>> l_bound;
};

template<typename Tint>
struct FullBoundary {
  std::vector<ListBoundEntry<Tint>> ll_bound;
};


template<typename T, typename Tint, typename Tgroup>
struct ResultStepEnumeration {
  FacesPerfectComplex<T,Tint,Tgroup> level;
  std::optional<FullBoundary<Tint>> boundary;
};

template<typename Tint, typename Telt>
struct ElementMapper {
private:
  MyMatrix<Tint> EXT;
  ContainerMatrix<Tint> cont;
public:
  ElementMapper(MyMatrix<Tint> const& _EXT) : EXT(_EXT), cont(EXT) {
  }
  Telt map_elt(MyMatrix<Tint> const& M) const {
    using Tidx = typename Telt::Tidx;
    int n_ext = EXT.rows();
    std::vector<Tidx> v(n_ext);
    for (int u=0; u<n_ext; u++) {
      MyVector<Tint> V1 = GetMatrixRow(EXT, u);
      MyVector<Tint> V2 = M.transpose() * V1;
      std::optional<size_t> opt = cont.GetIdx_v(V2);
      size_t pos = unfold_opt(opt, "PERFCOMPL: Failed to find the index");
      v[u] = pos;
    }
    Telt elt(v);
    return elt;
  }
};

template<typename Tint>
struct FaceEntry {
  int iOrb;
  Tint value;
  MyMatrix<Tint> M;
};

template<typename T, typename Tint, typename Tgroup>
ResultStepEnumeration<T,Tint,Tgroup> compute_next_level(PerfectComplexTopDimInfo<T,Tint,Tgroup> const& pctdi, FacesPerfectComplex<T,Tint,Tgroup> const& level, std::ostream & os) {
  using Telt = typename Tgroup::Telt;
  int n = pctdi.LinSpa.n;
#ifdef DEBUG_PERFECT_COMPLEX
  os << "PERFCOMP: compute_next_level, start, n=" << n << "\n";
#endif
  std::vector<FacePerfectComplex<T,Tint,Tgroup>> l_faces;
  std::vector<ListBoundEntry<Tint>> ll_bound;
  auto find_matching_entry=[&](triple<Tint> const& t) -> std::optional<std::pair<int,MyMatrix<Tint>>> {
    int i_domain = 0;
    for (auto & face1: l_faces) {
      std::optional<MyMatrix<Tint>> opt =
        test_triple_in_listtriple(pctdi.l_perfect, face1.l_triple, t, os);
      if (opt) {
        MyMatrix<Tint> const& M = *opt;
        std::pair<int,MyMatrix<Tint>> p{i_domain, M};
        return p;
      }
      i_domain += 1;
    }
    return {};
  };
  using Tfull_triple = std::pair<std::vector<triple<Tint>>, std::vector<MyMatrix<Tint>>>;
  auto need_opt_t=[&]([[maybe_unused]] PerfectBoundednessProperty const& pbp, [[maybe_unused]] triple<Tint> const& t) -> bool {
#ifdef SANITY_CHECK_PERFECT_COMPLEX_EXTENSIVE
    std::cerr << "Returning true because we want to compute whether the group is finite\n";
    Tfull_triple pair = get_spanning_list_triple(pctdi.l_perfect, t, os);
    std::vector<MyMatrix<Tint>> l_matr = pair.second;
    std::string Prefix = "FinitenessMatrixGroupTest_" + std::to_string(l_matr.size()) + "_n" + std::to_string(n) + "_";
    std::string FileOut = FindAvailableFileFromPrefix(Prefix);
    bool write_test_case = false;
    if (write_test_case) {
      bool is_finite = test_finiteness_group<T,Tint>(l_matr, os);
      std::ofstream os_out(FileOut);
      os_out << "return rec(GRPmatr:=";
      WriteListMatrixGAP(os_out, l_matr);
      os_out << ", is_finite:=" << GAP_logical(is_finite) << ");\n";
    } else {
      WriteListMatrixFile(FileOut, l_matr);
    }
    return true;
#else
    if (pbp.bounded_self_dual) {
      return false;
    }
    if (pbp.bounded_spanning) {
      return false;
    }
    return true;
#endif
  };
  auto f_is_well_rounded=[&](triple<Tint> const& t, std::optional<Tfull_triple> & opt_t, PerfectBoundednessProperty & pbp) -> bool {
    if (need_opt_t(pbp, t)) {
      Tfull_triple pair = get_spanning_list_triple(pctdi.l_perfect, t, os);
      bool is_finite = test_finiteness_group<T,Tint>(pair.second, os);
      pbp.bounded_finite_stabilizer = is_finite;
      opt_t = pair;
    }
    //#ifdef SANITY_CHECK_PERFECT_COMPLEX
    // That check is really cheap. Do it all the time.
    check_pbp(pbp);
    //#endif
    return get_result(pbp);
  };
  auto is_insertable=[&](bool is_well_rounded) -> bool {
    if (pctdi.only_well_rounded) {
      return is_well_rounded;
    }
    return true;
  };
  auto get_initial_triple=[&](FacePerfectComplex<T,Tint,Tgroup> const& face, Face const& eIncd) -> triple<Tint> {
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: get_initial_triple, step 1 |face.l_triple|=" << face.l_triple.size() << "\n";
#endif
    triple<Tint> const& t_big = face.l_triple[0];
    size_t iCone = t_big.iCone;
    int n_ext = pctdi.l_perfect[iCone].EXT.rows();
    Face f(n_ext);
    size_t index = 0;
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: get_initial_triple, step 6\n";
    os << "PERFCOMP: get_initial_triple, |t_big.f_ext|=" << t_big.f_ext.size() << " / " << t_big.f_ext.count() << "\n";
    os << "PERFCOMP: get_initial_triple, |eIncd|=" << eIncd.size() << " / " << eIncd.count() << "\n";
#endif
    for (int i=0; i<n_ext; i++) {
      if (t_big.f_ext[i] == 1) {
        if (eIncd[index] == 1) {
          f[i] = 1;
        }
        index += 1;
      }
    }
    triple<Tint> t{iCone, f, t_big.eMat};
    triple<Tint> t_can = canonicalize_triple(pctdi.l_perfect, t, os);
    return t_can;
  };
  auto f_insert=[&](triple<Tint> const& t, std::optional<Tfull_triple> & opt_t, bool const& is_well_rounded) -> std::pair<int, MyMatrix<Tint>> {
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: f_insert, step 1\n";
#endif
    std::optional<std::pair<int,MyMatrix<Tint>>> opt = find_matching_entry(t);
    if (opt) {
      return *opt;
    }
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: f_insert, step 3\n";
#endif
    if (!opt_t) {
      opt_t = get_spanning_list_triple(pctdi.l_perfect, t, os);
    }
    Tfull_triple const& pair = *opt_t;
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: f_insert, step 4\n";
#endif
    int iCone = t.iCone;
    int n_ext = t.f_ext.count();
    int n_ext_big = pctdi.l_perfect[iCone].EXT.rows();
    MyMatrix<Tint> EXT(n_ext, n);
    int pos = 0;
    for (int i_big=0; i_big<n_ext_big; i_big++) {
      if (t.f_ext[i_big] == 1) {
        MyVector<Tint> V1 = GetMatrixRow(pctdi.l_perfect[iCone].EXT, i_big);
        MyVector<Tint> V2 = t.eMat.transpose() * V1;
        AssignMatrixRow(EXT, pos, V2);
        pos += 1;
      }
    }
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: f_insert, step 5\n";
#endif
    ElementMapper<Tint,Telt> elt_mapper(EXT);
    std::vector<Telt> l_gens;
    for (auto & eMatrGen: pair.second) {
      Telt elt = elt_mapper.map_elt(eMatrGen);
      l_gens.emplace_back(std::move(elt));
    }
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: f_insert, step 7\n";
#endif
    Tgroup GRP_ext(l_gens, n_ext);
    OrientationInfo or_info = get_orientation_info(EXT, pctdi.LinSpa.ListMat);
    bool is_orientable = is_face_orientable(EXT, pctdi.LinSpa.ListMat, or_info, pair.second);
    FacePerfectComplex<T,Tint,Tgroup> face{
      pair.first,
      pair.second,
      EXT,
      GRP_ext,
      is_well_rounded,
      or_info,
      is_orientable};
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: f_insert, step 8\n";
#endif
    int i_domain = l_faces.size();
    l_faces.push_back(face);
    MyMatrix<Tint> M = IdentityMat<Tint>(n);
    std::pair<int, MyMatrix<Tint>> p{i_domain, M};
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: f_insert, step 9\n";
#endif
    return p;
  };
  MyVector<T> interior_pt;
  std::vector<Telt> l_gens_perm;
  std::vector<int> l_status;
#ifdef SANITY_CHECK_PERFECT_COMPLEX
  std::unordered_set<MyVector<Tint>> set_EXT;
#endif
  auto get_sign=[&](MyVector<T> const& interior_pt, FacePerfectComplex<T,Tint,Tgroup> const& face, std::pair<int, MyMatrix<Tint>> const& p) -> int {
    int dim = face.or_info.ListRowSelect.size();
    MyMatrix<T> M(dim, dim);
    for (int i=0; i<dim; i++) {
      M(dim-1, i) = interior_pt(i);
    }
    int iOrb = p.first;
    for (int u=0; u<dim-1; u++) {
      int i_row = l_faces[iOrb].or_info.ListRowSelect[u];
      MyVector<Tint> V = GetMatrixRow(l_faces[iOrb].EXT, i_row);
      MyVector<Tint> Vimg = p.second.transpose() * V;
#ifdef SANITY_CHECK_PERFECT_COMPLEX
      if (set_EXT.count(Vimg) != 1) {
        std::cerr << "PERFCOMP: The vector does not belong to set\n";
        throw TerminalException{1};
      }
#endif
      for (int i_col=0; i_col<dim; i_col++) {
        int i_mat = face.or_info.ListColSelect[i_col];
        T scal = EvaluationQuadForm<T,Tint>(pctdi.LinSpa.ListMat[i_mat], Vimg);
        M(u, i_col) = scal;
      }
    }
    T det = DeterminantMat(M);
    if (det > 0) {
      return face.or_info.sign;
    } else {
      return -face.or_info.sign;
    }
  };
  auto initial_boundary_setup=[&](size_t const& i, RyshkovGRP<T, Tgroup> const& cone) -> void {
    FacePerfectComplex<T,Tint,Tgroup> const& face = level.l_faces[i];
    interior_pt = face.get_interior_point(pctdi.LinSpa.ListMat);
    l_gens_perm.clear();
    ElementMapper<Tint,Telt> elt_mapper(face.EXT);
    for (auto & eMatrGen: face.l_gens) {
      Telt elt1 = elt_mapper.map_elt(eMatrGen);
      Telt elt2 = cone.map_elt(elt1);
      int sign = get_face_orientation(face.EXT, pctdi.LinSpa.ListMat, face.or_info, eMatrGen);
      l_gens_perm.push_back(elt2);
      l_status.push_back(sign);
    }
#ifdef SANITY_CHECK_PERFECT_COMPLEX
    for (int i_row=0; i_row<face.EXT.rows(); i_row++) {
      MyVector<Tint> V = GetMatrixRow(face.EXT, i_row);
      set_EXT.insert(V);
    }
#endif
  };
  auto append_boundary=[&](size_t const& i, std::pair<int,MyMatrix<Tint>> const& p, Face const& incd_sma, std::vector<BoundEntry<Tint>> & l_bound) -> void {
    FacePerfectComplex<T,Tint,Tgroup> const& face = level.l_faces[i];
    int sign = get_sign(interior_pt, face, p);
    std::unordered_set<Face> set_faces_gen;
    std::vector<Face> l_faces_gen;
    std::vector<int> l_sign;
    std::vector<MyMatrix<Tint>> l_mat;
    auto insert_entry=[&](Face const& new_incd_sma, int const& new_sign, MyMatrix<Tint> const& new_mat) -> void {
      if (set_faces_gen.count(new_incd_sma) == 0) {
        set_faces_gen.insert(new_incd_sma);
        l_faces_gen.push_back(new_incd_sma);
        l_sign.push_back(new_sign);
        l_mat.push_back(new_mat);
        BoundEntry<Tint> be{p.first, sign, new_mat};
#ifdef SANITY_CHECK_PERFECT_COMPLEX
        MyMatrix<Tint> EXTimg = l_faces[be.iOrb].EXT * be.M;
        for (int iRow=0; iRow<EXTimg.rows(); iRow++) {
          MyVector<Tint> Vimg = GetMatrixRow(EXTimg, iRow);
          if (set_EXT.count(Vimg) != 1) {
            std::cerr << "PERFCOMP: The vector should belong to the image\n";
            throw TerminalException{1};
          }
        }
#endif
        l_bound.push_back(be);
      }
    };
    insert_entry(incd_sma, sign, p.second);
    std::vector<MyMatrix<Tint>> const& l_gens = face.l_gens;
    size_t n_gen = l_gens.size();
    size_t start=0;
    while(true) {
      size_t len = l_faces_gen.size();
      for (size_t u=start; u<len; u++) {
        for (size_t i_gen=0; i_gen<n_gen; i_gen++) {
          Face new_incd_sma = OnFace(l_faces_gen[u], l_gens_perm[i_gen]);
          int new_sign = l_sign[u] * l_status[i_gen];
          MyMatrix<Tint> new_mat = l_mat[u] * l_gens[i_gen];
          insert_entry(new_incd_sma, new_sign, new_mat);
        }
      }
      start = len;
      if (start == l_faces_gen.size()) {
        break;
      }
    }
  };
  for (size_t i=0; i<level.l_faces.size(); i++) {
    FacePerfectComplex<T,Tint,Tgroup> const& face = level.l_faces[i];
    std::vector<BoundEntry<Tint>> l_bound;
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T,Tint>(face.EXT);
    RyshkovGRP<T, Tgroup> cone =
      GetNakedPerfectCone_GRP<T, Tgroup>(pctdi.LinSpa, EXT_T, face.GRP_ext, os);
    if (pctdi.compute_boundary) {
      initial_boundary_setup(i, cone);
    }
    for (auto& incd_sma: cone.ListIncd) {
      Face incd_big = get_big_incd(cone, incd_sma);
      MyMatrix<Tint> EXTincd = SelectRow(face.EXT, incd_big);
      PerfectBoundednessProperty pbp = initial_bounded_property(pctdi.LinSpa, EXTincd, os);
      std::optional<Tfull_triple> opt_t;
      triple<Tint> t = get_initial_triple(face, incd_big);
      bool is_well_rounded = f_is_well_rounded(t, opt_t, pbp);
#ifdef DEBUG_PERFECT_COMPLEX
      os << "PERFCOMP: pbp=" << to_string(pbp) << "\n";
#endif
      if (is_insertable(is_well_rounded)) {
        std::pair<int, MyMatrix<Tint>> p = f_insert(t, opt_t, is_well_rounded);
        if (pctdi.compute_boundary) {
          append_boundary(i, p, incd_sma, l_bound);
        }
      }
    }
    if (pctdi.compute_boundary) {
      ListBoundEntry<Tint> lbe{l_bound};
      ll_bound.push_back(lbe);
    }
  }
  FacesPerfectComplex<T,Tint,Tgroup> pfc{l_faces};
  if (pctdi.compute_boundary) {
    FullBoundary<Tint> boundary{ll_bound};
    return {pfc, boundary};
  } else {
    return {pfc, {}};
  }
}

template<typename T, typename Tint, typename Tgroup>
struct FullComplexEnumeration {
  PerfectComplexTopDimInfo<T,Tint,Tgroup> pctdi;
  std::vector<FacesPerfectComplex<T,Tint,Tgroup>> levels;
  std::vector<FullBoundary<Tint>> boundaries;
};

template<typename T, typename Tint, typename Tgroup>
FullComplexEnumeration<T,Tint,Tgroup> full_perfect_complex_enumeration(std::vector<DatabaseEntry_Serial<PerfectTspace_Obj<T,Tint,Tgroup>,PerfectTspace_AdjO<Tint>>> const& l_tot, LinSpaceMatrix<T> const& LinSpa, bool const& only_well_rounded, bool const& compute_boundary, std::ostream& os) {
  PerfectComplexTopDimInfo<T,Tint,Tgroup> pctdi = generate_perfect_complex_top_dim_info(l_tot, LinSpa, only_well_rounded, compute_boundary, os);
#ifdef PRINT_PERFECT_COMPLEX_DIMENSION
  os << "We have pctdi, |pctdi.l_perfect|=" << pctdi.l_perfect.size() << "\n";
#endif
  int dim_spa = LinSpa.ListMat.size();
  FacesPerfectComplex<T,Tint,Tgroup> level = get_first_step_perfect_complex_enumeration(pctdi, os);
#ifdef PRINT_PERFECT_COMPLEX_DIMENSION
  os << "We have the first level |level|=" << level.l_faces.size() << "\n";
#endif
  std::vector<FacesPerfectComplex<T,Tint,Tgroup>> levels{level};
  std::vector<FullBoundary<Tint>> boundaries;
  for (int i=1; i<dim_spa; i++) {
    ResultStepEnumeration<T,Tint,Tgroup> result = compute_next_level(pctdi, levels[i-1], os);
#ifdef PRINT_PERFECT_COMPLEX_DIMENSION
    os << "i=" << i << " |result.level|=" << result.level.l_faces.size() << "\n";
#endif
    levels.push_back(result.level);
    if (result.boundary) {
      boundaries.push_back(*result.boundary);
    }
  }
  return {pctdi, levels, boundaries};
}

template<typename T>
MyMatrix<T> tot_set(MyMatrix<T> const& EXTin) {
  std::set<MyVector<T>> set;
  int dim = EXTin.cols();
  auto f_insert=[&](int i_row) -> void {
    for (int i=0; i<dim; i++) {
      if (EXTin(i_row, i) > 0) {
        MyVector<T> V = GetMatrixRow(EXTin, i_row);
        set.insert(V);
        return;
      }
    }
  };
  for (int i_row=0; i_row<EXTin.rows(); i_row++) {
    f_insert(i_row);
  }
  MyMatrix<T> M(EXTin.rows(), EXTin.cols());
  int pos = 0;
  for (auto & V: set) {
    AssignMatrixRow(M, pos, V);
    pos += 1;
  }
  return M;
}



template<typename T, typename Tint, typename Tgroup>
struct ChainBuilder {
private:
  std::unordered_map<MyMatrix<Tint>, size_t> map_ext_set;
  int idim;
  FullComplexEnumeration<T,Tint,Tgroup> const& fce;
  std::ostream& os;
public:
  std::vector<FaceEntry<Tint>> chain;
  ChainBuilder(int _idim, FullComplexEnumeration<T,Tint,Tgroup> const& _fce, std::ostream& _os) : idim(_idim), fce(_fce), os(_os) {
  }
  void f_insert(Tint const& value, int const& iOrb, MyMatrix<Tint> const& M) {
    FacePerfectComplex<T,Tint,Tgroup> const& face = fce.levels[idim].l_faces[iOrb];
    MyMatrix<Tint> const& EXT1 = face.EXT;
    OrientationInfo const& or_info = face.or_info;
    std::vector<MyMatrix<T>> const& ListMat = fce.pctdi.LinSpa.ListMat;
    MyMatrix<Tint> EXT2 = EXT1 * M;
    MyMatrix<Tint> EXT3 = tot_set(EXT2);
    size_t& idx = map_ext_set[EXT3];
    if (idx == 0) {
      idx = map_ext_set.size();
      FaceEntry<Tint> fe{iOrb, value, M};
      chain.emplace_back(std::move(fe));
    } else {
      int sign = 1;
      if (face.is_orientable) {
        MyMatrix<Tint> t = chain[idx - 1].M * Inverse(M);;
        sign = get_face_orientation(EXT1, ListMat, or_info, t);
      }
      chain[idx - 1].value += sign * value;
    }
  }
  std::vector<FaceEntry<Tint>> get_faces() const {
    std::vector<FaceEntry<Tint>> chain_ret;
    for (auto & fe: chain) {
      if (fe.value != 0) {
        chain_ret.push_back(fe);
      }
    }
    return chain_ret;
  }
  bool is_zero_chain() const {
    for (auto & fe: chain) {
      if (fe.value != 0) {
        return false;
      }
    }
    return true;
  }
};


template<typename T, typename Tint, typename Tgroup>
bool is_equal_chain(std::vector<FaceEntry<Tint>> const& chain1, std::vector<FaceEntry<Tint>> const& chain2, int const& idim, FullComplexEnumeration<T,Tint,Tgroup> const& fce, std::ostream& os) {
  ChainBuilder<T,Tint,Tgroup> chain_builder(idim, fce, os);
  for (auto & fe: chain1) {
    chain_builder(fe.value, fe.iOrb, fe.M);
  }
  for (auto & fe: chain2) {
    chain_builder(-fe.value, fe.iOrb, fe.M);
  }
  return chain_builder.is_zero_chain();
}




// Compute the boundary d(c) for c a part of the full complex.
template<typename T, typename Tint, typename Tgroup>
std::vector<FaceEntry<Tint>> compute_boundary(std::vector<FaceEntry<Tint>> const& c_idim, int const& idim, FullComplexEnumeration<T,Tint,Tgroup> const& fce, std::ostream& os) {
  ChainBuilder<T,Tint,Tgroup> chain_builder(idim+1, fce, os);
  FullBoundary<Tint> const& bnd = fce.boundaries[idim];
  for (auto & fe: c_idim) {
    int iOrb = fe.iOrb;
    for (auto & ebnd: bnd.ll_bound[iOrb].l_bound) {
      Tint value = fe.value * ebnd.sign;
      MyMatrix<Tint> M = ebnd.M * fe.M;
      chain_builder.f_insert(value, ebnd.iOrb, M);
    }
  }
  return chain_builder.get_faces();
}

// We should have d_{i+1}(d_i(c)) = 0 consistently.
template<typename T, typename Tint, typename Tgroup>
bool is_product_zero(int const& idim, FullComplexEnumeration<T,Tint,Tgroup> const& fce, std::ostream& os) {
  int dim = fce.levels[idim].l_faces.size();
  int n = fce.pctdi.LinSpa.n;
  for (int iOrb=0; iOrb<dim; iOrb++) {
    FaceEntry<Tint> fe{iOrb, Tint(1), IdentityMat<Tint>(n)};
    std::vector<FaceEntry<Tint>> chain1{fe};
    std::vector<FaceEntry<Tint>> chain2 = compute_boundary(chain1, idim, fce, os);
    std::vector<FaceEntry<Tint>> chain3 = compute_boundary(chain2, idim+1, fce, os);
    if (chain3.size() > 0) {
      return false;
    }
  }
  return true;
}

template<typename T, typename Tint, typename Tgroup>
bool are_product_zeros(FullComplexEnumeration<T,Tint,Tgroup> const& fce, std::ostream& os) {
  int dim_spa = fce.pctdi.LinSpa.ListMat.size();
  for (int idim=0; idim<dim_spa-2; idim++) {
    bool test = is_product_zero(idim, fce, os);
    if (!test) {
      return false;
    }
  }
  return true;
}

template<typename T, typename Tint, typename Tgroup>
std::vector<FaceEntry<Tint>> contracting_homotopy_kernel(std::vector<FaceEntry<Tint>> const& c_idim, int const& idim, FullComplexEnumeration<T,Tint,Tgroup> const& fce, std::ostream& os) {

}



template<typename T, typename Tint, typename Tgroup>
std::vector<FaceEntry<Tint>> contracting_homotopy(std::vector<FaceEntry<Tint>> const& c_idim, int const& idim, FullComplexEnumeration<T,Tint,Tgroup> const& fce, std::ostream& os) {
#ifdef SANITY_CHECK_COMPLEX_CONTRACTION
  std::vector<FaceEntry<Tint>> chain3 = compute_boundary(chain2, idim+1, fce, os);
  if (chain3.size() > 0) {
    std::cerr << "COMPCONT: The chain should have a zero boundary\n";
    throw TerminalException{1};
  }
#endif
  std::vector<FaceEntry<Tint>> one_sol = contracting_homotopy_kernel(c_idim, idim, fce, os);
#ifdef SANITY_CHECK_COMPLEX_CONTRACTION
  std::vector<FaceEntry<Tint>> one_sol_img = compute_boundary(one_sol, idim - 1, fce, os);
  if (is_equal_chain(one_sol_img, c_idim, idim, pctdi, rse, os)) {
    std::cerr << "COMPCONT: The proposed preimage is not a solution\n";
    throw TerminalException{1};
  }
#endif
  return one_sol;
}

// clang-format off
#endif  // SRC_PERFECT_PERFECT_COMPLEX_H_
// clang-format on

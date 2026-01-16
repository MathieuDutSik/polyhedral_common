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
 */

//
// The top dimensional complex
//

template<typename T, typename Tint_impl, typename Tgroup_impl>
struct PerfectFormInfoForComplex {
  using Tint = Tint_impl;
  using Tgroup = Tgroup_impl;
  MyMatrix<T> gram;
  MyMatrix<Tint> EXT; // Not necessarily full-dimensional
  MyMatrix<Tint> EXT_fd; // This is full-dimensional or empty.
  Tgroup GRP_ext; // Group acting on the shortest vectors
  std::vector<sing_adj<Tint>> l_sing_adj;
  MyMatrix<Tint> find_matrix(typename Tgroup::Telt const& x) const {
    if (EXT_fd.rows() == 0) {
      return FindTransformation(EXT, EXT, x);
    }
    return FindTransformation(EXT_fd, EXT_fd, x);
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
    // The vector family for the determination.
    MyMatrix<Tint> EXT_fd = get_fulldim_shv_tint<T,Tint>(ePerf.x.Gram, ePerf.x.rec_shv, os);
    if (RankMat(EXT) == EXT.cols()) { // EXT is already full dimensional, no need for EXT_fd
      EXT_fd = MyMatrix<Tint>(0, EXT.cols());
    }
    Tgroup const& GRP_ext = ePerf.x.GRP;
    PerfectFormInfoForComplex<T,Tint,Tgroup> perfect{ePerf.x.Gram, std::move(EXT), std::move(EXT_fd), GRP_ext, std::move(l_sing_adj)};
    l_perfect.emplace_back(std::move(perfect));
  }
  return {std::move(l_perfect), LinSpa, only_well_rounded, compute_boundary};
}

//
// The intermediate dimensional faces.
//

template<typename T, typename Tint, typename Tgroup>
struct FacePerfectComplex {
  std::vector<triple<Tint>> l_triple; // The containing triples.
  std::vector<MyMatrix<Tint>> l_gens; // generating set of the stabilizer
  MyMatrix<Tint> EXT;
  Tgroup GRP_ext; // Group acting on the vectors.
  bool is_well_rounded;
};


template<typename T, typename Tint, typename Tgroup>
struct FacesPerfectComplex {
  std::vector<FacePerfectComplex<T,Tint,Tgroup>> l_faces;
};

template<typename T, typename Tint, typename Tgroup>
FacesPerfectComplex<T,Tint,Tgroup> get_first_step_perfect_complex_enumeration(PerfectComplexTopDimInfo<T,Tint,Tgroup> const& pctdi, [[maybe_unused]] std::ostream & os) {
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
      MyMatrix<Tint> eMatrGen = top.find_matrix(ePermGen);
      l_gens.push_back(eMatrGen);
    }
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: Returning l_gens\n";
#endif
    return l_gens;
  };
  for (size_t i_domain=0; i_domain<pctdi.l_perfect.size(); i_domain++) {
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: get_first_step i_domain=" << i_domain << "\n";
#endif
    std::vector<triple<Tint>> l_triple = get_l_triple(i_domain);
    std::vector<MyMatrix<Tint>> l_gens = get_l_gens(i_domain);
    MyMatrix<Tint> const& EXT = pctdi.l_perfect[i_domain].EXT;
    Tgroup const& GRP_ext = pctdi.l_perfect[i_domain].GRP_ext;
    bool is_well_rounded = true; // Yes, the top dimensional cells are well rounded.
    FacePerfectComplex<T,Tint,Tgroup> face{l_triple, l_gens, EXT, GRP_ext, is_well_rounded};
    l_faces.push_back(face);
    i_domain += 1;
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

template<typename T, typename Tint, typename Tgroup>
ResultStepEnumeration<T,Tint,Tgroup> compute_next_level(PerfectComplexTopDimInfo<T,Tint,Tgroup> const& pctdi, FacesPerfectComplex<T,Tint,Tgroup> const& level, std::ostream & os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n = pctdi.LinSpa.n;
#ifdef DEBUG_PERFECT_COMPLEX
  os << "PERFCOMP: compute_next_level, start, n=" << n << "\n";
#endif
  std::vector<FacePerfectComplex<T,Tint,Tgroup>> l_faces;
  std::vector<ListBoundEntry<Tint>> ll_bound;
  auto find_matching_entry=[&](triple<Tint> const& t) -> std::optional<BoundEntry<Tint>> {
    int i_domain = 0;
    for (auto & face1: l_faces) {
      std::optional<MyMatrix<Tint>> opt =
        test_triple_in_listtriple(pctdi.l_perfect, face1.l_triple, t);
      if (opt) {
        MyMatrix<Tint> const& M = *opt;
        int sign = 1;
        BoundEntry<Tint> be{i_domain, sign, M};
        return be;
      }
      i_domain += 1;
    }
    return {};
  };
  auto need_opt_t=[&]([[maybe_unused]] PerfectBoundednessProperty const& pbp) -> bool {
#ifdef SANITY_CHECK_PERFECT_COMPLEX
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
  using Tfull_triple = std::pair<std::vector<triple<Tint>>, std::vector<MyMatrix<Tint>>>;
  auto f_is_well_rounded=[&](triple<Tint> const& t, std::optional<Tfull_triple> & opt_t, PerfectBoundednessProperty & pbp) -> bool {
    if (need_opt_t(pbp)) {
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
    return t;
  };
  auto f_insert=[&](triple<Tint> const& t, std::optional<Tfull_triple> & opt_t, bool const& is_well_rounded) -> BoundEntry<Tint> {
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: f_insert, step 1\n";
#endif
    std::optional<BoundEntry<Tint>> opt = find_matching_entry(t);
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
    ContainerMatrix<Tint> cont(EXT);
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: f_insert, step 6\n";
#endif
    std::vector<Telt> l_gens;
    for (auto & eMatrGen: pair.second) {
      std::vector<Tidx> v(n_ext);
      for (int u=0; u<n_ext; u++) {
        MyVector<Tint> V1 = GetMatrixRow(EXT, u);
        MyVector<Tint> V2 = eMatrGen.transpose() * V1;
        std::optional<size_t> opt = cont.GetIdx_v(V2);
        size_t pos = unfold_opt(opt, "PERFCOMPL: Failed to find the index");
        v[u] = pos;
      }
      Telt elt(v);
      l_gens.push_back(elt);
    }
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: f_insert, step 7\n";
#endif
    Tgroup GRP_ext(l_gens, n_ext);
    FacePerfectComplex<T,Tint,Tgroup> face{
      pair.first,
      pair.second,
      EXT,
      GRP_ext,
      is_well_rounded};
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: f_insert, step 8\n";
#endif
    int i_domain = l_faces.size();
    l_faces.push_back(face);
    MyMatrix<Tint> M = IdentityMat<Tint>(n);
    int sign = 1;
    BoundEntry<Tint> be{i_domain, sign, M};
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: f_insert, step 9\n";
#endif
    return be;
  };
  for (size_t i=0; i<level.l_faces.size(); i++) {
    FacePerfectComplex<T,Tint,Tgroup> const& face = level.l_faces[i];
    std::vector<BoundEntry<Tint>> l_bound;
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T,Tint>(face.EXT);
    RyshkovGRP<T, Tgroup> eCone =
      GetNakedPerfectCone_GRP<T, Tgroup>(pctdi.LinSpa, EXT_T, face.GRP_ext, os);
    for (auto& incd_sma: eCone.ListIncd) {
      Face incd_big = get_big_incd(eCone, incd_sma);
      MyMatrix<Tint> EXTincd = SelectRow(face.EXT, incd_big);
      PerfectBoundednessProperty pbp = initial_bounded_property(pctdi.LinSpa, EXTincd, os);
      std::optional<Tfull_triple> opt_t;
      triple<Tint> t = get_initial_triple(face, incd_big);
      bool is_well_rounded = f_is_well_rounded(t, opt_t, pbp);
#ifdef DEBUG_PERFECT_COMPLEX
      os << "PERFCOMP: pbp=" << to_string(pbp) << "\n";
#endif
      if (is_insertable(is_well_rounded)) {
        BoundEntry<Tint> be = f_insert(t, opt_t, is_well_rounded);
        if (pctdi.compute_boundary) {
          l_bound.push_back(be);
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

// clang-format off
#endif  // SRC_PERFECT_PERFECT_COMPLEX_H_
// clang-format on

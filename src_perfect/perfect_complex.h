// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_PERFECT_COMPLEX_H_
#define SRC_PERFECT_PERFECT_COMPLEX_H_

// clang-format off>
#include "perfect_tspace.h"
#include "triples.h"
#include "GampMatlab.h"
#include "boost_serialization.h"
#include <boost/dynamic_bitset/serialization.hpp>
#include <fstream>
// clang-format on

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

  Computation of the contracting homotopy.
  
  
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
        std::cerr << "PERFCOMP: Inconsistency in the pre_image computation\n";
        throw TerminalException{1};
      }
#endif
      return M;
    }
    return FindTransformation(EXT, EXT, x);
  }
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void save(Archive &ar, PerfectFormInfoForComplex<T, Tint, Tgroup> const &val,
                 [[maybe_unused]] const unsigned int version) {
  bool has_pre_imager = val.opt_pre_imager.has_value();
  ar &make_nvp("gram", val.gram);
  ar &make_nvp("EXT", val.EXT);
  ar &make_nvp("GRP_ext", val.GRP_ext);
  ar &make_nvp("l_sing_adj", val.l_sing_adj);
  ar &make_nvp("has_pre_imager", has_pre_imager);
  if (has_pre_imager) {
    std::vector<MyMatrix<Tint>> l_matr =
        val.opt_pre_imager->get_list_matr_gens();
    ar &make_nvp("l_matr", l_matr);
  }
}

template <class Archive, typename T, typename Tint, typename Tgroup>
inline void load(Archive &ar, PerfectFormInfoForComplex<T, Tint, Tgroup> &val,
                 [[maybe_unused]] const unsigned int version) {
  bool has_pre_imager = false;
  ar &make_nvp("gram", val.gram);
  ar &make_nvp("EXT", val.EXT);
  ar &make_nvp("GRP_ext", val.GRP_ext);
  ar &make_nvp("l_sing_adj", val.l_sing_adj);
  ar &make_nvp("has_pre_imager", has_pre_imager);
  if (has_pre_imager) {
    std::vector<MyMatrix<Tint>> l_matr;
    ar &make_nvp("l_matr", l_matr);
    if (l_matr.empty()) {
      std::cerr << "PERFCOMP: Missing l_matr for preimage reconstruction\n";
      throw TerminalException{1};
    }
    using Telt = typename Tgroup::Telt;
    using TintGroup = typename Tgroup::Tint;
    std::vector<Telt> l_perm =
        get_list_elt_from_list_matrices<Tint, Telt>(l_matr, val.EXT, std::cerr);
    MyMatrix<Tint> id = IdentityMat<Tint>(val.EXT.cols());
    val.opt_pre_imager = permutalib::PreImagerElement<Telt, MyMatrix<Tint>, TintGroup>(
        l_matr, l_perm, id);
  } else {
    val.opt_pre_imager = {};
  }
}

template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, PerfectFormInfoForComplex<T, Tint, Tgroup> &val,
                      const unsigned int version) {
  split_free(ar, val, version);
}
}  // namespace boost::serialization

struct PerfectComplexOptions {
  bool only_well_rounded;
  bool compute_boundary;
  bool compute_contracting_homotopy;
};

namespace boost::serialization {
template <class Archive>
inline void serialize(Archive &ar, PerfectComplexOptions &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("only_well_rounded", val.only_well_rounded);
  ar &make_nvp("compute_boundary", val.compute_boundary);
  ar &make_nvp("compute_contracting_homotopy", val.compute_contracting_homotopy);
}
}  // namespace boost::serialization


template<typename T, typename Tint, typename Tgroup>
struct PerfectComplexTopDimInfo {
  std::vector<PerfectFormInfoForComplex<T,Tint,Tgroup>> l_perfect;
  LinSpaceMatrix<T> LinSpa;
  PerfectComplexOptions pco;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, PerfectComplexTopDimInfo<T, Tint, Tgroup> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("l_perfect", val.l_perfect);
  ar &make_nvp("LinSpa", val.LinSpa);
  ar &make_nvp("pco", val.pco);
}
}  // namespace boost::serialization



template<typename T, typename Tint, typename Tgroup>
PerfectComplexTopDimInfo<T,Tint,Tgroup> generate_perfect_complex_top_dim_info(std::vector<DatabaseEntry_Serial<PerfectTspace_Obj<T,Tint,Tgroup>,PerfectTspace_AdjO<Tint>>> const& l_tot, LinSpaceMatrix<T> const& LinSpa, PerfectComplexOptions const& pco, std::ostream& os) {
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
  return {std::move(l_perfect), LinSpa, pco};
}

struct OrientationInfo {
  int sign;
  std::vector<int> ListRowSelect;
  std::vector<int> ListColSelect;
};

namespace boost::serialization {
template <class Archive>
inline void serialize(Archive &ar, OrientationInfo &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("sign", val.sign);
  ar &make_nvp("ListRowSelect", val.ListRowSelect);
  ar &make_nvp("ListColSelect", val.ListColSelect);
}
}  // namespace boost::serialization


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

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, FacePerfectComplex<T, Tint, Tgroup> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("l_triple", val.l_triple);
  ar &make_nvp("l_gens", val.l_gens);
  ar &make_nvp("EXT", val.EXT);
  ar &make_nvp("GRP_ext", val.GRP_ext);
  ar &make_nvp("is_well_rounded", val.is_well_rounded);
  ar &make_nvp("or_info", val.or_info);
  ar &make_nvp("is_orientable", val.is_orientable);
}
}  // namespace boost::serialization



template<typename T, typename Tint, typename Tgroup>
struct FacesPerfectComplex {
  std::vector<FacePerfectComplex<T,Tint,Tgroup>> l_faces;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, FacesPerfectComplex<T, Tint, Tgroup> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("l_faces", val.l_faces);
}
}  // namespace boost::serialization



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
#ifdef SANITY_CHECK_PERFECT_COMPLEX
  if (det == 0) {
    std::cerr << "PERFCOMP: determinant is 0, big error A\n";
    throw TerminalException{1};
  }
#endif
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
#ifdef SANITY_CHECK_PERFECT_COMPLEX
  if (det == 0) {
    std::cerr << "PERFCOMP: determinant is 0, big error B\n";
    throw TerminalException{1};
  }
#endif
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

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, BoundEntry<Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("iOrb", val.iOrb);
  ar &make_nvp("sign", val.sign);
  ar &make_nvp("M", val.M);
}
}  // namespace boost::serialization


template<typename Tint>
struct ListBoundEntry {
  std::vector<BoundEntry<Tint>> l_bound;
};

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, ListBoundEntry<Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("l_bound", val.l_bound);
}
}  // namespace boost::serialization


template<typename Tint>
struct FullBoundary {
  std::vector<ListBoundEntry<Tint>> ll_bound;
};

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, FullBoundary<Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("ll_bound", val.ll_bound);
}
}  // namespace boost::serialization



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

template<typename T, typename Tint>
struct PerfectFaceEntry {
  int iOrb;
  T value;
  MyMatrix<Tint> M;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint>
inline void serialize(Archive &ar, PerfectFaceEntry<T, Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("iOrb", val.iOrb);
  ar &make_nvp("value", val.value);
  ar &make_nvp("M", val.M);
}
}  // namespace boost::serialization


template<typename T, typename Tint, typename Tgroup>
ResultStepEnumeration<T,Tint,Tgroup> compute_next_level(PerfectComplexTopDimInfo<T,Tint,Tgroup> const& pctdi, FacesPerfectComplex<T,Tint,Tgroup> const& level, std::ostream & os) {
  using Telt = typename Tgroup::Telt;
  int n = pctdi.LinSpa.n;
#ifdef DEBUG_PERFECT_COMPLEX
  os << "PERFCOMP: compute_next_level, start, n=" << n << "\n";
  int i_level = 0;
  if (level.l_faces.size() > 0) {
    int dim_spa = pctdi.LinSpa.ListMat.size();
    int dim_ext = level.l_faces[0].or_info.ListRowSelect.size();
    i_level = dim_spa - dim_ext;
  }
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
    std::cerr << "PERFCOMP: Returning true because we want to compute whether the group is finite\n";
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
    // That check is really cheap. Do it all the time.
    check_pbp(pbp);
    return get_result(pbp);
  };
  auto is_insertable=[&](bool is_well_rounded) -> bool {
    if (pctdi.pco.only_well_rounded) {
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
  std::vector<int> l_gens_sign;
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
#ifdef SANITY_CHECK_PERFECT_COMPLEX
    if (det == 0) {
      std::cerr << "PERFCOMP: determinant is 0, big error C\n";
      throw TerminalException{1};
    }
#endif
    int sign = face.or_info.sign;
    if (det > 0) {
      return sign;
    } else {
      return -sign;
    }
  };
  auto initial_boundary_setup=[&](size_t const& i, RyshkovGRP<T, Tgroup> const& cone) -> void {
    FacePerfectComplex<T,Tint,Tgroup> const& face = level.l_faces[i];
    interior_pt = face.get_interior_point(pctdi.LinSpa.ListMat);
    l_gens_perm.clear();
    l_gens_sign.clear();
    ElementMapper<Tint,Telt> elt_mapper(face.EXT);
    for (auto & eMatrGen: face.l_gens) {
      Telt elt1 = elt_mapper.map_elt(eMatrGen);
      Telt elt2 = cone.map_elt(elt1);
      int sign = get_face_orientation(face.EXT, pctdi.LinSpa.ListMat, face.or_info, eMatrGen);
      l_gens_perm.push_back(elt2);
      l_gens_sign.push_back(sign);
    }
#ifdef SANITY_CHECK_PERFECT_COMPLEX
    set_EXT.clear();
    for (int i_row=0; i_row<face.EXT.rows(); i_row++) {
      MyVector<Tint> V = GetMatrixRow(face.EXT, i_row);
      set_EXT.insert(V);
    }
#endif
  };
  auto append_boundary=[&](size_t const& i, std::pair<int,MyMatrix<Tint>> const& p,
                           RyshkovGRP<T, Tgroup> const& cone,
                           Face const& incd_sma, std::vector<BoundEntry<Tint>> & l_bound) -> void {
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
        BoundEntry<Tint> be{p.first, new_sign, new_mat};
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
    size_t start = 0;
    while(true) {
      size_t len = l_faces_gen.size();
      for (size_t u=start; u<len; u++) {
        for (size_t i_gen=0; i_gen<n_gen; i_gen++) {
          Face new_incd_sma = OnFace(l_faces_gen[u], l_gens_perm[i_gen]);
          int new_sign = l_sign[u] * l_gens_sign[i_gen];
          MyMatrix<Tint> new_mat = l_mat[u] * l_gens[i_gen];
          insert_entry(new_incd_sma, new_sign, new_mat);
        }
      }
      start = len;
      if (start == l_faces_gen.size()) {
        break;
      }
    }
#ifdef SANITY_CHECK_PERFECT_COMPLEX
    using TintGroup = typename Tgroup::Tint;
    Tgroup stab = cone.GRPsub.Stabilizer_OnSets(incd_sma);
    TintGroup size1 = cone.GRPsub.size() / stab.size();
    TintGroup size2 = l_faces_gen.size();
    if (size1 != size2) {
      std::cerr << "PERFCOMP: Inconsistent orbit sizes size1=" << size1 << " size2=" << size2 << "\n";
      throw TerminalException{1};
    }
#endif
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: append_boundary l_sign=";
    for (auto &esign: l_sign) {
      os << " " << esign;
    }
    os << "\n";
#endif
  };
  for (size_t i=0; i<level.l_faces.size(); i++) {
    FacePerfectComplex<T,Tint,Tgroup> const& face = level.l_faces[i];
    std::vector<BoundEntry<Tint>> l_bound;
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T,Tint>(face.EXT);
    RyshkovGRP<T, Tgroup> cone =
      GetNakedPerfectCone_GRP<T, Tgroup>(pctdi.LinSpa, EXT_T, face.GRP_ext, os);
    if (pctdi.pco.compute_boundary) {
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
        if (pctdi.pco.compute_boundary) {
          append_boundary(i, p, cone, incd_sma, l_bound);
        }
      }
    }
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: i_level=" << i_level << " i=" << i
       << " |cone|=" << cone.PerfDomEXT.rows()
       << " |GRPsub|=" << cone.GRPsub.size()
       << " |ListIncd|=" << cone.ListIncd.size()
       << " |l_bound|=" << l_bound.size() << "\n";
    os << "PERFCOMP: i_level=" << i_level << " i=" << i << " status=";
    for (auto & e_sign: l_gens_sign) {
      os << " " << e_sign;
    }
    os << "\n";
#endif
#ifdef SANITY_CHECK_PERFECT_COMPLEX
    if (pctdi.pco.compute_boundary && !pctdi.pco.only_well_rounded) {
      std::vector<Telt> gens;
      Tgroup GRPtriv(gens, cone.PerfDomEXT.rows());
      vectface ListIncd_triv =
        DualDescriptionStandard<T, Tgroup>(cone.PerfDomEXT, GRPtriv, os);
      size_t size1 = ListIncd_triv.size();
      size_t size2 = l_bound.size();
      if (size1 != size2) {
        std::cerr << "PERFCOMP: inconsistent size size1=" << size1 << " size2=" << size2 << "\n";
        throw TerminalException{1};
      }
    }
#endif
    if (pctdi.pco.compute_boundary) {
      ListBoundEntry<Tint> lbe{l_bound};
      ll_bound.push_back(lbe);
    }
  }
  FacesPerfectComplex<T,Tint,Tgroup> pfc{l_faces};
  if (pctdi.pco.compute_boundary) {
    FullBoundary<Tint> boundary{ll_bound};
    return {pfc, boundary};
  } else {
    return {pfc, {}};
  }
}

template<typename Tint>
struct PerfectFace {
  int iOrb;
  MyMatrix<Tint> M;
};

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, PerfectFace<Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("iOrb", val.iOrb);
  ar &make_nvp("M", val.M);
}
}  // namespace boost::serialization


template<typename Tint>
struct TopDimensional {
  std::vector<std::vector<PerfectFace<Tint>>> ll_faces;
};

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, TopDimensional<Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("ll_faces", val.ll_faces);
}
}  // namespace boost::serialization



template<typename T, typename Tint, typename Tgroup>
struct FullComplexEnumeration {
  PerfectComplexTopDimInfo<T,Tint,Tgroup> pctdi;
  std::vector<FacesPerfectComplex<T,Tint,Tgroup>> levels;
  std::vector<FullBoundary<Tint>> boundaries;
  std::vector<TopDimensional<Tint>> l_topdims;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, FullComplexEnumeration<T, Tint, Tgroup> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("pctdi", val.pctdi);
  ar &make_nvp("levels", val.levels);
  ar &make_nvp("boundaries", val.boundaries);
  ar &make_nvp("l_topdims", val.l_topdims);
}
}  // namespace boost::serialization

template<typename T, typename Tint, typename Tgroup>
void WriteFullComplexEnumeration(std::string const &file,
                                 FullComplexEnumeration<T, Tint, Tgroup> const &fce) {
  std::ofstream ofs(file);
  boost::archive::text_oarchive oa(ofs);
  oa << fce;
}

template<typename T, typename Tint, typename Tgroup>
FullComplexEnumeration<T, Tint, Tgroup> ReadFullComplexEnumeration(
    std::string const &file) {
  FullComplexEnumeration<T, Tint, Tgroup> fce;
  std::ifstream ifs(file);
  boost::archive::text_iarchive ia(ifs);
  ia >> fce;
  return fce;
}


template<typename T, typename Tint, typename Tgroup>
FullComplexEnumeration<T,Tint,Tgroup> full_perfect_complex_enumeration(std::vector<DatabaseEntry_Serial<PerfectTspace_Obj<T,Tint,Tgroup>,PerfectTspace_AdjO<Tint>>> const& l_tot, LinSpaceMatrix<T> const& LinSpa, PerfectComplexOptions const& pco, std::ostream& os) {
  using Telt = typename Tgroup::Telt;
  PerfectComplexTopDimInfo<T,Tint,Tgroup> pctdi = generate_perfect_complex_top_dim_info(l_tot, LinSpa, pco, os);
#ifdef PRINT_PERFECT_COMPLEX_DIMENSION
  os << "PERFCOMP: We have pctdi, |pctdi.l_perfect|=" << pctdi.l_perfect.size() << "\n";
#endif
  int dim_spa = LinSpa.ListMat.size();
  FacesPerfectComplex<T,Tint,Tgroup> level = get_first_step_perfect_complex_enumeration(pctdi, os);
#ifdef PRINT_PERFECT_COMPLEX_DIMENSION
  os << "PERFCOMP: We have the first level |level|=" << level.l_faces.size() << "\n";
#endif
  std::vector<FacesPerfectComplex<T,Tint,Tgroup>> levels{level};
  std::vector<FullBoundary<Tint>> boundaries;
  for (int i=1; i<dim_spa; i++) {
    ResultStepEnumeration<T,Tint,Tgroup> result = compute_next_level(pctdi, levels[i-1], os);
#ifdef PRINT_PERFECT_COMPLEX_DIMENSION
    int n_orientable = 0;
    for (auto & ent: result.level.l_faces) {
      if (ent.is_orientable) {
        n_orientable += 1;
      }
    }
    os << "PERFCOMP: i=" << i << " |result.level|=" << result.level.l_faces.size() << " n_orientable=" << n_orientable << "\n";
#endif
    levels.push_back(result.level);
    if (result.boundary) {
      boundaries.push_back(*result.boundary);
    }
  }
  std::vector<TopDimensional<Tint>> l_topdims;
  if (pco.compute_contracting_homotopy) {
    size_t n_perfect = pctdi.l_perfect.size();
    size_t dim = levels.size();
    std::vector<PerfectFace<Tint>> l_faces;
    TopDimensional<Tint> top{std::vector<std::vector<PerfectFace<Tint>>>(dim, l_faces)};
    l_topdims = std::vector<TopDimensional<Tint>>(n_perfect, top);
    struct PairPermMatr {
      Telt ePerm;
      MyMatrix<Tint> eMatr;
    };
    std::vector<std::vector<PairPermMatr>> ll_pair;
    for (auto & ePerfect: pctdi.l_perfect) {
      std::vector<PairPermMatr> l_pair;
      for (auto & ePerm: ePerfect.GRP_ext.GeneratorsOfGroup()) {
        MyMatrix<Tint> eMatr = ePerfect.find_matrix(ePerm, os);
        PairPermMatr pair{ePerm, eMatr};
        l_pair.push_back(pair);
      }
      ll_pair.push_back(l_pair);
    }
    for (size_t i_dim=0; i_dim<dim; i_dim++) {
      FacesPerfectComplex<T,Tint,Tgroup> const& level = levels[i_dim];
      int n_face = level.l_faces.size();
      for (int i_face=0; i_face<n_face; i_face++) {
        FacePerfectComplex<T,Tint,Tgroup> const& face = level.l_faces[i_face];
        for (auto & triple: face.l_triple) {
          size_t i_perfect = triple.iCone;
          std::unordered_set<Face> set;
          std::vector<Face> l_set;
          std::vector<PerfectFace<Tint>>& l_pf = l_topdims[i_perfect].ll_faces[i_dim];
          size_t n_exist = l_pf.size();
          auto g_insert=[&](Face const& f, MyMatrix<Tint> const& M) -> void {
            if (set.count(f) == 0) {
              set.insert(f);
              l_set.push_back(f);
              PerfectFace<Tint> pf{i_face, M};
              l_pf.push_back(pf);
            }
          };
          g_insert(triple.f_ext, triple.eMat);
          size_t start = 0;
          while(true) {
            size_t len = l_pf.size() - n_exist;
            for (size_t u=start; u<len; u++) {
              for (auto & pair: ll_pair[i_perfect]) {
                Face set_img = OnFace(l_set[u], pair.ePerm);
#ifdef DEBUG_PERFECT_COMPLEX
                os << "PERFCOMP: n_exist=" << n_exist << " |l_pf|=" << l_pf.size() << "\n";
                os << "PERFCOMP: l_pf[n_exist + u].M=\n";
                WriteMatrix(os, l_pf[n_exist + u].M);
                os << "PERFCOMP: pair.eMatr=\n";
                WriteMatrix(os, pair.eMatr);
#endif
                MyMatrix<Tint> M_img = l_pf[n_exist + u].M * pair.eMatr;
                g_insert(set_img, M_img);
              }
            }
            start = len;
            if (n_exist + start == l_pf.size()) {
              break;
            }
          }
        }
      }
    }
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: full_perfect_complex_enumeration, we have l_topdims\n";
#endif
  }
  return {pctdi, levels, boundaries, l_topdims};
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
  std::vector<PerfectFaceEntry<T, Tint>> chain;
  ChainBuilder(int _idim, FullComplexEnumeration<T,Tint,Tgroup> const& _fce, std::ostream& _os) : idim(_idim), fce(_fce), os(_os) {
  }
  void f_insert(T const& value, int const& iOrb, MyMatrix<Tint> const& M) {
    FacePerfectComplex<T,Tint,Tgroup> const& face = fce.levels[idim].l_faces[iOrb];
    MyMatrix<Tint> const& EXT1 = face.EXT;
    OrientationInfo const& or_info = face.or_info;
    std::vector<MyMatrix<T>> const& ListMat = fce.pctdi.LinSpa.ListMat;
    MyMatrix<Tint> EXT2 = EXT1 * M;
    MyMatrix<Tint> EXT3 = tot_set(EXT2);
    size_t& idx = map_ext_set[EXT3];
    if (idx == 0) {
      idx = map_ext_set.size();
      PerfectFaceEntry<T, Tint> fe{iOrb, value, M};
      chain.emplace_back(std::move(fe));
    } else {
      MyMatrix<Tint> t = chain[idx - 1].M * Inverse(M);;
      int sign = get_face_orientation(EXT1, ListMat, or_info, t);
      chain[idx - 1].value += sign * value;
    }
  }
  std::vector<PerfectFaceEntry<T, Tint>> get_faces() const {
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: ChainBuilder get_faces |chain|=" << chain.size() << "\n";
#endif
    std::vector<PerfectFaceEntry<T, Tint>> chain_ret;
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
bool is_equal_chain(std::vector<PerfectFaceEntry<T, Tint>> const& chain1, std::vector<PerfectFaceEntry<T, Tint>> const& chain2, int const& idim, FullComplexEnumeration<T,Tint,Tgroup> const& fce, std::ostream& os) {
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
std::vector<PerfectFaceEntry<T, Tint>> compute_boundary(std::vector<PerfectFaceEntry<T, Tint>> const& chain, int const& idim, FullComplexEnumeration<T,Tint,Tgroup> const& fce, std::ostream& os) {
  ChainBuilder<T,Tint,Tgroup> chain_builder(idim+1, fce, os);
  FullBoundary<Tint> const& bnd = fce.boundaries[idim];
  for (auto & fe: chain) {
    int iOrb = fe.iOrb;
    for (auto & ebnd: bnd.ll_bound[iOrb].l_bound) {
      T value = fe.value * ebnd.sign;
      MyMatrix<Tint> M = ebnd.M * fe.M;
#ifdef DEBUG_PERFECT_COMPLEX
      os << "PERFCOMP: compute_boundary fe.value=" << fe.value << " sign=" << ebnd.sign << "\n";
#endif
      chain_builder.f_insert(value, ebnd.iOrb, M);
    }
  }
  return chain_builder.get_faces();
}

// We should have d_{i+1}(d_i(c)) = 0 consistently.
template<typename T, typename Tint, typename Tgroup>
bool is_product_zero(int const& idim, FullComplexEnumeration<T,Tint,Tgroup> const& fce, std::ostream& os) {
  int nOrb = fce.levels[idim].l_faces.size();
  int n = fce.pctdi.LinSpa.n;
#ifdef DEBUG_PERFECT_COMPLEX
  {
    size_t tot_dim = fce.boundaries.size();
    os << "PERFCOMP: tot_dim=" << tot_dim << "\n";
    for (size_t i_dim=0; i_dim<tot_dim; i_dim++) {
      FullBoundary<Tint> const& bnd = fce.boundaries[i_dim];
      size_t n_orb = bnd.ll_bound.size();
      os << "PERFCOMP: i_dim=" << i_dim << " n_orb=" << n_orb << "\n";
      for (size_t i_orb=0; i_orb<n_orb; i_orb++) {
        os << "PERFCOMP:   i_orb=" << i_orb << " status=";
        for (auto & ebnd: bnd.ll_bound[i_orb].l_bound) {
          os << " " << ebnd.sign;
        }
        os << "\n";
      }
    }
  }
#endif
  for (int iOrb=0; iOrb<nOrb; iOrb++) {
    PerfectFaceEntry<T, Tint> fe{iOrb, T(1), IdentityMat<Tint>(n)};
    std::vector<PerfectFaceEntry<T, Tint>> chain1{fe};
    std::vector<PerfectFaceEntry<T, Tint>> chain2 = compute_boundary(chain1, idim, fce, os);
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: We have chain2\n";
#endif
    std::vector<PerfectFaceEntry<T, Tint>> chain3 = compute_boundary(chain2, idim+1, fce, os);
#ifdef DEBUG_PERFECT_COMPLEX
    os << "PERFCOMP: We have chain3\n";
#endif
    if (chain3.size() > 0) {
#ifdef DEBUG_PERFECT_COMPLEX
      os << "PERFCOMP: |chain1|=" << chain1.size() << "\n";
      os << "PERFCOMP: |chain2|=" << chain2.size() << "\n";
      os << "PERFCOMP: |chain3|=" << chain3.size() << " values=";
      for (auto& fe: chain3) {
        os << " " << fe.value;
      }
      os << "\n";
      os << "PERFCOMP: not zero at i=" << idim << " iOrb=" << iOrb << "\n";
#endif
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
struct ComplexBuilder {
private:
  std::unordered_map<MyMatrix<Tint>, size_t> map_ext_set;
  int idim;
  FullComplexEnumeration<T,Tint,Tgroup> const& fce;
  std::ostream& os;
public:
  std::vector<PerfectFace<Tint>> faces;
  ComplexBuilder(int _idim, FullComplexEnumeration<T,Tint,Tgroup> const& _fce, std::ostream& _os) : idim(_idim), fce(_fce), os(_os) {
  }
  void f_insert(int const& iOrb, MyMatrix<Tint> const& M) {
    FacePerfectComplex<T,Tint,Tgroup> const& face = fce.levels[idim].l_faces[iOrb];
    MyMatrix<Tint> const& EXT1 = face.EXT;
    MyMatrix<Tint> EXT2 = EXT1 * M;
    MyMatrix<Tint> EXT3 = tot_set(EXT2);
    size_t& idx = map_ext_set[EXT3];
    if (idx == 0) {
      idx = map_ext_set.size();
      PerfectFace<Tint> pf{iOrb, M};
      faces.emplace_back(std::move(pf));
    }
  }
  int dimension() const {
    return faces.size();
  }
  std::vector<PerfectFace<Tint>> const& get_faces() const {
    return faces;
  }
  std::pair<size_t, int> get_position(int const& iOrb, MyMatrix<Tint> const& M) const {
    FacePerfectComplex<T,Tint,Tgroup> const& face = fce.levels[idim].l_faces[iOrb];
    MyMatrix<Tint> const& EXT1 = face.EXT;
    OrientationInfo const& or_info = face.or_info;
    std::vector<MyMatrix<T>> const& ListMat = fce.pctdi.LinSpa.ListMat;
    MyMatrix<Tint> EXT2 = EXT1 * M;
    MyMatrix<Tint> EXT3 = tot_set(EXT2);
    size_t idx = map_ext_set.at(EXT3);
    MyMatrix<Tint> t = faces[idx - 1].M * Inverse(M);;
    int sign = get_face_orientation(EXT1, ListMat, or_info, t);
    size_t pos = idx - 1;
    return {pos, sign};
  }
};



template<typename Tint>
struct TopPerfectCone {
  int i_perfect;
  MyMatrix<Tint> M;
};


/*
  chain is of index idim and we are looking for a chain of index idim-1 which induces it.
  By idim we mean of course, the index in the list, not the dimension
  ----
  Here we try to solve the equation in a specific set of top dimensional cells.
 */
template<typename T, typename Tint, typename Tgroup>
std::optional<std::vector<PerfectFaceEntry<T, Tint>>> contracting_homotopy_specified(std::vector<PerfectFaceEntry<T, Tint>> const& chain, int const& idim, FullComplexEnumeration<T,Tint,Tgroup> const& fce, std::vector<TopPerfectCone<Tint>> const& l_top, std::ostream& os) {
  ComplexBuilder<T,Tint,Tgroup> cb1(idim-1, fce, os);
  ComplexBuilder<T,Tint,Tgroup> cb2(idim, fce, os);
  for (auto &top: l_top) {
    int i_perfect = top.i_perfect;
    for (auto & face1: fce.l_topdims[i_perfect].ll_faces[idim-1]) {
      int iOrb = face1.iOrb;
      MyMatrix<Tint> M = face1.M * top.M;
      cb1.f_insert(iOrb, M);
    }
    for (auto & face2: fce.l_topdims[i_perfect].ll_faces[idim]) {
      int iOrb = face2.iOrb;
      MyMatrix<Tint> M = face2.M * top.M;
      cb2.f_insert(iOrb, M);
    }
  }
  int dim1 = cb1.dimension();
  int dim2 = cb2.dimension();
  // Building the B vector which is solution of XA = B
  // We have nbRow (A) = dim1 and nbCol(A) = dim2;
  int nbRow = dim1;
  int nbCol = dim2;
  MyVector<T> Bvect = ZeroVector<T>(dim2);
  for (auto & fe: chain) {
    std::pair<size_t, int> p = cb2.get_position(fe.iOrb, fe.M);
    Bvect(p.first) = p.second * fe.value; // No duplication allowed in chain.
  }
  // Now building the differential operator
  std::vector<PerfectFace<Tint>> const& faces1 = cb1.get_faces();
  using T2 = Eigen::Triplet<T>;
  std::vector<T2> tripletList;
  for (int i1=0; i1<dim1; i1++) {
    PerfectFace<Tint> const& face1 = faces1[i1];
    int iRow = i1;
    for (auto & ebnd: fce.boundaries[idim-1].ll_bound[face1.iOrb].l_bound) {
      MyMatrix<T> M = ebnd.M * face1.M;
      std::pair<size_t, int> p = cb2.get_position(ebnd.iOrb, M);
      T value = p.second * ebnd.value;
      int iCol = p.first;
      T2 entry(iRow, iCol, value);
      tripletList.push_back(entry);
    }
  }
  MySparseMatrix<T> SpMat(nbRow, nbCol);
  SpMat.setFromTriplets(tripletList.begin(), tripletList.end());
  MyVector<T> eSol1 = AMP_SolutionSparseSystem(SpMat, Bvect, os);

  std::vector<PerfectFaceEntry<T, Tint>> chain_ret;
  for (int i1=0; i1<dim1; i1++) {
    T val = eSol1(i1);
    if (val != 0) {
      PerfectFace<Tint> const& face1 = faces1[i1];
      PerfectFaceEntry<T,Tint> ent{face1.iOrb, val, face1.M};
      chain_ret.push_back(ent);
    }
  }
  return chain_ret;
}






template<typename T, typename Tint, typename Tgroup>
std::vector<PerfectFaceEntry<T, Tint>> contracting_homotopy_kernel(std::vector<PerfectFaceEntry<T, Tint>> const& chain, int const& idim, FullComplexEnumeration<T,Tint,Tgroup> const& fce, std::ostream& os) {

}



template<typename T, typename Tint, typename Tgroup>
std::vector<PerfectFaceEntry<T, Tint>> contracting_homotopy(std::vector<PerfectFaceEntry<T, Tint>> const& chain, int const& idim, FullComplexEnumeration<T,Tint,Tgroup> const& fce, std::ostream& os) {
#ifdef SANITY_CHECK_PERFECT_COMPLEX
  std::vector<PerfectFaceEntry<T, Tint>> chain3 = compute_boundary(chain, idim+1, fce, os);
  if (chain3.size() > 0) {
    std::cerr << "PERFCOMP: The chain should have a zero boundary\n";
    throw TerminalException{1};
  }
#endif
  std::vector<PerfectFaceEntry<T, Tint>> one_sol = contracting_homotopy_kernel(chain, idim, fce, os);
#ifdef SANITY_CHECK_PERFECT_COMPLEX
  std::vector<PerfectFaceEntry<T, Tint>> one_sol_img = compute_boundary(one_sol, idim - 1, fce, os);
  if (is_equal_chain(one_sol_img, chain, idim, fce, os)) {
    std::cerr << "PERFCOMP: The proposed preimage is not a solution\n";
    throw TerminalException{1};
  }
#endif
  return one_sol;
}

// clang-format off
#endif  // SRC_PERFECT_PERFECT_COMPLEX_H_
// clang-format on

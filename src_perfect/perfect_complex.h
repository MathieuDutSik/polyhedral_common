// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_PERFECT_COMPLEX_H_
#define SRC_PERFECT_PERFECT_COMPLEX_H_

#include "perfect_tspace.h"

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
    can be done with the classic algorithm.
    It is the faster algorithm.
  + For the other faces we can apply the
    algorithm of using the top dimensional
    faces.

  Ideally we would like to use both. But this
  has serious implications:
  + To use the second algorithm, we have two choices:
    . Computes all the facets of the cone.
      That approach is used in
      - Lorentzian.g for the perfect cones.
      - SaturationAndStabilizer of ProjectiveSystem.g
        for the symplectic story.
      Notes:
      - We cannot use the lower dimensional cells
        for that equivalence realistically. First
        those not well rounded faces will start occurring
        in relatively high rank. Second, the computation
        of the adjacency is likely a relatively tricky
        matter.
    . Compute the orbit and do some orbit splitting
      when the opportunity arise.
      This is likely much more complicated and would
      require a lot of double coset computation.
      And it would be a gain only for the cones that
      are symmetric and have a lot of facets.
      So, let us forget about that method.
  + When we switch from one representation to the
    other, we need to know one triple.
    So, we need to have at least one triple.
  + Since we are interested in boundaries, we have to keep
    track of that in the computation.

 */

//
// The top dimensional complex
//

template<typename Tint>
struct PerfectAdjInfo {
  MyMatrix<Tint> e_mat;
  int i_adj;
  Face f;
};

template<typename T, typename Tint, typename Tgroup>
struct PerfectFormInfoForComplex {
  MyMatrix<T> gram;
  Tshortest<T, Tint> rec_shv;
  Tgroup GRP; // Group acting on the shortest vectors
  std::vector<PerfectAdjInfo<Tint>> list_adj_info;
};

template<typename T, typename Tint, typename Tgroup>
struct PerfectComplexTopDimInfo {
  std::vector<PerfectFormInfoForComplex<T,Tint,Tgroup>> l_perfect;
  bool only_well_rounded;
  bool compute_boundary;
};


template<typename T, typename Tint, typename Tgroup>
PerfectComplexTopDimInfo<T,Tint,Tgroup> generate_perfect_complex_top_dim_info(std::vector<DatabaseEntry_Serial<PerfectTspace_Obj<T,Tint,Tgroup>,PerfectTspace_AdjO<Tint>>> const& l_tot, bool const& only_well_rounded) {
  std::vector<PerfectFormInfoForComplex<T,Tint,Tgroup>> l_perfect;
  for (auto & ePerf: l_tot) {
    std::vector<PerfectAdjInfo<Tint>> list_adj_info;
    for (auto & eAdj: ePerf.ListAdj) {
      int i_adj = eAdj.iOrb;
      Face f = eAdj.x.eInc;
      MyMatrix<Tint> e_mat = eAdj.x.eBigMat;
      PerfectAdjInfo<Tint> pai{e_mat, i_adj, f};
      list_adj_info.emplace_back(std::move(pai));
    }
    PerfectFormInfoForComplex<T,Tint,Tgroup>> perfect{ePerf.x.Gram, ePerf.x.rec_shv, ePerf.x.grp, std::move(list_adj_info)};
    l_perfect.emplace_back(std::move(perfect));
  }
  return {std::move(l_perfect), only_well_rounded};
}

//
//
//

template<typename Tint>
struct PerfectTriple {
  MyMatrix<Tint> e_mat;
  int i_perfect;
  Face f;
};

template<typename T, typename Tint, typename Tgroup>
struct FacePerfectComplex {
  std::vector<PerfectTriple<Tint>> l_triple; // The containing triples.
  MyMatrix<T> gram;
  Tshortest<T, Tint> rec_shv;
  Tgroup GRP; // Group acting on the
  std::vector<MyMatrix<Tint>> l_gens; // generating set of the stabilizer
  bool is_well_rounded;
};

template<typename Tint>
struct BounEntry {
  MyMatrix<Tint> M;
  int iOrb;
}

template<typename Tint, typename Tgroup>
struct ListBounEntry {
  Tgroup RotationSubgroup;
  std::vector<BoundEntry<Tint>> l_bound;
};

template<typename Tint, typename Tgroup>
struct FullBoundary {
  std::vector<ListBoundEntry<Tint,Tgroup>> ll_bound;
};




template<typename Tint, typename Tgroup>
std::vector<PerfectAdjInfo<Tint>> generate_perfect_adj_info(Face const& f, std::vector<std::pair<typename Tgroup::Telt, MyMatrix<Tint>>> const& l_gens, int const& i_adj) {

}







// clang-format off
#endif  // SRC_PERFECT_PERFECT_COMPLEX_H_
// clang-format on

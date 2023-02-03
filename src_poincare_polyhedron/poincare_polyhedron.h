// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POINCARE_POLYHEDRON_TH_POINCARE_POLYHEDRON_H_
#define SRC_POINCARE_POLYHEDRON_TH_POINCARE_POLYHEDRON_H_
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "QuadField.h"
#include "Namelist.h"
#include "POLY_DirectDualDesc.h"
#include "POLY_PolytopeFct.h"
// clang-format on

//
// First part, the data set used
//

// The data structure for keeîng track of the group elements:
// ---Positive value (so element 0 correspond to 1, X to X+1, ...) are the
// elements themselves.
// ---Negative values correspond to their inverse (that is -1 correspond to
// inverse of generator 0)
// ---0 should never show up in the list.
// DI stands for "Direct or Inverse"
struct TrackGroup {
  std::vector<int> ListDI;
};

// We do operations but we can keep track of what is happening.
template<typename T>
struct PairElt {
  TrackGroup tg;
  MyMatrix<T> mat;
};

TrackGroup ProductTrack(TrackGroup const &tg1, TrackGroup const &tg2) {
  std::vector<int> ListDI = tg1.ListDI;
  ListDI.insert(ListDI.end(), tg2.ListDI.begin(), tg2.ListDI.end());
  return {ListDI};
}

TrackGroup InverseTrack(TrackGroup const &tg) {
  std::vector<int> ListDI_ret;
  size_t len = tg.ListDI.size();
  for (size_t u = 0; u < len; u++) {
    size_t v = len - 1 - u;
    ListDI_ret.push_back(-tg.ListDI[v]);
  }
  return {ListDI_ret};
}


template <typename T>
PairElt<T> ProductPair(PairElt<T> const &p1, PairElt<T> const &p2) {
  return {ProductTrack(p1.tg, p2.tg), p1.mat * p2.mat};
}

template <typename T>
PairElt<T> InversePair(PairElt<T> const &p) {
  return {InverseTrack(p.tg), Inverse(p.mat)};
}

void WriteTrackGroup(std::ofstream & os, TrackGroup const& tg) {
  size_t n_elt = tg.ListDI.size();
  os << n_elt;
  for (size_t i_elt=0; i_elt<n_elt; i_elt++) {
    os << ":" << tg.ListDI[i_elt];
  }
}

template <typename T>
bool operator==(PairElt<T> const& pe1, PairElt<T> const& pe2) {
  return pe1.mat == pe2.mat;
}

namespace std {
  template<typename T>
  struct hash<PairElt<T>> {
    std::size_t operator()(PairElt<T> const& pe) {
      return std::hash<MyMatrix<T>>()(pe.mat);
    }
  };
}

template<typename T>
PairElt<T> GenerateIdentity(int const& n) {
  TrackGroup tg;
  MyMatrix<T> mat = IdentityMat<T>(n);
  return {tg, mat};
}

//
// Second part, the finite group generation
//

template<typename T>
std::vector<PairElt<T>> GroupGeneration(std::vector<PairElt<T>> const& input_l_ent) {
  std::vector<PairElt<T>> l_ent = input_l_ent;
  while(true) {
    std::unordered_set<PairElt<T>> GenElt;
    size_t n_ent = l_ent.size();
    for (size_t i_ent=0; i_ent<n_ent; i_ent++) {
      PairElt<T> const& pe1 = l_ent[i_ent];
      for (size_t j_ent=0; j_ent<n_ent; j_ent++) {
        PairElt<T> const& pe2 = l_ent[j_ent];
        PairElt<T> prod = ProductPair(pe1, pe2);
        GenElt.insert(prod);
      }
    }
    if (GenElt.size() == n_ent) {
      return input_l_ent;
    }
    l_ent.clear();
    for (auto & e_ent : GenElt)
      l_ent.push_back(e_ent);
  }
}


// a right coset is of the form Ug
template<typename T>
std::vector<PairElt<T>> IdentifyRightCosets(std::vector<PairElt<T>> const& l_ent, std::vector<PairElt<T>> const& list_grp_elt) {
  std::unordered_set<PairElt<T>> s_coset;
  auto f_insert=[&](PairElt<T> const& pe) -> void {
    for (auto & e_grp_elt : list_grp_elt) {
      PairElt<T> prod = ProductPair(e_grp_elt, pe);
      if (s_coset.count(prod) == 1)
        break;
    }
    s_coset.insert(pe);
  };
  for (auto & pe : l_ent)
    f_insert(pe);
  std::vector<PairElt<T>> l_coset;
  for (auto & e_coset : s_coset)
    l_coset.push_back(e_coset);
  return l_coset;
}

//
// The elementary data structures for the Poincare stuff
//


// The initial data for the Poincare Polyhedron Theorem
// ---a point x
// ---a list of group element which is of course assumed to generate the group
template <typename T> struct DataPoincare {
  MyVector<T> x;
  std::vector<PairElt<T>> ListGroupElt;
};


template <typename T> struct StepEnum {
  std::vector<PairElt<T>> stabilizerElt;
  std::vector<PairElt<T>> ListNeighbor;
};

template <typename T>
StepEnum<T> BuildInitialStepEnum(DataPoincare<T> const& dp) {
  std::vector<PairElt<T>> l_stab_elt;
  MyVector<T> x = dp.x;
  int n = x.size();
  l_stab_elt.push_back(GenerateIdentity<T>(n));
  for (auto & e_elt : dp.ListGroupElt) {
    // TODO: Check that part
    MyVector<T> x_img = e_elt.mat.transpose() * x;
    if (x_img == x) {
      l_stab_elt.push_back(e_elt);
    }
  }
  std::vector<PairElt<T>> stabilizerElt = GroupGeneration(l_stab_elt);
  std::vector<PairElt<T>> ListNeighbor1 = IdentifyRightCosets(dp.ListGroupElt, stabilizerElt);
  std::vector<PairElt<T>> ListNeighbor2;
  for (auto & e_grp_elt : stabilizerElt) {
    for (auto & e_coset : ListNeighbor1) {
      PairElt<T> prod = ProductPair(e_grp_elt, e_coset);
      ListNeighbor2.push_back(prod);
    }
  }
  return {stabilizerElt, ListNeighbor2};
}

template <typename T>
MyMatrix<T> Contragredient(MyMatrix<T> const& M)
{
  return Inverse(TransposedMat(M));
}


//
// Now the polyhedral stuff
//


// This is for the single adjacency in the polyhedron
// * iFaceAdj is the index of the facet in the polyhedron
//   which is adjacent to it in the ridge.
// * iPolyAdj is the index in the adjacent iFaceAdj face
//   in the ridge.
// * iFaceOpp is the index of the facet in the mapped polyhedron
// * iPolyAdj is similarly the corresponding index of the ridge
// * EXTadj is the record of the vertices of the codimension 2 face
struct TsingAdj {
  size_t iFaceAdj;
  size_t iPolyAdj;
  size_t iFaceOpp;
  size_t iPolyOpp;
  Face IncdRidge;
};



struct Tfacet {
  Face IncdFacet;
  std::vector<TsingAdj> l_sing_adj;
};



template <typename T> struct AdjacencyInfo {
  MyMatrix<T> EXT;
  StepEnum<T> se_red;
  std::vector<Tfacet> ll_adj;
};

// The domain is defined originally as
// Tr(AX) <= Tr(PAP^T X)
// which we rewrite
// a.x <= phi(a).x        0 <= (phi(a) - a).x
//
// Matrixwise the scalar product a.x is rewritten as
// A X^T
// phi(a) is expressed as AQ
// 0 <= (AQ - A) X^T
// The mapping A ----> AQ maps to the adjacent domain.
// The mapping of the X ----> X c(Q)^T with c(Q) the
// contragredient representation.
template <typename T>
AdjacencyInfo<T> ComputeAdjacencyInfo(MyVector<T> const &x,
                                      StepEnum<T> const &se,
                                      std::string const &eCommand) {
  size_t miss_val = std::numeric_limits<size_t>::max();
  int n = x.size();
  int n_mat = se.ListNeighbor.size();
  MyMatrix<T> FAC(n_mat, n);
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    MyMatrix<T> eMat = se.ListNeighbor[i_mat].mat;
    MyVector<T> x_img = eMat.transpose() * x;
    MyVector<T> x_diff = x_img - x;
    AssignMatrixRow(FAC, i_mat, x_diff);
  }
  vectface vf = DualDescExternalProgram(FAC, eCommand, std::cerr);
  int n_ext = vf.size();
  MyMatrix<T> EXT(n_ext, n);
  int pos = 0;
  std::vector<Face> v_red(n_mat, Face(n_ext));
  for (auto &eInc : vf) {
    MyVector<T> eEXT = FindFacetInequality(FAC, eInc);
    AssignMatrixRow(EXT, pos, eEXT);
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      v_red[i_mat][pos] = eInc[i_mat];
    }
    pos++;
  }
  std::vector<PairElt<T>> ListNeighborRed;
  int n_mat_red = 0;
  std::vector<int> l_i_mat;
  std::map<Face, size_t> s_facet;
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    MyMatrix<T> EXT_red = SelectRow(EXT, v_red[i_mat]);
    int rnk = RankMat(EXT_red);
    if (rnk == n - 1) {
      ListNeighborRed.push_back(se.ListNeighbor[i_mat]);
      l_i_mat.push_back(n_mat_red);
      s_facet[v_red[i_mat]] = n_mat_red + 1;
      n_mat_red++;
    }
  }
  std::cerr << "n_mat_red=" << n_mat_red << " |s_facet|=" << s_facet.size() << "\n";
  std::cerr << "First part: adjacency structure within the polyhedron\n";
  std::vector<Tfacet> ll_adj;
  for (int i_mat_red = 0; i_mat_red < n_mat_red; i_mat_red++) {
    int i_mat = l_i_mat[i_mat_red];
    Face const &f1 = v_red[i_mat];
    std::vector<TsingAdj> l_adj;
    for (int j_mat_red = 0; j_mat_red < n_mat_red; j_mat_red++) {
      if (i_mat_red != j_mat_red) {
        int j_mat = l_i_mat[j_mat_red];
        Face const &f2 = v_red[j_mat];
        Face f(n_ext);
        for (int i_ext = 0; i_ext < n_ext; i_ext++)
          if (f1[i_ext] == 1 && f2[i_ext] == 1)
            f[i_ext] = 1;
        MyMatrix<T> EXT_red = SelectRow(EXT, f);
        int rnk = RankMat(EXT_red);
        if (rnk == n - 2) {
          l_adj.push_back({static_cast<size_t>(j_mat_red), miss_val, miss_val, miss_val, f});
        }
      }
    }
    Face IncdFacet = v_red[i_mat];
    ll_adj.push_back({IncdFacet, l_adj});
  }
  auto get_iPoly=[&](size_t iFace, Face const& f1) -> size_t {
    size_t n_adjB = ll_adj[iFace].l_sing_adj.size();
    for (size_t i_adjB=0; i_adjB<n_adjB; i_adjB++) {
      Face f2 = ll_adj[iFace].l_sing_adj[i_adjB].IncdRidge;
      if (f1 == f2)
        return i_adjB;
    }
    std::cerr << "Failed to find a matching for f1\n";
    throw TerminalException{1};
  };
  for (int i_mat_red = 0; i_mat_red < n_mat_red; i_mat_red++) {
    size_t n_adj = ll_adj[i_mat_red].l_sing_adj.size();
    for (size_t i_adj=0; i_adj<n_adj; i_adj++) {
      size_t iFaceAdj = ll_adj[i_mat_red].l_sing_adj[i_adj].iFaceAdj;
      Face f1 = ll_adj[i_mat_red].l_sing_adj[i_adj].IncdRidge;
      ll_adj[i_mat_red].l_sing_adj[i_adj].iPolyAdj = get_iPoly(iFaceAdj, f1);
    }
  }
  std::cerr << "Second part: computing the opposite facets\n";
  MyMatrix<T> VectorContain(1, n_ext);
  for (int i_mat_red = 0; i_mat_red < n_mat_red; i_mat_red++) {
    size_t n_adj = ll_adj[i_mat_red].l_sing_adj.size();
    int i_mat = l_i_mat[i_mat_red];
    MyMatrix<T> Q = se.ListNeighbor[i_mat].mat;
    MyMatrix<T> cQ = Contragredient(Q);
    MyMatrix<T> EXTimg = EXT * cQ;
    ContainerMatrix<T> Cont(EXTimg, VectorContain);
    Face MapFace(n_ext);
    std::vector<size_t> l_idx(n_ext,0);
    std::map<int,size_t> map_index;
    for (int i_ext=0; i_ext<n_ext; i_ext++) {
      if (ll_adj[i_mat_red].IncdFacet[i_ext] == 1) {
        for (int i=0; i<n; i++)
          VectorContain(i) = EXT(i_ext, i);
        std::optional<size_t> opt = Cont.GetIdx();
        if (opt) {
          size_t idx = *opt;
          MapFace[idx] = 1;
          map_index[i_ext] = idx;
        }
      }
    }
    size_t pos = s_facet[MapFace];
    if (pos == 0) {
      std::cerr << "Failed to find the MapFace in s_facet\n";
      throw TerminalException{1};
    }
    size_t iFaceOpp = pos - 1;
    for (size_t i_adj=0; i_adj<n_adj; i_adj++) {
      Face f_map(n_ext);
      for (int i_ext=0; i_ext<n_ext; i_ext++) {
        if (ll_adj[i_mat_red].l_sing_adj[i_ext].IncdRidge[i_ext] == 1) {
          size_t idx = map_index[i_ext];
          f_map[idx] = 1;
        }
      }
      size_t iPolyOpp = get_iPoly(iFaceOpp, f_map);
      ll_adj[i_mat_red].l_sing_adj[i_adj].iFaceOpp = iFaceOpp;
      ll_adj[i_mat_red].l_sing_adj[i_adj].iPolyOpp = iPolyOpp;
    }
  }
  StepEnum<T> se_red = {ListNeighborRed};
  return {EXT, se_red, ll_adj};
}

template <typename T>
std::optional<StepEnum<T>> ComputeMissingNeighbors(AdjacencyInfo<T> const &ai) {
  return {};
}

struct RecOption {
  std::string eCommand;
  std::string FileI;
  std::string FileO;
  std::string Arithmetic;
  int n_iter_max;
};



FullNamelist NAMELIST_GetPoincareInput() {
  std::map<std::string, SingleBlock> ListBlock;
  // METHOD
  std::map<std::string, std::string> ListIntValues2_doc;
  std::map<std::string, std::string> ListStringValues2_doc;
  ListIntValues2_doc["n_iter_max"] = "Default: -1\n\
The maximum number of iteration. If negative then infinite";
  ListStringValues2_doc["eCommand"] = "eCommand: lrs\n\
The serial program for computing the dual description. Possibilities: lrs, cdd";
  ListStringValues2_doc["FileI"] = "The input file of the computation";
  ListStringValues2_doc["FileO"] = "The output file of the computation";
  ListStringValues2_doc["Arithmetic"] = "Default: rational\n\
Other possibilities are Qsqrt2, Qsqrt5 and RealAlgebraic=FileDesc where FileDesc is the description";
  SingleBlock BlockPROC;
  BlockPROC.setListIntValues(ListIntValues2_doc);
  BlockPROC.setListStringValues(ListStringValues2_doc);
  ListBlock["PROC"] = BlockPROC;
  // Merging all data
  return {std::move(ListBlock), "undefined"};
}



template <typename T>
DataPoincare<T> ReadDataPoincare(std::string const& FileI) {
  IsExistingFileDie(FileI);
  std::ifstream is(FileI);
  MyVector<T> x = ReadVector<T>(is);
  size_t n_elt;
  is >> n_elt;
  std::vector<PairElt<T>> ListGroupElt;
  for (size_t i_elt=0; i_elt<n_elt; i_elt++) {
    MyVector<T> eElt = ReadMatrix<T>(is);
    TrackGroup tg{{int(i_elt+1)}};
    PairElt<T> pe{tg, eElt};
    ListGroupElt.push_back(pe);
  }
  return {x, ListGroupElt};
}

RecOption ReadInitialData(FullNamelist const& eFull) {
  SingleBlock BlockPROC = eFull.ListBlock.at("DATA");
  std::string eCommand = BlockPROC.ListStringValues.at("eCommand");
  std::string FileI = BlockPROC.ListStringValues.at("FileI");
  std::string FileO = BlockPROC.ListStringValues.at("FileO");
  std::string Arithmetic = BlockPROC.ListStringValues.at("arithmetic");
  int n_iter_max = BlockPROC.ListIntValues.at("n_iter_max");
  return {eCommand, FileI, FileO, Arithmetic, n_iter_max};
}


template <typename T>
AdjacencyInfo<T> IterativePoincareRefinement(DataPoincare<T> const &dp,
                                             RecOption const &rec_option) {
  MyVector<T> x = dp.x;
  StepEnum<T> se = BuildInitialStepEnum(dp);
  int n_iter = 0;
  while (true) {
    std::cerr << "IterativePoincareRefinement n_iter=" << n_iter << "\n";
    AdjacencyInfo<T> ai = ComputeAdjacencyInfo(x, se, rec_option.eCommand);
    std::optional<StepEnum<T>> opt = ComputeMissingNeighbors(ai);
    if (!opt) {
      return ai;
    }
    se = *opt;
    if (rec_option.n_iter_max > 0) {
      if (n_iter > rec_option.n_iter_max) {
        std::cerr << "Reached the maximum number of iterations\n";
        throw TerminalException{1};
      }
    }
    n_iter++;
  }
}


template <typename T>
void PrintAdjacencyInfo(AdjacencyInfo<T> const& ai, std::string const& FileO) {
  std::ofstream os(FileO);
  size_t n_neigh = ai.se_red.ListNeighbor.size();
  os << n_neigh << "\n";
  for (size_t i_neigh=0; i_neigh<n_neigh; i_neigh++) {
    WriteTrackGroup(os, ai.se_red.ListNeighbor[i_neigh].tg);
  }
}


template <typename T>
void full_process_type(RecOption const& rec_option) {
  DataPoincare<T> dp = ReadDataPoincare<T>(rec_option.FileI);
  AdjacencyInfo<T> ai = IterativePoincareRefinement(dp, rec_option);
  PrintAdjacencyInfo(ai, rec_option.FileO);
}





void Process_rec_option(RecOption const& rec_option) {
  std::string arith = rec_option.Arithmetic;
  if (arith == "rational") {
    using T = mpq_class;
    return full_process_type<T>(rec_option);
  }
  if (arith == "Qsqrt5") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return full_process_type<T>(rec_option);
  }
  if (arith == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return full_process_type<T>(rec_option);
  }
  std::optional<std::string> opt_realalgebraic = get_postfix(arith, "RealAlgebraic=");
  if (opt_realalgebraic) {
    using T_rat = mpq_class;
    std::string FileAlgebraicField = *opt_realalgebraic;
    if (!IsExistingFile(FileAlgebraicField)) {
      std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                << " is missing\n";
      throw TerminalException{1};
    }
    HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
    int const idx_real_algebraic_field = 1;
    insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
    using T = RealField<idx_real_algebraic_field>;
    return full_process_type<T>(rec_option);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}








// clang-format off
#endif  // SRC_POINCARE_POLYHEDRON_TH_POINCARE_POLYHEDRON_H_
// clang-format on

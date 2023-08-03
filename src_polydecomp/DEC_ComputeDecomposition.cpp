// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "Decompositions.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    if (argc != 5 && argc != 4) {
      std::cerr << "DEC_ComputeDecomposition [opt] [TheLev] [FileI] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "DEC_ComputeDecomposition [opt] [TheLev] [FileI]\n";
      std::cerr << "\n";
      std::cerr << "  ------ opt -------\n";
      std::cerr << "\n";
      std::cerr << "opt can be one of the following:\n";
      std::cerr << "strategy2_from_high : Compute from the high dimensional to the lower dimensional one\n";
      std::cerr << "    --- We span the lower dimensional by extreme rays\n";
      std::cerr << "    --- We use Plesken-Souvignier like algorithm to decide isomorphism\n";
      std::cerr << "strategy1_from_low  : Compute from the extreme rays to the full one\n";
      std::cerr << "    --- We span the upper dimensional by the facets\n";
      std::cerr << "    --- We use track of the faces to decide isomorphism\n";
      std::cerr << "strategy1_from_high : Not implemented right now\n";
      std::cerr << "\n";
      std::cerr << "  ------ TheLev -------\n";
      std::cerr << "\n";
      std::cerr << "TheLev is the number of depth we compute\n";
      std::cerr << "\n";
      std::cerr << "  ------ FileI -------\n";
      std::cerr << "\n";
      std::cerr << "FileI is the input file\n";
      std::cerr << "\n";
      std::cerr << "  ------ FileO -------\n";
      std::cerr << "\n";
      std::cerr << "FileO is the output file. If absent then the std::cerr is used\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx_value = int32_t;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    std::string opt = argv[1];
    //
    int TheLev;
    sscanf(argv[2], "%d", &TheLev);
    //
    std::string FileI = argv[3];
    std::ifstream is(FileI);
    //
    // The polyhedral cone.
    // So far, we encode all the EXT / FAC and the group
    std::vector<ConeDesc<T, Tint, Tgroup>> ListCones;
    //
    //
    MyMatrix<Tint> G = ReadMatrix<Tint>(is);
    std::cerr << "We have G\n";
    size_t n_domain;
    is >> n_domain;
    for (size_t i = 0; i < n_domain; i++) {
      std::cerr << "i=" << i << " / " << n_domain << "\n";
      MyMatrix<T> EXT = ReadMatrix<T>(is);
      size_t n_ext = EXT.rows();
      std::cerr << "We have read EXT, |EXT|=" << EXT.rows() << "/" << EXT.cols()
                << "\n";
      MyMatrix<Tint> EXT_i = UniversalMatrixConversion<Tint, T>(EXT);
      MyMatrix<T> FAC = ReadMatrix<T>(is);
      size_t n_fac = FAC.rows();
      std::cerr << "We have read FAC, |FAC|=" << FAC.rows() << "/" << FAC.cols()
                << "\n";
      Tgroup GRP_ext = ReadGroup<Tgroup>(is);
      Tgroup GRP_fac = ReadGroup<Tgroup>(is);
      //
      Face extfac_incd = Compute_extfac_incd(FAC, EXT);
      std::cerr << "List(|FAC|) =";
      for (size_t iFac = 0; iFac < n_fac; iFac++) {
        Face f = GetFacet_extfac(extfac_incd, n_fac, n_ext, iFac);
        std::cerr << " " << f.count();
      }
      std::cerr << "\n";
      Face facext_incd = Compute_facext(extfac_incd, FAC.rows(), EXT.rows());
      ConeDesc<T, Tint, Tgroup> eCone{EXT,         EXT_i,   FAC,    extfac_incd,
                                      facext_incd, GRP_ext, GRP_fac};
      ListCones.push_back(eCone);
    }
    //
    std::cerr << "Now computing ListListDomain opt=" << opt << "\n";
    std::vector<std::vector<FaceDesc>> ListListDomain;
    if (opt == "strategy2_from_high") {
      std::cerr << "Matching strategy2_from_high\n";
#ifdef DEBUG_POLYEDRAL_DECOMPOSITION
      std::vector<std::vector<sing_adj<Tint>>> ll_sing_adj =
          compute_adjacency_structure<T, Tint, Tgroup, Tidx_value>(ListCones,
                                                                   G);
#else
      std::vector<std::vector<sing_adj<Tint>>> ll_sing_adj;
#endif
      ListListDomain =
          Compute_ListListDomain_strategy2<T, Tint, Tgroup, Tidx_value>(
              ListCones, G, ll_sing_adj, TheLev);
    }
    if (opt == "strategy1_from_low" || opt == "strategy1_from_high") {
      std::cerr << "Matching strategy1\n";
      std::vector<std::vector<sing_adj<Tint>>> ll_sing_adj =
          compute_adjacency_structure<T, Tint, Tgroup, Tidx_value>(ListCones,
                                                                   G);
      std::cerr << "We have ll_sing_adj\n";
      if (opt == "strategy1_from_high") {
        std::cerr << "Matching strategy1_from_high\n";
        std::cerr << "This needs to be implemented as it does not appear "
                     "needed right now\n";
        // Implementation should use above code
        throw TerminalException{1};
      }
      if (opt == "strategy1_from_low") {
        std::cerr << "Matching strategy1_from_low\n";
        ListListDomain =
            Compute_ListListDomain_strategy1<T, Tint, Tgroup, Tidx_value>(
                ListCones, G, ll_sing_adj, TheLev);
      }
    }
    //
    auto do_print = [&](std::ostream &os) -> void {
      size_t n_lev = ListListDomain.size();
      os << "Number of levels = " << n_lev << "\n";
      for (size_t i_lev = 0; i_lev < n_lev; i_lev++) {
        size_t n_orbit = ListListDomain[i_lev].size();
        os << "Number of orbits at level " << i_lev << " = " << n_orbit << "\n";
        for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
          const FaceDesc &fd = ListListDomain[i_lev][i_orbit];
          const ConeDesc<T, Tint, Tgroup> &eC = ListCones[fd.iCone];
          Face f_ext = Compute_faceEXT_from_faceFAC(
              eC.extfac_incd, eC.FAC.rows(), eC.EXT.rows(), fd.f_fac);
          //
          MyMatrix<Tint> EXT_sel = SelectRow(eC.EXT_i, f_ext);
          os << "i_orbit=" << i_orbit << " |EXT|=" << EXT_sel.rows()
             << " iCone=" << fd.iCone << "\n";
          WriteMatrix(os, EXT_sel);
        }
      }
    };
    //
    if (argc == 4) {
      do_print(std::cerr);
    } else {
      std::string FileO = argv[4];
      std::ofstream os(FileO);
      do_print(os);
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}

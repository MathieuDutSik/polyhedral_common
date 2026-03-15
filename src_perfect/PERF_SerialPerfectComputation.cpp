// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "perfect_complex.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
FullComplexEnumeration<T, Tint, Tgroup> get_full_complex_enumeration_kernel(LinSpaceMatrix<T> const& LinSpa,
                                                                            PerfectComplexOptions const& pco,
                                                                            std::ostream& os) {
  using TintGroup = typename Tgroup::Tint;
  int n = LinSpa.n;
  int dimEXT = LinSpa.ListMat.size();
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, os);
  //
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);
  bool keep_generators = false;
  bool reduce_gram_matrix = false;
  DataPerfectTspace<T, Tint, Tgroup> data{
      LinSpa, OnlineHierarchicalMatrixReduction<Tint>(n, std::cerr),
      keep_generators, reduce_gram_matrix, std::move(rddo)};
  using Tdata = DataPerfectTspaceFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_Serial<Tobj, TadjO>;
  auto f_incorrect = [&]([[maybe_unused]] Tobj const &x) -> bool {
    return false;
  };
  int max_runtime_second = 0;
  std::vector<Tout> l_tot =
      EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
          data_func, f_incorrect, max_runtime_second);
  os << "|l_tot|=" << l_tot.size() << "\n";
  FullComplexEnumeration<T, Tint, Tgroup> fce = full_perfect_complex_enumeration(l_tot, LinSpa, pco, os);
  if (pco.compute_boundary) {
    bool test = are_product_zeros(fce, os);
    os << "are_product_zeros, test=" << GAP_logical(test) << "\n";
  } else {
    os << "no boundary available\n";
  }
  return fce;
}

template <typename T, typename Tint, typename Tgroup>
bool is_correct_fce(FullComplexEnumeration<T, Tint, Tgroup> const& fce,
                    LinSpaceMatrix<T> const& LinSpa,
                    PerfectComplexOptions const& pco) {
  if (!LinSpaceMatrixEqual(LinSpa, fce.pctdi.LinSpa)) {
    return false;
  }
  if (!PerfectComplexOptionsEqual(pco, fce.pctdi.pco)) {
    return false;
  }
  return true;
}


template <typename T, typename Tint, typename Tgroup>
FullComplexEnumeration<T, Tint, Tgroup> get_full_complex_enumeration(LinSpaceMatrix<T> const& LinSpa,
                                                                     PerfectComplexOptions const& pco,
                                                                     std::string const& CacheFile,
                                                                     std::ostream& os) {
  if (CacheFile == "none") {
    return get_full_complex_enumeration_kernel<T,Tint,Tgroup>(LinSpa, pco, os);
  }
  if (IsExistingFile(CacheFile)) {
    FullComplexEnumeration<T, Tint, Tgroup> fce = read_fce_from_file<T,Tint,Tgroup>(CacheFile);
    if (is_correct_fce(fce, LinSpa, pco)) {
      return fce;
    }
  }
  FullComplexEnumeration<T, Tint, Tgroup> fce = get_full_complex_enumeration_kernel<T,Tint,Tgroup>(LinSpa, pco, os);
  write_fce_to_file(CacheFile, fce);
  return fce;
}



template <typename T, typename Tint, typename Tgroup>
void process_A(FullNamelist const &eFull, std::ostream& os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockQUERIES = eFull.get_block("QUERIES");
  SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa =
      ReadTspace<T, Tint, Tgroup>(BlockTSPACE, os);
  //
  bool compute_boundary = BlockDATA.get_bool("ComputeBoundary");
  bool only_well_rounded = BlockDATA.get_bool("OnlyWellRounded");
  bool compute_contracting_homotopy = BlockDATA.get_bool("ComputeContractingHomotopy");
  PerfectComplexOptions pco{only_well_rounded, compute_boundary, compute_contracting_homotopy};
  //
  std::string CacheFile = BlockDATA.get_string("CacheFile");
  FullComplexEnumeration<T, Tint, Tgroup> fce = get_full_complex_enumeration<T,Tint,Tgroup>(LinSpa, pco, CacheFile, os);

  std::string FileGroupGenerators = BlockQUERIES.get_string("FileGroupGenerators");
  if (FileGroupGenerators != "null") {
    std::vector<MyMatrix<Tint>> list_gen =
        get_reduced_generators_fce<T, Tint, Tgroup>(fce, os);
    std::ofstream os_out(FileGroupGenerators);
    os_out << "return ";
    WriteListMatrixGAP(os_out, list_gen);
    os_out << ";\n";
  }

  /*
    The queries related to individual cells
   */
  auto f_cell_oper=[&](std::string const& file_name, std::function<void(std::ostream&,MyMatrix<Tint>)> const&f) -> void {
    std::string file_oper = BlockQUERIES.get_string(file_name);
    if (file_oper != "null") {
      std::vector<MyMatrix<Tint>> l_ext =
        ReadListMatrixFile<Tint>(file_oper);
      std::string OutFile = file_oper + ".output";
      std::ofstream os_out(OutFile);
      bool is_first = true;
      os_out << "return [";
      for (auto &EXT : l_ext) {
        if (!is_first) {
          os_out << ",\n";
        }
        is_first = false;
        MyMatrix<Tint> EXT_sat = vector_family_saturation(EXT, fce.pctdi.LinSpa.PtStabGens);
        f(os_out, EXT_sat);
      }
      os_out << "];\n";
    }
  };
  // Dimension of the cell
  std::function<void(std::ostream&,MyMatrix<Tint> const&)> f_dim=[&](std::ostream& os_out, MyMatrix<Tint> const& EXT) -> void {
    int n_mat = fce.pctdi.LinSpa.ListMat.size();
    MyMatrix<T> ScalMat = get_scal_mat(fce.pctdi.LinSpa.ListMat, EXT);
    int rnk = RankMat(ScalMat);
    int index = n_mat - rnk;
    os_out << "rec(n_mat:=" << n_mat << ", rnk:=" << rnk << ", index:=" << index << ")";
  };
  f_cell_oper("FileDimensions", f_dim);
  // Whether it is a face or not
  std::function<void(std::ostream&,MyMatrix<Tint> const&)> f_is_face=[&](std::ostream& os_out, MyMatrix<Tint> const& EXT) -> void {
    std::optional<MyMatrix<T>> opt =
      is_bounded_face_iterative<T, Tint>(LinSpa, EXT, os);
    bool is_face = opt.has_value();
    os_out << GAP_logical(is_face);
  };
  f_cell_oper("FileIsFace", f_is_face);
  // Whether it is well rounded or not
  std::function<void(std::ostream&,MyMatrix<Tint> const&)> f_wr=[&](std::ostream& os_out, MyMatrix<Tint> const& EXT) -> void {
    PerfectBoundednessProperty pbp = initial_bounded_property(fce.pctdi.LinSpa, EXT, os);
    bool is_well_rounded = get_result(pbp);
    os_out << GAP_logical(is_well_rounded);
  };
  f_cell_oper("FileIsWellRounded", f_wr);
  // Finding the index and position in the cell complex
  std::function<void(std::ostream&,MyMatrix<Tint> const&)> f_ffs=[&](std::ostream& os_out, MyMatrix<Tint> const& EXT) -> void {
    FceFaceSearch<Tint> ffs = fce_face_search(fce, EXT, os);
    int iOrb = ffs.iOrb + 1;
    os_out << "rec(index:=" << ffs.index << ", iOrb:=" << iOrb << ", M:=" << StringMatrixGAP(ffs.M) << ")";
  };
  f_cell_oper("FileFaceSearch", f_ffs);
  // Determining the stabilizer
  std::function<void(std::ostream&,MyMatrix<Tint> const&)> f_stab=[&](std::ostream& os_out, MyMatrix<Tint> const& EXT) -> void {
    auto result = compute_stabilizer_ext<T, Tint, Tgroup>(fce, EXT, os);
    os_out << "Group(";
    WriteListMatrixGAP(os_out, result.second);
    os_out << ")";
  };
  f_cell_oper("FileStabilizerQueries", f_stab);
  // Determining the lower boundary
  std::function<void(std::ostream&,MyMatrix<Tint> const&)> f_lower=[&](std::ostream& os_out, MyMatrix<Tint> const& EXT) -> void {
    FceFaceSearch<Tint> ffs = fce_face_search(fce, EXT, os);
    std::vector<BoundEntry<Tint>> const& l_bound = fce.boundaries[ffs.index].ll_bound[ffs.iOrb].l_bound;
    std::vector<MyMatrix<Tint>> l_extbnd;
    for (auto & entry: l_bound) {
      MyMatrix<Tint> const& EXT1 = fce.levels[ffs.index+1].l_faces[entry.iOrb].EXT;
      MyMatrix<Tint> EXT2 = EXT1 * entry.M * ffs.M;
      l_extbnd.emplace_back(std::move(EXT2));
    }
    WriteListMatrixGAP(os_out, l_extbnd);
  };
  f_cell_oper("FileCellLowerBoundary", f_lower);
  // Determining the upper boundary
  std::function<void(std::ostream&,MyMatrix<Tint> const&)> f_upper=[&](std::ostream& os_out, MyMatrix<Tint> const& EXT) -> void {
    FceFaceSearch<Tint> ffs = fce_face_search(fce, EXT, os);
    std::pair<std::vector<MyMatrix<Tint>>, std::vector<PerfectFace<Tint>>> pair =
      get_all_upper_faces(fce, ffs.index, ffs.iOrb, os);
    std::vector<MyMatrix<Tint>> l_ext_upper;
    for (auto & EXT1: pair.first) {
      MyMatrix<Tint> EXT2 = EXT1 * ffs.M;
      l_ext_upper.emplace_back(std::move(EXT2));
    }
    WriteListMatrixGAP(os_out, l_ext_upper);
  };
  f_cell_oper("FileCellUpperBoundary", f_upper);
  // Computing the equivalence
  std::string FileEquivalenceQueries = BlockQUERIES.get_string("FileEquivalenceQueries");
  if (FileEquivalenceQueries != "null") {
    std::vector<MyMatrix<Tint>> l_ext =
      ReadListMatrixFile<Tint>(FileEquivalenceQueries);
    std::string OutFile = FileEquivalenceQueries + ".output";
    std::ofstream os_out(OutFile);
    size_t n_case = l_ext.size() / 2;
    os_out << "return [";
    for (size_t i_case=0; i_case<n_case; i_case++) {
      if (i_case > 0) {
        os_out << ",\n";
      }
      MyMatrix<Tint> const& EXT1 = l_ext[2*i_case];
      MyMatrix<Tint> const& EXT2 = l_ext[2*i_case+1];
      std::optional<MyMatrix<Tint>> opt =
        find_equivalence_ext(fce, EXT1, EXT2, os);
      if (opt) {
        WriteMatrixGAP(os_out, *opt);
      } else {
        os_out << "fail";
      }
    }
    os_out << "];\n";
  }
  std::string FileCells = BlockQUERIES.get_string("FileCells");
  if (FileCells != "null") {
    int index = BlockQUERIES.get_int("IndexCell");
    int n_orb = fce.levels[index].l_faces.size();
    std::ofstream os_out(FileCells);
    os_out << "return [";
    for (int iOrb=0; iOrb<n_orb; iOrb++) {
      if (iOrb > 0) {
        os_out << ",\n";
      }
      MyMatrix<Tint> const& EXT = fce.levels[index].l_faces[iOrb].EXT;
      os_out << "rec(EXT:=";
      WriteMatrixGAP(os_out, EXT);
      os_out << ")";
    }
    os_out << "];\n";
  }
  std::string FileListUpperBoundary = BlockQUERIES.get_string("FileListUpperBoundary");
  if (FileListUpperBoundary != "null") {
    int index = BlockQUERIES.get_int("IndexUpperBoundary");
    int n_orb = fce.levels[index].l_faces.size();
    std::ofstream os_out(FileListUpperBoundary);
    os_out << "return [";
    for (int iOrb=0; iOrb<n_orb; iOrb++) {
      if (iOrb > 0) {
        os_out << ",\n";
      }
      os_out << "rec(iOrb:=" << (iOrb + 1);
      std::pair<std::vector<MyMatrix<Tint>>, std::vector<PerfectFace<Tint>>> pair =
        get_all_upper_faces(fce, index, iOrb, os);
      os_out << ", ListEXT:=";
      WriteListMatrixGAP(os_out, pair.first);
      os_out << ", ListMap:=[";
      int len = pair.second.size();
      for (int i=0; i<len; i++) {
        if (i>0) {
          os_out << ",\n";
        }
        int jOrb = pair.second[i].iOrb + 1;
        os_out << "rec(jOrb:=" << jOrb << ", M:=";
        WriteMatrixGAP(os_out, pair.second[i].M);
        os_out << ")";
      }
      os_out << "]";
      os_out << ")";
    }
    os_out << "];\n";
  }
  std::string FileListLowerBoundary = BlockQUERIES.get_string("FileListLowerBoundary");
  if (FileListLowerBoundary != "null") {
    int index = BlockQUERIES.get_int("IndexLowerBoundary");
    int n_orb = fce.levels[index].l_faces.size();
    std::ofstream os_out(FileListLowerBoundary);
    os_out << "return [";
    for (int iOrb=0; iOrb<n_orb; iOrb++) {
      if (iOrb > 0) {
        os_out << ",\n";
      }
      std::vector<BoundEntry<Tint>> const& l_bound = fce.boundaries[index].ll_bound[iOrb].l_bound;
      os_out << "rec(iOrb:=" << iOrb;
      MyMatrix<Tint> EXT = fce.levels[index].l_faces[iOrb].EXT;
      os_out << ", EXT:=";
      WriteMatrixGAP(os_out, EXT);
      //
      int n_bnd = l_bound.size();
      os_out << ", ListBnd:=[";
      for (int i_bnd=0; i_bnd<n_bnd; i_bnd++) {
        if (i_bnd > 0) {
          os_out << ",\n";
        }
        int jOrb = l_bound[i_bnd].iOrb;
        MyMatrix<Tint> const& M = l_bound[i_bnd].M;
        MyMatrix<Tint> EXT = fce.levels[index-1].l_faces[jOrb].EXT * M;
        os_out << "rec(sign:=" << l_bound[i_bnd].sign;
        os_out << ", jOrb:=" << jOrb;
        os_out << ", M:=";
        WriteMatrixGAP(os_out, l_bound[i_bnd].M);
        os_out << ", EXT:=";
        WriteMatrixGAP(os_out, EXT);
        os_out << ")";
      }
      os_out << "]";
      os_out << ")";
    }
    os_out << "];\n";
  }
  /*
    Now the operations on chains
   */
  using Tchain = std::vector<PerfectFaceEntry<T, Tint>>;
  auto f_chain_oper=[&](std::string const& file_name, std::function<Tchain(int,Tchain const&)> const&f) -> void {
    std::string file_oper = BlockQUERIES.get_string(file_name);
    if (file_oper != "null") {
      std::ifstream is(file_oper);
      int index, n_chains;
      is >> index;
      is >> n_chains;
#ifdef DEBUG_PERFECT_COMPLEX
      os << "PERF: f_chain_oper, index=" << index << " n_chains=" << n_chains << "\n";
#endif
      std::string OutFile = file_oper + ".output";
      std::ofstream os_out(OutFile);
      os_out << "return [";
      for (int i_chain=0; i_chain<n_chains; i_chain++) {
        if (i_chain > 0) {
          os_out << ",\n";
        }
        std::vector<PerfectFaceEntry<T, Tint>> chain1 = ReadListPerfectFaceEntry<T,Tint>(is);
#ifdef DEBUG_PERFECT_COMPLEX
        os << "PERF: f_chain_oper, i_chain=" << i_chain << " |chain1|=" << chain1.size() << "\n";
#endif
        std::vector<PerfectFaceEntry<T, Tint>> chain2 = f(index, chain1);
        WriteListPerfectFaceEntryGAP(os_out, chain2);
      }
      os_out << "];\n";
    }
  };
  // The contracting homotopy.
  std::function<Tchain(int,Tchain const&)> f_ch=[&](int index, Tchain const& chain) -> Tchain {
    return contracting_homotopy(index, chain, fce, os);
  };
  f_chain_oper("FileChainContractingHomotopy", f_ch);
  // The chain boundary
  std::function<Tchain(int,Tchain const&)> f_bnd=[&](int index, Tchain const& chain) -> Tchain {
    return chain_boundary(index, chain, fce, os);
  };
  f_chain_oper("FileChainBoundary", f_bnd);
  // The chain simplification
  std::function<Tchain(int,Tchain const&)> f_simp=[&](int index, Tchain const& chain) -> Tchain {
    return chain_simplification(index, chain, fce, os);
  };
  f_chain_oper("FileChainSimplification", f_simp);
}

template <typename T, typename Tint> void process_B(FullNamelist const &eFull) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  return process_A<T, Tint, Tgroup>(eFull, std::cerr);
}

template <typename T> void process_C(FullNamelist const &eFull) {
  std::string arithmetic_Tint =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic_Tint");
  if (arithmetic_Tint == "gmp_integer") {
    using Tint = mpz_class;
    return process_B<T, Tint>(eFull);
  }
  std::cerr
      << "PERF_SerialEnumeratePerfectCones: Failed to find a matching type for "
         "arithmetic_Tint="
      << arithmetic_Tint << "\n";
  std::cerr << "Available types: gmp_integer\n";
  throw TerminalException{1};
}

void process_D(FullNamelist const &eFull) {
  std::string arithmetic_T =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic_T");
  if (arithmetic_T == "gmp_rational") {
    using T = mpq_class;
    return process_C<T>(eFull);
  }
  std::cerr
      << "PERF_SerialEnumeratePerfectCones: Failed to find a matching type for "
         "arithmetic_T="
      << arithmetic_T << "\n";
  std::cerr << "Available types: gmp_rational\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_ENUMERATE_PERFECT_COMPLEX_TSPACE();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "PERF_SerialPerfectComputation [file.nml]\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_D(eFull);
    std::cerr << "Normal termination of PERF_SerialPerfectComputation\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PERF_SerialPerfectComputation\n";
    exit(e.eVal);
  }
  runtime(time);
}

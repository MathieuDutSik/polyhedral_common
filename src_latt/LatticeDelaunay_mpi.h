// Copyright (C) 2024 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEDELAUNAY_MPI_H_
#define SRC_LATT_LATTICEDELAUNAY_MPI_H_

// clang-format off
#include "LatticeDelaunay.h"
#include "POLY_MPI_AdjacencyScheme.h"
// clang-format on

#ifdef TIMINGS
#define TIMINGS_MPI_DELAUNAY_ENUMERATION
#endif

#ifdef DEBUG
#define DEBUG_MPI_DELAUNAY_ENUMERATION
#endif


template <typename T, typename Tvert, typename Tgroup>
void WriteFamilyDelaunay_Mpi(
    boost::mpi::communicator &comm,
    MyMatrix<T> const& GramMat,
    std::string const &OutFormat,
    std::ostream& os_out,
    std::vector<DatabaseEntry_MPI<
        typename DataLatticeFunc<T, Tvert, Tgroup>::Tobj,
        typename DataLatticeFunc<T, Tvert, Tgroup>::TadjO>> const &ListDel,
    std::ostream &os) {
  int i_rank = comm.rank();
  if (OutFormat == "nothing") {
    std::cerr << "No output\n";
    return;
  }
  using Tout =
      DatabaseEntry_Serial<typename DataLatticeFunc<T, Tvert, Tgroup>::Tobj,
                           typename DataLatticeFunc<T, Tvert, Tgroup>::TadjO>;
  if (OutFormat == "CheckMergedOutput") {
    int i_proc_out = 0;
    std::vector<Tout> l_ent = my_mpi_gather(comm, ListDel, i_proc_out);
    if (i_proc_out == i_rank) {
      DelaunayTesselation<Tvert, Tgroup> DT =
          DelaunayTesselation_From_DatabaseEntries_Serial<T, Tvert, Tgroup>(l_ent);
      check_delaunay_tessellation(DT, os);
    }
    std::cerr << "The Delaunay tesselation passed the adjacency check\n";
    return;
  }
  if (OutFormat == "GAP") {
    int i_proc_out = 0;
    std::vector<Tout> l_ent = my_mpi_gather(comm, ListDel, i_proc_out);
    if (i_proc_out == i_rank) {
      DelaunayTesselation<Tvert, Tgroup> DT =
          DelaunayTesselation_From_DatabaseEntries_Serial<T, Tvert, Tgroup>(l_ent);
      os_out << "return ";
      WriteEntryGAP(os_out, DT);
      os_out << ";\n";
    }
    std::cerr << "The Delaunay tesselation has been written to file\n";
    return;
  }
  if (OutFormat == "RAW") {
    int nbDel = ListDel.size();
    os_out << "nbDel=" << nbDel << "\n";
    for (int iDel = 0; iDel < nbDel; iDel++) {
      os_out << "iDel=" << iDel << "/" << nbDel << "\n";
      WriteMatrix(os_out, ListDel[iDel].x.EXT);
    }
  }
  if (OutFormat == "GAP_Covering") {
    T max_radius(0);
    for (auto & eDel : ListDel) {
      MyMatrix<Tvert> const& EXT = eDel.x.EXT;
      MyMatrix<T> EXT_T = UniversalMatrixConversion<T,Tvert>(EXT);
      CP<T> cp = CenterRadiusDelaunayPolytopeGeneral<T>(GramMat, EXT_T);
      T SquareRadius = cp.SquareRadius;
      if (SquareRadius > max_radius) {
        max_radius = SquareRadius;
      }
    }
    int i_proc_out = 0;
    std::vector<T> l_ent = my_mpi_gather(comm, max_radius, i_proc_out);
    if (i_proc_out == i_rank) {
      int TheDim = GramMat.rows();
      T TheDet = DeterminantMat(GramMat);
      T TheCov(0);
      for (auto & rad: l_ent) {
        if (rad > TheCov) {
          TheCov = rad;
        }
      }
      ResultCov<T> x = ComputeCoveringDensityFromDimDetCov<T>(TheDim, TheDet, TheCov);
      os_out << "return " << to_stringGAP(x) << ";\n";
    }
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

template <typename T, typename Tint, typename Tgroup>
void ComputeDelaunayPolytope_MPI(boost::mpi::communicator &comm,
                             FullNamelist const &eFull) {
  using TintGroup = typename Tgroup::Tint;
  std::unique_ptr<std::ofstream> os_ptr = get_mpi_log_stream(comm, eFull);
  std::ostream &os = *os_ptr;
  SingleBlock const& BlockDATA = eFull.get_block("DATA");
  std::string GRAMfile = BlockDATA.get_string("GRAMfile");
  MyMatrix<T> GramMat = ReadMatrixFile<T>(GRAMfile);
  int dimEXT = GramMat.rows() + 1;
  SingleBlock const& BlockSTORAGE = eFull.get_block("STORAGE");
  //
  bool STORAGE_Saving = BlockSTORAGE.get_bool("Saving");
  std::string STORAGE_Prefix = BlockSTORAGE.get_string("Prefix");
  CreateDirectory(STORAGE_Prefix);
  //
  std::string OutFormat = BlockDATA.get_string("OutFormat");
  std::string OutFile = BlockDATA.get_string("OutFile");
#ifdef DEBUG_MPI_DELAUNAY_ENUMERATION
  std::cerr << "MPI_DEL_ENUM: OutFormat=" << OutFormat << " OutFile=" << OutFile << "\n";
#endif
  int max_runtime_second = BlockDATA.get_int("max_runtime_second");
#ifdef DEBUG_MPI_DELAUNAY_ENUMERATION
  std::cerr << "MPI_DEL_ENUM: max_runtime_second=" << max_runtime_second << "\n";
#endif
  //
  std::string FileDualDesc =
    BlockDATA.get_string("FileDualDescription");
  PolyHeuristicSerial<TintGroup> AllArr =
      Read_AllStandardHeuristicSerial_File<T, TintGroup>(FileDualDesc, dimEXT,
                                                         os);
  DataLattice<T, Tint, Tgroup> data = get_data_lattice<T,Tint,Tgroup>(eFull, AllArr, os);
  using Tdata = DataLatticeFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_MPI<Tobj, TadjO>;
  //
  std::pair<bool, std::vector<Tout>> pair = EnumerateAndStore_MPI<Tdata>(
      comm, data_func, STORAGE_Prefix, STORAGE_Saving, max_runtime_second);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "MPI_DEL_ENUM: We now have IsFinished=" << pair.first << "\n";
  os << "MPI_DEL_ENUM: We now have ListDel |ListDel|=" << pair.second.size()
     << "\n";
#endif
  //
  if (pair.first) {
    auto f=[&](std::ostream& os_out) -> void {
      WriteFamilyDelaunay_Mpi<T, Tint, Tgroup>(comm, GramMat, OutFormat, os_out, pair.second,
                                               os);
    };
    print_stderr_stdout_file(OutFile, f);
  } else {
    os << "MPI_DEL_ENUM: The enumeration did not finish\n";
  }
}

// clang-format off
#endif  // SRC_LATT_LATTICEDELAUNAY_MPI_H_
// clang-format on

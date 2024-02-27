// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ISODELAUNAYDOMAINS_H_
#define SRC_LATT_ISODELAUNAYDOMAINS_H_

// clang-format off
#include "NumberTheory.h"
#include "CtypeMPI_types.h"
#include "Namelist.h"
#include "hash_functions.h"
#include "rational.h"
#include "sparse-map/include/tsl/sparse_map.h"
#include <boost/mpi.hpp>
#include <netcdf>
#include <unordered_map>
// clang-format on

template<typename Tvert, typename Tgroup>
struct DataCtype {
  int n;
};



/*
  We implement the hash of a Delaunay tessellationn. The constraint is that two different
  but equivalent tessellations, must have the same hash. Hopefully, this is not a problem
  since we have many invariants:
  * The number of vertices of the orbit representatives of Delaunay polytopes.
  * The size of their automorphism groups
  * The number of vertices of the orbit representative of their facets.
 */
template<typename Tvert, typename Tgroup>
size_t CTYP_ComputeInvariant(DataCtype<Tvert,Tgroup> const& data,
                             size_t const& seed,
                             TypeCtypeExch<T> const &eCtype,
                             [[maybe_unused]] std::ostream & os) {

}


struct AdjInfo {
  int iOrb;
  int i_adj;
};

template<typename T>
struct FullAdjInfo {
  MyVector<T> eIneq;
  std::vector<AdjInfo> ListAdjInfo;
};

FullNamelist NAMELIST_GetStandard_COMPUTE_LATTICE_IsoDelaunayDomains() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["arithmetic_T"] = "gmp_rational";
  ListStringValues1["arithmetic_Tint"] = "gmp_integer";
  ListStringValues1["OutFormat"] = "nothing";
  ListStringValues1["OutFile"] = "unset.out";
  ListStringValues1["FileDualDescription"] = "unset";
  ListIntValues1["max_runtime_second"] = 0;
  ListBoolValues1["ApplyStdUnitbuf"] = false;
  ListBoolValues1["Saving"] = false;
  ListStringValues1["Prefix"] = "/irrelevant/";
  ListStringValues1["CVPmethod"] = "SVexact";
  SingleBlock BlockDATA;
  BlockDATA.ListIntValues = ListIntValues1;
  BlockDATA.ListBoolValues = ListBoolValues1;
  BlockDATA.ListDoubleValues = ListDoubleValues1;
  BlockDATA.ListStringValues = ListStringValues1;
  BlockDATA.ListListStringValues = ListListStringValues1;
  ListBlock["DATA"] = BlockDATA;
  // TSPACE
  ListBlock["TSPACE"] = SINGLEBLOCK_Get_Tspace_Description();
  // Merging all data
  return {ListBlock, "undefined"};
}

template<typename T, typename Tint>
struct IsoEdgeDomain_MPI_AdjO {
  int iProc;
  int iOrb;
};

namespace boost::serialization {
  template <class Archive, typename T, typename Tint>
  inline void serialize(Archive &ar, IsoEdgeDomain_MPI_AdjO<T, Tint> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("iProc", eRec.iProc);
    ar &make_nvp("iOrb", eRec.iOrb);
  }
}

template<typename T, typename Tint, typename Tgroup>
struct IsoEdgeDomain_AdjI {
  MyVector<T> V;
};

namespace boost::serialization {
  template <class Archive, typename T, typename Tint, typename Tgroup>
  inline void serialize(Archive &ar, IsoEdgeDomain_AdjI<T, Tint, Tgroup> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("V", eRec.V);
  }
}

template<typename Tint>
struct IsoDelaunayDomain_AdjO {
  int iOrb;
};

namespace boost::serialization {
  template <class Archive, typename Tint>
  inline void serialize(Archive &ar, IsoDelaunayDomain_AdjO<Tint> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("iOrb", eRec.iOrb);
  }
}

template<typename T, typename Tint, typename Tgroup>
struct IsoDelaunayDomain_MPI_Entry {
  TypeCtypeExch<Tint> ctype_arr;
  StructuralInfo struct_info;
  std::vector<IsoEdgeDomain_MPI_AdjO<T, Tint>> ListAdj;
};

namespace boost::serialization {
  template <class Archive, typename T, typename Tint, typename Tgroup>
  inline void serialize(Archive &ar, IsoDelaunayDomain_MPI_Entry<T, Tint, Tgroup> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("ctype_arr", eRec.ctype_arr);
    ar &make_nvp("struct_info", eRec.struct_info);
    ar &make_nvp("ListAdj", eRec.ListAdj);
  }
}

template<typename T, typename Tint, typename Tgroup>
IsoDelaunayDomain<T, Tint, Tgroup> GetInitialIsoDelaunayDomain(DataIsoDelaunayDomains<T,Tint,Tgroup> const& eData) {
  DelaunayTesselation<Tint, Tgroup> DT = GetInitialGenericDelaunayTesselation(eData);
  MyMatrix<T> GramMat = GetInteriorGramMatrix(eData.LinSpa, DT);
  return {DT, GramMat};
}




template<typename T, typename Tint, typename Tgroup>
std::vector<IsoDelaunayDomain_MPI_Entry<T,Tint,Tgroup>> MPI_EnumerationIsoEdgeDomains(boost::mpi::communicator &comm, DataIsoEdgeDomains<T,Tint,Tgroup> & eData, std::ostream & os) {
  using Tobj = TypeCtypeExch<Tvert>;
  using TadjI = IsoEdgeDomain_AdjI<T, Tint, Tgroup>;
  using TadjO = IsoEdgeDomain_MPI_AdjO<T, Tint>;
  auto f_init=[&]() -> Tobj {
    MyMatrix<Tint> M = GetPrincipalDomain<Tint>(eData.n);
    MyMatrix<Tint> CanM = LinPolytopeAntipodalIntegral_CanonicForm(M, os);
    return {CanM};
  };
  auto f_hash=[&](size_t const& seed, Tobj const& x) -> size_t {
    return Matrix_Hash(x.eMat, seed);
  };
  auto f_repr=[&](Tobj const& x, TadjI const& y, int const& i_rank, int const& i_orb) -> std::optional<TadjO> {
    if (x.eMat != y.ctype_arr.eMat) {
      return {};
    }
    TadjO ret{i_rank, i_orb};
    return ret;
  };
  auto f_spann=[&](TadjI const& x, int i_rank, int i_orb) -> std::pair<Tobj, TadjO> {
    Tobj IsoDel = x.DT_gram;
    MyMatrix<Tint> eBigMat = IdentityMat<Tint>(eData.LinSpa.n);
    TadjO ret{i_rank, i_orb, x.V, eBigMat};
    return {IsoDel, ret};
  };
  std::vector<IsoDelaunayDomain_MPI_Entry<T,Tint,Tgroup>> l_obj;
  std::vector<uint8_t> l_status;
  auto f_adj=[&](Tobj const& x, int i_orb) -> std::vector<TadjI> {

    int nb_free = CTYP_GetNumberFreeVectors(x);
    
    std::vector<FullAdjInfo<T>> ListIneq = ComputeDefiningIneqIsoDelaunayDomain<T,Tint,Tgroup>(x.DT, eData.LinSpa.ListLineMat);
    l_obj[i_orb].ListIneq = ListIneq;
    l_obj[i_orb].nb_free = nb_free;
    MyMatrix<T> FAC = GetFACineq(ListIneq);
    std::vector<int> ListIrred = cdd::RedundancyReductionClarkson(FAC);
    std::vector<TadjI> l_adj;
    for (auto & idxIrred : ListIrred) {
      FullAdjInfo<T> eRecIneq = ListIneq[idxIrred];
      DelaunayTesselation<Tint, Tgroup> DTadj = FlippingLtype<T,Tint,Tgroup>(x.DT, x.GramMat, eRecIneq.ListAdjInfo, eData.rddo);
      MyMatrix<T> GramMatAdj = GetInteriorGramMatrix(eData.LinSpa, DTadj);
      Tobj pair{DTadj, GramMatAdj};
      TadjI eAdj{eRecIneq.eIneq, pair};
      l_adj.push_back(eAdj);
    }
    return l_adj;
  };
  auto f_set_adj=[&](int const& i_orb, std::vector<TadjO> const& ListAdj) -> void {
    l_obj[i_orb].ListAdj = ListAdj;
  };
  auto f_obj=[&](TadjI const& x) -> Tobj {
    return x.DT_gram;
  };
  auto f_next=[&]() -> std::optional<std::pair<bool, Tobj>> {
    return {};
  };
  auto f_insert=[&](Tobj const& x) -> bool {
    l_obj.push_back({x, {}, {} });
    return false;
  };
  auto f_save_status=[&](size_t const& pos, bool const& val) -> void {
    uint8_t val_i = static_cast<uint8_t>(val);
    if (l_status.size() <= pos) {
      l_status.push_back(val_i);
    } else {
      l_status[pos] = val_i;
    }
  };
  bool test = compute_adjacency_mpi<Tobj,TadjI,TadjO,
    decltype(f_next),decltype(f_insert),decltype(f_obj),
    decltype(f_save_status),
    decltype(f_init),decltype(f_adj),decltype(f_set_adj),
    decltype(f_hash),decltype(f_repr),decltype(f_spann)>
    (comm, eData.max_runtime_second,
     f_next, f_insert, f_obj,
     f_save_status,
     f_init, f_adj, f_set_adj,
     f_hash, f_repr, f_spann, os);
  os << "Termination test=" << test << "\n";
  return l_obj;
}

template<typename T, typename Tint, typename Tgroup>
void WriteFamilyIsoDelaunayDomain(boost::mpi::communicator &comm, std::string const& OutFormat, std::string const& OutFile, std::vector<IsoDelaunayDomain_MPI_Entry<T,Tint,Tgroup>> const& ListIDD, [[maybe_unused]] std::ostream & os) {
  //  int i_rank = comm.rank();
  if (OutFormat == "nothing") {
    std::cerr << "No output\n";
    return;
  }
  if (OutFormat == "GAP") {
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
  throw TerminalException{1};
}


template<typename T, typename Tint, typename Tgroup>
void ComputeLatticeIsoDelaunayDomains(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  int i_rank = comm.rank();
  int n_proc = comm.size();
  std::string FileLog = "log_" + std::to_string(n_proc) + "_" + std::to_string(i_rank);
  std::ofstream os(FileLog);
  if (ApplyStdUnitbuf(eFull)) {
    os << std::unitbuf;
    os << "Apply UnitBuf\n";
  } else {
    os << "Do not apply UnitBuf\n";
  }
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  SingleBlock BlockTSPACE = eFull.ListBlock.at("TSPACE");
  //
  bool DATA_Saving = BlockDATA.ListBoolValues.at("Saving");
  std::string DATA_Prefix = BlockDATA.ListStringValues.at("Prefix");
  CreateDirectory(DATA_Prefix);
  //
  int max_runtime_second = BlockDATA.ListIntValues.at("max_runtime_second");
  std::cerr << "max_runtime_second=" <<	max_runtime_second << "\n";
  std::string OutFormat = BlockDATA.ListStringValues.at("OutFormat");
  std::string OutFile = BlockDATA.ListStringValues.at("OutFile");
  std::cerr << "OutFormat=" << OutFormat << " OutFile=" << OutFile << "\n";
  auto get_common=[&]() -> std::optional<MyMatrix<T>> {
    std::string CommonGramMat = BlockDATA.ListStringValues.at("CommonGramMat");
    if (CommonGramMat == "unset") {
      return {};
    }
    MyMatrix<T> eMat = ReadMatrixFile<T>(CommonGramMat);
    return eMat;
  };
  std::optional<MyMatrix<T>> CommonGramMat = get_common();
  //
  using TintGroup = typename Tgroup::Tint;
  PolyHeuristicSerial<TintGroup> AllArr = AllStandardHeuristicSerial<TintGroup>(os);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);

  LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint>(BlockTSPACE, os);

  DataIsoDelaunayDomains<T,Tint,Tgroup> eData{LinSpa,
    std::move(rddo),
    CommonGramMat,
    max_runtime_second,
    DATA_Saving,
    DATA_Prefix};

  std::vector<IsoDelaunayDomain_MPI_Entry<T, Tint, Tgroup>> ListIDD = MPI_EnumerationIsoDelaunayDomains<T,Tint,Tgroup>(comm, eData, os);
  WriteFamilyIsoDelaunayDomain(comm, OutFormat, OutFile, ListIDD, os);
}




// clang-format off
#endif  // SRC_LATT_ISODELAUNAYDOMAINS_H_
// clang-format on

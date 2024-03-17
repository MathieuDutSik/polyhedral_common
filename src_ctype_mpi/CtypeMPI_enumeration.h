// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_CTYPEMPI_ENUMERATION_H_
#define SRC_CTYPEMPI_ENUMERATION_H_

// clang-format off
#include "NumberTheory.h"
#include "CtypeMPI_types.h"
#include "POLY_AdjacencyScheme.h"
#include "Namelist.h"
#include "hash_functions.h"
#include "rational.h"
#include "sparse-map/include/tsl/sparse_map.h"
#include <boost/mpi.hpp>
#include <unordered_map>
// clang-format on

struct DataCtype {
  int n;
  int max_runtime_second;
  bool Saving;
  std::string Prefix;
  std::ostream& os;
};

FullNamelist NAMELIST_GetStandard_COMPUTE_LATTICE_IsoEdgeDomains() {
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
  ListIntValues1["n"] = -1;
  ListIntValues1["max_runtime_second"] = 0;
  ListBoolValues1["ApplyStdUnitbuf"] = false;
  ListBoolValues1["Saving"] = false;
  ListStringValues1["Prefix"] = "/irrelevant/";
  SingleBlock BlockDATA;
  BlockDATA.ListIntValues = ListIntValues1;
  BlockDATA.ListBoolValues = ListBoolValues1;
  BlockDATA.ListDoubleValues = ListDoubleValues1;
  BlockDATA.ListStringValues = ListStringValues1;
  BlockDATA.ListListStringValues = ListListStringValues1;
  ListBlock["DATA"] = BlockDATA;
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

template<typename Tint>
struct IsoEdgeDomain_AdjI {
  TypeCtypeExch<Tint> ctype_arr;
};

namespace boost::serialization {
  template <class Archive, typename Tint>
  inline void serialize(Archive &ar, IsoEdgeDomain_AdjI<Tint> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("ctype_arr", eRec.ctype_arr);
  }
}

struct IsoEdgeDomain_AdjO {
};

namespace boost::serialization {
  template <class Archive>
  inline void serialize([[maybe_unused]] Archive &ar, [[maybe_unused]] IsoEdgeDomain_AdjO &eRec,
                        [[maybe_unused]] const unsigned int version) {
  }
}

template<typename T, typename Tint, typename Tgroup>
struct IsoEdgeDomain_MPI_Entry {
  TypeCtypeExch<Tint> ctype_arr;
  StructuralInfo struct_info;
  std::vector<IsoEdgeDomain_MPI_AdjO<T, Tint>> ListAdj;
};

namespace boost::serialization {
  template <class Archive, typename T, typename Tint, typename Tgroup>
  inline void serialize(Archive &ar, IsoEdgeDomain_MPI_Entry<T, Tint, Tgroup> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("ctype_arr", eRec.ctype_arr);
    ar &make_nvp("struct_info", eRec.struct_info);
    ar &make_nvp("ListAdj", eRec.ListAdj);
  }
}

template<typename Tint>
struct IsoEdgeDomain_Obj {
  TypeCtypeExch<Tint> ctype_arr;
  StructuralInfo struct_info;
};

namespace boost::serialization {
  template <class Archive, typename Tint>
  inline void serialize(Archive &ar, IsoEdgeDomain_Obj<Tint> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("ctype_arr", eRec.ctype_arr);
    ar &make_nvp("struct_info", eRec.struct_info);
  }
}

template<typename T, typename Tint, typename Tgroup>
struct DataCtypeFunc {
  DataCtype data;
  using Tobj = IsoEdgeDomain_Obj<Tint>;
  using TadjI = IsoEdgeDomain_AdjI<Tint>;
  using TadjO = IsoEdgeDomain_AdjO;
  std::ostream& get_os() {
    return data.os;
  }
  Tobj f_init() {
    MyMatrix<Tint> M = GetPrincipalDomain<Tint>(data.n);
    MyMatrix<Tint> CanM = LinPolytopeAntipodalIntegral_CanonicForm(M, data.os);
    return {CanM};
  }
  size_t f_hash(size_t const& seed, Tobj const& x) {
    return Matrix_Hash(x.eMat, seed);
  };
  std::optional<TadjO> f_repr(Tobj const& x, TadjI const& y) {
    if (x.eMat != y.ctype_arr.eMat) {
      return {};
    }
    TadjO ret{};
    return ret;
  }
  std::pair<Tobj, TadjO> f_spann(TadjI const& x) {
    TypeCtypeExch<Tint> ctype_arr = x.ctype_arr;
    Tobj x_ret{ctype_arr, {} };
    TadjO ret{};
    return {x_ret, ret};
  }
  std::vector<TadjI> f_adj(Tobj & x_in) {
    using Tidx = typename Tgroup::Tidx;
    Tobj x = x_in.ctype_arr;
    int nb_free = CTYP_GetNumberFreeVectors(x);
    int nb_autom = CTYP_GetNbAutom<Tint, Tgroup>(x, data.os);
    DataCtypeFacet<Tint, Tidx> data = CTYP_GetConeInformation<Tint, Tidx>(x);
    int nb_triple = data.nb_triple;
    int nb_ineq = data.nb_ineq;
    int nb_ineq_after_crit = data.nb_ineq_after_crit;
    StructuralInfo struct_info{nb_triple, nb_ineq, nb_ineq_after_crit, nb_free, nb_autom};
    x_in.struct_info = struct_info;
    std::vector<TadjI> ListRet;
    for (auto &e_int : data.ListIrred) {
      MyMatrix<Tint> FlipMat = CTYP_TheFlipping(data.TheCtype, data.ListInformations[e_int]);
      MyMatrix<Tint> CanMat = LinPolytopeAntipodalIntegral_CanonicForm(FlipMat, data.os);
      TypeCtypeExch<Tint> x{std::move(CanMat)};
      TadjI x_adjI{x};
      ListRet.push_back(x_adjI);
    }
    return ListRet;
  }
  Tobj f_adji_obj=[&](TadjI const& x) {
    TypeCtypeExch<Tint> ctype_arr = x.ctype_arr;
    Tobj x_ret{ctype_arr, {} };
    return x_ret;
  };
};


template<typename T, typename Tint, typename Tgroup>
std::pair<bool, std::vector<IsoEdgeDomain_MPI_Entry<T,Tint,Tgroup>>> MPI_EnumerationIsoEdgeDomains(boost::mpi::communicator &comm, DataCtype const& eData, std::ostream & os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tover = IsoEdgeDomain_MPI_Entry<T,Tint,Tgroup>;
  using Tobj = TypeCtypeExch<Tint>;
  using TadjI = IsoEdgeDomain_AdjI<Tint>;
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
    Tobj IsoDel = x.ctype_arr;
    TadjO ret{i_rank, i_orb};
    return {IsoDel, ret};
  };
  std::vector<Tover> l_obj;
  std::vector<uint8_t> l_status;
  int i_rank = comm.rank();
  int n_proc = comm.size();
  std::string str_proc = "_nproc" + std::to_string(n_proc) + "_rank" + std::to_string(i_rank);
  PartialEnum_FullRead(eData.Prefix, str_proc, eData.Saving, l_obj, l_status, os);
  auto f=[](Tover const& x) -> Tobj {
    return x.ctype_arr;
  };
  //  auto next_iterator(l_obj, l_status, f);
  NextIterator<Tover,Tobj,decltype(f)> next_iterator(l_obj, l_status, f);
  auto f_next=[&]() -> std::optional<std::pair<bool, Tobj>> {
    return next_iterator.f_next();
  };
  auto f_adj=[&](int i_orb) -> std::vector<TadjI> {
    Tobj x = l_obj[i_orb].ctype_arr;
    int nb_free = CTYP_GetNumberFreeVectors(x);
    int nb_autom = CTYP_GetNbAutom<Tint, Tgroup>(x, os);
    DataCtypeFacet<Tint, Tidx> data = CTYP_GetConeInformation<Tint, Tidx>(x);
    int nb_triple = data.nb_triple;
    int nb_ineq = data.nb_ineq;
    int nb_ineq_after_crit = data.nb_ineq_after_crit;
    StructuralInfo struct_info{nb_triple, nb_ineq, nb_ineq_after_crit, nb_free, nb_autom};
    l_obj[i_orb].struct_info = struct_info;
    std::vector<TadjI> ListRet;
    for (auto &e_int : data.ListIrred) {
      MyMatrix<Tint> FlipMat = CTYP_TheFlipping(data.TheCtype, data.ListInformations[e_int]);
      MyMatrix<Tint> CanMat = LinPolytopeAntipodalIntegral_CanonicForm(FlipMat, os);
      TypeCtypeExch<Tint> x{std::move(CanMat)};
      TadjI x_adjI{x};
      ListRet.push_back(x_adjI);
    }
    return ListRet;
  };
  auto f_set_adj=[&](int const& i_orb, std::vector<TadjO> const& ListAdj) -> void {
    l_obj[i_orb].ListAdj = ListAdj;
  };
  auto f_adji_obj=[&](TadjI const& x) -> Tobj {
    return x.ctype_arr;
  };
  auto f_idx_obj=[&](size_t const& idx) -> Tobj {
    return l_obj[idx].ctype_arr;
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
    decltype(f_next),decltype(f_insert),decltype(f_adji_obj),
    decltype(f_idx_obj), decltype(f_save_status),
    decltype(f_init),decltype(f_adj),decltype(f_set_adj),
    decltype(f_hash),decltype(f_repr),decltype(f_spann)>
    (comm, eData.max_runtime_second,
     f_next, f_insert, f_adji_obj,
     f_idx_obj, f_save_status,
     f_init, f_adj, f_set_adj,
     f_hash, f_repr, f_spann, os);
  os << "Termination test=" << test << "\n";
  return {test, std::move(l_obj)};
}

template<typename T, typename Tint, typename Tgroup>
void WriteFamilyIsoEdgeDomain([[maybe_unused]] boost::mpi::communicator &comm, std::string const& OutFormat, [[maybe_unused]] std::string const& OutFile, [[maybe_unused]] std::vector<IsoEdgeDomain_MPI_Entry<T,Tint,Tgroup>> const& ListIDD, [[maybe_unused]] std::ostream & os) {
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
void ComputeLatticeIsoEdgeDomains(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  bool ApplyStdUnitbuf = BlockDATA.ListBoolValues.at("ApplyStdUnitbuf");
  int i_rank = comm.rank();
  int n_proc = comm.size();
  std::string FileLog = "log_" + std::to_string(n_proc) + "_" + std::to_string(i_rank);
  std::ofstream os(FileLog);
  if (ApplyStdUnitbuf) {
    os << std::unitbuf;
    os << "Apply UnitBuf\n";
  } else {
    os << "Do not apply UnitBuf\n";
  }
  //
  bool STORAGE_Saving = BlockDATA.ListBoolValues.at("Saving");
  std::string STORAGE_Prefix = BlockDATA.ListStringValues.at("Prefix");
  CreateDirectory(STORAGE_Prefix);
  //
  int n = BlockDATA.ListIntValues.at("n");
  int max_runtime_second = BlockDATA.ListIntValues.at("max_runtime_second");
  std::cerr << "max_runtime_second=" <<	max_runtime_second << "\n";
  std::string OutFormat = BlockDATA.ListStringValues.at("OutFormat");
  std::string OutFile = BlockDATA.ListStringValues.at("OutFile");
  std::cerr << "OutFormat=" << OutFormat << " OutFile=" << OutFile << "\n";
  //
  using TintGroup = typename Tgroup::Tint;

  DataCtype eData{n,
    max_runtime_second,
    STORAGE_Saving,
    STORAGE_Prefix, os};

  std::pair<bool, std::vector<IsoEdgeDomain_MPI_Entry<T, Tint, Tgroup>>> pair = MPI_EnumerationIsoEdgeDomains<T,Tint,Tgroup>(comm, eData, os);
  if (pair.first) {
    WriteFamilyIsoEdgeDomain(comm, OutFormat, OutFile, pair.second, os);
  }
}

// clang-format off
#endif  // SRC_CTYPEMPI_ENUMERATION_H_
// clang-format on

// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_CTYPE_MPI_CTYPEMPI_ENUMERATION_H_
#define SRC_CTYPE_MPI_CTYPEMPI_ENUMERATION_H_

// clang-format off
#include "NumberTheory.h"
#include "CtypeMPI_types.h"
#include "POLY_MPI_AdjacencyScheme.h"
#include "Namelist.h"
#include "hash_functions.h"
#include "rational.h"
#include "sparse-map/include/tsl/sparse_map.h"
#include <boost/mpi.hpp>
#include <unordered_map>
#include <map>
#include <string>
#include <utility>
#include <vector>
// clang-format on

struct DataCtype {
  int n;
  std::ostream &os;
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

template <typename Tint> struct IsoEdgeDomain_AdjI {
  TypeCtypeExch<Tint> ctype_arr;
};

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, IsoEdgeDomain_AdjI<Tint> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("ctype_arr", eRec.ctype_arr);
}
} // namespace boost::serialization

struct IsoEdgeDomain_AdjO {};

void WriteEntryGAP(std::ostream &os_out,
                   [[maybe_unused]] IsoEdgeDomain_AdjO const &_x) {
  os_out << "[]";
}

void WriteEntryPYTHON(std::ostream &os_out,
                   [[maybe_unused]] IsoEdgeDomain_AdjO const &_x) {
  os_out << "{}";
}

namespace boost::serialization {
template <class Archive>
inline void serialize([[maybe_unused]] Archive &ar,
                      [[maybe_unused]] IsoEdgeDomain_AdjO &eRec,
                      [[maybe_unused]] const unsigned int version) {}
} // namespace boost::serialization

template <typename Tint> struct IsoEdgeDomain_Obj {
  TypeCtypeExch<Tint> ctype_arr;
  StructuralInfo struct_info;
};

template <typename Tint>
void WriteEntryGAP(std::ostream &os_out, IsoEdgeDomain_Obj<Tint> const &entry) {
  os_out << "rec(Ctype:=";
  WriteMatrixGAP(os_out, entry.ctype_arr.eMat);
  os_out << ", struct_info:=rec(";
  os_out << "nb_triple:=" << entry.struct_info.nb_triple;
  os_out << ", nb_ineq:=" << entry.struct_info.nb_ineq;
  os_out << ", nb_ineq_after_crit:=" << entry.struct_info.nb_ineq_after_crit;
  os_out << ", nb_free:=" << entry.struct_info.nb_free;
  os_out << ", nb_autom:=" << entry.struct_info.nb_autom;
  os_out << "))";
}

template <typename Tint>
void WriteEntryPYTHON(std::ostream &os_out, IsoEdgeDomain_Obj<Tint> const &entry) {
  os_out << "{\"Ctype\":";
  WriteMatrixPYTHON(os_out, entry.ctype_arr.eMat);
  os_out << ", \"struct_info\":={";
  os_out << "\"nb_triple\":" << entry.struct_info.nb_triple;
  os_out << ", \"nb_ineq\":" << entry.struct_info.nb_ineq;
  os_out << ", \"nb_ineq_after_crit\":" << entry.struct_info.nb_ineq_after_crit;
  os_out << ", \"nb_free\":" << entry.struct_info.nb_free;
  os_out << ", \"nb_autom\":" << entry.struct_info.nb_autom;
  os_out << "}}";
}

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, IsoEdgeDomain_Obj<Tint> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("ctype_arr", eRec.ctype_arr);
  ar &make_nvp("struct_info", eRec.struct_info);
}
} // namespace boost::serialization

template <typename T, typename Tint, typename Tgroup> struct DataCtypeFunc {
  DataCtype data;
  using Tobj = IsoEdgeDomain_Obj<Tint>;
  using TadjI = IsoEdgeDomain_AdjI<Tint>;
  using TadjO = IsoEdgeDomain_AdjO;
  std::ostream &get_os() { return data.os; }
  Tobj f_init() {
    MyMatrix<Tint> M = GetPrincipalDomain<Tint>(data.n);
    MyMatrix<Tint> CanM = LinPolytopeAntipodalIntegral_CanonicForm(M, data.os);
    TypeCtypeExch<Tint> ctype_arr{CanM};
    Tobj x{ctype_arr, {}};
    return x;
  }
  size_t f_hash(size_t const &seed, Tobj const &x) {
    return Matrix_Hash(x.ctype_arr.eMat, seed);
  }
  std::optional<TadjO> f_repr(Tobj const &x, TadjI const &y) {
    if (x.ctype_arr.eMat != y.ctype_arr.eMat) {
      return {};
    }
    TadjO ret{};
    return ret;
  }
  std::pair<Tobj, TadjO> f_spann(TadjI const &x) {
    TypeCtypeExch<Tint> ctype_arr = x.ctype_arr;
    Tobj x_ret{ctype_arr, {}};
    TadjO ret{};
    return {x_ret, ret};
  }
  std::vector<TadjI> f_adj(Tobj &x_in) {
    using Tidx = typename Tgroup::Telt::Tidx;
    TypeCtypeExch<Tint> x = x_in.ctype_arr;
    int nb_free = CTYP_GetNumberFreeVectors(x);
    int nb_autom = CTYP_GetNbAutom<Tint, Tgroup>(x, data.os);
    DataCtypeFacet<Tint, Tidx> data_facet =
        CTYP_GetConeInformation<Tint, Tidx>(x, data.os);
    int nb_triple = data_facet.nb_triple;
    int nb_ineq = data_facet.nb_ineq;
    int nb_ineq_after_crit = data_facet.nb_ineq_after_crit;
    StructuralInfo struct_info{nb_triple, nb_ineq, nb_ineq_after_crit, nb_free,
                               nb_autom};
    x_in.struct_info = struct_info;
    std::vector<TadjI> ListRet;
    for (auto &e_int : data_facet.ListIrred) {
      MyMatrix<Tint> FlipMat = CTYP_TheFlipping(
          data_facet.TheCtype, data_facet.ListInformations[e_int]);
      MyMatrix<Tint> CanMat =
          LinPolytopeAntipodalIntegral_CanonicForm(FlipMat, data.os);
      TypeCtypeExch<Tint> x{std::move(CanMat)};
      TadjI x_adjI{x};
      ListRet.push_back(x_adjI);
    }
    return ListRet;
  }
  Tobj f_adji_obj(TadjI const &x) {
    TypeCtypeExch<Tint> ctype_arr = x.ctype_arr;
    Tobj x_ret{ctype_arr, {}};
    return x_ret;
  }
};

template <typename T, typename Tint, typename Tgroup>
void ComputeLatticeIsoEdgeDomains(boost::mpi::communicator &comm,
                                  FullNamelist const &eFull) {
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  bool ApplyStdUnitbuf = BlockDATA.ListBoolValues.at("ApplyStdUnitbuf");
  int i_rank = comm.rank();
  int n_proc = comm.size();
  std::string FileLog =
      "log_" + std::to_string(n_proc) + "_" + std::to_string(i_rank);
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
  std::cerr << "max_runtime_second=" << max_runtime_second << "\n";
  std::string OutFormat = BlockDATA.ListStringValues.at("OutFormat");
  std::string OutFile = BlockDATA.ListStringValues.at("OutFile");
  std::cerr << "OutFormat=" << OutFormat << " OutFile=" << OutFile << "\n";

  DataCtype data{n, os};
  using Tdata = DataCtypeFunc<T, Tint, Tgroup>;
  Tdata data_fct{data};
  using Tobj = typename DataCtypeFunc<T, Tint, Tgroup>::Tobj;
  using TadjO = typename DataCtypeFunc<T, Tint, Tgroup>::TadjO;
  using Tout = DatabaseEntry_MPI<Tobj, TadjO>;
  std::pair<bool, std::vector<Tout>> pair = EnumerateAndStore_MPI<Tdata>(
      comm, data_fct, STORAGE_Prefix, STORAGE_Saving, max_runtime_second);
  if (pair.first) {
    std::ofstream os_out(OutFile);
    bool result = WriteFamilyObjects_MPI<Tobj, TadjO>(comm, OutFormat, os_out, pair.second, os);
    if (result) {
      std::cerr << "CTYP_MPI: Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
      throw TerminalException{1};
    }
  }
}

// clang-format off
#endif  // SRC_CTYPE_MPI_CTYPEMPI_ENUMERATION_H_
// clang-format on

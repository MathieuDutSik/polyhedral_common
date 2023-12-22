// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEDELAUNAY_H_
#define SRC_LATT_LATTICEDELAUNAY_H_

// clang-format off
#include "boost_serialization.h"
#include "FundamentalDelaunay.h"
#include "GRP_DoubleCoset.h"
#include "MatrixGroup.h"
#include "POLY_RecursiveDualDesc.h"
#include "POLY_AdjacencyScheme.h"
#include <map>
#include <string>
#include <vector>
// clang-format on

#ifdef TIMINGS
#define TIMINGS_DELAUNAY_ENUMERATION
#endif

template <typename T, typename Tint> struct DataLattice {
  int n;
  MyMatrix<T> GramMat;
  MyMatrix<T> SHV;
  std::string CVPmethod;
  int max_runtime_second;
  bool Saving;
  std::string Prefix;
};

template <typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value>
GetWeightMatrixFromGramEXT(MyMatrix<T> const &EXT, MyMatrix<T> const &GramMat,
                           MyMatrix<T> const &SHV, std::ostream& os) {
  int n = GramMat.rows();
  CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral(GramMat, EXT);
  MyVector<T> TheCenter = eCP.eCent;
  std::vector<MyMatrix<T>> CharPair = CharacterizingPair(GramMat, TheCenter);
  MyMatrix<T> Qmat = CharPair[0];
  int nbVect = SHV.rows();
  int nbVert = EXT.rows();
  MyMatrix<T> EXText(nbVect + nbVert, n + 1);
  for (int iVert = 0; iVert < nbVert; iVert++)
    EXText.row(iVert) = EXT.row(iVert);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    EXText(nbVert + iVect, 0) = 0;
    for (int i = 0; i < n; i++)
      EXText(nbVert + iVect, i + 1) = SHV(iVect, i);
  }
  return GetSimpleWeightMatrix<T, Tidx_value>(EXText, Qmat, os);
}

template <typename T, typename Tgroup>
bool IsGroupCorrect(MyMatrix<T> const &EXT_T, Tgroup const &eGRP) {
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> LGen = eGRP.GeneratorsOfGroup();
  for (auto &eGen : LGen) {
    MyMatrix<T> eMat = FindTransformation(EXT_T, EXT_T, eGen);
    if (!IsIntegralMatrix(eMat))
      return false;
  }
  return true;
}

template <typename T, typename Tint, typename Tgroup>
Tgroup Delaunay_Stabilizer(DataLattice<T, Tint> const &eData,
                           MyMatrix<Tint> const &EXT, std::ostream& os) {
  using Tidx_value = int16_t;
  using Tgr = GraphListAdj;
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(EXT);
  //
  // Now extending with the SHV vector set
  //
  WeightMatrix<true, T, Tidx_value> WMat =
    GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT_T, eData.GramMat, eData.SHV, os);
  int nbVert = EXT_T.rows();
  int nbSHV = eData.SHV.rows();
  Face eFace(nbVert + nbSHV);
  for (int iVert = 0; iVert < nbVert; iVert++)
    eFace[iVert] = 1;
  for (int iSHV = 0; iSHV < nbSHV; iSHV++)
    eFace[nbVert + iSHV] = 0;
  Tgroup PreGRPisom =
    GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat, os);
  Tgroup GRPisom = ReducedGroupAction(PreGRPisom, eFace);
  return LinPolytopeIntegral_Stabilizer_Method8(EXT_T, GRPisom, os);
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
Delaunay_TestEquivalence(DataLattice<T, Tint> const &eData,
                         MyMatrix<Tint> const &EXT1,
                         MyMatrix<Tint> const &EXT2,
                         std::ostream & os) {
  using Telt = typename Tgroup::Telt;
  using Tgr = GraphListAdj;
  using Tidx_value = int16_t;
  os << "Begin Delaunay_TestEquivalence\n";
  MyMatrix<T> EXT1_T = UniversalMatrixConversion<T, Tint>(EXT1);
  MyMatrix<T> EXT2_T = UniversalMatrixConversion<T, Tint>(EXT2);
  //
  // Now extending by adding more vectors.
  //
  WeightMatrix<true, T, Tidx_value> WMat1 =
    GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT1_T, eData.GramMat, eData.SHV, os);
  WeightMatrix<true, T, Tidx_value> WMat2 =
    GetWeightMatrixFromGramEXT<T, Tidx_value>(EXT2_T, eData.GramMat, eData.SHV, os);
  std::optional<Telt> eRes =
    TestEquivalenceWeightMatrix<T, Telt, Tidx_value>(WMat1, WMat2, os);
  if (!eRes) {
    os << "Leaving Delaunay_TestEquivalence with false\n";
    return {};
  }
  Telt const& eElt = *eRes;
  MyMatrix<T> MatEquiv_T = FindTransformation(EXT1_T, EXT2_T, eElt);
  if (IsIntegralMatrix(MatEquiv_T)) {
    MyMatrix<Tint> MatEquiv_I = UniversalMatrixConversion<Tint, T>(MatEquiv_T);
    os << "Leaving Delaunay_TestEquivalence with true\n";
    return MatEquiv_I;
  }
  os << "Trying other strategies\n";
  Tgroup GRP1 = GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat1, os);
  std::optional<MyMatrix<T>> eEq = LinPolytopeIntegral_Isomorphism_Method8(EXT1_T, EXT2_T, GRP1, eElt, os);
  if (!eEq)
    return {};
  MyMatrix<Tint> eMat_I = UniversalMatrixConversion<Tint, T>(*eEq);
  return eMat_I;
}

template <typename T, typename Tint>
size_t ComputeInvariantDelaunay(DataLattice<T, Tint> const &eData,
                                size_t const& seed,
                                MyMatrix<Tint> const& EXT) {
  int nbVert = EXT.rows();
  int n = EXT.cols() - 1;
  Tint PreIndex = Int_IndexLattice(EXT);
  Tint eIndex = T_abs(PreIndex);
  MyVector<T> V(n);
  MyMatrix<T> EXT_fullT = UniversalMatrixConversion<T,Tint>(EXT);
  CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral(eData.GramMat, EXT_fullT);
  MyMatrix<T> EXT_T(nbVert, n);
  for (int iVert=0; iVert<nbVert; iVert++) {
    for (int i=0; i<n; i++) {
      T val = UniversalScalarConversion<T,Tint>(EXT(iVert, i+1));
      EXT_T(iVert, i) = val - eCP.eCent(i+1);
    }
  }
  std::map<T, size_t> ListDiagNorm;
  std::map<T, size_t> ListOffDiagNorm;
  for (int iVert = 0; iVert < nbVert; iVert++) {
    MyVector<T> V(n);
    for (int i=0; i<n; i++) {
      T eSum = 0;
      for (int j=0; j<n; j++) {
        eSum += eData.GramMat(i,j) * EXT_T(iVert, j);
      }
      V(i) = eSum;
    }
    T scal = 0;
    for (int i=0; i<n; i++) {
      scal += V(i) * EXT_T(iVert, i);
    }
    ListDiagNorm[scal] += 1;
    for (int jVert=iVert+1; jVert<nbVert; jVert++) {
      T scal = 0;
      for (int i=0; i<n; i++) {
        scal += V(i) * EXT_T(jVert, i);
      }
      ListOffDiagNorm[scal] += 1;
    }
  }
  auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
    seed ^= new_hash + 0x9e3779b8 + (seed << 6) + (seed >> 2);
  };
  size_t hash = seed;
  auto update_from_map=[&](std::map<T, size_t> const& map) -> void {
    for (auto & kv : map) {
      size_t hash1 = std::hash<T>()(kv.first);
      size_t hash2 = kv.second;
      combine_hash(hash, hash1);
      combine_hash(hash, hash2);
    }
  };
  update_from_map(ListDiagNorm);
  update_from_map(ListOffDiagNorm);
  return hash;
}

template<typename Tint>
struct Delaunay_AdjI {
  Face f;
  MyMatrix<Tint> obj;
};

namespace boost::serialization {
  template <class Archive, typename Tint>
  inline void serialize(Archive &ar, Delaunay_AdjI<Tint> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("f", eRec.f);
    ar &make_nvp("obj", eRec.obj);
  }
}

template<typename Tint>
struct Delaunay_MPI_AdjO {
  Face f;
  MyMatrix<Tint> P;
  int iProc;
  int iOrb;
};

namespace boost::serialization {
  template <class Archive, typename Tint>
  inline void serialize(Archive &ar, Delaunay_MPI_AdjO<Tint> &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("f", eRec.f);
    ar &make_nvp("P", eRec.P);
    ar &make_nvp("iProc", eRec.iProc);
    ar &make_nvp("iOrb", eRec.iOrb);
  }
}

template<typename Tint, typename Tgroup>
struct Delaunay_MPI_Entry {
  MyMatrix<Tint> obj;
  Tgroup GRP;
  std::vector<Delaunay_MPI_AdjO<Tint>> l_adj;
};


template <typename T, typename Tint, typename Tgroup>
std::vector<Delaunay_MPI_Entry<Tint, Tgroup>> EnumerationDelaunayPolytopes(boost::mpi::communicator &comm,
                                                                           std::ostream & os,
                                                                           DataLattice<T, Tint> const &eData,
                                                                           PolyHeuristic<mpz_class> const &AllArr) {
  using Tobj = MyMatrix<Tint>;
  using TadjI = Delaunay_AdjI<Tint>;
  using TadjO = Delaunay_MPI_AdjO<Tint>;
  auto f_init=[&]() -> Tobj {
    Tobj EXT = FindDelaunayPolytope<T, Tint>(
       eData.GramMat, eData.CVPmethod, os);
    os << "Creation of a Delaunay with |V|=" << EXT.rows() << " vertices\n";
    return EXT;
  };
  auto f_hash=[&](size_t const& seed, Tobj const& x) -> size_t {
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    MicrosecondTime time;
#endif
    size_t hash = ComputeInvariantDelaunay(eData, seed, x);
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    std::cerr << "|ComputeInvariantDelaunay|=" << time << "\n";
#endif
    return hash;
  };
  auto f_repr=[&](Tobj const& x, TadjI const& y, int const& i_rank, int const& i_orb) -> std::optional<TadjO> {
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    MicrosecondTime time;
#endif
    std::optional<MyMatrix<Tint>> opt = Delaunay_TestEquivalence<T, Tint, Tgroup>(eData, x, y.obj, os);
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    std::cerr << "|Delaunay_TestEquivalence|=" << time << "\n";
#endif
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const& P = *opt;
    TadjO ret{y.f, P, i_rank, i_orb};
    return ret;
  };
  auto f_spann=[&](TadjI const& x, int i_rank, int i_orb) -> std::pair<Tobj, TadjO> {
    Tobj EXT = x.obj;
    MyMatrix<Tint> P = IdentityMat<Tint>(eData.n);
    TadjO ret{x.f, P, i_rank, i_orb};
    return {EXT, ret};
  };
  std::vector<Delaunay_MPI_Entry<Tint,Tgroup>> l_obj;
  std::vector<int> l_status;
  auto f_adj=[&](Tobj const& x, int i_orb) -> std::vector<TadjI> {
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    MicrosecondTime time;
#endif
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T,Tint>(x);
    Tgroup GRPlatt = Delaunay_Stabilizer<T, Tint, Tgroup>(eData, x, os);
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    std::cerr << "|Delaunay_Stabilizer|=" << time << "\n";
#endif
    std::cerr << "i_orb=" << i_orb << " |l_obj|=" << l_obj.size() << "\n";
    l_obj[i_orb].GRP = GRPlatt;
    vectface TheOutput = DualDescriptionStandard(EXT_T, GRPlatt);
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    std::cerr << "|DualDescriptionStandard|=" << time << "\n";
#endif
    std::vector<TadjI> l_adj;
#ifdef TIMINGS_DELAUNAY_ENUMERATION
    std::cerr << "|l_adj|=" << l_adj.size() << "\n";
#endif
    for (auto &eOrbB : TheOutput) {
      MyMatrix<Tint> EXTadj = FindAdjacentDelaunayPolytope<T, Tint>(eData.GramMat, EXT_T, eOrbB, eData.CVPmethod);
#ifdef TIMINGS_DELAUNAY_ENUMERATION
      std::cerr << "|FindAdjacentDelaunayPolytope|=" << time << "\n";
#endif
      TadjI eAdj{eOrbB, EXTadj};
      l_adj.push_back(eAdj);
    }
    return l_adj;
  };
  auto f_set_adj=[&](int const& i_orb, std::vector<TadjO> const& l_adj) -> void {
    l_obj[i_orb].l_adj = l_adj;
  };
  auto f_exists=[&](int const& n_obj) -> bool {
    return false;
  };
  auto f_insert=[&](Tobj const& x) -> bool {
    Tgroup grp;
    l_obj.push_back({x, grp, {} });
    return false;
  };
  auto f_load=[&](size_t const& pos) -> Tobj {
    return l_obj[pos].obj;
  };
  auto f_save_status=[&](size_t const& pos, bool const& val) -> void {
    int val_i = static_cast<int>(val);
    if (l_status.size() <= pos) {
      l_status.push_back(val_i);
    } else {
      l_status[pos] = val_i;
    }
  };
  auto f_load_status=[&](size_t const& pos) -> bool {
    return static_cast<bool>(l_status[pos]);
  };
  compute_adjacency_mpi<Tobj,TadjI,TadjO,
    decltype(f_exists),decltype(f_insert),decltype(f_load),
    decltype(f_save_status),decltype(f_load_status),
    decltype(f_init),decltype(f_adj),decltype(f_set_adj),
    decltype(f_hash),decltype(f_repr),decltype(f_spann)>
    (comm, os, eData.max_runtime_second,
                        f_exists, f_insert, f_load,
                        f_save_status, f_load_status,
                        f_init, f_adj, f_set_adj,
                        f_hash, f_repr, f_spann);
  return l_obj;
}

FullNamelist NAMELIST_GetStandard_COMPUTE_DELAUNAY() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["arithmetic_T"] = "gmp_rational";
  ListStringValues1["arithmetic_Tint"] = "gmp_integer";
  ListStringValues1["GRAMfile"] = "unset.gram";
  ListStringValues1["SVRfile"] = "unset.svr";
  ListStringValues1["OutFormat"] = "nothing";
  ListStringValues1["OutFile"] = "unset.out";
  ListStringValues1["FileDualDescription"] = "unset";
  ListIntValues1["max_runtime_second"] = 0;
  SingleBlock BlockDATA;
  BlockDATA.ListIntValues = ListIntValues1;
  BlockDATA.ListBoolValues = ListBoolValues1;
  BlockDATA.ListDoubleValues = ListDoubleValues1;
  BlockDATA.ListStringValues = ListStringValues1;
  BlockDATA.ListListStringValues = ListListStringValues1;
  ListBlock["DATA"] = BlockDATA;
  // STORAGE
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListBoolValues2["Saving"] = false;
  ListStringValues2["Prefix"] = "/irrelevant/";
  ListStringValues2["CVPmethod"] = "SVexact";
  SingleBlock BlockSTORAGE;
  BlockSTORAGE.ListIntValues = ListIntValues2;
  BlockSTORAGE.ListBoolValues = ListBoolValues2;
  BlockSTORAGE.ListDoubleValues = ListDoubleValues2;
  BlockSTORAGE.ListStringValues = ListStringValues2;
  BlockSTORAGE.ListListStringValues = ListListStringValues2;
  ListBlock["STORAGE"] = BlockSTORAGE;
  // Merging all data
  return {ListBlock, "undefined"};
}

template<typename Tint, typename Tgroup>
void WriteFamilyDelaunay(std::string const& OutFormat, std::string const& OutFile, std::vector<Delaunay_MPI_Entry<Tint, Tgroup>> const& ListDel) {
  if (OutFormat == "nothing") {
    std::cerr << "No output\n";
    return;
  }
  if (OutFormat == "GAPformat") {
    std::ofstream OUTfs(OutFile);
    int nbDel = ListDel.size();
    OUTfs << "nbDel=" << nbDel << "\n";
    for (int iDel = 0; iDel < nbDel; iDel++) {
      OUTfs << "iDel=" << iDel << "/" << nbDel << "\n";
      WriteMatrix(OUTfs, ListDel[iDel].obj);
    }
  }
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
  throw TerminalException{1};
}


template<typename T, typename Tint, typename Tgroup>
void ComputeDelaunayPolytope(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  int i_rank = comm.rank();
  int n_proc = comm.size();
  std::string FileLog = "log_" + std::to_string(n_proc) + "_" + std::to_string(i_rank);
  std::ofstream os(FileLog);
  SingleBlock BlockSTORAGE = eFull.ListBlock.at("STORAGE");
  //
  bool STORAGE_Saving = BlockSTORAGE.ListBoolValues.at("Saving");
  std::string STORAGE_Prefix = BlockSTORAGE.ListStringValues.at("Prefix");
  CreateDirectory(STORAGE_Prefix);
  //
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  int max_runtime_second = BlockDATA.ListIntValues.at("max_runtime_second");
  std::cerr << "Reading DATA\n";
  std::string GRAMfile = BlockDATA.ListStringValues.at("GRAMfile");
  MyMatrix<T> GramMat = ReadMatrixFile<T>(GRAMfile);
  //
  std::string SVRfile = BlockDATA.ListStringValues.at("SVRfile");
  auto get_SVR=[&]() -> MyMatrix<T> {
    if (IsExistingFile(SVRfile)) {
      return ReadMatrixFile<T>(SVRfile);
    }
    int n = GramMat.rows();
    return ZeroMatrix<T>(0, n);
  };
  MyMatrix<T> SVR = get_SVR();
  std::cerr << "|SVR|=" << SVR.rows() << "\n";
  //
  std::string OutFormat = BlockDATA.ListStringValues.at("OutFormat");
  std::string OutFile = BlockDATA.ListStringValues.at("OutFile");
  std::cerr << "OutFile=" << OutFile << "\n";

  int n = GramMat.rows();
  std::string CVPmethod = "SVexact";
  DataLattice<T, Tint> eData{n,
                             GramMat,
                             SVR,
                             CVPmethod,
                             max_runtime_second,
                             STORAGE_Saving,
                             STORAGE_Prefix};
  //
  PolyHeuristic<mpz_class> AllArr = AllStandardHeuristic<mpz_class>();
  //
  std::vector<Delaunay_MPI_Entry<Tint, Tgroup>> ListDel =
    EnumerationDelaunayPolytopes<T,Tint,Tgroup>(comm, os, eData, AllArr);
  os << "We now have ListDel\n";
  //
  WriteFamilyDelaunay(OutFormat, OutFile, ListDel);
}

// clang-format off
#endif  // SRC_LATT_LATTICEDELAUNAY_H_
// clang-format on

// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_DATABANK_H_
#define SRC_DUALDESC_DATABANK_H_
// clang-format off
#include "basic_datafile.h"
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_DATABANK
#endif

#ifdef DISABLE_DEBUG_DATABANK
#undef DEBUG_DATABANK
#endif

template <typename Tkey, typename Tval>
std::pair<Tkey, Tval> Read_BankEntry(std::string const &Prefix) {
  using T = typename Tkey::value_type;
  using Tgroup = typename Tval::Tgroup;
  using Tint = typename Tgroup::Tint;
  std::string eFileEXT = Prefix + ".ext";
  std::string eFileGRP = Prefix + ".grp";
  std::string eFileOrbitSize = Prefix + ".orbitsize";
  std::string eFileNB = Prefix + ".nb";
  std::string eFileFF = Prefix + ".ff";
  //
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileEXT);
  size_t n_row = EXT.rows();
  //
  std::ifstream is_grp(eFileGRP);
  Tgroup GRP;
  is_grp >> GRP;
  //
  std::ifstream is_orbitsize(eFileOrbitSize);
  std::vector<Tint> ListPossOrbsize;
  is_orbitsize >> ListPossOrbsize;
  size_t n_factor = ListPossOrbsize.size();
  size_t n_bit_orbsize = get_matching_power(n_factor + 1);
  size_t delta = n_bit_orbsize + n_row;
  //
  FileNumber fn(eFileNB, false);
  size_t n_orbit = fn.getval();
  //
  FileFace ff(eFileFF, delta, n_orbit);
  vectface ListFace(delta);
  for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
    Face eFace = ff.getface(i_orbit);
    ListFace.push_back(eFace);
  }
#ifdef DEBUG_DATABANK
  std::cerr << "Read_BankEntry returning for Prefix=" << Prefix
            << " |EXT|=" << EXT.rows() << "/" << EXT.cols()
            << " |ListFace|=" << ListFace.size() << "\n";
#endif
  Tval eVal{std::move(GRP), std::move(ListPossOrbsize), std::move(ListFace)};
  return {std::move(EXT), std::move(eVal)};
}

template <typename T, typename Tgroup>
void Write_BankEntry(const std::string &Prefix, const MyMatrix<T> &EXT,
                     const TripleStore<Tgroup> &triple) {
  std::string eFileEXT = Prefix + ".ext";
  std::string eFileGRP = Prefix + ".grp";
  std::string eFileOrbitSize = Prefix + ".orbitsize";
  std::string eFileNB = Prefix + ".nb";
  std::string eFileFF = Prefix + ".ff";
  if (!FILE_IsFileMakeable(eFileEXT)) {
    std::cerr << "Error in Write_BankEntry: File eFileEXT=" << eFileEXT
              << " is not makeable\n";
    throw TerminalException{1};
  }
  //
  WriteMatrixFile(eFileEXT, EXT);
  size_t n_row = EXT.rows();
  //
  std::ofstream os_grp(eFileGRP);
  os_grp << triple.GRP;
  //
  std::ofstream os_orbitsize(eFileOrbitSize);
  os_orbitsize << triple.ListPossOrbsize;
  //
  FileNumber fn(eFileNB, true);
  size_t n_orbit = triple.ListFace.size();
  fn.setval(n_orbit);
  //
  FileFace ff(eFileFF, n_row);
  ff.direct_write(triple.ListFace.serial_get_std_vector_uint8_t());
}

template <typename Tkey, typename Tval>
void ReadingDatabaseFromPrefix(std::unordered_map<Tkey, Tval> &ListEnt,
                               bool const &Saving, std::string SavingPrefix,
                               [[maybe_unused]] std::ostream &os) {
  if (Saving) {
    size_t iOrbit = 0;
    while (true) {
      auto Prefix = SavingPrefix + "DualDesc" + std::to_string(iOrbit);
      std::string eFileBank = Prefix + ".ext";
      if (!IsExistingFile(eFileBank))
        break;
#ifdef DEBUG_DATABANK
      os << "Read iOrbit=" << iOrbit << " Prefix=" << Prefix << "\n";
#endif
      std::pair<Tkey, Tval> PairKV = Read_BankEntry<Tkey, Tval>(Prefix);
      ListEnt.emplace(std::move(PairKV.first), std::move(PairKV.second));
      iOrbit++;
    }
  }
}

// It is better to use std::unordered_map for the List of entries:
// This makes the check of equality rarer and instead uses the hash
// more strictly.
// Plus it is better because the tsl::sparse_map requires having the
// copy operator and we want to avoid that for the vectface.

template <typename Tkey, typename Tval> struct DataBank {
private:
  std::unordered_map<Tkey, Tval> ListEnt;
  bool Saving;
  std::string SavingPrefix;
  Tval TrivElt;
  std::ostream &os;

public:
  DataBank(const bool &_Saving, const std::string &_SavingPrefix,
           std::ostream &_os)
      : Saving(_Saving), SavingPrefix(_SavingPrefix), os(_os) {
    ReadingDatabaseFromPrefix(ListEnt, Saving, SavingPrefix, os);
  }
  void InsertEntry(Tkey &&eKey, Tval &&eVal) {
    // We have to face the situation that what we are trying to insert is
    // already present This can happen because of badly aligned heuristics
    if (ListEnt.contains(eKey)) {
#ifdef DEBUG_DATABANK
      os << "Exiting because the key is already present\n";
#endif
      return;
    }
    if (Saving) {
      size_t n_orbit = ListEnt.size();
      auto Prefix = SavingPrefix + "DualDesc" + std::to_string(n_orbit);
#ifdef DEBUG_DATABANK
      os << "Insert entry to file Prefix=" << Prefix << "\n";
#endif
      Write_BankEntry(Prefix, eKey, eVal);
    }
    ListEnt.emplace(std::move(eKey), std::move(eVal));
  }
  const Tval &GetDualDesc(const Tkey &eKey) const {
#ifdef DEBUG_DATABANK
    os << "Passing by GetDualDesc |ListEnt|=" << ListEnt.size() << "\n";
#endif
    auto iter = ListEnt.find(eKey);
    if (iter == ListEnt.end()) {
      // If returning empty then it means nothing has been found
      return TrivElt;
    }
    return iter->second;
  }
};

template <typename Tkey, typename Tval> struct TripleNKV {
  char nature;
  Tkey eKey;
  Tval eVal;
};

template <typename Tkey, typename Tval> struct PairKV {
  Tkey eKey;
  Tval eVal;
};

namespace boost::serialization {

template <class Archive, typename Tkey, typename Tval>
inline void serialize(Archive &ar, TripleNKV<Tkey, Tval> &triple,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("nature", triple.nature);
  ar &make_nvp("eKey", triple.eKey);
  ar &make_nvp("eVal", triple.eVal);
}

template <class Archive, typename Tkey, typename Tval>
inline void serialize(Archive &ar, PairKV<Tkey, Tval> &pair,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("eKey", pair.eKey);
  ar &make_nvp("eVal", pair.eVal);
}

// clang-format off
}  // namespace boost::serialization
// clang-format on

// clang-format off
#endif  // SRC_DUALDESC_DATABANK_H_
// clang-format on

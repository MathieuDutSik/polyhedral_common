// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_DATABANK_H_
#define SRC_DUALDESC_DATABANK_H_
// clang-format off
#include "basic_datafile.h"
#include <boost/asio.hpp>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

size_t get_matching_power(size_t const &val) {
  size_t pow = 1;
  size_t pos = 0;
  while (true) {
    if (pow >= val)
      return pos;
    pow *= 2;
    pos++;
  }
}

template <typename Tgroup_impl> struct TripleStore {
  using Tgroup = Tgroup_impl;
  using Tint = typename Tgroup::Tint;
  Tgroup GRP;
  std::vector<Tint> ListPossOrbsize;
  vectface ListFace;
};

namespace boost::serialization {

template <class Archive, typename Tgroup>
inline void serialize(Archive &ar, TripleStore<Tgroup> &triple,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("GRP", triple.GRP);
  ar &make_nvp("ListPossOrbsize", triple.ListPossOrbsize);
  ar &make_nvp("ListFace", triple.ListFace);
}

// clang-format off
}  // namespace boost::serialization
// clang-format on

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
  std::cerr << "Read_BankEntry returning for Prefix=" << Prefix
            << " |EXT|=" << EXT.rows() << "/" << EXT.cols()
            << " |ListFace|=" << ListFace.size() << "\n";
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
                               std::ostream &os) {
  if (Saving) {
    size_t iOrbit = 0;
    while (true) {
      std::string Prefix = SavingPrefix + "DualDesc" + std::to_string(iOrbit);
      std::string eFileBank = Prefix + ".ext";
      if (!IsExistingFile(eFileBank))
        break;
      os << "Read iOrbit=" << iOrbit << " Prefix=" << Prefix << "\n";
      std::pair<Tkey, Tval> PairKV = Read_BankEntry<Tkey, Tval>(Prefix);
      ListEnt.emplace(std::make_pair<Tkey, Tval>(std::move(PairKV.first),
                                                 std::move(PairKV.second)));
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
    if (ListEnt.count(eKey) > 0) {
      os << "Exiting because the key is already present\n";
      return;
    }
    if (Saving) {
      size_t n_orbit = ListEnt.size();
      std::string Prefix = SavingPrefix + "DualDesc" + std::to_string(n_orbit);
      os << "Insert entry to file Prefix=" << Prefix << "\n";
      Write_BankEntry(Prefix, eKey, eVal);
    }
    ListEnt.emplace(
        std::make_pair<Tkey, Tval>(std::move(eKey), std::move(eVal)));
  }
  const Tval &GetDualDesc(const Tkey &eKey) const {
    os << "Passing by GetDualDesc |ListEnt|=" << ListEnt.size() << "\n";
    typename std::unordered_map<Tkey, Tval>::const_iterator iter =
        ListEnt.find(eKey);
    if (iter == ListEnt.end()) {
      // If returning empty then it means nothing has been found
      return TrivElt;
    }
    return iter->second;
  }
};

template <typename T> T read_data(boost::asio::ip::tcp::socket &socket) {
  boost::asio::streambuf buf;
  boost::asio::read_until(socket, buf, "\n");
  std::string data = boost::asio::buffer_cast<const char *>(buf.data());
  T data_o;
  std::stringstream iss(data);
  boost::archive::text_iarchive ia(iss);
  ia >> data_o;
  return data_o;
}

template <typename T>
void send_data(boost::asio::ip::tcp::socket &socket, const T &data) {
  std::ostringstream archive_stream;
  boost::archive::text_oarchive archive(archive_stream);
  archive << data;
  std::string outbound_data = archive_stream.str() + "\n";
  boost::system::error_code error;
  boost::asio::write(socket, boost::asio::buffer(outbound_data), error);
  if (error) {
    std::cout << "send failed: " << error.message() << "\n";
    exit(1);
  }
}

template <typename T>
void send_data_atomic(const boost::asio::ip::tcp::endpoint &endpoint,
                      const T &data) {
  boost::asio::io_service io_service;
  boost::asio::ip::tcp::socket socket(io_service);
  socket.connect(endpoint);
  send_data<T>(socket, data);
  socket.close();
}

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

// ASIO client and server

template <typename Tkey, typename Tval>
void DataBankAsioServer(const bool &Saving, const std::string &SavingPrefix,
                        const uint16_t port, std::ostream &os) {
  std::unordered_map<Tkey, Tval> ListEnt;
  ReadingDatabaseFromPrefix(ListEnt, Saving, SavingPrefix, os);
  //
  boost::asio::io_service io_service;
  boost::asio::ip::tcp::endpoint endpoint(boost::asio::ip::tcp::v4(), port);
  boost::asio::ip::tcp::acceptor acceptor(io_service, endpoint);
  size_t n_iter = 0;
  while (true) {
    boost::asio::ip::tcp::socket socket = acceptor.accept();
    TripleNKV<Tkey, Tval> eTriple = read_data<TripleNKV<Tkey, Tval>>(socket);
    if (eTriple.nature == 'e') {
      os << "Terminating the Databank process\n";
      break;
    }
    if (eTriple.nature == 'i') {
      if (ListEnt.count(eTriple.eKey) > 0) {
        os << "The entry is already present\n";
      } else {
        if (Saving) {
          size_t n_orbit = ListEnt.size();
          std::string Prefix =
              SavingPrefix + "DualDesc" + std::to_string(n_orbit);
          os << "Insert entry to file Prefix=" << Prefix << "\n";
          Write_BankEntry(Prefix, eTriple.eKey, eTriple.eVal);
        }
        ListEnt.emplace(std::make_pair<Tkey, Tval>(std::move(eTriple.eKey),
                                                   std::move(eTriple.eVal)));
      }
    }
    if (eTriple.nature == 'g') {
      os << "Passing by GetDualDesc |ListEnt|=" << ListEnt.size() << "\n";
      typename std::unordered_map<Tkey, Tval>::const_iterator iter =
          ListEnt.find(eTriple.eKey);
      if (iter == ListEnt.end()) {
        // If returning empty then it means nothing has been found.
        send_data<Tval>(socket, Tval());
      } else {
        send_data<Tval>(socket, iter->second);
      }
    }
    os << "---------- " << n_iter << " -----------------------\n";
    n_iter++;
  }
}

template <typename Tkey, typename Tval> struct DataBankAsioClient {
private:
  uint16_t port;
  boost::asio::ip::tcp::endpoint endpoint;

public:
  DataBankAsioClient(const uint16_t &_port)
      : port(_port), endpoint(boost::asio::ip::tcp::v4(), port) {}
  void InsertEntry(Tkey &&eKey, Tval &&eVal) {
    send_data_atomic<TripleNKV<Tkey, Tval>>(
        endpoint, {'i', std::move(eKey), std::move(eVal)});
  }
  Tval GetDualDesc(const Tkey &eKey) const {
    boost::asio::io_service io_service;
    boost::asio::ip::tcp::socket socket(io_service);
    socket.connect(endpoint);
    send_data<TripleNKV<Tkey, Tval>>(socket, {'g', std::move(eKey), Tval()});
    return read_data<Tval>(socket);
  }
};

// clang-format off
#endif  // SRC_DUALDESC_DATABANK_H_
// clang-format on

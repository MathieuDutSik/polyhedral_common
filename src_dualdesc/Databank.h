#ifndef INCLUDE_DATABANK_H
#define INCLUDE_DATABANK_H

#include <boost/asio.hpp>


template<typename Tgroup_impl>
struct PairStore {
  using Tgroup = Tgroup_impl;
  Tgroup GRP;
  vectface ListFace;
};



namespace boost::serialization {

  template<class Archive, typename Tgroup>
  inline void serialize(Archive & ar, PairStore<Tgroup> & pair, const unsigned int version)
  {
    ar & make_nvp("GRP", pair.GRP);
    ar & make_nvp("ListFace", pair.ListFace);
  }

}




template<typename Tkey, typename Tval>
std::pair<Tkey, Tval> Read_BankEntry(std::string const& Prefix)
{
  using T = typename Tkey::value_type;
  using Tgroup = typename Tval::Tgroup;
  std::string eFileEXT = Prefix + ".ext";
  std::string eFileGRP = Prefix + ".grp";
  std::string eFileNB = Prefix + ".nb";
  std::string eFileFF = Prefix + ".ff";
 //
  std::ifstream is_ext(eFileEXT);
  MyMatrix<T> EXT = ReadMatrix<T>(is_ext);
  size_t n_row = EXT.rows();
  //
  std::ifstream is_grp(eFileGRP);
  Tgroup GRP;
  is_grp >> GRP;
  //
  FileNumber fn(eFileNB, false);
  size_t n_orbit = fn.getval();
  //
  FileFace ff(eFileFF, n_row, n_orbit);
  vectface ListFace(n_row);
  for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++) {
    Face eFace = ff.getface(i_orbit);
    ListFace.push_back(eFace);
  }
  std::cerr << " |EXT|=" << EXT.rows() << " |ListFace|=" << ListFace.size() << "\n";
  Tval eVal{std::move(GRP), std::move(ListFace)};
  return {std::move(EXT), std::move(eVal)};
}




template<typename T, typename Tgroup>
void Write_BankEntry(const std::string& Prefix, const MyMatrix<T>& EXT, const PairStore<Tgroup>& ePair)
{
  std::string eFileEXT = Prefix + ".ext";
  std::string eFileGRP = Prefix + ".grp";
  std::string eFileNB = Prefix + ".nb";
  std::string eFileFF = Prefix + ".ff";
  if (!FILE_IsFileMakeable(eFileEXT)) {
    std::cerr << "Error in Write_BankEntry: File eFileEXT=" << eFileEXT << " is not makeable\n";
    throw TerminalException{1};
  }
  //
  std::ofstream os_ext(eFileEXT);
  WriteMatrix(os_ext, EXT);
  size_t n_row = EXT.rows();
  //
  std::ofstream os_grp(eFileGRP);
  os_grp << ePair.GRP;
  //
  FileNumber fn(eFileNB, true);
  size_t n_orbit = ePair.ListFace.size();
  fn.setval(n_orbit);
  //
  FileFace ff(eFileFF, n_row);
  for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++)
    ff.setface(i_orbit, ePair.ListFace[i_orbit]);
}









// It is better to use std::unordered_map for the List of entries:
// This makes the check of equality rarer and instead uses the hash
// more strictly.
// Plus it is better because the tsl::sparse_map requires having the
// copy operator and we want to avoid that for the vectface.


template<typename Tkey, typename Tval>
struct DataBank {
private:
  std::unordered_map<Tkey, Tval> ListEnt;
  bool Saving;
  std::string SavingPrefix;
  Tval TrivElt;
public:
  DataBank(const bool& _Saving, const std::string& _SavingPrefix) : Saving(_Saving), SavingPrefix(_SavingPrefix)
  {
    if (Saving) {
      size_t iOrbit=0;
      while(true) {
        std::string eFileBank = SavingPrefix + "Dual_Desc_" + std::to_string(iOrbit) + ".nc";
        if (!IsExistingFile(eFileBank))
          break;
        std::cerr << "Read iOrbit=" << iOrbit << " FileBank=" << eFileBank << "\n";
        std::pair<Tkey, Tval> PairKV = Read_BankEntry<Tkey,Tval>(eFileBank);
        ListEnt.emplace(std::make_pair<Tkey, Tval>(std::move(PairKV.first), std::move(PairKV.second)));
        iOrbit++;
      }
    }
  }
  void InsertEntry(Tkey && eKey, Tval && eVal)
  {
    if (Saving) {
      size_t n_orbit = ListEnt.size();
      std::string Prefix = SavingPrefix + "DualDesc" + std::to_string(n_orbit);
      std::cerr << "Insert entry to file Prefix=" << Prefix << "\n";
      Write_BankEntry(Prefix, eKey, eVal);
    }
    ListEnt.emplace(std::make_pair<Tkey, Tval>(std::move(eKey), std::move(eVal)));
  }
  const Tval& GetDualDesc(const Tkey& eKey) const
  {
    std::cerr << "Passing by GetDualDesc |ListEnt|=" << ListEnt.size() << "\n";
    typename std::unordered_map<Tkey, Tval>::const_iterator iter = ListEnt.find(eKey);
    if (iter == ListEnt.end())
      return TrivElt; // If returning empty then it means nothing has been found.
    return iter->second;
  }
};




template<typename T>
T read_data(boost::asio::ip::tcp::socket & socket)
{
  boost::asio::streambuf buf;
  boost::asio::read_until( socket, buf, "\n" );
  std::string data = boost::asio::buffer_cast<const char*>(buf.data());
  T data_o;
  std::stringstream iss(data);
  boost::archive::text_iarchive ia(iss);
  ia >> data_o;
  return data_o;
}


template<typename T>
void send_data(boost::asio::ip::tcp::socket & socket, const T& data) {
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

template<typename T>
void send_data_atomic(const boost::asio::ip::tcp::endpoint& endpoint, const T& data) {
  boost::asio::io_service io_service;
  boost::asio::ip::tcp::socket socket(io_service);
  socket.connect(endpoint);
  send_data<T>(socket, data);
  socket.close();
}






template<typename Tkey, typename Tval>
struct TripleNKV {
  char nature;
  Tkey eKey;
  Tval eVal;
};


namespace boost::serialization {

  template<class Archive, typename Tkey, typename Tval>
  inline void serialize(Archive & ar, TripleNKV<Tkey,Tval> & triple, const unsigned int version)
  {
    ar & make_nvp("nature", triple.nature);
    ar & make_nvp("eKey", triple.eKey);
    ar & make_nvp("eVal", triple.eVal);
  }

}




template<typename Tkey, typename Tval>
struct DataBankServer {
private:
  std::unordered_map<Tkey, Tval> ListEnt;
  bool Saving;
  std::string SavingPrefix;
  short unsigned int port;
public:
DataBankServer(const bool& _Saving, const std::string& _SavingPrefix, const short unsigned int _port) : Saving(_Saving), SavingPrefix(_SavingPrefix), port(_port)
  {
    if (Saving) {
      size_t iOrbit=0;
      while(true) {
        std::string eFileBank = SavingPrefix + "Dual_Desc_" + std::to_string(iOrbit) + ".nc";
        if (!IsExistingFile(eFileBank))
          break;
        std::cerr << "Read iOrbit=" << iOrbit << " FileBank=" << eFileBank << "\n";
        std::pair<Tkey, Tval> PairKV = Read_BankEntry<Tkey,Tval>(eFileBank);
        ListEnt.emplace(std::make_pair<Tkey, Tval>(std::move(PairKV.first), std::move(PairKV.second)));
        iOrbit++;
      }
    }
    //
    boost::asio::io_service io_service;
    boost::asio::ip::tcp::endpoint endpoint(boost::asio::ip::tcp::v4(), port);
    boost::asio::ip::tcp::acceptor acceptor(io_service, endpoint);
    size_t n_iter=0;
    while (true) {
      boost::asio::ip::tcp::socket socket = acceptor.accept();
      TripleNKV<Tkey, Tval> eTriple = read_data<TripleNKV<Tkey,Tval>>(socket);
      if (eTriple.nature == 'e') {
        std::cerr << "Terminating the Databank process\n";
        break;
      }
      if (eTriple.nature == 'i') {
        if (Saving) {
          size_t n_orbit = ListEnt.size();
          std::string Prefix = SavingPrefix + "DualDesc" + std::to_string(n_orbit);
          std::cerr << "Insert entry to file Prefix=" << Prefix << "\n";
          Write_BankEntry(Prefix, eTriple.eKey, eTriple.eVal);
        }
        ListEnt.emplace(std::make_pair<Tkey, Tval>(std::move(eTriple.eKey), std::move(eTriple.eVal)));
      }
      if (eTriple.nature == 'g') {
        std::cerr << "Passing by GetDualDesc |ListEnt|=" << ListEnt.size() << "\n";
        typename std::unordered_map<Tkey, Tval>::const_iterator iter = ListEnt.find(eTriple.eKey);
        if (iter == ListEnt.end())
          send_data<Tval>(socket, Tval()); // If returning empty then it means nothing has been found.
        send_data<Tval>(socket, iter->second);
      }
      std::cerr << "------------------------------ " << n_iter << " ------------------------------\n";
      n_iter++;
    }
  }
};



template<typename Tkey, typename Tval>
struct DataBankClient {
private:
  std::unordered_map<Tkey, Tval> ListEnt;
  short unsigned int port;
  boost::asio::ip::tcp::endpoint endpoint;
public:
  DataBankClient(const short unsigned int& _port) : port(_port), endpoint(boost::asio::ip::tcp::v4(), port)
  {
  }
  void InsertEntry(Tkey && eKey, Tval && eVal)
  {
    send_data_atomic<TripleNKV<Tkey,Tval>>(endpoint, {'i', std::move(eKey), std::move(eVal)});
  }
  Tval GetDualDesc(const Tkey& eKey) const
  {
    boost::asio::io_service io_service;
    boost::asio::ip::tcp::socket socket(io_service);
    socket.connect(endpoint);
    send_data<TripleNKV<Tkey,Tval>>(socket, {'g', std::move(eKey), Tval()});
    return read_data<Tval>(socket);
  }
};

#endif

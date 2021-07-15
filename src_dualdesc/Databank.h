#ifndef INCLUDE_DATABANK_H
#define INCLUDE_DATABANK_H


#include <boost/asio.hpp>


#include <boost/archive/tmpdir.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>



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
      std::string eFile = SavingPrefix + "DualDesc" + std::to_string(n_orbit) + ".nc";
      std::cerr << "Insert entry to file eFile=" << eFile << "\n";
      Write_BankEntry(eFile, eKey, eVal);
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


template<typename Tkey, typename Tval>
struct TripleNKV {
  char nature;
  Tkey eKey;
  Tval eVal;
};


template<typename Tkey, typename Tval>
struct DataBankServer {
private:
  std::unordered_map<Tkey, Tval> ListEnt;
  bool Saving;
  std::string SavingPrefix;
  short unsigned int port;
  Tval TrivElt;
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
          std::string eFile = SavingPrefix + "DualDesc" + std::to_string(n_orbit) + ".nc";
          std::cerr << "Insert entry to file eFile=" << eFile << "\n";
          Write_BankEntry(eFile, eTriple.eKey, eTriple.eVal);
        }
        ListEnt.emplace(std::make_pair<Tkey, Tval>(std::move(eTriple.eKey), std::move(eTriple.eVal)));
      }
      if (eTriple.nature == 'g') {
        std::cerr << "Passing by GetDualDesc |ListEnt|=" << ListEnt.size() << "\n";
        typename std::unordered_map<Tkey, Tval>::const_iterator iter = ListEnt.find(eTriple.eKey);
        if (iter == ListEnt.end())
          send_data<Tval>(socket, TrivElt); // If returning empty then it means nothing has been found.
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
  Tval TrivElt;
  short unsigned int port;
public:
  DataBankClient(const short unsigned int& _port) : port(_port)
  {
  }
  void InsertEntry(Tkey && eKey, Tval && eVal)
  {
    boost::asio::io_service io_service;
    boost::asio::ip::tcp::endpoint endpoint(boost::asio::ip::tcp::v4(), port);
    boost::asio::ip::tcp::socket socket(io_service);
    socket.connect(endpoint);
    send_data<TripleNKV>(socket, {'i', eKey, eVal});
    socket.close();
  }
  const Tval& GetDualDesc(const Tkey& eKey) const
  {
    boost::asio::io_service io_service;
    boost::asio::ip::tcp::endpoint endpoint(boost::asio::ip::tcp::v4(), port);
    boost::asio::ip::tcp::socket socket(io_service);
    socket.connect(endpoint);
    send_data<TripleNKV>(socket, {'g', eKey, TrivElt});
    Tval eVal = read_data<Tval>(socket);
    socket.close();
  }
};

#endif

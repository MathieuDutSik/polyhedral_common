// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_DATABANK_ASIO_H_
#define SRC_DUALDESC_DATABANK_ASIO_H_
// clang-format off
#include "basic_datafile.h"
#include <boost/asio.hpp>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

template <typename T> T read_data(boost::asio::ip::tcp::socket &socket) {
  boost::asio::streambuf buf;
  boost::asio::read_until(socket, buf, "\n");
  const char *raw_data = static_cast<const char *>(buf.data().data());
  std::string data(raw_data, buf.size());
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
    std::cerr << "send failed: " << error.message() << "\n";
    throw TerminalException{1};
  }
}

template <typename T>
void send_data_atomic(const boost::asio::ip::tcp::endpoint &endpoint,
                      const T &data) {
  boost::asio::io_context io_context;
  boost::asio::ip::tcp::socket socket(io_context);
  socket.connect(endpoint);
  send_data<T>(socket, data);
  socket.close();
}

// ASIO client and server

template <typename Tkey, typename Tval>
void DataBankAsioServer(const bool &Saving, const std::string &SavingPrefix,
                        const uint16_t port, std::ostream &os) {
  std::unordered_map<Tkey, Tval> ListEnt;
  ReadingDatabaseFromPrefix(ListEnt, Saving, SavingPrefix, os);
  //
  boost::asio::io_context io_context;
  boost::asio::ip::tcp::endpoint endpoint(boost::asio::ip::tcp::v4(), port);
  boost::asio::ip::tcp::acceptor acceptor(io_context, endpoint);
#ifdef DEBUG_DATABANK
  size_t n_iter = 0;
#endif
  while (true) {
    boost::asio::ip::tcp::socket socket = acceptor.accept();
    TripleNKV<Tkey, Tval> eTriple = read_data<TripleNKV<Tkey, Tval>>(socket);
    if (eTriple.nature == 'e') {
#ifdef DEBUG_DATABANK
      os << "Terminating the Databank process\n";
#endif
      break;
    }
    if (eTriple.nature == 'i') {
      if (ListEnt.contains(eTriple.eKey)) {
#ifdef DEBUG_DATABANK
        os << "The entry is already present\n";
#endif
      } else {
        if (Saving) {
          size_t n_orbit = ListEnt.size();
          std::string Prefix =
              SavingPrefix + "DualDesc" + std::to_string(n_orbit);
#ifdef DEBUG_DATABANK
          os << "Insert entry to file Prefix=" << Prefix << "\n";
#endif
          Write_BankEntry(Prefix, eTriple.eKey, eTriple.eVal);
        }
        ListEnt.emplace(std::make_pair<Tkey, Tval>(std::move(eTriple.eKey),
                                                   std::move(eTriple.eVal)));
      }
    }
    if (eTriple.nature == 'g') {
#ifdef DEBUG_DATABANK
      os << "Passing by GetDualDesc |ListEnt|=" << ListEnt.size() << "\n";
#endif
      auto iter = ListEnt.find(eTriple.eKey);
      if (iter == ListEnt.end()) {
        // If returning empty then it means nothing has been found.
        send_data<Tval>(socket, Tval());
      } else {
        send_data<Tval>(socket, iter->second);
      }
    }
#ifdef DEBUG_DATABANK
    os << "---------- " << n_iter << " -----------------------\n";
    n_iter++;
#endif
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
    boost::asio::io_context io_context;
    boost::asio::ip::tcp::socket socket(io_context);
    socket.connect(endpoint);
    send_data<TripleNKV<Tkey, Tval>>(socket, {'g', std::move(eKey), Tval()});
    return read_data<Tval>(socket);
  }
};

// clang-format off
#endif  // SRC_DUALDESC_DATABANK_ASIO_H_
// clang-format on

// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_DATABANK_MPI_H_
#define SRC_DUALDESC_DATABANK_MPI_H_
// clang-format off
#include "Databank.h"
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
// clang-format on

static const int val_mpi_bank_end = 2045;
static const int tag_mpi_bank_end = 977;
static const int tag_mpi_bank_insert = 983;
static const int tag_mpi_bank_request = 991;
static const int tag_mpi_bank_reception = 997;


// MPI client and server

template <typename Tkey, typename Tval> struct DataBankMpiClient {
private:
  boost::mpi::communicator &comm;
  int dest_iproc;

public:
  DataBankMpiClient(boost::mpi::communicator &comm)
      : comm(comm), dest_iproc(comm.size() - 1) {}
  void InsertEntry(Tkey &&eKey, Tval &&eVal) {
    PairKV<Tkey, Tval> pair{std::move(eKey), std::move(eVal)};
    comm.send(dest_iproc, tag_mpi_bank_insert, pair);
  }
  Tval GetDualDesc(const Tkey &eKey) const {
    comm.send(dest_iproc, tag_mpi_bank_request, eKey);
    Tval val;
    comm.recv(dest_iproc, tag_mpi_bank_reception, val);
    return val;
  }
};

template <typename Tkey, typename Tval>
void DataBankMpiServer(boost::mpi::communicator &comm, const bool &Saving,
                       const std::string &SavingPrefix, std::ostream &os) {
  std::unordered_map<Tkey, Tval> ListEnt;
  ReadingDatabaseFromPrefix(ListEnt, Saving, SavingPrefix, os);
  //
  while (true) {
    boost::mpi::status msg = comm.probe();
    if (msg.tag() == tag_mpi_bank_end) {
      os << "Receiving termination message\n";
      int val_i;
      comm.recv(msg.source(), msg.tag(), val_i);
      if (val_i != val_mpi_bank_end) {
        std::cerr << "Error in the received value\n";
        std::cerr << "val_i=" << val_i
                  << " val_mpi_bank_end=" << val_mpi_bank_end << "\n";
        throw TerminalException{1};
      }
      os << "Terminating the DataBankMpiServer\n";
      break;
    }
    if (msg.tag() == tag_mpi_bank_insert) {
      os << "Receiving insert (key,val) from " << msg.source() << "\n";
      PairKV<Tkey, Tval> pair;
      comm.recv(msg.source(), msg.tag(), pair);
      if (ListEnt.count(pair.eKey) > 0) {
        os << "The entry is already present\n";
        os << "Not stopping but suboptimal. Maybe two threads computed the "
              "same thing. Otherwise bug.\n";
      } else {
        if (Saving) {
          size_t n_orbit = ListEnt.size();
          std::string Prefix =
              SavingPrefix + "DualDesc" + std::to_string(n_orbit);
          os << "Insert entry to file Prefix=" << Prefix << "\n";
          Write_BankEntry(Prefix, pair.eKey, pair.eVal);
        }
        ListEnt.emplace(std::make_pair<Tkey, Tval>(std::move(pair.eKey),
                                                   std::move(pair.eVal)));
      }
    }
    if (msg.tag() == tag_mpi_bank_request) {
      os << "Receiving request key from " << msg.source()
         << " |ListEnt|=" << ListEnt.size() << "\n";
      Tkey key;
      comm.recv(msg.source(), msg.tag(), key);
      typename std::unordered_map<Tkey, Tval>::const_iterator iter =
          ListEnt.find(key);
      if (iter == ListEnt.end()) {
        os << "Not found\n";
        comm.send(msg.source(), tag_mpi_bank_reception, Tval());
      } else {
        os << "Existing entry found. Returning it\n";
        comm.send(msg.source(), tag_mpi_bank_reception, iter->second);
      }
    }
  }
}

// clang-format off
#endif  // SRC_DUALDESC_DATABANK_MPI_H_
// clang-format on

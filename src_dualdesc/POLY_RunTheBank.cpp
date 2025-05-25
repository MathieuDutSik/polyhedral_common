// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "Permutation.h"
#include "Group.h"
#include "POLY_RecursiveDualDesc.h"
// clang-format on

boost::asio::ip::tcp::endpoint endpoint_bank;

template <typename Tkey, typename Tval>
void signal_callback_handler_bank(int signum) {
  std::cout << "Caught signal " << signum << "\n";
  std::cout << "We are going to exit from the bank hopefully\n";
  TripleNKV<Tkey, Tval> triple{'e', Tkey(), Tval()};
  send_data_atomic<TripleNKV<Tkey, Tval>>(endpoint_bank, triple);
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_BankingSystem();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Temp_ThreadedADM [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    using T = mpq_class;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    using Tkey = MyMatrix<T>;
    using Tval = TripleStore<Tgroup>;
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    SingleBlock const& BlockPROC = eFull.get_block("PROC");
    bool Saving = BlockPROC.get_bool("Saving");
    std::string SavingPrefix = BlockPROC.get_string("Prefix");
    int port_i = BlockPROC.get_int("port");
    std::cerr << "port_i=" << port_i << "\n";
    uint16_t port = port_i;
    //
    endpoint_bank =
        boost::asio::ip::tcp::endpoint(boost::asio::ip::tcp::v4(), port);
    auto process_signal = [](int signum) {
      signal_callback_handler_bank<Tkey, Tval>(signum);
    };
    signal(SIGINT, process_signal);
    DataBankAsioServer<Tkey, Tval>(Saving, SavingPrefix, port, std::cerr);
    std::cerr << "Normal termination of POLY_RunTheBank\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_RunTheBank\n";
    exit(e.eVal);
  }
  runtime(time1);
}

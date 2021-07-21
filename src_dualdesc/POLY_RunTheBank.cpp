#include "Permutation.h"
#include "Group.h"
#include "Group_serialization.h"
#include "POLY_RecursiveDualDesc.h"

int main(int argc, char *argv[])
{
  try {
    FullNamelist eFull=NAMELIST_GetStandard_BankingSystem();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Temp_ThreadedADM [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    using T = mpq_class;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt,Tint>;
    using Tkey = MyMatrix<T>;
    using Tval = PairStore<Tgroup>;
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    SingleBlock BlockPROC = eFull.ListBlock.at("PROC");
    bool Saving=BlockPROC.ListBoolValues.at("Saving");
    std::string SavingPrefix=BlockPROC.ListStringValues.at("SavingPrefix");
    int port_i=BlockPROC.ListIntValues.at("port");
    std::cerr << "port_i=" << port_i << "\n";
    short unsigned int port = port_i;
    //
    auto process_signal=[](int signum) {
      signal_callback_handler_bank<Tkey,Tval>(signum);
    };
    signal(SIGINT, process_signal);
    DataBankServer<Tkey,Tval> BankServer(Saving, SavingPrefix, port);
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}

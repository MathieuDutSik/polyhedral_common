#include "POLY_RecursiveDualDesc.h"
#include "Permutation.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "DD_ConvertGAP_to_NC [DATA_EXT] [DATA_FAC] [DATA_NC]\n";
      std::cerr << "\n";
      std::cerr << "DATA_EXT (in) : The polytope vertices and group\n";
      std::cerr << "DATA_FAC (in) : The list of orbits of facets\n";
      std::cerr << "DATA_NC (out): The ext, grp, ListOrbit in a netcdf file\n";
      return -1;
    }
    std::string File_EXT = argv[1];
    std::string File_FAC = argv[2];
    std::string File_NC  = argv[3];
    //
    using T = mpq_class;
    using Tint = mpz_class;
    using Telt = permutalib::DoubleSidedPerm;
    using Tgroup=permutalib::Group<Telt,Tint>;
    using Telt = typename Tgroup::Telt;
    datagap::DataGAP<T,Telt> dataEXT = ParseGAPFile(File_EXT);
    datagap::DataGAP<T,Telt> dataFAC = ParseGAPFile(File_FAC);
    EquivariantDualDescription<T,Tgroup> RecEXT_GRP_LOrb = ConvertGAPread_EquivDualDesc(dataEXT, dataFAC);
    Write_EquivDualDesc(RecEXT_GRP_LOrb, File_NC);
  }
}

// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "POLY_RecursiveDualDesc.h"
#include "Permutation.h"

template <typename T, typename Tgroup> struct EquivariantDualDescription {
  MyMatrix<T> EXT;
  Tgroup GRP;
  vectface ListFace;
};

template <typename T, typename Tgroup>
EquivariantDualDescription<T, Tgroup> ConvertGAPread_EquivDualDesc(
    datagap::DataGAP<T, typename Tgroup::Telt> const &dataEXT,
    datagap::DataGAP<T, typename Tgroup::Telt> const &dataFAC) {
  if (dataEXT.Nature != datagap::int_record) {
    std::cerr << "For EquivDualDesc, we need to have a record as entry\n";
    throw TerminalException{1};
  }
  int pos_EXT = -1;
  int pos_GRP = -1;
  int n_pos = dataEXT.ListRec.size();
  for (int pos = 0; pos < n_pos; pos++) {
    if (dataEXT.ListRec[pos].first == "EXT")
      pos_EXT = pos;
    if (dataEXT.ListRec[pos].first == "Group")
      pos_GRP = pos;
  }
  if (pos_EXT == -1) {
    std::cerr << "Failed to find entry EXT in the record\n";
    throw TerminalException{1};
  }
  if (pos_GRP == -1) {
    std::cerr << "Failed to find entry Group in the record\n";
    throw TerminalException{1};
  }
  MyMatrix<T> EXT =
      datagap::ConvertGAPread_MyMatrixT(dataEXT.ListRec[pos_EXT].second);
  int n_rows = EXT.rows();
  Tgroup GRP = datagap::ConvertGAPread_PermutationGroup<T, Tgroup>(
      dataEXT.ListRec[pos_GRP].second, n_rows);
  //
  vectface ListFace = ConvertGAPread_ListFace(dataFAC, n_rows);
  //
  return {std::move(EXT), std::move(GRP), std::move(ListFace)};
}

int main(int argc, char *argv[]) {
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
    std::string File_NC = argv[3];
    //
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Telt = permutalib::DoubleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    datagap::DataGAP<T, Telt> dataEXT =
        datagap::ParseGAPFile<T, Telt>(File_EXT);
    datagap::DataGAP<T, Telt> dataFAC =
        datagap::ParseGAPFile<T, Telt>(File_FAC);
    EquivariantDualDescription<T, Tgroup> RecEXT_GRP_LOrb =
        ConvertGAPread_EquivDualDesc<T, Tgroup>(dataEXT, dataFAC);
    Write_EquivDualDesc(RecEXT_GRP_LOrb, File_NC);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}

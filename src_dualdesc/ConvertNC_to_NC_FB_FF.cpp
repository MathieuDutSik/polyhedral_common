#include "Permutation.h"
#include "Group.h"
#include "POLY_RecursiveDualDesc.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "ConvertNC_to_NC_FB_FF [Prefix]\n";
      std::cerr << "\n";
      std::cerr << "Prefix (inout) : The prefix for the conversion\n";
      return -1;
    }
    std::string Prefix = argv[1];
    std::string eFileNC = Prefix + ".nc";
    std::string eFileFB = Prefix + ".fb";
    std::string eFileFF = Prefix + ".ff";
    if (!IsExistingFile(eFileNC)) {
      std::cerr << "The file eFileNC=" << eFileNC << " should be present\n";
      throw TerminalException{1};
    }
    if (IsExistingFile(eFileFB) || IsExistingFile(eFileFF)) {
      std::cerr << "The file eFileFB=" << eFileFB << " should be absent\n";
      std::cerr << "The file eFileFF=" << eFileFF << " should be absent\n";
      throw TerminalException{1};
    }
    //
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = int16_t;
    using Telt = permutalib::DoubleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt,Tint>;
    using Torbsize=uint16_t;
    netCDF::NcFile dataFile;
    //
    dataFile.open(eFileNC, netCDF::NcFile::read);
    MyMatrix<T> EXT = POLY_NC_ReadPolytope<T>(dataFile);
    Tgroup GRP = POLY_NC_ReadGroup<Tgroup>(dataFile);
    netCDF::NcDim eDim = dataFile.getDim("n_orbit");
    size_t n_orbit = eDim.getSize();
    std::vector<SingleEntryStatus<Tint>> ListOrb;
    for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++) {
      SingleEntryStatus<Tint> eEnt = POLY_NC_ReadSingleEntryStatus<Tint>(dataFile, i_orbit);
      ListOrb.push_back(eEnt);
    }
    dataFile.close();
    RemoveFileIfExist(eFileNC);
    //
    std::map<Tidx, int> LFact = GRP.factor_size();
    size_t n_factor = 1;
    for (auto & kv : LFact) {
      n_factor *= (1 + kv.second);
    }
    size_t n_act = GRP.n_act();
    size_t n_bit_orbsize = get_matching_power(n_factor + 1);
    std::vector<Tint> ListPossOrbsize = GetAllPossibilities<Tidx,Tint>(LFact);
    //
    UNORD_MAP<Tint,Torbsize> OrbSize_Map;
    auto GetOrbSizeIndex=[&](Tint const& orbSize) ->  Torbsize {
      Torbsize & idx = OrbSize_Map[orbSize];
      if (idx == 0) { // A rare case. The linear loop should be totally ok
        auto set=[&]() -> int {
          for (size_t u=0; u<ListPossOrbsize.size(); u++)
            if (ListPossOrbsize[u] == orbSize) {
              return u + 1;
            }
          return 0;
        };
        idx = set();
      }
      return idx - 1;
    };
    struct SingEnt {
      Face face;
      Torbsize idx_orb;
    };
    auto SingEntToFace=[&](SingEnt const& eEnt) -> Face {
      Face f(n_act + n_bit_orbsize);
      for (size_t i=0; i<n_act; i++)
        f[i] = eEnt.face[i];
      size_t work_idx = eEnt.idx_orb;
      size_t i_acc = n_act;
      for (size_t i=0; i<n_bit_orbsize; i++) {
        bool val = work_idx % 2;
        f[i_acc] = val;
        i_acc++;
        work_idx = work_idx / 2;
      }
      return f;
    };
    //
    FileBool fb(eFileFB);
    FileFace ff(eFileFF, n_act + n_bit_orbsize);
    dataFile.open(eFileNC, netCDF::NcFile::replace, netCDF::NcFile::nc4);
    POLY_NC_WritePolytope(dataFile, EXT);
    bool orbit_setup = false;
    bool orbit_status = false;
    POLY_NC_WriteGroup(dataFile, GRP, orbit_setup, orbit_status);
    POLY_NC_SetNbOrbit(dataFile);
    POLY_NC_WriteNbOrbit(dataFile, n_orbit);
    for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++) {
      SingleEntryStatus<Tint> eEnt = ListOrb[i_orbit];
      Torbsize idx_orb = GetOrbSizeIndex(eEnt.OrbSize);
      SingEnt eEntS{eEnt.face, idx_orb};
      Face f = SingEntToFace(eEntS);
      fb.setbit(i_orbit, eEnt.status);
      ff.setface(i_orbit, f);
    }
    //
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}

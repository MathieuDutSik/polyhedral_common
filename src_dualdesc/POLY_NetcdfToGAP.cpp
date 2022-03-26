#include "Group.h"
#include "POLY_RecursiveDualDesc.h"
#include "Permutation.h"
int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "POLY_NetcdfToGAP [file1] [file2]\n";
      return -1;
    }
    std::string FileI = argv[1];
    std::string FileO = argv[2];
    //
    using T = mpq_class;
    netCDF::NcFile dataFile(FileI, netCDF::NcFile::read);
    MyMatrix<T> EXT = POLY_NC_ReadPolytope<T>(dataFile);
    //
    WriteMatrixGAPfile(FileO, EXT);
  } catch (netCDF::exceptions::NcInvalidCoords &e) {
    std::cerr << "e complaint=" << e.what() << "\n";
    exit(1);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}

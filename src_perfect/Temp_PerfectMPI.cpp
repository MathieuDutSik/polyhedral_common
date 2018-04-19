#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
namespace mpi = boost::mpi;

#include "MAT_Matrix.h"

int main()
{
  using T=mpq_class;
  using Tint=mpz_class;
  //
  mpi::environment env;
  mpi::communicator world;
  int irank=world.rank();
  int size=world.size();
  std::string eFileO="LOG_" + irank;
  std::ofstream log(eFileO);
  log << "Initial log entry\n";
  //
  std::string FileMatrix = "ListMatrix";
  std::ifstream is(FileMatrix);
  int nbMatrixStart;
  is >> nbMatrixStart;
  struct KeyData {
    int StatusTreatedForm; // 0: untreated, 1: treated but status not written on disk, 2: done and treated
    MyMatrix<Tint> TheMat;
  };
  //
  std::vector<KeyData> ListCases;
  auto fInsert=[&](MyMatrix<Tint> const& NewMat) -> void {
    for (auto & fCase : ListCases) {
      if (fCace.TheMat == NewMat)
	return;
    }
    ListCases.push_back({0,NewMat});
    os << "Inserting new form, now we have |ListCases|=" << ListCases.size() << "\n";
  };
  for (int iMatStart=0; iMatStart<nbMatrixStart; iMatStart++) {
    int eStatus;
    is >> eStatus;
    MyMatrix<Tint> TheMat = ReadMatrix<Tint>(is);
    KeyData eData{eStatus,TheMat};
    int KeyInv=IntegerDiscriminantInvariant(TheMat);
    int res=KeyInv % size;
    if (res == irank)
      ListCases.push_back(eData);
  }
  os << "Reading finished, we have |ListCases|=" << ListCases.size() << "\n";
  //
  while(true) {
    boost::optional<mpi::status> prob = world.iprobe();
    if (prob) {
      if (prob->tag() == new_form) {
	
      }
    }
    else {
      
    }
    
  }
  
}

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
namespace mpi = boost::mpi;

#include "MAT_Matrix.h"

static tag_new_form = 37;


int main()
{
  using T=mpq_class;
  using Tint=mpz_class;
  struct TypePerfectExch {
    MyMatrix<Tint> eMat;
    int incd; // the number of shortest vectors divided by 2
  };
  //
  FullNamelist eFull = NAMELIST_GetStandard_ENUMERATE_PERFECT();
  std::string eFileName = "perfectenum.nml";
  NAMELIST_ReadNamelistFile(eFileName, eFull);
  SingleBlock BlDATA = eFull.ListBlock["DATA"];
  int n=BlDATA.ListIntValues.at("n");
  int MaxNumberFlyingMessage=BlDATA.ListIntValues.at("MaxNumerFlyingMessage");
  //
  boost::mpi::environment env;
  boost::mpi::communicator world;
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
  };
  //
  std::map<TypePerfectExch,KeyData> ListCases;
  auto fInsert=[&](TypePerfectExch const& NewMat) -> void {
    auto it = ListCases.find(NewMat);
    if (it != ListCases.end())
      return;
    ListCases[NewMat] = 0;
    os << "Inserting new form, now we have |ListCases|=" << ListCases.size() << "\n";
  };
  auto fTrySend=[&](TypePerfectMatrix const& NewMat) -> void {
    
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
  std::vector<boost::mpi::request> ListRequest(MaxNumberFlyingMessage);
  while(true) {
    boost::optional<boost::mpi::status> prob = world.iprobe();
    if (prob) {
      if (prob->tag() == tag_new_form) {
	std::vector<MyMatrix<Tint>> data;
	world.recv(msg->source(), msg->tag(), data);
	for (auto & eMat : data)
	  fInsert(eMat);
      }
    }
    else {
      

      
      std::vector<MyMatrix<Tint>> ListAdjacent = 
    }
    
  }
  
}

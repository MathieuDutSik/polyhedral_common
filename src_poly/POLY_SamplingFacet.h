#ifndef TEMP_SAMPLING_FACET_H
#define TEMP_SAMPLING_FACET_H


#include "POLY_PolytopeFct.h"


struct recSamplingOption {
  int critlevel;
  int maxnbcall;
  int maxnbsize;
  std::string prog;
};



template<typename T>
std::vector<Face> Kernel_DUALDESC_SamplingFacetProcedure(MyMatrix<T> const& EXT, recSamplingOption const& eOption, int & nbCall)
{
  int dim=RankMat(EXT);
  int len=EXT.rows();
  std::string prog=eOption.prog;
  int critlevel=eOption.critlevel;
  int maxnbcall=eOption.maxnbcall;
  int maxnbsize=eOption.maxnbsize;
  std::cerr << "critlevel=" << critlevel << " prog=" << prog << " maxnbcall=" << maxnbcall << "\n";
  auto IsRecursive=[&]() -> bool {
    if (len < critlevel)
      return false;
    if (dim < 15)
      return false;
    return true;
  };
  bool DoRecur=IsRecursive();
  std::vector<Face> ListFace;
  std::vector<int> ListStatus;
  auto FuncInsert=[&](Face const& eFace) -> void {
    for (auto & fFace : ListFace) {
      if (fFace.count() == eFace.count())
	return;
    }
    ListFace.push_back(eFace);
    ListStatus.push_back(0);
  };
  std::cerr << "dim=" << dim << "  len=" << len << "\n";
  if (!DoRecur) {
    std::vector<Face> ListIncd;
    if (prog == "lrs")
      ListIncd=lrs::DualDescription_temp_incd(EXT);
    if (prog == "cdd") {
      ListIncd=cdd::DualDescription_incd(EXT);
    }
    for (auto & eFace : ListIncd)
      FuncInsert(eFace);
    std::cerr << "DirectDualDesc |ListFace|=" << ListFace.size() << "\n";
    nbCall++;
    return ListFace;
  }
  Face eInc=FindOneInitialVertex(EXT);
  FuncInsert(eInc);
  while(true) {
    int nbCases=ListFace.size();
    bool IsFinished=true;
    for (int iC=0; iC<nbCases; iC++)
      if (ListStatus[iC] == 0) {
	nbCall++;  // we liberally increase the value
	IsFinished=false;
	ListStatus[iC]=1;
	Face eFace=ListFace[iC];
	MyMatrix<T> EXTred=SelectRow(EXT, eFace);
	std::vector<Face> ListRidge=Kernel_DUALDESC_SamplingFacetProcedure(EXTred, eOption, nbCall);
	for (auto & eRidge : ListRidge) {
	  Face eFlip=ComputeFlipping(EXT, eFace, eRidge);
	  FuncInsert(eFlip);
	}
	if (maxnbsize != -1) {
	  int siz=ListFace.size();
	  if (maxnbsize > siz) {
	    std::cerr << "Ending by maxsize criterion\n";
	    std::cerr << "siz=" << siz << " maxnbsize=" << maxnbsize << "\n";
	    return ListFace;
	  }
	}
	if (maxnbcall != -1) {
	  if (nbCall > maxnbcall) {
	    std::cerr << "Ending by maxnbcall\n";
	    return ListFace;
	  }
	}
      }
    if (IsFinished)
      break;
  }
  std::cerr << "RecursiveDualDesc |ListFace|=" << ListFace.size() << "\n";
  return ListFace;
}


template<typename T>
std::vector<Face> DUALDESC_SamplingFacetProcedure(MyMatrix<T> const& EXT, std::vector<std::string> const& ListOpt)
{
  std::string prog="lrs";
  int critlevel=50;
  int maxnbcall=-1;
  int maxnbsize=20;
  for (auto & eOpt : ListOpt) {
    std::vector<std::string> ListStrB=STRING_Split(eOpt, "_");
    if (ListStrB.size() == 2) {
      if (ListStrB[0] == "prog")
	prog=ListStrB[1];
      if (ListStrB[0] == "critlevel")
	std::istringstream(ListStrB[1]) >> critlevel;
      if (ListStrB[0] == "maxnbcall")
	std::istringstream(ListStrB[1]) >> maxnbcall;
      if (ListStrB[0] == "maxnbsize")
	std::istringstream(ListStrB[1]) >> maxnbsize;
    }
  }
  if (prog != "lrs" && prog != "cdd") {
    std::cerr << "We have prog=" << prog << "\n";
    std::cerr << "but the only allowed input formats are lrs and cdd\n";
    throw TerminalException{1};
  }
  recSamplingOption eOption;
  eOption.maxnbcall=maxnbcall;
  eOption.prog=prog;
  eOption.critlevel=critlevel;
  eOption.maxnbsize=maxnbsize;
  int nbcall=0;
  return Kernel_DUALDESC_SamplingFacetProcedure(EXT, eOption, nbcall);
}



template<typename T>
std::vector<Face> DirectComputationInitialFacetSet(MyMatrix<T> const& EXT, std::string const& ansSamp)
{
  bool WeAreDone=false;
  std::vector<Face> ListIncd;
  std::vector<std::string> ListStr=STRING_Split(ansSamp, ":");
  std::string ansOpt=ListStr[0];
  if (ansOpt == "lp_cdd") {
    // So possible format is lp_cdd:iter_100
    int iter=10;
    if (ListStr.size() > 1) {
      std::vector<std::string> ListStrB=STRING_Split(ListStr[1], "_");
      if (ListStrB.size() == 2 && ListStrB[0] == "iter")
	std::istringstream(ListStrB[1]) >> iter;
    }
    for (int i=0; i<iter; i++) {
      Face eFace=FindOneInitialVertex(EXT);
      ListIncd.push_back(eFace);
    }
    WeAreDone=true;
  }
  if (ansOpt == "sampling") {
    std::vector<std::string> ListOpt;
    for (int i=1; i<int(ListStr.size()); i++)
      ListOpt.push_back(ListStr[i]);
    ListIncd=DUALDESC_SamplingFacetProcedure(EXT, ListOpt);
    WeAreDone=true;
  }
  if (ansOpt == "lrs_limited") {
    int upperlimit=100;
    // So possible format is lrs_limited:upperlimit_1000
    if (ListStr.size() > 1) {
      std::vector<std::string> ListStrB=STRING_Split(ListStr[1], "_");
      if (ListStrB.size() == 2 && ListStrB[0] == "upperlimit")
	std::istringstream(ListStrB[1]) >> upperlimit;
    }
    ListIncd=lrs::DualDescription_temp_incd_limited(EXT, upperlimit);
    WeAreDone=true;
  }
  if (!WeAreDone) {
    std::cerr << "No right program found\n";
    std::cerr << "Let us die\n";
    throw TerminalException{1};
  }
  return ListIncd;
}



#endif

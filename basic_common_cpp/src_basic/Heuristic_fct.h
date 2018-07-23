#ifndef HEURISTIC_INCLUDE
#define HEURISTIC_INCLUDE

#include "Temp_common.h"

template<typename T>
struct SingleCondition {
  std::string eCond;
  std::string eType;
  T NumValue;
};


template<typename T>
struct OneFullCondition {
  std::vector<SingleCondition<T>> TheConditions;
  std::string TheResult;
};

template<typename T>
struct TheHeuristic {
  std::vector<OneFullCondition<T>> AllTests;
  std::string DefaultResult;
};



template<typename T>
TheHeuristic<T> ReadHeuristic(std::istream &is)
{
  if (!is.good()) {
    std::cerr << "ReadHeuristic operation failed because stream is not valied\n";
    throw TerminalException{1};
  }
  TheHeuristic<T> TheHeu;
  int nbFullCond;
  is >> nbFullCond;
  if (nbFullCond < 0) {
    std::cerr << "We must have nbFullCond >= 0\n";
    throw TerminalException{1};
  }
  for (int iFullCond=0; iFullCond<nbFullCond; iFullCond++) {
    std::vector<SingleCondition<T>> ListSingleCond;
    int nbCond;
    is >> nbCond;
    if (nbCond <= 0) {
      std::cerr << "Error, we must have nbCond > 0\n";
      std::cerr << "nbCond=" << nbCond << "\n";
      throw TerminalException{1};
    }
    for (int iCond=0; iCond<nbCond; iCond++) {
      std::string eType, eCond;
      T eNum;
      is >> eCond;
      is >> eType;
      is >> eNum;
      std::vector<std::string> ListType{">", "<", "=", "<=", ">="};
      bool IsMatch=false;
      for (auto & eTypePos : ListType) {
	if (eTypePos == eType)
	  IsMatch=true;
      }
      if (!IsMatch) {
	std::cerr << "Only allowed possibilities for eType are <, >, =, <= and >=\n";
	throw TerminalException{1};
      }
      if (eCond.size() == 0) {
	std::cerr << "eCond must be nontrivial otherwise evaluation will be impossible\n";
	throw TerminalException{1};
      }
      SingleCondition<T> eSingCond{eCond, eType, eNum};
      ListSingleCond.push_back(eSingCond);
    }
    std::string eResult;
    is >> eResult;
    if (eResult.size() == 0) {
      std::cerr << "eResult must be nontrivial otherwise evaluation will be impossible\n";
      throw TerminalException{1};
    }
    OneFullCondition<T> OneCond{ListSingleCond, eResult};
    TheHeu.AllTests.push_back(OneCond);
  }
  std::string DefaultResult;
  is >> DefaultResult;
  if (DefaultResult.size() == 0) {
    std::cerr << "DefaultResult must be nontrivial otherwise evaluation will be impossible\n";
    throw TerminalException{1};
  }
  TheHeu.DefaultResult=DefaultResult;
  return TheHeu;
}

template<typename T>
std::string HeuristicEvaluation(std::map<std::string, T> const& TheCand, TheHeuristic<T> const& TheHeu)
{
  for (auto const& eFullCond : TheHeu.AllTests) {
    bool IsOK=true;
    for (auto const& eSingCond : eFullCond.TheConditions) {
      std::string eCond=eSingCond.eCond;
      std::string eType=eSingCond.eType;
      T eNum=eSingCond.NumValue;
      auto search=TheCand.find(eCond);
      T eValue=0;
      if (search != TheCand.end()) {
	eValue=search->second;
      }
      else {
	std::cerr << "Entry " << eCond << " is required by heuristic\n";
	std::cerr << "Yet it is missing in the Candidate\n";
	std::cerr << "Please correct\n";
	throw TerminalException{1};
	//	exit(1);
      }
      bool WeMatch=false;
      if (eValue > eNum && eType == ">")
	WeMatch=true;
      if (eValue >= eNum && eType == ">=")
	WeMatch=true;
      if (eValue == eNum && eType == "=")
	WeMatch=true;
      if (eValue < eNum && eType == "<")
	WeMatch=true;
      if (eValue <= eNum && eType == "<=")
	WeMatch=true;
      if (!WeMatch)
	IsOK=false;
    }
    if (IsOK)
      return eFullCond.TheResult;
  }
  return TheHeu.DefaultResult;
}


template<typename T>
TheHeuristic<T> HeuristicFromListString(std::vector<std::string> const& ListString)
{
  size_t lenString = 30;
  std::string eRandString = random_string(lenString);
  std::string ePrefix = "/tmp/Std_adm";
  std::string TheFile=ePrefix + eRandString;
  std::ofstream OUTfs(TheFile);
  for (auto const& eStr : ListString)
    OUTfs << eStr << "\n";
  OUTfs.close();
  // Now reading it
  std::ifstream INfs(TheFile);
  TheHeuristic<T> TheHeu=ReadHeuristic<T>(INfs);
  std::remove(TheFile.c_str());
  return TheHeu;
}


#endif

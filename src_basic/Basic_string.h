#ifndef BASIC_STRING_INCLUDE
#define BASIC_STRING_INCLUDE

#include "Temp_common.h"
std::string STRING_GETENV(std::string const& eStr)
{
  //  std::string ePrefixAlti;
  //  std::cerr << "STRING_GETENV, step 1\n";
  //  char eStrEnvVar[]="ALTIMETER_DIRECTORY";
  //  std::cerr << "STRING_GETENV, step 2\n";
  char *ePre=std::getenv(eStr.c_str());
  //  std::cerr << "STRING_GETENV, step 3\n";
  if (ePre == NULL) {
    std::cerr << "Error in reading the environment variable : " << eStr << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "ePre=" << (void*)ePre << "\n";
  std::string eStrRet=ePre;
  //  std::cerr << "STRING_GETENV, step 4\n";
  return eStrRet;
}





bool STRING_IsStringReduceToSpace(std::string const& eStr)
{
  int len=eStr.length();
  std::string eChar=" ";
  for (int i=0; i<len; i++) {
    std::string eSubChar=eStr.substr(i,1);
    if (eSubChar != eChar)
      return false;
  }
  return true;
}



int STRING_GetCharPositionInString(std::string const& eStr, std::string const& eChar)
{
  int len=eStr.length();
  for (int i=0; i<len; i++) {
    std::string eSubChar=eStr.substr(i,1);
    if (eSubChar == eChar)
      return i;
  }
  return -1;
}



bool IsFullyNumeric(std::string const& eStr)
{
  std::string eLS="0123456789.";
  int nbChar=eStr.size();
  for (int uChar=0; uChar<nbChar; uChar++) {
    std::string eChar=eStr.substr(uChar,1);
    auto TheFind=[&](std::string const& eCharIn) -> bool {
      int nbCase=eLS.size();
      for (int iCase=0; iCase<nbCase; iCase++) {
	std::string eCase=eLS.substr(iCase,1);
	if (eCase == eCharIn)
	  return true;
      }
      return false;
    };
    if (TheFind(eChar) == false)
      return false;
  }
  return true;


}




std::string DoubleTo4dot2f(double const& x)
{
  char buffer[150];
  int n=sprintf(buffer, "%4.2f", x);
  if (n == 0) {
    std::cerr << "Clear error in DoubleTo4dot2f\n";
    throw TerminalException{1};
  }
  return std::string(buffer);
}


std::string DoubleToString(double const& x)
{
  std::stringstream s;
  s << std::fixed;
  s << std::setprecision(9);
  s << x;
  std::string converted(s.str());
  return converted;
}


std::string IntToString(int const & x)
{
  std::stringstream s;
  s << x;
  std::string converted(s.str());
  return converted;
}





std::string LongToString(long const & x)
{
  std::stringstream s;
  s << x;
  std::string converted(s.str());
  return converted;
}

int GetNumberDigit(int const& eVal)
{
  int nbDigit=1;
  while(true) {
    if (eVal < pow(10, nbDigit))
      return nbDigit;
    nbDigit++;
  }
}


std::string StringNumber(int const& nb, int const& nbDigit)
{
  if (nb > pow(10, nbDigit) ) {
    std::stringstream s;
    s << "Critical error in StringNumber\n";
    s << "nb=" << nb << "\n";
    s << "nbDigit=" << nbDigit << "\n";
    std::string eStr(s.str());
    throw eStr;
  }
  int idx=1;
  while(true) {
    if (nb < pow(10, idx) ) {
      std::string TheStr="";
      for (int i=0; i<nbDigit-idx; i++) {
	TheStr += "0";
      }
      TheStr += IntToString(nb);
      return TheStr;
    }
    idx++;
  }
}

std::string UpperCaseToLowerCase(std::string const& dataIn)
{
  std::string dataRet = dataIn;
  std::transform(dataRet.begin(), dataRet.end(), dataRet.begin(), ::tolower);
  return dataRet;
}




std::string STRING_RemoveSpacesBeginningEnd(std::string const& eStr)
{
  int len=eStr.size();
  std::vector<int> ListIsSpace(len,0);
  std::string eSpace=" ";
  for (int i=0; i<len; i++) {
    std::string eChar=eStr.substr(i, 1);
    if (eChar == eSpace)
      ListIsSpace[i]=1;
  }
  int PosLow=-1;
  for (int i=0; i<len; i++)
    if (PosLow == -1)
      if (ListIsSpace[i] == 0)
	PosLow=i;
  int PosUpp=-1;
  for (int i=0; i<len; i++) {
    int j=len-1-i;
    if (PosUpp == -1)
      if (ListIsSpace[j] == 0)
	PosUpp=j;
  }
  std::string RetStr;
  if (PosLow == -1) {
    return RetStr;
  }
  for (int iPos=PosLow; iPos<PosUpp+1; iPos++) {
    RetStr += eStr.at(iPos);
  }
  return RetStr;
}




// Example of use eStrA="A B C D" and eStrB=" "
std::vector<std::string> STRING_Split(std::string const& eStrA, std::string const& eStrB)
{
  int lenA=eStrA.length();
  int lenB=eStrB.length();
  std::vector<int> ListStatus(lenA,1);
  for (int iA=0; iA<lenA - lenB; iA++)
    if (ListStatus[iA] == 1) {
      bool IsMatch=true;
      for (int iB=0; iB<lenB; iB++) {
	std::string eCharA=eStrA.substr(iA+iB,1);
	std::string eCharB=eStrB.substr(iB,1);
	if (eCharA != eCharB)
	  IsMatch=false;
      }
      if (IsMatch)
	for (int iB=0; iB<lenB; iB++)
	  ListStatus[iA + iB]=0;
    }
  std::vector<std::string> RetList;
  std::string eFound;
  for (int iA=0; iA<lenA; iA++) {
    std::string eChar=eStrA.substr(iA, 1);
    if (ListStatus[iA] == 1)
      eFound += eChar;
    if (ListStatus[iA] == 0) {
      int siz=eFound.length();
      if (siz > 0)
	RetList.push_back(eFound);
      eFound="";
    }
  }
  int siz=eFound.size();
  if (siz > 0)
    RetList.push_back(eFound);
  return RetList;
}


std::vector<int> STRING_Split_Int(std::string const& eStrA, std::string const& eStrB)
{
  std::vector<std::string> LStr=STRING_Split(eStrA, eStrB);
  std::vector<int> LInt;
  for (auto & eStr : LStr) {
    int eVal;
    std::istringstream(eStr) >> eVal;
    LInt.push_back(eVal);
  }
  return LInt;
}



// Example of use eStrA="A B C D" and eStrB=" "
// Difference with the above is that splitting ";;" gives you a list
// of entries as {"", "", ""}, i.e. last and first entry and entry in the middle
std::vector<std::string> STRING_Split_Strict(std::string const& eStrA, std::string const& eStrB)
{
  int lenA=eStrA.length();
  int lenB=eStrB.length();
  std::vector<int> ListStatus(lenA,0);
  int idx=0;
  for (int iA=0; iA<lenA - lenB; iA++) {
    int sumEnt=0;
    for (int iB=0; iB<lenB; iB++)
      sumEnt += ListStatus[iA+iB];
    if (sumEnt == 0) {
      bool IsMatch=true;
      for (int iB=0; iB<lenB; iB++) {
	std::string eCharA=eStrA.substr(iA+iB,1);
	std::string eCharB=eStrB.substr(iB,1);
	if (eCharA != eCharB)
	  IsMatch=false;
      }
      if (IsMatch) {
	idx++;
	for (int iB=0; iB<lenB; iB++)
	  ListStatus[iA + iB]=idx;
      }
    }
  }
  int nbEnt=idx+1;
  std::vector<std::string> RetList(nbEnt);
  for (int iEnt=0; iEnt<nbEnt; iEnt++) {
    int posFirst=-1, posLast=-1;
    if (iEnt == 0) {
      posFirst=0;
    }
    else {
      bool WeFound=false;
      for (int iA=0; iA<lenA; iA++)
	if (!WeFound)
	  if (ListStatus[iA] == iEnt) {
	    WeFound=true;
	    posFirst=iA+lenB;
	  }
    }
    if (iEnt == nbEnt-1) {
      posLast=lenA-1;
    }
    else {
      bool WeFound=false;
      for (int iA=0; iA<lenA; iA++)
	if (!WeFound)
	  if (ListStatus[iA] == iEnt+1) {
	    WeFound=true;
	    posLast=iA-1;
	  }
    }
    if (posFirst == -1 || posLast == -1) {
      std::cerr << "posFirst = " << posFirst << "  posLast = " << posLast << "\n";
      std::cerr << "Positions have not been found\n";
      throw TerminalException{1};
    }
    std::string str;
    for (int i=posFirst; i<=posLast; i++) {
      std::string eChar=eStrA.substr(i,1);
      str += eChar;
    }
    RetList[iEnt]=str;
  }
  return RetList;
}


std::string STRING_Replace(std::string const&eStrA, std::string const& eStrB, std::string const& eStrC)
{
  //  std::cerr << "Before STRING_Split_Strict eStrA=" << eStrA << " eStrB=" << eStrB << "\n";
  std::vector<std::string> LStr=STRING_Split(eStrA, eStrB);
  //  std::cerr << " After STRING_Split_Strict\n";
  std::string str=LStr[0];
  int len=LStr.size();
  //  std::cerr << " len=" << len << "\n";
  for (int i=1; i<len; i++)
    str += eStrC + LStr[i];
  return str;
}





std::vector<std::string> STRING_SplitCharNb(std::string const& str)
{
  auto IsNumber=[](std::string const& eChar) {
    std::vector<std::string> ListCharNb{"-","0","1","2","3","4","5","6","7","8","9"};
    for (auto & fCharNb : ListCharNb)
      if (fCharNb == eChar)
	return true;
    return false;
  };
  int len=str.size();
  std::vector<bool> ListStat(len);
  for (int i=0; i<len; i++) {
    std::string eChar=str.substr(i,1);
    ListStat[i]=IsNumber(eChar);
  }
  std::string TotStr;
  for (int i=0; i<len; i++) {
    TotStr += str.substr(i,1);
    if (i<len-1) {
      if (ListStat[i] != ListStat[i+1]) {
	TotStr += "WRK";
      }
    }
  }
  //  std::cerr << "TotStr=" << TotStr << "\n";
  return STRING_Split(TotStr, "WRK");
}



std::string FILE_GetExtension(std::string const& eFile)
{
  std::vector<std::string> LStr=STRING_Split(eFile, "/");
  std::string eFinal=LStr[LStr.size()-1];
  std::vector<std::string> LBlck=STRING_Split(eFile, ".");
  return LBlck[LBlck.size()-1];
}


#endif

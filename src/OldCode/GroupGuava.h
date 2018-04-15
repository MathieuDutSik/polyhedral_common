#ifndef TEMP_GROUP_GUAVA_FUNCTIONS
#define TEMP_GROUP_GUAVA_FUNCTIONS

/*
  This code contains an attempt to access to the
  GUAVA program authored by Jeffrey Leon.
  The problem is that this code is deeply broken
  (memory leaks, etc.) and so should not be used.
  The right way seems to be to use the code authored
  by Thomas Rehn.
*/

#include "Temp_common.h"

void GUAVA_WriteGroup(std::string const& eFile, TheGroupFormat const& eGRP)
{
  int n=eGRP.n;
  std::string eFileTot="/tmp/" + eFile;
  std::ofstream os(eFileTot);
  std::vector<mpz_class> LFact_z=FactorsInt(eGRP.size);
  os << "LIBRARY " + eFile + ";\n";
  os << "\"group input for setStab and others\"\n";
  os << eFile + ": permutation group(" << n << ");\n";
  os << eFile + ".forder:";
  std::set<mpz_class> Set_z;
  for (auto& eVal_z : LFact_z)
    Set_z.insert(eVal_z);
  bool IsFirst=true;
  for (auto& ePrimeFact_z : Set_z) {
    int eMult=0;
    for (auto & eVal_z : LFact_z)
      if (eVal_z == ePrimeFact_z)
	eMult++;
    if (IsFirst == false) {
      os << " *";
    }
    IsFirst=false;
    if (eMult == 1) {
      os << " " << ePrimeFact_z;
    }
    else {
      os << " " << ePrimeFact_z << "^" << eMult;
    }
  }
  os << ";\n";
  //
  std::vector<permlib::dom_int> TheBase=GetBaseGroup(eGRP);
  os << eFile + ".base:  seq(";
  IsFirst=true;
  for (auto & eVal : TheBase) {
    if (IsFirst == false)
      os << ",";
    IsFirst=false;
    int eValP=eVal+1;
    os << eValP;
  }
  os << ");\n";
  os << eFile + ".strong generators:  [\n";
  int idx=0;
  IsFirst=true;
  for (auto & ePerm : eGRP.group->S) {
    idx++;
    if (IsFirst == false)
      os << ",\n";
    IsFirst=false;
    os << "x" << idx << " = " << *ePerm;
  }
  os << "];\n";
  os << "FINISH;\n";
}

void GUAVA_WritePointSet(std::string const& eFile, Face const& f)
{
  int nb=f.count();
  int siz=f.size();
  std::string eFileTot="/tmp/" + eFile;
  std::ofstream os(eFileTot);
  os << "LIBRARY " + eFile + ";\n";
  os << "\" A subset of 1,...," << siz << " of size " << nb << ".\"\n";
  os << eFile + " = [";
  int aRow=f.find_first();
  for (int i=0; i<nb; i++) {
    if (i>0)
      os << ",";
    int eVal=aRow+1;
    os << eVal;
    aRow=f.find_next(aRow);    
  }
  os << "];\n";
  os << "FINISH;\n";
}

TheGroupFormat GUAVA_ReadGroup(std::string const& eFile)
{
  std::string eFileTot="/tmp/" + eFile;
  std::vector<std::string> ListLines=ReadFullFile(eFile);
  int nbLine=ListLines.size();
  //
  int iLinestrong=-1;
  int iLinefinish=-1;
  int iLinesize=-1;
  for (int iLine=0; iLine<nbLine; iLine++) {
    std::string eLine=ListLines[iLine];
    std::vector<std::string> LSpltA=STRING_Split(eLine, "strong");
    if (LSpltA.size() > 1) {
      if (iLinestrong != -1) {
	std::cerr << "Error, two lines with strong in their expression\n";
	exit(1);
      }
      iLinestrong=iLine;
    }
    std::vector<std::string> LSpltB=STRING_Split(eLine, "INISH");
    if (LSpltB.size() > 1) {
      if (iLinefinish != -1) {
	std::cerr << "Error, two lines with INISH in their expression\n";
	exit(1);
      }
      iLinefinish=iLine;
    }
    std::vector<std::string> LSpltC=STRING_Split(eLine, "Permutation group");
    if (LSpltC.size() > 1) {
      if (iLinesize != -1) {
	std::cerr << "Error, two lines with Permutation group in their expression\n";
	exit(1);
      }
      iLinesize=iLine;
    }
  }
  //  std::cerr << "iLinestrong=" << iLinestrong << "\n";
  //  std::cerr << "iLinefinish=" << iLinefinish << "\n";
  //  std::cerr << "iLinesize=" << iLinesize << "\n";
  if (iLinestrong == -1 || iLinefinish == -1 || iLinesize == -1) {
    std::cerr << "Error with unset values\n";
    exit(1);
  }
  std::vector<std::string> LSplSizA=STRING_Split(ListLines[iLinesize], "Permutation group (");
  std::string eStrSizA=LSplSizA[1];
  std::vector<std::string> LSplSizB=STRING_Split(eStrSizA, ")");
  std::string eStrSizB=LSplSizB[0];
  int nbVert=atoi(eStrSizB.c_str());
  //  std::cerr << "nbVert=" << nbVert << "\n";
  //
  // Concatenating all the lines together for the relevant information
  //
  std::string eStringConcat="";
  for (int iLine=iLinestrong+1; iLine<iLinefinish; iLine++) {
    std::string eLine=ListLines[iLine];
    std::string eLineAppend;
    if (iLine == iLinefinish-1) {
      std::string eLastChar="]";
      int siz=eLine.size();
      for (int i=0; i<siz; i++) {
	std::string eSubStr=eLine.substr(i,1);
	if (eSubStr == "]")
	  eLineAppend=eLine.substr(0,i);
      }
    }
    else {
      eLineAppend=eLine;
    }
    eStringConcat += eLineAppend;
  }
  //
  // Determining the , separating the generators
  //
  int nbTotalLen=eStringConcat.size();
  std::vector<int> ListPos;
  int TheLevel=0;
  for (int i=0; i<nbTotalLen; i++) {
    std::string eChar=eStringConcat.substr(i,1);
    if (eChar == "(")
      TheLevel++;
    if (eChar == ")")
      TheLevel--;
    if (TheLevel != 0 && TheLevel != 1) {
      std::cerr << "Inconsistency in reading generator\n";
      exit(1);
    }
    if (TheLevel == 0 && eChar == ",")
      ListPos.push_back(i);
  }
  if (TheLevel != 0) {
    std::cerr << "Error, we should have TheLevel=0\n";
    exit(1);
  }
  //
  // separating the generators.
  //
  std::vector<std::string> ListStringGenerator;
  int nbGen=ListPos.size()+1;
  //  std::cerr << "nbGen=" << nbGen << "\n";
  for (int iGen=0; iGen<nbGen; iGen++) {
    int ePosFirst, ePosLast, len;
    if (iGen == 0) {
      ePosFirst=0;
    }
    else {
      ePosFirst=ListPos[iGen-1]+1;
    }
    if (iGen == nbGen-1) {
      ePosLast=nbTotalLen-1;
    }
    else {
      ePosLast=ListPos[iGen]-1;
    }
    len=ePosLast - ePosFirst + 1;
    std::string eStringGen=eStringConcat.substr(ePosFirst,len);
    ListStringGenerator.push_back(eStringGen);
  }
  //
  // Creating the permutation
  //
  std::vector<permlib::Permutation> ListPerm;
  for (auto & eStringGen : ListStringGenerator) {
    //    std::cerr << "eStringGen=" << eStringGen << "\n";
    std::vector<permlib::dom_int> v(nbVert);
    for (int iVert=0; iVert<nbVert; iVert++)
      v[iVert]=iVert;
    std::vector<int> ListPosParenFirst, ListPosParenLast;
    int sizP=eStringGen.size();
    for (int j=0; j<sizP; j++) {
      std::string eChar=eStringGen.substr(j,1);
      if (eChar == "(")
	ListPosParenFirst.push_back(j);
      if (eChar == ")")
	ListPosParenLast.push_back(j);
    }
    int nbCycle=ListPosParenFirst.size();
    if (size_t(nbCycle) != ListPosParenLast.size()) {
      std::cerr << "Error on the parenthesis counting\n";
      exit(1);
    }
    //    std::cerr << "  nbCycle=" << nbCycle << "\n";
    for (int iCycle=0; iCycle<nbCycle; iCycle++) {
      int ePosFirst=ListPosParenFirst[iCycle];
      int ePosLast=ListPosParenLast[iCycle];
      if (ePosFirst >= ePosLast) {
	std::cerr << "Error on the ePosFirst / ePosLast\n";
	exit(1);
      }
      int ePos=ePosFirst+1;
      int len=ePosLast - ePosFirst-1;
      std::string eStrPart=eStringGen.substr(ePos,len);
      std::vector<std::string> LSplGen=STRING_Split(eStrPart, ",");
      int lenCycle=LSplGen.size();
      //      std::cerr << "    iCycle=" << iCycle << " lenCycle=" << lenCycle << "\n";
      std::vector<int> ListVal(lenCycle);
      for (int i=0; i<lenCycle; i++) {
	int eVal=atoi(LSplGen[i].c_str());
	ListVal[i]=eVal-1;
      }
      for (int i=0; i<lenCycle; i++) {
	int iNext=NextIdx(lenCycle, i);
	v[ListVal[i]]=ListVal[iNext];
      }
    }
    permlib::Permutation ePerm(v);
    ListPerm.push_back(ePerm);
  }
  return GetPermutationGroup(nbVert, ListPerm);
}


TheGroupFormat GUAVA_GetStabilizer(TheGroupFormat const& TheGRP, Face const& eList)
{
  std::string eStrRand=random_string_restricted(20);
  std::string eFileGRP="grp_" + eStrRand;
  std::string eFileSET="set_" + eStrRand;
  std::string eFileSTAB="stab_" + eStrRand;
  std::cerr << "eFileGRP = " << eFileGRP << "\n";
  std::cerr << "eFileSET = " << eFileSET << "\n";
  std::string eFileERR="/tmp/err_" + eStrRand;
  std::string eFileERR2="/tmp/err2_" + eStrRand;
  GUAVA_WriteGroup(eFileGRP, TheGRP);
  GUAVA_WritePointSet(eFileSET, eList);
  /*
  std::string eFileGRPgap="/tmp/GRPgap_" + eStrRand;
  std::string eFileSETgap="/tmp/SETgap_" + eStrRand;
  std::cerr << "eFileGRPgap = " << eFileGRPgap << "\n";
  std::cerr << "eFileSETgap = " << eFileSETgap << "\n";
  std::ofstream os1(eFileGRPgap);
  WriteGroupGAP(os1, TheGRP);
  os1.close();
  std::ofstream os2(eFileSETgap);
  os2 << "return ";
  WriteFaceGAP(os2, eList);
  os2 << ";\n";
  os2.close();
  */
  std::cerr << "Files have been created\n";
  std::string order="(cd /tmp && setstab " + eFileGRP + " " + eFileSET + " " + eFileSTAB + " > " + eFileERR + " 2> " + eFileERR2 + ")";
  int iret=system(order.c_str());
  if (iret == -1) {
    std::cerr << "Error in running setStab from guava\n";
    exit(1);
  }
  std::string eFileGRP_tot="/tmp/" + eFileGRP;
  std::string eFileSET_tot="/tmp/" + eFileSET;
  std::string eFileSTAB_tot="/tmp/" + eFileSTAB;
  //  std::cerr << "eFileGRP_tot = " << eFileGRP_tot << "\n";
  //  std::cerr << "eFileSET_tot = " << eFileSET_tot << "\n";
  //  std::cerr << "eFileSTAB_tot = " << eFileSTAB_tot << "\n";
  TheGroupFormat eStab=GUAVA_ReadGroup(eFileSTAB_tot);
  //  std::cerr << "After GUAVA_ReadGroup\n";
  RemoveFileIfExist(eFileGRP_tot);
  RemoveFileIfExist(eFileSET_tot);
  RemoveFileIfExist(eFileSTAB_tot);
  RemoveFileIfExist(eFileERR);
  RemoveFileIfExist(eFileERR2);
  return eStab;
}



#endif

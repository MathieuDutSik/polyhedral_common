#ifndef BASIC_FILE_INCLUDE
#define BASIC_FILE_INCLUDE

#include "Temp_common.h"
#include "Basic_string.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>


void CopyOperation(std::string const& SrcFile, std::string const& DstFile)
{
  std::string eComm="cp " + SrcFile + " " + DstFile;
  int iret=system(eComm.c_str() );
  if (iret != 0) {
    std::cerr << "Error in copy operation\n";
    std::cerr << "SrcFile=" << SrcFile << "\n";
    std::cerr << "DstFile=" << DstFile << "\n";
    throw TerminalException{1};
  }
}


bool IsExistingFile(std::string const& eFile)
{
  std::ifstream f(eFile.c_str());
  if (f.good()) {
    f.close();
    return true;
  } else {
    f.close();
    return false;
  }   
}

void IsExistingFileDie(std::string const& eFile)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "The file eFile = " << eFile << "\n";
    std::cerr << "Is missing. DIE\n";
    throw TerminalException{1};
  }
}





std::vector<std::string> ReadFullFile(std::string const& eFile)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "ReadFullFile eFile=" << eFile << "\n";
    std::cerr << "Missing file\n";
    throw TerminalException{1};
  }
  std::ifstream is(eFile);
  std::string line;
  std::vector<std::string> ListLines;
  while (std::getline(is, line))
    ListLines.push_back(line);
  return ListLines;
}



std::string FILE_GetNakedFilename(std::string const& eFileFull)
{
  std::vector<std::string> LStr=STRING_Split(eFileFull, "/");
  std::string eLast=LStr[LStr.size() -1];
  return eLast;
}



std::vector<std::string> FILE_GetDirectoryListFile(std::string const& eDir)
{
  std::string ePath=eDir + ".";
  DIR* dirp=opendir(ePath.c_str());
  if (dirp == NULL) {
    std::cerr << "Error in routine FILE_GetDirectoryListFile\n";
    std::cerr << "Error in call to opendir\n";
    std::cerr << "eDir = " << eDir << "\n";
    throw TerminalException{1};
  }
  struct dirent *dp;
  std::vector<std::string> ListFile;
  while ((dp = readdir(dirp)) != NULL) {
    std::string eName=dp->d_name;
    //    free(dp); // not sure this is portable
    if (eName != ".." && eName != ".")
      ListFile.push_back(eName);
  }
  int err=closedir(dirp);
  if (err != 0) {
    std::cerr << "err=" << err << "\n";
    printf("Oh dear, something went wrong with ls! %s\n", strerror(errno));
    throw TerminalException{1};
  }
  return ListFile;
}

bool FILE_IsDirectoryEmpty(std::string const& eDir)
{
  std::vector<std::string> TheList = FILE_GetDirectoryListFile(eDir);
  if (TheList.size() == 0)
    return true;
  return false;
}



bool FILE_CheckFinalShashDirectory(std::string const& eDir)
{
  int len=eDir.size();
  std::string eChar=eDir.substr(len-1,1);
  if (eChar == "/")
    return true;
  return false;
}



bool FILE_IsRegularFile(std::string const& eFile)
{
  int status;
  struct stat st_buf;  
  status = stat(eFile.c_str(), &st_buf);
  if (status != 0) {
    std::cerr << "Problem in FILE_IsRegularFile\n";
    std::cerr << "Error, errno = " << errno << "\n";
    std::cerr << "eFile=" << eFile << "\n";
    throw TerminalException{1};
  }
  if (S_ISREG (st_buf.st_mode)) {
    return true;
  }
  return false;
}



std::vector<std::string> FILE_GetDirectoryFilesRecursively(std::string const& eDir)
{
  //  std::cerr << "Beginning of FILE_GetDirectoryFilesRecursively\n";
  std::vector<std::string> ListDir={eDir};
  std::vector<std::string> ListFile;
  while(true) {
    std::vector<std::string> NewListDir;
    for (auto & fDir : ListDir) {
      std::vector<std::string> LocalListFile=FILE_GetDirectoryListFile(fDir);
      for (auto & eFile : LocalListFile) {
	std::string NewEnt=fDir + eFile;
	if (FILE_IsRegularFile(NewEnt)) {
	  ListFile.push_back(NewEnt);
	}
	else {
	  std::string NewDir=NewEnt + "/";
	  NewListDir.push_back(NewDir);
	}
      }
    }
    if (NewListDir.size() == 0)
      break;
    ListDir=NewListDir;
  }
  //  std::cerr << "Ending of FILE_GetDirectoryFilesRecursively\n";
  return ListFile;
}

std::vector<std::string> FILE_DirectoryFilesSpecificExtension(std::string const& ePrefix, std::string const& eExtension)
{
  std::vector<std::string> ListFile=FILE_GetDirectoryFilesRecursively(ePrefix);
  std::vector<std::string> RetListFile;
  int lenExt=eExtension.size();
  for (auto & eFile : ListFile) {
    int len=eFile.size();
    if (len > lenExt) {
      std::string eEnd=eFile.substr(len-lenExt,lenExt);
      if (eEnd == eExtension)
	RetListFile.push_back(eFile);
    }
  }
  return RetListFile;
}





bool IsExistingDirectory(std::string const& ThePrefix)
{
  if (0 != access(ThePrefix.c_str(), F_OK)) {
    if (ENOENT == errno) {
      // does not exist
      return false;
    }
    if (ENOTDIR == errno) {
      return false;
      // not a directory
    }
    std::cerr << "Should not happen a priori\n";
    throw TerminalException{1};
  }
  return true;
}


void RemoveEmptyDirectory(std::string const& eDir)
{
  std::string eOrder="rm -d " + eDir;
  int iret=system(eOrder.c_str() );
  if (iret == -1) {
    std::cerr << "Error in RemoveEmptyDirectory\n";
    std::cerr << "eDir=" << eDir << "\n";
    std::cerr << "eOrder=" << eOrder << "\n";
    std::cerr << "unable to run the process\n";
    throw TerminalException{1};
  }
}




void RemoveFile(std::string const& eFile)
{
  std::remove(eFile.c_str());
}

void RemoveFileIfExist(std::string const& eFile)
{
  if (IsExistingFile(eFile))
    RemoveFile(eFile);
}


std::string FILE_RemoveEndingExtension(std::string const& FileName, std::string const& TheExtension)
{
  int len=FileName.size();
  int iCharLast=-1;
  for (int iChar=0; iChar<len; iChar++) {
    std::string eChar=FileName.substr(iChar,1);
    if (eChar == ".")
      iCharLast=iChar;
  }
  if (iCharLast == -1)
    return FileName;
  std::string FileNameRed=FileName.substr(0,iCharLast);
  std::string eExtension=FileName.substr(iCharLast+1,len-1-iCharLast);
  if (eExtension == TheExtension)
    return FileNameRed;
  return FileName;
}


void RemoveFileSpecificExtension(std::string const& ThePrefix, std::string const& TheExtension)
{
  bool test=IsExistingDirectory(ThePrefix);
  if (!test)
    return;
  std::vector<std::string> ListFile=FILE_GetDirectoryListFile(ThePrefix);
  int nbCharEnd=TheExtension.size();
  for (auto & eFile : ListFile) {
    int len=eFile.size();
    if (len > nbCharEnd) {
      std::string TheEnd=eFile.substr(len-nbCharEnd, nbCharEnd);
      if (TheEnd == TheExtension) {
	std::string eFileTot=ThePrefix + eFile;
	RemoveFile(eFileTot);
      }
    }
  }
}





int FILE_GetNumberLine(std::string const& eFile)
{
  int number_of_lines = 0;
  std::string line;
  std::ifstream myfile(eFile);
  while (std::getline(myfile, line))
    ++number_of_lines;
  return number_of_lines;
}




#ifndef WINDOWS
std::string GetCurrentDirectory()
{
  int size = pathconf(".", _PC_PATH_MAX);
  char *buf;
  char *ptr;
  if ((buf = (char *)malloc((size_t)size)) != NULL) {
    ptr=getcwd(buf, (size_t)size);
    if (ptr == NULL && errno != ERANGE) {
      std::cerr << "Error while trying to use getcwd\n";
      throw TerminalException{1};
    }
    std::string eRet = buf;
    eRet=eRet + "/";
    free(buf);
    if (ptr != NULL) {
      if (ptr != buf) {
	std::cerr << "Before ptr freeing\n";
	free(ptr);
	std::cerr << "After ptr freeing\n";
      }
    }
    return eRet;
  }
  else {
    std::cerr << "Not enough memory\n";
    throw TerminalException{1};
  }
}
#endif


#ifndef WINDOWS
std::string FILE_GetAbsoluteDirectory(std::string const& ePrefix)
{
  std::string FirstChar=ePrefix.substr(0, 1);
  if (FirstChar == "/") {
    return ePrefix;
  }
  else {
    std::string ePWD=GetCurrentDirectory();
    return ePWD + ePrefix;
  }
}
#endif

bool FILE_CheckPrefix(std::string const& ePrefix)
{
  int len;
  len=ePrefix.size();
  std::string LastChar=ePrefix.substr(len-1,1);
  //  std::cerr << "LastChar=" << LastChar << "\n";
  if (LastChar == "/")
    return true;
  return false;
}

std::string ExtractDirectoryFromFileString(std::string const& eFile)
{
  int len=eFile.size();
  int iCharFinal=-1;
  for (int iChar=0; iChar<len; iChar++) {
    std::string eChar=eFile.substr(iChar,1);
    if (eChar == "/")
      iCharFinal=iChar;
  }
  if (iCharFinal == -1) {
    std::cerr << "Error in ExtractDirectoryFromFileString\n";
    throw TerminalException{1};
  }
  return eFile.substr(0,iCharFinal+1);
}


bool FILE_IsFileMakeable(std::string const& eFile)
{
  std::string eDir=ExtractDirectoryFromFileString(eFile);
  if (!IsExistingFile(eDir))
    return false;
  return true;
}





void CreateDirectory(std::string const& eDir)
{
  //  std::cerr << "eDir=" << eDir << "\n";
  const char *dir=eDir.c_str();
  char tmp[256];
  char *p = NULL;
  size_t len;
  snprintf(tmp, sizeof(tmp),"%s",dir);
  len = strlen(tmp);
  if(tmp[len - 1] == '/')
    tmp[len - 1] = 0;
  for(p = tmp + 1; *p; p++)
    if(*p == '/') {
      *p = 0;
      //      std::cerr << "Before mkdir, tmp=" << tmp << "\n";
      mkdir(tmp, S_IRWXU);
      *p = '/';
    }
  //  std::cerr << "Before mkdir, tmp=" << tmp << "\n";
  mkdir(tmp, S_IRWXU);
}


void CreateDirectory_V1(std::string const& eDir)
{
  std::string eOrder="/bin/mkdir -p " + eDir;
  int iret=system(eOrder.c_str() );
  if (iret == -1) {
    std::cerr << "Error in CreateDirectory\n";
    std::cerr << "eDir=" << eDir << "\n";
    std::cerr << "eOrder=" << eOrder << "\n";
    std::cerr << "unable to run the process\n";
    throw TerminalException{1};
  }
}


std::vector<std::string> ls_operation(std::string const& ThePrefix)
{
  std::cerr << "Doing ls_operation\n";
  std::string strRand=random_string(20);
  std::string TmpFile="/tmp/file" + strRand;
  std::string ErrFile="/tmp/file" + strRand + ".err";
  std::string eOrder="ls " + ThePrefix + " > " + TmpFile + " 2> " + ErrFile;
  int iret=system(eOrder.c_str() );
  std::cerr << "iret=" << iret << "\n";
  if (iret != -1) {
    std::cerr << "Error in ls_operation\n";
    std::cerr << "ThePrefix=" << ThePrefix << "\n";
    std::cerr << "unable to run the process\n";
    throw TerminalException{1};
  }
  std::ifstream iserr(ErrFile);
  int nbCharErr=0;
  while(!iserr.eof()) {
    std::string PreStr;
    std::getline(iserr, PreStr);
    nbCharErr += PreStr.size();
  }
  if (nbCharErr > 0) {
    std::cerr << "Error in ls_operation\n";
    std::cerr << "We have nbCharErr = " << nbCharErr << "\n";
    std::cerr << "TmpFile=" << TmpFile << "\n";
    std::cerr << "ErrFile=" << ErrFile << "\n";
    std::cerr << "ThePrefix=" << ThePrefix << "\n";
    throw TerminalException{1};
  }
  //
  std::ifstream is(TmpFile);
  std::vector<std::string> ListFile;
  while(true) {
    if (is.eof()) {
      is.close();
      RemoveFile(TmpFile);
      RemoveFile(ErrFile);
      return ListFile;
    }
    std::string eFile;
    is >> eFile;
    ListFile.push_back(eFile);
  }
}



struct TempFile {
private:
  std::string FileName;
public:
  TempFile() = delete;
  TempFile(std::string const& eFile)
  {
    FileName=eFile;
  }
  TempFile(char* eFile)
  {
    FileName=eFile;
  }
  ~TempFile()
  {
    if (IsExistingFile(FileName)) {
      RemoveFile(FileName);
    }
  }
  //
  bool IsExisting() const
  {
    return IsExistingFile(FileName);
  }
  std::string string() const
  {
    return FileName;
  }
};



// This does not provide all the facility that we are after
// The problem is that when we have an exit(1)
// the destructors of the variables are not called.
//
// This means that we need to put a throw statement.
// in order to have the temporary directory eliminated.
// So eliminate all the exit(1) statement and replace by throw
// and catch the throws.
struct TempDirectory {
private:
  std::string DirName;
public:
  TempDirectory()
  {
    DirName="unset_and_irrelevant";
  }
  TempDirectory(std::string const& eDir)
  {
    DirName=eDir;
    CreateDirectory(DirName);
  }
  TempDirectory& operator=(TempDirectory&& eTemp)
  {
    DirName=eTemp.str();
    eTemp.DirName="unset_and_irrelevant";
    return *this;
  }
  TempDirectory(TempDirectory && eTemp) : DirName(eTemp.str())
  {
  }
  ~TempDirectory()
  {
    //    std::cerr << "Calling destructor\n";
    if (DirName != "unset_and_irrelevant") {
      //      std::cerr << "  Destructor is really needed\n";
      //      std::cerr << "  DirName=" << DirName << "\n";
      if (IsExistingDirectory(DirName)) {
	//	std::cerr << "    Directory really exist\n";
	if (!FILE_IsDirectoryEmpty(DirName)) {
	  std::cerr << "Keeping " << DirName << " since it is not empty\n";
	}
	else {
	  RemoveFile(DirName);
	}
      }
    }
  }
  //
  bool IsExisting() const
  {
    return IsExistingDirectory(DirName);
  }
  std::string str() const
  {
    return DirName;
  }
};


struct CondTempDirectory {
private:
  bool used;
  std::string DirName;
public:
  CondTempDirectory()
  {
    used=false;
    DirName="unset_and_irrelevant";
  }
  CondTempDirectory(bool const& eUsed, std::string const& eDir)
  {
    used=eUsed;
    if (used) {
      DirName=eDir;
      CreateDirectory(DirName);
    }
    else {
      DirName="unset_and_irrelevant";
    }
  }
  CondTempDirectory& operator=(CondTempDirectory&& eTemp)
  {
    used=eTemp.usedness();
    DirName=eTemp.str();
    eTemp.DirName="unset_and_irrelevant";
    return *this;
  }
  CondTempDirectory(CondTempDirectory && eTemp) : used(eTemp.usedness()), DirName(eTemp.str())
  {
  }
  ~CondTempDirectory()
  {
    //    std::cerr << "Calling destructor\n";
    if (used && DirName != "unset_and_irrelevant") {
      //      std::cerr << "  Destructor is really needed\n";
      //      std::cerr << "  DirName=" << DirName << "\n";
      if (IsExistingDirectory(DirName)) {
	//	std::cerr << "    Directory really exist\n";
	if (!FILE_IsDirectoryEmpty(DirName)) {
	  std::cerr << "Keeping " << DirName << " since it is not empty\n";
	}
	else {
	  RemoveFile(DirName);
	}
      }
    }
  }
  //
  bool IsExisting() const
  {
    return IsExistingDirectory(DirName);
  }
  bool usedness() const
  {
    return used;
  }
  std::string str() const
  {
    return DirName;
  }
};



#endif

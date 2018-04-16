#ifndef INCLUDE_FACE_BITSET
#define INCLUDE_FACE_BITSET

// Boost libraries

#include "Temp_common.h"
#include <bitset>
#include <boost/dynamic_bitset.hpp>

typedef boost::dynamic_bitset<> Face;

std::vector<int> FaceToVector(Face const& eSet)
{
  int nbVert=eSet.count();
  std::vector<int> eList(nbVert);
  int aRow=eSet.find_first();
  for (int i=0; i<nbVert; i++) {
    eList[i]=aRow;
    aRow=eSet.find_next(aRow);
  }
  return eList;
}


std::vector<int> FaceTo01vector(Face const& eSet)
{
  int nbVert=eSet.size();
  int siz=eSet.count();
  std::vector<int> eList(nbVert,0);
  int aRow=eSet.find_first();
  for (int i=0; i<siz; i++) {
    eList[aRow]=1;
    aRow=eSet.find_next(aRow);
  }
  return eList;
}







void WriteFace(std::ostream & os, Face const& eList)
{
  int len;
  len=eList.size();
  os << len;
  for (int i=0; i<len; i++) {
    int eVal=eList[i];
    os << " " << eVal;
  }
  os << "\n";
}

Face ReadFace(std::istream & is)
{
  if (!is.good()) {
    std::cerr << "ReadFace operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  int len, eVal;
  is >> len;
  Face eFace(len);
  for (int i=0; i<len; i++) {
    is >> eVal;
    eFace[i]=eVal;
  }
  return eFace;
}


std::vector<Face> ReadListFace(std::istream & is)
{
  if (!is.good()) {
    std::cerr << "ReadListFace operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  int nbFace;
  is >> nbFace;
  std::vector<Face> ListFace(nbFace);
  for (int iFace=0; iFace<nbFace; iFace++)
    ListFace[iFace]=ReadFace(is);
  return ListFace;
}

void WriteListFace(std::ostream & os, std::vector<Face> const& ListFace)
{
  int nbFace=ListFace.size();
  os << nbFace << "\n";
  for (int iFace=0; iFace<nbFace; iFace++)
    WriteFace(os, ListFace[iFace]);
}


void WriteFaceGAP(std::ostream &os, Face const& f)
{
  int nb=f.count();
  //  int siz=f.size();
  os << "[";
  int aPos=f.find_first();
  for (int i=0; i<nb; i++) {
    if (i>0)
      os << ",";
    int eVal=aPos+1;
    os << eVal;
    aPos=f.find_next(aPos);
  }
  os << "]";
}


void WriteListFaceGAP(std::ostream & os, std::vector<Face> const& ListFace)
{
  os << "[";
  bool IsFirst=true;
  for (auto & eFace : ListFace) {
    if (!IsFirst)
      os << ",";
    IsFirst=false;
    WriteFaceGAP(os, eFace);
  }
  os << "]";
}

void WriteListFaceGAPfile(std::string const& eFile, std::vector<Face> const& ListFace)
{
  std::ofstream os(eFile);
  os << "return ";
  WriteListFaceGAP(os, ListFace);
  os << ";\n";
}




// We require x and y to be of the same size
bool operator<(Face const& x, Face const& y)
{
  int len=x.size();
  for (int i=0; i<len; i++) {
    if (x[i] == 0 && y[i] == 1)
      return true;
    if (x[i] == 1 && y[i] == 0)
      return false;
  }
  return false;
}




void PrintVectInt(std::ostream &os, Face const& eList)
{
  int len, i;
  len=eList.size();
  for (i=0; i<len; i++)
    if (eList[i] == 1)
      os << " " << i;
  os << "\n";
}




struct GraphBitset {
public:
  GraphBitset() = delete;
  GraphBitset(int const& inpNbVert) : nbVert(inpNbVert)
  {
    LLAdj=boost::dynamic_bitset<>(nbVert*nbVert);
    HasVertexColor=false;
  }
  ~GraphBitset()
  {
    ~LLAdj;
    if (HasVertexColor)
      ListVertexColor.clear();
  }
  GraphBitset(GraphBitset const& eG)
  {
    nbVert=eG.GetNbVert();
    LLAdj=eG.GetBitset();
    HasVertexColor=eG.GetHasVertexColor();
    if (HasVertexColor)
      ListVertexColor=eG.GetListVertexColor();
  }
  GraphBitset operator=(GraphBitset const& eG)
  {
    nbVert=eG.GetNbVert();
    LLAdj=eG.GetBitset();
    HasVertexColor=eG.GetHasVertexColor();
    if (HasVertexColor)
      ListVertexColor=eG.GetListVertexColor();
    return *this;
  }
  // lighter stuff
  int GetNbVert() const
  {
    return nbVert;
  }
  boost::dynamic_bitset<> GetBitset() const
  {
    return LLAdj;
  }
  bool GetHasVertexColor() const
  {
    return HasVertexColor;
  }
  std::vector<int> GetListVertexColor() const
  {
    return ListVertexColor;
  }
  //
  void SetHasColor(bool const& TheVal)
  {
    if (TheVal == HasVertexColor) {
      return;
    }
    HasVertexColor=TheVal;
    if (TheVal) 
      ListVertexColor=std::vector<int>(nbVert);
    if (!TheVal)
      ListVertexColor.clear();
  }
  void SetColor(int const& iVert, int const& eColor)
  {
    ListVertexColor[iVert]=eColor;
  }
  std::vector<int> Adjacency(int const& iVert) const
  {
    std::vector<int> retList;
    for (int jVert=0; jVert<nbVert; jVert++) {
      int idxMat=iVert + nbVert*jVert;
      if (LLAdj[idxMat])
	retList.push_back(jVert);
    }
    return retList;
  }
  void AddAdjacent(int const& iVert, int const& jVert)
  {
    int idxMat=iVert + nbVert*jVert;
    LLAdj.set(idxMat,1);
  }
  void RemoveAdjacent(int const& iVert, int const& jVert)
  {
    int idxMat=iVert + nbVert*jVert;
    LLAdj.set(idxMat,1);
  }
  bool IsAdjacent(int const& iVert, int const& jVert) const
  {
    int idxMat=iVert + nbVert*jVert;
    if (LLAdj[idxMat])
      return true;
    return false;
  }
  int GetColor(int const& iVert) const
  {
    if (!HasVertexColor) {
      std::cerr << "Call to GetColor while HasVertexColor=false\n";
      throw TerminalException{1};
    }
    return ListVertexColor[iVert];
  }
private:
  int nbVert;
  boost::dynamic_bitset<> LLAdj;
  bool HasVertexColor=false;
  std::vector<int> ListVertexColor;
};

template<>
struct is_graphbitset_class<GraphBitset> {
  static const bool value = true;
};


template<typename Tret, typename Tgr>
inline typename std::enable_if<is_graphbitset_class<Tret>::value,Tret>::type InducedSubgraph(Tgr const& GR, std::vector<int> const& eList)
{
  int nbVert=GR.GetNbVert();
  std::vector<int> ListStat(nbVert,-1);
  int nbVertRed=eList.size();
  GraphBitset GRred(nbVertRed);
  for (int iVertRed=0; iVertRed<nbVertRed; iVertRed++) {
    int iVert=eList[iVertRed];
    ListStat[iVert]=iVertRed;
  }
  for (int iVertRed=0; iVertRed<nbVertRed; iVertRed++) {
    int iVert=eList[iVertRed];
    std::vector<int> LLadj=GR.Adjacency(iVert);
    for (int & jVert : LLadj) {
      int jVertRed=ListStat[jVert];
      if (jVertRed != -1) {
	GRred.AddAdjacent(iVertRed, jVertRed);
	GRred.AddAdjacent(jVertRed, iVertRed);
      }
    }
  }
  return GRred;
}



Face FullFace(int const& len)
{
  Face eFace(len);
  for (int u=0; u<len; u++)
    eFace[u]=1;
  return eFace;
}



ulong FaceToUnsignedLong(Face const& f)
{
  int len=f.size();
  if (len > 32) {
    std::cerr << "Too large value, conversion impossible";
    throw TerminalException{1};
  }
  ulong pos=0;
  ulong pow=1;
  for (int i=0; i<len; i++) {
    pos += pow*f[i];
    pow *= 2;
  }
  return pos;
}

Face UnsignedLongToFace(int const& len, ulong const& eVal)
{
  if (len > 32) {
    std::cerr << "length error\n";
    throw TerminalException{1};
  }
  ulong eWork = eVal;
  Face eFace(len);
  ulong pow=1;
  for (int i=0; i<len; i++) {
    ulong Pow2=pow * 2;
    ulong res=eWork % Pow2;
    if (res == pow) {
      eFace[i]=1;
      eWork -= pow;
    }
    pow=Pow2;
  }
  return eFace;
}

void VectVectInt_Magma_Print(std::ostream &os, std::vector<Face> const&ListOrbit)
{
  int nbOrbit=ListOrbit.size();
  os << "[";
  for (int iOrbit=0; iOrbit<nbOrbit; iOrbit++) {
    if (iOrbit>0)
      os << ",\n";
    Face eRepr=ListOrbit[iOrbit];
    int siz=eRepr.count();
    os << "[";
    int eVal=eRepr.find_first();
    for (int i=0; i<siz; i++) {
      if (i>0)
	os << ",";
      os << eVal;
      eVal=eRepr.find_next(eVal);
    }
    os << "]";
  }
  os << "]\n";
}



#endif

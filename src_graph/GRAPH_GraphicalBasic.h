#ifndef GRAPHICAL_BASIC_INCLUDE
#define GRAPHICAL_BASIC_INCLUDE

#include "Temp_common.h"





struct GraphListAdj {
public:
  GraphListAdj() = delete;
  GraphListAdj(int const& inpNbVert) : nbVert(inpNbVert)
  {
    ListListAdj.clear();
    std::vector<int> eList;
    for (int iVert=0; iVert<nbVert; iVert++) {
      ListListAdj.push_back(eList);
    }
    HasVertexColor=false;
  }
  ~GraphListAdj()
  {
  }
  GraphListAdj(GraphListAdj const& eG)
  {
    nbVert=eG.GetNbVert();
    ListListAdj=eG.GetListListAdj();
    HasVertexColor=eG.GetHasVertexColor();
    ListVertexColor=eG.GetListVertexColor();
  }
  GraphListAdj operator=(GraphListAdj const& eG)
  {
    nbVert=eG.GetNbVert();
    ListListAdj=eG.GetListListAdj();
    HasVertexColor=eG.GetHasVertexColor();
    ListVertexColor=eG.GetListVertexColor();
    return *this;
  }
  // lighter stuff
  int GetNbVert() const
  {
    return nbVert;
  }
  std::vector<std::vector<int> > GetListListAdj() const
  {
    return ListListAdj;
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
    if (TheVal == HasVertexColor)
      return;
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
    return ListListAdj[iVert];
  }
  void AddAdjacent(int const& iVert, int const& jVert)
  {
    ListListAdj[iVert].push_back(jVert);
  }
  void RemoveAdjacent(int const& iVert, int const& jVert)
  {
    // Not sure to find the right simple code for removing entry in std::vector
    std::vector<int> NewList;
    for (auto & eVert : ListListAdj[iVert])
      if (eVert != jVert)
	NewList.push_back(eVert);
    ListListAdj[iVert]=NewList;
  }
  bool IsAdjacent(int const& iVert, int const& jVert) const
  {
    for (auto & eVert : ListListAdj[iVert])
      if (eVert == jVert)
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
  std::vector<std::vector<int> > ListListAdj;
  bool HasVertexColor;
  std::vector<int> ListVertexColor;
};






struct GraphSparseImmutable {
public:
  GraphSparseImmutable() = delete;
GraphSparseImmutable(int const& _nbVert, std::vector<int> const& _ListStart, std::vector<int> const& _ListListAdj) : nbVert(_nbVert), ListStart(_ListStart), ListListAdj(_ListListAdj)
  {
    HasVertexColor=false;
  }
  ~GraphSparseImmutable()
  {
  }
  GraphSparseImmutable(MyMatrix<int> const& ListEdge, int const& _nbVert) : nbVert(_nbVert)
  {
    int nbEdge=ListEdge.rows();
    std::vector<int> ListDeg(nbVert,0);
    for (int iEdge=0; iEdge<nbEdge; iEdge++)
      for (int i=0; i<2; i++) {
	int eVert=ListEdge(iEdge,i);
	ListDeg[eVert]++;
      }
    ListStart.resize(nbVert+1,0);
    int nbAdj=0;
    for (int iVert=0; iVert<nbVert; iVert++) {
      int eDeg=ListDeg[iVert];
      ListStart[iVert+1] = ListStart[iVert] + eDeg;
      nbAdj += eDeg;
    }
    std::vector<int> ListPos(nbVert,0);
    ListListAdj.resize(nbAdj,0);
    for (int iEdge=0; iEdge<nbEdge; iEdge++)
      for (int i=0; i<2; i++) {
	int j=1-i;
	int eVert=ListEdge(iEdge,i);
	int fVert=ListEdge(iEdge,j);
	int eStart=ListStart[eVert];
	int pos=ListPos[eVert];
	ListListAdj[eStart + pos]=fVert;
	ListPos[eVert]++;
      }
  }
  GraphSparseImmutable(GraphSparseImmutable const& eG)
  {
    nbVert=eG.GetNbVert();
    ListStart=eG.GetListStart();
    ListListAdj=eG.GetListListAdj();
    HasVertexColor=eG.GetHasVertexColor();
    ListVertexColor=eG.GetListVertexColor();
  }
  GraphSparseImmutable operator=(GraphSparseImmutable const& eG)
  {
    nbVert=eG.GetNbVert();
    ListStart=eG.GetListStart();
    ListListAdj=eG.GetListListAdj();
    HasVertexColor=eG.GetHasVertexColor();
    ListVertexColor=eG.GetListVertexColor();
    return *this;
  }
  // lighter stuff
  int GetNbVert() const
  {
    return nbVert;
  }
  std::vector<int> GetListListAdj() const
  {
    return ListListAdj;
  }
  std::vector<int> GetListStart() const
  {
    return ListStart;
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
    int eStart=ListStart[iVert];
    int eEnd=ListStart[iVert+1];
    std::vector<int> TheRet;
    for (int i=eStart; i<eEnd; i++)
      TheRet.push_back(ListListAdj[i]);
    return TheRet;
  }
  bool IsAdjacent(int const& iVert, int const& jVert) const
  {
    int eStart=ListStart[iVert];
    int eEnd=ListStart[iVert+1];
    for (int i=eStart; i<eEnd; i++)
      if (ListListAdj[i] == jVert)
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
  std::vector<int> ListStart;
  std::vector<int> ListListAdj;
  bool HasVertexColor;
  std::vector<int> ListVertexColor;
};


template<>
struct is_graphsparseimmutable_class<GraphSparseImmutable> {
  static const bool value = true;
};





template<typename Tgr>
std::vector<int> ConnectedComponents_vector(Tgr const& GR)
{
  int nbVert=GR.GetNbVert();
  std::vector<int> ListStatus(nbVert,-1);
  int iStatus=0;
  auto Assignment=[&](int const & iVert) -> void {
    std::vector<int> TheSet{iVert};
    ListStatus[iVert]=iStatus;
    while(true) {
      std::vector<int> NewSet;
      for (int & eVal : TheSet) {
        for (int & fVal : GR.Adjacency(eVal)) {
          if (ListStatus[fVal] == -1) {
            ListStatus[fVal]=iStatus;
            NewSet.push_back(fVal);
          }
        }
      }
      if (NewSet.size() == 0)
        break;
      TheSet=NewSet;
    }
    iStatus++;
  };
  for (int iVert=0; iVert<nbVert; iVert++)
    if (ListStatus[iVert] == -1)
      Assignment(iVert);
  return ListStatus;
}



template<typename Tgr>
std::vector<std::vector<int>> ConnectedComponents_set(Tgr const& GR)
{
  int nbVert=GR.GetNbVert();
  std::vector<int> ListStatus=ConnectedComponents_vector(GR);
  int nbConn=VectorMax(ListStatus)+1;
  std::vector<std::vector<int> > ListConn(nbConn);
  for (int iVert=0; iVert<nbVert; iVert++) {
    int iStatus=ListStatus[iVert];
    ListConn[iStatus].push_back(iVert);
  }
  return ListConn;
}



template<typename Tret, typename Tgr>
inline typename std::enable_if<is_graphsparseimmutable_class<Tret>::value,Tret>::type InducedSubgraph(Tgr const& GR, std::vector<int> const& eList)
{
  int nbVert=GR.GetNbVert();
  std::vector<int> ListStat(nbVert,-1);
  int nbVertRed=eList.size();
  for (int iVertRed=0; iVertRed<nbVertRed; iVertRed++) {
    int iVert=eList[iVertRed];
    ListStat[iVert]=iVertRed;
  }
  std::vector<MyVector<int>> PreListEdge;
  for (int iVertRed=0; iVertRed<nbVertRed; iVertRed++) {
    int iVert=eList[iVertRed];
    std::vector<int> LLadj=GR.Adjacency(iVert);
    for (int & jVert : LLadj) {
      int jVertRed=ListStat[jVert];
      if (jVertRed != -1 && jVertRed > iVertRed) {
	MyVector<int> eVect(2);
	eVect(0)=iVertRed;
	eVect(1)=jVertRed;
	PreListEdge.push_back(eVect);
      }
    }
  }
  MyMatrix<int> ListEdge=MatrixFromVectorFamily(PreListEdge);
  return GraphSparseImmutable(ListEdge, nbVertRed);
}




#endif

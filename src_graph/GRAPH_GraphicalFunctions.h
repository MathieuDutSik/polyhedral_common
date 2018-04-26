#ifndef TEMP_GRAPHICAL_FUNCTION
#define TEMP_GRAPHICAL_FUNCTION

#include "COMB_Combinatorics.h"
#include "GRAPH_GraphicalBasic.h"
#include "MAT_Matrix.h"


struct GraphFunctional {
public:
  GraphFunctional() = delete;
  GraphFunctional(int const& inpNbVert, std::function<bool(int const&,int const&)> const& eFct) : nbVert(inpNbVert), f(eFct)
  {
    HasVertexColor=false;
  }
  ~GraphFunctional()
  {
  }
  GraphFunctional(GraphFunctional const& eG)
  {
    nbVert=eG.GetNbVert();
    f=eG.GetFCT();
    HasVertexColor=eG.GetHasVertexColor();
    fColor=eG.GetFColor();
  }
  GraphFunctional operator=(GraphFunctional const& eG)
  {
    nbVert=eG.GetNbVert();
    f=eG.GetFCT();
    HasVertexColor=eG.GetHasVertexColor();
    fColor=eG.GetFColor();
    return *this;
  }
  // lighter stuff
  int GetNbVert() const
  {
    return nbVert;
  }
  std::function<bool(int const&,int const&)> GetFCT() const
  {
    return f;
  }
  bool GetHasVertexColor() const
  {
    return HasVertexColor;
  }
  std::function<int(int const&)> GetFColor() const
  {
    return fColor;
  }
  //
  void SetFColor(std::function<int(int const&)> const& inpFColor)
  {
    HasVertexColor=true;
    fColor=inpFColor;
  }
  std::vector<int> Adjacency(int const& iVert) const
  {
    std::vector<int> retList;
    for (int jVert=0; jVert<nbVert; jVert++)
      if (f(iVert, jVert))
	retList.push_back(jVert);
    return retList;
  }
  bool IsAdjacent(int const& iVert, int const& jVert) const
  {
    return f(iVert, jVert);
  }
  int GetColor(int const& iVert) const
  {
    return fColor(iVert);
  }
private:
  int nbVert;
  std::function<bool(int const&,int const&)> f;
  bool HasVertexColor;
  std::function<int(int const&)> fColor;
};


template <typename T>
struct is_functional_graph_class {
  static const bool value = false;
};

template <>
struct is_functional_graph_class<GraphFunctional> {
  static const bool value = true;
};









template<typename Tgr>
void GRAPH_PrintOutput(std::ostream &os, Tgr const& GR)
{
  int nbVert=GR.GetNbVert();
  os << "nbVert=" << nbVert << "\n";
  for (int iVert=0; iVert<nbVert; iVert++) {
    std::vector<int> LVert=GR.Adjacency(iVert);
    os << LVert.size();
    for (auto &eVert : LVert)
      os << " " << eVert;
    os << "\n";
  }
}

template<typename Tgr>
void GRAPH_PrintOutputGAP_vertex_colored(std::string const& eFile, Tgr const& GR)
{
  int nbVert=GR.GetNbVert();
  std::vector<std::vector<int> > ListEdges;
  std::ofstream os(eFile);
  os << "local ListAdjacency, ThePartition;\n";
  os << "ListAdjacency:=[";
  for (int iVert=0; iVert<nbVert; iVert++) {
    if (iVert > 0)
      os << ",\n";
    os << "[";
    std::vector<int> LVert=GR.Adjacency(iVert);
    WriteStdVectorGAP(os, LVert);
  }
  os << "];\n";
  std::vector<int> ListColor=GR.GetListVertexColor();
  int nbBlock=VectorMax(ListColor)+1;
  std::vector<std::vector<int>> ListBlock(nbBlock);
  for (int iVert=0; iVert<nbVert; iVert++) {
    int eColor=ListColor[iVert];
    ListBlock[eColor].push_back(iVert+1);
  }
  os << "ThePartition:=[";
  for (int iBlock=0; iBlock<nbBlock; iBlock++) {
    if (iBlock>0)
      os << ",";
    WriteStdVectorGAP(os, ListBlock[iBlock]);
  }
  os << "];\n";
  os << "return [ListAdjacency, ThePartition];\n";
}



template<typename Tgr>
void GRAPH_PrintOutputGAP(std::ostream &os, Tgr const& GR)
{
  int nbVert=GR.GetNbVert();
  std::vector<std::vector<int> > ListEdges;
  for (int iVert=0; iVert<nbVert; iVert++) {
    std::vector<int> LVert=GR.Adjacency(iVert);
    for (auto & eVert : LVert)
      if (eVert > iVert)
	ListEdges.push_back({iVert+1, eVert+1});
  }
  int nbEdge=ListEdges.size();
  os << "local ListEdges, eEdge, GRA;\n";
  os << "ListEdges:=[";
  for (int iEdge=0; iEdge<nbEdge; iEdge++) {
    if (iEdge > 0)
      os << ",";
    os << "[" << ListEdges[iEdge][0] << "," << ListEdges[iEdge][1] << "]";
  }
  os << "];\n";
  os << "GRA:=NullGraph(Group(()), " << nbVert << ");\n";
  os << "for eEdge in ListEdges\n";
  os << "do\n";
  os << "  AddEdgeOrbit(GRA, eEdge);\n";
  os << "  AddEdgeOrbit(GRA, Reversed(eEdge));\n";
  os << "od;\n";
  os << "return GRA;\n";
}



template<typename Tgr>
Tgr GRAPH_Read(std::istream & is)
{
  if (!is.good()) {
    std::cerr << "GRAPH_Read operation failed because stream is not valied\n";
    throw TerminalException{1};
  }
  int nbVert;
  is >> nbVert;
  Tgr eGR(nbVert);
  for (int iVert=0; iVert<nbVert; iVert++) {
    int eDeg;
    is >> eDeg;
    for (int i=0; i<eDeg; i++) {
      int eAdj;
      is >> eAdj;
      eGR.AddAdjacent(iVert, eAdj);
    }
  }
  return eGR;
}






template<typename Tgr>
MyMatrix<int> ShortestPathDistanceMatrix(Tgr const& GR)
{
  int nbVert=GR.GetNbVert();
  MyMatrix<int> StatusMat(nbVert,nbVert);
  MyMatrix<int> DistMat(nbVert,nbVert);
  for (int iVert=0; iVert<nbVert; iVert++)
    for (int jVert=0; jVert<nbVert; jVert++)
      StatusMat(iVert,jVert)=0;
  for (int iVert=0; iVert<nbVert; iVert++) {
    StatusMat(iVert,iVert)=2;
    DistMat(iVert,iVert)=0;
  }
  int UpperBoundDiam=nbVert+2;
  int iter=0;
  while(true) {
    bool IsFinished=true;
    iter++;
    for (int iVert=0; iVert<nbVert; iVert++)
      for (int jVert=0; jVert<nbVert; jVert++)
	if (StatusMat(iVert,jVert) == 2) {
	  IsFinished=false;
	  StatusMat(iVert,jVert)=1;
	  std::vector<int> LLAdj=GR.Adjacency(jVert);
	  for (int & eAdj : LLAdj)
	    if (StatusMat(iVert,eAdj) == 0) {
	      StatusMat(iVert,eAdj)=-1;
	      DistMat(iVert,eAdj)=DistMat(iVert,jVert) + 1;
	    }
	}
    for (int iVert=0; iVert<nbVert; iVert++)
      for (int jVert=0; jVert<nbVert; jVert++)
	if (StatusMat(iVert,jVert) == -1)
	  StatusMat(iVert,jVert)=2;
    if (iter > UpperBoundDiam) {
      for (int iVert=0; iVert<nbVert; iVert++)
	DistMat(iVert,iVert)=-1;
      break;
    }
    if (IsFinished)
      break;
  }
  return DistMat;
}


template<typename Tgr>
std::vector<std::vector<int>> GRAPH_FindAllShortestPath(Tgr const& GR, int const& x, int const& y)
{
  MyMatrix<int> DistMat=ShortestPathDistanceMatrix(GR);
  int eDist=DistMat(x,y);
  if (eDist == -1)
    return {};
  std::vector<std::vector<int>> ListPath{{x}};
  for (int i=1; i<=eDist; i++) {
    std::vector<std::vector<int>> NewListPath;
    for (auto & ePath : ListPath) {
      int CurrentVert=ePath[ePath.size() - 1];
      std::vector<int> LAdj = GR.Adjacency(CurrentVert);
      for (auto & u : LAdj) {
	if (DistMat(u, y) == eDist - i) {
	  std::vector<int> fPath = ePath;
	  fPath.push_back(u);
	  NewListPath.push_back(fPath);
	}
      }
    }
    ListPath = NewListPath;
  }
  return ListPath;
}



template<typename Tgr>
int Diameter(Tgr const& GR)
{
  MyMatrix<int> DistMat = ShortestPathDistanceMatrix(GR);
  if (DistMat(0,0) == -1)
    return 1;
  int TheDiam=0;
  int nbVert=GR.GetNbVert();
  for (int iVert=0; iVert<nbVert; iVert++)
    for (int jVert=0; jVert<nbVert; jVert++) {
      int eDist=DistMat(iVert, jVert);
      if (eDist > TheDiam)
	TheDiam=eDist;
    }
  return TheDiam;
}

template<typename Tgr>
bool IsSimpleGraph(Tgr const& GR)
{
  int nbVert=GR.GetNbVert();
  for (int iVert=0; iVert<nbVert; iVert++)
    if (GR.IsAdjacent(iVert, iVert))
      return false;
  return true;
}

template<typename Tgr>
std::vector<std::vector<int>> GetEdgeSet(Tgr const& GR)
{
  int nbVert=GR.GetNbVert();
  /*
  int nbEdge=0;
  for (int iVert=0; iVert<nbVert-1; iVert++)
    for (int jVert=iVert+1; jVert<nbVert; jVert++)
      if (GR.IsAdjacent(iVert,jVert) == 1)
      nbEdge++;
  MyMatrix<int> Edges(nbEdge,2);*/
  std::vector<std::vector<int> > Edges;
  for (int iVert=0; iVert<nbVert-1; iVert++)
    for (int jVert=iVert+1; jVert<nbVert; jVert++)
      if (GR.IsAdjacent(iVert,jVert))
	Edges.push_back({iVert, jVert});
  return Edges;
}





template<typename Tgr>
std::vector<std::vector<int>> GRAPH_FindAllCycles(Tgr const& GR)
{
  std::vector<std::vector<int>> Edges = GetEdgeSet(GR);
  int nbEdge=Edges.size();
  int nbVert=GR.GetNbVert();
  std::set<std::vector<int>> ListCycle;
  auto FuncInsertCycle=[&](std::vector<int> const& eL) -> void {
    std::vector<int> eLcan = MinimumDihedralOrbit(eL);
    ListCycle.insert(eLcan);
  };
  for (int iEdge=0; iEdge<nbEdge; iEdge++) {
    GraphBitset GRred(nbVert);
    //    std::cerr << "nbVert=" << nbVert << "\n";
    for (int jEdge=0; jEdge<nbEdge; jEdge++)
      if (iEdge != jEdge) {
	int eVert=Edges[jEdge][0];
	int fVert=Edges[jEdge][1];
	GRred.AddAdjacent(eVert, fVert);
	GRred.AddAdjacent(fVert, eVert);
	//	std::cerr << "eVert=" << eVert << " fVert=" << fVert << "\n";
      }
    int x=Edges[iEdge][0];
    int y=Edges[iEdge][1];
    std::cerr << "iEdge=" << iEdge << " / " << nbEdge << " x=" << x << " y=" << y << "\n";
    std::vector<std::vector<int>> ListPath = GRAPH_FindAllShortestPath(GRred, x, y);
    std::cerr << "|ListPath|=" << ListPath.size() << "\n";
    for (auto & eCycle : ListPath)
      FuncInsertCycle(eCycle);
  }
  std::vector<std::vector<int>> ListCycleR;
  for (auto & eCycle : ListCycle)
    ListCycleR.push_back(eCycle);
  return ListCycleR;
}











template<typename Tgr>
bool IsClique(Tgr const& GR, std::vector<int> const& eList)
{
  int len=eList.size();
  for (int i=0; i<len-1; i++)
    for (int j=i+1; j<len; j++)
      if (!GR.IsAdjacent(eList[i],eList[j]))
	return false;
  return true;
}




template<typename Tgr>
std::vector<int> StartingCell(Tgr const& GR)
{
  int r;
  std::vector<int> Adj;
  //
  std::vector<int> A=GR.Adjacency(0);
  int eC=A[0];
  //  std::cerr << "eC=" << eC << "\n";
  Adj=IntersectionVect(GR.Adjacency(eC), A);
  r=Adj.size();
  std::vector<int> TT;
  //  std::cerr << "r=" << r << "\n";
  //  std::cerr << "1 : Adj=";
  //  WriteVectorInt_GAP(std::cerr, Adj);
  //  std::cerr << "\n";
  if (r == 0) {
    return {0, eC};
  }
  else {
    if (r == 1) {
      int Elt=Adj[0];
      int h=IntersectionVect(GR.Adjacency(Elt), GR.Adjacency(0)).size();
      int k=IntersectionVect(GR.Adjacency(Elt), GR.Adjacency(eC)).size();
      if (h == 1 && k == 1) {
	return {0, eC, Elt};
      }
      else {
	if (h>1) {
	  TT={Elt,0};
	}
	else {
	  TT={Elt,eC};
	}
      }
    }
    else {
      TT={0,eC};
      //      std::cerr << "TT=";
      //      WriteVectorInt_GAP(std::cerr, TT);
      //      std::cerr << "\n";
    }
  }
  
  Adj=IntersectionVect(GR.Adjacency(TT[0]), GR.Adjacency(TT[1]));
  //  std::cerr << "2 : Adj=";
  //  WriteVectorInt_GAP(std::cerr, Adj);
  //  std::cerr << "\n";
  r=Adj.size();
  std::vector<int> Lodd;
  int nbVert=GR.GetNbVert();
  for (int & i : Adj) {
    //    std::cerr << "i=" << i << "\n";
    std::vector<int> VertSet(nbVert);
    for (int iVert=0; iVert<nbVert; iVert++)
      VertSet[iVert]=iVert;
    std::vector<int> hSet{i,TT[0],TT[1]};
    //    std::cerr << "  hSet=";
    //    WriteVectorInt_GAP(std::cerr, hSet);
    //    std::cerr << "\n";
    std::vector<int> SE=DifferenceVect(VertSet, hSet);
    //    std::cerr << "  SE=";
    //    WriteVectorInt_GAP(std::cerr, SE);
    //    std::cerr << "\n";
    bool test=true;
    for (int & j : SE) {
      int oddness=IntersectionVect<int>(GR.Adjacency(j), hSet).size();
      //      std::cerr << " j=" << j << " oddness=" << oddness << "\n";
      if (oddness == 1 || oddness == 3)
	test=false;
    }
    if (!test)
      Lodd.push_back(i);
  }
  int s=Lodd.size();
  //  std::cerr << "r=" << r << " s=" << s << "\n";
  if (r == 2 && s == 0) {
    return {Adj[0], TT[0], TT[1]};
  }
  else {
    if (s == r || s == r-1) {
      Lodd.push_back(TT[0]);
      Lodd.push_back(TT[1]);
      if (!IsClique(GR, Lodd))
	return {-1};
      return Lodd;
    }
    else {
      return {-1};
    }
  }
}

template<typename Tgr>
std::vector<std::vector<int> > SpanningTree(Tgr const& GR)
{
  int nbVert=GR.GetNbVert();
  std::vector<int> ListStatus(nbVert,0);
  ListStatus[0]=1;
  std::vector<int> ListActiveVert{0};
  std::vector<std::vector<int> > TheSpann;
  while(true) {
    std::vector<int> NewListActiveVert;
    for (auto & eVert : ListActiveVert) {
      //      std::cerr << "Treating eVert=" << eVert << "\n";
      std::vector<int> LLadj=GR.Adjacency(eVert);
      for (auto & fVert : LLadj) {
	if (ListStatus[fVert] == 0) {
	  ListStatus[fVert]=1;
	  TheSpann.push_back(VectorAsSet<int>({eVert, fVert}));
	  NewListActiveVert.push_back(fVert);
	}
      }
    }
    if (NewListActiveVert.size() == 0)
      break;
    ListActiveVert=NewListActiveVert;
  }
  return TheSpann;
}



template<typename Tgr>
std::vector<std::vector<int> > InverseLineGraphConnected(Tgr const& GR)
{
  int nbVert=GR.GetNbVert();
  if (nbVert == 1) {
    return {{0,1}};
  }
  std::vector<int> eCell=StartingCell(GR);
  std::cerr << "eCell=";
  WriteVectorInt_GAP(std::cerr, eCell);
  std::cerr << "\n";
  std::vector<int> eCellS=VectorAsSet(eCell);
  if (eCell[0] == -1) {
    //    std::cerr << "Leaving InverseLineGraphConnected 1\n";
    return {{-1}};
  }
  std::vector<std::vector<int> > P={eCell};
  std::vector<std::vector<int> > TotalEdge(nbVert);
  for (int iVert=0; iVert<nbVert; iVert++)
    TotalEdge[iVert]=GR.Adjacency(iVert);
  for (auto & eVert : eCell)
    TotalEdge[eVert]=DifferenceVect(TotalEdge[eVert], eCell);
  int PreTot=0;
  for (int iVert=0; iVert<nbVert; iVert++) {
    int eSize=TotalEdge[iVert].size();
    //    std::cerr << "iVert=" << iVert << " |Adj|=" << eSize << "  Adj=";
    //    WriteVectorInt_GAP(std::cerr, TotalEdge[iVert]);
    //    std::cerr << "\n";
    PreTot += eSize;
  }
  int Tot=PreTot/2;
  std::cerr << "Tot=" << Tot << "\n";
  auto GetAdding=[&]() -> int {
    for (int & iVert : eCellS)
      if (TotalEdge[iVert].size() > 0)
	return iVert;
    std::cerr << "Error in GetAdding\n";
    throw TerminalException{1};
  };
  while(true) {
    std::cerr << "Tot=" << Tot << "\n";
    if (Tot == 0)
      break;
    int Adding=GetAdding();
    std::cerr << "Adding=" << Adding << "\n";
    eCell=UnionVect(TotalEdge[Adding], {Adding});
    //    std::cerr << "eCell=";
    //    WriteVectorInt_GAP(std::cerr, eCell);
    //    std::cerr << "\n";
    for (auto & eVert : eCell)
      for (auto & fVert : eCell)
	if (eVert != fVert)
	  if (PositionVect(TotalEdge[eVert], fVert) == -1)
	    return {{-1}};
    //    bool test=IsClique(GR, eCell);
    //std::cerr << "test=" << test << "\n";
    //    if (!test) {
    //      std::cerr << "Leaving InverseLineGraphConnected 2\n";
    //      return {{-1}};
    //    }
    P.push_back(eCell);
    eCellS=VectorAsSet(UnionVect(eCellS, eCell));
    //    std::cerr << "eCellS=";
    //    WriteVectorInt_GAP(std::cerr, eCellS);
    //    std::cerr << "\n";
    for (int & iV : eCell)
      TotalEdge[iV]=DifferenceVect(TotalEdge[iV], eCell);
    int Csiz=eCell.size();
    int nbEdge=Csiz*(Csiz - 1)/2;
    std::cerr << "Before Csiz=" << Csiz << " nbEdge=" << nbEdge << " Tot=" << Tot << "\n";
    Tot -= nbEdge;
    std::cerr << "After Tot=" << Tot << "\n";
  }
  for (int iVert=0; iVert<nbVert; iVert++) {
    int nb=0;
    for (auto & eP : P) {
      int pos=PositionVect(eP, iVert);
      if (pos != -1)
	nb++;
    }
    if (nb == 1)
      P.push_back({iVert});
  }
  std::vector<std::vector<int> > Label(nbVert);
  int Plen=P.size();
  for (int iVert=0; iVert<nbVert; iVert++) {
    for (int iP=0; iP<Plen; iP++) {
      int pos=PositionVect(P[iP], iVert);
      if (pos != -1)
	Label[iVert].push_back(iP);
    }
  }
  return Label;
}


template<typename T>
void Print_VectorVector(std::ostream &os, std::vector<std::vector<T> > const& TheList)
{
  int siz=TheList.size();
  os << "|TheList|=" << siz << "\n";
  for (int i=0; i<siz; i++) {
    os << "i=" << i << " V=";
    for (auto & eVal : TheList[i])
      os << " " << eVal;
    os << "\n";
  }
}




template<typename Tgr>
std::vector<std::vector<int> > InverseLineGraph(Tgr const& GR)
{
  int nbVert=GR.GetNbVert();
  std::vector<std::vector<int> > ListConn=ConnectedComponents_set(GR);
  int nbConn=ListConn.size();
  auto Hsize=[&](std::vector<std::vector<int> > const& ListListSet) -> int {
    int eMax=0;
    for (auto & eListSet : ListListSet)
      for (auto & eVal : eListSet)
	if (eVal > eMax)
	  eMax=eVal;
    return eMax+1;
  };
  std::vector<std::vector<int> > ListLabel(nbVert);
  int TotShift=0;
  for (int iConn=0; iConn<nbConn; iConn++) {
    std::vector<int> eConn=ListConn[iConn];
    int sizConn=eConn.size();
    //    std::cerr << "iConn=" << iConn << " sizConn=" << sizConn << "\n";
    GraphBitset GRind=InducedSubgraph<GraphBitset,Tgr>(GR, eConn);
    //    std::cerr << "-----------------------------\n";
    //    std::cerr << "GRind=";
    //    GRAPH_PrintOutput(std::cerr, GRind);
    //    GRAPH_PrintOutputGAP(std::cerr, GRind);
    std::vector<std::vector<int> > PartLabel=InverseLineGraphConnected(GRind);
    //    std::cerr << "GRind : PartLabel=\n";
    //    Print_VectorVector(std::cerr, PartLabel);
    int eFirstFirst=PartLabel[0][0];
    if (eFirstFirst == -1)
      return {{-1}};
    for (int iVertConn=0; iVertConn<sizConn; iVertConn++) {
      int eVert=eConn[iVertConn];
      std::vector<int> eLabel;
      for (int & eVal : PartLabel[iVertConn])
	eLabel.push_back(eVal + TotShift);
      ListLabel[eVert]=eLabel;
    }
    int eShift=Hsize(PartLabel);
    TotShift += eShift;
  }
  return ListLabel;
}

MyMatrix<int> CreateEmbedding(int const& StartPoint, std::vector<std::vector<int> > const& ListEdge, std::vector<std::vector<int> > const& ListLabel)
{
  int nbEdge=ListEdge.size();
  //
  int nbVert=0;
  for (auto & eEdge : ListEdge)
    for (auto & eVal : eEdge)
      if (eVal > nbVert)
	nbVert=eVal;
  nbVert++;
  //
  int nbLabel=0;
  for (auto & eLabel : ListLabel)
    for (auto & eVal : eLabel)
      if (eVal > nbLabel)
	nbLabel=eVal;
  nbLabel++;
  //
  MyMatrix<int> Embedding(nbVert, nbLabel);
  for (int i=0; i<nbLabel; i++)
    Embedding(StartPoint,i)=0;
  std::vector<int> ListStatus(nbVert,0);
  ListStatus[StartPoint]=1;
  while(true) {
    bool IsFinished=true;
    for (int iEdge=0; iEdge<nbEdge; iEdge++) {
      std::vector<int> eEdge=ListEdge[iEdge];
      int eStat1=ListStatus[eEdge[0]];
      int eStat2=ListStatus[eEdge[1]];
      int RelPos=-1;
      int eVert1=-1, eVert2=-1;
      if (eStat1 == 0 && eStat2 == 1) {
	RelPos=0;
	eVert1=eEdge[1];
	eVert2=eEdge[0];
      }
      if (eStat2 == 0 && eStat1 == 1) {
	RelPos=0;
	eVert1=eEdge[0];
	eVert2=eEdge[1];
      }
      if (RelPos != -1) {
	IsFinished=false;
	for (int i=0; i<nbLabel; i++)
	  Embedding(eVert2,i) = Embedding(eVert1,i);
	for (auto & eVal : ListLabel[iEdge]) {
	  Embedding(eVert2,eVal) = 1 - Embedding(eVert1,eVal);
	}
	ListStatus[eVert2]=1;
      }
    }
    if (IsFinished)
      break;
  }
  return Embedding;
}

template<typename Tgr>
std::vector<std::vector<int> > EnumerationClique(Tgr const& GR, int const& kSiz)
{
  int nbVert=GR.GetNbVert();
  std::vector<std::vector<int> > eList;
  if (kSiz == 1) {
    for (int iVert=0; iVert<nbVert; iVert++)
      eList.push_back({iVert});
  }
  else {
    for (auto & eLowList : EnumerationClique(GR, kSiz-1) ) {
      int eFirst=eLowList[kSiz-2]+1;
      for (int i=eFirst; i<nbVert; i++) {
	bool test=true;
	for (auto & eVal : eLowList)
	  if (!GR.IsAdjacent(eVal,i))
	    test=false;
	if (test) {
	  std::vector<int> NewList=UnionVect(eLowList, {i});
	  eList.push_back(NewList);
	}
      }
    }
  }
  return eList;
}


int L1_distance(std::vector<int> const& V1, std::vector<int> const& V2)
{
  size_t siz=V1.size();
  assert(siz == V2.size());
  int dist=0;
  for (size_t i=0; i<siz; i++)
    if (V1[i] != V2[i])
      dist++;
  return dist;
}


template<typename Tgr>
std::vector<MyMatrix<int> > GRAPH_S_Embedding(Tgr const& GR, int const& s, long const& MaxIter, long& iter)
{
  iter=0;
  if (!IsSimpleGraph(GR)) {
    std::cerr << "The graph should be simple\n";
    return {};
  }
  MyMatrix<int> M=ShortestPathDistanceMatrix(GR);
  std::cerr << "M=\n";
  WriteMatrix(std::cerr, M);
  if (M(0,0) == -1) {
    std::cerr << "The graph should be connected\n";
    return {};
  }
  std::vector<int> SetZero{0};
  std::vector<int> SetOne{1};
  std::vector<int> SetTwo{2};
  std::cerr << "****Start treating " << s << "-embedding\n";
  int n=GR.GetNbVert();
  std::vector<std::vector<int> > ListEdge=GetEdgeSet(GR);
  int nbEdge=ListEdge.size();
  std::cerr << "nbEdge=" << nbEdge << "\n";
  int DiamGraph=Diameter(GR);
  std::cerr << "DiamGraph=" << DiamGraph << "\n";
  auto EdgeDist=[&](int const& iEdge, int const& jEdge) -> int {
    int a, b, c, d, swp;
    a=ListEdge[iEdge][0];
    b=ListEdge[iEdge][1];
    c=ListEdge[jEdge][0];
    d=ListEdge[jEdge][1];
    if (M(a,c) > s || M(a,d) > s || M(b,c) > s || M(b,d) > s)
      return -400;
    std::vector<int> ListDist;
    for (int u=0; u<n; u++)
      if (T_abs(M(a,u) - M(b,u)) == 1 &&
	  T_abs(M(c,u) - M(d,u)) == 1 &&
	  M(a,u) <= s && M(b,u) <= s && 
	  M(c,u) <= s && M(d,u) <= s) {
	if (M(u,b) == M(u,a) -1) {
          swp=a;
          a=b;
          b=swp;
	}
	if (M(u,d) == M(u,c) -1) {
          swp=c;
          c=d;
          d=swp;
	}
	int DeltaE=1 + M(u,a) - M(u,b);
	int DeltaF=1 + M(u,c) - M(u,d);
	if (DeltaE != 0 || DeltaF != 0) {
	  std::cerr << "Not correct u distances\n";
	  throw TerminalException{1};
	}
	int eDist=-M(b,d) + M(a,d) + M(b,c) - M(a,c);
        ListDist.push_back(eDist);
      }
    std::vector<int> ListDistRed=VectorAsSet(ListDist);
    if (ListDistRed.size() == 0) {
      return -400;
    }
    if (ListDistRed.size() > 1) {
      std::cerr << "We find |ListDist|=" << ListDistRed.size() << "ListDist=";
      WriteVectorInt_GAP(std::cerr, ListDistRed);
      std::cerr << "\n";
    }
    for (auto & eVal : ListDistRed)
      if (eVal >= 0)
	return eVal;
    return -400;
  };
  std::vector<std::vector<int> > ConSet(nbEdge);
  for (int iEdge=0; iEdge<nbEdge; iEdge++)
    ConSet[iEdge]={iEdge};
  for (int iEdge=0; iEdge<nbEdge; iEdge++) {
    std::cerr << "iEdge=" << iEdge << " e=" << ListEdge[iEdge][0] << "," << ListEdge[iEdge][1] << "\n";
  }
  MyMatrix<std::vector<int> > MCE(nbEdge,nbEdge);
  int NbPrev=0;
  int NbUnsolved=0;
  std::vector<std::vector<int> > ListPairNegative;
  for (int iEdge=0; iEdge<nbEdge-1; iEdge++)
    for (int jEdge=iEdge+1; jEdge<nbEdge; jEdge++) {
      int val=EdgeDist(iEdge, jEdge);
      std::vector<int> LVal;
      std::cerr << "iEdge=" << iEdge << " jEdge=" << jEdge << " val=" << val << "\n";
      if (val == -400) {
	LVal={0,1,2};
	NbUnsolved++;
      }
      else {
	LVal={val};
	if (val > 2 || val < -2) {
	  std::cerr << "Logical error for the distances\n";
	  throw TerminalException{1};
	}
	if (val < 0) {
	  ListPairNegative.push_back({iEdge, jEdge});
	}
      }
      MCE(iEdge,jEdge)=LVal;
      MCE(jEdge,iEdge)=LVal;
      NbPrev+=2;
    }
  std::cerr << "NbUnsolved=" << NbUnsolved << "\n";
  std::cerr << "NbPrev    =" << NbPrev << "\n";
  std::cerr << "DiamGraph=" << DiamGraph << "\n";
  int nbNegative=ListPairNegative.size();
  if (nbNegative > 0) {
    std::cerr << "Found some pairs of edges with intersection of negative size\n";
    std::cerr << "|ListPairNegative|=" << nbNegative << "\n";
    for (int iNeg=0; iNeg<nbNegative; iNeg++) {
      int iEdge=ListPairNegative[iNeg][0];
      int jEdge=ListPairNegative[iNeg][1];
      int eDist=EdgeDist(iEdge,jEdge);
      int a=ListEdge[iEdge][0];
      int b=ListEdge[iEdge][1];
      int c=ListEdge[jEdge][0];
      int d=ListEdge[jEdge][1];
      std::cerr << "i=" << iNeg << " e={" << a << "," << b << "} f={" << c << "," << d << "} eDist=" << eDist << "\n";
    }
    throw TerminalException{1};
  }
  if (NbPrev == NbUnsolved) {
    std::cerr << "Nothing computable. Error actually\n";
    return {};
  }
  auto GetGraphIdentityEdge=[&ListEdge](std::vector<std::vector<int> > const& LEdge, MyMatrix<std::vector<int> > const& MCEwork) -> GraphBitset {
    int nbEdgeLoc=LEdge.size();
    GraphBitset Aspec(nbEdgeLoc);
    for (int iEdge=0; iEdge<nbEdgeLoc-1; iEdge++)
      for (int jEdge=iEdge+1; jEdge<nbEdgeLoc; jEdge++) {
	//	std::cerr << "iEdge=" << iEdge << " jEdge=" << jEdge << "\n";
	int iPos=PositionVect(ListEdge, LEdge[iEdge]);
	int jPos=PositionVect(ListEdge, LEdge[jEdge]);
	//	std::cerr << "iPos=" << iPos << " jPos=" << jPos << "\n";
	//	std::cerr << "nbRow=" << MCE.rows() << " nbCol=" << MCE.cols() << "\n";
	if (MCEwork(iPos,jPos) == std::vector<int>({2})) {
	  Aspec.AddAdjacent(iEdge, jEdge);
	  Aspec.AddAdjacent(jEdge, iEdge);
	}
      }
    return Aspec;
  };
  auto IsCoherentEdgeValues=[&ListEdge](std::vector<std::vector<int> > const& Cspec, std::vector<std::vector<int> > const& LEdge, MyMatrix<std::vector<int> > const & MCEwork) ->bool {
    int nbComp=Cspec.size();
    for (int iComp=0; iComp<nbComp-1; iComp++)
      for (int jComp=iComp+1; jComp<nbComp; jComp++) {
	std::vector<int> eComp=Cspec[iComp];
	std::vector<int> fComp=Cspec[jComp];
	int siz1=eComp.size();
	int siz2=fComp.size();
	bool IsFirst=true;
	std::vector<int> CommonVal;
	for (int iElt=0; iElt<siz1; iElt++)
	  for (int jElt=0; jElt<siz2; jElt++) {
	    std::vector<int> eEdge=LEdge[eComp[iElt]];
	    std::vector<int> fEdge=LEdge[fComp[jElt]];
	    int iPos=PositionVect(ListEdge, eEdge);
	    int jPos=PositionVect(ListEdge, fEdge);
	    std::vector<int> eVal=MCEwork(iPos,jPos);
	    if (IsFirst) {
	      CommonVal=eVal;
	      IsFirst=false;
	    }
	    else {
	      if (CommonVal != eVal)
		return false;
	    }
	  }

      }
    return true;
  };
  auto GraphIntersectionOne=[&ListEdge](std::vector<std::vector<int> > const& Cspec, std::vector<std::vector<int> > const& LEdge, MyMatrix<std::vector<int> > const & MCEwork) -> GraphBitset {
    int nbComp=Cspec.size();
    GraphBitset eGR(nbComp);
    for (int iComp=0; iComp<nbComp-1; iComp++)
      for (int jComp=iComp+1; jComp<nbComp; jComp++) {
	std::vector<int> eComp=Cspec[iComp];
	std::vector<int> fComp=Cspec[jComp];
	std::vector<int> eEdge=LEdge[eComp[0]];
	std::vector<int> fEdge=LEdge[fComp[0]];
	int iPos=PositionVect(ListEdge, eEdge);
	int jPos=PositionVect(ListEdge, fEdge);
	std::vector<int> eVal=MCEwork(iPos,jPos);
	if (eVal == std::vector<int>({1}) ) {
	  eGR.AddAdjacent(iComp,jComp);
	  eGR.AddAdjacent(jComp,iComp);
	};
      }
    return eGR;
  };
  auto ExtendListLabel=[](std::vector<std::vector<int> > const& Cspec, std::vector<std::vector<int> > const& PartListLabel) -> std::vector<std::vector<int> > {
    int nbVert=0;
    for (auto & eList : Cspec)
      for (auto & eVal : eList)
	if (eVal > nbVert)
	  nbVert=eVal;
    nbVert++;
    //
    std::vector<std::vector<int> > ListLabel(nbVert);
    int nbComp=Cspec.size();
    for (int iComp=0; iComp<nbComp; iComp++) {
      std::vector<int> eLabel=PartListLabel[iComp];
      for (auto& iVert : Cspec[iComp])
	ListLabel[iVert]=eLabel;
    }
    return ListLabel;
  };
  /*
  auto GetMCEfromEmbedding=[&nbEdge, &ListEdge](MyMatrix<int> const& eEmbed) -> MyMatrix<std::vector<int> > {
    int nbLabel=eEmbed.cols();
    MyMatrix<std::vector<int> > MCEret(nbEdge, nbEdge);
    for (int iEdge=0; iEdge<nbEdge-1; iEdge++)
      for (int jEdge=iEdge+1; jEdge<nbEdge; jEdge++) {
	std::vector<int> eEdge=ListEdge[iEdge];
	std::vector<int> fEdge=ListEdge[jEdge];
	MyVector<int> VEC1(nbLabel);
	MyVector<int> VEC2(nbLabel);
	for (int i=0; i<nbLabel; i++) {
	  int eVal1=(eEmbed(eEdge[0], i) - eEmbed(eEdge[1],i)) % 2;
	  int eVal2=(eEmbed(fEdge[0], i) - eEmbed(fEdge[1],i)) % 2;
	  VEC1(i)=eVal1;
	  VEC2(i)=eVal2;
	}
	int eDist=0;
	for (int i=0; i<nbLabel; i++)
	  if (VEC1(i) != VEC2(i))
	    eDist++;
	int eVal=(4-eDist)/2;
	MCEret(iEdge,jEdge)={eVal};
	MCEret(jEdge,iEdge)={eVal};
      }
    return MCEret;
    };*/
  auto CheckEmbedding=[&](MyMatrix<int> const& eEmbed) -> bool {
    int TheDim=eEmbed.cols();
    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++) {
	std::vector<int> V1(TheDim);
	std::vector<int> V2(TheDim);
	for (int iC=0; iC<TheDim; iC++) {
	  V1[iC]=eEmbed(i,iC);
	  V2[iC]=eEmbed(j,iC);
	}
	int dist1=M(i,j);
	int dist2=L1_distance(V1, V2);
	if (dist1 <= s) {
	  if (2*dist1 != dist2) {
	    std::cerr << "Fail embedding at i=" << i << " j=" << j << "\n";
	    return false;
	  }
	}
      }
    return true;
  };
  if (DiamGraph == s) {
    std::cerr << "Direct embed, step 0\n";
    std::vector<std::vector<int> > ET=SpanningTree(GR);
    std::cerr << "Direct embed, step 1\n";
    GraphBitset Aspec=GetGraphIdentityEdge(ET, MCE);
    std::cerr << "Direct embed, step 2\n";
    std::vector<std::vector<int> > Cspec=ConnectedComponents_set(Aspec);
    std::cerr << "Direct embed, step 3\n";
    for (std::vector<int> & eList : Cspec)
      if (!IsClique(Aspec, eList)) {
	std::cerr << "Not a clique. So no embedding\n";
	return {};
      }
    std::cerr << "Direct embed, step 4\n";
    if (!IsCoherentEdgeValues(Cspec, ET, MCE)) {
      std::cerr << "Non-coherency of values\n";
      return {};
    }
    std::cerr << "Direct embed, step 5\n";
    GraphBitset Dspec=GraphIntersectionOne(Cspec, ET, MCE);
    std::cerr << "Direct embed, step 6\n";
    std::vector<std::vector<int> > PartListLabel=InverseLineGraph(Dspec);
    std::cerr << "Direct embed, step 7\n";
    if (PartListLabel[0][0] == -1) {
      std::cerr << "Not a line graph\n";
      return {};
    }
    std::cerr << "Direct embed, step 8\n";
    std::vector<std::vector<int> > ListLabel=ExtendListLabel(Cspec, PartListLabel);
    std::cerr << "Direct embed, step 9\n";
    MyMatrix<int> Embedding=CreateEmbedding(0, ET, ListLabel);
    bool res=CheckEmbedding(Embedding);
    if (!res) {
      std::cerr << "Maybe some bug in the code\n";
      throw TerminalException{1};
    }
    std::vector<MyMatrix<int> > ListEmbedding{Embedding};
    std::cerr << "|ListEmbedding|=" << ListEmbedding.size() << "\n";
    return ListEmbedding;
  }
  auto RefinementBlock1=[&ListEdge,&GetGraphIdentityEdge](MyMatrix<std::vector<int> > & MCEwork, int &nbOper) -> bool {
    /*
    std::cerr << "ListEdge=\n";
    int nbEdge=ListEdge.size();
    for (int iEdge=0; iEdge<nbEdge; iEdge++) {
      std::cerr << "iEdge=" << iEdge << " e=" << ListEdge[iEdge][0] << "," << ListEdge[iEdge][1] << "\n";
    }
    for (int iEdge=0; iEdge<nbEdge; iEdge++)
      for (int jEdge=0; jEdge<nbEdge; jEdge++) {
	std::cerr << "i/jEdge=" << iEdge << "/" << jEdge << " MCE=";
	for (auto & eVal : MCE(iEdge,jEdge))
	  std::cerr << " " << eVal;
	std::cerr << "\n";
	}*/
    GraphBitset Aspec=GetGraphIdentityEdge(ListEdge, MCEwork);
    //    std::cerr << "Aspec=\n";
    //    GRAPH_PrintOutput(std::cerr, Aspec);
    std::vector<std::vector<int> > Cspec=ConnectedComponents_set(Aspec);
    //    std::cerr << "|Cspec|=" << Cspec.size() << "\n";
    //    std::cerr << "Cspec=\n";
    //    Print_VectorVector(std::cerr, Cspec);
    for (auto & eConn : Cspec) {
      int lenConn=eConn.size();
      //      std::cerr << "|eConn|=" << lenConn << "\n";
      for (int i=0; i<lenConn-1; i++)
	for (int j=i+1; j<lenConn; j++) {
	  int iEdge=eConn[i];
	  int jEdge=eConn[j];
	  std::vector<int> eVal=MCEwork(iEdge, jEdge);
	  int pos=PositionVect(eVal,2);
	  if (pos == -1) {
	    return false;
	  }
	  if (eVal != std::vector<int>({2})) {
	    MCEwork(iEdge,jEdge)={2};
	    MCEwork(jEdge,iEdge)={2};
	    nbOper++;
	  }
	}
    }
    return true;
  };
  auto RefinementBlock2=[&ListEdge,&GetGraphIdentityEdge](MyMatrix<std::vector<int> > & MCEwork, int & nbOper) -> bool {
    GraphBitset Aspec=GetGraphIdentityEdge(ListEdge, MCEwork);
    std::vector<std::vector<int> > Cspec=ConnectedComponents_set(Aspec);
    int nbConn=Cspec.size();
    for (int iConn=0; iConn<nbConn-1; iConn++)
      for (int jConn=iConn+1; jConn<nbConn; jConn++) {
	std::vector<int> ListPoss{0,1,2};
	std::vector<int> eConn=Cspec[iConn];
	std::vector<int> fConn=Cspec[jConn];
	int iSize=eConn.size();
	int jSize=fConn.size();
	for (int i=0; i<iSize; i++)
	  for (int j=0; j<jSize; j++) {
	    int iEdge=eConn[i];
	    int jEdge=fConn[j];
	    std::vector<int> eVal=MCEwork(iEdge,jEdge);
	    ListPoss=IntersectionVect(ListPoss, eVal);
	  }
	if (ListPoss.size() == 0) {
	  return false;
	}
	for (int i=0; i<iSize; i++)
	  for (int j=0; j<jSize; j++) {
	    int iEdge=eConn[i];
	    int jEdge=fConn[j];
	    std::vector<int> eVal=MCEwork(iEdge,jEdge);
	    if (eVal != ListPoss) {
	      MCEwork(iEdge,jEdge)=ListPoss;
	      MCEwork(jEdge,iEdge)=ListPoss;
	      nbOper++;
	    }
	  }
      }
    return true;
  };
  // Iterative refinement triangles
  // If M(a,b) = M(a,c) = M(b,c) = 1
  // TYPE 1:
  //   We have a=UV, b=UW, c=VW
  //   If we have a fourth edge with M(a,d)=1 then
  //   d=Ux with x not equal to V.
  //   ---If x is not V or W then M(b,d)=1 and M(b,c)=0 and there can be more than one such
  //   ---If x=W then M(b,d)=2 and M(c,d)=1
  //   So possible configurations: [1,1,0] and [1,1,2] or [0,0,0]
  // TYPE 2:
  //   We have a=UR, b=US, c=UT
  //   If we have a fourth edge with M(a,d)=1 then
  //   If d=Ux with x not R,S,T then pattern is [1,1,1]
  //   If d=Rx with x not U,S,T then pattern is [1,0,0]
  //   If d=RS or combination   then pattern is [1,1,0]
  //   Other possible patterns: [2,1,1] and [0,0,0]
  // In summary:
  //    TYPE 1: [1,1,0],                   [2,1,1], [0,0,0]
  //    TYPE 2: [1,1,1], [1,0,0], [1,1,0], [2,1,1], [0,0,0]
  //                    Pattern [1,1,0] can occur at most 3 times.
  auto RefinementTriangles=[&ListEdge,&GetGraphIdentityEdge,&GraphIntersectionOne,&SetOne,&SetZero,&SetTwo](MyMatrix<std::vector<int> > & MCEwork, int & nbOper) -> bool {
    GraphBitset Aspec=GetGraphIdentityEdge(ListEdge, MCEwork);
    //    std::cerr << "|Aspec|=" << Aspec.GetNbVert() << "\n";
    //    GRAPH_PrintOutput(std::cerr, Aspec);
    std::vector<std::vector<int> > Cspec=ConnectedComponents_set(Aspec);
    //    std::cerr << "|Cspec|=" << Cspec.size() << "\n";
    GraphBitset Dspec=GraphIntersectionOne(Cspec, ListEdge, MCEwork);
    int nbVertD=Dspec.GetNbVert();
    std::vector<int> VertSet(nbVertD);
    for (int iVert=0; iVert<nbVertD; iVert++)
      VertSet[iVert]=iVert;
    std::vector<std::vector<int> > ListTriple=EnumerationClique(Dspec, 3);
    for (auto & eTriple : ListTriple) {
      bool Type1feasible=true;
      bool Type2feasible=true;
      int nb110=0;
      for (auto & d : DifferenceVect(VertSet, eTriple) ) {
	int dEdge=Cspec[d][0];
	int nb1=0;
	int nb0=0;
	for (int i=0; i<3; i++) {
	  int jEdge=Cspec[eTriple[i]][0];
	  if (MCEwork(dEdge, jEdge) == SetOne)
	    nb1++;
	  if (MCEwork(dEdge, jEdge) == SetZero)
	    nb0++;
	}
	if (nb1 == 2 && nb0 == 1)
	  nb110++;
      }
      if (nb110 > 3)
	Type2feasible=false;
      for (int i=0; i<3; i++) {
	int iNext=NextIdx(3,i);
	int iPrev=PrevIdx(3,i);
	int iEdge=Cspec[eTriple[i]][0];
	int iEdgeNext=Cspec[eTriple[iNext]][0];
	int iEdgePrev=Cspec[eTriple[iPrev]][0];
	/*
	std::cerr << "iEdge=" << iEdge << " iEdgeNext=" << iEdgeNext << " iEdgePrev=" << iEdgePrev << "\n";
	std::cerr << "eEdge=" << ListEdge[iEdge][0] << " , " << ListEdge[iEdge][1] << "\n";
	std::cerr << "eEdgePrev=" << ListEdge[iEdgePrev][0] << " , " << ListEdge[iEdgePrev][1] << "\n";
	std::cerr << "eEdgeNext=" << ListEdge[iEdgeNext][0] << " , " << ListEdge[iEdgeNext][1] << "\n";
	std::cerr << "MCE01=";
	WriteVectorInt_GAP(std::cerr, MCE(iEdge, iEdgeNext));
	std::cerr << " MCE02=";
	WriteVectorInt_GAP(std::cerr, MCE(iEdge, iEdgePrev));
	std::cerr << " MCE12=";
	WriteVectorInt_GAP(std::cerr, MCE(iEdgeNext, iEdgePrev));
	std::cerr << "\n";
	std::cerr << "----------------------------------------------------------\n";
	*/
	for (auto & d : DifferenceVect(VertSet, eTriple) ) {
	  int dEdge=Cspec[d][0];
	  /*
	  std::cerr << "dEdge=" << dEdge << "\n";
	  std::cerr << "dEdge=" << ListEdge[dEdge][0] << " , " << ListEdge[dEdge][1] << "\n";
	  std::cerr << "MCEd0=";
	  WriteVectorInt_GAP(std::cerr, MCE(dEdge, iEdgePrev));
	  std::cerr << "\n";
	  std::cerr << "MCEd1=";
	  WriteVectorInt_GAP(std::cerr, MCE(dEdge, iEdgeNext));
	  std::cerr << "\n";
	  std::cerr << "MCEd2=";
	  WriteVectorInt_GAP(std::cerr, MCE(dEdge, iEdge));
	  std::cerr << "\n";*/
	  //
	  //	  int pos0=PositionVect(MCE(dEdge,iEdge), 0);
	  int pos1=PositionVect(MCEwork(dEdge,iEdge), 1);
	  int pos2=PositionVect(MCEwork(dEdge,iEdge), 2);
	  if (MCEwork(dEdge, iEdgeNext) == SetOne && MCEwork(dEdge, iEdgePrev) == SetOne) {
	    if (pos1 != -1 && !Type2feasible) {
	      std::vector<int> NewVal=DifferenceVect(MCEwork(dEdge,iEdge), SetOne);
	      MCEwork(dEdge, iEdge)=NewVal;
	      MCEwork(iEdge, dEdge)=NewVal;
	      nbOper++;
	    }
	    if (pos2 == -1) { // value is necessarily 1 or 0 and this forbids Type1
	      //	      std::cerr << "Now Type1feasible=false\n";
	      Type1feasible=false;
	    }
	  }
	  for (int j=0; j<2; j++) {
	    int iEdge1, iEdge2;
	    if (j == 0) {
	      iEdge1=iEdgePrev;
	      iEdge2=iEdgeNext;
	    }
	    else {
	      iEdge1=iEdgeNext;
	      iEdge2=iEdgePrev;
	    }
	    if (MCEwork(dEdge, iEdge1) == SetOne && MCEwork(dEdge,iEdge2) == SetZero) {
	      if (pos2 != -1) {
		std::vector<int> NewVal=DifferenceVect(MCEwork(dEdge,iEdge), SetTwo);
		MCEwork(dEdge, iEdge)=NewVal;
		MCEwork(iEdge, dEdge)=NewVal;
		nbOper++;
	      }
	      if (pos1 == -1 && !Type2feasible) {
		std::cerr << "Finding false 1\n";
		return false;
	      }
	      if (MCEwork(dEdge,iEdge) != SetOne && !Type2feasible) {
		MCEwork(dEdge,iEdge)=SetOne;
		MCEwork(iEdge,dEdge)=SetOne;
		nbOper++;
	      }
	    }
	  }
	  if (MCEwork(dEdge, iEdgePrev) == SetZero && MCEwork(dEdge,iEdgeNext) == SetZero) {
	    if (pos2 != -1) {
	      std::vector<int> NewVal=DifferenceVect(MCEwork(dEdge,iEdge), SetTwo);
	      MCEwork(dEdge, iEdge)=NewVal;
	      MCEwork(iEdge, dEdge)=NewVal;
	      nbOper++;
	    }
	    if (!Type2feasible && pos1 != -1) {
	      std::vector<int> NewVal=DifferenceVect(MCEwork(dEdge,iEdge), SetOne);
	      MCEwork(dEdge, iEdge)=NewVal;
	      MCEwork(iEdge, dEdge)=NewVal;
	      nbOper++;
	    }
	    if (MCEwork(dEdge,iEdge) == SetOne) {
	      //	      std::cerr << "Now Type1feasible=false\n";
	      Type1feasible=false;
	    }
	  }
	}
      }
      if (!Type1feasible && !Type2feasible)
	return false;
    }
    return true;
  };
  auto IsAprioriFeasible=[&nbEdge](MyMatrix<std::vector<int> > const& eMat) -> bool {
    for (int iEdge=0; iEdge<nbEdge-1; iEdge++)
      for (int jEdge=iEdge+1; jEdge<nbEdge; jEdge++) {
	int siz=eMat(iEdge,jEdge).size();
	if (siz == 0)
	  return false;
      }
    return true;
  };
  auto IterativeRefinementMCE=[&](MyMatrix<std::vector<int> > & MCEwork) -> bool {
    bool test;
    while(true) {
      int nbOper=0;
      test=RefinementBlock1(MCEwork, nbOper);
      std::cerr << "RefinementBlock1 test=" << test << " nbOper=" << nbOper << "\n";
      if (!test)
	return false;
      test=RefinementBlock2(MCEwork, nbOper);
      std::cerr << "RefinementBlock2 test=" << test << " nbOper=" << nbOper << "\n";
      if (!test)
	return false;
      test=RefinementTriangles(MCEwork, nbOper);
      std::cerr << "RefinementTriangles test=" << test << " nbOper=" << nbOper << "\n";
      if (!test)
	return false;
      test=IsAprioriFeasible(MCEwork);
      std::cerr << "IsAprioriFeasible test=" << test << "\n";
      if (!test)
	return false;
      if (nbOper == 0)
	break;
    }
    return true;
  };
  auto GetFirstNonZero=[&nbEdge](MyMatrix<std::vector<int> > const& MCEwork) -> std::vector<int> {
    for (int iEdge=0; iEdge<nbEdge-1; iEdge++)
      for (int jEdge=iEdge+1; jEdge<nbEdge; jEdge++) {
	int siz=MCEwork(iEdge,jEdge).size();
	if (siz > 1)
	  return {iEdge,jEdge};
      }
    std::cerr << "Fail to find an index\n";
    throw TerminalException{1};
  };
  auto IsFinalState=[&nbEdge](MyMatrix<std::vector<int> > const& MCEwork) -> bool {
    for (int iEdge=0; iEdge<nbEdge-1; iEdge++)
      for (int jEdge=iEdge+1; jEdge<nbEdge; jEdge++) {
	int siz=MCEwork(iEdge,jEdge).size();
	if (siz > 1)
	  return false;
      }
    return true;
  };
  auto GetEmbeddingMatrix=[&nbEdge,&ListEdge,&GetGraphIdentityEdge,&GraphIntersectionOne,&ExtendListLabel,&CheckEmbedding](MyMatrix<std::vector<int> > const& MCEwork, MyMatrix<int> & FinalEmbed) -> bool {
    std::cerr << "GetEmbeddingMatrix, step 1\n";
    /*
    for (int iEdge=0; iEdge<nbEdge; iEdge++)
      for (int jEdge=0; jEdge<nbEdge; jEdge++) {
	std::cerr << "i/jEdge=" << iEdge << "/" << jEdge << " MCE=";
	for (auto & eVal : MCE(iEdge,jEdge))
	  std::cerr << " " << eVal;
	std::cerr << "\n";
	}*/
    std::cerr << "GetEmbeddingMatrix, step 2\n";
    GraphBitset Aspec=GetGraphIdentityEdge(ListEdge, MCEwork);
    std::cerr << "GetEmbeddingMatrix, step 3\n";
    //    std::cerr << "Aspec=";
    //    GRAPH_PrintOutput(std::cerr, Aspec);
    std::vector<std::vector<int> > Cspec=ConnectedComponents_set(Aspec);
    std::cerr << "GetEmbeddingMatrix, step 4 |Cspec|=" << Cspec.size() << "\n";
    //    std::cerr << "Cspec=\n";
    //    Print_VectorVector(std::cerr, Cspec);
    GraphBitset Dspec=GraphIntersectionOne(Cspec, ListEdge, MCEwork);
    //    std::cerr << "--------------------------------------------\n";
    //    std::cerr << "GetEmbeddingMatrix, step 5\n";
    //    std::cerr << "Dspec=";
    //    GRAPH_PrintOutput(std::cerr, Dspec);
    //    GRAPH_PrintOutputGAP(std::cerr, Dspec);
    std::vector<std::vector<int> > PartListLabel=InverseLineGraph(Dspec);
    std::cerr << "GetEmbeddingMatrix, step 6\n";
    //    std::cerr << "PartListLabel=\n";
    //    Print_VectorVector(std::cerr, PartListLabel);
    //    std::cerr << "--------------------------------------------\n";
    if (PartListLabel[0][0] == -1) {
      std::cerr << "Embedding failure\n";
      return false;
    }
    std::vector<std::vector<int> > ListLabel=ExtendListLabel(Cspec, PartListLabel);
    std::cerr << "GetEmbeddingMatrix, step 7\n";
    //    std::cerr << "ListLabel=";
    //    Print_VectorVector(std::cerr, ListLabel);
    FinalEmbed=CreateEmbedding(0, ListEdge, ListLabel);
    bool res=CheckEmbedding(FinalEmbed);
    return res;
  };
  struct OneLevel {
    std::vector<int> PairEdge;
    std::vector<int> ListPoss;
    int idx;
    MyMatrix<std::vector<int> > MCE;
  };
  std::vector<OneLevel> TheTree;
  auto GoUpNextInTree=[&]() -> bool {
    int idx;
    while(true) {
      if (TheTree.size() == 0)
	return false;
      int len=TheTree.size();
      std::cerr << "TheTree len=" << len << "\n";
      int lenPoss=TheTree[len-1].ListPoss.size();
      std::cerr << "lenPoss=" << lenPoss << "\n";
      if (len == 1) {
	std::cerr << "ListPoss=";
	for (auto & eVal : TheTree[len-1].ListPoss)
	  std::cerr << " " << eVal;
	std::cerr << "\n";
      }
      idx=TheTree[len-1].idx;
      std::cerr << "idx=" << idx << "\n";
      if (idx == lenPoss -1) {
	TheTree.pop_back();
      }
      else {
	idx++;
	int iEdge=TheTree[len-1].PairEdge[0];
	int jEdge=TheTree[len-1].PairEdge[1];
	TheTree[len-1].idx=idx;
	MyMatrix<std::vector<int> > MCEnew;
	if (len == 1) {
	  MCEnew=MCE;
	}
	else {
	  MCEnew=TheTree[len-2].MCE;
	}
	if (len == 1) {
	  int eVal=TheTree[len-1].ListPoss[idx];
	  std::cerr << "1: Now operating with ListPoss=" << eVal << "\n";
	}
	int eVal=TheTree[len-1].ListPoss[idx];
	MCEnew(iEdge,jEdge)={eVal};
	bool test=IterativeRefinementMCE(MCEnew);
	TheTree[len-1].MCE=MCEnew;
	if (test)
	  return true;
      }
    }
  };
  auto GetLatestMCE=[&TheTree,&MCE]() -> MyMatrix<std::vector<int> > {
    int len=TheTree.size();
    if (len == 0)
      return MCE;
    return TheTree[len-1].MCE;
  };
  auto NextInTree=[&]() -> bool {
    MyMatrix<std::vector<int> > MCEstart=GetLatestMCE();
    bool test=IsFinalState(MCEstart);
    //    std::cerr << "test=" << test << "\n";
    if (test)
      return GoUpNextInTree();
    std::vector<int> ePairEdge=GetFirstNonZero(MCEstart);
    int iEdge=ePairEdge[0];
    int jEdge=ePairEdge[1];
    std::cerr << "Assign iEdge=" << iEdge << " jEdge=" << jEdge << "\n";
    std::vector<int> ListPoss=MCEstart(iEdge,jEdge);
    int idx=0;
    MCEstart(iEdge,jEdge)={ListPoss[idx]};
    MCEstart(jEdge,iEdge)={ListPoss[idx]};
    TheTree.push_back({ePairEdge, ListPoss, idx, MCEstart});
    int len=TheTree.size();
    bool test2=IterativeRefinementMCE(TheTree[len-1].MCE);
    if (!test2)
      return GoUpNextInTree();
    return true;
  };
  bool test=IterativeRefinementMCE(MCE);
  std::cerr << "IterativeRefinementMCE test=" << test << "\n";
  if (!test)
    return {};
  std::vector<MyMatrix<int> > ListEmbedding;
  int nbChoice=nbEdge*(nbEdge-1)/2;
  while(true) {
    MyMatrix<std::vector<int> > MCEstart=GetLatestMCE();
    bool testF=IsFinalState(MCEstart);
    std::cerr << "IsFinalState test=" << testF << "\n";
    if (testF) {
      MyMatrix<int> FinalEmbed;
      bool test2=GetEmbeddingMatrix(MCEstart, FinalEmbed);
      std::cerr << "test2=" << test2 << "\n";
      if (test2)
	ListEmbedding.push_back(FinalEmbed);
    }
    std::cerr << "Passing by NextInTree nbEmbedding=" << ListEmbedding.size() << "\n";
    std::cerr << "|TheTree|=" << TheTree.size() << " nbChoice=" << nbChoice << " iter=" << iter << "\n";
    if (iter > MaxIter)
      break;
    bool test2=NextInTree();
    if (!test2)
      break;
    iter++;
  }
  std::cerr << "|ListEmbedding|=" << ListEmbedding.size() << "\n";
  return ListEmbedding;
}


#endif

#ifndef DEFINE_BLOCK_SYSTEM_H
#define DEFINE_BLOCK_SYSTEM_H

#include "MAT_Matrix.h"
#include "Face_basic.h"
#include "GRAPH_GraphicalBasic.h"
#include "PermGroup.h"


namespace gap {

template<typename Telt>
std::pair<std::vector<std::vector<int>>,std::vector<Face>> Blocks_Kernel(std::vector<Telt> const& ListGen, std::vector<int> const& Omega, std::vector<int> const& seed)
{
  int nbMax=VectorMax(Omega);
  int n=Omega.size();
  std::vector<int> OmegaRev(nbMax, -1);
  for (int u=0; u<n; u++) {
    int ePt=Omega[u];
    OmegaRev[ePt] = u;
  }
  std::vector<Telt> ListGenRed;
  for (auto & eGen : ListGen) {
    std::vector<int> eList(n);
    for (int u=0; u<n; u++) {
      int ePt=Omega[u];
      int ePtImg=eGen.at(ePt);
      int uImg=OmegaRev[ePtImg];
      eList[u] = uImg;
    }
    Telt eGenRed(eList);
    ListGenRed.push_back(eGenRed);
  }
  // Building the orbit
  Face eFace(n);
  for (auto & ePt : seed) {
    int u=OmegaRev[ePt];
    eFace[u]=1;
  }
  std::function<Face(Face const&,Telt const&)> act=[&](Face const& x, Telt const& eElt) -> Face {
    Face eImg(n);
    int siz=x.count();
    int ePt=x.find_first();
    for (int u=0; u<siz; u++) {
      int ePtImg=eElt.at(ePt);
      eImg[ePtImg]=1;
      ePt = x.find_next(ePt);
    }
    return eImg;
  };
  std::vector<Face> SeedOrbit = Orbit(ListGenRed, eFace, act);
  // Building the adjacency graph
  auto IsIntersecting=[&](Face const& face1, Face const& face2) -> bool {
    int siz=face1.count();
    int ePt=face1.find_first();
    for (int u=0; u<siz; u++) {
      if (face2[ePt] == 1)
	return true;
      ePt = face1.find_next(ePt);
    }
    return false;
  };
  int sizOrb=SeedOrbit.size();
  std::vector<int> ListA;
  for (int u=0; u<sizOrb; u++)
    for (int v=u+1; v<sizOrb; v++) {
      if (IsIntersecting(SeedOrbit[u], SeedOrbit[v])) {
	ListA.push_back(u);
	ListA.push_back(v);
      }
    }
  int nbEdge=ListA.size() / 2;
  MyMatrix<int> ListEdge(nbEdge,2);
  for (int iEdge=0; iEdge<nbEdge; iEdge++) {
    ListEdge(iEdge,0) = ListA[2*iEdge];
    ListEdge(iEdge,1) = ListA[2*iEdge+1];
  }
  GraphSparseImmutable eGR(ListEdge, sizOrb);
  std::vector<std::vector<int>> ListConn = ConnectedComponents_set(eGR);
  std::vector<std::vector<int>> ListBlocks;
  for (auto & eConn : ListConn) {
    Face eFaceComb(n);
    for (auto & iElt : eConn) {
      Face eFace = SeedOrbit[iElt];
      int siz=eFace.count();
      int ePt=eFace.find_first();
      for (int u=0; u<siz; u++) {
	eFaceComb[ePt]=1;
	ePt = eFace.find_next(ePt);
      }
    }
    std::vector<int> eBlock;
    int blkSiz=eFaceComb.count();
    int ePt=eFaceComb.find_first();
    for (int u=0; u<blkSiz; u++) {
      eBlock.push_back(Omega[ePt]);
      ePt = eFaceComb.find_next(ePt);
    }
    ListBlocks.push_back(eBlock);
  }
  return {ListBlocks,SeedOrbit};
}


template<typename Telt>
std::vector<std::vector<int>> Blocks_from_seed(std::vector<Telt> const& ListGen, std::vector<int> const& Omega, std::vector<int> const& seed)
{
  return Blocks_Kernel(ListGen, Omega, seed).first;
}



// Maybe we can do better by using std::vector for the edges, but
// that is a secondary optimization
template<typename Telt>
std::vector<std::vector<int>> Blocks_without_seed(std::vector<Telt> const& ListGen, std::vector<int> const& Omega)
{
  int n=Omega.size();
  int nbTotalEdge=n * (n-1) / 2;
  Face UsedEdge(nbTotalEdge);
  int ThePos=0;
  auto UpdatePosition=[&]() -> bool {
    while(true) {
      if (UsedEdge[ThePos] == 0)
	return true;
      ThePos++;
      if (ThePos == nbTotalEdge)
	break;
    }
    return false;
  };
  auto GetEdgeFromPosition=[&](int const& pos) -> std::vector<int> {
    int posTot=pos;
    int firstElement=0;
    while(true) {
      int allowedSiz = n - 1 - firstElement;
      if (posTot < allowedSiz) {
	int secondElement=firstElement + 1 + posTot;
	return {firstElement,secondElement};
      }
      posTot -= allowedSiz;
      firstElement++;
    }
  };
  auto GetPosFromEdge=[&](std::vector<int> const& eEdge) -> int {
    int eFirst=eEdge[0];
    int eSecond=eEdge[1];
    int pos=0;
    for (int u=0; u<eFirst; u++) {
      int allowedSiz = n - 1 - eFirst;
      pos += allowedSiz;
    }
    return pos + eSecond - 1 - eFirst;
  };
  while(true) {
    std::vector<int> eEdge=GetEdgeFromPosition(ThePos);
    std::pair<std::vector<std::vector<int>>,std::vector<Face>> ePair = Blocks_Kernel(ListGen, Omega, eEdge);
    if (ePair.first.size() > 1) {
      return ePair.first;
    }
    for (auto & eEdgeFace : ePair.second) {
      std::vector<int> eEdge = FaceToVector(eEdgeFace);
      int pos=GetPosFromEdge(eEdge);
      UsedEdge[pos]=1;
    }
    bool test = UpdatePosition();
    if (!test)
      break;
  }
  return {Omega};
}

template<typename Telt>
std::function<Telt(Telt const&)> MapElementToSetPlusBlocks(std::vector<std::vector<int>> const& blks, int const& n)
{
  int TheMax=0;
  for (auto & eBlock : blks)
    for (auto & ePt : eBlock)
      if (ePt > TheMax)
	TheMax = ePt;
  std::vector<int> VectStatus(TheMax);
  int nbBlock=blks.size();
  for (int iBlk=0; iBlk<nbBlock; iBlk++)
    for (auto & ePt : blks[iBlk])
      VectStatus[ePt] = iBlk;
  std::function<Telt(Telt const&)> f=[=](Telt const& u) -> Telt {
    std::vector<int> eList(n+nbBlock);
    for (int i=0; i<n; i++)
      eList[i] = u.at(i);
    for (int iBlk=0; iBlk<nbBlock; iBlk++) {
      int ePt=blks[iBlk][0];
      int ePtImg=u.at(ePt);
      int iBlkImg=VectStatus[ePtImg];
      eList[n+iBlk] = n+iBlkImg;
    }
    return Telt(eList);
  };
  return f;
}

std::vector<int> GetBlock(std::vector<std::vector<int>> const& ListBlock, int const& ePt)
{
  for (auto & eBlock : ListBlock)
    if (PositionVect(eBlock, ePt) != -1)
      return eBlock;
  return {-1};
}


int GetNonTrivialPointInBlock(std::vector<int> const& eBlock, int const& ePt)
{
  for (auto & uPt : eBlock)
    if (uPt != ePt)
      return uPt;
  return -1;
}

}

#endif


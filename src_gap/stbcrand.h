#ifndef DEFINE_INCLUDE_STAB_CHAIN_RANDOM
#define DEFINE_INCLUDE_STAB_CHAIN_RANDOM

// Functions return type
// SCRSift           : returns permutation
// VerifyStabilizer  : returns permutation
// SCRStrongGenTest  : returns permutation
// SCRStrongGenTest2 : returns permutation
// VerifySGS         : returns permutation or result



#include "StabChain.h"
#include "BlockSystem.h"
#include "pseudorandom.h"
#include "COMB_Vectors.h"

namespace gap {

  
template<typename Telt>
void SCRExtend(std::vector<Telt> & labels, std::vector<int> & orb, std::vector<int> & transversal, std::vector<Telt> & treegen, std::vector<Telt> & treegeninv, int & len)
{
  int previous=len;
  len=orb.size();
  int nbTree=treegen.size();
  std::vector<int> ListPosGen(nbTree);
  std::vector<int> ListPosGeninv(nbTree);
  for (int iTree=0; iTree<nbTree; iTree++) {
    ListPosGen[iTree] = GetLabelIndex(labels, treegen[iTree]);
    ListPosGeninv[iTree] = GetLabelIndex(labels, treegeninv[iTree]);
  }
  for (int i=previous; i<len; i++)
    for (int j=0; j<nbTree; j++) {
      int img1=PowAct(orb[i], treegen[j]);
      int img2=PowAct(orb[i], treegeninv[j]);
      if (transversal[img1] == -1) {
	transversal[img1]=ListPosGeninv[j];
	orb.push_back(img1);
      }
      if (transversal[img2] == -1) {
	transversal[img2]=ListPosGen[j];
	orb.push_back(img2);
      }
    }
}
  

struct noticeType {
  bool res;
  int i;
  int j;
};

template<typename Telt>
noticeType SCRNotice_A(std::vector<int> const& orb, std::vector<int> const& transversal, std::vector<Telt> const& genlist)
{
  int siz1=orb.size();
  int siz2=genlist.size();
  for (int i=0; i<siz1; i++)
    for (int j=0; j<siz2; j++) {
      int ePt=PowAct(orb[i], genlist[j]);
      if (transversal[ePt] == -1)
	return {false, i, j};
    }
  return {true, -1,-1};
}


template<typename Telt>
noticeType SCRNotice_B(std::vector<int> const& orb, std::vector<int> const& transversal, std::vector<int> const& genlistIdx, std::vector<Telt> const& labels)
{
  int siz1=orb.size();
  int siz2=genlistIdx.size();
  for (int i=0; i<siz1; i++)
    for (int j=0; j<siz2; j++) {
      int posGen=genlistIdx[j];
      int ePt=PowAct(orb[i], labels[posGen]);
      if (transversal[ePt] == -1)
	return {false, i, j};
    }
  return {true, -1,-1};
}



 
int QuoInt(int const& a, int const& b)
{
  int res=a % b;
  return (a - res)/b;
}

template<typename T>
T LogInt(T const& TheNB, T const& expo)
{
  int result=0;
  T ThePow=1;
  while(true) {
    ThePow *= expo;
    if (ThePow > TheNB)
      return result;
    result++;
  }
}

 
Face SCRRandomString(int const& n)
{
  InfoPseudoRandom* R = GetPseudoRandom();
  Face string(n);
  int k=QuoInt(n-1, 28);
  for (int i=0; i<=k; i++) {
    RandomShift(R);
    Face f=Extract01vector(R);
    for (int j=0; j<28; j++) {
      int pos=28*i + j;
      if (pos < n)
	string[28*i + j] = f[j];
    }
  }
  return string;
}


template<typename Telt>
Telt SCRRandomPerm(int const& d )
{
  std::vector<int> rnd = ClosedInterval(0, d);
  for (int i=0; i<d; i++) {
    int u=d+1-i;
    int k = rand() % u;
    if (u != k)
      std::swap(rnd[u], rnd[k]);
  }
  return Telt(rnd);
}



 
template<typename Telt>
std::vector<Telt> CosetRepAsWord(std::vector<Telt> const& labels, int const& x, int const& y, std::vector<int> const& transversal)
{
  if (transversal[y] == -1)
    return {};
  int pnt=y;
  std::vector<Telt> word;
  while(true) {
    if (pnt == x)
      break;
    int pos=transversal[pnt];
    Telt eElt=labels[pos];
    word.push_back(eElt);
    pnt=PowAct(pnt, eElt);
  }
  return word;
}


template<typename Telt>
std::vector<Telt> InverseAsWord(std::vector<Telt> const& word, std::vector<Telt> const& list, std::vector<Telt> const& inverselist)
{
  int siz=word.size();
  if (siz == 1) {
    if (word[0].isIdentity())
      return word;
  }
  int sizList=list.size();
  std::vector<Telt> inverse(siz);
  for (int i=0; i<siz; i++) {
    auto GetPosition=[&](Telt const& x) -> int {
      for (int j=0; j<sizList; j++) {
	if (list[j] == x)
	  return j;
      }
      return -1;
    };
    inverse[i]=inverselist[GetPosition(word[siz-1-i])];
  }
  return inverse;
}
 
template<typename Telt>
int ImageInWord(int const& x, std::vector<Telt> const& word)
{
  int value=x;
  for (auto & eElt : word)
    value = PowAct(value, eElt);
  return value;
}


template<typename Telt>
std::pair<std::vector<Telt>,int> SiftAsWord(StabChain<Telt> const& Stot, int const& eLev, std::vector<Telt> const& perm)
{
  int len=Stot.stabilizer.size();
  int index=0;
  std::vector<Telt> word = perm;
  for (int iLev=eLev; iLev<len; iLev++) {
    index++;
    int pnt=Stot.stabilizer[iLev].orbit[0];
    int y=ImageInWord(pnt,word);
    if (Stot.stabilizer[iLev].transversal[y] == -1)
      return {word, index};
    std::vector<Telt> coset=CosetRepAsWord(Stot.labels, pnt, y, Stot.stabilizer[iLev].transversal);
    for (auto & eElt : coset)
      word.push_back(eElt);
  }
  index=0;
  return {word,index};
}


template<typename Telt>
std::vector<Telt> RandomElmAsWord(StabChain<Telt> const& Stot, int const& eLev)
{
  std::vector<Telt> word;
  int len=Stot.stabilizer.size();
  for (int iLev=eLev; iLev<len; iLev++) {
    int sizOrb=Stot.stabilizer[iLev].orbit.size();
    int pos=rand() % sizOrb;
    int ePt=Stot.stabilizer[iLev].orbit[0];
    int fPt=Stot.stabilizer[iLev].orbit[pos];
    std::vector<Telt> coset = CosetRepAsWord(Stot.labels, ePt, fPt, Stot.stabilizer[iLev].transversal);
    word.insert(word.end(), coset.begin(), coset.end());
  }
  return word;
}
 


template<typename Telt>
Telt SCRSift(StabChain<Telt> const& Stot, int const& eLev, Telt const& g)
{
  int len=Stot.stabilizer.size();
  Telt gRet=g;
  for (int iLev=eLev; iLev<len; iLev++) {
    int bpt=Stot.stabilizer[iLev].orbit[0];
    if (Stot.stabilizer[iLev].transversal[PowAct(bpt, gRet)] == -1)
      return gRet;
    while(true) {
      int img=PowAct(bpt, gRet);
      if (bpt == img)
	break;
      int pos=Stot.stabilizer[iLev].transversal[img];
      gRet=gRet * Stot.labels[pos];
    }
  }
  return gRet;
}

/*
  param[1] = number of pairs of random subproducts from generators in
             first checking phase
  param[2] = (number of random elements from created set)/S.diam
             in first checking phase
  param[3] = number of pairs of random subproducts from generators in
             second checking phase
  param[4] = (number of random elements from created set)/S.diam
             in second checking phase
  param[5] = maximum size of orbits in  which we evaluate words on all
             points of orbit
  param[6] = minimum number of random points from orbit to plug in to check
             whether given word is identity on orbit
*/
struct paramOpt {
  int param1;
  int param2;
  int param3;
  int param4;
  int param5;
  int param6;
};


template<typename Telt>
Telt Product(std::vector<Telt> const& eList)
{
  int siz=eList.size();
  Telt retVal=eList[0];
  for (int i=1; i<siz; i++)
    retVal = retVal * eList[i];
  return retVal;
}





  
template<typename Telt>
void SCRSchTree(StabChain<Telt> & Stot, int const& eLev, std::vector<Telt> const& newgens )
{
  int n=Stot.n;
  StabLevel<Telt> & S = Stot.stabilizer[eLev];
  noticeType l= SCRNotice_A(S.orbit, S.transversal, newgens);
  if (l.res)
    return;
  int i = l.i;
  int j = l.j;
  Telt witness = newgens[j];
  while(true) {
    std::vector<Telt> word = CosetRepAsWord(Stot.labels, S.orbit[0], S.orbit[i], S.transversal);
    Telt g = Product(word);
    Telt eGen = Inverse(g) * witness;
    Telt eGenInv = Inverse(witness) * g;
    S.treegen.push_back(eGen);
    S.treegeninv.push_back(eGenInv);
    S.orbit = {S.orbit[0]};
    S.transversal = std::vector<int>(n, -1);
    S.transversal[S.orbit[0]] = GetLabelIndex(Stot.labels, Stot.identity);
    S.treedepth = 0;
    int list5=0;
    while(true) {
      if (S.treedepth >= 2*int(S.treegen.size()))
	break;
      SCRExtend(Stot.labels, S.orbit, S.transversal, S.treegen, S.treegeninv, list5);
      if (int(S.orbit.size()) == list5) {
	break;
      }
      else {
	S.treedepth++;
      }
    }
    l = SCRNotice_B(S.orbit, S.transversal, S.genlabels, Stot.labels);
    if (l.res)
      break;
    i = l.i;
    j = l.j;
    int posGen=S.genlabels[j];
    witness = Stot.labels[posGen];
  }
  S.aux  = Concatenation(S.treegen, S.treegeninv, Stot.stabilizer[eLev+1].aux);
  S.diam = S.treedepth + Stot.stabilizer[eLev+1].diam;
}


template<typename Telt>
void SCRMakeStabStrong(StabChain<Telt> & Stot, int const& eLev, std::vector<Telt> const& newgens, paramOpt const& param, std::vector<std::vector<int>> const& orbits, std::vector<int> const& where, std::vector<int> & basesize, std::vector<int> const& base, bool const& correct, std::vector<int> & missing, bool const& top)
{
  int nbLev=Stot.stabilizer.size();
  int n=Stot.n;
  StabLevel<Telt> & S = Stot.stabilizer[eLev];
  if (newgens.size() > 0) {
    auto GetFirst=[&]() -> int {
      for (auto & eBas : base) {
	for (auto & eGen : newgens)
	  if (PowAct(eBas, eGen) != eBas)
	    return eBas;
      }
      return -1;
    };
    int firstmove = GetFirst();
    if (eLev == nbLev-1) {
      S.orbit = {firstmove};
      S.transversal = std::vector<int>(n, -1);
      S.transversal[S.orbit[0]]  = GetLabelIndex_const(Stot.labels, Stot.identity);
      S.genlabels = {};
      S.treegen     = {};
      S.treegeninv               = {};
      StabLevel<Telt> TheLevel;
      TheLevel.diam=0;
      TheLevel.treedepth=0;
      Stot.stabilizer.push_back(TheLevel);
      if (!correct)
	basesize[where[S.orbit[0]]]++;
      missing = DifferenceVect( missing, {firstmove});
    }
    else {
      if (PositionVect(base,firstmove) < PositionVect(base,S.orbit[0])) {
	StabLevel<Telt> TheLevel;
	TheLevel.diam=0;
	TheLevel.treedepth=0;
	TheLevel.transversal = std::vector<int>(n, -1);
	TheLevel.transversal[firstmove]  = GetLabelIndex_const(Stot.labels, Stot.identity);
	TheLevel.orbit = {firstmove};
	TheLevel.genlabels = S.genlabels;
	//
	auto iter=Stot.stabilizer.begin();
	iter += eLev;
	Stot.stabilizer.insert(iter, TheLevel);
	if (!correct)
	  basesize[where[S.orbit[0]]]++;
	missing = DifferenceVect( missing, {firstmove} );
      }
    }
    if (!top || S.genlabels.size() == 0) {
      for (auto & eGen : newgens) {
	int posGen=GetLabelIndex(Stot.labels, eGen);
	S.genlabels.push_back(posGen);
      }
    }
    SCRSchTree(Stot, eLev, newgens);
    int nbGen=newgens.size();
    for (int iGen=0; iGen<nbGen; iGen++) {
      int jGen=nbGen-1-iGen;
      Telt g = SCRSift(Stot, eLev, newgens[jGen]);
      if (!g.isIdentity()) {
	SCRMakeStabStrong(Stot, eLev+1, {g}, param, orbits, where, basesize, base, correct, missing, false);
	S.diam=S.treedepth + Stot.stabilizer[eLev+1].diam;
	S.aux=Concatenation(S.treegen, S.treegeninv, Stot.stabilizer[eLev+1].aux);
      }
    }
  }
  std::vector<Telt> gen = Concatenation(S.treegen,S.treegeninv,std::vector<Telt>({Stot.identity}));
  std::vector<Telt> inv = Concatenation(S.treegeninv,S.treegen,std::vector<Telt>({Stot.identity}));
  int len = S.aux.size();
  Telt w = Stot.identity;
  if (len > 1) {
    Telt ran1 = SCRRandomPerm<Telt>(len);
    Face string = SCRRandomString(len);
    for (int x=0; x<len; x++) {
      int ximg=PowAct(x, ran1);
      if (string[x] == 1)
	w = w*S.aux[ximg];
    }
  }
  else {
    w=S.aux[1];
  }
  int mlimit=1;
  int m=0;
  while (m < mlimit) {
    m++;
    int ran = rand() % (S.orbit.size());
    std::vector<Telt> coset = CosetRepAsWord(Stot.labels, S.orbit[0], S.orbit[ran], S.transversal);
    coset = InverseAsWord(coset,gen,inv);
    if (!w.isIdentity()) {
      coset.push_back(w);
      std::pair<std::vector<Telt>, int> residue = SiftAsWord(Stot, eLev, coset);
      if (residue.second > 0) {
	Telt g = Product(residue.first);
	SCRMakeStabStrong(Stot, eLev+1, {g}, param, orbits, where, basesize, base, correct, missing, false);
	S.diam = S.treedepth + Stot.stabilizer[eLev+1].diam;
	S.aux = Concatenation(S.treegen, S.treegeninv, Stot.stabilizer[eLev+1].aux);
	m = 0;
      }
      else if (correct) {
	int l = 0;
	while (l < int(missing.size())) {
	  l++;
	  if (ImageInWord(missing[l],residue.first) != missing[l]) {
	    Telt g = Product(residue.first);
	    SCRMakeStabStrong(Stot, eLev+1, {g}, param, orbits, where, basesize, base, correct, missing, false);
	    S.diam = S.treedepth + Stot.stabilizer[eLev+1].diam;
	    S.aux = Concatenation(S.treegen, S.treegeninv, Stot.stabilizer[eLev+1].aux);
	    m = 0;
	    l = missing.size();
	  }
	}
      }
      else {
	int l=0;
	while (l < int(orbits.size())) {
	  l++;
	  if (int(orbits[l].size()) > param.param5) {
	    int j=0;
	    int jlimit=std::max(param.param6, basesize[l]);
	    while (j < jlimit) {
	      j++;
	      int ran=rand() % (orbits[l].size());
	      if (ImageInWord(orbits[l][ran],residue.first) != orbits[l][ran]) {
		Telt g = Product(residue.first);
		SCRMakeStabStrong(Stot, eLev+1, {g}, param, orbits, where, basesize, base, correct, missing, false);
		S.diam = S.treedepth + Stot.stabilizer[eLev+1].diam;
		S.aux = Concatenation(S.treegen,S.treegeninv,Stot.stabilizer[eLev+1].aux);
		m = 0;
		j = jlimit;
		l = orbits.size();
	      }
	    }
	  }
	  else {
	    int j = 0;
	    while (j < int(orbits[l].size())) {
	      if (ImageInWord(orbits[l][j],residue.first) != orbits[l][j]) {
		Telt g = Product(residue.first);
		SCRMakeStabStrong(Stot, eLev+1, {g}, param, orbits, where, basesize, base, correct, missing, false);
		S.diam = S.treedepth + Stot.stabilizer[eLev+1].diam;
		S.aux = Concatenation(S.treegen, S.treegeninv, Stot.stabilizer[eLev+1].aux);
		m = 0;
		j = orbits[l].size();
		l = orbits.size();
	      }
	      j++;
	    }
	  }
	}
      }
    }
  }
}

template<typename Telt>
std::vector<Telt> GetWpair(std::vector<Telt> const& Saux, int const& k, paramOpt const& param, Telt const& TheId)
{
  int len = Saux.size();
  if (len > 2*param.param3) {
    Telt ran1 = SCRRandomPerm<Telt>(len);
    Telt ran2 = SCRRandomPerm<Telt>(len);
    int len2 = rand() % (QuoInt(len,2));
    Face string = SCRRandomString(len+len2);
    Telt w1 = TheId;
    for (int x=0; x<len; x++) {
      int ximg=PowAct(x, ran1);
      if (string[x] == 1)
	w1 = w1 * Saux[ximg];
    }
    Telt w2 = TheId;
    for (int x=0; x<len2; x++) {
      int ximg=PowAct(x, ran2);
      if (string[x+len] == 1)
	w2 = w2 * Saux[ximg];
    }
    return {w1, w2};
  }
  else {
    if (len >= 2*k)
      return {Saux[2*k-2], Saux[2*k-1]};
    else
      return {Saux[2*k-2], TheId};
  }
  
}


template<typename Telt>
Telt SCRStrongGenTest(StabChain<Telt> const& Stot, int const& eLev, paramOpt const& param, std::vector<std::vector<int>> const& orbits, std::vector<int> const& basesize, std::vector<int> const& base, bool const& correct, std::vector<int> const& missing)
{
  StabLevel<Telt> const& S = Stot.stabilizer[eLev];
  int k = 0;
  int l, j;
  while (k < param.param1) {
    k++;
    int len=S.aux.size();
    std::vector<Telt> w = GetWpair(S.aux, k, param, Stot.identity);
    int m = 0;
    int mlimit = param.param2*S.diam;
    while (m < mlimit) {
      m++;
      std::vector<Telt> ranword = RandomElmAsWord(Stot, eLev);
      int i = 0;
      while (i < 2) {
        i++;
	if (!w[i-1].isIdentity()) {
	  ranword.push_back(w[i-1]);
	  std::pair<std::vector<Telt>, int> residue = SiftAsWord(Stot, eLev, ranword);
	  if (residue.second > 0) {
	    return Product(residue.first);
	  }
	  else {
	    if (correct) {
	      l=0;
	      while (l < int(missing.size())) {
		if (ImageInWord(missing[l],residue.first) != missing[l])
		  return Product(residue.first);
		l++;
	      }
	    }
	    else {
	      l=0;
	      while (l < int(orbits.size())) {
	        l++;
		if (int(orbits[l].size()) > param.param5) {
		  j = 0;
		  int jlimit = std::max(param.param6, basesize[l]);
		  while (j < jlimit) {
		    j++;
		    int ran = rand() % (orbits[l].size());
		    if (ImageInWord(orbits[l][ran],residue.first) != orbits[l][ran])
		      return Product(residue.first);
		  }
		}
		else {
		  j = 0;
		  while (j < int(orbits[l].size())) {
		    j++;
		    if (ImageInWord(orbits[l][j],residue.first) != orbits[l][j])
		      return Product(residue.first);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    if (len <= 2*k)
      k = param.param1;
  }
  return Stot.identity;
}



template<typename Telt>
Telt SCRStrongGenTest2(StabChain<Telt> const& Stot, int const& eLev, paramOpt const& param)
{
  int k = 0;
  StabLevel<Telt> const& S = Stot.stabilizer[eLev];
  while (k < param.param3) {
    k++;
    int len = S.aux.size();
    std::vector<Telt> w = GetWpair(S.aux, k, param, Stot.identity);
    int m = 0;
    int mlimit = param.param4 * S.diam;
    while (m < mlimit) {
      m++;
      Telt ranelement = Stot.identity;
      int uLev=eLev;
      while (Stot.stabilizer[uLev].genlabels.size() > 0) {
        int p = Random(Stot.stabilizer[uLev].orbit);
	while (p != Stot.stabilizer[uLev].orbit[0]) {
	  int pos= Stot.stabilizer[uLev].transversal[p];
	  Telt eElt = Stot.labels[pos];
	  ranelement = LeftQuotient(eElt, ranelement );
	  p = PowAct(p, eElt);
	}
	uLev++;
      }
      int i = 0;
      while (i < 2) {
	i++;
	if (!w[i-1].isIdentity()) {
	  ranelement = ranelement*w[i-1];
	  Telt residue = SCRSift(Stot, eLev, ranelement);
	  if (!residue.isIdentity())
	    return residue;
	}
      }
    }
    if (len <= 2*k)
      k = param.param3;
  }
  return Stot.identity;
}




template<typename Telt>
Telt VerifyStabilizer(StabChain<Telt> const& Stot, int const& eLev, Telt const& z, std::vector<int> const& missing, bool const& correct)
{
  StabLevel<Telt> const& S = Stot.stabilizer[eLev];
  int pt1 = S.orbit[0];
  Telt zinv = ~z;
  int pt2 = PowAct(pt1, zinv);
  Telt result = Stot.identity;
  StabLevel<Telt> const& stab = Stot.stabilizer[eLev+1];
  bool flag;
  if (PowAct(pt2, zinv) == pt1) {
    flag= true;
    result = SCRSift(Stot, eLev+1, z*z);
  }
  else {
    flag=false;
  }
  int n=Stot.n;
  std::vector<int> where1(n, -1);
  std::vector<int> leaders{pt2};
  int orbit1count = 1;
  std::vector<int> transversal(n,-1);
  transversal[pt2] = GetLabelIndex_const(Stot.labels, Stot.identity);
  where1[pt2] = 1;
  std::vector<int> orb{pt2};
  int j = 1;
  while (j <= int(orb.size()) ) {
    for (auto & iGen : stab.genlabels) {
      int k = SlashAct(orb[j], Stot.labels[iGen]);
      if (transversal[k] == -1) {
	transversal[k] = iGen;
	orb.push_back(k);
	where1[k] = orbit1count;
      }
    }
    j++;
  }
  for (auto &i : Stot.stabilizer[eLev].orbit) {
    if (where1[i] == -1) {
      orbit1count++;
      leaders.push_back(i);
      orb = {i}; 
      where1[i] = orbit1count;
      transversal[i] = GetLabelIndex_const(Stot.labels, Stot.identity);
      int j = 1;
      while (j <= int(orb.size()) ) {
	for (auto & iGen : stab.genlabels) {
	  int k = SlashAct(orb[j], Stot.labels[iGen]);
	  if (transversal[k] == -1) {
	    transversal[k] = iGen;
	    orb.push_back(k);
	    where1[k] = orbit1count;
	  }
	}
        j++;
      }
    }
  }
  StabChain<Telt> chain = RestrictedStabChain(Stot, eLev+1);
  int nbLeader=leaders.size();
  int missSiz=missing.size();
  for (int j=0; j<nbLeader; j++) {
    if (result.isIdentity()) {
      int i = leaders[nbLeader - 1 - j];
      ChangeStabChain(chain, 0, {i}, false);
      Telt w1, w1inv;
      if (i == pt2) {
        w1 = z;
	w1inv = zinv;
      }
      else {
	std::vector<Telt> w1_w = CosetRepAsWord(Stot.labels, pt1, i, Stot.stabilizer[eLev].transversal);
        w1 = Product(w1_w);
        w1inv = Inverse(w1);
      }
      for (auto & iGen : chain.stabilizer[1].genlabels) {
	Telt g = chain.labels[iGen];
	if (result.isIdentity()) {
	  if (correct) {
	    // There is a possible source of problems here 
	    std::pair<std::vector<Telt>, int> residue = SiftAsWord(Stot, eLev+1, {w1inv,g,w1});
	    if (residue.second != 0) {
	      result = Product(residue.first);
	    }
	    else {
	      int l = 0;
	      while (l < missSiz && result.isIdentity()) {
		if (ImageInWord(missing[l],residue.first) != missing[l])
		  result = Product(residue.first);
		l++;
	      }
	    }
	  }
	  else {
	    result = SCRSift(Stot, eLev+1, w1inv*g*w1);
	  }
	}
      }
    }
  }
  if (result.isIdentity()) {
    Face where2(where1.size());
    for (auto &i : Stot.stabilizer[eLev].orbit) {
      if (result.isIdentity() && !where2[i]) {
        orb = {i};
	where2[i]=1;
	for (auto & pnt : orb) {
	  for (auto & iGen : chain.stabilizer[1].genlabels) {
	    Telt gen = chain.labels[iGen];
	    int img = PowAct(pnt, gen);
	    if (!where2[img]) {
	      orb.push_back(img);
	      where2[img]=1;
	    }
	  }
	}
	int ipZ=PowAct(i, z);
	if (flag && !where2[ipZ]) {
	  orb = {ipZ};
	  where2[ipZ]=1;
	  for (auto & pnt : orb) {
	    for (auto & iGen : chain.stabilizer[1].genlabels) {
	      Telt gen = chain.labels[iGen];
	      int img = PowAct(pnt, gen);
	      if (!where2[img]) {
		orb.push_back(img);
		where2[img]=1;
	      }
	    }
	  }
	}
	std::vector<Telt> w1 = CosetRepAsWord(Stot.labels, pt1, leaders[where1[i]], Stot.stabilizer[eLev].transversal);
	std::vector<Telt> w2 = CosetRepAsWord(Stot.labels, leaders[where1[i]], i, transversal); 
	std::vector<Telt> w3 = CosetRepAsWord(Stot.labels, leaders[where1[ipZ]], ipZ, transversal);
	std::vector<Telt> w4;
	if (where1[i] != where1[ipZ]) {
	  w4 = CosetRepAsWord(Stot.labels, pt1, leaders[where1[ipZ]], S.transversal);
	}
	else {
	  w4 = w1;
	}
        Telt schgen = Inverse(Product(w1)) * Inverse(Product(w2));
	if (correct) {
	  std::vector<Telt> schgen_w = Concatenation(std::vector<Telt>({schgen,z}),w3,w4);
	  int eLevLoc=0; // almost certainly wrong
	  std::pair<std::vector<Telt>,int> residue = SiftAsWord(Stot, eLevLoc, schgen_w);
	  if (residue.second != 0) {
	    result = Product(residue.first);
	  }
	  else {
	    int l = 0;
	    while (l < missSiz && result.isIdentity()) {
	      if (ImageInWord(missing[l],residue.first) != missing[l]) {
	        result = Product(residue.first);
	      }
	      l++;
	    }
	  }
	}
	else {
	  schgen = schgen*z*Product(w3)*Product(w4);
	  result = SCRSift(Stot, eLev+1, schgen);
	}
      }
    }
  }
  return result;
}



  

template<typename Telt, typename Tint>
StabChain<Telt> StabChainRandomPermGroup(std::vector<Telt> const& gens, Telt const& id, StabChainOptions<Tint> const& options)
{
  int k;
  if (options.random == 1000) {
    k = 1;
  }
  else {
    double eFrac=1 - double(options.random) / double(1000);
    double Expo=double(3) / double(5);
    k=0;
    double ThePow=1;
    while(true) {
      if (ThePow < eFrac)
	break;
      k++;
      ThePow *= Expo;
    }
  }
  //
  int degree = LargestMovedPoint( gens);
  paramOpt param;
  std::vector<int> UsedKnownBase;
  if (options.knownBase.size() > 0 && int(options.knownBase.size()) < 4 + LogInt(degree,10)) {
    param = {k, 4, 0, 0, 0, 0};
    UsedKnownBase = options.knownBase;
  }
  else {
    param = {QuoInt(k,2), 4, QuoInt(k+1,2), 4, 50, 5};
  }
  if (options.random <= 200) {
    param.param2 = 2;
    param.param4 = 2;
  }
  std::vector<int> givenbase = options.base;
  //
  Tint order, limit;
  int warning;
  if (options.size == 0) {
    order = options.size;
    warning = 0;
    limit = 0;
  }
  else {
    order = 0;
    limit = options.limit;
  }
  bool correct=false;
  if (UsedKnownBase.size() > 0)
    correct = true;
  //
  std::vector<int> missing;
  std::vector<int> base;
  std::vector<int> basesize;
  std::vector<int> where;
  std::vector<std::vector<int>> orbits;
  if (correct) {
    base=Concatenation(givenbase,DifferenceVect(UsedKnownBase,givenbase));
    missing = VectorAsSet(UsedKnownBase);
  }
  else {
    std::vector<int> TheBase(degree);
    for (int u=0; u<degree; u++)
      TheBase[u]=u;
    base=Concatenation(givenbase, DifferenceVect(TheBase, givenbase));
    std::vector<std::vector<int>> orbits2 = OrbitsPerms(gens,TheBase);
    for (auto & eOrb : orbits2)
      if (eOrb.size() > 1)
	orbits.push_back(eOrb);
    basesize = std::vector<int>(orbits.size(), 0);
    where = std::vector<int>(degree, -1);
    for (int iOrb=0; iOrb<int(orbits.size()); iOrb++)
      for (auto & ePt : orbits[iOrb])
	where[ePt] = iOrb;
    if (int(orbits.size()) * 40 > degree) {
      param.param1 = 0;
      param.param3 = k;
    }
  }

  bool ready=false;
  int n=degree;
  StabChain<Telt> Stot = EmptyStabChain<Telt>(n);
  std::vector<Telt> eNew=gens;
  Telt result;
  while(!ready) {
    int eLevLoc=0; // almost certainly wrong
    SCRMakeStabStrong(Stot, eLevLoc, eNew, param, orbits, where, basesize, base, correct, missing, true);
    Tint size = 1;
    for (auto const& eStab : Stot.stabilizer)
      size *= eStab.orbit.size();
    if (size == order || size == limit) {
      ready = true;
    }
    else {
      if (size < order) {
        result = Telt(n);
	if (options.random == 1000) {
	  Telt result;
	  int eLev=0; // certainly wrong
	  if (correct) {
	    paramOpt param{1, 10 / Stot.stabilizer[eLev].diam, 0, 0, 0, 0};
	    result = SCRStrongGenTest(Stot, eLev, param, orbits, basesize, base, correct, missing);
	  }
	  else {
	    paramOpt param{0, 0, 1, 10 / Stot.stabilizer[eLev].diam, 0, 0};
	    result = SCRStrongGenTest2(Stot, eLev, param);
	  }
	  std::pair<bool,Telt> resultComp = {true, result};
	  if (result.isIdentity()) {
	    resultComp = VerifySGS(Stot,missing,correct);
	    result=resultComp.second;
	  }
	  if (resultComp.first && result.isIdentity()) {
	    // There are mysterious things going on here
	    std::cerr << "Warning, computed and given size differ\n";
	    ready = true;
	  }
	  else {
	    if (!resultComp.first) {
	      while(true) {
		int eLev=0; // certainly wrong
		paramOpt param{0, 0, 1, 10 / Stot.stabilizer[eLev].diam, 0, 0};
		result = SCRStrongGenTest2(Stot, eLev, param);
		if (!result.isIdentity())
		  break;
	      }
	    }
	    eNew = {result};
	  }
	}
	else {
	  warning = 0;
	  if (correct) {
	    while (result.isIdentity()) {
	      warning++;
	      if (warning > 5)
		std::cerr << "Warning, computed and given size differ\n";
	      int eLev=0; // certainly wrong
	      result = SCRStrongGenTest(Stot, eLev, param, orbits, basesize, base, correct, missing);
	    }
	  }
	  else {
	    while (result.isIdentity()) {
	      warning++;
	      if (warning > 5)
		std::cerr << "Warning, computed and given size differ\n";
	      int eLev=0; // certainly wrong
	      result=SCRStrongGenTest(Stot, eLev, param, orbits, basesize, base, correct, missing);
	      if (result.isIdentity()) 
		result=SCRStrongGenTest2(Stot, eLev, param);
	    }
	  }
	  eNew = {result};
	}
      }
      else {
	if (options.random == 1000) {
	  int eLev=0; // certainly wrong
	  if (correct) {
	    paramOpt param{1,10/Stot.stabilizer[eLev].diam,0,0,0,0};
	    result = SCRStrongGenTest(Stot, eLev, param, orbits, basesize, base, correct, missing);
	  }
	  else {
	    paramOpt param{0,0,1,10/Stot.stabilizer[eLev].diam,0,0};
	    result = SCRStrongGenTest2(Stot, eLev, param);
	  }
	  std::pair<bool,Telt> resultComp = {true, result};
	  if (resultComp.first && result.isIdentity()) {
	    resultComp = VerifySGS(Stot,missing,correct);
	    result = resultComp.second;
	  }
	  if (resultComp.first && result.isIdentity()) {
	    ready = true;
	  }
	  else {
	    if (!resultComp.first) {
	      while(true) {
		int eLev=0; // certainly wrong
		paramOpt param{0,0,1,10/Stot.stabilizer[eLev].diam,0,0};
		result = SCRStrongGenTest2(Stot, eLev, param);
		if (!result.isIdentity())
		  break;
	      }
	    }
	    eNew = {result};
	  }
	}
	else {
	  int eLev=0; // almost certainly wrong
	  result =SCRStrongGenTest(Stot, eLev, param, orbits, basesize, base, correct, missing);
	  if (!result.isIdentity()) {
	    eNew = {result};
	  }
	  else {
	    if (correct) {
	      ready = true;
	    }
	    else {
	      int eLev=0; // certainly wrong
	      result=SCRStrongGenTest2(Stot, eLev, param);
	      if (result.isIdentity()) {
		ready =true;
	      }
	      else {
		eNew = {result};
	      }
	    }
	  }
	}
      }
    }
  }
  return Stot;
}



// Is the Stot really a const function? Not sure really.
template<typename Telt>
std::pair<bool,Telt> VerifySGS(StabChain<Telt> const& Stot, std::vector<int> const& missing, bool const& correct) 
{
  int n = Stot.n;
  int len = Stot.stabilizer.size();
  // list := ListStabChain(S);
  std::pair<bool, Telt> result = {true, Stot.identity};
  //
  std::vector<int> TotSet(n);
  for (int i=0; i<n; i++)
    TotSet[i] = i;
  //
  int i = 0;
  while (i < len && result.second.isIdentity()) {
    StabChain<Telt> temp = RestrictedStabChain(Stot, len-1-i);
    InsertTrivialStabilizer(temp, len-i, Stot.stabilizer[len -i].orbit[0]);
    int gencount = 0;
    while (gencount < int(Stot.stabilizer[len - i].genlabels.size()) && result.second.isIdentity()) {
      gencount++;
      int posgen=Stot.stabilizer[len - i].genlabels[gencount];
      Telt gen = Stot.labels[posgen];
      std::vector<int> set = VectorAsSet(temp.stabilizer[0].orbit);
      if (set == OnSets(set,gen)) {
	if (correct) {
	  int eLevLoc=0; // almost certainly wrong
	  std::pair<std::vector<Telt>, int> residue = SiftAsWord(temp, eLevLoc, {gen});
	  if (residue.second != 0) {
	    result.second = Product(residue.first);
	  }
	  else  {
	    int l = 0;
	    while ( l < int(missing.size()) && result.second.isIdentity()) {
	      if (ImageInWord(missing[l],residue.first) != missing[l]) {
		result.second = Product(residue.first);
	      }
	      l++;
	    }
	  }
	}
	else {
	  result.second = SCRSift(temp, 0, gen);
	}
      }
      else {
	StabChain<Telt> temp2;
	Telt newgen;
	if (set.size() == 1) {
	  temp2 = temp; 
	  newgen = gen;
	}
	else {
	  std::vector<Telt> longer;
	  for (auto & eIdx : temp.stabilizer[0].genlabels)
	    longer.push_back(temp.labels[eIdx]);
	  longer.push_back(gen);
	  std::vector<int> orbit = OrbitPerms(longer, temp.stabilizer[0].orbit[0]);
	  std::vector<std::vector<int>> blks = Blocks_from_seed(longer, orbit, set);
	  if (blks.size() * blks.size() != orbit.size()) {
	    result = {false, Telt(n)};
	  }
	  else {
	    int pos = PositionVect(blks, set);
	    std::function<Telt(Telt const&)> f = MapElementToSetPlusBlocks<Telt>(blks, n);
	    temp2 = HomomorphismMapping(temp, f);
	    newgen = f(gen);
	    InsertTrivialStabilizer(temp2, len-i, n + pos);
	  }
	}
	if (i == 0 && set.size() == 1) {
	  int eLen=CycleLength(newgen, temp2.stabilizer[0].orbit[0]);
	  result.second = PowerGroupElement(newgen, eLen);
	  int eLevLoc=0;
	  AddGeneratorsExtendSchreierTree(temp2, eLevLoc, {newgen});
	}
	else {
	  if (result.second.isIdentity()) {
	    int eLevLoc=0;
	    AddGeneratorsExtendSchreierTree(temp2, eLevLoc, {newgen} );
	    std::vector<std::vector<int>> blks = Blocks_without_seed(temp2.labels, temp2.stabilizer[0].orbit);
	    if (blks.size() > 1) {
	      int leader = temp2.stabilizer[0].orbit[0];
	      std::vector<int> block = GetBlock(blks, leader);
	      int point = GetNonTrivialPointInBlock(block, leader);
	      newgen = Product(CosetRepAsWord(Stot.labels, leader, point, temp2.stabilizer[0].transversal)); 
	      temp2 = RestrictedStabChain(temp2, 1);
	      InsertTrivialStabilizer(temp2, len-i, leader);
	      AddGeneratorsExtendSchreierTree(temp2, eLevLoc, {newgen});
	      result.second = VerifyStabilizer(temp2, eLevLoc, newgen, missing, correct);
	      if (leader > n) {
		AddGeneratorsExtendSchreierTree(temp, eLevLoc, {RestrictedPermNC(newgen, TotSet)});
	      }
	      else {
		temp = temp2;
	      }
	      gencount--;
	    }
	    else {
	      result.second = VerifyStabilizer(temp2, eLevLoc, newgen, missing, correct); 
	      if (temp2.stabilizer[0].orbit[0] >= n) {
		AddGeneratorsExtendSchreierTree(temp, eLevLoc, {gen});
	      }
	    }
	    int eLev=0; // almost certainly wrong
	    if (!result.second.isIdentity() && temp2.stabilizer[eLev].orbit[0] >= n) {
	      result.second = RestrictedPermNC(result.second, TotSet);
	    }
	  }
	}
      }
    }
    i++;
  }
  return result;
}


} 

#endif

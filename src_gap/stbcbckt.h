#ifndef DEFINE_STBCBCKT_H
#define DEFINE_STBCBCKT_H

#include "StabChainMain.h"
/*
#############################################################################
##
#W  stbcbckt.gi                 GAP library                    Heiko Theißen
##
##
#Y  Copyright (C)  1997,  Lehrstuhl D für Mathematik,  RWTH Aachen, Germany
#Y  (C) 1998 School Math and Comp. Sci., University of St Andrews, Scotland
#Y  Copyright (C) 2002 The GAP Group
##
*/

namespace gap {


template<typename Telt>
struct permPlusBool {
  int status; // values in {int_false, int_true, int_perm}
  Telt val;
};


template<typename Telt>
permPlusBool<Telt> ExtendedT(Telt const& t, int const& pnt, int& img, int const& simg, StabChain<Telt> const& Stot, int const& eLev)
{
  if (simg == -1)
    img = SlashAct(img, t);
  else
    img = simg;
    
  // If <G> fixes <pnt>, nothing more can  be changed, so test whether <pnt>
  // = <img>.
  int bpt = BasePoint(Stot, eLev);
  if (bpt != pnt) {
    if (pnt != img)
      return {int_false,{}};
  }
  else {
    if (Stot.stabilizer[eLev].translabels[img] == -1) 
      return {int_false,{}};
    else
      t = LeftQuotient(InverseRepresentative(Stot, eLev, img ), t);
  }
  return {int_true, t};
}

template<typename Telt>
struct StabChainPlusLev {
  int status; // possible values in {int_false, int_true, int_int, int_stablev} 
  int value_int;
  StabChain<Telt> Stot;
  int eLev;
};


struct singStrat {
  int p;
  int s;
  int i;
};


struct Refinement {
public:
  Refinement(int const& val1, int const& val2) {
    nature = 0;
    inputProcessfix = {val1, val2};
  }
  Refinement(Partition const& ePart, std::vector<singStrat> const& strat) {
    nature = 1;
    inputIntersection = {ePart, strat};
  }
  int nature; // 0 for PROCESSFIX, 1 for INTERSECTION
  std::pair<int,int> inputProcessfix;
  std::pair<Partition,std::vector<singStrat>> inputIntersection;
};


// The underscore nature of a function can be seen in stbcbckt top.
// Since we did not implement all the algorithms of stbcbckt, the value is always 0.
int UnderscoreNature(int const& nature)
{
  return 0;
}



template<typename Telt>
struct dataType {
  Partition P;
};


template<typename Telt>
struct rbaseType {
  std::vector<int> domain;
  std::vector<int> base;
  std::vector<int> where;
  //
  dataType<Telt> data;
  std::vector<std::vector<int>> fix;
  //
  std::vector<Refinement> rfm;
  Partition partition;
  std::vector<StabChainPlusLev<Telt>> lev;
  StabChainPlusLev<Telt> level;
  //
  std::vector<StabChainPlusLev<Telt>> lev2;
  StabChainPlusLev<Telt> level2;
};


bool ProcessFixpoint_rbase(rbaseType<Telt> & rbase, int const& pnt)
{
  if (rbase.level2.status != 0 && rbase.level2.status != 1) {
    ChangeStabChain(rbase.level2, {pnt});
    if (BasePoint(rbase.level2) == pnt) 
      rbase.level2.eLev++;
  }
  if (rbase.level_status == 2) {
    rbase.level.value_int--;
  }
  else {
    ChangeStabChain(rbase.level, {pnt} );
    if (BasePoint(rbase.level) == pnt) {
      rbase.level.eLev++;
    }
    else {
      return false;
    }
  }
  return true;
}



template<typename Telt>
struct imageType {
  int depth;
  Partition partition;
  permPlusBool<Telt> perm;
  StabChainPlusLev<Telt> level;
  std::vector<int> bimg;
  //
  permPlusBool<Telt> perm2;
  StabChainPlusLev<Telt> level2;
};


template<typename Telt>
bool ProcessFixpoint_image(imageType<Telt> & image, int const& pnt, int const& img, int const& simg)
{
  if (image.perm.status != 1) {
    permPlusBool<Telt> t = ExtendedT(image.perm, pnt, img, simg, image.level);
    if (!t.status)
      return false;
    else {
      if (BasePoint(image.level ) == pnt)
        image.level.eLev++;
    }
    image.perm = t;
  }
  if (image.level2.status != 0) {
    permPlusBool<Telt> t = ExtendedT(image.perm2, pnt, img, 0, image.level2);
    if (!t.status)
      return false;
    else {
      if (BasePoint(image.level2 ) == pnt)
        image.level2.eLev++;
    }
    image.perm2 = t;
  }
  return true;
}


template<typename Telt>
bool IsTrivialRBase(rbaseType<Telt> const& rbase)
{
  if (rbase.level.status == 2) {
    if (rbase.level.value_int <= 1)
      return true;
  }
  if (rbase.level.status == 3) {
    int eLev=base.level.eLev;
    if (rbase.level.stabilizer[eLev].genlabels.size() == 0)
      return true;
  }
  return false;
}



template<typename Telt>
rbaseType<Telt> EmptyRBase(std::vector<StabChain<Telt>> const& G, bool const& IsId, std::vector<int> const& Omega, Partition const& P)
{
  rbase<Telt> rbase;
  rbase.domain = Omega;
  rbase.base = {};
  rbase.where = {};
  rbase.rfm = {};
  rbase.partition = P;
  rbase.lev = {};
  if (F.size() == 2) {
    if (IsId) {
      rbase.NeedLevel2=false;
      rbase.SetLevelStabChain2=false;
    }
    else {
      rbase.SetLevelStabChain2=true;
      rbase.level2 = {G[1], 0};
      rbase.lev2 = {};
    }
  }
  else {
    rbase.NeedLevel2=false;
    rbase.SetLevelStabChain2=false;
  }
  rbase.level = {F[0], 0};
  for (auto & pnt : Fixcells(P))
    ProcessFixpoint_rbase(rbase, pnt);
  return rbase;
}




template<typename Telt>
bool MeetPartitionStrat(imageType<Telt> const& image, Partition const& S, Telt const& g, std::vector<singStrat> const& strat)
{
  if (strat.size() == 0)
    return false;
  for (auto & pRec : strat) {
    if ((pRec.p == -1 && !ProcessFixpoint_image(image, p.s, FixpointCellNo(P, p.i), -1)) ||
	(pRec.p != -1 && SplitCell(image.partition, pRec.p, S, pRec.s, g, pRec.i ) != pRec.i))
      return false;
  }
  return true;
}

/*
#############################################################################
##
#F  StratMeetPartition( <rbase>, <P>, <S>, <g> )  . construct a meet strategy
##
##  Entries in <strat> have the following meaning:
##    [p,s,i] (p<>0) means that `0 < |P[p]\cap S[s]/g| = i < |P[p]|',
##            i.e., a new cell with <i> points was appended to <P>
##                  (and these <i> have been taken out of `P[p]'),
##    [0,a,p] means that fixpoint <a> was mapped to fixpoint in `P[p]',
##            i.e., `P[p]' has become a one-point cell.
##
*/
template<typenam Telt>
std::vector<singStrat> StratMeetPartition(rbaseType<Telt> const& rbase, Partition const& P, Partition const& S, Telt const& g)
{
  std::vector<singStrat> strat;
  std::vector<int> cellsP;
  if (g.isIdentity()) {
    cellsP = P.cellno;
  }
  else {
    cellsP = std::vector<int>(P.cellno.size(), 0);
    for (int i=0; i<NumberCells(P); i++) {
      std::vector<int> cell = Cell( P, i );
      for (auto & eVal : cell) {
        int img=PowAct(eVal, g);
	cellsP[img] = i;
      }
    }
  }
  // If <S> is just a set, it is interpreted as partition ( <S>|<S>^compl ).
  int nrcells = NumberCells(S) - 1;

  for (int s=0; s<nrcells; s++) {
    // now split with cell number s of S.
    std::vector<int> p=Cell(S, s);
    std::vector<int> p2;
    for (auto & eVal : p)
      p2.push_back(cellsP[eVal]);
    CollectedResult<int> p3=Collected(p);
    std::vector<int> splits;
    for (int h=0; h<int(p3.LVal.size()); h++) {
      // a cell will split iff it contains more points than are in the s-cell
      if (P.lengths[p3.LVal[h]] > p3.LMult[h])
        splits.push_back(p3.LVal[h]);
    }
    for (auto & pVal : splits) {
      // Last argument true means that the cell will split.
      int i = SplitCell(P, pVal, S, s, g, true);
      if (!g.isIdentity()) {
	std::vector<int> cell = Cell(P, NumberCells(P));
	for (auto & eVal : cell) {
	  int img=PowAct(eVal, g);
	  cellsP[img] = NumberCells(P);
	}
      }
      strat.push_back({p, s, i});
      // If  we have one  or two  new fixpoints, put  them  into the base.
      if (i == 0) {
        int pnt = FixpointCellNo(P, NumberCells(P));
	ProcessFixpoint_rbase(rbase, pnt);
	strat.push_back({-1, pnt, NumberCells(P)});
	if (IsTrivialRBase(rbase)) 
	  return strat;
      }
      if (P.lengths[p] == 1) {
        int pnt = FixpointCellNo(P, p);
	ProcessFixpoint_rbase(rbase, pnt);
	strat.push_back({-1, pnt, p});
	if (IsTrivialRBase(rbase))
	  return strat;
      }
    }
  }
  return strat;
}


template<typename Telt>
void RegisterRBasePoint(Partition const& P, rbaseType<Telt> & rbase, int const& pnt)
{
  if (rbase.level2.status != 0 && rbase.level2.status != 1)
    rbase.lev2.push_back(rbase.level2);
  rbase.lev.push_back(rbase.level);
  rbase.base.push_back(pnt);
  int k = IsolatePoint(P, pnt);
  rbase.where.push_back(k);
  int len=rbase.rfm.size();
  rbase.rfm.push_back({});
  if (P.lengths[k] == 1) {
    int pnt = FixpointCellNo(P, k);
    ProcessFixpoint_rbase(rbase, pnt);
    rbase.rfm[len].push_back(Refinement({pnt,k}));
  }
  if (rbase.level2.status != 0) {
    auto MainInsert=[&](StabChainPlusLev const& lev) -> void {
      if (lev.status != 2) {
	Partition O = OrbitsPartition(lev, rbase.domain);
	std::vector<singStrat> strat = StratMeetPartition(rbase, P, O);
	rbase.rfm[len].push_back(Refinement({O,strat}));
      }
    };
    if (rbase.level2.status == 1)
      MainInsert(rbase.level);
    else
      MainInsert(rbase.level2);
  }
}




template<typename Telt>
void NextRBasePoint(Partition & P, rbaseType<Telt> & rbase)
{
  std::vector<int> lens = P.lengths;
  std::vector<int> order = ClosedInterval(0, NumberCells(P));
  SortParallel(lens, order);
  int k = PositionProperty(lens, [](int const& x) -> int {return x != 1;});
  int l = -1;
  while(l == -1) {
    if (rbase.level.status == 2) {
      l = 0;
    }
    else {
      l = PositionProperty(ClosedInterval(0, lens[k]), [](int const& i) -> int {return !IsFixedStabilizer(rbase.level, P.points[i+P.firsts[order[k]]]);});
    }
    k++;
  }
  int p = P.points[ P.firsts[ order[k] ] + l ];
  RegisterRBasePoint(P, rbase, p);
}

template<typename Telt>
bool Refinements_ProcessFixpoint(rbaseType<Telt> & rbase, imageType<Telt> & image, int const& pnt, int const& cellnum)
{
 int img = FixpointCellNo(image.partition, cellnum);
 return ProcessFixpoint_image(image, pnt, img, -1);
}


template<typename Telt>
bool Refinements_Intersection(rbaseType<Telt> & rbase, imageType<Telt> & image, Partition const& Q, std::vector<singStrat> const& strat)
{
  Telt t;
  if (image.level2.status == 0) {
    t = image.perm.val;
  }
  else {
    t = image.perm2.val;
  }
  Telt tinv =Inverse(t);
  return MeetPartitionStrat(rbase, image, Q, t, strat);
}


template<typename Telt>
int RRefine(rbaseType<Telt> & rbase, imageType<Telt> & image, bool const& uscore)
{
  auto BoolToInt=[&](bool const& val) -> int {
    if (val)
      return int_true;
    return int_false;
  };
  auto Evaluation=[&](Refinement const& eRef) -> bool {
    if (eRef.nature == int_false)
      return Refinements_ProcessFixpoint(rbase, image, eRef.inputProcessfix.first, eRef.inputProcessfix.second);
    if (eRef.nature == int_true)
      return Refinements_Intersection(rbase, image, eRef.inputIntersection.first, eRef.inputIntersection.second);
    return true;
  };
  if (!uscore) {
    for (auto & Rf : rbase.rfm[image.depth]) {
      bool t = Evaluation(Rf);
      if (!t) {
	return int_fail;
      }
      else {
	if (!t) {
	  return BoolToInt(t);
	}
      }
    }
    return int_true;
  }
  else {
    for (auto & Rf : rbase.rfm[image.depth]) {
      if (UnderscoreNature(Rf.nature) {
	bool t = Evaluation(Rf);
	if (!t) {
	  return int_fail;
	}
	else {
	  if (!t) {
	    return BoolToInt(t);
	  }
	}
      }
    }
    return int_true;
  }
}		       }

		       
template<typename Telt>
bool PBIsMinimal(std::vector<int> const& range, int const& a, int const& b, StabChain<Telt> const& Stot, int const& eLev)
{
  if (IsInBasicOrbit(Stot, eLev, b)) {
    for (auto & pVal : Stot.stabilizer[eLev].orbit) {
      if (a > pVal)
        return false;
    }
    return true;
  }
  if (b < a)
    return false;
  if (IsFixedStabilizer(Stot, eLev, b))
    return true;

  std::vector<int> orb{b};
  int pos=0;
  Face old = BlistList(range, orb);
  while(true) {
    int siz=orb.size();
    if (pos == siz)
      break;
    for (int i=pos; i<siz; i++) {
      int pnt=orb[i];
      for (auto & lVal : Stot.stabilizer[eLev].genlabels) {
        int img = PowAct(pnt, Stot.labels[lVal]);
        if (!old[img]) {
          if (img < a)
            return false;
          old[img]=true;
          orb.push_back(img);
        }
      }
    }
    pos=siz;
  }
  return true;
}

template<typename Telt>
void SubtractBlistOrbitStabChain(Face & blist, std::vector<Telt> const& LGen, int const& pnt)
{
  std::vector<int> orb{pnt};
  std::vector<int> orb = {pnt};
  blist[pnt]=false;
  int pos=0;
  int PrevPos=0;
  while(true) {
    int siz = orb.size();
    if (pos == siz) {
      break;
    }
    for (int ePos=0; ePos<siz; ePos++) {
      pnt=orb[ePos];
      for (auto& eGen : LGen) {
        int img = PowAct(pnt, gen);
        if (blist[img]) {
          blist[img]=false;
          orb.push_back(img);
        }
      }
    }
    pos=siz;
  }
}


template<typename Telt>
struct ResultPBT {
  int nature; // Allowed values in {int_group, int_fail, int_perm}.
              // int_group for group
              // int_fail for fail
              // int_perm for equivalence element
  StabChain<Telt> stab;
  Telt res;
};


template<typename Telt>
Telt MappingPermListList(n, std::vector<int> const& src, int const& dst)
{
  std::vector<int> ListImage(n);
  Face StatusSrc(n);
  Face StatusDst(n);
  for (int i=0; i<n; i++) {
    StatusSrc[i]=1;
    StatusDst[i]=1;
  }
  int len = src.size();
  for (int j=0; j<len; j++) {
    ListImage[src[i]] = dst[i];
    StatusSrc[src[i]] = 0;
    StatusDst[dst[i]] = 0;
  }
  int sizRemain = n - len;
  dynamic_bitset::size_type posSrc=StatusSrc.find_first();
  dynamic_bitset::size_type posDst=StatusDst.find_first();
  for (int u=0; i<sizRemain; i++) {
    ListImage[posSrc] = posDst;
    posSrc=StatusSrc.find_next(posSrc);
    posDst=StatusDst.find_next(posDst);
  }
  return Telt(ListImage);
}




template<typename Telt>
ResultPBT<Telt> PartitionBacktrack(StabChain<Telt> onst& G, std::function<bool(Telt const&)> const& Pr, bool const& repr, rbaseType<Telt> & rbase, dataType<Telt> const& data, std::vector<StabChain<Telt>> & L, std::vector<StabChain<Telt>> & R)
{
  imageType<Telt> image;
  Face orB_sing; // backup of <orb>. We take a single entry. Not sure it is correct
  int nrback;
  permPlusBool<Telt> rep;
  std::vector<Face> orb;
  std::vector<Face> org; // intersected (mapped) basic orbits of <G>
  Tplusinfinity<int> blen(true, 0);
  int dd, branch; // branch is level where $Lstab\ne Rstab$ starts
  std::vector<int> range;    // range for construction of <orb>
  Partition oldcel;       // old value of <image.partition.cellno>
  auto PBEnumerate = [&](int const& d, bool const & wasTriv) -> permPlusBool<Telt> {
    permPlusBool<Telt> oldprm, oldprm2;
    int a;                // current R-base point
    permPlusBool<Telt> t; // group element constructed, to be handed upwards
    int m;                // initial number of candidates in <orb>
    int max;              // maximal number of candidates still needed
    dynamic_bitset::size_type  b;        // image of base point currently being considered
               
    if (image.perm.status != int_true) 
      return {int_fail, {}};
    image.depth = d;

    // Store the original values of <image.*>.
    int undoto = NumberCells(image.partition);
    if (image.perm.status == int_true) {
      oldcel = image.partition;
    }
    else {
      oldcel = image.partition.cellno;
      oldprm = image.perm;
    }
    if (image.level2.status != int_false)
      oldprm2 = image.perm2;
    else
      oldprm2.status = int_false;
    // Recursion comes to an end  if all base  points have been prescribed
    // images.
    if (d >= rbase.base.size()) {
      if (IsTrivialRBase(rbase)) {
	blen = rbase.base.size();
	// Do     not  add the   identity    element  in the  subgroup
	// construction.
	if (wasTriv) {
	  // In the subgroup case, assign to  <L> and <R> stabilizer
	  // chains when the R-base is complete.
	  StabChainOptions<Tint> options = GetStandardOptions<Tint>();
	  options.base = rbase.base;
	  options.reduced = false;
	  L = StabChainOp(L, options);
	  R = L;
	  return {int_fail,{}};
	}
	else {
	  permPlusBool<Telt> prm;
	  if (image.perm.status == int_true)
	    prm = {int_perm, MappingPermListList(rbase.fix[rbase.base.size()-1], Fixcells(image.partition))};
	  else
	    prm = image.perm;
	  if (image.level2.status != int_false) {
	    if (SiftedPermutation(image.level2, prm.val * Inverse(image.perm2.val)).isIdentity()) 
	      return {int_perm, prm};
	  }
	  else {
	    if (Pr(prm.val))
	      return {int_perm, prm.val};
	  }
	  return {int_fail, {}};
	}
	// Construct the   next refinement  level. This  also  initializes
	// <image.partition> for the case ``image = base point''.
      }
      else {
	if (!repr) {
	  oldcel = StructuralCopy( oldcel );
	}
	NextRBasePoint(rbase.partition, rbase);
	if (image.perm.status == int_true)
	  rbase.fix.push_back(Fixcells(rbase.partition));
	Face eNewF(range.size());
	org.push_back(eNewF);
	if (repr) {
	  // In  the representative  case,  change  the   stabilizer
	  // chains of <L> and <R>.
	  ChangeStabChain(L, d, {rbase.base[d]}, int_false);
	  //	  L[ d + 1 ] := L[ d ].stabilizer;
	  ChangeStabChain(R, d, {rbase.base[d]}, int_false);
	  //	  R[ d + 1 ] := R[ d ].stabilizer;
	}
      }
    }
    a = rbase.base[d];
    
    // Intersect  the current cell of <P>  with  the mapped basic orbit of
    // <G> (and also with the one of <H> in the intersection case).
    if (image.perm.status == int_true) {
      orb[ d ] = BlistList(range, Cell(oldcel, rbase.where[d]) );
      if (image.level2.status != int_false) {
	b = orb[d].find_first();
	while (b != dynamic_bitset::npos) {
	  if (!IsInBasicOrbit(rbase.lev2[d], SlashAct(b, image.perm2.val))) 
	    orb[d][b] = false;
	  b = orb[d].find_next(b);
	}
      }
    }
    else {
      orb[d] = BlistList(range, {});
      for (auto & pVal : rbase.lev[d].orbit) {
	b = PowAct(p, image.perm.val);
	if (oldcel[b] == rbase.where[d] && (image.level2.status == int_false || IsInBasicOrbit(rbase.lev2[d], SlashAct(b,image.perm2.val)))) {
	  orb[d][b] = true;
	  org[d][b] = p;
	}
      }
    }
    if (d == 1 && ForAll(G.labels, [](Telt const& x){return PowAct(a, x) == a;})) {
      orb[d][a]=true; // ensure a is a possible image (can happen if acting on permutations with more points)
    }
    orB_sing = orb[d];
    
    // Loop  over the candidate images  for the  current base point. First
    // the special case image = base up to current level.
    if (wasTriv) {
      image.bimg[ d ] = a;
      // Refinements that start with '_' must be executed even when base
      // = image since they modify image.data, etc.
      RRefine(rbase, image, true);
      // Recursion.
      PBEnumerate(d + 1, true);
      image.depth = d;
      // Now we  can  remove  the  entire   <R>-orbit of <a>  from   the
      // candidate list.
      SubtractBlist(orb[d], BlistList(range, L.stabilizer[d].orbit));
    }
    
    // Only the early points of the orbit have to be considered.
    m = SizeBlist( orB_sing );
    if (m < L.stabilizer[d].orbit.size())
      return {int_fail,{}};
    max = PositionNthTrueBlist(orB_sing, m - L.stabilizer[d].orbit.size());
    
    if (wasTriv && a > max) {
      m--;
      if (m < L.stabilizer[d].orbit.size() )
	return {int_fail,{}};
      max = PositionNthTrueBlist( orB_sing, m - L.stabilizer[d].orbit.size());
    }
    // Now the other possible images.
    b = orb[d].find_first();
    while (b != dynamic_bitset::npos) {
      // Try to prune the node with prop 8(ii) of Leon paper.
      if (!repr && !wasTriv) {
	dd = branch;
	while (dd < d) {
	  if (IsInBasicOrbit(L, dd, a) && !PBIsMinimal(range, R.stabilizer[dd].orbit[0], b, R, d))
	    dd = d + 1;
	  else
	    dd = dd + 1;
	}
      }
      else {
	dd = d;
      }
      if (dd == d) {
	// Undo the  changes made to  <image.partition>, <image.level>
	// and <image.perm>.
	for (int i=undoto+1; i<NumberCells(image.partition); i++) 
	  UndoRefinement(image.partition);
	if (image.perm.status != 1) {
	  image.level = rbase.lev[ d ];
	  image.perm = oldprm;
	}
	if (image.level2.status != 0) {
	  image.level2 = rbase.lev2[ d ];
	  image.perm2  = oldprm2;
	}
	// If <b> could not be prescribed as image for  <a>, or if the
	// refinement was impossible, give up for this image.
	image.bimg[ d ] = b;
	IsolatePoint( image.partition, b );
	
	if (ProcessFixpoint_image(image, a, b, org[d][b]))
	  t = RRefine(rbase, image, false);
	else
	  t = int_fail;
	
	if (t != int_fail) {
	  // Subgroup case, base <> image   at current level:   <R>,
	  //   which until now is identical to  <L>, must be changed
	  //   without affecting <L>, so take a copy.
	  if (wasTriv && TestEqualityAtLevel(L, R, d)) {
	    SetStabChainFromLevel(R, L, d);
	    branch = d;
	  }
	  if (2 * d <= blen) {
	    ChangeStabChain(R, d, {b}, int_false);
	    //	    R[ d + 1 ] = R[ d ].stabilizer;
	  }
	  else {
	    std::vector<Telt> LGen = StrongGeneratorsStabChain( R, d);
	    std::vector<Telt> LGenB = Filtered(LGen, [&](Telt const& gen) -> bool {return PowAct(b, gen) == b;});
	    //	    R[ d + 1 ] := rec( generators := Filtered( R[ d + 1 ], gen -> b ^ gen = b ) );
	    int largMov=LargestMovedPoint(LGenB);
	    StabChainOptions<Tint> options = GetStandardOptions<Tint>();
	    options.base = ClosedInterval(0, largMov);
	    StabChainStrong(R, d+1, LGenB, options);
	  }
	}
        
	// Recursion.
	if (t.status == int_true) {
	  t = PBEnumerate(d + 1, false);
	  nrback++;
	  image.depth = d;
	}
                    
	// If   <t>   =   fail, either   the   recursive   call  was
	//   unsuccessful,  or all new  elements   have been added  to
	//   levels  below  the current one   (this happens if  base =
	//   image up to current level).
	if (t.status != int_fail) {
	  // Representative case, element found: Return it.
	  // Subgroup case, base <> image  before current level:  We
	  //   need  only find  a representative  because we already
	  //   know the stabilizer of <L> at an earlier level.
	  if (repr || !wasTriv)
	    return t;
	  else {
	    // Subgroup case, base  <> image at current level: Enlarge
	    //   <L>    with  <t>. Decrease <max>     according to the
	    //   enlarged <L>. Reset <R> to the enlarged <L>.
	    //	    for (int dd=0; dd<d; dd++)
	    //	      AddGeneratorsExtendSchreierTree( L[ dd ], {t});
	    AddGeneratorsExtendSchreierTree(L, 0, {t});
	    if (m < int(L.stabilizer[d].orbit.size()))
	      return {int_fail,{}};
	    max = PositionNthTrueBlist( orB_sing, m - L.stabilizer[d].orbit.size());
	    SetStabChainFromLevel(R, L, d);
	  }
	}
        
	// Now  we can remove the   entire <R>-orbit  of <b> from  the
	// candidate list.
	if  (R.stabilizer[d].translabels[b] != -1)
	  SubtractBlist(orb[d], BlistList(range, R.stabilizer[d].orbit));
	else
	  SubtractBlistOrbitStabChain(orb[d], StrongGeneratorsStabChain(R, d), b);
	b = orb[d].find_next(b);
      }
      
    }
    return {int_fail, {}};
  }

  nrback=0; // count the number of times we jumped up

  // Trivial cases first.
  if (IsTrivial(G)) {
    if (!repr)
      return {int_group, G, {}};
    if (Pr(G.identity))
      return {int_perm,{},F.identity};
    else
      return {int_fail,{},{}};
  }
    
  // Construct the <image>.
  imageType<Telt> image;
  image.data=data;
  image.depth=1;
  if (repr) {
    image.partition = data[0];
  }
  else {
    image.partition = rbase.partition;
  }
  if (IsBool(rbase.level2)) {
    image.level2 = 0;
  }
  else {
    image.level2 = rbase.level2;
    image.perm2  = rbase.level2.identity;
  }
    
  // If  <Pr> is  function,   multiply  permutations. Otherwise, keep   them
  // factorized.
  image.perm = G.identity;
  //  image.level = rbase.chain;
    
  if (repr) {
    // In the representative case, map the  fixpoints of the partitions at
    // the root of the search tree.
    if (rbase.partition.lengths != image.partition.lengths) {
      image.perm.status=0;
    }
    else {
      std::vector<int> fix  = Fixcells(rbase.partition);
      std::vector<int> fixP = Fixcells(image.partition);
      for (int i=0; i<int(fix.size()); i++)
	ProcessFixpoint_image(image, fix[i], fixP[i]);
    }
    // In   the representative case,   assign  to <L>  and <R>  stabilizer
    // chains.
    //    L := ListStabChain( CopyStabChain( StabChainMutable( L ) ) );
    //    R := ListStabChain( CopyStabChain( StabChainMutable( R ) ) );
  }
    
  int lenD=rbase.domain[rbase.domain.size()-1];
  for (int i=0; i<lenD; i++)
    range.push_back(i);
  permPlusBool<Telt> rep = PBEnumerate(0, !repr);
  if (!repr) {
    return {int_group, L, {}};
  }
  else {
    if (rep.status == int_perm)
      return {int_perm, {}, rep.val};
    else
      return {int_fail, {}, {}};
  }
}

bool IsSubset(Face const& f, std::vector<int> const& g)
{
  for (auto & eVal : g)
    if (f[eVal] == 0)
      return false;
  return true;
}


Face Difference_face(Face const& Phi, std::vector<int> const& Omega)
{
  Face fRet = Phi;
  for (auto & eVal : Omega)
    fRet[eVal] = 0;
  return fRet;
}


template<typename Telt>
Face OnSets(Face const& f, Telt const& g)
{
  int n=f.size();
  Face fRet(n);
  int b = f.find_first();
  while (b != dynamic_bitset::npos) {
    int eImg=g.at(b);
    fRet[eImg]=1;
    b = f.find_next(b);
  }
  return fRet;
}


template<typename Telt>
ResultPBT<Telt> RepOpSetsPermGroup(StabChain<Telt> const& G, bool const& repr, Face const& Phi, Face const& Psi)
{
  std::string<int> Omega = MovedPoints(G.labels);
  if (repr && Phi.size() != Psi.size())
    return {int_fail, {}, {}};
  if (IsSubset(Phi, Omega) || ForAll(Omega, [&](int const &p) -> bool {return !Phi[p];})) {
    if (repr) {
      if (Difference(Phi, Omega) != Difference(Psi, Omega))
	return {int_fail, {}, {}};
      else
	return {int_perm,{},F.identity};
    }
    else {
      return {int_group, G, {}};
    }
  }
  else {
    if (repr && (IsSubset(Psi, Omega) || ForAll(Omega, [&](int const& p) -> bool {return !Psi[p]})))
      return {int_fail, {}, {}};
  }
  auto GetPartitionFromPair=[&](std::vector<int> const& Omega, Face const& Ph) -> std::vector<std::vector<int>> {
    std::vector<int> IntVect;
    for (auto & eVal : Omega)
      if (Face[eVal] == 1)
	IntVect.push_back(eVal);
    std::vector<int> UnionVect;
    Face Ph_copy = Ph;
    for (auto & eVal : Omega)
      Ph_copy[eVal]=1;
    return Partition({IntVect, FaceToVector(Ph_copy)});
  };



  
  Partition P = GetPartitionFromPair(Phi);
  Partition Q = GetPartitionFromPair(Psi);


  auto GetSubgroup=[&](Face const& Ph) -> StabChain<Telt> {
    std::vector<Telt> sgs=Filtered(StrongGeneratorsStabChain(G), [&](Telt const& g)->bool{return OnSets(Ph, g) == Ph;});
    return MinimalStabChain(sgs);
  };
  
  
  StabChain<Telt> L = GetSubgroup(Phi);
  StabChain<Telt> R;
  if (repr)
    R = GetSubgroup(Psi);
  else
    R=L;
  rbaseType<Telt> rbase = EmptyRBase({G, G}, true, Omega, P);
  std::vector<int> Phi_vect = FaceToVector(Phi);
  std::function<bool(Telt const&)> Pr=[&](Telt const& gen) -> bool {
    for (auto & i : Phi_vect) {
      int iImg=gen.at(i);
      if (Psi[iImg] == 0)
	return false;
    }
    return true;
  };
  return PartitionBacktrack( G, Pr, repr, rbase, [ Q ], L, R );
}

}

#endif

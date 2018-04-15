#ifndef INCLUDE_PERMUTATION_GROUP
#define INCLUDE_PERMUTATION_GROUP

#include <vector>

namespace perm {


struct DoubleSidedPerm {
public:
  DoubleSidedPerm ()
  {
    siz=0;
    ListVal = {};
    ListRev = {};
  }
  DoubleSidedPerm (int const& n)
  {
    siz=n;
    ListVal = std::vector<int>(n);
    ListRev = std::vector<int>(n);
    for (int i=0; i<n; i++) {
      ListVal[i]=i;
      ListRev[i]=i;
    }
  }
  DoubleSidedPerm(std::vector<int> const& v)
  {
    ListVal=v;
    siz=v.size();
    ListRev.resize(siz);
    for (int i=0; i<siz; i++)
      ListRev[v[i]]=i;
  }
  DoubleSidedPerm(std::vector<int> const& v1, std::vector<int> const& v2)
  {
    siz=v1.size();
    ListVal=v1;
    ListRev=v2;
  }
  ~DoubleSidedPerm() 
  {
  }
  DoubleSidedPerm(DoubleSidedPerm const& ePerm)
  {
    siz=ePerm.siz;
    ListVal=ePerm.ListVal;
    ListRev=ePerm.ListRev;
  }
  DoubleSidedPerm operator=(DoubleSidedPerm const& ePerm)
  {
    siz=ePerm.siz;
    ListVal=ePerm.ListVal;
    ListRev=ePerm.ListRev;
    return *this;
  }
  bool isIdentity() const
  {
    for (int i=0; i<siz; i++)
      if (ListVal[i] != i)
	return false;
    return true;
  }
  int at(int const& i) const
  {
    return ListVal[i];
  }
  int atRev(int const& i) const
  {
    return ListRev[i];
  }
  std::vector<int> getListVal() const
  {
    return ListVal;
  }
  std::vector<int> getListRev() const
  {
    return ListRev;
  }
  int operator[](int const& i) const
  {
    return ListVal[i];
  }
  int size() const
  {
    return siz;
  }
  //
private:
  std::vector<int> ListVal;
  std::vector<int> ListRev;
  int siz;
};


bool operator==(DoubleSidedPerm const& v1, DoubleSidedPerm const& v2)
{
  int siz=v1.size();
  if (siz != v2.size() )
    return false;
  for (int i=0; i<siz; i++)
    if (v1.at(i) != v2.at(i))
      return false;
  return true;
}

bool operator<(DoubleSidedPerm const& v1, DoubleSidedPerm const& v2)
{
  int siz1=v1.size();
  int siz2=v2.size();
  if (siz1 != siz2) 
    return siz1<siz2;
  int siz=siz1;
  for (int i=0; i<siz; i++) {
    if (v1.at(i) != v2.at(i))
      return v1.at(i) < v2.at(i);
  }
  return false;
}
 
DoubleSidedPerm operator~(DoubleSidedPerm const& ePerm)
{
  return DoubleSidedPerm(ePerm.getListRev(), ePerm.getListVal());
}


// Form the product v1 * v2
DoubleSidedPerm operator*(DoubleSidedPerm const& v1, DoubleSidedPerm const& v2)
{
  int siz=v1.size();
#ifdef DEBUG
  if (siz != v2.size() ) {
    std::cerr << "Error in the DoubleSidedPerm product\n";
    throw TerminalException{1};
  }
#endif  
  std::vector<int> vVal(siz), vRev(siz);
  for (int i=0; i<siz; i++) {
    int j=v1.at(i);
    int k=v2.at(j);
    vVal[i]=k;
    //
    int j2=v2.atRev(i);
    int k2=v1.atRev(j2);
    vRev[i]=k2;
  }
  return DoubleSidedPerm(vVal, vRev);
}

DoubleSidedPerm Conjugation(DoubleSidedPerm const& v1, DoubleSidedPerm const& v2)
{
  int siz=v1.size();
#ifdef DEBUG
  if (siz != v2.size() ) {
    std::cerr << "Error in the DoubleSidedPerm conjugation\n";
    throw TerminalException{1};
  }
#endif  
  std::vector<int> v(siz);
  for (int i=0; i<siz; i++) {
    int j=v1[i];
    int i2=v2[i];
    int j2=v2[j];
    v[i2]=j2;
  }
  return DoubleSidedPerm(v);
}

int PowAct(int const& i, DoubleSidedPerm const& g)
{
  return g.at(i);
}

int SlashAct(int const& i, DoubleSidedPerm const& g)
{
  return g.atRev(i);
}


// LeftQuotient(x,y) = x^{-1}*y in the list.gi file
DoubleSidedPerm LeftQuotient(DoubleSidedPerm const& a, DoubleSidedPerm const& b)
{
  int siz=a.size();
  std::vector<int> ListVal(siz), ListRev(siz);
  for (int i=0; i<siz; i++) {
    int i1=a.atRev(i);
    int j1=b.at(i1);
    ListVal[i]=j1;
    int i2=b.atRev(i);
    int j2=a.at(i2);
    ListRev[i]=j2;
  }
  return DoubleSidedPerm(ListVal, ListRev);
}


DoubleSidedPerm SCRandomPerm(int const& d)
{
  std::vector<int> rnd(d);
  for (int i=0; i<d; i++)
    rnd[i]=i;
  for (int i=0; i<d; i++) {
    int idx=d-i;
    int res=d-i;
    int k=rand() % res;
    if (k != idx) {
      int tmp=rnd[idx];
      rnd[idx]=rnd[k];
      rnd[k]=tmp;
    }
  }
  return DoubleSidedPerm(rnd);
}


DoubleSidedPerm PermList(std::vector<int> const& eList)
{
  return DoubleSidedPerm(eList);
}

DoubleSidedPerm Inverse(DoubleSidedPerm const& ePerm)
{
  return ~ePerm;
}

std::string GapStyleString(DoubleSidedPerm const& ePerm)
{
  int n=ePerm.size();
  std::vector<int> ListStat(n,1);
  std::string eRet;
  /*
  for (int i=0; i<n; i++) {
    eRet += " " + std::to_string(ePerm.at(i));
  }
  eRet += "\n";*/

  
  for (int i=0; i<n; i++) {
    if (ListStat[i] == 1) {
      int eFirst=i;
      int eCurr=i;
      std::string ePart = "(";
      bool IsFirst=true;
      int len=0;
      while(true) {
	if (!IsFirst)
	  ePart += ",";
	IsFirst=false;
	ePart += std::to_string(eCurr);
	ListStat[eCurr]=0;
	int eNext = ePerm.at(eCurr);
	len++;
	if (eNext == eFirst)
	  break;
	eCurr = eNext;
      }
      ePart += ")";
      if (len > 1)
	eRet += ePart;
    }
  }
  if (eRet.size() > 0)
    return eRet;
  return "()";
}

 

  
}

#endif

#ifndef INCLUDE_FACET_FUNCTIONALITY
#define INCLUDE_FACET_FUNCTIONALITY

#include <locale>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>


// Code obtained from Stackoverflow
// http://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c
std::string currentDateTime() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
  // for more information about date/time format
  strftime(buf, sizeof(buf), "%Y-%m-%d %X:", &tstruct);
  return buf;
}

// Code obtained from Stackoverflow
// http://stackoverflow.com/questions/1391746/how-to-easily-indent-output-to-ofstream
class IndentFacet: public std::codecvt<char,char,std::mbstate_t>
{
 public:
  explicit IndentFacet(size_t ref = 0): std::codecvt<char,char,std::mbstate_t>(ref)    {}

  typedef std::codecvt_base::result               result;
  typedef std::codecvt<char,char,std::mbstate_t>  parent;
  typedef parent::intern_type                     intern_type;
  typedef parent::extern_type                     extern_type;
  typedef parent::state_type                      state_type;

  int&    state(state_type& s) const          {return *reinterpret_cast<int*>(&s);}
 protected:
  virtual result do_out(state_type& tabNeeded,
			const intern_type* rStart, const intern_type*  rEnd, const intern_type*&   rNewStart,
			extern_type*       wStart, extern_type*        wEnd, extern_type*&         wNewStart) const
  {
    result  res = std::codecvt_base::noconv;

    for(;(rStart < rEnd) && (wStart < wEnd);++rStart,++wStart)
      {
	// 0 indicates that the last character seen was a newline.
	// thus we will print a tab before it. Ignore it the next
	// character is also a newline
	if ((state(tabNeeded) == 0) && (*rStart != '\n'))
	  {
	    res                 = std::codecvt_base::ok;
	    state(tabNeeded)    = 1;
	    *wStart             = '\t';
	    ++wStart;
	    if (wStart == wEnd)
	      {
		res     = std::codecvt_base::partial;
		break;
	      }
	  }
	// Copy the next character.
	*wStart         = *rStart;

	// If the character copied was a '\n' mark that state
	if (*rStart == '\n')
	  {
	    state(tabNeeded)    = 0;
	  }
      }

    if (rStart != rEnd)
      {
	res = std::codecvt_base::partial;
      }
    rNewStart   = rStart;
    wNewStart   = wStart;

    return res;
  }

  // Override so the do_out() virtual function is called.
  virtual bool do_always_noconv() const noexcept
  {
    return false;   // Sometime we add extra tabs
  }
};




/*
class LogFacet: public std::codecvt<char,char,std::mbstate_t>
{
 public:
  explicit LogFacet(size_t ref = 0) : std::codecvt<char,char,std::mbstate_t>(ref) {}
  typedef std::codecvt_base::result               result;
  typedef std::codecvt<char,char,std::mbstate_t>  parent;
  typedef parent::intern_type                     intern_type;
  typedef parent::extern_type                     extern_type;
  typedef parent::state_type                      state_type;
  int&    state(state_type& s) const          {return *reinterpret_cast<int*>(&s);}
 protected:
  virtual result do_out(state_type& tabNeeded,
			const intern_type* rStart, const intern_type*  rEnd, const intern_type*&   rNewStart,
			extern_type*       wStart, extern_type*        wEnd, extern_type*&         wNewStart) const
  {
    result  res = std::codecvt_base::noconv;
    for(;(rStart < rEnd) && (wStart < wEnd);++rStart,++wStart) {
      // 0 indicates that the last character seen was a newline.
      // thus we will print a tab before it. Ignore it the next
      // character is also a newline
      if ((state(tabNeeded) == 0) && (*rStart != '\n')) {
	res                 = std::codecvt_base::ok;
	state(tabNeeded)    = 1;
	std::string strDT=currentDateTime();
	size_t siz=strDT.size();
	for (size_t i=0; i<siz; i++) {
	  *wStart = strDT.at(i);
	  ++wStart;
	}
	*wStart             = '\t';
	++wStart;
	if (wStart == wEnd) {
	  res     = std::codecvt_base::partial;
	  break;
	}
      }
      // Copy the next character.
      *wStart         = *rStart;
      // If the character copied was a '\n' mark that state
      if (*rStart == '\n')
	state(tabNeeded) = 0;
    }
    if (rStart != rEnd)
      res = std::codecvt_base::partial;
    rNewStart   = rStart;
    wNewStart   = wStart;
    return res;
  }
  // Override so the do_out() virtual function is called.
  virtual bool do_always_noconv() const throw()
  {
    return false;   // Sometime we add extra tabs
  }
};
*/




class LogFacet: public std::codecvt<char,char,std::mbstate_t>
{
 public:
  explicit LogFacet(size_t ref = 0) : std::codecvt<char,char,std::mbstate_t>(ref) {}
  typedef std::codecvt_base::result               result;
  typedef std::codecvt<char,char,std::mbstate_t>  parent;
  typedef parent::intern_type                     intern_type;
  typedef parent::extern_type                     extern_type;
  typedef parent::state_type                      state_type;
  int&    state(state_type& s) const          {return *reinterpret_cast<int*>(&s);}
 protected:
  result do_out(state_type& tabNeeded,
			const intern_type* rStart, const intern_type*  rEnd, const intern_type*&   rNewStart,
			extern_type*       wStart, extern_type*        wEnd, extern_type*&         wNewStart) const
  {
    result  res = std::codecvt_base::noconv;
    for(;(rStart < rEnd) && (wStart < wEnd);++rStart,++wStart) {
      // 0 indicates that the last character seen was a newline.
      // thus we will print a tab before it. Ignore it the next
      // character is also a newline
      if ((state(tabNeeded) == 0) && (*rStart != '\n')) {
	res                 = std::codecvt_base::ok;
	state(tabNeeded)    = 1;
	std::string strDT=currentDateTime();
	size_t siz=strDT.size();
	for (size_t i=0; i<siz; i++) {
	  *wStart = strDT.at(i);
	  ++wStart;
	}
	*wStart             = '\t';
	++wStart;
	if (wStart == wEnd) {
	  res     = std::codecvt_base::partial;
	  break;
	}
      }
      // Copy the next character.
      *wStart         = *rStart;
      // If the character copied was a '\n' mark that state
      if (*rStart == '\n')
	state(tabNeeded) = 0;
    }
    if (rStart != rEnd)
      res = std::codecvt_base::partial;
    rNewStart   = rStart;
    wNewStart   = wStart;
    return res;
  }
  // Override so the do_out() virtual function is called.
  virtual bool do_always_noconv() const noexcept
  {
    return false;   // Sometime we add extra tabs
  }
};




#endif

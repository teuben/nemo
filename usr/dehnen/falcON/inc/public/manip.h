// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// manip.h                                                                     |
//                                                                             |
// Copyright (C) 2004, 2005  Walter Dehnen                                     |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
// version 0.0  17/09/2004 WD  created based on acceleration.h v3.2            |
// version 1.0  19/05/2005 WD  allowed for manips given in single file         |
// version 1.1  20/05/2005 WD  added manpath as 4th argument to Manipulator    |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_manip_h
#define falcON_included_manip_h

#ifndef falcON_included_body_h
#  include <body.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::manipulator                                                //
  //                                                                          //
  // abstract base class                                                      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class manipulator {
  public:
    virtual const char* name    () const = 0;      // name of manipulator       
    virtual const char* describe() const = 0;      // description ----          
    virtual fieldset    need    () const = 0;      // which data are needed     
    virtual fieldset    change  () const = 0;      // which data get changed    
    virtual fieldset    provide () const = 0;      // which data are provided   
    virtual bool manipulate(const snapshot*)       // R: stop simulation?       
      const = 0;
    virtual operator bool() const { return true; }
    virtual ~manipulator() {}
    bool operator() (const snapshot*s) const { return manipulate( s); }
    bool operator() (snapshot const&s) const { return manipulate(&s); }
  };
#ifdef falcON_NEMO
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::Manipulator                                                //
  //                                                                          //
  //                                                                          //
  // NOTE on "concatination" of manipulators                                  //
  // We allow several (unlimited) manipulators to "line up". Manipulators are //
  // "concatinated" by exectuting them in the order listed. The constructor   //
  // will issue a warning if the needs of the first Manipulator exceed mxvapq.//
  //                                                                          //
  //                                                                          //
  // NOTE on construction of Manipulator                                      //
  // The constructor Manipulator::Manipulator() has 4 arguments: NAME, PARS,  //
  // FILE, and, optionally, a PATH. There two ways to pass the names,         //
  // parameter sets, and files.                                               //
  //                                                                          //
  //                                                                          //
  // 1 By argument                                                            //
  // If NAME is given and its first character is not a '+', then it is        //
  // assumed that it refers to the names of the to be concatinated            //
  // manipulators in the format                                               //
  //                                                                          //
  //   name1[+name2[+mane3 ... ]]                                             //
  //                                                                          //
  // Instead of '+' also ',' is allowed as separator. If given by PARS, the   //
  // parameter sets must match in number the names and must be in the format  //
  //                                                                          //
  //   pars1[#pars2[#pars3 ... ]]                                             //
  //                                                                          //
  // with each parameter set being a comma separated list of numbers (no      //
  // space) which may be empty. Instead of '#', also ';' is allowed as        //
  // separator. If given by FILE, the files must match in number the names    //
  // and must be in the format                                                //
  //                                                                          //
  //   file1[#file2[#file3 ... ]]                                             //
  //                                                                          //
  // Instead of '#' also ';' is allowed as separator.                         //
  //                                                                          //
  //                                                                          //
  // 2 By file                                                                //
  // If NAME is empty or starts with '+', then it is assumed that FILE refers //
  // to a file which contains the information about manipulator names,        //
  // parameters, and files. Each line of that file whose first non-space      //
  // character is not a '#' is assumed to be of the format                    //
  //                                                                          //
  //   name [pars] [file] [# comment]                                         //
  //                                                                          //
  // where pars is a comma separated list of numbers (no spaces) and the      //
  // first character of file must be a letter or '/', BUT NOT a '.' (these    //
  // rules are to distinguish a parameter set from a file name).              //
  //                                                                          //
  // Finally, if the optional argument PATH is given, the the file "name.so"  //
  // will be searched for in that path before looking in "." and              //
  // "$FALCON/manip/".                                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class Manipulator : public manipulator {
    typedef const manipulator *cp_manip;
    //--------------------------------------------------------------------------
    static const int NMAX=100;
    int       N;
    cp_manip  MANIP[NMAX];
    char     *NAME,*DSCR;
    fieldset  NEED, CHNG, PRVD;
  public:
    //--------------------------------------------------------------------------
    const char* name    () const { return N? NAME : "empty"; }
    const char* describe() const { return N? DSCR : "empty"; }
    fieldset    need    () const { return NEED; }
    fieldset    change  () const { return CHNG; }
    fieldset    provide () const { return PRVD; }
    bool manipulate(const snapshot*s) const {
      bool r = false;
      for(int i=0; i!=N; ++i)
	r = r || MANIP[i]->manipulate(s);
      return r;
    }
    operator bool() const { return N != 0; }
    //--------------------------------------------------------------------------
    Manipulator(const char*,                     // I: name(s) of manipulator(s)
		const char*,                     // I: set(s) of parameters     
		const char*,                     // I: data file(s)             
		const char* =0) falcON_THROWING; //[I: path to search]          
    //--------------------------------------------------------------------------
    ~Manipulator() {
      if(N) {
	falcON_DEL_A(NAME);
	falcON_DEL_A(DSCR);
	for(int i=0; i!=N; ++i) delete MANIP[i];
      }
      N = 0;
    }
    //--------------------------------------------------------------------------
  };
#endif // falcON_NEMO
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_manip_h

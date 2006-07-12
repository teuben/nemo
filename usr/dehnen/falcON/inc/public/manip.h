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
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::manipulator                                                //
  //                                                                          //
  /// abstract base class                                                     //
  //                                                                          //
  // ///////////////////////////////////////////////////////////////////////////
  class manipulator {
  public:
    /// name of manipulator
    virtual const char* name() const = 0;
    /// description of manipulator
    virtual const char* describe() const = 0;
    /// which data are needed by manipulator
    virtual fieldset need() const = 0;
    /// which data are changed by manipulator
    virtual fieldset change() const = 0;
    /// which data are provided by manipulator
    virtual fieldset provide() const = 0;
    /// manipulation of snapshot
    /// \return shall the simulation be stopped?
    virtual bool manipulate(const snapshot*) const = 0;
    /// conversion to bool: is this manipulator doing anything?
    virtual operator bool() const { return true; }
    /// destruction
    virtual ~manipulator() {}
    /// manipulation via function call
    /// \return shall the simulation be stopped?
    bool operator() (const snapshot*s) const { return manipulate( s); }
    /// manipulation via function call
    /// \return shall the simulation be stopped?
    bool operator() (snapshot const&s) const { return manipulate(&s); }
  };
#ifdef falcON_NEMO
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::Manipulator                                                  
  //                                                                            
  /// A manipulator which on construction loads other manipulators              
  ///                                                                           
  /// Several (unlimited) manipulator may be concatinated, i.e. a call to       
  /// Manipulator::manipulate() will call the individual manipulate() in the    
  /// order given at construction. The constructor will issue a warning if the  
  /// needs of the first manipulator exceed mxvapq. See also details for        
  /// Manipulator::Manipulator().                                               
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
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
    /// name of manipulator(s) loaded
    const char* name    () const { return N? NAME : "empty"; }
    //--------------------------------------------------------------------------
    /// descriptions of manipulator(s) loaded
    const char* describe() const { return N? DSCR : "empty"; }
    //--------------------------------------------------------------------------
    /// combined need of manipulators loaded
    fieldset    need    () const { return NEED; }
    //--------------------------------------------------------------------------
    /// data changed by manipulator(s) loaded
    fieldset    change  () const { return CHNG; }
    //--------------------------------------------------------------------------
    /// data provided by manipulator(s) loaded
    fieldset    provide () const { return PRVD; }
    //--------------------------------------------------------------------------
    /// manipulate: call individual manipulation in order
    /// \return: stop simulation?
    bool manipulate(const snapshot*s) const {
      int r = 0;
      for(int i=0; i!=N; ++i)
	r |= MANIP[i]->manipulate(s);
      return r;
    }
    //--------------------------------------------------------------------------
    /// are there any manipulators loaded?
    operator bool() const { return N != 0; }
    //--------------------------------------------------------------------------
    /// construction from character strings                                     
    ///                                                                         
    /// \param names name(s) of manipulator(s) to be loaded                     
    /// \param parss parameter set(s) for manipulator(s) to be loaded           
    /// \param files data file(s) for the manipulator(s) to be loaded           
    /// \param search path to search for .so of manipulator(s)                  
    ///                                                                         
    /// There are two different ways to pass the names, parameter sets and      
    /// file names.                                                             
    ///                                                                         
    /// First, by argument: If \a names is given (non-null) and its first       
    /// character is not an '+', then it is assumes that it refers to the names 
    /// of the to concatinated manipulators in the format                       
    ///                                                                         
    /// name1[+name2[+name3 ... ]]                                              
    ///                                                                         
    /// Instead of '+' as separator, also ',' is allowed. The corresponding     
    /// parameter sets must match in number the names and must be in the format 
    ///                                                                         
    /// pars1[#pars2[#pars3 ... ]]                                              
    ///                                                                         
    /// with each parameter set being a comma seperated list of numbers (no     
    /// space) with may be empty. Instead of '#', also ';' is allowed as        
    /// separator. The data files must match in number the names                
    /// and must be in the format                                               
    ///                                                                         
    /// file1[#file2[#file3 ... ]]                                              
    ///                                                                         
    /// where again ';' is also allowed as a separator.                         
    ///                                                                         
    /// Second, the names, parameter sets, and data files can be read from      
    /// a single file. This is done if either \a names is null or starts with   
    /// an '+'. In this case, \a files is interpreted as the name of a file     
    /// containing lines in the format                                          
    ///                                                                         
    /// name [pars] [file] [#comment]                                           
    ///                                                                         
    /// (lines the first non-white character of which is '#' are ignored).      
    /// Here, pars must be a comma separated list of values and file must       
    /// start with a letter or a '/', but not a '.' (these rules are to         
    /// distinguish a parameter set from a file name).                          
    ///                                                                         
    /// Finally, if the optional argument PATH is given, the the file "name.so" 
    /// will be searched for in that path before looking in "." and             
    /// "$FALCON/manip/".                                                       
    // /////////////////////////////////////////////////////////////////////////
    Manipulator(const char*names,
		const char*parss,
		const char*files,
		const char*path =0) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// destructor: destruct manipulators loaded
    ~Manipulator() {
      if(N) {
	falcON_DEL_A(NAME);
	falcON_DEL_A(DSCR);
	for(int i=0; i!=N; ++i) falcON_DEL_O(MANIP[i]);
      }
      N = 0;
    }
    //--------------------------------------------------------------------------
    /// auxiliary static function for parsing parameters to array of double     
    ///                                                                         
    /// \return number of parsed parameters                                     
    /// \param params (input) string with parameters                            
    /// \param pars   (output) array with parsed parameters                     
    /// \param maxp   (input) maximum number of parsed parameters               
    static int parse(const char*params, double*pars, int maxp);
#define MANIP_PARSE_AT_INIMANIP
  };
#endif // falcON_NEMO
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_manip_h

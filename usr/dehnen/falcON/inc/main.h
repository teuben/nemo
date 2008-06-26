// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// main.h                                                                      |
//                                                                             |
// Copyright (C) 2002-2008  Walter Dehnen                                      |
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
//                                                                             |
// based on nemomain.c by Peter Teuben                                         |
//                                                                             |
// to be included by mains only                                                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// to run your code with MPI, you must                                         |
//                                                                             |
// have falcON_MPI       #defined                                              |
// have falcON_PARALLEL  #defined                                              |
//                                                                             |
// falcON_MPI       shall be #defined (via a compiler option), if MPI is       |
//                  installed and to be enabled                                |
// falcON_PARALLEL  shall be #defined in the including .cc file, if the code   |
//                  is supposed to be a parallel one. Otherwise it will become |
//                  an ordinary serial code.                                   |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// to run your code with NEMO, you must                                        |
//                                                                             |
// have     falcON_NEMO      #defined                                          |
// have not falcON_NONEMO    #defined                                          |
//                                                                             |
// falcON_NEMO      shall be #defined (via a compiler option), if NEMO is      |
//                  installed and to be enabled                                |
// falcON_NONEMO    shall be #defined in the including .cc file, if the code   |
//                  is not to be using nemo features.                          |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifdef falcON_included_main_h
# error "main.h must only be included by a falcON::main"
#else
#define falcON_included_main_h 1

#if ! (defined(__cplusplus) || defined(c_plusplus) )
# error "main.h must be used as C++ main only"
#endif

#ifdef unix
#  ifndef falcON_included_unistd_h
#    include <unistd.h>
#    define falcON_included_unistd_h
#  endif
#else
#  ifndef falcON_included_fstream
#    include <fstream>
#    define falcON_included_fstream
#  endif
#endif
//------------------------------------------------------------------------------
// In case, this is a MPI application                                           
//                                                                              
// - set macro "falcON_USE_MPI"                                                 
// - include MPI basic stuff                                                    
//------------------------------------------------------------------------------

#if (defined(falcON_MPI) && defined(falcON_PARALLEL) ) 
                                                   // have MPI and want it?     
#  ifndef falcON_USE_MPI
#    define falcON_USE_MPI                         //   then use it             
#  endif
#  ifndef falcON_included_mpi_falcON_h
#    include <parallel/mpi_falcON.h>
#  endif

#else                                              // else                      

#  undef falcON_USE_MPI                            //   don't use it            

#endif

//------------------------------------------------------------------------------
// define macro indicating which compiler flags have been set                   
//------------------------------------------------------------------------------
#ifdef falcON_PROPER
#  define falcON_Pflag  "P"
#else
#  define falcON_Pflag
#endif

#ifdef falcON_SSE
#  define falcON_Sflag  "S"
#else
#  define falcON_Sflag
#endif

#ifdef falcON_INDI
#  ifdef falcON_ADAP
#    define falcON_Iflag  "IA"
#  else
#    define falcON_Iflag  "I"
#  endif
#else
#  define falcON_Iflag
#endif

#ifdef falcON_USE_MPI
#  define falcON_Mflag  "M"
#else
#  define falcON_Mflag
#endif

#define falcON_PSIFLAG falcON_Pflag falcON_Sflag falcON_Iflag falcON_Mflag

#if  !defined(falcON_PROPER) \
  && !defined(falcON_SSE) \
  && !defined(falcON_INDI) \
  && !defined(falcON_USE_MPI)
#undef falcON_PSIFLAG
#endif
//------------------------------------------------------------------------------
// include falcON stuff                                                         
//------------------------------------------------------------------------------

#ifndef falcON_included_body_h
#  include <body.h>
#endif

//------------------------------------------------------------------------------
// implement basic.h: compile_info                                              
//------------------------------------------------------------------------------
namespace falcON { namespace compile_info {
  bool __set = 0;      bool const&is_set  () { return __set; }
  char __version[100]; const char*version () { return __version; }
  char __origin [100]; const char*origin  () { return __origin; }
  char __time    [30]; const char*time    () { return __time; }
  char __compiler[30]; const char*compiler() { return __compiler; }
  void init() {
    if(! __set ) {
      SNprintf(__compiler,10,
#if   defined (__INTEL_COMPILER)
	       "icc-%d",__INTEL_COMPILER
#elif defined (__GNUC__)
	       "gcc-%d.%d.%d",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__
#else
	       "???"
#endif
	       );
#ifdef falcON_VERSION_D
      SNprintf(__origin,100,falcON_VERSION_D);
#else
      __origin[0] = 0;
#endif
#if defined(falcON_VERSION) || defined(falcON_PSIFLAG)
      SNprintf(__version,100,
#ifdef falcON_VERSION
	       falcON_VERSION
#endif
#ifdef falcON_PSIFLAG
	       falcON_PSIFLAG
#endif
	       );
#endif
      SNprintf(__time,30,__DATE__ ", " __TIME__);
      __set = 1;
    }
  }
} }
//------------------------------------------------------------------------------
// In case, this is a NEMO application                                          
//                                                                              
// - set macro "falcON_USE_NEMO"                                                
// - include some NEMO stuff                                                    
// - make some extern declarations (to be satisfied by the user)                
// - define defv_info with various informations for possible use in defv[]      
//                                                                              
//------------------------------------------------------------------------------

#if (defined(falcON_NEMO) && !defined(falcON_NONEMO) )
#  ifndef falcON_USE_NEMO
#    define falcON_USE_NEMO                        //   then use it             
#  endif
#  include <stdinc.h>                              // NEMO basic stuff          
#  include <getparam.h>                            // NEMO parameter handling   
#  include <history.h>                             // NEMO history              
#  include <public/io.h>                           // my NEMO I/O & in/output   
#  undef local
extern string defv[];                              // MUST be supplied by user  
extern string usage;	                           // MUST be supplied by user  
namespace falcON { namespace defv_info {
  char version [100];                              //   "VERSION=..."           
  char compiled[100];                              //   "COMPILED=..."          
  void init() {
    SNprintf(version,100,"VERSION=%s\n%s",
	     compile_info::version (),
	     compile_info::origin  ());
    SNprintf(compiled,100,"COMPILED=\n%s, with %s"
	     "                                   ",
	     compile_info::time(),
	     compile_info::compiler());
    compiled[61] = 0;
  };
} }
#  if (defined(falcON_PROPER) && !defined(falcON_NOT_USING_PROPER))     \
   ||  defined(falcON_USING_PROPER)
#    define falcON_DEFV							\
     falcON::defv_info::version,					\
     falcON::defv_info::compiled,					\
     "STATUS=\nproprietary version; usage restricted              "
#  else
#    define falcON_DEFV							\
     falcON::defv_info::version,					\
     falcON::defv_info::compiled,					\
     "STATUS=\npublic version                                     "
#  endif
#endif

//------------------------------------------------------------------------------
// define main() in namespace falcON                                            
//------------------------------------------------------------------------------
// MUST be supplied by user in the #including file                              
namespace falcON {
#ifdef falcON_USE_NEMO
  extern void main(void) falcON_THROWING;
#else
  extern void main(int argc, char *argv[]) falcON_THROWING;
#endif
//   //----------------------------------------------------------------------------
//   // define falcON::ERROR() to call nemo::error() if NEMO application           
//   //----------------------------------------------------------------------------
//   inline void ERROR(const char*m) {
// #ifndef falcON_USE_NEMO
//     falcON_ErrorN(m);
// #else
//     ::error(const_cast<char*>(m));
// #endif
//   }
}

//------------------------------------------------------------------------------
// define global main(), which calls nemo::main()                               
//------------------------------------------------------------------------------

int main(int argc, char *argv[])                   // global main               
{

  try {                                            // TRY:                      

#ifdef falcON_USE_MPI
    falcON::MPI::Init(&argc,&argv);                // start MPI: spawm processes
    falcON::set_exit(&falcON::MPI::Exit);          // make sure MPI_Abort() is  
    set_nemo_exit(&falcON::MPI::Exit);             //   called on exit()        
    falcON::RunInfo::set_mpi_proc(falcON::MPI::World.rank()); // inform RunInfo 
#ifdef falcON_NEMO
    ::set_mpi_rank(falcON::MPI::World.rank());
#endif
#endif

    falcON::CheckAgainstLibrary(falcON::CurrentStatus(),
				falcON::RunInfo::name_known() ?
				falcON::RunInfo::name() : "executable");
                                                   // assert status matches     

    falcON::compile_info::init();                  // initialize compile_info   

#ifdef falcON_USE_NEMO
    falcON::defv_info::init();                     // initialize defv_info      
#ifndef falcON_NoHelpDefault
    if(argc == 1) {                                // IF no command line args   
      char *argV[3];
      argV[0] = argv[0];
      argV[1] = "help=h";
      argV[2] = 0;
      std::cout << '\n' << usage << "\n\noption summary:\n";
      initparam(argV,defv);
    } else
#endif
      initparam(argv,defv);
    {
      int d = 100;
      while(!nemo_debug(d--));
      falcON::RunInfo::set_debug_level(++d);
    }
#  if defined(falcON_PROPER) &&			\
      defined(falcON_RepAction) && (falcON_RepAction==1)

    falcON::report::open_file(argv[0],*(ask_history()));

#  endif

    try {                                          // TRY:                      
      falcON::main();                              //   call user program       
    } catch(falcON::exception E) {                 // CATCH falcON errors       
      falcON_ErrorN(text(E));
    }
    finiparam();                                   // finish NEMO               

#else // falcON_USE_NEMO

#  if defined(falcON_PROPER) &&			\
      defined(falcON_RepAction) && (falcON_RepAction==1)

    falcON::report::open_file(argv[0]);

#  endif

    falcON::main(argc,argv);                       //   call user program       

#endif

#if defined(falcON_PROPER) &&			\
  defined(falcON_RepAction) &&			\
  (falcON_RepAction==1)
    falcON::report::close_file();
#endif

  } catch(falcON::exception E) {                   // CATCH falcON errors       
    falcON_ErrorN(text(E));
  }

#ifdef falcON_USE_MPI
  falcON::MPI::Finish();                           // finish MPI                
#endif
}
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
#ifdef falcON_USE_NEMO
#undef nemoinpr
#ifdef falcON_REAL_IS_FLOAT
#  define nemoinpr nemoinpf
#else
#  define nemoinpr nemoinpd
#endif
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // some utilities for getting NEMO parameters                               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  //----------------------------------------------------------------------------
  // like NEMO::hasvalue(), except that we take a const char*
  inline bool hasvalue(const char*param) {
    return ::hasvalue(const_cast<char*>(param));
  }
  //----------------------------------------------------------------------------
  // read vect from nemo parameter                                              
  vect getvparam(const char* option) falcON_THROWING {
    vect X;
    int  N = nemoinpr(getparam(const_cast<char*>(option)),
		      static_cast<real*>(X),Ndim);
    if(N>Ndim)
      warning("option \"%s\" requires %d values, but %d given\n",option,Ndim,N);
    else if(N<0)  error("parse error: processing option \"%s\"\n",option);
    else if(N==0) error("no data for option \"%s\"\n",option);
    else if(N<Ndim)
      falcON_THROW("option \"%s\" requires %d values, but %d given\n",
		   option,Ndim,N);
    return X;
  }
  //----------------------------------------------------------------------------
  // read vect from nemo parameter into 2nd argument, return pointer            
  vect* getvparam_z(const char* option, vect&X) falcON_THROWING {
    if(!hasvalue(const_cast<char*>(option))) return 0;
    int N = nemoinpr(getparam(const_cast<char*>(option)),
		     static_cast<real*>(X),Ndim);
    if(Ndim != N) {
      if(N<0) falcON_THROW("parse error: processing option \"%s\"\n",option);
      if(N>0) warning("option \"%s\" requires %d values, but %d given\n",
		      option,Ndim,N);
      return 0;
    }
    return &X;
  }
  //----------------------------------------------------------------------------
  // read array of type T
  template<typename Type> struct __getA;
  template<> struct __getA<double> {
    static int get(const char*o, double*a, int m) falcON_THROWING {
      return nemoinpd(getparam(const_cast<char*>(o)),a,m);
    } };
  template<> struct __getA<float> {
    static int get(const char*o, float*a, int m) falcON_THROWING {
      return nemoinpf(getparam(const_cast<char*>(o)),a,m);
    } };
  template<> struct __getA<int> {
    static int get(const char*o, int*a, int m) falcON_THROWING {
      return nemoinpi(getparam(const_cast<char*>(o)),a,m);
    } };
  template<> struct __getA<unsigned> {
    static int get(const char*o, unsigned*a, int m) falcON_THROWING {
      return nemoinpi(getparam(const_cast<char*>(o)),
		      static_cast<int*>(static_cast<void*>(a)),m);
    } };
  template<> struct __getA<bool> {
    static int get(const char*o, bool*a, int m) falcON_THROWING {
      return nemoinpb(getparam(const_cast<char*>(o)),a,m);
    } };
  /// read array of type T
  /// \param o  name of nemo option
  /// \param a  array
  /// \param m  physical size of array
  /// \return   number of elements actually read.
  template<typename Type>
  int getaparam(const char*o, Type*a, int m) falcON_THROWING {
    return __getA<Type>::get(o,a,m);
  }
  /// read array of type T, but allow for non-existence
  /// \param o  name of nemo option
  /// \param a  array
  /// \param m  physical size of array
  /// \return   number of elements actually read.
  template<typename Type>
  int getaparam_z(const char*o, Type*a, int m) falcON_THROWING {
    if(!hasvalue(const_cast<char*>(o))) {
      for(int i=0; i!=m; ++i) a[i] = Type(0);
      return 0;
    } else 
      return __getA<Type>::get(o,a,m);
  }
  //----------------------------------------------------------------------------
  // read float                                                                 
  inline float getfparam(const char* option) {
    return float(getdparam(const_cast<char*>(option)));
  }
  //----------------------------------------------------------------------------
  // read real                                                                  
#ifdef getrparam
#  undef getrparam
#endif
  inline real getrparam(const char* option) { 
    return real(getdparam(const_cast<char*>(option)));
  }
  //----------------------------------------------------------------------------
  // read unsigned                                                              
  inline unsigned getuparam(const char* option) {
    int i = getiparam(const_cast<char*>(option));
    if(i < 0)
      error("getuparam(%s): expected a positive integer, got %d\n", option,i);
    return unsigned(i);
  }
  //----------------------------------------------------------------------------
  // read fieldset                                                              
  inline fieldset getioparam(const char* option) {
    return fieldset(getparam(const_cast<char*>(option)));
  }
  //----------------------------------------------------------------------------
  // define getparam_z(arg), and related, which will return 0 if !hasvalue(arg).
  //----------------------------------------------------------------------------
  inline char* getparam_z(const char* option) {
    return hasvalue(const_cast<char*>(option))?
      getparam(const_cast<char*>(option)) : 0;
  }
  //----------------------------------------------------------------------------
  inline int getiparam_z(const char* option) {
    return hasvalue(const_cast<char*>(option))?
      getiparam(const_cast<char*>(option)) : 0;
  }
  //----------------------------------------------------------------------------
  inline unsigned getuparam_z(const char* option) {
    return hasvalue(const_cast<char*>(option))?
      getuparam(const_cast<char*>(option)) : 0u;
  }
  //----------------------------------------------------------------------------
  inline long getlparam_z(const char* option) {
    return hasvalue(const_cast<char*>(option))?
      getlparam(const_cast<char*>(option)) : 0;
  }
  //----------------------------------------------------------------------------
  inline bool getbparam_z(const char* option) {
    return hasvalue(const_cast<char*>(option))?
      getbparam(const_cast<char*>(option)) : false;
  }
  //----------------------------------------------------------------------------
  inline double getdparam_z(const char* option) {
    return hasvalue(const_cast<char*>(option))?
      getdparam(const_cast<char*>(option)) : 0.;
  }
  //----------------------------------------------------------------------------
  inline float getfparam_z(const char* option) {
    return hasvalue(const_cast<char*>(option))?
      getfparam(const_cast<char*>(option)) : 0.f;
  }
  //----------------------------------------------------------------------------
  inline real getrparam_z(const char* option) {
    return hasvalue(const_cast<char*>(option))?
      getrparam(const_cast<char*>(option)) : zero;
  }
  //----------------------------------------------------------------------------
  inline fieldset getioparam_z(const char* option) {
    return hasvalue(const_cast<char*>(option))?
      getioparam(const_cast<char*>(option)) : fieldset::empty;
  }
  //----------------------------------------------------------------------------
  // read fieldset, if not given return fieldset::all                           
  inline fieldset getioparam_a(const char* option) {
    return hasvalue(const_cast<char*>(option))?
      getioparam(const_cast<char*>(option)) : fieldset::all;
  }
  //----------------------------------------------------------------------------
  // read PotExp::symmetry
#ifdef falcON_included_PotExp_h
  inline PotExp::symmetry getsymparam(const char*symm) {
    int __sym (getiparam(const_cast<char*>(symm)));
    return
      __sym==4? PotExp::spherical   :
      __sym==3? PotExp::cylindrical :
      __sym==2? PotExp::triaxial    :
      __sym==1? PotExp::reflexion   : PotExp::none;
  }
#endif
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // check for file name to relate to a real output (is given && != ".")      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  inline bool file_for_output(const char* option) {
    return
      hasvalue(const_cast<char*>(option)) &&
      strcmp  (getparam(const_cast<char*>(option)),".");
  }
#endif // falcON_USE_NEMO
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // check fieldset for completeness                                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  inline void check_sufficient(fieldset const&read, fieldset const&need)
    throw(falcON::exception) {
    if(! read.contain(need)) {
      fieldset::wlist wneed(&need);
      fieldset::wlist wread(&read);
      throw exception("insufficient data: need \'%s\' but got only \'%s\'",
		      static_cast<const char*>(wneed),
 		      static_cast<const char*>(wread) );
    }
  }
  //----------------------------------------------------------------------------
  inline bool warn_insufficient(fieldset const&read, fieldset const&need) {
    if(! read.contain(need)) {
      fieldset::wlist wneed(&need);
      fieldset::wlist wread(&read);
      falcON_Warning("insufficient data: need \'%s\' but got only \'%s\'",
		     static_cast<const char*>(wneed),
		     static_cast<const char*>(wread) );
      return true;
    }
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // check for existence of a file                                            //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  inline bool file_exists(const char*file) {
#ifdef unix
    // simply use the access system call. F_OK asks for existence only.         
    return 0 == access(file, F_OK);
#else
    // in fact, this test is different, as existence and read permission are    
    // required to open an ifstream.                                            
    std::ifstream IN;
    IN.open(file);
    if(IN.is_open()) {
      IN.close();
      return true;
    } else
      return false;
#endif
  }
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {

////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_main_h

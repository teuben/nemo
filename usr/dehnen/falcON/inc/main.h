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
#  ifndef falcON_included_parallel_h
#    include <parallel/parallel.h>
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

#ifdef falcON_ADAP
#  define falcON_Iflag  "IA"
#else
#  define falcON_Iflag  "I"
#endif

#ifdef falcON_USE_MPI
#  define falcON_Mflag  "M"
#else
#  define falcON_Mflag
#endif

#define falcON_PSIFLAG falcON_Pflag falcON_Sflag falcON_Iflag falcON_Mflag

#if  !defined(falcON_PROPER) \
  && !defined(falcON_SSE) \
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
      SNprintf(__compiler,15,
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
#  include <public/nemo++.h>                       // my NEMO interface
#  undef local
extern const char* defv[];                         // MUST be supplied by user  
extern const char* usage;	                   // MUST be supplied by user  
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
#endif // falcON_NEMO

//------------------------------------------------------------------------------
// define main() in namespace falcON                                            
//------------------------------------------------------------------------------
// MUST be supplied by user in the #including file                              
namespace falcON {
#ifdef falcON_USE_NEMO
  extern void main(void) falcON_THROWING;
#else
  extern void main(int argc, const char**argv) falcON_THROWING;
#endif
}

//------------------------------------------------------------------------------
// define global main(), which calls nemo::main()                               
//------------------------------------------------------------------------------

int main(int argc, const char**argv)               // global main               
{

  try {                                            // TRY:                      

#ifdef falcON_USE_MPI
    falcON::MPI::Init(argc,argv);                  // start MPI: spawm processes
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
      const char *argV[3];
      argV[0] = argv[0];
      argV[1] = "help=h";
      argV[2] = 0;
      std::cout << '\n' << usage << "\n\noption summary:\n";
      falcON::initparam(argV,defv);
    } else
#endif
      falcON::initparam(argv,defv);
      falcON::RunInfo::set_debug_level(falcON::nemo_debug_level());
#  if defined(falcON_PROPER) &&			\
      defined(falcON_RepAction) && (falcON_RepAction==1)
      falcON::report::open_file(argv[0],falcON::RunInfo::cmd());
#  endif

    try {                                          // TRY:                      
      falcON::main();                              //   call user program       
    } catch(falcON::exception& E) {                 // CATCH falcON errors       
      falcON_ErrorN(text(E));
    }
    falcON::finiparam();                           // finish NEMO               

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

  } catch(falcON::exception& E) {                   // CATCH falcON errors       
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
  // check for file name to relate to a real output (is given && != ".")      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  inline bool file_for_output(const char*file) {
    return hasvalue(file) && strcmp(getparam(file),".");
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

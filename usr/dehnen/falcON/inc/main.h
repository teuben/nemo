// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// main.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
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
# error "main.h must only be included by a nbdy::main"
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
#  ifndef falcON_included_mpiu_h
#    include <walter/mpiu.h>
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

//------------------------------------------------------------------------------
// include other necessary header files                                         
//------------------------------------------------------------------------------

#ifndef falcON_included_auxx_h
#  include <public/auxx.h>                         // falcON types etc          
#endif

#ifndef falcON_included_exit_h
#  include <public/exit.h>                         // falcON exit & error       
#endif

#ifndef falcON_included_nbio_h
#  include <public/nbio.h>                         // falcON body I/O support   
#endif

//------------------------------------------------------------------------------
// implement exit.h: compile_info                                               
//------------------------------------------------------------------------------
namespace nbdy { namespace compile_info {
  bool __set = 0;      bool const&is_set  () { return __set; }
  char __version[100]; const char*version () { return __version; }
  char __origin [100]; const char*origin  () { return __origin; }
  char __time    [30]; const char*time    () { return __time; }
  char __compiler[10]; const char*compiler() { return __compiler; }
  void init() {
    if(! __set ) {
      snprintf(__compiler,10,
#if   defined (__INTEL_COMPILER)
	       "icc-%d",__INTEL_COMPILER
#elif defined (__GNUC__)
	       "gcc-%d.%d",__GNUC__,__GNUC_MINOR__
#else
	       "???"
#endif
	       );
#ifdef falcON_VERSION_D
      snprintf(__origin,100,falcON_VERSION_D);
#else
      __origin[0] = 0;
#endif
      snprintf(__version,100,
#ifdef falcON_VERSION
	       falcON_VERSION
#endif
	       falcON_PSIFLAG);
      snprintf(__time,30,__DATE__ ", " __TIME__);
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

extern string defv[];                              // MUST be supplied by user  
extern string usage;	                           // MUST be supplied by user  
namespace nbdy { namespace defv_info {
  char version [100];                              //   "VERSION=..."           
  char compiled[100];                              //   "COMPILED=..."          
  void init() {
    compile_info::init();
    run_info::init();
    snprintf(version,100,"VERSION=%s%s\n%s",
	     compile_info::version (),
	     compile_info::compiler(),
	     compile_info::origin  ());
    snprintf(compiled,100,"COMPILED=\n%s, with %s",
	     compile_info::time(),
	     compile_info::compiler());
  };
} }
#  define falcON_DEFV nbdy::defv_info::version, nbdy::defv_info::compiled
#endif

//------------------------------------------------------------------------------
// define main() in namespace nbdy                                              
//------------------------------------------------------------------------------

namespace nbdy {
#ifdef falcON_USE_NEMO
  extern void main(void);                          // MUST be supplied by user  
#else
  extern void main(int argc, char *argv[]);        // MUST be supplied by user  
#endif
}

//------------------------------------------------------------------------------
// define global main(), which calls nemo::main()                               
//------------------------------------------------------------------------------

int main(int argc, char *argv[])                   // global main               
{
#ifdef falcON_USE_MPI
  nbdy::mpi_init(&argc,&argv);                     // start MPI: spawm processes
#endif

  nbdy::set_name(argv[0]);                         // get name of application   
  nbdy::compile_info::init();                      // initialize compile_info   
  nbdy::run_info::init();                          // initialize run_info       

#ifdef falcON_USE_NEMO
  nbdy::defv_info::init();                         // initialize defv_info      
  initparam(argv,defv);                            // start  NEMO               
#  if defined(falcON_RepAction) && (falcON_RepAction==1)
  nbdy::report::open_file(argv[0],*(ask_history()));
#  endif

  nbdy::main();                                    //   call user program       
  finiparam();                                     // finish NEMO               

#else // falcON_USE_NEMO

#  if defined(falcON_RepAction) && (falcON_RepAction==1)
  nbdy::report::open_file(argv[0]);
#  endif

  nbdy::main(argc,argv);                           //   call user program       

#endif

#if defined(falcON_RepAction) && (falcON_RepAction==1)
  nbdy::report::close_file();
#endif

#ifdef falcON_USE_MPI
  nbdy::mpi_finalize();                            // finish MPI                
#endif
}
//------------------------------------------------------------------------------
// define getvparam_z(arg,vect&) for getting a vect                             
//------------------------------------------------------------------------------
#ifdef falcON_USE_NEMO
#undef nemoinpr
#ifdef falcON_REAL_IS_FLOAT
#  define nemoinpr nemoinpf
#else
#  define nemoinpr nemoinpd
#endif
namespace nbdy {                                   // define in namespace nbdy  
  //----------------------------------------------------------------------------
  // read vect from nemo parameter                                              
  vect getvparam(char* option)
  {
    vect X;
    int  N = nemoinpr(getparam(option),static_cast<real*>(X),Ndim);
    if(N>Ndim)
      warning("option \"%s\" requires %d values, but %d given\n",option,Ndim,N);
    else if(N<0)  error("processing option \"%s\"\n",option);
    else if(N==0) error("no data for option \"%s\"\n",option);
    else if(N<Ndim)
      error  ("option \"%s\" requires %d values, but %d given\n",option,Ndim,N);
    return X;
  }
  //----------------------------------------------------------------------------
  // read vect from nemo parameter into 2nd argument, return pointer            
  vect* getvparam_z(char* option, vect&X)
  {
    if(!hasvalue(option)) return 0;
    int N = nemoinpr(getparam(option),static_cast<real*>(X),Ndim);
    if(Ndim != N) {
      if(N<0) error("processing option \"%s\"\n",option);
      if(N>0) warning("option \"%s\" requires %d values, but %d given\n",
		      option,Ndim,N);
      return 0;
    }
    return &X;
  }
//------------------------------------------------------------------------------
// define getparam_z(arg), and related, which will return 0 if !hasvalue(arg).  
//------------------------------------------------------------------------------
#ifdef getrparam
#  undef getrparam
#endif
  //----------------------------------------------------------------------------
  inline float getfparam(const char* option) {
    return float(getdparam(const_cast<char*>(option)));  }
  //----------------------------------------------------------------------------
  inline float getrparam(const char* option) { 
    return real(getdparam(const_cast<char*>(option))); }
  //----------------------------------------------------------------------------
  inline io getioparam(const char* option) {
    return io(getparam(const_cast<char*>(option))); }
//------------------------------------------------------------------------------
#  define GET_Z(TYPE,NAME,NULL)				\
  inline TYPE NAME##_z(char* option)			\
  { return hasvalue(option)? NAME(option) : NULL; }
  GET_Z(char*,    getparam, 0)
  GET_Z(int,      getiparam,0)
  GET_Z(long,     getlparam,0)
  GET_Z(bool,     getbparam,false)
  GET_Z(float,    getfparam,0.f)
  GET_Z(double,   getdparam,0.)
  GET_Z(real,     getrparam,zero)
  GET_Z(io,       getioparam, io::o)
#  undef GET_Z
  //////////////////////////////////////////////////////////////////////////////
}                                                  // END: namespace nbdy       
#endif                                             // falcON_USE_NEMO           

////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_main_h    

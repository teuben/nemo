// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// main.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
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
// In case, this is a NEMO application                                          
//                                                                              
// - set macro "falcON_USE_NEMO"                                                
// - include some NEMO stuff                                                    
// - make some extern declarations (to be satisfied by the user)                
//------------------------------------------------------------------------------

#if (defined(falcON_NEMO) && !defined(falcON_NONEMO) ) 
#  ifndef falcON_USE_NEMO
#    define falcON_USE_NEMO                        //   then use it             
#  endif
#  include <nemo.h>                                // NEMO basic stuff          

extern string defv[];                              // MUST be supplied by user  
extern string usage;	                           // MUST be supplied by user  
#endif

//------------------------------------------------------------------------------
// include other necessary header files                                         
//------------------------------------------------------------------------------

#include <public/exit.h>                           // falcON exit & error       

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
  nbdy::set_name(argv[0]);

#ifdef falcON_USE_NEMO
  initparam(argv,defv);                            //   start  NEMO             
#  if defined(falcON_RepAction) && (falcON_RepAction==1)
  nbdy::report::open_file(argv[0],*(ask_history()));
#  endif
  nbdy::main();                                    //     call user program     
  finiparam();                                     //   finish NEMO             

#else

#  if defined(falcON_RepAction) && (falcON_RepAction==1)
  nbdy::report::open_file(argv[0]);
#  endif
  nbdy::main(argc,argv);                           //     call user program     

#endif

#if defined(falcON_RepAction) && (falcON_RepAction==1)
  nbdy::report::close_file();
#endif
#ifdef falcON_USE_MPI
  nbdy::mpi_finalize();                            // finish MPI                
#endif
}
//------------------------------------------------------------------------------
// define getparam_z(arg), and related, which will return 0 if !hasvalue(arg).  
//------------------------------------------------------------------------------
#ifdef falcON_USE_NEMO
#  define GET_Z(TYPE,NAME) inline TYPE NAME##_z(char* option)	       	\
                         { return hasvalue(option)? NAME(option) : 0; }
namespace nbdy {                                   // define in namespace nbdy  
  GET_Z(char*,  getparam)
  GET_Z(int,    getiparam)
  GET_Z(long,   getlparam)
  GET_Z(bool,   getbparam)
  GET_Z(double, getdparam)
}
#  undef GET_Z
#endif                                             // falcON_included_main_h    
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_main_h    

// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// nemomain.h                                                                  |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002                                               |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// based on nemomain.c by Peter Teuben                                         |
//                                                                             |
// to be included by nemo_main only                                            |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// to run your code with MPI, you must                                         |
// have ALLOW_MPI  #defined                                                    |
// have PARALLEL   #defined                                                    |
//                                                                             |
// ALLOW_MPI    shall be #defined (via a compiler option), if MPI is installed |
// PARALLEL     shall be #defined in the including .cc file, if the code is    |
//              supposed to be a parallel one. Otherwise it will become an     |
//              ordinary serial code.                                          |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifdef included_nemomain_h
# error "nemomain.h must only be included by a nemo::main"
#else
#define included_nemomain_h 1

#if ! (defined(__cplusplus) || defined(c_plusplus) )
# error "nemomain.h must be used as C++ main only"
#endif

#if (defined(ALLOW_MPI) && defined(PARALLEL) )     // have MPI and want it?     
# ifndef USE_MPI
#  define USE_MPI                                  //   then use it             
# endif
#else                                              // else                      
# undef  USE_MPI                                   //   don't use it            
#endif

#ifdef ALLOW_NEMO                                  // this is a NEMO application
#ifdef USE_MPI                                     // use MPI                   
#  ifndef included_mpi_h
#    include <mpi.h>
#    define included_mpi_h
#  endif
#endif
// include some NEMO stuff                                                      
// (we cannot wrap it into namespace nemo here, 'cause they may #include others)
#include <stdinc.h>                                // NEMO basic stuff          
#include <getparam.h>                              // NEMO paramater input      
////////////////////////////////////////////////////////////////////////////////
extern string defv[];                              // MUST be supplied by user  
extern string usage;	                           // MUST be supplied by user  
namespace nemo {                                   // use namespace nemo        
  extern void main(void);                          // MUST be supplied by user  
}                                                  // CLOSE: namespace nemo     
//------------------------------------------------------------------------------
// define global main(), which calls nemo::main()                               
//------------------------------------------------------------------------------
int main(int argc,char *argv[])                    // global main               
{
#ifdef USE_MPI
  MPI_Init(&argc,&argv);                           // spawn MPI processes       
#endif
  initparam(argv,defv);                            //   start  NEMO             
  nemo::main();                                    //     call user program     
  finiparam();                                     //   finish NEMO             
#ifdef USE_MPI
  MPI_Finalize();                                  // finish MPI                
#endif
}
//------------------------------------------------------------------------------
// define getparam_z(arg), and related, which will return 0 if !hasvalue(arg).  
//------------------------------------------------------------------------------
#define GET_Z(TYPE,NAME) inline TYPE NAME##_z(char* option)	       	\
                         { return hasvalue(option)? NAME(option) : 0; }
namespace nemo {                                   // define in namespace nemo  
  GET_Z(char*,  getparam)
  GET_Z(int,    getiparam)
  GET_Z(long,   getlparam)
  GET_Z(bool,   getbparam)
  GET_Z(double, getdparam)
}
#undef GET_Z
////////////////////////////////////////////////////////////////////////////////
#endif                                             // ALLOW_NEMO                
#endif                                             // included_nemomain_h       

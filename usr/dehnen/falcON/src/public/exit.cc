//-----------------------------------------------------------------------------+
//                                                                             |
// exit.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002                                               |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/exit.h>

#ifdef ALLOW_MPI                                   // compiler option           
#  define  MPICH_SKIP_MPICXX
# include <mpi.h>                                  // C implement of MPI        
#endif

#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <cstring>

//------------------------------------------------------------------------------
void nbdy::exit(const int signal) {                // I: error signal           
#ifdef ALLOW_MPI                                   // IF MPI exists            >
  register int MPI_running;
  MPI_Initialized(&MPI_running);                   //   does MPI run?           
  if(MPI_running)                                  //   IF MPI runs            >
    MPI_Abort(MPI_COMM_WORLD,signal);              //     MPI-abort            <
  else                                             //   < ELSE(no MPI running) >
#endif                                             // < OR no MPI existing      
    std::exit(signal);                             //   ordinary exit()         
}
//------------------------------------------------------------------------------
void nbdy::error(const char* fmt,                  // I: error message          
                 ...             ) {               //[I: parameters]            
  va_list  ap;
  va_start(ap,fmt);
  fprintf(stderr,"@@@ Fatal error: ");
#ifdef ALLOW_MPI
  register int MPI_running;
  MPI_Initialized(&MPI_running);
  if(MPI_running) {
    register int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    fprintf(stderr,"on node %d",rank);
  }
#endif // ALLOW_MPI
  vfprintf(stderr,fmt,ap);
  if (fmt[strlen(fmt)-1] != '\n')
    fprintf(stderr,"\n");
  fflush(stderr);
  va_end(ap);
#ifdef ALLOW_MPI
  if(MPI_running)
    MPI_Abort(MPI_COMM_WORLD,1);                   //     MPI-abort            <
  else
#endif // ALLOW_MPI
  std::exit(1);
}
//------------------------------------------------------------------------------
void nbdy::warning(const char* fmt,                // I: error message          
                   ...             ) {             //[I: parameters]            
  va_list  ap;
  va_start(ap,fmt);
  fprintf(stderr,"@@@ Warning: ");
#ifdef ALLOW_MPI
  register int MPI_running;
  MPI_Initialized(&MPI_running);
  if(MPI_running) {
    register int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    fprintf(stderr,"on node %d",rank);
  }
#endif // ALLOW_MPI
  vfprintf(stderr,fmt,ap);
  if (fmt[strlen(fmt)-1] != '\n')
    fprintf(stderr,"\n");
  fflush(stderr);
  va_end(ap);
}
//------------------------------------------------------------------------------

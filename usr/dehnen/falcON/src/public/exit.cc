//-----------------------------------------------------------------------------+
//                                                                             |
// exit.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/exit.h>

#ifdef falcON_MPI                                  // compiler option           
#  define  MPICH_SKIP_MPICXX
#  include <mpi.h>                                 // C implement of MPI        
#endif

#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <cstring>
static char main_name[200] = {0};
//------------------------------------------------------------------------------
void nbdy::set_name(const char* name)
{
  strncpy(main_name,name,200);
}
//------------------------------------------------------------------------------
void nbdy::exit(const int signal) {                // I: error signal           
#ifdef falcON_MPI                                  // IF MPI exists            >
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
  if(main_name[0]) fprintf(stderr,"### falcON Fatal error [%s] : ",main_name);
  else             fprintf(stderr,"### falcON Fatal error : ");
  va_list  ap;
  va_start(ap,fmt);
#ifdef falcON_MPI
  register int MPI_running;
  MPI_Initialized(&MPI_running);
  if(MPI_running) {
    register int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    fprintf(stderr,"on node %d",rank);
  }
#endif // falcON_MPI
  vfprintf(stderr,fmt,ap);
  if (fmt[strlen(fmt)-1] != '\n')
    fprintf(stderr,"\n");
  fflush(stderr);
  va_end(ap);
#ifdef falcON_MPI
  if(MPI_running)
    MPI_Abort(MPI_COMM_WORLD,1);                   //     MPI-abort            <
  else
#endif // falcON_MPI
  std::exit(1);
}
//------------------------------------------------------------------------------
void nbdy::warning(const char* fmt,                // I: error message          
                   ...             ) {             //[I: parameters]            
  if(main_name[0]) fprintf(stderr,"### falcON Warning [%s] : ",main_name);
  else             fprintf(stderr,"### falcON Warning : ");
  va_list  ap;
  va_start(ap,fmt);
#ifdef falcON_MPI
  register int MPI_running;
  MPI_Initialized(&MPI_running);
  if(MPI_running) {
    register int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    fprintf(stderr,"on node %d",rank);
  }
#endif // falcON_MPI
  vfprintf(stderr,fmt,ap);
  if (fmt[strlen(fmt)-1] != '\n')
    fprintf(stderr,"\n");
  fflush(stderr);
  va_end(ap);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// support for reporting current actions to a temporary file                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_RepAction
#include <ctime>
#include <unistd.h>
//------------------------------------------------------------------------------
static struct report_data {
  char    FNAME[200];
  FILE   *STREAM;
  int     LEVEL;
  fpos_t  START;
  clock_t LAST;
  double  SECONDS;
} *REPORT = 0;
//------------------------------------------------------------------------------
const char* nbdy::report::file_name()
{
  return REPORT? REPORT->FNAME : 0;
}
//------------------------------------------------------------------------------
void nbdy::report::open_file(const char*prog, const char*com)
{
  if(REPORT) delete REPORT;
  REPORT = new report_data();
  REPORT->LEVEL   = 0;
  REPORT->LAST    = clock();
  REPORT->SECONDS = 0.;
  if(prog) snprintf(REPORT->FNAME,200,"%s.%d",prog,getpid());
  else     snprintf(REPORT->FNAME,200,"%d",getpid());
  REPORT->STREAM = fopen(REPORT->FNAME,"w");
  if(REPORT->STREAM) {
    time_t now = ::time(0);
    fprintf(REPORT->STREAM,
	    "\nfile \"%s\"\n"
	    "created by \"%s\"\n"
	    "        at %s"
	    "to report actions.\n"
	    "File will be deleted on succesful execution, "
	    "but not in case of abort.\n"
	    "This allows to narrow down the code position of mysterious "
	    "aborts\n\n",
	    REPORT->FNAME, com? com:prog, ctime(&now));
    fflush(REPORT->STREAM);
    fgetpos(REPORT->STREAM,&(REPORT->START));
  }
}
//------------------------------------------------------------------------------
void nbdy::report::close_file()
{
  if(REPORT) {
    if(REPORT->STREAM) fclose(REPORT->STREAM);
    if(unlink(REPORT->FNAME) != 0)
      warning("cannot delete file %s",REPORT->FNAME);
    delete REPORT;
    REPORT = 0;
  }
}
//------------------------------------------------------------------------------
nbdy::report::report(const char* fmt, ... )
{
  if(REPORT && REPORT->STREAM) {
    register clock_t NOW = clock();
    REPORT->SECONDS += (NOW - REPORT->LAST)/double(CLOCKS_PER_SEC);
    REPORT->LAST = NOW;
    fprintf(REPORT->STREAM,"%10.2f: ",REPORT->SECONDS);
    for(register int l=0; l!=REPORT->LEVEL; ++l) fprintf(REPORT->STREAM,"  ");
    fprintf(REPORT->STREAM,"begin: ");
    va_list  ap;
    va_start(ap,fmt);
    vfprintf(REPORT->STREAM,fmt,ap);
    va_end(ap);
    if (fmt[strlen(fmt)-1] != '\n')
      fprintf(REPORT->STREAM,"\n");
    fflush(REPORT->STREAM);
    fpos_t here;
    fgetpos(REPORT->STREAM,&here);
    fprintf(REPORT->STREAM,"<==\n");
    fflush (REPORT->STREAM);
    fsetpos(REPORT->STREAM,&here);
    ++(REPORT->LEVEL);
  }
}
//------------------------------------------------------------------------------
void nbdy::report::info(const char* fmt, ... )
{
  if(REPORT && REPORT->STREAM) {
    register clock_t NOW = clock();
    REPORT->SECONDS += (NOW - REPORT->LAST)/double(CLOCKS_PER_SEC);
    REPORT->LAST = NOW;
    fprintf(REPORT->STREAM,"%10.2f: ",REPORT->SECONDS);
    for(register int l=0; l!=REPORT->LEVEL; ++l) fprintf(REPORT->STREAM,"  ");
    fprintf(REPORT->STREAM,"doing: ");
    va_list  ap;
    va_start(ap,fmt);
    vfprintf(REPORT->STREAM,fmt,ap);
    va_end(ap);
    if (fmt[strlen(fmt)-1] != '\n')
      fprintf(REPORT->STREAM,"\n");
    fflush(REPORT->STREAM);
    fpos_t here;
    fgetpos(REPORT->STREAM,&here);
    fprintf(REPORT->STREAM,"<==\n");
    fflush (REPORT->STREAM);
    fsetpos(REPORT->STREAM,&here);
  }
}
//------------------------------------------------------------------------------
nbdy::report::~report()
{
  if(REPORT && REPORT->STREAM) {
    register clock_t NOW = clock();
    REPORT->SECONDS += (NOW - REPORT->LAST)/double(CLOCKS_PER_SEC);
    REPORT->LAST = NOW;
    --(REPORT->LEVEL);
    fprintf(REPORT->STREAM,"%10.2f: ",REPORT->SECONDS);
    for(register int l=0; l!=REPORT->LEVEL; ++l) fprintf(REPORT->STREAM,"  ");
    fprintf(REPORT->STREAM,"done\n");
    fflush(REPORT->STREAM);
    fpos_t here;
    fgetpos(REPORT->STREAM,&here);
    fprintf(REPORT->STREAM,"<==\n");
    fflush (REPORT->STREAM);
    if(REPORT->LEVEL == 0) fsetpos(REPORT->STREAM,&(REPORT->START));
    else                   fsetpos(REPORT->STREAM,&here);
  }
}
////////////////////////////////////////////////////////////////////////////////
#endif

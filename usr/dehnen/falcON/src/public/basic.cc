// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    src/public/basic.cc                                                
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2002-2006                                                          
///                                                                             
/// \brief   implements methods declared in inc/public/basic.h                  
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2002-2006  Walter Dehnen                                       
//                                                                              
// This program is free software; you can redistribute it and/or modify         
// it under the terms of the GNU General Public License as published by         
// the Free Software Foundation; either version 2 of the License, or (at        
// your option) any later version.                                              
//                                                                              
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//                                                                              
// You should have received a copy of the GNU General Public License            
// along with this program; if not, write to the Free Software                  
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#include <public/basic.h>

#ifdef falcON_MPI
#  include <proper/mpi_falcON.h>                   // C implement of MPI        
#endif

#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <cstring>
#include <ctime>
#ifdef unix
extern "C" {
#  include <unistd.h>
#  include <pwd.h>
}
#endif

////////////////////////////////////////////////////////////////////////////////
namespace {
  char lib_name  [300] = {0};
}
//------------------------------------------------------------------------------
const char* falcON::libdir()
{
  if(0 == *lib_name) {
    const char* d = directory();
    if(d) {
      strcpy(lib_name, d);
      strcat(lib_name,"/");
      strcat(lib_name,getenv("MACHTYPE"));
      strcat(lib_name,"_");
      strcat(lib_name,getenv("OSTYPE"));
    }
  }
  return lib_name;
}
//------------------------------------------------------------------------------
void falcON::error(const char* fmt,                // I: error message          
		   ...             ) {             //[I: parameters]            
  va_list  ap;
  va_start(ap,fmt);
  WDutils::printerr("### falcON Error: ", fmt, ap);
  va_end(ap);
  WDutils::exit();
}
//------------------------------------------------------------------------------
void falcON::warning(const char* fmt,              // I: warning message        
		     ...             ) {           //[I: parameters]            
  va_list  ap;
  va_start(ap,fmt);
  WDutils::printerr("### falcON Warning: ", fmt, ap);
  va_end(ap);
}
//------------------------------------------------------------------------------
void falcON::debug_info(int         deb,           // I: level for reporting    
			const char* fmt,           // I: debugging information  
			...             ) {        //[I: parameters]            
  if(RunInfo::debug(deb)) {
    va_list  ap;
    va_start(ap,fmt);
    printerr("### falcON Debug Info: ", fmt, ap, false);
    va_end(ap);
  }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// mechanism to avoid status mismatch                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
falcON::Status falcON::LibraryStatus()
{
  return CurrentStatus();
}
//------------------------------------------------------------------------------
void falcON::CheckAgainstLibrary(falcON::Status Current) falcON_THROWING
{
  Status Library = CurrentStatus();
  if( Current != Library ) {
    // check proprietary versus public
    if( Current&proper_version && !(Library&proper_version) )
      falcON_THROW("STATUS mismatch: proprietary executable, "
		   "but public-version library.\n");
    else if( Library&proper_version && !(Current&proper_version) )
      falcON_THROW("STATUS mismatch: public-version executable, "
		   "but proprietary library.\n");
    // check nemo
    if( Library&nemo_version && !(Current&nemo_version) )
      falcON_THROW("STATUS mismatch: "
		   "executable was not compiled with NEMO, "
		   "but library was.\n");
    else if( Current&nemo_version && !(Library&nemo_version) )
      falcON_THROW("STATUS mismatch: "
		   "executable was compiled with NEMO, "
		   "but library was not.\n");
    // check SPH
    if( Current&sph_version && !(Library&sph_version) )
      falcON_THROW("STATUS mismatch: "
		   "executable was compiled for SPH, "
		   "but library was not.\n");
    else if( Library&sph_version && !(Current&sph_version) )
      falcON_THROW("STATUS mismatch: "
		   "executable was not compiled for SPH, "
		   "but library was.\n");
    // check MPI
    if( Current&mpi_version && !(Library&mpi_version) )
      falcON_THROW("STATUS mismatch: "
		   "executable was compiled for MPI, "
		   "but library was not.\n");
    else if( Library&mpi_version && !(Current&mpi_version) )
      falcON_THROW("STATUS mismatch: "
		   "executable was not compiled for MPI, "
		   "but library was.\n");
    // check real
    if( Current&real_is_double && !(Library&real_is_double) )
      falcON_THROW("STATUS mismatch: "
		   "executable was compiled with real=double, "
		   "but library with real=float.\n");
    else if( Library&real_is_double && !(Current&real_is_double) )
      falcON_THROW("STATUS mismatch: "
		   "executable was compiled with real=float, "
		   "but library with real=double.\n");
    falcON_THROW("STATUS mismatch "
		 "between executable and library\n");
  }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// support for reporting current actions to a temporary file                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_RepAction
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
const char* falcON::report::file_name()
{
  return REPORT? REPORT->FNAME : 0;
}
//------------------------------------------------------------------------------
void falcON::report::open_file(const char*prog, const char*com)
{
  if(REPORT) delete REPORT;
  REPORT = new report_data();
  REPORT->LEVEL   = 0;
  REPORT->LAST    = clock();
  REPORT->SECONDS = 0.;
  if(prog) snprintf(REPORT->FNAME,200,"/tmp/%s.%d",prog,getpid());
  else     snprintf(REPORT->FNAME,200,"/tmp/%d",getpid());
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
void falcON::report::close_file()
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
falcON::report::report(const char* fmt, ... )
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
void falcON::report::info(const char* fmt, ... )
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
falcON::report::~report()
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
#endif // falcON_RepAction
////////////////////////////////////////////////////////////////////////////////

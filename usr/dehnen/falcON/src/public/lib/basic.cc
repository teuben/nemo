// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    src/public/basic.cc                                                
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2002-2008                                                          
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
  if(0 == *lib_name)
    strcpy(lib_name,getenv("FALCONLIB"));
  return lib_name;
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
void falcON::CheckAgainstLibrary(falcON::Status Current,
				 const char    *Program) falcON_THROWING
{
  Status Library = CurrentStatus();
  if( Current != Library ) {
    DebugInfo(5,"CheckAgainstLibrary(): Current=%d Library=%d\n",
	      Current, Library);
    // check proprietary versus public
    if( Current&proper_version && !(Library&proper_version) )
      falcON_THROW("STATUS mismatch: proprietary %s, "
		   "but public-version library.\n",Program);
    else if( Library&proper_version && !(Current&proper_version) )
      falcON_THROW("STATUS mismatch: public-version %s, "
		   "but proprietary library.\n",Program);
    // check nemo
    if( Library&nemo_version && !(Current&nemo_version) )
      falcON_THROW("STATUS mismatch: "
		   "%s was not compiled with NEMO, "
		   "but library was.\n",Program);
    else if( Current&nemo_version && !(Library&nemo_version) )
      falcON_THROW("STATUS mismatch: "
		   "%s was compiled with NEMO, "
		   "but library was not.\n",Program);
    // check SPH
    if( Current&sph_version && !(Library&sph_version) )
      falcON_THROW("STATUS mismatch: "
		   "%s was compiled for SPH, "
		   "but library was not.\n",Program);
    else if( Library&sph_version && !(Current&sph_version) )
      falcON_THROW("STATUS mismatch: "
		   "%s was not compiled for SPH, "
		   "but library was.\n",Program);
    // check real
    if( Current&real_is_double && !(Library&real_is_double) )
      falcON_THROW("STATUS mismatch: "
		   "%s was compiled with real=double, "
		   "but library with real=float.\n",Program);
    else if( Library&real_is_double && !(Current&real_is_double) )
      falcON_THROW("STATUS mismatch: "
		   "%s was compiled with real=float, "
		   "but library with real=double.\n",Program);
    falcON_THROW("STATUS mismatch "
		 "between %s and library\n",Program);
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
  if(REPORT) falcON_DEL_O(REPORT);
  REPORT = new report_data();
  REPORT->LEVEL   = 0;
  REPORT->LAST    = clock();
  REPORT->SECONDS = 0.;
  if(prog) SNprintf(REPORT->FNAME,200,"/tmp/%s.%d",prog,getpid());
  else     SNprintf(REPORT->FNAME,200,"/tmp/%d",getpid());
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
      falcON_Warning("cannot delete file %s",REPORT->FNAME);
    falcON_DEL_O(REPORT);
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

// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    src/exception.cc                                                   
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2000-2006                                                          
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2000-2006  Walter Dehnen                                       
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
#include <exception.h>
#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <cstring>
#include <ctime>
#include <iostream>
#include <fstream>
#ifdef unix
extern "C" {
#  include <unistd.h>
#  include <pwd.h>
}
#endif
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class RunInfo                                                                
//                                                                              
////////////////////////////////////////////////////////////////////////////////
WDutils::RunInfo::RunInfo()
  : __host_known(0),
    __user_known(0),
    __pid_known(0),
    __name_known(0),
    __is_mpi_proc(0),
    __debug(0)
{
  try {
  // set time
    time_t now = ::time(0);
    SNprintf(__time,100,ctime(&now));
    __time[24] = 0;
    // set host name, user name, and pid
#ifdef unix
    gethostname(__host,100);
    SNprintf(__user,100,(getpwuid(geteuid())->pw_name));
    SNprintf(__pid,20,"%d",getpid());
    __host_known = 1;
    __user_known = 1;
    __pid_known  = 1;
    char file[64];
    SNprintf(file,64,"/proc/%s/cmdline",__pid);
    std::ifstream in(file);
    if(file) {
      int i,e=0;
      for(i=0; i!=1024; ++i) __cmd[i]=0;
      in.getline(__cmd,1023);
      for(i=1023; i!=0; --i)
	if(__cmd[i]==0 || isspace(__cmd[i])) __cmd[i] = ' ';
	else if(e==0) e=i;
      __cmd[e+1] = 0;
      for(i=0; !isspace(__cmd[i]); ++i)
	__name[i] = __cmd[i];
      __name[i] = 0;
      __name[i] = 0;
      __cmd_known  = 1;
      __name_known = 1;
    }
#else
    SNprintf(__host,100,"unknown.host");
    SNprintf(__user,100,"unknown.user");
    SNprintf(__user,100,"unknown.main");
#endif
  } catch(exception E) {
    WDutils_RETHROW(E);
  }
}
WDutils::RunInfo WDutils::RunInfo::Info;
////////////////////////////////////////////////////////////////////////////////
WDutils::exiter WDutils::exit = &std::exit;
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace WDutils;
  inline void printerr(const char*header,
		       const char*fmt,
		       va_list   &ap,
		       const char*file = 0,
		       int        line = 0,
		       bool       name = true)
  {
    fprintf(stderr,header);
    if(name && RunInfo::name_known())
      fprintf(stderr,"[%s]: ",RunInfo::name());
    if(RunInfo::is_mpi_proc())
      fprintf(stderr,"@%d: ",RunInfo::mpi_proc());
    if(file) fprintf(stderr,"[%s:%d]: ",file,line);
    vfprintf(stderr,fmt,ap);
    if (fmt[strlen(fmt)-1] != '\n')
      fprintf(stderr,"\n");
    fflush(stderr);
  }
}
//------------------------------------------------------------------------------
void WDutils::Error::operator()(const char* fmt, ...) const
{
  char header[35];
  SNprintf(header,35,"### %s Error: ",lib);
  va_list  ap;
  va_start(ap,fmt);
  printerr(header, fmt, ap, file, line);
  va_end(ap);
  WDutils::exit(1);
}
//------------------------------------------------------------------------------
void WDutils::Warning::operator()(const char* fmt, ...) const
{
  char header[37];
  SNprintf(header,37,"### %s Warning: ",lib);
  va_list  ap;
  va_start(ap,fmt);
  printerr(header, fmt, ap, file, line);
  va_end(ap);
}
//------------------------------------------------------------------------------
void WDutils::__DebugInfo::operator()(const char* fmt, ...) const
{
  char header[40];
  SNprintf(header,40,"### %s Debug Info: ",lib);
  va_list  ap;
  va_start(ap,fmt);
  printerr(header, fmt, ap, file, line, false);
  va_end(ap);
}
//------------------------------------------------------------------------------
void WDutils::__DebugInfo::operator()(int deb, const char* fmt, ...) const
{
  if(RunInfo::debug(deb)) {
    char header[40];
    SNprintf(header,40,"### %s Debug Info: ",lib);
    va_list  ap;
    va_start(ap,fmt);
    printerr(header, fmt, ap, file, line, false);
    va_end(ap);
  }
}
//------------------------------------------------------------------------------
WDutils::exception WDutils::Thrower::operator()(const char*fmt, ...) const
{
  size_t size = 1024;
  char   buffer[1024], *buf=buffer;
  if(file) {
    int len = SNprintf(buf,size,"in [%s:%d]: ",file,line);
    buf  += len;
    size -= len;
  }
  va_list  ap;
  va_start(ap,fmt);
  vsnprintf(buf, size, fmt, ap);
  va_end(ap);
  return exception(buffer);
}
//------------------------------------------------------------------------------
WDutils::exception::exception(const char*fmt, ...)
{
  const int size=1024;
  char __text[size];
  va_list  ap;
  va_start(ap,fmt);
  int w = vsnprintf(__text,size,fmt,ap);
  if(w>=size) {
    WDutils_Warning("WDutils::exception::exception(): "
		    "string size of %d characters exceeded\n",size);
    __text[size-1]=0;
  }
  if(w<0)
    WDutils_Warning("WDutils::exception::exception(): formatting error\n");
  va_end(ap);
  std::string::operator= (__text);
}
//------------------------------------------------------------------------------
WDutils::message::message(const char*fmt, ...) throw(exception)
{
  va_list  ap;
  va_start(ap,fmt);
  int w = vsnprintf(__text,size,fmt,ap);
  if(w>=size) throw exception("WDutils::message::message(): "
			      "string size of %d characters exceeded\n",size);
  if(w <   0) throw exception("WDutils::message::message(): "
			      "formatting error\n");
  va_end(ap);
}
//------------------------------------------------------------------------------
int WDutils::snprintf(char*str, size_t l, const char* fmt, ...)
  WDutils_THROWING
{
  va_list  ap;
  va_start(ap,fmt);
  int w = std::vsnprintf(str,l,fmt,ap);
  va_end(ap);
  if(w==l) WDutils_THROW("snprintf(): trailing 0 lost");
  if(w >l) WDutils_THROW("snprintf(): string size exceeded [%d:%d]",w,l);
  if(w <0) WDutils_THROW("snprintf(): formatting error");
  return w;
}
//------------------------------------------------------------------------------
int WDutils::snprintf__::operator()(char*str, size_t l, const char* fmt, ...)
  WDutils_THROWING
{
  va_list  ap;
  va_start(ap,fmt);
  int w = std::vsnprintf(str,l,fmt,ap);
  va_end(ap);
  if(w==l) WDutils_THROW("[%s:%d]: snprintf(): trailing 0 lost",file,line);
  if(w >l) WDutils_THROW("[%s:%d]: snprintf(): string size exceeded [%d:%d]",
			 w,l,file,line);
  if(w <0) WDutils_THROW("[%s:%d]: snprintf(): formatting error",file,line);
  return w;
}
////////////////////////////////////////////////////////////////////////////////

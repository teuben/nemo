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
    __debug(0)
{
  // set time
  time_t now = ::time(0);
  snprintf(__time,20,ctime(&now));
  // set host name, user name, and pid
#ifdef unix
  gethostname(__host,100);
  snprintf(__user,100,(getpwuid(geteuid())->pw_name));
  snprintf(__pid,10,"%d",getpid());
  __host_known = 1;
  __user_known = 1;
  __pid_known  = 1;
  char file[64];
  snprintf(file,64,"/proc/%s/cmdline",__pid);
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
  snprintf(__host,100,"unknown.host");
  snprintf(__user,100,"unknown.user");
  snprintf(__user,100,"unknown.main");
#endif
}
WDutils::RunInfo WDutils::RunInfo::Info;
////////////////////////////////////////////////////////////////////////////////
void WDutils::exit(int signal) {                   // I: error signal           
#ifdef _MPI_INCLUDE                                // IF MPI exists             
  register int MPI_running;
  MPI_Initialized(&MPI_running);                   //   does MPI run?           
  if(MPI_running)                                  //   IF MPI runs             
    MPI_Abort(MPI_COMM_WORLD,signal);              //     MPI-abort             
  else                                             //   ELSE(no MPI running)    
#endif                                             // ELSE no MPI existing      
    std::exit(signal);                             //   ordinary exit()         
}
////////////////////////////////////////////////////////////////////////////////
template<typename VA_LIST>
void WDutils::printerr(const char*header,
		       const char*fmt,
		       VA_LIST   &ap,
		       bool       give_name)
{
  fprintf(stderr,header);
  if(give_name && RunInfo::name_known())
    fprintf(stderr,"[%s]: ",RunInfo::name());
#ifdef _MPI_INCLUDE
  register int MPI_running;
  MPI_Initialized(&MPI_running);
  if(MPI_running) {
    register int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    fprintf(stderr,"on node %d",rank);
  }
#endif
  vfprintf(stderr,fmt,ap);
  if (fmt[strlen(fmt)-1] != '\n')
    fprintf(stderr,"\n");
  fflush(stderr);
}
template 
void WDutils::printerr<va_list>(const char*,  const char*, va_list&, bool);
//------------------------------------------------------------------------------
void WDutils::error(const char* fmt,               // I: error message          
		    ...             )              //[I: parameters]            
{
  va_list  ap;
  va_start(ap,fmt);
  printerr("### WDutils Error: ", fmt, ap);
  va_end(ap);
  WDutils::exit();
}
//------------------------------------------------------------------------------
void WDutils::warning(const char* fmt,             // I: warning message        
		      ...             )            //[I: parameters]            
{
  va_list  ap;
  va_start(ap,fmt);
  printerr("### WDutils Warning: ", fmt, ap);
  va_end(ap);
}
//------------------------------------------------------------------------------
void WDutils::debug_info(const char* fmt,          // I: debugging information  
			 ...             )         //[I: parameters]            
{
  va_list  ap;
  va_start(ap,fmt);
  printerr("### WDutils Debug Info: ", fmt, ap, false);
  va_end(ap);
}
//------------------------------------------------------------------------------
void WDutils::debug_info(int         deb,          // I: level for reporting    
			 const char* fmt,          // I: debugging information  
			 ...             )         //[I: parameters]            
{
  if(RunInfo::debug(deb)) {
    va_list  ap;
    va_start(ap,fmt);
    printerr("### WDutils Debug Info: ", fmt, ap, false);
    va_end(ap);
  }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// WDutils::message                                                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
WDutils::message::message(const char*fmt, ...)
{
  va_list  ap;
  va_start(ap,fmt);
  vsnprintf(__text,size,fmt,ap);
  va_end(ap);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// WDutils::exception                                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
WDutils::exception::exception(const char*fmt, ...)
{
  char __text[1024];
  va_list  ap;
  va_start(ap,fmt);
  vsnprintf(__text,1024,fmt,ap);
  va_end(ap);
  std::string::operator= (__text);
}
////////////////////////////////////////////////////////////////////////////////

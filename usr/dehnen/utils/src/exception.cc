// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    src/exception.cc
///
/// \author  Walter Dehnen
///
/// \date    2000-14
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-14  Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
#include <exception.h>
#include <cerrno>
#include <cstdio>
#include <cstdarg>
#include <cstring>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#ifdef _OPENMP
#  include <omp.h>
#endif
extern "C" {
#if defined(__unix) || defined(__DARWIN_UNIX03)
#  include <unistd.h>
#  include <sys/time.h>
#endif
#ifdef WIN32
#  include <windows.h>
#endif
}

#ifdef __INTEL_COMPILER
#pragma warning (disable:981) /* operands are evaluated in unspecified order */
#pragma warning (disable:161) /* unrecognized pragma */
#endif

#if __cplusplus >= 201103L && defined(WDutilsDevel)
#  include <map>
#  include <thread>
#  include <mutex>
#endif
#ifdef WDutilsTBB
#  ifdef __clang__
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wundef"
#    pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#  endif
#  include <tbb/task_scheduler_init.h>
#  ifdef __clang__
#    pragma clang diagnostic pop
#  endif
namespace {
  typedef tbb::task_scheduler_init _tbb_ts_init;
}
#endif
//                                                                              
// class RunInfo                                                                
//                                                                              
WDutils::RunInfo::RunInfo()
  : _m_host_known(0)
  , _m_user_known(0)
  , _m_pid_known(0)
  , _m_name_known(0)
  , _m_is_mpi_proc(0)
  , _m_debug(0)
  , _m_tbb_init(0)
{
  try {
    // set wall-clock time
    {
#if defined(__unix) || defined(__DARWIN_UNIX03)
      timeval now;
      gettimeofday(&now, NULL);
      _m_sec = now.tv_sec;
      _m_usec = now.tv_usec;
#elif defined(WIN32)
      LARGE_INTEGER tmp;
      QueryPerformanceCounter(&tmp);
      QueryPerformanceFrequency(&tmp);
      _m_timecount = tmp.QuadPart;
      _m_timetick = 1.0/double(tmp.QuadPart);
#endif
    }
    // set run time
    {
      time_t now = ::time(0);
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wformat-security"
#endif
      SNprintf(_m_time,104,ctime(&now));
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
      _m_time[24] = 0;
    }
#if defined(__unix) || defined(__DARWIN_UNIX03)
    // set host name
    {
      gethostname(_m_host,104);
      _m_host_known = 1;
    }
    // set user name
    {
      const char*user_ = getenv("USER");
      if(user_) {
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wformat-security"
#endif
	SNprintf(_m_user,104,user_);
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
	_m_user_known = 1;
      } else
	SNprintf(_m_user,104,"unknown.user");
    }
    // set pid
    {
      _m_pid_num = getpid();
      SNprintf(_m_pid,24,"%d",_m_pid_num);
      _m_pid_known  = 1;
    }
    // set command and name of executable
    {
      char file[64];
      if(_m_pid_known) {
	SNprintf(file,64,"/proc/%s/cmdline",_m_pid);
	std::ifstream in(file);
	if(in) {
	  int i,e=0;
	  for(i=0; i!=1024; ++i) _m_cmd[i]=0;
	  in.getline(_m_cmd,1023);
	  for(i=1023; i!=0; --i)
	    if(_m_cmd[i]==0 || isspace(_m_cmd[i])) _m_cmd[i] = ' ';
	    else if(e==0) e=i;
	  _m_cmd[e+1] = 0;
	  for(i=0; !isspace(_m_cmd[i]); ++i)
	    _m_name[i] = _m_cmd[i];
	  _m_name[i] = 0;
	  _m_name[i] = 0;
	  _m_cmd_known  = 1;
	  _m_name_known = 1;
	}
      }
    }
#else // __unix
    {
      SNprintf(_m_host,104,"unknown.host");
      SNprintf(_m_user,104,"unknown.user");
      SNprintf(_m_name,104,"unknown.name");
    }
#endif
    // set # proc available for openMP
    {
#ifdef _OPENMP
      if(omp_in_parallel())
	WDutils_ErrorF("called inside OMP parallel region\n");
      _m_omp_proc = omp_get_num_procs();
#else
      _m_omp_proc = 1;
#endif
      _m_omp_size = _m_omp_proc;
    }
    // set # threads used by TBB
    {
#ifdef WDutilsTBB
      _m_tbb_proc = unsigned(_tbb_ts_init::default_num_threads());
#else
      _m_tbb_proc = 1u;
#endif 
      _m_tbb_size = _m_tbb_proc;
    }
  } 
  catch(WDutils::exception ex) { WDutils_RETHROW(ex); }
}
//
void WDutils::RunInfo::set_omp(int 
#ifdef _OPENMP
			       n
#endif
			       )
{
#ifdef _OPENMP
  Info._m_omp_size = n;
  if(Info._m_omp_size < 1) {
    Info._m_omp_size = 1;
    WDutils_WarningN("RunInfo::set_omp('%d') assume '1'\n",n);
  }
  omp_set_num_threads(Info._m_omp_size);
#else
  Info._m_omp_size = 1;
#endif
}
//
void WDutils::RunInfo::set_omp(const char*
#ifdef _OPENMP
			       arg
#endif
			       )
{
#ifdef _OPENMP
  if(arg==0 || arg[0]==0 || arg[0]=='t')
    Info._m_omp_size = Info._m_omp_proc;
  else if(arg[0] == 'f')
    Info._m_omp_size = 1;
  else if(arg && arg[0]) {
    Info._m_omp_size = strtol(arg,0,10);
    if(errno == EINVAL)
      WDutils_THROWN("RunInfo::set_omp('%s') (errno=EINVAL)\n",arg,errno);
    if(errno == ERANGE)
      WDutils_THROWN("RunInfo::set_omp('%s') (errno=ERANGE)\n",arg,errno);
    if(Info._m_omp_size < 1) {
      Info._m_omp_size = 1;
      WDutils_WarningN("RunInfo::set_omp('%s') assume '1'\n",arg);
    }
  }
  omp_set_num_threads(Info._m_omp_size);
#else
  Info._m_omp_size = 1;
#endif
}
//
WDutils::RunInfo::~RunInfo()
{
#ifdef WDutilsTBB
  if(Info._m_tbb_init)
    delete static_cast<_tbb_ts_init*>(Info._m_tbb_init);
#endif
  Info._m_tbb_init = 0;
}
//
#ifdef WDutilsTBB
void WDutils::RunInfo::set_tbb(unsigned n)
{
  if(Info._m_tbb_init==0)
    Info._m_tbb_init = n==0?
      new _tbb_ts_init(_tbb_ts_init::automatic) :
      new _tbb_ts_init(int(n))                  ;
  else if(n>=1 && n!=Info._m_tbb_size)
    static_cast<_tbb_ts_init*>(Info._m_tbb_init)
      -> initialize(int(n));
  else if(n==0 && Info._m_tbb_proc != Info._m_tbb_size)
    static_cast<_tbb_ts_init*>(Info._m_tbb_init)
      -> initialize(int(Info._m_tbb_proc));
  Info._m_tbb_size = n==0? Info._m_tbb_proc : n;
}
//
void WDutils::RunInfo::set_tbb(const char* arg)
{
  unsigned n=0;
  if(arg==0 || arg[0]==0 || arg[0]=='t')
    n = Info._m_tbb_proc;
  else if(arg[0] == 'f')
    n = 1u;
  else if(arg && arg[0]) {
    n = unsigned(strtol(arg,0,10));
    if(errno == EINVAL)
      WDutils_THROWN("RunInfo::set_tbb('%s') (errno=EINVAL)\n",arg,errno);
    if(errno == ERANGE)
      WDutils_THROWN("RunInfo::set_tbb('%s') (errno=ERANGE)\n",arg,errno);
  }
  set_tbb(n);
}
//
#elif defined(WDutilsDevel)
#  warning not implementing TBB stuff
#endif // WDutilsTBB
//
#if __cplusplus >= 201103L && defined(WDutilsDevel)
unsigned WDutils::RunInfo::thread_id()
{
#  ifdef __clang__
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wexit-time-destructors"
#  endif
  // yes, I know that ids will be destroyed at exit time. this is not problem.
  static unsigned nextindex = 0;
  static std::mutex my_mutex;
  static std::map<std::thread::id, unsigned> ids;
  auto id = std::this_thread::get_id();
  std::lock_guard<std::mutex> lock(my_mutex);
  if(ids.find(id) == ids.end())
    ids[id] = nextindex++;
  return ids[id];
#  ifdef __clang__
#    pragma clang diagnostic pop
#  endif
}
#endif
//
void WDutils::RunInfo::header(std::ostream&out)
{
  if(out) {
    if(RunInfo::cmd_known())
      out<<"# \""<<RunInfo::cmd()<<"\"\n#\n";
    out<<"# run at  "  <<RunInfo::time()<<"\n";
    if(RunInfo::user_known())  out<<"#     by  \""<<RunInfo::user()<<"\"\n";
    if(RunInfo::host_known())  out<<"#     on  \""<<RunInfo::host()<<"\"\n";
    if(RunInfo::pid_known())   out<<"#     pid  " <<RunInfo::pid() <<"\n";
    if(RunInfo::is_mpi_proc()) out<<"#     mpi  " <<RunInfo::mpi_size()<<"\n";
    out<<"#\n";
  }
}
//
#if defined(__unix) || defined(__DARWIN_UNIX03)
double WDutils::RunInfo::WallClock()
{
  timeval now;
  gettimeofday(&now, NULL);
  return (now.tv_sec - Info._m_sec) + (now.tv_usec - Info._m_usec) * 1.e-6;
}
void WDutils::RunInfo::WallClock(unsigned&sec, unsigned&usec)
{
  timeval now;
  gettimeofday(&now, NULL);
  if(now.tv_usec > Info._m_usec) {
    sec  = unsigned(now.tv_sec  - Info._m_sec);
    usec = unsigned(now.tv_usec - Info._m_usec);
  } else {
    sec  = unsigned(now.tv_sec  - Info._m_sec - 1);
    usec = unsigned((1000000 + now.tv_usec) - Info._m_usec);
  }
}
#elif defined(WIN32)
double WDutils::RunInfo::WallClock()
{
  LARGE_INTEGER now;
  QueryPerformanceCounter(&now);
  return (now.QuadPart - Info._m_timecount) * Info._m_timetick;
}
#endif
//
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wglobal-constructors"
#  pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif
WDutils::RunInfo WDutils::RunInfo::Info;
#ifdef __clang__
#  pragma clang diagnostic pop
#endif

//

namespace {
  using namespace WDutils;

#if defined(__PGI) || defined(__SUNPRO_CC)
  using ::snprintf;
  using ::vsnprintf;
#else
  using std::snprintf;
  using std::vsnprintf;
#endif

  void printerr(const char*lib,
		const char*issue,
		const char*fmt,
		va_list   &ap,
		int        depth,
		const char*func,
		const char*file,
		unsigned   line,
		unsigned   
#ifdef _OPENMP
		flag
#endif
		,
		bool       name)
  {
    const int size=1024;
    int s=size, w=0;
    char ffmt[size], *t=ffmt;
    char dpth[21] = "                    ";
    if(depth>20) depth=20;
    dpth[depth]=0;
    if(lib) {
      w=snprintf(t,size_t(s),"# %s %s",lib,issue);
      t+=w; s-=w;
    } else if(issue) {
      w=snprintf(t,size_t(s),"# %s",issue);
      t+=w; s-=w;
    }
    if(name && RunInfo::name_known()) {
      w=snprintf(t,size_t(s)," [%s]",RunInfo::name());
      t+=w; s-=w;
    }
    if(RunInfo::is_mpi_proc()) {
      w=snprintf(t,size_t(s)," @%2d",RunInfo::mpi_proc());
      t+=w; s-=w;
#ifdef _OPENMP
    } else if( (flag&1) && omp_in_parallel()) {
      w=snprintf(t,size_t(s)," @%2d",omp_get_thread_num());
      t+=w; s-=w;
#endif
    }
    if(file) {
      w=snprintf(t,size_t(s)," [%s:%d]",file,line);
      t+=w; s-=w;
    }
    if(func) {
      w=snprintf(t,size_t(s)," in %s",func);
      t+=w; s-=w;
    }      
    if (fmt[strlen(fmt)-1] != '\n')
      /*w=*/snprintf(t,size_t(s),": %s%s\n",dpth,fmt);
    else
      /*w=*/snprintf(t,size_t(s),": %s%s",dpth,fmt);
    //t+=w; s-=w;
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wformat-nonliteral"
#  ifdef __unix
#    pragma clang diagnostic ignored "-Wdisabled-macro-expansion"
#  endif
#endif
    vfprintf(stderr,ffmt,ap);
    fflush(stderr);
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
  }
}
//
#if __cplusplus >= 201103L
template<typename Traits>
void WDutils::Reporting<Traits>::print_header(std::ostringstream&ostr,
					      int depth) const
{
  char dpth[21] = "                    ";
  if(depth>20) depth=20;
  dpth[depth] = 0;
  ostr     << '#';
  if(library)
    ostr   << ' ' << library;
  if(Traits::issue())
    ostr   << ' ' << Traits::issue();
  bool mpi = RunInfo::is_mpi_proc();
  bool thr = flag&1;
  if(mpi || thr) {
    ostr   << " @";
    if(mpi)
      ostr << RunInfo::mpi_proc();
#if defined(WDutilsDevel)
    if(mpi && thr)
      ostr << ':';
    if(thr)
      ostr << RunInfo::thread_id();
#endif
  }
  if(depth)
    ostr   << dpth;
  if(file)
    ostr   << " [" << file << ':' << line << ']';
  if(func)
    ostr   << " in " << func;
  ostr     << ": ";
}
#endif
//
template<typename Traits>
void WDutils::Reporting<Traits>::operator()(int lev, const char* fmt, ...) const
{
  if(Traits::condition(lev)) {
    va_list  ap;
    va_start(ap,fmt);
    printerr(library, Traits::issue(), fmt, ap, lev, func, file, line, flag,
	     false);
    va_end(ap);
    Traits::after();
  }
}
//
template<typename Traits>
void WDutils::Reporting<Traits>::operator()(const char* fmt, ...) const
{
  va_list  ap;
  va_start(ap,fmt);
  printerr(library, Traits::issue(), fmt, ap, 0, func, file, line, flag, false);
  va_end(ap);
  Traits::after();
}
//
template struct WDutils::Reporting<WDutils::DebugInfoTraits>;
template struct WDutils::Reporting<WDutils::ErrorTraits>;
template struct WDutils::Reporting<WDutils::WarningTraits>;
//
WDutils::Thrower::handler WDutils::Thrower::InsteadOfThrow=0;
//
WDutils::exception WDutils::Thrower::operator()(const char*fmt, ...) const
{
  size_t size = 1024, len;
  char   buffer[1024], *buf=buffer;
  bool   error = false 
#ifdef _OPENMP
    || (omp_get_level() && InsteadOfThrow)
#endif
    ;
  if(!error && file) {
    len  = size_t(SNprintf(buf,size,"[%s:%d]",file,line));
    buf += len;
    size-= len;
  }
  if(func) {
    len  = size_t(file?
		  SNprintf(buf,size," in %s",func) :
		  SNprintf(buf,size, "in %s",func) ) ;
    buf += len;
    size-= len;
  }
  len  = size_t(SNprintf(buf,size,": "));
  buf += len;
  size-= len;
  va_list  ap;
  va_start(ap,fmt);
  vsnprintf(buf, size, fmt, ap);
  va_end(ap);
  if(error)
    InsteadOfThrow(file,line,buffer);
  return exception(buffer);
}
//
WDutils::exception::exception(const char*fmt, ...)
  : std::runtime_error(std::string())
{
  const int msize=1024;
  char _m_text[msize];
  va_list  ap;
  va_start(ap,fmt);
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wformat-nonliteral"
#endif
  int w = vsnprintf(_m_text,msize,fmt,ap);
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
  if(w>=msize) {
    WDutils_WarningF("string size of %d characters exceeded\n",msize);
    _m_text[msize-1]=0;
  }
  if(w<0)
    WDutils_WarningF("formatting error\n");
  va_end(ap);
  std::runtime_error::operator= (std::runtime_error(_m_text));
}
//
WDutils::message::message(const char*fmt, ...)
{
  va_list  ap;
  va_start(ap,fmt);
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wformat-nonliteral"
#endif
  int w = vsnprintf(_m_text,size,fmt,ap);
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
  if(w>=static_cast<int>(size))
    WDutils_THROW("string size of %ld characters exceeded\n",size);
  if(w < 0)
    WDutils_THROW("formatting error\n");
  va_end(ap);
}
//
int WDutils::snprintf(char*str, size_t l, const char* fmt, ...)
{
  va_list ap;
  va_start(ap,fmt);
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wformat-nonliteral"
#endif
  int w = vsnprintf(str,l,fmt,ap);
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
  va_end(ap);
  if(w==static_cast<int>(l))
    WDutils_THROWF("trailing 0 lost");
  if(w >static_cast<int>(l))
    WDutils_THROWF("string size exceeded [%d:%lu]",w,l);
  if(w <0)
    WDutils_THROWF("formatting error");
  return w;
}
//
int WDutils::snprintf__::operator()(char*str, size_t l, const char* fmt, ...)
{
  va_list ap;
  va_start(ap,fmt);
  int w = vsnprintf(str,l,fmt,ap);
  va_end(ap);
  if(w==static_cast<int>(l))
    WDutils_THROWER("snprintf()",file,line)
      ("trailing 0 lost");
  if(w >static_cast<int>(l))
    WDutils_THROWER("snprintf()",file,line)
      ("string size exceeded [%d:%lu]",w,l);
  if(w <0)
    WDutils_THROWER("snprintf()",file,line)
      ("formatting error");
  return w;
}
////////////////////////////////////////////////////////////////////////////////

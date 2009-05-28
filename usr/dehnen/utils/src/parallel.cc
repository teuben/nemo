// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    utils/src/parallel.cc
///
/// \author  Walter Dehnen
///
/// \date    2008,2009
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008,2009  Walter Dehnen
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
#include <parallel.h>
#include <io.h>
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <fstream>

// /////////////////////////////////////////////////////////////////////////////
// auxiliary stuff
namespace {
  using namespace WDutils;
#ifdef WDutils_EXCEPTIONS
  //----------------------------------------------------------------------------
  /// for generating exceptions from errors in MPI library functions
  struct Thrower {
    const char*file;          ///< file name
    const int  line;          ///< line number
    /// constructor: get file name, and line number
    Thrower(const char*__file, int __line) : file(__file), line(__line) {}
    /// generate an exception
    /// \param[in] err  error code as returned by MPI routines (C-binding)
    /// \param[in] fnc  (optional) name of calling function
    WDutils::exception operator()(int errcode, const char*func = 0) const
    {
      const int msize = 512+MPI_MAX_ERROR_STRING;
      int size = msize, len;
      char buffer[msize], *buf=buffer;
      if(file) {
	len   = SNprintf(buf,size,"in [%s:%d]: ",file,line);
	buf  += len;
	size -= len;
      }
      len = func?
	SNprintf(buf,size,"MPI Error in %s: ",func) :
	SNprintf(buf,size,"MPI Error: ") ;
      buf += len;
      MPI_Error_string(errcode, buf, &len);
      return WDutils::exception(buffer);
    }
  };
#  define MPI_THROWER  throw ::Thrower
#else
  //----------------------------------------------------------------------------
  /// to be used for reporting errors from the MPI library
  struct Error {
    const char*file;          ///< file name
    const int  line;          ///< line number
    /// default constructor: set data to NULL
    Error() : file(0), line(0) {}
    /// constructor: get file name, and line number
    Error(const char*__file, int __line) : file(__file), line(__line) {}
    /// print error message to stderr, report [file:line] if known.
    /// \param[in] err  error code as returned by MPI routines (C-binding)
    /// \param[in] fnc  (optional) name of calling function
    void operator() (int errcode, const char*func = 0) const
    {
      const int size = MPI_MAX_ERROR_STRING+512;
      char buffer[size], *buf=buffer;
      int  len = func?
	SNprintf(buf,size,"MPI Error in %s: ",func) :
	SNprintf(buf,size,"MPI Error: ") ;
      buf += len;
      MPI_Error_string(errcode, buf, &len);
      WDutils::Error(file,line)(buffer);
    }
  };
#  define MPI_THROWER  ::Error
#endif
#define MPI_THROW    MPI_THROWER(__FILE__,__LINE__)
  //
  double Itime = 0.;
  extern "C" void AbortMPI() {
    if(MPI::Initialized())
      MPI_Abort(MPI::World.Comm(),0);
  }
  // data type support
  const unsigned MaxNewTypes = 32, MaxTypeNameLen=16; 
  unsigned NewTypes = 0;
  MPI::DataType NewType[MaxNewTypes] = {0};
  char NewTypeName[MaxNewTypes][MaxTypeNameLen] = {{0}};
  //
  inline const char* name(MPI::DataType T) {
    if(T > MPI::MaxPresetType) {
      for(unsigned t=0; t!=NewTypes; ++t)
	if(T == NewType[t]) return NewTypeName[t];
    } else {
      switch(T) {
      case MPI::Byte:   return "Byte";
      case MPI::Bool:   return "Bool";
      case MPI::Int8:   return "Int8";
      case MPI::Int16:  return "Int16";
      case MPI::Int32:  return "Int32";
      case MPI::Int64:  return "Int64";
      case MPI::Uint8:  return "Uint8";
      case MPI::Uint16: return "Uint16";
      case MPI::Uint32: return "Uint32";
      case MPI::Uint64: return "Uint64";
      case MPI::Float:  return "Float";
      case MPI::Double: return "Double";
      default:          return "unknown";
      }
    }
    return "unknown";
  }
  //
  inline MPI_Datatype type(MPI::DataType T) {
    switch(T) {
    case MPI::Byte:   return MPI_BYTE;
    case MPI::Bool:   return MPI_CHAR;
    case MPI::Int8:   return MPI_CHAR;
    case MPI::Int16:  return MPI_SHORT;
    case MPI::Int32:  return MPI_INT;
    case MPI::Int64:  return MPI_LONG_LONG_INT;
    case MPI::Uint8:  return MPI_UNSIGNED_CHAR;
    case MPI::Uint16: return MPI_UNSIGNED_SHORT;
    case MPI::Uint32: return MPI_UNSIGNED;
    case MPI::Uint64: return MPI_LONG_LONG_INT;       // not quite correct
    case MPI::Float:  return MPI_FLOAT;
    case MPI::Double: return MPI_DOUBLE;
    default:          return T;
    }
  }
  //
  inline int oper(MPI::Operator O) {
    switch(O) {
    case MPI::Max: return MPI_MAX;
    case MPI::Min: return MPI_MIN;
    case MPI::Sum: return MPI_SUM;
    case MPI::Prod: return MPI_PROD;
    case MPI::And: return MPI_LAND;
    case MPI::Or: return MPI_LOR;
    case MPI::XOr: return MPI_LXOR;
    case MPI::BAnd: return MPI_BAND;
    case MPI::BOr: return MPI_BOR;
    case MPI::BXOr: return MPI_BXOR;
    default: WDutils_Error("unknown MPI::Operator\n"); // to make compiler happy
      return MPI_MAX;
    }
  }
  // ensure that:
  // -- MPI_Comm is int
  // -- MPI_Datatype is MPI::DataType
  // -- MPI::Status and MPI_Status are of the same size
  // -- MPI::Request and MPI_Request are of the same size
  struct __TMP {
    WDutilsStaticAssert((meta::TypeCompare<int,MPI_Comm>::identical  &&
			 meta::TypeCompare<MPI::DataType,MPI_Datatype>::identical &&
			 sizeof(MPI::Status)  == sizeof(MPI_Status)  &&
			 sizeof(MPI::Request) == sizeof(MPI_Request) ));
  };
} // namespace {
// /////////////////////////////////////////////////////////////////////////////
namespace WDutils {
namespace MPI {
  // ///////////////////////////////////////////////////////////////////////////
  // static data with external linkage
  Communicator Communicator::world;
  Communicator const& World = Communicator::world;
  int Communicator::DEBUG_LEVEL = 1;
  int Undefined = MPI_UNDEFINED;
  // ///////////////////////////////////////////////////////////////////////////
  // MPI::ContiguousType()
  DataType ContiguousType(int num, DataType etype, const char*tname)
			  WDutils_THROWING
  {
    // check whether there is still space for a new DataType
    if(NewTypes == MaxNewTypes)
      WDutils_THROW("MPI::ContiguousType(): cannot support more than %du"
		    "additional DataTypes\n",MaxNewTypes);
    // construct and commit the new type
    DataType ntype;
    int err;
    err = MPI_Type_contiguous(num,type(etype),&ntype);
    if(err != MPI_SUCCESS) MPI_THROW(err,"MPI::ContiguousType()");
    err = MPI_Type_commit    (&ntype);
    if(err != MPI_SUCCESS) MPI_THROW(err,"MPI::ContiguousType()");
    // check consistency
    if(ntype <= MaxPresetType)
      WDutils_THROW("MPI::ContiguousType(): internal problem\n");
    // add to list of new DataTypes
    NewType[NewTypes] = ntype;
    std::strncpy(NewTypeName[NewTypes],tname,MaxTypeNameLen-1);
    NewTypes++;
    return ntype;
  }
  // ///////////////////////////////////////////////////////////////////////////
  // MPI::Status
  unsigned Status::count(DataType t) const WDutils_THROWING
  {
    int c, err = MPI_Get_count(reinterpret_cast<MPI_Status*>
			       (const_cast<Status*>(this)), type(t), &c);
    if(err != MPI_SUCCESS) MPI_THROW(err,"Status::count()");
    if(c < 0) WDutils_THROW("Status::count(%s): c=%d<0\n",name(t),c);
    return c;
  }
  //
  unsigned Status::source() const
  {
    return reinterpret_cast<const MPI_Status*>(this)->MPI_SOURCE;
  }
  //
  int Status::tag() const
  {
    return reinterpret_cast<const MPI_Status*>(this)->MPI_TAG;
  }
  //
  int Status::errcode() const
  {
    return reinterpret_cast<const MPI_Status*>(this)->MPI_ERROR;
  }
  //
  unsigned Status::write_error(char*buf) const WDutils_THROWING
  {
    int len=0,err=errcode();
    if(err!=MPI_SUCCESS) MPI_Error_string(err, buf, &len);
    if(len < 0) WDutils_THROW("Status::write_error(): len=%d<0\n",len);
    return len;
  }
  // ///////////////////////////////////////////////////////////////////////////
  // MPI::Init()
  void Init(int&argc, const char**&argv,
	    const char**VAR, const char*file) WDutils_THROWING
  {
    if(!Initialized()) {

      // 1 start MPI and set initial wall-clock time
      int err;
      err = MPI_Init(&argc,const_cast<char***>(&argv));
      if(err != MPI_SUCCESS) MPI_THROW(err,"MPI::Init()");
      Itime = MPI_Wtime();
      // 2 set World and initialize it
      const_cast<Communicator*>(&World)->COMM = MPI_COMM_WORLD;
      const_cast<Communicator*>(&World)->init();
      // 3 set exit() functions and tell about MPI
      std::atexit(&AbortMPI);
      RunInfo::set_mpi_proc(World.rank(),World.size());
      // 4 ensure environment variables are set
      //   this is more complicated than you may think
      if(VAR && VAR[0]) {
	const unsigned NMAX=128;
	unsigned NENV=0;
	for(NENV=0; VAR[NENV] && NENV!=NMAX; ++NENV) ;
	if(NENV == NMAX)
	  WDutils_THROW("MPI::Init(): # environments vars exceeds %d",NMAX);
	DebugInfo(8,"MPI::Init(): trying to set %du environment vars\n",NENV);
	// 4.1 Do we have first environment variable on all processes?
	char haveEnv = getenv(VAR[0])!=0;
	char*HaveEnv = WDutils_NEW(char,World.size());
	err = MPI_Allgather(&haveEnv,1,MPI_CHAR,HaveEnv,1,
			    MPI_CHAR,MPI_COMM_WORLD);
	if(err != MPI_SUCCESS) MPI_THROW(err,"MPI::Init()");
	bool AllHaveEnv = true;
	for(unsigned p=0; p!=World.size(); ++p)
	  if(HaveEnv[p]==0) { AllHaveEnv=false; break; }
	if(AllHaveEnv) {
	  //     YES! all have VAR[0]: we are done
	  if(World.rank() == 0)
	    DebugInfo(8,"MPI::Init(): all processes have full environment\n");
	} else if(HaveEnv[0]) {
	  const size_t SENV=1024;
	  // 4.2 NO, only root has VAR[0]: getenv() on root, then broadcast
	  if(World.rank() == 0) {
	    DebugInfo(10,"MPI::Init(): only root as full environment "
		      "but will broadcast it now...\n");
	    for(unsigned i=0; i!=NENV; ++i) {
	      char  *envv=getenv(VAR[i]);
	      if(envv==0)
		WDutils_THROW("MPI::Init(): cannot getenv(\"%s\")\n",VAR[i]);
	      size_t lenv=strlen(envv);
	      if(lenv>=SENV)
		WDutils_THROW("MPI::Init(): strlen(\"$%s\")=%lu >= %lu\n",
			      VAR[i],lenv,SENV);
	      err = MPI_Bcast(envv,strlen(envv)+1,MPI_CHAR,0,MPI_COMM_WORLD);
	      if(err != MPI_SUCCESS) MPI_THROW(err,"MPI::Init()");
	    }
	    DebugInfo(8,"MPI::Init(): done: broadcasting environment\n");
	  } else {
	    char ENVV[SENV];
	    for(unsigned i=0; i!=NENV; ++i) {
	      err = MPI_Bcast(ENVV,SENV,MPI_CHAR,0,MPI_COMM_WORLD);
	      if(err != MPI_SUCCESS) MPI_THROW(err,"MPI::Init()");
	      setenv(VAR[i],ENVV,0);
	    }
	    DebugInfo(8,
		      "MPI::Init(): successfully set environment variables\n");
	  }
	} else {
	  // 4.3 NO, nobody has VAR[0]: try to get env vars from file
	  if(file==0)
	    WDutils_THROW("MPI::Init(): "
			  "unable to set environment variables: "
			  "please provide file argument to read them from\n");
	  std::ifstream IN(file);
	  if(!IN)
	    WDutils_THROW("MPI::Init(): "
			  "unable to set environment variables: "
			  "cannot open file \"%s\"\n",file);
	  DebugInfo(10,
		    "MPI::Init(): will try to get environment variables from "
		    "file \"%s\"\n",file);
	  char     ReadEnv[NMAX]={0};
	  unsigned Read=0;
	  std::string EVAR,ENVV;
	  while(IN && Read!=NENV) {
	    while(eat_line(IN,'#')) ;
	    IN >> EVAR >> ENVV;
	    for(unsigned i=0; i!=NENV; ++i)
	      if(0==strcmp(EVAR.c_str(),VAR[i])) {
		if(ReadEnv[i])
		  WDutils_THROW("MPI::Init(): found environment variable \"%s\""
				" twice in file \"%s\"\n",VAR[i],file);
		DebugInfo(10,"MPI::Init(): setting $%s to \"%s\"\n",
			  VAR[i], ENVV.c_str());
		setenv(VAR[i], ENVV.c_str(),1);
		ReadEnv[i] = 1;
		++Read;
		break;
	      }
	  }
	  if(Read != NENV) {
	    DebugInfo(0,"MPI::Init(): only found %d instead of %d environment "
		      "variables in file \"%s\"\n",Read,NENV,file);
	    for(unsigned i=0; i!=NENV; ++i)
	      if(ReadEnv[i]==0)
		WDutils_THROW("MPI::Init(): couldn't find environment variable "
			      "\"%s\" in \"%s\"\n",VAR[i],file);
	  }
	  DebugInfo(8,"MPI::Init(): successfully set environment variables\n");
	}
	WDutils_DEL_A(HaveEnv);
      }
    }
  }
  // ///////////////////////////////////////////////////////////////////////////
  // MPI::Finish()
  void Finish() WDutils_THROWING
  {
    int err;
    err = MPI_Finalize();
    if(err != MPI_SUCCESS) goto MPIError;
    const_cast<Communicator*>(&World)->COMM = 0;
    const_cast<Communicator*>(&World)->SIZE = 0;
    const_cast<Communicator*>(&World)->RANK = 0;
    return;
  MPIError: MPI_THROW(err,"MPI::Finish()");
  }
  // ///////////////////////////////////////////////////////////////////////////
  // MPI::WallClock()
  double WallClock() {
    return MPI_Wtime() - Itime;
  }
  // ///////////////////////////////////////////////////////////////////////////
  // MPI::WallClockTick()
  double WallClockTick() {
    return MPI_Wtick();
  }
  // ///////////////////////////////////////////////////////////////////////////
  // MPI::Communicator
#define DEBUGINFO if(debug(DEBUG_LEVEL)) DebugInformation(FILE,LINE)
#undef  MPI_THROW
#define MPI_THROW MPI_THROWER(FILE,LINE)
  //
  void Communicator::init() WDutils_THROWING
  {
    int err,tmp;
    err = MPI_Errhandler_set(COMM,MPI_ERRORS_RETURN);
    if(err != MPI_SUCCESS) goto MPIError;
    err = MPI_Comm_size(COMM,&tmp); SIZE = tmp;
    if(err != MPI_SUCCESS) goto MPIError;
    err = MPI_Comm_rank(COMM,&tmp); RANK = tmp;
    if(err != MPI_SUCCESS) goto MPIError;
    return;
  MPIError: MPI_THROW(err,"MPI::Communicator()");
  }
  //
  void Communicator::Barrier() const {
    DEBUGINFO("MPI::Communicator::Barrier()\n");
    MPI_Barrier(COMM);
    FILE = 0;
  };
  //
  void Communicator::BroadCast(unsigned root, void*buf, unsigned count,
			       DataType t, const char*caller)
    const WDutils_THROWING
  {
    if(!caller) DEBUGINFO("MPI::Communicator::BroadCast(): %d -> %d %s\n",
			  root,count,name(t));
    int err = MPI_Bcast(buf,static_cast<int>(count),type(t),
			static_cast<int>(root),COMM);
    if(err != MPI_SUCCESS) MPI_THROW(err,caller? caller : "BroadCast()");
    FILE = 0;
  }
  //
  void Communicator::Gather(unsigned root, const void*send, unsigned count,
			    void*recv, DataType t, const char*caller)
    const WDutils_THROWING
  {
    if(!caller) DEBUGINFO("MPI::Communicator::Gather(): %d %s at %d\n",
			  count,name(t),root);
    if(count == 0)
      Warning(FILE,LINE)("MPI::Communicator::Gather(): count=0\n");
    else {
      if(recv == 0 && root == RANK)
	WDutils_THROWER(FILE,LINE)("MPI::Communicator::Gather(): "
				  "recv=0 at root\n");
      int err =
	MPI_Gather(const_cast<void*>(send), static_cast<int>(count), type(t),
		   recv, static_cast<int>(count), type(t),
		   static_cast<int>(root), COMM);
      if(err != MPI_SUCCESS) MPI_THROW(err,caller? caller : "Gather()");
    }
  }
  //
  void Communicator::GatherVar(unsigned root,
			       const void*send, unsigned sendcount,
			       void*recv, const unsigned*recvcounts,
			       DataType t, const char*caller) const
    WDutils_THROWING
  {
    if(!caller) DEBUGINFO("MPI::Communicator::GatherVar(): %u %s at %u\n",
			  sendcount,name(t),root);
    int err;
    if(RANK == root) {
      if(recv == 0)
	WDutils_THROWER(FILE,LINE)("MPI::Communicator::GatherVar(): "
				  "recv=0 at root\n");
      if(recvcounts == 0)
	WDutils_THROWER(FILE,LINE)("MPI::Communicator::GatherVar(): "
				  "recvcounts=0 at root\n");
      if(sendcount != recvcounts[root])
	WDutils_THROWER(FILE,LINE)("MPI::Communicator::GatherVar(): at root: "
				  "sendcount=%u != recvcounts[root]=%u\n",
				  sendcount, recvcounts[root]);
      Array<int> displs(SIZE);
      Array<int> counts(SIZE);
      displs[0] = 0;
      counts[0] = recvcounts[0];
      for(unsigned p=1; p!=SIZE; ++p) {
	counts[p] = recvcounts[p];
	displs[p] = displs[p-1] + recvcounts[p-1];
      }
      err =
	MPI_Gatherv(const_cast<void*>(send), static_cast<int>(sendcount),
		    type(t), recv, counts.array(),
		    displs.array(), type(t), static_cast<int>(root), COMM);
    } else {
      err = MPI_Gatherv(const_cast<void*>(send), static_cast<int>(sendcount),
			type(t), 0,0,0, type(t), static_cast<int>(root), COMM);
    }
    if(err != MPI_SUCCESS) MPI_THROW(err,caller? caller : "GatherVar()");
    FILE = 0;
  }
  //
  void Communicator::AllGather(const void*send, unsigned count, void*recv,
			       DataType t, const char*caller)
    const WDutils_THROWING
  {
    if(!caller) DEBUGINFO("MPI::Communicator::AllGather(): %d %s\n",
			  count,name(t));
    if(count == 0)
      Warning(FILE,LINE)("MPI::Communicator::AllGather(): no data to send\n");
    else {
      int err = MPI_Allgather(const_cast<void*>(send), static_cast<int>(count),
			      type(t), recv,  static_cast<int>(count), type(t),
			      COMM);
      if(err != MPI_SUCCESS) MPI_THROW(err,caller? caller: "Allgather()");
    }
    FILE = 0;
  }
  //
  void Communicator::AllToAll(const void*send, unsigned count, void*recv,
			      DataType t, const char*caller)
    const WDutils_THROWING
  {
    if(!caller) DEBUGINFO("MPI::Communicator::AllToAll() %d %s\n",
			  count,name(t));
    int err = MPI_Alltoall(const_cast<void*>(send), static_cast<int>(count),
			   type(t), recv, static_cast<int>(count), type(t),
			   COMM);
    if(err != MPI_SUCCESS) MPI_THROW(err,caller? caller: "AllToAll()");
  }
  //
  void Communicator::Reduce(Operator op, unsigned root, const void*send,
			    void*recv, unsigned count, DataType t,
			    const char* caller) const WDutils_THROWING
  {
    if(!caller) DEBUGINFO("MPI::Communicator::Reduce(%s) %u %s\n",
			  OpName(op),count,name(t));
    int err = MPI_Reduce(const_cast<void*>(send), recv, static_cast<int>(count),
			 type(t), oper(op), static_cast<int>(root), COMM);
    if(err != MPI_SUCCESS)
      MPI_THROW(err, message(caller? caller : "Reduce(%s)",OpName(op)));
    FILE = 0;
  }
  //
  void Communicator::AllReduce(Operator op, const void*send, void*recv,
			       unsigned count, DataType t, const char*caller)
    const WDutils_THROWING
  {
    if(!caller) DEBUGINFO("MPI::Communicator::AllReduce(%s) %d %s\n",
			  OpName(op),count,name(t));
    int err = MPI_Allreduce(const_cast<void*>(send), recv,
			    static_cast<int>(count),type(t), oper(op),COMM);
    if(err != MPI_SUCCESS)
      MPI_THROW(err,message(caller? caller : "AllReduce(%s)",OpName(op)));
    FILE = 0;
  }
  //
  void Communicator::Send(unsigned dest, const void*buf, unsigned count,
			  int tag, DataType t, const char*caller)
    const WDutils_THROWING
  {
    if(!caller) DEBUGINFO("MPI::Communicator::Send() %u %s -> %u\n",
			  count,name(t),dest);
    int err = MPI_Send(const_cast<void*>(buf),static_cast<int>(count),type(t),
		       static_cast<int>(dest),tag,COMM);
    if(err != MPI_SUCCESS) MPI_THROW(err,caller? caller: "Send()");
    FILE = 0;
  }
  //
  unsigned Communicator::Recv(unsigned srce, void*buf, unsigned count, int tag,
			      DataType t, bool warn, const char*caller)
    const WDutils_THROWING
  {
    if(!caller) DEBUGINFO("MPI::Communicator::Recv() %d %s <- %d\n",
			  count,name(t),srce);
    MPI_Status S;
    int err, actual_count;
    err = MPI_Recv(buf, static_cast<int>(count), type(t),
		   static_cast<int>(srce), tag, COMM, &S);
    if(err != MPI_SUCCESS) MPI_THROW(err,"Recv()");
    err = MPI_Get_count(&S, type(t), &actual_count);
    if(err != MPI_SUCCESS) MPI_THROW(err,"Recv()");
    if(actual_count<0)
      WDutils_THROWER(FILE,LINE)("MPI::Communicator::Recv(%s): "
				"actual_count=%d<0\n", name(t),actual_count);
    if(warn && count != static_cast<unsigned>(actual_count))
      Warning(FILE,LINE)("MPI::Communicator::Recv(): "
			 "expected %u %s but got %d "
			 "(set last arg to false to suppress this warning\n",
			 count,name(t),actual_count);
    FILE = 0;
    return actual_count;
  }
  //
  Request Communicator::IssueSend(unsigned dest, const void*buf, unsigned count,
				  int tag, DataType t, const char*caller)
    const WDutils_THROWING
  {
    if(!caller) DEBUGINFO("MPI::Communicator::IssueSend() %d %s -> %d\n",
			  count,name(t),dest);
    MPI_Request REQ;
    int err = MPI_Isend(const_cast<void*>(buf),static_cast<int>(count),type(t),
			static_cast<int>(dest),tag,COMM,&REQ);
    if(err != MPI_SUCCESS) MPI_THROW(err,caller? caller : "IssueSend()");
    FILE = 0;
    return Request(REQ);
  }
  //
  Request Communicator::IssueRecv(unsigned srce, void*buf, unsigned count,
				  int tag, DataType t, const char*caller)
    const WDutils_THROWING
  {
    if(!caller) DEBUGINFO("MPI::Communicator::IssueRecv() %d %s <- %d\n",
			  count,name(t),srce);
    MPI_Request REQ;
    int err = MPI_Irecv(buf,static_cast<int>(count),type(t),
			static_cast<int>(srce),tag,COMM,&REQ);
    if(err != MPI_SUCCESS) MPI_THROW(err,caller? caller : "IssueRecv()");
    FILE = 0;
    return Request(REQ);
  }
  //
  void Communicator::WaitAll(unsigned num, Request*request, Status*status)
    const WDutils_THROWING
  {
    DEBUGINFO("calling MPI::Communicator::WaitAll()\n");
    Array<Status> S;
    if(0==status) {
      S.reset(num);
      status = S.array();
    }
    int err = MPI_Waitall(static_cast<int>(num),
			  reinterpret_cast<MPI_Request*>(request),
			  reinterpret_cast<MPI_Status *>(status));
    if(err != MPI_SUCCESS) MPI_THROW(err,"WaitAll()");
    FILE = 0;
  }
  //
  bool Communicator::Test(Request&request, Status&status) const WDutils_THROWING
  {
    DEBUGINFO("calling MPI::Communicator::Test()\n");
    int flag;
    int err = MPI_Test(reinterpret_cast<MPI_Request*>(&request),&flag,
		       reinterpret_cast<MPI_Status *>(&status));
    if(err != MPI_SUCCESS) MPI_THROW(err,"Test()");
    FILE = 0;
    return flag;
  }
  //
  bool Communicator::TestAny(unsigned num, Request*request, int&index,
			     Status&status)
    const WDutils_THROWING
  {
    DEBUGINFO("calling MPI::Communicator::TestAny()\n");
    int flag;
    int err = MPI_Testany(static_cast<int>(num),
			  reinterpret_cast<MPI_Request*>(request),
			  &index, &flag,
			  reinterpret_cast<MPI_Status *>(&status));
    if(err != MPI_SUCCESS) MPI_THROW(err,"TestAny()");
    FILE = 0;
    return flag;
  }
  //
  unsigned Communicator::TestSome(unsigned num, Request*request,
				  unsigned*indices,
				  Status*status) const WDutils_THROWING
  {
    DEBUGINFO("calling MPI::Communicator::TestSome()\n");
    int outcount;
    Array<int> ind(num);
    int err = MPI_Testsome(static_cast<int>(num),
			   reinterpret_cast<MPI_Request*>(request),
			   &outcount,ind.array(),
			   reinterpret_cast<MPI_Status *>(status));
    if(err != MPI_SUCCESS) MPI_THROW(err,"TestSome()");
    if(outcount<0)
      Warning(FILE,LINE)("MPI::Communicator::TestSome(): "
			 "none of the requests appears active\n");
    else 
      for(int i=0; i!=outcount; ++i)
	indices[i]=ind[i];
    FILE = 0;
    return outcount<0? 0:outcount;
  }
  //
} }
////////////////////////////////////////////////////////////////////////////////

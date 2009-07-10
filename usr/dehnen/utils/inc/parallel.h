// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/parallel.h
///
/// \author Walter Dehnen
///                                                                             
/// \date   2008,2009
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008,2009 Walter Dehnen
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
#ifndef WDutils_included_parallel_h
#define WDutils_included_parallel_h

#ifndef WDutils_included_memory_h
#  include <memory.h>
#endif
// /////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  /// \brief   defines a partial C++ interface for MPI.
  ///
  /// \details
  /// Here, we provide a partial C++ interface to some MPI functionality.
  /// The reason for not using the "official" C or C++ binding for MPI are
  /// as follows.\n
  /// 1. Both C and C++ binding are very basic, which has the benefit of
  ///    great flexibility but little convenience.\n
  /// 2. The "official" C++ binding is only a straight "translation" into C++
  ///    without any try to make it a proper C++ (STL-style) code.\n
  /// 3. Our implementation here is driven by our needs in falcON and provides
  ///    basic two-point and collective communications of all built-in data
  ///    types, plus contiguous types built from these, either as array or
  ///    Array<> (which allows some additional checks on buffer size).\n
namespace MPI {
  /// \name MPI::DataType and related functionality
  //@{
  typedef int DataType;
  enum {
    Byte  = 1,   ///< refers to a single byte of raw memory
    Bool  = 2,   ///< refers to a C++ type bool
    Int8  = 3,   ///< refers to a signed 8bit integer, usually char
    Int16 = 4,   ///< refers to a signed 16bit integer, usually short
    Int32 = 5,   ///< refers to a signed 32bit integer, usually int
    Int64 = 6,   ///< refers to a signed 64bit integer
    Uint8 = 7,   ///< refers to a unsigned 8bit integer, usually unsigned char
    Uint16= 8,   ///< refers to a unsigned 16bit integer, usually unsigned short
    Uint32= 9,   ///< refers to a unsigned 32bit integer, usually unsigned int
    Uint64= 10,  ///< refers to a unsigned 64bit integer
    Float = 11,  ///< refers to C type float: 32bit 
    Double= 12   ///< refers to C type double: 64bit
  };
  const static DataType MaxPresetType = 12;
  // auxiliary struct template
  // 1 generic template definition, must be incomplete so it cannot be used
  template<typename T> struct __DataType;
  // 2 specialisations for supported types
  template<> struct __DataType<bool> {
    static DataType type() { return Bool; }
  };
  template<> struct __DataType<int8> {
    static DataType type() { return Int8; }
  };
  template<> struct __DataType<int16> {
    static DataType type() { return Int16; }
  };
  template<> struct __DataType<int32> {
    static DataType type() { return Int32; }
  };
  template<> struct __DataType<int64> {
    static DataType type() { return Int64; }
  };
  template<> struct __DataType<uint8> {
    static DataType type() { return Uint8; }
  };
  template<> struct __DataType<uint16> {
    static DataType type() { return Uint16; }
  };
  template<> struct __DataType<uint32> {
    static DataType type() { return Uint32; }
  };
  template<> struct __DataType<uint64> {
    static DataType type() { return Uint64; }
  };
  template<> struct __DataType<float> {
    static DataType type() { return Float; }
  };
  template<> struct __DataType<double> {
    static DataType type() { return Double; }
  };
  /// return MPI::DataType given a type.
  /// use as in \code Type<vect>(); \endcode
  template<typename T> inline DataType Type() WDutils_THROWING {
    return __DataType<T>::type();
  }
  /// construct a contiguous DataType
  /// \param[in] num  number of elements
  /// \param[in] type datatype of elements
  /// \param[in] name name for new DataType
  /// \return         handle to be used for new type
  /// \note May only be called a limited number of times, so use with care!
  DataType ContiguousType(int num, DataType type, const char*name)
    WDutils_THROWING;
  //@}
  /// implementing a status
  class Status {
    unsigned X[4];
  public:
    /// number of data of given type
    /// \param[in] t data type of communication
    unsigned count(DataType t) const WDutils_THROWING;
    /// source of communication
    unsigned source() const;
    /// tag for communication
    int tag() const;
    /// error code, if any, for communication
    int errcode() const;
    /// write error message
    /// \param[out] buf buffer for error message, must hold up to 512
    /// \return actual size of error message
    unsigned write_error(char*buf) const WDutils_THROWING;
  };
  /// implementing a request
  class Request {
    void*R;
  public:
    Request() : R(0) {}
    explicit Request(void*r) : R(r) {}
    Request(Request const&r) : R(r.R) {}
    operator bool() const { return R != 0; }
  };
  /// return for undefined index values etc.
  extern int    Undefined;
  //----------------------------------------------------------------------------
  /// \name some global functions which do not need a communicator
  //@{
  /// is MPI initialised?
  /// \return true if MPI::Init() has been called
  inline bool Initialized();
  /// \brief Initialise MPI and more.
  /// \details
  /// Here, the following things are done (in this order).
  /// - MPI_Init() is called spawning the MPI processes and distributing the
  ///   command-line arguments.\n
  /// - the initial wall-clock time is remembered.\n
  /// - MPI::World::init() is called.\n
  /// - std::atexit() is called to allow MPI_Abort() on std::exit().\n
  /// - RunInfo is informed about MPI and proc number.\n
  /// - optionally some environment variables are spawned.\n
  /// \param[in,out] argc  number of command-line arguments
  /// \param[in,out] argv  values of command-line arguments
  /// \param[in]     env   null-terminated list of environment variables
  /// \param[in]     file  file to look for environment variables
  /// \note arguments \a argc and \a argv are input on the master node and
  ///       output on all others: this call spawns the command-line arguments
  ///       via MPI_Init()
  /// \note argument \a env must be the same on all nodes. We will try to
  ///       ensure that all nodes have the environment variables in the list
  ///       env. If they don't initially, we try to find them at root and,
  ///       failing that, read them from the \a file provided, if any. It is
  ///       assumed that if the first environment variable is present, so are
  ///       all the others.
  void Init(int&argc, const char**&argv,
	    const char**env=0, const char*file=0) WDutils_THROWING;
  /// Finalise MPI and destruct World
  void Finish() WDutils_THROWING;
  /// wall-clock time in seconds since call to MPI::Init()
  double WallClock();
  /// tick size in seconds of wall-clock in Wtime()
  double WallClockTick();
  //@}
  //----------------------------------------------------------------------------
  /// \brief operators for collective reductions
  enum Operator {
    Max,    ///< maximum
    Min,    ///< minimum
    Sum,    ///< sum
    Prod,   ///< product
    And,    ///< logical and
    Or,     ///< logical or
    XOr,    ///< logical exclusive or
    BAnd,   ///< bit-wise and
    BOr,    ///< bit-wise or
    BXOr    ///< bit-wise exclusive or
  };
  /// \brief name for reduction operator
  inline const char* OpName(Operator Op) {
    switch(Op) {
    case Max:    return "Max";
    case Min:    return "Min";
    case Sum:    return "Sum";
    case Prod:   return "Prod";
    case And:    return "And";
    case Or:     return "Or";
    case XOr:    return "XOr";
    case BAnd:   return "BAnd";
    case BOr:    return "BOr";
    case BXOr:   return "BXOr";
    default:     return "Unknown";     // to make compiler happy
    }
  }
  //----------------------------------------------------------------------------
#define NonTemplateCommunications
#undef  NonTemplateCommunications
  /// \brief
  /// implements some functionality of MPI communicators as C++ class.
  /// \details
  /// We provide the point-to-point and collective communications as templates
  /// over datatype with the same type for send and receive.
  /// \note 
  /// We deviate from MPI standard by setting the error handler to 
  /// MPI_ERRORS_RETURN (MPI standard: MPI_ERRORS_ARE_FATAL), to enable us
  /// to throw exceptions in the member methods.
  /// \note
  /// This is not complete: I will add stuff as I need it
  class Communicator {
    int                 COMM;          ///< MPI handle for communicator
    unsigned            SIZE;          ///< size of communicator
    unsigned            RANK;          ///< rank within communicator
  protected:
    virtual ~Communicator() {}
    mutable const char *FILE;          ///< name of file from which called
    mutable int         LINE;          ///< line in file from which called
    static int          DEBUG_LEVEL;
#define DEBUGINFO if(debug(DEBUG_LEVEL)) DebugInformation(FILE,LINE)
    /// sets SIZE and RANK as well as the error handler to MPI_ERRORS_RETURN
    void init() WDutils_THROWING;
    // Init() must call init() on World, for it is constructed before Init()
    friend void Init(int&, const char**&, const char**, const char*)
      WDutils_THROWING;
    // Finish() will undo the call to init() and set COMM to NULL
    friend void Finish() WDutils_THROWING;
    //--------------------------------------------------------------------------
  public:
    int Comm() const { return COMM; }
    static Communicator world;
    static void SetDebugLevel(int d) { DEBUG_LEVEL = d; }
    static int const&DebugLevel() { return DEBUG_LEVEL; }
    /// default constructor: World
    /// \note This constructor is (implicitly) called to construct the static
    /// member world, which is visible outside the class as MPI::World. Since
    /// this is done \b before MPI::Init() is called in main(), we cannot do
    /// any MPI related stuff here. This is deferred to member init(), which
    /// is called for Communicator::world from MPI::Init(), which in turn also
    /// sets member World.COMM to MPI_COMM_WORLD.
    /// \note There is only one World, so a second call is errorneous.
    Communicator() WDutils_THROWING : COMM(0)
    {
      if(Initialized())
	throw exception("trying to construct another MPI::World");
    }
    /// size of communicator
    unsigned const&size() const { return SIZE; }
    /// rank within communicator
    unsigned const&rank() const { return RANK; }
    /// equality
    bool operator==(Communicator const&C) const { return COMM == C.COMM; }
    /// inequality
    bool operator!=(Communicator const&C) const { return COMM != C.COMM; }
    /// set file and line information
    virtual Communicator* set(const char*f, int l)
    {
      FILE = f;
      LINE = l;
      return this;
    }
    /// set file and line information
    virtual const Communicator* set(const char*f, int l) const
    {
      FILE = f;
      LINE = l;
      return this;
    }
#define COMMUN(C) (C)->set(__FILE__,__LINE__)
    /// \name collective communications (not complete)
    //@{
    /// Barrier
    void Barrier() const;
    // --- BroadCast ----------------------------------------------------------
    void BroadCast(unsigned root, void*buf, unsigned count, DataType type,
		   const char* =0) const WDutils_THROWING;
    /// broadcast a buffer
    /// \param[in]     root      sender
    /// \param[in,out] buf       buffer both for send and receive
    /// \param[in]     count     amount of data to be send
    template<typename T>
    void BroadCast(unsigned root, T*buf, unsigned count) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::BroadCast(): %d -> %d %s\n",
		root,count,nameof(T));
      BroadCast(root,buf,count,Type<T>(),"BroadCast()");
    }
    /// broadcast a single datum
    /// \param[in]      root      sender
    /// \param[in,out]  buf       buffer both for send and receive
    template<typename T>
    void BroadCast(unsigned root, T&buf) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::BroadCast(): %d -> 1 %s\n",
		root,nameof(T));
      BroadCast(root,&buf,1,Type<T>(),"BroadCast()");
    }
    /// broadcast an Array<>
    /// \param[in]     root      sender
    /// \param[in,out] buf       buffer both for send and receive
    /// \note This routine is a specialisation of BCast with a single datum
    template<typename T>
    void BroadCast(unsigned root, Array<T>&buf) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::BroadCast(): %d -> Array<%s>(%d)\n",
		root,nameof(T),buf.size());
      BroadCast(root,buf.array(),buf.size(),Type<T>(),"BroadCast()");
    }
    // --- Gather --------------------------------------------------------------
    void Gather(unsigned root, const void*send, unsigned count,
		void*recv, DataType type, const char* =0)
      const WDutils_THROWING;
    /// gather multiple data at root
    /// \param[in]  root   recipient
    /// \param[in]  send   buffer with data to be send
    /// \param[in]  count  amount of data to be send
    /// \param[out] recv   at root only: buffer for data to be received
    /// \note count must be the same on all processes and at root recv must
    ///       hold memory for count*this->size() elements
    template<typename T>
    void Gather(unsigned root, const T*send, unsigned count, T*recv=0)
      const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::Gather(): %d %s at %d\n",
		count,nameof(T),root);
      Gather(root,send,count,recv,Type<T>(),"Gather()");
    }
    /// gather a single datum into Array<>
    /// \param[in]  root  recipient
    /// \param[in]  send  datum to be gather by root
    /// \param[out] recv  buffer for data to be received at root
    /// \note If at root recv.size() < this->size(), we reset recv
    template<typename T>
    void Gather(unsigned root, T const&send, Array<T>&recv=Array<T>())
      const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::Gather(): 1 %s at %d\n",
		nameof(T),root);
      if(root == RANK && recv.size() < SIZE) {
	Warning(FILE,SIZE)("MPI::Communicator::Gather(): at root: "
			   "recv.size()=%d < this->size()=%d: "
			   "will reset recv.size() to %d\n",
			   recv.size(),SIZE,SIZE);
	recv.reset(SIZE);
      }
      Gather(root,&send,1,recv.array(),Type<T>(),"Gather()");
    }
    /// gather multiple data using Array<> arguments
    /// \param[in]  root  recipient
    /// \param[in]  send  data to be send to root
    /// \param[out] recv  at root only: buffer for data to be received
    /// \note send.size() must be the same on each process (not checked)
    /// \note We enforce that at root recv has correct size to hold the data
    template<typename T>
    void Gather(unsigned root, Array<T> const&send,
		Array<T,2>&recv=(Array<T,2>()))
      const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::Gather(): Array<%s>(%d) at %d\n",
		nameof(T),send.size(),root);
      if(root == RANK && (recv.size(0)!=SIZE || recv.size(1)!=send.size())) {
	Warning(FILE,SIZE)("MPI::Communicator::Gather(): at root:"
			   "recv.size()=[%d,%d] != "
			   "this->size(),send.size()=[%d,%d]: "
			   "will reset recv.size() to [%d,%d]\n",
			   recv.size(0),recv.size(1),SIZE,send.size(),
			   SIZE,send.size());
	unsigned n[2] = {SIZE,send.size()};
	recv.reset(n);
      }
      Gather(root,send.array(),send.size(),recv.array(),Type<T>(),"Gather()");
    }
    // --- GatherVar -----------------------------------------------------------
    void GatherVar(unsigned root, const void*send, unsigned sendcount,
		   void*recv, const unsigned*recvcounts, DataType type,
		   const char* =0) const WDutils_THROWING;
    /// contiguously gather multiple data of variable size: template over type
    /// \param[in]  root       recipient
    /// \param[in]  send       buffer with data to be send
    /// \param[in]  sendcount  amount of data to be send
    /// \param[out] recv       at root only: buffer for data to be received
    /// \param[in]  recvcounts at root only: # data to be received per proc
    template<typename T>
    void GatherVar(unsigned root, const T*send, unsigned sendcount,
		   T*recv=0, const unsigned*recvcounts=0) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::GatherVar(): %d %s at %d\n",
		sendcount,nameof(T),root);
      GatherVar(root,send,sendcount,recv,recvcounts,Type<T>(),"GatherVar()");
    }
    /// contiguously gather multiple data of variable size, using Array<> args
    /// \param[in]  root       recipient
    /// \param[in]  send       Array with data to be send from this proc
    /// \param[in]  counts     Array: # data gathered from proc p
    /// \param[out] recv       at root only: Array for data to be received
    /// \note Unlike the non-Array<> version, we require counts to be correct
    ///       not just at root, but at all processes, allowing more stringent
    ///       error checking.
    /// \note It is not permissible that send.size() != counts[RANK], but
    ///       it is permissible that recv.size() >= Sum counts[p]
    template<typename T>
    void GatherVar(unsigned root, Array<T> const&send,
		   Array<unsigned> const&counts,
		   Array<T>&recv = (Array<T>())) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::GatherVar(): Array<%s>(%d) at %d\n",
		nameof(T),send.size(),root);
      if(counts.size() < SIZE)
	WDutils_THROWER(FILE,LINE)("MPI::Communicator::GatherVar(): "
				  "counts.size()=%d < this->size()=%d\n",
				  counts.size(),SIZE);
      if(send.size() != counts[RANK])
	WDutils_THROWER(FILE,LINE)("MPI::Communicator::GatherVar(): "
				  "send.size()=%d != counts[RANK]=%d\n",
				  send.size(),counts[RANK]);
      if(RANK == root) {
	unsigned total_count=0;
	for(unsigned p=0; p!=SIZE; ++p) total_count += counts[p];
	if(total_count > recv.size()) {
	  Warning(FILE,SIZE)("MPI::Communicator::GatherVar(): at root: "
			     "recv.size()=%d < Sum counts[p]=%d: "
			     "will reset recv.size() to %d\n",
			     recv.size(),total_count,total_count);
	  recv.reset(total_count);
	}
      }
      GatherVar(root,send.array(),send.size(),
		recv.array(),counts.array(),Type<T>(),"GatherVar()");
    }
    // --- AllGather -----------------------------------------------------------
    void AllGather(const void*send, unsigned count, void*recv, DataType type,
		   const char* =0) const WDutils_THROWING;
    /// gather data at all processes: template over type
    /// \param[in]  send   buffer with data to be send
    /// \param[in]  count  amount of data to be send and received
    /// \param[out] recv   buffer for data to be received
    /// \note count must be the same on all processes (not checked)
    /// \note recv must hold memory for count*this->size() elements
    template<typename T>
    void AllGather(const T*send, unsigned count, T*recv) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::AllGather(): %d %s\n", count,nameof(T));
      AllGather(send,count,recv,Type<T>(),"AllGather()");
    }
    /// gather single datum at all processes
    /// \param[in]  send  datum to be send
    /// \param[out] recv  buffer for data to be received
    /// \note If needed, we reset recv to match size requirements
    template<typename T>
      void AllGather(const T&send, Array<T>&recv) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::AllGather(): 1 %s\n", nameof(T));
      if(recv.size() != SIZE) {
	Warning(FILE,SIZE)("MPI::Communicator::AllGather(): "
			   "recv.size() mismatch ([%d] vs [%d]:"
			   "will reset sizes\n",recv.size(),SIZE);
	recv.reset(SIZE);
      }
      AllGather(&send,1,recv.array(),Type<T>(),"AllGather()");
    }
    /// gather data at all processes, using Array<> arguments
    /// \param[in]  send   buffer with data to be send
    /// \param[out] recv   buffer for data to be received
    /// \note send.size() must be the same on all processes (not checked)
    /// \note If needed, we reset recv to match size requirements
    template<typename T>
      void AllGather(const Array<T>&send, Array<T,2>&recv)
      const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::AllGather(): Array<%s>(%d)\n",
		nameof(T),send.size());
      if(send.size() == 0)
	Warning(FILE,LINE)("MPI::Communicator::AllGather(): no data to send\n");
      else {
	if(recv.size(0) != SIZE || recv.size(1) != send.size()) {
	  Warning(FILE,LINE)("MPI::Communicator::AllGather(): "
			     "recv.size() mismatch ([%d,%d] vs [%d,%d]:"
			     "will reset recv.size() to [%d,%d]\n",
			     recv.size(0),recv.size(1),
			     SIZE,send.size(),SIZE,send.size());
	  unsigned n[2] = {SIZE,send.size()};
	  recv.reset(n);
	}
	Allgather(send.array(),send.size(),recv.array(),Type<T>(),
		  "AllGather()");
      }
    }
    // --- AllToAll -----------------------------------------------------------
    void AllToAll(const void*send, unsigned count, void*recv, DataType type,
		  const char* =0) const WDutils_THROWING;
    /// send data between all processes
    /// \param[in]  send   data to be send to each process
    /// \param[in]  count  # data to be send to each process
    /// \param[out] recv   buffer for data to be received from all processes
    /// \note send must contain count*this->size() data and, equally,
    ///       recv must be able to hold the same amount of data
    template<typename T>
    void AllToAll(const T*send, unsigned count, T*recv) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::AllToAll() %d %s\n",count,nameof(T));
      AllToAll(send,count,recv,Type<T>(),"AllToAll()");
    }
    /// send data between all processes using Array<> arguments
    /// \param[in]  send   data to be send to each process
    /// \param[out] recv   buffer for data to be received from all processes
    /// \note send.size(0) must equal this->size(); the same holds for recv.
    ///       Similarly, send.size(1) must equal recv.size(1). We reset recv
    ///       if necessary.
    template<typename T>
    void AllToAll(Array<T,2>const&send, Array<T,2>&recv) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::AllToAll() Array<%s>(%d)\n",
		nameof(T),send.size(1));
      if(send.size(0) != SIZE)
	WDutils_THROWER(FILE,LINE)("MPI::Communicator::AllToAll(): "
				  "send.size(0)=%d != this->size=%d\n",
				  send.size(0),SIZE);
      if(recv.size(0) != SIZE || recv.size(1) != send.size(1)) {
	Warning(FILE,LINE)("MPI::Communicator::AllToAll(): "
			   "recv.size() mismatch ([%d,%d] vs [%d,%d]:"
			   "will reset recv.size() to [%d,%d]\n",
			   recv.size(0),recv.size(1),
			   SIZE,send.size(1),SIZE,send.size(1));
	unsigned n[2] = {SIZE,send.size(1)};
	recv.reset(n);
      }
      AllToAll(send.array(),send.size(1),recv.array(),Type<T>(),"AllToAll()");
    }
    // --- Reduce -------------------------------------------------------------
    void Reduce(Operator op, unsigned root, const void*send, void*recv,
		unsigned count, DataType type, const char* =0)
      const WDutils_THROWING;
    /// reduce data at one process: template over type and reduction operator
    /// \param[in]  root   recipient
    /// \param[in]  send   buffer with data to be reduced
    /// \param[out] recv   at root only: buffer for data to be received
    /// \param[in]  count  amount of data to be reduced
    /// \note count must be the same on all processes (not checked)
    template<Operator O, typename T>
    void Reduce(unsigned root, const T*send, T*recv, unsigned count)
      const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::Reduce<%s> %d %s\n",
		OpName(O),count,nameof(T));
      Reduce(O,root,send,recv,count,Type<T>(),"Reduce<%s>");
    }
    /// reduce a single datum at one process
    /// \param[in]  root      recipient
    /// \param[in]  send      datum to be reduced
    /// \param[out] recv      datum to be received
    template<Operator O, typename T>
    void Reduce(unsigned root, T const&send, T&recv) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::Reduce<%s> 1 %s\n",
		OpName(O),nameof(T));
      Reduce(O,root,&send,&recv,1,Type<T>(),"Reduce<%s>");
    }
    /// reduce data at one process, using Array<> arguments
    /// \param[in]  root   recipient
    /// \param[in]  send   buffer with data to be reduced
    /// \param[out] recv   at root only: buffer for data to be received
    /// \note send.size() must be the same on all processes (not checked)
    /// \note At root, we enforce that recv.size() >= send.size()
    /// \note This routine is a specialisation of Reduce with a single datum
    template<Operator O, typename T>
    void Reduce(unsigned root, Array<T>const&send, Array<T>&recv=(Array<T>()))
      const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::Reduce<%s> Array<%s>(%d)\n",
		OpName(O),nameof(T),send.size());
      if(RANK==root && send.size() > recv.size()) {
	Warning(FILE,LINE)("MPI::Communicator::Reduce<%s> :"
			   "recv.size()=%d < send.size()=%d: "
			   "will reset recv.size() to %d",
			   OpName(O),recv.size(),send.size(),send.size());
	recv.reset(send.size());
      }
      Reduce(O,root,send.array(),recv.array(),send.size(),Type<T>(),
	     "Reduce<%s>");
    }
    // --- ReduceInPlace ------------------------------------------------------
    /// reduce data at one process in place: template over type & reduction op
    /// \param[in]     root  recipient
    /// \param[in,out] buf   buffer with data to be reduced
    /// \param[in]     count amount of data to be reduced
    /// \note at root buf will on output contain the reduce data
    template<Operator O, typename T>
    void ReduceInPlace(unsigned root, T*buf, unsigned count)
      const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::ReduceInPlace<%s> %d %s\n",
		OpName(O),count,nameof(T));
      if(RANK==root) {
	T*send = WDutils_NEW(T,count);
	memcpy(send,buf,count*sizeof(T));
	Reduce(O,root,send,buf,count,Type<T>(),"ReduceInPlace<%s>");
	WDutils_DEL_A(send);
      } else 
	Reduce(O,root,buf,   0,count,Type<T>(),"ReduceInPlace<%s>");
    }
    /// reduce a single datum in place
    /// \param[in]     root     recipient
    /// \param[in,out] buf      datum to be reduced (output only at root)
    template<Operator O, typename T>
    void ReduceInPlace(unsigned root, T&buf) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::ReduceInPlace<%s> 1 %s\n",
		OpName(O),nameof(T));
      if(RANK==root) {
	T send(buf);
	Reduce(O,root,&send,&buf,1,Type<T>(),"ReduceInPlace<%s>");
      } else
	Reduce(O,root,&buf ,   0,1,Type<T>(),"ReduceInPlace<%s>");
    }
    /// reduce data at one process in place, using Array<> argument
    /// \param[in]     root  recipient
    /// \param[in,out] buf   buffer with data to be reduced
    /// \note at root buf will on output contain the reduce data
    /// \note This is a specialisation of ReduceInPlace with a single datum
    template<Operator O, typename T>
    void ReduceInPlace(unsigned root, Array<T,1>&buf) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::ReduceInPlace<%s> Array<%s>(%d)\n",
		OpName(O),nameof(T),buf.size());
      if(RANK==root) {
	T*send = WDutils_NEW(T,buf.size());
	memcpy(send,buf.array(),buf.size()*sizeof(T));
	Reduce(O,root,send,buf.array(),buf.size(),Type<T>(),
	       "ReduceInPlace<%s>");
	WDutils_DEL_A(send);
      } else 
	Reduce(O,root,buf.array(),0,buf.size(),Type<T>(),"ReduceInPlace<%s>");
    }
    // --- AllReduce ----------------------------------------------------------
    void AllReduce(Operator op, const void*send, void*recv, unsigned count,
		   DataType type, const char* =0) const WDutils_THROWING;
    /// reduce data at all processes: template over type and reduction operator
    /// \param[in]  send   buffer with data to be reduced
    /// \param[out] recv   buffer for data to be received
    /// \param[in]  count  amount of data to be reduced
    template<Operator O, typename T>
    void AllReduce(const T*send, T*recv, unsigned count) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::AllReduce<%s> %d %s\n",
		OpName(O),count,nameof(T));
      AllReduce(O,send,recv,count,Type<T>(),"AllReduce<%s>");
    }
    /// reduce a single datum at all processes
    /// \param[in]  send   datum to be reduced
    /// \param[out] recv   reduced datum
    template<Operator O, typename T>
    void AllReduce(T const&send, T&recv) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::AllReduce<%s> 1 %s\n", OpName(O),nameof(T));
      AllReduce(O,&send,&recv,1,Type<T>(),"AllReduce<%s>");
    }
    /// reduce data at all processes, using Array<> arguments
    /// \param[in]  send   buffer with data to be reduced
    /// \param[out] recv   buffer for data to be received
    /// \note we enforce recv to hold enough memory to receive data
    /// \note This is a specialisation of AllReduce with a single datum
    template<Operator O, typename T>
    void AllReduce(Array<T> const&send, Array<T> &recv) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::AllReduce(%s) Array<%s>(%d)\n",
		OpName(O),nameof(T),send.size());
      if(send.size() > recv.size()) {
	Warning(FILE,LINE)("MPI::Communicator::AllReduce<%s>: "
			   "recv.size()=%d < send.size()=%d: "
			   "will reset recv.size() to %d\n",
			   OpName(O),recv.size(),send.size(),send.size());
	recv.reset(send.size());
      }
      AllReduce(O,send.array(),recv.array(),send.size(),Type<T>(),
		"AllReduce<%s>");
    }
    // --- AllReduceInPlace ---------------------------------------------------
    /// reduce data at all processes in place
    /// \param[in,out] buf    buffer with data to be reduced
    /// \param[in]     count  amount of data to be reduced
    /// \note on output buf will contain the reduced data
    template<Operator O, typename T>
    void AllReduceInPlace(T*buf, unsigned count) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::AllReduceInPlace<%s> %d %s\n",
		OpName(O),count,nameof(T));
      T*send = WDutils_NEW(T,count);
      memcpy(send,buf,count*sizeof(T));
      AllReduce(O,send,buf,count,Type<T>(),"AllReduceInPlace<%s>");
      WDutils_DEL_A(send);
    }
    /// reduce a single datum at all processes in place
    /// \param[in,out]  buf    datum to be reduced
    template<Operator O, typename T>
    void AllReduceInPlace(T&buf) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::AllReduceInPlace<%s> 1 %s\n",
		OpName(O),nameof(T));
      T send(buf);
      AllReduce(O,&send,&buf,1,Type<T>(),"AllReduceInPlace<%s>");
    }
    /// reduce data at all processes in place, using Array<> argument
    /// \param[in,out] buf    buffer with data to be reduced
    /// \note on output buf will contain the reduce data
    /// \note This is a specialisation of AllReduceInPlace with a single datum
    template<Operator O, typename T>
    void AllReduceInPlace(Array<T>&buf) const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::AllReduceInPlace<%s> Array<%s>(%d)\n",
		OpName(O),nameof(T),buf.size());
      T*send = WDutils_NEW(T,buf.size());
      memcpy(send,buf.array(),buf.size()*sizeof(T));
      AllReduce(O,send,buf.array(),buf.size(),Type<T>(),"AllReduceInPlace<%s>");
      WDutils_DEL_A(send);
    }
    //@}
    /// \name blocking point to point communications (not complete)
    //@{
    // --- Send ---------------------------------------------------------------
    void Send(unsigned dest, const void*buf, unsigned count, int tag, DataType type,
	      const char* =0) const WDutils_THROWING;
    /// send data
    /// \param[in]  dest  destination process
    /// \param[in]  buf   buffer with data to be send
    /// \param[in]  count amount of data to be send
    /// \param[in]  tag   identifier
    template<typename T>
    void Send(unsigned dest, const T*buf, unsigned count, int tag)
      const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::Send() %d %s -> %d\n",count,nameof(T),dest);
      Send(dest,buf,count,tag,Type<T>(),"Send()");
    }
    // --- Recv ---------------------------------------------------------------
    unsigned Recv(unsigned srce, void*buf, unsigned count, int tag,
		  DataType type, bool warn, const char* =0)
      const WDutils_THROWING;
    /// receive data
    /// \return            # data actually received
    /// \param[in]  srce   source process: only receive matching message
    /// \param[in]  buf    buffer to receive data
    /// \param[in]  count  size of receive buffer
    /// \param[in]  tag    identifier: only receive matching message
    /// \param[in]  warn   optional: warn about count mismatch (default: true)
    template<typename T>
    unsigned Recv(unsigned srce, T*buf, unsigned count, int tag, bool warn=true)
      const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::Recv() %d %s <- %d\n",
		count,nameof(T),srce);
      return Recv(srce,buf,count,tag,Type<T>(),warn,"Recv()");
    }
    //@}
    /// \name non-blocking point to point communications (not complete)
    //@{
    // --- IssueSend ----------------------------------------------------------
    Request IssueSend(unsigned dest, const void*buf, unsigned count, int tag,
		      DataType type, const char* =0) const WDutils_THROWING;
    /// issue send request
    /// \param[in]  dest  destination process
    /// \param[in]  buf   buffer with data to be send
    /// \param[in]  count amount of data to be send
    /// \param[in]  tag   identifier
    /// \return           request issued
    template<typename T>
    Request IssueSend(unsigned dest, const T*buf, unsigned count, int tag)
      const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::IssueSend() %d %s -> %d\n",
		count,nameof(T),dest);
      return IssueSend(dest,buf,count,tag,Type<T>(),"IssueSend()");
    }
    // --- IssueRecv ----------------------------------------------------------
    Request IssueRecv(unsigned srce, void*buf, unsigned count, int tag,
		      DataType type, const char* =0) const WDutils_THROWING;
    /// issue receive request
    /// \param[in]  srce  source process: only receive matching message
    /// \param[in]  buf   buffer to receive data
    /// \param[in]  count size of receive buffer
    /// \param[in]  tag   identifier: only receive matching message
    /// \return           request issued
    template<typename T>
    Request IssueRecv(unsigned srce, T*buf, unsigned count, int tag)
      const WDutils_THROWING
    {
      DEBUGINFO("MPI::Communicator::IssueRecv() %d %s <- %d\n",
		count,nameof(T),srce);
      return IssueRecv(srce,buf,count,tag,Type<T>(),"IssueRecv()");
    }
    // --- Waits and Tests ----------------------------------------------------
    /// blocking: wait for a completion of a operations
    /// \param[in,out] request  Request as issued by IssueSend() or IssueRecv()
    /// \param[out]    status   status of operation
    /// \note the request is void after return.
    void Wait(Request&request, Status&status) const WDutils_THROWING;
    /// blocking: wait for completion of any of several requests
    /// \param[in]     num     number of operations to be completed
    /// \param[in,out] reqs    array of requests to be completed
    /// \param[out]    status  optional: status of completed operation
    /// \return                index of completed operation
    /// \note The request[index] is void after return.
    /// \note If none of the requests was valid, MPI::Undefined is returned
    unsigned WaitAny(unsigned num, Request*reqs, Status*status=0)
      const WDutils_THROWING;
    /// blocking: wait for completion of all of several operations
    /// \param[in]     num     number of operations to be completed
    /// \param[in,out] reqs    array of requests to be completed
    /// \param[out]    status  optional: status of completed operations
    /// \note The request[] are void after return.
    /// \note If none of the requests was valid, MPI::Undefined is returned
    void WaitAll(unsigned num, Request*reqs, Status*status=0)
      const WDutils_THROWING;
    /// non-blocking: test for completion of a single operation
    /// \param[in,out] request Request as issued by IssueSend() or IssueRecv()
    /// \param[out]    status  if success: status of operation
    /// \return        true    if operation is complete
    /// \note if true is returned, the Request is void after the call
    bool Test(Request&request, Status&status) const WDutils_THROWING;
    /// non-blocking: test for completion of one or none of severeal operations
    /// \param[in]     num     number of operations to be tested
    /// \param[in,out] request array of requests to be tested
    /// \param[out]    index   index of completed operation, if any
    /// \note index is set to MPI::Undefined if no operation completed.
    /// \param[out]    status  status of completed operation, if any
    /// \return true if either one operation completed or no request was active
    /// \note If one operation finished, the associated request becomes void.
    /// \note If none of the requests was valid, true is returned and index is
    ///       MPI::Undefined.
    bool TestAny(unsigned num, Request*request, int&index, Status&status)
      const WDutils_THROWING;
    /// non-blocking: test for completion of any of several operations
    /// \param[in]     num      number of operations to be tested
    /// \param[in,out] request  array: requests to be tested
    /// \param[in,out] indices  indices of operations completed, if any
    /// \param[out]    status   array: statuses of completed operations, if any
    /// \return                 number of completed operations, if any
    /// \note if no active request was in the list, 0 is returned but a
    ///       warning is issued
    unsigned TestSome(unsigned num, Request*request, unsigned*indices,
		      Status*status) const WDutils_THROWING;
    //@}
    /// shift rank right
    /// \param[in,out] r  will be increased by one, modulo SIZE
    void to_right(unsigned&r) const {
      if(++r == SIZE) r=0u;
    }
    /// shift rank left
    /// \param[in,out] r  will be decreased by one, modulo SIZE
    void to_left(unsigned&r) const {
      if(r == 0u) r = SIZE;
      --r;
    }
    /// right of given rank
    /// \param[in] r  given rank
    /// \return    r+1 mod SIZE
    unsigned right(unsigned r) const {
      to_right(r);
      return r;
    }
    /// left of given rank
    /// \param[in] r  given rank
    /// \return    r-1 mod SIZE
    unsigned left(unsigned r) const {
      to_left(r);
      return r;
    }
    /// right of us
    unsigned right() const { return right(RANK); }
    /// left of us
    unsigned left() const { return left(RANK); }
#undef DEBUGINFO
  };// class Communicator
  /// Communicator World
  extern Communicator const&World;
  //----------------------------------------------------------------------------
  inline bool Initialized() { return World.Comm() != 0; }
} } // namespace WDutils { namespace MPI {
// /////////////////////////////////////////////////////////////////////////////
#ifndef LoopProcLeft
#  define LoopProcLeft(COMM,LEFT)			\
  for(unsigned LEFT=COMM.left();			\
      LEFT!=COMM.rank();				\
      COMM.to_left(LEFT))
#endif

#ifndef LoopProcRight
#  define LoopProcRight(COMM,RIGHT)			\
  for(unsigned RIGHT=COMM.right();			\
      RIGHT!=COMM.rank();				\
      COMM.to_right(RIGHT))
#endif

#ifndef LoopProcLeftRight
#  define LoopProcLeftRight(COMM,LEFT,RIGHT)		\
  for(unsigned LEFT=COMM.left(), RIGHT=COMM.right();	\
      LEFTT!=COMM.rank();				\
      COMM.to_left(LEFT), COMM.to_right(RIGHT))
#endif
// /////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_parallel_h

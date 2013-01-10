// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/memory.h
///
/// \author Walter Dehnen
///                                                                             
/// \date   2000-2012
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2012 Walter Dehnen
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
#ifndef WDutils_included_memory_h
#define WDutils_included_memory_h

#ifndef WDutils_included_cstddef
#  include <cstddef>
#  define WDutils_included_cstddef
#endif
#ifndef WDutils_included_cstring
#  include <cstring>                           // for memcpy
#  define WDutils_included_cstring
#endif
#ifndef WDutils_included_iostream
#  include <iostream>
#  define WDutils_included_iostream
#endif
#ifndef WDutils_included_traits_h
#  include <traits.h>
#endif
#ifndef WDutils_included_inline_h
#  include <inline.h>
#endif
#ifndef WDutils_included_cstdlib
#  include <cstdlib>
#endif

#if __cplusplus < 201103L
# define noexcept
# define constexpr
#endif

namespace WDutils {
#ifndef WDutilsAllocDebugLevel
# define WDutilsAllocDebugLevel 8
#endif
  // ///////////////////////////////////////////////////////////////////////////
  //
  /// array allocation giving useful info in case of allocation error; mostly
  /// used from macro WDutils_NEW.
  ///
  /// In case of allocation error we abort or throw an exception (depending on
  /// the WDutils error handling settings).  If the debugging level exceeds
  /// WDutilsAllocDebugLevel (default 8), we always print debugging
  /// information about memory allocation.
  ///
  /// \return         a valid pointer (unless an error occurs) 
  /// \param[in] num  number of array elements to allocate
  /// \param[in] file name of the source file where this routines is called
  /// \param[in] line number of the line in that file 
  /// \param[in] lib  (optional) name of calling library (default: "WDutils")
  template<typename T> inline
  T* NewArray(size_t num, const char*file, int line, const char*lib="WDutils")
    WDutils_THROWING
  {
    T*t;
    bool failed=0;
    try {
      t = new T[num];
    } catch(std::bad_alloc E) {
      t = 0;
      failed = 1;
    }
    if(failed || (num && t==0))
      throw Thrower(file,line)("allocation of %u '%s' (%u bytes) failed\n",
			       uint32_t(num),nameof(T),uint32_t(num*sizeof(T)));
    DebugInformation(file,line,lib)(WDutilsAllocDebugLevel,
				    "allocated %u %s = %u bytes @ %p\n",
				    uint32_t(num),nameof(T),
				    uint32_t(num*sizeof(T)),
				    static_cast<void*>(t));
    return t;
  }
  // ///////////////////////////////////////////////////////////////////////////
  //
  /// C MACRO to be used for array allocation
  ///
  /// Calling WDutils::NewArray<TYPE>(), which in case of an error generates
  /// an error message detailing the source file and line of the call. In case
  /// the debugging level exceeds 10, we always print debugging information
  /// about memory allocation.
  ///
  /// \param  TYPE name of the element type
  /// \param  SIZE number of elements
#define WDutils_NEW(TYPE,SIZE)				\
  WDutils::NewArray<TYPE>(SIZE,__FILE__,__LINE__)
  template<typename T> struct _report
  {
    static bool report(int n) { return n>0; }
    static int  bytes (int n) { return n*sizeof(T); }
  };
  template<> struct _report<void>
  {
    static bool report(int) { return 0; }
    static int  bytes (int) { return 0; }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //
  /// array de-allocation giving useful info in case of error; mostly used
  /// from macro WDutils_DEL_A.
  ///
  /// In case of a de-allocation error (if the pointer provided was not valid)
  /// an error is generated (or an exception thrown, depending on the WDutils
  /// error settings). If the debugging level exceeds WDutilsAllocDebugLevel
  /// (default 8), we always print debugging information about memory
  /// de-allocation.
  ///                                                                           
  /// \param[in] array pointer previously allocated with WDutils::NewArray<>()
  ///                  or ::operator new[].
  /// \param[in] file  name of the source file where this routines is called
  /// \param[in] line  number of the line in that file
  /// \param[in] num   (optional) number of elements de-allocated
  /// \param[in] lib   (optional) name of calling library (default: "WDutils")
  template<typename T> inline
  void DelArray(const T*array, const char*file, int line, int num=0,
		const char*lib = "WDutils")
    WDutils_THROWING
  {
    if(0==array) {
      Warning(file,line,lib)
	("trying to delete zero pointer to array of '%s'", nameof(T));
      return;
    }
    try {
      delete[] array;
    } catch(...) {
      throw Thrower(file,line)
	("de-allocating array of '%s' @ %p failed\n", nameof(T),array);
    }
    if(_report<T>::report(num))
      DebugInformation(file,line,lib)
	(WDutilsAllocDebugLevel, "de-allocated array of %d %s [%d byes] @ %p\n",
	 num, nameof(T), _report<T>::bytes(num), array);
    else
      DebugInformation(file,line,lib)
	(WDutilsAllocDebugLevel, "de-allocated array of %s @ %p\n",
	 nameof(T), array);
  }
  // ///////////////////////////////////////////////////////////////////////////
  //
  /// C MACRO to be used for array de-allocation
  ///
  /// Calling WDutilsN::DelArray<TYPE>(), which in case of an error generates
  /// an error message detailing the source file and line of the call. In case
  /// the debugging level exceeds 10, we always print debugging information
  /// about memory de-allocation.
  ///
  /// \param P  pointer to be de-allocated
#define WDutils_DEL_A(P)			\
  WDutils::DelArray(P,__FILE__,__LINE__)
#define WDutils_DEL_AN(P,N)			\
  WDutils::DelArray(P,__FILE__,__LINE__,N)
  // ///////////////////////////////////////////////////////////////////////////
  //
  /// Object de-allocation giving useful info in case of error; mostly used
  /// from macro WDutils_DEL_O.
  ///
  /// In case of a de-allocation error (if the pointer provided was not valid)
  /// an error is generated (or an exception thrown, depending on the WDutils
  /// error settings).
  ///
  /// \param[in] pobj  pointer previously allocated with ::operator new().
  /// \param[in] file  name of the source file where this routines is called
  /// \param[in] line  number of the line in that file
  /// \param[in] lib   (optional) name of calling library (default: "WDutils")
  template<typename T> inline
  void DelObject(const T*pobj, const char*file, int line,
		 const char*lib="WDutils") WDutils_THROWING
  {
    if(0==pobj) {
      Warning(file,line,lib)
	("trying to delete zero pointer to object '%s'", nameof(T));
      return;
    }
    try {
      delete pobj;
    } catch(...) {
      throw Thrower(file,line)
	("de-allocating object '%s' @ %p failed\n", nameof(T),pobj);
    }
    DebugInformation(file,line,lib)
      (WDutilsAllocDebugLevel,"de-allocated %s object @ %p\n", nameof(T),pobj);
  }
  // ///////////////////////////////////////////////////////////////////////////
  ///
  /// C MACRO to be used for object de-allocation
  ///
  /// should be used for all object de-allocation in WDutils.
  ///
  /// Calling WDutils::DelObject<TYPE>(), which in case of an error generates
  /// an error message detailing the source file and line of the call.
  ///
  /// \param P pointer to object to be de-allocated
#define WDutils_DEL_O(P) \
  WDutils::DelObject(P,__FILE__,__LINE__)
}
//
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#include <mm_malloc.h>
#endif
//
namespace WDutils {
  //////////////////////////////////////////////////////////////////////////////
  /// \defgroup  MemAligned  memory alignment to K bytes

  /// Macro enforcing memory alignment to K bytes
  /// \ingroup MemAligned
  ///
  /// Forces the corresponding variable/type to be aligned to K bytes; Works
  /// with icc (icpc) and gcc (g++) [versions > 3]; Use it like \code
  ///    struct WDutils__align_to(K) name { ... };              \endcode
#if defined (__INTEL_COMPILER)
#  define WDutilsAlignTo(K) __declspec(align(K))
#elif (defined (__GNUC__) && __GNUC__ > 2) || defined(__PGI)
#  define WDutilsAlignTo(K) __attribute__ ((aligned(K)))
#else
#  define WDutilsAlignTo(K)
#endif
  /// struct holding aligned datum
  template<int alignment, typename data_type> struct AlignedDatum
  {
    WDutilsStaticAssert(alignment<=1);
    data_type Datum;
  };
  //
  template<typename data_type> struct AlignedDatum<2,data_type>
  { WDutilsAlignTo(2) data_type Datum; };
  //
  template<typename data_type> struct AlignedDatum<4,data_type>
  { WDutilsAlignTo(4) data_type Datum; };
  //
  template<typename data_type> struct AlignedDatum<8,data_type>
  { WDutilsAlignTo(8) data_type Datum; };
  //
  template<typename data_type> struct AlignedDatum<16,data_type>
  { WDutilsAlignTo(16) data_type Datum; };
  //
  template<typename data_type> struct AlignedDatum<32,data_type>
  { WDutilsAlignTo(32) data_type Datum; };

  /// Macro enforcing memory alignment to 16 bytes
  /// \ingroup MemAligned
#define WDutils__align16 WDutilsAlignTo(16)
  /// Macro enforcing memory alignment to 32 bytes
  /// \ingroup MemAligned
#define WDutils__align32 WDutilsAlignTo(32)

  ///
  /// is a given memory address aligned to K bytes?
  /// \ingroup MemAligned
  /// \param p  memory address to be tested
  /// \param al alignemt to a bytes will be tested
  template<int K>
  inline constexpr bool is_aligned(const void*p)
  { return size_t(p) % K == 0; }
  ///
  /// find the smallest multiple of K not smaller than @a n
  /// \note @a K must be power of 2
  template<int K>
  inline constexpr size_t next_aligned(size_t n)
  {
    static_assert(K>0 && (K&(K-1))==0,"K not a power of 2");
    return (n+K-1)&(~(K-1));
  }
  ///
  /// find the smallest K-byte aligned address not smaller than @a p
  template<int K, typename T>
  inline T* next_aligned(T*p)
  { return reinterpret_cast<T*>(next_aligned<K>(reinterpret_cast<size_t>(p))); }
#if(0)
  ///
  /// is a given memory address aligned to a 16 bytes memory location?
  /// \param p  memory address to be tested
  inline bool is_aligned16(const void*p)
  { return size_t(p) % 16 == 0; }
  ///
  /// find the smallest multiple of 16 not smaller than @a n
  inline size_t next_aligned16(size_t n)
  { return (n+15)&(~15); }
  ///
  /// find the smallest 16-byte aligned address not smaller than @a p
  template<typename T>
  inline T* next_aligned16(T*p)
  { return reinterpret_cast<T*>(next_aligned16(reinterpret_cast<size_t>(p))); }
#endif
  ///
  /// Allocate memory at a address aligned to a K byte memory location
  /// \ingroup MemAligned
  ///
  /// \return a newly allocated memory address at a K byte memory location
  /// \param[in] nobj  number of objects to allocate
  /// \param[in] file  name of the source file where this routines is called
  /// \param[in] line  number of the line in that file
  /// \param[in] lib   (optional) name of calling library (default: "WDutils")
  /// \version  debugged 02-09-2004 WD
  /// \note Unlike NewArray<>, we do not call the default ctor for each
  ///       allocated object!
  template<int K, typename T> inline
  T* NewArrayAligned(size_t nobj, const char*file, int line,
		     const char*lib = "WDutils") WDutils_THROWING
  {
    static_assert(K>0 && (K&(K-1))==0,"K not power of 2");
    size_t nbytes = nobj*sizeof(T);
#if defined(__GNUC__) || defined (__INTEL_COMPILER)
    void*t;
    bool failed=0;
    try {
      t = _mm_malloc(nbytes,K);
    } catch(...) {
      t = 0;
      failed = 1;
    }
    if(failed || (nbytes && t==0))
      throw Thrower(file,line)
	("NewArrayAligned<%d,%s>(%u): allocation of %u bytes failed\n",
	 K,nameof(T),uint32_t(nobj),uint32_t(nbytes));
    DebugInformation(file,line,lib)
      (WDutilsAllocDebugLevel,
       "allocated %u %s = %u bytes aligned to %d @ %p\n",
       uint32_t(nobj),nameof(T),uint32_t(nbytes),K,t);
    return static_cast<T*>(t);
#else // __GNUC__ or __INTEL_COMPILER
    static_assert(K==16,"only implemented for K=16 (for other than gcc & icc)");
    // linear memory model:                                                     
    // ^    = 16byte alignment points                                           
    // S    = sizeof(void*) (assumed 4 in this sketch)                          
    // PPPP = memory where p is stored (needed in deletion)                     
    // def0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef01    
    //    ^               ^               ^               ^               ^     
    //        |-S-|       |  ->  at least n bytes                               
    //        |   |       |                                                     
    //        p   q   ->  q                                                     
    //    |--off--|-16-off|                                                     
    //                |-S-|                                                     
    //                PPPP|                                                     
    //                                                                          
    // the original allocation gave p, we return the shifted q and remember     
    // the original allocation address at PPPP.
    char*p = NewArray<char>(nbytes+16+sizeof(void*),func,line,lib);
    char*q = p + sizeof(void*);                     // go sizeof pointer up     
    size_t off = size_t(q) % 16;                    // offset from 16b alignment
    if(off) q += 16-off;                            // IF offset, shift         
    *((void**)(q-sizeof(void*))) = p;               // remember allocation point
    return static_cast<T*>(q);                      // return aligned address   
#endif
  }
  //
  /// C MACRO to be used for array allocation aligned to K bytes
  /// \ingroup MemAligned
  ///
  /// Calling WDutils::NewArrayAligned<K,TYPE>(), which in case of an error
  /// generates an error message detailing the source file and line of the
  /// call. In case the debugging level exceeds 10, we always print debugging
  /// information about memory allocation.
  ///
  /// \param  TYPE name of the element type
  /// \param  SIZE number of elements
  /// \note   Unlike WDutils_NEW(TYPE,SIZE), we do not call the default ctor for
  ///         the objects allocated!
#define WDutils_NEW_aligned(K,TYPE,SIZE)			\
  WDutils::NewArrayAligned<K,TYPE>(SIZE,__FILE__,__LINE__)
  /// for backwards compatibility
#define WDutils_NEW16(TYPE,SIZE) WDutils_NEW_aligned(16,TYPE,SIZE)
  ///
  /// de-allocate memory previously allocated with WDutils::NewArrayAligned<K>()
  /// \ingroup MemAligned
  ///
  /// \param[in] array pointer previously allocated with
  ///                  WDutils::NewArrayAligned<K>()
  /// \param[in] file  name of the source file where this routines is called
  /// \param[in] line  number of the line in that file
  /// \param[in] num   (optional) number of elements de-allocated
  /// \param[in] lib   (optional) name of calling library (default: "WDutils")
  ///
  /// \note This routine \b must be used to properly de-allocate memory that
  ///       has been previously allocated by WDutils::NewArrayAligned<K>();
  ///       other de-allocation may result in a run-time \b error!
  template<int K, typename T> inline
  void DelArrayAligned(const T*array, const char*file, int line,
		       int num=0, const char*lib="WDutils") WDutils_THROWING
  {
#if defined(__GNUC__) || defined (__INTEL_COMPILER)
    if(0==array) {
      Warning(file,line,lib)("WDutils::DelArrayAligned<%d,%s>(0x0)\n",
			     K,nameof(T));
      return;
    }
    if(!is_aligned<K>(array)) {
      throw Thrower(file,line)
	("WDutils::DelArrayAligned<%d,%s>(%p): not aligned",K,nameof(T),array);
    }
    try {
      _mm_free(const_cast<T*>(array));
    } catch(...) {
      throw Thrower(file,line)
	("WDutils::DelArrayAligned<%d,%s>(%p): de-allocation failed\n",
	 K,nameof(T),array);
    }
    if(_report<T>::report(num))
      DebugInformation(file,line,lib)
	(WDutilsAllocDebugLevel,
	 "de-allocated %d-byte aligned array of %d '%s' [%d bytes] @ %p\n", 
	 K,num,nameof(T),_report<T>::bytes(num),array);
    else
      DebugInformation(file,line,lib)
	(WDutilsAllocDebugLevel,
	 "de-allocated %d-byte aligned array of '%s' @ %p\n",
	 K,nameof(T),array);
#else
    static_assert(K==16,"only implemented for K=16 (for other than gcc & icc)");
    DelArray((char*)(*((void**)(((char*)array)-sizeof(void*)))),
	     file,line,num,lib);
#endif
  }
  // ///////////////////////////////////////////////////////////////////////////
  //
  /// C MACRO to be used for array de-allocation of K-byte aligned stuff
  /// \ingroup MemAligned
  ///
  /// Calling WDutilsN::DelArrayAlgined<K,TYPE>(), which in case of an error
  /// generates an error message detailing the source file and line of the
  /// call. In case the debugging level exceeds 10, we always print debugging
  /// information about memory de-allocation.
  ///
  /// \param P  pointer to be de-allocated
#define WDutils_DEL_aligned(K,P)			\
  WDutils::DelArrayAligned<K>(P,__FILE__,__LINE__)
#define WDutils_DEL_alignedN(K,P,N)			\
  WDutils::DelArrayAligned<K>(P,__FILE__,__LINE__,N)
  /// for backward compatibility
#define WDutils_DEL16(P) WDutils_DEL_aligned(16,P)
#define WDutils_DEL16N(P,N) WDutils_DEL_alignedN(16,P,N)
  // ///////////////////////////////////////////////////////////////////////////
  //
  /// a simple one-dimensional array of data aligned to 16 bytes
  /// \note  sizeof(T) must be either a multiple or a dividor of 16.
  /// \ingroup MemAligned
  template<typename T>
  class Array16 {
    /// ensure sizeof(T) is either multiple or dividor of 16
    WDutilsStaticAssert( 0 == (sizeof(T) % 16)  ||  0 == (16 % sizeof(T)) );
    /// # objects to allocate for n
    static unsigned Nalloc(unsigned n)
    { return sizeof(T)>=16? n : ((n*sizeof(T)+15)&(~15))/sizeof(T); }
    //  data
    unsigned N; ///< # allocated data
    T       *A; ///< allocates array
    //  no copy ctor and no operator=
    Array16           (const Array16&) WDutilsCXX11Delete;
    Array16& operator=(const Array16&) WDutilsCXX11Delete;
  public:
    /// default ctor
    Array16()
      : N(0), A(0) {}
    /// ctor from size
    explicit Array16(unsigned n)
      : N(Nalloc(n)), A(N?WDutils_NEW16(T,N):0) {}
    /// dtor
    ~Array16()
    {
      if(A) WDutils_DEL16N(A,N);
      N = 0;
      A = 0;
    }
    /// will reset() allocate new data (and delete any old data)?
    bool reset_will_allocate(unsigned n) const
    {
      n = Nalloc(n);
      return n>N || (3*n<2*N && sizeof(T)*n>16);
    }
    /// reset(): only re-allocate if n>N or n<2N/3
    void reset(unsigned n)
    {
      n = Nalloc(n);
      if(n>N || (3*n<2*N && sizeof(T)*n>16) ) {
	if(A) WDutils_DEL16N(A,N);
	N = n;
	A = WDutils_NEW16(T,N);
      }
    }
    /// grow: increase size by n, but keep old data
    /// \param[in] n  grow by this much, default: double size
    void grow(unsigned n=0)
    {
      n = n? N+n : N+N;
      if(n) {
	T* newA = WDutils_NEW16(T,n);
	if(N) {
	  memcpy(newA,A,sizeof(T)*N);
	  WDutils_DEL16N(A,N);
	}
	N = n;
	A = newA;
      }
    }
    /// # allocated elements
    unsigned nalloc() const { return N; }
    /// const data access
    template<typename integer>
    T const&operator[](integer i) const { return A[i]; }
    /// non-const data access
    template<typename integer>
    T&operator[](integer i) { return A[i]; }
    /// const array
    const T*array() const { return A; }
    /// non-const array
    T*array() { return A; }
  };
  //
  //  WDutils::block_alloc<T>
  //
  /// allocator for elements of type T
  ///
  /// Elements are allocated in blocks and given out either singly or in small
  /// contiguous chunks. A forward iterator for sequential access to all
  /// elements given out sofar.
  ///
  /// In order to (i) guarantee the correct amount of elements (not known a
  /// priori) will be allocated and (ii) not to waste memory by allocating too
  /// many, we use the following strategy, similar to that of std::string.
  /// Elements are allocated in blocks, which in turn are organized in a
  /// linked list. When a number of new elements is to be allocated and the
  /// last block in the list cannot provide them, we allocated a new block and
  /// add it to the list. The number of elements allocated in this new block
  /// is taken to be a function of the total number of elements used
  /// sofar. This function must be provided by the user, otherwise the same
  /// number as in the last block is used.
  ///
  /// \warning Since using 16-byte aligned memory for the elements, the
  ///          default constructor is not called for each of them, rather
  ///          uninitialised memory is returned. The user should therefore
  ///          initialise elements given out, rather than trust the default
  ///          constructor to have been called.
  ///
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T, int Align=16>
  class block_alloc {
  public:
    WDutilsStaticAssert(0== (Align & (Align-1)));
    /// \name some public typedefs
    //@{
    typedef T         value_type;                ///< type of elements          
    typedef size_t    size_type;                 ///< type of number of elements
    typedef ptrdiff_t difference_type;           ///< type of pointer difference
    typedef T        *pointer;                   ///< type of pointer to element
    typedef const T  *const_pointer;             ///< type of const pointer     
    typedef T        &reference;                 ///< type of reference to elem 
    typedef const T  &const_reference;           ///< type of const reference   
    //@}
    //  no copy ctor and no operator=
    block_alloc(const block_alloc&) WDutilsCXX11Delete;
    block_alloc& operator=(const block_alloc&) WDutilsCXX11Delete;
  private:
    /// allocates and manages a contiguous chunk of elements                    
    class block {
    private:
      /// \name data of class WDutils::block_alloc::block
      //@{
      block      *NEXT;    ///< next block in linked list
      value_type *FIRST;   ///< front element
      value_type *END;     ///< end of active elements
      value_type *ENDTOT;  ///< end of all elements
      //@}
      block();
    public:
      /// constructor
      /// \param n number of elements to allocate
      explicit
      block(size_type const&n) :
	NEXT    ( 0 ),
	FIRST   ( WDutils_NEW_aligned(Align,value_type,n) ),
	END     ( FIRST ),
	ENDTOT  ( FIRST + n ) {}
      /// destructor: de-allocate
      ~block()
      { WDutils_DEL_aligned(Align,FIRST); }
      /// link block to next block
      void link(block*__n)
      { NEXT = __n; }
      /// give out: another element
      pointer new_element()
      { return END++; }
      /// give out: N new elements
      pointer new_elements(size_type n)
      {
	pointer OLD_END=END;
	END += n;
	return OLD_END;
      }
      /// can we give out up to \a n elements?
      bool has_free(size_type n) const
      { return END + n <= ENDTOT; }
      /// is block full (no free elements)?
      bool is_full() const
      { return END>=ENDTOT; }
      /// return next block in linked list of blocks
      block  *next () const
      { return NEXT; }
      /// return first element allocated
      pointer front() const
      { return FIRST; }
      /// return last element used
      pointer back () const
      { return END-1; }
      /// return end of elements used
      pointer end  () const
      { return END; }
      /// return number of elements used
      difference_type N_used () const
      { return END-FIRST; }
      /// return number of elemnets allocated
      difference_type N_alloc() const
      { return ENDTOT-FIRST; }
    };
  public:
    /// for forward sequential iteration through all elements used              
    class iterator {
    private:
      /// \name data of WDutils::block_alloc::iterator
      //@{
      block  *B;   ///< current block
      pointer E;   ///< current element
      //@}
    public:
      /// \name construction and assignment
      //@{
      /// construct from block pointer and element pointer
      iterator(block*b, pointer e)
	: B(b), E(e) {}
      /// copy constructor
      iterator(iterator const&I)
	: B(I.B), E(I.E) {}
      /// copy assignment
      iterator& operator= (iterator const&I)
      {
	B=I.B;
	E=I.E;
	return *this;
      }
      //@}
      /// prefix ++: increment pointer; if out of block: get next block
      iterator& operator++()
      {
	++E;
	if(E == B->end()) {
	  B = B->next();
	  E = B? B->front() : 0;
	}
	return *this;
      }
      /// postfix ++: return temporary equal original, but increment this
      iterator operator++(int)
      {
	iterator tmp(*this);
	this->operator++();
	return tmp;
      }
      /// equality: refer to the same element?
      bool operator==(const iterator&I) const
      { return E == I.E; }
      /// inequality: not equal
      bool operator!=(const iterator&I) const
      { return ! operator==(I); }
      /// conversion to bool: is iterator valid?
      operator bool() const
      { return B!=0; }
      /// \name conversions to pointer to elment or reference to element
      //@{
      /// type conversion to pointer to element
      operator pointer()
      { return E; }
      /// type conversion to pointer to const element
      operator const_pointer() const
      { return E; }
      /// reference operator
      reference operator* ()
      { return*E; }
      /// const reference operator
      const_reference operator* () const
      { return*E; }
      /// dereference operator
      pointer operator->()
      { return E; }
      /// const dereference operator
      const_pointer operator->() const
      { return E; }
      //@}
    };
  private:
    /// \name member data of class WDutils::block_alloc
    //@{
    block       *FIRST;  ///< first block in linked list
    block       *LAST;   ///< last block in linked list
    size_type    NTOT;   ///< # elements allocated
    size_type    NUSED;  ///< # elements used
    size_type    NBLCK;  ///< # blocks
    //@}
  public:
    /// \name member functions of class WDutils::block_alloc
    //@{
    /// constructor: allocate first block
    /// \param[in] Ns number of elements in 1st block
    explicit
    block_alloc(size_type Ns)
      : FIRST ( new block(Ns) ),
	LAST  ( FIRST ),
	NTOT  ( Ns ),
	NUSED ( 0 ),
	NBLCK ( 1 ) {}
    /// destructor: delete all blocks
    ~block_alloc() WDutils_THROWING;
    /// give out: another element
    /// \return pointer to allocated element
    pointer new_element()
    {
      if(LAST->is_full()) {
	size_type New = LAST->N_alloc();
	LAST->link(new block(New));
	LAST = LAST->next();
	NTOT+= New;
	++NBLCK;
      }
      ++NUSED;
      return LAST->new_element();
    }
    /// give out: another element
    /// \return pointer to allocated element
    /// \param F estimator = function object returning N_required(N_used yet)
    template<class estimator>
    pointer new_element(estimator const&F)
    {
      if(LAST->is_full()) {
	size_type New = F(NUSED);
	LAST->link(new block(New));
	LAST = LAST->next();
	NTOT+= New;
	++NBLCK;
      }
      ++NUSED;
      return LAST->new_element();
    }
    /// give out: more elements
    /// \return pointer to \a Ne allocated elements
    /// \param Ne number of elements to return pointer to
    pointer new_elements(size_type Ne)
    {
      if(!LAST->has_free(Ne)) {
	size_type New=max(Ne,static_cast<size_type>(LAST->N_alloc()));
	LAST->link(new block(New));
	LAST = LAST->next();
	NTOT+= New;
	++NBLCK;
      }
      NUSED+=Ne;
      return LAST->new_elements(Ne);
    }
    /// give out: more elements
    /// \return pointer to \a Ne allocated elements
    /// \param Ne number of elements to return pointer to
    /// \param F estimator = function object returning N_required(N_used yet)
    template<class estimator>
    pointer new_elements(size_type Ne, estimator const&F)
    {
      if(!LAST->has_free(Ne)) {
	size_type New=max(Ne,F(NUSED));
	LAST->link(new block(New));
	LAST = LAST->next();
	NTOT+= New;
	++NBLCK;
      }
      NUSED+=Ne;
      return LAST->new_elements(Ne);
    }
    /// an invalid number for an element
    static const unsigned invalid = ~0;
    /// the running number of a given element
    ///
    /// If @a E does not to point to an element, we return @a invalid
    ///
    /// \return the running number of element pointed to by @a E
    /// \param[in]  E  pointer to element
    unsigned number_of_element(const_pointer E) const
    {
      if(E==0) return invalid;
      unsigned num=0;
      for(block*B=FIRST; B; B=B->next())
	if(E >= B->front() && E < B->end())
	  return num + unsigned(E-B->front());
	else
	  num += B->N_used();
      return invalid;
    }
    /// element of a given running number
    ///
    /// \param[in]  n  running number, must be in [0, N_used[
    /// \return     pointer to element # @a n, NULL is @a n outside range
    pointer element(unsigned n) const
    {
      for(block*B=FIRST; B; B=B->next()) {
	unsigned nb = B->N_used();
	if(n<nb) return B->front() + n;
	n -= nb;
      }
      return 0;
    }
    /// was a given element given out by this?
    /// \return true if a given pointer to element was given out by this
    /// \param E pointer to element
    bool is_element(const_pointer E) const
    {
      if(E==0) return 0;
      for(register block*B=FIRST; B; B=B->next())
	if(E >= B->front() && E < B->end())
	  return 1;
      return 0;
    }
    /// return pointer to first element
    pointer first() const
    { return FIRST->front(); }
    /// return iterator referring to first element
    iterator front() const
    { return iterator(FIRST,FIRST->front()); }
    /// return iterator referring to first element
    iterator begin() const
    { return iterator(FIRST,FIRST->front()); }
    /// return iterator referring to end of elements
    iterator end() const
    { return iterator(0,0); }
    /// return number of elements given out sofar
    size_type N_used() const
    { return NUSED; }
    /// return number of elements allocated sofar
    size_type N_allocated() const
    { return NTOT; }
    /// return number of blocks
    size_type N_blocks() const
    { return NBLCK; }
  };// class WDutils::block_alloc
  //
  template<typename T, int K> struct traits< block_alloc<T,K> > {
    static const char  *name () {
      return message("block_alloc<%s,%d>",traits<T>::name(),K);
    }
  };
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T, int K> inline
  block_alloc<T,K>::~block_alloc() WDutils_THROWING {
    block *A=FIRST, *N;
    while(A) {
      N = A->next();
      WDutils_DEL_O(A);
      A = N;
    }
  }
  //
  //  class WDutils::pool
  //
  /// allocates blocks of K bytes (=elements) in chunks of N elements
  ///
  /// elements are defined solely by their size K in bytes. They are allocated
  /// in chunks of N, which are organized as linked list. Single elements can be
  /// handed out (allocate) or freed (de-allocated). Free elements are kept in a
  /// linked list. The actual number of bytes per element is at least the size
  /// of a pointer. Thus, for K < sizeof(void*), this class is inefficient.
  ///
  class pool {
    //  no copy ctor and no operator=
    pool(const pool&) WDutilsCXX11Delete;
    pool& operator=(const pool&) WDutilsCXX11Delete;
  public:
    /// \name sub-types of class WDutils::pool
    //@{
    typedef size_t    size_type;         ///< type of number of elements
    typedef ptrdiff_t difference_type;   ///< type of pointer difference
  public:
    /// elementary of a linked list
    struct link
    {
      link *NEXT;   ///< pter to next link
    };
#define LINK(NAME) reinterpret_cast<link*>(NAME)
    /// a chunk of elements
    struct chunk {
      char   *DATA; ///< pter to allocated memory
      chunk  *NEXT; ///< pter to next chunk
      /// constructor
      /// \param[in] N  number of element in chunk
      /// \param[in] Kp sizeof(elements)
      chunk(size_type N, size_type Kp)
	: DATA ( WDutils_NEW16(char,N*Kp) ),
	  NEXT ( 0 )
      {
	const char *END=DATA+N*Kp;
	char*l,*n;
	for(l=DATA, n=DATA+Kp; n!=END; l+=Kp,n+=Kp)
	  LINK(l)->NEXT = LINK(n);
	LINK(l)->NEXT = 0;
      }
      /// destructor: de-allocate memory
      ~chunk()
      { WDutils_DEL16(DATA); }
    };// struct pool::chunk
    //@}
    /// \name data
    //@{
  private:
    const size_type N;                    ///< # elements / chunk
    const size_type Kp;                   ///< sizeof(element)
    unsigned        NC;                   ///< # chunks
    size_type       Na, Nmax;             ///< # elements given out
    chunk          *CHUNKS;               ///< pter to 1st chunk
    link           *HEAD;                 ///< pter to 1st free element
    //@}
    /// \name member methods
    //@{
    /// grow: add another chunk
    void grow()
    {
      chunk *c = new chunk(N,Kp);
      c->NEXT  = CHUNKS;
      CHUNKS   = c;
      HEAD     = LINK(CHUNKS->DATA);
      ++NC;
    }
  public:
    /// ctor
    /// \param[in] n desired number of elements
    /// \param[in] k sizeof elements
    pool(size_type n, size_type k)
      : N      ( n<1? 1 : n ),
	Kp     ( sizeof(link)<k? k : sizeof(link) ),
	NC     ( 1 ),
	Na     ( 0 ),
	Nmax   ( 0 ),
	CHUNKS ( new chunk(N,Kp) ),
	HEAD   ( LINK(CHUNKS->DATA) )
    {}
#undef LINK
    /// dtor
    ~pool()
    {
      chunk *a=CHUNKS, *n;
      while(a) {
	n = a->NEXT;
	WDutils_DEL_O(a);
	a = n;
      }
    }
    /// hand out an element = pointer to Kp (=at least K) free bytes
    void *alloc()
    {
      if(HEAD==0) grow();
      link*p = HEAD;
      HEAD = p->NEXT;
      Na++;
      if(Na>Nmax) Nmax=Na;
      return p;
    }
    /// take back an element (freeing)
    /// \param e pointer to element
    /// \note the pointer \a e MUST have been previously allocated from this
    void free(void* e)
    {
      link *p = static_cast<link*>(e);
      p->NEXT = HEAD;
      HEAD    = p;
      Na--;
    }
    /// # chunks used
    unsigned N_chunks() const
    { return NC; }
    /// # elements handed out
    size_type N_alloc() const
    { return Na; }
    /// # elements allocated
    size_type N_alloc_max() const
    { return Nmax; }
  };// class WDutils::pool
  // ///////////////////////////////////////////////////////////////////////////
  template<> struct traits< pool >
  { static const char*name () { return "pool"; } };
  //
  template<> struct traits< pool::chunk >
  { static const char*name () { return "pool::chunk"; } };
  //
  //  class WDutils::Pool<>
  //
  /// template class, based on WDutils::pool, for allocating elements of type T
  //
  template<typename T>
  class Pool : private pool {
    //  no copy ctor and no operator=
    Pool(const Pool&) WDutilsCXX11Delete;
    Pool& operator=(const Pool&) WDutilsCXX11Delete;
  public:
    /// ctor: allocate 1st chunk
    /// \param[in] n number of elements in first chunk
    explicit Pool(size_type n)
      : pool(n,sizeof(T)) {}
    /// allocation: hand out a single element
    T*alloc()
    { return static_cast<T*>(pool::alloc()); }
    /// freeing: take back a single element
    void free (T*e)
    { pool::free(e); }
    using pool::N_chunks;    ///< # chunks used
    using pool::N_alloc;     ///< # elements handed out
    using pool::N_alloc_max; ///< # elements allocated
  };// class WDutils::Pool<>
  //
  template<typename T> struct traits< Pool<T> > {
    static const char  *name ()
    { return message("Pool<%s>",traits<T>::name()); }
  };
  //
  // Arrays<T,D> of type T and arbitrary dimension D 
  //
  template<typename T, unsigned D> class SubArray;
  /// \brief    a const array of arbitrary type in D dimensions
  /// \details  used as a return type to support A[i0][i1][i2][i3] (in this case
  ///           for an Array<T,4>.
  template<typename T, unsigned D> class ConstSubArray {
    WDutilsStaticAssert((D>2));
    friend class SubArray<T,D+1>;
    friend class ConstSubArray<T,D+1>;
    const unsigned*const N; // sizes in all dimensions
    const unsigned*const K; // offsets in all dimensions
    const T       *const A; // pointer to first element
    //  no copy ctor and no operator=
    ConstSubArray           (const ConstSubArray&) WDutilsCXX11Delete;
    ConstSubArray& operator=(const ConstSubArray&) WDutilsCXX11Delete;
    // private constructor: accessible from friends only
    ConstSubArray(const unsigned*n, const unsigned*k, const T*a)
      : N(n), K(k), A(a) {}
  public:
    /// rank: number of dimensions
    static const unsigned rank = D;
    /// return size in dimension @a d
    unsigned size(unsigned d=0) const
    { return N[d]; }
    /// return product of size in all dimensions
    unsigned prod_sizes() const
    { return N[0]*K[0]; }
    /// type resulting from a const [] operation
    typedef ConstSubArray<T,D-1> ConstSub;
    /// acts like the operator[] on a const pointer
    ConstSub operator[] (unsigned i) const
    { return ConstSub(N+1, K+1, A+i*K[0]); }
    /// type conversion to const pointer: return first element
    const T*array() const { return A; }
  };
  /// a const array of arbitrary type in 2 dimensions
  template<typename T> class ConstSubArray<T,2> {
    friend class SubArray<T,3>;
    friend class ConstSubArray<T,3>;
    const unsigned*const N; ///< sizes in all dimensions
    const unsigned*const K; ///< offsets in all dimensions
    const T       *const A; ///< pointer to first element
    //  no copy ctor and no operator=
    ConstSubArray           (const ConstSubArray&) WDutilsCXX11Delete;
    ConstSubArray& operator=(const ConstSubArray&) WDutilsCXX11Delete;
    /// private constructor: accessible from friends only
    ConstSubArray(const unsigned*n, const unsigned*k, const T*a)
      : N(n), K(k), A(a) {}
  public:
    /// rank: number of dimensions
    static const unsigned rank = 2;
    /// return size in dimension @a d, default: in lowest dimension
    unsigned size(unsigned d=0) const
    { return N[d]; }
    /// return product of size in all dimensions
    unsigned prod_sizes() const
    { return N[0]*K[0]; }
    /// type resulting from a const [] operation
    typedef const T* ConstSub;
    /// acts like the operator[] on a const pointer
    ConstSub operator[] (unsigned i) const
    { return A+i*K[0]; }
    /// type conversion to const pointer: return first element
    const T*array() const { return A; }
  };
  //
  /// an array of arbitrary type in D dimensions.
  /// used as base and return type for class Array<T,D>
  template<typename T, unsigned D> class SubArray {
    WDutilsStaticAssert((D>2));
    friend class SubArray<T,D+1>;
  protected:
    const unsigned*const N; // sizes in all dimensions
    const unsigned*const K; // offsets in all dimensions
    T             *const A; // pointer to allocated memory 
    //  no copy ctor and no operator=
    SubArray           (const SubArray&) WDutilsCXX11Delete;
    SubArray& operator=(const SubArray&) WDutilsCXX11Delete;
    // protected constructor: accessible from friends and derived
    SubArray(const unsigned*n, const unsigned*k, T*a) : N(n), K(k), A(a) {}
  public:
    /// rank: number of dimensions
    static const unsigned rank = D;
    /// return size in dimension @a d, default: in lowest dimension
    /// \note We assume 0 <= @a d < @a rank. For efficiency, this is not tested.
    unsigned size(unsigned d=0) const
    { return N[d]; }
    /// return product of size in all dimensions
    unsigned prod_sizes() const
    { return N[0]*K[0]; }
    /// type resulting from a non-const [] operation
    typedef SubArray<T,D-1> Sub;
    /// type resulting from a const [] operation
    typedef ConstSubArray<T,D-1> ConstSub;
    /// acts like the operator[] on a pointer
    /// \param[in] i  index in lowest dimension
    Sub operator[] (unsigned i)
    { return Sub(N+1,K+1,A+i*K[0]); }
    /// acts like the operator[] on a const pointer
    /// \param[in] i  index in lowest dimension
    ConstSub operator[] (unsigned i) const
    { return ConstSub(N+1,K+1,A+i*K[0]); }
    /// type conversion to pointer: return first element
    T      *array()       { return A; }
    /// type conversion to const pointer: return first element
    const T*array() const { return A; }
  };
  /// an array of arbitrary type in 2 dimensions
  template<typename T> class SubArray<T,2> {
    friend class SubArray<T,3>;
  protected:
    const unsigned*const N; // sizes in all dimensions
    const unsigned*const K; // offsets in all dimensions
    T             *const A; // pointer to allocated memory 
    //  no copy ctor and no operator=
    SubArray           (const SubArray&) WDutilsCXX11Delete;
    SubArray& operator=(const SubArray&) WDutilsCXX11Delete;
    // protected constructor: accessible from friends and derived
    SubArray(const unsigned*n, const unsigned*k, T*a) : N(n), K(k), A(a) {}
  public:
    /// rank: number of dimensions
    static const unsigned rank = 2;
    /// return size in dimension @a d, default: in lowest dimension
    /// \note We assume 0 <= @a d < @a rank. For efficiency, this is not tested.
    unsigned size(unsigned d=0) const
    { return N[d]; }
    /// return product of size in all dimensions
    unsigned prod_sizes() const
    { return N[0]*K[0]; }
    /// type resulting from a non-const [] operation
    typedef T* Sub;
    /// type resulting from a const [] operation
    typedef const T* ConstSub;
    /// acts like the operator[] on a pointer
    /// \param[in] i  index in lowest dimension
    Sub operator[] (unsigned i)
    { return A+i*K[0]; }
    /// acts like the operator[] on a const pointer
    /// \param[in] i  index in lowest dimension
    ConstSub operator[] (unsigned i) const
    { return A+i*K[0]; }
    /// type conversion to pointer: return first element
    T      *array()       { return A; }
    /// type conversion to const pointer: return first element
    const T*array() const { return A; }
  };
  //
  // class Array<T,D>
  //
  /// \brief   A @a D dimensional array on the heap.
  /// \details Acts like a @a T[N1]...[ND], but @a N1 ...@a ND can be chosen at
  ///          run-time and only one chunk of memory is ever allocated. The
  ///          operator[] behaves as for an array and is implemented via types
  ///          SubArray and ConstSubArray.
  template<typename T, unsigned D=1> class Array : public SubArray<T,D>
  {
    typedef SubArray<T,D> Base;
    using Base::N;
    using Base::K;
    using Base::A;
    /// \name data
    //@{
    unsigned __N[D];  ///< __N[d]: size in dimension d 
    unsigned __K[D];  ///< __K[d] = Prod_i>d __N[i]      
    //@}
    /// set N[d] and K[d]
    /// \param[in] n size of array in each dimension
    void set(const unsigned*n) {
      for(unsigned d=0; d!=D; ++d)
	__N[d] = n? n[d] : 0;
      __K[D-1] = 1;
      for(unsigned d=D-1; d!=0; --d)
	__K[d-1] = K[d] * N[d];
    }
    /// is a set @a n of sizes equal to ours?
    /// \param[in] n  size of array in each dimension
    bool equal(const unsigned n[D]) const {
      for(unsigned d=0; d!=D; ++d)
	if(N[d] != n[d]) return false;
      return true;
    }
    //  no copy ctor and no operator=
    Array           (const Array&) WDutilsCXX11Delete;
    Array& operator=(const Array&) WDutilsCXX11Delete;
  public:
    using Base::Sub;
    using Base::ConstSub;
    /// default constructor: sizes are all equal to 0
    Array()
      : Base(__N,__K,0) { set(0); }
    /// construction from sizes
    /// \param[in] n size of array in each dimension
    explicit Array(const unsigned n[D]) WDutils_THROWING
      : Base(__N,__K,0) { reset(n); }
    /// construction from sizes for D=2 (compile-time error otherwise)
    /// \param[in] n0 size of array in dimension 0
    /// \param[in] n1 size of array in dimension 1
    Array(unsigned n0, unsigned n1) WDutils_THROWING
    : Base(__N,__K,0)
    {
      WDutilsStaticAssert(D==2);
      const unsigned n[2] = {n0,n1};
      reset(n);
    }
    /// construction from sizes for D=3 (compile-time error otherwise)
    /// \param[in] n0 size of array in dimension 0
    /// \param[in] n1 size of array in dimension 1
    /// \param[in] n2 size of array in dimension 2
    Array(unsigned n0, unsigned n1, unsigned n2) WDutils_THROWING
    : Base(__N,__K,0)
    {
      WDutilsStaticAssert(D==3);
      const unsigned n[3] = {n0,n1,n2};
      reset(n);
    }
    /// construction from sizes for D=4 (compile-time error otherwise)
    /// \param[in] n0 size of array in dimension 0
    /// \param[in] n1 size of array in dimension 1
    /// \param[in] n2 size of array in dimension 2
    /// \param[in] n3 size of array in dimension 3
    Array(unsigned n0, unsigned n1, unsigned n2, unsigned n3) WDutils_THROWING
    : Base(__N,__K,0)
    {
      WDutilsStaticAssert(D==4);
      const unsigned n[4] = {n0,n1,n2,n3};
      reset(n);
    }
    /// construction from sizes for D=5 (compile-time error otherwise)
    /// \param[in] n0 size of array in dimension 0
    /// \param[in] n1 size of array in dimension 1
    /// \param[in] n2 size of array in dimension 2
    /// \param[in] n3 size of array in dimension 3
    /// \param[in] n4 size of array in dimension 4
    Array(unsigned n0, unsigned n1, unsigned n2, unsigned n3, unsigned n4)
    WDutils_THROWING
    : Base(__N,__K,0)
    {
      WDutilsStaticAssert(D==5);
      const unsigned n[5] = {n0,n1,n2,n3,n4};
      reset(n);
    }
    /// construction from sizes and initial value
    /// \param[in] n size of array in each dimension
    /// \param[in] x initialize each element with this value
    Array(const unsigned n[D], T const&x) WDutils_THROWING
    : Base(__N,__K,0) { reset(n,x); }
    /// destruction: de-allocate memory
    ~Array() WDutils_THROWING
    {
      if(A) {
	WDutils_DEL_AN(A,K[0]*N[0]);
	const_cast<T*&>(A) = 0;
      }
      set(0);
    }
    /// set all values to given constant
    /// \param[in] x initialize each element with this value
    void setval(T const&x = T(0) )
    { for(unsigned i=0; i!=K[0]*N[0]; ++i) A[i]=x; }
    /// reset: destruct and construct again
    /// \param[in] n new size of array in each dimension
    void reset(const unsigned n[D]) WDutils_THROWING
    {
      if(A==0 || !equal(n) ) {
	if(A) WDutils_DEL_AN(A,K[0]*N[0]);
	set(n);
	const_cast<T*&>(A) = K[0]*N[0] ? WDutils_NEW(T,K[0]*N[0]) : 0;
      }
    }
    /// reset for D=2 (compile-time error otherwise)
    /// \param[in] n0 new size of array in dimension 0
    /// \param[in] n1 new size of array in dimension 1
    void reset(unsigned n0, unsigned n1) WDutils_THROWING
    {
      WDutilsStaticAssert(D==2);
      const unsigned n[D] = {n0,n1};
      reset(n);
    }
    /// reset for D=3 (compile-time error otherwise)
    /// \param[in] n0 new size of array in dimension 0
    /// \param[in] n1 new size of array in dimension 1
    /// \param[in] n2 new size of array in dimension 2
    void reset(unsigned n0, unsigned n1, unsigned n2) WDutils_THROWING
    {
      WDutilsStaticAssert(D==3);
      const unsigned n[D] = {n0,n1,n2};
      reset(n);
    }
    /// reset for D=4 (compile-time error otherwise)
    /// \param[in] n0 new size of array in dimension 0
    /// \param[in] n1 new size of array in dimension 1
    /// \param[in] n2 new size of array in dimension 2
    /// \param[in] n3 new size of array in dimension 3
    void reset(unsigned n0, unsigned n1, unsigned n2, unsigned n3)
      WDutils_THROWING
    {
      WDutilsStaticAssert(D==4);
      const unsigned n[D] = {n0,n1,n2,n3};
      reset(n);
    }
    /// reset for D=5 (compile-time error otherwise)
    /// \param[in] n0 new size of array in dimension 0
    /// \param[in] n1 new size of array in dimension 1
    /// \param[in] n2 new size of array in dimension 2
    /// \param[in] n3 new size of array in dimension 3
    /// \param[in] n4 new size of array in dimension 4
    void reset(unsigned n0, unsigned n1, unsigned n2, unsigned n3, unsigned n4)
      WDutils_THROWING
    {
      WDutilsStaticAssert(D==5);
      const unsigned n[D] = {n0,n1,n2,n3,n4};
      reset(n);
    }
    /// reset: destruct and construct again
    /// \param[in] n new size of array in each dimension
    /// \param[in] x initialize each element with this value
    void reset(const unsigned n[D], T const&x) WDutils_THROWING {
      reset(n);
      setval(x);
    }
  };// class Array<T,D>
  //
  /// \brief   Specialisation for D=1
  /// \details Acts like a simple T[N], except that memory is on the heap not
  ///          the stack.
  template<typename T> class Array<T,1> {
    /// \name data
    //@{
    unsigned  N;     ///< N[d]: size in dimension d 
    T        *A;     ///< pointer to allocated memory 
    //@}
    //  no copy ctor and no operator=
    Array           (const Array&) WDutilsCXX11Delete;
    Array& operator=(const Array&) WDutilsCXX11Delete;
  public:
#if __cplusplus >= 201103L

    /// move ctor
    Array(Array&&a)
      : N(a.N), A(a.A)
    { a.N=0; a.A=0; }
    /// move operator
    Array& operator=(Array&&a)
    {
      N = a.N; a.N=0;
      A = a.A; a.A=0;
      return*this;
    }
#endif
    /// rank: number of dimensions
    static const unsigned rank = 1;
    /// return size of array
    unsigned const&size() const
    { return N; }
    /// return size of array
    unsigned const&size(unsigned) const
    { return N; }
    /// return product of size in all dimensions
    unsigned prod_sizes() const
    { return N; }
    /// set all values to given constant
    /// \param[in] x initialize each element with this value
    void setval(T const&x = T(0) )
    { for(unsigned i=0; i!=N; ++i) A[i] = x; }
    /// grow: increase size by n, but keep old data
    /// \param[in] n  grow by this much, default: old size ->
    void grow(unsigned n=0) WDutils_THROWING
    {
      n = N + (n?n:N);
      if(n) {
	T*newA = WDutils_NEW(T,n);
	if(N) {
	  memcpy(newA,A,sizeof(T)*N);
	  WDutils_DEL_A(A);
	}
	N = n;
	A = newA;
      }
    }
    /// reset: destruct and construct again
    /// \param[in] n new size of array
    void reset(unsigned n) WDutils_THROWING
    {
      if(n!=N || (n && A==0)) {
	if(A) WDutils_DEL_AN(A,N);
	N = n;
	A = N>0? WDutils_NEW(T,N) : 0;
      }
    }
    /// reset: destruct and construct again
    /// \param[in] n new size of array
    void reset(const unsigned n[1]) WDutils_THROWING
    { reset(n[0]); }
    /// reset: destruct and construct again
    /// \param[in] n new size of array
    /// \param[in] x initial value for each element
    void reset(unsigned n, T const&x) WDutils_THROWING {
      reset(n);
      for(unsigned i=0; i!=N; ++i) A[i] = x;
    }
    /// reset: destruct and construct again
    /// \param[in] n new size of array
    /// \param[in] x initial value for each element
    void reset(const unsigned n[1], T const&x) WDutils_THROWING {
      reset(n[0],x);
    }
    /// default constructor: size equal to 0
    Array() : N(0), A(0) {}
    /// construction from sizes
    /// \param[in] n size of array in each dimension
    explicit Array(unsigned n) WDutils_THROWING
      : N(0), A(0) { reset(n); }
    /// construction from sizes
    /// \param[in] n size of array in each dimension
    explicit Array(const unsigned n[1]) WDutils_THROWING
    : N(0), A(0) { reset(n); }
    /// construction from size and initial value
    /// \param[in] n size of array in each dimension
    /// \param[in] x initialize each element with this value
    Array(const unsigned n, T const&x) WDutils_THROWING
    : N(0), A(0) { reset(n,x); }
    /// construction from size and initial value
    /// \param[in] n size of array in each dimension
    /// \param[in] x initialize each element with this value
    Array(const unsigned n[1], T const&x) WDutils_THROWING
    : N(0), A(0) { reset(n,x); }
    /// destruction: de-allocate memory
    ~Array() WDutils_THROWING { 
      if(A) {
	WDutils_DEL_AN(A,N);
	A = 0;
      }
      N = 0;
    }
    /// type conversion to pointer: return first element
    T      *array()       { return A; }
    /// type conversion to const pointer: return first element
    const T*array() const { return A; }
    /// type resulting from a non-const [] operation
    typedef T& Sub;
    /// type resulting from a const [] operation
    typedef T const& ConstSub;
    /// non-const array sub-scription: return reference to element
    Sub operator[] (unsigned i)
    { return A[i]; }
    /// const array sub-scription: return const reference to element
    ConstSub operator[] (unsigned i) const
    { return A[i]; }
  };// class Array<T,1>
  template<int D> class Array<void,D> {};
  template<> class Array<void,1> {};
  //
  /// \brief  Specialisation for D=0, provided for completeness only
  template<typename T> class Array<T,0> {
    T A;
  public:
    static const unsigned rank = 0;
    unsigned size() const { return 1; }
    unsigned size(unsigned) const { return 1; }
    explicit Array(unsigned*) {}
    Array(unsigned*, T const&x) : A(x) {}
    typedef T& Sub;
    typedef T const& ConstSub;
    Sub      operator[] (unsigned)       { return A; }
    ConstSub operator[] (unsigned) const { return A; }
  };
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T, unsigned D> struct traits< Array<T,D> > {
    static const char*name () {
      return message("Array<%s,%d>",traits<T>::name(),D);
    }
  };
  // ///////////////////////////////////////////////////////////////////////////
  // ///////////////////////////////////////////////////////////////////////////
  class BitArray : private Array<uint64_t>
  {
    typedef Array<uint64_t> Base;
    static const uint64_t null = 0;
    static const uint64_t full = ~null;
    static unsigned rsize(unsigned n)
    { return (n>>6) + (n&full)? 1 : 0; }
    //  no copy ctor and no operator=
    BitArray           (const BitArray&) WDutilsCXX11Delete;
    BitArray& operator=(const BitArray&) WDutilsCXX11Delete;
  public:
    /// return size of array
    unsigned size() const
    { return Base::size()<<6; }
    /// set all bits to either off (default) or on
    /// \param[in] b  initialise each bit to this value
    void setval(bool b=false)
    { Base::setval(b? full:null); }
    /// reset: destruct and construct again
    /// \param[in] n  new size of array
    void reset(unsigned n)
    { Base::reset(rsize(n)); }
    /// reset: destruct and construct again
    /// \param[in] n  new size of array
    /// \param[in] b  initialise each bit to this value
    void reset(unsigned n, bool b)
    { Base::reset(rsize(n),b? full:null); }
    /// default constructor: size equal to 0
    BitArray()
      : Base() {}
    /// construction from size
    /// \param[in] n  size of array
    explicit BitArray(unsigned n) WDutils_THROWING
      : Base(rsize(n)) {}
    /// construction from size and initial value
    /// \param[in] n  size of array
    /// \param[in] b  initialise each bit to this value
    BitArray(unsigned n, bool b) WDutils_THROWING
      : Base(rsize(n),b? full:null) {}
    /// const data access
    /// \param[in] i  index of bit to access
    bool operator[] (unsigned i) const
    { return Base::operator[](i>>6) & (1<<(i&full)); }
    /// switch bit on
    /// \param[in] i  index of bit to switch on
    void on(unsigned i)
    { Base::operator[](i>>6) |= 1<<(i&full); }
    /// switch bit off
    /// \param[in] i  index of bit to switch off
    void off(unsigned i)
    { Base::operator[](i>>6) &= ~(1<<(i&full)); }
    /// set a bit
    /// \param[in] i  index of bit to set
    /// \param[in] b  value to set bit to
    void set(unsigned i, bool b)
    { if(b) on(i); else off(i); }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //
  //  class Stack<X>
  //
  /// a simple stack of elements of type X
  ///
  // ///////////////////////////////////////////////////////////////////////////
  template<typename X>
  class Stack {
    X *S;        ///< begin of stack
    X *P;        ///< top of stack
    X *const SN; ///< end of stack
    //  no copy ctor and no operator=
    Stack           (const Stack&) WDutilsCXX11Delete;
    Stack& operator=(const Stack&) WDutilsCXX11Delete;
  public:
    /// ctor: allocate memory, empty stack
    explicit
    Stack(unsigned n) : S(n? WDutils_NEW(X,n) : 0), P(S), SN(S+n) {}
    /// ctor: allocate memory, put one element on stack
    Stack(unsigned n, const X&a) : S(WDutils_NEW(X,n? n:1)), P(S), SN(S+n)
    { push(a); }
    /// dtor: de-allocate memory
    ~Stack() { if(S) WDutils_DEL_AN(S,unsigned(SN-S)); S=0; }
    /// is stack empty?
    bool is_empty () const { return P<=S; }
    /// is there space for more to stack?
    bool has_space() const { return P>=SN; }
    /// push another X onto the stack
    /// \note we use the operator=(X,X), which must be defined and accessible
    void push(const X&a) WDutils_THROWING {
      if(P>=SN) WDutils_THROW("Stack<%s>::push(): exceeding stack\n",nameof(X));
      *(P++) = a;
    }
    /// pop top element off the stack
    X pop() { return *(--P); }
    /// return top element, but don't pop it
    X&top() { return *P; }
    /// empty stack
    void reset() { P=S; }
    /// capacity of stack
    unsigned capacity() const
    { return unsigned(SN-S); }
    /// reset stack capacity
    void reset_capacity(unsigned new_capacity)
    {
      unsigned old_size = (S && P>S) ? unsigned(P-S) : 0;
      if(new_capacity < old_size)
	new_capacity = old_size;
      X *newS = WDutils_NEW(X,new_capacity);
      if(old_size)
	std::memcpy(newS,S,sizeof(X)*old_size);
      if(S)
	WDutils_DEL_AN(S,unsigned(SN-S));
      S = newS;
      P = S+old_size;
      SN= S+new_capacity;
    }
  };
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T> struct traits< Stack<T> > {
    static const char  *name () {
      return message("Stack<%s>",traits<T>::name());
    }
  };
  ///
  /// allocate and de-allocate aligned memory
  /// \ingroup MemAligned
  ///
  template<std::size_t alignment>
  struct static_allocator {
    WDutilsStaticAssert(alignment>1);
    static void*allocate(std::size_t n)
    {
      if(n == 0) return 0;
      DebugInfoN(WDutilsAllocDebugLevel+2,
		 "static_allocator<%lu>: trying to allocate %lu bytes\n",
		 alignment,n);
      if(n > max_size())
	throw std::bad_alloc();
      void*ret =
#if defined(__GNUC__) || defined (__INTEL_COMPILER)
	_mm_malloc
#else
	_aligned_malloc
#endif
	(n,alignment);
      if(!ret)
	throw std::bad_alloc();
      DebugInfoN(WDutilsAllocDebugLevel,
		 "static_allocator<%lu>: allocated %lu bytes @ %p\n",
		 alignment,n,ret);
      return ret;
    }
    static void deallocate(void*p)
    {
      DebugInfoN(WDutilsAllocDebugLevel+2,
		 "static_allocator<%lu>: trying to de-allocated memory @ %p\n",
		 alignment,p);
#if defined(__GNUC__) || defined (__INTEL_COMPILER)
      _mm_free
#else
      _aligned_free
#endif
	(p);
      DebugInfoN(WDutilsAllocDebugLevel,
		 "static_allocator<%lu>: de-allocated memory @ %p\n",
		 alignment,p);
    }
    static std::size_t max_size () noexcept
    { return std::numeric_limits<std::size_t>::max(); }
  };
  ///
  /// allocate and de-allocate unaligned memory
  /// \ingroup MemAligned
  ///
  template<>
  struct static_allocator<1> {
    static std::size_t max_size () noexcept
    { return std::numeric_limits<std::size_t>::max(); }
    static void*allocate(std::size_t n)
    { 
      if(n == 0) return 0;
      DebugInfoN(WDutilsAllocDebugLevel+2,
		 "static_allocator<1>: trying to allocate %lu bytes\n",n);
      void*ret = new char[n];
      DebugInfoN(WDutilsAllocDebugLevel,
		 "static_allocator<1>: allocated %lu bytes @ %p\n",n,ret);
      return ret;
    }
    static void deallocate(void*p)
    { 
      DebugInfoN(WDutilsAllocDebugLevel+2,
		 "static_allocator<1>: trying to de-allocated memory @ %p\n",p);
      delete[] static_cast<char*>(p);
      DebugInfoN(WDutilsAllocDebugLevel,
		 "static_allocator<1>: de-allocated memory @ %p\n",p);
    }
  };
  template<> struct static_allocator<0>;
  ///
  /// allocator with explicit alignment
  /// \ingroup MemAligned
  ///
  template<typename _Tp, std::size_t alignment = 16>
  class AlignmentAllocator
  {
    WDutilsStaticAssert(alignment);
    typedef static_allocator<alignment> static_alloc;
  public:
    typedef size_t     size_type;
    typedef ptrdiff_t  difference_type;
    typedef _Tp*       pointer;
    typedef const _Tp* const_pointer;
    typedef _Tp&       reference;
    typedef const _Tp& const_reference;
    typedef _Tp        value_type;

    template <typename _Tp1>
    struct rebind
    { typedef AlignmentAllocator<_Tp1, alignment> other; };

    AlignmentAllocator() noexcept {}

    AlignmentAllocator (const AlignmentAllocator&) noexcept {}

    template <typename _Tp1>
    AlignmentAllocator (const AlignmentAllocator<_Tp1, alignment> &) noexcept {}

    ~AlignmentAllocator () noexcept {}

    pointer address (reference __x) const noexcept
    {
#if __cplusplus >= 201103L
      return std::addressof(__x);
#else
      return reinterpret_cast<_Tp*>(&reinterpret_cast<char&>(__x));
#endif
    }

    const_pointer address (const_reference __x) const noexcept
    {
#if __cplusplus >= 201103L
      return std::addressof(__x);
#else
      return reinterpret_cast<const _Tp*>(&reinterpret_cast<const char&>(__x));
#endif
    }

    pointer allocate (size_type __n, const void* = 0)
    {
      return
	static_cast<pointer>(static_alloc::allocate(__n*sizeof(value_type)));
    }

    void deallocate (pointer __p, size_type)
    { static_alloc::deallocate(__p); }

    size_type max_size () const noexcept
    { return static_alloc::max_size() / sizeof (value_type); }

#if __cplusplus >= 201103L

    template<typename _Up, typename... _Args>
    void construct(_Up* __p, _Args&&... __args)
    { ::new(static_cast<void*>(__p)) _Up(std::forward<_Args>(__args)...); }
    
    template<typename _Up>
    void destroy(_Up* __p)
    { __p->~_Up(); }

#else

    void construct (pointer __p, const_reference __val)
    { ::new(static_cast<void*>(__p)) value_type(__val); }

    void destroy (pointer __p)
    { __p->~value_type (); }

#endif

    bool operator!=(const AlignmentAllocator&) const 
    { return false; }
  
    // Returns true if and only if storage allocated from *this
    // can be deallocated from other, and vice versa.
    // Always returns true for stateless allocators.
    bool operator==(const AlignmentAllocator&) const 
    { return true; }

  };// class AlignmentAllocator<>

  /// AlignmentAllocator<void> specialization.
  /// \ingroup MemAligned
  template<std::size_t alignment>
  class AlignmentAllocator<void, alignment>
  {
  public:
    typedef size_t      size_type;
    typedef ptrdiff_t   difference_type;
    typedef void*       pointer;
    typedef const void* const_pointer;
    typedef void        value_type;

    template<typename _Tp1>
    struct rebind
    { typedef AlignmentAllocator<_Tp1, alignment> other; };
  };
  ///
  /// managing raw memory
  /// \ingroup MemAligned
  ///
  template<std::size_t alignment = 0>
  class raw_memory
  {
    typedef static_allocator<alignment> static_alloc;
    char       *_mem;
    std::size_t _num;
  public:
    /// default ctor: no memory
    raw_memory() noexcept
      : _mem(0)
      , _num(0) {}
#if __cplusplus >= 201103L
    /// move ctor: steal memory
    raw_memory(raw_memory&&other) noexcept
      : _mem(other._mem)
      , _num(other._num)
    { other._mem = 0; other._num = 0; }
    /// move operator: steal memory
    raw_memory&operator=(raw_memory&&other) noexcept
    {
      _mem = other._mem;
      _num = other._num;
      other._mem = 0;
      other._num = 0;
      return*this;
    }
#endif// c++11
    /// ctor: allocate @a n bytes
    explicit raw_memory(std::size_t n)
      : _mem(static_cast<char*>(static_alloc::allocate(n)))
      , _num(n)
    {}
    /// dtor
    ~raw_memory()
    { if(_mem) static_alloc::deallocate(_mem); }
    /// release memory
    char*release() noexcept
    {
      char*tmp = _mem;
      _mem = 0;
      _num = 0;
      return tmp;
    }
    /// reset capacity
    void reset(std::size_t n)
    {
      if(n != _num) {
	if(_mem) static_alloc::deallocate(_mem);
	_num = n;
	_mem = _num? static_cast<char*>(static_alloc::allocate(_num)) : 0;
      }
    }
    /// reset capacity, but shrink only if  n*q < p*capacity()
    template<unsigned q, unsigned p>
    void reset_conditional(std::size_t n)
    {
      static_assert(q>=p,"shrink not allowed");
      if(n>_num || n*q < p*_num) reset(n);
    }
    /// number of bytes allocated
    std::size_t capacity() const noexcept
    { return _num; }
    /// access to raw memory
    const char*get() const noexcept
    { return _mem; }
    /// access to raw memory
    char*get() noexcept
    { return _mem; }
    /// begin of raw memory
    const char*begin() const noexcept
    { return _mem; }
    /// begin of raw memory
    char*begin() noexcept
    { return _mem; }
    /// end of raw memory
    const char*end() const noexcept
    { return _mem + _num; }
    /// end of raw memory
    char*end() noexcept
    { return _mem + _num; }
    //  disable copy from lvalue
#if __cplusplus >= 201103L
    raw_memory(raw_memory const&) = delete;
    raw_memory&operator=(raw_memory const&) = delete;
#else // c++11
  private:
    raw_memory(raw_memory const&);
    raw_memory&operator=(raw_memory const&);
#endif// c++11
  };
} // namespace WDutils {
////////////////////////////////////////////////////////////////////////////////
#if __cplusplus < 201103L
# undef noexcept
# undef constexpr
#endif
#endif // WDutils_included_memory_h

// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/memory.h
///
/// \author Walter Dehnen
///                                                                             
/// \date   2000-2008
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2008 Walter Dehnen
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
#ifndef WDutils_included_iostream
#  include <iostream>
#  define WDutils_included_iostream
#endif
#ifndef WDutils_included_exception_h
#  include <exception.h>
#endif
#ifndef WDutils_included_traits_h
# include <traits.h>
#endif
#ifndef WDutils_included_inline_h
#  include <inline.h>
#endif

namespace WDutils {
#ifndef WDutilsAllocDebugLevel
#define WDutilsAllocDebugLevel 8
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
  /// \return    a valid pointer (unless an error occurs) 
  /// \param[in] n number of array elements
  /// \param[in] f name of the source file where this routines is called
  /// \param[in] l number of the line in that file 
  /// \param[in] lib (optional) name of calling library (default: "WDutils")
  template<typename T> inline
  T* NewArray(size_t n, const char*f, int l, const char*lib = "WDutils")
    WDutils_THROWING
  {
    T*t;
    try {
      t = new T[n];
    } catch(std::bad_alloc E) {
      t = 0;
#ifdef WDutils_EXCEPTIONS
      throw Thrower
#else
      Error
#endif
	(f,l)("caught std::bad_alloc\n");
    }
    if(debug(WDutilsAllocDebugLevel))
      DebugInformation(f,l,lib)("allocated %lu %s = %lu bytes @ %p\n",
				n,nameof(T),n*sizeof(T),static_cast<void*>(t));
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
#define WDutils_NEW(TYPE,SIZE) WDutils::NewArray<TYPE>(SIZE,__FILE__,__LINE__)

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
  /// \param[in] a  pointer previously allocated with WDutils::NewArray<>()
  ///               or ::operator new[].
  /// \param[in] f  name of the source file where this routines is called
  /// \param[in] l  number of the line in that file
  /// \param[in] lib (optional) name of calling library (default: "WDutils")
  template<typename T> inline
  void DelArray(T* a, const char*f, int l, const char*lib = "WDutils")
    WDutils_THROWING {
    if(0==a) {
      Warning(f,l)("trying to delete zero pointer to array of '%s'", nameof(T));
      return;
    }
    try {
      delete[] a;
    } catch(...) {
#ifdef WDutils_EXCEPTIONS
      throw Thrower
#else
      Error
#endif
	(f,l)("de-allocating array of '%s' @ %p' failed\n", nameof(T),a);
    }
    if(debug(WDutilsAllocDebugLevel))
      DebugInformation(f,l,lib)("de-allocated array of %s @ %p\n",
				nameof(T), static_cast<void*>(a));
  }
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T> inline
  void DelArray(const T* a, const char*f, int l, const char*lib = "WDutils")
    WDutils_THROWING {
    if(0==a) {
      Warning(f,l)("trying to delete zero pointer to array of '%s'", nameof(T));
      return;
    }
    try {
      delete[] a;
    } catch(...) {
#ifdef WDutils_EXCEPTIONS
      throw Thrower
#else
      Error
#endif
	(f,l)("de-allocating array of '%s' @ %p' failed\n", nameof(T),a);
    }
    if(debug(WDutilsAllocDebugLevel))
      DebugInformation(f,l,lib)("de-allocated array of %s @ %p\n",
				nameof(T), static_cast<const void*>(a));
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
#define WDutils_DEL_A(P) WDutils::DelArray(P,__FILE__,__LINE__)
  // ///////////////////////////////////////////////////////////////////////////
  //
  /// Object de-allocation giving useful info in case of error; mostly used
  /// from macro WDutils_DEL_O.
  ///
  /// In case of a de-allocation error (if the pointer provided was not valid)
  /// an error is generated (or an exception thrown, depending on the WDutils
  /// error settings).
  ///
  /// \param[in] a  pointer previously allocated with ::operator new().
  /// \param[in] f  name of the source file where this routines is called
  /// \param[in] l  number of the line in that file
  /// \param[in] lib (optional) name of calling library (default: "WDutils")
  template<typename T> inline
  void DelObject(T* a, const char*f, int l, const char*lib="WDutils")
    WDutils_THROWING {
    if(0==a) {
      Warning(f,l)("trying to delete zero pointer to object '%s'", nameof(T));
      return;
    }
    try {
      delete a;
    } catch(...) {
#ifdef WDutils_EXCEPTIONS
      throw Thrower
#else
      Error
#endif
	(f,l)("de-allocating object '%s' @ %p failed\n", nameof(T),a);
    }
    if(debug(WDutilsAllocDebugLevel))
      DebugInformation(f,l,lib)("de-allocated %s object @ %p\n",
				nameof(T), static_cast<void*>(a));
  }
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T> inline
  void DelObject(const T* a, const char*f, int l, const char*lib="WDutils")
    WDutils_THROWING {
    if(0==a) {
      Warning(f,l)("trying to delete zero pointer to object '%s'", nameof(T));
      return;
    }
    try {
      delete a;
    } catch(...) {
#ifdef WDutils_EXCEPTIONS
      throw Thrower
#else
      Error
#endif
	(f,l)("de-allocating object '%s' @ %p failed\n",nameof(T),a);
    }
    if(debug(WDutilsAllocDebugLevel))
      DebugInformation(f,l,lib)("de-allocated %s object @ %p\n",
				nameof(T), static_cast<const void*>(a));
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
#define WDutils_DEL_O(P) WDutils::DelObject(P,__FILE__,__LINE__)
  // ///////////////////////////////////////////////////////////////////////////
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
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T>
  class block_alloc {
  public:
    ////////////////////////////////////////////////////////////////////////////
    //
    /// \name some public typedefs
    //@{
    ////////////////////////////////////////////////////////////////////////////
    typedef T         value_type;                ///< type of elements          
    typedef size_t    size_type;                 ///< type of number of elements
    typedef ptrdiff_t difference_type;           ///< type of pointer difference
    typedef T*        pointer;                   ///< type of pointer to element
    typedef const T*  const_pointer;             ///< type of const pointer     
    typedef T&        reference;                 ///< type of reference to elem 
    typedef const T&  const_reference;           ///< type of const reference   
    //@}
    // /////////////////////////////////////////////////////////////////////////
    //
    //  sub-type WDutils::block_alloc::block
    //
    /// allocates and manages a contiguous chunk of elements                    
    ///
    // /////////////////////////////////////////////////////////////////////////
  private:
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
      //------------------------------------------------------------------------
    public:
      /// constructor
      /// \param n number of elements to allocate
      explicit
      block(size_type const&n) :
	NEXT    ( 0 ),
	FIRST   ( WDutils_NEW(value_type,n) ),
	END     ( FIRST ),
	ENDTOT  ( FIRST + n ) {}
      /// destructor: de-allocate
      ~block() {
	WDutils_DEL_A(FIRST);
      }
      /// link block to next block
      void link(block* next) {
	NEXT = next;
      }
      /// give out: another element
      pointer new_element() {
	return END++;
      }
      /// give out: N new elements
      pointer new_elements(size_type n) {
	register pointer OLD_END=END;
	END += n;
	return OLD_END;
      }
      /// can we give out up to \a n elements?
      bool has_free(size_type n) const {
	return END + n <= ENDTOT;
      }
      /// is block full (no free elements)?
      bool is_full() const {
	return END>=ENDTOT;
      }
      /// return next block in linked list of blocks
      block  *next () const {
	return NEXT;
      }
      /// return first element allocated
      pointer front() const {
	return FIRST;
      }
      /// return last element used
      pointer back () const {
	return END-1;
      }
      /// return end of elements used
      pointer end  () const {
	return END;
      }
      /// return number of elements used
      difference_type N_used () const {
	return END-FIRST;
      }
      /// return number of elemnets allocated
      difference_type N_alloc() const {
	return ENDTOT-FIRST;
      }
    };
  public:
    // /////////////////////////////////////////////////////////////////////////
    //
    //  sub-type WDutils::block_alloc::iterator                                 
    //
    /// for forward sequential iteration through all elements used              
    ///
    // /////////////////////////////////////////////////////////////////////////
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
	: B(b),   E(e)   {}
      /// copy constructor
      iterator(iterator const&I)
	: B(I.B), E(I.E) {}
      /// copy assignment
      iterator& operator= (iterator const&I) {
	B=I.B;
	E=I.E;
	return *this;
      }
      //@}
      /// prefix ++: increment pointer; if out of block: get next block
      iterator& operator++() {
	++E;
	if(E == B->end()) {
	  B = B->next();
	  E = B? B->front() : 0;
	}
	return *this;
      }
      /// postfix ++: return temporary equal original, but increment this
      iterator  operator++(int) {
	iterator tmp(*this);
	this->operator++();
	return tmp;
      }
      /// equality: refer to the same element?
      bool operator==(const iterator&I) {
	return B == I.B  &&  E == I.E;
      }
      /// inequality: not equal
      bool operator!=(const iterator&I) {
	return ! operator==(I);
      }
      /// \name conversions to pointer to elment or reference to element
      //@{
      /// type conversion to pointer to element
      operator pointer() {
	return E;
      }
      /// type conversion to pointer to const element
      operator const_pointer() const {
	return E;
      }
      /// reference operator
      reference operator* () {
	return*E;
      }
      /// const reference operator
      const_reference operator* () const {
	return*E;
      }
      /// dereference operator
      pointer operator->() {
	return E;
      }
      /// const dereference operator
      const_pointer operator->() const {
	return E;
      }
      //@}
    };
  private:
    // /////////////////////////////////////////////////////////////////////////
    //
    /// \name member data of class WDutils::block_alloc
    //@{
    // /////////////////////////////////////////////////////////////////////////
    block       *FIRST;  ///< first block in linked list
    block       *LAST;   ///< last block in linked list
    size_type    NTOT;   ///< # elements allocated
    size_type    NUSED;  ///< # elements used
    //@}
  public:
    // /////////////////////////////////////////////////////////////////////////
    //
    /// \name member functions of class WDutils::block_alloc
    //@{
    // /////////////////////////////////////////////////////////////////////////
    /// constructor: allocate first block
    /// \param[in] Ns number of elements in 1st block
    explicit
    block_alloc(size_type Ns)
      : FIRST ( new block(Ns) ),
	LAST  ( FIRST ),
	NTOT  ( Ns ),
	NUSED ( 0 )  {}
    /// destructor: delete all blocks
    ~block_alloc() WDutils_THROWING;
    /// give out: another element
    /// \return pointer to allocated element
    pointer new_element() {
      if(LAST->is_full()) {
	size_type New = LAST->N_alloc();
	LAST->link(new block(New));
	LAST = LAST->next();
	NTOT+= New;
      }
      ++NUSED;
      return LAST->new_element();
    }
    /// give out: another element
    /// \return pointer to allocated element
    /// \param F estimator = function object returning N_required(N_used yet)
    template<class estimator>
    pointer new_element(estimator const&F) {
      if(LAST->is_full()) {
	size_type New = F(NUSED);
	LAST->link(new block(New));
	LAST = LAST->next();
	NTOT+= New;
      }
      ++NUSED;
      return LAST->new_element();
    }
    /// give out: more elements
    /// \return pointer to \a Ne allocated elements
    /// \param Ne number of elements to return pointer to
    pointer new_elements(size_type Ne) {
      if(!LAST->has_free(Ne)) {
	size_type New=max(Ne,static_cast<size_type>(LAST->N_alloc()));
	LAST->link(new block(New));
	LAST = LAST->next();
	NTOT+= New;
      }
      NUSED+=Ne;
      return LAST->new_elements(Ne);
    }
    /// give out: more elements
    /// \return pointer to \a Ne allocated elements
    /// \param Ne number of elements to return pointer to
    /// \param F estimator = function object returning N_required(N_used yet)
    template<class estimator>
    pointer new_elements(size_type       Ne,
			 estimator const&F) {
      if(!LAST->has_free(Ne)) {
	size_type New=max(Ne,F(NUSED));
	LAST->link(new block(New));
	LAST = LAST->next();
	NTOT+= New;
      }
      NUSED+=Ne;
      return LAST->new_elements(Ne);
    }
    /// the running number of a given element
    ///
    /// if \a E does appear not to point to an element given out by this, we
    /// return a negative number
    /// \return the running number of element pointed to by \a E
    /// \param[in]  E  pointer to element
    int number_of_element(const_pointer E) const {
      if(E==0) return -99;
      int num=0;
      for(block*B=FIRST; B; B=B->next())
	if(E >= B->front() && E < B->end())
	  return num + int(E-B->front());
	else
	  num += B->N_used();
      return -99;
    }
    /// was a given element given out by this?
    /// \return true if a given pointer to element was given out by this
    /// \param E pointer to element
    bool is_element(const_pointer E) const {
      if(E==0) return 0;
      for(register block*B=FIRST; B; B=B->next())
	if(E >= B->front() && E < B->end())
	  return 1;
      return 0;
    }
    /// return pointer to first element
    pointer first() const {
      return FIRST->front();
    }
    /// return iterator referring to first element
    iterator front() const {
      return iterator(FIRST,FIRST->front());
    }
    /// return iterator referring to first element
    iterator begin() const {
      return iterator(FIRST,FIRST->front());
    }
    /// return iterator referring to end of elements
    iterator end() const {
      return iterator(0,0);
    }
    /// return number of elements given out sofar
    size_type N_used() const {
      return NUSED;
    }
    /// return number of elements allocated sofar
    size_type N_allocated() const {
      return NTOT;
    }
  };// class WDutils::block_alloc
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T> struct traits< block_alloc<T> > {
    static const char  *name () {
      return message("block_alloc<%s>",traits<T>::name());
    }
  };
  // ///////////////////////////////////////////////////////////////////////////
  // does not compile with gcc 4.3.1, which seems a compiler bug
#if defined(__GNUC__) && ( __GNUC__ < 4 || __GNUC_MINOR__ < 3)
  template<typename T> struct traits< typename block_alloc<T>::block > {
    static const char  *name () {
      return message("block_alloc<%s>::block",traits<T>::name());
    }
  };
#endif
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T> inline
  block_alloc<T>::~block_alloc() WDutils_THROWING {
    block *A=FIRST, *N;
    while(A) {
      N = A->next();
      WDutils_DEL_O(A);
      A = N;
    }
  }
  // ///////////////////////////////////////////////////////////////////////////
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
  // ///////////////////////////////////////////////////////////////////////////
  class pool {
  public:
    // /////////////////////////////////////////////////////////////////////////
    //
    /// \name sub-types of class WDutils::pool
    //@{
    ////////////////////////////////////////////////////////////////////////////
    typedef size_t    size_type;         ///< type of number of elements
    typedef ptrdiff_t difference_type;   ///< type of pointer difference
  public:
    /// elementary of a linked list
    struct link {
      link *NEXT;   ///< pter to next link
    };
#define LINK(NAME) static_cast<link*>(static_cast<void*>(NAME))
    /// a chunk of elements
    struct chunk {
      char   *DATA; ///< pter to allocated memory
      chunk  *NEXT; ///< pter to next chunk
      /// constructor
      /// \param[in] N  number of element in chunk
      /// \param[in] Kp sizeof(elements)
      chunk(size_type N, size_type Kp)
	: DATA ( WDutils_NEW(char,N*Kp) ),
	  NEXT ( 0 ) {
	const    char *END=DATA+N*Kp;
	register char *l,*n;
	for(l=DATA, n=DATA+Kp; n!=END; l+=Kp,n+=Kp)
	  LINK(l)->NEXT = LINK(n);
	LINK(l)->NEXT = 0;
      }
      /// destructor: de-allocate memory
      ~chunk() {
	WDutils_DEL_A(DATA);
      }
    };// struct pool::chunk
    //@}
    // /////////////////////////////////////////////////////////////////////////
    //
    /// \name data of class WDutils::pool
    //@{
    // /////////////////////////////////////////////////////////////////////////
  private:
    const size_type N;                    ///< # elements / chunk
    const size_type Kp;                   ///< sizeof(element)
    unsigned        NC;                   ///< # chunks
    unsigned        Na, Nmax;             ///< # elements given out
    chunk          *CHUNKS;               ///< pter to 1st chunk
    link           *HEAD;                 ///< pter to 1st free element
    //@}
    // /////////////////////////////////////////////////////////////////////////
    //
    /// \name member functions of class WDutils::pool
    //@{
    // /////////////////////////////////////////////////////////////////////////
    /// grow: add another chunk
    void grow() {
      chunk *c = new chunk(N,Kp);
      c->NEXT  = CHUNKS;
      CHUNKS   = c;
      HEAD     = LINK(CHUNKS->DATA);
      ++NC;
    }
  public:
    /// construction
    /// \param[in] n desired number of elements
    /// \param[in] k sizeof elements
    pool(size_type n, size_type k)
      : N      ( n<1? 1 : n ),
	Kp     ( sizeof(link)<k? k : sizeof(link) ),
	NC     ( 1u ),
	Na     ( 0u ),
	Nmax   ( 0u ),
	CHUNKS ( new chunk(N,Kp) ),
	HEAD   ( LINK(CHUNKS->DATA) )
	  {}
#undef LINK
    /// destruction: delete all chunks
    ~pool() {
      register chunk *a=CHUNKS, *n;
      while(a) {
	n = a->NEXT;
	WDutils_DEL_O(a);
	a = n;
      }
    }
    /// hand out an element = pointer to Kp (=at least K) free bytes
    void *alloc() {
      if(HEAD==0) grow();
      register link*p = HEAD;
      HEAD = p->NEXT;
      Na++;
      if(Na>Nmax) Nmax=Na;
      return p;
    }
    /// take back an element (freeing)
    /// \param e pointer to element
    /// \note the pointer \a e MUST have been previously allocated from this
    void free(void* e) {
      register link *p = static_cast<link*>(e);
      p->NEXT = HEAD;
      HEAD    = p;
      Na--;
    }
    /// return number of chunks used
    unsigned const&N_chunks   () const {
      return NC;
    }
    /// return number of bytes handed out
    unsigned const&N_alloc    () const {
      return Na;
    }
    /// return number of bytes allocated
    unsigned const&N_alloc_max() const {
      return Nmax;
    }
  };// class WDutils::pool
  // ///////////////////////////////////////////////////////////////////////////
  template<> struct traits< pool > {
    static const char    *name () { return "pool"; }
  };
#if(0) // gcc v 3.3.5 is buggy and doesn't like this
  template<> struct traits< pool::chunk > {
    static const char    *name () { return "pool::chunk"; }
  };
#endif
  // ///////////////////////////////////////////////////////////////////////////
  //
  //  class WDutils::Pool<>
  //
  /// template class, based on WDutils::pool, for allocating elemnts of type T
  //
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T>
  class Pool : private pool {
  public:
    /// constructor: allocate 1st chunk
    /// \param[in] n number of elements in first chunk
    explicit
    Pool(size_type n)
      : pool(n,sizeof(T)) {}
    /// allocation: hand out a single element
    T*   alloc()    {
      return static_cast<T*>(pool::alloc());
    }
    /// freeing: take back a single element
    void free (T*e) {
      pool::free(e);
    }
    pool::N_chunks;    ///< return number of chunks used
    pool::N_alloc;     ///< return number of bytes handed out
    pool::N_alloc_max; ///< return number of bytes allocated
  };// class WDutils::Pool<>
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T> struct traits< Pool<T> > {
    static const char  *name () {
      return message("Pool<%s>",traits<T>::name());
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //
  // Arrays<T,D> of type T and arbitrary dimension D 
  //
  //////////////////////////////////////////////////////////////////////////////
#if defined(DEBUG) || defined(EBUG)
#  define THROW_BAD WDutils_THROWING
#  define CHECK_BAD(M,R)						\
  if(i< 0) {								\
    WDutils_Warning("i=%d <  0 in Array subscript #%d",i,R);		\
    WDutils_THROW("i=%d <  0 in Array subscript #%d",i,R);		\
  }									\
  if(i>=M) {								\
    WDutils_Warning("i=%d >= N=%d in Array subscript #%d",i,M,R);	\
    WDutils_THROW("i=%d >= N=%d in Array subscript #%d",i,M,R);		\
  }
#else
#  define THROW_BAD
#  define CHECK_BAD(M,R)
#endif
  //////////////////////////////////////////////////////////////////////////////
  template<typename, unsigned> class ConstPseudoArray;
  template<typename, unsigned> class PseudoArray;
  template<typename, unsigned> class Array;
  // ///////////////////////////////////////////////////////////////////////////
  //
  //  class ConstPseudoArray<T,D>
  //
  /// used as return type for Array<>::operator[] const
  ///
  /// Apart from a function returning the size in the first dimension, the only
  /// public member functions for D>1 are the operator[] and element(), which
  /// both return a ConstPseudoArray<T,D-1>.  For D=1, the operator[] returns a
  /// 'T const&'. Moreover for D=1, there is a type conversion to const T*.
  // 
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T, unsigned D> class ConstPseudoArray {
    friend class ConstPseudoArray<T,D+1>;
    friend class PseudoArray<T,D+1>;
    friend class Array<T,D+1>;
    /// \name data
    //@{
    const unsigned*const N; ///< sizes in all dimensions
    const unsigned*const K; ///< offsets in all dimensions
    const T       *const A; ///< pointer to first element
    //@}
    /// constructor is private: accessible only from friends
    ConstPseudoArray(const T*a, const unsigned*n, const unsigned*k)
      : N(n), K(k), A(a) {}
  public:
    /// return size in 1st dimension
    unsigned const& size() const {
      return N[0];
    }
    /// acts like the operator[] on a pointer const T*...*
    ConstPseudoArray<T,D-1> operator[] (unsigned i) const THROW_BAD {
      CHECK_BAD(N[0],D);
      return ConstPseudoArray<T,D-1>(A+i*K[0],N+1,K+1);
    }
    /// same as operator[]
    ConstPseudoArray<T,D-1> element (unsigned i) const THROW_BAD {
      CHECK_BAD(N[0],D);
      return ConstPseudoArray<T,D-1>(A+i*K[0],N+1,K+1);
    }
  };// class ConstPseudoArray<T,D>
  // ///////////////////////////////////////////////////////////////////////////
  /// special case D=1
  template<typename T> class ConstPseudoArray<T,1> {
    friend class ConstPseudoArray<T,2>;
    friend class PseudoArray<T,2>;
    friend class Array<T,2>;
    /// \name data
    //@{
    const unsigned*const N; ///< N[0] = length
    const T       *const A; ///< pointer to first element
    //@}
    /// constructor is private: accessible only from friends
    ConstPseudoArray(const T*a, const unsigned*n, const unsigned*)
      : N(n), A(a) {}
  public:
    /// return size in 1st dimension
    unsigned const& size() const {
      return N[0];
    }
    /// conversion to pointer to constant element
    operator const T* () const {
      return A;
    }
    /// acts like the operator[] on a pointer const T*
    T const & operator[](unsigned i) const THROW_BAD {
      CHECK_BAD(N[0],1);
      return A[i];
    }
    /// same as operator[]
    T const & element(unsigned i) const THROW_BAD {
      CHECK_BAD(N[0],1);
      return A[i];
    }
  };// class ConstPseudoArray<T,D=1>
  // ///////////////////////////////////////////////////////////////////////////
  //
  //  class PseudoArray<T,D>
  //
  /// used as return type for Array<>::operator[]                               
  ///
  /// Apart from a function returning the size in the first dimension, the only
  /// public member functions for D>1 are the operator[], which returns a
  /// PseudoArray<T,D-1>, and the operator[] const, which returns a
  /// ConstPseudoArray<T,D-1>.\n
  /// For D=1, the operator[] returns a 'T&' and the operator[] const a 'T
  /// const&'. Moreover for D=1, there are type conversions to const T* and T*.
  //
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T, unsigned D> class PseudoArray {
    friend class PseudoArray<T,D+1>;
    friend class Array<T,D+1>;
    /// \name data
    //@{
    const unsigned*const N; ///< sizes in all dimensions
    const unsigned*const K; ///< offsets in all dimensions
    T             *const A; ///< pointer to first element
    //@}
    /// constructor is private: accessible only from friends
    PseudoArray(T*a, const unsigned*n, const unsigned*k) : N(n), K(k), A(a) {}
  public:
    /// return size in 1st dimension
    unsigned const& size() const {
      return N[0];
    }
    /// acts like the operator[] on a pointer T*...*
    PseudoArray<T,D-1> operator[] (unsigned i) THROW_BAD {
      CHECK_BAD(N[0],D);
      return PseudoArray<T,D-1>(A+i*K[0],N+1,K+1);
    }
    /// acts like the operator[] on a pointer const T*...*
    ConstPseudoArray<T,D-1> operator[] (unsigned i) const THROW_BAD {
      CHECK_BAD(N[0],D);
      return ConstPseudoArray<T,D-1>(A+i*K[0],K+1);
    }
    /// same as operator[]
    PseudoArray<T,D-1> element (unsigned i) THROW_BAD {
      CHECK_BAD(N[0],D);
      return PseudoArray<T,D-1>(A+i*K[0],N+1,K+1);
    }
    /// same as operator[]
    ConstPseudoArray<T,D-1> element (unsigned i) const THROW_BAD {
      CHECK_BAD(N[0],D);
      return ConstPseudoArray<T,D-1>(A+i*K[0],K+1);
    }
  };// class PseudoArray<T,D>
  // ///////////////////////////////////////////////////////////////////////////
  /// special case D=1
  template<typename T> class PseudoArray<T,1> {
    friend class PseudoArray<T,2>;
    friend class Array<T,2>;
    /// \name data
    //@{
    const unsigned*const N; ///< sizes in all dimensions
    T             *const A; ///< pointer to first element
    //@}
    /// constructor is private: accessible only from friends
    PseudoArray(T*a, const unsigned*n, const unsigned*) : N(n), A(a) {}
  public:
    /// return size in 1st dimension
    unsigned const& size() const {
      return N[0];
    }
    /// conversion to pointer to element
    operator T* () {
      return A;
    }
    /// conversion to pointer to constant element
    operator const T* () const {
      return A;
    }
    /// acts like the operator[] on a pointer T*
    T       & operator[](unsigned i)       THROW_BAD {
      CHECK_BAD(N[0],1);
      return A[i];
    }
    /// acts like the operator[] on a pointer const T*
    T const & operator[](unsigned i) const THROW_BAD {
      CHECK_BAD(N[0],1);
      return A[i];
    }
    /// same as operator[]
    T       & element(unsigned i)       THROW_BAD {
      CHECK_BAD(N[0],1);
      return A[i];
    }
    /// same as operator[]
    T const & element(unsigned i) const THROW_BAD {
      CHECK_BAD(N[0],1);
      return A[i];
    }
  };// class PseudoArray<T,D=1>
  // ///////////////////////////////////////////////////////////////////////////
  //
  //  class Array<T,D>
  ///
  /// a D-dimensional array, operator[] as expected, specialisation for D=1.
  /// acts like a T[N_1]...[N_D], but allocates only one chunk of memory and no
  /// pointers.
  //
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T, unsigned D=1> class Array {
    /// \name data
    //@{
    unsigned N[D];     ///< N[d]: size in dimension d 
    unsigned K[D];     ///< K[d] = Prod_i>d N[i]      
    union {
      T      *A;       ///< pointer to allocated memory 
      const T*C;       ///< const pointer to allocated memory 
    };
    //@}
    /// set N[d] and K[d]
    /// \param[in] n size of array in each dimension
    void set(const unsigned*n) {
      for(unsigned d=0; d!=D; ++d)
	N[d] = n? n[d] : 0;
      K[D-1] = 1;
      for(unsigned d=D-1; d!=0; --d)
	K[d-1] = K[d] * N[d];
    }
    /// is a set \a n of sizes equal to ours?
    /// \param[in] n size of array in each dimension
    bool equal(const unsigned n[D]) const {
      for(unsigned d=0; d!=D; ++d)
	if(N[d] != n[d]) return false;
      return true;
    }
    /// copy constructor is disabled (private); use references instead of copies
    Array(Array const&);
  public:
    /// rank: number of dimensions
    static const unsigned rank = D;
    /// type resulting from a non-const [] operation
    typedef PseudoArray<T,D-1> Sub;
    /// type resulting from a const [] operation
    typedef ConstPseudoArray<T,D-1> ConstSub;
    /// return size in dimension \a d
    /// \param[in] d dimension to return size for
    unsigned const&size(unsigned d) const {
      return N[d];
    }
    /// set all values to given constant
    /// \param[in] x initialize each element with this value
    void setval(T const&x = T(0) ) WDutils_THROWING {
      for(unsigned i=0; i!=K[0]*N[0]; ++i) A[i] = x;
    }
    /// reset: destruct and construct again
    /// \param[in] n new size of array in each dimension
    void reset(const unsigned n[D]) WDutils_THROWING {
      if(A==0 || !equal(n) ) {
	if(A) WDutils_DEL_A(A);
	set(n);
	A = K[0]*N[0] ? WDutils_NEW(T,K[0]*N[0]) : 0;
      }
    }
    /// reset: destruct and construct again
    /// \param[in] n new size of array in each dimension
    /// \param[in] x initialize each element with this value
    void reset(const unsigned n[D], T const&x) WDutils_THROWING {
      reset(n);
      setval(x);
    }
    /// default ctor. construction from nothing: sizes are all equal to 0
    Array() : A(0) {
      set(0);
    }
    /// construction from sizes
    /// \param[in] n size of array in each dimension
    explicit
    Array(const unsigned n[D]) WDutils_THROWING
    : A(0) {
      reset(n);
    }
    /// construction from sizes for D=2
    /// \param[in] n0 size of array in dimension 0
    /// \param[in] n1 size of array in dimension 1
    Array(unsigned n0, unsigned n1) WDutils_THROWING
    : A(0) {
      WDutilsStaticAssert(D==2);
      const unsigned n[2] = {n0,n1};
      reset(n);
    }
    /// construction from sizes for D=3
    /// \param[in] n0 size of array in dimension 0
    /// \param[in] n1 size of array in dimension 1
    /// \param[in] n2 size of array in dimension 2
    Array(unsigned n0, unsigned n1, unsigned n2) WDutils_THROWING
    : A(0) {
      WDutilsStaticAssert(D==3);
      const unsigned n[3] = {n0,n1,n2};
      reset(n);
    }
    /// construction from sizes for D=4
    /// \param[in] n0 size of array in dimension 0
    /// \param[in] n1 size of array in dimension 1
    /// \param[in] n2 size of array in dimension 2
    /// \param[in] n3 size of array in dimension 3
    Array(unsigned n0, unsigned n1, unsigned n2, unsigned n3) WDutils_THROWING
    : A(0) {
      WDutilsStaticAssert(D==4);
      const unsigned n[4] = {n0,n1,n2,n3};
      reset(n);
    }
    /// construction from sizes for D=5
    /// \param[in] n0 size of array in dimension 0
    /// \param[in] n1 size of array in dimension 1
    /// \param[in] n2 size of array in dimension 2
    /// \param[in] n3 size of array in dimension 3
    /// \param[in] n4 size of array in dimension 4
    Array(unsigned n0, unsigned n1, unsigned n2, unsigned n3, unsigned n4)
    WDutils_THROWING
    : A(0) {
      WDutilsStaticAssert(D==5);
      const unsigned n[5] = {n0,n1,n2,n3,n4};
      reset(n);
    }
    /// construction from sizes and initial value
    /// \param[in] n size of array in each dimension
    /// \param[in] x initialize each element with this value
    Array(const unsigned n[D], T const&x) WDutils_THROWING : A(0) {
      reset(n,x);
    }
    /// destruction: de-allocate memory
    ~Array() WDutils_THROWING {
      if(A) {
	WDutils_DEL_A(A);
	A = 0;
      }
      set(0);
    }
    /// type conversion to pointer: return first element
    /// \return pointer to allocated memory
    T      *const&array()       { return A; }
    /// type conversion to const pointer: return first element
    /// \return const pointer to allocated memory
    const T*const&array() const { return C; }
    /// non-const array sub-scription: return PseudoArray
    /// \param[in] i index in first dimension (dimension 0)
    /// \return sub-array of D-1 dimensions
    PseudoArray<T,D-1> operator[] (unsigned i) THROW_BAD {
      CHECK_BAD(N[0],D);
      return PseudoArray<T,D-1>(A+i*K[0],N+1,K+1);
    }
    /// const array sub-scription: return ConstPseudoArray
    /// \param[in] i index in first dimension (dimension 0)
    /// \return sub-array of D-1 dimensions
    ConstPseudoArray<T,D-1> operator[] (unsigned i) const THROW_BAD {
      CHECK_BAD(N[0],D);
      return ConstPseudoArray<T,D-1>(C+i*K[0],N+1,K+1);
    }
    /// same as operator[]
    /// \param[in] i index in first dimension (dimension 0)
    /// \return sub-array of D-1 dimensions
    PseudoArray<T,D-1> element (unsigned i) THROW_BAD {
      CHECK_BAD(N[0],D);
      return PseudoArray<T,D-1>(A+i*K[0],N+1,K+1);
    }
    /// same as operator[]
    /// \param[in] i index in first dimension (dimension 0)
    /// \return sub-array of D-1 dimensions
    ConstPseudoArray<T,D-1> element (unsigned i) const THROW_BAD {
      CHECK_BAD(N[0],D);
      return ConstPseudoArray<T,D-1>(C+i*K[0],N+1,K+1);
    }
  };// class Array<T,D>
  // ///////////////////////////////////////////////////////////////////////////
  /// special case D=1 for Array<T,D>.
  /// acts like a simple T[N], except that memory is on the heap not the stack.
  template<typename T> class Array<T,1> {
    /// \name data
    //@{
    unsigned N;        ///< N[d]: size in dimension d 
    union {
      T      *A;  ///< pointer to allocated memory 
      const T*C;  ///< const pointer to allocated memory 
    };
    //@}
    /// copy constructor is disabled (private); use references instead of copies
    Array(Array const&) : N(0), A(0) {}
  public:
    /// rank: number of dimensions
    static const unsigned rank = 1;
    /// type resulting from a non-const [] operation
    typedef T& Sub;
    /// type resulting from a const [] operation
    typedef T const& ConstSub;
    /// return size of array
    unsigned const&size() const {
      return N;
    }
    /// return size of array
    unsigned const&size(unsigned) const {
      return N;
    }
    /// set all values to given constant
    /// \param[in] x initialize each element with this value
    void setval(T const&x = T(0) ) WDutils_THROWING {
      for(unsigned i=0; i!=N; ++i) A[i] = x;
    }
    /// reset: destruct and construct again
    /// \param[in] n new size of array
    void reset(unsigned n) WDutils_THROWING {
      if(n!=N || (n && A==0)) {
	if(A) WDutils_DEL_A(A);
	N = n;
	A = N>0? WDutils_NEW(T,N) : 0;
      }
    }
    /// reset: destruct and construct again
    /// \param[in] n new size of array
    void reset(const unsigned n[1]) WDutils_THROWING {
      reset(n[0]);
    }
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
    Array() : N(0), A(0){}
    /// construction from sizes
    /// \param[in] n size of array in each dimension
    explicit
    Array(unsigned n) WDutils_THROWING
    : N(0), A(0) {
      reset(n);
    }
    /// construction from sizes
    /// \param[in] n size of array in each dimension
    explicit
    Array(const unsigned n[1]) WDutils_THROWING
    : N(0), A(0) {
      reset(n);
    }
    /// construction from size and initial value
    /// \param[in] n size of array in each dimension
    /// \param[in] x initialize each element with this value
    Array(const unsigned n, T const&x) WDutils_THROWING
    : N(0), A(0) {
      reset(n,x);
    }
    /// construction from size and initial value
    /// \param[in] n size of array in each dimension
    /// \param[in] x initialize each element with this value
    Array(const unsigned n[1], T const&x) WDutils_THROWING
    : N(0), A(0) {
      reset(n,x);
    }
    /// destruction: de-allocate memory
    ~Array() WDutils_THROWING { 
      if(A) {
	WDutils_DEL_A(A);
	A = 0;
      }
      N = 0;
    }
    /// type conversion to pointer: return first element
    T      *const&array()       { return A; }
    /// type conversion to const pointer: return first element
    const T*const&array() const { return C; }
    /// non-const array sub-scription: return reference to element
    T      & operator[] (unsigned i)       THROW_BAD {
      CHECK_BAD(N,1);
      return A[i];
    }
    /// const array sub-scription: return const reference to element
    T const& operator[] (unsigned i) const THROW_BAD {
      CHECK_BAD(N,1);
      return C[i];
    }
    /// same as operator[]
    T      & element (unsigned i)       THROW_BAD {
      CHECK_BAD(N,1);
      return A[i];
    }
    /// same as operator[]
    T const& element (unsigned i) const THROW_BAD {
      CHECK_BAD(N,1);
      return C[i];
    }
  };// class Array<T,1>
  // ///////////////////////////////////////////////////////////////////////////
  /// special case D=0, a scalar, provided for completeness only
  template<typename T> class Array<T,0> {
    T A;
  public:
    static const unsigned rank = 0;
    unsigned const&size(unsigned) const { return 1; }
    explicit Array(unsigned*) {}
    Array(unsigned*, T const&x) : A(x) {}
    //--------------------------------------------------------------------------
    T      & operator[] (unsigned i)       THROW_BAD {
      CHECK_BAD(0,0);
      return A;
    }
    T const& operator[] (unsigned i) const THROW_BAD {
      CHECK_BAD(0,0);
      return A;
    }
    T      & element (unsigned i)       THROW_BAD {
      CHECK_BAD(0,0);
      return A;
    }
    T const& element (unsigned i) const THROW_BAD {
      CHECK_BAD(0,0);
      return A;
    }
  };
#undef THROW_BAD
#undef CHECK_BAD
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T, unsigned D> struct traits< Array<T,D> > {
    static const char  *name () {
      return message("Array<%s,%d>",traits<T>::name(),D);
    }
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
  public:
    /// ctor: allocate memory, empty stack
    Stack(unsigned n) : S(n? WDutils_NEW(X,n) : 0), P(S), SN(S+n) {}
    /// dtor: de-allocate memory
    ~Stack() { if(S) WDutils_DEL_A(S); S=0; }
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
  };
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T> struct traits< Stack<T> > {
    static const char  *name () {
      return message("Stack<%s>",traits<T>::name());
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  /// \defgroup  Mem16  memory alignment to 16 bytes
  /// \name memory alignment to 16 bytes
  //@{

  //----------------------------------------------------------------------------
  /// Macro enforcing memory alignment to 16 bytes
  /// \ingroup Mem16
  ///
  /// Forces the corresponding variable/type to be 16-byte aligned; Works with
  /// icc (icpc) and gcc (g++) [versions > 3]; Use it like \code 
  ///    struct falcON__align16 name { ... };              \endcode
#if defined (__INTEL_COMPILER)
#  define falcON__align16 __declspec(align(16)) 
#elif defined (__GNUC__) && __GNUC__ > 2
#  define falcON__align16 __attribute__ ((aligned(16)))
#else
#  define falcON__align16
#endif
  //----------------------------------------------------------------------------
  /// is a given memory address aligned?
  /// \ingroup Mem16
  /// \param p  memory address to be tested
  /// \param al alignemt to a bytes will be tested
  inline bool is_aligned(const void*p, int al)
  {
    return size_t(p) % al == 0;
  }
  //----------------------------------------------------------------------------
  /// is a given memory address aligned to a 16 bytes memory location?
  /// \param p  memory address to be tested
  inline bool is_aligned16(const void*p)
  {
    return size_t(p) % 16 == 0;
  }
  //----------------------------------------------------------------------------
  /// Allocate memory at a address aligned to a 16 byte memory location
  /// \ingroup Mem16
  ///
  /// Will allocate slightly more memory than required to ensure we can find
  /// the required amount at a 16b memory location. To de-allocate, you \b must
  /// use falcON::free16(), otherwise an error will occur!
  ///
  /// \return   a newly allocated memory address at a 16 byte memory location
  /// \param n  number of bytes to allocate
  /// \version  debugged 02-09-2004 WD
  inline void* malloc16(size_t n) WDutils_THROWING
  {
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
    char *p = WDutils_NEW(char,n+16+sizeof(void*)); // alloc: (n+16)b + pter    
    char *q = p + sizeof(void*);                    // go sizeof pointer up     
    size_t off = size_t(q) % 16;                    // offset from 16b alignment
    if(off) q += 16-off;                            // IF offset, shift         
    *((void**)(q-sizeof(void*))) = p;               // remember allocation point
    return static_cast<void*>(q);                   // return aligned address   
  }
  //----------------------------------------------------------------------------
  ///
  /// de-allocate memory previously allocated with WDutils::malloc16()
  /// \ingroup Mem16
  ///
  /// This routine \b must be used to properly de-allocate memory that has been
  /// previously allocated by WDutils::malloc16(); other de-allocation will
  /// inevitably result in a run-time \b error!
  ///
  /// \param q  pointer previously allocated by WDutils::malloc16()
  inline void  free16  (void*q) WDutils_THROWING
  {
    WDutils_DEL_A( (char*)( *( (void**) ( ( (char*)q )-sizeof(void*) ) ) ) );
  }
  //----------------------------------------------------------------------------
  ///
  /// allocate memory at a address alignged to at a 16 byte memory location
  /// \ingroup Mem16
  ///
  /// Will allocate slightly more memory than required to ensure we can find
  /// the required amount at a 16b memory location. To de-allocate, you \b must
  /// use WDutils::delete16<>(), otherwise an error will occur!
  ///
  /// \return   a newly allocated memory address at a 16 byte memory location
  /// \param n  number of objects to allocate
  template<typename T> inline T* new16(size_t n) WDutils_THROWING
  {
    return static_cast<T*>(malloc16(n * sizeof(T)));
  }
  //----------------------------------------------------------------------------
  ///
  /// de-allocate memory previously allocated with WDutils::new16().
  /// \ingroup Mem16
  ///
  /// This routine \b must be used to properly de-allocate memory that has been
  /// previously allocated by WDutils::new16(); other de-allocation will
  /// inevitably result in an error.
  ///
  /// \param q  pointer previously allocated by WDutils::new16()
  template<typename T> inline void delete16(T* q)WDutils_THROWING 
  {
    free16(static_cast<void*>(q));
  }
  //@}
  //////////////////////////////////////////////////////////////////////////////
} // namespace WDutils {
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_memory_h

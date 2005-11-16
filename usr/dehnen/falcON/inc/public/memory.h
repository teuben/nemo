// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// memory.h                                                                    |
//                                                                             |
// Copyright (C) 2000-2005  Walter Dehnen                                      |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// container/allocater classes of my own design                                |
//                                                                             |
// class block_alloc<>  - elements are allocated in large blocks.              |
//                      - elements can given out either singly or in chunks    |
//                        of a few                                             |
//                      - a forward iterator (not random access) is provided   |
//                                                                             |
// class pool           - allocates blocks of K bytes in chunks of N elements  |
//                      - single elements can be allocated and freed           |
//                      - chunks are organized as linked list                  |
//                      - free elements are kept in a linked list              |
//                                                                             |
// simple routines for allocating/deallocating 1,2,3D arrays                   |
//                                                                             |
// Array<T,D>           - arrays of type T and arbitrary dimension D           |
//                      - operator [] acts as on pointer with D *s             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_memory_h
#define falcON_included_memory_h

#ifndef falcON_included_cstddef
#  include <cstddef>
#  define falcON_included_cstddef
#endif

#ifndef falcON_included_basic_h
#  include <public/basic.h>
#endif
#ifndef falcON_included_inline_h
#  include <public/inline.h>
#endif

namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::block_alloc<T>                                             //
  //                                                                          //
  // - allocates elements of type T in blocks                                 //
  // - gives elements out either singly or in small amounts                   //
  // - provides a forward iterator for sequential access to all elements      //
  //   given out sofar                                                        //
  //                                                                          //
  // In order to (i) guarantee the correct amount of elements (not known a    //
  // priori) will be allocated and (ii) not to waste memory by allocating too //
  // many, we use the following strategy, similar to that of std::string.     //
  //                                                                          //
  // Elements are allocated in blocks, which in turn are organized in a       //
  // linked list. When a number of new elements is to be allocated and the    //
  // last block in the list cannot provide them, we allocated a new block     //
  // and add it to the list. The number of elements allocated in this new     //
  // block is taken to be a function of the total number of elements used     //
  // sofar. This function must be provided by the user, otherwise the same    //
  // number as in the last block is used.                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename T>
  class block_alloc {
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // 1. member types of class falcON::block_alloc                           //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // 1.0 some public typedefs                                               //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
  public:
    typedef T         value_type;                  // type of elements          
    typedef size_t    size_type;                   // type of number of elements
    typedef ptrdiff_t difference_type;             // type of pointer difference
    typedef T*        pointer;                     // type of pointer to element
    typedef const T*  const_pointer;               // type of const pointer     
    typedef T&        reference;                   // type of reference to elem 
    typedef const T&  const_reference;             // type of const reference   
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // 1.1 class falcON::block_alloc::block                                   //
    //                                                                        //
    //     - allocates a large block of elements                              //
    //     - gives elements out either singly or in small blocks              //
    //     - provides pointers to first and last element given out sofar      //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
  private:
    class block {
    private:
      block      *NEXT;                            // next block in linked list 
      value_type *FIRST;                           // front element             
      value_type *END;                             // end of active elements    
      value_type *ENDTOT;                          // end of all elements       
      block();                                     // not implemented           
      //------------------------------------------------------------------------
    public:
      explicit
      block(size_type const&n) :                   // constructor               
	NEXT    ( 0 ),                             //   no next block           
	FIRST   ( falcON_NEW(value_type,n) ),      //   first of n elements     
	END     ( FIRST ),                         //   no elements used yet    
	ENDTOT  ( FIRST + n ) {}                   //   end of all elements     
      //------------------------------------------------------------------------
      ~block() { falcON_DEL_A(FIRST); }            // destructor: de-allocate   
      //------------------------------------------------------------------------
      void link(block* next){ NEXT = next; }       // link to next block        
      //------------------------------------------------------------------------
      pointer new_element() { return END++; }      // give out: another element 
      //------------------------------------------------------------------------
      pointer new_elements(size_type n) {          // give out: N new elements  
	register pointer OLD_END=END;              //   reserve boxes           
	END += n;                                  //   increment next free     
	return OLD_END;                            //   return old end          
      }
      //------------------------------------------------------------------------
      bool has_free(size_type n) const {           // Do we have N free spaces ?
	return END + n <= ENDTOT;                  //   true if END+N <= ENDTOT 
      }
      //------------------------------------------------------------------------
      bool is_full() const { return END>=ENDTOT; } // is block full?            
      //------------------------------------------------------------------------
      block  *next () const { return NEXT; }       // next block in linked list 
      pointer front() const { return FIRST; }      // first element             
      pointer back () const { return END-1; }      // last element used         
      pointer end  () const { return END; }        // beyond last element used  
      //------------------------------------------------------------------------
      difference_type N_used () const { return END-FIRST; }    // # given out   
      difference_type N_alloc() const { return ENDTOT-FIRST;}  // # allocated   
    };
  public:
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // 1.2 class falcON::block_alloc::iterator                                //
    //                                                                        //
    //     for forward sequential iteration through all boxes                 //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    class iterator{
    private:
      block  *B;                                   // current block             
      pointer E;                                   // current element           
    public:
      iterator(block*b, pointer e)   : B(b),   E(e)   {}
      iterator(iterator const&I)     : B(I.B), E(I.E) {}
      iterator& operator= (iterator const&I) { B=I.B; E=I.E; return *this; }
      //------------------------------------------------------------------------
      iterator& operator++() {                     // prefix ++                 
	++E;                                       //   increment pointer       
	if(E == B->end()) {                        //   IF(out of block)        
	  B = B->next();                           //     get next block in list
	  E = B? B->front() : 0;                   //     get its first element 
	}                                          //   ENDIF                   
	return *this;                              //   return *this            
      }
      //------------------------------------------------------------------------
      iterator  operator++(int) {                  // postfix ++                
	iterator tmp(*this);                       //   temp copy *this         
	this->operator++();                        //   prefix ++ *this         
	return tmp;                                //   return temp copy        
      }
      //------------------------------------------------------------------------
      bool operator==(const iterator&I) { return B == I.B  &&  E == I.E; }
      bool operator!=(const iterator&I) { return B != I.B  ||  E != I.E; }
      //------------------------------------------------------------------------
      operator        pointer   ()       { return E; }
      operator  const_pointer   () const { return E; }
      reference       operator* ()       { return*E; }
      const_reference operator* () const { return*E; }
      pointer         operator->()       { return E; }
      const_pointer   operator->() const { return E; }
    };
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // 2 member data of class falcON::block_alloc                             //
    //                                                                        //
    //   hold & manage the blocks, keep track of elements used and allocated  //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
  private:
    block       *FIRST;                            // first block in linked list
    block       *LAST;                             // last block in linked list 
    size_type    NTOT;                             // # elements allocated      
    size_type    NUSED;                            // # elements used           
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // 3 member functions of class falcON::block_alloc                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
  public:
    explicit
    block_alloc(                                   // constructor:              
		size_type Ns) :                    // I: # elements in 1st block
      FIRST ( new block(Ns) ),                     //   allocate first block    
      LAST  ( FIRST ),                             //   last=first block        
      NTOT  ( Ns ),                                //   # elements allocated    
      NUSED ( 0 )  {}                              //   # elements used sofar   
    //--------------------------------------------------------------------------
    ~block_alloc() {                               // destructor:               
      register block *A=FIRST, *N;                 //   actual & next block     
      while(A) {                                   //   WHILE(actual is valid)  
	N = A->next();                             //     get next block        
	delete A;                                  //     delete actual block   
	A = N;                                     //     set actual = next     
      }                                            //   END WHILE               
    }
    //--------------------------------------------------------------------------
    pointer new_element() {                        // give out: another element 
      if(LAST->is_full()) {                        //   IF(last block is full)  
	register size_type New = LAST->N_alloc();  //     # elements to allocate
	LAST->link(new block(New));                //     allocate new block &  
	LAST = LAST->next();                       //     add it to linked list 
	NTOT+= New;                                //     update # allocated    
      }                                            //   ENDIF                   
      ++NUSED;                                     //   # elements used         
      return LAST->new_element();                  //   return new element      
    }
    //--------------------------------------------------------------------------
    template<class estimator>
    pointer new_element(                           // give out: another element 
			estimator const&F) {       // I: class returning an     
                                                   //    estimate for # elements
                                                   //    still needed given     
                                                   //    # elements used sofar  
      if(LAST->is_full()) {                        //   IF(last block is full)  
	register size_type New = F(NUSED);         //     # elements to allocate
	LAST->link(new block(New));                //     allocate new block &  
	LAST = LAST->next();                       //     add it to linked list 
	NTOT+= New;                                //     update # allocated    
      }                                            //   ENDIF                   
      ++NUSED;                                     //   # elements used         
      return LAST->new_element();                  //   return new element      
    }
    //--------------------------------------------------------------------------
    pointer new_elements(                          // give out: more elements   
			 size_type Ne) {           // I: # elements wanted      
      if(!LAST->has_free(Ne)) {                    //   IF(last block is full)  
	register size_type
	  New=max(Ne,LAST->N_alloc());             //     # elements to allocate
	LAST->link(new block(New));                //     allocate new block &  
	LAST = LAST->next();                       //     add it to linked list 
	NTOT+= New;                                //     update # allocated    
      }                                            //   ENDIF                   
      NUSED+=Ne;                                   //   # elements used         
      return LAST->new_elements(Ne);               //   return new elements     
    }
    //--------------------------------------------------------------------------
    template<class estimator>
    pointer new_elements(                          // give out: more elements   
			 size_type       Ne,       // I: # elements wanted      
			 estimator const&F) {      // I: class returning an     
                                                   //    estimate for # elements
                                                   //    still needed given     
                                                   //    # elements used sofar  
      if(!LAST->has_free(Ne)) {                    //   IF(last block is full)  
	register size_type New=max(Ne,F(NUSED));   //     # elements to allocate
	LAST->link(new block(New));                //     allocate new block &  
	LAST = LAST->next();                       //     add it to linked list 
	NTOT+= New;                                //     update # allocated    
      }                                            //   ENDIF                   
      NUSED+=Ne;                                   //   # elements used         
      return LAST->new_elements(Ne);               //   return new elements     
    }
    //--------------------------------------------------------------------------
    int number_of_element(                         // R: number of given element
			  const_pointer E) const { // I: element                
      if(E==0) return -99;                         //   IF null: return -99     
      register int num=0;                          //   reset counter num       
      for(register block*B=FIRST; B; B=B->next())  //   LOOP blocks             
	if(E >= B->front() && E < B->end())        //     IF in block           
	  return num + int(E-B->front());          //       return num + n(in b)
	else                                       //     ELSE                  
	  num += B->N_used();                      //       add to counter      
      return -99;                                  //   not found: return -99   
    }
    //--------------------------------------------------------------------------
    bool is_element(                               // R: element from this?     
		    const_pointer E) const {       // I: element                
      if(E==0) return 0;                           //   IF null: return false   
      for(register block*B=FIRST; B; B=B->next())  //   LOOP blocks             
	if(E >= B->front() && E < B->end())        //     IF in block           
	  return 1;                                //       return true         
      return 0;                                    //   return false            
    }
    //--------------------------------------------------------------------------
    pointer   first      () const { return FIRST->front(); }
    iterator  front      () const { return iterator(FIRST,FIRST->front()); }
    iterator  begin      () const { return iterator(FIRST,FIRST->front()); }
    iterator  end        () const { return iterator(0,0); }
    size_type N_used     () const { return NUSED; }
    size_type N_allocated() const { return NTOT; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::pool                                                       //
  //                                                                          //
  // - allocates blocks of K bytes (=elements) in chunks of N elements        //
  // - single elements can be allocated and freed                             //
  // - chunks are organized as linked list                                    //
  // - free elements are kept in a linked list                                //
  //                                                                          //
  // The actual number of bytes allocated per element is                      //
  //                                                                          //
  //   K' = max( K, sizeof(*void) ).                                          //
  //                                                                          //
  // Thus, for K < sizeof(void*) this code becomes inefficient and should     //
  // not be used.                                                             //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class pool {
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // 1. member types of class falcON::pool                                  //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // 1.0 some public typedefs                                               //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
  public:
    typedef size_t    size_type;                   // type of number of elements
    typedef ptrdiff_t difference_type;             // type of pointer difference
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // 1.1 struct falcON::pool::link                                          //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
  private:
    struct link { link *NEXT; };                   // pter to next link         
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // 1.2 struct falcON::pool::chunk                                         //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
#define LINK(NAME) static_cast<link*>(static_cast<void*>(NAME))
    struct chunk {
      char   *DATA;                                // pter to allocated memory  
      chunk  *NEXT;                                // pter to next chunk        
      //------------------------------------------------------------------------
      chunk(size_type N,                           // I: # elements             
	    size_type Kp) :                        // I: size of elements       
	DATA ( falcON_NEW(char,N*Kp) ),            // allocate memory           
	NEXT ( 0 )                                 // reset pter to next chunk  
      {
	const    char *END=DATA+N*Kp;              // beyond last byte          
	register char *l,*n;                       // now we link the elements  
	for(l=DATA, n=DATA+Kp; n!=END; l+=Kp,n+=Kp)// LOOP elems of Kp bytes    
	  LINK(l)->NEXT = LINK(n);                 //   link to next one        
	LINK(l)->NEXT = 0;                         // last one: null link       
      }
      //------------------------------------------------------------------------
      ~chunk() { falcON_DEL_A(DATA); }             // de-allocate memory        
    };
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // 2 member data of class falcON::pool                                    //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
  private:
    const size_type N;                             // # elements / chunk        
    const size_type Kp;                            // sizeof(element)           
    unsigned        NC;                            // # chunks                  
    unsigned        Na, Nmax;                      // # elements given out      
    chunk          *CHUNKS;                        // pter to 1st chunk         
    link           *HEAD;                          // pter to 1st free element  
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // 3 member functions of class falcON::pool                               //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
    void grow() {                                  // add another chunk         
      chunk *c = new chunk(N,Kp);                  //   allocate new chunk      
      c->NEXT  = CHUNKS;                           //   add to list of chunks   
      CHUNKS   = c;                                //   make it first chunk     
      HEAD     = LINK(CHUNKS->DATA);               //   add elements to list    
      ++NC;                                        //   increment chunks counter
    }
    //--------------------------------------------------------------------------
    // constructor                                                              
    //--------------------------------------------------------------------------
  public:
    pool(size_type n,                              // I: desired # elements     
	 size_type k) :                            // I: desired sizeof(elem)   
      N      ( n<1? 1 : n ),                       //  actual # elements        
      Kp     ( sizeof(link)<k? k : sizeof(link) ), //  actual sizeof(element)   
      CHUNKS ( new chunk(N,Kp) ),                  //  allocate 1st chunk       
      HEAD   ( LINK(CHUNKS->DATA) ),               //  get list of elements     
      NC     ( 1u ),                               //  we start with 1 chunk    
      Na     ( 0u ),                               //  reset counter            
      Nmax   ( 0u ) {}                             //  reset counter            
#undef LINK
    //--------------------------------------------------------------------------
    // destructor                                                               
    //--------------------------------------------------------------------------
    ~pool() {
      register chunk *a=CHUNKS, *n;                //  actual & next chunk      
      while(a) {                                   //  WHILE(actual is valid)   
	n = a->NEXT;                               //    get next chunk         
	delete a;                                  //    delete actual chunk    
	a = n;                                     //    set actual = next      
      }                                            //  END WHILE                
    }
    //--------------------------------------------------------------------------
    // hand out an element, i.e. a pter to Kp (=at least K) free bytes          
    //--------------------------------------------------------------------------
    void *alloc() {
      if(HEAD==0) grow();                          // no more elements: grow    
      register link*p = HEAD;                      // get 1st free element      
      HEAD   = p->NEXT;                            // reset pter to ----        
      Na++;                                        // increment counter         
      if(Na>Nmax) Nmax=Na;                         // update maximum            
      return p;                                    // return free element       
    }
    //--------------------------------------------------------------------------
    // take back an element                                                     
    //--------------------------------------------------------------------------
    void free(void* e) {                           // I: pter to element        
      register link *p = static_cast<link*>(e);    // cast into our format      
      p->NEXT = HEAD;                              // add to list of free       
      HEAD    = p;                                 // elements                  
      Na--;                                        // decrement counter         
    }
    //--------------------------------------------------------------------------
    // inform on the number of chunks                                           
    //--------------------------------------------------------------------------
    unsigned const&N_chunks   () const { return NC; }
    unsigned const&N_alloc    () const { return Na; }
    unsigned const&N_alloc_max() const { return Nmax; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::Pool<>                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename T>
  class Pool : private pool {
  public:
    explicit
    Pool(size_type n) : pool(n,sizeof(T)) {}
    T*   alloc()    { return static_cast<T*>(pool::alloc()); }
    void free (T*e) { pool::free(e); }
    pool::N_chunks;
    pool::N_alloc;
    pool::N_alloc_max;
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // simple methods for 1, 2 & 3D array allocation/deallocation               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template <typename T>
  inline void Alloc1D(T* &A, int N) falcON_THROWING {
    A = falcON_NEW(T,N); 
  }
  //----------------------------------------------------------------------------
  template <typename T>
  inline void Alloc1D(T* &A, int N, T const&X) falcON_THROWING {
    A = falcON_NEW(T,N); 
    for(T a=A; a != A+N; ++a) *a = X;
  }
  //----------------------------------------------------------------------------
  template <typename T>
  inline void Free1D(T* A) falcON_THROWING {
    falcON_DEL_A(A);
  }
  //----------------------------------------------------------------------------
  template <typename T>
  inline void Alloc2D(T** &A, int N[2]) falcON_THROWING {
    register int i, iN1;
    A    = falcON_NEW(T*,N[0]);
    A[0] = falcON_NEW(T ,N[0]*N[1]);
    for(int i=0,iN1=0; i!=N[0]; ++i, iN1+=N[1])
      A[i] = A[0] + iN1;
  }
  //----------------------------------------------------------------------------
  template <typename T>
  inline void Alloc2D(T** &A, int N[2], T const&X) falcON_THROWING {
    register int i, iN1;
    A    = falcON_NEW(T*,N[0]);
    A[0] = falcON_NEW(T ,N[0]*N[1]);
    for(int i=0,iN1=0; i!=N[0]; ++i, iN1+=N[1]) {
      A[i] = A[0] + iN1;
      for(int j=0; j!=N[1]; ++j)
	A[i][j] = X;
    }
  }
  //----------------------------------------------------------------------------
  template <typename T>
  inline void Free2D(T** A) falcON_THROWING {
    falcON_DEL_A(A[0]);
    falcON_DEL_A(A);
  }
  //----------------------------------------------------------------------------
  template <typename T>
  inline void Alloc3D(T*** &A, int N[3]) falcON_THROWING {
    register int i, iN1;
    A       = falcON_NEW(T**,N[0]);
    A[0]    = falcON_NEW(T* ,N[0]*N[1]);
    A[0][0] = falcON_NEW(T  ,N[0]*N[1]*N[2]);
    const int N12 = N[1]*N[2];
    for(int i=0,iN1=0,iN12=0; i!=N[0]; ++i, iN1+=N[1], iN12+=N12) {
      A[i]    = A[0]    + iN1;
      A[i][0] = A[0][0] + iN12;
      for(int j=0, jN2=0; j!=N[1]; ++j, jN2+=N[2])
	A[i][j] = A[i][0] + jN2;
    }
  }
  //----------------------------------------------------------------------------
  template <typename T>
  inline void Alloc3D(T*** &A, int N[3], T const&X) falcON_THROWING {
    register int i, iN1;
    A       = falcON_NEW(T**,N[0]);
    A[0]    = falcON_NEW(T* ,N[0]*N[1]);
    A[0][0] = falcON_NEW(T  ,N[0]*N[1]*N[2]);
    const int N12 = N[1]*N[2];
    for(int i=0,iN1=0,iN12=0; i!=N[0]; ++i, iN1+=N[1], iN12+=N12) {
      A[i]    = A[0]    + iN1;
      A[i][0] = A[0][0] + iN12;
      for(int j=0, jN2=0; j!=N[1]; ++j, jN2+=N[2]) {
	A[i][j] = A[i][0] + jN2;
	for(int k=0; k!=N[2]; ++k) A[i][j][k] = X;
      }
    }
  }
  //----------------------------------------------------------------------------
  template <typename T>
  inline void Free3D(T*** A) falcON_THROWING {
    falcON_DEL_A(A[0][0]);
    falcON_DEL_A(A[0]);
    falcON_DEL_A(A);
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // Arrays<T,D> of type T and arbitrary dimension D                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#if defined(DEBUG) || defined(EBUG)
#  define THROW_BAD falcON_THROWING
#  define CHECK_BAD(M,R)					\
  if(i< 0) {							\
    warning("i=%d <  0 in Array subscript #%d",i,R);		\
    falcON_THROW("i=%d <  0 in Array subscript #%d",i,R);	\
  }								\
  if(i>=M) {							\
    warning("i=%d >= N=%d in Array subscript #%d",i,M,R);	\
    falcON_THROW("i=%d >= N=%d in Array subscript #%d",i,M,R);	\
  }
#else
#  define THROW_BAD
#  define CHECK_BAD(M,R)
#endif
  //////////////////////////////////////////////////////////////////////////////
  template<typename, int> class ConstPseudoArray;
  template<typename, int> class PseudoArray;
  template<typename, int> class Array;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // ConstPseudoArray<T,D> acts like a pointer of type const P*...* with D *s //
  // only public member is the [] operator which returns a                    //
  // 'ConstPseudoArray<T,D-1>', or if D==1, a 'T const&'                      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename T, int D> class ConstPseudoArray {
    friend class ConstPseudoArray<T,D+1>;
    friend class PseudoArray<T,D+1>;
    friend class Array<T,D+1>;
    const int*const N;
    const int*const K;
    const T  *const A;
    ConstPseudoArray(const T*a, const int*n, const int*k) : N(n), K(k), A(a) {}
    //--------------------------------------------------------------------------
  public:
    int const& size() const { return N[0]; }
    ConstPseudoArray<T,D-1> operator[] (int i) const THROW_BAD {
      CHECK_BAD(N[0],D);
      return ConstPseudoArray<T,D-1>(A+i*K[0],N+1,K+1);
    }
  };
  //----------------------------------------------------------------------------
  template<typename T> class ConstPseudoArray<T,1> {
    friend class ConstPseudoArray<T,2>;
    friend class PseudoArray<T,2>;
    friend class Array<T,2>;
    const int*const N;
    const T  *const A;
    ConstPseudoArray(const T*a, const int*n, const int*) : N(n), A(a) {}
    //--------------------------------------------------------------------------
  public:
    T const & operator[](int i) const THROW_BAD {
      CHECK_BAD(N[0],1);
      return A[i];
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // PseudoArray<T,D> acts like a pointer of type P*...* with D *s            //
  // only public member are the [] operators which return a                   //
  // 'ConstPseudoArray<T,D-1>', or if D==1, a 'T const&'  and                 //
  // 'PseudoArray     <T,D-1>', or if D==1, a 'T      &'                      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename T, int D> class PseudoArray {
    friend class PseudoArray<T,D+1>;
    friend class Array<T,D+1>;
    const int*const N;
    const int*const K;
    T        *const A;
    PseudoArray(T*a, const int*n, const int*k) : N(n), K(k), A(a) {}
    //--------------------------------------------------------------------------
  public:
    PseudoArray<T,D-1> operator[] (int i) THROW_BAD {
      CHECK_BAD(N[0],D);
      return PseudoArray<T,D-1>(A+i*K[0],N+1,K+1);
    }
    ConstPseudoArray<T,D-1> operator[] (int i) const THROW_BAD {
      CHECK_BAD(N[0],D);
      return ConstPseudoArray<T,D-1>(A+i*K[0],K+1);
    }
  };
  //----------------------------------------------------------------------------
  template<typename T> class PseudoArray<T,1> {
    friend class PseudoArray<T,2>;
    friend class Array<T,2>;
    const int*const N;
    T        *const A;
    PseudoArray(T*a, const int*n, const int*) : N(n), A(a) {}
    //--------------------------------------------------------------------------
  public:
    T       & operator[](int i)       THROW_BAD {
      CHECK_BAD(N[0],1);
      return A[i];
    }
    T const & operator[](int i) const THROW_BAD {
      CHECK_BAD(N[0],1);
      return A[i];
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class Array<T,D> acts like an D-dimensional array T[N_1]...[N_D]         //
  // operator[] act as expected                                               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename T, int D=1> class Array {
    int N[D];                                      // N[d]: size in dimension d 
    int K[D];                                      // K[d] = Prod_i>d N[i]      
    T  *A;                                         // pointer: allocated memory 
    void set(const int*n) {                        // set N[d], K[d]            
      for(int d=0; d!=D; ++d)
	N[d] = n? n[d] : 0;
      K[D-1] = 1;
      for(int d=D-1; d!=0; --d)
	K[d-1] = K[d] * N[d];
    }
    bool equal(const int n[D]) const {
      for(int d=0; d!=D; ++d)
	if(N[d] != n[d]) return false;
      return true;
    }
  public:
    static const int rank = D;
    int const&size(int d) const { return N[d]; }
    //--------------------------------------------------------------------------
    void reset(const int n[D]) falcON_THROWING {
      if(A==0 || !equal(n) ) {
	if(A) falcON_DEL_A(A);
	set(n);
	A = falcON_NEW(T,K[0]*N[0]);
      }
    }
    void reset(const int n[D], T const&x) falcON_THROWING {
      reset(n);
      for(int i=0; i!=K[0]*N[0]; ++i) A[i] = x;
    }
    //--------------------------------------------------------------------------
    Array() : A(0) {
      set(0);
    }
    explicit
    Array(const int n[D]) falcON_THROWING : A(0) {
      reset(n);
    }
    Array(const int n[D], T const&x) falcON_THROWING : A(0) {
      reset(n,x);
    }
    //--------------------------------------------------------------------------
    ~Array() falcON_THROWING {
      if(A) {
	falcON_DEL_A(A);
	A = 0;
      }
      set(0);
    }
    //--------------------------------------------------------------------------
    PseudoArray<T,D-1> operator[] (int i) THROW_BAD {
      CHECK_BAD(N[0],D);
      return PseudoArray<T,D-1>(A+i*K[0],N+1,K+1);
    }
    //--------------------------------------------------------------------------
    ConstPseudoArray<T,D-1> operator[] (int i) const THROW_BAD {
      CHECK_BAD(N[0],D);
      return ConstPseudoArray<T,D-1>(A+i*K[0],N+1,K+1);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class Array<T,1> acts like an 1-dimensional array T[N]                   //
  // operator[] act as expected                                               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename T> class Array<T,1> {
    int N;
    T  *A;
  public:
    static const int rank = 1;
    int const&size() const { return N; }
    int const&size(int) const { return N; }
    //--------------------------------------------------------------------------
    void reset(int n) falcON_THROWING {
      if(A==0 || n != N) {
	if(A) falcON_DEL_A(A);
	N = n;
	A = falcON_NEW(T,N);
      }
    }
    void reset(int n, T const&x) falcON_THROWING {
      reset(n);
      for(int i=0; i!=N; ++i) A[i] = x;
    }
    void reset(const int n[1]) falcON_THROWING { reset(n[0]); }
    void reset(const int n[1], T const&x) falcON_THROWING { reset(n[0],x); }
    //--------------------------------------------------------------------------
    Array() : A(0), N(0) {}
    explicit
    Array(int n) falcON_THROWING : A(0) { reset(n); }
    explicit
    Array(const int n[1]) falcON_THROWING : A(0) { reset(n); }
    Array(const int n[1], T const&x) falcON_THROWING : A(0) { reset(n,x); }
    //--------------------------------------------------------------------------
    ~Array() falcON_THROWING { 
      if(A) {
	falcON_DEL_A(A);
	A = 0;
      }
      N = 0;
    }
    //--------------------------------------------------------------------------
    T      & operator[] (int i)       THROW_BAD {
      CHECK_BAD(N,1);
      return A[i];
    }
    //--------------------------------------------------------------------------
    T const& operator[] (int i) const THROW_BAD {
      CHECK_BAD(N,1);
      return A[i];
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class Array<T,0> is just a scalar                                        //
  // this class is provided for completeness only                             //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename T> class Array<T,0> {
    T A;
  public:
    static const int rank = 0;
    int const&size(int) const { return 1; }
    explicit Array(int*) {}
    Array(int*, T const&x) : A(x) {}
    //--------------------------------------------------------------------------
    T      & operator[] (int i)       THROW_BAD {
      CHECK_BAD(0,0);
      return A;
    }
    T const& operator[] (int i) const THROW_BAD {
      CHECK_BAD(0,0);
      return A;
    }
  };
#undef THROW_BAD
#undef CHECK_BAD
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_memory_h

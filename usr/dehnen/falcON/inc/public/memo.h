// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// memo.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
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
//-----------------------------------------------------------------------------+
#ifndef falcON_included_memo_h
#define falcON_included_memo_h

#ifndef falcON_included_cstddef
#  include <cstddef>
#  define falcON_included_cstddef
#endif

#ifndef falcON_included_exit_h
#  include <public/exith>
#endif
#ifndef falcON_included_inln_h
#  include <public/inln.h>
#endif

namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  // class nbdy::block_alloc<T>                                                 
  //                                                                            
  // - allocates elements of type T in blocks                                   
  // - gives elements out either singly or in small amounts                     
  // - provides a forward iterator for sequential access to all elements given  
  //   out sofar                                                                
  //                                                                            
  // In order to (i) guarantee the correct amount of elements (not known a      
  // priori) will be allocated and (ii) not to waste memory by allocating too   
  // many, we use the following strategy, similar to that of std::string.       
  //                                                                            
  // Elements are allocated in blocks, which in turn are organized in a linked  
  // list. When a number of new elements is to be allocated and the last block  
  // in the list cannot provide them, we allocated a new block and add it to    
  // the list. The number of elements allocated in this new block is taken to   
  // be a function of the total number of elements used sofar. This function    
  // must be provided by the user, otherwise the same number as in the last     
  // block is used.                                                             
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  template<typename T>
  class block_alloc {
    ////////////////////////////////////////////////////////////////////////////
    // 1. member types of class nbdy::block_alloc                               
    ////////////////////////////////////////////////////////////////////////////
    // 1.0 some public typedefs                                                 
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
    // 1.1 class nbdy::block_alloc::block                                       
    //     - allocates a large block of elements                                
    //     - gives elements out either singly or in small blocks                
    //     - provides pointers to first and last element given out sofar        
    ////////////////////////////////////////////////////////////////////////////
  private:
    class block {
    private:
      block      *NEXT;                            // next block in linked list 
      value_type *FIRST;                           // front element             
      value_type *END;                             // and of active elements    
      value_type *ENDTOT;                          // end of all elements       
      block();                                     // not implemented           
      //------------------------------------------------------------------------
    public:
      block(size_type const&n) :                   // constructor               
	NEXT    ( 0 ),                             //   no next block           
	FIRST   ( falcON_New(value_type,n) ),      //   first of n elements     
	END     ( FIRST ),                         //   no elements used yet    
	ENDTOT  ( FIRST + n ) {}                   //   end of all elements     
      //------------------------------------------------------------------------
      ~block() { delete[] FIRST; }                 // destructor: de-allocate   
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
    // 1.2 class nbdy::block_alloc::iterator                                    
    //     for forward sequential iteration through all boxes                   
    ////////////////////////////////////////////////////////////////////////////
    class iterator{
    private:
      block  *B;                                   // current block             
      pointer E;                                   // current element           
    public:
      iterator(block*b, pointer e)   : B(b), E(e) {}
      iterator(const iterator&I) : B(I.B), E(I.E) {}
      iterator& operator= (const iterator&I) { B=I.B; E=I.E; return *this; }
      //------------------------------------------------------------------------
      iterator& operator++() {                     // prefix ++                 
	++E;                                       //   increment pointer       
	if(E == B->end()) {                        //   IF(out of block) >      
	  B = B->next();                           //     get next block in list
	  E = B? B->front() : 0;                   //     get its first element 
	}                                          //   <                       
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
    // 2 member data of class nbdy::block_alloc                                 
    //   hold & manage the blocks, keep track of elements used and allocated    
    ////////////////////////////////////////////////////////////////////////////
  private:
    block       *FIRST;                            // first block in linked list
    block       *LAST;                             // last block in linked list 
    size_type    NTOT;                             // # elements allocated      
    size_type    NUSED;                            // # elements used           
    ////////////////////////////////////////////////////////////////////////////
    // 3 member functions of class nbdy::block_alloc                            
    ////////////////////////////////////////////////////////////////////////////
  public:
    block_alloc(                                   // constructor:              
		size_type Ns) :                    // I: # elements in 1st block
      FIRST ( falcON_Memory(new block(Ns)) ),      //   allocate first block    
      LAST  ( FIRST ),                             //   last=first block        
      NTOT  ( Ns ),                                //   # elements allocated    
      NUSED ( 0 )   {}                             //   # elements used sofar   
    //--------------------------------------------------------------------------
    ~block_alloc() {                               // destructor:               
      register block *A=FIRST, *N;                 //   actual & next block     
      while(A) {                                   //   WHILE(actual is valid) >
	N = A->next();                             //     get next block        
	delete A;                                  //     delete actual block   
	A = N;                                     //     set actual = next     
      }                                            //   <                       
    }
    //--------------------------------------------------------------------------
    pointer new_element() {                        // give out: another element 
      if(LAST->is_full()) {                        //   IF(last block is full) >
	register size_type New = LAST->N_alloc();  //     # elements to allocate
	LAST->link(falcON_Memory(new block(New))); //     allocate new block &  
	LAST = LAST->next();                       //     add it to linked list 
      }                                            //   <                       
      NUSED++;                                     //   # elements used         
      return LAST->new_element();                  //   return new element      
    }
    //--------------------------------------------------------------------------
    template<class estimator>
    pointer new_element(                           // give out: another element 
			estimator F) {             // I: class returning an     
                                                   //    estimate for # elements
                                                   //    still needed given     
                                                   //    # elements used sofar  
      if(LAST->is_full()) {                        //   IF(last block is full) >
	register size_type New = F(NUSED);         //     # elements to allocate
	LAST->link(falcON_Memory(new block(New))); //     allocate new block &  
	LAST = LAST->next();                       //     add it to linked list 
      }                                            //   <                       
      NUSED++;                                     //   # elements used         
      return LAST->new_element();                  //   return new element      
    }
    //--------------------------------------------------------------------------
    pointer new_elements(                          // give out: more elements   
			 size_type Ne) {           // I: # elements wanted      
      if(!LAST->has_free(Ne)) {                    //   IF(last block is full) >
	register size_type
	  New=max(Ne,LAST->N_alloc());             //     # elements to allocate
	LAST->link(falcON_Memory(new block(New))); //     allocate new block &  
	LAST = LAST->next();                       //     add it to linked list 
      }                                            //   <                       
      NUSED+=Ne;                                   //   # elements used         
      return LAST->new_elements(Ne);               //   return new elements     
    }
    //--------------------------------------------------------------------------
    template<class estimator>
    pointer new_elements(                          // give out: more elements   
			 size_type Ne,             // I: # elements wanted      
			 estimator F) {            // I: class returning an     
                                                   //    estimate for # elements
                                                   //    still needed given     
                                                   //    # elements used sofar  
      if(!LAST->has_free(Ne)) {                    //   IF(last block is full) >
	register size_type New=max(Ne,F(NUSED));   //     # elements to allocate
	LAST->link(falcON_Memory(new block(New))); //     allocate new block &  
	LAST = LAST->next();                       //     add it to linked list 
      }                                            //   <                       
      NUSED+=Ne;                                   //   # elements used         
      return LAST->new_elements(Ne);               //   return new elements     
    }
    //--------------------------------------------------------------------------
    int number_of_element(                         // R: number of given element
			  const_pointer E) const { // I: element                
      if(E==0) return -99;                         //   IF null: return -99     
      register int num=0;                          //   reset counter num       
      for(register block*B=FIRST; B; B=B->next())  //   loop blocks         >   
	if(E >= B->front() && E < B->end())        //     IF in block         > 
	  return num + int(E-B->front());          //       return num + n(in b)
	else                                       //     < ELSE >              
	  num += B->N_used();                      //   <   add to counter    < 
      return -99;                                  //   not found: return -99   
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
  //                                                                            
  // class nbdy::pool                                                           
  //                                                                            
  // - allocates blocks of K bytes in chunks of N elements                      
  // - single elements can be allocated and freed                               
  // - chunks are organized as linked list                                      
  // - free elements are kept in a linked list                                  
  //                                                                            
  // The actual number of bytes allocated per element is                        
  //                                                                            
  //   K' = max( K, sizeof(*void) ).                                            
  //                                                                            
  // Thus, for K < sizeof(void*) this code becomes inefficient and should not   
  // be used.                                                                   
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class pool {
    ////////////////////////////////////////////////////////////////////////////
    // 1. member types of class nbdy::pool                                      
    ////////////////////////////////////////////////////////////////////////////
    // 1.0 some public typedefs                                                 
    ////////////////////////////////////////////////////////////////////////////
  public:
    typedef size_t    size_type;                   // type of number of elements
    typedef ptrdiff_t difference_type;             // type of pointer difference
    ////////////////////////////////////////////////////////////////////////////
    // 1.1 struct nbdy::pool::link                                              
    ////////////////////////////////////////////////////////////////////////////
  private:
    struct link { link *NEXT; };                     // pter to next link       
    ////////////////////////////////////////////////////////////////////////////
    // 1.2 struct nbdy::pool::chunk                                             
    ////////////////////////////////////////////////////////////////////////////
#define LINK(NAME) static_cast<link*>(static_cast<void*>(NAME))
   struct chunk {
      char   *DATA;                                  // pter to allocated memory
      chunk  *NEXT;                                  // pter to next chunk      
      //------------------------------------------------------------------------
      chunk(size_type N,                             // I: # elements           
	    size_type Kp) :                          // I: size of elements     
	DATA ( falcON_New(char,N*Kp) ),              // allocate memory         
	NEXT ( 0 )                                   // reset pter to next chunk
      {
	const    char *END=DATA+N*Kp;                // beyond last byte        
	register char *l,*n;                         // now we link the elements
	for(l=DATA, n=DATA+Kp; n!=END; l+=Kp, n+=Kp) // loop bits of Kp bytes   
	  LINK(l)->NEXT = LINK(n);                   //   link to next one      
	LINK(l)->NEXT = 0;                           // last one: null link     
      }
      //------------------------------------------------------------------------
      ~chunk() { delete DATA; }                      // de-allocate memory      
    };
    ////////////////////////////////////////////////////////////////////////////
    // 2 member data of class nbdy::pool                                        
    ////////////////////////////////////////////////////////////////////////////
  private:
    const size_type N;                               // # elements              
    const size_type Kp;                              // sizeof(element)         
    chunk          *CHUNKS;                          // pter to 1st chunk       
    link           *HEAD;                            // pter to 1st free element
  
    ////////////////////////////////////////////////////////////////////////////
    // 3 member functions of class nbdy::pool                                   
    ////////////////////////////////////////////////////////////////////////////
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
    void grow() {                                    // add another chunk       
      chunk *c = falcON_Memory(new chunk(N,Kp));     // allocate new chunk      
      c->NEXT  = CHUNKS;                             //   add to list of chunks 
      CHUNKS   = c;                                  //   make it first chunk   
      HEAD     = LINK(CHUNKS->DATA);                 //   add elements to list  
    }
    //--------------------------------------------------------------------------
    // constructor                                                              
    //--------------------------------------------------------------------------
  public:
    pool(size_type n,                                // I: desired # elements   
	 size_type k) :                              // I: desired sizeof(elem) 
      N      ( n<1? 1 : n ),                         //  actual # elements      
      Kp     ( sizeof(link) < k? k : sizeof(link) ), //  actual sizeof(element) 
      CHUNKS ( falcON_Memory(new chunk(N,Kp)) ),     //  allocate 1st chunk     
      HEAD   ( LINK(CHUNKS->DATA) ) {}               //  get list of elements   
#undef LINK
    //--------------------------------------------------------------------------
    // destructor                                                               
    //--------------------------------------------------------------------------
    ~pool() {
      register chunk *a=CHUNKS, *n;                  //  actual & next chunk    
      while(a) {                                     //  WHILE(actual is valid)>
	n = a->NEXT;                                 //    get next chunk       
	delete a;                                    //    delete actual chunk  
	a = n;                                       //    set actual = next    
      }                                              //  <                      
    }
    //--------------------------------------------------------------------------
    // hand out an element, i.e. a pter to Kp (=at least K) free bytes          
    //--------------------------------------------------------------------------
    void *alloc() {
      if(HEAD==0) grow();                            // no more elements: grow  
      link*p = HEAD;                                 // get 1st free element    
      HEAD   = p->NEXT;                              // reset pter to ----      
      return p;                                      // return free element     
    }
    //--------------------------------------------------------------------------
    // take back an element                                                     
    //--------------------------------------------------------------------------
    void free(void* e) {                             // I: pter to element      
      link *p = static_cast<link*>(e);               // cast into our format    
      p->NEXT = HEAD;                                // add to list of free     
      HEAD    = p;                                   // elements                
    }
  };
}
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_memo_h

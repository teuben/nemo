// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// iact.h                                                                      |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_iact_h
#define falcON_included_iact_h
#ifndef falcON_included_exit_h
#  include <public/exit.h>
#endif
#ifdef _OPENMP
#  ifndef falcON_included_omp_h
#    include <omp.h>
#    define falcON_included_omp_h
#  endif
#endif

#if defined(WRITE_IACT_INFO) && !defined(DEBUG)
#  define DEBUG
#endif

namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // auxiliary stuff                                                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  //----------------------------------------------------------------------------
  // struct iaction<>                                                           
  //----------------------------------------------------------------------------
  template<typename A, typename B>
  struct iaction {
#ifdef WRITE_IACT_INFO
    int lev;                                       // level of iaction          
#endif
    A fst;                                         // first object              
    B snd;                                         // second object             
    iaction () : fst(0), snd(0)                    // constructor               
#ifdef WRITE_IACT_INFO
    ,lev(0)
#endif
    {}
    void set(A a, B b                              // set pair                  
#ifdef WRITE_IACT_INFO
	     , int l
#endif
	     ) {
      fst=a;  snd=b;
#ifdef WRITE_IACT_INFO
      lev=l;
#endif
    }
  };
  //----------------------------------------------------------------------------
  // struct iastack<>                                                           
  //----------------------------------------------------------------------------
  template<typename A, typename B>
  class iastack {
  private:
    typedef iaction<A,B> iact;                     // type of stack objects     
    iact    *IA, *pi;                              // first & current element   
#ifdef DEBUG
    iact    *IEND;
#endif
    //--------------------------------------------------------------------------
  public:
    iastack (unsigned const&M)                     // constructor               
      : IA   ( falcON_New(iact,M) ),               //   allocate memory         
	pi   ( IA - 1 )                            //   set pter to activ       
#ifdef DEBUG
      , IEND ( IA + M )
#endif
    {
//       falcON_MemoryCheck(IA = new iact[M]);
//       pi = IA-1;
// #ifdef DEBUG
//       IEND = IA + M;
// #endif
    }
    //--------------------------------------------------------------------------
    ~iastack () { delete[] IA; }                   // destructor: deallocate    
    bool is_empty() const   { return pi<IA;}       // if stack empty?           
    iact pop     ()         { return *(pi--); }    // give last:   pop          
    void push    (A a, B b                         // add element: push         
#ifdef WRITE_IACT_INFO
		  ,int l=0
#endif
		  ) {
      ++pi;
#ifdef DEBUG
      if(pi >= IEND) error("push()ing beyond end of iastack");
#endif
      pi->set(a,b
#ifdef WRITE_IACT_INFO
	      ,l
#endif
	      );
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::InteractionBase<>                                            //
  //                                                                          //
  // encodes the mutual interaction algorithm.                                //
  // the order of nodes (first,second) is preserved when splitting on node    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename INTERACTOR> class InteractionBase {
    //--------------------------------------------------------------------------
    const INTERACTOR *IA;                            // interactor              
    //--------------------------------------------------------------------------
  protected:
    InteractionBase(const INTERACTOR* const&ia) : IA(ia) {}
    //--------------------------------------------------------------------------
#ifdef __GNUC__
  private:
#endif
    typedef typename INTERACTOR::cell_iter cell_iter;// iterator over cells     
    typedef typename INTERACTOR::soul_iter soul_iter;// iterator over souls     
    //--------------------------------------------------------------------------
    typedef iaction<cell_iter,cell_iter> cx_iact;    // cell-self iaction       
    typedef iaction<cell_iter,cell_iter> cc_iact;    // cell-cell iaction       
    typedef iaction<cell_iter,soul_iter> cs_iact;    // cell-soul iaction       
    typedef iaction<soul_iter,cell_iter> sc_iact;    // soul-cell iaction       
    //--------------------------------------------------------------------------
    typedef iastack<cell_iter,cell_iter> cx_stack;   // stack: cell-self iaction
    typedef iastack<cell_iter,cell_iter> cc_stack;   // stack: cell-cell iaction
    typedef iastack<cell_iter,soul_iter> cs_stack;   // stack: cell-soul iaction
    typedef iastack<soul_iter,cell_iter> sc_stack;   // stack: soul-cell iaction
    //--------------------------------------------------------------------------
#define LoopCKids(C,A)    LoopCellKids(typename cell_iter,C,A)
#define LoopSKids(C,A)    LoopSoulKids(typename cell_iter,C,A)
#define LoopCPairs(C,A,B) LoopCellSecd(typename cell_iter,C,A,B)
#define LoopSPairs(C,A,B) LoopSKids(C,A) LoopSoulSecd(typename cell_iter,C,A+1,B)
    //--------------------------------------------------------------------------
  public:
    inline void soul_soul(soul_iter const&s1, soul_iter const&s2) const
    {
      IA->interact(s1,s2);
    }
    //--------------------------------------------------------------------------
  protected:
    void clear_soul_cell_stack(sc_stack&SC) const {
      // clear a stack of soul-cell interactions                                
      while(!SC.is_empty()) {                        // WHILE(SC non-empty)     
	register sc_iact  sc = SC.pop();             //   pop new S-C iaction   
#ifdef WRITE_IACT_INFO
	std::cerr<<" SC "<<setw(2)<<sc.lev<<" :"
		 <<std::setw(6)<< sc.snd.my_tree()->index(sc.fst)<<' '
		 <<std::setw(6)<< sc.snd.index();
#endif
	if(!IA->interact(sc.fst,sc.snd)) {           //   IF(not performed)     
#ifdef WRITE_IACT_INFO
	  std::cerr<<": splitting\n";
#endif
	  LoopCKids(sc.snd,c2) SC.push  (sc.fst,c2   //     push    sub S-C     
#ifdef WRITE_IACT_INFO
					 ,sc.lev+1
#endif
					 );
	  LoopSKids(sc.snd,s2) soul_soul(sc.fst,s2); //     perform sub S-S     
	}                                            //   ENDIF                 
#ifdef WRITE_IACT_INFO
	else std::cerr<<": DONE\n";
#endif
      }                                              // END WHILE               
    }
    //--------------------------------------------------------------------------
    void clear_cell_soul_stack(cs_stack&CS) const {
      // clear a stack of cell-soul interactions                                
      while(!CS.is_empty()) {                        // WHILE(CS non-empty)     
	register cs_iact cs = CS.pop();              //   pop new C-S iaction   
#ifdef WRITE_IACT_INFO
	std::cerr<<" CS "<<setw(2)<<cs.lev<<" :"
		 <<std::setw(6)<< cs.fst.index()<<' '
		 <<std::setw(6)<< cs.fst.my_tree()->index(cs.snd);
#endif
	if(!IA->interact(cs.fst,cs.snd)) {           //   IF(not performed)     
#ifdef WRITE_IACT_INFO
	  std::cerr<<": splitting\n";
#endif
	  LoopCKids(cs.fst,c1) CS.push  (c1,cs.snd   //     push    sub C-S     
#ifdef WRITE_IACT_INFO
					 ,cs.lev+1
#endif
					 );
	  LoopSKids(cs.fst,s1) soul_soul(s1,cs.snd); //     perform sub S-S     
	}                                            //   ENDIF                 
#ifdef WRITE_IACT_INFO
	else std::cerr<<": DONE\n";
#endif
      }                                              // END WHILE               
    }
    //--------------------------------------------------------------------------
    void split_cell_cell(cc_iact &cc,
			 cc_stack&CC,
			 cs_stack&CS,
			 sc_stack&SC) const {
      // split cell-cell interaction and perform resulting soul interactions    
      if(IA->split_first(cc.fst,cc.snd)) {           // IF(split 1st)           
	LoopCKids(cc.fst,c1) CC.push(c1,cc.snd       //   push sub C-C          
#ifdef WRITE_IACT_INFO
				     ,cc.lev+1
#endif
				     );
	if(has_soul_kids(cc.fst)) {                  //   IF(soul kids)         
	  LoopSKids(cc.fst,s1) SC.push(s1,cc.snd     //     push sub S-C        
#ifdef WRITE_IACT_INFO
				       ,cc.lev+1
#endif
				     );
	  clear_soul_cell_stack(SC);                 //     clear the SC stack  
	}                                            //   ENDIF                 
      } else {                                       // ELSE(split 2nd)         
	LoopCKids(cc.snd,c2) CC.push(cc.fst,c2       //   push sub C-C          
#ifdef WRITE_IACT_INFO
				     ,cc.lev+1
#endif
				     );
	if(has_soul_kids(cc.snd)) {                  //   IF(soul kids)         
	  LoopSKids(cc.snd,s2) CS.push(cc.fst,s2     //     push sub C-S        
#ifdef WRITE_IACT_INFO
				       ,cc.lev+1
#endif
				       );
	  clear_cell_soul_stack(CS);                 //     clear the CS stack  
	}                                            //   ENDIF                 
      }                                              // ENDIF                   
    }
    //--------------------------------------------------------------------------
    void clear_cell_cell_stack(cc_stack&CC,
			       cs_stack&CS,
			       sc_stack&SC) const {
      // clear an stack of cell-cell interactions                               
      while(!CC.is_empty()) {                        // WHILE(CC non-empty)     
	register cc_iact cc = CC.pop();              //   pop new C-C iaction   
#ifdef WRITE_IACT_INFO
	std::cerr<<" CC "<<setw(2)<<cc.lev<<" :"
		 <<std::setw(6)<< cc.fst.index()<<' '
		 <<std::setw(6)<< cc.snd.index();
#endif
	if(!IA->interact(cc.fst,cc.snd)) {           //   IF(not performed)     
#ifdef WRITE_IACT_INFO
	  std::cerr<<": splitting\n";
#endif
	  split_cell_cell(cc,CC,CS,SC);              //     split               
	}
#ifdef WRITE_IACT_INFO
	else std::cerr<<": DONE\n";
#endif
      }                                              // END WHILE               
    }
    //--------------------------------------------------------------------------
  private:
    void work_cell_cell_stack(cc_stack &CC,
			      cs_stack &CS,
			      sc_stack &SC,
			      int      &i,
			      int const&M) const {
      // perform up to M cell-cell, perform all occuring soul interactions      
      while(!CC.is_empty() && i<M) {                 // WHILE(CC non-empty)     
	register cc_iact cc = CC.pop();              //   pop new C-C iaction   
	if(IA->interact(cc.fst,cc.snd)) ++i;         //   IF(performed): count  
	else split_cell_cell(cc,CC,CS,SC);           //   ELSE:          split  
      }                                              // END WHILE               
    }
    //--------------------------------------------------------------------------
  protected:
    void work_cell_cell_stack(cc_stack &CC,
			      cs_stack &CS,
			      sc_stack &SC,
			      int const&M) const {
      // perform up to M cell-cell, perform all occuring soul interactions      
      register int i=0;                              // counter: performed C-C  
      work_cell_cell_stack(CC,CS,SC,i,M);            // work it out             
    }
    //--------------------------------------------------------------------------
    void split_cell_self(cx_iact &cx,
			 cx_stack&CX,
			 cc_stack&CC,
			 cs_stack&CS) const {
      // split cell-self interaction and perform resulting soul interactions    
      LoopCKids(cx.fst,c1) {                         // LOOP(cell kids)         
	CX.push(c1,c1                                //   push sub C-X          
#ifdef WRITE_IACT_INFO
		,cx.lev+1
#endif
		);
	LoopSKids (cx.fst,s2)      CS.push(c1,s2   //   push sub C-S          
#ifdef WRITE_IACT_INFO
					   ,cx.lev+1
#endif
					   );
	LoopCPairs(cx.fst,c1+1,c2) CC.push(c1,c2   //   push sub C-C          
#ifdef WRITE_IACT_INFO
					   ,cx.lev+1
#endif
					   );
      }                                              // END LOOP                
      LoopSPairs(cx.fst,s1,s2) soul_soul(s1,s2);     // perform sub S-S         
      clear_cell_soul_stack(CS);                     // clear the CS stack      
    }
    //--------------------------------------------------------------------------
    void clear_cell_self_stack(cx_stack&CX,
			       cc_stack&CC,
			       cs_stack&CS,
			       sc_stack&SC) const {
      // clear an stack of cell-self interactions                               
      while(!CX.is_empty()) {                        // WHILE(CX non-empty)     
	register cx_iact cx = CX.pop();              //   pop new self iaction  
#ifdef WRITE_IACT_INFO
	std::cerr<<" CX "<<setw(2)<<cx.lev<<" :"
		 <<std::setw(6)<< cx.fst.index() <<"       ";
#endif
	if(!IA->interact(cx.fst)) {                  //   IF(not performed)     
#ifdef WRITE_IACT_INFO
	  std::cerr<<": splitting\n";
#endif
	  split_cell_self(cx,CX,CC,CS);              //     split               
	  clear_cell_cell_stack(CC,CS,SC);           //     clear the CC stack  
	}                                            //   ENDIF                 
#ifdef WRITE_IACT_INFO
	else std::cerr<<": DONE\n";
#endif
      }                                              // END WHILE               
    }
    //--------------------------------------------------------------------------
    void work_cell_self_stack(cx_stack &CX,
			      cc_stack &CC,
			      cs_stack &CS,
			      sc_stack &SC,
			      int const&M) const {
      // perform up to M cell iactions, perform all occuring soul iactions      
      register int i=0;                              // counter: performed      
      work_cell_cell_stack(CC,CS,SC,i,M);            // work on C-C stack       
      while(!CX.is_empty() && i<M) {                 // WHILE(CX non-empty)     
	register cx_iact cx = CX.pop();              //   pop new self iaction  
	if(IA->interact(cx.fst))                     //   IF(not performed)     
	  ++i;                                       //     count               
	else {                                       //   ELSE:                 
	  split_cell_self(cx,CX,CC,CS);              //     split               
	  work_cell_cell_stack(CC,CS,SC,i,M);        //     work on CC stack    
	}                                            //   ENDIF                 
      }                                              // END WHILE               
    }    
    //--------------------------------------------------------------------------
    static unsigned n_x(unsigned const&d) { return Nsub*d+1; }
    static unsigned n_c(unsigned const&d) { return Nsub*(Nsub-1)/2+Nsub*d; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::MutualInteractor<>                                           //
  //                                                                          //
  // encodes the mutual interaction algorithm.                                //
  // the order of nodes (first,second) is preserved when splitting a node     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename INTERACTOR> class MutualInteractor :
    public InteractionBase<INTERACTOR> {
#ifdef __GNUC__
    //--------------------------------------------------------------------------
    // typedefs identical to those in InteractorBase<>                          
    // We repeat them here, because otherwise gcc 3.2 complains, eg.:           
    //                                                                          
    // "`typename nbdy::MutualInteractor<INTERACTOR>::cell_iter' is implicitly  
    //  a typename"                                                             
    //                                                                          
    // (which looks like a compiler bug to me).                                 
    //--------------------------------------------------------------------------
    typedef typename INTERACTOR::cell_iter cell_iter;// iterator over cells     
    typedef typename INTERACTOR::soul_iter soul_iter;// iterator over souls     
    //--------------------------------------------------------------------------
    typedef iaction<cell_iter,cell_iter> cx_iact;    // cell-self iaction       
    typedef iaction<cell_iter,cell_iter> cc_iact;    // cell-cell iaction       
    typedef iaction<cell_iter,soul_iter> cs_iact;    // cell-soul iaction       
    typedef iaction<soul_iter,cell_iter> sc_iact;    // soul-cell iaction       
    //--------------------------------------------------------------------------
    typedef iastack<cell_iter,cell_iter> cx_stack;   // stack: cell-self iaction
    typedef iastack<cell_iter,cell_iter> cc_stack;   // stack: cell-cell iaction
    typedef iastack<cell_iter,soul_iter> cs_stack;   // stack: cell-soul iaction
    typedef iastack<soul_iter,cell_iter> sc_stack;   // stack: soul-cell iaction
#endif
    //--------------------------------------------------------------------------
    mutable cx_stack CX;                             // stack: cell-self iaction
    mutable cc_stack CC;                             // stack: cell-cell iaction
    mutable cs_stack CS;                             // stack: cell-soul iaction
    mutable sc_stack SC;                             // stack: soul-cell iaction
    //--------------------------------------------------------------------------
  public:
    MutualInteractor(                                // constructor             
		     INTERACTOR*const&ia,            // I: interactor           
		     unsigned   const&d1) :          // I: depth of 1st tree    
      InteractionBase<INTERACTOR> ( ia ),	     //   initialize base       
      CX ( n_x(d1) ),                                //   initialize cell-self  
      CC ( n_c(2*d1) ),                              //   initialize cell-cell  
      CS ( n_c(2*d1) ),                              //   initialize cell-soul  
      SC ( n_c(2*d1) ) {}                            //   initialize soul-cell  
    //--------------------------------------------------------------------------
    MutualInteractor(                                // constructor             
		     INTERACTOR*const&ia,            // I: interactor           
		     unsigned   const&d1,            // I: depth of 1st tree    
		     unsigned   const&d2) :          // I: depth of 2nd tree    
      InteractionBase<INTERACTOR> ( ia ),	     //   initialize base       
      CX ( n_x(d1) ),                                //   initialize cell-self  
      CC ( n_c(d2? d1+d2 : 2*d1) ),                  //   initialize cell-cell  
      CS ( n_c(d2? d1+d2 : 2*d1) ),                  //   initialize cell-soul  
      SC ( n_c(d2? d1+d2 : 2*d1) ) {}                //   initialize soul-cell  
    //--------------------------------------------------------------------------
    void cell_self(cell_iter const&a) const {        // cell-self interaction   
      CX.push(a,a);                                  // initialize stack        
      clear_cell_self_stack(CX,CC,CS,SC);            // work until stack empty  
    }
    //--------------------------------------------------------------------------
    void cell_cell(cell_iter const&a,
		   cell_iter const&b) const {        // cell-cell interaction   
      if(a==b)falcON_ErrorF("self-interaction","MutualInteractor::cell_cell()");
      CC.push(a,b);                                  // initialize stack        
      clear_cell_cell_stack(CC,CS,SC);               // work until stack empty  
    }
    //--------------------------------------------------------------------------
    void cell_soul(cell_iter const&a,
		   soul_iter const&b) const {        // cell-soul interaction   
      CS.push(a,b);                                  // initialize stack        
      clear_cell_soul_stack(CS);                     // work until stack empty  
    }
    //--------------------------------------------------------------------------
    void soul_cell(soul_iter const&a,
		   cell_iter const&b) const {        // soul-cell interaction   
      SC.push(a,b);                                  // initialize stack        
      clear_soul_cell_stack(SC);                     // work until stack empty  
    }
    //--------------------------------------------------------------------------
    // these routines are to be used by the MPI parallel code                   
    // they allow to interrupt the clearing of the CX and/or CC stack, e.g., for
    // testing for the sending or receipt of an MPI message.                    
    // use init_..()   to add a new root-root interaction onto a stack          
    // use work_..()   to perform at most a pre-defined number of cell iactions 
    //                 (soul-cell and soul-soul iaction will are not counted)   
    //                 returns true if still some work is left.                 
    // use finish_..() to clear the stack                                       
    //--------------------------------------------------------------------------
    void init_cell_self(cell_iter const&a) const {   // add cell-self iaction   
      CX.push(a,a);                                  // initialize stack        
    }
    //--------------------------------------------------------------------------
    bool work_cell_self(                             // work on CX, performing  
			int const&M) const {         // I: max # CX & CC iaction
      work_cell_self_stack(CX,CC,CS,SC,M);           // work at most M CC       
      return !CX.is_empty();                         // still work left         
    }
    //--------------------------------------------------------------------------
    void finish_cell_self() const {                  // clear CX stack          
      clear_cell_self_stack(CX,CC,CS,SC);            // work until stack empty  
    }
    //==========================================================================
    void init_cell_cell(cell_iter const&a,
			cell_iter const&b) const {   // add cell-cell iaction   
      CC.push(a,b);                                  // initialize stack        
    }
    //--------------------------------------------------------------------------
    bool work_cell_cell(                             // work on CC, performing  
			int const&M) const {         // I: max # CC iaction     
      work_cell_cell_stack(CC,CS,SC,M);              // work at most M CC       
      return !CC.is_empty();                         // still work left         
    }
    //--------------------------------------------------------------------------
    void finish_cell_cell() const {                  // clear CC stack          
      clear_cell_cell_stack(CC,CS,SC);               // work until stack empty  
    }
    //--------------------------------------------------------------------------
  };
#ifdef _OPENMP
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // auxiliary stuff                                                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  //----------------------------------------------------------------------------
  // struct union_iter<>                                                        
  //----------------------------------------------------------------------------
  template<typename A, typename B>
  struct union_iter {
    A     C;
    B     S;
    //--------------------------------------------------------------------------
    union_iter() : C(), S(0) {}
    union_iter(A const&c) : C(c), S(0) {}
    union_iter(B const&s) : C(),  S(s) {}
  };
#define is_cell(N) (N.S==0)
  //----------------------------------------------------------------------------
  // struct piaction<>                                                          
  //----------------------------------------------------------------------------
  template<typename A, typename B>
  struct piaction {
    union_iter<A,B> n1,n2;
  };
  //----------------------------------------------------------------------------
  // struct liaction<>                                                          
  //----------------------------------------------------------------------------
  template<typename A, typename B>
  struct liaction : public piaction<A,B> {
    omp_lock_t *l1, *l2;
    liaction   *next;
    //--------------------------------------------------------------------------
    liaction() : l1(0), l2(0), next(0) {}
    //--------------------------------------------------------------------------
    void try_to_get(bool &get) {
      get = false;
      if(omp_test_lock(l1)) {
	if(l1==l2 || omp_test_lock(l2)) get = true;
	else                            omp_unset_lock(l1);
      }
    }
    //--------------------------------------------------------------------------
    void free_locks() {
      omp_unset_lock(l1);
      if(l1 != l2) omp_unset_lock(l2);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::SelfInteractorP<>                                            //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename IACTOR>
  class SelfInteractorP : public InteractionBase<IACTOR> {
    //--------------------------------------------------------------------------
    typedef union_iter<cell_iter,soul_iter> node_iter;
    typedef liaction <cell_iter,soul_iter> l_iaction;
    //--------------------------------------------------------------------------
    unsigned                       nX,nC;          // size of stacks            
    int                            Nnod,Ncut;
    l_iaction                     *iact,*head;
    omp_lock_t                    *lock;
    //--------------------------------------------------------------------------
    inline void add_iact_to_list(register l_iaction* &iac,
                                 register int i1, register int i2,
				 node_iter* const&node)
    {
      if(i1 < Nnod && i2 < Nnod) {
        iac->n1   = node[i1];
        iac->n2   = node[i2];
        iac->l1   = lock+i1;
        iac->l2   = lock+i2;
        iac->next = head;
        head      = iac++;
      }
    }
    //--------------------------------------------------------------------------
    void get_iaction(l_iaction* &i) {
      register l_iaction **p;                      // current iaction           
      register bool got=false;                     // get iaction?              
      do {                                         // UNTIL got an interaction: 
	i =  head;                                 //   current iaction = first 
	p = &head;                                 //   pointer to current      
	while(i && !got) {                         //   WHILE list non-empty   >
	  i->try_to_get(got);                      //     try to get iaction    
	  if(got) {                                //     IF got it             
#pragma omp critical(pialist_self)                 //       critical:           
	    *p = i->next;                          //       remove from list    
	  } else {                                 //     ELSE                  
	    p = &(i->next);                        //       update pointer: next
	    i =   i->next;                         //       try next            
	  }                                        //     ENDIF                 
	}                                          //   <                       
      } while(!got);                               // <                         
    }
    //--------------------------------------------------------------------------
    void split(cell_iter const&C, node_iter* const&node, int const&Nmax) {
      LoopSKids(C,s) {
	if(Nnod > Nmax)
	  falcON_ErrorF("exceeding max # nodes","SelfInteractorP::split()");
	node[Nnod++] = node_iter(s);
      }
      LoopCKids(C,c)
	if(number(c) > Ncut)
	  split(c,node,Nnod,Nmax);
	else {
	  if(Nnod > Nmax)
	    falcON_ErrorF("exceeding max # nodes","SelfInteractorP::split()");
	  node[Nnod++] = node_iter(c);
	}
    }
    //--------------------------------------------------------------------------
  public:
    SelfInteractorP(                               // constructor               
		    IACTOR*   const&ia,            // I: interactor             
		    cell_iter const&root,          // I: root of tree           
		    unsigned  const&d) :           // I: depth of tree          
      InteractionBase<IACTOR> ( ia ),		   //   initialize base         
      nX   ( n_x(d) ),                             //   initialize CX size      
      nC   ( n_c(d+d) ),                           //   initialize CC size      
      Nnod ( 0 ),
      head ( 0 )
    {
      const int f    = 1;
      const int pmax = omp_get_max_threads();
      const int fmax = f*pmax;
      const int Nmax = Nsub*Nsub*fmax;
      Ncut = number(root) / fmax;
//       node_iter *node = new node_iter[Nmax];    falcON_MemoryCheck(node);
      node_iter *node = falcON_New(node_iter,Nmax);
      split(root,node,Nmax);
//       iact = new l_iaction[(Nnod*(Nnod+1))/2];  falcON_MemoryCheck(iact);
      iact = falcON_New(l_iaction,(Nnod*(Nnod+1))/2);
//       lock = new omp_lock_t[Nnod];              falcON_MemoryCheck(lock);
      lock = falcON_New(omp_lock_t,Nnod);
      for(register int i=0; i!=Nnod; ++i) omp_init_lock(lock+i);
      register l_iaction *iac = iact;
      register int N = Nnod + (Nnod%2? 0 : 1);
      for(register int r=0; r!=N/2+1; ++r) {
        for(register int i1=0,i2=i1+r; i1!=N; ++i1, ++i2)
          add_iact_to_list(iac,i1,i2%N,node);
      }
      delete[] node;
    }
    //--------------------------------------------------------------------------
    ~SelfInteractorP() {
      delete[] iact;
      for(register int i=0; i!=Nnod; ++i)
        omp_destroy_lock(lock+i);
      delete[] lock;
    }
    //--------------------------------------------------------------------------
    void interact_parallel() {
#pragma omp parallel
      {
        iastack<cell_iter,cell_iter>   CX(nX);       // create CX stack         
        iastack<cell_iter,cell_iter>   CC(nC);       // create CC stack         
        iastack<cell_iter,soul_iter>   CS(nC);       // create CS stack         
        iastack<soul_iter,cell_iter>   SC(nC);       // create SC stack         
        register l_iaction *iac;                     // current interaction     
        while(head) {                                // while iaction on list:  
          get_iaction(iac);                          // get next iaction        
          if(is_cell(iac->n1))                       // 1     node 1 is cell    
            if(is_cell(iac->n2))                     // 1.1   node 2 is cell    
              if(iac->n1.C == iac->n2.C) {           // 1.1.1 self-interaction  
                CX.push(iac->n1.C,iac->n2.C);        //   initialize CX stack   
                clear_cell_self_stack(CX,CC,CS,SC);  //   work stack down       
              } else {                               // 1.1.2 mutual interaction
                CC.push(iac->n1.C,iac->n2.C);        //   initialize CC stack   
                clear_cell_cell_stack(CC,CS,SC);     //   work stack down       
              }                                      //                         
            else {                                   // 1.2   node 2 is soul    
              CS.push(iac->n1.C,iac->n2.S);          //   initialize CS stack   
              clear_cell_soul_stack(CS);             //   work stack down       
            }                                        //                         
          else                                       // 2     node 1 is soul    
            if(is_cell(iac->n2)) {                   // 2.1   node 2 is cell    
              SC.push(iac->n1.S,iac->n2.C);          //   initialize SC stack   
              clear_soul_cell_stack(SC);             //   work stack down       
            } else if(iac->n1.S != iac->n2.S)        // 2.2   node 2 is soul    
              soul_soul(iac->n1.S,iac->n2.S);        //   work it out           
          iac->free_locks();                         // free locks for nodes    
        }
      }
    }
    //--------------------------------------------------------------------------
  };
#undef  is_cell
#endif  // _OPENMP
  //////////////////////////////////////////////////////////////////////////////
#undef LoopCKids
#undef LoopSKids
#undef LoopCPairs
#undef LoopSPairs
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::basic_iactor<TREE>                                             
  //                                                                            
  // a base class for a possible template argument of MutualInteractor<>        
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  template<typename TREE> class basic_iactor {
    //--------------------------------------------------------------------------
    // types of class basic_iactor                                              
    //--------------------------------------------------------------------------
  public:
    typedef TREE                              tree_type;   // type of tree      
    typedef typename tree_type::cell_iterator cell_iter;   // iterator over cell
    typedef typename tree_type::soul_iterator soul_iter;   // iterator over soul
    //--------------------------------------------------------------------------
    // data of class basic_iactor                                               
    //--------------------------------------------------------------------------
  private:
    const int NCB, NCC, NCS;                         // params control direct   
    //--------------------------------------------------------------------------
    // abstract methods, MUST be provided by derived class                      
    //--------------------------------------------------------------------------
  protected:
    virtual void many       (bool const&, soul_iter const&,
			     soul_iter const&, soul_iter const&) const = 0;
    virtual void single     (soul_iter const&, soul_iter const&) const = 0;
    virtual bool discard    (cell_iter const&, soul_iter const&) const = 0;
    virtual bool discard    (cell_iter const&, cell_iter const&) const = 0;
  public:
    virtual bool split_first(cell_iter const&, cell_iter const&) const = 0;
    //--------------------------------------------------------------------------
    // other methods                                                            
    //--------------------------------------------------------------------------
  protected:
    basic_iactor(const int dir[4] = Default::direct) :
      NCB(dir[1]), NCC(dir[2]), NCS(dir[3]) {}
    //--------------------------------------------------------------------------
  private:
    bool many  (cell_iter const&A, soul_iter const&B) const {
      many(al_active(A) || is_active(B), B, A.begin_souls(), A.end_soul_desc());
      return true;
    }
    //--------------------------------------------------------------------------
    bool many  (cell_iter const&A, cell_iter const&B) const {
      LoopAllSouls(typename cell_iter,B,Bi) many(A,Bi);
      return true;
    }
    //--------------------------------------------------------------------------
    bool many  (cell_iter const&A) const {
      // debugged 12-may-2003                                                   
      // bug detected by Clayton Heller in falcON::make_iaction_list()          
      LoopLstSouls(typename cell_iter,A,Ai)
	many(is_active(Ai), Ai,Ai+1,A.end_soul_desc());
      return true;
    }
    //--------------------------------------------------------------------------
  public:
    bool interact(cell_iter const&A, soul_iter const&B) const {
      if(!is_active(A) && !is_active(B))return true;     // no actives -> no job
      if(discard(A,B))                  return true;     // try to discard      
      if(number(A) < NCB)               return many(A,B);// perform many single 
                                        return false;    // must split          
    }
    //--------------------------------------------------------------------------
    bool interact(soul_iter const&B, cell_iter const&A) const {
      if(!is_active(A) && !is_active(B))return true;     // no actives -> no job
      if(discard(A,B))                  return true;     // try to discard      
      if(number(A) < NCB)               return many(A,B);// perform many single 
                                        return false;    // must split          
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A, cell_iter const&B) const {
      if(!is_active(A) && !is_active(B))return true;     // no actives -> no job
      if(discard(A,B))                  return true;     // try to discard      
      if(number(A) < NCC &&
	 number(B) < NCC    )           return many(A,B);// perform many single 
                                        return false;    // must split          
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A) const {
      if(!is_active(A))              return true;        // no actives -> no job
      if(number(A) < NCS)            return many(A);     // perform many single 
                                     return false;       // must split          
    }
    //--------------------------------------------------------------------------
    void interact(soul_iter const&A, soul_iter const&B) const {
      single(A,B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::neighbour_counter                                            //
  //                                                                          //
  // for each active soul: count # souls with |R| < size(A);                  //
  // for individual sizes(soft_type=individual): requires:                    //
  //     - each active soul to hold its size and size^2.                      //
  //     - each active cell to hold its size                                  //
  // for global size, it is given as argument to the constructor              //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename TREE, soft_type SOFT>  class neighbour_counter;
  //////////////////////////////////////////////////////////////////////////////
  template<typename TREE>
  class neighbour_counter<TREE,individual> : public basic_iactor<TREE> {
    //--------------------------------------------------------------------------
    // nbdy::neighbour_counter<individual>                                      
    // MUST NOT be commented out via #ifdef falcON_INDI                         
    //--------------------------------------------------------------------------
  public:
    typedef typename TREE::soul_iterator soul_iter;
    typedef typename TREE::cell_iterator cell_iter;
  protected:
    //--------------------------------------------------------------------------
    void many(bool      const&all,
	      soul_iter const&A,
	      soul_iter const&B0,
	      soul_iter const&BN) const {
      register real      Rq;
      register soul_iter B;
      if(is_active(A))
	if(all)
	  for(B=B0; B<BN; B++) {
	    Rq = dist_sq(pos(A),pos(B));
	    if(Rq < sizeq(A)) A->inc();
	    if(Rq < sizeq(B)) B->inc();
	  }
	else
	  for(B=B0; B<BN; B++) {
	    Rq = dist_sq(pos(A),pos(B));
	    if(                Rq < sizeq(A)) A->inc();
	    if(is_active(B) && Rq < sizeq(B)) B->inc();
	  }
      else
	if(all)
	  for(B=B0; B<BN; B++) {
	    Rq = dist_sq(pos(A),pos(B));
	    if(Rq < sizeq(B)) B->inc();
	  }
	else
	  for(B=B0; B<BN; B++) if(is_active(B)) {
	    Rq = dist_sq(pos(A),pos(B));
	    if(Rq < sizeq(B)) B->inc();
	  }
    }
    //--------------------------------------------------------------------------
    void single(soul_iter const&A, soul_iter const&B) const {
      if(!is_active(A) && !is_active(B)) return;   // no actives -> no job      
      register real Rq = dist_sq(pos(A),pos(B));   // R^2 = distance^2          
      if(is_active(A) && Rq < sizeq(A)) A->inc();  // A is active && |R|<siza(A)
      if(is_active(B) && Rq < sizeq(B)) B->inc();  // B is active && |R|<siza(B)
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, soul_iter const&B) const
    {
      return dist_sq(pos(A),pos(B)) >              // |xA-xB| >                 
	square(max(rmax(A)+size(B),size(A)));      // max(sA,rA+sB)  ?          
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, cell_iter const&B) const
    {
      return dist_sq(pos(A),pos(B)) >              // |xA-xB| >                 
	square(max(rmax(A)+size(B),rmax(B)+size(A))); // max(rb+sA,rA+sB)  ?    
    }
  public:
    //--------------------------------------------------------------------------
    neighbour_counter() {}
    //--------------------------------------------------------------------------
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return nbdy::is_twig(B) || size(A) > size(B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  template<typename TREE>
  class neighbour_counter<TREE,global>
    : public basic_iactor<TREE> {
  public:
    typedef typename TREE::soul_iterator soul_iter;
    typedef typename TREE::cell_iterator cell_iter;
  private:
    const real EPS,EPQ;
  protected:
    //--------------------------------------------------------------------------
    void many(bool      const&all,
	      soul_iter const&A,
	      soul_iter const&B0,
	      soul_iter const&BN) const {
      register soul_iter B;
      if(is_active(A))
	if(all)
	  for(B=B0; B<BN; B++) {
	    if(dist_sq(pos(A),pos(B)) < EPQ) { A->inc();
	    B->inc(); }
	  }
	else
	  for(B=B0; B<BN; B++) {
	    if(dist_sq(pos(A),pos(B)) < EPQ) {                  A->inc();
	    if(is_active(B)) B->inc(); }
	  }
      else
	if(all)
	  for(B=B0; B<BN; B++) {
	    if(dist_sq(pos(A),pos(B)) < EPQ) { B->inc(); }
	  }
	else
	  for(B=B0; B<BN; B++) 
	    if(is_active(B) && dist_sq(pos(A),pos(B)) < EPQ) { B->inc(); }
    }
    //--------------------------------------------------------------------------
    void single(soul_iter const&A, soul_iter const&B) const {
      if(!is_active(A) && !is_active(B)) return;   // no actives -> no job      
      if(dist_sq(pos(A),pos(B)) < EPQ) {           // distance^2 < eps^2        
	if(is_active(A)) A->inc();                 // A is active && |R|<siza(A)
	if(is_active(B)) B->inc();                 // B is active && |R|<siza(B)
      }
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, soul_iter const&B) const
    {
      return dist_sq(pos(A),pos(B)) >                // |xA-xB| >               
	square(max(rmax(A)+EPS,size(A)));            // max(sA,rA+sB)  ?        
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, cell_iter const&B) const
    {
      return dist_sq(pos(A),pos(B)) >                // |xA-xB| >               
	square(max(rmax(A)+size(B),rmax(B)+size(A)));// max(rb+sA,rA+sB)  ?     
    }
  public:
    //--------------------------------------------------------------------------
    neighbour_counter(const real e) : EPS(e), EPQ(e*e) {}
    //--------------------------------------------------------------------------
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || size(A) > size(B); }
  };
  //////////////////////////////////////////////////////////////////////////////
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_iact_h    

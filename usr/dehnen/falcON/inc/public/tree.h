// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tree.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2002                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines                                                                     |
//                                                                             |
// class basic_soul                                                            |
// class basic_cell                                                            |
// class basic_tree<TREE,CELL>                                                 |
//                                                                             |
// macros for access to cells & souls of basic_tree<cell>                      |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef included_tree_h
#define included_tree_h

#ifndef included_auxx_h
#  include <public/auxx.h>
#endif
#ifndef included_flag_h
#  include <public/flag.h>
#endif
#ifndef included_iomanip
#  include <iomanip>
#  define included_iomanip
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  class sbodies;                                   // forward declaration       
#ifdef ALLOW_MPI
  class pbodies;                                   // forward declaration       
#endif
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::basic_soul                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class basic_soul : public flag {
  private:
    basic_soul           (const basic_soul&);      // not implemented           
    basic_soul& operator=(const basic_soul&);      // not implemented           
    //--------------------------------------------------------------------------
    // data of class basic_soul                                                 
    //--------------------------------------------------------------------------
    vect COFM;                                     // center of mass = position 
    uint LINK;                                     // associated body index     
    //--------------------------------------------------------------------------
    // constructor                                                              
    //--------------------------------------------------------------------------
  protected:
    basic_soul() {}                                // for optimization, we don't
                                                   // initialize COFM & LINK    
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
    void  link(unsigned  const&L) { LINK = L; }
    vect &cofm()                  { return COFM; }
    real &cofm(const int i)       { return COFM[i]; }
    //--------------------------------------------------------------------------
  public:
    void copy_link (const basic_soul*S) { LINK = S->LINK; }
    void copy_basic(const basic_soul*S) { COFM = S->COFM; LINK = S->LINK; }
    //--------------------------------------------------------------------------
    // data access                                                              
    //--------------------------------------------------------------------------
    vect const&cofm  () const { return COFM; }
    uint const&mybody() const { return LINK; }
    //--------------------------------------------------------------------------
    // flag manipulations to use with bodies                                    
    //--------------------------------------------------------------------------
    template<typename BODIES>
    void copy_body_flag(const BODIES*BB) {
      flag::set_to_part(&(BB->flg(mybody())), flag::BODY_FLAGS);
    }
    //--------------------------------------------------------------------------
    template<typename BODIES>
    void add_body_flag(const BODIES*BB) {
      flag::add_part(&(BB->flg(mybody())),flag::BODY_FLAGS);
    }
    //--------------------------------------------------------------------------
    template<typename BODIES>
    void set_body_flag(const BODIES*BB) {
      flag::set_part(&(BB->flg(mybody())),flag::BODY_FLAGS);
    }
    //--------------------------------------------------------------------------
    template<typename BODIES>
    void copy_sink_flag(const BODIES*BB) {
      flag::set_part(&(BB->flg(mybody())),flag::SINK);
    }
    //--------------------------------------------------------------------------
    // flag manipulations to use with arrays                                    
    //--------------------------------------------------------------------------
    void copy_body_flag(const int*F) {
      flag::set_to_part(F[mybody()],flag::BODY_FLAGS);
    }
    //--------------------------------------------------------------------------
    void add_body_flag(const int*F) {
      flag::add_part(F[mybody()],flag::BODY_FLAGS);
    }
    //--------------------------------------------------------------------------
    void copy_sink_flag(const int*F) {
      flag::set_part(F[mybody()],flag::SINK);
    }
    //--------------------------------------------------------------------------
    // other non-const methods                                                  
    //--------------------------------------------------------------------------
    void set_pos(const areal*x[NDIM]) {
      cofm()[0] = x[0][mybody()];
      cofm()[1] = x[1][mybody()];
#if NDIM==3
      cofm()[2] = x[2][mybody()];
#endif
    }
    //--------------------------------------------------------------------------
    template<typename BODIES>
    void set_pos(const BODIES*BB) {
      cofm() = BB->pos(mybody());
    }
    //--------------------------------------------------------------------------
    void set_link_and_pos(const uint&i, vect const&x) {
      link(i);
      cofm()=x;
    }
    //--------------------------------------------------------------------------
    void copy(const basic_soul*S)  {
      copy_basic(S);
      set_to_part(S,flag::BODY_FLAGS);
    }
    //--------------------------------------------------------------------------
    void copy_prune(const basic_soul*S, size_t const &L) { 
      COFM = S->COFM;
      LINK = L; 
      set_to_part(S,flag::BODY_FLAGS);
    }
    //--------------------------------------------------------------------------
    // dump soul basic data                                                     
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream&o) {
      o<<" flg mybody          position";
    }
    //--------------------------------------------------------------------------
    static void dump(std::ostream&o, const basic_soul* const&S) {
      o<<" "<<std::setw(3)<<flag(*S)
       <<" "<<std::setw(6)<<S->mybody();
      for(register int i=0; i!=NDIM; ++i)
	o<<" "<<std::setw(8)<<std::setprecision(4)<<S->cofm()[i];
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // nbdy::class basic_cell                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class basic_cell : public flag {
    //--------------------------------------------------------------------------
    // friendship                                                               
    //--------------------------------------------------------------------------
    friend class basic_cell_access;                // allows data access        
    //--------------------------------------------------------------------------
    // data of class basic_cell (write access through class basic_cell_access)  
    //--------------------------------------------------------------------------
  private:
    real  RADIUS;                                  // half edge length of cube  
    vect  CENTER;                                  // center of cube            
    int   NUMBER;                                  // # soul descendants        
    int   FCSOUL;                                  // index of fst soul child   
    int   FCCELL;                                  // index of fst cell child   
    short NSOULS;                                  // # soul children           
    short NCELLS;                                  // # cell children           
    //--------------------------------------------------------------------------
    // constructor                                                              
    //--------------------------------------------------------------------------
  protected:
    basic_cell() {}
    //--------------------------------------------------------------------------
    // data access                                                              
    //--------------------------------------------------------------------------
  public:
    real  const&radius() const { return RADIUS; }
    vect  const&center() const { return CENTER; }
    int   const&number() const { return NUMBER; }
    int   const&fcsoul() const { return FCSOUL; }
    int   const&fccell() const { return FCCELL; }
    short const&nsouls() const { return NSOULS; }
    short const&ncells() const { return NCELLS; }
    int   const ecsoul() const { return FCSOUL+NSOULS; }
    int   const ncsoul() const { return FCSOUL+NUMBER; }
    int   const eccell() const { return FCCELL+NCELLS; }
    //--------------------------------------------------------------------------
    // flag manipulations                                                       
    //--------------------------------------------------------------------------
    void reset_flag() {
      flag::reset();
      flag::add(flag::AL_SINK);
    }
    //--------------------------------------------------------------------------
    void reset_sink_flag() {
      flag::un_set(flag::SINK);
      flag::add(flag::AL_SINK);
    }
    //--------------------------------------------------------------------------
    void add_flag(int const&f) {
      flag::add(f);
    }
    //--------------------------------------------------------------------------
  protected:
    void add_sink_flag_from_soul(const flag* const&S) {
      flag::add_part(S,flag::SINK);
      if(!is_sink(S)) flag::un_set(flag::AL_SINK);
    }
    //--------------------------------------------------------------------------
    void add_sink_flag_from_cell(const flag* const&C) {
      flag::add_part(C,flag::SINK);
      if(!al_sink(C)) flag::un_set(flag::AL_SINK);
    }
    //--------------------------------------------------------------------------
    // boolean information                                                      
    //--------------------------------------------------------------------------
  public:
    bool has_cell_kids() const { return ncells()!=0; }
    bool has_soul_kids() const { return nsouls()!=0; }
    bool is_twig      () const { return ncells()==0; }
    bool is_branch    () const { return nsouls()==0; }
    //--------------------------------------------------------------------------
    // dump cell data                                                           
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream&o) {
      o<<" flg cells ncell souls nsoul number    radius            center";
    }
    //--------------------------------------------------------------------------
    static void dump(std::ostream&o, const basic_cell* const&C) {
      o<<" "<<std::setw(3)<<flag(*C);
      if(C->ncells())
	o<<" "<<std::setw(5)<<C->fccell();
      else
	o<<"     -";
      o<<" "<<std::setw(5)<<C->ncells()
       <<" "<<std::setw(5)<<C->fcsoul()
       <<" "<<std::setw(5)<<C->nsouls()
       <<" "<<std::setw(6)<<C->number()
       <<" "<<std::setw(9)<<C->radius();
      for(register int i=0; i<NDIM; ++i)
	o<<" "<<std::setw(8)<<std::setprecision(4)<<C->center()[i];
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  const int PRUNED      = 1<<10;
  const int PUREREP     = 1<<11;
  const int PRUNE_FLAGS = PRUNED | PUREREP;
  //////////////////////////////////////////////////////////////////////////////
  template<typename TREE>
  class tree_transport;                            // forward declaration     //
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::basic_tree<TREE,CELL>                                        //
  //                                                                          //
  // NOTE We cannot avoid the second template argument CELL. This is because  //
  //      in the usual application, TREE is derived from basic_tree<>, and    //
  //      cannot be used like a complete type within basic_tree. In particular//
  //      we cannot say                                                       //
  //        typedef typename TREE::cell_type cell_type;                       //
  //      This would create a compiler error, since type TREE is not complete //
  //      without basic_tree<TREE> and vice versa.                            //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename TREE, typename CELL>           // derived tree, cell of it  
  class basic_tree {
    //--------------------------------------------------------------------------
    friend class tree_transport<TREE>;             // for parallel code         
    //--------------------------------------------------------------------------
    basic_tree           (const basic_tree&);      // not implemented           
    basic_tree& operator=(const basic_tree&);      // not implemented           
    //--------------------------------------------------------------------------
    // types of class basic_tree                                                
    //--------------------------------------------------------------------------
  public:
    typedef CELL                         cell_type;// type of tree cell         
    typedef typename CELL::soul_type     soul_type;// type of tree soul         
    typedef basic_tree<TREE,CELL>        base_tree;// type of basic_tree        
    //--------------------------------------------------------------------------
    // type cell_iterator                                                       
    // Rationale: In order to access subcells in a parallel code, we cannot use 
    //            'pointer to cell' (or we re-adjust the pointers after each    
    //            MPI communication).                                           
    //            So, in class cell, not the pointers to the cell and soul      
    //            children are stored, but their indices. This means that in    
    //            order to access them, we do also need the tree that holds the 
    //            cell and soul arrays.                                         
    //            basic_tree::cell_iterator is protected to allow the derived   
    //            TREE to have its own cell_iterator, derived from this one.    
    //--------------------------------------------------------------------------
    class cell_iterator {
      friend class basic_tree;
      //........................................................................
      // types                                                                  
      //........................................................................
    public:
      typedef cell_iterator            cell_child; // type used in child access 
      typedef basic_tree::soul_type     soul_type; // type used in child access 
      typedef soul_type               *soul_child; // type used in child access 
      //........................................................................
      // data members                                                           
      //........................................................................
    private:
      const basic_tree*T;                          // pointer to basic_tree     
      CELL            *C;                          // pointer to derived cell   
      //........................................................................
      // construction                                                           
      //........................................................................
      cell_iterator(const basic_tree*const&t, CELL* const&c) : T(t),   C(c)   {}
    public:
      cell_iterator()                                        : T(0),   C(0)   {}
      cell_iterator(int)                                     : T(0),   C(0)   {}
      cell_iterator(cell_iterator const&I)                   : T(I.T), C(I.C) {}
      //........................................................................
      // assigning                                                              
      //........................................................................
      cell_iterator&operator= (cell_iterator const&I) {
	T = I.T;
	C = I.C;
	return *this;
      }
      //........................................................................
      // forward iteration                                                      
      //........................................................................
      cell_iterator&operator++()                 { ++C; return *this; }
      cell_iterator operator++(int)              { return cell_iterator(T,C++);}
      cell_iterator&operator+=(int const&k)      { C+=k; return *this; }
      cell_iterator operator+ (int const&k) const{ return cell_iterator(T,C+k);}
      //........................................................................
      // backward iteration                                                     
      //........................................................................
      cell_iterator&operator--()                 { --C; return *this; }
      cell_iterator operator--(int)              { return cell_iterator(T,C--);}
      cell_iterator&operator-=(int const&k)      { C-=k; return *this; }
      cell_iterator operator- (int const&k) const{ return cell_iterator(T,C-k);}
      //........................................................................
      // boolean methods                                                        
      //........................................................................
      bool operator==(cell_iterator const&I) const { return C == I.C; }
      bool operator!=(cell_iterator const&I) const { return C != I.C; }
      bool operator< (cell_iterator const&I) const { return C <  I.C; }
      bool operator<=(cell_iterator const&I) const { return C <= I.C; }
      bool operator> (cell_iterator const&I) const { return C >  I.C; }
      bool operator>=(cell_iterator const&I) const { return C >= I.C; }
      //........................................................................
      // tree and index                                                         
      //........................................................................
      const basic_tree*const&my_tree() const          { return T; }
      int              index  () const                { return C - T->root().C;}
      friend int       index  (cell_iterator const&I) { return I.index(); }
      //........................................................................
      // conversion to pointer to cell  and  dereferencing to cell*             
      //........................................................................
      operator CELL*const&()           const { return C; }
               CELL*const&operator->() const { return C; }
               CELL*const&c_pter()     const { return C; }
      //........................................................................
      // const access to basic_cell methods via friends                         
      //........................................................................
#define CICI cell_iterator const&I
      friend real  const&radius(CICI)        { return I->radius(); }
      friend vect  const&center(CICI)        { return I->center(); }
      friend int   const&number(CICI)        { return I->number(); }
      friend bool  const has_cell_kids(CICI) { return I->has_cell_kids();}
      friend bool  const has_soul_kids(CICI) { return I->has_soul_kids();}
      friend bool  const is_twig (CICI)      { return I->is_twig (); }
      friend bool  const is_branch(CICI)     { return I->is_branch(); }
      friend bool  const is_pruned(CICI)     { return I->is_set(PRUNED); }
      friend bool  const is_pure_rep(CICI)   { return I->is_set(PUREREP); }
#undef CICI
      //........................................................................
      // access to cell childreen                                               
      //........................................................................
      cell_child begin_cell_kids() const { 
	return cell_child(T, T->root() + C->fccell()); }
      cell_child end_cell_kids  () const {
	return cell_child(T, T->root() + C->eccell()); }
      //........................................................................
      // access to soul childreen and soul  descendants                         
      //........................................................................
      soul_child begin_souls   () const { return T->begin_souls()+ C->fcsoul();}
      soul_child end_soul_kids () const { return T->begin_souls()+ C->ecsoul();}
      soul_child end_soul_desc () const { return T->begin_souls()+ C->ncsoul();}
      soul_child last_soul_desc() const { return end_soul_desc()-1; }
      //........................................................................
      // output: give just the index (an iterator acts like a pointer)          
      //........................................................................
      friend std::ostream& operator<<(std::ostream&o, const cell_iterator&I) {
	return o<<std::setw(5)<<I.index();
      }
    };
    //--------------------------------------------------------------------------
    // data of class basic_tree                                                 
    //--------------------------------------------------------------------------
  private:
    const sbodies       *BSRCES;                   // pointer to sbodies        
#ifdef ALLOW_MPI
    const pbodies       *PSRCES;                   // pointer to pbodies        
#endif
    const int           *F;                        // array with flags          
    const areal        **X;                        // arrays with x,y,z         
    uint                 Nb,Ns,Nc;                 // # bodies/souls/cells      
    int                  Dp;                       // depth                     
    mutable soul_type   *S0, *SN;                  // begin & end of souls      
    mutable cell_type   *C0, *CN;                  // begin & end of cells      
    //--------------------------------------------------------------------------
    // construction and such                                                    
    //--------------------------------------------------------------------------
    basic_tree() :                                 // for friend tree_transport 
      BSRCES(0), 
#ifdef ALLOW_MPI
      PSRCES(0),
#endif
      F(0), X(0), Ns(0), Nc(0), S0(0), C0(0) {}
  public:
    //--------------------------------------------------------------------------
    // construct from sbodies                                                   
    //--------------------------------------------------------------------------
    basic_tree(                                    // openMP//ized              
	       const sbodies* const&,              // I: bodies' flags & pos's  
	       const int    = 6,                   //[I: N_crit]                
	       const int    = 50);                 //[I: max tree depth]        
    //--------------------------------------------------------------------------
    // construct from sbodies                                                   
    //--------------------------------------------------------------------------
    basic_tree(                                    // openMP//ized              
	       const sbodies* const&,              // I: bodies' flags & pos's  
	       vect           const&,              // I: x_min                  
	       vect           const&,              // I: x_max                  
	       const int    = 6,                   //[I: N_crit]                
	       const int    = 50);                 //[I: max tree depth]        
#ifdef ALLOW_MPI
    //--------------------------------------------------------------------------
    // construct from pbodies                                                   
    //--------------------------------------------------------------------------
    basic_tree(                                    // openMP//ized              
	       const pbodies* const&,              // I: bodies' flags & pos's  
	       const int    = 6,                   //[I: N_crit]                
	       const int    = 50);                 //[I: max tree depth]        
    //--------------------------------------------------------------------------
    // construct from pbodies                                                   
    //--------------------------------------------------------------------------
    basic_tree(                                    // openMP//ized              
	       const pbodies* const&,              // I: bodies' flags & pos's  
	       vect           const&,              // I: x_min                  
	       vect           const&,              // I: x_max                  
	       const int    = 6,                   //[I: N_crit]                
	       const int    = 50);                 //[I: max tree depth]        
#endif
    //--------------------------------------------------------------------------
    // construct from flags & pos                                               
    //--------------------------------------------------------------------------
    basic_tree(                                    // openMP//ized              
	       const int   *,                      // I: array with flags       
	       const areal *[NDIM],                // I: arrays with x,y,z      
	       const uint,                         // I: size of arrays         
	       const int   = 6,                    //[I: N_crit]                
	       const int   = 50);                  //[I: max tree depth]        
    //--------------------------------------------------------------------------
    // construct as sub-tree                                                    
    //--------------------------------------------------------------------------
    template<typename PARENT_TREE>                 //   type of parent-tree     
    basic_tree(                                    // construct as sub-tree     
	       const PARENT_TREE*const&);          // I: parent tree            
    //--------------------------------------------------------------------------
    // construct a pruned tree:                                                 
    //   a pruned tree gives a coarser representation of the original tree:     
    //   cells for which the function prune(cell_iter) returns true, will be    
    //   mapped to "purely representative" cells. Any cell that contains a      
    //   purely representative cell is called a "pruned" cell.                  
    //   A purely representative cell has no meaningful entries for sub-nodes   
    //   but represents them by their properties (N, ...).                      
    //   On a pruned cell, direct access to ALL soul descendants is impossible. 
    //   The souls and cells of the prunded tree hold the index of the corres-  
    //   ponding soul and cell of the original tree.                            
    //                                                                          
    // The purpose of a pruned tree is to represent the tree of an external     
    // domain in the MPI parallel version of the code: only those cells that    
    // need to be resolved even from outside will not be pruned.                
    //--------------------------------------------------------------------------
    basic_tree(                                    // construct as sub-tree     
	       const basic_tree*const&,            // I: parent tree            
	       bool(*)(const cell_type*const&));   // I: prune cell?            
    //--------------------------------------------------------------------------
    void reset_bodies(                             // exchange bodies           
		      const sbodies*B)             // I: sbodies                
    {
      BSRCES = B; 
#ifdef ALLOW_MPI
      PSRCES = 0;
#endif
    }
#ifdef ALLOW_MPI
    //--------------------------------------------------------------------------
    void reset_bodies(                             // exchange bodies           
		      const pbodies*P)             // I: pbodies                
    { 
      BSRCES = 0;
      PSRCES = P;
    }
#endif
    //--------------------------------------------------------------------------
    // tree-rebuild: use old tree structure                                     
    //--------------------------------------------------------------------------
    // use order of souls only                                                  
    //--------------------------------------------------------------------------
    void build  (                                  // openMP//ized              
		 const int   = 6,                  //[I: N_crit]                
		 const int   = 50);                //[I: max tree depth]        
    //--------------------------------------------------------------------------
    // use full tree structure                                                  
    //--------------------------------------------------------------------------
    void rebuild(                                  //                           
		 const int,                        // I: N_cut                  
		 const int   = 6,                  //[I: N_crit]                
		 const int   = 50);                //[I: max tree depth]        
    //--------------------------------------------------------------------------
    // destructor                                                               
    //--------------------------------------------------------------------------
    ~basic_tree() {
      if(S0) delete[] S0;
      if(C0) delete[] C0;
    }
    //--------------------------------------------------------------------------
    // other public methods                                                     
    //--------------------------------------------------------------------------
    void mark_for_subtree(uint&, uint&) const;
    const sbodies*const&my_sbodies () const { return BSRCES; }
#ifdef ALLOW_MPI
    const pbodies*const&my_pbodies () const { return PSRCES; }
    bool                use_pbodies() const { return PSRCES!=0; }
    bool                use_arrays () const { return BSRCES==0 && PSRCES==0;}
#else
    bool                use_arrays () const { return BSRCES==0;}
#endif
    const int    *const&my_flags   () const { return F; }
    const areal **const&my_pos     () const { return X; }
    bool                use_sbodies() const { return BSRCES!=0; }
    unsigned      const&N_souls    () const { return Ns; }
    unsigned      const&N_cells    () const { return Nc; }
    unsigned      const&N_array    () const { return Nb; }
    int           const&depth      () const { return Dp; }
    //--------------------------------------------------------------------------
    // access to cells via cell_iterators                                       
    //--------------------------------------------------------------------------
  public:
    cell_iterator root        () const { return cell_iterator(this,C0); }
    cell_iterator begin_cells () const { return cell_iterator(this,C0); }
    cell_iterator end_cells   () const { return cell_iterator(this,CN); }
    cell_iterator last_cells  () const { return cell_iterator(this,CN-1); }
    cell_iterator rbegin_cells() const { return cell_iterator(this,CN-1); }
    cell_iterator rend_cells  () const { return cell_iterator(this,C0-1); }
    cell_iterator cell_No(int const&i) const { return cell_iterator(this,C0+i);}
    //--------------------------------------------------------------------------
    // type soul_iterator                                                       
    //--------------------------------------------------------------------------
    typedef soul_type* soul_iterator;
    //--------------------------------------------------------------------------
    // access to souls via soul_iterators                                       
    //--------------------------------------------------------------------------
    soul_iterator const& begin_souls () const { return S0; }
    soul_iterator const& end_souls   () const { return SN; }
    soul_iterator const  soul_No     (int const&i) const { return S0+i; }
    size_t index(const soul_type* const&S) const { return S-S0; }
    //--------------------------------------------------------------------------
    // dump cell and soul data                                                  
    //--------------------------------------------------------------------------
  private:
    void dump_nodes(cell_iterator const&,
		    std::ostream* const&, std::ostream* const&) const;
  public:
    void dump_nodes(std::ostream*, std::ostream*) const;
    //--------------------------------------------------------------------------
  };
  //============================================================================
  static const int SUBTREE = 1<<7;                 // node = node of subtree    
  inline void flag_for_subtree(flag*F) { F->add    (SUBTREE); }
  inline void unflag_subtree  (flag*F) { F->un_set (SUBTREE); }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// macros for looping all souls and cells in a tree                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
// loop the cells down: root first;                                             
//------------------------------------------------------------------------------
#define LoopCellsDown(TREE_TYPE,                   /* type of parent tree    */\
		      TREE_PTER,                   /* pointer to parent tree */\
		      NAME)                        /* name for cells         */\
  for(register TREE_TYPE::cell_iterator	           /* tree::cell_iterator    */\
      NAME  = TREE_PTER->begin_cells();            /* from root cell         */\
      NAME != TREE_PTER->end_cells();              /* until beyond last cell */\
    ++NAME)                                        // get next cell             
//------------------------------------------------------------------------------
#define LoopMyCellsDown(NAME)                      /* name for cells         */\
  for(register cell_iterator	                   /* type of cell           */\
      NAME  = begin_cells();                       /* from root cell         */\
      NAME != end_cells();                         /* until beyond last cell */\
    ++NAME)                                        // get next cell             
//------------------------------------------------------------------------------
// loop the cells up: root last;                                                
//------------------------------------------------------------------------------
#define LoopCellsUp(TREE_TYPE,                     /* type of parent tree    */\
		    TREE_PTER,                     /* pointer to parent tree */\
	            NAME)                          /* name for cell* s       */\
  for(register TREE_TYPE::cell_iterator	           /* type of cell           */\
      NAME  = TREE_PTER->rbegin_cells();           /* from last cell         */\
      NAME != TREE_PTER->rend_cells();             /* until beyond root cell */\
    --NAME)                                        // get previous cell         
//------------------------------------------------------------------------------
#define LoopMyCellsUp(NAME)                        /* name for cells         */\
  for(register cell_iterator	                   /* type of cell           */\
      NAME  = rbegin_cells();                      /* from last cell         */\
      NAME != rend_cells();                        /* until beyond root cell */\
    --NAME)                                        // get previous cell         
//------------------------------------------------------------------------------
// loop the souls down;                                                         
//------------------------------------------------------------------------------
#define LoopSouls(TREE_TYPE,                       /* type of parent tree    */\
		  TREE_PTER,                       /* pointer to parent tree */\
	          NAME)                            /* name for cell* s       */\
  for(register TREE_TYPE::soul_iterator	           /* type of cell           */\
      NAME  = TREE_PTER->begin_souls();            /* from first soul        */\
      NAME != TREE_PTER->end_souls();              /* until beyond last soul */\
    ++NAME)                                        // get next soul             
//------------------------------------------------------------------------------
#define LoopMySouls(NAME)                          /* name for souls         */\
  for(register soul_iterator	                   /* type of soul           */\
      NAME  = begin_souls();                       /* from first soul        */\
      NAME != end_souls();                         /* until beyond last soul */\
    ++NAME)                                        // get next soul             
//------------------------------------------------------------------------------
#define LoopSoulsRange(TREE_TYPE,                  /* type of parent tree    */\
		       TREE_PTER,                  /* pointer to parent tree */\
                       BEGIN,                      /* index: first soul      */\
		       END,                        /* index: beyond last soul*/\
	               NAME)                       /* name for cell* s       */\
  for(register TREE_TYPE::soul_iterator	           /* type of cell           */\
      NAME  = TREE_PTER->soul_No(BEGIN);           /* from begin soul        */\
      NAME != TREE_PTER->soul_No(END);             /* until beyond end soul  */\
    ++NAME)                                        // get next soul             
//------------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// macros for looping the cell and soul children and descendants of a cell    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
// loop the cell children of a given cell                                       
//------------------------------------------------------------------------------
#define LoopCellKids(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
		     NAME)                         /* name for child cells   */\
  for(register CELL_TYPE::cell_child               /* type of child cell     */\
      NAME  = CELL.begin_cell_kids();              /* from first child       */\
      NAME != CELL.end_cell_kids();                /* until beyond last      */\
    ++NAME)                                        // get next child            
//------------------------------------------------------------------------------
// loop the cell children of a given cell starting somewhere                    
//------------------------------------------------------------------------------
#define LoopCellSecd(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
		     CELL_START,                   /* start loop here        */\
		     NAME)                         /* name for child cells   */\
  for(register CELL_TYPE::cell_child               /* type of child cell     */\
      NAME  = CELL_START;                          /* from start child       */\
      NAME != CELL.end_cell_kids();                /* until beyond last      */\
    ++NAME)                                        // get next child            
//------------------------------------------------------------------------------
// loop the soul children of a given cell                                       
//------------------------------------------------------------------------------
#define LoopSoulKids(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
		     NAME)                         /* name for child souls   */\
  for(register CELL_TYPE::soul_child	           /* pter to soul child     */\
      NAME  = CELL.begin_souls();                  /* from first child       */\
      NAME != CELL.end_soul_kids();                /* until beyond last      */\
    ++NAME)                                        // get next child            
//------------------------------------------------------------------------------
// loop the soul children of a given cell starting somewhere                    
//------------------------------------------------------------------------------
#define LoopSoulSecd(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
                     SOUL_START,                   /* start loop here        */\
		     NAME)                         /* name for child souls   */\
  for(register CELL_TYPE::soul_child	           /* type of child cell     */\
      NAME  = SOUL_START;                          /* from start child       */\
      NAME != CELL.end_soul_kids();                /* until beyond last      */\
    ++NAME)                                        // get next child            
//------------------------------------------------------------------------------
// loop the soul descendants of a given cell                                    
//------------------------------------------------------------------------------
#define LoopAllSouls(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
		     NAME)                         /* name for child souls   */\
  for(register CELL_TYPE::soul_child               /* type of child cell     */\
      NAME  = CELL.begin_souls();                  /* from first child       */\
      NAME != CELL.end_soul_desc();                /* until beyond last      */\
    ++NAME)                                        // get next child            
//------------------------------------------------------------------------------
// loop the soul descendants of a given cell starting somewhere                 
//------------------------------------------------------------------------------
#define LoopSecSouls(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
                     SOUL_START,                   /* start loop here        */\
		     NAME)                         /* name for child souls   */\
  for(register CELL_TYPE::soul_child               /* type of child cell     */\
      NAME  = SOUL_START;                          /* from start child       */\
      NAME != CELL.end_soul_desc();                /* until beyond last      */\
    ++NAME)                                        // get next child            
//------------------------------------------------------------------------------
// loop the all, except the last, souls descendants of a given cell             
//------------------------------------------------------------------------------
#define LoopLstSouls(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
		     NAME)                         /* name for child souls   */\
  for(register CELL_TYPE::soul_child               /* type of child cell     */\
      NAME  = CELL.begin_souls();                  /* from start child       */\
      NAME != CELL.last_soul_desc();               /* until last             */\
    ++NAME)                                        // get next child            
////////////////////////////////////////////////////////////////////////////////
#endif                                             // included_tree_h           

// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tree.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines                                                                     |
//                                                                             |
// class basic_leaf                                                            |
// class basic_cell                                                            |
// class oct_tree                                                              |
//                                                                             |
// macros for access to cells & leafs of a tree                                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_tree_h
#define falcON_included_tree_h

#ifndef falcON_included_auxx_h
#  include <public/auxx.h>
#endif
#ifndef falcON_included_flag_h
#  include <public/flag.h>
#endif
#ifndef falcON_included_deft_h
#  include <public/deft.h>
#endif
#ifndef falcON_included_iomanip
#  include <iomanip>
#  define falcON_included_iomanip
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  class  bodies;                                   // forward declaration      
  class ebodies;                                   // forward declaration      
#ifdef falcON_MPI
  class pbodies;                                   // forward declaration      
#endif
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::basic_leaf                                                   //
  //                                                                          //
  // IMPORTANT NOTE                                                           //
  //   Any derived leaf type must NOT define new data members, since the      //
  //   sizeof(leaf)=32 is supposed to be that of basic_leaf.                  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class basic_leaf : public flag {
  private:
    basic_leaf           (const basic_leaf&);      // not implemented           
    basic_leaf& operator=(const basic_leaf&);      // not implemented           
    //--------------------------------------------------------------------------
    // data of class basic_leaf which will be initialized by oct_tree           
    //--------------------------------------------------------------------------
    vect POS;                                      // center of mass = position 
    uint LNK;                                      // associated body index     
    //--------------------------------------------------------------------------
    // data of class basic_leaf which will not be dealt with by oct_tree        
    // the design is such that for the most common case of                      
    //      sizeof(real) = sizeof(void*) = 4                                    
    // we have                                                                  
    //      sizeof(basic_leaf) = 32 = 2*16                                      
    //--------------------------------------------------------------------------
  protected:
    real SCAL;                                     // any scalar (eg. mass)     
    void*PROP;                                     // pointer to more data      
    union {                                        // only of of these, please: 
      int  AUXI;                                   //   auxiliary integer       
      uint AUXU;                                   //   auxiliary unsigned      
      real AUXR;                                   //   auxiliary real          
      void*AUXP;                                   //   auxiliary pointer       
    };
    //--------------------------------------------------------------------------
    // constructor                                                              
    //--------------------------------------------------------------------------
    basic_leaf() {}                                // for optimization, we don't
                                                   // initialize POS & LNK      
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
    void  link    (unsigned  const&L) { LNK  = L; }
    //--------------------------------------------------------------------------
    // const data access                                                        
    //--------------------------------------------------------------------------
  protected:
    vect const&pos   () const { return POS; }
    uint const&mybody() const { return LNK; }
    real const&scalar() const { return SCAL; }
    real const&auxr  () const { return AUXR; }
    uint const&auxu  () const { return AUXU; }
    int  const&auxi  () const { return AUXI; }
    //--------------------------------------------------------------------------
    // non-const data access via members                                        
    //--------------------------------------------------------------------------
  public:
    vect&pos   () { return POS; }
    real&scalar() { return SCAL; }
    real&auxr  () { return AUXR; }
    uint&auxu  () { return AUXU; }
    int &auxi  () { return AUXI; }
    //--------------------------------------------------------------------------
    // const data access via friends                                            
    //--------------------------------------------------------------------------
    friend vect const&pos   (const basic_leaf*const&L) { return L->POS; }
    friend real const&scalar(const basic_leaf*const&L) { return L->SCAL; }
    friend real const&auxr  (const basic_leaf*const&L) { return L->AUXR; }
    friend uint const&auxu  (const basic_leaf*const&L) { return L->AUXU; }
    friend int  const&auxi  (const basic_leaf*const&L) { return L->AUXI; }
    friend uint const&mybody(const basic_leaf*const&L) { return L->LNK; }
    //--------------------------------------------------------------------------
    void set_prop  (real*            const&P) { PROP = P; }
    void copy_link (const basic_leaf*const&L) { LNK = L->LNK; }
    void copy_basic(const basic_leaf*const&L) { POS = L->POS; LNK = L->LNK; }
    //--------------------------------------------------------------------------
    // copy data from body to leaf                                              
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void copy_from_bodies_flag(const bodies_type*const&B) {
      flag::set_to(B->flg(mybody()));
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void copy_from_bodies_pos(const bodies_type*const&B) {
      pos() = B->pos(mybody());
    }
    //--------------------------------------------------------------------------
    // other non-const methods                                                  
    //--------------------------------------------------------------------------
    void set_link_and_pos(const uint&i, vect const&x) {
      link(i);
      pos()=x;
    }
    //--------------------------------------------------------------------------
    void copy(const basic_leaf*const&S)  {
      copy_basic(S);
      set_to_part(S,flag::BODY_FLAGS);
    }
    //--------------------------------------------------------------------------
    void copy_prune(const basic_leaf*const&S) { 
      POS = S->POS;
      LNK = S->LNK; 
      set_to_part(S,flag::BODY_FLAGS);
    }
    //--------------------------------------------------------------------------
    // dump leaf basic data                                                     
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream&o) {
      o<<"     # flg mybody          position";
    }
    //--------------------------------------------------------------------------
    void dump(std::ostream&o) const {
      o<<' '<<std::setw(3)<<flag(*this)
       <<' '<<std::setw(6)<<LNK;
      for(register int d=0; d!=Ndim; ++d)
	o<<' '<<std::setw(8)<<std::setprecision(4)<<POS[d];
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  const int PRUNED      = 1<<10;
  const int PUREREP     = 1<<11;
  const int PRUNE_FLAGS = PRUNED | PUREREP;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // nbdy::class basic_cell                                                   //
  //                                                                          //
  // IMPORTANT NOTE                                                           //
  //   Any derived cell type must NOT define new data members, since the      //
  //   sizeof(cell)=48 is supposed to be that of basic_cell.                  //
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
  private:                                         // 4b: flag                  
    indx  LEVEL;                                   // 2b: tree level            
    indx  OCTANT;                                  // 2b: octant of parent cell 
    indx  NLEAFS;                                  // 2b: # leaf children       
    indx  NCELLS;                                  // 2b: # cell children       
    int   NUMBER;                                  // 4b: # leaf descendants    
    int   FCLEAF;                                  // 4b: index of fst leaf desc
    int   FCCELL;                                  // 4b: index of fst cell kid 
    vect  CENTER;                                  //12b: center of cube        
    //--------------------------------------------------------------------------
    // data of class basic_cell which will not be dealt with by oct_tree        
    //--------------------------------------------------------------------------
  protected:
    union data {                                   // three times either of:    
      void *PTER;                                  //   generic pointer         
      real  SCAL;                                  //   real number             
      uint  NUMB;                                  //   unsigned integer number 
    } AUX1, AUX2, AUX3;                            //48b in basic_cell          
    //--------------------------------------------------------------------------
    // constructor                                                              
    //--------------------------------------------------------------------------
    basic_cell() {}
    //--------------------------------------------------------------------------
    // leaf type to be superseeded in any derived cell                          
    //--------------------------------------------------------------------------
  public:
    typedef basic_leaf leaf_type;
    //--------------------------------------------------------------------------
    // const data access via friends only                                       
    //--------------------------------------------------------------------------
#define CBCC const basic_cell*const&C
    friend indx const&level        (CBCC) { return C->LEVEL; }
    friend indx const&octant       (CBCC) { return C->OCTANT; }
    friend indx const&nleafs       (CBCC) { return C->NLEAFS; }
    friend indx const&ncells       (CBCC) { return C->NCELLS; }
    friend int  const&number       (CBCC) { return C->NUMBER; }
    friend int  const&fcleaf       (CBCC) { return C->FCLEAF; }
    friend int  const&fccell       (CBCC) { return C->FCCELL; }
    friend vect const&center       (CBCC) { return C->CENTER; }
    friend int        ecleaf       (CBCC) { return C->FCLEAF+C->NLEAFS; }
    friend int        ncleaf       (CBCC) { return C->FCLEAF+C->NUMBER; }
    friend int        eccell       (CBCC) { return C->FCCELL+C->NCELLS; }
    friend bool       has_cell_kids(CBCC) { return C->NCELLS != 0; }
    friend bool       has_leaf_kids(CBCC) { return C->NLEAFS != 0; }
    friend bool       is_twig      (CBCC) { return C->NCELLS == 0; }
    friend bool       is_branch    (CBCC) { return C->NLEAFS == 0; }
    friend bool       is_pruned    (CBCC) { return C->is_set(PRUNED); }
    friend bool       is_pure_rep  (CBCC) { return C->is_set(PUREREP); }
#undef CBCC
    //--------------------------------------------------------------------------
    // flag manipulations                                                       
    //--------------------------------------------------------------------------
    void reset_flag() {
      flag::reset();
      flag::add(flag::AL_ACTIVE);
    }
    //--------------------------------------------------------------------------
    void add_flag(int const&f) {
      flag::add(f);
    }
    //==========================================================================
    void reset_active_flag() {
      flag::un_set(flag::ACTIVE);
      flag::add(flag::AL_ACTIVE);
    }
    //--------------------------------------------------------------------------
    void add_active_flag(const basic_leaf*const&L) {
      flag::add_part(L,flag::ACTIVE);
      if(!is_active(L)) flag::un_set(flag::AL_ACTIVE);
    }
    //--------------------------------------------------------------------------
    void add_active_flag(const basic_cell*const&C) {
      flag::add_part(C,flag::ACTIVE);
      if(!al_active(C)) flag::un_set(flag::AL_ACTIVE);
    }
    //--------------------------------------------------------------------------
    // dump cell data                                                           
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream&o) {
      o<<"     # flg lev oct cells ncell leafs nleaf number            center";
    }
    //--------------------------------------------------------------------------
    void dump(std::ostream&o) const {
      o<<' '<<std::setw(3)<<flag(*this)
       <<' '<<std::setw(3)<<LEVEL
       <<' '<<std::setw(3)<<OCTANT;
      if(NCELLS)
	o<<' '<<std::setw(5)<<FCCELL;
      else
	o<<"     -";
      o<<' '<<std::setw(5)<<NCELLS
       <<' '<<std::setw(5)<<FCLEAF
       <<' '<<std::setw(5)<<NLEAFS
       <<' '<<std::setw(6)<<NUMBER;
      for(register int d=0; d!=Ndim; ++d)
	o<<' '<<std::setw(8)<<std::setprecision(4)<<CENTER[d];
    }
    //--------------------------------------------------------------------------
    // other non-const methods                                                  
    //--------------------------------------------------------------------------
  protected:
    void copy_sub(const basic_cell*C) {
      LEVEL  = C->LEVEL;
      OCTANT = C->OCTANT;
      CENTER = C->CENTER;
    }
    //--------------------------------------------------------------------------
    void copy_prune(const basic_cell*C) {
      copy_sub(C);
      NUMBER = C->NUMBER;
      set_to_part(C,flag::ACTIVE_FLAGS);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::oct_tree                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class oct_tree {
    friend class basic_cell_access;
    //--------------------------------------------------------------------------
    oct_tree           (const oct_tree&);          // not implemented           
    oct_tree& operator=(const oct_tree&);          // not implemented           
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
    enum state {
      fresh    = 0,
      re_grown = 1,
      re_used  = 2,
      sub_tree = 4,
      pruned   = 8,
      origins  = sub_tree | pruned
    };
    //--------------------------------------------------------------------------
    enum usage {                                  // tree is used by whom?     |
      un_used  = 0,                               //   not used currently      |
      grav_use = 1,                               //   used by grav_estimator  |
      stsp_use = 2,                               //   used by stsp_estimator  |
      sph_use  = 3                                //   used by sph_estimator   |
    };                                            //                           |
    //--------------------------------------------------------------------------
    // data of class oct_tree                                                   
    //--------------------------------------------------------------------------
  private:
    const  bodies       *BSRCES;                   // pointer to  bodies        
    const ebodies       *ESRCES;                   // pointer to ebodies        
#ifdef falcON_MPI
    const pbodies       *PSRCES;                   // pointer to pbodies        
#endif
    const int            SPFLAG;                   // specific body flag        
    uint                 Ns,Nc;                    // #leafs/cells              
    mutable basic_leaf  *LEAFS;                    // memory for leafs          
    mutable basic_cell  *CELLS;                    // memory for cells          
    mutable real        *RA;                       // table of cell radii       
    vect                 RCENTER;                  // root center               
    union {
      mutable char      *ALLOC;                    // allocated memory          
      mutable unsigned  *DUINT;                    // easy access to 1st few    
    };
    uint                 NALLOC;                   // # bytes allocated         
    state                STATE;                    // tree state                
    mutable  usage       USAGE;                    // tree usage                
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
    void allocate (uint const&, uint const&, uint const&, real const&);
    void set_depth(uint const&);
    //--------------------------------------------------------------------------
    // constructors                                                             
    //                                                                          
    // If the specific body flag is set, only bodies with this flag set will be 
    // used in the tree.                                                        
    //--------------------------------------------------------------------------
    // construct from bodies                                                    
    //--------------------------------------------------------------------------
  public:
    oct_tree(const bodies*const&,                  // I: bodies' flags & pos's  
	     int          const&,                  // I: N_crit                 
	     const vect*const& = 0,                //[I: pre-determined center] 
	     int        const& = Default::MaxDepth,//[I: max tree depth]        
	     int        const& = 0);               //[I: specific flag]         
    //--------------------------------------------------------------------------
    oct_tree(const bodies*const&,                  // I: bodies' flags & pos's  
	     vect         const&,                  // I: x_min                  
	     vect         const&,                  // I: x_max                  
	     int          const&,                  // I: N_crit                 
	     const vect*const& = 0,                //[I: pre-determined center] 
	     int        const& = Default::MaxDepth,//[I: max tree depth]        
	     int        const& = 0);               //[I: specific flag]         
#ifdef falcON_MPI
    //--------------------------------------------------------------------------
    // construct from pbodies                                                   
    //--------------------------------------------------------------------------
    oct_tree(const pbodies*const&,                 // I: bodies' flags & pos's  
	     int           const&,                 // I: N_crit                 
	     const vect*const& = 0,                //[I: pre-determined center] 
	     int        const& = Default::MaxDepth,//[I: max tree depth]        
	     int        const& = 0);               //[I: specific flag]         
    //--------------------------------------------------------------------------
    oct_tree(const pbodies*const&,                 // I: bodies' flags & pos's  
	     vect          const&,                 // I: x_min                  
	     vect          const&,                 // I: x_max                  
	     int           const&,                 // I: N_crit                 
	     const vect*const& = 0,                //[I: pre-determined center] 
	     int        const& = Default::MaxDepth,//[I: max tree depth]        
	     int        const& = 0);               //[I: specific flag]         
#endif
    //--------------------------------------------------------------------------
    // construct from ebodies                                                   
    //--------------------------------------------------------------------------
    oct_tree(const ebodies*const&,                 // I: abodies' flags & pos's 
	     int           const&,                 // I: N_crit                 
	     const vect*const& = 0,                //[I: pre-determined center] 
	     int        const& = Default::MaxDepth,//[I: max tree depth]        
	     int        const& = 0);               //[I: specific flag]         
    //--------------------------------------------------------------------------
    oct_tree(const ebodies*const&,                 // I: abodies' flags & pos's 
	     vect          const&,                 // I: x_min                  
	     vect          const&,                 // I: x_max                  
	     int           const&,                 // I: N_crit                 
	     const vect*const& = 0,                //[I: pre-determined center] 
	     int        const& = Default::MaxDepth,//[I: max tree depth]        
	     int        const& = 0);               //[I: specific flag]         
    //--------------------------------------------------------------------------
    // construct as sub-tree                                                    
    //--------------------------------------------------------------------------
    oct_tree(const oct_tree*const&,                // I: parent tree            
	     int            const&,                // I: flag specif'ing subtree
	     int            const&);               // I: N_crit                 
    //--------------------------------------------------------------------------
    // construct a pruned tree:                                                 
    //   a pruned tree gives a coarser representation of the original tree:     
    //   cells for which the function prune(cell_iter) returns true, will be    
    //   mapped to "purely representative" cells. Any cell that contains a      
    //   purely representative cell is called a "pruned" cell.                  
    //   A purely representative cell has no meaningful entries for sub-nodes   
    //   but represents them by their properties (N, ...).                      
    //   On a pruned cell, direct access to ALL leaf descendants is impossible. 
    //                                                                          
    // The purpose of a pruned tree is to represent the tree of an external     
    // domain in the MPI parallel version of the code: only those cells that    
    // need to be resolved even from outside will not be pruned.                
    //--------------------------------------------------------------------------
    oct_tree(                                      // construct as sub-tree     
	     const oct_tree*const&,                // I: parent tree            
	     bool(*)(const basic_cell*const&));    // I: prune cell?            
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
    basic_leaf*const&FstLeaf()                         const { return LEAFS; }
    basic_cell*const&FstCell()                         const { return CELLS; }
    basic_leaf*      EndLeaf()                         const { return LEAFS+Ns;}
    basic_cell*      EndCell()                         const { return CELLS+Nc;}
    size_t           NoCell (const basic_cell*const&C) const { return C-CELLS; }
    size_t           NoLeaf (const basic_leaf*const&L) const { return L-LEAFS; }
    basic_cell*      CellNo (int const&i)              const { return CELLS+i; }
    basic_leaf*      LeafNo (int const&i)              const { return LEAFS+i; }
    //--------------------------------------------------------------------------
    bool is_re_used    () const { return STATE & re_used; }
    bool is_re_grown   () const { return STATE & re_grown; }
    bool is_fresh      () const { return STATE & fresh; }
    bool is_sub_tree   () const { return STATE & sub_tree; }
    bool is_pruned_tree() const { return STATE & pruned; }
    //--------------------------------------------------------------------------
    bool is_used_for_grav() const { return USAGE == grav_use; }
    bool is_used_for_stsp() const { return USAGE == stsp_use; }
    bool is_used_for_sph () const { return USAGE == sph_use; }
    bool is_not_used     () const { return USAGE == un_used; }
    void mark_grav_usage () const { USAGE = grav_use; }
    void mark_stsp_usage () const { USAGE = stsp_use; }
    void mark_sph_usage  () const { USAGE = sph_use; }
    void mark_un_used    () const { USAGE = un_used; }
    //--------------------------------------------------------------------------
    // tree-rebuild: use old tree structure                                     
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    // use order of leafs only                                                  
    //--------------------------------------------------------------------------
    void                                           //                           
    build(int        const&,                       // I: N_crit                 
	  const vect*const& = 0,                   //[I: pre-determined center] 
	  int        const& = Default::MaxDepth);  //[I: max tree depth]        
    //--------------------------------------------------------------------------
    // use full tree structure                                                  
    //--------------------------------------------------------------------------
#if(0) // code currently broken
    void rebuild(                                  //                           
		 int const&,                       // I: N_cut                  
		 int const&,                       // I: N_crit                 
		 int const& = Default::MaxDepth);  //[I: max tree depth]        
#endif
    //--------------------------------------------------------------------------
    // tree-re-use: just update the leafs' positions                            
    //--------------------------------------------------------------------------
    void reuse();                                  //                           
    //--------------------------------------------------------------------------
    // destructor                                                               
    //--------------------------------------------------------------------------
    ~oct_tree();
    //--------------------------------------------------------------------------
    // other public methods                                                     
    //--------------------------------------------------------------------------
    void mark_for_subtree(int const&, int const&, uint&, uint&) const;
    bool                use_bodies () const { return BSRCES!=0; }
    const  bodies*const&my_bodies  () const { return BSRCES; }
    bool                use_ebodies() const { return ESRCES!=0;}
    const ebodies*const&my_ebodies () const { return ESRCES; }
#ifdef falcON_MPI
    bool                use_pbodies() const { return PSRCES!=0; }
    const pbodies*const&my_pbodies () const { return PSRCES; }
#endif
    int           const&SP_flag    () const { return SPFLAG; }
    unsigned      const&N_leafs    () const { return Ns; }
    unsigned      const&N_cells    () const { return Nc; }
    unsigned      const&depth      () const { return DUINT[2]; }
    unsigned      const&Maxdepth   () const { return DUINT[3]; }
    real const&rad (int  const&l) const { return RA[l]; }
    real const&rad (indx const&l) const { return RA[l]; }
    real const&root_rad   () const { return RA[level(CellNo(0))]; }
    real const&root_radius() const { return RA[level(CellNo(0))]; }
    vect const&root_center() const { return RCENTER; }
    //--------------------------------------------------------------------------
    // type CellIter<>                                                          
    //--------------------------------------------------------------------------
    template<typename CELL> class CellIter {
      //........................................................................
      // types                                                                  
      //........................................................................
#define pC(P) static_cast<CELL*>(P)
#define tI template<typename Cell_Iter>
    public:
      typedef CellIter                 cell_child; // type used in child access 
      typedef typename CELL::leaf_type leaf_type;  // type used in child access 
      typedef leaf_type               *leaf_child; // type used in child access 
      //........................................................................
      // data members                                                           
      //........................................................................
    private:
      const oct_tree*T;                            // pointer to oct_tree       
      CELL          *C;                            // pointer to derived cell   
      //........................................................................
      // construction                                                           
      //........................................................................
    public:
      CellIter(const oct_tree*const&t,
	       CELL          *const&c) : T(t),   C(c)   {}
      //........................................................................
      template<typename cell_type>
      CellIter(const oct_tree*const&t,
	       cell_type     *const&c) : T(t),   C(pC(c)) {}
      //........................................................................
      tI CellIter(Cell_Iter const&I)     : T(I.my_tree()), 
					   C(static_cast<CELL*>
					     (I.c_pter())) {}
      //........................................................................
      CellIter()                         : T(0),   C(0)   {}
      CellIter(int)                      : T(0),   C(0)   {}
      CellIter(CellIter const&I)         : T(I.T), C(I.C) {}
      //........................................................................
      // assigning                                                              
      //........................................................................
      CellIter&operator= (CellIter const&I) {
	T = I.T;
	C = I.C;
	return *this;
      }
      //........................................................................
      tI CellIter&operator= (Cell_Iter const&I) {
	T = I.T;
	C = pC(I.C);
	return *this;
      }
      //........................................................................
      // forward iteration                                                      
      //........................................................................
      CellIter&operator++()                  { ++C; return *this; }
      CellIter operator++(int)               { return CellIter(T,C++); }
      CellIter&operator+=(int const&k)       { C+=k; return *this; }
      CellIter operator+ (int const&k) const { return CellIter(T,C+k); }
      //........................................................................
      // backward iteration                                                     
      //........................................................................
      CellIter&operator--()                  { --C; return *this; }
      CellIter operator--(int)               { return CellIter(T,C--); }
      CellIter&operator-=(int const&k)       { C-=k; return *this; }
      CellIter operator- (int const&k) const { return CellIter(T,C-k); }
      //........................................................................
      // boolean methods                                                        
      //........................................................................
      bool operator==(CellIter  const&I) const { return C == I.C; }
      bool operator!=(CellIter  const&I) const { return C != I.C; }
      bool operator< (CellIter  const&I) const { return C <  I.C; }
      bool operator<=(CellIter  const&I) const { return C <= I.C; }
      bool operator> (CellIter  const&I) const { return C >  I.C; }
      bool operator>=(CellIter  const&I) const { return C >= I.C; }
      tI bool operator==(Cell_Iter const&I) const { return C ==pC(I.c_pter()); }
      tI bool operator!=(Cell_Iter const&I) const { return C !=pC(I.c_pter()); }
      tI bool operator< (Cell_Iter const&I) const { return C < pC(I.c_pter()); }
      tI bool operator<=(Cell_Iter const&I) const { return C <=pC(I.c_pter()); }
      tI bool operator> (Cell_Iter const&I) const { return C > pC(I.c_pter()); }
      tI bool operator>=(Cell_Iter const&I) const { return C >=pC(I.c_pter()); }
#undef pI
#undef pC
      //........................................................................
      // tree and index                                                         
      //........................................................................
      const oct_tree*const&my_tree() const { return T; }
      size_t index() const { return T->NoCell(C); }
      friend size_t index  (CellIter const&I) { return I.index(); }
      //........................................................................
      // conversion to pointer to cell  and  dereferencing to cell*             
      //........................................................................
      operator CELL*const&()           const { return C; }
               CELL*const&operator->() const { return C; }
               CELL*const&c_pter()     const { return C; }
      //........................................................................
      // const access to basic_cell methods via friends                         
      //........................................................................
      friend real const&radius(CellIter const&I) {return I.T->rad(level(I)); }
      //........................................................................
      // access to cell childreen                                               
      //........................................................................
      cell_child begin_cell_kids() const { 
	return cell_child(T, static_cast<CELL*>(T->CellNo(fccell(C)))); }
      cell_child end_cell_kids  () const {
	return cell_child(T, static_cast<CELL*>(T->CellNo(eccell(C)))); }
      //........................................................................
      // access to leaf childreen and leaf  descendants                         
      //........................................................................
      leaf_child begin_leafs   () const {
	return static_cast<leaf_child>(T->LeafNo(fcleaf(C))); }
      leaf_child end_leaf_kids () const {
	return static_cast<leaf_child>(T->LeafNo(ecleaf(C))); }
      leaf_child end_leaf_desc () const {
	return static_cast<leaf_child>(T->LeafNo(ncleaf(C))); }
      leaf_child last_leaf_desc() const { return end_leaf_desc()-1; }
      //........................................................................
      // output: give just the index (an iterator acts like a pointer)          
      //........................................................................
      friend std::ostream& operator<<(std::ostream&o, const CellIter&I) {
	return o<<std::setw(5)<<I.index();
      }
    };
    //--------------------------------------------------------------------------
    // public types, to be superseeded by any derived tree                      
    //--------------------------------------------------------------------------
    typedef basic_cell           cell_type;        // type of tree cell         
    typedef basic_leaf           leaf_type;        // type of tree leaf         
    typedef CellIter<basic_cell> cell_iterator;    // oct_tree::cell_iterator   
    typedef leaf_type*           leaf_iterator;    // oct_tree::leaf_iterator   
    //--------------------------------------------------------------------------
    // access to cells via cell_iterators, to be superseeded by derived tree    
    //--------------------------------------------------------------------------
    cell_iterator root        ()const{ return cell_iterator(this,FstCell()); }
    cell_iterator begin_cells ()const{ return cell_iterator(this,FstCell()); }
    cell_iterator end_cells   ()const{ return cell_iterator(this,EndCell()); }
    cell_iterator last_cells  ()const{ return cell_iterator(this,EndCell()-1); }
    cell_iterator rbegin_cells()const{ return cell_iterator(this,EndCell()-1); }
    cell_iterator rend_cells  ()const{ return cell_iterator(this,FstCell()-1); }
    cell_iterator cell_No     (int const&i) const {
                                       return cell_iterator(this,CellNo(i)); }
    //--------------------------------------------------------------------------
    // access to leafs via leaf_iterators, to be superseeded by derived tree    
    //--------------------------------------------------------------------------
    leaf_iterator const&begin_leafs () const { return FstLeaf(); }
    leaf_iterator       end_leafs   () const { return EndLeaf(); }
    leaf_iterator       leaf_No     (int const&i) const { return LeafNo(i); }
    size_t index(const leaf_type* const&L) const { return NoLeaf(L); }
    //--------------------------------------------------------------------------
    // dump cell and leaf data                                                  
    //--------------------------------------------------------------------------
    template<typename cell_type> void dump_cells(std::ostream&) const;
    template<typename leaf_type> void dump_leafs(std::ostream&) const;
    //--------------------------------------------------------------------------
    // perform some manipulations on the bodies                                 
    //--------------------------------------------------------------------------
    template<typename BodiesManip>
    uint UseBodies(BodiesManip const&BM) const {
      if(BSRCES) return BM(BSRCES);
#ifdef falcON_MPI
      if(PSRCES) return BM(PSRCES);
#endif
      if(ESRCES) return BM(ESRCES);
      error("tree has neither bodies nor arrays");
      return 0u;
    }
    //--------------------------------------------------------------------------
  private:
    uint mark_sub (int const&, int const&, cell_iterator const&, uint &) const;
  };
  //============================================================================
  static const int SUBTREE = 1<<7;                 // node = node of subtree    
  inline void flag_for_subtree(flag*F) { F->add    (SUBTREE); }
  inline void unflag_subtree  (flag*F) { F->un_set (SUBTREE); }
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// macros for looping all leafs and cells in a tree                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
// loop the cells down: root first;                                             
//------------------------------------------------------------------------------
#define LoopCellsDown(CELL_ITER_TYPE,              /* type of cell iterator  */\
		      TREE_PTER,                   /* pointer to parent tree */\
		      NAME)                        /* name for cells         */\
  for(register CELL_ITER_TYPE	                   /* tree::cell_iterator    */\
      NAME((TREE_PTER)->begin_cells());            /* from root cell         */\
      NAME!=(TREE_PTER)->end_cells();              /* until beyond last cell */\
    ++NAME)                                        // get next cell             
//------------------------------------------------------------------------------
// loop the cells up: root last;                                                
//------------------------------------------------------------------------------
#define LoopCellsUp(CELL_ITER_TYPE,                /* type of cell iterator  */\
		    TREE_PTER,                     /* pointer to parent tree */\
	            NAME)                          /* name for cell* s       */\
  for(register CELL_ITER_TYPE	                   /* type of cell           */\
      NAME((TREE_PTER)->rbegin_cells());           /* from last cell         */\
      NAME!=(TREE_PTER)->rend_cells();             /* until beyond root cell */\
    --NAME)                                        // get previous cell         
//------------------------------------------------------------------------------
// loop the leafs down;                                                         
//------------------------------------------------------------------------------
#define LoopLeafs(LEAF_TYPE,                       /* type of leaf           */\
		  TREE_PTER,                       /* pointer to parent tree */\
	          NAME)                            /* name for cell* s       */\
  for(register LEAF_TYPE*	                   /* type of leaf           */\
      NAME =static_cast<LEAF_TYPE*>((TREE_PTER)->begin_leafs()); /*from begin*/\
      NAME!=static_cast<LEAF_TYPE*>((TREE_PTER)->end_leafs());   /*until end */\
    ++NAME)                                        // get next leaf             
//------------------------------------------------------------------------------
#define LoopLeafsRange(LEAF_TYPE,                  /* type of leaf           */\
		       TREE_PTER,                  /* pointer to parent tree */\
                       BEGIN,                      /* index: first leaf      */\
		       END,                        /* index: beyond last leaf*/\
	               NAME)                       /* name for cell* s       */\
  for(register LEAF_TYPE*	                   /* type of leaf           */\
      NAME =static_cast<LEAF_TYPE*>((TREE_PTER)->leaf_No(BEGIN));/*from begin*/\
      NAME!=static_cast<LEAF_TYPE*>((TREE_PTER)->leaf_No(END));  /*until end */\
    ++NAME)                                        // get next leaf             
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// macros for looping the cell and leaf children and descendants of a cell    //
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
// loop the leaf children of a given cell                                       
//------------------------------------------------------------------------------
#define LoopLeafKids(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
		     NAME)                         /* name for child leafs   */\
  for(register CELL_TYPE::leaf_child	           /* pter to leaf child     */\
      NAME  = CELL.begin_leafs();                  /* from first child       */\
      NAME != CELL.end_leaf_kids();                /* until beyond last      */\
    ++NAME)                                        // get next child            
//------------------------------------------------------------------------------
// loop the leaf children of a given cell starting somewhere                    
//------------------------------------------------------------------------------
#define LoopLeafSecd(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
                     LEAF_START,                   /* start loop here        */\
		     NAME)                         /* name for child leafs   */\
  for(register CELL_TYPE::leaf_child	           /* type of child cell     */\
      NAME  = LEAF_START;                          /* from start child       */\
      NAME != CELL.end_leaf_kids();                /* until beyond last      */\
    ++NAME)                                        // get next child            
//------------------------------------------------------------------------------
// loop the leaf descendants of a given cell                                    
//------------------------------------------------------------------------------
#define LoopAllLeafs(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
		     NAME)                         /* name for child leafs   */\
  for(register CELL_TYPE::leaf_child               /* type of child cell     */\
      NAME  = CELL.begin_leafs();                  /* from first child       */\
      NAME != CELL.end_leaf_desc();                /* until beyond last      */\
    ++NAME)                                        // get next child            
//------------------------------------------------------------------------------
// loop the leaf descendants of a given cell starting somewhere                 
//------------------------------------------------------------------------------
#define LoopSecLeafs(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
                     LEAF_START,                   /* start loop here        */\
		     NAME)                         /* name for child leafs   */\
  for(register CELL_TYPE::leaf_child               /* type of child cell     */\
      NAME  = LEAF_START;                          /* from start child       */\
      NAME != CELL.end_leaf_desc();                /* until beyond last      */\
    ++NAME)                                        // get next child            
//------------------------------------------------------------------------------
// loop the all, except the last, leafs descendants of a given cell             
//------------------------------------------------------------------------------
#define LoopLstLeafs(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
		     NAME)                         /* name for child leafs   */\
  for(register CELL_TYPE::leaf_child               /* type of child cell     */\
      NAME  = CELL.begin_leafs();                  /* from start child       */\
      NAME != CELL.last_leaf_desc();               /* until last             */\
    ++NAME)                                        // get next child            
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  template<typename CELL_TYPE> inline
  void oct_tree::dump_cells(std::ostream &o) const {
    CELL_TYPE::dump_head(o);
    o <<'\n';
    for(register cell_iterator Ci=begin_cells(); Ci!=end_cells(); ++Ci) {
      o <<' '<< std::setw(5)<<nbdy::index(Ci);
      static_cast<CELL_TYPE*>
	(static_cast<basic_cell*>(Ci))->dump(o);
      o <<'\n';
    }
    o.flush();
  }
  //----------------------------------------------------------------------------
  template<typename LEAF_TYPE> inline
  void oct_tree::dump_leafs(std::ostream&o) const {
    LEAF_TYPE::dump_head(o);
    o <<'\n';
    for(register leaf_iterator Li=begin_leafs(); Li!=end_leafs(); ++Li) {
      o <<' '<< std::setw(5)<<index(Li);
      static_cast<LEAF_TYPE*>(Li)->dump(o);
      o <<'\n';
    }
    o.flush();
  }
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_tree_h    

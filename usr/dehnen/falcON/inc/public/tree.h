// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/tree.h                                                   
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2000-2005                                                           
///                                                                             
/// \brief  contains definition of classes \a BasicLeaf, \a BasicCell, and      
///         \a OctTree as well as macros for access to cells & leafs of a tree  
///                                                                             
/// \todo   complete doxygen documentation, support parallel code               
///                                                                             
// /////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2000-2005  Walter Dehnen                                       
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
// /////////////////////////////////////////////////////////////////////////////
#ifndef falcON_included_tree_h
#define falcON_included_tree_h

#ifndef falcON_included_body_h
#  include <body.h>
#endif
#ifndef falcON_included_default_h
#  include <public/default.h>
#endif
#ifndef falcON_included_iomanip
#  include <iomanip>
#  define falcON_included_iomanip
#endif
#ifdef falcON_MPI
#  ifndef falcON_included_peano_h
#    include <proper/peano.h>
#  endif
#else
namespace falcON {
  typedef uint8 PeanoMap;
}
#endif
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::BasicLeaf                                                  //
  //                                                                          //
  // IMPORTANT NOTE                                                           //
  //   Any derived leaf type must NOT define new data members, since the      //
  //   sizeof(Leaf)=32 is supposed to be that of BasicLeaf.                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class BasicLeaf {
  private:
    BasicLeaf           (BasicLeaf const&);        // not implemented           
    BasicLeaf& operator=(BasicLeaf const&);        // not implemented           
    //--------------------------------------------------------------------------
    // data of class BasicLeaf which will be initialized by OctTree             
    //--------------------------------------------------------------------------
    vect POS;                                      // center of mass = position 
  protected:
    union {                                        // only one of these, please:
      int         AUXI;                            //   auxiliary integer       
      unsigned    AUXU;                            //   auxiliary unsigned      
      real        AUXR;                            //   auxiliary real          
      void       *AUXP;                            //   auxiliary pointer       
    };
  private:
    flag          FLG;                             // body flag                 
    bodies::index LNK;                             // associated body index     
    //--------------------------------------------------------------------------
    // data of class BasicLeaf which will not be dealt with by OctTree          
    // the design is such that for the most common case of                      
    //      sizeof(real) = sizeof(void*) = 4                                    
    // we have                                                                  
    //      sizeof(BasicLeaf) = 32 = 2*16                                       
    //--------------------------------------------------------------------------
  protected:
    real SCAL;                                     // any scalar (eg. mass)     
    void*PROP;                                     // pointer to more data      
    //--------------------------------------------------------------------------
    // constructor                                                              
    //--------------------------------------------------------------------------
    BasicLeaf() 
#ifdef EBUG
      : SCAL(0), PROP(0), AUXP(0)
#endif
    {}
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
    void  link (bodies::index L) { LNK  = L; }
    //--------------------------------------------------------------------------
    // const data access                                                        
    //--------------------------------------------------------------------------
  protected:
    vect          const&pos   () const { return POS; }
    bodies::index const&mybody() const { return LNK; }
    real          const&scalar() const { return SCAL; }
    real          const&auxr  () const { return AUXR; }
    unsigned      const&auxu  () const { return AUXU; }
    int           const&auxi  () const { return AUXI; }
    flag          const&flg   () const { return FLG; }
    //--------------------------------------------------------------------------
    // non-const data access via members                                        
    //--------------------------------------------------------------------------
  public:
    vect    &pos   () { return POS; }
    real    &scalar() { return SCAL; }
    real    &auxr  () { return AUXR; }
    unsigned&auxu  () { return AUXU; }
    int     &auxi  () { return AUXI; }
    flag    &flg   () { return FLG; }
    //--------------------------------------------------------------------------
    // const data access via friends                                            
    //--------------------------------------------------------------------------
    friend vect          const&pos   (const BasicLeaf*L) { return L->POS; }
    friend real          const&scalar(const BasicLeaf*L) { return L->SCAL; }
    friend real          const&auxr  (const BasicLeaf*L) { return L->AUXR; }
    friend unsigned      const&auxu  (const BasicLeaf*L) { return L->AUXU; }
    friend int           const&auxi  (const BasicLeaf*L) { return L->AUXI; }
    friend bodies::index const&mybody(const BasicLeaf*L) { return L->LNK; }
    friend flag          const&flg   (const BasicLeaf*L) { return L->FLG; }
    //--------------------------------------------------------------------------
    void set_prop  (real*            P) { PROP = P; }
    void copy_link (const BasicLeaf*L) { LNK = L->LNK; }
    void copy_basic(const BasicLeaf*L) { POS = L->POS; LNK = L->LNK; }
    //--------------------------------------------------------------------------
    // flag information                                                         
    //--------------------------------------------------------------------------
    friend bool is_active(const BasicLeaf*L) { return is_active(L->FLG); }
    friend bool is_sph   (const BasicLeaf*L) { return is_sph(L->FLG); }
    friend bool is_sticky(const BasicLeaf*L) { return is_sticky(L->FLG); }
    friend bool is_set   (const BasicLeaf*L, int const&T) {
      return is_set(L->FLG,T); }
    //--------------------------------------------------------------------------
    // copy data from body to leaf                                              
    //--------------------------------------------------------------------------
    void copy_from_bodies_flag(const bodies*B) {
      FLG.set_to_part(B->flg(mybody()), flag::LEAF_FLAGS);
    }
    //--------------------------------------------------------------------------
    void copy_from_bodies_pos(const bodies*B) {
      pos() = B->pos(mybody());
    }
    //--------------------------------------------------------------------------
    // other non-const methods                                                  
    //--------------------------------------------------------------------------
    void set_link_and_pos(bodies::index i, vect const&x) {
      link(i);
      pos()=x;
    }
    //--------------------------------------------------------------------------
    void copy(const BasicLeaf*L)  {
      copy_basic(L);
      FLG.set_to_part(L->FLG,flag::BODY_FLAGS);
    }
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    // dump leaf basic data                                                     
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream&o) {
      o<<"     # flg mybody          position";
    }
    //--------------------------------------------------------------------------
    void dump(std::ostream&o) const {
      o<<' '<<std::setw(3)<<FLG
       <<' '<<std::setw(6)<<LNK;
      for(register int d=0; d!=Ndim; ++d)
	o<<' '<<std::setw(8)<<std::setprecision(4)<<POS[d];
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // falcON::class BasicCell                                                  //
  //                                                                          //
  // IMPORTANT NOTE                                                           //
  //   Any derived cell type must NOT define new data members, since the      //
  //   sizeof(cell)=48 is supposed to be that of BasicCell.                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class BasicCell : public flag {
    //--------------------------------------------------------------------------
    // friendship                                                               
    //--------------------------------------------------------------------------
    friend class BasicCellAccess;                  // allows data access        
    //--------------------------------------------------------------------------
    // data of class BasicCell (write access through class BasicCellAccess)     
    //--------------------------------------------------------------------------
  private:                                         // 4b: flag                  
    uint8    LEVEL;                                // 1b: level in tree         
    uint8    OCTANT;                               // 1b: octant of parent cell 
    PeanoMap PEANO;                                // 1b: Peano-Hilbert map     
    uint8    KEY;                                  // 1b: local Peano key       
    indx     NLEAFS;                               // 2b: # leaf children       
    indx     NCELLS;                               // 2b: # cell children       
    int      NUMBER;                               // 4b: # leaf descendants    
    int      FCLEAF;                               // 4b: index of fst leaf desc
    int      FCCELL;                               // 4b: index of fst cell kid 
    vect     CENTER;                               //12b: center of cube        
    //--------------------------------------------------------------------------
    // data of class BasicCell which will not be dealt with by OctTree          
    //--------------------------------------------------------------------------
  protected:
    vect  POS;                                     //12b: position/cofm         
    real  RAD;                                     // 4b: radius/size/rcrit     
    union data {                                   // three times either of:    
      void     *PTER;                              //   generic pointer         
      real      SCAL;                              //   real number             
      unsigned  NUMB;                              //   unsigned integer number 
    } AUX1, AUX2, AUX3;                            //64b in BasicCell           
    //--------------------------------------------------------------------------
    // constructor                                                              
    //--------------------------------------------------------------------------
    BasicCell() {}
    //--------------------------------------------------------------------------
    // leaf type to be superseeded in any derived cell                          
    //--------------------------------------------------------------------------
  public:
    typedef BasicLeaf leaf_type;
    //--------------------------------------------------------------------------
    // const data access via friends only                                       
    //--------------------------------------------------------------------------
#define CBCC const BasicCell*C
    friend uint8    const&level        (CBCC) { return C->LEVEL; }
    friend uint8    const&octant       (CBCC) { return C->OCTANT; }
#ifdef falcON_MPI
    friend PeanoMap const&peano        (CBCC) { return C->PEANO; }
    friend uint8    const&localkey     (CBCC) { return C->KEY; }
#endif
    friend indx     const&nleafs       (CBCC) { return C->NLEAFS; }
    friend indx     const&ncells       (CBCC) { return C->NCELLS; }
    friend int      const&number       (CBCC) { return C->NUMBER; }
    friend int      const&fcleaf       (CBCC) { return C->FCLEAF; }
    friend int      const&fccell       (CBCC) { return C->FCCELL; }
    friend vect     const&center       (CBCC) { return C->CENTER; }
    friend int            ecleaf       (CBCC) { return C->FCLEAF+C->NLEAFS; }
    friend int            ncleaf       (CBCC) { return C->FCLEAF+C->NUMBER; }
    friend int            eccell       (CBCC) { return C->FCCELL+C->NCELLS; }
    friend bool           has_cell_kids(CBCC) { return C->NCELLS != 0; }
    friend bool           has_leaf_kids(CBCC) { return C->NLEAFS != 0; }
    friend bool           is_twig      (CBCC) { return C->NCELLS == 0; }
    friend bool           is_branch    (CBCC) { return C->NLEAFS == 0; }
#undef CBCC
    //--------------------------------------------------------------------------
    // flag manipulations                                                       
    //--------------------------------------------------------------------------
    void reset_flag() {
      flag::reset();
      flag::add(flag::AL_ACTIVE);
    }
    //--------------------------------------------------------------------------
    void add_flag(int f) {
      flag::add(f);
    }
    //==========================================================================
    void reset_active_flag() {
      flag::un_set(flag::ACTIVE);
      flag::add(flag::AL_ACTIVE);
    }
    //--------------------------------------------------------------------------
    void add_active_flag(const BasicLeaf*L) {
      flag::add_part(flg(L),flag::ACTIVE);
      if(!is_active(flg(L))) flag::un_set(flag::AL_ACTIVE);
    }
    //--------------------------------------------------------------------------
    void add_active_flag(const BasicCell*C) {
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
      o<<' '<<std::setw(3)<< flag(*this)
       <<' '<<std::setw(3)<< falcON::level(this)
       <<' '<<std::setw(3)<< falcON::octant(this);
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
    void copy_sub(const BasicCell*C) {
      LEVEL  = C->LEVEL;
      OCTANT = C->OCTANT;
      CENTER = C->CENTER;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::OctTree                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class OctTree {
    friend class BasicCellAccess;
    //--------------------------------------------------------------------------
    OctTree           (const OctTree&);            // not implemented           
    OctTree& operator=(const OctTree&);            // not implemented           
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
    enum state {
      fresh    = 0,
      re_grown = 1,
      re_used  = 2,
      sub_tree = 4,
      origins  = sub_tree
    };
    //--------------------------------------------------------------------------
    enum usage {                                  // tree is used by whom?     |
      un_used  = 0,                               //   not used currently      |
      grav_use = 1,                               //   used by GravEstimator   |
      stsp_use = 2,                               //   used by stsp_estimator  |
      sph_use  = 3                                //   used by SphEstimator    |
    };                                            //                           |
    //--------------------------------------------------------------------------
    // data of class OctTree                                                    
    //--------------------------------------------------------------------------
  private:
    const   bodies      *BSRCES;                   // pointer to  bodies        
    const   int          SPFLAG;                   // specific body flag        
    unsigned             Ns,Nc;                    // #leafs/cells              
    mutable BasicLeaf   *LEAFS;                    // memory for leafs          
    mutable BasicCell   *CELLS;                    // memory for cells          
    mutable real        *RA;                       // table of cell radii       
    vect                 RCENTER;                  // root center               
    union {
      mutable char      *ALLOC;                    // allocated memory          
      mutable unsigned  *DUINT;                    // easy access to 1st few    
    };
    unsigned             NALLOC;                   // # bytes allocated         
    state                STATE;                    // tree state                
    mutable  usage       USAGE;                    // tree usage                
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
    void allocate (unsigned , unsigned , unsigned , real);
    void set_depth(unsigned);
    //--------------------------------------------------------------------------
    // constructors                                                             
    //                                                                          
    // If the specific body flag is set, only bodies with this flag set will be 
    // used in the tree.                                                        
    //--------------------------------------------------------------------------
    // construct from bodies                                                    
    //--------------------------------------------------------------------------
  public:
    OctTree(const bodies*const&,                   // I: bodies' flags & pos's  
	    int          const&,                   // I: N_crit                 
	    const vect*const& = 0,                 //[I: pre-determined center] 
	    int        const& = Default::MaxDepth, //[I: max tree depth]        
	    int        const& = 0);                //[I: specific flag]         
    //--------------------------------------------------------------------------
    OctTree(const bodies*const&,                   // I: bodies' flags & pos's  
	    vect         const&,                   // I: x_min                  
	    vect         const&,                   // I: x_max                  
	    int          const&,                   // I: N_crit                 
	    const vect*const& = 0,                 //[I: pre-determined center] 
	    int        const& = Default::MaxDepth, //[I: max tree depth]        
	    int        const& = 0);                //[I: specific flag]         
    //--------------------------------------------------------------------------
    // construct as sub-tree                                                    
    //--------------------------------------------------------------------------
    OctTree(const OctTree*const&,                  // I: parent tree            
	    int           const&,                  // I: flag specif'ing subtree
	    int           const&);                 // I: N_crit                 
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
    BasicLeaf*const&FstLeaf()                  const { return LEAFS; }
    BasicCell*const&FstCell()                  const { return CELLS; }
    BasicLeaf*      EndLeaf()                  const { return LEAFS+Ns;}
    BasicCell*      EndCell()                  const { return CELLS+Nc;}
    size_t          NoCell (const BasicCell*C) const { return C-CELLS; }
    size_t          NoLeaf (const BasicLeaf*L) const { return L-LEAFS; }
    BasicCell*      CellNo (int i)             const { return CELLS+i; }
    BasicLeaf*      LeafNo (int i)             const { return LEAFS+i; }
    //--------------------------------------------------------------------------
    bool is_re_used    () const { return STATE & re_used; }
    bool is_re_grown   () const { return STATE & re_grown; }
    bool is_fresh      () const { return STATE & fresh; }
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
    // tree-re-use: just update the leafs' positions                            
    //--------------------------------------------------------------------------
    void reuse();                                  //                           
    //--------------------------------------------------------------------------
    // destructor                                                               
    //--------------------------------------------------------------------------
    ~OctTree();
    //--------------------------------------------------------------------------
    // other public methods                                                     
    //--------------------------------------------------------------------------
    void mark_for_subtree(int, int, unsigned&, unsigned&) const;
    const  bodies*const&my_bodies  () const { return BSRCES; }
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
      const OctTree*T;                             // pointer to OctTree        
      CELL         *C;                             // pointer to derived cell   
      //........................................................................
      // construction                                                           
      //........................................................................
    public:
      CellIter(const OctTree*const&t,
	       CELL         *const&c) : T(t),   C(c)   {}
      //........................................................................
      template<typename cell_type>
      CellIter(const OctTree*const&t,
	       cell_type    *const&c) : T(t),   C(pC(c)) {}
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
      const OctTree*const&my_tree() const { return T; }
      size_t index() const { return T->NoCell(C); }
      friend size_t index  (CellIter const&I) { return I.index(); }
      //........................................................................
      // conversion to pointer to cell  and  dereferencing to cell*             
      //........................................................................
      operator CELL*const&()           const { return C; }
               CELL*const&operator->() const { return C; }
               CELL*const&c_pter()     const { return C; }
      //........................................................................
      // const access to BasicCell methods via friends                          
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
    typedef BasicCell            cell_type;        // type of tree cell         
    typedef BasicLeaf            leaf_type;        // type of tree leaf         
    typedef CellIter<BasicCell>  cell_iterator;    // OctTree::cell_iterator    
    typedef leaf_type*           leaf_iterator;    // OctTree::leaf_iterator    
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
    unsigned UseBodies(BodiesManip const&BM) const {
      return BM(BSRCES);
    }
    //--------------------------------------------------------------------------
  private:
    unsigned mark_sub (int, int, cell_iterator const&, unsigned&) const;
  };
  //////////////////////////////////////////////////////////////////////////////
  static const int SUBTREE = 1<<7;                 // node = node of subtree    
  inline void flag_for_subtree(flag*F) { F->add    (SUBTREE); }
  inline void unflag_subtree  (flag*F) { F->un_set (SUBTREE); }
  inline void flag_for_subtree(BasicLeaf*F) { F->flg().add    (SUBTREE); }
  inline void unflag_subtree  (BasicLeaf*F) { F->flg().un_set (SUBTREE); }
  //////////////////////////////////////////////////////////////////////////////
  template<typename CELL_TYPE> inline
  void OctTree::dump_cells(std::ostream &o) const {
    CELL_TYPE::dump_head(o);
    o <<'\n';
    for(register cell_iterator Ci=begin_cells(); Ci!=end_cells(); ++Ci) {
      o <<' '<< std::setw(5)<<falcON::index(Ci);
      static_cast<CELL_TYPE*>
	(static_cast<BasicCell*>(Ci))->dump(o);
      o <<'\n';
    }
    o.flush();
  }
  //----------------------------------------------------------------------------
  template<typename LEAF_TYPE> inline
  void OctTree::dump_leafs(std::ostream&o) const {
    LEAF_TYPE::dump_head(o);
    o <<'\n';
    for(register leaf_iterator Li=begin_leafs(); Li!=end_leafs(); ++Li) {
      o <<' '<< std::setw(5)<<index(Li);
      static_cast<LEAF_TYPE*>(Li)->dump(o);
      o <<'\n';
    }
    o.flush();
  }
} // namespace falcON {
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
#endif // falcON_included_tree_h

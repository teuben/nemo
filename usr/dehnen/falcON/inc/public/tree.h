// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///
/// \file   inc/public/tree.h
///
/// \author Walter Dehnen
/// \date   2000-2007,2010,2012
///
/// \brief  contains definition of class \a OctTree and macros for access to
///         cells & leafs of a tree
///
// /////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2012  Walter Dehnen
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
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
#ifdef falcON_PROPER
#  ifndef falcON_included_peano_h
#    include <proper/peano.h>
#  endif
#else
namespace falcON {
  typedef uint8_t PeanoMap;
}
#endif
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::OctTree                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class OctTree {
    OctTree           (const OctTree&);            // not implemented           
    OctTree& operator=(const OctTree&);            // not implemented           
    //==========================================================================
    //                                                                          
    // sub-type OctTree::Leaf                                                   
    //                                                                          
    /// represents a body in the tree.                                          
    /// \note any derived type shall NOT add any data (must be of the same      
    /// size), to allow allocation here                                         
    ///                                                                         
    //==========================================================================
  public:
    class Leaf {
    private:
      Leaf           (Leaf const&);
      Leaf& operator=(Leaf const&);
      //------------------------------------------------------------------------
      /// \name data of class Leaf
      //@{
      vect POS;                                 ///< centre of mass = position
    protected:
      /// union to provide auxiliary storage for derived class; not initialized
      union {
	int         AUXI;                       ///< auxiliary integer
	unsigned    AUXU;                       ///< auxiliary unsigned
	real        AUXR;                       ///< auxiliary real
	void       *AUXP;                       ///< auxiliary pointer
      };
    private:
      flags         FLAGS;                      ///< body flags
      bodies::index LINK;                       ///< index of associated body
    protected:
      real SCAL;                                ///< any scalar (eg. mass)
      void*PROP;                                ///< pointer to more data
      //@}
      //------------------------------------------------------------------------
      /// construction: do nothing (unless C-macro EBUG is set
      Leaf() 
#if defined(EBUG) || defined(DEBUG)
	: SCAL(0), PROP(0), AUXP(0)
#endif
	{}
      //------------------------------------------------------------------------
      /// set the link to the associated body
      void link(bodies::index L) { LINK  = L; }
      //------------------------------------------------------------------------
      /// \name const data access via members
      //@{
      vect          const&pos   () const { return POS; }   ///< position
      bodies::index const&mybody() const { return LINK; }  ///< index of my body
      real          const&scalar() const { return SCAL; }  ///< scalar
      real          const&auxr  () const { return AUXR; }  ///< auxiliary scalar
      unsigned      const&auxu  () const { return AUXU; }  ///< auxi'ry unsigned
      int           const&auxi  () const { return AUXI; }  ///< auxiliary int
      flags         const&flag  () const { return FLAGS; } ///< body flags
      //@}
      //------------------------------------------------------------------------
    public:
      /// \name non-const data access via members
      //@{
      vect    &pos   () { return POS; }   ///< position
      real    &scalar() { return SCAL; }  ///< scalar
      real    &auxr  () { return AUXR; }  ///< auxiliary scalar
      unsigned&auxu  () { return AUXU; }  ///< auxiliary unsigned
      int     &auxi  () { return AUXI; }  ///< auxiliary int
      flags   &flag  () { return FLAGS; } ///< body flags
      //@}
      //------------------------------------------------------------------------
      /// \name const data access via friends
      //@{
      /// position
      friend vect const&pos(const Leaf*);
      /// scalar
      friend real const&scalar(const Leaf*);
      /// auxiliary scalar
      friend real const&auxr(const Leaf*);
      /// auxiliary unsigned
      friend unsigned const&auxu(const Leaf*);
      /// auxiliary int
      friend int const&auxi(const Leaf*);
      /// index of associated body
      friend bodies::index const&mybody(const Leaf*);
      /// body flag
      friend flags const&flag(const Leaf*);
      //@}
      //------------------------------------------------------------------------
      /// \name flag information via friends
      //@{
      /// does the body flag indicate activity?
      friend bool is_active(const Leaf*);
      /// does the body flag indicate a SINK particle?
      friend bool is_sink  (const Leaf*);
      /// does the body flag indicate a SPH particle?
      friend bool is_sph   (const Leaf*);
      /// does the body flag indicate a sticky particle?
      friend bool is_sticky(const Leaf*);
      /// is the flag \a f set in the body flag?
      friend bool is_set(const Leaf*, flags::single);
      /// are the flags \a f set in the body flag?
      friend bool are_set(const Leaf*, flags const&);
      //@}
      //------------------------------------------------------------------------
      /// set data member \a PROP
      void set_prop  (real*const&P) {
	PROP = P; }
      //------------------------------------------------------------------------
      void copy_link (const Leaf*L) {
	LINK = L->LINK; }
      void copy_basic(const Leaf*L) {
	POS = L->POS; LINK = L->LINK; }
      void copy      (const Leaf*L) {
	copy_basic(L);
	FLAGS.set_to_part(L->FLAGS,flags::body_flags);
      }
      //------------------------------------------------------------------------
      /// copy body flag from associated body
      void copy_from_bodies_flag(const bodies*B) {
	FLAGS.set_to_part(B->flag(mybody()),flags::leaf_flags);
      }
      //------------------------------------------------------------------------
      /// copy position from associated body
      void copy_from_bodies_pos(const bodies*B) {
	pos() = B->pos(mybody());
      }
      //------------------------------------------------------------------------
      /// set link to associated body and position
      void set_link_and_pos(bodies::index i, vect const&x) {
	link(i);
	pos()=x;
      }
      //------------------------------------------------------------------------
      static void dump_head(std::ostream&o) {
	o<<"#      flag blck in            position";
      }
      //------------------------------------------------------------------------
      void dump(std::ostream&o) const {
	o<<' '<<std::setw(3) << FLAGS
	 <<' '<<std::setw(2) << LINK.no()<<' '<<std::setw(6)<<LINK.in();
	for(int d=0; d!=Ndim; ++d)
	  o<<' '<<std::setw(9)<<std::setprecision(4)<<POS[d];
      }
    }; // class Leaf
    //==========================================================================
    //                                                                          
    // sub-type OctTree::Cell                                                   
    //                                                                          
    /// represents a cell in the tree, size: 64 bytes.                          
    /// \note any derived type shall NOT add any data (must be of the same      
    /// size), to allow allocation here                                         
    ///                                                                         
    //==========================================================================
  public:
    class Cell : public flags {
      friend class CellAccess;            // allows data access in tree building
      //------------------------------------------------------------------------
      /// \name data initialized in OctTree build
      //@{
    private:                                 ///< 4bytes: flag
      uint8_t  LEVEL;                        ///< 1byte : level in tree
      uint8_t  OCTANT;                       ///< 1byte : octant of parent cell
      PeanoMap PEANO;                        ///< 1byte : Peano-Hilbert map
      uint8_t  KEY;                          ///< 1byte : local Peano key
      indx     NLEAFS;                       ///< 2bytes: # leaf children
      indx     NCELLS;                       ///< 2bytes: # cell children
      unsigned NUMBER;                       ///< 4bytes: # leaf descendants
      unsigned FCLEAF;                       ///< 4bytes: index of fst leaf desc
      unsigned FCCELL;                       ///< 4bytes: index of fst cell kid
      unsigned PACELL;                       ///< 4bytes: index of parent cell
      vect     CENTRE;                       ///<12bytes: centre of cube
      // 40 bytes
      //@}
      //------------------------------------------------------------------------
      /// \name data not used by OctTree at all
      //@{
    protected:
      vect  POS;                             ///<12bytes: position/cofm         
      real  RAD;                             ///< 4bytes: radius/size/rcrit     
      // 56 bytes
      /// union \a data to hold either a pointer, a scalar, or an unsigned
      union data {
	void     *PTER;
	real      SCAL;
	unsigned  NUMB;
      } AUX1, AUX2;                          ///< 8bytes for 2 \a data
      // 64 bytes
      //@}
      //------------------------------------------------------------------------
      Cell() {}  ///< construction: do nothing
      //------------------------------------------------------------------------
    public:
      /// value used for invalid parent cell
      static const unsigned INVALID = static_cast<unsigned>(-1);
      /// leaf_type to be superseeded in derived classes
      typedef Leaf leaf_type;
      //------------------------------------------------------------------------
      /// \name const data access via friends only
      //@{
      friend uint8_t    const&level   (const Cell*);
      friend uint8_t    const&octant  (const Cell*);
      friend indx     const&nleafs  (const Cell*);
      friend indx     const&ncells  (const Cell*);
      friend unsigned const&number  (const Cell*);
      friend unsigned const&fcleaf  (const Cell*);
      friend unsigned const&fccell  (const Cell*);
      friend unsigned const&pacell  (const Cell*);
      friend vect     const&centre  (const Cell*);
      friend vect     const&center  (const Cell*);
      friend unsigned ecleaf        (const Cell*);
      friend unsigned ncleaf        (const Cell*);
      friend unsigned eccell        (const Cell*);
      friend bool     has_cell_kids (const Cell*);
      friend bool     has_leaf_kids (const Cell*);
      friend bool     is_twig       (const Cell*);
      friend bool     is_branch     (const Cell*);
#ifdef falcON_PROPER
      friend PeanoMap const&peano   (const Cell*);
      friend uint8_t    const&localkey(const Cell*);
#endif
      //@}
      //------------------------------------------------------------------------
      /// \name flag manipulations
      //@{
      /// reset the flag to flag::AL_ACTIVE
      void reset_flag() {
	flags::operator=(flags::all_active);
      }
      /// add a single flag
      void add_flag(flags::single f) {
	flags::add(f);
      }
      /// add a set of flags
      void add_flags(flags const&f) {
	flags::add(f);
      }
      /// set the acticity flag to flags::all_active
      void reset_active_flag() {
	flags::un_set(flags::active);
	flags::add   (flags::all_active);
      }
      /// set the sink flag to flags::all_sink
      void reset_sink_flag() {
	flags::un_set(flags::sink);
	flags::add   (flags::all_sink);
      }
      /// add the activity flag of a Leaf
      void add_active_flag(const Leaf*L) {
	if( is_active(flag(L)))
	  flags::add(flags::active);
	else
	  flags::un_set(flags::all_active);
      }
      /// add the activity flag of a Cell
      void add_active_flag(const Cell*C) {
	flags::add_part(*C,flags::active);
	if(!al_active(C)) flags::un_set(flags::all_active);
      }
      /// add the sink flag of a Leaf
      void add_sink_flag(const Leaf*L) {
	if( is_sink(flag(L)))
	  flags::add(flags::sink);
	else
	  flags::un_set(flags::all_sink);
      }
      /// add the sink flag of a Cell
      void add_sink_flag(const Cell*C) {
	flags::add_part(*C,flags::sink);
	if(!al_sink(C)) flags::un_set(flags::all_sink);
      }
      //@}
      //------------------------------------------------------------------------
      static void dump_head(std::ostream&o) {
	o<<"#      flag    lev oct paren cells ncell leafs nleaf number"
	 <<"            centre        ";
      }
      //------------------------------------------------------------------------
      void dump(std::ostream&o) const {
	o<<' '<<std::setw(7)<< flags(*this)
	 <<' '<<std::setw(3)<< int(LEVEL)
	 <<' '<<std::setw(3)<< int(OCTANT);
	if(PACELL==INVALID)
	  o<<"     -";
	else
	  o<<' '<<std::setw(5)<<PACELL;
	if(NCELLS)
	  o<<' '<<std::setw(5)<<FCCELL;
	else
	  o<<"     -";
	o<<' '<<std::setw(5)<<NCELLS
	 <<' '<<std::setw(5)<<FCLEAF
	 <<' '<<std::setw(5)<<NLEAFS
	 <<' '<<std::setw(6)<<NUMBER;
	for(int d=0; d!=Ndim; ++d)
	  o<<' '<<std::setw(8)<<std::setprecision(4)<<CENTRE[d];
      }
      //------------------------------------------------------------------------
    protected:
      void copy_sub(const Cell*C) {
	LEVEL  = C->LEVEL;
	OCTANT = C->OCTANT;
	CENTRE = C->CENTRE;
      }
    }; // class Cell
    //--------------------------------------------------------------------------
  private:
    enum state {
      fresh    = 0,
      re_grown = 1,
      re_used  = 2,
      sub_tree = 4,
      origins  = sub_tree
    };
    //--------------------------------------------------------------------------
    enum usage {                                  // tree is used by whom?      
      un_used  = 0,                               //   not used currently       
      grav_use = 1,                               //   used by GravEstimator    
      stsp_use = 2,                               //   used by PartnerEstimator 
      sph_use  = 3                                //   used by SphEstimator     
    };                                            //                            
    //--------------------------------------------------------------------------
    /// \name data
    //@{
    const   bodies    *BSRCES;                   ///< pointer to  bodies        
    const   flags      SPFLAG;                   ///< specific body flag        
    unsigned           Ns,Nc;                    ///< #leafs/cells              
    mutable Leaf      *LEAFS;                    ///< memory for leafs          
    mutable Cell      *CELLS;                    ///< memory for cells          
    mutable real      *RA;                       ///< table of cell radii       
    vect               RCENTRE;                  ///< root centre               
    union {
      mutable char    *ALLOC;                    ///< allocated memory          
      mutable unsigned*DUINT;                    ///< easy access to 1st few    
    };
    unsigned           NALLOC;                   ///< # bytes allocated         
    state              STATE;                    ///< tree state                
    mutable  usage     USAGE;                    ///< tree usage                
    //@}
    //--------------------------------------------------------------------------
    void allocate (unsigned , unsigned , unsigned , real);
    void set_depth(unsigned);
    void mark_for_subtree(flags, int, unsigned&, unsigned&) const;
    //--------------------------------------------------------------------------
  public:
    /// \name constructors
    //@{
    /// construct from scratch given bodies
    /// \param B  (input) bodies
    /// \param Nc (input) Ncrit: minimum # bodies/cell
    /// \param C  (optional input) centre of root cell
    /// \param D  (optional input) maximum tree depth
    /// \param S  (optional input) specific flag: only use bodies flagged
    /// \param xmin (optional input) minimum position in each dimension
    /// \param xmax (optional input) maximum position in each dimension
    /// \param putsinksout put leaf for sink particles in root cell
    OctTree(const bodies*B,
	    int          Nc,
	    const vect*  C= 0,
	    int          D= Default::MaxDepth,
	    flags        S= flags::empty,
	    const vect  *xmin = 0,
	    const vect  *xmax = 0,
	    bool putsinksout  = true);
    //--------------------------------------------------------------------------
    /// construct as sub-tree
    /// \param T  (input) parent tree
    /// \param S  (input) flag used to specify sub-tree
    /// \param Nc (input) Ncrit: minimum # bodies/cell
    OctTree(const OctTree*T,
	    flags         S,
	    int           Nc);
    //@}
    //--------------------------------------------------------------------------
    /// return pointer to 1st leaf
    Leaf*const&FstLeaf()             const { return LEAFS; }
    /// return pointer to 1st cell
    Cell*const&FstCell()             const { return CELLS; }
    /// return pointer to end of leafs
    Leaf*      EndLeaf()             const { return LEAFS+Ns;}
    /// return pointer to end of cells
    Cell*      EndCell()             const { return CELLS+Nc;}
    /// return running number of cell given its pointer
    size_t     NoCell (const Cell*C) const { return C-CELLS; }
    /// return running number of leaf given its pointer
    size_t     NoLeaf (const Leaf*L) const { return L-LEAFS; }
    /// return pointer to cell of given running number
    Cell*      CellNo (int i)        const { return CELLS+i; }
    /// return pointer to leaf of given running number
    Leaf*      LeafNo (int i)        const { return LEAFS+i; }
    //--------------------------------------------------------------------------
    bool is_re_used () const { return STATE & re_used; }  ///< is tree re-used?
    bool is_re_grown() const { return STATE & re_grown; } ///< is tree re-grown?
    bool is_fresh   () const { return STATE & fresh; }    ///< is tree fresh?
    //--------------------------------------------------------------------------
    /// mark for usage with GravEstimator
    void mark_grav_usage () const { USAGE = grav_use; }
    /// is used with GravEstimator?
    bool is_used_for_grav() const { return USAGE == grav_use; }
    //--------------------------------------------------------------------------
    /// mark for usage with PartnerEstimator
    void mark_stsp_usage () const { USAGE = stsp_use; }
    /// is used with PartnerEstimator?
    bool is_used_for_stsp() const { return USAGE == stsp_use; }
    //--------------------------------------------------------------------------
    /// mark for usage with SPHEstimator
    void mark_sph_usage  () const { USAGE = sph_use; }
    /// is used with SPHEstimator?
    bool is_used_for_sph () const { return USAGE == sph_use; }
    //--------------------------------------------------------------------------
    /// mark as unused
    void mark_un_used    () const { USAGE = un_used; }
    /// is unused?
    bool is_not_used     () const { return USAGE == un_used; }
    //--------------------------------------------------------------------------
    /// tree-rebuild: use old tree structure
    /// \param Nc (input) Ncrit: minimum # bodies/cell
    /// \param C  (optional input): root centre
    /// \param D  (optional input): maximum tree depth
    /// \param putsinksout put leaf for sink particles under root
    void build(int        const&Nc,
	       const vect*const&C= 0,
	       int        const&D= Default::MaxDepth,
	       bool putsinksout=true);
    //--------------------------------------------------------------------------
    /// tree-re-use: just update the leafs' positions
    void reuse();
    //--------------------------------------------------------------------------
    /// destructor
    ~OctTree();
    //--------------------------------------------------------------------------
    /// return pointer to bodies (const)
    const  bodies*const&my_bodies  () const { return BSRCES; }
    /// return flag used to select bodies on building
    flags         const&SP_flag    () const { return SPFLAG; }
    /// return number of leafs
    unsigned      const&N_leafs    () const { return Ns; }
    /// return number of cells
    unsigned      const&N_cells    () const { return Nc; }
    /// return tree depth
    unsigned      const&depth      () const { return DUINT[2]; }
    /// return maximum allowed tree depth
    unsigned      const&Maxdepth   () const { return DUINT[3]; }
    /// return radius of cell of given tree level
    real const&rad (int l) const { return RA[l]; }
    /// return radius of root cell
    real const&root_radius() const { return RA[level(CellNo(0))]; }
    /// return centre of root cell
    vect const&root_center() const { return RCENTRE; }
    /// return centre of root cell
    vect const&root_centre() const { return RCENTRE; }
    /// find smallest cell, if any, containing a given position
    /// \param  x position
    /// \return pointer to smallest cell containing x, 0 if x not in root
    const Cell*surrounding_cell(vect const&x) const;
    //==========================================================================
    //                                                                          
    // sub-type OctTree::CellIter<>                                             
    //                                                                          
    /// represents a cell in the tree.                                          
    /// \note to be used instead of a pointer to derived Cell                   
    ///                                                                         
    //==========================================================================
    template<typename CELL> class CellIter {
    public:
      typedef CellIter                 cell_child; ///< used in child access    
      typedef typename CELL::leaf_type leaf_type;  ///< used in child access    
      typedef leaf_type               *leaf_child; ///< used in child access    
      //------------------------------------------------------------------------
    private:
      const OctTree*T;                             ///< pointer to OctTree      
      CELL         *C;                             ///< pointer to derived cell 
      //------------------------------------------------------------------------
      /// constructors
      //@{
    public:
      /// from pointers to tree and cell
      CellIter(const OctTree*const&t, CELL*const&c) : T(t), C(c)   {}
      /// from CellIter<CELL>
      CellIter(CellIter const&I) : T(I.T), C(I.C) {}
      /// from CellIter<DERIVED_CELL>
      template<typename DerivedCell>
      CellIter(CellIter<DerivedCell> const&I) 
	: T( I.my_tree() ), C( static_cast<CELL*>(I.c_pter()) ) {}
      /// from nothing: set pointers to zero
      CellIter() : T(0), C(0) {}
      /// from any int: as from nothing
      explicit CellIter(int) : T(0), C(0) {}
      //@}
      //------------------------------------------------------------------------
      /// assignment
      template<typename SomeCell>
      CellIter&operator= (CellIter<SomeCell> const&I) {
	T = I.my_tree();
	C = static_cast<CELL*>(I.c_pter());
	return *this;
      }
      //------------------------------------------------------------------------
      /// \name forward iteration
      //@{
      /// pre-fix
      CellIter&operator++() {
	++C;
	return *this;
      }
      /// post-fix
      CellIter operator++(int) {
	return CellIter(T,C++);
      }
      /// increment by any amount
      CellIter&operator+=(int k) {
	C+=k; return *this;
      }
      /// '+' operator: get CellIter incremented by any amount
      CellIter operator+ (int k) const {
	return CellIter(T,C+k);
      }
      //@}
      //------------------------------------------------------------------------
      /// \name forward iteration
      //@{
      /// pre-fix
      CellIter&operator--() {
	--C;
	return *this;
      }
      /// post-fix
      CellIter operator--(int) {
	return CellIter(T,C--);
      }
      /// decrement by any amount
      CellIter&operator-=(int k) {
	C-=k;
	return *this;
      }
      /// '-' operator: get CellIter decremented by any amount
      CellIter operator- (int k) const {
	return CellIter(T,C-k);
      }
      //@}
      //------------------------------------------------------------------------
      /// \name boolean methods
      //@{
      /// same cell?
      template<typename SomeCell>
      bool operator==(CellIter<SomeCell> const&I) const {
	return C == static_cast<CELL*>(I.c_pter());
      }
      /// different cell?
      template<typename SomeCell>
      bool operator!=(CellIter<SomeCell> const&I) const {
	return C != static_cast<CELL*>(I.c_pter());
      }
      /// earlier cell?
      template<typename SomeCell>
      bool operator< (CellIter<SomeCell> const&I) const {
	return C <  static_cast<CELL*>(I.c_pter());
      }
      /// earlier or same cell?
      template<typename SomeCell>
      bool operator<=(CellIter<SomeCell> const&I) const {
	return C <= static_cast<CELL*>(I.c_pter());
      }
      /// later cell?
      template<typename SomeCell>
      bool operator> (CellIter<SomeCell> const&I) const {
	return C >  static_cast<CELL*>(I.c_pter());
      }
      /// later or same cell?
      template<typename SomeCell>
      bool operator>=(CellIter<SomeCell> const&I) const {
	return C >= static_cast<CELL*>(I.c_pter());
      }
      /// are we refering to a valid cell?
      bool is_valid() const {
	return T != 0;
      }
      /// are we the root cell
      bool is_root() const {
	return level(C) == 0;
      }
      /// conversion to bool: are we referring to a valid cel?
      operator bool () const {
	return is_valid();
      }
      //@}
      //------------------------------------------------------------------------
      /// \name tree and index of cell in tree
      //@{
      /// return tree of cell
      const OctTree*const&my_tree() const {
	return T;
      }
      /// return index of cell within its tree
      size_t index() const {
	return T->NoCell(C);
      }
      //@}
      //------------------------------------------------------------------------
      /// \name conversion and dereferencing to *cell
      //@{
      operator
      CELL*const&()           const { return C; } ///< conversion
      CELL*const&operator->() const { return C; } ///< dereferencing
      CELL*const&c_pter()     const { return C; } ///< dereferencing
      //@}
      //------------------------------------------------------------------------
      /// radius of cell
      real const&radius() const {return T->rad(level(C)); }
      //------------------------------------------------------------------------
      /// \name access to cell relatives
      //@{
      /// parent cell (if not root)
      CellIter parent() const {
	return pacell(C) == CELL::INVALID?
	  CellIter(0,0) : 
	  CellIter(T, static_cast<CELL*>(T->CellNo(pacell(C)))); }
      /// promote to parent cell
      CellIter& to_parent() {
	if(pacell(C) == CELL::INVALID) { T=0; C=0; }
	else C = static_cast<CELL*>(T->CellNo(pacell(C)));
	return *this;
      }
      /// first cell kid (if any)
      cell_child begin_cell_kids() const { 
	return cell_child(T, static_cast<CELL*>(T->CellNo(fccell(C)))); }
      /// end of cell kids (if any)
      cell_child end_cell_kids  () const {
	return cell_child(T, static_cast<CELL*>(T->CellNo(eccell(C)))); }
      /// cell of given index in the same tree
      CellIter CellNo(unsigned i) const {
	return CellIter(T, static_cast<CELL*>(T->CellNo(i))); }
      //@}
      //------------------------------------------------------------------------
      /// \name access to leaf kids and descendants
      //@{
      /// first leaf kid (always exists)
      leaf_child begin_leafs   () const {
	return static_cast<leaf_child>(T->LeafNo(fcleaf(C))); }
      /// end of leaf kids
      leaf_child end_leaf_kids () const {
	return static_cast<leaf_child>(T->LeafNo(ecleaf(C))); }
      /// end of leaf descendants
      leaf_child end_leaf_desc () const {
	return static_cast<leaf_child>(T->LeafNo(ncleaf(C))); }
      /// last leaf descendant
      leaf_child last_leaf_desc() const { return end_leaf_desc()-1; }
      //@}
      //------------------------------------------------------------------------
      /// does this cell contain a given position?
      bool contains(const vect& pos) {
	const real r = radius();
	return abs(C->centre()[0]-pos[0]) <= r
	  &&   abs(C->centre()[1]-pos[1]) <= r
	  &&   abs(C->centre()[2]-pos[2]) <= r ;
      }
      //------------------------------------------------------------------------
    };// class CellIter<CELL>
    //--------------------------------------------------------------------------
    /// \name public types, to be superseeded by any derived tree
    //@{
    typedef Cell            cell_type;        ///< type of tree cell         
    typedef Leaf            leaf_type;        ///< type of tree leaf         
    typedef CellIter<Cell>  cell_iterator;    ///< OctTree::cell_iterator    
    typedef leaf_type*      leaf_iterator;    ///< OctTree::leaf_iterator    
    //@}
    //--------------------------------------------------------------------------
    /// \name access to cells via cell_iterators
    //@{
    /// root cell
    cell_iterator root        () const {
      return cell_iterator(this,FstCell());
    }
    /// 1st cell == root cell
    cell_iterator begin_cells () const {
      return cell_iterator(this,FstCell());
    }
    /// end of cells
    cell_iterator end_cells   () const {
      return cell_iterator(this,EndCell());
    }
    /// last cell
    cell_iterator last_cells  () const {
      return cell_iterator(this,EndCell()-1);
    }
    /// begin of cells for backward == last cell
    cell_iterator rbegin_cells() const {
      return cell_iterator(this,EndCell()-1);
    }
    /// end of cells for backward == cell before 1st
    cell_iterator rend_cells  () const {
      return cell_iterator(this,FstCell()-1);
    }
    /// cell of given index
    cell_iterator cell_No     (int i) const {
      return cell_iterator(this,CellNo(i));
    }
    //@}
    //--------------------------------------------------------------------------
    /// \name access to leafs via leaf_iterators
    //@{
    /// 1st leaf
    leaf_iterator const&begin_leafs () const {
      return FstLeaf();
    }
    /// end of leafs
    leaf_iterator       end_leafs   () const {
      return EndLeaf();
    }
    /// leaf of given index
    leaf_iterator       leaf_No     (int i) const {
      return LeafNo(i);
    }
    /// index of given leaf
    size_t index(const leaf_type*L) const {
      return NoLeaf(L);
    }
    //@}
    //--------------------------------------------------------------------------
    template<typename cell_type> void dump_cells(std::ostream&) const;
    template<typename leaf_type> void dump_leafs(std::ostream&) const;
    //--------------------------------------------------------------------------
  private:
    unsigned mark_sub (flags, int, cell_iterator const&, unsigned&) const;
  }; // class OctTree {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // inline definitions of friends of class OctTree::Leaf                     //
  // also serve to inject these functions into namespace falcON               //
  //                                                                          //
  // ///////////////////////////////////////////////////////////////////////////
  inline vect const&pos(const OctTree::Leaf*L) { return L->POS; }
  inline real const&scalar(const OctTree::Leaf*L) { return L->SCAL; }
  inline real const&auxr(const OctTree::Leaf*L) { return L->AUXR; }
  inline unsigned const&auxu(const OctTree::Leaf*L) { return L->AUXU; }
  inline int const&auxi(const OctTree::Leaf*L) { return L->AUXI; }
  inline bodies::index const&mybody(const OctTree::Leaf*L) { return L->LINK; }
  inline flags const&flag(const OctTree::Leaf*L) { return L->FLAGS; }
  inline bool is_active(const OctTree::Leaf*L) { return is_active(L->FLAGS); }
  inline bool is_sink(const OctTree::Leaf*L) { return is_sink(L->FLAGS); }
  inline bool is_sph(const OctTree::Leaf*L) { return is_sph(L->FLAGS); }
  inline bool is_sticky(const OctTree::Leaf*L) { return is_sticky(L->FLAGS); }
  inline bool is_set(const OctTree::Leaf*L, flags::single f) {
    return L->FLAGS.is_set(f);
  }
  inline bool are_set(const OctTree::Leaf*L, flags const&f) {
    return L->FLAGS.are_set(f);
  }
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // inline definitions of friends of class OctTree::Cell                     //
  // also serve to inject these functions into namespace falcON               //
  //                                                                          //
  // ///////////////////////////////////////////////////////////////////////////
#ifdef falcON_PROPER
  inline PeanoMap const&peano(const OctTree::Cell*C) { return C->PEANO; }
  inline uint8_t const&localkey(const OctTree::Cell*C) { return C->KEY; }
#endif
  inline uint8_t const&level(const OctTree::Cell*C) { return C->LEVEL; }
  inline uint8_t const&octant(const OctTree::Cell*C) { return C->OCTANT; }
  inline indx const&nleafs(const OctTree::Cell*C) { return C->NLEAFS; }
  inline indx const&ncells(const OctTree::Cell*C) { return C->NCELLS; }
  inline unsigned const&number(const OctTree::Cell*C) { return C->NUMBER; }
  inline unsigned const&fcleaf(const OctTree::Cell*C) { return C->FCLEAF; }
  inline unsigned const&fccell(const OctTree::Cell*C) { return C->FCCELL; }
  inline unsigned const&pacell(const OctTree::Cell*C) { return C->PACELL; }
  inline vect const&center(const OctTree::Cell*C) { return C->CENTRE; }
  inline vect const&centre(const OctTree::Cell*C) { return C->CENTRE; }
  inline unsigned ecleaf(const OctTree::Cell*C) { return C->FCLEAF+C->NLEAFS; }
  inline unsigned ncleaf(const OctTree::Cell*C) { return C->FCLEAF+C->NUMBER; }
  inline unsigned eccell(const OctTree::Cell*C) { return C->FCCELL+C->NCELLS; }
  inline bool has_cell_kids(const OctTree::Cell*C) { return C->NCELLS != 0; }
  inline bool has_leaf_kids(const OctTree::Cell*C) { return C->NLEAFS != 0; }
  inline bool is_twig(const OctTree::Cell*C) { return C->NCELLS == 0; }
  inline bool is_branch(const OctTree::Cell*C) { return C->NLEAFS == 0; }
  // ///////////////////////////////////////////////////////////////////////////
  /// \related falcON::OctTree::CellIter
  /// \name functions taking OctTree::CellIter arguments
  //@{
  /// return index of cell within its tree
  template<typename CELL> inline
  size_t index(OctTree::CellIter<CELL> const&I) { return I.index(); }
  /// radius of cell
  template<typename CELL> inline
  real const&radius(OctTree::CellIter<CELL> const&I) { return I.radius(); }
  /// formatted output: just write index
  template<typename CELL> inline
  std::ostream& operator<<(std::ostream&o, const OctTree::CellIter<CELL>&I) {
    return o<<std::setw(5)<<I.index();
  }
  //@}
  //////////////////////////////////////////////////////////////////////////////
  inline void flag_for_subtree(OctTree::Cell*C) {
    C->add         (flags::subtree); }
  inline void unflag_subtree  (OctTree::Cell*C) {
    C->un_set      (flags::subtree); }
  inline void flag_for_subtree(OctTree::Leaf*L) {
    L->flag().add   (flags::subtree); }
  inline void unflag_subtree  (OctTree::Leaf*L) {
    L->flag().un_set(flags::subtree); }
  //////////////////////////////////////////////////////////////////////////////
  template<typename CELL_TYPE> inline
  void OctTree::dump_cells(std::ostream &o) const {
    CELL_TYPE::dump_head(o);
    o <<'\n';
    for(cell_iterator Ci=begin_cells(); Ci!=end_cells(); ++Ci) {
      o <<' '<< std::setw(5)<<falcON::index(Ci);
      static_cast<CELL_TYPE*>(static_cast<OctTree::Cell*>(Ci))->dump(o);
      o <<'\n';
    }
    o.flush();
  }
  //----------------------------------------------------------------------------
  template<typename LEAF_TYPE> inline
  void OctTree::dump_leafs(std::ostream&o) const {
    LEAF_TYPE::dump_head(o);
    o <<'\n';
    for(leaf_iterator Li=begin_leafs(); Li!=end_leafs(); ++Li) {
      o <<' '<< std::setw(5)<<index(Li);
      static_cast<LEAF_TYPE*>(Li)->dump(o);
      o <<'\n';
    }
    o.flush();
  }
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
falcON_TRAITS(falcON::OctTree,"OctTree");
falcON_TRAITS(falcON::OctTree::Leaf,"OctTree::Leaf");
falcON_TRAITS(falcON::OctTree::Cell,"OctTree::Cell");
namespace WDutils {
  template<typename CELL>
  struct traits<typename falcON::OctTree::CellIter<CELL> > {
    static const char* name() {
      return message("OctTree::CellIter<%s>",traits<CELL>::name());
    }
  };
}
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
  for(CELL_ITER_TYPE	                           /* tree::cell_iterator    */	\
      NAME((TREE_PTER)->begin_cells());            /* from root cell         */\
      NAME!=(TREE_PTER)->end_cells();              /* until beyond last cell */\
    ++NAME)                                        // get next cell             
//------------------------------------------------------------------------------
// loop the cells up: root last;                                                
//------------------------------------------------------------------------------
#define LoopCellsUp(CELL_ITER_TYPE,                /* type of cell iterator  */\
		    TREE_PTER,                     /* pointer to parent tree */\
	            NAME)                          /* name for cell* s       */\
  for(CELL_ITER_TYPE	                           /* type of cell           */\
      NAME((TREE_PTER)->rbegin_cells());           /* from last cell         */\
      NAME!=(TREE_PTER)->rend_cells();             /* until beyond root cell */\
    --NAME)                                        // get previous cell         
//------------------------------------------------------------------------------
// loop the leafs down;                                                         
//------------------------------------------------------------------------------
#define LoopLeafs(LEAF_TYPE,                       /* type of leaf           */\
		  TREE_PTER,                       /* pointer to parent tree */\
	          NAME)                            /* name for cell* s       */\
  for(LEAF_TYPE*	                           /* type of leaf           */\
      NAME =static_cast<LEAF_TYPE*>((TREE_PTER)->begin_leafs()); /*from begin*/\
      NAME!=static_cast<LEAF_TYPE*>((TREE_PTER)->end_leafs());   /*until end */\
    ++NAME)                                        // get next leaf             
//------------------------------------------------------------------------------
#define LoopLeafsRange(LEAF_TYPE,                  /* type of leaf           */\
		       TREE_PTER,                  /* pointer to parent tree */\
                       BEGIN,                      /* index: first leaf      */\
		       END,                        /* index: beyond last leaf*/\
	               NAME)                       /* name for cell* s       */\
  for(LEAF_TYPE*	                           /* type of leaf           */\
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
  for(CELL_TYPE::cell_child                        /* type of child cell     */\
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
  for(CELL_TYPE::cell_child                        /* type of child cell     */\
      NAME  = CELL_START;                          /* from start child       */\
      NAME != CELL.end_cell_kids();                /* until beyond last      */\
    ++NAME)                                        // get next child            
//------------------------------------------------------------------------------
// loop the leaf children of a given cell                                       
//------------------------------------------------------------------------------
#define LoopLeafKids(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
		     NAME)                         /* name for child leafs   */\
  for(CELL_TYPE::leaf_child	                   /* pter to leaf child     */\
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
  for(CELL_TYPE::leaf_child	                   /* type of child cell     */\
      NAME  = LEAF_START;                          /* from start child       */\
      NAME != CELL.end_leaf_kids();                /* until beyond last      */\
    ++NAME)                                        // get next child            
//------------------------------------------------------------------------------
// loop the leaf descendants of a given cell                                    
//------------------------------------------------------------------------------
#define LoopAllLeafs(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
		     NAME)                         /* name for child leafs   */\
  for(CELL_TYPE::leaf_child                        /* type of child cell     */\
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
  for(CELL_TYPE::leaf_child                        /* type of child cell     */\
      NAME  = LEAF_START;                          /* from start child       */\
      NAME != CELL.end_leaf_desc();                /* until beyond last      */\
    ++NAME)                                        // get next child            
//------------------------------------------------------------------------------
// loop all, except the last, leaf descendants of a given cell                  
//------------------------------------------------------------------------------
#define LoopLstLeafs(CELL_TYPE,                    /* type of parent cell    */\
		     CELL,                         /* parent cell            */\
		     NAME)                         /* name for child leafs   */\
  for(CELL_TYPE::leaf_child                        /* type of child cell     */\
      NAME  = CELL.begin_leafs();                  /* from start child       */\
      NAME != CELL.last_leaf_desc();               /* until last             */\
    ++NAME)                                        // get next child            
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_tree_h

// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tree.cc                                                                     |
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
// tree-building  is done in three steps:                                      |
// - root construction                                                         |
// - building of a box-leaf tree                                               |
// - linking to a cell-soul tree                                               |
//                                                                             |
// root constructions:                                                         |
// In each dimension, the mininum and maximum position is found and from them  |
// the center and size of an appropriate root box computed.                    |
//                                                                             |
// building of a box-leaf tree                                                 |
// We first construct a box-leaf tree. The leaf-adding algorithm is used, ie.  |
// the leafs are added one-by-one to the root box (the alternative would be    |
// the box-dividing algorithm, which we found to be slightly less efficient).  |
// The boxes are allocated in blocks, using block_alloc of public/memo.h.      |
// Boxes with less than Ncrit leafs are not divided (ie. we wait until a box   |
// has Ncrit leafs before splitting it).                                       |
//                                                                             |
// linking to a cell-soul tree                                                 |
// The box-leaf tree is mapped to a cell-soul tree, such that all cells that   |
// are contained in a given cell are contiguous in memory. The same holds for  |
// the souls.                                                                  |
//                                                                             |
// Notes                                                                       |
// There are several reasons that make the two-step process of first building  |
// a box-leaf tree and then mapping it to a cell-soul tree worth our while:    |
// - we can arrange sub-cells and sub-souls to be contiguous in memory;        |
//   this implies that looping over sub-souls is very fast (no linked lists    |
//   spawming randomly through memory are used), the immediate child souls     |
//   as well as all the soul descendants may be addressed easily.              |
// - we can build the tree with memory-minimal entities (boxes are smaller     |
//   then cells, leafs smaller than souls), saving CPU time;                   |
// - we can allocate EXACTLY the correct number of cells;                      |
//                                                                             |
// Variants                                                                    |
// When an old tree is already existent, we may employ the fact that the order |
// of the new tree may differ only little. There are two ways to exploit that: |
// - We may just add the leafs to the new tree in the same order as they are   |
//   in the old tree. This ensures that subsequent leafs will fall into the    |
//   same box for the most part, reducing random memory access on the boxes.   |
//   This simple method reduces the costs for tree-building by 50% or more for |
//   large N.                                                                  |
// - We may actually take the sorting of the old tree directly. If we still    |
//   want to have a cell-soul tree with contiguous sub-nodes (and we do), then |
//   we must still go via a box-leaf tree. The resulting code is not           |
//   significantly faster then the much simpler method above and in some cases |
//   actually much slower. It is NOT RECOMMENDED; we keep it for reference.    |
//                                                                             |
// Naming                                                                      |
// Throughout this file, we use the following names:                           |
// body:     iterator through either sbodies or pbodies                        |
// soul:     body-representative in the cell-soul tree.                        |
//           Here, we assume that the template parameter SOUL is a class that  |
//           has been derived from class basic_soul of tree.h .                |
// leaf:     a minimal copy of a body/soul, defined below. A leaf is used      |
//           only in this file for tree construction                           |
// cell:     cell of the tree to be built.                                     |
//           However, here, we only assume that the template parameter CELL    |
//           is a class that has been derived from basic_cell of tree.h.       |
// box:      tree cell reduced to the tree-building specific needs, defined    |
//           below; only used in this file for tree construction               |
// node:     either a leaf or a box; actually, node is base of leaf and box    |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/tree.h>
#include <public/memo.h>
#include <public/grat.h>
#include <public/stic.h>
#ifdef falcON_MPI
#  include <walter/grap.h>
#endif
#ifdef falcON_SPH
#  include <sph/spht.h>
#endif

using namespace nbdy;
////////////////////////////////////////////////////////////////////////////////
#define TEST_TIMING
#undef  TEST_TIMING

#ifdef TEST_TIMING
#  include <ctime>
using std::clock_t;
using std::clock;

static clock_t C_0, C_i;
static double  Tini, Tgro, Tlnk, Ttot;
#endif

////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
#define LoopDims for(register int d=0; d<Ndim; d++)
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // auxiliarty constants                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  const int BIT[8]        = {1,2,4,8,16,32,64,128};// 2^i, i=1,..2^Ndim         
  const int SUBTREECELL   = 1<<8;                  // cell = cell of subtree    
  const int SUBTREEZOMBIE = 1<<9;                  // cell with 1 subtree subcel
  const int SUBTREE_FLAGS = SUBTREE | SUBTREECELL | SUBTREEZOMBIE;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // auxiliarty functions                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#define CFCF const flag*const&F
  inline void flag_as_subtreecell  (flag*F) { F->add   (SUBTREECELL); }
  inline void flag_as_subtreezombie(flag*F) { F->add   (SUBTREEZOMBIE);}
  inline void unflag_subtreecell   (flag*F) { F->un_set(SUBTREECELL); }
  inline void unflag_subtreezombie (flag*F) { F->un_set(SUBTREEZOMBIE);}
  inline void unflag_subtree_flags (flag*F) { F->un_set(SUBTREE_FLAGS);}
  inline bool in_subtree           (CFCF)   { return F->is_set(SUBTREE); }
  inline bool is_subtreecell       (CFCF)   { return F->is_set(SUBTREECELL); }
  inline bool is_subtreezombie     (CFCF)   { return F->is_set(SUBTREEZOMBIE);}
#undef CFCF
  //----------------------------------------------------------------------------
  // This routines returns, in each dimension, the nearest integer to x         
  //----------------------------------------------------------------------------
  inline vect Integer(const vect& x) {
    register vect c = zero;                        // reset return position     
    LoopDims c[d]=int(x[d]+half);                  // find center position      
    return c;                                      // and return it             
  }
  //----------------------------------------------------------------------------
  // in which octant of the cube centred on cen is pos?                         
  //----------------------------------------------------------------------------
  inline int octant(const vect& cen, const vect& pos) {
    register int 
      oct=(pos[0] > cen[0])? 1 : 0;
    if    (pos[1] > cen[1]) oct |= 2;
#if falcON_NDIM==3
    if    (pos[2] > cen[2]) oct |= 4;
#endif
    return oct;                                    // return octant             
  }
  //----------------------------------------------------------------------------
  // is pos inside the cube centred on cen and with radius (=half size) rad?    
  //----------------------------------------------------------------------------
  inline bool contains(const vect& cen,
		       const real& rad,
		       const vect& pos) {
    return abs(cen[0]-pos[0]) <= rad
      &&   abs(cen[1]-pos[1]) <= rad
#if falcON_NDIM==3
      &&   abs(cen[2]-pos[2]) <= rad
#endif
      ;
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // auxiliarty types                                                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // nbdy::class basic_cell_access                                            //
  //                                                                          //
  // any class derived from this one has write access to the tree-specific    //
  // entries of tree cells, which are otherwise not writable.                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class basic_cell_access {
    //--------------------------------------------------------------------------
    // protected types and methods                                              
    //--------------------------------------------------------------------------
  protected:
    static vect &center_(basic_cell*C) { return C->CENTER; }
    static real &radius_(basic_cell*C) { return C->RADIUS; }
    static int  &number_(basic_cell*C) { return C->NUMBER; }
    static short&nsouls_(basic_cell*C) { return C->NSOULS; }
    static short&ncells_(basic_cell*C) { return C->NCELLS; }
    static int  &fcsoul_(basic_cell*C) { return C->FCSOUL; }
    static int  &fccell_(basic_cell*C) { return C->FCCELL; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::sub_tree_builder<TREE,PARENT_TREE>                           //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename TREE, typename PARENT_TREE>
  class sub_tree_builder : private basic_cell_access {
    //--------------------------------------------------------------------------
    // some types                                                               
    //--------------------------------------------------------------------------
    typedef typename PARENT_TREE::cell_iterator       parent_c_iter;
    typedef typename PARENT_TREE::soul_type          *parent_s_pter;
    typedef typename TREE::cell_iterator          cell_iter;
    typedef typename TREE::soul_type             *soul_pter;
    //--------------------------------------------------------------------------
    // all methods are static                                                   
    //--------------------------------------------------------------------------
    static parent_s_pter find_soul(const parent_c_iter&C) {
      // parent cell with one soul in subtree: find that one                    
      LoopCellKids(typename parent_c_iter,C,c)
	if(in_subtree(c)) return find_soul(c);
      LoopSoulKids(typename parent_c_iter,C,s)
	if(in_subtree(s)) return s;
      error("blind ending in subtree cutting");
      return 0;
    }
    //--------------------------------------------------------------------------
    static parent_c_iter find_cell(const parent_c_iter&C) {
      // parent cell with one cell descendant in subtree: find that one         
      LoopCellKids(typename parent_c_iter,C,c) if(in_subtree(c)) {
	if(is_subtreecell  (c)) return c;
	if(is_subtreezombie(c)) return find_cell(c);
      }
      error("blind ending in subtree cutting");
      return C;
    }
    //--------------------------------------------------------------------------
    // public method                                                            
    //--------------------------------------------------------------------------
  public:
    static
    int link(parent_c_iter const&P,                // I:   parent cell to link  
	     cell_iter     const&C,                // I:   current cell         
	     cell_iter          &Cf,               // I/O: index: next free cell
	     soul_pter     const&S0,               // I:   first sub-tree soul  
	     int                &sf)               // I/O: index: next free soul
    {
      register int dep=0;                          // depth                     
      radius_(C) = P->radius();                    // copy radius               
      center_(C) = P->center();                    // copy center               
      nsouls_(C) = 0;                              // reset cell: # soul kids   
      ncells_(C) = 0;                              // reset cell: # cell kids   
      fcsoul_(C) = sf;                             // set cell: sub-souls       
      LoopSoulKids(typename parent_c_iter,P,ps)    // LOOP(soul kids of Pcell) >
	if(in_subtree(ps)) {                       //   IF(soul == subt soul)   
	  S0[sf++].copy(ps);                       //     copy link to body etc 
	  nsouls_(C)++;                            //     increment # soul kids 
	}                                          //   ENDIF                  <
      LoopCellKids(typename parent_c_iter,P,pc)    // LOOP(cell kids of Pcell) >
	if(is_subtreecell(pc) ||                   //   IF(cell == subt cell)   
	   is_subtreezombie(pc) ) ncells_(C)++;    //     count # cell kids     
	else if(in_subtree(pc)) {                  //   ELIF(cell == subt node) 
	  S0[sf++].copy(find_soul(pc));            //     find soul & copy it   
	  nsouls_(C)++;                            //     increment # soul kids 
	}                                          //   ENDIF                  <
      number_(C) = nsouls_(C);                     // # souls >= # soul kids    
      if(ncells_(C)) {                             // IF(cell has cell kids)    
	register int de;                           //   depth of sub-cell       
	register cell_iter Ci=Cf;                  //   remember free cells     
	fccell_(C) = Ci.index();                   //   set cell children       
	Cf += ncells_(C);                          //   reserve children cells  
	LoopCellKids(typename parent_c_iter,P,pc)  //   LOOP(c kids of Pcell)  >
	  if(is_subtreecell(pc)) {                 //     IF(cell == subt cell) 
	    de = link(pc,Ci,Cf,S0,sf);             //       link sub cells      
	    if(de>dep) dep=de;                     //       update depth        
	    number_(C) += number_(Ci++);           //       count soul descends 
	  } else if(is_subtreezombie(pc)) {        //     ELIF(zombie cell)     
	    de = link(find_cell(pc),Ci,Cf,S0,sf);  //       link sub cells      
	    if(de>dep) dep=de;                     //       update depth        
	    number_(C) += number_(Ci++);           //       count soul descends 
	  }                                        //     ENDIF                <
      } else                                       // ELSE                      
	fccell_(C) =-1;                            //   set pter to cell kids   
      return dep+1;                                // return cells depth        
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::tree_pruner<TREE,PARENT_TREE>                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename TREE, typename PARENT_TREE>
  class tree_pruner : private basic_cell_access {
    tree_pruner(tree_pruner const&);               // not implemented           
    //--------------------------------------------------------------------------
    typedef typename TREE::cell_type            cell_type;
    typedef typename TREE::cell_iterator        cell_iter;
    typedef typename TREE::soul_type           *soul_pter;
    typedef typename PARENT_TREE::cell_iterator parent_ci;
    typedef typename PARENT_TREE::soul_type    *parent_sp;
    typedef bool(*pruner)(parent_ci const&);
    //--------------------------------------------------------------------------
    const pruner        prune;
    const PARENT_TREE  *parent;
    int                 nC,nS;
    //--------------------------------------------------------------------------
    void count_nodes(parent_ci const&C) {
      ++nC;                                        // count C                   
      if(!(*prune)(C)) {                           // IF(C is not pruned)       
	nS += C->nsouls();                         //   count soul kids         
	LoopCellKids(typename parent_ci,C,c)       //   LOOP(cell kids)         
	  count_nodes(c);                          //     RECURSIVE call        
      }                                            // ENDIF                     
    }
    //--------------------------------------------------------------------------
    int link(                                      // R:   depth                
	     parent_ci const&P,                    // I:   parent cell to link  
	     cell_iter const&C,                    // I:   cell to link with    
	     cell_iter      &Cf,                   // I/O: next free cell       
	     soul_pter const&S0,                   // I:   souls                
	     int            &Sf)                   // I/O: next free soul       
    {
      register int dep=0;                          // depth                     
      C->copy_prune(P);                            // copy relevant data        
      if((*prune)(P)) {                            // IF(prune cell P)          
	C->add_flag(PRUNE_FLAGS);                  //   mark C as purely repr.  
	nsouls_(C) = 0;                            //   reset # subsouls        
	ncells_(C) = 0;                            //   reset # subcells        
	fcsoul_(C) =-1;                            //   subsouls meaningless    
	fccell_(C) =-1;                            //   subcells meaningless    
      } else {                                     // ELSE(retain original)     
	nsouls_(C) = P->nsouls();                  //   copy # subsouls         
	ncells_(C) = P->ncells();                  //   copy # subcells         
	fcsoul_(C) = Sf;                           //   set sub-souls           
	LoopSoulKids(typename parent_ci,P,ps)      //   LOOP(sub-souls)         
	  S0[Sf++].copy_prune(ps);                 //     copy properties       
	if(P->ncells()) {                          //   IF(P has cell kids)     
	  register int       de,pruned=0;          //     sub-depth             
	  register cell_iter Ci=Cf;                //     remember free cells   
	  fccell_(C) = Ci.index();                 //     set cell: 1st sub-cell
	  Cf += P->ncells();                       //     reserve sub-cells     
	  LoopCellKids(typename parent_ci,P,pc) {  //     LOOP(sub-cells)       
	    de = link(pc,Ci,Cf,S0,Sf);             //       RECURSIVE link      
	    if(de>dep) dep=de;                     //       update depth        
// #if defined(__GNUC__) && __GNUC__ < 3
// 	    pruned |= Ci->is_set(PRUNED);          //       pruning info        
// #else
	    pruned |= is_pruned(Ci);               //       pruning info        
// #endif
	    ++Ci;                                  //       next subcell        
	  }                                        //     END LOOP              
	  if(pruned) C->add_flag(PRUNED);          //     IF(pruned) flag so    
	} else                                     //   ELSE                    
	  fccell_(C) = -1;                         //     set cell: 1st sub-cell
      }                                            // ENDIF                     
      return dep+1;                                // return cell's depth       
    }
    //--------------------------------------------------------------------------
  public:
    tree_pruner(pruner            const&p,
		const PARENT_TREE*const&t) : prune(p), parent(t), nC(0), nS(0)
    {
      count_nodes(t->root());
    }
    //--------------------------------------------------------------------------
    int const& N_cells() const { return nC; }
    int const& N_souls() const { return nS; }
    //--------------------------------------------------------------------------
    int link(                                      // R:   depth                
	     cell_iter const&C0,                   // I:   cell to link with    
	     soul_pter const&S0)                   // I:   souls                
    {
      register cell_iter Cf=C0+1;
      register int       Sf=0;
      return link(parent->root(),C0,Cf,S0,Sf);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::node                                                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class node {
  private:
    vect POS;                                      // position                  
    node(const node&);                             // not implemented           
  public:
    typedef node *node_pter;                       // pointer to node           
    node() {}                                      // default constructor       
    node(const vect&X) : POS(X) {}                 // constructor from position 
    vect&       pos()                  { return POS; }
    vect const &pos()            const { return POS; }
    real const &pos(const int i) const { return POS[i]; }
    friend vect const &pos(const node*N) { return N->POS; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::leaf                                                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  struct leaf : public node {
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
    mutable leaf*NEXT;
    unsigned     LINK;
    //--------------------------------------------------------------------------
    // methods                                                                  
    //--------------------------------------------------------------------------
    void add_to_list(leaf* &LIST, int& COUNTER) {
      NEXT = LIST;
      LIST = this;
      ++COUNTER;
    }
    //--------------------------------------------------------------------------
    void add_to_list(node* &LIST, int& COUNTER) {
      NEXT = static_cast<leaf*>(LIST);
      LIST = this;
      ++COUNTER;
    }
    //--------------------------------------------------------------------------
    template<typename soul_type>
    void  set_up  (const soul_type* const&S) {
      pos() = nbdy::pos(S);
      LINK  = mybody(S);
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void  set_up  (const bodies_type* const&BB, int const&i) {
      LINK  = i;
      pos() = BB->pos(i);
    }
    //--------------------------------------------------------------------------
    void  set_up  (const barrays* const&BB, int const&i) {
      LINK = i;
      pos()[0] = BB->pos_x(i);
      pos()[1] = BB->pos_y(i);
#if falcON_NDIM==3
      pos()[2] = BB->pos_z(i);
#endif
    }
    //--------------------------------------------------------------------------
    template<typename soul_type>
    void  set_soul(soul_type &S) {
      S.set_link_and_pos(LINK,pos());
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::leaf_list                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  struct leaf_list {
    leaf *HEAD;                                    // head of list              
    int   SIZE;                                    // size of list              
    //--------------------------------------------------------------------------
    leaf_list() : HEAD(0), SIZE(0) {}
    leaf_list(leaf* const&L, int const&N) : HEAD(L), SIZE(N) {}
  private:
    leaf_list(const leaf_list&L);                  // not implemented           
    leaf_list& operator=(const leaf_list&);        // not implemented           
    //--------------------------------------------------------------------------
  public:
    bool is_empty() { return SIZE == 0; }          // R: list empty?            
    //--------------------------------------------------------------------------
    void add_leaf(leaf* const&L) {                 // I: leaf to be added       
      L->NEXT = HEAD;                              // set L' next to our list   
      HEAD    = L;                                 // update head of list       
      ++SIZE;                                      // increment size of list    
    }
    //--------------------------------------------------------------------------
    void append(leaf_list const&L) {               // I: list to be appended    
      register leaf* T=L.HEAD;                     // take head of list L       
      while(T->NEXT) T=T->NEXT;                    // find tail of list L       
      T->NEXT = HEAD;                              // let it point to our list  
      HEAD    = L.HEAD;                            // update head of our list   
      SIZE   += L.SIZE;                            // update size               
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  // this macro requires you to close the curly bracket or use macro EndLeafList
#define BeginLeafList(LIST,NAME)		   /* loop elements of list   */\
  for(register leaf*NAME=LIST.HEAD, *NEXT_NAME;    /* register 2 pointers:    */\
      NAME;					   /* loop until current=0    */\
      NAME = NEXT_NAME) { 			   /* set current = next      */\
  NEXT_NAME = NAME->NEXT;                          // get next elemetnt         
#define EndLeafList }                              // close curly brackets      
  // this macro fails if leaf.NEXT is manipulated within the loop               
#define LoopLeafListS(LIST,NAME)		   /* loop elements of list   */\
  for(register leaf* NAME = LIST.HEAD;	           /* currentleaf             */\
      NAME;					   /* loop until current=0    */\
      NAME = NAME->NEXT)                           // set current = next        
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::box                                                          //
  //                                                                          //
  // basic link structure in a box-leaf tree.                                 //
  // a box represents a cube centered on center() with half-size RADIUS       //
  // if N <= Ncrit, it only contains sub-leafs, which are in a linked list    //
  // pointed to by LEAFS.                                                     //
  // if N >  Ncrit, LEAFS=0 and the sub-nodes are in the array OCT of octants //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  struct box : public node {
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
    real       RADIUS;                             // half size of the cube     
    node      *OCT[Nsub];                          // octants                   
    int        TYPE;                               // bitfield: 1=cell, 0=leaf  
    int        NUMBER;                             // number of leafs           
    leaf      *LEAFS;                              // linked list of leafs      
    //--------------------------------------------------------------------------
    // const methods                                                            
    //--------------------------------------------------------------------------
    bool marked_as_box  (int const&i) const {return TYPE & BIT[i]; }
    bool marked_as_leaf (int const&i) const {return !marked_as_box(i); }
    bool octant_is_empty(int const&i) const {return OCT[i]==0; }
    bool octant_is_box  (int const&i) const {return OCT[i] && marked_as_box(i);}
    bool octant_is_leaf (int const&i) const {return OCT[i] && marked_as_leaf(i);}
    vect const&center   ()            const {return pos(); }
    real const&center   (int const&i) const {return pos(i); }
    //--------------------------------------------------------------------------
    int octant(const vect&x) const {               // I: pos  (must be in box)  
      return nbdy::octant(pos(),x);
    }
    //--------------------------------------------------------------------------
    int octant(const leaf*const &L) const {        // I: leaf (must be in box)  
      return nbdy::octant(pos(),L->pos());
    }
    //--------------------------------------------------------------------------
    int octant(const box* const &P) const {        // I: box (must be in box)   
      return nbdy::octant(pos(),P->pos());
    }
    //--------------------------------------------------------------------------
    template<typename cell_type>
    int octant(const cell_type* const&C) const {   // I: cell (must be in box)  
      return nbdy::octant(pos(),center(C));
    }
    //--------------------------------------------------------------------------
    bool contains(const vect&x) const {
      return nbdy::contains(center(),RADIUS,x);
    }
    //--------------------------------------------------------------------------
    bool contains(const leaf* const &L) const {
      return nbdy::contains(center(),RADIUS,L->pos());
    }
    //--------------------------------------------------------------------------
    bool oct_is_empty(int const&i) const {
      return OCT[i] == 0;
    }
    //--------------------------------------------------------------------------
    bool is_twig() const {
      return LEAFS != 0;
    }
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
    box() {}
    //--------------------------------------------------------------------------
    template<typename cell_iter>
    void set_up(cell_iter const&C) {
      pos()     = nbdy::center(C);
      RADIUS    = nbdy::radius(C);
      TYPE      = 0;
      NUMBER    = 0;
      LEAFS     = 0;
      for(register int b=0; b!=Nsub; ++b) OCT[b] = 0;
    }
    //--------------------------------------------------------------------------
    void mark_as_box    (int const&i) { TYPE |=  BIT[i]; }
    void mark_as_leaf   (int const&i) { TYPE &= ~BIT[i]; }
    vect&      center   ()            { return pos(); }
    //--------------------------------------------------------------------------
    box& reset_octants() {
      for(register node**P=OCT; P!=OCT+Nsub; ++P) *P = 0;
      return *this;
    }
    //--------------------------------------------------------------------------
    box& reset() {
      TYPE      = 0;
      NUMBER    = 0;
      LEAFS     = 0;
      reset_octants();
      return *this;
    }
    //--------------------------------------------------------------------------
    void addleaf_to_list(leaf* const&L) {          // add leaf L to linked list 
      L->add_to_list(LEAFS, NUMBER);               // add to linked list        
    }
    //--------------------------------------------------------------------------
    void addleaf_to_octs(leaf* const &L) {         // add leaf L to octants     
      register int b = octant(L);                  // find appropriate octant   
      OCT[b] = L;                                  // fill into octant          
      mark_as_leaf(b);                             // mark octant as leaf       
      ++NUMBER;                                    // increment number          
    }
    //--------------------------------------------------------------------------
    void addbox_to_octs(box *const &P) {           // add box P to octants      
      register int b = octant(P);                  // find appropriate octant   
      OCT[b] = P;                                  // fill into octant          
      mark_as_box(b);                              // mark octant as box        
      NUMBER += P->NUMBER;                         // increment number          
    }
    //--------------------------------------------------------------------------
    void shrink_to_octant(int const&i) {           // shrink to octant i        
      RADIUS *= half;                              // reduce radius             
      if(i&1) center()[0]+= RADIUS; else center()[0]-= RADIUS;
      if(i&2) center()[1]+= RADIUS; else center()[1]-= RADIUS;
#if falcON_NDIM==3
      if(i&4) center()[2]+= RADIUS; else center()[2]-= RADIUS;
#endif
    }
    //--------------------------------------------------------------------------
    void expand_to_octant(int const&i) {           // expand in direction i     
      if(i&1) center()[0]+= RADIUS; else center()[0]-= RADIUS;
      if(i&2) center()[1]+= RADIUS; else center()[1]-= RADIUS;
#if falcON_NDIM==3
      if(i&4) center()[2]+= RADIUS; else center()[2]-= RADIUS;
#endif
      RADIUS *= two;                               // increase radius           
    }
    //--------------------------------------------------------------------------
    // const friends                                                            
    //--------------------------------------------------------------------------
    friend vect      &center (      box*const&B) {return B->center();  }
    friend vect const&center (const box*const&B) {return B->center();  }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::estimate_N_alloc                                             //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class estimate_N_alloc {
  private:
    const size_t &Nleafs;
    const size_t &Nsofar;
  public:
    estimate_N_alloc(const size_t&a, const size_t&b) : Nleafs(a), Nsofar(b) {}
    size_t operator() (const size_t&Nused) {
      register real x = Nused*(real(Nleafs)/real(Nsofar)-one);
      return size_t(x+4*sqrt(x)+16);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
#define PLEAF static_cast<leaf*>
#define PBOX  static_cast<box* >
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::box_leaf_tree                                                //
  //                                                                          //
  // for building of a box-leaf tree by addleaf()                             //
  // for linking of the box-leaf tree to a cell-soul tree by link_cells()     //
  // does not itself allocate the leafs.                                      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class box_leaf_tree : private basic_cell_access {
    box_leaf_tree           (const box_leaf_tree&);// not implemented           
    box_leaf_tree& operator=(const box_leaf_tree&);// not implemented           
    //--------------------------------------------------------------------------
    // data of class box_leaf_tree                                              
    //--------------------------------------------------------------------------
    int                NCRIT;                      // Ncrit                     
    int                NCUT, NCUT_FAC, FACH;       // Ncut, Ncut*fac/2          
    int                DMAX, DEPTH;                // max/actual tree depth     
    size_t             NLEAFS;                     // # leafs (to be) added     
    block_alloc<box>  *BM;                         // allocator for boxes       
  protected:
    box               *P0;                         // root of box-leaf tree     
    //--------------------------------------------------------------------------
    // protected methods are all inlined                                        
    //--------------------------------------------------------------------------
    inline box* new_box(size_t const&nl) {
      return &(BM->new_element(estimate_N_alloc(NLEAFS,nl))->reset());
    }
    //--------------------------------------------------------------------------
    inline box* make_parbox(                       // R: new box                
			    box*   const&B,        // I: daughter box           
			    int    const&i,        // I: direction of extension 
			    size_t const&nl)       // I: # leafs added sofar    
      // provides a new box, whose i th octant is box B                         
    {
      register box* P = new_box(nl);               // get box off the stock     
      P->RADIUS   = B->RADIUS;                     // copy radius of parent     
      P->center() = B->center();                   // copy center of parent     
      P->expand_to_octant(i);                      // expand in direction       
      P->addbox_to_octs(B);                        // add B as sub-box          
      return P;                                    // return new box            
    }
    //--------------------------------------------------------------------------
    inline box* make_subbox(                       // R: new box                
			    const box*const&B,     // I: parent box             
			    int       const&i,     // I: parent box's octant    
			    size_t    const&nl)    // I: # leafs added sofar    
      // provides a new empty box in the i th octant of B                       
    {
      register box* P = new_box(nl);               // get box off the stock     
      P->RADIUS   = B->RADIUS;                     // copy radius of parent     
      P->center() = B->center();                   // copy center of parent     
      P->shrink_to_octant(i);                      // shrink to correct octant  
      return P;                                    // return new box            
    }
    //--------------------------------------------------------------------------
    inline box* make_subbox_1(                     // R: new box                
			      const box*const&B,   // I: parent box             
			      int       const&i,   // I: parent box's octant    
			      leaf*     const&L,   // I: leaf of octant         
			      size_t    const&nl)  // I: # leafs added sofar    
      // provides a new box in the i th octant of B and containing leaf L       
      // requires that NCRIT == 1                                               
    {
      register box* P = make_subbox(B,i,nl);       // make new sub-box          
      P->addleaf_to_octs(L);                       // add leaf to its octant    
      return P;                                    // return new box            
    }
    //--------------------------------------------------------------------------
    inline box* make_subbox_N(                     // R: new box                
			      const box*const&B,   // I: parent box             
			      int       const&i,   // I: parent box's octant    
			      leaf*     const&L,   // I: leaf of octant         
			      size_t    const&nl)  // I: # leafs added sofar    
      // provides a new box in the i th octant of B and containing leaf L       
      // requires that NCRIT >  1                                               
    {
      register box* P = make_subbox(B,i,nl);       // make new sub-box          
      P->addleaf_to_list(L);                       // add old leaf to list      
      return P;                                    // return new box            
    }
    //--------------------------------------------------------------------------
    template<typename cell_iter>
    inline box* make_cellbox(                      // R: new box                
			     cell_iter const&C,    // I: associated cell        
			     size_t    const&nl)   // I: # leafs added sofar    
      // provides a new empty box at the place of cell C                        
    {
      register box* P = new_box(nl);               // get box off the stock     
      P->set_up(C);                                // copy cell properties      
      return P;                                    // return new box            
    }
    //--------------------------------------------------------------------------
    inline box* make_copybox(                      // R: new box                
			     box    const&Po,      // I: box to copy            
			     size_t const&nl) {    // I: # leafs added sofar    
      // provides a new box with same properties as box Po                      
      return &((new_box(nl))->operator= (Po));
    }
    //--------------------------------------------------------------------------
#define TREE_TOO_DEEP							       \
    error("exceeding maxium tree depth of %d\n                 "	       \
	  "(presumably more than Ncrit=%d bodies have a common position)",     \
	  DMAX,NCRIT);
    //--------------------------------------------------------------------------
    inline void split_box(box*         P,          // I: box to be splitted     
			  register int l,          // I: level                  
			  size_t const&nl)         // I: # leafs added so far   
      // This routines splits a box:                                            
      // The leafs in the linked list are sorted into octants. Octants with one 
      // leaf will just hold that leaf, octants with many leafs will be boxes   
      // with the leafs in the linked list.                                     
      // If all leafs happen to be in just one octant, the process is repeated  
      // recursively on the box of this octant.                                 
    {
      register int   NUM[Nsub];                    // array with number of leafs
      register int   b,ne;                         // octant & counters         
      register box  *sub=0;                        // current sub-box           
      register leaf *Li, *Ln;                      // current & next leaf       
      do {                                         // DO until # octants > 1    
	if(l>DMAX) TREE_TOO_DEEP                   // IF level>DMAX ERROR       
	for(b=0; b<Nsub; b++) NUM[b] = 0;          //   reset counters          
	for(Li=P->LEAFS; Li; Li=Ln) {              //   LOOP linked list        
	  Ln = Li->NEXT;                           //     next leaf in list     
	  b  = P->octant(Li);                      //     octant of current leaf
	  Li->add_to_list(P->OCT[b], NUM[b]);      //     add leaf to list [b]  
	}                                          //   END LOOP                
	P->LEAFS = 0;                              //   reset list of sub-leafs 
	for(ne=b=0; b<Nsub; b++) if(NUM[b]) {      //   LOOP non-empty octs     
	  ne++;                                    //     count them            
	  if(NUM[b]>1) {                           //     IF many leafs         
	    sub = make_subbox(P,b,nl);             //       make sub-box        
	    sub->LEAFS = PLEAF(P->OCT[b]);         //       assign sub-box's    
	    sub->NUMBER= NUM[b];                   //       leaf list & number  
	    P->OCT[b]= sub;                        //       set octant=sub-box  
	    P->mark_as_box(b);                     //       mark octant as box  
	  }                                        //     ENDIF                 
	}                                          //   END LOOP                
	P = sub;                                   //   set current box=sub-box 
	++l;                                       //   increment depth counter 
      } while(ne==1);                              // WHILE only 1 octant       
    }
    //--------------------------------------------------------------------------
    inline void addleaf_1(                         // add single leaf to box    
			  box*   const&base,       // I: base box to add to     
			  leaf*        Li,         // I: leaf to add            
			  size_t const&nl,         // I: # leafs added sofar    
			  int    const&l0)         // I: level of base box      
      // This routine will leave the box-leaf tree such that twig boxes contain 
      // at most 1 leafs                                                        
    {
      register int    b,l;                         // octant, level             
      register leaf  *Lo;                          // actual leaf loaded        
      register box   *P=base;;                     // actual box                
      register node **oc;                          // actual sub-node           
      for(l=l0; l<DMAX; l++) {                     // LOOP over boxes           
	b  = P->octant(Li);                        //   leaf's octant           
	oc = P->OCT+b;                             //   pointer to octant       
	P->NUMBER++;                               //   increment number        
	if((*oc)==0) {                             //   IF octant empty         
	  *oc = Li;                                //     assign leaf to it     
	  P->mark_as_leaf(b);                      //     mark octant as leaf   
	  return;                                  // <-  DONE with this leaf   
	} else if(P->marked_as_leaf(b)) {          //   ELIF octant=leaf        
	  Lo = PLEAF(*oc);                         //     get old leaf          
	  P->mark_as_box(b);                       //     mark octant as box    
	  P = make_subbox_1(P,b,Lo,nl);            //     create sub-box        
	  *oc = P;                                 //     assign sub-box to oc  
	} else                                     //   ELSE octant=box         
	  P = PBOX(*oc);                           //     set current box       
      }                                            // END LOOP                  
      if(l>=DMAX) TREE_TOO_DEEP                    // IF level>DMAX ERROR       
    }
    //--------------------------------------------------------------------------
    inline void addleaf_N(                         // add single leaf to box    
			  box*   const&base,       // I: base box to add to     
			  leaf*        Li,         // I: leaf to add            
			  size_t const&nl,         // I: # leafs added sofar    
			  int    const&l0)         // I: level of base box      
      // This routine will leave the box-leaf tree such that twig boxes contain 
      // at most NCRIT leafs, where NCRIT > 1 is required (else use addleaf_1)  
    {
      register int    b,l;                         // octant, level             
      register leaf  *Lo;                          // actual leaf loaded        
      register box   *P=base;;                     // actual box                
      register node **oc;                          // actual sub-node           
      for(l=l0; l<DMAX; l++) {                     // LOOP over boxes           
	if(P->is_twig()) {                         //   IF box == twig          
	  P->addleaf_to_list(Li);                  //     add leaf to list      
	  if(P->NUMBER>NCRIT) split_box(P,l,nl);   //     IF(N > NCRIT) split   
	  return;                                  //     DONE with this leaf   
	} else {                                   //   ELIF box==branch        
	  b  = P->octant(Li);                      //     leaf's octant         
	  oc = P->OCT+b;                           //     pointer to octant     
	  P->NUMBER++;                             //     increment number      
	  if((*oc)==0) {                           //     IF octant empty       
	    *oc = Li;                              //       assign leaf to it   
	    P->mark_as_leaf(b);                    //       mark octant as leaf 
	    return;                                // <-    DONE with this leaf 
	  } else if(P->marked_as_leaf(b)) {        //     ELIF octant=leaf      
	    Lo = PLEAF(*oc);                       //       get old leaf        
	    P->mark_as_box(b);                     //       mark octant as box  
	    P = make_subbox_N(P,b,Lo,nl);          //       create sub-box      
	    *oc = P;                               //       assign sub-box to oc
	  } else                                   //     ELSE octant=box       
	    P = PBOX(*oc);                         //       set current box     
	}                                          //   ENDIF                   
      }                                            // END LOOP                  
      if(l>=DMAX) TREE_TOO_DEEP                    // IF level>DMAX  ERROR      
    }
    //--------------------------------------------------------------------------
    inline void addleaf(                           // add single leaf to box    
			box*   const&base,         // I: base box to add to     
			leaf*        Li,           // I: leaf to add            
			size_t const&nl,           // I: # leafs added sofar    
			int    const&l0)           // I: level of base box      
      // This routine will leave the box-leaf tree such that twig boxes contain 
      // at most NCRIT leafs.                                                   
    {
      if(NCRIT > 1) addleaf_N(base,Li,nl,l0);
      else          addleaf_1(base,Li,nl,l0);
    }
    //--------------------------------------------------------------------------
    template<typename cell_iter>
    inline box* unlink(                            // R: box; 0 if box empty    
		       cell_iter const &C,         // I: actual cell            
		       int       const &de,        // I: depth of cell          
		       size_t          &Na,        // I/O: # leafs added sofar  
		       leaf*           &Lf,        // I/O: free leafs           
		       leaf_list       &Lp)        // I/O: list: 'lost' leafs   
    {
      // map cell-soul tree back to box-leaf tree:                              
      // cells with N <= N_cut are grown from scratch of their containing souls 
      // larger cells accumulate a list of "lost" souls                         
      // NOTE:     this version ASSUMES that all souls are still in tree        
      register box* P = make_cellbox(C,Na);        // get new box = cell        
      if(number(C) > NCUT_FAC) {                   // IF N > N_cut * fac/2      
	register leaf_list Lo;                     //   list: leafs lost from C 
	LoopSoulKids(typename cell_iter,C,Si) {    //   LOOP soul kids of C     
	  Lf->set_up(Si);                          //     map soul -> leaf      
	  Lo.add_leaf(Lf);                         //     add to list Lo        
	  ++Lf;                                    //     increment leaf pter   
	}                                          //   END LOOP                
	register box* Pi;                          //   pter to child box       
	LoopCellKids(typename cell_iter,C,Ci) {    //   LOOP cell kids of C     
	  Pi = unlink(Ci,de+1,Na,Lf,Lo);           //     RECURSIVE call        
	  if(Pi) P->addbox_to_octs(Pi);            //     add box to octants    
	}                                          //   END LOOP                
	if(Lo.SIZE   > NCUT        &&              //   IF enough lost leafs    
	   number(C) > FACH*Lo.SIZE    ) {         //      but much less than N 
	  BeginLeafList(Lo,Li)                     //     LOOP list of lost L   
	    if(P->contains(Li))                    //       IF leaf in box      
	      addleaf(P,Li,Na++,de);               //         add leaf to box   
	    else                                   //       ELSE(still lost)    
	      Lp.add_leaf(Li);                     //         add leaf to Lp    
	  EndLeafList                              //     END LOOP              
	} else if(!Lo.is_empty())                  //   ELSE (few lost leafs)   
	  Lp.append(Lo);                           //     append list Lo to Lp  
      } else if(number(C) > NCUT) {                // ELIF N > N_cut            
	LoopSoulKids(typename cell_iter,C,Si) {    //   LOOP soul kids of C     
	  Lf->set_up(Si);                          //     map soul -> leaf      
	  Lp.add_leaf(Lf);                         //     add to list Lp        
	  ++Lf;                                    //     increment leaf pter   
	}                                          //   END LOOP                
	register box* Pi;                          //   pter to child box       
	LoopCellKids(typename cell_iter,C,Ci) {    //   LOOP cell kids of C     
	  Pi = unlink(Ci,de+1,Na,Lf,Lp);           //     RECURSIVE call        
	  if(Pi) P->addbox_to_octs(Pi);            //     add box to octants    
	}                                          //   END LOOP                
      } else {                                     // ELSE(N <= N_cut)          
	LoopAllSouls(typename cell_iter,C,Si) {    //   LOOP all souls in cell  
	  Lf->set_up(Si);                          //     initialize leaf       
	  if(P->contains(Lf))                      //     IF(leaf in box)       
	    addleaf(P,Lf,Na++,de);                 //       add leaf to box     
	  else                                     //     ELSE(not in box)      
	    Lp.add_leaf(Lf);                       //       add to list Lp      
	  ++Lf;                                    //     incr current leaf     
	}                                          //   END LOOP                
      }                                            // ENDIF                     
      return P->NUMBER ? P : 0;                    // return this box           
    }
    //--------------------------------------------------------------------------
    template<typename cell_iter>
    inline int link_cells_1(                       // R:   tree depth of cell   
			    const box* const&P,    // I:   current box          
			    cell_iter  const&C,    // I:   current cell         
			    cell_iter       &Cf,   // I/O: free cells           
			    typename cell_iter::
			    soul_type* const&S0,   // I:   array with souls     
			    int             &sf)   // I/O: index: free souls    
      const
      // RECURSIVE                                                              
      // This routines transforms the box-leaf tree into the cell-soul tree,    
      // such that all the cells that are contained within some cell are conti- 
      // guous in memory, as are the souls.                                     
    {
      register int dep=0;                          // depth of cell             
      radius_(C) = P->RADIUS;                      // copy radius               
      center_(C) = P->center();                    // copy center               
      number_(C) = P->NUMBER;                      // copy number               
      fcsoul_(C) = sf;                             // set cell: soul kids       
      nsouls_(C) = 0;                              // reset cell: # soul kids   
      register int i,nsub=0;                       // index, # sub-boxes        
      register node*const*N;                       // pointer to current sub    
      for(i=0,N=P->OCT; i!=Nsub; ++i,++N)          // LOOP octants              
	if(P->TYPE & BIT[i]) ++nsub;               //   IF   sub-boxes: count   
	else if(*N) {                              //   ELIF sub-leafs:         
	  PLEAF(*N)->set_soul(S0[sf++]);           //     set soul              
	  nsouls_(C)++;                            //     inc # sub-souls       
	}                                          // END LOOP                  
      if(nsub) {                                   // IF sub-boxes              
	register int       de;                     //   sub-depth               
	register box      *Pi;                     //   sub-box pointer         
	register cell_iter Ci=Cf;                  //   remember free cells     
	fccell_(C) = Ci.index();                   //   set cell: 1st sub-cell  
	ncells_(C) = nsub;                         //   set cell: # sub-cells   
	Cf += nsub;                                //   reserve nsub cells      
	for(i=0,N=P->OCT; i!=Nsub; ++i,++N)        //   LOOP octants            
	  if(*N && P->marked_as_box(i)) {          //     IF oct==sub-box       
	    Pi = PBOX(*N);                         //       sub-box             
	    de = link_cells_1(Pi,Ci++,Cf,S0,sf);   //       recursive call      
	    if(de>dep) dep=de;                     //       update depth        
	  }                                        //   END LOOP                
      } else {                                     // ELSE (no sub-boxes)       
	fccell_(C) =-1;                            //   set cell: 1st sub-cell  
	ncells_(C) = 0;                            //   set cell: # sub-cells   
      }                                            // ENDIF                     
      dep++;                                       // increment depth           
      return dep;                                  // return cells depth        
    }
    //--------------------------------------------------------------------------
    template<typename cell_iter>
    inline int link_cells_N(                       // R:   tree depth of cell   
			    const box* const&P,    // I:   current box          
			    cell_iter  const&C,    // I:   current cell         
			    cell_iter       &Cf,   // I/O: free cells           
			    typename cell_iter::
			    soul_type* const&S0,   // I:   array with souls     
			    int             &sf)   // I/O: index: free souls    
      const
      // RECURSIVE                                                              
      // This routines transforms the box-leaf tree into the cell-soul tree,    
      // such that all the cells that are contained within some cell are conti- 
      // guous in memory, as are the souls.                                     
    {
      register int dep=0;                          // depth of cell             
      radius_(C) = P->RADIUS;                      // copy radius               
      center_(C) = P->center();                    // copy center               
      number_(C) = P->NUMBER;                      // copy number               
      fcsoul_(C) = sf;                             // set cell: soul kids       
      nsouls_(C) = 0;                              // reset cell: # soul kids   
      if(P->is_twig()) {                           // IF box==twig              
	fccell_(C) =-1;                            //   set cell: sub-cells     
	ncells_(C) = 0;                            //   set cell: # sub-cells   
	register leaf*Li=P->LEAFS;                 //   sub-leaf pointer        
	for(; Li; Li=Li->NEXT) {                   //   LOOP sub-leafs          
	  Li->set_soul(S0[sf++]);                  //     set soul              
	  nsouls_(C)++;                            //     increment # sub-souls 
	}                                          //   END LOOP                
      } else {                                     // ELSE (box==branch)        
	register int   i,nsub=0;                   //   index, # sub-boxes      
	register node*const*N;                     //   pointer to current sub  
	for(i=0,N=P->OCT; i!=Nsub; ++i,++N)        //   LOOP octants            
	  if(P->TYPE & BIT[i]) ++nsub;             //     IF is sub-boxe: count 
	  else if(*N) {                            //     ELIF is sub-leafs     
	    PLEAF(*N)->set_soul(S0[sf++]);         //       set soul            
	    nsouls_(C)++;                          //       inc # sub-souls     
	  }                                        //   END LOOP                
	if(nsub) {                                 //   IF has sub-boxes        
	  register int       de;                   //     sub-depth             
	  register box      *Pi;                   //     sub-box pointer       
	  register cell_iter Ci=Cf;                //     remember free cells   
	  fccell_(C) = Ci.index();                 //     set cell: 1st sub-cell
	  ncells_(C) = nsub;                       //     set cell: # sub-cells 
	  Cf += nsub;                              //     reserve nsub cells    
	  for(i=0,N=P->OCT; i!=Nsub; ++i,++N)      //     LOOP octants          
	    if(*N && P->marked_as_box(i)) {        //       IF oct==sub-box     
	      Pi = PBOX(*N);                       //         sub-box           
	      de = link_cells_N(Pi,Ci++,Cf,S0,sf); //         RECURSIVE call    
	      if(de>dep) dep=de;                   //         update depth      
	    }                                      //     END LOOP              
	} else {                                   //   ELSE (no sub-boxes)     
	  fccell_(C) =-1;                          //     set cell: 1st sub-cell
	  ncells_(C) = 0;                          //     set cell: # sub-cells 
	}                                          //   ENDIF                   
      }                                            // ENDIF                     
      dep++;                                       // increment depth           
      return dep;                                  // return cell's depth       
    }
    //--------------------------------------------------------------------------
    box_leaf_tree()
      : BM(0), P0(0) {}
    //--------------------------------------------------------------------------
    void reset(int    const&nc,                    // I: N_crit                 
	       int    const&nu,                    // I: N_cut                  
	       int    const&dm,                    // I: D_max                  
	       size_t const&nl,                    // I: N_leafs                
	       vect   const&x0,                    // I: root center            
	       real   const&sz,                    // I: root radius            
	       size_t const nb = 0)                //[I: #boxes initially alloc]
    {
      NCRIT  = nc;
      NCUT   = nu;
      DMAX   = dm;
      NLEAFS = nl;
      if(BM) delete BM;
      BM     = falcON_Memory(new block_alloc<box>(nb>0? nb : 1+NLEAFS/4));
      P0     = new_box(1);
      P0->center() = x0;
      P0->RADIUS   = sz;
    }
    //--------------------------------------------------------------------------
    box_leaf_tree(int    const&nc,                 // I: N_crit                 
		  int    const&nu,                 // I: N_cut                  
		  int    const&dm,                 // I: D_max                  
		  size_t const&nl,                 // I: N_leafs                
		  vect   const&x0,                 // I: root center            
		  real   const&sz,                 // I: root radius            
		  size_t const&nb = 0) :           //[I: #boxes initially alloc]
      NCRIT  ( nc ),
      NCUT   ( nu ),
      DMAX   ( dm ),
      NLEAFS ( nl ),
      BM     ( falcON_Memory(new block_alloc<box>(nb>0? nb : 1+NLEAFS/4))),
      P0     ( new_box(1) )
    {
      P0->center() = x0;
      P0->RADIUS   = sz;
    }
    //--------------------------------------------------------------------------
    ~box_leaf_tree()
    {
      delete BM;
    }
    //--------------------------------------------------------------------------
    // const public methods (all inlined)                                       
    //--------------------------------------------------------------------------
  public:
    inline size_t       N_allocated() const { return BM->N_allocated(); }
    inline size_t       N_used     () const { return BM->N_used(); }
    inline size_t       N_boxes    () const { return BM->N_used(); }
    inline size_t       N_free     () const { return N_allocated()-N_used(); }
    inline int    const&depth      () const { return DEPTH; }
    inline int    const&maxdepth   () const { return DMAX; }
    inline int    const&Ncrit      () const { return NCRIT; }
    inline size_t const&N_leafs    () const { return NLEAFS; }
    inline box   *const&root       () const { return P0; }
    //--------------------------------------------------------------------------
    // non-const public methods                                                 
    //--------------------------------------------------------------------------
    template<typename cell_iter>
    void link(cell_iter const&C0, typename cell_iter::soul_type*const&S0)
    {
      report REPORT("box_leaf_tree::link()");
      register int       Sf=0;
      register cell_iter Cf=C0+1;
      DEPTH = NCRIT > 1?
	link_cells_N(P0,C0,Cf,S0,Sf) : 
	link_cells_1(P0,C0,Cf,S0,Sf) ;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::tree_builder<tree_type>                                      //
  //                                                                          //
  // for serial tree-building                                                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename tree_type>
  class tree_builder : public box_leaf_tree {
    tree_builder           (const tree_builder&);  // not implemented           
    tree_builder& operator=(const tree_builder&);  // not implemented           
    //--------------------------------------------------------------------------
    // types of class tree_builder                                              
    //--------------------------------------------------------------------------
    typedef typename tree_type::cell_iterator cell_iter; // for cell access     
    typedef typename tree_type::soul_type    *soul_pter; // for soul access     
    typedef typename tree_type::cell_type    *cell_pter; // pointer to cell     
    //--------------------------------------------------------------------------
    // data of class tree_builder                                               
    //--------------------------------------------------------------------------
    vect               XAVE, XMIN, XMAX;           // extreme positions         
    leaf              *L0, *LN;                    // begin/end of leafs        
    const tree_type   *TREE;                       // tree for re-build         
    //--------------------------------------------------------------------------
    // methods of class tree_builder                                            
    //--------------------------------------------------------------------------
#if falcON_NDIM == 2
#define SET_XMIN_XMAX					\
    if     (Li->pos(0) > XMAX[0]) XMAX[0]=Li->pos(0);	\
    else if(Li->pos(0) < XMIN[0]) XMIN[0]=Li->pos(0);	\
    if     (Li->pos(1) > XMAX[1]) XMAX[1]=Li->pos(1);	\
    else if(Li->pos(1) < XMIN[1]) XMIN[1]=Li->pos(1)
#else
#define SET_XMIN_XMAX					\
    if     (Li->pos(0) > XMAX[0]) XMAX[0]=Li->pos(0);	\
    else if(Li->pos(0) < XMIN[0]) XMIN[0]=Li->pos(0);	\
    if     (Li->pos(1) > XMAX[1]) XMAX[1]=Li->pos(1);	\
    else if(Li->pos(1) < XMIN[1]) XMIN[1]=Li->pos(1);	\
    if     (Li->pos(2) > XMAX[2]) XMAX[2]=Li->pos(2);	\
    else if(Li->pos(2) < XMIN[2]) XMIN[2]=Li->pos(2)
#endif
    //--------------------------------------------------------------------------
    inline vect root_center() {
      // This routines returns the root center nearest to the mean position     
      return Integer(XAVE);
    }
    //--------------------------------------------------------------------------
    inline real root_radius(const vect& X) {
      // This routines returns the half-size R of the smallest cube, centered   
      // on X that contains the points xmin and xmax.                           
      register real R,D=zero;                      // R, distance to leaf       
      LoopDims {                                   // loop dimensions          >
	R = max(abs(XMAX[d]-X[d]),
		abs(XMIN[d]-X[d]));                //   distance to xmin, xmax  
	update_max(D,R);                           //   update maximum distance 
      }                                            //                          <
      return pow(two,int(one+log(D)/M_LN2));       // M_LN2 == log(2) (math.h)  
    }
    //--------------------------------------------------------------------------
    inline void build_from_scratch()
      // This routine simply calls addleafs() with the root box as argument.    
    {
      report REPORT("tree_builder::build_from_scratch()");
      register size_t nl=0;                        // counter: # leafs added    
      register leaf  *Li;                          // actual leaf loaded        
      if(Ncrit() > 1)                              // IF(N_crit > 1)            
	for(Li=L0; Li<LN; Li++,nl++)               //   LOOP(leafs)            >
	  addleaf_N(P0,Li,nl,0);                   //     add leafs            <
      else                                         // ELSE                      
	for(Li=L0; Li<LN; Li++,nl++)               //   LOOP(leafs)            >
	  addleaf_1(P0,Li,nl,0);                   //     add leafs            <
    }
    //--------------------------------------------------------------------------
    inline void build_from_tree()
    {
      report REPORT("tree_builder::build_from_tree()");
      // - un_link the tree                                                     
      // - ensure the lost leafs are boxed by root                              
      // - add lost leafs to root                                               
      // IF all leafs were lost from root, we have to abandon our way and build 
      // the tree from scratch.                                                 
      register size_t Na=0;                        // # leafs added             
      register leaf*  Li=L0;                       // leaf: free                
      register leaf_list Lo;                       // list of lost leafs        
      P0 = unlink(TREE->root(),0,Na,Li,Lo);        // un-link tree              
      if(P0) {                                     // IF(root okay)             
	if(Lo.is_empty())                          //   IF(no lost leafs)       
	  return;                                  // <-  DONE                  
	XMIN = XMAX = Lo.HEAD->pos();              //   reset X_min/max         
	LoopLeafListS(Lo,Li) { SET_XMIN_XMAX; }    //   update XMIN & XMAX      
	while(! P0->contains(XMIN))                //   WHILE(Xmin ! in root)  >
	  P0= make_parbox(P0,P0->octant(XMIN),Na); //     enlarge root         <
	while(! P0->contains(XMAX))                //   WHILE(Xmax ! in root)  >
	  P0= make_parbox(P0,P0->octant(XMAX),Na); //     enlarge root         <
	if(Ncrit() > 1) {                          //   IF(N_crit > 1)          
	  BeginLeafList(Lo,Li)                     //     LOOP(lost leafs)     >
	    addleaf_N(P0,Li,Na++,0);               //       add them to root    
	  EndLeafList                              //                          <
	} else {                                   //   ELSE                    
	  BeginLeafList(Lo,Li)                     //     LOOP(lost leafs )    >
	    addleaf_1(P0,Li,Na++,0);               //       add them to root    
	  EndLeafList                              //                          <
	}                                          //   ENDIF                   
      } else {                                     // ELSE (all souls lost)     
	XAVE = XMIN = XMAX = L0->pos();            //   reset X_ave/min/max     
	for(Li=L0+1; Li!=LN; ++Li) {               //   LOOP(leafs)            >
	  SET_XMIN_XMAX;                           //     update X_min/max      
	  XAVE += Li->pos();                       //     sum up X              
	}                                          //                          <
	XAVE/= (LN-L0);                            //   average position        
	register vect X0 = root_center();
	box_leaf_tree::reset(Ncrit(),0,maxdepth(),LN-L0,X0,root_radius(X0));
	build_from_scratch();                      //   and start from scratch  
      }                                            // ENDIF                     
    }
    //==========================================================================
    inline void setup_from_scratch_bodies(const sbodies*const&BB,
					  uint          const&b0 = 0u,
					  uint                nb = 0u,
					  int           const&SP = 0)
    {
      if(nb == 0) nb=BB->N_bodies();               // number of bodies          
      L0 = falcON_New(leaf,nb);                    // allocate leafs            
      register leaf* Li=L0;                        // current leaf              
      XAVE = zero;                                 // reset X_ave               
      XMAX = XMIN = BB->pos(b0);                   // reset X_min/max           
      const uint bn=b0+nb;                         // end bodies                
      for(register uint b=b0; b!=bn; ++b)          // LOOP body flags & pos's   
	if(is_in_tree(BB->flg(b))                  //   IF body is in tree      
	   && SP==0 | is_set(BB->flg(b),SP)) {     //      AND specified        
	  Li->set_up(BB,b);                        //     initialize leaf       
	  SET_XMIN_XMAX;                           //     update XMIN & XMAX    
	  XAVE += Li->pos();                       //     sum up X              
	  Li++;                                    //     incr current leaf     
	}                                          // END LOOP                  
      LN    = Li;                                  // set: beyond last leaf     
      XAVE /= LN-L0;                               // set: X_ave                
    }
    //--------------------------------------------------------------------------
    inline void setup_from_scratch_bodies(const sbodies*const&BB,
					  vect          const&xmin,
					  vect          const&xmax,
					  uint          const&b0 = 0u,
					  uint                nb = 0u,
					  int           const&SP = 0)
    {
      if(nb == 0u) nb=BB->N_bodies();              // number of bodies          
      L0 = falcON_New(leaf,nb);                    // allocate leafs            
      register leaf* Li=L0;                        // current leaf              
      XAVE = zero;                                 // reset X_ave               
      XMIN = xmin;                                 // believe delivered x_min   
      XMAX = xmax;                                 // believe delivered x_max   
      const uint bn=b0+nb;                         // end bodies                
      for(register uint b=b0; b!=bn; ++b)          // LOOP body flags & pos's   
	if(is_in_tree(BB->flg(b))                  //   IF body is in tree      
	   && SP==0 | is_set(BB->flg(b),SP)) {     //      AND specified        
	  Li->set_up(BB,b);                        //     initialize leaf       
	  XAVE += Li->pos();                       //     sum up X              
	  Li++;                                    //     incr current leaf     
	}                                          // END LOOP                  
      LN    = Li;                                  // set: beyond last leaf     
      XAVE /= LN-L0;                               // set: X_ave                
    }
#ifdef falcON_MPI
    //--------------------------------------------------------------------------
    inline void setup_from_scratch_bodies(const pbodies*const&BB,
					  uint          const&b0 = 0u,
					  uint                nb = 0u,
					  int           const&SP = 0)
    {
      if(nb == 0u) nb=BB->N_bodies();              // number of bodies          
      L0 = falcON_New(leaf,nb);                    // allocate leafs            
      register leaf* Li=L0;                        // current leaf              
      XAVE = zero;                                 // reset X_ave               
      XMAX = XMIN = BB->pos(b0);                   // reset x_min & x_max       
      const uint bn=b0+nb;                         // end bodies                
      for(register uint b=b0; b!=bn; ++b)          // LOOP body flags & pos's   
	if(is_in_tree(BB->flg(b))                  //   IF body is in tree      
	   && SP==0 | is_set(BB->flg(b),SP)) {     //      AND specified        
	  Li->set_up(BB,b);                        //     initialize leaf       
	  SET_XMIN_XMAX;                           //     update XMIN & XMAX    
	  XAVE += Li->pos();                       //     sum up X              
	  Li++;                                    //     incr current leaf     
	}                                          // END LOOP                  
      LN    = Li;                                  // set: beyond last leaf     
      XAVE /= LN-L0;                               // set: X_ave                
    }
    //--------------------------------------------------------------------------
    inline void setup_from_scratch_bodies(const pbodies*const&BB,
					  vect          const&xmin,
					  vect          const&xmax,
					  uint          const&b0 = 0u,
					  uint                nb = 0u,
					  int           const&SP = 0)
    {
      if(nb == 0u) nb=BB->N_bodies();              // number of bodies          
      L0 = falcON_New(leaf,nb);                    // allocate leafs            
      register leaf* Li=L0;                        // current leaf              
      XAVE = zero;                                 // reset X_ave               
      XMIN = xmin;                                 // believe delivered x_min   
      XMAX = xmax;                                 // believe delivered x_max   
      const uint bn=b0+nb;                         // end bodies                
      for(register uint b=b0; b!=bn; ++b)          // LOOP body flags & pos's   
	if(is_in_tree(BB->flg(b))                  //   IF body is in tree      
	   && SP==0 | is_set(BB->flg(b),SP)) {     //      AND specified        
	  Li->set_up(BB,b);                        //     initialize leaf       
	  XAVE += Li->pos();                       //     sum up X              
	  Li++;                                    //     incr current leaf     
	}                                          // END LOOP                  
      LN    = Li;                                  // set: beyond last leaf     
      XAVE /= LN-L0;                               // set: X_ave                
    }
#endif
    //--------------------------------------------------------------------------
    inline void setup_from_scratch_arrays(const barrays*const&BA,
					  uint          const&b0 = 0u,
					  uint                nb = 0u,
					  int           const&SP = 0)
    {
      if(nb == 0u) nb=BA->N_bodies();              // number of bodies          
      L0 = falcON_New(leaf,nb);                    // allocate leafs            
      register leaf* Li=L0;                        // current          leaf     
      XAVE = zero;                                 // reset X_ave               
      LoopDims XMAX[d] = XMIN[d] = BA->pos(d)[b0]; // reset x_min & x_max       
      const uint bn=b0+nb;                         // end bodies                
      for(register uint b=b0; b!=bn; ++b)          // LOOP body flags & pos's   
	if(is_in_tree(BA->flg(b))                  //   IF body is in tree      
	   && SP==0 | is_set(BA->flg(b),SP)) {     //      AND specified        
	  Li->set_up(BA,b);                        //     initialize leaf       
	  SET_XMIN_XMAX;                           //     update XMIN & XMAX    
	  XAVE += Li->pos();                       //     sum up X              
	  Li++;                                    //     incr current leaf     
	}                                          // END LOOP                  
      LN    = Li;                                  // set: beyond last leaf     
      XAVE /= LN-L0;                               // set: X_ave                
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    inline void setup_soul_order_bodies(const tree_type  *const&T,
					const bodies_type*const&BB,
					uint              const&b0 = 0u,
					uint                    nb = 0u) {
      if(nb==0u) nb=T->N_souls();                  // number of souls to add    
      L0 = falcON_New(leaf,nb);                    // allocate leafs            
      register leaf*Li = L0;                       // current leaf              
      XAVE = zero;                                 // reset X_ave               
      XMAX = XMIN = BB->pos(T->soul_No(b0)->mybody());  // reset x_min & x_max  
      const uint bn=b0+nb;                         // end souls                 
      LoopSoulsRange(typename tree_type,T,b0,bn,S) {
	Li->set_up(BB,S->mybody());                //   initialize leaf         
	SET_XMIN_XMAX;                             //   update XMIN & XMAX      
	XAVE += Li->pos();                         //     sum up X              
	Li++;                                      //   incr current leaf       
      }                                            // <                         
      LN    = Li;                                  // set: beyond last leaf     
      XAVE /= LN-L0;                               // set: X_ave                
    }
    //--------------------------------------------------------------------------
    inline bool setup_soul_order_arrays(const tree_type*const&T,
					const barrays*  const&BB,
					uint            const&b0 = 0u,
					register uint         nb = 0u) {
      if(nb==0u) nb=T->N_souls();                  // number of souls to add    
      L0 = falcON_New(leaf,nb);                    // allocate leafs            
      register leaf* Li = L0;                      // current leaf              
      XAVE = zero;                                 // reset X_ave               
      LoopDims                                     // reset x_min & x_max       
	XMAX[d] = XMIN[d] = BB->pos(d)[T->soul_No(b0)->mybody()];
      const uint bn=b0+nb;                         // end souls                 
      LoopSoulsRange(typename tree_type,T,b0,bn,S) // LOOP souls of old tree    
	if(is_in_tree(BB->flg(S->mybody()))) {     //   IF(body is in tree)     
	  Li->set_up(BB,S->mybody());              //     initialize leaf       
	  SET_XMIN_XMAX;                           //     update XMIN & XMAX    
	  XAVE += Li->pos();                       //     sum up X              
	  Li++;                                    //     incr current leaf     
	} else {                                   //   ELSE(body not in tree)  
	  delete L0;                               //     delete allocated leafs
	  return false;                            //     return unsuccessful   
	}                                          //   ENDIF                   
      LN    = Li;                                  // set: beyond last leaf     
      XAVE /= LN-L0;                               // set: X_ave                
      return true;                                 // successful                
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    inline void setup_old_tree_bodies(const tree_type  *const&T,
				      const bodies_type*const&BB) {
      LoopSouls(typename tree_type,T,Si)           // loop souls of tree       >
	Si->set_pos(BB);                           //   set soul from body     <
      TREE = T;                                    // tree is still useful      
      L0   = falcON_New(leaf,TREE->N_souls());     // allocate leafs            
      LN   = L0 + TREE->N_souls();                 // set: beyond last leaf     
    }
    //--------------------------------------------------------------------------
    inline bool setup_old_tree_arrays(const tree_type  *const&T,
				      const barrays    *const&BB) {
      LoopSouls(typename tree_type,T,S)            // LOOP souls of old tree    
	if(is_in_tree(BB->flg(S->mybody())))       //   IF(body is active)      
	  S->set_pos(BB);                          //     set new position      
	else {                                     //   ELSE(body not in tree)  
	  TREE = 0;                                //     tree is useless       
	  return false;                            //     return unsuccessful   
	}                                          //   ENDIF                   
      TREE = T;                                    // tree is still useful      
      L0   = falcON_New(leaf,TREE->N_souls());     // allocate leafs            
      LN   = L0 + TREE->N_souls();                 // set: beyond last leaf     
      return true;                                 // successful                
    }
  public:
    //--------------------------------------------------------------------------
    // non-const public methods (almost all non-inline)                         
    //--------------------------------------------------------------------------
    void build();                                  // build box-leaf tree       
    //--------------------------------------------------------------------------
    // constructors of class tree_builder                                       
    //--------------------------------------------------------------------------
    // 1   completely from scratch                                              
    //--------------------------------------------------------------------------
    // 1.1 for the use with bodies with and without input of X_min/max          
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    tree_builder(int               const&,         // I: Ncrit                  
		 int               const&,         // I: Dmax                   
		 const bodies_type*const&,         // I: body sources           
		 int               const& = 0,     //[I: flag specifying bodies]
		 uint              const& = 0u,    //[I: first body]            
		 uint              const& = 0u);   //[I: number of bodies]      
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    tree_builder(int               const&,         // I: Ncrit                  
		 int               const&,         // I: Dmax                   
		 const bodies_type*const&,         // I: body sources           
		 vect              const&,         // I: x_min                  
		 vect              const&,         // I: x_max                  
		 int               const& = 0,     //[I: flag specifying bodies]
		 uint              const& = 0u,    //[I: first body]            
		 uint              const& = 0u);   //[I: number of bodies]      
    //--------------------------------------------------------------------------
    // 1.2 for the use with barrays                                             
    //--------------------------------------------------------------------------
    tree_builder(int               const&,         // I: Ncrit                  
		 int               const&,         // I: Dmax                   
		 const barrays    *const&,         // I: body source arrays     
		 int               const& = 0,     //[I: flag specifying bodies]
		 uint              const& = 0u,    //[I: first body]            
		 uint              const& = 0u);   //[I: number of bodies]      
    //--------------------------------------------------------------------------
    // 2   from scratch, but aided by old tree                                  
    //     we put the leafs to be added in the same order as the souls of the   
    //     old tree. This reduces random memory access, yielding a significant  
    //     speed-up.                                                            
    //                                                                          
    // NOTE  In order to make the code more efficient, we no longer check for   
    //       any potential changes in the tree usage flags (in particular for   
    //       arrays). Thus, if those have changed, don't re-build the tree!     
    //--------------------------------------------------------------------------
    tree_builder(int               const&,         // I: Ncrit                  
		 int               const&,         // I: Dmax                   
		 const tree_type  *const&,         // I: old tree               
		 uint              const& = 0u,    //[I: first soul]            
		 uint              const& = 0u);   //[I: number of bodies/souls]
    //--------------------------------------------------------------------------
    // 3   aided by old tree:                                                   
    //     we build only the lower parts of the old tree anew -- except, of     
    //     course, for bodies that have left the cell. Here, we recognize a     
    //     parameter, N_cut, such that cells with N <= N_cut are re-build from  
    //     scratch. For optimal efficiency (slightly better than the first      
    //     technique above), N_cut should be such that about 90% of bodies      
    //     are still in the cell.                                               
    //                                                                          
    //     if Ncut <= Ncrit, we restore to method 2 above.                      
    //                                                                          
    // NOTE  In order to make the code more efficient, we no longer check for   
    //       any potential changes in the tree usage flags (in particular for   
    //       arrays). Thus, if those have changed, don't re-build the tree!     
    //--------------------------------------------------------------------------
    tree_builder(const tree_type  *const&,         // I: old tree               
		 int               const&,         // I: Ncrit                  
		 int               const&,         // I: Dmax                   
		 int               const&,         // I: Ncut                   
		 uint              const& = 0u,    //[I: first soul]            
		 uint              const& = 0u);   //[I: number of bodies/souls]
    //--------------------------------------------------------------------------
    // destructor                                                               
    //--------------------------------------------------------------------------
    inline ~tree_builder()  {
      delete[] L0;                                 // de-allocate leafs         
    }
    //--------------------------------------------------------------------------
  };
  //============================================================================
#undef LoopDims
#undef PLEAF
#undef PBOX
  //============================================================================
  // non-inline methods of tree_builder                                         
  //============================================================================
  template<typename tree_type>
  void tree_builder<tree_type>::build()
  {
    if(TREE) build_from_tree();
    else     build_from_scratch();
  }
  //----------------------------------------------------------------------------
  template<typename tree_type>
  template<typename bodies_type>
  tree_builder<tree_type>::tree_builder(int               const&nc,
					int               const&dm,
					const bodies_type*const&bb,
					int               const&sp,
					uint              const&b0,
					uint              const&nb)
    : TREE ( 0 )
  {
    report REPORT("tree_builder::tree_builder(): 1");
    setup_from_scratch_bodies(bb,b0,nb,sp);
    register vect X0 = root_center();
    box_leaf_tree::reset(nc,0,dm,LN-L0,X0,root_radius(X0));
  }
  //----------------------------------------------------------------------------
  template<typename tree_type>
  template<typename bodies_type>
  tree_builder<tree_type>::tree_builder(int               const&nc,
					int               const&dm,
					const bodies_type*const&bb,
					vect              const&xmin,
					vect              const&xmax,
					int               const&sp,
					uint              const&b0,
					uint              const&nb)
    : TREE ( 0 )
  {
    report REPORT("tree_builder::tree_builder(): 2");
    setup_from_scratch_bodies(bb,xmin,xmax,b0,nb,sp);
    register vect X0 = root_center();
    box_leaf_tree::reset(nc,0,dm,LN-L0,X0,root_radius(X0));
  }
  //----------------------------------------------------------------------------
  template<typename tree_type>
  tree_builder<tree_type>::tree_builder(int               const&nc,
					int               const&dm,
					const barrays    *const&BA,
					int               const&sp,
					uint              const&b0,
					uint              const&nb)
    : TREE ( 0 )
  {
    report REPORT("tree_builder::tree_builder(): 3");
    setup_from_scratch_arrays(BA,b0,nb,sp);
    register vect X0 = root_center();
    box_leaf_tree::reset(nc,0,dm,LN-L0,X0,root_radius(X0));
  }
  //----------------------------------------------------------------------------
  template<typename tree_type>
  tree_builder<tree_type>::tree_builder(int               const&nc,
					int               const&dm,
					const tree_type  *const&T,
					uint              const&b0,
					uint              const&nb)
    : TREE ( 0 )
  {
    report REPORT("tree_builder::tree_builder(): 4");
    if       (T->use_sbodies()) {                  // CASE 1: use sbodies       
      const sbodies* BB=T->my_sbodies();           //   get bodies              
      if(BB->changes_in_tree_usage_flags())        //   IF(tree-usage changed)  
	setup_from_scratch_bodies(BB,b0,nb);       //     build from scratch    
      else                                         //   ELSE                    
	setup_soul_order_bodies(T,BB,b0,nb);       //     use soul order        
#ifdef falcON_MPI
    } else if(T->use_pbodies()) {                  // CASE 2: use pbodies       
      const pbodies* BB=T->my_pbodies();           //   get bodies              
      if(BB->changes_in_tree_usage_flags())        //   IF(tree-usage changed)  
	setup_from_scratch_bodies(BB,b0,nb);       //     build from scratch    
      else                                         //   ELSE                    
	setup_soul_order_bodies(T,BB,b0,nb);       //     use soul order        
#endif
    } else if(T->use_barrays()) {                  // CASE 3: use barrays       
      const barrays* BA=T->my_barrays();           //   get bodies              
      if(!setup_soul_order_arrays(T,BA,b0,nb))     //   use soul order          
	setup_from_scratch_arrays(BA,b0,nb,T->SP_flag());
    } else 
      error("cannot build from old tree");
    register vect X0 = root_center();
    box_leaf_tree::reset(nc,0,dm,LN-L0,X0,root_radius(X0));
  }
  //----------------------------------------------------------------------------
  template<typename tree_type>
  tree_builder<tree_type>::tree_builder(const tree_type  *const&T,
					int               const&nc,
					int               const&dm,
					int               const&nu,
					uint              const&b0,
					uint              const&nb) 
    : TREE ( 0 )
  {
    report REPORT("tree_builder::tree_builder(): 5");
    if       (T->use_sbodies()) {                  // CASE 1: use sbodies       
      const sbodies* BB=T->my_sbodies();           //   get bodies              
      if(BB->changes_in_tree_usage_flags())        //   IF(tree-usage changed)  
	setup_from_scratch_bodies(BB,b0,nb,T->SP_flag()); // build from scratch 
      else if(nu <= nc)                            //   ELIF(Ncut<=Ncrit)       
	setup_soul_order_bodies(T,BB,b0,nb);       //     use soul order        
      else if(nb != 0u)
	error("cannot partially rebuild");
      else                                         //   ELSE(Ncut >Ncrit)       
	setup_old_tree_bodies(T,BB);               //     use old tree fully   <
#ifdef falcON_MPI
    } else if(T->use_pbodies()) {                  // CASE 2: use pbodies       
      const pbodies* BB=T->my_pbodies();           //   get bodies              
      if(BB->changes_in_tree_usage_flags())        //   IF(tree-usage changed)  
	setup_from_scratch_bodies(BB,b0,nb,T->SP_flag()); // build from scratch 
      else if(nu <= nc)                            //   ELIF(Ncut<=Ncrit)       
	setup_soul_order_bodies(T,BB,b0,nb);       //     use soul order        
      else if(nb != 0u)
	error("cannot partially rebuild");
      else                                         //   ELSE(Ncut >Ncrit)       
	setup_old_tree_bodies(T,BB);               //     use old tree fully   <
#endif
    } else if(T->use_barrays()) {                  // CASE 3: use arrays        
      const barrays* BA=T->my_barrays();           //   get bodies              
      if(nu <= nc) {                               //   IF   (Ncut<=Ncrit)     >
	if(!setup_soul_order_arrays(T,BA,b0,nb))   //   use soul order          
	  setup_from_scratch_arrays(BA,b0,nb,T->SP_flag());
      } else if(nb != 0u)
	error("cannot partially rebuild");
      else {                                       //   ELSE (Ncut >Ncrit)     >
	if(!setup_old_tree_arrays(T,BA))           //     use old tree fully   <
	  setup_from_scratch_arrays(BA,b0,nb,T->SP_flag());
      }
    } else 
      error("cannot build from old tree");
    register vect X0 = root_center();
    if(TREE)
      box_leaf_tree::reset(nc,nu,dm,LN-L0,X0,root_radius(X0),
			   TREE->N_cells()+10*int(sqrt(TREE->N_cells())));
    else
      box_leaf_tree::reset(nc,nu,dm,LN-L0,X0,root_radius(X0));
  }
  //============================================================================
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::basic_tree<TREE,CELL>                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
template<typename TREE, typename CELL> void basic_tree<TREE,CELL>::
mark_for_subtree(                                  // sub-tree mark my nodes    
		 int const&F,                      // I: flag for subtree       
		 uint &Nsubc,                      // O: # subtree cells        
		 uint &Nsubs) const                // O: # subtree souls        
{
  register int  nsun,nsuc;                         // # subt nodes per cell     
  register uint ns=0u, nc=0u;                      // # subt nodes in total     
  LoopMyCellsUp(Ci) {                              // LOOP tree cells backward  
    unflag_subtree_flags(Ci);                      //   reset subtree flags     
    nsun=0;                                        //   counter: subt subnodes  
    nsuc=0;                                        //   counter: subt subcells  
    LoopSoulKids(typename cell_iterator,Ci,s)      //   LOOP child souls        
      if(is_set(s,F)) {                            //     IF flag F is set      
	flag_for_subtree(s);                       //       flag for subtree    
	++nsun;                                    //       count               
      }                                            //   END LOOP                
    ns += nsun;                                    //   count total subt souls  
    LoopCellKids(typename cell_iterator,Ci,c)      //   LOOP child cells        
      if(in_subtree(c)) {                          //     IF cell is in subt    
	++nsun;                                    //       count subt nodes    
	if(is_subtreecell(c) ||                    //       IF is subt cell     
	   is_subtreezombie(c) )                   //       OR is subt zombie   
	  ++nsuc;                                  //         count subcells    
      }                                            //   END LOOP                
    if(nsun > 0) {                                 //   IF cell is subt node    
      flag_for_subtree(Ci);                        //     flag it as such       
      if(nsun > 1) {                               //     IF tcell subt cell    
	flag_as_subtreecell(Ci);                   //       flag it as such     
	++nc;                                      //       count subt cells    
      } else if(nsuc > 0)                          //     ELSE tcell=zombie     
	flag_as_subtreezombie(Ci);                 //       so flag it as such  
    }                                              //   ENDIF                   
  }                                                // END LOOP                  
  Nsubs = ns;                                      // # subt souls              
  Nsubc = nc;                                      // # subt cells              
}
//------------------------------------------------------------------------------
// construction and helpers                                                     
//------------------------------------------------------------------------------
#ifdef TEST_TIMING
#  define SET_I        C_i  = clock();
#  define SET_C        C_0  = clock();
#  define SET_T(TIME)  TIME = (clock() - C_0)/double(CLOCKS_PER_SEC);
#  define SET_F        Ttot = (clock() - C_i)/double(CLOCKS_PER_SEC);	\
      std::cerr<<" time for tree_builder::tree_builder(): "<<Tini<<"\n"	\
               <<" time for tree_builder::build():        "<<Tgro<<"\n"	\
               <<" time for tree_builder::link():         "<<Tlnk<<"\n"	\
               <<" time for tree::tree()                  "<<Ttot<<std::endl;
#else
#  define SET_I        {}
#  define SET_C        {}
#  define SET_T(TIME)  {}
#  define SET_F        {}
#endif
//------------------------------------------------------------------------------
template<typename TREE, typename CELL> inline void basic_tree<TREE,CELL>::
allocate (uint const&ns, uint const&nc) {
  if(ALLOC) delete[] ALLOC;
  ALLOC = falcON_New(char,   3*sizeof(uint)        // Ns, Nc, Dp                
		          + ns*sizeof(soul_type)   // souls                     
		          + nc*sizeof(cell_type)); // cells                     
  DUINT[0] = Ns = ns;
  DUINT[1] = Nc = nc;
  S0 = static_cast<soul_type*>(static_cast<void*>(DUINT+3));
  SN = S0+Ns;
  C0 = static_cast<cell_type*>(static_cast<void*>(SN));
  CN = C0+Nc;
}
//------------------------------------------------------------------------------
template<typename TREE, typename CELL> inline void basic_tree<TREE,CELL>::
set_depth(uint const&dp) {
  DUINT[2] = dp;
}
//------------------------------------------------------------------------------
// construction from sbodies                                                    
//------------------------------------------------------------------------------
template<typename TREE, typename CELL> basic_tree<TREE,CELL>::
basic_tree(const sbodies*const&bb,                 // I: body sources           
	   int           const&nc,                 //[I: N_crit]                
	   int           const&dm,                 //[I: max tree depth]        
	   int           const&sp) :               //[I: flag specifying bodies]
  BSRCES(bb),  ASRCES(0),
#ifdef falcON_MPI
  PSRCES(0),
#endif
  SPFLAG(sp), S0(0), C0(0), ALLOC(0)
{
  SET_I
  SET_C
  tree_builder<base_tree> TB(nc,dm,bb,sp);         // initialize tree_builder   
  SET_T(Tini)
  if(TB.N_leafs()) {                               // IF(leafs in tree)         
    SET_C
    TB.build();                                    //   build box-leaf tree     
    SET_T(Tgro)
    SET_C
    allocate(TB.N_leafs(),TB.N_boxes());           //   allocate souls & cells  
    TB.link(root(),S0);                            //   box-leaf -> cell-soul   
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(Tlnk)
  } else {                                         // ELSE                      
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0);                                 //   reset souls & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // ENDIF                     
  BSRCES->after_tree_growth();                     // reset change flag         
  SET_F
}  
//------------------------------------------------------------------------------
// construction from sbodies with X_min/max known already                       
//------------------------------------------------------------------------------
template<typename TREE, typename CELL> basic_tree<TREE,CELL>::
basic_tree(const sbodies*const&bb,                 // I: body sources           
	   vect          const&xi,                 // I: x_min                  
	   vect          const&xa,                 // I: x_max                  
	   int           const&nc,                 //[I: N_crit]                
	   int           const&dm,                 //[I: max tree depth]        
	   int           const&sp) :               //[I: flag specifying bodies]
  BSRCES(bb), ASRCES(0),
#ifdef falcON_MPI
  PSRCES(0),
#endif
  SPFLAG(sp), S0(0), C0(0), ALLOC(0)
{
  SET_I
  SET_C
  tree_builder<base_tree> TB (nc,dm,bb,xi,xa,sp);  // initiliaze tree_builder   
  SET_T(Tini)
  if(TB.N_leafs()) {                               // IF(leafs in tree)         
    SET_C
    TB.build();                                    //   build box-leaf tree     
    SET_T(Tgro)
    SET_C
    allocate(TB.N_leafs(),TB.N_boxes());           //   allocate souls & cells  
    TB.link(root(),S0);                            //   box-leaf -> cell-soul   
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(Tlnk)
  } else {                                         // < ELSE                  > 
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0);                                 //   reset souls & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // <                         
  BSRCES->after_tree_growth();                     // reset change flag         
  SET_F
}  
#ifdef falcON_MPI
//------------------------------------------------------------------------------
// construction from pbodies                                                    
//------------------------------------------------------------------------------
template<typename TREE, typename CELL> basic_tree<TREE,CELL>::
basic_tree(const pbodies*const&bb,                 // I: body sources           
	   int           const&nc,                 //[I: N_crit]                
	   int           const&dm,                 //[I: max tree depth]        
	   int           const&sp) :               //[I: flag specifying bodies]
  BSRCES(0), ASRCES(0), PSRCES(bb), SPFLAG(sp), S0(0), C0(0), ALLOC(0)
{
  SET_I
  SET_C
  tree_builder<base_tree> TB (nc,dm,bb,sp);        // initiliaze tree_builder   
  SET_T(Tini)
  if(TB.N_leafs()) {                               // IF(leafs in tree)         
    SET_C
    TB.build();                                    //   build box-leaf tree     
    SET_T(Tgro)
    SET_C
    allocate(TB.N_leafs(),TB.N_boxes());           //   allocate souls & cells  
    TB.link(root(),S0);                            //   box-leaf -> cell-soul   
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(Tlnk)
  } else {                                         // < ELSE                  > 
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0);                                 //   reset souls & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // <                         
  PSRCES->after_tree_growth();                     // reset change flag         
  SET_F
}  
//------------------------------------------------------------------------------
// construction from pbodies with X_min/max known already                       
//------------------------------------------------------------------------------
template<typename TREE, typename CELL> basic_tree<TREE,CELL>::
basic_tree(const pbodies*const&bb,                 // I: body sources           
	   vect          const&xi,                 // I: x_min                  
	   vect          const&xa,                 // I: x_max                  
	   int           const&nc,                 //[I: N_crit]                
	   int           const&dm,                 //[I: max tree depth]        
	   int           const&sp) :               //[I: flag specifying bodies]
  BSRCES(0), ASRCES(0), PSRCES(bb), SPFLAG(sp), S0(0), C0(0), ALLOC(0)
{
  SET_I
  SET_C
  tree_builder<base_tree> TB (nc,dm,bb,xi,xa,sp);  // initiliaze tree_builder   
  SET_T(Tini)
  if(TB.N_leafs()) {                               // IF(leafs in tree)         
    SET_C
    TB.build();                                    //   build box-leaf tree     
    SET_T(Tgro)
    SET_C
    allocate(TB.N_leafs(),TB.N_boxes());           //   allocate souls & cells  
    TB.link(root(),S0);                            //   box-leaf -> cell-soul   
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(Tlnk)
  } else {                                         // < ELSE                  > 
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0);                                 //   reset souls & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // <                         
  PSRCES->after_tree_growth();                     // reset change flag         
  SET_F
}
#endif  
//------------------------------------------------------------------------------
// construction from barrays                                                    
//------------------------------------------------------------------------------
template<typename TREE, typename CELL> basic_tree<TREE,CELL>::
basic_tree(const barrays*const&ba,                 // I: body arrays            
	   int           const&nc,                 //[I: N_crit]                
	   int           const&dm,                 //[I: max tree depth]        
	   int           const&sp) :               //[I: flag specifying bodies]
  BSRCES(0),  ASRCES(ba),
#ifdef falcON_MPI
  PSRCES(0),
#endif
  SPFLAG(sp), S0(0), C0(0), ALLOC(0)
{
  SET_I
  SET_C
  tree_builder<base_tree> TB(nc,dm,ba,sp);         // initiliaze tree_builder   
  SET_T(Tini)
  if(TB.N_leafs()) {                               // IF(leafs in tree)         
    SET_C
    TB.build();                                    //   build box-leaf tree     
    SET_T(Tgro)
    SET_C
    allocate(TB.N_leafs(),TB.N_boxes());           //   allocate souls & cells  
    TB.link(root(),S0);                            //   box-leaf -> cell-soul   
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(Tlnk)
  } else {                                         // < ELSE                  > 
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0);                                 //   reset souls & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // <                         
  SET_F
}
//------------------------------------------------------------------------------
// construction as sub-tree from another tree                                   
//------------------------------------------------------------------------------
template<typename TREE, typename CELL>             // type of sub-tree & -cell  
template<typename PARENT_TREE>                     // type of parent-tree       
basic_tree<TREE, CELL>::
basic_tree(const PARENT_TREE* const&parent,        // I: parent tree            
	   int                const&F     ) :      // I: flag specif'ing subtree
  BSRCES(parent->my_sbodies()),                    // copy parent's sbodies     
  ASRCES(parent->my_barrays()),                    // copy parent's barrays     
#ifdef falcON_MPI
  PSRCES(parent->my_pbodies()),                    // copy parent's pbodies     
#endif
  SPFLAG(parent->SP_flag()),                       // copy body specific flag   
  S0(0), C0(0), ALLOC(0)                           // reset some data           
{
  parent->mark_for_subtree(F,Nc,Ns);               // mark parent tree          
  if(Ns==0 || Nc==0) {                             // IF no nodes marked       >
    warning("empty subtree");                      //   issue warning and       
    allocate(0,0);                                 //   reset souls & cells     
    set_depth(0);                                  //   set tree depth to zero  
  } else {                                         // < ELSE                   >
    allocate(Ns,Nc);                               //   allocate souls & cells  
    register cell_iterator Cr=root(),Cf=Cr+1;      //   root, first empty cell  
    register int sf=0;                             //   index: free cell,soul   
    set_depth(                                     //   set tree depth          
	      sub_tree_builder<base_tree,typename PARENT_TREE::base_tree>::
	      link(parent->root(),Cr,Cf,S0,sf));   //   link sub-tree           
  }                                                // <                         
}
//------------------------------------------------------------------------------
// construction as pruned version of another tree                               
//------------------------------------------------------------------------------
template<typename TREE, typename CELL>             // type of sub-tree & -cell  
template<typename PARENT_TREE>                     // type of parent-tree       
basic_tree<TREE,CELL>::
basic_tree(const PARENT_TREE* const&parent,
	   bool(*prune)(typename PARENT_TREE::cell_iterator const&)) :
  BSRCES(0), ASRCES(0),
#ifdef falcON_MPI
  PSRCES(0),
#endif
  SPFLAG(parent->SP_flag()),
  ALLOC(0)
{
  tree_pruner<base_tree,PARENT_TREE> TP(prune,parent);
  allocate(TP.N_souls(),TP.N_cells());
  set_depth(TP.link(root(),S0));
}
//------------------------------------------------------------------------------
// building using the soul-order of the old tree structure                      
//------------------------------------------------------------------------------
template<typename TREE, typename CELL> void basic_tree<TREE,CELL>::
build(int const&nc,                                //[I: N_crit]                
      int const&dm)                                //[I: max tree depth]        
{
  report REPORT("basic_tree::build(%d,%d)",nc,dm);
  SET_I
  SET_C
  tree_builder<base_tree> TB(nc,dm,this);          // initialize tree_builder   
  if(ALLOC) { delete[] ALLOC; ALLOC=0; }           // de-allocate old memory    
  SET_T(Tini)
  if(TB.N_leafs()) {                               // IF(leafs in tree)         
    SET_C
    TB.build();                                    //   build box-leaf tree     
    SET_T(Tgro)
    SET_C
    allocate(TB.N_leafs(),TB.N_boxes());           //   allocate souls & cells  
    TB.link(root(),S0);                            //   box-leaf -> cell-soul   
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(Tlnk)
  } else {                                         // < ELSE                  > 
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0);                                 //   reset souls & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // <                         
  if     (BSRCES) BSRCES->after_tree_growth();     // reset change flag         
#ifdef falcON_MPI
  else if(PSRCES) PSRCES->after_tree_growth();     // reset change flag         
#endif
  SET_F
}
//------------------------------------------------------------------------------
// rebuilding the old tree structure                                            
//------------------------------------------------------------------------------
template<typename TREE, typename CELL> void basic_tree<TREE,CELL>::
rebuild(int const&nu,                              // I: N_cut                  
	int const&nc,                              //[I: N_crit]                
	int const&dm)                              //[I: max tree depth]        
{
  report REPORT("basic_tree::rebuild(%d,%d,%d)",nu,nc,dm);
  SET_I
  SET_C
  tree_builder<base_tree> TB(this,nc,dm,nu);       // initialize tree_builder   
  SET_T(Tini)
  if(TB.N_leafs()) {                               // IF(leafs in tree)         
    SET_C
    TB.build();                                    //   build box-leaf tree     
    SET_T(Tgro)
    SET_C
    allocate(TB.N_leafs(),TB.N_boxes());           //   allocate souls & cells  
    TB.link(root(),S0);                            //   box-leaf -> cell-soul   
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(Tlnk)
  } else {                                         // < ELSE                  > 
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0);                                 //   reset souls & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // <                         
  if     (BSRCES) BSRCES->after_tree_growth();     // reset change flag         
#ifdef falcON_MPI
  else if(PSRCES) PSRCES->after_tree_growth();     // reset change flag         
#endif
  SET_F
}
//------------------------------------------------------------------------------
template<typename TREE, typename CELL>
inline void basic_tree<TREE,CELL>::dump_nodes(cell_iterator const&Ci,
					      std::ostream* const&oc,
					      std::ostream* const&os) const
{
  if(oc) { cell_type::dump(*oc,Ci,Ci.index()); *oc <<"\n"; }
  LoopCellKids(typename cell_iterator,Ci,c) dump_nodes(c,oc,os);
  if(os)
    LoopSoulKids(typename cell_iterator,Ci,s) {
      soul_type::dump(*os,s,int(s-S0));
      *os <<"\n";
    }
}
//------------------------------------------------------------------------------
template<typename TREE, typename CELL>
void basic_tree<TREE,CELL>::dump_nodes(std::ostream*oc, std::ostream*os) const
{
  if(oc == 0 && os == 0) return;
  if(oc) { cell_type::dump_head(*oc); *oc <<"\n"; }
  if(os) { soul_type::dump_head(*os); *os <<"\n"; }
  dump_nodes(root(),oc,os);
  if(oc) oc->flush();
  if(os) os->flush();
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// specialiazations                                                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
template class basic_tree<grav_tree,grav_cell>;

template class basic_tree<sticky_tree,sticky_cell>;
template basic_tree<sticky_tree,sticky_cell>::
  basic_tree(const grav_tree* const&, int const&);
#ifdef falcON_USE_MPI
template class basic_tree<grap_tree,grap_cell>;
template basic_tree<grap_tree,grap_cell>::
  basic_tree(const grav_tree* const&,
	     bool(*)(grav_tree::cell_iterator const&));
#endif
#ifdef falcON_SPH
template class basic_tree<sph_tree,sph_cell>;
template basic_tree<sph_tree,sph_cell>::
  basic_tree(const grav_tree* const&, int const&);
#endif
////////////////////////////////////////////////////////////////////////////////

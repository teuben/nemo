// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tree.cc                                                                     |
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
// tree-building  is done in three steps:                                      |
// - root construction                                                         |
// - building of a box-dot tree                                                |
// - linking to a cell-leaf tree                                               |
//                                                                             |
// root constructions:                                                         |
// In each dimension, the mininum and maximum position is found and from them  |
// the center and size of an appropriate root box computed.                    |
//                                                                             |
// building of a box-dot tree                                                  |
// We first construct a box-dot tree. The dot-adding algorithm is used, ie.    |
// the dots are added one-by-one to the root box (the alternative would be     |
// the box-dividing algorithm, which we found to be slightly less efficient).  |
// The boxes are allocated in blocks, using block_alloc of public/memo.h.      |
// Boxes with less than Ncrit dots are not divided (ie. we wait until a box    |
// has Ncrit dots before splitting it).                                        |
//                                                                             |
// linking to a cell-leaf tree                                                 |
// The box-dot tree is mapped to a cell-leaf tree, such that all cells that    |
// are contained in a given cell are contiguous in memory. The same holds for  |
// the leafs.                                                                  |
//                                                                             |
// Notes                                                                       |
// There are several reasons that make the two-step process of first building  |
// a box-dot tree and then mapping it to a cell-leaf tree worth our while:     |
// - we can arrange sub-cells and sub-leafs to be contiguous in memory;        |
//   this implies that looping over sub-leafs is very fast (no linked lists    |
//   spawming randomly through memory are used), the immediate child leafs     |
//   as well as all the leaf descendants may be addressed easily.              |
// - we can build the tree with memory-minimal entities (boxes are smaller     |
//   then cells, dots smaller than leafs), saving CPU time;                    |
// - we can allocate EXACTLY the correct number of cells;                      |
//                                                                             |
// Variants                                                                    |
// When an old tree is already existent, we may employ the fact that the order |
// of the new tree may differ only little. There are two ways to exploit that: |
// - We may just add the dots to the new tree in the same order as they are    |
//   in the old tree. This ensures that subsequent dots will fall into the     |
//   same box for the most part, reducing random memory access on the boxes.   |
//   This simple method reduces the costs for tree-building by 50% or more for |
//   large N.                                                                  |
// - We may actually take the sorting of the old tree directly. If we still    |
//   want to have a cell-leaf tree with contiguous sub-nodes (and we do), then |
//   we must still go via a box-dot tree. The resulting code is not            |
//   significantly faster then the much simpler method above and in some cases |
//   actually much slower. It is NOT RECOMMENDED (retained for reference only).|
//                                                                             |
// Naming                                                                      |
// Throughout this file, we use the following names:                           |
// body:     iterator through either bodies or ebodies                         |
// leaf:     body-representative in the cell-leaf tree.                        |
//           Here, we assume that the template parameter LEAF is a class that  |
//           has been derived from class basic_leaf of tree.h .                |
// dot:      a minimal copy of a body/leaf, defined below. A dot is used       |
//           only in this file for tree construction                           |
// cell:     cell of the tree to be built.                                     |
//           However, here, we only assume that the template parameter CELL    |
//           is a class that has been derived from basic_cell of tree.h.       |
// box:      tree cell reduced to the tree-building specific needs, defined    |
//           below; only used in this file for tree construction               |
// node:     either a dot or a box; actually, node is base of dot and box      |
// level:    the level of a cell is its 'distance' from root. Usually, root has|
//           level zero.                                                       |
// depth:    the depth of a cell is equal to the maximum level of any of its   |
//           descendants minus its own level.                                  |
//                                                                             |
//-----------------------------------------------------------------------------+
// #define TEST_TIMING

#include <public/tree.h>
#include <public/memo.h>
#include <body.h>
#ifdef falcON_MPI
#  include <walter/pody.h>
#endif

#ifdef  falcON_PROPER
#  define falcON_track_bug
#  undef  falcON_track_bug
#endif

using namespace nbdy;
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::basic_cell_access                                            //
  //                                                                          //
  // any class derived from this one has write access to the tree-specific    //
  // entries of tree cells, which are otherwise not writable.                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class basic_cell_access {
    //--------------------------------------------------------------------------
    // public types                                                             
    //--------------------------------------------------------------------------
  public:
    typedef       oct_tree   *BT;
    typedef       basic_cell *BC;
    typedef       basic_leaf *BL;
    typedef const oct_tree  *CBT;
    typedef const basic_cell*CBC;
    typedef const basic_leaf*CBL;
    //--------------------------------------------------------------------------
    // protected types and methods                                              
    //--------------------------------------------------------------------------
  protected:
    static indx&level_ (BC const&C) { return C->LEVEL; }
    static indx&octant_(BC const&C) { return C->OCTANT; }
    static indx&nleafs_(BC const&C) { return C->NLEAFS; }
    static indx&ncells_(BC const&C) { return C->NCELLS; }
    static int &number_(BC const&C) { return C->NUMBER; }
    static int &fcleaf_(BC const&C) { return C->FCLEAF; }
    static int &fccell_(BC const&C) { return C->FCCELL; }
    static vect&center_(BC const&C) { return C->CENTER; }
    //--------------------------------------------------------------------------
    static void copy_sub  (BC const&C, CBC const&P) { C->copy_sub(P); }
    static void copy_prune(BC const&C, CBC const&P) { C->copy_prune(P); }
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    static size_t   NoLeaf (CBT const&T, CBL const&L) { return T->NoLeaf(L); }
    static size_t   NoCell (CBT const&T, CBC const&C) { return T->NoCell(C); }
    static BC const&FstCell(CBT const&T)              { return T->FstCell(); }
    static BL const&FstLeaf(CBT const&T)              { return T->FstLeaf(); }
    static BC       EndCell(CBT const&T)              { return T->EndCell(); }
    static BL       EndLeaf(CBT const&T)              { return T->EndLeaf(); }
    static BC       CellNo (CBT const&T, int const&I) { return T->CellNo(I); }
    static BL       LeafNo (CBT const&T, int const&I) { return T->LeafNo(I); }
  };
}                                                  // end: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
namespace {
  using nbdy::uint;
  using namespace nbdy;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // auxiliarty macros                                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#define LoopDims for(register int d=0; d<Ndim; ++d)
  //----------------------------------------------------------------------------
#define __LoopLeafKids(TREE,CELL,KID)			\
    for(register BL KID = LeafNo(TREE,fcleaf(CELL));	\
	KID != LeafNo(TREE,ecleaf(CELL)); ++KID)
#define __LoopAllLeafs(TREE,CELL,KID)			\
    for(register BL KID = LeafNo(TREE,fcleaf(CELL));	\
	KID != LeafNo(TREE,ncleaf(CELL)); ++KID)
#define __LoopCellKids(TREE,CELL,KID)			\
    for(register BC KID = CellNo(TREE,fccell(CELL));	\
	KID != CellNo(TREE,eccell(CELL)); ++KID)
#define __LoopLeafs(TREE,LEAF)				\
    for(register BL LEAF = FstLeaf(TREE);		\
        LEAF != EndLeaf(TREE); ++LEAF)
#define __LoopMyCellsUp(CELL)				\
    for(register BC CELL=EndCells()-1;			\
	CELL != FstCell()-1; --CELL)
#define __LoopLeafsRange(TREE,FIRST,END,LEAF)		\
    for(register BL LEAF=LeafNo(TREE,FIRST);		\
	LEAF != LeafNo(TREE,END); ++LEAF)
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // auxiliarty constants                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  const int BIT[8]        = {1,2,4,8,16,32,64,128};// 2^i, i=1,..2^Ndim         
  const int SUBTREECELL   = 1<<8;                  // cell = cell of subtree    
  const int SUBTREE_FLAGS = SUBTREE | SUBTREECELL;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // auxiliarty functions                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#define CFCF const flag*const&F
  inline void flag_as_subtreecell  (flag*F) { F->add   (SUBTREECELL); }
  inline void unflag_subtree_flags (flag*F) { F->un_set(SUBTREE_FLAGS);}
  inline bool in_subtree           (CFCF)   { return F->is_set(SUBTREE); }
  inline bool is_subtreecell       (CFCF)   { return F->is_set(SUBTREECELL); }
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
  // class nbdy::sub_tree_builder                                             //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class sub_tree_builder : private basic_cell_access {
    //--------------------------------------------------------------------------
    // public static method                                                     
    //--------------------------------------------------------------------------
  public:
    static int link(CBT const&,                    // I:   parent tree          
		    CBC const&,                    // I:   current parent cell  
		    CBT const&,                    // I:   daughter tree        
		    BC  const&,                    // I:   current cell to link 
		    BC       &,                    // I/O: next free cell       
		    BL       &);                   // I/O: next free leaf       
    //--------------------------------------------------------------------------
    static int link(                               // R:   depth of tree        
		    CBT const&PT,                  // I:   parent tree          
		    CBT const&DT)                  // I:   daughter tree        
    {
      register BL Lf=FstLeaf(DT);
      register BC Cf=FstCell(DT)+1;
      return link(PT,FstCell(PT), DT, FstCell(DT), Cf, Lf);
    }
  };
  //----------------------------------------------------------------------------
  int sub_tree_builder::link(                      // R:   depth of tree        
			     CBT const&PT,         // I:   parent tree          
			     CBC const&P,          // I:   current parent cell  
			     CBT const&T,          // I:   daughter tree        
			     BC  const&C,          // I:   current cell to link 
			     BC       &Cf,         // I/O: next free cell       
			     BL       &Lf)         // I/O: next free leaf       
  {
    register int dep=0;                            // depth                     
    copy_sub(C,P);                                 // copy level, octant, center
    nleafs_(C) = 0;                                // reset cell: # leaf kids   
    ncells_(C) = 0;                                // reset cell: # cell kids   
    fcleaf_(C) = NoLeaf(T,Lf);                     // set cell: sub-leafs       
    __LoopLeafKids(PT,P,pl)                        // LOOP(leaf kids of Pcell)  
      if(in_subtree(pl)) {                         //   IF(leaf == subt leaf)   
	(Lf++)->copy(pl);                          //     copy link to body etc 
	nleafs_(C)++;                              //     increment # leaf kids 
      }                                            //   ENDIF                   
    __LoopCellKids(PT,P,pc)                        // LOOP(cell kids of Pcell)  
      if(is_subtreecell(pc))                       //   IF(cell==subt cell)     
	ncells_(C)++;                              //     count # subt cell kids
      else if(in_subtree(pc))                      //   ELIF(cell==subt node)   
	__LoopAllLeafs(PT,pc,pl)                   //     LOOP sub cell's leafs 
	  if(in_subtree(pl)) {                     //       IF(leaf==subt leaf) 
	    (Lf++)->copy(pl);                      //         copy link etc     
	    nleafs_(C)++;                          //         incr # leaf kids  
	  }                                        //   ENDIF                   
    number_(C) = nleafs_(C);                       // # leafs >= # leaf kids    
    if(ncells_(C)) {                               // IF(cell has cell kids)    
      register int de;                             //   depth of sub-cell       
      register BC  Ci=Cf;                          //   remember free cells     
      fccell_(C) = NoCell(T,Ci);                   //   set cell children       
      Cf += ncells_(C);                            //   reserve children cells  
      __LoopCellKids(PT,P,pc)                      //   LOOP(c kids of Pcell)   
	if(is_subtreecell(pc)) {                   //     IF(cell == subt cell) 
	  de =link(PT,pc,T,Ci,Cf,Lf);              //       link sub cells      
	  if(de>dep) dep=de;                       //       update depth        
	  number_(C)+= number_(Ci++);              //       count leaf descends 
	}                                          //     ENDIF                 
    } else                                         // ELSE                      
      fccell_(C) =-1;                              //   set pter to cell kids   
    return dep+1;                                  // return cells depth        
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::tree_pruner                                                  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class tree_pruner : private basic_cell_access {
    tree_pruner(tree_pruner const&);               // not implemented           
    //--------------------------------------------------------------------------
    typedef bool(*pruner)(const basic_cell*const&);
    //--------------------------------------------------------------------------
    const pruner prune;
    CBT          parent;
    CBT          pruned;
    int          nC,nS;
    //--------------------------------------------------------------------------
    void count_nodes(CBC const&C) {                // I: cell of parent treee   
      ++nC;                                        // count C                   
      if(!(*prune)(C)) {                           // IF(C is not pruned)       
	nS += nleafs(C);                           //   count leaf kids         
	__LoopCellKids(parent,C,c)                 //   LOOP(cell kids)         
	  count_nodes(c);                          //     RECURSIVE call        
      }                                            // ENDIF                     
    }
    //--------------------------------------------------------------------------
    int link(                                      // R:   depth                
	     CBC const&,                           // I:   parent cell to link  
	     BC  const&,                           // I:   cell to link with    
	     BC       &,                           // I/O: next free cell       
	     BL       &);                          // I/O: next free leaf       
    //--------------------------------------------------------------------------
  public:
    tree_pruner(pruner const&p,
		CBT    const&t)
      : prune(p), parent(t), pruned(0), nC(0), nS(0)
    {
      count_nodes(FstCell(parent));
    }
    //--------------------------------------------------------------------------
    int  const& N_cells () const { return nC; }
    int  const& N_leafs () const { return nS; }
    //--------------------------------------------------------------------------
    int link(                                      // R:   depth                
	     CBT const&t)                          // I:   tree to link         
    {
      pruned = t;
      register BC  C0=FstCell(pruned), Cf=C0+1;
      register BL  Lf=FstLeaf(pruned);
      return link(FstCell(parent),C0,Cf,Lf);
    }
  };
  //----------------------------------------------------------------------------
  int tree_pruner::link(                           // R:   depth                
			CBC const&P,               // I:   parent cell to link  
			BC  const&C,               // I:   cell to link with    
			BC       &Cf,              // I/O: next free cell       
			BL       &Lf)              // I/O: next free leaf       
  {
    register int dep=0;                            // depth                     
    copy_prune(C,P);                               // copy relevant data        
    if((*prune)(P)) {                              // IF(prune cell P)          
      C->add_flag(PRUNE_FLAGS);                    //   mark C as purely repr.  
      nleafs_(C) = 0;                              //   reset # subleafs        
      ncells_(C) = 0;                              //   reset # subcells        
      fcleaf_(C) =-1;                              //   subleafs meaningless    
      fccell_(C) =-1;                              //   subcells meaningless    
    } else {                                       // ELSE(retain original)     
      nleafs_(C) = nleafs(P);                      //   copy # subleafs         
      ncells_(C) = ncells(P);                      //   copy # subcells         
      fcleaf_(C) = NoLeaf(pruned,Lf);              //   set sub-leafs           
      __LoopLeafKids(parent,P,ps)                  //   LOOP(sub-leafs)         
	(Lf++)->copy_prune(ps);                    //     copy properties       
      if(ncells(P)) {                              //   IF(P has cell kids)     
	register int de,prund=0;                   //     sub-depth             
	register BC  Ci=Cf;                        //     remember free cells   
	fccell_(C) = NoCell(pruned,Ci);            //     set cell: 1st sub-cell
	Cf += ncells(P);                           //     reserve sub-cells     
	__LoopCellKids(parent,P,pc) {              //     LOOP(sub-cells)       
	  de = link(pc,Ci,Cf,Lf);                  //       RECURSIVE link      
	  if(de>dep) dep=de;                       //       update depth        
	  prund |= is_pruned(Ci++);                //       pruning info        
	}                                          //     END LOOP              
	if(prund) C->add_flag(PRUNED);             //     IF(pruned) flag so    
      } else                                       //   ELSE                    
	fccell_(C) = -1;                           //     set cell: 1st sub-cell
    }                                              // ENDIF                     
    return dep+1;                                  // return cell's depth       
  }
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
    vect             &pos()             { return POS; }
    vect        const&pos() const       { return POS; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::dot                                                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  struct dot : public node {
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
    mutable dot*NEXT;
    unsigned    LINK;
    //--------------------------------------------------------------------------
    // methods                                                                  
    //--------------------------------------------------------------------------
    void add_to_list(dot* &LIST, int& COUNTER) {
      NEXT = LIST;
      LIST = this;
      ++COUNTER;
    }
    //--------------------------------------------------------------------------
    void add_to_list(node* &LIST, int& COUNTER) {
      NEXT = static_cast<dot*>(LIST);
      LIST = this;
      ++COUNTER;
    }
    //--------------------------------------------------------------------------
    void  set_up  (const basic_leaf* const&L) {
      pos() = nbdy::pos(L);
      LINK  = nbdy::mybody(L);
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void  set_up  (const bodies_type* const&BB, int const&i) {
      LINK  = i;
      pos() = BB->pos(i);
    }
    //--------------------------------------------------------------------------
    void  set_leaf(basic_leaf* const&L) {
      L->set_link_and_pos(LINK,pos());
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::dot_list                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  struct dot_list {
    dot *HEAD;                                     // head of list              
    int  SIZE;                                     // size of list              
    //--------------------------------------------------------------------------
    dot_list() : HEAD(0), SIZE(0) {}
//     dot_list(dot* const&L, int const&N) : HEAD(L), SIZE(N) {}
  private:
    dot_list(const dot_list&L);                    // not implemented           
    dot_list& operator=(const dot_list&);          // not implemented           
    //--------------------------------------------------------------------------
  public:
    bool is_empty() { return SIZE == 0; }          // R: list empty?            
    //--------------------------------------------------------------------------
    void add_dot(dot* const&L) {                   // I: dot to be added        
      L->NEXT = HEAD;                              // set L' next to our list   
      HEAD    = L;                                 // update head of list       
      ++SIZE;                                      // increment size of list    
    }
    //--------------------------------------------------------------------------
    void append(dot_list const&L) {                // I: list to be appended    
      register dot* T=L.HEAD;                      // take head of list L       
      while(T->NEXT) T=T->NEXT;                    // find tail of list L       
      T->NEXT = HEAD;                              // let it point to our list  
      HEAD    = L.HEAD;                            // update head of our list   
      SIZE   += L.SIZE;                            // update size               
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  // this macro requires you to close the curly bracket or use macro EndDotList 
#define BeginDotList(LIST,NAME)		           /* loop elements of list  */\
  for(register dot*NAME=LIST.HEAD, *NEXT_NAME;     /* pointers: current, next*/\
      NAME;					   /* loop until current=0   */\
      NAME = NEXT_NAME) { 			   /* set current = next     */\
  NEXT_NAME = NAME->NEXT;                          /* get next elemetnt      */
#define EndDotList }                               // close curly brackets      
  // this macro fails if dot.NEXT is manipulated within the loop                
#define LoopDotListS(LIST,NAME)		           /* loop elements of list  */\
  for(register dot* NAME = LIST.HEAD;	           /* current dot            */\
      NAME;					   /* loop until current=0   */\
      NAME = NAME->NEXT)                           // set current = next        
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::box                                                          //
  //                                                                          //
  // basic link structure in a box-dot tree.                                  //
  // a box represents a cube centered on center() with half size ("radius")   //
  // equal to box_dot_tree::RA[LEVEL].                                        //
  // if N <= Ncrit, it only contains sub-dots, which are in a linked list     //
  // pointed to by DOTS.                                                      //
  // if N >  Ncrit, DOTS=0 and the sub-nodes are in the array OCT of octants  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  struct box : public node {
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
    indx  LEVEL;                                   // tree level of box         
    indx  TYPE;                                    // bitfield: 1=cell, 0=dot   
    node *OCT[Nsub];                               // octants                   
    int   NUMBER;                                  // number of dots            
    dot  *DOTS;                                    // linked list of dots       
    //--------------------------------------------------------------------------
    // const methods                                                            
    //--------------------------------------------------------------------------
    bool marked_as_box  (int const&i) const{return TYPE & BIT[i]; }
    bool marked_as_dot  (int const&i) const{return !marked_as_box(i); }
    vect const&center   ()            const{return pos(); }
    //--------------------------------------------------------------------------
    int octant(const vect&x) const {               // I: pos  (must be in box)  
      return ::octant(pos(),x);
    }
    //--------------------------------------------------------------------------
    int octant(const dot*const &D) const {         // I: dot (must be in box)   
      return ::octant(pos(),D->pos());
    }
    //--------------------------------------------------------------------------
    int octant(const box* const &P) const {        // I: box (must be in box)   
      return ::octant(pos(),P->pos());
    }
    //--------------------------------------------------------------------------
//     int octant(const basic_cell* const&C) const {  // I: cell (must be in box)  
//       return ::octant(pos(),nbdy::center(C));
//     }
    //--------------------------------------------------------------------------
    bool is_twig() const {
      return DOTS != 0;
    }
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
    box() {}
    //--------------------------------------------------------------------------
    // copy cell properties to box                                              
    //--------------------------------------------------------------------------
    void set_up(const basic_cell*const&C) {
      pos()  = nbdy::center(C);
      LEVEL  = nbdy::level (C);
      TYPE   = 0;
      NUMBER = 0;
      DOTS   = 0;
      for(int b=0; b!=Nsub; ++b) OCT[b] = 0;
    }
    //--------------------------------------------------------------------------
    void mark_as_box(int const&i) { TYPE |=  BIT[i]; }
    void mark_as_dot(int const&i) { TYPE &= ~BIT[i]; }
    vect&center     ()            { return pos(); }
    //--------------------------------------------------------------------------
    box& reset_octants() {
      for(register node**P=OCT; P!=OCT+Nsub; ++P) *P = 0;
      return *this;
    }
    //--------------------------------------------------------------------------
    box& reset() {
      TYPE   = 0;
      NUMBER = 0;
      DOTS   = 0;
      reset_octants();
      return *this;
    }
    //--------------------------------------------------------------------------
    void adddot_to_list(dot* const&L) {            // add dot L to linked list  
      L->add_to_list(DOTS, NUMBER);                // add to linked list        
    }
    //--------------------------------------------------------------------------
    void adddot_to_octs(dot* const &L) {           // add dot L to octants      
      register int b = octant(L);                  // find appropriate octant   
      OCT[b] = L;                                  // fill into octant          
      mark_as_dot(b);                              // mark octant as dot        
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
    const size_t &Ndots;
    const size_t &Nsofar;
  public:
    estimate_N_alloc(const size_t&a, const size_t&b) : Ndots(a), Nsofar(b) {}
    size_t operator() (const size_t&Nused) const {
      register real x = Nused*(real(Ndots)/real(Nsofar)-one);
      return size_t(x+4*sqrt(x)+16);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
#define PDOT static_cast<dot*>
#define PBOX  static_cast<box* >
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::box_dot_tree                                                 //
  //                                                                          //
  // for building of a box-dot tree by adddot()                               //
  // for linking of the box-dot tree to a cell-leaf tree by link_cells()      //
  // does not itself allocate the dots.                                       //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class box_dot_tree : protected basic_cell_access {
    box_dot_tree           (const box_dot_tree&);  // not implemented           
    box_dot_tree& operator=(const box_dot_tree&);  // not implemented           
    //--------------------------------------------------------------------------
    // data of class box_dot_tree                                               
    //--------------------------------------------------------------------------
    static const int   FAC = 100;                  // factor used in unlink()   
    int                NCRIT;                      // Ncrit                     
    int                NCUT, NCUT_FAC;             // Ncut, Ncut*fac            
    int                DMAX, DEPTH;                // max/actual tree depth     
    size_t             NDOTS;                      // # dots (to be) added      
    block_alloc<box>  *BM;                         // allocator for boxes       
  protected:
    CBT                TREE;                       // tree to link              
    real              *RA;                         // array with radius(level)  
    box               *P0;                         // root of box-dot tree      
#ifdef falcON_track_bug
    BL                 LEND;                       // beyond leaf pter range    
    BC                 CEND;                       // beyond cell pter range    
#endif
    //--------------------------------------------------------------------------
    // protected methods are all inlined                                        
    //--------------------------------------------------------------------------
    // radius of box                                                            
    inline real const& radius(const box* const&B) const {
      return RA[B->LEVEL];
    }
    //--------------------------------------------------------------------------
    // does box contain a given position?                                       
    inline bool contains(const box* const&B, const vect&x) const {
      return ::contains(center(B),radius(B),x);
    }
    //--------------------------------------------------------------------------
    // does box contain a given dot?                                            
    inline bool contains(const box* const&B, const dot* const &D) const {
      return contains(B,D->pos());
    }
    //--------------------------------------------------------------------------
#define TREE_TOO_DEEP							       \
    error("exceeding maxium tree depth of %d\n                 "	       \
	  "(presumably more than Ncrit=%d bodies have a common position)",     \
	  DMAX,NCRIT);
    //--------------------------------------------------------------------------
    // shrink box to its octant i                                               
    inline void shrink_to_octant(box* const&B, int const&i) {
      register int l = ++(B->LEVEL);
      if(l > DMAX) TREE_TOO_DEEP
      register real rad=RA[l];
      if(i&1) center(B)[0] += rad;  else  center(B)[0] -= rad;
      if(i&2) center(B)[1] += rad;  else  center(B)[1] -= rad;
#if falcON_NDIM==3
      if(i&4) center(B)[2] += rad;  else  center(B)[2] -= rad;
#endif
    }
    //--------------------------------------------------------------------------
    // expand box in direction i                                                
    inline void expand_to_octant(box* const&B, int const&i) {
      if(B->LEVEL==0) falcON_Error("cannot expand box with level 0");
      register real rad=RA[B->LEVEL];
      if(i&1) center(B)[0] += rad;  else  center(B)[0] -= rad;
      if(i&2) center(B)[1] += rad;  else  center(B)[1] -= rad;
#if falcON_NDIM==3
      if(i&4) center(B)[2] += rad;  else  center(B)[2] -= rad;
#endif
      --(B->LEVEL);
    }
    //--------------------------------------------------------------------------
    inline box* new_box(size_t const&nl) {
      return &(BM->new_element(estimate_N_alloc(NDOTS,nl))->reset());
    }
    //--------------------------------------------------------------------------
    // provides a new (parent) box, whose i th octant is box B                  
    inline box* make_parbox(                       // R: new box                
			    box*   const&B,        // I: daughter box           
			    int    const&i,        // I: direction of extension 
			    size_t const&nl)       // I: # dots added sofar     
    {
      register box* P = new_box(nl);               // get box off the stock     
      P->LEVEL    = B->LEVEL;                      // copy level of parent      
      P->center() = B->center();                   // copy center of parent     
      expand_to_octant(P,i);                       // expand in direction       
      P->addbox_to_octs(B);                        // add B as sub-box          
      return P;                                    // return new box            
    }
    //--------------------------------------------------------------------------
    // provides a new empty (daughter) box in the i th octant of B              
    inline box* make_subbox(                       // R: new box                
			    const box*const&B,     // I: parent box             
			    int       const&i,     // I: parent box's octant    
			    size_t    const&nl)    // I: # dots added sofar     
    {
      register box* P = new_box(nl);               // get box off the stock     
      P->LEVEL    = B->LEVEL;                      // copy level of parent      
      P->center() = B->center();                   // copy center of parent     
      shrink_to_octant(P,i);                       // shrink to correct octant  
      return P;                                    // return new box            
    }
    //--------------------------------------------------------------------------
    // provides a new (daughter) box in the i th octant of B containing dot L   
    // requires that NCRIT == 1                                                 
    inline box* make_subbox_1(                     // R: new box                
			      const box*const&B,   // I: parent box             
			      int       const&i,   // I: parent box's octant    
			      dot      *const&L,   // I: dot of octant          
			      size_t    const&nl)  // I: # dots added sofar     
    {
      register box* P = make_subbox(B,i,nl);       // make new sub-box          
      P->adddot_to_octs(L);                        // add dot to its octant     
      return P;                                    // return new box            
    }
    //--------------------------------------------------------------------------
    // provides a new (daughter) box in the i th octant of B containing dot L   
    // requires that NCRIT >  1                                                 
    inline box* make_subbox_N(                     // R: new box                
			      const box*const&B,   // I: parent box             
			      int       const&i,   // I: parent box's octant    
			      dot      *const&L,   // I: dot of octant          
			      size_t    const&nl)  // I: # dots added sofar     
    {
      register box* P = make_subbox(B,i,nl);       // make new sub-box          
      P->adddot_to_list(L);                        // add old dot to list       
      return P;                                    // return new box            
    }
    //--------------------------------------------------------------------------
    // provides a new empty box at the place of cell C                          
    inline
    box* make_cellbox(                             // R: new box                
		      const basic_cell*const&C,    // I: associated cell        
		      size_t           const&nl)   // I: # dots added sofar     
    {
      register box* P = new_box(nl);               // get box off the stock     
      P->set_up(C);                                // copy cell properties      
      return P;                                    // return new box            
    }
    //--------------------------------------------------------------------------
    // provides a new box with same properties as box Po                        
//     inline box* make_copybox(                      // R: new box                
// 			     box    const&Po,      // I: box to copy            
// 			     size_t const&nl) {    // I: # dots added sofar     
//       return &((new_box(nl))->operator= (Po));
//     }
    //--------------------------------------------------------------------------
    // This routines splits a box:                                              
    // The dots in the linked list are sorted into octants. Octants with one    
    // dot will just hold that dot, octants with many dots will be boxes with   
    // the dots in the linked list.                                             
    // If all dots happen to be in just one octant, the process is repeated on  
    // the box of this octant.                                                  
    void split_box(box*         P,                 // I: box to be splitted     
		   size_t const&nl)                // I: # dots added so far    
    {
      register int NUM[Nsub];                      // array with number of dots 
      register int b,ne;                           // octant & counters         
      register box*sub=0;                          // current sub-box           
      register dot*Di,*Dn;                         // current & next dot        
      do {                                         // DO until # octants > 1    
	for(b=0; b!=Nsub; ++b) NUM[b] = 0;         //   reset counters          
	for(Di=P->DOTS; Di; Di=Dn) {               //   LOOP linked list        
	  Dn = Di->NEXT;                           //     next dot in list      
	  b  = P->octant(Di);                      //     octant of current dot 
	  Di->add_to_list(P->OCT[b], NUM[b]);      //     add dot to list [b]   
	}                                          //   END LOOP                
	P->DOTS = 0;                               //   reset list of sub-dots  
	for(ne=b=0; b!=Nsub; ++b) if(NUM[b]) {     //   LOOP non-empty octs     
	  ne++;                                    //     count them            
	  if(NUM[b]>1) {                           //     IF many dots          
	    sub = make_subbox(P,b,nl);             //       make sub-box        
	    sub->DOTS = PDOT(P->OCT[b]);           //       assign sub-box's    
	    sub->NUMBER = NUM[b];                  //       dot list & number   
	    P->OCT[b] = sub;                       //       set octant=sub-box  
	    P->mark_as_box(b);                     //       mark octant as box  
	  }                                        //     ENDIF                 
	}                                          //   END LOOP                
	P = sub;                                   //   set current box=sub-box 
      } while(ne==1);                              // WHILE only 1 octant       
    }
    //--------------------------------------------------------------------------
    // This routine makes twig boxes contain at most 1 dot                      
    void adddot_1(                                 // add single dot to box     
		  box   *const&base,               // I: base box to add to     
		  dot   *const&Di,                 // I: dot to add             
		  size_t const&nl)                 // I: # dots added sofar     
    {
      for(register box*P=base;;) {                 // LOOP over boxes           
	register int    b  = P->octant(Di);        //   dot's octant            
	register node **oc = P->OCT+b;             //   pointer to octant       
	P->NUMBER++;                               //   increment number        
	if((*oc)==0) {                             //   IF octant empty         
	  *oc = Di;                                //     assign dot to it      
	  P->mark_as_dot(b);                       //     mark octant as dot    
	  return;                                  // <=  DONE with this dot    
	} else if(P->marked_as_dot(b)) {           //   ELIF octant=dot         
	  register dot*Do = PDOT(*oc);             //     get old dot           
	  P->mark_as_box(b);                       //     mark octant as box    
	  P = make_subbox_1(P,b,Do,nl);            //     create sub-box        
	  *oc = P;                                 //     assign sub-box to oc  
	} else                                     //   ELSE octant=box         
	  P = PBOX(*oc);                           //     set current box       
      }                                            // END LOOP                  
    }
    //--------------------------------------------------------------------------
    // This routine makes twig boxes contain at most NCRIT > 1 dots             
    void adddot_N(                                 // add single dot to box     
		  box   *const&base,               // I: base box to add to     
		  dot   *const&Di,                 // I: dot to add             
		  size_t const&nl)                 // I: # dots added sofar     
    {
      for(register box*P=base;;) {                 // LOOP over boxes           
	if(P->is_twig()) {                         //   IF box == twig          
	  P->adddot_to_list(Di);                   //     add dot to list       
	  if(P->NUMBER>NCRIT) split_box(P,nl);     //     IF(N > NCRIT) split   
	  return;                                  //     DONE with this dot    
	} else {                                   //   ELIF box==branch        
	  register int    b  = P->octant(Di);      //   dot's octant            
	  register node **oc = P->OCT+b;           //   pointer to octant       
	  P->NUMBER++;                             //     increment number      
	  if((*oc)==0) {                           //     IF octant empty       
	    *oc = Di;                              //       assign dot to it    
	    P->mark_as_dot(b);                     //       mark octant as dot  
	    return;                                // <-    DONE with this dot  
	  } else if(P->marked_as_dot(b)) {         //     ELIF octant=dot       
	    register dot*Do = PDOT(*oc);           //       get old dot         
	    P->mark_as_box(b);                     //       mark octant as box  
	    P = make_subbox_N(P,b,Do,nl);          //       create sub-box      
	    *oc = P;                               //       assign sub-box to oc
	  } else                                   //     ELSE octant=box       
	    P = PBOX(*oc);                         //       set current box     
	}                                          //   ENDIF                   
      }                                            // END LOOP                  
    }
    //--------------------------------------------------------------------------
    // This routine makes twig boxes contain at most NCRIT dots                 
    void adddot(                                   // add single dot to box     
		box   *const&base,                 // I: base box to add to     
		dot   *const&Di,                   // I: dot to add             
		size_t const&nl)                   // I: # dots added sofar     
    {
      if(NCRIT > 1) adddot_N(base,Di,nl);
      else          adddot_1(base,Di,nl);
    }
    //--------------------------------------------------------------------------
    // map cell-leaf tree back to box-dot tree:                                 
    // cells with N <= N_cut are grown from scratch of their containing leafs   
    // larger cells accumulate a list of "lost" leafs                           
    // NOTE:     this version ASSUMES that all leafs are still in tree          
    box* unlink(                                   // R:   box; 0 if box empty  
		basic_cell*const&,                 // I:   current cell         
		size_t          &,                 // I/O: # dots added sofar   
		dot       *     &,                 // I/O: free dots            
		dot_list        &);                // I/O: list: 'lost' dots    
    //--------------------------------------------------------------------------
    // RECURSIVE                                                                
    // This routines transforms the box-dot tree into the cell-leaf tree,       
    // such that all the cells that are contained within some cell are conti-   
    // guous in memory, as are the leafs.                                       
    int link_cells_1(                              // R:   tree depth of cell   
		     const box *const&,            // I:   current box          
		     int        const&,            // I:   octant of current box
		     basic_cell*const&,            // I:   current cell         
		     basic_cell*     &,            // I/O: index: free cells    
		     basic_leaf*     &) const;     // I/O: index: free leafs    
    //--------------------------------------------------------------------------
    // RECURSIVE                                                                
    // This routines transforms the box-dot tree into the cell-leaf tree,       
    // such that all the cells that are contained within some cell are conti-   
    // guous in memory, as are the leafs.                                       
    int link_cells_N(                              // R:   tree depth of cell   
		     const box* const&,            // I:   current box          
		     int        const&,            // I:   octant of current box
		     basic_cell*const&,            // I:   current cell         
		     basic_cell*     &,            // I/O: index: free cells    
		     basic_leaf*     &) const;     // I/O: index: free leafs    
    //--------------------------------------------------------------------------
    box_dot_tree()
      : BM(0), TREE(0), RA(0), P0(0) {}
    //--------------------------------------------------------------------------
    // to be called before adding any dots.                                     
    // - allocates boxes                                                        
    // - initializes root cell                                                  
    void reset(CBT    const&t,                     // I: tree to be build       
	       int    const&nc,                    // I: N_crit                 
	       int    const&nu,                    // I: N_cut                  
	       int    const&dm,                    // I: D_max                  
	       size_t const&nl,                    // I: N_dots                 
	       vect   const&x0,                    // I: root center            
	       real   const&sz,                    // I: root radius            
	       size_t const nb = 0,                //[I: #boxes initially alloc]
	       indx   const&l0 = 0)                //[I: level of root box]     
    {
      NCRIT    = nc;
      NCUT     = nu;
      NCUT_FAC = NCUT * FAC;
      DMAX     = dm;
      NDOTS    = nl;
      if(BM) delete BM;
      BM       = falcON_Memory(new block_alloc<box>(nb>0? nb : 1+NDOTS/4));
      TREE     = t;
      if(RA) delete RA;
      RA       = falcON_Memory(new real[DMAX+1]);
      P0       = new_box(1);
      RA[l0]   = sz;
      for(register int l=l0; l!=DMAX; ++l) RA[l+1] = half * RA[l];
      for(register int l=l0; l!=0;    --l) RA[l-1] = two  * RA[l];
      P0->LEVEL    = l0;
      P0->center() = x0;
    }
    //--------------------------------------------------------------------------
#if (0) // "...was defined but never used ...
    box_dot_tree(CBT    const&t,                   // I: tree to be build       
		 int    const&nc,                  // I: N_crit                 
		 int    const&nu,                  // I: N_cut                  
		 int    const&dm,                  // I: D_max                  
		 size_t const&nl,                  // I: N_dots                 
		 vect   const&x0,                  // I: root center            
		 real   const&sz,                  // I: root radius            
		 size_t const&nb = 0,              //[I: #boxes initially alloc]
		 indx   const&l0 = 0) :            //[I: level of root box]     
      NCRIT    ( nc ),
      NCUT     ( nu ),
      NCUT_FAC ( nu * FAC ),
      DMAX     ( dm ),
      NDOTS    ( nl ),
      BM       ( falcON_Memory(new block_alloc<box>(nb>0? nb : 1+NDOTS/4))),
      TREE     ( t ),
      RA       ( falcON_Memory(new real[DMAX+1]) ),
      P0       ( new_box(1) )
    {
      RA[l0]       = sz;
      for(register int l=l0; l!=DMAX; ++l) RA[l+1] = half * RA[l];
      for(register int l=l0; l!=0;    --l) RA[l-1] = two  * RA[l];
      P0->center() = x0;
      P0->LEVEL    = l0;
    }
#endif
    //--------------------------------------------------------------------------
    ~box_dot_tree()
    {
      if(BM) delete   BM;
      if(RA) delete[] RA;
    }
    //--------------------------------------------------------------------------
    // const public methods (all inlined)                                       
    //--------------------------------------------------------------------------
  public:
//     inline size_t       N_allocated() const { return BM->N_allocated(); }
//     inline size_t       N_used     () const { return BM->N_used(); }
    inline size_t       N_boxes    () const { return BM->N_used(); }
//     inline size_t       N_free     () const { return N_allocated()-N_used(); }
    inline int    const&depth      () const { return DEPTH; }
    inline int    const&maxdepth   () const { return DMAX; }
    inline int    const&Ncrit      () const { return NCRIT; }
    inline size_t const&N_dots     () const { return NDOTS; }
    inline int          N_levels   () const { return DMAX - P0->LEVEL; }
//     inline box   *const&root       () const { return P0; }
    inline real   const&root_rad   () const { return RA[P0->LEVEL]; }
    //--------------------------------------------------------------------------
    // non-const public methods                                                 
    //--------------------------------------------------------------------------
    void link()
    {
      report REPORT("box_dot_tree::link()");
#ifdef falcON_track_bug
      LEND = EndLeaf(TREE);
      if(LEND != LeafNo(TREE,N_dots()))
	error("box_dot_tree::link(): leaf number mismatch");
      CEND = EndCell(TREE);
      if(CEND != CellNo(TREE,N_boxes()))
	error("box_dot_tree::link(): cell number mismatch");
#endif
      register BC C0 = FstCell(TREE), Cf=C0+1;
      register BL Lf = FstLeaf(TREE);
      DEPTH = NCRIT > 1?
	link_cells_N(P0,0,C0,Cf,Lf) :
	link_cells_1(P0,0,C0,Cf,Lf) ;
    }
  };
  //============================================================================
  // non-inline methods of class box_dot_tree                                   
  //============================================================================
  box* box_dot_tree::                              // R:   box; 0 if box empty  
  unlink(basic_cell*const&C,                       // I:   current cell         
	 size_t          &Na,                      // I/O: # dots added sofar   
	 dot       *     &Df,                      // I/O: free dots            
	 dot_list        &Dp)                      // I/O: list: 'lost' dots    
  {
    register box* P = make_cellbox(C,Na);          // get new box = cell        
    if(level(C) == 0 || number(C) > NCUT_FAC) {    // IF N > N_cut * fac/2      
      register dot_list Do;                        //   list: dots lost from C  
      __LoopLeafKids(TREE,C,Li) {                  //   LOOP leaf kids of C     
	Df->set_up(Li);                            //     map leaf -> dot       
	Do.add_dot(Df);                            //     add to list Do        
	++Df;                                      //     increment dot pter    
      }                                            //   END LOOP                
      register box* Pi;                            //   pter to child box       
      __LoopCellKids(TREE,C,Ci) {                  //   LOOP cell kids of C     
	Pi = unlink(Ci,Na,Df,Do);                  //     RECURSIVE call        
	if(Pi) P->addbox_to_octs(Pi);              //     add box to octants    
      }                                            //   END LOOP                
      if(level(C) == 0  || (                       //   IF root OR              
	 Do.SIZE   > NCUT        &&                //      enough lost dots     
	 number(C) > FAC*Do.SIZE   )) {            //      but much less than N 
	BeginDotList(Do,Di)                        //     LOOP list of lost L   
	  if(contains(P,Di))                       //       IF dot in box       
	    adddot(P,Di,Na++);                     //         add dot to box    
	  else                                     //       ELSE(still lost)    
	    Dp.add_dot(Di);                        //         add dot to Dp     
	EndDotList                                 //     END LOOP              
      } else if(!Do.is_empty())                    //   ELSE (few lost dots)    
	Dp.append(Do);                             //     append list Do to Dp  
    } else if(number(C) > NCUT) {                  // ELIF N > N_cut            
      __LoopLeafKids(TREE,C,Li) {                  //   LOOP leaf kids of C     
	Df->set_up(Li);                            //     map leaf -> dot       
// 	if(contains(P,Df))                         //     IF(dot in box)        
// 	  adddot(P,Df,Na++);                       //       add dot to box      
// 	else                                       //     ELSE(not in box)      
	  Dp.add_dot(Df);                          //       add to list Dp      
	++Df;                                      //     increment dot pter    
      }                                            //   END LOOP                
      register box* Pi;                            //   pter to child box       
      __LoopCellKids(TREE,C,Ci) {                  //   LOOP cell kids of C     
	Pi = unlink(Ci,Na,Df,Dp);                  //     RECURSIVE call        
	if(Pi) P->addbox_to_octs(Pi);              //     add box to octants    
      }                                            //   END LOOP                
    } else {                                       // ELSE(N <= N_cut)          
      __LoopAllLeafs(TREE,C,Li) {                  //   LOOP all leafs in cell  
	Df->set_up(Li);                            //     initialize dot        
	if(contains(P,Df))                         //     IF(dot in box)        
	  adddot(P,Df,Na++);                       //       add dot to box      
	else                                       //     ELSE(not in box)      
	  Dp.add_dot(Df);                          //       add to list Dp      
	++Df;                                      //     incr current dot      
      }                                            //   END LOOP                
    }                                              // ENDIF                     
    return P->NUMBER ? P : 0;                      // return this box           
  }
  //----------------------------------------------------------------------------
  int box_dot_tree::                               // R:   tree depth of cell   
  link_cells_1(const box *const&P,                 // I:   current box          
	       int        const&op,                // I:   octant of current box
	       basic_cell*const&C,                 // I:   current cell         
	       basic_cell*     &Cf,                // I/O: index: free cells    
	       basic_leaf*     &Lf) const          // I/O: index: free leafs    
  {
    register int dep=0;                            // depth of cell             
    level_ (C) = P->LEVEL;                         // copy level                
    octant_(C) = op;                               // set octant                
    center_(C) = P->center();                      // copy center               
    number_(C) = P->NUMBER;                        // copy number               
    fcleaf_(C) = NoLeaf(TREE,Lf);                  // set cell: leaf kids       
    nleafs_(C) = 0;                                // reset cell: # leaf kids   
    register int i,nsub=0;                         // index, # sub-boxes        
    register node*const*N;                         // pointer to current sub    
    for(i=0,N=P->OCT; i!=Nsub; ++i,++N)            // LOOP octants              
      if(P->TYPE & BIT[i]) ++nsub;                 //   IF   sub-boxes: count   
      else if(*N) {                                //   ELIF sub-dots:          
	PDOT(*N)->set_leaf(Lf++);                  //     set leaf              
#ifdef falcON_track_bug
	if(Lf > LEND) report::info("tree_builder::link_cells_1(): "
				   "exceeding max # free leafs");
#endif
	nleafs_(C)++;                              //     inc # sub-leafs       
      }                                            // END LOOP                  
    if(nsub) {                                     // IF sub-boxes              
      register int de;                             //   sub-depth               
      register box*Pi;                             //   sub-box pointer         
      register BC  Ci=Cf;                          //   remember free cells     
      fccell_(C) = NoCell(TREE,Ci);                //   set cell: 1st sub-cell  
      ncells_(C) = nsub;                           //   set cell: # sub-cells   
      Cf += nsub;                                  //   reserve nsub cells      
#ifdef falcON_track_bug
      if(Cf > CEND) report::info("tree_builder::link_cells_1(): "
				 "exceeding max # free cells");
#endif
      for(i=0,N=P->OCT; i!=Nsub; ++i,++N)          //   LOOP octants            
	if(*N && P->marked_as_box(i)) {            //     IF oct==sub-box       
	  Pi = PBOX(*N);                           //       sub-box             
	  de = link_cells_1(Pi,i,Ci++,Cf,Lf);      //       recursive call      
	  if(de>dep) dep=de;                       //       update depth        
	}                                          //   END LOOP                
    } else {                                       // ELSE (no sub-boxes)       
      fccell_(C) =-1;                              //   set cell: 1st sub-cell  
      ncells_(C) = 0;                              //   set cell: # sub-cells   
    }                                              // ENDIF                     
    dep++;                                         // increment depth           
    return dep;                                    // return cells depth        
  }
  //----------------------------------------------------------------------------
  int box_dot_tree::                               // R:   tree depth of cell   
  link_cells_N(const box* const&P,                 // I:   current box          
	       int        const&op,                // I:   octant of current box
	       basic_cell*const&C,                 // I:   current cell         
	       basic_cell*     &Cf,                // I/O: index: free cells    
	       basic_leaf*     &Lf) const          // I/O: index: free leafs    
  {
    register int dep=0;                            // depth of cell             
    level_ (C) = P->LEVEL;                         // copy level                
    octant_(C) = op;                               // set octant                
    center_(C) = P->center();                      // copy center               
    number_(C) = P->NUMBER;                        // copy number               
    fcleaf_(C) = NoLeaf(TREE,Lf);                  // set cell: leaf kids       
    if(P->is_twig()) {                             // IF box==twig              
      fccell_(C) =-1;                              //   set cell: sub-cells     
      ncells_(C) = 0;                              //   set cell: # cell kids   
      nleafs_(C) = P->NUMBER;                      //   set cell: # leaf kids   
      register dot*Di=P->DOTS;                     //   sub-dot pointer         
      for(; Di; Di=Di->NEXT) {                     //   LOOP sub-dots           
	Di->set_leaf(Lf++);                        //     set leaf              
#ifdef falcON_track_bug
	if(Lf > LEND) report::info("tree_builder::link_cells_N(): "
				   "exceeding max # free leafs in twig");
#endif
      }                                            //   END LOOP                
    } else {                                       // ELSE (box==branch)        
      nleafs_(C) = 0;                              //   reset cell: # leaf kids 
      register int i,nsub=0;                       //   index, # sub-boxes      
      register node*const*N;                       //   pointer to current sub  
      for(i=0,N=P->OCT; i!=Nsub; ++i,++N) if(*N)   //   LOOP non-empty octants  
	if(P->marked_as_box(i)) ++nsub;            //     IF is sub-box: count  
	else {                                     //     ELSE (sub-dot)        
	  PDOT(*N)->set_leaf(Lf++);                //       set leaf            
#ifdef falcON_track_bug
	  if(Lf > LEND) report::info("tree_builder::link_cells_N(): "
				     "exceeding max # free leafs");
#endif
	  nleafs_(C)++;                            //       inc # sub-leafs     
	}                                          //   END LOOP                
      if(nsub) {                                   //   IF has sub-boxes        
	register int de;                           //     sub-depth             
	register box*Pi;                           //     sub-box pointer       
	register BC  Ci=Cf;                        //     remember free cells   
	fccell_(C) = NoCell(TREE,Ci);              //     set cell: 1st sub-cel 
	ncells_(C) = nsub;                         //     set cell: # cell kids 
	Cf += nsub;                                //     reserve nsub cells    
#ifdef falcON_track_bug
	if(Cf > CEND) report::info("tree_builder::link_cells_N(): "
				   "exceeding max # free cells");
#endif
	for(i=0,N=P->OCT; i!=Nsub; ++i,++N)        //     LOOP octants          
	  if(*N && P->marked_as_box(i)) {          //       IF oct==sub-box     
	    Pi = PBOX(*N);                         //         sub-box           
	    de = link_cells_N(Pi,i,Ci++,Cf,Lf);    //         RECURSIVE call    
	    if(de>dep) dep=de;                     //         update depth      
	  }                                        //     END LOOP              
      } else {                                     //   ELSE (no sub-boxes)     
	fccell_(C) =-1;                            //     set cell: 1st sub-cell
	ncells_(C) = 0;                            //     set cell: # sub-cells 
      }                                            //   ENDIF                   
    }                                              // ENDIF                     
    dep++;                                         // increment depth           
    return dep;                                    // return cell's depth       
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::tree_builder                                                 //
  //                                                                          //
  // for serial tree-building                                                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class tree_builder : public box_dot_tree {
    tree_builder           (const tree_builder&);  // not implemented           
    tree_builder& operator=(const tree_builder&);  // not implemented           
    //--------------------------------------------------------------------------
    // data of class tree_builder                                               
    //--------------------------------------------------------------------------
    const vect *ROOTCENTER;                        // pre-determined root center
    vect        XAVE, XMIN, XMAX;                  // extreme positions         
    dot        *D0, *DN;                           // begin/end of dots         
    bool        REBUILD;                           // re-build old tree?        
    //--------------------------------------------------------------------------
    // This routines returns the root center nearest to the mean position       
    inline vect root_center() {
      return ROOTCENTER? *ROOTCENTER : Integer(XAVE);
    }
    //--------------------------------------------------------------------------
    // This routines returns the half-size R of the smallest cube, centered     
    // on X that contains the points xmin and xmax.                             
    inline real root_radius(const vect& X) {
      register real R,D=zero;                      // R, distance to dot        
      LoopDims {                                   // LOOP dimensions           
	R=max(abs(XMAX[d]-X[d]),abs(XMIN[d]-X[d]));//   distance to xmin, xmax  
	update_max(D,R);                           //   update maximum distance 
      }                                            // END LOOP                  
      return pow(two,int(one+log(D)/M_LN2));       // M_LN2 == log(2) (math.h)  
    }
    //--------------------------------------------------------------------------
    // This routine simply calls adddots() with the root box as argument.       
    void build_from_scratch();
    //--------------------------------------------------------------------------
    // - un_link the tree                                                       
    // - ensure the lost dots are boxed by root                                 
    // - add lost dots to root                                                  
    // IF all dots were lost from root, we have to abandon our way and build    
    // the tree from scratch.                                                   
    void build_from_tree();
    //--------------------------------------------------------------------------
#define CBB const bodies_type*
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void setup_from_scratch(CBB  const&,
			    uint const& = 0u,
			    uint        = 0u,
			    int  const& = 0);
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void setup_from_scratch(CBB  const&,
			    vect const&,
			    vect const&,
			    uint const& = 0u,
			    uint        = 0u,
			    int  const& = 0);
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void setup_leaf_order(CBB  const&,
			  uint const& = 0u,
			  uint        = 0u);
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void setup_old_tree(CBB const&);
    //--------------------------------------------------------------------------
  public:
    //--------------------------------------------------------------------------
    // non-const public methods (almost all non-inline)                         
    //--------------------------------------------------------------------------
    void build()                                   // build box-dot tree        
    {
      if(REBUILD) build_from_tree();
      else        build_from_scratch();
    }
    //--------------------------------------------------------------------------
    // constructors of class tree_builder                                       
    //--------------------------------------------------------------------------
    // 1   completely from scratch                                              
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    tree_builder(CBT        const&,                // I: tree to be build       
		 const vect*const&,                // I: pre-determined center  
		 int        const&,                // I: Ncrit                  
		 int        const&,                // I: Dmax                   
		 CBB        const&,                // I: body sources           
		 int        const& = 0,            //[I: flag specifying bodies]
		 uint       const& = 0u,           //[I: first body]            
		 uint       const& = 0u);          //[I: number of bodies]      
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    tree_builder(CBT        const&,                // I: tree to be build       
		 const vect*const&,                // I: pre-determined center  
		 int        const&,                // I: Ncrit                  
		 int        const&,                // I: Dmax                   
		 CBB        const&,                // I: body sources           
		 vect       const&,                // I: x_min                  
		 vect       const&,                // I: x_max                  
		 int        const& = 0,            //[I: flag specifying bodies]
		 uint       const& = 0u,           //[I: first body]            
		 uint       const& = 0u);          //[I: number of bodies]      
    //--------------------------------------------------------------------------
    // 2   from scratch, but aided by old tree                                  
    //     we put the dots to be added in the same order as the leafs of the    
    //     old tree. This reduces random memory access, yielding a significant  
    //     speed-up.                                                            
    //                                                                          
    // NOTE  In order to make the code more efficient, we no longer check for   
    //       any potential changes in the tree usage flags (in particular for   
    //       arrays). Thus, if those have changed, don't re-build the tree!     
    //--------------------------------------------------------------------------
    tree_builder(CBT        const&,                // I: old/new tree           
		 const vect*const&,                // I: pre-determined center  
		 int        const&,                // I: Ncrit                  
		 int        const&,                // I: Dmax                   
		 uint       const& = 0u,           //[I: first leaf]            
		 uint       const& = 0u);          //[I: number of bodies/leafs]
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
    tree_builder(CBT  const&,                      // I: old/new tree           
		 int  const&,                      // I: Ncrit                  
		 int  const&,                      // I: Dmax                   
		 int  const&,                      // I: Ncut                   
		 uint const& = 0u,                 //[I: first leaf]            
		 uint const& = 0u);                //[I: number of bodies/leafs]
    //--------------------------------------------------------------------------
    // destructor                                                               
    //--------------------------------------------------------------------------
    inline ~tree_builder()  {
      delete[] D0;                                 // de-allocate dots          
    }
    //--------------------------------------------------------------------------
  };
  //============================================================================
#undef PDOT
#undef PBOX
  //============================================================================
  // non-inline routines of class tree_builder<>                                
  //============================================================================
  void tree_builder::build_from_scratch()
  {
    report REPORT("tree_builder::build_from_scratch()");
    register size_t nl=0;                          // counter: # dots added     
    register dot   *Di;                            // actual dot loaded         
    if(Ncrit() > 1)                                // IF(N_crit > 1)            
      for(Di=D0; Di!=DN; ++Di,++nl)                //   LOOP(dots)              
	adddot_N(P0,Di,nl);                        //     add dots              
    else                                           // ELSE                      
      for(Di=D0; Di!=DN; ++Di,++nl)                //   LOOP(dots)              
	adddot_1(P0,Di,nl);                        //     add dots              
  }
  //----------------------------------------------------------------------------
  void tree_builder::build_from_tree()
  {
    report REPORT("tree_builder::build_from_tree()");
    register size_t Na=0;                          // # dots added              
    register dot*   Di=D0;                         // dot: free                 
    register dot_list Do;                          // list of lost dots         
    P0 = unlink(FstCell(TREE),Na,Di,Do);           // un-link tree              
    if(P0) {                                       // IF(root okay)             
      if(Do.is_empty())                            //   IF(no lost dots)        
	return;                                    // <=  DONE                  
      XMIN = XMAX = Do.HEAD->pos();                //   reset X_min/max         
      LoopDotListS(Do,Di) {                        //   LOOP dots in list       
	XMAX.up_max(Di->pos());                    //     update X_max          
	XMIN.up_min(Di->pos());                    //     update X_min          
      }                                            //   END LOOP                
      while(! contains(P0,XMIN))                   //   WHILE(Xmin ! in root)   
	P0= make_parbox(P0,P0->octant(XMIN),Na);   //     enlarge root          
      while(! contains(P0,XMAX))                   //   WHILE(Xmax ! in root)   
	P0= make_parbox(P0,P0->octant(XMAX),Na);   //     enlarge root          
      if(Ncrit() > 1) {                            //   IF(N_crit > 1)          
	BeginDotList(Do,Di)                        //     LOOP(lost dots)       
	  adddot_N(P0,Di,Na++);                    //       add them to root    
	EndDotList                                 //     END LOOP              
      } else {                                     //   ELSE                    
	BeginDotList(Do,Di)                        //     LOOP(lost dots )      
	  adddot_1(P0,Di,Na++);                    //       add them to root    
	EndDotList                                 //     END LOOP              
      }                                            //   ENDIF                   
    } else {                                       // ELSE (all leafs lost)     
      XAVE = XMIN = XMAX = D0->pos();              //   reset X_ave/min/max     
      for(Di=D0+1; Di!=DN; ++Di) {                 //   LOOP(dots)              
	XMAX.up_max(Di->pos());                    //     update X_max          
	XMIN.up_min(Di->pos());                    //     update X_min          
	XAVE += Di->pos();                         //     sum up X              
      }                                            //   END LOOP                
      XAVE /= real(DN-D0);                         // set: X_ave                
      register vect X0 = root_center();
      box_dot_tree::reset(TREE,Ncrit(),0,maxdepth(),DN-D0,X0,root_radius(X0));
      build_from_scratch();                        //   and start from scratch  
    }                                              // ENDIF                     
  }
  //----------------------------------------------------------------------------
  template<typename bodies_type> void tree_builder::
  setup_from_scratch(CBB  const&BB,
		     uint const&b0,
		     uint       nb,
		     int  const&SP)
  {
    if(nb == 0) nb=BB->N_bodies();                 // number of bodies          
    D0 = falcON_New(dot,nb);                       // allocate dots             
    register dot* Di=D0;                           // current dot               
    XAVE = zero;                                   // reset X_ave               
    XMAX = XMIN = BB->pos(b0);                     // reset X_min/max           
    const uint bn=b0+nb;                           // end bodies                
    for(register uint b=b0; b!=bn; ++b)            // LOOP body flags & pos's   
      if(is_in_tree(BB->flg(b))                    //   IF body is in tree      
	 && SP==0 | is_set(BB->flg(b),SP)) {       //      AND specified        
	Di->set_up(BB,b);                          //     initialize dot        
	XMAX.up_max(Di->pos());                    //     update X_max          
	XMIN.up_min(Di->pos());                    //     update X_min          
	XAVE += Di->pos();                         //     sum up X              
	Di++;                                      //     incr current dot      
      }                                            // END LOOP                  
    DN    = Di;                                    // set: beyond last dot      
    XAVE /= real(DN-D0);                           // set: X_ave                
  }
  //----------------------------------------------------------------------------
  template<typename bodies_type> void tree_builder::
  setup_from_scratch(CBB  const&BB,
		     vect const&xmin,
		     vect const&xmax,
		     uint const&b0,
		     uint       nb,
		     int  const&SP)
  {
    if(nb == 0u) nb=BB->N_bodies();                // number of bodies          
    D0 = falcON_New(dot,nb);                       // allocate dots             
    register dot* Di=D0;                           // current dot               
    XAVE = zero;                                   // reset X_ave               
    XMIN = xmin;                                   // believe delivered x_min   
    XMAX = xmax;                                   // believe delivered x_max   
    const uint bn=b0+nb;                           // end bodies                
    for(register uint b=b0; b!=bn; ++b)            // LOOP body flags & pos's   
      if(is_in_tree(BB->flg(b))                    //   IF body is in tree      
	 && SP==0 | is_set(BB->flg(b),SP)) {       //      AND specified        
	Di->set_up(BB,b);                          //     initialize dot        
	XAVE += Di->pos();                         //     sum up X              
	Di++;                                      //     incr current dot      
      }                                            // END LOOP                  
    DN    = Di;                                    // set: beyond last dot      
    XAVE /= real(DN-D0);                           // set: X_ave                
  }
  //----------------------------------------------------------------------------
  template<typename bodies_type> void tree_builder::
  setup_leaf_order(CBB  const&BB,
		   uint const&b0,
		   uint       nb)
  {
    if(nb==0u) nb=TREE->N_leafs();                 // number of leafs to add    
    D0 = falcON_New(dot,nb);                       // allocate dots             
    register dot*Di = D0;                          // current dot               
    XAVE =zero;                                    // reset X_ave               
    XMAX =XMIN =BB->pos(mybody(LeafNo(TREE,b0)));  // reset x_min & x_max       
    const uint bn=b0+nb;                           // end leafs                 
    __LoopLeafsRange(TREE,b0,bn,Li) {              // LOOP leaf                 
      Di->set_up(BB,mybody(Li));                   //   initialize dot          
      XMAX.up_max(Di->pos());                      //   update X_max            
      XMIN.up_min(Di->pos());                      //   update X_min            
      XAVE += Di->pos();                           //   sum up X                
      Di++;                                        //   incr current dot        
    }                                              // END LOOP                  
    DN    = Di;                                    // set: beyond last dot      
    XAVE /= real(DN-D0);                           // set: X_ave                
  }
  //----------------------------------------------------------------------------
  template<typename bodies_type> void tree_builder::
  setup_old_tree(CBB const&BB)
  {
    __LoopLeafs(TREE,Li)                           // LOOP leafs of tree        
      Li->copy_from_bodies_pos(BB);                //   copy pos: body -> leaf  
    REBUILD = true;                                // we may rebuild            
    D0      = falcON_New(dot,TREE->N_leafs());     // allocate dots             
    DN      = D0 + TREE->N_leafs();                // set: beyond last dot      
  }
  //----------------------------------------------------------------------------
  // constructor 1.1.1                                                          
  template<typename bodies_type>
  tree_builder::tree_builder(CBT        const&t,
			     const vect*const&x0,
			     int        const&nc,
			     int        const&dm,
			     CBB        const&bb,
			     int        const&sp,
			     uint       const&b0,
			     uint       const&nb) :
    ROOTCENTER(x0), REBUILD(false)
  {
    report REPORT("tree_builder::tree_builder(): 1.1.1");
    setup_from_scratch(bb,b0,nb,sp);
    register vect X0 = root_center();
    box_dot_tree::reset(t,nc,0,dm,DN-D0,X0,root_radius(X0));
  }
  //----------------------------------------------------------------------------
  // constructor 1.1.2                                                          
  template<typename bodies_type>
  tree_builder::tree_builder(CBT        const&t,
			     const vect*const&x0,
			     int        const&nc,
			     int        const&dm,
			     CBB        const&bb,
			     vect       const&xmin,
			     vect       const&xmax,
			     int        const&sp,
			     uint       const&b0,
			     uint       const&nb) :
    ROOTCENTER(x0), REBUILD(false)
  {
    report REPORT("tree_builder::tree_builder(): 1.1.2");
    setup_from_scratch(bb,xmin,xmax,b0,nb,sp);
    register vect X0 = root_center();
    box_dot_tree::reset(t,nc,0,dm,DN-D0,X0,root_radius(X0));
  }
#undef CBB
  //----------------------------------------------------------------------------
  // constructor 2                                                              
  tree_builder::tree_builder(CBT        const&t,
			     const vect*const&x0,
			     int        const&nc,
			     int        const&dm,
			     uint       const&b0,
			     uint       const&nb) : 
    ROOTCENTER(x0), REBUILD (false)
  {
    TREE = t;                                      // set tree                  
    report REPORT("tree_builder::tree_builder(): 2");
    if       (TREE->use_bodies()) {                // CASE 1: use bodies        
      const bodies* BB=TREE->my_bodies();          //   get bodies              
      if(BB->changes_in_tree_usage_flags())        //   IF(tree-usage changed)  
	setup_from_scratch(BB,b0,nb);              //     build from scratch    
      else                                         //   ELSE                    
	setup_leaf_order(BB,b0,nb);                //     use leaf order        
#ifdef falcON_MPI
    } else if(TREE->use_pbodies()) {               // CASE 2: use pbodies       
      const pbodies* BB=TREE->my_pbodies();        //   get bodies              
      if(BB->changes_in_tree_usage_flags())        //   IF(tree-usage changed)  
	setup_from_scratch(BB,b0,nb);              //     build from scratch    
      else                                         //   ELSE                    
	setup_leaf_order(BB,b0,nb);                //     use leaf order        
#endif
    } else if(TREE->use_ebodies()) {               // CASE 3: use abodies       
      const ebodies* BB=TREE->my_ebodies();        //   get bodies              
      if(BB->changes_in_tree_usage_flags())        //   IF(tree-usage changed)  
	setup_from_scratch(BB,b0,nb);              //     build from scratch    
      else                                         //   ELSE                    
	setup_leaf_order(BB,b0,nb);                //     use leaf order        
    } else 
      error("cannot build from old tree");
    register vect X0 = root_center();
    box_dot_tree::reset(t,nc,0,dm,DN-D0,X0,root_radius(X0));
  }
  //----------------------------------------------------------------------------
  // constructor 3                                                              
  tree_builder::tree_builder(CBT  const&t,
			     int  const&nc,
			     int  const&dm,
			     int  const&nu,
			     uint const&b0,
			     uint const&nb) : ROOTCENTER(0), REBUILD (false)
  {
    TREE = t;                                      // set tree                  
    report REPORT("tree_builder::tree_builder(): 3");
    if       (TREE->use_bodies()) {                // CASE 1: use bodies        
      const  bodies* BB=TREE->my_bodies();         //   get bodies              
      if(BB->changes_in_tree_usage_flags())        //   IF(tree-usage changed)  
	setup_from_scratch(BB,b0,nb,TREE->SP_flag());//   build frm scratch
      else if(nu <= nc)                            //   ELIF(Ncut<=Ncrit)       
	setup_leaf_order(BB,b0,nb);                //     use leaf order        
      else if(nb != 0u)                            //   ELIF(not all)           
	error("cannot partially rebuild");         //     fatal error           
      else                                         //   ELSE(Ncut >Ncrit)       
	setup_old_tree(BB);                        //     use old tree fully    
#ifdef falcON_MPI
    } else if(TREE->use_pbodies()) {               // CASE 2: use pbodies       
      const pbodies* BB=TREE->my_pbodies();        //   get bodies              
      if(BB->changes_in_tree_usage_flags())        //   IF(tree-usage changed)  
	setup_from_scratch(BB,b0,nb,TREE->SP_flag());//   build frm scratch
      else if(nu <= nc)                            //   ELIF(Ncut<=Ncrit)       
	setup_leaf_order(BB,b0,nb);                //     use leaf order        
      else if(nb != 0u)                            //   ELIF(not all)           
	error("cannot partially rebuild");         //     fatal error           
      else                                         //   ELSE(Ncut >Ncrit)       
	setup_old_tree(BB);                        //     use old tree fully    
#endif
    } else if(TREE->use_ebodies()) {               // CASE 3: use arrays        
      const ebodies* BB=TREE->my_ebodies();        //   get bodies              
      if(BB->changes_in_tree_usage_flags())        //   IF(tree-usage changed)  
	setup_from_scratch(BB,b0,nb,TREE->SP_flag());//   build frm scratch
      else if(nu <= nc)                            //   ELIF(Ncut<=Ncrit)       
	setup_leaf_order(BB,b0,nb);                //     use leaf order        
      else if(nb != 0u)                            //   ELIF(not all)           
	error("cannot partially rebuild");         //     fatal error           
      else                                         //   ELSE(Ncut >Ncrit)       
	setup_old_tree(BB);                        //     use old tree fully    
    }
    if(REBUILD) {
      box_dot_tree::reset(t,nc,nu,dm,DN-D0,t->root_center(),t->root_radius(),
			  TREE->N_cells()+10*int(sqrt(TREE->N_cells())),4);
    } else {
      register vect X0 = root_center();
      box_dot_tree::reset(t,nc,nu,dm,DN-D0,X0,root_radius(X0));
    }
  }
  //////////////////////////////////////////////////////////////////////////////
}                                                  // END: empty namespace      
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::oct_tree                                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
inline unsigned                                    // R: # subtree leafs        
oct_tree::mark_sub(int           const&F,          // I: subtree flag           
		   int           const&Ncr,        // I: Ncrit                  
		   cell_iterator const&C,          // I: cell                   
		   uint               &nc) const   // O: # subtree cells        
  // RECURSIVE                                                                  
  // - count leafs in the subtree                                               
  // - flag cells with any subtree leafs as 'subtree'                           
  // - flag cells with more than Ncrit subtree leafs as 'subtree cells'         
{
  unflag_subtree_flags(C);                         // reset subtree flags       
  register int ns=0;                               // counter: subtree dots     
  LoopLeafKids(cell_iterator,C,Li)                 // LOOP leaf kids            
    if(is_set(Li,F)) {                             //   IF flag F is set        
      flag_for_subtree(Li);                        //     flag for subtree      
      ++ns;                                        //     count                 
    }                                              // END LOOP                  
  LoopCellKids(cell_iterator,C,Ci)                 // LOOP cell kids            
    ns += mark_sub(F,Ncr,Ci,nc);                   //   RECURSIVE call          
  if(ns) {                                         // IF any subtree leafs      
    flag_for_subtree(C);                           //   mark for subtree        
    if(ns >= Ncr) {                                //   IF >=Ncrit subree leafs 
      flag_as_subtreecell(C);                      //     mark as subtree cell  
      ++nc;                                        //     count subtree cells   
    }                                              //   ENDIF                   
  }                                                // ENDIF                     
  return ns;                                       // return: # subtree dots    
}
//------------------------------------------------------------------------------
void oct_tree::mark_for_subtree(int  const&F,      // I: flag for subtree       
				int  const&Ncr,    // I: Ncrit for subtree      
				uint      &Nsubc,  // O: # subtree cells        
				uint      &Nsubs)  // O: # subtree leafs        
  const 
{
  if(Ncr > 1) {                                    // IF Ncrit > 1              
    Nsubc = 0u;                                    //   reset subt cell counter 
    Nsubs = mark_sub(F,Ncr,root(),Nsubc);          //   set flags: subtree_cell 
  } else {                                         // ELSE( Ncrit == 1)         
    register uint subs=0,subc=0;                   //   counter: # subt nodes   
    LoopCellsUp(cell_iterator,this,Ci) {           //   LOOP cells up           
      unflag_subtree_flags(Ci);                    //     reset subtree flags   
      register int ns=0;                           //     # subt dots in cell   
      LoopLeafKids(cell_iterator,Ci,l)             //     LOOP child leafs      
	if(is_set(l,F)) {                          //       IF flag F is set    
	  flag_for_subtree(l);                     //         flag for subtree  
	  ++ns;                                    //         count             
	}                                          //     END LOOP              
      if(ns) {                                     //     IF any subt dots      
	subs += ns;                                //       count # subt dots   
	subc ++;                                   //       count # subt cells  
	flag_for_subtree(Ci);                      //       mark for subtree    
	flag_as_subtreecell(Ci);                   //       mark as subtree cell
      } else                                       //     ELSE (no subt dots)   
	LoopCellKids(cell_iterator,Ci,c)           //       LOOP child cells    
	  if(in_subtree(c)) {                      //         IF cell is in subt
	    flag_for_subtree(Ci);                  //           mark C: subtree 
	    flag_as_subtreecell(Ci);               //           mark C: subtcell
	    break;                                 //           break this loop 
	  }                                        //       END LOOP            
    }                                              //   END LOOP                
    Nsubc = subc;                                  //   set # subtree cells     
    Nsubs = subs;                                  //   set # subtree leafs     
  }                                                // ENDIF                     
}
//------------------------------------------------------------------------------
// construction and helpers                                                     
//------------------------------------------------------------------------------
inline void oct_tree::allocate (uint const&ns, uint const&nc,
				uint const&dm, real const&r0) {
  const uint need =   4*sizeof(uint)               // Ns, Nc, Dp, Dm            
		    +ns*sizeof(basic_leaf)         // leafs                     
		    +nc*sizeof(basic_cell)         // cells                     
                    +(dm+1)*sizeof(real);          // radii of cells            
  if((need > NALLOC) || (need+need < NALLOC)) {
    if(ALLOC) delete16(ALLOC);
    ALLOC  = falcON_New16(char,need);
//     if(ALLOC) delete[] ALLOC;
//     ALLOC  = falcON_New(char,need);
    NALLOC = need;
  }
  DUINT[0] = Ns = ns;
  DUINT[1] = Nc = nc;
  DUINT[3] = dm;
  LEAFS    = static_cast<basic_leaf*>(
	     static_cast<void*>( DUINT+4 ));       // offset of 16bytes         
  CELLS    = static_cast<basic_cell*>(
             static_cast<void*>( LEAFS+Ns ));
  RA       = static_cast<real*>(
             static_cast<void*>( CELLS+Nc ));
  RA[0]    = r0;
  for(register uint l=0; l!=dm; ++l) RA[l+1] = half * RA[l];
}
//------------------------------------------------------------------------------
inline void oct_tree::set_depth(uint const&dp) {
  DUINT[2] = dp;
}
//------------------------------------------------------------------------------
// construction from bodies                                                     
//------------------------------------------------------------------------------
oct_tree::oct_tree(const bodies*const&bb,          // I: body sources           
		   int          const&nc,          // I: N_crit                 
		   const vect  *const&x0,          // I: pre-determined center  
		   int          const&dm,          // I: max tree depth         
		   int          const&sp) :        // I: flag specifying bodies 
  BSRCES(bb), ESRCES(0),
#ifdef falcON_MPI
  PSRCES(0),
#endif
  SPFLAG(sp), LEAFS(0), CELLS(0), ALLOC(0), NALLOC(0u),
  STATE(fresh), USAGE(un_used)
{
  SET_I
  tree_builder TB(this,x0,nc,dm,bb,sp);            // initialize tree_builder   
  SET_T(" time for tree_builder::tree_builder(): ");
  if(TB.N_dots()) {                                // IF(dots in tree)          
    TB.build();                                    //   build box-dot tree      
    SET_T(" time for tree_builder::build():        ");
    allocate(TB.N_dots(),TB.N_boxes(),             //   allocate leafs & cells  
	     TB.N_levels(),TB.root_rad());         //   & set up table: radii   
    TB.link();                                     //   box-dot -> cell-leaf    
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(" time for tree_builder::link():         ");
  } else {                                         // ELSE                      
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0,0,zero);                          //   reset leafs & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // ENDIF                     
  BSRCES->after_tree_growth();                     // reset change flag         
  RCENTER = center(root());                        // set root center           
}  
//------------------------------------------------------------------------------
// construction from bodies with X_min/max known already                        
//------------------------------------------------------------------------------
oct_tree::oct_tree(const bodies*const&bb,          // I: body sources           
		   vect         const&xi,          // I: x_min                  
		   vect         const&xa,          // I: x_max                  
		   int          const&nc,          // I: N_crit                 
		   const vect  *const&x0,          // I: pre-determined center  
		   int          const&dm,          // I: max tree depth         
		   int          const&sp) :        // I: flag specifying bodies 
  BSRCES(bb), ESRCES(0),
#ifdef falcON_MPI
  PSRCES(0),
#endif
  SPFLAG(sp), LEAFS(0), CELLS(0), ALLOC(0), NALLOC(0u),
  STATE(fresh), USAGE(un_used)
{
  SET_I
  tree_builder TB(this,x0,nc,dm,bb,xi,xa,sp);      // initialize tree_builder   
  SET_T(" time for tree_builder::tree_builder(): ");
  if(TB.N_dots()) {                                // IF(dots in tree)          
    TB.build();                                    //   build box-dot tree      
    SET_T(" time for tree_builder::build():        ");
    allocate(TB.N_dots(),TB.N_boxes(),             //   allocate leafs & cells  
	     TB.N_levels(),TB.root_rad());         //   & set up table: radii   
    TB.link();                                     //   box-dot -> cell-leaf    
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(" time for tree_builder::link():         ");
  } else {                                         // ELSE                      
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0,0,zero);                          //   reset leafs & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // ENDIF                     
  BSRCES->after_tree_growth();                     // reset change flag         
  RCENTER = center(root());                        // set root center           
}
#ifdef falcON_MPI
//------------------------------------------------------------------------------
// construction from pbodies                                                    
//------------------------------------------------------------------------------
oct_tree::oct_tree(const pbodies*const&bb,         // I: body sources           
		   int           const&nc,         // I: N_crit                 
		   const vect   *const&x0,         // I: pre-determined center  
		   int           const&dm,         // I: max tree depth         
		   int           const&sp) :       // I: flag specifying bodies 
  BSRCES(0), ESRCES(0), PSRCES(bb),
  SPFLAG(sp), LEAFS(0), CELLS(0), ALLOC(0), NALLOC(0u),
  STATE(fresh), USAGE(un_used)
{
  SET_I
  tree_builder TB(this,x0,nc,dm,bb,sp);            // initialize tree_builder   
  SET_T("time for tree_builder::tree_builder():");
  if(TB.N_dots()) {                                // IF(dots in tree)          
    TB.build();                                    //   build box-dot tree      
    SET_T(" time for tree_builder::build():        ");
    allocate(TB.N_dots(),TB.N_boxes(),             //   allocate leafs & cells  
	     TB.N_levels(),TB.root_rad());         //   & set up table: radii   
    TB.link();                                     //   box-dot -> cell-leaf    
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(" time for tree_builder::link():         ");
  } else {                                         // ELSE                      
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0,0,zero);                          //   reset leafs & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // ENDIF                     
  PSRCES->after_tree_growth();                     // reset change flag         
  RCENTER = center(root());                        // set root center           
}  
//------------------------------------------------------------------------------
// construction from pbodies with X_min/max known already                       
//------------------------------------------------------------------------------
oct_tree::oct_tree(const pbodies*const&bb,         // I: body sources           
		   vect          const&xi,         // I: x_min                  
		   vect          const&xa,         // I: x_max                  
		   int           const&nc,         // I: N_crit                 
		   const vect   *const&x0,         // I: pre-determined center  
		   int           const&dm,         // I: max tree depth         
		   int           const&sp) :       // I: flag specifying bodies 
  BSRCES(0), ESRCES(0), PSRCES(bb),
  SPFLAG(sp), LEAFS(0), CELLS(0), ALLOC(0), NALLOC(0u),
  STATE(fresh), USAGE(un_used)
{
  SET_I
  tree_builder TB(this,x0,nc,dm,bb,xi,xa,sp);      // initialize tree_builder   
  SET_T(" time for tree_builder::tree_builder(): ");
  if(TB.N_dots()) {                                // IF(dots in tree)          
    TB.build();                                    //   build box-dot tree      
    SET_T(" time for tree_builder::build():        ");
    allocate(TB.N_dots(),TB.N_boxes(),             //   allocate leafs & cells  
	     TB.N_levels(),TB.root_rad());         //   & set up table: radii   
    TB.link();                                     //   box-dot -> cell-leaf    
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(" time for tree_builder::link():         ");
  } else {                                         // ELSE                      
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0,0,zero);                          //   reset leafs & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // ENDIF                     
  PSRCES->after_tree_growth();                     // reset change flag         
  RCENTER = center(root());                        // set root center           
}
#endif
//------------------------------------------------------------------------------
// construction from ebodies                                                    
//------------------------------------------------------------------------------
oct_tree::oct_tree(const ebodies*const&bb,         // I: body arrays            
		   int           const&nc,         //[I: N_crit]                
		   const vect   *const&x0,         // I: pre-determined center  
		   int           const&dm,         //[I: max tree depth]        
		   int           const&sp) :       //[I: flag specifying bodies]
  BSRCES(0), ESRCES(bb),
#ifdef falcON_MPI
  PSRCES(0),
#endif
  SPFLAG(sp), LEAFS(0), CELLS(0), ALLOC(0), NALLOC(0u),
  STATE(fresh), USAGE(un_used)
{
  SET_I
  tree_builder TB(this,x0,nc,dm,bb,sp);            // initiliaze tree_builder   
  SET_T(" time for tree_builder::tree_builder(): ");
  if(TB.N_dots()) {                                // IF(dots in tree)          
    TB.build();                                    //   build box-dot tree      
    SET_T(" time for tree_builder::build():        ");
    allocate(TB.N_dots(),TB.N_boxes(),             //   allocate leafs & cells  
	     TB.N_levels(),TB.root_rad());         //   & set up table: radii   
    TB.link();                                     //   box-dot -> cell-leaf    
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(" time for tree_builder::link():         ");
  } else {                                         // ELSE                      
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0,0,zero);                          //   reset leafs & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // ENDIF                     
  RCENTER = center(root());                        // set root center           
}
//------------------------------------------------------------------------------
// construction from ebodies  with X_min/max known already                      
//------------------------------------------------------------------------------
oct_tree::oct_tree(const ebodies*const&bb,         // I: body arrays            
		   vect          const&xi,         // I: x_min                  
		   vect          const&xa,         // I: x_max                  
		   int           const&nc,         //[I: N_crit]                
		   const vect   *const&x0,         // I: pre-determined center  
		   int           const&dm,         //[I: max tree depth]        
		   int           const&sp) :       //[I: flag specifying bodies]
  BSRCES(0), ESRCES(bb),
#ifdef falcON_MPI
  PSRCES(0),
#endif
  SPFLAG(sp), LEAFS(0), CELLS(0), ALLOC(0), NALLOC(0u),
  STATE(fresh), USAGE(un_used)
{
  SET_I
  tree_builder TB(this,x0,nc,dm,bb,xi,xa,sp);      // initiliaze tree_builder   
  SET_T(" time for tree_builder::tree_builder(): ");
  if(TB.N_dots()) {                                // IF(dots in tree)          
    TB.build();                                    //   build box-dot tree      
    SET_T(" time for tree_builder::build():        ");
    allocate(TB.N_dots(),TB.N_boxes(),             //   allocate leafs & cells  
	     TB.N_levels(),TB.root_rad());         //   & set up table: radii   
    TB.link();                                     //   box-dot -> cell-leaf    
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(" time for tree_builder::link():         ");
  } else {                                         // ELSE                      
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0,0,zero);                          //   reset leafs & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // ENDIF                     
  RCENTER = center(root());                        // set root center           
}
//------------------------------------------------------------------------------
// construction as sub-tree from another tree                                   
//------------------------------------------------------------------------------
oct_tree::oct_tree(const oct_tree*const&par,       // I: parent tree            
		   int            const&F,         // I: flag specif'ing subtree
		   int            const&Ncrit) :   //[I: N_crit]                
  BSRCES(par->my_bodies()),                        // copy parent's  bodies     
  ESRCES(par->my_ebodies()),                       // copy parent's ebodies     
#ifdef falcON_MPI
  PSRCES(par->my_pbodies()),                       // copy parent's pbodies     
#endif
  SPFLAG(par->SP_flag()),                          // copy body specific flag   
  LEAFS(0), CELLS(0), ALLOC(0), NALLOC(0u),        // reset some data           
  STATE ( state( par->STATE | sub_tree) ),         // set state                 
  USAGE ( un_used )                                // set usage                 
{
  par->mark_for_subtree(F,Ncrit,Nc,Ns);            // mark parent tree          
  if(Ns==0 || Nc==0) {                             // IF no nodes marked        
    warning("empty subtree");                      //   issue warning and       
    allocate(0,0,0,zero);                          //   reset leafs & cells     
    set_depth(0);                                  //   set tree depth to zero  
  } else {                                         // ELSE                      
    allocate(Ns,Nc,par->depth(),par->root_rad());  //   allocate leafs & cells  
    set_depth(                                     //   set tree depth          
	      sub_tree_builder::link(par,this));   //   link sub-tree           
  }                                                // ENDIF                     
  RCENTER = center(root());                        // set root center           
}
//------------------------------------------------------------------------------
// construction as pruned version of another tree                               
//------------------------------------------------------------------------------
oct_tree::oct_tree(const oct_tree* const&par,
		   bool(*prune)(const basic_cell*const&)) :
  BSRCES(0), ESRCES(0),
#ifdef falcON_MPI
  PSRCES(0),
#endif
  SPFLAG(par->SP_flag()),
  ALLOC(0), NALLOC(0u), STATE(state(par->STATE | pruned)), USAGE(un_used)
{
  tree_pruner TP(prune,par);
  allocate(TP.N_leafs(),TP.N_cells(),              //   allocate leafs & cells  
	   par->depth(),par->root_rad());
  set_depth(TP.link(this));
  RCENTER = center(root());                        // set root center           
}
//------------------------------------------------------------------------------
// building using the leaf-order of the old tree structure                      
//------------------------------------------------------------------------------
void oct_tree::build(int        const&nc,          //[I: N_crit]                
		     const vect*const&x0,          //[I: pre-determined center] 
		     int        const&dm)          //[I: max tree depth]        
{
  report REPORT("oct_tree::build(%d,%d)",nc,dm);
  SET_I
  tree_builder TB(this,x0,nc,dm);                  // initialize tree_builder   
  SET_T(" time for tree_builder::tree_builder(): ");
  if(TB.N_dots()) {                                // IF(dots in tree)          
    TB.build();                                    //   build box-dot tree      
    SET_T(" time for tree_builder::build():        ");
    allocate(TB.N_dots(),TB.N_boxes(),             //   allocate leafs & cells  
	     TB.N_levels(),TB.root_rad());         //   & set up table: radii   
    TB.link();                                     //   box-dot -> cell-leaf    
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(" time for tree_builder::link():         ");
  } else {                                         // ELSE                      
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0,0,zero);                          //   reset leafs & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // ENDIF                     
  if     (BSRCES) BSRCES->after_tree_growth();     // reset bodies' change flag 
#ifdef falcON_MPI
  else if(PSRCES) PSRCES->after_tree_growth();     // reset bodies' change flag 
#endif
  STATE = state((STATE & origins) | re_grown);     // reset state               
  USAGE = un_used;                                 // reset usage flag          
  RCENTER = center(root());                        // set root center           
}
//------------------------------------------------------------------------------
// rebuilding the old tree structure                                            
//------------------------------------------------------------------------------
#if(0) // code currently broken
void oct_tree::rebuild(int const&nu,               // I: N_cut                  
		       int const&nc,               //[I: N_crit]                
		       int const&dm)               //[I: max tree depth]        
{
  report REPORT("oct_tree::rebuild(%d,%d,%d)",nu,nc,dm);
  SET_I
  tree_builder TB(this,nc,dm,nu);                  // initialize tree_builder   
  SET_T(" time for tree_builder::tree_builder(): ");
  if(TB.N_dots()) {                                // IF(dots in tree)          
    TB.build();                                    //   build box-dot tree      
    SET_T(" time for tree_builder::build():        ")
    allocate(TB.N_dots(),TB.N_boxes(),             //   allocate leafs & cells  
	     TB.N_levels(),TB.root_rad());         //   & set up table: radii   
    TB.link();                                     //   box-dot -> cell-leaf    
    set_depth(TB.depth());                         //   set tree depth          
    SET_T(" time for tree_builder::link():         ");
  } else {                                         // ELSE                      
    warning("nobody in tree");                     //   issue a warning         
    allocate(0,0,0,zero);                          //   reset leafs & cells     
    set_depth(0);                                  //   set tree depth to zero  
  }                                                // ENDIF                     
  if     (BSRCES) BSRCES->after_tree_growth();     // reset change flag         
#ifdef falcON_MPI
  else if(PSRCES) PSRCES->after_tree_growth();     // reset change flag         
#endif
  STATE = state((STATE & origins) | re_grown);     // reset state               
  USAGE = un_used;                                 // reset usage flag          
  RCENTER = center(root());                        // set root center           
}
#endif
//------------------------------------------------------------------------------
// re-using old tree structure                                                  
//------------------------------------------------------------------------------
void oct_tree::reuse() {
  if     (BSRCES)
    for(register leaf_iterator Li=begin_leafs(); Li!=end_leafs(); ++Li)
      Li->copy_from_bodies_pos(BSRCES);
#ifdef falcON_MPI
  else if(PSRCES)
    for(register leaf_iterator Li=begin_leafs(); Li!=end_leafs(); ++Li)
      Li->copy_from_bodies_pos(PSRCES);
#endif
  else if(ESRCES)
    for(register leaf_iterator Li=begin_leafs(); Li!=end_leafs(); ++Li)
      Li->copy_from_bodies_pos(ESRCES);
  else
    error("oct_tree::reuse(): without bodies");
  STATE = state((STATE & origins) | re_used);      // reset state               
  USAGE = un_used;                                 // reset usage flag          
}
////////////////////////////////////////////////////////////////////////////////
oct_tree::~oct_tree()
{
  if(ALLOC) { delete16(ALLOC); }
}
////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------+
//                                                                             |
// stic.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// support for sticky particles                                                |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/stic.h>
#include <public/tree.h>
#include <public/iact.h>

namespace {
  using namespace nbdy; using nbdy::uint;
  //////////////////////////////////////////////////////////////////////////////
#define LoopStSpCellKids(CELLITER,CELL,NAME,STSP)		\
  LoopCellKids(CELLITER,CELL,NAME) if(is_##STSP(NAME))
#define LoopStSpLeafKids(CELLITER,CELL,NAME,STSP)		\
  LoopLeafKids(CELLITER,CELL,NAME) if(is_##STSP(NAME))
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // struct take_sph                                                          //
  // struct take_sticky                                                       //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  struct take_sph {
    static bool take    (stsp_estimator::leaf_iterator const&A) 
    { return is_sph(A); }
    static bool take    (stsp_estimator::cell_iterator const&A)
    { return is_sph(A); }
    static bool take_all(stsp_estimator::cell_iterator const&A)
    { return al_sph(A); }
  };
  //////////////////////////////////////////////////////////////////////////////
  struct take_sticky {
    static bool take    (stsp_estimator::leaf_iterator const&A) 
    { return is_sticky(A); }
    static bool take    (stsp_estimator::cell_iterator const&A)
    { return is_sticky(A); }
    static bool take_all(stsp_estimator::cell_iterator const&A)
    { return al_sticky(A); }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class basic_finder                                                       //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename taker>
  class basic_finder : 
    public taker,
    public basic_iactor<stsp_estimator> {
    //--------------------------------------------------------------------------
    // abstract methods                                                         
    //--------------------------------------------------------------------------
    virtual void check_pair(leaf_iter const&, leaf_iter const&) const = 0;
    //--------------------------------------------------------------------------
    // private and protected methods                                            
    //--------------------------------------------------------------------------
  private:
    void many(bool      const&take_all,
	      bool      const&all_active,
	      leaf_iter const&A,
	      leaf_iter const&B0,
	      leaf_iter const&BN) const
    {
      if(take_all) {
	if(all_active) {
	  for(register leaf_iter B=B0; B!=BN; ++B)
	    check_pair(A,B);
	} else {
	  for(register leaf_iter B=B0; B!=BN; ++B)
	    if(is_active(B)) check_pair(A,B);
	}
      } else {
	if(all_active) {
	  for(register leaf_iter B=B0; B!=BN; ++B)
	    if(take(B)) check_pair(A,B);
	} else {
	  for(register leaf_iter B=B0; B!=BN; ++B)
	    if(take(B) && is_active(B)) check_pair(A,B);
	}
      }
    }
    //--------------------------------------------------------------------------
  protected:
    bool many(cell_iter const&A, leaf_iter const&B) const
    {
      many(take_all(A),
	   al_active(A) || is_active(B), B, A.begin_leafs(), A.end_leaf_desc());
      return true;
    }
    //--------------------------------------------------------------------------
    bool many(cell_iter const&A, cell_iter const&B) const
    {
      if(take_all(A)) {
	if(take_all(B))
	  LoopAllLeafs(cell_iter,B,Bi)
	    many(1,           al_active(A)||is_active(Bi),
		 Bi,A.begin_leafs(),A.end_leaf_desc());
	else
	  LoopAllLeafs(cell_iter,B,Bi) if(take(Bi))
	    many(1,           al_active(A)||is_active(Bi),
		 Bi,A.begin_leafs(),A.end_leaf_desc());
      } else
	LoopAllLeafs(cell_iter,A,Ai) if(take(Ai))
	  many(take_all(B),al_active(B)||is_active(Ai),
	       Ai,B.begin_leafs(),B.end_leaf_desc());
      return true;
    }
    //--------------------------------------------------------------------------
    bool many(cell_iter const&A) const {
      // debugged 12-may-2003                                                   
      // bug detected by Clayton Heller in falcON::make_iaction_list()          
      if(take_all(A)) {
	if(al_active(A)) {
	  LoopLstLeafs(cell_iter,A,Ai)
	    many(1,1, Ai,Ai+1,A.end_leaf_desc());
	} else {
	  LoopLstLeafs(cell_iter,A,Ai)
	    many(1,is_active(Ai), Ai,Ai+1,A.end_leaf_desc());
	}
      } else {
	if(al_active(A)) {
	  LoopLstLeafs(cell_iter,A,Ai) if(take(Ai))
	    many(0,1, Ai,Ai+1,A.end_leaf_desc());
	} else {
	  LoopLstLeafs(cell_iter,A,Ai) if(take(Ai))
	    many(0,is_active(Ai), Ai,Ai+1,A.end_leaf_desc());
	}
      }
      return true;
    }
    //--------------------------------------------------------------------------
    void single(leaf_iter const&A, leaf_iter const&B) const
    {
      if(!is_active(A) && !is_active(B)) return;   // no actives -> no job      
      check_pair(A,B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class basic_lister                                                       //
  //                                                                          //
  // for finding leaf pairs in a sticky_tree;                                 //
  // derived from basic_iactor of iact.h, which satisfies the requirements    //
  // for an INTERACTOR template parameter to class MutualInteractor<>;        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename taker>
  class basic_lister : public basic_finder<taker> {
    //--------------------------------------------------------------------------
    // types of class basic_lister                                              
    //--------------------------------------------------------------------------
  public:
    typedef uint elem_pair[2];                     // pair of indices           
    //--------------------------------------------------------------------------
    // data of class basic_lister                                               
    //--------------------------------------------------------------------------
  private:
    const uint    MAX;                             // maximal size of list      
    elem_pair    *BL;                              // list of interaction pairs 
    mutable uint  N;                               // actual size of list       
    //--------------------------------------------------------------------------
    // other methods                                                            
    //--------------------------------------------------------------------------
  protected:
    void add_pair(stsp_estimator::leaf_iterator const&A, 
		  stsp_estimator::leaf_iterator const&B) const
    {
      if(N<MAX)
	if(mybody(A) < mybody(B)) {
	  BL[N][0] = mybody(A);
	  BL[N][1] = mybody(B);
	} else {
	  BL[N][0] = mybody(B);
	  BL[N][1] = mybody(A);
	}
      N++;
      if(N==MAX) warning("interaction list overflow");
    }
    //--------------------------------------------------------------------------
    // abstract methods                                                         
    //--------------------------------------------------------------------------
    virtual void check_pair(stsp_estimator::leaf_iterator const&, 
			    stsp_estimator::leaf_iterator const&) const = 0;
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    uint const &actual_size_of_list() const { return N; }
    //--------------------------------------------------------------------------
    basic_lister(const uint n, elem_pair*l) : MAX(n),BL(l),N(0){}
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class sticky_finder                                                      //
  //                                                                          //
  // for finding stsp_leaf pairs which satisfy:                               //
  // - both are flagged as sticky                                             //
  // - at least one is flagged as active                                      //
  // - there is a time t in [0,tau] such that:                                //
  //   | (x_i+v_i*t) - (x_j+v+j*t) | < size_i+size_j                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class sticky_finder : public basic_lister<take_sticky> {
    const real    TAU;                             // time period               
    //--------------------------------------------------------------------------
  protected:
    void check_pair(leaf_iter const&A, leaf_iter const&B) const
    {
      register vect R = pos(A)-pos(B);             // distance vector           
      register real Sq= square(size(A)+size(B));   // (combined size)^2         
      if(norm(R) < Sq) return add_pair(A,B);       // IF(overlap at t=0): add   
      if(TAU==zero) return;                        // IF(no time > 0) DONE      
      register vect V = vel(A)-vel(B);             // velocity difference       
      register real RV= R*V;                       // scalar product            
      if(RV > zero) return;                        // IF(diverging orbits) DONE 
      register real t = min(TAU,-RV/norm(V));      // time of min distance      
      if(norm(R+t*V) < Sq) add_pair(A,B);          // IF(overlap) add           
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, leaf_iter const&B) const
      // true  if there cannot possibly be an interaction between any A and B   
      // false if an interaction between any A and B is possible                
    {
      register vect R = pos(A)-pos(B);             // R   = distance vector     
      register real
	Rq = norm(R),                              // R^2 = distance^2          
	x  = size(A)+size(B);                      // x   = combined size       
      if(Rq < x*x) return false;                   // overlap at t=0            
      if(TAU==zero) return true;                   // no times t>0              
      register vect V = vel(A)-vel(B);             // V   = velocity diff       
      register real
	w  = vrad(A),                              // w   = v-size              
	wq = w*w,                                  // w^2 = (v-size)^2          
	RV = R*V,                                  // R*V = scalar product      
	RVq= RV*RV;                                // (R*V)^2                   
      if(RV>zero && RVq>wq*Rq) return true;        // R*V/|R| > w: diverging    
      register real
	Vq = norm(V),                              // V^2 = (v-distance)^2      
	t  = (wq>=Vq)? TAU :                       // v-size too large          
	min(TAU, (w*sqrt((Rq*Vq-RVq)/(Vq-wq))-RV)/Vq );// t of min dist         
      if(norm(R+t*V) < square(x+t*w)) return false;// possible overlap          
      return true;                                 // no overlap -> discard     
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, cell_iter const&B) const
      // true  if there cannot possibly be an interaction between any A and B   
      // false if an interaction between any A and B is possible                
    {
      register vect R = pos(A)-pos(B);             // R   = distance vector     
      register real
	Rq = norm(R),                              // R^2 = distance^2          
	x  = size(A)+size(B);                      // x   = combined size       
      if(Rq < x*x) return false;                   // overlap at t=0            
      if(TAU==zero) return true;                   // no times t>0              
      register vect V = vel(A)-vel(B);             // V   = velocity diff       
      register real
	w  = vrad(A)+vrad(B),                      // w   = combined v-size     
	wq = w*w,                                  // w^2 = (comb v-size)^2     
	RV = R*V,                                  // R*V = scalar product      
	RVq= RV*RV;                                // (R*V)^2                   
      if(RV>zero && RVq>wq*Rq) return true;        // R*V/|R| > w: diverging    
      register real
	Vq = norm(V),                              // V^2 = (v-distance)^2      
	t  = (wq>=Vq)? TAU :                       // v-size too large          
	min(TAU, (w*sqrt((Rq*Vq-RVq)/(Vq-wq))-RV)/Vq );// t of min dist         
      if(norm(R+t*V) < square(x+t*w)) return false;// possible overlap          
      return true;                                 // no overlap -> discard     
    }
    //--------------------------------------------------------------------------
  public:
    sticky_finder(const real t, const uint n, elem_pair*l) :
      basic_lister<take_sticky>(n,l), TAU(t) {}
    //--------------------------------------------------------------------------
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || size(A)+TAU*vrad(A) > size(B)+TAU*vrad(B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class neighbour_finder                                                   //
  //                                                                          //
  // for finding stsp_leaf pairs which satisfy:                               //
  // - both are flagged as sph                                                //
  // - at least one is flagged as active                                      //
  // - | x_i - x_j | < max(size_i,size_j)                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class neighbour_finder : public basic_lister<take_sph> {
  protected:
    void check_pair(leaf_iter const&A, leaf_iter const&B) const
    {
      register real Rq = dist_sq(pos(A),pos(B));
      if(Rq < sizeq(A) || Rq < sizeq(B)) add_pair(A,B);
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, leaf_iter const&B) const
    {
      return dist_sq(pos(A),pos(B)) > square(max(rmax(A)+size(B),size(A)));
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, cell_iter const&B) const
    {
      return dist_sq(pos(A),pos(B)) >
	     square( max( rmax(A)+size(B) , rmax(B)+size(A) ) );
    }
  public:
    //--------------------------------------------------------------------------
    neighbour_finder(const uint n, elem_pair*l) : basic_lister<take_sph>(n,l) {}
    //--------------------------------------------------------------------------
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || size(A) > size(B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class partner_finder                                                     //
  //                                                                          //
  // for finding stsp_leaf pairs which satisfy:                               //
  // - both are flagged as sph                                                //
  // - at least one is flagged as active                                      //
  // - | x_i - x_j | < size_i + size_j                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class partner_finder : public basic_lister<take_sph> {
  protected:
    void check_pair(leaf_iter const&A, leaf_iter const&B) const
    {
      if(dist_sq(pos(A),pos(B)) < square(size(A)+size(B)))
	add_pair(A,B);
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, leaf_iter const&B) const
    {
      return dist_sq(pos(A),pos(B)) > square(size(A)+size(B));
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, cell_iter const&B) const
    {
      return dist_sq(pos(A),pos(B)) > square(size(A)+size(B));
    }
  public:
    //--------------------------------------------------------------------------
    partner_finder(const uint n, elem_pair*l) : basic_lister<take_sph>(n,l) {}
    //--------------------------------------------------------------------------
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || size(A) > size(B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class partner_counter                                                    //
  //                                                                          //
  // for counting stsp_leaf pairs which satisfy:                              //
  // - both are flagged as sph                                                //
  // - | x_i - x_j | < size_i + size_j                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class partner_counter : public basic_finder<take_sph> {
  protected:
    void check_pair(leaf_iter const&A, leaf_iter const&B) const
    {
      if(dist_sq(pos(A),pos(B)) < square(size(A)+size(B))) {
	A->inc();
	B->inc();
      }
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, leaf_iter const&B) const
    {
      return dist_sq(pos(A),pos(B)) > square(size(A)+size(B));
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, cell_iter const&B) const
    {
      return dist_sq(pos(A),pos(B)) > square(size(A)+size(B));
    }
  public:
    //--------------------------------------------------------------------------
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || size(A) > size(B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  typedef stsp_estimator::cell_iterator StSpCellIter;
  //----------------------------------------------------------------------------
  inline stsp_leaf*
  first_sph_leaf_kid(stsp_estimator::cell_iterator const&C) {
    LoopStSpLeafKids(StSpCellIter,C,l,sph) return l;
    return 0;
  }
  //----------------------------------------------------------------------------
  inline stsp_cell*
  first_sph_cell_kid(stsp_estimator::cell_iterator const&C) {
    LoopStSpCellKids(StSpCellIter,C,c,sph) return c;
    return 0;
  }
  //----------------------------------------------------------------------------
  inline vect first_sph_kid_pos(stsp_estimator::cell_iterator const&C) {
    register stsp_cell*c = first_sph_cell_kid(C);
    if(c) return pos(c);
    register stsp_leaf*l = first_sph_leaf_kid(C);
    if(l) return pos(l);
    return center(C);
  }
  //============================================================================
  inline stsp_leaf*
  first_sticky_leaf_kid(stsp_estimator::cell_iterator const&C) {
    LoopStSpLeafKids(StSpCellIter,C,l,sticky) return l;
    return 0;
  }
  //----------------------------------------------------------------------------
  inline stsp_cell*
  first_sticky_cell_kid(stsp_estimator::cell_iterator const&C) {
    LoopStSpCellKids(StSpCellIter,C,c,sticky) return c;
    return 0;
  }
  //----------------------------------------------------------------------------
  inline 
  void first_sticky_kid_pos_and_vel(stsp_estimator::cell_iterator const&C,
				    vect&X, vect&V) {
    register stsp_cell*c = first_sticky_cell_kid(C);
    if(c) { X=pos(c); V=vel(c); return; }
    register stsp_leaf*l = first_sticky_leaf_kid(C);
    if(l) { X=pos(l); V=vel(l); return; }
  }
  //////////////////////////////////////////////////////////////////////////////
  // struct CountBodiesSticky                                                 //
  //////////////////////////////////////////////////////////////////////////////
  struct CountBodiesSticky {
    template<typename bodies_type>
    uint operator() (const bodies_type*const&B) const {
      register uint N = 0;
      for(register int i=0; i!=B->N_bodies(); ++i)
	if(is_in_tree(B->flg(i)) && is_sticky(B->flg(i))) ++N;
      return N;
    } };
  //////////////////////////////////////////////////////////////////////////////
  // struct UpdateLeafsSticky                                                 //
  //////////////////////////////////////////////////////////////////////////////
  struct UpdateLeafsSticky {
    const oct_tree              *T;
    mutable uint                 NS,NA;
    mutable stsp_leaf::leaf_data*Di;
    UpdateLeafsSticky(const oct_tree      *const&t,
		      stsp_leaf::leaf_data*const&d) : T(t), Di(d) {}
    template<typename bodies_type>
    uint operator() (const bodies_type*const&B) const {
      NS = 0;  NA = 0;
      LoopLeafs(stsp_leaf,T,Li) {
	Li->copy_from_bodies_flag(B);
	if(is_sticky(Li)) {
	  ++NS;
	  Li->set_data(Di++);
	  Li->set_sticky(B);
	  if(is_active(Li)) ++NA;
	}
      }
      return 0;
    } };
  //////////////////////////////////////////////////////////////////////////////
  // struct CountBodiesSPH                                                    //
  //////////////////////////////////////////////////////////////////////////////
  struct CountBodiesSPH {
    template<typename bodies_type>
    uint operator() (const bodies_type*const&B) const {
      return B->N_sph();
    } };
  //////////////////////////////////////////////////////////////////////////////
  // struct UpdateLeafsSPH                                                    //
  //////////////////////////////////////////////////////////////////////////////
  struct UpdateLeafsSPH {
    const oct_tree              *T;
    mutable uint                 NS,NA;
    mutable stsp_leaf::leaf_data*Di;
    UpdateLeafsSPH(const oct_tree      *const&t,
		   stsp_leaf::leaf_data*const&d) : T(t), Di(d) {}
    template<typename bodies_type>
    uint operator() (const bodies_type*const&B) const {
      NS = 0;  NA = 0;
      LoopLeafs(stsp_leaf,T,Li) {
	Li->copy_from_bodies_flag(B);
	if(is_sph(Li)) {
	  ++NS;
	  Li->set_data(Di++);
	  Li->set_sph(B);
	  if(is_active(Li)) ++NA;
	}
      }
      return 0;
    } };
  //////////////////////////////////////////////////////////////////////////////
}                                                  // END: unamed namespace     
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class nbdy::stsp_estimator                                                   
//                                                                              
////////////////////////////////////////////////////////////////////////////////
void stsp_estimator::update_leafs_sticky()
{
  if(TREE==0) error("stsp_estimator: no tree");    // IF no tree, FATAL ERROR   
  if(! TREE->is_used_for_stsp() )                  // IF tree not used by stsp  
    reset();                                       //   reset allocation & flags
  if(! STC_UPTODATE ) {                            // IF not up to date         
    // 1. count sticky bodies                                                   
    NL = TREE->UseBodies(CountBodiesSticky());
    ALL_STSP = NL == TREE->N_leafs();
    if(NL) {
      // 2. allocate leaf data for this number of sticky bodies                 
      if(LEAF_DATA) delete[] LEAF_DATA;
      LEAF_DATA = falcON_New(stsp_leaf::leaf_data,NL);
      // 3. Loop leafs, copy body flag; for sticky leafs: set mem, copy data    
      UpdateLeafsSticky ULS(TREE,LEAF_DATA);
      TREE->UseBodies(ULS);
      NL = ULS.NS;                                 // NL could be < B->N_sph()  
      ALL_STSP   = NL == TREE->N_leafs();
      ALL_ACTIVE = NL == ULS.NA;
    }
  }
  SPH_UPTODATE = 0;
}
//------------------------------------------------------------------------------
void stsp_estimator::update_leafs_sph() {
  if(TREE==0) error("stsp_estimator: no tree");    // IF no tree, FATAL ERROR   
  if(! TREE->is_used_for_stsp() )                  // IF tree not used by stsp  
    reset();                                       //   reset allocation & flags
  if(! SPH_UPTODATE) {
    NL = TREE->UseBodies(CountBodiesSPH());
    if(NL) {
      if(LEAF_DATA) delete[] LEAF_DATA;
      LEAF_DATA = falcON_New(stsp_leaf::leaf_data,NL);
      UpdateLeafsSPH ULS(TREE,LEAF_DATA);
      TREE->UseBodies(ULS);
      NL = ULS.NS;                                 // NL could be < B->N_sph()  
      ALL_STSP   = NL == TREE->N_leafs();
      ALL_ACTIVE = NL == ULS.NA;
    }
  }
  STC_UPTODATE = 0;
}
//------------------------------------------------------------------------------
void stsp_estimator::prepare_sph()
{
  // 1. allocate leaf memory, copy leaf data from bodies etc                    
  update_leafs_sph();                              // update flags & leafs      
  if(SPH_UPTODATE) return;                         // IF up to date: DONE       
  // 2. loop cells: pass flags & number of sph leafs up; count sph cells        
  register uint nc = 0;                            // counter: # sph cells      
  LoopCellsUp(cell_iterator,TREE,Ci) {             // LOOP cells upwards        
    register uint ns = 0;                          //   counter: # sph leafs    
    Ci->reset_active_flag();                       //   reset activity flag     
    LoopStSpLeafKids(cell_iterator,Ci,l,sph) {     //   LOOP stsp sub-leafs     
      Ci->add_active_flag(l);                      //     add in activity flag  
      ++ns;                                        //     count sph leafs       
    }                                              //   END LOOP                
    LoopStSpCellKids(cell_iterator,Ci,c,sph) {     //   LOOP stsp sub-cells c   
      Ci->add_active_flag(c);                      //     add in activity flag  
      ns += numb(c);                               //     count sph leafs       
    }                                              //   END LOOP                
    Ci->numb() = ns;                               //   set: # stsp leaf descs  
    if     (ns == 0u)                              //   IF no stsp leaf descs   
      Ci->un_set(flag::SPH);                       //     set flag: no stsp     
    else {                                         //   ELSE                    
      ++nc;                                        //     count stsp cells      
      Ci->add(flag::SPH);                          //     flag as stsp cell     
      if(ns == number(Ci)) Ci->add(flag::AL_SPH);  //     IF all stsp: flag     
    }                                              //   ENDIF                   
  }                                                // END LOOP                  
  NC = nc;                                         // # stsp cells              
  // 3. allocate memory for stsp cells                                          
  if(CELL_SRCE) delete[] CELL_SRCE;
  CELL_SRCE = falcON_New(stsp_cell::srce_data,NC);
  // 4. loop cells: give memory to stsp cells, pass up pos, size, rmax          
  register stsp_cell::srce_data*ci=CELL_SRCE;      // pter to cell's source     
  LoopCellsUp(cell_iterator,TREE,Ci)
  if(is_sph(Ci)) {                                 // LOOP sph ells upwards     
    Ci->set_srce(ci++);                            //   give srce memory        
    register vect x=zero;                          //   mean position           
    LoopStSpLeafKids(cell_iterator,Ci,l,sph)       //   LOOP stsp sub-leafs     
      x += pos(l);                                 //     add up: mean pos      
    LoopStSpCellKids(cell_iterator,Ci,c,sph)       //   LOOP stsp sub-cells     
      x.add_times(pos(c),numb(c));                 //     add up: mean pos      
    x /= real(numb(Ci));                           //   mean position           
    Ci->pos() = x;                                 //   set cell: mean position 
    register real s=zero,r=zero;                   //   size, r_max             
    LoopStSpLeafKids(cell_iterator,Ci,l,sph) {     //   LOOP stsp sub-leafs     
      register real R = dist(x,pos(l));            //     distance to sub leaf  
      update_max(r,R);                             //     update r_max          
      update_max(s,R+size(l));                     //     update size           
    }                                              //   END LOOP                
    LoopStSpCellKids(cell_iterator,Ci,c,sph) {     //   LOOP stsp sub-cells     
      register real R = dist(x,pos(c));            //     distance to sub cell  
      update_max(r,R+rmax(c));                     //     update r_max          
      update_max(s,R+size(c));                     //     update size           
    }                                              //   END LOOP                
    Ci->size() = s;                                //   set cell: size          
    Ci->rmax() = r;                                //   set cell: rmax          
  }                                                // END LOOP                  
  SPH_UPTODATE = 1;                                // tree up to date for sph   
}
//------------------------------------------------------------------------------
void stsp_estimator::prepare_sticky()
{
  // 1. allocate leaf memory, copy leaf data from bodies etc                    
  update_leafs_sticky();                           // update flags & leafs      
  if(STC_UPTODATE) return;                         // IF up to date: DONE       
  // 2. loop cells: pass flags & number of sticky leafs up; count sticky cells  
  register uint nc = 0;                            // counter: # sticky cells   
  LoopCellsUp(cell_iterator,TREE,Ci) {             // LOOP cells upwards        
    register uint ns = 0;                          //   counter: # sticky leafs 
    Ci->reset_active_flag();                       //   reset activity flag     
    LoopStSpLeafKids(cell_iterator,Ci,l,sticky) {  //   LOOP stsp sub-leafs     
      Ci->add_active_flag(l);                      //     add in activity flag  
      ++ns;                                        //     count sticky leafs    
    }                                              //   END LOOP                
    LoopStSpCellKids(cell_iterator,Ci,c,sticky) {  //   LOOP stsp sub-cells     
      Ci->add_active_flag(c);                      //     add in activity flag  
      ns += numb(c);                               //     count sticky leafs    
    }                                              //   END LOOP                
    Ci->numb() = ns;                               //   set: # stsp leaf descs  
    if     (ns == 0u)                              //   IF no stsp leaf descs   
      Ci->un_set(flag::STICKY);                    //     set flag: no stsp     
    else {                                         //   ELSE                    
      ++nc;                                        //     count stsp cells      
      Ci->add(flag::STICKY);                       //     flag as stsp cell     
      if(ns == number(Ci)) Ci->add(flag::AL_STICKY); //   IF all stsp: flag     
    }                                              //   ENDIF                   
  }                                                // END LOOP                  
  NC = nc;                                         // # stsp cells              
  // 3. allocate memory for stsp cells                                          
  if(CELL_SRCE) delete[] CELL_SRCE;
  CELL_SRCE = falcON_New(stsp_cell::srce_data,NC);
  // 4. loop cells: give memory to stsp cells, pass up pos, vel, size, vrad     
  register stsp_cell::srce_data*ci=CELL_SRCE;      // pter to cell's source     
  LoopCellsUp(cell_iterator,TREE,Ci)
  if(is_sticky(Ci)) {                              // LOOP sticky ells upwards  
    Ci->set_srce(ci++);                            //   give srce memory        
    register vect x=zero;                          //   mean position           
    register vect v=zero;                          //   mean velocity           
    LoopStSpLeafKids(cell_iterator,Ci,l,sticky) {  //   LOOP stsp sub-leafs     
      x += pos(l);                                 //     add up: mean pos      
      v += vel(l);                                 //     add up: mean vel      
    }                                              //   END LOOP                
    LoopStSpCellKids(cell_iterator,Ci,c,sticky) {  //   LOOP stsp sub-cells     
      x.add_times(pos(c),numb(c));                 //     add up: mean pos      
      v.add_times(vel(c),numb(c));                 //     add up: mean vel      
    }                                              //   END LOOP                
    register real iN=one/real(numb(Ci));           //   1/N_sticky              
    x *= iN;                                       //   mean position           
    v *= iN;                                       //   mean velocity           
    Ci->pos() = x;                                 //   set cell: mean position 
    Ci->vel() = v;                                 //   set cell: mean position 
    register real s=zero,r=zero;                   //   size, vrad              
    LoopStSpLeafKids(cell_iterator,Ci,l,sticky) {  //   LOOP stsp sub-leafs     
      update_max(s,dist(x,pos(l))+size(l));        //     update size           
      update_max(r,dist(v,vel(l)));                //     update vrad           
    }                                              //   END LOOP                
    LoopStSpCellKids(cell_iterator,Ci,c,sticky) {  //   LOOP stsp sub-cells     
      update_max(s,dist(x,pos(c))+size(c));        //     update size           
      update_max(r,dist(v,vel(c))+vrad(c));        //     update vrad           
    }                                              //   END LOOP                
    Ci->size() = s;                                //   set cell: size          
    Ci->vrad() = r;                                //   set cell: vrad          
  }                                                // END LOOP                  
  STC_UPTODATE = 1;                                // tree up to date for stc   
}
//------------------------------------------------------------------------------
void stsp_estimator::make_sticky_list (elem_pair *bl,
				       uint const&nl,
				       uint      &na,
				       real const&tau)
{
  prepare_sticky();
  sticky_finder sfind(tau,nl,bl);
  MutualInteractor<sticky_finder> MI(&sfind,TREE->depth());
  MI.cell_self(root());
  na = sfind.actual_size_of_list();
  TREE->mark_stsp_usage();
}
//------------------------------------------------------------------------------
void stsp_estimator::make_sph_list (elem_pair *bl,
				    uint const&nl,
				    uint      &na,
				    bool const&Max)
{
  prepare_sph();
  if(Max) {
    neighbour_finder sfind(nl,bl);
    MutualInteractor<neighbour_finder> MI(&sfind,TREE->depth());
    MI.cell_self(root());
    na = sfind.actual_size_of_list();
  } else {
    partner_finder sfind(nl,bl);
    MutualInteractor<partner_finder> MI(&sfind,TREE->depth());
    MI.cell_self(root());
    na = sfind.actual_size_of_list();
  }
  TREE->mark_stsp_usage();
}
//------------------------------------------------------------------------------
template<typename bodies_type>
inline void stsp_estimator::copy_to_bodies_num(const bodies_type*const&B) const
{
  if(ALL_STSP)
    LoopLeafs(leaf_type,TREE,Li)
      Li->copy_to_bodies_num(B);
  else
    LoopLeafs(leaf_type,TREE,Li) if(is_sph(Li))
      Li->copy_to_bodies_num(B);
}
//------------------------------------------------------------------------------
void stsp_estimator::count_sph_partners()
{
  prepare_sph();
  partner_counter scount;
  MutualInteractor<partner_counter> MI(&scount,TREE->depth());
  MI.cell_self(root());
  if     (TREE->use_sbodies()) copy_to_bodies_num(TREE->my_sbodies());
#ifdef falcON_MPI
  else if(TREE->use_pbodies()) copy_to_bodies_num(TREE->my_pbodies());
#endif
  else if(TREE->use_abodies()) copy_to_bodies_num(TREE->my_abodies());
  else falcON_Error("tree has neither bodies nor arrays for data");
  TREE->mark_stsp_usage();
}
////////////////////////////////////////////////////////////////////////////////

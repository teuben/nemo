//------------------------------------------------------------------------------
//                                                                              
/// \file src/public/partner.cc                                                 
//                                                                              
// Copyright (C) 2000-2008  Walter Dehnen                                       
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
//------------------------------------------------------------------------------
#include <public/partner.h>
#include <public/interact.h>

namespace {
  using namespace falcON;
  typedef PartnerEstimator::indx_pair indx_pair;
  //////////////////////////////////////////////////////////////////////////////
#define LoopPartnerCKids(CELLITER,CELL,NAME,STSP)		\
  LoopCellKids(CELLITER,CELL,NAME) if(is_##STSP(NAME))
#define LoopPartnerLKids(CELLITER,CELL,NAME,STSP)		\
  LoopLeafKids(CELLITER,CELL,NAME) if(is_##STSP(NAME))
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // struct take_sph                                                            
  // struct take_sticky                                                         
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  struct take_sph {
    static bool take    (PartnerEstimator::leaf_iterator const&A) 
    { return is_sph(A); }
    static bool take    (PartnerEstimator::cell_iterator const&A)
    { return is_sph(A); }
    static bool take_all(PartnerEstimator::cell_iterator const&A)
    { return al_sph(A); }
  };
  //////////////////////////////////////////////////////////////////////////////
  struct take_sticky {
    static bool take    (PartnerEstimator::leaf_iterator const&A) 
    { return is_sticky(A); }
    static bool take    (PartnerEstimator::cell_iterator const&A)
    { return is_sticky(A); }
    static bool take_all(PartnerEstimator::cell_iterator const&A)
    { return al_sticky(A); }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class BasicFinder                                                          
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  template<typename taker>
  class BasicFinder : 
    public taker,
    public BasicIactor<PartnerEstimator> {
    //--------------------------------------------------------------------------
  public:
    taker::take;
    taker::take_all;
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
	  for(leaf_iter B=B0; B!=BN; ++B)
	    check_pair(A,B);
	} else {
	  for(leaf_iter B=B0; B!=BN; ++B)
	    if(is_active(B)) check_pair(A,B);
	}
      } else {
	if(all_active) {
	  for(leaf_iter B=B0; B!=BN; ++B)
	    if(take(B)) check_pair(A,B);
	} else {
	  for(leaf_iter B=B0; B!=BN; ++B)
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
    //--------------------------------------------------------------------------
    explicit
    BasicFinder(const unsigned dir[4] = Default::direct)
      : BasicIactor<PartnerEstimator>(dir) {}
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class BasicLister                                                          
  //                                                                            
  // for finding leaf pairs in a sticky_tree;                                   
  // derived from BasicFinder above, which satisfies the requirements           
  // for an INTERACTOR template parameter to class MutualInteractor<>;          
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  template<typename taker>
  class BasicLister : public BasicFinder<taker> {
    //--------------------------------------------------------------------------
    // data of class BasicLister                                                
    //--------------------------------------------------------------------------
  protected:
    const bodies    *BODIES;                       // bodies (to get index)     
  private:
    const unsigned   MAX;                          // maximal size of list      
    indx_pair       *BL;                           // list of interaction pairs 
    mutable unsigned N;                            // actual size of list       
    //--------------------------------------------------------------------------
    // other methods                                                            
    //--------------------------------------------------------------------------
  protected:
    void add_pair(PartnerEstimator::leaf_iterator const&A, 
		  PartnerEstimator::leaf_iterator const&B) const
    {
      if(N<MAX) {
	if(BODIES->is_less(mybody(A), mybody(B))) {
	  BL[N][0] = mybody(A);
	  BL[N][1] = mybody(B);
	} else {
	  BL[N][0] = mybody(B);
	  BL[N][1] = mybody(A);
	}
      }
      N++;
      if(N==MAX) falcON_Warning("interaction list overflow");
    }
    //--------------------------------------------------------------------------
    // abstract methods                                                         
    //--------------------------------------------------------------------------
    virtual void check_pair(PartnerEstimator::leaf_iterator const&, 
			    PartnerEstimator::leaf_iterator const&) const = 0;
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    unsigned const &actual_size_of_list() const { return N; }
    //--------------------------------------------------------------------------
    BasicLister(const bodies  *b,
		const unsigned n,
		indx_pair     *l,
		const unsigned dir[4] = Default::direct)
      : BasicFinder<taker>(dir), BODIES(b), MAX(n), BL(l), N(0) {}
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // struct Counter<bool>                                                       
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  template<bool> struct Counter;
  template<> struct Counter<true> {
    static void count(PartnerEstimator::leaf_iterator const&A,
		      PartnerEstimator::leaf_iterator const&B) {
      if(is_active(A)) A->inc();
      if(is_active(B)) B->inc();
    } };
  template<> struct Counter<false> {
    static void count(PartnerEstimator::leaf_iterator const&A,
		      PartnerEstimator::leaf_iterator const&B) {}
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class StickyFinder<>                                                       
  //                                                                            
  // for finding stsp_leaf pairs which satisfy:                                 
  // - both are flagged as sticky                                               
  // - at least one is flagged as active                                        
  // - there is a time t in [0,tau] such that:                                  
  //   | (x_i+v_i*t) - (x_j+v+j*t) | < size_i+size_j                            
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  template<bool COUNT>
  class StickyFinder : public BasicLister<take_sticky> {
    const real    TAU;                             // time period               
    //--------------------------------------------------------------------------
  protected:
    void check_pair(leaf_iter const&A, leaf_iter const&B) const
    {
      vect R = pos(A)-pos(B);                      // distance vector           
      real Sq= square(size(A)+size(B));            // (combined size)^2         
      if(norm(R) < Sq) {                           // IF(overlap at t=0):       
	add_pair(A,B);                             //   add                     
	Counter<COUNT>::count(A,B);                //   count collision partners
	return;                                    //   DONE                    
      }                                            // ENDIF                     
      if(TAU==zero) return;                        // IF(no time > 0) DONE      
      vect V = vel(A)-vel(B);                      // velocity difference       
      real RV= R*V;                                // scalar product            
      if(RV > zero) return;                        // IF(diverging orbits) DONE 
      real t = min(TAU,-RV/norm(V));               // time of min distance      
      if(norm(R+t*V) < Sq) {                       // IF(overlap)               
	add_pair(A,B);                             //   add                     
	Counter<COUNT>::count(A,B);                //   count collision partners
      }                                            // ENDIF                     
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, leaf_iter const&B) const
      // true  if there cannot possibly be an interaction between any A and B   
      // false if an interaction between any A and B is possible                
    {
      vect R = pos(A)-pos(B);                      // R   = distance vector     
      real
	Rq = norm(R),                              // R^2 = distance^2          
	x  = size(A)+size(B);                      // x   = combined size       
      if(Rq < x*x) return false;                   // overlap at t=0            
      if(TAU==zero) return true;                   // no times t>0              
      vect V = vel(A)-vel(B);                      // V   = velocity diff       
      real
	w  = vrad(A),                              // w   = v-size              
	wq = w*w,                                  // w^2 = (v-size)^2          
	RV = R*V,                                  // R*V = scalar product      
	RVq= RV*RV;                                // (R*V)^2                   
      if(RV>zero && RVq>wq*Rq) return true;        // R*V/|R| > w: diverging    
      real
	Vq = norm(V),                              // V^2 = (v-distance)^2      
	t  = (wq>=Vq)? TAU :                       // v-size too large          
	min(TAU, (w*real(sqrt((Rq*Vq-RVq)/(Vq-wq)))-RV)/Vq );// t of min dist   
      if(norm(R+t*V) < square(x+t*w)) return false;// possible overlap          
      return true;                                 // no overlap -> discard     
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, cell_iter const&B) const
      // true  if there cannot possibly be an interaction between any A and B   
      // false if an interaction between any A and B is possible                
    {
      vect R = pos(A)-pos(B);                      // R   = distance vector     
      real
	Rq = norm(R),                              // R^2 = distance^2          
	x  = size(A)+size(B);                      // x   = combined size       
      if(Rq < x*x) return false;                   // overlap at t=0            
      if(TAU==zero) return true;                   // no times t>0              
      vect V = vel(A)-vel(B);                      // V   = velocity diff       
      real
	w  = vrad(A)+vrad(B),                      // w   = combined v-size     
	wq = w*w,                                  // w^2 = (comb v-size)^2     
	RV = R*V,                                  // R*V = scalar product      
	RVq= RV*RV;                                // (R*V)^2                   
      if(RV>zero && RVq>wq*Rq) return true;        // R*V/|R| > w: diverging    
      real
	Vq = norm(V),                              // V^2 = (v-distance)^2      
	t  = (wq>=Vq)? TAU :                       // v-size too large          
	min(TAU, (w*real(sqrt((Rq*Vq-RVq)/(Vq-wq)))-RV)/Vq );// t of min dist   
      if(norm(R+t*V) < square(x+t*w)) return false;// possible overlap          
      return true;                                 // no overlap -> discard     
    }
    //--------------------------------------------------------------------------
  public:
    StickyFinder(const bodies  *b,
		 const real     t,
		 const unsigned n,
		 indx_pair     *l,
		 const unsigned dir[4] = Default::direct) :
      BasicLister<take_sticky>(b,n,l,dir), TAU(t) {}
    //--------------------------------------------------------------------------
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || size(A)+TAU*vrad(A) > size(B)+TAU*vrad(B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class NeighbourCounter                                                     
  //                                                                            
  // for counting PartnerEstimator::Leaf pairs which satisfy:                   
  // - both are flagged as sph                                                  
  // - at least one is flagged as active                                        
  // - | x_i - x_j | < max(size_i,size_j)                                       
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class NeighbourCounter : public BasicFinder<take_sph> {
  protected:
    void check_pair(leaf_iter const&A, leaf_iter const&B) const
    {
      real Rq = dist_sq(pos(A),pos(B));
      if(Rq < sizeq(A) || Rq < sizeq(B))
	Counter<1>::count(A,B);
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
    //--------------------------------------------------------------------------
  public:
    NeighbourCounter(const unsigned dir[4] = Default::direct)
      : BasicFinder<take_sph>(dir) {}
    //--------------------------------------------------------------------------
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || size(A) > size(B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class NeighbourLister<>                                                    
  //                                                                            
  // for finding and counting PartnerEstimator::Leaf pairs which satisfy:       
  // - both are flagged as sph                                                  
  // - at least one is flagged as active                                        
  // - | x_i - x_j | < max(size_i,size_j)                                       
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  template<bool COUNT>
  class NeighbourLister : public BasicLister<take_sph> {
  protected:
    void check_pair(leaf_iter const&A, leaf_iter const&B) const
    {
      real Rq = dist_sq(pos(A),pos(B));
      if(Rq < sizeq(A) || Rq < sizeq(B)) {
	add_pair(A,B);
	Counter<COUNT>::count(A,B);
      }
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
    //--------------------------------------------------------------------------
  public:
    NeighbourLister(const bodies  *b,
		    const unsigned n,
		    indx_pair     *l,
		    const unsigned dir[4] = Default::direct)
      : BasicLister<take_sph>(b,n,l,dir) {}
    //--------------------------------------------------------------------------
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || size(A) > size(B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class PartnerCounter                                                       
  //                                                                            
  // for counting PartnerEstimator::Leaf pairs which satisfy:                   
  // - both are flagged as sph                                                  
  // - | x_i - x_j | < size_i + size_j                                          
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class PartnerCounter : public BasicFinder<take_sph> {
  protected:
    void check_pair(leaf_iter const&A, leaf_iter const&B) const
    {
      if(dist_sq(pos(A),pos(B)) < square(size(A)+size(B)))
	Counter<1>::count(A,B);
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
    //--------------------------------------------------------------------------
    explicit PartnerCounter(const unsigned dir[4] = Default::direct)
      : BasicFinder<take_sph>(dir) {}
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class PartnerLister<>                                                      
  //                                                                            
  // for finding and counting PartnerEstimator::Leaf pairs which satisfy:       
  // - both are flagged as sph                                                  
  // - at least one is flagged as active                                        
  // - | x_i - x_j | < size_i + size_j                                          
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  template<bool COUNT>
  class PartnerLister : public BasicLister<take_sph> {
  protected:
    void check_pair(leaf_iter const&A, leaf_iter const&B) const
    {
      if(dist_sq(pos(A),pos(B)) < square(size(A)+size(B))) {
	add_pair(A,B);
	Counter<COUNT>::count(A,B);
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
    //--------------------------------------------------------------------------
  public:
    PartnerLister(const bodies  *b,
		  const unsigned n,
		  indx_pair     *l,
		  const unsigned dir[4] = Default::direct)
      : BasicLister<take_sph>(b,n,l,dir) {}
    //--------------------------------------------------------------------------
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || size(A) > size(B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace {
////////////////////////////////////////////////////////////////////////////////
falcON_TRAITS(StickyFinder<0>,"StickyFinder<0>");
falcON_TRAITS(StickyFinder<1>,"StickyFinder<1>");
falcON_TRAITS(NeighbourCounter,"NeighbourCounter");
falcON_TRAITS(NeighbourLister<0>,"NeighbourLister<0>");
falcON_TRAITS(NeighbourLister<1>,"NeighbourLister<1>");
falcON_TRAITS(PartnerCounter,"PartnerCounter");
falcON_TRAITS(PartnerLister<0>,"PartnerLister<0>");
falcON_TRAITS(PartnerLister<1>,"PartnerLister<1>");
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class falcON::PartnerEstimator                                               
//                                                                              
////////////////////////////////////////////////////////////////////////////////
void PartnerEstimator::update_leafs_sticky()
{
  if(TREE==0)
    falcON_Error("PartnerEstimator: no tree");     // IF no tree, FATAL ERROR   
  if(! TREE->is_used_for_stsp() )                  // IF tree not used by stsp  
    reset();                                       //   reset allocation & flags
  if(! STC_UPTODATE ) {                            // IF not up to date         
    NL = TREE->my_bodies()->N_sph();
    ALL_STSP = NL == TREE->N_leafs();
    if(NL) {
      if(LEAF_DATA) falcON_DEL_A(LEAF_DATA);
      LEAF_DATA = falcON_NEW(Leaf::leaf_data,NL);
      Leaf::leaf_data* Di = LEAF_DATA;
      unsigned NS=0, NA=0;
      LoopLeafs(Leaf,TREE,Li) {
	Li->copy_from_bodies_flag(TREE->my_bodies());
	if(is_sticky(Li)) {
	  ++NS;
	  Li->set_data(Di++);
	  Li->set_sticky(TREE->my_bodies());
	  if(is_active(Li)) ++NA;
	}
      }
      if(NS > NL) falcON_Error("PartnerEstimator: too many sticky leafs");
      NL = NS;
      ALL_STSP   = NL == TREE->N_leafs();
      ALL_ACTIVE = NL == NA;
    }
  }
  SPH_UPTODATE = 0;
}
//------------------------------------------------------------------------------
void PartnerEstimator::update_leafs_sph() {
  if(TREE==0)
    falcON_Error("PartnerEstimator: no tree");     // IF no tree, FATAL ERROR   
  if(! TREE->is_used_for_stsp() )                  // IF tree not used by stsp  
    reset();                                       //   reset allocation & flags
  if(! SPH_UPTODATE) {
    NL = TREE->my_bodies()->N_sph();
    ALL_STSP = NL == TREE->N_leafs();
    if(NL) {
      if(LEAF_DATA) falcON_DEL_A(LEAF_DATA);
      LEAF_DATA = falcON_NEW(Leaf::leaf_data,NL);
      Leaf::leaf_data* Di = LEAF_DATA;
      unsigned NS=0, NA=0;
      LoopLeafs(Leaf,TREE,Li) {
	Li->copy_from_bodies_flag(TREE->my_bodies());
	if(is_sph(Li)) {
	  ++NS;
	  Li->set_data(Di++);
	  Li->set_sph(TREE->my_bodies());
	  if(is_active(Li)) ++NA;
	}
      }
      if(NS > NL) falcON_Error("PartnerEstimator: too many sticky leafs");
      NL = NS;
      ALL_STSP   = NL == TREE->N_leafs();
      ALL_ACTIVE = NL == NA;
    }
  }
  STC_UPTODATE = 0;
}
//------------------------------------------------------------------------------
void PartnerEstimator::prepare_sph()
{
  // 1. allocate leaf memory, copy leaf data from bodies etc                    
  update_leafs_sph();                              // update flags & leafs      
  if(SPH_UPTODATE) return;                         // IF up to date: DONE       
  // 2. loop cells: pass flags & number of sph leafs up; count sph cells        
  unsigned nc = 0;                                 // counter: # sph cells      
  LoopCellsUp(cell_iterator,TREE,Ci) {             // LOOP cells upwards        
    unsigned ns = 0;                               //   counter: # sph leafs    
    Ci->reset_active_flag();                       //   reset activity flag     
    LoopPartnerLKids(cell_iterator,Ci,l,sph) {     //   LOOP partner sub-leafs  
      Ci->add_active_flag(l);                      //     add in activity flag  
      ++ns;                                        //     count sph leafs       
    }                                              //   END LOOP                
    LoopPartnerCKids(cell_iterator,Ci,c,sph) {     //   LOOP partner sub-cells c
      Ci->add_active_flag(c);                      //     add in activity flag  
      ns += numb(c);                               //     count sph leafs       
    }                                              //   END LOOP                
    Ci->numb() = ns;                               //   set: # partner leaf desc
    if(ns == 0u)                                   //   IF no partner leaf descs
      Ci->un_set(flags::sph);                      //     set flag: no partner  
    else {                                         //   ELSE                    
      ++nc;                                        //     count partner cells   
      Ci->add(flags::sph);                         //     flag as partner cell  
      if(ns == number(Ci))                         //     IF all partner:       
	Ci->add(flags::all_sph);                   //       flag as such        
    }                                              //   ENDIF                   
  }                                                // END LOOP                  
  NC = nc;                                         // # partner cells           
  // 3. allocate memory for Cell::srce_data                                     
  if(CELL_SRCE) falcON_DEL_A(CELL_SRCE);
  CELL_SRCE = falcON_NEW(Cell::srce_data,NC);
  // 4. loop cells: give memory to PartnerCells, pass up pos, size, rmax        
  Cell::srce_data*ci=CELL_SRCE;                    // pter to cell's source     
  LoopCellsUp(cell_iterator,TREE,Ci)
  if(is_sph(Ci)) {                                 // LOOP sph ells upwards     
    Ci->set_srce(ci++);                            //   give srce memory        
    vect x(zero);                                  //   mean position           
    LoopPartnerLKids(cell_iterator,Ci,l,sph)       //   LOOP partner sub-leafs  
      x += pos(l);                                 //     add up: mean pos      
    LoopPartnerCKids(cell_iterator,Ci,c,sph)       //   LOOP partner sub-cells  
      x += numb(c) * pos(c);                       //     add up: mean pos      
    x /= real(numb(Ci));                           //   mean position           
    Ci->pos() = x;                                 //   set cell: mean position 
    real s=zero,r=zero;                            //   size, r_max             
    LoopPartnerLKids(cell_iterator,Ci,l,sph) {     //   LOOP partner sub-leafs  
      real R = dist(x,pos(l));                     //     distance to sub leaf  
      update_max(r,R);                             //     update r_max          
      update_max(s,R+size(l));                     //     update size           
    }                                              //   END LOOP                
    LoopPartnerCKids(cell_iterator,Ci,c,sph) {     //   LOOP partner sub-cells  
      real R = dist(x,pos(c));                     //     distance to sub cell  
      update_max(r,R+rmax(c));                     //     update r_max          
      update_max(s,R+size(c));                     //     update size           
    }                                              //   END LOOP                
    Ci->size() = s;                                //   set cell: size          
    Ci->rmax() = r;                                //   set cell: rmax          
  }                                                // END LOOP                  
  SPH_UPTODATE = 1;                                // tree up to date for sph   
}
//------------------------------------------------------------------------------
void PartnerEstimator::prepare_sticky()
{
  // 1. allocate leaf memory, copy leaf data from bodies etc                    
  update_leafs_sticky();                           // update flags & leafs      
  if(STC_UPTODATE) return;                         // IF up to date: DONE       
  // 2. loop cells: pass flags & number of sticky leafs up; count sticky cells  
  unsigned nc = 0;                                 // counter: # sticky cells   
  LoopCellsUp(cell_iterator,TREE,Ci) {             // LOOP cells upwards        
    unsigned ns = 0;                               //   counter: # sticky leafs 
    Ci->reset_active_flag();                       //   reset activity flag     
    LoopPartnerLKids(cell_iterator,Ci,l,sticky) {  //   LOOP partner sub-leafs  
      Ci->add_active_flag(l);                      //     add in activity flag  
      ++ns;                                        //     count sticky leafs    
    }                                              //   END LOOP                
    LoopPartnerCKids(cell_iterator,Ci,c,sticky) {  //   LOOP partner sub-cells  
      Ci->add_active_flag(c);                      //     add in activity flag  
      ns += numb(c);                               //     count sticky leafs    
    }                                              //   END LOOP                
    Ci->numb() = ns;                               //   set: # partner leaf desc
    if     (ns == 0u)                              //   IF no partner leaf descs
      Ci->un_set(flags::sticky);                   //     set flag: no partner  
    else {                                         //   ELSE                    
      ++nc;                                        //     count partner cells   
      Ci->add(flags::sticky);                      //     flag as partner cell  
      if(ns == number(Ci))                         //     IF all partner:       
	Ci->add(flags::all_sticky);                //       flag as such        
    }                                              //   ENDIF                   
  }                                                // END LOOP                  
  NC = nc;                                         // # partner cells           
  // 3. allocate memory for partner cells                                       
  if(CELL_SRCE) falcON_DEL_A(CELL_SRCE);
  CELL_SRCE = falcON_NEW(Cell::srce_data,NC);
  // 4. loop cells: give memory to partner cells, pass up pos, vel, size, vrad  
  Cell::srce_data*ci=CELL_SRCE;                    // pter to cell's source     
  LoopCellsUp(cell_iterator,TREE,Ci)
  if(is_sticky(Ci)) {                              // LOOP sticky ells upwards  
    Ci->set_srce(ci++);                            //   give srce memory        
    vect x(zero);                                  //   mean position           
    vect v(zero);                                  //   mean velocity           
    LoopPartnerLKids(cell_iterator,Ci,l,sticky) {  //   LOOP partner sub-leafs  
      x += pos(l);                                 //     add up: mean pos      
      v += vel(l);                                 //     add up: mean vel      
    }                                              //   END LOOP                
    LoopPartnerCKids(cell_iterator,Ci,c,sticky) {  //   LOOP partner sub-cells  
      x += numb(c) * pos(c);                       //     add up: mean pos      
      v += numb(c) * vel(c);                       //     add up: mean vel      
    }                                              //   END LOOP                
    real iN=one/real(numb(Ci));                    //   1/N_sticky              
    x *= iN;                                       //   mean position           
    v *= iN;                                       //   mean velocity           
    Ci->pos() = x;                                 //   set cell: mean position 
    Ci->vel() = v;                                 //   set cell: mean position 
    real s=zero,r=zero;                            //   size, vrad              
    LoopPartnerLKids(cell_iterator,Ci,l,sticky) {  //   LOOP partner sub-leafs  
      update_max(s,dist(x,pos(l))+size(l));        //     update size           
      update_max(r,dist(v,vel(l)));                //     update vrad           
    }                                              //   END LOOP                
    LoopPartnerCKids(cell_iterator,Ci,c,sticky) {  //   LOOP partner sub-cells  
      update_max(s,dist(x,pos(c))+size(c));        //     update size           
      update_max(r,dist(v,vel(c))+vrad(c));        //     update vrad           
    }                                              //   END LOOP                
    Ci->size() = s;                                //   set cell: size          
    Ci->vrad() = r;                                //   set cell: vrad          
  }                                                // END LOOP                  
  STC_UPTODATE = 1;                                // tree up to date for stc   
}
//------------------------------------------------------------------------------
inline void PartnerEstimator::copy_to_bodies_num(bool sph) const
{
  if(sph)
    LoopLeafs(leaf_type,TREE,Li) if(is_sph(Li)    && is_active(Li))
      Li->copy_to_bodies_num(TREE->my_bodies());
  else
    LoopLeafs(leaf_type,TREE,Li) if(is_sticky(Li) && is_active(Li))
      Li->copy_to_bodies_num(TREE->my_bodies());
}
//------------------------------------------------------------------------------
template<bool COUNT>
void PartnerEstimator::make_st_list (indx_pair*bl,
				     unsigned  nl,
				     unsigned &na,
				     real      tau) falcON_THROWING
{
  prepare_sticky();
  StickyFinder<COUNT> sfind(TREE->my_bodies(),tau,nl,bl);
  MutualInteractor< StickyFinder<COUNT> > MI(&sfind,TREE->depth());
  MI.cell_self(root());
  na = sfind.actual_size_of_list();
  TREE->mark_stsp_usage();
}
//------------------------------------------------------------------------------
void PartnerEstimator::make_sticky_list (indx_pair*bl,
					 unsigned  nl,
					 unsigned &na,
					 real      tau,
					 bool      count) falcON_THROWING
{
  if(count && ! (TREE->my_bodies()->have(fieldbit::N))) {
    falcON_Warning("PartnerEstimator: cannot count: field 'N' not supported\n");
    count = false;
  }
  if(count) {
    make_st_list<1>(bl,nl,na,tau);
    copy_to_bodies_num(false);
  } else
    make_st_list<0>(bl,nl,na,tau);
}
//------------------------------------------------------------------------------
template<bool COUNT>
void PartnerEstimator::make_sp_list (indx_pair*bl,
				     unsigned  nl,
				     unsigned &na,
				     bool      Max) falcON_THROWING
{
  prepare_sph();
  if(Max) {
    NeighbourLister<COUNT> sfind(TREE->my_bodies(),nl,bl);
    MutualInteractor<NeighbourLister<COUNT> > MI(&sfind,TREE->depth());
    MI.cell_self(root());
    na = sfind.actual_size_of_list();
  } else {
    PartnerLister<COUNT> sfind(TREE->my_bodies(),nl,bl);
    MutualInteractor<PartnerLister<COUNT> > MI(&sfind,TREE->depth());
    MI.cell_self(root());
    na = sfind.actual_size_of_list();
  }
  TREE->mark_stsp_usage();
}
//------------------------------------------------------------------------------
void PartnerEstimator::make_sph_list(indx_pair*bl,
				     unsigned  nl,
				     unsigned &na,
				     bool      Max,
				     bool      count) falcON_THROWING
{
  if(count && ! (TREE->my_bodies()->have(fieldbit::N))) {
    falcON_Warning("PartnerEstimator: cannot count: field 'N' not supported\n");
    count = false;
  }
  if(count) {
    make_sp_list<1>(bl,nl,na,Max);
    copy_to_bodies_num(true);
  } else
    make_sp_list<0>(bl,nl,na,Max);
}
//------------------------------------------------------------------------------
void PartnerEstimator::count_sph_partners(bool Max) falcON_THROWING
{
  prepare_sph();
  if(Max) {
    NeighbourCounter scount;
    MutualInteractor<NeighbourCounter> MI(&scount,TREE->depth());
    MI.cell_self(root());
  } else {
    PartnerCounter scount;
    MutualInteractor<PartnerCounter> MI(&scount,TREE->depth());
    MI.cell_self(root());
  }
  copy_to_bodies_num(true);
  TREE->mark_stsp_usage();
}
////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------+
//                                                                             |
// stic.cc                                                                     |
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
// support for sticky particles                                                |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/stic.h>
#include <public/tree.h>
#include <public/iact.h>

namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::basic_finder                                                   
  //                                                                            
  // for finding soul pairs in a sticky_tree;                                   
  // derived from basic_iactor of iact.h, which satisfies the requirements for  
  // an INTERACTOR template parameter to class MutualInteractor<>;              
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class basic_finder : public basic_iactor<sticky_tree> {
    //--------------------------------------------------------------------------
    // types of class basic_finder                                              
    //--------------------------------------------------------------------------
  public:
    typedef uint elem_pair[2];                     // pair of indices           
    //--------------------------------------------------------------------------
    // data of class basic_finder                                               
    //--------------------------------------------------------------------------
  private:
    const uint    MAX;                             // maximal size of list      
    elem_pair    *BL;                              // list of interaction pairs 
    mutable uint  N;                               // actual size of list       
    //--------------------------------------------------------------------------
    // other methods                                                            
    //--------------------------------------------------------------------------
  protected:
    void add_pair(soul_iter const&A, soul_iter const&B) const
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
  public:
    uint const &actual_size_of_list() const { return N; }
    //--------------------------------------------------------------------------
    // constructors                                                             
    //--------------------------------------------------------------------------
    basic_finder(const uint n, elem_pair*l) : MAX(n),BL(l),N(0){}
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::sticky_finder                                                  
  //                                                                            
  // for finding sticky_soul pairs which satisfy:                               
  // - at least one is flagged as active                                        
  // - there is a time t in [0,tau] such that:                                  
  //   | (x_i+v_i*t) - (x_j+v+j*t) | < size_i+size_j                            
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
#define CHECK_PAIR							       \
  R = pos(A)-pos(B);                               /* distance vector       */ \
  xq= square(size(A)+size(B));                     /* (combined size)^2     */ \
  if(norm(R) < xq) { add_pair(A,B); continue; }    /* overlap at t=0        */ \
  if(TAU==zero) continue;                          /* no times t>0          */ \
  V = vel(A)-vel(B);                               /* velocity difference   */ \
  RV= R*V;                                         /* scalar product        */ \
  if(RV>zero) continue;                            /* diverging orbits      */ \
  t = min(TAU,-RV/norm(V));                        /* time of min distance  */ \
  if(norm(R+t*V) < xq) add_pair(A,B);              /* overlap: add to list  */
  //////////////////////////////////////////////////////////////////////////////
  class sticky_finder : public basic_finder {
  private:
    const real    TAU;                             // time period               
    //--------------------------------------------------------------------------
  protected:
    void many(bool const&all, soul_iter const&A,
	      soul_iter const&B0, soul_iter const&BN) const
    {
      register vect      R,V;
      register real      xq,RV,t;
      if(all) for(register soul_iter B=B0; B<BN; B++)
	{ CHECK_PAIR }
      else    for(register soul_iter B=B0; B<BN; B++) if(is_active(B))
	{ CHECK_PAIR }
    }
#undef CHECK_PAIR
    //--------------------------------------------------------------------------
    void single(soul_iter const&A, soul_iter const&B) const
    {
      if(!is_active(A) && !is_active(B)) return;   // no actives -> no job      
      register vect R = pos(A)-pos(B);             // R   = distance vector     
      register real xq= square(size(A)+size(B));   // x^2 = (combined size)^2   
      if(norm(R) < xq) return add_pair(A,B);       // |R|<x: overlap at t=0     
      if(TAU==zero) return;                        // no times t>0              
      register vect V = vel(A)-vel(B);             // V   = velocity diff       
      register real RV= R*V;                       // R*V = scalar product      
      if(RV>zero) return;                          // R*V>0: diverging orbits   
      register real t = min(TAU, -RV/norm(V));     // t of minimum distance     
      if(norm(R+t*V) < xq) return add_pair(A,B);   // d_min<x: overlapping      
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, soul_iter const&B) const
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
      basic_finder(n,l), TAU(t) {}
    //--------------------------------------------------------------------------
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || size(A)+TAU*vrad(A) > size(B)+TAU*vrad(B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::neighbour_finder                                               
  //                                                                            
  // for finding sticky_soul pairs which satisfy:                               
  // - at least one is flagged as active                                        
  // - | x_i - x_j | < max(size_i,size_j)                                       
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
#define CHECK_PAIR							       \
  Rq = dist_sq(pos(A),pos(B));                      /* R^2 = distance^2      */\
  if(Rq < sizeq(A) || Rq < sizeq(B)) add_pair(A,B); /* |R| < max(s_A,s_B) ?  */
  //////////////////////////////////////////////////////////////////////////////
  class neighbour_finder : public basic_finder {
  protected:
    //--------------------------------------------------------------------------
    void many(bool const&all, soul_iter const&A,
	      soul_iter const&B0, soul_iter const&BN) const
    {
      register real      Rq;
      if(all) for(register soul_iter B=B0; B<BN; B++)
	{ CHECK_PAIR }
      else    for(register soul_iter B=B0; B<BN; B++) if(is_active(B))
	{ CHECK_PAIR }
    }
#undef CHECK_PAIR
    //--------------------------------------------------------------------------
    void single(soul_iter const&A, soul_iter const&B) const
    {
      if(!is_active(A) && !is_active(B)) return;   // no actives -> no job      
      register real Rq = dist_sq(pos(A),pos(B));   // R^2 = distance^2          
      if(Rq < sizeq(A) || Rq < sizeq(B))           // |R| < max(s_A,s_B) ?      
	return add_pair(A,B);
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, soul_iter const&B) const
    {
      register real
	Rq = dist_sq(pos(A),pos(B)),               // R^2 = distance^2          
	x  = max(rmax(A)+size(B),size(A));         // max size distance         
      if(Rq < square(x)) return false;             // |R| < max(s_A,r_A+s_B)    
      return true;
    }
    //--------------------------------------------------------------------------
    bool discard(cell_iter const&A, cell_iter const&B) const
    {
      register real
	Rq = dist_sq(pos(A),pos(B)),               // R^2 = distance^2          
	x  = max(rmax(A)+size(B),rmax(B)+size(A)); // max size distance         
      if(Rq < square(x)) return false;             // |R| < max(rB+sA,rA+sB)    
      return true;
    }
  public:
    //--------------------------------------------------------------------------
    neighbour_finder(const uint n, elem_pair*l) : basic_finder(n,l) {}
    //--------------------------------------------------------------------------
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || size(A) > size(B);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
}                                                  // END: namespace nbdy       
using namespace nbdy;
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// auxiliaries for class nbdy::sticky_tree                                      
//                                                                              
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  typedef sticky_tree::cell_iterator cell_iter;
  typedef sticky_tree::soul_iterator soul_iter;
  //////////////////////////////////////////////////////////////////////////////
  inline void prepare_sticky(const sticky_tree*const&T)
  {
    if(T->begin_souls()==0) return;                // no souls -> no job        
    // 1. update souls with size and vel (body flag & pos known from tree souls)
    if     (T->use_sbodies())
      LoopSouls(sticky_tree,T,Si) Si->set_sticky(T->my_sbodies());
#ifdef falcON_MPI
    else if(T->use_pbodies())
      error("[sticky_tree::prepare_sticky()]: use of pbodies not implemented");
#endif
    else if(T->use_barrays())
      LoopSouls(sticky_tree,T,Si) Si->set_sticky(T->my_barrays());
    else
      falcON_Error("tree has neither bodies nor array data");
    // 2. for all cells: compute Sum x_i and Sum v_i, pass active flag up       
    register vect spos,svel;                       // vel & pos sum             
    LoopCellsUp(sticky_tree,T,Ci) {
      spos = zero;
      svel = zero;
      Ci->reset_flag();
      LoopSoulKids(cell_iter,Ci,s) {
	spos+= pos(s);
	svel+= vel(s);
	Ci->add_active_flag(s);
      }
      LoopCellKids(cell_iter,Ci,c) {
	spos+= pos(c);
	svel+= vel(c);
	Ci->add_active_flag(c);
      }
      Ci->pos() = spos;
      Ci->vel() = svel;
    }
    // 3. for all cells: normalize pos & vel; compute size & vrad               
    register real iN,dmax,vmax;
    LoopCellsUp(sticky_tree,T,Ci) {
      iN = 1./real(number(Ci));
      Ci->pos() *= iN;
      Ci->vel() *= iN;
      dmax = zero;
      vmax = zero;
      LoopSoulKids(cell_iter,Ci,s) {
	update_max(dmax,sqrt(dist_sq(pos(Ci),pos(s))) + size(s));
	update_max(vmax,sqrt(dist_sq(vel(Ci),vel(s))));
      }
      LoopCellKids(cell_iter,Ci,c) {
	update_max(dmax,sqrt(dist_sq(pos(Ci),pos(c))) + size(c));
	update_max(vmax,sqrt(dist_sq(vel(Ci),vel(c))) + vrad(c));
      }
      Ci->size() = dmax;
      Ci->vrad() = vmax;
    }
  }
  //----------------------------------------------------------------------------
  inline void prepare_sph(const sticky_tree*const&T)
  {
    if(T->begin_souls()==0) return;                  // no souls -> no job      
    // 1. update souls with size and vel (body flag & pos known from tree souls)
    if     (T->use_sbodies())
      LoopSouls(sticky_tree,T,Si) Si->set_sph(T->my_sbodies());
#ifdef falcON_MPI
    else if(T->use_pbodies())
      error("[sticky_tree::prepare_sph()]: use of pbodies not implemented");
#endif
    else if(T->use_barrays())
      LoopSouls(sticky_tree,T,Si) Si->set_sph(T->my_barrays());
    else
      falcON_Error("tree has neither bodies nor array data");
    // 2. for all cells: compute Sum x_i, pass active flag up
    register vect spos;
    LoopCellsUp(sticky_tree,T,Ci) {
      spos = zero;
      Ci->reset_flag();
      LoopSoulKids(cell_iter,Ci,s) {
	spos+= pos(s);
	Ci->add_active_flag(s);
      }
      LoopCellKids(cell_iter,Ci,c) {
	spos+= pos(c);
	Ci->add_active_flag(c);
      }
      Ci->pos() = spos;
    }
    // 3. for all cells: normalize pos, compute rmax & size
    register real dmax,smax,d;
    LoopCellsUp(sticky_tree,T,Ci) {
      Ci->pos() /= real(number(Ci));
      dmax = zero;
      smax = zero;
      LoopSoulKids(cell_iter,Ci,s) {
	d = sqrt(dist_sq(pos(Ci),pos(s)));
	update_max(dmax, d + size(s));
	update_max(smax, d);
      }
      LoopCellKids(cell_iter,Ci,c) {
	d = sqrt(dist_sq(pos(Ci),pos(c)));
	update_max(dmax, d + size(c));
	update_max(smax, d + rmax(c));
      }
      Ci->size() = dmax;
      Ci->rmax() = smax;
    }
  }
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class nbdy::sticky_tree                                                      
//                                                                              
////////////////////////////////////////////////////////////////////////////////
void sticky_tree::make_iaction_list(elem_pair *bl,
				    uint const&nl,
				    uint      &na,
				    real const&tau) const
{
  if(tau < zero) {             // SPH
    prepare_sph(this);
    neighbour_finder sfind(nl,bl);
    MutualInteractor<neighbour_finder> MI(&sfind,depth());
    MI.cell_self(root());
    na = sfind.actual_size_of_list();
  } else {                     // sticky particles
    prepare_sticky(this);
     sticky_finder sfind(tau,nl,bl);
   MutualInteractor<sticky_finder> MI(&sfind,depth());
    MI.cell_self(root());
    na = sfind.actual_size_of_list();
  }
}
//------------------------------------------------------------------------------
void sticky_tree::count_neighbours() const
{
  prepare_sph(this);
  neighbour_counter<sticky_tree,individual> scount;
  MutualInteractor<neighbour_counter<sticky_tree,individual> > 
    MI(&scount,depth());
  MI.cell_self(root());
  if     (use_sbodies())
    LoopMySouls(Si) { if(is_active(Si)) my_sbodies()->num(mybody(Si))=num(Si); }
#ifdef falcON_MPI
  else if(use_pbodies())
    falcON_ErrorF("use of pbodies not implemented","sticky_tree::set_num()");
#endif
  else if(use_barrays())
    LoopMySouls(Si) { if(is_active(Si)) my_barrays()->num(mybody(Si))=num(Si); }
  else
    falcON_Error("tree has neither bodies nor array data");
}
////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------+
//                                                                             |
// stic.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2001                                          |
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
#include <body.h>
#include <public/tree.h>

using namespace nbdy;
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class nbdy::basic_finder                                                     
//                                                                              
////////////////////////////////////////////////////////////////////////////////
inline void basic_finder::add_pair(soul_iter const&A, soul_iter const&B) const
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
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class nbdy::sticky_finder                                                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
void sticky_finder::single(soul_iter const&A, soul_iter const&B) const
{
  if(!is_sink(A) && !is_sink(B)) return;             // no sinks -> no job      
  register vect R = pos(A)-pos(B);                   // R   = distance vector   
  register real xq= square(size(A)+size(B));         // x^2 = (combined size)^2 
  if(norm(R) < xq) return add_pair(A,B);             // |R|<x: overlap at t=0   
  if(TAU==zero) return;                              // no times t>0            
  register vect V = vel(A)-vel(B);                   // V   = velocity diff     
  register real RV= R*V;                             // R*V = scalar product    
  if(RV>zero) return;                                // R*V>0: diverging orbits 
  register real t = min(TAU, -RV/norm(V));           // t of minimum distance   
  if(norm(R+t*V) < xq) return add_pair(A,B);         // d_min<x: overlapping    
}
//------------------------------------------------------------------------------
#define CHECK_PAIR							       \
  R = pos(A)-pos(B);                                 /* distance vector     */ \
  xq= square(size(A)+size(B));                       /* (combined size)^2   */ \
  if(norm(R) < xq) { add_pair(A,B); continue; }      /* overlap at t=0      */ \
  if(TAU==zero) continue;                            /* no times t>0        */ \
  V = vel(A)-vel(B);                                 /* velocity difference */ \
  RV= R*V;                                           /* scalar product      */ \
  if(RV>zero) continue;                              /* diverging orbits    */ \
  t = min(TAU,-RV/norm(V));                          /* time of min distance*/ \
  if(norm(R+t*V) < xq) add_pair(A,B);                /* overlap: add to list*/

void sticky_finder::many(bool const&all, soul_iter const&A,
			 soul_iter const&B0, soul_iter const&BN) const
{
  register vect      R,V;
  register real      xq,RV,t;
  if(all) for(register soul_iter B=B0; B<BN; B++)                { CHECK_PAIR }
  else    for(register soul_iter B=B0; B<BN; B++) if(is_sink(B)) { CHECK_PAIR }
}
#undef CHECK_PAIR
//------------------------------------------------------------------------------
bool sticky_finder::discard(cell_iter const&A, soul_iter const&B) const
  // true  if there cannot possibly be an interaction between any A and B       
  // false if an interaction between any A and B is possible                    
{
  register vect R = pos(A)-pos(B);                   // R   = distance vector   
  register real
    Rq = norm(R),                                    // R^2 = distance^2        
    x  = size(A)+size(B);                            // x   = combined size     
  if(Rq < x*x) return false;                         // overlap at t=0          
  if(TAU==zero) return true;                         // no times t>0            
  register vect V = vel(A)-vel(B);                   // V   = velocity diff     
  register real
    w  = vrad(A),                                    // w   = v-size            
    wq = w*w,                                        // w^2 = (v-size)^2        
    RV = R*V,                                        // R*V = scalar product    
    RVq= RV*RV;                                      // (R*V)^2                 
  if(RV>zero && RVq>wq*Rq) return true;              // R*V/|R| > w: diverging  
  register real
    Vq = norm(V),                                    // V^2 = (v-distance)^2    
    t  = (wq>=Vq)? TAU :                             // v-size too large        
    min(TAU, (w*sqrt((Rq*Vq-RVq)/(Vq-wq))-RV)/Vq );  // t of min dist           
  if(norm(R+t*V) < square(x+t*w)) return false;      // possible overlap        
  return true;                                       // no overlap -> discard   
}
//------------------------------------------------------------------------------
bool sticky_finder::discard(cell_iter const&A, cell_iter const&B) const
  // true  if there cannot possibly be an interaction between any A and B       
  // false if an interaction between any A and B is possible                    
{
  register vect R = pos(A)-pos(B);                   // R   = distance vector   
  register real
    Rq = norm(R),                                    // R^2 = distance^2        
    x  = size(A)+size(B);                            // x   = combined size     
  if(Rq < x*x) return false;                         // overlap at t=0          
  if(TAU==zero) return true;                         // no times t>0            
  register vect V = vel(A)-vel(B);                   // V   = velocity diff     
  register real
    w  = vrad(A)+vrad(B),                            // w   = combined v-size   
    wq = w*w,                                        // w^2 = (comb v-size)^2   
    RV = R*V,                                        // R*V = scalar product    
    RVq= RV*RV;                                      // (R*V)^2                 
  if(RV>zero && RVq>wq*Rq) return true;              // R*V/|R| > w: diverging  
  register real
    Vq = norm(V),                                    // V^2 = (v-distance)^2    
    t  = (wq>=Vq)? TAU :                             // v-size too large        
    min(TAU, (w*sqrt((Rq*Vq-RVq)/(Vq-wq))-RV)/Vq );  // t of min dist           
  if(norm(R+t*V) < square(x+t*w)) return false;      // possible overlap        
  return true;                                       // no overlap -> discard   
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class nbdy::neighbour_finder                                                 
//                                                                              
////////////////////////////////////////////////////////////////////////////////
void neighbour_finder::single(soul_iter const&A, soul_iter const&B) const
{
  if(!is_sink(A) && !is_sink(B)) return;             // no sinks -> no job      
  register real Rq = dist_sq(pos(A),pos(B));         // R^2 = distance^2        
  if(Rq < sizeq(A) || Rq < sizeq(B))                 // |R| < max(s_A,s_B) ?    
    return add_pair(A,B);
}
//------------------------------------------------------------------------------
#define CHECK_PAIR							       \
  Rq = dist_sq(pos(A),pos(B));                        /* R^2 = distance^2    */\
  if(Rq < sizeq(A) || Rq < sizeq(B)) add_pair(A,B);   /* |R| < max(s_A,s_B) ?*/

void neighbour_finder::many(bool const&all, soul_iter const&A,
			    soul_iter const&B0, soul_iter const&BN) const
{
  register real      Rq;
  if(all) for(register soul_iter B=B0; B<BN; B++)                { CHECK_PAIR }
  else    for(register soul_iter B=B0; B<BN; B++) if(is_sink(B)) { CHECK_PAIR }
}
#undef CHECK_PAIR
//------------------------------------------------------------------------------
bool neighbour_finder::discard(cell_iter const&A, soul_iter const&B) const
{
  register real
    Rq = dist_sq(pos(A),pos(B)),                     // R^2 = distance^2        
    x  = max(rmax(A)+size(B),size(A));               // max size distance       
  if(Rq < square(x)) return false;                   // |R| < max(s_A,r_A+s_B)  
  return true;
}
//------------------------------------------------------------------------------
bool neighbour_finder::discard(cell_iter const&A, cell_iter const&B) const
{
  register real
    Rq = dist_sq(pos(A),pos(B)),                     // R^2 = distance^2        
    x  = max(rmax(A)+size(B),rmax(B)+size(A));       // max size distance       
  if(Rq < square(x)) return false;                   // |R| < max(rB+sA,rA+sB)  
  return true;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class nbdy::neighbour_counter                                                
//                                                                              
////////////////////////////////////////////////////////////////////////////////
template<typename TREE> void neighbour_counter<TREE,individual>::
single(soul_iter const&A, soul_iter const&B) const {
  if(!is_sink(A) && !is_sink(B)) return;             // no sinks -> no job      
  register real Rq = dist_sq(pos(A),pos(B));         // R^2 = distance^2        
  if(is_sink(A) && Rq < sizeq(A)) A->inc();          // A is sink && |R|<siza(A)
  if(is_sink(B) && Rq < sizeq(B)) B->inc();          // B is sink && |R|<siza(B)
}
//------------------------------------------------------------------------------
template<typename TREE>
void neighbour_counter<TREE,individual>::many(bool      const&all,
					      soul_iter const&A,
					      soul_iter const&B0,
					      soul_iter const&BN) const {
  register real      Rq;
  register soul_iter B;
  if(is_sink(A))
    if(all)
      for(B=B0; B<BN; B++) {
        Rq = dist_sq(pos(A),pos(B));
        if(Rq < sizeq(A)) A->inc();
        if(Rq < sizeq(B)) B->inc();
      }
    else
      for(B=B0; B<BN; B++) {
        Rq = dist_sq(pos(A),pos(B));
        if(              Rq < sizeq(A)) A->inc();
        if(is_sink(B) && Rq < sizeq(B)) B->inc();
      }
  else
    if(all)
      for(B=B0; B<BN; B++) {
        Rq = dist_sq(pos(A),pos(B));
        if(Rq < sizeq(B)) B->inc();
      }
    else
      for(B=B0; B<BN; B++) if(is_sink(B)) {
        Rq = dist_sq(pos(A),pos(B));
        if(Rq < sizeq(B)) B->inc();
      }
}
//------------------------------------------------------------------------------
template<typename TREE> bool neighbour_counter<TREE,individual>::
discard(cell_iter const&A, soul_iter const&B) const
{
  return dist_sq(pos(A),pos(B)) >                    // |xA-xB| >               
    square(max(rmax(A)+size(B),size(A)));            // max(sA,rA+sB)  ?        
}
//------------------------------------------------------------------------------
template<typename TREE> bool neighbour_counter<TREE,individual>::
discard(cell_iter const&A, cell_iter const&B) const
{
  return dist_sq(pos(A),pos(B)) >                    // |xA-xB| >               
    square(max(rmax(A)+size(B),rmax(B)+size(A)));    // max(rb+sA,rA+sB)  ?     
}
//------------------------------------------------------------------------------
template<typename TREE> void neighbour_counter<TREE,global>::
single(soul_iter const&A, soul_iter const&B) const {
  if(!is_sink(A) && !is_sink(B)) return;             // no sinks -> no job      
  if(dist_sq(pos(A),pos(B)) < EPQ) {                 // distance^2 < eps^2      
    if(is_sink(A)) A->inc();                         // A is sink && |R|<siza(A)
    if(is_sink(B)) B->inc();                         // B is sink && |R|<siza(B)
  }
}
//------------------------------------------------------------------------------
template<typename TREE>
void neighbour_counter<TREE,global>::many(bool      const&all,
				      soul_iter const&A,
				      soul_iter const&B0,
				      soul_iter const&BN) const {
  register soul_iter B;
  if(is_sink(A))
    if(all)
      for(B=B0; B<BN; B++) {
	if(dist_sq(pos(A),pos(B)) < EPQ) { A->inc();
	                                   B->inc(); }
      }
    else
      for(B=B0; B<BN; B++) {
	if(dist_sq(pos(A),pos(B)) < EPQ) {                A->inc();
	                                   if(is_sink(B)) B->inc(); }
      }
  else
    if(all)
      for(B=B0; B<BN; B++) {
	if(dist_sq(pos(A),pos(B)) < EPQ) { B->inc(); }
      }
    else
      for(B=B0; B<BN; B++) 
	if(is_sink(B) && dist_sq(pos(A),pos(B)) < EPQ) { B->inc(); }
}
//------------------------------------------------------------------------------
template<typename TREE> bool neighbour_counter<TREE,global>::
discard(cell_iter const&A, soul_iter const&B) const
{
  return dist_sq(pos(A),pos(B)) >                    // |xA-xB| >               
    square(max(rmax(A)+EPS,size(A)));                // max(sA,rA+sB)  ?        
}
//------------------------------------------------------------------------------
template<typename TREE> bool neighbour_counter<TREE,global>::
discard(cell_iter const&A, cell_iter const&B) const
{
  return dist_sq(pos(A),pos(B)) >                    // |xA-xB| >               
    square(max(rmax(A)+size(B),rmax(B)+size(A)));    // max(rb+sA,rA+sB)  ?     
}
//------------------------------------------------------------------------------
#ifdef ALLOW_INDI
template class neighbour_counter<grav_tree,individual>;
#endif
template class neighbour_counter<sticky_tree,individual>;
template class neighbour_counter<grav_tree,global>;
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class nbdy::sticky_tree                                                      
//                                                                              
////////////////////////////////////////////////////////////////////////////////
void sticky_tree::prepare_sticky()
{
  if(begin_souls()==0) return;                     // no souls -> no job        
  // 1. update souls with size and vel (body flag & pos known from tree souls)  
  if     (use_sbodies())
    LoopMySouls(Si) Si->set_sticky(my_sbodies());
#ifdef ALLOW_MPI
  else if(use_pbodies())
    error("[sticky_tree::prepare_sticky()]: use of pbodies not implemented");
#endif
  else 
    LoopMySouls(Si) Si->set_sticky(S,V);
  // 2. for all cells: compute Sum x_i and Sum v_i, pass sink flag up           
  register vect spos,svel;                         // vel & pos sum             
  LoopMyCellsUp(Ci) {
    spos = zero;
    svel = zero;
    Ci->reset_flag();
    LoopSoulKids(cell_iterator,Ci,s) {
      spos+= pos(s);
      svel+= vel(s);
      Ci->add_sink_flag(s);
    }
    LoopCellKids(cell_iterator,Ci,c) {
      spos+= pos(c);
      svel+= vel(c);
      Ci->add_sink_flag(c);
    }
    Ci->pos() = spos;
    Ci->vel() = svel;
  }
  // 3. for all cells: normalize pos & vel; compute size & vrad                 
  register real iN,dmax,vmax;
  LoopMyCellsUp(Ci) {
    iN = 1./real(number(Ci));
    Ci->pos() *= iN;
    Ci->vel() *= iN;
    dmax = zero;
    vmax = zero;
    LoopSoulKids(cell_iterator,Ci,s) {
      update_max(dmax,sqrt(dist_sq(pos(Ci),pos(s))) + size(s));
      update_max(vmax,sqrt(dist_sq(vel(Ci),vel(s))));
    }
    LoopCellKids(cell_iterator,Ci,c) {
      update_max(dmax,sqrt(dist_sq(pos(Ci),pos(c))) + size(c));
      update_max(vmax,sqrt(dist_sq(vel(Ci),vel(c))) + vrad(c));
    }
    Ci->size() = dmax;
    Ci->vrad() = vmax;
  }
}
//------------------------------------------------------------------------------
void sticky_tree::prepare_sph()
{
  if(begin_souls()==0) return;                     // no souls -> no job        
  // 1. update souls with size and vel (body flag & pos known from tree souls)  
  if     (use_sbodies())
    LoopMySouls(Si) Si->set_sph(my_sbodies());
#ifdef ALLOW_MPI
  else if(use_pbodies())
    error("[sticky_tree::prepare_sph()]: use of pbodies not implemented");
#endif
  else 
    LoopMySouls(Si) Si->set_sph(S);
  // 2. for all cells: compute Sum x_i, pass sink flag up
  register vect spos;
  LoopMyCellsUp(Ci) {
    spos = zero;
    Ci->reset_flag();
    LoopSoulKids(cell_iterator,Ci,s) {
      spos+= pos(s);
      Ci->add_sink_flag(s);
    }
    LoopCellKids(cell_iterator,Ci,c) {
      spos+= pos(c);
      Ci->add_sink_flag(c);
    }
    Ci->pos() = spos;
  }
  // 3. for all cells: normalize pos, compute rmax & size
  register real iN,dmax,smax,d;
  LoopMyCellsUp(Ci) {
    iN = 1./real(number(Ci));
    Ci->pos() *= iN;
    dmax = zero;
    smax = zero;
    LoopSoulKids(cell_iterator,Ci,s) {
      d = sqrt(dist_sq(pos(Ci),pos(s)));
      update_max(dmax, d + size(s));
      update_max(smax, d);
    }
    LoopCellKids(cell_iterator,Ci,c) {
      d = sqrt(dist_sq(pos(Ci),pos(c)));
      update_max(dmax, d + size(c));
      update_max(smax, d + rmax(c));
    }
    Ci->size() = dmax;
    Ci->rmax() = smax;
  }
}
//------------------------------------------------------------------------------
void sticky_tree::set_num(grav_soul*S0) const
{
  if     (use_sbodies())
    LoopMySouls(Si) { if(is_sink(Si)) my_sbodies()->num(mybody(Si)) = num(Si); }
#ifdef ALLOW_MPI
  else if(use_pbodies())
    error("[sticky_tree::set_num()]: use of pbodies not implemented");
#endif
  else error("[sticky_tree::set_num()]: bodies/arrays mismatch");
}
//------------------------------------------------------------------------------
void sticky_tree::set_num(int*NUM) const
{
  if(! use_arrays()) error("[sticky_tree::set_num()]: bodies/arrays mismatch");
  LoopMySouls(Si) if(is_sink(Si)) NUM[mybody(Si)] = num(Si);
}
////////////////////////////////////////////////////////////////////////////////

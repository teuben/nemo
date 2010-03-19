// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    src/public/lib/gravity.cc
///
/// \author  Walter Dehnen
///
/// \date    2000-2010
///
/// \brief   implements inc/public/gravity.h
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2010  Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
#include <public/gravity.h>
#include <body.h>
#include <public/interact.h>
#include <public/kernel.h>
#include <numerics.h>

using namespace falcON;
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // auxiliary stuff for class GravMAC                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::InvertZ                                                    //
  //                                                                          //
  // methods for inverting                                                    //
  //                                                                          //
  //         theta^(P+2)/(1-theta)^2 * y^a = 1                                //
  //                                                                          //
  // for  1/theta(y)                                                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class InvertZ {
  private:
    static const unsigned N = 1000, N1=N-1;        // size of tables            
    const unsigned P;                              // expansion order           
    const     real A,hA,sA;                        // parameters                
    real          *Z,*Y;                           // tables                    
    //--------------------------------------------------------------------------
    real z(real y) const {                         // z(y) = 1/theta - 1        
      if(y < Y[ 0]) return std::pow(y,hA);
      if(y > Y[N1]) return std::pow(y,sA);
      return Polev(y,Y,Z,N);
    }    
  public:
    //--------------------------------------------------------------------------
    InvertZ(real     a,                            // I: power a                
	    unsigned p) :                          // I: order P                
      P   ( p ),
      A   ( a ),
      hA  ( half * A ),
      sA  ( A/(P+2.) ),
      Z   ( falcON_NEW(real,N) ),
      Y   ( falcON_NEW(real,N) )
    {
      double _z,iA=1./A,
	zmin = 1.e-4,
	zmax = 1.e4,
	lmin = log(zmin),
	dlz  = (log(zmax)-lmin)/double(N1);
      for(unsigned i=0; i!=N; ++i) {
	_z   = std::exp(lmin+i*dlz);
	Z[i] = _z;
	Y[i] = pow(_z*_z*pow(1+_z,P),iA);
      }
    }
    //--------------------------------------------------------------------------
    ~InvertZ() {
      falcON_DEL_A(Z);
      falcON_DEL_A(Y);
    }
    //--------------------------------------------------------------------------
    real invtheta(real y) const {
      return one + z(y);
    }    
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::GravMAC                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
GravMAC::GravMAC(MAC_type mc,
		 real     t0,
		 unsigned p) :
  MAC  ( mc ),
  P    ( p ),
  TH0  ( min(one, falcON::abs(t0)) ),
  iTH0 ( one/TH0)
{
  switch(MAC) {
  case const_theta:
    IZ = 0;
    break;
  case theta_of_M:
    // th^(p+2)    M  (d-2)/d   th0^(p+2)
    // --------  (---)        = ---------
    // (1-th)^2   M0            (1-th0)^2
    IZ = new InvertZ(third,P);
    break;
  case theta_of_M_ov_r:
    // th^(p+2)    Q  (d-2)/(d-1)   th0^(p+2)               M  
    // --------  (---)            = ---------  with  Q := -----
    // (1-th)^2   Q0                (1-th0)^2             r_max
    IZ = new InvertZ(half,P);
    break;
  case theta_of_M_ov_rq:
    // th^(p+2)    S     th0^(p+2)                M   
    // --------  (---) = ---------  with  S := -------
    // (1-th)^2   S0     (1-th0)^2             r_max^2
    IZ = new InvertZ(one,P);
    break;
  }
}
//------------------------------------------------------------------------------
void GravMAC::reset(MAC_type mc,
		    real     t0, 
		    unsigned p) {
  TH0  = min(one,abs(t0));
  iTH0 = one/TH0;
  if(MAC != mc || P != p) {
    if(IZ) falcON_DEL_O(IZ);
    MAC  = mc;
    P    = p;
    switch(MAC) {
    case const_theta:
      IZ = 0;
      break;
    case theta_of_M:
      IZ = new InvertZ(third,P);
      break;
    case theta_of_M_ov_r:
      IZ = new InvertZ(half,P);
      break;
    case theta_of_M_ov_rq:
      IZ = new InvertZ(one,P);
      break;
    }
  }
}
//------------------------------------------------------------------------------
void GravMAC::reset_theta(real t0)
{
  TH0  = min(one,abs(t0));
  iTH0 = one/TH0;
}
//------------------------------------------------------------------------------
void GravMAC::set_rcrit(const GravEstimator*G) const {
  switch(MAC) {
  case const_theta:
    LoopCellsDown(grav::cell_iter,G->my_tree(),Ci)
      Ci->set_rcrit(iTH0);
    break;
  case theta_of_M: {
    real 
      M0 = mass(G->root()),
      iF = pow(square(1-TH0)/pow(TH0,P+2u), 3u) / M0;
    LoopCellsDown(grav::cell_iter,G->my_tree(),Ci)
      Ci->set_rcrit(IZ->invtheta(mass(Ci)*iF));
  } break;
  case theta_of_M_ov_r: {
    int  i  = 0;
    real Q0 = mass(G->root()) / rmax(G->root());
    real *Q = falcON_NEW(real,G->my_tree()->N_cells());
    LoopCellsDown(grav::cell_iter,G->my_tree(),Ci) {
      Q[i] = mass(Ci)/rmax(Ci);
      if(Q[i] > Q0) Q0 = Q[i];
      ++i;
    }
    real iF = square(square(1-TH0)/pow(TH0,P+2u)) / Q0;
    i = 0;
    LoopCellsDown(grav::cell_iter,G->my_tree(),Ci)
      Ci->set_rcrit(IZ->invtheta(iF*Q[i++]));
    falcON_DEL_A(Q);
  } break;
  case theta_of_M_ov_rq: {
    int  i  = 0;
    real S0 = mass(G->root()) / square(rmax(G->root()));
    real *S = falcON_NEW(real,G->my_tree()->N_cells());
    LoopCellsDown(grav::cell_iter,G->my_tree(),Ci) {
      S[i] = mass(Ci)/square(rmax(Ci));
      if(S[i] > S0) S0 = S[i];
      ++i;
    }
    real iF = square(1-TH0)/pow(TH0,P+2u) / S0;
    i = 0;
    LoopCellsDown(grav::cell_iter,G->my_tree(),Ci)
      Ci->set_rcrit(IZ->invtheta(iF*S[i++]));
    falcON_DEL_A(S);
  } break;
  }
}
//------------------------------------------------------------------------------
GravMAC::~GravMAC()
{
  if(IZ) falcON_DEL_O(IZ);
}
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace grav;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // auxiliary stuff for class GravEstimator                                    
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  inline real bmax(vect const&com, cell_iter const&C)
    // This routines returns the distance from the cell's cofm (com)            
    // to its most distant corner.                                              
  {
    return sqrt( square(radius(C)+abs(com[0]-center(C)[0])) +
		 square(radius(C)+abs(com[1]-center(C)[1])) +
		 square(radius(C)+abs(com[2]-center(C)[2])) );
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class GravIactBase                                                         
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class GravIactBase {
    //--------------------------------------------------------------------------
    // types required in interact.h                                             
    //--------------------------------------------------------------------------
  public:
    typedef grav::leaf_iter leaf_iter;
    typedef grav::cell_iter cell_iter;
    //--------------------------------------------------------------------------
    // static methods                                                           
    //--------------------------------------------------------------------------
    static bool take(grav::leaf_iter const&) { return true; }
    static bool take(grav::cell_iter const&) { return true; }
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
  private:
    unsigned              N_PRE[3], N_POST[3];     // direct sums control       
  protected:
    GravStats* const      STAT;                    // statistics                
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
    bool do_direct_pre (cell_iter const&A, leaf_iter const&) const
    {
      return number(A) < N_PRE[0];
    }
    //--------------------------------------------------------------------------
    bool do_direct_post(cell_iter const&A, leaf_iter const&) const
    {
      return is_twig(A) || number(A) < N_POST[0];
    }
    //--------------------------------------------------------------------------
    bool do_direct_post(cell_iter const&A, cell_iter const&B) const
    {
      return (is_twig(A) && is_twig(B)) ||
	     (number(A) < N_POST[1] && number(B) < N_POST[1]);
    }
    //--------------------------------------------------------------------------
    bool do_direct(cell_iter const&A) const
    { 
      return is_twig(A) || number(A) < N_PRE[2];
    }
    //--------------------------------------------------------------------------
    static bool well_separated(cell_iter const&A, cell_iter const&B, real Rq)
    { 
      return Rq > square(rcrit(A)+rcrit(B));
    }
    //--------------------------------------------------------------------------
    static bool well_separated(cell_iter const&A, leaf_iter const&, real Rq)
    { 
      return Rq > rcrit2(A);
    }
    //--------------------------------------------------------------------------
  protected:
    GravIactBase(GravStats* t, unsigned const nd[4]= Default::direct)
      : STAT ( t )
    {
      N_PRE [0] = nd[0];                           // C-B direct before try     
      N_PRE [1] = 0u;                              // C-C direct before try     
      N_PRE [2] = nd[3];                           // C-S direct before try     
      N_POST[0] = nd[1];                           // C-B direct after fail     
      N_POST[1] = nd[2];                           // C-C direct after fail     
      N_POST[2] = 0u;                              // C-S direct after fail     
    }
    //--------------------------------------------------------------------------
  public:
    bool split_first(cell_iter const&A, cell_iter const&B) const
    {
      return is_twig(B) || (!is_twig(A) && rmax(A) > rmax(B));
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::GravIact                                                     
  //                                                                            
  // This class is at the heart of the algorithm. It serves as INTERACTOR in    
  // the template class MutualInteract<>, defined in interact.h, which          
  // encodes the interaction phase in a most general way.                       
  // class falcON::GravIact has member functions for leaf-leaf,                 
  // leaf-cell, cell-leaf, cell-cell, and cell-self interactions (methods       
  // interact()), as well as a function for the evaluation phase, method        
  // evaluate_grav().                                                           
  //                                                                            
  // NOTE. We organize the cells' Taylor coefficient: at a cell's first         
  // interaction, memory for its coefficients is taken from a pre-allocated     
  // pool and returned to the pool when the cell eventually passes through      
  // evaluation phase.                                                          
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class GravIact : 
    public GravIactBase,
    public GravKern
  {
    GravIact           (GravIact const&);          // not implemented           
    GravIact& operator=(GravIact const&);          // not implemented           
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    GravIact(kern_type k,                           // I: type of kernel       
	     GravStats*t,                           // I: statistics           
	     real      e,                           // I: softening length     
	     unsigned  np,                          // I: initial pool size    
	     bool      s    = Default::soften,      //[I: use individual eps?] 
	     unsigned  const nd[4]= Default::direct)//[I: direct sum control]  
	: GravIactBase  ( t,nd ),
	  GravKern      ( k,e,s,np ) {}
    //--------------------------------------------------------------------------
    // interaction phase                                                        
    //--------------------------------------------------------------------------
    void direct_summation(cell_iter const&A) const { 
      direct(A);
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A) const {
      if(!is_active(A)) return true;               // no interaction -> DONE    
      if(do_direct(A)) {                           // IF(suitable)              
	direct(A);                                 //   perform BB iactions     
	STAT->record_direct_CX(A);                 //   record stats            
#ifdef WRITE_IACT_INFO
	std::cerr<<": direct\n";
#endif
	return true;                               //   DONE                    
      }                                            // ELSE                      
      return false;                                //   must be splitted        
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A, cell_iter const&B) const {
      if(!(is_active(A)||is_active(B)))return true;// no interaction -> DONE    
      vect dX = cofm(A)-cofm(B);                   // compute dX = X_A - X_B    
      real Rq = norm(dX);                          // and dX^2                  
      if(well_separated (A,B,Rq)) {                // IF(well separated)        
	approx(A,B,dX,Rq);                         //   interact                
	STAT->record_approx_CC(A,B);               //   record stats            
#ifdef WRITE_IACT_INFO
	std::cerr<<": approx\n";
#endif
	return true;                               //   DONE                    
      }                                            // ENDIF                     
      if(do_direct_post(A,B)) {                    // IF(suitable)              
	direct(A,B);                               //   perform BB iactions     
	STAT->record_direct_CC(A,B);               //   record stats            
#ifdef WRITE_IACT_INFO
	std::cerr<<": direct\n";
#endif
	return true;                               //   DONE                    
      }                                            // ELSE                      
      return false;                                //   must be splitted        
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A, leaf_iter const&B) const {
      if(!(is_active(A)||is_active(B)))return true;// no interaction -> DONE    
      if(do_direct_pre(A,B)) {                     // IF(suitable)              
	direct(A,B);                               //   perform BB iactions     
	STAT->record_direct_CB(A,B);               //   record statistics       
#ifdef WRITE_IACT_INFO
	std::cerr<<": direct\n";
#endif
	return true;                               //   DONE                    
      }                                            // ENDIF                     
      vect dX = cofm(A)-cofm(B);                   // compute R = x_A-x_B       
      real Rq = norm(dX);                          // compute R^2               
      if(well_separated(A,B,Rq)) {                 // IF(well separated)        
	approx(A,B,dX,Rq);                         //   interact                
	STAT->record_approx_CB(A,B);               //   record statistics       
#ifdef WRITE_IACT_INFO
	std::cerr<<": approx\n";
#endif
	return true;                               //   DONE                    
      }                                            // ENDIF                     
      if(do_direct_post(A,B)) {                    // IF(suitable)              
	direct(A,B);                               //   perform BB iactions     
	STAT->record_direct_CB(A,B);               //   record statistics       
#ifdef WRITE_IACT_INFO
	std::cerr<<": direct\n";
#endif
	return true;                               //   DONE                    
      }                                            // ELSE                      
      return false;                                //   cell must be splitted   
    }
    //--------------------------------------------------------------------------
    bool interact(leaf_iter const&A, cell_iter const&B) const {
      return interact(B,A);
    }
    //--------------------------------------------------------------------------
    void interact(leaf_iter const&A, leaf_iter const&B) const {
      if(!(is_active(A) || is_active(B))) return;  // no interaction -> DONE    
      single(A,B);                                 // perform interaction       
      STAT->record_BB(A,B);                        // record statistics         
    }  
    //--------------------------------------------------------------------------
    void evaluate(cell_iter const&C) const {       // evaluation phase          
      flush_buffers();                             // finish interactions       
      eval_grav(C,TaylorSeries(cofm(C)));          // start recursion           
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::GravIactAll                                                  
  //                                                                            
  // Like GravIact, except that all cells and leafs are assumed active.         
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class GravIactAll : 
    public GravIactBase,
    public GravKernAll
  {
    GravIactAll           (GravIactAll const&);    // not implemented           
    GravIactAll& operator=(GravIactAll const&);    // not implemented           
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    GravIactAll(kern_type    k,                     // I: type of kernel        
		GravStats*   t,                     // I: statistics            
		real         e,                     // I: softening length      
		unsigned     np,                    // I: initial pool size     
		bool         s    =Default::soften, //[I: use individual eps?]  
		unsigned const nd[4]=Default::direct  //[I: direct sum control]
		) :
      GravIactBase ( t,nd ),
      GravKernAll  ( k,e,s,np ) {}
    //--------------------------------------------------------------------------
    // interaction phase                                                        
    //--------------------------------------------------------------------------
    void direct_summation(cell_iter const&A) const { 
      direct(A);
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A) const {
      if(do_direct(A)) {                           // IF(suitable)              
	direct(A);                                 //   perform BB iactions     
	STAT->record_direct_CX(A);                 //   record stats            
#ifdef WRITE_IACT_INFO
	std::cerr<<": direct\n";
#endif
	return true;                               //   DONE                    
      }                                            // ELSE                      
      return false;                                //   must be splitted        
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A, cell_iter const&B) const {
      vect dX = cofm(A)-cofm(B);                   // compute dX = X_A - X_B    
      real Rq = norm(dX);                          // and dX^2                  
      if(well_separated (A,B,Rq)) {                // IF(well separated)        
	approx(A,B,dX,Rq);                         //   interact                
	STAT->record_approx_CC(A,B);               //   record stats            
#ifdef WRITE_IACT_INFO
	std::cerr<<": approx\n";
#endif
	return true;                               //   DONE                    
      }                                            // ENDIF                     
      if(do_direct_post(A,B)) {                    // IF(suitable)              
	direct(A,B);                               //   perform BB iactions     
	STAT->record_direct_CC(A,B);               //   record stats            
#ifdef WRITE_IACT_INFO
	std::cerr<<": direct\n";
#endif
	return true;                               //   DONE                    
      }                                            // ELSE                      
      return false;                                //   SPLIT <                 
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A, leaf_iter const&B) const {
      if(do_direct_pre(A,B)) {                     // IF(suitable)              
	direct(A,B);                               //   perform BB iactions     
	STAT->record_direct_CB(A,B);               //   record statistics       
#ifdef WRITE_IACT_INFO
	std::cerr<<": direct\n";
#endif
	return true;                               //   DONE                    
      }                                            // ENDIF                     
      vect dX = cofm(A)-cofm(B);                   // compute R = x_A-x_B       
      real Rq = norm(dX);                          // compute R^2               
      if(well_separated(A,B,Rq)) {                 // IF(well separated)        
	approx(A,B,dX,Rq);                         //   interact                
	STAT->record_approx_CB(A,B);               //   record statistics       
#ifdef WRITE_IACT_INFO
	std::cerr<<": approx\n";
#endif
	return true;                               //   DONE                    
      }                                            // ENDIF                     
      if(do_direct_post(A,B)) {                    // IF(suitable)              
	direct(A,B);                               //   perform BB iactions     
	STAT->record_direct_CB(A,B);               //   record statistics       
#ifdef WRITE_IACT_INFO
	std::cerr<<": direct\n";
#endif
	return true;                               //   DONE                    
      }                                            // ELSE                      
      return false;                                //   cell must be splitted   
    }
    //--------------------------------------------------------------------------
    bool interact(leaf_iter const&A, cell_iter const&B) const {
      return interact(B,A);
    }
    //--------------------------------------------------------------------------
    void interact(leaf_iter const&A, leaf_iter const&B) const {
      single(A,B);                                 // perform interaction       
      STAT->record_BB(A,B);                        // record statistics         
    }  
    //--------------------------------------------------------------------------
    void evaluate(cell_iter const&C) const {       // evaluation phase          
      flush_buffers();                             // finish interactions       
      eval_grav_all(C,TaylorSeries(cofm(C)));      // start recursion           
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // UpdateLeafs()                                                            //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  unsigned UpdateLeafs(const OctTree*tree, bool i_soft)
  {
    unsigned n=0;
    if(i_soft) {
      CheckMissingBodyData(tree->my_bodies(),
			   fieldset::m|fieldset::e|fieldset::f);
      if(debug(1))
	LoopLeafs(grav::leaf,tree,Li) {
	  Li->copy_from_bodies_mass(tree->my_bodies());
	  Li->copy_from_bodies_eph (tree->my_bodies());
	  Li->copy_from_bodies_flag(tree->my_bodies());
	  if(is_active(Li)) ++n;
	  if(mass(Li) <= zero)
	    falcON_THROW("GravEstimator: mass of body #%d=%f "
			 "but falcON requires positive masses\n",
			 tree->my_bodies()->bodyindex(mybody(Li)),
			 mass(Li));
	} 
      else
	LoopLeafs(grav::leaf,tree,Li) {
	  Li->copy_from_bodies_mass(tree->my_bodies());
	  Li->copy_from_bodies_eph (tree->my_bodies());
	  Li->copy_from_bodies_flag(tree->my_bodies());
	  if(is_active(Li)) ++n;
	} 
    } else {
      CheckMissingBodyData(tree->my_bodies(),fieldset::m|fieldset::f);
      if(debug(1))
	LoopLeafs(grav::leaf,tree,Li) {
	  Li->copy_from_bodies_mass(tree->my_bodies());
	  Li->copy_from_bodies_flag(tree->my_bodies());
	  if(is_active(Li)) ++n;
	  if(mass(Li) <= zero)
	    falcON_THROW("GravEstimator: mass of body #%d=%f "
			 "but falcON requires positive masses\n",
			 tree->my_bodies()->bodyindex(mybody(Li)),
			 mass(Li));
	}
      else
	LoopLeafs(grav::leaf,tree,Li) {
	  Li->copy_from_bodies_mass(tree->my_bodies());
	  Li->copy_from_bodies_flag(tree->my_bodies());
	  if(is_active(Li)) ++n;
	}
    }
    return n;
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // is_act<bool>()                                                           //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<bool> struct __IsAct;
  template<> struct __IsAct<0> {
    template<typename T> static bool is(T const&t) { return is_active(t); }
  };
  template<> struct __IsAct<1> {
    template<typename T> static bool is(T const&) { return 1; }
  };

  template<bool ALL> inline
  bool is_act(const grav::leaf*L) { return __IsAct<ALL>::is(L); }
  template<bool ALL> inline
  bool is_act(flags F) { return __IsAct<ALL>::is(F); }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // UpdateBodiesGrav<ALL_ACT>()                                              //
  //                                                                          //
  // for active OR all leafs:                                                 //
  // - copy pot & acc to their associated bodies                              //
  // - optionally copy eps                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<bool ALL_ACT> 
  void UpdateBodiesGrav(const OctTree*T,
			real          G
#ifdef falcON_ADAP
		      , bool          U
#endif
			)
  {
#ifdef falcON_ADAP
    if(U) {
      CheckMissingBodyData(T->my_bodies(),
			   fieldset::e|fieldset::a|fieldset::p);
      if(G!=one) {
	LoopLeafs(grav::leaf,T,Li) if(is_act<ALL_ACT>(Li)) {
	  Li->copy_to_bodies_eps (T->my_bodies());
	  Li->copy_to_bodies_grav(T->my_bodies(),G);
	}
      } else { 
	LoopLeafs(grav::leaf,T,Li) if(is_act<ALL_ACT>(Li)) {
	  Li->copy_to_bodies_eps (T->my_bodies());
	  Li->copy_to_bodies_grav(T->my_bodies());
	}
      }
    } else {
#endif
      CheckMissingBodyData(T->my_bodies(),fieldset::a|fieldset::p);
      if(G!=one) {
	LoopLeafs(grav::leaf,T,Li) if(is_act<ALL_ACT>(Li))
	  Li->copy_to_bodies_grav(T->my_bodies(),G);
      } else {
	LoopLeafs(grav::leaf,T,Li) if(is_act<ALL_ACT>(Li))
	  Li->copy_to_bodies_grav(T->my_bodies());
      }
#ifdef falcON_ADAP
    }
#endif
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // ResetBodiesGrav<ALL_ACT>()                                               //
  //                                                                          //
  // - reset pot & acc of active OR all bodies to zero                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<bool ALL_ACT> void ResetBodiesGrav(const bodies*B)
  {
    CheckMissingBodyData(B,fieldset::a|fieldset::p);
    LoopAllBodies(B,b)
      if(is_act<ALL_ACT>(flag(b))) {
	b.pot() = zero;
	b.acc() = zero;
      }
  } 
}                                                  // END: unnamed namespace    
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::GravEstimator                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void GravEstimator::update_leafs()
{
  if(TREE==0)
    falcON_Error("GravEstimator: no tree");        // IF no tree, FATAL ERROR   
  if(! TREE->is_used_for_grav() )                  // IF tree not used by grav  
    reset();                                       //   reset allocation & flags
  if( TREE->my_bodies()->srce_data_changed() )     // IF body source are changed
    LEAFS_UPTODATE = 0;                            //   leafs are out of date   
  if(! LEAFS_UPTODATE ) {                          // IF leafs are out of date  
    NLA_needed = UpdateLeafs(TREE, INDI_SOFT);     //   update leafs            
    LEAFS_UPTODATE = 1;                            //   leafs are up to date now
    CELLS_UPTODATE = 0;                            //   but cells are not       
    TREE->my_bodies()->mark_srce_data_read();      //   mark bodies: srce read  
  }                                                // ENDIF                     
}
//------------------------------------------------------------------------------
#ifdef falcON_ADAP
# include <proper/gravity_ind.cc>                 // GravEstimator::adjust_eph()
#endif
//------------------------------------------------------------------------------
GravEstimator::~GravEstimator() {
  if(CELL_SRCE) falcON_DEL_A(CELL_SRCE);
  if(LEAF_ACPN) falcON_DEL_A(LEAF_ACPN);
}
//------------------------------------------------------------------------------
unsigned GravEstimator::pass_up(const GravMAC*MAC,
				bool          REUSE)
{
  // passes up: flag, mass, cofm, rmax[, eph], multipoles; sets rcrit           
  report REPORT("GravEstimator::pass_up_for_approx()");
  int n=0;                                         // counter: active cells     
  if(INDI_SOFT) {                                  // IF(individual eps_i)      
    // 1    with eps_i: pass flag, mass, N*eps/2, cofm, rmax, multipoles
#define MASS_WEIGHTED_SOFTENING 1
    LoopCellsUp(grav::cell_iter,TREE,Ci) {         //   LOOP cells upwards      
      Ci->reset_active_flag();                     //     reset activity flag   
      Ci->reset_sink_flag();                       //     reset sink flag       
      real eh (zero);                              //     reset eps/2           
      real mon(zero);                              //     reset monopole        
      vect com(zero);                              //     reset dipole          
      LoopCellKids(cell_iter,Ci,c) {               //     LOOP sub-cells c      
	eh  += eph (c) *
#ifdef MASS_WEIGHTED_SOFTENING
	    mass(c)                                //       sum up M * eps/2    
#else
	    number(c)
#endif
	    ;
	mon += mass(c);                            //       sum up monopole     
	com += mass(c) * cofm(c);                  //       sum up dipole       
	Ci->add_active_flag(c);                    //       add in activity flag
	Ci->add_sink_flag(c);                      //       add in sink flag    
      }                                            //     END LOOP              
      LoopLeafKids(cell_iter,Ci,l) {               //     LOOP sub-leafs s      
#ifdef MASS_WEIGHTED_SOFTENING
	eh  += mass(l) * eph(l);                   //       sum up M * eps/2    
#else
	eh  += eph (l);                            //       sum up eps/2        
#endif
	mon += mass(l);                            //       sum up monopole     
	com += mass(l) * cofm(l);                  //       sum up dipole       
	Ci->add_active_flag(l);                    //       add in activity flag
	Ci->add_sink_flag(l);                      //       add in sink flag    
      }                                            //     END LOOP              
      if(is_active(Ci)) n++;                       //     count active cells    
      Ci->mass() = mon;                            //     set mass              
      mon        = (mon==zero)? zero:one/mon;      //     1/mass                
      com       *= mon;                            //     cofm = dipole/mass    
#ifdef MASS_WEIGHTED_SOFTENING
      eh        /= mass(Ci);                       //     mean eps/2            
#else
      eh        /= number(Ci);                     //     mean eps/2            
#endif
      Ci->eph()  = eh;                             //     set eps/2             
      Mset P(zero);                                //     reset multipoles      
      real dmax(zero);                             //     reset d_max           
      LoopLeafKids(cell_iter,Ci,l) {               //     LOOP sub-leafs s      
	vect Xi = cofm(l); Xi-= com;               //       distance vector     
	update_max(dmax,norm(Xi));                 //       update d_max^2      
	P.add_body(Xi,mass(l));                    //       add multipoles      
      }                                            //     END LOOP              
      if(has_leaf_kids(Ci)) dmax = sqrt(dmax);     //     d_max due to sub-leafs
      LoopCellKids(cell_iter,Ci,c) {               //     LOOP sub-cells c      
	vect Xi = cofm(c); Xi-= com;               //       distance vector     
	real Xq = norm(Xi);                        //       distance^2          
	real x  = dmax - rmax(c);                  //       auxiliary           
	if(zero>x || Xq>square(x))                 //       IF(d>d_max)         
	  dmax = sqrt(Xq) + rmax(c);               //         set d_max = d     
	P.add_cell(Xi,mass(c),poles(c));           //       add multipoles      
      }                                            //     END LOOP              
      Ci->rmax() = REUSE? dmax :                   //     r_max=d_max           
	           min(dmax,bmax(com,Ci));         //     r_max=min(d_max,b_max)
      Ci->cofm() = com;                            //     set dipole = mass*cofm
      Ci->poles()= P;                              //     assign multipoles     
    }                                              //   END LOOP                
  } else {                                         // ELSE (no individual eps)  
    // 2    without eps_i: pass flag, mass, cofm, rmax, multipoles              
    LoopCellsUp(grav::cell_iter,TREE,Ci) {         //   LOOP cells upwards      
      Ci->reset_active_flag();                     //     reset activity flag   
      Ci->reset_sink_flag();                       //     reset sink flag       
      real mon(zero);                              //     reset monopole        
      vect com(zero);                              //     reset dipole          
      LoopCellKids(cell_iter,Ci,c) {               //     LOOP sub-cells c      
	mon += mass(c);                            //       sum up monopole     
	com += mass(c) * cofm(c);                  //       sum up dipole       
	Ci->add_active_flag(c);                    //       add in activity flag
	Ci->add_sink_flag(c);                      //       add in sink flag    
      }                                            //     END LOOP              
      LoopLeafKids(cell_iter,Ci,l) {               //     LOOP sub-leafs s      
	mon += mass(l);                            //       sum up monopole     
	com += mass(l) * cofm(l);                  //       sum up dipole       
	Ci->add_active_flag(l);                    //       add in activity flag
	Ci->add_sink_flag(l);                      //       add in sink flag    
      }                                            //     END LOOP              
      if(is_active(Ci)) n++;                       //     count active cells    
      Ci->mass() = mon;                            //     set mass              
      mon        = (mon==zero)? zero:one/mon;      //     1/mass                
      com       *= mon;                            //     cofm = dipole/mass    
      Mset P(zero);                                //     reset multipoles      
      real dmax(zero);                             //     reset d_max           
      LoopLeafKids(cell_iter,Ci,l) {               //     LOOP sub-leafs s      
	vect Xi = cofm(l); Xi-= com;               //       distance vector     
	update_max(dmax,norm(Xi));                 //       update d_max^2      
	P.add_body(Xi,mass(l));                    //       add multipoles      
      }                                            //     END LOOP              
      if(has_leaf_kids(Ci)) dmax = sqrt(dmax);     //     d_max due to sub-leafs
      LoopCellKids(cell_iter,Ci,c) {               //     LOOP sub-cells c      
	vect Xi = cofm(c); Xi-= com;               //       distance vector     
	real Xq = norm(Xi);                        //       distance^2          
	real x  = dmax - rmax(c);                  //       auxiliary           
	if(zero>x || Xq>square(x))                 //       IF(d>d_max)         
	  dmax = sqrt(Xq) + rmax(c);               //         set d_max = d     
	P.add_cell(Xi,mass(c),poles(c));           //       add multipoles      
      }                                            //     END LOOP              
      Ci->rmax() = REUSE? dmax :                   //     r_max=d_max           
	           min(dmax,bmax(com,Ci));         //     r_max=min(d_max,b_max)
      Ci->cofm() = com;                            //     set dipole = mass*cofm
      Ci->poles()= P;                              //     assign multipoles     
    }                                              //   END LOOP                
  }                                                // ENDIF                     
  // 3  normalize multipoles                                                    
  for(unsigned i=0; i!=TREE->N_cells(); ++i)       // LOOP cell sources         
    CELL_SRCE[i].normalize_poles();                //   normalize multipoles    
  // 4  set rcrit                                                               
  if(MAC) MAC->set_rcrit(this);                    // set r_crit for all cells  
  return n;                                        // return # active cells     
}
//------------------------------------------------------------------------------
bool GravEstimator::prepare(const GravMAC*MAC,
			    bool          al)
{
  SET_I
  if(al) NLA_needed = TREE->N_leafs();             // all leafs are active      
  if(NLA_needed==0) {
    falcON_Warning("in GravEstimator::prepare(): no body active");
    return 1;
  }
  //  - allocate memory for leaf acc/pot/num properties for active leafs        
  if(NLA!=NLA_needed) {                            // IF #active leafs changed  
    if(LEAF_ACPN) falcON_DEL_A(LEAF_ACPN);         //   delete old allocation   
    NLA = NLA_needed;                              //   # new allocation        
    LEAF_ACPN=falcON_NEW(Leaf::acpn_data,NLA);     //   allocate memory         
  }                                                // ENDIF                     
  const bool all = al || NLA==TREE->N_leafs();     // are all active?           
  Leaf::acpn_data*si=LEAF_ACPN;                    // pter to leafs' acpn data  
  if(all)                                          // IF all leafs              
    LoopLeafs(Leaf,TREE,Li) {                      //   LOOP leafs              
      si->reset();                                 //     reset acpn data       
      Li->set_acpn(si++);                          //     set leaf: acpn        
    }                                              //   END LOOP                
  else                                             // ELSE (only active)        
    LoopLeafs(Leaf,TREE,Li)                        //   LOOP leafs              
      if(is_active(Li)) {                          //     IF(leaf is active)    
	si->reset();                               //       reset acpn data     
	Li->set_acpn(si++);                        //       set leaf: acpn      
      } else                                       //     ELSE                  
	Li->set_acpn(0);                           //       leaf has no acpn    
  // update cell source properties                                              
  if(!CELLS_UPTODATE ||                            // IF sources out-of-date    
     NCT != TREE->N_cells()) {                     //    OR changed in number   
    //    - allocate memory for cell source properties & reset their coeff pter 
    if(NCT     < TREE->N_cells() ||                //   IF #cells too large     
       NCT+NCT > TREE->N_cells()   ) {             //   OR #cells too small     
      if(CELL_SRCE) falcON_DEL_A(CELL_SRCE);       //     delete old allocation 
      NCT = TREE->N_cells();                       //     # new allocation      
      CELL_SRCE=falcON_NEW(Cell::srce_data,NCT);   //     allocate memory       
    }                                              //   ENDIF                   
    Cell::srce_data*ci=CELL_SRCE;                  //   pter to cell's source   
    LoopCellsDown(grav::cell_iter,TREE,Ci) {       //   LOOP cells              
      Ci->set_srce(ci++);                          //     give memory to cell   
      Ci->resetCoeffs();                           //     reset cell: Coeffs    
    }                                              //   END LOOP                
    //    - pass source properties up the tree, count active cells              
    NCA = pass_up(MAC,TREE->is_re_used());         //   pass source data up tree
    if(debug(11)) {
      std::ofstream dump;
      dump.open("/tmp/leafs");
      TREE->dump_leafs<Leaf>(dump);
      dump.open("/tmp/cells");
      TREE->dump_cells<Cell>(dump);
      if(debug(11)) DebugInfo("GravEstimator::prepare(): "
			      "leafs dumped to file \"/tmp/leafs\" "
			      "and cells to file \"/tmp/cells\"\n");
    }
    CELLS_UPTODATE = 1;                            //   update up-to-date flag  
  } else {                                         // ELSE                      
    Cell::srce_data*ci=CELL_SRCE;                  //   pter to cell's source   
    LoopCellsDown(grav::cell_iter,TREE,Ci)         //   LOOP cells              
      Ci->set_srce(ci++);                          //     give memory to cell   
  }                                                // ENDIF                     
  return all;                                      // return all                
}
//------------------------------------------------------------------------------
void GravEstimator::exact(bool       al
#ifdef falcON_ADAP
			 ,real       Nsoft,
			  unsigned   Nref,
			  real       emin,
			  real       efac
#endif
			  )
{
  if(GRAV==zero) {
    falcON_Warning("GravEstimator::exact(): G=0\n");
    if(al) ResetBodiesGrav<1>(TREE->my_bodies());
    else   ResetBodiesGrav<0>(TREE->my_bodies());
    return;
  }
  update_leafs();
#ifdef falcON_ADAP
  adjust_eph(al,Nsoft,emin,EPS,Nref,efac);
#endif
  const bool all = prepare(0,al);
  if(N_active_cells()==0)
    return falcON_Warning("GravEstimator::exact(): nobody active");
  STATS->reset(
#ifdef WRITE_IACTION_INFO
	       TREE
#endif
               );
  if(TREE->my_bodies()->N_bodies(bodytype::sink) && EPSSINK != EPS)
    falcON_Warning("GravEstimator::exact(): will ignore eps_sink\n");
  if(all) {
    GravIactAll K(KERNEL,STATS,EPS,0,INDI_SOFT);
    K.direct_summation(root());
    LoopLeafs(Leaf,TREE,Li) Li->normalize_grav();
  } else {
    GravIact K(KERNEL,STATS,EPS,0,INDI_SOFT);
    K.direct_summation(root());
    LoopLeafs(Leaf,TREE,Li) if(is_active(Li)) Li->normalize_grav();
  }
#ifdef falcON_ADAP
  if(all) UpdateBodiesGrav<1>(TREE,GRAV,INDI_SOFT && Nsoft);
  else    UpdateBodiesGrav<0>(TREE,GRAV,INDI_SOFT && Nsoft);
#else
  if(all) UpdateBodiesGrav<1>(TREE,GRAV);
  else    UpdateBodiesGrav<0>(TREE,GRAV);
#endif
  TREE->mark_grav_usage();
}
//------------------------------------------------------------------------------
void GravEstimator::approx(const GravMAC*GMAC,
			   bool          al
#ifdef falcON_ADAP
			   ,
			   real          Nsoft,
			   unsigned      Nref,
			   real          emin,
			   real          efac
#endif
			   )
{
  if(GRAV==zero) {
    falcON_Warning("[GravEstimator::approx()]: G=0\n");
    if(al) ResetBodiesGrav<1>(TREE->my_bodies());
    else   ResetBodiesGrav<0>(TREE->my_bodies());
    return;
  }
  SET_I
  report REPORT("GravEstimator::approx()");
  // update leafs' source fields (mass, flags, epsh) & count active
  update_leafs();
  SET_T(" time: GravEstimator::update_leafs():  ");
#ifdef falcON_ADAP
  // adjust epsh of leafs
  adjust_eph(al,Nsoft,emin,EPS,Nref,efac);
  SET_T(" time: GravEstimator::adjust_eph():    ");
#endif
  // prepare tree: allocate memory (leafs & cells), pass up source, count active
  const bool all = prepare(GMAC,al);
  if(!all && N_active_cells()==0)
    return falcON_Warning("[GravEstimator::approx()]: nobody active");
  SET_T(" time: GravEstimator::prepare():       ");
  report REPORT2("interaction & evaluation");
  STATS->reset(
#ifdef WRITE_IACTION_INFO
	       TREE
#endif
               );
  Ncsize = 4+(all? TREE->N_cells() : N_active_cells())/16;
  if(all) {                                        // IF all are active         
    GravIactAll GK(KERNEL,STATS,EPS,Ncsize,INDI_SOFT,DIR);
                                                   //   init gravity kernel     
    MutualInteractor<GravIactAll> MI(&GK,TREE->depth()-1);
                                                   //   init mutual interactor  
    LoopCellKids(cell_iter,root(),c1) {            //   LOOP cell kids c1     
      MI.cell_self(c1);                            //     self-iaction c1     
      LoopCellSecd(cell_iter,root(),c1+1,c2)       //     LOOP kids c2>c1     
	MI.cell_cell(c1,c2);                       //       interaction c1,2  
      LoopLeafKids(cell_iter,root(),s2) {          //     LOOP leaf kids s    
	if(is_sink(s2)) GK.reset_eps(EPSSINK);     //       sink: switch eps
	MI.cell_leaf(c1,s2);                       //       interaction c1,s  
	if(is_sink(s2)) GK.reset_eps(EPS);         //       sink: switch back
      }                                            //     END LOOP
      GK.evaluate(c1);                             //     evaluation phase    
    }                                              //   END LOOP              
    LoopLeafKids(cell_iter,root(),s1) {            //   LOOP leaf kids s1     
      LoopLeafSecd(cell_iter,root(),s1+1,s2) {     //     LOOP kids s2>s1     
	if(is_sink(s1) || is_sink(s2))             //       IF either is sink
	  GK.reset_eps(EPSSINK);                   //         switch eps
	GK.interact(s1,s2);                        //       interaction s1,s2 
	if(is_sink(s1) || is_sink(s2))             //       IF either is sink
	  GK.reset_eps(EPS);                       //         switch eps back
      }                                            //     END LOOP
      s1->normalize_grav();                        //     evaluation phase    
    }                                              //   END LOOP              
    Ncoeffs = GK.coeffs_used();                    //   remember # coeffs used  
    Nchunks = GK.chunks_used();                    //   remember # chunks used  
  } else {                                         // ELSE: not all are active  
    GravIact GK(KERNEL,STATS,EPS,Ncsize,INDI_SOFT,DIR);
                                                   //   init gravity kernel     
    MutualInteractor<GravIact> MI(&GK,TREE->depth()-1);
                                                   //   init mutual interactor  
    LoopCellKids(cell_iter,root(),c1) {            //   LOOP cell kids c1     
      if(is_active(c1)) {                          //     IF active c1:        
	MI.cell_self(c1);                          //       self-iaction c1     
	LoopCellSecd(cell_iter,root(),c1+1,c2)     //       LOOP kids c2>c1     
	  MI.cell_cell(c1,c2);                     //         interaction c1,2  
	LoopLeafKids(cell_iter,root(),s2) {        //       LOOP leaf kids s
	  if(is_sink(s2)) GK.reset_eps(EPSSINK);   //         sink: switch eps
	  MI.cell_leaf(c1,s2);                     //         interaction c1,s
	  if(is_sink(s2)) GK.reset_eps(EPS);       //         sink: switch back
	}                                          //       END LOOP
	GK.evaluate(c1);                           //       evaluation phase    
      } else {                                     //     ELSE: inactive c1    
	LoopCellSecd(cell_iter,root(),c1+1,c2)     //       LOOP kids c2>c1    
	  if(is_active(c2))MI.cell_cell(c1,c2);    //         interaction c1,2 
	LoopLeafKids(cell_iter,root(),s2)          //       LOOP leaf kids s
	  if(is_active(s2)) {
	    if(is_sink(s2)) GK.reset_eps(EPSSINK); //         sink: switch eps
	    MI.cell_leaf(c1,s2);                   //         interaction c1,s
	    if(is_sink(s2)) GK.reset_eps(EPS);     //         sink: switch back
	  }                                        //       END LOOP
      }                                            //     ENDIF                
    }                                              //   END LOOP              
    LoopLeafKids(cell_iter,root(),s1) {            //   LOOP leaf kids s1     
      if(is_active(s1)) {                          //     IF active s1:        
	LoopLeafSecd(cell_iter,root(),s1+1,s2) {   //       LOOP kids s2>s1     
	  if(is_sink(s1) || is_sink(s2))           //         IF either is sink
	    GK.reset_eps(EPSSINK);                 //           switch eps
	  GK.interact(s1,s2);                      //         interaction s1,2  
	  if(is_sink(s1) || is_sink(s2))           //         IF either is sink
	    GK.reset_eps(EPS);                     //           switch eps back
	}                                          //       END LOOP
	s1->normalize_grav();                      //       evaluation phase    
      } else {                                     //     ELSE: inactive s1    
	LoopLeafSecd(cell_iter,root(),s1+1,s2)     //       LOOP kids s2>s1     
	  if(is_active(s2)) {                      //         active s2 only
	    if(is_sink(s1) || is_sink(s2))         //         IF either is sink
	      GK.reset_eps(EPSSINK);               //           switch eps
	    GK.interact(s1,s2);                    //         interaction s1,2  
	    if(is_sink(s1) || is_sink(s2))         //         IF either is sink
	      GK.reset_eps(EPS);                   //           switch eps back
	  }                                        //       END LOOP
      }                                            //     ENDIF                
    }                                              //   END LOOP              
    Ncoeffs = GK.coeffs_used();                    //   remember # coeffs used  
    Nchunks = GK.chunks_used();                    //   remember # chunks used  
  }                                                // ENDIF                     
  SET_T(" time: interaction & evaluation:        ");
#ifdef falcON_ADAP
  if(all) UpdateBodiesGrav<1>(TREE,GRAV,INDI_SOFT && Nsoft);
  else    UpdateBodiesGrav<0>(TREE,GRAV,INDI_SOFT && Nsoft);
#else
  if(all) UpdateBodiesGrav<1>(TREE,GRAV);
  else    UpdateBodiesGrav<0>(TREE,GRAV);
#endif
  TREE->mark_grav_usage();
  SET_T(" time: updating bodies gravity:         ");
}
//------------------------------------------------------------------------------
namespace {
  using namespace falcON;
  //============================================================================
  unsigned NX;
  real pdim(real const&x) { return cube(x); }
  //----------------------------------------------------------------------------
  struct number_density {
    static real dens(grav::cell_iter const&C)
    { return number(C)/(Nsub*pdim(radius(C))); }
  };
  //----------------------------------------------------------------------------
  struct surface_density {
    static real dens(grav::cell_iter const&C)
    { return mass(C)/(4*square(radius(C))); }
  };
  //----------------------------------------------------------------------------
  struct mass_density {
    static real dens(grav::cell_iter const&C)
    { return mass(C)/(Nsub*pdim(radius(C))); }
  };
  //============================================================================
  template<typename, bool> class guess;
  //----------------------------------------------------------------------------
  template<typename density> class guess<density,1> {
  public:
    static void do_it(cell_iter const&C, real d) {
      if(number(C)>NX || d==zero) d = density::dens(C);
      LoopLeafKids(grav::cell_iter,C,l) l->rho() = d;
      LoopCellKids(grav::cell_iter,C,c) do_it(c,d);
    }
  };
  //----------------------------------------------------------------------------
  template<typename density> class guess<density,0> {
  public:
    static void do_it(cell_iter const&C, real d) {
      if(number(C)>NX || d==zero) d = density::dens(C);
      LoopLeafKids(grav::cell_iter,C,l)
	if(is_active(l)) l->rho() = d;
      LoopCellKids(grav::cell_iter,C,c)
	if     (al_active(c)) guess<density,1>::do_it(c,d);
	else if(is_active(c))                   do_it(c,d);
    }
  };
  //----------------------------------------------------------------------------
  void UpdateBodiesRho(const OctTree*T,
		       bool          all)
  {
    if(all)
      LoopLeafs(grav::leaf,T,Li)
	Li->copy_to_bodies_rho(T->my_bodies());
    else
      LoopLeafs(grav::leaf,T,Li) if(is_active(Li))
	Li->copy_to_bodies_rho(T->my_bodies());
  }
}                                                  // END: unnamed namespace    
//------------------------------------------------------------------------------
void GravEstimator::estimate_nd(bool al, unsigned Nx) const
{
  NX = Nx;
  if(al) guess<number_density,1>::do_it(root(),zero);
  else   guess<number_density,0>::do_it(root(),zero);
  UpdateBodiesRho(TREE,al);
}
//------------------------------------------------------------------------------
void GravEstimator::estimate_sd(bool al, unsigned Nx)
{
  update_leafs();
  prepare(0,al);
  NX = Nx;
  if(al) guess<surface_density,1>::do_it(root(),zero);
  else   guess<surface_density,0>::do_it(root(),zero);
  UpdateBodiesRho(TREE,al);
  TREE->mark_grav_usage();
}
//------------------------------------------------------------------------------
void GravEstimator::estimate_md(bool al, unsigned Nx)
{
  update_leafs();
  prepare(0,al);
  NX = Nx;
  if(al) guess<mass_density,1>::do_it(root(),zero);
  else   guess<mass_density,0>::do_it(root(),zero);
  UpdateBodiesRho(TREE,al);
  TREE->mark_grav_usage();
}
//------------------------------------------------------------------------------
void GravEstimator::dump_cells(std::ostream&o) const
{
  if(CELL_SRCE) TREE->dump_cells<Cell>(o);
  else          TREE->dump_cells<OctTree::Cell>(o);
}
//------------------------------------------------------------------------------
void GravEstimator::dump_leafs(std::ostream&o) const
{
  TREE->dump_leafs<Leaf>(o);
}
////////////////////////////////////////////////////////////////////////////////

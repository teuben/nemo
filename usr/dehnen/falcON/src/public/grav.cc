//-----------------------------------------------------------------------------+
//                                                                             |
// grav.cc                                                                     |
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
// implementing nbdy/grav.h                                                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/grav.h>
#include <body.h>
#ifdef falcON_MPI
#  include <walter/pody.h>
#endif
#include <public/iact.h>
#include <public/kern.h>
#include <public/Pi.h>
#include <public/nums.h>

using namespace nbdy;
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // auxiliary stuff for class grav_mac                                         
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  //===========================================================================#
  // class nbdy::InvertZ                                                       |
  //                                                                           |
  // methods for inverting                                                     |
  //                                                                           |
  //         theta^(P+2)/(1-theta)^2 * y^a = 1                                 |
  //                                                                           |
  // for  1/theta(y)                                                           |
  //                                                                           |
  //===========================================================================#
  class InvertZ {
  private:
    static const uint N = 1000, N1=N-1;            // size of tables            
    const    unsigned P;                           // expansion order           
    const        real A,hA,sA;                     // parameters                
    real             *Z,*Y;                        // tables                    
    //--------------------------------------------------------------------------
    real z(const real y) const {                   // z(y) = 1/theta - 1        
      if(y < Y[ 0]) return std::pow(y,hA);
      if(y > Y[N1]) return std::pow(y,sA);
      return polev(y,Y,Z,N);
    }    
  public:
    //--------------------------------------------------------------------------
    InvertZ(real const&a,                          // I: power a                
	    uint const&p) :                        // I: order P                
      P   ( p ),
      A   ( a ),
      hA  ( half * A ),
      sA  ( A/(P+2.) ),
      Z   ( falcON_New(real,N) ),
      Y   ( falcON_New(real,N) )
    {
      register double z,iA=1./A,
	zmin = 1.e-4,
	zmax = 1.e4,
	lmin = log(zmin),
	dlz  = (log(zmax)-lmin)/double(N1);
      for(register int i=0; i!=N; ++i) {
	z    = std::exp(lmin+i*dlz);
	Z[i] = z;
	Y[i] = pow(z*z*pow(1+z,P),iA);
      }
    }
    //--------------------------------------------------------------------------
    ~InvertZ() {
      delete[] Z;
      delete[] Y;
    }
    //--------------------------------------------------------------------------
    real invtheta(const real y) const {
      return one + z(y);
    }    
  };
  //============================================================================
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class nbdy::grav_mac                                                         
//                                                                              
////////////////////////////////////////////////////////////////////////////////
grav_mac::grav_mac(const MAC_type mc,
		   const real     t0,
		   const uint     p) :
  MAC  ( mc ),
  P    ( p ),
  TH0  ( min(one,abs(t0)) ),
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
#if falcON_NDIM==2
    IZ = 0;
#else
    IZ = falcON_Memory(new InvertZ(third,P));
#endif
    break;
  case theta_of_M_ov_r:
    // th^(p+2)    Q  (d-2)/(d-1)   th0^(p+2)               M  
    // --------  (---)            = ---------  with  Q := -----
    // (1-th)^2   Q0                (1-th0)^2             r_max
#if falcON_NDIM==2
    IZ = 0;
#else
    IZ = falcON_Memory(new InvertZ(half,P));
#endif
    break;
  case theta_of_M_ov_rq:
    // th^(p+2)    S     th0^(p+2)                M   
    // --------  (---) = ---------  with  S := -------
    // (1-th)^2   S0     (1-th0)^2             r_max^2
    IZ = falcON_Memory(new InvertZ(one,P));
    break;
  }
}
//------------------------------------------------------------------------------
void grav_mac::reset(const MAC_type mc,
		     const real     t0, 
		     const uint     p) {
  TH0  = min(one,abs(t0));
  iTH0 = one/TH0;
  if(MAC != mc || P != p) {
    if(IZ) delete IZ;
    MAC  = mc;
    P    = p;
    switch(MAC) {
    case const_theta:
      IZ = 0;
      break;
    case theta_of_M:
#if falcON_NDIM==2
      IZ = 0;
#else
      IZ = falcON_Memory(new InvertZ(third,P));
#endif
      break;
    case theta_of_M_ov_r:
#if falcON_NDIM==2
      IZ = 0;
#else
      IZ = falcON_Memory(new InvertZ(half,P));
#endif
      break;
    case theta_of_M_ov_rq:
      IZ = falcON_Memory(new InvertZ(one,P));
      break;
    }
  }
}
//------------------------------------------------------------------------------
void grav_mac::reset_theta(const real t0)
{
  TH0  = min(one,abs(t0));
  iTH0 = one/TH0;
}
//------------------------------------------------------------------------------
void grav_mac::set_rcrit(const grav_estimator*G) const {
  switch(MAC) {
  case const_theta:
    LoopCellsDown(grav::cell_iter,G->my_tree(),Ci)
      Ci->set_rcrit(iTH0);
    break;
  case theta_of_M: {
#if falcON_NDIM==2
    LoopCellsDown(grav::cell_iter,G->my_tree(),Ci)
      Ci->set_rcrit(iTH0);
#else
    register real 
      M0 = mass(G->root()),
      iF = pow(square(1-TH0)/pow(TH0,P+2u), 3u) / M0;
    LoopCellsDown(grav::cell_iter,G->my_tree(),Ci)
      Ci->set_rcrit(IZ->invtheta(mass(Ci)*iF));
#endif
  } break;
  case theta_of_M_ov_r: {
#if falcON_NDIM==2
    LoopCellsDown(grav::cell_iter,G->my_tree(),Ci)
      Ci->set_rcrit(iTH0);
#else
    register int  i  = 0;
    register real Q0 = mass(G->root()) / rmax(G->root());
    register real *Q = falcON_New(real,G->my_tree()->N_cells());
    LoopCellsDown(grav::cell_iter,G->my_tree(),Ci) {
      Q[i] = mass(Ci)/rmax(Ci);
      if(Q[i] > Q0) Q0 = Q[i];
      ++i;
    }
    register real iF = square(square(1-TH0)/pow(TH0,P+2u)) / Q0;
    i = 0;
    LoopCellsDown(grav::cell_iter,G->my_tree(),Ci)
      Ci->set_rcrit(IZ->invtheta(iF*Q[i++]));
    delete[] Q;
#endif
  } break;
  case theta_of_M_ov_rq: {
    register int  i  = 0;
    register real S0 = mass(G->root()) / square(rmax(G->root()));
    register real *S = falcON_New(real,G->my_tree()->N_cells());
    LoopCellsDown(grav::cell_iter,G->my_tree(),Ci) {
      S[i] = mass(Ci)/square(rmax(Ci));
      if(S[i] > S0) S0 = S[i];
      ++i;
    }
    register real iF = square(1-TH0)/pow(TH0,P+2u) / S0;
    i = 0;
    LoopCellsDown(grav::cell_iter,G->my_tree(),Ci)
      Ci->set_rcrit(IZ->invtheta(iF*S[i++]));
    delete[] S;
  } break;
  }
}
//------------------------------------------------------------------------------
grav_mac::~grav_mac()
{
  if(IZ) delete IZ;
}
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace nbdy; using nbdy::uint;
  using namespace grav;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // auxiliary stuff for class grav_estimator                                   
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  inline real bmax(vect const&com, cell_iter const&C)
    // This routines returns the distance from the cell's cofm (com)            
    // to its most distant corner.                                              
  {
    register real bq=zero;                         // initialize return variable
    for(register int d=0; d!=Ndim; ++d)            // loop dimensions           
      bq += square(radius(C)+abs(com[d]-center(C)[d])); // add in square        
    return sqrt(bq);                               // return sqrt of total      
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class grav_iact_base                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_iact_base {
    //--------------------------------------------------------------------------
    // types required in iact.h                                                 
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
    int                   N_PRE[3], N_POST[3];     // direct sums control       
  protected:
    grav_stat* const      STAT;                    // statistics                
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
    bool do_direct_pre (cell_iter const&A, leaf_iter const&B) const
    {
      return number(A) < N_PRE[0];
    }
    //--------------------------------------------------------------------------
    bool do_direct_post(cell_iter const&A, leaf_iter const&B) const
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
    static bool well_separated(cell_iter const&A, cell_iter const&B,
			       real      const&Rq)
    { 
      return Rq > square(rcrit(A)+rcrit(B));
    }
    //--------------------------------------------------------------------------
    static bool well_separated(cell_iter const&A, leaf_iter const&B,
			       real      const&Rq)
    { 
      return Rq > rcrit2(A);
    }
    //--------------------------------------------------------------------------
  protected:
    grav_iact_base(
		   grav_stat*const&t,                 // I: statistics          
		   int const nd[4]= Default::direct): //[I: direct sum control] 
      STAT       ( t )
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
  //                                                                          //
  // class nbdy::grav_iact                                                    //
  //                                                                          //
  // This class is at the heart of the algorithm. It serves as INTERACTOR in  //
  // the template class MutualInteract<>, defined in iact.h, which encodes    //
  // the interaction phase in a most general way.                             //
  // class nbdy::grav_iact has member functions for leaf-leaf, leaf-cell,     //
  // cell-leaf, cell-cell, and cell-self interactions (methods interact()),   //
  // as well as a function for the evaluation phase, method evaluate_grav().  //
  //                                                                          //
  // NOTE. We organize the cells' Taylor coefficient: at a cell's first       //
  // interaction, memory for its coefficients is taken from a pre-allocated   //
  // pool and returned to the pool when the cell eventually passes through    //
  // evaluation phase.                                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_iact : 
    public grav_iact_base,
    public grav_kern
  {
    grav_iact           (grav_iact const&);        // not implemented           
    grav_iact& operator=(grav_iact const&);        // not implemented           
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    grav_iact(kern_type const&k,                      // I: type of kernel      
	      grav_stat*const&t,                      // I: statistics          
	      real      const&e,                      // I: softening length    
	      uint      const&np,                     // I: initial pool size   
#ifdef falcON_INDI
	      bool      const&s    = Default::soften, //[I: use individual eps?]
#endif
	      int       const nd[4]= Default::direct, //[I: direct sum control] 
	      bool      const&fp   = false) :         //[I: use Pth pole in pot]
      grav_iact_base  ( t,nd ),
      grav_kern       ( k,e,
#ifdef falcON_INDI
			s,
#endif
			np ) {}
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
      register vect dX = cofm(A)-cofm(B);          // compute dX = X_A - X_B    
      register real Rq = norm(dX);                 // and dX^2                  
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
      register vect dX = cofm(A)-cofm(B);          // compute R = x_A-x_B       
      register real Rq = norm(dX);                 // compute R^2               
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
      STAT->record_BB();                           // record statistics         
    }  
    //--------------------------------------------------------------------------
    void evaluate(cell_iter const&C) const {       // evaluation phase          
      flush_buffers();                             // finish interactions       
      eval_grav(C,TaylorSeries(cofm(C)));          // start recursion           
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::grav_iact_all                                                //
  //                                                                          //
  // Like grav_iact, except that all cells and leafs are assumed active.      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_iact_all : 
    public grav_iact_base,
    public grav_kern_all
  {
    grav_iact_all           (grav_iact_all const&);// not implemented           
    grav_iact_all& operator=(grav_iact_all const&);// not implemented           
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    grav_iact_all(kern_type const&k,                  // I: type of kernel      
		  grav_stat*const&t,                  // I: statistics          
		  real      const&e,                  // I: softening length    
		  uint      const&np,                 // I: initial pool size   
#ifdef falcON_INDI
		  bool      const&s =Default::soften, //[I: use individual eps?]
#endif
		  int const    nd[4]=Default::direct  //[I: direct sum control] 
		  ) :
      grav_iact_base ( t,nd ),
      grav_kern_all  ( k, e,
#ifdef falcON_INDI
		       s,
#endif
		       np ) {}
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
      register vect dX = cofm(A)-cofm(B);          // compute dX = X_A - X_B    
      register real Rq = norm(dX);                 // and dX^2                  
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
      register vect dX = cofm(A)-cofm(B);          // compute R = x_A-x_B       
      register real Rq = norm(dX);                 // compute R^2               
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
      STAT->record_BB();                           // record statistics         
    }  
    //--------------------------------------------------------------------------
    void evaluate(cell_iter const&C) const {       // evaluation phase          
      flush_buffers();                             // finish interactions       
      eval_grav_all(C,TaylorSeries(cofm(C)));      // start recursion           
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class UpdateLeafs                                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  struct UpdateLeafs {
    const oct_tree*tree;
#ifdef falcON_INDI
    const bool     i_soft;
    UpdateLeafs(const oct_tree*const&t,
		bool           const&i) : tree(t), i_soft(i) {}
#else
    UpdateLeafs(const oct_tree*const&t) : tree(t) {}
#endif
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    uint operator() (const bodies_type*const&B) const {
      register uint n=0;
#ifdef falcON_INDI
      if(i_soft) {
	LoopLeafs(grav::leaf_type,tree,Li) {
	  Li->copy_from_bodies_mass(B);
	  Li->copy_from_bodies_eph (B);
	  Li->copy_from_bodies_flag(B);
	  if(is_active(Li)) ++n;
	} 
      } else
#endif
	LoopLeafs(grav::leaf_type,tree,Li) {
	  Li->copy_from_bodies_mass(B);
	  Li->copy_from_bodies_flag(B);
	  if(is_active(Li)) ++n;
	}
      return n;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class UseAct<>                                                           //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<bool> struct IsAct;
  template<> struct IsAct<0> {
    static bool is_act(const grav_leaf*const&L) { return is_active(L); }
    static bool is_act(flag            const&F) { return is_active(F); }
  };
  template<> struct IsAct<1> {
    static bool is_act(const grav_leaf*const&L) { return 1; }
    static bool is_act(flag            const&F) { return 1; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class UpdateBodiesGrav                                                   //
  //                                                                          //
  // for active OR all leafs:                                                 //
  // - copy pot & acc to their associated bodies                              //
  // - optionally copy eps                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<bool ALL_ACT> struct UpdateBodiesGrav : 
    protected IsAct<ALL_ACT>
  {
    const real     G;
    const oct_tree*T;
#ifdef falcON_ADAP
    const bool U;                                  // also update eps?          
    UpdateBodiesGrav(const oct_tree*const&t,
		     real           const&g,
		     bool           const&u) : T(t), G(g), U(u) {}
#else
    UpdateBodiesGrav(const oct_tree*const&t,
		     real           const&g) : T(t), G(g) {}
#endif
    template<typename bodies_type>
    uint operator()(const bodies_type*const&B) const {
#ifdef falcON_ADAP
      if(U) {
	if(G!=one) {
	  LoopLeafs(grav_leaf,T,Li) if(is_act(Li)) {
	    Li->copy_to_bodies_eps (B);
	    Li->copy_to_bodies_grav(B,G);
	  }
	} else { 
	  LoopLeafs(grav_leaf,T,Li) if(is_act(Li)) {
	    Li->copy_to_bodies_eps (B);
	    Li->copy_to_bodies_grav(B);
	  }
	}
      } else
#endif
	if(G!=one) {
	  LoopLeafs(grav_leaf,T,Li) if(is_act(Li))
	    Li->copy_to_bodies_grav(B,G);
	} else {
	  LoopLeafs(grav_leaf,T,Li) if(is_act(Li))
	    Li->copy_to_bodies_grav(B);
	}
      return 0;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class ResetBodiesGrav                                                    //
  //                                                                          //
  // - reset pot & acc of active OR all bodies to zero                        //
  //   NOTE: we cannot assume leaf flags to be up to date!                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<bool ALL_ACT> struct ResetBodiesGrav : 
    protected IsAct<ALL_ACT>
  {
    template<typename bodies_type>
    uint operator()(const bodies_type*const&B) const {
      for(register int b=0; b!=B->N_bodies(); ++b) 
	if(is_act(B->flg(b))) {
	  B->pot(b) = zero;
	  B->acc(b) = zero;
	}
      return 0;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class MarkBodiesRead                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  struct MarkBodiesRead {
    template<typename bodies_type>
    uint operator() (const bodies_type*const&B) const {
      B->mark_srce_data_read();
      return 0;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class CheckBodiesChanged                                                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  struct CheckBodiesChanged {
    template<typename bodies_type>
    uint operator() (const bodies_type*const&B) const {
      return B->srce_data_changed();
    }
  };
}                                                  // END: unnamed namespace    
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::grav_estimator                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_INDI
#  define __I_SOFT INDI_SOFT
#else
#  define __I_SOFT
#endif
void grav_estimator::update_leafs()
{
  if(TREE==0) error("grav_estimator: no tree");    // IF no tree, FATAL ERROR   
  if(! TREE->is_used_for_grav() )                  // IF tree not used by grav  
    reset();                                       //   reset allocation & flags
  if( TREE->UseBodies(CheckBodiesChanged()) )      // IF body source are changed
    LEAFS_UPTODATE = 0;                            //   leafs are out of date   
  if(! LEAFS_UPTODATE ) {                          // IF leafs are out of date  
    NLA_needed = TREE->UseBodies(UpdateLeafs(TREE,__I_SOFT)); // update leafs   
    LEAFS_UPTODATE = 1;                            //   leafs are up to date now
    CELLS_UPTODATE = 0;                            //   but cells are not       
    TREE->UseBodies(MarkBodiesRead());             //   mark bodies: srce read  
  }                                                // ENDIF                     
}
//------------------------------------------------------------------------------
#ifdef falcON_ADAP
# include <proper/grav_ind.cc>                // grav_estimator::adjust_eph()   
#endif
//------------------------------------------------------------------------------
nbdy::uint grav_estimator::pass_up(const grav_mac*const&MAC,
				   bool           const&REUSE)
{
  // passes up: flag, mass, cofm, rmax[, eph], multipoles; sets rcrit           
  report REPORT("grav_estimator::pass_up_for_approx()");
  register int n=0;                                // counter: active cells     
#ifdef falcON_INDI
  if(INDI_SOFT) {                                  // IF(individual eps_i)      
    // 1    with eps_i: pass flag, mass, N*eps/2, cofm, rmax, multipoles        
    LoopCellsUp(grav::cell_iter,TREE,Ci) {         //   LOOP cells upwards      
      Ci->reset_active_flag();                     //     reset activity flag   
      register real eh (zero);                     //     reset eps/2           
      register real mon(zero);                     //     reset monopole        
      register vect com(zero);                     //     reset dipole          
      LoopCellKids(cell_iter,Ci,c) {               //     LOOP sub-cells c      
	eh  += eph (c) * number(c);                //       sum up N * eps/2    
	mon += mass(c);                            //       sum up monopole     
	com.add_times(cofm(c),mass(c));            //       sum up dipole       
	Ci->add_active_flag(c);                    //       add in activity flag
      }                                            //     END LOOP              
      LoopLeafKids(cell_iter,Ci,l) {               //     LOOP sub-leafs s      
	eh  += eph (l);                            //       sum up eps/2        
	mon += mass(l);                            //       sum up monopole     
	com.add_times(cofm(l),mass(l));            //       sum up dipole       
	Ci->add_active_flag(l);                    //       add in activity flag
      }                                            //     END LOOP              
      if(is_active(Ci)) n++;                       //     count active cells    
      Ci->mass() = mon;                            //     set mass              
      mon        = (mon==zero)? zero:one/mon;      //     1/mass                
      com       *= mon;                            //     cofm = dipole/mass    
      eh        /= number(Ci);                     //     mean eps/2            
      Ci->eph()  = eh;                             //     set eps/2             
      Mset P(zero);                                //     reset multipoles      
      register real dmax(zero);                    //     reset d_max           
      LoopLeafKids(cell_iter,Ci,l) {               //     LOOP sub-leafs s      
	register vect Xi = cofm(l); Xi-= com;      //       distance vector     
	update_max(dmax,norm(Xi));                 //       update d_max^2      
	P.add_body(Xi,mass(l));                    //       add multipoles      
      }                                            //     END LOOP              
      if(has_leaf_kids(Ci)) dmax = sqrt(dmax);     //     d_max due to sub-leafs
      LoopCellKids(cell_iter,Ci,c) {               //     LOOP sub-cells c      
	register vect Xi = cofm(c); Xi-= com;      //       distance vector     
	register real Xq = norm(Xi);               //       distance^2          
	register real x  = dmax - rmax(c);         //       auxiliary           
	if(zero>x || Xq>square(x))                 //       IF(d>d_max)         
	  dmax = sqrt(Xq) + rmax(c);               //         set d_max = d     
	P.add_cell(Xi,mass(c),poles(c));           //       add multipoles      
      }                                            //     END LOOP              
      Ci->rmax() = REUSE? dmax :                   //     r_max=d_max           
	           min(dmax,bmax(com,Ci));         //     r_max=min(d_max,b_max)
      Ci->cofm() = com;                            //     set dipole = mass*cofm
      Ci->poles()= P;                              //     assign multipoles     
    }                                              //   END LOOP                
  } else                                           // ELSE (no individual eps)  
#endif
  { 
    // 2    without eps_i: pass flag, mass, cofm, rmax, multipoles              
    LoopCellsUp(grav::cell_iter,TREE,Ci) {         //   LOOP cells upwards      
      Ci->reset_active_flag();                     //     reset activity flag   
      register real mon(zero);                     //     reset monopole        
      register vect com(zero);                     //     reset dipole          
      LoopCellKids(cell_iter,Ci,c) {               //     LOOP sub-cells c      
	mon += mass(c);                            //       sum up monopole     
	com.add_times(cofm(c),mass(c));            //       sum up dipole       
	Ci->add_active_flag(c);                    //       add in activity flag
      }                                            //     END LOOP              
      LoopLeafKids(cell_iter,Ci,l) {               //     LOOP sub-leafs s      
	mon += mass(l);                            //       sum up monopole     
	com.add_times(cofm(l),mass(l));            //       sum up dipole       
	Ci->add_active_flag(l);                    //       add in activity flag
      }                                            //     END LOOP              
      if(is_active(Ci)) n++;                       //     count active cells    
      Ci->mass() = mon;                            //     set mass              
      mon        = (mon==zero)? zero:one/mon;      //     1/mass                
      com       *= mon;                            //     cofm = dipole/mass    
      Mset P(zero);                                //     reset multipoles      
      register real dmax(zero);                    //     reset d_max           
      LoopLeafKids(cell_iter,Ci,l) {               //     LOOP sub-leafs s      
	register vect Xi = cofm(l); Xi-= com;      //       distance vector     
	update_max(dmax,norm(Xi));                 //       update d_max^2      
	P.add_body(Xi,mass(l));                    //       add multipoles      
      }                                            //     END LOOP              
      if(has_leaf_kids(Ci)) dmax = sqrt(dmax);     //     d_max due to sub-leafs
      LoopCellKids(cell_iter,Ci,c) {               //     LOOP sub-cells c      
	register vect Xi = cofm(c); Xi-= com;      //       distance vector     
	register real Xq = norm(Xi);               //       distance^2          
	register real x  = dmax - rmax(c);         //       auxiliary           
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
  for(register int i=0; i!=TREE->N_cells(); ++i)   // LOOP cell sources         
    CELL_SRCE[i].normalize_poles();                //   normalize multipoles    
  // 4  set rcrit                                                               
  if(MAC) MAC->set_rcrit(this);                    // set r_crit for all cells  
  return n;                                        // return # active cells     
}
//------------------------------------------------------------------------------
#ifdef falcON_MPI
inline void grav_estimator::set_cell_sink(bool const&all) {
  register real*ci=CELL_COEF;                      // pter to cell's sink       
  if(all)                                          // IF all leafs              
    LoopCellsDown(grav::cell_iter,TREE,Ci)         //   LOOP cells              
      Ci->setCoeffs(ci++);                         //     set cell: pter->sink  
  else                                             // ELSE (only active)        
    LoopCellsDown(grav::cell_iter,TREE,Ci)         //   LOOP cells              
      if(is_active(Ci))                            //     IF(leaf is active)    
	Ci->setCoeffs(ci++);                       //     set cell: pter->sink  
      else                                         //     ELSE                  
	Ci->resetCoeffs();                         //     set cell: pter->sink  
}
#endif
//------------------------------------------------------------------------------
bool grav_estimator::prepare(const grav_mac*const&MAC,
			     bool           const&al,
			     bool           const&alloc_cell_coeffs)
{
  SET_I
  if(al) NLA_needed = TREE->N_leafs();             // all leafs are sink        
  if(NLA_needed==0) {
    falcON_WarningF("no body active","grav_estimator");
    return 1;
  }
  //  - allocate memory for leaf sink properties for active leafs               
  if(NLA!=NLA_needed) {                            // IF #active leafs changed  
    if(LEAF_SINK) delete[] LEAF_SINK;              //   delete old allocation   
    NLA = NLA_needed;                              //   # new allocation        
    LEAF_SINK=falcON_New(grav_leaf::sink_data,NLA);//   allocate memory         
  }                                                // ENDIF                     
  const bool all = al || NLA==TREE->N_leafs();     // are all active?           
  register grav_leaf::sink_data*si=LEAF_SINK;      // pter to leafs' sink data  
  if(all)                                          // IF all leafs              
    LoopLeafs(grav::leaf_type,TREE,Li) {           //   LOOP leafs              
      si->reset();                                 //     reset sink data       
      Li->set_sink(si++);                          //     set leaf: sink        
    }                                              //   END LOOP                
  else                                             // ELSE (only active)        
    LoopLeafs(grav::leaf_type,TREE,Li)             //   LOOP leafs              
      if(is_active(Li)) {                          //     IF(leaf is active)    
	si->reset();                               //       reset sink data     
	Li->set_sink(si++);                        //       set leaf: sink      
      } else                                       //     ELSE                  
	  Li->set_sink(0);                         //       leaf has no sink    
  // update cell source properties                                              
  if(!CELLS_UPTODATE ||                            // IF sources out-of-date    
     NCT != TREE->N_cells()) {                     //    OR changed in number   
    //    - allocate memory for cell source properties & reset their coeff pter 
    if(NCT     < TREE->N_cells() ||                //   IF #cells too large     
       NCT+NCT > TREE->N_cells()   ) {             //   OR #cells too small     
      if(CELL_SRCE) delete[] CELL_SRCE;            //     delete old allocation 
      NCT = TREE->N_cells();                       //     # new allocation      
      CELL_SRCE=falcON_New(grav_cell::srce_data,NCT);//   allocate memory       
    }                                              //   ENDIF                   
    register grav_cell::srce_data*ci=CELL_SRCE;    //   pter to cell's source   
    LoopCellsDown(grav::cell_iter,TREE,Ci) {       //   LOOP cells              
      Ci->set_srce(ci++);                          //     give memory to cell   
      Ci->resetCoeffs();                           //     reset cell: Coeffs    
    }                                              //   END LOOP                
    //      - pass source properties up the tree, count active cells            
    NCA = pass_up(MAC,TREE->is_re_used());         //   pass source data up tree
    CELLS_UPTODATE = 1;                            //   update up-to-date flag  
  } else {                                         // ELSE                      
    register grav_cell::srce_data*ci=CELL_SRCE;    //   pter to cell's source   
    LoopCellsDown(grav::cell_iter,TREE,Ci)         //   LOOP cells              
      Ci->set_srce(ci++);                          //     give memory to cell   
  }                                                // ENDIF                     
#ifdef falcON_MPI
  if(alloc_cell_coeffs) {                          // IF desired                
    if(CELL_COEF) delete[] CELL_COEF;              //   delete old cell coeffs  
    CELL_COEF = falcON_New(grav::Cset,NCA);        //   allocate memory         
    set_cell_sink(all);                            //   set cells' coeffs fields
  }                                                // ENDIF                     
#endif
  return all;                                      // return all                
}
//------------------------------------------------------------------------------
void grav_estimator::exact(bool       const&al
#ifdef falcON_ADAP
			  ,real       const&Nsoft,
			   uint       const&Nref,
			   real       const&emin,
			   real       const&efac
#endif
			   )
{
  if(GRAV==zero) {
    warning("[grav_estimator::exact()]: G=0\n");
    if(al) TREE->UseBodies(ResetBodiesGrav<1>());
    else   TREE->UseBodies(ResetBodiesGrav<0>());
    return;
  }
  update_leafs();
#ifdef falcON_ADAP
  adjust_eph(al,Nsoft,emin,EPS,Nref,efac);
#endif
  const bool all = prepare(0,al,0);
  if(N_active_cells()==0)
    return warning("[grav_estimator::exact()]: nobody active");
  STATS->reset();
#ifdef falcON_INDI
#  define ARGS KERNEL,STATS,EPS,0,INDI_SOFT
#else
#  define ARGS KERNEL,STATS,EPS,0
#endif
  if(all) {
    grav_iact_all K(ARGS);
    K.direct_summation(root());
    LoopLeafs(grav_leaf,TREE,Li) Li->normalize_grav();
  } else {
    grav_iact K(ARGS);
    K.direct_summation(root());
    LoopLeafs(grav_leaf,TREE,Li) if(is_active(Li)) Li->normalize_grav();
  }
#undef ARGS
#ifdef falcON_ADAP
  if(all) TREE->UseBodies(UpdateBodiesGrav<1>(TREE,GRAV,INDI_SOFT && Nsoft));
  else    TREE->UseBodies(UpdateBodiesGrav<0>(TREE,GRAV,INDI_SOFT && Nsoft));
#else
  if(all) TREE->UseBodies(UpdateBodiesGrav<1>(TREE,GRAV));
  else    TREE->UseBodies(UpdateBodiesGrav<0>(TREE,GRAV));
#endif
  TREE->mark_grav_usage();
}
//------------------------------------------------------------------------------
void grav_estimator::approx(const grav_mac*const&GMAC,
			    bool           const&al,
			    bool           const&split
#ifdef falcON_ADAP
			    ,
			    real           const&Nsoft,
			    uint           const&Nref,
			    real           const&emin,
			    real           const&efac
#endif
			    )
{
  if(GRAV==zero) {
    warning("[grav_estimator::approx()]: G=0\n");
    if(al) TREE->UseBodies(ResetBodiesGrav<1>());
    else   TREE->UseBodies(ResetBodiesGrav<0>());
    return;
  }
  SET_I
  report REPORT("grav_estimator::approx()");
  update_leafs();
  SET_T(" time: grav_estimator::update_leafs():  ");
#ifdef falcON_ADAP
  adjust_eph(al,Nsoft,emin,EPS,Nref,efac);         //[adjust eps_i/2 of leafs]  
  SET_T(" time: grav_estimator::adjust_eph():    ");
#endif
  const bool all = prepare(GMAC,al,0);             // prepare tree              
  if(!all && N_active_cells()==0)
    return warning("[grav_estimator::approx()]: nobody active");
  SET_T(" time: grav_estimator::prepare():       ");
  report REPORT2("interaction & evaluation");
  STATS->reset();
#ifdef falcON_INDI
#  define ARGS KERNEL,STATS,EPS,Ncsize,INDI_SOFT,DIR
#else
#  define ARGS KERNEL,STATS,EPS,Ncsize,DIR
#endif
  Ncsize = 4+(all? TREE->N_cells() : N_active_cells())/(split? 16:4);
  if(all) {                                        // IF all are active         
    grav_iact_all GK(ARGS);                        //   init gravity kernel     
    MutualInteractor<grav_iact_all> MI(&GK,split? TREE->depth()-1 : 
				                  TREE->depth());
                                                   //   init mutual interactor  
    if(split) {                                    //   IF splitting            
      LoopCellKids(cell_iter,root(),c1) {          //     LOOP cell kids c1     
	MI.cell_self(c1);                          //       self-iaction c1     
	LoopCellSecd(cell_iter,root(),c1+1,c2)     //       LOOP kids c2>c1     
	  MI.cell_cell(c1,c2);                     //         interaction c1,2  
	LoopLeafKids(cell_iter,root(),s2)          //       LOOP leaf kids s    
	  MI.cell_leaf(c1,s2);                     //         interaction c1,s  
	GK.evaluate(c1);                           //       evaluation phase    
      }                                            //     END LOOP              
      LoopLeafKids(cell_iter,root(),s1) {          //     LOOP leaf kids s1     
	LoopLeafSecd(cell_iter,root(),s1+1,s2)     //       LOOP kids s2>s1     
	  GK.interact(s1,s2);                      //         interaction s1,2  
	s1->normalize_grav();                      //       evaluation phase    
      }                                            //     END LOOP              
    } else {                                       //   ELSE                    
      MI.cell_self(root());                        //     interaction phase     
      GK.evaluate(root());                         //     evaluation phase      
    }                                              //   ENDIF                   
    Ncoeffs = GK.coeffs_used();                    //   remember # coeffs used  
    Nchunks = GK.chunks_used();                    //   remember # chunks used  
  } else {                                         // ELSE: not all are active  
    grav_iact GK(ARGS);                            //   init gravity kernel     
#undef ARGS
    MutualInteractor<grav_iact> MI(&GK,split? TREE->depth()-1 :
					      TREE->depth());
                                                   //   init mutual interactor  
    if(split) {                                    //   IF splitting            
      LoopCellKids(cell_iter,root(),c1) {          //     LOOP cell kids c1     
	if(is_active(c1)) {                        //      IF active s1:        
	  MI.cell_self(c1);                        //       self-iaction c1     
	  LoopCellSecd(cell_iter,root(),c1+1,c2)   //       LOOP kids c2>c1     
	    MI.cell_cell(c1,c2);                   //         interaction c1,2  
	  LoopLeafKids(cell_iter,root(),s2)        //       LOOP leaf kids s    
	    MI.cell_leaf(c1,s2);                   //         interaction c1,s  
	  GK.evaluate(c1);                         //       evaluation phase    
	} else {                                   //      ELSE: inactive c1    
	  LoopCellSecd(cell_iter,root(),c1+1,c2)   //       LOOP kids c2>c1     
	    if(is_active(c2))MI.cell_cell(c1,c2);  //         interaction c1,2  
	  LoopLeafKids(cell_iter,root(),s2)        //       LOOP leaf kids s    
	    if(is_active(s2))MI.cell_leaf(c1,s2);  //         interaction c1,s  
	}                                          //      ENDIF                
      }                                            //     END LOOP              
      LoopLeafKids(cell_iter,root(),s1) {          //     LOOP leaf kids s1     
	if(is_active(s1)) {                        //      IF active s1:        
	  LoopLeafSecd(cell_iter,root(),s1+1,s2)   //       LOOP kids s2>s1     
	    GK.interact(s1,s2);                    //         interaction s1,2  
	  s1->normalize_grav();                    //       evaluation phase    
	} else {                                   //      ELSE: inactive s1    
	  LoopLeafSecd(cell_iter,root(),s1+1,s2)   //       LOOP kids s2>s1     
	    if(is_active(s2)) GK.interact(s1,s2);  //         interaction s1,2  
	}                                          //      ENDIF                
      }                                            //     END LOOP              
    } else {                                       //   ELSE                    
      MI.cell_self(root());                        //     interaction phase     
      GK.evaluate(root());                         //     evaluation phase      
    }                                              //   ENDIF                   
    Ncoeffs = GK.coeffs_used();                    //   remember # coeffs used  
    Nchunks = GK.chunks_used();                    //   remember # chunks used  
  }                                                // ENDIF                     
  SET_T(" time: interaction & evaluation:        ");
#ifdef falcON_ADAP
  if(all) TREE->UseBodies(UpdateBodiesGrav<1>(TREE,GRAV,INDI_SOFT && Nsoft));
  else    TREE->UseBodies(UpdateBodiesGrav<0>(TREE,GRAV,INDI_SOFT && Nsoft));
#else
  if(all) TREE->UseBodies(UpdateBodiesGrav<1>(TREE,GRAV));
  else    TREE->UseBodies(UpdateBodiesGrav<0>(TREE,GRAV));
#endif
  TREE->mark_grav_usage();
  SET_T(" time: updating bodies gravity:         ");
}
//------------------------------------------------------------------------------
namespace {
  using namespace nbdy; using nbdy::uint;
  //============================================================================
  uint NX;
#if falcON_NDIM==3
  real pdim(real const&x) { return cube(x); }
#else
  real pdim(real const&x) { return square(x); }
#endif
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
  class UpdateBodiesRho {
    const bool     all;
    const oct_tree*T;
  public:
    UpdateBodiesRho(const oct_tree*const&t,
		    bool           const&a) : T(t), all(a) {}
    template<typename bodies_type>
    uint operator() (const bodies_type*const&B) const {
      if(all)
	LoopLeafs(grav_leaf,T,Li)
	  Li->copy_to_bodies_rho(B);
      else
	LoopLeafs(grav_leaf,T,Li) if(is_active(Li))
	  Li->copy_to_bodies_rho(B);
      return 0;
    }
  };
}                                                  // END: unnamed namespace    
//------------------------------------------------------------------------------
void grav_estimator::estimate_nd(bool const&al, uint const&Nx) const
{
  NX = Nx;
  if(al) guess<number_density,1>::do_it(root(),zero);
  else   guess<number_density,0>::do_it(root(),zero);
  TREE->UseBodies(UpdateBodiesRho(TREE,al));
}
//------------------------------------------------------------------------------
void grav_estimator::estimate_sd(bool const&al, uint const&Nx)
{
  update_leafs();
  prepare(0,al,0);
  NX = Nx;
  if(al) guess<surface_density,1>::do_it(root(),zero);
  else   guess<surface_density,0>::do_it(root(),zero);
  TREE->UseBodies(UpdateBodiesRho(TREE,al));
  TREE->mark_grav_usage();
}
//------------------------------------------------------------------------------
void grav_estimator::estimate_md(bool const&al, uint const&Nx)
{
  update_leafs();
  prepare(0,al,0);
  NX = Nx;
  if(al) guess<mass_density,1>::do_it(root(),zero);
  else   guess<mass_density,0>::do_it(root(),zero);
  TREE->UseBodies(UpdateBodiesRho(TREE,al));
  TREE->mark_grav_usage();
}
//------------------------------------------------------------------------------
void grav_estimator::dump_cells(std::ostream&o) const
{
  if(CELL_SRCE) TREE->dump_cells<grav_cell>(o);
  else          TREE->dump_cells<basic_cell>(o);
}
//------------------------------------------------------------------------------
void grav_estimator::dump_leafs(std::ostream&o) const
{
  TREE->dump_leafs<grav_leaf>(o);
}
////////////////////////////////////////////////////////////////////////////////

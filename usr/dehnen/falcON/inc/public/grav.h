// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// grav.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2002                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines                                                                     |
//                                                                             |
// class grav_iact                                                             |
// class grav_iact_s                                                           |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef included_grav_h
#define included_grav_h

#ifndef included_grat_h
#  include <public/grat.h>
#endif
#ifndef included_deft_h
#  include <public/deft.h>
#endif

#if defined(REAL_IS_FLOAT) && defined(USE_SSE)
#  define  SSE_CODE
#  include <proper/kern_sse.h>
#  define  grav_kern grav_kern_sse
#else
#  undef   SSE_CODE
#  include <public/kern.h>
#endif

////////////////////////////////////////////////////////////////////////////////
#define P_ORDER 3                                  // expansion order is fixed !
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::grav_iact                                                      
  //                                                                            
  // This class is at the heart of the algorithm. It serves as INTERACTOR in    
  // the template class MutualInteract<>, defined in iact.h, which encodes      
  // the interaction phase in a most general way.                               
  // class nbdy::grav_iact has member functions for soul-soul, soul-cell,       
  // cell-soul, cell-cell, and cell-self interactions (methods interact()), as  
  // well as a function for the evaluation phase, method evaluate_grav().       
  //                                                                            
  // NOTE. We do NOT allocate memory for the cell's Taylor coefficients. Rather 
  // this has to be done before the first call to any interact().               
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class grav_iact : public grav_kern {
    grav_iact           (grav_iact const&);        // not implemented           
    grav_iact& operator=(grav_iact const&);        // not implemented           
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
  protected:
    grav_stat* const      STAT;                    // statistics                
    int                   N_PRE[3], N_POST[3];     // direct sums control       
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
  public:
    typedef grav_tree::cell_iterator  cell_iter;   // cell iterator             
    typedef grav_tree::soul_iterator  soul_iter;   // soul iterator             
    //--------------------------------------------------------------------------
  protected:
    struct TaylorSeries {
              vect X;                              // expansion center          
      mutable real C[N_COEFF];                     // coefficients              
      TaylorSeries(vect const&);
      TaylorSeries(TaylorSeries const&);
//       TaylorSeries(vect const&x) : X(x)
//       { for(register int i=0; i!=N_COEFF; ++i) C[i] = zero; }
//       TaylorSeries(TaylorSeries const&T) : X(T.X)
//       { for(register int i=0; i!=N_COEFF; ++i) C[i] = T.C[i]; }
    };
    //--------------------------------------------------------------------------
    real&      c0(real*         const&T) const { return *(T); }
    real const&c0(TaylorSeries  const&T) const { return *(T.C); }
    ten1       c1(real*         const&T) const { return ten1(T+N_C1); }
    ten1       c1(TaylorSeries  const&T) const { return c1(T.C); }
    ten2       c2(real*         const&T) const { return ten2(T+N_C2); }
    ten2       c2(TaylorSeries  const&T) const { return c2(T.C); }
    ten3       c3(real*         const&T) const { return ten3(T+N_C3); }
    ten3       c3(TaylorSeries  const&T) const { return c3(T.C); }
#if P_ORDER > 3
    ten4       c4(real*         const&T) const { return ten4(T+N_C4); }
    ten4       c4(TaylorSeries  const&T) const { return c4(T.C); }
#if P_ORDER > 4
    ten5       c5(real*         const&T) const { return ten5(T+N_C5); }
    ten5       c5(TaylorSeries  const&T) const { return c5(T.C); }
#if P_ORDER > 5
    ten6       c6(real*         const&T) const { return ten6(T+N_C6); }
    ten6       c6(TaylorSeries  const&T) const { return c6(T.C); }
#endif
#endif
#endif
    //--------------------------------------------------------------------------
    void shift_and_add  (TaylorSeries&, vect const&,
			 const real* const&, real const&) const;
    void extract_grav   (TaylorSeries const&, soul_iter const&) const;
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
  protected:
    //--------------------------------------------------------------------------
    bool do_direct_pre (cell_iter const&A, soul_iter const&B) const {
      return number(A) < N_PRE[0];
    }
    //--------------------------------------------------------------------------
    bool do_direct_post(cell_iter const&A, soul_iter const&B) const {
      return is_twig(A) || number(A) < N_POST[0];
    }
    //--------------------------------------------------------------------------
    bool do_direct_post(cell_iter const&A, cell_iter const&B) const {
      return (is_twig(A) && is_twig(B)) ||
	     (number(A) < N_POST[1] && number(B) < N_POST[1]);
    }
    //--------------------------------------------------------------------------
    bool do_direct(cell_iter const&A) const { 
      return is_twig(A) || number(A) < N_PRE[2];
    }
    //--------------------------------------------------------------------------
    static bool well_separated(cell_iter const&A, cell_iter const&B,
			       real      const&Rq)
    { return Rq > square(rcrit(A)+rcrit(B)); }
    //--------------------------------------------------------------------------
    static bool well_separated(cell_iter const&A, soul_iter const&B,
			       real      const&Rq)
    { return Rq > rcrit2(A); }
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
  private:
    void shift_by (real*  const&, vect const&) const;
    void shift_to (TaylorSeries&, vect const&) const;
    void eval_grav(cell_iter const&, TaylorSeries const&) const;
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    grav_iact(
	      grav_stat*const&t,                      // I: statistics          
	      real            e,                      // I: softening length    
	      kern_type       k    = Default::kernel, //[I: type of kernel]     
#ifdef ALLOW_INDI
	      soft_type       s    = Default::soften, //[I: type of softening]  
#endif
	      int const       nd[4]= Default::direct, //[I: direct sum control] 
	      bool            fp   = false) :         //[I: use Pth pole in pot]
#ifdef ALLOW_INDI
      grav_kern( k, e, s, fp ),
#else
      grav_kern( k, e, fp ),
#endif
      STAT     ( t )                               // set statistics            
    {
      N_PRE [0] = nd[0];                           // C-B direct before try     
      N_PRE [1] = 0u;                              // C-C direct before try     
      N_PRE [2] = nd[3];                           // C-S direct before try     
      N_POST[0] = nd[1];                           // C-B direct after fail     
      N_POST[1] = nd[2];                           // C-C direct after fail     
      N_POST[2] = 0u;                              // C-S direct after fail     
    }
    //--------------------------------------------------------------------------
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      return is_twig(B) || (!is_twig(A) && rmax(A) > rmax(B));
    }
    //--------------------------------------------------------------------------
    // interaction phase                                                        
    //--------------------------------------------------------------------------
    void direct_summation(cell_iter const&A) const { 
      grav_kern::many(A);
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A) const {
      if(!is_sink(A)) return true;                 // no interaction -> DONE    
      if(do_direct(A)) {                           // IF(suitable)             >
	grav_kern::many(A);                        //   perform BB iactions     
	STAT->record_direct_CS(A);                 //   record stats            
	return true;                               //   DONE                    
      }                                            // < ELSE >                  
      return false;                                //   must be splitted        
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A, cell_iter const&B) const {
      if(!(is_sink(A) || is_sink(B))) return true; // no interaction -> DONE    
      register vect dX = cofm(A)-cofm(B);          // compute dX = X_A - X_B    
      register real Rq = norm(dX);                 // and dX^2                  
      if(well_separated (A,B,Rq)) {                // IF(well separated)      > 
	grav_kern::grav(A,B,dX,Rq);                //   interact                
	STAT->record_taylor_CC(A,B);               //   record stats            
	return true;                               //   DONE                    
      }                                            // <                         
      if(do_direct_post(A,B)) {                    // IF(suitable)            > 
	grav_kern::many(A,B);                      //   perform BB iactions     
	STAT->record_direct_CC(A,B);               //   record stats            
	return true;                               //   DONE                    
      }                                            // < ELSE  >                 
      return false;                                //   SPLIT <                 
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A, soul_iter const&B) const {
      if(!(is_sink(A) || is_sink(B))) return true; // no interaction -> DONE    
      if(do_direct_pre(A,B)) {                     // IF(suitable)             >
	many(A,B);                                 //   perform BB iactions     
	STAT->record_direct_CB(A,B);               //   record statistics       
	return true;                               //   DONE                    
      }                                            // <                         
      register vect dX = cofm(A)-cofm(B);          // compute R = x_A-x_B       
      register real Rq = norm(dX);                 // compute R^2               
      if(well_separated(A,B,Rq)) {                 // IF(well separated)       >
	grav(A,B,dX,Rq);                           //   interact                
	STAT->record_taylor_CB(A,B);               //   record statistics       
	return true;                               //   DONE                    
      }                                            // <                         
      if(do_direct_post(A,B)) {                    // IF(suitable)             >
	many(A,B);                                 //   perform BB iactions     
	STAT->record_direct_CB(A,B);               //   record statistics       
	return true;                               //   DONE                    
      }                                            // < ELSE                    
      return false;                                //   cell must be splitted   
    }
    //--------------------------------------------------------------------------
    bool interact(soul_iter const&A, cell_iter const&B) const {
      return interact(B,A);
    }
    //--------------------------------------------------------------------------
    void interact(soul_iter const&A, soul_iter const&B) const {
      if(!(is_sink(A) || is_sink(B))) return;      // no interaction -> DONE    
      grav_iact::single(A,B);                      // perform interaction       
      STAT->record_BB();                           // record statistics         
    }  
    //--------------------------------------------------------------------------
    void evaluate(cell_iter const&C) const;        // evaluation phase          
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::grav_iact_s                                                    
  //                                                                            
  // This class is similar to grav_iact. It differs only in that it organizes   
  // the cells' Taylor coefficients by itself: at a cell's first interaction,   
  // memory for its coefficients is taken from a pre-allocated pool and returned
  // to the pool when the cell eventually passes through evaluation phase.      
  // Together with the interwaeving of interaction and evalutation phase, this  
  // saves about 20b/body on typical applications.                              
  //                                                                            
  // NOTE. This code is useful predominantly for serial code.                   
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class grav_iact_s : public grav_iact {
  private:
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
    nbdy::pool   *POOL;                             // pool for TaylorCoeffs    
    mutable int   NC, MAXNC;                        // # TaylorCoeffs ever used 
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
    void give_coeffs(cell_iter const&C) const {
      if(coeffs(C)==0) {
	register real*X = static_cast<real*>(POOL->alloc());
	for(register int i=0; i!=N_COEFF; ++i) X[i] = zero;
	C->coeffs() = X;
	++NC;
      }
    }
    //--------------------------------------------------------------------------
    void take_coeffs(cell_iter const&C) const {
      if(coeffs(C)) {
	POOL->free(C->coeffs());
	C->coeffs() = 0;
	--NC;
      }
    }
    //--------------------------------------------------------------------------
    void eval_grav  (cell_iter const&, TaylorSeries const&) const;
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    grav_iact_s(
	      grav_stat*const&t,                      // I: statistics          
	      real            e,                      // I: softening length    
	      uint            np,                     // I: initial pool size   
	      kern_type       k    = Default::kernel, //[I: type of kernel]     
#ifdef ALLOW_INDI
	      soft_type       s    = Default::soften, //[I: type of softening]  
#endif
	      int const       nd[4]= Default::direct, //[I: direct sum control] 
	      bool            fp   = false) :         //[I: use Pth pole in pot]
#ifdef ALLOW_INDI
      grav_iact ( t,e,k,s,nd,fp ),
#else
      grav_iact ( t,e,k,nd,fp ),
#endif
      NC        ( 0 ), MAXNC( 0 )
    {
      MemoryCheck( POOL=new nbdy::pool(np,N_COEFF*sizeof(real)) );
    }
    //--------------------------------------------------------------------------
    ~grav_iact_s() { delete POOL; }
    //--------------------------------------------------------------------------
    const  int& coeffs_used  () const { return MAXNC; }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A) const {
      return grav_iact::interact(A);
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A, cell_iter const&B) const {
      if(!(is_sink(A) || is_sink(B))) return true; // no interaction -> DONE    
      register vect dX = cofm(A)-cofm(B);          // compute dX = X_A - X_B    
      register real Rq = norm(dX);                 // and dX^2                  
      if(well_separated (A,B,Rq)) {                // IF(well separated)      > 
	if(is_sink(A)) give_coeffs(A);             //   ensure coeffs exist     
	if(is_sink(B)) give_coeffs(B);             //   ensure coeffs exist     
	grav_kern::grav(A,B,dX,Rq);                //   interact                
	STAT->record_taylor_CC(A,B);               //   record stats            
	return true;                               //   DONE                    
      }                                            // <                         
      if(do_direct_post(A,B)) {                    // IF(suitable)            > 
	grav_kern::many(A,B);                      //   perform BB iactions     
	STAT->record_direct_CC(A,B);               //   record stats            
	return true;                               //   DONE                    
      }                                            // < ELSE  >                 
      return false;                                //   SPLIT <                 
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A, soul_iter const&B) const {
      if(!(is_sink(A) || is_sink(B))) return true; // no interaction -> DONE    
      if(do_direct_pre(A,B)) {                     // IF(suitable)             >
	many(A,B);                                 //   perform BB iactions     
	STAT->record_direct_CB(A,B);               //   record statistics       
	return true;                               //   DONE                    
      }                                            // <                         
      register vect dX = cofm(A)-cofm(B);          // compute R = x_A-x_B       
      register real Rq = norm(dX);                 // compute R^2               
      if(well_separated(A,B,Rq)) {                 // IF(well separated)       >
	if(is_sink(A)) give_coeffs(A);             //   ensure coeffs exist     
	grav(A,B,dX,Rq);                           //   interact                
	STAT->record_taylor_CB(A,B);               //   record statistics       
	return true;                               //   DONE                    
      }                                            // <                         
      if(do_direct_post(A,B)) {                    // IF(suitable)             >
	many(A,B);                                 //   perform BB iactions     
	STAT->record_direct_CB(A,B);               //   record statistics       
	return true;                               //   DONE                    
      }                                            // < ELSE                    
      return false;                                //   cell must be splitted   
    }
    //--------------------------------------------------------------------------
    bool interact(soul_iter const&A, cell_iter const&B) const {
      return interact(B,A);
    }
    //--------------------------------------------------------------------------
    void interact(soul_iter const&A, soul_iter const&B) const {
      return grav_iact::interact(A,B);
    }
    //--------------------------------------------------------------------------
    void evaluate(cell_iter const&) const;         // evaluation phase          
    //--------------------------------------------------------------------------
  };
}
////////////////////////////////////////////////////////////////////////////////
#endif                                             // included_grav_h           

// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// step.h                                                                      |
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
// defines                                                                     |
//                                                                             |
// class LeapFrog            fully inline                                      |
// class BlockStep           fully inline                                      |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_step_h
#define falcON_included_step_h

#ifndef falcON_included_iostream
#  include <iostream>
#  define falcON_included_iostream
#endif

#ifndef falcON_included_iomanip
#  include <iomanip>
#  define falcON_included_iomanip
#endif

#ifndef falcON_included_body_h
#  include <body.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::LeapFrog                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class LeapFrog {
    LeapFrog(const LeapFrog&);                     // not implemented           
    LeapFrog operator= (const LeapFrog&);          // not implemented           
    //--------------------------------------------------------------------------
    // data members                                                             
    //--------------------------------------------------------------------------
  private:
    real         TAU, TAUH;                        // time step                 
    mutable real TIME;                             // simulation time           
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    LeapFrog(const int f, const real t)
      : TAU (pow(half,f)), TAUH (half*TAU), TIME(t) {}
    const real& time() const { return TIME; }
    const real& tau () const { return TAU; }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void predict    (const bodies_type* const&BB) const
    {
      LoopBodies(bodies,BB,Bi) {
	Bi.vel().add_mul(acc(Bi),TAUH);            // v_1/2 = v_0 +a_0 * dt/2   
	Bi.pos().add_mul(vel(Bi),TAU );            // x_1   = x_0 +v_1/2 * dt   
      }
      TIME += TAU;                                 // t_1   = t_0 + dt          
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void accelerate (const bodies_type* const&BB) const
    {
      LoopBodies(bodies,BB,Bi)
	Bi.vel().add_mul(acc(Bi),TAUH);            // v_1   = v_1/2 * a_1 * dt/2
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::BlockStep                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class BlockStep {
    BlockStep(const BlockStep&);                   // not implemented           
    BlockStep& operator= (const BlockStep&);       // not implemented           
    //--------------------------------------------------------------------------
    // private types                                                            
    //--------------------------------------------------------------------------
  private:
    struct TimeStep {
      mutable uint  N;                             // # bodies in step          
      real          TAU, TAUH, TAUSQ;              // tau, tau/2, tau^2         
      void re_set(const int f=0)        { TAU   = pow(half,f);
                                          TAUH  = half*TAU;
                                          TAUSQ = TAU*TAU; }
      TimeStep(const int f=0)           { N=0u; re_set(f); }
      void reset_count ()         const { N=0u; }
      void dump (std::ostream&to) const { to<<N<<" "<<TAU<<" "<<"\n"; }
      void join ()                const { ++N; }
      void leave()                const { --N; }
    };
    //--------------------------------------------------------------------------
    // data members                                                             
    //--------------------------------------------------------------------------
  private:
    int            H0;                             // sets longest time step    
    int            NSTEPS, HIGHEST;                // # steps, #steps-1         
    TimeStep      *TS, *SHORTEST;                  // arrays with time steps    
    mutable real   TIME;                           // simulation time           
  protected:
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
    void        to_lower      (body&B, int const&low) const
    {
      if(level(B) > low) {                         // IF(above lowest active)  >
	TS[level(B)].leave();                      //   leave old level         
	B.level()--;                               //   set new level           
	TS[level(B)].join();                       //   join new level          
      }                                            // <                         
    }
    //--------------------------------------------------------------------------
    void        to_higher     (body&B) const
    {
      if(level(B) < highest_level()) {             // IF(below highest level)  >
	TS[level(B)].leave();                      //   leave old level         
	B.level()++;                               //   set new level           
	TS[level(B)].join();                       //   join new level          
      }                                            // <                         
    }
    //--------------------------------------------------------------------------
    void        reset_steps   (int  const&h0,      // I: h0, tau_max = 2^-h0    
			       int  const&Ns,      // I: #steps                 
			       real const&t =zero) // I: actual time            
    {
      if(TS) delete[] TS;
      NSTEPS   = abs(Ns);
      HIGHEST  = Ns-1;
      TS       = falcON_New(TimeStep,NSTEPS);
      SHORTEST = TS+HIGHEST;
      TIME     = t;
      H0       = h0;
      register int h,i;
      for(i=0, h=H0; i!=NSTEPS; ++i,++h) (TS+i)->re_set(h);
    }
    //--------------------------------------------------------------------------
    // construction & destruction                                               
    //--------------------------------------------------------------------------
    BlockStep    () : TS(0) {}
    //--------------------------------------------------------------------------
    BlockStep    (int  const&h0,                   // I: h0, tau_max = 2^-h0    
		  int  const&Ns,                   // I: #steps                 
		  real const&t)                    // I: actual time            
      : TS ( 0 )
    {
      reset_steps (h0,Ns,t);
    }
    //--------------------------------------------------------------------------
    ~BlockStep                  () { if(TS) delete[] TS;}
    //--------------------------------------------------------------------------
    // other protected methods                                                  
    //--------------------------------------------------------------------------
    void        dump          (std::ostream&to)  const
    {
      for(register int i=0; i!=NSTEPS; ++i) {
	to<<i<<": ";
	(TS+i)->dump(to);
      }
    }
    //--------------------------------------------------------------------------
    void        clock_on      ()             const
    {
      TIME += tau_min();
    }
    //--------------------------------------------------------------------------
    void        reset_counts  () const
    {
      for(register TimeStep* s=TS; s<=SHORTEST; s++) s->reset_count();
    }
    //--------------------------------------------------------------------------
    bool        need_move     (int const&low) const
    {
      for(register TimeStep* s=TS+low; s<=SHORTEST; s++) if(s->N) return true;
      return false;
    }
    //--------------------------------------------------------------------------
    void        change_level  (body&B, indx const&l) const
    {
      TS[level(B)].leave();                        // leave old level           
      B.level() = l;                               // set new level             
      TS[level(B)].join ();                        // join new level            
    }
    //--------------------------------------------------------------------------
    int         longest_moving(uint const&T) const
    {
      register int i=HIGHEST, t=T;
      for(; !(t&1) && i>0; t>>=1, --i);
      return i;
    }
    //--------------------------------------------------------------------------
    void        add_to_level  (const indx l) const { (TS+l)->join (); }
    bool        is_leap_frog  ()             const {return NSTEPS == 1; }
    const uint& number        (const indx l) const {return(TS+l)->N; }
    const real& tau           (const indx l) const {return(TS+l)->TAU; }
    const real& tau           (const body B) const {return(TS+level(B))->TAU; }
    const real& tau_h         (const body B) const {return(TS+level(B))->TAUH; }
    const real& tausq         (const body B) const {return(TS+level(B))->TAUSQ;}
    const real& time          ()             const {return TIME; }
    //--------------------------------------------------------------------------
    void        half_kick     (body&B) const
    {
      B.vel().add_mul(acc(B),tau_h(B));              // v -> v + h/2 * a        
    }
    //--------------------------------------------------------------------------
    void        tiny_drift    (body&B) const
    {
      B.pos().add_mul(vel(B),tau_min());             // x -> x + h_min * v      
    }
    //--------------------------------------------------------------------------
    void        drift_by      (body&B, real const&dt) const
    {
      B.pos().add_mul(vel(B),dt);                    // x -> x + dt * v         
    }
    //--------------------------------------------------------------------------
    void        acce_half     (body&B) const { half_kick(B); }
    void        move_tiny     (body&B) const { tiny_drift(B); }
    void        move_by       (body&B, real const&dt) const { drift_by(B,dt); }
    //--------------------------------------------------------------------------
    int         highest_level ()             const {return HIGHEST; }
    const real& tau_min       ()             const {return SHORTEST->TAU; }
    const real& tau_max       ()             const {return TS->TAU; }
    const int & h0            ()             const {return H0; }
    const uint& No_in_level   (const indx l) const {return(TS+l)->N; }
    //--------------------------------------------------------------------------
    void        short_stats   (std::ostream&to, int const&w=7) const
    {
      for(register int l=0; l!=NSTEPS; ++l) to<<std::setw(w)<<number(l)<<" ";
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::GravBlockStep                                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class GravBlockStep : public BlockStep {
    GravBlockStep(const GravBlockStep&);            // not implemented          
    GravBlockStep& operator= (const GravBlockStep&);// not implemented          
    //--------------------------------------------------------------------------
    // private types                                                            
    //--------------------------------------------------------------------------
  private:
    enum {
      use_a    = 1,
      use_p    = 2,
      use_ap   = 3,
      use_c    = 4,
      use_ac   = 5,
      use_pc   = 6,
      use_apc  = 7,
      use_e    = 8,
      use_ae   = 9,
      use_pe   = 10,
      use_ce   = 12,
      use_ape  = 11,
      use_ace  = 13,
      use_pce  = 14,
      use_apce = 15,
    };
    //--------------------------------------------------------------------------
    // data members                                                             
    //--------------------------------------------------------------------------
    int            SCH;                            // stepping scheme           
    real           FAQ,FPQ,FCQ,FEQ;                // factors^2 for stepping    
  protected:
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
    real        step_sq       (body const&B, real const&eps=one)   const
    {
      if(SCH == 0)     return zero;
      if(SCH == use_p) return FPQ/square(pot(B));
      else {
	register real ia = one/norm(acc(B)), tq=1.e7;
	if(SCH & use_a) update_min(tq, FAQ*ia);
	if(SCH & use_p) update_min(tq, FPQ/square(pot(B)));
	if(SCH & use_c) update_min(tq, FCQ*abs(pot(B))*ia);
	if(SCH & use_e) update_min(tq, FEQ*eps*sqrt(ia));
	return tq;
      }
    }
    //--------------------------------------------------------------------------
    void        assign_level  (body&B, real const&eps=one)           const
    {
      register real t = sqrt(two*step_sq(B,eps)); // tau * 2^(1/2)              
      register indx l = 0;                        // try longest step           
      while(tau(l)          > t  &&               // WHILE shorter step desired 
	    highest_level() > l) ++l;             // AND possible: decrease     
      add_to_level(B.level()=l);                  // account for in BlockStep   
    }
    //--------------------------------------------------------------------------
    void        adjust_level  (body&B, int const&low, real const&eps=one) const
      // implementing scheme 2 of my notes: change only by factors of two       
    {
      const double root_half=0.7071067811865475244;// sqrt(1/2)                 
      register real tq=step_sq(B,eps)* root_half;  // tau^2 * sqrt(1/2)         
      if     (tausq(B) < tq)    to_lower (B,low);  // change to longer  if poss 
      else if(tausq(B) > tq+tq) to_higher(B);      // change to shorter if poss 
    }
    //--------------------------------------------------------------------------
    void        reset_scheme  (real const&fa,
			       real const&fp=zero,
			       real const&fc=zero,
			       real const&fe=zero)
    {
      SCH = 0;
      FAQ = fa*fa; if(FAQ) SCH |= use_a;
      FPQ = fp*fp; if(FPQ) SCH |= use_p;
      FCQ = fc*fc; if(FCQ) SCH |= use_c;
      FEQ = fe*fe; if(FEQ) SCH |= use_e;
      if(SCH==0) 
	error("[%s.%d]: in GravBlockStep::reset_scheme(%f,%f,%f,%f): "
	      "time step control factors=0",__FILE__,__LINE__,fa,fp,fc,fe);
    }
    //--------------------------------------------------------------------------
    // construction & destruction                                               
    //--------------------------------------------------------------------------
    GravBlockStep() : BlockStep() {}
    //--------------------------------------------------------------------------
    GravBlockStep(int  const&h0,                   // O: h0, tau_max = 2^-h0    
		  int  const&Ns,                   // I: #steps                 
		  real const&t,                    // I: actual time            
		  real const&fa,                   // I: f_acc                  
		  real const&fp = zero,            // I: f_phi                  
		  real const&fc = zero,            // I: f_pa                   
		  real const&fe = zero)            // I: f_ea                   
      : BlockStep(h0,Ns,t)
    {
      reset_scheme(fa,fp,fc,fe);
    }
    //--------------------------------------------------------------------------
  };
}
////////////////////////////////////////////////////////////////////////////////
#endif   // falcON_included_step_h

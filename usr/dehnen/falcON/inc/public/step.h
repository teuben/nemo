// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// step.h                                                                      |
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
// class LeapFrog            fully inline                                      |
// class BlockStep           fully inline                                      |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef included_step_h
#define included_step_h

#ifndef included_iostream
#  include <iostream>
#  define included_iostream
#endif

#ifndef included_iomanip
#  include <iomanip>
#  define included_iomanip
#endif

#ifndef included_body_h
#  include <body.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  // class nbdy::LeapFrog                                                       
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
    void predict    (const bodies*) const;         // before force calculation  
    void accelerate (const bodies*) const;         // after force calculation   
  };
  //////////////////////////////////////////////////////////////////////////////
  // class nbdy::BlockStep                                                      
  //////////////////////////////////////////////////////////////////////////////
  class BlockStep {
    BlockStep(const BlockStep&);                   // not implemented           
    BlockStep& operator= (const BlockStep&);       // not implemented           
    //--------------------------------------------------------------------------
    // friendships                                                              
    //--------------------------------------------------------------------------
    friend class Stepping;
    //--------------------------------------------------------------------------
    // private types                                                            
    //--------------------------------------------------------------------------
  private:
    enum {
      use_a   = 1,
      use_p   = 2,
      use_c   = 4,
      use_ap  = 3,
      use_ac  = 5,
      use_pc  = 6,
      use_apc = 7
    };
    //--------------------------------------------------------------------------
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
    int            SCH;                            // stepping scheme           
    int            H0;                             // sets longest time step    
    int            NSTEPS, HIGHEST;                // # steps, #steps-1         
    TimeStep      *TS, *SHORTEST;                  // arrays with time steps    
    mutable real   TIME;                           // simulation time           
    real           FAQ,FPQ,FCQ;                    // factors^2 for stepping    
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
    void        to_lower      (body&, const int) const;
    void        to_higher     (body&) const;
  protected:
    //--------------------------------------------------------------------------
    // construction & destruction                                               
    //--------------------------------------------------------------------------
    BlockStep    () : TS(0) {}
    //--------------------------------------------------------------------------
    BlockStep    (const int,                       // O: h0, tau_max = 2^-h0    
		  const int,                       // I: #steps                 
		  const real,                      // I: actual time            
		  const real,                      // I: f_acc                  
		  const real = zero,               // I: f_phi                  
		  const real = zero);              // I: f_pa                   
    //--------------------------------------------------------------------------
    ~BlockStep                  () { if(TS) delete[] TS;}
    //--------------------------------------------------------------------------
    // other protected methods                                                  
    //--------------------------------------------------------------------------
    void        dump          (std::ostream&)const;
    void        clock_on      ()             const;
    real        step_sq       (const body)   const;
    void        assign_level  (body&)        const;
    void        reset_counts  () const;
    bool        need_move     (const int)    const;
    void        reset_steps   (const int, const int, const real=zero);
    void        reset_scheme  (const real, const real=zero, const real=zero);
    void        change_level  (body&, const indx) const;
    void        adjust_level  (body&, const int) const;
    int         longest_moving(const uint)   const;
    void        add_to_level  (const indx l) const { (TS+l)->join (); }
    bool        is_leap_frog  ()             const {return NSTEPS == 1; }
    const uint& number        (const indx l) const {return(TS+l)->N; }
    const real& tau           (const indx l) const {return(TS+l)->TAU; }
    const real& tau           (const body B) const {return(TS+level(B))->TAU; }
    const real& tau_h         (const body B) const {return(TS+level(B))->TAUH; }
    const real& tausq         (const body B) const {return(TS+level(B))->TAUSQ;}
    const real& time          ()             const {return TIME; }
    void        acce_half     (body&)        const;
    void        move_tiny     (body&)        const;
    void        move_by       (body&, const real) const;
    int         highest_level ()             const {return HIGHEST; }
    const real& tau_min       ()             const {return SHORTEST->TAU; }
    const real& tau_max       ()             const {return TS->TAU; }
    const int & h0            ()             const {return H0; }
    const uint& No_in_level   (const indx l) const {return(TS+l)->N; }
    void        short_stats   (std::ostream&,
			       const int=7)  const;
  };
  //MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM//
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW//
  //                                                                          //
  //                      INLINE FUNCTION DEFINITIONS                         //
  //                                                                          //
  //MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM//
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW//

  //////////////////////////////////////////////////////////////////////////////
  // class nbdy::LeapFrog                                                       
  //////////////////////////////////////////////////////////////////////////////
  inline void LeapFrog::predict(const bodies *BODIES) const {
    LoopBodies(bodies,BODIES,B) {
      B.vel().add_mul(acc(B),TAUH);                // v_1/2 = v_0 +a_0 * dt/2   
      B.pos().add_mul(vel(B),TAU );                // x_1   = x_0 +v_1/2 * dt   
    }
    TIME += TAU;                                   // t_1   = t_0 + dt          
  }
  //----------------------------------------------------------------------------
  inline void LeapFrog::accelerate(const bodies *BODIES) const {
    LoopBodies(bodies,BODIES,B)
      B.vel().add_mul(acc(B),TAUH);                // v_1   = v_1/2 * a_1 * dt/2
  }
  //////////////////////////////////////////////////////////////////////////////
  // class nbdy::BlockStep                                                      
  //////////////////////////////////////////////////////////////////////////////
  inline void BlockStep::reset_steps(const int hz, const int Ns, const real t)
  {
    if(TS) delete[] TS;
    NSTEPS   = abs(Ns);
    HIGHEST  = Ns-1;
    MemoryCheck(TS = new TimeStep[NSTEPS]);
    SHORTEST = TS+HIGHEST;
    TIME     = t;
    H0       = hz;
    register int h,i;
    for(i=0, h=H0; i!=NSTEPS; ++i,++h) (TS+i)->re_set(h);
  }
  //----------------------------------------------------------------------------
  inline void BlockStep::reset_counts() const
  {
    for(register TimeStep* s=TS; s<=SHORTEST; s++) s->reset_count();
  }
  //----------------------------------------------------------------------------
  inline bool BlockStep::need_move(const int low) const
  {
    for(register TimeStep* s=TS+low; s<=SHORTEST; s++) if(s->N) return true;
    return false;
  }
  //----------------------------------------------------------------------------
  inline void BlockStep::reset_scheme(const real fa,
				      const real fp,
				      const real fc)
  {
    SCH = 0;
    FAQ = fa*fa; if(FAQ) SCH |= use_a;
    FPQ = fp*fp; if(FPQ) SCH |= use_p;
    FCQ = fc*fc; if(FCQ) SCH |= use_c;
  }
  //----------------------------------------------------------------------------
  inline BlockStep::BlockStep(const int  h0, const int  Ns, const real t,
			      const real fa, const real fp, const real fc)
    : TS(0)
  {
    reset_steps (h0,Ns,t);
    reset_scheme(fa,fp,fc);
  }
  //----------------------------------------------------------------------------
  inline void BlockStep::dump(std::ostream&to)  const
  {
    for(register int i=0; i!=NSTEPS; ++i) {
      to<<i<<": ";
      (TS+i)->dump(to);
    }
  }
  //----------------------------------------------------------------------------
  inline int BlockStep::longest_moving(const uint T) const
  {
    register int i=HIGHEST, t=T;
    for(; !(t&1) && i>0; t>>=1, --i);
    return i;
  }
  //----------------------------------------------------------------------------
  inline void BlockStep::clock_on() const
  {
    TIME += tau_min();
  }
  //----------------------------------------------------------------------------
  inline void BlockStep::acce_half(body&B) const {
    B.vel().add_mul(acc(B),tau_h(B));                // v -> v + h/2 * a        
  }
  //----------------------------------------------------------------------------
  inline void BlockStep::move_tiny(body&B) const {
    B.pos().add_mul(vel(B),tau_min());               // x -> x + h_min * v      
  }
  //----------------------------------------------------------------------------
  inline void BlockStep::move_by(body&B, const real dt) const {
    B.pos().add_mul(vel(B),dt);                      // x -> x + dt * v         
  }
  //----------------------------------------------------------------------------
  inline real BlockStep::step_sq(const body B) const {
    switch(SCH) {
    case use_a:   return FAQ/norm(acc(B));
    case use_p:   return FPQ/square(pot(B));
    case use_c:   return FCQ*abs(pot(B))/norm(acc(B));
    case use_ap:  return min(FAQ/norm(acc(B)), FPQ/square(pot(B)));
    case use_ac:  return min(FAQ, FCQ*abs(pot(B))) / norm(acc(B));
    case use_pc:  return min(FPQ/square(pot(B)),
			     FCQ*abs(pot(B))/norm(acc(B)));
    case use_apc: return min(FPQ/square(pot(B)),
			     FAQ/norm(acc(B)),
			     FPQ/square(pot(B)));
    }
    return zero;
  }
  //----------------------------------------------------------------------------
  inline void BlockStep::change_level(body&B, const indx l) const
  {
    TS[level(B)].leave();                          // leave old level           
    B.level() = l;                                 // set new level             
    TS[level(B)].join ();                          // join new level            
  }
  //----------------------------------------------------------------------------
  inline void BlockStep::to_lower(body&B, const int low) const
  {
    if(level(B) > low) {                           // IF(above lowest active)  >
      TS[level(B)].leave();                        //   leave old level         
      B.level()--;                                 //   set new level           
      TS[level(B)].join();                         //   join new level          
    }                                              // <                         
  }
  //----------------------------------------------------------------------------
  inline void BlockStep::to_higher(body&B) const
  {
    if(level(B) < highest_level()) {               // IF(below highest level)  >
      TS[level(B)].leave();                        //   leave old level         
      B.level()++;                                 //   set new level           
      TS[level(B)].join();                         //   join new level          
    }                                              // <                         
  }
  //----------------------------------------------------------------------------
  inline void BlockStep::adjust_level(body&B, const int low) const
    // implementing scheme 2 of my notes: change only by factors of two         
  {
    const double root_half=0.7071067811865475244;  // sqrt(1/2)                 
    register real tq=step_sq(B)* root_half;        // tau^2 * sqrt(1/2)         
    if     (tausq(B) < tq)    to_lower (B,low);    // change to longer  if poss 
    else if(tausq(B) > tq+tq) to_higher(B);        // change to shorter if poss 
  }
  //----------------------------------------------------------------------------
  inline void BlockStep::short_stats(std::ostream& to, const int w) const
  {
    for(register int l=0; l!=NSTEPS; ++l) to<<std::setw(w)<<number(l)<<" ";
  }
  //----------------------------------------------------------------------------
  inline void BlockStep::assign_level(body&B) const {
    register real t = sqrt(two*step_sq(B));       // tau * 2^(1/2)              
    register indx l = 0;                          // try longest step           
    while(tau(l)          > t  &&                 // WHILE shorter step desired 
	  highest_level() > l) ++l;               // AND possible: decrease     
    add_to_level(B.level()=l);                    // account for in BlockStep   
  }
  //----------------------------------------------------------------------------
}
////////////////////////////////////////////////////////////////////////////////
#endif   // included_step_h

// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// step.h                                                                      |
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
  template<typename bodies_type = bodies>
  class LeapFrog {
    LeapFrog(const LeapFrog&);                     // not implemented           
    LeapFrog operator= (const LeapFrog&);          // not implemented           
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
  protected:
    typedef typename bodies_type::iterator body;   // type of body              
    //--------------------------------------------------------------------------
    // data members                                                             
    //--------------------------------------------------------------------------
    const   real   TAU, TAUH;                      // time step                 
    mutable double TIME;                           // simulation time           
    const   uint REUSE;                            // # re-uses allowed         
    uint         REUSED;                           // # re-uses done            
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    LeapFrog(int    const&f,                       // I: -log_2(tau)            
	     double const&t,                       // I: initial time           
	     uint   const&r) :                     // I: # re-uses allowed      
      TAU (pow(half,f)), TAUH (half*TAU), TIME(t), REUSE(r), REUSED(0u) {}
    //--------------------------------------------------------------------------
    const double &time() const { return TIME; }
    const real   &tau () const { return TAU; }
    //--------------------------------------------------------------------------
    void predict(const bodies_type* const&B) const
    {
      LoopBodies(typename bodies_type,B,Bi) {
	Bi.vel().add_times(acc(Bi),TAUH);
	Bi.pos().add_times(vel(Bi),TAU );
      }
      TIME += double(TAU);
    }
    //--------------------------------------------------------------------------
    void accelerate(const bodies_type* const&B) const
    {
      LoopBodies(typename bodies_type,B,Bi)
	Bi.vel().add_times(acc(Bi),TAUH);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::BlockStep                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename bodies_type = bodies>
  class BlockStep {
    BlockStep(const BlockStep&);                   // not implemented           
    BlockStep& operator= (const BlockStep&);       // not implemented           
    //--------------------------------------------------------------------------
    // private types                                                            
    //--------------------------------------------------------------------------
    typedef typename bodies_type::iterator body;   // type of body              
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
    mutable double TIME;                           // simulation time           
  protected:
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
    void to_lower(body&B, int const&low) const {
      if(level(B) > low) {                         // IF(above lowest active)  >
	TS[level(B)].leave();                      //   leave old level         
	B.level()--;                               //   set new level           
	TS[level(B)].join();                       //   join new level          
      }                                            // <                         
    }
    //--------------------------------------------------------------------------
    void to_higher(body&B) const {
      if(level(B) < highest_level()) {             // IF(below highest level)  >
	TS[level(B)].leave();                      //   leave old level         
	B.level()++;                               //   set new level           
	TS[level(B)].join();                       //   join new level          
      }                                            // <                         
    }
    //--------------------------------------------------------------------------
    void reset_steps(int    const&h0,              // I: h0, tau_max = 2^-h0    
		     int    const&Ns,              // I: #steps                 
		     double const&t =0.)           // I: actual time            
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
    BlockStep    (int    const&h0,                 // I: h0, tau_max = 2^-h0    
		  int    const&Ns,                 // I: #steps                 
		  double const&t)                  // I: actual time            
      : TS ( 0 )
    {
      reset_steps (h0,Ns,t);
    }
    //--------------------------------------------------------------------------
    ~BlockStep                  () { if(TS) delete[] TS;}
    //--------------------------------------------------------------------------
    // other protected methods                                                  
    //--------------------------------------------------------------------------
    void dump(std::ostream&to)  const {
      for(register int i=0; i!=NSTEPS; ++i) {
	to<<i<<": ";
	(TS+i)->dump(to);
      }
    }
    //--------------------------------------------------------------------------
    void clock_on() const {
      TIME += double(tau_min());
    }
    //--------------------------------------------------------------------------
    void reset_counts() const {
      for(register TimeStep* s=TS; s<=SHORTEST; s++) s->reset_count();
    }
    //--------------------------------------------------------------------------
    bool need_move(int const&low) const {
      for(register TimeStep* s=TS+low; s<=SHORTEST; s++) if(s->N) return true;
      return false;
    }
    //--------------------------------------------------------------------------
    void change_level(body const&B, indx const&l) const {
      TS[level(B)].leave();                        // leave old level           
      B.level() = l;                               // set new level             
      TS[level(B)].join ();                        // join new level            
    }
    //--------------------------------------------------------------------------
    int longest_moving(uint const&T) const {
      register int i=HIGHEST, t=T;
      for(; !(t&1) && i>0; t>>=1, --i);
      return i;
    }
    //--------------------------------------------------------------------------
    void          add_to_level(const indx l) const { (TS+l)->join (); }
    bool          is_leap_frog()             const { return NSTEPS == 1; }
    const uint   &number      (const indx l) const { return(TS+l)->N; }
    const real   &tau         (const indx l) const { return(TS+l)->TAU; }
    const double &time        ()             const { return TIME; }
    //--------------------------------------------------------------------------
    const real& tau         (body const&B) const { return(TS+level(B))->TAU; }
    const real& tau_h       (body const&B) const { return(TS+level(B))->TAUH; }
    const real& tausq       (body const&B) const { return(TS+level(B))->TAUSQ;}
    //--------------------------------------------------------------------------
    void half_kick(body&B) const {
      B.vel().add_times(acc(B),tau_h(B));            // v -> v + h/2 * a        
    }
    //--------------------------------------------------------------------------
    void tiny_drift(body&B) const {
      B.pos().add_times(vel(B),tau_min());           // x -> x + h_min * v      
    }
    //--------------------------------------------------------------------------
    void drift_by(body&B, real const&dt) const {
      B.pos().add_times(vel(B),dt);                  // x -> x + dt * v         
    }
    //--------------------------------------------------------------------------
    void acce_half(body&B) const {
      half_kick(B); }
    //--------------------------------------------------------------------------
    void move_tiny(body&B) const {
      tiny_drift(B); }
    //--------------------------------------------------------------------------
    void move_by(body&B, real const&dt) const {
      drift_by(B,dt); }
    //--------------------------------------------------------------------------
    int         highest_level ()             const { return HIGHEST; }
    const real& tau_min       ()             const { return SHORTEST->TAU; }
    const real& tau_max       ()             const { return TS->TAU; }
    const int & h0            ()             const { return H0; }
    const uint& No_in_level   (indx const&l) const { return(TS+l)->N; }
    //--------------------------------------------------------------------------
    void short_stats(std::ostream&to, int const&w=7) const {
      for(register int l=0; l!=NSTEPS; ++l) to<<std::setw(w)<<number(l)<<" ";
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::GravBlockStep                                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename bodies_type = bodies>
  class GravBlockStep : public BlockStep<bodies_type> {
    GravBlockStep(const GravBlockStep&);            // not implemented          
    GravBlockStep& operator= (const GravBlockStep&);// not implemented          
    //--------------------------------------------------------------------------
    // private types                                                            
    //--------------------------------------------------------------------------
  private:
    typedef typename bodies_type::iterator body;   // type of body              
    //--------------------------------------------------------------------------
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
    bool           UPX;                            // use pot+pex for potential 
    int            SCH;                            // stepping scheme           
    real           FAQ,FPQ,FCQ,FEQ;                // factors^2 for stepping    
  protected:
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
    real fpot(body const&B) const
    {
      return UPX? pex(B)+pot(B) : pot(B);
    }
    //--------------------------------------------------------------------------
    real step_sq(body const&B, real const&eps=one)   const
    {
      if(SCH == 0)     return zero;
      if(SCH == use_p) return FPQ/square(fpot(B));
      else {
	register real ia = one/norm(acc(B)), tq=1.e7;
	if(SCH & use_a) update_min(tq, FAQ*ia);
	if(SCH & use_p) update_min(tq, FPQ/square(fpot(B)));
	if(SCH & use_c) update_min(tq, FCQ*abs(fpot(B))*ia);
	if(SCH & use_e) update_min(tq, real(FEQ*eps*sqrt(ia)));
	return tq;
      }
    }
    //--------------------------------------------------------------------------
    void assign_level(body&B, real const&eps=one)
      const {
      register real t = sqrt(two*step_sq(B,eps)); // tau * 2^(1/2)              
      register indx l = 0;                        // try longest step           
      while(tau(l)          > t  &&               // WHILE shorter step desired 
	    highest_level() > l) ++l;             // AND possible: decrease     
      add_to_level(B.level()=l);                  // account for in BlockStep   
    }
    //--------------------------------------------------------------------------
    void adjust_level(body&B, int const&low, real const&eps=one) const {
      // implementing scheme 2 of my notes: change only by factors of two       
      const double root_half=0.7071067811865475244;// sqrt(1/2)                 
      register real tq=step_sq(B,eps)* root_half;  // tau^2 * sqrt(1/2)         
      if     (tausq(B) < tq)    to_lower (B,low);  // change to longer  if poss 
      else if(tausq(B) > tq+tq) to_higher(B);      // change to shorter if poss 
    }
    //--------------------------------------------------------------------------
    void reset_scheme(real const&fa,
		      real const&fp,
		      real const&fc,
		      real const&fe,
		      bool const&up)
    {
      UPX = up;
      SCH = 0;
      FAQ = fa*fa; if(FAQ) SCH |= use_a;
      FPQ = fp*fp; if(FPQ) SCH |= use_p;
      FCQ = fc*fc; if(FCQ) SCH |= use_c;
      FEQ = fe*fe; if(FEQ) SCH |= use_e;
      if(SCH==0) 
	error("[%s.%d]: in GravBlockStep::reset_scheme(%f,%f,%f,%f,%d): "
	      "time step control factors=0",__FILE__,__LINE__,fa,fp,fc,fe,up);
    }
    //--------------------------------------------------------------------------
    // construction & destruction                                               
    //--------------------------------------------------------------------------
    GravBlockStep() : BlockStep<bodies_type>(), UPX(0) {}
    //--------------------------------------------------------------------------
    GravBlockStep(int    const&h0,                 // O: h0, tau_max = 2^-h0    
		  int    const&Ns,                 // I: #steps                 
		  double const&t,                  // I: actual time            
		  real   const&fa,                 // I: f_acc                  
		  real   const&fp,                 // I: f_phi                  
		  real   const&fc,                 // I: f_pa                   
		  real   const&fe,                 // I: f_ea                   
		  bool   const&up)                 // I: use external pot       
      : BlockStep<bodies_type> ( h0,Ns,t )
    {
      reset_scheme(fa,fp,fc,fe,up);
    }
    //--------------------------------------------------------------------------
  };
}
////////////////////////////////////////////////////////////////////////////////
#endif   // falcON_included_step_h

// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// nbdy.h                                                                      |
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
// See the detailled comments for how to use this code. Some abbreviations:    |
//                                                                             |
// R:       return value                                                       |
// I:       function argument is input                                         |
// O:       function argument is output                                        |
// I/O:     function argument is input & output                                |
// [I]:     function argument is optional                                      |
//                                                                             |
// Note that in case of optional function arguments, one may only omit the     |
// last one (and the second but last if the last is already omitted ...).      |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// Note that it is often more convenient to use the wrapper provided by the    |
// classes defined in files yanc.h and ssgs.h                                  |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
//  ON THE TIME-STEPPING SCHEMES                                               |
//  ============================                                               |
//                                                                             |
//  As of February 2002, we have abondoned the old schemes and introduced a    |
//  new scheme instead. It is based on the following prescription (Tremaine):  |
//                                                                             |
//    x  ->  x + h/2 * v    and compute the acceleration a(x)                  |
//    v  ->  v + h/2 * a(x)                                                    |
//    h  ->  h'             based on s(h,h') = tau(x,v,Phi,a,...)              |
//    v  ->  v + h/2 * a(x)                                                    |
//    x  ->  x + h/2 * v                                                       |
//                                                                             |
//  where s(h,h') is a symmetric function that satisfies s(h,h)=h, while tau() |
//  gives the ideal time step. This scheme can also be written as              |
//                                                                             |
//    v  ->  v + h/2 * a(x)                                                    |
//    x  ->  x + h   * v    and compute the acceleration a(x)                  |
//    h  ->  h'                                                                |
//                                                                             |
//  Currently, we are allowing for the following options                       |
//                                                                             |
//    tau = f_a / |a|                (as in GADGET, Springel et al)            |
//    tau = f_p / |Phi|                                                        |
//    tau = f_c * sqrt(|Phi|) /|a|   (as in PKDGRAV, Stadel)                   |
//    tau = f_e * sqrt(eps/|a|)      (as in PKDGRAV, Stadel)                   |
//                                                                             |
//  and combinations of these (whereby the minimum tau is used).               |
//                                                                             |
//  Note that with this method the time step h may not be a monotonic function |
//  of time even when tau is. In particular, if tau=const but h<tau, we have   |
//  an oscillation between a value h<tau and a value h'>tau.                   |
//                                                                             |
//                                                                             |
//  Using the Block Step Scheme                                                |
//  ---------------------------                                                |
//                                                                             |
//  In our application of the above to the blockstep scheme we specify s(h,h') |
//  to be the geometric mean and replace tau by the integer power of sqrt(2)   |
//  nearest (in the log) to tau. This results in the following scheme.         |
//                                                                             |
//  h' = h/2  if  h > tau*2^( 1/4)                                             |
//  h' = 2 h  if  h < tau*2^(-1/4) and  h'= 2 h is possible with the blockstep |
//  h' = h    otherwise                                                        |
//                                                                             |
//  Note that our restriction to change the time step by no more than a factor |
//  of two results in violation of the time reversibility. However, for        |
//  carefully chosen time steps, the necessity to change by a factor of 4 or   |
//  more should hardly ever occur.                                             |
//                                                                             |
//  For the initial h we take the integer power of 2 that is nearest to tau.   |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_nbdy_h
#define falcON_included_nbdy_h

#ifndef falcON_included_ctime
#  include <ctime>
#  define falcON_included_ctime
#endif
#ifndef falcON_included_body_h
#  include <body.h>
#endif
#ifndef falcON_included_falcON_h
#  include <falcON.h>
#endif
#ifndef falcON_included_pext_h
#  include <pext.h>
#endif
#ifndef falcON_included_step_h
#  include <public/step.h>
#endif
#if defined(falcON_NEMO) && !defined(falcON_included_nmio_h)
#  include <public/nmio.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace meta {
  using namespace nbdy;
  //////////////////////////////////////////////////////////////////////////////
  template<int N, int I=0> struct __tr {
    template<typename scalar> static scalar a(scalar t[N+1][N+1]) {
      return t[I][I] + __tr<N,I+1>::a(t); } };
  template<int N> struct __tr<N,N> {
    template<typename scalar> static scalar a(scalar t[N+1][N+1]) {
      return t[N][N]; } };
  //////////////////////////////////////////////////////////////////////////////
  template<int N, int I=0, int J=0> struct __addt {
    template<typename scalar> static 
    void a(scalar t[N+1][N+1], const scalar*x, const real*y) {
      t[I][J] += x[I] * y[J];
      __addt<N,I,J+1>::a(t,x,y); } };
  template<int N,int I> struct __addt<N,I,N> {
    template<typename scalar> static 
    void a(scalar t[N+1][N+1], const scalar*x, const real*y) {
      t[I][N] += x[I] * y[N];
      __addt<N,I+1,0>::a(t,x,y); } };
  template<int N> struct __addt<N,N,N> {
    template<typename scalar> static 
    void a(scalar t[N+1][N+1], const scalar*x, const real*y) {
      t[N][N] += x[N] * y[N]; } };
  //////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::nbody_base                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class nbody_base : public falcON
  {
    nbody_base(const nbody_base&);                 // not implemented           
    nbody_base& operator= (const nbody_base&);     // not implemented           
    //--------------------------------------------------------------------------
    // type tensor                                                              
    //--------------------------------------------------------------------------
  protected:
    typedef real tensor[Ndim][Ndim];               // not necessarily symmetric 
    //--------------------------------------------------------------------------
    // trace of tensor                                                          
    //--------------------------------------------------------------------------
    static double tr(double t[Ndim][Ndim]) { return meta::__tr<Ndim-1>::a(t); }
    static float  tr(float  t[Ndim][Ndim]) { return meta::__tr<Ndim-1>::a(t); }
    //--------------------------------------------------------------------------
    // tensor computation as outer product                                      
    //--------------------------------------------------------------------------
    static void AddTensor(double t[Ndim][Ndim], vect_d const&x, vect const&y) {
      meta::__addt<Ndim-1>::a(t,
			      static_cast<const double*> (x),
			      static_cast<const real  *> (y));
    }
    static void AddTensor(real t[Ndim][Ndim], vect const&x, vect const&y) {
      meta::__addt<Ndim-1>::a(t,
			      static_cast<const real*> (x),
			      static_cast<const real*> (y));
    }
    //--------------------------------------------------------------------------
    // record a CPU timing                                                      
    //--------------------------------------------------------------------------
    static void record_cpu(clock_t& c0, real& CPU) {
      register clock_t c1 = clock();
      CPU += (c1-c0)/real(CLOCKS_PER_SEC);
      c0 = c1;
    }
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
    mutable clock_t C_OLD;
    const   sbodies*BODIES;                        // bodies                    
    const   bool    SELF_GRAV;                     // is G != 0?                
    const   vect   *ROOTCENTER;                    // pre-determined root center
    const   gravity*PEX;                           // external potential if any 
    const   real    TINI;                          // initial simulation time   
    const   int     NCRIT;                         // max # bodies in cell      
    mutable real    CPU_BUILD,                     // time for falcON::grow etc 
                    CPU_GRAV,                      // time for gravity etc      
                    CPU_PEX,                       // time for ext potential    
                    CPU_STEP,                      // total time for a longstep 
                    CPU_TOTAL;                     // total time so far         
    mutable real    M,T,Vin,Vex,W,TW;              // mass, kin & pot E         
    mutable amom    L;                             // total angular momentum    
    mutable tensor  KT,WT;                         // kin & pot energy          
    mutable vect    CMX,CMV;                       // center of mass pos & vel  
    mutable bool    DIAG;                          // diagnose up-to-date ?     
    //--------------------------------------------------------------------------
    // abstract protected methods                                               
    //--------------------------------------------------------------------------
    virtual void diagnose  () const=0;             // update E,T,V,W etc        
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
    void reset_cpus() const
    {
      CPU_BUILD = zero;
      CPU_GRAV  = zero;
      CPU_PEX   = zero;
      CPU_STEP  = zero;
    }
    //--------------------------------------------------------------------------
    void update_diags() const {
      if(! DIAG) {
	diagnose();
	DIAG = true;
      }
    }
    //--------------------------------------------------------------------------
    void update_cpu_total() const
    {
      register clock_t C_NEW = clock();
      CPU_TOTAL += (C_NEW-C_OLD)/real(CLOCKS_PER_SEC);
      C_OLD = C_NEW;
    }
    //--------------------------------------------------------------------------
    nbody_base (const sbodies*const&b,             // I: bodies                 
		real          const&e,             // I: eps                    
		real          const&ti,            // I: initial time           
		real          const&th,            // I: tolerance parameter    
		int           const&nc,            // I: N_crit                 
		const vect*   const&x0,            // I: pre-set root center    
		kern_type     const&ke,            // I: softening kernel       
		real          const&g,             // I: Newton's G             
#ifdef falcON_INDI
		bool          const&is,            // I: individual softening?  
#endif
		const gravity*const&px,            // I: P_ex                   
		const int           gd[4]          // I: direct sum: gravity    
#ifdef falcON_SPH
		,
		const int           sd[3]          // I: direct sum: SPH        
		= Default::SPHdirect
#endif
		) : 
      falcON    ( b, abs(e) , abs(th), ke, 
#ifdef falcON_INDI
		  is,
#endif
		  g, th<0? const_theta : theta_of_M, gd
#ifdef falcON_SPH
		  ,sd
#endif
		  ),
      BODIES    ( b ),
      SELF_GRAV ( g != zero ),
      ROOTCENTER( x0 ),
      PEX       ( px ),
      TINI      ( ti ),
      C_OLD     ( clock() ), 
      NCRIT     ( max(1,nc) ),
      CPU_TOTAL ( zero ), 
      DIAG      ( false )
    {
      if(PEX != 0 && !BODIES->has(io::q))
	falcON_ErrorF("external potential desired, but nobody has memory",
		      "nbody_base");
    }
    //--------------------------------------------------------------------------
    // abstract public methods                                                  
    //--------------------------------------------------------------------------
  public:
    virtual real const&time() const=0;             // actual time               
    virtual void full_step () = 0;                 // one full time step        
    virtual void stats     (std::ostream&) const=0;// statistics output         
    virtual void stats_head(std::ostream&) const=0;// header for stats output   
    virtual void stats_line(std::ostream&) const=0;// line for stats output     
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
    real const&t_ini          () const { return TINI; }
    int  const&Ncrit          () const { return NCRIT; }
    //--------------------------------------------------------------------------
    real const&kin_energy     () const { update_diags(); return T; }
    real const&pot_self_energy() const { update_diags(); return Vin; }
    real const&pot_ext_energy () const { update_diags(); return Vex; }
    real const pot_energy     () const { update_diags(); return Vin+Vex; }
    real const pot_energy_acc () const { update_diags(); return W; }
    real       total_energy   () const { update_diags(); return T+Vin+Vex;}
    real const&virial_ratio   () const { update_diags(); return TW; }
    amom const&total_angmom   () const { update_diags(); return L; }
    vect const&total_momentum () const { update_diags(); return CMV; }
    vect const&center_of_mass () const { update_diags(); return CMX; }
    real const&kin_energy     (int i, int j) const { update_diags();
                                                     return KT[i][j]; }
    real const&pot_energy     (int i, int j) const { update_diags();
                                                     return WT[i][j]; }
    bool       using_extpot   () const { return PEX && !PEX->is_empty(); }
    //--------------------------------------------------------------------------
    real       as_angmom     (int const&i, int const&j) const 
    {
      update_diags();
#if falcON_NDIM == 3
      switch(i) {
      case 0:
	switch(j) {
	case 0:  return  zero;
	case 1:  return  L[2];
	default: return -L[1];
	}  
      case 1:
	switch(j) {
	case 0:  return -L[2];
	case 1:  return  zero;
	default: return  L[0];
	}  
      default:
	switch(j) {
	case 0:  return  L[1];
	case 1:  return -L[0];
	default: return  zero;
	}  
      }
#else
      return (i==j) ? zero : real(L);
#endif
    }
    //--------------------------------------------------------------------------
    real const&cpu_build       () const { return CPU_BUILD; }
    real const&cpu_self_grav   () const { return CPU_GRAV; }
    real const&cpu_ext_grav    () const { return CPU_PEX; }
    real       cpu_grav        () const { return CPU_GRAV+CPU_PEX; }
    real const&cpu_longstep    () const { return CPU_STEP; }
    real const&cpu_total       () const { return CPU_TOTAL; }
    //--------------------------------------------------------------------------
    void  reset_cpu_total       (const real sec) const { CPU_TOTAL=sec; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::basic_nbody                                                  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class basic_nbody : public nbody_base
  {
    basic_nbody(const basic_nbody&);               // not implemented           
    basic_nbody& operator= (const basic_nbody&);   // not implemented           
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
#ifdef falcON_INDI
  public:
    enum soft_type {                               // type of softening         
      global_fixed        = 0,                     // globally time-constant    
      individual_fixed    = 1,                     // individual but time-const 
#ifdef falcON_ADAP
      individual_adaptive = 2                      // individual and adaptive   
#endif
    };
#endif
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
  protected:
#ifdef falcON_INDI
    soft_type       SOFTENING;                     // type of softening method  
#endif
  private:
#ifdef falcON_ADAP
    real            NSOFT;                         // #/sphere for eps_i setting
    uint            NREF;                          // #/cell for eps_i setting  
    real            EMIN;                          // lower limit for eps_i     
    real            EFAC;                          // max change factor for epsi
#endif
  protected:
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
    // - compute KT, WT, CMX, CMV, T, Vin, Vex, W                               
    void diagnose          () const;
    //--------------------------------------------------------------------------
    // - IF Grav != 0:                                                          
    //   - (re-)grows or re-uses the tree (the latter if arg1 = true)           
    //   - optionally adjusts eps_i of active bodies  (individual_adaptive)     
    //   - upates self-gravity of active bodies                                 
    //   ELSE (Grav == 0):                                                      
    //   - reset accelerations & potentials of active bodies                    
    // - add external gravity                                                   
    void set_gravity    (bool const&,              // I: re-use old tree?       
			 bool const&               // I: all or only active?    
#ifdef falcON_ADAP
			,bool const& =true         //[I: indi&adap: adjust epsi]
#endif
			);
    //--------------------------------------------------------------------------
#ifdef falcON_ADAP
    void estimate_mass_density(const bool);
#endif
    //--------------------------------------------------------------------------
    void stats_front     (std::ostream&) const;    // front part of stats       
    void stats_back      (std::ostream&) const;    // back  part of stats       
    void stats_head_front(std::ostream&) const;    // front part of stats_head  
    void stats_head_back (std::ostream&) const;    // back  part of stats_head  
    void stats_line_front(std::ostream&) const;    // front part of stats_line  
    void stats_line_back (std::ostream&) const;    // back  part of stats_line  
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
#ifdef falcON_ADAP
    void estimate_mass_densities(                  // estimate every bodies rho 
				 const bool=true); //[I: add contrib of Pot_ext]
    void estimate_surf_densities();                // estimate every bodies SD  
#endif
    //--------------------------------------------------------------------------
    // - sets Nsoft, Nref, eps, kernel, but NOT soft_type!                      
    void reset_softening(                          // resets softening params   
			 kern_type const&,         // I: softening kernel       
			 real      const&          // I: eps                    
#ifdef falcON_ADAP
			,real      const&,         // I: Nsoft                  
			 unsigned  const&,         // I: Nref                   
			 real      const&,         // I: eps_min                
			 real      const&          // I: eps_fac                
#endif
			 );
    //--------------------------------------------------------------------------
    void reset_opening(const real) const;
    //--------------------------------------------------------------------------
    inline
    basic_nbody(const sbodies*const&,              // I: bodies                 
		real          const&,              // I: eps/eps_max            
		real          const&,              // I: initial time           
		real          const&,              // I: tolerance parameter    
		int           const&,              // I: N_crit                 
		const vect*   const&,              // I: pre-set root center    
		kern_type     const&,              // I: softening kernel       
		real          const&,              // I: Newton's G             
#ifdef falcON_INDI
		soft_type     const&,              // I: softening type         
#ifdef falcON_ADAP
		real          const&,              // I: N_soft                 
		uint          const&,              // I: N_ref                  
		real          const&,              // I: eps_min                
		real          const&,              // I: eps_fac                
#endif
#endif
		const gravity*const&,              // I: P_ex                   
		const int[4]                       // I: direct sum: gravity    
#ifdef falcON_SPH
		,
		const int[3] = Default::SPHdirect  // I: direct sum: SPH        
#endif
		);
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::LeapFrogCode                                                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class LeapFrogCode :
    public    basic_nbody,
    protected LeapFrog<sbodies>
  {
    LeapFrogCode(const LeapFrogCode&);             // not implemented           
    LeapFrogCode& operator= (const LeapFrogCode&); // not implemented           
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    LeapFrogCode(const sbodies*const&,             // I: bodies                 
		 real          const&,             // I: eps/eps_max            
		 real          const&,             // I: initial time           
		 int           const&,             // I: h0                     
		 int           const&,             // I: h_grow                 
		 real          const&,             // I: tolerance param        
		 int           const&,             // I: N_crit                 
		 const vect*   const&,             // I: root center            
		 kern_type     const&,             // I: softening kernel       
		 real          const&,             // I: Newton's G             
#ifdef falcON_INDI
		 soft_type     const&,             // I: softening type         
#ifdef falcON_ADAP
		 real          const&,             // I: N_soft                 
		 uint          const&,             // I: N_ref                  
		 real          const&,             // I: eps_min                
		 real          const&,             // I: eps_fac                
#endif
#endif
		 const gravity*const&,             // I: P_ex                   
		 const         int[4]);            // I: direct sum control     
    //--------------------------------------------------------------------------
    void full_step       ();                       // a single leap-frog step   
    const real& time() const { return LeapFrog<sbodies>::time(); }
    void stats     (std::ostream&) const;          // statistic of long_step    
    void stats_head(std::ostream&) const;          // table head for short_stats
    void stats_line(std::ostream&) const;          // line of proper size       
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::BlockStepCode                                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class BlockStepCode :
    public    basic_nbody,
    protected GravBlockStep<sbodies>
  {
    BlockStepCode(const BlockStepCode&);            // not implemented          
    BlockStepCode& operator= (const BlockStepCode&);// not implemented          
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
    const indx LGROW;
    int   W;                                       // width for output #levels  
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
  private:
    void set_active_flags(int const&);
    void elementary_step (int const&);
    void prepare         (int const&,              // I: h0                     
			  int const&);             // I: N_steps                
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    BlockStepCode(const sbodies*const&,            // I: bodies                 
		  real          const&,            // I: eps/eps_max            
		  real          const&,            // I: initial time           
		  int           const&,            // I: h0                     
		  int           const&,            // I: # levels               
		  real          const&,            // I: f_a: for stepping      
		  real          const&,            // I: f_p: for stepping      
		  real          const&,            // I: f_c: for stepping      
		  real          const&,            // I: f_e: for stepping      
		  int           const&,            // I: h_grow                 
		  real          const&,            // I: tolerance param        
		  int           const&,            // I: N_crit                 
		  const vect*   const&,            // I: root center            
		  kern_type     const&,            // I: softening kernel       
		  real          const&,            // I: Newton's G             
#ifdef falcON_INDI
		  soft_type     const&,            // I: softening type         
#ifdef falcON_ADAP
		  real          const&,            // I: N_soft                 
		  uint          const&,            // I: N_ref                  
		  real          const&,            // I: eps_min                
		  real          const&,            // I: eps_fac                
#endif
#endif
		  const gravity*const&,            // I: P_external             
		  const         int[4]);           // I: direct sum control     
    //--------------------------------------------------------------------------
    void reset_stepping (real const&,              // I: f_a                    
			 real const& =zero,        //[I: f_p]                   
			 real const& =zero,        //[I: f_c]                   
			 real const& =zero);       //[I: f_e]                   
    void full_step      ();                        // do one blockstep          
    //--------------------------------------------------------------------------
    const real& time    () const {
      return BlockStep<sbodies>::time(); }
    //--------------------------------------------------------------------------
    void dump_steps     (std::ostream& to) const {
      BlockStep<sbodies>::dump(to); }
    //--------------------------------------------------------------------------
    void stats     (std::ostream&) const;          // statistic of long_step    
    void stats_head(std::ostream&) const;          // statistic of long_step    
    void stats_line(std::ostream&) const;          // line of proper size       
  };
#ifdef falcON_NEMO
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::NbodyCode                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class NbodyCode {
    NbodyCode           (NbodyCode const&);        // not implemented           
    NbodyCode& operator=(NbodyCode const&);        // not implemented           
    //--------------------------------------------------------------------------
  private:
    real         const EFAC;                       // max factor to change eps_i
    int          const ND;                         // # NEMO devices            
    io           const SUPPORT, INPUT;             // what: supported & given   
    nemo_out*    const NEMO;                       // NEMO output devices       
    sbodies*     const BODIES;                     // bodies themselves         
    std::string  const IFILE;                      // input file                
    basic_nbody*       NBODY;                      // actual N-body code        
    //--------------------------------------------------------------------------
    // construction                                                             
    //--------------------------------------------------------------------------
  public:
    NbodyCode(// data input                                                     
	      const char         *,                //   I: input file           
	      bool          const&,                //   I: resume old (if nemo) 
	      // time-integration related                                       
	      int           const&,                //   I: hmin: t_min=2^(-hmin)
	      int           const&,                //   I: # time-step levels   
	      real          const&,                //   I: f_a                  
	      real          const&,                //   I: f_p                  
	      real          const&,                //   I: f_c                  
	      real          const&,                //   I: f_e                  
	      // tree related                                                   
	      int           const&,                //   I: N_crit               
	      int           const&,                //   I: h_grow               
	      const vect   *const&,                //   I: pre-set root center  
	      // gravity related                                                
	      real          const&,                //   I: eps                  
	      kern_type     const&,                //   I: softening kernel     
	      const gravity*const&,                //   I: P_ex                 
	      real          const&,                //   I: tolerance parameter  
	      real          const&,                //   I: Newton's G           
	      // default arguments                                              
#ifdef falcON_INDI
#ifdef falcON_ADAP
	      real      const& = zero,             //  [I: Nsoft]               
	      uint      const& = 32,               //  [I: Nref]                
	      real      const& = zero,             //  [I: emin]                
#endif
	      basic_nbody::soft_type
	      const& =basic_nbody::global_fixed,   //  [I: softening type]      
#endif
	      int       const& = 2,                //  [I: # nemo output devs]  
	      io        const& = io::o,            //  [I: what else to read?]  
	      const int[4] = Default::direct);     //  [I: direct sum: gravity] 
    //--------------------------------------------------------------------------
    // destruction                                                              
    //--------------------------------------------------------------------------
    ~NbodyCode() {                                 // destructor                
      delete BODIES;
      delete NBODY;
    }
    //--------------------------------------------------------------------------
    // nemo outputs                                                             
    //--------------------------------------------------------------------------
    void  describe_nemo(                           // describe simulation       
			std::ostream&,             // I: output stream          
			const char*);              // I: command line           
    //--------------------------------------------------------------------------
    void  open_nemo    (                           // open NEMO output stream   
			int const&  =0,            //[I: index of nemo stream]  
			const char* =0,            //[I: file name ]            
			bool const& =0);           //[I" resume old sim?]       
    //--------------------------------------------------------------------------
    void  close_nemo   (int const&d=0) {           //[I: index of nemo stream]  
      if(d<ND) (NEMO+d)->close(); }
    //--------------------------------------------------------------------------
    void  write_nemo   (io   const&,               //[I: what to write out]     
			int  const& =0,            //[I: index of nemo stream]  
			bool const& =true);        //[I: output diagnostics?]   
    //--------------------------------------------------------------------------
    int   nemo_devices () const {                  // R: # nemo output streams  
      return ND; }
    //--------------------------------------------------------------------------
    bool  nemo_is_open (int const&d=0) const {     //[I: index of nemo stream]  
      return d < ND && (NEMO+d)->is_open(); }
    //--------------------------------------------------------------------------
    // time integration                                                         
    //--------------------------------------------------------------------------
    void  full_step    () { NBODY->full_step(); }  // perform one full time step
    //--------------------------------------------------------------------------
    // statistic outputs                                                        
    //--------------------------------------------------------------------------
    void  stats       (std::ostream&o) const {
      o<<' '; NBODY->stats(o); }
    //--------------------------------------------------------------------------
    void  stats_head  (std::ostream&o) const {
      o<<'#'; NBODY->stats_head(o);
      o<<'#'; NBODY->stats_line(o); }
    //--------------------------------------------------------------------------
    // data access                                                              
    //--------------------------------------------------------------------------
    int         const&No_nemo_devices() const { return ND; }
    real        const&initial_time   () const { return NBODY->t_ini(); }
    sbodies    *const&bodies         () const { return BODIES; }
    std::string const&input_file     () const { return IFILE; }
    real        const&time           () const { return NBODY->time(); }
    real        const&Grav           () const { return NBODY->NewtonsG(); }
    bool              okay           () const { return BODIES!=0 && NBODY!=0; }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_NEMO               
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_nbdy_h    


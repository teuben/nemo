// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// nbdy.h                                                                      |
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
//    tau = f_c * sqrt(|Phi|) / a    (as in PKDGRAV, Stadel)                   |
//                                                                             |
//  and combinations of these (whereby the minimum tau is used).               |
//                                                                             |
//  Note that with this method the time step h may not be a monotonic function |
//  of time even when tau is. In particular, if tau=const but h<tau, we have   |
//  an oscillation between a value h<tau and a value h'>tau.                   |
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
//  more should hardly occur.                                                  |
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
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::basic_nbody                                                  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class basic_nbody : protected falcON
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
      individual_adaptive = 2                      // individual and adaptive   
    };
#endif
    //--------------------------------------------------------------------------
  protected:
    typedef real tensor[Ndim][Ndim];               // not necessarily symmetric 
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
    const bodies   *BODIES;
    const extpot   *PEX;
#ifdef falcON_INDI
    soft_type       SOFTENING;                     // type of softening method  
#endif
    real            TINI;                          // initial simulation time   
  private:
    mutable clock_t C_OLD;
#ifdef falcON_INDI
    real            NSOFT;                         // #/sphere for eps_i setting
    uint            NREF;                          // #/cell for eps_i setting  
    real            EMIN;                          // lower limit for eps_i     
    real            EFAC;                          // max change factor for epsi
#endif
    const int       NCRIT;                         // max # bodies in cell      
    const int       NCUT;                          // for tree re-build         
    int             DIR[4];                        // direct summation          
    mutable real    M,Ktot,Uin,Uex,TU;             // mass, kin & pot E, vir rat
    mutable amom    L;                             // total angular momentum    
    mutable tensor  KT,WT;                         // kin & pot energy, AM tens 
    mutable vect    CMX,CMV;                       // center of mass pos & vel  
  protected:
    mutable bool    DIAG;                          // diagnose up-to-date ?     
    mutable real    CPU_BUILD,                     // time for falcON::grow etc 
                    CPU_GRAV,                      // time for gravity etc      
                    CPU_PEX,                       // time for ext potential    
                    CPU_DENS,                      // time for density estimates
                    CPU_STEP,                      // total time for a longstep 
                    CPU_TOTAL;                     // total time so far         
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
  protected:
    //--------------------------------------------------------------------------
    // - compute KT, WT, CMX, CMV, Ktot, Uin, Uex, TU                           
    void diagnose          () const;
    //--------------------------------------------------------------------------
    // - IF out of date: compute KT, WT, CMX, CMV, Ktot, Uin, Uex, TU           
    void update_diags() const { if(! DIAG) { diagnose(); DIAG = true; } }
    //--------------------------------------------------------------------------
    // - (re-)grows or re-uses the tree (the latter if arg = true)              
    void set_tree          (const bool);
    //--------------------------------------------------------------------------
    // - adjusts eps_i of active bodies       (individual_adaptive)             
    // - upates acceleration of active bodies                                   
    void eps_and_acc    (bool const&               // I: all or only active?    
#ifdef falcON_INDI
			,bool const& =true         //[I: indi&adap: adjust epsi]
#endif
			);
    //--------------------------------------------------------------------------
    inline void reset_cpus           () const;
    void update_cpu_total     () const;
#ifdef falcON_INDI
    void estimate_mass_density(const bool);
#endif
    //--------------------------------------------------------------------------
    // abstract methods                                                         
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
#ifdef falcON_INDI
    void estimate_mass_densities(                  // estimate every bodies rho 
				 const bool=true); //[I: add contrib of Pot_ext]
    void estimate_surf_densities();                // estimate every bodies SD  
#endif
    //--------------------------------------------------------------------------
    real const&kin_energy     () const { update_diags(); return Ktot; }
    real const&pot_self_energy() const { update_diags(); return Uin; }
    real const&pot_ext_energy () const { update_diags(); return Uex; }
    real const pot_energy     () const { update_diags(); return Uin+Uex; }
    real       total_energy   () const { update_diags(); return Ktot+Uin+Uex;}
    real const&virial_ratio   () const { update_diags(); return TU; }
    amom const&total_angmom   () const { update_diags(); return L; }
    vect const&total_momentum () const { update_diags(); return CMV; }
    vect const&center_of_mass () const { update_diags(); return CMX; }
    real const&kin_energy     (int const&i, int const&j) const {
                                         update_diags(); return KT[i][j]; }
    real const&pot_energy     (int const&i, int const&j) const {
                                         update_diags(); return WT[i][j]; }
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
    real const&cpu_tree_force  () const { return CPU_GRAV; }
    real const&cpu_ext_force   () const { return CPU_PEX; }
    real       cpu_force       () const { return CPU_GRAV+CPU_PEX; }
    real const&cpu_density     () const { return CPU_DENS; }
    real const&cpu_longstep    () const { return CPU_STEP; }
    real const&cpu_total       () const { return CPU_TOTAL; }
    //--------------------------------------------------------------------------
    void  reset_cpu_total       (const real sec) const { CPU_TOTAL=sec; }
    //--------------------------------------------------------------------------
    // - sets Nsoft, Nref, eps, kernel, but NOT soft_type!                      
    void reset_softening(                          // resets softening params   
			 kern_type const&,         // I: softening kernel       
			 real      const&          // I: eps                    
#ifdef falcON_INDI
			,real      const&,         // I: Nsoft                  
			 unsigned  const&,         // I: Nref                   
			 real      const&,         // I: eps_min                
			 real      const&          // I: eps_fac                
#endif
			 );
    //--------------------------------------------------------------------------
    void reset_opening(const real) const;
    //--------------------------------------------------------------------------
    const real t_ini() const { return TINI; }
    //--------------------------------------------------------------------------
    basic_nbody(const bodies*const&,               // I: bodies                 
		real         const&,               // I: eps/eps_max            
		real         const&,               // I: initial time           
		real         const&,               // I: tolerance parameter    
		int          const&,               // I: N_crit                 
		kern_type    const&,               // I: softening kernel       
#ifdef falcON_INDI
		soft_type    const&,               // I: softening type         
		real         const&,               // I: N_soft                 
		uint         const&,               // I: N_ref                  
		real         const&,               // I: eps_min                
		real         const&,               // I: eps_fac                
#endif
		const extpot*const&,               // I: P_ex                   
		const int[4]           );          // I: direct sum control     
    //--------------------------------------------------------------------------
    virtual ~basic_nbody() {}
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::LeapFrogCode                                                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class LeapFrogCode :
    public    basic_nbody,
    protected LeapFrog
  {
    LeapFrogCode(const LeapFrogCode&);              // not implemented          
    LeapFrogCode& operator= (const LeapFrogCode&);  // not implemented          
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
    const uint REUSE;
          uint REUSED;
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    LeapFrogCode(const bodies*const&,                  // I: bodies             
		 real         const&,                  // I: eps/eps_max        
		 real         const&,                  // I: initial time       
		 int          const&,                  // I: h0                 
		 int          const& =0,               //[I: h_grow]            
		 real         const& =Default::theta,  //[I: tolerance param]   
		 int          const& =Default::Ncrit,  //[I: N_crit]            
		 kern_type    const& =Default::kernel, //[I: softening kernel]  
#ifdef falcON_INDI
		 soft_type    const& =global_fixed,    //[I: softening type]    
		 real         const& =zero,            //[I: N_soft]            
		 uint         const& =32,              //[I: N_ref]             
		 real         const& =zero,            //[I: eps_min]           
		 real         const& =two,             //[I: eps_fac]           
#endif
		 const extpot*const& =0,               //[I: P_ex]              
		 const int[4]        =Default::direct);//[I: direct sum control]
    //--------------------------------------------------------------------------
    void full_step       ();                       // a single leap-frog step   
    const real& time() const { return LeapFrog::time(); }
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
    protected GravBlockStep
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
    BlockStepCode(const bodies*const&,                 // I: bodies             
		  real         const&,                 // I: eps/eps_max        
		  real         const&,                 // I: initial time       
		  int          const&,                 // I: h0                 
		  int          const&,                 // I: # levels           
		  real         const&,                 // I: f_a: for stepping  
		  real         const& =zero,           //[I: f_p: for stepping] 
		  real         const& =zero,           //[I: f_c: for stepping] 
		  real         const& =zero,           //[I: f_e: for stepping] 
		  int          const& =0,              //[I: h_grow]            
		  real         const& =Default::theta, //[I: tolerance param]   
		  int          const& =Default::Ncrit, //[I: N_crit]            
		  kern_type    const& =Default::kernel,//[I: softening kernel]  
#ifdef falcON_INDI
		  soft_type    const& =global_fixed,   //[I: softening type]    
		  real         const& =zero,           //[I: N_soft]            
		  uint         const& =32,             //[I: N_ref]             
		  real         const& =zero,           //[I: eps_min]           
		  real         const& =two,            //[I: eps_fac]           
#endif
		  const extpot*const& =0,              //[I: P_ex]              
		  const int[4]       =Default::direct);//[I: direct sum control]
    //--------------------------------------------------------------------------
    void reset_stepping (real const&,              // I: f_a                    
			 real const& =zero,        //[I: f_p]                   
			 real const& =zero,        //[I: f_c]                   
			 real const& =zero);       //[I: f_e]                   
    void full_step      ();                        // do one blockstep          
    const real& time    () const { return BlockStep::time(); }
    void dump_steps     (std::ostream& to) const { BlockStep::dump(to); }
    void stats     (std::ostream&) const;          // statistic of long_step    
    void stats_head(std::ostream&) const;          // statistic of long_step    
    void stats_line(std::ostream&) const;          // line of proper size       
  };
}
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_nbdy_h

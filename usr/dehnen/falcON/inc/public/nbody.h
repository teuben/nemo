// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/nbody.h                                                  
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2000-2008                                                           
///                                                                             
/// \brief  provides classes for the efficient implementation of N-body codes   
///                                                                             
/// \todo   complete doxygen documentation                                      
///                                                                             
/// Here, we disentangle the integration from force computation and diagnostics.
/// The latter two are so tightly connected that it makes no sense to force them
/// to divorce.                                                               \n
/// Any N-body code fitting in our scheme can be put together with class        
/// NBodyCode. For practical applications see class FalcONCode below or         
/// class PotExpCode in pexp.h (proprietary only).                              
///                                                                             
/// \version 14-01-2005 WD integration of H,R moved Integrator -> ForceSPH      
/// \version 01-04-2005 WD soft_type -> enum.h as falcON::soft_type             
///                     WD last 4 args of NbodyCode::init() default to          
///                        fieldset::empty                                      
/// \version 05-04-2005 WD SELF_GRAV removed from ForceALCON (in ForceDiagGrav) 
///                     WD added BlockStepCode::StepLevels::always_adjust() to  
///                        allow for the new schemes in DirectCode.cc           
///                     WD increase precision of log output if real==double     
/// \version 13-07-2005 WD adapt for new falcON                                 
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2000-2008 Walter Dehnen                                        
//                                                                              
// This program is free software; you can redistribute it and/or modify         
// it under the terms of the GNU General Public License as published by         
// the Free Software Foundation; either version 2 of the License, or (at        
// your option) any later version.                                              
//                                                                              
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//                                                                              
// You should have received a copy of the GNU General Public License            
// along with this program; if not, write to the Free Software                  
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#ifndef falcON_included_nbody_h
#define falcON_included_nbody_h

#ifndef falcON_included_body_h
#  include <body.h>
#endif
#ifndef falcON_included_externacc_h
#  include <externacc.h>
#endif
#ifndef falcON_included_forces_h
#  include <forces.h>
#endif
#ifndef falcON_included_iomanip
#  include <iomanip>
#  define falcON_included_iomanip
#endif
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  namespace meta {
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
    void a(scalar t[N+1][N+1], const scalar*x, const falcON::real*y) {
      t[I][J] += x[I] * y[J];
      __addt<N,I,J+1>::a(t,x,y); } };
  template<int N,int I> struct __addt<N,I,N> {
    template<typename scalar> static 
    void a(scalar t[N+1][N+1], const scalar*x, const falcON::real*y) {
      t[I][N] += x[I] * y[N];
      __addt<N,I+1,0>::a(t,x,y); } };
  template<int N> struct __addt<N,N,N> {
    template<typename scalar> static 
    void a(scalar t[N+1][N+1], const scalar*x, const falcON::real*y) {
      t[N][N] += x[N] * y[N]; } };
  } // namespace meta {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::ForceAndDiagnose                                             
  //                                                                            
  /// abstract base class for N-body force computation and diagnosis            
  /// provides some useful tools related to diagnosis                           
  ///                                                                           
  /// \note      the force solver may also need to assist in the integration:   
  /// evolution of quantities not fitting in the integration scheme below MUST  
  /// be dealt with by the force solver, for instance, adaption of the softening
  /// or smoothing lengths (which really is a force rather than an integration  
  /// issue).                                                                   
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class ForceAndDiagnose {
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
    snapshot*          const SNAPSHOT;
    const acceleration*const ACCEXTERN;
    //--------------------------------------------------------------------------
    /// \name type tensor and related operations
    //@{
  public:
    /// 3x3 matrix, not necessarily symmetric
    typedef real tensor[Ndim][Ndim];
    /// trace of tensor                                                         
    static double tr(double t[Ndim][Ndim]) { return meta::__tr<Ndim-1>::a(t); }
    /// trace of tensor                                                         
    static float  tr(float  t[Ndim][Ndim]) { return meta::__tr<Ndim-1>::a(t); }
    /// tensor computation as outer product or two vectors                      
    static void AddTensor(double t[Ndim][Ndim], vect_d const&x, vect const&y) {
      meta::__addt<Ndim-1>::a(t,
			      static_cast<const double*> (x),
			      static_cast<const real  *> (y));
    }
    /// tensor computation as outer product or two vectors                      
    static void AddTensor(float t[Ndim][Ndim], vect_f const&x, vect const&y) {
      meta::__addt<Ndim-1>::a(t,
			      static_cast<const float*> (x),
			      static_cast<const real *> (y));
    }
    /// asymmetric angular momentum tensor (funny idea, not mine anyway)        
    static real as_angmom     (vect const&L,
			       int  const&i,
			       int  const&j) {
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
    }
    //@}
    //--------------------------------------------------------------------------
    /// construction                                                            
    ForceAndDiagnose(snapshot          *s,
		     const acceleration*a) :
      SNAPSHOT(s), ACCEXTERN(a) {}
    //--------------------------------------------------------------------------
    /// \name abstract methods                                                  
    //@{                                                                        
    //--------------------------------------------------------------------------
    ///  compute forces (and other time derivatives if applicable)
    ///
    ///  IF \e dia is true, prepare for \c diagnose() to be called later.
    ///  \e tau may be needed for the time integration of auxiliary quantities,
    ///  such as SPH smoothing lengths.
    ///
    /// \param all (input) compute forces for all bodies or only the active?
    /// \param dia (input) diagnostics will be done (prepare for it)
    /// \param tau (input) size of elementary time step
    virtual void setforces  (bool   all,
			     bool   dia,
			     double tau) const = 0;
    //--------------------------------------------------------------------------
    /// perform a diagnose; only called after \c setforces(all,true,tau);
    virtual void     diagnose   ()   const = 0;
    //--------------------------------------------------------------------------
    /// which quantities are required in force computation and diagnosis?
    virtual fieldset requires   ()   const = 0;
    //--------------------------------------------------------------------------
    /// which SPH quantities are required in force computation and diagnosis?
    virtual fieldset requiresSPH()   const = 0;
    //--------------------------------------------------------------------------
    /// which quantities are computed by setforces()
    virtual fieldset computes   ()   const = 0;
    //--------------------------------------------------------------------------
    /// which SPH quantities are computed by setforces()
    virtual fieldset computesSPH()   const = 0;
    //--------------------------------------------------------------------------
    /// is a particular quantity computed?
    bool computes(fieldbit b) const {
      return computes().contain(fieldset(b));
    }
    //--------------------------------------------------------------------------
    /// is a particular quantity required?
    bool requires(fieldbit b) const {
      return requires().contain(fieldset(b));
    }
    //--------------------------------------------------------------------------
    /// is a particular SPH quantity computed?
    bool computesSPH(fieldbit b) const {
      return computesSPH().contain(fieldset(b));
    }
    //--------------------------------------------------------------------------
    /// is a particular SPH quantity required?
    bool requiresSPH(fieldbit b) const {
      return requiresSPH().contain(fieldset(b));
    }
    //--------------------------------------------------------------------------
    /// write diagnostic statistics
    virtual void dia_stats_body(std::ostream&) const=0;
    //--------------------------------------------------------------------------
    /// write CPU statistics
    virtual void cpu_stats_body(std::ostream&) const=0;
    //--------------------------------------------------------------------------
    /// write diagnostic statistics header
    virtual void dia_stats_head(std::ostream&) const=0;
    //--------------------------------------------------------------------------
    /// write CPU statistics header
    virtual void cpu_stats_head(std::ostream&) const=0;
    //--------------------------------------------------------------------------
    /// write diagnostic statistics line
    virtual void dia_stats_line(std::ostream&) const=0;
    //--------------------------------------------------------------------------
    /// write CPU statistics line
    virtual void cpu_stats_line(std::ostream&) const=0;
    //@}
    //--------------------------------------------------------------------------
    /// \name public data access
    //@{
    /// pointer to the N-body data
    snapshot*          const&snap_shot() const { return SNAPSHOT; }
    /// pointer to an external acceleration field (may be NULL)
    const acceleration*const&acc_ext  () const { return ACCEXTERN; }
    //@}
    //--------------------------------------------------------------------------
  };// class falcON::ForceAndDiagnose
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::Integrator                                                   
  //                                                                            
  /// abstract base class for N-body integrator                                 
  /// provides support for implementations (derived classes)                    
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class Integrator {
    Integrator(const Integrator&);
  private:
    const fieldset
      predALL,        ///< properties to be predicted for all bodies
      kickALL,        ///< properties to be predicted for SPH bodies only
      rembALL,        ///< properties to be remembered for all bodies
      predSPH,        ///< properties to be remembered for SPH bodies only
      kickSPH,        ///< properties to be kicked for all bodies
      rembSPH;        ///< properties to be kicked for SPH bodies only
    fieldset
      requALL,        ///< properties required for all bodies
      requSPH;        ///< properties required for SPH bodies only
  protected:
    /// pointer to the force solver
    const ForceAndDiagnose *const SOLVER;
  private:
    mutable clock_t    C_OLD;           
    mutable double     CPU_STEP,                   // total time for a longstep 
                       CPU_TOTAL;                  // total time so far         
    //--------------------------------------------------------------------------
  public:
    /// \name record CPU timings
    //@{
    /// add elapsed CPU time to CPU time record \c CPU
    ///
    /// \param c0  input: last CPU clock reading; output: current CPU clock
    /// \param CPU (output) time since \c c0 is added to CPU in seconds
    static void record_cpu(clock_t& c0, double& CPU)
    {
      register clock_t c1 = clock();
      CPU += (c1-c0)/real(CLOCKS_PER_SEC);
      c0   = c1;
    }
    //--------------------------------------------------------------------------
    /// print a CPU time
    static void print_cpu(double const&x, std::ostream&to)
    {
      if(x < 100)
	to<<std::setw(2)<<std::setfill(' ')<<int(x)<<'.'
	  <<std::setw(2)<<std::setfill('0')<<int(100*(x-int(x)));
      else if(x<1000)
	to<<std::setw(3)<<std::setfill(' ')<<int(x)<<'.'
	  <<std::setw(1)<<std::setfill('0')<<int(10*(x-int(x)));
      else
	to<<std::setw(5)<<std::setfill(' ')<<int(x+0.5);
    }
    //--------------------------------------------------------------------------
    /// print a CPU time in hhh:mm:ss.cc format
    static void print_cpu_hms(double t, std::ostream&to)
    {
      int    h,m,s,c;
      h = int(t/3600); t-= 3600*h;
      m = int(t/60);   t-= 60*m;
      s = int(t);      t-= s;
      c = int(100*t);
      to<<std::setw(3)<<std::setfill(' ')<<h<<':'
	<<std::setw(2)<<std::setfill('0')<<m<<':'
	<<std::setw(2)<<s<<'.'<<std::setw(2)<<c
	<<std::setfill(' ');
    }
  protected:
    //--------------------------------------------------------------------------
    /// record a CPU timing and cumulative CPU time after a full step
    void add_to_cpu_step() const {
      register clock_t c1 = clock();
      CPU_STEP    += (c1-C_OLD)/real(CLOCKS_PER_SEC);
      CPU_TOTAL   += (c1-C_OLD)/real(CLOCKS_PER_SEC);
      C_OLD        = c1;
    }
    //--------------------------------------------------------------------------
    /// reset CPU timing records
    void reset_CPU() const {
      CPU_STEP = 0.;
    }
    //--------------------------------------------------------------------------
    /// print CPU statistics, implements abstract method of base class
    void cpu_stats_body(std::ostream&) const;
    //--------------------------------------------------------------------------
    /// print CPU statistics header, implements abstract method of base class
    void cpu_stats_head(std::ostream&to) const {
      SOLVER->cpu_stats_head(to);
      to<< " step  accumulated";
    }
    //--------------------------------------------------------------------------
    /// print CPU statistics line, implements abstract method of base class
    void cpu_stats_line(std::ostream&to) const {
      SOLVER->cpu_stats_line(to);
      to<< "------------------";
    }
    //@}
    //--------------------------------------------------------------------------
    /// protected access to N-body data
    snapshot* const&snap_shot() const { return SOLVER->snap_shot(); }
    //--------------------------------------------------------------------------
    /// compute forces & SPH time derivatives, records CPU_GRAV, AEX, SPH
    void set_time_derivs(bool all,
			 bool dia,
			 double t) const {
      SOLVER->setforces(all,dia,t);
    }
    //--------------------------------------------------------------------------
    /// finish the diagnostics;
    /// previous call to \c set_time_derivs() had \c dia true
    void finish_diagnose() const { SOLVER->diagnose(); }
    //--------------------------------------------------------------------------
    /// specify which quantities will be predicted by this scheme
    fieldset const&predicted   () const { return predALL; }
    /// specify which SPH quantities will be predicted by this scheme
    fieldset const&predictedSPH() const { return predSPH; }
    /// specify which quantities will be required by this scheme
    fieldset const&required    () const { return requALL; }
    /// specify which SPH quantities will be required by this scheme
    fieldset const&requiredSPH () const { return requSPH; }
    //--------------------------------------------------------------------------
    /// drift by \c dt body properties in \a predicted() and \a predictedSPH()  
    ///                                                                         
    /// for all bodies: properties \a predicted() are moved by \e dt          \n
    /// for SPH bodies: properties \a predictedSPH() are moved by \e dt       \n
    /// currently we support the following drifts:                            \n
    /// <b> x += v * dt  </b>  position; for all bodies                       \n
    /// <b> w += a * dt  </b>  predicted velocities; for all bodies           \n
    /// <b> Y += I * dt  </b>  predicted SPH internal energy                    
    /// \param dt  (input) time step to drift                                   
    /// \param all (input) drift all or active bodies only?                     
    void drift   (double dt, bool all = true) const;
    //--------------------------------------------------------------------------
    /// kick by \c dt body properties to in \a kickALL and \a kickSPH           
    ///                                                                         
    /// currently we support the following kicks:                             \n
    /// <b> v += a * dt  </b>  velocity; for all bodies                       \n
    /// <b> U += I * dt  </b>  internal gas energy; for SPH bodies only         
    /// \param dt  (input) time step to kick                                    
    /// \param all (input) kick all or active bodies only?                      
    void kick    (double dt, bool all = true) const;
    //--------------------------------------------------------------------------
    /// similar to \a kick(), except that bodies are kicked by their individual 
    /// time step, given by \e dt[level(B)]                                     
    /// \param dt  (input) table with time step per level                       
    /// \param all (input) kick all or active bodies only?                      
    void kick_i  (const double*dt, bool all = true) const;
    //--------------------------------------------------------------------------
    /// remember properties to be predicted: in \a rembALL and \a rembSPH       
    ///                                                                         
    /// currently we support the following remembrances:                      \n
    /// <b> w = v  </b>  remember velocity                                    \n
    /// <b> Y = U  </b>  remember SPH internal energy                           
    /// \param all (input) remember for all or active bodies only?              
    void remember(bool all = true) const;
    //--------------------------------------------------------------------------
    /// construction: set properties for \a drift(), \a kick(), \a remember()   
    ///                                                                         
    /// we will double check that the input is sensible                         
    /// \param solver  pointer to solver for time derivatives                   
    /// \param p_all   -> \a predALL                                            
    /// \param k_all   -> \a kickALL                                            
    /// \param r_all   -> \a rembALL                                            
    /// \param p_sph   -> \a predSPH                                            
    /// \param k_sph   -> \a kickSPH                                            
    /// \param r_sph   -> \a rembSPH                                            
    Integrator(const ForceAndDiagnose* solver,
	       fieldset p_all,
	       fieldset k_all,
	       fieldset r_all,
	       fieldset p_sph,
	       fieldset k_sph,
	       fieldset r_sph) falcON_THROWING;
  public:
    //--------------------------------------------------------------------------
    /// \name virtual and pure virtual methods
    //@{
    //--------------------------------------------------------------------------
    /// perform a full time step
    virtual void fullstep() const = 0;
    //--------------------------------------------------------------------------
    /// print full statistic
    virtual void stats_body(std::ostream&to) const {
      SOLVER -> dia_stats_body(to);
      cpu_stats_body(to);
      to<<std::endl;
    }
    //--------------------------------------------------------------------------
    /// print statistic header
    virtual void stats_head(std::ostream&to) const {
      SOLVER -> dia_stats_head(to);
      cpu_stats_head(to);
      to<<std::endl;
    }
    //--------------------------------------------------------------------------
    /// print statistic line
    virtual void stats_line(std::ostream&to) const {
      SOLVER -> dia_stats_line(to);
      cpu_stats_line(to);
      to<<std::endl;
    }
    //@}
    //--------------------------------------------------------------------------
    // non-virtual public methods:                                              
    //--------------------------------------------------------------------------
    /// print some details 
    void describe(std::ostream&) const;
#ifdef falcON_NEMO
    //--------------------------------------------------------------------------
    /// write snapshot to NEMO output
    ///
    /// \param o NEMO output stream
    /// \param f data to be written out; default: mass, position, velocity
    void write(nemo_out const&o,
	       fieldset       f = fieldset::basic) const;
    //--------------------------------------------------------------------------
#endif // falcON_NEMO
  };// class falcON::Integrator
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::LeapFrogCode                                                 
  //                                                                            
  /// implements kick-drift-kick leap-frog                                      
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class LeapFrogCode : public Integrator, public bodies::TimeSteps {
    void account_new() const;
  public:
    //--------------------------------------------------------------------------
    /// construction
    ///
    /// \param k       tau = 0.5^k
    /// \param solver  pointer to solver for time derivatives
    /// \param p_all   -> \a predALL
    /// \param k_all   -> \a kickALL
    /// \param r_all   -> \a rembALL
    /// \param p_sph   -> \a predSPH
    /// \param k_sph   -> \a kickSPH
    /// \param r_sph   -> \a rembSPH
    LeapFrogCode(int k,
		 const ForceAndDiagnose* solver,
		 fieldset p_all,
		 fieldset k_all,
		 fieldset r_all,
		 fieldset p_sph,
		 fieldset k_sph,
		 fieldset r_sph) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// perform a full kick-drift-kick step
    void fullstep() const;
    //--------------------------------------------------------------------------
  };// class falcON::LeapFrogCode
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::BlockStepCode                                                
  //                                                                            
  //  implements a hierarchical blockstep scheme with kick-drift-kick leap-frog 
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class BlockStepCode : public Integrator, public bodies::TimeSteps {
  public:
    //--------------------------------------------------------------------------
    /// abstract sub-type: interface for assigning/adjusting time-step levels   
    struct StepLevels {
      typedef double *const c_pdouble;
      /// pure virtual: assign time-step level to body
      /// \param B  body to assign level to
      /// \param N  table for # bodies / step
      /// \param H  highest table index
      virtual void assign_level(body&B, unsigned *N, int H) const =0;
      /// pure virtual: adjust time-step level of body
      /// \param B  body whose level to adjust
      /// \param N  table for # bodies / step
      /// \param L  lowest allowed level
      /// \param H  highest allowed level
      virtual void adjust_level(body&B, unsigned *N, int L, int H) const =0;
      /// virtual function: shall adjust_level() always be called,
      /// even if there is only one level?
      virtual bool always_adjust () const { return false; }
    };
    //--------------------------------------------------------------------------
    // data members                                                             
    //--------------------------------------------------------------------------
  private:
    unsigned              *N;                      // table:  #                 
    int                    W;                      // width in stats fields     
    const StepLevels*const ST;                     // for adjusing levels       
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
  private:
    void assign_levels() const;
    void adjust_levels(int, bool) const;
    void update_Nlev(const bodies*);
    void account_del() const;
    void account_new() const;
    void elementary_step(int) const;
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    const unsigned&No_in_level   (int l)  const { return N[l]; }
    //--------------------------------------------------------------------------
    // to satisfy class Integrator::fullstep()                                  
    void fullstep() const;
    //--------------------------------------------------------------------------
    // statistics output                                                        
    // we have to superseed Integrator::stats_etc.., in order to add step stats 
    void stats_head(std::ostream&) const;
    //--------------------------------------------------------------------------
    void stats_line(std::ostream&to) const {
      SOLVER -> dia_stats_line(to);
      if(highest_level())
	for(int l=0; l!=Nsteps(); ++l) 
	  for(int i=0; i<=W; ++i) to<<'-';
      cpu_stats_line(to);
      to<<std::endl;
    }
    //--------------------------------------------------------------------------
    void stats_body(std::ostream&to) const {
      SOLVER -> dia_stats_body(to);
      if(highest_level())
	for(int l=0; l!=Nsteps(); ++l)
	  to<<std::setw(W)<<N[l]<<' ';
      cpu_stats_body(to);
      to<<std::endl;
    }
    //--------------------------------------------------------------------------
    // construction & destruction                                               
    //--------------------------------------------------------------------------
    BlockStepCode (int,                            // I: kmax, tau_max = 2^-kmax
		   unsigned,                       // I: #steps                 
		   const ForceAndDiagnose*,
		   const StepLevels      *,
		   fieldset, fieldset, fieldset, fieldset, fieldset, fieldset,
		   int) falcON_THROWING;
    //--------------------------------------------------------------------------
    ~BlockStepCode () { 
      if(N) { falcON_DEL_A(N); N=0; }
    }
  };// class falcON::BlockStepCode   
#ifdef falcON_NEMO
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::NBodyCode                                                  //
  //                                                                          //
  // 10-11-2004: added parameter times to constructor                         //
  // 06-10-2005: moved TINI (initial time) to snapshot                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class nemo_out;
  class NBodyCode {
  protected:
    std::string        FILE;
    snapshot           SHOT;
    const Integrator  *CODE;
    fieldset           READ;
    //--------------------------------------------------------------------------
    // construction & related                                                   
    //                                                                          
    // NOTE                                                                     
    //      for technical reasons, we must follow the following order:          
    //      1. load initial snapshot (& allocate required memory)               
    //      2. initialize the forcesolver (it needs the snapshot)               
    //      3. initialize the integrator                                        
    //      Since the force solver is dealt with in a derived class, we CANNOT  
    //      follow this procedure in the constructor below. Therefore, we merely
    //      load the initial snapshot and postpone the initialisation of the    
    //      integrator, which shall be done via the routine init() below.       
    //--------------------------------------------------------------------------
    // read initial snapshot                                                    
    NBodyCode(const char*,                         // I: input file             
	      bool       ,                         // I: resume old (if nemo)   
	      fieldset   ,                         // I: read after mxv         
	      const char* =0,                      // I: read 1st snapshot in   
	                                           //    1st file matching      
	      fieldset    =fieldset::empty)        // I: try to read these data 
      falcON_THROWING;
    //--------------------------------------------------------------------------
    // initialize Integrator, see NOTE above                                    
    void init(const ForceAndDiagnose         *,    // I: e.g. ForceALCON        
	      int                             ,    // I: kmax: t_max=2^(-kmax)  
	      int                             ,    // I: # time-step levels     
	      const BlockStepCode::StepLevels*,
	      fieldset, fieldset,
	      fieldset=fieldset::empty,
	      fieldset=fieldset::empty,
	      fieldset=fieldset::empty,
	      fieldset=fieldset::empty)
      falcON_THROWING;
    //--------------------------------------------------------------------------
    // destruction                                                              
    //--------------------------------------------------------------------------
    ~NBodyCode() {                                 // destructor                
      falcON_DEL_O(CODE);
    }
    //--------------------------------------------------------------------------
    // nemo outputs                                                             
    //--------------------------------------------------------------------------
  public:
    void  describe(                                // describe simulation       
		   std::ostream&o) const {         // I: output stream          
      CODE->describe(o);
    }
    //--------------------------------------------------------------------------
    void  write   (nemo_out&o,                     // I: nemo output stream     
		   fieldset w) const {             //[I: what to write out]     
      CODE->write(o,w);
    }
    //--------------------------------------------------------------------------
    // time integration                                                         
    //--------------------------------------------------------------------------
    void  full_step    () {
      CODE->fullstep();
    }
    //--------------------------------------------------------------------------
    // statistic outputs                                                        
    //--------------------------------------------------------------------------
    void  stats     (std::ostream&o) const { 
      o <<' ';
      CODE->stats_body(o);
    }
    //--------------------------------------------------------------------------
    void  stats_head(std::ostream&o) const {
      o<<'#'; CODE->stats_head(o);
      o<<'#'; CODE->stats_line(o);
    }
    //--------------------------------------------------------------------------
    // data access                                                              
    //--------------------------------------------------------------------------
    double      const&initial_time   () const { return SHOT.initial_time(); }
    const bodies     *my_bodies      () const { return&SHOT; }
    const snapshot   *my_snapshot    () const { return&SHOT; }
    std::string const&input_file     () const { return FILE; }
    double      const&time           () const { return SHOT.time(); }
    //--------------------------------------------------------------------------
  };// class falcON::NBodyCode
#endif // #ifdef falcON_NEMO
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::ForceDiagGrav                                              //
  //                                                                          //
  // derived from ForceAndDiagnose                                            //
  // provides basis diagnostics for gravity only codes                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class ForceDiagGrav : public ForceAndDiagnose {
    //--------------------------------------------------------------------------
    // diagnostics data                                                         
    //--------------------------------------------------------------------------
  protected:
    const   bool         SELF_GRAV;                // compute V_in?             
    mutable double       TIME;                     // time of diagnose          
    mutable double       M,T,Vin,Vex,W,TW;         // mass, kin & pot E         
    mutable vect_d       L;                        // total angular momentum    
    mutable double       DVDT;                     // dV/dt = - dT/dt           
    mutable tensor       KT,WT;                    // kin & pot energy          
    mutable vect_d       CMX,CMV;                  // center of mass pos & vel  
    //--------------------------------------------------------------------------
    void diagnose_grav() const;
    void diagnose_vels() const falcON_THROWING;
    void diagnose_full() const;
    //--------------------------------------------------------------------------
    ForceDiagGrav(snapshot          *s,
		  const acceleration*a,
		  bool               g) : ForceAndDiagnose(s,a), SELF_GRAV(g) {}
    //--------------------------------------------------------------------------
  public:
    virtual fieldset requires() const {
      return fieldset(fieldset::m |
		      fieldset::x |
		      fieldset::v |
		      fieldset::a |
		      fieldset::p |
		      (acc_ext()? fieldset::q : fieldset::empty) );
    }
    //--------------------------------------------------------------------------
    virtual double Ekin() const { return T; }
    virtual double Epot() const { return Vin + Vex; }
    virtual double Etot() const { return Ekin() + Epot(); }
    //--------------------------------------------------------------------------
    double const  &Vrat() const { return TW; }
    double         Wvir() const { return W; }
    double const  &dVdt() const { return DVDT; }
    vect_d const  &Ltot() const { return L; }
    vect_d const  &Xave() const { return CMX; }
    vect_d const  &Vave() const { return CMV; }
    //--------------------------------------------------------------------------
    virtual void dia_stats_body (std::ostream&) const;
    virtual void dia_stats_head (std::ostream&) const;
    virtual void dia_stats_line (std::ostream&) const;
    //--------------------------------------------------------------------------
    bool const&self_grav() const { return SELF_GRAV; }

  };// class falcON::ForceDiagGrav
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::GravStepper                                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class GravStepper {
    //--------------------------------------------------------------------------
    // enum used in switch() statements                                         
    //--------------------------------------------------------------------------
    enum {
      use_a = 1,
      use_p = 2,
      use_g = 4,
      use_e = 8
    };
    //--------------------------------------------------------------------------
    typedef BlockStepCode::StepLevels::c_pdouble c_pdouble;
    //--------------------------------------------------------------------------
    // data members                                                             
    //--------------------------------------------------------------------------
    bool    UPX;                                   // use pot+pex for potential 
    bool    UGE;                                   // use global soften'g length
    int     SCH;                                   // stepping scheme           
    int     SINKL;                                 // minimum level for sinks   
    double  EPS;                                   // global softening length   
    double  FAQ,FPQ,FGQ,FEQ;                       // factors^2 for stepping    
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
    real fpot(body const&B) const {
      return UPX? pex(B)+pot(B) : pot(B);
    }
    //--------------------------------------------------------------------------
    real soft(body const&B) const {
      return UGE? EPS : eps(B);
    }
    //--------------------------------------------------------------------------
    int minlevel(body const&B) const {
      return is_sink(B)? SINKL : 0;
    }
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
  protected:
    double tq_grav(body const&B) const {
      if(SCH == 0)     return zero;
      if(SCH == use_p) return FPQ/square(fpot(B));
      else {
	double ia = 1./norm(acc(B));
	double tq = 1.e7;
	if(SCH & use_a) update_min(tq, FAQ*ia);
	if(SCH & use_p) update_min(tq, FPQ/square(fpot(B)));
	if(SCH & use_g) update_min(tq, FGQ*abs(fpot(B))*ia);
	if(SCH & use_e) update_min(tq, FEQ*soft(B)*std::sqrt(ia));
	return tq;
      }
    }
    //--------------------------------------------------------------------------
    void assign_level(body        &Bi,             // I: body                   
		      unsigned    *N,              // I: table: # / step        
		      int          H) const        // I: highest table index    

    {
      double tq=twice(tq_grav(Bi));
      for(Bi.level()=minlevel(Bi); tauq(Bi)>tq && level(Bi)<H; ++(Bi.level()));
      N[level(Bi)]++;
    }
    //--------------------------------------------------------------------------
    void adjust_level(body        &Bi,             // I: bodies                 
		      unsigned    *N,              // I: table: # / step        
		      int          L,              // I: lowest  allowed level  
		      int          H) const        // I: highest allowed level  
    {
      const double root_half=0.7071067811865475244;
      double tq=tq_grav(Bi)*root_half;             // tau^2 * sqrt(1/2)         
      L = max(L,minlevel(Bi));                     // set lowest level for sinks
      if       (tauq(Bi) < tq) {                   // IF too short              
	if(level(Bi) > L) {                        //   IF(above lowest active) 
	  N[level(Bi)] --;                         //     leave old level       
	  Bi.level()   --;                         //     set new level         
	  N[level(Bi)] ++;                         //     join new level        
	}                                          //   ENDIF                   
      } else if(tauq(Bi) > tq+tq) {                // ELIF too long             
	if(level(Bi) < H) {                        //   IF(below highest level) 
	  N[level(Bi)] --;                         //     leave old level       
	  Bi.level()   ++;                         //     set new level         
	  N[level(Bi)] ++;                         //     join new level        
	}                                          //   ENDIF                   
      }                                            // ENDIF                     
    }
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    GravStepper(int  sl,
		real fa,
		real fp,
		real fg,
		real fe,
		bool up,
		real ep = zero,
		bool ue = 0) : UPX(up), UGE(ue), EPS(ep), SCH(0), SINKL(sl)
    {
      FAQ = fa*fa; if(FAQ) SCH |= use_a;
      FPQ = fp*fp; if(FPQ) SCH |= use_p;
      FGQ = fg*fg; if(FGQ) SCH |= use_g;
      FEQ = fe*fe; if(FEQ) SCH |= use_e;
    }
    //--------------------------------------------------------------------------
    int const&scheme() const { return SCH; }
  };// class falcON::GravStepper
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::GravSteps                                                  //
  //                                                                          //
  // non-abstract, derived from abstract class BlockStepCode::StepLevels      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class GravSteps : 
    public BlockStepCode::StepLevels,
    public GravStepper 
  {
    //--------------------------------------------------------------------------
    typedef BlockStepCode::StepLevels::c_pdouble c_pdouble;
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
  protected:
    void assign_level(body        &Bi,             // I: body                   
		      unsigned    *N,              // I: table: # / step        
		      int          H) const        // I: highest table index    
    { GravStepper::assign_level(Bi,N,H); }
    //--------------------------------------------------------------------------
    void adjust_level(body        &Bi,             // I: bodies                 
		      unsigned    *N,              // I: table: # / step        
		      int          L,              // I: lowest  allowed level  
		      int          H) const        // I: highest allowed level  
    { GravStepper::adjust_level(Bi,N,L,H); }
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    GravSteps(int  sl,
	      real fa,
	      real fp,
	      real fg,
	      real fe,
	      bool up,
	      real ep = zero,
	      bool ue = 0) : GravStepper(sl,fa,fp,fg,fe,up,ep,ue) {}
    //--------------------------------------------------------------------------
  };// class falcON::GravSteps
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::ForceALCON                                                 //
  //                                                                          //
  // standard force solver based on forces, providing gravity only            //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class ForceALCON : public ForceDiagGrav {
    //--------------------------------------------------------------------------
  protected:
#ifdef falcON_INDI
    soft_type         SOFTENING;
#ifdef falcON_ADAP
    real              NSOFT;                       // #/sphere for eps_i setting
    unsigned          NREF;                        // #/cell for eps_i setting  
    real              EMIN;                        // lower limit for eps_i     
    real              EFAC;                        // max change factor for epsi
#endif
#endif
    const   vect* const ROOTCENTRE;
    const   int         NCRIT, REUSE;
    mutable forces      FALCON;
    mutable int         REUSED;
    mutable double      CPU_TREE, CPU_GRAV, CPU_AEX;
    //--------------------------------------------------------------------------
    void set_tree_and_forces(bool,                 // I: all or only active?    
			     bool) const;          // I: build tree anyway?     
    //--------------------------------------------------------------------------
    // - sets Nsoft, Nref, eps, kernel, but NOT soft_type!                      
  public:
    void reset_softening(                          // resets softening params   
			 kern_type,                // I: softening kernel       
			 real                      // I: eps                    
#ifdef falcON_ADAP
			,real     ,                // I: Nsoft                  
			 unsigned ,                // I: Nref                   
			 real     ,                // I: eps_min                
			 real                      // I: eps_fac                
#endif
			 );
    //--------------------------------------------------------------------------
    // construction                                                             
    ForceALCON(snapshot          *,                // I: snapshot: time & bodies
	       real               ,                // I: global softening length
	       real               ,                // I: tolerance parameter    
	       int                ,                // I: N_crit                 
	       const vect        *,                // I: pre-set root centre    
	       kern_type          ,                // I: softening kernel       
	       real               ,                // I: Newton's G             
	       real               ,                // I: theta_sink/theta       
	       int                ,                // I: # reused of tree       
	       const acceleration*,                // I: external acceleration  
	       const int          [4]              // I: direct sum: gravity    
#ifdef falcON_INDI
	       ,soft_type                          // I: softening type         
#ifdef falcON_ADAP
	       ,real                               // I: N_soft                 
	       ,unsigned                           // I: N_ref                  
	       ,real                               // I: eps_min                
	       ,real                               // I: eps_fac                
#endif
#endif
#ifdef falcON_SPH
	       ,const int         [3]              //[I: direct sum: SPH]       
	       = Default::SPHdirect
#endif
	       ) falcON_THROWING;
    //--------------------------------------------------------------------------
    // satisfy all the abstract functions of ForceDiagGrav                      
    //                                                                          
    // functions are virtual to allow a class derived from ForceALCON to serve  
    // in NBodyCode (otherwise the functions below are taken).                  
    //--------------------------------------------------------------------------
    virtual void cpu_stats_head(std::ostream&to) const {
      if(SELF_GRAV) to << "l2R  D  tree  grav ";
      if(acc_ext()) to << " pext ";
    }
    //--------------------------------------------------------------------------
    virtual void cpu_stats_line(std::ostream&to) const {
      if(SELF_GRAV) to << "-------------------";
      if(acc_ext()) to << "------";
    }
    //--------------------------------------------------------------------------
    virtual void cpu_stats_body(std::ostream&) const;
    //--------------------------------------------------------------------------
    virtual fieldset computes() const { 
      return fieldset(fieldset::p |
		      fieldset::a |
		      (acc_ext()? fieldset::q : fieldset::empty) );
    }
    //--------------------------------------------------------------------------
    // note: we don't require eps_i, but use them as part of force computation  
    virtual fieldset requires() const {
      return
	( fieldset(fieldset::m | fieldset::x) | ForceDiagGrav::requires() )
	& ~computes();
    }
    //--------------------------------------------------------------------------
    virtual fieldset requiresSPH () const { return fieldset::empty; } 
    virtual fieldset computesSPH () const { return fieldset::empty; } 
    //--------------------------------------------------------------------------
    virtual void setforces (bool all, bool, double) const {
      set_tree_and_forces(all,false);
      debug_info(5,"ForceALCON::setforces(): done\n");
    }
    virtual void diagnose  () const { ForceDiagGrav::diagnose_full(); }
    //--------------------------------------------------------------------------
  };// class falcON::ForceALCON
#ifdef falcON_NEMO
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::FalcONCode                                                 //
  //                                                                          //
  // completely inline                                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class FalcONCode :
    public  NBodyCode,
    private ForceALCON {
  private:
    GravSteps GS;
    //--------------------------------------------------------------------------
    /// constructor
    ///
    /// \param file   data input file for N-body data
    /// \param resume resume old simulation?
    /// \param kmax   tau_max=(1/2)^kmax
    /// \param Nlev   # time-step levels
    /// \param ksink  minimum time step for sinks: (1/2)^ksink
    /// \param fa     time step factor f_acc
    /// \param fp     time step factor f_phi
    /// \param fc     time step factor f_c
    /// \param fe     time step factor f_eps
    /// \param Ncrit  N_crit for tree building
    /// \param hgrow  build tree only every 2^hgrow shortest steps
    /// \param croot  if non-null: pter to pre-set root centre position
    /// \param eps    if > 0: global softening length; if < 0: use eps_i
    /// \param kernel softening kernel
    /// \param aex    pter to external acceleration field
    /// \param theta  opening angle
    /// \param Grav   Newton's constant of gravity
    /// \param sfac   theta_sink/theta
#ifdef falcON_INDI
#ifdef falcON_ADAP
    /// \param Nsoft  # of bodies in softening shpere
    /// \param Nref   # of bodies in smallest cell to estimate eps_i
    /// \param emin   minimum eps_i
#endif
    /// \param soft   softening type
#endif
    /// \param time   time to read input for (take first if null)
    /// \param read   data to read in addition to required for integration
    /// \param dir    number controlling direct summation
  public:
    FalcONCode(// data input                                                    
	       const char *       file,
	       bool               resume,
	       // time-integration related                                      
	       int                kmax,
	       int                Nlev,
	       int                ksink,
	       real               fa,
	       real               fp,
	       real               fc,
	       real               fe,
	       // tree related                                                  
	       int                Ncrit,
	       int                hgrow,
	       const vect        *croot,
	       // gravity related                                               
	       real               eps,
	       kern_type          kernel,
	       const acceleration*aex,
	       real               theta,
	       real               Grav,
	       real               sfac,
	       // default arguments                                             
#ifdef falcON_INDI
#ifdef falcON_ADAP
	       real               Nsoft = zero,
	       unsigned           Nref = 32,
	       real               emin = zero,
#endif
	       soft_type          soft = global_fixed, 
#endif
	       const char        *time = 0,
	       fieldset           read = fieldset::empty,
	       const int          dir[4] = Default::direct) falcON_THROWING :
      NBodyCode ( file, resume, read | fieldset(
#ifdef falcON_INDI
		  soft != global_fixed ? fieldset::e :
#endif
		  fieldset::empty), time ),
      ForceALCON( &SHOT, eps, theta, Ncrit, croot, kernel, Grav, sfac,
		  (1<<hgrow)-1, aex, dir
#ifdef falcON_INDI
		  ,soft
#ifdef falcON_ADAP
		  ,Nsoft, Nref, emin, 1.1
#endif
#endif
		  ),
      GS         ( ksink-kmax, fa,fp,fc,fe, aex!=0, abs(eps), eps>=zero)
    {
#ifdef falcON_ADAP
      if(soft == individual_adaptive && hgrow)
	falcON_THROW("inidividual adaptive softening lengths "
		     "must not be used with hgrow != 0");
#endif
      if(ksink > kmax+Nlev-1)
	falcON_THROW("FalcONCode: ksink=%d > kmin=%d\n",ksink,kmax+Nlev-1);
      if(ksink < kmax)
	falcON_THROW("FalcONCode: ksink=%d < kmax=%d\n",ksink,kmax);
      if(Nlev > 1 && GS.scheme() == 0)
	falcON_THROW("FalcONCode: time stepping factors=0\n");
      NBodyCode::init(this, kmax, Nlev, &GS, fieldset::x, fieldset::v);
    }
  };// class falcON::FalcONCode
  //////////////////////////////////////////////////////////////////////////////
#endif // #ifdef falcON_NEMO
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // included_falcON_nbody_h

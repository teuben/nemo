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
#ifndef falcON_included_io_h
#  include <public/io.h>
#endif
#ifndef falcON_included_nemopp_h
#  include <public/nemo++.h>
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
    /// \name data
    //{@
    snapshot*          const SNAPSHOT;   ///< particle snapshot
    const acceleration*const ACCEXTERN;  ///< external potential/acceleration
    //}@
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
    /// compute forces (and other time derivatives if applicable).
    /// IF \e dia is true, prepare for \c diagnose() to be called later.
    /// \e tau may be needed for the time integration of auxiliary quantities,
    /// such as SPH smoothing lengths.
    /// \param all (input) compute forces for all bodies or only the active?
    /// \param dia (input) diagnostics will be done (prepare for it)
    /// \param tau (input) size of elementary time step
    virtual void setforces  (bool   all,
			     bool   dia,
			     double tau) const = 0;
    /// perform a diagnose; only called after \c setforces(all,true,tau);
    virtual void     diagnose   ()   const = 0;
    /// which quantities are required in force computation and diagnosis?
    virtual fieldset requires   ()   const = 0;
    /// which SPH quantities are required in force computation and diagnosis?
    virtual fieldset requiresSPH()   const = 0;
    /// which quantities are computed by setforces()
    virtual fieldset computes   ()   const = 0;
    /// which SPH quantities are computed by setforces()
    virtual fieldset computesSPH()   const = 0;
    /// write diagnostic statistics
    virtual void dia_stats_body(output&) const=0;
    /// write CPU statistics
    virtual void cpu_stats_body(output&) const=0;
    /// write diagnostic statistics header
    virtual void dia_stats_head(output&) const=0;
    /// write CPU statistics header
    virtual void cpu_stats_head(output&) const=0;
    /// write diagnostic statistics line
    virtual void dia_stats_line(output&) const=0;
    /// write CPU statistics line
    virtual void cpu_stats_line(output&) const=0;
    //@}
    //--------------------------------------------------------------------------
    /// is a particular quantity computed?
    bool computes(fieldbit b) const {
      return computes().contain(fieldset(b));
    }
    /// is a particular quantity required?
    bool requires(fieldbit b) const {
      return requires().contain(fieldset(b));
    }
    /// is a particular SPH quantity computed?
    bool computesSPH(fieldbit b) const {
      return computesSPH().contain(fieldset(b));
    }
    /// is a particular SPH quantity required?
    bool requiresSPH(fieldbit b) const {
      return requiresSPH().contain(fieldset(b));
    }
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
    /// \param c0  input: last CPU clock reading; output: current CPU clock
    /// \param CPU (output) time since \c c0 is added to CPU in seconds
    static void record_cpu(clock_t& c0, double& CPU)
    {
      register clock_t c1 = clock();
      CPU += (c1-c0)/real(CLOCKS_PER_SEC);
      c0   = c1;
    }
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
    /// record a CPU timing and cumulative CPU time after a full step
    void add_to_cpu_step() const {
      register clock_t c1 = clock();
      CPU_STEP    += (c1-C_OLD)/real(CLOCKS_PER_SEC);
      CPU_TOTAL   += (c1-C_OLD)/real(CLOCKS_PER_SEC);
      C_OLD        = c1;
    }
    /// reset CPU timing records
    void reset_CPU() const {
      CPU_STEP = 0.;
    }
    /// print CPU statistics, implements abstract method of base class
    void cpu_stats_body(output&) const;
    /// print CPU statistics header, implements abstract method of base class
    void cpu_stats_head(output&to) const {
      SOLVER->cpu_stats_head(to);
      if(to) to<< " step  accumulated";
    }
    /// print CPU statistics line, implements abstract method of base class
    void cpu_stats_line(output&to) const {
      SOLVER->cpu_stats_line(to);
      if(to) to<< "------------------";
    }
    //@}
    //--------------------------------------------------------------------------
    /// protected access to N-body data
    snapshot* const&snap_shot() const { return SOLVER->snap_shot(); }
    /// compute forces & SPH time derivatives, records CPU_GRAV, AEX, SPH
    void set_time_derivs(bool all,
			 bool dia,
			 double t) const {
      SOLVER->setforces(all,dia,t);
    }
    /// finish the diagnostics;
    /// previous call to \c set_time_derivs() had \c dia true
    void finish_diagnose() const { SOLVER->diagnose(); }
    /// specify which quantities will be predicted by this scheme
    fieldset const&predicted   () const { return predALL; }
    /// specify which SPH quantities will be predicted by this scheme
    fieldset const&predictedSPH() const { return predSPH; }
    /// specify which quantities will be required by this scheme
    fieldset const&required    () const { return requALL; }
    /// specify which SPH quantities will be required by this scheme
    fieldset const&requiredSPH () const { return requSPH; }
    /// drift by \c dt body properties in \a predicted() and \a predictedSPH().
    /// for all bodies: properties \a predicted() are moved by \e dt          \n
    /// for SPH bodies: properties \a predictedSPH() are moved by \e dt       \n
    /// currently we support the following drifts:                            \n
    /// <b> x += v * dt  </b>  position; for all bodies                       \n
    /// <b> w += a * dt  </b>  predicted velocities; for all bodies           \n
    /// <b> Y += I * dt  </b>  predicted SPH internal energy                    
    /// \param dt  (input) time step to drift                                   
    /// \param all (input) drift all or active bodies only?                     
    void drift   (double dt, bool all = true) const;
    /// kick by \c dt body properties to in \a kickALL and \a kickSPH.
    /// currently we support the following kicks:                             \n
    /// <b> v += a * dt  </b>  velocity; for all bodies                       \n
    /// <b> U += I * dt  </b>  internal gas energy; for SPH bodies only         
    /// \param dt  (input) time step to kick                                    
    /// \param all (input) kick all or active bodies only?                      
    void kick    (double dt, bool all = true) const;
    /// similar to \a kick(), except that bodies are kicked by their individual 
    /// time step, given by \e dt[level(B)]                                     
    /// \param dt  (input) table with time step per level                       
    /// \param all (input) kick all or active bodies only?                      
    void kick_i  (const double*dt, bool all = true) const;
    /// remember properties to be predicted: in \a rembALL and \a rembSPH.
    /// currently we support the following remembrances:                      \n
    /// <b> w = v  </b>  remember velocity                                    \n
    /// <b> Y = U  </b>  remember SPH internal energy                           
    /// \param all (input) remember for all or active bodies only?              
    void remember(bool all = true) const;
    /// construction: set properties for \a drift(), \a kick(), \a remember().
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
    /// perform a full time step
    virtual void fullstep() const = 0;
    /// print full statistic
    virtual void stats_body(output&to) const {
      SOLVER -> dia_stats_body(to);
      cpu_stats_body(to);
      if(to) to<<std::endl;
    }
    /// print statistic header
    virtual void stats_head(output&to) const {
      SOLVER -> dia_stats_head(to);
      cpu_stats_head(to);
      if(to) to<<std::endl;
    }
    /// print statistic line
    virtual void stats_line(output&to) const {
      SOLVER -> dia_stats_line(to);
      cpu_stats_line(to);
      if(to) to<<std::endl;
    }
    //@}
    //--------------------------------------------------------------------------
    /// print some details 
    void describe(output&) const;
#ifdef falcON_NEMO
    /// write snapshot to NEMO output.
    /// \param[in] o NEMO output stream
    /// \param[in] f data to be written out; default: mass, position, velocity
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
    /// construction.
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
    /// perform a full kick-drift-kick step
    void fullstep() const;
    //--------------------------------------------------------------------------
  };// class falcON::LeapFrogCode
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::BlockStepCode                                                
  //                                                                            
  /// implements a hierarchical blockstep scheme with kick-drift-kick leap-frog 
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class BlockStepCode : public Integrator, public bodies::TimeSteps {
  public:
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
  private:
    //--------------------------------------------------------------------------
    /// \name data members
    //@{
    unsigned              *N;                    ///< table: # bodies per level
    int                    W;                    ///< width in stats fields     
    const StepLevels*const ST;                   ///< for adjusing levels       
    //@}
    //--------------------------------------------------------------------------
  private:
    void assign_levels() const;
    void adjust_levels(int, bool) const;
    void update_Nlev(const bodies*);
    void account_del() const;
    void account_new() const;
    void elementary_step(int) const;
  public:
    //--------------------------------------------------------------------------
    /// how many bodies are in level l?
    const unsigned&No_in_level   (int l)  const { return N[l]; }
    /// implements Integrator::fullstep()
    void fullstep() const;
    /// implements Integrator::stats_head
    void stats_head(output&) const;
    /// implements Integrator::stats_line
    void stats_line(output&to) const {
      SOLVER -> dia_stats_line(to);
      if(to && highest_level())
	for(int l=0; l!=Nsteps(); ++l) 
	  for(int i=0; i<=W; ++i) to<<'-';
      cpu_stats_line(to);
      if(to) to<<std::endl;
    }
    /// implements Integrator::stats_head
    void stats_body(output&to) const {
      SOLVER -> dia_stats_body(to);
      if(to && highest_level())
	for(int l=0; l!=Nsteps(); ++l)
	  to<<std::setw(W)<<N[l]<<' ';
      cpu_stats_body(to);
      if(to) to<<std::endl;
    }
    /// construction
    /// \param[in] kmax  \f$ \tau_{\mathrm{max}}=2^{-\mathrm{kmax}} \f$
    /// \param[in] nlev  number of time step levels
    /// \param[in] F     code for force evaluation etc
    /// \param[in] S     code for time-stepping levels and criteria
    /// \param[in] p_all \a Integrator::predALL
    /// \param[in] k_all \a Integrator::kickALL
    /// \param[in] r_all \a Integrator::rembALL
    /// \param[in] p_sph \a Integrator::predSPH
    /// \param[in] k_sph \a Integrator::kickSPH
    /// \param[in] r_sph \a Integrator::rembSPH
    /// \param[in] w     width in stats fields
    BlockStepCode (int kmax, unsigned nlev, const ForceAndDiagnose*F,
		   const StepLevels*S,
		   fieldset p_all, fieldset k_all, fieldset r_all,
		   fieldset p_sph, fieldset k_sph, fieldset r_sph,
		   int w) falcON_THROWING;
    /// destruction
    ~BlockStepCode () { 
      if(N) { falcON_DEL_A(N); N=0; }
    }
  };// class falcON::BlockStepCode   
} // namespace falcON {
//------------------------------------------------------------------------------
falcON_TRAITS(falcON::Integrator,"Integrator");
falcON_TRAITS(falcON::LeapFrogCode,"LeapFrogCode");
falcON_TRAITS(falcON::BlockStepCode,"BlockStepCode");
//------------------------------------------------------------------------------
namespace falcON {
#ifdef falcON_NEMO
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::NBodyCode                                                    
  //                                                                            
  /// puts together a (virtual) N-body code.                                    
  /// \version 10-11-2004: added parameter times to constructor                 
  /// \version 06-10-2005: moved TINI (initial time) to snapshot                
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class nemo_out;
  class NBodyCode {
  protected:
    std::string        FILE;   ///< input file
    snapshot          *SHOT;   ///< snapshot data
    ParallelSnapshot  *PSHT;   ///< if MPI job: body data including remote
    const Integrator  *CODE;   ///< time integrator
    fieldset           READ;   ///< which data have been read?
    /// construction
    /// \note For technical reasons, we must follow the following order.\n
    ///       1. load initial snapshot (& allocate required memory);\n
    ///       2. initialize the forcesolver (it needs the snapshot);\n
    ///       3. initialize the integrator\n
    ///       Since the force solver is dealt with in a derived class, we 
    ///       cannot follow this procedure in the constructor below. Therefore,
    ///       we merely load the initial snapshot and postpone the 
    ///       initialisation of the integrator, which shall be done via the
    ///       routine init() below.
    /// \param[in] file   name of data input file (in NEMO snapshot format)
    /// \param[in] resume resume old simulation: read last snapshot from file
    /// \param[in] read_more which data to read in addition to mxv
    /// \param[in] time (optional) if given, read snapshot matching time
    /// \param[in] read_try (optional) if given, try to read these data
    /// \version 29-10-2008 enabled MPI parallelism
    NBodyCode(const char* file, bool resume, fieldset read_more,
	      const char* time=0, fieldset read_try=fieldset::empty)
      falcON_THROWING;
    //--------------------------------------------------------------------------
    /// initialize Integrator, see note with constructor.
    /// \param[in] F force solver, e.g. ForceALCON
    /// \param[in] kmax  \f$ \tau_{\mathrm{max}}=2^{-\mathrm{kmax}} \f$
    /// \param[in] nlev  number of time step levels
    /// \param[in] S     time step levels (needed if nlev>1)
    /// \param[in] p_all \a Integrator::predALL
    /// \param[in] k_all \a Integrator::kickALL
    /// \param[in] r_all \a Integrator::rembALL
    /// \param[in] p_sph \a Integrator::predSPH
    /// \param[in] k_sph \a Integrator::kickSPH
    /// \param[in] r_sph \a Integrator::rembSPH
    void init(const ForceAndDiagnose*F, int kmax, int nlev,
	      const BlockStepCode::StepLevels*S,
	      fieldset p_all, fieldset k_all,
	      fieldset r_all=fieldset::empty,
	      fieldset p_sph=fieldset::empty,
	      fieldset k_sph=fieldset::empty,
	      fieldset r_sph=fieldset::empty)
      falcON_THROWING;
    /// destruction
    ~NBodyCode();
  public:
    //--------------------------------------------------------------------------
    /// describe simulation
    void describe(output&o) const { CODE->describe(o);  }
    /// write N-body data to NEMO snapshot file
    /// \param[in] o stream to NEMO snapshot file
    /// \param[in] w what to write?
    void write(nemo_out&o, fieldset w) const { CODE->write(o,w); }
    /// time integration
    void full_step() { CODE->fullstep(); }
    /// print statistic outputs
    /// \param[in] to ostream to print to
    void stats(output&to) const { 
      if(to) to <<' ';
      CODE->stats_body(to);
    }
    /// header for statistics output
    /// \param[in] to ostream to print to
    void  stats_head(output&to) const {
      if(to) to<<'#';
      CODE->stats_head(to);
      if(to) to<<'#';
      CODE->stats_line(to);
    }
    /// \name data access
    //@{
    /// simulation time of initial snapshot
    double      const&initial_time   () const { return SHOT->initial_time(); }
    /// pointer to body data
    const bodies     *my_bodies      () const { return SHOT; }
    /// pointer to snapshot data
    const snapshot   *my_snapshot    () const { return SHOT; }
    /// input file name
    std::string const&input_file     () const { return FILE; }
    /// current simulation time
    double      const&time           () const { return SHOT->time(); }
    //@}
    //--------------------------------------------------------------------------
  };// class falcON::NBodyCode
#endif // #ifdef falcON_NEMO
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::ForceDiagGrav                                                
  //                                                                            
  /// provides basis diagnostics for gravity only codes                         
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class ForceDiagGrav : public ForceAndDiagnose {
  protected:
    /// \name diagnostics data
    //@{
    const   bool         SELF_GRAV;           ///< compute V_in?
    mutable double       TIME;                ///< time of diagnose
    mutable double       M,T,Vin,Vex,W,TW;    ///< mass, kin & pot E
    mutable vect_d       L;                   ///< total angular momentum
    mutable tensor       KT,WT;               ///< kin & pot energy
    mutable vect_d       CMX,CMV;             ///< center of mass pos & vel
    //@}
    /// perform gravity part of diagnose
    void diagnose_grav() const;
    /// perform velocity part of diagnose
    void diagnose_vels() const falcON_THROWING;
    /// perform full diagnose
    void diagnose_full() const;
    /// construction
    /// \param[in] s  pointer to snapshot to use
    /// \param[in] a  external acceleration field (may be NULL)
    /// \param[in] g  doing self-gravity?
    ForceDiagGrav(snapshot*s, const acceleration*a, bool g)
      : ForceAndDiagnose(s,a), SELF_GRAV(g) {}
  public:
    /// which body data are required?
    virtual fieldset requires() const {
      return fieldset(fieldset::m |
		      fieldset::x |
		      fieldset::v |
		      fieldset::a |
		      fieldset::p |
		      (acc_ext()? fieldset::q : fieldset::empty) );
    }
    /// kinetic energy
    virtual double Ekin() const { return T; }
    /// potential energy
    virtual double Epot() const { return Vin + Vex; }
    /// total energy
    virtual double Etot() const { return Ekin() + Epot(); }
    /// virial ratio
    double const  &Vrat() const { return TW; }
    /// virial W
    double         Wvir() const { return W; }
    /// total angular momentum
    vect_d const  &Ltot() const { return L; }
    /// centre of mass (=mass-averaged) position
    vect_d const  &Xave() const { return CMX; }
    /// centre of mass (=mass-averaged) velocity
    vect_d const  &Vave() const { return CMV; }
    /// print out body of diagnostics statistics
    virtual void dia_stats_body (output&) const;
    /// print out header for diagnostics statistics
    virtual void dia_stats_head (output&) const;
    /// print out line for diagnostics statistics
    virtual void dia_stats_line (output&) const;
    /// do we do self-gravity?
    bool const&self_grav() const { return SELF_GRAV; }
  };// class falcON::ForceDiagGrav
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::GravStepper                                                  
  //                                                                            
  /// used to implement class GravSteps below                                   
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class GravStepper {
    /// enum used in switch() statements
    enum {
      use_a = 1, ///< using criterion  f/acc
      use_p = 2, ///< using criterion  f/pot^2
      use_g = 4, ///< using criterion  f*pot/acc
      use_e = 8  ///< using criterion  f*eps/sqrt(acc)
    };
    typedef BlockStepCode::StepLevels::c_pdouble c_pdouble;
    /// \name data members
    //@{
    bool    UPX;                                 ///< use pot+pex for potential 
    bool    UGE;                                 ///< use global soften'g length
    int     SCH;                                 ///< stepping scheme           
    int     SINKL;                               ///< minimum level for sinks   
    double  EPS;                                 ///< global softening length   
    double  FAQ,FPQ,FGQ,FEQ;                     ///< factors^2 for stepping    
    //@}
    /// full potential (internal and external) of body
    real fpot(body const&B) const {
      return UPX? pex(B)+pot(B) : pot(B);
    }
    /// softening length (individual or global, whatever applies)
    real soft(body const&B) const {
      return UGE? EPS : eps(B);
    }
    /// minimum time-step level of given body
    int minlevel(body const&B) const {
      return is_sink(B)? SINKL : 0;
    }
    //--------------------------------------------------------------------------
  protected:
    /// for a given body compute actual time-step squared
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
    /// assign time-step level to body
    /// \param[in]     Bi body to assign level to
    /// \param[in,out] N  table of number of body per level, to be updated
    /// \param[in]     H  highest time step level
    void assign_level(body&Bi, unsigned*N, int H) const

    {
      double tq=twice(tq_grav(Bi));
      for(Bi.level()=minlevel(Bi); tauq(Bi)>tq && level(Bi)<H; ++(Bi.level()));
      N[level(Bi)]++;
    }
    /// adjust time-step level for given body
    /// \param[in]     Bi body to adjust time-step level for
    /// \param[in,out] N  table of number of body per level, to be updated
    /// \param[in]     L  lowest allowed level at present time
    /// \param[in]     H  highest allowedtime step level
    void adjust_level(body&Bi, unsigned*N, int L, int H) const
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
  public:
    /// constructor
    /// \param[in] sl  minimum time-step level for sinks
    /// \param[in] fa  factor for criterion fa/|acc|
    /// \param[in] fp  factor for criterion fp/pot^2
    /// \param[in] fg  factor for criterion fg*|pot|/|acc|
    /// \param[in] fe  factor for criterion fe*eps/sqrt(|acc|)
    /// \param[in] up  use internal and external potentials
    /// \param[in] ep  global softening length
    /// \param[in] ue  (optional) use global softening length
    /// \note if any of the factors is given as 0, we don't use the
    /// corresponding criterium.
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
    /// full scheme: bit sum of enumeration values.
    int const&scheme() const { return SCH; }
  };// class falcON::GravStepper
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::GravSteps                                                    
  //                                                                            
  /// implements BlockStepCode::StepLevels for gravity-only codes               
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class GravSteps : 
    public BlockStepCode::StepLevels,
    public GravStepper 
  {
    typedef BlockStepCode::StepLevels::c_pdouble c_pdouble;
  protected:
    /// assign time-step level to body
    /// \param[in]     Bi body to assign level to
    /// \param[in,out] N  table of number of body per level, to be updated
    /// \param[in]     H  highest time step level
    void assign_level(body        &Bi,             // I: body                   
		      unsigned    *N,              // I: table: # / step        
		      int          H) const        // I: highest table index    
    { GravStepper::assign_level(Bi,N,H); }
    /// adjust time-step level for given body
    /// \param[in]     Bi body to adjust time-step level for
    /// \param[in,out] N  table of number of body per level, to be updated
    /// \param[in]     L  lowest allowed level at present time
    /// \param[in]     H  highest allowedtime step level
    void adjust_level(body        &Bi,             // I: bodies                 
		      unsigned    *N,              // I: table: # / step        
		      int          L,              // I: lowest  allowed level  
		      int          H) const        // I: highest allowed level  
    { GravStepper::adjust_level(Bi,N,L,H); }
  public:
    /// constructor
    /// \param[in] sl  minimum time-step level for sinks
    /// \param[in] fa  factor for criterion fa/|acc|
    /// \param[in] fp  factor for criterion fp/pot^2
    /// \param[in] fg  factor for criterion fg*|pot|/|acc|
    /// \param[in] fe  factor for criterion fe*eps/sqrt(|acc|)
    /// \param[in] up  use internal and external potentials
    /// \param[in] ep  global softening length
    /// \param[in] ue  (optional) use global softening length
    /// \note if any of the factors is given as 0, we don't use the
    /// corresponding criterium.
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
  //                                                                            
  // class falcON::ForceALCON                                                   
  //                                                                            
  /// standard force solver based on the falcON self-gravity,                   
  /// providing gravity only                                                    
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class ForceALCON : public ForceDiagGrav {
  protected:
    /// \name data are protected
    //@{
#ifdef falcON_INDI
    soft_type         SOFTENING;           ///< type of softening
#ifdef falcON_ADAP
    real              NSOFT;               // #/sphere for eps_i setting
    unsigned          NREF;                // #/cell for eps_i setting  
    real              EMIN;                // lower limit for eps_i     
    real              EFAC;                // max change factor for epsi
#endif
#endif
    const   vect* const ROOTCENTRE;        ///< centre for root, if non-NULL
    const   int         NCRIT, REUSE;      ///< N_crit for tree, re-using
    mutable forces      FALCON;            ///< force algorithm(s)
    mutable int         REUSED;            ///< how often has tree been re-used
    mutable double      CPU_TREE, CPU_GRAV, CPU_AEX; ///< CPU timings
    //@}
    /// build tree and compute forces
    /// \param[in] all   for all bodies (or active only)?
    /// \param[in] build build the tree in any case?
    void set_tree_and_forces(bool all, bool build) const;
  public:
    /// reset softening parameters
    /// \param[in] kern  softening kernel
    /// \param[in] eps   softening length
    void reset_softening(kern_type kern, real eps
#ifdef falcON_ADAP
			,real     ,                // I: Nsoft                  
			 unsigned ,                // I: Nref                   
			 real     ,                // I: eps_min                
			 real                      // I: eps_fac                
#endif
			 );
    /// construction
    /// \param[in] s  snapshot: time & bodies
    /// \param[in] e  global softening length
    /// \param[in] th tree opening parameters
    /// \param[in] nc N_crit for tree build
    /// \param[in] rc pter to root-centre (if non-NULL)
    /// \param[in] k  softening kernel
    /// \param[in] G  Newton's constant of gravity
    /// \param[in] fs theta_sink/theta
    /// \param[in] nr tree re-using (not recommended)
    /// \param[in] ae external acceleration
    /// \param[in] nd direct-summation control for gravity
#ifdef falcON_INDI
    /// \param[in] st softening type
#endif
#ifdef falcON_SPH
    /// \param[in] sd direct-summation control for SPH
#endif
    ForceALCON(snapshot*s, real e, real th, int nc, const vect*rc,
	       kern_type k, real G, real fs, int nr, const acceleration*ae,
	       const int nd[4]
#ifdef falcON_INDI
	       , soft_type st
#ifdef falcON_ADAP
	       , real, unsigned, real, real
#endif
#endif
#ifdef falcON_SPH
	       ,const int sd[3]= Default::SPHdirect
#endif
	       ) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// \name implementing the abstract functions of ForceDiagGrav
    /// \note functions are virtual to allow a class derived from ForceALCON to
    ///       serve in NBodyCode (otherwise the functions below are taken).
    //@{
    /// print header for CPU stats
    virtual void cpu_stats_head(output&) const;
    /// print line for CPU stats
    virtual void cpu_stats_line(output&) const;
    /// print CPU time statistics
    virtual void cpu_stats_body(output&) const;
    /// which fields do we compute
    virtual fieldset computes() const { 
      return fieldset(fieldset::p |
		      fieldset::a |
		      (acc_ext()? fieldset::q : fieldset::empty) );
    }
    /// which fields to we require
    /// \note: we don't require eps_i, but use them as part of force
    ///        computation  
    virtual fieldset requires() const {
      return
	( fieldset(fieldset::m | fieldset::x) | ForceDiagGrav::requires() )
	& ~computes();
    }
    virtual fieldset requiresSPH () const { return fieldset::empty; } 
    virtual fieldset computesSPH () const { return fieldset::empty; } 
    /// compute the forces.
    /// this routine grows the tree; computes self-gravity (acc & pot), unless
    /// self-gravity is not desired; and adds any external gravity, if any.
    /// \param all  set forces for all bodies (or only the active ones)?
    virtual void setforces (bool all, bool, double) const {
      set_tree_and_forces(all,false);
      DebugInfo(8,"ForceALCON::setforces(): done\n");
    }
    /// diagnose: compute total pot & kin energy, etc.
    virtual void diagnose  () const { ForceDiagGrav::diagnose_full(); }
  };// class falcON::ForceALCON
#ifdef falcON_NEMO
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::FalcONCode                                                   
  //                                                                            
  /// N-body code using the falcON force solver; fully inline                   
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class FalcONCode :
    public  NBodyCode,
    private ForceALCON {
  private:
    GravSteps GS;
    //--------------------------------------------------------------------------
    static fieldset to_read(fieldset need, bool soft, bool grav)
    {
      if(soft)
	need += fieldset::e;
      if(grav) {
	need -= fieldset::a;
	need -= fieldset::p;
      }
      return need;
    }
  public:
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
      NBodyCode ( file, resume, to_read(read,
#ifdef falcON_INDI
					soft!=global_fixed ? 1 :
#endif
					0, Grav || aex), time ),
      ForceALCON( SHOT, eps, theta, Ncrit, croot, kernel, Grav, sfac,
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

//-----------------------------------------------------------------------------+
//                                                                             |
// nbdy.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2002                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/nbdy.h>
#include <public/inln.h>
#include <public/Pi.h>
#include <fstream>
#include <iomanip>
#include <ctime>
using std::ostream;
using std::endl;
using std::ios;
using std::setw;
using std::cerr;
using std::setprecision;

using namespace nbdy;
////////////////////////////////////////////////////////////////////////////////
// useful macro looping bodies                                                  
////////////////////////////////////////////////////////////////////////////////
#define LoopMyBodies LoopBodies(bodies,BODIES,Bi)
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  // auxiliary inline function                                                  
  //////////////////////////////////////////////////////////////////////////////
  inline void record_cpu(clock_t& c0, real& CPU) {
    register clock_t c1 = clock();
    CPU += (c1-c0)/real(CLOCKS_PER_SEC);
    c0 = c1;
  }
}
////////////////////////////////////////////////////////////////////////////////
// class nbdy::basic_nbody                                                      
////////////////////////////////////////////////////////////////////////////////
inline real basic_nbody::trace(const tensor&T)
{
  register real tr=zero;
  for(register int i=0; i<NDIM; i++) tr += T[i][i];
  return tr;
}
//------------------------------------------------------------------------------
inline void basic_nbody::reset_cpus() const
{
  CPU_BUILD   = zero;
  CPU_GRAV    = zero;
  CPU_PEX     = zero;
  CPU_DENS    = zero;
  CPU_STEP    = zero;
}
//------------------------------------------------------------------------------
inline void basic_nbody::update_cpu_total() const
{
  register clock_t C_NEW = clock();
  CPU_TOTAL += (C_NEW-C_OLD)/real(CLOCKS_PER_SEC);
  C_OLD = C_NEW;
}
//------------------------------------------------------------------------------
inline void basic_nbody::set_tree(const bool grow_tree)
{
  register clock_t cpu = clock();                  // get current CPU time      
  if(grow_tree)   falcON::grow(NCRIT);             // grow() new tree OR        
  else            falcON::reuse();                 // reuse() old tree          
  record_cpu(cpu,CPU_BUILD);                       // record CPU time required  
}
//------------------------------------------------------------------------------
// - adjusts eps_i of sink bodies       (individual_adaptive)                   
// - upates acceleration of sink bodies                                         
inline void basic_nbody::eps_and_acc(
#ifdef ALLOW_INDI
				     const bool adjust
#endif
				     )
{
  register clock_t cpu = clock();                  // CPU clock                 
#ifdef ALLOW_INDI
  if(SOFTENING==individual_adaptive) {
    if(adjust) falcON::approximate_gravity(true,NSOFT,NREF,EFAC,DIR);
    else       falcON::approximate_gravity(true,NSOFT,NREF,zero,DIR);
  } else 
    falcON::approximate_gravity(true,zero,0,EFAC,DIR);
#else
    falcON::approximate_gravity(true,DIR);
#endif
  record_cpu(cpu,CPU_GRAV);                        // record CPU                
  if(using_extpot()) {                             // IF external potential    >
    register vect A;                               //   external acceleration   
    LoopMyBodies if(is_sink(Bi)) {                 //   loop body sinks        >
      Bi.pot() += PEX->pot_f(A,pos(Bi),time());    //     add external pot      
      Bi.acc() += A;                               //     add external acc      
    }                                              //   <                       
    record_cpu(cpu,CPU_PEX);                       //   record CPU              
  }                                                // <                         
}
//------------------------------------------------------------------------------
void basic_nbody::reset_softening(
				  kern_type const&ker,
				  real      const&e
#ifdef ALLOW_INDI
				  ,
				  real      const&ns,
				  uint      const&nr,
				  real      const&ef
#endif
				  )
{
#ifdef ALLOW_INDI
  NSOFT = ns;
  NREF  = nr;
  EFAC  = abs(ef);
  if(SOFTENING == individual_adaptive && EFAC  == zero)
    NbdyErrorF("using individual adaptive softening, but eps_fac=0",
	       "basic_nbody::reset_softening()")
  if(SOFTENING == individual_adaptive && NSOFT == zero)
    NbdyErrorF("using individual adaptive softening, but Nsoft=0",
	       "basic_nbody::reset_softening()")
#endif
  falcON::reset_softening(abs(e),ker);
}
//------------------------------------------------------------------------------
inline
basic_nbody::basic_nbody(const bodies*const&b,     // I: bodies                 
			 real         const&e,     // I: eps/eps_max            
			 real         const&ti,    // I: t_initial              
			 real         const&th,    // I: tolerance parameter    
			 int          const&nc,    // I: N_crit                 
			 kern_type    const&ker,   // I: softening kernel       
#ifdef ALLOW_INDI
			 soft_type    const&sf,    // I: softening type         
			 real         const&ns,    // I: N_soft                 
			 uint         const&nr,    // I: N_ref                  
			 real         const&ef,    // I: eps_fac                
#endif
			 const extpot*const&p,     // I: P_ex                   
			 const int          nd[4]):// I: direct sum control     
  falcON   ( b,  abs(e), abs(th), ker,
#ifdef ALLOW_INDI
	     sf==global_fixed? global : individual,
#endif
	     th<0? const_theta : theta_of_M ),
  BODIES   ( b ),
  PEX      ( p ),
#ifdef ALLOW_INDI
  SOFTENING( sf ),
#endif
  TINI     ( ti ),
  C_OLD    ( clock() ), 
  NCRIT    ( max(1,nc) ),
  NCUT     ( 0 ), 
  DIAG     ( false ),
  CPU_TOTAL( zero )
{
  DIR[0] = nd[0];
  DIR[1] = nd[1];
  DIR[2] = nd[2];
  DIR[3] = nd[3];
#ifdef ALLOW_INDI
  if(SOFTENING==individual_fixed && !BODIES->has_eps()) 
    NbdyErrorF("individual fixed softening, but no epsi given","basic_nbody")
  reset_softening(ker,e,ns,nr,ef);
#else
  reset_softening(ker,e);
#endif
}
//------------------------------------------------------------------------------
void basic_nbody::reset_opening(const real th) const {
  falcON::reset_opening(abs(th),
			th<zero? const_theta : theta_of_M);
}
//------------------------------------------------------------------------------
#ifdef ALLOW_INDI
inline void basic_nbody::estimate_mass_density(const bool ext)
{
  register clock_t  cpu = clock();                 // cpu time                  
  falcON::estimate_rho(NREF);                      // estimate rho in tree      
  if(PEX && !PEX->is_empty() && ext)               // if external potential     
    LoopMyBodies if(is_sink(Bi))                   //   loop body sinks         
      Bi.rho()+=PEX->rho(pos(Bi),time());          //     add ext density       
  record_cpu(cpu,CPU_DENS);                        // record CPU                
}
//------------------------------------------------------------------------------
void basic_nbody::estimate_mass_densities(const bool ext)
{
  LoopMyBodies Bi.flag_as_sink();
  estimate_mass_density(ext);
}
//------------------------------------------------------------------------------
void basic_nbody::estimate_surf_densities()
{
  LoopMyBodies Bi.flag_as_sink();
  falcON::estimate_sd(NREF);
}
#endif
//------------------------------------------------------------------------------
void basic_nbody::diagnose() const
#define LoopTensor for(i=0; i<NDIM; i++) for(j=0; j<NDIM; j++)
{
  register indx   i,j;
  register double u=zero, m=zero, t=zero, ue=zero;
  register vect   x=zero, v=zero;
  register amom   l=zero;
  register tensor T,W;
  LoopTensor T[i][j] = W[i][j] = zero;
  if(PEX && !PEX->is_empty()) {
    register real uer=zero, mr=zero;
    PEX->energies(uer,W,T,x,v,l,mr);
    ue = uer;
    m  = mr;
  }
  LoopMyBodies {
    m += double(mass(Bi));                         // total mass                
    u += double(mass(Bi)) * pot(Bi);               // internal potential energy 
    t += double(mass(Bi)) * 
      ( square(double(vel(Bi)[0])) +
        square(double(vel(Bi)[1])) +
        square(double(vel(Bi)[2])) );
    LoopTensor {
      T[i][j] += mass(Bi) * vel(Bi)[i]*vel(Bi)[j]; // kin energy tensor         
      W[i][j] += mass(Bi) * pos(Bi)[i]*acc(Bi)[j]; // pot energy tensor         
    }
    x += mass(Bi) * pos(Bi);                       // dipole                    
    v += mass(Bi) * vel(Bi);                       // momentum                  
    l += mass(Bi) * pos(Bi) ^ vel(Bi);             // total angular momentum    
  }
  M    = m;                                        // total mass                
  Ktot = half*t;                                   // total kin energy          
  Utot = half*u+ue;                                // total pot energy          
  TU   =-half*trace(T)/trace(W);                   // virial ratio              
  L    = l;                                        // total angular momentum    
  CMX  = x/m;                                      // center of mass            
  CMV  = v/m;                                      // center of mass velocity   
  LoopTensor {
    KT[i][j] = half * T[i][j];                     // kin energy tensor         
    WT[i][j] =        W[i][j];                     // pos energy tensor         
  }
}
#undef LoopTensor
////////////////////////////////////////////////////////////////////////////////
// class fasstree::LeapFrogCode                                                 
////////////////////////////////////////////////////////////////////////////////
LeapFrogCode::LeapFrogCode(const bodies*const&b,   // I: bodies                 
			   real         const&e,   // I: eps/eps_max            
			   real         const&ti,  // I: t_initial              
			   int          const&h0,  // I: h0                     
			   int          const&hg,  //[I: h_grow]                
			   real         const&th,  //[I: tolerance parameter]   
			   int          const&nc,  //[I: N_crit]                
			   kern_type    const&ker, //[I: softening kernel]      
#ifdef ALLOW_INDI
			   soft_type    const&sf,  //[I: softening type]        
			   real         const&ns,  //[I: N_soft]                
			   uint         const&nr,  //[I: N_ref]                 
			   real         const&ef,  //[I: eps_fac]               
#endif
			   const extpot*const&p,   //[I: P_ex]                  
			   const int       nd[4]): //[I: direct sum control]    
  basic_nbody(b,e,ti,th,nc,ker,
#ifdef ALLOW_INDI
	      sf,ns,nr,ef,
#endif
	      p,nd),
  LeapFrog   (h0,ti),
  REUSE      ((1<<hg)-1),
  REUSED     (0u)
{
  register clock_t  cpu = clock();                 // cpu time                  
  reset_cpus();                                    // reset CPU time counters   
  set_tree(true);                                  // create tree               
  LoopMyBodies Bi.flag_as_sink();                  // flag everybody as sink    
  eps_and_acc(                                     // set eps & get acc         
#ifdef ALLOW_INDI
	      false
#endif
	      );
  record_cpu(cpu,CPU_STEP);                        // record CPU time           
  update_cpu_total();                              // record accumulated time   
}
//------------------------------------------------------------------------------
void LeapFrogCode::full_step()
{
  register clock_t  cpu = clock();                 // cpu time                  
  reset_cpus();                                    // reset CPU time counters   
  predict(BODIES);                                 // move = drift              
  if(REUSED < REUSE) {                             // IF(re-using old tree)    >
    set_tree(false);                               //   create tree             
    REUSED++;                                      //   count # re-using        
  } else {                                         // < ELSE                   >
    set_tree(true);                                //   create tree             
    REUSED=0;                                      //   reset # re-using        
  }                                                // <                         
  eps_and_acc();                                   // asjust eps & get acc      
  accelerate(BODIES);                              // accelerate = kick         
  DIAG = false;                                    // diagnose out-of-date      
  record_cpu(cpu,CPU_STEP);                        // record CPU time           
  update_cpu_total();                              // record accumulated time   
}
//------------------------------------------------------------------------------
void LeapFrogCode::stats(ostream& to) const
{
  to.setf(ios::left, ios::adjustfield);
  update_diagnostics();
  to<<setprecision(5)<<setw(10)<<time()         <<" "
    <<setprecision(7)<<setw(13)<<total_energy() <<" "
    <<setprecision(min(5,max(1,5-int(-log10(virial_ratio())))))
    <<setw(7)<<virial_ratio()<<" ";
#if NDIM==3
  to<<setprecision(5)<<setw(10)<<sqrt(norm(total_angmom()))  <<" ";
#else
  to<<setprecision(5)<<setw(10)<<total_angmom()              <<" ";
#endif
  to<<setprecision(2)<<setw(8)<<sqrt(norm(total_momentum()))  <<" "
    <<setprecision(4)<<setw(6)<<cpu_build()    <<" "
    <<setprecision(4)<<setw(7)<<cpu_force()    <<" ";
  to.setf(ios::right, ios::adjustfield);
  to<<setprecision(4)<<setw(7)<<cpu_longstep() <<" "
    <<setprecision(6)<<setw(9)<<cpu_total()    <<endl;
}
//------------------------------------------------------------------------------
void LeapFrogCode::stats_head(ostream& to) const
{
  to<<"    time   "
    <<"   energy     "
    <<"  -T/U  "
#if NDIM==3
    <<"   |L|     "
#else
    <<"   L_z     "
#endif
    <<"  |v_cm| "
    <<" build   force    step     accum"<<endl;
}
//------------------------------------------------------------------------------
void LeapFrogCode::stats_line(ostream &to) const
{
  to<<" ----------------------------------------------------"
    "--------------------------------"<<endl;
}
////////////////////////////////////////////////////////////////////////////////
// class nbdy::BlockStepCode                                                    
////////////////////////////////////////////////////////////////////////////////
void BlockStepCode::reset_stepping(const real fa, const real fp, const real fc)
{
  BlockStep::reset_scheme(fa,fp,fc);
}
//------------------------------------------------------------------------------
#define ForAll       LoopMyBodies
#define ForActive    LoopMyBodies if(is_sink(Bi))
//------------------------------------------------------------------------------
inline void BlockStepCode::set_sink_flags(const indx low)
{
  LoopMyBodies                                     // loop bodies in tree       
    if(level(Bi) >= low) Bi.flag_as_sink();        //   IF(level>=low): sink    
    else                 Bi.unflag_sink ();        //   ELSE          : no sink 
}
//------------------------------------------------------------------------------
void BlockStepCode::
elementary_step(const indx low)                    // I: lowest level moving    
{
  static uint m = 0u;                              // # of tiny moves omitted   
  clock_on();                                      // increase actual clock     
  if(need_move(low)) {                             // IF anybody is sink now   >
    register real dt = tau_min() * (m+1);          //   dt = (m+1)*tau_min      
    ForAll move_by(Bi,dt);                         //   all:    x -> x+v*dt     
    m = 0;                                         //   reset m = 0             
    set_tree(low<LGROW);                           //   initialize tree         
    set_sink_flags(low);                           //   IF(l>=low) -> active    
    eps_and_acc();                                 //   active: a = acc(x)      
    ForActive acce_half(Bi);                       //   active: v -> v+a*h/2    
    if(low != highest_level())                     //   IF(levels may change)  >
      ForActive adjust_level(Bi,low);              //     active: h -> h'      <
    if(low)                                        //   IF(not last step)      >
      ForActive acce_half(Bi);                     //     active: v -> v+a*h/2 <
  } else m++;                                      // < ELSE count move omission
}
//------------------------------------------------------------------------------
void BlockStepCode::full_step()
{
  register clock_t  cpu = clock();                 // cpu time                  
  const    unsigned nn  = 1 << highest_level();    // number of small steps     
  register unsigned tn;                            // index of small step       
  reset_cpus();                                    // reset CPU time counters   
  ForAll acce_half(Bi);                            // v -> v + a*h/2            
  for(tn=1; tn<nn; tn++)                           // loop over all but last    
    elementary_step(longest_moving(tn));           //   elementary steps        
  elementary_step(0);                              // last step: grow tree      
  DIAG = false;                                    // diagnose out-of-date      
  record_cpu(cpu,CPU_STEP);                        // record CPU time           
  update_cpu_total();                              // record accumulated time   
}
#undef ForAll
#undef ForActive
//------------------------------------------------------------------------------
void BlockStepCode::prepare(const int h0,          // I: h0                     
			    const int nl)          // I: N_levels               
{
  register clock_t cpu = clock();                  // local CPU counter         
  reset_cpus();                                    // reset CPU counters        
  BlockStep::reset_steps(h0,nl,TINI);              // set time steps            
  set_tree(true);                                  // initialize tree           
  LoopMyBodies Bi.flag_as_sink();                  // flag every body as sink   
  eps_and_acc(                                     // set eps & get acc         
#ifdef ALLOW_INDI
	      false
#endif
	      );
  diagnose();                                      // energy et al. at t=t0     
  DIAG = true;                                     // diagnose up-to-date       
  BlockStep::reset_counts();                       // reset counts              
  LoopMyBodies                                     // loop bodies               
    assign_level(Bi);                              //   assign step level       
  record_cpu(cpu,CPU_STEP);                        // record total CPU          
  update_cpu_total();                              // record accumulated time   
}
//------------------------------------------------------------------------------
BlockStepCode::BlockStepCode(const bodies*const&b,   // I: bodies               
			     real         const&e,   // I: eps/eps_max          
			     real         const&ti,  // I: t_initial            
			     int          const&h0,  // I: h0                   
			     int          const&nl,  // I: # levels             
			     real         const&fa,  // I: f_a: for stepping    
			     real         const&fp,  //[I: f_p: for stepping]   
			     real         const&fc,  //[I: f_c: for stepping]   
			     int          const&hg,  //[I: h_grow]              
			     real         const&th,  //[I: tolerance parameter] 
			     int          const&nc,  //[I: N_crit]              
			     kern_type    const&ker, //[I: softening kernel]    
#ifdef ALLOW_INDI
			     soft_type    const&sf,  //[I: softening type]      
			     real         const&ns,  //[I: N_soft]              
			     uint         const&nr,  //[I: N_ref]               
			     real         const&ef,  //[I: eps_fac]             
#endif
			     const extpot*const&p,   //[I: P_ex]                
			     const int       nd[4]): //[I: direct sum control]  
  basic_nbody(b,e,ti,th,nc,ker,
#ifdef ALLOW_INDI
	      sf,ns,nr,ef,
#endif
	      p,nd),
  BlockStep  (),
  LGROW      ((nl <= hg)? 1 : nl-hg),
  W          (max(6, int(log10(double(b->N_bodies())))))
{
  BlockStep::reset_scheme(fa,fp,fc);
  if(!BODIES->has_lev())
    NbdyErrorF("bodies::has_lev() == false","BlockStepCode")
  if(fa==zero && fp==zero) 
    NbdyErrorF("time step control factors = 0","BlockStepCode")
  prepare(h0,nl);
}
//------------------------------------------------------------------------------
void BlockStepCode::stats(ostream& to) const
{
  to.setf(ios::left, ios::adjustfield);
  if(time() != TINI) update_diagnostics();
  to<<setprecision(5)<<setw(10)<<time()         <<" "
    <<setprecision(7)<<setw(13)<<total_energy() <<" "
    <<setprecision(min(5,max(1,5-int(-log10(virial_ratio())))))
    <<setw(7)<<virial_ratio()<<" ";
#if NDIM==3
  to<<setprecision(5)<<setw(10)<<sqrt(norm(total_angmom()))  <<" ";
#else
  to<<setprecision(5)<<setw(10)<<total_angmom()              <<" ";
#endif
  to<<setprecision(2)<<setw(8)<<sqrt(norm(total_momentum()))  <<" ";
  to.setf(ios::right, ios::adjustfield);
  if(highest_level()) BlockStep::short_stats(to,W);
  to.setf(ios::left, ios::adjustfield);
  to<<setprecision(4)<<setw(6)<<cpu_build()    <<" "
    <<setprecision(4)<<setw(7)<<cpu_force()    <<" "
    <<setprecision(4)<<setw(7)<<cpu_longstep() <<" ";
  to.setf(ios::right, ios::adjustfield);
  to<<setprecision(6)<<setw(9)<<cpu_total()    <<endl;
}
//------------------------------------------------------------------------------
namespace nbdy {
  inline ostream& put_char(ostream&o, const char c, const int s)
  {
    if(s>0) for(register int i=0; i<s; ++i) o<<c;
    return o;
  }
}
//------------------------------------------------------------------------------
void BlockStepCode::stats_head(ostream& to) const
{
  to<<"    time   "
    <<"   energy     "
    <<"  -T/U  "
#if NDIM==3
    <<"   |L|     "
#else
    <<"    L      "
#endif
    <<"  |v_cm| ";
  if(highest_level())
    for(register int i=0, h=-h0(); i<=highest_level(); i++, h--)
           if(h>13)  put_char(to,' ',W-5)<<" 2^" <<     h  <<" ";
      else if(h> 9)  put_char(to,' ',W-5)<<" "   << (1<<h) <<" ";
      else if(h> 6)  put_char(to,' ',W-5)<<"  "  << (1<<h) <<" ";
      else if(h> 3)  put_char(to,' ',W-5)<<"   " << (1<<h) <<" ";
      else if(h>=0)  put_char(to,' ',W-5)<<"    "<< (1<<h) <<" ";
      else if(h==-1) put_char(to,' ',W-5)<<"  1/2 ";
      else if(h==-2) put_char(to,' ',W-5)<<"  1/4 ";
      else if(h==-3) put_char(to,' ',W-5)<<"  1/8 ";
      else if(h==-4) put_char(to,' ',W-5)<<" 1/16 ";
      else if(h==-5) put_char(to,' ',W-5)<<" 1/32 ";
      else if(h==-6) put_char(to,' ',W-5)<<" 1/64 ";
      else if(h>-10) put_char(to,' ',W-5)<<" 2^" <<     h  <<" ";
      else           put_char(to,' ',W-5)<<"2^"  <<     h  <<" ";
  to<<" build   force    step     accum"<<endl;
}
//------------------------------------------------------------------------------
void BlockStepCode::stats_line(ostream &to) const
{
  to<<" ----------------------------------------------------";
  if(highest_level())
    for(register int i=0; i<=highest_level(); i++) put_char(to,'-',W+1);
  to<<"--------------------------------"<<endl;
}
//------------------------------------------------------------------------------
#undef LoopMyBodies

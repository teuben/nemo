//-----------------------------------------------------------------------------+
//                                                                             |
// nbdy.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
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
  //                                                                          //
  // auxiliary inline function                                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  inline void record_cpu(clock_t& c0, real& CPU) {
    register clock_t c1 = clock();
    CPU += (c1-c0)/real(CLOCKS_PER_SEC);
    c0 = c1;
  }
  //----------------------------------------------------------------------------
  template<typename scalar, int N>
  inline scalar tr(scalar T[N][N])
  {
    register scalar x=scalar(0);
    for(register int i=0; i<N; i++) x += T[i][i];
    return x;
  }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::basic_nbody                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
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
  report REPORT("basic_nbody::set_tree(%d)",grow_tree);
  register clock_t cpu = clock();                  // get current CPU time      
  if(grow_tree)   falcON::grow(NCRIT);             // grow() new tree OR        
  else            falcON::reuse();                 // reuse() old tree          
  record_cpu(cpu,CPU_BUILD);                       // record CPU time required  
}
//------------------------------------------------------------------------------
// - adjusts eps_i of active bodies       (individual_adaptive)                 
// - upates acceleration of active bodies                                       
inline void basic_nbody::eps_and_acc(bool const&all
#ifdef falcON_INDI
				    ,bool const&adjust
#endif
				     )
{
#ifdef falcON_INDI
  report REPORT("basic_nbody::eps_and_acc(%d,%d)",all,adjust);
#else
  report REPORT("basic_nbody::eps_and_acc(%d)",all);
#endif
  register clock_t cpu = clock();                  // CPU clock                 
#ifdef falcON_INDI
  if(SOFTENING==individual_adaptive) {
    if(adjust) falcON::approximate_gravity(true,all,NSOFT,NREF,EMIN,EFAC,DIR);
    else       falcON::approximate_gravity(true,all,NSOFT,NREF,zero,zero,DIR);
  } else 
    falcON::approximate_gravity(true,all,zero,0u,EMIN,zero,DIR);
#else
    falcON::approximate_gravity(true,all,DIR);
#endif
  record_cpu(cpu,CPU_GRAV);                        // record CPU                
  if(using_extpot()) {                             // IF external potential     
    report REPORT2("adding external potential");
    register vect A;                               //   external acceleration   
    LoopMyBodies if(is_active(Bi)) {               //   LOOP active bodies      
      Bi.pex()  = PEX->pot_f(A,pos(Bi),time());    //     get external pot      
      Bi.acc() += A;                               //     add external acc      
    }                                              //   END LOOP                
    record_cpu(cpu,CPU_PEX);                       //   record CPU              
  }                                                // ENDIF                     
}
//------------------------------------------------------------------------------
void basic_nbody::reset_softening(
				  kern_type const&ker,
				  real      const&e
#ifdef falcON_INDI
				 ,real      const&ns,
				  uint      const&nr,
				  real      const&em,
				  real      const&ef
#endif
				  )
{
#ifdef falcON_INDI
  NSOFT = ns;
  NREF  = nr;
  EMIN  = abs(em);
  EFAC  = abs(ef);
  if(SOFTENING == individual_adaptive && EFAC  == zero)
    falcON_ErrorF("using individual adaptive softening, but eps_fac=0",
		  "basic_nbody::reset_softening()");
  if(SOFTENING == individual_adaptive && NSOFT == zero)
    falcON_ErrorF("using individual adaptive softening, but Nsoft=0",
		  "basic_nbody::reset_softening()");
#endif
  falcON::reset_softening(abs(e),ker);
}
//------------------------------------------------------------------------------
#ifdef TESTING_
static std::ofstream test[10];
#endif
inline
basic_nbody::basic_nbody(const bodies*const&b,     // I: bodies                 
			 real         const&e,     // I: eps/eps_max            
			 real         const&ti,    // I: t_initial              
			 real         const&th,    // I: tolerance parameter    
			 int          const&nc,    // I: N_crit                 
			 kern_type    const&ker,   // I: softening kernel       
#ifdef falcON_INDI
			 soft_type    const&sf,    // I: softening type         
			 real         const&ns,    // I: N_soft                 
			 uint         const&nr,    // I: N_ref                  
			 real         const&em,    // I: eps_min                
			 real         const&ef,    // I: eps_fac                
#endif
			 const extpot*const&p,     // I: P_ex                   
			 const int          nd[4]):// I: direct sum control     
  falcON   ( b,  abs(e), abs(th), ker,
#ifdef falcON_INDI
	     sf==global_fixed? global : individual,
#endif
	     th<0? const_theta : theta_of_M ),
  BODIES   ( b ),
  PEX      ( p ),
#ifdef falcON_INDI
  SOFTENING( sf ),
#endif
  TINI     ( ti ),
  C_OLD    ( clock() ), 
  NCRIT    ( max(1,nc) ),
  NCUT     ( 0 ), 
  DIAG     ( false ),
  CPU_TOTAL( zero )
{
#ifdef TESTING_
  // TEST
  char file[10];
  for(register int i=0; i!=10; ++i) {
    sprintf(file,"test.%03d",i+1);
    test[i].open(file);
  }
  // TSET
#endif
  DIR[0] = nd[0];
  DIR[1] = nd[1];
  DIR[2] = nd[2];
  DIR[3] = nd[3];
#ifdef falcON_INDI
  if(SOFTENING==individual_fixed && !BODIES->has(io::e)) 
    falcON_ErrorF("individual fixed softening, but no epsi given",
		  "basic_nbody");
  reset_softening(ker,e,ns,nr,em,ef);
#else
  reset_softening(ker,e);
#endif
  if(PEX != 0 && !BODIES->has(io::P))
    falcON_ErrorF("external potential desired, but nobody has memory",
		  "basic_nbody");
}
//------------------------------------------------------------------------------
void basic_nbody::reset_opening(const real th) const {
  falcON::reset_opening(abs(th),
			th<zero? const_theta : theta_of_M);
}
//------------------------------------------------------------------------------
#ifdef falcON_INDI
inline void basic_nbody::estimate_mass_density(const bool ext)
{
  register clock_t  cpu = clock();                 // cpu time                  
  falcON::estimate_rho(NREF);                      // estimate rho in tree      
  if(PEX && !PEX->is_empty() && ext)               // IF external potential     
    LoopMyBodies if(is_active(Bi))                 //   LOOP active bodies      
      Bi.rho()+=PEX->rho(pos(Bi),time());          //     add ext density       
  record_cpu(cpu,CPU_DENS);                        // record CPU                
}
//------------------------------------------------------------------------------
void basic_nbody::estimate_mass_densities(const bool ext)
{
  LoopMyBodies Bi.flag_as_active();
  estimate_mass_density(ext);
}
//------------------------------------------------------------------------------
void basic_nbody::estimate_surf_densities()
{
  LoopMyBodies Bi.flag_as_active();
  falcON::estimate_sd(NREF);
}
#endif
//------------------------------------------------------------------------------
void basic_nbody::diagnose() const
#define LoopTensor							\
  for(register int i=0; i!=Ndim; ++i) for(register int j=0; j!=Ndim; ++j)
{
  report REPORT("basic_nbody::diagnose()");
  register double  m=0., ui=0., ue=0.;
  register vect_d  x=0., v=0., mx, mv;
  register amom_d  l=0.;
  register double  K[Ndim][Ndim]={0.}, W[Ndim][Ndim]={0.};
  if(PEX && !PEX->is_empty()) {                    // IF have external pot      
    PEX->energies(ue,K,W,x,v,l,m);                 //   get external contrib    
    LoopMyBodies {                                 //   LOOP bodies             
      m += mass(Bi);                               //     add: total mass       
      ui+= mass(Bi) * pot(Bi);                     //     add: int pot energy   
      ue+= mass(Bi) * pex(Bi);                     //     add: ext pot energy   
      mx = mass(Bi) * pos(Bi);                     //     m * x                 
      mv = mass(Bi) * vel(Bi);                     //     m * v                 
      LoopTensor {                                 //     LOOP i,j              
	K[i][j] += mv[i] * vel(Bi)[j];             //       add: K_ij           
	W[i][j] += mx[i] * acc(Bi)[j];             //       add: W_ij           
      }                                            //     END LOOP              
      x += mx;                                     //     add: dipole           
      v += mv;                                     //     add: total momentum   
      l += mx ^ vel(Bi);                           //     add: total ang mom    
    }                                              //   END LOOP                
  } else {                                         // ELSE: no external pot     
    LoopMyBodies {                                 //   LOOP bodies             
      m += mass(Bi);                               //     add: total mass       
      ui+= mass(Bi) * pot(Bi);                     //     add: int pot energy   
      mx = mass(Bi) * pos(Bi);                     //     m * x                 
      mv = mass(Bi) * vel(Bi);                     //     m * v                 
      LoopTensor {                                 //     LOOP i,j              
	K[i][j] += mv[i] * vel(Bi)[j];             //       add: K_ij           
	W[i][j] += mx[i] * acc(Bi)[j];             //       add: W_ij           
      }                                            //     END LOOP              
      x += mx;                                     //     add: dipole           
      v += mv;                                     //     add: total momentum   
      l += mx ^ vel(Bi);                           //     add: total ang mom    
    }                                              //   END LOOP                
  }                                                // ENDIF                     
  M    = m;                                        // total mass                
  Ktot = 0.5*tr(K);                                // total kin energy          
  Uin  = 0.5*ui;                                   // total int pot energy      
  Uex  = ue;                                       // total ext pot energy      
  TU   =-half*tr(K)/tr(W);                         // virial ratio              
  L    = l;                                        // total angular momentum    
  m    = 1./m;
  CMX  = x*m;                                      // center of mass            
  CMV  = v*m;                                      // center of mass velocity   
  LoopTensor {
    KT[i][j] = half * K[i][j];                     // kin energy tensor         
    WT[i][j] = half *(W[i][j]+W[j][i]);            // pot energy tensor         
  }
}
#undef LoopTensor
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::LeapFrogCode                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
LeapFrogCode::LeapFrogCode(const bodies*const&b,   // I: bodies                 
			   real         const&e,   // I: eps/eps_max            
			   real         const&ti,  // I: t_initial              
			   int          const&h0,  // I: h0                     
			   int          const&hg,  //[I: h_grow]                
			   real         const&th,  //[I: tolerance parameter]   
			   int          const&nc,  //[I: N_crit]                
			   kern_type    const&ker, //[I: softening kernel]      
#ifdef falcON_INDI
			   soft_type    const&sf,  //[I: softening type]        
			   real         const&ns,  //[I: N_soft]                
			   uint         const&nr,  //[I: N_ref]                 
			   real         const&em,  //[I: eps_min]               
			   real         const&ef,  //[I: eps_fac]               
#endif
			   const extpot*const&p,   //[I: P_ex]                  
			   const int       nd[4]): //[I: direct sum control]    
  basic_nbody(b,e,ti,th,nc,ker,
#ifdef falcON_INDI
	      sf,ns,nr,em,ef,
#endif
	      p,nd),
  LeapFrog   (h0,ti),
  REUSE      ((1<<hg)-1),
  REUSED     (0u)
{
  register clock_t  cpu = clock();                 // cpu time                  
  reset_cpus();                                    // reset CPU time counters   
  set_tree(true);                                  // create tree               
  LoopMyBodies Bi.flag_as_active();                // flag everybody as active  
  eps_and_acc(true                                 // set eps & get acc         
#ifdef falcON_INDI
	     ,false
#endif
	      );
  record_cpu(cpu,CPU_STEP);                        // record CPU time           
  update_cpu_total();                              // record accumulated time   
}
//------------------------------------------------------------------------------
void LeapFrogCode::full_step()
{
  report REPORT("LeapFrogCode::full_step()");
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
#ifdef TESTING_
  // TEST
  for(register int i=0; i!=10; ++i)
    test[i] << time() <<"  "
	    << BODIES->pos(i) << "  "
	    << BODIES->vel(i) << "  "
	    << BODIES->acc(i) << "  "
	    << BODIES->pot(i) << "  "
	    << BODIES->eps(i) << std::endl;
  // TSET
#endif
  eps_and_acc(true);                               // asjust eps & get acc      
  accelerate(BODIES);                              // accelerate = kick         
  DIAG = false;                                    // diagnose out-of-date      
  record_cpu(cpu,CPU_STEP);                        // record CPU time           
  update_cpu_total();                              // record accumulated time   
}
//------------------------------------------------------------------------------
void LeapFrogCode::stats(ostream& to) const
{
  report REPORT("LeapFrogCode::stats()");
  to.setf(ios::left, ios::adjustfield);
  update_diags();
  to<<setprecision(5)<<setw(10)<<time()         <<" "
    <<setprecision(7)<<setw(13)<<total_energy() <<" "
    <<setprecision(min(5,max(1,5-int(-log10(virial_ratio())))))
    <<setw(7)<<virial_ratio()<<" ";
#if falcON_NDIM==3
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
#if falcON_NDIM==3
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
//                                                                            //
// class nbdy::BlockStepCode                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void BlockStepCode::reset_stepping(real const&fa,
				   real const&fp,
				   real const&fc,
				   real const&fe)
{
  GravBlockStep::reset_scheme(fa,fp,fc,fe);
}
//------------------------------------------------------------------------------
#define ForAll       LoopMyBodies
#define ForActive    LoopMyBodies if(is_active(Bi))
//------------------------------------------------------------------------------
inline void BlockStepCode::set_active_flags(int const&low)
{
  report REPORT("BlockStepCode::set_active_flags(%d)",low);
  LoopMyBodies                                     // LOOP bodies in tree       
    if(level(Bi) >= low) Bi.flag_as_active();      //   IF(level>=low): active  
    else                 Bi.unflag_active ();      //   ELSE          : inactive
}
//------------------------------------------------------------------------------
void BlockStepCode::
elementary_step(int const&low)                     // I: lowest level moving    
{
  report REPORT("BlockStepCode::elementary_step(%d)",low);
#ifdef TESTING_
  // TEST
  for(register int i=0; i!=10; ++i) if(low==0 || is_active(BODIES->flg(i)))
    test[i] << time() <<"  "
	    << BODIES->pos(i) << "  "
	    << BODIES->vel(i) << "  "
	    << BODIES->acc(i) << "  "
	    << BODIES->pot(i) << "  "
	    << BODIES->eps(i) << "  "
	    << BODIES->level(i) << std::endl;
  // TSET
#endif
  static uint m = 0u;                              // # of tiny moves omitted   
  clock_on();                                      // increase actual clock     
  if(need_move(low)) {                             // IF anybody is active now  
    register real dt = tau_min() * (m+1);          //   dt = (m+1)*tau_min      
    ForAll move_by(Bi,dt);                         //   all:    x -> x+v*dt     
    m = 0;                                         //   reset m = 0             
    set_tree(low<LGROW);                           //   initialize tree         
    set_active_flags(low);                         //   IF(l>=low) -> active    
    eps_and_acc(low==0);                           //   active: a = acc(x)      
    ForActive acce_half(Bi);                       //   active: v -> v+a*h/2    
    if(low != highest_level())                     //   IF(levels may change)   
#ifdef falcON_INDI
      if(SOFTENING==global_fixed) {                //     IF fixed eps          
	ForActive adjust_level(Bi,low,eps());      //       active: h -> h'     
      } else                                       //     ELSE: individual eps_i
#endif
	ForActive adjust_level(Bi,low,::eps(Bi));  //       active: h -> h'     
    if(low)                                        //   IF(not last step)       
      ForActive acce_half(Bi);                     //     active: v -> v+a*h/2  
  } else m++;                                      // ELSE count move omission  
}
//------------------------------------------------------------------------------
void BlockStepCode::full_step()
{
  report REPORT("BlockStepCode::full_step()");
  register clock_t  cpu = clock();                 // cpu time                  
  const    unsigned nn  = 1 << highest_level();    // number of small steps     
  register unsigned tn;                            // index of small step       
  reset_cpus();                                    // reset CPU time counters   
  ForAll acce_half(Bi);                            // v -> v + a*h/2            
  for(tn=1; tn<nn; tn++)                           // LOOP all but last         
    elementary_step(longest_moving(tn));           //   elementary steps        
  elementary_step(0);                              // last step: grow tree      
  DIAG = false;                                    // diagnose out-of-date      
  record_cpu(cpu,CPU_STEP);                        // record CPU time           
  update_cpu_total();                              // record accumulated time   
}
#undef ForAll
#undef ForActive
//------------------------------------------------------------------------------
void BlockStepCode::prepare(int const&h0,          // I: h0                     
			    int const&nl)          // I: N_levels               
{
  report REPORT("BlockStepCode::prepare(%d,%d)",h0,nl);
  register clock_t cpu = clock();                  // local CPU counter         
  reset_cpus();                                    // reset CPU counters        
  GravBlockStep::reset_steps(h0,nl,TINI);          // set time steps            
  set_tree(true);                                  // initialize tree           
  LoopMyBodies Bi.flag_as_active();                // flag every body as active 
  eps_and_acc(true                                 // set eps & get acc         
#ifdef falcON_INDI
	     ,false
#endif
	      );
  diagnose();                                      // energy et al. at t=t0     
  DIAG = true;                                     // diagnose up-to-date       
  GravBlockStep::reset_counts();                   // reset counts              
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
			     real         const&fe,  //[I: f_e: for stepping]   
			     int          const&hg,  //[I: h_grow]              
			     real         const&th,  //[I: tolerance parameter] 
			     int          const&nc,  //[I: N_crit]              
			     kern_type    const&ker, //[I: softening kernel]    
#ifdef falcON_INDI
			     soft_type    const&sf,  //[I: softening type]      
			     real         const&ns,  //[I: N_soft]              
			     uint         const&nr,  //[I: N_ref]               
			     real         const&em,  //[I: eps_min]             
			     real         const&ef,  //[I: eps_fac]             
#endif
			     const extpot*const&p,   //[I: P_ex]                
			     const int       nd[4]): //[I: direct sum control]  
  basic_nbody  (b,e,ti,th,nc,ker,
#ifdef falcON_INDI
	        sf,ns,nr,em,ef,
#endif
	        p,nd),
  GravBlockStep (),
  LGROW         ((nl <= hg)? 1 : nl-hg),
  W             (max(6, 1+int(log10(double(b->N_bodies())))))
{
  GravBlockStep::reset_scheme(fa,fp,fc,fe);
  if(!BODIES->has(io::l))
    falcON_ErrorF("bodies::has(io::l) == false","BlockStepCode");
  prepare(h0,nl);
}
//------------------------------------------------------------------------------
void BlockStepCode::stats(ostream& to) const
{
  report REPORT("BlockStepCode::stats()");
  to.setf(ios::left, ios::adjustfield);
  if(time() != TINI) update_diags();
  to<<setprecision(5)<<setw(10)<<time()         <<" "
    <<setprecision(7)<<setw(13)<<total_energy() <<" "
    <<setprecision(min(5,max(1,5-int(-log10(virial_ratio())))))
    <<setw(7)<<virial_ratio()<<" ";
#if falcON_NDIM==3
  to<<setprecision(5)<<setw(10)<<sqrt(norm(total_angmom()))  <<" ";
#else
  to<<setprecision(5)<<setw(10)<<total_angmom()              <<" ";
#endif
  to<<setprecision(2)<<setw(8)<<sqrt(norm(total_momentum()))  <<" ";
  to.setf(ios::right, ios::adjustfield);
  if(highest_level()) GravBlockStep::short_stats(to,W);
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
#if falcON_NDIM==3
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

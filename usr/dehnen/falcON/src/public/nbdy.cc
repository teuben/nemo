//-----------------------------------------------------------------------------+
//                                                                             |
// nbdy.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <nbdy.h>
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
using std::setfill;

using namespace nbdy;
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// useful macro looping bodies                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#define LoopMyBodies       LoopBodies(bodies,BODIES,Bi)
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::basic_nbody                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void basic_nbody::set_gravity(bool const&grow_tree,
			      bool const&all
#ifdef falcON_ADAP
			     ,bool const&adjust
#endif
			      )
  // - adjusts eps_i of active bodies       (individual_adaptive)               
  // - upates acceleration of active bodies                                     
  // we deal with four cases, depending on whether or not self-gravity is to be 
  // computed and whether or not external gravity is to be added                
{
#ifdef falcON_ADAP
  report REPORT("basic_nbody::set_gravity(%d,%d,%d)",grow_tree,all,adjust);
#else
  report REPORT("basic_nbody::set_gravity(%d,%d)",grow_tree,all);
#endif
  register clock_t cpu = clock();                  // CPU clock                 
  if(SELF_GRAV) {                                  // IF self-gravitating (G!=0)
    report::info("growing tree");
    if(grow_tree) falcON::grow(NCRIT,ROOTCENTER);  //   grow new tree OR        
    else          falcON::reuse();                 //   re-use old tree         
    record_cpu(cpu,CPU_BUILD);                     //   record CPU time required
    report::info("compute gravity");
#ifdef falcON_ADAP
    if(SOFTENING==individual_adaptive) {           //   IF adaptive softening   
      if(adjust) falcON::approximate_gravity(true,all,NSOFT,NREF,EMIN,EFAC);
      else       falcON::approximate_gravity(true,all,NSOFT,NREF,zero,zero);
    } else                                         //   ELIF: fixed softening   
      falcON::approximate_gravity(true,all,zero,0u,EMIN,zero);
#else
    falcON::approximate_gravity(true,all);         //   ELIF: fixed softening   
#endif
    record_cpu(cpu,CPU_GRAV);                      //   record CPU for gravity  
    if(using_extpot()) {                           //   IF external potential   
      report::info("adding external potential");
      PEX->sad(time(),                             //   get external pot, acc   
	       BODIES->N_bodies(),
	       BODIES->mass_s(),
	       BODIES->pos_s(),
	       BODIES->flg_s(),
	       BODIES->pex_s(),
	       BODIES->acc_s(), all);              //     set pex, add acc      
      record_cpu(cpu,CPU_PEX);                     //     record CPU            
    }                                              //   ENDIF                   
  } else if(using_extpot()) {                      // ELIF external potential   
    report::info("setting external potential");
    PEX->set(time(),                               //   get external pot, acc   
	     BODIES->N_bodies(),
	     BODIES->mass_s(),
	     BODIES->pos_s(),
	     BODIES->flg_s(),
	     BODIES->pex_s(),
	     BODIES->acc_s(), all);                //   set pex, set acc        
    record_cpu(cpu,CPU_PEX);                       //   record CPU              
    if(all)                                        //   IF all are active       
      LoopMyBodies Bi.pot()= zero;                 //     reset internal pot    
    else                                           //   ELSE                    
      LoopMyBodies                                 //     LOOP active bodies    
	if(is_active(Bi)) Bi.pot()=zero;           //       reset internal pot  
  } else {                                         // ELSE (no gravity at all)  
    if(all)                                        //   IF all are active       
      LoopMyBodies {                               //     LOOP all bodies       
        Bi.pot() = zero;                           //     reset potential       
        Bi.acc() = zero;                           //     reset acceleration    
      }                                            //   END LOOP                
    else                                           //   ELSE                    
      LoopMyBodies if(is_active(Bi)) {             //     LOOP active bodies    
	Bi.pot() = zero;                           //       reset potential     
	Bi.acc() = zero;                           //       reset acceleration  
      }                                            //     END LOOP              
    record_cpu(cpu,CPU_GRAV);                      //   record CPU for gravity  
  }                                                // ENDIF                     
}
//------------------------------------------------------------------------------
void basic_nbody::reset_softening(
				  kern_type const&ker,
				  real      const&e
#ifdef falcON_ADAP
				 ,real      const&ns,
				  uint      const&nr,
				  real      const&em,
				  real      const&ef
#endif
				  )
{
#ifdef falcON_ADAP
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
inline
basic_nbody::basic_nbody(const bodies *const&b,    // I: bodies                 
			 real          const&e,    // I: eps/eps_max            
			 double        const&ti,   // I: t_initial              
			 real          const&th,   // I: tolerance parameter    
			 int           const&nc,   // I: N_crit                 
			 const vect*   const&x0,   // I: pre-determined center  
			 kern_type     const&ke,   // I: softening kernel       
			 real          const&g,    // I: Newton's G             
#ifdef falcON_INDI
			 soft_type     const&sf,   // I: softening type         
#ifdef falcON_ADAP
			 real          const&ns,   // I: N_soft                 
			 uint          const&nr,   // I: N_ref                  
			 real          const&em,   // I: eps_min                
			 real          const&ef,   // I: eps_fac                
#endif
#endif
			 const gravity*const&px,   // I: P_ex                   
			 const int           gd[4] // I: direct sum: gravity    
#ifdef falcON_SPH
			 ,
			 const int           sd[3] // I: direct sum: SPH        
#endif
			 ) :
  nbody_base( b, e, ti, th, nc, x0, ke, g,
#ifdef falcON_INDI
	      sf != global_fixed,
#endif
	      px, gd
#ifdef falcON_SPH
	      ,sd
#endif
	      )
#ifdef falcON_INDI
  , SOFTENING ( sf )
#endif
{
#ifdef falcON_INDI
  if(SOFTENING==individual_fixed && !BODIES->has(io::e)) 
    falcON_ErrorF("individual fixed softening, but no epsi given",
		  "basic_nbody");
#endif
#ifdef falcON_ADAP
  reset_softening(ke,e,ns,nr,em,ef);
#else
  reset_softening(ke,e);
#endif
}
//------------------------------------------------------------------------------
void basic_nbody::reset_opening(const real th) const {
  falcON::reset_opening(abs(th),
			th<zero? const_theta : theta_of_M);
}
//------------------------------------------------------------------------------
#ifdef falcON_ADAP
inline void basic_nbody::estimate_mass_density(const bool ext)
{
  register clock_t  cpu = clock();                 // cpu time                  
  falcON::estimate_rho(NREF);                      // estimate rho in tree      
  if(PEX && !PEX->is_empty() && ext)               // IF external potential     
    PEX->add_rho(time(),
		 BODIES->N_bodies(),
		 BODIES->mass_s(),
		 BODIES->pos_s(),
		 BODIES->flg_s(),
		 BODIES->rho_s());
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
{
  report REPORT("basic_nbody::diagnose()");
  double m=0., vin=0.;
  vect_d x=0., v=0.;
  amom_d l=0.;
  double k[Ndim][Ndim]={0.}, w[Ndim][Ndim]={0.};
  if(PEX && !PEX->is_empty()) {                    // IF have external pot      
    double vex=0.;
    LoopMyBodies {                                 //   LOOP bodies             
      register double mi = mass(Bi);               //     m                     
      m  += mi;                                    //     add: total mass       
      vin+= mi * pot(Bi);                          //     add: int pot energy   
      vex+= mi * pex(Bi);                          //     add: ext pot energy   
      register vect_d mx = mi * pos(Bi);           //     m * x                 
      register vect_d mv = mi * vel(Bi);           //     m * v                 
      AddTensor(k,mv,vel(Bi));                     //     add to K_ij           
      AddTensor(w,mx,acc(Bi));                     //     add to W_ij           
      x += mx;                                     //     add: dipole           
      v += mv;                                     //     add: total momentum   
      l += mx ^ vect_d(vel(Bi));                   //     add: total ang mom    
    }                                              //   END LOOP                
    Vex = vex;                                     // total ext pot energy      
  } else {                                         // ELSE: no external pot     
    LoopMyBodies {                                 //   LOOP bodies             
      register double mi = mass(Bi);               //     m                     
      m  += mi;                                    //     add: total mass       
      vin+= mi * pot(Bi);                          //     add: int pot energy   
      register vect_d mx = mi * pos(Bi);           //     m * x                 
      register vect_d mv = mi * vel(Bi);           //     m * v                 
      AddTensor(k,mv,vel(Bi));                     //     add to K_ij           
      AddTensor(w,mx,acc(Bi));                     //     add to W_ij           
      x += mx;                                     //     add: dipole           
      v += mv;                                     //     add: total momentum   
      l += mx ^ vect_d(vel(Bi));                   //     add: total ang mom    
    }                                              //   END LOOP                
    Vex = zero;                                    // total ext pot energy      
  }                                                // ENDIF                     
  M    = m;                                        // total mass                
  T    = half*tr(k);                               // total kin energy          
  Vin  = half*vin;                                 // total int pot energy      
  W    = tr(w);                                    // total pot energy from acc 
  TW   =-T/W;                                      // virial ratio              
  L    = l;                                        // total angular momentum    
  m    = 1./m;
  CMX  = x*m;                                      // center of mass            
  CMV  = v*m;                                      // center of mass velocity   
  for(register int i=0; i!=Ndim; ++i)
    for(register int j=0; j!=Ndim; ++j) {
      KT[i][j] = half * k[i][j];                   // kin energy tensor         
      WT[i][j] = half *(w[i][j]+w[j][i]);          // pot energy tensor         
    }
}
//------------------------------------------------------------------------------
inline void basic_nbody::stats_front(ostream& to) const
{
  ios::fmtflags old = to.flags();
  to.setf(ios::left | ios::showpoint);
  update_diags();
  to<<setprecision(5)<<setw(10)<<time()
    <<' '
    <<setprecision(7)<<setw(13)<<total_energy()
    <<' '
    <<setprecision(5)<<setw(10)<<kin_energy()
    <<' ';
  if(SELF_GRAV)
    to<<setprecision(5)<<setw(11)<<pot_self_energy()
      <<' ';
  if(PEX && !PEX->is_empty())
    to<<setprecision(5)<<setw(11)<<pot_ext_energy()
      <<' ';
  to<<setprecision(5)<<setw(11)<<pot_energy_acc()  
    <<' '
    <<setprecision(min(5,max(1,5-int(-log10(twice(virial_ratio()))))))
    <<setw(7)<<twice(virial_ratio())
    <<' '
    <<setprecision(5)<<setw(10)
#if falcON_NDIM==3
    <<sqrt(norm(total_angmom()))
#else
    <<total_angmom()
#endif
    <<' '
    <<setprecision(2)<<setw(8)<<sqrt(norm(total_momentum()));
  to.flags(old);
}
//------------------------------------------------------------------------------
inline void basic_nbody::stats_head_front(ostream& to) const
{
  to<< "    time   "
    << "    E=T+V     "
    << "   T       ";
  if(SELF_GRAV)               to<<"   V_in     ";
  if(PEX && !PEX->is_empty()) to<<"   V_ex     ";
  to<< "   W        "
    << " -2T/W  "
#if falcON_NDIM==3
    << "   |L|     "
#else
    << "   L_z     "
#endif
    << " |v_cm| ";
}
//------------------------------------------------------------------------------
inline void basic_nbody::stats_line_front(ostream& to) const
{
  to<<" ---------------------------------------"
    <<"------------------------------------";
  if(SELF_GRAV)               to<<"------------";
  if(PEX && !PEX->is_empty()) to<<"------------";
}
//------------------------------------------------------------------------------
namespace {
  inline void print5(double const&x, ostream&to) 
  {
    if(x < 100)
      to<<setw(2)<<setfill(' ')<<int(x)<<'.'
	<<setw(2)<<setfill('0')<<int(100*(x-int(x)));
    else if(x<1000)
      to<<setw(3)<<setfill(' ')<<int(x)<<'.'
	<<setw(1)<<setfill('0')<<int(10*(x-int(x)));
    else
      to<<setw(5)<<setfill(' ')<<int(x+0.5);
  }
}
//------------------------------------------------------------------------------
inline void basic_nbody::stats_back(ostream& to) const
{
  int    h,m,s,c;
  double t = cpu_total();
  h = int(t/3600); t-= 3600*h;
  m = int(t/60);   t-= 60*m;
  s = int(t);      t-= s;
  c = int(100*t);
  print5(cpu_build(),   to); to<<' ';
  print5(cpu_grav(),    to); to<<' ';
  print5(cpu_longstep(),to); to<<' ';
  to<<setw(3)<<setfill(' ')<<h<<':'
    <<setw(2)<<setfill('0')<<m<<':'
    <<setw(2)<<s<<'.'<<setw(2)<<c
    <<setfill(' ');
}
//------------------------------------------------------------------------------
inline void basic_nbody::stats_head_back(ostream& to) const
{
//to<< "ss.cc ss.cc ss.cc hhh:mm:ss.cc";
  to<< " tree  grav  step  accumulated";
}
//------------------------------------------------------------------------------
inline void basic_nbody::stats_line_back(ostream& to) const
{
  to<<"------------------------------";
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::LeapFrogCode                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
LeapFrogCode::LeapFrogCode(const bodies *const&b,  // I: bodies                 
			   real          const&e,  // I: eps/eps_max            
			   double        const&ti, // I: t_initial              
			   int           const&h0, // I: h0                     
			   int           const&hg, // I: h_grow                 
			   real          const&th, // I: tolerance parameter    
			   int           const&nc, // I: N_crit                 
			   const vect*   const&x0, // I: root center            
			   kern_type     const&ker,// I: softening kernel       
			   real          const&g,  // I: Newton's G             
#ifdef falcON_INDI
			   soft_type     const&sf, // I: softening type         
#ifdef falcON_ADAP
			   real          const&ns, // I: N_soft                 
			   uint          const&nr, // I: N_ref                  
			   real          const&em, // I: eps_min                
			   real          const&ef, // I: eps_fac                
#endif
#endif
			   const gravity*const&p,  // I: P_ex                   
			   const int        nd[4]):// I: direct sum control     
  basic_nbody       ( b,e,ti,th,nc,x0,ker,g,
#ifdef falcON_INDI
		      sf,
#ifdef falcON_ADAP
		      ns,nr,em,ef,
#endif
#endif
		      p,nd ),
  LeapFrog<bodies> ( h0, ti, (1<<hg)-1 )
{
  register clock_t  cpu = clock();                 // cpu time                  
  reset_cpus();                                    // reset CPU time counters   
  LoopMyBodies Bi.flag_as_active();                // flag everybody as active  
  set_gravity(true,true                            // set eps & get acc         
#ifdef falcON_ADAP
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
  if(REUSED < REUSE) {                             // IF may re-use old tree    
    set_gravity(false, true);                      //   adjust eps & get acc    
    REUSED++;                                      //   increment # re-using    
  } else {                                         // ELSE                      
    set_gravity(true, true);                       //   adjust eps & get acc    
    REUSED=0;                                      //   reset # re-using        
  }                                                // ENDIF                     
  accelerate(BODIES);                              // accelerate = kick         
  DIAG = false;                                    // diagnose out-of-date      
  record_cpu(cpu,CPU_STEP);                        // record CPU time           
  update_cpu_total();                              // record accumulated time   
}
//------------------------------------------------------------------------------
void LeapFrogCode::stats(ostream& to) const
{
  stats_front(to); to << ' ';
  stats_back (to); to << endl;
}
//------------------------------------------------------------------------------
void LeapFrogCode::stats_head(ostream& to) const
{
  stats_head_front(to);
  to<<' ';
  stats_head_back(to);
  to<<endl;
}
//------------------------------------------------------------------------------
void LeapFrogCode::stats_line(ostream &to) const
{
  stats_line_front(to);
  stats_line_back (to);
  to << endl;
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
  GravBlockStep<bodies>::reset_scheme(fa,fp,fc,fe,using_extpot());
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
  static uint m = 0u;                              // # of tiny moves omitted   
  clock_on();                                      // increase actual clock     
  if(need_move(low)) {                             // IF anybody is active now  
    register real dt = tau_min() * (m+1);          //   dt = (m+1)*tau_min      
    ForAll move_by(Bi,dt);                         //   all:    x -> x+v*dt     
    m = 0;                                         //   reset m = 0             
    set_active_flags(low);                         //   IF(l>=low) -> active    
    set_gravity(low<LGROW,low==0);                 //   active: a = acc(x)      
    ForActive acce_half(Bi);                       //   active: v -> v+a*h/2    
    if(low != highest_level())                     //   IF(levels may change)   
#ifdef falcON_INDI
      if(SOFTENING != global_fixed) {              //     IF individual eps     
	ForActive adjust_level(Bi,low,::eps(Bi));  //       active: h -> h'     
      } else                                       //     ELSE: individual eps_i
#endif
	ForActive adjust_level(Bi,low,eps());      //       active: h -> h'     
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
  GravBlockStep<bodies>::reset_steps(h0,nl,TINI);  // set time steps            
  LoopMyBodies Bi.flag_as_active();                // flag every body as active 
  set_gravity(true,                                // set eps & get acc         
	      true                                 // for all leafs             
#ifdef falcON_ADAP
	     ,false                                // set, don't adjust eps_i   
#endif
	      );
  diagnose();                                      // energy et al. at t=t0     
  DIAG = true;                                     // diagnose up-to-date       
  GravBlockStep<bodies>::reset_counts();           // reset counts              
  LoopMyBodies                                     // loop bodies               
    assign_level(Bi);                              //   assign step level       
  record_cpu(cpu,CPU_STEP);                        // record total CPU          
  update_cpu_total();                              // record accumulated time   
}
//------------------------------------------------------------------------------
BlockStepCode::
BlockStepCode(const bodies *const&b,               // I: bodies                 
	      real          const&e,               // I: eps/eps_max            
	      double        const&ti,              // I: t_initial              
	      int           const&h0,              // I: h0                     
	      int           const&nl,              // I: # levels               
	      real          const&fa,              // I: f_a: for stepping      
	      real          const&fp,              // I: f_p: for stepping      
	      real          const&fc,              // I: f_c: for stepping      
	      real          const&fe,              // I: f_e: for stepping      
	      int           const&hg,              // I: h_grow                 
	      real          const&th,              // I: tolerance parameter    
	      int           const&nc,              // I: N_crit                 
	      const vect*   const&x0,              // I: root center            
	      kern_type     const&ker,             // I: softening kernel       
	      real          const&g,               // I: Newton's G             
#ifdef falcON_INDI
	      soft_type     const&sf,              // I: softening type         
#ifdef falcON_ADAP
	      real          const&ns,              // I: N_soft                 
	      uint          const&nr,              // I: N_ref                  
	      real          const&em,              // I: eps_min                
	      real          const&ef,              // I: eps_fac                
#endif
#endif
	      const gravity*const&p,               // I: P_ex                   
	      const int        nd[4]):             // I: direct sum control     
  basic_nbody            (b,e,ti,th,nc,x0,ker,g,
#ifdef falcON_INDI
			  sf,
#ifdef falcON_ADAP
			  ns,nr,em,ef,
#endif
#endif
			  p,nd),
  GravBlockStep<bodies>  (),
  LGROW                  ((nl <= hg)? 1 : nl-hg),
  W                      (max(5, 1+int(log10(double(b->N_bodies())))))
{
  GravBlockStep<bodies>::reset_scheme(fa,fp,fc,fe,using_extpot());
  if(!BODIES->has(io::l))
    falcON_ErrorF("bodies::has(io::l) == false","BlockStepCode");
  prepare(h0,nl);
}
//------------------------------------------------------------------------------
void BlockStepCode::stats(ostream& to) const
{
  stats_front(to); to<<' ';
  if(highest_level()) {
    ios::fmtflags old = to.setf(ios::right);
    GravBlockStep<bodies>::short_stats(to,W);
    to.flags(old);
  }
  stats_back (to); to << endl;
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
  stats_head_front(to);
  to<<' ';
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
  stats_head_back(to);
  to<<endl;
}
//------------------------------------------------------------------------------
void BlockStepCode::stats_line(ostream &to) const
{
  stats_line_front(to);
  if(highest_level())
    for(register int i=0; i<=highest_level(); i++) put_char(to,'-',W+1);
  stats_line_back (to);
  to << endl;
}
//------------------------------------------------------------------------------
#undef LoopMyBodies
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::NbodyCode                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
NbodyCode::NbodyCode(const char         *file,     //   I: input file           
		     bool          const&resume,   //   I: resume old (if nemo) 
		     int           const&hmin,     //   I: hmin: t_min=2^(-hmin)
		     int           const&Nlev,     //   I: # time-step levels   
		     real          const&fac,      //   I: f_a                  
		     real          const&fph,      //   I: f_p                  
		     real          const&fap,      //   I: f_c                  
		     real          const&fep,      //   I: f_e                  
		     int           const&Ncrit,    //   I: N_crit               
		     int           const&hgrow,    //   I: h_grow               
		     const vect*   const&proot,    //   I: pre-set root center  
		     real          const&eps,      //   I: eps                  
		     kern_type     const&kernel,   //   I: softening kernel     
		     const gravity*const&pex,      //   I: P_ex                 
		     real          const&theta,    //   I: tolerance parameter  
		     real          const&Grav,     //   I: Newton's G           
#ifdef falcON_INDI
#ifdef falcON_ADAP
		     real          const&Nsoft,    //  [I: Nsoft]               
		     uint          const&Nref,     //  [I: Nref]                
		     real          const&emin,     //  [I: emin]                
#endif
		     basic_nbody::soft_type
		                   const&soften,   //  [I: softening type]      
#endif
		     int           const&Nout,     //  [I: # nemo output devs]  
		     io            const&read_more,//  [I: what else to read?]  
		     const         int direct[4]) ://  [I: direct sum: gravity] 
  EFAC   ( 1.1 ),
  ND     ( Nout ),
  SUPPORT( io::gravity                            |
	   (Nlev   > 1           ? io::l : io::o) |
	   (pex   != 0           ? io::q : io::o) |
	   (soften!=basic_nbody::global_fixed ? io::e : io::o) ),
  INPUT  ( io::mxv | read_more |
	   (soften==basic_nbody::individual_fixed? io::e : io::o) ),
  NEMO   ( ND>0? new nemo_out[ND] : 0 ),
  BODIES ( new bodies(0, SUPPORT) ),
  IFILE  ( file ),
  NBODY  ( 0 )
{
  nemo_in NemoIn(IFILE.c_str());                   // open nemo input           
  io      got;                                     // what did we get?          
  double  tini, cini;                              // initial time, CPU time    
  if(resume) {                                     // IF (resuming old sim)     
    do {                                           //   DO: read particles      
      NemoIn.open_set(nemo_io::snap);              //     open nemo snapshot    
      BODIES->read_nemo_particles(NemoIn,got,&tini,INPUT,0,0); // read bodies   
      if(NemoIn.is_present(nemo_io::diags)) {      //     IF(diags set)         
	NemoIn.open_set(nemo_io::diags);           //       open diags set      
	if(NemoIn.is_present(nemo_io::cputime))    //       IF(cpu time there)  
	  cini =NemoIn.read(nemo_io::cputime)*60;  //         read CPU [min]    
	NemoIn.close_set(nemo_io::diags);          //       close diags set     
      }                                            //     ENDIF                 
      NemoIn.close_set(nemo_io::snap);             //     close nemo snapshot   
    } while(NemoIn.is_present(nemo_io::snap));     //   WHILE more to be read   
  } else {                                         // ELSE (not resuming)       
    NemoIn.open_set(nemo_io::snap);                //   open nemo snapshot      
    BODIES->read_nemo_particles(NemoIn,got,&tini,INPUT,0,0);  // read bodies    
    NemoIn.close_set(nemo_io::snap);               //   close nemo snapshot     
  }                                                // ENDIF                     
  if(!got.contains(INPUT)) {                       // IF some data missing      
    char w[20];                                    //   string for missing data 
    got.missing(INPUT).make_word(w);               //   which data are missing  
    error("couldn't read body data: %s",w);        //   issue error: data amiss 
  }                                                // ENDIF                     
#ifdef falcON_ADAP
  if(soften == individual_adaptive && hgrow)
    error("inidividual adaptive softening lengths must not be used"
	  "with hgrow != 0");
#endif
  NBODY = Nlev > 1 ?
    static_cast<basic_nbody*>( new BlockStepCode(BODIES,eps,tini,hmin+1-Nlev,
						 Nlev,fac,fph,fap,fep,
						 hgrow,theta,Ncrit,
						 proot,kernel,Grav,
#ifdef falcON_INDI
						 soften,
#ifdef falcON_ADAP
						 Nsoft,Nref,emin,EFAC,
#endif
#endif
						 pex,direct) ) :
    static_cast<basic_nbody*>( new LeapFrogCode (BODIES,eps,tini,hmin,
						 hgrow,theta,Ncrit,
						 proot,kernel,Grav,
#ifdef falcON_INDI
						 soften,
#ifdef falcON_ADAP
						 Nsoft,Nref,emin,EFAC,
#endif
#endif
						 pex,direct) );
  if(resume) NBODY->reset_cpu_total(cini);
}
//------------------------------------------------------------------------------
void NbodyCode::describe_nemo(std::ostream&out,    // I: output stream          
			      const char  *com)    // I: command line           
{
  if(!okay())
    error("something wrong with N-body data, perhaps not read yet?");
  out<<"#"; NBODY->stats_line(out);
  out<<"# \""<<com<<"\"\n#\n";
  out<<"# run at  "  <<run_info::time()<<"\n";
  if(run_info::user_known()) out<<"#     by  \""<<run_info::user()<<"\"\n";
  if(run_info::host_known()) out<<"#     on  \""<<run_info::host()<<"\"\n";
  if(run_info::pid_known())  out<<"#     pid  " <<run_info::pid() <<"\n";
  out<<"#\n";
  out.flush();
}
//------------------------------------------------------------------------------
void NbodyCode::open_nemo(int  const&d,            //[I: index of nemo stream]  
			  const char*file,         //[I: file name ]            
			  bool const&resume)       //[I" resume old sim?]       
{
  if(d >= ND) error("NbodyCode::open_nemo(): nemo device %d does not exist",d);
  if(resume) {
    if(file==0 || IFILE == std::string(file))
      (NEMO+d)->open_to_append(IFILE.c_str());
    else
      (NEMO+d)->open(file);
  } else if(file) {
    if(IFILE == std::string(file) && IFILE != "-")
      error("out==in; use option resume instead");
    (NEMO+d)->open(file);
  } else
    warning("cannot open nemo output");
}
//------------------------------------------------------------------------------
void NbodyCode::write_nemo(io   const&w,           //[I: what to write out]     
			   int  const&d,           //[I: index of nemo stream]  
			   bool const&diag)        //[I: output diagnostics?]   
{
  if(d >= ND)
    error("NbodyCode::write_nemo(): nemo device %d does not exist\n",d);
  nemo_out*OUT=NEMO+d;
  if(!OUT->is_open())
    error("NbodyCode::write_nemo(): nemo device %d not open\n",d);
  register uint        i,j;
  OUT->open_set(nemo_io::snap);                    // OPEN a new nemo snapshot  
  double t_now = time();
  BODIES->write_nemo_particles(*OUT, &t_now,       //   write out particles     
#ifdef falcON_INDI
			       (NBODY->use_individual_eps() ? io::e : io::o) |
#endif
			       w);
  if(diag) {                                       //   IF diags output wantd   
    OUT->open_set(nemo_io::diags);                 //     open diagnostics set  
    OUT->single_vec(0)=NBODY->total_energy();      //       copy total energy   
    OUT->single_vec(1)=NBODY->kin_energy();        //       copy kin energy     
    OUT->single_vec(2)=NBODY->pot_energy();        //       copy pot energy     
    OUT->write(nemo_io::energy);                   //       write energies      
    for(i=0;i!=Ndim;++i) for(j=0;j!=Ndim;++j)      //     LOOP dims twice       
      OUT->single_mat(i,j) = NBODY->kin_energy(i,j);//      copy K_ij           
    OUT->write(nemo_io::KinT);                     //     write K_ij            
    for(i=0;i!=Ndim;++i) for(j=0;j!=Ndim;++j)      //     LOOP dims twice       
      OUT->single_mat(i,j) = NBODY->pot_energy(i,j);//      copy W_ij           
    OUT->write(nemo_io::PotT);                     //     write W_ij            
    for(i=0;i!=Ndim;++i) for(j=0;j!=Ndim;++j)      //     LOOP dimensions       
      OUT->single_mat(i,j) = NBODY->as_angmom(i,j); //      copy A_ij           
    OUT->write(nemo_io::AmT);                      //     write A_ij            
    for(i=0;i!=Ndim;++i) {                         //     LOOP dims             
      OUT->single_phs(0,i)=NBODY->center_of_mass()[i];//    copy c-of-m pos     
      OUT->single_phs(1,i)=NBODY->total_momentum()[i];//    copy c-of-m vel     
    }                                              //     END LOOP              
    OUT->write(nemo_io::cofm);                     //     write c-of-m (x,v)    
    OUT->write(nemo_io::cputime,                   //     write accum CPU [min] 
	       NBODY->cpu_total()/60.);            //                           
    OUT->close_set(nemo_io::diags);                //     close diagnostics set 
  }                                                //   ENDIF                   
  OUT->close_set(nemo_io::snap);                   // CLOSE nemo snapshot       
  OUT->reset();                                    // reset nemo output         
}
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_NEMO               
  

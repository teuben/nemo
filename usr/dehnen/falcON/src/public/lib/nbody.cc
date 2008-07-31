//-----------------------------------------------------------------------------+
//                                                                             |
// nbody.cc                                                                    |
//                                                                             |
// Copyright (C) 2000-2008  Walter Dehnen                                      |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/nbody.h>
#ifdef falcON_MPI
#  include <parallel/snapshot.h>
#endif
#include <iomanip>
using namespace falcON;
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::Integrator                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
Integrator::Integrator(const ForceAndDiagnose*S,
		       fieldset p, fieldset k, fieldset r,
		       fieldset P, fieldset K, fieldset R) falcON_THROWING :
  // all quantities remembered must also be predicted
  predALL     ((S->requires()&fieldset::w ? p | fieldset::w : p) | r),
  kickALL     ( S->requires()&fieldset::w ? k | fieldset::v : k ),
  rembALL     ( S->requires()&fieldset::w ? r | fieldset::w : r ),
  // don't do to SPH what has already been done for all
  predSPH     ((P|R) &~ predALL ),
  kickSPH     ( K    &~ kickALL ),
  rembSPH     ( R    &~ rembALL ),
  SOLVER      ( S ),
  C_OLD       ( clock() ),
  CPU_TOTAL   ( 0. )
{
  // sanity check prediction, kick, and rembembrance settings
  char     comp[32];
  fieldset test;
  if( (test=predALL & ~fieldset(fieldset::x|fieldset::w)) ) 
    falcON_Warning("Integration: will not predict '%s'", test.make_word(comp));
  if( (test=kickALL & ~fieldset(fieldset::v)) )
    falcON_Warning("Integration: will not kick '%s'", test.make_word(comp));
  if( (test=rembALL & ~fieldset(fieldset::w)) )
    falcON_Warning("Integration: will not remember '%s'", test.make_word(comp));
  if( predALL & fieldset::w && !(kickALL & fieldset::v) )
    falcON_THROW("Integration: cannot predict w without kicking v");
  if( predALL & fieldset::x && !(kickALL & fieldset::v) )
    falcON_THROW("Integration: request to predict x without kicking v");
  requALL = 
    kickALL & fieldset::v ? fieldset::a : fieldset::empty |
    predALL & fieldset::w ? fieldset::a : fieldset::empty;
#ifdef falcON_SPH
  if( (test=predSPH & ~fieldset(fieldset::x|fieldset::H|fieldset::R|
				fieldset::Y|fieldset::V)) ) 
    falcON_Warning("Integration: will not predict '%s'", test.make_word(comp));
  if( (test=kickSPH & ~fieldset(fieldset::U|fieldset::v)) )
    falcON_Warning("Integration: will not kick '%s'", test.make_word(comp));
  if( (test=rembSPH & ~fieldset(fieldset::V|fieldset::Y)) )
    falcON_Warning("Integration: will not remember '%s'", test.make_word(comp));
  if( predSPH & fieldset::V && !((kickSPH|kickALL) & fieldset::v) )
    falcON_THROW("Integration: cannot predict V without kicking v");
  if( predSPH & fieldset::Y && !(kickSPH & fieldset::U) )
    falcON_THROW("Integration: cannot predict Y without kicking U");
  if( predSPH & fieldset::x && !((kickSPH|kickALL) & fieldset::v) )
    falcON_THROW("Integration: request to predict x without kicking v");
  requSPH =
//     predSPH & fieldset::H ? fieldset::J : fieldset::empty |
    predSPH & fieldset::V ? fieldset::a : fieldset::empty |
//     predSPH & fieldset::R ? fieldset::D : fieldset::empty |
    kickSPH & fieldset::U ? fieldset::I : fieldset::empty;
#endif
  // sanity check requirements and delivieries
  reset_CPU();
  if(! SOLVER->computes().contain(required()))
    falcON_THROW("Integrator requires '%s', but ForceSolver computes '%s'",
		 word(required()),word(S->computes()));
  fieldset full = SOLVER->computes() | predicted() | fieldset::m |
                  rembALL | kickALL;
  if(! full.contain(SOLVER->requires()) )
    falcON_THROW("ForceAndDiagnose requires '%s', but code delivers only '%s'",
		 word(S->requires()), word(full));
  if(! SOLVER->computesSPH().contain(requiredSPH()))
    falcON_THROW("SPH: Integrator requires '%s', but ForceSolver computes '%s'",
		 word(requiredSPH()), word(S->computesSPH()));
  full |= SOLVER->computesSPH() | predictedSPH() | rembSPH | kickSPH;
  if(! full.contain(SOLVER->requiresSPH()) )
    falcON_THROW("SPH: ForceAndDiagnose requires '%s', "
		 "but code delivers only '%s'",
		 word(S->requiresSPH()), word(full));
  // make sure snapshot supports all data required
  fieldset need = p | k | r | P | K | R | fieldset::f;
  need |= SOLVER->computes() | SOLVER->computesSPH();
  snap_shot()->add_fields(need);
}
////////////////////////////////////////////////////////////////////////////////
namespace {
  //----------------------------------------------------------------------------
  template<int TYPE, int DERIV>
  inline void move_all(const bodies*B, double dt, bool all) {
    if(all)
      LoopAllBodies(B,b) {
	b. template datum<TYPE>() += dt * const_datum<DERIV>(b);
      }
    else
      LoopAllBodies(B,b) if(is_active(b)) {
	b. template datum<TYPE>() += dt * const_datum<DERIV>(b);
      }
  }
  //----------------------------------------------------------------------------
  template<int TYPE, int DERIV>
  inline void move_all_i(const bodies*B, const double*dt, bool all) {
    if(all)
      LoopAllBodies(B,b) {
	b. template datum<TYPE>() += dt[level(b)] * const_datum<DERIV>(b);
      }
    else
      LoopAllBodies(B,b) if(is_active(b)) {
	b. template datum<TYPE>() += dt[level(b)] * const_datum<DERIV>(b);
      }
  }
  //----------------------------------------------------------------------------
  template<int TYPE, int DERIV>
  inline void move_sph(const bodies*B, double dt, bool all) {
    if(all)
      LoopSPHBodies(B,b)
	b. template datum<TYPE>() += dt * const_datum<DERIV>(b);
    else
      LoopSPHBodies(B,b) if(is_active(b))
	b. template datum<TYPE>() += dt * const_datum<DERIV>(b);
  }
  //----------------------------------------------------------------------------
  template<int TYPE, int DERIV>
  inline void move_sph_i(const bodies*B, const double*dt, bool all) {
    if(all)
      LoopSPHBodies(B,b)
	b. template datum<TYPE>() += dt[level(b)] * const_datum<DERIV>(b);
    else
      LoopSPHBodies(B,b) if(is_active(b))
	b. template datum<TYPE>() += dt[level(b)] * const_datum<DERIV>(b);
  }
  //----------------------------------------------------------------------------
  template<int TYPE, int COPY>
  inline void copy_all(const bodies*B, bool all) {
    if(all)
      LoopAllBodies(B,b)
	b. template datum<TYPE>() = const_datum<COPY>(b);
    else
      LoopAllBodies(B,b) if(is_active(b))
	b. template datum<TYPE>() = const_datum<COPY>(b);
  }
  //----------------------------------------------------------------------------
  template<int TYPE, int COPY>
  inline void copy_sph(const bodies*B, bool all) {
    if(all)
      LoopSPHBodies(B,b)
	b. template datum<TYPE>() = const_datum<COPY>(b);
    else
      LoopSPHBodies(B,b) if(is_active(b))
	b. template datum<TYPE>() = const_datum<COPY>(b);
  }
}
//------------------------------------------------------------------------------
void Integrator::drift(double dt, bool all) const
{
  const snapshot*const&B(SOLVER->snap_shot());
  B->advance_time_by(dt);
  if(predALL & fieldset::x) move_all<fieldbit::x,fieldbit::v>(B,dt,all);
  if(predALL & fieldset::w) move_all<fieldbit::w,fieldbit::a>(B,dt,all);
#ifdef falcON_SPH
  if(predSPH & fieldset::x) move_sph<fieldbit::x,fieldbit::v>(B,dt,all);
  if(predSPH & fieldset::V) move_sph<fieldbit::V,fieldbit::a>(B,dt,all);
  if(predSPH & fieldset::Y) move_sph<fieldbit::Y,fieldbit::I>(B,dt,all);
#endif
}
//------------------------------------------------------------------------------
void Integrator::kick(double dt, bool all) const
{
  const snapshot*const&B(SOLVER->snap_shot());
  if(kickALL & fieldset::v) move_all<fieldbit::v,fieldbit::a>(B,dt,all);
#ifdef falcON_SPH
  if(kickSPH & fieldset::v) move_sph<fieldbit::v,fieldbit::a>(B,dt,all);
  if(kickSPH & fieldset::U) move_sph<fieldbit::U,fieldbit::I>(B,dt,all);
#endif
}
//------------------------------------------------------------------------------
void Integrator::kick_i(const double*dt, bool all) const
{
  const snapshot*const&B(SOLVER->snap_shot());
  if(kickALL & fieldset::v) move_all_i<fieldbit::v,fieldbit::a>(B,dt,all);
#ifdef falcON_SPH
  if(kickSPH & fieldset::v) move_sph_i<fieldbit::v,fieldbit::a>(B,dt,all);
  if(kickSPH & fieldset::U) move_sph_i<fieldbit::U,fieldbit::I>(B,dt,all);
#endif
}
//------------------------------------------------------------------------------
void Integrator::remember(bool all) const
{
  const snapshot*const&B(SOLVER->snap_shot());
  if(rembALL & fieldset::w) copy_all<fieldbit::w,fieldbit::v>(B,all);
#ifdef falcON_SPH
  if(rembSPH & fieldset::V) copy_sph<fieldbit::V,fieldbit::v>(B,all);
  if(rembSPH & fieldset::Y) copy_sph<fieldbit::Y,fieldbit::U>(B,all);
#endif
}
//------------------------------------------------------------------------------
void Integrator::cpu_stats_body(output&to) const
{
  SOLVER->cpu_stats_body(to);
#ifdef falcON_MPI
  // need to add code to sum cpu timings over all processes
#endif
  if(to) {
    print_cpu(CPU_STEP,to); to<<' ';
    print_cpu_hms(CPU_TOTAL,to);
  } 
}
//------------------------------------------------------------------------------
void Integrator::describe(output&out)        // I: output stream          
const {
  out<<"#"; stats_line(out);
  if(RunInfo::cmd_known())
    out<<"# \""<<RunInfo::cmd()<<"\"\n#\n";
  out<<"# run at  "  <<RunInfo::time()<<"\n";
  if(RunInfo::user_known()) out<<"#     by  \""<<RunInfo::user()<<"\"\n";
  if(RunInfo::host_known()) out<<"#     on  \""<<RunInfo::host()<<"\"\n";
  if(RunInfo::pid_known())  out<<"#     pid  " <<RunInfo::pid() <<"\n";
  out<<"#\n";
  out.flush();
}
//------------------------------------------------------------------------------
#ifdef falcON_NEMO
//------------------------------------------------------------------------------
void Integrator::write(nemo_out const&o,           // I: nemo output            
		       fieldset       w) const     //[I: what to write]         
{
  if( o.is_sink()) return;
  if(!o.is_open()) 
    falcON_THROW("Integrator::write(): nemo device not open\n");
  snap_shot()->write_nemo(o,w);
} 
#endif
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::LeapFrogCode                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
LeapFrogCode::LeapFrogCode(int kstep, const ForceAndDiagnose *F,
			   fieldset p, fieldset k, fieldset r,
			   fieldset P, fieldset K, fieldset R) falcON_THROWING :
  Integrator(F,p,k,r,P,K,R),
  bodies::TimeSteps(kstep,1)
{
  snap_shot()->set_steps(this);                    // set time steps in bodies  
  remember();                                      // eg: w = v                 
  set_time_derivs(1,1,0.);                         // eg: a = F(x,w)            
  finish_diagnose();                               // finish diagnosis          
  add_to_cpu_step();                               // record CPU time           
  DebugInfo(4,"LeapFrogCode constructed\n");
}
//------------------------------------------------------------------------------
void LeapFrogCode::account_new() const {
  if(snap_shot()->N_new()) {
    LoopAllBodies(snap_shot(),b) 
      if(is_new(b)) b.flag_as_active();
      else          b.unflag_active ();
    set_time_derivs(0,0,0.);
    LoopAllBodies(snap_shot(),b) 
      if(is_new(b)) b.unflag_new();
    snap_shot()->reset_Nnew();
  }
}
//------------------------------------------------------------------------------
void LeapFrogCode::fullstep() const {
  reset_CPU();                                     // reset cpu timers          
  account_new();                                   // account for new bodies    
  kick(tauh(0));                                   // eg: v+= a*tau/2           
  drift(tau(0));                                   // eg: x+= v*tau;  w+= a*tau 
  set_time_derivs(1,1,tau(0));                     // eg: a = F(x,w)            
  kick(tauh(0));                                   // eg: v+= a*tau/2           
  remember();                                      // eg: w = v                 
  finish_diagnose();                               // finish diagnosis          
  add_to_cpu_step();                               // record CPU time           
}    
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::BlockStepCode                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void BlockStepCode::elementary_step(int t) const { // I: number of step         
  static int m = 0;                                // # of tiny moves omitted   
  ++m;                                             // count shortest steps      
  ++t;                                             // add one to t              
  int l=highest_level();                           // find lowest level moving  
  for(; !(t&1) && l; t>>=1, --l);                  // l: lowest level moving    
  bool move=false;                                 // need to do anything?      
  for(int i=l; i!=Nsteps(); ++i)                   // LOOP levels up to highest 
    if(N[i]) move = true;                          //   IF any non-empty: move  
  if(!move) return;                                // IF not moving: DONE       
  bool all=true;                                   // are all active?           
  for(int i=0; i!=l; ++i)                          // LOOP lower levels         
    if(N[i]) all  = false;                         //   IF all empty: all active
  double dt=tau_min() * m;                         // dt = m*tau_min            
  drift(dt);                                       // predict @ new time        
  m = 0;                                           // reset m = 0               
  LoopAllBodies(snap_shot(),b)                     // LOOP bodies               
    if(level(b) >= l) b.flag_as_active();          //   IF(level>=l): active    
    else              b.unflag_active ();          //   ELSE        : inactive  
  set_time_derivs(all,l==0,dt);                    // set accelerations etc     
  kick_i(tauh(),all);                              // kick velocity etc         
  if(l != highest_level() ||                       // IF(levels may change OR   
     ST->always_adjust() )                         //    always adjusting       
    adjust_levels(l, l==0);                        //   h -> h'                 
  if(l) {                                          // IF not last step          
    remember(all);                                 //   remember to be predicted
    kick_i(tauh(),all);                            //   kick by half a step     
  }                                                // ENDIF                     
}
//------------------------------------------------------------------------------
inline void BlockStepCode::account_del() const {
  if(snap_shot()->N_del()) {
    for(int l=0; l!=Nsteps(); ++l)
      N[l] = 0u;
    LoopAllBodies(snap_shot(),b)
      if(!is_new(b)) ++(N[level(b)]);
    snap_shot()->reset_Ndel();
  }
}
//------------------------------------------------------------------------------
inline void BlockStepCode::account_new() const {
  if(snap_shot()->N_new()) {
    LoopAllBodies(snap_shot(),b) 
      if(is_new(b)) b.flag_as_active();
      else          b.unflag_active ();
    set_time_derivs(0,0,0.);
    LoopAllBodies(snap_shot(),b) if(is_new(b)) {
      b.unflag_new();
      ST->assign_level(b, N, highest_level());
    }
    snap_shot()->reset_Nnew();
  }
}
//------------------------------------------------------------------------------
void BlockStepCode::assign_levels() const {
  if(!snap_shot()->have_steps())
    falcON_Error("BlockStepCode::assign_levels(): steps not set\n");
  LoopAllBodies(snap_shot(),b)
    ST->assign_level(b, N, highest_level());
}
//------------------------------------------------------------------------------
void BlockStepCode::adjust_levels(int low, bool all) const {
  if(all)
    LoopAllBodies(snap_shot(),b)
      ST->adjust_level(b, N, low, highest_level());
  else
    LoopAllBodies(snap_shot(),b) if(is_active(b))
      ST->adjust_level(b, N, low, highest_level());
}
//------------------------------------------------------------------------------
void BlockStepCode::update_Nlev(const bodies*B) {
  for(int l=0; l!=Nsteps(); ++l) N[l] = 0;
  LoopAllBodies(B,b)
    ++(N[level(b)]);
}
//------------------------------------------------------------------------------
void BlockStepCode::fullstep() const {
  reset_CPU();                                     // reset cpu timers          
  account_new();                                   // account for new bodies    
  account_del();                                   // account for removed bodies
  remember(true);                                  // remember to be predicted  
  kick_i(tauh(),true);                             // kick by half a step       
  for(int t=0; t!=1<<highest_level(); ++t)         // LOOP elementary steps     
    elementary_step(t);                            //   elementary step         
  finish_diagnose();                               // finish diagnosis          
  add_to_cpu_step();                               // record CPU time           
}
//------------------------------------------------------------------------------
BlockStepCode::BlockStepCode(int      km,          // I: tau_max = 2^-kmax      
			     unsigned Ns,          // I: #steps                 
			     const ForceAndDiagnose *F,
			     const StepLevels       *S,
			     fieldset p, fieldset k, fieldset r,
			     fieldset P, fieldset K, fieldset R, int w)
  falcON_THROWING :
  Integrator        ( F,p,k,r,P,K,R ),
  bodies::TimeSteps ( km, Ns),
  N                 ( Ns? falcON_NEW(unsigned,Ns) : 0 ),
  W                 ( (kmax()+highest_level())>9? max(5,w) : max(4,w) ),
  ST                ( S )
{
  snap_shot()->set_steps(this);                    // set time steps in bodies  
  snap_shot()->add_fields(fieldset::l);            // make sure we have levels  
  for(int n=0; n!=Nsteps(); ++n) N[n] = 0;
  remember();                                      // to be predicted quantities
  set_time_derivs(1,1,0.);                         // set initial forces        
  assign_levels();                                 // get bodies into levels    
  finish_diagnose();                               // finish diagnosis          
  add_to_cpu_step();                               // record CPU time           
  DebugInfo(4,"BlockStepCode constructed\n");
}
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  inline std::ostream& put_char(std::ostream&o, const char c, const int s) {
    if(s>0) for(register int i=0; i<s; ++i) o<<c;
    return o;
  }
}
//------------------------------------------------------------------------------
void BlockStepCode::stats_head(output&to) const {
  SOLVER -> dia_stats_head(to);
  if(to && highest_level())
    for(int i=0, h=-kmax(); i!=Nsteps(); i++, h--)
      if     (h>13)  put_char(to,' ',W-4)<<"2^" <<     h  <<' ';
      else if(h> 9)  put_char(to,' ',W-4)       << (1<<h) <<' ';
      else if(h> 6)  put_char(to,' ',W-4)<<' '  << (1<<h) <<' ';
      else if(h> 3)  put_char(to,' ',W-4)<<"  " << (1<<h) <<' ';
      else if(h>=0)  put_char(to,' ',W-4)<<"   "<< (1<<h) <<' ';
      else if(h==-1) put_char(to,' ',W-4)<<" 1/2 ";
      else if(h==-2) put_char(to,' ',W-4)<<" 1/4 ";
      else if(h==-3) put_char(to,' ',W-4)<<" 1/8 ";
      else if(h==-4) put_char(to,' ',W-4)<<"1/16 ";
      else if(h==-5) put_char(to,' ',W-4)<<"1/32 ";
      else if(h==-6) put_char(to,' ',W-4)<<"1/64 ";
      else if(h>-10) put_char(to,' ',W-4)<<"2^" <<     h  <<' ';
      else           put_char(to,' ',W-5)<<"2^" <<     h  <<' ';
  cpu_stats_head(to);
  if(to) to<<std::endl;
}
#ifdef falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::NBodyCode                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
NBodyCode::NBodyCode(const char*file,              // I: input file             
		     bool       resume,            // I: resume old (if nemo)   
		     fieldset   read_more,         // I: further data to read   
		     const char*time,              // I: time for initial data  
		     fieldset   read_try)          // I: data to try to read    
  falcON_THROWING :
  FILE ( file ),
  SHOT ( fieldset::gravity | read_more ),
  CODE ( 0 ),
  READ ( fieldset::empty )
{
  const fieldset must(fieldset::basic | read_more);
  const fieldset read(must | read_try);
  nemo_in  In(file);                               // open nemo input           
  if(resume) {                                     // IF resuming: last snapshot
    do   SHOT.read_nemo(In,READ,read,0,0);         //   DO:  read bodies        
    while(In.has_snapshot());                      //   WHILE more to be read   
  } else {                                         // ELIF:                     
    bool gotit=false;                              //   read snapshot?          
    do   gotit=SHOT.read_nemo(In,READ,read,time,0);//   DO:  try to read them   
    while(!gotit && In.has_snapshot());            //   WHILE snapshots present 
    if(!gotit)                                     //   didn't read any -> ERROR
      falcON_THROW("NBodyCode: no snapshot matching \"time=%s\""
		   "found in file \"%s\"",time? time:"  ", file);
  }                                                // ENDIF                     
  if(!READ.contain(fieldset::f))                   // UNLESS flags just read    
    SHOT.reset_flags();                            //   reset them              
  if(!READ.contain(must))                          // IF some data missing      
    falcON_THROW("NBodyCode: couldn't read body data: %s",
		 word(READ.missing(must)));
  DebugInfo(4,"NBodyCode constructed\n");
}
//------------------------------------------------------------------------------
void NBodyCode::init(const ForceAndDiagnose         *FS,
		     int                             kmax,
		     int                             Nlev,
		     const BlockStepCode::StepLevels*St,
		     fieldset p, fieldset k, fieldset r,
		     fieldset P, fieldset K, fieldset R) falcON_THROWING
{
  DebugInfo(5,"NBodyCode::init(): called ... \n");
  try {
    if(FS->acc_ext()) SHOT.add_fields(fieldset::q);
    if(Nlev <= 1 || St == 0)
      CODE = static_cast<const Integrator*>
	( new LeapFrogCode(kmax,FS,p,k,r,P,K,R) );
    else
      CODE = static_cast<const Integrator*>
	( new BlockStepCode(kmax,Nlev,FS,St,p,k,r,P,K,R,
			    int(1+std::log10(double(SHOT.N_bodies())))));
  } catch(falcON::exception E) {
    DebugInfo(2,"NBodyCode::init(): caught error \"%s\"\n",E.text());
    falcON_RETHROW(E);
  }
  DebugInfo(4,"NBodyCode::init(): done\n");
}
#endif // falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::ForceDiagGrav                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void ForceDiagGrav::diagnose_grav() const
{
  double m(0.), vin(0.), vex(0.), w[Ndim][Ndim]={0.};
  vect_d x(0.);
  if(snap_shot()->have(fieldbit::q)) {             // IF have external pot      
    LoopAllBodies(snap_shot(),b) {                 //   LOOP bodies             
      register double mi = mass(b);                //     m                     
      m  += mi;                                    //     add: total mass       
      vin+= mi * pot(b);                           //     add: int pot energy   
      vex+= mi * pex(b);                           //     add: ext pot energy   
      register vect_d mx = mi * pos(b);            //     m * x                 
      AddTensor(w,mx,acc(b));                      //     add to W_ij           
      x += mx;                                     //     add: dipole           
    }                                              //   END LOOP                
  } else {                                         // ELSE: no external pot     
    LoopAllBodies(snap_shot(),b) {                 //   LOOP bodies             
      register double mi = mass(b);                //     m                     
      m  += mi;                                    //     add: total mass       
      vin+= mi * pot(b);                           //     add: int pot energy   
      register vect_d mx = mi * pos(b);            //     m * x                 
      AddTensor(w,mx,acc(b));                      //     add to W_ij           
      x += mx;                                     //     add: dipole           
    }                                              //   END LOOP                
  }                                                // ENDIF                     
#ifdef falcON_MPI
  if(snap_shot()->parallel()) {
    const int Num=Ndim*(Ndim+1)+3;
    double Loc[Num];
    double Tmp[Num];
    int p=0;
    Loc[p++] = m;
    Loc[p++] = vin;
    Loc[p++] = vex;
    for(int i=0; i!=Ndim; ++i) {
      Loc[p++] = x[i];
      for(int j=0; j!=Ndim; ++j) 
	Loc[p++] = w[i][j];
    }
    COMMUN(Comm(snap_shot()))->AllReduce(Loc,Tmp,Num,MPI::Sum);
    m    = Tmp[p=0];
    vin  = Tmp[++p];
    vex  = Tmp[++p];
    for(int i=0; i!=Ndim; ++i) {
      x[i] = Tmp[++p];
      for(int j=0; j!=Ndim; ++j) 
	w[i][j] = Tmp[++p];
    }
  }
#endif
  M   = m;                                         // total mass                
  Vin = half*vin;                                  // total int pot energy      
  Vex = vex;                                       // total ext pot energy      
  CMX = x/m;                                       // center of mass            
  for(int i=0; i!=Ndim; ++i)
    for(int j=0; j!=Ndim; ++j) 
      WT[i][j] = half *(w[i][j]+w[j][i]);          // pot energy tensor         
  W    = tr(WT);
  TIME = snap_shot()->time();
}
////////////////////////////////////////////////////////////////////////////////
void ForceDiagGrav::diagnose_vels() const falcON_THROWING
{
  if(snap_shot()->time() != TIME)
    falcON_THROW("ForceDiagGrav::diagnose_vels(): time mismatch");
  double m(0.), k[Ndim][Ndim]={0.};
  vect_d v(0.), l(0.);
  LoopAllBodies(snap_shot(),b) {                   // LOOP bodies               
    register double mi = mass(b);                  //   m                       
    m  += mi;                                      //   add: total mass         
    register vect_d mv = mi * vel(b);              //   m * v                   
    AddTensor(k,mv,vel(b));                        //   add to K_ij             
    v += mv;                                       //   add: total momentum     
    l += vect_d(pos(b)) ^ mv;                      //   add: total ang mom      
  }                                                // END LOOP                  
#ifdef falcON_MPI
  if(snap_shot()->parallel()) {
    const int Num=Ndim*(Ndim+2);
    double Loc[Num];
    double Tmp[Num];
    for(int i=0,p=0; i!=Ndim; ++i) {
      Loc[p++] = l[i];
      Loc[p++] = v[i];
      for(int j=0; j!=Ndim; ++j) 
	Loc[p++] = k[i][j];
    }
    COMMUN(Comm(snap_shot()))->AllReduce(Loc,Tmp,Num,MPI::Sum);
    for(int i=0,p=0; i!=Ndim; ++i) {
      l[i] = Tmp[p++];
      v[i] = Tmp[p++];
      for(int j=0; j!=Ndim; ++j) 
	k[i][j] = Tmp[p++];
    }
  }
#endif
  L   = l;                                         // total angular momentum    
  CMV = v/m;                                       // center of mass velocity   
  for(int i=0; i!=Ndim; ++i)
    for(int j=0; j!=Ndim; ++j) 
      KT[i][j] = half * k[i][j];                   // kin energy tensor         
  T   = tr(KT);                                    // total kin energy          
  TW  =-T/W;                                       // virial ratio              
}
////////////////////////////////////////////////////////////////////////////////
void ForceDiagGrav::diagnose_full() const
{
  double m(0.), vin(0.), vex(0.), w[Ndim][Ndim]={0.}, k[Ndim][Ndim]={0.};
  vect_d x(0.), v(0.), l(0.);
  if(snap_shot()->have(fieldbit::q)) {             // IF have external pot      
    LoopAllBodies(snap_shot(),b) {                 //   LOOP bodies             
      register double mi = mass(b);                //     m                     
      m  += mi;                                    //     add: total mass       
      vin+= mi * pot(b);                           //     add: int pot energy   
      vex+= mi * pex(b);                           //     add: ext pot energy   
      register vect_d mx = mi * pos(b);            //     m * x                 
      register vect_d mv = mi * vel(b);            //     m * v                 
      AddTensor(w,mx,acc(b));                      //     add to W_ij           
      AddTensor(k,mv,vel(b));                      //     add to K_ij           
      x += mx;                                     //     add: dipole           
      v += mv;                                     //     add: total momentum   
      l += mx ^ vect_d(vel(b));                    //     add: total ang mom    
    }                                              //   END LOOP                
  } else {                                         // ELSE: no external pot     
    LoopAllBodies(snap_shot(),b) {                 //   LOOP bodies             
      register double mi = mass(b);                //     m                     
      m  += mi;                                    //     add: total mass       
      vin+= mi * pot(b);                           //     add: int pot energy   
      register vect_d mx = mi * pos(b);            //     m * x                 
      register vect_d mv = mi * vel(b);            //     m * v                 
      AddTensor(w,mx,acc(b));                      //     add to W_ij           
      AddTensor(k,mv,vel(b));                      //     add to K_ij           
      x += mx;                                     //     add: dipole           
      v += mv;                                     //     add: total momentum   
      l += mx ^ vect_d(vel(b));                    //     add: total ang mom    
    }                                              //   END LOOP                
  }                                                // ENDIF                     
#ifdef falcON_MPI
  if(snap_shot()->parallel()) {
    const int Num=Ndim*(2*Ndim+3)+3;
    double Loc[Num];
    double Tmp[Num];
    int p=0;
    Loc[p++] = m;
    Loc[p++] = vin;
    Loc[p++] = vex;
    for(int i=0; i!=Ndim; ++i) {
      Loc[p++] = x[i];
      Loc[p++] = v[i];
      Loc[p++] = l[i];
      for(int j=0; j!=Ndim; ++j) {
	Loc[p++] = w[i][j];
	Loc[p++] = k[i][j];
      }
    }
    COMMUN(Comm(snap_shot()))->AllReduce(Loc,Tmp,Num,MPI::Sum);
    m   = Tmp[p=0];
    vin = Tmp[++p];
    vex = Tmp[++p];
    for(int i=0; i!=Ndim; ++i) {
      x[i] = Tmp[++p];
      v[i] = Tmp[++p];
      l[i] = Tmp[++p];
      for(int j=0; j!=Ndim; ++j) {
	w[i][j] = Tmp[++p];
	k[i][j] = Tmp[++p];
      }
    }
  }
#endif
  M   = m;                                         // total mass                
  Vin = half*vin;                                  // total int pot energy      
  Vex = vex;                                       // total ext pot energy      
  W   = tr(w);                                     // total pot energy from acc 
  T   = half*tr(k);                                // total kin energy          
  TW  =-T/W;                                       // virial ratio              
  L   = l;                                         // total angular momentum    
  CMX = x/m;                                       // center of mass            
  CMV = v/m;                                       // center of mass velocity   
  for(int i=0; i!=Ndim; ++i)
    for(int j=0; j!=Ndim; ++j) {
      WT[i][j] = half *(w[i][j]+w[j][i]);          // pot energy tensor         
      KT[i][j] = half * k[i][j];                   // kin energy tensor         
    }
  TIME = snap_shot()->time();
}
////////////////////////////////////////////////////////////////////////////////
#if 0
void ForceDiagGrav::write_diag_nemo(nemo_out const&out,
				   double         cpu) const
{
  out.open_set(nemo_io::diags);                    // OPEN diagnostics set      
  out.single_vec(0) = Ekin() + Epot();             //     copy total energy     
  out.single_vec(1) = Ekin();                      //     copy kin energy       
  out.single_vec(2) = Epot();                      //     copy pot energy       
  out.write(nemo_io::energy);                      //     write energies        
  for(int i=0;i!=Ndim;++i) for(int j=0;j!=Ndim;++j)//   LOOP dims twice         
    out.single_mat(i,j) = KT[i][j];                //     copy K_ij             
  out.write(nemo_io::KinT);                        //   write K_ij              
  for(int i=0;i!=Ndim;++i) for(int j=0;j!=Ndim;++j)//   LOOP dims twice         
    out.single_mat(i,j) = WT[i][j];                //     copy W_ij             
  out.write(nemo_io::PotT);                        //   write W_ij              
  for(int i=0;i!=Ndim;++i) for(int j=0;j!=Ndim;++j)//   LOOP dimensions         
    out.single_mat(i,j) = as_angmom(L,i,j);        //     copy A_ij             
  out.write(nemo_io::AmT);                         //   write A_ij              
  for(int i=0;i!=Ndim;++i) {                       //   LOOP dims               
    out.single_phs(0,i) = CMX[i];                  //     copy c-of-m pos       
    out.single_phs(1,i) = CMV[i];                  //     copy c-of-m vel       
  }                                                //   END LOOP                
  out.write(nemo_io::cofm);                        //   write c-of-m (x,v)      
  out.write(nemo_io::cputime, (cpu)/60.);          //   write accum CPU [min]   
  out.close_set(nemo_io::diags);                   // CLOSE diagnostics set     
}
#endif
////////////////////////////////////////////////////////////////////////////////
void ForceDiagGrav::dia_stats_head (output& to) const {
  if(to) {
    const char *space = sizeof(real)==4? " " : "     ";
    to  << "      time  "<<space
	<< "    E=T+V    "<<space
	<< "   T     "<<space;
    if(SELF_GRAV)
      to<< "   V_in   "<<space;
    if(acc_ext())
      to<< "   V_ex   "<<space;
    if(SELF_GRAV || acc_ext())
      to<< "   W      "<<space
	<< " -2T/W"<<space;
    to  << "   |L| "<<space
	<< " |v_cm|"<<space;
  }
}
////////////////////////////////////////////////////////////////////////////////
void ForceDiagGrav::dia_stats_line (output&to) const {
  if(to) {
    const char *space = sizeof(real)==4? "-" : "-----";
    to  << " -----------"<<space
	<< "-------------"<<space
	<< "---------"<<space;
    if(SELF_GRAV)
      to<< "----------"<<space;
    if(acc_ext())
      to<< "----------"<<space;
    if(SELF_GRAV || acc_ext())
      to<< "----------"<<space;
    to<< "------"<<space;
    to  << "-------"<<space
	<< "-------"<<space;
  }
}
////////////////////////////////////////////////////////////////////////////////
void ForceDiagGrav::dia_stats_body(output&o) const
{
  if(o) {
    std::ostream&to(o);
    int ACC = 1+sizeof(real);
    std::ios::fmtflags old = to.flags();
    to.setf(std::ios::left | std::ios::showpoint);
    to  << print(TIME,ACC+7,ACC+2) << ' '
	<< print(T+Vin+Vex,ACC+8,ACC+2) << ' '
	<< print(T,ACC+4,ACC-1) << ' ';
    if(SELF_GRAV)
      to<< print(Vin,ACC+5,ACC-1) << ' ';
    if(acc_ext())
      to<< print(Vex,ACC+5,ACC-1) << ' ';
    if(SELF_GRAV || acc_ext())
      to<< print(W,ACC+5,ACC-1) << ' '
	<< print(twice(TW),ACC+1,1) << ' ';
    to  << print(std::sqrt(norm(L)),ACC+2,ACC-3) << ' '
	<< print(std::sqrt(norm(CMV)),ACC+2,ACC-3) << ' ';
    to.flags(old);
  }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::ForceALCON                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void ForceALCON::reset_softening(
				  kern_type ker,
				  real      e
#ifdef falcON_ADAP
				 ,real      ns,
				  unsigned  nr,
				  real      em,
				  real      ef
#endif
				  )
{
#ifdef falcON_ADAP
  NSOFT = ns;
  NREF  = nr;
  EMIN  = abs(em);
  EFAC  = abs(ef);
  if(SOFTENING == individual_adaptive && EFAC  == zero)
    falcON_THROW("ForceALCON: using individual adaptive softening, "
		 " but eps_fac=0\n");
  if(SOFTENING == individual_adaptive && NSOFT == zero)
    falcON_THROW("ForceALCON: using individual adaptive softening, "
		 "but Nsoft=0\n");
#endif
  FALCON.reset_softening(abs(e),ker);
}
////////////////////////////////////////////////////////////////////////////////
ForceALCON::ForceALCON(snapshot          *s,       // I: snapshot: time & bodies
		       real               e,       // I: global softening length
		       real               th,      // I: tolerance parameter    
		       int                nc,      // I: N_crit                 
		       const vect        *x0,      // I: pre-set root centre    
		       kern_type          ke,      // I: softening kernel       
		       real               g,       // I: Newton's G             
		       real               f,       // I: theta_sink/theta       
		       int                ru,      // I: # reused of tree       
		       const acceleration*ae,      // I: external acceleration  
		       const int          gd[4]    // I: direct sum: gravity    
#ifdef falcON_INDI
		       ,soft_type         sf       // I: softening type         
#ifdef falcON_ADAP
		       ,real              ns       // I: N_soft                 
		       ,unsigned          nr       // I: N_ref                  
		       ,real              em       // I: eps_min                
		       ,real              ef       // I: eps_fac                
#endif
#endif
#ifdef falcON_SPH
		       ,const int         sd[3]    //[I: direct sum: SPH]       
#endif
		       ) falcON_THROWING :
  ForceDiagGrav ( s, ae, g!=zero ),
#ifdef falcON_INDI
  SOFTENING     ( sf ),
#endif
  ROOTCENTRE    ( x0 ),
  NCRIT         ( max(1,nc) ),
  REUSE         ( ru ),
  FALCON        ( s, abs(e), abs(th), ke,
#ifdef falcON_INDI
		  SOFTENING != global_fixed,
#endif
		  g, th < zero? const_theta : theta_of_M, f, gd
#ifdef falcON_SPH
		  , sd
#endif
		  ),
  REUSED        ( ru ),
  CPU_TREE      ( 0. ),
  CPU_GRAV      ( 0. ),
  CPU_AEX       ( 0. )
{
#ifdef falcON_MPI
  if(SELF_GRAV && MPI::Initialized())
    falcON_THROW("ForceALCON: cannot (yet) do parallel self-gravity\n");
#endif
#ifdef falcON_INDI
  if(SOFTENING==individual_fixed && !snap_shot()->have(fieldbit::e)) 
    falcON_THROW("ForceALCON: individual fixed softening, but no eps_i given");
#endif
#ifdef falcON_ADAP
  reset_softening(ke,e,ns,nr,em,ef);
#else
  reset_softening(ke,e);
#endif
  DebugInfo(4,"ForceALCON constructed\n");
}
//------------------------------------------------------------------------------
void ForceALCON::set_tree_and_forces(bool all, bool build_tree) const
{
  clock_t cpu = clock();
  // 1. build tree if required
  if(SELF_GRAV || build_tree) {
    if(REUSED < REUSE) {                           // IF may re-use old tree    
      REUSED++;                                    //   increment # re-using    
      FALCON.reuse();                              //   re-use old tree         
    } else {                                       // ELSE                      
      FALCON.grow(NCRIT,ROOTCENTRE);               //   grow new tree           
      REUSED=0;                                    //   reset # re-using        
    }                                              // ENDIF                     
    Integrator::record_cpu(cpu,CPU_TREE);          // record CPU consumption    
  }
  // 2.  deal with self-gravity
  if(SELF_GRAV) {
    // 2.1 compute self-gravity
#ifdef falcON_ADAP
    if(SOFTENING==individual_adaptive)             // IF adaptive softening     
      FALCON.approximate_gravity(true,all,NSOFT,NREF,EMIN,EFAC);
    else                                           // ELIF: fixed softening     
#endif
      FALCON.approximate_gravity(true,all);        // ELIF: fixed softening     
    Integrator::record_cpu(cpu,CPU_GRAV);          // record CPU consumption    
  } else {
    // 2.2 no self-gravity: //   reset pot and acc
    if(acc_ext()) {
      LoopAllBodies(snap_shot(),b)
	if(all || is_active(b))
	  b.pot() = zero;
    } else {
      LoopAllBodies(snap_shot(),b)
	if(all || is_active(b)) {
	  b.pot() = zero;
	  b.acc() = zero;
	}
    }
  }
  // 3 compute external gravity
  if(acc_ext()) {
    acc_ext() -> set(snap_shot(), all, SELF_GRAV? 2 : 0);
    Integrator::record_cpu(cpu,CPU_AEX);           // record CPU consumption    
  }
}
//------------------------------------------------------------------------------
void ForceALCON::cpu_stats_head(output&to) const {
  if(to) {
    if(SELF_GRAV) to << "l2R  D  tree  grav ";
    if(acc_ext()) to << " pext ";
  }
}
//------------------------------------------------------------------------------
void ForceALCON::cpu_stats_line(output&to) const {
  if(to) {
    if(SELF_GRAV) to << "-------------------";
    if(acc_ext()) to << "------";
  }
}
//------------------------------------------------------------------------------
void ForceALCON::cpu_stats_body(output&to) const
{
#ifdef falcON_MPI
  if(snap_shot()->parallel()) {
    double loc[3]={CPU_TREE,CPU_GRAV,CPU_AEX},cpu[3];
    DebugInfo(4,"ForceALCON::cpu_stats_body(): "
	      "calling Communicator::Reduce()\n");
    COMMUN(Comm(snap_shot()))->Reduce(0,loc,cpu,3,MPI::Sum);
    CPU_TREE = cpu[0];
    CPU_GRAV = cpu[1];
    CPU_AEX  = cpu[2];
  }
#endif
  if(to) {
    if(SELF_GRAV) {
      to << std::setw(3) << int(log(FALCON.root_radius())/M_LN2) <<' '
	 << std::setw(2) << FALCON.root_depth() <<' ';
      Integrator::print_cpu(CPU_TREE, to);
      to<<' ';
      Integrator::print_cpu(CPU_GRAV, to);
      to<<' ';
    }
    if(acc_ext()) {
      Integrator::print_cpu(CPU_AEX, to);
      to<<' ';
    }
  }
  CPU_TREE = 0.;
  CPU_GRAV = 0.;
  CPU_AEX  = 0.;
}
////////////////////////////////////////////////////////////////////////////////

// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    src/public/lib/nbody.cc
///
/// \author  Walter Dehnen
///
/// \date    2000-2010
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2010  Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
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
  predALL     ((S->requires()&fieldset::u ? p | fieldset::u : p) | r),
  kickALL     ( S->requires()&fieldset::u ? k | fieldset::v : k ),
  rembALL     ( S->requires()&fieldset::u ? r | fieldset::u : r ),
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
  if( (test=predALL & ~fieldset(fieldset::x|fieldset::u)) ) 
    falcON_Warning("Integration: will not predict '%s'", test.make_word(comp));
  if( (test=kickALL & ~fieldset(fieldset::v)) )
    falcON_Warning("Integration: will not kick '%s'", test.make_word(comp));
  if( (test=rembALL & ~fieldset(fieldset::u)) )
    falcON_Warning("Integration: will not remember '%s'", test.make_word(comp));
  if( predALL & fieldset::u && !(kickALL & fieldset::v) )
    falcON_THROW("Integration: cannot predict w without kicking v");
  if( predALL & fieldset::x && !(kickALL & fieldset::v) )
    falcON_THROW("Integration: request to predict x without kicking v");
  requALL = 
    ((kickALL & fieldset::v) ? fieldset::a : fieldset::empty) |
    ((predALL & fieldset::u) ? fieldset::a : fieldset::empty) ;
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
    (predSPH & fieldset::V ? fieldset::a : fieldset::empty) |
//     predSPH & fieldset::R ? fieldset::D : fieldset::empty |
    (kickSPH & fieldset::U ? fieldset::I : fieldset::empty) ;
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
  inline void move_all(const bodies*B, real dt, bool all) {
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
	b. template datum<TYPE>() += real(dt[level(b)]) * const_datum<DERIV>(b);
      }
    else
      LoopAllBodies(B,b) if(is_active(b)) {
	b. template datum<TYPE>() += real(dt[level(b)]) * const_datum<DERIV>(b);
      }
  }
  //----------------------------------------------------------------------------
  template<int TYPE, int DERIV>
  inline void move_sph(const bodies*B, real dt, bool all) {
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
	b. template datum<TYPE>() += real(dt[level(b)]) * const_datum<DERIV>(b);
    else
      LoopSPHBodies(B,b) if(is_active(b))
	b. template datum<TYPE>() += real(dt[level(b)]) * const_datum<DERIV>(b);
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
  if(predALL & fieldset::u) move_all<fieldbit::u,fieldbit::a>(B,dt,all);
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
  if(rembALL & fieldset::u) copy_all<fieldbit::u,fieldbit::v>(B,all);
#ifdef falcON_SPH
  if(rembSPH & fieldset::V) copy_sph<fieldbit::V,fieldbit::v>(B,all);
  if(rembSPH & fieldset::Y) copy_sph<fieldbit::Y,fieldbit::U>(B,all);
#endif
}
//------------------------------------------------------------------------------
void Integrator::cpu_stats_body(output&to) const
{
  SOLVER->cpu_stats_body(to);
  double C_S=CPU_STEP, C_T=CPU_TOTAL;
#ifdef falcON_MPI
  if(snap_shot()->parallel()) {
    COMMUN(snap_shot()->parallel()->Comm())->Reduce<MPI::Sum>(0,CPU_STEP ,C_S);
    COMMUN(snap_shot()->parallel()->Comm())->Reduce<MPI::Sum>(0,CPU_TOTAL,C_T);
  }
#endif
  if(to) {
    print_cpu    (C_S,to);
    print_cpu_hms(C_T,to<<' ');
#ifdef falcON_MPI
    if(snap_shot()->parallel()) print_cpu_hms(MPI::WallClock(), to<<' ');
#endif
  } 
}
//------------------------------------------------------------------------------
void Integrator::describe(output&to) const {
  if(to) {
    to<<"#"; stats_line(to);
    RunInfo::header(to);
    to.flush();
  }
}
//------------------------------------------------------------------------------
#ifdef falcON_NEMO
//------------------------------------------------------------------------------
void Integrator::write(nemo_out const&o,           // I: nemo output            
		       fieldset       w) const     //[I: what to write]         
{
#ifdef falcON_MPI
  if(snap_shot()->parallel())
    snap_shot()->parallel()->write_nemo(o,w);
  else
#endif
  {
    if( o.is_sink()) return;
    if(!o.is_open()) 
      falcON_THROW("Integrator::write(): nemo device not open\n");
    snap_shot()->write_nemo(o,w);
  } 
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
  }
}
//------------------------------------------------------------------------------
void LeapFrogCode::fullstep(bool rf) const {
  reset_CPU();                                     // reset cpu timers          
  account_new();                                   // account for new bodies    
  if(rf) set_time_derivs(1,1,0.);                  // re-compute initial forces 
  kick(tauh(0));                                   // eg: v+= a*tau/2           
  drift(tau(0));                                   // eg: x+= v*tau;  w+= a*tau 
  set_time_derivs(1,1,tau(0));                     // eg: a = F(x,w)            
  kick(tauh(0));                                   // eg: v+= a*tau/2           
  remember();                                      // eg: w = v                 
  finish_diagnose();                               // finish diagnosis          
  snap_shot()->reset_Nnew();
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
  for(; !(t&1) && l; t>>=1, --l) ;                 // l: lowest level moving    
  bool move = false;                               // need to do anything?      
  for(unsigned i=l; i!=Nsteps(); ++i)              // LOOP levels up to highest 
    if(N[i]) move = true;                          //   IF any non-empty: move  
#ifdef falcON_MPI
  if(snap_shot()->parallel())
    COMMUN(snap_shot()->parallel()->Comm())->
      AllReduceInPlace<MPI::And>(move);
#endif
  if(!move) return;                                // none moving anywhere: DONE
  bool all=true;                                   // are all active?           
  for(int i=0; i!=l; ++i)                          // LOOP lower levels         
    if(N[i]) all = false;                          //   IF all empty: all active
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
    for(unsigned l=0; l!=Nsteps(); ++l)
      N[l] = 0u;
    LoopAllBodies(snap_shot(),b)
      if(!is_new(b)) ++(N[level(b)]);
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
  for(unsigned l=0; l!=Nsteps(); ++l) N[l] = 0;
  LoopAllBodies(B,b)
    ++(N[level(b)]);
}
//------------------------------------------------------------------------------
void BlockStepCode::fullstep(bool rf) const {
  reset_CPU();                                     // reset cpu timers          
  account_new();                                   // account for new bodies    
  account_del();                                   // account for removed bodies
  if(rf) set_time_derivs(1,1,0.);                  // re-compute initial forces 
  remember(true);                                  // remember to be predicted  
  kick_i(tauh(),true);                             // kick by half a step       
  for(int t=0; t!=1<<highest_level(); ++t)         // LOOP elementary steps     
    elementary_step(t);                            //   elementary step         
  finish_diagnose();                               // finish diagnosis          
  snap_shot()->reset_Nnew();
  snap_shot()->reset_Ndel();
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
  for(unsigned n=0; n!=Nsteps(); ++n) N[n] = 0;
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
    if(s>0) for( int i=0; i<s; ++i) o<<c;
    return o;
  }
}
//------------------------------------------------------------------------------
void BlockStepCode::stats_head(output&to) const {
  SOLVER -> dia_stats_head(to);
  if(to && highest_level())
    for(int i=0, h=-kmax(); static_cast<unsigned>(i)!=Nsteps(); i++, h--)
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
NBodyCode::NBodyCode(const char*file,
		     bool       resume,
		     fieldset   read_more,
		     const char*trange,
		     fieldset   read_try) falcON_THROWING :
  FILE ( file ),
  PSHT ( 
#ifdef falcON_MPI
	 MPI::Initialized()? new ParallelSnapshot :
#endif
	 0 ),	     
  SHOT ( 
#ifdef falcON_MPI
	 PSHT? PSHT->local() : 
#endif
	 new snapshot ), 
  CODE ( 0 ),
  READ ( fieldset::empty )
{
  SHOT->add_fields(fieldset::gravity | read_more);
  const fieldset must(fieldset::basic | (read_more-fieldset::k));
  const fieldset read(must | read_try | (read_more&fieldset::k));
  nemo_in In;
#ifdef falcON_MPI
  if(!PSHT || Comm(PSHT)->rank() == 0)
#endif
    In.open(file);
  bool more, gotT=true;
  do {
    gotT =
#ifdef falcON_MPI
      PSHT ?
      PSHT->read_nemo(In,READ,read,resume? 0:trange, 0) :
#endif
      SHOT->read_nemo(In,READ,read,resume? 0:trange, 0) ;
    more = In.has_snapshot();
#ifdef falcON_MPI
    if(PSHT) COMMUN(Comm(PSHT))->BroadCast(0,more);
#endif
    DebugInfo(3,"NBodyCode::NBodyCode: more=%d, resume=%d, gotT=%d\n",
	      more,resume,gotT);
  } while(more && (resume || !gotT));
  if(!gotT)
    falcON_THROW("NBodyCode: no snapshot matching \"time=%s\""
		 "found in file \"%s\"",trange? trange:"  ", file);
  if(!READ.contain(must))
    falcON_THROW("NBodyCode: couldn't read body data: %s",
		 word(READ.missing(must)));
  if(!READ.contain(fieldset::f))
    SHOT->reset_flags();
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
    if(FS->acc_ext()) SHOT->add_fields(fieldset::q);
    if(Nlev <= 1 || St == 0)
      CODE = static_cast<const Integrator*>
	( new LeapFrogCode(kmax,FS,p,k,r,P,K,R) );
    else
      CODE = static_cast<const Integrator*>
	( new BlockStepCode(kmax,Nlev,FS,St,p,k,r,P,K,R,
			    int(1+std::log10(double(SHOT->N_bodies())))));
  } catch(falcON::exception& E) {
    DebugInfo(2,"NBodyCode::init(): caught error \"%s\"\n",E.what());
    falcON_RETHROW(E);
  }
  DebugInfo(4,"NBodyCode::init(): done\n");
}
//------------------------------------------------------------------------------
NBodyCode::~NBodyCode() {
  if(CODE)
    falcON_DEL_O(CODE);
#ifdef falcON_MPI
  if(PSHT)
    falcON_DEL_O(PSHT);
  else
#endif
  if(SHOT)
    falcON_DEL_O(SHOT);
  CODE = 0;
  PSHT = 0;
  SHOT = 0;
}
#endif // falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::ForceDiagGrav                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void ForceDiagGrav::diagnose_grav() const
{
  double m(0.), vin(0.), vex(0.), w[Ndim][Ndim]={{0.}};
  vect_d x(0.);
  if(snap_shot()->have(fieldbit::q)) {             // IF have external pot      
    LoopAllBodies(snap_shot(),b) {                 //   LOOP bodies             
       double mi = mass(b);                //     m                     
      m  += mi;                                    //     add: total mass       
      vin+= mi * pot(b);                           //     add: int pot energy   
      vex+= mi * pex(b);                           //     add: ext pot energy   
      vect_d mx = mi * vect_d(pos(b));             //     m * x                 
      AddTensor(w,mx,acc(b));                      //     add to W_ij           
      x += mx;                                     //     add: dipole           
    }                                              //   END LOOP                
  } else {                                         // ELSE: no external pot     
    LoopAllBodies(snap_shot(),b) {                 //   LOOP bodies             
       double mi = mass(b);                //     m                     
      m  += mi;                                    //     add: total mass       
      vin+= mi * pot(b);                           //     add: int pot energy   
      vect_d mx = mi * vect_d(pos(b));             //     m * x                 
      AddTensor(w,mx,acc(b));                      //     add to W_ij           
      x += mx;                                     //     add: dipole           
    }                                              //   END LOOP                
  }                                                // ENDIF                     
#ifdef falcON_MPI
  if(snap_shot()->parallel()) {
    const int Num=Ndim*(Ndim+1)+3;
    double Tmp[Num];
    int p=0;
    Tmp[p++] = m;
    Tmp[p++] = vin;
    Tmp[p++] = vex;
    for(int i=0; i!=Ndim; ++i) {
      Tmp[p++] = x[i];
      for(int j=0; j!=Ndim; ++j) 
	Tmp[p++] = w[i][j];
    }
    COMMUN(Comm(snap_shot()))->AllReduceInPlace<MPI::Sum>(Tmp,Num);
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
  double m(0.), k[Ndim][Ndim]={{0.}};
  vect_d v(0.), l(0.);
  LoopAllBodies(snap_shot(),b) {                   // LOOP bodies               
     double mi = mass(b);                  //   m                       
    m  += mi;                                      //   add: total mass         
    vect_d mv = mi * vect_d(vel(b));               //   m * v                   
    AddTensor(k,mv,vel(b));                        //   add to K_ij             
    v += mv;                                       //   add: total momentum     
    l += vect_d(pos(b)) ^ mv;                      //   add: total ang mom      
  }                                                // END LOOP                  
#ifdef falcON_MPI
  if(snap_shot()->parallel()) {
    const int Num=Ndim*(Ndim+2);
    double Tmp[Num];
    for(int i=0,p=0; i!=Ndim; ++i) {
      Tmp[p++] = l[i];
      Tmp[p++] = v[i];
      for(int j=0; j!=Ndim; ++j) 
	Tmp[p++] = k[i][j];
    }
    COMMUN(Comm(snap_shot()))->AllReduceInPlace<MPI::Sum>(Tmp,Num);
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
  double m(0.), vin(0.), vex(0.), w[Ndim][Ndim]={{0.}}, k[Ndim][Ndim]={{0.}};
  vect_d x(0.), v(0.), l(0.);
  if(snap_shot()->have(fieldbit::q)) {             // IF have external pot      
    LoopAllBodies(snap_shot(),b) {                 //   LOOP bodies             
       double mi = mass(b);                //     m                     
      m  += mi;                                    //     add: total mass       
      vin+= mi * pot(b);                           //     add: int pot energy   
      vex+= mi * pex(b);                           //     add: ext pot energy   
      vect_d mx = mi * vect_d(pos(b));             //     m * x                 
      vect_d mv = mi * vect_d(vel(b));             //     m * v                 
      AddTensor(w,mx,acc(b));                      //     add to W_ij           
      AddTensor(k,mv,vel(b));                      //     add to K_ij           
      x += mx;                                     //     add: dipole           
      v += mv;                                     //     add: total momentum   
      l += mx ^ vect_d(vel(b));                    //     add: total ang mom    
    }                                              //   END LOOP                
  } else {                                         // ELSE: no external pot     
    LoopAllBodies(snap_shot(),b) {                 //   LOOP bodies             
       double mi = mass(b);                //     m                     
      m  += mi;                                    //     add: total mass       
      vin+= mi * pot(b);                           //     add: int pot energy   
      vect_d mx = mi * vect_d(pos(b));             //     m * x                 
      vect_d mv = mi * vect_d(vel(b));             //     m * v                 
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
    double Tmp[Num];
    int p=0;
    Tmp[p++] = m;
    Tmp[p++] = vin;
    Tmp[p++] = vex;
    for(int i=0; i!=Ndim; ++i) {
      Tmp[p++] = x[i];
      Tmp[p++] = v[i];
      Tmp[p++] = l[i];
      for(int j=0; j!=Ndim; ++j) {
	Tmp[p++] = w[i][j];
	Tmp[p++] = k[i][j];
      }
    }
    COMMUN(Comm(snap_shot()))->AllReduceInPlace<MPI::Sum>(Tmp,Num);
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
void ForceDiagGrav::dia_stats_head (output&to) const {
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
ForceALCON::ForceALCON(snapshot          *s,       // I: snapshot: time & bodies
		       real               e,       // I: global softening length
		       real               th,      // I: tolerance parameter    
		       unsigned           nc,      // I: N_crit                 
		       const vect        *x0,      // I: pre-set root centre    
		       kern_type          ke,      // I: softening kernel       
		       real               g,       // I: Newton's G             
		       real               es,      // I: eps for sink particles
		       real               fs,      // I: theta_sink/theta
		       unsigned           ru,      // I: # reused of tree       
		       const acceleration*ae,      // I: external acceleration  
		       const unsigned     gd[4],   // I: direct sum: gravity    
		       soft_type          sf       // I: softening type         
#ifdef falcON_ADAP
		       ,real              ns       // I: N_soft                 
		       ,unsigned          nr       // I: N_ref                  
		       ,real              em       // I: eps_min                
		       ,real              ef       // I: eps_fac                
#endif
#ifdef falcON_SPH
		       ,const unsigned    sd[3]    //[I: direct sum: SPH]       
#endif
		       ) falcON_THROWING :
  ForceDiagGrav ( s, ae, g!=zero ),
  SOFTENING     ( sf ),
  ROOTCENTRE    ( x0 ),
  NCRIT         ( max(1u,nc) ),
  REUSE         ( ru ),
  FALCON        ( s, abs(e), abs(th), ke,
		  SOFTENING != global_fixed,
		  g, th < zero? const_theta : theta_of_M, abs(es), abs(fs), gd
#ifdef falcON_SPH
		  , sd
#endif
		  ),
  REUSED        ( ru ),
  CPU_TREE      ( 0. ),
  CPU_GRAV      ( 0. ),
  CPU_AEX       ( 0. ),
  _EPS          ( e ),
  _EPSSINK      ( es? es:e ),
  _KERN         ( ke )
{
#ifdef falcON_MPI
  if(SELF_GRAV && MPI::Initialized())
    falcON_THROW("ForceALCON: cannot (yet) do parallel self-gravity\n");
#endif
  if(SOFTENING==individual_fixed && !snap_shot()->have(fieldbit::e)) 
    falcON_THROW("ForceALCON: individual fixed softening, but no eps_i given");
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
  s->add_pointer(&_EPS,"eps");
  s->add_pointer(&_EPSSINK,"epssink");
  s->add_pointer(&_KERN,"kernel");
  s->add_pointer(&FALCON,"forces");
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
      FALCON.approximate_gravity(all,NSOFT,NREF,EMIN,EFAC);
    else                                           // ELIF: fixed softening     
#endif
      FALCON.approximate_gravity(all);             // ELIF: fixed softening     
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
    double cpu[3]={CPU_TREE,CPU_GRAV,CPU_AEX};
    COMMUN(Comm(snap_shot()))->ReduceInPlace<MPI::Sum>(0,cpu,3);
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

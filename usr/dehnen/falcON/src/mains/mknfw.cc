// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// mknfw.cc                                                                    |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002                                               |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// creates N-body initial conditions from a truncated NFW model                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0    02/04/2002  WD  created for nemo 3.0.12                            |
// v 1.0.a  30/08/2002  WD  adapted this file for usage of MPI otherwise       |
// v 1.1    13/11/2002  WD  adapted for various changes in orther files        |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef ALLOW_NEMO                                 // this is a NEMO program    
#error You need "ALLOW_NEMO" to compile __FILE__
#endif
#include <public/Pi.h>
#include <proper/tnfw.h>
#include <body.h>
#include <public/nmio.h>
#include <iomanip>
#include <nemomain.h>
extern "C" {
#include <nemo.h>                                  // NEMO basics & main        
#include <mathfns.h>                               // NEMO random numbers       
}
////////////////////////////////////////////////////////////////////////////////
string defv[] = {
  "out=???\n         output file                                  ",
  "nbody=???\n       number of bodies                             ",
  "v=1\n             maximum circular speed                       ",
  "a=1\n             scale radius:      rho_0 = C/r(r+a)^2        ",
  "b=10\n            truncation radius: rho   = rho_0 sech(r/b)   ",
  "seed=0\n          seed for the randum number generator         ",
  "time=0\n          simulation time of snapshot                  ",
  "zerocm=t\n        centrate velocities of snapshot?             ",
  "WD_units=f\n      input:  kpc, km/s\n"
  "                   output: kpc, kpc/Gyr, G=1 (-> mass unit)     ",
  "give_pot=f\n      give potentials with snapshot                ",
  "outputs=f\n       give some global quantities to stderr        ",
  "VERSION=1.0.a\n   30/August/2002 WD                            ",
  NULL};
string usage =
"mknfw -- construct a truncated NFW model with density\n\n"
"                       sech(r/b)\n"
"           rho(r) = C -----------\n"
"                       r (r+a)^2\n\n";
////////////////////////////////////////////////////////////////////////////////
inline double xrandom_() { return xrandom(0.,1.); }
////////////////////////////////////////////////////////////////////////////////
void nemo::main()
{
  //----------------------------------------------------------------------------
  // 1. create truncated NFW model, i.e. integrate to get Phi(r) etc.           
  //----------------------------------------------------------------------------
  const    double vf  = 0.977775320024919, mf = 2.2228847e5;
  register double vcm = getdparam("v");
  if(getbparam("WD_units")) vcm /= vf;
  nbdy::tnfw_model NFW(getdparam("a"), getdparam("b"), vcm);
  //----------------------------------------------------------------------------
  // 2. create initial conditions from truncated NFW model                      
  //----------------------------------------------------------------------------
  int N = getiparam("nbody");
  nbdy::bodies BB(N);
  init_xrandom(getparam("seed"));
  const    double     m = NFW.total_mass()/double(N);
  register double     rad,vel,cth,R,phi;
  LoopBodies(nbdy::bodies,(&BB),Bi) {
    //                                                                          
    // 2.1 set mass & potential and get r & v                                   
    //                                                                          
    Bi.mass() = m;
    Bi.pot()  = NFW.random(xrandom_, rad, vel);
    //                                                                          
    // 2.2 get x,y,z                                                            
    //                                                                          
    cth       = xrandom(-1.,1.);
    R         = rad * sqrt(1.-cth*cth);
    phi       = xrandom(0.,nbdy::TPi);
    Bi.pos(0) = R * cos(phi);
    Bi.pos(1) = R * sin(phi);
    Bi.pos(2) = rad * cth;
    //                                                                          
    // 2.3 get vx,vy,vz                                                         
    //                                                                          
    cth       = xrandom(-1.,1.);
    R         = vel * sqrt(1.-cth*cth);
    phi       = xrandom(0.,nbdy::TPi);
    Bi.vel(0) = R * cos(phi);
    Bi.vel(1) = R * sin(phi);
    Bi.vel(2) = vel * cth;
  }
  //----------------------------------------------------------------------------
  // 3. centrate velocity on origin                                             
  //----------------------------------------------------------------------------
  if(getbparam("zerocm")) {
    register nbdy::real iN = nbdy::one/real(N);
    register nbdy::vect cvel(nbdy::zero);
    LoopBodies(nbdy::bodies,(&BB),Bi) cvel += nbdy::vel(Bi);
    cvel *= iN;
    LoopBodies(nbdy::bodies,(&BB),Bi) Bi.vel() -= cvel;
  }
  //----------------------------------------------------------------------------
  // 4. output of snapshot                                                      
  //----------------------------------------------------------------------------
  nbdy::nemo_out out(getparam("out"));
  out.write_history();
  nbdy::real time = getdparam("time");
  BB.write_nemo_snapshot(out,&time,
			 getbparam("give_pot")? nbdy::io::mxv | nbdy::io::p :
			                        nbdy::io::mxv);
  //----------------------------------------------------------------------------
  // 5. optional outputs of global quantities                                   
  //----------------------------------------------------------------------------
  if(getbparam("outputs")) {
    std::clog<<
      "# ------------------------------------------------------\n"
      "# truncated NFW model:\n"
      "#\n"
      "#                       sech(r/b)\n"
      "#           rho(r) = C -----------\n"
      "#                       r (r+a)^2\n"
      "#";
    if(getbparam("WD_units"))
      std::clog
	<<"\n# C      = "<<std::setw(13)<<NFW.mass_normal()
	<<" = "<<std::setw(13)<<mf*NFW.mass_normal()   <<" M_sun"
	<<"\n# M_tot  = "<<std::setw(13)<<NFW.total_mass()
	<<" = "<<std::setw(13)<<mf*NFW.total_mass()    <<" M_sun"
	<<"\n# M(r<a) = "<<std::setw(13)<<NFW.cum(NFW.break_radius())
	<<" = "<<std::setw(13)<<mf*NFW.cum(NFW.break_radius()) <<" M_sun"
	<<"\n# M(r<b) = "<<std::setw(13)<<NFW.cum(NFW.trunc_radius())
	<<" = "<<std::setw(13)<<mf*NFW.cum(NFW.trunc_radius()) <<" M_sun"
	<<"\n# a      = "<<std::setw(13)<<NFW.break_radius()
	<<"                 kpc"
	<<"\n# b      = "<<std::setw(13)<<NFW.trunc_radius()
	<<"                 kpc"
	<<"\n# Phi_0  = "<<std::setw(13)<<NFW.central_pot()
	<<" = "<<std::setw(13)<<vf*vf*NFW.central_pot()<<" (km/s)^2"
	<<"\n# vc_max = "<<std::setw(13)<<sqrt(NFW.vcirc_square())
	<<" = "<<std::setw(13)<<vf*sqrt(NFW.vcirc_square())<<" km/s"
	<<"\n# E_tot  = "<<std::setw(13)<<NFW.total_energy()
	<<" = "<<std::setw(13)<<vf*vf*mf*NFW.total_energy()<<" M_sun(km/s)^2";
    else
      std::clog
	<<"\n# C      = "<<     NFW.mass_normal()
	<<"\n# M_tot  = "<<     NFW.total_mass()
	<<"\n# M(r<a) = "<<     NFW.cum(NFW.break_radius())
	<<"\n# M(r<b) = "<<     NFW.cum(NFW.trunc_radius())
	<<"\n# a      = "<<     NFW.break_radius()
	<<"\n# b      = "<<     NFW.trunc_radius()
	<<"\n# Phi_0  =" <<     NFW.central_pot()
	<<"\n# vc_max = "<<sqrt(NFW.vcirc_square())
	<<"\n# E_tot  =" <<     NFW.total_energy();
    std::clog<<"\n# ------------------------------------------------------"
	     <<std::endl;
  }
};

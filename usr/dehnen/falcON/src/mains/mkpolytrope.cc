// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// mkpolytrope.cc                                                              |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2004                                               |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// creates SPH initial conditions from a polytropic gas sphere                 |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 0.0   09/02/2004  WD created (only test for initial poly.h)               |
// v 0.1   10/02/2004  WD setup of masses & positions added & tested           |
// v 1.0   11/02/2004  WD first running version                                |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "1.0"
#define falcON_VERSION_D "11-feb-2004 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#error You need NEMO to compile "src/mains/mkpolytrope.cc"
#endif
#define falcON_RepAction 0                         // no action reporting       
#include <body.h>                                  // bodies                    
#include <public/Pi.h>                             // NEMO file I/O             
#include <public/nmio.h>                           // NEMO file I/O             
#include <proper/poly.h>                           // polytropic sphere         
#include <proper/pack.h>                           // dense sphere packing      
#include <public/tree.h>                           // oct-tree                  
#include <sph/spht.h>                              // SPH stuff                 
#include <main.h>                                  // main & NEMO stuff         
#include <fstream>                                 // C++ file I/O              
#include <iomanip>                                 // C++ I/O formatting        
////////////////////////////////////////////////////////////////////////////////
string defv[] = {
  "out=???\n          output file                                        ",
  "nbody=???\n        number of bodies                                   ",
  "gamma=???\n        polytropic index (for gas)                         ",
  "Nh=42\n            # neighbours -> h_i (must be in {12,18,42,54,78})  ",
  "h_max=0.4\n        max h                                              ",
  "mass_ratio=1\n     m(core)/m(edge) in ]0,infinity[                    ",
  "mass=1\n           total mass                                         ",
  "r_l=1\n            limiting radius                                    ",
  "r_c=\n             core radius (don't use r_l)                        ",
  "time=0\n           simulation time of snapshot                        ",
  "write=mxvHU\n      what to write (you may add N=#partners or omit v)  ",
  "WD_units=f\n       input:  kpc, solar masses\n"
  "                   output: kpc, kpc/Gyr, G=1 (-> mass unit)           ",
  "table=\n           file for table of r,Phi,rho,M(<r)                  ",
  falcON_DEFV, NULL };
string usage = "mkpolytrope -- construct a polytropic sphere";
////////////////////////////////////////////////////////////////////////////////
void nbdy::main()
{
  //----------------------------------------------------------------------------
  // 1. create polytropic sphere, i.e. integrate to get Phi(r) etc.             
  //----------------------------------------------------------------------------
  polytrope PS(getdparam("gamma"));
  //----------------------------------------------------------------------------
  // 2. rescale model to mass and core/limiting radius wanted                   
  //----------------------------------------------------------------------------
  const double vf = 0.977775320024919, mf = 2.2228847e5;
  double mass = getdparam("mass");
  if(getbparam("WD_units")) mass /= mf;
  if(hasvalue("r_c")) PS.reset_scales_core (mass, getdparam("r_c"));
  else                PS.reset_scales_limit(mass, getdparam("r_l"));
  if(hasvalue("table")) PS.write_table(getparam("table"));
  if(0==strcmp(getparam("out"),".")) return;
  //----------------------------------------------------------------------------
  // 3. generate unit sphere densely packed with spheres                        
  //----------------------------------------------------------------------------
  const int N    = getiparam("nbody");
  const io write(getparam("write"));
  const io data(io::mxvf | io::U  | io::H | io::N);
  bodies   BB(N, data, N);
  PackDenseBodies(&BB);
  //----------------------------------------------------------------------------
  // 4. set body masses, re-scale positions, set rho, Uin and sizes             
  //----------------------------------------------------------------------------
  const double m  = mass / double(N);
  const int    Nh = getiparam("Nh");
  double h;
  if     (Nh == 12) h = 0.5*(sqrt(0.5) + 1        );
  else if(Nh == 18) h = 0.5*(1.        + sqrt(1.5));
  else if(Nh == 42) h = 0.5*(sqrt(1.5) + sqrt(2.) ); 
  else if(Nh == 54) h = 0.5*(sqrt(2.)  + 0.5*sqrt(10.));
  else if(Nh == 78) h = 0.5*(0.5*sqrt(10.) + sqrt(3.) );
  else error(" Nh=%f not in [12,18,42,54,78]\n",Nh);
  h *=  
#ifdef falcON_SPH_spline_kernel
    0.5 *
#endif
    cbrt(4*mass/double(N));
  const double a = (getdparam("mass_ratio")-1)/(getdparam("mass_ratio")+1);
  if(a <= -1.) error("mass_ratio <= 0\n");
  const    real   hmax = getdparam("h_max");
  const    double a1   = a+1;
  register double M    = 0.;
  LoopBodies(bodies,&BB,Bi) {
    const double x  = sqrt(norm(pos(Bi)));
    const double xc = cube(x);
    const double mu = a1 - twice(a*xc);
    M        += m * mu;
    Bi.mass() = m * mu;
    const double Mx = mass * xc * (a1 - a*xc);
    Bi.pos() *= x>0.? PS.rad(Mx,polytrope::mass) / x : 1.;
    const double rh = PS.rho(Mx,polytrope::mass);
    Bi.uin()  = PS.Uin(rh,polytrope::density);
    const double z = min(1.,(1-x)/h);
    const double e = z<1.? 4./(z*(3-z*z)+2) : 1.;
    Bi.size() = rh>0.? h * cbrt(e*mu/rh) : hmax;
    Bi.vel()  = zero;
  }
  M = mass/M;
  LoopBodies(bodies,&BB,Bi) {
    Bi.mass() *= M;
    if(size(Bi) > hmax) Bi.size() = hmax;
  }
  //----------------------------------------------------------------------------
  // 5. adjust h_i                                                              
  //----------------------------------------------------------------------------
  oct_tree TREE(&BB,Default::SPHNcrit);
  sph_estimator SPH(&TREE);
  SPH.adjust_sizes(Nh,hmax,1,(4*Nh)/5,(5*Nh)/4);
  //----------------------------------------------------------------------------
  // 6. output of snapshot                                                      
  //----------------------------------------------------------------------------
  nemo_out out(getparam("out"));
  double time = getdparam("time");
  BB.write_nemo_snapshot(out,&time,write);
}

// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// mkcdm.cc                                                                    |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// creates N-body initial conditions from a truncated CDM model                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0    11/06/2002  WD created for nemo 3.0.12                             |
// v 1.1    19/06/2002  WD added constant anisotropy                           |
// v 1.1.a  30/08/2002  WD adapted this file for usage of MPI otherwise        |
// v 2.0    24/09/2002  WD enabled MPI                                         |
// v 2.0.1  13/11/2002  WD typo in tabfile output corrected                    |
// v 2.1    19/04/2004  WD made it work again using nsam.h                     |
// v 2.2    04/05/2004  WD happy icc 8.0; new body.h; new make                 |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "2.2"
#define falcON_VERSION_D "04-may-2004 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile mkcdm
#endif
#define falcON_RepAction 0                         // no action reporting       
#include <public/nsam.h>                           // my N-body sampler         
#include <public/nmio.h>                           // my NEMO file I/O          
#include <public/ionl.h>                           // my I/O utilities          
#include <proper/tcdm.h>                           // my truncated CDM models   
#include <fstream>                                 // C++ file I/O              
#include <iomanip>                                 // C++ I/O formatting        
#include <main.h>                                  // main & NEMO stuff         
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace nbdy;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class TruncCDMSampler                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class TruncCDMSampler :
    private TruncCDMModel,
    public  SphericalSampler
  {
  public:
    TruncCDMSampler(double const&__g,              // I: inner power-law slope  
		    double const&__a,              // I: scale radius           
		    double const&__b,              // I: truncation radius      
		    double const&__v,              // I: maximum v_circ         
		    double const&__beta,           // I: anisotropy             
		    const double*__r,              // I: mass adaption: radii   
		    int    const&__n,              // I: mass adaption: # --    
		    double const&__f,              // I: mass adaption: factor  
		    bool   const&__p) :            // I: mass adaption: R_-/Re  
      TruncCDMModel   ( __g, __a,__b,__v, __beta ),
      SphericalSampler( total_mass(), 0., __beta, __r,__n,__f,__p ) {}
    //--------------------------------------------------------------------------
    TruncCDMModel const&TCDM() const { return *this; }
    //--------------------------------------------------------------------------
    double DF(double q)           const { return TruncCDMModel::g_E(q); }
    double Ps(double r)           const { return TruncCDMModel::psi(r); }
    double rM(double m)           const { return TruncCDMModel::rM(m); }
    double Re(double e)           const { return TruncCDMModel::RcE(e); } 
    double Rp(double e, double l) const { return TruncCDMModel::Rp(e,l*l); }
  };
}                                                  // END: unnamed namespace    
////////////////////////////////////////////////////////////////////////////////
string defv[] = {
  "out=???\n          output file                                        ",
  "nbody=???\n        number of bodies                                   ",
  "g=???\n            inner density exponent                             ",
  "v=1\n              maximum circular speed                             ",
  "a=1\n              scale radius:      rho_0 = C/[r^g(r+a)^(3-g)]      ",
  "b=10\n             truncation radius: rho   = rho_0 sech(r/b)         ",
  "beta=0\n           anisotropy; -1.5 <= beta <= g/2                    ",
//   "r_a=\n             Ossipkov-Merritt anisotropy radius                 ",
  "seed=0\n           seed for the randum number generator               ",
  "q-ran=f\n          use quasi- instead of pseudo-random numbers        ",
  "time=0\n           simulation time of snapshot                        ",
  "f_pos=0.5\n        fraction of bodies with positive sense of rotation ",
  "Rp=\n              for mass adaption: list of R in increasing order   ",
  "fac=1.2\n          for mass adaption: factor between mass bins        ",
  "peri=f\n           for mass adaption: R_peri(E,L) rather than R_c(E)  ",
  "epar=\n            if given, set eps_i = epar * sqrt(m_i/M_tot)       ",
  "WD_units=f\n       input:  kpc, km/s\n"
  "                   output: kpc, kpc/Gyr, G=1 (-> mass unit)           ",
  "outputs=f\n        give some global quantities to stderr              ",
  "tabfile=\n         write table with r, pot, vc, sigma to file         ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage =
"mkcdm -- construct a truncated CDM model with density\n"
"\n"
"                          sech(r/b)\n"
"           rho(r) = C ----------------\n"
"                       r^g (r+a)^(2-g)\n"
"\n"
"         and constant anisotropy beta = 1-(sigma_theta/sigma_r)^2\n"
"\n";
//------------------------------------------------------------------------------
void nbdy::main()
{
  const double vf = 0.977775320024919;             // km/s  in WD_units         
  const double mf = 2.2228847e5;                   // M_sun in WD_units         
  //----------------------------------------------------------------------------
  // 1. set some parameters                                                     
  //----------------------------------------------------------------------------
  if(getdparam("v")<= 0.) error("max v_circ <= 0\n");
  const bool   WD (getbparam("WD_units"));         // using WD_units?           
  const Random Ran(getparam("seed"),6);
  const io data= hasvalue("epar")? io::mxv | io::e : io::mxv;
  const int    nbmax(100);
  double       Rad[nbmax];
  int          nb=0;
  if(hasvalue("Rp")) {
    nb=nemoinpd(getparam("Rp"),Rad,nbmax);
    if(nb+1>nbmax)
      ::error("exceeding expected number of radii in mass adaption");
    Rad[nb] = Rad[nb-1] * 1.e20;
    nb++;
  }
  //----------------------------------------------------------------------------
  // 2. create initial conditions from a Dehnen model using mass adaption       
  //----------------------------------------------------------------------------
  TruncCDMSampler TS(getdparam("g"),
		     getdparam("a"),
		     getdparam("b"),
		     WD? getdparam("v")/vf : getdparam("v"),
		     getdparam("beta"),
		     Rad,nb,
		     getdparam("fac"),
		     getbparam("peri"));
  bodies BB(getiparam("nbody"), data);
  TS.sample(BB,getbparam("q-ran"),Ran,getdparam("f_pos"),getdparam("epar"));
  //----------------------------------------------------------------------------
  // 3. output of snapshot                                                      
  //----------------------------------------------------------------------------
  nemo_out out(getparam("out"));
  double time = getdparam("time");
  BB.write_nemo_snapshot(out,&time,data);
  //----------------------------------------------------------------------------
  // 4. optional outputs of global quantities                                   
  //----------------------------------------------------------------------------
  if(getbparam("outputs")) {
    TruncCDMModel const&CDM (TS.TCDM());
    std::clog<<
      "# ------------------------------------------------------\n"
      "# truncated CDM model:\n"
      "#\n"
      "#                          sech(r/b)\n"
      "#           rho(r) = C ----------------\n"
      "#                       r^g (r+a)^(3-g)\n"
      "#";
    if(getbparam("WD_units"))
      std::clog
	<<"\n# g      = "<<setw(13)<<CDM.inner_slope()
	<<"\n# C      = "<<setw(13)<<CDM.mass_normal()
	<<" = "<<setw(13)<<mf*CDM.mass_normal()   <<" M_sun"
	<<"\n# M_tot  = "<<setw(13)<<CDM.total_mass()
	<<" = "<<setw(13)<<mf*CDM.total_mass()    <<" M_sun"
	<<"\n# M(r<a) = "<<setw(13)<<CDM.cum(CDM.break_radius())
	<<" = "<<setw(13)<<mf*CDM.cum(CDM.break_radius()) <<" M_sun"
	<<"\n# M(r<b) = "<<setw(13)<<CDM.cum(CDM.trunc_radius())
	<<" = "<<setw(13)<<mf*CDM.cum(CDM.trunc_radius()) <<" M_sun"
	<<"\n# a      = "<<setw(13)<<CDM.break_radius()
	<<"                 kpc"
	<<"\n# b      = "<<setw(13)<<CDM.trunc_radius()
	<<"                 kpc"
	<<"\n# Phi_0  = "<<setw(13)<<CDM.central_pot()
	<<" = "<<setw(13)<<vf*vf*CDM.central_pot()<<" (km/s)^2"
	<<"\n# vc_max = "<<setw(13)<<sqrt(CDM.vcirc_square())
	<<" = "<<setw(13)<<vf*sqrt(CDM.vcirc_square())<<" km/s"
	<<"\n# W_tot  =" <<setw(13)<<CDM.pot_energy()
	<<" =" <<setw(13)<<vf*vf*mf*CDM.pot_energy()<<" M_sun(km/s)^2"
	<<"\n# K_tot  = "<<setw(13)<<CDM.kin_energy()
	<<" = "<<setw(13)<<vf*vf*mf*CDM.kin_energy()<<" M_sun(km/s)^2"
	<<"\n# E_tot  =" <<setw(13)<<CDM.total_energy()
	<<" =" <<setw(13)<<vf*vf*mf*CDM.total_energy()<<" M_sun(km/s)^2";
    else
      std::clog
	<<"\n# g      = "<<     CDM.inner_slope()
	<<"\n# C      = "<<     CDM.mass_normal()
	<<"\n# M_tot  = "<<     CDM.total_mass()
	<<"\n# M(r<a) = "<<     CDM.cum(CDM.break_radius())
	<<"\n# M(r<b) = "<<     CDM.cum(CDM.trunc_radius())
	<<"\n# a      = "<<     CDM.break_radius()
	<<"\n# b      = "<<     CDM.trunc_radius()
	<<"\n# Phi_0  =" <<     CDM.central_pot()
	<<"\n# vc_max = "<<sqrt(CDM.vcirc_square())
	<<"\n# W_tot  =" <<     CDM.pot_energy()
	<<"\n# K_tot  = "<<     CDM.kin_energy()
	<<"\n# E_tot  =" <<     CDM.total_energy();
    std::clog<<"\n# ------------------------------------------------------"
	     <<std::endl;
  }
  //----------------------------------------------------------------------------
  // 5. optional output of table                                                
  //----------------------------------------------------------------------------
  if(hasvalue("tabfile")) {
    TruncCDMModel const&CDM (TS.TCDM());
    std::ofstream tab;
    open_error(tab,getparam("tabfile"));
    tab<<"#\n"
       <<"# file created by mkcdm\n"
       <<"# run at  "  <<run_info::time()<<"\n";
    if(run_info::user_known()) tab<<"#     by  \""<<run_info::user()<<"\"\n";
    if(run_info::host_known()) tab<<"#     on  \""<<run_info::host()<<"\"\n";
    if(run_info::pid_known())  tab<<"#     pid  " <<run_info::pid() <<"\n";
    tab<<"#\n"
       <<"# table with properties of isotropic spherical model with density\n"
       <<"#\n"
       <<"#                         sech(r/b)\n"
       <<"#           rho(r) = C ----------------\n"
       <<"#                       r^g (r+a)^(3-g)\n"
       <<"#\n"
       <<"# with g      = "<<getdparam("g")<<"\n"
       <<"#      a      = "<<getdparam("a")<<"\n"
       <<"#      b      = "<<getdparam("b")<<"\n"
       <<"#      vc_max = "<<getdparam("v")<<"\n"
       <<"#      beta   = "<<getdparam("beta")<<"\n"
       <<"#\n"
       <<"#          r       rho(r)       M(r)        "
       <<" psi(r)        vc(r)     sig_r(r)      ln(g(E))\n";
    for(int i=0; i!=CDM.N_grid(); ++i)
      tab<<std::setw(12)<<     CDM.rad_grid(i) <<" "
	 <<std::setw(12)<<     CDM.rho(CDM.rad_grid(i))  <<" "
	 <<std::setw(12)<<     CDM.cum_grid(i) <<" "
	 <<std::setw(12)<<     CDM.psi_grid(i) <<" "
	 <<std::setw(12)<<sqrt(CDM.vcq_grid(i))<<" "
	 <<std::setw(12)<<sqrt(CDM.sqr_grid(i))<<" "
	 <<std::setw(12)<<     CDM.lng_grid(i) <<"\n";
    tab.close();
  }
};
////////////////////////////////////////////////////////////////////////////////

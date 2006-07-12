// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// WD99disk.cc                                                                 |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2005                                          |
//           Paul McMillan, 2004-2005                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
//           paul.mcmillan@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/basic.h>
#include <pjm/WD99disk.h>
#include <utils/Pi.h>
#include <utils/inline_io.h>
#include <public/io.h>
#include <utils/tupel.h>
#include <utils/numerics.h>
#include <utils/random.h>
#include <utils/spline.h>
#include <utils/WDMath.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>

using std::pow;
using std::exp;
using std::log;

using namespace falcON;

//----------------------------------------------------------------------------
// static splines (required throughout)                                   
//----------------------------------------------------------------------------
spline<double> *SPLPOT,*SPLACC;                   // global splines

//---------------------------------------------------------------------------+
// WD99disk::PlanarOrbit                                                     |
//                                                                           |
// A C++ sub-class (technically a structure as the whole thing is "public")  |
// used by WD99disk to do all calculations for an individual orbit in the disk
//                                                                           |
// It is given two inputs specific for any particular orbit, R (the radius)  |
// and Xi (used to determine the eccentricity), and integrates position and  |
// velocity in R over a full period of radial motion                         |
//                                                                           |
// Once this has been done, the function WD99disk::PlanarOrbit::sample       |
// provides a randomly sampled point on the orbit (with randomised phi)      |
//                                                                           |
//---------------------------------------------------------------------------+

WD99disk::PlanarOrbit::PlanarOrbit(double R,
				   double Xi,
				   double Rs,
				   double S0,
				   double Qmin,
				   double Sdens,
				   double sigcorr)
  : re(R),
    xi(Xi),
    rs(Rs),
    s0(S0),
    ire(1./re),
    Q(Qmin),
    sdens(Sdens)
{
  double d2p;
  POT=(*SPLPOT)(log(re));
  CENACC=(*SPLACC)(log(re),&d2p);
  d2p *= ire;
  double rootr=sqrt(re),
    roota=sqrt(CENACC);
  // All the derived parameters needed:
  Ec= 0.5*re*CENACC + POT;               // Energy        (of circular orbit)
  Lc=re*rootr*roota;                     // Ang. Momentum (ditto)
  Omc=roota/rootr;                       // Omega: Circular frequency
  kap=sqrt(d2p+3*CENACC* ire);           // kappa: Epicycle frequency
  gam=2.* Omc/kap;                       // gamma (2*Omega/kappa)
  sigre=    (rs)? s0*exp(-re/rs)  :      // sigma: Velocity disp. at Re 
              3.36*sdens*Q/kap;
  sigre *= sigcorr;
  Lorb=Lc+(sigre*sigre/Omc)*log(xi);     // Ang. Mom. of THIS orbit.
  range=(Lc+Lorb<0)? false : true;       // Try again? 
  
  vectd wo,wn;
  W = new  vectd[30000];  // Use CashKarp integrator to integrate orbit
  dT = 0.001*TPi/kap; 
  W[0][0] = 0.;                                   // Time
  W[0][1] = re;                                   // Radius
  W[0][2] = re*CENACC - Lorb*Lorb/(re*re);        // v_r
  W[0][2] = (W[0][2] >0.)? sqrt(W[0][2]) : 0.;

  wn=W[0];
  if (range){
  for(i=1;(wo[1]-re)*(wn[1]-re) >= 0. || wn[2] < 0.;i++){
    wo=wn;
    CashKarpStep(wn,1.e-14); // integrate over full orbit (plus a tiny bit)
    W[i]=wn;  
    Npoints=i;
  }
  step_back_to_R_equal(re,wn); 
  // Step back so table covers exactly one period.
  } else {Npoints=0;}
  W[Npoints]=wn;
  Npoints++;
 
  ttable  = new double[Npoints];
  Rtable  = new double[Npoints];
  vRtable = new double[Npoints];
  if(range){ 
    for (i=0;i!=Npoints;i++){
      ttable[i]  = W[i][0];
      Rtable[i]  = W[i][1];
      vRtable[i] = W[i][2];
    }
  }
  delete[] W;

  Tr=wn[0];                             // Radial period
  omr=TPi/Tr;                           // Radial frequency 
  g2=kap/omr;                           // Correction factor 
  
}

WD99disk::PlanarOrbit::~PlanarOrbit(){
  if(ttable) delete[] ttable;
  if(Rtable) delete[] Rtable;
  if(vRtable) delete[] vRtable;
}

//----------------------------------------------------------------------------+
// function WD99disk::PlanarOrbit::sample                                     |
// samples a point on the orbit                                               |
//----------------------------------------------------------------------------+

void WD99disk::PlanarOrbit::sample( // samples a point on the orbit
		  Random const&RNG,            // I: random
		  bool   const&q,              // I: quasi?
		  double      &rad,            // O: radius
		  double      &vrad,           // O: radial velocity
		  double      &phi,            // O: azimuth
		  double      &vphi) const    // 0: azimuthal velocity
{

  phi= q? TPi*(RNG(3,0.,1.)) : TPi*(RNG(0.,1.));      // Azimuth

  double t    = q? RNG(4,0.,1.) : RNG(0.,1.);
  t *=Tr;

 rad  = polev(t, ttable,  Rtable, Npoints);
 vrad = polev(t, ttable, vRtable, Npoints);
 vphi = Lorb/rad;

}

//----------------------------------------------------------------------------+
// Integrator for orbits                                                      |
//----------------------------------------------------------------------------+

void WD99disk::PlanarOrbit::CashKarpStep(vectd& W, 
					  const double  eps)
{
  vectd W0=W;
  register double h;
  do {
    h   = CashKarp(W0,dT,eps);
    dT *= h;
  } while(h<0.8);
  W = W0;
}

double WD99disk::PlanarOrbit::CashKarp(vectd& W, 
					const double dt, 
					const double eps)
{
  const double
    b10=1./5,
    b20=3./40,       b21=9./40,
    b30=3./10,       b31=-9./10,   b32=6./5,
    b40=-11./54,     b41=5./2,     b42=-70./27,    b43=35./27,
    b50=1631./55296, b51=175./512, b52=575./13824, b53=44275./110592, 
    b54=253./4096,
    c0 =37./378,     c2 =250./621,     c3 =125./594, c5=512./1771,
    cs0=2825./27648, cs2=18575./48384, cs3=13525./55296,
    cs4=277./14336,  cs5=1./4;
  
  vectd Kx0, Kx1, Kx2, Kx3, Kx4, Kx5, Ws(W);
  
  Kx0=dt*Dt(W); 
  Kx1=dt*Dt(W+b10*Kx0); 
  Kx2=dt*Dt(W+b20*Kx0+b21*Kx1); 
  Kx3=dt*Dt(W+b30*Kx0+b31*Kx1+b32*Kx2); 
  Kx4=dt*Dt(W+b40*Kx0+b41*Kx1+b42*Kx2+b43*Kx3); 
  Kx5=dt*Dt(W+b50*Kx0+b51*Kx1+b52*Kx2+b53*Kx3+b54*Kx4); 
  
  Ws+= cs0*Kx0+cs2*Kx2+cs3*Kx3+cs4*Kx4+cs5*Kx5;
  W += c0*Kx0+c2*Kx2+c3*Kx3+c5*Kx5;
  
  register double h = maxnorm(W-Ws);
  if(h>1.e-19*maxnorm(W) ) h = pow(eps/h,0.2);
  else                            h = 1.e3;
  return min(h,2.);
}

vectd WD99disk::PlanarOrbit::Dt(const vectd& w)
{
  vectd dwdt;
  dwdt[0] = 1;
  dwdt[1] = w[2];
  dwdt[2] =  Lorb*Lorb/(w[1]*w[1]*w[1]) - (*SPLACC)(log(w[1]));
  return dwdt;
}

void WD99disk::PlanarOrbit::step_back_to_R_equal(const double R, vectd& y)
{
  register double dR=R-y[1];
  vectd dy, y0(y);
  y+=( dy=dR*DR(y0)        )/6.;
  y+=( dy=dR*DR(y0+0.5*dy) )/3.;
  y+=( dy=dR*DR(y0+0.5*dy) )/3.;
  y+=( dy=dR*DR(y0+dy)     )/6.;
}

vectd WD99disk::PlanarOrbit::DR(const vectd& y)	// dy/dR
{
  vectd dydR = Dt(y);			      	        // dy/dt
  dydR   *= 1./dydR[2];			   		// dy/dt * (dt/dR)
  return dydR;
}



//-------------------------------------------------------------------------+
// CONSTRUCTOR FOR CLASS WD99disk                                          |
//                                                                         |
// Creates table of Potential and acceleration in the z=0 plane of the     |
// given external field as a function of ln(R)                             |
//                                                                         |
//-------------------------------------------------------------------------+

WD99disk::WD99disk(int    no,                 // # particles/orbit (approx)
		   double rd,                 // Disk scale radius
		   double dens0,              // Disk scale surface density
		   double rsig,               // Vdisp scale radius
		   double qmin,               // Toomre's Q
		   double z0,                 // scale height
		   double eps,                // particle smoothing length
		   const acceleration *pe) :             
  No(no),
  Rd(rd),
  iRd(1./rd),
  Dens0(dens0),
  Qmin(qmin),  // Not a minimum any more. Too lazy to change all the references
  Rsig(rsig),
  sig0(0),
  Z0(z0),
  Eps(eps),
  Mt(TPi*Dens0*Rd*Rd),
  Pe(pe),
  Disk(0, Rd),
  rmin(0.00001*Rd),
  rmax(1000*Rd)
{
  // output tab("RPA.txt");
  n1= int(200*log10(rmax/rmin));
  n=n1+1;
  const double dlr= log(rmax/rmin)/double(n1);
  lr    = new double[n];
  pot   = new double[n];
  dpdr  = new double[n];
  lr[0] = log(rmin);
  for(int i=1; i!=n; ++i)  lr[i] = lr[0] + i * dlr;

 vectd *pos_e  = new vectd[n];
 vectd *acc_e  = new vectd[n];
 for(int i=0; i!=n; ++i) {
   pos_e[i]    = zero;
   pos_e[i][0] = exp(lr[i]);
 }
 // Find Psi and Acc from external potential
 Pe->set(0.,n,0,pos_e,0,0,pot,acc_e,0);
 for(int i=0; i!=n; ++i) { dpdr[i] = -acc_e[i][0]; } 
 double s0  = rmin * dpdr[0],                     // Gradients for pot in ln(r) 
        sn1 = rmax * dpdr[n1];                    // at first and last point 
                                                  // (wanted for spline)         
 // The Splines: -------------------------------------------------------
 SPLPOT = new spline<double>(n,lr,pot,&s0,&sn1);
 SPLACC = new spline<double>(n,lr,dpdr);
 //---------------------------------------------------------------------
 delete[] pos_e;
 delete[] acc_e;
 // If needed, define radial dispersion:
 if (Rsig){
   double d2pdr2,gradp,lrsig,epfreq;
   lrsig=log(Rsig);
   gradp = (*SPLACC)(lrsig,&d2pdr2);
   d2pdr2/=Rsig;
   epfreq = sqrt(d2pdr2+3*gradp/Rsig);
   sig0=3.36*Dens0*Qmin*exp(1-(rsig*iRd))/epfreq;
   warning("FYI, vel. disp. at R=0 is %lf",sig0);
 } 
}
//
// Well that's those tables built. Now what.
//------------------------------------------------------------------------------



//-----------------------------------------------------------------------------+
// function WD99disk::sample                                                   |
// Populates the disk (except in the case where Q=0). See Dehnen 1999 for      |
// details, but in brief: Pick a radius R (in such a way that the final disk   |
// has the surface density desired), define an orbit with 1) energy equal to   |
// that of a circular orbit at radius R, z=0, 2) Angular momentum picked from a| 
// distribution defined so that sig_R is as defined. Place 'Nbpero' bodies on  |
// that orbit. Repeat.                                                         |
//                                                                             |
// Since this is effectively sampling particles in energy, rather than radius, |
// the N-body distribution found will not have exactly the desired density     |
// distribution or velocity dispersion.                                        |
//                                                                             |
// An iterative scheme is used to ensure that the density distribution and     |
// velocity dispersion of the resulting N-body distribution is as close to the |
// required as possible.                                                       |
//-----------------------------------------------------------------------------+

void WD99disk::sample( bodies const&B,            // I/O: bodies to sample
		       int const&NI,              // I: No. iterations
		       bool   const&q,            // I: quasi random?          
		       Random const&RNG) const    // I: pseudo & quasi RNG 
{
  
  int iset=0,Tsize=3,counter=1,lasti=0;
  int Np=B.N_bodies();
  double mpart=Mt/double(Np);
  double Mcum=0.;
  double Re,Xi,z,vz,temporary = 0.1*Rd;

  //--------------------------------------------------------------------------//
  // Iterate to find better target mass function and target velocity disp.    //
  //--------------------------------------------------------------------------//

  for(;temporary<6*Rd;++Tsize) // find No. points in table
    temporary += 0.01*Rd*Rd*exp(temporary*iRd)/temporary;

  // Create tables for R, M, Sigma (aimed, found and used)
  double *RTar,*MInp,*STar,*SInp,*SOut,*VOut;
  RTar = new double[Tsize];  // Radii for all tables
  MInp = new double[Tsize];  // Table of M(R) to pick orbit radii from.
  STar = new double[Tsize];  // Table of target vel disp (sigma_R)
  SInp = new double[Tsize];  // Table of sigma_R correction factor
  SOut = new double[Tsize];  // Table of sigma_R for N-body system given
                             // input SInp*STar
  VOut = new double[Tsize];  // Used in calculation of SOut
  RTar[0] = 0.;
  RTar[1] = 0.05 * Rd;
  RTar[2] = 0.1 * Rd;  //Put in first values by hand, trust me it's easier
  // Find other radii for table, roughly evenly spaced in mass
  for (int i=0;i!=Tsize;i++){ 
    if(i>2) RTar[i] = RTar[i-1] + 0.01*Rd*Rd*exp(RTar[i-1]*iRd)/RTar[i-1];
    MInp[i] = 0.;
    STar[i] = 0.;
    SOut[i] = 0.;
    SInp[i] = 1.;
  }
  RTar[Tsize-1] = 6.*Rd;  //Truncate to avoid numerical trouble from polev

  //--- Create a table of the target sigma_R ------------------//
  for(int i=1; i!=5000000;) {                                //
    Re= Disk.radius(double(i)/5000000.);                     //
    if(Re<RTar[counter]) {                                   //
      double sigtmp = (Rsig)? sig0*exp(-Re/Rsig)     :       //
	              3*(*SPLACC)(log(Re),&temporary);       //
      sigtmp = (Rsig)? sigtmp :                              //
	       3.36*Qmin*Dens0*exp(-Re*iRd)/sqrt((temporary+sigtmp)/Re); 
      STar[counter]+=sigtmp;                                 //
      i++;                                                   //
    } else {                                                 //
      if(i!=lasti) STar[counter]/=double(i-lasti);           //
      lasti=i;                                               //
      counter++;}                                            //
    if (Re>6.*Rd) { break;}                                  //
  }                                                          //
  //---------------------------------------------------------//

  for (int i=0;i != 5*Np;){
    // Sample and bin up. Sample more points than you need to better accuracy.
    do{
      Re= Disk.radius(q? RNG(0,0.,1.) : RNG(0.,1.));  // Radius from exp disk
      Xi= q? RNG(1,0.,1.) : RNG(0.,1.);               // Used for sigma_R
    } while(Re>rmax || Re<rmin);                      // If in pot table range 

    double Densit=Dens0*exp(-Re/Rd);
    PlanarOrbit PO(Re,Xi,Rsig,sig0,Qmin,Densit,1.);  
    // sample an orbit (else try again)
    
    if(PO.range){                           // Then if ang mom is in range
      double corrit=No*(PO.g2);             // Get correction factor
      int Nospecit=int(corrit);
      // get No. points per orbit
      if(q)
	{ Nospecit+=(RNG(2,0.,1.)+Nospecit-corrit>0)? 1: 0;}
      else   Nospecit+=(RNG(0.,1.)+Nospecit-corrit>0)? 1: 0; 
      
      // loop to find points on orbit
      for(int j=0;j!=Nospecit && i!=5*Np;++j,++i) {
	double R,vR,ph,vph;
	PO.sample(RNG, q,R,vR,ph,vph);
	for (int k=1;k<Tsize;++k){
	  if(R<RTar[k]) {
	    VOut[k] +=vR;
	    SOut[k] +=vR*vR;
	    for(int k2=k;k2!=Tsize;++k2) MInp[k2] += 1.;
	    k=Tsize;
	  }
	}
      }
    }
  }

  for(int i=1;i!=Tsize;i++) {
    temporary = 1./(MInp[i]-MInp[i-1]);
    SOut[i] = sqrt((SOut[i]-(VOut[i]*VOut[i]*temporary))*temporary);
    SInp[i] *= STar[i]/SOut[i];}
  for(int i=1;i!=Tsize;i++)
    MInp[i]=(MInp[i])? pow((1 - exp(-(RTar[i])/Rd) - 
			    ((RTar[i])/Rd)*(exp(-(RTar[i])/Rd))),2)/MInp[i]:
      0.;

  temporary = (1.0 - 7.*exp(-6.))/MInp[Tsize-1];
  for(int i=0;i!=Tsize;i++) {
    MInp[i] *= temporary;
    if(!(MInp[i])) 
      MInp[i]=TPi*Rd*Rd*Dens0*(1 - exp(-(RTar[i])/Rd) 
			       - ((RTar[i])/Rd)*(exp(-(RTar[i])/Rd)));
	}

  // We have just used information from sampling to improve accuracy. 
  // See Dehnen '99 for details 
  //                   (though this is not that implementation)

  for (int i=1;i!=NI;++i) iterate(Tsize,Np,q,RNG,RTar,MInp,STar,SInp);

  // Function that iterates to improve this still further

  // Iteration finished
  //----------------------------------------------------------------------
  
  // Now

  for(body Bi=B.begin_all_bodies(); Bi; ) { // until all bodies are sampled
    do{
      double rando=q? RNG(0,0.,1.) : RNG(0.,1.);
      Re= (rando< (1- 7*exp(-6.)))? polev(rando,MInp,RTar,Tsize): 
	// use polynomial interpolator
	Disk.radius(rando);
      Xi= q? RNG(1,0.,1.) : RNG(0.,1.);
    } while(Re>rmax || Re<rmin);                 // If in table range 
    
    double Dens=Dens0*exp(-Re/Rd);
    if(Re>RTar[Tsize-1]) temporary=1.;
    else 
      for(int m=1;m!=Tsize;++m) if(Re<RTar[m] && Re>=RTar[m-1]) temporary= SInp[m];
    
    PlanarOrbit PO(Re,Xi,Rsig,sig0,Qmin,Dens,temporary); // sample an orbit
 
    if(PO.range){                              // Then if ang mom is in range
      double corr=No*(PO.g2);                  // Get correction factor
      int Nospec=int(corr);
      if(q)
	{ Nospec+=(RNG(2,0.,1.)+Nospec-corr>0)? 1: 0;}
      else   Nospec+=(RNG(0.,1.)+Nospec-corr>0)? 1: 0;// get #points per orbit
      
      for(int j=0; j!=Nospec && Bi; ++j,++Bi) {    // loop points on orbit
	double R,vR,ph,vph;
	
	PO.sample(RNG, q,R,vR,ph,vph); // sample planar phase-space coords

	//----------------------------------------------------------------------
	// The vertical component
	//----------------------------------------------------------------------
	if(Z0){
	  double tmp1=q? RNG(5,0.,1.) : RNG(0.,1.);
	  z=0.5*Z0*log((2.*tmp1)/(2.-2.*tmp1)); // sample z from sech^2 disk
	  double tmp2,tmp3,spare;
	  if (iset) {
	    iset=0;
	    vz=spare;
	  } else { 
	    do {
	      tmp1=q? RNG(6,0.,1.) : RNG(0.,1.);
	      tmp2=q? RNG(7,0.,1.) : RNG(0.,1.);     // Sample vz from Gaussian
	      tmp1= 2. * tmp1 - 1.;                  // using trick from rand.cc
	      tmp2= 2. * tmp2 - 1.;                  // (analytical in 2-D)
	      tmp3= tmp1*tmp1 + tmp2*tmp2;
	    } while (tmp3>=1. || tmp3 <=0. );
	    
	    spare = tmp1 * sqrt(-2*log(tmp3)/tmp3);
	    vz = tmp2 * sqrt(-2*log(tmp3)/tmp3);
	  }
	  
	  vz *= sqrt(Pi*Z0*Dens0*exp(-R*iRd));      // Sigma for isothermal disk
	  } else {
	     z = 0.;
	     vz = 0.;}
	//----------------------------------------------------------------------
	
	Bi.mass() = mpart;                         //     set mass           
	Mcum += mpart;                             //     cumulate total mass
	register double cph=cos(ph), sph=sin(ph);  //
	Bi.pos()[0] = R * cph;                     //     x=r*sin(th)*cos(ph)   
	Bi.pos()[1] = R * sph;                     //     y=r*sin(th)*sin(ph)   
	Bi.pos()[2] = z;                           //     z=r*cos(th)   
	
	Bi.vel()[0] = vR * cph - vph * sph;        //     v_x=...   
	Bi.vel()[1] = vR * sph + vph * cph;        //     v_y=...   
	Bi.vel()[2] = vz;                          //     v_z=...   
	
	if(Eps)
	  Bi.eps() = Eps;
      }
      
    }
    
  }
  if(Mcum != Mt) {
    const double mfac = Mt/Mcum;
    LoopAllBodies(&B,Bi)
      Bi.mass()*= mfac;
  }

  delete[] RTar;
  delete[] MInp;
  delete[] STar;
  delete[] SInp;
  delete[] SOut;
  delete[] VOut;
}





///--------------------------------------------------------------------------
// This bit only used in the special case of a completely cold disk


void WD99disk::coldsample( bodies const&B,            // I/O: bodies to sample
			bool   const&q,           // I: quasi random?          
			Random const&RNG) const    // I: pseudo & quasi RNG 
		       

{

  int iset=0;
  int Np=B.N_bodies();
  double mpart=Mt/double(Np);
  double Mcum=0.;
  double Re,z,vz;
  for(body Bi=B.begin_all_bodies(); Bi; ++Bi) { // until all bodies are sampled
    do{
      Re= Disk.radius(q? RNG(0,0.,1.) : RNG(0.,1.));
    }while(Re>rmax || Re<rmin);                 // If in table range 
      double vph=sqrt(Re*(*SPLACC)(log(Re)));
      double phi= q? TPi*(RNG(3,0.,1.)) : TPi*(RNG(0.,1.));      // Azimuth

      //----------------------------------------------------------------------
      // The vertical component
      //----------------------------------------------------------------------
      if(Z0){
	double tmp1=q? RNG(5,0.,1.) : RNG(0.,1.);
	z=0.5*Z0*log((2.*tmp1)/(2.-2.*tmp1)); // sample z from sech^2 disk
	double tmp2,tmp3,spare;
	if (iset) {
	  iset=0;
	  vz=spare;
	} else { 
	  do {
	    tmp1=q? RNG(6,0.,1.) : RNG(0.,1.);
	    tmp2=q? RNG(7,0.,1.) : RNG(0.,1.);       // Sample vz from Gaussian
	    tmp1= 2. * tmp1 - 1.;                    // using trick from rand.cc
	    tmp2= 2. * tmp2 - 1.;                    // (analytical in 2-D)
	    tmp3= tmp1*tmp1 + tmp2*tmp2;
	  } while (tmp3>=1. || tmp3 <0. );
	  
	  spare = tmp1 * sqrt(-2*log(tmp3)/tmp3);
	  vz = tmp2 * sqrt(-2*log(tmp3)/tmp3);
	}
	
	vz *= sqrt(Pi*Z0*Dens0*exp(-Re*iRd));      // Sigma for isothermal disk
      } else {
	z = 0.;
	vz = 0.;

	//----------------------------------------------------------------------
	
	  Bi.mass() = mpart;                          //     set mass           
	  Mcum += mpart;                              //     cumulate total mass
	  register double cph=cos(phi), sph=sin(phi); //
	  Bi.pos()[0] = Re * cph;                     //    x=r*sin(th)*cos(ph) 
	  Bi.pos()[1] = Re * sph;                     //    y=r*sin(th)*sin(ph) 
	  Bi.pos()[2] = z;                            //    z=r*cos(th)   
	  Bi.vel()[0] = -vph * sph;                   //     v_x=...   
	  Bi.vel()[1] = vph * cph;                    //     v_y=...   
	  Bi.vel()[2] = vz;                           //     v_z=...   
      }   
  }
  if(Mcum != Mt) {
    const double mfac = Mt/Mcum;
    LoopAllBodies(&B,Bi)
      Bi.mass()*= mfac;
  }
}

//-----------------------------------------------------------------------------


void WD99disk::iterate(int    const&tsize,
			 int    const&Np,
			 bool   const&q,           // I: quasi random?          
			 Random const&RNG,
			 double *rtar,
			 double *minp,
			 double *star,
			 double *sinp) const
{
  //output tab("RPA.txt");
  // Function called to iterate mass distribution still closer to required
  // Same idea as above.
  
  double rando,Re,Xi,temporary,tempminp[tsize],tmpso[tsize],tmpvo[tsize];

  for(int i=0;i!=tsize;++i){
    tempminp[i]=0.;
    tmpso[i]=0.;
    tmpvo[i]=0.;
  }

  for (int i=0;i != 5*Np;){
    do{
      rando= q? RNG(0,0.,1.) : RNG(0.,1.);
      Re= (rando< (1- 7.*exp(-6.)))? polev(rando,minp,rtar,tsize):
	Disk.radius(rando);
      Xi= q? RNG(1,0.,1.) : RNG(0.,1.);
    } while(Re>rmax || Re<rmin);                 // If in potential table range 
    double Densit=Dens0*exp(-Re/Rd);
    
    if(Re>rtar[tsize-1]) temporary=1.;
    else 
      for(int m=1;m!=tsize;++m) if(Re<rtar[m] && Re>=rtar[m-1]) temporary= sinp[m];
    
     PlanarOrbit PO(Re,Xi,Rsig,sig0,Qmin,Densit,temporary); // sample an orbit 
 
    if(PO.range){                           // Then if ang mom is in range
      double corrit=No*(PO.g2);             // Get correction factor
      int Nospecit=int(corrit);
      // get No. points per orbit
      if(q) { Nospecit+=(RNG(2,0.,1.)+Nospecit-corrit>0)? 1: 0;}
      else    Nospecit+=(RNG(0.,1.)+Nospecit-corrit>0)? 1: 0;
       // loop points on orbit     
      for(int j=0; j!=Nospecit && i != 5*Np; ++j,++i) {
	double R,vR,ph,vph;
	PO.sample(RNG, q,R,vR,ph,vph);
	for (int k=1;k<tsize;++k){
	  if(R<rtar[k]) {
	    tmpvo[k] +=vR;
	    tmpso[k] +=vR*vR;
	    for(int k2=k;k2!=tsize;++k2) tempminp[k2] += 1.;
	    k=tsize;}}
      }
    }
  }
  for(int i=1;i!=tsize;i++) {
    temporary = 1./(tempminp[i]-tempminp[i-1]);
    tmpso[i] = sqrt((tmpso[i]-(tmpvo[i]*tmpvo[i]*temporary))*temporary);
    sinp[i] *= star[i]/tmpso[i];}

  for(int i=1;i!=tsize;i++)
    tempminp[i]=(tempminp[i])? minp[i] * (1 - exp(-(rtar[i])/Rd) - 
					  ((rtar[i])/Rd)*(exp(-(rtar[i])/Rd)))/tempminp[i]:
      0.;

  temporary = (1.0 - 7.*exp(-6.))/tempminp[tsize-1];
  for(int i=0;i!=tsize;i++)  {
    minp[i] = tempminp[i]*temporary;
    if(!(minp[i])) 
      minp[i]=TPi*Rd*Rd*Dens0*(1 - exp(-(rtar[i])/Rd) 
			       - ((rtar[i])/Rd)*(exp(-(rtar[i])/Rd)));
  }
//   tab  <<std::endl;
//   tab.close();
}
//------------------------------------------------------------------------------
//  DESTRUCTOR

WD99disk::~WD99disk()
{

  if(lr)      delete[] lr;
  if(pot)     delete[] pot;
  if(dpdr)    delete[] dpdr;
  if (SPLPOT) delete  SPLPOT;
  if (SPLACC) delete  SPLACC;

}

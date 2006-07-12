// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// anyhalo.cc                                                                  |
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
#include <pjm/anyhalo.h>
#include <utils/Pi.h>
#include <utils/inline_io.h>
#include <public/io.h>
#include <utils/tupel.h>
#include <utils/numerics.h>
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

namespace {
  typedef tupel<4,double> quad_d;                  // quad of doubles          

  using namespace falcON;

  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class HaloModel                                                          //
  //                                                                          //
  // model for any density profile requiring it, it's derivatives, and certain//
  //              knowledge of the very innermost region                      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class HaloModel {
    const double b,g,rs,irs,rt;
    double rsont;
    //--------------------------------------------------------------------------
  public:
    HaloModel(const double _g,
	      const double _b,
	      const double _rs,
	      const double _rt):
      b   (_b),                 // Outer density exponent
      g   (_g),                 // Inner density exponent
      rs  (_rs),                // Scale radius
      irs (1./rs),
      rt  (_rt)                 // Truncation radius (if required)
    {
      rsont=(rt)? rs/rt : 0;
    }
    //--------------------------------------------------------------------------
    double Minner(const double r0) const
    {
      return r0*r0*r0*pow(r0*irs,-g)/(3-g)
	+ r0*r0*r0*(g-b)*pow(r0*irs,1.-g)/(4-g); // M(very small r) / (4*Pi) !!!!!!
    }
    
    double operator() (const double r) const
    {
      register double x=r*irs,
	x1=1.+ x;
      if(rt){return pow(x,-g)*pow(x1,(g-b))*sech(x*rsont);}
      else  return pow(x,-g)*pow(x1,(g-b)); // put density here
    }
    //--------------------------------------------------------------------------
    double operator() (const double r, double &rh1) const
    {
      register double x=r*irs,x1=1.+x, f2=pow(x,-g), f3=pow(x1,(g-b)),
	fp2=-g*pow(x,-g-1), fp3=(g-b)*pow(x1,g-b-1);
      if(rt) {
	rh1 =irs*((fp3*f2 + fp2*f3)*sech(x*rsont) 
		  - rsont*f2*f3*sech(x*rsont)*tanh(x*rsont));  // First deriv.
	return f2*f3*sech(x*rsont);}                    // DENSITY
      else {
	rh1=irs*(fp3*f2 + fp2*f3 ); // First deriv.
	return f2*f3;  }            // DENSITY
    }
    //--------------------------------------------------------------------------
    double operator() (const double r, double &rh1, double &rh2) const
    {
      register double x=r*irs,x1=1.+x, f2=pow(x,-g), f3=pow(x1,(g-b)),
	fp2=-g*pow(x,-g-1), fp3=(g-b)*pow(x1,g-b-1), fpp2=g*(g+1)*pow(x,-g-2),
	fpp3=(g-b)*(g-b-1)*pow(x1,g-b-2);
      if(rt) {
	rh2 =irs*irs*((fpp2*f3 + 2*fp2*fp3 + f2*fpp3)*sech(x*rsont)
		      - 2*(fp3*f2 + fp2*f3)*rsont*sech(x*rsont)*tanh(x*rsont)
		      + f2*f3*rsont*rsont*sech(x*rsont)*(2*tanh(x*rsont) - 1.));
	                                             // Second deriv.
	rh1 = irs*((fp3*f2 + fp2*f3)*sech(x*rsont) 
		  - rsont*f2*f3*sech(x*rsont)*tanh(x*rsont));  // First deriv.
	return f2*f3*sech(x*rsont);}                           // DENSITY
      else{
	rh2 = irs*irs*(fpp2*f3 + 2*fp2*fp3 + f2*fpp3);             // Second deriv.
	rh1 = irs*(fp3*f2 + fp2*f3);                               // First deriv.
	return f2*f3;  }                                         // DENSITY
    }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // auxiliary variables                                                      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  HaloModel               *RHO_;            // full density    
  double                   B_,R_A;          // beta, anisotropy rad.

  inline double red_rho(                                 
			 double x) 
  { 
    if (R_A){

      return pow((1.+(x*x/(R_A*R_A))),-B_+1) * pow(x,2*B_) * (*RHO_)(x);
      
    } else {
      
      return pow(x,2.*B_) * (*RHO_)(x);}   // REDUCED density
  }
  //----------------------------------------------------------------------------
  inline double red_rho(                                     
			 double x,                                
 			 double&d1)                    // First deriv.      
  {
    register double 
      rh1,rh=(*RHO_)(x,    rh1);
  if (R_A){
    register double tmp=1.+(x*x/(R_A*R_A));
    d1 = pow(tmp,-B_+1) * (pow(x,2.*B_) * rh1 + 2.*B_*pow(x,2.*B_-1)*rh)
      + 2.*(-B_+1.)*pow(tmp,-B_)*pow(x,1.+2.*B_)*rh/(R_A*R_A);
    return pow(tmp,-B_+1) * pow(x,2*B_) * rh;
    
  } else {
  
    d1 = pow(x,2*B_) * rh1 + 2*B_ * rh * pow(x,2*B_-1 ); 
    return pow(x,2*B_)*rh;}
  
  }
  //----------------------------------------------------------------------------
  inline double red_rho(                                     
			 double x,                               
			 double&d1,                          
			 double&d2)          // second deriv.             
  {
    register double 
      rh1,rh2,rh=(*RHO_)(x,    rh1,rh2);

 if (R_A){
    register double tmp=1.+(x*x/(R_A*R_A));
  d1 = pow(tmp,-B_+1) * (pow(x,2.*B_) * rh1 + 2.*B_*pow(x,2.*B_-1)*rh)
      + 2.*(-B_+1.)*pow(tmp,-B_)*pow(x,1.+2.*B_)*rh/(R_A*R_A);

  d2 =  pow(tmp,-B_+1)*(pow(x,2.*B_) * rh2 + 4*B_*pow(x,2.*B_-1)*rh1 
		       - 2*B_*(-2*B_+1)*pow(x,2.*B_-2)*rh);
  d2 += 2*(-B_+1)*pow(tmp,-B_)*(pow(x,1+2.*B_)*rh1 + (1+2*B_)*pow(x,+2.*B_)*rh)/(R_A*R_A);
  d2 -= 4*B_*(-B_+1)*pow(tmp,-B_-1)*pow(x,2+2.*B_)*rh/pow(R_A,4);
  d2 += 2*(-B_+1)*pow(tmp,-B_)*(pow(x,1+2.*B_)*rh1 + 2*B_*pow(x,2.*B_)*rh)/(R_A*R_A);

  return pow(tmp,-B_+1) * pow(x,2*B_) * rh;

 } else {
  
    d1 = pow(x,2*B_)*rh1 + 2*B_*rh*pow(x,2*B_-1);
    d2 = -2*B_*(-2*B_+1)*rh*pow(x,2*B_-2)+4*B_*rh1*pow(x,(2*B_-1))+rh2*pow (x,2*B_);
    return pow(x,2*B_)*rh;}
  }


  //============================================================================
  inline double dM(double const&lr, double const&M)// dM/dln r                  
  {
    register double r = exp(lr);                   // ln r as independent var   
    return cube(r)*(*RHO_)(r);
  }
  //----------------------------------------------------------------------------
  // static spline and function returning it.                                   
  //----------------------------------------------------------------------------
  spline<double> *SPLINE, *SPLIN2;                 // global spline             
  inline double glob_spline(const double M)
  {
    return (*SPLINE)(M);
  }
  //----------------------------------------------------------------------------
  // integrand for dPsi/dlnr and d(r^2beta rho sigma_r^2) / dlnr                
  //----------------------------------------------------------------------------
  double TBET;
  //----------------------------------------------------------------------------
  inline quad_d dP(double const&l, quad_d const&F)
  {
    const double x = exp(l), rd=red_rho(x), rh=(*RHO_)(x);
    register quad_d D;
    D[0] =-(*SPLINE)(l) / x;                       // dP       = M/r * dln r    
    D[1] = -rd * (*SPLIN2)(l) / x;                 // d(rd*sq) = rd * M/r dln r 
    D[2] = rh * F[0] * cube(x);                    // dW       = rh * dP * r^3  
    D[3] = F[1] * pow(x,TBET) * TBET;              // dK       = Pr * r^3 (3-2b)
    return D;
  }
}                                                  // END: unnamed namespace    
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class AnyHalo                                                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
AnyHalo::AnyHalo(double MonM,            // Mdisc/Mhalo            
		 double inner,           // inner density exponent
		 double outer,           // outer density exponent
		 double r_s,             // scale radius
		 double r_t,             // trunc. radius
		 double be,              // anisotropy parameter      
		 double ra,              // anisotropy radius
	         const acceleration *MONO) :  //Monopole
  MM   (MonM),
  A    (inner),
  B    ( be ),
  RA   (ra),
  RHO  ( new HaloModel(A,outer,r_s,r_t) ),
  external (MONO)
{
//   output tab("Huh.txt");
  R_A=RA;
  B_=B;
  RHO_ = RHO;
  atb = A-B-B;
  if(B) {
    if( atb < 0 ) error("AnyHalo: beta > 2*alpha: unphysical model\n");
    if( B > 1.0 ) error("AnyHalo: beta > 1: unphysical model\n");
    if( B <-0.5 ) error("AnyHalo: beta < -1/2 not supported\n");
  }
  if(!external)
    error("This is going to go horribly wrong without an external field");

  double extA=3;
  // 1. initialize grids: lr,r,rh,m;                                            
  const double 
    scale = (r_t)? r_t : r_s,
    rmin = 0.001*min(scale,1.),
    rmax = 200.* max(scale,1.);
  n1 = int(200*log10(rmax/rmin));
  n  = 1+n1;
  const double dlr= log(rmax/rmin)/double(n1);
  
  //--------------------------------------------------------------------------
  double M  = RHO->Minner(rmin);
  // Mass at rmin
  //--------------------------------------------------------------------------
  
  double *rd= new double[n];
  lr    = new double[n];
  r     = new double[n];
  mh    = new double[n];
  r [0] = rmin;
  lr[0] = log(rmin);
  rd[0] = red_rho(rmin);
  mh[0] = FPi * M;
  vqm  = mh[0] / r[0];
  for(int i=1; i!=n; ++i) {
    lr[i] = lr[0] + i * dlr;
    r [i] = exp(lr[i]);
    rd[i] = red_rho(r[i]);
    M     = rk4(M,lr[i-1],dlr,dM); //Integrate mass outwards
    mh[i]  = FPi * M;
    if(mh[i] > vqm * r[i]) vqm = mh[i] / r[i];
  }
  //  Scale factors changed to fit the new rules 
  //---------------------------------------------------------------------------
  // SECTION FOR CASE OF EXTERNAL MONOPOLE
  ps = new double[n];
  mg = new double[n];
  rhe = new double[n];

  if(!external) {
    for(int i=0; i!=n; ++i) {
      ps[i] = 0.;
      mg[i] = 0.;
      rhe[i] = 0.;
    }
  } else {
    double *pot_e = new double[n];
    vectd *pos_e  = new vectd[n];
    vectd *acc_e  = new vectd[n];
    for(int i=0; i!=n; ++i) {
      pos_e[i]    = zero;
      pos_e[i][0] = r[i];
    }
    // Find Psi and Acc from external potential
    external->set(0.,n,0,pos_e,0,0,pot_e,acc_e,0);
    // get M_e(<r)
    for(int i=0; i!=n; ++i) {
      ps   [i]  =-pot_e[i];                      // de-scale psi_e
      mg   [i]  =-square(r[i])*acc_e[i][0];      // get M_e(<r)
    }
    delete[] pot_e;
    delete[] pos_e;
    delete[] acc_e;
    extA=(log(mg[1])-log(mg[0]))/dlr;       // 3-inner density exp for external
    double dmdlrin= mg[0] * extA; // assuming power law

    // find density generating external potential
    spline<double> Sext(n,lr,mg,&dmdlrin);
    for(int i=0; i!=n; ++i){
      Sext(lr[i], rhe+i);                  // G*rho = d(G*M)/dlnr               
      rhe[i] /= FPi * cube(r[i]);          //         / (4 * Pi * r^3)          

    }
  }
  //---------------------------------------------------------------------------

  // Descale Disc Mass

  c=mg[n1]/(mh[n1]*MM);
  ic=1/c;
  // 3. potential & radial velocity dispersion, Ec                              

  // ADD psi_h to descaled ps[]  and  mh[i] to mg[i] (which holds external mass)     
  for(int i=0; i!=n; ++i){
    ps[i]  *=ic;
    rhe[i] *= ic;
    mg[i]  *= ic;
    mg[i]  += mh[i];
}

 double
   s0h     = FPi * power<3>(r[0] ) * (*RHO)(r[ 0]),
   sNh     = FPi * power<3>(r[n1]) * (*RHO)(r[n1]),
   s0g     = FPi * power<3>(r[0] ) * ((*RHO)(r[ 0])+rhe[ 0]),
   sNg     = FPi * power<3>(r[n1]) * ((*RHO)(r[n1])+rhe[n1]);

 // Gradients for end points of below
  SPLINE = new spline<double>(n,lr,mh,&s0h,&sNh);   // static spline: m_halo(ln r)
  SPLIN2 = new spline<double>(n,lr,mg,&s0g,&sNg);   // static spline: m_grav(ln r)
  sq     = new double[n];
  ec     = new double[n];
  quad_d y; // Four doubles in a neat wrapper
  y[0]   = mh[n1] / r[n1]; //Psi(halo only) at outer point

  y[1]   = mg[n1]*pow(1,2*B-4)*exp(LogGamma(2*B-4,r[n1])); // Original code
  y[2]   = y[3] = 0.;             // contained these for finding sigma, KE & PE
                                  //   I retained them (while not actually
                                  // adapting them) to save fuss.
  ps[n1]+= y[0];
  sq[n1] = rd[n1]? y[1] / rd[n1] : 0;
  TBET = 3-2*B;

  for(int i=n-2; i!=-1; --i) {
    y     = rk4(y,lr[i+1],-dlr,dP);
    ps[i] += y[0];
    sq[i] = rd[i]? y[1] / rd[i] : 0.;
    ec[i] = ps[i] - 0.5*mg[i]/r[i];
  }
  ps0  = A<2 && extA>1 ? ps[0] + mh[0]/(r[0]*(2-A)) + (mg[0]-mh[0])/(r[0]*(extA-1)): 0.;
  sq0  = sq[0]*rd[0] + FPi/(3-A)*Ipow(r[0],1-A-atb);/// WRONG, too much fuss to fix
  wtot = TPi * y[2]; // WRONG, too much fuss to fix
  ktot =-TPi * y[3]; //       Ditto
  delete   SPLINE;
  delete   SPLIN2;
  delete[] rd;
  for(nm=0; nm!=n1; ++nm) if(mg[nm+1] == mg[nm]) break;
  nm++;
  // 4. compute distibution function                                            
  setgE();
}
////////////////////////////////////////////////////////////////////////////////
namespace {

  //----------------------------------------------------------------------------
  //  Clever way of integrating the equation for the distribution function
  //  N.B. the part for when psi<last point on spline. It IS used, and works
  //  for all density profiles providing red_rho is done right
  //----------------------------------------------------------------------------

  double E_,M_;                             // Eps, Mtot, Mtot*a/b,beta  
  double nu;                                       // 1/(1-p)                   
  inline double intge_of_q(double q)
  {
    //   d Psi       Q^(1-p)                                                    
    // ----------- = ------- dq    with  Psi = Q (1-q^[1/(1-p)])                
    // (Q - Psi)^p    1-p                                                       
    // 
    // N.B. Q and p are both constants in this integration, making life very easy
    register double psi = E_*(1-pow(q,nu));        // Psi = Eps * (1 - q^nu)    
    if     (psi <= 0.)
      return 0.;
    else if(psi < SPLINE->last_X()) {
      if(B_ > 0.5) {
	double tmp1,tmp0=red_rho((M_/psi),tmp1);
	return -M_/(psi*psi)*tmp1;}
      if(B_ >-0.5){
	double tmp2,tmp1,tmp0=red_rho((M_/psi),tmp1,tmp2);
 	return M_/(psi*psi*psi)*(2*tmp1 + M_*tmp2/psi);
      }
    } else
      return (*SPLINE)(psi);
  }

  //----------------------------------------------------------------------------
  inline double intgE(double z)                    // q = z^4; dq = 4 z^3 dz    
  {
    register double f=cube(z);
    return times<4>(f)*intge_of_q(z*f);
  }
}                                                  // END: unnamed namespace   
//------------------------------------------------------------------------------
inline void AnyHalo::setgE()
{
  //output tab("g.tab"); 
  // 1. consider possible cases for beta, tabulate integrand on grid            
  double *in = new double[n];
  if       (B >= 0.5) {        //  0.5 <= beta <  1  :  tabulate d red/ dpsi    
    double rd,rd1,ps1;
    for(int i=0; i!=n; ++i) {
      rd    = red_rho(r[i],rd1);                 // red, red'                  
      ps1   =-mg[i]/square(r[i]);                  // psi'                       
      in[i] = rd1/ps1;                            // dred/dpsi                  
    }
  } else if(B >=-0.5) {        // -0.5 <= beta <  0.5:  tabulate d^2 red/ dpsi^2
    double rh,rd,rd1,rd2,ps1,ps2;
    for(int i=0; i!=n; ++i) {
      rh    = (*RHO)(r[i]);                     // rho                        
      rd    = red_rho(r[i],rd1,rd2);             // red, red', red"            
      ps1   =-mg[i]/square(r[i]);                 // psi'                       
      ps2   =-2*ps1/r[i]-FPi*(rh+rhe[i]);         // psi"                       
      in[i] = (rd2*ps1 - rd1*ps2)/cube(ps1);      // d^2red/dpsi^2              
    }
  }
   else
    error("AnyHalo: beta < -0.5 not supported\n");
  // 2. compute g(E) on grid                                                    
  lg = new double[n];

  const double
    Logalfa = (1.5-B) * LogofTwo + LogBeta(0.5,1-B) - LogofPi +
              log( (B >= 0.5 ? 1 : ( 0.5-B)) *
	           (B >=-0.5 ? 1 : (-0.5-B)) );
  if(B == 0.5  ||  B ==-0.5 ||  B ==-1.5 )
    for(int i=0; i!=n; ++i) lg[i] = log(in[i]) - Logalfa;
  else {
    const int mm = int(1.5-B);
    const double
      p      = 1.5-B-mm,
      p1     = 1-p,
      lfc    = log(sin(p1*Pi)) - LogofPi - Logalfa,
      s0     = (atb==0 && A <2)? (mm+1-1/(2-A))*in[0]*r[0]/mg[0]   :
             (atb==0 && A==2)? -in[0]/FPi                       :
             (A <2)          ? (mm+1-atb/(2-A))*in[0]*r[0]/mg[0] :
             (A==2)          ? atb * in[0]/FPi                    :
                               (mm+1-atb/(2-A))*in[0]/ps[0]     ;
    nu     = 1./p1;
    SPLINE = new spline<double>(n,ps,in,&s0);
    M_     = mg[n1];
    B_     = B;
    tester=0;
    brokeps=ps[n1]-1.;
    double gprev = 1E10000;
    for(register int i=0; i!=n; ++i) {
      E_       = ps[i];
      double err, g = qbulir(intgE,0.,1.,1.e-7,&err,0,50);
      // FUDGE unavoidable
      if(g < 0. )
	{
	  g=gprev;
	  if(!tester){
	    warning("AnyHalo: g(E) < 0 at %lf,", r[i]);
	    warning("distribution function inexact around this radius outwards");
	    brokeps=ps[i];
	    tester++;
	  }
	}
      if (tester!=0) ++tester;
      if(err>1.e-3)
	warning("AnyHalo: very inaccurate integration for g(E)\n");
      lg[i] = lfc + p1*log(E_) + log(nu*g);
      gprev=g;
      // END FUDGE
    }
    delete SPLINE;
  }

  delete[] in;
  //tab  <<std::endl;
  //tab.close();
  fac    = c / pow(c,1.5);
  ln_fac = log(fac);

}
//------------------------------------------------------------------------------
AnyHalo::~AnyHalo()
{
  if(r)   delete[] r;
  if(lr)  delete[] lr;
  if(rhe) delete[] rhe;
  if(mh)  delete[] mh;
  if(mg)  delete[] mg;
  if(ps)  delete[] ps;
  if(sq)  delete[] sq;
  if(lg)  delete[] lg;
  if(ec)  delete[] ec;
  if(RHO) delete   RHO;
}


////////////////////////////////////////////////////////////////////////////////
//          The functions below are those actually called by mkanyhalo        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
inline double AnyHalo::redx(double x) const {
  return red_rho(x);
}

//------------------------------------------------------------------------------
inline double AnyHalo::rhox_h(double x) const {
    return (*RHO)(x);
}
//------------------------------------------------------------------------------
double AnyHalo::rho_h(double rad) const {
  return c * rhox_h(rad);
}
//------------------------------------------------------------------------------
inline double AnyHalo::rhox_e(double x) const {
  // BEWARE points beyond the grid. Consider finding an alternative method
  return external? polev(log(x),lr,rhe,n) : 0.;
}
//------------------------------------------------------------------------------
inline double AnyHalo::rhox_g(double x) const {
  return rhox_h(x) + rhox_e(x);
}
//------------------------------------------------------------------------------
double AnyHalo::rho_g(double rad) const {
  return c * rhox_g(rad);
}
//------------------------------------------------------------------------------
inline double AnyHalo::psix(double x) const {
  if(x<= 0.)   return ps0;
  if(x< r[0] ) {
    if(A< 2.) return  ps[0]-mg[0]/(r[0]*(2-A))*pow(x/r[0],2-A); ///// NEEDS
    if(A==2.) return  FPi * log(x);                             /////   A
    else      return  mg[0]/(r[0]*(A-2))*pow(x/r[0],2-A);       /////  FIX
  }
  if(x> r[n1]) return mg[n1]/x;
  else         return polev(log(x),lr,ps,n);
}
//------------------------------------------------------------------------------
double AnyHalo::psi(double rad) const {
  return c*psix(rad);
}
//------------------------------------------------------------------------------
inline double AnyHalo::cumx_g(double x) const {
  if(x<= 0.)  return 0.;
  if(x<r[0] ) return mg[0]*pow(x/r[0],3-A); // FIX
  if(x>r[n1]) return mg[n1];
  else        return polev(log(x),lr,mg,n);
}
//------------------------------------------------------------------------------
double AnyHalo::cum_g(double rad) const {
  return c*cumx_g(rad);
}
//------------------------------------------------------------------------------
inline double AnyHalo::cumx_h(double x) const {
  if(x<= 0.)  return 0.;
  if(x<r[0] ) return mh[0]*pow(x/r[0],3-A); // Fix
  if(x>r[n1]) return mh[n1];
  else        return polev(log(x),lr,mh,n);
}
//------------------------------------------------------------------------------
double AnyHalo::cum_h(double rad) const {
  return c*cumx_h(rad);
}
//------------------------------------------------------------------------------
inline double AnyHalo::xcum_h(double mr) const {
  if(mr<= 0.)   return 0.;
  if(mr<= mh[0]) return r[0]*pow(mr/mh[0],1/(3-A)); //FIX
  if(mr>  mh[n1]) error("AnyHalo: M>Mtot\n");
  return exp(polev(mr,mh,lr,n));
}
//------------------------------------------------------------------------------
double AnyHalo::rM_h(double M) const {
  return xcum_h(M/c);
}
//------------------------------------------------------------------------------
double AnyHalo::vcq(double rad) const {
  if(rad <= 0.) return 0.;
  return cum_g(rad)/rad;//changed from using cumx
}
//------------------------------------------------------------------------------
double AnyHalo::omq(double rad) const {
  if(rad <= 0.) return 0.;
  return cum_g(rad)/cube(rad);
}
//------------------------------------------------------------------------------
double AnyHalo::kpq(double rad) const {
  if(rad <= 0.) return 0.;
  return omq(rad) + FPi * rho_g(rad);
}
//------------------------------------------------------------------------------
inline double AnyHalo::gamx(double x) const {
  if(x < r[0] ) return 2./sqrt(4-A);  //FIX
  if(x > r[n1]) return 1.;
  double oq = cumx_g(x)/cube(x);
  return 2*sqrt(oq/(oq + FPi*rhox_g(x)));
}
//------------------------------------------------------------------------------
double AnyHalo::gam(double rad) const {
  return gamx(rad);
}
//------------------------------------------------------------------------------
inline double AnyHalo::epcx(double x) const {
  if(x<=0.)    return ps0;
  if(x< r[0])  return psix(x) - 0.5*vcq(x);
  if(x> r[n1]) return 0.5*mg[n1]/x;
  else         return polev(log(x),lr,ec,n);
}
//------------------------------------------------------------------------------
double AnyHalo::Epc(double rad) const {
  return c*epcx(rad);
}
//------------------------------------------------------------------------------
inline double AnyHalo::xepc(double e) const {
  if(e>=ps0)    return 0.;
  if(e> ec[0]) {
    if(A< 2.)   return r[0]*pow((ps[0]-e)/(ps[0]-ec[0]),1/(2-A));// FIX
    if(A==2.)   return exp(0.5+e/FPi);
    else        return r[0]*pow(e/ec[0],1/(2-A));
  }
  if(e< ec[n1]) return 0.5*mg[n1]/e;
  else          return exp(polev(e,ec,lr,n));
}
//------------------------------------------------------------------------------
double AnyHalo::RcE(double e) const {
  return xepc(e/c);
}
//------------------------------------------------------------------------------
inline double AnyHalo::xap(double e, double lq, double ce) const {
  register double
    xc  = xepc(e),
    ecc = sqrt(1.-lq/(cumx_g(xc)*xc));
  return xc * pow(1+ecc*ce, 0.5*gamx(xc));
}
//------------------------------------------------------------------------------
double AnyHalo::Rp(double e, double q) const {
  return xap(e/c, q/c, -1.);
}
//------------------------------------------------------------------------------
double AnyHalo::Ra(double e, double q) const {
  return xap(e/c, q/c, 1.);
}
//------------------------------------------------------------------------------
inline double AnyHalo::sqrx(double x) const {
  if(x <= 0.) return 0.;
  if(x<r[0] ) return (sq0-FPi/(3-A)*Ipow(x,1-A-atb))/redx(x);
  if(x>r[n1]) return mg[n1]*pow(x,2*B-4)*exp(LogGamma(2*B-4,x))/redx(x);//////
              return polev(log(x),lr,sq,n);
}
//------------------------------------------------------------------------------
double AnyHalo::sqr(double rad) const {
  return c*sqrx(rad);
}
//------------------------------------------------------------------------------
double AnyHalo::sqt(double rad) const {
  return c*sqrx(rad)*(1-B);
}
//------------------------------------------------------------------------------
inline double AnyHalo::lngE(double e) const {
  if(e <= ps[n1] || A>2 && e>ps0)
    return -1.e100;
  if(e >ps[0]) {
    if(atb==0)
      return lg[0] -0.5*((A-5)*A+4)/(2-A)*log((ps0-e)/(ps0-ps[0])); // NEEDS A FIX
    if(A  < 2)
      return lg[0] + (B-0.5*(6-A-4*B)/(2-A))*log((ps0-e)/(ps0-ps[0]));
    if(A  ==2)
      return lg[0] + 2*(1-B)*(e-ps[0]);
    if(A  > 2)
      return lg[0] + (B+0.5*(6-A-4*B)/(A-2))*log(e/ps[0]);
  }
  if(e<brokeps && e>ps[n1]){
    return polev(e,ps,lg,n,tester+6);}
  return polev(e,ps,lg,n); // Danger Will Robinson, for small psi.
}
//------------------------------------------------------------------------------
inline double AnyHalo::gofE(double e) const {
  if(e <= ps[n1] || A>2 && e>ps0) return 0.; //FIX?
  return exp(lngE(e));
}
//------------------------------------------------------------------------------
double AnyHalo::g_E(double Eps) const {
   return fac * gofE(Eps/c);
}
//------------------------------------------------------------------------------
double AnyHalo::fEL(double Eps, double Lq) const {
  if(B) return fac * gofE(Eps/c) * pow(Lq/c,-B);
  else  return fac * gofE(Eps/c);
}
////////////////////////////////////////////////////////////////////////////////

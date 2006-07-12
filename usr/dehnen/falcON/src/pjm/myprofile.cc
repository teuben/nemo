//-----------------------------------------------------------------------------+
//                                                                             |
// profile.cc                                                                  |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// This is a non-public part of the code.                                      |
// It is property of its author and not to be made public without his written  |
// consent.                                                                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <pjm/myprofile.h>
#include <body.h>
#include <utils/Pi.h>
#include <utils/numerics.h>
#include <cmath>
#include <assert.h>

using std::log;
using std::pow;
using std::exp;

#define TESTING
#undef  TESTING

#ifdef TESTING
#  include <iostream>                                 // C++ I/O                
#  include <fstream>                                  // C++ file I/O           
#  include <iomanip>
using std::setw;
using std::ofstream;
using std::endl;
#endif

using namespace falcON;
//------------------------------------------------------------------------------
namespace {
  real radius(body const&b) { return abs(pos(b)); }
}
//------------------------------------------------------------------------------
void spherical_model::setup_mass(const bodies      *B,// I: bodies              
				 Array<real>  const&R,// I: body radii, sorted  
				 Array<index> const&I)// I: body indices, sorted
  falcON_THROWING
{
  // 1 set radii, cumulative mass, etc..
  r  = falcON_NEW(double,n);
  lr = falcON_NEW(double,n);
  mr = falcON_NEW(double,n);
  oq = falcON_NEW(double,n);
  lq = falcON_NEW(double,n);
  mt = 0.;
  int jb = 0;
  for(int i=0,ib=k; i!=n; ++i,ib+=k) {
    ib    = min(ib,I.size()-1);
    r [i] = R[ib];
    lr[i] = log(r[i]);
    for(; jb!=ib; ++jb)
      mt += B->mass(I[jb]);   
    mr[i] = mt + 0.5*B->mass(I[ib]);       // Mr excactly at  r_i
    oq[i] = mr[i] / cube(r[i]);
    lq[i] = mr[i]*r[i];
  }
  for(; jb!=I.size(); ++jb)
    mt += B->mass(I[jb]);
  // 3 compute density and kappa^2
  kq = falcON_NEW(double,n);
  rh = falcON_NEW(double,n);
  for(int i=0; i!=n-1; ++i) {
    rh[i] = i==0? 3* mr[i+1] / cube(r[i+1]) :
                  3*(mr[i+1]-mr[i-1])/(cube(r[i+1])-cube(r[i-1]));
    kq[i] = rh[i] + oq[i];
    rh[i]/= FPi;
  }
  rh[n-1] = rh[n-2] * pow(r[n-1]/r[n-2],(log(rh[n-2]/rh[n-3]))/
			                (log(r [n-2]/r [n-3])));
  kq[n-1] = FPi*rh[n-1] + oq[n-1];
  // 4 compute psi, ec
  ps = falcON_NEW(double,n);
  ec = falcON_NEW(double,n);
  ps[n-1] = mr[n-1] / r[n-1];
  ec[n-1] = 0.5*ps[n-1];
  register double dr,a,b=mr[n-1]/square(r[n-1]);
  for(register int i=n-2; i>=0; --i) {
    a     = mr[i] / square(r[i]);
    ps[i] = ps[i+1] + (r[i+1]-r[i]) * 0.5 * (a+b);
    ec[i] = ps[i] - 0.5 * mr[i]/r[i];
    b     = a;
  }
#ifdef TESTING
  ofstream test("smod.mass");
  for(register int i=0; i!=n; ++i)
    test<<setw( 5)<<   i     <<' '
	<<setw(10)<<r [i]    <<' '
	<<setw(10)<<mr[i]    <<' '
	<<setw( 8)<<oq[i]    <<' '
	<<setw( 8)<<rh[i]    <<' '
	<<setw( 8)<<ps[i]    <<' '
	<<setw( 8)<<ec[i]    <<'\n';
  test.close();
#endif
}
//------------------------------------------------------------------------------
void spherical_model::setup_vels(const bodies      *B,// I: bodies              
				 Array<real>  const&R,// I: body radii, sorted  
				 Array<index> const&I)// I: body indices, sorted
  falcON_THROWING
{
  // 1 set cumulative quantities:
  double *mvr  = falcON_NEW(double,n);  // sum_i=0..j m_i * vr_i,    j taken 1/2
  double *mvrq = falcON_NEW(double,n);  // sum_i=0..j m_i * vr_i^2,  j taken 1/2
  double *mvt  = falcON_NEW(double,n);  // sum_i=0..j m_i * vt_i,    j taken 1/2
  double *mvtq = falcON_NEW(double,n);  // sum_i=0..j m_i * vt_i^2,  j taken 1/2
  double *mvp  = falcON_NEW(double,n);  // sum_i=0..j m_i * vp_i,    j taken 1/2
  double *mvpq = falcON_NEW(double,n);  // sum_i=0..j m_i * vp_i^2,  j taken 1/2
  vect_d *mam  = falcON_NEW(vect_d,n);  // sum_i=0..j m_i * L_i,     j taken 1/2
  vect_d Mam(0.), ami;
  double Mvr=0,Mvrq=0,Mvt=0,Mvtq=0,Mvp=0,Mvpq=0, vri,vti,vpi,Ri,mi;
  register int jb = 0;
  vri  = B->pos(I[jb]) * B->vel(I[jb]) / R[jb];
  ami  = B->pos(I[jb]) ^ B->vel(I[jb]);
  Ri   = sqrt(square(B->pos(I[jb])[0])+square(B->pos(I[jb])[1]));
  vpi  = ami[2] / Ri;
  vti  = (B->pos(I[jb])[0]*B->vel(I[jb])[0]+
	  B->pos(I[jb])[1]*B->vel(I[jb])[1])/Ri;
  vti  =(B->pos(I[jb])[2] * vti - Ri * B->vel(I[jb])[2])/R[jb];
  mi   = B->mass(I[jb]);
#ifdef TESTING
  ofstream tcum("smod.cums");
#endif
  for(int i=0,ib=k; i!=n; ++i,ib+=k) {
    ib  = min(ib,I.size()-1);
    for(; jb!=ib; ) {
      Mvr += mi*vri;
      Mvrq+= mi*vri*vri;
      Mvt += mi*vti;
      Mvtq+= mi*vti*vti;
      Mvp += mi*vpi;
      Mvpq+= mi*vpi*vpi;
      Mam += mi*ami;
//       // TEST
//       std::cerr<<" jb="<<jb<<" mi="<<mi
// 	       <<" vri="<<vri<<" Mvr="<<Mvr<<'\n';
//       // TSET
      ++jb;
      vri  = B->pos(I[jb]) * B->vel(I[jb]) / R[jb];
      ami  = B->pos(I[jb]) ^ B->vel(I[jb]);
      Ri   = sqrt(square(B->pos(I[jb])[0])+square(B->pos(I[jb])[1]));
      vpi  = ami[2] / Ri;
      vti  = (B->pos(I[jb])[0]*B->vel(I[jb])[0]+
	      B->pos(I[jb])[1]*B->vel(I[jb])[1])/Ri;
      vti  =(B->pos(I[jb])[2] * vti - Ri * B->vel(I[jb])[2])/R[jb];
      mi   = B->mass(I[jb]);
    }
    mvr [i]= Mvr + 0.5*mi*vri;
    mvrq[i]= Mvrq+ 0.5*mi*vri*vri;
    mvt [i]= Mvt + 0.5*mi*vti;
    mvtq[i]= Mvtq+ 0.5*mi*vti*vti;
    mvp [i]= Mvp + 0.5*mi*vpi;
    mvpq[i]= Mvpq+ 0.5*mi*vpi*vpi;
    mam [i]= Mam + 0.5*mi*ami;
#ifdef TESTING
    tcum<<setw( 5)<<     i     <<' '
	<<setw(10)<<mvr [i]    <<' '
	<<setw(10)<<mvrq[i]    <<' '
	<<setw(10)<<mvt [i]    <<' '
	<<setw(10)<<mvtq[i]    <<' '
	<<setw(10)<<mam [i]    <<'\n';
#endif    
  }
#ifdef TESTING
  tcum.close();
#endif
  for(;;) {
    Mvr += mi*vri;
    Mvrq+= mi*vri*vri;
    Mvt += mi*vti;
    Mvtq+= mi*vti*vti;
    Mvp += mi*vpi;
    Mvpq+= mi*vpi*vpi;
    Mam += mi*ami;
//     // TEST
//     std::cerr<<" jb="<<jb<<" mi="<<mi<<" vri="<<vri<<" Mvr="<<Mvr<<'\n';
//     // TSET
    if(++jb == I.size()) break;
    vri  = B->pos(I[jb]) * B->vel(I[jb]) / R[jb];
    ami  = B->pos(I[jb]) ^ B->vel(I[jb]);
    Ri   = sqrt(square(B->pos(I[jb])[0])+square(B->pos(I[jb])[1]));
    vpi  = ami[2] / Ri;
    vti  = (B->pos(I[jb])[0]*B->vel(I[jb])[0]+
	    B->pos(I[jb])[1]*B->vel(I[jb])[1])/Ri;
    vti  =(B->pos(I[jb])[2] * vti - Ri * B->vel(I[jb])[2])/R[jb];
    mi   = B->mass(I[jb]);
  }
  // 2 compute mean vr, sigma_r, sigma_t, mean angular momentum
  vr = falcON_NEW(double,n);
  vt = falcON_NEW(double,n);
  vp = falcON_NEW(double,n);
  sr = falcON_NEW(double,n);
  st = falcON_NEW(double,n);
  sp = falcON_NEW(double,n);
  am = falcON_NEW(vect_d,n);
  register double
    dM  = mr[1],
    idM = 1./dM;
  vr[0] = mvr[1] * idM;
  vt[0] = mvt[1] * idM;
  vp[0] = mvp[1] * idM;
  am[0] = mam[1] * idM;
  sr[0] = sqrt(dM*mvrq[1] - square(mvr[1])) * idM;
  st[0] = sqrt(dM*mvtq[1] - square(mvt[1])) * idM;
  sp[0] = sqrt(dM*mvpq[1] - square(mvp[1])) * idM;
  for(int i=1; i!=n-1; ++i) {
    dM    = mr[i+1]-mr[i-1];
    idM   = 1./dM;
    vr[i] = (mvr[i+1]-mvr[i-1]) * idM;
    vt[i] = (mvt[i+1]-mvt[i-1]) * idM;
    vp[i] = (mvp[i+1]-mvp[i-1]) * idM;
    am[i] = (mam[i+1]-mam[i-1]) * idM;
    sr[i] = sqrt(dM*(mvrq[i+1]-mvrq[i-1]) - square(mvr[i+1]-mvr[i-1])) * idM;
    st[i] = sqrt(dM*(mvtq[i+1]-mvtq[i-1]) - square(mvt[i+1]-mvt[i-1])) * idM;
    sp[i] = sqrt(dM*(mvpq[i+1]-mvpq[i-1]) - square(mvp[i+1]-mvp[i-1])) * idM;
  }
  dM      = mt-mr[n-2];
  idM     = 1./dM;
  vr[n-1] = (Mvr-mvr[n-2]) * idM;
  vt[n-1] = (Mvt-mvt[n-2]) * idM;
  vp[n-1] = (Mvp-mvp[n-2]) * idM;
  am[n-1] = (Mam-mam[n-2]) * idM;
  sr[n-1] = sqrt(dM*(Mvrq-mvrq[n-2]) - square((Mvr-mvr[n-2]))) * idM;
  st[n-1] = sqrt(dM*(Mvtq-mvtq[n-2]) - square((Mvt-mvt[n-2]))) * idM;
  sp[n-1] = sqrt(dM*(Mvpq-mvpq[n-2]) - square((Mvp-mvp[n-2]))) * idM;
  falcON_DEL_A(mvr);
  falcON_DEL_A(mvrq);
  falcON_DEL_A(mvt);
  falcON_DEL_A(mvtq);
  falcON_DEL_A(mvp);
  falcON_DEL_A(mvpq);
  falcON_DEL_A(mam);
#ifdef TESTING
  ofstream test("smod.vels");
  for(register int i=0; i!=n; ++i)
    test<<setw( 5)<<   i     <<' '
	<<setw(10)<<r [i]    <<' '
	<<setw(10)<<vr[i]    <<' '
	<<setw(10)<<vt[i]    <<' '
	<<setw(10)<<vp[i]    <<' '
	<<setw(10)<<sr[i]    <<' '
	<<setw(10)<<st[i]    <<' '
	<<setw(10)<<sp[i]    <<' '
	<<setw(10)<<am[i]    <<'\n';
  test.close();
#endif
}
//------------------------------------------------------------------------------
spherical_model::spherical_model(const bodies*B,
				 int          W,
				 bool         V) falcON_THROWING
  : k(W/2), n(B->N_bodies()/k)
{
  if(!B->have_all(fieldset(fieldset::m|fieldset::x)))
    error("spherical_model: need body masses & positions\n");
  Array<real>  rad;
  Array<index> ind;
  B->sorted(ind,rad,radius);
  setup_mass(B,rad,ind);
  if( V && B->have(fieldbit::v) )
    setup_vels(B,rad,ind);
  else
    vr = 0;
}
//------------------------------------------------------------------------------
inline const double* spherical_model::table(const arg a) const
{
  switch(a) {
  case radial: return r;
  case lograd: return lr;
  case psient: return ps;
  case epcirc: return ec;
  case lcircq: return lq;
  case omegaq: return oq;
  }
  return 0;
}
//------------------------------------------------------------------------------
double spherical_model::R_c (const double&x, const arg a) const {
  if(a == radial) return x;
  return polev(x,table(a),r,n);
}
//------------------------------------------------------------------------------
double spherical_model::Eps_c  (const double&x, const arg a) const {
  if(a == epcirc) return x;
  return polev(x,table(a),ec,n);
}
//------------------------------------------------------------------------------
double spherical_model::L_c_sq (const double&x, const arg a) const {
  if(a == lcircq) return x;
  return polev(x,table(a),lq,n);
}
//------------------------------------------------------------------------------
double spherical_model::Mofr(const double&x, const arg a) const {
  return polev(x,table(a),mr,n);
}
//------------------------------------------------------------------------------
double spherical_model::Psi(const double&x, const arg a) const {
  if(a == psient) return x;
  return polev(x,table(a),ps,n);
}
//------------------------------------------------------------------------------
double spherical_model::v_c_sq (const double&x, const arg a) const {
  return polev(x,table(a),mr,n)/R_c(x,a);
}
//------------------------------------------------------------------------------
double spherical_model::density(const double&x, const arg a) const {
  return polev(x,table(a),rh,n);
}
//------------------------------------------------------------------------------
double spherical_model::Omega_sq(const double&x, const arg a) const {
  if(a == omegaq) return x;
  return polev(x,table(a),oq,n);
}
//------------------------------------------------------------------------------
double spherical_model::kappa_sq(const double&x, const arg a) const {
  return polev(x,table(a),kq,n);
}
//------------------------------------------------------------------------------
double spherical_model::gamma_sq(const double&x, const arg a) const {
  return 4*Omega_sq(x,a)/kappa_sq(x,a);
}
//------------------------------------------------------------------------------
double spherical_model::mean_vr(const double&x, const arg a) const {
  return polev(x,table(a),vr,n);
}
//------------------------------------------------------------------------------
double spherical_model::mean_vth(const double&x, const arg a) const {
  return polev(x,table(a),vt,n);
}
//------------------------------------------------------------------------------
double spherical_model::mean_vphi(const double&x, const arg a) const {
  return polev(x,table(a),vp,n);
}
//------------------------------------------------------------------------------
double spherical_model::sigma_r(const double&x, const arg a) const {
  return polev(x,table(a),sr,n);
}
//------------------------------------------------------------------------------
double spherical_model::sigma_th(const double&x, const arg a) const {
  return polev(x,table(a),st,n);
}
//------------------------------------------------------------------------------
double spherical_model::sigma_phi(const double&x, const arg a) const {
  return polev(x,table(a),sp,n);
}
//------------------------------------------------------------------------------
vect_d spherical_model::angmom(const double&x, const arg a) const {
  return polev(x,table(a),am,n);
}
//------------------------------------------------------------------------------
inline
void spherical_model::binup_v(int b, int w, const double*v, const double*s,
			      double&V, double&S) const
{
  double m=0.,mv=0.,mvq=0.;
  for(int i=b-w; i<=b+w; i+=2) {
    register double mi = mr[i+1]-mr[i-1];
    m   += mi;
    mv  += mi * v[i];
    mvq += mi * (s[i]*s[i] + v[i]*v[i]);
  }
  V = mv/m;
  S = sqrt(mvq*m-mv*mv)/m;
}
//------------------------------------------------------------------------------
void spherical_model::binup_vr(int b, int w, double&v, double&s) const {
  binup_v(b,w,vr,sr,v,s);
}
//------------------------------------------------------------------------------
void spherical_model::binup_vt(int b, int w, double&v, double&s) const {
  binup_v(b,w,vt,st,v,s);
}
//------------------------------------------------------------------------------
void spherical_model::binup_vp(int b, int w, double&v, double&s) const {
  binup_v(b,w,vp,sp,v,s);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::spherical_profile                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
spherical_profile::spherical_profile(const bodies*B,
				     int          __nmin,
				     double       __dmax,
				     bool         do_vels) falcON_THROWING :
  kmin(__nmin/2), dmax(0.5*__dmax), vr(0)
{
  if(!B->have_all(fieldset(fieldset::m|fieldset::x)))
    falcON_THROW("spherical_profile: need \"mx\" to estimate density\n");
  // 1. sort radii of all bodies                                                
  Array<real>  R;
  Array<index> I;
  B->sorted(I,R,radius);
  const int nb = I.size();
  // 2. step through the bodies to make radial bins                             
  n = nb/kmin;                                     // maximum possible #bins    
  int m=0;                                         // index of bin              
  int*ir= falcON_NEW(int,n);                       // rank of largest radius    
  const double fac=exp(dmax);                      // radius factor between bins
  ir[0] = kmin;                                    // first rank = kmin         
  for(m=1; ; ++m) {                                // LOOP radial bins          
    assert(m<n);
    int j = ir[m-1];                               //   begin of bin            
    const double rm = R[j] * fac;                  //   maximum radius in bin   
    j    += kmin;                                  //   minimum rank in bin     
    while(R[j++] < rm && j < nb);                  //   increase rank to r=rmax 
    if(j + kmin >= nb) j = nb;                     //   enough left for next?   
    ir[m] = min(nb-1, j-1);                        //   highest rank in bin     
    if(ir[m] == nb-1) break;                       //   all bodies used? DONE!  
  }                                                // END LOOP                  
  n = m-1;                                         // this is the true #bins    
  // 3. establish tables with median radii and cumulative masses at them        
  rr = falcON_NEW(double,n);                       // median radius of bin      
  mr = falcON_NEW(double,n);                       // M(<=r) at median radius   
  double mcum=0.;                                  // cumulate mass             
  int j=0;                                         // body index                
  for(int i=0; i!=n; ++i) {                        // LOOP bins                 
    register int 
      s  = i==0? ir[i+1]+1 : ir[i+1]+ir[i-1],      //   lower + upper rank      
      sh = s/2;                                    //   mean rank               
    for(; j!=sh; ++j)                              //   LOOP ranks up to mean   
      mcum += B->mass(I[j]);                       //     cumulate mass         
    if(s&1) {                                      //   IF #radii is odd        
      rr[i] = 0.5*(R[j]+R[j+1]);                   //     median radius         
      mr[i] = mcum +     B->mass(I[j]);            //     cumulative mass at ---
    } else {                                       //   ELSE (#radii is even)   
      rr[i] = R[j];                                //     median radius         
      mr[i] = mcum + 0.5*B->mass(I[j]);            //     cumulative mass at ---
    }                                              //   ENDIF                   
  }                                                // END LOOP                  
  for(; j!=nb; ++j)                                // LOOP ranks up to end      
    mcum += B->mass(I[j]);                         //   cumulate mass           
  mt = mcum;                                       // total mass                
  // 4. for the density we need to sum over the bin contents ...                
  rh = falcON_NEW(double,n);                       // density in bin            
  if(do_vels && B->have(fieldbit::v)) {
    vr = falcON_NEW(double,n);
    vt = falcON_NEW(double,n);
    vp = falcON_NEW(double,n);
    sr = falcON_NEW(double,n);
    st = falcON_NEW(double,n);
    sp = falcON_NEW(double,n);
    skr= falcON_NEW(double,n);
    skt= falcON_NEW(double,n);
    skp= falcON_NEW(double,n);
    kur= falcON_NEW(double,n);
    kut= falcON_NEW(double,n);
    kup= falcON_NEW(double,n);   
    am = falcON_NEW(vect_d,n);
    for(int i=0; i!=n; ++i) {                      // LOOP bins                 
      double M    = 0.;
      double Mvr  = 0.;
      double Mvrq = 0.;
      double Mvrc = 0.;
      double Mvrf = 0.;
      double Mvp  = 0.;
      double Mvpq = 0.;
      double Mvpc = 0.;
      double Mvpf = 0.;
      double Mvt  = 0.;
      double Mvtq = 0.;
      double Mvtc = 0.;
      double Mvtf = 0.;
      vect_d Mam   (0.);
      for(j=i? ir[i-1]:0; j<=ir[i+1]; ++j) {       // LOOP bodies in bin        
	vect_d
	  ami = B->pos(I[j]) ^ B->vel(I[j]);
	double
	  mi  = (i!=0 && j==ir[i-1]  ||  j==ir[i+1]) ?
	        0.5 * B->mass(I[j]) : B->mass(I[j]),
	  vri = B->pos(I[j]) * B->vel(I[j]) / R[j],
	  Ri  = sqrt(square(B->pos(I[j])[0])+square(B->pos(I[j])[1])),
	  vpi = ami[2] / Ri,
	  vti = (B->pos(I[j])[2] * (B->pos(I[j])[0]*B->vel(I[j])[0]+
				    B->pos(I[j])[1]*B->vel(I[j])[1]) / Ri
		 - Ri * B->vel(I[j])[2])/R[j];
	M    += mi;
	Mvr  += mi*vri;
	Mvrq += mi*vri*vri;
	Mvrc += mi*vri*vri*vri;
	Mvrf += mi*vri*vri*vri*vri;
	Mvp  += mi*vpi;
	Mvpq += mi*vpi*vpi;
	Mvpc += mi*vpi*vpi*vpi;
	Mvpf += mi*vpi*vpi*vpi*vpi;
	Mvt  += mi*vti;
	Mvtq += mi*vti*vti;
	Mvtc += mi*vti*vti*vti;
	Mvtf += mi*vti*vti*vti*vti;
	Mam  += mi*ami;
      }
      rh[i] = iFPit*M/(i? cube(R[ir[i+1]])-cube(R[ir[i-1]]) :
		          cube(R[ir[i+1]]) );
      vr[i] = Mvr/M;
      vt[i] = Mvt/M;
      vp[i] = Mvp/M;
      am[i] = Mam/M;
      sr[i] = sqrt(M*Mvrq-Mvr*Mvr)/M;
      st[i] = sqrt(M*Mvtq-Mvt*Mvt)/M;
      sp[i] = sqrt(M*Mvpq-Mvp*Mvp)/M;
      skr[i]= (M*M*Mvrc-3*M*Mvr*Mvrq+2*Mvr*Mvr*Mvr)/(M*M*M*sr[i]*sr[i]*sr[i]);
      skt[i]= (M*M*Mvtc-3*M*Mvt*Mvtq+2*Mvt*Mvt*Mvt)/(M*M*M*st[i]*st[i]*st[i]);
      skp[i]= (M*M*Mvpc-3*M*Mvp*Mvpq+2*Mvp*Mvp*Mvp)/(M*M*M*sp[i]*sp[i]*sp[i]);
      kur[i]= (M*M*M*Mvrf-4*M*M*Mvrc*Mvr+6*M*Mvrq*Mvr*Mvr-3*pow(Mvr,4))/
	      (M*M*M*M*sr[i]*sr[i]*sr[i]*sr[i]);
      kut[i]= (M*M*M*Mvtf-4*M*M*Mvtc*Mvt+6*M*Mvtq*Mvt*Mvt-3*pow(Mvt,4))/
	      (M*M*M*M*st[i]*st[i]*st[i]*st[i]);
      kup[i]= (M*M*M*Mvpf-4*M*M*Mvpc*Mvp+6*M*Mvpq*Mvp*Mvp-3*pow(Mvp,4))/
	      (M*M*M*M*sp[i]*sp[i]*sp[i]*sp[i]);
    }
  } else {
    for(int i=0; i!=n; ++i) {                      // LOOP bins                 
      register double M = 0.;
      for(j=i? ir[i-1]:0; j<=ir[i+1]; ++j)         // LOOP bodies in bin        
	M  += (i!=0 && j==ir[i-1]  ||  j==ir[i+1]) ?
	       0.5 * B->mass(I[j]) : B->mass(I[j]);
      rh[i] = iFPit*M/(i? cube(R[ir[i+1]])-cube(R[ir[i-1]]) :
		          cube(R[ir[i+1]]) );
    }
  }
  falcON_DEL_A(ir);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::DistributionProfile                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  class Minus {
    const double* T;
  public:
    Minus(const double*t) : T(t) {}
    double operator[] (int i) const { return -T[i]; }
  };
  class Mass {
    const double* M;
  public:
    Mass(const double*m) : M(m) {}
    double operator[] (int i) const { return M? M[i] : 1.; }
  };
}
//------------------------------------------------------------------------------
DistributionProfile::DistributionProfile(int          N,
					 const double*F,
					 const double*__M,
					 int          __nmin,
					 double       __dmax) :
  kmin(__nmin/2), dmax(0.5*__dmax)
{
  // 1. sort F of all bodies in DESCENDING order
  int *I = new int[N];
  HeapIndex<double,Minus>(Minus(F),N,I);
  // 2. step through the bodies to make bins in F
  n = N/kmin;                                      // maximum possible #bins    
  int m=0;                                         // index of bin              
  int*r = new int[n];                              // rank of smallest F  in bin
  const double fac=exp(-dmax);                     // F factor between bins     
  r[0] = kmin;                                     // first rank = kmin         
  for(m=1; ; ++m) {                                // LOOP F bins               
    int j = r[m-1];                                //   begin of bin            
    const double fm = F[I[j]] * fac;               //   minimum F in bin        
    j    += kmin;                                  //   minimum rank in bin     
    while(F[I[j++]] > fm && j < N);                //   increase rank to F=Fmin 
    if(j + kmin >= N) j = N;                       //   enough left for next?   
    r[m] = min(N-1, j-1);                          //   highest rank in bin     
    if(r[m] == N-1) break;                         //   all bodies used? DONE!  
  }                                                // END LOOP                  
  n = m-1;                                         // this is the true #bins    
  // 3. establish tables with median F and cumulative masses at them
  Mass M(__M);                                     // mass_i or 1               
  fm = falcON_NEW(double,n);
  Mf = falcON_NEW(double,n);
  vf = falcON_NEW(double,n);
  double mcum=0.;                                  // cumulate mass             
  int j=0;                                         // body index                
  for(int i=0; i!=n; ++i) {                        // LOOP bins                 
    // we combine two adjacent small bins
    register int 
      s  = i==0? r[i+1]+1 : r[i+1]+r[i-1],         //   lower + upper rank      
      sh = s/2;                                    //   mean rank -> median F   
    // cumulate the mass up to the median F of combined bin
    for(; j!=sh; ++j)
      mcum += M[I[j]];
    // set median F and cumulated mass; be careful about even/odd sized bins
    if(s&1) {
      fm[i] = 0.5*(F[I[j]]+F[I[j+1]]);
      Mf[i] = mcum +     M[I[j]];
    } else {
      fm[i] = F[I[j]];
      Mf[i] = mcum + 0.5*M[I[j]];
    }
    // compute v(f) = (1/f) * dM/dF;
    register double dM = (i==0? M[I[0]] : 0.5*M[I[r[i-1]]]) + 0.5*M[I[r[i+1]]];
    for(int k=i? r[i-1]+1:1; k<r[i+1]; ++k) dM += M[I[k]];
    const double dF = (i==0? F[I[0]] : F[I[r[i-1]]]) - F[I[r[i+1]]];
    vf[i] = dM / (dF * fm[i]);
  }
  for(; j!=N; ++j) mcum += M[I[j]];
  Mt = mcum;                                       // total mass                
  Vf = falcON_NEW(double,n);
  Vf[0] = 0.5* vf[0] * (F[I[0]]-fm[1]);
  for(int i=1; i!=n; ++i)
    Vf[i] = Vf[i-1] + 0.5*(vf[i]+vf[i-1]) * abs(fm[i-1]-fm[i]);
  delete[] r;
  delete[] I;
}
////////////////////////////////////////////////////////////////////////////////

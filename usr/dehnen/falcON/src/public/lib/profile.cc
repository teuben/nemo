// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/profile.cc                                               
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2000-2008                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2000-2008 Walter Dehnen                                        
//                                                                              
// This program is free software; you can redistribute it and/or modify         
// it under the terms of the GNU General Public License as published by         
// the Free Software Foundation; either version 2 of the License, or (at        
// your option) any later version.                                              
//                                                                              
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//                                                                              
// You should have received a copy of the GNU General Public License            
// along with this program; if not, write to the Free Software                  
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#include <public/profile.h>
#include <utils/numerics.h>
#include <cmath>
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
  typedef tupel<2,real>   vect2;
  typedef tupel<2,double> vect2_d;
  vect centre;
  vect project;
  real radius        (body const&b) { return abs(pos(b)); }
  real radius_centred(body const&b) { return dist(centre,pos(b)); }
  real Radius        (body const&b) { return abs(pos(b)^project); }
  real Radius_centred(body const&b) { return abs((pos(b)-centre)^project); }
  //----------------------------------------------------------------------------
  template<typename X> inline
  void add_outer_product3(X p[3][3], vect const&x, real const&m) {
    p[0][0] += m * x[0] * x[0];
    p[0][1] += m * x[0] * x[1];
    p[0][2] += m * x[0] * x[2];
    p[1][1] += m * x[1] * x[1];
    p[1][2] += m * x[1] * x[2];
    p[2][2] += m * x[2] * x[2];
  }
  //----------------------------------------------------------------------------
  template<typename X> inline
  void symmetrize3(X p[3][3]) {
    p[1][0] = p[0][1];
    p[2][0] = p[0][2];
    p[2][1] = p[1][2];
  }
  //----------------------------------------------------------------------------
  template<typename X> inline
  void add_outer_product2(X p[2][2], vect2 const&x, real const&m) {
    p[0][0] += m * x[0] * x[0];
    p[0][1] += m * x[0] * x[1];
    p[1][1] += m * x[1] * x[1];
  }
  //----------------------------------------------------------------------------
  template<typename X> inline
  void symmetrize2(X p[2][2]) {
    p[1][0] = p[0][1];
  }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::spherical_profile                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
spherical_profile::spherical_profile(const bodies*B,
				     unsigned     __nmin,
				     double       __dmax,
				     const vect  *x0,
				     const vect  *v0)
  falcON_THROWING
: kmin(__nmin/2), dmax(0.5*__dmax), rr(0), mr(0), rh(0), ps(0), vr(0), sr(0),
  st(0), sp(0), ca(0), ba(0), am(0), vp(0), xa(0), xi(0)
{
  if(!B->have_all(fieldset(fieldset::m|fieldset::x)))
    falcON_THROW("spherical_profile: need \"mx\" to estimate density\n");
  // 1. sort radii of all bodies                                                
  Array<real>  R;
  Array<index> I;
  if(x0) centre = *x0;
  B->sorted(I,R,x0? radius_centred : radius);
  nb = I.size();
  if(nb < kmin)
    falcON_THROW("spherical_profile: too few (%d) bodies in_subset()\n",nb);
  // 2. step through the bodies to make radial bins                             
  n = nb/kmin;                                     // maximum possible #bins    
  int m=0;                                         // index of bin              
  int*ir= falcON_NEW(int,n);                       // rank of largest radius    
  const double fac=exp(dmax);                      // radius factor between bins
  ir[0] = kmin;                                    // first rank = kmin         
  for(m=1; ; ++m) {                                // LOOP radial bins          
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
  // 4. for the density etc we need to sum over the bin contents ...            
  rh = falcON_NEW(double,n);                       // density in bin            
  ca = falcON_NEW(double,n);
  ba = falcON_NEW(double,n);
  xa = falcON_NEW(vect_d,n);
  xi = falcON_NEW(vect_d,n);
  if(B->have(fieldbit::v)) {
    vr = falcON_NEW(double,n);
    vp = falcON_NEW(vect_d,n);
    am = falcON_NEW(vect_d,n);
    sr = falcON_NEW(double,n);
    st = falcON_NEW(double,n);
    sp = falcON_NEW(double,n);
    for(int i=0; i!=n; ++i) {                      // LOOP bins                 
      double M  (0.), Mvr(0.), Mvrq(0.);
      vect_d Mvp(0.), Mam(0.);
      double Mxx[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
      for(j=i? ir[i-1]:0; j<=ir[i+1]; ++j) {       // 1st LOOP of bodies in bin 
	double mi  = (i!=0 && j==ir[i-1]  ||  j==ir[i+1]) ?
	             0.5 * B->mass(I[j]) : B->mass(I[j]);
	vect_d ri  = x0? B->pos(I[j])-(*x0) : B->pos(I[j]);
	vect_d vi  = v0? B->vel(I[j])-(*v0) : B->vel(I[j]);
	if(R[j]==zero) {
	  double vri = abs(vi);
	  M   += mi;
	  Mvr += mi*vri;
	  Mvrq+= mi*vri*vri;
	} else {
	  vect_d er  = ri/R[j];
	  double vri = er*vi;
	  M   += mi;
	  Mvr += mi*vri;
	  Mvrq+= mi*vri*vri;
	  Mam += mi*(ri^vi);
	  Mvp += mi*(er^vi);
	}
	add_outer_product3(Mxx,ri,mi);
      }
      rh[i] = iFPit*M/(i? cube(R[ir[i+1]])-cube(R[ir[i-1]]) : cube(R[ir[i+1]]));
      vr[i] = Mvr/M;
      vp[i] = Mvp/M;
      am[i] = Mam/M;
      sr[i] = sqrt(M*Mvrq-Mvr*Mvr)/M;
      symmetrize3(Mxx);
      double IV[3][3], ID[3];
      int    IR;
      EigenSymJacobiSorted<3,double>(Mxx,IV,ID,IR);
      ca[i] = sqrt(ID[2]/ID[0]);
      ba[i] = sqrt(ID[1]/ID[0]);
      xa[i] = vect_d(IV[0][0],IV[1][0],IV[2][0]);
      xi[i] = vect_d(IV[0][2],IV[1][2],IV[2][2]);
      vect_d erot = norm(Mvp)>0.? normalized(Mvp) : vect_d(0.,0.,1.);
      double Mvpq(0.), Mvtq(0.);
      for(j=i? ir[i-1]:0; j<=ir[i+1]; ++j) {       // 2nd LOOP of bodies in bin 
	double mi  = (i!=0 && j==ir[i-1]  ||  j==ir[i+1]) ?
	             0.5 * B->mass(I[j]) : B->mass(I[j]);
	vect_d ri  = x0? B->pos(I[j])-(*x0) : B->pos(I[j]);
	vect_d vi  = v0? B->vel(I[j])-(*v0) : B->vel(I[j]);
	if(R[j]>zero) {
	  vect_d er  = ri/R[j];
	  vect_d ep  = normalized(er^erot);
	  vect_d et  = normalized(er^ep);
	  Mvpq += mi * square(vi*ep);
	  Mvtq += mi * square(vi*et);
	}
      }
      st[i] = sqrt(Mvtq/M);
      sp[i] = sqrt(M*Mvpq-Mvp*Mvp)/M;
    }
  } else {
    for(int i=0; i!=n; ++i) {                      // LOOP bins                 
      double M  (0.);
      double Mxx[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
      for(j=i? ir[i-1]:0; j<=ir[i+1]; ++j) {       // 1st LOOP of bodies in bin 
	double mi  = (i!=0 && j==ir[i-1]  ||  j==ir[i+1]) ?
	             0.5 * B->mass(I[j]) : B->mass(I[j]);
	vect_d ri  = x0? B->pos(I[j])-(*x0) : B->pos(I[j]);
	M   += mi;
	add_outer_product3(Mxx,ri,mi);
      }
      rh[i] = iFPit*M/(i? cube(R[ir[i+1]])-cube(R[ir[i-1]]) : cube(R[ir[i+1]]));
      symmetrize3(Mxx);
      double IV[3][3], ID[3];
      int    IR;
      EigenSymJacobiSorted<3,double>(Mxx,IV,ID,IR);
      ca[i] = sqrt(ID[2]/ID[0]);
      ba[i] = sqrt(ID[1]/ID[0]);
      xa[i] = vect_d(IV[0][0],IV[1][0],IV[2][0]);
      xi[i] = vect_d(IV[0][2],IV[1][2],IV[2][2]);
    }
  }
  falcON_DEL_A(ir);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::projected_profile                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
projected_profile::projected_profile(const bodies*B,
				     vect const  &__proj,
				     unsigned     __nmin,
				     double       __dmax,
				     const vect  *x0,
				     const vect  *v0)
  falcON_THROWING :
kmin(__nmin/2), dmax(0.5*__dmax), elos(normalized(__proj)),
mr(0), rr(0), sd(0), vl(0), vr(0), sl(0), ba(0), ph(0), al(0)
{
  if(elos == 0.)
    falcON_THROW("projected_profile: projection along null vector\n");
  vect_d eY = elos == vect_d(1.,0.,0.) ?
    normalized(elos ^ vect_d(0.,0.,1.)) : normalized(elos ^ vect_d(1.,0.,0.));
  vect_d eX = eY ^ elos;
  if(!B->have_all(fieldset(fieldset::m|fieldset::x)))
    falcON_THROW("projected_profile: need \"mx\" to estimate density\n");
  // 1. sort radii of all bodies                                                
  Array<real>  R;
  Array<index> I;
  if(x0) centre = *x0;
  project = elos;
  B->sorted(I, R, x0? Radius_centred : Radius);
  nb = I.size();
  if(nb < kmin)
    falcON_THROW("projected_profile: too few (%d) bodies in_subset()\n",nb);
  // 2. step through the bodies to make radial bins                             
  n = nb/kmin;                                     // maximum possible #bins    
  int m=0;                                         // index of bin              
  int*ir= falcON_NEW(int,n);                       // rank of largest radius    
  const double fac=exp(dmax);                      // radius factor between bins
  ir[0] = kmin;                                    // first rank = kmin         
  for(m=1; ; ++m) {                                // LOOP radial bins          
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
  // 4. for the rest we need to sum over the bin contents ...                   
  sd = falcON_NEW(double,n);
  ba = falcON_NEW(double,n);
  ph = falcON_NEW(double,n);
  if(B->have(fieldbit::v)) {
    vl = falcON_NEW(double,n);
    vr = falcON_NEW(double,n);
    sl = falcON_NEW(double,n);
    al = falcON_NEW(double,n);
    for(int i=0; i!=n; ++i) {                      // LOOP bins                 
      double  M  (0.), Mvl(0.), Mvlq(0.);
      vect2_d Mvr(0.);
      double  Mxx[2][2] = {{0.,0.}, {0.,0.}};
      for(j=i? ir[i-1]:0; j<=ir[i+1]; ++j) {       // LOOP of bodies in bin 
	double  mi = (i!=0 && j==ir[i-1]  ||  j==ir[i+1]) ?
	             0.5 * B->mass(I[j]) : B->mass(I[j]);
	vect_d  ri = x0? B->pos(I[j])-(*x0) : B->pos(I[j]);
	vect2_d xi = vect2_d(eX*ri, eY*ri);
	vect_d  vi = v0? B->vel(I[j])-(*v0) : B->vel(I[j]);
	double  vl = elos * vi;
	vect2_d ei = normalized(xi);
	M    += mi;
	Mvl  += mi * vl;
	Mvlq += mi * vl*vl;
	Mvr  += xi * (mi*vl);
	add_outer_product2(Mxx,xi,mi);
      }
      sd[i] = M / (Pi * i? square(R[ir[i+1]]) - square(R[ir[i-1]])
		        :  square(R[ir[i+1]]) );
      vl[i] = Mvl/M;
      vr[i] = abs(Mvr)/M;
      sl[i] = sqrt(M*Mvlq-Mvl*Mvl)/M;
      al[i] = atan2(Mvr[1],Mvr[0]);
      symmetrize2(Mxx);
      double apb = Mxx[0][0]+Mxx[1][1];
      double amb = Mxx[0][0]-Mxx[1][1];
      double det = sqrt(square(amb)+4*square(Mxx[0][1]));
      double l1  = 0.5*(apb+det);
      ba[i] = (apb - det) / (apb + det);
      ph[i] = acos(Mxx[0][1] / sqrt(square(l1-Mxx[0][0])+square(Mxx[0][1])));
    }
  } else {
    for(int i=0; i!=n; ++i) {                      // LOOP bins                 
      double  M  (0.);
      double  Mxx[2][2] = {{0.,0.}, {0.,0.}};
      for(j=i? ir[i-1]:0; j<=ir[i+1]; ++j) {       // LOOP of bodies in bin 
	double  mi = (i!=0 && j==ir[i-1]  ||  j==ir[i+1]) ?
	             0.5 * B->mass(I[j]) : B->mass(I[j]);
	vect_d  ri = x0? B->pos(I[j])-(*x0) : B->pos(I[j]);
	vect2_d xi = vect2_d(eX*ri, eY*ri);
	M    += mi;
	add_outer_product2(Mxx,xi,mi);
      }
      sd[i] = M / (Pi * i? square(R[ir[i+1]]) - square(R[ir[i-1]])
		        :  square(R[ir[i+1]]) );
      symmetrize2(Mxx);
      double ab  = Mxx[0][0]-Mxx[1][1];
      double det = sqrt(square(ab)+4*square(Mxx[0][1]));
      double l1  = 0.5*(ab+det);
      ba[i] = (ab - det) / (ab + det);
      ph[i] = acos(Mxx[0][1] / sqrt(square(l1-Mxx[0][0])+square(Mxx[0][1])));
    }
  }
  falcON_DEL_A(ir);
}
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_PROPER
#  include <proper/profile.cc>
#endif

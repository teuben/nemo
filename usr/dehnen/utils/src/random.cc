// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    src/random.cc
///
/// \author  Walter Dehnen
/// \author  Paul McMillan
///
/// \date    1994-2005, 2008-2011
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 1994-2005, 2008,2009,2011  Walter Dehnen
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
#include <random.h>
#include <numerics.h>
#include <exception.h>

#ifdef __INTEL_COMPILER
#pragma warning (disable:981) /* "operands are evaluated in unspecified order" */
#pragma warning (disable:1418) /* intel: "this warning can safely be ignored" */
#endif

using namespace WDutils;
////////////////////////////////////////////////////////////////////////////////
//
// class WDutils::Random3
//
// as in Numerical Recipes
//
////////////////////////////////////////////////////////////////////////////////
namespace {
  const long mbig = 1000000000, mseed = 161803398;
}
//------------------------------------------------------------------------------
Random3::Random3(long idum) : inext(0), inextp(31)
{
  long  mj,mk;
  int   i,ii,k;
  mj     = mseed - (idum<0 ? -idum : idum);
  mj    %= mbig;
  ma[55] = mj;
  mk     = 1;
  for(i=1; i<=54; i++) {
    ii     = (21*i) % 55;
    ma[ii] = mk;
    mk     = mj-mk;
    mj     = ma[ii];
    if(mk<0) mk += mbig;
  }
  for(k=1; k<=4; k++)
    for(i=1; i<=55; i++) {
      ma[i] -= ma[1+(i+30) % 55];
      if(ma[i]<0) ma[i] += mbig;
    }
}
//------------------------------------------------------------------------------
double Random3::RandomDouble() const
{
  static double fac = 1./double(mbig);
  long   mj;
  double r;
  do {
    if(++inext  >= 56) inext  = 1;
    if(++inextp >= 56) inextp = 1;
    mj = ma[inext] - ma[inextp];
    while(mj<0) mj+=mbig;
    ma[inext] = mj;
    r=mj*fac;
  } while (r<0. || r>1.);
  return r;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class WDutils::Sobol                                                       //
//                                                                            //
// with this installation not more than 52 objects (Sobol sequences) can be   //
// created simultaneously.                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  static unsigned sobol_setb  = 30;
  const unsigned  sobol_MO = 52;
  static char sobol_f[sobol_MO]
  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  const unsigned sobol_d[sobol_MO]
  = {1,2,3,3,4,4,5,5,5,5,5,5,6,6,6,6,6,6,
     7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
     8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8};
  const unsigned sobol_p[sobol_MO]
  = {0,1,1,2,1,4,2,4,7,11,13,14,1,13,16,19,22,25,
     1,4,7,8,14,19,21,28,31,32,37,41,42,50,55,56,59,62,
     14,21,22,38,47,49,50,52,56,67,70,84,97,103,115,122};
}
//------------------------------------------------------------------------------
void Sobol::set_bits(const unsigned BITS)
{
  sobol_setb = (BITS<=0) ? 30 : BITS;
}
//------------------------------------------------------------------------------
Sobol::Sobol(int ACTL, unsigned BITS)
{
  // set degree and polynomial from the static tables. It's important, that     
  // no other object with the same degree and polynomial is currently existing, 
  // because its quasi random numbers would be the same as ours. If ACTL on     
  // input is set >=0 we don't care and construct the sequence No. ACTL.        
  if(ACTL>=0 && ACTL < int(sobol_MO))
    actl=unsigned(ACTL);
  else {
    for(actl=0; actl!=sobol_MO && sobol_f[actl]; ++actl) {}
    if(actl>=sobol_MO)
      WDutils_Error("in Sobol::Sobol(): trying to create the 53th object");
  }
  sobol_f[actl] += 1;
  if(BITS==0)
    bits = sobol_setb;
  else if((bits=BITS)<10)
    WDutils_Warning("in Sobol::Sobol(): "
		    "creating object with less than 10 bits");
  in  = 0;
  ix  = 0;
  fac = 1./(1L<<bits);
  unsigned degs = sobol_d[actl];
  unsigned poly = sobol_p[actl];
  // seed initial Mi; i=1,...,degs                                              
  // these must be odd integer numbers less than 2^i.                           
  // Finally the direction numbers are Vi = 2^(bits-i) * Mi                     
  unsigned i,i2,ip,l;
  unsigned long vi;
  v_allocated = WDutils_NEW(unsigned long,bits);
  v = v_allocated - 1;
  for(i=1,i2=2; i<=degs; i++,i2<<=1) {
    if(i2<=poly) 
      vi = 1;
    else {
      vi = unsigned(i2-poly);
      if(!(vi&1)) vi-= 1;
    }
    if(bits>i)
      v[i] = vi << (bits-i);
  }
  // now use the recurrence (Press et al. 1992, eq. 7.7.2) to create            
  // the remaining direction numbers. With Vi = 2^(bits-i) Mi it reads          
  // V[i] = (a[1]*V[i-1]) XOR (a[2]*V[i-2]) XOR ... XOR (a[q-1]*V[i-q+1])       
  //   XOR ( V[i-q] XOR V[i-q]/2^q )                                            
  for(i=degs+1; i<=bits; i++) {
    ip = poly;
    vi = v[i-degs];
    vi^= (vi>>degs);
    for(l=degs-1; l>0; l--) {
      if(ip&1) vi ^= v[i-l];
      ip >>= 1;
    }
    v[i] = vi;
  }
}
//------------------------------------------------------------------------------
Sobol::~Sobol() {
  WDutils_DEL_A(v_allocated);
  sobol_f[actl]  = 0;
}
//------------------------------------------------------------------------------
double Sobol::RandomDouble () const {
  unsigned long im=in++, j;
  for(j=1; j<=bits; j++) {
    if( !(im&1) ) break;
    im >>= 1;
  }
  if(j>bits)
    WDutils_Error("in Sobol::RandomDouble(): "
		  "trying to call more than 2^BITS times");
  ix^= v[j];
  return double(ix)*fac;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class WDutils::Normal                                                        
//                                                                              
////////////////////////////////////////////////////////////////////////////////
Normal::Normal(const RandomNumberGenerator*r1,
	       const RandomNumberGenerator*r2)
  : iset(0), R1(r1), R2(r2? r2:r1)
{
  if(R1==R2 && !R1->is_pseudo())
    WDutils_THROW("trying to construct \"Normal\" with a "
		  "single quasi-random number generator\n");
}
//------------------------------------------------------------------------------
double Normal::operator() () const
{
  if(iset) {
    iset = 0;
    return gset;
  } else {
    double v1,v2,rsq,fac;
    do {
      v1  = 2 * R1->RandomDouble() - 1.;
      v2  = 2 * R2->RandomDouble() - 1.;
      rsq = v1*v1 + v2*v2;
    } while (rsq>=1. || rsq <=0. );
    fac  = sqrt(-2.*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class WDutils::ExpDisk                                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
ExpDisk::ExpDisk(const RandomNumberGenerator*r, const double H) : 
  R(r), h(H), hi(1./h), hqi(hi*hi)
{
  Y[0] = P[0] = 0.;
  Y[N] = P[N] = 1.;
  for(unsigned i=1; i<N; ++i) {
    Y[i] = i/double(N);
    double y1 = 1.-Y[i];
    P[i] = 1.-exp(-Y[i]/y1)/y1;
  }
}
//------------------------------------------------------------------------------
double ExpDisk::ranvar() const
{
  double y, p=(*R)();
  while(p>=1.) p=(*R)();
  y = polev(p,P,Y,N1);
  return h*y/(1.-y);
}
//------------------------------------------------------------------------------
double ExpDisk::value(double x) const
{
  return x>=0? hqi*x*exp(-hi*x) : 0.;
}
//------------------------------------------------------------------------------
double ExpDisk::radius(double p) const
{
  double y;
  y = polev(p,P,Y,N1);
  return h*y/(1.-y);
}

//////////////////////////////////////////////////////////////////////////////
//
// class WDutils::PowerLawDist
//
////////////////////////////////////////////////////////////////////////////////
PowerLawDist::PowerLawDist(const RandomNumberGenerator*rng, double alpha,
			   double x1, double x2)
  : R(rng), xmin(x1), xmax(x2),
    p(alpha), p1(1+alpha), ip1(1/p1), islog(abs(p1)<1.e-14), 
    ranfc(islog? std::log(xmax/xmin) : (std::pow(xmax/xmin,p1)-1)),
    pnorm(islog? (1/ranfc) : (p1/(std::pow(xmax,p1)-std::pow(xmin,p1))))
{
  if(p1<=1.e14? xmin<=0 : xmin<0)
    WDutils_THROW("PowerLawDist: xmin=%g is too small\n",xmin);
  if(xmin>=xmax)
    WDutils_THROW("PowerLawDist: xmin=%g > xmax=%g\n",xmin,xmax);
}
//
double PowerLawDist::ranvar() const
{
  double ran = ranfc*(*R)();
  if(islog) ran = std::exp(ran);
  else      ran = std::pow(ran+1,ip1);
  return xmin*ran;
}
////////////////////////////////////////////////////////////////////////////////


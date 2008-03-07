// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    src/random.cc                                                      
///                                                                             
/// \author  Walter Dehnen                                                      
/// \author  Paul McMillan                                                      
///                                                                             
/// \date    1994-2005, 2008                                                    
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1994-2005, 2008  Walter Dehnen                                 
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
#include <random.h>
#include <numerics.h>
#include <exception.h>

using namespace WDutils;
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class WDutils::Random3                                                      //
//                                                                            //
// as in Numerical Recipes                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  const long   mbig  = 1000000000, mseed = 161803398;
  const double fac   = 1./double(mbig);
}
//------------------------------------------------------------------------------
Random3::Random3(long idum) : inext(0), inextp(31)
{
  register long  mj,mk;
  register int   i,ii,k;
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
double Random3::RandomDouble() const {
  register long   mj;
  register double r;
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
  static int  setb  = 30;
  const  int  MO    = 52;
  static char f[MO] ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  const  int  d[MO] ={1,2,3,3,4,4,5,5,5,5,5,5,6,6,6,6,6,6,
		      7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
		      8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8};
  const  int  p[MO] ={0,1,1,2,1,4,2,4,7,11,13,14,1,13,16,19,22,25,
		      1,4,7,8,14,19,21,28,31,32,37,41,42,50,55,56,59,62,
		      14,21,22,38,47,49,50,52,56,67,70,84,97,103,115,122};
}
//------------------------------------------------------------------------------
void set_bits(const int BITS)
{
  setb = (BITS<=0) ? 30 : BITS;
}
//------------------------------------------------------------------------------
Sobol::Sobol(int ACTL, int BITS)
{
  // set degree and polynomial from the static tables. It's important, that     
  // no other object with the same degree and polynomial is currently existing, 
  // because its quasi random numbers would be the same as ours. If ACTL on     
  // input is set >=0 we don't care and construct the sequence No. ACTL.        
  if(ACTL>=0 && ACTL < MO)
    actl=ACTL;
  else {
    for(actl=0; actl!=MO && f[actl]; ++actl) {}
    if(actl>=MO)
      WDutils_ErrorF("trying to create the 53th object","Sobol::Sobol()");
  }
  f[actl] += 1;
  if(BITS==0)
    bits = setb;
  else if((bits=BITS)<10)
    WDutils_WarningF("creating object with less than 10 bits","Sobol::Sobol()");
  in       = 0;
  ix       = 0;
  fac      = 1./(1L<<bits);
  int degs = d[actl];
  int poly = p[actl];
  // seed initial Mi; i=1,...,degs                                              
  // these must be odd integer numbers less than 2^i.                           
  // Finally the direction numbers are Vi = 2^(bits-i) * Mi                     
  register int i,i2,ip,l;
  register unsigned long vi;
  v = WDutils_NEW(unsigned long,bits)-1;
  for(i=1,i2=2; i<=degs; i++,i2<<=1) {
    if(i2<=poly) 
      vi = 1;
    else {
      vi = i2 - poly;
      if(!(vi&1)) vi-= 1;
    }
    v[i] = vi << (bits-i);
  }
  // now use the recurrence (Press et al. 1992, eq. 7.7.2) to create            
  // the remaining direction numbers. With Vi = 2^(bits-i) Mi it reads          
  // V[i] = (a[1]*V[i-1]) XOR (a[2]*V[i-2]) XOR ... XOR (a[q-1]*V[i-q+1])       
  //   XOR ( V[i-q] XOR V[i-q]/2^q )                                            
  for(i=degs+1; i<=int(bits); i++) {
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
  WDutils_DEL_A(v+1);
  f[actl]  = 0;
}
//------------------------------------------------------------------------------
double Sobol::RandomDouble () const {
  register unsigned long im=in++, j;
  for(j=1; j<=bits; j++) {
    if( !(im&1) ) break;
    im >>= 1;
  }
  if(j>bits)
    WDutils_ErrorF("trying to call more than 2^BITS times",
		  "Sobol::RandomDouble()");
  ix^= v[j];
  return double(ix)*fac;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class WDutils::Normal                                                        
//                                                                              
////////////////////////////////////////////////////////////////////////////////
Normal::Normal(const RandomNumberGenerator*r1,
	       const RandomNumberGenerator*r2) :
  iset(0), R1(r1), R2(r2? r2:r1)
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
    register double v1,v2,rsq,fac;
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
ExpDisk::ExpDisk(RandomNumberGenerator* r, const double H) : 
  N(256), N1(N+1), R(r), h(H), hi(1./h), hqi(hi*hi)
{
  register int    i;
  register double y1;
  Y = WDutils_NEW(double,N1);
  P = WDutils_NEW(double,N1);
  Y[0] = P[0] = 0.;
  Y[N] = P[N] = 1.;
  for(i=1; i<N; i++) {
    Y[i] = i/double(N);
    y1   = 1.-Y[i];
    P[i] = 1.-exp(-Y[i]/y1)/y1;
  }
}
//------------------------------------------------------------------------------
ExpDisk::~ExpDisk()
{
  WDutils_DEL_A(Y);
  WDutils_DEL_A(P);
}
//------------------------------------------------------------------------------
double ExpDisk::ranvar() const
{
  register double y, p=(*R)();
  while(p>=1.) p=(*R)();
  y = polev(p,P,Y,N1,4);
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
  register double y;
  y = polev(p,P,Y,N1,4);
  return h*y/(1.-y);
}

////////////////////////////////////////////////////////////////////////////////


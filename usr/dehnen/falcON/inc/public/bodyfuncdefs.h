// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/bodyfuncdefs.h                                           
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2004-2006                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2004-2006 Walter Dehnen                                        
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
//                                                                              
// This file is to included in source code written, compiled and loaded at      
// run time by routines in bodyfunc.cc. There is no other sensible usage of     
// this file. It defines macros that are resolved as body properties.           
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#if !defined(body_func) && !defined(bodies_func)
#error "bodyfuncdefs.h must only be included in a bodyfunc or bodies_func application"
#endif

#if defined(BD_TEST)
#  include<utils/random.h>
#endif
namespace {
  template<typename T> struct make_bool {
    static bool cond(const T&__x) { return bool(__x); } };
  template<> struct make_bool<vect> {
    static bool cond(const vect&__x) { return __x!=zero; } };
  // cond(x): type conversion to bool
  template<typename T> inline
  bool cond(const T&__x) { return make_bool<T>::cond(__x); }
}

#if defined(BD_TEST)
// this switch is used to create functions that merely compute need and type.

#define real double
#define vect tupel<3,double>

namespace {

  fieldset need(fieldset::o);
  WDutils::Random3 RNG(1);

  inline int  __index() { return int(RNG()); }
  inline int  __key  () { need |= fieldset::k; return int(RNG()); }
  inline int  __ntot () { return int(RNG()); }
  inline real __vcom () { need |= fieldset::x;
                          need |= fieldset::v; return real(RNG()); }
  inline real __mass () { need |= fieldset::m; return real(RNG()); }
  inline real __pot  () { need |= fieldset::p; return real(RNG()); }
  inline real __aux  () { need |= fieldset::y; return real(RNG()); }
  inline real __eps  () { need |= fieldset::e; return real(RNG()); }
  inline vect __pos  () { need |= fieldset::x; return vect(RNG()); }
  inline vect __vel  () { need |= fieldset::v; return vect(RNG()); }
  inline vect __acc  () { need |= fieldset::a; return vect(RNG()); }
  inline vect __jrk  () { need |= fieldset::j; return vect(RNG()); }

#ifdef falcON_SPH

  inline real __size () { need |= fieldset::H; return real(RNG()); }
  inline real __snum () { need |= fieldset::N; return real(RNG()); }
  inline real __uin  () { need |= fieldset::U; return real(RNG()); }
  inline real __udot () { need |= fieldset::I; return real(RNG()); }
  inline real __udex () { need |= fieldset::E; return real(RNG()); }
  inline real __entr () { need |= fieldset::S; return real(RNG()); }
  inline real __srho () { need |= fieldset::R; return real(RNG()); }
  inline real __fact () { need |= fieldset::F; return real(RNG()); }
  inline real __csnd () { need |= fieldset::C; return real(RNG()); }

#endif

}

#define i     __index()
#define m     __mass()
#define pos   __pos()
#define vel   __vel()
#define pot   __pot()
#define acc   __acc()
#define jrk   __jrk()
#define aux   __aux()
#define key   __key()
#define eps   __eps()
#define vr    __vcom()
#define vt    __vcom()
#define vphi  __vcom()
#define N     __ntot()
#define M     __mass()
#define Pi    3.141592653589793238
#define vector(X) vect(real(X))

#ifdef falcON_SPH

#  define h   __size()
#  define n   __snum()
#  define U   __uin()
#  define Ui  __udot()
#  define Ue  __udex()
#  define S   __entr()
#  define rho __srho()
#  define f   __fact()
#  define c   __csnd()

#endif

#else // defined(BD_TEST)
// this gives the bodyfunc macros

namespace {

  typedef bodies::iterator body;

  inline int ntot(bodies const&B) {
    return B.N_bodies();
  }

  inline real mtot(bodies const&B) {
    static const bodies*_B = 0;
    static       double _M = 0.;
    if( &B != _B || _M == 0.) {
      _B = &B;
      _M = 0.;
      LoopAllBodies(_B,Bi) _M += double(mass(Bi));
    }
    return _M;
  }

  inline int  ntot(body const&b) { return ntot(*(b.my_bodies())); }
  inline real mtot(body const&b) { return mtot(*(b.my_bodies())); }

  inline real vrad(body const&b) {
    register real r = abs(pos(b));
    if(r==zero) return zero;
    return (pos(b)*vel(b))/r;
  }
  inline real vtan(body const&b) {
    register real r = abs(pos(b));
    if(r==zero) return zero;
    return abs(pos(b)^vel(b))/r;
  }
  inline real vphi(body const&b) {
    register real R = sqrt(square(pos(b)[0])+square(pos(b)[1]));
    if(R==zero) return zero;
    return (pos(b)[0]*vel(b)[1] - pos(b)[1]*vel(b)[0])/R;
  }
}

// global properties of all bodies
#if   defined(body_func)
# define N    ntot(b)                     /* number N of bodies              */
# define M    mtot(b)                     /* total mass of bodies            */
#elif defined(bodies_func)
# define N    ntot(B)                     /* number N of bodies              */
# define M    mtot(B)                     /* total mass of bodies            */
#endif
#define Pi    3.141592653589793238        /* Pi                              */
#define vector(X) vect(real(X))           /* type conversion to vector       */

// basic body properties
#define i     bodyindex(b)                /* index                           */
#define m     mass(b)                     /* mass                            */
#define pos   pos(b)                      /* position vector                 */
#define vel   vel(b)                      /* velocity vector                 */
#define pot   pot(b)                      /* gravitational potential         */
#define acc   acc(b)                      /* acceleration vector             */
#define jrk   jerk(b)                     /* jerk vector                     */
#define aux   aux(b)                      /* auxiliary float                 */
#define key   key(b)                      /* auxiliary integer               */
#define eps   eps(b)                      /* softening length                */

#ifdef falcON_SPH

#  define h   size(b)                     /* SPH: smoothing length           */
#  define n   snum(b)                     /* # SPH neihbours                 */
#  define U   uin(b)                      /* SPH: internal energy            */
#  define Ud  udot(b)                     /* SPH: (dU/dt)_total              */
#  define Ue  udex(b)                     /* SPH: (dU/dt)_external           */
#  define S   entr(b)                     /* SPH: entropy (function)         */
#  define rho srho(b)                     /* SPH: gas density                */
#  define f   fact(b)                     /* SPH: factor f_i                 */
#  define c   csnd(b)                     /* SPH: sound speed                */

#endif

// components of and compounds from basic body properties
#define vt    vtan(b)                     /* |tangential velocity|           */
#define vphi  vphi(b)                     /* azimuthal component of vel      */
#define vr    vrad(b)                     /* radial velocity                 */

#endif // defined(BD_TEST)

#define x     pos[0]                      /* x-position                      */
#define y     pos[1]                      /* y-position                      */
#define z     pos[2]                      /* z-position                      */
#define vx    vel[0]                      /* x-velocity                      */
#define vy    vel[1]                      /* y-velocity                      */
#define vz    vel[2]                      /* z-velocity                      */
#define ax    acc[0]                      /* x-acceleration                  */
#define ay    acc[1]                      /* y-acceleration                  */
#define az    acc[2]                      /* z-acceleration                  */
#define r     abs(pos)                    /* radius                          */
#define R     sqrt(x*x+y*y)               /* cylindrical radius              */
#define v     abs(vel)                    /* total velocity                  */
#define a     abs(acc)                    /* total acceleration              */
#define l     (pos^vel)                   /* specific angular momentum       */
#define lx    (y*vz-z*vy)                 /* x-component of ---              */
#define ly    (z*vx-x*vz)                 /* y-component of ---              */
#define lz    (x*vy-y*vx)                 /* z-component of ---              */
#define ltot  abs(l)                      /* total specific angular momentum */
#define jtot  abs(l)                      /* --- for NEMO compatibility      */
#define L     (m*l)                       /* angular momentum                */
#define Lx    (m*lx)                      /* x-component of ---              */
#define Ly    (m*ly)                      /* y-component of ---              */
#define Lz    (m*lz)                      /* z-component of ---              */
#define Ltot  abs(L)                      /* total angular momentum          */
#define E     (half*norm(vel)+phi)        /* total specific energy           */
#define etot  (half*norm(vel)+phi)        /* --- for NEMO compatibility      */
#define phi   pot                         /* potential for NEMO ---          */

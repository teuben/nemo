// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/bodyfuncdefs.h                                           
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2004-2008                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2004-2008 Walter Dehnen                                        
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

#if defined(BD_TEST)
  // dummy functions that aid computing need and type

  template<typename T> struct TypeShift { typedef T      type; };
  template<> struct TypeShift<flags>    { typedef int    type; };
  template<> struct TypeShift<unsigned> { typedef int    type; };
  template<> struct TypeShift<indx>     { typedef int    type; };
  template<> struct TypeShift<peanokey> { typedef int    type; };
  template<> struct TypeShift<real>     { typedef double type; };
  template<> struct TypeShift<vect>     { typedef tupel<3,double> type; };
  // TestType<BIT>::type is
  // int, double,tupel<3,double>  for integers, scalars, vectors
  template<int BIT> struct TestType:
  public TypeShift<typename field_traits<BIT>::type> {};

  fieldset need(fieldset::empty);  // to accumulate fields needed
  // provides return value for dummy functions,
  // preventing compiler from optimising dummy function calls away
  WDutils::Random3 RNG(1);

  // dummy functions which are allowed but not field functions
  inline int    bodyindex() {                      return int   (RNG()); }
  inline int    ntot     () {                      return int   (RNG()); }
  inline double mtot     () { need |= fieldset::m; return double(RNG()); }
  inline double vrad     () { need |= fieldset::phases; return double(RNG()); }
  inline double vtan     () { need |= fieldset::phases; return double(RNG()); }
  inline double vphi     () { need |= fieldset::phases; return double(RNG()); }
  inline bool   is_sph   () { return RNG() > 0.5; }
  inline bool   is_std   () { return RNG() > 0.5; }
  inline bool   is_sink  () { return RNG() > 0.5; }
  // dummy field functions
#define DEF_DUMMY(BIT,NAME)			\
  inline TestType<BIT>::type NAME()		\
  {						\
    need |= fieldset(fieldbit(BIT));		\
    return TestType<BIT>::type(RNG());		\
  }

  DEF_NAMED(DEF_DUMMY);
#undef DEF_DUMMY
  // define 'b' to be empty, such that 'pos(b)' --> 'pos()'
#define b
  // define 'B' to be empty, such that 'mtot(B)' --> 'mtot()'
#define B

#else // defined(BD_TEST)

  typedef bodies::iterator body;

  // functions which are not field functions
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

#endif // defined(BD_TEST)

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
#  define K   entr(b)                     /* SPH: entropy (function)         */
#  define rho srho(b)                     /* SPH: gas density                */
#  define f   fact(b)                     /* SPH: factor f_i                 */
#  define c   csnd(b)                     /* SPH: sound speed                */

#endif

// components of and compounds from basic body properties
#define vt    vtan(b)                     /* |tangential velocity|           */
#define vphi  vphi(b)                     /* azimuthal component of vel      */
#define vr    vrad(b)                     /* radial velocity                 */

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
#define k     key                         /* body key                        */
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
#define phi   pot                         /* potential for NEMO ---          */
#define E     (half*norm(vel)+phi)        /* total specific energy           */
#define etot  (half*norm(vel)+phi)        /* --- for NEMO compatibility      */

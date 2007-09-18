// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/simd.h                                                   
///                                                                             
/// \brief  defined fcev4, a 4-vector of floats aligned to 16 byte address      
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2002-2004                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2002-2004  Walter Dehnen, Paul McMillan                        
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
#ifndef falcON_included_simd_h
#define falcON_included_simd_h

#ifndef falcON_included_iostream
#  include <iostream>
#  define falcON_included_iostream
#endif

#ifndef falcON_included_cmath
#  include <cmath>
#  define falcON_included_cmath
#endif

#ifndef falcON_included_basic_h
#  include <public/basic.h>
#endif

#ifdef rsqrt          // NEMO has the funny idea to #define rsqrt               
#  undef rsqrt
#endif

////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // struct falcON::fvec4                                                       
  //                                                                            
  // for 4 floating point values that are (supposed to be) alligned to 16 bytes 
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  struct falcON__align16 fvec4 {
    float W,X,Y,Z;
    //--------------------------------------------------------------------------
    // constructors                                                             
    //--------------------------------------------------------------------------
    fvec4() {}
    fvec4(fvec4 const&v) : W(v.W), X(v.X), Y(v.Y), Z(v.Z) {}
    fvec4(int   const&t) : W(t),   X(t),   Y(t),   Z(t)   {}
    fvec4(float const&t) : W(t),   X(t),   Y(t),   Z(t)   {}
    fvec4(float const&w, float const&x) : W(w), X(x), Y(x), Z(x) {}
    fvec4(float const&w, float const&x, float const&y, float const&z)
      : W(w), X(x), Y(y), Z(z) {}
    fvec4(const float*const&v) : W(v[0]), X(v[1]), Y(v[2]), Z(v[3]) {}
    //--------------------------------------------------------------------------
    // assignment                                                               
    //--------------------------------------------------------------------------
    fvec4& operator=(fvec4 const&v)       { W=v.W;  X=v.X;  Y=v.Y;  Z=v.Z;
                                            return *this; }
    fvec4& operator=(float const&t)       { W=t;    X=t;    Y=t;    Z=t;
                                            return *this; }
    fvec4& operator=(const float*const&v) { W=v[0]; X=v[1]; Y=v[2]; Z=v[3];
                                            return *this; }
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
    fvec4& operator*=(fvec4 const&v)      { W*=v.W;  X*=v.X;  Y*=v.Y;  Z*=v.Z;
                                            return *this; }
    //--------------------------------------------------------------------------
    fvec4& operator/=(fvec4 const&v) {
#ifdef falcON_SSE
      __asm__ ("movaps %0,     %%xmm0\n\t"	\
               "divps  %1,     %%xmm0\n\t"	\
               "movaps %%xmm0, %0    \n\t"	\
               : "=m" (W)			\
	       : "m"  (v) );
#else
      W/=v.W;
      X/=v.X;
      Y/=v.Y;
      Z/=v.Z;
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& operator+=(fvec4 const&v)      { W+=v.W;  X+=v.X;  Y+=v.Y;  Z+=v.Z;
                                            return *this; }
    //--------------------------------------------------------------------------
    fvec4& operator-=(fvec4 const&v)      { W-=v.W;  X-=v.X;  Y-=v.Y;  Z-=v.Z;
                                            return *this; }
    //--------------------------------------------------------------------------
    fvec4& operator*=(float const&t)      { W*=t;    X*=t;    Y*=t;    Z*=t;
                                            return *this; }
    //--------------------------------------------------------------------------
    fvec4& operator/=(float const&t)      { register float x=1./t;
                                            return operator*=(x);  }
    //--------------------------------------------------------------------------
    fvec4& ass_sum (fvec4 const&v, fvec4 const&w) {
      W = v.W + w.W;
      X = v.X + w.X;
      Y = v.Y + w.Y;
      Z = v.Z + w.Z;
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& ass_times (fvec4 const&v, float const&t) {
      W = v.W * t;
      X = v.X * t;
      Y = v.Y * t;
      Z = v.Z * t;
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& add_times (fvec4 const&v, float const&t) {
      W += v.W * t;
      X += v.X * t;
      Y += v.Y * t;
      Z += v.Z * t;
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& sub_times (fvec4 const&v, float const&t) {
      W -= v.W * t;
      X -= v.X * t;
      Y -= v.Y * t;
      Z -= v.Z * t;
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& ass_times (fvec4 const&v, fvec4 const&w) {
      W = v.W * w.W;
      X = v.X * w.X;
      Y = v.Y * w.Y;
      Z = v.Z * w.Z;
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& add_times (fvec4 const&v, fvec4 const&w) {
      W += v.W * w.W;
      X += v.X * w.X;
      Y += v.Y * w.Y;
      Z += v.Z * w.Z;
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& sub_times (fvec4 const&v, fvec4 const&w) {
      W -= v.W * w.W;
      X -= v.X * w.X;
      Y -= v.Y * w.Y;
      Z -= v.Z * w.Z;
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& ass_times (fvec4 const&v, fvec4 const&w, float const&t) {
      W = v.W * w.W * t;
      X = v.X * w.X * t;
      Y = v.Y * w.Y * t;
      Z = v.Z * w.Z * t;
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& add_times (fvec4 const&v, fvec4 const&w, float const&t) {
      W += v.W * w.W * t;
      X += v.X * w.X * t;
      Y += v.Y * w.Y * t;
      Z += v.Z * w.Z * t;
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& sub_times (fvec4 const&v, fvec4 const&w, float const&t) {
      W -= v.W * w.W * t;
      X -= v.X * w.X * t;
      Y -= v.Y * w.Y * t;
      Z -= v.Z * w.Z * t;
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& negate() { W=-W; X=-X; Y=-Y; Z=-Z; return *this; }
    //--------------------------------------------------------------------------
    fvec4& rcp() {    // w,x,y,z  ->  1/w, 1/x, 1/y, 1/z                        
#ifdef falcON_SSE
      __asm__ ("movaps %0,     %%xmm0\n\t"	\
               "rcpps  %%xmm0, %%xmm1\n\t"	\
               "movaps %%xmm1, %0    \n\t"	\
               : "=X" (W));
#else
      W=1./W;
      X=1./X;
      Y=1./Y;
      Z=1./Z;
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& rcp (fvec4 const&v) {
#ifdef falcON_SSE
      __asm__ ("movaps %1,     %%xmm0\n\t"	\
               "rcpps  %%xmm0, %%xmm1\n\t"	\
               "movaps %%xmm1, %0    \n\t"	\
               : "=m" (W)			\
               : "m"  (v));
#else
      W=1./v.W;
      X=1./v.X;
      Y=1./v.Y;
      Z=1./v.Z;
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& sqrt() {   // w,x,y,z  ->  sqrt(w), sqrt(x), sqrt(y), sqrt(z)        
#ifdef falcON_SSE
      __asm__ ("movaps %0,     %%xmm0\n\t"	\
               "sqrtps %%xmm0, %%xmm1\n\t"	\
               "movaps %%xmm1, %0    \n\t"	\
               : "=m" (W));
#else
      W=std::sqrt(W);
      X=std::sqrt(X);
      Y=std::sqrt(Y);
      Z=std::sqrt(Z); 
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& sqrt (fvec4 const&v) {
#ifdef falcON_SSE
      __asm__ ("movaps %1,     %%xmm0\n\t"	\
               "sqrtps %%xmm0, %%xmm1\n\t"	\
               "movaps %%xmm1, %0    \n\t"	\
               : "=m" (W)			\
               : "m"  (v));
#else
      W=std::sqrt(v.W);
      X=std::sqrt(v.X);
      Y=std::sqrt(v.Y);
      Z=std::sqrt(v.Z);
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& rsqrt() {  // w,x,y,z  ->  1/sqrt(w), 1/sqrt(x), 1/sqrt(y), 1/sqrt(z)
#ifdef falcON_SSE
      __asm__ ("movaps  %0,     %%xmm0\n\t"	\
               "rsqrtps %%xmm0, %%xmm1\n\t"	\
               "movaps  %%xmm1, %0    \n\t"	\
               : "=m" (W));
#else
      W=1./std::sqrt(W);
      X=1./std::sqrt(X);
      Y=1./std::sqrt(Y);
      Z=1./std::sqrt(Z); 
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    fvec4& rsqrt (fvec4 const&v) {
#ifdef falcON_SSE
      __asm__ ("movaps  %1,     %%xmm0\n\t"	\
               "rsqrtps %%xmm0, %%xmm1\n\t"	\
               "movaps  %%xmm1, %0    \n\t"	\
               : "=m" (W)			\
               : "m"  (v));
#else
      W=1./std::sqrt(v.W);
      X=1./std::sqrt(v.X);
      Y=1./std::sqrt(v.Y);
      Z=1./std::sqrt(v.Z);
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    // binary operations (the use of which should be minimized)                 
    //--------------------------------------------------------------------------
    fvec4 operator+ (fvec4 const&v) const { register fvec4 y(*this);
                                            return y+=v; }
    //--------------------------------------------------------------------------
    fvec4 operator- (fvec4 const&v) const { register fvec4 y(*this);
                                            return y-=v; }
    //--------------------------------------------------------------------------
    fvec4 operator* (fvec4 const&v) const { register fvec4 y(*this);
                                            return y*=v; }
    //--------------------------------------------------------------------------
    fvec4 operator* (float const&t) const { register fvec4 y(*this);
                                            return y*=t; }
    //--------------------------------------------------------------------------
    fvec4 operator/ (float const&t) const { register fvec4 y(*this);
                                            return y*=1./t; }
    //--------------------------------------------------------------------------
    // unary minus (preferrably use negate())                                   
    //--------------------------------------------------------------------------
    fvec4 operator- () const { register fvec4 y(*this); return y.negate(); }
    //--------------------------------------------------------------------------
    // const methods                                                            
    //--------------------------------------------------------------------------
    float norm() const {
      return W*W + X*X + Y*Y + Z*Z;
    }
    //--------------------------------------------------------------------------
    float dot(fvec4 const&v) const { return W*v.W + X*v.X + Y*v.Y + Z*v.Z; }
    //--------------------------------------------------------------------------
    // friends versions of some methods                                         
    //--------------------------------------------------------------------------
    friend float norm(fvec4 const&v) { return v.norm(); }
    //--------------------------------------------------------------------------
    friend fvec4 rcp (fvec4 const&v) {
      register fvec4 w;
      return w.rcp(v);
    }
    //--------------------------------------------------------------------------
    friend fvec4 sqrt (fvec4 const&v) {
      register fvec4 w;
      return w.sqrt(v);
    }
    //--------------------------------------------------------------------------
    friend fvec4 rsqrt (fvec4 const&v) {
      register fvec4 w;
      return w.rsqrt(v);
    }
    //--------------------------------------------------------------------------
    // conversion and data access                                               
    //--------------------------------------------------------------------------
    operator float*        ()       { return (      float*)(this); }
    operator const float*  () const { return (const float*)(this); }
    //--------------------------------------------------------------------------
    float&       operator[] (int i)       { return ((      float*)(this))[i]; }
    float const& operator[] (int i) const { return ((const float*)(this))[i]; }
    //--------------------------------------------------------------------------
    // allocation of many fvec4                                                 
    //--------------------------------------------------------------------------
    void* operator new   [](        size_t n) { return malloc16(n); }
    void  operator delete[](void*q, size_t n) { return free16(q); }
    //--------------------------------------------------------------------------
    // output to std::ostream; there is no input from std::istream              
    //--------------------------------------------------------------------------
    friend std::ostream& operator<<(std::ostream&s, fvec4 const&x) {
      return s<<x.W<<' '<<x.X<<' '<<x.Y<<' '<<x.Z;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_simd_h

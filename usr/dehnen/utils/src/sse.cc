// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/src/sse.cc
///
/// \brief  implements sse.h
///
/// \author Walter Dehnen
///
/// \date   2009
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Walter Dehnen
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
#include <sse.h>
#include <cmath>
#include <cstring>

using namespace WDutils;

#ifdef __SSE__
# define __FloatTmp							\
  __m128 Tm = _mm_setzero_ps(); float tm(0.f);
# define __FloatRet							\
  __attribute__ ((aligned(16))) float tmp[4];				\
  _mm_store_ps(tmp,Tm);							\
  return tmp[0]+tmp[1]+tmp[2]+tmp[3]+tm;
# define __FloatTmp16							\
  __m128 Tm = _mm_setzero_ps();
# define __FloatRet16							\
  __attribute__ ((aligned(16))) float tmp[4];				\
  _mm_store_ps(tmp,Tm);							\
  return tmp[0]+tmp[1]+tmp[2]+tmp[3]
# define __FloatSSE(__NAME,__INIT,__WORKS,__WORKU,__WORKA)		\
  size_t i=size_t(f)&16;						\
  if(n<4 || i&4) {		      /* small or unaligned     */	\
    if(i&4)								\
      WDutils_Warning("SSE::%s(): "					\
		      "array of floats @ %p not 4-byte aligned: "	\
		      "unaligned (slow) SSE instructions forced\n",	\
		      __NAME,f);					\
    i=0;								\
    size_t n4=SSE::bottom<4>(n);					\
    __INIT;								\
    for(; i!=n4; i+=4) { __WORKU; }   /* unaligned blocks work  */	\
    for(; i!=n ;  ++i) { __WORKS; }   /* final few elements     */	\
  } else {								\
    i >>= 2;                          /* begin: shifted pointer */	\
    n  += i;                          /* end:   shifted pointer */	\
    f  -= i;                          /* shift pointer          */	\
    size_t i4=SSE::top<4>(i);         /* begin: aligned blocks  */	\
    size_t n4=SSE::bottom<4>(n);      /* end:   aligned blocks  */	\
    __INIT;								\
    for(; i!=i4;  ++i) { __WORKS; }   /* STD work @ begin       */	\
    for(; i!=n4; i+=4) { __WORKA; }   /* SSE work in blocks     */	\
    for(; i!=n ;  ++i) { __WORKS; }   /* STD work @ end         */	\
  }
# define __FloatSSE16(__INIT,__WORKS,__WORKA)				\
  __INIT;								\
  for(size_t i=0; i!=n; i+=4) { __WORKA; }
#else
# warning SSE instructions not supported: will implement simple code
# ifdef __GNUC__
# warning perhaps using the -msse compiler option helps
# endif
# define __FloatTmp							\
  float tm(0.f);
# define __FloatRet							\
  return tm;
# define __FloatTmp16 __FloatTmp
# define __FloatRet16 __FloatRet
# define __FloatSSE(__NAME,__INIT,__WORKS,__WORKU,__WORKA)		\
  for(size_t i=0; i!=n; ++i) { __WORKS; }
# define __FloatSSE16(__INIT,__WORKS,__WORKA)				\
  for(size_t i=0; i!=n; ++i) { __WORKS; }
#endif

#ifdef __SSE2__
# define __DoubleTmp							\
  __m128d Tm = _mm_setzero_pd(); double tm(0.0);
# define __DoubleRet							\
  __attribute__ ((aligned(16))) double tmp[2];				\
  _mm_store_pd(tmp,Tm);							\
  return tmp[0]+tmp[1]+tm;
# define __DoubleTmp16							\
  __m128d Tm = _mm_setzero_pd();
# define __DoubleRet16							\
  __attribute__ ((aligned(16))) double tmp[2];				\
  _mm_store_pd(tmp,Tm);							\
  return tmp[0]+tmp[1];
# define __DoubleSSE(__NAME,__INIT,__WORKS,__WORKU,__WORKA)		\
  size_t i=size_t(f)&16;						\
  if(n<2 || i&8) {		      /* small or unaligned     */	\
    if(i&8) 								\
      WDutils_Warning("SSE::%s(): "					\
		      "array of doubles @ %p not 8-byte aligned: "	\
		      "unaligned (slow) SSE instructions forced\n",	\
		      __NAME,f);					\
    i=0;								\
    size_t n2=SSE::bottom<2>(n); 					\
    __INIT;								\
    for(; i!=n2; i+=2) { __WORKU; }   /* unaligned blocks work  */	\
    for(; i!=n ;  ++i) { __WORKS; }   /* final few elements     */	\
  } else {								\
    i >>= 3;                          /* begin: shifted pointer */	\
    n  += i;                          /* end:   shifted pointer */	\
    f  -= i;                          /* shift pointer          */	\
    size_t i2=SSE::top<2>(i);         /* begin: aligned blocks  */	\
    size_t n2=SSE::bottom<2>(n);      /* end:   aligned blocks  */	\
    __INIT;								\
    for(; i!=i2;  ++i) { __WORKS; }   /* STD work @ begin       */	\
    for(; i!=n2; i+=2) { __WORKA; }   /* SSE work in blocks     */	\
    for(; i!=n ;  ++i) { __WORKS; }   /* STD work @ end         */	\
  }
# define __DoubleSSE16(__INIT,__WORKS,__WORKA)				\
  __INIT;								\
  for(size_t i=0; i!=n; i+=4) { __WORKA; }
#else
# warning SSE2 instructions not supported: will implement simple code
# ifdef __GNUC__
# warning perhaps using the -msse2 compiler option helps
# endif
# define __DoubleTmp							\
  double tm(0.0);
# define __DoubleRet							\
  return tm;
# define __DoubleTmp16 __DoubleTmp
# define __DoubleRet16 __DoubleRet
# define __DoubleSSE(__NAME,__INIT,__WORKS,__WORKU,__WORKA)		\
  for(size_t i=0; i!=n; ++i) { __WORKS; }
# define __DoubleSSE16(__INIT,__WORKS,__WORKA)	\
  for(size_t i=0; i!=n; ++i) { __WORKS; }
#endif

// SSE::Ass, SSE::Ass16
#define __Works f[i] = v;
#define __Init  __m128 V = _mm_set1_ps(v);
#define __Worku _mm_storeu_ps(f+i,V);
#define __Worka _mm_store_ps(f+i,V);
void SSE::Ass(float*f, size_t n, float v)
{ __FloatSSE("Ass",__Init,__Works,__Worku,__Worka); }
void SSE::Ass16(float*f, size_t n, float v)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128d V = _mm_set1_pd(v);
#define __Worku _mm_storeu_pd(f+i,V);
#define __Worka _mm_store_pd(f+i,V);
void SSE::Ass(double*f, size_t n, double v)
{ __DoubleSSE("Ass",__Init,__Works,__Worku,__Worka); }
void SSE::Ass16(double*f, size_t n, double v)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka
#undef  __Works
// SSE::Neg, SSE::Neg16
#define __Works f[i] =-f[i];
#define __Init  __m128 V = _mm_set1_ps(0.f);
#define __Worku _mm_storeu_ps(f+i,_mm_sub_ps(V,_mm_loadu_ps(f+i)));
#define __Worka _mm_store_ps(f+i,_mm_sub_ps(V,_mm_load_ps(f+i)));
void SSE::Neg(float*f, size_t n)
{ __FloatSSE("Ass",__Init,__Works,__Worku,__Worka); }
void SSE::Neg16(float*f, size_t n)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128d V = _mm_set1_pd(0.f);
#define __Worku _mm_storeu_pd(f+i,_mm_sub_pd(V,_mm_loadu_pd(f+i)));
#define __Worka _mm_store_pd(f+i,_mm_sub_pd(V,_mm_load_pd(f+i)));
void SSE::Neg(double*f, size_t n)
{ __DoubleSSE("Ass",__Init,__Works,__Worku,__Worka); }
void SSE::Neg16(double*f, size_t n)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka
#undef  __Works
// SSE::Add, SSE::Add16
#define __Works f[i] += v;
#define __Init  __m128 V = _mm_set1_ps(v);
#define __Worku _mm_storeu_ps(f+i,_mm_add_ps(_mm_loadu_ps(f+i),V));
#define __Worka _mm_store_ps(f+i,_mm_add_ps(_mm_load_ps(f+i),V));
void SSE::Add(float*f, size_t n, float v)
{ __FloatSSE("Add",__Init,__Works,__Worku,__Worka); }
void SSE::Add16(float*f, size_t n, float v)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128d V = _mm_set1_pd(v);
#define __Worku _mm_storeu_pd(f+i,_mm_add_pd(_mm_loadu_pd(f+i),V));
#define __Worka _mm_store_pd(f+i,_mm_add_pd(_mm_load_pd(f+i),V));
void SSE::Add(double*f, size_t n, double v)
{ __DoubleSSE("Add",__Init,__Works,__Worku,__Worka); }
void SSE::Add16(double*f, size_t n, double v)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka
#undef  __Works

#define __Init
#define __Works a[i] += b[i];
#define __Worka _mm_store_ps(a+i,_mm_add_ps(_mm_load_ps(a+i),_mm_load_ps(b+i)));
void SSE::Add16(float*a, size_t n, const float*b)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Worka
#define __Worka _mm_store_pd(a+i,_mm_add_pd(_mm_load_pd(a+i),_mm_load_pd(b+i)));
void SSE::Add16(double*a, size_t n, const double*b)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Works
#undef  __Worka

#define __Init  __m128 W = _mm_set1_ps(w);
#define __Works a[i] += w*b[i];
#define __Worka _mm_store_ps(a+i,_mm_add_ps(_mm_load_ps(a+i), \
					    _mm_mul_ps(W,_mm_load_ps(b+i))));
void SSE::Add16(float*a, size_t n, float w, const float*b)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worka
#define __Init  __m128d W = _mm_set1_pd(w);
#define __Worka _mm_store_pd(a+i,_mm_add_pd(_mm_load_pd(a+i), \
					    _mm_mul_pd(W,_mm_load_pd(b+i))));
void SSE::Add16(double*a, size_t n, double w, const double*b)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Works
#undef  __Worka
// SSE::Sum, SSE::Sum16
#define __Init
#define __Works tm += f[i];
#define __Worka Tm = _mm_add_ps(Tm,_mm_load_ps(f+i));
#define __Worku Tm = _mm_add_ps(Tm,_mm_loadu_ps(f+i));
float SSE::Sum(const float*f, size_t n)
{
  __FloatTmp;
  __FloatSSE("Sum",__Init,__Works,__Worku,__Worka);
  __FloatRet;
}
float SSE::Sum16(const float*f, size_t n)
{
  __FloatTmp16;
  __FloatSSE16(__Init,__Works,__Worka);
  __FloatRet16;
}
#undef  __Worku
#undef  __Worka
#define __Worka Tm = _mm_add_pd(Tm,_mm_load_pd(f+i));
#define __Worku Tm = _mm_add_pd(Tm,_mm_loadu_pd(f+i));
double SSE::Sum(const double*f, size_t n)
{
  __DoubleTmp;
  __DoubleSSE("Sum",__Init,__Works,__Worku,__Worka);
  __DoubleRet;
}
double SSE::Sum16(const double*f, size_t n)
{
  __DoubleTmp16;
  __DoubleSSE16(__Init,__Works,__Worka);
  __DoubleRet16;
}
#undef  __Init
#undef  __Worku
#undef  __Works
#undef  __Worka
// SSE::Dot16
#define __Init
#define __Works tm += a[i]*b[i];
#define __Worka Tm = _mm_add_ps(Tm,_mm_mul_ps(_mm_load_ps(a+i),		\
					      _mm_load_ps(b+i)));
float SSE::Dot16(const float*a, size_t n, const float*b)
{
  __FloatTmp16;
  __FloatSSE16(__Init,__Works,__Worka);
  __FloatRet16;
}
#undef  __Worka
#define __Worka Tm = _mm_add_pd(Tm,_mm_mul_pd(_mm_load_pd(a+i),		\
					      _mm_load_pd(b+i)));
double SSE::Dot16(const double*a, size_t n, const double*b)
{
  __DoubleTmp16;
  __DoubleSSE16(__Init,__Works,__Worka);
  __DoubleRet16;
}
#undef  __Init
#undef  __Works
#undef  __Worka
// SSE::Sub, SSE::Sub16
#define __Works f[i] -= v;
#define __Init  __m128 V = _mm_set1_ps(v);
#define __Worku _mm_storeu_ps(f+i,_mm_sub_ps(_mm_loadu_ps(f+i),V));
#define __Worka _mm_store_ps(f+i,_mm_sub_ps(_mm_load_ps(f+i),V));
void SSE::Sub(float*f, size_t n, float v)
{ __FloatSSE("Sub",__Init,__Works,__Worku,__Worka); }
void SSE::Sub16(float*f, size_t n, float v)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128d V = _mm_set1_pd(v);
#define __Worku _mm_storeu_pd(f+i,_mm_sub_pd(_mm_loadu_pd(f+i),V));
#define __Worka _mm_store_pd(f+i,_mm_sub_pd(_mm_load_pd(f+i),V));
void SSE::Sub(double*f, size_t n, double v)
{ __DoubleSSE("Sub",__Init,__Works,__Worku,__Worka); }
void SSE::Sub16(double*f, size_t n, double v)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka
#undef  __Works

#define __Init
#define __Works a[i] -= b[i];
#define __Worka _mm_store_ps(a+i,_mm_sub_ps(_mm_load_ps(a+i),_mm_load_ps(b+i)));
void SSE::Sub16(float*a, size_t n, const float*b)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Worka
#define __Worka _mm_store_pd(a+i,_mm_sub_pd(_mm_load_pd(a+i),_mm_load_pd(b+i)));
void SSE::Sub16(double*a, size_t n, const double*b)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Works
#undef  __Worka

#define __Init  __m128 W = _mm_set1_ps(w);
#define __Works a[i] -= w*b[i];
#define __Worka _mm_store_ps(a+i,_mm_sub_ps(_mm_load_ps(a+i), \
					    _mm_mul_ps(W,_mm_load_ps(b+i))));
void SSE::Sub16(float*a, size_t n, float w, const float*b)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worka
#define __Init  __m128d W = _mm_set1_pd(w);
#define __Worka _mm_store_pd(a+i,_mm_sub_pd(_mm_load_pd(a+i), \
					    _mm_mul_pd(W,_mm_load_pd(b+i))));
void SSE::Sub16(double*a, size_t n, double w, const double*b)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Works
#undef  __Worka
// SSE::Mul, SSE::Mul16
#define __Works f[i] *= v;
#define __Init  __m128 V = _mm_set1_ps(v);
#define __Worku _mm_storeu_ps(f+i,_mm_mul_ps(_mm_loadu_ps(f+i),V));
#define __Worka _mm_store_ps(f+i,_mm_mul_ps(_mm_load_ps(f+i),V));
void SSE::Mul(float*f, size_t n, float v)
{ __FloatSSE("Mul",__Init,__Works,__Worku,__Worka); }
void SSE::Mul16(float*f, size_t n, float v)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128d V = _mm_set1_pd(v);
#define __Worku _mm_storeu_pd(f+i,_mm_mul_pd(_mm_loadu_pd(f+i),V));
#define __Worka _mm_store_pd(f+i,_mm_mul_pd(_mm_load_pd(f+i),V));
void SSE::Mul(double*f, size_t n, double v)
{ __DoubleSSE("Mul",__Init,__Works,__Worku,__Worka); }
void SSE::Mul16(double*f, size_t n, double v)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka
#undef  __Works

#define __Init
#define __Works a[i] -= b[i];
#define __Worka _mm_store_ps(a+i,_mm_mul_ps(_mm_load_ps(a+i),_mm_load_ps(b+i)));
void SSE::Mul16(float*a, size_t n, const float*b)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Worka
#define __Worka _mm_store_pd(a+i,_mm_mul_pd(_mm_load_pd(a+i),_mm_load_pd(b+i)));
void SSE::Mul16(double*a, size_t n, const double*b)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Works
#undef  __Worka
// SSE::Div16
#define __Init
#define __Works a[i] -= b[i];
#define __Worka _mm_store_ps(a+i,_mm_div_ps(_mm_load_ps(a+i),_mm_load_ps(b+i)));
void SSE::Div16(float*a, size_t n, const float*b)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Worka
#define __Worka _mm_store_pd(a+i,_mm_div_pd(_mm_load_pd(a+i),_mm_load_pd(b+i)));
void SSE::Div16(double*a, size_t n, const double*b)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Works
#undef  __Worka
// SSE::Inv, SSE::Inv16
#define __Works f[i] = x/f[i];
#define __Init  __m128 X = _mm_set1_ps(x);
#define __Worku _mm_storeu_ps(f+i,_mm_div_ps(X,_mm_loadu_ps(f+i)));
#define __Worka _mm_store_ps(f+i,_mm_div_ps(X,_mm_load_ps(f+i)));
void SSE::Inv(float*f, size_t n, float x)
{ __FloatSSE("Inv",__Init,__Works,__Worku,__Worka); }
void SSE::Inv16(float*f, size_t n, float x)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128d X = _mm_set1_pd(x);
#define __Worku _mm_storeu_pd(f+i,_mm_div_pd(X,_mm_loadu_pd(f+i)));
#define __Worka _mm_store_pd(f+i,_mm_div_pd(X,_mm_load_pd(f+i)));
void SSE::Inv(double*f, size_t n, double x)
{ __DoubleSSE("Inv",__Init,__Works,__Worku,__Worka); }
void SSE::Inv16(double*f, size_t n, double x)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka
#undef  __Works

#define __Init  __m128 X = _mm_set1_ps(x);
#define __Works a[i] = x/b[i];
#define __Worka _mm_store_ps(a+i,_mm_div_ps(X,_mm_load_ps(b+i)));
void SSE::Inv16(float*a, size_t n, float x, const float*b)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worka
#define __Init  __m128d X = _mm_set1_pd(x);
#define __Worka _mm_store_pd(a+i,_mm_div_pd(X,_mm_load_pd(b+i)));
void SSE::Inv16(double*a, size_t n, double x, const double*b)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Works
#undef  __Worka
// SSE::Sqrt, SSE::Sqrt16
#define __Works f[i] = std::sqrt(f[i]);
#define __Init
#define __Worku _mm_storeu_ps(f+i,_mm_sqrt_ps(_mm_loadu_ps(f+i)));
#define __Worka _mm_store_ps(f+i,_mm_sqrt_ps(_mm_load_ps(f+i)));
void SSE::Sqrt(float*f, size_t n)
{ __FloatSSE("Sqrt",__Init,__Works,__Worku,__Worka); }
void SSE::Sqrt16(float*f, size_t n)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Worku
#undef  __Worka

#define __Worku _mm_storeu_pd(f+i,_mm_sqrt_pd(_mm_loadu_pd(f+i)));
#define __Worka _mm_store_pd(f+i,_mm_sqrt_pd(_mm_load_pd(f+i)));
void SSE::Sqrt(double*f, size_t n)
{ __DoubleSSE("Sqrt",__Init,__Works,__Worku,__Worka); }
void SSE::Sqrt16(double*f, size_t n)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka
#undef  __Works

#define __Init
#define __Works a[i] = std::sqrt(b[i]);
#define __Worka _mm_store_ps(a+i,_mm_sqrt_ps(_mm_load_ps(b+i)));
void SSE::Sqrt16(float*a, size_t n, const float*b)
{ __FloatSSE16(__Init,__Works,__Worka); }
#undef  __Worka
#define __Worka _mm_store_pd(a+i,_mm_sqrt_pd(_mm_load_pd(b+i)));
void SSE::Sqrt16(double*a, size_t n, const double*b)
{ __DoubleSSE16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Works
#undef  __Worka
//
#undef __FloatTmp
#undef __FloatSSE
#undef __FloatRet
#undef __FloatTmp16
#undef __FloatSSE16
#undef __FloatRet16
#undef __DoubleTmp
#undef __DoubleSSE
#undef __DoubleRet
#undef __DoubleTmp16
#undef __DoubleSSE16
#undef __DoubleRet16
// SSE::Array16
template<typename _F>
SSE::Array16<_F>&SSE::Array16<_F>::reset(size_t n)
{
  if(n==_N) return*this;
  size_t s = Top<_F>(n);
  if(s != _S) {
    if(_A) delete16(_A);
    const_cast<size_t&>(_S) = s;
    const_cast<_F*&>(_A) = new16<_F>(_S);
  }
  const_cast<size_t&>(_N) = n;
  return*this;
}

template<typename _F>
SSE::Array16<_F>&SSE::Array16<_F>::assign(Array16 const&a)
{
  if(a._S != _S) {
    if(_A) delete16(_A);
    const_cast<size_t&>(_S) = a._S;
    const_cast<_F*&>(_A) = new16<_F>(_S);
  }
  const_cast<size_t&>(_N) = a._N;
  std::memcpy(_A,a._A,_N*sizeof(_F));
  return*this;
}

template<typename _F>
SSE::Array16<_F>&SSE::Array16<_F>::operator=(Array16 const&B) WDutils_THROWING
{
  check_size(B,"operator=(Array16&)");
  std::memcpy(_A,B._A,_N*sizeof(_F));
  return*this;
}

template class WDutils::SSE::Array16<float>;
template class WDutils::SSE::Array16<double>;

//

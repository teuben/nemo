// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/src/sse.cc
///
/// \brief  implements sse.h
///
/// \author Walter Dehnen
///
/// \date   2010
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Walter Dehnen
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


#ifndef __SSE__
# warning SSE instructions not supported: will implement simple code
# ifdef __GNUC__
# warning perhaps using the -msse compiler option helps
# endif
#endif

#ifndef __SSE2__
# warning SSE2 instructions not supported: will implement simple code
# ifdef __GNUC__
# warning perhaps using the -msse2 compiler option helps
# endif
#endif

#ifdef __INTEL_COMPILER
#pragma warning (disable:981) /* operands are evaluated in unspecified order */
#endif

namespace WDutils {
  namespace SSE {
    // loading and storing 32-bit integers
#define _mm_load_epi32(A) \
    (__m128i)(_mm_load_ps(reinterpret_cast<float const*>(A)))
#define _mm_loadu_epi32(A) \
    (__m128i)(_mm_loadu_ps(reinterpret_cast<float const*>(A)))
#define _mm_store_epi32(M,A) \
    _mm_store_ps(reinterpret_cast<float*>(M),(__m128)(A))
#define _mm_storeu_epi32(M,A) \
    _mm_storeu_ps(reinterpret_cast<float*>(M),(__m128)(A))

#ifdef __SSE2__
# define __Loop2SSE(__NAME,__INIT,__WORKS,__WORKU,__WORKA,__TYPE)	\
    size_t i=size_t(f)&16;						\
    if(n<2 || i&8) {		      /* small or unaligned     */	\
      if(i&8) 								\
	WDutils_Warning("SSE::%s(): "					\
			"array of %ss @ %p not 8-byte aligned: "	\
			"unaligned (slow) SSE instructions forced\n",	\
			__NAME,__TYPE,f);				\
      i=0;								\
      size_t n2=bottom<2>(n);						\
      __INIT;								\
      for(; i!=n2; i+=2) { __WORKU; }   /* unaligned blocks work  */	\
      for(; i!=n ;  ++i) { __WORKS; }   /* final few elements     */	\
    } else {								\
      i >>= 3;                          /* begin: shifted pointer */	\
      n  += i;                          /* end:   shifted pointer */	\
      f  -= i;                          /* shift pointer          */	\
      size_t i2=top<2>(i);              /* begin: aligned blocks  */	\
      size_t n2=bottom<2>(n);           /* end:   aligned blocks  */	\
      __INIT;								\
      for(; i!=i2;  ++i) { __WORKS; }   /* STD work @ begin       */	\
      for(; i!=n2; i+=2) { __WORKA; }   /* SSE work in blocks     */	\
      for(; i!=n ;  ++i) { __WORKS; }   /* STD work @ end         */	\
    }
# define __Loop2SSE16(__INIT,__WORKA)					\
    __INIT;								\
    for(size_t i=0; i!=n; i+=2) { __WORKA; }
#endif 
#ifdef __SSE__
# define __Loop4SSE(__NAME,__INIT,__WORKS,__WORKU,__WORKA,__TYPE)	\
    size_t i=size_t(f)&16;						\
    if(n<4 || i&4) {		        /* small or unaligned     */	\
      if(i&4)								\
	WDutils_Warning("SSE::%s(): "					\
			"array of %ss @ %p not 4-byte aligned: "	\
			"unaligned (slow) SSE instructions forced\n",	\
			__NAME,__TYPE,f);				\
      i=0;								\
      size_t n4=bottom<4>(n);						\
      __INIT;								\
      for(; i!=n4; i+=4) { __WORKU; }   /* unaligned blocks work  */	\
      for(; i!=n ;  ++i) { __WORKS; }   /* final few elements     */	\
    } else {								\
      i >>= 2;                          /* begin: shifted pointer */	\
      n  += i;                          /* end:   shifted pointer */	\
      f  -= i;                          /* shift pointer          */	\
      size_t i4=top<4>(i);              /* begin: aligned blocks  */	\
      size_t n4=bottom<4>(n);           /* end:   aligned blocks  */	\
      __INIT;								\
      for(; i!=i4;  ++i) { __WORKS; }   /* STD work @ begin       */	\
      for(; i!=n4; i+=4) { __WORKA; }   /* SSE work in blocks     */	\
      for(; i!=n ;  ++i) { __WORKS; }   /* STD work @ end         */	\
    }
# define __Loop4SSE16(__INIT,__WORKA)					\
    __INIT;								\
    for(size_t i=0; i!=n; i+=4) { __WORKA; }
#endif 

#ifdef __SSE__
# define __LoopFloat(__NAME,__INIT,__WORKS,__WORKU,__WORKA)		\
    __Loop4SSE(__NAME,__INIT,__WORKS,__WORKU,__WORKA,"float")
# define __LoopFloat16(__INIT,__WORKS,__WORKA)				\
    __Loop4SSE16(__INIT,__WORKA)
#else
# define __LoopFloat(__NAME,__INIT,__WORKS,__WORKU,__WORKA)		\
    for(size_t i=0; i!=n; ++i) { __WORKS; }
# define __LoopFloat16(__INIT,__WORKS,__WORKA)				\
    for(size_t i=0; i!=n; ++i) { __WORKS; }
#endif

#ifdef __SSE2__
# define __LoopDouble(__NAME,__INIT,__WORKS,__WORKU,__WORKA)		\
    __Loop2SSE(__NAME,__INIT,__WORKS,__WORKU,__WORKA,"double")
# define __LoopDouble16(__INIT,__WORKS,__WORKA)				\
    __Loop2SSE16(__INIT,__WORKA)
# define __LoopInt(__NAME,__INIT,__WORKS,__WORKU,__WORKA)		\
    __Loop4SSE(__NAME,__INIT,__WORKS,__WORKU,__WORKA,"int")
# define __LoopInt16(__INIT,__WORKS,__WORKA)				\
    __Loop4SSE16(__INIT,__WORKA)
#else
# define __LoopDouble(__NAME,__INIT,__WORKS,__WORKU,__WORKA)		\
    for(size_t i=0; i!=n; ++i) { __WORKS; }
# define __LoopDouble16(__INIT,__WORKS,__WORKA)				\
    for(size_t i=0; i!=n; ++i) { __WORKS; }
# define __LoopInt(__NAME,__INIT,__WORKS,__WORKU,__WORKA)		\
    for(size_t i=0; i!=n; ++i) { __WORKS; }
# define __LoopInt16(__INIT,__WORKS,__WORKA)				\
    for(size_t i=0; i!=n; ++i) { __WORKS; }
#endif

#ifdef __SSE4_1__
# define __LoopXInt(__NAME,__INIT,__WORKS,__WORKU,__WORKA)		\
    __Loop4SSE(__NAME,__INIT,__WORKS,__WORKU,__WORKA,"int")
# define __LoopXInt16(__INIT,__WORKS,__WORKA)				\
    __Loop4SSE16(__INIT,__WORKA)
#else
# define __LoopXInt(__NAME,__INIT,__WORKS,__WORKU,__WORKA)		\
    for(size_t i=0; i!=n; ++i) { __WORKS; }
# define __LoopXInt16(__INIT,__WORKS,__WORKA)				\
    for(size_t i=0; i!=n; ++i) { __WORKS; }
#endif
    // for(i=0; i!=n; ++i) f[i]=v
#define __Works f[i]=v;

#define __Init  __m128i V = _mm_set1_epi32(v);
#define __Worku _mm_storeu_epi32(f+i,V);
#define __Worka _mm_store_epi32(f+i,V);
    void UnAligned::Ass(int*f, size_t n, int v)
    { __LoopInt("Ass",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Ass(int*f, size_t n, int v)
    { __LoopInt16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128 V = _mm_set1_ps(v);
#define __Worku _mm_storeu_ps(f+i,V);
#define __Worka _mm_store_ps(f+i,V);
    void UnAligned::Ass(float*f, size_t n, float v)
    { __LoopFloat("Ass",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Ass(float*f, size_t n, float v)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128d V = _mm_set1_pd(v);
#define __Worku _mm_storeu_pd(f+i,V);
#define __Worka _mm_store_pd(f+i,V);
    void UnAligned::Ass(double*f, size_t n, double v)
    { __LoopDouble("Ass",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Ass(double*f, size_t n, double v)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#undef  __Works
    // for(i=0; i!=n; ++i) a[i]=w*b[i]
#define __Works a[i]=w*b[i];

#define __Init  __m128i W = _mm_set1_epi32(w);
#define __Worka _mm_store_epi32(a+i,_mm_mul_epi32(W,_mm_load_epi32(b+i)));
    void Aligned::Ass(int*a, size_t n, int w, const int*b)
    { __LoopXInt16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worka

#define __Init  __m128 W = _mm_set1_ps(w);
#define __Worka _mm_store_ps(a+i,_mm_mul_ps(W,_mm_load_ps(b+i)));
    void Aligned::Ass(float*a, size_t n, float w, const float*b)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worka

#define __Init  __m128d W = _mm_set1_pd(w);
#define __Worka _mm_store_pd(a+i,_mm_mul_pd(W,_mm_load_pd(b+i)));
    void Aligned::Ass(double*a, size_t n, double w, const double*b)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worka

#undef  __Works
    // for(i=0; i!=n; ++i) f[i]=-f[i];
#define __Works f[i]=-f[i];

#define __Init  __m128i V = _mm_set1_epi32(0);
#define __Worku _mm_storeu_epi32(f+i,_mm_sub_epi32(V,_mm_loadu_epi32(f+i)));
#define __Worka _mm_store_epi32(f+i,_mm_sub_epi32(V,_mm_load_epi32(f+i)));
    void UnAligned::Neg(int*f, size_t n)
    { __LoopInt("Ass",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Neg(int*f, size_t n)
    { __LoopInt16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128 V = _mm_set1_ps(0.f);
#define __Worku _mm_storeu_ps(f+i,_mm_sub_ps(V,_mm_loadu_ps(f+i)));
#define __Worka _mm_store_ps(f+i,_mm_sub_ps(V,_mm_load_ps(f+i)));
    void UnAligned::Neg(float*f, size_t n)
    { __LoopFloat("Ass",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Neg(float*f, size_t n)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128d V = _mm_set1_pd(0.f);
#define __Worku _mm_storeu_pd(f+i,_mm_sub_pd(V,_mm_loadu_pd(f+i)));
#define __Worka _mm_store_pd(f+i,_mm_sub_pd(V,_mm_load_pd(f+i)));
    void UnAligned::Neg(double*f, size_t n)
    { __LoopDouble("Ass",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Neg(double*f, size_t n)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#undef  __Works
    // for(i=0; i!=n; ++i) f[i]+=v
#define __Works f[i]+=v;

#define __Init  __m128i V = _mm_set1_epi32(v);
#define __Worku _mm_storeu_epi32(f+i,_mm_add_epi32(V,_mm_loadu_epi32(f+i)));
#define __Worka _mm_store_epi32 (f+i,_mm_add_epi32(V,_mm_load_epi32 (f+i)));
    void UnAligned::Add(int*f, size_t n, int v)
    { __LoopInt("Add",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Add(int*f, size_t n, int v)
    { __LoopInt16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128 V = _mm_set1_ps(v);
#define __Worku _mm_storeu_ps(f+i,_mm_add_ps(V,_mm_loadu_ps(f+i)));
#define __Worka _mm_store_ps(f+i,_mm_add_ps(V,_mm_load_ps(f+i)));
    void UnAligned::Add(float*f, size_t n, float v)
    { __LoopFloat("Add",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Add(float*f, size_t n, float v)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128d V = _mm_set1_pd(v);
#define __Worku _mm_storeu_pd(f+i,_mm_add_pd(V,_mm_loadu_pd(f+i)));
#define __Worka _mm_store_pd(f+i,_mm_add_pd(V,_mm_load_pd(f+i)));
    void UnAligned::Add(double*f, size_t n, double v)
    { __LoopDouble("Add",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Add(double*f, size_t n, double v)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#undef  __Works
    // for(i=0; i!=n; ++i) f[i]+=b[i]
#define __Works a[i] += b[i];
#define __Init

#define __Worka _mm_store_epi32(a+i,_mm_add_epi32(_mm_load_epi32(a+i),	\
						  _mm_load_epi32(b+i)));
    void   Aligned::Add(int*a, size_t n, const int*b)
    { __LoopInt16(__Init,__Works,__Worka); }
#undef  __Worka

#define __Worka _mm_store_ps(a+i,_mm_add_ps(_mm_load_ps(a+i),	\
					    _mm_load_ps(b+i)));
    void   Aligned::Add(float*a, size_t n, const float*b)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Worka

#define __Worka _mm_store_pd(a+i,_mm_add_pd(_mm_load_pd(a+i),	\
					    _mm_load_pd(b+i)));
    void   Aligned::Add(double*a, size_t n, const double*b)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Worka

#undef  __Init
#undef  __Works
    // for(i=0; i!=n; ++i) f[i]+=w*b[i]
#define __Works a[i]+=w*b[i];

#define __Init  __m128i W = _mm_set1_epi32(w);
#define __Worka _mm_store_epi32(a+i,_mm_add_epi32(_mm_load_epi32(a+i),	\
		_mm_mul_epi32(W,_mm_load_epi32(b+i))));
    void   Aligned::Add(int*a, size_t n, int w, const int*b)
    { __LoopXInt16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worka

#define __Init  __m128 W = _mm_set1_ps(w);
#define __Worka _mm_store_ps(a+i,_mm_add_ps(_mm_load_ps(a+i),		 \
					    _mm_mul_ps(W,_mm_load_ps(b+i))));
    void   Aligned::Add(float*a, size_t n, float w, const float*b)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worka

#define __Init  __m128d W = _mm_set1_pd(w);
#define __Worka _mm_store_pd(a+i,_mm_add_pd(_mm_load_pd(a+i),		 \
					    _mm_mul_pd(W,_mm_load_pd(b+i))));
    void   Aligned::Add(double*a, size_t n, double w, const double*b)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worka

#undef  __Works
    // for(i=0; i!=n; ++i) f[i]-=v;
#define __Works f[i]-=v;

#define __Init  __m128i V = _mm_set1_epi32(v);
#define __Worku _mm_storeu_epi32(f+i,_mm_sub_epi32(_mm_loadu_epi32(f+i),V));
#define __Worka _mm_store_epi32(f+i,_mm_sub_epi32(_mm_load_epi32(f+i),V));
    void UnAligned::Sub(int*f, size_t n, int v)
    { __LoopInt("Sub",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Sub(int*f, size_t n, int v)
    { __LoopInt16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128 V = _mm_set1_ps(v);
#define __Worku _mm_storeu_ps(f+i,_mm_sub_ps(_mm_loadu_ps(f+i),V));
#define __Worka _mm_store_ps(f+i,_mm_sub_ps(_mm_load_ps(f+i),V));
    void UnAligned::Sub(float*f, size_t n, float v)
    { __LoopFloat("Sub",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Sub(float*f, size_t n, float v)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128d V = _mm_set1_pd(v);
#define __Worku _mm_storeu_pd(f+i,_mm_sub_pd(_mm_loadu_pd(f+i),V));
#define __Worka _mm_store_pd(f+i,_mm_sub_pd(_mm_load_pd(f+i),V));
    void UnAligned::Sub(double*f, size_t n, double v)
    { __LoopDouble("Sub",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Sub(double*f, size_t n, double v)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#undef  __Works
    // for(i=0; i!=n; ++i) f[i]-=b[i];
#define __Works a[i]-=b[i];
#define __Init

#define __Worka _mm_store_epi32(a+i,_mm_sub_epi32(_mm_load_epi32(a+i),	\
						  _mm_load_epi32(b+i)));
    void   Aligned::Sub(int*a, size_t n, const int*b)
    { __LoopInt16(__Init,__Works,__Worka); }
#undef  __Worka

#define __Init
#define __Worka _mm_store_ps(a+i,_mm_sub_ps(_mm_load_ps(a+i),	\
					    _mm_load_ps(b+i)));
    void   Aligned::Sub(float*a, size_t n, const float*b)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Worka

#define __Worka _mm_store_pd(a+i,_mm_sub_pd(_mm_load_pd(a+i),	\
					    _mm_load_pd(b+i)));
    void   Aligned::Sub(double*a, size_t n, const double*b)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Worka

#undef  __Init
#undef  __Works
    // for(i=0; i!=n; ++i) f[i]-=w*b[i];
#define __Works a[i] -= w*b[i];

#define __Init  __m128i W = _mm_set1_epi32(w);
#define __Worka _mm_store_epi32(a+i,_mm_sub_epi32(_mm_load_epi32(a+i),	\
	        _mm_mul_epi32(W,_mm_load_epi32(b+i))));
    void   Aligned::Sub(int*a, size_t n, int w, const int*b)
    { __LoopXInt16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worka

#define __Init  __m128 W = _mm_set1_ps(w);
#define __Worka _mm_store_ps(a+i,_mm_sub_ps(_mm_load_ps(a+i),		\
					    _mm_mul_ps(W,_mm_load_ps(b+i))));
    void   Aligned::Sub(float*a, size_t n, float w, const float*b)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worka

#define __Init  __m128d W = _mm_set1_pd(w);
#define __Worka _mm_store_pd(a+i,_mm_sub_pd(_mm_load_pd(a+i),		\
					    _mm_mul_pd(W,_mm_load_pd(b+i))));
    void   Aligned::Sub(double*a, size_t n, double w, const double*b)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worka

#undef  __Works
    // for(i=0; i!=n; ++i) f[i]*=v;
#define __Works f[i]*=v;

#define __Init  __m128 V = _mm_set1_epi32(v);
#define __Worku _mm_storeu_epi32(f+i,_mm_mul_epi32(_mm_loadu_epi32(f+i),V));
#define __Worka _mm_store_epi32(f+i,_mm_mul_epi32(_mm_load_epi32(f+i),V));
    void UnAligned::Mul(int*f, size_t n, int v)
    { __LoopXInt("Mul",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Mul(int*f, size_t n, int v)
    { __LoopXInt16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128 V = _mm_set1_ps(v);
#define __Worku _mm_storeu_ps(f+i,_mm_mul_ps(_mm_loadu_ps(f+i),V));
#define __Worka _mm_store_ps(f+i,_mm_mul_ps(_mm_load_ps(f+i),V));
    void UnAligned::Mul(float*f, size_t n, float v)
    { __LoopFloat("Mul",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Mul(float*f, size_t n, float v)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128d V = _mm_set1_pd(v);
#define __Worku _mm_storeu_pd(f+i,_mm_mul_pd(_mm_loadu_pd(f+i),V));
#define __Worka _mm_store_pd(f+i,_mm_mul_pd(_mm_load_pd(f+i),V));
    void UnAligned::Mul(double*f, size_t n, double v)
    { __LoopDouble("Mul",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Mul(double*f, size_t n, double v)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#undef  __Works
    // for(i=0; i!=n; ++i) f[i]*=b[i];
#define __Works a[i]*=b[i];
#define __Init

#define __Worka _mm_store_epi32(a+i,_mm_mul_epi32(_mm_load_epi32(a+i),	\
						  _mm_load_epi32(b+i)));
    void   Aligned::Mul(int*a, size_t n, const int*b)
    { __LoopXInt16(__Init,__Works,__Worka); }
#undef  __Worka

#define __Worka _mm_store_ps(a+i,_mm_mul_ps(_mm_load_ps(a+i),_mm_load_ps(b+i)));
    void   Aligned::Mul(float*a, size_t n, const float*b)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Worka

#define __Worka _mm_store_pd(a+i,_mm_mul_pd(_mm_load_pd(a+i),_mm_load_pd(b+i)));
    void   Aligned::Mul(double*a, size_t n, const double*b)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Worka

#undef  __Init
#undef  __Works
    // for(i=0; i!=n; ++i) a[i]/=b[i];
#define __Works a[i]/=b[i];
#define __Init

#define __Worka _mm_store_ps(a+i,_mm_div_ps(_mm_load_ps(a+i),_mm_load_ps(b+i)));
    void   Aligned::Div(float*a, size_t n, const float*b)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Worka

#define __Worka _mm_store_pd(a+i,_mm_div_pd(_mm_load_pd(a+i),_mm_load_pd(b+i)));
    void   Aligned::Div(double*a, size_t n, const double*b)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Worka

#undef  __Init
#undef  __Works
    // for(i=0; i!=n; ++i) f[i]=x/f[i];
#define __Works f[i]=x/f[i];

#define __Init  __m128 X = _mm_set1_ps(x);
#define __Worku _mm_storeu_ps(f+i,_mm_div_ps(X,_mm_loadu_ps(f+i)));
#define __Worka _mm_store_ps(f+i,_mm_div_ps(X,_mm_load_ps(f+i)));
    void UnAligned::Inv(float*f, size_t n, float x)
    { __LoopFloat("Inv",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Inv(float*f, size_t n, float x)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#define __Init  __m128d X = _mm_set1_pd(x);
#define __Worku _mm_storeu_pd(f+i,_mm_div_pd(X,_mm_loadu_pd(f+i)));
#define __Worka _mm_store_pd(f+i,_mm_div_pd(X,_mm_load_pd(f+i)));
    void UnAligned::Inv(double*f, size_t n, double x)
    { __LoopDouble("Inv",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Inv(double*f, size_t n, double x)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worku
#undef  __Worka

#undef  __Works
    // for(i=0; i!=n; ++i) a[i]=x/b[i];
#define __Works a[i]=x/b[i];

#define __Init  __m128 X = _mm_set1_ps(x);
#define __Worka _mm_store_ps(a+i,_mm_div_ps(X,_mm_load_ps(b+i)));
    void   Aligned::Inv(float*a, size_t n, float x, const float*b)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worka

#define __Init  __m128d X = _mm_set1_pd(x);
#define __Worka _mm_store_pd(a+i,_mm_div_pd(X,_mm_load_pd(b+i)));
    void   Aligned::Inv(double*a, size_t n, double x, const double*b)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Init
#undef  __Worka

#undef  __Works
    // for(i=0; i!=n; ++i) f[i]=std::sqrt(f[i])
#define __Works f[i]=std::sqrt(f[i]);
#define __Init

#define __Worku _mm_storeu_ps(f+i,_mm_sqrt_ps(_mm_loadu_ps(f+i)));
#define __Worka _mm_store_ps(f+i,_mm_sqrt_ps(_mm_load_ps(f+i)));
    void UnAligned::Sqrt(float*f, size_t n)
    { __LoopFloat("Sqrt",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Sqrt(float*f, size_t n)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Worku
#undef  __Worka

#define __Worku _mm_storeu_pd(f+i,_mm_sqrt_pd(_mm_loadu_pd(f+i)));
#define __Worka _mm_store_pd(f+i,_mm_sqrt_pd(_mm_load_pd(f+i)));
    void UnAligned::Sqrt(double*f, size_t n)
    { __LoopDouble("Sqrt",__Init,__Works,__Worku,__Worka); }
    void   Aligned::Sqrt(double*f, size_t n)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Worku
#undef  __Worka

#undef  __Init
#undef  __Works
    // for(i=0; i!=n; ++i) a[i]=std::sqrt(b[i])
#define __Works a[i]=std::sqrt(b[i]);
#define __Init

#define __Worka _mm_store_ps(a+i,_mm_sqrt_ps(_mm_load_ps(b+i)));
    void   Aligned::Sqrt(float*a, size_t n, const float*b)
    { __LoopFloat16(__Init,__Works,__Worka); }
#undef  __Worka

#define __Worka _mm_store_pd(a+i,_mm_sqrt_pd(_mm_load_pd(b+i)));
    void   Aligned::Sqrt(double*a, size_t n, const double*b)
    { __LoopDouble16(__Init,__Works,__Worka); }
#undef  __Worka

#undef  __Init
#undef  __Works
    // S=0; for(i=0; i!=n; ++i) S+=f[i]; return S;
#define __Works tm+=f[i];

#ifdef __SSE__
# define __FloatTmp  __m128 Tm = _mm_setzero_ps();float tm(0.f)
# define __FloatRet  return _mm_getsum_ps(Tm) + tm
#else
# define __FloatTmp  float tm(0.f);	
# define __FloatRet  return tm;
#endif

#ifdef __SSE2__

# define __IntTmp    __m128i Tm = _mm_setzero_si128(); int tm(0)
# define __IntRet    return _mm_getsum_epi32(Tm) + tm

# define __DoubleTmp __m128d Tm = _mm_setzero_pd(); double tm(0.0)
# define __DoubleRet							\
    __attribute__ ((aligned(16))) double tmp[2];			\
    _mm_store_pd(tmp,Tm);						\
    return tmp[0]+tmp[1] + tm

#else

# define __IntTmp    int tm(0)
# define __IntRet    return tm

# define __DoubleTmp double tm(0.0)
# define __DoubleRet return tm

#endif

#define __Init

#define __Worka Tm = _mm_add_epi32(Tm,_mm_load_epi32(f+i));
#define __Worku Tm = _mm_add_epi32(Tm,_mm_loadu_epi32(f+i));
    int UnAligned::Sum(const int*f, size_t n)
    {
      __IntTmp;
      __LoopInt("Sum",__Init,__Works,__Worku,__Worka);
      __IntRet;
    }
    int   Aligned::Sum(const int*f, size_t n)
    {
      __IntTmp;
      __LoopInt16(__Init,__Works,__Worka);
      __IntRet;
    }
#undef  __Worku
#undef  __Worka

#define __Worka Tm = _mm_add_ps(Tm,_mm_load_ps(f+i));
#define __Worku Tm = _mm_add_ps(Tm,_mm_loadu_ps(f+i));
    float UnAligned::Sum(const float*f, size_t n)
    {
      __FloatTmp;
      __LoopFloat("Sum",__Init,__Works,__Worku,__Worka);
      __FloatRet;
    }
    float   Aligned::Sum(const float*f, size_t n)
    {
      __FloatTmp;
      __LoopFloat16(__Init,__Works,__Worka);
      __FloatRet;
    }
#undef  __Worku
#undef  __Worka

#define __Worka Tm = _mm_add_pd(Tm,_mm_load_pd(f+i));
#define __Worku Tm = _mm_add_pd(Tm,_mm_loadu_pd(f+i));
    double UnAligned::Sum(const double*f, size_t n)
    {
      __DoubleTmp;
      __LoopDouble("Sum",__Init,__Works,__Worku,__Worka);
      __DoubleRet;
    }
    double   Aligned::Sum(const double*f, size_t n)
    {
      __DoubleTmp;
      __LoopDouble16(__Init,__Works,__Worka);
      __DoubleRet;
    }
#undef  __Worku
#undef  __Worka

#undef  __Works
    // S=0; for(i=0; i!=n; ++i) S+=f[i]*f[i]; return S;
#define __Works tm+=f[i]*f[i];

#define __Worka __m128i iF=_mm_load_epi32(f+i); \
    Tm = _mm_add_epi32(Tm,_mm_mul_epi32(iF,iF));
#define __Worku __m128i iF=_mm_loadu_epi32(f+i); \
    Tm = _mm_add_epi32(Tm,_mm_mul_epi32(iF,iF));
    int UnAligned::Norm(const int*f, size_t n)
    {
      __IntTmp;
      __LoopXInt("Norm",__Init,__Works,__Worku,__Worka);
      __IntRet;
    }
    int   Aligned::Norm(const int*f, size_t n)
    {
      __IntTmp;
      __LoopXInt16(__Init,__Works,__Worka);
      __IntRet;
    }
#undef  __Worku
#undef  __Worka

#define __Worka __m128 iF=_mm_load_ps(f+i); \
    Tm = _mm_add_ps(Tm,_mm_mul_ps(iF,iF));
#define __Worku __m128 iF=_mm_loadu_ps(f+i); \
    Tm = _mm_add_ps(Tm,_mm_mul_ps(iF,iF));
    float UnAligned::Norm(const float*f, size_t n)
    {
      __FloatTmp;
      __LoopFloat("Norm",__Init,__Works,__Worku,__Worka);
      __FloatRet;
    }
    float   Aligned::Norm(const float*f, size_t n)
    {
      __FloatTmp;
      __LoopFloat16(__Init,__Works,__Worka);
      __FloatRet;
    }
#undef  __Worku
#undef  __Worka

#define __Worka __m128d iF=_mm_load_pd(f+i); \
    Tm = _mm_add_pd(Tm,_mm_mul_pd(iF,iF));
#define __Worku __m128d iF=_mm_loadu_pd(f+i); \
    Tm = _mm_add_pd(Tm,_mm_mul_pd(iF,iF));
    double UnAligned::Norm(const double*f, size_t n)
    {
      __DoubleTmp;
      __LoopDouble("Norm",__Init,__Works,__Worku,__Worka);
      __DoubleRet;
    }
    double   Aligned::Norm(const double*f, size_t n)
    {
      __DoubleTmp;
      __LoopDouble16(__Init,__Works,__Worka);
      __DoubleRet;
    }
#undef  __Worku
#undef  __Worka

#undef  __Works
    // S=0; for(i=0; i!=n; ++i) S+=a[i]*b[i]; return S;
#define __Works tm += a[i]*b[i];

#define __Worka Tm = _mm_add_epi32(Tm,_mm_mul_epi32(_mm_load_epi32(a+i), \
						    _mm_load_epi32(b+i)));
    int   Aligned::Dot(const int*a, size_t n, const int*b)
    {
      __IntTmp;
      __LoopXInt16(__Init,__Works,__Worka);
      __IntRet;
    }
#undef  __Worka

#define __Worka Tm = _mm_add_ps(Tm,_mm_mul_ps(_mm_load_ps(a+i),		\
					      _mm_load_ps(b+i)));
    float   Aligned::Dot(const float*a, size_t n, const float*b)
    {
      __FloatTmp;
      __LoopFloat16(__Init,__Works,__Worka);
      __FloatRet;
    }
#undef  __Worka

#define __Worka Tm = _mm_add_pd(Tm,_mm_mul_pd(_mm_load_pd(a+i),		\
					      _mm_load_pd(b+i)));
    double   Aligned::Dot(const double*a, size_t n, const double*b)
    {
      __DoubleTmp;
      __LoopDouble16(__Init,__Works,__Worka);
      __DoubleRet;
    }
#undef  __Worka

#undef  __Init
#undef  __IntTmp
#undef  __IntRet
#undef  __FloatTmp
#undef  __FloatRet
#undef  __DoubleTmp
#undef  __DoubleRet
#undef  __Works
    //
#undef __Loop2SSE
#undef __Loop2SSE16
#undef __Loop4SSE
#undef __Loop4SSE16
#undef __LoopInt
#undef __LoopInt16
#undef __LoopFloat
#undef __LoopFloat16
#undef __LoopDouble
#undef __LoopDouble16
    // SSE::Array16
    template<typename _F>
    Array16<_F>&Array16<_F>::reset(size_t n)
    {
      if(n==_N) return*this;
      size_t s = Top<_F>(n);
      if(s != _S) {
	if(_A) WDutils_DEL16(_A);
	const_cast<size_t&>(_S) = s;
	const_cast<_F*&>(_A) = WDutils_NEW16(_F,_S);
      }
      const_cast<size_t&>(_N) = n;
      return*this;
    }

    template<typename _F>
    Array16<_F>&Array16<_F>::assign(Array16 const&a)
    {
      if(a._S != _S) {
	if(_A) WDutils_DEL16(_A);
	const_cast<size_t&>(_S) = a._S;
	const_cast<_F*&>(_A) = WDutils_NEW16(_F,_S);
      }
      const_cast<size_t&>(_N) = a._N;
      std::memcpy(_A,a._A,_N*sizeof(_F));
      return*this;
    }

    template<typename _F>
    Array16<_F>&Array16<_F>::operator=(Array16 const&B) WDutils_THROWING
    {
      check_size(B,"operator=(Array16&)");
      std::memcpy(_A,B._A,_N*sizeof(_F));
      return*this;
    }

    template class Array16<float>;
    template class Array16<double>;

    /// SSE::Swap16
#ifdef __SSE__
    void __swap16(float*a, float*b, size_t n)
    {
      for(float*an=a+n; a<an; a+=4,b+=4) {
	__m128 t = _mm_load_ps(a);
	_mm_store_ps(a,_mm_load_ps(b));
	_mm_store_ps(b,t);
      }
    }
#endif
} } // namespace WDutils::SSE
//

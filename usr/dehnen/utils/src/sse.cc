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

using namespace WDutils;

#ifdef __SSE__
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
#else
# warning SSE instructions not supported: will implement simple code
# define __FloatSSE(__NAME,__INIT,__WORKS,__WORKU,__WORKA)		\
  for(size_t i=0; i!=n; ++i) { __WORKS; }
#endif

#ifdef __SSE2__
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
#else
# warning SSE2 instructions not supported: will implement simple code
# define __DoubleSSE(__NAME,__INIT,__WORKS,__WORKU,__WORKA)		\
  for(size_t i=0; i!=n; ++i) { __WORKS; }
#endif

// SSE::Assign
void SSE::Assign(float*f, size_t n, float v) {
#define __Works f[i] = v;
#define __Init  __m128 V = _mm_set1_ps(v);
#define __Worku _mm_storeu_ps(f+i,V);
#define __Worka _mm_store_ps(f+i,V);
  __FloatSSE("Assign",__Init,__Works,__Worku,__Worka);
#undef  __Init
#undef  __Worku
#undef  __Worka
}
void SSE::Assign(double*f, size_t n, double v) {
#define __Init  __m128d V = _mm_set1_pd(v);
#define __Worku _mm_storeu_pd(f+i,V);
#define __Worka _mm_store_pd(f+i,V);
  __DoubleSSE("Assign",__Init,__Works,__Worku,__Worka);
#undef  __Init
#undef  __Worku
#undef  __Worka
#undef  __Works
}
// SSE::Add
void SSE::Add(float*f, size_t n, float v) {
#define __Works f[i] += v;
#define __Init  __m128 V = _mm_set1_ps(v);
#define __Worku _mm_storeu_ps(f+i,_mm_add_ps(_mm_loadu_ps(f+i),V));
#define __Worka _mm_store_ps(f+i,_mm_add_ps(_mm_load_ps(f+i),V));
  __FloatSSE("Add",__Init,__Works,__Worku,__Worka);
#undef  __Init
#undef  __Worku
#undef  __Worka
}
void SSE::Add(double*f, size_t n, double v) {
#define __Init  __m128d V = _mm_set1_pd(v);
#define __Worku _mm_storeu_pd(f+i,_mm_add_pd(_mm_loadu_pd(f+i),V));
#define __Worka _mm_store_pd(f+i,_mm_add_pd(_mm_load_pd(f+i),V));
  __DoubleSSE("Add",__Init,__Works,__Worku,__Worka);
#undef  __Init
#undef  __Worku
#undef  __Worka
#undef  __Works
}
// SSE::Subtract
void SSE::Subtract(float*f, size_t n, float v) {
#define __Works f[i] -= v;
#define __Init  __m128 V = _mm_set1_ps(v);
#define __Worku _mm_storeu_ps(f+i,_mm_sub_ps(_mm_loadu_ps(f+i),V));
#define __Worka _mm_store_ps(f+i,_mm_sub_ps(_mm_load_ps(f+i),V));
  __FloatSSE("Subtract",__Init,__Works,__Worku,__Worka);
#undef  __Init
#undef  __Worku
#undef  __Worka
}
void SSE::Subtract(double*f, size_t n, double v) {
#define __Init  __m128d V = _mm_set1_pd(v);
#define __Worku _mm_storeu_pd(f+i,_mm_sub_pd(_mm_loadu_pd(f+i),V));
#define __Worka _mm_store_pd(f+i,_mm_sub_pd(_mm_load_pd(f+i),V));
  __DoubleSSE("Subtract",__Init,__Works,__Worku,__Worka);
#undef  __Init
#undef  __Worku
#undef  __Worka
#undef  __Works
}
// SSE::Multiply
void SSE::Multiply(float*f, size_t n, float v) {
#define __Works f[i] *= v;
#define __Init  __m128 V = _mm_set1_ps(v);
#define __Worku _mm_storeu_ps(f+i,_mm_mul_ps(_mm_loadu_ps(f+i),V));
#define __Worka _mm_store_ps(f+i,_mm_mul_ps(_mm_load_ps(f+i),V));
  __FloatSSE("Multiply",__Init,__Works,__Worku,__Worka);
#undef  __Init
#undef  __Worku
#undef  __Worka
}
void SSE::Multiply(double*f, size_t n, double v) {
#define __Init  __m128d V = _mm_set1_pd(v);
#define __Worku _mm_storeu_pd(f+i,_mm_mul_pd(_mm_loadu_pd(f+i),V));
#define __Worka _mm_store_pd(f+i,_mm_mul_pd(_mm_load_pd(f+i),V));
  __DoubleSSE("Multiply",__Init,__Works,__Worku,__Worka);
#undef  __Init
#undef  __Worku
#undef  __Worka
#undef  __Works
}
// SSE::Invert
void SSE::Invert(float*f, size_t n, float v) {
#define __Works f[i] = v/f[i];
#define __Init  __m128 V = _mm_set1_ps(v);
#define __Worku _mm_storeu_ps(f+i,_mm_div_ps(V,_mm_loadu_ps(f+i)));
#define __Worka _mm_store_ps(f+i,_mm_div_ps(V,_mm_load_ps(f+i)));
  __FloatSSE("Invert",__Init,__Works,__Worku,__Worka);
#undef  __Init
#undef  __Worku
#undef  __Worka
}
void SSE::Invert(double*f, size_t n, double v) {
#define __Init  __m128d V = _mm_set1_pd(v);
#define __Worku _mm_storeu_pd(f+i,_mm_div_pd(V,_mm_loadu_pd(f+i)));
#define __Worka _mm_store_pd(f+i,_mm_div_pd(V,_mm_load_pd(f+i)));
  __DoubleSSE("Invert",__Init,__Works,__Worku,__Worka);
#undef  __Init
#undef  __Worku
#undef  __Worka
#undef  __Works
}
// SSE::SquareRoot
void SSE::SquareRoot(float*f, size_t n) {
#define __Works f[i] = std::sqrt(f[i]);
#define __Init
#define __Worku _mm_storeu_ps(f+i,_mm_sqrt_ps(_mm_loadu_ps(f+i)));
#define __Worka _mm_store_ps(f+i,_mm_sqrt_ps(_mm_load_ps(f+i)));
  __FloatSSE("SquareRoot",__Init,__Works,__Worku,__Worka);
#undef  __Worku
#undef  __Worka
}
void SSE::SquareRoot(double*f, size_t n) {
#define __Worku _mm_storeu_pd(f+i,_mm_sqrt_pd(_mm_loadu_pd(f+i)));
#define __Worka _mm_store_pd(f+i,_mm_sqrt_pd(_mm_load_pd(f+i)));
  __DoubleSSE("SquareRoot",__Init,__Works,__Worku,__Worka);
#undef  __Init
#undef  __Worku
#undef  __Worka
#undef  __Works
}
//

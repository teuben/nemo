// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/meta.h
///
/// \brief  support for metaprogramming
///
/// \author Walter Dehnen
///                                                                             
/// \date   2008-2009
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008-2009 Walter Dehnen
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
#ifndef WDutils_included_meta_h
#define WDutils_included_meta_h

#ifndef WDutils_included_exception_h
# include <exception.h>
#endif

namespace WDutils {
  namespace meta {
    ///
    /// support for operation count
    class OpCounting {
      static int M, A, D, S;
    public:
      static void reset() { M=0u; A=0u; D=0u; S=0u; }
      static void inc_mul() { ++ M; }
      static void inc_add() { ++ A; }
      static void inc_div() { ++ D; }
      static void inc_sqr() { ++ S; }
      static void count_ma(int m, int a) { M+=m; A+=a; }
      static void count_ds(int d, int s) { D+=d; S+=s; }
      static int mul() { return M; }
      static int add() { return A; }
      static int div() { return D; }
      static int sqr() { return S; }
    };
#ifdef WDutilsMetaCountOperations
# define __ResetC()     WDutils::meta::OpCounting::reset()
# define __CountMA(M,A) WDutils::meta::OpCounting::count_ma(M,A)
# define __CountDS(D,S) WDutils::meta::OpCounting::count_ds(D,S)
# define __IncMul()     WDutils::meta::OpCounting::inc_mul()
# define __IncAdd()     WDutils::meta::OpCounting::inc_add()
# define __IncDiv()     WDutils::meta::OpCounting::inc_div()
# define __IncSqr()     WDutils::meta::OpCounting::inc_sqr()
#else
# define __ResetC();
# define __CountMA(M,A);
# define __CountDS(D,S);
# define __IncMul();
# define __IncAdd();
# define __IncDiv();
# define __IncSqr();
#endif
    /// Inverse of an integer: replace integer divisions of real numbers by
    /// multiplications to generate faster code.
    /// \note instantinated for all integers up to 100
    template<int N> struct IntegerInverse
#ifndef WDutilsMetaNoDefaultIntegerInverse
    {
      WDutilsStaticAssert((N!=0 && N!=1));
      /// inverse of integer
      template<typename Real> Real Inverse() {
	__IncDiv();
	return Real(1)/Real(N);
      }
      /// divide scalar
      /// \param[in,out] X scalar, replaced by X/N on output
      template<typename Real>
      static void Divide(Real&X) {
	__IncMul();
	X *= Inverse<Real>();
      }
      /// ratio of scalar and N
      /// \param[in] X scalar
      /// \return    X/N
      template<typename Real>
      static Real Ratio(Real X) {
	__IncMul();
	return X*Inverse<Real>();
      }
    }
#endif
    ;
    template<> struct IntegerInverse<0> {};
    template<> struct IntegerInverse<1> {
      template<typename Real> static Real Inverse() { return Real(1); }
      template<typename Real> static void Divide(Real& ) {}
      template<typename Real> static Real Ratio(Real X) { return X; }
    };
#define DEFINVERSE(NUM,INVERSE)					\
    template<> struct IntegerInverse<NUM> {			\
      template<typename Real> static Real Inverse() {		\
	return Real(INVERSE); }					\
      template<typename Real> static void Divide(Real&X) {	\
	__IncMul(); X *= Real(INVERSE); }			\
      template<typename Real> static Real Ratio(Real X) {	\
	__IncMul(); return X * Real(INVERSE); }			\
    }
    // all numbers up to 100
    DEFINVERSE(   2, 0.5);
    DEFINVERSE(   3, 0.333333333333333333333);
    DEFINVERSE(   4, 0.25);
    DEFINVERSE(   5, 0.2);
    DEFINVERSE(   6, 0.166666666666666666667);
    DEFINVERSE(   7, 0.142857142857142857143);
    DEFINVERSE(   8, 0.125);
    DEFINVERSE(   9, 0.111111111111111111111);
    DEFINVERSE(  10, 0.1);
    DEFINVERSE(  11, 0.090909090909090909091);
    DEFINVERSE(  12, 0.083333333333333333333);
    DEFINVERSE(  13, 0.076923076923076923077);
    DEFINVERSE(  14, 0.071428571428571428571);
    DEFINVERSE(  15, 0.066666666666666666667);
    DEFINVERSE(  16, 0.0625);
    DEFINVERSE(  17, 0.058823529411764705882);
    DEFINVERSE(  18, 0.055555555555555555556);
    DEFINVERSE(  19, 0.052631578947368421053);
    DEFINVERSE(  20, 0.05);
    DEFINVERSE(  21, 0.047619047619047619048);
    DEFINVERSE(  22, 0.045454545454545454545);
    DEFINVERSE(  23, 0.043478260869565217391);
    DEFINVERSE(  24, 0.041666666666666666667);
    DEFINVERSE(  25, 0.04);
    DEFINVERSE(  26, 0.038461538461538461538);
    DEFINVERSE(  27, 0.037037037037037037037);
    DEFINVERSE(  28, 0.035714285714285714286);
    DEFINVERSE(  29, 0.034482758620689655172);
    DEFINVERSE(  30, 0.033333333333333333333);
    DEFINVERSE(  31, 0.032258064516129032258);
    DEFINVERSE(  32, 0.03125);
    DEFINVERSE(  33, 0.030303030303030303030);
    DEFINVERSE(  34, 0.029411764705882352941);
    DEFINVERSE(  35, 0.028571428571428571429);
    DEFINVERSE(  36, 0.027777777777777777778);
    DEFINVERSE(  37, 0.027027027027027027027);
    DEFINVERSE(  38, 0.026315789473684210526);
    DEFINVERSE(  39, 0.025641025641025641026);
    DEFINVERSE(  40, 0.025);
    DEFINVERSE(  41, 0.024390243902439024390);
    DEFINVERSE(  42, 0.023809523809523809524);
    DEFINVERSE(  44, 0.022727272727272727273);
    DEFINVERSE(  45, 0.022222222222222222222);
    DEFINVERSE(  46, 0.021739130434782608696);
    DEFINVERSE(  47, 0.021276595744680851064);
    DEFINVERSE(  48, 0.020833333333333333333);
    DEFINVERSE(  49, 0.020408163265306122449);
    DEFINVERSE(  50, 0.02);
    DEFINVERSE(  51, 0.019607843137254901961);
    DEFINVERSE(  52, 0.019230769230769230769);
    DEFINVERSE(  53, 0.018867924528301886792);
    DEFINVERSE(  54, 0.018518518518518518519);
    DEFINVERSE(  55, 0.018181818181818181818);
    DEFINVERSE(  56, 0.017857142857142857143);
    DEFINVERSE(  57, 0.017543859649122807018);
    DEFINVERSE(  58, 0.017241379310344827586);
    DEFINVERSE(  59, 0.016949152542372881356);
    DEFINVERSE(  60, 0.016666666666666666667);
    DEFINVERSE(  61, 0.016393442622950819672);
    DEFINVERSE(  62, 0.016129032258064516129);
    DEFINVERSE(  63, 0.015873015873015873016);
    DEFINVERSE(  64, 0.015625);
    DEFINVERSE(  65, 0.015384615384615384615);
    DEFINVERSE(  66, 0.015151515151515151515);
    DEFINVERSE(  67, 0.014925373134328358209);
    DEFINVERSE(  68, 0.014705882352941176471);
    DEFINVERSE(  69, 0.014492753623188405797);
    DEFINVERSE(  70, 0.014285714285714285714);
    DEFINVERSE(  71, 0.014084507042253521127);
    DEFINVERSE(  72, 0.013888888888888888889);
    DEFINVERSE(  73, 0.013698630136986301370);
    DEFINVERSE(  74, 0.013513513513513513514);
    DEFINVERSE(  75, 0.013333333333333333333);
    DEFINVERSE(  76, 0.013157894736842105263);
    DEFINVERSE(  77, 0.012987012987012987013);
    DEFINVERSE(  78, 0.012820512820512820513);
    DEFINVERSE(  79, 0.012658227848101265823);
    DEFINVERSE(  80, 0.0125);
    DEFINVERSE(  81, 0.012345679012345679012);
    DEFINVERSE(  82, 0.012195121951219512195);
    DEFINVERSE(  83, 0.012048192771084337349);
    DEFINVERSE(  84, 0.011904761904761904762);
    DEFINVERSE(  85, 0.011764705882352941176);
    DEFINVERSE(  86, 0.011627906976744186047);
    DEFINVERSE(  87, 0.011494252873563218391);
    DEFINVERSE(  88, 0.011363636363636363636);
    DEFINVERSE(  89, 0.011235955056179775281);
    DEFINVERSE(  90, 0.011111111111111111111);
    DEFINVERSE(  91, 0.010989010989010989011);
    DEFINVERSE(  92, 0.010869565217391304348);
    DEFINVERSE(  93, 0.010752688172043010753);
    DEFINVERSE(  94, 0.010638297872340425532);
    DEFINVERSE(  95, 0.010526315789473684211);
    DEFINVERSE(  96, 0.010416666666666666667);
    DEFINVERSE(  97, 0.010309278350515463918);
    DEFINVERSE(  98, 0.010204081632653061224);
    DEFINVERSE(  99, 0.010101010101010101010);
    DEFINVERSE( 100, 0.01);
    // all products of two numbers up to 16 and some more
    DEFINVERSE( 104, 0.0096153846153846153846);
    DEFINVERSE( 105, 0.0095238095238095238095);
    DEFINVERSE( 108, 0.0092592592592592592593);
    DEFINVERSE( 110, 0.0090909090909090909091);
    DEFINVERSE( 112, 0.0089285714285714285714);
    DEFINVERSE( 115, 0.0086956521739130434783);
    DEFINVERSE( 117, 0.0085470085470085470085);
    DEFINVERSE( 119, 0.0084033613445378151260);
    DEFINVERSE( 120, 0.0083333333333333333333);
    DEFINVERSE( 121, 0.0082644628099173553719);
    DEFINVERSE( 125, 0.008);
    DEFINVERSE( 126, 0.0079365079365079365079);
    DEFINVERSE( 128, 0.0078125);
    DEFINVERSE( 130, 0.0076923076923076923077);
    DEFINVERSE( 131, 0.0076335877862595419847);
    DEFINVERSE( 132, 0.0075757575757575757576);
    DEFINVERSE( 133, 0.0075187969924812030075);
    DEFINVERSE( 135, 0.0074074074074074074074);
    DEFINVERSE( 140, 0.0071428571428571428571);
    DEFINVERSE( 143, 0.0069930069930069930070);
    DEFINVERSE( 144, 0.0069444444444444444444);
    DEFINVERSE( 147, 0.0068027210884353741497);
    DEFINVERSE( 150, 0.0066666666666666666667);
    DEFINVERSE( 153, 0.0065359477124183006536);
    DEFINVERSE( 154, 0.0064935064935064935065);
    DEFINVERSE( 156, 0.0064102564102564102564);
    DEFINVERSE( 160, 0.00625);
    DEFINVERSE( 161, 0.0062111801242236024845);
    DEFINVERSE( 165, 0.0060606060606060606061);
    DEFINVERSE( 168, 0.0059523809523809523810);
    DEFINVERSE( 169, 0.0059171597633136094675);
    DEFINVERSE( 171, 0.0058479532163742690058);
    DEFINVERSE( 175, 0.0057142857142857142857);
    DEFINVERSE( 176, 0.0056818181818181818182);
    DEFINVERSE( 180, 0.0055555555555555555556);
    DEFINVERSE( 182, 0.0054945054945054945055);
    DEFINVERSE( 187, 0.0053475935828877005348);
    DEFINVERSE( 189, 0.0052910052910052910053);
    DEFINVERSE( 192, 0.0052083333333333333333);
    DEFINVERSE( 195, 0.0051282051282051282051);
    DEFINVERSE( 196, 0.0051020408163265306122);
    DEFINVERSE( 200, 0.005);
    DEFINVERSE( 207, 0.0048309178743961352657);
    DEFINVERSE( 208, 0.0048076923076923076923);
    DEFINVERSE( 209, 0.0047846889952153110048);
    DEFINVERSE( 210, 0.0047619047619047619048);
    DEFINVERSE( 216, 0.0046296296296296296296);
    DEFINVERSE( 224, 0.0044642857142857142857);
    DEFINVERSE( 220, 0.0045454545454545454545);
    DEFINVERSE( 221, 0.0045248868778280542986);
    DEFINVERSE( 225, 0.0044444444444444444444);
    DEFINVERSE( 231, 0.0043290043290043290043);
    DEFINVERSE( 240, 0.0041666666666666666667);
    DEFINVERSE( 247, 0.0040485829959514170040);
    DEFINVERSE( 250, 0.004);
    DEFINVERSE( 252, 0.0039682539682539682540);
    DEFINVERSE( 255, 0.0039215686274509803922);
    DEFINVERSE( 256, 0.00390625);
    // further powers of two
    DEFINVERSE( 512, 0.001953125);
    DEFINVERSE(1024, 0.0009765625);
    DEFINVERSE(2048, 0.00048828125);
    DEFINVERSE(4096, 0.000244140625);
    DEFINVERSE(8192, 0.0001220703125);
    // some low factorials
    DEFINVERSE( 720, 0.0013888888888888888889);
    DEFINVERSE(5040, 0.00019841269841269841270);
    DEFINVERSE(40320, 0.000024801587301587301587);
#undef DEFINVERSE
    /// meta programming arithmetic with small non-negative integers
    /// \note specialisations for N=0,1,2,3
    template<int N> struct Integer : public IntegerInverse<N> {
      IntegerInverse<N>::Inverse;
      IntegerInverse<N>::Divide;
      IntegerInverse<N>::Ratio;
      /// modulus with 2
      static const int Odd = N & 2;
      /// half rounded down to integer
      static const int Half = N>>1;
      /// largest power of two not greater than N
      static const int LargestPow2 = 2*Integer<Half>::LargestPow2;
      /// logarithm of two rounded down to integer
      static const int Log2 = 1+Integer<Half>::Log2;
      /// multiply scalar
      /// \param[in,out] X scalar, replaced by N*X on output
      template<typename Real>
      static void Multiply(Real&X) { __IncMul(); X *= N; }
      /// product of scalar with N
      /// \param[in] X scalar
      /// \return    X*N
      template<typename Real>
      static Real Product(Real X) { __IncMul(); return X*N; }
    };
    //
    template<> struct Integer<0> {
      static const int Odd = 0;
      static const int Half = 0;
      static const int Factorial = 1;
      template<typename Real> static void Multiply    (Real&X) { X=0; }
      template<typename Real> static Real Product     (Real  ) { return 0; }
      template<typename Real> static void Exponentiate(Real&X) { X=1; }
      template<typename Real> static Real Power       (Real  ) { return 1; }
      template<typename Real> static Real PowerProduct(Real, Real Y) {
	return Y; }
    };
    //
    template<> struct Integer<1> : IntegerInverse<1> {
      static const int Odd = 1;
      static const int Half = 0;
      static const int LargestPow2 = 1;
      static const int Log2 = 0;
      static const int Factorial = 1;
      IntegerInverse<1>::Inverse;
      IntegerInverse<1>::Divide;
      IntegerInverse<1>::Ratio;
      template<typename Real> static void Multiply    (Real& ) {}
      template<typename Real> static Real Product     (Real X) { return X; }
      template<typename Real> static void Exponentiate(Real& ) {}
      template<typename Real> static Real Power       (Real X) { return X; }
      template<typename Real> static Real PowerProduct(Real X, Real Y) {
	__IncMul(); return X*Y; }
    };
    //
    template<> struct Integer<2> : IntegerInverse<2> {
      static const int Odd = 0;
      static const int Half = 1;
      static const int LargestPow2 = 2;
      static const int Log2 = 1;
      static const int Factorial = 2;
      IntegerInverse<2>::Inverse;
      IntegerInverse<2>::Divide;
      IntegerInverse<2>::Ratio;
      template<typename Real> static void Multiply    (Real&X) { __IncAdd(); X+=X; }
      template<typename Real> static Real Product     (Real X) { __IncAdd(); return X+X; }
      template<typename Real> static void Exponentiate(Real&X) { __IncMul(); X*=X; }
      template<typename Real> static Real Power       (Real X) { __IncMul(); return X*X; }
    };
    //
    template<> struct Integer<3> : IntegerInverse<3> {
      static const int Odd = 1;
      static const int Half = 1;
      static const int LargestPow2 = 2;
      static const int Log2 = 1;
      static const int Factorial = 6;
      IntegerInverse<3>::Inverse;
      IntegerInverse<3>::Divide;
      IntegerInverse<3>::Ratio;
      template<typename Real> static void Multiply    (Real&X) { __CountMA(0,2); X+=X+X; }
      template<typename Real> static Real Product     (Real X) { __CountMA(0,2); return X+X+X; }
      template<typename Real> static void Exponentiate(Real&X) { __CountMA(2,0); X*=X*X; }
      template<typename Real> static Real Power       (Real X) { __CountMA(2,0); return X*X*X; }
    };
    /// Ratio with integer: convert into real-valued multiplication
    template<int N, typename Real> inline
    Real Over(Real x) { return Integer<N>::Ratio(x); }
    /// Product with integer: convert to sum for N=0,1,2,3
    template<int N, typename Real> inline
    Real Times(Real x) { return Integer<N>::Product(x); }
    ////////////////////////////////////////////////////////////////////////////
    // add or subtract, depending on whether template parameter is even or odd
    template<int> struct __AddSubOdd;
    template<> struct __AddSubOdd<0> {
      template<typename Real> static void neg (Real&) {}
      template<typename Real> static void ass (Real&x, Real y) { x =y; }
      template<typename Real> static void add (Real&x, Real y) { __IncAdd(); x+=y; }
      template<typename Real> static void sub (Real&x, Real y) { __IncAdd(); x-=y; }
      template<typename Real> static void add2(Real&x, Real y) { __CountMA(0,2); x+=y; x+=y; }
      template<typename Real> static void sub2(Real&x, Real y) { __CountMA(0,2); x-=y; x-=y; }
      template<typename Real> static Real sum (Real x, Real y) { __IncAdd(); return x+y; }
      template<typename Real> static Real diff(Real x, Real y) { __IncAdd(); return x-y; }
    };
    template<> struct __AddSubOdd<1> {
      template<typename Real> static void neg (Real&x) { x=-x; }
      template<typename Real> static void ass (Real&x, Real y) { x =-y; }
      template<typename Real> static void add (Real&x, Real y) { __IncAdd(); x-=y; }
      template<typename Real> static void sub (Real&x, Real y) { __IncAdd(); x+=y; }
      template<typename Real> static void add2(Real&x, Real y) { __CountMA(0,2); x-=y; x-=y; }
      template<typename Real> static void sub2(Real&x, Real y) { __CountMA(0,2); x+=y; x+=y; }
      template<typename Real> static Real sum (Real x, Real y) { __IncAdd(); return x-y; }
      template<typename Real> static Real diff(Real x, Real y) { __IncAdd(); return x+y; }
    };
    /// assigns positive or negative, depending on whether N is even or odd
    /// \param[in,out] x lvalue, replaced by y if N is even, -y if N is odd
    /// \param[in]     y rvalue
    template<int N, typename Real> inline
    void Assign(Real&x, Real y) { __AddSubOdd<N&1>::ass(x,y); }
    /// switches sign of x if N is odd
    /// \param[in,out] x lvalue, replaced by -x if N is odd
    template<int N, typename Real> inline
    void Negate(Real&x) { __AddSubOdd<N&1>::neg(x); }
    /// adds or subtracts, depending on whether N is even or odd
    /// \param[in,out] x lvalue, replaced by x+y if N is even, x-y if N is odd
    /// \param[in]     y rvalue
    template<int N, typename Real> inline
    void Add(Real&x, Real y) { __AddSubOdd<N&1>::add(x,y); }
    /// adds or subtracts twice, depending on whether N is even or odd
    /// \param[in,out] x lvalue, replaced by x+2y if N is even, x-2y if N is odd
    /// \param[in]     y rvalue
    template<int N, typename Real> inline
    void AddTwice(Real&x, Real y) { __AddSubOdd<N&1>::add2(x,y); }
    /// subtracts or adds, depending on whether N is even or odd
    /// \param[in,out] x lvalue, replaced by x-y if N is even, x+y if N is odd
    /// \param[in]     y rvalue
    template<int N, typename Real> inline
    void Sub(Real&x, Real y) { __AddSubOdd<N&1>::sub(x,y); }
    /// subtracts or adds twice, depending on whether N is even or odd
    /// \param[in,out] x lvalue, replaced by x-2y if N is even, x+2y if N is odd
    /// \param[in]     y rvalue
    template<int N, typename Real> inline
    void SubTwice(Real&x, Real y) { __AddSubOdd<N&1>::sub2(x,y); }
    /// sum or difference, depending on whether N is even or odd
    /// \param[in] x rvalue
    /// \param[in] y rvalue
    /// \return x+y if N is even, x-y if N is odd
    template<int N, typename Real> inline
    Real Sum(Real x, Real y) { return __AddSubOdd<N&1>::sum(x,y); }
    /// difference or sum, depending on whether N is even or odd
    /// \param[in] x rvalue
    /// \param[in] y rvalue
    /// \return x-y if N is even, x+y if N is odd
    template<int N, typename Real> inline
    Real Diff(Real x, Real y) { return __AddSubOdd<N&1>::diff(x,y); }
    ////////////////////////////////////////////////////////////////////////////
    /// sign of an integer
    template<int N> struct Sign {
      const static int S = N<0? -1 : N>0? 1 : 0;
    };
    ////////////////////////////////////////////////////////////////////////////
  }
}

#endif

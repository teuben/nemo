/*
 * N2H+ 93.17 GHz triplet fitting 7 component fitting
 * 
 *
 * 93171.619 (10) N2H+         1-0 F1=1-1 F=0-1              0.5   L134N           NRAO     11m Sny79  Caz86
 * 93171.947 (10) N2H+         1-0 F1=1-1 F=2-2              0.7   L134N           NRAO     11m Sny77  Caz86
 * 93172.098 (10) N2H+         1-0 F1=1-1 F=1-0              0.8   L134N           NRAO     11m Sny77  Caz86
 * 93173.505 (10) N2H+         1-0 F1=2-1 F=2-1              0.9   L134N           NRAO     11m Sny77  Caz86
 * 93173.809 (10) N2H+         1-0 F1=2-1 F=3-2              0.9   L134N           NRAO     11m Sny77  Caz86
 * 93174.016 (10) N2H+         1-0 F1=2-1 F=1-1              0.6   L134N           NRAO     11m Sny77  Caz86
 * 93176.310 (10) N2H+         1-0 F1=0-1 F=1-2              0.7   L134N           NRAO     11m Sny77  Caz86
 *
 * f=[93171.619,93171.947,93172.098,93173.505,93173.809,93174.016,93176.310]
 * v=[    6.068,    5.013,    4.527,    0.0,     -0.978,   -1.644,   -9.025]
 * a=[    0.5  ,    0.7,      0.8,      0.9,      0.9,      0.6,      0.7]
 *
 * splatalogue:
 *
 * / 93171.6103, 93171.9051, 93172.0423, 93173.4677, 93173.7642, 93173.9587, 93176.2543 /
 * / 0.1429, 0.7143, 0.4286, 0.7143, 1.0, 0.4286, 0.4286 /
 *   5.976 5.028 4.586 0.0 -0.954 -1.580 -8.966
                            ^^^^^
                 todo:    center one
 *
 *
 *
 * gauss7_n2h+: (22 parameters)
 *	f = a + b1*exp(-(x-c1)^2/(2*d1^2)) + 
 *              b2*exp(-(x-c2)^2/(2*d2^2))
 *              b3*exp(-(x-c3)^2/(2*d3^2))
 *              b4*exp(-(x-c4)^2/(2*d4^2))
 *              b5*exp(-(x-c5)^2/(2*d5^2))
 *              b6*exp(-(x-c6)^2/(2*d6^2))
 *              b7*exp(-(x-c7)^2/(2*d7^2))
 *  
 * gauss7_n2h+_1: (16 parameters)
 *   with fixing the delta V's between the 7 hyperfines
 *        par=a,c5,b1,d1,.....
 *
 * gauss7_n2h+_2: (10 parameters)
 *   with fixing the line widths to be all the same
 *        par=a,c5,d,b1,b2,b3,...
 *
 * 28-dec-2011	created for N1333
 */

/*  delta wings in km/s */

#define F01  6.068
#define F22  5.013
#define F10  4.527
#define F21  0.000
#define F32 -0.978
#define F11 -1.644
#define F12 -9.025

#include <stdinc.h>


static real dv1 = F01;
static real dv2 = F22;
static real dv3 = F10;
static real dv4 = F21;
static real dv5 = F32;
static real dv6 = F11;
static real dv7 = F12;

static debug_first = 1;


real func_gauss7_n2hp(real *x, real *p, int np)
{
  real a,b,arg1,arg2,arg3,arg4,arg5,arg6,arg7;

  if (debug_first) {
    dprintf(0,"gauss7_n2h+: %g %g %g %g %g %g %g\n",
	    dv1,dv2,dv3,dv4,dv5,dv6,dv7);
    debug_first = 0;
  }

  a = p[2]-x[0];
  b = p[3];
  arg1 = a*a/(2*b*b);

  a = p[5]-x[0];
  b = p[6];
  arg2 = a*a/(2*b*b);

  a = p[8]-x[0];
  b = p[9];
  arg3 = a*a/(2*b*b);

  a = p[11]-x[0];
  b = p[12];
  arg4 = a*a/(2*b*b);

  a = p[14]-x[0];
  b = p[15];
  arg5 = a*a/(2*b*b);

  a = p[17]-x[0];
  b = p[18];
  arg6 = a*a/(2*b*b);

  a = p[20]-x[0];
  b = p[21];
  arg7 = a*a/(2*b*b);

  return p[0] + p[1]*exp(-arg1) + p[4]*exp(-arg2) + p[7]*exp(-arg3) +
    p[10]*exp(-arg4) + p[13]*exp(-arg5) + p[16]*exp(-arg6) + p[19]*exp(-arg7);

}
 
void derv_gauss7_n2hp(real *x, real *p, real *e, int np)
{
  real a1,b1,arg1,a2,b2,arg2,a3,b3,arg3,a4,b4,arg4,a5,b5,arg5,a6,b6,arg6,a7,b7,arg7;

  a1 = p[2]-x[0];
  b1 = p[3];
  arg1 = a1*a1/(2*b1*b1);

  a2 = p[5]-x[0];
  b2 = p[6];
  arg2 = a2*a2/(2*b2*b2);

  a3 = p[8]-x[0];
  b3 = p[9];
  arg3 = a3*a3/(2*b3*b3);

  a4 = p[11]-x[0];
  b4 = p[12];
  arg4 = a3*a3/(2*b3*b3);

  a5 = p[14]-x[0];
  b5 = p[15];
  arg5 = a3*a3/(2*b3*b3);

  a6 = p[17]-x[0];
  b6 = p[18];
  arg6 = a3*a3/(2*b3*b3);

  a7 = p[20]-x[0];
  b7 = p[21];
  arg7 = a3*a3/(2*b3*b3);
  
  e[0] = 1.0;

  e[1] = exp(-arg1);
  e[2] = -p[1]*e[1] * a1 / (b1*b1);
  e[3] = p[1] * e[1] * a1 * a1 / (b1*b1*b1);

  e[4] = exp(-arg2);
  e[5] = -p[4]*e[4] * a2 / (b2*b2);
  e[6] = p[4] * e[4] * a2 * a2 / (b2*b2*b2);

  e[7] = exp(-arg3);
  e[8] = -p[7]*e[7] * a3 / (b3*b3);
  e[9] = p[7] * e[7] * a3 * a3 / (b3*b3*b3);

  e[10] = exp(-arg4);
  e[11] = -p[10]*e[10] * a4 / (b4*b4);
  e[12] = p[10] * e[10] * a4 * a4 / (b4*b4*b4);

  e[13] = exp(-arg5);
  e[14] = -p[13]*e[13] * a5 / (b5*b5);
  e[15] = p[13] * e[13] * a5 * a5 / (b5*b5*b5);

  e[16] = exp(-arg6);
  e[17] = -p[16]*e[16] * a6 / (b6*b6);
  e[18] = p[16] * e[16] * a6 * a6 / (b6*b6*b6);

  e[19] = exp(-arg7);
  e[20] = -p[19]*e[19] * a7 / (b7*b7);
  e[21] = p[19] * e[19] * a7 * a7 / (b7*b7*b7);

}

/* parameters:  fix the velocity offsets
 *      a,c,b1,d1,b2,d2,b3,d3,b4,d4,b5,d5,b6,d6,b7,d7
 *      0 1 2  3  4  5  6  7  8  9  10 11 12 13 14 15
 */


real func_gauss7_n2hp_1(real *x, real *p, int np)
{
  real a,b,arg1,arg2,arg3,arg4,arg5,arg6,arg7;
  a = p[1]-x[0]+dv1;
  b = p[3];
  arg1 = a*a/(2*b*b);

  a = p[1]-x[0]+dv2;
  b = p[5];
  arg2 = a*a/(2*b*b);

  a = p[1]-x[0]+dv3;
  b = p[7];
  arg3 = a*a/(2*b*b);

  a = p[1]-x[0]+dv4;
  b = p[9];
  arg4 = a*a/(2*b*b);

  a = p[1]-x[0]+dv5;
  b = p[11];
  arg5 = a*a/(2*b*b);

  a = p[1]-x[0]+dv6;
  b = p[13];
  arg6 = a*a/(2*b*b);

  a = p[1]-x[0]+dv7;
  b = p[15];
  arg7 = a*a/(2*b*b);

  return p[0] + p[2]*exp(-arg1) + p[4]*exp(-arg2) + p[6]*exp(-arg3) +
    p[8]*exp(-arg4) + p[10]*exp(-arg5) + p[12]*exp(-arg6) + p[14]*exp(-arg7);

}

void derv_gauss7_n2hp_1(real *x, real *p, real *e, int np)
{
  real a1,b1,arg1,a2,b2,arg2,a3,b3,arg3,a4,b4,arg4,a5,b5,arg5,a6,b6,arg6,a7,b7,arg7;

  a1 = p[1]-x[0]+dv1;
  b1 = p[3];
  arg1 = a1*a1/(2*b1*b1);

  a2 = p[1]-x[0]+dv2;
  b2 = p[5];
  arg2 = a2*a2/(2*b2*b2);

  a3 = p[1]-x[0]+dv3;
  b3 = p[7];
  arg3 = a3*a3/(2*b3*b3);

  a4 = p[1]-x[0]+dv4;
  b4 = p[9];
  arg4 = a4*a4/(2*b4*b4);

  a5 = p[1]-x[0]+dv5;
  b5 = p[11];
  arg5 = a5*a5/(2*b5*b5);

  a6 = p[1]-x[0]+dv6;
  b6 = p[13];
  arg6 = a6*a6/(2*b6*b6);

  a7 = p[1]-x[0]+dv7;
  b7 = p[15];
  arg7 = a7*a7/(2*b7*b7);

  
  e[0] = 1.0;

  e[2]  = exp(-arg1);
  e[4]  = exp(-arg2);
  e[6]  = exp(-arg3);
  e[8]  = exp(-arg4);
  e[10] = exp(-arg5);
  e[12] = exp(-arg6);
  e[14] = exp(-arg7);

  e[1] = -p[2]*e[2]   * a1 / (b1*b1)
         -p[4]*e[4]   * a2 / (b2*b2)
         -p[6]*e[6]   * a3 / (b3*b3)
         -p[8]*e[8]   * a4 / (b4*b4)
         -p[10]*e[10] * a5 / (b5*b5)
         -p[12]*e[12] * a6 / (b6*b6)
         -p[14]*e[14] * a7 / (b7*b7);

  e[3]  = p[2]  * e[2]  * a1 * a1 / (b1*b1*b1);
  e[5]  = p[4]  * e[4]  * a2 * a2 / (b2*b2*b2);
  e[7]  = p[6]  * e[6]  * a3 * a3 / (b3*b3*b3);
  e[9]  = p[8]  * e[8]  * a4 * a4 / (b4*b4*b4);
  e[11] = p[10] * e[10] * a5 * a5 / (b5*b5*b5);
  e[13] = p[12] * e[12] * a6 * a6 / (b6*b6*b6);
  e[15] = p[14] * e[14] * a7 * a7 / (b7*b7*b7);
}


/* parameters: fix all widths to be the same as well as DV's.
 *      a,c,d,b1,b2,b3,b4,b5,b6,b7
 *      0 1 2 3  4  5  6  7  8  9
 */
real func_gauss7_n2hp_2(real *x, real *p, int np)
{
  real a,b,arg1,arg2,arg3,arg4,arg5,arg6,arg7;
  a = p[1]-x[0]+dv1;
  b = p[2];
  arg1 = a*a/(2*b*b);

  a = p[1]-x[0]+dv2;
  b = p[2];
  arg2 = a*a/(2*b*b);

  a = p[1]-x[0]+dv3;
  b = p[2];
  arg3 = a*a/(2*b*b);

  a = p[1]-x[0]+dv4;
  b = p[2];
  arg4 = a*a/(2*b*b);

  a = p[1]-x[0]+dv5;
  b = p[2];
  arg5 = a*a/(2*b*b);

  a = p[1]-x[0]+dv6;
  b = p[2];
  arg6 = a*a/(2*b*b);

  a = p[1]-x[0]+dv7;
  b = p[2];
  arg7 = a*a/(2*b*b);

  return p[0] + p[3]*exp(-arg1) + p[4]*exp(-arg2) + p[5]*exp(-arg3) +
    p[6]*exp(-arg4) + p[7]*exp(-arg5) + p[8]*exp(-arg6) + p[9]*exp(-arg7);

}

/* parameters:
 *      a,c,d,b1,b2,b3
 *      0 1 2 3  4  5
 */

void derv_gauss7_n2hp_2(real *x, real *p, real *e, int np)
{
  real a1,b1,arg1,a2,b2,arg2,a3,b3,arg3,a4,b4,arg4,a5,b5,arg5,a6,b6,arg6,a7,b7,arg7;

  a1 = p[1]-x[0]+dv1;
  b1 = p[2];
  arg1 = a1*a1/(2*b1*b1);

  a2 = p[1]-x[0]+dv2;
  b2 = p[2];
  arg2 = a2*a2/(2*b2*b2);

  a3 = p[1]-x[0]+dv3;
  b3 = p[2];
  arg3 = a3*a3/(2*b3*b3);

  a4 = p[1]-x[0]+dv4;
  b4 = p[2];
  arg4 = a4*a4/(2*b4*b4);

  a5 = p[1]-x[0]+dv5;
  b5 = p[2];
  arg5 = a5*a5/(2*b5*b5);

  a6 = p[1]-x[0]+dv6;
  b6 = p[2];
  arg6 = a6*a6/(2*b6*b6);

  a7 = p[1]-x[0]+dv7;
  b7 = p[2];
  arg7 = a7*a7/(2*b7*b7);
  
  e[0] = 1.0;

  e[3] = exp(-arg1);
  e[4] = exp(-arg2);
  e[5] = exp(-arg3);
  e[6] = exp(-arg4);
  e[7] = exp(-arg5);
  e[8] = exp(-arg6);
  e[9] = exp(-arg7);

  e[1] = -p[3]*e[3] * a1 / (b1*b1)
         -p[4]*e[4] * a2 / (b2*b2)
         -p[5]*e[5] * a3 / (b3*b3)
         -p[6]*e[6] * a4 / (b4*b4)
         -p[7]*e[7] * a5 / (b5*b5)
         -p[8]*e[8] * a6 / (b6*b6)
         -p[9]*e[9] * a7 / (b7*b7);

  e[2] = p[3] * e[3] * a1 * a1 / (b1*b1*b1) +
         p[4] * e[4] * a2 * a2 / (b2*b2*b2) +
         p[5] * e[5] * a3 * a3 / (b3*b3*b3) +
         p[6] * e[6] * a4 * a4 / (b4*b4*b4) +
         p[7] * e[7] * a5 * a5 / (b5*b5*b5) +
         p[8] * e[8] * a6 * a6 / (b6*b6*b6) +
         p[9] * e[9] * a7 * a7 / (b7*b7*b7);
}

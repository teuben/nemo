/*
 * N2H+ 93.17 GHz 3-3-1 triplet fitting (7) component fitting
 * 
 *
 * 93171.619 (10) N2H+    1-0 F1=1-1 F=0-1   0.5   L134N      NRAO     11m Sny79  Caz86
 * 93171.947 (10) N2H+    1-0 F1=1-1 F=2-2   0.7   L134N      NRAO     11m Sny77  Caz86
 * 93172.098 (10) N2H+    1-0 F1=1-1 F=1-0   0.8   L134N      NRAO     11m Sny77  Caz86
 * 93173.505 (10) N2H+    1-0 F1=2-1 F=2-1   0.9   L134N      NRAO     11m Sny77  Caz86
 * 93173.809 (10) N2H+    1-0 F1=2-1 F=3-2   0.9   L134N      NRAO     11m Sny77  Caz86 *center*
 * 93174.016 (10) N2H+    1-0 F1=2-1 F=1-1   0.6   L134N      NRAO     11m Sny77  Caz86
 * 93176.310 (10) N2H+    1-0 F1=0-1 F=1-2   0.7   L134N      NRAO     11m Sny77  Caz86
 *
 * f=[93171.619,93171.947,93172.098,93173.505,93173.809,93174.016,93176.310]
 * v=[    7.046,    5.991,    5.505,    0.978,    0.000,   -0.666,   -8.047]
 * a=[    0.5  ,    0.7,      0.8,      0.9,      0.9,      0.6,      0.7]
 *                                               center

 *   93171.62,  93171.92,  93172.05,  93173.48,  93173.78,  93173.97,  93176.27    (shaye,splat)
         6.950      5.985      5.566      0.965      0.000     -0.611     -8.012
 * splatalogue:
 * / 93171.6103, 93171.9051, 93172.0423, 93173.4677, 93173.7642, 93173.9587, 93176.2543 /
 * / 0.1429, 0.7143, 0.4286, 0.7143, 1.0,   0.4286, 0.4286 /
     6.930   5.982   5.540   0.954   0.000 -0.626  -8.012
                                     ^^^^^
 *
 *
 *
 * gauss7_n2h+: (22 parameters)   [see also gaussn.c]
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
 * gauss7_n2h+_3: (4 parameters)
 *   with fixing linewidths, delta V's and relative amps
 *        par=a,b5,c5,d
 *
 * 28-dec-2011	created for N1333
 * 10-oct-2013  redefined some functions, gauss7_n2hp is now the 1+N*3 parameter 7-let with
 *              equal sigma's and fixed delta-v's. 
 * 11-oct-2013  remove the _22, since gaussn will do this.    Order lines by increasing velocity
 */

/*  delta wings in km/s */

#define F01  6.950
#define F22  5.985 
#define F10  5.566
#define F21  0.965
#define F32  0.000
#define F11 -0.611
#define F12 -8.012

#if 0
  /* lovas */
# define A01  0.556
# define A22  0.778 
# define A10  0.889
# define A21  1.0
# define A32  1.0
# define A11  0.667
# define A12  0.778
#else
  /* splatalogue */
# define A01  0.1429
# define A22  0.7143
# define A10  0.4286
# define A21  0.7143
# define A32  1.0
# define A11  0.4286
# define A12  0.4286
#endif

#include <stdinc.h>

/* lines are ordered from low to high vel */

static real dv1 = F12;
static real dv2 = F11;
static real dv3 = F32;  /* center */
static real dv4 = F21;
static real dv5 = F10;
static real dv6 = F22;
static real dv7 = F01;

static real av1 = A12;
static real av2 = A11;
static real av3 = A32;  /* center */
static real av4 = A21;
static real av5 = A10;
static real av6 = A22;
static real av7 = A01;


static debug_first = 1;


/* parameters: fix all widths to be the same as well as DV's as well as all Amp's
 *             to a pre-described ratio
 * 
 *      a,b5,c5,d
 *      0 1  2  3 
 *
 */
real func_gauss7_n2hp(real *x, real *p, int np)
{
  real a,b,arg1,arg2,arg3,arg4,arg5,arg6,arg7,val;
  int ic,nc, ia,iv,id;


  if ( (np-1)%3) error("Need 1+3*Ncomponents parameters: %d",(np-1)%3);
  nc = (np-1)/3;

  if (debug_first) {
    dprintf(0,"gauss7_n2hp: np=%d nc=%d\n",np,nc);
    dprintf(0,"gauss7_n2hp: a,b5,c5,d: %g %g %g %g\n",p[0],p[1],p[2],p[3]);
    dprintf(0,"gauss7_n2hp: dv=%g,%g,%g,%g,%g,%g,%g\n",
	    dv1,dv2,dv3,dv4,dv5,dv6,dv7);
    dprintf(0,"gauss7_n2hp: av=%g,%g,%g,%g,%g,%g,%g\n",
	    av1,av2,av3,av4,av5,av6,av7);
    debug_first = 0;
  }

  val = p[0];
  for (ic=0; ic<nc; ic++) {
    ia = 1+ic*3;            /* offsets into parameter array */
    iv = 2+ic*3;
    id = 3+ic*3;
    
    a = p[iv]-x[0]+dv1;
    b = p[id];
    arg1 = a*a/(2*b*b);

    a = p[iv]-x[0]+dv2;
    b = p[id];
    arg2 = a*a/(2*b*b);

    a = p[iv]-x[0]+dv3;
    b = p[id];
    arg3 = a*a/(2*b*b);

    a = p[iv]-x[0]+dv4;
    b = p[id];
    arg4 = a*a/(2*b*b);
    
    a = p[iv]-x[0]+dv5;
    b = p[id];
    arg5 = a*a/(2*b*b);
    
    a = p[iv]-x[0]+dv6;
    b = p[id];
    arg6 = a*a/(2*b*b);
    
    a = p[iv]-x[0]+dv7;
    b = p[id];
    arg7 = a*a/(2*b*b);
    
    val += p[ia]*(av1*exp(-arg1)+av2*exp(-arg2)+av3*exp(-arg3)+
                  av4*exp(-arg4)+av5*exp(-arg5)+av6*exp(-arg6)+av7*exp(-arg7));
  }
  return val;

}

/* parameters:
 *      a,c,d,b1,b2,b3
 *      0 1 2 3  4  5
 */

void derv_gauss7_n2hp(real *x, real *p, real *e, int np)
{
  real a1,b1,arg1,a2,b2,arg2,a3,b3,arg3,a4,b4,arg4,a5,b5,arg5,a6,b6,arg6,a7,b7,arg7;
  real earg1,earg2,earg3,earg4,earg5,earg6,earg7;
  int ic,nc, ia,iv,id;

  if ( (np-1)%3) error("Need 1+3*Ncomponents parameters: %d",(np-1)%3);
  nc = (np-1)/3;

  e[0] = 1.0;
  for (ic=1; ic<np; ic++) e[ic] = 0.0;

  for (ic=0; ic<nc; ic++) {
    ia = 1+ic*3;            /* offsets into parameter array */
    iv = 2+ic*3;
    id = 3+ic*3;    

    a1 = p[iv]-x[0]+dv1;
    b1 = p[id];
    arg1 = a1*a1/(2*b1*b1);
    earg1 = exp(-arg1);

    a2 = p[iv]-x[0]+dv2;
    b2 = p[id];
    arg2 = a2*a2/(2*b2*b2);
    earg2 = exp(-arg2);
    
    a3 = p[iv]-x[0]+dv3;
    b3 = p[id];
    arg3 = a3*a3/(2*b3*b3);
    earg3 = exp(-arg3);
    
    a4 = p[iv]-x[0]+dv4;
    b4 = p[id];
    arg4 = a4*a4/(2*b4*b4);
    earg4 = exp(-arg4);
    
    a5 = p[iv]-x[0]+dv5;
    b5 = p[id];
    arg5 = a5*a5/(2*b5*b5);
    earg5 = exp(-arg5);
    
    a6 = p[iv]-x[0]+dv6;
    b6 = p[id];
    arg6 = a6*a6/(2*b6*b6);
    earg6 = exp(-arg6);
    
    a7 = p[iv]-x[0]+dv7;
    b7 = p[id];
    arg7 = a7*a7/(2*b7*b7);
    earg7 = exp(-arg7);
      
    e[ia] +=  av1*earg1 + av2*earg2 + av3*earg3 +
              av4*earg4 + av5*earg5 + av6*earg6 + av7*earg7;


    e[iv] += -p[ia]*av1*earg1*a1/(b1*b1)
             -p[ia]*av2*earg2*a2/(b2*b2)
             -p[ia]*av3*earg3*a3/(b3*b3)
             -p[ia]*av4*earg4*a4/(b4*b4)
             -p[ia]*av5*earg5*a5/(b5*b5)
             -p[ia]*av6*earg6*a6/(b6*b6)
             -p[ia]*av7*earg7*a7/(b7*b7);


    e[id] += p[ia] * av1*earg1 * a1*a1 / (b1*b1*b1) +
             p[ia] * av2*earg2 * a2*a2 / (b2*b2*b2) +
             p[ia] * av3*earg3 * a3*a3 / (b3*b3*b3) +
             p[ia] * av4*earg4 * a4*a4 / (b4*b4*b4) +
             p[ia] * av5*earg5 * a5*a5 / (b5*b5*b5) +
             p[ia] * av6*earg6 * a6*a6 / (b6*b6*b6) +
             p[ia] * av7*earg7 * a7*a7 / (b7*b7*b7);
  }
}


/* parameters:  fix the velocity offsets
 *      a,c,b1,d1,b2,d2,b3,d3,b4,d4,b5,d5,b6,d6,b7,d7
 *      0 1 2  3  4  5  6  7  8  9  10 11 12 13 14 15
 */


real func_gauss7_n2hp_16(real *x, real *p, int np)
{
  real a,b,arg1,arg2,arg3,arg4,arg5,arg6,arg7;

  if (debug_first) {
    dprintf(0,"gauss7_n2hp_16: a,c,b1,d1,b2,d2,...b7,d7\n");
    dprintf(0,"gauss7_n2hp_16: %g,%g,%g,%g,%g,%g,%g\n",
	    dv1,dv2,dv3,dv4,dv5,dv6,dv7);
    debug_first = 0;
  }
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

void derv_gauss7_n2hp_16(real *x, real *p, real *e, int np)
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
real func_gauss7_n2hp_10(real *x, real *p, int np)
{
  real a,b,arg1,arg2,arg3,arg4,arg5,arg6,arg7;

  if (debug_first) {
    dprintf(0,"gauss7_n2hp_10: a,c,d,b1,b2,...b7: %g %g %g %g %g %g %g %g %g %g\n",
	    p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10]);
    dprintf(0,"gauss7_n2hp_10: dv=%g,%g,%g,%g,%g,%g,%g\n",
	    dv1,dv2,dv3,dv4,dv5,dv6,dv7);
    debug_first = 0;
  }

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

void derv_gauss7_n2hp_10(real *x, real *p, real *e, int np)
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


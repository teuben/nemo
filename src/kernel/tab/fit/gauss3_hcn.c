/*
 * HCN 88.63 GHz triplet fitting 3 (constrained) gaussians
 *
 * 88630.4157(10) HCN          1-0 F=1-1                     9.6   OriMC-1         NRAO     11m Uli76  DeL69
 * 88631.8473(10) HCN          1-0 F=2-1                    17.2   OriMC-1         NRAO     11m Uli76  DeL69
 * 88633.9360(10) HCN          1-0 F=0-1                     6.8   OriMC-1         NRAO     11m Uli76  DeL69
 *
 * f=[88630.4157,88631.8473,88633.9360]
 * v=[    4.842,     0.0,      -7.065 ]
 * / 88630.41560, 88631.84750, 88633.93570 /
 * / 0.6, 1.0, 0.2 /
 * 
 * splatalogue:
 * f=[88.63042,88.63185,88.63394]
 * s=[9.6,17.2,6.8]
 * id=11,21,01
 *
 *
 * gauss3_hcn: (10 parameters)
 *	f = a + b1*exp(-(x-c1)^2/(2*d1^2)) + 
 *              b2*exp(-(x-c2)^2/(2*d2^2))
 *              b3*exp(-(x-c3)^2/(2*d3^2))
 *  
 * gauss3_hcn_1: (8 parameters)
 *   with
 *        c3=c2+4.8457  1-0 F=1-1  middle
 *        c2            1-0 F=2-1  strongest
 *        c1=c2-7.0698  1-0 F=0-1  weakest
 *
 * gauss3_hcn_2: (6 parameters)
 *        d1=d2=d3 
 *
 * 31-mar-2011	created for N1333
 */

/*  delta wings in km/s */

#define F11 +4.8457
#define F21  0.0000
#define F01 -7.0698

#include <stdinc.h>

static real dv1 = F01;
static real dv2 = F21;
static real dv3 = F11;

static debug_first = 1;


real func_gauss3_hcn(real *x, real *p, int np)
{
  real a,b,arg1,arg2,arg3;

  if (debug_first) {
    dprintf(0,"gauss3_hcn: a,b1,c1,d1,b2,c2,d2,b3,c3,d3\n");
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

  return p[0] + p[1]*exp(-arg1) + p[4]*exp(-arg2) + p[7]*exp(-arg3);

}
 
void derv_gauss3_hcn(real *x, real *p, real *e, int np)
{
  real a1,b1,arg1,a2,b2,arg2,a3,b3,arg3;

  a1 = p[2]-x[0];
  b1 = p[3];
  arg1 = a1*a1/(2*b1*b1);

  a2 = p[5]-x[0];
  b2 = p[6];
  arg2 = a2*a2/(2*b2*b2);

  a3 = p[8]-x[0];
  b3 = p[9];
  arg3 = a3*a3/(2*b3*b3);
  
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
}

/* parameters:
 *      a,c,b1,d1,b2,d2,b3,d3
 *      0 1 2  3  4  5  6  7
 */

real func_gauss3_hcn_1(real *x, real *p, int np)
{
  real a,b,arg1,arg2,arg3;

  if (debug_first) {
    dprintf(0,"gauss3_hcn_1: a,c,b1,d1,b2,d2,b3,d3\n");
    dprintf(0,"gauss3_hcn_1: dv=%g,%g,%g\n",dv1,dv2,dv3);
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

  return p[0] + p[2]*exp(-arg1) + p[4]*exp(-arg2) + p[6]*exp(-arg3);

}

void derv_gauss3_hcn_1(real *x, real *p, real *e, int np)
{
  real a1,b1,arg1,a2,b2,arg2,a3,b3,arg3;

  a1 = p[1]-x[0]+dv1;
  b1 = p[3];
  arg1 = a1*a1/(2*b1*b1);

  a2 = p[1]-x[0]+dv2;
  b2 = p[5];
  arg2 = a2*a2/(2*b2*b2);

  a3 = p[1]-x[0]+dv3;
  b3 = p[7];
  arg3 = a3*a3/(2*b3*b3);
  
  e[0] = 1.0;

  e[2] = exp(-arg1);
  e[4] = exp(-arg2);
  e[6] = exp(-arg3);

  e[1] = -p[2]*e[2] * a1 / (b1*b1)
         -p[4]*e[4] * a2 / (b2*b2)
         -p[6]*e[6] * a3 / (b3*b3);

  e[3] = p[2] * e[2] * a1 * a1 / (b1*b1*b1);
  e[5] = p[4] * e[4] * a2 * a2 / (b2*b2*b2);
  e[7] = p[6] * e[6] * a3 * a3 / (b3*b3*b3);
}


/* parameters:
 *      a,c,d,b1,b2,b3
 *      0 1 2 3  4  5
 */
real func_gauss3_hcn_2(real *x, real *p, int np)
{
  real a,b,arg1,arg2,arg3;

  if (debug_first) {
    dprintf(0,"gauss3_hcn_2: a,c,d,b1,b2,b3\n");
    dprintf(0,"gauss3_hcn_2: dv=%g,%g,%g\n",dv1,dv2,dv3);
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

  return p[0] + p[3]*exp(-arg1) + p[4]*exp(-arg2) + p[5]*exp(-arg3);

}

/* parameters:
 *      a,c,d,b1,b2,b3
 *      0 1 2 3  4  5
 */

void derv_gauss3_hcn_2(real *x, real *p, real *e, int np)
{
  real a1,b1,arg1,a2,b2,arg2,a3,b3,arg3;

  a1 = p[1]-x[0]+dv1;
  b1 = p[2];
  arg1 = a1*a1/(2*b1*b1);

  a2 = p[1]-x[0]+dv2;
  b2 = p[2];
  arg2 = a2*a2/(2*b2*b2);

  a3 = p[1]-x[0]+dv3;
  b3 = p[2];
  arg3 = a3*a3/(2*b3*b3);
  
  e[0] = 1.0;

  e[3] = exp(-arg1);
  e[4] = exp(-arg2);
  e[5] = exp(-arg3);

  e[1] = -p[3]*e[3] * a1 / (b1*b1)
         -p[4]*e[4] * a2 / (b2*b2)
         -p[5]*e[5] * a3 / (b3*b3);

  e[2] = p[3] * e[3] * a1 * a1 / (b1*b1*b1) +
         p[4] * e[4] * a2 * a2 / (b2*b2*b2) +
         p[5] * e[5] * a3 * a3 / (b3*b3*b3) ;
}

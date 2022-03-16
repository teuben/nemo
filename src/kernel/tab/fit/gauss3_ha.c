/*
 * NII-a, Ha, NII-b triplet fitting 3 (constrained) gaussians
 * WARNING:  untested code
 *
 * 6549.86    NII-a
 * 6564.614   Halpha
 * 6585.27    NII-b
 *
 * WARNING: some hcn notation from the clone left over.....
 *
 * gauss3_ha_10: (10 parameters - see gaussn, not fitted here anymore
 *	f = a + b1*exp(-(x-c1)^2/(2*d1^2)) + 
 *              b2*exp(-(x-c2)^2/(2*d2^2))
 *              b3*exp(-(x-c3)^2/(2*d3^2))
 *  
 * gauss3_ha_1: (8 parameters)
 *   with
 *        c3=c2+4.8457  1-0 F=1-1  middle
 *        c2            1-0 F=2-1  strongest
 *        c1=c2-7.0698  1-0 F=0-1  weakest
 *
 * gauss3_ha_2: (6 parameters)
 *        d1=d2=d3 
 *
 * gauss3_ha_2: (1+3N parameters, where N is number of lines
 *        d1=d2=d3 
 *
 * 14-mar-2022  cloned off gauss3_hcn
 */

/*  delta wings in AA */

#define F11 +14.754
#define F21  0.0000
#define F01 -20.656

// ??
#define A11 1
#define A21 1
#define A01 1


#include <stdinc.h>

/* the lines are ordered from low to high velocity */
static real dv1 = F01;
static real dv2 = F21;
static real dv3 = F11;


static real av1 = A01;
static real av2 = A21;
static real av3 = A11;


static int debug_first = 1;


/* parameters: a simple NII-Halpha-NII triplet, only continuum, middle amp, middle vel
 *             and dispersion given (4 parameters per component)
 *
 */
real func_gauss3_ha(real *x, real *p, int np)
{
  real a,b,arg1,arg2,arg3,val;
  int ic,nc, ia,iv,id;


  if ( (np-1)%3) error("gauss3_ha: need 1+3*Ncomp parameters: %d",(np-1)%3);
  nc = (np-1)/3;

  if (debug_first) {
    dprintf(1,"gauss3_ha: np=%d nc=%d\n",np,nc);
    dprintf(1,"gauss3_ha: a,b,c,d:");
    for (ic=0; ic<np; ic++) dprintf(1," %g",p[ic]);
    dprintf(1,"\n");
    dprintf(1,"gauss3_ha: dv=%g,%g,%g\n",dv1,dv2,dv3);
    dprintf(1,"gauss3_ha: av=%g,%g,%g\n",av1,av2,av3);
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

    val += p[ia]*(av1*exp(-arg1)+av2*exp(-arg2)+av3*exp(-arg3));
  }
  return val;

}

void derv_gauss3_ha(real *x, real *p, real *e, int np)
{
  real a1,b1,arg1,a2,b2,arg2,a3,b3,arg3;
  real earg1,earg2,earg3;
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
    
    e[ia] +=  av1*earg1 + av2*earg2 + av3*earg3 ;

    e[iv] += -p[ia]*av1*earg1*a1/(b1*b1)
             -p[ia]*av2*earg2*a2/(b2*b2)
             -p[ia]*av3*earg3*a3/(b3*b3);

    e[id] += p[ia] * av1*earg1 * a1*a1 / (b1*b1*b1) +
             p[ia] * av2*earg2 * a2*a2 / (b2*b2*b2) +
             p[ia] * av3*earg3 * a3*a3 / (b3*b3*b3);
  }
}

/*  the remainder of this code is for less restricted fits, where
 *    amplitudes or velocity dispersions are kept free for each wing
 *    for all free, use gaussn
 */

/* parameters:
 *      a,c,b1,d1,b2,d2,b3,d3
 *      0 1 2  3  4  5  6  7
 */

real func_gauss3_ha_1(real *x, real *p, int np)
{
  real a,b,arg1,arg2,arg3;

  if (debug_first) {
    dprintf(0,"gauss3_ha_1: a,c,b1,d1,b2,d2,b3,d3\n");
    dprintf(0,"gauss3_ha_1: dv=%g,%g,%g\n",dv1,dv2,dv3);
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

void derv_gauss3_ha_1(real *x, real *p, real *e, int np)
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
real func_gauss3_ha_2(real *x, real *p, int np)
{
  real a,b,arg1,arg2,arg3;

  if (np!=6) error("Needs 6 parameters");

  if (debug_first) {
    dprintf(0,"gauss3_ha_2: a,c,d,b1,b2,b3\n",
	    p[0],p[1],p[2],p[3],p[4],p[5]);
    dprintf(0,"gauss3_ha_2: dv=%g,%g,%g\n",dv1,dv2,dv3);
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

void derv_gauss3_ha_2(real *x, real *p, real *e, int np)
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

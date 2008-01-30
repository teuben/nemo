
#include <nemo.h>
#include <rotcurshape.h>



/* a bunch of rotation curves and parameter derivatives */


real rotcur_flat(real r, int n, real *p, real *d)
{
  d[0] = 1.0;
  return p[0];
}

real rotcur_linear(real r, int n, real *p, real *d)
{
  d[0] = r;
  return p[0] * r;
}

real rotcur_poly(real r, int np, real *p, real *d)
{
  real v, dp, x = r/p[1];
  int i;

  i = np-1;
  v = 0;
  dp = 0;
  d[1] = p[0] * x;   /* fake placeholder for recursion coming up */
  while (i > 1) {    /* p[0] and p[1] are special, p[2] last one in loop */
    v = v*x + p[i];
    dp = dp*x + i*p[i];
    d[np+1-i] = d[np-i] * x;
    i--;
  }
  v  = x*(1+x*v);
  dp = x*(1+x*dp);

  d[0] = v;
  d[1] = -p[0]*dp/p[1];

  return p[0] * v;
}

/*
 *  v = p0 * x / (1+x)
 */

real rotcur_core1(real r, int n, real *p, real *d)     /* power, with c=1 fixed */
{
  real x = r / p[1];
  d[0] = x/(1+x);
  d[1] = -p[0]*d[0]/(p[1]*(1+x));
  return p[0] * d[0];
}

/*
 *  v = p0 * x / sqrt(1+x^2)
 */

real rotcur_core2(real r, int n, real *p, real *d)     /* power, with c=2 fixed */
{
  real x = r/p[1];
  real y1 = 1+x*x;
  real y2 = sqrt(y1);
  real v = x / y2;

  d[0] = v;
  d[1] = -p[0]*v/(y1*p[1]);
  return p[0] * v;
}

real rotcur_core(real r, int np, real *p, real *d)
{
  real x = r/p[1];
  real c = p[2];
  real q1 = pow(x,c);
  real q = 1+q1;
  real lnx = log(x);
  real lnq = log(q);
  real y = pow(q,1/c);

  d[0] = x / y;
  d[1] = -p[0]*d[0]/(p[1]*q);
  d[2] = (-((q1*lnx)/(c*q)) + lnq/(c*c))/y;     /* CForm[D[(1+x^c)^(-1/c),c]]  */
  d[2] *= p[0] * x;
  return p[0] * d[0];
}

real rotcur_plummer(real r, int np, real *p, real *d)
{
  real x = r/p[1];
  real y = pow(1+x*x,-0.75);
  d[0] = y;
  d[1] = -x*p[0]/p[1]*(1-x*x/2)/(1+x*x)/y;
  return p[0] * x * y;
}

#if 0
real rotcur_tanh(real r, int np, real *p, real *d)
{
  v = tanh(x);
  dvdx = sqr(sech(x));
}
#endif

/*
 * softened iso-thermal sphere: (a.k.a. pseudo-isothermal)
 *    rho  = rho0/(1+x^2)                    x = r/r0
 *    vrot = vrot0*(1-atan(x)/x)^(1/2)   ,   vrot0 = sqrt(4.pi.G.rho0*r0^2)
 *    

In[1]:= Integrate[x^2/(1+x^2),x]

Out[1]= x - ArcTan[x]


In[7]:=D[Sqrt[1-ArcTan[x]/x],x]

              1         ArcTan[x]
        -(----------) + ---------
                  2         2
          x (1 + x )       x
Out[7]= -------------------------
                     ArcTan[x]
          2 Sqrt[1 - ---------]
                         x
*/


real rotcur_iso(real r, int np, real *p, real *d)
{
  real x = r/p[1];
  real v = sqrt(1-atan(x)/x);
  d[0] = v;
  d[1] = -p[0]/p[1]*(1/(1+x*x) - v);
  return p[0] * v;
}

/*
 * used in van Moorsel & Wells, AJ 90, 1038 (1985)
 *
 *    V/Vmax = 1 - e^{-ln{100) R/Rmax}
 *  or as we write:
 *    V = Vmax ( 1 - e^{-R/Rmax} )
 *
 */

real rotcur_exp(real r, int np, real *p, real *d)
{
  real x = r/p[1];
  real y = exp(-x);
  d[0] = 1-y;
  d[1] = -p[0]*y/p[1]*x;
  return p[0] * d[0];
}

/*
 *  NFW profile:  pars = V_200,R_200,c
 *  V_c^2(x)=V_{200}^2 \frac{\ln(1+cx)-cx(1+cx)^{-1}}
 *           {x[\ln(1+c)-c(1+c)^{-1}]}
 *
 *  In[7]:=D[ Sqrt[(Log[1+c*x]-c*x/(1+c*x))/x] ,x]
 */

real rotcur_nfw(real r, int np, real *p, real *d)
{
  real x = r/p[1];
  real c = p[2];
  real cx = x*c;
  real lncx = log(1+cx);
  real a = -cx/(1+cx) + lncx;
  real v2 = a/(log(1+c)-c/(1+c))/x;
  real v=sqrt(v2);
  d[0] = v;
  d[1] = -p[0]/p[1]*x  * (sqr(c/(1+cx)) - a/x)/2/sqrt(a);
  d[2] = 0.0;  /* Note, don't allow derivatives w.r.t. c -- otherwise degenerate */
  return p[0] * d[0];
}

/*
 * Moore et al (1999, MNRAS 310, 1147:
 *
 * rho \propto  r^{-3/2}
 * i.e.
 * v   \propto  r^{1/4}      ->   see 'power' with P3=0.25
 */

real rotcur_moore(real r, int np, real *p, real *d)
{
  return 0.0;
}


/*
 * Brandt, J.C. (1960, ApJ 131, 293)
 * Brandt, J.C. & Scheer, L.S.  1965 AJ 70, 471
 *
 * v = v_0   x /   (1/3 + 3/2*x^n)^(3/2n)
 * 
 * 
 */

real rotcur_brandt(real r, int np, real *p, real *d)
{
  real x = r/p[1];
  real n = p[2];
  real dn = pow(1.0/3.0+1.5*pow(x,n), 1.5/n);
  real v = x/dn;
  d[0] = v;
  d[1] = 0.0;   /* to do !! */
  d[2] = 0.0;   /* to do !! */
  return p[0] * v;
}

/* 
 * simple Power Law rotation curve
 *
 *  v = v_0 x^a
 */

real rotcur_power(real r, int np, real *p, real *d)
{
  real x = r/p[1];
  real a = p[2];
  real v = pow(x,a);
  d[0] = v;
  d[1] = -a*p[0]*v/p[1];
  d[2] = p[0]*v*log(x);
  return p[0] * d[0];
}

/*
 * some kind of toy disk for max disk degeneracy simulations
 *
 */

real rotcur_disk1(real r, int np, real *p, real *d)
{
  /* this #if 0 is needed for the 2.96 compiler on mdk81 :-) */
#if 0
  error("disk1 not implemented yet");
#endif
}

/*
 * gauss.c:  (spherical) gauss potential
 *
 *      19-mar-2021    spring break in covid time... while reading Cappelleri 2002
 *
 */

/*CTEX
 *	{\bf potname=gauss
 *       potpars={\it $\Omega,M,\sigma,p,q$}}
 *
 *  Gauss potential (e.g. Cappellari 2002 from the MGE formalism).  For the spherical case
 *  (cases with P and q < 1 will be considered another time)
 *
 * $$
 *    \rho = { M \over {(\sigma \sqrt{2\pi})}^3  }   e^{- {R^2}\over{2\sigma^2}}
 * $$
 *
 * $$
 *    \Phi = -  {  M  \over R } erf({ {R}\over{\sigma\sqrt{2}}})
 * $$
 */


// benchmark:
// /usr/bin/time potccd ccd1.pot gauss 0,1,1,1,0.1 x=-5:5:0.1 y=-5:5:0.1 z=-5:5:0.1
// #   6.87user 0.01system 0:06.88elapsed 99%CPU
// /usr/bin/time potccd ccd1.den gauss 0,1,1,1,0.1 x=-5:5:0.1 y=-5:5:0.1 z=-5:5:0.1 mode=den dr=0.001 nder=2
// #  48.23user 0.01system 0:48.28elapsed 99%CPU 
  
#include <stdinc.h>
 
local double omega = 0.0;
local double g_mass = 1.0;
local double g_sigma = 1.0;
local double g_q = 1.0;         // discarded for now
local double g_p = 1.0;         // discarded for now
local double accuracy = 1e-7;
local int    g_debug = 0;
local int    need_warning = 1;

local double r2,z2,s2,s3,eps,delta;

#define MAX_STEPS 32

local double romberg(double (*f)(double), double a, double b, size_t max_steps, double acc);		     

local double oblate(double t) {
  return exp(-t*t*(r2+z2/(1-eps*t*t))/(2*s2)) / sqrt(1-eps*t*t);
}

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega   = par[0];
    if (n>1) g_mass  = par[1];
    if (n>2) g_sigma = par[2];
    if (n>3) g_p     = par[3];
    if (n>4) g_q     = par[4];
    if (n>5) warning("gauss: npar=%d only up to 5 parameters accepted",n);

    dprintf (1,"INIPOTENTIAL Gauss: \n");
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  mass, sigma, p, q = %f %f %f %f\n",g_mass,g_sigma,g_p,g_q);
	
    par[0] = omega;
    s2 = sqrt(2);
    s3 = sqrt(PI);
    delta = 1 - g_p*g_p;
    eps = 1 - g_q*g_q;

    if (eps > 0 && g_p < 1) error("potname gauss: triaxial case not implemented");

}


void potential_double(int *ndim,double *pos,double *acc,double *pot,double *time)
{
  double r, tmp, x;

  if (eps>0) {   // oblate case in the z=0 plane
    r2 = pos[0]*pos[0] + pos[1]*pos[1];
    z2 = pos[2]*pos[2];
    *pot = -g_mass / g_sigma * sqrt(2/PI) * romberg(oblate,0.0,1.0,30,accuracy);
    if (need_warning) {
      warning("gauss.c: oblate q=%g case has no forces computed yet, only potential",g_q);
      need_warning = 0;
    }
    acc[0] = acc[1] = acc[2] = 0;
    return;
  }

  r2 = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];

  if (r2 == 0) {
    *pot = -g_mass / g_sigma * sqrt(2/PI);
    acc[0] = acc[1] = acc[2] = 0;
    return;
  }

  r = sqrt(r2);
  x = r/(s2*g_sigma);

  tmp = erf(x) / r;
    
  *pot = -g_mass * tmp;

  tmp = tmp - s2*exp(-x*x)/g_sigma/s3;
  tmp = -tmp * g_mass/r2;
  
  acc[0] = tmp*pos[0];
  acc[1] = tmp*pos[1];
  acc[2] = tmp*pos[2];
}

/* adapted from https://en.wikipedia.org/wiki/Romberg's_method */
 
local void dump_row(size_t i, double *R) {
  printf("R[%2zu] = ", i);
  for (size_t j = 0; j <= i; ++j)
    printf("%f ", R[j]);
  printf("\n");
}
 

local double romberg(double (*f)(double), double a, double b, size_t max_steps, double acc) {
  double R1[MAX_STEPS], R2[MAX_STEPS];  // buffers
  double *Rp = &R1[0], *Rc = &R2[0];    // Rp is previous row, Rc is current row
  double h = (b-a);                     // step size
  Rp[0] = (f(a) + f(b))*h*.5;           // first trapezoidal step

  if (g_debug) dump_row(0, Rp);
  
  for (size_t i = 1; i < max_steps; ++i) {
    h /= 2.;
    double c = 0;
    size_t ep = 1 << (i-1); //2^(n-1)
    for (size_t j = 1; j <= ep; ++j) {
      c += f(a+(2*j-1)*h);
    }
    Rc[0] = h*c + .5*Rp[0]; //R(i,0)

    for (size_t j = 1; j <= i; ++j) {
      double n_k = pow(4, j);
      Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); // compute R(i,j)
    }

    // Dump ith row of R, R[i,i] is the best estimate so far
    if (g_debug) dump_row(i, Rc);

    if (i > 1 && fabs(Rp[i-1]-Rc[i]) < acc) {
      return Rc[i-1];
    }

    // swap Rn and Rc as we only need the last row
    double *rt = Rp;
    Rp = Rc;
    Rc = rt;
  }
  return Rp[max_steps-1]; // return our best guess
}

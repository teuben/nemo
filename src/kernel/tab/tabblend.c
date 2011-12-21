/*
 * TABBLEND: create (blined) lines, grid them,output table

 *
 *      21-dec-2011 V1.0    Created - shortest day of the year
 */

#include <stdinc.h>	
#include <getparam.h>
#include <moment.h>


#define MAXL       100
#define MAXP    100000

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
  "x=0\n		 X coordinates of the line(s)",
  "a=1\n                 Amplitudes",
  "d=1\n                 FWHM of the lines (FWHM=2.355*sigma)",
  "s=-16:16:0.01\n       Points to sample",
  "g=-10:10:1\n          Points to grid",
  "fwhm=0\n              If non-zero, smooth points with this beam",
  "hanning=f\n           Optional hanning (not implemented)",
  "fft=f\n               Use FFT to compute spectrum",
  "VERSION=0.1\n	 21-dec-2011 PJT",
  NULL
};

string usage = "create (blended) lines and grid them";

local real x[MAXL], a[MAXL], d[MAXL], d1[MAXL];
local real s[MAXP], g[MAXP], y[MAXP];
local real fwhm;



nemo_main()
{
  int i, j, nx, na, nd, ns, ng;
  real arg;
  real fac1 = 2.3548;   /* 2sqrt(2ln2): fwhm to sigma factor */
  real fac2 = 1.6651;   /* 2sqrt(ln2)  */

  nx = nemoinpr(getparam("x"),x,MAXL);
  if (nx < 1) error("Error parsing x=%s",getparam("x"));

  na = nemoinpr(getparam("a"),a,MAXL);
  if (na < 1) error("Error parsing a=%s",getparam("a"));
  if (na > nx) error("Too many a=%s",getparam("a"));
  for (i=na; i<nx; i++) a[i] = a[i-1];
  

  nd = nemoinpr(getparam("d"),d,MAXL);
  if (nd < 1) error("Error parsing d=%s",getparam("d"));
  if (nd > nx) error("Too many d=%s",getparam("d"));
  for (i=nd; i<nx; i++) d[i] = d[i-1];
  for (i=0; i<nx; i++) d1[i] = fac2 / d[i];

  ns = nemoinpr(getparam("s"),s,MAXP);
  if (ns < 1) error("Error parsing s=%s",getparam("s"));

  ng = nemoinpr(getparam("g"),g,MAXP);
  if (ng < 1) error("Error parsing g=%s",getparam("g"));

  fwhm = getdparam("fwhm");

  /* write out a header reminding which lines used */

  printf("#  X    A   D\n");
  for (i=0; i<nx; i++) {
    printf("#  %g %g %g\n",x[i], a[i], d[i]);
  }

  for (i=0; i<ns; i++) {     /* loop over all sample points */
    y[i] = 0;
    for (j=0; j<nx; j++) {   /* add in all lines */
      arg = (s[i]-x[j])*d1[j];
      y[i] += a[j]*exp(-arg*arg);
    }
    printf("%g  %g\n",s[i],y[i]);
  }

  
}


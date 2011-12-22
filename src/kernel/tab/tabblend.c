/*
 * TABBLEND: create (blended) lines, grid them, output table(s)
 *
 *
 *      21-dec-2011 V0.1    Created - shortest day of the year
 *
 *  @TODO:   gridding shift
 */

#include <stdinc.h>	
#include <getparam.h>
#include <moment.h>
#include <grid.h>


#define MAXL       100
#define MAXP    100000

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
  "x=0\n		 X coordinates of the line(s)",
  "a=1\n                 Amplitudes",
  "d=1\n                 FWHM of the lines (FWHM=2.355*sigma)",
  "s=-16:16:0.01\n       Points to sample",
  "g=-16:16:1\n          Points to grid",
  "fwhm=0\n              If non-zero, convolve points with this beam",
  "hanning=f\n           Optional hanning",
  "fft=f\n               Use FFT to compute spectrum (not implemented)", 
  "rms=0\n               Add gaussian noise",
  "seed=0\n              seed for random number generator",
  "mode=1\n              0 = output raw   1=output final",
  "VERSION=0.3\n	 22-dec-2011 PJT",
  NULL
};

string usage = "create (blended) lines and grid them";

string cvsid="$Id$";


local real x[MAXL], a[MAXL], d[MAXL], d1[MAXL];
local real s[MAXP], g[MAXP], y[MAXP], z[MAXP], z1[MAXP];
local real fwhm;



nemo_main()
{
  Grid G;
  int i, j, il, ir, nx, na, nd, ns, ng;
  real arg, ds, dg;
  real fac1 = 2.3548;   /* 2sqrt(2ln2): fwhm to sigma factor */
  real fac2 = 1.6651;   /* 2sqrt(ln2)                        */
  bool Qraw;
  bool Qfft = getbparam("fft");
  bool Qhan = getbparam("hanning");
  int mode = getiparam("mode");
  int seed = init_xrandom(getparam("seed"));
  real rms = getdparam("rms");

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
  ds = s[1]-s[0];  /* better be a uniform grid */

  ng = nemoinpr(getparam("g"),g,MAXP);
  if (ng < 1) error("Error parsing g=%s",getparam("g"));
  dg = g[1]-g[0];  /* better be a uniform grid */
  inil_grid(&G, ng, g[0], g[ng-1]);

  fwhm = getdparam("fwhm");

  /* write out a header reminding which lines used */

  printf("#  X    A   D\n");
  for (i=0; i<nx; i++) {
    printf("#  %g %g %g\n",x[i], a[i], d[i]);
  }
  if (mode==0) printf("# raw spectrum\n");

  /* compute raw sampled spectrum */

  for (i=0; i<ns; i++) {     /* loop over all sample points */
    y[i] = 0;
    for (j=0; j<nx; j++) {   /* add in all lines */
      arg = (s[i]-x[j])*d1[j];
      arg = arg*arg;
      if (arg > 100) continue;   
      y[i] += a[j]*exp(-arg);
    }
    if (mode==0)
      printf("%g  %g\n",s[i],y[i]);
  }

  /* smooth the raw spectrum */

  if (fwhm > 0.0) {
    warning("fwhm not implemented yet");
  }

  /* use an FFT engine ?  */

  if (Qfft) {
    warning("fft not implemented yet");
  }

  /* grid the 's' to 'g' grid */
  if (ds < dg) {
    for (i=0; i<ng; i++) z[i] = z1[i] = 0.0;
    for (i=0; i<ns; i++) {
      /* @TODO: gridding shift ?  */
      j = index_grid(&G, s[i]);
      if (j<0) continue;
      z[j]  += y[i];
      z1[j] += 1.0;
    }
    for (i=0; i<ng; i++) {
      if (z1[i] > 0.0) z[i] /= z1[i];
    }
  } else
    error("No gridding possible?");

  /* Hanning */
  if (Qhan) {
    warning("Hanning");
    for (i=0; i<ng; i++) {
      il = MAX(i-1,0);
      ir = MIN(i+1,ng-1);
      z1[i] = 0.25 * z[il] + 0.5 * z[i] + 0.25 * z[ir];
    }
    for (i=0; i<ng; i++) z[i] = z1[i];
  } 


  /* add noise */
  if (rms > 0.0) {
    for (i=0; i<ng; i++)
      z[i] += grandom(0.0,rms);
  }

  if (mode==1) {
    printf("# final grid\n");
    for (i=0; i<ng; i++)
      printf("%g  %g\n",g[i],z[i]);
  }

  
}


/*
 * POTROT: query a potential to derive a rotation curve
 *
 *	Peter Teuben		7-feb-05       Created
 */

#include <stdinc.h>
#include <getparam.h>
#include <potential.h>
#include <vectmath.h>

string defv[] = {
    "potname=???\n  Name of potential",
    "potpars=\n     Parameters for potential (1st one is pattern speed)",
    "potfile=\n     Any optional data file associated with potential",
    "r=0:2:0.1\n    Radii to sample",
    "p=0\n          Angles to sample",
    "t=0.0\n        Time to test potential at, if relevant",
    "omega=\n       Use this instead of any returned pattern speed",
    "format=%g\n    Format used to print numbers",
    "VERSION=0.1\n  7-feb-05 PJT",
    NULL,
};

string usage = "query a NEMO potential in the XY plane, and derive a rotation curve ";

string cvsid = "$Id$";

#ifndef MAXPT
#define MAXPT 10001
#endif

local potproc_double mypotd;     /* pointer to potential calculator function : double */
local potproc_float  mypotf;     /* pointer to potential calculator function : float */
local void do_potential(bool,int *, double *, double *, double *, double *);

local double report_q(int,double);

void nemo_main(void)
{
  int    i, j, k, nr,np, nstep;
  double pos[3],acc[3],pot, time;
  real   radii[MAXPT],phi[MAXPT],xarr[MAXPT],yarr[MAXPT],zarr[MAXPT];
  double dr, dphi, phi;
  double fourpi = 4*PI;
  double omega, qt, qt0, dq;
  char *fmt, s[20], pfmt[256];
  int niter = getiparam("niter");
  double eps = getdparam("eps");
  
  nr = nemoinpr(getparam("r"), radii, MAXPT);  /* get radii */
  np = nemoinpr(getparam("p"), phi,   MAXPT);  /* get angles */
  
  time = getdparam("t");
  
  mypotd = get_potential_double(getparam("potname"), 
				getparam("potpars"), 
				getparam("potfile"));
  
  if (hasvalue("omega"))
    omega = getdparam("omega");
  else
    omega = get_pattern();
  dprintf(1,"Found pattern speed = %g\n",omega);
  
  fmt = getparam("format");
  strcpy(s,fmt);    /* use format from command line */
  if (strchr(fmt,' ')==NULL && strchr(fmt,',')==NULL)
    strcat(s," ");      /* append separator if none specified */
  
  for (i=0;i<nr;i++) {
    reset_q();
    add_q(radii[i], 0);
    add_q(radii[i], HALF_PI);
    add_q(radii[i], HALF_PI/2.0);
    qt0 = report_q(1,-1.0);
    
    dphi  = HALF_PI/2.0;
    nstep = 2;
    for (j=0; j<niter; j++) {
      for (k=0, phi = dphi/2.0; k<nstep-1; k++, phi += dphi) {
	add_q(radii[i], phi);
      }
      qt0 = qt;
      qt = report_q(1,qt0);
      dq = (qt-qt0)/qt;
      dq = ABS(dq);
      if (dq < eps && j>5) {
	// warning("Early convergence after %d iterations",j);
	break;
      }
      nstep *= 2;
      dphi  /= 2.0;
    }
    qt = report_q(0,qt0);
  }
}


double rad_q;
double sum_fr;
double max_ft;
int    nsum;
int    iter;

reset_q()
{
  sum_fr = max_ft = 0.0;
  nsum = 0;
  iter = 0;
}

add_q (double rad, double phi)
{
  double pos[3], acc[3], pot, r, fr, ft, time = 0.0;
  int ndim = 3;

  pos[0] = rad * cos(phi);
  pos[1] = rad * sin(phi);
  pos[2] = 0.0;

  (*mypotd)(&ndim, pos, acc, &pot, &time);

  r  = pos[0]*pos[0] + pos[1]*pos[1];
  fr = (pos[0]*acc[0] + pos[1]*acc[1])/sqrt(r);
  ft = sqrt( (acc[0]*acc[0] + acc[1]*acc[1]) - fr*fr);
  dprintf(2,"r,fr,ft=%g %g %g\n",r,fr,ft);

  sum_fr += fr;
  nsum++;
  if (ft > max_ft) max_ft = ft;

  rad_q = rad;
  iter++;
  
}

double report_q(int debug, double old)
{
  double x, eps;

  x = nsum * max_ft/sum_fr;
  x = ABS(x);

  if (old > 0) {
    eps = (x-old)/x;
    eps = ABS(eps);
  } else {
    eps = 0.0;
  }

  if (debug)
    dprintf(1,"%g %d %g %g = %g %g\n",rad_q,iter,sum_fr/nsum,max_ft,x,eps);
  else
    printf("%g %g %g\n",rad_q,x,eps);
  return x;
}

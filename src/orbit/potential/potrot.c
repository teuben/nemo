/*
 * POTROT: query a potential to derive a rotation curve
 *
 *	Peter Teuben		7-feb-05       cloned off potq
 */

#include <stdinc.h>
#include <getparam.h>
#include <potential.h>
#include <vectmath.h>
#include <moment.h>

string defv[] = {
    "potname=???\n  Name of potential",
    "potpars=\n     Parameters for potential (1st one is pattern speed)",
    "potfile=\n     Any optional data file associated with potential",
    "r=0:2:0.1\n    Radii to sample",
    "p=0\n          Angles to sample (degrees)",
    "t=0.0\n        Time to test potential at, if relevant",
    "omega=\n       Use this instead of any returned pattern speed",
    "format=%g\n    Format used to print numbers",
    "niter=10\n     Iteration max",
    "eps=0.001\n    Accuracy",
    "VERSION=0.2\n  7-feb-05 PJT",
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
local double force (double rad, double phi);

local double report_q(int,double);

void nemo_main(void)
{
  int    i, j, k, nr,np, nstep;
  double pos[3],acc[3],pot, time;
  real   radii[MAXPT],phis[MAXPT],xarr[MAXPT],yarr[MAXPT],zarr[MAXPT];
  double dr, dphi, phi, fr, f1, f0, df, vrot;
  double fourpi = 4*PI;
  double omega, qt, qt0, dq;
  char *fmt, s[20], pfmt[256];
  int niter = getiparam("niter");
  double eps = getdparam("eps");
  Moment m;
  
  nr = nemoinpr(getparam("r"), radii, MAXPT);  /* get radii */
  np = nemoinpr(getparam("p"), phis,  MAXPT);  /* get angles */
  for (i=0; i<np; i++) phis[i] *= PI/180.0;    /* convert to radian */
  
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
  ini_moment(&m,2,0);

  if (niter>0 && np == 2) {
    warning("Testing an iterative procedure: niter=%d eps=%f",niter,eps);
    for (i=0;i<nr;i++) {
      reset_moment(&m);
      accum_moment(&m, force(radii[i],phis[0]),1.0);
      accum_moment(&m, force(radii[i],phis[1]),1.0);
      accum_moment(&m, force(radii[i],0.5* (phis[0]+phis[1])),1.0);
      f0 = mean_moment(&m);
      dphi  = (phis[1]-phis[0])/2.0;
      nstep = 2;
      for (j=0; j<niter; j++) {
	for (k=0, phi = dphi/2.0; k<nstep-1; k++, phi += dphi) {
	  accum_moment(&m, force(radii[i],phi),1.0);
	}
	f1 = mean_moment(&m);
	df = (f1-f0)/f0;
	df = ABS(df);
	if (df < eps && j>5) {
	  // warning("Early convergence after %d iterations",j);
	  break;
	}
	nstep *= 2;
	dphi  /= 2.0;
      }
      vrot = sqrt(-f1*radii[i]);
      printf("%g %g %g %d\n",radii[i],fr,vrot,j);
    }
  } else {
    for (i=0;i<nr;i++) {
      reset_moment(&m);
      for (j=0; j<np; j++) {
	fr = force(radii[i],phis[j]);
	accum_moment(&m,fr,1.0);
      }
      fr = mean_moment(&m);
      vrot = sqrt(-fr*radii[i]);            /* force better be inward bound ! */
      printf("%g %g %g\n",radii[i],fr,vrot); 
    }
  }
}

/* radial force */

double force (double rad, double phi)
{
  double pos[3], acc[3], pot, r, fr, ft, time = 0.0;
  int ndim = 3;

  pos[0] = rad * cos(phi);
  pos[1] = rad * sin(phi);
  pos[2] = 0.0;

  (*mypotd)(&ndim, pos, acc, &pot, &time);

  r  = pos[0]*pos[0] + pos[1]*pos[1];
  fr = (pos[0]*acc[0] + pos[1]*acc[1])/sqrt(r);

  return fr;
#if 0
  ft = sqrt( (acc[0]*acc[0] + acc[1]*acc[1]) - fr*fr);
  dprintf(2,"r,fr,ft=%g %g %g\n",r,fr,ft);

  sum_fr += fr;
  nsum++;
  if (ft > max_ft) max_ft = ft;

  rad_q = rad;
  iter++;
#endif
}


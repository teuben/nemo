/* 
 * PVCORR:   PV diagram correlation to find multiple weaker lines that
 *           match a template.
 *
 *     29-may-2013  V0.1    for astute, cloned off pvtrace
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
  "in=???\n       Input image file - must be an PV diagram",
  "clip=0\n       Clip value to search for profile in template area",
  "y0=0\n         First row to search for template",
  "y1=30\n        Last row to search for template",
  "region=\n      Miriad style region definition of the template (**not implemented**)",
  "yscale=1\n     Scale the Dy values to more presentable numbers",
  "VERSION=0.2\n  29-may-2013 PJT",
  NULL,
};

string usage="PV correlation search";

string cvsid="$Id$";


#define RPD (PI/180.0)
#ifndef HUGE
#define HUGE 1.0e20
#endif

local void pv_corr(imageptr iptr, real clip, int y0, int y1, real yscale);

local real gfit1(int n, real *s, real *clip, int ipeak, real sig, int mode);

local real qabs(real x, bool Q);


nemo_main()
{
  stream  instr;
  imageptr iptr = NULL;
  real clip, yscale;
  int y0, y1;
  
  instr = stropen(getparam("in"), "r");     /* get file name and open file */
  read_image( instr, &iptr);                /* read image */
  strclose(instr);                          /* close file */
  
  clip = getrparam("clip");
  y0 = getiparam("y0");
  y1 = getiparam("y1");
  yscale = getrparam("yscale");
  pv_corr(iptr, clip, y0, y1, yscale);

}

/*
 * CORR:  in PV diagram, only in V direction, but search for template first
 *
 */

#define MAXV  32
#define MV(i,j)   MapValue(iptr,i,j)

local void pv_corr(imageptr iptr, real clip, int ymin, int ymax, real yscale)
{
  int  ix, iy, nx, ny, dy, delta;
  real mv, sum, imax;
  int *y0, *y1, y0min, y1max;

    
  nx = Nx(iptr);     /* assumed to be position for now */
  ny = Ny(iptr);     /* assumed to be velocity for now */
  y0 = (int *) allocate(nx*sizeof(int));
  y1 = (int *) allocate(nx*sizeof(int));

  imax = MapMax(iptr);
  if (clip>imax) error("mapmax=%g too large for clip=%g\n",imax,clip);

  y0min = ny;
  y1max = -1;
  for (ix=0; ix<nx; ix++) { 
    y0[ix] = y1[ix] = -1;
    for (iy=ymin; iy<ymax; iy++) {
      mv = MapValue(iptr,ix,iy);
      if (mv < clip && y0[ix]<0) continue;
      if (y0[ix]<0 && mv > clip) {
	y0[ix] = iy;
	if (iy < y0min) y0min = iy;
	continue;
      }
      if (y0[ix]>=0 && mv < clip) {
	y1[ix] = iy;
	if (iy > y1max) y1max = iy;
	break;
      }
    }
    dprintf(1,"%d %d %d\n",ix,y0[ix],y1[ix]);
  }
  printf("# y0min=%d y1max=%d\n",y0min,y1max);

  delta = y1max-y0min+1;
  for (dy=-y0min; dy<ny; dy++) {
    sum = 0.0;
    for (ix=0; ix<nx; ix++) {
      if (y0[ix] < 0) continue;
      for (iy=y0[ix]; iy<=y1[ix]; iy++) {
	if(iy+dy<0 || iy+dy>=ny) continue;
	sum += MV(ix,iy) * MV(ix,iy+dy);
      }
    }
    printf("%g %g %d\n",dy*Dy(iptr)*yscale,sum,dy);
  }


  
}

local real qabs(real x, bool Q)
{
  if (Q && x < 0) return -x;
  return x;
}


/*
 * fit GAUSS1d:       y = a + b * exp( - (x-c)^2/(2*d^2) )
 *
 */


static real func_gauss1d(real *x, real *p, int np)
{
  real a,b,arg;
  a = p[2]-x[0];
  b = p[3];
  arg = a*a/(2*b*b);
  return p[0] + p[1] * exp(-arg);
}

static void derv_gauss1d(real *x, real *p, real *e, int np)
{
  real a,b,arg;
  a = p[2]-x[0];
  b = p[3];
  arg = a*a/(2*b*b);
  e[0] = 1.0;
  e[1] = exp(-arg);
  e[2] = -p[1]*e[1] * a   / (b*b);
  e[3] =  p[1]*e[1] * a*a / (b*b*b);
}

#define MAXP  1024

/*
 *  gfit fits a gauss, returns the central velocity
 *  in integer coordinates, which must be 0..n-1
 *  returns -1 if bad fit, no data etc.
 *  mode = 0: full gaussian
 *  mode = 1: half gaussian
 *  mode = 2: half gaussian + 0.5*FWHM (needed for MET)
 *
 */

local real gfit1(int n, real *s, real *clip, int ipeak, real sig, int mode)
{
  real x[MAXP], y[MAXP], fpar[4], epar[4], tol, lab;
  int i, npt, mpar[4], nrt, itmax;

  if (ipeak+1 >= n) return -1.0;

  /* accumulate the points; for half-gaussian, take one point before the peak
   * and then all after, as long as it's not inside the clip range 
   */

  if (mode == 0) {         /* gaussian */
    for (i=0, npt=0; i<n;  i++) {
      if (s[i] <= clip[0] || s[i] >= clip[1]) {
	if (npt==MAXP) error("Too many points, MAXP=%d",MAXP);
	x[npt] = i;
	y[npt] = s[i];
	npt++;
      }
    }
  } else {                /* half gaussian */
    for (i=ipeak-1, npt=0; i<n;  i++) {
      if (s[i] <= clip[0] || s[i] >= clip[1]) {
	if (npt==MAXP) error("Too many points, MAXP=%d",MAXP);
	x[npt] = i;
	y[npt] = s[i];
	npt++;
      }
    }
  } 

  if (npt < 3) return -1.0;    /* not enuf points for fit */


  while (sig > 0.5) {
    /* determine some reasonable initial estimates */
    fpar[0] = 0.0;      /* we'll fix this, there better be no continuum here */
    fpar[1] = s[ipeak]; /* peak value */
    fpar[2] = ipeak;    /* fitting is done in integer channel space */
    fpar[3] = sig;      /* this is the tricky one : we iterate on it */

    mpar[0] = 0;
    mpar[1] = mpar[2] = mpar[3] = 1;

    tol = 0.001;
    lab = 0.001;
    itmax = 50;
    
    nrt = nllsqfit(x,1,y,NULL,NULL,npt,fpar,epar,mpar,4,tol,itmax,lab, 
		 func_gauss1d, derv_gauss1d);
    if (nrt < 0) {
      sig = sig - 1.0;  /* lower sig (units are channel #s) for next iter */
      continue;
    }

    dprintf(1,"@ %d/%d %d %d %g: %g %g    %g %g  %d\n",ipeak,n,npt,mode,sig,
	    fpar[2],epar[2], fpar[3],epar[3],nrt);

    if (nrt > 0) {
      if (mode==2) return fpar[2] + 1.1774*fpar[3];   /* MET ;  1.1774 = sqrt(2ln2) */
      return fpar[2];
    } else {
#if 0
      printf("# %g %g %d %g\n",0.0,s[ipeak],ipeak,sig);
      for(i=0; i<npt; i++)
	printf("%g %g\n",x[i],y[i]);
#endif
      return -1;
    }
  }
  return -1;
}

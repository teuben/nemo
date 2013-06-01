/* 
 * PVCORR:   PV diagram correlation to find multiple weaker lines that
 *           match a template.
 *
 *     29-may-2013  V0.1    for astute, cloned off pvtrace, Q&D
 *     30-may-2013   0.4    extension to cubes
 *     31-may-2013   0.5    mode=0 implemented (1=default), and a simple rounded mode=2
 *      1-jun-2013   0.6    auto-peak segment searching phase 1 implemented for PV, fixed y1[] 
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <mdarray.h>

string defv[] = {
  "in=???\n       Input image file - must be an PV map or XYV cube",
  "clip=???\n     Clip value to search for profile and template area",
  "v0=\n          First row/plane (1 based) to search for template (default uses peak line)",
  "v1=\n          Last row/plane to search for template",
  "vscale=1\n     Scale the velocity values to more presentable numbers",
  "region=\n      Miriad style region definition of the template (**not implemented**)",
  "rms=\n         Use RMS to scale the cross correlations (**experimental**)",
  "mode=1\n       0=template cross-corr no weight  1=intensity weight   2=sample in VMEAN",
  "out=\n         Optional output. For PV this is a table. For XYV this is a 3 plane map",
  "VERSION=0.6\n  1-jun-2013 PJT",
  NULL,
};

string usage="PV map or XYV cube cross correlation search for spectral lines";

string cvsid="$Id$";


local void pv_corr(imageptr iptr, real clip, int y0, int y1, real yscale, int mode, stream outstr);
local void xyv_corr(imageptr iptr, real clip, int z0, int z1, real zscale, int mode, stream outstr);
local real gfit1(int n, real *s, real *clip, int ipeak, real sig, int mode);

local real rms = -1.0;

local real psqrt(real x) { return x<0 ? -sqrt(-x) : sqrt(x); }

nemo_main()
{
  stream  instr, outstr = NULL;
  imageptr iptr = NULL;
  real clip, vscale;
  int v0, v1, mode;
  
  instr = stropen(getparam("in"), "r");     /* get file name and open file */
  read_image( instr, &iptr);                /* read image */
  strclose(instr);                          /* close file */
  if (hasvalue("out")) outstr = stropen(getparam("out"),"w");
  
  clip = getrparam("clip");
  if (hasvalue("v0") && hasvalue("v1")) {
    v0 = getiparam("v0") - 1;   /* user supplies 1-based numbers */
    v1 = getiparam("v1") - 1;   /* internally we use 0-based V's */
    if (v0<0 || v1<0 || v1<=v0) error("bad values v0,v1=%d,%d",v0+1,v1+1);
  } else
    v0 = v1 = -1;
  vscale = getrparam("vscale");
  mode = getiparam("mode");
  if (mode<0 || mode>2) error("invalid mode=%d",mode);
  if (hasvalue("rms")) rms = getrparam("rms");
  if (Nz(iptr) > 1)
    xyv_corr(iptr, clip, v0, v1, vscale, mode, outstr);
  else
    pv_corr(iptr, clip, v0, v1, vscale, mode, outstr);

}

/* 
 * helper macros and functions for PV diagram 
 */

#define MV(i,j)   MapValue(iptr,i,j)

/* find the peak , if not found, allow to grow a bit */

int pv_peak(imageptr iptr, int ix, int iy0, int iy1, real clip, int grow)
{
  int iy, iymax = iy0;
  real peak = MV(ix,iy0);

  if (iy0<0 || iy1<0) return -1;   /* need to find solution for this */
  
  for (iy=iy0; iy<=iy1; iy++) {
    if (MV(ix,iy) > peak) {
      iymax = iy;
      peak = MV(ix,iy);
    }
  }
  if (peak > clip) return iymax;
   
  if (grow) warning("pv_peak: growing not implemented yet");

  /* if not found by now, return failure -1 */
  return -1;
}

/* find lower clip valued boundary, v0 */

int pv_y0(imageptr iptr, int ix, int iymax, real clip)
{
  int iy;

  for (iy=iymax; iy>=0; iy--) {
    if (MV(ix,iy) < clip) return iy+1;
  }
  warning("edge-0 > clip");
  return 0;  /* bad, this means edge still > clip */
}



/* find upper clip valued boundary, v1 */

int pv_y1(imageptr iptr, int ix, int iymax, real clip)
{
  int iy;

  for (iy=iymax; iy<Ny(iptr); iy++) {
    if (MV(ix,iy) < clip) return iy-1;
  }
  warning("edge-1 > clip");
  return Ny(iptr)-1;   /* bad, this means edge still > clip */

}


/*
 * CORR:  in PV diagram, only in V direction, but search for template first
 *
 */


local void pv_corr(imageptr iptr, real clip, int ymin, int ymax, real yscale, int mode, stream outstr)
{
  int  ix, iy, nx, ny, np, dy, delta, imax_x, imax_y;
  real mv, sum, sum0, imax, sumi, sumiv, *ym;
  int *y0, *y1, y0min, y1max, iylm;

    
  nx = Nx(iptr);     /* assumed to be position for now */
  ny = Ny(iptr);     /* assumed to be velocity for now */
  y0 = (int *) allocate(nx*sizeof(int));
  y1 = (int *) allocate(nx*sizeof(int));
  ym = (real *) allocate(nx*sizeof(real));

  imax = MapMax(iptr);
  if (clip>imax) error("mapmax=%g too large for clip=%g\n",imax,clip);
  if (ymax >=ny) error("invalid v1=%d",ymax);

  y0min = ny;
  y1max = -1;
  np = 0;
  if (ymin>=0 && ymax>=0) {       /* if (v0,v1) were given, force looking between those */
    for (ix=0; ix<nx; ix++) { 
      if (MV(ix,ymin) > clip) error("Clip too large, PV(%d,%d)=%g",ix,ymin,MV(ix,ymin));
      if (MV(ix,ymax) > clip) error("Clip too large, PV(%d,%d)=%g",ix,ymax,MV(ix,ymax));
      y0[ix] = y1[ix] = -1;
      for (iy=ymin; iy<=ymax; iy++) {
	mv = MV(ix,iy);
	if (mv < clip && y0[ix]<0) continue;
	if (y0[ix]<0 && mv > clip) {
	  y0[ix] = iy;
	  if (y0[ix] < y0min) y0min = y0[ix];
	  np++;
	  continue;
	}
	if (y0[ix]>=0 && mv < clip) {
	  y1[ix] = iy-1;
	  if (y1[ix] > y1max) y1max = y1[ix];
	  break;
	}
      }
    }
  } else {                          /* else find peak, and look around the peak */
    imax = MV(0,0);
    for (ix=0; ix<nx; ix++) {       /* loop over PV and find the location of peak */
      y0[ix] = y1[ix] = -1;
      for (iy=0; iy<ny; iy++) {    
	if (MV(ix,iy)>imax) {
	  imax   = MV(ix,iy);
	  imax_x = ix;
	  imax_y = iy;
	}
      }
    }
    dprintf(0,"Peak in map: %g @ (%d,%d)\n",imax,imax_x,imax_y);
    if (imax<clip) error("clip and imax bad");

    y0[imax_x] = pv_y0(iptr,imax_x,imax_y,clip);
    y1[imax_x] = pv_y1(iptr,imax_x,imax_y,clip);
    dprintf(0,"MaxPeak clip segment: %d %d %d %d\n",imax_x,imax_y,y0[imax_x],y1[imax_x]); 
    delta =  y1[imax_x]  - y0[imax_x] + 1;          
    for (ix=imax_x+1; ix<nx; ix++) {
      iylm = pv_peak(iptr, ix, y0[ix-1], y1[ix-1], clip, 0);
      if (iylm < 0) break;  /* for now - need solition to peak up new blobs */
      y0[ix] = pv_y0(iptr,ix,iylm,clip);
      y1[ix] = pv_y1(iptr,ix,iylm,clip);
      np += y1[ix]-y0[ix]+1;
      if (y0[ix] < y0min) y0min = y0[ix];
      if (y1[ix] > y1max) y1max = y1[ix];
      dprintf(0,"Up clip segment: %d %d %d %d\n",ix,iylm,y0[ix],y1[ix]); 
    } 
    for (ix=imax_x-1; ix>=0; ix--) {
      iylm = pv_peak(iptr, ix, y0[ix+1]-delta, y1[ix+1]+delta, clip, 0);
      if (iylm < 0) break;  /* for now - need solition to peak up new blobs */
      y0[ix] = pv_y0(iptr,ix,iylm,clip);
      y1[ix] = pv_y1(iptr,ix,iylm,clip);
      np += y1[ix]-y0[ix]+1;
      if (y0[ix] < y0min) y0min = y0[ix];
      if (y1[ix] > y1max) y1max = y1[ix];
      dprintf(0,"Down clip segment: %d %d %d %d\n",ix,iylm,y0[ix],y1[ix]); 
    }
  }
  printf("# y0min=%d y1max=%d np=%d\n",y0min,y1max,np);

  /*  
   *  compute some potentially normalizing sum , and
   *  grab the mean velocity in case we need that as
   *  well
   *  @todo   sum0 for mode=2 needs adjusting
   */
  sum = 0.0;
  for (ix=0; ix<nx; ix++) {
    ym[ix] = 0;
    sumi = sumiv = 0.0;
    if (y0[ix] < 0) continue;
    for (iy=y0[ix]; iy<=y1[ix]; iy++) {
      sumi  = sumi  + MV(ix,iy);
      sumiv = sumiv + MV(ix,iy) * iy;
      if (rms > 0)
	sum0 += rms*rms;
      else
	sum0 += MV(ix,iy)*MV(ix,iy);
    }
    ym[ix] = sumiv/sumi;
    if (outstr) fprintf(outstr,"%d %d %d %g\n",ix,y0[ix],y1[ix],ym[ix]);
  }
  printf("# sum0=%g\n",sum0);

  delta = y1max-y0min+1;
  for (dy=-y1max; dy<ny-y0min; dy++) {
    sum = 0.0;
    for (ix=0; ix<nx; ix++) {
      if (y0[ix] < 0) continue;
      if (mode==0 || mode==1) {
	for (iy=y0[ix]; iy<=y1[ix]; iy++) {
	  if(iy+dy<0 || iy+dy>=ny) continue;
	  if (mode == 0)
	    sum += MV(ix,iy+dy);
	  else /* mode == 1 */
	    sum += MV(ix,iy) * MV(ix,iy+dy);
	}
      } else { /* mode == 2 */
	iy = dy + ym[ix];
	if (iy<0 || iy>=ny) continue;
	sum += MV(ix,iy);  /* @todo need to properly interpolate */
      }
    }
    printf("%g %g %d\n",dy*Dy(iptr)*yscale,(sum/sum0),dy);
  }
}



#define CV(i,j,k)   CubeValue(iptr,i,j,k)

/* find the peak , if not found, allow to grow a bit */

int xyv_peak(imageptr iptr, int ix, int iy, int iv0, int iv1, real clip, int grow)
{
  /* if not found by now, return failure -1 */
  return -1;
}

/* find lower clip valued boundary, v0 */

int xyv_v0(imageptr iptr, int ix, int ivmax, real clip)
{
  return 0;  /* bad, this means edge still > clip */
}



/* find upper clip valued boundary, v1 */

int xyv_v1(imageptr iptr, int ix, int ivmax, real clip)
{
  return Nz(iptr)-1;   /* bad, this means edge still > clip */
}


local void xyv_corr(imageptr iptr, real clip, int vmin, int vmax, real vscale, int mode, stream outstr)
{
  int  ix, iy, iz, nx, ny, nz, nxy, dv, delta, imax_x, imax_y;
  real cv, sum, imax;
  int v0min, v1max;
  mdarray2 v0, v1;

  nx = Nx(iptr);     /* assumed to be position for now */
  ny = Ny(iptr);     /* assumed to be position for now */
  nz = Nz(iptr);     /* assumed to be velocity for now */
  v0 = allocate_mdarray2(ny,nx);
  v1 = allocate_mdarray2(ny,nx);

  imax = MapMax(iptr);
  if (clip>imax) error("mapmax=%g too large for clip=%g\n",imax,clip);
  if (vmax >=nz) error("invalid v1=%d",vmax);

  nxy = 0;
  v0min = nz;
  v1max = -1;
  if (vmin>=0 && vmax>=0) {
    for (ix=0; ix<nx; ix++) { 
      for (iy=0; iy<ny; iy++) { 
	v0[iy][ix] = v1[iy][ix] = -1;
	for (iz=vmin; iz<=vmax; iz++) {
	  cv = CV(ix,iy,iz);
	  if (cv < clip && v0[iy][ix]<0) continue;
	  if (v0[iy][ix]<0 && cv > clip) {
	    v0[iy][ix] = iz;
	    if (iz < v0min) v0min = iz;
	    nxy++;
	    continue;
	  }
	  if (v0[iy][ix]>=0 && cv < clip) {
	    v1[iy][ix] = iz;
	    if (iz > v1max) v1max = iz;
	    break;
	  }
	}
	dprintf(1,"%d %d %g %g\n",ix,iy,v0[iy][ix],v1[iy][ix]);
      }
    }
    printf("# v0min=%d v1max=%d nxy=%d\n",v0min,v1max,nxy);
  } else {
    error("auto peak and template finding code not finished here");
  }

  delta = v1max-v0min+1;
  for (dv=-v1max; dv<nz-v0min; dv++) {
    sum = 0.0;
    for (ix=0; ix<nx; ix++) {
      for (iy=0; iy<ny; iy++) {
	if (v0[iy][ix] < 0) continue;
	for (iz=v0[iy][ix]; iz<=v1[iy][ix]; iz++) {
	  if(iz+dv<0 || iz+dv>=nz) continue;
	  sum += CV(ix,iy,iz) * CV(ix,iy,iz+dv);
	}
      }
    }
    printf("%g %g %d\n",dv*Dz(iptr)*vscale,sum,dv);
  }
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

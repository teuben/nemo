/* 
 * CCDGEN:   create 2D 'astronomical' type objects, modeled after MIRIAD'd imgen
 *
 *      4-jan-05        V0.1: for Kartik's bars modeling experiment 
 *                      (IMGEN is just tooooo cumbersome, and so is ccdmath) 
 *      6-jan-05        V0.7: added many more models and features. added factor=
 *      8-jan-05        V0.8: add m!=2 multi-arm spirals
 *      5-aug-11        V0.9: object = test
 *     15-oct-2014      V0.10:  object=noise now can do 3D
 *     12-mar-2020      V1.0    object=blobs in 3D
 * TODO:
 *    - find out why the normalization was PI, and not TWO_PI, which worked before.
 *       (this happened when I changed from 1 to 1/3600 scaling factor in the examples)
 *    - add a few more of the imgen models
 *    - replace crpix/crval/cdelt with cell= and maybe radec= ? as in imgen
 *    - if A=0, then factor most like won't work (hard to fix the way it's written)
 *
 */

#include <stdinc.h>
#include <ctype.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <strlib.h>
#include <image.h>

string defv[] = {
  "in=\n           Input file (optional) to be added to the new file",
  "out=???\n       Output file",
  "object=flat\n   Object type (test,flat,exp,gauss,bar,spiral,....)",
  "spar=\n         Parameters for this object",
  "center=\n       Center of object (defaults to map(reference)center) in pixels 0..size-1",
  "size=10,10,1\n  2- or 3D dimensions of map/cube (only if no input file given)",
  "cell=1\n        Cellsize",
  "pa=0\n          Position Angle of disk in which object lives",
  "inc=0\n         Inclination Angle of disk in which object lives",
  "totflux=f\n     Interpret peak values (spar[0]) as total flux instead",
  "factor=1\n      Factor to multiple input image by before adding object",
  "crpix=\n        Override/Set crpix (1,1,1) // ignored ** mapcenter if left blank **",
  "crval=\n        Override/Set crval (0,0,0) // ignored",
  "cdelt=\n        Override/Set cdelt (1,1,1) // ignored",
  "seed=0\n        Random seed",
  "headline=\n     Random veriage for the history",
  "VERSION=1.0\n   11-mar-2020 PJT",
  NULL,
};

string usage = "image creation/modification with objects";

string cvsid = "$Id$";

#ifndef HUGE
# define HUGE 1.0e20
#endif

#define MAXIMAGE 20
#define MAXPAR  32768           /* max in nemoinp */

imageptr iptr;	                /* pointers to (input) image */

#define MAXNAX 3

double crval[MAXNAX], crpix[MAXNAX], cdelt[MAXNAX];
int nwcs = 0;

real spar[MAXPAR];
int npar;

real pa,inc,center[2];
real sinp, cosp, sini, cosi;

bool Qtotflux;
real factor;
real surface=1.0;


local void do_create(int nx, int ny, int nz);

local void object_test(int npars, real *pars);
local void object_flat(int npars, real *pars);
local void object_exp(int npars, real *pars);
local void object_gauss(int npars, real *pars);
local void object_bar(int npars, real *pars);
local void object_ferrers(int npars, real *pars);
local void object_spiral(int npars, real *pars);
local void object_noise(int npars, real *pars);
local void object_j1x(int npars, real *pars);
local void object_isothermal(int npars, real *pars);
local void object_comet(int npars, real *pars);
local void object_jet(int npars, real *pars);
local void object_shell(int npars, real *pars);
local void object_point(int npars, real *pars);
local void object_blobs(int ndim, int npars, real *pars);

extern string *burststring(string,string);


void nemo_main ()
{
  string  fnames;
  stream  instr;                      /* input file (optional) */
  stream  outstr;                     /* output file */
  int     size[3],nx, ny, nz;         /* size of scratch map */
  int     ncen;
  string  object;
  string  headline;
  int seed = init_xrandom(getparam("seed"));
  int ndim;

  object = getparam("object");
  npar = nemoinpr(getparam("spar"),spar,MAXPAR);
  if (npar < 0) error("Syntax error %s",getparam("spar"));

  Qtotflux = getbparam("totflux");
  factor = getdparam("factor");

  pa = getdparam("pa");
  inc = getdparam("inc");
  sinp = sin(pa*PI/180.0);
  cosp = cos(pa*PI/180.0);
  sini = sin(inc*PI/180.0);
  cosi = cos(inc*PI/180.0);
  dprintf(0,"%s: disk with pa=%g inc=%g\n",object,pa,inc);
  dprintf(1,"pa(%g %g) inc(%g %g)\n",sinp,cosp,sini,cosi);

  ncen = nemoinpr(getparam("center"),center,2);
  if (ncen<0) error("Syntax error %s",getparam("center"));
  if (ncen==1) center[1] = center[0];

  nwcs = nemorinpd(getparam("crval"),crval,MAXNAX,0.0,FALSE);
  nwcs = nemorinpd(getparam("cdelt"),cdelt,MAXNAX,1.0,FALSE);
  nwcs = nemorinpd(getparam("crpix"),crpix,MAXNAX,1.0,FALSE);
  if (hasvalue("in")) {
    instr = stropen(getparam("in"),"r");    /* open file */
    read_image (instr, &iptr);
  } else {
    dprintf(0,"Generating a map from scratch\n");
    switch (nemoinpi(getparam("size"),size,3)) {
    case 1:			/*  nx[,nx,1] */
      size[1] = size[0];
      size[2] = 1;
      ndim = 2;
      break;
    case 2:			/*  nx,ny[,1] */
      size[2] = 1;
      ndim = 2;
      break;
    case 3:			/*  nx,ny,nz  */
      ndim = 3;
      break;
    case 0:			/*  [10,10,1] */
      dprintf(0,"Cannot have no size, default 10 assumed\n");
      size[0] = size[1] = 10;
      size[2] = 1;
      ndim = 2;
      break;
    default:			/* --- some error --- */
      error("Syntax error in size keyword\n");
    }
    nx = size[0];
    ny = size[1];
    nz = size[2];
    do_create(nx,ny,nz);
#if 0
    if (nwcs == 0) {    /* last one for crpix */
      warning("going to set the center of the map in new NEMO convention");
      crpix[0] = (nx-1)/2.0;
      crpix[1] = (ny-1)/2.0;
      crpix[2] = (nz-1)/2.0;
    }
#endif

  }
  if (ncen==0) {   /* fix center if it has not been set yet */
    center[0] = (Nx(iptr)-1)/2.0;     /* 0 based center= */
    center[1] = (Ny(iptr)-1)/2.0;
  }
  dprintf(0,"%s: center pixel: %g %g\n",object,center[0],center[1]);
  surface = Dx(iptr)*Dy(iptr);
  surface = ABS(surface);

  if (streq(object,"flat"))
    object_flat(npar,spar);
  else if (streq(object,"exp"))
    object_exp(npar,spar);
  else if (streq(object,"gauss"))
    object_gauss(npar,spar);
  else if (streq(object,"bar"))
    object_bar(npar,spar);
  else if (streq(object,"ferrers"))
    object_ferrers(npar,spar);
  else if (streq(object,"spiral"))
    object_spiral(npar,spar);
  else if (streq(object,"noise"))
    object_noise(npar,spar);
  else if (streq(object,"jet"))
    object_jet(npar,spar);
  else if (streq(object,"j1x"))
    object_j1x(npar,spar);
  else if (streq(object,"isothermal"))
    object_isothermal(npar,spar);
  else if (streq(object,"comet"))
    object_comet(npar,spar);
  else if (streq(object,"shell"))
    object_shell(npar,spar);
  else if (streq(object,"point"))
    object_point(npar,spar);
  else if (streq(object,"test"))
    object_test(npar,spar);
  else if (streq(object,"blobs"))
    object_blobs(ndim,npar,spar);
  else
    error("Unknown object %g",object);
  
  outstr = stropen (getparam("out"),"w");  /* open output file first ... */
  if (hasvalue("headline"))
    set_headline(getparam("headline"));
  write_image (outstr,iptr);         /* write image to file */
  strclose(outstr);
}

/*
 *  create new map from scratch, using %x and %y as position parameters 
 *		0..nx-1 and 0..ny-1
 */
local void do_create(int nx, int ny, int nz)
{
    double m_min, m_max, total;
    real   fin[5], fout;
    int    ix, iy, iz;
    int    badvalues;
    
    m_min = HUGE; m_max = -HUGE;
    total = 0.0;		/* count total intensity in new map */
    badvalues = 0;		/* count number of bad operations */

    if (nz > 0) {
      // warning("cube");
      if (!create_cube (&iptr, nx, ny, nz))	/* create default empty image */
        error("Could not create 3D image from scratch");
#if 0      
      Axis(iptr) = 1;      /* set linear axistype with a fits-style reference pixel */
#endif
      wcs_f2i(3,crpix,crval,cdelt,iptr);

      for (iz=0; iz<nz; iz++) {
        fin[2] = iz-crpix[2]+1;   /* crpix is 1 for first pixel (FITS convention) */
        for (iy=0; iy<ny; iy++) {
	  fin[1] = iy-crpix[1]+1;
	  for (ix=0; ix<nx; ix++) {
	    fout = 0.0;
	    CubeValue(iptr,ix,iy,iz) = fout;
	    m_min = MIN(m_min,fout);         /* and check for new minmax */
	    m_max = MAX(m_max,fout);
	    total += fout;                   /* add up totals */
	  }
        }
      }
    } else {
      warning("2d-map");
      if (!create_image (&iptr, nx, ny))	
        error("Could not create 2D image from scratch");
      wcs_f2i(2,crpix,crval,cdelt,iptr);

      for (iy=0; iy<ny; iy++) {
	fin[1] = iy;
	for (ix=0; ix<nx; ix++) {
	  fout = 0.0;
	  MapValue(iptr,ix,iy) = fout;
	  m_min = MIN(m_min,fout);         /* and check for new minmax */
	  m_max = MAX(m_max,fout);
	  total += fout;                   /* add up totals */
	}
      }
    } 
    
    MapMin(iptr) = m_min;
    MapMax(iptr) = m_max;

    dprintf(1,"New min and max in map are: %f %f\n",m_min,m_max);
    dprintf(1,"New total brightness/mass is %f\n",
			total*Dx(iptr)*Dy(iptr));
    if (badvalues)
    	warning ("There were %d bad operations in dofie",badvalues);
}


local void object_test(int npars, real *pars)
{
  int i,j,k;
  int nx = Nx(iptr);
  int ny = Ny(iptr);
  int nz = Ny(iptr);
  real A = 1.0;

  for (k=0; k<nz; k++)
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++)
	CubeValue(iptr,i,j,k) = i + 10*j + 100*k;
}


local void object_flat(int npars, real *pars)
{
  int i,j;
  int nx = Nx(iptr);
  int ny = Ny(iptr);
  real A = 1.0;

  if (A==0) return;

  if (npar > 0) A = pars[0];

  if (Qtotflux) A /= (nx*ny);

  for (j=0; j<ny; j++)
    for (i=0; i<nx; i++)
      MapValue(iptr,i,j) = factor*MapValue(iptr,i,j) + A;

}

local void object_exp(int npars, real *pars)
{
  int i,nx = Nx(iptr);
  int j,ny = Ny(iptr);
  real A = 1.0;
  real h = 1.0;
  real x1,y1,x2,y2,r,arg,value;

  if (npar > 0) A = pars[0];
  if (npar > 1) h = pars[1];
  dprintf(0,"exp: %g %g\n",A,h);

  if (A==0) return;

  if (Qtotflux) {
    A /= (PI*h*h/surface);
    dprintf(0,"exp: A->%g\n",A);
  }

  for (j=0; j<ny; j++) {
    y1 = (j-center[1])*Dy(iptr) + Ymin(iptr);
    for (i=0; i<nx; i++) {
      x1 = (i-center[0])*Dx(iptr) + Xmin(iptr);
      x2 =  -x1*sinp - y1*cosp;
      y2 =   x1*cosp - y1*sinp;
      r = sqrt(x2*x2 + sqr(y2/cosi));
      arg = r/h;
      dprintf(2,"%d %d : %g %g %g %g  %g\n",i,j,x1,y1,x2,y2,arg);
      value = (arg < 80) ? A * exp(-arg) : 0.0;
      MapValue(iptr,i,j) = factor*MapValue(iptr,i,j) + value;
    }
  }
}

local void object_gauss(int npars, real *pars)
{
  int i,nx = Nx(iptr);
  int j,ny = Ny(iptr);
  real A = 1.0;
  real h = 1.0;
  real x1,y1,x2,y2,r,arg,value;

  if (npar > 0) A = pars[0];
  if (npar > 1) h = pars[1];
  dprintf(0,"gauss: %g %g\n",A,h);

  if (A==0) return;

  if (Qtotflux) {
    A /= (PI*h*h/surface);
    dprintf(0,"gauss: A->%g\n",A);
  }

  for (j=0; j<ny; j++) {
    y1 = (j-center[1])*Dy(iptr) + Ymin(iptr);
    for (i=0; i<nx; i++) {
      x1 = (i-center[0])*Dx(iptr) + Xmin(iptr);
      x2 =  -x1*sinp - y1*cosp;
      y2 =   x1*cosp - y1*sinp;
      r = sqrt(x2*x2 + sqr(y2/cosi));
      arg = r/h;
      dprintf(2,"%d %d : %g %g %g %g   %g\n",i,j,x1,y1,x2,y2,arg);
      value = (arg < 13) ?  A * exp(-0.5*arg*arg) : 0.0;
      MapValue(iptr,i,j)  = factor*MapValue(iptr,i,j) + value;
    }
  }
}

local void object_bar(int npars, real *pars)
{
  int i,nx = Nx(iptr);
  int j,ny = Ny(iptr);
  real A = 1.0;
  real h = 1.0;
  real e = 0.0;
  real b = 0.0;
  real x1,y1,x2,y2,x3,y3,r,arg,value;
  real sinb,cosb,omba;

  if (npar > 0) A = pars[0];
  if (npar > 1) h = pars[1];
  if (npar > 2) e = pars[2];
  if (npar > 3) b = pars[3];
  dprintf(0,"bar: %g %g %g %g\n",A,h,e,b);

  sinb = sin(b*PI/180.0);
  cosb = cos(b*PI/180.0);

  if (A==0) return;

  if (Qtotflux) {
    A /= (PI*h*h*(1-e)/surface);
    dprintf(0,"bar: A->%g\n",A);
  }
  dprintf(1,"bar b=%g\n",b);

  for (j=0; j<ny; j++) {
    y1 = (j-center[1])*Dy(iptr) + Ymin(iptr);
    for (i=0; i<nx; i++) {
      x1 = (i-center[0])*Dx(iptr) + Xmin(iptr);

      x2 =  -x1*sinp - y1*cosp;
      y2 =   (x1*cosp - y1*sinp)/cosi;

      x3 =   x2*cosb - y2*sinb;
      y3 = (x2*sinb  + y2*cosb)/(1-e);

      r = sqrt(x3*x3+y3*y3);
      arg = r/h;
      value = (arg < 100) ? A * exp(-arg) : 0.0;
      MapValue(iptr,i,j) = factor*MapValue(iptr,i,j) + value;
    }
  }
}


local void object_ferrers(int npars, real *pars)
{
  int i,nx = Nx(iptr);
  int j,ny = Ny(iptr);
  int l, lmax;
  real A = 1.0;   /* peak */
  real h = 1.0;   /* size */
  real e = 0.0;   /* 1-b/a */
  real b = 0.0;   /* phi */
  real p = 1.0;   /* power factor */
  real x1,y1,x2,y2,x3,y3,r,arg,value;
  real amp,phi,sum;
  real sinb,cosb;

  if (npar > 0) A = pars[0];
  if (npar > 1) h = pars[1];
  if (npar > 2) e = pars[2];
  if (npar > 3) b = pars[3];
  if (npar > 4) p = pars[4];
  dprintf(0,"ferrers: %g %g %g %g %g\n",A,h,e,b,p);

  if (A==0) return;

  lmax = Qtotflux ? 2 : 1;

  sinb = sin(b*PI/180.0);
  cosb = cos(b*PI/180.0);

  for (l=0; l<lmax; l++) {     /* 1st loop: sum up the flux   if in 2nd loop: normalize */
    if (l==0) 
      sum = 0.0;
    else {
      A /= sum;
      dprintf(0,"ferrers: A->%g\n",A);
    } 
    for (j=0; j<ny; j++) {
      y1 = (j-center[1])*Dy(iptr) + Ymin(iptr);
      for (i=0; i<nx; i++) {
	x1 = (i-center[0])*Dx(iptr) + Xmin(iptr);
	
	x2 =  -x1*sinp - y1*cosp;
	y2 =  (x1*cosp - y1*sinp)/cosi;

	x3 =   x2*cosb - y2*sinb;
	y3 = (x2*sinb  + y2*cosb)/(1-e);

	r = (x3*x3+y3*y3)/(h*h);
	value = r < 1 ? A * pow(1.0-r, p) : 0.0;

	if (Qtotflux) {   
	  if (l==0) 
	    sum += value;
	  else
	    MapValue(iptr,i,j) = factor*MapValue(iptr,i,j) + value;
	} else
	  MapValue(iptr,i,j) = factor*MapValue(iptr,i,j) + value;
      } /* i */
    } /* j */
  } /* k */
}


local void object_spiral(int npars, real *pars)
{
  int i,nx = Nx(iptr);
  int j,ny = Ny(iptr);
  int l, lmax;
  real A = 1.0;   /* Amplitude */
  real h = 1.0;   /* radial scale length */
  real k = 1.0;   /* wave number */
  real p = 1.0;   /* 1/2 power of cos */
  real r0 = 0.0;  /* starting radius */
  real p0 = 0.0;  /* starting angle */
  int m = 2;      /* number of arms */
  real x1,y1,x2,y2,x3,y3,r,arg,value;
  real amp,phi,sum;
  bool Qint, pint;

  if (npar > 0) A = pars[0];
  if (npar > 1) h = pars[1];
  if (npar > 2) k = pars[2];
  if (npar > 3) p = pars[3];
  if (npar > 4) m =  (int) pars[4];  /* note rounding ? */
  if (npar > 5) r0 = pars[5];
  if (npar > 6) p0 = pars[6];

  dprintf(0,"spiral: %g %g   %g %g %d   %g %g\n",A,h,k,p,m,r0,p0);

  if (A==0) return;

  p0 *= PI/180.0;  /* convert from degrees to radians */
  k *= TWO_PI;     /* angles are 2.PI.k.r , so absorb 2.PI.k in one */
  pint = (int) p;
  Qint = (p-pint == 0);
  if (Qint)
    warning("Integer power p = %d\n",pint);

  lmax = Qtotflux ? 2 : 1;

  for (l=0; l<lmax; l++) {     /* 1st loop: sum up the flux   if in 2nd loop: normalize */
    if (l==0) 
      sum = 0.0;
    else {
      A /= sum;
      dprintf(0,"spiral: A->%g\n",A);
    } 
    for (j=0; j<ny; j++) {
      y1 = (j-center[1])*Dy(iptr) + Ymin(iptr);
      for (i=0; i<nx; i++) {
	x1 = (i-center[0])*Dx(iptr) + Xmin(iptr);
	
	x2 =  -x1*sinp - y1*cosp;
	y2 =  (x1*cosp - y1*sinp)/cosi;
	r = sqrt(x2*x2+y2*y2);
	if (r < r0) 
	  value = 0.0;
	else {
	  /* ? should match this up better so we can connect bar and spiral ? */
	  phi = atan2(y2,x2) + k*(r-r0) + p0;
	  if (Qint) {
	    amp = powi(cos(m*phi),pint);   /* these can come out negative for odd p's !! */
	  } else
	    amp = powd(cos(m*phi),(double)p);
	  arg = r/h;
	  value = (arg < 80) ? A * amp * exp(-arg) :  0.0;
	}
	if (Qtotflux) {   
	  if (l==0) 
	    sum += value;
	  else
	    MapValue(iptr,i,j) = factor*MapValue(iptr,i,j) + value;
	} else
	  MapValue(iptr,i,j) = factor*MapValue(iptr,i,j) + value;
      } /* i */
    } /* j */
  } /* k */
}

local void object_noise(int npars, real *pars)
{
  int i,j,k;
  int nx = Nx(iptr);
  int ny = Ny(iptr);
  int nz = Nz(iptr);
  double m = 1.0;
  double s = 1.0;

  if (npar > 0) m = pars[0];
  if (npar > 1) s = pars[1];

  dprintf(0,"noise:%g %g\n",m,s);

  if (s==0) return;

  if (Qtotflux) {
    m /= (nx*ny);
    s /= (nx*ny);
    dprintf(0,"noise: m->%g  s->%g\n",m,s);
  }
  
  for (k=0; k<nz; k++)
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++)
	CubeValue(iptr,i,j,k) = factor*CubeValue(iptr,i,j,k) + grandom(m,s);
}



local void object_j1x(int npars, real *pars)
{
  warning("alas, j1x not supported yet");
}

local void object_isothermal(int npars, real *pars)
{
  int i,nx = Nx(iptr);
  int j,ny = Ny(iptr);
  int l, lmax;
  real A = 1.0;   /* peak */
  real h = 1.0;   /* size */
  real p = -0.5;  /* power factor */
  real x1,y1,x2,y2,x3,y3,r,arg,value;
  real amp,phi,sum;
  real sinb,cosb;

  if (npar > 0) A = pars[0];
  if (npar > 1) h = pars[1];
  if (npar > 2) p = pars[2];
  dprintf(0,"isothermal: %g %g %g \n",A,h,p);

  if (A==0) return;

  lmax = Qtotflux ? 2 : 1;

  for (l=0; l<lmax; l++) {     /* 1st loop: sum up the flux   if in 2nd loop: normalize */
    if (l==0) 
      sum = 0.0;
    else {
      A /= sum;
      dprintf(0,"isothermal: A->%g\n",A);
    } 
    for (j=0; j<ny; j++) {
      y1 = (j-center[1])*Dy(iptr) + Ymin(iptr);
      for (i=0; i<nx; i++) {
	x1 = (i-center[0])*Dx(iptr) + Xmin(iptr);
	
	x2 =  -x1*sinp - y1*cosp;
	y2 =  (x1*cosp - y1*sinp)/cosi;

	r = (x2*x2+y2*y2)/(h*h);
	value = A * pow(1.0+r, p);

	if (Qtotflux) {   
	  if (l==0) 
	    sum += value;
	  else
	    MapValue(iptr,i,j) = factor*MapValue(iptr,i,j) + value;
	} else
	  MapValue(iptr,i,j) = factor*MapValue(iptr,i,j) + value;
      } /* i */
    } /* j */
  } /* k */
}

local void object_comet(int npars, real *pars)
{
  warning("alas, comet not supported yet");
}

local void object_jet(int npars, real *pars)
{
  warning("alas, jet not supported yet");
}

local void object_shell(int npars, real *pars)
{
  warning("alas, shell not supported yet");
}

local void object_point(int npars, real *pars)
{
  int nx = Nx(iptr);
  int ny = Ny(iptr);
  int i = nx/2;
  int j = ny/2;

  if (npars == 0) error("object=point requires a value");

  MapValue(iptr,i,j) = factor*MapValue(iptr,i,j) + pars[0];
}


local void object_blobs(int ndim, int npars, real *pars)
{
  int x0, ix, nx = Nx(iptr);
  int y0, iy, ny = Ny(iptr);
  int z0, iz, nz = Ny(iptr);
  int i, j, k;
  int npar = 5;  // peak, x, y, z, s
  int ns, ib, nblobs = npars/npar;
  real pk, s0, dx, dy, dz, arg;

  if (npars == 0) error("object=blobs requires a value");
  if (npars % npar != 0) error("object=blobs not commensurate to npar=%d",npar);
  if (ndim != 3) error("object=blobs not yet coded for 2d");

  dprintf(1,"Found %d blobs\n",nblobs);

  for (ib=0; ib<nblobs; ib++) {
    pk =      pars[ib*npar + 0];
    x0 = (int)pars[ib*npar + 1];
    y0 = (int)pars[ib*npar + 2];
    z0 = (int)pars[ib*npar + 3];
    s0 =      pars[ib*npar + 4];
    ns = (int) (5*s0);
    dprintf(1,"%d  %g  %d %d %d %g\n",ib,pk,x0,y0,z0,s0);

    for (ix=x0-ns; ix<=x0+ns; ix++) {
      if (ix<0 || ix>=nx) continue;
      dx = ix-x0;
      for (iy=y0-ns; iy<=y0+ns; iy++) {
	if (iy<0 || iy>=ny) continue;
	dy = iy-y0;
	for (iz=z0-ns; iz<=z0+ns; iz++) {
	  if (iz<0 || iz>=nz) continue;	  
	  dz = iz-z0;
	  arg = (dx*dx + dy*dy + dz*dz)/(s0*s0);
	  CubeValue(iptr,ix,iy,iz) += pk * exp(-arg);
	}//iz
      }//iy
    }//iz
  }//ib
}

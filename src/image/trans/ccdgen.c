/* 
 * CCDGEN:   create 2D 'astronomical' type objects, modeled after MIRIAD'd imgen
 *
 *      4-jan-05        V0.1: for Kartik's bars modeling experiment (IMGEN is just tooooo cumbersome)
 *                      
 */

#include <stdinc.h>
#include <ctype.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <strlib.h>
#include <image.h>

string defv[] = {
  "in=\n           Input file (optional) to be added to the new file",
  "out=???\n       Output file",
  "object=flat\n   Object type (flat,exp,gauss,bar,spiral,....)",
  "spar=\n         Parameters for this object",
  "center=\n       Center of object (defaults to mapcenter) in pixels 0..size-1",
  "size=10,10,1\n  2- or 3D dimensions of map/cube (only if no input file given)",
  "cell=1\n        Cellsize",
  "pa=0\n          Position Angle of disk in which object lives",
  "inc=0\n         Inclination Angle of disk in which object lives",
  "totflux=f\n     Interpret peak values (spar[0]) as total flux instead",
  "crpix=\n        Override/Set crpix (1,1,1) // ignored",
  "crval=\n        Override/Set crval (0,0,0) // ignored",
  "cdelt=\n        Override/Set cdelt (1,1,1) // ignored",
  "seed=0\n        Random seed",
  "VERSION=0.6\n   5-jan-05 PJT",
  NULL,
};

string usage = "image creation with objects";

string cvsid = "$Id$";

#ifndef HUGE
# define HUGE 1.0e20
#endif

#define MAXIMAGE 20
#define MAXPAR  20

imageptr iptr;	                /* pointers to (input) image */

#define MAXNAX 3

double crval[MAXNAX], crpix[MAXNAX], cdelt[MAXNAX];
int nwcs = 0;

real spar[MAXPAR];
int npar;

real pa,inc,center[2];
real sinp, cosp, sini, cosi;

bool Qtotflux;
real surface=1.0;


local int set_axis(string var, int n, double *xvar, double defvar);
local void do_create(int nx, int ny, int nz);
local void wcs_f2i(             int ndim, double *crpix, double *crval, double *cdelt, image *iptr);
local void wcs_i2f(image *iptr, int ndim, double *crpix, double *crval, double *cdelt);

local void object_flat(int npars, real *pars);
local void object_exp(int npars, real *pars);
local void object_gauss(int npars, real *pars);
local void object_bar(int npars, real *pars);
local void object_spiral(int npars, real *pars);
local void object_noise(int npars, real *pars);


extern string *burststring(string,string);


void nemo_main ()
{
  string  fnames;
  stream  instr;                      /* input file (optional) */
  stream  outstr;                     /* output file */
  int     size[3],nx, ny, nz;         /* size of scratch map */
  int     ncen;
  string  object;
  int seed = init_xrandom(getparam("seed"));

  object = getparam("object");
  npar = nemoinpr(getparam("spar"),spar,MAXPAR);
  if (npar < 0) error("Syntax error %s",getparam("spar"));

  Qtotflux = getbparam("totflux");

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

  init_xrandom(getparam("seed"));
  nwcs += set_axis(getparam("crpix"),MAXNAX,crpix,1.0);
  nwcs += set_axis(getparam("crval"),MAXNAX,crval,0.0);
  nwcs += set_axis(getparam("cdelt"),MAXNAX,cdelt,1.0);
  if (hasvalue("in")) {
    instr = stropen(getparam("in"),"r");    /* open file */
    read_image (instr, &iptr);
  } else {
    dprintf(0,"Generating a map from scratch\n");
    switch (nemoinpi(getparam("size"),size,3)) {
    case 1:			/*  nx[,nx,1] */
      size[1] = size[0];
      size[2] = 1;
      break;
    case 2:			/*  nx,ny[,1] */
      size[2] = 1;
      break;
    case 3:			/*  nx,ny,nz  */
      break;
    case 0:			/*  [10,10,1] */
      dprintf(0,"Cannot have no size, default 10 assumed\n");
      size[0] = size[1] = 10;
      size[2] = 1;
      break;
    default:			/* --- some error --- */
      error("Syntax error in size keyword\n");
    }
    nx = size[0];
    ny = size[1];
    nz = size[2];
    do_create(nx,ny,nz);
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
  else if (streq(object,"spiral"))
    object_spiral(npar,spar);
  else if (streq(object,"noise"))
    object_noise(npar,spar);
  else
    error("Unknown object %g",object);
  
  outstr = stropen (getparam("out"),"w");  /* open output file first ... */
  write_image (outstr,iptr);         /* write image to file */
  strclose(outstr);
}


local int set_axis(string var, int n, double *xvar, double defvar)
{
  int i, nret;
  if (var == 0 || *var == 0) {
    for (i=0; i<n; i++) 
      xvar[i] = defvar;
    return 0;
  } else {
    nret = nemoinpd(var,xvar,n);
    if (nret < 0) error("Parsing error %d in %s",nret,var);
    for (i=nret; i<n; i++)
      xvar[i] = defvar;
    return 1;
  }
  return 999;
}

/*
 *    FITS: x = (i-crpix)*cdelt + crval        lower/left is 1 (i=1...naxis)
 *    NEMO: x = i*Dx + Xmin                    lower/left is 0 (i=0...naxis-1)
 */

local void wcs_f2i(int ndim, double *crpix, double *crval, double *cdelt, image *iptr)
{
  int i;
  if (ndim<1) return;

  for (i=0; i<ndim; i++)
    dprintf(1,"axis %d: %g %g %g\n",i+1,crpix[i],crval[i],cdelt[i]);
  
  Dx(iptr) = cdelt[0];
  Xmin(iptr) = (1.0-crpix[0])*cdelt[0] + crval[0];
  if (ndim==1) return;
  
  Dy(iptr) = cdelt[1];
  Ymin(iptr) = (1.0-crpix[1])*cdelt[1] + crval[1];
  if (ndim==2) return;

  Dz(iptr) = cdelt[2];
  Zmin(iptr) = (1.0-crpix[2])*cdelt[2] + crval[2];

  dprintf(1,"XYZMin/Dxyz: %g %g %g %g 5g 5g\n",
	  Xmin(iptr),Ymin(iptr),Zmin(iptr),Dx(iptr),Dy(iptr),Dz(iptr));

}

local void wcs_i2f(image *iptr, int ndim, double *crpix, double *crval, double *cdelt)
{
  int i;
  if (ndim<1) return;

  dprintf(1,"XYZMin/Dxyz: %g %g %g %g 5g 5g\n",
	  Xmin(iptr),Ymin(iptr),Zmin(iptr),Dx(iptr),Dy(iptr),Dz(iptr));

  crpix[0] = 1.0;
  crval[0] = Xmin(iptr);
  cdelt[0] = Dx(iptr);
  if (ndim==1) return;
  
  crpix[1] = 1.0;
  crval[1] = Ymin(iptr);
  cdelt[1] = Dy(iptr);
  if (ndim==2) return;

  crpix[2] = 1.0;
  crval[2] = Ymin(iptr);
  cdelt[2] = Dy(iptr);

  for (i=0; i<ndim; i++)
    dprintf(1,"axis %d: %g %g %g\n",i+1,crpix[i],crval[i],cdelt[i]);

}


/*
 *  create new map from scratch, using %x and %y as position parameters 
 *		0..nx-1 and 0..ny-1
 */
local void do_create(int nx, int ny,int nz)
{
    double m_min, m_max, total;
    real   fin[5], fout;
    int    ix, iy, iz;
    int    badvalues;
    
    m_min = HUGE; m_max = -HUGE;
    total = 0.0;		/* count total intensity in new map */
    badvalues = 0;		/* count number of bad operations */

    if (nz > 0) {
      if (!create_cube (&iptr, nx, ny, nz))	/* create default empty image */
        error("Could not create 3D image from scratch");
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
      MapValue(iptr,i,j) += A;
}

local void object_exp(int npars, real *pars)
{
  int i,nx = Nx(iptr);
  int j,ny = Ny(iptr);
  real A = 1.0;
  real h = 1.0;
  real x1,y1,x2,y2,r,arg;

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
      if (arg < 100) MapValue(iptr,i,j) += A * exp(-arg);
    }
  }
}

local void object_gauss(int npars, real *pars)
{
  int i,nx = Nx(iptr);
  int j,ny = Ny(iptr);
  real A = 1.0;
  real h = 1.0;
  real x1,y1,x2,y2,r,arg;

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
      if (arg < 100) MapValue(iptr,i,j) += A * exp(-0.5*arg*arg);
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
  real x1,y1,x2,y2,x3,y3,r,arg;
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
      if (arg < 100) MapValue(iptr,i,j) += A * exp(-arg);
    }
  }
}

local void object_spiral(int npars, real *pars)
{
  int i,nx = Nx(iptr);
  int j,ny = Ny(iptr);
  int l, lmax;
  real A = 1.0;
  real h = 1.0;
  real k = 1.0;   /* wave number */
  real p = 1.0;   /* 1/2 power of cos */
  real r0 = 0.0;  /* starting radius */
  real p0 = 0.0;  /* starting angle */
  real x1,y1,x2,y2,x3,y3,r,arg;
  real amp,phi,value,sum;

  if (npar > 0) A = pars[0];
  if (npar > 1) h = pars[1];
  if (npar > 2) k = pars[2];
  if (npar > 3) p = pars[3];
  if (npar > 4) r0 = pars[4];
  if (npar > 5) p0 = pars[5];
  dprintf(0,"spiral: %g %g %g %g %g %g\n",A,h,k,p,r0,p0);

  if (A==0) return;

  p0 *= PI/180.0;  /* convert from degrees to radians */
  k *= TWO_PI;     /* angles are 2.PI.k.r , so absorb 2.PI.k in one */

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
	if (r < r0) continue;
	/* ? should match this up better so we can connect bar and spiral ? */
	phi = atan2(y2,x2) + k*(r-r0) + p0;
	amp = pow(cos(phi),2*p);
	arg = r/h;
	value = (arg < 100) ? A * amp * exp(-arg) :  0.0;
	if (Qtotflux) {   
	  if (l==0) 
	    sum += value;
	  else
	    MapValue(iptr,i,j) += value;
	} else
	  MapValue(iptr,i,j) += value;
      } /* i */
    } /* j */
  } /* k */
}

local void object_noise(int npars, real *pars)
{
  int i,j;
  int nx = Nx(iptr);
  int ny = Ny(iptr);
  double m = 1.0;
  double s = 1.0;

  if (npar > 0) m = pars[0];
  if (npar > 1) s = pars[1];

  dprintf(0,"noise:%g %g\n",m,s);

  if (m==0) return;

  if (Qtotflux) {
    m /= (nx*ny);
    s /= (nx*ny);
    dprintf(0,"noise: m->%g  s->%g\n",m,s);
  }
  

  for (j=0; j<ny; j++)
    for (i=0; i<nx; i++)
      MapValue(iptr,i,j) += grandom(m,s);
}




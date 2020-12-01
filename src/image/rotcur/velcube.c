/* 
 * VELCUBE: create a cube from a velocity field map, with optional
 *          density and velocity dispersion maps
 *
 *
 *  28-may-04   V1.0  first version (see also Miriad's velimage program)
 *  30-nov-20   V1.1  fix header
 *
 *  TODO:
 *    - normalization of exp(-x^2/2) :
 *         Recall: Integrate[Exp[-x^2],{x,-Infinity,Infinity}] = Pi
 *      is currently enforced by adding the flux and equating it to input
 *      density
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <image.h>

string defv[] = {
  "out=???\n            Output file name (a cube)",
  "invel=???\n          Velocity field",
  "inden=\n             Surface density (optionally)",
  "insig=\n             Velocity dispersion map (optionally)",
  "zrange=-2:2\n        Range in Z velocity",
  "nz=64\n              Number of pixels in Z",
  "sigdefault=0\n       Default velocity dispersion",
  "VERSION=1.1\n        30-nov-2020 PJT",
  NULL,
};

string usage="create a cube from a velocity field";

#define RPD  PI/180
#ifndef HUGE
#define HUGE 1.0e20
#endif

local stream  outstr;                    /* output file */

local int  nz;                     /* 3D size of cube */
local real zrange[3];              /* min, max and cellsize */

local imageptr denptr=NULL;
local imageptr velptr=NULL;
local imageptr sigptr=NULL;


local real undef = 0.0;            /* could also use IEEE NaN ??? set_fblank ? */
local real sigdef = 0.0;

local void cube_create(stream);

extern double grandom(double,double);
void setparams(void);
void setrange(real *rval, string rexp);


void nemo_main(void)
{
    setparams();
    outstr = stropen(getparam("out"),"w");
    cube_create(outstr);
    strclose(outstr);
}

void setparams(void)
{
  stream instr;
  string intpol;
  int i,j,n, nbadvel,nbadsig,nbadden;
  real vel, sig, den, minvel, maxvel, minsig, maxsig;
  real f=4.0;  /* f*sigma is the assumed contributions of signal */
  
  
  nz = getiparam("nz");
  setrange(zrange, getparam("zrange"));
  if (zrange[0] >= zrange[1])
    error("Bad zrange");
  zrange[2] = (zrange[1]-zrange[0])/nz;

  instr = stropen(getparam("invel"),"r");
  read_image(instr,&velptr);
  strclose(instr);

  if (hasvalue("inden")) {
    instr = stropen(getparam("inden"),"r");
    read_image(instr,&denptr);
    strclose(instr);
  }

  if (hasvalue("insig")) {
    instr = stropen(getparam("insig"),"r");
    read_image(instr,&sigptr);
    strclose(instr);
  } else
    sigdef = getdparam("sigdefault");

  n = nbadvel = nbadsig = nbadden = 0;
  for (j=0; j<Ny(velptr); j++)
    for (i=0; i<Nx(velptr); i++) {
      vel = MapValue(velptr,i,j);
      if (denptr) {
	den = MapValue(denptr,i,j);
	if (den < 0) {
	  nbadden++;
	  continue;
	} else if (den == 0)
	  continue;
      } 
      if (sigptr) {
	sig = MapValue(sigptr,i,j);
	if (sig < 0) {
	  nbadsig++;
	  continue;
	}
      } else
	sig = sigdef;
      if (vel-f*sig < zrange[0]) nbadvel++;
      if (vel+f*sig > zrange[1]) nbadvel++;
      if (n==0) {
	minvel = vel-f*sig;
	maxvel = vel+f*sig;
	minsig = maxsig = sig;
      } else {
	minvel = MIN(minvel,vel-f*sig);
	maxvel = MAX(maxvel,vel+f*sig);
	minsig = MIN(minsig,sig);
	maxsig = MAX(maxsig,sig);
      }	
      n++;
    }
  if (nbadden) warning("%d bad density pixels (density < 0)",nbadden);
  if (nbadsig) warning("%d bad sigma pixels (sigma < 0)",nbadsig);
  if (nbadvel) warning("%d/%d pixels outside your zrange %g %g\n",
		       nbadvel,Nx(velptr)*Ny(velptr),zrange[0],zrange[1]);
  dprintf(0,"Minvel=%g Maxvel=%g (using %g sigma) Minsig=%g Maxsig=%g\n",
	  minvel,maxvel,f,minsig,maxsig);
  
}

void setrange(real *rval, string rexp)
{
    char *cptr, *tmpstr;
    double dpar;
                                                                                                               
    cptr = strchr(rexp, ':');
    if (cptr != NULL) {
        tmpstr = allocate(cptr-rexp+1);
        strncpy(tmpstr,rexp,cptr-rexp);
        if (nemoinpd(tmpstr,&dpar,1) != 1)
            error("setrange: parsing error %s",tmpstr);
        free(tmpstr);
        rval[0] = dpar;
                                                                                                               
        if (nemoinpd(cptr+1,&dpar,1) != 1)
            error("setrange: parsing error %s",cptr+1);
        rval[1] = dpar;
    } else {
        rval[0] = 0.0;
        if (nemoinpd(rexp,&dpar,1) != 1)
            error("setrange: parsing error %s",rexp);
        rval[1] = dpar;
    }
    rval[2] = rval[1] - rval[0];
}



/*
 * create a velocity field (GIPSY method)
 *              0..nx-1 and 0..ny-1 
 *      start from pixel, work back to gal plane and interpolate
 *      (a.k.a. retracing method)
 *
 */
local void cube_create(stream outstr)
{
  int  i, j, k, n, nx, ny;
  real m_min, m_max, sum;
  real den, vel, sig, velk, x;
  imageptr vptr;
  real f = 4.0;
  
  m_min = HUGE; m_max = -HUGE;
  nx = Nx(velptr);
  ny = Ny(velptr);
  
  if (!create_cube(&vptr, nx, ny, nz))   /* output data cube */
    error("Could not create cube from scratch");
  
  for (j=0; j<ny; j++)            /* Loop over all pixels */
    for (i=0; i<nx; i++) {
      for (k=0; k<nz; k++)
	CubeValue(vptr,i,j,k) = undef;       /* first set all to 'undefined' */
      sum = 0.0;
      vel = MapValue(velptr,i,j);

      if (denptr) {
	den = MapValue(denptr,i,j);
	if (den <= 0.0) continue;
      } else
	den = 1.0;

      if (sigptr) {
	sig = MapValue(sigptr,i,j);
	if (sig < 0.0) continue;
      } else
	sig = sigdef;
      if (sig == 0.0) {      /* special case, only populate 1 cell */
	k = (int) floor(    (vel-zrange[0])/zrange[2]   );
	if (k<0 || k>nz-1) continue;
	CubeValue(vptr,i,j,k) = den;
	sum = den;
      } else {               /* else use a gaussian profile */
	for (k=0, velk=zrange[0]+0.5*zrange[2]; k<nz; k++, velk += zrange[2]) {
	  x = (velk-vel)/sig;
	  x = ABS(x);
	  if (x > f) continue;
	  CubeValue(vptr,i,j,k) = den * exp(-0.5*x*x);
	  sum += CubeValue(vptr,i,j,k);
	}
      }
      sum = den/sum;              /* now normalize the spectrum so the sum is 'den' */
      for (k=0; k<nz; k++) {
	if (CubeValue(vptr,i,j,k) != undef)
	  CubeValue(vptr,i,j,k) *= sum;
      }
    }
      
  n=0;
  for (k=0; k<nz; k++)        /* get min and max in map */
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++) 
	if (CubeValue(vptr,i,j,k) != undef) {
	  m_min = MIN(CubeValue(vptr,i,j,k),m_min);
	  m_max = MAX(CubeValue(vptr,i,j,k),m_max);
	} else
	  n++;
  MapMin(vptr) = m_min;
  MapMax(vptr) = m_max;
  Dx(vptr) = Dx(velptr);
  Dy(vptr) = Dy(velptr);
  Dz(vptr) = zrange[2];
  Xmin(vptr) = Xmin(velptr);
  Ymin(vptr) = Ymin(velptr);
  Zmin(vptr) = zrange[0];
  Namex(vptr) = strdup("X");
  Namey(vptr) = strdup("Y");
  Namez(vptr) = strdup("Z");
  
  if (n>0) warning("%d/%d cells with no signal",n,nx*ny*nz);
  printf("Min and max in map: %g %g\n",m_min,m_max);
  write_image (outstr,vptr);  
}





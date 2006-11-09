/*
 * MKCOSMO: set up a cube from a density grid, for cosmology
 *	
 *	 1-nov-06  V0.1  Created     - Peter Teuben / Ed Shaya
 *       7-nov-06  V0.4  add rhob=, a=
 *
 * todo:
 *  - first point is always 0,0,0
 *  - do rescale_image()
 *  - implement (1+z) scaling in rescale_image
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

#include <image.h>

string defv[] = {	/* DEFAULT INPUT PARAMETERS */
  "in=???\n       Input density (fluctuation) cube",
  "out=???\n      Output file name",
  "z=\n           Redshift, for growth function 1/(1+x)",
  "D=\n           Growth function value, given explicitly (D<=1)",
  "a=1\n          a scaling?",
  "rhob=\n        Mean density to use",
  "absrho=t\n     Is map absolute density or relative d(rho)/rho?",
  "sigma=0\n      Also perturb distances by gaussian sigma",
  "seed=0\n       Initial seed",
  "rejection=f\n  Use rejection technique to seed the 'grid'",
  "nbody=\n       Use this instead of NX*NY*NZ if rejection is used",
  "headline=\n    Random verbiage",
  "VERSION=0.5\n  9-nov-06 PJT",
  NULL,
};

string usage = "create a cosmology cube of equal massive stars";

string cvsid="$Id$";

local Body *btab;
local int nbody;
local real sigma;
local real cube_aver, cube_max;

local imageptr iptr=NULL;

extern double xrandom(double,double), grandom(double,double);

void check_image(void);
void rescale_image(bool Qabs);
void rescale_image_d(real d, real a, real rb);
void write_snap(string name, string headline);
void mkcube(void), mkcube_reject(void);
void fiddle_x(void),  fiddle_y(void),  fiddle_z(void), drifter(void);

void nemo_main()
{
  int seed;
  stream instr;
  real z = getdparam("z");
  bool Qabs = getbparam("absrho");
  bool Qreject = getbparam("rejection");

  instr = stropen (getparam("in"),"r");
  read_image(instr,&iptr); 
  strclose(instr);      
  dprintf(0,"MinMax = %g %g \n",MapMin(iptr),MapMax(iptr));

  if (hasvalue("D") && hasvalue("rhob") && hasvalue("a"))
    rescale_image_d(getdparam("D"), getdparam("a"), getdparam("rhob"));
  else {
    warning("Currently program needs D= a= and rhob=; doing no scaling now");
  }
  check_image();


  
  if (Nx(iptr) != Ny(iptr) || Nx(iptr) != Nz(iptr))
    warning("Input data is not a cube: %d x %d x %d",Nx(iptr),Ny(iptr),Nz(iptr));
  nbody = Nx(iptr)*Ny(iptr)*Nz(iptr);

  sigma = getdparam("sigma");
  seed = init_xrandom(getparam("seed"));
  
  if (Qreject) {
    if (hasvalue("nbody")) nbody = getiparam("nbody");
    mkcube_reject();
  } else {
    mkcube();

    fiddle_x();
    fiddle_y();
    fiddle_z();
    drifter();
  }

  write_snap(getparam("out"), getparam("headline"));
  free(btab);
}

void check_image(void)
{  
  int ix,iy,iz;
  int nx=Nx(iptr), ny=Ny(iptr), nz=Nz(iptr);
  real sum, cmax;

  /* first get the total number in densities (mass) */

  sum = 0.0;
  cmax = CubeValue(iptr,0,0,0);
  for (iz=0; iz<nz; iz++)
    for (iy=0; iy<ny; iy++)
      for (ix=0; ix<nx; ix++) {
	sum += CubeValue(iptr,ix,iy,iz);
	if (cmax < CubeValue(iptr,ix,iy,iz)) cmax = CubeValue(iptr,ix,iy,iz);
      }
  cube_aver = sum/(nx*ny*nz);
  cube_max  = cmax;
  dprintf(0,"Average pixel value in cube: %g\n",cube_aver);
  dprintf(0,"Max pixel value in cube: %g\n",cube_max);
}


void rescale_image(bool Qabs)
{  
  int ix,iy,iz;
  int nx=Nx(iptr), ny=Ny(iptr), nz=Nz(iptr);
  real sum;

  /* first get the total number in densities (mass) */

  sum = 0.0;
  for (iz=0; iz<nz; iz++)
    for (iy=0; iy<ny; iy++)
      for (ix=0; ix<nx; ix++)
	sum += CubeValue(iptr,ix,iy,iz);

  if (Qabs) {
    /* if Absolute densities: */
    for (iz=0; iz<nz; iz++)
      for (iy=0; iy<ny; iy++)
	for (ix=0; ix<nx; ix++)
	  CubeValue(iptr,ix,iy,iz) /= sum;
  } else {
    /* if Relative densities: */
    if (sum != 0.0) warning("Relative densities, but sum=%g should be 0.0",sum);
    for (iz=0; iz<nz; iz++)
      for (iy=0; iy<ny; iy++)
	for (ix=0; ix<nx; ix++)
	  CubeValue(iptr,ix,iy,iz) += 1.0;
  }
}

/*
 * rescale:  for given growth function D (D=1 if none, and <<1 (typically 1/(1+z))
 *           rhobar = mean density
 *
 */

void rescale_image_d(real d, real a, real rhobar)
{  
  int ix,iy,iz;
  int nx=Nx(iptr), ny=Ny(iptr), nz=Nz(iptr);
  real sum, offset;

  /* first get the total number in densities (mass) */

  sum = 0.0;
  offset = (1/(a*a*a)-d)*rhobar;
  for (iz=0; iz<nz; iz++)
    for (iy=0; iy<ny; iy++)
      for (ix=0; ix<nx; ix++) {
	sum += 	CubeValue(iptr,ix,iy,iz);
	CubeValue(iptr,ix,iy,iz) = d*CubeValue(iptr,ix,iy,iz) + offset;
      }
  sum /= nx*ny*nz;
  dprintf(0,"Average density before scaling=%g\n",sum);
  dprintf(0,"Offset=%g\n",offset);
}


/*
 * WRITEGALAXY: write galaxy model to output.
 */

void write_snap(string name, string headline)
{
  stream outstr;
  real tsnap = 0.0;
  int bits = MassBit | PhaseSpaceBit;
  
  if (! streq(headline, ""))
    set_headline(headline);
  outstr = stropen(name, "w");
  put_history(outstr);
  put_snap(outstr, &btab, &nbody, &tsnap, &bits);
  strclose(outstr);
}

/*
 * MKCUBE: initial homogenous cube based on the input image
 */

void mkcube(void)
{
  Body *bp;
  real x,y,z,mass_i;
  int ix,iy,iz;

  bp = btab = (Body *) allocate(nbody * sizeof(Body));
  mass_i = 1.0/nbody;
  for (iz=0, z=Zmin(iptr); iz<Nz(iptr); iz++, z+=Dz(iptr))
    for (iy=0, y=Ymin(iptr); iy<Ny(iptr); iy++, y+=Dy(iptr))
      for (ix=0, x=Xmin(iptr); ix<Nx(iptr); ix++, x+=Dx(iptr)) {
	Mass(bp) = mass_i;
	Phase(bp)[0][0] = x;
	Phase(bp)[0][1] = y;
	Phase(bp)[0][2] = z;
	Phase(bp)[1][0] = 0.0;
	Phase(bp)[1][1] = 0.0;
	Phase(bp)[1][2] = 0.0;
	if (sigma > 0) {
	  Phase(bp)[0][0] += grandom(0.0,sigma);
	  Phase(bp)[0][1] += grandom(0.0,sigma);
	  Phase(bp)[0][2] += grandom(0.0,sigma);
	}
	bp++;
      }
}

void mkcube_reject(void)
{
  Body *bp;
  real x,y,z,v,mass_i;
  int ix,iy,iz, redo=0;
  real xmin, xmax, ymin, ymax, zmin, zmax;

  warning("using rejection to create a cube"); 

  bp = btab = (Body *) allocate(nbody * sizeof(Body));
  mass_i = 1.0/nbody;

  xmin = Xmin(iptr) - 0.5*Dx(iptr);
  xmax = xmin       + Nx(iptr)*Dx(iptr);
  ymin = Ymin(iptr) - 0.5*Dy(iptr);
  ymax = ymin       + Ny(iptr)*Dy(iptr);
  zmin = Zmin(iptr) - 0.5*Dz(iptr);
  zmax = zmin       + Nz(iptr)*Dz(iptr);

  for (bp=btab; bp<btab+nbody; bp++) {
    Mass(bp) = mass_i;
    while(1) {
      x = xrandom(xmin,xmax);
      y = xrandom(ymin,ymax);
      z = xrandom(zmin,zmax);
      ix = floor((x-xmin)/Dx(iptr));
      iy = floor((y-ymin)/Dy(iptr));
      iz = floor((z-zmin)/Dz(iptr));
      v = xrandom(0,cube_max);
      if (v<CubeValue(iptr,ix,iy,iz)) break;
      redo++;
    }
    Phase(bp)[0][0] = x;
    Phase(bp)[0][1] = y;
    Phase(bp)[0][2] = z;
    Phase(bp)[1][0] = 0.0;
    Phase(bp)[1][1] = 0.0;
    Phase(bp)[1][2] = 0.0;
    if (sigma > 0) {
      Phase(bp)[0][0] += grandom(0.0,sigma);
      Phase(bp)[0][1] += grandom(0.0,sigma);
      Phase(bp)[0][2] += grandom(0.0,sigma);
    }
  }
  dprintf(0,"Redrawn %d times for %d particles\n",redo,nbody);

}

/* 
 * fiddle_x: sweep in X and accumulate mass as such, re-adjust
 *           X position from center of cell as to agree with
 *           the thus obtained cumulative M(<x)
 *
 * fiddle_y: --same, but then in y--
 * fiddle_z: --same, but then in z--
 */


void fiddle_x(void)
{
  int ix,iy,iz,i,nerr=0;
  int nx=Nx(iptr),ny=Ny(iptr),nz=Nz(iptr);
  real dx, xold, xnew, xhit, mold, mnew, mhit, err;
  Body *bp;

  mold = mhit = 0.0;
  xold = 0.0;
  dx = Dx(iptr);
  for (iz=0; iz<nz; iz++) {
    for (iy=0; iy<ny; iy++) {
      xold = Xmin(iptr);
      for (ix=0; ix<nx; ix++) {
	mnew = mold + CubeValue(iptr,ix,iy,iz);
	xnew = xold + dx;
	xhit = xold + dx*(mhit-mold)/(mnew-mold);
	i = ix + ny*(iy+nz*iz);
	bp = btab + i;
	err = (Phase(bp)[0][0]-xhit)/dx;
	if (err<-0.5 || err>0.5) {
	  dprintf(1,"X: %d %d %d -> xold=%g xhit=%g   ****\n",ix,iy,iz,xold,xhit);
	  nerr++;
	} else {
	  dprintf(1,"X: %d %d %d -> xold=%g xhit=%g\n",ix,iy,iz,xold,xhit);
	}
	Phase(bp)[0][0] = xhit;
	mhit += cube_aver;
	mold = mnew;
	xold = xnew;
      }
    }
  }
  if (nerr) warning("Found %d non-linear deviations in X",nerr);
}

void fiddle_y(void)
{
  int ix,iy,iz,i,nerr=0;
  int nx=Nx(iptr),ny=Ny(iptr),nz=Nz(iptr);
  real dy, yold, ynew, yhit, mold, mnew, mhit, err;
  Body *bp;
  
  mold = mhit = 0.0;
  yold = 0.0;
  dy = Dy(iptr);
  for (ix=0; ix<nx; ix++) {
    for (iz=0; iz<nz; iz++) {
      yold = Ymin(iptr);
      for (iy=0; iy<ny; iy++) {
	mnew = mold + CubeValue(iptr,ix,iy,iz);
	ynew = yold + dy;
	yhit = yold + dy*(mhit-mold)/(mnew-mold);
	dprintf(1,"Y: %d %d %d -> yold=%g yhit=%g\n",ix,iy,iz,yold,yhit);
	i = ix + ny*(iy+nz*iz);
	bp = btab + i;
	err = (Phase(bp)[0][1]-yhit)/dy;
	if (err<-0.5 || err>0.5) nerr++;
	Phase(bp)[0][1] = yhit;
	mhit += cube_aver;
	mold = mnew;
	yold = ynew;
      }
    }
  }
  if (nerr) warning("Found %d non-linear deviations in Y",nerr);
}

void fiddle_z(void)
{
  int ix,iy,iz,i,nerr=0;
  int nx=Nx(iptr),ny=Ny(iptr),nz=Nz(iptr);
  real dz, zold, znew, zhit, mold, mnew, mhit, err;
  Body *bp;

  
  mold = mhit = 0.0;
  zold = 0.0;
  dz = Dz(iptr);
  for (iy=0; iy<ny; iy++) {
    for (ix=0; ix<nx; ix++) {
      zold = Zmin(iptr);
      for (iz=0; iz<nz; iz++) {
	mnew = mold + CubeValue(iptr,ix,iy,iz);
	znew = zold + dz;
	zhit = zold + dz*(mhit-mold)/(mnew-mold);
	dprintf(1,"Z: %d %d %d -> zold=%g zhit=%g\n",ix,iy,iz,zold,zhit);
	i = ix + ny*(iy+nz*iz);
	bp = btab + i;
	err = (Phase(bp)[0][2]-zhit)/dz;
	if (err<-0.5 || err>0.5) nerr++;
	Phase(bp)[0][2] = zhit;
	mhit += cube_aver;
	mold = mnew;
	zold = znew;
      }
    }
  }
  if (nerr) warning("Found %d non-linear deviations in Z",nerr);
}

void drifter(void)
{
  warning("no drifting done yet");
}

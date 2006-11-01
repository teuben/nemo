/*
 * MKCOSMO: set up a cube from a density grid, for cosmology
 *	
 *	 1-nov-06  V0.1  Created     - Peter Teuben / Ed Shaya
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
  "z=0.0\n        Size of cube",
  "absrho=t\n     Is map absolute density or relative d(rho)/rho?",
  "sigma=0\n      Perturb distances by gaussian sigma",
  "seed=0\n       Initial seed",
  "headline=\n    Random verbiage",
  "VERSION=0.1\n  1-nov-06 PJT",
  NULL,
};

string usage = "create a cosmology cube of equal massive stars";

local Body *btab;
local int nbody;
local real sigma;

local imageptr iptr=NULL;
local real *mcum;

extern double xrandom(double,double), grandom(double,double);


void nemo_main()
{
  int seed;
  stream instr;
  real z = getdparam("z");
  bool Qabs = getbparam("absrho");

  instr = stropen (getparam("in"),"r");
  read_image (instr,&iptr); 
  strclose(instr);      

  /* if z > 0 rescale by 1+z ; if relative, make map absolute
   *  but after rescaling ?
   */
  
  
  if (Nx(iptr) != Ny(iptr) || Nx(iptr) != Nz(iptr))
    warning("Input data is not a cube: %d x %d x %d",Nx(iptr),Ny(iptr),Nz(iptr));
  nbody = Nx(iptr)*Ny(iptr)*Nz(iptr);
  mcum = (real *) allocate((nbody+1)*sizeof(real));

  sigma = getdparam("sigma");
  seed = init_xrandom(getparam("seed"));
  
  mkcube();
  fiddle_x();
  fiddle_y();
  fiddle_z();
  writegalaxy(getparam("out"), getparam("headline"));
  free(btab);
  free(mcum);
}

/*
 * WRITEGALAXY: write galaxy model to output.
 */

writegalaxy(string name, string headline)
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

mkcube()
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



fiddle_x()
{
  int ix,iy,iz,i;
  int nx=Nx(iptr),ny=Ny(iptr),nz=Nz(iptr);

  mcum[0] = 0.0;
  i = 1;
  
  for (iz=0; iz<nz; iz++)
    for (iy=0; iy<ny; iy++)
      for (ix=0; ix<nx; ix++)
	mcum[i] = mcum[i-1] + CubeValue(iptr,ix,iy,iz);
}

fiddle_y()
{
  int ix,iy,iz,i;
  int nx=Nx(iptr),ny=Ny(iptr),nz=Nz(iptr);

  mcum[0] = 0.0;
  i = 1;
  
  for (ix=0; ix<nx; ix++)
    for (iz=0; iz<nz; iz++)
      for (iy=0; iy<ny; iy++)
	mcum[i] = mcum[i-1] + CubeValue(iptr,ix,iy,iz);
}

fiddle_z()
{
  int ix,iy,iz,i;
  int nx=Nx(iptr),ny=Ny(iptr),nz=Nz(iptr);

  mcum[0] = 0.0;
  i = 1;
  
  for (iy=0; iy<ny; iy++)
    for (ix=0; ix<nx; ix++)
      for (iz=0; iz<nz; iz++)
	mcum[i] = mcum[i-1] + CubeValue(iptr,ix,iy,iz);
}

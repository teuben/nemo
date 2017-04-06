/* 
 * MOM2CUBE: Use moment maps to reconstruct a cube
 *
 *       5-apr-2017  V0.1   drafted for EDGE simulations        PJT
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
  "mom0=???\n         Input mom0 image file",
  "mom1=???\n         Input mom1 image file",
  "mom2=\n            Optional mom2 image file",
  "out=???\n          Output cube",
  "sigma=1.0\n        Optional mom2 (if constant)",
  "cube=\n            Optional cube to describe 3rd axis",
  "z=-1.0,1.0,0.1\n   Zmin,Zmax,dZ (or use cube)",
  "clip=0.0\n         Clipvalue for mom0",
  "norm=t\n           Normalize integral to mom0?",
  "clone=f\n          Clone the cube for output",
  "VERSION=0.4\n      6-apr-2017 PJT",
  NULL,
};

string usage = "Use moments to reconstruct a cube";
string cvsid = "$Id$";



void nemo_main()
{
    stream  mom0str, mom1str, mom2str, cubestr, outstr;
    real    sigma;    /* replaces mom2str */
    real    z[3];     /* replaces cubestr */
    bool    Qnorm, Qclone;
    int     i,j,k,nx,ny,nz;
    real    m0, m1, m2, v, clip, sfactor, sum0, sumc;
    imageptr mom0=NULL, mom1=NULL, mom2=NULL, cube=NULL;

    outstr  = stropen(getparam("out"), "w");
    mom0str = stropen(getparam("mom0"), "r");
    mom1str = stropen(getparam("mom1"), "r");
    if (hasvalue("mom2")) {
      dprintf(0,"Using mom2, skipping sigma\n");
      mom2str = stropen(getparam("mom2"), "r");
      sigma = -1;
    } else
      sigma = getrparam("sigma");
    nz    = nemoinpd(getparam("z"),z,3);    /* or cubestr */
    Qnorm = getbparam("norm");
    Qclone = getbparam("clone");
    clip  = getrparam("clip");

    read_image( mom0str, &mom0);
    read_image( mom1str, &mom1);
    if (sigma < 0) read_image( mom2str, &mom2);

    nx = Nx(mom0);	
    ny = Ny(mom1);

    if (hasvalue("cube")) {
      dprintf(0,"Using cube\n");
      cubestr = stropen(getparam("cube"),"r");
      read_image(cubestr, &cube);
      nz   = Nz(cube);
      z[0] = Zmin(cube);
      z[1] = Zmin(cube) + (nz-1) * Dz(cube);
      z[2] = Dz(cube);
      if (!Qclone) {
	free_image(cube);
	create_cube(&cube,nx,ny,nz);
      }
    } else {
      nz = round((z[1]-z[0])/z[2]);
      create_cube(&cube,nx,ny,nz);
    }
    dprintf(0,"Cube with %d planes from %g to %g in steps %g\n",nz,z[0],z[1],z[2]);

    sumc = sum0 = 0.0;
    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
	m0 = MapValue(mom0,i,j);
	m1 = MapValue(mom1,i,j);
	if (mom2)
	  m2 = MapValue(mom2,i,j);
	else
	  m2 = sigma;
	if (m0 > clip)sum0 += m0;
	if (Qnorm)
	  sfactor = 1.0 / (sqrt(TWO_PI) * m2);
	else
	  sfactor = 1.0;
	for (k=0; k<nz; k++) {
	  if (m0 > clip && m2 > 0) {
	    v = z[0] + k*z[2] - m1;
	    v = 0.5*v*v/(m2*m2);
	    CubeValue(cube,i,j,k) = sfactor * m0 * exp(-v);
	    sumc += CubeValue(cube,i,j,k);
	  } else
	    CubeValue(cube,i,j,k) = 0.0;
	}
      }
    }
    dprintf(0,"MOM0 sum=%g  CUBE sum=%g\n",sum0,sumc*z[2]);
    
    write_image(outstr, cube);
    
}



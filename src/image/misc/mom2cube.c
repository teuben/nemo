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
#include <moment.h>

string defv[] = {
  "mom0=???\n         Input mom0 image file",
  "mom1=???\n         Input mom1 image file",
  "mom2=\n            Optional mom2 image file",
  "out=???\n          Output cube",
  "sigma=1.0\n        Optional mom2 (if constant)",
  "cube=\n            Optional cube to describe 3rd axis",
  "z=-1.0,1.0,0.1\n   Zmin,Zmax,dZ (or use cube)",
  "norm=t\n           Normalization?",
  "clip=0.0           Clipvalue for mom0",
  "VERSION=0.1\n      5-apr-2017 PJT",
  NULL,
};

string usage = "Use moments to reconstruct a cube";
string cvsid="$Id$";

local real peak_spectrum(int n, real *spec, int p);
local real peak_mom(int n, real *spec, int *smask, int peak, int mom);
local real peak_axis(imageptr iptr, int i, int j, int k, int axis);
local int  peak_find(int n, real *data, int *mask, int npeak);
local int  peak_assign(int n, real *data, int *mask);
local bool out_of_range(real);
local void image_oper(imageptr ip1, string oper, imageptr ip2);



void nemo_main()
{
    stream  mom0str, mom1str, mom2str, cubestr, outstr;
    real    sigma;    /* replaces mom2str */
    real    z[3];     /* replaces cubestr */
    bool    Qnorm;
    int     i,j,k,nx,ny,nz;
    real    m0, m1, m2, v;
    

    outstr  = stropen(getparam("out"), "w");
    mom0str = stropen(getparam("mom0"), "r");
    mom1str = stropen(getparam("mom1"), "r");    
    sigma = getrparam("sigma");              /* or mom2str */
    nz    = nemoinpd(getparam("z"),z,3);    /* or cubestr */
    Qnorm = getbparam("norm");
    imageptr mom0=NULL, mom1=NULL, mom2=NULL, cube=NULL;

    read_image( mom0str, &mom0);
    read_image( mom1str, &mom1);
    
    nx = Nx(mom0);	
    ny = Ny(mom1);
    nz = round((z[1]-z[0])/z[2]);
    dprintf(0,"Cube with %d planes from %g to %g in steps %g\n",nz,z[0],z[1],z[2]);

    create_cube(&cube,nx,ny,nz);

    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
	m0 = MapValue(mom0,i,j);
	m1 = MapValue(mom1,i,j);
	if (mom2)
	  m2 = MapValue(mom2,i,j);
	else
	  m2 = sigma;
	for (k=0; k<nz; k++) {
	  v = z[0] + k*z[2] - m1;
	  v = 0.5*v*v/(m2*m2);
	  CubeValue(cube,i,j,k) = m0 * exp(-v); 
	}
      }
    }
    
    write_image(outstr, cube);
    
}



/* 
 * CCDWRITE: write a zero image
 * 
 *    17-sep-2003    Created for tutorials           PJT
 *                      
 */


#include <nemo.h>

#include <image.h>

string defv[] = {
  "out=???\n      Output image file",
  "nx=10\n        Size in X",
  "ny=10\n        Size in Y",
  "nz=1\n         Size in Z",
  "VERSION=1.0\n  17-sep-03 PJT",
  NULL,
};

string usage = "tutorial: write an image";


void nemo_main()
{
  stream  outstr;
  int     nx, ny, nz;  
  int     ix, iy, iz;
  imageptr optr=NULL;        /* pointer to image, needs to be NULL to force new */
  real    tmp, sum = 0.0;

  outstr = stropen(getparam("out"), "w");

  nx = getiparam("nx");
  ny = getiparam("ny");
  nz = getiparam("nz");

  dprintf(0,"Creating image cube (size : %d x %d x %d)\n",nx,ny,nz);

  create_cube(&optr,nx,ny,nz);

  for (iz=0; iz<nz; iz++) {
    for (iy=0; iy<ny; iy++) {
      for (ix=0; ix<nx/2; ix++) {
        CubeValue(optr,ix,iy,iz) = 0.0;
      }
    }
  }

  write_image(outstr, optr);



}

/* 
 * CCDREAD: read an image
 * 
 *    17-sep-2003    Created - for tutorials          PJT
 *                      
 */


#include <nemo.h>

#include <image.h>

string defv[] = {
  "in=???\n       Input image file",
  "VERSION=1.0\n  17-sep-03 PJT",
  NULL,
};

string usage = "tutorial: read an image";


void nemo_main()
{
  stream  instr;
  int     nx, ny, nz;  
  int     ix, iy, iz;
  imageptr iptr=NULL;        /* pointer to image, needs to be NULL to force new */
  real    tmp, sum = 0.0;

  instr = stropen(getparam("in"), "r");

  read_image( instr, &iptr);

  nx = Nx(iptr);	
  ny = Ny(iptr);
  nz = Nz(iptr);
  dprintf(0,"Found image cube (size : %d x %d x %d)\n",nx,ny,nz);

  for (iz=0; iz<nz; iz++) {
    for (iy=0; iy<ny; iy++) {
      for (ix=0; ix<nx/2; ix++) {
	sum += CubeValue(iptr,ix,iy,iz);
      }
    }
  }
  dprintf(0,"Sum = %g\n",sum);
}

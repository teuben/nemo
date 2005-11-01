/* 
 * CCDMERGE: merge all images into a cube
 *	quick and dirty: 1-nov-05
 *                      
 */


#include <stdinc.h>		/* also gets <stdio.h>	*/
#include <getparam.h>
#include <vectmath.h>		/* otherwise NDIM undefined */
#include <filestruct.h>
#include <image.h>

string defv[] = {
  "in=???\n       Input image file",
  "out=???\n      Output file",
  "VERSION=0.1\n  1-nov-05 PJT",
  NULL,
};

string usage = "merge all input images into one big cube - memory intensive";

string cvsid="$Id$";


#define MAXIM   512

void nemo_main()
{
    stream  instr, outstr;
    int     nx, ny, nz;        /* size of scratch map */
    int     nx1, ny1, nz1, nz2;
    int     ni, i, ix, iy, iz, iz1;
    real    dmin, dmax;
    imageptr iptr[MAXIM], optr;        /* pointer to image */
    string  flipmode;

    instr = stropen(getparam("in"), "r");
    outstr = stropen(getparam("out"), "w");

    for (i=0; i<MAXIM; i++) {               /* loop over all to gather data */
      iptr[i] = 0;
      if (read_image( instr, &iptr[i]) == 0) break;
      nx1 = Nx(iptr[i]);	
      ny1 = Ny(iptr[i]);
      nz1 = Nz(iptr[i]);
      dprintf(1,"Image %d: %d x %d x %d\n",i,nx1,ny1,nz1);
      if (i==0) {
	nx = nx1;
	ny = ny1;
	nz = nz1;
	dmin = MapMin(iptr[i]);
	dmax = MapMax(iptr[i]);
      } else {
	if (nx != nx1) error("size nx: %d != %d",nx,nx1);
	if (ny != ny1) error("size ny: %d != %d",ny,ny1);
	nz += nz1;
	dmin = MIN(dmin,MapMin(iptr[i]));
	dmax = MAX(dmax,MapMax(iptr[i]));
      }
    }
    ni = i;
    dprintf(0,"Final cube: %d x %d x %d\n",nx,ny,nz);
    dprintf(0,"Data min/max: %g %g\n",dmin,dmax);
    create_cube(&optr,nx,ny,nz);
    MapMin(optr) = dmin;
    MapMax(optr) = dmax;
    Xmin(optr) = Xmin(iptr[0]);
    Ymin(optr) = Ymin(iptr[0]);
    Zmin(optr) = Zmin(iptr[0]);
    Xref(optr) = Xref(iptr[0]);
    Yref(optr) = Yref(iptr[0]);
    Zref(optr) = Zref(iptr[0]);
    Dx(optr) = Dx(iptr[0]);
    Dy(optr) = Dy(iptr[0]);
    Dz(optr) = Dz(iptr[0]);

    for (i=0, iz=0; i<ni; i++) {       /* grab all data in output cube */
      nz1 = Nz(iptr[i]);
      for (iz1=0; iz1< nz1; iz1++, iz++) {
	for (iy=0; iy<ny; iy++) {
	  for (ix=0; ix<nx; ix++) {
	    CubeValue(optr,ix,iy,iz) = CubeValue(iptr[i],ix,iy,iz1);
	  }
	}
      }
    }

    write_image(outstr, optr);
}

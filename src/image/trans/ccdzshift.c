/* 
 * CCDZSHIFT:   shift profiles in a cube to re-align their zero point
 *
 *     30-nov-2020   0.1    drafted - for a stacking problem
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <mdarray.h>

string defv[] = {
  "in=???\n       Input image cube",
  "map=???\n      Z reference map",
  "out=???\n      Re-aligned output cube",
  "wcs=2\n        0=use 0-based pixel   1=use 1-based pixel  2=use WCS",
  "nearest=t\n    Shift with nearest pixel, or interpolate?",
  "center=-1\n    If given, use this pixel to move the centroids to",
  "VERSION=0.1\n  30-nov-2020 PJT",
  NULL,
};

string usage="re-align a cube Z dimension based on a ZMAP (for stacking)";

string cvsid="$Id$";


void nemo_main(void)
{
  stream  instr, zstr, outstr;
  imageptr iptr = NULL, zptr = NULL;
  real zval, *spec;
  int i,j,k,z,kz, nx,ny,nz, ns;
  int wcs = getiparam("wcs");
  int Qnear = getbparam("nearest");
  int iz2 = getiparam("center");
  
  instr = stropen(getparam("in"), "r");     /* get file name and open file */
  read_image( instr, &iptr);                /* read image */
  strclose(instr);                          /* close file */

  zstr = stropen(getparam("map"), "r");     /* get file name and open file */
  read_image( zstr, &zptr);                 /* read image */
  strclose(zstr);                           /* close file */

  outstr = stropen(getparam("out"),"w");
  
  spec = (real *) allocate(2*Nz(iptr)*sizeof(real));
  nx = Nx(iptr);
  ny = Ny(iptr);
  nz = Nz(iptr);
  if (iz2 < 0)  iz2 = nz/2;
  dprintf(0,"Shifting spectra to z pixel = %d (zval=%g)\n",iz2,  (iz2-Zref(iptr))*Dz(iptr) + Zmin(iptr));
  dprintf(0,"WCS mode %d\n",wcs);

  if (!Qnear) warning("nearest=false interpolation not implemented yet");

  ns = 0;
  for (j=0; j<Ny(iptr); j++)
    for (i=0; i<Nx(iptr); i++) {
      zval = MapValue(zptr,i,j);
      if (wcs==0) {         // 0-based pixels (should be default)
	z = rintl(zval);
	for (k=0; k<Nz(iptr); k++) {
	  spec[k] = CubeValue(iptr,i,j,k);
	  CubeValue(iptr,i,j,k) = 0.0;
	}
	for (k=0; k<Nz(iptr); k++) {
	  kz = k + z - iz2;
	  if (kz < 0 || kz >= Nz(iptr)) continue;
	  CubeValue(iptr,i,j,k) = spec[kz];
	  ns++;
	}
      } else if (wcs==1) {  // 1 based pixels (e.g. if you have FITS based values)
	z = (int) zval;	
      } else if (wcs==2) {  // wcs based (will result in 0 based pixels locally)

	// zval = (k-Zmin(iptr))*Dz(iptr) + Zref(iptr);

        z  = rintl((zval - Zmin(iptr))/Dz(iptr) + Zref(iptr));
	dprintf(1,"i,j,zval,z: %d %d %g %d\n",i,j,zval,z);
	
	for (k=0; k<Nz(iptr); k++) {
	  spec[k] = CubeValue(iptr,i,j,k);
	  CubeValue(iptr,i,j,k) = 0.0;
	}
	for (k=0; k<Nz(iptr); k++) {
	  kz = k + z - iz2;
	  if (kz < 0 || kz >= Nz(iptr)) continue;
	  CubeValue(iptr,i,j,k) = spec[kz];
	  ns++;
	}
      } 
    }

  dprintf(0,"%d/%d pixels shifted\n", ns, nx*ny*nz);

  write_image(outstr, iptr);
  
}

/* 
 * helper macros and functions for PV diagram 
 */

#define MV(i,j)   MapValue(iptr,i,j)


/*
 *  wcsio:  Convert a classic FITS wcs into image(5NEMO) notation, and vice versa
 *
 *    FITS: x = (i-crpix)*cdelt + crval        lower/left is 1 (i=1...naxis)
 *    NEMO: x = i*Dx + Xmin                    lower/left is 0 (i=0...naxis-1)
 *
 *  History:
 *    6-jan-2005    pulled out of ccdmath/ccdgen and stuck in library
 */


#include <stdinc.h>
#include <image.h>


/* for now we use the classes version 1 in NEMO, i.e. crpix=1 */

void wcs_f2i(int ndim, double *crpix, double *crval, double *cdelt, 
	     image *iptr)
{
  int i;
  if (ndim<1) return;

  for (i=0; i<ndim; i++)
    dprintf(1,"axis %d: %g %g %g\n",i+1,crpix[i],crval[i],cdelt[i]);
  
  Dx(iptr) = cdelt[0];
  Xmin(iptr) = (1.0-crpix[0])*cdelt[0] + crval[0];
  if (Axis(iptr)==1) {
    Xref(iptr) = crval[0];
  }
  if (ndim==1) return;
  
  Dy(iptr) = cdelt[1];
  Ymin(iptr) = (1.0-crpix[1])*cdelt[1] + crval[1];
  if (Axis(iptr)==1) {
    Yref(iptr) = crval[1];
  }
  if (ndim==2) return;

  Dz(iptr) = cdelt[2];
  Zmin(iptr) = (1.0-crpix[2])*cdelt[2] + crval[2];
  if (Axis(iptr)==1) {
    Zref(iptr) = crval[2];
  }

  dprintf(1,"XYZMin/Dxyz: %g %g %g %g %g %g\n",
	  Xmin(iptr),Ymin(iptr),Zmin(iptr),Dx(iptr),Dy(iptr),Dz(iptr));

}


void wcs_i2f(image *iptr, 
	     int ndim, double *crpix, double *crval, double *cdelt)
{
  int i;
  if (ndim<1) return;

  dprintf(1,"XYZMin/Dxyz: %g %g %g %g %g %g\n",
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

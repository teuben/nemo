/*
 * write a simple 3D FITS image, taken from the 2D write_image example in 
 * cookbook.c , which comes distributed with cfitsio
 *
 * In NEMO the command "bake writeimage" should suffice to compile/link,
 * although make sure $NEMOLIB/makedefs contains CFITSIO_LIB = -lcfitsio
 *
 * Otherwise something like this will probably work:
 *
 * gcc -I/usr/local/include -L/usr/local/lib \
 *     -o writeimage writeimage.c -lcfitsio -lm
 *
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "fitsio.h"

void printerror( int status);

#define NX 300
#define NY 200
#define NZ 100

int main( void )
{
  fitsfile *fptr;  
  int status, ii, jj, kk;
  long  fpixel, nelements, ncode=10;
  float eps=0.05, tol=0.7, *array[NY];

  char filename[] = "cube.fits";                    /* name for output FITS file */
  int bitpix      =  FLOAT_IMG;                 /* single precision pixel values */
  long naxis      =   3;                             /* 3-dimensional image cube */  
  long naxes[3]   = { NX, NY, NZ };      /* X=300 cols, Y=200 rows, Z=100 planes */
  
  array[0] = (float *)malloc( naxes[0] * naxes[1] * sizeof(float) );  /* a plane */

  for( ii=1; ii<naxes[1]; ii++ )                            /* array of pointers */
    array[ii] = array[ii-1] + naxes[0];

  remove(filename);                      /* Delete old file if it already exists */

  status = 0;                /* initialize status before calling fitsio routines */

  if (fits_create_file(&fptr, filename, &status))        /* create new FITS file */
    printerror( status );

  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    printerror( status );

  fits_write_key_dbl(fptr,"CRPIX1",1.0,10,"X reference pixel",&status);   /* WCS */
  fits_write_key_dbl(fptr,"CRPIX2",1.0,10,"Y reference pixel",&status);
  fits_write_key_dbl(fptr,"CRPIX3",1.0,10,"Z reference pixel",&status);

  fits_write_key_dbl(fptr,"CRVAL1",2.0,10,"X reference value",&status);
  fits_write_key_dbl(fptr,"CRVAL2",1.0,10,"Y reference value",&status);
  fits_write_key_dbl(fptr,"CRVAL3",0.0,10,"Z reference value",&status);

  fits_write_key_dbl(fptr,"CDELT1",1.0,10,"X pixel increment",&status);
  fits_write_key_dbl(fptr,"CDELT2",0.1,10,"Y pixel increment",&status);
  fits_write_key_dbl(fptr,"CDELT3",1.0,10,"Z pixel increment",&status);

  fits_write_key_str(fptr,"CTYPE1","X----CAR","X axis label",&status);
  fits_write_key_str(fptr,"CTYPE2","Y----CAR","Y axis label",&status);
  fits_write_key_str(fptr,"CTYPE3","Z----CAR","Z axis label",&status);

  fits_write_comment(fptr, "Code parameters:", & status);     /* some other pars */
  fits_write_key_flt(fptr, "EPS", eps, 10, "Softening", &status);
  fits_write_key_flt(fptr, "TOL", tol, 10, "Tolerance", &status);
  fits_write_key_lng(fptr, "NCODE", ncode, "Code Version", &status);
  fits_write_comment(fptr, "Note: code is experimental", & status);
  fits_write_history(fptr, "Code History:", & status);
  fits_write_history(fptr, " - none", & status);

  nelements = naxes[0] * naxes[1];                  /* number of pixels to write */

  for (kk = 0; kk < naxes[2]; kk++) { 
    for (jj = 0; jj < naxes[1]; jj++) {        /* create some simple linear ramp */
      for (ii = 0; ii < naxes[0]; ii++) {                  /* foreach plane 'kk' */
	array[jj][ii] = ii + jj + kk;
      }
    }

    fpixel = 1 + kk*nelements;                                         /* offset */

    if ( fits_write_img(fptr, TFLOAT, fpixel, nelements, array[0], &status) )
      printerror( status );
  }
  free( array[0] ); 


  if ( fits_close_file(fptr, &status) ) 
    printerror( status );           

}

void printerror( int status)
{
  if (status) {
    fits_report_error(stderr, status);
    exit( status );  
  }
  return;
}

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
  long  fpixel, nelements, ncode;
  float eps, tol, *array[NY];

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

  if ( fits_write_comment(fptr, "Code parameters:", & status) )
    printerror( status );

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

  eps = 0.05;
  if ( fits_update_key(fptr, TFLOAT, "EPS", &eps, "Softening", &status) )
    printerror( status );           

  tol = 0.7;
  if ( fits_update_key(fptr, TFLOAT, "TOL", &tol, "Tolerance", &status) )
    printerror( status );           
  ncode = 10;
  if ( fits_update_key(fptr, TLONG, "NCODE", &ncode, "Code Version", &status) )
    printerror( status );           

  if ( fits_write_comment(fptr, "Note: code is experimental", & status) )
      printerror( status );           

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

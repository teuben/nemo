/*
 * write a simple 3D image, taken from the 2D write_image example in 
 * cookbook.c , which comes distributed with cfitsio
 *
 * In NEMO the command "bake writeimage" should suffice to compile/link,
 * but otherwise something like this will probably do (where you'll
 * have to find
 *
 * gcc -I/usr/local/include -L/usr/local/lib -o writeimage \
 *     writeimage.c -lcfitsio -lm
 *
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "fitsio.h"

void printerror( int status);


int main( void )
{
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  int status, ii, jj, kk;
  long  fpixel, nelements, ncode;
  float eps, tol, *array[200];

  char filename[] = "cube.fits";                        /* name for new FITS file */
  int bitpix      =  FLOAT_IMG;            /* single precision pixel values       */
  long naxis      =   3;                              /* 2-dimensional image      */    
  long naxes[3]   = { 300, 200, 100 };    /* X=300 wide, Y=200 rows, Z=100 planes */
  
  array[0] = (float *)malloc( naxes[0] * naxes[1] * sizeof(float) );

  for( ii=1; ii<naxes[1]; ii++ )
    array[ii] = array[ii-1] + naxes[0];

  remove(filename);               /* Delete old file if it already exists */

  status = 0;         /* initialize status before calling fitsio routines */

  if (fits_create_file(&fptr, filename, &status)) /* create new FITS file */
    printerror( status );

  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    printerror( status );

  nelements = naxes[0] * naxes[1];          /* number of pixels to write */

  for (kk = 0; kk < naxes[2]; kk++) { 
    for (jj = 0; jj < naxes[1]; jj++) {      /* create some simple linear ramp */
      for (ii = 0; ii < naxes[0]; ii++) {    /* foreach plane 'kk' */
	array[jj][ii] = ii + jj + kk;
      }
    }

    fpixel = 1 + kk*nelements;               /* offset */

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
  
  if ( fits_close_file(fptr, &status) ) 
    printerror( status );           

}

void printerror( int status)
{
  if (status) {
    fits_report_error(stderr, status); /* print error report */
    exit( status );    /* terminate the program, returning error status */
  }
  return;
}

/*
 * wrapper so that NEMO can call cfitsio, since both are using <fitsio.h>
 *
 *   - this file is not used yet, using <cfitsio/fitsio.h> is the current
 *     proposed solution
 */

#include <stdio.h>
#include <fitsio.h>


static fitsfile *fptr;

void ccd_cfitsio_open(char *filename)
{
  int status;
  int bitpix = USHORT_IMG;
  long naxis = 2;
  long naxes[2] = {128, 128};

  if (fits_create_file(&fptr, filename, &status))
    printf("%s\n",status);

  if (fits_create_img(fptr, bitpix, naxis, naxes, &status))
    printf("%s\n",status);
}

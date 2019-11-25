/*
 * sample program to use some basic routines from cfitsio
 *
 * 12-sep-01     written
 * 29-jul-04     added to NEMO
 */
#include <nemo.h>

string defv[] = {
    "out=???\n      Output fits file",
    "VERSION=0.2a\n 23-nov-2019 PJT",
    NULL,
};

string usage = "Trial with cfitsio";

#if defined(HAVE_LIBCFITSIO)
# include <fitsio.h>  
# include <longnam.h>
#else
#error CFITSIO not available or properly configured
#endif

int nemo_main(void)
{
    fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
    int status, ii, jj;
    long  fpixel = 1, naxis = 2, nelements, exposure;
    long naxes[2] = { 300, 200 };   /* image is 300 pixels wide by 200 rows */
    short array[200][300];
    string fname = getparam("out");
    float time4 = PI;
    double time8 = PI;

    status = 0;         /* initialize status before calling fitsio routines */
    fits_create_file(&fptr, fname, &status);   /* create new file */

    /* Create the primary array image (16-bit short integer pixels */
    fits_create_img(fptr, SHORT_IMG, naxis, naxes, &status);

    /* Write a keyword; must pass the ADDRESS of the value */
    exposure = 1500.;
    fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,
         "Total Exposure Time", &status);
    dprintf(1,"exposure %d\n",status);
    fits_update_key(fptr, TFLOAT, "TIME4", &time4,
		    "Sample REAL*4", &status);    
    dprintf(1,"time4    %d\n",status);
    fits_update_key(fptr, TDOUBLE, "TIME8", &time8,
		    "Sample REAL*8", &status);
    dprintf(1,"time8    %d\n",status);
    fits_update_key(fptr, TSTRING, "DATEOBS", "2001/09/30",
		    "Sample STRING", &status);
    dprintf(1,"string   %d\n",status);

    /* Initialize the values in the image with a linear ramp function */
    for (jj = 0; jj < naxes[1]; jj++)
        for (ii = 0; ii < naxes[0]; ii++)
            array[jj][ii] = ii + jj;

    nelements = naxes[0] * naxes[1];          /* number of pixels to write */

    /* Write the array of integers to the image */
    fits_write_img(fptr, TSHORT, fpixel, nelements, array[0], &status);

    fits_close_file(fptr, &status);            /* close the file */

    fits_report_error(stderr, status);  /* print out any error messages */
    return status ;
}

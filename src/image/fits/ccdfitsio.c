/* 
 *      CCDFITSIO: convert CCD image to a fits file, using cfitsio library
 *
 *      26-feb-01  PJT   0.1 quick and dirty compilation example
 *                           DOESN'T WORK .... file_size() conflict with cfitsio library
 *       8-dec-01        using HAVE_LIBCFITSIO, but the image is written transposed
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <history.h>

#ifdef HAVE_LIBCFITSIO
#include <cfitsio/fitsio.h>
#endif

string defv[] = {
        "in=???\n        Input image filename",
        "out=???\n       Output fits filename",
	"bitpix=-64\n    Bitpix",
        "VERSION=0.1\n   8-nov-01 PJT",
        NULL,
};
string usage = "convert image to a fits file using cfitsio";


#ifdef HAVE_LIBCFITSIO

stream  instr;    /* file streams */
fitsfile *fptr;   /* output file */


imageptr iptr=NULL;                     /* image, allocated dynamically */

double scale[3];        /* scale conversion for FITS (CDELT) */
double iscale[2];	/* intensity rescale */
char *object;           /* name of object in FITS header */
char *comment;          /* extra comments */
char *headline;         /* optional NEMO headline, added as COMMENT */
bool Qcdmatrix;         /* writing out new-style cdmatrix ? */
bool Qradecvel;         /* fake astronomy WCS header */

void setparams(void);
void write_fits(string,imageptr);
void fits_error( int status);


void nemo_main()
{
    setparams();                               /* set cmdln par's */
    instr = stropen (getparam("in"), "r");     /* open image file */
    read_image (instr,&iptr);                  /* read image */
    headline = ask_headline();                 /* possible headline */
    strclose(instr);                           /* close image file */

    write_fits(getparam("out"),iptr);          /* write fits file */
    free_image(iptr);
}

void setparams(void)
{
}

void write_fits(string name,imageptr iptr)
{
    int  nx, ny, nz, i, j, k, bitpix, npix;
    char *cp;
    string *hitem;
    float *buffer, *bp;
    int ndim;
    long naxis[3];   /* at most 3D cubes for now */
    int status;
    double bscale, bzero;
    
    nx = naxis[1] = Nx(iptr);
    ny = naxis[0] = Ny(iptr);
    nz = naxis[2] = Nz(iptr);   if (nz <= 0) nz = 1;
    npix = nx*ny*nz;


    dprintf(1,"NEMO Image file written to FITS disk file\n");

    ndim = (nz > 1 ? 3 : 2);    /* Make it 2D or real 3D map/cube */
    bitpix = getiparam("bitpix");

    status = 0;
    fptr = 0;
    if (fits_create_file(&fptr, name, &status))
      fits_error(status);
    if (fits_create_img(fptr,bitpix,ndim,naxis, &status))
      fits_error(status);
    if (ffpprd(fptr,0,1,npix, Frame(iptr), &status))
      fits_error(status);
    if (fits_close_file(fptr, &status))
      fits_error(status);
}

void fits_error( int status)     
{
  /*****************************************************/
  /* Print out cfitsio error messages and exit program */
  /*****************************************************/

  char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];

  if (status==0) return;                  /* 0=ok, no need to do more checking */
  
  fits_get_errstatus(status, status_str);   /* get the error description */
  warning("fits_error: status = %d: %s\n", status, status_str);


  if ( fits_read_errmsg(errmsg) ) {          /* get first message; null if stack is empty */
    dprintf(0, "Error message stack:\n");
    dprintf(0, " %s\n", errmsg);
    while ( fits_read_errmsg(errmsg) )  /* get remaining messages */
      dprintf(0, " %s\n", errmsg);
  }
  stop(status);       /* terminate the program, returning error status to parent */
}



#else
nemo_main()
{
    error("CFITSIO library not available");
}
#endif

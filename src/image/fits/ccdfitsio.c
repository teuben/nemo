/* 
 *      CCDFITSIO: convert CCD image to a fits file, using cfitsio library
 *
 *      26-feb-10  PJT   0.1 quick and dirty compilation example
 *                           DOESN"T WORK .... file_size() conflict with cfitsio library
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <history.h>


#include "/usr/local/include/fitsio.h"   /* cfitsio */




string defv[] = {
        "in=???\n        Input image filename",
        "out=???\n       Output fits filename",
	"bitpix=-32\n    Bitpix",
        "VERSION=0.1\n   26-feb-01 PJT",
        NULL,
};

string usage = "convert image to a fits file using cfitsio";

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
    int  nx, ny, nz, i, j, k, bitpix;
    char *cp;
    string *hitem;
    float *buffer, *bp;
    int ndim, naxis[3];   /* at most 3D cubes for now */
    int status;
    double bscale, bzero;
    
    nx = naxis[0] = Nx(iptr);
    ny = naxis[1] = Ny(iptr);
    nz = naxis[2] = Nz(iptr);   if (nz <= 0) nz = 1;


    dprintf(1,"NEMO Image file written to FITS disk file\n");

    ndim = (nz > 1 ? 3 : 2);    /* Make it 2D or real 3D map/cube */
    bitpix = getiparam("bitpix");

    if (fits_create_file(&fptr, name, &status))
      warning("cfitsio: %d",status);
    if (fits_create_img(fptr,bitpix,ndim,naxis, &status))
      warning("cfitsio: %d",status);
    /* should write the array now */
    if (fits_close_file(fptr, &status))
      warning("cfitsio: %d",status);
}


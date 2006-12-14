/*
   Wrapper routines for the public routines in the NEMO fitsio.c file,
   which call the CFITSIO library.
   This code also includes the old MIRIAD-style code if CFITSIO is not
   available

   16-dec-2001     Original version sketched                 Bill Pence
   18-dec-2001     finalized with the new fitsio_nemo.h      Peter Teuben
   19-dec-2001     shift over comment/history cards          PJT
*/

#include <nemo.h>
#include "fitsio_nemo.h" 

/*#ifndef HAVE_LIBCFITSIO*/
#if 1
#include "fitsio.c"       /* old self-coded MIRIAD-style interface */
#else

/* there is some old outstanding bug in the cfitsio interface... 
 * so as of 14-dec-2006 this version has been hardcoded out
 */

static int w_bitpix = -32;               /* see: fit_setbitpix()    */
static FLOAT w_bscale = 1.0;             /* see: fit_setscale()     */
static FLOAT w_bzero = 0.0;              /* see: fit_setscale()     */
static char *cvs_id="$Id$ fitsio_nemo.c";

/**********************************************************************/
FITS *fitopen (char *name, char *status, int naxis, int *nsize)
/*
  This opens a FITS file and readies it for i/o.

  Inputs:
    name        A string giving the name of the file to be opened.
    status      Either "old" or "new".
    naxis       The number of dimensions of the image. When opening
                an "old" file, fitopen makes sure that the number of
                dimensions does not exceed this. For a "new" file,
                the output file is created with this many dimensions.
  Input or Output:
    nsize       This is input for a "new" file, and output for an "old"
                file. It is of size "naxis". For "new" files it gives
                the number of pixels along each dimension. For "old"
                files, it returns the number of pixels along each
                axis. On output, it is filled with 1s for dummy axes.
  Output:
    fitopen     This is a pointer to a FITS structure
                which describes the file. This pointer is used in all
                subsequent calls to the fitsio routines.
----------------------------------------------------------------------*/
{
    FITS *fptr;
    int tstatus = 0, ndims, hdutype, ii;
    long naxes[MAXNAX];
    
    nemo_dprintf(1,"[fitopen: using CFITSIO wrapper routine]\n");

    if (!strncmp(status, "old", 3) || *status == 'r') { 
         /* open existing image */

         fptr = 0;
         if (fits_open_file(&fptr, name, READONLY, &tstatus) )
             return 0;

         fits_get_hdu_type (fptr, &hdutype, &tstatus);
         if (hdutype != IMAGE_HDU)
             tstatus = 1;  /* not an image HDU */

         fits_get_img_dim(fptr, &ndims, &tstatus);
         if (ndims > naxis)
	   tstatus = 1;   /* no space for all these dimensions */

         fits_get_img_size(fptr, MAXNAX, naxes, &tstatus);
         for (ii = 0; ii < naxis; ii++) {
	   if (ii < ndims)
	     nsize[ii] = naxes[ii];
	   else
	     nsize[ii] = 1;
         }

    }  else if (strstr(status,"new") || *status == 'w') {
        /* create new FITS image file */
        fptr = 0;
        if (fits_create_file(&fptr, name, &tstatus) )
            return(0);

        for (ii = 0; ii < naxis; ii++) 
            naxes[ii] = nsize[ii];   /* copy from int to long array */

        fits_create_img(fptr, w_bitpix, naxis, naxes, &tstatus);
    }

    if (tstatus)  {
        fits_close_file(fptr, &tstatus);
        return 0;      /* error occured */
    }
    else
         return(fptr);    /* OK */
}

/**********************************************************************/
void fitclose (FITS *fptr)
/*
  This closes a FITS file, and deletes any memory associated with it.

  Input:
    file        This is the pointer returned by fitopen.
*/
{
    int status = 0;
    fits_close_file(fptr, &status);
    return;
}
/**********************************************************************/
void fitsetpl (FITS *file, int n, int *nsize)
{
  int i;
  for (i=0; i<n; i++)
    warning("skipping fitsetpl(%d)  nsize[%d]=%d",n,i,nsize[i]);
}
/**********************************************************************/
void fitread(FITS *file, int j, FLOAT *data)
/*
  This reads a row of a FITS image.

  Input:
    file        The pointer to the data structure returned by the fitopen
                routine.
    j           The row number to be read. This varies from 0 to naxis2-1.
  Output:
    data        A FLOAT array of naxis1 elements, being the pixel values
                read.
*/
{
/*  NOTE:  This assumes a 2D image */

    int naxis, status = 0;
    long naxes[MAXNAX], fpixel[MAXNAX], nelements;

    fits_get_img_param(file, MAXNAX, NULL, &naxis, naxes, &status);
    nelements = naxes[0];  /* number of pixels in a row */
    fpixel[0] = 1;
    fpixel[1] = j+1;

    fits_read_pix(file, TFLOAT, fpixel, nelements, NULL, data, 
                   NULL, &status);
    return;
}
/**********************************************************************/
void fitwrite(FITS *file, int j, FLOAT *data)
/*
  This writes a row of a FITS image. 
                    
  Inputs:
   file        The pointer returned by fitopen.
   j           The row number to be written. This varies from 0 to naxis2-1.
   data        A FLOAT array of naxis1 elements, being the pixel values  
                to write.
*/
{
/*  NOTE:  This assumes a 2D image */

    int naxis, status = 0;
    long naxes[MAXNAX], fpixel[MAXNAX], nelements;

    fits_get_img_param(file, MAXNAX, NULL, &naxis, naxes, &status);
    nelements = naxes[0];  /* number of pixels in a row */
    fpixel[0] = 1;   /* cfitsio is 1-based, we are 0-based .... */
    fpixel[1] = j+1;

    fits_write_pix(file, TFLOAT, fpixel, nelements, data, &status);
    return;
}
/**********************************************************************/
void fitrdhdr (FITS *file, char *keyword, FLOAT *value, FLOAT def)
/*
  This reads the value of a real-valued FITS keyword from the file header.

  Input:
    file        The pointer returned by fitopen.
    keyword     The keyword to search for.
    def         If the keyword is not found, this "default value" is
                returned.
  Output:
    value       The value read from the FITS header.
*/
{
    int status = 0;

    if (fits_read_key(file, TFLOAT, keyword, value, NULL, &status) )
        *value = def;   /* keyword not found, so return default value */

    return;
}
/**********************************************************************/
void fitrdhdi (FITS *file, char *keyword, int *value, int def)
/*
  This reads the value of a integer-valued FITS keyword from the file header.

  Input:
    file        The pointer returned by fitopen.
    keyword     The keyword to search for.
    def         If the keyword is not found, this "default value" is
                returned.
  Output:
    value       The value read from the FITS header.
*/
{
    int status = 0;

    if (fits_read_key(file, TINT, keyword, value, NULL, &status) )
        *value = def;   /* keyword not found, so return default value */

    return;
}
/**********************************************************************/
void fitrdhda (FITS *file, char *keyword, char *value, char *def)
/*
  This reads the value of a character-valued FITS keyword from the header.

  Input:
    file        The pointer returned by fitopen.
    keyword     The keyword to search for.
    def         If the keyword is not found, this "default value" is
                returned.
  Output:
    value       The value read from the FITS header.
*/
{
    int status = 0;

    if (fits_read_key(file, TSTRING, keyword, value, NULL, &status) )
        strcpy(value, def);   /* keyword not found, so return default value */

    return;
}
/**********************************************************************/
void fitwrhdr(FITS *file, char *keyword, FLOAT value)
/*
  This writes the value of a real-valued FITS keyword to the file header.

  Input:
    file        The pointer returned by the fitopen routine.
    keyword     The name of the keyword to write.
    value       The value of the keyword.
*/
{
    int status = 0;

    /* write 9 decimals of precision */
    fits_write_key_flt(file, keyword, value, 9, NULL, &status);
    return;
}

/**********************************************************************/
void fitwrhdi(FITS *file, char *keyword, int value)
/*
  This writes the value of a integer-valued FITS keyword to the file header.

  Input:
    file        The pointer returned by the fitopen routine.
    keyword     The name of the keyword to write.
    value       The value of the keyword.
*/
{
    int status = 0, temp;

    temp = value;
    fits_write_key(file, TINT, keyword, &temp, NULL, &status);
    return;
}
/**********************************************************************/
void fitwrhdl(FITS *file, char *keyword, int value)
/*
  This writes the value of a logical-valued FITS keyword to the file header.

  Input:
    file        The pointer returned by the fitopen routine.
    keyword     The name of the keyword to write.
    value       The value of the keyword.
*/
{
    int status = 0, temp;

    temp = value;
    fits_write_key(file, TLOGICAL, keyword, &temp, NULL, &status);
    return;
}
/**********************************************************************/
void fitwrhda(FITS *file, char *keyword, char *value)
/*
  This writes the value of a character-valued FITS keyword to the file header.

  Input:
    file        The pointer returned by the fitopen routine.
    keyword     The name of the keyword to write.
    value       The value of the keyword.
*/
{
    int status = 0;

    fits_write_key(file, TSTRING, keyword, value, NULL, &status);
    return;
}
/**********************************************************************/
void fitwrhd  (FITS *file, char *keyword, char *value)
/*
  This writes a character-value to a FITS keyword with '='
*/
{
    char tmp[100], card[100];
    int status = 0, dummy;

    strcpy(tmp, keyword);   /* construct template */
    strcat(tmp, " ");
    strncat(tmp, value, 70);
    fits_parse_template(tmp, card, &dummy, &status); /* make formated card */
    fits_write_record(file, card, &status);
    return;
}
/**********************************************************************/
void fitwra(FITS *file, string keyword, string value)
/*
  This writes a character-value to a FITS keyword without '='
  as e.g. used by COMMENT and HISTORY.
*/
{
    char tmp[100], card[100];
    int status = 0, dummy;
#if 0
    strcpy(tmp, keyword);   /* construct template */
    strcat(tmp, " ");
    strncat(tmp, value, 70);
    fits_parse_template(tmp, card, &dummy, &status); /* make formatted card */
    card[8] = ' ';   /* erase the '=' from the card, if present */
#else
    if ((int)strlen(value)>71) {
      strncpy(tmp,value,80);  /* chop */
      sprintf(card,"%-8s %-71s",keyword,tmp);
    } else
      sprintf(card,"%-8s %s",keyword,value);
#endif
    fits_write_record(file, card, &status);
    return;
}
/**********************************************************************/

void fit_setbitpix(int bp)
/*
    fit_setbitpix:  set the bitpix value for subsequent output
    Note: this routine must be called before fitopen

    Input:  
        bp      bitpix value; supported are: 8, 16, 32, -32, -64
*/
{
    switch(bp) {
      case -64:
      case -32:
      case  64:
      case  32:
      case  16:
      case   8:
        dprintf(1,"BITPIX reset to %d\n",bp);
        w_bitpix = bp;    break;
      default:
        warning("illegal value of bitpix, %d taken", w_bitpix);
    }
}
/**********************************************************************/
void fit_setscale(FLOAT bs, FLOAT bz)
/*
    fit_setscale:   set bscale and bzero scaling factor for subsequent
                    output. Normally only used for positive bitpix,
                    but is also supported for negative bitpix.
    Note: this routine must be called before fitopen
    
    Input:
        bs      bscale factor
        bz      bzero factor

    (bscale,bzero) such that
                image_value = fits_value * bscale + bzero
*/
{
    w_bscale = bs;
    w_bzero = bz;
}
/**********************************************************************/
void fit_setblocksize(int n)
/*
    fit_setblocksize:  set a new blocksize
                    Although the default is 2880 for standard fits,
                    some files are now written with a blocking factor
                    > 1, in which case the blocksize should be 
                    2880*blocking_factor.
    
    Input:
        n       blocksize to be set
*/
{
    return;   /* this function is not needed with CFITSIO */
}
/**********************************************************************/
int  fitexhd (FITS *fptr, char *keyword)
/*
  Public routine to check for a FITS keyword existence in a file.
  Returns 1 if found, 0 if not.
*/
{
    char card[81];
    int status = 0;

    if (fits_read_card(fptr, keyword, card, &status) )
        return 0;   /* keyword not found */
    else
        return 1;   /* keyword does exist */
}
/**********************************************************************/
#endif

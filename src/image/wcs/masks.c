/** @file masks.c
 *  @brief Read fits images and queries masking
 *
 * A set of routines to read in a set of -presumably related- FITS 
 *   images and query if a certain pixel or (ra,dec) location 
 *   is masked or not.
 *
 * requires: nemo, cfitsio
 *
 * public routines:
 *        void ini_mask(char *filename, int slot);
 *        int  is_mask(double ra_deg, double dec_deg, int slot);
 *        int  is_maski(int ix, int iy, int slot);
 *        void end_mask(int slot);
 *
 * NOTE that is_mask() returns 1 if the pixel is good!
 *
 *  12-jun-2003    Initial version for the new YIELD bandmerging program - PJT
 *  12-sep-2003    Added is_edge routine, similar to is_mask - NLC
 *  19-feb-2004    more verbose error messages
 */

#include <nemo.h>
#include "fitsio.h"

/* -------------------------------------------------------------------------------- */

/**
 *  @def MAXSLOT
 *  @brief number of image slots this program can hold
 *
 *  @todo why do i see MAXSLOT+1 ?
 */

#define MAXSLOT 16

/**
 * @typedef axis 
 * @brief structure that hold information on a single 'fits'-like axis 
 * 
 * An 'axis' hold information for typical astronomical axes, as found
 * in for example FITS files.
 */

typedef struct {           
  int n;                        /**< length */
  double crpix;                 /**< reference pixel (does not have to be 1..n) */
  double crval;                 /**< value at the reference pixel */
  double cdelt;                 /**< pixel increment (at the reference pixel) */
  char ctype[16];               /**< axis type, e.g. 'RA---TAN' */
} axis;

/**
 * @typedef mask_image
 * @brief structure holds all info for an image 
 *
 * mask_image 
 *
 * @todo  explain the rot (CROTA?) stuff, and what about cd matrix
 */

typedef struct {      
  int used;            /**< 0=free 1=used */
  char *filename;      /**< filename where the image came from (NOT USED) */
  float **data;        /**< the image data, structured as "data[ny][nx]"   */
  float *buffer;       /**< buffer[nx*ny] */
  char **cdata;        /**< the mask itself */
  char *cbuffer;       /**< the mask itself */
  axis  x;             /**< X axis info */
  axis  y;             /**< Y axis info */
  double rot;          /**< rotation angle (HOW DEFINED,is this like CROTA??) */
  char proj[10];       /**< projection type for wcs() routines, e.g. '-TAN' */
} mask_image;

static mask_image mask[MAXSLOT+1];           /**< a simple array of masks */
static int nmask = 0;                        /**< counter how many masks used */

/* -------------------------------------------------------------------------------- */

/** @fn static void check_slot(int slot) 
 *
 * check_slot:  should be called by *every* routine, for sanity checking 
 *
 *  @param slot  slot number (0..MAXSLOT-1)
 *  @warning    should be called by *every* routine, for sanity checking 
 */

static void check_slot(int slot) 
{
  int i;

  if (slot < 0 || slot > MAXSLOT)
    error("Illegal slot number %d\n",slot);

  if (nmask == 0) {             /* initialize 'used' slots */
    for (i=0; i<=MAXSLOT; i++)
      mask[i].used = 0;
  }

  return;
}

/** @fn static int get_slot(int slot) 
 *
 */

static int get_slot(int slot) 
{
  check_slot(slot);
  return mask[slot].used;
}

/** @fn static void set_slot(int slot) 
 */
static void set_slot(int slot) 
{
  check_slot(slot);
  if (mask[slot].used == 0) nmask++;
  mask[slot].used = 1;
  
}


/** @fn static void fits_error( int status, string mesg)
 */

static void fits_error(int status, string msg)
{
  char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
                                                                                
  if (status==0) return;                  /* 0=ok, no need to do more checking */

  fits_get_errstatus(status, status_str);   /* get the error description */
  warning("fits_error: status = %d msg = %s: %s\n", status, msg, status_str);


  if ( fits_read_errmsg(errmsg) ) {          /* get first message; null if stack is empty */
    nemo_dprintf(0, "Error message stack:\n");
    nemo_dprintf(0, " %s\n", errmsg);
    while ( fits_read_errmsg(errmsg) )  /* get remaining messages */
      nemo_dprintf(0, " %s\n", errmsg);
  }
  stop(status);       /* terminate the program, NEMO style, return error status to parent */
}



/* -------------------------------------------------------------------------------- */

/** @fn ini_mask(char *filename, int slot)
 *  @brief initializes the mask
 *
 *  ini_mask reads in a FITS file using CFITSIO routines
 *  the data is actually not used and discarded in this
 *  routine. only the mask is retained.
 *
 *  @param filename    filename of the FITS file
 *  @param slot        requested slot
 */


void ini_mask(char *filename, int slot)
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status,  nfound, anynull, nulls;
    long naxes[2], fpixel, npixels, ii;
    float datamin, datamax, *buffer, **data;
    double ra_min, ra_max, dec_min, dec_max;
    char *cbuffer, **cdata;
    mask_image *mp;
    axis *ra, *dec;

    check_slot(slot);
    if (get_slot(slot))
      warning("Reusing slot %d with %s",slot,filename);

    status = 0;
    fits_open_file(&fptr, filename, READONLY, &status);
    fits_error(status,"open_file");

    mp = &mask[slot];
    ra  = &mp->x;
    dec = &mp->y;


    fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);
    fits_error(status,"read_keys_lng");
    /* this routine would crash and burn if the CD matrix is used ... */
    /* hence the special check for error code 506 (APPROX_WCS_KEY)    */
    /* it would be better to have an opaque structure to handle the   */
    /* world_to_pix and pix_to_world conversions  " EXACTLY "         */
    fits_read_img_coord(fptr,
			&mp->x.crval, &mp->y.crval,
			&mp->x.crpix, &mp->y.crpix,
			&mp->x.cdelt, &mp->y.cdelt,
			&mp->rot, mp->proj, &status);
    if (status) {
      if (status == APPROX_WCS_KEY) {    /* for now, ignore minor skews in CD matrix */
	warning("Ignoring minor WCS skew in %s due to CD matrix",filename);
	status = 0;
      } else
	fits_error(status,"read_img_coord");
    }
    ra->n = naxes[0];
    dec->n = naxes[1];

    if (ra->cdelt < 0) {    /* the normal astronomical convention */
      ra_min = (ra->n - ra->crpix)*ra->cdelt/cos(PI/180*dec->crval) + ra->crval;
      ra_max = (1     - ra->crpix)*ra->cdelt/cos(PI/180*dec->crval) + ra->crval;
    } else {
      ra_min = (1     - ra->crpix)*ra->cdelt/cos(PI/180*dec->crval) + ra->crval;
      ra_max = (ra->n - ra->crpix)*ra->cdelt/cos(PI/180*dec->crval) + ra->crval;
    }
    dec_min = (1      - dec->crpix)*dec->cdelt + dec->crval;
    dec_max = (dec->n - dec->crpix)*dec->cdelt + dec->crval;
    dprintf(0,"RA:  %lf - %lf [%d]\n", ra_min,  ra_max,  mp->x.n);
    dprintf(0,"DEC: %lf - %lf [%d]\n", dec_min, dec_max, mp->y.n);
    dprintf(0,"ROT: %lf\n",mp->rot);
    dprintf(0,"PROJ: %s\n",mp->proj);

    npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
    fpixel   = 1;                          /* first pixel */
    datamin  = 1.0E30;
    datamax  = -1.0E30;
    mask[slot].buffer = (float *) allocate(npixels*sizeof(float));   /* temp */
    mask[slot].data   = (float **)allocate(naxes[1]*sizeof(float *));/* temp */
    mask[slot].cbuffer= (char *)  allocate(npixels*sizeof(char));
    mask[slot].cdata  = (char **) allocate(naxes[1]*sizeof(char *));

    buffer = mask[slot].buffer;
    data   = mask[slot].data;
    cbuffer = mask[slot].cbuffer;
    cdata   = mask[slot].cdata;

    for (ii = 0; ii < naxes[1]; ii++) {
      mask[slot].data[ii] = buffer;
      mask[slot].cdata[ii] = cbuffer;
      buffer += naxes[0];
      cbuffer += naxes[0];
    }
    buffer = mask[slot].buffer;
    cbuffer = mask[slot].cbuffer;

    fits_read_imgnull(fptr, TFLOAT, fpixel, npixels, buffer, cbuffer, &anynull, &status);
    fits_error(status,"read_imgnull");

#if 1
    nulls = 0;
    for (ii = 0; ii < npixels; ii++)  {
      if (!cbuffer[ii]) {
	datamin = MIN(datamin, buffer[ii]);
	datamax = MAX(datamax, buffer[ii]);
      } else
	nulls++;
    }
#else
    nulls = 0;
    for (j=0; i<naxes[1]; j++) {
      for (i=0; i<naxes[0]; i++) {
	if (!cdata[j][i]) {
	  datamin = MIN(datamin, data[j][i]);
	  datamax = MAX(datamax, data[j][i]);
	} else
	  nulls++;
      }
    }
#endif

    printf("\nMin and max image pixels =  %g, %g  ; found %d/%ld null values\n",
	   datamin, datamax, nulls, npixels);

    fits_close_file(fptr, &status);
    fits_error(status,"close_file");

#if 1
    /* we actually don't need the data anymore, just the mask file */
    free(mask[slot].buffer);
    free(mask[slot].data);
#endif
    set_slot(slot);
    return;
}

/** @fn int is_maski(int ix, int iy, int slot)
 *
 * @retval is_maski      1 if the pixel is bad, 0 if it's good
 */

int is_maski(int ix, int iy, int slot)
{
  check_slot(slot);
  mask_image *mp = &mask[slot];
  double x,y;
  int status = 0;

  if (ix < 0 || ix >= mp->x.n) return 0;   /* outside the grid */
  if (iy < 0 || iy >= mp->y.n) return 0;
#if 1
  fits_pix_to_world( (double) ix, (double) iy, mp->x.crval, mp->y.crval,
		     mp->x.crpix, mp->y.crpix, mp->x.cdelt, mp->y.cdelt,
		     mp->rot, mp->proj, &x, &y, &status);
  fits_error(status,"pix_to_world");
  dprintf(1,"%d %d -> %g %g\n",ix,iy, x, y);
#endif

  return mp->cdata[iy][ix] == 0;           /* cdata=0 means good data */
}

/** @fn int is_mask(double ra_deg, double dec_deg, int slot)
 *
 *  @retval is_mask     returns 1 if the pixel is bad, 0 if it's good 
 */

int is_mask(double ra_deg, double dec_deg, int slot)
{ /* returns 1 if outside the grid (i.e. yes, masked)
     and a 0 if inside the fits file (no, not masked) */
  check_slot(slot);
  mask_image *mp = &mask[slot];
  double x, y;
  int ix, iy, status=0;

  fits_world_to_pix(ra_deg, dec_deg, mp->x.crval, mp->y.crval,
		    mp->x.crpix, mp->y.crpix, mp->x.cdelt, mp->y.cdelt,
		    mp->rot, mp->proj, &x, &y, &status);
  fits_error(status,"world_to_pix");
  dprintf(1,"masks.c:is_mask: %g %g -> %g %g\n",ra_deg, dec_deg, x, y);

  if (x < 0 || x >= mp->x.n) return 1;   /* outside the grid */
  if (y < 0 || y >= mp->y.n) return 1;
  ix = x;
  iy = y;
  return mp->cdata[iy][ix] == 1; /* cdata=0 means good data */
}

/** @fn int is_edge(double ra_deg, double dec_deg, int slot)
 *
 * @retval is_edge    1 if on the edge(i.e. yes, near edge)
 *          and a 0 not on an edge (no, not near an edge)
 */

int is_edge(double ra_deg, double dec_deg, int slot)
{
  const int edge_dist=4; /* num pixel distance to an edge */
  check_slot(slot);
  mask_image *mp = &mask[slot];
  double x, y;
  int ix, iy, status=0;
  int ixmin,ixmax,iymin,iymax;

  fits_world_to_pix(ra_deg, dec_deg, mp->x.crval, mp->y.crval,
		    mp->x.crpix, mp->y.crpix, mp->x.cdelt, mp->y.cdelt,
		    mp->rot, mp->proj, &x, &y, &status);
  fits_error(status,"world_to_pix");
  dprintf(1,"%g %g -> %g %g\n",ra_deg, dec_deg, x, y);
  ix = x;
  iy = y;
  ixmin=ix-edge_dist; /* ixmin */
  ixmax=ix+edge_dist; /* ixmax */
  iymin=iy-edge_dist; /* iymin */
  iymax=iy+edge_dist; /* iymax */

  if (ixmin < 0 || ixmax >= mp->x.n) return 1;   /* outside the grid */
  if (iymin < 0 || iymax >= mp->y.n) return 1;

  if (mp->cdata[iy][ixmin] == 1)
     return 1; /* pixel is near left edge */
  if (mp->cdata[iy][ixmax] == 1)
     return 1; /* pixel is near right edge */
  if (mp->cdata[iymin][ix] == 1)
     return 1; /* pixel is near bottom edge */
  if (mp->cdata[iymax][ix] == 1)
     return 1; /* pixel is near top edge */
  return 0; /* default return value.  Means not near an edge */
}

/** @fn void end_mask(int slot)
 *  end mask (dummy routine)
 */

void end_mask(int slot)
{
  warning("end_mask: this routine still needs to be implemented");
}



#ifdef TESTBED

string defv[] = {
  "in1=\n    fits file 1",
  "in2=\n    fits file 2",
  "in3=\n    fits file 3",
  "in4=\n    fits file 4",
  "ix=\n   pixel",
  "iy=\n   pixel",
  "x=\n     wcs",
  "y=\n     wcs",
  "VERSION=2\n  pjt",
  NULL,
};

string usage="test";

nemo_main()
{
  char key[10];
  int i, ix,iy;
  double x,y;

  for (i=1; i<=4; i++) {
    sprintf(key,"in%d",i);
    if (!hasvalue(key)) continue;
    ini_mask(getparam(key),i);
  }
  if (hasvalue("ix")) {
    ix = getiparam("ix");
    iy = getiparam("iy");
    printf("im_maski(%d,%d,%d) = %d\n",ix,iy,1,is_maski(ix,iy,1));
  }
  if (hasvalue("x")) {
    x = getdparam("x");
    y = getdparam("y");
    printf("im_mask(%g,%g,%d) = %d\n",x,y,1,is_mask(x,y,1));
  }


}

#endif

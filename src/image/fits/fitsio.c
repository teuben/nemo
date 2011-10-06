/************************************************************************/
/*                                                                      */
/*      This module consists of a number of routines to perform         */
/*      I/O on a FITS image file.                                       */
/*                                                                      */
/*      originally written by Bob Sault while at:                       */
/*        Giant Metrewave Telescope Project, Tata Institute             */
/*        Poona University, Pune, India.                                */
/*	adapted for NEMO by Peter Teuben                                */
/*                                                                      */
/*  This should be reasonably easy adapted to any system which supports */
/*  a standard C compiler and UNIX-like open,read,write and close       */
/*  system subroutines. However, note the following:                    */
/*  * Currently only machines which support FITS integers and IEEE      */
/*    floating point data                                               */
/*  * endianism is taken care of by WORDS_BIGENDIAN, normally set via   */
/*    autoconf                                                          */
/*  * A #define of BSD should be inserted on Berkeley UNIX machines.    */
/*  * A typedef FLOAT is used for the standard floating point data type.*/
/*    This can be either "float" or "double".                           */
/*                                                                      */
/*  History:                                                            */
/*    15-jul-90 Original version. rjs.                                  */
/*    28-jul-90 hurah - replace werong in NEMO. pjt                     */
/*    30-jul-90 Fixed bug in fitsetpl. Added comments. Improved some    */
/*              error messages. rjs.                                    */
/*     8-aug-90 Fixed bug in fitread, when reading 4-byte integers. rjs.*/
/*     1-oct-90 merged PJT and RJS versions - (PJT) - check fitsio.h    */
/*              kludged sprintf() when %-#s when strlen()>>#            */
/*              see fitwrhda and fitwra                                 */
/*    10-oct-90 made it reset fits->type according to bitpix at output  */
/*		and allow different blocksize from 2880			*/
/*    12-sep-91 added routine fitexhd() to return existence of keyword  */
/*    16-sep-91 dummy comments                                          */
/*    23-sep-91 some support for bitpix=-64; fixed fitexhd() bug        */
/*    23-dec-91 support for bitpix=8 images                             */
/*     7-mar-92 happy gcc2.0 - bad strcmp(,,) -> strncmp()              */
/*    21-may-92 fixed close bug when wrong bitpix..  pjt		*/
/*    26-may-93 fixed reading bitpix=-64 data                           */
/*    20-jan-94 solaris fixes						*/
/*    22-feb-94 ansi - nemo - stdio (fopen instead of open)             */
/*    10-jun-94 more ansifying - seems like fixing compiler bugs        */
/*    13-mar-95 added support for -DNEED_SWAP (linux, alpha etc.)       */
/*    20-mar-95 fixed bugs when card is exactly 80 chars long           */
/*     5-jan-99 added fitswrkv						*/
/*    18-may-99 added support to read CD_i_j matrices                   */
/*    15-oct-99 string.h into stdinc.h					*/
/*    21-mar-00 fixed offset bug for raw mode - using WORDS_BIGENDIAN   */
/*     9-jul-00 added message in first fitopen reporting endianism      */
/*              fixed offset bug when fitsetpl called before fitwrite   */
/*     7-aug-01 fixed bug rounding bug getting offset in reading        */
/*              (redefined f->ncards now to be 1-based upon reading too)*/
/*    28-sep-01 added bitpix=64 support - portability not solved yet    */
/*              seems bitpix=-64 is working ok                          */
/*    12-oct-01 new standard added FITS reference                       */
/*    23-jul-02 add fitresize                                           */
/*    14-nov-02 add dummy DATASUM and CHECKSUM                          */
/*     8-nov-05 also recognize XTENSION = 'IMAGE'                       */
/*    11-dec-06 store cvsID in output                                   */
/* ToDo:                                                                */
/*  - BLANK substitution                                                */
/************************************************************************/

#include <stdinc.h>
#include <ctype.h>
#include <fitsio_nemo.h>

static char *cvs_id="$Id$ fitsio.c";

/* 
 * this next #def's obviously needs to be refined. 
 * This file can currently only deal with simple endian swappings
 * There is no support for non-IEEE or non-twos complement
 * we used to use NEED_SWAP, now we use WORDS_BIGENDIAN (from $NEMO's
 * autoconf generated include files)
 */

#if 0
#define DEBUG 1
#endif

#if SIZEOF_LONG_LONG==8
typedef long long int8;         /* e.g. i386; sparc <= sol7; ppc ? */
#elif SIZEOF_LONG==8
typedef long int8;              /* e.g. alpha 64's */
#elif SIZEOF_INT==8
typedef int int8;               /* will never happen ? */
#else
typedef long int8;		/* some stupid fallback, probably wrong */
#endif

local int  fitsrch    (FITS *, char *, char *);
local void fitpad     (FITS *, int, char),
           fitput     (FITS *, char *),
           fitalloc   (FITS *),
           fitcvt_8i  (char *, byte *, int),
           fitcvt_16i (char *, short int *, int),
           fitcvt_32i (char *, int *, int),
           fitcvt_64i (char *, int8 *, int),
           fitcvt_fr  (char *, FLOAT *, int),
           fitcvt_dr  (char *, DOUBLE *, int),
           fitcvt_i8  (byte *, char *, int),
           fitcvt_i16 (short int *, char *, int),
           fitcvt_i32 (int *, char *, int),
           fitcvt_i64 (int8 *, char *, int),
           fitcvt_rf  (FLOAT *, char *, int),
           fitcvt_rd  (DOUBLE *, char *, int);

local char *buf1=NULL, *buf2=NULL;	/* local conversion buffers */
local int maxdim=0;
local int w_bitpix = -32;               /* see: fit_setbitpix()    */
local FLOAT w_bscale = 1.0;             /* see: fit_setscale()     */
local FLOAT w_bzero = 0.0;              /* see: fit_setscale()     */
local int blocksize= 2880;	        /* See: fit_setblocksize() */
local int first_message = 1;		/* See: fitopen */

local string cfits1="FITS (Flexible Image Transport System) format is defined in 'Astronomy";
local string cfits2="and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H";
local string cfits3="nemo::fitsio.c $Id$";

/**********************************************************************/
FITS *fitopen(string name,string status,int naxis,int *nsize)
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
  FITS *f;
  int n,t,i,size,bitpix;
  char keyword[9],line[81];

  if (first_message) {
#ifdef WORDS_BIGENDIAN
    dprintf(1,"fitopen: Big-endian machine; no need to swap bytes\n");
#else
    dprintf(1,"fitopen: Little-endian machine; swapping bytes for FITSIO\n");
#endif
    first_message = 0;
  }
#if 0  
  if (sizeof(int8) != 8)
    warning("Cannot read BITPIX=64 items with this code, sizeof(int8)=%d",sizeof(int8));
  else
#endif
    dprintf(1,"This version has experimental handling of BITPIX 64 data\n");

  f = (FITS *)allocate(sizeof(FITS));

/* Handle a new file. */

  if(streq(status,"new") || streq(status,"w")){
    f->status = STATUS_NEW;
    if(naxis <= 0 || naxis > MAXNAX)
      error("Bad value of NAXIS, when opening %s",name);
    for(i=0; i < naxis; i++){
      if(nsize[i] <= 0)
        error("Bad image dimensions for file %s",name);
    }
    f->fd = stropen(name,"w");
    switch (w_bitpix) {
      case   8:    f->type = TYPE_8INT;   break;
      case  16:    f->type = TYPE_16INT;  break;
      case  32:    f->type = TYPE_32INT;  break;
      case  64:    f->type = TYPE_64INT;  break;  /* experimental */
      case -32:    f->type = TYPE_FLOAT;  break;
      case -64:    f->type = TYPE_DOUBLE; break;
      default:
          error("Illegal bitpix %d",w_bitpix);
    }
    f->bytepix = ABS(w_bitpix)/8;

    f->ncards = 0;
    fitwrhdl(f,"SIMPLE",TRUE);
    fitwrhdi(f,"BITPIX",w_bitpix);
    fitwrhdi(f,"NAXIS",naxis);
    for(i=0; i < naxis; i++){
      sprintf(keyword,"NAXIS%d",i+1);
      fitwrhdi(f,keyword,nsize[i]);
    }
    if (w_bitpix>0 || w_bscale!=1 || w_bzero!=0) {
        fitwrhdr(f,"BSCALE",w_bscale);
        fitwrhdr(f,"BZERO",w_bzero);
    }
    f->bscale = w_bscale;
    f->bzero  = w_bzero;
    fitwra(f,"COMMENT",cfits1);
    fitwra(f,"COMMENT",cfits2);
    fitwra(f,"COMMENT",cfits3);
    fitwrhda(f,"DATASUM", "0000000000000000");
    fitwrhda(f,"CHECKSUM","0000000000000000");

/* Handle an old file. */

  } else if(streq(status,"old") || streq(status,"r")) {
    f->status = STATUS_OLD;
    f->fd = stropen(name,"r");

/* Check it has SIMPLE and END keywords. Calculate the byte offset to
   the start of the image data. -- ncards is now 1 based ! */

    if(fitsrch(f,"SIMPLE  ",line) == 1) {
      dprintf(1,"Detected SIMPLE=\n");
    } else if (fitsrch(f,"XTENSION",line) == 1) {
      dprintf(1,"Detected XTENSION=\n");
    } else
      error("FITS file with no SIMPLE= or XTENSION=");
    f->ncards = fitsrch(f,"END     ",line);
    if (f->ncards < 0)
      error("no END found: File \'%s\' does not appear to be FITS",name);
    f->offset = blocksize*((80*f->ncards + (blocksize-1))/blocksize); /* round up */
    f->skip = f->offset;		/* offset can change !!! */
    dprintf(1,"END found at card %d, offset=%d\n",f->ncards,f->offset);

/* Determine things about the file. */

    fitrdhdi(f,"BITPIX",&bitpix,0);
    if(bitpix ==   8)      f->type = TYPE_8INT;
    else if(bitpix == 16)  f->type = TYPE_16INT;
    else if(bitpix ==  32) f->type = TYPE_32INT;
    else if(bitpix ==  64) f->type = TYPE_64INT;  /* experimental */
    else if(bitpix == -32) f->type = TYPE_FLOAT;
    else if(bitpix == -64) f->type = TYPE_DOUBLE;
    else
      error("Unsupported value of BITPIX (%d) when opening %s",bitpix,name);
    f->bytepix = ABS(bitpix)/8;
    fitrdhdi(f,"NAXIS",&n,0);
    if(n <= 0)
      error("Bad value of NAXIS (%d) when opening %s",n,name);
    size = 1;
    for(i=0; i < n; i++){
      sprintf(keyword,"NAXIS%d",i+1);
      fitrdhdi(f,keyword,&t,0);
      if(t <= 0 || (i >= naxis && t != 1))
        error("Cannot handle dimension %d of %s being %d\n",i+1,name,t);
      if(i < naxis)nsize[i] = t;
      size *= t;
    }
    for(i=n; i < naxis; i++) nsize[i] = 1;
    fitrdhdr(f,"BSCALE",&(f->bscale),1.0);
    fitrdhdr(f,"BZERO",&(f->bzero),0.0);

/* Check that the file is of the right size. */

    size *= f->bytepix;
    size += f->offset;
    if(fseek(f->fd,0,2))
      error("Cannot seek - file %s appears too small",name);
/* Neither old nor new - an error for now */

  } else
    error("Bad status argument \'%s\', in fitopen",status);

/* Save dimension info. */

  f->naxis = naxis;
  for(i=0; i < naxis; i++)f->axes[i] = nsize[i];
  for(i=naxis; i < MAXNAX; i++) f->axes[i] = 1;

  return f; /* Return an opaque (although typed) pointer */
}

/* Warning: should not use this for the miriad based version......
 */

#if 0
void fitresize(FITS *file,int naxis,int *nsize)
{
  int i;
  char keyword[9];

  for(i=0; i < naxis; i++) {
      file->axes[i] = nsize[i];
      sprintf(keyword,"NAXIS%d",i+1);
      fitwrhdi(file,keyword,nsize[i]);
  }
}
#endif
/**********************************************************************/
void fitclose(FITS *file)
/*
  This closes a FITS file, and deletes any memory associated with it.

  Input:
    file        This is the pointer returned by fitopen.
----------------------------------------------------------------------*/
{
  FITS *f;
  int i,offset;

  f = file;
  if(f->status == STATUS_NEW){
    warning("No image data written to FITS file");
    fitpad(f,80*f->ncards,' ');
  } else if(f->status == STATUS_NEW_WRITE){
    offset = f->bytepix;
    for(i=0; i < f->naxis; i++)offset *= f->axes[i];
    offset += blocksize*((80*f->ncards + (blocksize-1))/blocksize);
    fitpad(f,offset,0);
  }
  strclose(f->fd);
  free((char *)f);
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

  Data is read into 'buf1', then converted into 'buf2', (which also swaps
  the bytes in buf1 as a side effect), from where it
  is converted into a FLOAT array that was passed to this routine
----------------------------------------------------------------------*/
{
  int i, offset,length,bytes;
  FITS *f;
  FLOAT bscale,bzero;
  short int *idat;
  int  *jdat;
  int8 *kdat;
  byte *sdat;
  double *ddat;

  f = file;
  fitalloc(f);
  bytes = f->bytepix;
  offset = bytes*j*f->axes[0] + f->offset;
  length = bytes * f->axes[0];
  fseek(f->fd,offset,0);
  		dprintf(2,"fitread: offset(%d)=%d => %d\n",j,f->offset,offset);
  if(j > f->axes[1])
    error("Attempt to read beyond image boundaries, in fitread");
  else if(f->status != STATUS_OLD)
    error("Attempt to read from a new file, in fitread");
  else if(length != fread(buf1,1,length,f->fd))
    error("I/O read error in fitread");

/* We have the (possibly swapped) raw data now. Convert and scale it. */

  bscale = f->bscale; bzero = f->bzero;
  if(f->type == TYPE_8INT){
    sdat = (byte *)buf2;
    fitcvt_8i(buf1,sdat,f->axes[0]);
    for(i=0; i < f->axes[0]; i++) *data++ = bscale * *sdat++ + bzero;/* ansi warning */
  } else if(f->type == TYPE_16INT){
    idat = (short int *)buf2;
    fitcvt_16i(buf1,idat,f->axes[0]);
    for(i=0; i < f->axes[0]; i++) *data++ = bscale * *idat++ + bzero;
  } else if(f->type == TYPE_32INT){
    jdat = (int *)buf2;
    fitcvt_32i(buf1,jdat,f->axes[0]);
    for(i=0; i < f->axes[0]; i++) *data++ = bscale * *jdat++ + bzero;
  } else if(f->type == TYPE_64INT){
#ifdef DEBUG
    kdat = (int8 *)buf1;
    /* print out raw data, will be in the wrong endian on e.g. i386 */
    /* also note %lld format may not be supported on all versions of printf.. */
    for(i=0; i < f->axes[0]; i++) dprintf(3,"DEBUG1: %d -> %lld\n",i,kdat[i]);
#endif
    kdat = (int8 *)buf2;
    fitcvt_64i(buf1,kdat,f->axes[0]);
#ifdef DEBUG
    kdat = (int8 *)buf2;
    for(i=0; i < f->axes[0]; i++) dprintf(3,"DEBUG2: %d -> %lld\n",i,kdat[i]);
#endif
    for(i=0; i < f->axes[0]; i++) *data++ = bscale * *kdat++ + bzero;
  } else if(f->type == TYPE_DOUBLE){
    ddat = (double *)buf1;
    fitcvt_dr(buf1,ddat,f->axes[0]);
    for(i=0; i < f->axes[0]; i++) *data++ = bscale * *ddat++ + bzero;
  } else if(f->type == TYPE_FLOAT) {
    fitcvt_fr(buf1,data,f->axes[0]);
    if(bscale != 1 || bzero != 0) {
      for(i=0; i < f->axes[0]; i++) {
        *data = bscale * *data + bzero;
        data++;
      }
    }
  } else {
    error("fitread: Illegal datatype %d",f->type);    /* NEVER REACHED */
  }
}
/**********************************************************************/
void fitwrite(FITS *file, int j, FLOAT *data)
/*
  This writes a row of a FITS image. Note that this should not be called
  until the programmer is done writing to the FITS header (routines fitwrhd).
                    
  Inputs:
   file        The pointer returned by fitopen.
   j           The row number to be written. This varies from 0 to naxis2-1.
   data        A FLOAT array of naxis1 elements, being the pixel values  
                to write.

   FLOAT array data is written into buf1,  then converted to buf2 as FITS
   format, which is then written to the fits file
----------------------------------------------------------------------*/
{
  int i,offset,length,bytes;
  FITS *f;
  FLOAT bscale,bzero, *fdat;
  int *jdat;
  short int *idat;
  int8 *kdat;
  byte *sdat;
  double *ddat;

  f = file;
  fitalloc(f);
/* Finish off the header, if this is the first write. */

  if(f->status == STATUS_NEW){
    fitput(f,"END");
    fitpad(f,80*f->ncards,' ');
    f->skip = f->offset = blocksize*((80*f->ncards + (blocksize-1))/blocksize);
    f->status = STATUS_NEW_WRITE;
  } else if(f->status != STATUS_NEW_WRITE) {
    error("Illegal operation, in fitwrite");
  } if(j >= f->axes[1]){
    error("Attempt to write beyond image boundaries, in fitwrite");
  }

  bytes = f->bytepix;
  offset = bytes*j*f->axes[0] + f->offset;
  length = bytes * f->axes[0];
  fseek(f->fd,offset,0);

/* Convert the data. */

  bscale = f->bscale; bzero = f->bzero;
  if(f->type == TYPE_8INT) {
    sdat = (byte *)buf2;
    for(i=0; i < f->axes[0]; i++)
        *sdat++ = (byte) ( (*data++ - bzero) / bscale);
    sdat = (byte *)buf2;	/* reset pointer to buffer */
    fitcvt_i8(sdat,buf1,f->axes[0]);
  } else if(f->type == TYPE_16INT){
    idat = (short int *)buf2;   /* set pointer to conversion buffer */
    for(i=0; i < f->axes[0]; i++)
        *idat++ = (short int) ( (*data++ - bzero) / bscale);
    idat = (short int *)buf2;	/* reset pointer to buffer */
    fitcvt_i16(idat,buf1,f->axes[0]);
  } else if(f->type == TYPE_32INT){
    jdat = (int *)buf2;         /* set pointer to conversion buffer */
    for(i=0; i < f->axes[0]; i++)
        *jdat++ = (int) ( (*data++ - bzero) / bscale);
    jdat = (int *)buf2;		/* reset pointer to buffer */
    fitcvt_i32(jdat,buf1,f->axes[0]);
  } else if(f->type == TYPE_64INT){
    kdat = (int8 *)buf2;         /* set pointer to conversion buffer */
    for(i=0; i < f->axes[0]; i++)
        *kdat++ = (int8) ( (*data++ - bzero) / bscale);
    kdat = (int8 *)buf2;		/* reset pointer to buffer */
    fitcvt_i64(kdat,buf1,f->axes[0]);
  } else if(f->type == TYPE_DOUBLE){
    ddat = (double *)buf2;	/* set pointer to conversion buffer */
    for(i=0; i < f->axes[0]; i++)
        *ddat++ = (double)( (*data++ - bzero) / bscale);
    ddat = (double *)buf2;	/* reset pointer to buffer */
    fitcvt_rd(ddat,buf1,f->axes[0]);
  } else if(f->type == TYPE_FLOAT){    
    if(bscale != 1 || bzero != 0) {     /* only convert when needed */
      fdat = data;                      /* BUT DATA IS MODIFIED !!! */
      for(i=0; i < f->axes[0]; i++) {
        *fdat = ( (*fdat - bzero) / bscale);
        fdat++;
      }
    }
    fitcvt_rf(data,buf1,f->axes[0]);
  } else {
    error("fitwrite: Illegal datatype %d",f->type);    /* NEVER REACHED */
  }

/* More checking, then do the output. */

  if(length != fwrite(buf1,1,length,f->fd))
    error("I/O write error in fitwrite");
}
/**********************************************************************/
void fitsetpl(FITS *file, int n, int *nsize)
/*
  This sets the plane to be accessed in a FITS file which has more than
  two dimensions.

  Input:
    file        The pointer returned by fitopen.
    n           This gives the size of the nsize array.
    nsize       This gives the indices of the higher dimensions of the
                FITS image. They are zero-relative. nsize[0] gives the
                index along the 3rd dimension, nsize[1] is the indice along
                the 4th dimension, etc.
----------------------------------------------------------------------*/
{
  FITS *f;
  int i;
  size_t offset;

  f = file;
  if(f->status == STATUS_NEW){
    fitput(f,"END");
    fitpad(f,80*f->ncards,' ');
    f->skip = blocksize*((80*f->ncards + (blocksize-1))/blocksize);    
    f->status = STATUS_NEW_WRITE;
  }
  offset = 0;
  for(i=n-1; i >= 0; i--){
    if(nsize[i] >= f->axes[i+2]){
      error("Illegal coordinate index, in fitsetpl, aborting ...");
    }
    offset = offset * f->axes[i+2] + nsize[i];
  }
  offset *= f->bytepix * f->axes[0] * f->axes[1];
  offset += f->skip;
  f->offset = offset;
  if (f->skip == 0)
  	warning("fitsetpl: f->skip is 0, should be multiple of 2880");
  if (offset < 0)
	error("fitsetpl: bad offset=%ld (%d,...)\n",offset,nsize[0]);
  dprintf(4,"fitsetpl: offset=%ld (%d,...)\n",offset,nsize[0]);
}
/**********************************************************************/
void fit_setbitpix(int bp)
/*
    fit_setbitpix:  set the bitpix value for subsequent output
    Note: this routine must be called before fitopen

    Input:  
        bp      bitpix value; supported are: 8, 16, 32, -32, -64
----------------------------------------------------------------------*/
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
----------------------------------------------------------------------*/
{
    dprintf(1,"BSCALE and BZERO reset to %g and %g\n",bs,bz);
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
    Note: this routine must be called before fitopen
    
    Input:
        n       blocksize to be set
----------------------------------------------------------------------*/
{
    dprintf(1,"BLOCKSIZE reset to %d\n",n);
    blocksize = n;
}
/**********************************************************************/
void fitrdhdr(FITS *file, string keyword, FLOAT *value, FLOAT def)
/*
  This reads the value of a real-valued FITS keyword from the file header.

  Input:
    file        The pointer returned by fitopen.
    keyword     The keyword to search for.
    def         If the keyword is not found, this "default value" is
                returned.
  Output:
    value       The value read from the FITS header.
----------------------------------------------------------------------*/
{
  char card[81],*s,*s1;
  double d;

  if(fitsrch(file,keyword,card) < 0) *value = def;
  else {
    card[80] = 0;
    s = card + strlen(keyword);
    while(*s == ' ' && *s != 0)s++;
    while(*s == '=' && *s != 0)s++;
    while(*s == ' ' && *s != 0)s++;
    s1 = s;
    while(isdigit(*s) || *s == '+' || *s == '-' || *s == '.' || *s == 'E' ||
      *s == 'D' || *s == 'e' || *s == 'd')s++;
    *s = 0;
    if(sscanf(s1,"%lg",&d) != 1) *value = def;
    else *value = d;
  }
}
/**********************************************************************/
void fitrdhdi(FITS *file, string keyword, int *value, int def)
/*
  This reads the value of a integer-valued FITS keyword from the file header.
  It's done rather liberally, but reading it has real-valued, and converting
  back to an integer.

  Input:
    file        The pointer returned by fitopen.
    keyword     The keyword to search for.
    def         If the keyword is not found, this "default value" is
                returned.
  Output:
    value       The value read from the FITS header.
----------------------------------------------------------------------*/
{
  FLOAT temp;
  fitrdhdr(file,keyword,&temp,(FLOAT)def);
  *value = (int) temp;
}
/***********************************************************************/
void fitrdhda(FITS *file, string keyword, char *value, char *def)
/*
  This reads the value of a character-valued FITS keyword from the file header.

  Input:
    file        The pointer returned by fitopen.
    keyword     The keyword to search for.
    def         If the keyword is not found, this "default value" is
                returned.
  Output:
    value       The value read from the FITS header.
----------------------------------------------------------------------*/
{
  char card[81],*s,*s1;

  if(fitsrch(file,keyword,card) < 0) strcpy(value,def);
  else {
    card[80] = 0;
    s = card + strlen(keyword);
    while(*s == ' ' && *s != 0)s++;
    while(*s == '=' && *s != 0)s++;
    while(*s == ' ' && *s != 0)s++;
    if(*s == '\'') s++;
    s1 = s;
    while(*s != '\'') s++;
    *s = 0;
    strcpy(value,s1);
  }
}
/**********************************************************************/
void fitwrhdr(FITS *file, string keyword, FLOAT value)
/*
  This writes the value of a real-valued FITS keyword to the file header.

  Input:
    file        The pointer returned by the fitopen routine.
    keyword     The name of the keyword to write.
    value       The value of the keyword.
----------------------------------------------------------------------*/
{
  char line[81];
  sprintf(line,"%-8s= %20.9E /",keyword,value);
  fitput(file,line);
}
/**********************************************************************/
void fitwrhdi(FITS *file, string keyword, int value)
/*
  This writes the value of a integer-valued FITS keyword to the file header.

  Input:
    file        The pointer returned by the fitopen routine.
    keyword     The name of the keyword to write.
    value       The value of the keyword.
----------------------------------------------------------------------*/
{
  char line[81];
  sprintf(line,"%-8s= %20d /",keyword,value);
  fitput(file,line);
}
/**********************************************************************/
void fitwrhdl(FITS *file, string keyword, int value)
/*
  This writes the value of a logical-valued FITS keyword to the file header.

  Input:
    file        The pointer returned by the fitopen routine.
    keyword     The name of the keyword to write.
    value       The value of the keyword.
----------------------------------------------------------------------*/
{
  char line[81];
  sprintf(line,"%-8s=                    %c /",keyword,(value ? 'T' : 'F'));
  fitput(file,line);
}
/**********************************************************************/
void fitwrhda(FITS *file, string keyword, string value)
/*
  This writes the value of a character-valued FITS keyword to the file header.
----------------------------------------------------------------------*/
{
  char line[81], tmp[81];
  if ((int)strlen(value)<8)
    sprintf(line,"%-8s= '%-8s'           /",keyword,value);  /* 32 */
  else if ((int)strlen(value)>68) {
    strncpy(tmp,value,80);
    sprintf(line,"%-8s= '%-68s'",keyword,tmp);
  } else
    sprintf(line,"%-8s= '%s'",keyword,value);
  fitput(file,line);
}
/**********************************************************************/
void fitwra(FITS *file, string keyword, string value)
/*
  This writes a character-value to a FITS keyword without '='
  as e.g. used by COMMENT and HISTORY.
----------------------------------------------------------------------*/
{
  char line[81], tmp[81];
  if ((int)strlen(value)>71) {
    strncpy(tmp,value,80);
    sprintf(line,"%-8s %-71s",keyword,tmp);
  } else
    sprintf(line,"%-8s %s",keyword,value);
  fitput(file,line);
}
/**********************************************************************/
void fitwrhd(FITS *file, string keyword, string value)
/*
  This writes a character-value to a FITS keyword with '='
----------------------------------------------------------------------*/
{
  char line[81], tmp[81];
  if ((int)strlen(value)>70) {
    strncpy(tmp,value,80);
    sprintf(line,"%-8s= %-70s",keyword,tmp);
  } else
    sprintf(line,"%-8s= %s",keyword,value);
  fitput(file,line);
}


/*====================================================================*/
/**********************************************************************/
local void fitpad(FITS *f, int offset, char pad)
/*
  This pads a FITS file from location 'offset' (0-based) up to the 
  next 2880 block boundary.
----------------------------------------------------------------------*/
{
#define MAXLEN 512
  char buf[MAXLEN];
  int k,ktot,i,length;
  for(i=0; i < MAXLEN; i++) buf[i] = pad;
  k = offset;
  ktot = blocksize*((k + (blocksize-1))/blocksize);
  fseek(f->fd,k,0);
  while(k < ktot){
    length = ktot - k;
    if(length > MAXLEN) length = MAXLEN;
    if(length != fwrite(buf,1,length,f->fd)){
      error("I/O error while padding FITS file, aborting ...");
    }
    k += length;
  }
}
/**********************************************************************/
local void fitput(FITS *f,string card)
/*
  This adds any FITS card to the header; it therefore also keeps track 
  of the number of cards written to the header (ncards)
----------------------------------------------------------------------*/
{
  int i;
  char line[81],*s;

  if(f->status != STATUS_NEW){
    error("Illegal call to a fitwrhd routine, aborting ...");
  }
  s = line;
  while(*card != 0)*s++ = *card++;
  for(i = s - line; i < 80; i++) *s++ = ' ';
  fseek(f->fd,80*(f->ncards++),0);
  if(80 != fwrite(line,1,80,f->fd)){
    error("Error writing a FITS header card, aborting ...");
  }
}
/**********************************************************************/
int fitexhd(FITS *f, string keyword)
/*
  Public routine to check for a FITS keyword existence in a file.
  Returns 1 if found, 0 if not.
  Does not work on NEW files !!!
----------------------------------------------------------------------*/
{
  char card[81];
  if (fitsrch(f,keyword,card)<0)
    return(0);                      /* card not found */
  else
    return(1);                      /* card found */
}
/**********************************************************************/
local int fitsrch(FITS *f,string keyword,char *card)
/*
  This searches for a FITS keyword in a file.
  Returns cardnumber index (1 based) :
  -1 if not found, and 1 or higher if found
----------------------------------------------------------------------*/
{
  int length,ncard;

  length = strlen(keyword);
  ncard = 0;
  fseek(f->fd,0,0);
  while(fread(card,1,80,f->fd) == 80){
    if((card[length] == ' ' || card[length] == '=') &&
       !strncmp(card,keyword,length))return(ncard+1);
    else if(!strncmp(card,"END     ",8)) return(-1);
    ncard++;
  }
  return(-1);
}
/**********************************************************************/
local void fitcvt_8i(char *in,byte *out,int n)
/*
  This converts an array of FITS 8 bit integers into byte

  Input:
    in          An array of 8 bit FITS integers.
    n           Number of integers to convert.
  Output:
    out         The array of host bytes.
----------------------------------------------------------------------*/
/* The following assumes that 8 bit FITS integers are equivalent to
   host bytes */
{
  memcpy((char *)out,in,sizeof(byte)*n);
}
/**********************************************************************/
local void fitcvt_16i(char *in,short int *out, int n)
/*
  This converts an array of FITS 16 bit integers into host short integers.

  Input:
    in          An array of 16 bit FITS integers.
    n           Number of integers to convert.
  Output:
    out         The array of host short integers.
----------------------------------------------------------------------*/
/* The following assumes that 16 bit FITS integers are equivalent to
   host short integers. */
{
#ifndef WORDS_BIGENDIAN 
  bswap(in, 2, n);
#endif
  memcpy((char *)out,in,sizeof(short int)*n);
}
/**********************************************************************/
local void fitcvt_32i(char *in, int *out, int n)
/*
  This converts an array of FITS 32 bit integers into host integers.

  Input:
    in          An array of 32 bit FITS integers.
    n           Number of integers to convert.
  Output:
    out         The array of host integers.
----------------------------------------------------------------------*/
/* The following assumes that 32 bit FITS integers are equivalent to
   host integers. */
{
#ifndef WORDS_BIGENDIAN 
  bswap(in, 4, n);
#endif  
  memcpy((char *)out,in,sizeof(int)*n);
}
/**********************************************************************/
local void fitcvt_64i(char *in, int8 *out, int n)
/*
  This converts an array of FITS 64 bit integers into host int8's

  Input:
    in          An array of 64 bit FITS integers.
    n           Number of integers to convert.
  Output:
    out         The array of host int8's
----------------------------------------------------------------------*/
/* The following assumes that 64 bit FITS integers are equivalent to
   host int8's */
{
#ifndef WORDS_BIGENDIAN 
  bswap(in, 8, n);
#endif  
  memcpy((char *)out,in,sizeof(int8)*n);
}
/**********************************************************************/
local void fitcvt_fr(char *in,FLOAT *out, int n)
/*
  This converts an array of IEEE 32 bit floating point numbers into
  the hosts "FLOAT" format.

  Input:
    in          An array of 32 bit FITS integers.
    n           Number of integers to convert.
  Output:
    out         The array of host integers.
----------------------------------------------------------------------*/
/* The following assumes that the host "float" format is the IEEE
   32 bit format. It does not assume that a "FLOAT" is a "float". */
{
  int i;
  float *f;

#ifndef WORDS_BIGENDIAN 
  bswap(in, 4, n);
#endif  
  f = (float *)in;
  for(i=0; i < n; i++) *out++ = *f++;
}
/**********************************************************************/
local void fitcvt_dr(char *in, DOUBLE *out, int n)
/*
  This converts an array of IEEE 64 bit floating point numbers into
  the hosts "DOUBLE" format.

  Input:
    in          An array of 32 bit FITS integers.
    n           Number of integers to convert.
  Output:
    out         The array of host integers.
----------------------------------------------------------------------*/
/* The following assumes that the host "double" format is the IEEE
   64 bit format. It does not assume that a "DOUBLE" is a "double". */
{
  int i;
  double *f;

#ifndef WORDS_BIGENDIAN 
  bswap(in, 8, n);
#endif  
  f = (double *)in;
  for(i=0; i < n; i++) *out++ = *f++;
}
/**********************************************************************/
local void fitcvt_i8(byte *in,char *out, int n)
/*
  This converts an array of host byte integers into FITS 8 bit integers

  Input:
    in          The array of host bytes
    n           Number of integers to convert.
  Output:
    out         An array of 8 bit FITS integers.
----------------------------------------------------------------------*/
/* The following assumes that 8 bit FITS integers are equivalent to
   host bytes */
{
  memcpy(out,(char *)in,sizeof(byte)*n);
}
/**********************************************************************/
local void fitcvt_i16(short int *in,char *out, int n)
/*
  This converts an array of host short integers into FITS 16 bit integers

  Input:
    in          The array of host short integers.
    n           Number of integers to convert.
  Output:
    out         An array of 16 bit FITS integers.
----------------------------------------------------------------------*/
/* The following assumes that 16 bit FITS integers are equivalent to
   host short integers. */
{
  memcpy(out,(char *)in,sizeof(short int)*n);
#ifndef WORDS_BIGENDIAN 
  bswap(out, 2, n);
#endif  

}
/**********************************************************************/
local void fitcvt_i32(int *in,char *out,int n)
/*
  This converts an array of host integers into FITS 32 bit integers.

  Input:
    in          The array of host integers.
    n           Number of integers to convert.
  Output:
    out         An array of 32 bit FITS integers.
----------------------------------------------------------------------*/
/* The following assumes that 32 bit FITS integers are equivalent to
   host integers. */
{
  memcpy(out,(char *)in,sizeof(int)*n);
#ifndef WORDS_BIGENDIAN 
  bswap(out, 4, n);
#endif  

}
/**********************************************************************/
local void fitcvt_i64(int8 *in,char *out,int n)
/*
  This converts an array of host integers into FITS 64 bit integers.

  Input:
    in          The array of host integers.
    n           Number of integers to convert.
  Output:
    out         An array of 64 bit FITS integers.
----------------------------------------------------------------------*/
/* The following assumes that 64 bit FITS integers are equivalent to
   host integers. */
{
  memcpy(out,(char *)in,sizeof(int8)*n);
#ifndef WORDS_BIGENDIAN 
  bswap(out, 8, n);
#endif  

}
/**********************************************************************/
local void fitcvt_rf(FLOAT *in,char *out,int n)
/*
  This converts an array of host FLOATs into IEEE 32 bit floating
  point numbers.

  Input:
    in          An array of host FLOATs.
    n           Number of integers to convert.
  Output:
    out         The array of IEEE 32 bit floating point numbers.
----------------------------------------------------------------------*/
/* The following assumes that the host "float" format is the IEEE
   32 bit format. It does not assume that a "FLOAT" is a "float". */
{
  int i;
  float *f;

  f = (float *)out;
  for(i=0; i < n; i++) *f++ = *in++;
#ifndef WORDS_BIGENDIAN 
  bswap(out, 4, n);
#endif  

}
/**********************************************************************/
local void fitcvt_rd(DOUBLE *in,char *out,int n)
/*
  This converts an array of host DOUBLEs into IEEE 64 bit floating
  point numbers.

  Input:
    in          An array of host FLOATs.
    n           Number of integers to convert.
  Output:
    out         The array of IEEE 64 bit floating point numbers.
----------------------------------------------------------------------*/
/* The following assumes that the host "double" format is the IEEE
   64 bit format. It does not assume that a "DOUBLE" is a "double". */
{
  int i;
  double *f;

  f = (double *)out;
  for(i=0; i < n; i++) *f++ = *in++;
#ifndef WORDS_BIGENDIAN 
  bswap(out, 8, n);
#endif  

}

local void fitalloc(FITS *f)
{

/* Make sure we have enough memory to deal with this file. */
/* buf1 and buf2 are static buffers used for conversion of */
/* from one type to the other                              */
/* Reading goes into buf1, conversion into ..              */

  if(maxdim < f->axes[0]){
    maxdim = f->axes[0];
    buf1 = (buf1 == NULL ? (char *) allocate(sizeof(double)*maxdim) : 
			  (char *) reallocate(buf1,sizeof(double)*maxdim));
    buf2 = (buf2 == NULL ? (char *) allocate(8*maxdim) : 
			  (char *) reallocate(buf2,8*maxdim));
  }
}

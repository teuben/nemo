# include "fits.h"
# include <stdio.h>
extern char *sprintf();
# include <math.h>
# include <sys/file.h>
# include <ctype.h>
# include <strings.h>
#include <varargs.h>

#define FITSRECL 2880

struct of_imagehdr *imageHdr[MAXIMAGES];

static char EMPTY[] = "";
static char histKey[] = "HISTORY COMB";
static double sixty = 60.;
static double one = 1.0;
static float FDRAGON = 1.697655329e37;
static struct of_fitskey fitsKey[] = {

#include "fitskeywords.h"

};
#define FITSKEYLEN ( sizeof(fitsKey) / sizeof(struct of_fitskey) )

# ifndef PI
# define PI 3.14159265358979323846
# endif

/* Array to hold units conversion factors */
double cunits[] = { PI/180.,1.,PI/12.,PI/10800.};

#if 0
static char *fitsBuf;		/* pointer to buffer for fits data */
static char *headPos;		/* pointer to next card position */
static int ffd;			/* file descriptor for fits file or tape */
#endif 0

extern char *malloc();
extern void free();
long lseek();
int round();


/* fits.c*/

double delta(/*cdelt*/);
int imagelinked(/*imnumber*/);
int readhdr(/*ihp*/);
struct of_fitskey *findkey(/*cardp*/);
void sethdr(/*fk, ih, cp*/);
struct AXISNAMES FitsAxisUnits(/*ctype*/);
char *ChkAlloc(/*size,name*/);

/*
 * Open a FITS file and read its header into an imageHdr.  Take care of any
 * existing image at this imnumber.  If asked to reopen the same file which
 * is currently open at the same imnumber, no action is taken.
 */
void OpenFits(imnumber,fname)
register int imnumber;
char *fname;
{
	register struct of_imagehdr *ihp;

	ImNumErr(imnumber);
	if(access(fname, 4) < 0)
		error("Can't access %s",fname);
	NewImage(imnumber);
	ihp = imageHdr[imnumber];
	ihp->ftype = EXTERNAL | FITS;
	ihp->cdelt1 = 1;	/* Defaults for simple headers */
	ihp->cdelt2 = -1;	/* Don't default to LH coordinate system */
	ihp->blank = FDRAGON;
	strcpy(ihp->fname, fname);
	if( (ihp->fd = open(fname, O_RDONLY, 0664) ) < 0) {
		CloseImage(imnumber);
		fprintf(stderr, "Error trying to open %s because", fname);
		perror(EMPTY);
		error(EMPTY);
	}
	if(readhdr(ihp)) {
		CloseImage(imnumber);
		error("Error reading header of %s", fname);
	}
	if(ihp->bscale == 0) {
		static float dmax[3] = {255,  32767,  2147483647};
		static float dmin[3] = {  0, -32767, -2147483647};
		ihp->bscale = 1;
		switch(ihp->bitpix) {
		case 16:
			ihp->datamax = 32767;
			ihp->datamin = -32767;
			break;
		case 32:
			ihp->datamax = 2147483647;
			ihp->datamin = -2147483647;
			break;
		default:
			ihp->datamax = 255;
			ihp->datamin = 0;
		}
	}
}
#if 0
/*
 * Start a standard image with the header values which cm and vc will need
 */
void StartStdImage(imnumber)
int imnumber;
{
# include <sys/time.h>
	register struct of_imagehdr *ihp;
	long clock;
	struct tm *tm;

	long time();
	struct tm *localtime();

	NewImage(imnumber);
	ihp = imageHdr[imnumber];
	ihp->simple = 1;
	ihp->bitpix = 16;
	clock = time(0);
	tm = localtime(&clock);
	sprintf(ihp->date_map,"%02d/%02d/%02d",tm->tm_mday,tm->tm_mon + 1,
		tm->tm_year);
}
#endif 0

/*
 * Set up for a new image in the given image number.  Allocate and clear
 * the header.
 */
void NewImage(imnumber)
register int imnumber;
{
	register struct of_imagehdr *ihp;

	CloseImage(imnumber);
	imageHdr[imnumber] = ihp = (struct of_imagehdr *)ChkAlloc(
		sizeof(struct of_imagehdr), "Image Header");
	bzero(ihp, sizeof(struct of_imagehdr));
	ihp->ftype = TMP;
	ihp->fd = -1;
}

#if 0
/*
 * Link an existing image to a new image number
 */
void LinkImage(imSrc, imDest)
register imSrc, imDest;
{
	ImNumErr(imSrc);
	ImNumErr(imDest);
	if( ! imageHdr[imSrc] )
		error("Image %d is not open", imSrc);
	CloseImage(imDest);
	imageHdr[imDest] = imageHdr[imSrc];
}
#endif 0

/*
 * Release the data buffer associated with an image.  If the imnumber is 
 * < 0, release the buffers for all images.
 */
void ReleaseImageBuff(imnumber)
register int imnumber;
{
	register int end;
	register struct of_imagehdr *ihp;

	if(imnumber < 0) {
		imnumber = 1;
		end = MAXIMAGES - 1;
	} else {
		ImNumErr(imnumber);
		end = imnumber;
	}
	for(;imnumber <= end; imnumber++) {
		if(imageHdr[imnumber] && (ihp = imageHdr[imnumber])->buf) {
			free( (char *)ihp->buf);
			ihp->buf = 0;
		}
	}
}

/*
 * Close a given image number.  If the same imageHdr is attached to another
 * imnumber, just clear the pointer.  Otherwise free the data buffer if any,
 * unlink the file if it is temporary, close the file, and release the
 * imageHdr.  Close all images if imnumber is < 0.
 */
void CloseImage(imnumber)
register int imnumber;
{
	register int end;
	register struct of_imagehdr *ih;

	if(imnumber < 0) {
		imnumber = 0;
		end = MAXIMAGES - 1;
	} else {
		ImNumErr(imnumber);
		end = imnumber;
	}
	for(; imnumber <= end; imnumber++) {
		if( (ih = imageHdr[imnumber]) == 0 )
			continue;
		if( imagelinked(imnumber) ) {
			imageHdr[imnumber] = 0;
			continue;
		}
		ReleaseImageBuff(imnumber);
		if(ih->fd > 0) {
			close(ih->fd);
		}
		if(ih->fname[0] && !(ih->ftype & EXTERNAL) ) {
			unlink(ih->fname);
		}
		free((char *)ih);
		imageHdr[imnumber] = (struct of_imagehdr *)0;
	}
}

#if 0
/*
 * Check to see that the two given images are similar enough to combine or
 * make a scatter plot.
 */
int ImageSimilar(im1, im2)
int im1, im2;
{
	register struct of_imagehdr *ihp1 = imageHdr[im1];
	register struct of_imagehdr *ihp2 = imageHdr[im2];
	double d1,d2;
	static double tol = 0.1;

	ChkImage(im1);
	ChkImage(im2);
	if( strcmp(ihp1->ctype1, ihp2->ctype1) || strcmp(ihp1->ctype2,
	    ihp2->ctype2) ) {
		fprintf(stderr, "%s has coordinates '%s' and '%s', but %s has '%s' and '%s'\n",
			ihp1->fname, ihp1->ctype1, ihp1->ctype2, ihp2->fname,
			ihp2->ctype1, ihp2->ctype2);
		warn(EMPTY);
	}
	d1 = fabs(ihp1->cdelt1);
	d2 = fabs(ihp1->cdelt2);
	if( cifdif(d1, fabs(ihp2->cdelt1), tol * d1 / ihp1->naxis1 ) ||
	    cifdif(d2, fabs(ihp2->cdelt2), tol * d2 / ihp1->naxis2 ) )
		error("%s and %s have different grids", ihp1->fname,
			ihp2->fname);
	if( cifdif(ihp1->crval1, ihp2->crval1, tol * d1) ||
	    cifdif(ihp1->crval2, ihp2->crval2, tol * d2) ) 
		error("%s and %s have different center positions", ihp1->fname,
			ihp2->fname);
	if( cifdif(ihp1->crota1, ihp2->crota1, 1.0) ||
	    cifdif(ihp1->crota2, ihp2->crota2, 1.0)  )
		error("%s and %s have different rotations", ihp1->fname, 
			ihp2->fname);
}
#endif 0

/*
 * Return the address of the requested line of pixels after filling the data
 * buffer if necessary.  The pointer will point to the pixel with the smallest
 * x in the line and x will increase as the pointer is incremented.  This
 * works with either an image in internal format or a FITS file.
 */
float * GetImageLine(imnumber,y,nlines,ydir)
int imnumber;
int y;				/* line # relative to reference position */
int nlines;			/* required number of lines in buffer */
int ydir;			/* extend buffer to larger or smaller y */
{
	register struct of_imagehdr *ihp = imageHdr[imnumber];
	register int i;
	int reqLine;		/* requested line number in image 0 -> bottom */

	if(nlines > ihp->naxis2)
		nlines = ihp->naxis2;
	if( !ihp->buf) {
		ihp->buflen = nlines;
		ihp->buf = (float *)ChkAlloc(ihp->buflen * ihp->naxis1 *
			sizeof(float) , "Image buffer");
		/* preposterous bufline to force a read */
		ihp->bufline = 1000000000;
	}
	/* If the file has cdelt2 < 0, crpix will also be expressed from the
	 * other end */
	if( ihp->cdelt2 > 0 ) {
		reqLine = round(ihp->crpix2) - 1 + y;
	} else {
		reqLine = ihp->naxis2 - round(ihp->crpix2) + y;
	}
	if( reqLine < ihp->bufline || reqLine >= ihp->bufline +
		ihp->buflen) {

		/* decide what part of the file to read in */
		if(ydir > 0) {
			if(reqLine  > (i = ( ihp->naxis2 - ihp->buflen)))
				ihp->bufline = i;
			else
				ihp->bufline = reqLine;
		} else {
			if( (i = reqLine - ihp->buflen + 1) < 0)
				ihp->bufline = 0;
			else
				ihp->bufline = i;
		}
		if(ihp->ftype & FITS) {
			int fitsLine;	/* Line in Fits array of first line of
					 * buffer (beg of data is 0) */
			int fitsWordLen;/* Len of an int in the fits file */
			int fitsBufLen;	/* Len of fits buffer */
			int line;	/* loop counter */
			register int pixel;	/* int value of Fits pixel */
			char *fb;	/* pointer to fits buffer */
			/* pointer to next byte in fitsbuffer */
			register char *ip;
			/* pointer to next loc in image buf */
			register float *op;
			/* Note, op will always advance regularly, but ip may
			 * have to move backward through lines or pixels if
			 * cdelt2 or cdelt1 is negative.  NextWord and/or
			 * nextLine may be zero.
			 */
			int nextWord;	/* Add to ip after each word */
			int nextLine;	/* Add to ip after each line */

			fitsWordLen = ihp->bitpix >> 3;
			fitsBufLen = ihp->buflen * ihp->naxis1 * fitsWordLen;
			fb = ChkAlloc(fitsBufLen, "FITS buffer");
			if(ihp->cdelt2 > 0)
				fitsLine = ihp->bufline;
			else
				fitsLine = ihp->naxis2 - ihp->bufline - nlines;
			if( lseek(ihp->fd, (long)(ihp->dataStart + fitsLine *
					ihp->naxis1 * fitsWordLen), 0) < 0) {
				sprintf(stderr, "Trouble seeking in %s because "
					,ihp->fname);
				perror(EMPTY);
				error(EMPTY);
			}
			if(read(ihp->fd, fb, fitsBufLen) != fitsBufLen)
				error("Trouble reading FITS file %s",
					ihp->fname);
			op = ihp->buf;
			if(ihp->cdelt2 > 0)
				ip = fb;
			else
				ip = fb + fitsBufLen - ihp->naxis1 *
					fitsWordLen;
			/* the common case - cdelt1 < 0 && cdelt2 > 0 */
			nextLine = 2 * ihp->naxis1 * fitsWordLen;
			if(ihp->cdelt1 < 0) {
				/* Start at line end */
				ip += (ihp->naxis1 - 1) * fitsWordLen;
				/* compensate for the increment reading bytes */
				nextWord = -2 * fitsWordLen;
				if(ihp->cdelt2 < 0)
					nextLine = 0;	/* No increment */
			} else {
				nextWord = 0;
				if(ihp->cdelt2 < 0)
					nextLine = -nextLine;
				else
					nextLine = 0;
			}
			for(line = 0; line < ihp->buflen; line++) {
			    for(i = 0; i < ihp->naxis1; i++) {
				switch(fitsWordLen) {
				case 4:
#if BYTEREVERSED
					pixel = *ip++ << 24;
					pixel |= (*ip++ & 0xff) << 16;
					pixel |= (*ip++ & 0xff) << 8;
					pixel |= *ip++ & 0xff;
#else
					pixel = *(long *)ip;
					ip += 4;
#endif
					break;
				case 2:
#if BYTEREVERSED
					pixel = *ip++ << 8;
					pixel |= *ip++ & 0xff;
#else
					pixel = *(short *)ip;
					ip += 2;
#endif
					break;
				case 1:
					pixel = *(unsigned char *)ip++;
					break;
				}
				if(pixel == ihp->blank)
					*op++ = FDRAGON;
				else
					*op++ = pixel * ihp->bscale + ihp->fbzero;
				ip += nextWord;
			    }
			    ip += nextLine;
			}
			free(fb);
		} else {
		    i = ihp->dataStart + ihp->bufline * ihp->naxis1
			* sizeof(float);
		    if( lseek(ihp->fd, (long)i, 0) < 0 ) {
			fprintf(stderr, "Trouble seeking in %s because ",
				ihp->fname);
			perror(EMPTY);
			error(EMPTY);
		    }
		    i = ihp->buflen *ihp->naxis1 * sizeof(float);
		    if( read(ihp->fd, (char *)ihp->buf, i) != i)
			error("Trouble reading %d bytes in %s",i,ihp->fname);
		}
	}
	return( &ihp->buf[ (reqLine - ihp->bufline) * ihp->naxis1 ]);
}

/*
 * Compute deltax or deltay given the corresponding cdelt.  Attempt to round
 * delta to an even number of Arc min.
 */
static double delta(cdelt)
double cdelt;
{
	register int arcMinPerPix = round(one / cdelt);
	if(cifdif(one / (double) arcMinPerPix, cdelt, .001 * cdelt))
		return(cdelt * sixty);
	else
		return(sixty / (double) arcMinPerPix);
}

/*
 * Get the x grid values.  Lowx and highx are in pixel numbers relative to the
 * center of the file.  They increase with increasing x coordinate.  Deltax
 * is the positive pixel spacing (in FitsAxisUnits if the coordinate is an
 *  angle), so lowx * delx is the most neg x offset from the center position.
 */
void ImageXGrid(imnumber,lowx,highx,deltax,xunit)
int imnumber;
int *lowx, *highx;		/* min and max pixel numbers */
double *deltax;			/* Absolute value of pixel spacing */
int *xunit;			/* Units of deltax (ARCMINUTES, etc.) */
{
	register struct of_imagehdr *ihp = imageHdr[imnumber];
	register int refpix;
	struct AXISNAMES units;

	refpix = round(ihp->crpix1);
	if( !(ihp->ftype & FITS) || ihp->cdelt1 > 0) {
		*lowx = 1 - refpix;
		*highx = ihp->naxis1 - refpix;
		*deltax = ihp->cdelt1;
	} else {
		*lowx = refpix - ihp->naxis1;
		*highx = refpix - 1;
		*deltax = -ihp->cdelt1;
	}
	units = FitsAxisUnits(ihp->ctype1);
	*xunit = units.unit;
	if( (units.type & 0xf) == SPATIAL )
		*deltax = delta( *deltax) * cunits[ARCMINUTES] /
			cunits[units.unit];
}

/*
 * Get the y grid values.  See coments above.
 */
void ImageYGrid(imnumber,lowy,highy,deltay,yunit)
int imnumber;
int *lowy, *highy;		/* min and max pixel (line) numbers */
double *deltay;			/* Absolute value of pixel spacing */
int *yunit;			/* units of deltay (ARCMINUTES, etc.) */
{
	register struct of_imagehdr *ihp = imageHdr[imnumber];
	register int refpix;
	struct AXISNAMES units;

	refpix = round(ihp->crpix2);
	if( !(ihp->ftype & FITS) || ihp->cdelt2 > 0) {
		*lowy = 1 - refpix;
		*highy = ihp->naxis2 - refpix;
		*deltay = ihp->cdelt2;
	} else {
		*lowy = refpix - ihp->naxis2;
		*highy = refpix - 1;
		*deltay = -ihp->cdelt2;
	}
	units = FitsAxisUnits(ihp->ctype2);
	*yunit = units.unit;
	if( (units.type & 0xf) == SPATIAL )
		*deltay = delta( *deltay) * cunits[ARCMINUTES] /
			cunits[units.unit];
}

/*
 * Check for a valid imnumber
 */
void ImNumErr(imnumber)
int imnumber;
{
	if(imnumber < 0 || imnumber >= MAXIMAGES)
		error("%d is an illegal image number",imnumber);
}

/*
 * Local routines
 * Check whether the imagehdr corresponding to imnumber is aka another imnumber
 */
static int imagelinked(imnumber)
register int imnumber;
{
	register int i;

	for(i = 1; i < MAXIMAGES; i++) {
		if(i != imnumber && imageHdr[i] == imageHdr[imnumber]) {
printf("imagelinked(%d) found linked\n",imnumber);
			return(1);
		}
	}
	return(0);
}

static int readhdr(ihp)
register struct of_imagehdr *ihp;
{
	register char *cardp;
	register struct of_fitskey *fk;
	register int started = 0;
	char hb[FITSRECL];
	char key[9];

	key[9] = '\0';
	lseek(ihp->fd, 0L, 0);
	for(;;) {
		if(read(ihp->fd, hb, FITSRECL) != FITSRECL) {
			printf("Trouble with read of header\n");
			return(1);
		}
		
		for(cardp = hb ; cardp < hb + FITSRECL; cardp +=80) {
			if( (fk = findkey(cardp)) != 0 ) {
				if(fk->type == LOGICAL && fk->index == END) {
					ihp->dataStart = lseek(ihp->fd, 0L, 1);
					return(0);
				}
				sethdr(fk, ihp, cardp);
				if( !started ) {
				    if( !ihp->simple) {
					printf("No SIMPLE card found\n");
					printf("first card was %.60s\n",cardp);
					return(1);
				    }
				    started++;
				}

			} else if( !started) {
				printf("First card not recognized");
				printf("first card was %.60s\n",cardp);
				return(1);
			}
		}
	}
}

static struct of_fitskey *findkey(cardp)
char *cardp;
{
	char key[9];
	register int i, high, low;
	register int c;

	if(strncmp(cardp, histKey, sizeof(histKey) - 1) == 0)
		cardp += sizeof(histKey);
	for(i = 0; i < 8; i++) {
		if(cardp[i] == ' ') {
			key[i] = '\0';
			break;
		} else {
			key[i] = cardp[i];
		}
	}
	key[8] = '\0';
	high = FITSKEYLEN - 1;
	for(low = 0; high >= low; ) {
		i = (high + low ) / 2;
		if( (c = strcmp(key, fitsKey[i].keyword )) == 0)
			return( &fitsKey[i]);
		if( c < 0)
			high = i - 1;
		else
			low = i + 1;
	}	
/* printf("%s not found\n", key); */
	return( (struct of_fitskey *)0 );
}

static void sethdr(fk, ih, cp)
register struct of_fitskey *fk;
struct of_imagehdr *ih;
char *cp;
{
	register char *ip = cp + 10;
	register char *end;
	register char *op;

	switch(fk->type) {
	case LOGICAL:
		ih->l[fk->index] = (ip[19] == 'T');
		break;
	case INT:
		sscanf(ip,"%d", &ih->i[fk->index]);
		break;
	case STR:
		op = ih->s[fk->index];
		for(end = ip + 10; *ip++ != '\''; )
			if(ip >= end)
				goto strerr;
		for(end = op + 11; op < end; ) {
			if(*ip == '\'') {
				*op = '\0';
				return;
			}
			*op++ = *ip++;
		}
strerr:
		*op = '\0';
/* error flagging when the string for the keyword is too long. IRAF does
 * this alot. eg CTYPE1 = 'GLAT--TAN         ' instead of
 *		 CTYPE1 = 'GLAT--TAN' 
 *  -- mwp 8/15/90 
 */
		fprintf(stderr,"Non-standard FITS card: %.50s\n", cp);
		break;
	case HISTORY:
		ip = index(cp + sizeof(histKey), ' ');
		while(*++ip == ' ')
			;
		for(end = cp + 80; *--end == ' '; )
			;
		if(end <= ip)
			goto strerr;
		bcopy(ip, ih->h[fk->index], end + 1 - ip);
		*++end = '\0';
	case DOUBLE:
		sscanf(ip,"%lf", &ih->d[fk->index]);
	}
}

#if 0
/*
 * The Routines which follow are for writing a FITS format file of an image
 */

/*
 * Write an image as a fits file assuming that the image header has been
 * completely filled in including DATAMAX, etc.
 */
void WriteFits(imnumber, fname)
int imnumber;
char *fname;
{
	register struct of_imagehdr * ihp = imageHdr[imnumber];
	register int i;

	ChkImage(imnumber);
	OpenFitsOut(fname);
	if(ihp->ftype & FITS) {	/* If it is a FITS file, just copy it */
		while( (i = read(ihp->fd,fitsBuf, FITSRECL)) == FITSRECL) {
			write(ffd, fitsBuf, FITSRECL);
		}
		if(i > 0)
			error("Incorrect FITS file length");
	} else {
		WriteFitsHdr(imnumber);
		WriteFitsData(imnumber);
	}
	CloseFitsOut();
}

/*
 * SetDataScale finds datamax and datamin by looking through the buffer and
 * sets FITS header values BLANK, BSCALE, BZERO, DATAMAX, and DATAMIN.  It
 * assumes that imageHdr[imnumber]->buf contains the whole image and that only
 * the first two dimansions are > 1.
 */
void SetDataScale(imnumber)
int imnumber;
{
	register struct of_imagehdr * ihp = imageHdr[imnumber];
	register float *pix, *endPix;

	pix = ihp->buf;
	if(ihp->naxis >= 2)
		endPix = pix + ihp->naxis1 * ihp->naxis2;
	else
		endPix = pix + ihp->naxis1;
	while(*pix == FDRAGON)
		if(++pix >= endPix)
			error("No data in image %d\n",imnumber);
	ihp->datamax = ihp->datamin = *pix;
	while(++pix < endPix) {
	    if(*pix != FDRAGON)
		if(*pix > ihp->datamax)
			ihp->datamax = *pix;
		else if(*pix < ihp->datamin)
			ihp->datamin = *pix;
	}
	ihp->blank = (-1) << ihp->bitpix - 1;
	ihp->fbzero = (ihp->datamax + ihp->datamin) / 2.;
	if(ihp->bitpix == 16)
		ihp->bscale = (ihp->datamax - ihp->datamin) / 65534.;
	else if(ihp->bitpix == 32)
		ihp->bscale = (ihp->datamax - ihp->datamin) / 4294967294.;
	else if(ihp->bitpix == 8)
		ihp->bscale = (ihp->datamax - ihp->datamin) / 254.;
	else
		error("Can't use BITPIX = %d\n", ihp->bitpix);
}

static void OpenFitsOut(fname)
char *fname;
{
	if(strncmp(fname, "/dev/rmt", 8) == 0) {
		ffd = SafeOpen("Image File",fname,O_WRONLY , 0666);
	} else {
	    if( access(fname, F_OK) == 0) {
		if( access(fname, W_OK) == 0)
			warn("%s exists and will be overwritten\n", fname);
		else
			error("%s exists and can't be written", fname);

	    }
	    ffd = SafeOpen("Image File", fname, O_WRONLY | O_CREAT | O_TRUNC,
		0666);
	}
	headPos = fitsBuf = SafeAlloc(FITSRECL + 1,"Fits Buffer");
	cfill(fitsBuf,' ',FITSRECL);
}

static void CloseFitsOut()
{
	SafeClose(ffd);
	ffd = -1;
	SafeFree(fitsBuf);
	fitsBuf = 0;
}

/*
 * Write header data from the given imnumber into the fits file which is
 * presumed open.
 */
static void WriteFitsHdr(imnumber)
int imnumber;
{
	struct of_fitskey *list[FITSKEYLEN];
	register struct of_imagehdr *ihp = imageHdr[imnumber];
	register struct of_fitskey *fkp;
	register int i;

	for(fkp = fitsKey; fkp < fitsKey + FITSKEYLEN; fkp++)
		list[fkp->sequence] = fkp;
	for(i = 0; i < FITSKEYLEN; i++) {
		fkp = list[i];
		if(fkp->special) {
			switch(fkp->special) {
			case 1:
			case 2:
			case 3:
			case 4:
				if(ihp->naxis < fkp->special) {
					i += MAXAXES - fkp->special;
					continue;
				}
				break;
			case 5:
			case 6:
			case 7:
			case 8:
				if(ihp->naxis < fkp->special - 4) {
					i += 24 + (MAXAXES - fkp->special) * 5;
					continue;
				}
				break;
			case 9:
			/* This is the END card which has no value and should
			 * be preceeded by the extra cards */
				FitsCard(fkp->keyword);
				WriteFitsBuf();
				return;
			/* case 10 - output only if non zero is handled later */
			}
		}
		switch(fkp->type) {
		case LOGICAL:
			LFitsCard(fkp->keyword, ihp->l[fkp->index]);
			break;
		case INT:
			IFitsCard(fkp->keyword, ihp->i[fkp->index]);
			break;
		case STR:
			if(fkp->special != 10 || *(ihp->s[fkp->index]) != 0)
				SFitsCard(fkp->keyword, ihp->s[fkp->index]);
			break;
		case HISTORY:
			if(fkp->special != 10 || *(ihp->h[fkp->index]) != 0)
				HFitsCard(fkp->keyword, ihp->h[fkp->index]);
			break;
		case DOUBLE:
			if(fkp->special != 10 || ihp->d[fkp->index] != 0)
				DFitsCard(fkp->keyword, ihp->d[fkp->index]);
			break;
		default:
			error("Internal error: Bad FITS keyword type %d",
				fkp->type);
		}
	}
}

static void WriteFitsData(imnumber)
int imnumber;
{
	register struct of_imagehdr *ihp = imageHdr[imnumber];
	register int i;				/* scaled in value of pixel */
	register char *op;			/* output pointer */
	register float *pix, *endPix;		/* current pixel, last pixel of
						 * current line in image buf */
	char *endBuf;				/* end of fits buffer */
	int xdir, ydir;				/* order of pix and lines in
						 * FITS file */
	int y, endy;				/* range of lines in image */
	double deltay;				/* extra output of ImageYGRid */
	int yunit;				/* more extra output */

	ImageYGrid(imnumber, &y, &endy, &deltay, &yunit);
	if(ihp->cdelt2 > 0) {
		ydir = 1;
		endy++;
	} else {
		ydir = -1;
		i = y;
		y = endy;
		endy = i - 1;
	}
	xdir = (ihp->cdelt1 > 0)? 1: -1;
	op = fitsBuf;
	endBuf = op + FITSRECL;

	for( ; y != endy; y += ydir) {
	    pix = GetImageLine(imnumber, y, 5, ydir);
	    if(ihp->cdelt1 > 0) {
		endPix = pix + ihp->naxis1;
	    } else {
		endPix = --pix;
		pix += ihp->naxis1;
	    }
	    for( ;pix != endPix; pix += xdir) {
		if(*pix == FDRAGON)
		    i = ihp->blank;
		else
		    i = (*pix - ihp->fbzero) / ihp->bscale;

		if(ihp->bitpix == 16) {
# if BYTEREVERSED
		    *op++ = i >> 8;
		    *op++ = i;
# else
		    *(short *)op  = i;
		    op += 2;
# endif BYTEREVERSED
		} else if(ihp->bitpix == 32) {
# if BYTEREVERSED
		    *op++ = i >> 24;
		    *op++ = i >> 16;
		    *op++ = i >> 8;
		    *op++ = i;
# else
		    *(long *)op  = i;
		    op += 4;
# endif BYTEREVERSED
		} else {
		    *op++ = i;
		}
		if(op >= endBuf) {
		    WriteFitsBuf();
		    op = fitsBuf;
		}
	    }
	}
	if(op != fitsBuf) {
	    bzero(op, endBuf - op);
	    WriteFitsBuf();
	}
}

void FFitsCard(name,value)
char *name;
float value;
{
	sprintf(headPos,"%-8s= %#20.9G /",name,value);
	headPos[32] = ' ';
	NextCrd();
}

void DFitsCard(name,value)
char *name;
double value;
{
	sprintf(headPos,"%-8s= %#20.9G /",name,value);
	headPos[32] = ' ';
	NextCrd();
}

void IFitsCard(name,value)
char *name;
int value;
{
	sprintf(headPos,"%-8s= %20d /",name,value);
	headPos[32] = ' ';
	NextCrd();
}

void SFitsCard(name,value)
char *name;
char *value;
{
	char tempstr[80];

	sprintf(tempstr,"%-8s= '%-8s'%11s/",name,value,EMPTY);
	FitsCard(tempstr);
}

void HFitsCard(name,value)
char *name;
char *value;
{
	sprintf(headPos,"%s %s %.57s",histKey, name, value);
	headPos[strlen(headPos)] = ' ';
	NextCrd();
}

void LFitsCard(name,value)
char *name;
int value;
{
	value = (value)?'T':'F';
	sprintf(headPos,"%-8s= %20c /",name,value);
	headPos[32] = ' ';
	NextCrd();
}

/* copy up to 80 chars into the current header position */
void FitsCard(card)
char *card;
{
	register char *ip = card, *op = headPos;
	register char *end = op + 80;
	while(*ip && op < end)
		if(islower(*ip))
			*op++ = toupper(*ip++);
		else
			*op++ = *ip++;
	NextCrd();
}

static void NextCrd()
{
	if( (headPos += 80) > fitsBuf + FITSRECL) {
		WriteFitsBuf();
		cfill(fitsBuf,' ',FITSRECL);
	}
}

static void WriteFitsBuf()
{
	if(write(ffd,fitsBuf,FITSRECL) != FITSRECL)
		error("Trouble writing FITS file");
	headPos = fitsBuf;
}

char *FitsAxisLabel(axis)
int axis;
{
	static char name[3][2][2][5] = {
		{{"RA  ","DEC"},{"RA--","DEC-"}},
		{{"GLON","GLAT"},{"GLON","GLAT"}},
		{{"DX  ","DY  "},{"DX--","DY--"}}
	};
	static char suffix[4][2][5] = {
		{"    ","    "},{"-SIN","-SIN"},
		{"---X","---X"},{"-TAN","-TAN"}
	};
	static char label[9];
	register char *suf;

	suf = suffix[cproj.type][axis];
	strcpy(label,name[csys.type][*suf != ' '][axis]);
	strcpy(label + 4, suf);
	return(label);
}
#endif 0

struct AXISNAMES FitsAxisUnits(ctype)
register char *ctype;
{
	static struct AXISNAMES axisNames[] = {
		{60.,ARCMINUTES,SPATIAL | OFFSETTYPE,"RA",2,"R.A. Offset"},
		{60.,ARCMINUTES,SPATIAL | OFFSETTYPE | LATTYPE,"DEC",3,
			"Dec. Offset"},
		{1.,DEGREES,SPATIAL,"GLON",4,"l"},
		{1.,DEGREES,SPATIAL | LATTYPE,"GLAT",4,"b"},
		{1e-3,KMPERSEC,VELOCITY,"VELO",4,"Velocity"},
		{1e-3,KMPERSEC,VELOCITY,"FELO",4,"Velocity"},
		{1,SECONDS,TIME,"TIME",4,"Time"},
		{1e-6,MEGAHERTZ,FREQUENCY,"FREQ",4,"mHz"},
		{1e6,MICROMETERS,WAVELENGTH,"LAMBD",5,"Microns"},
		{1.,NONE,UNDEF,"PIXEL",5,"Pix"},
		{1.,NONE,UNDEF,"     ",5,"Pix"},
	};
	struct AXISNAMES tempName;
	register int n;

	for(n = 0; n < sizeof(axisNames)/sizeof(struct AXISNAMES); n++){
		if( !strncmp(ctype, axisNames[n].fname, axisNames[n].len) )
			return(axisNames[n]);
	}
	/* if we didn't find it, assume degrees */
	tempName = axisNames[2];
	for(n = 0; n < 8; n++) {
		if(ctype[n] == ' ')
			break;
		tempName.lname[n] = ctype[n];
		if(n < tempName.len)
			tempName.fname[n] = ctype[n];
	}
	tempName.lname[n] = '\0';
	return(tempName);
}

/************************************************************
 * ChkAlloc - call malloc and check for no available memory *
 ************************************************************/

char *ChkAlloc(size,name)
int size;
char *name;
{
	register char *cp;

	if( !(cp = malloc((unsigned int)size)))
		error("Can't allocate %d bytes for %s",size,name);
	return(cp);
}

/*VARARGS*/
error(va_alist)
va_dcl
{
	register char *s;
	va_list ap;

	va_start(ap);
	s = va_arg(ap, char*);
	_doprnt(s, ap, stderr);
	va_end(ap);
	putc('\n', stderr);
	exit(77);
}
static double half = 0.5;
round(d)
double d;
{
	return((int)((d < 0.0)? d - half:d + half));
}

/*
 * This subroutine compares a and b. If they differ by more than fuzz
 * it returns 1, otherwise 0. A,b,and fuzz are double.
 */

cifdif(a,b,fuzz)
double a,b,fuzz;
{
	register double diff;
	diff = a - b;
	if(diff < 0)
		diff = - diff;
	return(diff > fuzz);
}

#include "ql.h"

int colorMapSize = 256;
extern Frame frame;
extern Canvas canvas;

static int *blob = (int *)0;		/* Blob points to a null terminated set
				 * of offsets from the current pixel for
				 * filling a zoom x zoom size blob in a
				 * pixrect
				 */
static void ShowAngScale(/*pr, x, y, width, xdegrees*/);
static double ChooseStep(/*minStep*/);

/* Read a fits file into a memory pixrect in p. Make the memory pixrect
 * to hold the data if needed.  If maxv != minv, scale the image so that
 * pixels between minv and maxv use the range 0-(colorMapSize - 1).  Otherwise
 * use DATAMAX and DATAMIN from fits.  If this is a 32 bit frame buffer,
 * offset controls where the 8 bit data is put.  Offsets of 0, 1, 2, or
 * 3 write grey, blue, green, or red values.
 */
void LoadImage(p, use24Bits, offset)
PICTURE *p;
int use24Bits, offset;
{
	double scale, /* deltax, */ deltay, /* xunit, */ yunit;
	float *xp;
	Pixrect *pr;
	struct of_imagehdr *ihp;
	int cmLimit = colorMapSize - 1;
	int imnumber = 0;
	int width, height;	/* image width and height in pixrect pixels */
	int val;
	int y1, y2;
	int ix, iy;
	char *pix0, *pixLine, *pix;
	int linebytes;		/* # of bytes in a scan line of the pixrect */
	int zoomlinebytes;	/* # of bytes in a line of the original image */
	int pixelbytes;		/* # of bytes in a pixel of the pixrect */
	int zoompixelbytes;	/* # of bytes in zoom pixels of the pixrect */
	int widthbytes;		/* # of bytes of the pixrect which are used
				 * in a scan line <= linebytes */
	int *blobp;

	OpenFits(imnumber, p->fn);
	ihp = imageHdr[imnumber];
	if(ihp->naxis < 2 || ihp->naxis1 < 2 || ihp->naxis2 < 2)
		error("Not enough axes or points in FITS file\n");
	width = ihp->naxis1 * zoom;
	height = ihp->naxis2 * zoom;
	if(ihp->cdelt1 == 0)
		ihp->cdelt1 = 1;
	if(ihp->cdelt2 == 0)
		ihp->cdelt2 = 1;

	if(minv == maxv) {
		maxv = ihp->datamax;
		minv = ihp->datamin;
	}
	if(!p->pr)
		p->pr = mem_create(width, height + 1 + WEDGE_HEIGHT,
			(use24Bits)? 32: 8);
	pr = p->pr;
	if(! pr)
		error("Can't create a %d x %d pixrect", width, height + zoom +
			 1 + WEDGE_HEIGHT);
	if(pr->pr_depth != ((use24Bits)? 32: 8))
		error("LoadImage: mis-matched pr_depth");
	linebytes = mpr_d(pr)->md_linebytes; 
	zoomlinebytes = linebytes * zoom; 
	pixelbytes = (use24Bits)? 4: 1;
	zoompixelbytes = pixelbytes * zoom;
	widthbytes = pixelbytes * width;
	/* The blob is to the upper right of the base pixel */
	if(!blob) {
		/* allocate enough memory for 32 bit pixrects with 24 bits used
		 * so we don't have to keep track of it.
		 */
		blob = (int *)ChkAlloc((zoom * zoom * 3 + 1) * sizeof(int),
			"blob");
	}
	blobp = blob;
	if(pixelbytes == 1)
		offset = 0;
	else if(offset == 0)
		offset = -1;
	if(offset < -1 || offset > 3)
		error("Bad offset");
	for(iy = 0; iy > - zoomlinebytes; iy -= linebytes)
	    for(ix = 0; ix < zoompixelbytes; ix += pixelbytes)
		if(offset >= 0) {
			*blobp++ = ix + iy + offset;
		} else {
			*blobp++ = ix + iy + 1;
			*blobp++ = ix + iy + 2;
			*blobp++ = ix + iy + 3;
		}
	*blobp = 0;
	scale = cmLimit / (maxv - minv);
/*	ImageXGrid(imnumber, &x1, &x2, &deltax, &xunit); */
	ImageYGrid(imnumber, &y1, &y2, &deltay, &yunit);
	pix0 =  (char *)mpr_d(pr)->md_image;
	for(iy = y1, pixLine = pix0 + (height - 1) * linebytes;
		iy <= y2; iy++, pixLine -= zoomlinebytes) {
	    xp = GetImageLine(imnumber,iy, 32, 1);
	    if(ihp->cdelt1 < 0) {
		for(pix = pixLine + widthbytes - zoompixelbytes;
			pix >= pixLine; pix -= zoompixelbytes) {
		    if(*xp < 1e36) {
			val = scale * (*xp++ - minv);
			if(val < 0)
				val = 0;
			if(val > cmLimit)
				val = cmLimit;
			blobp = blob;
			do {
				pix[*blobp] = val;
			} while(*++blobp);
		    } else {
			xp++;
		    }
		}
	    } else {
		for(pix = pixLine; pix < pixLine + widthbytes; pix +=
			zoompixelbytes) {
		    if(*xp < 1e36) {
			val = scale * (*xp++ - minv);
			if(val < 0)
				val = 0;
			if(val > cmLimit)
				val = cmLimit;
			blobp = blob;
			do {
				pix[*blobp] = val;
			} while(*++blobp);
		    } else {
			xp++;
		    }
		}
	    }
	}
	ReleaseImageBuff(imnumber);
	make_wedge(pr, 0, height + 1, width, WEDGE_HEIGHT, minv, maxv);
	ShowAngScale(pr, 0, height + 1, width, (ihp->naxis1 - 1) *
		fabs(ihp->cdelt1));
}

/* Make and label a ramp of colors */
void make_wedge(pr, x, y, w, h, minLabel, maxLabel)
Pixrect *pr;
int x, y;				/* Origin in pixrect */
int w, h;
double minLabel, maxLabel;		/* value of min, max color */
{
	int i;				/* Current color in wedge */
	int ix, iy;			/* Current x, y in wedge */
	int dx, dy;			/* width, hieght of each grey value */
	char pbuf[32];			/* print buffer */

	/* Clear space for ramp and its labels */
	pr_rop(pr, x, y, w, h, PIX_CLR, (Pixrect *)0, 0, 0);
	iy = y + pf->pf_defaultsize.y;
	sprintf(pbuf, "%.3g", minLabel);
	pr_text(pr, x, iy, PIX_SRC | PIX_COLOR(colorMapSize-1), pf, pbuf);
	pr_text(pr, x, iy, PIX_SRC , pf, pbuf);
	
	sprintf(pbuf, "%.3g", maxLabel);
	ix = x + w - pf->pf_defaultsize.x * strlen(pbuf);
	pr_text(pr, ix, iy, PIX_SRC | PIX_COLOR(colorMapSize-1), pf, pbuf);
	pr_text(pr, ix, iy, PIX_SRC , pf, pbuf);
	dx = (w + colorMapSize - 1) / colorMapSize;
	if(dx <= 0)
		dx = 1;
	iy += pf->pf_defaultsize.y / 3;
	dy = h + y - iy;
	for(i = 0; i < colorMapSize; i++) {
		ix = i * w / colorMapSize;
		pr_rop(pr, i * w / colorMapSize, iy, dx, dy,
			PIX_SRC | PIX_COLOR((i | (i << 8) | (i << 16))),
			(Pixrect *)0, 0, 0);
	}
}

static void ShowAngScale(pr, x, y, width, xdegrees)
Pixrect *pr;
int x, y;				/* Origin in pixrect */
int width;				/* width in pixels */
double xdegrees;			/* width in degrees on sky */
{
#define MINSIZ 0.35
	double size;			/* mark width in Deg or 'Arc */
	int useMinutes;
	int ix, iy, center;
	int markwid;			/* width of the mark in pixels */
	int markstart;
	int textwid;
	int halfheight;
	char pbuf[32];			/* print buffer */

/*printf("ShowAngScale: xdeg = %g\n", xdegrees);*/
	if(xdegrees > 0.5 /MINSIZ) {
		useMinutes = 0;
	} else {
		useMinutes = 1;
		xdegrees *= 60.;
	}
	size = ChooseStep(xdegrees * MINSIZ);
	markwid = width * size / xdegrees;
	center = x + width / 2;
	iy = y + pf->pf_defaultsize.y;
	sprintf(pbuf, (useMinutes)? "%.3g 'Arc": "%.3g Deg.", size);
	textwid = pf->pf_defaultsize.x * strlen(pbuf);
	ix = center - textwid / 2;
	halfheight = pf->pf_defaultsize.y / 2;

	pr_text(pr, ix, iy, PIX_SRC , pf, pbuf);
	/*pr_text(pr, ix, iy, PIX_SRC | PIX_COLOR(colorMapSize-1), pf, pbuf);*/
	markstart = center - markwid / 2;
	pr_vector(pr, markstart, iy - halfheight, markstart, iy, PIX_SET, 0);
	/*pr_vector(pr, markstart, iy - halfheight, markstart, iy,
                PIX_SRC|PIX_COLOR(colorMapSize-1), 0);*/
	markstart += markwid;
	pr_vector(pr, markstart, iy - halfheight, markstart, iy, PIX_SET, 0);
	/*pr_vector(pr, markstart, iy - halfheight, markstart, iy,
		PIX_SRC|PIX_COLOR(colorMapSize-1), 0);*/
	if(markwid -textwid > 20) {
		iy -= halfheight / 2;
		/*pr_vector(pr, center - markwid / 2, iy, ix, iy, 
                          PIX_SRC|PIX_COLOR(colorMapSize-1), 0);*/
		pr_vector(pr, center - markwid / 2, iy, ix, iy, PIX_SET, 0);
		/*pr_vector(pr, center + markwid / 2, iy, ix + textwid, iy,
			  PIX_SRC|PIX_COLOR(colorMapSize-1), 0);*/
		pr_vector(pr, center + markwid / 2, iy,ix+textwid, iy, PIX_SET, 0);
	}
}

static double ChooseStep(minStep)
double minStep;
{
	static double steps[] = { 1.0, 2.0, 2.5, 4.0, 5.0, 10.0};
	double fact, *stepp;
	int expon;

	expon = (int)floor(log10(minStep));
/*printf("minStep = %g, expon = %d\n",  minStep, expon);*/
	minStep /= (fact = exp10((double)expon));
	for(stepp = steps; *stepp <= minStep; stepp++)
		;
	return( *stepp * fact);
}

void StartDisplay()
{
	window_set(frame, WIN_SHOW, TRUE, 0);
	notify_do_dispatch();
}

int Is24Bit()
{
	/* check for read/write permission and existance */
	/* this isn't the recommended way, but Sun's way would be hard */
	return(access("/dev/cgeight0", 6) == 0);
}

Pixrect *CreatePr24()
{
	Pixrect *pr;

	if( !(pr = mem_create(xsize, ysize, 32)) )
		error("Can't create %d by %d pixrect 32 bits deep");
	return(pr);
}

ShowCurrent(pw)
Pixwin *pw;
{
	PICTURE *dPic = (  (is24Bit)? &pic24: &pic[curPic]);
	Pixrect *pr = dPic->pr;

	window_set(frame, FRAME_LABEL, dPic->label, 0);
        if(invert_pixels==-1)
           pw_rop(pw,0,0, pr->pr_size.x, pr->pr_size.y, PIX_NOT(PIX_DST), pr,0,0); 
        else
           pw_rop(pw, 0, 0, pr->pr_size.x, pr->pr_size.y, PIX_SRC, pr, 0, 0); 
}

/* Put the picture in the srcPic in the 32 bit pixrect (pr24).  How determines
 * what will be done.  If how ==:
 * 0 - Make a grey image.  Duplicate the 8 bit pixrect in red, green ,and blue
 * 1, 2, or 3 - Put the 8 bit pixrect in Blue, green, or red.
 * 4 - Shift green and red into blue and green and load red from src
 * 5 - Shift blue and green into green and red and load blue from src
 */
ConvertToPr24(srcPic, how)
PICTURE *srcPic;
int how;
{
	PICTURE *dstPic = &pic24;
	unsigned char *srcPix, *endPix, *endLinePix;
	unsigned int *dstPix;
	int srclinebytes, dstLineExtra;
	int width;
	int i;

	if((width = srcPic->pr->pr_size.x) != dstPic->pr->pr_size.x ||
	    srcPic->pr->pr_size.y != dstPic->pr->pr_size.y)
		error("Size mismatch between 32 bit and 8 bit pixrects");

	strcpy(dstPic->label, srcPic->label);

	srclinebytes = mpr_d(srcPic->pr)->md_linebytes; 
	dstLineExtra = mpr_d(dstPic->pr)->md_linebytes / 4 - width;
	srcPix =  (unsigned char *)mpr_d(srcPic->pr)->md_image;
	endPix = srcPix + srcPic->pr->pr_size.y * srclinebytes;
	dstPix =  (unsigned int *)mpr_d(dstPic->pr)->md_image;

	for( ; srcPix < endPix; srcPix += srclinebytes - width,
		dstPix += dstLineExtra) {

	    endLinePix = srcPix + width;
	    switch(how) {
	    case 0:
		for( ; srcPix < endLinePix; srcPix++, dstPix++) {
			i = *srcPix;
			*dstPix = i | i << 8 | i << 16;
		}
		break;
	    case 1:
	    case 2:
	    case 3:
		for( ; srcPix < endLinePix; srcPix++, dstPix++) {
			( (char *)dstPix)[how] = *srcPix;
		}
		break;
	    case 4:
		for( ; srcPix < endLinePix; srcPix++, dstPix++) {
			*dstPix = *srcPix | (*dstPix << 8);
		}
		break;
	    case 5:
		for( ; srcPix < endLinePix; srcPix++, dstPix++) {
			*dstPix = (*srcPix << 16) | (65535 & (*dstPix >> 8));
		}
		break;
	    }
	}
}

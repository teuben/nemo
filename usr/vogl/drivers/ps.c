
#undef VOGLE
/*
 * The next bit is site and UNIX specific...
 */
#undef LASERWRITER

/*
 * Some more bug fixes from  ralf@physik3.gwdg.de (Ralf Fassel)
 * regarding %%Pages and devname .... 28/03/94
 */

#include <stdio.h>
#ifdef VOGLE
#include "vogle.h"
#else
#include "vogl.h"
#endif

typedef struct {
	float	r, g, b;
} Coltab;

static Coltab	*coltab;
extern char	*vallocate();

extern FILE	*_voutfile();

static int	ps_first_time = 1, drawn = 0,
		curcol = 0,			/* black */
		pslstx = -1, pslsty = -1,	/* last (x, y) drawn */
		colour = 0, page;

/*
 * gray scale map for our standard colours
 */
static float	graymap[8] = {
			0.0,
			0.30,
			0.59,
			0.89,
			0.11,
			0.41,
			0.70,
			0.99
};

static FILE	*fp;

/*
 * noop
 *
 *	do nothing but return -1
 */
static int
noop()
{
	return(-1);
}
/*
 * PS_color
 *
 *	change the grey value of the ink
 */
static int
PS_color(col)
	int	col;
{
	curcol = col;
	if (colour) {
		curcol %= 256;
		fprintf(fp, "%3.2f %3.2f %3.2f c\n", coltab[curcol].r, coltab[curcol].g, coltab[curcol].b);
		return(0);
	}

	if (col > 7)
		return(0);


#ifdef GREY_LINES
	fprintf(fp, "%3.2f g\n", graymap[curcol]);
#endif

	return(0);
}

/*
 * PS_mapcolor
 *
 *	Set our values in our pseudo colour map.
 */
static int
PS_mapcolor(indx, r, g, b)
	int	indx, r, g, b;
{
	if (colour && indx < 256 && indx >= 0) {
		coltab[indx].r = r / 255.0;
		coltab[indx].g = g / 255.0;
		coltab[indx].b = b / 255.0;
	}

	return(0);
}

/*
 * PS_common_init
 *
 *	 Initialization that is common to both layouts
 */
static int
PS_common_init()
{

	vdevice.depth = colour ? 8 : 1;

	/*	Set other line drawing parameters	*/

	fprintf(fp, "2 setlinewidth\n1 setlinejoin\n1 setlinecap\n");

	/*	Speed up symbol font handling	*/

	fprintf(fp, "/sf /Courier findfont def\n");

	/*	Move	*/

	fprintf(fp, "/m /moveto load def\n");

	/*	Draw	*/

	fprintf(fp, "/d { lineto currentpoint stroke moveto } def\n");

	/*	Polygon Draw	*/

	fprintf(fp, "/p /lineto load def\n");

	/*	Set character height	*/

	fprintf(fp, "/h { sf exch scalefont setfont } def\n");

	/*	Show character string	*/

	fprintf(fp, "/s /show load def\n");

	/*	Set gray scale	*/

	fprintf(fp, "/g /setgray load def\n");

	if (colour) {
		fprintf(fp, "/c /setrgbcolor load def\n");
		coltab = (Coltab *)vallocate(256 * sizeof(Coltab));
		PS_mapcolor(0, 0, 0, 0);
		PS_mapcolor(1, 255, 0, 0);
		PS_mapcolor(2, 0, 255, 0);
		PS_mapcolor(3, 255, 255, 0);
		PS_mapcolor(4, 0, 0, 255);
		PS_mapcolor(5, 255, 0, 255);
		PS_mapcolor(6, 0, 255, 255);
		PS_mapcolor(7, 255, 255, 255);
	}

	/*	Set a default font height	*/
	
	fprintf(fp, "45 h\n");

	return(1);
}

/*
 * PS_init
 *
 *	set up the postcript environment. Returns 1 on success.
 */
static int
PS_init()
{
	fp = _voutfile();

	if (!ps_first_time)
		return(1);

	page = 1;
	fputs("%!PS-Adobe-2.0 EPSF-1.2\n", fp);
	fputs("%%BoundingBox: 74 96 528 728\n", fp);
	fprintf(fp, "%%%%Page: %d %d\n", page, page);
	fputs("%%EndComments\n", fp);
	fprintf(fp, "72 300 div dup scale\n90 rotate\n400 -2200 translate\n");

	vdevice.sizeSy = 1890; 
	vdevice.sizeSx = 2634; 
	vdevice.sizeX = vdevice.sizeY = 1890; 

	PS_common_init();

	return (1);
}

/*
 * PSP_init
 *
 *	set up the postscript (Portrait) environment. Returns 1 on success.
 */
static int
PSP_init()
{
	fp = _voutfile();

	if (!ps_first_time)
		return(1);

	page = 1;
	fputs("%!PS-Adobe-2.0 EPSF-1.2\n", fp);
	fputs("%%BoundingBox: 72 96 526 728\n", fp);
	fprintf(fp, "%%%%Page: %d %d\n", page, page);
	fputs("%%EndComments\n", fp);

	fprintf(fp, "72 300 div dup scale\n300 400 translate\n");

	vdevice.sizeSy = 2634; 
	vdevice.sizeSx = 1890; 
	vdevice.sizeX = vdevice.sizeY = 1890; 

	PS_common_init();

	return (1);
}

/*
 * PS_exit
 *
 *	do a showpage and close the output file if neccessary.
 */
static int
PS_exit()
{
	fputs("showpage\n", fp);
	fputs("%%Trailer\n", fp);
	fflush(fp);

	if (fp != stdout)
		fclose(fp);

	/*
	 * Bug fix from brett@kirk.es.go.dlr.de (Bernward Bretthauer)
	 */
	ps_first_time = 1;
	drawn = 0;
	curcol = 0;			/* black */
	pslstx = pslsty = -1;
	colour = 0;

	return(0);
}

/*
 * PS_draw
 *
 *	draw to an x, y point.
 */
static int
PS_draw(x, y)
	int	x, y;
{
	if (pslstx != vdevice.cpVx || pslsty != vdevice.cpVy)
		fprintf(fp, "%d %d m\n", vdevice.cpVx, vdevice.cpVy);

	fprintf(fp, "%d %d d\n", x, y);
	pslstx = x;
	pslsty = y;
	drawn = 1;

	return(0);
}

static int
PS_pnt(x, y)
	int	x, y;
{
	fprintf(fp, "%d %d m\n", x, y);
	fprintf(fp, "%d %d d\n", x, y);

	return(0);
}


/*
 * PS_font
 *
 * load in small or large - could be improved.
 */
static int
PS_font(font)
	char	*font;
{
	if (strcmp(font, "small") == 0) {
		vdevice.hwidth = 22.0;
		vdevice.hheight = vdevice.hwidth * 1.833;
		fprintf(fp, "%d h\n", (int)vdevice.hheight);
	} else if (strcmp(font, "large") == 0) {
		vdevice.hwidth = 35.0;
		vdevice.hheight = vdevice.hwidth * 1.833;
		fprintf(fp, "%d h\n", (int)vdevice.hheight);

	} else
		return(0);

	return(1);
}

/*
 * PS_clear
 *
 *	flush the current page without resetting the graphics state of the
 * laser printer.
 */
static int
PS_clear()
{
	if (drawn) {
		fprintf(fp, "gsave showpage grestore\n");
		/* This is the end of the page, not of the document. */
		/*  ralf@physik3.gwdg.de (Ralf Fassel) */
		fputs("%%PageTrailer\n", fp);
		page++;
		fprintf(fp, "%%%%Page: %d %d\n", page, page);
	}

	drawn = 0;

	return(0);
}

	
/*
 * PS_char
 *
 *	output a character making sure that a '\' is sent first when
 * appropriate.
 */
static int
PS_char(c)
	char	c;
{
	if (pslstx != vdevice.cpVx || pslsty != vdevice.cpVy)
		fprintf(fp, "%d %d m\n", vdevice.cpVx, vdevice.cpVy);

	fprintf(fp, "(");

	switch(c) {
	case '(':
		fprintf(fp, "\\(");
		break;
	case ')':
		fprintf(fp, "\\)");
		break;
	case '\\':
		fprintf(fp, "\\");
		break;
	default:
		fprintf(fp, "%c",c);
	}

	fprintf(fp,") s \n");

	drawn = 1;
	pslstx = pslsty = -1;

	return(0);
}

/*
 * PS_string
 *
 *	output a string one char at a time.
 */
static int
PS_string(s)
	char	*s;
{
	char	c;

	if (pslstx != vdevice.cpVx || pslsty != vdevice.cpVy)
		fprintf(fp, "%d %d m\n", vdevice.cpVx, vdevice.cpVy);

	fprintf(fp, "(");
	while ((c = *s++))
		switch(c) {
		case '(':
			fprintf(fp, "\\(");
			break;
		case ')':
			fprintf(fp, "\\)");
			break;
		case '\\':
			fprintf(fp, "\\");
			break;
		default:
		fprintf(fp, "%c",c);
		}

	fprintf(fp,") s \n");
	drawn = 1;
	pslstx = pslsty = -1;

	return(0);
}

/*
 * PS_fill
 *
 *      fill a polygon
 */
static int
PS_fill(n, x, y)
	int     n, x[], y[];
{
	int     i;


	fprintf(fp, "newpath \n");

	fprintf(fp, "%d %d m\n", x[0], y[0]);

	for (i = 1; i < n; i++)
		fprintf(fp, "%d %d p\n", x[i], y[i]);

	fprintf(fp, "closepath\n");

	if (!colour)
		fprintf(fp, "%3.2f g\n", graymap[curcol]);

	fprintf(fp, "fill\n");

	if (!colour)
		fprintf(fp, "0 g\n");

	vdevice.cpVx = x[n - 1];
	vdevice.cpVy = y[n - 1];

	pslstx = pslsty = -1;		/* fill destroys current path */

	return(0);
}

#ifndef VOGLE
/*
 * Set the line width...
 */
static int
PS_setlw(w)
	int	w;
{
	fprintf(fp, "%d setlinewidth\n", w * 2 + 1);

	return(0);
}

/*
 * Set the line style...
 */
static int
PS_setls(lss)
	int	lss;
{
	unsigned ls = lss;
	int	i, d, a, b, offset;

	if (ls == 0xffff) {
		fprintf(fp, "[] 0 setdash\n");
		return;
	}

	fputc('[', fp);

	for (i = 0; i < 16; i++)	/* Over 16 bits */
		if ((ls & (1 << i)))
			break;

	offset = i;

#define	ON	1
#define	OFF	0
		
	a = b = OFF;
	if (ls & (1 << 0))
		a = b = ON;

	d = 0;
	for (i = 0; i < 16; i++) {	/* Over 16 bits */
		if (ls & (1 << i))
			a = ON;
		else
			a = OFF;

		if (a != b) {
			b = a;
			fprintf(fp, "%d ", d * 2 + 1);
			d = 0;
		}

		d++;
	}

	fprintf(fp, "] %d setdash\n", offset);

	return(0);
}

#else
/*
 * Set the line width...
 */
static int
PS_setlw(w)
	int	w;
{
	if (w == 0)
		w = 2;
	else if (w == 1)
		w = 4;

	fprintf(fp, "%d setlinewidth\n", w);

	return(0);
}
#endif

static DevEntry psdev = {
	"postscript",
	"large",
	"small",
	noop,
	PS_char,
	noop,
	PS_clear,
	PS_color,
	PS_draw,
	PS_exit,
	PS_fill,
	PS_font,
	noop,
	noop,
	PS_init,
	noop,
	PS_mapcolor,
#ifndef VOGLE
	PS_setls,
#endif
	PS_setlw,
	PS_string,
	noop,
	noop
};

/*
 * _CPS_devcpy
 *
 *	copy the postscript device into vdevice.dev.
 * 	Set it so we can use colours.
 */
int
_CPS_devcpy()
{
	vdevice.dev = psdev;
	vdevice.dev.devname = "cps";
	colour = 1;

	return(0);
}

/*
 * _PS_devcpy
 *
 *	copy the postscript device into vdevice.dev.
 */
int
_PS_devcpy()
{
	vdevice.dev = psdev;

	return(0);
}

/*
 * _PSP_devcpy
 *
 *	copy the postscript portrait device into vdevice.dev.
 */
int
_PSP_devcpy()
{
	vdevice.dev = psdev;
	vdevice.dev.Vinit = PSP_init;
	vdevice.dev.devname = "ppostscript";

	return(0);
}

/*
 * _PCPS_devcpy
 *
 *	copy the portrait postscript device into vdevice.dev.
 * 	Set it so we can use colours.
 */
int
_PCPS_devcpy()
{
	vdevice.dev = psdev;
	vdevice.dev.Vinit = PSP_init;
	vdevice.dev.devname = "pcps";
	colour = 1;

	return(0);
}


#ifdef LASERWRITER

/*
 * LASER_init
 *
 *	set up the postcript environment. Returns 1 on success.
 *	Opens pipe direct to lpr -Plw
 */
static int
LASER_init()
{
	fp = popen("lpr -Plw", "w");
	if (!fp) {
		fprintf(stderr, "Couldn't open pipe to lpr command.\n");
		exit(1);
	}

	page = 1;

	if (!ps_first_time)
		return(1);

	fputs("%!PS-Adobe-2.0 EPSF-1.2\n", fp);
	fputs("%%BoundingBox: 74 96 528 728\n", fp);
	fprintf(fp, "%%%%Page: %d %d\n", page, page);
	fputs("%%EndComments\n", fp);
	fprintf(fp, "72 300 div dup scale\n90 rotate\n400 -2200 translate\n");

	vdevice.sizeSy = 1890; 
	vdevice.sizeSx = 2634; 
	vdevice.sizeX = vdevice.sizeY = 1890; 

	PS_common_init();

	return (1);
}

/*
 * LASER_exit
 *
 *	do a showpage and close the output file if neccessary.
 */
static int
LASER_exit()
{
	fprintf(fp, "showpage\n");
	fprintf(fp, "%%Trailer\n");
	fflush(fp);

	pclose(fp);

	/*
	 * Bug fix from brett@kirk.es.go.dlr.de (Bernward Bretthauer)
	 */
	ps_first_time = 1;
	drawn = 0;
	curcol = 0;			/* black */
	pslstx = pslsty = -1;
	colour = 0;

	return(0);
}

/*
 * _LASER_devcpy
 *
 *	copy the postscript portrait device into vdevice.dev.
 */
int
_LASER_devcpy()
{
	vdevice.dev = psdev;
	vdevice.dev.Vinit = LASER_init;
	vdevice.dev.Vexit = LASER_exit;
	vdevice.dev.devname = "laser";
	colour = 1;	/* If you have a colour printer */

	return(0);
}
#endif

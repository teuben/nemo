/*
 * YAPP: Yet Another Plotting Package.  - Computer Graphics Metafile (CGM)
 *
 *      19-sep-89       Created for Pittsburgh UNICOS
 *	 9-dec-90       pl_matrix
 */

#include <stdinc.h>

extern int    debug_level;      /* see getparam.c */
extern string yapp_string;      /* see getparam.c */


local double dxymax;	            /* size of user window */
local double ax,bx,ay,by;           /* scaling factors user -> NDC */
local string def_device = "ps";     /* some default */

local int    ierr;                  /* error code to/from cgm routines */
local float xcur=0.0, ycur=0.0;     /* current location in VDC (0..1) */
local ncolors = 1;                  /* number of colors (1=B/W) */

float convx(), convy();             /* conversion functions */

/*
 * PLINIT: initalize the plotting package.
 */

plinit(pltdev, xmin, xmax, ymin, ymax)
char *pltdev;		/* output device name (ignored) */
double xmin, xmax;	/* user plotting area */
double ymin, ymax;
{
    float  width, height;
    char *device;

    if (pltdev==NULL || *pltdev==NULL || strlen(pltdev)==0) {
        dprintf(0,"YAPP_CGM: No plotdevice given, default %s used\n",
            def_device);
        device = def_device;
    } else if (strncmp(pltdev,"***")==0) {
        dprintf(0,"YAPP_CGM (***) : Default plotdevice %s used\n",
            def_device);
        device = def_device;
    }

    dprintf(1,"yapp_cgm: device = %s\n",yapp_string);

    csetdev(device,&ierr);              /* set driver */
    dprintf(1,"csetdev: ierr=%d\n",ierr);
    
    wrcopn("cgm.meta",&ierr);           /* open meta file */
    dprintf(1,"wrcopn: ierr=%d\n",ierr);

    wrmxci_(&ncolors,&ierr);             /* set number of colors */

    wrbegp_(&ierr);                     /* begin picture */

    wrbgpb_( &ierr);                     /* begin picture body */

    if (ymax - ymin < xmax - xmin) {
	dxymax = xmax - xmin;
	width = 1.0;
	height = (ymax - ymin) / dxymax;
    } else {
	dxymax = ymax - ymin;
	width = (xmax - xmin) / dxymax;
	height = 1.0;
    }
    
    ax = 1/(xmax-xmin);      bx = -1.0*ax;
    ay = 1/(ymax-ymin);      by = -1.0*ay;
}

/*
 * PLSWAP: does nothing.
 */

plswap() { }

/*
 * PLXSCALE, PLYSCALE: transform from usr to plotting coordinates.
 * At present, these do nothing (identity transformation).
 */

double plxscale(x, y)
double x, y;		/* user coordinates */
{
    return (x);
}

double plyscale(x, y)
double x, y;		/* user coordinates */
{
    return (y);
}

/*
 * PLLTYPE: select line width and dot-dash pattern.
 */

plltype(lwid, lpat)
int lwid;		/* line width */
int lpat;		/* line pattern */
{
    int lw, lp;

    if (lwid > 0) {
	lw = lwid;
    }
    if (lpat > 0) {
	lp = lpat-1;
    }
}

/*
 * PLLINE, PLMOVE, PLPOINT: plot lines, moves, and points.
 */

plline(x, y)
double x, y;		/* user coordinates */
{
    float xp[2],yp[2];
    int   np = 2;

    xp[0] = xcur;          yp[0] = ycur; 
    xcur = convx(x);        ycur = convy(y);
    xp[1] = xcur;          yp[1] = ycur; 
    wrplin_(xp,yp,&np,&ierr);               /* polyline */
    
}

plmove(x, y)
double x, y;		/* user coordinates */
{
   xcur = convx(x);        ycur = convy(y);
}

plpoint(x, y)
double x, y;		/* user coordinates */
{
    int n=0,istyle=0;		/* just a dot */
    float xp,yp;
    
    xp=x; yp=y;		/* RECALC !! */
}

/*
 * PLCIRCLE, PLCROSS, PLBOX: plot various symbols.
 */

plcircle(x, y, r)
double x, y;
double r;
{
    int npnts, i;
    double theta, sin(), cos();

    npnts = MAX(2400 * r / dxymax, 6.0);
    plmove(x + r, y);
    for (i = 1; i <= npnts; i++) {
	theta = TWO_PI * ((double) i) / ((double) npnts);
	plline(x + r * cos(theta), y + r * sin(theta));
    }
}

plcross(x, y, s)
double x, y;
double s;
{
    if (s > 0.0) {
	plmove(x - s, y);
	plline(x + s, y);
	plmove(x, y - s);
	plline(x, y + s);
    } else {
	s = s / 1.4142;
	plmove(x - s, y - s);
	plline(x + s, y + s);
	plmove(x - s, y + s);
	plline(x + s, y - s);
    }
}

plbox(x, y, s)
double x, y;
double s;
{
    if (s > 0.0) {
	plmove(x - s, y - s);
	plline(x + s, y - s);
	plline(x + s, y + s);
	plline(x - s, y + s);
	plline(x - s, y - s);
    } else {
	s = s * 1.4142;
	plmove(x - s, y);
	plline(x, y - s);
	plline(x + s, y);
	plline(x, y + s);
	plline(x - s, y);
    }
}

/*
 * PLJUST: specify justification of strings and numbers.
 * Imports: jus: 
 */

static int textjust = -1;

pljust(jus)
int jus;		/* -1, 0, 1 for left, mid, right just */
{
    textjust = (jus < -1 ? -1 : (jus > 1 ? 1 : jus));
}

/*
 * PLTEXT: plot a text string.
 */

pltext(msg, x, y, hgt, ang)
char *msg;		/* message to plot, with NULL termination */
double x, y;		/* user coordinates (modified by justification) */
double hgt;		/* height of characters in user coordinates */
double ang;		/* angle of message, counterclockwise in degrees */
{
    double c, cos(), s, sin();
    float dx, dy, xp, yp, ap;
    float newsize, sl, sh;
    int   n;

#if 1
	if (1)
		return;
#endif
    xp=x; yp=y; ap=ang;
    c = cos(ang / 57.296);
    s = sin(ang / 57.296);

}

/*
 * PLFLUSH: output any pending graphics.
 */

plflush() { }

/*
 * PLFRAME: advance to next frame.
 */

plframe()
{
    wrendp_(&ierr);     /* end picture */
    wrbegp_(&ierr);     /* begin picture */
    wrbgpb_(&ierr);     /* begin picture body */
}

/*
 * PLSTOP: finish up a display.
 */

#define ONEMIN (60 * 1000 * 1000)

plstop()
{
    wrendp_(&ierr);     /* end picture */
    wrtend_(&ierr);     /* end metafile */
}

pl_matrix(frame,nx,ny,xmin,ymin,cell,fmin,fmax,findex)
double *frame, xmin, ymin, cell, fmin, fmax, findex;
int nx, ny;
{
    double x,y,f,grayscale,ds,pow();
    int ix,iy;

#if 1
    dprintf(0,"PL_MATRIX is unsupported in cgm version of YAPP\n");
    return(0);
#endif
}


local float convx(x)
double x;
{
    return( (float) (ax*x+bx) );
}

local float convy(y)
double y;
{
    return( (float) (ay*y+by) );
}


#ifdef TESTBED

int debug_level = 1;
char *yapp_string = "";

main(argc, argv)
int argc;
string argv[];
{
    int i, j;


    plinit(argv[1], 0.0, 20.0, 0.0, 20.0);
    plmove(0.0, 0.0);
    plline(20.0, 0.0);
    plline(20.0, 20.0);
    plline(0.0, 20.0);
    plline(0.0, 0.0);
    plline(20.0, 20.0);
    plmove(20.0, 0.0);
    plline(0.0, 20.0);
    plltype(12, 0);
    plmove(4.0, 18.0);
    plline(16.0, 18.0);
    plltype(-6, 0);
    plmove(6.0, 18.0);
    plline(14.0, 18.0);
    for (i = 1; i <= 4; i++) {
	plltype(i, 1);
        plmove(1.0, 13.0 - i);
        plline(3.0, 13.0 - i);
        plpoint(3.5, 13.0 - i);
	plltype(1, i);
	for (j = 1; j <= 4; j++) {
	    plmove(1.5, 13.0 - i - 0.2*j);
	    plline(1.5 + j, 13.0 - i - 0.2*j);
	}
    }
    plltype(1, 1);
    plcircle(15.0, 9.0, 0.5);
    plcircle(16.0, 9.0, 0.25);
    plcircle(17.0, 9.0, 0.125);
    plcircle(18.0, 9.0, 0.0625);
    plbox(16.0, 8.0, 0.4);
    plbox(17.0, 8.0, 0.2);
    plbox(18.0, 8.0, -0.2);
    plcross(16.0, 7.0, 0.4);
    plcross(17.0, 7.0, 0.2);
    plcross(18.0, 7.0, -0.2);
    pltext("Foo Bar!", 8.0, 5.0, 0.5, 0.0);
    pltext("Fum Bar!", 8.0, 3.0, 0.25, 0.0);
    for (i = 0; i <= 4; i++)
	pltext(" testing angles", 16.0, 10.0, 0.2, 45.0*i);
    plmove(10.0, 8.5);
    plline(10.0, 11.5);
    pljust(-1);
    pltext("left justified",  10.0,  9.0, 0.25, 0.0);
    pljust(0);
    pltext("centered",        10.0, 10.0, 0.25, 0.0);
    pljust(1);
    pltext("right justified", 10.0, 11.0, 0.25, 0.0);
    plstop();
}

#endif

pl_screendump(fname)
char *fname;
{
  dprintf(0,"pl_screendump(%s): Not implemented for yapp_cgm\n",fname);
}

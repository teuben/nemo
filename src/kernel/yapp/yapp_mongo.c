/*
 * YAPP: Yet Another Plotting Package.   :   mongo-87
 *
 *	xx-Mar-88:	mongo87 interface version - Peter Teuben
 *	 1-feb-89:      yapp==0 is now an error
 *	17-apr-89	yapp==0 will run without graphics output though	
 *	26-oct-90       fixed some bugs:
 *                      corrected some 'real' to 'float' in mongopar_
 *			I cannot imagine this ever worked.....
 *                      also calling mgogstring_() was wrong.....
 *
 *		Normally link cc program with:
 *	-lmongo -lsuntool -lsunwindow -lpixrect -lF77 -lI77 -lU77 -lm
 *		some of the '77' libs might be not necessary in SUNOS >=4.1
 *	-lmongo can also be replaced with $(MONGOFILES)/libmongo.a
 *		in case this library is not in search path
 *
 * This yapp_mongo links fortran and C and hence, as written, will only work
 * on machines which conform to the BSD convention. The /mongopar/ common block 
 * defined in C as a structure, also likely to be machine dependant.
 */

#include <stdinc.h>

extern int   yapp_dev;	/* see getparam.c */
extern struct {		/* The /COMMON/MONGOPAR from fortran */
	float x1, x2, y1, y2;   /* User coordinates of plot region */
	float gx1,gx2,gy1,gy2;	/* Physical coordinates of user plot region */
				/* GX1<GX2,GY1<GY2 */
 	int   lx1,lx2,ly1,ly2;	/* Limit physical coordinates of device */
	float xp,yp;		/* Present plot position (phys.coordinates) */
	float cxp,cyp;		/* Present plot position (user coordinates)  */
	float expand; 		/* expansion factor of characters */
	float angle;		/* Rotation of characters (etc) */
				/*	 (counterclockwise) */
	int   ltype;		/* type of output lines */
	int   lweight;		/* weight: 0,1 for single, 2 for double, etc */
	float cheight;		/* Expand 1.0 height of characters */
				/*		(physical units) */
	float cwidth;		/* Expand 1.0 width of characters */
				/*		(physical units) */
	float cxdef,cydef;	/* Expand 1.0 expansion of characters drawn */
				/*		by vector (X dir) */
	float pxdef,pydef;	/* Expand 1.0 radius of points, (X dir) */
	float coff;		/* Terminal vertical offset of characters */
				/*			(physical) */
        int   termout;		/* True for a terminal, false for hardcopy */
	int   xyswapped;	/* True for an inverted picture x <-> y */
	int   numdev;		/* Device number */
        float pi;		/* 3.1415... */
	float uservar[10];	/* Array of user set-able variables */
	int   autodot;		/* True if device has automatic dotted/dashed */
				/*		line generation */
} mongopar_;

local real   dxymax;	/* size of user window */
local int    iterm;	/* terminal number */

/*
 * PLINIT: initalize the plotting package.
 */

plinit(pltdev, xmin, xmax, ymin, ymax)
char *pltdev;		/* output device name (ignored) */
real xmin, xmax;	/* user plotting area */
real ymin, ymax;
{
    float width, height, x1,x2,y1,y2;
    bool  Qsquare;

    if (pltdev==NULL || *pltdev==NULL)
        iterm = 0;
    else
        iterm = atoi(pltdev);

    if (iterm==0)		/* take default from environment */
        iterm = yapp_dev;       /* or as initialized from commandline 'yapp=' */
    if (iterm==0) {
        warning("No YAPP device specified, use 'setenv YAPP #' or 'yapp=#'");
        warning("All graphics output will be suppressed");
	dprintf(0,"[This is yapp_mongo]\n");  
        dprintf(0,"N       YAPP                     N      YAPP\n\n");
	dprintf(0,"1       Retrographics 640       -1      Versatec (vertical)\n");
        dprintf(0,"2       DEC VT125               -2      Versatec (horizontal)\n");
        dprintf(0,"3       Tektronix 4010          -3      Printronix (vertical)\n");
        dprintf(0,"4       Grinnell 270            -4      Printronix (horizontal)\n");
        dprintf(0,"5       HP2648A                 -5      Laser (vertical)\n");
        dprintf(0,"6       SUN windows             -6      Laser (horizontal)\n");
	dprintf(0,"                                -7      Laser (square)\n");							
	return;
    } else
        dprintf(1,"yapp_mongo: term (yapp) = %d\n",iterm);

    mgoinit_();
    mgosetup_(&iterm);
    mgoerase_();

    if (ymax - ymin < xmax - xmin) {
	dxymax = xmax - xmin;
	width = 1.0;
	height = (ymax - ymin) / dxymax;
    } else {
	dxymax = ymax - ymin;
	width = (xmax - xmin) / dxymax;
	height = 1.0;
    }
    
    x1=mongopar_.lx1; 	/* use the g__ or l__ ?? */
    y1=mongopar_.ly1;
    x2=mongopar_.lx2;
    y2=mongopar_.ly2;

    x1 = 0.0;
    y1 = 0.0;
    x2 = MIN(x2,y2);    /* make square; ignore xmin, xmax etc for now */
    y2 = x2;

    mgosetloc_(&x1,&y1,&x2,&y2); /* use full page (l__) or partial (g__)?*/

    x1=xmin; y1=ymin;
    x2=xmax; y2=ymax;
    mgosetlim_(&x1,&y1,&x2,&y2);	/* set user coordinates as 'cm' */


    dprintf (2,"Physical dev: %d %d %d %d\n",mongopar_.lx1,mongopar_.ly1,
        mongopar_.lx2,mongopar_.ly2);
    dprintf (2,"Physical plt: %f %f %f %f\n",mongopar_.gx1,mongopar_.gy1,
        mongopar_.gx2,mongopar_.gy2);
    dprintf (2,"    User plt: %f %f %f %f\n",mongopar_.x1,mongopar_.y1,
        mongopar_.x2,mongopar_.y2);

}

/*
 * PLSWAP: does nothing.
 */

plswap() { }

/*
 * PLXSCALE, PLYSCALE: transform from usr to plotting coordinates.
 * At present, these do nothing (identity transformation).
 */

real plxscale(x, y)
real x, y;		/* user coordinates */
{
    if (iterm==0) return;       /* no graphics output requested */
    return (x);
}

real plyscale(x, y)
real x, y;		/* user coordinates */
{
    if (iterm==0) return;       /* no graphics output requested */
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

    if (iterm==0) return;       /* no graphics output requested */
    if (lwid > 0) {
	lw = lwid;
	mgosetlweight_(&lw);
    }
    if (lpat > 0) {
	lp = lpat-1;
	mgosetltype_(&lp);
    }
}

/*
 * PLLINE, PLMOVE, PLPOINT: plot lines, moves, and points.
 */

plline(x, y)
real x, y;		/* user coordinates */
{
    float xp,yp;

    if (iterm==0) return;       /* no graphics output requested */

    xp=x; yp=y;		/* RECALC !! */
    mgodraw_(&xp, &yp);
}

plmove(x, y)
real x, y;		/* user coordinates */
{
   float xp,yp;

   if (iterm==0) return;       /* no graphics output requested */

   xp=x; yp=y;		/* RECALC !! */
   mgorelocate_(&xp, &yp);
}

plpoint(x, y)
real x, y;		/* user coordinates */
{
    int n=0,istyle=0;		/* just a dot */
    float xp,yp;
    
    if (iterm==0) return;       /* no graphics output requested */

    xp=x; yp=y;		/* RECALC !! */
    mgorelocate_(&xp, &yp);
    mgopoint_(&n, &istyle);
}

/*
 * PLCIRCLE, PLCROSS, PLBOX: plot various symbols.
 */

plcircle(x, y, r)
real x, y;
real r;
{
    int npnts, i;
    double sin(), cos();
    real theta;

    if (iterm==0) return;       /* no graphics output requested */

    npnts = MAX(2400 * r / dxymax, 6.0);
    plmove(x + r, y);
    for (i = 1; i <= npnts; i++) {
	theta = TWO_PI * ((real) i) / ((real) npnts);
	plline(x + r * cos(theta), y + r * sin(theta));
    }
}

plcross(x, y, s)
real x, y;
real s;
{
    if (iterm==0) return;       /* no graphics output requested */

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
real x, y;
real s;
{
    if (iterm==0) return;       /* no graphics output requested */

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
    if (iterm==0) return;       /* no graphics output requested */

    textjust = (jus < -1 ? -1 : (jus > 1 ? 1 : jus));
}

/*
 * PLTEXT: plot a text string.
 */

pltext(msg, x, y, hgt, ang)
char *msg;		/* message to plot, with NULL termination */
real x, y;		/* user coordinates (modified by justification) */
real hgt;		/* height of characters in user coordinates */
real ang;		/* angle of message, counterclockwise in degrees */
{
    float dx, dy, xp, yp, ap;
    float sl, sh, newsize;
    int   n, loc;

    if (iterm==0) return;       /* no graphics output requested */

    xp=x; yp=y; 
    ap=ang;
    mgosetangle_(&ap);			/* set angle */
    newsize = (3 * hgt);		/* new expansion size */
    mgosetexpand_(&newsize);
    mgorelocate_(&xp,&yp);              /* go to this position */
    loc = -textjust + 5;		/* mongo uses loc */
    n = strlen(msg);                    /* how long is string */
    mgoputlabel_(&n,msg,&loc,n);        /* put the label with the righ 'loc' */
    ap=0.0; mgosetangle_(&ap);		/* set back to default ? */
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
    int nvec;

    if (iterm==0) return;       /* no graphics output requested */

    if (iterm<0) {
	mgoprntplot_(&nvec);
	dprintf (0,"%d vectors plotted\n",nvec);
    } else {
	mgotidle_();
   }
	
   mgoerase_();
}

/*
 * PLSTOP: finish up a display.
 */

#define ONEMIN (60 * 1000 * 1000)

plstop()
{
    int nvec, bell=0x07;

    if (iterm==0) return;       /* no graphics output requested */

    if (iterm<0) {
	mgoprntplot_(&nvec);
	dprintf (0,"%d vectors plotted\n",nvec);
    } else {
	mgotidle_();
	putchar(bell); fflush(stdout);		/* wait for keyboard */
	getchar();
   }

}

pl_matrix(frame,nx,ny,xmin,ymin,cell,fmin,fmax,findex)
real *frame, xmin, ymin, cell, fmin, fmax, findex;
int nx, ny;
{
    real x,y,f,grayscale,ds,pow();
    int ix,iy;
    
    if (iterm==0) return;       /* no graphics output requested */

    warning("pl_matrix: Not implemented in yapp_mongo");
    return(0);
}

pl_screendump(fname)
string fname;
{
  warning("pl_screendump(%s): Not implemented in yapp_mongo",fname);
}

pl_getpoly(x,y,n)
float x[], y[];
int n;
{
  warning("pl_getpoly: Not implemented in yapp_mongo");
  return(0);
}


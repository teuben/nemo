/*
 * YAPP: Yet Another Plotting Package.
 *	
 *      YAPP_X11:    YAPP for X11:
 *		Also designed for movies within Xwindows
 *		Has support through pl_screendump
 *
 *		 7-mar-90	Peter Teuben
 *		20-jan-94  comments on portability in code
 *
 */
#include <stdinc.h>
#include <getparam.h>
#include <ctype.h>


#include <X11/Xlib.h>		/* X11 stuff */
#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#include <X11/Intrinsic.h>
#ifdef linux
#include <X11/Xaw/Viewport.h>	/* <====== see also <X11/Xaw> */
#else
#include <X11/Viewport.h>	/* <====== see also <X11/Xaw> */
#endif
#include <X11/StringDefs.h>

#ifndef linux
#include "Canvas.h"
#endif

/*   #include "mongo.h"	*/

typedef struct {
   unsigned char red,
                 green,
                 blue;
} COLOR;
				     


#define MAX_POLYGON 10		/* maximum corners for PTYPE 3 points */
#define START_NVEC 200		/* starting value for max_nvec */

static int Wwidth, Wheight;             /* size of window */
static double ax, bx, ay, by;           /* scale factors */
static double xp, yp;                   /* current pen position */
static double width, height, dxymax;    /* window stuff */
static int nvec,			/* number of vectors drawn */
	   max_nvec;			/* max number of vectors allocated */

static Widget toplevel = 0,
	      frame,
	      canvas;

static Arg arglist[] = {
    {XtNwidth, (XtArgVal) 512},
    {XtNheight, (XtArgVal) 512},
};

static GC graphics, erasegc, redrawgc;

static XSegment *xvec;

static Cursor curs;			/* graphics cursor to use */

static Pixmap backing;

static Colormap cmap;

static XColor xcolor;


/* local functions: */
static x11_redraw(), size_11window();

char *malloc();

plinit(pltdev, xmin, xmax, ymin, ymax)
char *pltdev;					/*  s is old name in x11 */
double xmin, xmax, ymin, ymax;
{
    int argc;
    char *argv[2], *s, *argv0;
    XGCValues gcvalues, evalues, revalues;
    
    if (strncmp(pltdev,"***",3)==0)
        s = NULL;               /* no device speified */
    else
        s = pltdev;
    argv0 = getargv0();
    argv[0] = argv0;
    argv[1] = s;
    if(toplevel == 0) {
	if(s != NULL) {
	    while(*s != '\0' && isspace(*s)) s++;   /* skip whitestuff */
	    if(s[0] == '\0') {
		s = NULL;
		argc = 1;
	    } else argc = 2;
	} else argc = 1;
	toplevel = XtInitialize(argv0, argv0, NULL, 0, &argc, argv);
	frame = XtCreateManagedWidget(argv0, viewportWidgetClass, toplevel,
				    arglist, XtNumber(arglist));
	canvas = XtCreateManagedWidget("canvas", canvasWidgetClass, frame,
				       arglist, XtNumber(arglist));
	XtRealizeWidget(toplevel);
	backing = XCreatePixmap(XtDisplay(canvas), XtWindow(canvas), 512, 512,
			    DefaultDepth(XtDisplay(canvas),0));
	gcvalues.function = GXset;
	gcvalues.fill_style = FillSolid;
	gcvalues.subwindow_mode = IncludeInferiors;
	graphics = XCreateGC(XtDisplay(canvas), XtWindow(canvas),
			     GCFunction|GCFillStyle|GCSubwindowMode,
			     &gcvalues);
	evalues.function = GXclear;
	evalues.subwindow_mode = IncludeInferiors;
	erasegc = XCreateGC(XtDisplay(canvas), XtWindow(canvas),
			     GCFunction|GCSubwindowMode,
			     &evalues);
	revalues.function = GXcopy;
	revalues.subwindow_mode = IncludeInferiors;
	redrawgc = XCreateGC(XtDisplay(canvas), XtWindow(canvas),
			     GCFunction|GCSubwindowMode,
			     &revalues);
	curs = XCreateFontCursor(XtDisplay(canvas), XC_crosshair);
	max_nvec = START_NVEC;
        if((xvec = (XSegment *)malloc((unsigned)max_nvec*
					 sizeof(XSegment))) == NULL) {
	    fprintf(stderr,"Can't allocate vectors in x_setup\n");
	    max_nvec = 0;
	    return(-1);
	}
	nvec = 0;
        cmap = DefaultColormap(XtDisplay(canvas), 
                               DefaultScreen(XtDisplay(canvas)));
    }
    size_11window();
    if (ymax - ymin < xmax - xmin) {        /* make a square area */
	dxymax = xmax - xmin;
	width = 1.0;
	height = (ymax - ymin) / dxymax;
    } else {
	dxymax = ymax - ymin;
	width = (xmax - xmin) / dxymax;
	height = 1.0;
    }
    ax = (Wwidth) / (xmax-xmin);
    bx = -ax*xmin;
    ay = (Wheight) / (ymax-ymin);
    by = -ay*ymin;
    x11_redraw();

#if 0
/* finally, set some variables for SM */
   default_ctype("black");		/* is this correct? */
   ldef = 32;
#endif

}
plswap() {}

double plxscale(x,y)
double x,y;
{
    return(x);
}

double plyscale(x,y)
double x,y;
{
    return(y);
}

pl_line(xf, yf, xt, yt)         /* NON-STANDARD */
double xf, yf, xt, yt;          /*  from .. to */
{
    int x1,x2,y1,y2;

    xp = xt;                    /* save end point */
    yp = yt;
    x1 = (int) (ax*xf+bx);
    x2 = (int) (ax*xt+bx);    
    y1 = Wheight - (int) (ay*yf+by);
    y2 = Wheight - (int) (ay*yt+by);
    xvec[nvec].x1 = x1;
    xvec[nvec].y1 = y1;
    xvec[nvec].x2 = x2;
    xvec[nvec].y2 = y2;
    if(++nvec >= max_nvec) {
	XDrawSegments(XtDisplay(canvas), XtWindow(canvas), graphics, xvec,
		      nvec);
	XDrawSegments(XtDisplay(canvas), backing, graphics, xvec,
		      nvec);
	nvec = 0;
    }
}

/*
 * PLLTYPE: select line width and dot-dash pattern.
 */

plltype(lwid, lpat)
int lwid;		/* line width */
int lpat;		/* line pattern */
{
}

plmove(x, y)
double x, y;   
{
    xp = x;
    yp = y;
}

plline(x, y)
double x, y;
{
    pl_line(xp, yp, x, y);
}

plpoint(x, y)
double x, y;		/* user coordinates */
{
    printf("plpoint: %g %g no implemented\n",x,y);
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
 */

static int textjust = -1;

pljust(jus)
int jus;		/* -1, 0, 1 for left, mid, right just */
{
    textjust = (jus < -1 ? -1 : (jus > 1 ? 1 : jus));
}

pltext(s, xt, yt, hgt, ang)
char *s;
double xt, yt;
double hgt, ang;            /* ignored for now */
{
    int x,y;

    if(s == NULL) return(0);		/* we do hard characters */

    x = (int) (ax*xt+bx);
    y = Wheight - (int) (ay*yt+by);
    XDrawString(XtDisplay(canvas),XtWindow(canvas),graphics,x,y,s,strlen(s));
    XDrawString(XtDisplay(canvas),backing,graphics,x,y,s,strlen(s));
}

pl_erase()          /* NON-STANDARD */
{
    Dimension width;
    Dimension height;
    Arg arglist[2];

    XtSetArg(arglist[0], XtNwidth, (XtArgVal) &width);
    XtSetArg(arglist[1], XtNheight, (XtArgVal) &height);
    XtGetValues(canvas, arglist, XtNumber(arglist));
    XFillRectangle(XtDisplay(canvas), XtWindow(canvas), erasegc, 0, 0, width,
		   height);
    XFillRectangle(XtDisplay(canvas), backing, erasegc, 0, 0, width,
		   height);
    nvec = 0;
}

/*  Dummy for x11  */
pl_set_ctype(colors, ncolors)
COLOR *colors;
int ncolors;
{
      return(-1);
}

/*
 * Set colour of line 
 */
pl_ctype(r,g,b)
int r,g,b;
{
   xcolor.red = 65535*r/255;
   xcolor.green = 65535*g/255;
   xcolor.blue = 65535*b/255;
   xcolor.flags = DoRed & DoGreen & DoBlue;
   if (XAllocColor(XtDisplay(canvas), cmap, &xcolor) == 0) {
      dprintf(0,"couldn't allocate specified color");
      return(-1);
    }
    XSetForeground(XtDisplay(canvas), graphics, xcolor.pixel);
    return(0);
}

pl_idle()           /* non-standard YAPP */
{
   if(XtPending() > 0) {
      size_11window();
      x11_redraw();
   }
}

plflush()
{
    if(nvec) {
	XDrawSegments(XtDisplay(canvas), XtWindow(canvas), graphics, xvec,
		      nvec);
	XDrawSegments(XtDisplay(canvas), backing, graphics, xvec,
		      nvec);
	nvec = 0;
    }
}

#if 0
pl_fill_pt(n)           /* non-standard YAPP */
int n;                          /* number of sides */
{
   float dtheta, theta;
   static float old_ang;	/* old values of angle */
   int i,
       xpsize,			/* scale for points == g_dx*PDEF*eexpand */
       ypsize;
   static int num = -1,         /* number of vertices used last time */
	      old_xp,old_yp,	/* old values of xpsize, ypsize */
	      x0,y0;		/* constant part of vertex[0].x, .y */
   static XPoint vlist[MAX_POLYGON + 1];  /* vertices describing point */

   if(n < 2) {
      x11_line(xp,yp,xp,yp);
      return;
   }

   dtheta = 2*PI/n;
   xpsize = 2*PDEF*sin(dtheta/2)*eexpand*g_dx;
   ypsize = 2*PDEF*sin(dtheta/2)*eexpand*g_dy;
   if(n != num || aangle != old_ang || xpsize != old_xp || ypsize != old_yp) {
      if(n > MAX_POLYGON) num = MAX_POLYGON;
      else num = n;

      theta = 3*PI/2 + dtheta/2 + aangle*PI/180;

      old_ang = aangle;
      old_xp = xpsize;
      old_yp = ypsize;

      x0 = PDEF*g_dx*eexpand*cos(theta);
      y0 = height - PDEF*g_dy*eexpand*sin(theta);

      theta += dtheta/2;
      for(i = 1;i <= num + 1;i++) {
	 vlist[i].x = -xpsize*sin(theta);
	 vlist[i].y = -ypsize*cos(theta);	/* screen is upside down */
	 theta += dtheta;
      }
   }
   vlist[0].x = g_dx*xp + x0;
   vlist[0].y = y0 - g_dy*yp;
   XFillPolygon(XtDisplay(canvas), XtWindow(canvas), graphics, vlist, num+1,
		Convex, CoordModePrevious);
   XFillPolygon(XtDisplay(canvas), backing, graphics, vlist, num+1,
		Convex, CoordModePrevious);
}
#endif

#if 0
int
pl_cursor(x, y)     /* non-standard YAPP */
int *x, *y;
{
   XEvent rep;
   char buf[5];

   XDefineCursor(XtDisplay(canvas), XtWindow(canvas), curs);
   						/* define graphics cursor */
   XSelectInput(XtDisplay(canvas), XtWindow(canvas), ExposureMask
		| KeyPressMask);

   for(rep.type = 0; rep.type != KeyPress;)
       XtNextEvent(&rep);		/* sense cursor button */

   XSelectInput(XtDisplay(canvas), XtWindow(canvas), ExposureMask);
   *x = rep.xkey.x/g_dx;
   *y = (height - rep.xkey.y)/g_dy;
   XLookupString(&rep, buf, sizeof(buf), 0, (XComposeStatus *) NULL);
   return((int) buf[0]);
}
#endif

plstop()
{
   ;
}

static
x11_redraw()
{
    XEvent rep;

    while(XtPending()) {
	XtNextEvent(&rep);
    }
    XCopyArea(XtDisplay(canvas), backing, XtWindow(canvas), redrawgc, 0, 0,
	      512, 512, 0, 0);
}

static
size_11window()
{
    
    Dimension height;		/* height of window */
    Dimension width;            /* width of window */
    Arg arglist[2];

    XtSetArg(arglist[0], XtNwidth, (XtArgVal) &width);
    XtSetArg(arglist[1], XtNheight, (XtArgVal) &height);
    XtGetValues(canvas, arglist, XtNumber(arglist));
    Wheight = height;               /* save the global parameter for yapp */
    Wwidth = width;
}

#if defined(TESTBED)

main(argc, argv)
int argc;
string argv[];
{
    int i, j;
    char *getenv();

    plinit(getenv("DISPLAY"), 0.0, 20.0, 0.0, 20.0);
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


/*
 * YAPP: Yet Another Plotting Package.
 *	
 *      YAPP_SV:    YAPP for SunView
 *
 *	Specifically designed for movies within sunwindows, does not
 *	plot text, uses pl_screendump
 *
 *		Has support through pl_screendump
 *		|yapp| is size of screen in pixels
 *		yapp>0:  at plstop() immediately quits
 *		yapp<0:  waits for mouse press
 *
 *		 9-dec-88	written - Peter Teuben
 *		14-mar-90	trouble with GCC ...
 *		22-jun-90	appease GCC - but still does not compile
 *		11-jul-91       extra warning when |yapp|<50
 */

#include <stdinc.h>

#include <suntool/sunview.h>
#include <suntool/canvas.h>

local Frame  base_frame;              /* base frame structure  */
local Canvas canvas;                  /* the canvas itself */
local Pixwin *pw;                     /* pixwin for screen handling */
local struct pixrect *pr;             /* memory pixrect for screendumps */


local short icon_image[] = {
#include <yapp.icon>
};

DEFINE_ICON_FROM_IMAGE(yapp_icon, icon_image);

local double dxymax;		/* size of user window */
local double ax,bx,ay,by;   	/* scaling between NDC and WC */

local int ix, iy;               /* current CP */
local int npx, npy;             /* number of pixels in X and Y */
local int pcolor;               /* color of pixels to plot */
   
local int interact;             /* 0: no text, no interaction (movies) */

local int plix(), pliy(), next();

/*
 * PLINIT: initalize the plotting package.
 */

plinit(pltdev, xmin, xmax, ymin, ymax)
char *pltdev;		/* output device name (ignored in this yapp) */
double xmin, xmax;	/* user plotting area */
double ymin, ymax;
{
    extern int yapp_dev;
    float  width, height;
    int    x0,x1,y0,y1, xsize, ysize;

    dprintf (1,"YAPP_SV: pltdev=%s\n",pltdev);

    if (streq(pltdev,"***") || yapp_dev==0 || strlen(pltdev)==0)
        xsize = ysize = 256;		/* default size */
    else if (strlen(pltdev) > 0)
        yapp_dev = atoi(pltdev);

    if (yapp_dev != 0)
        xsize = ysize = ABS(yapp_dev);
    interact = (yapp_dev < 0);  /* allow interaction if yapp<0 */
    if (xsize<50 || ysize<50) 
        warning("YAPP_SV has a really small window of %d X %d pixels",
            xsize, ysize);

    dprintf (1,"              size = %d  interact=%d\n",xsize,interact);

    base_frame = window_create(NULL, FRAME,
                        FRAME_SHOW_LABEL, FALSE,
                        FRAME_ICON,    &yapp_icon,
                        WIN_ERROR_MSG, "Can't create window",
                        0);
    window_set (base_frame,
                        WIN_SHOW,    TRUE,      /* show it anyhow */
                        WIN_WIDTH,   xsize+10,  /* size, including border */
                        WIN_HEIGHT,  ysize+10,  /* which seems to be 10 pixels */
                        0);
                        
    canvas = window_create(base_frame, CANVAS,
                        0);
                        
    pw = canvas_pixwin(canvas);
    
    npx = (int) window_get(canvas, CANVAS_WIDTH);   /* pixels in X */
    npy = (int) window_get(canvas, CANVAS_HEIGHT);  /* pixels in Y */
    printf ("Window is %d by %d pixels\n",npx,npy);
    
    pr =  mem_create(npx,npy,1);        /* create 1 bit deep mem pixrect */
                                        /* used in pl_screendump() */   
        
    ax = (npx - 1.0)/(xmax-xmin);
    bx = -xmin * ax;
    ay = (npy - 1.0)/(ymax-ymin);
    by = -ymin * ay;        
    
    ix = iy = 0;

    if (ymax - ymin < xmax - xmin) {        /* make a square area */
	dxymax = xmax - xmin;
	width = 1.0;
	height = (ymax - ymin) / dxymax;
    } else {
	dxymax = ymax - ymin;
	width = (xmax - xmin) / dxymax;
	height = 1.0;
    }

   (void) notify_do_dispatch();     /* implicit dispatching (see manual) */
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
}

/*
 * PLLINE, PLMOVE, PLPOINT: plot lines, moves, and points.
 */

plline(x, y)
double x, y;		/* user coordinates */
{
    int ixx,iyy;
    
    ixx = plix(x);
    iyy = pliy(y);
    
    pw_vector(pw, ix,iy,ixx,iyy, PIX_SRC, 1);
    ix = ixx;
    iy = iyy;
    
}

plmove(x, y)
double x, y;		/* user coordinates */
{
    ix = plix(x);
    iy = pliy(y);

}

plpoint(x, y)
double x, y;		/* user coordinates */
{
    ix = plix(x);
    iy = pliy(y);

    pw_put(pw, ix, iy, 1);
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

/*
 * PLTEXT: plot a text string.
 */

pltext(msg, x, y, hgt, ang)
char *msg;		/* message to plot, with NULL termination */
double x, y;		/* user coordinates (modified by justification) */
double hgt;		/* height of characters in user coordinates */
double ang;		/* angle of message, counterclockwise in degrees */
{
    if (interact) {
       ix = plix(x);
       iy = pliy(y);
       pw_text(pw, ix, iy, PIX_SRC | PIX_DST, 0, msg);
    }
}
 
local int plix(x)
double x;
{
    int i;

    i = ax*x+bx;
    if (i<0) i=0;
    if (i>=npx) i=npx-1;
    return(i);
}

local int pliy(y)
double y;
{
    int i;

    i = ay*y+by;

    if (i<0) i=0;
    if (i>=npy) i=npy-1;
    return(npy-1-i);
}


/*
 * PLFLUSH: output any pending graphics.
 */

plflush() 
{
   window_set (base_frame,
                        WIN_SHOW,    TRUE,      /* show it anyhow */
                        0);
}

/*
 * PLFRAME: advance to next frame, with optional pause??
 */

plframe()
{
   int my_done, n;

   next();
   pw_writebackground(pw,0,0,npx,npy,PIX_SRC);      /* clear old stuff */
}

/*
 * PLSTOP: finish up a display.
 */

plstop()
{
    next();
}


local int next()
{
   int c;

   if (interact) {
     printf ("PLSTOP: enter any key to advance to next plot or quit: \n");
     fflush(stdout);
     fflush(stdin);
     read(0,&c,1);
   }
}
    

pl_matrix(frame,nx,ny,xmin,ymin,cell,fmin,fmax,findex)
double *frame, xmin, ymin, cell, fmin, fmax, findex;
int nx, ny;
{
    dprintf(0,"PLMATRIX is unsupported in YAPP_SV\n");
    return(0);
}

pl_screendump(fname)
char *fname;
{
    FILE *fscrdmp, *fopen();
    int err;
    
    printf ("PL_SCREENDUMP: %s\n",fname);
    plflush();                                  /* flush just to be sure */
    
    if ((fscrdmp = fopen(fname,"w"))==NULL)
        error("pl_screendump: Cannot create dump file %s [yapp_sv]\n",
                fname);
        
    pw_read (pr,0,0,npx,npy,PIX_SRC,pw,0,0);    /* put pixwin -> pixrect */
    
    if (pr_dump (pr, fscrdmp, NULL, RT_STANDARD, TRUE)) /* save in file */
        error("pl_screendump: error dumping rasterfile %s [yapp_sv]\n",
                fname);
    
    fclose(fscrdmp);
}

pl_getpoly(x,y,n)
float *x, *y;
int n;
{
  dprintf(0,"Warning: pl_getpoly not implemented\n");
  return(0);
}


/*
 * YAPP: Yet Another Plotting Package.
 * Joshua Barnes  Dec 1986  I. A. S.  Princeton, NJ.
 *
 *    29-jun-87  added pl_matrix() for yapp_ps consistency   PJT
 *    13-aug-88	 added code to autodetect screen mode       PJT
 *               added mouse functions for polygons         PJT
 *    16-dec-88  added pl_screendump() function             PJT
 *
 * This version works with the SunCore graphics library.
 *	*** in new SUN OS 4.x "usercore.h" is not in /usr/include;
 *	    anymore; we have put an old copy in $NEMO/inc
 *	    in SUN OS 4.1.2 the libcore.a has disappeared from /usr/lib
 *	    and that will be the end of this yapp...
 */

#include <stdinc.h>
#include <usercore.h>

/*			Possible compiler defines to expand functionality
 * #define MOUSE
 * #define COLOR
 * #define RETAINED
 */


#if defined(COLOR)
  int cgpixwindd();	/* Color window */
  local struct vwsurf winsurf = DEFAULT_VWSURF(cgpixwindd);
  int cg2dd();          /* Color full screen */
  local struct vwsurf rawsurf = DEFAULT_VWSURF(cg2dd);
#else
  int pixwindd();	/* BW window */
  local struct vwsurf winsurf = DEFAULT_VWSURF(pixwindd);
  int bw2dd();          /* BW full screen */
  local struct vwsurf rawsurf = DEFAULT_VWSURF(bw2dd);
#endif

local double dxymax;		/* size of user window */
local double ax,bx,ay,by;   	/* scaling between NDC and WC */
local struct vwsurf *surface;	/* this will be our surface */
#if defined(RETAINED)
local int keep_old_segments=0;	/* 1=TRUE, keep'm  0=delete'm */
local int segment_name=0;	/* name of segment, 0=unused */
#endif

local struct vwsurf *get_surface(); 
local bell();


/*
 * PLINIT: initalize the plotting package.
 */

plinit(pltdev, xmin, xmax, ymin, ymax)
char *pltdev;		/* output device name (ignored in this yapp) */
double xmin, xmax;	/* user plotting area */
double ymin, ymax;
{
    float width, height;

    initialize_core(BASIC, SYNCHRONOUS, TWOD);
    surface = get_surface();
    initialize_view_surface(surface, FALSE);
    select_view_surface(surface);
    if (ymax - ymin < xmax - xmin) {        /* make a square area */
	dxymax = xmax - xmin;
	width = 1.0;
	height = (ymax - ymin) / dxymax;
    } else {
	dxymax = ymax - ymin;
	width = (xmax - xmin) / dxymax;
	height = 1.0;
    }
    set_ndc_space_2(width, height);
    set_viewport_2(0.0, width, 0.0, height);
    set_window(xmin, xmax, ymin, ymax);
    ax = (xmax-xmin)/width;                 /* conversion NDC to WC */
    ay = (ymax-ymin)/height;                /* x_NDC = ax.x_WC+bx   */
    bx = xmin;                              /* y_NDC = ay.y_WC+by   */
    by = ymin;
    set_font(ROMAN);
    set_charprecision(CHARACTER);
    initialize_device(KEYBOARD, 1);
    set_echo(KEYBOARD, 1, 0);		/* no echo */
#if defined(MOUSE)
    initialize_device(BUTTON, 1);
    initialize_device(BUTTON, 2);
    initialize_device(BUTTON, 3);
    set_echo_surface(BUTTON,1,surface);
    set_echo_surface(BUTTON,2,surface);
    set_echo_surface(BUTTON,3,surface);
    initialize_device(LOCATOR, 1);
    set_echo_surface(LOCATOR,1,surface);
    set_echo(LOCATOR,  1, 1);           /* show the mouse as a finger */
#endif

#if defined(RETAINED)
    create_retained_segment(++segment_name);
#else
    create_temporary_segment();
#endif
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
    if (lwid > 0)
	set_linewidth(0.1 * (lwid - 1));
    if (lpat > 0)
	set_linestyle(lpat - 1);
}

/*
 * PLLINE, PLMOVE, PLPOINT: plot lines, moves, and points.
 */

plline(x, y)
double x, y;		/* user coordinates */
{
    line_abs_2(x, y);
}

plmove(x, y)
double x, y;		/* user coordinates */
{
    move_abs_2(x, y);
}

plpoint(x, y)
double x, y;		/* user coordinates */
{
    move_abs_2(x, y);
    line_rel_2(0.0, 0.0);
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
    double c, cos(), s, sin();
    float dx, dy;

    c = cos(ang / 57.296);
    s = sin(ang / 57.296);
    set_charpath_2(c, s);
    set_charup_2(- s, c);
    set_charsize(0.75 * hgt, hgt);
    inquire_text_extent_2(msg, &dx, &dy);
    move_abs_2(x - 0.58 * (1 + textjust) * dx,
               y - 0.55 * (1 + textjust) * dy);
				/* note fudge factors in offsets */
    text(msg);
}

/*
 * PLFLUSH: output any pending graphics.
 */

plflush() { }

/*
 * PLFRAME: advance to next frame, with optional pause
 */

plframe()
{
#if defined(RETAINED)
    if (keep_old_segments)
        close_retained_segment();
    else
        delete_retained_segment(segment_name);
#else
    close_temporary_segment();
#endif
    new_frame();                        /* clear screen and that stuff */
#if defined(RETAINED)
    create_retained_segment(++segment_name);
#else
    create_temporary_segment();
#endif
}

/*
 * PLSTOP: finish up a display.
 */

#define ONEMIN (60 * 1000 * 1000)

plstop()
{
    char msg[80];
    int len;

#if defined(RETAINED)
    close_retained_segment();
#else
    close_temporary_segment();
#endif

    bell();         /* always beep */

#if defined(MOUSE)
    printf ("Press any button to finish\n");
    len = 0;
    do {
        await_any_button(1,&len);
    } while (len==0);
    terminate_device(BUTTON,1);
    terminate_device(BUTTON,2);
    terminate_device(BUTTON,3);
    terminate_device(LOCATOR,1);
#else
    printf ("Press RETURN to finish\n");
    await_keyboard(ONEMIN, 1, msg, &len);
#endif
    deselect_view_surface(surface);
    terminate_core();
}

pl_screendump (fname)
char *fname;
{
    char cmd[64];

    sprintf (cmd,"screendump %s",fname);
    system(cmd);    /* execute */
}


pl_matrix(frame,nx,ny,xmin,ymin,cell,fmin,fmax,findex)
double *frame, xmin, ymin, cell, fmin, fmax, findex;
int nx, ny;
{
    double x,y,f,grayscale,ds,pow();
    int ix,iy;
    
#if TRUE
    printf("PL_MATRIX is unsupported in suncore version of YAPP\n");
    return(0);
#else
    /*
     * this is code from PS- must be converted to SunCore stuff
     * see:  set_fill_index
     *       define_color_indices
     */
    ds = cell*PTPCM;			/* convert cm -> pixels */
    printf("GRAY_MATRIX cm: %d * %d matrix  ll corner: %f %f cell %f\n",
	   nx,ny,xmin,ymin,cell);
    printf("GRAY_MATRIX DC: %d * %d matrix  ll corner: %f %f cell %f\n",
	   nx,ny,convx(xmin),convy(ymin),ds);	
    grayscale = 1.0/(fmax-fmin);
    				/* normally positive for a negative image */
    for (ix = 0, x = xmin; ix<nx; ix++, x += cell) {
	for (iy = 0, y = ymin; iy<ny; iy++, y+=cell) {
	    f = *(frame + ix*ny + iy);
					/* apply linear grayscale + cutoff */
	    if (grayscale>0.0)
		f *= grayscale;
	    else
		f = f*grayscale + 1.0;
	    if (f>1.0)			/* upper cutoff */
		f = 1.0;
	    if (f<0.0)			/* lower cutoff */
		f = 0.0;
	    f = pow(f,findex);		/* transfer function */
	    vec_paint (convx(x),convy(y),ds,f);		/* paint this square */
	}  
    }
#endif
}

local struct vwsurf *get_surface()
{
    if (getenv("WINDOW_ME"))
        return(&winsurf);           /* window mode */
    else
        return(&rawsurf);           /* full screen mode */
}

local bell()
{
    int ring=7;
#if 0
    set_echo(KEYBOARD,1,1);  
#endif
    putchar(ring);
    putchar('\n');		/* send a line feed to flush the buffer */
#if 0
    set_echo(KEYBOARD,1,0);   
#endif
}

/* Unsupported YAPP functions for mouse locator polygonic stuff in MOUSE mode */

pl_getpoly(x,y,n)
float x[];	/* x positions of polygon */
float y[];	/* y positions of polygon */
int  n;		/* maximum number in x,y to  be held */
{
    int nn, delay, loc, k;
    float xold,yold, xnew,ynew;

    bell();

#if defined(MOUSE)
    printf("Define a polygon:\n");
    printf("   LEFT   = define first/next vertex of polygon\n");
    printf("   MIDDLE = delete previous vertex point from polygon\n");
    printf("	    (this option cannot remove already plotted vertex lines\n");
    printf("   RIGHT  = close polygon\n");
    set_locator_2 (1,0.5,0.5);        /* set to dummy position */
    nn = 0;                           /* count points in polygon */
    while(1) {
        do {
            await_any_button_get_locator_2 (1, 1, &k, &xnew, &ynew);
        } while (k==0);
        if (k==1) {                         /* left button = DEFINE POLYGON */
            xnew=xnew*ax+bx;                   /* back to WC */
            ynew=ynew*ay+by;
            xold = xnew;
            yold = ynew;
            dprintf (2,"Button A at x=%f y=%f\n",xold,yold);
            if (nn==0) {
                move_abs_2(xold,yold);
                set_marker_symbol(43);      /* a 'plus' symbol */
                marker_abs_2(xold,yold);
            } else
                line_abs_2(xold,yold);
            if (nn<n) {
               x[nn] = xold;
               y[nn] = yold;
               nn++;
            } else
                printf ("Warning, no more space to store this data\n");
        } else if (k==2) {                  /* middle  = DELETE PREVIOUS POINT */
            if (nn>0)
                nn--;
            if (nn!=0)
               move_abs_2(x[nn-1],y[nn-1]);    /* reset */
        } else {                          /* right = QUIT */
            break;
        }
    }   /* while */
    dprintf (2,"Button C pressed to finish up polygon (%d)\n",nn);
    if (nn>0) {
        line_abs_2(x[0],y[0]);
    }
    return(nn<3 ? 0 : nn);        
#else
    error("No mouse functions compiled into YAPP\n");
#endif
}

#if 0
/* 	This function can be called to test the algorithm */

pl_testpoly (x,y,n)
float x[], y[];
int   n;
{
    int k,i, b;
    float xnew,ynew;

    printf("Testing points in/out of polygon\n");
    printf("   LEFT   = define point to test\n");
    printf("   MIDDLE,RIGHT = quit\n");
    set_locator_2 (1,0.5,0.5);        /* set to dummy position */
    x[n] = x[0];
    y[n] = y[0];
    while(1) {
        do {
            await_any_button_get_locator_2 (1, 1, &k, &xnew, &ynew);
        } while (k==0);
        if (k==1) {                         /* left button = DEFINE POLYGON */
            xnew=xnew*ax+bx;                   /* back to WC */
            ynew=ynew*ay+by;
            printf ("Button A at x=%f y=%f ",xnew,ynew);
            printf ("found inside = %d\n",inpolygon(n,x,y,xnew,ynew));
        } else
            break;
    }   /* while */
}

/*	Early version, for updates see snapplot.c or some library */
/*  The current algorithm has some flaws, may not be the fastest,
 *  (I have seen a faster one somewhere, but can't recall where),
 *  works. It's from CACM aug 1962, algorithm 112 by M. Shimrat (p434)
 *  but see also remark in CACM dec 1962, p606)
 */

/*		Code for this is now in editplot/snapplot */
inpolygon (n, x, y, x0, y0)
int n;
float x[], y[], x0, y0;
{
    int i,b;              /* b=0:outside b=1:inside */
        
    x[n] = x[0];        /* make sure polygon is closed */
    y[n] = y[0];

    b = 0;              /* set initially to false */
    for (i=0; i<n; i++) {
        if ( ((y0<y[i])==(y0>y[i+1])) &&
             ((x0-x[i]-(y0-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i])) < 0 ))
        b = !b;
    }
    return(b);
}
#endif

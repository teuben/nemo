/*
 * YAPP: Yet Another Plotting Package.
 * Joshua Barnes  Dec 1986  I. A. S.  Princeton, NJ.
 *
 *    29-jun-87  added pl_matrix() for yapp_ps consistency   PJT
 *    26-oct-90	 added pl_screendump for consistency	    PJT
 *    22-jun-91  include stdinc  BEFORE usercode
 *
 * This version works with the SunCore graphics library.
 */

#include <stdinc.h>
#include <usercore.h>


#ifndef FullScreen
  int pixwindd();		/* B&W window driver */
  local struct vwsurf vwsurf = DEFAULT_VWSURF(pixwindd);
#else
  int bw2dd();                  /* B&W full-screen driver */
  local struct vwsurf vwsurf = DEFAULT_VWSURF(bw2dd);
#endif

local double dxymax;	/* size of user window */
local int retval;       /* if retval != 0, no graphics calls processed */

/*
 * PLINIT: initalize the plotting package.
 */

plinit(pltdev, xmin, xmax, ymin, ymax)
char *pltdev;		/* output device name (ignored) */
double xmin, xmax;	/* user plotting area */
double ymin, ymax;
{
    float width, height;

    retval=initialize_core(BASIC, SYNCHRONOUS, TWOD);
    if (retval) {
        warning("Error_1 initializing graphics screen:\n\tare you on the console?");
        return;
    }
    retval= initialize_view_surface(&vwsurf, FALSE);
    if (retval) {
        warning("Error_2 initializing graphics screen:\n\tare you on the console?");
        return;
    }
    retval=select_view_surface(&vwsurf);
    if (retval) {
        warning("Error_3 initializing graphics screen:\n\tare you on the console?");
        return;
    }
    if (ymax - ymin < xmax - xmin) {
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
    set_font(ROMAN);
    set_charprecision(CHARACTER);
    initialize_device(KEYBOARD, 1);
    set_echo(KEYBOARD, 1, 0);		/* no echo */
    create_temporary_segment();
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
    if (retval) return;
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
    if (retval) return;
    line_abs_2(x, y);
}

plmove(x, y)
double x, y;		/* user coordinates */
{
    if (retval) return;
    move_abs_2(x, y);
}

plpoint(x, y)
double x, y;		/* user coordinates */
{
    if (retval) return;
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

    if (retval) return;
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
    if (retval) return;
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
    if (retval) return;
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

    if (retval) return;
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
 * PLFRAME: advance to next frame.
 */

plframe()
{
    if (retval) return;
    close_temporary_segment();
    new_frame();
    create_temporary_segment();
}

/*
 * PLSTOP: finish up a display.
 */

#define ONEMIN (60 * 1000 * 1000)

plstop()
{
    char msg[80];
    int len;

    if (retval) return;
    close_temporary_segment();
    await_keyboard(ONEMIN, 1, msg, &len);
    deselect_view_surface(&vwsurf);
    terminate_core();
}

pl_matrix(frame,nx,ny,xmin,ymin,cell,fmin,fmax,findex,blank)
double *frame, xmin, ymin, cell, fmin, fmax, findex, blank;
int nx, ny;
{
    double x,y,f,grayscale,ds,pow();
    int ix,iy;

    if (retval) return;    
    warning("pl_matrix: unsupported in YAPP_CORE");
    return(0);
}

pl_screendump(fname)
string fname;
{
    if (retval) return;
    warning("pl_screendump: not implemented in yapp_core");
}

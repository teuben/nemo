/*
 * YAPP-CG: Yet Another Plotting Package -- Color Graphics.
 * Joshua Barnes	Dec 1986  I. A. S.  Princeton, NJ.
 * Peter Wisnovsky	May 1987  I. A. S.  Princeton, NJ.
 *
 * This version works with the SunCore graphics library,
 * and drives a color display through a simple interface.
 *
 *  2-jul-87  added dummy call to pl_matrix() for consistency
 */

#include <usercore.h>
#include <stdinc.h>

#include "get_view_surface.c"

local struct vwsurf vwsurf;

local float red[256];
local float blue[256];
local float green[256];

/*
 * Definition of the colormap segment CMS_RGB, a collection of rgb
 * (red-green-blue) binary combinations.
 */

#define	BLACK	0
#define	RED	1
#define	YELLOW	2
#define	GREEN	3
#define	CYAN	4
#define	BLUE	5
#define	MAGENTA	6
#define	WHITE	7

#define	cms_rgbsetup(r,g,b)						\
{									\
    (r)[BLACK] = 0.0;	(g)[BLACK] = 0.0;	(b)[BLACK] = 0.0;	\
    (r)[RED] = 1.0;	(g)[RED] = 0.0;		(b)[RED] = 0.0;		\
    (r)[YELLOW] = 1.0;	(g)[YELLOW] = 1.0;	(b)[YELLOW] = 0.0;	\
    (r)[GREEN] = 0.0;	(g)[GREEN] = 1.0;	(b)[GREEN] = 0.0;	\
    (r)[CYAN] = 0.0;	(g)[CYAN] = 1.0;	(b)[CYAN] = 1.0;	\
    (r)[BLUE] = 0.75;	(g)[BLUE] = 0.75;	(b)[BLUE] = 1.0;	\
    (r)[MAGENTA] = 1.0;	(g)[MAGENTA] = 0.75;	(b)[MAGENTA] = 1.0;	\
    (r)[WHITE] = 1.0;	(g)[WHITE] = 1.0;	(b)[WHITE] = 1.0;	\
}
  
local double dxymax;	/* size of user window */

/*
 * PLINIT: initalize the plotting package.
 */

plinit(pltdev, xmin, xmax, ymin, ymax)
char *pltdev;		/* output device name (ignored) */
double xmin, xmax;	/* user plotting area */
double ymin, ymax;
{
    float width, height;
    string nullargv[1];

    nullargv[0] = NULL;
    initialize_core(BASIC, SYNCHRONOUS, TWOD);
    my_get_view_surface(&vwsurf, nullargv, 8);
    initialize_view_surface(&vwsurf, FALSE);
    select_view_surface(&vwsurf);
    cms_rgbsetup(red, green, blue);
    define_color_indices(&vwsurf, 0, 8, red, green, blue);
    plcolor(2.0);
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
 * PLCOLOR: simple function to specify plotting color: if less then
 * 0.0 then black, if within 0.0 to 1.0 (inclusive) then a cycle of
 * saturated colors, if greater than 1.0 then white.
 */

plcolor(color)
double color;
{
    int icol;

    if (color < 0.0)
	icol = BLACK;
    else if (color > 1.0)
	icol = WHITE;
    else
	icol = 1 + (int) (6 * color);
    set_line_index(icol);
    set_text_index(icol);
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
    float dx, dy;

    c = cos(ang / 57.296);
    s = sin(ang / 57.296);
    set_charpath_2(c, s);
    set_charup_2(- s, c);
    set_charsize(0.75 * hgt, hgt);
    inquire_text_extent_2(msg, &dx, &dy);
    move_abs_2(x - 0.55 * (1 + textjust) * dx,
               y - 0.55 * (1 + textjust) * dy);
				/* note fudge factor of 1.1 in offsets */
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

    close_temporary_segment();
    await_keyboard(ONEMIN, 1, msg, &len);
    deselect_view_surface(&vwsurf);
    terminate_core();
}


pl_matrix (frame,nx,ny,xmin,ymin,cell,fmin,fmax,findex)
double *frame, xmin, ymin, cell, fmin, fmax, findex;
int nx, ny;
{
	double x,y,f,grayscale,ds,pow();
	int    ix,iy;
	if (1) {
		printf ("PL_MATRIX is still unsupported in suncore version of YAPP_CG\n");
		return(0);
	}
}

#ifdef TESTBED

main(argc, argv)
int argc;
string argv[];
{
    int i, j;

    plinit("***", 0.0, 20.0, 0.0, 20.0);
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
    for (i = 0; i < 32; i++ ) {
	plcolor(i / 31.0);
	plmove(12.0, 2.0 + i / 31.0);
	plline(13.0, 2.0 + i / 31.0);
    }
    plcolor(2.0);
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

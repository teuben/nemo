/*
 * YAPP_GL.C: Yet Another Plotting Package -- SiliconGraphics GL - OpenGL - MESA 
 * Version 1.0	Joshua Barnes	14 Oct  89	CITA, Toronto.
 *         1.1  pjt              9 jul  91  added dummy pl_screendump() for
 *					    testyapp.c - NEMO V2
 *	   1.2  pjt             30 dec  92  testing - using VOGL library on SUN
 *			                    minor NEMO V2.x mods - -DVOGL
 *         1.3  pjt              1 apr  97  added plswap,plxscale,plyscale
 *					    made it SINGLEPREC friendly
 *         1.3a pjt              9 may  97  fixed bug for non-JOSH code
 *                                          variables not initialized
 *                                          fix dec alpha plstop()
 *	   1.3b pjt             18 nov  99  Added MESA support (-DMESA)
 *                                          --but not working yet--
 *					    see separate _mesa version?
 */

#include <stdinc.h>
#include <getparam.h>
#ifdef SGI
# include <gl.h>
# include <device.h>
  static string yapp_gl_id = "yapp_gl=SGI";
#elif defined(YGL)
# include <ygl/Ygl.h>
  static string yapp_gl_id = "yapp_gl=YGL";
#elif defined(MESA)
# include <GL/gl.h>
  static string yapp_gl_id = "yapp_gl=MESA";
# if 1
#  define BLACK			0
#  define WHITE 		1
#  define KEYBD 		257
#  define INPUTCHANGE		262
#  define ESCKEY		'\033'
# endif
#else
# include <vogl/gl.h>                           /* mostly:  NEMOINC/vogl */
# include <vogl/device.h>
  static string yapp_gl_id = "yapp_gl=VOGL";
#endif

#define INACTIVE   0				/* not yet initialized	    */
#define ACTIVE     1				/* ready to draw graphics   */

local int yappstat = INACTIVE;			/* status of yapp	    */
local long yappwind;				/* id. of yapp window       */
local real pixscale;				/* pixel size in user c.s.  */
local bool twobuf;				/* true if double-buffered  */
local bool gfxpnd;				/* true if graphics pending */
local int ncolors;				/* number of colors defined */
local int just;					/* text justification param */
local string inmsg = "yapp: not initialized\n";	/* a common error message   */
local short l_styles[] = {
    0xFFFF, 0xFFF0, 0xCFF0, 0xF0F0, 0xFF66, 0};



/*
 * PLINIT: initialize YAPP_GT.
 */

plinit(string opt, real x0, real x1, real y0, real y1)
{
    real dx, dy;
    int idx, idy;
    long xsz, ysz;

    if (yappstat != INACTIVE)				/* check: 2nd call? */
	error("yapp: already initialized\n");
    yappstat = ACTIVE;					/* make yapp active */
    dx = x1 - x0;					/* find width and   */
    dy = y1 - y0;					/*  height of area  */
    idx = 1024 * dx / MIN(dx, dy);			/* compute integer  */
    idy = 1024 * dy / MIN(dx, dy);			/*  width, height   */
    keepaspect(idx, idy);				/* set window shape */
    foreground();					/* for user input   */
#ifdef JOSH
    yappwind = winopen(getargv0());			/* open yapp window */
    getsize(&xsz, &ysz);				/* get window size  */
#else
    prefsize(512L,512L);
    winopen(getargv0());
    /* hfont("times.rb"); */
    /* htextsize(0.1,0.2); */
    xsz = 512;
    ysz = 512;
#endif
    pixscale = dx / xsz;				/* find pixel size  */
    ortho2(x0, x1, y0, y1);				/* set user coords  */
    if (scanopt(opt, "twobuf")) {			/* double-buf mode? */
	twobuf = TRUE;					/*   set mode flag  */
	doublebuffer();					/*   set IRIS mode  */
	gconfig();					/*   and reconfig   */
    } else						/* single-buf mode? */
	twobuf = FALSE;					/* clear mode flag  */
    gfxpnd = FALSE;					/* no gfx pending   */
    color(BLACK);					/* set bkgnd color  */
    clear();						/* and erase screen */
    color(WHITE);					/* default to white */
    ncolors = 1 + WHITE;				/* for now...	    */
    for (idx=0; l_styles[idx]; idx++)
        deflinestyle((short)(idx+1), l_styles[idx]);    /* define ltype's   */
}

/*
 * PLFLUSH, PLSWAP, PLXSCALE, PLYSCALE: currently do nothing.
 */

plflush()
{}

plswap()
{}

real plxscale(real x, real y)
{}

real plyscale(real x, real y)
{}



/*
 * PLFRAME: advance to next frame.
 */

plframe()
{
    if (yappstat == INACTIVE)				/* skipped plinit?  */
	error(inmsg);
    if (twobuf)						/* double-buf mode? */
	swapbuffers();					/*   show back buf  */
    gfxpnd = FALSE;					/* no gfx pending   */
    color(BLACK);					/* set bkgnd color  */
    clear();						/* and erase screen */
    color(ncolors-1);					/* default to white */
}

/*
 * PLSTOP: finish up and quit.
 */

plstop()
{
    short dev, val;

    if (yappstat == INACTIVE)
	error(inmsg);
    if (twobuf && gfxpnd)				/* any pending gfx? */
	swapbuffers();					/*   then show it   */
#ifdef JOSH
    qdevice(ESCKEY);					/* listen for <esc> */
    while ((dev = qread(&val)) != ESCKEY)		/* loop till read   */
	;						/*   doing nothing  */
    winclose(yappwind);					/* close the window */
#else
    qdevice(KEYBD);
    unqdevice(INPUTCHANGE);
#if 0
    qread(&val);   
#else
    /* digital unix  may 1997 */
    /* bug: it doesn't read the mouse keys ???? */
    for (;;) {
      qread(&val);   
      if (val == ESCKEY) break;
    }
#endif
    gexit();
#endif
    yappstat = INACTIVE;				/* and shut down    */
}

/*
 * PLCOLOR: specify plotting color.
 */

plcolor(int col)
{
    if (yappstat == INACTIVE)
	error(inmsg);
    color(MAX(0, MIN(col, ncolors-1)));
}

/*
 * PLNCOLORS: return number of colors defined; valid colors range
 * range from zero to ncolors-1.
 */

int plncolors()
{
    if (yappstat == INACTIVE)
	error(inmsg);
    return (ncolors);
}

/*
 * PLPALETTE: specify new color table.
 */

plpalette(real *r, real *g, real *b, int nc)
{
    int i, ir, ig, ib;

    if (yappstat == INACTIVE)
	error(inmsg);
    for (i = 0; i < nc; i++) {
	ir = 255 * MAX(0.0, MIN(r[i], 1.0));
	ig = 255 * MAX(0.0, MIN(g[i], 1.0));
	ib = 255 * MAX(0.0, MIN(b[i], 1.0));
	mapcolor(i, ir, ig, ib);
    }
    ncolors = nc;
}

/*
 * PLLTYPE: select line width and pattern.
 */

plltype(int lwid, int lpat)
{
    if (yappstat == INACTIVE)
	error(inmsg);
    if (lwid)  linewidth((short)lwid);
    if (lpat)  setlinestyle((short)lpat);
}

/*
 * PLMOVE, PLLINE, PLPOINT: plot lines, moves, and points.
 */

plmove(real x, real y)
{
    if (yappstat == INACTIVE)
	error(inmsg);
    move2(x, y);
    gfxpnd = TRUE;					/* flag pending gfx */
}

plline(real x, real y)
{
    if (yappstat == INACTIVE)
	error(inmsg);
    draw2(x, y);
    gfxpnd = TRUE;
}

plpoint(real x, real y)
{
    if (yappstat == INACTIVE)
	error(inmsg);
    pnt2(x, y);
    gfxpnd = TRUE;
}

/*
 *  PLCIRCLE, PLCROSS, PLBOX: draw useful marking symbols.
 */
 
plcircle(real x, real y, real r)
{
    if (yappstat == INACTIVE)
	error(inmsg);
    if (r > 0.0)
	circf(x, y, r);
    else
	circ(x, y, -r);
    gfxpnd = TRUE;
}

plcross(real x, real y, real s)
{
    if (s > 0.0) {
	plmove(x-s, y);
	plline(x+s, y);
	plmove(x, y-s);
	plline(x, y+s);
    } else {
	s = s / 1.4142;
	plmove(x-s, y-s);
	plline(x+s, y+s);
	plmove(x-s, y+s);
	plline(x+s, y-s);
    }
}

plbox(real x, real y, real s)
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
 * PLTEXT: plot text string with specified position, size and angle.
 *	note: hgt and ang are currently ignored !!!
 */

pltext(string msg, real x, real y, real hgt, real ang)
{
    real sw, sh;

    if (yappstat == INACTIVE)
	error(inmsg);
    sw = pixscale * strwidth(msg);			/* find msg width   */
#ifdef JOSH
    sh = pixscale * (getheight() - getdescender());	/* find msg height  */
#else
    sh = pixscale * (getheight());	/* find msg height  */
#endif
    cmov2(x - sw * (just + 1) / 2.0, y - 0.3 * sh);	/* set char origin  */
    rot((float)ang,'z');
    charstr(msg);					/* plot text string */
    rot(0.0,'z');
    gfxpnd = TRUE;
}

/*
 * PLJUST: specify justification parameter.
 */

pljust(int j)
{
    if (yappstat == INACTIVE)
	error(inmsg);
    just = (j <= -1 ? -1 : (j >= 1 ? 1 : j));
}

pl_screendump(string fname)
{
    warning("pl_screendump: not implemented for YAPP_GT");
}

#ifdef TESTBED_JOSH

string defv[] = {
    "opt=",
    "x0=0.0",
    "x1=20.0",
    "y0=-2.5",
    "y1=22.5",
    "VERSION=1.0",
    NULL,
};

main(argc, argv)
int argc;
string argv[];
{
    int i, j, nc;

    initparam(argv, defv);
    plinit(getparam("opt"),
	   getdparam("x0"), getdparam("x1"),
	   getdparam("y0"), getdparam("y1"));
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
    plltype(1, 0);
    nc = plncolors();
    for (i = 1; i < nc; i++) {
	plcolor(i);
	plmove(6.0 + 2 * (i - 1.0) / nc, 16.0);
	plline(6.5 + 2 * (i - 1.0) / nc, 17.0);
    }
    for (i = 1; i <= 5; i++) {
	plltype(2*i, 1);
        plmove(1.0, 13.5 - i);
        plline(3.0, 13.5 - i);
        plpoint(3.5, 13.5 - i);
	plltype(1, i);
	for (j = 1; j <= 4; j++) {
	    plmove(1.5, 13.5 - i - 0.2*j);
	    plline(1.5 + j, 13.5 - i - 0.2*j);
	}
    }

    plltype(1, 1);
    for (i = 0; i < 8; i++) {
	plcross(  7.5, 2.0 + 0.5*i, -0.4 / (1 + i));
	plbox(    8.5, 2.0 + 0.5*i, -0.4 / (1 + i));
	plcircle( 9.5, 2.0 + 0.5*i, -0.4 / (1 + i));
	plcircle(10.5, 2.0 + 0.5*i,  0.4 / (1 + i));
	plbox(   11.5, 2.0 + 0.5*i,  0.4 / (1 + i));
	plcross( 12.5, 2.0 + 0.5*i,  0.4 / (1 + i));
    }
    pltext("Foo Bar!", 8.0, 16.0, 0.5, 0.0);
    pltext("Fum Bar!", 8.0, 15.0, 0.25, 0.0);
    for (i = 0; i <= 4; i++)
	pltext(" testing angles", 16.0, 10.0, 0.32, 45.0*i);
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



/*
 * GRX driver for VOGL c1993 by Gary Murphy (garym@virtual.rose.utoronto.ca)
 * 
 * To compile:
 * 
 * 1) add GRX to device.c and mash-up your makefiles for MsDOS 2) compile with
 * DOBJ=-DPOSTSCRIPT -DHPGL -DGRX and MFLAGS=-O2
 * 
 * To run:
 * 
 * set VDEVICE=grx
 * 
 * grateful thanks to Lance Norskog (thinman@netcom.com) and Bernie Kirby
 * (bernie@ecr.mu.oz.au) --- should either of you be in my neighbourhood, my
 * offer of an Ice Beer is still open! (Things in #ifdef BLART disabled by
 * bernie...)
 */
#undef VOGLE

#include <stdio.h>
#undef DBG
#ifdef DBG
FILE	*dfp = NULL;
#endif
#include <assert.h>

#include <stdlib.h>
#include <memory.h>
#include <grx.h>
#include <mousex.h>

#define MSG( m ) fprintf(stderr, "\n%s: %d: %s", __FILE__, __LINE__, (m))
#define ERROR1( m, p ) fprintf(stderr, "\n%s: %d: " m, __FILE__, __LINE__, (p))

#ifdef VOGLE
#include	"vogle.h"
#else
#include	"vogl.h"
#endif

#ifndef TRUE
#define TRUE	1
#endif

#ifndef FALSE
#define FALSE	0
#endif

#define MAXCOLOR 256

static struct {

	GR_graphics_modes old_mode;

	int             width, height, planes;
	unsigned        scrsize;/* size of buffer in long words */
	GrContext      *cbuf;	/* current context */
	GrContext      *fbuf;
	GrContext      *bbuf;

	int             palette[8];

	GrLineOption    lopt;	/* pen drawing options */
	int             fg;	/* foreground/background colours */
	int             bg;

	int             has_mouse;

	GrFont         *font;	/* Current font */
	GrFont         *lfont;	/* Preloaded small font */
	GrFont         *sfont;	/* Preloaded large font */
	char           *fname;	/* Fontname */
	GrTextOption   to;	/* Other text stuff */

	int             cx;
	int             cy;

} grx;

/*
 * access functions: *
 * 
/* I'm going to need this to fudge in stereo graphics ...
 */

GrContext      *
setBackBuffer(GrContext * newBB)
{
	GrContext      *oldBB = grx.bbuf;
	assert(newBB != NULL);

	grx.bbuf = newBB;
	return oldBB;
}

static int
grx_init()
{
#ifdef DBG
	dfp = fopen("grx.dbg", "w");
#endif
	grx.old_mode = GrCurrentMode();
	GrSetMode(GR_default_graphics);

#ifndef VOGLE
	vdevice.devname = "Grx";
#endif

	/* set the VOGL device */
	vdevice.sizeX = GrSizeY();	/* square max, was GrScreenX(); */
	vdevice.sizeY = GrSizeY();

	grx.width = vdevice.sizeSx = GrScreenX();
	grx.height = vdevice.sizeSy = GrScreenY();
	grx.planes = vdevice.depth = GrNumPlanes();

	grx.scrsize = (GrPlaneSize(grx.width, grx.height) * grx.planes) / sizeof(long);

	/* setup default palette */
	GrSetRGBcolorMode();
	grx.lopt.lno_color = grx.fg = GrWhite();
	grx.bg = GrBlack();

	grx.palette[BLACK] = GrAllocColor(0, 0, 0);
	grx.palette[RED] = GrAllocColor(255, 0, 0);
	grx.palette[GREEN] = GrAllocColor(0, 255, 0);
	grx.palette[YELLOW] = GrAllocColor(255, 255, 0);
	grx.palette[BLUE] = GrAllocColor(0, 0, 255);
	grx.palette[MAGENTA] = GrAllocColor(255, 0, 255);
	grx.palette[CYAN] = GrAllocColor(0, 255, 255);
	grx.palette[WHITE] = GrAllocColor(255, 255, 255);

	/*
	 * setup back/front buffers: frontbuffer is the current screen, back
	 * is a ram context
	 */
	grx.cbuf = grx.fbuf = GrSaveContext(NULL);
	grx.bbuf = NULL;

	/* initialize mouse */
	if ((grx.has_mouse = MouseDetect()) == TRUE) {
		/* dare I do interrupts? ... */
		MouseEventMode(1);
		MouseInit();

		/* no keyboard (use getch) */
		MouseEventEnable(0, 1);

		/* cheezy mouse speed algorithm (blame Lance for the pun) */
		if (grx.width * grx.height < 100000)
			MouseSetSpeed(6);
		else if (grx.width * grx.height < 200000)
			MouseSetSpeed(4);
		else if (grx.width * grx.height < 500000)
			MouseSetSpeed(3);
		else
			MouseSetSpeed(2);

		MouseWarp(1, 1);
		MouseDisplayCursor();
	};

	/* initial drawing style to thin solid lines */
	grx.lopt.lno_width = 1;
	grx.lopt.lno_pattlen = 0;
	grx.lopt.lno_dashpat = NULL;
	/* load initial fonts */
	if (getenv("GRXFONT") == NULL)
		GrSetFontPath("fonts");

	grx.font = grx.sfont = GrLoadFont(vdevice.dev.small);
#ifdef DBG
	if (!grx.font) {
		fprintf(dfp, "GrLoadFont failed");
		fflush(dfp);
	}
#endif

	grx.lfont = GrLoadFont(vdevice.dev.large);
#ifdef DBG
	if (!grx.lfont) {
		fprintf(dfp, "GrLoadFont failed");
		fflush(dfp);
	}
#endif

	grx.to.txo_font = grx.font;
	grx.to.txo_xmag = grx.to.txo_ymag = 1;
	grx.to.txo_direct = GR_TEXT_RIGHT;
	grx.to.txo_xalign = grx.to.txo_xalign = GR_ALIGN_DEFAULT;
	grx.to.txo_fgcolor = 0;

	grx.fname = vdevice.dev.small;
	vdevice.hwidth = 8.0;
	vdevice.hheight = 8.0;

	return (1);
}

/*
 * grx_frontbuffer, grx_backbuffer, grx_swapbuffers
 * 
 */
static
int 
grx_frontbuffer()
{
	grx.cbuf = grx.fbuf;
	GrSetContext(grx.fbuf);
	return (0);
}

static
int 
grx_backbuffer()
{
	/* if they want a backbuffer, we'd better make one ... */

	if (grx.bbuf == NULL)
		grx.bbuf = GrCreateContext(GrSizeX(), GrSizeY(), NULL, NULL);

	assert(grx.bbuf != NULL);

	grx.cbuf = grx.bbuf;
	GrSetContext(grx.bbuf);
	return (0);
}

static
int 
grx_swapbuffers()
{
	if (grx.cbuf == grx.fbuf)
		grx_backbuffer();
	else {
		/*
		 * there are rumours of a portable VGA backbuffer using VESA
		 * but I've yet to track it down.
		 * 
		 * the following copies by long words from back to front buffer
		 * modify this for regions by triming the first x-limit and
		 * y-limit and upping the pointers to the start of your
		 * subcontext
		 */

		MouseEraseCursor();

		/* WARNING WILL ROBINSON - WARNING WILL ROBINSON */
		/*
		 * We're using the NC version so I can copy a 2W by H/2
		 * backbuffer in my stereo graphics 'interlaced' mode
		 */

		GrBitBltNC(grx.fbuf, 0, 0,
			   grx.bbuf, 0, 0,
			   grx.bbuf->gc_xmax, grx.bbuf->gc_ymax, GrWRITE);

		MouseDisplayCursor();


	}

	return (0);
}

#ifdef VOGLE
/*
 * grx_vclear
 * 
 * Clear the screen to current colour
 */
grx_vclear()
{
	grx.to.txo_bgcolor = grx.fg;
	GrClearContext(grx.fg);

}

#else

/*
 * grx_vclear
 * 
 * Clear the viewport to current colour
 */
static int
grx_vclear()
{
	unsigned int    vw = vdevice.maxVx - vdevice.minVx;
	unsigned int    vh = vdevice.maxVy - vdevice.minVy;

	grx.to.txo_bgcolor = grx.fg;
	if ((vdevice.sizeSx == vw) && (vdevice.sizeSy == vh)) {
		GrClearContext(grx.fg);	/* full screen */
	} else
		GrFilledBox(
			    vdevice.minVx,
			    vdevice.sizeSy - vdevice.maxVy,
			    grx.width,
			    grx.height,
			    grx.fg);

	return (1);
}

#endif

/*
 * grx_exit
 * 
 * Sets the display back to text mode.
 */
static
grx_exit()
{
	MouseUnInit();		/* disable mouse/keyboard interrupts */

	GrSetMode(grx.old_mode);
	GrDestroyContext(grx.bbuf);

	return (1);
}

static int 
grx_sync()
{
};

static int
noop()
{
	return (-1);
}

/*
 * grx_font : load either of the fonts
 */
static int
grx_font(char *name)
{
#ifdef DBG
	fprintf(dfp, "fontname=%s\n", name);
	fflush(dfp);
#endif
	/*
	 * Hacky way to quicky test for small or large font
	 * ... see of they are the same pointers.... this
	 * assumes that they have been called from the main
	 * library routine with *vdevice.Vfont(vdevice.smallfont);
	 */
	if (name == vdevice.dev.small) {
		grx.font = grx.sfont;
		grx.to.txo_font = grx.font;
		vdevice.hheight = grx.font->fnt_height;
		vdevice.hwidth = grx.font->fnt_width;
		grx.fname = name;
#ifdef DBG
		fprintf(dfp, "w, h: %f %f\n", vdevice.hheight, vdevice.hwidth);
		fflush(dfp);
#endif
		return(1);
	} else if (name == vdevice.dev.large) {
		grx.font = grx.lfont;
		grx.to.txo_font = grx.font;
		vdevice.hheight = grx.font->fnt_height;
		vdevice.hwidth = grx.font->fnt_width;
		grx.fname = name;
#ifdef DBG
		fprintf(dfp, "w, h: %f %f\n", vdevice.hheight, vdevice.hwidth);
		fflush(dfp);
#endif
		return(1);
	} else

	/* 
	 * It must be a completely different font (ala vogle possibility).
	 */
	if (strcmp(name, grx.fname)) {
		if (grx.font != grx.sfont && grx.font != grx.lfont)
			GrUnloadFont(grx.font);

		if (grx.fname) {
			free(grx.fname);
			grx.fname = (char *)malloc(strlen(name) + 1);
			strcpy(grx.fname, name);
		}

		if (!(grx.font = GrLoadFont(name))) {
			return(0);
		}

		grx.to.txo_font = grx.font;

		vdevice.hheight = grx.font->fnt_height;
		vdevice.hwidth = grx.font->fnt_width;
#ifdef DBG
		fprintf(dfp, "w, h: %f %f\n", vdevice.hheight, vdevice.hwidth);
		fflush(dfp);
#endif
	}

	return (1);
}

static
int 
grx_char(int c)
{
	GrDrawChar(c, vdevice.cpVx, vdevice.sizeSy - vdevice.cpVy - vdevice.hheight, &grx.to);
	vdevice.cpVx += vdevice.hwidth;

	return (1);
};

static int
grx_string(char *s)
{
	int	len = strlen(s);
	GrDrawString(s, len, vdevice.cpVx, vdevice.sizeSy - vdevice.cpVy - vdevice.hheight, &grx.to);
	vdevice.cpVx += vdevice.hwidth * len;
	return (1);
}


/*
 * Everything is supposed to have been through the higher up clippers in
 * vogl.. so no need to clip here..
 * 
 * Draw a solid 1 pixel wide line... libgrx checks for horizontal and vertical
 * lines for us.
 */
static int
grx_solid(int x, int y)
{

	GrLineNC(x, vdevice.sizeSy - y,
		 vdevice.cpVx, vdevice.sizeSy - vdevice.cpVy,
		 grx.fg
	);

	vdevice.cpVx = x;
	vdevice.cpVy = y;

	return (0);
}

/*
 * Draw a patterned and/or > 1 pixel wide line. (Waiting for libgrx to
 * actually implement this...)
 */
static int
grx_pattern(int x, int y)
{
	GrCustomLine(vdevice.cpVx, vdevice.sizeSy - vdevice.cpVy,
		     x, vdevice.sizeSy - y,
		     &grx.lopt
	);

	vdevice.cpVx = x;
	vdevice.cpVy = y;

	return (0);
};

static int
grx_colour(int i)
{

	if (i < MAXCOLOR)
		grx.fg = grx.palette[i];	/* for now */
	else
		grx.fg = GrBlack();

	grx.lopt.lno_color = grx.fg;

	grx.to.txo_fgcolor = grx.fg;
	return (0);
};

/*
 * grx_mapcolor
 * 
 * change index i in the color map to the appropriate r, g, b, value.
 */
static int
grx_mapcolor(int c, int r, int g, int b)
{
	int             j;

	if (c >= MAXCOLOR || vdevice.depth == 1)
		return (-1);

	grx.palette[c] = GrAllocColor(r, g, b);

}


static int
grx_fill(int sides, int *x, int *y)
{
	int             i, j;
	int             points[sides][2];

	for (i = 0; i < sides; i++) {
		points[i][0] = x[i];
		points[i][1] = grx.height - y[i];
	}

	GrFilledPolygon(sides, points, grx.fg);

	return (0);
};

static int
grx_checkkey()
{
	char            c;

	if (kbhit()) {
		if ((c = getkey()) == 3) {	/* control-c */
			grx_exit();
			/* don't call vexit(), avoid back-refs */
			exit(0);
		} else
			return c;
	} else
		return 0;
}

static int
grx_locator(int *x, int *y)
{
	MouseEvent      mEv;
	static          ox = 0, oy = 0, obuttons = 0;

	if (!grx.has_mouse) {
		*x = *y = 0;
		return (-1);
	}
	/*
	 * if (MousePendingEvent()) {
	 */
	MouseGetEvent(M_MOTION | M_BUTTON_CHANGE | M_POLL, &mEv);

	if (mEv.flags & M_MOTION) {
		ox = mEv.x;
		oy = vdevice.sizeSy - mEv.y;
	}
	/*
	 * HACK... the RIGHT button is the second button and we want it to
	 * return 2...
	 */

	if (mEv.flags & M_BUTTON_CHANGE) {
		obuttons = ((mEv.buttons & M_LEFT) ? 1 : 0) |
			((mEv.buttons & M_MIDDLE) ? 2 : 0) |
			((mEv.buttons & M_RIGHT) ? 2 : 0);
	}
	/*
	 * }
	 */

	*x = ox;
	*y = oy;

	return (obuttons);
}

static
int 
grx_lwidth(int w)
{

	grx.lopt.lno_width = w;
	if (w == 1 && grx.lopt.lno_pattlen == 0)
		vdevice.dev.Vdraw = grx_solid;
	else
		vdevice.dev.Vdraw = grx_pattern;
}

static
int 
grx_lstyle(int s)
{

	static unsigned char dashes[16];
	unsigned        ls = s;
	int             i, n, a, b;

	if (grx.lopt.lno_width == 1 && (ls == 0 || ls == 0xffff)) {
		vdevice.dev.Vdraw = grx_solid;
		grx.lopt.lno_pattlen = 0;
		return;
	}
	for (i = 0; i < 16; i++)
		dashes[i] = 0;

#define	ON	1
#define	OFF	0

	a = b = OFF;
	if (ls & (1 << 0))
		a = b = ON;

	n = 0;
	for (i = 0; i < 16; i++) {	/* Over 16 bits */
		if (ls & (1 << i))
			a = ON;
		else
			a = OFF;

		if (a != b) {
			b = a;
			n++;
		}
		dashes[n]++;
	}

	grx.lopt.lno_pattlen = 16;
	grx.lopt.lno_dashpat = dashes;
	vdevice.dev.Vdraw = grx_pattern;
}

static DevEntry grxdev = {
	"grx",
	"@:pc8x16.fnt",		/* Large font */
	"@:pc8x8.fnt",		/* Small font */
	grx_backbuffer,		/* backbuffer */
	grx_char,		/* hardware char */
	grx_checkkey,		/* keyhit */
	grx_vclear,		/* clear viewport to current colour */
	grx_colour,		/* set current colour */
	grx_solid,		/* draw line */
	grx_exit,		/* close graphics & exit */
	grx_fill,		/* fill polygon */
	grx_font,		/* set hardware font */
	grx_frontbuffer,	/* front buffer */
	getkey,			/* wait for and get key */
	grx_init,		/* begin graphics */
	grx_locator,		/* get mouse position */
	grx_mapcolor,		/* map colour (set indices) */
#ifndef VOGLE
	grx_lstyle,		/* set linestyle */
#endif
	grx_lwidth,		/* set line width */
	grx_string,		/* draw string of chars */
	grx_swapbuffers,	/* swap buffers */
	grx_sync		/* sync display */
};

/*
 * _grx_devcpy
 * 
 * copy the pc device into vdevice.dev. (as listed in drivers.c)
 */
_grx_devcpy()
{
	vdevice.dev = grxdev;
	return (0);
}

/*
 * VOGLE driver for X11 under the DEC window manager.
 *
 *	ok so we admit that we don't know heaps about X11, but come
 * on guys, this shouldn't be neccessary (wild comment from DEC about
 * why our standard X11 driver won't work is welcome.)
 */
#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "vogl.h"

#define LARGEX11R2	"courier12f.snf"
#define SMALLX11R2	"courier10f.snf"

#define LARGEX11R3	"-adobe-courier-medium-r-normal--24-240-75-75-m-150-iso8859-1"
#define SMALLX11R3	"-adobe-courier-medium-r-normal--10-100-75-75-m-60-iso8859-1"

#define MIN(x,y)	((x) < (y) ? (x) : (y))
#define	CMAPSIZE	256
#define	EV_MASK		KeyPressMask|ButtonReleaseMask|ExposureMask|ButtonPressMask

static	Window		winder;
static	Display		*display;
static	int		screen;
static	unsigned long	carray[CMAPSIZE];
static	Colormap	colormap;

static	Drawable	theDrawable;
static	GC		theGC;
static	Pixmap		bbuf;		/* Back buffer pixmap */
static	int		back_used;	/* Have we backbuffered ? */

static	XFontStruct	*font_id;
XEvent			event;

static	int		size;
static	unsigned long	colour;
static	unsigned int	h, w;

/*
 * DECX11_init
 *
 *	initialises X11 display.
 */
DECX11_init()
{
	int		i;
	int		x, y, prefx, prefy, prefxs, prefys;
	unsigned int	bw, depth, mask;
	Window		rootw, childw, topwinder, *kiddies;
	char		*av[2], name[50];

	XSetWindowAttributes    theWindowAttributes;
        XSizeHints              theSizeHints;
        unsigned long           theWindowMask;
	XWMHints                theWMHints;


	av[0] = "vogl.X11";
	av[1] = (char *)NULL;

	if ((display = XOpenDisplay((char *)NULL)) == (Display *)NULL) {
		fprintf(stderr,"vogl: DECX11_init: can't connect to X server\n");
		exit(1);
	}

	winder = DefaultRootWindow(display);
	screen = DefaultScreen(display);
	vdevice.depth = DefaultDepth(display, screen);
	colormap = DefaultColormap(display, screen);

	/*
	 * Set our standard colors...
	 */
	if (vdevice.depth == 1) {
		/*
		 * Black and white - anything that's not black is white.
		 */
		carray[0] = BlackPixel(display, screen);
		for (i = 1; i < CMAPSIZE; i++)
			carray[i] = WhitePixel(display, screen);
	} else {
		/*
		 * Color, try to get our colors close to what's in the
		 * default colormap.
		 */
		DECX11_mapcolor(0, 0, 0, 0);
		DECX11_mapcolor(1, 255, 0, 0);
		DECX11_mapcolor(2, 0, 255, 0);
		DECX11_mapcolor(3, 255, 255, 0);
		DECX11_mapcolor(4, 0, 0, 255);
		DECX11_mapcolor(5, 255, 0, 255);
		DECX11_mapcolor(6, 0, 255, 255);
		DECX11_mapcolor(7, 255, 255, 255);
	}

	getprefposandsize(&prefx, &prefy, &prefxs, &prefys);

	XQueryPointer(display, winder, &rootw, &childw, &x, &y, &x, &y, &mask);

	if (childw == None)
		childw = rootw;

	/*
	 * there is something very weird about dec's window manager
	 * as expressed on a dec station, to get the details for the window
	 * we are actually in we have to get the root window of the childw
	 * which gives us the root window of the real window stack, and then
	 * we use XQueryPointer to find the real child window.
	 */

	XQueryTree(display, childw, &rootw, &rootw, &kiddies, &i);

	topwinder = winder;
	winder = kiddies[0];
	XQueryPointer(display, winder, &rootw, &childw, &x, &y, &x, &y, &mask);

	XGetGeometry(display, childw, &rootw, &x, &y, &w, &h, &bw, &depth);

	if (prefx > -1) {
	        x = prefx;
	        y = prefy;
	}

	if (prefxs > -1) {
	        w = prefxs;
	        h = prefys;
	}

	x += bw;
	y += bw;

	w -= 2 * bw;
	h -= 2 * bw;

	theWindowAttributes.override_redirect = False;

        /*theWindowMask = CWBackPixel|CWBorderPixel|CWOverrideRedirect;*/

        theWindowMask = CWOverrideRedirect;

        winder = XCreateWindow(display,
                                topwinder,
                                x, y,
                                w, h,
                                bw,
                                (int)depth,
                                InputOutput,
                                CopyFromParent,
                                theWindowMask,
                                &theWindowAttributes
                        );

        theWMHints.initial_state = NormalState;
        theWMHints.flags = StateHint;
        XSetWMHints(display, winder, &theWMHints);
 
        theSizeHints.flags = PPosition|PSize;
        theSizeHints.x = x;
        theSizeHints.y = y;
        theSizeHints.width = w;
        theSizeHints.height = h;
 
        XSetNormalHints(display, winder, &theSizeHints);

	sprintf(name, "vogl %d", getpid());

	XSetStandardProperties(display,
		winder,
		name,
		name,
		None,
		av,
		1,
		&theSizeHints
	);

	XSelectInput(display, winder, EV_MASK);

	theDrawable = (Drawable)winder;

	/*
	 *  Let VOGLE know about the window size.
	 */
	vdevice.sizeX = vdevice.sizeY = MIN(h, w) - 1;
	vdevice.sizeSx = w - 1;
	vdevice.sizeSy = h - 1;

	/*
	 * Create Graphics Context and Drawable
	 */
	theGC = XDefaultGC(display, screen);
	theDrawable = (Drawable)winder;
	DECX11_color(0);

	XMapRaised(display, winder);
	XFlush(display);

	/*
	 * Wait for Exposure event.
	 */
	do {
		XNextEvent(display, &event);
	} while (event.type != Expose);

	if (prefx == -1 && prefxs == -1)
		XSetInputFocus(display, winder, RevertToParent, CurrentTime);

	back_used = 0;

	return(1);
}

/*
 * DECX11_exit
 *
 *	cleans up before returning the window to normal.
 */
DECX11_exit()
{
	XFreeGC(display, theGC);

	if (back_used) 
		XFreePixmap(display, bbuf);

	XUnmapWindow(display, winder);

	XDestroyWindow(display, winder);

	return(1);
}

/*
 * DECX11_draw
 *
 *	draws a line from the current graphics position to (x, y).
 *
 * Note: (0, 0) is defined as the top left of the window in X (easy
 * to forget).
 */
DECX11_draw(x, y)
	int	x, y;
{
	XDrawLine(display,
		theDrawable,
		theGC,
		vdevice.cpVx, vdevice.sizeSy - vdevice.cpVy,
		x, vdevice.sizeSy - y
	);

	XFlush(display);
}

/*
 * DECX11_getkey
 *
 *	grab a character from the keyboard - blocks until one is there.
 */
int
DECX11_getkey()
{
	char	c;

	do {
		XNextEvent(display, &event);
		if (event.type == KeyPress) {
			if (XLookupString(&event, &c, 1, NULL, NULL) > 0)
				return((int)c);
			else
				return(0);
		}
	} while (event.type != KeyPress);
}

/*
 * DECX11_checkkey
 *
 *	Check if there has been a keyboard key pressed.
 *	and return it if there is.
 */
int
DECX11_checkkey()
{
	char	c;
	int	i;

	if (!XCheckWindowEvent(display, winder, KeyPressMask, &event))
		return(0);

	if (event.type == KeyPress)
		if (XLookupString(&event, &c, 1, NULL, NULL) > 0)
			return((int)c);

	return(0);
}

/*
 * DECX11_locator
 *
 *	return the window location of the cursor, plus which mouse button,
 * if any, is been pressed.
 */
int
DECX11_locator(wx, wy)
	int	*wx, *wy;
{
	Window	rootw, childw;
	int	x, y, mask;

	XQueryPointer(display, winder, &rootw, &childw, &x, &y, wx, wy, &mask);

	*wy = (int)vdevice.sizeSy - *wy;

	return(mask >> 8);
}

/*
 * DECX11_clear
 *
 * Clear the screen (or current buffer )to current colour
 */
DECX11_clear()
{
	XSetBackground(display, theGC, colour);
	XFillRectangle(display,
		theDrawable,
		theGC,
		vdevice.minVx,
		vdevice.minVy, 
		(unsigned int)vdevice.maxVx,
		(unsigned int)vdevice.maxVy
	);
}

/*
 * DECX11_color
 *
 *	set the current drawing color index.
 */
DECX11_color(ind)
        int	ind;
{
	colour = carray[ind];
	XSetForeground(display, theGC, colour);
}

/*
 * DECX11_mapcolor
 *
 *	change index i in the color map to the appropriate r, g, b, value.
 */
DECX11_mapcolor(i, r, g, b)
	int	i;
	int	r, g, b;
{
	int	stat;
	XColor	tmp;

	if (i >= CMAPSIZE)
		return(-1);


	/*
	 * For Black and White.
	 * If the index is 0 and r,g,b != 0 then we are remapping black.
	 * If the index != 0 and r,g,b == 0 then we make it black.
	 */
	if (vdevice.depth == 1) {
		if (i == 0 && (r != 0 || g != 0 || b != 0)) 
			carray[i] = WhitePixel(display, screen);
		else if (i != 0 && r == 0 && g == 0 && b == 0)
			carray[i] = BlackPixel(display, screen);
	} else {
		tmp.red = (unsigned short)(r / 255.0 * 65535);
		tmp.green = (unsigned short)(g / 255.0 * 65535);
		tmp.blue = (unsigned short)(b / 255.0 * 65535);
		tmp.flags = 0;
		tmp.pixel = (unsigned long)i;

		if ((stat = XAllocColor(display, colormap, &tmp)) == 0) {
			fprintf(stderr, "XAllocColor failed (status = %d)\n", stat);
			exit(1);
		}
		carray[i] = tmp.pixel;
	}

	XFlush(display);
	return(0);
}
	
/*
 * DECX11_font
 *
 *   Set up a hardware font. Return 1 on success 0 otherwise.
 *
 */
DECX11_font(fontfile)
        char	*fontfile;
{
	XGCValues	xgcvals;

	if (font_id != (XFontStruct *)NULL)
		XFreeFont(display, font_id);

	if (strcmp(fontfile, "small") == 0) {
		if ((font_id = XLoadQueryFont(display, SMALLX11R2)) == (XFontStruct *)NULL) {		/* X11 R2 */
			if ((font_id = XLoadQueryFont(display, SMALLX11R3)) == (XFontStruct *)NULL)	 	/* X11 R3 */
				return(0);
			else
				fontfile = SMALLX11R3;
		} else
			fontfile = SMALLX11R2;
	} else if (strcmp(fontfile, "large") == 0) {
		if ((font_id = XLoadQueryFont(display, LARGEX11R2)) == (XFontStruct *)NULL) {		/* X11 R2 */
			if ((font_id = XLoadQueryFont(display, LARGEX11R3)) == (XFontStruct *)NULL)	 	/* X11 R3 */
				return(0);
			else
				fontfile = LARGEX11R3;
		} else
			fontfile = LARGEX11R2;
	} else if ((font_id = XLoadQueryFont(display, fontfile)) == (XFontStruct *)NULL)
		return(0);

	vdevice.hheight = font_id->max_bounds.ascent + font_id->max_bounds.descent;
	vdevice.hwidth = font_id->max_bounds.width;

	xgcvals.font = XLoadFont(display, fontfile);
	XChangeGC(display, theGC, GCFont, &xgcvals);

	return(1);
}

/* 
 * DECX11_char
 *
 *	 outputs one char - is more complicated for other devices
 */
DECX11_char(c)
	char	c;
{
	char	*s = " ";

	s[0] = c;
	XDrawString(display, theDrawable, theGC, vdevice.cpVx, (int)(vdevice.sizeSy - vdevice.cpVy), s, 1);
	XFlush(display);
}

/*
 * DECX11_string
 *
 *	Display a string at the current drawing position.
 */
DECX11_string(s)
        char	s[];
{
	XDrawString(display, theDrawable, theGC, vdevice.cpVx, (int)(vdevice.sizeSy - vdevice.cpVy), s, strlen(s));
	XSync(display, 0);
}

/*
 * DECX11_fill
 *
 *	fill a polygon
 */
DECX11_fill(n, x, y)
	int	n, x[], y[];
{
	XPoint	plist[128];
	int	i;

	if (n > 128)
		verror("vogl: more than 128 points in a polygon");

	for (i = 0; i < n; i++) {
		plist[i].x = x[i];
		plist[i].y = vdevice.sizeSy - y[i];
	}

	XFillPolygon(display, theDrawable, theGC, plist, n, Nonconvex, CoordModeOrigin);

	vdevice.cpVx = x[n-1];
	vdevice.cpVy = y[n-1];

	XFlush(display);
}

#define	GC_COPY_MASK	~0

/*
 * DECX11_backbuf
 *
 *	Set up double buffering by allocating the back buffer and
 *	setting drawing into it.
 */
DECX11_backbuf()
{
	if (!back_used)
		bbuf = XCreatePixmap(display,
			(Drawable)winder,
			(unsigned int)vdevice.sizeSx + 1,
			(unsigned int)vdevice.sizeSy + 1,
			(unsigned int)vdevice.depth
		);

	theDrawable = (Drawable)bbuf;

	back_used = 1;

	return(1);
}

/*
 * DECX11_swapbuf
 *
 *	Swap the back and from buffers. (Really, just copy the
 *	back buffer to the screen).
 */
DECX11_swapbuf()
{
	XCopyArea(display,
		theDrawable,
		winder,
		theGC,
		0, 0,
		(unsigned int)vdevice.sizeSx + 1,
		(unsigned int)vdevice.sizeSy + 1,
		0, 0
	);

	XSync(display, 0);	/* Not XFlush */
}

/*
 * DECX11_frontbuf
 *
 *	Make sure we draw to the screen.
 */
DECX11_frontbuf()
{
	theDrawable = (Drawable)winder;
}

/*
 * the device entry
 */
static DevEntry DECX11dev = {
	"decX11",
	"large",
	"small",
	DECX11_backbuf,
	DECX11_char,
	DECX11_checkkey,
	DECX11_clear,
	DECX11_color,
	DECX11_draw,
	DECX11_exit,
	DECX11_fill,
	DECX11_font,
	DECX11_frontbuf,
	DECX11_getkey,
	DECX11_init,
	DECX11_locator,
	DECX11_mapcolor,
	noop,
	noop,
	DECX11_string,
	DECX11_swapbuf,
	noop
};

/*
 * _DECX11_devcpy
 *
 *	copy the decX11 device into vdevice.dev.
 */
_DECX11_devcpy()
{
	vdevice.dev = DECX11dev;
}

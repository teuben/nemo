/*
 * VOGL/VOGLE driver for X11.
 * 
 * Define VOGLE if this driver is really for the VOGLE Libarary.
 *
 */
#undef VOGLE

#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#ifdef VOGLE

#include "vogle.h"
static	char	*me = "vogle";
#define LARGEFONT       "-adobe-courier-medium-r-normal--24-240-75-75-m-150-iso8859-1"
#define SMALLFONT       "-adobe-courier-medium-r-normal--10-100-75-75-m-60-iso8859-1"

#else

#include "vogl.h"
static	char	*me = "vogl";
#define LARGEFONT	"9x15bold"
#define SMALLFONT	"6x13bold"

#endif

#define MIN(x,y)	((x) < (y) ? (x) : (y))
#define MAX(x,y)	((x) > (y) ? (x) : (y))
#define	CMAPSIZE	256
#define	EV_MASK		KeyPressMask|ButtonReleaseMask|ExposureMask|ButtonPressMask

static	int		maxw = -1, maxh = -1;
static	Window		winder;
static	Display		*display;
static	int		screen;
static	unsigned long	carray[CMAPSIZE];
static	Colormap	colormap;

static	Drawable	theDrawable = 0xffffffff;	/* (unsigned)-1 */
static	GC		theGC;
static	XGCValues	theGCvalues;
static	Pixmap		bbuf;		/* Back buffer pixmap */
static	int		back_used = 0;	/* Have we backbuffered ? */

static	XFontStruct	*font_id = (XFontStruct *)NULL;
XEvent			event;

static	unsigned long	colour;
static	unsigned int	h, w;
static	char		*smallf, *largef;
static	char		use_toolkit_win = 0;

/*
 * vo_xt_set_win
 *
 *	Just sets the drawable to the partucular window.
 */
int
vo_xt_set_win(dis, win, xw, xh)
	Display		*dis;
	Drawable	win;
	int		xw, xh;
{
	int	backb;

	backb = (theDrawable == bbuf);

	winder = win;

	vdevice.sizeX = vdevice.sizeY = MIN(xh, xw);
	vdevice.sizeSx = xw;
	vdevice.sizeSy = xh;

        if (xw > maxw || xh > maxh) {
		if (back_used) {
			back_used = 0;
			XFreePixmap(display, bbuf);
			X11_backbuf();
		}
        }

	display = dis;
	if (backb)
		theDrawable = bbuf;
	else
		theDrawable = win;

	return(1);
}

/*
 * vo_xt_window
 *
 *	Tells VOGL/VOGLE to use a window from an X11 toolkit (eg xview)
 *	and not to make it's own window.
 */
int
vo_xt_window(dis, win, xw, xh)
	Display	*dis;
	Window	win;
	int	xw, xh;
{
	int	backb, i, depth;

	backb = (theDrawable == bbuf);

	display = dis;
	winder = win;
	screen = DefaultScreen(display);
	colormap = DefaultColormap(display, screen);
	depth = vdevice.depth = DefaultDepth(display, screen);
	theDrawable = winder;

	use_toolkit_win = 1;
	w = xw;
	h = xh;

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
		X11_mapcolor(0, 0, 0, 0);
		X11_mapcolor(1, 255, 0, 0);
		X11_mapcolor(2, 0, 255, 0);
		X11_mapcolor(3, 255, 255, 0);
		X11_mapcolor(4, 0, 0, 255);
		X11_mapcolor(5, 255, 0, 255);
		X11_mapcolor(6, 0, 255, 255);
		X11_mapcolor(7, 255, 255, 255);
	}

	if ((smallf = XGetDefault(display, me, "smallfont")) == (char *)NULL)
		smallf = SMALLFONT;

	if ((largef = XGetDefault(display, me, "largefont")) == (char *)NULL)
		largef = LARGEFONT;

	/*
	 * Create Graphics Context and Drawable
	 */
	theGC = XDefaultGC(display, screen);
	theGCvalues.graphics_exposures = False;
	theGCvalues.cap_style = CapButt;
	XChangeGC(display, theGC, GCGraphicsExposures|GCCapStyle, &theGCvalues);
	X11_color(0);

	vdevice.sizeX = vdevice.sizeY = MIN(xh, xw);
	vdevice.sizeSx = xw;
	vdevice.sizeSy = xh;

        if (back_used && (xw > maxw || xh > maxh)) {
                back_used = 0;
		XFreePixmap(display, bbuf);
                X11_backbuf();
        }

	if (backb)
		theDrawable = bbuf;
	else
		theDrawable = win;


#ifndef VOGLE
	vdevice.devname = "X11";
#endif

	return(1);
}

/*
 *	vo_xt_win_size
 *
 * If the X toolkit has changed the window size, then
 * you might wish to call this routine to tell vogl/vogle about it.
 */
void
vo_xt_win_size(xw, xh)
	int	xw, xh;
{
	char	backb;

	w = xw;
	h = xh;

	vdevice.sizeX = vdevice.sizeY = MIN(h, w);
	vdevice.sizeSx = w;
	vdevice.sizeSy = h;

	backb = (theDrawable == bbuf);

	if (back_used) {

		/* Have to re allocate the back buffer */

		XFreePixmap(display, bbuf);

		bbuf = XCreatePixmap(display,
			(Drawable)winder,
			(unsigned int)vdevice.sizeSx,
			(unsigned int)vdevice.sizeSy,
			(unsigned int)vdevice.depth
		);
	}
	if (backb)
		theDrawable = (Drawable)bbuf;
}

/*
 * return the X display in use.
 */
Display *
vo_xt_get_display()
{
	return(display);
}

/*
 * return the X Window in use.
 */
Window
vo_xt_get_window()
{
	return(winder);
}

/*
 * return the Graphics Context in use.
 */
GC
vo_xt_get_GC()
{
	return(theGC);
}

/*
 * Set the Graphics Context to use.
 */
void
vo_xt_set_GC(gc)
	GC	gc;
{
	theGC = gc;
}


/*
 * X11_init
 *
 *	initialises X11 display.
 */
X11_init()
{
	int		i;
	int		x, y, prefx, prefy, prefxs, prefys;
	unsigned int	bw, depth, mask;
	Window		rootw, childw;
	char		*av[2], name[128], *geom;

	XSetWindowAttributes    theWindowAttributes;
	XWindowAttributes	retWindowAttributes;
        XSizeHints              theSizeHints;
        unsigned long           theWindowMask;
	XWMHints                theWMHints;


	if (use_toolkit_win)
		return(1);

	av[0] = me;
	av[1] = (char *)NULL;

	if ((display = XOpenDisplay((char *)NULL)) == (Display *)NULL) {
		fprintf(stderr,"%s: X11_init: can't connect to X server\n", me);
		exit(1);
	}

	screen = DefaultScreen(display);
	winder = RootWindow(display, screen);
	colormap = DefaultColormap(display, screen);
	depth = vdevice.depth = DefaultDepth(display, screen);

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
		X11_mapcolor(0, 0, 0, 0);
		X11_mapcolor(1, 255, 0, 0);
		X11_mapcolor(2, 0, 255, 0);
		X11_mapcolor(3, 255, 255, 0);
		X11_mapcolor(4, 0, 0, 255);
		X11_mapcolor(5, 255, 0, 255);
		X11_mapcolor(6, 0, 255, 255);
		X11_mapcolor(7, 255, 255, 255);
	}

	getprefposandsize(&prefx, &prefy, &prefxs, &prefys);

	/*
	 * NEED TO USE XGRABPOINTER here???
	 */
	XQueryPointer(display, winder, &rootw, &childw, &x, &y, &x, &y, &mask);

	if (childw == None)
		childw = rootw;

/*
	if (!XGetWindowAttributes(display, childw, &retWindowAttributes)) {
		fprintf(stderr,"Can't get window attributes.");
		exit(1);
	}

	x = retWindowAttributes.x;
	y = retWindowAttributes.y;
	w = retWindowAttributes.width;
	h = retWindowAttributes.height;
	bw = retWindowAttributes.border_width;
	depth = vdevice.depth = retWindowAttributes.depth;

	XTranslateCoordinates(display,
			childw, retWindowAttributes.root,
			0, 0,
			&x, &y,
			&rootw
	);
*/

	XGetGeometry(display, childw, &rootw, &x, &y, &w, &h, &bw, &depth);

        theWindowAttributes.backing_store = WhenMapped;
        theWindowAttributes.save_under = True;
        theWindowAttributes.border_pixel = carray[1];


	/*
	 * See if there is something in the .Xdefaults file regarding
	 * VOGL/VOGLE.
	 */

	if ((smallf = XGetDefault(display, me, "smallfont")) == (char *)NULL)
		smallf = SMALLFONT;

	if ((largef = XGetDefault(display, me, "largefont")) == (char *)NULL)
		largef = LARGEFONT;

	geom = XGetDefault(display, me, "Geometry");

	theSizeHints.flags = PPosition | PSize;

	if (geom != (char *)NULL) {

		theSizeHints.flags = 0;

		mask = XParseGeometry(geom, &x, &y, &w, &h);

		if (mask & XValue)
			theSizeHints.flags |= USPosition;

		if (mask & YValue)
			theSizeHints.flags |= USPosition;

		if (mask & WidthValue)
			theSizeHints.flags |= USSize;

		if (mask & HeightValue)
			theSizeHints.flags |= USSize;

		if (mask & XNegative)
			 x = DisplayWidth(display, screen) - 2*bw - w + x;

		if (mask & YNegative)
			y = DisplayHeight(display, screen) - 2*bw - h + y;
	}

	if (prefx > -1) {
	        x = prefx;
	        y = prefy;
	}

	if (prefxs > -1) {
	        w = prefxs;
	        h = prefys;
	}

	if (bw == 0)
		bw = 4;

	x -= bw;
	y -= bw;

	if (x <= 0)
		x = 0;

	if (y <= 0)
		y = 0;

	w -= 4 * bw;
	h -= 4 * bw;

        theWindowMask = CWBorderPixel|CWBackingStore;

        winder = XCreateWindow(display,
                                winder,
                                x, y,
                                w, h,
                                bw,
                                (int)vdevice.depth,
                                InputOutput,
                                CopyFromParent,
                                theWindowMask,
                                &theWindowAttributes
                        );

	XSetWindowColormap(display, winder, colormap);
 
        theSizeHints.x = x;
        theSizeHints.y = y;
        theSizeHints.width = w;
        theSizeHints.height = h;

#ifndef VOGLE
	if (vdevice.wintitle)
		strcpy(name, vdevice.wintitle);
	else
		sprintf(name, "%s %d (win id 0x%x)", me, getpid(), winder);
#else
	sprintf(name, "%s %d (win id 0x%x)", me, getpid(), winder);
#endif

	XSetStandardProperties(display,
		winder,
		name,
		name,
		None,
		av,
		1,
		&theSizeHints
	);

        theWMHints.initial_state = NormalState;
        theWMHints.input = True;
        theWMHints.flags = StateHint | InputHint;
        XSetWMHints(display, winder, &theWMHints);

	XSelectInput(display, winder, EV_MASK);

	theDrawable = (Drawable)winder;

	/*
	 * Create Graphics Context and Drawable
	 */
	theGC = XDefaultGC(display, screen);
	theGCvalues.graphics_exposures = False;
	theGCvalues.cap_style = CapButt;
	XChangeGC(display, theGC, GCGraphicsExposures|GCCapStyle, &theGCvalues);
	theDrawable = (Drawable)winder;
	X11_color(0);

	XMapRaised(display, winder);
	XFlush(display);

	/*
	 * Wait for Exposure event.
	do {
		XNextEvent(display, &event);
	} while (event.type != Expose);
	 */
	XWindowEvent(display, winder, ExposureMask, &event);

	/*
	 * Set the input Focus to us.

        if (prefx == -1 && prefxs == -1)
                XSetInputFocus(display, winder, RevertToParent, CurrentTime);
	 */

	/*
	 *  Let VOGL/VOGLE know about the window size.
	 *  (We may have been resized..... )
	 */
	if (!XGetWindowAttributes(display, winder, &retWindowAttributes)) {
		fprintf(stderr,"Can't get window attributes.");
		exit(1);
	}

	x = retWindowAttributes.x;
	y = retWindowAttributes.y;
	w = retWindowAttributes.width;
	h = retWindowAttributes.height;

	XTranslateCoordinates(display,
			winder, retWindowAttributes.root,
			0, 0,
			&x, &y,
			&rootw
	);

	vdevice.sizeX = vdevice.sizeY = MIN(h, w);
	vdevice.sizeSx = w;
	vdevice.sizeSy = h;

	if (back_used && (maxw < w || maxh < h)) {
		back_used = 0;
		XFreePixmap(display, bbuf);
		X11_backbuf();
	}

	return(1);
}

/*
 * X11_exit
 *
 *	cleans up before returning the window to normal.
 */
X11_exit()
{
	if (back_used || bbuf != 0xffffffff) 
		XFreePixmap(display, bbuf);

	back_used = 0;
	bbuf = 0xffffffff;
	maxw = maxh = -1;

	if (font_id != (XFontStruct *)NULL)
		XFreeFont(display, font_id);

	font_id = (XFontStruct *)NULL;

#ifdef NEWCMAP
	if (colormap != DefaultColormap(display, screen))
		XFreeColormap(display, colormap);

	colormap = 0;
#endif

	if (use_toolkit_win) {
		use_toolkit_win = 0;
		return(1);
	}

	XDestroyWindow(display, winder);

	XSync(display, 0);

	XCloseDisplay(display);

	display = (Display *)NULL;
	winder = 0;

	return(1);
}

/*
 * X11_draw
 *
 *	draws a line from the current graphics position to (x, y).
 *
 * Note: (0, 0) is defined as the top left of the window in X (easy
 * to forget).
 */
X11_draw(x, y)
	int	x, y;
{
	if (x == vdevice.cpVx && y == vdevice.cpVy)
		/*
		 * Hack for some X servers... my MIT X11 R5 manual states:
		 * CapButt	the results are device-dependent, but the
		 *		desired effect is that a single pixel is
		 *		drawn.
		 *
		 * This does work on a real MIT R5 server... but on some
		 * machines, we have to do this XDrawPoint thing.
		 * (It's probably faster this way anyway).
		 */
		XDrawPoint(display, theDrawable, theGC, x, vdevice.sizeSy - y);
	else
		XDrawLine(display,
			theDrawable,
			theGC,
			vdevice.cpVx, vdevice.sizeSy - vdevice.cpVy,
			x, vdevice.sizeSy - y
		);

	if (vdevice.sync)
		XSync(display, 0);
}

X11_pnt(x, y)
	int	x, y;
{
	XDrawPoint(display,
		theDrawable,
		theGC,
		x, vdevice.sizeSy - y
	);

	if (vdevice.sync)
		XSync(display, 0);
}

/*
 * X11_getkey
 *
 *	grab a character from the keyboard - blocks until one is there.
 */
int
X11_getkey()
{
	char	c;

	do {
		XNextEvent(display, &event);
		if (event.type == KeyPress) {
			if (XLookupString((XKeyEvent *)&event, &c, 1, NULL, NULL) > 0)
				return((int)c);
			else
				return(0);
		}
	} while (event.type != KeyPress);
}

/*
 * X11_checkkey
 *
 *	Check if there has been a keyboard key pressed.
 *	and return it if there is.
 */
int
X11_checkkey()
{
	char	c;

	if (!XCheckWindowEvent(display, winder, KeyPressMask, &event))
		return(0);

	if (event.type == KeyPress)
		if (XLookupString((XKeyEvent *)&event, &c, 1, NULL, NULL) > 0)
			return((int)c);

	return(0);
}

/*
 * X11_locator
 *
 *	return the window location of the cursor, plus which mouse button,
 * if any, is been pressed.
 */
int
X11_locator(wx, wy)
	int	*wx, *wy;
{
	Window		rootw, childw;
	int		x, y;
	unsigned int	mask;

	XQueryPointer(display, winder, &rootw, &childw, &x, &y, wx, wy, &mask);

	*wy = (int)vdevice.sizeSy - *wy;

	return(mask >> 8);
}

#ifdef VOGLE
/*
 * X11_clear
 *
 * Clear the screen (or current buffer )to current colour
 */
X11_clear()
{
	XSetBackground(display, theGC, colour);
	XFillRectangle(display,
		theDrawable,
		theGC,
		0,
		0,
		(unsigned int)vdevice.sizeSx,
		(unsigned int)vdevice.sizeSy
	);

	if (vdevice.sync)
		XFlush(display);
}

#else 

/*
 * X11_clear
 *
 * Clear the screen (or current buffer )to current colour
 */
X11_clear()
{
	unsigned int	w = vdevice.maxVx - vdevice.minVx;
	unsigned int	h = vdevice.maxVy - vdevice.minVy;

	XSetBackground(display, theGC, colour);

	XFillRectangle(display,
		theDrawable,
		theGC,
		vdevice.minVx,
		vdevice.sizeSy - vdevice.maxVy - 1, 
		w + 1, 
		h + 1
	);

	if (vdevice.sync)
		XFlush(display);
}
#endif

/*
 * X11_color
 *
 *	set the current drawing color index.
 */
X11_color(ind)
        int	ind;
{
	colour = carray[ind];
	XSetForeground(display, theGC, colour);
}

/*
 * X11_mapcolor
 *
 *	change index i in the color map to the appropriate r, g, b, value.
 */
X11_mapcolor(i, r, g, b)
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
#ifdef NEWCMAP
			colormap = XCopyColormapAndFree(display, colormap);
			XSetWindowColormap(display, winder, colormap);
#else
			fprintf(stderr, "XAllocColor failed (status = %d)\n", stat);

			exit(1);
#endif
		}
		carray[i] = tmp.pixel;
	}

	XFlush(display);
	return(0);
}
	
/*
 * X11_font
 *
 *   Set up a hardware font. Return 1 on success 0 otherwise.
 *
 */
X11_font(fontfile)
        char	*fontfile;
{
	XGCValues	xgcvals;
	char	*name = fontfile;

	if (font_id != (XFontStruct *)NULL)
		XFreeFont(display, font_id);

	if (strcmp(fontfile, "small") == 0) {
		if ((font_id = XLoadQueryFont(display, smallf)) == (XFontStruct *)NULL) {
			fprintf(stderr, "%s X11.c couldn't open small font '%s'\n", me, smallf);
			fprintf(stderr, "You'll have to redefine it....\n");
			return(0);
		} else
			name = smallf;
		
	} else if (strcmp(fontfile, "large") == 0) {
		if ((font_id = XLoadQueryFont(display, largef)) == (XFontStruct *)NULL) {
			fprintf(stderr, "%s X11.c couldn't open large font '%s'\n", me, largef);
			fprintf(stderr, "You'll have to redefine it....\n");
			return(0);
		}
			name = largef;
	} else {
		if ((font_id = XLoadQueryFont(display, fontfile)) == (XFontStruct *)NULL) {
			fprintf(stderr, "%s X11.c couldn't open fontfile '%s'\n", me, fontfile);
			return(0);
		}
	}

	/*
	vdevice.hheight = font_id->max_bounds.ascent + font_id->max_bounds.descent;
	*/
	vdevice.hheight = font_id->ascent + font_id->descent;
	vdevice.hwidth = font_id->max_bounds.width;


	xgcvals.font = XLoadFont(display, name);
	XChangeGC(display, theGC, GCFont, &xgcvals);

	return(1);
}

/* 
 * X11_char
 *
 *	 outputs one char - is more complicated for other devices
 */
X11_char(c)
	char	c;
{
	XDrawString(display, theDrawable, theGC, vdevice.cpVx, (int)(vdevice.sizeSy - vdevice.cpVy), &c, 1);

	if (vdevice.sync)
		XFlush(display);
}

/*
 * X11_string
 *
 *	Display a string at the current drawing position.
 */
X11_string(s)
        char	s[];
{
	XDrawString(display, theDrawable, theGC, vdevice.cpVx, (int)(vdevice.sizeSy - vdevice.cpVy), s, strlen(s));
	if (vdevice.sync)
		XFlush(display);
}

/*
 * X11_fill
 *
 *	fill a polygon
 */
X11_fill(n, x, y)
	int	n, x[], y[];
{
	char	buf[BUFSIZ];
	XPoint	plist[128];
	int	i;

	if (n > 128) {
		sprintf(buf, "%s: more than 128 points in a polygon", me);
		verror(buf);
	}

	for (i = 0; i < n; i++) {
		plist[i].x = x[i];
		plist[i].y = vdevice.sizeSy - y[i];
	}

	XFillPolygon(display, theDrawable, theGC, plist, n, Nonconvex, CoordModeOrigin);

	vdevice.cpVx = x[n-1];
	vdevice.cpVy = y[n-1];

	if (vdevice.sync)
		XFlush(display);
}

/*
 * X11_backbuf
 *
 *	Set up double buffering by allocating the back buffer and
 *	setting drawing into it.
 */
int
X11_backbuf()
{
	if (!back_used) {
		bbuf = XCreatePixmap(display,
			(Drawable)winder,
			(unsigned int)vdevice.sizeSx,
			(unsigned int)vdevice.sizeSy,
			(unsigned int)vdevice.depth
		);

		maxw = MAX(vdevice.sizeSx, maxw);
		maxh = MAX(vdevice.sizeSy, maxh);
	}

	theDrawable = (Drawable)bbuf;

	back_used = 1;

	return(1);
}

/*
 * X11_swapbuf
 *
 *	Swap the back and from buffers. (Really, just copy the
 *	back buffer to the screen).
 */
X11_swapbuf()
{
	XCopyArea(display,
		theDrawable,
		winder,
		theGC,
		0, 0,
		(unsigned int)vdevice.sizeSx,
		(unsigned int)vdevice.sizeSy,
		0, 0
	);
	XSync(display, 0);
}

/*
 * X11_frontbuf
 *
 *	Make sure we draw to the screen.
 */
X11_frontbuf()
{
	theDrawable = (Drawable)winder;
}

/*
 * Syncronise the display with what we think has been sent to it...
 */
X11_sync()
{
	XSync(display, 0);
}

#undef VORTDUMP
#ifdef VORTDUMP
/*
 * HACK
 * Dump the contents of the current buffer to a VORT file....
 * ONLY WORKS WITH 8Bit Drawables!
 */
#include "vort.h"

X11_dump_pixmap(filename, dx, dy, dw, dh)
	char	*filename;
	int	dx, dy, dw, dh;
{
	XImage	*ximage;
	image	*im;
	unsigned char	*line, *rm, *gm, *bm;
	XColor	*cols;
	int	i;

	if (dw > vdevice.sizeSx || dw < 0)
		dw = vdevice.sizeSx;
	if (dh > vdevice.sizeSy || dh < 0)
		dh = vdevice.sizeSy;

	if (dx > vdevice.sizeSx || dx < 0)
		dx = 0;
	if (dy > vdevice.sizeSy || dy < 0)
		dy = 0;

	ximage = XGetImage(display, 
			theDrawable, 
			dx, dy,
			(unsigned int)dw,
			(unsigned int)dh,
			AllPlanes,
			ZPixmap
		);

	if (!ximage) {
		fprintf(stderr, "X11_dump_pixmap: can't do XGetImage\n");
		exit(1);
	}

	if ((im = openimage(filename, "w")) == (image *)NULL) {
		fprintf(stderr, "X11_dump_pixmap: can't open %s\n", filename);
		exit(1);
	}

	if (!(rm = (unsigned char *)malloc(256))) {
		fprintf(stderr, "X11_dump_pixmap: can't alloc rm\n");
		exit(1);
	}
	if (!(gm = (unsigned char *)malloc(256))) {
		fprintf(stderr, "X11_dump_pixmap: can't alloc gm\n");
		exit(1);
	}
	if (!(bm = (unsigned char *)malloc(256))) {
		fprintf(stderr, "X11_dump_pixmap: can't alloc bm\n");
		exit(1);
	}
	if (!(cols = (XColor *)malloc(256 * sizeof(XColor)))) {
		fprintf(stderr, "X11_dump_pixmap: can't alloc cols\n");
		exit(1);
	}

	/*
	 * Get our colormap...
	 */
	for (i = 0; i < 256; i++) {
		cols[i].pixel = (unsigned long)i;
		cols[i].red = cols[i].green = cols[i].blue = 0;
		cols[i].flags = DoRed | DoGreen | DoBlue;
	}

	XQueryColors(display, colormap, cols, 256);

	for (i = 0; i < 256; i++) {
		rm[i] = (unsigned char)(cols[i].red >> 8);
		gm[i] = (unsigned char)(cols[i].green >> 8);
		bm[i] = (unsigned char)(cols[i].blue >> 8);
	}

	imagetype(im) = PIX_RLECMAP;
	imageheight(im) = dh;
	imagewidth(im) = dw;
	imagedate(im) = time(0);
	titlelength(im) = 0;
	setcmap(im, 256, rm, gm, bm);

	writeheader(im);

	line = (unsigned char *)ximage->data;
	for (i = 0; i < dh; i++) {
		writemappedline(im, line);
		line += ximage->bytes_per_line;
	}
	
	closeimage(im); 

	free(rm);
	free(gm);
	free(bm);
	free(cols);
}

#endif

#ifndef VOGLE
/*
 * X11_setlw
 *
 *	Set the line width....
 */
X11_setlw(w)
	int	w;
{
	XGCValues vals;

	vals.line_width = w;
	XChangeGC(display, theGC, GCLineWidth, &vals);
}

/*
 * X11_setls
 *
 *	Set the line style....
 */

X11_setls(lss)
	int	lss;
{
	unsigned ls = lss;
	char	dashes[16];
	int	i, n, a, b, offset;

	if (ls == 0xffff) {
		XSetLineAttributes(display, theGC, vdevice.attr->a.lw, LineSolid, CapButt, JoinMiter);
		return;
	}

	for (i = 0; i < 16; i++)
		dashes[i] = 0;

	for (i = 0; i < 16; i++)	/* Over 16 bits */
		if ((ls & (1 << i)))
			break;

	offset = i;

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
	n++;

	XSetLineAttributes(display, theGC, vdevice.attr->a.lw, LineOnOffDash, CapButt, JoinMiter);
	XSetDashes(display, theGC, offset, dashes, n);
}

#else
/*
 * X11_setlw (this one for VOGLE only)
 *
 *	Set the line width....THICK or THIN
 */
X11_setlw(w)
	int	w;
{
	XGCValues vals;

	if (w == 0)
		w = 1;
	else if (w == 1)
		w = 2;

	vals.line_width = w;
	XChangeGC(display, theGC, GCLineWidth, &vals);
}

#endif 

/*
 * the device entry
 */
static DevEntry X11dev = {
	"X11",
	"large",
	"small",
	X11_backbuf,
	X11_char,
	X11_checkkey,
	X11_clear,
	X11_color,
	X11_draw,
	X11_exit,
	X11_fill,
	X11_font,
	X11_frontbuf,
	X11_getkey,
	X11_init,
	X11_locator,
	X11_mapcolor,
#ifndef VOGLE
	X11_setls,
#endif
	X11_setlw,
	X11_string,
	X11_swapbuf,
	X11_sync
};

/*
 * _X11_devcpy
 *
 *	copy the X11 device into vdevice.dev.
 */
_X11_devcpy()
{
	vdevice.dev = X11dev;
}

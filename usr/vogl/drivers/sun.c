
/*
 * Vogl/Vogle driver for Sun using sunview.
 *
 */
/*
 * define VOGLE for VOGLE library or leave blank for VOGL library
#define VOGLE 1
 */
#undef VOGLE

#include <stdio.h>

#include <suntool/sunview.h>
#include <suntool/canvas.h>
#include <fcntl.h>
#include <errno.h>

#ifdef SUN_3_5
#define event_action	event_id
#endif

#ifdef VOGLE
#include "vogle.h"
#else
#include "vogl.h"
#endif

#define	CMAPSIZE	256
#define	STDFONTDIR	"/usr/lib/fonts/fixedwidthfonts/"

#define MIN(x,y)	((x) < (y) ? (x) : (y))

#define OP_WHITE	(PIX_SRC | PIX_COLOR(7) | PIX_DONTCLIP)
#define OP_BLACK	(PIX_SRC | PIX_COLOR(0) | PIX_DONTCLIP)
#define COL_WHITE	7
#define COL_BLACK	0

static Pixwin	*pw_tmp, *pw;

#ifdef BASE_COL_SET
static Pixwin	*pw0;
#endif

static Pixrect	*backb;
static Pixfont	*font_id;
static int	wfd, blanket_win;
static int	oldflags;
static int	pwfd, h, w;
static int	colour;
static Rect	wrect;

static char	use_sunview_canvas = 0;

#define	DASHED		1
#define	FATLINES	2
static unsigned char	lineflags = 0;

static Pr_brush	brush = {1};
static Pr_texture	*t = (Pr_texture *)NULL;
static Pr_texture	tex = {
				(short *)NULL,	
				0,
				{1, 1, 1, 0},
				0,
			};

/*
 * default colour map
 */

static int	colnum = 8;
static u_char   red[CMAPSIZE] = {
			0, 255, 0, 255, 0, 255, 0, 255, 0,
		},
		green[CMAPSIZE] = {
			0, 0, 255, 255, 0, 0, 255, 255, 0,
		},
		blue[CMAPSIZE] = {
			0, 0, 0, 0, 255, 255, 255, 255, 0,
		};

/*
 * redisplay
 *
 *	redisplay the window.
 */
static void
redisplay()
{
	pw_damaged(pw);

	pw_repairretained(pw);

	pw_donedamaged(pw);
}

/*
 * To be called from a resize procedure from withing sunview
 */
vo_sunview_canvas_size(sw, sh)
	int	sw, sh;
{
	w = sw;
        h = sh;

        vdevice.sizeX = vdevice.sizeY = MIN(h, w);
        vdevice.sizeSx = w;
        vdevice.sizeSy = h;

	if (pw->pw_prretained) {
		mem_destroy(pw->pw_prretained);
		backb = pw->pw_prretained = mem_create(w, h, vdevice.depth);
	}
}

/*
 * vo_sunview_canvas
 *
 * 	Tells VOGLE/VOGL to use a supplied sunview pixwin
 */
vo_sunview_canvas(canvas, cw, ch)
	Canvas	canvas;
	int	cw, ch;
{
	pw = canvas_pixwin(canvas);
	use_sunview_canvas = 1;

	w = cw;
	h = ch;

	vdevice.sizeX = vdevice.sizeY = MIN(h, w);
	vdevice.sizeSx = w;
	vdevice.sizeSy = h;

	vdevice.depth = pw->pw_pixrect->pr_depth;
	if (!pw->pw_prretained)	/* Make us retained */
		backb = pw->pw_prretained = mem_create(w, h, vdevice.depth);

	/* 
	 *  Set up the color map.  
	 */

	if (vdevice.depth > 1) {
		pw_setcmsname(pw, "vogle");
		pw_putcolormap(pw, 0, colnum, red, green, blue);
	}

	wfd = (int)window_get(canvas, WIN_FD);

	/*
	 * Set non-blocking input for window.
	oldflags = fcntl(wfd, F_GETFL, 0);
	if (fcntl(wfd, F_SETFL, FNDELAY) < 0) {
		perror("F_SETFL");
		exit(1);
	}
	 */


	pw_batch_on(pw);

#ifndef VOGLE
	vdevice.devname = "sun";
#endif

	return(1);
}

/*
 * SUN_init
 *
 *	initialises drawing canvas to occupy current window
 */
SUN_init()
{
	int		i, prefx, prefy, prefxs, prefys, bw;
	char		name[WIN_NAMESIZE];
	Inputmask	imk, imp;
	int		rootfd;

	if (use_sunview_canvas)
		return(1);

	pw = (Pixwin *)NULL;

#ifdef BASE_COL_SET
	pw0 = pw;
#endif

	/*
	 * get the gfx win so we have some default sizes to use
	 */
	we_getgfxwindow(name);
	pwfd = open(name, 2);
	win_getrect(pwfd, &wrect);
	
        /*
         * Get the input masks of the base window...
         */
        win_get_pick_mask(pwfd, &imp);
        win_get_kbd_mask(pwfd, &imk);
 

	/*
	 * get a new window (either going to be a blanket window or
	 * a window in it's own right)
	 */
	if ((wfd = win_getnewwindow()) == -1) {
		fprintf(stderr, "No new windows!\n");
		exit(1);
	}

	getprefposandsize(&prefx, &prefy, &prefxs, &prefys);

	if (prefx > -1) {
		wrect.r_left = prefx;
		wrect.r_top = prefy;
	}

	if (prefxs > -1) {
		wrect.r_width = prefxs;
		wrect.r_height = prefys;
	}

	w = wrect.r_width;
	h = wrect.r_height;

	bw = 3;
	if (prefx > -1 || prefxs > -1) {
		/*
		 * Make room for a 3 pixel border
		 */
		if (wrect.r_left <= 2)
			wrect.r_left = 0;
		else
			wrect.r_left -= bw;

		if (wrect.r_top <= 2)
			wrect.r_top = 0;
		else
			wrect.r_top -= bw;

		wrect.r_width += 2 * bw;
		wrect.r_height += 2 * bw;

		win_setrect(wfd, &wrect);

		/*
		 * get the parent (probably full screen window)
		 * so we can size our window to any size we like
		 * on the screen.
		 */
		we_getparentwindow(name);
		rootfd = open(name, 2);

		win_setlink(wfd, WL_PARENT, win_fdtonumber(rootfd));
		win_setlink(wfd, WL_COVERED, WIN_NULLLINK);
		win_insert(wfd);

		wmgr_top(wfd, rootfd);

		pw = pw_open(wfd);

#ifdef BASE_COL_SET
		/*
		 * Get the pixrect for the window that we started in
		 * so we can set it's colourmap as well
		 */
		pw0 = pw_open(pwfd);
#endif

		close(rootfd);
		blanket_win = 0;
	} else {
		win_insertblanket(wfd, pwfd);
		pw = pw_region(pw_open(wfd), 0, 0, w, h);
		blanket_win = 1;
	}

	/*
	 * Set non-blocking input for window.
	 */
	oldflags = fcntl(wfd, F_GETFL, 0);
	if (fcntl(wfd, F_SETFL, FNDELAY) < 0) {
		perror("F_SETFL");
		exit(1);
	}

	/*
	 * Setup the input masks for window.
	 */

	win_set_kbd_mask(wfd, &imk);
        win_set_pick_mask(wfd, &imp);

	vdevice.depth = pw->pw_pixrect->pr_depth;

	/* 
	 *  Set up the color map.  
	 */

	if (vdevice.depth > 1) {
		pw_setcmsname(pw, "vogle");
		pw_putcolormap(pw, 0, colnum, red, green, blue);
#ifdef BASE_COL_SET
		if (pw0 != (Pixwin *)NULL) {
			pw_setcmsname(pw0, "vogle");
			pw_putcolormap(pw0, 0, colnum, red, green, blue);
		}
#endif
	}

	if (prefx > -1 || prefxs > -1) {
		/*
		 * Draw the border...
		 */
		int	x0, y0, x1, y1;

		x0 = y0 = 0;
		x1 = wrect.r_width - 1;
		y1 = 0;
		pw_vector(pw, x0, y0, x1, y1, OP_WHITE, COL_WHITE);
		pw_vector(pw, x0 + 2, y0 + 2, x1 - 2, y1 + 2, OP_WHITE, COL_WHITE);
		pw_vector(pw, x0 + 1, y0 + 1, x1 - 1, y1 + 1, OP_BLACK, COL_BLACK);
		x0 = x1;
		y0 = y1;
		x1 = x0;
		y1 = wrect.r_height - 1;
		pw_vector(pw, x0, y0, x1, y1, OP_WHITE, COL_WHITE);
		pw_vector(pw, x0 - 2, y0 + 2, x1 - 2, y1 - 2, OP_WHITE, COL_WHITE);
		pw_vector(pw, x0 - 1, y0 + 1, x1 - 1, y1 - 1, OP_BLACK, COL_BLACK);
		x0 = x1;
		y0 = y1;
		x1 = 0;
		y1 = wrect.r_height - 1;
		pw_vector(pw, x0, y0, x1, y1, OP_WHITE, COL_WHITE);
		pw_vector(pw, x0 - 2, y0 - 2, x1 + 2, y1 - 2, OP_WHITE, COL_WHITE);
		pw_vector(pw, x0 - 1, y0 - 1, x1 + 1, y1 - 1, OP_BLACK, COL_BLACK);
		x0 = x1;
		y0 = y1;
		x1 = 0;
		y1 = 0;
		pw_vector(pw, x0, y0, x1, y1, OP_WHITE, COL_WHITE);
		pw_vector(pw, x0 + 2, y0 - 2, x1 + 2, y1 + 2, OP_WHITE, COL_WHITE);
		pw_vector(pw, x0 + 1, y0 - 1, x1 + 1, y1 + 1, OP_BLACK, COL_BLACK);
		pw_tmp = pw;
		pw = pw_region(pw_tmp, 3, 3, w, h);
		pw_close(pw_tmp);
	}


	backb = pw->pw_prretained = mem_create(w, h, vdevice.depth);

	signal(SIGWINCH, redisplay);

	/*
	 *  Let VOGLE/VOGL know about the window size.
	 */
        vdevice.sizeX = vdevice.sizeY = MIN(w, h);
	vdevice.sizeSx = w;
	vdevice.sizeSy = h;

	/*
	 * Set up batching.....(for speed and "pseudo" double buffering)
	 */
	pw_batch_on(pw);

	return(1);
}

/*
 * SUN_exit
 *
 *	cleans up before returning the window to normal.
 */
SUN_exit()
{
	long	nbytes;
	int	i;
	Event	event;

	if (use_sunview_canvas)
		return(1);

	/*
	 * Flush all events for this window.
	 *
	 * While doing non-blocking input input_readevent returns -1 and
	 * errno == EWOULDBLOCK when everything has been read, so if 
	 * errno != EWOULDBLOCK then something naughty happened...
	 */
	while (input_readevent(wfd, &event) >= 0)
		;

	if (errno != EWOULDBLOCK) {
		perror("SUN_exit(flushing), input_readevent");
		exit();
	}

	/*
	 * reset wfd to blocking input.
	 */
	if (fcntl(wfd, F_SETFL, oldflags) < 0) {
		perror("oldflags, F_SETFL");
		exit(1);
	}

	if (blanket_win)
		win_removeblanket(wfd);
	else 
		win_remove(wfd);

	signal(SIGWINCH, SIG_DFL);

	return(1);
}

/*
 * SUN_draw
 *
 *	draws a line from the current graphics position to (x, y).
 *
 * Note: (0, 0) is defined as the top left of the window on a sun (easy
 * to forget).
 */
SUN_draw(x, y)
	int	x, y;
{
	if (!lineflags)	/* If thin and solid */
		pw_vector(pw, vdevice.cpVx, vdevice.sizeSy - vdevice.cpVy, x, vdevice.sizeSy - y, PIX_SRC | PIX_COLOR(colour), colour);
	else	/* Fat and/or dashed */
		pw_line(pw, vdevice.cpVx, vdevice.sizeSy - vdevice.cpVy, x, vdevice.sizeSy - y, &brush, t, PIX_SRC | PIX_COLOR(colour));

	if (vdevice.sync)
		pw_show(pw);
}

/*
 * SUN_getkey
 *
 *	grab a character from the keyboard.
 */
int
SUN_getkey()
{
	Event	event;

	pw_show(pw);

	do {
		while ((input_readevent(wfd, &event) < 0) && (errno == EWOULDBLOCK))
		;	/* Nothing to read - wait for something */
		
	} while (!event_is_ascii(&event));	/* Wait for a key press */

	return(event_action(&event));
}

/*
 * SUN_checkkey
 *
 *	Check if a keyboard key has been hit. If so return it.
 */
static Event *saved_event;

int
SUN_checkkey()
{
	static Event	event;
	
	if(saved_event == &event)
	    saved_event = NULL;

	if (saved_event && event_is_ascii(saved_event)) {
		Event *tmp = saved_event;

		saved_event = NULL;
		return(event_action(tmp));
	}

	saved_event = NULL;

	if (input_readevent(wfd, &event) < 0) {
		if (errno == EWOULDBLOCK) {
			return(0);
		} else {
			perror("SUN_checkkey, input_readevent");
			exit(1);
		}
	} else if (event_is_ascii(&event))
		return(event_action(&event));
	else
		saved_event = &event;

	return(0);
}
		

/*
 * SUN_locator
 *
 *	return the window location of the cursor, plus which mouse button,
 * if any, is been pressed.
 */
int
SUN_locator(wx, wy)
	int	*wx, *wy;
{
	int	but;
	static Event	event;
	
	if (vdevice.sync)
		pw_show(pw);

	but = 0;

	*wx = win_get_vuid_value(wfd, LOC_X_ABSOLUTE);
	*wy = (int)vdevice.sizeSy - win_get_vuid_value(wfd, LOC_Y_ABSOLUTE);

	/* This used to work under 4.0 but not longer ..... Maybe
	 * SUN don't want us to use sunview anymore?
	 */
	if (win_get_vuid_value(wfd, BUT(1)))
		but |= 1;

	if (win_get_vuid_value(wfd, BUT(2)))
		but |= 2;

	if (win_get_vuid_value(wfd, BUT(3)))
		but |= 3;

	
	if (use_sunview_canvas)
		return(but);

	if (saved_event == &event) /*Don't Re-use a previously rejected event*/
	    saved_event = NULL;

	if (saved_event && event_is_button(saved_event) && event_is_down(saved_event)) {
		if (event_action(saved_event) == MS_LEFT)
			but |= 1;
		if (event_action(saved_event) == MS_MIDDLE)
			but |= 2;
		if (event_action(saved_event) == MS_RIGHT)
			but |= 4;

		saved_event = NULL;
		return(but);
	}


	if (input_readevent(wfd, &event) < 0) {
		if (errno == EWOULDBLOCK) {
			return(0);
		} else {
			perror("SUN_locator, input_readevent");
			exit(1);
		}
	} else if (event_is_button(&event) && event_is_down(&event)) {
			if (event_action(&event) == MS_LEFT)
				but |= 1;
			if (event_action(&event) == MS_MIDDLE)
				but |= 2;
			if (event_action(&event) == MS_RIGHT)
				but |= 4;
	} else
		saved_event = &event;

	return(but);
}

#ifdef VOGLE
/*
 * SUN_clear
 *
 *	Clear the screen to current colour
 */
SUN_clear()
{
	pw_writebackground(pw, 0, 0, w, h, PIX_SRC | PIX_COLOR(colour) | PIX_DONTCLIP);

	if (vdevice.sync)
		pw_show(pw);
}

#else

/*
 * SUN_clear
 *
 *	Clear the viewport to current colour
 */
SUN_clear()
{
        unsigned int    vw = vdevice.maxVx - vdevice.minVx;
        unsigned int    vh = vdevice.maxVy - vdevice.minVy;

	pw_writebackground(pw, 
		vdevice.minVx, vdevice.sizeSy - vdevice.maxVy, 
		vw, vh,
		PIX_SRC | PIX_COLOR(colour) | PIX_DONTCLIP
	);

	if (vdevice.sync)
		pw_show(pw);
}

#endif

/*
 * SUN_color
 *
 *	set the current drawing color index.
 */
SUN_color(ind)
        int	ind;
{
	colour = ind;
}

/*
 * SUN_mapcolor
 *
 *	change index i in the color map to the appropriate r, g, b, value.
 */
SUN_mapcolor(i, r, g, b)
	int	i;
	int	r, g, b;
{
	int	j;

	if (i >= 255 || vdevice.depth == 1)
		return(-1);

	if (i >= colnum)
		colnum = i;

	red[i] = r;
	green[i] = g;
	blue[i] = b;

	red[255] = (u_char)~red[0];
	green[255] = (u_char)~green[0];
	blue[255] = (u_char)~blue[0];

	pw_putcolormap(pw, 0, colnum, red, green, blue);
#ifdef BASE_COL_SET
	pw_putcolormap(pw0, 0, colnum, red, green, blue);
#endif
}
	
/*
 * SUN_font
 *
 *   Set up a hardware font. Return 1 on success 0 otherwise.
 *
 */
SUN_font(fontfile)
        char	*fontfile;
{
	char	name[BUFSIZ];

	if (font_id != (Pixfont *)NULL)
		pf_close(font_id);

	if ((font_id = pf_open(fontfile)) == NULL)
		if (*fontfile != '/') {
			strcpy(name, STDFONTDIR);
			strcat(name, fontfile);
			if ((font_id = pf_open(name)) == NULL)
				return(0);
		} else 
			return(0);

	vdevice.hheight = font_id->pf_defaultsize.y;
	vdevice.hwidth = font_id->pf_defaultsize.x;

	return(1);
}

/* 
 * SUN_char
 *
 *	 outputs one char - is more complicated for other devices
 */
SUN_char(c)
	char	c;
{
	char	*s = " ";

	s[0] = c;
	pw_ttext(pw, vdevice.cpVx, (int)(vdevice.sizeSy - vdevice.cpVy), PIX_SRC | PIX_COLOR(colour), font_id, s);

	if (vdevice.sync)
		pw_show(pw);
}

/*
 * SUN_string
 *
 *	Display a string at the current drawing position.
 */
SUN_string(s)
        char	s[];
{
	pw_ttext(pw, vdevice.cpVx, (int)(vdevice.sizeSy - vdevice.cpVy), PIX_SRC | PIX_COLOR(colour), font_id, s);

	if (vdevice.sync)
		pw_show(pw);
}

/*
 * SUN_fill
 *
 *	fill a polygon
 */
SUN_fill(n, x, y)
	int	n, x[], y[];
{
	struct	pr_pos	vlist[128];
	int	i, npnts;

	if (n > 128)
		verror("vogle: more than 128 points in a polygon");

	npnts = n;

	for (i = 0; i < n; i++) {
		vlist[i].x = x[i];
		vlist[i].y = vdevice.sizeSy - y[i];
	}

	pw_polygon_2(pw, 0, 0, 1, &npnts, vlist, PIX_SRC | PIX_COLOR(colour) | PIX_DONTCLIP, (Pixwin *)NULL, 0, 0);

	vdevice.cpVx = x[n-1];
	vdevice.cpVy = y[n-1];

/*
	if (vdevice.sync)
		pw_show(pw);
*/
}

/*
 * SUN_backb
 *
 *	swap to memory only drawing (backbuffer) - a little slow but it
 * works on everything. 
 */
SUN_backb()
{
	/*
	 * Batching is already on....
	 */
	return(0);
}

/*
 * SUN_swapb
 *
 *	swap the front and back buffers.
 */
SUN_swapb()
{
	if (vdevice.inbackbuffer)
		pw_show(pw);

	return(0);
}

/*
 * SUN_frontb
 *
 *	draw in the front buffer
 */
SUN_frontb()
{
	/*
	 * Make it visible anyway....
	 */
	pw_show(pw);
}

/*
 * SUN_sync
 *
 *	Syncronise the display with what we thing has been output
 */
SUN_sync()
{
	pw_show(pw);
}

/*
 * SUN_setls
 *
 *	Set the line style....
 */
SUN_setls(lss)
	int	lss;
{
	
	unsigned ls = lss;
	static short	dashes[17];

	int	i, n, a, b, offset;

	if (ls == 0xffff) {
		lineflags &= ~DASHED;
		t = (Pr_texture *)NULL;
		return;
	}

	lineflags |= DASHED;

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
	dashes[n] = 0;

	tex.pattern = &dashes[0];
	tex.offset = offset;
	t = &tex;
}

/*
 * SUN_setlw
 *
 *	Set the line style....
 */
SUN_setlw(lw)
	int	lw;
{
	if (lw > 1) {
		lineflags |= FATLINES;
		brush.width = lw;
	} else
		lineflags &= ~FATLINES;
}

/*
 * the device entry
 */
static DevEntry sundev = {
	"sun",
	"screen.b.16",
	"screen.b.12",
	SUN_backb,
	SUN_char,
	SUN_checkkey,
	SUN_clear,
	SUN_color,
	SUN_draw,
	SUN_exit,
	SUN_fill,
	SUN_font,
	SUN_frontb,
	SUN_getkey,
	SUN_init,
	SUN_locator,
	SUN_mapcolor,
#ifndef VOGLE
	SUN_setls,
	SUN_setlw,
#endif
	SUN_string,
	SUN_swapb,
	SUN_sync
};

/*
 * _SUN_devcpy
 *
 *	copy the sun device into vdevice.dev.
 */
_SUN_devcpy()
{
	vdevice.dev = sundev;
}

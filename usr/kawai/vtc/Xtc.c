/*
 * 1992 Author: Takanori Nagae
 *
 *
 * xtc : Turbo C like graphics library on X
 * compile: use makefile
 * PseudoColor only.
 */

#include <stdio.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "Xtc.h"

#define FONT_NAME	 "-adobe-courier-bold-r-normal-*-24-*-*-*-*-*-*-*"

static int op_table[5] ={
	GXcopy, GXxor, GXor, GXand, GXcopyInverted
};

static unsigned char color_table[256];
static unsigned long pixels[NCOLORS];
static int fgcolor, bkcolor, CPx, CPy;

static Display *disp;
static int screen, depth;
static Window win, canvas;
static Colormap cmap;
static GC gc, gc_putpixel, gc_putimage, clearGC, pmGC;
static XFontStruct *font_info;
static Pixmap pixmap;
static int window_size_x=SCREEN_WIDTH,window_size_y=SCREEN_HEIGHT;
static int window_size_pixmap_x=SCREEN_WIDTH,window_size_pixmap_y=SCREEN_HEIGHT;

static void PolypointsToXPoints(npoints, polypoints, xpoints)
int npoints, *polypoints;
XPoint *xpoints;
{
	for(; 0 < npoints; npoints--){
		xpoints->x = *polypoints++;
		xpoints->y = *polypoints++;
		xpoints++;
	}
}

void xtcmainloop(button)
int button;
{
	XEvent ev;
	char c;

	XFlush(disp);
	while(1){
		XNextEvent(disp, &ev);
		cur_pointer_x = ev.xbutton.x;
		cur_pointer_y = ev.xbutton.y;
		switch(ev.type){
			case ButtonPress:
					if(button == ev.xbutton.button)
					return;
				break;
			case KeyPress:
				if(1 == XLookupString(&ev, &c, 1, NULL, NULL)){
					switch(c){
						case 'q':
						case 'Q':
						case '\03':
							return;
					}
				}
				break;
		}
	}
}

void xflush()
{
	XFlush(disp);
}

int xgetbutton()
{
	XEvent ev;

	while(1){
		XNextEvent(disp, &ev);

		cur_pointer_x = ev.xbutton.x;
		cur_pointer_y = ev.xbutton.y;
		if(ButtonPress == ev.type)
			return ev.xbutton.button;
	}
}

int xgetbutton_now()
{
	XEvent ev;

	XNextEvent(disp, &ev);

	cur_pointer_x = ev.xbutton.x;
	cur_pointer_y = ev.xbutton.y;
	if(ButtonPress == ev.type){
		return ev.xbutton.button;
	}

        return 0;
}

void arc(x, y, stangle, endangle, radius)
int x, y, stangle, endangle, radius;
{
	XDrawArc(disp, win, gc, x - radius, y - radius,
		2 * radius, 2 * radius, 64 * stangle, 64 *(endangle - stangle));
}

void bar(left, top, right, bottom)
int left, top, right, bottom;
{
	XFillRectangle(disp, win, gc, left, top,
		right - left + 1, bottom - top + 1);
}

void circle(x, y, radius)
int x, y, radius;
{
	XDrawArc(disp, win, gc, x - radius, y - radius,
		2 * radius, 2 * radius, 0, 23040);
}

void circle_pixmap(x, y, radius)
int x, y, radius;
{
	XDrawArc(disp, pixmap, pmGC, x - radius, y - radius,
		2 * radius, 2 * radius, 0, 23040);
}

void cleardevice()
{
	XClearWindow(disp, win);
	CPx = CPy = 0;
}

void cleardevice_pixmap()
{
        XFillRectangle(disp,pixmap,clearGC,0,0,window_size_pixmap_x,window_size_pixmap_y);
/*	XClearWindow(disp, pixmap);*/
	CPx = CPy = 0;
}

void closegraph()
{
	XUnmapWindow(disp, pixmap);
	XUnmapWindow(disp, win);
}

void drawpoly(numpoints, polypoints)
int numpoints, *polypoints;
{
	XPoint *xpoints =(XPoint *)malloc(numpoints * sizeof(XPoint));

	PolypointsToXPoints(numpoints, polypoints, xpoints);
	XDrawLines(disp, win, gc, xpoints, numpoints, CoordModeOrigin);
	free(xpoints);
}

void ellipse(x, y, stangle, endangle, xradius, yradius)
int x, y, stangle, endangle, xradius, yradius;
{
	XDrawArc(disp, win, gc, x - xradius, y - yradius,
		2 * xradius, 2 * yradius, 64 * stangle, 64 *(endangle - stangle)); 
}

void fillellipse(x, y, xradius, yradius)
int x, y, xradius, yradius;
{
	XFillArc(disp, win, gc, x - xradius, y - yradius,
		2 * xradius, 2 * yradius, 0, 23040); 
}

void fillellipse_pixmap(x, y, xradius, yradius)
int x, y, xradius, yradius;
{
	XFillArc(disp, pixmap, pmGC, x - xradius, y - yradius,
		2 * xradius, 2 * yradius, 0, 23040); 
}

void fillpoly(numpoints, polypoints)
int numpoints, *polypoints;
{
	XPoint *xpoints =(XPoint *)malloc(numpoints * sizeof(XPoint));

	PolypointsToXPoints(numpoints, polypoints, xpoints);
	XFillPolygon(disp, win, gc, xpoints, numpoints, Complex,
		CoordModeOrigin);
	free(xpoints); 
}

int getbkcolor(){ return bkcolor; }

int getcolor(){ return fgcolor; }

int getmaxcolor() { return NCOLORS - 1; }

int getmaxx() { return window_size_x - 1; }

int getmaxy() { return window_size_y - 1; }

int getmaxx_pixmap() { return window_size_pixmap_x - 1; }

int getmaxy_pixmap() { return window_size_pixmap_y - 1; }

void getimage(left, top, right, bottom, bitmap)
int left, top, right, bottom;
void *bitmap;
{
	XImage *image;
	int width, height, i, j, k; 
	unsigned char *data, *buffer =(unsigned char *)bitmap;

	width = right - left + 1;
	height = bottom - top + 1;
	image = XGetImage(disp, win, left, top, width, height,
				AllPlanes, ZPixmap);
	data =(unsigned char *)image->data;
	*buffer++ = width;
	*buffer++ = height;
	for(i = 0; height > i; i++){
		for(j = k = 0; image->bytes_per_line > j; j++, data++){
			if(j < image->xoffset) continue;
			k++;
			if(width < k) continue;
			*buffer++ = color_table[*data];
		}
	}
}

int getpixel(x, y)
int x, y;
{
	XImage *image;
	int val;

	image = XGetImage(disp, win, x, y, 1, 1, AllPlanes, ZPixmap);
	val = color_table[image->data[image->xoffset]];
	XDestroyImage(image);
	return(unsigned)val;
}

int getx() { return CPx; }

int gety() { return CPy; }

int imagesize(left, top, right, bottom)
int left, top, right, bottom;
{
	return(unsigned)(2 +(right - left + 1)*(bottom - top + 1));
}

void set_window_size(x,y)
int x,y;
{
    window_size_x=x;
    window_size_y=y;
}

void set_window_size_pixmap(x,y)
int x,y;
{
    window_size_pixmap_x=x;
    window_size_pixmap_y=y;
}

void initgraph(graphdriver, graphmode, pathtodriver)
int *graphdriver, *graphmode;
char *pathtodriver;
{
	static unsigned char def_rgb[NCOLORS][3] = {
		{0, 0, 0},		/* black */
		{0, 0, 7},		/* blue */
		{0, 7, 0},		/* green */
		{0, 7, 7},		/* cyan */
		{7, 0, 0},		/* red */
		{0, 7, 7},		/* magenta */
		{7, 7, 0},		/* brown */
		{11, 11, 11},	/* light gray */
		{5, 5, 5},		/* gray */
		{0, 0, 15},		/* light blue */
		{0, 15, 0},		/* light green */
		{0, 15, 15},	/* light cyan */
		{15, 0, 0},		/* right red */
		{15, 0, 15},	/* right magenta */ 
		{15, 15, 0},	/* yellow */
		{15, 15, 15},	/* white */
	};

	XEvent	ev;
	int i;
	unsigned long plane_mask[1];

	disp = XOpenDisplay(NULL);
	screen = DefaultScreen(disp);
	win = XCreateSimpleWindow(disp, RootWindow(disp, screen) , 0, 0,
			window_size_x, window_size_y, 1, 0, BlackPixel(disp, screen));
	depth = DefaultDepth(disp, screen);
	cmap = DefaultColormap(disp, screen);
	XSelectInput(disp, win, ExposureMask | ButtonPressMask | KeyPressMask);
	XStoreName(disp, win, "Xtc");

	XAllocColorCells(disp, cmap, False, plane_mask, 0, pixels, NCOLORS);
	for(i = 0; NCOLORS > i; i++) color_table[pixels[i]] = i;

	font_info = XLoadQueryFont(disp, FONT_NAME);

	gc = XCreateGC(disp, win, 0, NULL);
	XSetForeground(disp, gc, pixels[WHITE]); 
	XSetFont(disp, gc, font_info->fid);

	gc_putpixel = XCreateGC(disp, win, 0, NULL);
	gc_putimage = XCreateGC(disp, win, 0, NULL);

	for(i = 0; NCOLORS > i; i++){
		setrgbpalette(i, def_rgb[i][0], def_rgb[i][1], def_rgb[i][2]);
	}
	bkcolor = BLACK;
	fgcolor = WHITE;

	XMapWindow(disp, win);

	XFlush(disp);

	do{ XNextEvent(disp, &ev); } while(Expose != ev.type);

	canvas=XCreateSimpleWindow(disp,win,0,0,window_size_pixmap_x,window_size_pixmap_y,1,
				   BlackPixel(disp,DefaultScreen(disp)),
				   WhitePixel(disp,DefaultScreen(disp)));
	pixmap=XCreatePixmap(disp,canvas,window_size_pixmap_x,window_size_pixmap_y,
			     DefaultDepth(disp,DefaultScreen(disp)));
/*	XMapWindow(disp, pixmap);---------------------- */
	clearGC=XCreateGC(disp,pixmap,0,0);
	XSetForeground(disp,clearGC,BlackPixel(disp,DefaultScreen(disp)));
	XFillRectangle(disp,pixmap,clearGC,0,0,window_size_pixmap_x,window_size_pixmap_y);
	pmGC=XCreateGC(disp,pixmap,0,0);
	XSetForeground(disp,pmGC,WhitePixel(disp,DefaultScreen(disp)));
	XSetBackground(disp,pmGC,BlackPixel(disp,DefaultScreen(disp)));
	
	XSetFont(disp, pmGC, font_info->fid);

/*	XFillRectangle(disp,pixmap,pmGC,0,0,100,100);*/
	

}

void copy_from_pixmap(des_x,des_y)
int des_x,des_y;
{
/*        XSetForeground(disp,pmGC,pixels[12]);
	XFillRectangle(disp,pixmap,pmGC,0,0,50,50);*/
	XCopyArea(disp,pixmap,win,pmGC,0,0,window_size_pixmap_x,window_size_pixmap_y,des_x,des_y);
}


void line(x1, y1, x2, y2)
int x1, y1, x2, y2;
{
	XDrawLine(disp, win, gc, x1, y1, x2, y2);
}

void linerel(dx, dy)
int dx, dy;
{
	XDrawLine(disp, win, gc, CPx, CPy, CPx + dx, CPy + dy);
	CPx += dx; CPy += dy;
}

void lineto(x, y)
int x, y;
{
	XDrawLine(disp, win, gc, CPx, CPy, x, y);
	CPx = x; CPy = y;
}

void moverel(dx, dy)
int dx, dy;
{
	CPx += dx; CPy += dy;
}

void moveto(x, y)
int x, y;
{
	CPx = x; CPy = y;
}

void outtextxy(x, y, textstring)
int x, y;
char *textstring;
{
	XDrawString(disp, win, gc, x, font_info->max_bounds.ascent + y,
		textstring, strlen(textstring));
}

void outtextxy_pixmap(x, y, textstring)
int x, y;
char *textstring;
{
	XDrawString(disp, pixmap, pmGC, x, font_info->max_bounds.ascent + y,
		textstring, strlen(textstring));
}

void pieslice(x, y, stangle, endangle, radius)
int x, y, stangle, endangle, radius;
{
	XFillArc(disp, win, gc, x - radius, y - radius,
		2 * radius, 2 * radius, 64 * stangle, 64 *(endangle - stangle));
}

void putimage(left, top, bitmap, op)
int left, top, op;
void *bitmap;
{
	static int sta_op = COPY_PUT;
	int width, height;
	XImage *image;
	unsigned char *b, *buffer =(unsigned char *)bitmap;
	int	i;

	width = *buffer++;
	height = *buffer++;
	for(b = buffer, i = 0; width * height > i; i++) *b++ = pixels[*b];
	if(op != sta_op){
		sta_op = op;
		XSetFunction(disp, gc_putimage, op_table[sta_op]);
	}
	image = XCreateImage(disp, DefaultVisual(disp, screen), depth,
				ZPixmap, 0,(char *)buffer, width, height, 8, 0);
	XPutImage(disp, win, gc_putimage, image, 0, 0, left, top, 
		width, height);
	XDestroyImage(image);
}

void putpixel(x, y, color)
int x, y, color;
{
	static int sta_color = WHITE;

	if(sta_color != color){
		sta_color = color;
		XSetForeground(disp, gc_putpixel, pixels[sta_color]);
	}
	XDrawPoint(disp, win, gc_putpixel, x, y);
}

void rectangle(left, top, right, bottom)
int left, top, right, bottom;
{
	XDrawRectangle(disp, win, gc, left, top, right - left, bottom - top);
}

void sector(x, y, stangle, endangle, xradius, yradius)
int x, y, stangle, endangle, xradius, yradius;
{
	XFillArc(disp, win, gc, x - xradius, y - yradius,
		2 * xradius, 2 * yradius, 64 * stangle, 64 *(endangle - stangle));
}

void setbkcolor(color)
int color;
{
	XSetWindowBackground(disp, win, pixels[bkcolor = color]);
	XClearWindow(disp, win);
}

void setcolor(color)
int color;
{
	XSetForeground(disp, gc, pixels[fgcolor = color]);
}

void setcolor_pixmap(color)
int color;
{
	XSetForeground(disp, pmGC, pixels[fgcolor = color]);
}

void setlinestyle(linestyle, upattern, thickness)
int linestyle, thickness;
unsigned upattern;
{
	static char dotted_line[] = {1, 3};
	static char center_line[] = {1, 3, 4, 3};
	static char dashed_line[] = {7, 4};

	switch(linestyle){
		case DOTTED_LINE:
			XSetDashes(disp, gc, 0, dotted_line, 2);
			break;
		case CENTER_LINE:
			XSetDashes(disp, gc, 0, center_line, 4);
			break;
		case DASHED_LINE:
			XSetDashes(disp, gc, 0, dashed_line, 2);
	}	
	XSetLineAttributes(disp, gc, thickness,
		((0 == linestyle)? LineSolid: LineOnOffDash),
		CapButt, JoinBevel);
} 

void setrgbpalette(colornum, red, green, blue)
int colornum, red, green, blue;
{
	XColor color;

	color.pixel = pixels[colornum];	
	color.red = 4096 * red;
	color.green = 4096 * green;
	color.blue = 4096 * blue;
	color.flags = DoRed | DoGreen | DoBlue;
#if 1 /* AK */
	XAllocColor(disp, cmap, &color);
	pixels[colornum] = color.pixel;
	/* look for closest color from current colormap,
	   insted of sorting a new color */
#else
	XStoreColor(disp, cmap, &color); 
#endif
}

void setwritemode(mode)
int mode;
{
	XSetFunction(disp, gc, op_table[mode]);
}

int textheight(textstring)
char *textstring;
{
	return font_info->max_bounds.ascent + font_info->max_bounds.descent;
}

int textwidth(textstring)
char *textstring;
{
	return XTextWidth(font_info, textstring, strlen(textstring));
}

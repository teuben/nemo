/*
 *   Xutil.c - Part of my_x and yapp_x (NEMO) graphics driver programs
 *         Provides generic X11-based plotting functions
 *         
 *   Author:     Patrick (AdM) Frisch ( frisch@walhall.uni-sw.gwdg.de )
 *   Date:       sometimes in 1992
 *   Revised:    March 93
 *   Version:    1.01
 *   Module
 *   Remarks:    Definition of user-callable functions, easy interface
 *               between complicated X stuff and user.
 *
 *
 *   Global Functions:
 *
 *               x_space(real,real,real,real)
 *               x_erase(void)
 *               x_point(real,real)
 *               x_line(real,real,real,real)
 *               x_move(real, real)
 *               x_drawto(real,real)
 *               x_circle(real,real,real)
 *               x_arc(real,real,real,real,real,real)
 *               x_label(char *)
 *               x_linemod(char *, int)
 *               x_color_by_number(real,real,real)
 *               x_color_by_name(char *)
 *               x_number_of_colors()
 *               x_calc_size_in_pixel()
 *               x_screendump()
 *               x_doplot()
 *               x_buttonwait()   (*this one really bad here*)
 *
 */

#include "config.h"
#include <X11/Xlib.h>

extern Display *x_dpy;
extern Window x_win;
extern Pixmap x_pixmap;
extern GC x_gc;
extern XWindowAttributes x_wa;

extern int rwidth, rheight;
extern int fgcolor, bgcolor, initfg;

extern int x_filling;

real x_xmin, x_xmax, x_ymin, x_ymax;
real x_xfactor, x_yfactor;

#define MAPX(x) (int) (((x) - x_xmin) * x_xfactor)
#define MAPY(y) (int) ((x_ymax - (y)) * x_yfactor)

local int spaceset = 0;    

real x_curpt[] = {0,0};

int (*event_function)() = NULL;


typedef struct {
    char type;
    float x,y;
    float tx,ty;
    void *opt;
} plobjrec, *plobjrecptr;

static plobjrecptr plobj=NULL;
static int nplobj = 0;
static int do_replot = 0;
extern int x_no_resize;

local Window_Dump(FILE *);
local Image_Size(XImage *);
local Get_XColors(XWindowAttributes *, XColor **);
local _swapshort(char *, unsigned);
local _swaplong(char *, unsigned);

x_init_plobj()
{
    plobj = (plobjrecptr)malloc(sizeof(plobj)*1000);
    if(plobj == NULL)
    {
	fprintf(stderr,"Can't allocate objects ...\n");
	exit(1);
    }
    x_no_resize = 0;
    x_set_color(initfg);
    x_linemod("solid",1);
}

static reset_plobj()
{
    int i;

    for(i = 0; i < nplobj; i++)
	if(plobj[i].opt)
	    free(plobj[i].opt);

    nplobj = 0;
}

x_close_plobj()
{
    reset_plobj();
    free(plobj);
    x_no_resize = 1;
}

x_insert_object(type, xp, yp, xt, yt, opt)
char type;
real xp, yp, xt, yt;
char *opt;
{
    if(do_replot || x_no_resize) /* don't insert during replot */
	return;
    
    plobj[nplobj].type = type;
    plobj[nplobj].x = (float)xp;
    plobj[nplobj].y = (float)yp;
    plobj[nplobj].tx = (float)xt;
    plobj[nplobj].ty = (float)yt;
    if(opt)
    {
	plobj[nplobj].opt = (void *)malloc(strlen(opt)+1);
	if(plobj[nplobj].opt == NULL)
	{
	    fprintf(stderr,"Can't allocate string: %s ...\n", opt);
	    exit(1);
	}
	strcpy((char *)plobj[nplobj].opt, opt);
    }
    else
	plobj[nplobj].opt = (void *) NULL;
    
    if(++nplobj >= 1000)
    {
	fprintf(stderr,"Too many plot objects ... sorry\n");
	exit(1);
    }
	
}

x_replot()
{
    int i;

    if(x_no_resize)
	return;
    
    do_replot = 1;
    x_clear_window();
    for(i = 0; i < nplobj; i++)
    {
	switch(plobj[i].type)
	{
	  case 'p':
	    x_point((real)plobj[i].x, (real)plobj[i].y);
	    break;
	  case 'd':
	    x_move((real)plobj[i].x, (real)plobj[i].y);
	    break;
	  case 'l':
	    x_line((real)plobj[i].x, (real)plobj[i].y,
		   (real)plobj[i].tx, (real)plobj[i].ty);
	    break;
	  case 'm':
	    x_linemod(plobj[i].opt, (int)plobj[i].x);
	    break;
	  case 'c':
	    x_set_color((int)plobj[i].x);
	  case 't':
	    x_label(plobj[i].opt);
	    break;
	  case 's':
	    x_alabel((char)plobj[i].x, (char)plobj[i].y,
		     plobj[i].opt, (real)plobj[i].tx);
	    break;
	  case 'f':
	    if(plobj[i].opt)
		x_font_by_name(plobj[i].opt);
	    else
		x_fontsize((real)plobj[i].x);
	    break;
	  case 'a':
	    x_filling = (int)plobj[i].ty;
	    x_circle((real)plobj[i].x, (real)plobj[i].y,
		     (real)plobj[i].tx);
	    x_filling = 0;
	    break;
	}
    }
    do_replot = 0;
    x_doplot();
}

x_space(x0,y0,x1,y1)
real x0,y0,x1,y1;
{
    spaceset = 1;
    x_xmin = x0;
    x_ymin = y0;
    x_xmax = x1;
    x_ymax = y1;
    x_xfactor = rwidth/(x_xmax - x_xmin);
    x_yfactor = rheight/(x_ymax-x_ymin);
}

static x_drawellipse(x,y,r1,r2)
int x,y,r1,r2;
{
#if 1
    /* nice feature, keep it for later use (AdM) */   
    if (x_filling) {
	XFillArc(x_dpy, x_pixmap, x_gc, x - r1, y - r2, 2*r1, 2*r2,
		 0, 64*180);
	XFillArc(x_dpy, x_pixmap, x_gc, x - r1, y - r2, 2*r1, 2*r2,
		 64*180, 64*180);
    }
#endif    
    XDrawArc(x_dpy, x_pixmap, x_gc, x - r1, y - r2, 2*r1, 2*r2,
	     0, 64*180);
    XDrawArc(x_dpy, x_pixmap, x_gc, x - r1, y - r2, 2*r1, 2*r2,
	     64*180, 64*180);
}

static x_drawarc(x,y,x0,y0,x1,y1)
int x,y,x0,y0,x1,y1;
{
    int a0,b0,a1,b1;
    int a02,b02,a12,b12;
    real ar,br,ar2,br2;
    real theta0,theta1;

    a0 = x0 - x;
    a02 = a0*a0;
    a1 = x1 - x;
    a12 = a1*a1;
    b0 = y0 - y;
    b02 = b0*b0;
    b1 = y1 - y;
    b12 = b1*b1;
    if (b12 == b02) return 0;
    ar2 = (a02*b12 - a12*b02)/(b12 - b02);
    if (ar2 < 0) return 0;
    ar = sqrt(ar2);
    br2 = (b02*a12 - b12*a02)/(a12 - a02);
    if (br2 < 0) return 0;
    br = sqrt(br2);
    theta0 = -atan2(b0/br,a0/ar);
    theta1 = -atan2(b1/br,a1/ar);
    if (theta0 > theta1) theta1 += 2 * M_PI;
    if (theta0 < 0){
	theta0 += 2 * M_PI;
	theta1 += 2 * M_PI;
    }
    XDrawArc(x_dpy,x_pixmap,x_gc,(int) (x - ar), (int) (y - br),
	     (int)(2 * ar), (int) (2 * br), (int) (64* theta0
	     * 180/M_PI), (int) (64 * (theta1 - theta0) * 180/M_PI));
}

x_clear_window()
{
    XSetForeground(x_dpy,x_gc,bgcolor);
    XFillRectangle(x_dpy,x_pixmap,x_gc,0,0,rwidth,rheight);
    XSetForeground(x_dpy,x_gc,fgcolor);
}

x_erase()
{
    x_clear_window();
    reset_plobj();
}

x_point(x,y)
real x,y;
{
    x_curpt[0] = x;
    x_curpt[1] = y;
    if (spaceset == 0) return;
    XDrawPoint(x_dpy,x_pixmap,x_gc, MAPX(x), MAPY(y));
    x_insert_object('p',x,y,0.,0.,NULL);
}

x_line(x1,y1,x2,y2)
real x1,y1,x2,y2;
{
    x_curpt[0] = x2;
    x_curpt[1] = y2;
    if (spaceset == 0) return;
    XDrawLine(x_dpy,x_pixmap,x_gc, MAPX(x1), MAPY(y1), MAPX(x2), MAPY(y2));
    x_insert_object('l',x1,y1,x2,y2,NULL);
}

x_move(x,y)
real x,y;
{
    x_curpt[0] = x;
    x_curpt[1] = y;
    x_insert_object('d',x,y,0.,0.,NULL);
}

x_drawto(x, y)
real x,y;
{
    if (spaceset == 0) return;
    XDrawLine(x_dpy,x_pixmap,x_gc, MAPX(x_curpt[0]), MAPY(x_curpt[1]),
	      MAPX(x), MAPY(y));
    x_insert_object('l',x_curpt[0],x_curpt[1],x,y,NULL);
    x_curpt[0] = x;
    x_curpt[1] = y;
}

x_circle(x,y,r)
real x,y,r;
{
    if (spaceset == 0) return;
    x_drawellipse(MAPX(x), MAPY(y), (int)(r * x_xfactor), (int)(r*x_yfactor));
    x_insert_object('a', x, y, r, (real)x_filling, NULL);
}

x_arc(x,y,x0,y0,x1,y1)
real x,y,x0,y0,x1,y1;
{
    if (spaceset == 0) return;
    x_drawarc(MAPX(x), MAPY(y), MAPX(x0), MAPY(y0), MAPX(x1), MAPY(y1));
}

x_label(s)
char *s;
{
    XDrawString(x_dpy,x_pixmap,x_gc,MAPX(x_curpt[0]),MAPY(x_curpt[1]),s,
		strlen(s));
    x_insert_object('t',0.,0.,0.,0.,s);
}

/* Now to control line styles */

static char *linemode_labels[] ={
    "solid",
    "dotted",
    "shortdashed",
    "dotdashed",
    "longdashed",
    "disconnected",
};

static int no_of_linemodes = 6;

static char dashes[][5] = {
    {0},                    /* solid */
    {1,1,0},                /* dotted */
    {2,2,0},                /* shortdashed */
    {4,2,1,2,0},            /* dotdashed */
    {4,4,0},                /* longdashed */
    {1,8,0}                 /* disconnected */
};

x_linemod(s, lwid)
char *s;
int lwid;
{
    int i=0;
    /* search through the linemode names for a match.  */
    while ((i < no_of_linemodes) && (strcmp(s, linemode_labels[i]) != 0))
    i++;

    if (dashes[i][0])
        XSetDashes(x_dpy, x_gc, 0, dashes[i], strlen (dashes[i]));
    XSetLineAttributes (x_dpy, x_gc, lwid /* width */,
			(i == 0) ? LineSolid : LineOnOffDash,
			CapButt, JoinBevel);
    x_insert_object('m',(real)lwid,0.,0.,0.,s);
}

x_set_color(color)
int color;
{
    XSetForeground(x_dpy, x_gc, color);
    x_insert_object('c',(real)color,0.,0.,0.,NULL);
}

/* Color range is 0...65536 */

x_color_by_number(red, green, blue)
unsigned short red, green, blue;
{
    XColor col;
    col.red = red;
    col.green = green;
    col.blue = blue;
    col.flags = DoRed | DoGreen | DoBlue;
    if (XAllocColor (x_dpy, x_wa.colormap, &col)) {
	fgcolor = col.pixel;
	x_set_color(fgcolor);
    } else {
	fprintf(stderr, "Warning: Cannot allocate color %d %d %d\
                        (RGB).\n", red, green, blue);
    }
}

x_color_by_name(colname)
char *colname;
{
    if (nameToPixel(colname, &fgcolor)) {
	x_set_color(fgcolor);
	return 1;
    }
    return 0;
}

x_doplot()
{
    XClearArea(x_dpy,x_win,0,0,0,0,True); /* Show the image! */
    XFlush(x_dpy);
}

int default_event_function(ev)
XEvent *ev;
{
    if(ev->type == ButtonPress)
        return(1);
    if (ev->type == ConfigureNotify) /* Resize event */
	x_resize_window(ev->xconfigure.width, ev->xconfigure.height);
  
    return(0);
}

x_buttonwait()
{
    XEvent event;
    int ev;
    int (*event_handler)();
    long xeventmask = ButtonPressMask;

    event_handler = default_event_function;
    
    while (1)
    {
        ev = x_ask_event(&event);
	if (event_handler(&event))
	{
	    return(0);
	}
    }
}

x_number_of_colors()
{
    int planes;
    double pow();

    planes = DefaultDepth(x_dpy, DefaultScreen(x_dpy));

    return((int)(pow(2.0,planes)));
}

int x_calc_size_in_pixel(size, axis)
real size;                          
char axis;
{
    if(axis == 'x')
        return((int)(x_xfactor*size + .5));
    else if (axis == 'y')
        return((int)(x_yfactor*size + .5));
    else
        return(0);
}

/* very new stuff, handle with care ... */

x_screendump(dump_name)
char *dump_name;
{
    FILE *fptr, *fopen();

#if DEBUG
    printf("dumping screen ... ");
#endif
    fptr = fopen(dump_name, "w");
    if(fptr == NULL)
    {
       fprintf(stderr,"cannot open file %s\n", dump_name);
       exit(1);
    }
    Window_Dump(fptr);
#if DEBUG
    printf("ready !\n");
#endif
}
/*
 * Window_Dump: dump a window to a file which must already be open for
 *              writting.
 */
#ifndef alpha
extern char *calloc();
#endif
typedef unsigned long Pixel;

#include "X11/XWDFile.h"

static Window_Dump(out)
     FILE *out;
{
    unsigned long swaptest = 1;
    XColor *colors;
    unsigned buffer_size;
    int win_name_size;
    int header_size;
    int ncolors, i;
    char *win_name;
    Bool got_win_name;
    XWindowAttributes win_info;
    XImage *image;
    int absx, absy, x, y;
    unsigned width, height;
    int dwidth, dheight;
    int bw;
    Window dummywin;
    XWDFileHeader xwdheader;
    
    
    /*
     * Get the parameters of the window being dumped.
     */
    if(!XGetWindowAttributes(x_dpy, x_win, &win_info))
    {
	fprintf(stderr,"Can't get target window attributes.");
	exit(1);
    }

    /* handle any frame window */
    if (!XTranslateCoordinates (x_dpy, x_win, RootWindow (x_dpy, DefaultScreen(x_dpy)), 0, 0,
				&absx, &absy, &dummywin)) {
	fprintf (stderr, 
		 "unable to translate window coordinates (%d,%d)\n",
		 absx, absy);
	exit (1);
    }
    win_info.x = absx;
    win_info.y = absy;
    width = win_info.width;
    height = win_info.height;
    bw = 0;

    dwidth = DisplayWidth (x_dpy, DefaultScreen(x_dpy));
    dheight = DisplayHeight (x_dpy, DefaultScreen(x_dpy));


    /* clip to window */
    if (absx < 0) width += absx, absx = 0;
    if (absy < 0) height += absy, absy = 0;
    if (absx + width > dwidth) width = dwidth - absx;
    if (absy + height > dheight) height = dheight - absy;

    XFetchName(x_dpy, x_win, &win_name);
    if (!win_name || !win_name[0]) {
	win_name = "xwdump";
	got_win_name = False;
    } else {
	got_win_name = True;
    }

    /* sizeof(char) is included for the null string terminator. */
    win_name_size = strlen(win_name) + sizeof(char);

    /*
     * Snarf the pixmap with XGetImage.
     */

    x = absx - win_info.x;
    y = absy - win_info.y;

    image = XGetImage (x_dpy, x_win, x, y, width, height, AllPlanes, ZPixmap);

    if (!image) {
	fprintf (stderr, "unable to get image at %dx%d+%d+%d\n",
		 width, height, x, y);
	exit (1);
    }

    /*
     * Determine the pixmap size.
     */
    buffer_size = Image_Size(image);

    ncolors = Get_XColors(&win_info, &colors);

    XFlush(x_dpy);

    /*
     * Calculate header size.
     */
    header_size = sizeof(xwdheader) + win_name_size;

    /*
     * Write out header information.
     */
    xwdheader.header_size = (CARD32) header_size;
    xwdheader.file_version = (CARD32) XWD_FILE_VERSION;
    xwdheader.pixmap_format = (CARD32) ZPixmap;
    xwdheader.pixmap_depth = (CARD32) image->depth;
    xwdheader.pixmap_width = (CARD32) image->width;
    xwdheader.pixmap_height = (CARD32) image->height;
    xwdheader.xoffset = (CARD32) image->xoffset;
    xwdheader.byte_order = (CARD32) image->byte_order;
    xwdheader.bitmap_unit = (CARD32) image->bitmap_unit;
    xwdheader.bitmap_bit_order = (CARD32) image->bitmap_bit_order;
    xwdheader.bitmap_pad = (CARD32) image->bitmap_pad;
    xwdheader.bits_per_pixel = (CARD32) image->bits_per_pixel;
    xwdheader.bytes_per_line = (CARD32) image->bytes_per_line;
    xwdheader.visual_class = (CARD32) win_info.visual->class;
    xwdheader.red_mask = (CARD32) win_info.visual->red_mask;
    xwdheader.green_mask = (CARD32) win_info.visual->green_mask;
    xwdheader.blue_mask = (CARD32) win_info.visual->blue_mask;
    xwdheader.bits_per_rgb = (CARD32) win_info.visual->bits_per_rgb;
    xwdheader.colormap_entries = (CARD32) win_info.visual->map_entries;
    xwdheader.ncolors = ncolors;
    xwdheader.window_width = (CARD32) win_info.width;
    xwdheader.window_height = (CARD32) win_info.height;
    xwdheader.window_x = absx;
    xwdheader.window_y = absy;
    xwdheader.window_bdrwidth = (CARD32) win_info.border_width;

    if (*(char *) &swaptest) {
	_swaplong((char *) &xwdheader, sizeof(xwdheader));
	for (i = 0; i < ncolors; i++) {
	    _swaplong((char *) &colors[i].pixel, sizeof(long));
	    _swapshort((char *) &colors[i].red, 3 * sizeof(short));
	}
    }

    (void) fwrite((char *)&xwdheader, sizeof(xwdheader), 1, out);
    (void) fwrite(win_name, win_name_size, 1, out);

    /*
     * Write out the color maps, if any
     */

    (void) fwrite((char *) colors, sizeof(XColor), ncolors, out);

    /*
     * Write out the buffer.
     */

    /*
     *    This copying of the bit stream (data) to a file is to be replaced
     *  by an Xlib call which hasn't been written yet.  It is not clear
     *  what other functions of xwd will be taken over by this (as yet)
     *  non-existant X function.
     */
    (void) fwrite(image->data, (int) buffer_size, 1, out);

    /*
     * free the color buffer.
     */

    if(ncolors > 0) free(colors);

    /*
     * Free window name string.
     */
    if (got_win_name) XFree(win_name);

    /*
     * Free image
     */
    XDestroyImage(image);

    fclose(out);
}

/*
 * Determine the pixmap size.
 */

static int Image_Size(image)
     XImage *image;
{
    if (image->format != ZPixmap)
      return(image->bytes_per_line * image->height * image->depth);

    return(image->bytes_per_line * image->height);
}

#define lowbit(x) ((x) & (~(x) + 1))

/*
 * Get the XColors of all pixels in image - returns # of colors
 */
static int Get_XColors(win_info, colors)
     XWindowAttributes *win_info;
     XColor **colors;
{
    int i, ncolors;
    Colormap cmap = win_info->colormap;
    Pixel red, green, blue, red1, green1, blue1;

    if (!cmap)
	return(0);

    ncolors = win_info->visual->map_entries;
    if (!(*colors = (XColor *) malloc (sizeof(XColor) * ncolors)))
    {
	fprintf(stderr,"Out of memory!");
	exit(1);
    }

    if ((win_info->visual->class == DirectColor) ||
	(win_info->visual->class == TrueColor)) {
	red = green = blue = 0;
	red1 = lowbit(win_info->visual->red_mask);
	green1 = lowbit(win_info->visual->green_mask);
	blue1 = lowbit(win_info->visual->blue_mask);
	for (i=0; i<ncolors; i++) {
	  (*colors)[i].pixel = red|green|blue;
	  (*colors)[i].pad = 0;
	  red += red1;
	  if (red > win_info->visual->red_mask)
	    red = 0;
	  green += green1;
	  if (green > win_info->visual->green_mask)
	    green = 0;
	  blue += blue1;
	  if (blue > win_info->visual->blue_mask)
	    blue = 0;
	}
    } else {
	for (i=0; i<ncolors; i++) {
	  (*colors)[i].pixel = i;
	  (*colors)[i].pad = 0;
	}
    }

    XQueryColors(x_dpy, cmap, *colors, ncolors);
    
    return(ncolors);
}

static _swapshort (bp, n)
    register char *bp;
    register unsigned n;
{
    register char c;
    register char *ep = bp + n;

    while (bp < ep) {
	c = *bp;
	*bp = *(bp + 1);
	bp++;
	*bp++ = c;
    }
}

static _swaplong (bp, n)
    register char *bp;
    register unsigned n;
{
    register char c;
    register char *ep = bp + n;
    register char *sp;

    while (bp < ep) {
	sp = bp + 3;
	c = *sp;
	*sp = *bp;
	*bp++ = c;
	sp = bp + 1;
	c = *sp;
	*sp = *bp;
	*bp++ = c;
	bp += 2;
    }
}














/*
 *   X.c - Part of my_x and yapp_x (NEMO) graphics driver programs
 *         Provides generic X11-based plotting functions
 *
 *   Author:     Patrick (AdM) Frisch ( frisch@walhall.uni-sw.gwdg.de )
 *   Date:       sometimes in 1992
 *   Revised:    March 93
 *   Version:    1.01
 *   Module
 *   Remarks:    Provides basic window stuff, not very elegant, but it
 *               seems to work quite stable
 *
 *   15-jan-95  double -> real for yapp		PJT
 *
 *   Global Functions:
 *
 *
 *               init_X(int, char **)
 *               close_X(void)
 *               nameToPixel(char *,unsignedlong)  - a bit misplaced here ???
 *               setgeom(int, int)
 *
 */

#include "config.h"
#include <X11/Xlib.h>
#include <X11/Intrinsic.h>
#include <X11/Xatom.h>

#define DEFAULTFONT "fixed"

extern char *getenv();

char *x_bgname = NULL;
char *x_displayname = NULL;
char *x_fgname = NULL;
char *x_fontname = NULL;
char *x_geom = NULL;
int x_reverse=0;
int x_filling=0;
int x_no_resize=1;

Display *x_dpy;
Window x_win;
Pixmap x_pixmap;
GC x_gc;
GC fillgc;
XGCValues x_gcv;
XWindowAttributes x_wa;
XFontStruct *x_font = NULL;
int xeventmask = 0;

int rwidth = 500, rheight = 500;
int xwin = 0, ywin = 0;
long x_geomflag;

int fgcolor = 1; /* Default pixel values... maybe not desirable? */
int bgcolor = 0;
int initfg, initbg;
XFontStruct *initfont;

local getgeom(string);
local usage(void);

init_X(argc, argv)
int argc;
char *argv[];
{
    int i, tmpcolor;
    Window rwindow;
    XSetWindowAttributes swa;
    XSizeHints      sh;
    int reverse_denied = 0;
#if 1
    char *ident = "my_x - another fine product from AdM - Software";
#else
    char ident[256];
    char *getparam();

    sprintf(ident,"%s: yapp_x",getparam("argv0"));
#endif


    /* first handle arguments, beginning with 0 against common
       rules because program name is not included in this argv[0] */

    for (i=0; i < argc; i++)
    {
	/* Not using getopt: Be somewhat inflexible... */
	if (argv[i][0] == '+') {
	    if ((!strcmp(argv[i],"+rv")) ||
		(!strcmp(argv[i], "+reverse"))) {
		reverse_denied = 1;
		continue;
	    }
	}
	if (argv[i][0] == '-') {
	    if ((!strcmp(argv[i],"-bg")) ||
		(!strcmp(argv[i], "-background"))) {
		if ((x_bgname = argv[++i]) == NULL) {
		    usage();
		    exit(-1);
		}
		continue; /* Next arg */
	    }
	    if ((!strcmp(argv[i], "-d")) ||
		(!strcmp(argv[i],"-display"))) {
		/* Set display */
		x_displayname = argv[++i];
		continue;
	    }
	    if ((!strcmp(argv[i],"-fg")) ||
		(!strcmp(argv[i], "-foreground"))) {
		if ((x_fgname = argv[++i]) == NULL) {
		    usage();
		    exit(-1);
		}
		continue; /* Next arg */
	    }
	    if ((!strcmp(argv[i], "-fn")) ||
		(!strcmp(argv[i], "-font"))) {
		if ((x_fontname = argv[++i]) == NULL) {
		    usage();
		    exit(-1);
		}
		continue; /* Next arg */
	    }
	    if ((!strcmp(argv[i], "-g")) ||
		(!strcmp(argv[i],"-geometry"))) {
		if ((x_geom = argv[++i]) == NULL) {
		    usage();
		    exit(-1);
		}
		continue; /* Next arg */
	    }
	    if ((!strcmp(argv[i],"-rv")) ||
		(!strcmp(argv[i], "-reverse"))) {
		/* Reverse bg and fg */
		x_reverse = 1;
		continue;
	    }
	    /* Default */
	    usage();
	    exit(0);
	}
    }
    
    if (x_displayname == NULL)
        x_displayname = getenv("DISPLAY");

    if ((x_dpy = XOpenDisplay(x_displayname)) == NULL) {
	fprintf(stderr, "Can't open display %s\n",
	        XDisplayName(x_displayname));
	exit(1);
    }

    rwindow = XDefaultRootWindow(x_dpy);

    XGetWindowAttributes(x_dpy, rwindow, &x_wa);
    if (x_wa.colormap == 0)
        x_wa.colormap = DefaultColormap(x_dpy, DefaultScreen(x_dpy));

    if (x_geom == NULL) {
	if (x_geom = XGetDefault(x_dpy,"yapp","geometry")) {
	    getgeom(x_geom);
	} /* Otherwise just use our default geometry */
    } else
        getgeom(x_geom);

    if (x_fgname == NULL) {
	if (x_fgname = XGetDefault(x_dpy, "yapp", "foreground")){
	    nameToPixel(x_fgname, &fgcolor);
	}
    } else
        nameToPixel(x_fgname, &fgcolor);

    if (x_bgname == NULL) {
	if (x_bgname = XGetDefault(x_dpy, "yapp", "background")){
	    nameToPixel(x_bgname, &bgcolor);
	}
    } else
        nameToPixel(x_bgname,&bgcolor);

    if (!x_reverse && !reverse_denied)
        if (XGetDefault(x_dpy,"yapp","reverseVideo"))
	    x_reverse = 1;
    
/* my own font routine needs initialized gc !! */
    
    x_pixmap = XCreatePixmap(x_dpy,DefaultRootWindow(x_dpy), rwidth, rheight,
			   DefaultDepth(x_dpy, DefaultScreen(x_dpy)));
    x_gc = XCreateGC(x_dpy, x_pixmap, 0L, &x_gcv);
    
    if (x_fontname == NULL) {
	if (x_fontname = XGetDefault(x_dpy, "yapp", "font")){
	    if ((x_font = XLoadQueryFont(x_dpy, x_fontname)) == NULL) {
		/* try own routine */
		if (!x_font_by_name(x_fontname))
		{
		    fprintf(stderr, "Font %s does not exist\n",
			    x_fontname);
		    exit(0);
		}
	    }
	} else {
	    if ((x_font = XLoadQueryFont(x_dpy, DEFAULTFONT)) == NULL) {
		/* try own routine */
		if (!x_font_by_name(x_fontname))
		{
		    fprintf(stderr, "Font %s does not exist\n",
			    x_fontname);
		    exit(0);
		}
	    }
	}
    } else {
	if ((x_font = XLoadQueryFont(x_dpy, x_fontname)) == NULL) {
	    /* try own routine */
	    if (!x_font_by_name(x_fontname))
	    {
		fprintf(stderr, "Font %s does not exist\n",
			x_fontname);
		exit(0);
	    }
	}
    }
/*    
    x_pixmap = XCreatePixmap(x_dpy,DefaultRootWindow(x_dpy), rwidth, rheight,
			   DefaultDepth(x_dpy, DefaultScreen(x_dpy)));
    x_gc = XCreateGC(x_dpy, x_pixmap, 0L, &x_gcv);
*/
/* not necessary, if font_by_name is used */
    
    XSetFont(x_dpy, x_gc, x_font->fid); /* Font must be set by now */
    if (x_reverse){
	tmpcolor = fgcolor;
	fgcolor = bgcolor;
	bgcolor = tmpcolor;
    }
    XSetForeground(x_dpy,x_gc,fgcolor);
    XSetBackground(x_dpy,x_gc,bgcolor);

    initfg = fgcolor; /* Need these for proper resizing */
    initbg = bgcolor;
    initfont = x_font;
	
    x_erase();
    swa.background_pixmap = x_pixmap;
    x_win = XCreateWindow(x_dpy, rwindow, (unsigned int) xwin,
			(unsigned int) ywin, (unsigned int) rwidth,
			(unsigned int) rheight, 0, CopyFromParent,
			InputOutput, CopyFromParent,
			CWBackPixmap, &swa);
    XChangeProperty(x_dpy, x_win, XA_WM_NAME, XA_STRING, 8,
		    PropModeReplace, ident, strlen(ident));
		    
    sh.width = rwidth;
    sh.height = rheight;
    sh.flags = USSize;
    if ((x_geomflag&XValue)&&(x_geomflag&YValue)) {
	sh.x = xwin;
	sh.y = ywin;
	sh.flags |= USPosition;
    }
    XSetNormalHints(x_dpy, x_win, &sh);
    XMapWindow(x_dpy, x_win);
    XFlush(x_dpy);
}

setgeom(rw, rh)
int rw, rh;
{
    rwidth = rw;
    rheight = rh;
}

static getgeom(ptr)
char *ptr;
{
    x_geomflag = XParseGeometry(ptr, &xwin, &ywin,
			      (unsigned int *) &rwidth,
			      (unsigned int *) &rheight);
    if (x_geomflag&XNegative) {
	xwin = DisplayWidth(x_dpy, DefaultScreen(x_dpy)) - xwin - rwidth;
    }
    if (x_geomflag&YNegative) {
	ywin = DisplayHeight(x_dpy, DefaultScreen(x_dpy)) - ywin -rheight;
    }
}

int nameToPixel(name, pixel)
char *name;
unsigned long *pixel;
{
    XColor color;

    if (!XParseColor(x_dpy, x_wa.colormap, name, &color)) {
#ifdef DEBUG
	fprintf(stderr, "Unknown color '%s'\n", name);
#endif
	return 0;
    }
    if (!XAllocColor(x_dpy, x_wa.colormap, &color)) {
#ifdef DEBUG
	fprintf(stderr, "Cannot allocate color '%s'\n", name);
#endif
	return 0;
    }
    *pixel= color.pixel;
    return 1;
}

close_X()
{
    XFreeGC(x_dpy,x_gc);
    XFreePixmap(x_dpy,x_pixmap);
    XDestroyWindow(x_dpy,x_win);

    if(!x_no_resize)
	x_close_plobj();
}

x_resize_window(w_new, h_new)
int w_new, h_new;
{
    extern real x_xmin, x_ymin, x_xmax, x_ymax;

    if(x_no_resize)
	return;
    
    if((w_new==rwidth)&&(h_new==rheight))
	return;

    rwidth = w_new;
    rheight = h_new;

    XFreePixmap(x_dpy,x_pixmap);
    x_pixmap = XCreatePixmap(x_dpy,DefaultRootWindow(x_dpy), rwidth, rheight,
			   DefaultDepth(x_dpy, DefaultScreen(x_dpy)));
    x_clear_window();
    XSetWindowBackgroundPixmap(x_dpy, x_win, x_pixmap);
    /* set up new scaling factors */
    x_space(x_xmin, x_ymin, x_xmax, x_ymax);
    x_replot();
}

static usage()
{
    fprintf(stderr, "-bg colorname [background color]\n");
    fprintf(stderr, "-display dispname\n");
    fprintf(stderr, "-fg colorname [foreground color]\n");
    fprintf(stderr, "-fn fontname [font for labels]\n");
    fprintf(stderr, "-geometry geom [Standard Xgeometry string]\n");
    fprintf(stderr, "-rv [reverse foreground and background]\n");
}


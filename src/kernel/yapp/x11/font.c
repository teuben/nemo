/*
 *   font.c - Part of my_x and yapp_x (NEMO) graphics driver programs
 *         Provides generic X11-based plotting functions
 *         
 *   Author:     Patrick (AdM) Frisch ( frisch@walhall.uni-sw.gwdg.de )
 *   Date:       sometimes in 1992
 *   Revised:    March 93
 *   Version:    1.01
 *   Module
 *   Remarks:    font handling
 *
 *   15-jan-95  changed double -> real				pjt
 *
 *   Global Functions:
 *
 *               x_alabel(char, char, char *, real)
 *               x_font_by_name(char *)
 *               x_fontsize(real)
 *               set_font_and_size(char *, int) (*do not use this one !*)
 *
 */


#include "config.h"
#include <X11/Xlib.h>

static real font_size = 14.;
static char *font_name = NULL;
static int font_name_length = 0;
static char *user_font_name = NULL;    /* user specified font name */
static int user_font_name_length = 0;
static char *x_name = NULL;        /* x font name */
static int x_name_length = 0;

extern Display *x_dpy;
extern Window x_win;
extern Pixmap x_pixmap;
extern GC x_gc;
extern XFontStruct *x_font;
extern XGCValues x_gcv;
extern real x_curpt[];
extern real x_xmin, x_ymax;
extern real x_xfactor, x_yfactor;

#define DEFAULTFONT "fixed"

#define   RAD            0.017453293

static struct name_info_struct
{
      char *plot;
      char *x;
};

static struct name_info_struct name_info[] = {
    {
	    "courier-bold",
	        "courier-bold-r"},
    {
	    "courier-boldoblique",
	        "courier-bold-o"},
    {
	    "courier-oblique",
	        "courier-medium-o"},
    {
	    "courier",
	        "courier-medium-r"},
    {
	    "helvetica-bold",
	        "helvetica-bold-r"},
    {
	    "helvetica-boldoblique",
	        "helvetica-bold-o"},
    {
	    "helvetica-oblique",
	        "helvetica-medium-o"},
    {
	    "helvetica",
	        "helvetica-medium-r"},
    {
	    "symbol",
	        "symbol-medium-r"},
    {
	    "times-bold",
	        "times-bold-r"},
    {
	    "times-bolditalic",
	        "times-bold-i"},
  {
          "times-italic",
	      "times-medium-i"},
    {
	    "times-roman",
	        "times-medium-r"},
    {
	    NULL,
	        NULL}
};

/* font structures necessary for finding the bounding boxes */
static int direction_return;
static int font_ascent_return, font_descent_return;
static XCharStruct overall_return;

/* locals */
local int text_rotate(string, int, int, real, int, int);
local int transform_pixel(int *, int *, real, real);

int
x_alabel (x_justify, y_justify, s, angle)
         char x_justify, y_justify;
         char *s;
         real angle;
{
    int alabel_width, alabel_height, alabel_decent;
    int xpos, ypos,hot_x, hot_y;
    real x_char_offset, y_char_offset, hot_x_offset, hot_y_offset;
    real sphi, cphi;

    switch (x_justify)        /* placement of label with respect
				 to x coordinate */
   {
	case 'l':              /* left justified */
	x_char_offset = 0.0;
	hot_x_offset = 0.0;
	break;
	
	case 'c':              /* centered */
	x_char_offset = -0.5;
	hot_x_offset = 0.5;
	break;
	
	case 'r':              /* right justified */
	x_char_offset = -1.0;
	hot_x_offset = 1.0;
	break;
	
	default:               /* None of the above? */
	x_char_offset = 0.0;
	hot_x_offset = 0.0;
#ifdef DEBUG
	warning("x_alabel: Unknown x justification: %c", x_justify);
#endif
	break;
    }
    
    switch (y_justify)        /* placement of label with respect
				 to y coordinate */
    {
	case 'b':              /* bottom */
	y_char_offset = -1.0;
	hot_y_offset = 1.0;
	break;
	
	case 'c':              /* centered */
	y_char_offset = -0.5;
	hot_y_offset = 0.5;
	break;

	case 't':              /* top */
	y_char_offset = 0.0;
	hot_y_offset = 0.0;
	break;

	default:            /* None of the above? */
	y_char_offset = 0.0;
	hot_y_offset = 0.0;
#ifdef DEBUG
	warning("x_alabel: Unknown y justification: %c", y_justify);
#endif
	break;
    }
    alabel_width = XTextWidth (x_font, s, strlen (s));
    XTextExtents (x_font, s, strlen (s), &direction_return,
		  &font_ascent_return, &font_descent_return, &overall_return);
    alabel_height = font_ascent_return + font_descent_return;
    alabel_decent = font_descent_return;
    
    xpos =  (int)((x_curpt[0] - x_xmin) * x_xfactor
		  + x_char_offset * alabel_width);
    ypos =  (int)((x_ymax - x_curpt[1]) * x_yfactor
		  + y_char_offset * alabel_height);

    hot_x = (int)(hot_x_offset*alabel_width);
    hot_y = (int)(hot_y_offset*alabel_height);
#ifdef DEBUG
    warning("x_alabel: drawstring %d %d %d %s", x_gc, xpos, ypos, s);
#endif
/*
    XDrawString (x_dpy, x_pixmap, x_gc, xpos, ypos, s, strlen (s));
*/
    text_rotate(s, xpos, ypos, angle, hot_x, hot_y);

    x_insert_object('s', (real)x_justify, (real)y_justify, angle, 0., s);
    
    return 0;
}

int set_font_and_size (name, size)
char *name;
int size;
{
    int name_size;
    int i;
    XFontStruct *newfont;
#ifdef DEBUG
        warning("set_font_and_size: fontname `%s' size %d", 
		name ? name : "", size);
#endif

    /* If the name is null or zero length, don't change the font */
    
    if ((name == NULL)
	|| (strlen(name) == 0))
    {
	return 0;
    }
    /* If size = 0 then name should be an X font name */
    if (size == 0)
    {
	newfont = XLoadQueryFont (x_dpy, name);
	if (newfont != NULL)       /* if the font does exist */
	{
	    if (x_font != NULL)
	        XFreeFont (x_dpy, x_font);
	    x_font = newfont;
	    XSetFont(x_dpy, x_gc, x_font->fid);
#ifdef DEBUG
	    warning("set_font_and_size: Font `%s' %d", name, x_font);
#endif
	    return 1;          /* font found */
	}
	return 0;          /* font not found */
    }
    
    /* Try buildin an x font name from the name and size: */
    
    /* save the user specified name for later use */
    if (user_font_name_length < strlen (name) + 256)
    {
	user_font_name_length = strlen (name) + 256;
	if (user_font_name)
	    user_font_name = (char *) realloc (user_font_name,
					   user_font_name_length);
	else
	    user_font_name = (char *) malloc (user_font_name_length);
    }
    strcpy (user_font_name, name);

    /* search for the X font name correxponding to the requested name */
    i = -1;
    while (name_info[++i].plot)
    {
	if (strcmp (name_info[i].plot, name) == 0)
	    break;
    }

    /* if the X font name is found, use it.  Otherwise, assume the supplied
       name is an X font name itself. */
    if (name_info[i].plot)
    {
	if (!x_name)
	{
	    x_name_length = strlen (name_info[i].x) + 24;
	    x_name = malloc (x_name_length);
	}
	else if (strlen (name_info[i].x) + 24 >= x_name_length)
	{
	    x_name_length = strlen (name_info[i].x) + 24;
	    x_name = realloc (x_name, x_name_length);
	}
	sprintf (x_name, "*-%s-*-%d-*", name_info[i].x, size);
	if (set_font_and_size (x_name, 0))
	    return 1;
    }
    else
    {
	/* Try using the name as an X font name: */
	if (set_font_and_size (name, 0))
	    return 1;
    }
    
    return 0; /* Font not changed */
}

int
x_font_by_name (name)
char *name;
{
    static char *oldfont = NULL;

    if(oldfont && (strcmp(oldfont, name) == 0))
        return 1;

    x_insert_object('f', 0., 0., 0., 0., name);
    
    if (!set_font_and_size (name, (int) font_size)) {
	fprintf(stderr, "warning: font `%s' not found.\n", user_font_name);
	return 0;
    }
    if(oldfont != NULL)
       free(oldfont);
    if((oldfont = malloc(strlen(name)+1))==NULL) {
	error("cannot allocate ...");
    }
    strcpy(oldfont, name);
    
    return 1;
}

int
x_fontsize (rpoints)
real rpoints;
{
    static real oldsize = 0.0;
    int isize, iup, idown;

    if(oldsize == rpoints)
        return 1;

    x_insert_object('f', rpoints, 0., 0., 0., NULL);
    
    font_size = x_xfactor * rpoints;

#if 0 /* old method, look up a size and skip if it is not found */
    if (!set_font_and_size (user_font_name, (int) font_size))
    {
#ifdef DEBUG
	fprintf(stderr, "warning: font `%s' at size %2f not found.\n",
		user_font_name, font_size);
#endif
	return 0;
    }
#else /* new method, try to find the neighbouring font */
    isize = iup = idown = (int) (font_size+.5); /* round !!*/
    do
    {
	if (set_font_and_size (user_font_name, idown))
	{
	    oldsize = rpoints;
	    return 1;
	}
#if DEBUG
	else
	    fprintf(stderr, "warning: font `%s' at size %d not found.\n",
		    user_font_name, idown);
#endif
	if(iup != idown)
	    if (set_font_and_size (user_font_name, iup))
	    {
	        oldsize = rpoints;
	        return 1;
	    }
#if DEBUG
	    else
	        fprintf(stderr, "warning: font `%s' at size %d not found.\n",
		        user_font_name, iup);
#endif
	idown--;
	iup++;
    }
    while((idown > 0)&&(idown > isize-5)&&(iup < isize+5)); 
#endif
    return 0;
}

/* new stuff (still in development, handle with care ) */

static int
text_rotate(s, xpos, ypos, angle, hot_x, hot_y)
char *s;
int xpos, ypos;
real angle;
int hot_x, hot_y;
{
    Pixmap textBackingPixmap = NULL;
    GC rotateGC = NULL;
    XGCValues values;
    XImage *from_image;
    int i, j, w, h;
    int direction_return;
    int font_ascent_return, font_descent_return;
    XCharStruct overall_return;
    int dx, dy, xoff, yoff;
    real sphi, cphi;

    w = XTextWidth (x_font, s, strlen (s));
    XTextExtents (x_font, s, strlen (s), &direction_return,
		  &font_ascent_return, &font_descent_return, &overall_return);
    h = font_ascent_return + font_descent_return;

    textBackingPixmap = XCreatePixmap(x_dpy, x_win, MAX(w,h), MAX(w,h),
			            DefaultDepth(x_dpy, DefaultScreen(x_dpy)));
    rotateGC = XCreateGC(x_dpy, x_win, 0L, &x_gcv);

    XSetForeground(x_dpy, rotateGC, 0);
    XFillRectangle(x_dpy, textBackingPixmap, rotateGC, 0,0,MAX(w,h),MAX(w,h));
    XSetForeground(x_dpy, rotateGC, 1);
    XSetFont(x_dpy, rotateGC, x_font->fid);

    XDrawString(x_dpy, textBackingPixmap, rotateGC, 0, font_ascent_return,
		s, strlen(s));

    from_image = XGetImage(x_dpy, textBackingPixmap, 0, 0, w, h, 1, ZPixmap);

    if(from_image == NULL)
    {
	error("text_rotate: cannot get text as image ??");
    }
    
    sphi = sin(-RAD*angle); cphi = cos(-RAD*angle);
    xoff = hot_x; yoff = hot_y;
    transform_pixel(&hot_x, &hot_y, sphi, cphi);
    xoff -= hot_x; yoff -= hot_y;

    for(i = 0; i < w; i++)
        for(j = 0; j < h; j++)
	    if(XGetPixel(from_image, i, j) == 1)
	    {
		dx = i; dy = j;
		transform_pixel(&dx, &dy, sphi, cphi);
	        XDrawPoint(x_dpy, x_pixmap, x_gc, xpos+dx+xoff, ypos+dy+yoff);
	    }
    
    XFreeGC(x_dpy, rotateGC);
    XFreePixmap(x_dpy, textBackingPixmap);
}

static int
transform_pixel(x, y, sphi, cphi)
int *x, *y;
real sphi, cphi;
{
    register real xi, yi;

    xi = (real) *x; yi = (real) *y;

    *x = (int)(xi*cphi - yi*sphi + .5);
    *y = (int)(xi*sphi + yi*cphi + .5);
}












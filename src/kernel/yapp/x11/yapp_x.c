/***********************************************************
 ** YAPP_X     - a NEMO-driver for X11 R5 (should be run ***
 **              on other versions too).                 ***
 ***********************************************************
 ** Based on xplot from standard MIT distribution        ***
 ** Author : Patrick (AdM) Frisch                        ***
 **          Department of Earth Science & Astronomy     ***
 **          College of Arts & Sciences                  ***
 **          University of Tokyo                         ***
 **          e-mail: patrick@kyohou.c.u-tokyo.ac.jp      ***
 **                                                      ***
 ** Version: 1.0                                         ***
 ** Date   : 10/23/92                                    ***
 ** Remarks:                                             ***
 **          10/23/92 initial coding and debugging (AdM) ***
 **          10/26/92 adjust plcolor for snapplot (AdM)  ***
 **          10/30/92 rotating text implemented (AdM)    ***
 **          10/31/92 code cleaning (AdM)                ***
 **          11/05/92 bug in plcircle (AdM)              ***
 **          11/06/92 fiddling with more colors (PGPLOT) ***
 **                   by PJT                             ***
 **          11/07/92 color order changed (AdM)          ***
 **          11/10/92 global vars renamed +              ***
 **                   code cleaning (AdM)                ***
 **          11/15/92 reverse video now foolproof        ***
 **                   (hopely)+screendump implemented    ***
 **          08/26/93 using pltdev in the correct way    ***
 **          01/14/95 made it work with 'real' double    ***
 **          05/05/95 no more special first plframe      ***
 **                                                      ***
 **                                                      ***
 ** Many thank's to:                                     ***
 **                   Authors of xplot (unknown)         ***
 **                   Jun Makino (for discussion, hints, ***
 **                               ...)                   ***
 **********************************************************/

#include "config.h"
#include <signal.h>

#define PGPLOT

/* line styles  */

#define SOLID          0
#define DOTTED         1
#define SHORTDASHED    2
#define DOTDASHED      3
#define LONGDASHED     4
#define DISCON         5

/* colors */

#if 1
	/* PGPLOT colors */
#define BLACK           0
#define RED             1
#define ORANGE          2
#define YELLOW          3       /* (Red + Green) */
#define GREEN_YELLOW    4
#define GREEN           5
#define GREEN_CYAN      6
#define CYAN            7       /* (Green + Blue) */
#define BLUE            8
#define BLUE_CYAN       9
#define RED_MAGENTA     10
#define MAGENTA         11      /* (Red + Blue) */  
#define BLUE_MAGENTA    12
#define DARK_GRAY       13
#define LIGHT_GRAY      14
#define WHITE           15


#else

#define BLACK          0
#define RED            1
#define YELLOW         2
#define GREEN          3
#define CYAN           4
#define BLUE           5
#define MAGENTA        6
#define WHITE          7

#endif

#define MAXCOLORS      256     /* 2^8 , e.g. 8 Biplanes maximum */

/* local used variables */

local int ncolors;

local real red[MAXCOLORS],
           green[MAXCOLORS],
           blue[MAXCOLORS];

local cp_str2arg(string, int *, string **);
local rgb_setup(void);
/*
 * PLINIT: initalize the plotting package.
 */

plinit(pltdev, xmin, xmax, ymin, ymax)
char *pltdev;		/* output device name (ignored) */
real xmin, xmax;	/* user plotting area */
real ymin, ymax;
{
    int i, my_argc = 0;
    char **my_argv;
    extern char *yapp_string;
    extern int fgcolor,bgcolor;
    
    dprintf(1, "**************************\n");
    dprintf(1, "***    YAPP FOR X     ****\n");
    dprintf(1, "**************************\n\n");
    dprintf(1, "   written (92) by AdM\n\n");

    if(strcmp(pltdev,"***"))
      cp_str2arg(pltdev, &my_argc, &my_argv);
    else if(yapp_string != NULL)
      cp_str2arg(yapp_string, &my_argc, &my_argv);

    init_X(my_argc, my_argv);
    x_space(xmin, ymin, xmax, ymax);
    rgb_setup();

    dprintf(2,"%d %d\n",fgcolor,bgcolor);
}

static cp_str2arg(s, argc, argv)
char *s;
int *argc;
char ***argv;
{
    int numargs=1;
    char *tmp;
    char **argtmp;
    
    /* simply count whitespaces */
    tmp = s;
    while(*tmp)
      if(*tmp++ == ' ') numargs+=1;

    dprintf(2, "numargs = %d\n", numargs);
    *argc = numargs;
    
    /* allocate space for pointers */
    *argv = argtmp = (char **)malloc(sizeof(char **)*numargs);

    /* copy pointers into argv by simply exchange space by '\n' */
    while(numargs--)
    {
	dprintf(2, "Count: %d s= %s ", numargs, s);
	tmp = s;
	while(*tmp && (*tmp++ != ' '));
	if(*(tmp-1) == ' ')
	   *(tmp-1) = '\0';
	*argtmp++ = s;
	dprintf(2, " arg= %s\n", s);
	s = tmp;
    }
}

static rgb_setup()
{
#if 1
  ncolors = WHITE+1;
  /* default PGPLOT  colors: although, defined here, plpalette is not called */
  red[BLACK]=0.00;         green[BLACK]=0.00;          blue[BLACK]=0.00;  
  red[WHITE]=1.00;         green[WHITE]=1.00;          blue[WHITE]=1.00;  
  red[RED]=1.00;           green[RED]=0.00;            blue[RED]=0.00;    
  red[GREEN]=0.00;         green[GREEN]=1.00;          blue[GREEN]=0.00;  
  red[BLUE]=0.00;          green[BLUE]=0.00;           blue[BLUE]=1.00;   
  red[CYAN]=0.00;          green[CYAN]=1.00;           blue[CYAN]=1.00;   
  red[MAGENTA]=1.00;       green[MAGENTA]=0.00;        blue[MAGENTA]=1.00;
  red[YELLOW]=1.00;        green[YELLOW]=1.00;         blue[YELLOW]=0.00;
  red[ORANGE]=1.00;        green[ORANGE]=0.50;         blue[ORANGE]=0.00;
  red[GREEN_YELLOW]=0.50;  green[GREEN_YELLOW]=1.00;   blue[GREEN_YELLOW]=0.00;
  red[GREEN_CYAN]=0.00;    green[GREEN_CYAN]=1.00;     blue[GREEN_CYAN]=0.50;  
  red[BLUE_CYAN]=0.00;     green[BLUE_CYAN]=0.50;      blue[BLUE_CYAN]=1.00;   
  red[BLUE_MAGENTA]=0.50;  green[BLUE_MAGENTA]=0.00;   blue[BLUE_MAGENTA]=1.00;
  red[RED_MAGENTA]=1.00;   green[RED_MAGENTA]=0.00;    blue[RED_MAGENTA]=0.50;
  red[DARK_GRAY] =0.33;    green[DARK_GRAY]=0.33;      blue[DARK_GRAY]=0.33; 
  red[LIGHT_GRAY]=0.66;    green[LIGHT_GRAY]=0.66;     blue[LIGHT_GRAY]=0.66;
#else
    ncolors = WHITE+1;
    red[BLACK]    = 0.0;  green[BLACK]    = 0.0;   blue[BLACK]     = 0.0;
    red[RED]      = 1.0;  green[RED]      = 0.0;   blue[RED]       = 0.0;
    red[YELLOW]   = 1.0;  green[YELLOW]   = 1.0;   blue[YELLOW]    = 0.0;
    red[GREEN]    = 0.0;  green[GREEN]    = 1.0;   blue[GREEN]     = 0.0;
    red[CYAN]     = 0.0;  green[CYAN]     = 0.8;   blue[CYAN]      = 0.8;
    red[BLUE]     = 0.2;  green[BLUE]     = 0.6;   blue[BLUE]      = 1.0;
    red[MAGENTA]  = 1.0;  green[MAGENTA]  = 0.4;   blue[MAGENTA]   = 1.0;
    red[WHITE]    = 1.0;  green[WHITE]    = 1.0;   blue[WHITE]     = 1.0;
#endif
}

/*
 * PLSWAP: does nothing.
 */

plswap() { }

/*
 * PLXSCALE, PLYSCALE: transform from usr to plotting coordinates.
 * At present, these do nothing (identity transformation).
 */

real plxscale(x, y)
real x, y;		/* user coordinates */
{
    return (x);
}

real plyscale(x, y)
real x, y;		/* user coordinates */
{
    return (y);
}

/*
 * PLLTYPE: select line width and dot-dash pattern.
 */

plltype(lwid, lpat)
int lwid;		/* line width */
int lpat;		/* line pattern */
{
    static int oldlwid=0;

    if(lwid <= 0)
        lwid = oldlwid;
    switch(lpat)
    {
	case SOLID:       x_linemod("solid",lwid);        break;
	case DOTTED:      x_linemod("dotted",lwid);       break;
	case SHORTDASHED: x_linemod("shortdashed",lwid);  break;
	case DOTDASHED:   x_linemod("dotdashed",lwid);    break;
	case LONGDASHED:  x_linemod("longdashed",lwid);   break;
	case DISCON:      x_linemod("disconnected",lwid); break;
	default:
	                  x_linemod("solid",lwid);
    }
    oldlwid = lwid;
}

/*
 * PLLINE, PLMOVE, PLPOINT: plot lines, moves, and points.
 */

real cur_x, cur_y;

plline(x, y)
real x, y;		/* user coordinates */
{
    x_drawto(x, y);
}

plmove(x, y)
real x, y;		/* user coordinates */
{
    x_move(x,y);
}

plpoint(x, y)
real x, y;		/* user coordinates */
{
    x_point(x, y);
}

/*
 * PLCIRCLE, PLCROSS, PLBOX: plot various symbols.
 */

plcircle(x, y, r)
real x, y;
real r;
{
    extern int x_filling;
    
    if(r == 0)
      return;
    
    if(r > 0)
      x_circle(x, y, r); /* use standard X - routine instead of self written */
    else
    {
      x_filling = 1;
      x_circle(x,y,-r);
      x_filling = 0;
  }
/*
    int npnts, i;
    real theta;

    npnts = MAX(2400 * r / dxymax, 6.0);
    plmove(x + r, y);
    for (i = 1; i <= npnts; i++) {
	theta = TWO_PI * ((real) i) / ((real) npnts);
	plline(x + r * cos(theta), y + r * sin(theta));
    }
*/
}

plcross(x, y, s)
real x, y;
real s;
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
real x, y;
real s;
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
 */

static int textjust = -1;
static char xjus = 'l',    /* default: left justified  */
            yjus = 'c';    /* default: centered        */

pljust(jus)
int jus;		/* -1, 0, 1 for left, mid, right just */
{
    textjust = (jus < -1 ? -1 : (jus > 1 ? 1 : jus));

    switch(textjust)
    {
	case -1: xjus = 'l'; break;
	case  0: xjus = 'c'; break;
	case  1: xjus = 'r'; break;
    }

    /* no idea how to calculate the vertical justification now */
    /* maybe bitwise coded in jus ???                          */
}

/*
 * PLTEXT: plot a text string.
 */

pltext(msg, x, y, hgt, ang)
char *msg;		/* message to plot, with NULL termination */
real x, y;		/* user coordinates (modified by justification) */
real hgt;		/* height of characters in user coordinates */
real ang;		/* angle of message, counterclockwise in degrees */
{
    plmove(x, y);
    x_fontsize(2.0*hgt);
    x_alabel(xjus, yjus, msg, ang); 
}

/*
 * PLFLUSH: output any pending graphics.
 */

plflush()
{
    x_doplot();
}

plerase()
{
    x_erase();
}

/*
 * PLFRAME: advance to next frame.
 */

plframe()
{
    static int first = TRUE;
    
    plflush();	/* flush graphics */
    first = FALSE;      /* no more special first time, */
    if(!first) /* most NEMO - programs send first a plframe,        */
    {          /* so we have an empty screen at first .. avoid this */
/*	dprintf(0, "\nClick Mouse Button ... "); */
	x_buttonwait();
/*	dprintf(0, " continuing ...\n"); */
	plerase();
    }
}

/*
 * PLSTOP: finish up a display.
 */

#define ONEMIN (60 * 1000 * 1000)

plstop()
{
    plframe();
/*    dprintf(0, "\nClick Mouse Button ... "); */
    x_buttonwait();
/*    dprintf(1, "ok, shutting down plot now !\n"); */
    close_X();
}

pl_matrix(frame,nx,ny,xmin,ymin,cell,fmin,fmax,findex,blankval)
real *frame, xmin, ymin, cell, fmin, fmax, findex, blankval;
int nx, ny;
{
    warning("pl_matrix: unsupported in yapp_x");
    return(0);
}

pl_screendump(fname)
string fname;
{
#if 0
    string command[128];
    
    dprintf(1, "pl_screendump: only for testing reasons !!!\n");
    dprintf(1, "\n Click in the plotting window, three beeps should appear,\n\
               after that, the picture should be saved !\n");

    plflush();
    sprintf(command,"xwd > %s", fname);
    system(command);
#else
    x_screendump(fname);
#endif    
}

/*
 *  specify color directly by it's name in rgb.txt (Standard X11)
 */

pl_color_by_name(color_name)
char *color_name;
{
    x_color_by_name(color_name);
}


/*
 *  specify plotting colors by RGB numbers in the range 0.0 ... 1.0
 */

#define RANGE(x) !(((x)>=0.0)&&((x)<=1.0))

plcolor_rgb(r,g,b)
real r,g,b;
{
    unsigned short ri,gi,bi;

    if(RANGE(r) || RANGE(g) || RANGE(b)) {
	warning("(yapp) colors not in range !!! Colors unchanged");
    } else {
	ri = (unsigned short) (r*65535);
	gi = (unsigned short) (g*65535);
	bi = (unsigned short) (b*65535);
	
	x_color_by_number(ri, gi, bi);
    }
}

/*
 *  for compatibility reasons:
 *  specify plotting color as integer between 0 and ncolors-1;
 *  values outside this range are mapped to the nearest endpoint.
 *
 *  Addition: The most important program (for that this is written
 *  actually), sends a big number hardwired for setting the foreground
 *  color. This leads to big troubles.
 */

void plcolor(color)
int color;
{
    extern int x_reverse, initbg;
    
    if(color < 0)
        color = ((x_reverse||(initbg==0)) ? 0 : ncolors-1);
    else if (color > ncolors-1)
        color = ((x_reverse||(initbg==0)) ? ncolors-1 : 0);

    plcolor_rgb(red[color], green[color], blue[color]);
}

plncolors()
{
    return(ncolors);
}

void plpalette(r, g, b, nc)
real r[], g[], b[];
int nc;
{
    int i;

    if(nc > MAXCOLORS)
        error("plpalette: cannot define more than %d colors\n", MAXCOLORS);
    ncolors = nc;
    for(i = 0; i < ncolors; i++)
    {
	red[i]    = r[i];
	green[i]  = g[i];
	blue[i]   = b[i];
    }
    red[ncolors] = green[ncolors] = blue[ncolors] = 0.0;
#if 0
    plcolor(1);              /* black */
#else
    plcolor(32768);
#endif
}

char *pl_getfont()
{
    extern char *x_fontname;

    return(x_fontname);
}

pl_setfont(fn)
char *fn;
{
    x_font_by_name(fn);
}

int plot_to_pixel(size, axis)
real size;
char axis;
{
    return(x_calc_size_in_pixel(size, axis));
}








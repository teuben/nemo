/*
 * YAPP: Yet Another Plotting Package. - calling PLPLOT package 
 *       See also:  http://plplot.sourceforge.net/ 
 *                  https://github.com/PLplot/PLplot   (a mirror)
 *
 *	8-oct-94  Created - note: need -DYAPP_LONG since plinit() is
 *		  used by both NEMO and the PLPLOT library
 *	4-nov-94  Abandoned the -DYAPP_LONG and use c_ directly
 *                but need to declare them here ... can't use plplot.h
 *	7-feb-95  proper real=float version
 *     12-mar-97  non COLOR version didn't link properly
 *  
 *     25-oct-03  various fixed for the new CVS 5.x series of plplot
 *     11-oct-07  compile option to use the double instead of float
 *   24-dec-2019  Trying PLplot-5.15.0 (june 2019) - fixed plwid -> plwidth
 *
 *  ToDo: colors
 *        fix the --without-double requirement
 *        fix aspect ratio problems
 *        presistent plots on the screen
 */

 
#include <stdinc.h>
#include <extstring.h>
#include <yapp.h>

/*   we don't want to use NEMO's real, since we may want to mix & match */
#if 1
 typedef double pl_real;
 local char *pl_nemoreal = "nemo: pl_real = double";
#else
 typedef float  pl_real;
l ocal char *pl_nemoreal = "nemo: pl_real = float";
#endif

/* we use the c_ versions since NEMO uses the same namespace ... */
/* this means we cannot include the plplot.h header as is ...... */
/* and be warned about double vs. float mis-use                  */

extern string yapp_string;	/* a kludge, see: getparam.c */

extern void c_plinit(void);
extern void c_plsdev(string);
extern void c_plsfnam(string);
extern void c_pladv(int);
extern void c_plvpor(pl_real, pl_real, pl_real, pl_real);
extern void c_plwind(pl_real, pl_real, pl_real, pl_real);
extern void c_pljoin(pl_real, pl_real, pl_real, pl_real);
extern void c_plsym(int, pl_real *, pl_real *, int);	/* hershey */
extern void c_plpoin(int, pl_real *, pl_real *, int);	/* ascii */
extern void c_plptex(pl_real, pl_real, pl_real, pl_real, pl_real, char *);
extern void c_plschr(pl_real, pl_real);
extern void c_plend(void);
extern void c_plfontld(int);
extern void c_pllsty(int);
extern void c_plwidth(int);
extern void c_plcol0(int);
extern void c_plcol1(pl_real);
extern void c_plsori(int);

extern void c_plgver(char *);
extern void c_plgpage(int *,int *,int *,int *,int *,int *);

extern string *burststring(string, string);

/* make COLOR and PG6 now the defaults */

#define COLOR  1

local real dxymax;    /* size of user window */
local pl_real xcur=0.0, ycur=0.0;

#ifdef COLOR
#define MAXCOLOR 256
local int ncolors=0;
local pl_real red[MAXCOLOR];		/* RGB color tables		    */
local pl_real blue[MAXCOLOR];
local pl_real green[MAXCOLOR];
local int cms_rgbsetup(void);
#else
#define MAXCOLOR 0
#endif

/*
 * PLINIT: initalize the plotting package.
 */

int plinit(string pltdev, real xmin, real xmax, real ymin, real ymax)
{
    pl_real width, height, x1,x2,y1,y2, zero, one;
    int   dummy, nx, ny, npldev;
    int   nxp, nyp, xoff, yoff, xlen, ylen;
    string *pldev;
    char ver[80];

    dprintf(1,"plinit (PLPLOT): %s\n",pl_nemoreal);
    
    if (yapp_string == NULL || *yapp_string == 0)
	if (pltdev != NULL && *pltdev != 0)
	    yapp_string = pltdev;

    nx = ny = 1;        /* only one window on the page */
    if (yapp_string) {
        dprintf(1,"plinit (PLPLOT): %s\n",yapp_string);
        pldev = burststring(yapp_string,",");
        npldev = xstrlen(pldev,sizeof(string))-1;
        c_plsdev(pldev[0]);
        if (npldev>1)
            c_plsfnam(pldev[1]);
    } else
        warning("No device name specified - plplot running interactive");

    c_plsori(0);     /* 0=landscape 1=portrait upside 3=portrait ?*/
    
    c_plinit();
    c_pladv(0);
    c_plvpor(0.0, 1.0, 0.0, 1.0);
    c_plwind((pl_real)xmin, (pl_real)xmax, (pl_real)ymin, (pl_real)ymax);
    c_plfontld(1);
    c_plcol0(15);    /* white color will be the default */

    c_plgver(ver);
    dprintf(1,"Plplot library version: %s\n", ver);
    c_plgpage(&nxp,&nyp,&xlen,&ylen,&xoff,&yoff);
    dprintf(1,"Plplot: plpgpage %d %d  %d %d  %d %d\n",
            nxp,nyp,xlen,ylen,xoff,yoff);
     

    x1=xmin;  x2=xmax;
    y1=ymin;  y2=ymax;

    if (ymax - ymin < xmax - xmin) {		/* not used for now */
        dxymax = xmax - xmin;
        width = 1.0;
        height = (ymax - ymin) / dxymax;
    } else {
        dxymax = ymax - ymin;
        width = (xmax - xmin) / dxymax;
        height = 1.0;
    }
#if defined(COLOR)
    cms_rgbsetup();
#endif
}

/*
 * PLSWAP: does nothing.
 */

int plswap() { }

/*
 * PLXSCALE, PLYSCALE: transform from usr to plotting coordinates.
 * At present, these do nothing (identity transformation).
 */

real plxscale(real x, real y)
{
    return x;
}

real plyscale(real x, real y)
{
    return y;
}

/*
 * PLLTYPE: select line width and dot-dash pattern.
 */

int plltype(int lwid, int lpat)
{
    if (lwid > 0) 
	c_plwidth(lwid);		/* set line width */
    if (lpat > 0)
	c_pllsty(lpat);			/* set line style */
}

/*
 * PLLINE, PLMOVE, PLPOINT: plot lines, moves, and points.
 */

int plline(real x, real y)
{
    pl_real xp,yp;

    xp=x; yp=y;
    c_pljoin(xcur,ycur,xp,yp);
    xcur=xp;
    ycur=yp;

}

int plmove(real x, real y)
{
   pl_real xp,yp;

   xcur=x; ycur=y;          /* RECALC !! */
}

int plpoint(real x, real y)
{
    int n=0,istyle=0;           /* just a dot */
    pl_real xp,yp;
    
    xp=x; yp=y;
    c_plpoin(1,&xp,&yp,1);
}

/*
 * PLCIRCLE, PLCROSS, PLBOX: plot various symbols.
 */

int plcircle(real x, real y, real r)
{
    int npnts, i;
    real theta;

    npnts = MAX(2400 * r / dxymax, 6.0);
    plmove(x + r, y);
    for (i = 1; i <= npnts; i++) {
        theta = TWO_PI * ((real) i) / ((real) npnts);
        plline(x + r * cos(theta), y + r * sin(theta));
    }
}

int plcross(real x, real y, real s)
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

int plbox(real x, real y, real s)
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
 * Imports: jus: 
 */

static pl_real fjust = 0.0;   /* pgplot default: left justified */

int pljust(int jus)          /* -1, 0, 1 for left, mid, right just */
{
    fjust = (jus < -1 ? 0.0 : (jus > 1 ? 1.0 : 0.5));
}

/*
 * PLTEXT: plot a text string.
 */

int pltext(string msg, real x, real y, real hgt, real ang)
{
    pl_real dx, dy, xp, yp;

    xp=x; yp=y;
    ang *= PI/180.0;
    dx = cos(ang);
    dy = sin(ang);
    dprintf(1,"just %f: %s\n",fjust,msg);
    c_plschr((pl_real)10.0*hgt,1.0);
    c_plptex(xp,yp,dx,dy,fjust,msg);
}

/*
 * PLFLUSH: output any pending graphics.
 */

int plflush() 
{ 

}

/*
 * PLFRAME: advance to next frame.
 */

int plframe()
{
    c_pladv(0);
}

/*
 * PLSTOP: finish up a display.
 */

#define ONEMIN (60 * 1000 * 1000)

int plstop()
{
    c_plend();
}

int pl_matrix(real *frame,int nx,int ny,
	  real xmin,real ymin,real cell,real fmin,real fmax,real findex, real blankval)
{
    real x,y,f,grayscale,ds;
    int ix,iy;
    
    printf("PL_MATRIX is unsupported in plplot version of YAPP\n");
    return(0);
}

int pl_screendump(string fname)
{
  warning("pl_screendump(%s): Not implemented for yapp_plplot",fname);
}

local void bell(void)
{
    int ring=7;
    putchar(ring);
    putchar('\n');		/* send a line feed to flush the buffer */
}

int pl_getpoly(float *x, float *y, int n)
{
    int nn, delay, loc, k, symbol;
    float xold,yold, xnew,ynew;
    char ch[10];

    bell();

    printf("Define a polygon:\n");
    printf("   LEFT   = define first/next vertex of polygon\n");
    printf("   MIDDLE = delete previous vertex point from polygon\n");
    printf("	    (this option cannot remove already plotted vertex lines\n");
    printf("   RIGHT  = close polygon\n");
#if 0
    nn = 0;                           /* count points in polygon */
    for(;;) {
        k = pgcurse_(&xnew, &ynew, ch, 1);
        if(k==0) {
            warning("Device has no cursor...");
            return;
        }
        if (ch[0]=='A') {                     /* left button = DEFINE POLYGON */
            xold = xnew;
            yold = ynew;
            dprintf (2,"Button A at x=%f y=%f\n",xold,yold);
            if (nn==0) {
                pgmove_(&xold,&yold);
                k=1;
                symbol=2;
                pgpoint_(&k,&xold,&yold,&symbol);
#if 0
                set_marker_symbol(43);      /* a 'plus' symbol */
                marker_abs_2(xold,yold);
        
#endif
            } else
                pgdraw_(&xold,&yold);
            if (nn<n) {
               x[nn] = xold;
               y[nn] = yold;
               nn++;
            } else
                warning("No more space to store this data");
        } else if (ch[0]=='D') {                  /* middle  = DELETE PREVIOUS POINT */
            if (nn>0)
                nn--;
            if (nn!=0)
               pgmove_(&x[nn-1],&y[nn-1]);    /* reset */
        } else if (ch[0]=='X') {                          /* right = QUIT */
            break;
        } else
            warning("Unknown return char %c from PGCURSE (expected A,D,X)",ch[0]);
    }   /* for(;;) */
    dprintf (2,"Button C pressed to finish up polygon (%d)\n",nn);
    if (nn>0) {
        pgdraw_(&x[0],&y[0]);       /* close polygon */
    }
    return(nn<3 ? 0 : nn);
#else
    return 0;    
#endif
    
}


/*  The rest of this file is for optional color support */
#if defined(COLOR)
/*
 * CMS_RGBSETUP: default color table initialization.
 */

#define	BLACK		0
#define RED		1
#define GREEN		2
#define BLUE		3
#define CYAN 		4	/* (Green + Blue) */	
#define MAGENTA 	5	/* (Red + Blue) */	
#define YELLOW		6	/*   (Red + Green) */	
#define ORANGE		7
#define GREEN_YELLOW	8
#define GREEN_CYAN	9
#define BLUE_CYAN	10
#define BLUE_MAGENTA	11
#define RED_MAGENTA	12
#define DARK_GRAY	13
#define LIGHT_GRAY	14
#define	WHITE   	15

local int cms_rgbsetup(void)
{
  ncolors = WHITE+1;
  /* default PGPLOT  colors: although, defined here, plpalette is not called */
  red[BLACK]=0.00;         green[BLACK]=0.00;          blue[BLACK]=0.00;
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
  red[WHITE]=1.00;         green[WHITE]=1.00;          blue[WHITE]=1.00;
}

/*
 * PLCOLOR: specify new plotting color as an integer between 0 and ncolors-1;
 * values outside this range are mapped to the nearest endpoint.
 */

void plcolor(int color)
{
    if (color < 0)
	color = 0;
    else if (color > ncolors - 1)
	color = ncolors - 1;
    c_plcol0(color);
}

/*
 * PLNCOLORS: return current value of local variable ncolors.
 */

int plncolors()
{
    return ncolors;
}

/*
 * PLPALETTE: re-initialize color table from user-supplied values.
 */

void plpalette(real *r, real *g, real *b, int nc)
{
    int i;

    if (nc > MAXCOLOR)
	error("plpalette: cannot define more than %d colors", MAXCOLOR);
    ncolors = nc;
    for (i = 0; i < ncolors; i++) {
	red[i] = r[i];
	green[i] = g[i];
	blue[i] = b[i];
#if 0
        pgscr_(&i,&red[i],&green[i],&blue[i]);
#endif
    }
    red[ncolors] = green[ncolors] = blue[ncolors] = 0.0;    /* terminate */
    c_plcol0(15);			/* reset default color to foreground */
    dprintf(0,"Setting default color to: %g %g %g\n",
            red[1], green[1], blue[1]);
    warning("plpalette not implemented yet for PLPLOT");
}

/*
 * PLLUT:  read a new RGB palette from a LUT color table
 *	   The lut must be a simple ascii file, with 1st, 2nd and 3rd column being
 *	   the R, G and B response, a real number between 0.0 and 1.0
 */
#define NOCOLOR(x)  ((x)<0||(x)>1)

void pllut(string fname, bool compress)
{
    stream cstr;
    int ncolors, nskip=0;
    real red[MAXCOLOR], green[MAXCOLOR], blue[MAXCOLOR];
    pl_real r, g, b, r_old, g_old, b_old;
    float fr, fg, fb;
    char line[256];

    if (fname==NULL || *fname == 0) return;
    cstr = stropen(fname, "r");
    ncolors = 0;
    while (fgets(line,256,cstr)) {
       if(ncolors>=MAXCOLOR) error("(%d/%d): Too many colors in %s",
              ncolors+1, MAXCOLOR, fname);
       sscanf(line,"%f %f %f",&fr, &fg, &fb);
       r=fr;  g=fg; b=fb;
       if (NOCOLOR(r) || NOCOLOR(g) || NOCOLOR(b)) {
          warning("Skipping RGB=%f %f %f",r,g,b);
	  continue;
       }
       if (ncolors>0 && r==r_old && g==g_old && b==b_old && compress) {
          nskip++;
       } else {
          dprintf(1,"LUT(%d): %f %f %f\n",ncolors,r,g,b);
          red[ncolors]   = r;   r_old = r;
          green[ncolors] = g;   g_old = g;
          blue[ncolors]  = b;   b_old = b;
          ncolors++;
       }
    }
    strclose(cstr);
    dprintf(0,"Colortable %s: %d entries; %d redundant entries skipped\n",
            fname, ncolors,nskip);
    plpalette(red, green, blue, ncolors);
}

#endif

int pl_cursor(real *x, real *y, char *c)
{
  return 0;
}

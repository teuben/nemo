/*
 * YAPP: Yet Another Plotting Package: Lupton & Monger's super mongo (SM)
 *
 *	21-May-90:  sm interface version - Peter Teuben
 *	26-oct-90:  stripped TESTBED for new NEMO
 *      12-nov-90   more debugging, almost all works, except angle and size
 *                  of fonts seems to be almost arbitrary...
 *	 7-jul-91   sm functions now have leading sm_ function names
 *	25-jul-91   using -DOLDSM can still link with older sm libraries??
 *	18-feb-93   SM2_2_1E (0th Nov 1992) testing again - worked !!!!
 *	16-jan-95   SM2.3.10 (alpha release) ; also fixed 'real'
 *      28-dec-95   prototyped interface
 */

#include <stdinc.h>

extern int    yapp_dev;	    /* see getparam.c */
extern string yapp_string;  /* see getparam.c */

local real   dxymax;	/* size of user window */
local int    iterm;	/* terminal number */

#if defined(OLDSM)
	        /* if you would need this, you have a really really old SM */
#define sm_lweight lweight
#define sm_curs curs
#define sm_draw draw
#define sm_relocate relocate
#define sm_putlabel putlabel
#define sm_location location
#define sm_graphics graphics
#define sm_hardcopy hardcopy
#define sm_alpha alpha
#define sm_erase erase
#define sm_angle angle
#define sm_ptype ptype
#define sm_ltype ltype
#define sm_expand expand
#define sm_device device
#define sm_defvar defvar
#define sm_points points
#define sm_limits limits
#endif

/* SIGH: from sm's sm_declare.h: */

extern int sm_array_to_vector( float *, int, char * );
extern int sm_axis( double, double, double, double, int, int, int, int, int);
extern void sm_box( int, int, int, int );
extern int sm_conn( float [], float [], int );
extern int sm_contour( void );
extern void sm_ctype( char * );
extern void sm_curs( float *, float *, int *); /* callable SM */
extern void sm_cursor( int );
extern int sm_defimage( float **, double, double, double, double, int, int );
extern void sm_delimage( void );
extern void sm_disable( void );
extern void sm_dot( void );
extern void sm_draw( double, double );
extern void sm_draw_point( double, double, int, int );
extern void sm_draw_surface( int, double, double, float *, int, float *,int );
extern void sm_enable( void );
extern int sm_errorbar( float [], float [], float [], int, int );
extern void sm_format( char *, char * );
extern void sm_gflush( void );
extern void sm_grelocate( double, double );
extern void sm_grid( int, int );
extern int sm_histogram( float [], float [], int );
extern void sm_label( char * );
extern void sm_levels( float [], int );
extern int sm_limits( double, double, double, double );
extern void sm_line( double, double, double, double );
extern int sm_location( int, int, int, int );
extern void sm_ltype( int );
extern void sm_lweight( double );
extern void sm_minmax( float *, float * );
extern void sm_notation( double, double, double, double );
extern int sm_peek( void );
extern void sm_plotsym( float [], float [], int, int [], int );
extern void sm_points( float [], float [], int );
extern void sm_ptype( float *, int );
extern void sm_putlabel( int, char * );
extern void sm_relocate( double, double );
extern void sm_set_viewpoint( double, double, double );
extern void sm_shade( int, float *, float *, int );
extern int sm_main( int, char ** );
extern int sm_parser( char * );
extern int sm_suspend( void );
extern void sm_ticksize( double, double, double, double );
extern void sm_window( int, int, int, int );
extern void sm_xlabel( char * );
extern void sm_ylabel( char * );
extern void sm_alpha( void );
extern void sm_angle( double );
extern void sm_defvar( char *, char * );
extern int sm_device( char * );
extern void sm_erase( void );
extern void sm_expand( double );
extern void sm_identification( char * );
extern void sm_graphics( void );
extern void sm_hardcopy( void );
extern void sm_redraw( int );
extern void sm_toscreen( double, double, int *, int * );

local int more(void);


/*
 * PLINIT: initalize the plotting package.
 */

plinit(string pltdev, real xmin, real xmax, real ymin, real ymax)
{
    float width, height, x1,x2,y1,y2;

    if (yapp_string == NULL) {
	  dprintf(0,"YAPP device is not given, use setenv YAPP or yapp=\n");
	  dprintf(0,"This is yapp_sm:");
	  dprintf(0,"all graphics output suppressed\n");
	  iterm = 0;
	  return;
    } else
	  iterm = 1;
    dprintf(1,"YAPP_SM: %s iterm=%d\n",yapp_string,iterm);

    sm_device(yapp_string);
    sm_location(1000,32000,1000,32000);
    sm_graphics();
    sm_defvar("TeX_strings","1");

    if (ymax - ymin < xmax - xmin) {
	dxymax = xmax - xmin;
	width = 1.0;
	height = (ymax - ymin) / dxymax;
    } else {
	dxymax = ymax - ymin;
	width = (xmax - xmin) / dxymax;
	height = 1.0;
    }
    
    x1=xmin; y1=ymin;
    x2=xmax; y2=ymax;
    sm_limits(x1,x2,y1,y2);                /* set user coordinates as 'cm' */

}

/*
 * PLSWAP: does nothing.
 */

plswap() { }

/*
 * PLXSCALE, PLYSCALE: transform from usr to plotting coordinates.
 * At present, these do nothing (identity transformation).
 */

real plxscale(real x, real y)
{
    if (iterm==0) return;       /* no graphics output requested */
    return x;
}

real plyscale(real x, real y)
{
    if (iterm==0) return;       /* no graphics output requested */
    return y;
}

/*
 * PLLTYPE: select line width and dot-dash pattern.
 */

plltype(int lwid, int lpat)
{
    int lp;
    float lw;

    if (iterm==0) return;       /* no graphics output requested */
    if (lwid > 0) {
	lw = lwid;
        sm_lweight(lw);
    }
    if (lpat > 0) {
	lp = lpat-1;
        sm_ltype(lp);
    }
}

/*
 * PLLINE, PLMOVE, PLPOINT: plot lines, moves, and points.
 */

plline(real x, real y)
{
    float xp,yp;

    if (iterm==0) return;       /* no graphics output requested */

    xp=x; yp=y;		/* RECALC !! */
    sm_draw(xp,yp);
}

plmove(real x, real y)
{
   float xp,yp;

   if (iterm==0) return;       /* no graphics output requested */

   xp=x; yp=y;		/* RECALC !! */
   sm_relocate(xp,yp);
}

plpoint(real x, real y)
{
    int n=0,istyle=0;		/* just a dot */
    float xp,yp,pp=0.0;
    
    if (iterm==0) return;       /* no graphics output requested */

    xp=x; yp=y;		/* RECALC !! */

/*    ptype(n,istyle);	*/
    sm_ptype(&pp,1);
    sm_points(&xp,&yp,1);
}

/*
 * PLCIRCLE, PLCROSS, PLBOX: plot various symbols.
 */

plcircle(real x, real y, real r)
{
    int npnts, i;
    real theta;

    if (iterm==0) return;       /* no graphics output requested */

    npnts = MAX(2400 * r / dxymax, 6.0);
    plmove(x + r, y);
    for (i = 1; i <= npnts; i++) {
	theta = TWO_PI * ((real) i) / ((real) npnts);
	plline(x + r * cos(theta), y + r * sin(theta));
    }
}

plcross(real x, real y, real s)
{
    if (iterm==0) return;       /* no graphics output requested */

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

plbox(real x, real y, real s)
{
    if (iterm==0) return;       /* no graphics output requested */

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
 * Imports: jus:  -1, 0, 1 for left, mid, right just 
 */

local int textjust = -1;

pljust(int jus)
{
    if (iterm==0) return;       /* no graphics output requested */

    textjust = (jus < -1 ? -1 : (jus > 1 ? 1 : jus));
}

/*
 * PLTEXT: plot a text string.
 */

pltext(string msg, real x, real y, real hgt, real ang)
{
    float dx, dy, xp, yp, ap;
    float newsize, sl, sh, slen;
    int   n, loc;

    if (iterm==0) return;       /* no graphics output requested */

    xp=x; yp=y; 
    ap=ang;
    sm_angle(ap);			/* set angle */
    newsize = 2*hgt;		/* new expansion size */
    sm_expand(newsize);
    dprintf(1,"pltext(%s) at hgt=%f ang=%f\n",msg,newsize,ap);

    sm_relocate(xp,yp);
    loc = -textjust + 5;        /* sm uses 4,5, or 6 for loc of string */
    sm_putlabel(loc,msg);

    ap=0.0; sm_angle(ap);		/* set back to default ? */
}

/*
 * PLFLUSH: output any pending graphics.
 */

plflush() 
{
/*    sm_gflush();	*/
}

/*
 * PLFRAME: advance to next frame.
 */

plframe()
{

#if 1
    more();
    sm_hardcopy();
    sm_erase();
#endif
}

/*
 * PLSTOP: finish up a display.
 */

plstop()
{
#if 0
    sm_hardcopy();
    sm_alpha();
    dprintf(0,"[Hit RETURN to continue] ");
    sm_getchar();	/* doesn't exist anymore ??? */
#else
    more();
    sm_alpha();
#endif

}

local int more()
{
    float xp, yp;
    int n;

    dprintf(0,"[YAPP_SM(pl_stop) Cursor]\n");
    sm_curs(&xp,&yp,&n);
    dprintf(1,"Cursor: %f %f %d\n",xp,yp,n);
    sm_hardcopy();
}

pl_matrix(frame,nx,ny,xmin,ymin,cell,fmin,fmax,findex,blankval)
real *frame, xmin, ymin, cell, fmin, fmax, findex, blankval;
int nx, ny;
{
    real x,y,f,grayscale,ds;
    int ix,iy;
    
    if (iterm==0) return;       /* no graphics output requested */

    warning("pl_matrix is unsupported in yapp_sm");
    return(0);
}

pl_screendump(fname)
string fname;
{
  warning("pl_screendump(%s): Not implemented for yapp_sm",fname);
}


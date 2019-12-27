/*
 * YAPP_PS.C: Yet Another Plotting Package -- PostScript (version 2).
 * Version 1	Peter Teuben	 8 July 87	IAS, Princeton.
 * Version 2	Josh Barnes	20 Jan 88	IAS, Princeton.
 *	   2.1  PJT		12 Apr 89  plcirlce()
 *	   2.2  PJT		27 sep 89  added dummy pl_screendump
 *	   2.3  PJT		26 nov 89  added timestamp back
 *	   2.4  PJT		14 mar 90  fixed GCC warnings 
 *
 *         3.0  PJT             22 nov 93  use yapp= to set ps out name
 *                                         use stropen, allow yapp=lpr
 *	   3.0a			 1 mar 94  some ansi fixes
 *	   3.1  PJT             22 jan 95  fixed ANSI 'real=float' problems
 *	   3.2  PJT	         5 may 95  put plflush() back, this is now
 *					   the official YAPP_PS 
 *                               5 dec 99  stub pl_contour
 *                              21 jan 00  stub pl_cursor
 *				 5 apr 01  added plswap/plxscale/plyscale for completeness
 *					   and plcolor (but this ought to do color too) + friends
 *                                         allow yapp.ps to be overwritten; others not
 *                              20 jun 01  debug level now 1
 *         3.3  PJT             15-jun-2018 prototypes
 */

#define VERSIONID "Version 3.3 15-jun-2018 PJT"

#include <stdinc.h>
//#include <yapp.h>

extern string yapp_string;	/* a kludge, see: getparam.c */
extern int debug_level;         /* see dprintf.c   from DEBUG env.var. */

/* Scaling and layout parameters. */

#define PTCM  (72.0 / 2.54)			/* points per cm */

#define XRANGE	(20.0 * PTCM)			/* width of user area */
#define YRANGE  (25.0 * PTCM)			/* height of user area */

#define XOFFSET 26				/* offset of user area */
#define YOFFSET 42

/* YAPP state variables. */

#define INACTIVE  0				/* not yet initialized */
#define OPENPAGE  1				/* ready to open a page */
#define OPENPATH  2				/* ready to open a path */
#define DRAWPATH  3				/* drawing a path */

local int yappstat = 0;				/* status of yapp; see above */

local stream psst;				/* postscript output stream */
local string yappfile="yapp.ps";                /* output ps filename */
local bool lpr = FALSE;                         /* flag to print & delete ? */

local int npage = 0;				/* count of pages output */

#define MAXPNT  1024				/* max. points per path */

local int npnt;					/* point count this path */

local real uscale;				/* user to pt scaling */
local real uxmin, uymin;			/* origin of user c.s. */

#define DEFWIDTH  0.5				/* default line width */
#define DEFLTYPE  1				/* default line type */
#define DEFJUST  -1				/* default text just'n */
#define YAPPFONT  "Times-Roman"			/* font used by pltext */

local real width;				/* line width in pts */
local int ltype;				/* line type: solid, etc */
local int just;					/* text centering parameter */
local real fntsz;				/* current size of font */
local string did;                               /* date id string */


#define MAXCOLOR 256
local int ncolors=2;				/* black and white ??? */
local float red[MAXCOLOR];              /* RGB color tables                 */
local float blue[MAXCOLOR];
local float green[MAXCOLOR];

static real xscale(real x);
static real yscale(real y);
static void prolog(real dx, real dy);
static void trailer(void);
static void begpage(void);
static void endpage(void);
static void begpath(void);
static void endpath(void);

extern string date_id(void);                    /* date id */


/*
 * PLINIT: initialize YAPP_PS.
 */

plinit(string opt, real x0, real x1, real y0, real y1)
{
    real dx, dy;

    if (yappstat != INACTIVE)
	error("plinit: called twice\n");
    if (yapp_string != 0 && *yapp_string != 0)
        yappfile = yapp_string;
    if (streq(yappfile,"lpr")) {
      yappfile = "yapp.tmp.ps";
      lpr = TRUE;
    } else {
      dprintf(1,"YAPP_PS %s:\n do \'lpr %s\' to print\n", VERSIONID, yappfile);
    }
    if (streq(yappfile,"yapp.ps"))
        psst = stropen(yappfile, "w!");
    else
        psst = stropen(yappfile, "w");
    if (psst == NULL)
	error("plinit: can't open %s", yappfile);
    did = date_id();
    dprintf(1,"Date_id=%s\n",did);
    dx = x1 - x0;
    dy = y1 - y0;
    if (dx / dy < XRANGE / YRANGE)
	uscale = YRANGE / dy;
    else
	uscale = XRANGE / dx;
    uxmin = x0;
    uymin = y0;
    prolog(dx, dy);
    width = DEFWIDTH;
    ltype = DEFLTYPE;
    just = DEFJUST;
    fntsz = 0.0;
}

/*
 * PLFRAME: finish current page, make ready for next.
 */

plframe()
{
    endpath();
    endpage();
    width = DEFWIDTH;
    ltype = DEFLTYPE;
    just = DEFJUST;
    fntsz = 0.0;
}

/*
 * PLFLUSH: output pending graphics
 */

plflush()
{
        /*      does nothing yet */
}  



/*
 * PLSTOP: finish up current page and quit.
 */

plstop()
{
    char cmd[256];

    endpath();
    endpage();
    trailer();
    strclose(psst);
    if (lpr) { 
      sprintf(cmd,"lpr %s; rm -f %s",yappfile, yappfile);
      dprintf(0,"%s\n",cmd);
      system(cmd);
    } 
}

/*
 * PLLTYPE: set line width and pattern.
 */

#define NPAT  5

local string dashtab[NPAT] = {
    "[] 0",			/* solid */
    "[1 3] 0",			/* dotted */
    "[3] 0",			/* short dash */
    "[6 3] 0",			/* long dash */
    "[3 3 6 3] 0",		/* short long */
};

plltype(int lwd, int ltp)
{
    real w;

    begpage();
    endpath();
    w = 0.5 * lwd;
    if (lwd != 0 && w != width) {
	fprintf(psst, "%.1f setlinewidth\n", ABS(w));
	if (w * width < 0.0)
	    fprintf(psst, "%1d setgray\n", lwd > 0 ? 0 : 1);
	width = w;
    }
    if (ltp > 0 && ltp != ltype) {
	fprintf(psst, "%s setdash\n", dashtab[(ltp - 1) % NPAT]);
	ltype = ltp;
    }
}

/*
 * PLMOVE, PLLINE, PLPOINT: low-level graphics commands.
 */

plmove(real x, real y)
{
    begpage();
    begpath();
    fprintf(psst, "%.1f %.1f mt\n", xscale(x), yscale(y));
    npnt++;
    if (npnt >= MAXPNT) {
	endpath();
	plmove(x, y);
    }
}

plline(real x, real y)
{
    begpage();
    begpath();
    fprintf(psst, "%.1f %.1f lt\n", xscale(x), yscale(y));
    npnt++;
    if (npnt >= MAXPNT) {
	endpath();
	plmove(x, y);
    }
}

plpoint(real x, real y)
{
    begpage();
    endpath();
    fprintf(psst, "%.1f %.1f pt\n", xscale(x), yscale(y));
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
    /* now set some color here ??? */
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
        dprintf(1,"->PGSCR_(%d,%g,%g,%g) \n",i,red[i],green[i],blue[i]);
	/* pgscr_(&i,&red[i],&green[i],&blue[i]); */
    }
    red[ncolors] = green[ncolors] = blue[ncolors] = 0.0;    /* terminate */
    plcolor(1);                 /* reset default color to foreground */
    dprintf(0,"Setting default color to: %g %g %g\n",
            red[1], green[1], blue[1]);
}
 

/*
 *  PLCIRCLE, PLCROSS, PLBOX: draw useful marking symbols.
 */
 
plcircle(real x, real y, real r)
{
    begpage();
    endpath();
    fprintf(psst, "%1d %.1f %.1f %.1f circ\n",
	    r > 0.0 ? 1 : 0, xscale(x), yscale(y), uscale * ABS(r));
}

plcross(real x, real y, real s)
{
    if (s > 0.0) {
	plmove(x-s, y);
	plline(x+s, y);
	plmove(x, y-s);
	plline(x, y+s);
    } else {
	s = s / 1.4142;
	plmove(x-s, y-s);
	plline(x+s, y+s);
	plmove(x-s, y+s);
	plline(x+s, y-s);
    }
}

plbox(real x, real y, real s)
{
    begpage();
    endpath();
    fprintf(psst, "%.1f %d %.1f %.1f box\n",
	    uscale * ABS(s), s > 0.0 ? 0 : 45, xscale(x), yscale(y));
}

/*
 * PLTEXT: plot text string with specified position, size and angle.
 */

#define FUDGE1  1.50		/* fudge factor to make pltext() larger by */
#define FUDGE2  0.50		/* fudge factor to offset text coords by */

pltext(string msg, real x, real y, real hgt, real ang)
{
    real fs, x1, y1;

    begpage();
    endpath();
    fs = FUDGE1 * uscale * hgt;
    if (fs != fntsz) {
	fprintf(psst, "/%s findfont %1f scalefont setfont\n", YAPPFONT, fs);
	fntsz = fs;
    }
    x1 = x + FUDGE2 * hgt * sin((PI / 180.0) * ang);
    y1 = y - FUDGE2 * hgt * cos((PI / 180.0) * ang);
    fprintf(psst, "(%s) %.1f %.1f %.1f ", msg, ang, xscale(x1), yscale(y1));
    switch (just) {
      case -1:
	fprintf(psst, "lsh\n");
	break;
      case 0:
	fprintf(psst, "csh\n");
	break;
      case 1:
	fprintf(psst, "rsh\n");
	break;
    }
}

/*
 * PLJUST: specify justification parameter.
 */

pljust(int j)
{
    just = (j < -1 ? -1 : (j > 1 ? 1 : j));
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
    return x;
}

real plyscale(real x, real y)
{
    return y;
}

/*
 * PL_MATRIX: paint a matrix (someday).
 */

pl_matrix(
    real *frame,		/* image array, stored in 1-D */
    int nx, int ny,		/* size of image */
    real xmin, real ymin,	/* lower-left corner */
    real cell,		/* pixel size in user units */
    real fmin, real fmax,	/* values mapped to white, black */
    real findex,
    real blankval)
{
    static bool virgin = TRUE;

    if (virgin) {
	warning("pl_matrix: not yet implemented");
	virgin = FALSE;
    }
}

pl_contour (real *frame, int nx, int ny, int ncntval, real *cntval)
{
    warning("pl_contour: not yet implemented for yapp_ps");
}

int pl_cursor(float *x,float *y, char *c)      
{
 return 0;
}


pl_screendump (string fname)
{
    warning("yapp_ps_new: pl_screendump not implemented");
}
/*
 * XSCALE, YSCALE: user to pt transformations.
 */

local real xscale(real x)
{
    return (uscale * (x - uxmin));
}

local real yscale(real y)
{
    return (uscale * (y - uymin));
}

/*
 * PROLOG: output prolog, make ready to begin 1st page.
 */

local void prolog(real dx, real dy)
{
    int llx, lly, urx, ury;

    llx = XOFFSET + 0.5 * (XRANGE - uscale * dx);
    lly = YOFFSET + 0.5 * (YRANGE - uscale * dy);
    urx = XOFFSET + 0.5 * (XRANGE + uscale * dx);
    ury = YOFFSET + 0.5 * (YRANGE + uscale * dy);
    fprintf(psst, "%%!PS-Adobe-1.0\n");
    fprintf(psst, "%%%%Title: yapp.ps\n");
    fprintf(psst, "%%%%Pages: (atend)\n");
    fprintf(psst, "%%%%BoundingBox: %d %d %d %d\n", llx, lly, urx, ury);
    fprintf(psst, "%%%%DocumentsFonts: %s\n", YAPPFONT);
    fprintf(psst, "%%%%EndComments\n");
    fprintf(psst, "/np {newpath} def\n");
    fprintf(psst, "/st {gsave stroke grestore} def\n");
    fprintf(psst, "/mt {moveto} def\n");
    fprintf(psst, "/lt {lineto} def\n");
    fprintf(psst, "/pt {0.33 newpath 0 360 arc closepath fill} def\n");
    fprintf(psst, "/circ {newpath 0 360 arc closepath\n");
    fprintf(psst, "       gsave setgray fill grestore stroke} def\n");
    fprintf(psst, "/boxdict 1 dict def\n");
    fprintf(psst, "/box {gsave translate rotate\n");
    fprintf(psst, "      boxdict begin /s exch def newpath\n");
    fprintf(psst, "        s s mt s s neg lt s neg s neg lt\n");
    fprintf(psst, "        s neg s lt closepath end\n");
    fprintf(psst, "      gsave 1 setgray fill grestore stroke\n");
    fprintf(psst, "      grestore} def\n");
    fprintf(psst, "/lsh {gsave translate rotate np\n");
    fprintf(psst, "      0 0 mt show grestore} def\n");
    fprintf(psst, "/csh {gsave translate rotate np\n");
    fprintf(psst, "      dup stringwidth pop 2 div neg 0 mt\n");
    fprintf(psst, "      show grestore} def\n");
    fprintf(psst, "/rsh {gsave translate rotate np\n");
    fprintf(psst, "      dup stringwidth pop neg 0 mt\n");
    fprintf(psst, "      show grestore} def\n");
    fprintf(psst, "%.1f setlinewidth\n", DEFWIDTH);
    fprintf(psst, "[] 0 setdash\n");
    fprintf(psst, "/%s findfont %1f scalefont setfont\n", YAPPFONT, 10.0);
    fprintf(psst, "(%s) 0 548 685 rsh %%%%%% date_id\n",did);   /* date id */
    fprintf(psst, "%d %d translate\n", llx, lly);
    fprintf(psst, "%%%%EndProlog\n");
    yappstat = OPENPAGE;
}

/*
 * TRAILER: output trailer.
 */

local void trailer(void)
{
    fprintf(psst, "%%%%Trailer\n");
    fprintf(psst, "%%%%Pages: %d\n", npage);
}

/*
 * BEGPAGE: start new output page.
 */

local void begpage(void)
{
    if (yappstat == OPENPAGE) {
	npage++;
	fprintf(psst, "%%%%Page: %d %d\n", npage, npage);
	fprintf(psst, "save\n");
	yappstat = OPENPATH;
    }
}

/*
 * ENDPAGE: finish output page.
 */

local void endpage(void)
{
    if (yappstat == OPENPATH) {
	fprintf(psst, "showpage\n");
	fprintf(psst, "restore\n");
	yappstat = OPENPAGE;
    }
}

/*
 * BEGPATH: start path specification.
 */

local void begpath(void)
{
    if (yappstat == OPENPATH) {
	fprintf(psst, "np\n");
	yappstat = DRAWPATH;
	npnt = 0;
    }
}

/*
 * ENDPATH: terminate path specification.
 */

local void endpath(void)
{
    if (yappstat == DRAWPATH) {
	fprintf(psst, "st\n");
	yappstat = OPENPATH;
    }
}

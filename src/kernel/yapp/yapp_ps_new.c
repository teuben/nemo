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
 */

#include <stdinc.h>

extern string yapp_string;

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

local real xscale(), yscale();			/* user to pt mappings */
local prolog(), trailer(), begpage(), endpage(), begpath(), endpath();

#define DEFWIDTH  0.5				/* default line width */
#define DEFLTYPE  1				/* default line type */
#define DEFJUST  -1				/* default text just'n */
#define YAPPFONT  "Times-Roman"			/* font used by pltext */
#define VERSIONID "Version 3.0 22-nov-92 PJT"

local real width;				/* line width in pts */
local int ltype;				/* line type: solid, etc */
local int just;					/* text centering parameter */
local real fntsz;				/* current size of font */
local string did;                               /* date id string */

/*
 * PLINIT: initialize YAPP_PS.
 */

plinit(opt, x0, x1, y0, y1)
string opt;				/* mostly ignored name string */
real x0, x1, y0, y1;
{
    real dx, dy;
    string date_id();                    /* date id */

    if (yappstat != INACTIVE)
	error("plinit: called twice\n");
    if (yapp_string != 0 && *yapp_string != 0)
        yappfile = yapp_string;
    if (streq(yappfile,"lpr")) {
      yappfile = "yapp.tmp.ps";
      lpr = TRUE;
    } else {
      dprintf(0,"YAPP_PS %s:\n do \'lpr %s\' to print\n", VERSIONID, yappfile);
    }
    psst = stropen(yappfile, "w");
    if (psst == NULL)
	error("plinit: can't open %s", yappfile);
    did = date_id();
    dprintf(0,"Date_id=%s\n",did);
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

plltype(lwd, ltp)
int lwd;
int ltp;
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

plmove(x, y)
real x, y;
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

plline(x, y)
real x, y;
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

plpoint(x, y)
real x, y;
{
    begpage();
    endpath();
    fprintf(psst, "%.1f %.1f pt\n", xscale(x), yscale(y));
}

/*
 *  PLCIRCLE, PLCROSS, PLBOX: draw useful marking symbols.
 */
 
plcircle(x, y, r)
real x, y, r;
{
    begpage();
    endpath();
    fprintf(psst, "%1d %.1f %.1f %.1f circ\n",
	    r > 0.0 ? 1 : 0, xscale(x), yscale(y), uscale * ABS(r));
}

plcross(x, y, s)
real x, y, s;
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

plbox(x, y, s)
real x, y, s;
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

pltext(msg, x, y, hgt, ang)
string msg;
real x, y;
real hgt;
real ang;
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

pljust(j)
int j;
{
    just = (j < -1 ? -1 : (j > 1 ? 1 : j));
}

/*
 * PL_MATRIX: paint a matrix (someday).
 */

pl_matrix(frame, nx, ny, xmin, ymin, cell, fmin, fmax, findex, blank)
real *frame;		/* image array, stored in 1-D */
int nx, ny;		/* size of image */
real xmin, ymin;	/* lower-left corner */
real cell;		/* pixel size in user units */
real fmin, fmax;	/* values mapped to white, black */
real findex;
real blank;
{
    static bool virgin = TRUE;

    if (virgin) {
	printf("pl_matrix: not yet implemented\n");
	virgin = FALSE;
    }
}

pl_screendump (fname)
string fname;
{
    printf("yapp_ps_new: pl_screendump not implemented\n");
}
/*
 * XSCALE, YSCALE: user to pt transformations.
 */

local real xscale(x)
real x;
{
    return (uscale * (x - uxmin));
}

local real yscale(y)
real y;
{
    return (uscale * (y - uymin));
}

/*
 * PROLOG: output prolog, make ready to begin 1st page.
 */

local prolog(dx, dy)
real dx, dy;
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

local trailer()
{
    fprintf(psst, "%%%%Trailer\n");
    fprintf(psst, "%%%%Pages: %d\n", npage);
}

/*
 * BEGPAGE: start new output page.
 */

local begpage()
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

local endpage()
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

local begpath()
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

local endpath()
{
    if (yappstat == DRAWPATH) {
	fprintf(psst, "st\n");
	yappstat = OPENPATH;
    }
}

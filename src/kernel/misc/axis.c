/*
 * AXIS.C: routines to draw various kinds of axies.
 * 
 *	xx-xxx-86   created...					 JEB
 *	10-sep-90   trimval() declared as void in both places    PJT
 *	19-oct-90   index() -> strchr()                          PJT
 *	15-apr-95   ansi prototypes, no more ARGS		 PJT
 *      17-feb-97   proc prototype (for SINGLEPREC)		 pjt
 *       2-jun-13   3 digits,finally?                            pjt
 */

#include <stdinc.h>
#include <yapp.h>
#include <axis.h>

struct axisvars xaxisvar = { 3, 0.2, 0.0, 0.4, 0.32, 1.0, 0.32 };
struct axisvars yaxisvar = { 3, 0.2, 0.0, 0.3, 0.32, 1.0, 0.32 };

/*
 * Accessors for old-fashioned routines.
 */

#define nxdig  xaxisvar.ndig
#define xtikup xaxisvar.tikup
#define xtikdn xaxisvar.tikdn
#define xnumdn xaxisvar.numdn
#define xsznum xaxisvar.sznum
#define xlabdn xaxisvar.labdn
#define xszlab xaxisvar.szlab

#define nydig  yaxisvar.ndig
#define ytikrt yaxisvar.tikup
#define ytiklf yaxisvar.tikdn
#define ynumlf yaxisvar.numdn
#define ysznum yaxisvar.sznum
#define ylablf yaxisvar.labdn
#define yszlab yaxisvar.szlab

local void trimval(string);

/*
 * FORMALAXIS: if TRUE, draw axies in formal style.
 */

bool formalaxis;

/* #define FUDGEAXIS  TRUE */

/*
 * XAXIS: draw an x axis.
 */
 
void xaxis(
    real x0, real y0,	    /* axis starting point */
    real xl,		    /* axis length */
    real tick[],	    /* tick values or limits */
    int nticks,		    /* number of ticks */
    axis_proc xtrans,	    /* plotting transformation */
    string label	    /* label for axis */
)
{
    int i;
    real t, x;
    char val[32];

    plmove(x0, y0);
    plline(x0 + xl, y0);
    for (i = 0; i < abs(nticks); i++) {
	if (nticks > 0)
	    t = tick[i];
	else
	    t = tick[0] + (i + 1) * (tick[1] - tick[0]) / (1.0 - nticks);
	x = (*xtrans)(t);
	if (label != NULL) {
	    plmove(x, y0 + xtikup);
	    plline(x, y0 - xtikdn);
	    sprintf(val, "%-.*f", nxdig, t);
	    trimval(val);
	    pljust(0);
	    pltext(val, x, y0 - xnumdn, xsznum, 0.0);
	} else {
	    plmove(x, y0 - xtikup);
	    plline(x, y0 + xtikdn);
	}
    }
    if (label != NULL && *label != 0) {
	if (! formalaxis) {
	    pljust(0);
	    pltext(label, x0 + xl / 2, y0 - xlabdn, xszlab, 0.0);
	} else {
	    pljust(1);
#ifdef FUDGEAXIS
	    pltext(label, x0 + xl, y0 - xlabdn + xszlab/4.0, xszlab, 0.0);
#else
	    pltext(label, x0 + xl, y0 - xlabdn, xszlab, 0.0);
#endif
	}
    }
    pljust(-1);
}

/*
 * YAXIS: draw a y axis.
 */

void yaxis(
    real x0,
    real y0,    		/* axis starting point */
    real yl,		        /* axis length */
    real tick[],		/* tick values or limits */
    int nticks,		        /* number of ticks */
    axis_proc ytrans,		/* plotting transformation */
    string label		/* label for axis */
)
{
    int i;
    real t, y;
    char val[32];

    plmove(x0, y0);
    plline(x0, y0 + yl);
    for (i = 0; i < abs(nticks); i++) {
	if (nticks > 0)
	    t = tick[i];
	else
	    t = tick[0] + (i + 1) * (tick[1] - tick[0]) / (1.0 - nticks);
	y = (*ytrans)(t);
	if (label != NULL) {
	    plmove(x0 + ytikrt, y);
	    plline(x0 - ytiklf, y);
	    sprintf(val, "%-.*f", nydig, t);
	    trimval(val);
	    if (! formalaxis) {
		pljust(0);
		pltext(val, x0 - ynumlf, y, ysznum, 90.0);
	    } else {
		pljust(1);
		pltext(val, x0 - ynumlf, y, ysznum, 0.0);
	    }
	} else {
	    plmove(x0 - ytikrt, y);
	    plline(x0 + ytiklf, y);
	}
    }
    if (label != NULL && *label != 0) {
	if (! formalaxis) {
	    pljust(0);
	    pltext(label, x0 - ylablf, y0 + yl / 2, yszlab, 90.0);
	} else {
	    pljust(1);
#ifdef FUDGEAXIS
	    pltext(label, x0 - ylablf, y0 + yl - yszlab/4.0, yszlab, 0.0);
#else
	    pltext(label, x0 - ylablf, y0 + yl - yszlab/2.0, yszlab, 0.0);
#endif
	}
    }
    pljust(-1);
}

local void trimval(string val)
{
    char *ep;

    if (strchr(val, '.') != NULL) {
	ep = val + strlen(val);
	while (*--ep == '0')
	    *ep = 0;
	if (*ep == '.')
	    *ep = 0;
    }
}

/*
 * AXVAR: specify general parameters for axis subroutines.
 */

void axvar(real sn, real sl, real scale)
{
    if (sn > 0.0)
	xsznum = ysznum = sn;
    if (sl > 0.0)
	xszlab = yszlab = sl;
    if (scale > 0.0) {
	xtikup *= scale;
	xtikdn *= scale;
	xnumdn *= scale;
	xsznum *= scale;
	xlabdn *= scale;
	xszlab *= scale;
	ytikrt *= scale;
	ytiklf *= scale;
	ynumlf *= scale;
	ysznum *= scale;
	ylablf *= scale;
	yszlab *= scale;
    }
}

/*
 * XAXVAR: specify parameters for xaxis.
 */

void xaxvar(int nxd, real xtu, real xtd, real xnd, real xld)
{
    if (nxd >= 0)
	nxdig = nxd;
    if (xtu >= 0.0)
	xtikup = xtu;
    if (xtd >= 0.0)
	xtikdn = xtd;
    if (xnd >= 0.0) 
	xnumdn = xnd;
    if (xld >= 0.0)
	xlabdn = xld;
}

/*
 * YAXVAR: specify parameters for yaxis.
 * Imports: as for xaxvar.
 */

void yaxvar(int nyd, real ytr, real ytl, real ynl, real yll)
{
    if (nyd >= 0)
	nydig = nyd;
    if (ytr >= 0.0)
	ytikrt = ytr;
    if (ytl >= 0.0)
	ytiklf = ytl;
    if (ynl >= 0.0)
	ynumlf = ynl;
    if (yll >= 0.0)
	ylablf = yll;
}

#ifdef TESTBED

#include <getparam.h>

string defv[] = {
    "formalaxis=false",
    "xlabdn=1.0",
    "xszlab=0.32",
    "ylablf=1.0",
    "yszlab=0.32",
    NULL,
};

real xrange[] = { -2.0, 2.0 }, yrange[] = { 0.0, 4.0 };

nemo_main()
{
    real xtrans(real), ytrans(real);

    formalaxis = getbparam("formalaxis");
    xlabdn = getdparam("xlabdn");
    xszlab = getdparam("xszlab");
    ylablf = getdparam("ylablf");
    yszlab = getdparam("yszlab");
    plinit("***", 0.0, 20.0, 0.0, 20.0);
    xaxis( 2.0,  2.0, 16.0, xrange, -3, xtrans, "x");
    xaxis( 2.0, 18.0, 16.0, xrange, -3, xtrans, NULL);
    yaxis( 2.0,  2.0, 16.0, yrange, -7, ytrans, "y");
    yaxis(18.0,  2.0, 16.0, yrange, -7, ytrans, NULL);
    plstop();
}

real xtrans(real x)
{
    return (10.0 + 4.0 * x);
}

real ytrans(real y)
{
    return (2.0 + 4.0 * y);
}

#endif

/*
 * SNAPPLOT3D.C: plot particle positions from a snapshot output file;

 *      V0.1  14-feb-11  Created, cloned off snapplot3, using s2plot library
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <filefn.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <loadobj.h>
#include <yapp.h>
#include <axis.h>
#include "s2plot.h"

#ifdef HAVE_LIBPGPLOT
#define COLOR
#endif

string defv[] = {
    "in=???\n                     input file name",
    "times=all\n		  range of times to plot",

    "xvar=x\n			  x-axis plotting variable",
    "xlabel=\n			  x-axis label; defaults to xvar",
    "xrange=-2.0:2.0\n		  x-axis variable range",

    "yvar=y\n			  y-axis plotting variable",
    "ylabel=\n			  y-axis label; defaults to yvar",
    "yrange=-2.0:2.0\n		  y-axis variable range",

    "zvar=z\n			  y-axis plotting variable",
    "zlabel=\n			  y-axis label; defaults to yvar",
    "zrange=-2.0:2.0\n		  y-axis variable range",

    "visib=1\n			  determine visibility of points",
    "psize=0\n			  point type and size",
    "fill_circle=t\n		  fill points plotted as circles",
    "formal=false\n		  produce more publication-style plots",
    "nobox=false\n		  draw axis, ticks, labels",

    "xbox=2.0:10.0\n		  extent of x-y frame in x direction",
    "ybox=2.0:10.0\n		  extent of x-y frame in y direction",
    "zbox=11.0:19.0\n		  extent of zy and xz frame in x resp. y direction",

    "nxticks=7\n		  number of ticks on x axis",
    "nyticks=7\n		  number of ticks on y axis",
    "nzticks=7\n		  number of ticks on y axis",

    "xticks=\n                    x-tickmark positions, if not default",
    "yticks=\n		          y-",
    "zticks=\n		          z-",

#ifdef COLOR
    "color=0\n			  determine color of points",
    "color_table=\n		  specify new color table to use",
    "crange=0:1\n                 range in colors to map",
#endif
    "VERSION=0.2\n		  2-nov-05 PJT",
    NULL,
};

string usage = "plot particle positions in a 3d s2plot window";

string cvsid="$Id$";




#define MAXTICKS  32

local bool trakflag;			/* if TRUE, plot trajectories */
local string input, times;
local stream instr;
local string xvar, yvar, zvar;
local btrproc xfunc, yfunc, zfunc;
local string xlabel, ylabel, zlabel;
local string visib, psize, color;
local btiproc vfunc;
local btrproc pfunc, cfunc;
local bool fillcircle;
local bool formal;
local bool nobox;
local real xbox[3], ybox[3], zbox[3];
local real xrange[3], yrange[3], zrange[3], crange[3];
local int nxticks, nyticks, nzticks;
local real xticks[MAXTICKS], yticks[MAXTICKS], zticks[MAXTICKS];

/*
 * Data read from input file.
 */

local real *timeptr = NULL;
local int nbody = 0;
local real *massptr = NULL;
local real *phaseptr = NULL;
local real *phiptr = NULL;
local real *accptr = NULL;
local real *auxptr = NULL;

/* Layout:
 *   BOX3  N/A
 *   BOX1  BOX2
 *
 * Variables labeled 1,2,3 belong to box 1,2,3
 *                   x,y,z belong to the plotted x,y,z variable
 */



real xtrans1(real), ytrans1(real), xtrans2(real), ytrans2(real), xtrans3(real), ytrans3(real);

local bool scansnap(void);

#ifndef FRAMEDELAY
#  define FRAMEDELAY 1
#endif

nemo_main()
{
    permanent bool first=TRUE;
    int argc = 1;
    static char *argv[] = {"snaps2plot", NULL};
    float size = 1.0;

    s2opend("/?", argc, argv);
    s2swin(-size,size, -size,size, -size,size);
    s2show(1);


    setparams();
    instr = stropen(input, "r");
    get_history(instr);
    compfuncs();
    plinit("", 0.0, 20.0, 0.0, 20.0);
#ifdef COLOR
    setcolors();
#endif
    while (scansnap()) {
	if (! trakflag) {
	  if (first) {
	    first = FALSE;
	  } else {
	    sleep(FRAMEDELAY);
	    plframe();
	  } 
	}
	plotbox();
	plltype(1, 0);
	plotsnap();
    }
    plstop();
}

setparams()
{
    trakflag = FALSE;
    input = getparam("in");
    times = getparam("times");

    xvar = getparam("xvar");
    if (hasvalue("xlabel"))
        xlabel = getparam("xlabel");
    else
	xlabel = xvar;
    setrange(xrange, getparam("xrange"));

    yvar = getparam("yvar");
    if (hasvalue("ylabel"))
        ylabel = getparam("ylabel");
    else
	ylabel = yvar;
    setrange(yrange, getparam("yrange"));
 
    zvar = getparam("zvar");
    if (hasvalue("zlabel"))
        zlabel = getparam("zlabel");
    else
	zlabel = zvar;
    setrange(zrange, getparam("zrange"));
 
    visib = getparam("visib");
    psize = getparam("psize");
    fillcircle = getbparam("fill_circle");
    formal = getbparam("formal");
    nobox = getbparam("nobox");
    setrange(xbox, getparam("xbox"));
    setrange(ybox, getparam("ybox"));
    setrange(zbox, getparam("zbox"));

    if (hasvalue("xticks"))
      setticks(xticks, &nxticks, getparam("xticks"));
    else {
      xticks[0] = xrange[0];
      xticks[1] = xrange[1];
      nxticks = - getiparam("nxticks");
    } 
    if (hasvalue("yticks"))
      setticks(yticks, &nyticks, getparam("yticks"));
    else {
      yticks[0] = yrange[0];
      yticks[1] = yrange[1];
      nyticks = - getiparam("nyticks");
    } 
    if (hasvalue("zticks"))
      setticks(zticks, &nzticks, getparam("zticks"));
    else {
      zticks[0] = zrange[0];
      zticks[1] = zrange[1];
      nzticks = - getiparam("nzticks");
    } 
    
#ifdef COLOR
    color = getparam("color");
    setrange(crange, getparam("crange"));
#endif
}

setrange(real *rval, string rexp)
{
    char *cptr, *tmpstr;
    double dpar;

    cptr = strchr(rexp, ':');
    if (cptr != NULL) {
        tmpstr = allocate(cptr-rexp+1);
        strncpy(tmpstr,rexp,cptr-rexp);
        if (nemoinpd(tmpstr,&dpar,1) != 1)
            error("setrange: parsing error %s",tmpstr);
        free(tmpstr);
        rval[0] = dpar;

        if (nemoinpd(cptr+1,&dpar,1) != 1)
            error("setrange: parsing error %s",cptr+1);
	rval[1] = dpar;
    } else {
        rval[0] = 0.0;
        if (nemoinpd(rexp,&dpar,1) != 1)
            error("setrange: parsing error %s",rexp);
	rval[1] = dpar;
    }
    rval[2] = rval[1] - rval[0];
}


extern btrproc btrtrans(string);	/* in reality: rproc */
extern btiproc btitrans(string);	/* in reality: iproc */

compfuncs()
{
    xfunc = btrtrans(xvar);
    yfunc = btrtrans(yvar);
    zfunc = btrtrans(zvar);

    vfunc = btitrans(visib);
    pfunc = btrtrans(psize);
#ifdef COLOR
    cfunc = btrtrans(color);
#endif
}

#ifdef COLOR

#define MAXCOL  256

setcolors()
{
    stream cstr;
    int ncolors;
    real red[MAXCOL], green[MAXCOL], blue[MAXCOL];
    char line[256];

    if (hasvalue("color_table")) {
	cstr = stropen(getparam("color_table"), "r");
        if(qsf(cstr)) {
	    get_data(cstr, "ncolors", IntType, &ncolors, 0);
	    get_data(cstr, "red", RealType, red, ncolors, 0);
	    get_data(cstr, "green", RealType, green, ncolors, 0);
	    get_data(cstr, "blue", RealType, blue, ncolors, 0);
	    plpalette(red, green, blue, ncolors);
        } else {
            rewind(cstr);
            ncolors = 0;
            while (fgets(line,256,cstr)) {
                sscanf(line,"%g %g %g",&red[ncolors], &green[ncolors], &blue[ncolors]);
                ncolors++;
            }
	    plpalette(red, green, blue, ncolors);
        }
	strclose(cstr);
        dprintf(0,"New colortable with %d entries\n",ncolors);
    }
}

#endif

#ifndef TIMEFUZZ
#  define TIMEFUZZ  0.001
#endif

bool scansnap(void)
{
    bool success;
    int i;
    real *pptr, *xptr;
    

    success = FALSE;
    while (! success) {
	get_history(instr);
	if (! get_tag_ok(instr, SnapShotTag))
	    return (FALSE);
	get_set(instr, SnapShotTag);
	get_set(instr, ParametersTag);
	if (get_tag_ok(instr, TimeTag)) {
	    if (timeptr == NULL)
		timeptr = (real *) allocate(sizeof(real));
	    get_data_coerced(instr, TimeTag, RealType, timeptr, 0);
	}
	if (get_tag_ok(instr, NobjTag))
	    get_data(instr, NobjTag, IntType, &nbody, 0);
	get_tes(instr, ParametersTag);
	if (get_tag_ok(instr, ParticlesTag) &&
	      (timeptr == NULL || streq(times, "all") ||
	         within(*timeptr, times, TIMEFUZZ))) {
	    get_set(instr, ParticlesTag);
	    if (get_tag_ok(instr, MassTag)) {
		if (massptr == NULL)
		    massptr = (real *) allocate(sizeof(real) * nbody);
		get_data_coerced(instr, MassTag, RealType, massptr,
				 nbody, 0);
	    }
	    if (get_tag_ok(instr, PhaseSpaceTag)) {
		if (phaseptr == NULL)
		    phaseptr = (real *) allocate(sizeof(real) * 2*NDIM * nbody);
		get_data_coerced(instr, PhaseSpaceTag, RealType, phaseptr,
				 nbody, 2, NDIM, 0);
		success = TRUE;
	    } else if (get_tag_ok(instr, PosTag)) {
	        real *ptmp = (real *) allocate(sizeof(real)*nbody*NDIM);
		if (phaseptr == NULL)
		    phaseptr = (real *) allocate(sizeof(real) * 2*NDIM * nbody);
		get_data_coerced(instr, PosTag, RealType, ptmp, nbody, NDIM, 0);

		for (i=0, pptr=phaseptr, xptr=ptmp; i<nbody; i++) {
		  *pptr++ = *xptr++;
		  *pptr++ = *xptr++;
		  if (NDIM==3) *pptr++ = *xptr++;
		  pptr += NDIM;
		}
		get_data_coerced(instr, VelTag, RealType, ptmp, nbody, NDIM, 0);
		for (i=0, pptr=phaseptr+NDIM, xptr=ptmp; i<nbody; i++) {
		  *pptr++ = *xptr++;
		  *pptr++ = *xptr++;
		  if (NDIM==3) *pptr++ = *xptr++;
		  pptr += NDIM;
		}
		free(ptmp);
		success = TRUE;
	    }
	    if (get_tag_ok(instr, PotentialTag)) {
		if (phiptr == NULL)
		    phiptr = (real *) allocate(sizeof(real) * nbody);
		get_data_coerced(instr, PotentialTag, RealType, phiptr,
				 nbody, 0);
	    }
	    if (get_tag_ok(instr, AccelerationTag)) {
		if (accptr == NULL)
		    accptr = (real *) allocate(sizeof(real) * NDIM * nbody);
		get_data_coerced(instr, AccelerationTag, RealType, accptr,
				 nbody, NDIM, 0);
	    }
	    if (get_tag_ok(instr, AuxTag)) {
		if (auxptr == NULL)
		    auxptr = (real *) allocate(sizeof(real) * nbody);
		get_data_coerced(instr, AuxTag, RealType, auxptr,
				 nbody, 0);
	    }
	    get_tes(instr, ParticlesTag);
	}
	get_tes(instr, SnapShotTag);
    }
    return TRUE;
}

plotbox()
{
    char msg[128];

    if (formal) {
	formalaxis = TRUE;
	xaxisvar.labdn = 0.44;
	xaxisvar.szlab = 0.40;
	yaxisvar.numdn = 0.24;
	yaxisvar.labdn = 0.24;
	yaxisvar.szlab = 0.40;
    }
    if (! nobox) {
        xaxis(xbox[0], ybox[0], xbox[2], xticks, nxticks, xtrans1, xlabel);
	xaxis(xbox[0], ybox[1], xbox[2], xticks, nxticks, xtrans1, NULL);
        yaxis(xbox[0], ybox[0], ybox[2], yticks, nyticks, ytrans1, ylabel);
	yaxis(xbox[1], ybox[0], ybox[2], yticks, nyticks, ytrans1, NULL);

        xaxis(zbox[0], ybox[0], zbox[2], zticks, nzticks, xtrans2, zlabel);
	xaxis(zbox[0], ybox[1], zbox[2], zticks, nzticks, xtrans2, NULL);
        yaxis(zbox[0], ybox[0], ybox[2], yticks, nyticks, ytrans2, "");
	yaxis(zbox[1], ybox[0], ybox[2], yticks, nyticks, ytrans2, NULL);

        xaxis(xbox[0], zbox[0], xbox[2], xticks, nxticks, xtrans3, "");
	xaxis(xbox[0], zbox[1], xbox[2], xticks, nxticks, xtrans3, NULL);
        yaxis(xbox[0], zbox[0], ybox[2], zticks, nzticks, ytrans3, zlabel);
	yaxis(xbox[1], zbox[0], ybox[2], zticks, nzticks, ytrans3, NULL);

	if (! formal) {
	    sprintf(msg, "File: %s", input);
	    pltext(msg, xbox[0], zbox[1] + 0.4, 0.32, 0.0);
	    if (timeptr != NULL) {
	      sprintf(msg, "Time: %8.3f", *timeptr);
	      pltext(msg, xbox[0] + 10.0, zbox[1] + 0.4, 0.32, 0.0);
	    }
	}


    } else {
	if (! formal && timeptr != NULL) {
	    sprintf(msg, "%.2f", *timeptr);
	    pltext(msg, xbox[1] - 1.0, zbox[1] - 1.0, 0.24, 0.0);
	}
    }
}

setticks(real *tiks, int *ntik, string tikstr)
{
    *ntik = nemoinpr(tikstr, tiks, MAXTICKS);
    if (*ntik < 0) error("(5d) Parsing %s",*ntik, tikstr);
}

plotsnap()
{
    real t, *mp, *psp, *pp, *ap, *acp;
    int vismax, visnow, i, vis, icol;
    real psz, col, x, y, z;
    real x1,x2,x3,y1,y2,y3;
    bool v1,v2,v3;
    Body b;

    t = (timeptr != NULL ? *timeptr : 0.0);	/* get current time value   */
    CLRV(Acc(&b));				/* zero unsupported fields  */
    Key(&b) = 0;
    visnow = vismax = 0;
    do {					/* loop painting layers     */
	visnow++;				/*   make next layer visib. */
	mp  = massptr;				/*   (re)set data pointers  */
	psp = phaseptr;
	pp  = phiptr;
	ap  = auxptr;
	acp = accptr;
	for (i = 0; i < nbody; i++) {		/*   loop over all bodies   */
	    Mass(&b) = (mp != NULL ? *mp++ : 0.0);
						/*     set mass if supplied */
	    SETV(Pos(&b), psp);			/*     always set position  */
	    psp += NDIM;			/*     and advance p.s. ptr */
	    SETV(Vel(&b), psp);			/*     always set velocity  */
	    psp += NDIM;			/*     and advance ptr      */
	    Phi(&b) = (pp != NULL ? *pp++ : 0.0);	
	    Aux(&b) = (ap != NULL ? *ap++ : 0.0);
	    if (acp) {				
	    	SETV(Acc(&b),acp);		/*     set accel's          */	
	    	acp += NDIM;			/*     and advance ptr      */
	    }
	    					/*     set phi,aux if given */
	    vis = (*vfunc)(&b, t, i);		/*     evaluate visibility  */
	    vismax = MAX(vismax, vis);		/*     remember how hi to go*/
	    if (vis == visnow) {		/*     if body is visible   */
	        x = (*xfunc)(&b, t, i);	        /*     evaluate x,y,z coords*/
		y = (*yfunc)(&b, t, i);
		z = (*zfunc)(&b, t, i);

		x1 = xtrans1(x);
		y1 = ytrans1(y);

		x2 = xtrans2(z);
		y2 = ytrans2(y);

		x3 = xtrans3(x);
		y3 = ytrans3(z);
		
		v1 = (xbox[0] < x1 && x1 < xbox[1] && ybox[0] < y1 && y1 < ybox[1]);
		v2 = (zbox[0] < x2 && x2 < zbox[1] && ybox[0] < y2 && y2 < ybox[1]);
		v3 = (xbox[0] < x3 && x3 < xbox[1] && zbox[0] < y3 && y3 < zbox[1]);

		if (v1 || v2 || v3) {
		    psz = (*pfunc)(&b, t, i);
#ifdef COLOR
		    col = (*cfunc)(&b, t, i);
                    col = (col - crange[0])/(crange[1] - crange[0]);
		    icol = 1 + (plncolors() - 2) *
			         MAX(0.0, MIN(1.0, col));
		    plcolor(icol);
#endif
		    if (v1) {
		      if (psz == 0.0)
			plpoint(x1, y1);
		      else if (psz < 0.0) {
			if (fillcircle)
			    plcross(x1, y1, - psz);
                        else
                            plcross(x1, y1, psz);
		      } else if (psz > 0.0) {
			if (fillcircle)
			    plcircle(x1, y1, -psz);
			else
			    plcircle(x1, y1, psz);
		      }
		    }
		    if (v2) {
		      if (psz == 0.0)
			plpoint(x2, y2);
		      else if (psz < 0.0) {
			if (fillcircle)
			    plcross(x2, y2, - psz);
                        else
                            plcross(x2, y2, psz);
		      } else if (psz > 0.0) {
			if (fillcircle)
			    plcircle(x2, y2, -psz);
			else
			    plcircle(x2, y2, psz);
		      }
		    }
		    if (v3) {
		      if (psz == 0.0)
			plpoint(x3, y3);
		      else if (psz < 0.0) {
			if (fillcircle)
			    plcross(x3, y3, - psz);
                        else
                            plcross(x3, y3, psz);
		      } else if (psz > 0.0) {
			if (fillcircle)
			    plcircle(x3, y3, -psz);
			else
			    plcircle(x3, y3, psz);
		      }
		    }

		}
	    }
	}
    } while (visnow < vismax);			/* until final layer done   */
#ifdef COLOR
    plcolor(32767);				/* reset to white */
#endif
}



/* lower left box */
real xtrans1( real x )
{
    return xbox[0] + xbox[2] * (x - xrange[0]) / xrange[2];
}

real ytrans1( real y )
{
    return ybox[0] + ybox[2] * (y - yrange[0]) / yrange[2];
}

/* lower right box */
real xtrans2( real x )
{
    return zbox[0] + zbox[2] * (x - zrange[0]) / zrange[2];
}

real ytrans2( real y )
{
    return ybox[0] + ybox[2] * (y - yrange[0]) / yrange[2];
}

/* uppper left box */
real xtrans3( real x )
{
    return xbox[0] + xbox[2] * (x - xrange[0]) / xrange[2];
}

real ytrans3( real y )
{
    return zbox[0] + zbox[2] * (y - zrange[0]) / zrange[2];
}


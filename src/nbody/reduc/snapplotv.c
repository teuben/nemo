/*
 * SNAPPLOTV.C: plot particle positions from a snapshot output file
 *              including vectors for velocity
 *      V1.0  cloned off snapplot V3.1                                    pjt
 *       1.1  9-oct-03  support new snapshot option pos/vel               pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <loadobj.h>
#include <yapp.h>
#include <axis.h>

#define VECTOR  1

string defv[] = {
    "in=???\n                     input file name",
    "times=all\n		  range of times to plot",
    "xvar=x\n			  x-axis plotting variable",
    "xlabel=\n			  x-axis label; defaults to xvar",
    "xrange=-2.0:2.0\n		  x-axis variable range",
    "yvar=y\n			  y-axis plotting variable",
    "ylabel=\n			  y-axis label; defaults to yvar",
    "yrange=-2.0:2.0\n		  y-axis variable range",
#ifdef VECTOR    
    "vxvar=vx\n                   vx-variable",
    "vyvar=vy\n                   vy-variable",
    "scale=1\n                    mult. factor for vx/vy to make 1 cm vectors",
#endif    
    "visib=1\n			  determine visibility of points",
    "psize=0\n			  point type and size",
    "fill_circle=t\n		  fill points plotted as circles",
    "formal=false\n		  produce publication-style plots",
    "xbox=2.0:18.0\n		  extent of frame in x direction",
    "ybox=2.0:18.0\n		  extent of frame in y direction",
    "nobox=false\n		  draw axis, ticks, labels",
    "nxticks=7\n		  number of ticks on x axis",
    "nyticks=7\n		  number of ticks on y axis",
    "xticks=\n                    x-tickmark positions, if not default",
    "yticks=\n		          y-",
    "nxy=1,1\n                    Number of X and Y panels per page",
#ifdef COLOR
    "color=0\n			  determine color of points",
    "color_table=\n		  specify new color table to use",
#endif
    "frame=\n			  base filename for rasterfiles(5)",
    "VERSION=1.1\n		  9-oct-03 PJT",
    NULL,
};

string usage = "plot particle positions from a snapshot file";



#define MAXTICKS  32

local bool trakflag;			/* if TRUE, plot trajectories */
local string input, times;
local stream instr;
local string xvar, yvar, vxvar, vyvar;
local btrproc xfunc, yfunc, vxfunc, vyfunc;
local string xlabel, ylabel;
local string visib, psize, color;
local btiproc vfunc;
local btrproc pfunc, cfunc;
local string frame;
local bool fillcircle;
local bool formal;
local bool nobox;
local real xbox[3], ybox[3];
local real xrange[3], yrange[3];
local real scale;
local int  nxy[2], ix, iy;

/*
 * Data read from input file.
 */

local real *timeptr = NULL;
local int nbody = 0;
local real *massptr = NULL;
local real *phaseptr = NULL;
local real *phiptr = NULL;
local real *auxptr = NULL;

real xtrans(real),  ytrans(real);
void arrow(real, real, real, real);

#ifndef FRAMEDELAY
#  define FRAMEDELAY 1
#endif

nemo_main()
{
    permanent bool first=TRUE;
    int frameno;
    bool scansnap();

    setparams();
    instr = stropen(input, "r");
    get_history(instr);
    compfuncs();
    plinit("", 0.0, 20.0, 0.0, 20.0);
#ifdef COLOR
    setcolors();
#endif
    frameno = 0;	/* BUG for nx,ny>1: each window counts as a frame */
    ix = 0;		/* top left corner */
    iy = 0;		/* if nx>1 and/or ny> 1 */
    while (scansnap()) {
	if (! trakflag) {
            if (ix==0 && iy==0) {
                if (first) {
                    first = FALSE;
                } else {
    	            sleep(FRAMEDELAY);
	            plframe();
                } 
	    }
	}
	if (frameno == 0 || ! trakflag)
	    plotbox();
	plltype(1, 0);
	plotsnap();
#ifdef HACKTRACK
	if (frameno == 0)
	    plottrack();
#endif
	if (frame)
	    scrdump(frameno);
	frameno++;

	/* the remainder is to keep track of the sub-window on the page */

	ix++;
	xbox[0] += xbox[2];
	xbox[1] += xbox[2];
	if (ix==nxy[0]) {		/* select a new row: go down */
	    ix=0;
	    xbox[0] -= nxy[0]*xbox[2];
	    xbox[1] -= nxy[0]*xbox[2];
	    iy++;
	    ybox[0] -= ybox[2];
	    ybox[1] -= ybox[2];
	    if (iy==nxy[1]) {		/* full page, need to flush graphics */
	        iy=0;
                ybox[0] += nxy[1]*ybox[2];
                ybox[1] += nxy[1]*ybox[2];
	    }
	}
    }
    plstop();
}

setparams()
{
    string tail();

    trakflag = (strncmp(tail(getargv0()), "trak", 4) == 0);
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
    vxvar = getparam("vxvar");
    vyvar = getparam("vyvar");
    scale = getdparam("scale");
    setrange(yrange, getparam("yrange"));
    visib = getparam("visib");
    psize = getparam("psize");
    fillcircle = getbparam("fill_circle");
    formal = getbparam("formal");
    nobox = getbparam("nobox");
    setrange(xbox, getparam("xbox"));
    setrange(ybox, getparam("ybox"));
    switch (nemoinpi(getparam("nxy"),nxy,2)) {
        case 0: nxy[0] = nxy[1] = 1;            break;
        case 1: nxy[1] = nxy[0];                break;
        case 2:                                 break;
        default: error("Parsing error nxy=; max two integer values"); break;
    }
    xbox[2] /= nxy[0];
    ybox[2] /= nxy[1];
    xbox[1] = xbox[0] + xbox[2];
    ybox[0] = ybox[1] - ybox[2];
#ifdef COLOR
    color = getparam("color");
#endif
    if (hasvalue("frame"))
        frame = getparam("frame");
    else
        frame = NULL;
}

setrange(rval, rexp)
real rval[];
string rexp;
{
    char *cptr;

    cptr = strchr(rexp, ':');
    if (cptr != NULL) {
        rval[0] = atof(rexp);
	rval[1] = atof(cptr+1);
    } else {
        rval[0] = 0.0;
	rval[1] = atof(rexp);
    }
    rval[2] = rval[1] - rval[0];
}

extern btrproc btrtrans(string);
extern btiproc btitrans(string);

compfuncs()
{

    xfunc = btrtrans(xvar);
    yfunc = btrtrans(yvar);
    vfunc = btitrans(visib);
    pfunc = btrtrans(psize);
#ifdef COLOR
    cfunc = btrtrans(color);
#endif
#ifdef VECTOR
    vxfunc = btrtrans(vxvar);
    vyfunc = btrtrans(vyvar);
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

bool scansnap()
{
    bool success;
    int i;
    real *pptr, *xptr;

    success = FALSE;
    while (! success) {
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
    return (TRUE);
}

plotbox()
{
    char msg[128];
    int nticks;
    real ticks[MAXTICKS];

    if (formal) {
	formalaxis = TRUE;
	xaxisvar.labdn = 0.44;
	xaxisvar.szlab = 0.40;
	yaxisvar.numdn = 0.24;
	yaxisvar.labdn = 0.24;
	yaxisvar.szlab = 0.40;
    }
    if (! nobox) {
        if (hasvalue("xticks"))
	    setticks(ticks, &nticks, getparam("xticks"));
        else {
	    ticks[0] = xrange[0];
	    ticks[1] = xrange[1];
	    nticks = - getiparam("nxticks");
	} 
	
        xaxis(xbox[0], ybox[0], xbox[2], ticks, nticks, xtrans,
              iy==nxy[1]-1 ? xlabel : NULL);
	xaxis(xbox[0], ybox[1], xbox[2], ticks, nticks, xtrans, NULL);
        if (hasvalue("yticks"))
	    setticks(ticks, &nticks, getparam("yticks"));
	else {
	    ticks[0] = yrange[0];
	    ticks[1] = yrange[1];
	    nticks = - getiparam("nyticks");
	} 
        yaxis(xbox[0], ybox[0], ybox[2], ticks, nticks, ytrans,
	      ix==0 ? ylabel : NULL);
	yaxis(xbox[1], ybox[0], ybox[2], ticks, nticks, ytrans, NULL);
	if (! formal && ix==0 && iy==0) {
	    sprintf(msg, "File: %s", input);
	    pltext(msg, xbox[0], ybox[1] + 0.4, 0.32, 0.0);
	    if (timeptr != NULL) {
		if (trakflag || nxy[0]>1 || nxy[1]>1)
		    sprintf(msg, "Times: %s", times);
		else
		    sprintf(msg, "Time: %8.3f", *timeptr);
		pltext(msg, xbox[0] + 10.0, ybox[1] + 0.4, 0.32, 0.0);
	    }
	}
    } else {
	if (! formal && timeptr != NULL) {
	    sprintf(msg, "%.2f", *timeptr);
	    pltext(msg, xbox[1] - 1.0, ybox[1] - 1.0, 0.24, 0.0);
	}
    }
}

setticks(tiks, ntik, tikstr)
real tiks[];
int *ntik;
string tikstr;
{
#if 0
    /* this section has memory leak; allocate, but never free's */
    string *burststring(), *tptr;

    tptr = burststring(tikstr, ", ");
    *ntik = 0;
    while (*tptr != NULL)
	tiks[(*ntik)++] = atof(*tptr++);
#else
    *ntik = nemoinpr(tikstr, tiks, MAXTICKS);
    if (*ntik < 0) error("(5d) Parsing %s",*ntik, tikstr);
#endif
}

plotsnap()
{
    real t, *mp, *psp, *pp, *ap;
    int vismax, visnow, i, vis, icol;
    real psz, col, x, y, vx, vy;
    Body b;

    t = (timeptr != NULL ? *timeptr : 0.0);	/* get current time value   */
    CLRV(Acc(&b));				/* zero unsupported fields  */
    Key(&b) = 0;
    visnow = vismax = 0;
    do {					/* loop painting layers     */
	visnow++;				/*   make next layer visib. */
	mp = massptr;				/*   (re)set data pointers  */
	psp = phaseptr;
	pp = phiptr;
	ap = auxptr;
	for (i = 0; i < nbody; i++) {		/*   loop over all bodies   */
	    Mass(&b) = (mp != NULL ? *mp++ : 0.0);
						/*     set mass if supplied */
	    SETV(Pos(&b), psp);			/*     always set position  */
	    psp += NDIM;			/*     and advance p.s. ptr */
	    SETV(Vel(&b), psp);			/*     always set velocity  */
	    psp += NDIM;			/*     and advance ptr      */
	    Phi(&b) = (pp != NULL ? *pp++ : 0.0);	
	    Aux(&b) = (ap != NULL ? *ap++ : 0.0);	
	    					/*     set phi,aux if given */
	    vis = (*vfunc)(&b, t, i);		/*     evaluate visibility  */
	    vismax = MAX(vismax, vis);		/*     remember how hi to go*/
	    if (vis == visnow) {		/*     if body is visible   */
		x = xtrans((*xfunc)(&b, t, i));	/*       evaluate x,y coords*/
		y = ytrans((*yfunc)(&b, t, i));
		if (xbox[0] < x && x < xbox[1] && ybox[0] < y && y < ybox[1]) {
		    psz = (*pfunc)(&b, t, i);	/*         eval point size  */
#ifdef COLOR
		    col = (*cfunc)(&b, t, i);
		    icol = 1 + (plncolors() - 2) *
			         MAX(0.0, MIN(1.0, col));
		    plcolor(icol);
#endif
		    if (psz == 0.0)
			plpoint(x, y);
		    else if (psz < 0.0) {
			if (fillcircle)
			    plcross(x, y, - psz);
                        else
                            plcross(x, y, psz);
		    } else if (psz > 0.0) {
			if (fillcircle)
			    plcircle(x, y, -psz);
			else
			    plcircle(x, y, psz);
		    }
#ifdef VECTOR
		    vx = (*vxfunc)(&b, t, i); /*  evaluate vector */
		    vy = (*vyfunc)(&b, t, i);
		    vx *= scale;
		    vy *= scale;
		    arrow(x,y,vx,vy);
#endif		    
		}
	    }
	}
    } while (visnow < vismax);			/* until final layer done   */
#ifdef COLOR
    plcolor(32767);				/* reset to white */
#endif
}

#ifdef HACKTRACK

#define X0  -0.61
#define X1   1.22

#define Y0   0.1
#define Y2  -2.5

#define L0   0.05

plottrack()
{
    int i;
    real x, y, dx, dy;

#ifdef COLOR
    plcolor(32767);				/* reset to white */
#endif
    for (i = 0; i <= 128; i++) {
	x = X0 + X1 * i / 128.0;
	y = Y0 + Y2 * x*x;
	if (i == 0)
	    plmove(xtrans(y), ytrans(- x));
	else
	    plline(xtrans(y), ytrans(- x));
    }
    dx = - L0 / sqrt(1.0 + 4 * Y2*Y2 * x*x);
    dy = 2 * Y2 * x * dx;
    plmove(xtrans(y + 0.94*dy - 0.34*dx), ytrans(- x - 0.94*dx - 0.34*dy));
    plline(xtrans(y), ytrans(- x));
    plline(xtrans(y + 0.94*dy + 0.34*dx), ytrans(- x - 0.94*dx + 0.34*dy));
    for (i = 0; i <= 128; i++) {
	x = X0 + X1 * i / 128.0;
	y = Y0 + Y2 * x*x;
	if (i == 0)
	    plmove(xtrans(- y), ytrans(x));
	else
	    plline(xtrans(- y), ytrans(x));
    }
    dx = - L0 / sqrt(1.0 + 4 * Y2*Y2 * x*x);
    dy = 2 * Y2 * x * dx;
    plmove(xtrans(- y - 0.94*dy + 0.34*dx), ytrans(x + 0.94*dx + 0.34*dy));
    plline(xtrans(- y), ytrans(x));
    plline(xtrans(- y - 0.94*dy - 0.34*dx), ytrans(x + 0.94*dx - 0.34*dy));
}

#endif

real xtrans(real x)
{
    return xbox[0] + xbox[2] * (x - xrange[0]) / xrange[2];
}

real ytrans(real y)
{
    return ybox[0] + ybox[2] * (y - yrange[0]) / yrange[2];
}

/* draw an arrow, anchored at (xcm,ycm) with physical units vx,vy 
 *   with a head, at 30deg
 */

void arrow(real xcm, real ycm, real vx, real vy)
{
    real dx, dy;
    
    plmove(xcm,ycm);
    xcm += vx;
    ycm += vy;
    plline(xcm, ycm);
#if 1
    dx = -vx *0.866 + vy*0.5;
    dy = -vx *0.5 - vy*0.866;
    plline(xcm + 0.2*dx, ycm + 0.2*dy);

    plmove(xcm,ycm);
    dx = -vx*0.866 -vy*0.5;
    dy = vx*0.5 - vy*0.866;
    plline(xcm + 0.2*dx, ycm + 0.2*dy);
#endif
}


scrdump(frameno)
int frameno;
{
    char s[64];

#if 1
    /* new method: let yapp figure it out how to make movie frames */
    sprintf(s, "%s.%d", frame, frameno);
    pl_screendump(s);
#else
    /* Old method: suntools only */
    sprintf(s, "screendump %s.%d", frame, frameno);
    dprintf(0,"%s\n",s);
    system(s);
#endif
}

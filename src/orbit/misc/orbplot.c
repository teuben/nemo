/*
 * ORBPLOT.C: read a orbit and plot particle trajectory
 *
 *	14-jul-87	V1.0: cloned off SNAPPLOT-PJT
 *	28-jul-87	V2.0: new orbit(5)              	PJT
 *	 1-jun-88	V2.1: new filestruct, no code mod's 	PJT
 *      22-may-90       V2.2: read multiple orbits              PJT
 *      19-jun-91       V2.3: handle multiple files too         PJT
 *	 7-mar-92	V2.3a: gcc happy			pjt
 *	24-may-92           b: included <potential.h> now	Pjt
 *	23-mar-95	    c: index->strchr			pjt
 *	17-apr-95           d: -DNO_NBODY, simple orbrans	pjt
 *      20-feb-97           d: fixed for SINGLEPREC             pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <orbit.h>
#include <vectmath.h>
#include <yapp.h>
#include <axis.h>

string defv[] = {               /* DEFAULT INPUT PARAMETERS */
    "in=???\n                      input file name ",
    "times=all\n		   range of times to plot ",
    "nplot=1\n			   skip intermaediates? ",
    "xvar=x\n			   x-axis plotting variable ",
    "xlabel=\n			   x-axis label; defaults to xvar ",
    "xrange=-2.0:2.0\n		   x-axis variable range ",
    "yvar=y\n			   y-axis plotting variable ",
    "ylabel=\n			   y-axis label; defaults to yvar ",
    "yrange=-2.0:2.0\n		   y-axis variable range ",
    "visib=1\n			   determine visibility of points ",
    "psize=0\n			   point type and size ",
#ifdef COLOR
    "color=0\n			   determine color of points ",
#endif
    "maxsteps=5000\n               max. number of steps per orbit",
    "VERSION=2.3e\n                20-feb-97 PJT",
    NULL,
};

string usage = "Plot an orbit trajectory";

#ifndef HUGE
# define	HUGE  1.0e20
#endif

string input, otimes;
orbitptr optr;
int    nplot, maxsteps;
string xvar, yvar;
string xlabel, ylabel;
string visib, psize, color;
int    xvar_idx, yvar_idx;      /* non-bodytrans code */

real xrange[2], yrange[2], trange[2];
real xtrans(real),  ytrans(real);

int    get_idx(string);
void   setparams(void), setrange(real *, string), compfuncs(void);
void   plot_path(orbitptr, real, real, int);

#ifndef NO_NBODY
rproc xfunc;                    /* only used if bodytrans enabled */
rproc yfunc;
iproc vfunc;
rproc pfunc;
rproc cfunc;
extern rproc btrtrans();
extern iproc btitrans();
#endif

extern string *burststring(string, string);

void nemo_main()
{
    stream istr;
    char msg[128];
    int  i, ndim;
    string *fn;
    
    setparams();
    compfuncs();	/* btrrtans snapplot-like interface */

    plinit("***", 0.0, 20.0, 0.0, 20.0);

    xaxis( 2.0,  2.0, 16.0, xrange, -7, xtrans, xlabel);
    xaxis( 2.0, 18.0, 16.0, xrange, -7, xtrans, NULL);
    yaxis( 2.0,  2.0, 16.0, yrange, -7, ytrans, ylabel);
    yaxis(18.0,  2.0, 16.0, yrange, -7, ytrans, NULL);
    sprintf(msg, "File: %s", input);
    pltext(msg, 2.0, 18.4, 0.32, 0.0);

    optr = NULL;
    ndim = NDIM;
    allocate_orbit(&optr,ndim,maxsteps);    /* allocate a large orbit */

    fn = burststring(input," ,");
    for (i=0; fn[i]; i++) {             /* loop over all files */
        istr = stropen(fn[i], "r");             /* open orbit file */
        while (read_orbit(istr,&optr)) {         /* while an orbit found....*/
            dprintf(0,"Read orbit with %d phase-points\n",Nsteps(optr));
            plot_path(optr,trange[0],trange[1],nplot);
            Nsteps(optr) = maxsteps;             /* reset for next orbit */
        }
        strclose(istr);
    }

    plstop();
}

void setparams()
{
    input = getparam("in");
    otimes = getparam("times");
    if (strcmp(otimes,"all")==0) {
        trange[0] = -HUGE;
        trange[1] =  HUGE;
    } else
    	setrange(trange, otimes);

    nplot=getiparam("nplot");    

    xvar = getparam("xvar");
    xvar_idx = get_idx(xvar);
    if (hasvalue("xlabel"))
        xlabel = getparam("xlabel");
    else
        xlabel = xvar;
    setrange(xrange, getparam("xrange"));

    yvar = getparam("yvar");
    yvar_idx = get_idx(yvar);
    if (hasvalue("ylabel"))
        ylabel = getparam("ylabel");
    else
        ylabel = yvar;
    setrange(yrange, getparam("yrange"));

    visib = getparam("visib");
    psize = getparam("psize");
#ifdef COLOR
    color = getparam("color");
#endif
    maxsteps = getiparam("maxsteps");
}

void setrange(real *rval, string rexp)
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
}

int get_idx(string var)
{
    /* a poor man's bodytrans */

    if (streq(var,"x"))
        return 0;
    else if (streq(var,"y"))
        return 1;
    else if (streq(var,"z"))
        return 2;
    else if (streq(var,"vx"))		/* remainder only if ndim=NDIM=3 */
        return 3;
    else if (streq(var,"vy"))
        return 4;
    else if (streq(var,"vz"))
        return 5;
    else
        error("this version does not support var=%s",var);
    return 0;
}

void compfuncs()
{
#ifndef NO_NBODY

    xfunc = btrtrans(xvar);
    yfunc = btrtrans(yvar);
    vfunc = btitrans(visib);
    pfunc = btrtrans(psize);
#ifdef COLOR
    cfunc = btrtrans(color);
#endif
#endif
}

void plot_path(orbitptr o, real tstart, real tend, int n)
{
	int i, plotfirst = 1;
	real x,y;

	for (i=0; i<(Nsteps(o)-n); i += n) {
	    if ((tstart<Torb(optr,i)) && (Torb(optr,i)<tend)) {
		x = xtrans(Posorb(o,i,xvar_idx));
		y = ytrans(Posorb(o,i,yvar_idx));
		if (plotfirst) {
		    plpoint (x,y);			/* mark first point */
	 	    plotfirst = 0;
		}
		plmove (x,y);
		x = xtrans(Posorb(o,i+n,xvar_idx));
		y = ytrans(Posorb(o,i+n,yvar_idx));
		plline (x,y);
	     }
	}
}

#if 0
		/* THIS NEEDS TO BE IMPLEMENTED: see new snapplot code */
var_plot_path()
{
    int i, key, vis;
    real *pos, *vel, *acc, avec[NDIM], pot, mass, aux;
    real psz, col, x, y;

    CLRV(avec);
    mass = aux = 0.0;
    key = 0;
    for (i = 0; i < nobj; i++) {
	pos = phase[i][0];
	vel = phase[i][1];
	acc = avec;
	pot = (phi != NULL ? *(phi + i) : 0.0);
	vis = (*vfunc)(pos,vel,acc,pot,mass,aux,key,t,i);
	if (vis != 0) {
	    x = xtrans((*xfunc)(pos,vel,acc,pot,mass,aux,key,t,i));
	    y = ytrans((*yfunc)(pos,vel,acc,pot,mass,aux,key,t,i));
	    psz = (*pfunc)(pos,vel,acc,pot,mass,aux,key,t,i);
#ifdef COLOR
	    col = (*cfunc)(pos,vel,acc,pot,mass,aux,key,t,i);
	    plcolor(vis == 1 ? col : 2.0);	/* highlight if vis > 1 */
#endif
	    if (2.0 < x && x < 18.0 && 2.0 < y && y < 18.0) {
		if (psz == 0.0)
		    plpoint(x, y);
	        else if (psz > 0.0)
		    plcircle(x, y, psz);
	        else if (psz < 0.0)
		    plcross(x, y, psz);
	    }
	}
    }
#ifdef COLOR
    plcolor(2.0);				/* reset to white */
#endif
}
#endif

real xtrans(real x)
{
    return (2.0 + 16.0 * (x - xrange[0]) / (xrange[1] - xrange[0]));
}

real ytrans(real y)
{
    return (2.0 + 16.0 * (y - yrange[0]) / (yrange[1] - yrange[0]));
}

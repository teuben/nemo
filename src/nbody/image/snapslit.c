/* 
 *  SNAPSLIT: put a slit accross an Nbody system and compute 3 moments
 *
 *          modeled after snapplot/snapgrid
 *
 *
 *	27-May-87  V1.0 original program   	P.J. Teuben
 *	16-Jun-87  V1.0a recreated after disk-scratch   PJT
 *	 8-jun-88  V1.1  new filestruct                 PJT
 *	21-feb-89  V1.4  tab=t does not plot anymore	PJT
 *	22-apr-92  V2.0  NEMO V2.x; major overhaul: xvar=yvar= etc.  PJT
 *	21-jul-95  V2.0a bugfix - calling sprintf for strcpy   Dave Shone
 *	 4-mar-97      b fix for SINGLEPREC
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

#include <yapp.h>
#include <axis.h>

string defv[] = {
    "in=???\n		    input filename (snapshot)",
    "xvar=x\n		    X variable to grid",
    "yvar=y\n		    Y variable to grid",
    "zvar=-vz\n		    Z variable to use in moment calculation",
    "evar=m\n		    Emission variable to grid",	
    "origin=0,0\n	    origin of slit on 'sky' (xvar,yvar)",
    "pa=0\n		    position angle ; N=0 E=90 ",
    "width=0.1\n	    width of slit 'sky' units ",
    "length=10.0\n	    length of slit ",
    "cell=0.1\n		    cellsize ",
    "smooth=0.25,0.5,0.25\n symmetric smoothing array along slit",
    "nsmooth=1\n	    number of smoothings",
    "gauss=\n		    If supplied, this is the gaussian FWHM beam",
    "mmax=0.0\n		    max emission to plot",
    "vmin=0.0\n	            min mean to plot",
    "vmax=0.0\n             max mean to plot",
    "smax=0.0\n             max dispersion to plot",
    "tab=f\n                Need a table? If yes, no plot",
    "VERSION=2.0b\n         4-mar-97 PJT",
    NULL,
};

string usage="slit spectra of N-body systems: 3 moments";

#define RPD (3.141592/180.0)

string	infile;				/* input file */
stream  instr;

Body *btab;                             /* snapshot */
int    nobj;
real tsnap;
string headline;

#define MSLIT 4096			/* maximum # pixels in a slit */
real v0star[MSLIT];			/* zeroth moment*/
real v1star[MSLIT];			/* first moment */
real v2star[MSLIT];			/* second moment */
int    nslit;				/* actual slit length */

real origin[2];			/* the slit in 'sky' coordinates */
real pa;
real slit_width;
real slit_cell;
real slit_length;

#define MSMOOTH 101

real smooth[MSMOOTH];                 /* full symmetric smoothing array */
int    lsmooth;				/* actual smoothing length */
int    nsmooth;				/* number of smoothings */
real gauss;                           /* gaussian width, if gauss is used */

bool   Qtab;                            /* need a table output instead ? */

rproc xvar, yvar, zvar, evar;           /* gridding variables */

real  mmax, vmin, vmax, smax;             /* plot min and max */
real xplot[2],   yplot[2];		    /* ranges in plot variables */
char   xlabel[80], ylabel[80], plabel[132];  /* labels */

local real xtrans(real), ytransm(real), ytransv1(real), ytransv2(real);

extern rproc btrtrans(string);


nemo_main()
{
    setparams();

    instr = stropen (infile, "r");
    read_snap();
    strclose(instr);

    out_slit();        /* no plot when Qtab=t */
}

setparams()
{
    string  tmpstr;
    int tmpint, i;
    real  x, fx, sum, sigma;

    infile = getparam ("in");

    xvar = btrtrans(getparam("xvar"));
    yvar = btrtrans(getparam("yvar"));
    zvar = btrtrans(getparam("zvar"));
    evar = btrtrans(getparam("evar"));

    if (hasvalue("origin")) {
        if (nemoinpr(getparam("origin"),origin,2) != 2)
            error ("origin= must have two values, or no value");
    }

    pa=getdparam("pa")*RPD;         /* input must be in degrees */
    slit_width = getdparam("width");
    slit_cell = getdparam("cell");
    slit_length = getdparam("length");
    nslit = slit_length/slit_cell;
    dprintf (1,"#pixels along slit = %d\n",nslit);
    if (nslit>MSLIT)
        error ("too many slit cells, maximum %s",MSLIT);

    if(hasvalue("gauss")) {         /* gaussian beam */
                gauss = getdparam("gauss");
                sigma = gauss/2.355;    /* FWHM = 2 sqrt(2 ln(2)) * sigma */
                lsmooth=1;              /* count how many we need */
                x=0;
                sum=1.0;                /* unnormalized value at x=0 */
                do {
                        x += slit_cell;
                        lsmooth += 2;
                        fx = exp(-0.5*sqr(x/sigma));            
                        dprintf (1,"making beam %f %f\n",x,fx);
                        sum += 2*fx;
                } while ((fx>0.01) && (lsmooth<MSMOOTH));  /* cutoff at 1% */
                sum = 1.0/sum;
                dprintf (0,"Gaussean beam will contain %d points \n",lsmooth);
                if ((lsmooth==MSMOOTH) && (fx>0.01)) 
                        dprintf(0,"Warning: beam cutoff at %d points\n",MSMOOTH);
                for (i=0; i<lsmooth; i++) {
                        x = (i-(lsmooth-1)/2)*slit_cell/sigma;
                        smooth[i] = exp(-0.5*x*x) * sum;   /* normalize */
                }
                nsmooth = getiparam("nsmooth");
                if (nsmooth<=0) nsmooth=1;
    } else {                        /* beam by hand */
                if ( (lsmooth=nemoinpr(getparam("smooth"),smooth,MSMOOTH)) < 0)
                    error("Parameter smooth=%s",getparam("smooth"));
                nsmooth = getiparam("nsmooth");
                if (nsmooth<=0) nsmooth=1;
    }
        
    mmax=getdparam("mmax");
    vmin=getdparam("vmin");
    vmax=getdparam("vmax");
    smax=getdparam("smax");
    Qtab = getbparam("tab");
}

read_snap()
{                               
    int    bits;

    get_history(instr);    
    get_snap(instr,&btab,&nobj,&tsnap,&bits);
}

#define  SYMBOLSIZE  0.1		/* size of plotsymbol in 'cm' */


out_slit()
{
    real xsky, ysky, vrad, inv_surden, sigma, mass;
    real xslit, yslit, xplt, yplt, sinpa, cospa;
    real m_max, v_min, v_max, s_max;	      /* local min/max */
    int    i, islit;
    Body *bp;

    for (islit=0; islit<nslit; islit++)     /* reset local variables */
	v0star[islit] = v1star[islit] = v2star[islit] = 0.0;
    m_max = v_min = v_max = s_max = 0.0;
    inv_surden = 1.0 / (slit_width*slit_cell);
    sinpa = sin(pa); cospa = cos(pa);

    for(bp=btab, i=0; i<nobj; bp++, i++) {      /* loop over all particles */
 	xsky = xvar(bp,tsnap,i);
 	ysky = yvar(bp,tsnap,i);
 	vrad = zvar(bp,tsnap,i);
 	mass = evar(bp,tsnap,i) * inv_surden;
 	
	xsky -= origin[0];			/* translate to slit origin */
	ysky -= origin[1];
	xslit = -cospa*ysky + sinpa*xsky;	/* and rotate to slit frame */
	yslit =  sinpa*ysky + cospa*xsky;	/* !!! check signs !!! */

	if (fabs(yslit) > 0.5*slit_width) 
	   continue;			/* not in slit */

	islit =  (xslit+0.5*slit_length)/slit_cell;
	if (islit<0 || islit>=nslit)
	   continue;			/* not in slit */

	v0star[islit] += mass;
	v1star[islit] += vrad * mass;
	v2star[islit] += sqr(vrad) * mass;
     } /*-- end particles loop --*/


     while (nsmooth-- > 0) {            	/* convolution */
     	dprintf (0,"Convolving with %d-length beam: ",lsmooth);
     	for (i=0; i<lsmooth; i++) 
            dprintf (0,"%f ",smooth[i]);
    	convolve (v0star, nslit, smooth, lsmooth);
    	convolve (v1star, nslit, smooth, lsmooth);
    	convolve (v2star, nslit, smooth, lsmooth);
     	dprintf (0,"\n");
    }
    	
    for (islit=0; islit<nslit; islit++) { /* moment analysis: M, MV and MV^2 */
	if (v0star[islit]==0.0)
	    continue;		/* no data - skip to next pixel */
	v1star[islit] /= v0star[islit];
	sigma = v2star[islit]/v0star[islit] - sqr(v1star[islit]);
	if (sigma<0.0) {        /* should never happen */
	    warning("islit=%d sigma^2=%e < 0 !!!\n",islit,sigma);
	    v2star[islit] = 0.0;
	    continue;		/* something really wrong */
	}
	v2star[islit] = sqrt(sigma);

	if (v0star[islit] > m_max)  m_max = v0star[islit];
	if (v1star[islit] < v_min)  v_min = v1star[islit];
	if (v1star[islit] > v_max)  v_max = v1star[islit];
	if (v2star[islit] > s_max)  s_max = v2star[islit];
	if (Qtab) {
	    xslit = islit*slit_cell;
	    printf ("%g %g %g %g\n",
			xslit,v0star[islit], v1star[islit], v2star[islit]);
	}
    } /* for(islit) */
 
    if (Qtab)
        return(0);

    plinit ("***", 0.0, 20.0, 0.0, 20.0);

	/* reset default autoscales to user supplied if necessary */
    if (mmax==0.0) mmax=m_max;
    if (vmin==0.0) vmin=v_min;
    if (vmax==0.0) vmax=v_max;
    if (smax==0.0) smax=s_max;
    dprintf (0,"mmax=%f vmin=%f vmax=%f smax=%f  reset to:\n",m_max,v_min,v_max,s_max);
    if (mmax==0) mmax=1;
    if (vmin==0 && vmax==0) vmax=1;
    if (smax==0) smax=1;
    dprintf (0,"mmax=%f vmin=%f vmax=%f smax=%f           \n",mmax, vmin, vmax, smax);
    
	/* general plot header */
    sprintf (plabel,"File: %s; var{%s,%s,%s,%s} slit{%s %s %s %s}",
		infile,getparam("xvar"),getparam("yvar"),
		getparam("zvar"),getparam("evar"),
		getparam("origin"),getparam("pa"),getparam("width"),
		getparam("length"),getparam("cell"));
		
    pltext (plabel,2.0,18.4, 0.32, 0.0);
#if 0
    if (*headline!=NULL)				/* identification */
	pltext (headline,2.0,19.0,0.25,0.0);
#endif

    xplot[0] = -0.5*slit_length;            /* PLOT1: upper panel */
    xplot[1] =  0.5*slit_length;
    sprintf(xlabel,"slit: {x=%s,y=%s}",getparam("xvar"),getparam("yvar"));

    yplot[0]=0.0; 
    yplot[1]=mmax;  
    strcpy (ylabel,"mass surface density");

    xaxis ( 2.0,12.0, 16.0, xplot, -7, xtrans,  NULL);
    xaxis ( 2.0,17.0, 16.0, xplot, -7, xtrans,  NULL);
    yaxis ( 2.0,12.0,  5.0, yplot, -3, ytransm, ylabel);
    yaxis (18.0,12.0,  5.0, yplot, -3, ytransm, NULL);

    for (islit=0; islit<nslit; islit++) {
	xplt = xtrans (-0.5*slit_length + (islit+0.5)*slit_cell);
	yplt = ytransm (v0star[islit]);
	plbox (xplt, yplt, SYMBOLSIZE);
    }
                                           
    yplot[0]=vmin;                      /* PLOT2: middle panel */
    yplot[1]=vmax;  
    strcpy (ylabel,"velocity");

    xaxis (2.0, 7.0, 16.0, xplot, -7, xtrans,   NULL);	/* line ?? */
    yaxis (2.0, 7.0,  5.0, yplot, -3, ytransv1, ylabel);
    yaxis (18.0,7.0,  5.0, yplot, -3, ytransv1, NULL);

    for (islit=0; islit<nslit; islit++) {
	xplt = xtrans (-0.5*slit_length + (islit+0.5)*slit_cell);
	yplt = ytransv1 (v1star[islit]);
	plcross (xplt, yplt, SYMBOLSIZE);
    }
    if (vmin<0.0 || vmax>0.0) {
       plltype (1,2);	/* dashed line at v=0 */
       plmove (xtrans(xplot[0]), ytransv1(0.0));
       plline (xtrans(xplot[1]), ytransv1(0.0));
       plltype (1,1);
    }

    yplot[0]=0.0;                       /* PLOT3: bottom panel */
    yplot[1]=smax; 
    strcpy (ylabel,"velocity dispersion");
    xaxis (2.0, 2.0, 16.0, xplot, -7, xtrans,   xlabel);
    yaxis (2.0, 2.0,  5.0, yplot, -3, ytransv2, ylabel);
    yaxis (18.0,2.0,  5.0, yplot, -3, ytransv2, NULL);

    for (islit=0; islit<nslit; islit++) {
	xplt = xtrans (-0.5*slit_length + (islit+0.5)*slit_cell);
	yplt = ytransv2 (v2star[islit]);
	plcross (xplt, yplt, -SYMBOLSIZE);
    }

    plstop();
}

/*********** SOME LOCAL UTILITIES , for plotting ***********/

real xtrans(x)
real x;
{
    return (2.0 + 16.0*(x-xplot[0])/(xplot[1]-xplot[0]) );
}


real ytransm(y)	/* upper panel: mass surface density */
real y;
{
    return (12.0 + 5.0*(y-yplot[0])/(yplot[1]-yplot[0]) );
}

real ytransv1(y)	/* middle panel: velocity */
real y;
{
    return (7.0 + 5.0*(y-yplot[0])/(yplot[1]-yplot[0]) );
}

real ytransv2(y)	/* lower panel: velocity dispersion */
real y;
{
    return (2.0 + 5.0*(y-yplot[0])/(yplot[1]-yplot[0]) );
}


#define MSIZE MSLIT

/* CONVOLVE: 1D, the naive approach */

convolve (a, na, b, nb)
real a[], b[];
int    na,  nb;
{
    real c[MSIZE];
    int    i, j, k;
	
    if (na>MSIZE) {
	warning("convolve: MSIZE=%d array too small for input %d\n",MSIZE,na);
	return(0);
    }
	
    for (i=0; i<na; i++) {
		c[i] = a[i];		/* copy array */
		a[i] = 0.0;		/* reset for accumulation */
    }
		
    for (i=0; i<na; i++) 
        for (j=0; j<nb; j++) {
            k = i + j - (nb-1)/2;
            if (k>=0)
                if (k<na)
                    a[k] += b[j]*c[i];
                else
                    continue;
        }
			
    return(1);
}

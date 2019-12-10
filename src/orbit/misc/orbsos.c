/*
 * ORBSOS:	 tabulate an orbit as a surface of section
 *
 *	9-jun-88  PJT  Created - (still named mksos in a beta version)
 *     26-mar-92  V1.1 put table version into system - named orbsos.c   PJT
 *     24-may-92  V1.1a added <potential.h>                             PJT
 *     30-dec-93  V1.2 added pabs,vabs=t|f                              PJT
 *    9-dec-2019  V1.3 symm option to deal with e.g. hh64               PJT
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <potential.h>
#include <orbit.h>

string defv[] = {
    "in=???\n		Input orbit",
    "mode=x\n           Mode (x|y|xy)",
    "pabs=f\n           Show absolute values of positions?",
    "vabs=f\n           Show absolute values of velocities?",
    "symm=t\n           Symmetric potential?",
    "VERSION=1.3\n  	9-dec-2019 PJT",
    NULL,
};

string usage="Tabulate the SOS coordinates of an orbit";


string	infile;			/* file names */
stream  instr;			/* file streams */

orbitptr optr;
int xmode=0, ymode=0;

local void plot_sos2(orbitptr optr, bool symm, bool pabs,bool vabs);

void nemo_main(void)
{
    string mode;
    bool symm,pabs,vabs;

    infile = getparam("in");
    instr = stropen (infile,"r");
    mode = getparam("mode");
    xmode = scanopt(mode,"y");
    ymode = scanopt(mode,"x");
    if (xmode==0 && ymode==0) 
        error("mode must contain at least one of 'x' or 'y'");
    if (xmode) dprintf(1,"SOS X-VX will be computed\n");
    if (ymode) dprintf(1,"SOS Y-VY will be computed\n");
    pabs = getbparam("pabs");
    vabs = getbparam("vabs");
    symm = getbparam("symm");
    
    optr=NULL;
    read_orbit(instr,&optr);
    plot_sos2(optr,symm,pabs,vabs);
    strclose(instr);
}


/* 	
 * plots a surface of section from a 2D orbit : classical method
 *
 */
 
#define TEST 1			/* if defined, no plotting */

static int sos_count=0;
static real *pos_sos=NULL;
static real *vel_sos=NULL;

int lzs(real x, real y, real u, real v)
{
  real lz = x*v - y*u;
  return SGN(lz);
}

void plot_sos2(orbitptr optr, bool symm, bool pabs, bool vabs)
{
	real xold,xnew, yold,ynew,f, xtest, ytest;
	real yplot,vplot,tplot;
	int  i, lz;
	int lz0 = lzs(Xorb(optr,0),Yorb(optr,0),Uorb(optr,0),Vorb(optr,0));
	
	dprintf(1,"Orbit with SGN(lz) = %d\n",lz0);
	
#ifndef TEST
		/* allocate memory if necessary */
	if (pos_sos==NULL) {
		pos_sos = (real *) allocate(Nsteps(optr)/20*sizeof(real));
		vel_sos = (real *) allocate(Nsteps(optr)/20*sizeof(real));
	}
#endif
		/* walk along orbit path and determine SOS datapoints */
        xtest = 1;    /* this will change sign if the X axis was crossed */
        ytest = 1;    /* this will change sign if the Y axis was crossed */
	for (i=1; i<Nsteps(optr); ++i) {
		xold=Xorb(optr,i-1);
		xnew=Xorb(optr,i);
		yold=Yorb(optr,i-1);
		ynew=Yorb(optr,i);
                if(xmode) xtest=xold*xnew;
                if(ymode) ytest=yold*ynew;

                if (xtest<0) {  /* crossing X axis: get the Y-VY's */
                    f = xold/(xnew-xold);
                    tplot = (1+f)*Torb(optr,i-1) - f*Torb(optr,i);
                    yplot = (1+f)*Yorb(optr,i-1) - f*Yorb(optr,i);
                    vplot = (1+f)*Vorb(optr,i-1) - f*Vorb(optr,i);
		    lz = lzs(0,yplot,Uorb(optr,i-1),0);
                    if (lz != lz0) {
		        if (symm) {
			    yplot = -yplot;
			    vplot = -vplot;
			} else
   			    continue;
                    }
                    if (pabs) yplot = ABS(yplot);
                    if (vabs) vplot = ABS(vplot);
                    printf ("%d %f %f %f\n",sos_count, tplot, yplot, vplot);
                    sos_count++;
                } 
                if (ytest<0) {  /* crossing Y axis: get the X-VX */
                    f = yold/(ynew-yold);
                    tplot = (1+f)*Torb(optr,i-1) - f*Torb(optr,i);
                    yplot = (1+f)*Xorb(optr,i-1) - f*Xorb(optr,i);
                    vplot = (1+f)*Uorb(optr,i-1) - f*Uorb(optr,i);
		    lz = lzs(yplot,0,0,Vorb(optr,i-1));
                    if (lz != lz0) {
		        if (symm) {
                            yplot = -yplot;
                            vplot = -vplot;
			} else
			    continue;
                    }
                    if (pabs) yplot = ABS(yplot);
                    if (vabs) vplot = ABS(vplot);
                    printf ("%d %f %f %f\n",sos_count, tplot, yplot, vplot);
                    sos_count++;
                }
        }
	
}

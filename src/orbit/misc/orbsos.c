/*
 * ORBSOS:	 tabulate an orbit as a surface of section
 *
 *	9-jun-88  PJT  Created - (still named mksos in a beta version)
 *     26-mar-92  V1.1 put table version into system - named orbsos.c   PJT
 *     24-may-92  V1.1a added <potential.h>                             PJT
 *     30-dec-93  V1.2 added pabs,vabs=t|f                              PJT
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
    "VERSION=1.2\n  	30-dec-93 PJT",
    NULL,
};

string usage="Generate SOS from an orbit";


string	infile;			/* file names */
stream  instr;			/* file streams */

orbitptr optr;
int xmode=0, ymode=0;

nemo_main ()
{
    string mode;
    bool pabs,vabs;

    infile = getparam("in");
    instr = stropen (infile,"r");
    mode = getparam("mode");
    xmode = scanopt(mode,"y");
    ymode = scanopt(mode,"x");
    if (xmode==0 && ymode==0) 
        error("mode must contain at least one of 'x' or 'y'");
    pabs = getbparam("pabs");
    vabs = getbparam("vabs");
    
    optr=NULL;
    read_orbit (instr,&optr);
    plot_sos2 (optr,pabs,vabs);
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

plot_sos2(optr,pabs,vabs)
orbitptr optr;
bool pabs,vabs;
{
	real xold,xnew, yold,ynew,f, xtest, ytest;
	real yplot,vplot,tplot;
	int  i;
	
#ifndef TEST
		/* allocate memory if necessary */
	if (pos_sos==NULL) {
		pos_sos = (real *) malloc(Nsteps(optr)/20*sizeof(real));
		vel_sos = (real *) malloc(Nsteps(optr)/20*sizeof(real));
		if (vel_sos==NULL) {
			printf ("No memory\n");
			exit(1);
		}
	}
#endif
		/* walk along orbit path and determine SOS datapoints */
        xtest = 1;
        ytest = 1;
	for (i=1; i<Nsteps(optr); i++) {
		xold=Xorb(optr,i-1);
		xnew=Xorb(optr,i);
		yold=Yorb(optr,i-1);
		ynew=Yorb(optr,i);
                if(xmode) xtest=xold*xnew;
                if(ymode) ytest=yold*ynew;

		     /* determine rest of coord's by linear interpolation */
                if (xtest<0) {
                    f = xold/(xnew-xold);
                    tplot = (1+f)*Torb(optr,i-1) - f*Torb(optr,i);
                    yplot = (1+f)*Yorb(optr,i-1) - f*Yorb(optr,i);
                    vplot = (1+f)*Vorb(optr,i-1) - f*Vorb(optr,i);
                    if (Uorb(optr,i-1)>0) {
                        yplot = -yplot;
                        vplot = -vplot;
                    }
                    if (pabs) yplot = ABS(yplot);
                    if (vabs) vplot = ABS(vplot);
                    printf ("%d %f %f %f\n",sos_count, tplot, yplot, vplot);
                    sos_count++;
                } 
                if (ytest<0) {
                    f = yold/(ynew-yold);
                    tplot = (1+f)*Torb(optr,i-1) - f*Torb(optr,i);
                    yplot = (1+f)*Xorb(optr,i-1) - f*Xorb(optr,i);
                    vplot = (1+f)*Uorb(optr,i-1) - f*Uorb(optr,i);
                    if (Vorb(optr,i-1)<0) {
                        yplot = -yplot;
                        vplot = -vplot;
                    }
                    if (pabs) yplot = ABS(yplot);
                    if (vabs) vplot = ABS(vplot);
                    printf ("%d %f %f %f\n",sos_count, tplot, yplot, vplot);
                    sos_count++;
                }
        }
	
}

/*
 * ORBSOS:	 tabulate an orbit as a surface of section
 *
 *	9-jun-88  PJT  Created - (still named mksos in a beta version)
 *     26-mar-92  V1.1 put table version into system - named orbsos.c   PJT
 *     24-may-92  V1.1a added <potential.h>                             PJT
 *     30-dec-93  V1.2 added pabs,vabs=t|f                              PJT
 *    9-dec-2019  V1.3 symm option to deal with e.g. hh64               PJT
 *    2-jan-2020  V1.4 process all orbits, not just the first           PJT
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
    "VERSION=1.4\n  	2-jan-2020 PJT",
    NULL,
};

string usage="Tabulate the SOS coordinates of an orbit";


string	infile;			/* file names */
stream  instr;			/* file streams */

orbitptr optr = NULL;
int xmode=0;                    // crossing X axis?
int ymode=0;                    // crossing Y axis?

local int lzs(real x, real y, real u, real v);
local void sos2(orbitptr optr, bool symm, bool pabs, bool vabs);

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
    if (xmode) dprintf(1,"SOS Y-VY will be computed\n");
    if (ymode) dprintf(1,"SOS X-VX will be computed\n");
    pabs = getbparam("pabs");
    vabs = getbparam("vabs");
    symm = getbparam("symm");
    
    while (read_orbit(instr,&optr))
      sos2(optr,symm,pabs,vabs);
    strclose(instr);
}


/* 	
 * compute a surface of section from a 2D orbit : classical method
 *
 */
 
static int sos_count=0;

int lzs(real x, real y, real u, real v)
{
  real lz = x*v - y*u;    // counter clock wise is +1
  return SGN(lz);
}

void sos2(orbitptr optr, bool symm, bool pabs, bool vabs)
{
  real xold,xnew, yold,ynew,f, xtest, ytest;
  real yplot,vplot,tplot;
  int  i, lz;
  int lz0 = lzs(Xorb(optr,0),Yorb(optr,0),Uorb(optr,0),Vorb(optr,0));
  
  dprintf(1,"Orbit with SGN(lz) = %d\n",lz0);
	
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
      if (Uorb(optr,i-1) > 0) {
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
      if (Vorb(optr,i-1) < 0) {
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

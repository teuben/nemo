/**********************************************************************/
/* PLOT.C: routine to plot snapshot.                                  */
/* Copyright (C) 2000 by Jin Koda, Tokyo, JAPAN                       */
/**********************************************************************/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "treecode.h"
#include "/usr/local/pgplot/cpgplot.h"

/*
 * Plot snapshot.
 */

void plot(bodyptr btab, int nbody)
{
    
    int i, pinc;
    bodyptr p;
    real prange, omgnow, xbar[2], ybar[2], tmp1, tmp2;
    vector xpos, xvel;

    pinc = 5;                                     /* ptcl incli. to plot  */
    prange = rinit / rscale * 1.5;                /* plot range           */
    omgnow = (1.02273e-9 * omgb) *(tnow * tscale);/* speed in sys.unit    */

    xbar[0] = 0.7 * prange * rcos(omgnow);        /* set bar direction    */
    ybar[0] = 0.7 * prange * rsin(omgnow);
    xbar[1] = 0.9 * prange * rcos(omgnow);
    ybar[1] = 0.9 * prange * rsin(omgnow);

    cpgbeg(0, "/xterm", 1, 1);                    /* begin pgplot         */
    cpgsubp(2,1);                                 /* display is two parts */

    cpgenv(-prange, prange, -prange, prange, 1, 0); 
    cpglab("x", "y", " ");                        /* set plot coordinate  */
    for (p = btab; p < btab + nbody; p+=pinc) {   /* plot particles       */
    	SETV(xpos, Pos(p));
	cpgpt(1, &xpos[0], &xpos[1], -2);
    }
    cpgsch(2.0);
    cpgarro(xbar[0], ybar[0], xbar[1], ybar[1]);
    MULVS(xbar, xbar, -1.0);
    MULVS(ybar, ybar, -1.0);
    cpgarro(xbar[0], ybar[0], xbar[1], ybar[1]);
    cpgsch(1.0);

    cpgenv(-prange, prange, -prange, prange, 1, 0);
    cpglab("x", "y", " ");
    cpgsch(0.3);
    for (p = btab; p < btab + nbody; p+=pinc) {
    SETV(xpos, Pos(p));
	SETV(xvel, xpos);
	ADDMULVS(xvel, Vel(p), 0.05);
	cpgarro(xpos[0], xpos[1], xvel[0], xvel[1]);
    }
    cpgsch(2.0);
    cpgarro(xbar[0], ybar[0], xbar[1], ybar[1]);
    MULVS(xbar, xbar, -1.0);
    MULVS(ybar, ybar, -1.0);
    cpgarro(xbar[0], ybar[0], xbar[1], ybar[1]);
    cpgsch(1.0);

    cpgend();                                     /* end pgplot           */
}

/*
 *  DIFFUSE:   diffuse orbits a bit, but conserve energy : XY plane only
 *
 *
 *       20-jan-04     special version for flowcode, where Acc = Vel
 *
 */

#include "defs.h"
 

#define ECONS  1		/* flag energy conservation */

diffuse (btab, nb, ndim, sigma)
Body *btab;
int   nb;
int   ndim;
real  sigma;
{
    real p, s, c, vx, vy;
    Body *b;
    bool Qspin = ome != 0.0;
    
    if (sigma<=0.0) return 1;          /* no work to do ... */
    dprintf(0,"diffuse %d\n",nb);

    for (b=btab; b<btab+nb; b++) {              /* loop over all bodies */
        p = grandom(0.0,sigma);         /* get random angle +/- sigma */
        s = sin(p);
        c = cos(p);
#ifdef SPIN
	if (Qspin) {                       /* correct for rot. frame */
	    Acc(b)[0] -= ome*Pos(b)[1];
            Acc(b)[1] += ome*Pos(b)[0];
	}
#endif
        vx = Acc(b)[0] * c  -  Acc(b)[1] * s;   /* rotate the vector */
        vy = Acc(b)[0] * s  +  Acc(b)[1] * c;
	if (b==btab+50) dprintf(1,"#50: %20.16g\n",
		dotvp(Acc(b),Acc(b))/(vx*vx+vy*vy));
        Acc(b)[0] = vx;
        Acc(b)[1] = vy;
#ifdef SPIN
        if (Qspin) {                       /* correct back to rot. frame */
            Acc(b)[0] += ome*Pos(b)[1];
            Acc(b)[1] -= ome*Pos(b)[0];
        }
#endif
    }

    return 1;                              /* success */
}	

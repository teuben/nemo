/*
 *  DIFFUSE:   diffuse orbits a bit, but conserve energy : XY plane only
 *
 *		Peter Teuben - june 1992
 *		20-jun-92   fiddled with -DSPIN ... ???
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
    double cos(), sin(), grandom(), sqrt();
    
    if (sigma==0.0) return 1;          /* no work to do ... */

    for (b=btab; b<btab+nb; b++) {              /* loop over all bodies */
        p = grandom(0.0,sigma);         /* get random angle +/- sigma */
        s = sin(p);
        c = cos(p);
#ifdef SPIN
	if (ome != 0.0) {                       /* correct for rot. frame */
	    Vel(b)[0] -= ome*Pos(b)[1];
            Vel(b)[1] += ome*Pos(b)[0];
	}
#endif
        vx = Vel(b)[0] * c  -  Vel(b)[1] * s;   /* rotate the vector */
        vy = Vel(b)[0] * s  +  Vel(b)[1] * c;
	if (b==btab+50) dprintf(1,"#50: %20.16g\n",
		dotvp(Vel(b),Vel(b))/(vx*vx+vy*vy));
        Vel(b)[0] = vx;
        Vel(b)[1] = vy;
#ifdef SPIN
        if (ome != 0.0) {                       /* correct back to rot. frame */
            Vel(b)[0] += ome*Pos(b)[1];
            Vel(b)[1] -= ome*Pos(b)[0];
        }
#endif
    }

    return 1;                              /* success */
}	

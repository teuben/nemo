/*
 *  DIFFUSE:   diffuse orbits a bit, but conserve energy : XY plane only
 *
 *             this routine can optionally just compute the angle
 *             and diffuse later with the same diffusion angle
 *
 *
 *       20-jan-04     special version for flowcode, where Acc = Vel
 *        6-feb-04     store deflection sigma in Aux(), optional rotation here
 *
 */

#include "defs.h"


/* 
 * ROTATE_AUX:  rotate velocities (stored in Acc()) by angle defined in Aux()
 *              and stuff them back into Acc() for later retrieval
 */

void rotate_aux(bodyptr btab, int nb)
{
  real c,s, vx,vy;
  bodyptr p;
  static int count;

  count += nb;
  dprintf(1,"rotate_aux: %d %d\n",nb,count);

  for (p = btab; p < btab+nb; p++) {		/* loop over bodies */
    dprintf(1," angle %g\n",Aux(p));
    s = sin(Aux(p));
    c = cos(Aux(p));
    
    vx = Acc(p)[0] * c  -  Acc(p)[1] * s;
    vy = Acc(p)[0] * s  +  Acc(p)[1] * c;
    
    Acc(p)[0] = vx;
    Acc(p)[1] = vy;
  }
}


/* 
 *  DIFFUSE:   find, and optionally apply, a diffusion angle
 */


void diffuse (Body *btab, int   nb, int   ndim, real  sigma, bool Qrotate)
{
    Body *b;
    
    if (sigma<=0.0) return;          /* no work to do ... */
    dprintf(2,"diffuse %d %g\n",nb,sigma);

    for (b=btab; b<btab+nb; b++) {              /* loop over all bodies */
	Aux(b) = grandom(0.0,sigma);         /* get random angle +/- sigma */
    }
    if (Qrotate) rotate_aux(btab, nb);
}	

 

void diffuse_old (Body *btab, int   nb, int   ndim, real  sigma, bool Qrotate)
{
    real p, s, c, vx, vy;
    Body *b;
#ifdef SPIN
    bool Qspin = ome != 0.0;
#endif
    
    if (sigma<=0.0) return;          /* no work to do ... */
    dprintf(2,"diffuse %d %g\n",nb,sigma);

    for (b=btab; b<btab+nb; b++) {              /* loop over all bodies */
        p = grandom(0.0,sigma);         /* get random angle +/- sigma */
	Aux(b) = p;                     /* save it in Aux() for potential re-use */
        s = sin(p);
        c = cos(p);
	if (!Qrotate) continue;
#ifdef SPIN
	if (Qspin) {                       /* correct for rot. frame */
	    Acc(b)[0] -= ome*Pos(b)[1];
            Acc(b)[1] += ome*Pos(b)[0];
	}
#endif
        vx = Acc(b)[0] * c  -  Acc(b)[1] * s;   /* rotate the vector */
        vy = Acc(b)[0] * s  +  Acc(b)[1] * c;
	if (b==btab+5) dprintf(1,"#5: %20.16g\n",
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
}	

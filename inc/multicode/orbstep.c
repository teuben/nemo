/*
 * ORBSTEP.C: general-purpose orbit-integration routines.
 * Defines: initstep(), orbstep().
 * Requires: MBODY, Body, Pos(), Vel(), Acc().
 */

#include "quaddefs.h"

/*
 * ABAK0, ..., ABAK3: saved accelerations, latest to oldest.
 */

#define MCOR  (NDIM * MBODY)

local real abak0[MCOR];
local real abak1[MCOR];
local real abak2[MCOR];
local real abak3[MCOR];

/*
 * INITSTEP: initialize the orbit-integrator.
 */

local int nstep;	/* integration step counter */

initstep(btab, nb, tptr, force)
Body *btab;		/* array of bodies */
int nb;			/* number of bodies */
real *tptr;		/* initial time */
force_proc force;	/* acceleration calculation */
{
    nstep = 0;					/* start counting steps     */
    (*force)(btab, nb, *tptr);			/* compute (t-dep) force    */
    moveaccel(btab, nb);			/* save resulting accel     */
}

/*
 * ORBSTEP: integrate the orbit of a set of bodies.
 */

orbstep(btab, nb, tptr, force, dt, mode)
Body *btab;		/* array of bodies */
int nb;			/* number of bodies */
real *tptr;		/* current time */
force_proc force;	/* acceleration calculation */
real dt;		/* integration time step */
int mode;		/* select integration algorithm */
{
    if (mode == 1 || nstep < 3) {		/* RK algorithm required?   */
        rkstep(btab, nb, tptr, force, dt, abak3);
						/*   take RK time-step      */
	(*force)(btab, nb, *tptr);		/*   compute final force    */
    } else if (mode == 2 || mode == 3) {	/* PC algorithm possible?   */
        pcstep(btab, nb, tptr, force, dt);	/*   take PC time-step      */
	if (mode == 2)				/*   final force requested? */
	    (*force)(btab, nb, *tptr);		/*   compute final force    */
    } else
	error("orbstep: unknown mode %d\n", mode);
    moveaccel(btab, nb);			/* save final accel         */
    nstep++;					/* count another time-step  */
}

rkstep(btab, nb, tptr, force, dt, atmp1)
Body *btab;		/* array of bodies */
int nb;			/* number of bodies */
real *tptr;		/* current time */
force_proc force;	/* acceleration calculation */
real dt;		/* integration time step */
real atmp1[];		/* scratch accelerations */
{
    Body *p;
    int i, k;
    register real *pptr, *vptr, *aptr;
    real dt2, dts4, dt6, dts6;
    static real atmp2[MCOR];			/* dont alloc on stack */

    dt2 = dt / 2;
    for (p = btab; p < btab+nb; p++) {		/* loop over bodies */
        pptr = Pos(p);				/*   get body coords */
	vptr = Vel(p);
	for (k = 0; k < NDIM; k++)		/*   loop over coords */
	    *pptr++ += dt2 * (*vptr++);		/*     set position x_1 */
    }
    (*force)(btab, nb, *tptr + dt2);		/* get accel a_1 */
    dts4 = dt2 * dt2;
    for (p = btab, i = 0; p < btab+nb; p++) {	/* loop over bodies */
        pptr = Pos(p);				/*   get body coords */
	aptr = Acc(p);
	for (k = 0; k < NDIM; k++, i++) {	/*   loop over coords */
	    *pptr++ += dts4 * abak0[i];		/*     set position x_2 */
	    atmp1[i] = *aptr++;			/*     save accel a_1 */
	}
    }
    (*force)(btab, nb, *tptr + dt2);		/* get accel a_2 */
    for (p = btab, i = 0; p < btab+nb; p++) {	/* loop over bodies */
        pptr = Pos(p);				/*   get body coords */
	vptr = Vel(p);
	aptr = Acc(p);
	for (k = 0; k < NDIM; k++, i++) {	/*   loop over coords */
	    *pptr++ += dt2 * (*vptr++) + dts4 * (2*atmp1[i] - abak0[i]);
	    					/*     set position x_3 */
	    atmp2[i] = *aptr++;			/*     save accel a_2 */
	}
    }
    (*force)(btab, nb, *tptr + dt);		/* get accel a_3 */
    dt6 = dt / 6;
    dts6 = dt * dt6;
    for (p = btab, i = 0; p < btab+nb; p++) {	/* loop over bodies */
        pptr = Pos(p);				/*   get body coords */
	vptr = Vel(p);
	aptr = Acc(p);
	for (k = 0; k < NDIM; k++, i++) {	/*   loop over coords */
	    *vptr++ += dt6 * (abak0[i] + 2*atmp1[i] + 2*atmp2[i] + *aptr++);
	    *pptr++ += dts6 * (abak0[i] - 2*atmp1[i] + atmp2[i]);
	}
    }
    *tptr += dt;
}

pcstep(
       Body *btab,		/* array of bodies */
       int nb,			/* number of bodies */
       real *tptr,		/* current time */
       force_proc force,	/* acceleration calculation */
       real dt)		        /* integration time step */
{
    Body *p;
    int i, k;
    register real *pptr, *vptr, *aptr;
    real dt360, dts32, dt720, app, acp, acv;

    dt360 = dt / 360;
    for (p = btab, i = 0; p < btab+nb; p++) {	/* loop over bodies */
        pptr = Pos(p);				/*   get body coords */
	vptr = Vel(p);
	for (k = 0; k < NDIM; k++, i++) {	/*   loop over coords */
	    app = 323 * abak0[i] - 264 * abak1[i] +
		    159 * abak2[i] -  38 * abak3[i];
	    *pptr++ += dt * (*vptr++ + dt360 * app);
						/*     predict new position */
	}
    }
    (*force)(btab, nb, *tptr + dt);		/* find force at pred. pos. */
    dts32 = dt*dt / 32;
    dt720 = dt / 720;
    for (p = btab, i = 0; p < btab+nb; p++) {	/* loop over bodies */
        pptr = Pos(p);				/*   get body coords */
	vptr = Vel(p);
	aptr = Acc(p);
	for (k = 0; k < NDIM; k++, i++) {	/*   loop over coords */
	    acp = 3 * (*aptr) - 12 * abak0[i] +
		    18 * abak1[i] - 12 * abak2[i] + 3 * abak3[i];
	    *pptr++ += dts32 * acp;		/*     correct position */
	    acv = 251 * (*aptr++) + 646 * abak0[i] -
		    264 * abak1[i] + 106 * abak2[i] - 19 * abak3[i];
	    *vptr++ += dt720 * acv;		/*     advance velocity */
	}
    }
    *tptr += dt;				/* advance time */
}


moveaccel(Body *btab, int nb)
{
    Body *p;
    int i, k;
    register real *aptr;

    for (p = btab, i = 0; p < btab+nb; p++) {	/* loop over bodies */
        aptr = Acc(p);				/*   get body accel */
	for (k = 0; k < NDIM; k++, i++) {	/*   loop over coords */
	    abak3[i] = abak2[i];		/*     copy history back */
	    abak2[i] = abak1[i];
	    abak1[i] = abak0[i];
	    abak0[i] = *aptr++;			/*     copy latest accel */
	}
    }
}

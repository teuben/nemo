
/*
 * ORBSTEP.C: general-purpose orbit-integration routines for flowcode
 *            
 *
 * Defines: initstep(), orbstep().
 * Requires: MBODY, body, bodyptr, Pos(), Vel(), Acc().
 *
 *  10-jun-92  Added the 'rk4' method, but this
 *             now uses VECTMATH and assumes particles are
 *             not 'interacting' - build for rotating potentials  
 *             Note that RK, PC and PC1 don't work for rotating
 *             potentials.					PJT
 *  11-apr-96  adapted for flowcode from the potcode version    PJT
 *   6-feb-04  overhauled the code and defined diffusion angles in Aux()  PJT
 *
 */

#include "defs.h"

/*
 * ABAK0, ..., ABAK3: saved accelerations, latest to oldest.
 */

#define MCOR  (NDIM * MBODY)

local real abak0[MCOR];
local real abak1[MCOR];
local real abak2[MCOR];
local real abak3[MCOR];


/*
 * RKSTEP: (some kind of) Runge-Kutta step - mode=1
 */
static void rkstep(btab, nb, tptr, force, dt, atmp1)
bodyptr btab;		/* array of bodies */
int nb;			/* number of bodies */
real *tptr;		/* current time */
fproc force;		/* acceleration calculation */
real dt;		/* integration time step */
real atmp1[];		/* scratch accelerations */
{
    bodyptr p;
    int i, k;
    register real *pptr, *vptr, *aptr;
    real dt2, dts4, dt6, dts6;
    static real atmp2[MCOR];			/* dont alloc on stack */

    dt2 = dt / 2;
    for (p = btab; p < btab+nb; p++) {		/* loop over bodies */
        pptr = Pos(p);				/*   get body coords */
	vptr = Acc(p);
	for (k = 0; k < NDIM; k++)		/*   loop over coords */
	    *pptr++ += dt2 * (*vptr++);		/*     set position x_1 */
    }
    (*force)(btab, nb, *tptr + dt2, FALSE);	/* get accel a_1 */
    dts4 = dt2 * dt2;
    for (p = btab, i = 0; p < btab+nb; p++) {	/* loop over bodies */
        pptr = Pos(p);				/*   get body coords */
	aptr = Acc(p);
	for (k = 0; k < NDIM; k++, i++) {	/*   loop over coords */
	    *pptr++ += dts4 * abak0[i];		/*     set position x_2 */
	    atmp1[i] = *aptr++;			/*     save accel a_1 */
	}
    }
    (*force)(btab, nb, *tptr + dt2, FALSE);	/* get accel a_2 */
    for (p = btab, i = 0; p < btab+nb; p++) {	/* loop over bodies */
        pptr = Pos(p);				/*   get body coords */
	vptr = Acc(p);
	aptr = Acc(p);
	for (k = 0; k < NDIM; k++, i++) {	/*   loop over coords */
	    *pptr++ += dt2 * (*vptr++);         /*     set position x_3 */
	    atmp2[i] = *aptr++;			/*     save accel a_2 */
	}
    }
    (*force)(btab, nb, *tptr + dt, TRUE);	/* get accel a_3 */
    dt6 = dt / 6;
    dts6 = dt * dt6;
#if 0
    for (p = btab, i = 0; p < btab+nb; p++) {	/* loop over bodies */
        pptr = Pos(p);				/*   get body coords */
	vptr = Acc(p);
	aptr = Acc(p);
	for (k = 0; k < NDIM; k++, i++) {	/*   loop over coords */
	    *pptr++ += dts6 * (abak0[i] - 2*atmp1[i] + atmp2[i]);
	}
    }
#endif
    *tptr += dt;
}

/* 
 * PCSTEP: Predictor-Corrector step  (mode=2,3)
 */
static void pcstep(btab, nb, tptr, force, dt)
bodyptr btab;		/* array of bodies */
int nb;			/* number of bodies */
real *tptr;		/* current time */
fproc force;		/* acceleration calculation */
real dt;		/* integration time step */
{
    bodyptr p;
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
    (*force)(btab, nb, *tptr + dt,TRUE);	/* find force at pred. pos. */
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

/*
 *   EULER: a simple first order integration method (mode=0)
 */
static void eulerstep(btab, nb, tptr, force, dt)
bodyptr btab;		/* array of bodies */
int nb;			/* number of bodies */
real *tptr;		/* current time */
fproc force;		/* acceleration calculation */
real dt;		/* integration time step */
{
    bodyptr p;
    int k;
    register real *pptr, *vptr, *aptr;

    for (p=btab; p<btab+nb; p++) {
        pptr = Pos(p);
        vptr = Vel(p);
        aptr = Acc(p);
        for (k=0; k<NDIM; k++) {
            *pptr++ += dt * (*aptr);
            *vptr++ = *aptr++;
        }
    }
    *tptr += dt;
}
/*
 *   LEAPFROG:   mode=5
 */

static void leapfrogstep(btab, nb, tptr, force, dt)
bodyptr btab;		/* array of bodies */
int nb;			/* number of bodies */
real *tptr;		/* current time */
fproc force;		/* acceleration calculation */
real dt;		/* integration time step */
{
    error("leapfrog stepping not implemented");
}


/*
 * RK4: mode=4
 */
static void rk4step(btab, nb, tptr, force, dt)
bodyptr btab;		/* array of bodies */
int nb;			/* number of bodies */
real *tptr;		/* current time */
fproc force;		/* acceleration calculation */
real dt;		/* integration time step */
{
    body tmp1;
    bodyptr p, pt=&tmp1;
    int k;
    vector dpos1, dpos2;
    real dth, dt6;

    dth = 0.5*dt;
    dt6 = dt/6.0;

    (*force)(btab,nb,*tptr,TRUE);	/* needed only when fake physics */

    for (p=btab; p<btab+nb; p++) {

        dprintf(1,"RK4/0: %g %g %g\n",Aux(p),Acc(p)[0],Acc(p)[1]);
        for (k=0; k<NDIM; k++) {
            Pos(pt)[k] = Pos(p)[k] + dth*Acc(p)[k];
        }
	Aux(pt) = Aux(p);   /* update sigma for this particle */

        (*force)(pt,1,*tptr+dth,FALSE);
	dprintf(1,"RK4/1: %g %g %g\n",Aux(pt),Acc(pt)[0],Acc(pt)[1]);

        SETV(dpos1, Acc(pt));                           /* save 'dyt' */
        for (k=0; k<NDIM; k++) {
            Pos(pt)[k] = Pos(p)[k] + dth*Acc(pt)[k];
        }

        (*force)(pt,1,*tptr+dth,FALSE);
	dprintf(1,"RK4/2: %g %g %g\n",Aux(pt),Acc(pt)[0],Acc(pt)[1]);
        
        SETV(dpos2, Acc(pt));                           /* save 'dym' */
        for (k=0; k<NDIM; k++) {
            Pos(pt)[k] = Pos(p)[k] + dt*Acc(pt)[k];
        }

        ADDV(dpos2, dpos2, dpos1);

        (*force)(pt,1,*tptr+dt,FALSE);
	dprintf(1,"RK4/3: %g %g %g\n",Aux(pt),Acc(pt)[0],Acc(pt)[1]);
        
        for (k=0; k<NDIM; k++) {
            Pos(p)[k] += dt6*(Acc(p)[k]+Acc(pt)[k]+2*dpos2[k]);
        }
	
	/* and compute the correct velocity in case somebody needs it */
	(*force)(p,1,*tptr+dt,FALSE);
	dprintf(1,"RK4/4: %g %g %g\n",Aux(p),Acc(p)[0],Acc(p)[1]);
	SETV(Vel(p),Acc(p));
    }
    *tptr += dt;
}


/* 
 * MOVEACCEL: local utility to stack back old values of the
 *	      forces. Only used in RK, PC and PC1.
 *
 *  I/O bodyptr btab;		array of bodies
 *  I   int nb;			number of bodies
 */
local void moveaccel(bodyptr btab, int nb)
{
    bodyptr p;
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


/*
 * INITSTEP: initialize the orbit-integrator.
 */

local int nstep;	/* integration step counter */

void initstep(btab, nb, tptr, force)
bodyptr btab;		/* array of bodies */
int nb;			/* number of bodies */
real *tptr;		/* initial time */
fproc force;		/* acceleration calculation */
{
    nstep = 0;					/* start counting steps */
    (*force)(btab, nb, *tptr, TRUE);		/* compute (t-dep) force */
    moveaccel(btab, nb);			/* save resulting accel */
}
 


/*
 * ORBSTEP: integrate the orbit of a set of bodies.
 */

void orbstep(btab, nb, tptr, force, dt, mode)
bodyptr btab;		/* array of bodies */
int nb;			/* number of bodies */
real *tptr;		/* current time */
fproc force;		/* acceleration calculation */
real dt;		/* integration time step */
int mode;		/* select integration algorithm */
{
  if (mode == 0) {                              /* simplest Eulerian */
    eulerstep(btab, nb, tptr, force, dt);       /* take Euler step */
    (*force)(btab, nb, *tptr,TRUE);             /* compute new force */
  } else if (mode == 5) {			/* use leapfrog ? */
    leapfrogstep(btab, nb, tptr, force, dt);    /* take step */
    (*force)(btab, nb, *tptr,TRUE);             /* compute final force ? */
  } else if (mode == 4) {			/* use 4th order RK ? */
    rk4step(btab, nb, tptr, force, dt);	        /* take RK4 step */
    (*force)(btab, nb, *tptr, TRUE);	        /* compute final force */
  } else if (mode == 1 || nstep < 3) { 	        /* RK algorithm required? */
    rkstep(btab, nb, tptr, force, dt, abak3);   /*   take RK time-step */
    (*force)(btab, nb, *tptr, TRUE);  	        /*   compute final force */
  } else if (mode == 2 || mode == 3) {	        /* PC algorithm possible? */
    pcstep(btab, nb, tptr, force, dt);	        /*   take PC time-step */
    if (mode == 2)				/*   final force requested? */
      (*force)(btab, nb, *tptr,TRUE);		/*   compute final force */
  } else 
    error("orbstep: unknown mode %d", mode);
  moveaccel(btab, nb);			        /* save final accel */
  nstep++;					/* count another time-step */
}

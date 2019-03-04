/* NBODY00.C: handedited f2c translated version of NBODY0.F 
 *
 * 	S.J. Aarseth's Standard N-Body Program (nbody0) 
 *
 * 	Adapted from the source code in Binney & Tremaine's (1987)  
 *       book 'Galactic Dynamics',  
 *         by Peter Teuben - June '89:  
 *       - all variables to be declared and a choice of 'real' or  
 *         'double precision'  
 * 	- all I/O is done through subroutines in order to allow easy  
 * 	  interface with development different systems, e.g. NEMO89.  
 * 	- PARAMETERS for NMAX and NDIM via an include file (see Makefile)  
 * 	  Note that NDIM should NEVER be reset from 3, the code does  
 * 	  not handle 2D stuff (yet).  
 *
 *          jun-89   Created  
 * 	 23-jan-90   Back to double precision  - PJT  
 * 	  5-apr-90   started to declare all variables - 	PJT  
 *       14-nov-91   finished(!) declaring all variables;  
 *                   tested using IMPLICIT NONE and compile with -u  
 *       15-nov-91   Connection Machine testing in nbody0.fcm file  
 *	 10-may-92   f2c conversion - manually optimized
 *                   body index (i,j) are 0 based
 *                   dim index (k) is 0 based
 *	 11-feb-98   Fixed array index for f2dot(k)  a(16)->a(15)
 *	  6-jan-00   new F77_FUNC macros, renamed to nbody00.c since
 *		     compiled name could clash with nbody0.f
 *       21-jan-00   added debug output and reset= option
 *       21-feb-04   stop if f2dot is 0
 *       24-feb-04   compute steps from fdot/f3dot if f2dot=0
 *       13-mar-04   nreset was not initialized
 *       20-apr-04   was never using the header file. now it does.
 */

#include "nbody0.h"

/* extract X, Y and Z components of a 3-vector */
#define X(a,i)	(a[i*NDIM])
#define Y(a,i)  (a[i*NDIM+1])
#define Z(a,i)  (a[i*NDIM+2])

#define SQR(x)  ((x)*(x))

/* Table of constant values */

static int c_nmax = NMAX;

#include <nemo.h>
#include "proto.h"

/* *********************************************************************** */
void nbody0(void)
{
    /* Initialized data */
    double time = 0.0, tnext = 0.0, e0;
    int nsteps = 0, makesure, nout = 0, nreset=0;
    double time0, time1=0.0;        /* counters for my timing/debugging */

    /* System (f2c) generated locals */
    double d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    double fdot[NDIM*NMAX], body[NMAX], step[NMAX], 
           f1dot[NDIM], f2dot[NDIM], f3dot[NDIM], 
           x0dot[NDIM*NMAX], a[17], f[NDIM*NMAX],
           e, s, x[NDIM*NMAX], tcrit, 
           d1[NDIM*NMAX], d2[NDIM*NMAX], d3[NDIM*NMAX], f1[NDIM], 
	   t0[NMAX], t1[NMAX], t2[NMAX], t3[NMAX], x0[NDIM*NMAX], 
	   dt, deltat, dt1, dt2, dt3, eta, eps2, t1pr, t2pr, t3pr;
    double sum_f, sum_fdot, sum_f2dot, sum_f3dot;
    int    i, j, k, n, itmp, reset, use3dot;


    inpars(&c_nmax, &n, &eta, &deltat, &tcrit, &eps2, &reset, &use3dot);
    inbods(&n, body, x0, x0dot);

    /* obtain total forces and first derivative for each body */
    for (i = 0; i < n; ++i) {
	for (k = 0; k < NDIM; ++k) {
	    f[k+i*3] = 0.0;
	    fdot[k+i*3] = 0.0;
	    d1[k+i*3] = 0.0;
	    d2[k+i*3] = 0.0;
	    d3[k+i*3] = 0.0;
	}
	for (j = 0; j < n; ++j) {
	    if (j == i) continue;
	    for (k = 0; k < NDIM; ++k) {
		a[k    ] = x0[k+j*3] - x0[k+i*3];
		a[k + 3] = x0dot[k+j*3] - x0dot[k+i*3];
	    }
	    d__1 = a[0];
	    d__2 = a[1];
	    d__3 = a[2];
	    a[6] = 1.0 / (SQR(d__1) + SQR(d__2) + SQR(d__3) + eps2);
	    a[7] = body[j] * a[6] * sqrt(a[6]);
	    a[8] = (a[0]*a[3] + a[1]*a[4] + a[2]*a[5]) * 3.0 * a[6];
	    for (k = 0; k < NDIM; ++k) {
		f[k+i*3] += a[k] * a[7];
		fdot[k+i*3] += (a[k + 3] - a[k] * a[8]) * a[7];
	    }
	}
    } /* for(i) */
    /*  form second and third derivative */
    for (i = 0; i < n; ++i) {
	for (j = 0; j < n; ++j) {
	    if (i == j) continue;
	    for (k = 0; k < NDIM; ++k) {
		a[k    ] = x0[k+j*3] - x0[k+i*3];
		a[k + 3] = x0dot[k+j*3] - x0dot[k+i*3];
		a[k + 6] = f[k+j*3] - f[k+i*3];
		a[k + 9] = fdot[k+j*3] - fdot[k+i*3];
	    }
	    d__1 = a[0];
	    d__2 = a[1];
	    d__3 = a[2];
	    a[12] = 1.0 / (SQR(d__1) + SQR(d__2) + SQR(d__3) + eps2);
	    a[13] = body[j] * a[12] * sqrt(a[12]);
	    a[14] = (a[0]*a[3] + a[1]*a[4] + a[2]*a[5]) * a[12];
	    d__1 = a[3];
	    d__2 = a[4];
	    d__3 = a[5];
	    d__4 = a[14];
	    a[15] = (SQR(d__1) + SQR(d__2) + SQR(d__3) + 
	             a[0]*a[6] + a[1]*a[7] + a[2]*a[8]) * a[12] + SQR(d__4);
	    d__1 = a[14];
	    a[16] = ((a[3]*a[6] + a[4]*a[7] + a[5]*a[8]) * 9 +
		     (a[0]*a[9] + a[1]*a[10] + a[2]*a[11]) * 3) * 
		    a[12] + a[14] * (a[15] * 9.0 - SQR(d__1) * 12.);
	    for (k = 0; k < NDIM; ++k) {
		f1dot[k] = a[k + 3] - a[14] * 3 * a[k];
		f2dot[k] = (a[k + 6] - a[14] * 6 * f1dot[k] - a[15] * 
			3 * a[k]) * a[13];
		f3dot[k] = (a[k + 9] - a[15] * 9 * f1dot[k] - a[16] * 
			a[k]) * a[13];
		d2[k+i*3] += f2dot[k];
		d3[k+i*3] = d3[k+i*3] + f3dot[k] - a[14] * 9 * f2dot[k];
		d1[k+i*3] += f3dot[k];  /* PJT */
	    }
	}
    } /* for (i) */

    /* initialize integration steps and convert to force diffences          */
    /* STEP = sqrt(ETA * (F*F2DOT + FDOT*FDOT)/(F2DOT*F2DOT + FDOT*F3DOT))  */
    for (i = 0; i < n; ++i) {
	d__1 = X(f,i);
	d__2 = Y(f,i);
	d__3 = Z(f,i);
	d__4 = X(d2,i);
	d__5 = Y(d2,i);
	d__6 = Z(d2,i);
	sum_f     = SQR(d__1)+SQR(d__2)+SQR(d__3);
	sum_f2dot = SQR(d__4)+SQR(d__5)+SQR(d__6);
	d__1 = X(fdot,i);
	d__2 = Y(fdot,i);
	d__3 = Z(fdot,i);
	d__4 = X(d1,i);
	d__5 = Y(d1,i);
	d__6 = Z(d1,i);
	sum_fdot  = SQR(d__1)+SQR(d__2)+SQR(d__3);
	sum_f3dot = SQR(d__4)+SQR(d__5)+SQR(d__6);	

	if (sum_f2dot == 0.0) {
	  warning("You have hit upon an F2DOT=0 case, correcting....");
	  step[i] = sqrt(eta * sqrt(sum_fdot/sum_f3dot));	  
	  dprintf(1,"INT0: i=%d time=%g sum_fdot=%g sum_f3dot=%g -> %g\n",
		  i,time,sum_fdot,sum_f3dot,step[i]);
	} else {
	  step[i] = sqrt(eta * sqrt(sum_f/sum_f2dot));
	  dprintf(1,"INT0: i=%d time=%g sum_f=%g sum_f2dot=%g -> %g \n",
		  i,time,sum_f,sum_f2dot,step[i]);
	  dprintf(1,"INT0: i=%d time=%g sum_fdot=%g sum_f3dot=%g -> %g [***]\n",
		  i,time,sum_fdot,sum_f3dot,
		  sqrt(eta*(sqrt(sum_f)*sqrt(sum_f2dot)+sum_fdot)/
		       (sum_f2dot + sqrt(sum_fdot)*sqrt(sum_f2dot))));
	}
	t0[i] = time;
	t1[i] = time - step[i];
	t2[i] = time - step[i] * 2;
	t3[i] = time - step[i] * 3;
	for (k = 0; k < NDIM; ++k) {
	    d1[k+i*3] = (d3[k+i*3] * step[i] / 6 - 
	                         d2[k+i*3] * 0.5) * step[i] + 
	                         fdot[k+i*3];
	    d2[k+i*3] = d2[k+i*3] * 0.5 - d3[k+i*3] * 0.5 * step[i];
	    d3[k+i*3] /= 6;
	    f[k+i*3] *= 0.5;
	    fdot[k+i*3] /= 6;
	}
    } /* for (i) */

/* ----------------------------------------------------------------------- */
/* ==============> entry point when major output + diagnostics to be done */

/*              energy check and output */
L100:
    e = 0.0;
    for (i = 0; i < n; ++i) {
	dt = tnext - t0[i];
	for (k = 0; k < NDIM; ++k) {
	    f2dot[k] = d3[k+i*3] * (t0[i] - t1[i] + 
	        (t0[i] - t2[i])) + d2[k+i*3];
	    x[k+i*3] = ((((d3[k+i*3] * 0.05 * dt + f2dot[k] / 12.0) * dt  +
		    fdot[k+i*3]) * dt + f[k+i*3]) * dt + 
		    x0dot[k+i*3]) * dt + x0[k+i*3];
	    a[k] = (((d3[k+i*3] * 0.25 * dt + f2dot[k] /
		     3.0) * dt + fdot[k+i*3] * 3.0) * dt 
		    + f[k+i*3] * 2.0) * dt + x0dot[k+i*3];
	}
	itmp=i+1;
	//  this suffers from the Heisenbug, output depends on the timestep
	outbods(&body[i], &x[i*3], a, &step[i], &itmp);
	d__1 = a[0];
	d__2 = a[1];
	d__3 = a[2];
	e += 0.5 * body[i] * (SQR(d__1) + SQR(d__2) + SQR(d__3));
    }
    for (i = 0; i < n; ++i) {
	for (j = 0; j < n; ++j) {
	    if (j == i) continue;
	    d__1 = X(x,i) - X(x,j);
	    d__2 = Y(x,i) - Y(x,j);
	    d__3 = Z(x,i) - Z(x,j);
	    e -= body[i] * 0.5 * body[j] / 
		    sqrt(SQR(d__1) + SQR(d__2) + SQR(d__3) + eps2);
	}
    }
    if (nout ==0) e0 = e;
    nout++;
    outene(&tnext, &nsteps, &e);
    if (time > tcrit) {
        dprintf(0,"Time spent in searching for next advancement: %g\n",
                time1*60.0);
        dprintf(0,"Energy conservation: %g / %g = %g\n",e-e0,e0,(e-e0)/e);
	if (reset) 
            dprintf(0,"Time resets needed %d times / %d dumps\n",nreset,nout);
        return;
    } 
    tnext += deltat;
    makesure = 1;


/* ----------------------------------------------------------------------- */
/*   ===============>  Normal entry point for next timstep */


/*      find next body (i) to be advanced and set new time 
 *      this linear search takes about 5% of the total CPU
 *      see sift() to speed up this process 
 */
L200:
    time0=cputime();
    time = 1e10;
    for (j = 0; j < n; ++j) {
	if (time > t0[j] + step[j]) {
	    i = j;
	    time = t0[j] + step[j];
	}
    }
    time1 += (cputime()-time0);

    /* however, the time needs to be reset if the next datadump is
     * coming up earlier than the next particle. Why? Because
     * if we don't this next particle will then be integrated beyond
     * tnext ??  When Sverre and I were discussing this in January 2000,
     * two solutions were proposed:
     */
#if 1
                                            /* Peter : only do it after dump */
    dprintf(3,"time=%g i=%d\n",time,i);
    if (makesure) {
        dprintf(2,"TIME?: %g %g\n",time,tnext);
        if (reset && (time > tnext)) {
            dprintf(1,"RESET: %g > %g\n",time,tnext);
	    time = tnext;
	    nreset++;
        }
        makesure=0;
    }
#else
                                            /* Sverre : do it always */
	/* bug: integrates too far */
    if (time > tnext) {
            dprintf(1,"RESET: %g > %g\n",time,tnext);
	    time = tnext;
	    nreset++;
    }
#endif
    
    /*  predict all coordinates to first order in force derivative */
    for (j = 0; j < n; ++j) {
	s = time - t0[j];
	X(x,j) = ((X(fdot,j) * s + X(f,j)) * s + X(x0dot,j)) * s + X(x0,j);
	Y(x,j) = ((Y(fdot,j) * s + Y(f,j)) * s + Y(x0dot,j)) * s + Y(x0,j);
	Z(x,j) = ((Z(fdot,j) * s + Z(f,j)) * s + Z(x0dot,j)) * s + Z(x0,j);
    }
    /*  include 2nd and 3rd order and obtain the velocity */
    dt = time - t0[i];
    for (k = 0; k < NDIM; ++k) {
	f2dot[k] = d3[k+i*3] * (t0[i] - t1[i] + (t0[i]
		 - t2[i])) + d2[k+i*3];
	d__1 = dt, d__1 *= d__1;
	x[k+i*3] = (d3[k+i*3] * 0.05 * dt + f2dot[k]
		 / 12.) * SQR(d__1) + x[k+i*3];
	x0dot[k+i*3] = (((d3[k+i*3] * 0.25 * dt + f2dot[k] / 3.0) * 
	        dt + fdot[k+i*3] * 3.0) * 
		dt + f[k+i*3] * 2.0) * dt + x0dot[k+i*3];
	f1[k] = 0.0;
    }
    /*  obtain current force on i-th body */
    for (j = 0; j < n; ++j) {
	if (j == i) continue;
	d__1 = a[0] = X(x,j) - X(x,i);
	d__2 = a[1] = Y(x,j) - Y(x,i);
	d__3 = a[2] = Z(x,j) - Z(x,i);
	a[3] = 1.0 / (SQR(d__1) + SQR(d__2) + SQR(d__3) + eps2);
	a[4] = body[j] * a[3] * sqrt(a[3]);
	f1[0] += a[0] * a[4];
	f1[1] += a[1] * a[4];
	f1[2] += a[2] * a[4];
    }

    /*  set time intervals for new difference and update the times */
    dt1 = time - t1[i];
    dt2 = time - t2[i];
    dt3 = time - t3[i];
    t1pr = t0[i] - t1[i];
    t2pr = t0[i] - t2[i];
    t3pr = t0[i] - t3[i];
    t3[i] = t2[i];
    t2[i] = t1[i];
    t1[i] = t0[i];
    t0[i] = time;
    /*  form new differences and include 4th order semi-iterative  */
    for (k = 0; k < NDIM; ++k) {
	a[k    ] = (f1[k   ] -  f[k+i*3 ] * 2.0) / dt;
	a[k + 3] = (a[k    ] - d1[k+i*3]) / dt1;
	a[k + 6] = (a[k + 3] - d2[k+i*3]) / dt2;
	a[k + 9] = (a[k + 6] - d3[k+i*3]) / dt3;
	d1[k+i*3] = a[k    ];
	d2[k+i*3] = a[k + 3];
	d3[k+i*3] = a[k + 6];
	f1dot[k] = t1pr * t2pr * t3pr * a[k + 9];
	f2dot[k] = (t1pr * t2pr + t3pr * (t1pr + t2pr)) * a[k + 9];
	f3dot[k] = (t1pr + t2pr + t3pr) * a[k + 9];
	d__1 = dt, d__2 = d__1;
	x0[k+i*3] = (((a[k + 9] * dt / 30. + f3dot[k] * 
		0.05) * dt + f2dot[k] / 12.) * dt + f1dot[k] / 6.) 
		* (d__2 * SQR(d__1)) + x[k+i*3];
	d__1 = dt;
	x0dot[k+i*3] = (((a[k + 9] * 0.2 * dt + f3dot[k] * 
		0.25) * dt + f2dot[k] / 3.0) * dt + f1dot[k] * 0.5) * 
		SQR(d__1) + x0dot[k+i*3];
    }
    /*  scale F and FDOT by factorials and set new integration step  */
    for (k = 0; k < NDIM; ++k) {
	f[k+i*3] = f1[k] * 0.5;
	fdot[k+i*3] = ((d3[k+i*3] * dt1 + d2[k+i*3]) *
		 dt + d1[k+i*3]) / 6.0;
	f2dot[k] = (d3[k+i*3] * (dt + dt1) + d2[k+i*3]) * 2.0;
    }
    d__1 = f1[0];
    d__2 = f1[1];
    d__3 = f1[2];
    d__4 = f2dot[0];
    d__5 = f2dot[1];
    d__6 = f2dot[2];
    sum_f     = SQR(d__1)+SQR(d__2)+SQR(d__3);
    sum_f2dot = SQR(d__4)+SQR(d__5)+SQR(d__6);
    /* should do the same here at T=0 in case f2dot=0 */
    if (sum_f2dot == 0.0) error("You have hit upon an F2DOT=0 case during integration");
    step[i] = sqrt(eta * sqrt(sum_f/sum_f2dot));
    dprintf(1,"INT1: i=%d time=%g sum_f=%g sum_f2dot=%g -> %g\n",
	    i,time,sum_f,sum_f2dot,step[i]);
    ++nsteps;
    if (time - tnext >= 0.0) {
      goto L100;               // energy check and output
    } else {
      goto L200;               // next timestep
    }
} /* nbody0 */


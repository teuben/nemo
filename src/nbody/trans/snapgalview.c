/*
 * SNAPGALVIEW.C: coordinate transformations for galactic viewing
 *
 *      16-feb-03   V1.0    Created for project w/ R.Olling     Peter Teuben
 *	
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n	    Input snapshot filename",
    "out=???\n	    Output snapshot filename",
    "pos=\n         Position vector of LSR",
    "vel=\n         Velocity vector of LSR",
    "center=\n      Star # to represent LSR (0=first)",
    "times=all\n    Times to select snapshots from",
    "VERSION=1.0\n  16-feb-03 PJT",
    NULL,
};

string usage="coordinate transformations for galactic viewing";

#define TIMEFUZZ	0.000001	/* tolerance in time comparisons */

void nemo_main()
{
    stream instr, outstr;
    string times, ctypei, ctypeo;
    vector pos, vel;
    real   tsnap, phi, r2, v2,sinp, cosp, xnew, ynew;
    int i, nbody, bits, mode, center;
    Body *btab = NULL, *bp;

    times = getparam("times");
    
    instr = stropen(getparam("in"), "r");   
    outstr = stropen(getparam("out"), "w");
    if (hasvalue("center")) {
      center = getiparam("center");
    } else {
      center = -1;
      if (!hasvalue("pos") || !hasvalue("vel")) error("Need pos= and vel=");
      if ((i=nemoinpr(getparam("pos"),pos,3)) != 3) error("Need 3 values for pos=, found %d",i);
      if ((i=nemoinpr(getparam("vel"),vel,3)) != 3) error("Need 3 values for vel=, found %d",i);
    }

    get_history(instr);
    put_history(outstr);		
    for (;;) {                /* loop over all snapshots */
    	get_history(instr);		/* skip over stuff we can forget */
        if (!get_tag_ok(instr, SnapShotTag))
		break;			/* done with work in loop */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0) {
	    continue;       /* just skip it's probably a diagnostics */
        }

        if ((bits & TimeBit) == 0)
	    tsnap = 0.0;
        else if (!streq(times,"all") && !within(tsnap, times, TIMEFUZZ))
            continue;
        dprintf (1,"Transforming snapshot at time= %f bits=0x%x\n",tsnap,bits);

	if (center >= 0 && center < nbody) {
	  SETV(pos,Pos(btab+center));
	  SETV(vel,Vel(btab+center));
	} else if (center >= nbody) {
	  error("center=%d is too large, nbody=%d",center,nbody);
	}

	r2 = sqrt(sqr(pos[0])+sqr(pos[1]));     /* rotation angle */
	if (r2 > 0) {
	  cosp = pos[0]/r2;
	  sinp = pos[1]/r2;
	} else
	  phi = 0.0;
	v2 = sqrt(sqr(vel[0])+sqr(vel[1]));    /* LSR planar velocity  */

	
	dprintf(0,"LSR(%d):  Pos:%g %g %g  |%g|    Vel: %g %g %g |%g|\n", center,
		pos[0],pos[1],pos[2],r2,vel[0],vel[1],vel[2],v2);

        for (bp = btab; bp < btab+nbody; bp++) {
	  for (i=0; i<NDIM; i++)                    /* shift POS */
	    Pos(bp)[i] -= pos[i];
	  for (i=0; i<NDIM; i++)                    /* shift VEL */
	    Vel(bp)[i] -= vel[i];
	  xnew = -cosp*Pos(bp)[0] - sinp*Pos(bp)[1];   /* rotate POS */
	  ynew =  sinp*Pos(bp)[0] - cosp*Pos(bp)[1];
	  Pos(bp)[0] = xnew;
	  Pos(bp)[1] = ynew;
	  xnew = -cosp*Vel(bp)[0] - sinp*Vel(bp)[1];   /* rotate VEL */
	  ynew =  sinp*Vel(bp)[0] - cosp*Vel(bp)[1];
	  Vel(bp)[0] = xnew;
	  Vel(bp)[1] = ynew;
        }
        put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    }
}


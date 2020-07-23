/*
 *  SNAPSPIN: gives a N-body snapshot a spin around the Z-axis.
 *            also corrects potential and forces, if they are present
 *
 *	18-Aug-87	V1.0 created				PJT
 *	 8-Mar-88	V1.1 data history added			PJT
 *	 6-jun-88	V1.2 new filstruct			PJT
 *      25-apr-89       V1.3 get/put_snap , omega bug?          PJT
 *	 1-nov-90	V1.4 helpvec				PJT
 *	 2-may-92	V1.4a  usage 
 *	22-jul-92	V1.5 added rotcur=, sign= 		PJT
 *       6-oct-92       V1.6 fix potential, if present          PJT
 *       1-feb-93       V1.7 Phi/Acc corrections too            PJT
 *	28-apr-98	V2.0 also allow radial inflow/outflow	PJT
 *	mar-94 ansi
 *      23-may-02       V2.0a  nemo_file_lines                  pjt
 *      23-nov-03           b  carry over Key and Aux if present pjt
 *      23-jul-20       V2.1 fix prototypes                     pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <spline.h>
#include <table.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n		input file (snapshot)",
    "out=???\n		output filename",
    "omega=0.0\n	rotation (counterclock)",
    "rotcur=\n          Optional rotation/outflow curve spin (table: r,v)",
    "sign=1\n           Sign of rotcur",
    "outflow=f\n        Outflow or Rotation/Spin",
    "times=all\n	times of snapshots to copy",
    "VERSION=2.1\n	23-jul-2020 PJT",
    NULL,
};

string usage="give a snapshot a spin around or outflow from the Z-axis";

string cvsid="$Id$";

#define TIMEFUZZ	0.0001	/* tolerance in time comparisons */

void get_rotcur(string);
real rotcur(real);

void nemo_main(void)
{
    stream instr, outstr;
    real   tsnap, omega, omega2, sign, r;
    string headline = NULL;
    string times;
    Body *btab = NULL, *bp;
    int i, nbody, bits, nrscale, nvscale;
    bool Qrotcur, Qphi, Qacc, Qflow, Qkey, Qaux;

    instr = stropen (getparam("in"), "r");
    omega=getdparam("omega");
    omega2 = omega*omega;
    Qrotcur = hasvalue("rotcur");
    sign = getdparam("sign");
    if (Qrotcur) {
        get_rotcur(getparam("rotcur"));
        if (omega != 0.0) warning("Ignoring spin omega=%g",omega);
    }
    Qflow = getbparam("outflow");
    outstr = stropen (getparam("out"), "w");
    times = getparam("times");    

    get_history(instr);
    put_history(outstr);					    					
    for(;;) {
    	get_history(instr);
        while (get_tag_ok(instr, HeadlineTag))
	    headline = get_string(instr, HeadlineTag);
        if (!get_tag_ok(instr, SnapShotTag))
	    break;			/* done with work in loop */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0) {
            unlink(getparam("out"));
            error("essential data missing  bits = 0x%x\n", bits);
        }
        if ((bits & TimeBit) == 0)
            tsnap = 0.0;
        else if (!streq(times,"all") && !within(tsnap, times, TIMEFUZZ))
            continue;		/* however skip this snapshot */
        dprintf (1,"Snapshot time=%f spun\n",tsnap);
        Qphi = bits & PotentialBit;
        Qacc = bits & AccelerationBit;
	Qkey = bits & KeyBit;
	Qaux = bits & AuxBit;
        for (bp = btab; bp < btab+nbody; bp++) {
            if (Qrotcur) {
                r = sqrt(sqr(Pos(bp)[0])+sqr(Pos(bp)[1]));
                omega = (r > 0.0 ? rotcur(r)/r * sign  : 0.0);
                dprintf(2,"r=%g omega=%g v=%g\n",r,omega,r*omega);
                omega2 = omega * omega;
            }
            if (Qflow) {                        /* radial in/out-flow */
                Vel(bp)[0] += omega*Pos(bp)[0];     /*  vx = vx + x*omega */
                Vel(bp)[1] += omega*Pos(bp)[1];     /*  vy = vy + y*omega */
                Vel(bp)[2] += omega*Pos(bp)[2];     /*  vz = vz + z*omega */
            } else {                            /* rotation */
                Vel(bp)[0] -= omega*Pos(bp)[1];     /*  vx = vx - y*omega */
                Vel(bp)[1] += omega*Pos(bp)[0];     /*  vy = vy + x*omega */
            }
#if 1
            if (Qphi && !Qrotcur && !Qflow)
                Phi(bp) += 0.5*omega2 * (sqr(Pos(bp)[0])+sqr(Pos(bp)[1]));
            if (Qacc && !Qrotcur && !Qflow) {
                Acc(bp)[0] += omega2 * Pos(bp)[0];
                Acc(bp)[1] += omega2 * Pos(bp)[1];
            }
#endif
        }
        bits = bits & (TimeBit | MassBit | PhaseSpaceBit | KeyBit | AuxBit);
        put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    }
}

local real *rad, *vel, *coef;
local int nrad;

/* 
 * The rotation curve must be an ascii table(5NEMO) with 
 * radius in the first column and velocity in the second.
 * Currently one cannot use pipes yet, since nemo_file_lines
 * would return a < 0 number.
 */

void get_rotcur(string name)
{
    stream instr;
    int  i,colnr[2], nmax;
    real *coldat[2];

    nmax = nemo_file_lines(name,0);
    if (nmax<1) error("No data to read from %s",name);
    rad = (real *) allocate(nmax*sizeof(real));
    vel = (real *) allocate(nmax*sizeof(real));
    coef= (real *) allocate(3*nmax*sizeof(real));   /* spline coeff's */

    colnr[0] = 1;       coldat[0] = rad;
    colnr[1] = 2;       coldat[1] = vel;

    instr = stropen(name,"r");
    nrad = get_atable(instr,2,colnr,coldat,nmax);
    strclose(instr);
    if (nrad > 2)
        spline(coef,rad,vel,nrad);
    else {
      //warning("Could only read %d lines from %s - lin intpol",nrad,name);
        if (nrad<1) error("Not enuf lines (%d) for rotcur",nrad);
    }
    for(i=0; i<nrad; i++)
    	dprintf(1,"rad=%g vel=%g\n",rad[i],vel[i]);

}

real rotcur(real r)
{
    if (r<=0.0)                 /* catch easy points */
        return 0.0;

    if (r<=rad[0])               /* catch case to linearly interpolate */
        return r*vel[0]/rad[0];

    if (r>=rad[nrad-1])          /* keplarian beyond last point */
        return vel[nrad-1]*sqrt(rad[nrad-1]/r);

    if (nrad==2) {          /* linear interpolate 2 points */
        return vel[0] + (r-rad[0])*(vel[1]-vel[0])/(rad[1]-rad[0]);
    }

    return seval(r,rad,vel,coef,nrad);      /* spline interpolation */
}


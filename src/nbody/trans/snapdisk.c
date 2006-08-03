/*
 *  SNAPDISK: gives an assumed disk-like snapshot rotation in line
 *            with radial acccellerations
 *
 *       3-aug-06       V1.0 written at GH2006, cloned off snapspin         PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <spline.h>

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
    "VERSION=1.0\n	3-aug-06 PJT",
    NULL,
};

string usage="assign rotation to a disk";

string cvsid="$Id$";

#define TIMEFUZZ	0.0001	/* tolerance in time comparisons */

void get_rotcur(string);
real rotcur(real);

nemo_main()
{
    stream instr, outstr;
    real   tsnap, omega, omega2, sign, r,v,f,v2, cost,sint;
    string headline = NULL;
    string times;
    Body *btab = NULL, *bp;
    int i, nbody, bits, nrscale, nvscale, nzero;
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
	if (!Qrotcur && !Qacc) error("no acc's");
	nzero = 0;
        for (bp = btab; bp < btab+nbody; bp++) {
	  r = sqrt(sqr(Pos(bp)[0])+sqr(Pos(bp)[1]));
	  if (r>0) {
	    cost = Pos(bp)[0]/r;
	    sint = Pos(bp)[1]/r;
	  } else
	    sint = cost = 0.0;
	  if (Qrotcur) {
	    v = rotcur(r);
	    dprintf(2,"r=%g v=%g\n",r,v);
	  } else {
	    v2 = -(Acc(bp)[0]*Pos(bp)[0] + Acc(bp)[1]*Pos(bp)[1]);
	    if (v2>0.0) 
	      v = sqrt(v2);
	    else {
	      v = 0.0;
	      nzero++;
	    }
	  }
	  Vel(bp)[0] = -v*sint*sign;
	  Vel(bp)[1] =  v*cost*sign;
	}
	if (nzero) warning("%d with 0 velocity",nzero);
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

    nmax = nemo_file_lines(name);
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
        warning("Could only read %d lines from %s - lin intpol",nrad,name);
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


/*
 * MKOMMOD: generate N-body realizations of a spherical system with a
 * spheroidal anisotropic distribution function a la Osipkov-Merritt.
 *
 *	22-jan-89   V1.1                    Josh Barnes
 *      15-nov-90   V1.2  helpvec, Nemo 2.x               PJT
 *	20-feb-92   V1.2a usage			          PJT
 *      22-nov-93   V1.2b cosmetic mods for deVauc table  PJT
 *                        and some sanity checks on table     
 *	11-jan-95   V1.2c proper declarations for Linux	  PJT
 *      29-mar-97   V1.2D SINGLEPREC                      PJT
 *       3-jul-97   V1.2e minor code cleanup              pjt
 *	 5-jul-97   V2.0  added epsilon=		  pjt
 *	21-jul-97       a MXTB is now 1024, from 512	  pjt
 *       1-apr-01       b fixed some compiler warnings    pjt 
 *       9-sep-01       c gsl/xrandom
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/barebody.h>
#include <snapshot/put_snap.c>

#include <spline.h>

string defv[] = {
    "in=???\n			  file with tables of mass, ... vs radius",
    "out=???\n			  snapshot output file name",
    "nbody=512\n		  number of particles to generate",
    "seed=0\n			  seed for random number generator",
    "ftrun=0.999\n		  fractional mass truncation paramete",
    "zerocm=true\n		  if true, zero center of mass",
    "nmodel=1\n			  number of copies to generate",
    "headline=\n		  random verbiage for output file",
    "epsilon=1.0e-10\n		  roundoff control parameter",
    "VERSION=2.0b\n		  1-apr-01 PJT",
    NULL,
};

string usage = 
 "spherical system with anisotropic distribution function a la Osipkov-Merritt";

/*
 * The structure of the model to generate is specified by the following.
 */

#define MXTB 1024		/* max number of values per table             */

local real anisorad;		/* anisotropy radius (if le 0, isotropic)     */
local int ntab;			/* count of tabulated values                  */

local real mtab[MXTB];		/* table of interior masses                   */
local real rtab[MXTB], rcof[3*MXTB]; /* table of radii, coefs for r(M)         */
local real ptab[MXTB], pcof[3*MXTB]; /* table of potentials, coefs for phi(r)  */
local real ftab[MXTB], fcof[3*MXTB]; /* table of distrib. func, coefs for f(Q) */

/*
 * The resulting N-body system is stored in the following structure:
 */

local int nbody;			/* number of bodies in system        */
local Body *btab = NULL;		/* pointer to array of bodies        */

local int count1=0, count2=0;		/* counters fo vnpick 		     */
local real epsilon;			/* roundoff controll		     */

extern double xrandom(double,double);

local real rad_m(real), phi_r(real), pickspeed(void);
local real g_v(real), f_e(real), vnpick(rproc, real, real, real, string);
local readmodel(string);
local anisomod(void);
local initgmax(void);
local snapcenter(void);
local snapwrite(stream);

void nemo_main()
{
    int nmodel, seed;
    stream outstr;
    char hisline[80];

    readmodel(getparam("in"));			/* read tabulated model     */

    nbody = getiparam("nbody");
    if (nbody < 1) error("nbody=%d is absurd",nbody);
    epsilon = getdparam("epsilon");
	
    btab = (Body *) allocate(nbody * sizeof(Body));
    						/* allocate body array      */
    seed = init_xrandom(getparam("seed"));
    nmodel = getiparam("nmodel");
    if (nmodel < 1) error("nmodel=%d is absurd",nmodel);
	
    outstr = stropen(getparam("out"), "w");
    sprintf(hisline,"init_xrandom: seed used %d",seed);
    put_string(outstr, HeadlineTag, hisline);
    if (hasvalue("headline"))
	set_headline(getparam("headline"));
    put_history(outstr);
    while (--nmodel >= 0) {
	anisomod();
	snapwrite(outstr);
    }
    strclose(outstr);
    dprintf(1,"vnpick: von Neumann rejection %d/%d = %g\n",
		count2,count1, (real)count2/(real)count1);
}

/*
 * READMODEL: initalize tables for radius, enclosed mass, potential,
 * anisotropic "energy" and distribution function from a file.
 */

local readmodel(string name)
{
    stream instr;
    int i;

    instr = stropen(name, "r");
    get_history(instr);
    get_set(instr, "OsipkovMerrittModel");
    get_data_coerced(instr, "AnisoRadius", RealType, &anisorad, 0);
    get_data(instr, "Ntab", IntType, &ntab, 0);
    if (ntab > MXTB) error("Not enough space for tables, MXTB=%d",MXTB);
    get_data_coerced(instr, "Radius", RealType, rtab, ntab, 0);
    get_data_coerced(instr, "Mass", RealType, mtab, ntab, 0);
    get_data_coerced(instr, "Potential", RealType, ptab, ntab, 0);
    get_data_coerced(instr, "DistFunc", RealType, ftab, ntab, 0);
    get_tes(instr, "OsipkovMerrittModel");
    strclose(instr);

    dprintf(0, "[readmode: ntab = %d]\n", ntab);
    for (i=1; i<ntab; i++) {
      if (rtab[i-1] > rtab[i]) warning("Radius %d not increasing?",i+1);
      if (mtab[i-1] > mtab[i]) warning("Mass %d not increasing?",i+1);
      if (ptab[i-1] > ptab[i]) warning("Potential %d not increasing?",i+1);
      if (ftab[i-1] < ftab[i]) warning("DistFunc %d not decreasing?",i+1);
    }
}

/*
 * ANISOMOD: generate a realization of an anisotropic model.
 * Note: To save time in the rejection phase, the realization is
 * truncated at the radius enclosing a fraction ftrun of the mass.
 */

local real phix;		/* for pickspeed and g_v */

local anisomod(void)
{
    real ftrun, x, rx, vx, svt, vr;
    vector vrad;
    Body *bp;

    ftrun = getdparam("ftrun");			/* get trunc. mass fraction */
    spline(rcof, mtab, rtab, ntab);		/* set up r = r(M)          */
    spline(pcof, rtab, ptab, ntab);		/* set up phi = phi(r)      */
    spline(fcof, ptab, ftab, ntab);		/* set up f = f(Q)          */
    initgmax();					/* initialize gmax(phix)    */
    for (bp = btab; bp < btab + nbody; bp++) {	/* loop over body array     */
	x = xrandom(mtab[0],ftrun * mtab[ntab-1]); /*   pick mass enclosed  */	
	rx = rad_m(x);				/*   pick radius from r(M)  */
	pickshell(Pos(bp), NDIM, rx);		/*   pick rand point at rx  */
	phix = phi_r(rx);			/*   find potential  at rx  */
	vx = pickspeed();			/*   pick speed from g(v)   */
	pickshell(Vel(bp), NDIM, vx);		/*   pick rand velocity vec */
	if (anisorad > 0.0) {			/*   sqeeze vel distrib?    */
	    svt = 1.0 / sqrt(1 + sqr(rx / anisorad));
	    					/*     find trans v. scale  */
	    vr = dotvp(Vel(bp), Pos(bp)) / rx;	/*     and radial vel comp  */
	    MULVS(vrad, Pos(bp), vr/rx);	/*     and radial proj of v */
	    MULVS(Vel(bp), Vel(bp), svt);	/*     scale transverse vel */
	    MULVS(vrad, vrad, 1.0 - svt);	/*     and radial vel       */
	    ADDV(Vel(bp), Vel(bp), vrad);	/*     put parts together   */
	}
	Mass(bp) = mtab[ntab-1] / nbody;	/*   assign equal masses    */
    }
    if (getbparam("zerocm"))			/* zero center of mass?     */
	snapcenter();
}

/*
 * INITGMAX: initalize the table used to find gmax(phi).
 */

#define NGMAX 64	/* size of gmax[phi] table */
#define NGPNT 24	/* points evaluated to find gmax */

local real pmax[NGMAX], gmax[NGMAX], gcoef[3*NGMAX];

local initgmax(void)
{
    int i, j;
    real vmax, v, g;

    for (i = 0; i < NGMAX; i++) {
	pmax[i] = ptab[0] + i * (ptab[ntab-1] - ptab[0]) / (NGMAX - 1.0);
	phix = pmax[i];
	vmax = sqrt(ABS(-2.0 * (pmax[i] - ptab[ntab-1])));
	gmax[i] = 0.0;
	for (j = 1; j <= NGPNT; j++) {
	    v = j * vmax / NGPNT;
	    g = g_v(v);
	    if (g > gmax[i])
		gmax[i] = g;
	}
    }
    spline(gcoef, pmax, gmax, NGMAX);
}

/*
 * PICKSPEED: chose speed distributed according to g(v) = v^2 f(v;phix).
 */

local real pickspeed()
{
    real vm, gm;

    vm = sqrt(-2.0 * (phix - ptab[ntab-1]));
    gm = seval(phix, pmax, gmax, gcoef, NGMAX);
    return vnpick(g_v, 0.0, vm, gm, "g_v");
}

/*
 * G_V: compute speed distribution (up to a factor of 4 PI).
 */

local real g_v(real v)
{
    return v*v * f_e(0.5 * v*v + phix);
}

/*
 * F_E: compute distribution as a function of E = 0.5 v^2 + phi.
 */

local real f_e(real e)
{
    /* check if valid energy, but allow for roundoff error */
    if (ptab[0] <= e && e <= ptab[ntab-1] + epsilon)
	return (seval(e, ptab, ftab, fcof, ntab));
    else
	error("f_e: e  = %g outside permitted range [%g,%g]",
		e,ptab[0],ptab[ntab-1]+epsilon);
}

/*
 * RAD_M: radius as a function of mass.
 */

local real rad_m(real m)
{
    return seval(m, mtab, rtab, rcof, ntab);
}

/*
 * PHI_R: potential as a function of radius.
 */

local real phi_r(real r)
{
    return seval(r, rtab, ptab, pcof, ntab);
}

/*
 * VNPICK: chose from a distribution by von Neumann technique.
 * Ref: von Neumann, J. 1963. Collected Works, Vol. 5.
 *
 * fun      :		distribution function
 * xl, xh   :		range of values to consider
 * fh       :		max value of f to generate
 * name     :		function name, for warning or error
 */

local real vnpick(rproc fun, real xl, real xh, real fh, string name)
{
    real f, fx, x;

    f = 1.0;
    fx = 0.0;
    count1++;
    while (f > fx) {
    	count2++;
	x = xrandom(xl, xh);
	f = xrandom(0.0, 1.1 * fh);
	fx = (*fun)(x);
	if (fx > 1.05 * fh) {
	    dprintf(0,"vnpick: x = %g\t%s(x) = %g\t%s_max = %g\n",
		   x, name, fx, name, fh);
	    if (fx > 1.1 * fh)
	        error("vnpick: %s(x) out of bounds",name);
	}
    }
    return x;
}

/*
 * SNAPCENTER: transform to center-of mass coordinates.
 */

local snapcenter(void)
{
    vector cmpos, cmvel, tmp;
    real mtot;
    Body *bp;

    CLRV(cmpos);
    CLRV(cmvel);
    mtot = 0.0;
    for (bp = btab; bp < btab + nbody; bp++) {
	mtot = mtot + Mass(bp);
	MULVS(tmp, Pos(bp), Mass(bp));
	ADDV(cmpos, cmpos, tmp);
	MULVS(tmp, Vel(bp), Mass(bp));
	ADDV(cmvel, cmvel, tmp);
    }
    DIVVS(cmpos, cmpos, mtot);
    DIVVS(cmvel, cmvel, mtot);
    for (bp = btab; bp < btab + nbody; bp++) {
	SUBV(Pos(bp), Pos(bp), cmpos);
	SUBV(Vel(bp), Vel(bp), cmvel);
    }
    dprintf(0, "[snapcenter in %s: mtot = %f",getargv0(), mtot);
    dprintf(0, "    cmpos = %f    cmvel = %f]\n",absv(cmpos), absv(cmvel));
}

/*
 * SNAPWRITE: output completed N-body model.
 */

local snapwrite(stream outstr)
{
    real tzero = 0.0;
    int bits = TimeBit | MassBit | PhaseSpaceBit;

    put_snap(outstr, &btab, &nbody, &tzero, &bits);
}

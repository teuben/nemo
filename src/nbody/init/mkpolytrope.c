/*
 * MKPOLYTROPE: generate an N-body realization of one of Henon's
 * generalized polytropes.
 * 
 *	22-jan-89   V1.1 				Josh Barnes
 *	15-nov-90   V1.2  Nemo 2.x                      PJT
 *	13-dec-94   V1.3  renamed 4th parameter nobj= to 2nd nbody=	PJT
 *                        default seed is now 0
 *	29-mar-97   V1.3a usage, SINGLEPREC		pjt
 *       9-sep-01       b gsl/xrandom
 *      19-feb-06    1.4  mdarray                       pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <mdarray.h>

#include <snapshot/snapshot.h>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "out=???\n			  output file name ",
    "nbody=512\n		  number of particles ",
    "n=1\n			  energy exponent ",
    "m=-1\n			  angular momentum exponent ",
    "seed=0\n			  seed for random numbers ",
    "headline=\n		  verbiage for output ",
    "VERSION=1.4\n		  19-feb-06 PJT",
    NULL,
};

string usage="N-body realization of one of Henon's generalized polytropes";

string cvsid="$Id$";

local string headline;		/* random text message */

/*
 * The model to be generated is characterized by two parameters:
 */

local real nval;			/* energy index: 0.5 <= nval <= inf. */
local real mval;			/* angmom index:  -1 <= mval <= inf. */

/*
 * The generated N-body model is stored in the following:
 */

local int nobj;			/* number of particles */
local real *mass;		/* masses of particles generated */
local mdarray2 phase;		/* phase-space coordinates of particles */

extern double xrandom(double, double);

extern void pickshell(real *,int,double);

nemo_main()
{
    string oname;		/* files for model, output */

    nobj = getiparam("nbody");
    mass = (real *) allocate(nobj*sizeof(real));
    phase = allocate_mdarray2(nobj,6);
    init_xrandom(getparam("seed"));
    nval = getdparam("n");
    mval = getdparam("m");
    oname = getparam("out");
    headline = getparam("headline");
    if (nval == 1.0 && mval == -1.0)	/* exceptional case n=1, m=-1? */
        polymod0();
    else
        error("general polytropic models not implemented yet");
    writesnap(oname);
}

#define SQRT2 1.414214

polymod0()
{
    int i, j;
    double ri, x, vi;

    for (i = 0; i < nobj; i++) {		/* loop generating bodies */
        ri = xrandom(0.0, 2.0);			/*   pick uniform dist */
	pickshell(&phase[i][0], 3, ri);		/*   pick point at ri */
	x = SQRT2 * cos(PI * (2.0 - xrandom(0.0, 1.0)) / 4.0);
	vi = (1.0 - x*x) * sqrt(log(2.0 / ri));	/*   pick speed (BGH A15) */
	if (xrandom(0.0, 1.0) < 0.5)		/*   pick sign at random */
	    vi = -vi;
	for (j = 0; j < 3; j++)			/*   loop setting velocity */
	    phase[i][j+3] = vi * phase[i][j] / ri;
	mass[i] = 1.0 / nobj;			/*   total mass M is unity */
    }
    zerocms(phase[0], 6, mass, nobj, nobj);
}

writesnap(string name)
{
    stream outstr;
    int coord = CSCode(Cartesian, 3, 2);

    if (! streq(headline, ""))
	set_headline(headline);
    outstr = stropen(name, "w");
    put_history(outstr);
    put_set(outstr, SnapShotTag);
    put_set(outstr, ParametersTag);
    put_data(outstr, NobjTag, IntType, &nobj, 0),
    put_data(outstr, "nval", RealType, &nval, 0),
    put_data(outstr, "mval", RealType, &mval, 0),
    put_tes(outstr, ParametersTag);
    put_set(outstr, ParticlesTag);
    put_data(outstr, CoordSystemTag, IntType, &coord, 0);
    put_data(outstr, MassTag, RealType, mass, nobj, 0);
    put_data(outstr, PhaseSpaceTag, RealType, phase[0], nobj, 2, 3, 0),
    put_tes(outstr, ParticlesTag);
    put_tes(outstr, SnapShotTag);
    strclose(outstr);
}

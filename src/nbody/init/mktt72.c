/* 
 * MKTT
 *
 *   19-nov-2002	Created
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "out=???\n		output file name",
    "nbody=100\n	number of particles per ring",
    "radius=1:6:1\n	radii of rings",
    "seed=0\n		seed for random numbers",
    "zerocm=false\n	if true, zero center of mass",
    "headline=\n	verbiage for output",
    "VERSION=0.1\n	19-nov-02 PJT",
    NULL,
};

string usage = "Create a Toomre & Toomre 1972 test disk";


#ifndef MOBJ
#  define MOBJ 4096
#endif

#define MAXRAD 1024


local int nobj;
local real mass[MOBJ];
local vector phase[MOBJ][2];
local double radius[MAXRAD];

local stream outstr;
local string headline;

extern double xrandom(double, double);



void nemo_main()
{
    int i, nrad;

    nrad = nemoinpd(getparam("radius"),radius,MAXRAD);
    nobj = getiparam("nbody");
    if (nobj > MOBJ)
	error("Too many particles requested: nbody > MOBJ [%d]", MOBJ);
    init_xrandom(getparam("seed"));
    headline = getparam("headline");

    for (i=0; i<nrad; i++) {
        makering(radius[i]);
        if (getbparam("zerocm"))
	    zerocms(phase, 2 * NDIM, mass, nobj, nobj);
        writesnap(i);
    }
    strclose(outstr);
}

makering(real radius)
{
    int i;
    real theta;

    for (i = 0; i < nobj; i++) {
	mass[i] = 1.0 / nobj;
	theta = TWO_PI * ((real) i) / nobj;
	CLRV(phase[i][0]);
	phase[i][0][0] = radius * sin(theta);
	phase[i][0][1] = radius * cos(theta);
	CLRV(phase[i][1]);
    }
}

writesnap(int i)
{

    int cs = CSCode(Cartesian, NDIM, 2);

    if (i==0) {
        if (! streq(headline, ""))
            set_headline(headline);
        outstr = stropen(getparam("out"), "w");
        put_history(outstr);
    }

    put_set(outstr, SnapShotTag);
     put_set(outstr, ParametersTag);
      put_data(outstr, NobjTag, IntType, &nobj, 0);
     put_tes(outstr, ParametersTag);
     put_set(outstr, ParticlesTag);
      put_data(outstr, CoordSystemTag, IntType, &cs, 0);
      put_data(outstr, MassTag, RealType, mass, nobj, 0);
      put_data(outstr, PhaseSpaceTag, RealType, phase, nobj, 2, NDIM, 0);
     put_tes(outstr, ParticlesTag);
    put_tes(outstr, SnapShotTag);

}

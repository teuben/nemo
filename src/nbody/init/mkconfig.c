/* 
 * MKCONFIG: make a static configuration of particles.
 *
 *	22-jan-89	JEB	???
 *	 8-jul-90  0.1	PJT	using <snapshot/...> macro's now	
 *	20-feb-92  0.2	PJT	usage, nemo_main
 *	23-mar-97  0.2b pjt     fixed protos
 *	27-mar-97  0.3  pjt	moved nbody= as 2nd keyword
 *       9-sep-01       a       gsl/xrandom
 *      10-feb-06  0.4  pjt     using mdarray, no more need for MOBJ, but slightly ugly syntax
 *      13-sep-09  0.4b pjt     added snapshot time
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <mdarray.h>
#include <history.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "out=???\n		output file name",
    "nbody=512\n	number of particles",
    "shape=shell\n	shape of configuration {shell,ball,cusp,line,ring}",
    "radius=1.0\n	radius of shell, ball or ring",
    "seed=123\n		seed for random numbers",
    "mass=1\n           total mass of system",
    "zerocm=false\n	if true, zero center of mass",
    "headline=\n	verbiage for output",
    "VERSION=1.0\n	5-may-10 PJT",
    NULL,
};

string usage = "make a static configuration of particles";

string cvsid="$Id$";


local int nobj;
local real *mass;
local mdarray3 phase;
local double total_mass;                
local double radius;	                /* must be double !! */

extern double xrandom(double, double);
extern void  zerocms(double *, int, double *, int, int);
void makeshell(void);
void makeball(void);
void makecusp(void);
void makeline(void);
void makering(void);
void writesnap(void);


void nemo_main()
{
    string shape = getparam("shape");

    total_mass = getdparam("mass");
    radius = getdparam("radius");
    nobj = getiparam("nbody");

    mass = (real *) allocate(nobj * sizeof(real));
    phase = allocate_mdarray3(nobj,2,NDIM);

    init_xrandom(getparam("seed"));
    if (streq(shape, "shell"))
	makeshell();
    else if (streq(shape, "ball"))
	makeball();
    else if (streq(shape, "cusp"))
	makecusp();
    else if (streq(shape, "line"))
	makeline();
    else if (streq(shape, "ring"))
	makering();
    else
	error("%s: shape %s unknown\n", getargv0(), shape);
    if (total_mass > 0 && getbparam("zerocm"))
	zerocms(phase[0][0], 2 * NDIM, mass, nobj, nobj);
    writesnap();
}

void makeshell()
{
    int i;

    for (i = 0; i < nobj; i++) {
	mass[i] = total_mass / nobj;
	pickshell(phase[i][0], NDIM, (double)radius);
	CLRV(phase[i][1]);
    }
}

void makeball()
{
    int i;

    for (i = 0; i < nobj; i++) {
	mass[i] = total_mass / nobj;
	pickball(phase[i][0], NDIM, radius);
	CLRV(phase[i][1]);
    }
}

void makecusp()
{
    int i;

    for (i = 0; i < nobj; i++) {
	mass[i] = total_mass / nobj;
	pickshell(phase[i][0], NDIM, xrandom(0.0, radius));
	CLRV(phase[i][1]);
    }
}

void makeline()
{
    int i;

    for (i = 0; i < nobj; i++) {
	mass[i] = total_mass / nobj;
	CLRV(phase[i][0]);
	phase[i][0][0] = radius * i / (nobj - 1.0);
	CLRV(phase[i][1]);
    }
}

void makering()
{
    int i;
    real theta;

    for (i = 0; i < nobj; i++) {
	mass[i] = total_mass / nobj;
	theta = TWO_PI * ((real) i) / nobj;
	CLRV(phase[i][0]);
	phase[i][0][0] = radius * sin(theta);
	phase[i][0][1] = radius * cos(theta);
	CLRV(phase[i][1]);
    }
}

void writesnap()
{
    stream outstr;
    string headline = getparam("headline");
    int cs = CSCode(Cartesian, NDIM, 2);
    real tsnap = 0.0;

    if (! streq(headline, ""))
	set_headline(headline);
    outstr = stropen(getparam("out"), "w");

    put_history(outstr);
    put_set(outstr, SnapShotTag);
     put_set(outstr, ParametersTag);
      put_data(outstr, NobjTag, IntType, &nobj, 0);
      put_data(outstr, TimeTag, RealType, &tsnap, 0);
     put_tes(outstr, ParametersTag);
     put_set(outstr, ParticlesTag);
      put_data(outstr, CoordSystemTag, IntType, &cs, 0);
      put_data(outstr, MassTag, RealType, mass, nobj, 0);
      put_data(outstr, PhaseSpaceTag, RealType, phase[0][0], nobj, 2, NDIM, 0);
     put_tes(outstr, ParticlesTag);
    put_tes(outstr, SnapShotTag);
    strclose(outstr);
}
